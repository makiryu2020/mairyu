/*
 * 電解質モニター 計算コア（純粋関数・UIなし）
 *
 * 方針（引き継ぎまとめ §5・§8）:
 *  - 数値解釈はこの決定論的コアだけが行う。UI や生成AIはここの結果を表示・説明するだけ。
 *  - 体液量は入力に含めない。
 *  - SIG は文献標準の SIDa − SIDe（正値＝未測定陰イオン, Kellum 1995）に統一。
 *
 * 入力単位（正準単位）:
 *  pH, pCO2 mmHg, HCO3 mmol/L, Na/K/Cl mmol/L, Lac mmol/L, Alb g/dL,
 *  iCa mmol/L, tCa mg/dL, Mg mg/dL, P mg/dL, Glu mg/dL, BUN mg/dL, Osm mOsm/kg
 *
 * ブラウザでは window.EMCore、Node では module.exports として使う。
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.EMCore = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ---------------------------------------------------------------- 定数
  const S_CO2 = 0.0307;          // mmol/L/mmHg, 37℃
  const PKA = 6.1;
  const MW_P = 30.97;            // mg/dL → mmol/L は ÷3.097
  const MW_MG = 24.305;          // mg/dL → mmol/L は ÷2.4305
  const MW_CA = 40.08;           // mg/dL → mmol/L は ÷4.008
  const H24 = 24 * 3600 * 1000;

  // 既定設定。基準範囲は施設・測定法依存なので UI から変更できる。
  // src: 出典。"memory" は一次文献を本実装で直接確認していない値。
  const DEFAULT_SETTINGS = {
    agRef: 12,            // ΔΔ計算の基準AG
    agUpper: 16,          // AG開大の判定上限（12±4）。ISE法では低めになるので施設値に合わせる
    glucoseFactor: 2.4,   // 補正Na係数（Hillier 1999: 2.4 / Katz 1973: 1.6）
    naLimit24h: 10,       // 低Na補正の24h上限（ESE/ERA-EDTA 2014: 初日10, 以降8）
    naHighRisk: false,    // 高リスク（ODS）では 8（Verbalis 2013）
    rcvZ: 1.96,           // 両側95%
    ranges: {
      pH:   { lo: 7.35, hi: 7.45, unit: '' },
      pCO2: { lo: 35, hi: 45, unit: 'mmHg' },
      HCO3: { lo: 22, hi: 26, unit: 'mmol/L' },
      Na:   { lo: 135, hi: 145, unit: 'mmol/L' },
      K:    { lo: 3.5, hi: 5.0, unit: 'mmol/L' },
      Cl:   { lo: 98, hi: 106, unit: 'mmol/L' },
      Lac:  { lo: 0.5, hi: 2.0, unit: 'mmol/L' },
      Alb:  { lo: 4.0, hi: 5.0, unit: 'g/dL' },
      iCa:  { lo: 1.15, hi: 1.30, unit: 'mmol/L' },
      tCa:  { lo: 8.6, hi: 10.2, unit: 'mg/dL', src: 'memory' },
      Mg:   { lo: 1.7, hi: 2.4, unit: 'mg/dL', src: 'memory' },
      P:    { lo: 2.5, hi: 4.5, unit: 'mg/dL' },
      Glu:  { lo: 70, hi: 140, unit: 'mg/dL', src: 'memory' },
    },
    // 個体内生物学的変動 CVI（%）: EFLM Biological Variation Database（検索要約経由・原典で要確認）
    // CVA（%）: 既定は「望ましい分析精度」0.5×CVI（Fraser 2001）。自施設の実測CVへの置換を推奨
    bv: {
      Na:  { cvi: 0.5, cva: 0.25 },
      K:   { cvi: 4.1, cva: 2.05 },
      Cl:  { cvi: 1.1, cva: 0.55 },
      tCa: { cvi: 1.8, cva: 0.9 },
      Mg:  { cvi: 2.9, cva: 1.45 },
      P:   { cvi: 7.8, cva: 3.9 },
    },
  };

  // ---------------------------------------------------------------- 小道具
  const isNum = (v) => typeof v === 'number' && Number.isFinite(v);
  const has = (s, ...keys) => keys.every((k) => isNum(s[k]));
  const r = (v, d = 1) => (isNum(v) ? Math.round(v * 10 ** d) / 10 ** d : null);

  function mergeSettings(user) {
    const s = JSON.parse(JSON.stringify(DEFAULT_SETTINGS));
    if (!user) return s;
    for (const k of Object.keys(user)) {
      if (k === 'ranges' || k === 'bv') {
        for (const a of Object.keys(user[k] || {})) s[k][a] = Object.assign({}, s[k][a], user[k][a]);
      } else s[k] = user[k];
    }
    return s;
  }

  // ---------------------------------------------------------------- 単位換算
  const conv = {
    pMgdlToMmol: (v) => v / (MW_P / 10),
    mgMgdlToMmol: (v) => v / (MW_MG / 10),
    caMgdlToMmol: (v) => v / (MW_CA / 10),
    albGdlToGl: (v) => v * 10,
  };

  // ---------------------------------------------------------------- 基本計算
  /** Henderson–Hasselbalch による HCO3⁻ */
  function hhHCO3(pH, pCO2) {
    return S_CO2 * pCO2 * 10 ** (pH - PKA);
  }

  /** SBE（Van Slyke, Hb 5 g/dL 相当の細胞外液） */
  function sbe(pH, hco3) {
    return 0.9287 * (hco3 - 24.4 + 14.83 * (pH - 7.4));
  }

  function anionGap(na, cl, hco3) {
    return na - (cl + hco3);
  }

  /** アルブミン補正AG（Figge 1998: 2.5 mEq/L per g/dL） */
  function correctedAG(ag, albGdl) {
    return isNum(albGdl) ? ag + 2.5 * (4.0 - albGdl) : ag;
  }

  /** グルコース補正Na */
  function correctedNa(na, glu, factor) {
    return glu > 100 ? na + factor * (glu - 100) / 100 : na;
  }

  /** Payne 補正Ca（mg/dL）。Alb < 4 のときのみ補正 */
  function correctedCa(tca, albGdl) {
    return isNum(albGdl) && albGdl < 4 ? tca + 0.8 * (4.0 - albGdl) : tca;
  }

  /** 計算浸透圧 */
  function calcOsm(na, glu, bun) {
    return 2 * na + (isNum(glu) ? glu / 18 : 0) + (isNum(bun) ? bun / 2.8 : 0);
  }

  // ---------------------------------------------------------------- Stewart
  /**
   * Figge 簡易式。Alb g/L、Pi mmol/L。
   * Ca/Mg は2価なので mmol/L ×2 で mEq/L に換算する。
   * iCa がなければ SIDa は算出しない（総Caで代用すると SIDa を過大評価するため）。
   * Mg は総Mgを使用（イオン化Mgは通常未測定。SIDa をわずかに過大評価する）。
   */
  function stewart(s) {
    if (!has(s, 'pH', 'pCO2', 'Na', 'K', 'Cl', 'Alb', 'P')) return null;
    const hco3 = isNum(s.HCO3) ? s.HCO3 : hhHCO3(s.pH, s.pCO2);
    const albGl = conv.albGdlToGl(s.Alb);
    const pi = conv.pMgdlToMmol(s.P);
    const albCharge = albGl * (0.123 * s.pH - 0.631);
    const piCharge = pi * (0.309 * s.pH - 0.469);
    const sidE = hco3 + albCharge + piCharge;
    const out = { albCharge, piCharge, atotCharge: albCharge + piCharge, sidE, sidA: null, sig: null, notes: [] };
    if (isNum(s.iCa)) {
      const mg = isNum(s.Mg) ? 2 * conv.mgMgdlToMmol(s.Mg) : 0;
      if (!isNum(s.Mg)) out.notes.push('Mg 未入力のため SIDa から除外');
      out.sidA = s.Na + s.K + 2 * s.iCa + mg - s.Cl - (isNum(s.Lac) ? s.Lac : 0);
      if (!isNum(s.Lac)) out.notes.push('乳酸未入力のため 0 として計算');
      out.sig = out.sidA - sidE;
    } else {
      out.notes.push('イオン化Ca 未入力のため SIDa・SIG は算出しない');
    }
    return out;
  }

  /** BE 分解（Story 2004 簡易法）。乳酸は残差から分けて示す */
  function beDecomposition(s, hco3) {
    if (!has(s, 'pH', 'Na', 'Cl', 'Alb')) return null;
    const total = sbe(s.pH, hco3);
    const naCl = s.Na - s.Cl - 38;
    const alb = 0.25 * (42 - conv.albGdlToGl(s.Alb));
    const lac = isNum(s.Lac) ? -s.Lac : 0;
    const other = total - naCl - alb - lac;
    return { sbe: total, naCl, alb, lac, other };
  }

  // ---------------------------------------------------------------- 酸塩基分類（Boston法）
  /**
   * ヒューリスティック分類。三重障害は「示唆」まで。最終判断は臨床医。
   * 代償式は Adrogué & Madias の慣用式（引き継ぎまとめ §6-2）。
   */
  function classifyAcidBase(s, hco3, cag, st) {
    const out = { primary: [], secondary: [], comp: null, flags: [] };
    const { pH, pCO2 } = s;
    const R = st.ranges;
    const acidemia = pH < R.pH.lo, alkalemia = pH > R.pH.hi;
    const lowH = hco3 < R.HCO3.lo, highH = hco3 > R.HCO3.hi;
    const highC = pCO2 > R.pCO2.hi, lowC = pCO2 < R.pCO2.lo;

    const metAcid = () => {
      const exp = 1.5 * hco3 + 8;
      out.comp = { label: "Winter式 予測pCO₂", expected: exp, range: [exp - 2, exp + 2], unit: 'mmHg' };
      if (pCO2 > exp + 2) out.secondary.push('呼吸性アシドーシス合併（代償不十分）');
      else if (pCO2 < exp - 2) out.secondary.push('呼吸性アルカローシス合併');
    };
    const metAlk = () => {
      const exp = 40 + 0.7 * (hco3 - 24);
      out.comp = { label: '代謝性アルカローシス 予測pCO₂', expected: exp, range: [exp - 5, exp + 5], unit: 'mmHg' };
      if (pCO2 > exp + 5) out.secondary.push('呼吸性アシドーシス合併');
      else if (pCO2 < exp - 5) out.secondary.push('呼吸性アルカローシス合併');
    };
    const respAcid = () => {
      const d = pCO2 - 40, a = 24 + 0.1 * d, c = 24 + 0.35 * d;
      out.comp = { label: '呼吸性アシドーシス 予測HCO₃⁻（急性〜慢性）', expected: null, range: [a, c], unit: 'mmol/L' };
      if (hco3 < a - 2) out.secondary.push('代謝性アシドーシス合併');
      else if (hco3 > c + 2) out.secondary.push('代謝性アルカローシス合併');
      else out.flags.push(hco3 < (a + c) / 2 ? '急性寄り' : '慢性寄り');
    };
    const respAlk = () => {
      const d = 40 - pCO2, a = 24 - 0.2 * d, c = 24 - 0.5 * d;
      out.comp = { label: '呼吸性アルカローシス 予測HCO₃⁻（急性〜慢性）', expected: null, range: [c, a], unit: 'mmol/L' };
      if (hco3 < c - 2) out.secondary.push('代謝性アシドーシス合併');
      else if (hco3 > a + 2) out.secondary.push('代謝性アルカローシス合併');
      else out.flags.push(hco3 > (a + c) / 2 ? '急性寄り' : '慢性寄り');
    };

    if (acidemia) {
      if (lowH) { out.primary.push('代謝性アシドーシス'); metAcid(); }
      if (highC) { out.primary.push('呼吸性アシドーシス'); if (!lowH) respAcid(); }
      if (!lowH && !highC) out.flags.push('酸血症だが HCO₃⁻・pCO₂ に一次障害を示す変化なし：入力値の整合性を確認');
    } else if (alkalemia) {
      if (highH) { out.primary.push('代謝性アルカローシス'); metAlk(); }
      if (lowC) { out.primary.push('呼吸性アルカローシス'); if (!highH) respAlk(); }
      if (!highH && !lowC) out.flags.push('アルカリ血症だが一次障害を示す変化なし：入力値の整合性を確認');
    } else {
      if (lowH && lowC) out.flags.push('pH正常で HCO₃⁻↓・pCO₂↓：代謝性アシドーシス＋呼吸性アルカローシスの混合を考慮');
      if (highH && highC) out.flags.push('pH正常で HCO₃⁻↑・pCO₂↑：代謝性アルカローシス＋呼吸性アシドーシスの混合を考慮');
    }

    // AG による代謝性アシドーシスの検出（pH・HCO₃⁻ が正常でも隠れている場合がある）
    if (isNum(cag) && cag > st.agUpper) {
      if (!out.primary.includes('代謝性アシドーシス')) out.secondary.push('AG開大性代謝性アシドーシスの併存');
      const dd = (cag - st.agRef) / (24 - hco3);
      if (hco3 < 24) {
        out.deltaRatio = dd;
        if (dd < 0.8) out.secondary.push('非AG開大性（高Cl性）代謝性アシドーシスの併存（ΔΔ<0.8）');
        else if (dd > 2) out.secondary.push('代謝性アルカローシス併存 または 慢性呼吸性アシドーシス（ΔΔ>2）');
      } else {
        out.secondary.push('HCO₃⁻≥24 で AG開大：代謝性アルカローシスの併存を考慮');
      }
    } else if (out.primary.includes('代謝性アシドーシス')) {
      out.secondary.push('非AG開大性（高Cl性）代謝性アシドーシス');
    }
    return out;
  }

  // ---------------------------------------------------------------- 電解質判定
  function band(v, rg) {
    if (!isNum(v) || !rg) return null;
    return v < rg.lo ? 'low' : v > rg.hi ? 'high' : 'normal';
  }

  /**
   * 重症度は出典のある区分だけを付ける。
   * level: 'ok' | 'watch' | 'warn' | 'crit'
   */
  function electrolyteFindings(s, d, st) {
    const R = st.ranges, f = [];
    const add = (key, level, text, src) => f.push({ key, level, text, src });

    if (isNum(s.Na)) {
      const na = isNum(d.correctedNa) ? d.correctedNa : s.Na;
      const lab = isNum(d.correctedNa) && d.correctedNa !== s.Na ? `補正Na ${r(na)}` : `Na ${s.Na}`;
      if (s.Na < 125) add('Na', 'crit', `高度（profound）低Na血症：Na ${s.Na}`, 'ESE/ERA-EDTA 2014');
      else if (s.Na < 130) add('Na', 'warn', `中等度低Na血症：Na ${s.Na}`, 'ESE/ERA-EDTA 2014');
      else if (s.Na < R.Na.lo) add('Na', 'watch', `軽度低Na血症：Na ${s.Na}`, 'ESE/ERA-EDTA 2014');
      else if (s.Na > R.Na.hi) add('Na', 'warn', `高Na血症：Na ${s.Na}`, '重症度区分の統一基準なし');
      if (s.Na < R.Na.lo && isNum(s.Glu) && s.Glu > 100 && na >= R.Na.lo) {
        add('Na', 'watch', `${lab}：高血糖による低Na（高張性/偽性の可能性）`, 'Hillier 1999');
      }
      if (isNum(d.osmGap) && d.osmGap > 10) add('Osm', 'warn', `浸透圧ギャップ ${r(d.osmGap)}（>10）：中毒性アルコール等を考慮`, '');
    }
    if (isNum(s.K)) {
      if (s.K >= 6.5) add('K', 'crit', `重度高K血症：K ${s.K}（≥6.5）`, 'UKKA 2023');
      else if (s.K >= 6.0) add('K', 'warn', `中等度高K血症：K ${s.K}（6.0–6.4）`, 'UKKA 2023（区分は記憶に基づく）');
      else if (s.K >= 5.5) add('K', 'watch', `軽度高K血症：K ${s.K}（5.5–5.9）`, 'UKKA 2023（区分は記憶に基づく）');
      else if (s.K > R.K.hi) add('K', 'watch', `基準上限超：K ${s.K}`, '');
      else if (s.K < 2.5) add('K', 'crit', `重度低K血症：K ${s.K}（<2.5）`, '慣用区分（記憶に基づく）');
      else if (s.K < R.K.lo) add('K', 'warn', `低K血症：K ${s.K}`, '');
      if (isNum(s.pH) && s.pH < 7.2 && s.K < 5.0) {
        add('K', 'warn', 'アシデミアで K が高くない：総体内K欠乏の可能性（インスリン・アルカリ投与で急低下に注意）', 'Kitabchi 2009（DKA）');
      }
    }
    if (isNum(s.Cl)) {
      if (s.Cl > R.Cl.hi) add('Cl', 'watch', `高Cl血症：Cl ${s.Cl}（Na補正Cl ${r(d.naCorrCl)}）`, '');
      else if (s.Cl < R.Cl.lo) add('Cl', 'watch', `低Cl血症：Cl ${s.Cl}（Na補正Cl ${r(d.naCorrCl)}）`, '');
    }
    if (isNum(s.iCa)) {
      const b = band(s.iCa, R.iCa);
      if (b === 'low') add('Ca', s.iCa < 0.9 ? 'crit' : 'warn', `イオン化Ca低値：${s.iCa} mmol/L`, s.iCa < 0.9 ? '<0.9 は慣用的な重症目安（記憶に基づく）' : '');
      if (b === 'high') add('Ca', 'warn', `イオン化Ca高値：${s.iCa} mmol/L`, '');
    } else if (isNum(d.correctedCa)) {
      const c = d.correctedCa;
      if (c > 14) add('Ca', 'crit', `重度高Ca血症：補正Ca ${r(c)}（>14）`, 'Minisola 2015（記憶に基づく）');
      else if (c >= 12) add('Ca', 'warn', `中等度高Ca血症：補正Ca ${r(c)}（12–14）`, 'Minisola 2015（記憶に基づく）');
      else if (c > R.tCa.hi) add('Ca', 'watch', `軽度高Ca血症：補正Ca ${r(c)}`, '');
      else if (c < R.tCa.lo) add('Ca', 'warn', `低Ca血症：補正Ca ${r(c)}（イオン化Caで確認推奨）`, 'Payne 補正は誤分類が多い');
    }
    if (isNum(s.iCa) && isNum(s.tCa)) {
      const ratio = s.tCa / s.iCa; // mg/dL ÷ mmol/L
      if (ratio > 10) add('Ca', 'watch', `総Ca/イオン化Ca 比 ${r(ratio)}（総Ca↑・iCa↓の乖離）：クエン酸蓄積を考慮`, '比の閾値は文献により異なる（要確認）');
    }
    if (isNum(s.Mg)) {
      const b = band(s.Mg, R.Mg);
      if (b === 'low') add('Mg', 'warn', `低Mg血症：Mg ${s.Mg}${isNum(s.K) && s.K < R.K.lo ? '（難治性低Kの原因となりうる）' : ''}`, '');
      if (b === 'high') add('Mg', 'watch', `高Mg血症：Mg ${s.Mg}`, '');
    }
    if (isNum(s.P)) {
      const b = band(s.P, R.P);
      if (b === 'low') add('P', s.P < 1.0 ? 'crit' : 'warn', `低P血症：P ${s.P}${s.P < 1.0 ? '（<1.0 重度）' : ''}`, s.P < 1.0 ? '慣用区分（記憶に基づく）' : '');
      if (b === 'high') add('P', 'watch', `高P血症：P ${s.P}`, '');
    }
    if (isNum(s.Lac) && s.Lac > R.Lac.hi) {
      add('Lac', s.Lac >= 4 ? 'crit' : 'warn', `高乳酸血症：${s.Lac} mmol/L`, s.Lac >= 4 ? 'Surviving Sepsis（≥4）' : '');
    }
    return f;
  }

  // ---------------------------------------------------------------- 単一時点の解析
  function analyze(s, userSettings) {
    const st = mergeSettings(userSettings);
    const d = { warnings: [] };

    if (has(s, 'pH', 'pCO2')) {
      d.hco3Calc = hhHCO3(s.pH, s.pCO2);
      if (isNum(s.HCO3)) {
        const diff = Math.abs(s.HCO3 - d.hco3Calc);
        if (diff > Math.max(2, 0.1 * d.hco3Calc)) {
          d.warnings.push(`入力 HCO₃⁻ ${s.HCO3} と pH/pCO₂ からの計算値 ${r(d.hco3Calc)} が乖離：転記・検体を確認`);
        }
      }
    }
    d.hco3 = isNum(s.HCO3) ? s.HCO3 : d.hco3Calc;

    if (has(s, 'Na', 'Cl') && isNum(d.hco3)) {
      d.ag = anionGap(s.Na, s.Cl, d.hco3);
      d.cag = correctedAG(d.ag, s.Alb);
    }
    if (has(s, 'Na', 'Cl')) d.naCorrCl = s.Cl * 140 / s.Na;
    if (isNum(s.Na)) d.correctedNa = isNum(s.Glu) ? correctedNa(s.Na, s.Glu, st.glucoseFactor) : s.Na;
    if (isNum(s.Na)) {
      d.calcOsm = calcOsm(s.Na, s.Glu, s.BUN);
      if (isNum(s.Osm) && has(s, 'Glu', 'BUN')) d.osmGap = s.Osm - d.calcOsm;
    }
    if (isNum(s.tCa)) d.correctedCa = correctedCa(s.tCa, s.Alb);

    if (has(s, 'pH', 'pCO2') && isNum(d.hco3)) {
      d.sbe = sbe(s.pH, d.hco3);
      d.acidBase = classifyAcidBase(s, d.hco3, d.cag, st);
      d.beParts = beDecomposition(s, d.hco3);
    }
    d.stewart = stewart(s);
    d.findings = electrolyteFindings(s, d, st);
    return d;
  }

  // ---------------------------------------------------------------- 時系列
  /** 対称RCV（%）= √2 × Z × √(CVA² + CVI²)（Fraser 2001） */
  function rcv(cvi, cva, z) {
    return Math.SQRT2 * z * Math.sqrt(cva * cva + cvi * cvi);
  }

  const TS_KEYS = ['Na', 'K', 'Cl', 'tCa', 'Mg', 'P'];

  /**
   * samples: [{ t: ISO文字列, ...値 }]
   * 返り値: 時刻順に並べた各時点の解析＋時系列アラート
   */
  function analyzeSeries(samples, userSettings) {
    const st = mergeSettings(userSettings);
    const rows = samples
      .map((s, i) => ({ s, i, ms: Date.parse(s.t) }))
      .filter((x) => Number.isFinite(x.ms))
      .sort((a, b) => a.ms - b.ms)
      .map((x) => ({ sample: x.s, ms: x.ms, index: x.i, result: analyze(x.s, st), deltas: {} }));
    const alerts = [];
    const naLimit = st.naHighRisk ? 8 : st.naLimit24h;

    for (let j = 1; j < rows.length; j++) {
      const cur = rows[j], prev = rows[j - 1];
      const hrs = (cur.ms - prev.ms) / 3.6e6;
      for (const k of TS_KEYS) {
        const a = prev.sample[k], b = cur.sample[k];
        if (!isNum(a) || !isNum(b) || a === 0) continue;
        const pct = (b - a) / a * 100;
        const bv = st.bv[k];
        const lim = bv ? rcv(bv.cvi, bv.cva, st.rcvZ) : null;
        cur.deltas[k] = { abs: b - a, pct, hrs, per24h: hrs > 0 ? (b - a) / hrs * 24 : null, rcv: lim, significant: lim != null && Math.abs(pct) > lim };
      }
    }

    // 低Na補正速度：直近24h窓内の最低値からの上昇量（窓内に低Na<135 がある場合のみ）
    // 高血糖では測定Na は血糖低下とともに上がるため、グルコース補正Na で評価する
    const naOf = (x) => x.result.correctedNa;
    for (let j = 0; j < rows.length; j++) {
      const cur = rows[j];
      if (!isNum(naOf(cur))) continue;
      const win = rows.slice(0, j + 1).filter((x) => isNum(naOf(x)) && cur.ms - x.ms <= H24);
      const nadir = win.reduce((m, x) => (naOf(x) < naOf(m) ? x : m), win[0]);
      if (naOf(nadir) >= 135 || nadir === cur) continue;
      const rise = naOf(cur) - naOf(nadir);
      const hrs = (cur.ms - nadir.ms) / 3.6e6;
      const corrected = win.some((x) => naOf(x) !== x.sample.Na);
      const lab = corrected ? '補正Na' : 'Na';
      cur.naRise24h = { rise, hrs, nadir: naOf(nadir), limit: naLimit, corrected };
      if (rise >= naLimit) {
        alerts.push({ ms: cur.ms, level: 'crit', key: 'Na', text: `低Na補正速度 超過：${r(hrs)}時間で ${lab} +${r(rise)} mmol/L（上限 ${naLimit}/24h）。過剰補正への対応（再低下）を検討`, src: st.naHighRisk ? 'Verbalis 2013（高リスク 8）' : 'ESE/ERA-EDTA 2014（10/24h）' });
      } else if (rise >= naLimit - 2) {
        alerts.push({ ms: cur.ms, level: 'warn', key: 'Na', text: `低Na補正速度 上限接近：${r(hrs)}時間で ${lab} +${r(rise)} mmol/L（上限 ${naLimit}/24h）`, src: '' });
      } else {
        // 直前区間の速度が続いた場合の24h窓での上昇見込み（線形外挿）
        const prev = rows[j - 1];
        if (prev && isNum(naOf(prev)) && cur.ms > prev.ms && hrs < 24) {
          const rate = (naOf(cur) - naOf(prev)) / ((cur.ms - prev.ms) / 3.6e6) * 24;
          const proj = rise + rate * (24 - hrs) / 24;
          if (rate > 0 && proj >= naLimit) {
            cur.naRise24h.projected = proj;
            alerts.push({ ms: cur.ms, level: 'watch', key: 'Na', text: `現在の上昇速度（${r(rate)}/24h）が続くと24h窓で ${lab} +${r(proj)} に達する見込み`, src: '線形外挿（参考）' });
          }
        }
      }
    }

    // 高Na血症の低下速度：表示のみ（成人の速度上限は議論中）
    for (const row of rows) {
      const d = row.deltas.Na;
      if (d && d.abs < 0 && isNum(row.sample.Na) && row.sample.Na - d.abs > 145 && d.per24h != null) {
        row.naFallRate = d.per24h;
      }
    }

    // K の急変（RCV超かつ 6時間以内に 1.0 mmol/L 以上）：DKA治療中の急落など
    for (const row of rows) {
      const d = row.deltas.K;
      if (d && d.significant && Math.abs(d.abs) >= 1.0 && d.hrs <= 6) {
        alerts.push({ ms: row.ms, level: 'warn', key: 'K', text: `K の急変：${r(d.hrs)}時間で ${d.abs > 0 ? '+' : ''}${r(d.abs)} mmol/L`, src: '閾値は本ツールの暫定値（根拠文献なし）' });
      }
    }

    // AG の閉鎖（DKA 等の治療反応）
    const agRows = rows.filter((x) => isNum(x.result.cag));
    if (agRows.length >= 2) {
      const first = agRows[0], last = agRows[agRows.length - 1];
      if (first.result.cag > st.agUpper && last.result.cag <= st.agUpper) {
        alerts.push({ ms: last.ms, level: 'ok', key: 'AG', text: `補正AG が ${r(first.result.cag)} → ${r(last.result.cag)} に低下（基準内）`, src: '' });
      }
    }

    // 最新時点の重大所見もアラートに載せる（過去時点の分は表と解析パネルで確認）
    const latest = rows[rows.length - 1];
    if (latest) for (const f of latest.result.findings) if (f.level === 'crit') alerts.push({ ms: latest.ms, level: 'crit', key: f.key, text: `最新値：${f.text}`, src: f.src });

    alerts.sort((a, b) => b.ms - a.ms);
    return { rows, alerts, settings: st };
  }

  return {
    DEFAULT_SETTINGS, mergeSettings, conv,
    hhHCO3, sbe, anionGap, correctedAG, correctedNa, correctedCa, calcOsm,
    stewart, beDecomposition, classifyAcidBase, electrolyteFindings,
    analyze, rcv, analyzeSeries, round: r,
  };
});
