// 実行: node --test electrolyte-monitor/core.test.js
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const C = require('./core.js');

const near = (a, b, tol = 0.05) => assert.ok(Math.abs(a - b) <= tol, `${a} ≉ ${b}`);

// 引き継ぎまとめ §7-1 の検証症例（Mg 0.5 mmol/L = 1.215 mg/dL, P 3.5 mg/dL）
const mvpCase = { pH: 7.32, pCO2: 32, HCO3: 16, Lac: 9, Alb: 4.0, Na: 140, K: 4.0, Cl: 104, iCa: 1.2, Mg: 1.215, P: 3.5 };

test('MVP症例: AG・補正AG・ΔΔ・Winter が引き継ぎ文書の値と一致', () => {
  const d = C.analyze(mvpCase);
  assert.equal(d.ag, 20);
  assert.equal(d.cag, 20);
  near(d.acidBase.deltaRatio, 1.0);
  near(d.acidBase.comp.expected, 32);
  assert.deepEqual(d.acidBase.primary, ['代謝性アシドーシス']);
  assert.equal(d.acidBase.secondary.length, 0, '単一AG開大性なので合併なし');
});

test('MVP症例: Stewart（SIDa 34.4 / Alb電荷 10.77 / SIDe 28.8 / SIG 5.6）', () => {
  const st = C.analyze(mvpCase).stewart;
  near(st.sidA, 34.4);
  near(st.albCharge, 10.77, 0.01);
  near(st.sidE, 28.8);
  near(st.sig, 5.6);
  assert.ok(st.sig > 0, 'SIG 正値＝未測定陰イオン（標準符号）');
});

test('Henderson–Hasselbalch: pH 7.40 / pCO2 40 → HCO3 ≈ 24.5', () => {
  near(C.hhHCO3(7.4, 40), 24.50, 0.01);
});

test('HCO3 入力と pH/pCO2 の乖離を警告する', () => {
  const d = C.analyze({ pH: 7.40, pCO2: 40, HCO3: 15, Na: 140, Cl: 100 });
  assert.equal(d.warnings.length, 1);
});

test('アルブミン補正AG（Alb 2.0 で +5）', () => {
  assert.equal(C.correctedAG(10, 2.0), 15);
});

test('ΔΔ < 0.8 で高Cl性代謝性アシドーシス併存を示す', () => {
  // AG 20 (ΔAG 8), HCO3 8 (ΔHCO3 16) → ΔΔ 0.5
  const d = C.analyze({ pH: 7.10, pCO2: 20, HCO3: 6, Na: 140, Cl: 114, Alb: 4.0 });
  assert.ok(d.acidBase.deltaRatio < 0.8);
  assert.ok(d.acidBase.secondary.some((x) => x.includes('高Cl')));
});

test('代償不十分な代謝性アシドーシスで呼吸性アシドーシス合併を示す', () => {
  const d = C.analyze({ pH: 7.05, pCO2: 40, HCO3: 11, Na: 140, Cl: 105, Alb: 4.0 });
  assert.ok(d.acidBase.secondary.some((x) => x.includes('呼吸性アシドーシス')));
});

test('pH正常でも AG 開大を検出する', () => {
  const d = C.analyze({ pH: 7.40, pCO2: 40, HCO3: 24, Na: 145, Cl: 100, Alb: 4.0 });
  assert.ok(d.cag > 16);
  assert.ok(d.acidBase.secondary.some((x) => x.includes('代謝性アルカローシス')));
});

test('補正Na: Glu 600 で Hillier 2.4 / Katz 1.6', () => {
  assert.equal(C.correctedNa(130, 600, 2.4), 142);
  assert.equal(C.correctedNa(130, 600, 1.6), 138);
});

test('浸透圧ギャップ', () => {
  const d = C.analyze({ Na: 140, Glu: 90, BUN: 14, Osm: 320 });
  near(d.calcOsm, 290);
  near(d.osmGap, 30);
  assert.ok(d.findings.some((f) => f.key === 'Osm'));
});

test('K の重症度区分（UKKA）', () => {
  const lv = (k) => C.analyze({ K: k }).findings.find((f) => f.key === 'K').level;
  assert.equal(lv(6.5), 'crit');
  assert.equal(lv(6.1), 'warn');
  assert.equal(lv(5.6), 'watch');
  assert.equal(lv(2.3), 'crit');
});

test('Na の重症度区分（ESE/ERA-EDTA 2014）', () => {
  const lv = (na) => C.analyze({ Na: na }).findings.find((f) => f.key === 'Na').level;
  assert.equal(lv(122), 'crit');
  assert.equal(lv(127), 'warn');
  assert.equal(lv(132), 'watch');
});

test('BE 分解の各成分の和は SBE に一致する', () => {
  const d = C.analyze(mvpCase);
  const p = d.beParts;
  near(p.naCl + p.alb + p.lac + p.other, p.sbe, 1e-9);
});

test('RCV: Na (CVI 0.5, CVA 0.25) ≈ 1.55%', () => {
  near(C.rcv(0.5, 0.25, 1.96), 1.55, 0.01);
});

test('時系列: 低Na補正 24h で +11 → 上限超過アラート', () => {
  const s = [
    { t: '2026-10-01T08:00', Na: 118 },
    { t: '2026-10-01T14:00', Na: 122 },
    { t: '2026-10-02T06:00', Na: 129 },
  ];
  const out = C.analyzeSeries(s);
  const a = out.alerts.find((x) => x.key === 'Na' && x.level === 'crit' && x.text.includes('超過'));
  assert.ok(a, '補正速度超過アラートが出る');
  assert.equal(out.rows[2].naRise24h.rise, 11);
});

test('時系列: 高リスク設定では上限 8', () => {
  const s = [{ t: '2026-10-01T08:00', Na: 118 }, { t: '2026-10-01T20:00', Na: 126 }];
  const out = C.analyzeSeries(s, { naHighRisk: true });
  assert.ok(out.alerts.some((x) => x.text.includes('上限 8')));
});

test('時系列: 24h 窓の外の最低値は使わない', () => {
  const s = [
    { t: '2026-10-01T08:00', Na: 118 },
    { t: '2026-10-02T09:00', Na: 126 },
    { t: '2026-10-02T20:00', Na: 130 },
  ];
  const out = C.analyzeSeries(s);
  assert.equal(out.rows[2].naRise24h.rise, 4);
});

test('時系列: RCV を超える変化に印が付く', () => {
  const out = C.analyzeSeries([{ t: '2026-10-01T00:00', K: 5.8 }, { t: '2026-10-01T04:00', K: 4.2 }]);
  const d = out.rows[1].deltas.K;
  assert.ok(d.significant);
  assert.ok(out.alerts.some((x) => x.key === 'K' && x.text.includes('急変')));
});

test('時系列: 入力順に依存せず時刻順に並べる', () => {
  const out = C.analyzeSeries([{ t: '2026-10-02T00:00', Na: 130 }, { t: '2026-10-01T00:00', Na: 125 }]);
  assert.equal(out.rows[0].sample.Na, 125);
});

test('時系列: DKA で血糖低下に伴う測定Na上昇は補正Na で評価し誤警告しない', () => {
  const s = [
    { t: '2026-10-05T02:00', Na: 131, Glu: 620 },
    { t: '2026-10-05T10:00', Na: 137, Glu: 260 },
    { t: '2026-10-06T02:00', Na: 139, Glu: 160 },
  ];
  const out = C.analyzeSeries(s);
  assert.ok(!out.alerts.some((x) => x.key === 'Na'));
});

test('時系列: 重大所見は最新時点のものだけをアラートに載せる', () => {
  const out = C.analyzeSeries([{ t: '2026-10-01T00:00', K: 6.8 }, { t: '2026-10-01T06:00', K: 6.6 }]);
  assert.equal(out.alerts.filter((x) => x.text.startsWith('最新値')).length, 1);
});
