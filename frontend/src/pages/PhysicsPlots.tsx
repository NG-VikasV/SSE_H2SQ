// Publication-quality physics plots — PRB / PRL figure style
// Okabe–Ito palette · Times New Roman · box frame · no grid · per-panel PNG save

import React, { useMemo, useEffect, useState, useRef } from 'react';
import {
  ComposedChart, Line, XAxis, YAxis,
  Tooltip, ResponsiveContainer, ReferenceLine, Label, Customized,
} from 'recharts';
import { Printer, BarChart2, Download, Settings } from 'lucide-react';

// ─── Okabe–Ito 7-color colorblind-safe palette ────────────────────────────────
const PALETTE = [
  '#0072B2', // blue
  '#D55E00', // vermillion
  '#009E73', // green
  '#CC79A7', // pink / magenta
  '#E69F00', // orange
  '#56B4E9', // sky blue
  '#000000', // black
];

// ─── Typographic constants ────────────────────────────────────────────────────
// "Times New Roman" is universally available and matches LaTeX / PRL output.
const RF = "'Times New Roman', Times, serif";

// ─── Chart geometry & theme ───────────────────────────────────────────────────
export interface ChartTheme {
  chartHeight: number;
  tickSize: number;
  tickFontSize: number;
  xAxisLabelOffset: number;
  yAxisLabelOffset: number;
  axisLabelFontSize: number;
  panelLabelSize: number;
  title: string;
  titleFontSize: number;
  titleX: number;
  titleY: number;
}

export const defaultTheme: ChartTheme = {
  chartHeight: 360,
  tickSize: 5,
  tickFontSize: 15,
  xAxisLabelOffset: -44,
  yAxisLabelOffset: 28,
  axisLabelFontSize: 18,
  panelLabelSize: 18,
  title: "H₂SQ Model — Physical Observables",
  titleFontSize: 20,
  titleX: 0,
  titleY: 0,
};

export const ChartThemeContext = React.createContext<ChartTheme>(defaultTheme);

const CHART_M = { top: 28, right: 32, bottom: 64, left: 85 };

// ─── Custom tick renderers (Outward Ticks) ────────────────────────────────────
function OutwardBottomTick(props: any) {
  const { x, y, payload } = props;
  const theme = React.useContext(ChartThemeContext);
  return (
    <g transform={`translate(${x},${y})`}>
      <line x1={0} y1={0} x2={0} y2={theme.tickSize} stroke="#1a1a1a" strokeWidth={1.2} />
      <text x={0} y={theme.tickSize + 19} textAnchor="middle" fill="#111" fontSize={theme.tickFontSize} fontFamily={RF} fontWeight={400}>
        {Number(payload.value).toFixed(1)}
      </text>
    </g>
  );
}

function OutwardLeftTick(props: any) {
  const { x, y, payload, yFmt } = props;
  const theme = React.useContext(ChartThemeContext);
  return (
    <g transform={`translate(${x},${y})`}>
      <line x1={0} y1={0} x2={-theme.tickSize} y2={0} stroke="#1a1a1a" strokeWidth={1.2} />
      <text x={-theme.tickSize - 8} y={5} textAnchor="end" fill="#111" fontSize={theme.tickFontSize} fontFamily={RF} fontWeight={400}>
        {yFmt(payload.value)}
      </text>
    </g>
  );
}

function OutwardTopTick(props: any) {
  const { x, y } = props;
  const theme = React.useContext(ChartThemeContext);
  return (
    <g transform={`translate(${x},${y})`}>
      <line x1={0} y1={0} x2={0} y2={-theme.tickSize} stroke="#1a1a1a" strokeWidth={1.2} />
    </g>
  );
}

function OutwardRightTick(props: any) {
  const { x, y } = props;
  const theme = React.useContext(ChartThemeContext);
  return (
    <g transform={`translate(${x},${y})`}>
      <line x1={0} y1={0} x2={theme.tickSize} y2={0} stroke="#1a1a1a" strokeWidth={1.2} />
    </g>
  );
}

// ─── Custom dot renderers ─────────────────────────────────────────────────────
function CircleDot(color: string) {
  return function Dot(props: any) {
    const { cx, cy } = props;
    if (cx == null || cy == null) return null;
    return (
      <circle cx={cx} cy={cy} r={5} fill="white" stroke={color} strokeWidth={2.2} />
    );
  };
}
function SquareDot(color: string) {
  return function Dot(props: any) {
    const { cx, cy } = props;
    if (cx == null || cy == null) return null;
    return <rect x={cx - 4.5} y={cy - 4.5} width={9} height={9} fill={color} />;
  };
}

// ─── PNG download — 3× scale for journal submission ──────────────────────────
function downloadAsPng(ref: React.RefObject<HTMLDivElement | null>, filename: string) {
  const wrapper = ref.current;
  if (!wrapper) return;
  const svg = wrapper.querySelector('svg');
  if (!svg) return;

  const { width, height } = svg.getBoundingClientRect();
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.setAttribute('width',  String(width));
  clone.setAttribute('height', String(height));
  clone.setAttribute('xmlns',  'http://www.w3.org/2000/svg');

  // Explicit white background
  const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  bg.setAttribute('x', '0'); bg.setAttribute('y', '0');
  bg.setAttribute('width', '100%'); bg.setAttribute('height', '100%');
  bg.setAttribute('fill', '#ffffff');
  clone.insertBefore(bg, clone.firstChild);

  // Frame + panel label are already in the SVG via the Customized component.
  // No manual rect injection needed here.

  const svgStr = new XMLSerializer().serializeToString(clone);
  const blob   = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
  const url    = URL.createObjectURL(blob);

  const DPI = 3; // ~300 DPI
  const canvas = document.createElement('canvas');
  canvas.width  = Math.round(width  * DPI);
  canvas.height = Math.round(height * DPI);
  const ctx = canvas.getContext('2d')!;
  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  const img = new Image();
  img.onload = () => {
    ctx.scale(DPI, DPI);
    ctx.drawImage(img, 0, 0);
    URL.revokeObjectURL(url);
    const a = document.createElement('a');
    a.href = canvas.toDataURL('image/png');
    a.download = filename + '.png';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  };
  img.onerror = () => URL.revokeObjectURL(url);
  img.src = url;
}

// ─── Physics data types ───────────────────────────────────────────────────────
interface PhysicsPoint {
  beta: number; N: number;
  ed_enrg?: number;  ed_enrg2?: number;
  ed_SMag2?: number; ed_SMag4?: number; ed_Magx?: number;
  ed_Cv?: number;    ed_Binder?: number;
  qmc_enrg?: number;  qmc_enrg2?: number;
  qmc_SMag2?: number; qmc_SMag4?: number; qmc_Magx?: number;
  qmc_rho_x?: number; qmc_rho_y?: number;
  qmc_Cv?: number; qmc_Binder?: number; qmc_rho?: number;
}
interface Group { label: string; systemLabel: string; pts: PhysicsPoint[]; }

// ─── Group results by (Lx,Ly,J0,J1,J2) ──────────────────────────────────────
function groupResults(results: any[]): Map<string, Group> {
  const map = new Map<string, Group>();
  for (const r of results) {
    const ed = r.ed as any, qmc = r.qmc as any, src = ed ?? qmc;
    if (!src) continue;
    const Lx = src.Lx ?? 4, Ly = src.Ly ?? 4, Nl = src.Nl ?? src.Lz ?? 1;
    const N = Lx * Ly * Nl, beta = Number(src.beta);
    const key         = `${Lx}x${Ly}_${src.J0}_${src.J1}_${src.J2}`;
    const label       = `L = ${Lx}×${Ly},  J₀ = ${src.J0},  J₁ = ${src.J1},  J₂ = ${src.J2}`;
    const systemLabel = `${Lx}×${Ly}`;
    if (!map.has(key)) map.set(key, { label, systemLabel, pts: [] });
    const grp = map.get(key)!;
    let pt = grp.pts.find(p => Math.abs(p.beta - beta) < 1e-9);
    if (!pt) { pt = { beta, N }; grp.pts.push(pt); }

    if (ed) {
      pt.ed_enrg  = ed.enrg;  pt.ed_enrg2 = ed.enrg2;
      pt.ed_SMag2 = ed.SMag_square; pt.ed_SMag4 = ed.SMag_four;
      pt.ed_Magx  = ed.Mag_x != null ? Math.abs(ed.Mag_x) : undefined;
      if (ed.enrg != null && ed.enrg2 != null)
        pt.ed_Cv = beta * beta * N * (ed.enrg2 - ed.enrg ** 2);
      if (ed.SMag_square != null && ed.SMag_four != null && ed.SMag_square > 1e-14)
        pt.ed_Binder = 1 - ed.SMag_four / (3 * ed.SMag_square ** 2);
    }
    if (qmc) {
      pt.qmc_enrg  = qmc.enrg;  pt.qmc_enrg2 = qmc.enrg2;
      pt.qmc_SMag2 = qmc.SMag_square; pt.qmc_SMag4 = qmc.SMag_four;
      pt.qmc_Magx  = qmc.Mag_x != null ? Math.abs(qmc.Mag_x) : undefined;
      pt.qmc_rho_x = qmc.rho_x; pt.qmc_rho_y = qmc.rho_y;
      if (qmc.enrg != null && qmc.enrg2 != null)
        pt.qmc_Cv = beta * beta * N * (qmc.enrg2 - qmc.enrg ** 2);
      if (qmc.SMag_square != null && qmc.SMag_four != null && qmc.SMag_square > 1e-14)
        pt.qmc_Binder = 1 - qmc.SMag_four / (3 * qmc.SMag_square ** 2);
      if (qmc.rho_x != null && qmc.rho_y != null)
        pt.qmc_rho = (qmc.rho_x + qmc.rho_y) / 2;
    }
  }
  for (const g of map.values()) g.pts.sort((a, b) => a.beta - b.beta);
  return map;
}

// ─── Flatten groups into recharts data array ──────────────────────────────────
function buildChartData(
  groups: Map<string, Group>,
  edKey: keyof PhysicsPoint | null,
  qmcKey: keyof PhysicsPoint | null,
) {
  const betaSet = new Set<number>();
  for (const { pts } of groups.values()) pts.forEach(p => betaSet.add(p.beta));
  return [...betaSet].sort((a, b) => a - b).map(beta => {
    const row: Record<string, number | null> = { beta };
    for (const [gkey, { pts }] of groups.entries()) {
      const pt = pts.find(p => Math.abs(p.beta - beta) < 1e-9);
      if (edKey)  row[`${gkey}|ed`]  = (pt && pt[edKey]  !== undefined) ? (pt[edKey]  as number) : null;
      if (qmcKey) row[`${gkey}|qmc`] = (pt && pt[qmcKey] !== undefined) ? (pt[qmcKey] as number) : null;
    }
    return row;
  });
}

// ─── Format Helpers ─────────────────────────────────────────────────────────────
function formatScientific(v: number, precision: number): string {
  if (v === 0) return '0';
  const str = v.toExponential(precision);
  const [base, expStr] = str.split('e');
  const exp = parseInt(expStr, 10);
  if (exp === 0) return base;
  
  const supMap: Record<string, string> = {
    '-': '⁻', '+': '⁺', '0': '⁰', '1': '¹', '2': '²', 
    '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'
  };
  const unicodeExp = exp.toString().split('').map(c => supMap[c]).join('');
  return `${base} × 10${unicodeExp}`;
}

// ─── Panel configuration ──────────────────────────────────────────────────────
interface PanelCfg {
  id: string; label: string;
  yAxis: string; yDesc: string;
  edKey: keyof PhysicsPoint | null; qmcKey: keyof PhysicsPoint | null;
  yFmt: (v: number) => string;
  refLine?: number; refLabel?: string;
  useScientific?: boolean;
}

const PANELS: PanelCfg[] = [
  {
    id: 'energy', label: '(a)',
    yAxis: '⟨E⟩ / N',            yDesc: 'Energy per site',
    edKey: 'ed_enrg',   qmcKey: 'qmc_enrg',   yFmt: v => v.toFixed(5),
  },
  {
    id: 'cv', label: '(b)',
    yAxis: 'Cᵥ / N',                  yDesc: 'Specific heat per site',
    edKey: 'ed_Cv',     qmcKey: 'qmc_Cv',     yFmt: v => v.toFixed(4),
  },
  {
    id: 'smag2', label: '(c)',
    yAxis: '⟨m²ᵣ⟩',    yDesc: 'Staggered magnetization squared',
    edKey: 'ed_SMag2',  qmcKey: 'qmc_SMag2',  yFmt: v => formatScientific(v, 2),
    useScientific: true,
  },
  {
    id: 'binder', label: '(d)',
    yAxis: 'Uᴸ',                       yDesc: 'Binder cumulant',
    edKey: 'ed_Binder', qmcKey: 'qmc_Binder', yFmt: v => v.toFixed(4),
    refLine: 2 / 3, refLabel: 'U* = 2/3',
  },
  {
    id: 'magx', label: '(e)',
    yAxis: '|⟨Mₓ⟩|',         yDesc: 'Transverse magnetization',
    edKey: 'ed_Magx',   qmcKey: 'qmc_Magx',   yFmt: v => formatScientific(v, 2),
    useScientific: true,
  },
  {
    id: 'rho', label: '(f)',
    yAxis: 'ρₛ',                  yDesc: 'Superfluid stiffness  (QMC only)',
    edKey: null,         qmcKey: 'qmc_rho',    yFmt: v => v.toFixed(5),
  },
];

// ─── Single panel component ───────────────────────────────────────────────────
function PlotPanel({ cfg, groups, colorIdx }: { cfg: PanelCfg, groups: Map<string, Group>, colorIdx: Map<string, number> }) {
  const theme = React.useContext(ChartThemeContext);
  const chartRef = useRef<HTMLDivElement>(null);
  const [hover, setHover] = useState(false);
  const data = useMemo(() => buildChartData(groups, cfg.edKey, cfg.qmcKey), [groups, cfg.edKey, cfg.qmcKey]);

  const maxVal = useMemo(() => {
    let max = 0;
    data.forEach(row => {
      Object.entries(row).forEach(([k, v]) => {
        if (v != null && (k.endsWith('|ed') || k.endsWith('|qmc'))) {
          if (Math.abs(v as number) > max) max = Math.abs(v as number);
        }
      });
    });
    return max;
  }, [data]);

  const exponent = (cfg.useScientific && maxVal > 0) ? Math.floor(Math.log10(maxVal)) : 0;
  const multiplier = Math.pow(10, -exponent);

  const hasData = data.some(row =>
    Object.keys(row).some(k => k !== 'beta' && row[k] != null),
  );



  return (
    <div
      style={{ background: '#fff', position: 'relative', paddingBottom: 4 }}
      onMouseEnter={() => setHover(true)}
      onMouseLeave={() => setHover(false)}
    >
      {/* ── Save PNG button (appears on hover) ── */}
      <button
        className="no-print"
        onClick={() => downloadAsPng(chartRef, `h2sq_${cfg.id}`)}
        style={{
          position: 'absolute', top: 8, right: 8, zIndex: 10,
          opacity: hover ? 1 : 0, transition: 'opacity 0.15s',
          background: '#111', color: '#fff', border: 'none',
          borderRadius: 3, padding: '4px 10px', fontSize: 10.5,
          fontFamily: RF, cursor: 'pointer',
          display: 'flex', alignItems: 'center', gap: 5,
          boxShadow: '0 1px 6px rgba(0,0,0,0.25)',
        }}
      >
        <Download size={10.5} /> Save PNG
      </button>

      {/* ── Chart — frame + panel letter rendered inside SVG via Customized ── */}
      <div ref={chartRef} style={{ position: 'relative', height: theme.chartHeight }}>
        <ResponsiveContainer width="100%" height="100%">
          <ComposedChart data={data} margin={CHART_M}>
            {/* NO CartesianGrid */}

            {/* Panel label drawn explicitly via Customized using standard recharts coords */}
            <Customized component={() => {
              const left = CHART_M.left;
              const top = CHART_M.top;
              const supMap: Record<string, string> = {
                '-': '⁻', '+': '⁺', '0': '⁰', '1': '¹', '2': '²', 
                '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'
              };
              const unicodeExp = exponent !== 0 ? exponent.toString().split('').map(c => supMap[c]).join('') : '';

              return (
                <g>
                  <text x={left + 9} y={top + 22} fontFamily={RF} fontSize={theme.panelLabelSize} fontWeight={700} fill="#111" textAnchor="start">
                    {cfg.label}
                  </text>
                  {cfg.useScientific && exponent !== 0 && (
                    <text x={left} y={top - 10} fontFamily={RF} fontSize={theme.tickFontSize} fontWeight={400} fill="#111" textAnchor="middle">
                      × 10{unicodeExp}
                    </text>
                  )}
                </g>
              );
            }} />

            {/* ── Top axis (box frame) ── */}
            <XAxis xAxisId="top" orientation="top" dataKey="beta" type="number" domain={['auto', 'auto']}
              padding={{ left: 15, right: 15 }} tickSize={0} tickMargin={0}
              tick={<OutwardTopTick />} tickLine={false} axisLine={{ stroke: '#1a1a1a', strokeWidth: 1.8 }} />

            {/* ── Bottom X axis ── */}
            <XAxis
              dataKey="beta" type="number" domain={['auto', 'auto']}
              padding={{ left: 15, right: 15 }} tickSize={0} tickMargin={0}
              tick={<OutwardBottomTick />}
              tickLine={false}
              axisLine={{ stroke: '#1a1a1a', strokeWidth: 1.8 }}
            >
              <Label
                value="β = 1/T"
                position="insideBottom" offset={theme.xAxisLabelOffset}
                style={{ fontFamily: RF, fontSize: theme.axisLabelFontSize, fill: '#111', fontWeight: 700, fontStyle: 'italic' }}
              />
            </XAxis>

            {/* ── Right Y axis (box frame) ── */}
            <YAxis yAxisId="right" orientation="right" domain={['auto', 'auto']} tickSize={0} tickMargin={0} tick={<OutwardRightTick />} tickLine={false} axisLine={{ stroke: '#1a1a1a', strokeWidth: 1.8 }} />

            {/* ── Left Y axis ── */}
            <YAxis
              domain={['auto', 'auto']} tickSize={0} tickMargin={0}
              tick={(props: any) => <OutwardLeftTick {...props} yFmt={(v: number) => cfg.useScientific ? (v * multiplier).toFixed(1) : cfg.yFmt(v)} />}
              tickLine={false}
              axisLine={{ stroke: '#1a1a1a', strokeWidth: 1.8 }}
              width={CHART_M.left - 2}
            >
              <Label
                value={cfg.yAxis}
                angle={-90} position="insideLeft" offset={theme.yAxisLabelOffset}
                style={{ fontFamily: RF, fontSize: theme.axisLabelFontSize, fill: '#111', fontWeight: 700, textAnchor: 'middle' }}
              />
            </YAxis>

            {/* ── Binder reference line ── */}
            {cfg.refLine !== undefined && (
              <ReferenceLine
                y={cfg.refLine}
                stroke="#999" strokeDasharray="5 3" strokeWidth={1.2}
                label={{
                  value: cfg.refLabel, position: 'insideTopRight',
                  fontSize: 10.5, fill: '#888', fontFamily: RF,
                }}
              />
            )}

            {/* ── Tooltip ── */}
            <Tooltip
              contentStyle={{
                background: '#fff', border: '1px solid #aaa', borderRadius: 2,
                fontSize: 11.5, fontFamily: RF, color: '#111',
                padding: '5px 10px', boxShadow: '0 2px 8px rgba(0,0,0,0.12)',
              }}
              itemStyle={{ color: '#333', fontSize: 11 }}
              labelFormatter={v => `β = ${Number(v).toFixed(4)}`}
              formatter={(v: any, name: any) => {
                const ns  = String(name ?? '');
                const isED  = ns.endsWith('|ed');
                const isQMC = ns.endsWith('|qmc');
                const gkey  = ns.replace(/\|(ed|qmc)$/, '');
                const grp   = groups.get(gkey);
                const tag   = isED ? ' (ED)' : isQMC ? ' (QMC)' : '';
                return [typeof v === 'number' ? cfg.yFmt(v) : v, (grp?.systemLabel ?? gkey) + tag];
              }}
            />

            {/* ── One ED + QMC pair per parameter group ── */}
            {[...groups.keys()].map(gkey => {
              const ci    = colorIdx.get(gkey) ?? 0;
              const color = PALETTE[ci % PALETTE.length];
              return (
                <React.Fragment key={gkey}>
                  {cfg.edKey && (
                    <Line
                      dataKey={`${gkey}|ed`} name={`${gkey}|ed`}
                      stroke={color} strokeWidth={1.8} strokeDasharray="8 4"
                      dot={CircleDot(color)}
                      activeDot={{ r: 6, fill: 'white', stroke: color, strokeWidth: 2.2 }}
                      connectNulls={false} type="linear"
                      isAnimationActive={false} legendType="none"
                    />
                  )}
                  {cfg.qmcKey && (
                    <Line
                      dataKey={`${gkey}|qmc`} name={`${gkey}|qmc`}
                      stroke={color} strokeWidth={1.8}
                      dot={SquareDot(color)}
                      activeDot={{ r: 6, fill: color }}
                      connectNulls={false} type="linear"
                      isAnimationActive={false} legendType="none"
                    />
                  )}
                </React.Fragment>
              );
            })}
          </ComposedChart>
        </ResponsiveContainer>

        {/* Frame is drawn inside the SVG via Customized — no HTML overlay needed */}

        {/* No-data message */}
        {!hasData && (
          <div style={{
            position: 'absolute', inset: 0, display: 'flex',
            alignItems: 'center', justifyContent: 'center',
            fontSize: 11, color: '#bbb', fontFamily: RF, fontStyle: 'italic',
          }}>
            no data
          </div>
        )}
      </div>

      {/* Panel descriptor */}
      <div style={{
        textAlign: 'center', marginTop: 3,
        fontSize: 10.5, color: '#555', fontFamily: RF, fontStyle: 'italic',
      }}>
        {cfg.yDesc}
      </div>
    </div>
  );
}

// ─── Legend SVG atoms ─────────────────────────────────────────────────────────
function LegendLine({ color }: { color: string }) {
  return (
    <svg width="34" height="14" style={{ display: 'inline', verticalAlign: 'middle' }}>
      <line x1="2" y1="7" x2="32" y2="7" stroke={color} strokeWidth={1.8} strokeDasharray="8 4" />
      <circle cx="17" cy="7" r="5" fill="white" stroke={color} strokeWidth={2.2} />
    </svg>
  );
}
function LegendSq({ color }: { color: string }) {
  return (
    <svg width="34" height="14" style={{ display: 'inline', verticalAlign: 'middle' }}>
      <line x1="2" y1="7" x2="32" y2="7" stroke={color} strokeWidth={1.8} />
      <rect x="12.5" y="2.5" width="9" height="9" fill={color} />
    </svg>
  );
}

// ─── Main page ────────────────────────────────────────────────────────────────
export default function PhysicsPlots() {
  const [results, setResults] = useState<any[]>([]);
  const [theme, setTheme] = useState<ChartTheme>(defaultTheme);
  const [showSettings, setShowSettings] = useState(false);

  const refresh = () => {
    try { setResults(JSON.parse(localStorage.getItem('results') ?? '[]')); }
    catch { setResults([]); }
  };
  useEffect(refresh, []);

  const groups   = useMemo(() => groupResults(results), [results]);
  const colorIdx = useMemo(() => {
    const m = new Map<string, number>(); let i = 0;
    for (const k of groups.keys()) m.set(k, i++);
    return m;
  }, [groups]);

  const hasData   = groups.size > 0;
  const paramSets = useMemo(() => {
    const seen = new Set<string>();
    return results.map(r => {
      const s = r.ed ?? r.qmc; if (!s) return null;
      const k = `${s.Lx}|${s.J0}|${s.J1}|${s.J2}`;
      if (seen.has(k)) return null; seen.add(k); return s;
    }).filter(Boolean) as any[];
  }, [results]);

  return (
    // Outer wrapper creates an isolated stacking context so the dark App
    // glow-bg elements cannot bleed through the physics-plots white background.
    <div style={{
      position: 'relative', zIndex: 1, isolation: 'isolate',
      background: '#ebebeb', minHeight: '100vh',
    }}>
      <style>{`
        @media print {
          .no-print { display: none !important; }
          body, html { background: #fff !important; }
        }
      `}</style>

      {/* ── Toolbar ─────────────────────────────────────────────────────────── */}
      <div
        className="no-print"
        style={{
          background: '#f5f5f5', borderBottom: '1px solid #d5d5d5',
          padding: '10px 28px',
          display: 'flex', alignItems: 'center', justifyContent: 'space-between',
          position: 'sticky', top: 0, zIndex: 20,
        }}
      >
        <div>
          <div style={{ fontWeight: 700, fontSize: 13.5, color: '#111', fontFamily: RF }}>
            Physics Plots
          </div>
          <div style={{ fontSize: 10.5, color: '#888', marginTop: 2, fontFamily: RF, fontStyle: 'italic' }}>
            Hover any panel → Save PNG (3× DPI) · Okabe–Ito palette · box-frame axes
          </div>
        </div>
        <div style={{ display: 'flex', gap: 8 }}>
          <button onClick={() => setShowSettings(s => !s)} style={{
            background: showSettings ? '#111' : '#eee', 
            color: showSettings ? '#fff' : '#444', 
            border: showSettings ? 'none' : '1px solid #ccc', 
            borderRadius: 4, padding: '5px 12px', fontSize: 11, fontFamily: RF,
            cursor: 'pointer', display: 'flex', alignItems: 'center', gap: 5,
          }}>
            <Settings size={12} /> {showSettings ? 'Close Settings' : 'Style Settings'}
          </button>
          <button onClick={refresh} style={{
            background: '#eee', border: '1px solid #ccc', borderRadius: 4,
            padding: '5px 12px', fontSize: 11, fontFamily: RF,
            cursor: 'pointer', color: '#444',
            display: 'flex', alignItems: 'center', gap: 5,
          }}>
            <BarChart2 size={12} /> Refresh
          </button>
          <button onClick={() => window.print()} style={{
            background: '#111', color: '#fff', border: 'none', borderRadius: 4,
            padding: '5px 14px', fontSize: 11, fontFamily: RF,
            cursor: 'pointer', display: 'flex', alignItems: 'center', gap: 5,
          }}>
            <Printer size={12} /> Print all
          </button>
        </div>
      </div>

      {/* ── Style Settings Panel ── */}
      {showSettings && (
        <div className="no-print" style={{ background: '#fff', borderBottom: '1px solid #ccc', padding: '16px 28px', display: 'flex', gap: 32, flexWrap: 'wrap', fontFamily: 'sans-serif', fontSize: 12 }}>
          <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }}>
            <h4 style={{ margin: 0, fontSize: 13, fontWeight: 600, color: '#111' }}>Plot Geometry</h4>
            <label style={{ display: 'flex', justifyContent: 'space-between', width: 200, color: '#444' }}>
              Height: <input type="number" value={theme.chartHeight} onChange={e => setTheme({...theme, chartHeight: +e.target.value})} style={{ width: 60 }} />
            </label>
            <label style={{ display: 'flex', justifyContent: 'space-between', width: 200, color: '#444' }}>
              X Label Offset: <input type="number" value={theme.xAxisLabelOffset} onChange={e => setTheme({...theme, xAxisLabelOffset: +e.target.value})} style={{ width: 60 }} />
            </label>
            <label style={{ display: 'flex', justifyContent: 'space-between', width: 200, color: '#444' }}>
              Y Label Offset: <input type="number" value={theme.yAxisLabelOffset} onChange={e => setTheme({...theme, yAxisLabelOffset: +e.target.value})} style={{ width: 60 }} />
            </label>
          </div>
          
          <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }}>
            <h4 style={{ margin: 0, fontSize: 13, fontWeight: 600, color: '#111' }}>Ticks & Fonts</h4>
            <label style={{ display: 'flex', justifyContent: 'space-between', width: 200, color: '#444' }}>
              Tick Length: <input type="number" value={theme.tickSize} onChange={e => setTheme({...theme, tickSize: +e.target.value})} style={{ width: 60 }} />
            </label>
            <label style={{ display: 'flex', justifyContent: 'space-between', width: 200, color: '#444' }}>
              Tick Font Size: <input type="number" value={theme.tickFontSize} onChange={e => setTheme({...theme, tickFontSize: +e.target.value})} style={{ width: 60 }} />
            </label>
            <label style={{ display: 'flex', justifyContent: 'space-between', width: 200, color: '#444' }}>
              Axis Label Size: <input type="number" value={theme.axisLabelFontSize} onChange={e => setTheme({...theme, axisLabelFontSize: +e.target.value})} style={{ width: 60 }} />
            </label>
            <label style={{ display: 'flex', justifyContent: 'space-between', width: 200, color: '#444' }}>
              Panel Label Size: <input type="number" value={theme.panelLabelSize} onChange={e => setTheme({...theme, panelLabelSize: +e.target.value})} style={{ width: 60 }} />
            </label>
          </div>

          <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }}>
            <h4 style={{ margin: 0, fontSize: 13, fontWeight: 600, color: '#111' }}>Global Title</h4>
            <label style={{ display: 'flex', flexDirection: 'column', gap: 4, width: 200, color: '#444' }}>
              Title Text: 
              <input type="text" value={theme.title} onChange={e => setTheme({...theme, title: e.target.value})} style={{ width: '100%', boxSizing: 'border-box' }} />
            </label>
            <label style={{ display: 'flex', justifyContent: 'space-between', width: 200, color: '#444' }}>
              Title Font Size: <input type="number" value={theme.titleFontSize} onChange={e => setTheme({...theme, titleFontSize: +e.target.value})} style={{ width: 60 }} />
            </label>
            <div style={{ display: 'flex', gap: 8 }}>
              <label style={{ display: 'flex', flexDirection: 'column', width: 96, color: '#444' }}>
                X Offset: <input type="number" value={theme.titleX} onChange={e => setTheme({...theme, titleX: +e.target.value})} />
              </label>
              <label style={{ display: 'flex', flexDirection: 'column', width: 96, color: '#444' }}>
                Y Offset: <input type="number" value={theme.titleY} onChange={e => setTheme({...theme, titleY: +e.target.value})} />
              </label>
            </div>
          </div>
        </div>
      )}

      {/* ── Figure area ──────────────────────────────────────────────────────── */}
      <ChartThemeContext.Provider value={theme}>
        <div style={{ padding: 32 }}>
        <div style={{
          maxWidth: 1020, margin: '0 auto',
          background: '#fff', border: '1px solid #bbb',
          padding: '28px 28px 22px',
          boxShadow: '0 2px 20px rgba(0,0,0,0.09)',
        }}>

          {/* Header */}
          <div style={{
            marginBottom: 18, paddingBottom: 14,
            borderBottom: '1.5px solid #ddd',
          }}>
            <div style={{
              fontSize: theme.titleFontSize, fontWeight: 700, color: '#111',
              fontFamily: RF, letterSpacing: '-0.2px',
              transform: `translate(${theme.titleX}px, ${theme.titleY}px)`,
            }}>
              {theme.title}
            </div>

            {hasData && (
              <div style={{
                marginTop: 11,
                display: 'flex', flexWrap: 'wrap', gap: '10px 22px', alignItems: 'center',
              }}>
                {[...groups.entries()].map(([gkey, { label }]) => {
                  const ci    = colorIdx.get(gkey) ?? 0;
                  const color = PALETTE[ci % PALETTE.length];
                  return (
                    <span key={gkey} style={{
                      display: 'flex', alignItems: 'center', gap: 6,
                      fontSize: 16, color: '#222', fontFamily: RF,
                    }}>
                      <span style={{
                        display: 'inline-block', width: 11, height: 11,
                        background: color, borderRadius: 1,
                      }} />
                      {label}
                    </span>
                  );
                })}
                {/* Symbol key */}
                <span style={{
                  display: 'flex', alignItems: 'center', gap: 16,
                  fontSize: 14, color: '#555', fontFamily: RF,
                  marginLeft: 10, paddingLeft: 14, borderLeft: '1px solid #ccc',
                }}>
                  <span style={{ display: 'flex', alignItems: 'center', gap: 5 }}>
                    <LegendLine color="#555" /> ED
                  </span>
                  <span style={{ display: 'flex', alignItems: 'center', gap: 5 }}>
                    <LegendSq color="#555" /> QMC (SSE)
                  </span>
                </span>
              </div>
            )}
          </div>

          {/* ── Panels ── */}
          {!hasData ? (
            <div style={{
              height: 260, display: 'flex', flexDirection: 'column',
              alignItems: 'center', justifyContent: 'center',
              color: '#bbb', fontFamily: RF, gap: 9,
            }}>
              <BarChart2 size={40} />
              <div style={{ fontSize: 16 }}>No simulation data available.</div>
              <div style={{ fontSize: 14, fontStyle: 'italic' }}>
                Run simulations from the Simulation tab, then return here.
              </div>
            </div>
          ) : (
            <>
              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 22 }}>
                {PANELS.map((cfg, i) => (
                  <PlotPanel key={i} cfg={cfg} groups={groups} colorIdx={colorIdx} />
                ))}
              </div>

              {/* ── Figure caption ── */}
              <div style={{
                marginTop: 22, paddingTop: 14, borderTop: '1px solid #ddd',
                fontSize: 14, lineHeight: 1.72, color: '#444',
                fontFamily: 'Georgia, "Times New Roman", serif', fontStyle: 'italic',
              }}>
                <strong style={{ fontStyle: 'normal' }}>FIG. 1.</strong>{' '}
                Physical observables of the H₂SQ model computed by stochastic series expansion
                quantum Monte Carlo (QMC-SSE, filled squares, solid lines) and exact diagonalization
                (ED, open circles, dashed lines) as functions of inverse temperature β = 1/T.
                (a) Energy per site ⟨E⟩/N.
                {' '}(b) Specific heat per site C<sub style={{ fontStyle: 'normal' }}>V</sub>/N
                  {' '}= β²N[⟨(E/N)²⟩ − ⟨E/N⟩²].
                {' '}(c) Staggered magnetization squared ⟨m²<sub style={{ fontStyle: 'normal' }}>z</sub>⟩.
                {' '}(d) Binder cumulant U<sub style={{ fontStyle: 'normal' }}>L</sub>
                  {' '}= 1 − ⟨m⁴<sub style={{ fontStyle: 'normal' }}>z</sub>⟩/(3⟨m²<sub style={{ fontStyle: 'normal' }}>z</sub>⟩²);
                  {' '}dashed line marks the 3D-Ising critical fixed point U* = 2/3.
                {' '}(e) Transverse magnetization |⟨M<sub style={{ fontStyle: 'normal' }}>x</sub>⟩|.
                {' '}(f) Superfluid stiffness ρ<sub style={{ fontStyle: 'normal' }}>s</sub>
                  {' '}= (ρ<sub style={{ fontStyle: 'normal' }}>x</sub>
                  + ρ<sub style={{ fontStyle: 'normal' }}>y</sub>)/2 (QMC only).
                {paramSets.length > 0 && (
                  <>
                    {' '}Parameters:
                    {' '}{paramSets.map((p, i) => (
                      <span key={i}>
                        {i > 0 && ';  '}
                        L = {p.Lx}×{p.Ly}, J₀ = {p.J0}, J₁ = {p.J1}, J₂ = {p.J2}
                      </span>
                    ))}.
                  </>
                )}
              </div>
            </>
          )}
        </div>
        </div>
      </ChartThemeContext.Provider>
    </div>
  );
}
