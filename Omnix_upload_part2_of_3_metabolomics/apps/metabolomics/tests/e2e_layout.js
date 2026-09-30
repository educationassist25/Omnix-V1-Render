// No-overlap check of the network view across every "pathways shown" setting.
//   node tests/e2e_layout.js report1.html [report2.html ...] [--shots=dir]
// Each file is an interactive report written by rewiring_ui.render_view_html.
const { chromium } = require('playwright');
const path = require('path');
const overlapProblems = require('./overlap_check.js');
const files = process.argv.slice(2).filter(a => !a.startsWith('--'));
const shots = (process.argv.find(a => a.startsWith('--shots=')) || '').slice(8);
(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
  await page.route(/^https?:/, r => r.abort());   // offline: no web fonts, fallback fonts are measured too
  let bad = 0;
  for (const f of files) {
    await page.goto('file://' + path.resolve(f), { waitUntil: 'domcontentloaded' });
    await page.waitForFunction(() => window.__rewiringReady === true);
    const hasPw = await page.$('#pwTop');
    const settings = hasPw ? [] : [[null, null, null]];
    if (hasPw) for (const top of ['5', '10', '15', '20', '0']) for (const [met, cut] of [['p', '0.05'], ['q', '0.05'], ['q', '1'], ['p', '1']]) settings.push([top, met, cut]);
    for (const [top, met, cut] of settings) {
      if (top !== null) { await page.selectOption('#pwTop', top); await page.selectOption('#pwMetric', met); await page.selectOption('#pwCut', cut); }
      for (const selectNode of [false, true]) {
        if (selectNode) await page.evaluate(() => { const n = [...document.querySelectorAll('#net .node')]; if (n.length) n[Math.floor(n.length / 2)].dispatchEvent(new MouseEvent('click', { bubbles: true })); });
        const probs = await page.evaluate(`(${overlapProblems.toString()})()`);
        const n = await page.evaluate(() => document.querySelectorAll('#net .node').length);
        const trunc = await page.evaluate(() => [...document.querySelectorAll('#net text')].filter(t => t.textContent.includes('…')).length);
        const tag = `${path.basename(f)} top=${top} ${met}<=${cut}${selectNode ? ' +selected' : ''}`;
        if (probs.length || trunc) { bad++; console.log('FAIL', tag, n, 'nodes', trunc ? `${trunc} truncated` : '', probs.slice(0, 6)); }
        else console.log('ok  ', tag, n, 'nodes');
      }
      if (shots && top !== null && met === 'p' && cut === '0.05' && ['10', '0'].includes(top)) {
        const el = await page.$('#net'); await el.screenshot({ path: `${shots}/${path.basename(f, '.html')}_top${top}.png` });
      }
    }
  }
  console.log(bad ? `LAYOUT FAILED (${bad})` : 'LAYOUT PASSED');
  process.exitCode = bad ? 1 : 0;
  await browser.close();
})();
