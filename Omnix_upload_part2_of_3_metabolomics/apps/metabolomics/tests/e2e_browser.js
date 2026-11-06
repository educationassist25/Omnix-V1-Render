// End-to-end browser check of the Rewiring Map inside a running Streamlit app.
//   streamlit run rewiring_app.py --server.headless true --server.port 8599
//   node tests/e2e_browser.js http://localhost:8599 [screenshot_dir]
const { chromium } = require('playwright');
const url = process.argv[2] || 'http://localhost:8599';
const out = (process.argv[3] && !process.argv[3].startsWith('--')) ? process.argv[3] : '.';
const fs = require('fs');
const overlapProblems = require('./overlap_check.js');
const fail = (m) => { console.error('FAIL: ' + m); process.exitCode = 1; };
(async () => {
  const browser = await chromium.launch();
  const context = await browser.newContext({ acceptDownloads: true, viewport: { width: 1440, height: 1000 } });
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('page: ' + e.message));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource|fonts\.g|metrics config/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto(url, { waitUntil: 'networkidle' });
  if (process.argv.includes('--main')) {
    // main MetaboAI Pro app: open Network Analysis and switch on the simulated demo
    await page.getByRole('button', { name: 'Network Analysis' }).click();
    await page.getByText('Explore with the simulated rewiring demo').click();
  }
  await page.getByRole('button', { name: /Run Rewiring Analysis/ }).click();
  await page.waitForSelector('iframe', { timeout: 60000 });
  // the Streamlit iframe that holds our view
  let frame = null;
  for (let t = 0; t < 60 && !frame; t++) {
    for (const f of page.frames()) { try { if (await f.evaluate(() => window.__rewiringReady === true)) frame = f; } catch (_) {} }
    if (!frame) await page.waitForTimeout(500);
  }
  if (!frame) { fail('view never became ready'); await browser.close(); return; }
  // default network: top 10 pathways by p-value <= 0.05, drawn without any overlap
  const def = await frame.evaluate(() => ({ note: document.getElementById('netNote').innerText, top: document.getElementById('pwTop').value,
    metric: document.getElementById('pwMetric').value, cut: document.getElementById('pwCut').value }));
  console.log('default network', JSON.stringify(def));
  if (!(def.top === '10' && def.metric === 'p' && def.cut === '0.05' && /Showing the top 10 pathways by p-value ≤ 0\.05/.test(def.note))) fail('default should be top 10 by p-value ≤ 0.05');
  for (const top of ['10', '0']) {
    await frame.selectOption('#pwTop', top);
    const probs = await frame.evaluate(`(${overlapProblems.toString()})()`);
    console.log(`overlap check top=${top}:`, probs.length ? probs.slice(0, 5) : 'clean');
    if (probs.length) fail('text overlaps a line, dot or other text');
  }
  await frame.selectOption('#pwTop', '0');   // all pathways for the structural checks below
  await page.waitForTimeout(200);
  const info = await frame.evaluate(() => ({
    edges: document.querySelectorAll('#net .hit:not(.rxhit)').length, rxEdges: document.querySelectorAll('#net .rxe').length, pwRows: document.querySelectorAll('#ptbl tbody tr').length,
    nodes: document.querySelectorAll('#net .node').length,
    cards: document.querySelectorAll('.card').length,
    rows: document.querySelectorAll('#tbl tbody tr').length,
    stats: document.getElementById('stats').innerText.replace(/\s+/g, ' '),
    detail: document.querySelector('#detail h2').innerText,
    verdicts: [...document.querySelectorAll('.verdict')].map(v => v.innerText),
    firstCard: document.querySelector('.card h3').innerText,
  }));
  console.log(JSON.stringify(info, null, 1));
  if (info.nodes !== 43) fail('expected 43 nodes, got ' + info.nodes);
  if (info.edges < 20) fail('too few rewired edges drawn: ' + info.edges);
  if (info.rows !== info.edges) fail(`table rows (${info.rows}) != rewired edges (${info.edges})`);
  if (info.cards < 3) fail('expected hypothesis cards');
  if (!info.verdicts.some(v => /Fragile|Mixed/.test(v))) fail('a flagged card should be shown');
  if (info.rxEdges < 20) fail('reaction layer should draw the measured reaction pairs: ' + info.rxEdges);
  if (info.pwRows < 5) fail('pathway-level table should list pathways: ' + info.pwRows);
  // pathway row click highlights its metabolites
  await (await frame.$('#ptbl tbody tr')).click();
  await page.waitForTimeout(200);
  const pwSel = await frame.evaluate(() => ({ active: document.querySelectorAll('#ptbl tr.active').length, dimNodes: document.querySelectorAll('#net .node.dim').length }));
  console.log('after pathway click', pwSel);
  if (pwSel.active !== 1 || pwSel.dimNodes === 0) fail('pathway row should highlight its metabolites');
  await (await frame.$('#ptbl tbody tr')).click();   // toggle off
  // reaction-only edge opens the reaction panel with a product/substrate ratio
  await frame.evaluate(() => document.querySelector('#net .rxhit').dispatchEvent(new MouseEvent('click', { bubbles: true })));
  await page.waitForTimeout(200);
  const rxPanel = await frame.evaluate(() => { const b = document.querySelector('#detail .rx-box'); return b ? b.innerText : ''; });
  console.log('reaction panel:', rxPanel.replace(/\s+/g, ' ').slice(0, 160));
  if (!/reaction|steps/i.test(rxPanel) || !/log2\(/.test(rxPanel)) fail('reaction edge should open the reaction panel');
  // iframe auto-height: the whole view should be visible without an inner scrollbar
  const heights = await page.evaluate(() => [...document.querySelectorAll('iframe')].map(f => f.getBoundingClientRect().height));
  console.log('iframe heights', heights);
  if (!heights.some(h => h > 1200)) fail('iframe did not size to content');
  // interactions
  const cards = await frame.$$('.card');
  await cards[1].click();
  await page.waitForTimeout(300);
  const afterCard = await frame.evaluate(() => ({ active: document.querySelectorAll('.card.active').length, dim: document.querySelectorAll('#net .edge.dim').length, detail: document.querySelector('#detail h2').innerText }));
  console.log('after card click', afterCard);
  if (afterCard.active !== 1 || afterCard.dim === 0) fail('card click should highlight a module');
  const lastRow = await frame.$('#tbl tbody tr:last-child');
  const rowText = await lastRow.evaluate(r => r.cells[0].innerText);
  await lastRow.click();
  await page.waitForTimeout(300);
  const detail = await frame.evaluate(() => document.querySelector('#detail h2').innerText);
  if (detail !== rowText) fail(`row click should select pair: ${detail} vs ${rowText}`);
  const chip = await frame.$('.chip[data-k="lost"]');
  const before = await frame.evaluate(() => document.querySelectorAll('#net .hit:not(.rxhit)').length);
  await chip.click();
  const after = await frame.evaluate(() => document.querySelectorAll('#net .hit:not(.rxhit)').length);
  console.log('edges before/after hiding lost', before, after);
  if (!(after < before)) fail('class chip should hide edges');
  await (await frame.$('.chip[data-k="lost"]')).click();   // chips re-render, so query again
  // hover tooltip on an edge
  // a curved edge's bounding-box centre is usually off the stroke, so fire the event directly
  await frame.evaluate(() => document.querySelector('#net .hit').dispatchEvent(new MouseEvent('mousemove', { clientX: 200, clientY: 200, bubbles: true })));
  const tip = await frame.evaluate(() => !document.getElementById('tip').hidden && document.getElementById('tip').innerText);
  if (!tip || !/↔/.test(tip)) fail('edge hover should show a tooltip');
  // network filter controls: top 5 shows fewer metabolites than all; 'no pathway' cutoff empties the ring
  const nodeCount = () => frame.evaluate(() => document.querySelectorAll('#net .node').length);
  await frame.selectOption('#pwTop', '0'); const nAll = await nodeCount();
  await frame.selectOption('#pwTop', '5'); const n5 = await nodeCount();
  await frame.selectOption('#pwMetric', 'p'); await frame.selectOption('#pwCut', '0.001'); const n0 = await nodeCount();
  console.log('nodes all/top5/none', nAll, n5, n0);
  // full names, nothing clipped, download buttons below the figures
  await frame.selectOption('#pwCut', '1'); await frame.selectOption('#pwTop', '0'); await page.waitForTimeout(300);
  const clip = await frame.evaluate(() => { const svg = document.getElementById('net'); const sb = svg.getBoundingClientRect();
    const out = []; svg.querySelectorAll('text').forEach(t => { const b = t.getBoundingClientRect();
      if (b.left < sb.left - 1 || b.right > sb.right + 1 || b.top < sb.top - 1 || b.bottom > sb.bottom + 1) out.push(t.textContent); });
    const netBtn = document.querySelector('.fig-foot [data-dl="net"]'); const scBtn = document.querySelector('[data-dl="scatter"]');
    return { clipped: out, truncated: [...svg.querySelectorAll('text')].filter(t => t.textContent.includes('…')).length,
      netBelow: !!netBtn && netBtn.getBoundingClientRect().top > sb.bottom - 2 && netBtn.getBoundingClientRect().top - sb.bottom < 40,
      scatterBelow: !!scBtn && scBtn.getBoundingClientRect().top > document.getElementById('scatter').getBoundingClientRect().bottom - 2
        && scBtn.getBoundingClientRect().top - document.getElementById('scatter').getBoundingClientRect().bottom < 40,
      formats: ['net', 'scatter'].map(w => [...document.querySelectorAll(`button[data-dl="${w}"]`)].map(b => b.dataset.fmt).join(',')),
      scatterClipped: (() => { const sv = document.getElementById('scatter'); const r = sv.getBoundingClientRect(); return [...sv.querySelectorAll('text')].filter(t => { const b = t.getBoundingClientRect(); return b.left < r.left - 1 || b.right > r.right + 1 || b.top < r.top - 1 || b.bottom > r.bottom + 1; }).map(t => t.textContent); })(),
      bareQ: /\bq\b(?!-value)/.test([...document.querySelectorAll('th,.kv .k,#net text')].map(e => e.textContent).join(' ')) }; });
  console.log('labels', JSON.stringify(clip));
  if (clip.clipped.length || clip.truncated || !clip.netBelow || !clip.scatterBelow || clip.bareQ || clip.scatterClipped.length) fail('labels / buttons / wording');
  if (clip.formats.some(f => f !== 'svg,png,jpeg,tiff,pdf')) fail('download buttons for SVG, PNG, JPEG, TIFF, PDF under each image');
  if (!(nAll === 43 && n5 > 0 && n5 < nAll && n0 === 0)) fail('pathway filter controls');
  await frame.selectOption('#pwCut', '1'); await frame.selectOption('#pwTop', '10');
  if ((await nodeCount()) === 0 || (await frame.evaluate(() => document.querySelectorAll('#net .pw-label tspan').length)) === 0) fail('wrapped pathway labels');
  const axis = await frame.evaluate(() => [...document.querySelectorAll('#scatter text')].map(t => t.textContent).filter(t => /Relative Abundance/.test(t)).length);
  if (axis !== 2) fail('correlation plot axes should read Relative Abundance (...)');
  // the buttons right under each image download from inside the Streamlit iframe (sandbox allows downloads)
  const magic = { svg: '<svg', png: '\x89PNG', jpeg: '\xff\xd8', tiff: 'II*\x00', pdf: '%PDF' };
  for (const [which, fmt] of [['net', 'svg'], ['net', 'png'], ['net', 'tiff'], ['net', 'pdf'], ['scatter', 'jpeg'], ['scatter', 'pdf']]) {
    try {
      const [dl] = await Promise.all([page.waitForEvent('download', { timeout: 15000 }), frame.click(`button[data-dl="${which}"][data-fmt="${fmt}"]`)]);
      const buf = fs.readFileSync(await dl.path()); const head = buf.slice(0, 5).toString('latin1');
      const ok = fmt === 'svg' ? buf.toString('utf8').includes('<svg') : head.startsWith(magic[fmt]);
      console.log(`in-view download ${which}.${fmt}:`, dl.suggestedFilename(), buf.length, 'bytes', ok ? 'valid' : 'INVALID');
      if (!ok || buf.length < 1000) fail(`in-view ${which} ${fmt} download invalid`);
    } catch (e) { fail(`in-view ${which} ${fmt} download did not start: ${e.message.slice(0, 80)}`); }
  }
  // the legend is part of the network image (and so of every download)
  const leg = await frame.evaluate(() => [...document.querySelectorAll('#net .net-legend text')].map(t => t.textContent));
  console.log('network legend:', leg.join(' | '));
  if (!['Coupling lost', 'Coupling gained', 'Sign flipped', 'Reaction (direct)', 'Reaction (2 steps)'].every(x => leg.includes(x)) || !leg.some(x => /Node colour: log2 FC/.test(x))) fail('network legend missing');
  if (await page.getByText('More figure formats').count()) fail('the separate figure-format section should be gone');
  await page.screenshot({ path: out + (process.argv.includes('--main') ? '/e2e_main.png' : '/e2e_standalone.png'), fullPage: true });
  const el = await page.$('iframe'); await el.screenshot({ path: out + (process.argv.includes('--main') ? '/e2e_main_view.png' : '/e2e_view.png') });
  if (errors.length) fail('JS errors: ' + errors.join(' | '));
  console.log(process.exitCode ? 'E2E FAILED' : 'E2E PASSED');
  await browser.close();
})();
