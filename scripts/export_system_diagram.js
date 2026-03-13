#!/usr/bin/env node
/**
 * Export ELM system diagram PNGs using Cytoscape.js + Puppeteer.
 *
 * Usage:
 *   node scripts/export_system_diagram.js              # both versions
 *   node scripts/export_system_diagram.js --paper      # paper version only (no cana→DNA edge)
 *   node scripts/export_system_diagram.js --deck       # deck version only (all edges)
 *
 * Outputs to docs/figures/:
 *   system_diagram_paper1.png   — Paper 1 Figure 1 (no canagliflozin antioxidant edge)
 *   system_diagram_deck.png     — Full deck slide 9 version (all edges)
 *
 * Requires: npm install puppeteer  (or npx will auto-fetch)
 */

const path = require('path');
const fs = require('fs');

const sleep = ms => new Promise(r => setTimeout(r, ms));

async function main() {
    const args = process.argv.slice(2);
    const doPaper = args.length === 0 || args.includes('--paper');
    const doDeck  = args.length === 0 || args.includes('--deck');

    const htmlPath = path.resolve(__dirname, 'generate_system_diagram_cy.html');
    const figDir   = path.resolve(__dirname, '..', 'docs', 'figures');

    if (!fs.existsSync(htmlPath)) {
        console.error('ERROR: HTML template not found at', htmlPath);
        process.exit(1);
    }
    if (!fs.existsSync(figDir)) {
        fs.mkdirSync(figDir, { recursive: true });
    }

    // Import puppeteer
    let puppeteer;
    try {
        puppeteer = require('puppeteer');
    } catch (e) {
        console.error('Puppeteer not found. Install with: npm install puppeteer');
        console.error('Or run via: npx -y puppeteer node scripts/export_system_diagram.js');
        process.exit(1);
    }

    const browser = await puppeteer.launch({
        headless: true,
        args: ['--no-sandbox', '--disable-setuid-sandbox'],
    });

    const page = await browser.newPage();
    await page.setViewport({ width: 1600, height: 1000 });

    // Load the HTML file
    const fileUrl = 'file:///' + htmlPath.replace(/\\/g, '/');
    await page.goto(fileUrl, { waitUntil: 'networkidle0', timeout: 30000 });

    // Wait for Cytoscape to be ready
    await page.waitForFunction('typeof window.getCyPNG === "function"', { timeout: 10000 });
    // Give Cytoscape a moment to settle layout
    await sleep(500);

    async function exportVariant(showCanaAox, filename) {
        // Set the edge visibility
        await page.evaluate((show) => {
            window.setShowCanaAox(show);
        }, showCanaAox);

        // Wait for render to settle
        await sleep(300);

        // Get PNG data URL from Cytoscape
        const pngDataUrl = await page.evaluate(() => {
            return window.getCyPNG(5);
        });

        // Convert data URL to buffer
        const base64 = pngDataUrl.replace(/^data:image\/png;base64,/, '');
        const buffer = Buffer.from(base64, 'base64');

        const outPath = path.join(figDir, filename);
        fs.writeFileSync(outPath, buffer);
        console.log('  ' + filename + ' (' + Math.round(buffer.length / 1024) + ' KB)');
    }

    console.log('Exporting system diagrams...');

    if (doPaper) {
        await exportVariant(true, 'system_diagram_paper1.png');
    }
    if (doDeck) {
        await exportVariant(true, 'system_diagram_deck.png');
    }

    await browser.close();
    console.log('Done.');
}

main().catch(err => {
    console.error('ERROR:', err.message);
    process.exit(1);
});
