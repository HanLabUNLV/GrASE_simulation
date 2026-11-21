#!/usr/bin/env python3
"""One SVG per manuscript figure, in figs/manuscript_svg/.

Script-drawn figures come from the figs/svg/ exports (87_figs_export.sh, run
it first). Figures 6 and 7 are cropped to what the manuscript shows. Figures assembled by
hand in the manuscript (graphical abstract, 1, 2 and 5) have no exported master, so their SVG is
taken from the manuscript docx itself; Word stores those with embedded raster
images, so they are not fully vector.

Usage: python3 -I scripts/88_manuscript_svg.py [manuscript.docx]
"""
import os, re, shutil, sys, zipfile
os.chdir("/mnt/data1/home/mirahan/GrASE_simulation")
DOCX = sys.argv[1] if len(sys.argv) > 1 else "/mnt/data1/home/mirahan/GrASE/Bioinformatics_manuscript.updated.docx"
S, OUT = "figs/svg", "figs/manuscript_svg"
os.makedirs(OUT, exist_ok=True)

def nested(f, tag, x, y):
    s = open(f).read()
    s = re.sub(r"<\?xml[^>]*\?>\s*", "", s)
    s = re.sub(r'id="([^"]+)"', lambda m: 'id="%s_%s"' % (tag, m.group(1)), s)
    s = re.sub(r'url\(#([^)]+)\)', lambda m: 'url(#%s_%s)' % (tag, m.group(1)), s)
    s = re.sub(r'href="#([^"]+)"', lambda m: 'href="#%s_%s"' % (tag, m.group(1)), s)
    m = re.search(r'<svg\b[^>]*\bwidth="([\d.]+)(?:pt)?"[^>]*\bheight="([\d.]+)(?:pt)?"', s)
    s = s.replace("<svg", '<svg x="%g" y="%g"' % (x, y), 1)
    return s, float(m.group(1)), float(m.group(2))

def crop(src, box, out, labels=(), fs=16):
    """Show only box = (x0, y0, x1, y1) of an export; labels are (text, y) in
    the export's own coordinates, drawn in a margin added at the left."""
    x0, y0, x1, y1 = box
    m = 24 if labels else 0
    inner, _, _ = nested(os.path.join(S, src), "c", m, 0)
    inner = re.sub(r'(<svg\b[^>]*?)\swidth="[^"]*"\s+height="[^"]*"\s+viewBox="[^"]*"',
                   r'\1 width="%g" height="%g" viewBox="%g %g %g %g"' % (x1 - x0, y1 - y0, x0, y0, x1 - x0, y1 - y0),
                   inner, count=1)
    W, H = x1 - x0 + m, y1 - y0
    t = "".join('<text x="2" y="%g" font-family="Arial, sans-serif" font-weight="bold" font-size="%d">%s</text>\n'
                % (y - y0 + fs, fs, lab) for lab, y in labels)
    open(os.path.join(OUT, out), "w").write(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        'width="%gpt" height="%gpt" viewBox="0 0 %g %g">\n<rect width="100%%" height="100%%" fill="white"/>\n'
        '%s%s\n</svg>\n' % (W, H, W, H, t, inner))

def copy(src, out):
    shutil.copy(os.path.join(S, src), os.path.join(OUT, out))

copy("AIF1_ev11.svg", "Figure9_AIF1.svg")
copy("ARHGAP9_ev73.svg", "Figure10_ARHGAP9.svg")
## Figures 6 and 7 are crops of an export, as placed in the manuscript:
## Figure 6 drops the Unit column, the title and the footnote and adds row
## labels; Figure 7 drops the title and the footnote. Boxes are in pt.
crop("pr_three_levels.GT_rule.svg", (258, 20, 1051, 527), "Figure6_pr_three_levels.svg",
     labels=[("a.", 22), ("b.", 276)])
crop("roc_partial.GT_rule.tolerant.svg", (0, 36, 784, 372), "Figure7_roc_partial.svg")
## hand-assembled figures: the docx's own SVG for each picture, in document order
z = zipfile.ZipFile(DOCX)
doc = z.read("word/document.xml").decode()
rels = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="(media/[^"]+)"', z.read("word/_rels/document.xml.rels").decode()))
svgs = [rels[r] for r in re.findall(r'svgBlip[^>]*r:embed="(rId\d+)"', doc)]
print("docx SVG pictures in order:", svgs)
## order: graphical abstract, Figure 1, Figure 2 (algorithm), Figure 5
for media, out in zip(svgs[0:4], ["graphical_abstract.docx.svg", "Figure1_lower_set.docx.svg", "Figure2_algorithm.docx.svg",
                                  "Figure5_lfc_filter.docx.svg"]):
    open(os.path.join(OUT, out), "wb").write(z.read("word/" + media))
for f in sorted(os.listdir(OUT)):
    print("  %-36s %8.0f KB" % (f, os.path.getsize(os.path.join(OUT, f)) / 1024))
