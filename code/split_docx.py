#!/usr/bin/env python3
"""
Split the Quarto-rendered Word manuscript into the two files Genome Medicine
wants as separate uploads:

    main_manuscript.docx    title page .. Additional files
    additional_file_1.docx  Supplementary Materials onward

Approach: a *surgical* split of the valid Quarto-generated index.docx at the
XML level. We trim the <w:body> children at the "Supplementary Materials"
heading and keep every other part of the package (styles, numbering, media,
relationships) exactly as Quarto wrote them. This avoids the lossy
docx -> markdown -> docx pandoc round-trip, which re-creates the figure tables
and produced files Word flagged as "unreadable content".

While trimming we also:
  * balance bookmarks (drop any bookmarkStart/End whose partner fell in the
    other half) so Word does not see dangling bookmarks, and
  * renumber supplementary floats to "S" form (Fig. S1, Table S1, ...) to
    match the PDF. Each figure/table label is a single text run
    ("Figure 2:", "Fig. 1") so this is a per-run regex. A float is
    "supplementary" iff its number exceeds the count of same-type floats whose
    captions appear in the main section.

Usage:  python3 code/split_docx.py <index.docx> <main_manuscript.docx> <additional_file_1.docx>
Requires: lxml (preserves the OOXML namespace prefixes on write).
"""
import os, re, sys, zipfile
from lxml import etree

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
def w(tag): return '{%s}%s' % (W, tag)

HEADING = 'Supplementary Materials'
SP = r'[\s ]+'          # ASCII or non-breaking space (Quarto uses NBSP)
FIG_RE = re.compile(r'(\b)(Figure|Fig\.)(' + SP + r')(\d+)')
TBL_RE = re.compile(r'(\b)(Table)(' + SP + r')(\d+)')
CAP_FIG = re.compile(r'Figure' + SP + r'(\d+):')
CAP_TBL = re.compile(r'Table' + SP + r'(\d+):')

SRC, MAIN_OUT, SUPP_OUT = sys.argv[1], sys.argv[2], sys.argv[3]
DOC_BYTES = zipfile.ZipFile(SRC).read('word/document.xml')


def ptext(p):
    return ''.join(t.text or '' for t in p.iter(w('t')))


def split_index(body):
    for i, c in enumerate(body):
        if c.tag == w('p') and ptext(c).strip() == HEADING:
            return i
    raise SystemExit("ERROR: '%s' heading not found in document.xml" % HEADING)


def build(keep_main):
    """Return a document root containing only the main- or supplementary-half body."""
    root = etree.fromstring(DOC_BYTES)
    body = root.find(w('body'))
    kids = list(body)
    idx = split_index(body)
    sectPr = kids[-1] if kids[-1].tag == w('sectPr') else None
    keep = (kids[:idx] + ([sectPr] if sectPr is not None else [])) if keep_main else kids[idx:]
    for c in kids:
        body.remove(c)
    for c in keep:
        body.append(c)
    balance_bookmarks(root)
    return root


def balance_bookmarks(root):
    starts = {b.get(w('id')) for b in root.iter(w('bookmarkStart'))}
    ends = {b.get(w('id')) for b in root.iter(w('bookmarkEnd'))}
    for b in list(root.iter(w('bookmarkStart'))):
        if b.get(w('id')) not in ends:
            b.getparent().remove(b)
    for b in list(root.iter(w('bookmarkEnd'))):
        if b.get(w('id')) not in starts:
            b.getparent().remove(b)


def caption_max(root, pat):
    mx = 0
    for t in root.iter(w('t')):
        if t.text:
            for m in pat.finditer(t.text):
                mx = max(mx, int(m.group(1)))
    return mx


def renumber(root, mf, mt):
    def fig(m):
        n = int(m.group(4))
        return m.group(1) + m.group(2) + m.group(3) + (('S%d' % (n - mf)) if n > mf else str(n))
    def tbl(m):
        n = int(m.group(4))
        return m.group(1) + m.group(2) + m.group(3) + (('S%d' % (n - mt)) if n > mt else str(n))
    for t in root.iter(w('t')):
        if not t.text:
            continue
        s = TBL_RE.sub(tbl, FIG_RE.sub(fig, t.text))
        if s != t.text:
            t.text = s


def write(out, root):
    new = etree.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)
    tmp = out + '.tmp'
    with zipfile.ZipFile(SRC) as zi, zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zo:
        for item in zi.infolist():
            data = new if item.filename == 'word/document.xml' else zi.read(item.filename)
            zo.writestr(item, data)
    os.replace(tmp, out)


main_root = build(keep_main=True)
mf = caption_max(main_root, CAP_FIG)
mt = caption_max(main_root, CAP_TBL)
sys.stderr.write("Main section defines %d figures and %d tables; "
                 "numbering the rest as S1, S2, ...\n" % (mf, mt))
renumber(main_root, mf, mt)
write(MAIN_OUT, main_root)

supp_root = build(keep_main=False)
renumber(supp_root, mf, mt)
write(SUPP_OUT, supp_root)

print("Wrote", MAIN_OUT)
print("Wrote", SUPP_OUT)
