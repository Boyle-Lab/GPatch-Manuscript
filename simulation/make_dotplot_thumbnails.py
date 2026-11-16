#!/usr/bin/env python3
"""Tile dot-plot PDFs as vector thumbnails onto US Letter pages.

Each set (label + directory) starts on a new page with a header line; thumbnails are
placed in a rows x cols grid inside the margins, in chromosome order, each captioned
with its chromosome. Requires pypdf.

Usage (from manuscript/simulation):
  ./make_dotplot_thumbnails.py -o dotplot_thumbnails.pdf \\
      "HG002 - GPatch=T2T-HG002/GPatch" "HG002 - RagTag scaffold=T2T-HG002/RagTag/SCAFFOLD" ...
"""

__version__ = "0.1.0"

import argparse
import glob
import os
import re

from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import (DecodedStreamObject, DictionaryObject, NameObject)

PT = 72.0


def chrom_key(path):
    c = re.search(r"\.(chr[^.]+)\.pdf$", path).group(1)[3:]
    return (0, int(c), "") if c.isdigit() else (1, 0, c)


def pdf_str(s):
    return s.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def add_text(writer, page, items):
    """Draw (x, y, size, text, bold, centered) items with standard Helvetica fonts."""
    fonts = DictionaryObject()
    for key, base in (("/FH", "/Helvetica"), ("/FB", "/Helvetica-Bold")):
        f = DictionaryObject({NameObject("/Type"): NameObject("/Font"),
                              NameObject("/Subtype"): NameObject("/Type1"),
                              NameObject("/BaseFont"): NameObject(base)})
        fonts[NameObject(key)] = writer._add_object(f)
    res = page.setdefault(NameObject("/Resources"), DictionaryObject())
    res[NameObject("/Font")] = fonts
    ops = []
    for x, y, size, text, bold, centered in items:
        if centered:  # approximate Helvetica width (0.5 em average) for centring
            x -= 0.5 * size * len(text) * (0.55 if bold else 0.5)
        ops.append("BT /%s %.1f Tf %.2f %.2f Td (%s) Tj ET"
                   % ("FB" if bold else "FH", size, x, y, pdf_str(text)))
    stream = DecodedStreamObject()
    stream.set_data("\n".join(ops).encode("latin-1"))
    page[NameObject("/Contents")] = writer._add_object(stream)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sets", nargs="+", help='"Label=directory" pairs, in page order.')
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("-g", "--glob", default="*.dotplot.*.chr*.pdf")
    ap.add_argument("--rows", type=int, default=4)
    ap.add_argument("--cols", type=int, default=4)
    ap.add_argument("--width", type=float, default=8.5, help="Page width (in).")
    ap.add_argument("--height", type=float, default=11.0, help="Page height (in).")
    ap.add_argument("--margin", type=float, default=0.5, help="Margin (in).")
    args = ap.parse_args()

    W, H, M = args.width * PT, args.height * PT, args.margin * PT
    header_h, caption_h, pad = 18.0, 10.0, 4.0
    cell_w = (W - 2 * M) / args.cols
    cell_h = (H - 2 * M - header_h) / args.rows
    size = min(cell_w - 2 * pad, cell_h - caption_h - 2 * pad)
    per_page = args.rows * args.cols

    writer = PdfWriter()
    for spec in args.sets:
        label, d = spec.split("=", 1)
        files = sorted(glob.glob(os.path.join(d, args.glob)), key=chrom_key)
        if not files:
            raise SystemExit("no files matching %s in %s" % (args.glob, d))
        npages = -(-len(files) // per_page)
        for p in range(npages):
            page = writer.add_blank_page(W, H)
            title = label if npages == 1 else "%s (%d/%d)" % (label, p + 1, npages)
            text = [(M, H - M - 12, 11, title, True, False)]
            for i, f in enumerate(files[p * per_page:(p + 1) * per_page]):
                r, c = divmod(i, args.cols)
                x0 = M + c * cell_w
                y_top = H - M - header_h - r * cell_h
                src = PdfReader(f).pages[0]
                bw, bh = float(src.mediabox.width), float(src.mediabox.height)
                s = size / max(bw, bh)
                tx = x0 + (cell_w - bw * s) / 2 - float(src.mediabox.left) * s
                ty = y_top - caption_h - pad - bh * s - float(src.mediabox.bottom) * s
                text.append((x0 + cell_w / 2, y_top - caption_h + 1, 8,
                             re.search(r"\.(chr[^.]+)\.pdf$", f).group(1), False, True))
                page.merge_transformed_page(src, Transformation().scale(s, s).translate(tx, ty))
            # Captions and header drawn on top in their own content stream.
            overlay = writer.add_blank_page(W, H)
            add_text(writer, overlay, text)
            page.merge_page(overlay)
            writer.remove_page(len(writer.pages) - 1)
    with open(args.out, "wb") as fh:
        writer.write(fh)


if __name__ == "__main__":
    main()
