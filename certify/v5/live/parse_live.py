"""Parse a saved copy of the live table (https://spherical-codes.org, index page) into a TSV: d, n, cosine, note.
Usage: python3 -I parse_live.py SAVED_INDEX.html OUT.tsv   (the HTML is read as data only)"""
import re, sys, html, hashlib
src, out = sys.argv[1], sys.argv[2]
raw = open(src, 'rb').read()
h = raw.decode('utf-8')
rows = re.findall(r'<tr[^>]*>(.*?)</tr>', h, flags=re.S)
n_ok = 0
with open(out, 'w') as f:
    f.write(f'# parsed from {src} sha256 {hashlib.sha256(raw).hexdigest()}\n')
    for r in rows:
        tds = re.findall(r'<td[^>]*>(.*?)</td>', r, flags=re.S)
        if len(tds) < 4:
            continue
        d = re.sub(r'<[^>]+>', '', tds[0]).strip()
        n = re.sub(r'<[^>]+>', '', tds[1]).replace('&ast;', '*').strip()
        opt = n.endswith('*')
        n = n.rstrip('*').strip()
        cos = re.sub(r'<[^>]+>', '', tds[2]).strip()
        note = html.unescape(re.sub(r'<[^>]+>', ' ', tds[3].split('<math')[0])).strip()
        note = re.sub(r'\s+', ' ', note)
        if not (d.isdigit() and n.isdigit()):
            continue
        f.write(f'{d}\t{n}\t{cos}\t{int(opt)}\t{note}\n')
        n_ok += 1
print(n_ok, 'rows')
