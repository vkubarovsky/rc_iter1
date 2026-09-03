"""Read an SVG produced by pdftocairo, carrying transforms down the group tree."""
import re
import numpy as np

TOK = re.compile(r'<g([^>]*)>|</g>|<path([^>]*)/>')
MAT = re.compile(r'matrix\(\s*([-\d.eE]+)[, ]+([-\d.eE]+)[, ]+([-\d.eE]+)[, ]+'
                 r'([-\d.eE]+)[, ]+([-\d.eE]+)[, ]+([-\d.eE]+)\s*\)')

def mul(A, B):
    a,b,c,d,e,f = A; p,q,r,s,t,u = B
    return (a*p+c*q, b*p+d*q, a*r+c*s, b*r+d*s, a*t+c*u+e, b*t+d*u+f)

def apply(M, x, y):
    a,b,c,d,e,f = M
    return a*x+c*y+e, b*x+d*y+f

def read(path):
    s = open(path).read()
    I = (1,0,0,1,0,0)
    stack = [(I, None, None)]          # transform, fill, stroke
    out = []
    for m in TOK.finditer(s):
        if m.group(0).startswith('<g'):
            a = m.group(1)
            M = stack[-1][0]
            t = MAT.search(a)
            if t: M = mul(M, tuple(float(z) for z in t.groups()))
            fl = re.search(r'fill="([^"]*)"', a); st = re.search(r'stroke="([^"]*)"', a)
            stack.append((M, fl.group(1) if fl else stack[-1][1],
                             st.group(1) if st else stack[-1][2]))
        elif m.group(0) == '</g>':
            if len(stack) > 1: stack.pop()
        else:
            a = m.group(2)
            d = re.search(r'd="([^"]+)"', a)
            if not d: continue
            dd = d.group(1)
            pts = [(float(x), float(y)) for _, x, y in
                   re.findall(r'([MLC])\s*([-\d.eE]+)[ ,]+([-\d.eE]+)', dd)]
            if not pts: continue
            M = stack[-1][0]
            t = MAT.search(a)
            if t: M = mul(M, tuple(float(z) for z in t.groups()))
            P = [apply(M, x, y) for x, y in pts]
            xs = [p[0] for p in P]; ys = [p[1] for p in P]
            fl = re.search(r'fill="([^"]*)"', a); st = re.search(r'stroke="([^"]*)"', a)
            out.append(dict(x0=min(xs), x1=max(xs), y0=min(ys), y1=max(ys), n=len(P),
                            curved='C' in dd,
                            fill=fl.group(1) if fl else stack[-1][1],
                            stroke=st.group(1) if st else stack[-1][2]))
    return out
