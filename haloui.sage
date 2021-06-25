"""
This file compares the results of the WeilPolynomials iterator with the main results of [Hal10] and [HS12].
These two papers give descriptions of Weil polynomials in dimension 3 and 4 in terms of explicit inequalities.

* [Hal10] - Safia Haloui.  The characteristic polynomials of abelian varieties of dimensions 3 over finite fields, J. of Number Theory 130 (12), 2745--2752.
* [Hal12] - Safia Haloui and Vijaykumar Singh.  The characteristic polynomials of abelian varieties of dimension 4 over finite fields, in Arithmetic, Geometry, Cryptography and Coding Theory, AMS, Providence, 2012.
"""

from sage.rings.polynomial.weil.all import WeilPolynomials
from collections import defaultdict

def wp(f):
    # polynomial -> reversed list
    return tuple(reversed(list(f)))

def check_weil3(q):
    #dbpolys = set(tuple(poly) for poly in db.av_fq_isog.search({"g":3, "q":q}, "poly"))
    dbpolys = set(wp(f) for f in WeilPolynomials(6, q))
    comppolys = set()
    R.<t> = ZZ[]
    upper0 = (2*RR(q).sqrt()).floor()
    if q.is_square():
        upper0 -= 1
    lower0 = -upper0
    for beta in [lower0..upper0]:
        comppolys.add(wp((t^2-q)^2*(t^2+beta*t+q)))
    upper1 = (6*RR(q).sqrt()).floor()
    if q.is_square():
        upper1 -= 1
    lower1 = -upper1
    for a1 in [lower1..upper1]:
        lower2 = (4*RR(q).sqrt()*a1.abs()).ceil() - 9*q
        #if q.is_square():
        #    lower2 += 1
        upper2 = (a1^2/3 + 3*q).floor()
        for a2 in [lower2..upper2]:
            lower3a = (-2*a1^3/27 + a1*a2/3 + q*a1 - 2/27*RR(a1^2-3*a2+9*q)^(3/2)).ceil()
            upper3a = (-2*a1^3/27 + a1*a2/3 + q*a1 + 2/27*RR(a1^2-3*a2+9*q)^(3/2)).floor()
            lower3b = (-2*q*a1 - 2*RR(q).sqrt()*a2 - 2*q*RR(q).sqrt()).ceil()
            upper3b = (-2*q*a1 + 2*RR(q).sqrt()*a2 + 2*q*RR(q).sqrt()).floor()
            #if q.is_square() or a2 == -q:
            #    lower3b += 1
            #    upper3b -= 1
            for a3 in [max(lower3a, lower3b)..min(upper3a, upper3b)]:
                comppolys.add((1, a1, a2, a3, q*a2, q^2*a1, q^3))
    if comppolys != dbpolys:
        print("Mismatch q=%s!"%q, len(comppolys), len(dbpolys))
        for f in dbpolys:
            if f not in comppolys:
                F = R(list(reversed(f)))
                if not F.is_weil_polynomial():
                    print("ERROR!")
                print("Extra", F.factor())
        for f in comppolys:
            if f not in dbpolys:
                print("Missing", f)
#sage: check_weil3(2)
#sage: check_weil3(3)
#sage: check_weil3(4)
#Mismatch q=4! 1637 1641
#Extra (t - 2)^2 * (t^2 + 3*t + 4)^2
#Extra (t + 2)^2 * (t^2 - 3*t + 4)^2
#Extra (t - 2)^4 * (t^2 + 3*t + 4)
#Extra (t + 2)^4 * (t^2 - 3*t + 4)
#sage: check_weil3(5)
#Mismatch q=5! 2949 2953
#Extra (t^2 + 4*t + 5) * (t^2 - 3*t + 5)^2
#Extra (t^2 - 4*t + 5) * (t^2 + 3*t + 5)^2
#Extra (t^2 + 3*t + 5) * (t^2 - 4*t + 5)^2
#Extra (t^2 - 3*t + 5) * (t^2 + 4*t + 5)^2
#sage: check_weil3(7)
#Mismatch q=7! 7973 7979
#Extra (t^2 - 4*t + 7) * (t^2 + 4*t + 7)^2
#Extra (t^2 + 4*t + 7) * (t^2 - 3*t + 7)^2
#Extra (t^2 - 4*t + 7) * (t^2 + 3*t + 7)^2
#Extra (t^2 - t + 7) * (t^2 + 3*t + 7)^2
#Extra (t^2 + t + 7) * (t^2 - 3*t + 7)^2
#Extra (t^2 + 4*t + 7) * (t^2 - 4*t + 7)^2
#sage: check_weil3(8)
#Mismatch q=8! 11821 11823
#Extra (t^2 - 4*t + 8) * (t^2 + 3*t + 8)^2
#Extra (t^2 + 4*t + 8) * (t^2 - 3*t + 8)^2
#sage: check_weil3(9)
#Mismatch q=9! 17115 17121
#Extra (t^2 - 4*t + 9) * (t^2 + 3*t + 9)^2
#Extra (t^2 + 4*t + 9) * (t^2 - 3*t + 9)^2
#Extra (t - 3)^2 * (t^2 + 5*t + 9)^2
#Extra (t - 3)^4 * (t^2 + 5*t + 9)
#Extra (t + 3)^4 * (t^2 - 5*t + 9)
#Extra (t + 3)^2 * (t^2 - 5*t + 9)^2

def check_weil4(q):
    dbpolys = set(wp(f) for f in WeilPolynomials(8, q))
    a4s = defaultdict(list)
    for f in dbpolys:
        a4s[f[1:4]].append(f[4])
    comppolys = set()
    R.<t> = ZZ[]
    i = CC.gen()
    j = CC.zeta(3)
    if q.is_square():
        for h in WeilPolynomials(4, q, polring=R):
            comppolys.add((t^2-q.sqrt())^2*h)
            comppolys.add((t^2+q.sqrt())^2*h)
    upper1 = (8*RR(q).sqrt()).floor()
    #if q.is_square():
    #    upper1 -= 1
    lower1 = -upper1
    for a1 in [lower1..upper1]:
        lower2 = (6*RR(q).sqrt()*a1.abs()).ceil() - 20*q
        #if q.is_square():
        #    lower2 += 1
        upper2 = (3*a1^2/8 + 4*q).floor()
        for a2 in [lower2..upper2]:
            lower3a = (-9*q*a1 - 4*RR(q).sqrt()*a2 - 16*q*RR(q).sqrt()).ceil()
            upper3a = (-9*q*a1 + 4*RR(q).sqrt()*a2 + 16*q*RR(q).sqrt()).floor()
            #if q.is_square():
            #    lower3a += 1
            #    upper3a -= 1
            lower3b = (-a1^3/8 + a1*a2/2 + q*a1 - RR(2/3 * (3*a1^2/8 - a2 + 4*q))^(3/2)).ceil()
            upper3b = (-a1^3/8 + a1*a2/2 + q*a1 + RR(2/3 * (3*a1^2/8 - a2 + 4*q))^(3/2)).floor()
            for a3 in [max(lower3a, lower3b)..min(upper3a, upper3b)]:
                lower4a = (2*RR(q).sqrt()*(q*a1+a3).abs() - 2*q*a2 - 2*q^2).ceil()
                #if q.is_square():
                #    lower4a += 1
                y = -3*a1^2/8 + a2 - 4*q
                z = a1^3/8 - q*a1 - a1*a2/2 + a3
                if y == z == 0:
                    w = CC(0)
                else:
                    if z^2 == -8/27*y^3:
                        w3 = CC(8*y^6 + 540*y^3*z^2 - 729*z^4)
                    else:
                        w3 = 8*y^6 + 540*y^3*z^2 - 729*z^4 + 9*i*z.abs()*CC(-z^2 - 8/27*y^3)^(3/2)
                    w = 1/24*w3.abs()^(1/3)*CC(e)^(i*w3.argument()/3)
                lower4b = (9*a1^4/256 - 3*a1^2*a2/16 + a1*a3/4 + a2^2/6 + 2*q*a2/3 + 2*q^2/3 + 2*w.real()).ceil()
                assert (j*w + j^2*w.conjugate()).imag() < 10^-8
                upper4b = (9*a1^4/256 - 3*a1^2*a2/16 + a1*a3/4 + a2^2/6 + 2*q*a2/3 + 2*q^2/3 + (j*w + j^2*w.conjugate()).real()).floor()
                print(max(lower4a, lower4b),"-", upper4b, a4s[a1,a2,a3])
                for a4 in [max(lower4a, lower4b)..upper4b]:
                    print(a1,a2,a3,a4)
                    comppolys.add((1, a1, a2, a3, a4, q*a3, q^2*a2, q^3*a1, q^4))
    if comppolys != dbpolys:
        print("Mismatch q=%s!"%q, len(comppolys), len(dbpolys))
        for f in dbpolys:
            if f not in comppolys:
                F = R(list(reversed(f)))
                if not F.is_weil_polynomial():
                    print("ERROR!")
                print("Extra", F.factor())
        for f in comppolys:
            if f not in dbpolys:
                print("Missing", f)
#sage: check_weil4(2)
#169 - 167 []
#129 - 128 []
#136 - 136 [136]
#-8 32 -80 136
#...
#97 - 95 [96]
#105 - 103 [104]
#...
#Mismatch q=2! 5 1645
#Extra (t^2 + t + 2) * (t^6 + 3*t^5 + 2*t^4 - t^3 + 4*t^2 + 12*t + 8)
#Extra (t^2 - 2*t + 2) * (t^6 - t^5 - t^4 + 5*t^3 - 2*t^2 - 4*t + 8)
#Extra (t^4 - t^3 + t^2 - 2*t + 4)^2
#...
