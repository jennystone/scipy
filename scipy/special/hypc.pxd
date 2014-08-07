cimport sf_error

cdef extern from "c_misc/misc.h":
    double poch(double x, double m) nogil

cdef extern int sgngam
cdef extern double MACHEP

cdef extern from "cephes.h":
    double Gamma ( double x ) nogil
    double lgam ( double x ) nogil
    double psi ( double x ) nogil
    
from _complexstuff cimport *
from libc.math cimport cos, sqrt, exp, fabs, round, M_PI as pi, INFINITY, pow, log, isnan
from libc.stdlib cimport abs

cdef extern from "complex.h":
     double cabs (double complex z) nogil
     double creal (double complex z) nogil
     double cimag (double complex z) nogil
     double complex cpow(double complex x, double complex y) nogil
     double complex ctan (double complex z) nogil
     double complex clog (double complex z) nogil
 

cdef inline double complex chyp2f1(double a, double b, double c, double complex x) nogil:

    cdef double d, e, f, g, h, m
    cdef double p, q, r, ax, k, umax
    cdef double complex u, s, sum, y, d1, d2
    cdef double  err, EL = 0.5772156649015329
    cdef double t1
    cdef double ETHRESH = 1.0e-12, EPS = 1.0e-13
    cdef int MAX_ITER = 10000
    cdef double complex sum1, sum2, sum3, sum4, sum5,sum6
    cdef double complex term1, term2, term3, term4, term5, term6, term7 
    cdef int i, aid, max, ia, ib, ic, id, j, mab, nca,ncb
    cdef int neg_int_a = 0, neg_int_b = 0, neg_int_c = 0
    cdef int neg_int_ca = 0, neg_int_cb = 0

    err = 0.0
    ax = cabs(x)
    s = 1.0 - x
    ia = <int>round(a)
    ib = <int>round(b)
    ic = <int>round(c)

    d = c - a - b
    id = <int>round(d)

    if x == 0.0:
        return 1.0

    if ((a == 0 or b == 0) and c != 0):
        return 1.0
    if (a <= 0 and fabs(a - ia) < EPS):  #a is a negative integer
        neg_int_a = 1

    if (b <= 0 and fabs(b - ib) < EPS):  #b is a negative integer
        neg_int_b = 1

    if (c <= 0 and fabs(c - ic) < EPS):
        neg_int_c = 1

    if neg_int_c and not (neg_int_a or neg_int_b):
        return INFINITY

    p = c - a
    ia = <int>round(p)

    r = c - b
    ib = <int>round(r)
   
    if fabs(1 - creal(x)) < EPS and cimag(x) == 0.0:
        if d < 0:
            return INFINITY
        if d > 0:
            return poch(p, a)*poch(r, -a)
 	
    if ((ia <= 0.0) and (fabs(p - ia) < EPS)):
        neg_int_ca = 1

    if ((ib <= 0.0) and (fabs(r - ib) < EPS)):
        neg_int_cb = 1

    if neg_int_a or neg_int_b:
        if neg_int_a:
            max = <int>fabs(a)
        else:
            max = <int>fabs(b)
        if neg_int_c:
            if fabs(c) < max:
                return INFINITY
        u = 1.0
        sum = 1.0
        for i in range(1, max+1):
            u *= (a+i-1.0)*(b+i-1.0)*x/(i*(c+i-1.0))
            sum += u
        return sum 

    if neg_int_ca or neg_int_cb:
        if neg_int_ca:
            max = <int>fabs(c - a)
        else:
            max = <int>fabs(c - b)
        if neg_int_c:
            if fabs(c) < max:
                return INFINITY
        u = 1.0
        sum = 1.0
        for i in range(1, max+1):
            u *= (c-a+i-1.0)*(c-b+i-1.0)*x/(i*(c+i-1.0))
            sum += u
        sum *= s**d
        return sum

    if fabs(b - c) < 10**-13:
        return cpow(s, -a)

    if fabs(a - c) < 10**-13:
        return cpow(s, -b)

    if ax > 1.0:
        EPS = 10**-8
        if a != b:
            mab = <int> (a - b + EPS*((a - b)/fabs(a - b)))
        else:
            mab = <int> (a - b)
        if (fabs(a - b - mab) < EPS) and (ax < 1.1):
            b += EPS
        if fabs(a - b - mab) > EPS:
            term1 = Gamma(c)*Gamma(b - a)/(Gamma(p)*Gamma(b)*cpow(-x, a))
            term2 = Gamma(c)*Gamma(a - b)/(Gamma(r)*Gamma(a)*cpow(-x, b))
            sum = term1 + term2
            m = 1
            u = sum
            while (cabs(u/sum) > EPS):
                term1*=(a + m - 1.0)*(-p + 1.0)/((a - b + m)*m*x)
                term2*=(b + m - 1.0)*(-r + 1.0)/((b - a + m)*m*x)
                sum += term1 + term2
                u = term1 + term2
                m += 1
                if  m > MAX_ITER:
                    return sum
            return sum
        else:
            if a-b < 0:
                a, b = b, a
            nca = <int> (p + EPS*p/fabs(p))
            ncb = <int> (r + EPS*r/fabs(r))
            #if fabs(p - nca) < EPS or fabs(r - ncb) < EPS:
            c += EPS

            mab = <int>(a - b + EPS)
            term1 = Gamma(c)/(Gamma(a)*cpow(-x, b))
            sum = Gamma(a - b)/Gamma(r)*term1
            for i in range(1, mab):
               term1*=(b + i -1.0)/(i*x)
               sum += term1*Gamma(a - b - i)/Gamma(r - i)

            if mab == 0:
               sum = 0
            term2 = Gamma(c)/(Gamma(a)*Gamma(r)*(-x)**a)
            term3 = -2.0*EL-psi(a)-psi(p)
            for i in range(1,mab+1):
                term3 += 1.0/i
            term4 = 1.0
            for i in range(1,mab+1):
                term4 *= (b + i - 1.0)*(b - c + i)/i
            sum1 = term4*(term3 + clog(-x))*term2
            term6 = 1
            sum3 = 0.0
            sum4 = 0.0
            for i in range(1, MAX_ITER):
                term2 = term2/x
                term6 *= (b + i - 1.0)*(b - c+ i)/(i*i)
                term7 = term6
                for j in range(i + 1, i+mab+1):
                    term7 *= (b + j -1.0)*( b - c + j)/j
                sum3 += (a-1.0)/(i*(a+i-1.0))+(-p - 1.0)/(i*(-p+i-1.0))
                sum5 = sum3
                for j in range(i + 1, i+mab+1):
                    sum5 += 1.0/j

                sum6 = -2.0*EL-psi(a)-psi(-p)+sum5-1.0/(i-p)-pi/ctan(pi*(i-p))+clog(-x)
                sum1+=sum6*term2*term7
            return sum1 + sum
                
    if d < 0.0:
        y = hys2f1(a, b, c, x, &err)
        if err < ETHRESH:
            return y
        err = 0.0
        aid = 2 - id
        e = c + aid
        d2 = chyp2f1(a, b, e, x)
        d1 = chyp2f1(a, b, e + 1.0, x)
        q = a + b + 1.0
        for i in range(1, aid):
            r = e - 1.0
            y = (e * (r - (2.0 * e - q) * x) * d2 +
                 (e - a) * (e - b) * x * d1) / (e * r * s)
            e = r
            d1 = d2
            d2 = y
        if err > ETHRESH:
            sf_error.error("chyp2f1", sf_error.LOSS, "Loss of precision")
        return y


    if neg_int_ca or neg_int_cb:
        y = cpow(s, d) * hys2f1(c - a, c - b, c, x, &err)
        if err > ETHRESH:
            sf_error.error("chyp2f1", sf_error.LOSS, "Loss of precision")
        return y


  
cdef inline double complex hys2f1(double a, double b, double c, double complex x, double *loss) nogil:
    cdef double ETHRESH = 1.0e-12, EPS = 1.0e-13
    cdef int MAX_ITER = 10000
    cdef double f, g, h, k, m, umax, prec
    cdef double complex u, s
    cdef int i
    cdef int ib, ic, intflag = 0

    if fabs(b) > fabs(a):
        a, b = b, a

    ib = <int>round(b);
    ic = <int>round(c);

    if fabs(b - ib) < EPS and ib <= 0 and fabs(b) < fabs(a):
        a, b = b, a
        intflag = 1

    if ((fabs(a) > fabs(c) + 1 or intflag) and fabs(c - a) > 2
	and fabs(a) > 2) :
        return hyp2f1ra(a, b, c, x, loss)

    i = 0
    umax = 0.0
    f = a
    g = b
    h = c
    s = 1.0
    u = 1.0
    k = 0.0
    prec = MACHEP

    if c < 0 and not(a < 0 or b < 0) and fabs( c - ic) > EPS: 
        prec = pow(10,-700)

    while (s == 0 or cabs(u/s) > prec): 
        if fabs(h) < EPS:
            loss[0] = 1.0 
            return INFINITY
        m = k + 1.0
        u = u*((f + k)*(g + k)*x/((h + k)*m))
        s += u
        k = cabs(u)
        if (k > umax):
            umax = k
        k = m
        i = i + 1
        if  i > MAX_ITER:
            loss[0] = 1.0
            return s
    loss[0] = (MACHEP * umax)/cabs(s) + (MACHEP*i)
    return (s)


cdef inline double complex hyp2f1ra(double a, double b, double c, double complex x, double *loss) nogil:
    cdef double complex f2, f1, f0
    cdef double ETHRESH = 1.0e-12, EPS = 1.0e-13
    cdef int MAX_ITER = 10000
    cdef int n, da
    cdef double t, err

    if (c < 0 and a <= c) or (c >= 0 and a >= c):
        da = <int>round(a - c)

    else: 
        da = <int>round(a)

    t = a - da
    loss[0] = 0
    if (fabs(da) > MAX_ITER):
        loss[0] = 1.0
        return nan

    if da < 0:
        f2 = 0
        f1 = hys2f1(t, b, c, x, &err)
        loss[0] += err
        f0 = hys2f1(t - 1, b, c, x, &err)
        loss[0] += err
        t -= 1
        for n in range(1, -da): 
            f2 = f1
            f1 = f0
            f0 = -(2*t - c - t*x + b*x)/(c - t)*f1 - t*(x - 1)/(c - t)*f2
            t -= 1
	
    else:
        f2 = 0
        f1 = hys2f1(t, b, c, x, &err)
        loss[0] += err
        f0 = hys2f1(t + 1, b, c, x, &err)
        loss[0] += err
        t += 1
        for n in range(1, da):
            f2 = f1
            f1 = f0
            f0 = -((2*t - c - t*x + b*x)*f1 + (c - t)*f2)/(t*(x - 1))
            t += 1
    return f0

