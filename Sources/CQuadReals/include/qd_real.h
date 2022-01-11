/*
 * include/qd_real.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Quad-double precision (>= 212-bit significand) floating point arithmetic
 * package, written in ANSI C++, taking full advantage of operator overloading.
 * Uses similar techniques as that of David Bailey's double-double package 
 * and that of Jonathan Shewchuk's adaptive precision floating point 
 * arithmetic package.  See
 *
 *   http://www.nersc.gov/~dhbailey/mpdist/mpdist.html
 *   http://www.cs.cmu.edu/~quake/robust.html
 *
 * for more details.
 *
 * Yozo Hida
 */
#ifndef _QD_QD_REAL_H
#define _QD_QD_REAL_H

//#include <iostream>
//#include <string>
//#include <limits>
#include <math.h>
#include "qd_config.h"
#include "dd_real.h"

typedef struct QD_API qd_real {
  double x[4];    /* The Components. */
} qd_real;

  /* Eliminates any zeros in the middle component(s). */
//  void zero_elim();
//  void zero_elim(double &e);
//
//  void renorm();
//  void renorm(double &e);
//
//  void quick_accum(double d, double &e);
//  void quick_prod_accum(double a, double b, double &e);
//
//  qd_real(double x0, double x1, double x2, double x3);
//  explicit qd_real(const double *xx);
//
  static const qd_real qd_2pi;
  static const qd_real qd_pi;
  static const qd_real qd_3pi4;
  static const qd_real qd_pi2;
  static const qd_real qd_pi4;
  static const qd_real qd_e;
  static const qd_real qd_log2;
  static const qd_real qd_log10c;
  static const qd_real qd_nan;
  static const qd_real qd_inf;

  static const double qd_eps;
  static const double qd_min_normalized;
  static const qd_real qd_max;
  static const qd_real qd_safe_max;
  static const int qd_ndigits;

inline void qd_quick_renorm(double *c0, double *c1, double *c2, double *c3, double *c4) {
  double t0, t1, t2, t3;
  double s;
  s  = quick_two_sum(*c3, *c4, &t3);
  s  = quick_two_sum(*c2, s , &t2);
  s  = quick_two_sum(*c1, s , &t1);
  *c0 = quick_two_sum(*c0, s , &t0);

  s  = quick_two_sum(t2, t3, &t2);
  s  = quick_two_sum(t1, s , &t1);
  *c1 = quick_two_sum(t0, s , &t0);

  s  = quick_two_sum(t1, t2, &t1);
  *c2 = quick_two_sum(t0, s , &t0);
  
  *c3 = t0 + t1;
}

inline void qd_renorm(double *c0, double *c1, double *c2, double *c3) {
  double s0, s1, s2 = 0.0, s3 = 0.0;

//  if (QD_ISINF(c0)) return;

  s0 = quick_two_sum(*c2, *c3, c3);
  s0 = quick_two_sum(*c1, s0, c2);
  *c0 = quick_two_sum(*c0, s0, c1);

  s0 = *c0;
  s1 = *c1;
  if (s1 != 0.0) {
    s1 = quick_two_sum(s1, *c2, &s2);
    if (s2 != 0.0)
      s2 = quick_two_sum(s2, *c3, &s3);
    else
      s1 = quick_two_sum(s1, *c3, &s2);
  } else {
    s0 = quick_two_sum(s0, *c2, &s1);
    if (s1 != 0.0)
      s1 = quick_two_sum(s1, *c3, &s2);
    else
      s0 = quick_two_sum(s0, *c3, &s1);
  }

  *c0 = s0;
  *c1 = s1;
  *c2 = s2;
  *c3 = s3;
}

inline void qd_renorm2(double *c0, double *c1, double *c2, double *c3, double *c4) {
  double s0, s1, s2 = 0.0, s3 = 0.0;

//  if (QD_ISINF(c0)) return;

  s0 = quick_two_sum(*c3, *c4, c4);
  s0 = quick_two_sum(*c2, s0, c3);
  s0 = quick_two_sum(*c1, s0, c2);
  *c0 = quick_two_sum(*c0, s0, c1);

  s0 = *c0;
  s1 = *c1;

  if (s1 != 0.0) {
    s1 = quick_two_sum(s1, *c2, &s2);
    if (s2 != 0.0) {
      s2 = quick_two_sum(s2, *c3, &s3);
      if (s3 != 0.0)
        s3 += *c4;
      else
        s2 = quick_two_sum(s2, *c4, &s3);
    } else {
      s1 = quick_two_sum(s1, *c3, &s2);
      if (s2 != 0.0)
        s2 = quick_two_sum(s2, *c4, &s3);
      else
        s1 = quick_two_sum(s1, *c4, &s2);
    }
  } else {
    s0 = quick_two_sum(s0, *c2, &s1);
    if (s1 != 0.0) {
      s1 = quick_two_sum(s1, *c3, &s2);
      if (s2 != 0.0)
        s2 = quick_two_sum(s2, *c4, &s3);
      else
        s1 = quick_two_sum(s1, *c4, &s2);
    } else {
      s0 = quick_two_sum(s0, *c3, &s1);
      if (s1 != 0.0)
        s1 = quick_two_sum(s1, *c4, &s2);
      else
        s0 = quick_two_sum(s0, *c4, &s1);
    }
  }

  *c0 = s0;
  *c1 = s1;
  *c2 = s2;
  *c3 = s3;
}
//}

//inline void qd_renorm() {
//  renorm(x[0], x[1], x[2], x[3]);
//}
//
//inline void qd_renorm(double *e) {
//  renorm(x[0], x[1], x[2], x[3], e);
//}


/********** Additions ************/
//namespace qd {

inline void three_sum(double *a, double *b, double *c) {
    double t1, t2, t3;
    t1 = two_sum(*a, *b, &t2);
    *a  = two_sum(*c, t1, &t3);
    *b  = two_sum(t2, t3, c);
}

inline void three_sum2(double *a, double *b, double *c) {
    double t1, t2, t3;
    t1 = two_sum(*a, *b, &t2);
    *a  = two_sum(*c, t1, &t3);
    *b = t2 + t3;
}

//}

/* quad-double + double */
inline qd_real qd_addd (const qd_real *a, double b) {
    double c0, c1, c2, c3;
    double e;
    
    c0 = two_sum(a->x[0], b, &e);
    c1 = two_sum(a->x[1], e, &e);
    c2 = two_sum(a->x[2], e, &e);
    c3 = two_sum(a->x[3], e, &e);
    
    qd_renorm2(&c0, &c1, &c2, &c3, &e);
    
    qd_real c = {c0, c1, c2, c3};
    return c;
}

/* quad-double + double-double */
inline qd_real qd_add (const qd_real *a, const dd_real *b) {
    
    double s0, s1, s2, s3;
    double t0, t1;
    
    s0 = two_sum(a->x[0], b->x[0], &t0);
    s1 = two_sum(a->x[1], b->x[1], &t1);
    
    s1 = two_sum(s1, t0, &t0);
    
    s2 = a->x[2];
    three_sum(&s2, &t0, &t1);
    
    s3 = two_sum(t0, a->x[3], &t0);
    t0 += t1;
    
    qd_renorm2(&s0, &s1, &s2, &s3, &t0);
    qd_real s = {s0, s1, s2, s3};
    return s;
}

//
///* double + quad-double */
//inline qd_real operator+(double a, const qd_real *b) {
//  return (b + a);
//}
//
///* double-double + quad-double */
//inline qd_real operator+(const dd_real *a, const qd_real *b) {
//  return (b + a);
//}

//namespace qd {

/* s = quick_three_accum(a, b, c) adds c to the dd-pair (a, b).
 * If the result does not fit in two doubles, then the sum is
 * output into s and (a,b) contains the remainder.  Otherwise
 * s is zero and (a,b) contains the sum. */
inline double quick_three_accum(double *a, double *b, double c) {
    double s;
    bool za, zb;
    
    s = two_sum(*b, c, b);
    s = two_sum(*a, s, a);
    
    za = (*a != 0.0);
    zb = (*b != 0.0);
    
    if (za && zb)
        return s;
    
    if (!zb) {
        *b = *a;
        *a = s;
    } else {
        *a = s;
    }
    
    return 0.0;
}

//}

inline qd_real qd_ieee_add(const qd_real *a, const qd_real *b) {
    int i, j, k;
    double s, t;
    double u, v;   /* double-length accumulator */
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    
    i = j = k = 0;
    if (fabs(a->x[i]) > fabs(b->x[j]))
        u = a->x[i++];
    else
        u = b->x[j++];
    if (fabs(a->x[i]) > fabs(b->x[j]))
        v = a->x[i++];
    else
        v = b->x[j++];
    
    u = quick_two_sum(u, v, &v);
    
    while (k < 4) {
        if (i >= 4 && j >= 4) {
            x[k] = u;
            if (k < 3)
                x[++k] = v;
            break;
        }
        
        if (i >= 4)
            t = b->x[j++];
        else if (j >= 4)
            t = a->x[i++];
        else if (fabs(a->x[i]) > fabs(b->x[j])) {
            t = a->x[i++];
        } else
            t = b->x[j++];
        
        s = quick_three_accum(&u, &v, t);
        
        if (s != 0.0) {
            x[k++] = s;
        }
    }
    
    /* add the rest. */
    for (k = i; k < 4; k++)
        x[3] += a->x[k];
    for (k = j; k < 4; k++)
        x[3] += b->x[k];
    
    qd_renorm(&x[0], &x[1], &x[2], &x[3]);
    qd_real s2 = {x[0], x[1], x[2], x[3]};
    return s2;
}

inline qd_real qd_sloppy_add(const qd_real *a, const qd_real *b) {
    /*
     double s0, s1, s2, s3;
     double t0, t1, t2, t3;
     
     s0 = two_sum(a[0], b[0], t0);
     s1 = two_sum(a[1], b[1], t1);
     s2 = two_sum(a[2], b[2], t2);
     s3 = two_sum(a[3], b[3], t3);
     
     s1 = two_sum(s1, t0, t0);
     three_sum(s2, t0, t1);
     three_sum2(s3, t0, t2);
     t0 = t0 + t1 + t3;
     
     renorm(s0, s1, s2, s3, t0);
     return qd_real(s0, s1, s2, s3, t0);
     */
    
    /* Same as above, but addition re-organized to minimize
     data dependency ... unfortunately some compilers are
     not very smart to do this automatically */
    double s0, s1, s2, s3;
    double t0, t1, t2, t3;
    
    double v0, v1, v2, v3;
    double u0, u1, u2, u3;
    double w0, w1, w2, w3;
    
    s0 = a->x[0] + b->x[0];
    s1 = a->x[1] + b->x[1];
    s2 = a->x[2] + b->x[2];
    s3 = a->x[3] + b->x[3];
    
    v0 = s0 - a->x[0];
    v1 = s1 - a->x[1];
    v2 = s2 - a->x[2];
    v3 = s3 - a->x[3];
    
    u0 = s0 - v0;
    u1 = s1 - v1;
    u2 = s2 - v2;
    u3 = s3 - v3;
    
    w0 = a->x[0] - u0;
    w1 = a->x[1] - u1;
    w2 = a->x[2] - u2;
    w3 = a->x[3] - u3;
    
    u0 = b->x[0] - v0;
    u1 = b->x[1] - v1;
    u2 = b->x[2] - v2;
    u3 = b->x[3] - v3;
    
    t0 = w0 + u0;
    t1 = w1 + u1;
    t2 = w2 + u2;
    t3 = w3 + u3;
    
    s1 = two_sum(s1, t0, &t0);
    three_sum(&s2, &t0, &t1);
    three_sum2(&s3, &t0, &t2);
    t0 = t0 + t1 + t3;
    
    /* renormalize */
    qd_renorm2(&s0, &s1, &s2, &s3, &t0);
    qd_real s = {s0, s1, s2, s3};
    return s;
}

///* quad-double + quad-double */
//inline qd_real operator+(const qd_real *a, const qd_real *b) {
//#ifndef QD_IEEE_ADD
//  return qd_sloppy_add(a, b);
//#else
//  return qd_ieee_add(a, b);
//#endif
//}
//
//
//
///********** Self-Additions ************/
///* quad-double += double */
//inline qd_real &qd_operator+=(double a) {
//  *this = *this + a;
//  return *this;
//}
//
///* quad-double += double-double */
//inline qd_real &qd_operator+=(const dd_real *a) {
//  *this = *this + a;
//  return *this;
//}
//
///* quad-double += quad-double */
//inline qd_real &qd_operator+=(const qd_real *a) {
//  *this = *this + a;
//  return *this;
//}
//
///********** Unary Minus **********/
inline qd_real qd_negate(const qd_real *a) {
    qd_real b = {-a->x[0], -a->x[1], -a->x[2], -a->x[3]};
    return b;
}

///********** Subtractions **********/
//inline qd_real operator-(const qd_real *a, double b) {
//  return (a + (-b));
//}
//
//inline qd_real operator-(double a, const qd_real *b) {
//  return (a + (-b));
//}
//
//inline qd_real operator-(const qd_real *a, const dd_real *b) {
//  return (a + (-b));
//}
//
//inline qd_real operator-(const dd_real *a, const qd_real *b) {
//  return (a + (-b));
//}
//
//inline qd_real operator-(const qd_real *a, const qd_real *b) {
//  return (a + (-b));
//}
//
///********** Self-Subtractions **********/
//inline qd_real &qd_operator-=(double a) {
//  return ((*this) += (-a));
//}
//
//inline qd_real &qd_operator-=(const dd_real *a) {
//  return ((*this) += (-a));
//}
//
//inline qd_real &qd_operator-=(const qd_real *a) {
//  return ((*this) += (-a));
//}
//
//
//inline qd_real operator*(double a, const qd_real *b) {
//  return (b * a);
//}
//
//inline qd_real operator*(const dd_real *a, const qd_real *b) {
//  return (b * a);
//}
//
//inline qd_real mul_pwr2(const qd_real *a, double b) {
//  return qd_real(a[0] * b, a[1] * b, a[2] * b, a[3] * b);
//}

/********** Multiplications **********/
inline qd_real qd_muldd(const qd_real *a, double b) {
    double p0, p1, p2, p3;
    double q0, q1, q2;
    double s0, s1, s2, s3, s4;
    
    p0 = two_prod(a->x[0], b, &q0);
    p1 = two_prod(a->x[1], b, &q1);
    p2 = two_prod(a->x[2], b, &q2);
    p3 = a->x[3] * b;
    
    s0 = p0;
    
    s1 = two_sum(q0, p1, &s2);
    
    three_sum(&s2, &q1, &p2);
    
    three_sum2(&q1, &q2, &p3);
    s3 = q1;
    
    s4 = q2 + p2;
    
    qd_renorm2(&s0, &s1, &s2, &s3, &s4);
    qd_real s = {s0, s1, s2, s3};
    return s;
}

/* quad-double * double-double */
/* a0 * b0                        0
        a0 * b1                   1
        a1 * b0                   2
             a1 * b1              3
             a2 * b0              4
                  a2 * b1         5
                  a3 * b0         6
                       a3 * b1    7 */
inline qd_real qd_mul (const qd_real *a, const dd_real *b) {
    double p0, p1, p2, p3, p4;
    double q0, q1, q2, q3, q4;
    double s0, s1, s2;
    double t0, t1;
    
    p0 = two_prod(a->x[0], b->x[0], &q0);
    p1 = two_prod(a->x[0], b->x[1], &q1);
    p2 = two_prod(a->x[1], b->x[0], &q2);
    p3 = two_prod(a->x[1], b->x[1], &q3);
    p4 = two_prod(a->x[2], b->x[0], &q4);
    
    three_sum(&p1, &p2, &q0);
    
    /* Five-Three-Sum */
    three_sum(&p2, &p3, &p4);
    q1 = two_sum(q1, q2, &q2);
    s0 = two_sum(p2, q1, &t0);
    s1 = two_sum(p3, q2, &t1);
    s1 = two_sum(s1, t0, &t0);
    s2 = t0 + t1 + p4;
    p2 = s0;
    
    p3 = a->x[2] * b->x[0] + a->x[3] * b->x[1] + q3 + q4;
    three_sum2(&p3, &q0, &s1);
    p4 = q0 + s2;
    
    qd_renorm2(&p0, &p1, &p2, &p3, &p4);
    qd_real p = {p0, p1, p2, p3};
    return p;
}

/* quad-double * quad-double */
/* a0 * b0                    0
        a0 * b1               1
        a1 * b0               2
             a0 * b2          3
             a1 * b1          4
             a2 * b0          5
                  a0 * b3     6
                  a1 * b2     7
                  a2 * b1     8
                  a3 * b0     9  */
inline qd_real qd_sloppy_mul(const qd_real *a, const qd_real *b) {
    double p0, p1, p2, p3, p4, p5;
    double q0, q1, q2, q3, q4, q5;
    double t0, t1;
    double s0, s1, s2;
    
    p0 = two_prod(a->x[0], b->x[0], &q0);
    
    p1 = two_prod(a->x[0], b->x[1], &q1);
    p2 = two_prod(a->x[1], b->x[0], &q2);
    
    p3 = two_prod(a->x[0], b->x[2], &q3);
    p4 = two_prod(a->x[1], b->x[1], &q4);
    p5 = two_prod(a->x[2], b->x[0], &q5);
    
    /* Start Accumulation */
    three_sum(&p1, &p2, &q0);
    
    /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
    three_sum(&p2, &q1, &q2);
    three_sum(&p3, &p4, &p5);
    /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5). */
    s0 = two_sum(p2, p3, &t0);
    s1 = two_sum(q1, p4, &t1);
    s2 = q2 + p5;
    s1 = two_sum(s1, t0, &t0);
    s2 += (t0 + t1);
    
    /* O(eps^3) order terms */
    s1 += a->x[0]*b->x[3] + a->x[1]*b->x[2] + a->x[2]*b->x[1] + a->x[3]*b->x[0] + q0 + q3 + q4 + q5;
    qd_renorm2(&p0, &p1, &s0, &s1, &s2);
    qd_real p = {p0, p1, s0, s1};
    return p;
}

inline qd_real qd_accurate_mul(const qd_real *a, const qd_real *b) {
    double p0, p1, p2, p3, p4, p5;
    double q0, q1, q2, q3, q4, q5;
    double p6, p7, p8, p9;
    double q6, q7, q8, q9;
    double r0, r1;
    double t0, t1;
    double s0, s1, s2;
    
    p0 = two_prod(a->x[0], b->x[0], &q0);
    
    p1 = two_prod(a->x[0], b->x[1], &q1);
    p2 = two_prod(a->x[1], b->x[0], &q2);
    
    p3 = two_prod(a->x[0], b->x[2], &q3);
    p4 = two_prod(a->x[1], b->x[1], &q4);
    p5 = two_prod(a->x[2], b->x[0], &q5);
    
    /* Start Accumulation */
    three_sum(&p1, &p2, &q0);
    
    /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
    three_sum(&p2, &q1, &q2);
    three_sum(&p3, &p4, &p5);
    /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5). */
    s0 = two_sum(p2, p3, &t0);
    s1 = two_sum(q1, p4, &t1);
    s2 = q2 + p5;
    s1 = two_sum(s1, t0, &t0);
    s2 += (t0 + t1);
    
    /* O(eps^3) order terms */
    p6 = two_prod(a->x[0], b->x[3], &q6);
    p7 = two_prod(a->x[1], b->x[2], &q7);
    p8 = two_prod(a->x[2], b->x[1], &q8);
    p9 = two_prod(a->x[3], b->x[0], &q9);
    
    /* Nine-Two-Sum of q0, s1, q3, q4, q5, p6, p7, p8, p9. */
    q0 = two_sum(q0, q3, &q3);
    q4 = two_sum(q4, q5, &q5);
    p6 = two_sum(p6, p7, &p7);
    p8 = two_sum(p8, p9, &p9);
    /* Compute (t0, t1) = (q0, q3) + (q4, q5). */
    t0 = two_sum(q0, q4, &t1);
    t1 += (q3 + q5);
    /* Compute (r0, r1) = (p6, p7) + (p8, p9). */
    r0 = two_sum(p6, p8, &r1);
    r1 += (p7 + p9);
    /* Compute (q3, q4) = (t0, t1) + (r0, r1). */
    q3 = two_sum(t0, r0, &q4);
    q4 += (t1 + r1);
    /* Compute (t0, t1) = (q3, q4) + s1. */
    t0 = two_sum(q3, s1, &t1);
    t1 += q4;
    
    /* O(eps^4) terms -- Nine-One-Sum */
    t1 += a->x[1] * b->x[3] + a->x[2] * b->x[2] + a->x[3] * b->x[1] + q6 + q7 + q8 + q9 + s2;
    
    qd_renorm2(&p0, &p1, &s0, &t0, &t1);
    qd_real p = {p0, p1, s0, t0};
    return p;
}

//inline qd_real operator*(const qd_real *a, const qd_real *b) {
//#ifdef QD_SLOPPY_MUL
//  return qd_sloppy_mul(a, b);
//#else
//  return qd_accurate_mul(a, b);
//#endif
//}

/* quad-double ^ 2  = (x0 + x1 + x2 + x3) ^ 2
                    = x0 ^ 2 + 2 x0 * x1 + (2 x0 * x2 + x1 ^ 2)
                               + (2 x0 * x3 + 2 x1 * x2)           */
inline qd_real qd_sqr(const qd_real *a) {
    double p0, p1, p2, p3, p4, p5;
    double q0, q1, q2, q3;
    double s0, s1;
    double t0, t1;
    
    p0 = two_sqr(a->x[0], &q0);
    p1 = two_prod(2.0 * a->x[0], a->x[1], &q1);
    p2 = two_prod(2.0 * a->x[0], a->x[2], &q2);
    p3 = two_sqr(a->x[1], &q3);
    
    p1 = two_sum(q0, p1, &q0);
    
    q0 = two_sum(q0, q1, &q1);
    p2 = two_sum(p2, p3, &p3);
    
    s0 = two_sum(q0, p2, &t0);
    s1 = two_sum(q1, p3, &t1);
    
    s1 = two_sum(s1, t0, &t0);
    t0 += t1;
    
    s1 = quick_two_sum(s1, t0, &t0);
    p2 = quick_two_sum(s0, s1, &t1);
    p3 = quick_two_sum(t1, t0, &q0);
    
    p4 = 2.0 * a->x[0] * a->x[3];
    p5 = 2.0 * a->x[1] * a->x[2];
    
    p4 = two_sum(p4, p5, &p5);
    q2 = two_sum(q2, q3, &q3);
    
    t0 = two_sum(p4, q2, &t1);
    t1 = t1 + p5 + q3;
    
    p3 = two_sum(p3, t0, &p4);
    p4 = p4 + q0 + t1;
    
    qd_renorm2(&p0, &p1, &p2, &p3, &p4);
    qd_real p = {p0, p1, p2, p3};
    return p;
}

///********** Self-Multiplication **********/
///* quad-double *= double */
//inline qd_real &qd_operator*=(double a) {
//  *this = (*this * a);
//  return *this;
//}
//
///* quad-double *= double-double */
//inline qd_real &qd_operator*=(const dd_real *a) {
//  *this = (*this * a);
//  return *this;
//}
//
///* quad-double *= quad-double */
//inline qd_real &qd_operator*=(const qd_real *a) {
//  *this = *this * a;
//  return *this;
//}
//
//inline qd_real operator/ (const qd_real *a, const dd_real *b) {
//#ifdef QD_SLOPPY_DIV
//  return qd_sloppy_div(a, b);
//#else
//  return qd_accurate_div(a, b);
//#endif
//}
//
//inline qd_real operator/(const qd_real *a, const qd_real *b) {
//#ifdef QD_SLOPPY_DIV
//  return qd_sloppy_div(a, b);
//#else
//  return qd_accurate_div(a, b);
//#endif
//}

///* double / quad-double */
//inline qd_real operator/(double a, const qd_real *b) {
//  return qd_real(a) / b;
//}
//
///* double-double / quad-double */
//inline qd_real operator/(const dd_real *a, const qd_real *b) {
//  return qd_real(a) / b;
//}
//
///********** Self-Divisions **********/
///* quad-double /= double */
//inline qd_real &qd_operator/=(double a) {
//  *this = (*this / a);
//  return *this;
//}
//
///* quad-double /= double-double */
//inline qd_real &qd_operator/=(const dd_real *a) {
//  *this = (*this / a);
//  return *this;
//}
//
///* quad-double /= quad-double */
//inline qd_real &qd_operator/=(const qd_real *a) {
//  *this = (*this / a);
//  return *this;
//}
//
//
///********** Exponentiation **********/
//inline qd_real qd_operator^(int n) const {
//  return pow(*this, n);
//}
//
///********** Miscellaneous **********/
//inline qd_real abs(const qd_real *a) {
//  return (a[0] < 0.0) ? -a : a;
//}
//
//inline qd_real fabs(const qd_real *a) {
//  return abs(a);
//}
//
///* Quick version.  May be off by one when qd is very close
//   to the middle of two integers.                         */
//inline qd_real quick_nint(const qd_real *a) {
//  qd_real r = qd_real(qd::nint(a[0]), qd::nint(a[1]),
//      qd::nint(a[2]), qd::nint(a[3]));
//  r.renorm();
//  return r;
//}
//
///*********** Assignments ************/
///* quad-double = double */
//inline qd_real &qd_operator=(double a) {
//  x[0] = a;
//  x[1] = x[2] = x[3] = 0.0;
//  return *this;
//}
//
///* quad-double = double-double */
//inline qd_real &qd_operator=(const dd_real *a) {
//  x[0] = a._hi();
//  x[1] = a._lo();
//  x[2] = x[3] = 0.0;
//  return *this;
//}
//
/////********** Equality Comparison **********/
//inline bool operator==(const qd_real *a, double b) {
//  return (a[0] == b && a[1] == 0.0 && a[2] == 0.0 && a[3] == 0.0);
//}
//
//inline bool operator==(double a, const qd_real *b) {
//  return (b == a);
//}
//
//inline bool operator==(const qd_real *a, const dd_real *b) {
//  return (a[0] == b._hi() && a[1] == b._lo() &&
//          a[2] == 0.0 && a[3] == 0.0);
//}
//
//inline bool operator==(const dd_real *a, const qd_real *b) {
//  return (b == a);
//}
//
//inline bool operator==(const qd_real &a, const qd_real &b) {
//  return (a[0] == b[0] && a[1] == b[1] &&
//          a[2] == b[2] && a[3] == b[3]);
//}
//
//
///********** Less-Than Comparison ***********/
//inline bool operator<(const qd_real &a, double b) {
//  return (a[0] < b || (a[0] == b && a[1] < 0.0));
//}
//
//inline bool operator<(double a, const qd_real &b) {
//  return (b > a);
//}
//
//inline bool operator<(const qd_real &a, const dd_real &b) {
//  return (a[0] < b._hi() ||
//          (a[0] == b._hi() && (a[1] < b._lo() ||
//                            (a[1] == b._lo() && a[2] < 0.0))));
//}
//
//inline bool operator<(const dd_real &a, const qd_real &b) {
//  return (b > a);
//}
//
//inline bool operator<(const qd_real &a, const qd_real &b) {
//  return (a[0] < b[0] ||
//          (a[0] == b[0] && (a[1] < b[1] ||
//                            (a[1] == b[1] && (a[2] < b[2] ||
//                                              (a[2] == b[2] && a[3] < b[3]))))));
//}
//
///********** Greater-Than Comparison ***********/
//inline bool operator>(const qd_real &a, double b) {
//  return (a[0] > b || (a[0] == b && a[1] > 0.0));
//}
//
//inline bool operator>(double a, const qd_real &b) {
//  return (b < a);
//}
//
//inline bool operator>(const qd_real &a, const dd_real &b) {
//  return (a[0] > b._hi() ||
//          (a[0] == b._hi() && (a[1] > b._lo() ||
//                            (a[1] == b._lo() && a[2] > 0.0))));
//}
//
//inline bool operator>(const dd_real &a, const qd_real &b) {
//  return (b < a);
//}
//
//inline bool operator>(const qd_real &a, const qd_real &b) {
//  return (a[0] > b[0] ||
//          (a[0] == b[0] && (a[1] > b[1] ||
//                            (a[1] == b[1] && (a[2] > b[2] ||
//                                              (a[2] == b[2] && a[3] > b[3]))))));
//}
//
//
///********** Less-Than-Or-Equal-To Comparison **********/
//inline bool operator<=(const qd_real &a, double b) {
//  return (a[0] < b || (a[0] == b && a[1] <= 0.0));
//}
//
//inline bool operator<=(double a, const qd_real &b) {
//  return (b >= a);
//}
//
//inline bool operator<=(const qd_real &a, const dd_real &b) {
//  return (a[0] < b._hi() ||
//          (a[0] == b._hi() && (a[1] < b._lo() ||
//                            (a[1] == b._lo() && a[2] <= 0.0))));
//}
//
//inline bool operator<=(const dd_real &a, const qd_real &b) {
//  return (b >= a);
//}
//
//inline bool operator<=(const qd_real &a, const qd_real &b) {
//  return (a[0] < b[0] ||
//          (a[0] == b[0] && (a[1] < b[1] ||
//                            (a[1] == b[1] && (a[2] < b[2] ||
//                                              (a[2] == b[2] && a[3] <= b[3]))))));
//}
//
///********** Greater-Than-Or-Equal-To Comparison **********/
//inline bool operator>=(const qd_real &a, double b) {
//  return (a[0] > b || (a[0] == b && a[1] >= 0.0));
//}
//
//inline bool operator>=(double a, const qd_real &b) {
//  return (b <= a);
//}
//
//inline bool operator>=(const qd_real &a, const dd_real &b) {
//  return (a[0] > b._hi() ||
//          (a[0] == b._hi() && (a[1] > b._lo() ||
//                            (a[1] == b._lo() && a[2] >= 0.0))));
//}
//
//inline bool operator>=(const dd_real &a, const qd_real &b) {
//  return (b <= a);
//}
//
//inline bool operator>=(const qd_real &a, const qd_real &b) {
//  return (a[0] > b[0] ||
//          (a[0] == b[0] && (a[1] > b[1] ||
//                            (a[1] == b[1] && (a[2] > b[2] ||
//                                              (a[2] == b[2] && a[3] >= b[3]))))));
//}
//
//
//
///********** Not-Equal-To Comparison **********/
//inline bool operator!=(const qd_real &a, double b) {
//  return !(a == b);
//}
//
//inline bool operator!=(double a, const qd_real &b) {
//  return !(a == b);
//}
//
//inline bool operator!=(const qd_real &a, const dd_real &b) {
//  return !(a == b);
//}
//
//inline bool operator!=(const dd_real &a, const qd_real &b) {
//  return !(a == b);
//}
//
//inline bool operator!=(const qd_real &a, const qd_real &b) {
//  return !(a == b);
//}
//
//
//
//inline qd_real aint(const qd_real &a) {
//  return (a[0] >= 0) ? floor(a) : ceil(a);
//}
//
inline bool qd_is_zero(const qd_real *a) {
  return (a->x[0] == 0.0);
}

inline bool qd_is_one(const qd_real *a) {
  return (a->x[0] == 1.0 && a->x[1] == 0.0 && a->x[2] == 0.0 && a->x[3] == 0.0);
}

inline bool qd_is_positive(const qd_real *a) {
  return (a->x[0] > 0.0);
}

inline bool qd_is_negative(const qd_real *a) {
  return (a->x[0] < 0.0);
}

//inline dd_real to_dd_real(const qd_real &a) {
//  return dd_real(a[0], a[1]);
//}
//
//inline double to_double(const qd_real &a) {
//  return a[0];
//}
//
//inline int to_int(const qd_real &a) {
//  return static_cast<int>(a[0]);
//}
//
//inline qd_real inv(const qd_real &qd) {
//  return 1.0 / qd;
//}
//
//inline qd_real max(const qd_real &a, const qd_real &b) {
//  return (a > b) ? a : b;
//}
//
//inline qd_real max(const qd_real &a, const qd_real &b,
//                   const qd_real &c) {
//  return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
//}
//
//inline qd_real min(const qd_real &a, const qd_real &b) {
//  return (a < b) ? a : b;
//}
//
//inline qd_real min(const qd_real &a, const qd_real &b,
//                   const qd_real &c) {
//  return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
//}
//
///* Random number generator */
//inline qd_real qd_rand() {
//  return qdrand();
//}
//
//inline qd_real ldexp(const qd_real &a, int n) {
//  return qd_real(std::ldexp(a[0], n), std::ldexp(a[1], n),
//                 std::ldexp(a[2], n), std::ldexp(a[3], n));
//}
//
//  qd_real();
//  qd_real(const char *s);
//  qd_real(const dd_real &dd);
//  qd_real(double d);
//  qd_real(int i);
//
//  double operator[](int i) const;
//  double &operator[](int i);
//
//  static void error(const char *msg);
//
//  bool isnan() const;
//  bool isfinite() const { return QD_ISFINITE(x[0]); }
//  bool isinf() const { return QD_ISINF(x[0]); }

//  static qd_real qd_ieee_add(const qd_real *a, const qd_real *b);
//  static qd_real qd_sloppy_add(const qd_real *a, const qd_real *b);

//  qd_real &operator+=(double a);
//  qd_real &operator+=(const dd_real &a);
//  qd_real &operator+=(const qd_real &a);
//
//  qd_real &operator-=(double a);
//  qd_real &operator-=(const dd_real &a);
//  qd_real &operator-=(const qd_real &a);

//  static qd_real qd_sloppy_mul(const qd_real *a, const qd_real *b);
//  static qd_real qd_accurate_mul(const qd_real *a, const qd_real *b);

//  qd_real &operator*=(double a);
//  qd_real &operator*=(const dd_real &a);
//  qd_real &operator*=(const qd_real &a);

//  static qd_real qd_sloppy_divd(const qd_real *a, const dd_real *b);
//  static qd_real qd_accurate_divd(const qd_real *a, const dd_real *b);
//  static qd_real qd_sloppy_div(const qd_real *a, const qd_real *b);
//  static qd_real qd_accurate_div(const qd_real *a, const qd_real *b);

//  qd_real &operator/=(double a);
//  qd_real &operator/=(const dd_real &a);
//  qd_real &operator/=(const qd_real &a);
//
//  qd_real operator^(int n) const;
//
//  qd_real operator-() const;
//
//  qd_real &operator=(double a);
//  qd_real &operator=(const dd_real &a);
//  qd_real &operator=(const char *s);

//  bool qd_is_zero(const qd_real *a);
//  bool qd_is_one(const qd_real *a);
//  bool qd_is_positive(const qd_real *a);
//  bool qd_is_negative(const qd_real *a);
//
//  static qd_real rand(void);

//  void to_digits(char *s, int &expn, int precision = _ndigits) const;
//  void write(char *s, int len, int precision = _ndigits,
//      bool showpos = false, bool uppercase = false) const;
//  std::string to_string(int precision = _ndigits, int width = 0,
//      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0),
//      bool showpos = false, bool uppercase = false, char fill = ' ') const;
//  static int read(const char *s, qd_real &a);

  /* Debugging methods */
//  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
//  void dump_bits(const std::string &name = "",
//                 std::ostream &os = std::cerr) const;
//
//  static qd_real debug_rand();

//};
//
//namespace std {
//  template <>
//  class numeric_limits<qd_real> : public numeric_limits<double> {
//  public:
//    inline static double epsilon() { return qd_real::_eps; }
//    inline static double min() { return qd_real::_min_normalized; }
//    inline static qd_real max() { return qd_real::_max; }
//    inline static qd_real safe_max() { return qd_real::_safe_max; }
//    static const int digits = 209;
//    static const int digits10 = 62;
//  };
//}

QD_API qd_real qd_polyeval(const qd_real *c, int n, const qd_real *x);
QD_API qd_real qd_polyroot(const qd_real *c, int n,
    const qd_real *x0, int max_iter, double thresh);

QD_API qd_real qd_qdrand(void);
QD_API qd_real qd_sqrt(const qd_real *a);

//QD_API inline bool isnan(const qd_real &a) { return a.isnan(); }
//QD_API inline bool isfinite(const qd_real &a) { return a.isfinite(); }
//QD_API inline bool isinf(const qd_real &a) { return a.isinf(); }

/* Computes  qd * d  where d is known to be a power of 2.
   This can be done component wise.                      */
QD_API qd_real qd_mul_pwr2(const qd_real *qd, double d);

//QD_API qd_real operator+(const qd_real *a, const qd_real *b);
//QD_API qd_real operator+(const dd_real *a, const qd_real *b);
//QD_API qd_real operator+(const qd_real *a, const dd_real *b);
//QD_API qd_real operator+(const qd_real *a, double b);
//QD_API qd_real operator+(double a, const qd_real *b);
//
//QD_API qd_real operator-(const qd_real &a, const qd_real &b);
//QD_API qd_real operator-(const dd_real &a, const qd_real &b);
//QD_API qd_real operator-(const qd_real &a, const dd_real &b);
//QD_API qd_real operator-(const qd_real &a, double b);
//QD_API qd_real operator-(double a, const qd_real &b);
//
//QD_API qd_real operator*(const qd_real &a, const qd_real &b);
//QD_API qd_real operator*(const dd_real &a, const qd_real &b);
//QD_API qd_real operator*(const qd_real &a, const dd_real &b);
//QD_API qd_real operator*(const qd_real &a, double b);
//QD_API qd_real operator*(double a, const qd_real &b);
//
//QD_API qd_real operator/(const qd_real &a, const qd_real &b);
//QD_API qd_real operator/(const dd_real &a, const qd_real &b);
//QD_API qd_real operator/(const qd_real &a, const dd_real &b);
//QD_API qd_real operator/(const qd_real &a, double b);
//QD_API qd_real operator/(double a, const qd_real &b);

//QD_API qd_real qd_sqr(const qd_real *a);
QD_API qd_real qd_sqrt(const qd_real *a);
QD_API qd_real qd_powi(const qd_real *a, int n);
QD_API qd_real qd_pow(const qd_real *a, const qd_real *b);
QD_API qd_real qd_npwr(const qd_real *a, int n);

QD_API qd_real qd_nroot(const qd_real *a, int n);

QD_API qd_real qd_rem(const qd_real *a, const qd_real *b);
QD_API qd_real qd_drem(const qd_real *a, const qd_real *b);
QD_API qd_real qd_divrem(const qd_real *a, const qd_real *b, qd_real *r);

//dd_real to_dd_real(const qd_real &a);
//double  to_double(const qd_real &a);
//int     to_int(const qd_real &a);
//
//QD_API bool operator==(const qd_real &a, const qd_real &b);
//QD_API bool operator==(const qd_real &a, const dd_real &b);
//QD_API bool operator==(const dd_real &a, const qd_real &b);
//QD_API bool operator==(double a, const qd_real &b);
//QD_API bool operator==(const qd_real &a, double b);
//
//QD_API bool operator<(const qd_real &a, const qd_real &b);
//QD_API bool operator<(const qd_real &a, const dd_real &b);
//QD_API bool operator<(const dd_real &a, const qd_real &b);
//QD_API bool operator<(double a, const qd_real &b);
//QD_API bool operator<(const qd_real &a, double b);
//
//QD_API bool operator>(const qd_real &a, const qd_real &b);
//QD_API bool operator>(const qd_real &a, const dd_real &b);
//QD_API bool operator>(const dd_real &a, const qd_real &b);
//QD_API bool operator>(double a, const qd_real &b);
//QD_API bool operator>(const qd_real &a, double b);
//
//QD_API bool operator<=(const qd_real &a, const qd_real &b);
//QD_API bool operator<=(const qd_real &a, const dd_real &b);
//QD_API bool operator<=(const dd_real &a, const qd_real &b);
//QD_API bool operator<=(double a, const qd_real &b);
//QD_API bool operator<=(const qd_real &a, double b);

//QD_API bool operator>=(const qd_real &a, const qd_real &b);
//QD_API bool operator>=(const qd_real &a, const dd_real &b);
//QD_API bool operator>=(const dd_real &a, const qd_real &b);
//QD_API bool operator>=(double a, const qd_real &b);
//QD_API bool operator>=(const qd_real &a, double b);
//
//QD_API bool operator!=(const qd_real &a, const qd_real &b);
//QD_API bool operator!=(const qd_real &a, const dd_real &b);
//QD_API bool operator!=(const dd_real &a, const qd_real *b);
//QD_API bool operator!=(double a, const qd_real &b);
//QD_API bool operator!=(const qd_real &a, double b);

QD_API qd_real qd_fabs(const qd_real *a);
QD_API qd_real qd_abs(const qd_real *a);    /* same as fabs */

QD_API qd_real qd_ldexp(const qd_real *a, int n);

QD_API qd_real qd_nint(const qd_real *a);
QD_API qd_real qd_quick_nint(const qd_real *a);
QD_API qd_real qd_floor(const qd_real *a);
QD_API qd_real qd_ceil(const qd_real *a);
QD_API qd_real qd_aint(const qd_real *a);

QD_API qd_real qd_sin(const qd_real *a);
QD_API qd_real qd_cos(const qd_real *a);
QD_API qd_real qd_tan(const qd_real *a);
QD_API void qd_sincos(const qd_real *a, qd_real *s, qd_real *c);

QD_API qd_real qd_asin(const qd_real *a);
QD_API qd_real qd_acos(const qd_real *a);
QD_API qd_real qd_atan(const qd_real *a);
QD_API qd_real qd_atan2(const qd_real *y, const qd_real *x);

QD_API qd_real qd_exp(const qd_real *a);
QD_API qd_real qd_log(const qd_real *a);
QD_API qd_real qd_log10(const qd_real *a);

QD_API qd_real qd_sinh(const qd_real *a);
QD_API qd_real qd_cosh(const qd_real *a);
QD_API qd_real qd_tanh(const qd_real *a);
QD_API void qd_sincosh(const qd_real *a, qd_real *sin_qd, qd_real *cos_qd);

QD_API qd_real qd_asinh(const qd_real *a);
QD_API qd_real qd_acosh(const qd_real *a);
QD_API qd_real qd_atanh(const qd_real *a);

QD_API qd_real qdrand(void);

//QD_API qd_real max(const qd_real *a, const qd_real *b);
//QD_API qd_real max(const qd_real *a, const qd_real *b, const qd_real *c);
//QD_API qd_real min(const qd_real *a, const qd_real *b);
//QD_API qd_real min(const qd_real *a, const qd_real *b, const qd_real *c);

QD_API qd_real qd_fmod(const qd_real *a, const qd_real *b);

//QD_API std::ostream &operator<<(std::ostream &s, const qd_real &a);
//QD_API std::istream &operator>>(std::istream &s, qd_real &a);
//#ifdef QD_INLINE
//#include "qd_inline.h"
//#endif

#endif /* _QD_QD_REAL_H */

