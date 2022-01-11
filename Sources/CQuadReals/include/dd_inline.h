/*
 * include/dd_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains small functions (suitable for inlining) in the double-double
 * arithmetic package.
 */
#ifndef _QD_DD_INLINE_H
#define _QD_DD_INLINE_H

//#include <cmath>
#include "inline.h"

#ifndef QD_INLINE
#define inline
#endif

/*********** Additions ************/
/* double-double = double + double */
inline dd_real dd_add(double a, double b) {
    double s, e;
    s = two_sum(a, b, &e);
    dd_real x = {s, e};
    return x;
}

/* double-double + double */
inline dd_real dd_addd(const dd_real *a, double b) {
    double s1, s2;
    s1 = two_sum(a->x[0], b, &s2);
    s2 += a->x[1];
    s1 =  quick_two_sum(s1, s2, &s2);
    dd_real x = {s1, s2};
    return x;
}

/* double-double + double-double */
inline dd_real dd_ieee_add(const dd_real *a, const dd_real *b) {
    /* This one satisfies IEEE style error bound,
     due to K. Briggs and W. Kahan.                   */
    double s1, s2, t1, t2;
    
    s1 = two_sum(a->x[0], b->x[0], &s2);
    t1 = two_sum(a->x[1], b->x[1], &t2);
    s2 += t1;
    s1 = quick_two_sum(s1, s2, &s2);
    s2 += t2;
    s1 = quick_two_sum(s1, s2, &s2);
    dd_real x = {s1, s2};
    return x;
}

inline dd_real dd_sloppy_add(const dd_real *a, const dd_real *b) {
    /* This is the less accurate version ... obeys Cray-style
     error bound. */
    double s, e;
    
    s = two_sum(a->x[0], b->x[0], &e);
    e += (a->x[1] + b->x[1]);
    s = quick_two_sum(s, e, &e);
    dd_real x = {s, e};
    return x;
}

//inline dd_real dd_add(const dd_real *a, const dd_real *b) {
//#ifndef QD_IEEE_ADD
//  return dd_sloppy_add(a, b);
//#else
//  return dd_ieee_add(a, b);
//#endif
//}

/* double + double-double */
//inline dd_real operator+(double a, const dd_real *b) {
//  return (b + a);
//}


/*********** Self-Additions ************/
/* double-double += double */
//inline void dd_plusEqual(dd_real *x, double a) {
//  double s1, s2;
//  s1 = two_sum(*x[0], a, &s2);
//  s2 += *x[1];
//  x[0] = quick_two_sum(s1, s2, &x[1]);
//}

/* double-double += double-double */
//inline dd_real dd_plusEqual(const dd_real *a) {
//#ifndef QD_IEEE_ADD
//  double s, e;
//  s = two_sum(x[0], a->x[0], e);
//  e += x[1];
//  e += a->x[1];
//  x[0] = quick_two_sum(s, e, x[1]);
//  return *this;
//#else
//  double s1, s2, t1, t2;
//  s1 = two_sum(x[0], a->x[0], s2);
//  t1 = two_sum(x[1], a->x[1], t2);
//  s2 += t1;
//  s1 = quick_two_sum(s1, s2, s2);
//  s2 += t2;
//  x[0] = quick_two_sum(s1, s2, x[1]);
//  return *this;
//#endif
//}

//inline dd_real dd_sloppy_sub(const dd_real *a, const dd_real *b) {
//    double s, e;
//    s = two_diff(a->x[0], b->x[0], &e);
//    e += a->x[1];
//    e -= b->x[1];
//    s = quick_two_sum(s, e, &e);
//    dd_real x = {s, e};
//    return x;
//}

inline dd_real dd_sqrd(double a) {
    double p1, p2;
    p1 = two_sqr(a, &p2);
    dd_real res = {p1, p2};
    return res;
}

/*********** Squaring **********/
inline dd_real dd_sqr(const dd_real *a) {
    double p1, p2;
    double s1, s2;
    p1 = two_sqr(a->x[0], &p2);
    p2 += 2.0 * a->x[0] * a->x[1];
    p2 += a->x[1] * a->x[1];
    s1 = quick_two_sum(p1, p2, &s2);
    dd_real x = {s1, s2};
    return x;
}

/*********** Subtractions ************/
/* double-double = double - double */
inline dd_real dd_sub(double a, double b) {
    double s, e;
    s = two_diff(a, b, &e);
    dd_real x = {s, e};
    return x;
}

/* double-double - double */
inline dd_real dd_subd(const dd_real *a, double b) {
    double s1, s2;
    s1 = two_diff(a->x[0], b, &s2);
    s2 += a->x[1];
    s1 = quick_two_sum(s1, s2, &s2);
    dd_real x = {s1, s2};
    return x;
}

/* double-double - double-double */
/* double-double + double-double */
inline dd_real dd_ieee_sub(const dd_real *a, const dd_real *b) {
    double s1, s2, t1, t2;
    s1 = two_diff(a->x[0], b->x[0], &s2);
    t1 = two_diff(a->x[1], b->x[1], &t2);
    s2 += t1;
    s1 = quick_two_sum(s1, s2, &s2);
    s2 += t2;
    s1 = quick_two_sum(s1, s2, &s2);
    dd_real x = {s1, s2};
    return x;
}

inline dd_real dd_sloppy_sub(const dd_real *a, const dd_real *b) {
    double s, e;
    s = two_diff(a->x[0], b->x[0], &e);
    e += a->x[1];
    e -= b->x[1];
    s = quick_two_sum(s, e, &e);
    dd_real x = {s, e};
    return x;
}

/* double - double-double */
inline dd_real dd_subd2(double a, const dd_real *b) {
    double s1, s2;
    s1 = two_diff(a, b->x[0], &s2);
    s2 -= b->x[1];
    s1 = quick_two_sum(s1, s2, &s2);
    dd_real x = {s1, s2};
    return x;
}

/*********** Self-Subtractions ************/
/* double-double -= double */
//inline dd_real &dd_operator-=(double a) {
//  double s1, s2;
//  s1 = two_diff(x[0], a, s2);
//  s2 += x[1];
//  x[0] = quick_two_sum(s1, s2, x[1]);
//  return *this;
//}
//
///* double-double -= double-double */
inline void minusEqual(dd_real *x, const dd_real *a) {
  double s1, s2, t1, t2;
  s1 = two_diff(x->x[0], a->x[0], &s2);
  t1 = two_diff(x->x[1], a->x[1], &t2);
  s2 += t1;
  s1 = quick_two_sum(s1, s2, &s2);
  s2 += t2;
  x->x[0] = quick_two_sum(s1, s2, &x->x[1]);
}

/*********** Unary Minus ***********/
//inline dd_real dd_operator-() const {
//  return dd_real(-x[0], -x[1]);
//}

/*********** Multiplications ************/
/* double-double = double * double */
inline dd_real dd_mul(double a, double b) {
    double p, e;
    p = two_prod(a, b, &e);
    dd_real x = {p, e};
    return x;
}

/* double-double * (2.0 ^ exp) */
//inline dd_real dd_ldexp(const dd_real *a, int exp) {
//  return dd_real(std::ldexp(a->x[0], exp), std::ldexp(a->x[1], exp));
//}
//
///* double-double * double,  where double is a power of 2. */
//inline dd_real dd_mul_pwr2(const dd_real *a, double b) {
//  return dd_real(a->x[0] * b, a->x[1] * b);
//}

/* double-double * double */
inline dd_real dd_muld(const dd_real *a, double b) {
    double p1, p2;
    
    p1 = two_prod(a->x[0], b, &p2);
    p2 += (a->x[1] * b);
    p1 = quick_two_sum(p1, p2, &p2);
    dd_real x = {p1, p2};
    return x;
}

/* double-double * double-double */
inline dd_real dd_muldd(const dd_real *a, const dd_real *b) {
    double p1, p2;
    
    p1 = two_prod(a->x[0], b->x[0], &p2);
    p2 += (a->x[0] * b->x[1] + a->x[1] * b->x[0]);
    p1 = quick_two_sum(p1, p2, &p2);
    dd_real x = {p1, p2};
    return x;
}

/* double * double-double */
//inline dd_real operator*(double a, const dd_real *b) {
//  return (b * a);
//}

/*********** Self-Multiplications ************/
/* double-double *= double */
//inline dd_real &dd_operator*=(double a) {
//  double p1, p2;
//  p1 = two_prod(x[0], a, p2);
//  p2 += x[1] * a;
//  x[0] = quick_two_sum(p1, p2, x[1]);
//  return *this;
//}
//
///* double-double *= double-double */
//inline dd_real &dd_operator*=(const dd_real *a) {
//  double p1, p2;
//  p1 = two_prod(x[0], a->x[0], p2);
//  p2 += a->x[1] * x[0];
//  p2 += a->x[0] * x[1];
//  x[0] = quick_two_sum(p1, p2, x[1]);
//  return *this;
//}

/*********** Divisions ************/
inline dd_real dd_div(double a, double b) {
    double q1, q2;
    double p1, p2;
    double s, e;
    
    q1 = a / b;
    
    /* Compute  a - q1 * b */
    p1 = two_prod(q1, b, &p2);
    s = two_diff(a, p1, &e);
    e -= p2;
    
    /* get next approximation */
    q2 = (s + e) / b;
    
    s = quick_two_sum(q1, q2, &e);
    
    dd_real x = {s, e};
    return x;
}

/* double-double / double */
inline dd_real dd_divd(const dd_real *a, double b) {

  double q1, q2;
  double p1, p2;
  double s, e;
  dd_real r;
  
  q1 = a->x[0] / b;   /* approximate quotient. */

  /* Compute  this - q1 * d */
  p1 = two_prod(q1, b, &p2);
  s = two_diff(a->x[0], p1, &e);
  e += a->x[1];
  e -= p2;
  
  /* get next approximation. */
  q2 = (s + e) / b;

  /* renormalize */
  r.x[0] = quick_two_sum(q1, q2, &r.x[1]);

  return r;
}

inline dd_real dd_sloppy_div(const dd_real *a, const dd_real *b) {
    double s1, s2;
    double q1, q2;
    dd_real r;
    
    q1 = a->x[0] / b->x[0];  /* approximate quotient */
    
    /* compute  this - q1 * dd */
    r = dd_muld(b, q1);
    //  r = b * q1;
    s1 = two_diff(a->x[0], r.x[0], &s2);
    s2 -= r.x[1];
    s2 += a->x[1];
    
    /* get next approximation */
    q2 = (s1 + s2) / b->x[0];
    
    /* renormalize */
    r.x[0] = quick_two_sum(q1, q2, &r.x[1]);
    return r;
}

inline dd_real dd_accurate_div(const dd_real *a, const dd_real *b) {
    double q1, q2, q3;
    dd_real r;
    
    q1 = a->x[0] / b->x[0];  /* approximate quotient */
    
    r = dd_muld(b, q1);
    r = dd_ieee_sub(a, &r);
    //  r = a - q1 * b;
    
    q2 = r.x[0] / b->x[0];
    dd_real t = dd_muld(b, q2);
    minusEqual(&r, &t);
    // r -= (q2 * b);
    
    q3 = r.x[0] / b->x[0];
    
    q1 = quick_two_sum(q1, q2, &q2);
    r.x[0] = q1; r.x[1] = q2;
    r = dd_addd(&r, q3);
    return r;
}

///* double-double / double-double */
//inline dd_real operator/(const dd_real *a, const dd_real *b) {
//#ifdef QD_SLOPPY_DIV
//  return dd_sloppy_div(a, b);
//#else
//  return dd_accurate_div(a, b);
//#endif
//}

///* double / double-double */
//inline dd_real operator/(double a, const dd_real *b) {
//  return dd_real(a) / b;
//}
//
//inline dd_real inv(const dd_real *a) {
//  return 1.0 / a;
//}

/*********** Self-Divisions ************/
/* double-double /= double */
//inline dd_real &dd_operator/=(double a) {
//  *this = *this / a;
//  return *this;
//}
//
///* double-double /= double-double */
//inline dd_real &dd_operator/=(const dd_real *a) {
//  *this = *this / a;
//  return *this;
//}
//
///********** Remainder **********/
//inline dd_real drem(const dd_real *a, const dd_real *b) {
//  dd_real n = nint(a / b);
//  return (a - n * b);
//}
//
//inline dd_real divrem(const dd_real *a, const dd_real *b, dd_real *r) {
//  dd_real n = nint(a / b);
//  r = a - n * b;
//  return n;
//}


///********** Exponentiation **********/
//inline dd_real dd_operator^(int n) {
//  return npwr(*this, n);
//}
//
//
///*********** Assignments ************/
///* double-double = double */
//inline dd_real &dd_operator=(double a) {
//  x[0] = a;
//  x[1] = 0.0;
//  return *this;
//}

/*********** Equality Comparisons ************/
/* double-double == double */
//inline bool operator==(const dd_real *a, double b) {
//  return (a->x[0] == b && a->x[1] == 0.0);
//}
//
///* double-double == double-double */
//inline bool operator==(const dd_real *a, const dd_real *b) {
//  return (a->x[0] == b->x[0] && a->x[1] == b->x[1]);
//}
//
///* double == double-double */
//inline bool operator==(double a, const dd_real *b) {
//  return (a == b->x[0] && b->x[1] == 0.0);
//}
//
///*********** Greater-Than Comparisons ************/
///* double-double > double */
//inline bool operator>(const dd_real *a, double b) {
//  return (a->x[0] > b || (a->x[0] == b && a->x[1] > 0.0));
//}
//
///* double-double > double-double */
//inline bool operator>(const dd_real *a, const dd_real *b) {
//  return (a->x[0] > b->x[0] || (a->x[0] == b->x[0] && a->x[1] > b->x[1]));
//}
//
///* double > double-double */
//inline bool operator>(double a, const dd_real *b) {
//  return (a > b->x[0] || (a == b->x[0] && b->x[1] < 0.0));
//}

/*********** Less-Than Comparisons ************/
/* double-double < double */
//inline bool operator<(const dd_real *a, double b) {
//  return (a->x[0] < b || (a->x[0] == b && a->x[1] < 0.0));
//}
//
///* double-double < double-double */
//inline bool operator<(const dd_real *a, const dd_real *b) {
//  return (a->x[0] < b->x[0] || (a->x[0] == b->x[0] && a->x[1] < b->x[1]));
//}
//
///* double < double-double */
//inline bool operator<(double a, const dd_real *b) {
//  return (a < b->x[0] || (a == b->x[0] && b->x[1] > 0.0));
//}
//
///*********** Greater-Than-Or-Equal-To Comparisons ************/
///* double-double >= double */
//inline bool operator>=(const dd_real *a, double b) {
//  return (a->x[0] > b || (a->x[0] == b && a->x[1] >= 0.0));
//}
//
///* double-double >= double-double */
//inline bool operator>=(const dd_real *a, const dd_real *b) {
//  return (a->x[0] > b->x[0] || (a->x[0] == b->x[0] && a->x[1] >= b->x[1]));
//}
//
///* double >= double-double */
//inline bool operator>=(double a, const dd_real *b) {
//  return (b <= a);
//}

/*********** Less-Than-Or-Equal-To Comparisons ************/
/* double-double <= double */
//inline bool operator<=(const dd_real *a, double b) {
//  return (a->x[0] < b || (a->x[0] == b && a->x[1] <= 0.0));
//}
//
///* double-double <= double-double */
//inline bool operator<=(const dd_real *a, const dd_real *b) {
//  return (a->x[0] < b->x[0] || (a->x[0] == b->x[0] && a->x[1] <= b->x[1]));
//}
//
///* double <= double-double */
//inline bool operator<=(double a, const dd_real *b) {
//  return (b >= a);
//}
//
///*********** Not-Equal-To Comparisons ************/
///* double-double != double */
//inline bool operator!=(const dd_real *a, double b) {
//  return (a->x[0] != b || a->x[1] != 0.0);
//}
//
///* double-double != double-double */
//inline bool operator!=(const dd_real *a, const dd_real *b) {
//  return (a->x[0] != b->x[0] || a->x[1] != b->x[1]);
//}
//
///* double != double-double */
//inline bool operator!=(double a, const dd_real *b) {
//  return (a != b->x[0] || b->x[1] != 0.0);
//}
//


///* Absolute value */
//inline dd_real abs(const dd_real *a) {
//  return (a->x[0] < 0.0) ? -a : a;
//}
//
//inline dd_real fabs(const dd_real *a) {
//  return abs(a);
//}
//
///* Round to Nearest integer */
//inline dd_real nint(const dd_real *a) {
//  double hi = nint(a->x[0]);
//  double lo;
//
//  if (hi == a->x[0]) {
//    /* High word is an integer already.  Round the low word.*/
//    lo = nint(a->x[1]);
//
//    /* Renormalize. This is needed if x[0] = some integer, x[1] = 1/2.*/
//    hi = quick_two_sum(hi, lo, lo);
//  } else {
//    /* High word is not an integer. */
//    lo = 0.0;
//    if (std::abs(hi-a->x[0]) == 0.5 && a->x[1] < 0.0) {
//      /* There is a tie in the high word, consult the low word
//         to break the tie. */
//      hi -= 1.0;      /* NOTE: This does not cause INEXACT. */
//    }
//  }
//
//  return dd_real(hi, lo);
//}
//
//inline dd_real floor(const dd_real *a) {
//  double hi = std::floor(a->x[0]);
//  double lo = 0.0;
//
//  if (hi == a->x[0]) {
//    /* High word is integer already.  Round the low word. */
//    lo = std::floor(a->x[1]);
//    hi = quick_two_sum(hi, lo, lo);
//  }
//
//  return dd_real(hi, lo);
//}
//
//inline dd_real ceil(const dd_real *a) {
//  double hi = std::ceil(a->x[0]);
//  double lo = 0.0;
//
//  if (hi == a->x[0]) {
//    /* High word is integer already.  Round the low word. */
//    lo = std::ceil(a->x[1]);
//    hi = quick_two_sum(hi, lo, lo);
//  }
//
//  return dd_real(hi, lo);
//}
//
//inline dd_real aint(const dd_real *a) {
//  return (a->x[0] >= 0.0) ? floor(a) : ceil(a);
//}
//
///* Cast to double. */
//inline double to_double(const dd_real *a) {
//  return a->x[0];
//}
//
///* Cast to int. */
//inline int to_int(const dd_real *a) {
//  return static_cast<int>(a->x[0]);
//}

//**

#endif /* _QD_DD_INLINE_H */
