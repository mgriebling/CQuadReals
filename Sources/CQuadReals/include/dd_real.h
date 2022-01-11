/*
 * include/dd_real.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Double-double precision (>= 106-bit significand) floating point
 * arithmetic package based on David Bailey's Fortran-90 double-double
 * package, with some changes. See  
 *
 *   http://www.nersc.gov/~dhbailey/mpdist/mpdist.html
 *   
 * for the original Fortran-90 version.
 *
 * Overall structure is similar to that of Keith Brigg's C++ double-double
 * package.  See  
 *
 *   http://www-epidem.plansci.cam.ac.uk/~kbriggs/doubledouble.html
 *
 * for more details.  In particular, the fix for x86 computers is borrowed
 * from his code.
 *
 * Yozo Hida
 */

#ifndef _QD_DD_REAL_H
#define _QD_DD_REAL_H

//#include <cmath>
//#include <iostream>
//#include <string>
//#include <limits>
#include "qd_config.h"
#include "inline.h"
//#include <qd/fpu.h>

// Some compilers define isnan, isfinite, and isinf as macros, even for
// C++ codes, which cause havoc when overloading these functions.  We undef
// them here.
#ifdef isnan
#undef isnan
#endif

#ifdef isfinite
#undef isfinite
#endif

#ifdef isinf
#undef isinf
#endif

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

#ifndef bool
typedef unsigned char bool;
#define true 1
#define false 0
#endif

typedef struct {
    double x[2];
} dd_real;

/*********** Additions ************/
/* double-double = double + double */
static inline dd_real dd_add(double a, double b) {
    double s, e;
    s = two_sum(a, b, &e);
    dd_real x = {s, e};
    return x;
}

/* double-double + double */
static inline dd_real dd_addd(const dd_real *a, double b) {
    double s1, s2;
    s1 = two_sum(a->x[0], b, &s2);
    s2 += a->x[1];
    s1 =  quick_two_sum(s1, s2, &s2);
    dd_real x = {s1, s2};
    return x;
}

/* double-double + double-double */
static inline dd_real dd_ieee_add(const dd_real *a, const dd_real *b) {
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

static inline dd_real dd_sloppy_add(const dd_real *a, const dd_real *b) {
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

static inline dd_real dd_sqrd(double a) {
    double p1, p2;
    p1 = two_sqr(a, &p2);
    dd_real res = {p1, p2};
    return res;
}

/*********** Squaring **********/
static inline dd_real dd_sqr(const dd_real *a) {
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
static inline dd_real dd_sub(double a, double b) {
    double s, e;
    s = two_diff(a, b, &e);
    dd_real x = {s, e};
    return x;
}

/* double-double - double */
static inline dd_real dd_subd(const dd_real *a, double b) {
    double s1, s2;
    s1 = two_diff(a->x[0], b, &s2);
    s2 += a->x[1];
    s1 = quick_two_sum(s1, s2, &s2);
    dd_real x = {s1, s2};
    return x;
}

/* double-double - double-double */
/* double-double + double-double */
static inline dd_real dd_ieee_sub(const dd_real *a, const dd_real *b) {
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

static inline dd_real dd_sloppy_sub(const dd_real *a, const dd_real *b) {
    double s, e;
    s = two_diff(a->x[0], b->x[0], &e);
    e += a->x[1];
    e -= b->x[1];
    s = quick_two_sum(s, e, &e);
    dd_real x = {s, e};
    return x;
}

/* double - double-double */
static inline dd_real dd_subd2(double a, const dd_real *b) {
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
static inline void minusEqual(dd_real *x, const dd_real *a) {
  double s1, s2, t1, t2;
  s1 = two_diff(x->x[0], a->x[0], &s2);
  t1 = two_diff(x->x[1], a->x[1], &t2);
  s2 += t1;
  s1 = quick_two_sum(s1, s2, &s2);
  s2 += t2;
  x->x[0] = quick_two_sum(s1, s2, &x->x[1]);
}

/*********** Unary Minus ***********/
//static inline dd_real dd_operator-() const {
//  return dd_real(-x[0], -x[1]);
//}

/*********** Multiplications ************/
/* double-double = double * double */
static inline dd_real dd_mul(double a, double b) {
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
static inline dd_real dd_muld(const dd_real *a, double b) {
    double p1, p2;
    
    p1 = two_prod(a->x[0], b, &p2);
    p2 += (a->x[1] * b);
    p1 = quick_two_sum(p1, p2, &p2);
    dd_real x = {p1, p2};
    return x;
}

/* double-double * double-double */
static inline dd_real dd_muldd(const dd_real *a, const dd_real *b) {
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
static inline dd_real dd_div(double a, double b) {
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
static inline dd_real dd_divd(const dd_real *a, double b) {

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

static inline dd_real dd_sloppy_div(const dd_real *a, const dd_real *b) {
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

static inline dd_real dd_accurate_div(const dd_real *a, const dd_real *b) {
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

//  dd_real(double hi, double lo) { x[0] = hi; x[1] = lo; }
//  dd_real() {x[0] = 0.0; x[1] = 0.0; }
//  dd_real(double h) { x[0] = h; x[1] = 0.0; }
//  dd_real(int h) {
//    x[0] = (static_cast<double>(h));
//    x[1] = 0.0;
//  }
//
//  dd_real (const char *s);
//  explicit dd_real (const double *d) {
//    x[0] = d[0]; x[1] = d[1];
//  }

//  static void error(const char *msg);

//  double _hi() { return x[0]; }
//  double _lo() { return x[1]; }

  static const dd_real dd_2pi;
  static const dd_real dd_pi;
  static const dd_real dd_3pi4;
  static const dd_real dd_pi2;
  static const dd_real dd_pi4;
  static const dd_real dd_e;
  static const dd_real dd_log2;
  static const dd_real dd_log10c;
  static const dd_real dd_nan;
  static const dd_real dd_inf;
  static const dd_real dd_zero;

  static const double dd_eps;
  static const double dd_min_normalized;
  static const dd_real dd_max;
  static const dd_real dd_safe_max;
  static const int dd_ndigits;

//  bool isnan() { return QD_ISNAN(x[0]) || QD_ISNAN(x[1]); }
//  bool isfinite() { return QD_ISFINITE(x[0]); }
//  bool isinf() { return QD_ISINF(x[0]); }

//  static dd_real dd_add(double a, double b);
//  static dd_real dd_addd(const dd_real *a, double b);
//  static dd_real dd_ieee_add(const dd_real *a, const dd_real *b);
//  static dd_real dd_sloppy_add(const dd_real *a, const dd_real *b);

//  dd_real &operator+=(double a);
//  dd_real &operator+=(const dd_real &a);

//  static dd_real dd_sub(double a, double b);
//  static dd_real dd_subd2(double a, const dd_real *b);
//  static dd_real dd_ieee_sub(const dd_real *a, const dd_real *b);
//  static dd_real dd_sloppy_sub(const dd_real *a, const dd_real *b);

//  dd_real &operator-=(double a);
//  dd_real &operator-=(const dd_real &a);
//
//  dd_real operator-() const;

//  static dd_real dd_mul(double a, double b);

//  dd_real &operator*=(double a);
//  dd_real &operator*=(const dd_real &a);

//  static dd_real dd_div(double a, double b);
//  static dd_real dd_sloppy_div(const dd_real *a, const dd_real *b);
//  static dd_real dd_accurate_div(const dd_real *a, const dd_real *b);
  
//  dd_real &operator/=(double a);
//  dd_real &operator/=(const dd_real *a);
//
//  dd_real &operator=(double a);
//  dd_real &operator=(const char *s);
//
//  dd_real operator^(int n);
//QD_API dd_real sqrt(const dd_real *a);
  
//  bool dd_is_zero(void);
//  bool dd_is_one(void);
//  bool dd_is_positive(void);
//  bool dd_is_negative(void);
///*********** Micellaneous ************/
///*  this == 0 */
//inline bool dd_is_zero(const dd_real *a) {
//  return (a->x[0] == 0.0);
//}
//
///*  this == 1 */
//inline bool dd_is_one(const dd_real *a) {
//  return (a->x[0] == 1.0 && a->x[1] == 0.0);
//}
//
///*  this > 0 */
//inline bool dd_is_positive(const dd_real *a) {
//  return (a->x[0] > 0.0);
//}
//
///* this < 0 */
//inline bool dd_is_negative(const dd_real *a) {
//  return (a->x[0] < 0.0);
//}

//  static dd_real dd_rand(void);

//  void to_digits(char *s, int *expn, int precision);
//  void write(char *s, int len, int precision,
//      bool showpos, bool uppercase );
//  std::string to_string(int precision = _ndigits, int width = 0,
//      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0),
//      bool showpos = false, bool uppercase = false, char fill = ' ');
//  int read(const char *s, dd_real *a);

  /* Debugging Methods */
//  void dump(const std::string &name = "", std::ostream &os = std::cerr);
//  void dump_bits(const std::string &name = "",
//                 std::ostream &os = std::cerr);

//  static dd_real debug_rand(void);


//namespace std {
//  template <>
//  class numeric_limits<dd_real> : public numeric_limits<double> {
//  public:
//    inline static double epsilon() { return dd_real::_eps; }
//    inline static dd_real max() { return dd_real::_max; }
//    inline static dd_real safe_max() { return dd_real::_safe_max; }
//    inline static double min() { return dd_real::_min_normalized; }
//    static const int digits = 104;
//    static const int digits10 = 31;
//  };
//}

QD_API dd_real ddrand(void);
QD_API dd_real dd_sqrt(const dd_real *a);

QD_API dd_real dd_polyeval(const dd_real *c, int n, const dd_real *x);
QD_API dd_real dd_polyroot(const dd_real *c, int n,
    const dd_real *x0, int max_iter, double thresh);

//QD_API inline bool isnan(const dd_real *a) { return a->isnan(); }
//QD_API inline bool isfinite(const dd_real *a) { return a->isfinite(); }
//QD_API inline bool isinf(const dd_real *a) { return a->isinf(); }

/* Computes  dd * d  where d is known to be a power of 2. */
QD_API dd_real dd_mul_pwr2(const dd_real *dd, double d);

//QD_API dd_real operator+(const dd_real *a, double b);
//QD_API dd_real operator+(double a, const dd_real *b);
//QD_API dd_real operator+(const dd_real *a, const dd_real *b);
//
//QD_API dd_real operator-(const dd_real *a, double b);
//QD_API dd_real operator-(double a, const dd_real *b);
//QD_API dd_real operator-(const dd_real *a, const dd_real *b);
//
//QD_API dd_real operator*(const dd_real *a, double b);
//QD_API dd_real operator*(double a, const dd_real *b);
//QD_API dd_real operator*(const dd_real *a, const dd_real *b);
//
//QD_API dd_real operator/(const dd_real *a, double b);
//QD_API dd_real operator/(double a, const dd_real *b);
//QD_API dd_real operator/(const dd_real *a, const dd_real *b);

//QD_API dd_real inv(const dd_real *a);

QD_API dd_real dd_rem(const dd_real *a, const dd_real *b);
QD_API dd_real dd_drem(const dd_real *a, const dd_real *b);
QD_API dd_real dd_divrem(const dd_real *a, const dd_real *b, dd_real *r);

QD_API dd_real dd_powi(const dd_real *a, int n);
QD_API dd_real dd_pow(const dd_real *a, const dd_real *b);
QD_API dd_real dd_npwr(const dd_real *a, int n);
//QD_API dd_real dd_sqr(const dd_real *a);

QD_API dd_real dd_sqrt(const dd_real *a);
QD_API dd_real dd_nroot(const dd_real *a, int n);

//QD_API bool operator==(const dd_real *a, double b);
//QD_API bool operator==(double a, const dd_real *b);
//QD_API bool operator==(const dd_real *a, const dd_real *b);
//
//QD_API bool operator<=(const dd_real *a, double b);
//QD_API bool operator<=(double a, const dd_real *b);
//QD_API bool operator<=(const dd_real *a, const dd_real *b);
//
//QD_API bool operator>=(const dd_real *a, double b);
//QD_API bool operator>=(double a, const dd_real *b);
//QD_API bool operator>=(const dd_real *a, const dd_real *b);
//
//QD_API bool operator<(const dd_real *a, double b);
//QD_API bool operator<(double a, const dd_real *b);
//QD_API bool operator<(const dd_real *a, const dd_real *b);
//
//QD_API bool operator>(const dd_real *a, double b);
//QD_API bool operator>(double a, const dd_real *b);
//QD_API bool operator>(const dd_real *a, const dd_real *b);
//
//QD_API bool operator!=(const dd_real *a, double b);
//QD_API bool operator!=(double a, const dd_real *b);
//QD_API bool operator!=(const dd_real *a, const dd_real *b);

QD_API dd_real dd_nint(const dd_real *a);
QD_API dd_real dd_floor(const dd_real *a);
QD_API dd_real dd_ceil(const dd_real *a);
QD_API dd_real dd_aint(const dd_real *a);

QD_API dd_real dd_ddrand(void);

double to_double(const dd_real *a);
int    to_int(const dd_real *a);

QD_API dd_real dd_exp(const dd_real *a);
QD_API dd_real dd_ldexp(const dd_real *a, int exp);
QD_API dd_real dd_log(const dd_real *a);
QD_API dd_real dd_log10(const dd_real *a);

QD_API dd_real dd_sin(const dd_real *a);
QD_API dd_real dd_cos(const dd_real *a);
QD_API dd_real dd_tan(const dd_real *a);
QD_API void sincos(const dd_real *a, dd_real *sin_a, dd_real *cos_a);

QD_API dd_real dd_asin(const dd_real *a);
QD_API dd_real dd_acos(const dd_real *a);
QD_API dd_real dd_atan(const dd_real *a);
QD_API dd_real dd_atan2(const dd_real *y, const dd_real *x);

QD_API dd_real dd_sinh(const dd_real *a);
QD_API dd_real dd_cosh(const dd_real *a);
QD_API dd_real dd_tanh(const dd_real *a);
QD_API void sincosh(const dd_real *a,
                      dd_real *sinh_a, dd_real *cosh_a);

QD_API dd_real dd_asinh(const dd_real *a);
QD_API dd_real dd_acosh(const dd_real *a);
QD_API dd_real dd_atanh(const dd_real *a);

QD_API dd_real dd_fabs(const dd_real *a);
QD_API dd_real dd_abs(const dd_real *a);   /* same as fabs */

QD_API dd_real dd_fmod(const dd_real *a, const dd_real *b);

//QD_API std::ostream& operator<<(std::ostream &s, const dd_real &a);
//QD_API std::istream& operator>>(std::istream &s, dd_real &a);
//#ifdef QD_INLINE
//#include "dd_inline.h"
//#endif

#endif /* _QD_DD_REAL_H */

