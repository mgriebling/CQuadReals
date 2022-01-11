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

typedef struct QD_API dd_real {
    double x[2];
} dd_real;

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

  static dd_real dd_add(double a, double b);
  static dd_real dd_addd(const dd_real *a, double b);
  static dd_real dd_ieee_add(const dd_real *a, const dd_real *b);
  static dd_real dd_sloppy_add(const dd_real *a, const dd_real *b);

//  dd_real &operator+=(double a);
//  dd_real &operator+=(const dd_real &a);

  static dd_real dd_sub(double a, double b);
  static dd_real dd_subd2(double a, const dd_real *b);
  static dd_real dd_ieee_sub(const dd_real *a, const dd_real *b);
  static dd_real dd_sloppy_sub(const dd_real *a, const dd_real *b);

//  dd_real &operator-=(double a);
//  dd_real &operator-=(const dd_real &a);
//
//  dd_real operator-() const;

  static dd_real dd_mul(double a, double b);

//  dd_real &operator*=(double a);
//  dd_real &operator*=(const dd_real &a);

  static dd_real dd_div(double a, double b);
  static dd_real dd_sloppy_div(const dd_real *a, const dd_real *b);
  static dd_real dd_accurate_div(const dd_real *a, const dd_real *b);
  
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
#ifdef QD_INLINE
#include "dd_inline.h"
#endif

#endif /* _QD_DD_REAL_H */

