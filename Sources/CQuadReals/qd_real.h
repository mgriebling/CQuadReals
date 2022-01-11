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

  static qd_real qd_ieee_add(const qd_real *a, const qd_real *b);
  static qd_real qd_sloppy_add(const qd_real *a, const qd_real *b);

//  qd_real &operator+=(double a);
//  qd_real &operator+=(const dd_real &a);
//  qd_real &operator+=(const qd_real &a);
//
//  qd_real &operator-=(double a);
//  qd_real &operator-=(const dd_real &a);
//  qd_real &operator-=(const qd_real &a);

  static qd_real qd_sloppy_mul(const qd_real *a, const qd_real *b);
  static qd_real qd_accurate_mul(const qd_real *a, const qd_real *b);

//  qd_real &operator*=(double a);
//  qd_real &operator*=(const dd_real &a);
//  qd_real &operator*=(const qd_real &a);

  static qd_real qd_sloppy_divd(const qd_real *a, const dd_real *b);
  static qd_real qd_accurate_divd(const qd_real *a, const dd_real *b);
  static qd_real qd_sloppy_div(const qd_real *a, const qd_real *b);
  static qd_real qd_accurate_div(const qd_real *a, const qd_real *b);

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
#ifdef QD_INLINE
#include "qd_inline.h"
#endif

#endif /* _QD_QD_REAL_H */

