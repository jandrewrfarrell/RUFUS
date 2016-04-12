/* DO NOT EDIT! GENERATED AUTOMATICALLY! */
/* A GNU-like <math.h>.

   Copyright (C) 2002-2003, 2007-2011 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

#ifndef _GL_MATH_H

#if __GNUC__ >= 3
#pragma GCC system_header
#endif


/* The include_next requires a split double-inclusion guard.  */
#include_next <math.h>

#ifndef _GL_MATH_H
#define _GL_MATH_H


/* The definitions of _GL_FUNCDECL_RPL etc. are copied here.  */
#ifndef _GL_CXXDEFS_H
#define _GL_CXXDEFS_H

/* The three most frequent use cases of these macros are:

   * For providing a substitute for a function that is missing on some
     platforms, but is declared and works fine on the platforms on which
     it exists:

       #if @GNULIB_FOO@
       # if !@HAVE_FOO@
       _GL_FUNCDECL_SYS (foo, ...);
       # endif
       _GL_CXXALIAS_SYS (foo, ...);
       _GL_CXXALIASWARN (foo);
       #elif defined GNULIB_POSIXCHECK
       ...
       #endif

   * For providing a replacement for a function that exists on all platforms,
     but is broken/insufficient and needs to be replaced on some platforms:

       #if @GNULIB_FOO@
       # if @REPLACE_FOO@
       #  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
       #   undef foo
       #   define foo rpl_foo
       #  endif
       _GL_FUNCDECL_RPL (foo, ...);
       _GL_CXXALIAS_RPL (foo, ...);
       # else
       _GL_CXXALIAS_SYS (foo, ...);
       # endif
       _GL_CXXALIASWARN (foo);
       #elif defined GNULIB_POSIXCHECK
       ...
       #endif

   * For providing a replacement for a function that exists on some platforms
     but is broken/insufficient and needs to be replaced on some of them and
     is additionally either missing or undeclared on some other platforms:

       #if @GNULIB_FOO@
       # if @REPLACE_FOO@
       #  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
       #   undef foo
       #   define foo rpl_foo
       #  endif
       _GL_FUNCDECL_RPL (foo, ...);
       _GL_CXXALIAS_RPL (foo, ...);
       # else
       #  if !@HAVE_FOO@   or   if !@HAVE_DECL_FOO@
       _GL_FUNCDECL_SYS (foo, ...);
       #  endif
       _GL_CXXALIAS_SYS (foo, ...);
       # endif
       _GL_CXXALIASWARN (foo);
       #elif defined GNULIB_POSIXCHECK
       ...
       #endif
*/

/* _GL_EXTERN_C declaration;
   performs the declaration with C linkage.  */
#if defined __cplusplus
# define _GL_EXTERN_C extern "C"
#else
# define _GL_EXTERN_C extern
#endif

/* _GL_FUNCDECL_RPL (func, rettype, parameters_and_attributes);
   declares a replacement function, named rpl_func, with the given prototype,
   consisting of return type, parameters, and attributes.
   Example:
     _GL_FUNCDECL_RPL (open, int, (const char *filename, int flags, ...)
                                  _GL_ARG_NONNULL ((1)));
 */
#define _GL_FUNCDECL_RPL(func,rettype,parameters_and_attributes) \
  _GL_FUNCDECL_RPL_1 (rpl_##func, rettype, parameters_and_attributes)
#define _GL_FUNCDECL_RPL_1(rpl_func,rettype,parameters_and_attributes) \
  _GL_EXTERN_C rettype rpl_func parameters_and_attributes

/* _GL_FUNCDECL_SYS (func, rettype, parameters_and_attributes);
   declares the system function, named func, with the given prototype,
   consisting of return type, parameters, and attributes.
   Example:
     _GL_FUNCDECL_SYS (open, int, (const char *filename, int flags, ...)
                                  _GL_ARG_NONNULL ((1)));
 */
#define _GL_FUNCDECL_SYS(func,rettype,parameters_and_attributes) \
  _GL_EXTERN_C rettype func parameters_and_attributes

/* _GL_CXXALIAS_RPL (func, rettype, parameters);
   declares a C++ alias called GNULIB_NAMESPACE::func
   that redirects to rpl_func, if GNULIB_NAMESPACE is defined.
   Example:
     _GL_CXXALIAS_RPL (open, int, (const char *filename, int flags, ...));
 */
#define _GL_CXXALIAS_RPL(func,rettype,parameters) \
  _GL_CXXALIAS_RPL_1 (func, rpl_##func, rettype, parameters)
#if defined __cplusplus && defined GNULIB_NAMESPACE
# define _GL_CXXALIAS_RPL_1(func,rpl_func,rettype,parameters) \
    namespace GNULIB_NAMESPACE                                \
    {                                                         \
      rettype (*const func) parameters = ::rpl_func;          \
    }                                                         \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#else
# define _GL_CXXALIAS_RPL_1(func,rpl_func,rettype,parameters) \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#endif

/* _GL_CXXALIAS_RPL_CAST_1 (func, rpl_func, rettype, parameters);
   is like  _GL_CXXALIAS_RPL_1 (func, rpl_func, rettype, parameters);
   except that the C function rpl_func may have a slightly different
   declaration.  A cast is used to silence the "invalid conversion" error
   that would otherwise occur.  */
#if defined __cplusplus && defined GNULIB_NAMESPACE
# define _GL_CXXALIAS_RPL_CAST_1(func,rpl_func,rettype,parameters) \
    namespace GNULIB_NAMESPACE                                     \
    {                                                              \
      rettype (*const func) parameters =                           \
        reinterpret_cast<rettype(*)parameters>(::rpl_func);        \
    }                                                              \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#else
# define _GL_CXXALIAS_RPL_CAST_1(func,rpl_func,rettype,parameters) \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#endif

/* _GL_CXXALIAS_SYS (func, rettype, parameters);
   declares a C++ alias called GNULIB_NAMESPACE::func
   that redirects to the system provided function func, if GNULIB_NAMESPACE
   is defined.
   Example:
     _GL_CXXALIAS_SYS (open, int, (const char *filename, int flags, ...));
 */
#if defined __cplusplus && defined GNULIB_NAMESPACE
  /* If we were to write
       rettype (*const func) parameters = ::func;
     like above in _GL_CXXALIAS_RPL_1, the compiler could optimize calls
     better (remove an indirection through a 'static' pointer variable),
     but then the _GL_CXXALIASWARN macro below would cause a warning not only
     for uses of ::func but also for uses of GNULIB_NAMESPACE::func.  */
# define _GL_CXXALIAS_SYS(func,rettype,parameters) \
    namespace GNULIB_NAMESPACE                     \
    {                                              \
      static rettype (*func) parameters = ::func;  \
    }                                              \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#else
# define _GL_CXXALIAS_SYS(func,rettype,parameters) \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#endif

/* _GL_CXXALIAS_SYS_CAST (func, rettype, parameters);
   is like  _GL_CXXALIAS_SYS (func, rettype, parameters);
   except that the C function func may have a slightly different declaration.
   A cast is used to silence the "invalid conversion" error that would
   otherwise occur.  */
#if defined __cplusplus && defined GNULIB_NAMESPACE
# define _GL_CXXALIAS_SYS_CAST(func,rettype,parameters) \
    namespace GNULIB_NAMESPACE                          \
    {                                                   \
      static rettype (*func) parameters =               \
        reinterpret_cast<rettype(*)parameters>(::func); \
    }                                                   \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#else
# define _GL_CXXALIAS_SYS_CAST(func,rettype,parameters) \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#endif

/* _GL_CXXALIAS_SYS_CAST2 (func, rettype, parameters, rettype2, parameters2);
   is like  _GL_CXXALIAS_SYS (func, rettype, parameters);
   except that the C function is picked among a set of overloaded functions,
   namely the one with rettype2 and parameters2.  Two consecutive casts
   are used to silence the "cannot find a match" and "invalid conversion"
   errors that would otherwise occur.  */
#if defined __cplusplus && defined GNULIB_NAMESPACE
  /* The outer cast must be a reinterpret_cast.
     The inner cast: When the function is defined as a set of overloaded
     functions, it works as a static_cast<>, choosing the designated variant.
     When the function is defined as a single variant, it works as a
     reinterpret_cast<>. The parenthesized cast syntax works both ways.  */
# define _GL_CXXALIAS_SYS_CAST2(func,rettype,parameters,rettype2,parameters2) \
    namespace GNULIB_NAMESPACE                                                \
    {                                                                         \
      static rettype (*func) parameters =                                     \
        reinterpret_cast<rettype(*)parameters>(                               \
          (rettype2(*)parameters2)(::func));                                  \
    }                                                                         \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#else
# define _GL_CXXALIAS_SYS_CAST2(func,rettype,parameters,rettype2,parameters2) \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#endif

/* _GL_CXXALIASWARN (func);
   causes a warning to be emitted when ::func is used but not when
   GNULIB_NAMESPACE::func is used.  func must be defined without overloaded
   variants.  */
#if defined __cplusplus && defined GNULIB_NAMESPACE
# define _GL_CXXALIASWARN(func) \
   _GL_CXXALIASWARN_1 (func, GNULIB_NAMESPACE)
# define _GL_CXXALIASWARN_1(func,namespace) \
   _GL_CXXALIASWARN_2 (func, namespace)
/* To work around GCC bug <http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43881>,
   we enable the warning only when not optimizing.  */
# if !__OPTIMIZE__
#  define _GL_CXXALIASWARN_2(func,namespace) \
    _GL_WARN_ON_USE (func, \
                     "The symbol ::" #func " refers to the system function. " \
                     "Use " #namespace "::" #func " instead.")
# elif __GNUC__ >= 3 && GNULIB_STRICT_CHECKING
#  define _GL_CXXALIASWARN_2(func,namespace) \
     extern __typeof__ (func) func
# else
#  define _GL_CXXALIASWARN_2(func,namespace) \
     _GL_EXTERN_C int _gl_cxxalias_dummy
# endif
#else
# define _GL_CXXALIASWARN(func) \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#endif

/* _GL_CXXALIASWARN1 (func, rettype, parameters_and_attributes);
   causes a warning to be emitted when the given overloaded variant of ::func
   is used but not when GNULIB_NAMESPACE::func is used.  */
#if defined __cplusplus && defined GNULIB_NAMESPACE
# define _GL_CXXALIASWARN1(func,rettype,parameters_and_attributes) \
   _GL_CXXALIASWARN1_1 (func, rettype, parameters_and_attributes, \
                        GNULIB_NAMESPACE)
# define _GL_CXXALIASWARN1_1(func,rettype,parameters_and_attributes,namespace) \
   _GL_CXXALIASWARN1_2 (func, rettype, parameters_and_attributes, namespace)
/* To work around GCC bug <http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43881>,
   we enable the warning only when not optimizing.  */
# if !__OPTIMIZE__
#  define _GL_CXXALIASWARN1_2(func,rettype,parameters_and_attributes,namespace) \
    _GL_WARN_ON_USE_CXX (func, rettype, parameters_and_attributes, \
                         "The symbol ::" #func " refers to the system function. " \
                         "Use " #namespace "::" #func " instead.")
# elif __GNUC__ >= 3 && GNULIB_STRICT_CHECKING
#  define _GL_CXXALIASWARN1_2(func,rettype,parameters_and_attributes,namespace) \
     extern __typeof__ (func) func
# else
#  define _GL_CXXALIASWARN1_2(func,rettype,parameters_and_attributes,namespace) \
     _GL_EXTERN_C int _gl_cxxalias_dummy
# endif
#else
# define _GL_CXXALIASWARN1(func,rettype,parameters_and_attributes) \
    _GL_EXTERN_C int _gl_cxxalias_dummy
#endif

#endif /* _GL_CXXDEFS_H */

/* The definition of _GL_ARG_NONNULL is copied here.  */
/* _GL_ARG_NONNULL((n,...,m)) tells the compiler and static analyzer tools
   that the values passed as arguments n, ..., m must be non-NULL pointers.
   n = 1 stands for the first argument, n = 2 for the second argument etc.  */
#ifndef _GL_ARG_NONNULL
# if (__GNUC__ == 3 && __GNUC_MINOR__ >= 3) || __GNUC__ > 3
#  define _GL_ARG_NONNULL(params) __attribute__ ((__nonnull__ params))
# else
#  define _GL_ARG_NONNULL(params)
# endif
#endif

/* The definition of _GL_WARN_ON_USE is copied here.  */
#ifndef _GL_WARN_ON_USE

# if 4 < __GNUC__ || (__GNUC__ == 4 && 3 <= __GNUC_MINOR__)
/* A compiler attribute is available in gcc versions 4.3.0 and later.  */
#  define _GL_WARN_ON_USE(function, message) \
extern __typeof__ (function) function __attribute__ ((__warning__ (message)))
# elif __GNUC__ >= 3 && GNULIB_STRICT_CHECKING
/* Verify the existence of the function.  */
#  define _GL_WARN_ON_USE(function, message) \
extern __typeof__ (function) function
# else /* Unsupported.  */
#  define _GL_WARN_ON_USE(function, message) \
_GL_WARN_EXTERN_C int _gl_warn_on_use
# endif
#endif

/* _GL_WARN_ON_USE_CXX (function, rettype, parameters_and_attributes, "string")
   is like _GL_WARN_ON_USE (function, "string"), except that the function is
   declared with the given prototype, consisting of return type, parameters,
   and attributes.
   This variant is useful for overloaded functions in C++. _GL_WARN_ON_USE does
   not work in this case.  */
#ifndef _GL_WARN_ON_USE_CXX
# if 4 < __GNUC__ || (__GNUC__ == 4 && 3 <= __GNUC_MINOR__)
#  define _GL_WARN_ON_USE_CXX(function,rettype,parameters_and_attributes,msg) \
extern rettype function parameters_and_attributes \
     __attribute__ ((__warning__ (msg)))
# elif __GNUC__ >= 3 && GNULIB_STRICT_CHECKING
/* Verify the existence of the function.  */
#  define _GL_WARN_ON_USE_CXX(function,rettype,parameters_and_attributes,msg) \
extern rettype function parameters_and_attributes
# else /* Unsupported.  */
#  define _GL_WARN_ON_USE_CXX(function,rettype,parameters_and_attributes,msg) \
_GL_WARN_EXTERN_C int _gl_warn_on_use
# endif
#endif

/* _GL_WARN_EXTERN_C declaration;
   performs the declaration with C linkage.  */
#ifndef _GL_WARN_EXTERN_C
# if defined __cplusplus
#  define _GL_WARN_EXTERN_C extern "C"
# else
#  define _GL_WARN_EXTERN_C extern
# endif
#endif

/* Helper macros to define a portability warning for the
   classification macro FUNC called with VALUE.  POSIX declares the
   classification macros with an argument of real-floating (that is,
   one of float, double, or long double).  */
#define _GL_WARN_REAL_FLOATING_DECL(func) \
static inline int                                                   \
rpl_ ## func ## f (float f)                                         \
{                                                                   \
  return func (f);                                                  \
}                                                                   \
static inline int                                                   \
rpl_ ## func ## d (double d)                                        \
{                                                                   \
  return func (d);                                                  \
}                                                                   \
static inline int                                                   \
rpl_ ## func ## l (long double l)                                   \
{                                                                   \
  return func (l);                                                  \
}                                                                   \
_GL_WARN_ON_USE (rpl_ ## func ## f, #func " is unportable - "       \
                 "use gnulib module " #func " for portability");    \
_GL_WARN_ON_USE (rpl_ ## func ## d, #func " is unportable - "       \
                 "use gnulib module " #func " for portability");    \
_GL_WARN_ON_USE (rpl_ ## func ## l, #func " is unportable - "       \
                 "use gnulib module " #func " for portability")
#define _GL_WARN_REAL_FLOATING_IMPL(func, value) \
  (sizeof (value) == sizeof (float) ? rpl_ ## func ## f (value)     \
   : sizeof (value) == sizeof (double) ? rpl_ ## func ## d (value)  \
   : rpl_ ## func ## l (value))


/* POSIX allows platforms that don't support NAN.  But all major
   machines in the past 15 years have supported something close to
   IEEE NaN, so we define this unconditionally.  We also must define
   it on platforms like Solaris 10, where NAN is present but defined
   as a function pointer rather than a floating point constant.  */
#if !defined NAN || 0
# undef NAN
  /* The Compaq (ex-DEC) C 6.4 compiler chokes on the expression 0.0 / 0.0.  */
# ifdef __DECC
static float
_NaN ()
{
  static float zero = 0.0f;
  return zero / zero;
}
#  define NAN (_NaN())
# else
#  define NAN (0.0f / 0.0f)
# endif
#endif

/* Solaris 10 defines HUGE_VAL, but as a function pointer rather
   than a floating point constant.  */
#if 0
# undef HUGE_VAL
# define HUGE_VAL (1.0 / 0.0)
#endif


/* Write x as
     x = mantissa * 2^exp
   where
     If x finite and nonzero: 0.5 <= |mantissa| < 1.0.
     If x is zero: mantissa = x, exp = 0.
     If x is infinite or NaN: mantissa = x, exp unspecified.
   Store exp in *EXPPTR and return mantissa.  */
#if 1
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   define frexp rpl_frexp
#  endif
_GL_FUNCDECL_RPL (frexp, double, (double x, int *expptr) _GL_ARG_NONNULL ((2)));
_GL_CXXALIAS_RPL (frexp, double, (double x, int *expptr));
# else
_GL_CXXALIAS_SYS (frexp, double, (double x, int *expptr));
# endif
_GL_CXXALIASWARN (frexp);
#elif defined GNULIB_POSIXCHECK
# undef frexp
/* Assume frexp is always declared.  */
_GL_WARN_ON_USE (frexp, "frexp is unportable - "
                 "use gnulib module frexp for portability");
#endif


#if 0
# if !1
_GL_EXTERN_C double logb (double x);
# endif
#elif defined GNULIB_POSIXCHECK
# undef logb
# if HAVE_RAW_DECL_LOGB
_GL_WARN_ON_USE (logb, "logb is unportable - "
                 "use gnulib module logb for portability");
# endif
#endif


#if 0
# if !1 || !1
_GL_FUNCDECL_SYS (acosl, long double, (long double x));
# endif
_GL_CXXALIAS_SYS (acosl, long double, (long double x));
_GL_CXXALIASWARN (acosl);
#elif defined GNULIB_POSIXCHECK
# undef acosl
# if HAVE_RAW_DECL_ACOSL
_GL_WARN_ON_USE (acosl, "acosl is unportable - "
                 "use gnulib module mathl for portability");
# endif
#endif


#if 0
# if !1 || !1
_GL_FUNCDECL_SYS (asinl, long double, (long double x));
# endif
_GL_CXXALIAS_SYS (asinl, long double, (long double x));
_GL_CXXALIASWARN (asinl);
#elif defined GNULIB_POSIXCHECK
# undef asinl
# if HAVE_RAW_DECL_ASINL
_GL_WARN_ON_USE (asinl, "asinl is unportable - "
                 "use gnulib module mathl for portability");
# endif
#endif


#if 0
# if !1 || !1
_GL_FUNCDECL_SYS (atanl, long double, (long double x));
# endif
_GL_CXXALIAS_SYS (atanl, long double, (long double x));
_GL_CXXALIASWARN (atanl);
#elif defined GNULIB_POSIXCHECK
# undef atanl
# if HAVE_RAW_DECL_ATANL
_GL_WARN_ON_USE (atanl, "atanl is unportable - "
                 "use gnulib module mathl for portability");
# endif
#endif


#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   define ceilf rpl_ceilf
#  endif
_GL_FUNCDECL_RPL (ceilf, float, (float x));
_GL_CXXALIAS_RPL (ceilf, float, (float x));
# else
#  if !1
_GL_FUNCDECL_SYS (ceilf, float, (float x));
#  endif
_GL_CXXALIAS_SYS (ceilf, float, (float x));
# endif
_GL_CXXALIASWARN (ceilf);
#elif defined GNULIB_POSIXCHECK
# undef ceilf
# if HAVE_RAW_DECL_CEILF
_GL_WARN_ON_USE (ceilf, "ceilf is unportable - "
                 "use gnulib module ceilf for portability");
# endif
#endif

#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   define ceil rpl_ceil
#  endif
_GL_FUNCDECL_RPL (ceil, double, (double x));
_GL_CXXALIAS_RPL (ceil, double, (double x));
# else
_GL_CXXALIAS_SYS (ceil, double, (double x));
# endif
_GL_CXXALIASWARN (ceil);
#endif

#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   define ceill rpl_ceill
#  endif
_GL_FUNCDECL_RPL (ceill, long double, (long double x));
_GL_CXXALIAS_RPL (ceill, long double, (long double x));
# else
#  if !1
_GL_FUNCDECL_SYS (ceill, long double, (long double x));
#  endif
_GL_CXXALIAS_SYS (ceill, long double, (long double x));
# endif
_GL_CXXALIASWARN (ceill);
#elif defined GNULIB_POSIXCHECK
# undef ceill
# if HAVE_RAW_DECL_CEILL
_GL_WARN_ON_USE (ceill, "ceill is unportable - "
                 "use gnulib module ceill for portability");
# endif
#endif


#if 0
# if !1 || !1
_GL_FUNCDECL_SYS (cosl, long double, (long double x));
# endif
_GL_CXXALIAS_SYS (cosl, long double, (long double x));
_GL_CXXALIASWARN (cosl);
#elif defined GNULIB_POSIXCHECK
# undef cosl
# if HAVE_RAW_DECL_COSL
_GL_WARN_ON_USE (cosl, "cosl is unportable - "
                 "use gnulib module mathl for portability");
# endif
#endif


#if 0
# if !1 || !1
_GL_FUNCDECL_SYS (expl, long double, (long double x));
# endif
_GL_CXXALIAS_SYS (expl, long double, (long double x));
_GL_CXXALIASWARN (expl);
#elif defined GNULIB_POSIXCHECK
# undef expl
# if HAVE_RAW_DECL_EXPL
_GL_WARN_ON_USE (expl, "expl is unportable - "
                 "use gnulib module mathl for portability");
# endif
#endif


#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   define floorf rpl_floorf
#  endif
_GL_FUNCDECL_RPL (floorf, float, (float x));
_GL_CXXALIAS_RPL (floorf, float, (float x));
# else
#  if !1
_GL_FUNCDECL_SYS (floorf, float, (float x));
#  endif
_GL_CXXALIAS_SYS (floorf, float, (float x));
# endif
_GL_CXXALIASWARN (floorf);
#elif defined GNULIB_POSIXCHECK
# undef floorf
# if HAVE_RAW_DECL_FLOORF
_GL_WARN_ON_USE (floorf, "floorf is unportable - "
                 "use gnulib module floorf for portability");
# endif
#endif

#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   define floor rpl_floor
#  endif
_GL_FUNCDECL_RPL (floor, double, (double x));
_GL_CXXALIAS_RPL (floor, double, (double x));
# else
_GL_CXXALIAS_SYS (floor, double, (double x));
# endif
_GL_CXXALIASWARN (floor);
#endif

#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   define floorl rpl_floorl
#  endif
_GL_FUNCDECL_RPL (floorl, long double, (long double x));
_GL_CXXALIAS_RPL (floorl, long double, (long double x));
# else
#  if !1
_GL_FUNCDECL_SYS (floorl, long double, (long double x));
#  endif
_GL_CXXALIAS_SYS (floorl, long double, (long double x));
# endif
_GL_CXXALIASWARN (floorl);
#elif defined GNULIB_POSIXCHECK
# undef floorl
# if HAVE_RAW_DECL_FLOORL
_GL_WARN_ON_USE (floorl, "floorl is unportable - "
                 "use gnulib module floorl for portability");
# endif
#endif


/* Write x as
     x = mantissa * 2^exp
   where
     If x finite and nonzero: 0.5 <= |mantissa| < 1.0.
     If x is zero: mantissa = x, exp = 0.
     If x is infinite or NaN: mantissa = x, exp unspecified.
   Store exp in *EXPPTR and return mantissa.  */
#if 1 && 0
# if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#  define frexpl rpl_frexpl
# endif
_GL_FUNCDECL_RPL (frexpl, long double,
                  (long double x, int *expptr) _GL_ARG_NONNULL ((2)));
_GL_CXXALIAS_RPL (frexpl, long double, (long double x, int *expptr));
#else
# if !1
_GL_FUNCDECL_SYS (frexpl, long double,
                  (long double x, int *expptr) _GL_ARG_NONNULL ((2)));
# endif
# if 1
_GL_CXXALIAS_SYS (frexpl, long double, (long double x, int *expptr));
# endif
#endif
#if 1 && !(0 && !1)
_GL_CXXALIASWARN (frexpl);
#endif
#if !1 && defined GNULIB_POSIXCHECK
# undef frexpl
# if HAVE_RAW_DECL_FREXPL
_GL_WARN_ON_USE (frexpl, "frexpl is unportable - "
                 "use gnulib module frexpl for portability");
# endif
#endif


/* Return x * 2^exp.  */
#if 0 && 0
# if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#  define ldexpl rpl_ldexpl
# endif
_GL_FUNCDECL_RPL (ldexpl, long double, (long double x, int exp));
_GL_CXXALIAS_RPL (ldexpl, long double, (long double x, int exp));
#else
# if !1
_GL_FUNCDECL_SYS (ldexpl, long double, (long double x, int exp));
# endif
# if 0
_GL_CXXALIAS_SYS (ldexpl, long double, (long double x, int exp));
# endif
#endif
#if 0
_GL_CXXALIASWARN (ldexpl);
#endif
#if !0 && defined GNULIB_POSIXCHECK
# undef ldexpl
# if HAVE_RAW_DECL_LDEXPL
_GL_WARN_ON_USE (ldexpl, "ldexpl is unportable - "
                 "use gnulib module ldexpl for portability");
# endif
#endif


#if 0
# if !1 || !1
_GL_FUNCDECL_SYS (logl, long double, (long double x));
# endif
_GL_CXXALIAS_SYS (logl, long double, (long double x));
_GL_CXXALIASWARN (logl);
#elif defined GNULIB_POSIXCHECK
# undef logl
# if HAVE_RAW_DECL_LOGL
_GL_WARN_ON_USE (logl, "logl is unportable - "
                 "use gnulib module mathl for portability");
# endif
#endif


#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   undef roundf
#   define roundf rpl_roundf
#  endif
_GL_FUNCDECL_RPL (roundf, float, (float x));
_GL_CXXALIAS_RPL (roundf, float, (float x));
# else
#  if !1
_GL_FUNCDECL_SYS (roundf, float, (float x));
#  endif
_GL_CXXALIAS_SYS (roundf, float, (float x));
# endif
_GL_CXXALIASWARN (roundf);
#elif defined GNULIB_POSIXCHECK
# undef roundf
# if HAVE_RAW_DECL_ROUNDF
_GL_WARN_ON_USE (roundf, "roundf is unportable - "
                 "use gnulib module roundf for portability");
# endif
#endif

#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   undef round
#   define round rpl_round
#  endif
_GL_FUNCDECL_RPL (round, double, (double x));
_GL_CXXALIAS_RPL (round, double, (double x));
# else
#  if !1
_GL_FUNCDECL_SYS (round, double, (double x));
#  endif
_GL_CXXALIAS_SYS (round, double, (double x));
# endif
_GL_CXXALIASWARN (round);
#elif defined GNULIB_POSIXCHECK
# undef round
# if HAVE_RAW_DECL_ROUND
_GL_WARN_ON_USE (round, "round is unportable - "
                 "use gnulib module round for portability");
# endif
#endif

#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   undef roundl
#   define roundl rpl_roundl
#  endif
_GL_FUNCDECL_RPL (roundl, long double, (long double x));
_GL_CXXALIAS_RPL (roundl, long double, (long double x));
# else
#  if !1
_GL_FUNCDECL_SYS (roundl, long double, (long double x));
#  endif
_GL_CXXALIAS_SYS (roundl, long double, (long double x));
# endif
_GL_CXXALIASWARN (roundl);
#elif defined GNULIB_POSIXCHECK
# undef roundl
# if HAVE_RAW_DECL_ROUNDL
_GL_WARN_ON_USE (roundl, "roundl is unportable - "
                 "use gnulib module roundl for portability");
# endif
#endif


#if 0
# if !1 || !1
_GL_FUNCDECL_SYS (sinl, long double, (long double x));
# endif
_GL_CXXALIAS_SYS (sinl, long double, (long double x));
_GL_CXXALIASWARN (sinl);
#elif defined GNULIB_POSIXCHECK
# undef sinl
# if HAVE_RAW_DECL_SINL
_GL_WARN_ON_USE (sinl, "sinl is unportable - "
                 "use gnulib module mathl for portability");
# endif
#endif


#if 0
# if !1 || !1
_GL_FUNCDECL_SYS (sqrtl, long double, (long double x));
# endif
_GL_CXXALIAS_SYS (sqrtl, long double, (long double x));
_GL_CXXALIASWARN (sqrtl);
#elif defined GNULIB_POSIXCHECK
# undef sqrtl
# if HAVE_RAW_DECL_SQRTL
_GL_WARN_ON_USE (sqrtl, "sqrtl is unportable - "
                 "use gnulib module mathl for portability");
# endif
#endif


#if 0
# if !1 || !1
_GL_FUNCDECL_SYS (tanl, long double, (long double x));
# endif
_GL_CXXALIAS_SYS (tanl, long double, (long double x));
_GL_CXXALIASWARN (tanl);
#elif defined GNULIB_POSIXCHECK
# undef tanl
# if HAVE_RAW_DECL_TANL
_GL_WARN_ON_USE (tanl, "tanl is unportable - "
                 "use gnulib module mathl for portability");
# endif
#endif


#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   define truncf rpl_truncf
#  endif
_GL_FUNCDECL_RPL (truncf, float, (float x));
_GL_CXXALIAS_RPL (truncf, float, (float x));
# else
#  if !1
_GL_FUNCDECL_SYS (truncf, float, (float x));
#  endif
_GL_CXXALIAS_SYS (truncf, float, (float x));
# endif
_GL_CXXALIASWARN (truncf);
#elif defined GNULIB_POSIXCHECK
# undef truncf
# if HAVE_RAW_DECL_TRUNCF
_GL_WARN_ON_USE (truncf, "truncf is unportable - "
                 "use gnulib module truncf for portability");
# endif
#endif

#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   define trunc rpl_trunc
#  endif
_GL_FUNCDECL_RPL (trunc, double, (double x));
_GL_CXXALIAS_RPL (trunc, double, (double x));
# else
#  if !1
_GL_FUNCDECL_SYS (trunc, double, (double x));
#  endif
_GL_CXXALIAS_SYS (trunc, double, (double x));
# endif
_GL_CXXALIASWARN (trunc);
#elif defined GNULIB_POSIXCHECK
# undef trunc
# if HAVE_RAW_DECL_TRUNC
_GL_WARN_ON_USE (trunc, "trunc is unportable - "
                 "use gnulib module trunc for portability");
# endif
#endif

#if 0
# if 0
#  if !(defined __cplusplus && defined GNULIB_NAMESPACE)
#   undef truncl
#   define truncl rpl_truncl
#  endif
_GL_FUNCDECL_RPL (truncl, long double, (long double x));
_GL_CXXALIAS_RPL (truncl, long double, (long double x));
# else
#  if !1
_GL_FUNCDECL_SYS (truncl, long double, (long double x));
#  endif
_GL_CXXALIAS_SYS (truncl, long double, (long double x));
# endif
_GL_CXXALIASWARN (truncl);
#elif defined GNULIB_POSIXCHECK
# undef truncl
# if HAVE_RAW_DECL_TRUNCL
_GL_WARN_ON_USE (truncl, "truncl is unportable - "
                 "use gnulib module truncl for portability");
# endif
#endif


#if 0
# if 0
_GL_EXTERN_C int gl_isfinitef (float x);
_GL_EXTERN_C int gl_isfinited (double x);
_GL_EXTERN_C int gl_isfinitel (long double x);
#  undef isfinite
#  define isfinite(x) \
   (sizeof (x) == sizeof (long double) ? gl_isfinitel (x) : \
    sizeof (x) == sizeof (double) ? gl_isfinited (x) : \
    gl_isfinitef (x))
# endif
#elif defined GNULIB_POSIXCHECK
# if defined isfinite
_GL_WARN_REAL_FLOATING_DECL (isfinite);
#  undef isfinite
#  define isfinite(x) _GL_WARN_REAL_FLOATING_IMPL (isfinite, x)
# endif
#endif


#if 0
# if 0
_GL_EXTERN_C int gl_isinff (float x);
_GL_EXTERN_C int gl_isinfd (double x);
_GL_EXTERN_C int gl_isinfl (long double x);
#  undef isinf
#  define isinf(x) \
   (sizeof (x) == sizeof (long double) ? gl_isinfl (x) : \
    sizeof (x) == sizeof (double) ? gl_isinfd (x) : \
    gl_isinff (x))
# endif
#elif defined GNULIB_POSIXCHECK
# if defined isinf
_GL_WARN_REAL_FLOATING_DECL (isinf);
#  undef isinf
#  define isinf(x) _GL_WARN_REAL_FLOATING_IMPL (isinf, x)
# endif
#endif


#if 0
/* Test for NaN for 'float' numbers.  */
# if 1
/* The original <math.h> included above provides a declaration of isnan macro
   or (older) isnanf function.  */
#  if __GNUC__ >= 4
    /* GCC 4.0 and newer provides three built-ins for isnan.  */
#   undef isnanf
#   define isnanf(x) __builtin_isnanf ((float)(x))
#  elif defined isnan
#   undef isnanf
#   define isnanf(x) isnan ((float)(x))
#  endif
# else
/* Test whether X is a NaN.  */
#  undef isnanf
#  define isnanf rpl_isnanf
_GL_EXTERN_C int isnanf (float x);
# endif
#endif

#if 0
/* Test for NaN for 'double' numbers.
   This function is a gnulib extension, unlike isnan() which applied only
   to 'double' numbers earlier but now is a type-generic macro.  */
# if 1
/* The original <math.h> included above provides a declaration of isnan
   macro.  */
#  if __GNUC__ >= 4
    /* GCC 4.0 and newer provides three built-ins for isnan.  */
#   undef isnand
#   define isnand(x) __builtin_isnan ((double)(x))
#  else
#   undef isnand
#   define isnand(x) isnan ((double)(x))
#  endif
# else
/* Test whether X is a NaN.  */
#  undef isnand
#  define isnand rpl_isnand
_GL_EXTERN_C int isnand (double x);
# endif
#endif

#if 0
/* Test for NaN for 'long double' numbers.  */
# if 1
/* The original <math.h> included above provides a declaration of isnan
   macro or (older) isnanl function.  */
#  if __GNUC__ >= 4
    /* GCC 4.0 and newer provides three built-ins for isnan.  */
#   undef isnanl
#   define isnanl(x) __builtin_isnanl ((long double)(x))
#  elif defined isnan
#   undef isnanl
#   define isnanl(x) isnan ((long double)(x))
#  endif
# else
/* Test whether X is a NaN.  */
#  undef isnanl
#  define isnanl rpl_isnanl
_GL_EXTERN_C int isnanl (long double x);
# endif
#endif

/* This must come *after* the snippets for GNULIB_ISNANF and GNULIB_ISNANL!  */
#if 0
# if 0
/* We can't just use the isnanf macro (e.g.) as exposed by
   isnanf.h (e.g.) here, because those may end up being macros
   that recursively expand back to isnan.  So use the gnulib
   replacements for them directly. */
#  if 1 && __GNUC__ >= 4
#   define gl_isnan_f(x) __builtin_isnanf ((float)(x))
#  else
_GL_EXTERN_C int rpl_isnanf (float x);
#   define gl_isnan_f(x) rpl_isnanf (x)
#  endif
#  if 1 && __GNUC__ >= 4
#   define gl_isnan_d(x) __builtin_isnan ((double)(x))
#  else
_GL_EXTERN_C int rpl_isnand (double x);
#   define gl_isnan_d(x) rpl_isnand (x)
#  endif
#  if 1 && __GNUC__ >= 4
#   define gl_isnan_l(x) __builtin_isnanl ((long double)(x))
#  else
_GL_EXTERN_C int rpl_isnanl (long double x);
#   define gl_isnan_l(x) rpl_isnanl (x)
#  endif
#  undef isnan
#  define isnan(x) \
   (sizeof (x) == sizeof (long double) ? gl_isnan_l (x) : \
    sizeof (x) == sizeof (double) ? gl_isnan_d (x) : \
    gl_isnan_f (x))
# elif __GNUC__ >= 4
#  undef isnan
#  define isnan(x) \
   (sizeof (x) == sizeof (long double) ? __builtin_isnanl ((long double)(x)) : \
    sizeof (x) == sizeof (double) ? __builtin_isnan ((double)(x)) : \
    __builtin_isnanf ((float)(x)))
# endif
/* Ensure isnan is a macro.  */
# ifndef isnan
#  define isnan isnan
# endif
#elif defined GNULIB_POSIXCHECK
# if defined isnan
_GL_WARN_REAL_FLOATING_DECL (isnan);
#  undef isnan
#  define isnan(x) _GL_WARN_REAL_FLOATING_IMPL (isnan, x)
# endif
#endif


#if 1
# if 1
#  undef signbit
   /* GCC 4.0 and newer provides three built-ins for signbit.  */
#  define signbit(x) \
   (sizeof (x) == sizeof (long double) ? __builtin_signbitl (x) : \
    sizeof (x) == sizeof (double) ? __builtin_signbit (x) : \
    __builtin_signbitf (x))
# endif
# if 0
#  undef signbit
_GL_EXTERN_C int gl_signbitf (float arg);
_GL_EXTERN_C int gl_signbitd (double arg);
_GL_EXTERN_C int gl_signbitl (long double arg);
#  if __GNUC__ >= 2 && !__STRICT_ANSI__
#   define _GL_NUM_UINT_WORDS(type) \
      ((sizeof (type) + sizeof (unsigned int) - 1) / sizeof (unsigned int))
#   if defined FLT_SIGNBIT_WORD && defined FLT_SIGNBIT_BIT && !defined gl_signbitf
#    define gl_signbitf_OPTIMIZED_MACRO
#    define gl_signbitf(arg) \
       ({ union { float _value;                                         \
                  unsigned int _word[_GL_NUM_UINT_WORDS (float)];       \
                } _m;                                                   \
          _m._value = (arg);                                            \
          (_m._word[FLT_SIGNBIT_WORD] >> FLT_SIGNBIT_BIT) & 1;          \
        })
#   endif
#   if defined DBL_SIGNBIT_WORD && defined DBL_SIGNBIT_BIT && !defined gl_signbitd
#    define gl_signbitd_OPTIMIZED_MACRO
#    define gl_signbitd(arg) \
       ({ union { double _value;                                        \
                  unsigned int _word[_GL_NUM_UINT_WORDS (double)];      \
                } _m;                                                   \
          _m._value = (arg);                                            \
          (_m._word[DBL_SIGNBIT_WORD] >> DBL_SIGNBIT_BIT) & 1;          \
        })
#   endif
#   if defined LDBL_SIGNBIT_WORD && defined LDBL_SIGNBIT_BIT && !defined gl_signbitl
#    define gl_signbitl_OPTIMIZED_MACRO
#    define gl_signbitl(arg) \
       ({ union { long double _value;                                   \
                  unsigned int _word[_GL_NUM_UINT_WORDS (long double)]; \
                } _m;                                                   \
          _m._value = (arg);                                            \
          (_m._word[LDBL_SIGNBIT_WORD] >> LDBL_SIGNBIT_BIT) & 1;        \
        })
#   endif
#  endif
#  define signbit(x) \
   (sizeof (x) == sizeof (long double) ? gl_signbitl (x) : \
    sizeof (x) == sizeof (double) ? gl_signbitd (x) : \
    gl_signbitf (x))
# endif
#elif defined GNULIB_POSIXCHECK
# if defined signbit
_GL_WARN_REAL_FLOATING_DECL (signbit);
#  undef signbit
#  define signbit(x) _GL_WARN_REAL_FLOATING_IMPL (signbit, x)
# endif
#endif


#endif /* _GL_MATH_H */
#endif /* _GL_MATH_H */
