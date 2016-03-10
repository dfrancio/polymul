#ifndef POLYMUL_H
#define POLYMUL_H

#include <cassert>
#include <iostream>

/*
  Multiply polynomials p1 and p2, both polynomials of nvar variables,
  of degree n1 and n2 respectively.  The result is a polynomial of
  degree n1+n2, which is _added_ to dst.

  A degree N polynomial of K variables has
  /N + K\
  \  N  /   terms.

  A polynomial P[K,N] is here organized as a sequence of monomials in
  K variables:

  M[K,0], M[K,1] .. M[K,N]  (for example 1, x y, x^2 xy y^2, ... )

  However, for each monomial all terms are of the same degree and the
  exponents of one variable can be inferred from the other
  exponents. We therefore have

  M[K,N] ~ P[K-1,N] (example: x^2 xy y^2 ~ 1 y y^2).

  To multiply monomials we can then use the ordinary polymul()
  function with lower nvar value.
*/

/*
 * Some compilers implement the C99 restrict also in C++. Need to find a way to differentiate
 * between newer g++ and other compilers that define __GNUC__ (like Intel icc).
 */
#ifdef __GNUC_AND_NOT_ICC__
#define POLYMUL_RESTRICT __restrict__
#else
#define POLYMUL_RESTRICT
#endif

namespace polymul
{
namespace detail
{
/*
 * compile-time factorial computation
 */
template <int n>
struct factorial
{
    enum
    {
        fac = n * factorial<n - 1>::fac
    };
};

template <>
struct factorial<0>
{
    enum
    {
        fac = 1
    };
};

template <>
struct factorial<-1>
{
    enum
    {
        fac = 1
    };
};

/*
 * Template to evaluate binomial coefficients at compile time
 */
template <int n, int k>
struct binomial
{
    enum
    {
        value = binomial<n - 1, k - 1>::value + binomial<n - 1, k>::value
    };
};

template <int n>
struct binomial<n, n>
{
    enum
    {
        value = 1
    };
};

template <>
struct binomial<0, 0>
{
    enum
    {
        value = 1
    };
};

template <int n>
struct binomial<n, 0>
{
    enum
    {
        value = n > 0 ? 1 : 0
    };
};

template <int k>
struct binomial<0, k>
{
    enum
    {
        value = k > 0 ? 1 : 0
    };
};

/*
 * number of coefficients of a polynomial
 */
template <int Nvar, int Ndeg>
struct polylen
{
    enum
    {
        len = binomial<Nvar + Ndeg, Ndeg>::value
    };
};

/*
 * term <-> degree
 */
template <int Nvar, int term>
struct term_deg
{
    enum
    {
        deg = term_deg<Nvar, term - 1>::deg
            + ((term == polylen<Nvar, term_deg<Nvar, term - 1>::deg>::len) ? 1 : 0)
    };
};

template <int Nvar>
struct term_deg<Nvar, 0>
{
    enum
    {
        deg = 0
    };
};

template <int term>
struct term_deg<1, term>
{
    enum
    {
        deg = term
    };
};

template <>
struct term_deg<1, 0>
{
    enum
    {
        deg = 0
    };
};

/*
 * term -> term within highest order
 */
template <int Nvar, int term>
struct term_sub
{
    enum
    {
        sub = term - polylen<Nvar, term_deg<Nvar, term>::deg - 1>::len
    };
};

/*
 * term multiplication table
 */
template <int Nvar, int t1, int t2>
struct term_prod
{
    enum
    {
        prod = polylen<Nvar, term_deg<Nvar, t1>::deg + term_deg<Nvar, t2>::deg - 1>::len
            + term_prod<Nvar - 1, t1 - polylen<Nvar, term_deg<Nvar, t1>::deg - 1>::len, t2 - polylen<Nvar, term_deg<Nvar, t2>::deg - 1>::len>::
                  prod
    };
};

template <int Nvar, int t2>
struct term_prod<Nvar, 0, t2>
{
    enum
    {
        prod = t2
    };
};

template <int Nvar, int t1>
struct term_prod<Nvar, t1, 0>
{
    enum
    {
        prod = t1
    };
};

template <int Nvar>
struct term_prod<Nvar, 0, 0>
{
    enum
    {
        prod = 0
    };
};

template <>
struct term_prod<0, 0, 0>
{
    enum
    {
        prod = 0
    };
};

/*
 * Determining the term exponents at compile time
 */
template <int Nvar, int term>
struct first_exponent
{
    enum
    {
        e = term_deg<Nvar, term>::deg - term_deg<Nvar - 1, term_sub<Nvar, term>::sub>::deg
    };
};

template <int term>
struct first_exponent<1, term>
{
    enum
    {
        e = term
    };
};

template <int Nvar>
struct first_exponent<Nvar, 0>
{
    enum
    {
        e = 0
    };
};

template <>
struct first_exponent<1, 0>
{
    enum
    {
        e = 0
    };
};

/*
 * Derivative factors n!m!k!.. (for term x^n y^m z^k..)
// WARNING: may overflow (do we get compiler warning from this?)
 */
template <int Nvar, int term>
struct deriv_fac
{
    enum
    {
        fac = factorial<first_exponent<Nvar, term>::e>::fac
            * deriv_fac<Nvar - 1, term_sub<Nvar, term>::sub>::fac
    };
};

template <int term>
struct deriv_fac<1, term>
{
    enum
    {
        fac = factorial<term>::fac
    };
};

template <int Nvar>
struct deriv_fac<Nvar, 0>
{
    enum
    {
        fac = 1
    };
};

template <>
struct deriv_fac<1, 0>
{
    enum
    {
        fac = 1
    };
};

/*
 * deriv_fac_multiplier
 */
template <class numtype, int Nvar, int term>
struct deriv_fac_multiplier
{
    static void mul_fac(numtype c[])
    {
        deriv_fac_multiplier<numtype, Nvar, term - 1>::mul_fac(c);
        c[term] *= deriv_fac<Nvar, term>::fac;
    }
};

template <class numtype, int Nvar>
struct deriv_fac_multiplier<numtype, Nvar, 0>
{
    static void mul_fac(numtype c[])
    {
    }
};

/*
 * Recursive template classes for multiplication
 */
template <class numtype, int Nvar, int Ndeg1, int Ndeg2>
struct polynomial_multiplier
{
    /*
     * _add_ the product between p1 and p2 to dst
     */
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        polynomial_multiplier<numtype, Nvar, Ndeg1, Ndeg2 - 1>::mul(dst, p1, p2);
        polynomial_multiplier<numtype, Nvar, Ndeg1, Ndeg2>::mul_monomial(
            dst, p1, p2 + binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value);
    }

    /*
     * m2 is a monomial in Nvar variables and of order Ndeg2.
     * _add_ the product of p1 and m2 to dst
     */
    static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        polynomial_multiplier<numtype, Nvar - 1, Ndeg1, Ndeg2>::mul(
            dst + binomial<Nvar + Ndeg1 + Ndeg2 - 1, Ndeg1 + Ndeg2 - 1>::value,
            p1 + binomial<Nvar + Ndeg1 - 1, Ndeg1 - 1>::value, m2);
        polynomial_multiplier<numtype, Nvar, Ndeg1 - 1, Ndeg2>::mul_monomial(dst, p1, m2);
    }

    /*
     * mul_set
     */
    static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        polynomial_multiplier<numtype, Nvar, Ndeg1, Ndeg2 - 1>::mul_set(dst, p1, p2);
        // now we just add, because lower part is already set.
        polynomial_multiplier<numtype, Nvar, Ndeg1 - 1, Ndeg2>::mul_monomial(
            dst, p1, p2 + binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value);
        // Set final two highest monomials
        polynomial_multiplier<numtype, Nvar - 1, Ndeg1, Ndeg2>::mul_set(
            dst + binomial<Nvar + Ndeg1 + Ndeg2 - 1, Ndeg1 + Ndeg2 - 1>::value,
            p1 + binomial<Nvar + Ndeg1 - 1, Ndeg1 - 1>::value,
            p2 + binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value);
    }

    /*
     * mul_monomial_set
     */
    static void mul_monomial_set(
        numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        // Unclear what this would mean, since the lower parts of dst are not touched.
        assert(0 && " I am not supposed to be here..");
    }

    /*
     * See polycontract below for a description of what this is for
     */
    static void antimul(const numtype dst[], numtype p1[], const numtype p2[])
    {
        polynomial_multiplier<numtype, Nvar, Ndeg1, Ndeg2 - 1>::antimul(dst, p1, p2);
        polynomial_multiplier<numtype, Nvar, Ndeg1, Ndeg2>::antimul_monomial(
            dst, p1, p2 + binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value);
    }

    /*
     * antimul_monomial
     */
    static void antimul_monomial(const numtype dst[], numtype p1[], const numtype m2[])
    {
        polynomial_multiplier<numtype, Nvar - 1, Ndeg1, Ndeg2>::antimul(
            dst + binomial<Nvar + Ndeg1 + Ndeg2 - 1, Ndeg1 + Ndeg2 - 1>::value,
            p1 + binomial<Nvar + Ndeg1 - 1, Ndeg1 - 1>::value, m2);
        polynomial_multiplier<numtype, Nvar, Ndeg1 - 1, Ndeg2>::antimul_monomial(dst, p1, m2);
    }
};

/*
 * polynomial_multiplier
 */
template <class numtype, int Ndeg1, int Ndeg2>
struct polynomial_multiplier<numtype, 1, Ndeg1, Ndeg2>
{
    /*
     * mul
     */
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        for (int i = 0; i <= Ndeg1; i++)
            for (int j = 0; j <= Ndeg2; j++)
                dst[i + j] += p1[i] * p2[j];
    }

    /*
     * mul_monomial
     */
    static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        for (int i = 0; i <= Ndeg1; i++)
            dst[i + Ndeg2] += m2[0] * p1[i];
    }

    /*
     * mul_set
     */
    static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        for (int i = 0; i <= Ndeg1; i++)
            dst[i] = p1[i] * p2[0];
        // Now set all where i = Ndeg1 and j > 0
        for (int j = 1; j <= Ndeg2; j++)
            dst[Ndeg1 + j] = p1[Ndeg1] * p2[j];
        // Add all inbetween
        for (int i = 0; i < Ndeg1; i++)
            for (int j = 1; j <= Ndeg2; j++)
                dst[i + j] += p1[i] * p2[j];
    }

    /*
     * mul_monomial_set
     */
    static void mul_monomial_set(
        numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        for (int i = 0; i < Ndeg1; i++)
            dst[i + Ndeg2] = m2[0] * p1[i];
    }

    /*
     * antimul
     * Implementation details go here...
     */
    static void antimul(const numtype dst[], numtype p1[], const numtype p2[])
    {
        for (int i = 0; i <= Ndeg1; i++)
            for (int j = 0; j <= Ndeg2; j++)
                p1[i] += p2[j] * dst[i + j];
    }

    /*
     * antimul_monomial
     * Implementation details go here...
     */
    static void antimul_monomial(const numtype dst[], numtype p1[], const numtype m2[])
    {
        for (int i = 0; i <= Ndeg1; i++)
            p1[i] += m2[0] * dst[i + Ndeg2];
    }
};

/*
 * 0 variables: this is needed as a limiting case to make other code simpler.
 * We have only the constant term.
 */
template <class numtype, int Ndeg1, int Ndeg2>
struct polynomial_multiplier<numtype, 0, Ndeg1, Ndeg2>
{
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        dst[0] += p1[0] * p2[0];
    }
    static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        assert(0 && "FIXME");
    }
    static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        dst[0] = p1[0] * p2[0];
    }
    static void mul_monomial_set(
        numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        assert(0 && "FIXME");
    }

    static void antimul(const numtype dst[], numtype p1[], const numtype p2[])
    {
        assert(0 && "FIXME");
    }
    static void antimul_monomial(const numtype dst[], numtype p1[], const numtype m2[])
    {
        assert(0 && "FIXME");
    }
};

template <class numtype, int Nvar, int Ndeg2>
struct polynomial_multiplier<numtype, Nvar, 0, Ndeg2>
{
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<Nvar + Ndeg2, Ndeg2>::value; i++)
            dst[i] += p1[0] * p2[i];
    }
    static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        for (int i = binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value;
             i < binomial<Nvar + Ndeg2, Ndeg2>::value; i++)
            dst[i] += p1[0] * m2[i - binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value];
    }
    static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<Nvar + Ndeg2, Ndeg2>::value; i++)
            dst[i] = p1[0] * p2[i];
    }
    static void mul_monomial_set(
        numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        for (int i = binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value;
             i < binomial<Nvar + Ndeg2, Ndeg2>::value; i++)
            dst[i] = p1[0] * m2[i - binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value];
    }
    static void antimul(const numtype dst[], numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<Nvar + Ndeg2, Ndeg2>::value; i++)
            p1[0] += dst[i] * p2[i];
    }
    static void antimul_monomial(const numtype dst[], numtype p1[], const numtype m2[])
    {
        for (int i = binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value;
             i < binomial<Nvar + Ndeg2, Ndeg2>::value; i++)
            p1[0] += dst[i] * m2[i - binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value];
    }
};

template <class numtype, int Ndeg2>
struct polynomial_multiplier<numtype, 1, 0, Ndeg2>
{
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<1 + Ndeg2, Ndeg2>::value; i++)
            dst[i] += p1[0] * p2[i];
    }
    static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        for (int i = binomial<1 + Ndeg2 - 1, Ndeg2 - 1>::value;
             i < binomial<1 + Ndeg2, Ndeg2>::value; i++)
            dst[i] += p1[0] * m2[i - binomial<1 + Ndeg2 - 1, Ndeg2 - 1>::value];
    }
    static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<1 + Ndeg2, Ndeg2>::value; i++)
            dst[i] = p1[0] * p2[i];
    }
    static void mul_monomial_set(
        numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        for (int i = binomial<1 + Ndeg2 - 1, Ndeg2 - 1>::value;
             i < binomial<1 + Ndeg2, Ndeg2>::value; i++)
            dst[i] = p1[0] * m2[i - binomial<1 + Ndeg2 - 1, Ndeg2 - 1>::value];
    }
    static void antimul(const numtype dst[], numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<1 + Ndeg2, Ndeg2>::value; i++)
            p1[0] += dst[i] * p2[i];
    }
    static void antimul_monomial(const numtype dst[], numtype p1[], const numtype m2[])
    {
        for (int i = binomial<1 + Ndeg2 - 1, Ndeg2 - 1>::value;
             i < binomial<1 + Ndeg2, Ndeg2>::value; i++)
            p1[0] += dst[i] * m2[i - binomial<1 + Ndeg2 - 1, Ndeg2 - 1>::value];
    }
};

template <class numtype, int Nvar, int Ndeg1>
struct polynomial_multiplier<numtype, Nvar, Ndeg1, 0>
{
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<Nvar + Ndeg1, Ndeg1>::value; i++)
            dst[i] += p1[i] * p2[0];
    }
    static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        for (int i = 0; i < binomial<Nvar + Ndeg1, Ndeg1>::value; i++)
            dst[i] += p1[i] * m2[0];
    }
    static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<Nvar + Ndeg1, Ndeg1>::value; i++)
            dst[i] = p1[i] * p2[0];
    }
    static void mul_monomial_set(
        numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        for (int i = 0; i < binomial<Nvar + Ndeg1, Ndeg1>::value; i++)
            dst[i] = p1[i] * m2[0];
    }
    static void antimul(const numtype dst[], numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<Nvar + Ndeg1, Ndeg1>::value; i++)
            p1[i] += dst[i] * p2[0];
    }
    static void antimul_monomial(const numtype dst[], numtype p1[], const numtype m2[])
    {
        for (int i = 0; i < binomial<Nvar + Ndeg1, Ndeg1>::value; i++)
            p1[i] += dst[i] * m2[0];
    }
};

template <class numtype, int Nvar>
struct polynomial_multiplier<numtype, Nvar, 0, 0>
{
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        dst[0] += p1[0] * p2[0];
    }
    static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        dst[0] += p1[0] * m2[0];
    }
    static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        dst[0] = p1[0] * p2[0];
    }
    static void mul_monomial_set(
        numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        dst[0] = p1[0] * m2[0];
    }
    static void antimul(const numtype dst[], numtype p1[], const numtype p2[])
    {
        p1[0] += dst[0] * p2[0];
    }
    static void antimul_monomial(const numtype dst[], numtype p1[], const numtype m2[])
    {
        p1[0] += dst[0] * m2[0];
    }
};

template <class numtype>
struct polynomial_multiplier<numtype, 1, 0, 0>
{
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        dst[0] += p1[0] * p2[0];
    }
    static void mul_monomial(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        dst[0] += p1[0] * m2[0];
    }
    static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        dst[0] = p1[0] * p2[0];
    }
    static void mul_monomial_set(
        numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype m2[])
    {
        dst[0] = p1[0] * m2[0];
    }
    static void antimul(const numtype dst[], numtype p1[], const numtype p2[])
    {
        p1[0] += dst[0] * p2[0];
    }
    static void antimul_monomial(const numtype dst[], numtype p1[], const numtype m2[])
    {
        p1[0] += dst[0] * m2[0];
    }
};

/*
 * Like polymul but truncates the result to Ndeg1 order.
 * Ndeg2 _must_ be less or equal to Ndeg1, otherwise the recursion will not terminate.
 */
template <class numtype, int Nvar, int Ndeg1, int Ndeg2>
struct taylor_multiplier
{
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        polynomial_multiplier<numtype, Nvar, Ndeg1 - Ndeg2, Ndeg2>::mul_monomial(
            dst, p1, p2 + binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value);
        taylor_multiplier<numtype, Nvar, Ndeg1, Ndeg2 - 1>::mul(dst, p1, p2);
    }
    static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        taylor_multiplier<numtype, Nvar, Ndeg1, Ndeg2 - 1>::mul_set(dst, p1, p2);
        polynomial_multiplier<numtype, Nvar, Ndeg1 - Ndeg2, Ndeg2>::mul_monomial(
            dst, p1, p2 + binomial<Nvar + Ndeg2 - 1, Ndeg2 - 1>::value);
    }
};

/*
 * (1)
 */
template <class numtype, int Nvar, int Ndeg1>
struct taylor_multiplier<numtype, Nvar, Ndeg1, 0>
{
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<Ndeg1 + Nvar, Ndeg1>::value; i++)
            dst[i] += p1[i] * p2[0];
    }
    static void mul_set(numtype POLYMUL_RESTRICT dst[], const numtype p1[], const numtype p2[])
    {
        for (int i = 0; i < binomial<Ndeg1 + Nvar, Ndeg1>::value; i++)
            dst[i] = p1[i] * p2[0];
    }
};

/*
 * (2)
 */
template <class numtype, int Nvar, int Ndeg, int Ndeg2, int i2>
struct taylor_inplace_multiplier
{
    static void mul(numtype p1[], const numtype p2[])
    {
        if (i2 <= Ndeg2)
        {
            // M1(Ndeg-i2)*M2(i2) -> M1(Ndeg)
            polynomial_multiplier<numtype, Nvar - 1, Ndeg - i2, i2>::mul(
                p1 + binomial<Nvar + Ndeg - 1, Ndeg - 1>::value,
                p1 + binomial<Nvar + Ndeg - i2 - 1, Ndeg - i2 - 1>::value,
                p2 + binomial<Nvar + i2 - 1, i2 - 1>::value);
        }
        // Back to (2), or to (3) when i2+1 == Ndeg
        taylor_inplace_multiplier<numtype, Nvar, Ndeg, Ndeg2, i2 + 1>::mul(p1, p2);
    }
};

/*
 * (3) final contribution to M1(Ndeg)
 */
template <class numtype, int Nvar, int Ndeg, int Ndeg2>
struct taylor_inplace_multiplier<numtype, Nvar, Ndeg, Ndeg2, Ndeg>
{
    static void mul(numtype POLYMUL_RESTRICT p1[], const numtype p2[])
    {
        if (Ndeg <= Ndeg2)
        {
            for (int i = binomial<Nvar + Ndeg - 1, Ndeg - 1>::value;
                 i < binomial<Nvar + Ndeg, Ndeg>::value; i++)
                p1[i] += p2[i] * p1[0];
        }
        taylor_inplace_multiplier<numtype, Nvar, Ndeg - 1, Ndeg2, 0>::mul(
            p1, p2); // Do lower degree terms, or go to (4)
    }
};

/*
 * (1) comes in here, sets M1(Ndeg)
 */
template <class numtype, int Nvar, int Ndeg, int Ndeg2>
struct taylor_inplace_multiplier<numtype, Nvar, Ndeg, Ndeg2, 0>
{
    static void mul(numtype POLYMUL_RESTRICT p1[], const numtype p2[])
    {
        for (int i = binomial<Nvar + Ndeg - 1, Ndeg - 1>::value;
             i < binomial<Nvar + Ndeg, Ndeg>::value; i++)
            p1[i] *= p2[0];
        taylor_inplace_multiplier<numtype, Nvar, Ndeg, Ndeg2, 1>::mul(p1, p2); // continue at (2)
    }
};

/*
 * (4), last coefficient.
 */
template <class numtype, int Nvar, int Ndeg2>
struct taylor_inplace_multiplier<numtype, Nvar, 0, Ndeg2, 0>
{
    static void mul(numtype POLYMUL_RESTRICT p1[], const numtype p2[])
    {
        p1[0] *= p2[0];
    }
};

/*
 * (4), last coefficient.
 */
template <class numtype, int Ndeg2>
struct taylor_inplace_multiplier<numtype, 1, 0, Ndeg2, 0>
{
    static void mul(numtype POLYMUL_RESTRICT p1[], const numtype p2[])
    {
        p1[0] *= p2[0];
    }
};

/*
 * Above code is only for Ndeg>1
 */
template <class numtype, int Ndeg, int Ndeg2>
struct taylor_inplace_multiplier<numtype, 1, Ndeg, Ndeg2, 0>
{
    static void mul(numtype POLYMUL_RESTRICT p1[], const numtype p2[])
    {
        p1[Ndeg] *= p2[0];
        if (Ndeg <= Ndeg2)
        {
            for (int i = 0; i < Ndeg; i++)
                p1[Ndeg] += p1[i] * p2[Ndeg - i];
        }
        else
        {
            for (int i = Ndeg - Ndeg2; i < Ndeg; i++)
                p1[Ndeg] += p1[i] * p2[Ndeg - i];
        }
        taylor_inplace_multiplier<numtype, 1, Ndeg - 1, Ndeg2, 0>::mul(p1, p2);
    }
};

/*
 * term_multiplier
 */
template <class numtype, int Nvar, int lterm, int rterm>
struct term_multiplier
{
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p[], const numtype &c)
    {
        term_multiplier<numtype, Nvar, lterm - 1, rterm>::mul(dst, p, c);
        dst[term_prod<Nvar, lterm, rterm>::prod] += p[lterm] * c;
    }
};

template <class numtype, int Nvar, int rterm>
struct term_multiplier<numtype, Nvar, 0, rterm>
{
    static void mul(numtype POLYMUL_RESTRICT dst[], const numtype p[], const numtype &c)
    {
        dst[rterm] += p[0] * c;
    }
};

/*
 * polynomial_evaluator
 */
template <class numtype, class vartype, int Nvar, int Ndeg>
struct polynomial_evaluator
{
    static vartype eval(const numtype p[], const vartype x[])
    {
        return polynomial_evaluator<numtype, vartype, Nvar, Ndeg>::eval_monomial(
                   p + binomial<Nvar + Ndeg - 1, Ndeg - 1>::value, x)
            + polynomial_evaluator<numtype, vartype, Nvar, Ndeg - 1>::eval(p, x);
    }

    // Evaluate monomial in Nvar variables.
    // M(N,K) = x[0]*M(N-1,K) + M(N,K-1)
    static vartype eval_monomial(const numtype p[], const vartype x[])
    {
        return x[0] * polynomial_evaluator<numtype, vartype, Nvar, Ndeg - 1>::eval_monomial(p, x)
            + polynomial_evaluator<numtype, vartype, Nvar - 1, Ndeg>::eval_monomial(
                   p + binomial<Nvar + Ndeg - 2, Ndeg - 1>::value, x + 1);
    }
    // cn is the monomial Nvar,Ndeg, so
    // cp is the monomial Nvar,Ndeg-1
    static void eval_terms_helper(vartype cn[], const vartype cp[], const vartype x[])
    {
        for (int i = 0; i < binomial<Nvar + Ndeg - 2, Ndeg - 1>::value; i++)
            cn[i] = cp[i] * x[0];
        polynomial_evaluator<numtype, vartype, Nvar - 1, Ndeg>::eval_terms_helper(
            cn + binomial<Nvar + Ndeg - 2, Ndeg - 1>::value,
            cp + binomial<Nvar + Ndeg - 3, Ndeg - 2>::value, x + 1);
    }
    // Evaluate all terms in a Nvar,Ndeg polymomial, given
    // the values of x1, x2 ..
    static void eval_terms(vartype p[], const vartype x[])
    {
        polynomial_evaluator<numtype, vartype, Nvar, Ndeg - 1>::eval_terms(p, x);
        polynomial_evaluator<numtype, vartype, Nvar, Ndeg>::eval_terms_helper(
            p + binomial<Nvar + Ndeg - 1, Ndeg - 1>::value,
            p + binomial<Nvar + Ndeg - 2, Ndeg - 2>::value, x);
    }
};

template <class numtype, class vartype, int Nvar>
struct polynomial_evaluator<numtype, vartype, Nvar, 1>
{
    /*
     * eval
     */
    static vartype eval(const numtype p[], const vartype x[])
    {
        return polynomial_evaluator<numtype, vartype, Nvar, 1>::eval_monomial(
                   p + binomial<Nvar + 1 - 1, 1 - 1>::value, x)
            + polynomial_evaluator<numtype, vartype, Nvar, 1 - 1>::eval(p, x);
    }

    /*
     * Evaluate monomial in Nvar variables. M(N,K) = x[0]*M(N-1,K) + M(N,K-1)
     */
    static vartype eval_monomial(const numtype p[], const vartype x[])
    {
        return x[0] * polynomial_evaluator<numtype, vartype, Nvar, 1 - 1>::eval_monomial(p, x)
            + polynomial_evaluator<numtype, vartype, Nvar - 1, 1>::eval_monomial(
                   p + binomial<Nvar + 1 - 2, 1 - 1>::value, x + 1);
    }

    /*
     * eval_terms_helper
     */
    static void eval_terms_helper(vartype cn[], const vartype cp[], const vartype x[])
    {
        for (int i = 0; i < Nvar; i++)
            cn[i] = x[i];
    }

    /*
     * Evaluate all terms in a Nvar,1 polymomial, given the values of x1, x2 ..
     */
    static void eval_terms(vartype p[], const vartype x[])
    {
        p[0] = 1;
        for (int i = 0; i < Nvar; i++)
            p[i + 1] = x[i];
    }
};

/*
 * 1 x y x^2 xy y^2 =
 *   1 [x] [y] x[x y] y^2
 */
template <class numtype, class vartype, int Ndeg>
struct polynomial_evaluator<numtype, vartype, 1, Ndeg>
{
    /*
     * a + bx + cx^2 = a + x(b + x(c))
     */
    static vartype eval(const numtype p[], const vartype x[])
    {
        // Horner scheme:
        vartype sum = p[Ndeg];
        for (int i = Ndeg - 1; i >= 0; i--)
            sum = sum * x[0] + p[i];
        return sum;
    }

    /*
     * eval_monomial
     */
    static vartype eval_monomial(const numtype p[], const vartype x[])
    {
        vartype xn = x[0];
        for (int i = 1; i < Ndeg; i++)
            xn *= x[0];
        return xn * p[0];
    }

    /*
     * eval_terms_helper
     */
    static void eval_terms_helper(vartype cn[], const vartype cp[], const vartype x[])
    {
        cn[0] = cp[0] * x[0];
    }

    /*
     * eval_terms
     */
    static void eval_terms(vartype p[], const vartype x[])
    {
        p[0] = 1;
        if (Ndeg > 0)
            p[1] = x[0];
        for (int i = 2; i <= Ndeg; i++)
            p[i] = x[0] * p[i - 1];
    }
};

template <class numtype, class vartype>
struct polynomial_evaluator<numtype, vartype, 1, 1>
{
    /*
     * a + bx + cx^2 = a + x(b + x(c))
     */
    static vartype eval(const numtype p[], const vartype x[])
    {
        // Horner scheme:
        vartype sum = p[1];
        for (int i = 1 - 1; i >= 0; i--)
            sum = sum * x[0] + p[i];
        return sum;
    }

    /*
     * eval_monomial
     */
    static vartype eval_monomial(const numtype p[], const vartype x[])
    {
        return x[0] * p[0];
    }

    /*
     * eval_terms_helper
     */
    static void eval_terms_helper(vartype cn[], const vartype cp[], const vartype x[])
    {
        cn[0] = cp[0] * x[0];
    }

    /*
     * Evaluate all terms in a 1,1 polymomial, given the values of x1, x2 ..
     */
    static void eval_terms(vartype p[], const vartype x[])
    {
        p[0] = 1;
        p[1] = x[0];
    }
};

template <class numtype, class vartype, int Nvar>
struct polynomial_evaluator<numtype, vartype, Nvar, 0>
{
    /*
     * eval
     */
    static vartype eval(const numtype p[], const vartype x[])
    {
        return p[0];
    }

    /*
     * eval_monomial
     */
    static vartype eval_monomial(const numtype p[], const vartype x[])
    {
        return p[0];
    }

    /*
     * eval_terms_helper
     */
    static void eval_terms_helper(vartype cn[], const vartype cp[], const vartype x[])
    {
        assert(0 && "BUG");
    }

    /*
     * eval_terms
     */
    static void eval_terms(vartype p[], const vartype x[])
    {
        p[0] = 1;
    }
};

/*
 * new_evaluator
 */
template <class numtype, class vartype, int Nvar, int Ndeg, int J, int N, int term>
struct new_evaluator
{
    static vartype eval(vartype tmp[Ndeg], const numtype x[], const numtype p[])
    {
        tmp[N - 1] = tmp[N - 2] * x[J];
        return tmp[N - 1] * p[term]
            + new_evaluator<numtype, vartype, Nvar, Ndeg, J, N + 1, term_prod<Nvar, J + 1, term>::prod>::
                  eval(tmp, x, p)
            + new_evaluator<numtype, vartype, Nvar, Ndeg, J + 1, N, term + 1>::eval(tmp, x, p);
    }
};

template <class numtype, class vartype, int Nvar, int Ndeg>
struct new_evaluator<numtype, vartype, Nvar, Ndeg, 0, 0, 0>
{
    static vartype eval(vartype tmp[Ndeg], const numtype x[], const numtype p[])
    {
        return p[0] + new_evaluator<numtype, vartype, Nvar, Ndeg, 0, 1, 1>::eval(tmp, x, p);
    }
};

template <class numtype, class vartype, int Nvar, int Ndeg, int J, int term>
struct new_evaluator<numtype, vartype, Nvar, Ndeg, J, 1, term>
{
    static vartype eval(vartype tmp[Ndeg], const numtype x[], const numtype p[])
    {
        tmp[0] = x[J];
        return tmp[0] * p[term]
            + new_evaluator<numtype, vartype, Nvar, Ndeg, J, 2, term_prod<Nvar, J + 1, term>::prod>::
                  eval(tmp, x, p)
            + new_evaluator<numtype, vartype, Nvar, Ndeg, J + 1, 1, term + 1>::eval(tmp, x, p);
    }
};

template <class numtype, class vartype, int Nvar, int Ndeg, int J, int term>
struct new_evaluator<numtype, vartype, Nvar, Ndeg, J, Ndeg, term>
{
    static vartype eval(vartype tmp[Ndeg], const numtype x[], const numtype p[])
    {
        tmp[Ndeg - 1] = tmp[Ndeg - 2] * x[J];
        return tmp[Ndeg - 1] * p[term]
            + new_evaluator<numtype, vartype, Nvar, Ndeg, J + 1, Ndeg, term + 1>::eval(tmp, x, p);
        return 0;
    }
};

template <class numtype, class vartype, int Nvar, int Ndeg, int N, int term>
struct new_evaluator<numtype, vartype, Nvar, Ndeg, Nvar, N, term>
{
    static vartype eval(vartype tmp[Ndeg], const numtype x[], const numtype p[])
    {
        return 0;
    }
};

template <class numtype, class vartype, int Nvar, int Ndeg, int term>
struct new_evaluator<numtype, vartype, Nvar, Ndeg, Nvar, Ndeg, term>
{
    static vartype eval(vartype tmp[Ndeg], const numtype x[], const numtype p[])
    {
        return 0;
    }
};

template <class numtype, class vartype, int Nvar, int Ndeg, int term>
struct new_evaluator<numtype, vartype, Nvar, Ndeg, Nvar, 1, term>
{
    static vartype eval(vartype tmp[Ndeg], const numtype x[], const numtype p[])
    {
        return 0;
    }
};

/*
 * Linear transformation of polynomials
 * Optimization opportunities:
 *   - do not put the highest order termporaries in tmp,
 *   - directly add them to dst with correct weight
 */
template <class numtype, int Nvar_dst, int Nvar_src, int Ndeg, int J, int N, int term>
struct transformer
{
    static void trans(
        numtype tmp[], numtype dst[], const numtype src[], const numtype T[Nvar_dst * Nvar_src])
    {
        // Calculat the N:th order monomial of tmp for this J value
        detail::polynomial_multiplier<numtype, Nvar_dst - 1, N - 1, 1>::mul_set(
            tmp + polylen<Nvar_dst, N - 1>::len, tmp + polylen<Nvar_dst, N - 2>::len,
            T + J * Nvar_dst);
        if (J == 0)
            for (int i = polylen<Nvar_dst, N - 1>::len; i < polylen<Nvar_dst, N>::len; i++)
                dst[i] = src[term] * tmp[i];
        else
            for (int i = polylen<Nvar_dst, N - 1>::len; i < polylen<Nvar_dst, N>::len; i++)
                dst[i] += src[term] * tmp[i];

        transformer<numtype, Nvar_dst, Nvar_src, Ndeg, J, N + 1, term_prod<Nvar_src, J + 1, term>::prod>::
            trans(tmp, dst, src, T);
        transformer<numtype, Nvar_dst, Nvar_src, Ndeg, J + 1, N, term + 1>::trans(tmp, dst, src, T);
    }
};

template <class numtype, int Nvar_dst, int Nvar_src, int Ndeg>
struct transformer<numtype, Nvar_dst, Nvar_src, Ndeg, 0, 0, 0>
{
    static void trans(
        numtype tmp[], numtype dst[], const numtype src[], const numtype T[Nvar_dst * Nvar_src])
    {
        dst[0] = src[0];
        transformer<numtype, Nvar_dst, Nvar_src, Ndeg, 0, 1, 1>::trans(tmp, dst, src, T);
    }
};

template <class numtype, int Nvar_dst, int Nvar_src>
struct transformer<numtype, Nvar_dst, Nvar_src, 0, 0, 0, 0>
{
    static void trans(
        numtype tmp[], numtype dst[], const numtype src[], const numtype T[Nvar_dst * Nvar_src])
    {
        dst[0] = src[0];
    }
};

/*
 * This case is here to make template resolution unambiguous for Ndeg=1.
 */
template <class numtype, int Nvar_dst, int Nvar_src>
struct transformer<numtype, Nvar_dst, Nvar_src, 1, 0, 0, 0>
{
    static void trans(
        numtype tmp[], numtype dst[], const numtype src[], const numtype T[Nvar_dst * Nvar_src])
    {
        dst[0] = src[0];
        for (int i = 0; i < Nvar_dst; i++)
            dst[1 + i] = src[1] * T[i];
        for (int j = 1; j < Nvar_src; j++)
            for (int i = 0; i < Nvar_dst; i++)
                dst[1 + i] += src[1 + j] * T[i + Nvar_dst * j];
    }
};

template <class numtype, int Nvar_dst, int Nvar_src, int Ndeg, int J, int term>
struct transformer<numtype, Nvar_dst, Nvar_src, Ndeg, J, 1, term>
{
    static void trans(
        numtype tmp[], numtype dst[], const numtype src[], const numtype T[Nvar_dst * Nvar_src])
    {
        for (int i = 0; i < Nvar_dst; i++)
            tmp[1 + i] = T[i + Nvar_dst * J];
        if (J == 0)
            for (int i = polylen<Nvar_dst, 0>::len; i < polylen<Nvar_dst, 1>::len; i++)
                dst[i] = src[term] * tmp[i];
        else
            for (int i = polylen<Nvar_dst, 0>::len; i < polylen<Nvar_dst, 1>::len; i++)
                dst[i] += src[term] * tmp[i];

        transformer<numtype, Nvar_dst, Nvar_src, Ndeg, J, 2, term_prod<Nvar_src, J + 1, term>::prod>::
            trans(tmp, dst, src, T);
        transformer<numtype, Nvar_dst, Nvar_src, Ndeg, J + 1, 1, term + 1>::trans(tmp, dst, src, T);
    }
};

template <class numtype, int Nvar_dst, int Nvar_src, int Ndeg, int J, int term>
struct transformer<numtype, Nvar_dst, Nvar_src, Ndeg, J, Ndeg, term>
{
    static void trans(
        numtype tmp[], numtype dst[], const numtype src[], const numtype T[Nvar_dst * Nvar_src])
    {
        detail::polynomial_multiplier<numtype, Nvar_dst - 1, Ndeg - 1, 1>::mul_set(
            tmp + polylen<Nvar_dst, Ndeg - 1>::len, tmp + polylen<Nvar_dst, Ndeg - 2>::len,
            T + J * Nvar_dst);
        if (J == 0)
            for (int i = polylen<Nvar_dst, Ndeg - 1>::len; i < polylen<Nvar_dst, Ndeg>::len; i++)
                dst[i] = src[term] * tmp[i];
        else
            for (int i = polylen<Nvar_dst, Ndeg - 1>::len; i < polylen<Nvar_dst, Ndeg>::len; i++)
                dst[i] += src[term] * tmp[i];

        transformer<numtype, Nvar_dst, Nvar_src, Ndeg, J + 1, Ndeg, term + 1>::trans(
            tmp, dst, src, T);
    }
};

template <class numtype, int Nvar_dst, int Nvar_src, int Ndeg, int N, int term>
struct transformer<numtype, Nvar_dst, Nvar_src, Ndeg, Nvar_src, N, term>
{
    static void trans(
        numtype tmp[], numtype dst[], const numtype src[], const numtype T[Nvar_dst * Nvar_src])
    {
    }
};

template <class numtype, int Nvar_dst, int Nvar_src, int Ndeg, int term>
struct transformer<numtype, Nvar_dst, Nvar_src, Ndeg, Nvar_src, Ndeg, term>
{
    static void trans(
        numtype tmp[], numtype dst[], const numtype src[], const numtype T[Nvar_dst * Nvar_src])
    {
    }
};

template <class numtype, int Nvar_dst, int Nvar_src, int Ndeg, int term>
struct transformer<numtype, Nvar_dst, Nvar_src, Ndeg, Nvar_src, 1, term>
{
    static void trans(
        numtype tmp[], numtype dst[], const numtype src[], const numtype T[Nvar_dst * Nvar_src])
    {
    }
};

template <class numtype, int Nvar, int Ndeg, int var>
struct differentiator
{
    /*
     * diff
     */
    static void diff(numtype dst[], const numtype src[])
    {
        differentiator<numtype, Nvar, Ndeg - 1, var>::diff(dst, src);
        differentiator<numtype, Nvar - 1, Ndeg, var - 1>::diff(
            dst + polylen<Nvar, Ndeg - 2>::len, src + polylen<Nvar, Ndeg - 1>::len);
    }

    /*
     * Multiply the highest part of dst with the appropriate factors, assuming it comes from a
     * differentiation wrt variable 0.
     *
     * var here indicated the factor, so first multiply the first part of the highest part of dst
     * with var, then the rest with lower values of var.
     *
     * The top monomial starts at polylen<Nvar,Ndeg-1>. It has length polylen<Nvar-1,Ndeg>
     * The x^m part of it has length polylen<Nvar-1,Ndeg-m>
     */
    static void dfacs(numtype dst[])
    {
        for (int i = polylen<Nvar, Ndeg - 1>::len + polylen<Nvar - 1, Ndeg - var - 1>::len;
             i < polylen<Nvar, Ndeg - 1>::len + polylen<Nvar - 1, Ndeg - var>::len; i++)
            dst[i] *= var + 1;
        differentiator<numtype, Nvar, Ndeg, var - 1>::dfacs(dst);
    }
};

template <class numtype, int Nvar, int var>
struct differentiator<numtype, Nvar, 0, var>
{
    static void diff(numtype dst[], const numtype src[])
    {
        dst[0] = numtype();
    }
    static void dfacs(numtype dst[])
    {
    }
};

template <class numtype, int Nvar>
struct differentiator<numtype, Nvar, 0, 0>
{
    static void diff(numtype dst[], const numtype src[])
    {
        dst[0] = numtype();
    }
    static void dfacs(numtype dst[])
    {
    }
};

template <class numtype, int Nvar, int Ndeg>
struct differentiator<numtype, Nvar, Ndeg, 0>
{
    /*
     * diff
     */
    static void diff(numtype dst[], const numtype src[])
    {
        differentiator<numtype, Nvar, Ndeg - 1, 0>::diff(dst, src);
        // copy first part of highest monomial to next highest
        for (int i = 0; i < polylen<Nvar - 1, Ndeg - 1>::len; i++)
            dst[polylen<Nvar, Ndeg - 2>::len + i] = src[polylen<Nvar, Ndeg - 1>::len + i];
        // Multiply by (old) x exponents
        differentiator<numtype, Nvar, Ndeg - 1, Ndeg - 1>::dfacs(dst);
    }

    /*
     * Multiply the highest part of dst with the appropriate factors
     */
    static void dfacs(numtype dst[])
    {
        // Multiply by 1, a NOP
    }
};

template <class numtype, int Ndeg>
struct differentiator<numtype, 1, Ndeg, 0>
{
    /*
     * diff
     */
    static void diff(numtype dst[], const numtype src[])
    {
        for (int i = 0; i < Ndeg; i++)
            dst[i] = (i + 1) * src[i + 1];
    }

    /*
     * dfacs
     */
    static void dfacs(numtype dst[])
    {
        assert(0 && "BUG");
    }
};

/*!
 * \brief Length of a nvar,neg polynomial = binomial(nvar+ndeg,ndeg). Evaluated at run time.
 * \param nvar
 * \param ndeg
 * \return
 */
inline int runtime_polylen(int nvar, int ndeg)
{
    int len = 1;
    for (int k = 1; k <= nvar; k++)
    {
        len *= ndeg + k;
        len /= k;
    }
    return len;
}

#ifdef POLYMUL_USE_TABLES
#include "polymul_tab.h"
#endif

} // End of namespace detail

/*!
 * \brief The basic_var_polynomial class
 */
template <class numtype, int Nvar, int Ndeg>
class basic_var_polynomial
{
public:
    /*!
     * \brief Number of coefficients.
     */
    enum
    {
        size = detail::polylen<Nvar, Ndeg>::len
    };

    /*!
     * \brief Coefficient array type
     */
    typedef numtype coeff_array[size];

    /*!
     * \brief operator []
     * \param i
     * \return
     */
    virtual numtype &operator[](int i) = 0;

    /*!
     * \brief get_coeffs
     * \return
     */
    virtual coeff_array &get_coeffs() = 0;
};

/*!
 * \brief Polynomial with Nvar variables and degree Ndeg. numtype is the type of the coefficients.
 * This version wraps a reference to a const coefficient array.
 */
template <class numtype, int Nvar, int Ndeg>
class const_polynomial
{
public:
    /*!
     * \brief Number of coefficients.
     */
    enum
    {
        size = detail::polylen<Nvar, Ndeg>::len
    };

    /*!
     * \brief Coefficient array type
     */
    typedef numtype coeff_array[size];

    /*!
     * \brief Construct a polynomial with coefficient array reference
     */
    const_polynomial(const coeff_array &coeffs) : c(coeffs)
    {
    }

    /*!
     * \brief Return coefficient i, in graded lexicographical order, for example 1 x y x^2 xy y^2 ..
     * \param i
     * \return
     */
    const numtype &operator[](int i) const
    {
        assert(i >= 0);
        assert(i < this->size);
        return c[i];
    }

    /*!
     * \brief get_coeffs
     * \return
     */
    const coeff_array &get_coeffs() const
    {
        return c;
    }

    /*!
     * \brief Convert to a polynomial of different degree and/or different scalar type.
     * Zeros are inserted if the new polynomial has higher degree than this.
     * \param dst
     * \todo Allow also change of Nvar
     */
    template <class T, int Ndeg2>
    void convert_to(basic_var_polynomial<T, Nvar, Ndeg2> &result) const
    {
        if (size < int(result.size))
        {
            for (int i = 0; i < size; i++)
                result[i] = c[i];
            for (int i = size; i < result.size; i++)
                result[i] = T();
        }
        else
        {
            for (int i = 0; i < result.size; i++)
                result[i] = c[i];
        }
    }

    /*!
     * \brief Calculate the exponents of the given term.
     * \note This is a slow function.
     * \param term
     * \param exponents
     */
    static void exponents(int term, int exponents[Nvar])
    {
        assert(term >= 0);
        if (Nvar == 1)
        {
            exponents[0] = term;
            return;
        }
        if (term >= size)
            assert(0 && "term > size");
        for (int i = 0; i < Nvar; i++)
            exponents[i] = 0;
        for (int i = 0; i < term; i++)
            const_polynomial<numtype, Nvar, Ndeg>::next_exponents(Nvar, exponents);
    }

    /*!
     * \brief Return the index of the term with the given exponents.
     * \param exponents
     * \return
     */
    static int term_index(const int exponents[Nvar])
    {
        int N = 0;
        for (int i = 0; i < Nvar; i++)
            N += exponents[i];
        int i = 0, idx = 0;
        N--;
        while (N >= 0)
        {
            idx += detail::runtime_polylen(Nvar - i, N);
            N -= exponents[i];
            i++;
        }
        return idx;
    }

    /*!
     * \brief Evaluate the polynomial at x.
     * \param x
     * \return
     */
    template <class T>
    T eval(const T x[Nvar]) const
    {
        return detail::polynomial_evaluator<numtype, T, Nvar, Ndeg>::eval(c, x);
    }

    /*!
     * \brief Calculate the derivative of this with respect to var.
     * \param derivative
     */
    template <int var>
    void diff(basic_var_polynomial<numtype, Nvar, Ndeg - 1> &derivative) const
    {
        detail::differentiator<numtype, Nvar, Ndeg, var>::diff(derivative.get_coeffs(), c);
    }

protected:
    /*!
     * \brief Next exponents
     * \param nvar
     * \param m
     */
    static void next_exponents(int nvar, int m[Nvar])
    {
        int k = 0;
        for (int i = 0; i < nvar - 1; i++)
            k += m[i];
        if (k == 0)
        {
            m[0] = m[nvar - 1] + 1;
            m[nvar - 1] = 0;
            return;
        }
        if (m[nvar - 2] > 0)
        {
            m[nvar - 1]++;
            m[nvar - 2]--;
        }
        else
        {
            next_exponents(nvar - 1, m);
            for (int i = nvar - 2; i >= 0; i--)
            {
                if (m[i] > 0)
                {
                    m[i] += m[nvar - 1];
                    break;
                }
            }
            m[nvar - 1] = 0;
        }
    }

private:
    /*!
     * \brief Coefficients
     */
    const coeff_array &c;
};

/*!
 * \brief Polynomial with Nvar variables and degree Ndeg. numtype is the type of the coefficients.
 * This version wraps a reference to a non-const coefficient array.
 */
template <class numtype, int Nvar, int Ndeg>
class var_polynomial : public basic_var_polynomial<numtype, Nvar, Ndeg>,
                       public const_polynomial<numtype, Nvar, Ndeg>
{
public:
    /*!
     * \brief Number of coefficients
     */
    using const_polynomial<numtype, Nvar, Ndeg>::size;

    /*!
     * \brief Coefficient array type
     */
    using typename const_polynomial<numtype, Nvar, Ndeg>::coeff_array;

    /*!
     * \brief Construct a polynomial with coefficient array reference
     */
    var_polynomial(coeff_array &coeffs) : const_polynomial<numtype, Nvar, Ndeg>(coeffs), c(coeffs)
    {
    }

    /*!
     * \brief Return coefficient i, in graded lexicographical order, for example 1 x y x^2 xy y^2 ..
     * \param i
     * \return
     */
    virtual numtype &operator[](int i)
    {
        assert(i >= 0);
        assert(i < this->size);
        return c[i];
    }

    /*!
     * \brief get_coeffs
     * \return
     */
    virtual coeff_array &get_coeffs()
    {
        return c;
    }

    /*!
     * \brief Copy assigment
     * \param other
     * \return
     */
    var_polynomial &operator=(const var_polynomial &other)
    {
        for (int i = 0; i < size; i++)
            c[i] = other.c[i];
        return *this;
    }

    /*!
     * \brief Assign the coefficients
     * \param coeffs
     * \return
     */
    var_polynomial &operator=(const coeff_array &coeffs)
    {
        for (int i = 0; i < size; i++)
            c[i] = coeffs[i];
        return *this;
    }

    /*!
     * \brief Assign the constant term the value c0 and zero the other terms
     * \param c0
     * \return
     */
    template <class T>
    var_polynomial &operator=(const T &c0)
    {
        c[0] = numtype(c0);
        for (int i = 1; i < size; i++)
            c[i] = numtype(0);
        return *this;
    }

private:
    /*!
     * \brief Coefficients
     */
    coeff_array &c;
};

/*!
 * \brief Polynomial with Nvar variables and degree Ndeg. numtype is the type of the coefficients.
 */
template <class numtype, int Nvar, int Ndeg>
class polynomial : public var_polynomial<numtype, Nvar, Ndeg>
{
public:
    /*!
     * \brief Number of coefficients
     */
    using var_polynomial<numtype, Nvar, Ndeg>::size;

    /*!
     * \brief Coefficient array type
     */
    using typename var_polynomial<numtype, Nvar, Ndeg>::coeff_array;

    /*!
     * \brief Construct a polynomial with constant coefficient c0 and other terms zero.
     * \param c0
     */
    explicit polynomial(const numtype &c0 = numtype()) : var_polynomial<numtype, Nvar, Ndeg>(c)
    {
        c[0] = c0;
        for (int i = 1; i < size; i++)
            c[i] = numtype();
    }

    /*!
     * \brief Using assigment operator from base class
     */
    using var_polynomial<numtype, Nvar, Ndeg>::operator=;

private:
    /*!
     * \brief Coefficients
     */
    coeff_array c;
};

/*!
 * \brief Multiplication
 * \param dst
 * \param p1
 * \param p2
 */
template <class numtype, int Nvar, int Ndeg1, int Ndeg2>
void mul(
    var_polynomial<numtype, Nvar, Ndeg1 + Ndeg2> &POLYMUL_RESTRICT dst,
    const const_polynomial<numtype, Nvar, Ndeg1> &p1,
    const const_polynomial<numtype, Nvar, Ndeg2> &p2)
{
    detail::polynomial_multiplier<numtype, Nvar, Ndeg1, Ndeg2>::mul_set(
        dst.get_coeffs(), p1.get_coeffs(), p2.get_coeffs());
}

/*!
 * \brief Taylor multiplication
 * \param dst
 * \param p1
 * \param p2
 */
template <class numtype, int Nvar, int Ndeg>
void taylormul(
    var_polynomial<numtype, Nvar, Ndeg> &POLYMUL_RESTRICT dst,
    const const_polynomial<numtype, Nvar, Ndeg> &p1,
    const const_polynomial<numtype, Nvar, Ndeg> &p2)
{
    detail::taylor_multiplier<numtype, Nvar, Ndeg, Ndeg>::mul_set(
        dst.get_coeffs(), p1.get_coeffs(), p2.get_coeffs());
}

/*!
 * \brief Taylor in-place multiplication
 * \param dst
 * \param p
 */
template <class numtype, int Nvar, int Ndeg, int Ndeg2>
void taylormul(
    var_polynomial<numtype, Nvar, Ndeg> &POLYMUL_RESTRICT dst,
    const const_polynomial<numtype, Nvar, Ndeg2> &p)
{
    assert(Ndeg2 <= Ndeg);
    detail::taylor_inplace_multiplier<numtype, Nvar, Ndeg, Ndeg2, 0>::mul(
        dst.get_coeffs(), p.get_coeffs());
}

/*!
 * \brief _Add_ the product of p and the single polynomial term rterm with coefficient c to dst
 * \param dst
 * \param p
 * \param c
 */
template <class numtype, int Nvar, int Ndeg, int rterm>
void termmul(
    var_polynomial<numtype, Nvar, Ndeg + detail::term_deg<Nvar, rterm>::deg> &POLYMUL_RESTRICT dst,
    const const_polynomial<numtype, Nvar, Ndeg> &p, const numtype &c)
{
    detail::term_multiplier<numtype, Nvar, detail::polylen<Nvar, Ndeg>::len - 1, rterm>::mul(
        dst.get_coeffs(), p.get_coeffs(), c);
}

/*!
 * \brief Number of coefficients
 * \param p
 * \param x
 */
template <class numtype, int Nvar, int Ndeg>
void terms(var_polynomial<numtype, Nvar, Ndeg> &p, const numtype x[])
{
    detail::polynomial_evaluator<numtype, numtype, Nvar, Ndeg>::eval_terms(p.get_coeffs(), x);
}

/*!
 * \brief Calculates p1 so that dot(P*p2,p3) = dot(P,p1)
 * P*p2 is polynomial multiplication, dot means summing the product of all coefficients).
 * \param p1
 * \param p2
 * \param p3
 */
template <class numtype, int Nvar, int Ndeg1, int Ndeg2>
void contract(
    var_polynomial<numtype, Nvar, Ndeg1> &POLYMUL_RESTRICT p1,
    const const_polynomial<numtype, Nvar, Ndeg2> &p2,
    const const_polynomial<numtype, Nvar, Ndeg1 + Ndeg2> &p3)
{
    for (int i = 0; i < p1.size; i++)
        p1[i] = numtype();
    detail::polynomial_multiplier<numtype, Nvar, Ndeg1, Ndeg2>::antimul(
        p3.get_coeffs(), p1.get_coeffs(), p2.get_coeffs());
}

/*!
 * \brief Multiply each term in the polynomial with appropriate factorials to generate derivatives,
 * i.e. x^n y^m z^k would be multiplied by n!m!k!
 * \note The factors are generated as enums, limiting the numeric range to rather small values.
 * \param p
 */
template <class numtype, int Nvar, int Ndeg>
void dfac(const_polynomial<numtype, Nvar, Ndeg> &p)
{
    detail::deriv_fac_multiplier<numtype, Nvar, detail::polylen<Nvar, Ndeg>::len - 1>::mul_fac(
        p.get_coeffs());
}

/*!
 * \brief Perform a linear change of variables, with x_i in src being equal to the linear
 * combination as x_i = sum_j T_ji y_j. T_ji is of the format T_00, T_10, T_20 ..
 * \param dst
 * \param src
 * \param T
 */
template <class numtype, int Nvar_src, int Nvar_dst, int Ndeg>
void trans(
    var_polynomial<numtype, Nvar_dst, Ndeg> &dst,
    const const_polynomial<numtype, Nvar_src, Ndeg> &src, const numtype T[Nvar_src * Nvar_dst])
{
    polynomial<numtype, Nvar_dst, Ndeg> tmp;
    detail::transformer<numtype, Nvar_dst, Nvar_src, Ndeg, 0, 0, 0>::trans(
        tmp.get_coeffs(), dst.get_coeffs(), src.get_coeffs(), T);
}

/*!
 * \brief mul_vec
 * \param size
 * \param dst
 * \param p1
 * \param p2
 */
template <class numtype, int Nvar, int Ndeg1, int Ndeg2>
inline void mul_vec(
    int size, polynomial<numtype, Nvar, Ndeg1 + Ndeg2> dst[],
    const polynomial<numtype, Nvar, Ndeg1> p1[], const polynomial<numtype, Nvar, Ndeg2> p2[])
{
    for (int i = 0; i < size; i++)
        mul(dst[i], p1[i], p2[i]);
}

/*!
 * \brief taylormul_vec
 * \param size
 * \param dst
 * \param p1
 * \param p2
 */
template <class numtype, int Nvar, int Ndeg>
inline void taylormul_vec(
    int size, polynomial<numtype, Nvar, Ndeg> dst[], const polynomial<numtype, Nvar, Ndeg> p1[],
    const polynomial<numtype, Nvar, Ndeg> p2[])
{
    for (int i = 0; i < size; i++)
        taylormul(dst[i], p1[i], p2[i]);
}

/*!
 * \brief taylormul_vec
 * \param size
 * \param p1
 * \param p2
 */
template <class numtype, int Nvar, int Ndeg, int Ndeg2>
inline void taylormul_vec(
    int size, polynomial<numtype, Nvar, Ndeg> p1[], const polynomial<numtype, Nvar, Ndeg2> p2[])
{
    for (int i = 0; i < size; i++)
        taylormul(p1[i], p2[i]);
}

/*!
 * \brief Print polynomial
 * \param stream
 * \param p
 */
template <class num, int Nvar, int Ndeg>
void print(const const_polynomial<num, Nvar, Ndeg> &p, std::ostream &stream = std::cout)
{
    int exps[Nvar] = {0};
    for (int i = 0; i < p.size; i++)
    {
        p.exponents(i, exps);
        for (int j = 0; j < Nvar; j++)
            stream << exps[j] << " ";
        stream << " " << p[i] << std::endl;
        assert(i == p.term_index(exps));
    }
}

} // namespace polymul

#endif
