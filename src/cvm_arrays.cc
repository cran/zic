//                  CVM Class Library
//                  http://cvmlib.com
//
//          Copyright Sergei Nikolaev 1992-2008
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include "cvm.h"
#include "cvm_blas.h"

#if defined (_MSC_VER)
#   pragma warning(disable:4786)
#endif


CVM_NAMESPACE_BEG

template<>
CVM_API void __scal<double, double>(double* mpD, int mnSize, int mnIncr, double dScal)
{
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(double))
    DSCAL (&mnSize, &dScal, mpD, &mnIncr);
}

template<>
CVM_API double __norm<double, double>(const double* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(double))
    return DNRM2 (&nSize, pD, &nIncr);
}

template<>
CVM_API double __norm<double, std::complex<double> >(const std::complex<double>* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(std::complex<double>))
    return DZNRM2 (&nSize, pD, &nIncr);
}

template<>
CVM_API int __idamax<double>(const double* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(double))
    return IDAMAX (&nSize, pD, &nIncr);
}

template<>
CVM_API int __idamin<double>(const double* pD, int nSize, int nIncr)
{
    CVM_ASSERT(pD, ((nSize - 1) * nIncr + 1) * sizeof(double))
    return IDAMIN (&nSize, pD, &nIncr);
}

template<>
CVM_API void __add<double>(double* mpD, int mnSize, int mnIncr, const double* pv, int nIncr)
{
    static const double one(1.);
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(double))
    CVM_ASSERT(pv, ((mnSize - 1) * nIncr + 1) * sizeof(double))
    DAXPY (&mnSize, &one, pv, &nIncr, mpD, &mnIncr);
}

template<>
CVM_API void __subtract<double>(double* mpD, int mnSize, int mnIncr, const double* pv, int nIncr)
{
    static const double mone(-1.F);
    CVM_ASSERT(mpD, ((mnSize - 1) * mnIncr + 1) * sizeof(double))
    CVM_ASSERT(pv, ((mnSize - 1) * nIncr + 1) * sizeof(double))
    DAXPY (&mnSize, &mone, pv, &nIncr, mpD, &mnIncr);
}

template<>
CVM_API void __randomize<double> (double* mpD, int mnSize, int mnIncr, double dFrom, double dTo)
{
    const int nSize = mnSize * mnIncr;
    for (int i = 0; i < nSize; i += mnIncr)
    {
        mpD[i] = Randomizer<double>::get (dFrom, dTo);
    }
}

CVM_NAMESPACE_END
