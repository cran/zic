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
#   pragma warning(disable:4503)
#endif

CVM_NAMESPACE_BEG

char Chars::mchars[15] = {'T','N','U','L','P','Q','B','E','R','A','S','V','O','I','C'};

// global error messages holder
CVM_API ErrMessages& ErrMessages::ErrMessagesInstance()
{
    static ErrMessages _ErrMessages;
    return _ErrMessages;
}

CVM_API ErrMessages::ErrMessages()
    : msUnknown ("Unknown exception")
{
    mmMsg.insert (pair_Msg (CVM_OK,                      "All OK"));
    mmMsg.insert (pair_Msg (CVM_OUTOFMEMORY,             "Out of memory"));
    mmMsg.insert (pair_Msg (CVM_OUTOFRANGE,              "Index %d is out of range"));
    mmMsg.insert (pair_Msg (CVM_OUTOFRANGE1,             "First index %d is out of range"));
    mmMsg.insert (pair_Msg (CVM_OUTOFRANGE2,             "Second index %d is out of range"));
    mmMsg.insert (pair_Msg (CVM_WRONGSIZE,               "Wrong size %d"));
    mmMsg.insert (pair_Msg (CVM_SIZESMISMATCH,           "Sizes mismatch"));
    mmMsg.insert (pair_Msg (CVM_WRONGMKLARG,             "Wrong argument passed to BLAS or LAPACK subroutine"));
    mmMsg.insert (pair_Msg (CVM_WRONGMKLARG2,            "Wrong argument %d passed to BLAS or LAPACK subroutine %s"));
    mmMsg.insert (pair_Msg (CVM_SINGULARMATRIX,          "The diagonal element (or main minor) %d of the matrix is zero (or singular)"));
    mmMsg.insert (pair_Msg (CVM_NOTPOSITIVEDEFINITE,     "The leading minor of order %d (and hence the matrix itself) is not positive-definite"));
    mmMsg.insert (pair_Msg (CVM_WRONGCHOLESKYFACTOR,     "The diagonal element %d of the Cholesky factor (and hence the factor itself) is zero"));
    mmMsg.insert (pair_Msg (CVM_WRONGBUNCHKAUFMANFACTOR, "The diagonal element %d of the Bunch-Kaufman factor (and hence the factor itself) is zero"));
    mmMsg.insert (pair_Msg (CVM_NOTPOSITIVEDIAG,         "The diagonal element %d of the matrix is nonpositive. Equilibration failed"));
    mmMsg.insert (pair_Msg (CVM_CONVERGENCE_ERROR,       "Method failed to converge"));
    mmMsg.insert (pair_Msg (CVM_DIVISIONBYZERO,          "Division by zero"));

#if defined (WIN32) || defined (_WIN32)
    mmMsg.insert (pair_Msg (CVM_SEMAPHOREERROR,          "Critical Section access error"));
#else
    mmMsg.insert (pair_Msg (CVM_SEMAPHOREERROR,          "Semaphore access error"));
#endif

    mmMsg.insert (pair_Msg (CVM_READ_ONLY_ACCESS,        "Attempt to change a read-only element"));
    mmMsg.insert (pair_Msg (CVM_SUBMATRIXACCESSERROR,    "Attempt to access non-continuous submatrix as continuous array, see programmer\'s reference for further details"));
    mmMsg.insert (pair_Msg (CVM_SUBMATRIXNOTAVAILABLE,   "Submatrix instantiation is not available for class \'%s\', see manual for details"));
    mmMsg.insert (pair_Msg (CVM_MATRIXNOTSYMMETRIC,      "The matrix passed doesn't appear to be symmetric"));
    mmMsg.insert (pair_Msg (CVM_MATRIXNOTHERMITIAN,      "The matrix passed doesn't appear to be hermitian"));
    mmMsg.insert (pair_Msg (CVM_BREAKS_HERMITIANITY,     "This operation could make the matrix non-hermitian. Use %s instead"));
    mmMsg.insert (pair_Msg (CVM_METHODNOTAVAILABLE,      "Function \'%s\' is not available for class \'%s\'. See programmer\'s reference for further details"));
    mmMsg.insert (pair_Msg (CVM_NOTIMPLEMENTED,          "Function \'%s\' is not implemented"));
}

CVM_API const std::string& ErrMessages::_get (int nException)
{
    citr_Msg i = mmMsg.size() > 0 ? mmMsg.find (nException) : mmMsg.end();
    return i == mmMsg.end() ? msUnknown : (*i).second;
}

CVM_API bool ErrMessages::_add (int nNewCause, const char* szNewMessage)
{
    bool bRes = true;
    itr_Msg i = mmMsg.find (nNewCause);
    if (i != mmMsg.end())
    {
        (*i).second = (*i).second + " | " + szNewMessage;       // Defenition is overlapped. This is not a good idea
        bRes = false;                                           // to do so, use CVM_THE_LAST_ERROR_CODE + 1 as an error code.
    }
    else
    {
        mmMsg.insert (pair_Msg (nNewCause, szNewMessage));      // new error definition
    }
    return bRes;
}

template <>
CVM_API void __copy<double> (int nSize, const double* pFrom, int nFromIncr, double* pTo, int nToIncr)
{
    CVM_ASSERT(pFrom, ((nFromIncr) * (nSize - 1) + 1) * sizeof(double))
    CVM_ASSERT(pTo,   ((nToIncr) * (nSize - 1) + 1) * sizeof(double))
    DCOPY (&nSize, pFrom, &nFromIncr, pTo, &nToIncr);
}

template <>
CVM_API void __copy<int> (int nSize, const int* pFrom, int nFromIncr, int* pTo, int nToIncr)
{
    CVM_ASSERT(pFrom, ((nFromIncr) * (nSize - 1) + 1) * sizeof(int))
    CVM_ASSERT(pTo,   ((nToIncr) * (nSize - 1) + 1) * sizeof(int))
    for (int i = 0; i < nSize; ++i)
    {
        pTo[i * nToIncr] = pFrom[i * nFromIncr];
    }
}

template <>
CVM_API void __swap<double> (int nSize, double* p1, int n1Incr, double* p2, int n2Incr)
{
    CVM_ASSERT(p1, ((n1Incr) * (nSize - 1) + 1) * sizeof(double))
    CVM_ASSERT(p2, ((n2Incr) * (nSize - 1) + 1) * sizeof(double))
    DSWAP (&nSize, p1, &n1Incr, p2, &n2Incr);
}

template <>
CVM_API void __swap<int> (int nSize, int* p1, int n1Incr, int* p2, int n2Incr)
{
    int n;
    CVM_ASSERT(p1, (n1Incr * (nSize - 1) + 1) * sizeof(int))
    CVM_ASSERT(p2, (n2Incr * (nSize - 1) + 1) * sizeof(int))
    for (int i = 0; i < nSize; ++i)
    {
        n = p1[i * n1Incr];
        p1[i * n1Incr] = p2[i * n2Incr];
        p2[i * n2Incr] = n;
    }
}

template<>
CVM_API void
__low_up<basic_srmatrix<double> > (basic_srmatrix<double>& m, int* nPivots) throw (cvmexception)
{
    int nOutInfo = 0;
    DGETRF (m._pm(), m._pn(), m, m._pld(), nPivots, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__low_up<basic_srbmatrix<double> >
    (basic_srbmatrix<double>& m, int* nPivots) throw (cvmexception)
{
    int nOutInfo = 0;
    const int nKL = m.lsize();
    const int nKU = m.usize();
    m.resize_lu (nKL, nKL + nKU);
    DGBTRF (m._pm(), m._pn(), &nKL, &nKU, m, m._pld(), nPivots, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API int
__cholesky<basic_srmatrix<double> >
    (basic_srmatrix<double>& m)                           // input is symmetric, output is triangular
{
    int nOutInfo = 0;
    DPOTRF (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            m._pm(), m, m._pld(), &nOutInfo);

    return nOutInfo;
}

template<>
CVM_API void
__bunch_kaufman<basic_srmatrix<double> >
    (basic_srmatrix<double>& m, int* nPivots) throw (cvmexception)        // input is symmetric, output is square
{
    int nOutInfo = 0;
    const int lwork = m.msize() * 64;
    basic_rvector<double> work (lwork);
    DSYTRF (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            m._pm(), m, m._pld(), nPivots, work, &lwork, &nOutInfo);

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_SINGULARMATRIX, nOutInfo);
}

template<>
CVM_API void
__ger<double, basic_rmatrix<double>, basic_rvector<double> >
    (basic_rmatrix<double>& m,
    const basic_rvector<double>& vCol,
    const basic_rvector<double>& vRow,
    double dAlpha)
{
    CVM_ASSERT(m.get(), vCol.size() * vRow.size() * sizeof(double))
    DGER (vCol._psize(), vRow._psize(), &dAlpha, vCol, vCol._pincr(), vRow, vRow._pincr(), m, m._pld());
}

template<>
CVM_API void
__poequ<double, basic_srsmatrix<double>, basic_rvector<double> >
    (const basic_srsmatrix<double>& m,
     basic_rvector<double>& vScalings, 
     double& dCond,
     double& dMax)
{
    int nOutInfo = 0;
    DPOEQU (m._pm(), m, m._pld(), vScalings, &dCond, &dMax, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDIAG, nOutInfo);
}

CVM_NAMESPACE_END
