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
CVM_API double __dot<double> (const double* mpD, int mnSize, int mnIncr, const double* pD, int nIncr)
{
    return DDOT (&mnSize, mpD, &mnIncr, pD, &nIncr);
}

template<>
CVM_API void
__gemv<double, basic_rmatrix<double>, basic_rvector<double> >
    (bool bLeft,
    const basic_rmatrix<double>& m,
    double dAlpha,
    const basic_rvector<double>& v,
    double dBeta,
    basic_rvector<double>& vRes)
{
    DGEMV (bLeft ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), m._pn(), &dAlpha, m._pd(), m._pldm(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__gbmv<double, basic_srbmatrix<double>, basic_rvector<double> >
    (bool bLeft,
    const basic_srbmatrix<double>& m,
    double dAlpha,
    const basic_rvector<double>& v,
    double dBeta,
    basic_rvector<double>& vRes)
{
    DGBMV (bLeft ? Chars::pT() : Chars::pN(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), m._pn(), m._pl(), m._pu(), &dAlpha, m, m._pld(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__symv<double, basic_srsmatrix<double>, basic_rvector<double> >
    (const basic_srsmatrix<double>& m,
    double dAlpha,
    const basic_rvector<double>& v,
    double dBeta,
    basic_rvector<double>& vRes)
{
    DSYMV (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
           1,
#endif
           m._pm(), &dAlpha, m, m._pld(), v, v._pincr(), &dBeta, vRes, vRes._pincr());
}

template<>
CVM_API void
__svd<double, basic_rmatrix<double>, basic_srmatrix<double> >
    (double* pD, int nSize, int nIncr, 
    const basic_rmatrix<double>& mArg,
    basic_srmatrix<double>* mU,
    basic_srmatrix<double>* mVH) throw (cvmexception)
{
    const bool bSimple  = (mU == NULL || mVH == NULL);
    const int  nM       = mArg.msize();
    const int  nN       = mArg.nsize();
    const int  m        = _cvm_min<int>(nM, nN);
    const int  M        = _cvm_max<int>(nM, nN);
    int lWork = -1; // to calculate lWork
    int nOutInfo = 0;

    if (m != nSize) throw cvmexception (CVM_SIZESMISMATCH);

    basic_rvector<double> mD (nSize);
    mD.assign (pD, nIncr);

    basic_rvector<double> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_rvector<double> vTauQ    (m);
    basic_rvector<double> vTauP    (m);
    basic_rmatrix<double> mA       (mArg);
    double dWork;

    DGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    basic_rvector<double> vWork(lWork);
    DGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    static const int zero(0);
    static const int one(1);
    if (bSimple)
    {
        basic_rvector<double> vWork2 (m * 4);
        DBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m,
                &zero, &zero, &zero,
                mD, vOffDiag,
                NULL, &one, NULL, &one, NULL, &one,
                vWork2, &nOutInfo);
    }
    else
    {
        basic_rmatrix<double> Q (mA);
        basic_rmatrix<double> P (mA);

        if (nM > nN) Q.resize(nM, nM);
        if (nM < nN) P.resize(nN, nN);

        {
            lWork = -1;
            DORGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                    1,
#endif
                    &nM, &nM, &nN,
                    Q, &nM,
                    vTauQ,
                    &dWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            lWork = static_cast<int>(dWork);
            basic_rvector<double> vWork3(lWork);
            DORGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                    1,
#endif
                    &nM, &nM, &nN,
                    Q, &nM,
                    vTauQ,
                    vWork3, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        }

        {
            lWork = -1;
            DORGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                    1,
#endif
                    &nN, &nN, &nM,
                    P, &M,
                    vTauP,
                    &dWork, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

            lWork = static_cast<int>(dWork);
            basic_rvector<double> vWork3(lWork);
            DORGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                    1,
#endif
                    &nN, &nN, &nM,
                    P, &M,
                    vTauP,
                    vWork3, &lWork, &nOutInfo);
            if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        }

        basic_rvector<double> vWork2 (m * 4);
        DBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m,
                &nN, &nM,
                &zero,
                mD, vOffDiag,
                P, &M, Q, &nM, NULL, &one,
                vWork2, &nOutInfo);

        (*mU)  = basic_srmatrix<double>(Q.resize(nM, nM));
        (*mVH) = basic_srmatrix<double>(P.resize(nN, nN));
    }

    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    __copy<double> (nSize, mD, mD.incr(), pD, nIncr);
}

template<>
CVM_API void
__svd<double, basic_srbmatrix<double>, basic_srmatrix<double> >
    (double* pD, int nSize, int nIncr,
    const basic_srbmatrix<double>& mArg,
    basic_srmatrix<double>* mU,
    basic_srmatrix<double>* mVH) throw (cvmexception)
{
    static const int zero(0);
    const int m = mArg.msize();
    if (m != nSize) throw cvmexception (CVM_SIZESMISMATCH);

    const bool bSimple  = (mU == NULL || mVH == NULL);
    int nOutInfo = 0;

    basic_rvector<double> mD (nSize);
    mD.assign (pD, nIncr);

    basic_rvector<double>   vOffDiag (_cvm_max<int>(1, m - 1));
    basic_srmatrix<double>  mQ       (bSimple ? 1 : m);
    basic_srmatrix<double>  mPT      (bSimple ? 1 : m);
    basic_srmatrix<double>  mC       (1);
    basic_rvector<double>   vWork    (2 * m);
    basic_srbmatrix<double> mA       (mArg);

    DGBBRD (bSimple ? Chars::pN() : Chars::pB(), 
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            mA._pm(), mA._pn(), &zero, mA._pl(), mA._pu(), mA, mA._pld(),
            mD, vOffDiag,
            mQ, mQ._pm(),
            mPT, mPT._pm(),
            mC, mC._pm(), 
            vWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    basic_rvector<double> vWork2 (m * 4);
    DBDSQR (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            bSimple ? &zero : &m,
            bSimple ? &zero : &m,
            &zero,
            mD, vOffDiag,
            mPT, mPT._pm(),
            mQ, mQ._pm(),
            mC, mC._pm(), 
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    if (!bSimple)
    {
        (*mU)  = mQ;
        (*mVH) = mPT;
    }

    __copy<double> (nSize, mD, mD.incr(), pD, nIncr);
}

static int _ssyevd_lwork (int n)
{
    int k = 1, nn = n - 1;
    while (nn >>= 1) ++k;
    return 3 * n * n + (5 + 2 * k) * n + 1;
}

template<>
CVM_API void
__eig<basic_rvector<double>, basic_srsmatrix<double>, basic_srmatrix<double> >
    (basic_rvector<double>& vRes,
    const basic_srsmatrix<double>& mArg,
    basic_srmatrix<double>* mEigVect,
    bool /*bRightVect*/) throw (cvmexception)
{
    const int nM = mArg.msize();
    if (nM != vRes.size()) throw cvmexception (CVM_SIZESMISMATCH);
    const bool bEigVect = (mEigVect != NULL);

    if (nM == 1)
    {
        vRes[1] = mArg(1,1);
        if (bEigVect)
        {
            static const double one(1.);
            mEigVect -> resize (1);
            (*mEigVect)[1].set (one);
        }
    }
    else
    {
        const char* pcJob = bEigVect ? Chars::pV() : Chars::pN();
        basic_srsmatrix<double> mA (mArg);
        int lwork = bEigVect ? _ssyevd_lwork (nM) : (2 * nM + 1);               // LAPACK wants (1 + 6*N + 2*N**2) though
        int liwork = bEigVect ? (5 * nM + 3) : 1;                               // MKL wants 5*n + 2 though
        basic_rvector<double> work (lwork);
        basic_array<int> iwork (liwork);
        int nOutInfo = 0;

        DSYEVD (pcJob,
    #if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
    #endif
                Chars::pU(),
    #if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
    #endif
                &nM, mA, mA._pld(), vRes, work, &lwork, iwork, &liwork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

        if (bEigVect)
        {
            (*mEigVect) << mA;
        }
    }
}

static int _cheevd_lwork (int n)
{
    int k = 1, nn = n - 1;
    while (nn >>= 1) ++k;
    return 3 * n * n + (4 + 2 * k) * n + 1;
}

template<>
CVM_API void
__pinv<double, basic_rmatrix<double>, basic_rmatrix<double> >
    (basic_rmatrix<double>& mX, 
    const basic_rmatrix<double>& mArg, double threshold) throw (cvmexception)
{
    const int nM = mArg.msize();
    const int nN = mArg.nsize();
    const int m  = _cvm_min<int>(nM, nN);
    int lWork    = -1; // to calculate lWork
    int nOutInfo = 0;

    basic_rvector<double> mD       (m);
    basic_rvector<double> vOffDiag (_cvm_max<int>(1, m - 1));
    basic_rvector<double> vTauQ    (m);
    basic_rvector<double> vTauP    (m);
    basic_rmatrix<double> mA       (mArg);
    double dWork;

    DGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, &dWork, &lWork, &nOutInfo);
    lWork = static_cast<int>(dWork);
    basic_rvector<double> vWork(lWork);
    DGEBRD (&nM, &nN, mA, mA._pld(), mD, vOffDiag, vTauQ, vTauP, vWork, &lWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    static const int zero(0);
    static const int one(1);
    
    // few words about economy:
    // for m > n case we care about m-by-n matrix U and n-by-n matrix VH
    // for m < n case we care about m-by-m matrix U and m-by-n matrix VH
    // however, the whole matrix A is needed to start computations
    basic_rmatrix<double> Q (mA);
    basic_rmatrix<double> P (mA);

    {
        lWork = -1;
        DORGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &m, &nN,
                Q, Q._pld(),
                vTauQ,
                &dWork, &lWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        lWork = static_cast<int>(dWork);
        basic_rvector<double> vWork3(lWork);
        DORGBR (Chars::pQ(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &nM, &m, &nN,
                Q, Q._pld(),
                vTauQ,
                vWork3, &lWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    }

    {
        lWork = -1;
        DORGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m, &nN, &nM,
                P, P._pld(),
                vTauP,
                &dWork, &lWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

        lWork = static_cast<int>(dWork);
        basic_rvector<double> vWork3(lWork);
        DORGBR (Chars::pP(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
                1,
#endif
                &m, &nN, &nM,
                P, P._pld(),
                vTauP,
                vWork3, &lWork, &nOutInfo);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    }

    basic_rvector<double> vWork2 (m * 4);
    DBDSQR (nM >= nN ? Chars::pU() : Chars::pL(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            &nN, &nM,
            &zero,
            mD, vOffDiag,
            P, P._pld(), Q, Q._pld(), NULL, &one,
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    if (nM > nN) P.resize(nN, nN);   // VH
    if (nM < nN) Q.resize(nM, nM);   // U
    for (int i = 1; i <= P.msize(); ++i) {
        if (mD[i] > threshold) {
            P[i] /= mD[i];  // multiplying V by S^-1
        }
        else {
            P[i].set(0.);
        }
    }
    mX.mult (~P, ~Q);
}

template<>
CVM_API void
__pinv<double, basic_srbmatrix<double>, basic_rmatrix<double> >
    (basic_rmatrix<double>& mX, 
    const basic_srbmatrix<double>& mArg, double threshold) throw (cvmexception)
{
    static const int zero(0);
    const int m = mArg.msize();
    int nOutInfo = 0;

    basic_rvector<double>   mD       (m);
    basic_rvector<double>   vOffDiag (_cvm_max<int>(1, m - 1));
    basic_srmatrix<double>  mQ       (m);
    basic_srmatrix<double>  mPT      (m);
    basic_srmatrix<double>  mC       (1);
    basic_rvector<double>   vWork    (2 * m);
    basic_srbmatrix<double> mA       (mArg);

    DGBBRD (Chars::pB(), 
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            mA._pm(), mA._pn(), &zero, mA._pl(), mA._pu(), mA, mA._pld(),
            mD, vOffDiag,
            mQ, mQ._pm(),
            mPT, mPT._pm(),
            mC, mC._pm(), 
            vWork, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);

    basic_rvector<double> vWork2 (m * 4);
    DBDSQR (Chars::pU(),
#if defined (CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES)
            1,
#endif
            &m,
            &m,
            &m,
            &zero,
            mD, vOffDiag,
            mPT, mPT._pm(),
            mQ, mQ._pm(),
            mC, mC._pm(), 
            vWork2, &nOutInfo);
    if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
    if (nOutInfo > 0) throw cvmexception (CVM_CONVERGENCE_ERROR);

    for (int i = 1; i <= m; ++i) {
        if (mD[i] > threshold) {
            mPT[i] /= mD[i];  // multiplying V by S^-1
        }
        else {
            mPT[i].set(0.);
        }
    }
    mX.mult (~mPT, ~mQ);
}

CVM_NAMESPACE_END
