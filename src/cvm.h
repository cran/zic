//                  CVM Class Library
//                  http://cvmlib.com
//
//          Copyright Sergei Nikolaev 1992-2008
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef _CVM_H
#define _CVM_H

#if defined (__INTEL_COMPILER)
#   pragma warning(disable:1744)
#   pragma warning(disable:1125)
#   pragma warning(disable:1195)
#   pragma warning(disable:383)
#   pragma warning(disable:810)
#   pragma warning(disable:981)
#   pragma warning(disable:1418)
#   pragma warning(disable:171)
#   pragma warning(disable:1684)
#   pragma warning(disable:1599)
#endif

#   if defined (_MT) && !defined (CVM_NO_MT)
#       define CVM_MT
#       if !defined (_PTHREADS)
#           define _PTHREADS
#       endif
#   endif

// MSVC++ 6.0 and higher settings
#if defined (_MSC_VER)
#   pragma once
#   define WIN32_LEAN_AND_MEAN        // Exclude rarely-used stuff from Windows headers
#   ifndef _WIN32_WINNT
#       define _WIN32_WINNT 0x500       // at least Win2000 is required for InitializeCriticalSectionAndSpinCount
#   endif
#   include <windows.h>
#   include <process.h>
#   include <time.h>
#   if (_MSC_VER < 1400)
#       error "Please use stable version 5.2 for older MSVC compilers"
#   endif
#   if (!defined(__INTEL_COMPILER) || !defined(_WIN64)) && !defined(CVM_ACML) && !(_MSC_VER >= 1500 && defined(_WIN64))
#       define CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES
#   endif
#   if defined(CVM_ACML)
#       define CVM_COMPLEX_NUMBER_RETURNED
#   endif
#   pragma warning(disable:4290)
#   pragma warning(push)
#   pragma warning(disable:4250)
#   pragma warning(disable:4251)
#   pragma warning(disable:4311)
#   if defined (CVM_FLOAT)
#       pragma warning(disable:4244)
#   endif
#   if defined (SRC_EXPORTS) && !defined (CVM_EXPORTS)
#       define CVM_EXPORTS
#   endif
#   include <hash_map>
#   if defined (CVM_USES_STLPORT)
#       define CVM_BLOCKS_MAP std::hash_map
#   else
#       define CVM_BLOCKS_MAP stdext::hash_map
#   endif
#   ifdef CVM_STATIC
#       define CVM_API
#   else
#       ifdef CVM_EXPORTS
#           define CVM_API __declspec(dllexport)
#       else
#           define CVM_API __declspec(dllimport)
#       endif
#   endif

#   include <limits>
    typedef __int64 CVM_LONGEST_INT;
#   define CVM_VSNPRINTF vsnprintf_s
#   define CVM_VSNPRINTF_S_DEFINED
#   define CVM_STRCPY_S_DEFINED
#endif

// GCC settings
#if defined (__GNUC__)
#   ifdef __MINGW32__               // Dev-C++ under Win32 assumed here
#       define WIN32_LEAN_AND_MEAN
#       include <windows.h>
#       include <process.h>
#   else
#       include <semaphore.h>       // Unix
#   endif
#   ifdef __stdcall
#       undef __stdcall
#   endif
#   define __stdcall
#   define CVM_API
#   define CVM_STDEXT stdext
#   define CVM_BLOCKS_MAP std::map

typedef long long CVM_LONGEST_INT;

#   if defined(__AMD64__)
        typedef unsigned long long CVM_PTR_WRAPPER;
#   elif defined(__MINGW64__)
        typedef unsigned long long CVM_PTR_WRAPPER;
#   else
        typedef unsigned long CVM_PTR_WRAPPER;
#   endif

#   define CVM_VSNPRINTF vsnprintf

#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <limits>
#include <string.h>

// fix for missing __builtin_clog functions
#if defined (__INTEL_COMPILER)
#   if defined (_GLIBCXX_USE_C99_COMPLEX)
#       undef _GLIBCXX_USE_C99_COMPLEX
#       define _GLIBCXX_USE_C99_COMPLEX 0
#   endif
#endif

#include <complex>
#include <algorithm>
#include <exception>
#include <new>

#if defined (STLPORT)
#   define CVM_USES_STLPORT
#endif

#if defined (_DEBUG) || defined (DEBUG)
#   define CVM_DEBUG
#   include <assert.h>
#   define CVM_ASSERT(p,n) _cvm_assert(p,n);
#else
#   define CVM_ASSERT(p,n)
#endif

// errors codes
#define CVM_OK                           0
#define CVM_OUTOFMEMORY                  1
#define CVM_OUTOFRANGE                   2
#define CVM_OUTOFRANGE1                  3
#define CVM_OUTOFRANGE2                  4
#define CVM_WRONGSIZE                    5
#define CVM_SIZESMISMATCH                6
#define CVM_WRONGMKLARG                  7
#define CVM_WRONGMKLARG2                 8
#define CVM_SINGULARMATRIX               9
#define CVM_NOTPOSITIVEDEFINITE          10
#define CVM_WRONGCHOLESKYFACTOR          11
#define CVM_WRONGBUNCHKAUFMANFACTOR      12
#define CVM_NOTPOSITIVEDIAG              13
#define CVM_CONVERGENCE_ERROR            14
#define CVM_DIVISIONBYZERO               15
#define CVM_SEMAPHOREERROR               16
#define CVM_READ_ONLY_ACCESS             17
#define CVM_SUBMATRIXACCESSERROR         18
#define CVM_SUBMATRIXNOTAVAILABLE        19
#define CVM_MATRIXNOTSYMMETRIC           20
#define CVM_MATRIXNOTHERMITIAN           21
#define CVM_BREAKS_HERMITIANITY          22
#define CVM_METHODNOTAVAILABLE           23
#define CVM_NOTIMPLEMENTED               24

#define CVM_THE_LAST_ERROR_CODE          25                    // use it in derived classes
#define CVM_MATRIX_ELEMENT_SEPARATOR     " "
#define CVM_EOL                          std::endl

typedef unsigned char tbyte;                                   // memory allocation quantum

#ifdef CVM_NO_NAMESPACE
#   define CVM_NAMESPACE_BEG
#   define CVM_NAMESPACE_END
#else
#   define CVM_NAMESPACE_BEG namespace cvm {
#   define CVM_NAMESPACE_END }
#endif


CVM_NAMESPACE_BEG

template <typename T>               class basic_array;
template <typename TR, typename TC> class Array;
template <typename TR, typename TC> class Matrix;
template <typename TR, typename TC> class SqMatrix;
template <typename TR, typename TC> class BandMatrix;
template <typename TR>              class basic_rvector;
template <typename TR>              class basic_rmatrix;
template <typename TR>              class basic_srmatrix;
template <typename TR>              class basic_srbmatrix;
template <typename TR>              class basic_srsmatrix;
template <typename T,  typename TR> class type_proxy;

// error messages holder
class ErrMessages
{
    // message string maps
    typedef std::map <int, std::string, std::less<int> > map_Msg;
    typedef map_Msg::iterator                            itr_Msg;
    typedef map_Msg::const_iterator                      citr_Msg;
    typedef std::pair <int, std::string>                 pair_Msg;

private:
    std::string msUnknown;
    map_Msg mmMsg;
    CVM_API ErrMessages();

public:
    CVM_API const std::string& _get (int nException);
    CVM_API bool _add (int nNewCause, const char* szNewMessage);

    static CVM_API ErrMessages& ErrMessagesInstance();
    ~ErrMessages() {}
};

class cvmexception : public std::exception
{
protected:
    int mnCause;
    mutable char mszMsg[256];

    virtual const char* _get_message(int nCause) const
    {
        return ErrMessages::ErrMessagesInstance()._get(nCause).c_str();
    }

public:
    cvmexception()
        : mnCause (CVM_OK)
    {
        mszMsg[0] = '\0';
    }

    CVM_API explicit cvmexception (int nCause, ...);
    CVM_API cvmexception (const cvmexception& e);

    virtual ~cvmexception() throw() {}

    int cause() const
    {
        return mnCause;
    }

    virtual const char* what() const throw()
    {
        return mszMsg;
    }

    static int getNextCause ()
    {
        return CVM_THE_LAST_ERROR_CODE;
    }

    static bool add (int nNewCause, const char* szNewMessage)
    {
        return ErrMessages::ErrMessagesInstance()._add (nNewCause, szNewMessage);
    }
};

// exported utilities - were recently refactored from specialization because of Borland's understanding of ANSI C++
template <typename T>
CVM_API void __copy (int nSize, const T* pFrom, int nFromIncr, T* pTo, int nToIncr);
template <typename T>
CVM_API void __swap (int nSize, T* p1, int n1Incr, T* p2, int n2Incr);

template <typename TC, typename TR>
CVM_API TR _real (const TC& mT);
template <typename TC, typename TR>
CVM_API TR _imag (const TC& mT);

template<typename TR>
CVM_API TR __dot (const TR* mpD, int mnSize, int mnIncr, const TR* pD, int nIncr);
template<typename TC>
CVM_API TC __dotu (const TC* mpD, int mnSize, int mnIncr, const TC* pD, int nIncr);
template<typename TC>
CVM_API TC __dotc (const TC* mpD, int mnSize, int mnIncr, const TC* pD, int nIncr);

template <typename TR, typename TC>
CVM_API TR __norm (const TC*  pD, int nSize, int nIncr);

template <typename TC>
CVM_API int __idamax (const TC* pD, int nSize, int nIncr);
template <typename TC>
CVM_API int __idamin (const TC* pD, int nSize, int nIncr);

template <typename TC>
CVM_API void __add (TC* mpD, int mnSize, int mnIncr, const TC* pv, int nIncr);

template <typename TC>
CVM_API void __subtract (TC* mpD, int mnSize, int mnIncr, const TC* pv, int nIncr);

template <typename TR, typename TC>
CVM_API void __scal (TC* mpD, int mnSize, int mnIncr, TR dScal);

template<typename TC>
CVM_API void __conj (TC* mpD, int mnSize, int mnIncr);

template <typename TR, typename TC>
CVM_API void __copy2 (TC* mpD, int mnSize, int mnIncr, const TR* pRe, const TR* pIm, int nReIncr = 1, int nImIncr = 1);

template <typename TC, typename TM, typename TV>
CVM_API void __gemv (bool bLeft, const TM& m, TC dAlpha, const TV& v, TC dBeta, TV& vRes);
template <typename TC, typename TM, typename TV>
CVM_API void __gbmv (bool bLeft, const TM& m, TC dAlpha, const TV& v, TC dBeta, TV& vRes);
template <typename TR, typename TM, typename TV>
CVM_API void __symv (const TM& m, TR dAlpha, const TV& v, TR dBeta, TV& vRes);
template <typename TC, typename TM, typename TV>
CVM_API void __shmv (const TM& m, TC dAlpha, const TV& v, TC dBeta, TV& vRes);
template <typename TC, typename TM>
CVM_API void __gemm (const TM& ml, bool bTrans1, const TM& mr, bool bTrans2, TC dAlpha, TM& mRes, TC dBeta);
template <typename TR, typename TSM, typename TM>
CVM_API void __symm (bool bLeft, const TSM& ml, const TM& mr, TR dAlpha, TM& mRes, TR dBeta);
template <typename TC, typename TSM, typename TM>
CVM_API void __hemm (bool bLeft, const TSM& ml, const TM& mr, TC dAlpha, TM& mRes, TC dBeta);

template <typename TC, typename TV>
CVM_API void __polynom (TC* mpD, int ldP, int mnM, const TC* pD, int ldA, const TV& v);
template <typename T>
CVM_API void __inv (T& m, const T& mArg) throw (cvmexception);
template <typename T, typename TR>
CVM_API void __exp (T& m, const T& mArg, TR tol) throw (cvmexception);
template <typename T, typename TR>
CVM_API void __exp_symm (T& m, const T& mArg, TR tol) throw (cvmexception);
template <typename TR, typename TC, typename TRM>
CVM_API void __solve (const TRM& m, int nrhs, const TC* pB, int ldB, TC* pX, int ldX, TR& dErr, const TC* pLU, const int* pPivots) throw (cvmexception);
template <typename TC, typename TM, typename TSM>
CVM_API void __svd (TC* pD, int nSize, int nIncr, const TM& mArg, TSM* mU, TSM* mVH) throw (cvmexception);
template <typename TR, typename TM, typename TX>
CVM_API void __pinv (TX& mX, const TM& mArg, TR threshold) throw (cvmexception);
template <typename TV, typename TSM, typename TSCM>
CVM_API void __eig (TV& vRes, const TSM& mArg, TSCM* mEigVect, bool bRightVect) throw (cvmexception);
template<typename TR, typename TM>
CVM_API void __cond_num (const TM& mArg, TR& dCondNum) throw (cvmexception);
template <typename TM>
void __low_up (TM& m, int* nPivots) throw (cvmexception);

template <typename TR, typename TM, typename TV>
CVM_API void __ger (TM& m, const TV& vCol, const TV& vRow, TR dAlpha);
template <typename TC, typename TM, typename TV>
CVM_API void __geru (TM& m, const TV& vCol, const TV& vRow, TC cAlpha);
template <typename TC, typename TM, typename TV>
CVM_API void __gerc (TM& m, const TV& vCol, const TV& vRow, TC cAlpha);

template <typename TC, typename TSM>
CVM_API void __syrk (bool bTransp, TC alpha, int k, const TC* pA, int ldA, TC beta, TSM& m);
template <typename TC, typename TSM>
CVM_API void __syr2k (bool bTransp, TC alpha, int k, const TC* pA, int ldA, const TC* pB, int ldB, TC beta, TSM& m);
template <typename TC, typename TSM>
CVM_API void __herk (bool bTransp, TC alpha, int k, const TC* pA, int ldA, TC beta, TSM& m);
template <typename TC, typename TSM>
CVM_API void __her2k (bool bTransp, TC alpha, int k, const TC* pA, int ldA, const TC* pB, int ldB, TC beta, TSM& m);

template <typename TM>
CVM_API int __cholesky (TM& m);
template <typename TM>
CVM_API void __bunch_kaufman (TM& m, int* nPivots) throw (cvmexception);
template<typename TR, typename TM, typename TV>
CVM_API void __poequ (const TM& m, TV& vScalings, TR& dCond, TR& dMax);

template <typename TM, typename TSM>
CVM_API void __qre (const TM& mA, TM& mQ, TSM& mR) throw (cvmexception);
template <typename TM, typename TSM>
CVM_API void __qrf (const TM& mA, TSM& mQ, TM& mR) throw (cvmexception);


template <typename T> inline
const T& _cvm_min (const T& x, const T& y)
{
    return x < y ? x : y;
}

template <typename T> inline
const T& _cvm_max (const T& x, const T& y)
{
    return x > y ? x : y;
}

template<class TR>
inline const TR* __get_real_p (const std::complex<TR>* c)
{
#if defined (CVM_USES_STLPORT)
    return &c->_M_re;
#else
    return reinterpret_cast<const TR*>(c);
#endif
}

template<class TR>
inline const TR* __get_imag_p (const std::complex<TR>* c)
{
#if defined (CVM_USES_STLPORT)
    return &c->_M_im;
#else
    return reinterpret_cast<const TR*>(c) + 1;
#endif
}

template<class TR>
inline TR* __get_real_p (std::complex<TR>* c)
{
#if defined (CVM_USES_STLPORT)
    return &c->_M_re;
#else
    return reinterpret_cast<TR*>(c);
#endif
}

template<class TR>
inline TR* __get_imag_p (std::complex<TR>* c)
{
#if defined (CVM_USES_STLPORT)
    return &c->_M_im;
#else
    return reinterpret_cast<TR*>(c) + 1;
#endif
}

template<class TR> inline const TR* __get_real_p (const TR* d) { return d; }
template<class TR> inline const TR* __get_imag_p (const TR* d) { return d; }
template<class TR> inline       TR* __get_real_p (      TR* d) { return d; }
template<class TR> inline       TR* __get_imag_p (      TR* d) { return d; }

// characters for fortran subroutines
class Chars
{
protected:
    static char mchars[15];

public:
    static const char* pT () {return mchars;}
    static const char* pN () {return mchars + 1;}
    static const char* pU () {return mchars + 2;}
    static const char* pL () {return mchars + 3;}
    static const char* pP () {return mchars + 4;}
    static const char* pQ () {return mchars + 5;}
    static const char* pB () {return mchars + 6;}
    static const char* pE () {return mchars + 7;}
    static const char* pR () {return mchars + 8;}
    static const char* pA () {return mchars + 9;}
    static const char* pS () {return mchars + 10;}
    static const char* pV () {return mchars + 11;}
    static const char* pO () {return mchars + 12;}
    static const char* pI () {return mchars + 13;}
    static const char* pC () {return mchars + 14;}
};

template <class TR>
inline const TR& basic_cvmMachMin()
{
    static const TR _min = (std::numeric_limits<TR>::min)();
    return _min;
}

template <class TR>
inline const TR& basic_cvmMachSp()
{
    static const TR _eps = (std::numeric_limits<TR>::epsilon)();
    return _eps;
}


// class of memory blocks
class CVM_API MemoryBlocks
{
    struct BlockProperty
    {
        size_t mnSize;
        int mnRefCount;
        BlockProperty (size_t nSize, int nRefCount) : mnSize (nSize), mnRefCount (nRefCount) {}
    };

#ifdef CVM_USE_POOL_MANAGER
    typedef std::map<tbyte*, BlockProperty, std::less<tbyte*> > map_Blocks;     // pointer -> {size, refcount}
    typedef map_Blocks::iterator itr_Blocks;

    typedef std::multimap<size_t, tbyte*> map_FreeBs;                           // size -> pointer
    typedef map_FreeBs::iterator itr_FreeBs;

    typedef std::map<tbyte*, itr_FreeBs, std::less<tbyte*> > map_FreeIt;        // pointer -> iterator to FreeBs
    typedef map_FreeIt::iterator itr_FreeIt;

    map_FreeBs mFreeBs;                                                         // currently free blocks by sizes
    map_FreeIt mFreeIt;                                                         // currently free blocks iterators by pointers
#else
    typedef CVM_BLOCKS_MAP<CVM_PTR_WRAPPER, BlockProperty> map_Blocks;          // pointer -> refcount
    typedef map_Blocks::iterator itr_Blocks;
#endif

    map_Blocks mBlocks;                                                         // currently occupied or freed blocks by pointers

public:
#ifdef CVM_USE_POOL_MANAGER
    void    AddBlock     (tbyte* pBlock, size_t nBytes, bool bOccupied);
    tbyte*  GetFreeBlock (size_t nBytes);
    void    AddPair      (tbyte* pBlock, size_t nBytes, size_t nRest);
#   ifdef CVM_DEBUG
    void    Assert       (const void* pvBlock, size_t nBytes);
#   endif
#else
    void    AddNew       (tbyte* pBlock, size_t nBytes);                        // for just allocated only, i.e. non-const
#endif

    tbyte*  AddRef       (const tbyte* pBlock);
    int     FreeBlock    (tbyte* pBlock);
};


// memory pool class
class CVM_API MemoryPool
{
#ifdef CVM_USE_POOL_MANAGER
    typedef std::list<tbyte*> list_blocks;

    struct DeletePtr {
        template<class T>
        void operator () (T* p) const
        {
            ::delete[] p;
        }
    };

    list_blocks  mOutBlocks;                                                    // outer memory blocks
#endif

    MemoryBlocks mMemoryBlocks;                                                 // currently existing blocks and their statuses

public:
    MemoryPool();
   ~MemoryPool();

    tbyte* Malloc (size_t nBytes) throw (cvmexception);
    tbyte* AddRef (const tbyte* pD);                                            // increases a reference counter
    int    FFFree   (tbyte*& pToFree) throw (cvmexception);                       // decreases a reference counter and
                                                                                // returns memory back to the pool if the counter is zeroed
#ifdef CVM_USE_POOL_MANAGER
#   ifdef CVM_DEBUG
    void Assert (const void* pvBlock, size_t nBytes)                            // synchronized outside
    {
        mMemoryBlocks.Assert (pvBlock, nBytes);
    }
#   endif
    void Clear();                                                               // destroys all outer blocks in reverse order
#endif
};


CVM_API tbyte* _cvmMalloc  (size_t nBytes) throw (cvmexception);
CVM_API tbyte* _cvmAddRef  (const tbyte* pD);
CVM_API int    _cvmFree    (tbyte*& pD);
CVM_API void   _cvm_assert (const void* pvBlock, size_t nBytes);
CVM_API void   cvmExit();

template <typename T>
inline T* cvmMalloc (size_t nEls) throw (cvmexception)
{
    return reinterpret_cast<T*>(_cvmMalloc (nEls * sizeof (T)));
}
template <typename T>
inline T* cvmAddRef (const T* pD)
{
    return reinterpret_cast<T*>(_cvmAddRef (reinterpret_cast<const tbyte*>(pD)));
}
template <typename T>
inline int cvmFree (T*& pD)
{
    return _cvmFree (reinterpret_cast<tbyte*&>(pD));
}



// inline utilities
template <typename T>
inline void CleanMemory (T* p, int nEls)
{
    memset (p, 0, nEls * sizeof(T));
}

inline float _abs (const float& v)
{
    return static_cast<float>(fabs (v));
}
inline double _abs (const double& v)
{
    return fabs(v);
}
inline long double _abs (const long double& v)
{
    return fabs(v);
}
inline float _abs (const std::complex<float>& v)
{
    return static_cast<float>(sqrt (v.real() * v.real() + v.imag() * v.imag()));
}
inline double _abs (const std::complex<double>& v)
{
    return sqrt (v.real() * v.real() + v.imag() * v.imag());
}
inline long double _abs (const std::complex<long double>& v)
{
    return sqrt (v.real() * v.real() + v.imag() * v.imag());
}
inline int _abs (const int& v)
{
    return abs(v);
}
inline float _sqrt (const float& v)
{
    return static_cast<float>(sqrt (v));
}
inline double _sqrt (const double& v)
{
    return sqrt(v);
}
inline std::complex<float> _conjugate (const std::complex<float>& v)
{
    return std::complex<float> (v.real(), - v.imag());
}
inline std::complex<double> _conjugate (const std::complex<double>& v)
{
    return std::complex<double> (v.real(), - v.imag());
}
template <typename TR>
inline bool _conjugated (const std::complex<TR>& v1, const std::complex<TR>& v2, TR tol)
{
    return _abs (v1.real() - v2.real()) <= tol &&
           _abs (v1.imag() + v2.imag()) <= tol;
}

template <typename T, typename TR>
inline std::ostream& operator << (std::ostream& os, const type_proxy<T,TR>& mOut)
{
#if defined (_MSC_VER) && !defined (CVM_USES_STLPORT)
    os.imbue (std::locale::empty());
#endif
    os << static_cast<T>(mOut);
    return os;
}

template <typename T, typename TR>
inline std::istream& operator >> (std::istream& is, type_proxy<T,TR>& v)
{
#if defined (_MSC_VER) && !defined (CVM_USES_STLPORT)
    is.imbue (std::locale::empty());
#endif
    T t;
    is >> t;
    v = t;
    return is;
}

template <typename T>
std::istream& operator >> (std::istream& is, basic_array<T>& aIn)
{
#if !defined (CVM_USES_STLPORT) && defined (_MSC_VER)
    is.imbue (std::locale::empty());
#endif
    CVM_ASSERT(aIn.mpD, aIn.mnSize * sizeof(T))
    for (int i = 0; i < aIn.mnSize && is.good(); ++i)
    {
        is >> aIn.mpD[i];
    }
    return is;
}

template <typename T>
std::ostream& operator << (std::ostream& os, const basic_array<T>& aOut)
{
#if !defined (CVM_USES_STLPORT) && defined (_MSC_VER)
    os.imbue (std::locale::empty());
#endif
    CVM_ASSERT(aOut.mpD, aOut.mnSize * sizeof(T))
    for (int i = 0; i < aOut.mnSize && os.good(); ++i)
    {
        os << aOut.mpD[i] << CVM_MATRIX_ELEMENT_SEPARATOR;
    }
    os << CVM_EOL;
    return os;
}

template <typename TR, typename TC>
std::istream& operator >> (std::istream& is, Array<TR,TC>& aIn)
{
#if !defined (CVM_USES_STLPORT) && defined (_MSC_VER)
    is.imbue (std::locale::empty());
#endif
    const int nSize = aIn.mnSize * aIn.mnIncr;
    CVM_ASSERT(aIn.mpD, ((aIn.mnSize - 1) * aIn.mnIncr + 1) * sizeof(TC))
    for (int i = 0; i < nSize && is.good(); i += aIn.mnIncr)
    {
        is >> aIn.mpD[i];
    }
    return is;
}

template <typename TR, typename TC>
std::ostream& operator << (std::ostream& os, const Array<TR,TC>& aOut)
{
#if !defined (CVM_USES_STLPORT) && defined (_MSC_VER)
    os.imbue (std::locale::empty());
#endif
    const int nSize = aOut.mnSize * aOut.mnIncr;
    CVM_ASSERT(aOut.mpD, ((aOut.mnSize - 1) * aOut.mnIncr + 1) * sizeof(TC))

    for (int i = 0; i < nSize && os.good(); i += aOut.mnIncr)
    {
        os << aOut.mpD[i] << CVM_MATRIX_ELEMENT_SEPARATOR;
    }
    os << CVM_EOL;
    return os;
}

template <typename TR, typename TC>
std::ostream& operator << (std::ostream& os, const Matrix<TR,TC>& mOut)
{
#if defined (_MSC_VER) && !defined (CVM_USES_STLPORT)
    os.imbue (std::locale::empty());
#endif
    for (int i = 0; i < mOut.mnM; ++i)
    {
        for (int j = 0; j < mOut.mnN && os.good(); ++j)
        {
            os << mOut._ij_val (i, j) << CVM_MATRIX_ELEMENT_SEPARATOR;
        }
        os << CVM_EOL;
    }
    return os;
}

template <typename TR, typename TC>
std::istream& operator >> (std::istream& is, Matrix<TR,TC>& mIn)
{
#if defined (_MSC_VER) && !defined (CVM_USES_STLPORT)
    is.imbue (std::locale::empty());
#endif
    TC v;
    for (int i = 0; i < mIn.mnM; ++i)
    {
        for (int j = 0; j < mIn.mnN && is.good(); ++j)
        {
            is >> v;
            mIn._ij_proxy_val (i, j) = v;
        }
    }
    return is;
}

template <typename TR, typename TC, typename RM, typename RBM>
inline void _copy_b_matrix (RM& m, RBM& mb, bool bLeftToRight)
{
    const int nM   = mb.msize();
    const int nN   = mb.nsize();
    const int nKL  = mb.lsize();
    const int nKU  = mb.usize();
    const int nCol = 1 + nKL + nKU;
    int nS, nShiftL, nShiftR;
    TC* pL;
    TC* pR;

    for (int i = 0; i < nN; ++i)
    {
        nS = nCol;
        nShiftL = 0;
        nShiftR = 0;
        if (i < nKU)
        {
            nShiftR = nKU - i;
            nS -= nShiftR;
        }
        else
        {
            nShiftL = i - nKU;
        }
        if (nN - i <= nKL)
        {
            nS -= nKL + 1 - (nN - i);
        }

        pL = m.get()  + i * nM + nShiftL;
        pR = mb.get() + i * nCol + nShiftR;

        __copy<TC> (nS,
                    bLeftToRight ? pL : pR,
                    1,
                    bLeftToRight ? pR : pL,
                    1);
    }
}


template <typename TR, typename TC>
inline void _sum (TC* pD, int nSize, int nIncr,
                  const TC* p1, int nIncr1, const TC* p2, int nIncr2)           // pD = a1 + a2
{
    if (pD == p1)
    {
        if (pD == p2)
        {
            static const TR two(2.);
            __scal<TR, TC> (pD, nSize, nIncr, two);
        }
        else
        {
            __add<TC> (pD, nSize, nIncr, p2, nIncr2);
        }
    }
    else
    {
        if (pD == p2)
        {
            __add<TC> (pD, nSize, nIncr, p1, nIncr1);
        }
        else
        {
            __copy<TC> (nSize, p1, nIncr1, pD, nIncr);
            __add<TC> (pD, nSize, nIncr, p2, nIncr2);
        }
    }
}

template <typename TR, typename TC>
inline void _diff (TC* pD, int nSize, int nIncr, 
                   const TC* p1, int nIncr1, const TC* p2, int nIncr2)          // this = v1 - v2
{
    if (pD == p1)
    {
        if (pD == p2)
        {
            static const TR zero(0.);
            __scal<TR, TC> (pD, nSize, nIncr, zero);
        }
        else
        {
            __subtract<TC> (pD, nSize, nIncr, p2, nIncr2);
        }
    }
    else
    {
        if (pD == p2)
        {
            static const TR mone(-1.);
            __subtract<TC> (pD, nSize, nIncr, p1, nIncr1);
            __scal<TR, TC> (pD, nSize, nIncr, mone);
        }
        else
        {
            __copy<TC> (nSize, p1, nIncr1, pD, nIncr);
            __subtract<TC> (pD, nSize, nIncr, p2, nIncr2);
        }
    }
}

template <typename TR, typename TC>
inline void _incr (TC* mpD, int mnSize, int mnIncr, const TC* pD, int nIncr)    // mpD = mpD + pD
{
    if (mpD == pD)
    {
        static const TR two(2.);
        __scal<TR, TC> (mpD, mnSize, mnIncr, two);
    }
    else
    {
        __add<TC> (mpD, mnSize, mnIncr, pD, nIncr);
    }
}

template <typename TR, typename TC>
inline void _decr (TC* mpD, int mnSize, int mnIncr, const TC* pD, int nIncr)    // mpD = mpD - pD
{
    if (mpD == pD)
    {
        static const TR zero(0.);
        __scal<TR, TC> (mpD, mnSize, mnIncr, zero);
    }
    else
    {
        __subtract<TC> (mpD, mnSize, mnIncr, pD, nIncr);
    }
}

template <typename TR, typename TC>
void _set_real (TC* mpD, int mnSize, int mnIncr, TR d)                          // fills real part
{
    const int nIncr2 = 2 * mnIncr;
    const int nSize = mnSize * nIncr2;
    TR* pD = __get_real_p<TR>(mpD);

    for (int i = 0; i < nSize; i += nIncr2)
    {
        CVM_ASSERT(pD, (i + 1) * sizeof(TR))
        pD[i] = d;
    }
}

template <typename TR, typename TC>
void _set_imag (TC* mpD, int mnSize, int mnIncr, TR d)                          // fills imaginary part
{
    const int nIncr2 = 2 * mnIncr;
    const int nSize = mnSize * nIncr2;
    TR* pD = __get_imag_p<TR>(mpD);

    for (int i = 0; i < nSize; i += nIncr2)
    {
        CVM_ASSERT(pD, (i + 1) * sizeof(TR))
        pD[i] = d;
    }
}


// this class provides read-write differentiation for a particular type
template <typename T, typename TR>
class type_proxy
{
    typedef type_proxy<T,TR> P;
    
protected:
    T&   mT;
    bool mbReadOnly;

public:
    type_proxy (T& ref, bool read_only) : mT(ref), mbReadOnly(read_only)
    {
    }

    type_proxy (const T& ref, bool read_only = true) :
            mT(const_cast<T&>(ref)), mbReadOnly(read_only)                      // read only by definition
    {
    }

    type_proxy (const type_proxy& p) : mT(p.mT), mbReadOnly(p.mbReadOnly)
    {
    }

    operator T () const
    {
        return mT;
    }
/*
    doesn't work - masks the operator above
    operator const T& () const
    {
        return mT;
    }

    operator T& ()
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        return mT;
    }
*/
    T val() const
    {
        return mT;
    }

    T& get() throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        return mT;
    }
    const T& get() const
    {
        return mT;
    }

    type_proxy& operator = (const type_proxy& p)
    {
        mT = p.mT;
        mbReadOnly = p.mbReadOnly;
        return *this;
    }
    type_proxy& operator = (const T& v) throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        mT = v;
        return *this;
    }

    T* operator & () throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        return &mT;
    }
    const T* operator & () const
    {
        return &mT;
    }

    template <typename U>
    T& operator += (const U& u) throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        mT += T(u);
        return mT;
    }
    template <typename U>
    T& operator -= (const U& u) throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        mT -= T(u);
        return mT;
    }
    template <typename U>
    T& operator *= (const U& u) throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        mT *= T(u);
        return mT;
    }
    template <typename U>
    T& operator /= (const U& u) throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        mT /= T(u);
        return mT;
    }

    template <typename U>
    T operator + (const U& u) const
    {
        return mT + T(u);
    }
    template <typename U>
    T operator - (const U& u) const
    {
        return mT - T(u);
    }
    template <typename U>
    T operator * (const U& u) const
    {
        return mT * T(u);
    }
    template <typename U>
    T operator / (const U& u) const
    {
        return mT / T(u);
    }

    T operator + (const P& u) const
    {
        return mT + u.mT;
    }
    T operator - (const P& u) const
    {
        return mT - u.mT;
    }
    T operator * (const P& u) const
    {
        return mT * u.mT;
    }
    T operator / (const P& u) const
    {
        return mT / u.mT;
    }

    T operator - () const
    {
        return - mT;
    }

    // specialized for std::complex<treal> only. link error would be received in other case
    TR real() const
    {
        return _real<T,TR>(mT);
    }
    TR imag() const
    {
        return _imag<T,TR>(mT);
    }
};


template <typename T, typename TR>
inline std::complex<TR> operator + (std::complex<TR> u, const type_proxy<T,TR>& p)
{
    return u + p.val();
}
template <typename T, typename TR>
inline std::complex<TR> operator - (std::complex<TR> u, const type_proxy<T,TR>& p)
{
    return u - p.val();
}
template <typename T, typename TR>
inline std::complex<TR> operator * (std::complex<TR> u, const type_proxy<T,TR>& p)
{
    return u * p.val();
}
template <typename T, typename TR>
inline std::complex<TR> operator / (std::complex<TR> u, const type_proxy<T,TR>& p)
{
    return u / p.val();
}

template <typename U, typename T, typename TR>
inline T operator + (U u, const type_proxy<T,TR>& p)
{
    return T(u) + p.val();
}
template <typename U, typename T, typename TR>
inline T operator - (U u, const type_proxy<T,TR>& p)
{
    return T(u) - p.val();
}
template <typename U, typename T, typename TR>
inline T operator * (U u, const type_proxy<T,TR>& p)
{
    return T(u) * p.val();
}
template <typename U, typename T, typename TR>
inline T operator / (U u, const type_proxy<T,TR>& p)
{
    return T(u) / p.val();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// basic_array
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
class basic_array
{
protected:
    int mnSize;                                                                 // number of array elements allocated
    T*  mpD;                                                                    // data pointer

public:
    // STL type definitions
    typedef T value_type;
    typedef value_type* pointer;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;

    basic_array()
        : mnSize(0), mpD(NULL)
    {
    }

    explicit basic_array (int nSize, bool bZeroMemory = true)
        : mnSize(nSize), mpD (cvmMalloc<T>(size_t(mnSize)))
    {
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        if (bZeroMemory) CleanMemory<T> (mpD, mnSize);
    }

    basic_array (const T* p, int nSize)
        : mnSize(nSize), mpD (cvmMalloc<T>(size_t(mnSize)))
    {
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        __copy<T> (mnSize, p, 1, mpD, 1);
    }

    basic_array (const T* first, const T* last)
        : mnSize(int(last - first)), mpD (cvmMalloc<T>(size_t(mnSize)))
    {
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        __copy<T> (mnSize, first, 1, mpD, 1);
    }

    basic_array (const basic_array& a)
        : mnSize(a.mnSize), mpD (cvmMalloc<T>(size_t(mnSize)))
    {
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        __copy<T> (mnSize, a.mpD, 1, mpD, 1);
    }

    virtual ~basic_array()
    {
        cvmFree<T>(mpD);
    }

    int size() const
    {
        return mnSize;
    }

    T* get()
    {
        return mpD;
    }

    const T* get() const
    {
        return mpD;
    }

    operator T* ()
    {
        return mpD;
    }

    operator const T* () const
    {
        return mpD;
    }

    // 1-based indexing operators
    T& operator () (int nFI) throw (cvmexception)                               // element access (returns l-value)
    {
        return this -> at(size_type(nFI - 1));
    }

    T operator () (int nFI) const throw (cvmexception)                          // element access (does not return l-value)
    {
        return this -> at(size_type(nFI - 1));
    }

    T& operator [] (size_type n) throw (cvmexception)                           // element access (returns l-value)
    {
        return this -> at(n - 1);
    }
    T operator [] (size_type n) const throw (cvmexception)                      // element access (does not return l-value)
    {
        return this -> at(n - 1);
    }
    T& operator [] (int n) throw (cvmexception)                                 // element access (returns l-value)
    {
        return this -> at(size_type(n - 1));
    }
    T operator [] (int n) const throw (cvmexception)                            // element access (does not return l-value)
    {
        return this -> at(size_type(n - 1));
    }

    basic_array& operator = (const basic_array& a) throw (cvmexception)
    {
        if (mnSize != a.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign (a.mpD, 1);
        return *this;
    }

    basic_array& assign (const T* p)
    {
        this -> _assign (p, 1);
        return *this;
    }

    basic_array& set (T x)
    {
        this -> _set (x);                                                       // fills the content
        return *this;
    }

    basic_array& resize (int nNewSize) throw (cvmexception)
    {
        this -> _resize (nNewSize);
        return *this;
    }

    // std stuff
    iterator       begin()        {return this -> mpD;}
    const_iterator begin() const  {return this -> mpD;}
    iterator       end()          {return this -> mpD + mnSize;}
    const_iterator end()   const  {return this -> mpD + mnSize;}

    reverse_iterator rbegin()             {return reverse_iterator(end());}
    const_reverse_iterator rbegin() const {return const_reverse_iterator(end());}
    reverse_iterator rend()               {return reverse_iterator(begin());}
    const_reverse_iterator rend()   const {return const_reverse_iterator(begin());}

    size_type max_size() const    {return size_type(-1) / sizeof(T);}
    size_type capacity() const    {return size_type(mnSize);}
    bool empty() const            {return mnSize > 0;}

    reference front()             {return *begin();}
    const_reference front() const {return *begin();}
    reference back()              {return *(end() - 1);}
    const_reference back()  const {return *(end() - 1);}

//    deprecated since 5.4.2 due to its incosistency with STL counterpart
//
//    void reserve (size_type n) throw (cvmexception)
//    {
//        this -> _resize (int(n));
//    }

    void assign (size_type n, const T& val) throw (cvmexception)
    {
        if (n > mnSize) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        for (int i = 0; i < n; ++i)
        {
            mpD[i] = val;
        }
    }

    void assign (const_iterator first, const_iterator last) throw (cvmexception)
    {
        const int n = last - first;
        if (n > mnSize) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        for (int i = 0; i < n; ++i)
        {
            mpD[i] = *(first + i);
        }
    }

    void resize (size_type nNewSize) throw (cvmexception)
    {
        this -> _resize (int(nNewSize));
    }

    void clear()
    {
        this -> _resize (0);
    }

    void swap (basic_array& v) throw (cvmexception)
    {
        if (mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        __swap<T> (mnSize, mpD, 1, v.mpD, 1);
    }

    // 0-based
    reference at (size_type n) throw (cvmexception)
    {
        return this -> _at(n);
    }

    const_reference at (size_type n) const throw (cvmexception)
    {
        return this -> _at(n);
    }

    // very very slow. provided for compatibility only
    void push_back (const T& x) throw (cvmexception)
    {
        this -> _resize (mnSize + 1);
        mpD[mnSize - 1] = x;
    }

    // very very slow. provided for compatibility only
    void pop_back () throw (cvmexception)
    {
        if (mnSize > 0) this -> _resize (mnSize - 1);
    }

    // very very slow. provided for compatibility only
    iterator insert (iterator position, const T& x) throw (cvmexception)
    {
        const int n = int (position - this -> begin());
        if (n > mnSize || n < 0) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        this -> _resize (mnSize + 1);
        for (int i = mnSize - 1; i > n; --i)
        {
            mpD[i] = mpD[i-1];
        }
        mpD[n] = x;
        return iterator (mpD + n);
    }

    // very very slow. provided for compatibility only
    iterator erase (iterator position) throw (cvmexception)
    {
        const int n = int (position - this -> begin());
        if (n > mnSize || n < 0) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        for (int i = n; i < mnSize - 1; ++i)
        {
            mpD[i] = mpD[i+1];
        }
        if (mnSize > 0) this -> _resize (mnSize - 1);
        return iterator (mpD + n);
    }
    // end of stl stuff

protected:
    // 0-based
    virtual T& _at (size_type n) throw (cvmexception)
    {
        if (n >= size_type(mnSize)) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, (n + 1) * sizeof(T))
        return mpD [n];
    }

    // 0-based
    virtual const T& _at (size_type n) const throw (cvmexception)
    {
        if (n >= size_type(mnSize)) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, (n + 1) * sizeof(T))
        return mpD [n];
    }

    virtual void _assign (const T* p, int)
    {
        if (mpD != p)
        {
            memcpy (mpD, p, mnSize * sizeof(T));
        }
    }

    virtual void _set (T d)                                                     // fills the content
    {
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        for (int i = 0; i < mnSize; ++i)
        {
            mpD[i] = d;
        }
    }

    virtual void _resize (int nNewSize) throw (cvmexception)
    {
        if (nNewSize < 0) throw cvmexception (CVM_WRONGSIZE, nNewSize);
        if (nNewSize == 0)                                                      // just let's free object memory in this case
        {
            cvmFree<T>(mpD);
            mnSize = 0;
        }
        else if (nNewSize != mnSize)
        {
            T* pD = cvmMalloc<T>(nNewSize);
            if (nNewSize > mnSize) CleanMemory<T> (pD, nNewSize);
            const int nMinSize = _cvm_min<int>(nNewSize, mnSize);

            if (nMinSize > 0)
            {
                __copy<T> (nMinSize, mpD, 1, pD, 1);
            }
            cvmFree<T>(mpD);
            mpD = pD;
            mnSize = nNewSize;
            CVM_ASSERT(mpD, mnSize * sizeof(T))
        }
    }

    friend std::istream& operator >> <> (std::istream& is, basic_array& aIn);
    friend std::ostream& operator << <> (std::ostream& os, const basic_array& aOut);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// basic_array
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// abstract array of numbers allocatable in pool
template <typename TR, typename TC>
class Array : public basic_array<TC>
{
    typedef size_t size_type;
    typedef basic_array<TC> BasicArray;

protected:
    int mnIncr;                                                                 // distance between array members (default is 1)

    Array()
        : mnIncr (0)
    {
    }

    explicit Array (int nSize, bool bZeroMemory = true)
        : BasicArray(nSize, bZeroMemory), mnIncr (1)
    {
    }

    Array (TC* pD, int nSize, int nIncr)
        : mnIncr (nIncr)
    {
        if (nSize <= 0) throw cvmexception (CVM_WRONGSIZE, nSize);
        this -> mnSize = nSize;
        this -> mpD = cvmAddRef<TC> (pD);
        CVM_ASSERT(this -> mpD, ((this -> mnSize - 1) * this -> mnIncr + 1) * sizeof(TC))
    }

public:
    int incr() const
    {
        return mnIncr;
    }

    int indofmax() const                                                        // index of max. module, undefined for matrices
    {
        return this -> _indofmax();
    }

    int indofmin() const                                                        // index of max. module, undefined for matrices
    {
        return this -> _indofmin();
    }

    virtual TR norm() const                                                     // Euclid norm
    {
        return __norm<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr);
    }

    virtual TR norminf() const
    {
        CVM_ASSERT(this -> mpD, ((this -> _indofmax() - 1) * this -> mnIncr + 1) * sizeof(TC))
        return _abs (this -> mpD[(this -> _indofmax() - 1) * this -> mnIncr]);
    }

    virtual TR norm1() const
    {
        TR dNorm(0.);
        const int nSize = this -> mnSize * this -> mnIncr;
        for (int i = 0; i < nSize; i += this -> mnIncr)
        {
            dNorm += _abs(this -> mpD[i]);
        }
        return dNorm;
    }

    virtual TR norm2() const
    {
        return this -> norm();
    }

    const int* _pincr() const
    {
        return &this -> mnIncr;
    }

    const int* _psize() const
    {
        return &this -> mnSize;
    }

    // to be redefined in classes with non-traditional storage, like band matrices etc.
    virtual TC* _pd()
    {
        return this -> mpD;
    }

    // to be redefined in classes with non-traditional storage, like band matrices etc.
    virtual const TC* _pd() const
    {
        return this -> mpD;
    }

    void _div (TR d) throw (cvmexception)                                       // this = this / d for real only
    {
        static const TR one(1.);
        if (_abs(d) <= basic_cvmMachMin<TR>()) throw cvmexception (CVM_DIVISIONBYZERO);
        this -> _scal (one / d);
    }

protected:
    // 0-based
    virtual TC& _at (size_type n) throw (cvmexception)
    {
        if (n >= size_type(this -> mnSize)) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(this -> mpD, (n + 1) * sizeof(TC))
        return this -> mpD [n * this -> mnIncr];
    }

    // 0-based
    virtual const TC& _at (size_type n) const throw (cvmexception)
    {
        if (n >= size_type(this -> mnSize)) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(this -> mpD, (n + 1) * sizeof(TC))
        return this -> mpD [n * this -> mnIncr];
    }

    bool _equals (const Array& a) const                                         // compares array elements
    {
        bool bRes = false;
        if (this -> mnSize == a.size())
        {
            bRes = true;
            if (this -> mpD != a.get())
            {
                for (int i = 0; i < this -> mnSize; ++i)
                {
                    if (_abs (this -> mpD[i * this -> mnIncr] - a.mpD[i * a.mnIncr]) > basic_cvmMachMin<TR>())
                    {
                        bRes = false;
                        break;
                    }
                }
            }
        }
        return bRes;
    }

    void _normalize()                                                           // array normalizing
    {
        const TR dNorm = this -> norm();
        if (dNorm > basic_cvmMachMin<TR>())
        {
            static const TR one(1.);
            this -> _scal (one / dNorm);
        }
    }

    void _replace (const Array& a) throw (cvmexception)                         // this = a
    {
        cvmFree<TC> (this -> mpD);
        this -> mpD = cvmMalloc<TC>(a.mnSize);
        this -> mnSize = a.mnSize;
        this -> mnIncr = 1;
        CVM_ASSERT(this -> mpD, ((this -> mnSize - 1) * this -> mnIncr + 1) * sizeof(TC))
    }

// virtual methods
    virtual void _scal (TR d)
    {
        __scal<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
    }

    virtual int _indofmax() const                                               // index of max. module, undefined for matrices
    {
        return __idamax<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
    }

    virtual int _indofmin() const                                               // index of max. module, undefined for matrices
    {
        return __idamin<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
    }

    virtual void _set (TC d)                                                    // fills the content
    {
        CVM_ASSERT(this -> mpD, ((this -> mnSize - 1) * this -> mnIncr + 1) * sizeof(TC))
        const int nSize = this -> mnSize * this -> mnIncr;
        for (int i = 0; i < nSize; i += this -> mnIncr)
        {
            this -> mpD[i] = d;
        }
    }

    virtual void _assign (const TC* pD, int nIncr)
    {
        if (this -> mpD != pD)
        {
            __copy<TC> (this -> mnSize, pD, nIncr, this -> mpD, this -> mnIncr);
        }
    }

    virtual void _assign_shifted (TC* pDshifted, const TC* pD, int nSize, int nIncr, int)
    {
        if (pDshifted != pD)
        {
            __copy<TC> (nSize, pD, nIncr, pDshifted, this -> mnIncr);
        }
    }

    virtual void _resize (int nNewSize) throw (cvmexception)
    {
        if (nNewSize < 0) throw cvmexception (CVM_WRONGSIZE, nNewSize);
        if (nNewSize == 0)                                                      // just let's free object memory in this case
        {
            cvmFree<TC>(this -> mpD);
            this -> mnSize = this -> mnIncr = 0;
        }
        else if (nNewSize != this -> mnSize)
        {
            TC* pD = cvmMalloc<TC>(nNewSize);
            if (nNewSize > this -> mnSize) CleanMemory<TC> (pD, nNewSize);
            const int nMinSize = _cvm_min<int>(nNewSize, this -> mnSize);

            if (nMinSize > 0)
            {
                __copy<TC> (nMinSize, this -> mpD, this -> mnIncr, pD, 1);
            }
            cvmFree<TC>(this -> mpD);
            this -> mpD    = pD;
            this -> mnSize = nNewSize;
            this -> mnIncr = 1;
            CVM_ASSERT(this -> mpD, ((this -> mnSize - 1) * this -> mnIncr + 1) * sizeof(TC))
        }
    }

    friend std::istream& operator >> <> (std::istream& is, Array& aIn);
    friend std::ostream& operator << <> (std::ostream& os, const Array& aOut);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// basic_rvector
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TR>
class basic_rvector : public Array<TR,TR>
{
    typedef std::complex<TR> TC;
    typedef Array<TR,TR> BaseArray;

public:
    basic_rvector()
    {
    }

    explicit basic_rvector (int nSize)
        : BaseArray (nSize)
    {
    }

    basic_rvector (int nSize, TR d)
        : BaseArray (nSize, false)
    {
        this -> _set (d);
    }

    // WARNING!
    // The following constructor does not allocate memory!
    // It just shares a memory allocated before.
    // It is intented to make possible the following syntax:
    //
    // basic_rmatrix m (10, 20);
    // basic_rvector v (20);
    // ...
    // m[1] = v;            // assigns v to the 1st row of m
    //
    //
    // And for example this code...
    //
    // basic_rmatrix m (10,20);
    // basic_rvector vRow = m[1];
    //
    // ...will also call THIS constructor, and memory will be shared!

    // If you need the code like this with memory allocation, use the following:
    //
    // basic_rmatrix m (10,20);
    // basic_rvector vRow (m.msize());
    // vRow = m[1];
    //
    basic_rvector (TR* pD, int nSize, int nIncr = 1)
        : BaseArray (pD, nSize, nIncr)
    {
    }

    basic_rvector (const basic_rvector& v)
        : BaseArray (v.size(), false)
    {
        __copy<TR> (this -> mnSize, v, v.incr(), this -> mpD, this -> mnIncr);
    }

    basic_rvector& operator = (const basic_rvector& v) throw (cvmexception)     // assignment (equal sizes)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign (v, v.incr());
        return *this;
    }

    basic_rvector& assign (const TR* pD, int nIncr = 1)                         // foreign array assignment
    {
        this -> _assign (pD, nIncr);
        return *this;
    }

    basic_rvector& assign (int n, const TR* pD, int nIncr = 1)                  // foreign array assignment to tail (till the end)
    {                                                                           // beginning from n-th element (1-based)
        if (n <= 0 || n > this -> mnSize) throw cvmexception (CVM_OUTOFRANGE, n);
        --n;
        this -> _assign_shifted (this -> mpD + this -> mnIncr * n, pD, this -> mnSize - n, nIncr, 0);
        return *this;
    }

    basic_rvector& assign (int n, const basic_rvector& v) throw (cvmexception)  // subvector assignment
    {
        if (v.mnSize + n > this -> mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        return assign (n, v, v.incr());
    }

    basic_rvector& set (TR d)                                                   // sets the content
    {
        this -> _set (d);
        return *this;
    }

    basic_rvector& resize (int nNewSize) throw (cvmexception)
    {
        this -> _resize (nNewSize);
        return *this;
    }

    bool operator == (const basic_rvector& v) const
    {
        return this -> _equals (v);
    }

    bool operator != (const basic_rvector& v) const
    {
        return !(this -> operator == (v));
    }

    // vector replacement
    basic_rvector& operator << (const basic_rvector& v) throw (cvmexception)
    {
        this -> _replace (v);
        __copy<TR> (this -> mnSize, v.mpD, v.mnIncr, this -> mpD, this -> mnIncr);
        return *this;
    }

    // v1 + v2
    basic_rvector operator + (const basic_rvector& v) const throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rvector vSum (*this);
        __add<TR>(vSum.mpD, vSum.mnSize, vSum.mnIncr, v._pd(), v.incr());
        return vSum;
    }

    // v1 - v2
    basic_rvector operator - (const basic_rvector& v) const throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rvector vDiff (*this);
        __subtract<TR>(vDiff.mpD, vDiff.mnSize, vDiff.mnIncr, v._pd(), v.incr());
        return vDiff;
    }

    // this = v1 + v2
    basic_rvector& sum (const basic_rvector& v1, const basic_rvector& v2) throw (cvmexception)
    {
        if (this -> mnSize != v1.mnSize || this -> mnSize != v2.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _sum<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, v1, v1.incr(), v2, v2.incr());
        return *this;
    }

    // this = v1 + v2
    basic_rvector& diff (const basic_rvector& v1, const basic_rvector& v2) throw (cvmexception)
    {
        if (this -> mnSize != v1.mnSize || this -> mnSize != v2.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _diff<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, v1, v1.incr(), v2, v2.incr());
        return *this;
    }

    basic_rvector& operator += (const basic_rvector& v) throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _incr<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, v, v.incr());
        return *this;
    }

    basic_rvector& operator -= (const basic_rvector& v) throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _decr<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, v, v.incr());
        return *this;
    }

    basic_rvector operator - () const throw (cvmexception)
    {
        static const TR mone(-1.);
        basic_rvector vRes (*this);
        vRes._scal (mone);
        return vRes;
    }

    basic_rvector operator * (TR dMult) const throw (cvmexception)
    {
        basic_rvector vRes (*this);
        vRes._scal (dMult);
        return vRes;
    }

    basic_rvector operator / (TR dDiv) const throw (cvmexception)
    {
        basic_rvector vRes (*this);
        vRes._div (dDiv);
        return vRes;
    }

    // this = this * d
    basic_rvector& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    // this = this / d
    basic_rvector& operator /= (TR dDiv) throw (cvmexception)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_rvector& normalize()
    {
        this -> _normalize();
        return *this;
    }

    TR operator * (const basic_rvector& v) const throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        return __dot<TR>(this -> mpD, this -> mnSize, this -> mnIncr, v.mpD, v.mnIncr);
    }

    basic_rvector operator * (const basic_rmatrix<TR>& m) const throw (cvmexception)
    {
        if (this -> mnSize != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rvector vRes (m.nsize());
        m._multiply (vRes, *this, true);
        return vRes;
    }

    // this = v * m
    basic_rvector& mult (const basic_rvector& v, const basic_rmatrix<TR>& m) throw (cvmexception)
    {
        if (this -> mnSize != m.nsize() || v.size() != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        m._multiply (*this, v, true);
        return *this;
    }

    // this = m * v
    basic_rvector& mult (const basic_rmatrix<TR>& m, const basic_rvector& v) throw (cvmexception)
    {
        if (this -> mnSize != m.msize() || v.size() != m.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        m._multiply (*this, v, false);
        return *this;
    }

    // v_col * v_row (rank-1 update)
    basic_rmatrix<TR> rank1update (const basic_rvector& v) const
    {
        basic_rmatrix<TR> mRes (this -> mnSize, v.mnSize);
        static const TR one(1.);
        __ger<TR, basic_rmatrix<TR>, basic_rvector>(mRes, *this, v, one);
        return mRes;
    }

    // linear solvers
    basic_rvector& solve (const basic_srmatrix<TR>& mA, const basic_rvector& vB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != this -> mnSize || mA.msize() != vB.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (vB, *this, dErr, NULL, NULL);
        return *this;
    }

    basic_rvector& solve (const basic_srmatrix<TR>& mA, const basic_rvector& vB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (mA, vB, dErr);
    }

    basic_rvector& solve_lu (const basic_srmatrix<TR>& mA, const basic_srmatrix<TR>& mLU, const int* pPivots,
                             const basic_rvector& vB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != this -> size() || mA.msize() != vB.size() || mA.msize() != mLU.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (vB, *this, dErr, mLU, pPivots);
        return *this;
    }

    basic_rvector& solve_lu (const basic_srmatrix<TR>& mA, const basic_srmatrix<TR>& mLU, const int* pPivots, 
                             const basic_rvector& vB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mA, mLU, pPivots, vB, dErr);
    }

    // singular value decomposition
    basic_rvector& svd (const basic_rmatrix<TR>& m) throw (cvmexception)
    {
        m._svd (*this, NULL, NULL);
        return *this;
    }

    basic_rvector& svd (const basic_rmatrix<TR>& m, basic_srmatrix<TR>& mU,
                                                    basic_srmatrix<TR>& mVH) throw (cvmexception)
    {
        m._svd (*this, &mU, &mVH);
        return *this;
    }

  

    // eigenvalues
    // we don't use _eig here since this is the special case - symmetric matrix
    basic_rvector& eig (const basic_srsmatrix<TR>& m) throw (cvmexception)
    {
        __eig<basic_rvector, basic_srsmatrix<TR>, basic_srmatrix<TR> > (*this, m, NULL, true);
        return *this;
    }

    basic_rvector& eig (const basic_srsmatrix<TR>& m, basic_srmatrix<TR>& mEigVect) throw (cvmexception)
    {
        __eig<basic_rvector, basic_srsmatrix<TR>, basic_srmatrix<TR> > (*this, m, &mEigVect, true);
        return *this;
    }

 
    // ?gemv routines perform a matrix-vector operation defined as
    // this = alpha*m*v + beta * this,
    basic_rvector& gemv (bool bLeft, const basic_rmatrix<TR>& m, TR dAlpha, const basic_rvector& v, TR dBeta) throw (cvmexception)
    {
        if ((bLeft ? m.msize() != v.mnSize : m.nsize() != v.mnSize) ||
            (bLeft ? m.nsize() != this -> mnSize : m.msize() != this -> mnSize)) throw cvmexception (CVM_SIZESMISMATCH);
        m._gemv (bLeft, dAlpha, v, dBeta, *this);
        return *this;
    }

    basic_rvector& gbmv (bool bLeft, const basic_srbmatrix<TR>& m, TR dAlpha, const basic_rvector& v, TR dBeta) throw (cvmexception)
    {
        if ((bLeft ? m.msize() != v.mnSize : m.nsize() != v.mnSize) ||
            (bLeft ? m.nsize() != this -> mnSize : m.msize() != this -> mnSize)) throw cvmexception (CVM_SIZESMISMATCH);
        m._gbmv (bLeft, dAlpha, v, dBeta, *this);
        return *this;
    }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Array
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// generalized matrix
template <typename TR, typename TC>
class Matrix : public Array<TR,TC>
{
    typedef Array<TR,TC> BaseArray;

protected:
    int mnM;                                                                    // number of rows
    int mnN;                                                                    // number of columns
    int mnLD;                                                                   // leading dimension

    Matrix()
        : mnM(0), mnN(0), mnLD(0)
    {
    }

    Matrix (int nM, int nN, int nLD, bool bZeroMemory)
        : BaseArray (nLD * nN, bZeroMemory), mnM(nM), mnN(nN), mnLD(nLD)
    {
    }

    Matrix (TC* pD, int nM, int nN, int nLD, int nSize)         // for submatrices
        : BaseArray (pD, nSize, 1), mnM(nM), mnN(nN), mnLD(nLD)
    {
    }

    Matrix (const BaseArray& v, bool bBeColumn)                 // true = column, false = row
        : BaseArray (v.size()), mnM (bBeColumn ? v.size() : 1), mnN (bBeColumn ? 1 : v.size()), mnLD(mnM)
    {
        __copy<TC> (this -> mnSize, v, v.incr(), this -> mpD, this -> mnIncr);
    }

public:
    int msize() const
    {
        return mnM;
    }

    int nsize() const
    {
        return mnN;
    }

    int ld() const
    {
        return mnLD;
    }

    int rowofmax() const
    {
        return (this -> _indofmax() - 1) % mnM + 1;
    }

    int rowofmin() const
    {
        return (this -> _indofmin() - 1) % mnM + 1;
    }

    int colofmax() const
    {
        return (this -> _indofmax() - 1) / mnM + 1;
    }

    int colofmin() const
    {
        return (this -> _indofmin() - 1) / mnM + 1;
    }

    virtual TR norm1() const
    {
        int i, j, k;
        TR  rSum, rNorm(0.);

        for (j = 0; j < mnN; ++j)
        {
            rSum = TR(0.);

            k = j * mnLD;
            for (i = 0; i < mnM; ++i)
            {
                CVM_ASSERT(this -> mpD, (k + i + 1) * sizeof(TC))
                rSum += _abs (this -> mpD [k + i]);
            }

            if (rSum > rNorm)
            {
                rNorm = rSum;
            }
        }
        return rNorm;
    }

    virtual TR norminf() const
    {
        int i, j;
        TR  rSum, rNorm(0.);

        for (i = 0; i < mnM; ++i)
        {
            rSum = TR(0.);

            for (j = 0; j < mnN; ++j)
            {
                CVM_ASSERT(this -> mpD, (j * this -> mnLD + i + 1) * sizeof(TC))
                rSum += _abs (this -> mpD [j * this -> mnLD + i]);
            }

            if (rSum > rNorm)
            {
                rNorm = rSum;
            }
        }
        return rNorm;
    }

    const int* _pm() const
    {
        return &mnM;
    }

    const int* _pn() const
    {
        return &mnN;
    }

    const int* _pld() const
    {
        return &mnLD;
    }

    TC* _sub_pointer_nocheck (int row, int col)
    {
        return this -> _pd() + ((col - 1) * this -> ld() + row - 1);
    }

    TC* _sub_pointer (int row, int col, int height, int width) throw (cvmexception)
    {
        if (row    <= 0) throw cvmexception (CVM_WRONGSIZE, row);
        if (col    <= 0) throw cvmexception (CVM_WRONGSIZE, col);
        if (height <= 0) throw cvmexception (CVM_WRONGSIZE, height);
        if (width  <= 0) throw cvmexception (CVM_WRONGSIZE, width);
        if (row + height - 1 > mnM || col + width - 1 > mnN) throw cvmexception (CVM_SIZESMISMATCH);
        return _sub_pointer_nocheck (row, col);
    }

    virtual int _ldm() const
    {
        return this -> ld();
    }

    virtual const int* _pldm() const
    {
        return this -> _pld();
    }

    virtual bool _continuous () const
    {
        return mnM == mnLD;
    }

    void _check_ld() const
    {
        if (!this -> _continuous())
        {
            throw cvmexception (CVM_SUBMATRIXACCESSERROR);
        }
    }

    virtual void _scal (TR d)
    {
        if (this -> _continuous())
        {
            __scal<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            __scal<TR, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, d);
        }
    }

protected:
    virtual int _indofmax() const
    {
        this -> _check_ld();
        return __idamax<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
    }

    virtual int _indofmin() const
    {
        this -> _check_ld();
        return __idamin<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
    }

    virtual void _assign (const TC* pD, int nIncr)
    {
        if (this -> mpD != pD)
        {
            if (this -> _continuous())
            {
                __copy<TC> (this -> mnSize, pD, nIncr, this -> mpD, this -> mnIncr);
            }
            else for (int i = 0; i < this -> mnN; ++i)
            {
                __copy<TC> (this -> mnM, pD + this -> mnM * i * nIncr, nIncr, this -> mpD + this -> mnLD * i, this -> mnIncr);
            }
        }
    }

    virtual void _assign_shifted (TC* pDshifted, const TC* pD, int nRows, int nCols, int nLD)  // reusing nSise and nIncr parameter
    {
        if (pDshifted != pD)
        {
            for (int i = 0; i < nCols; ++i)
            {
                __copy<TC> (nRows, pD + nLD * i, 1, pDshifted + this -> mnLD * i, this -> mnIncr);
            }
        }
    }

    virtual void _set (TC d)
    {
        CVM_ASSERT(this -> mpD, this -> mnSize * sizeof(TC))
        int i, j, k;
        for (j = 0; j < this -> mnN; ++j)
        {
            k = j * mnLD;
            for (i = 0; i < this -> mnM; ++i)
            {
                this -> mpD[k + i] = d;
            }
        }
    }

    virtual void _massign (const Matrix& m)
    {
        if (this -> mpD != m.mpD)
        {
            if (this -> _continuous() && m._continuous())
            {
                __copy<TC> (this -> mnSize, m._pd(), m.incr(), this -> mpD, this -> mnIncr);
            }
            else
            {
                const TC* p = m._pd();
                const int nLD = m._ldm();
                for (int i = 0; i < this -> mnN; ++i)
                {
                    __copy<TC> (this -> mnM, p + nLD * i, m.incr(), this -> mpD + this -> mnLD * i, this -> mnIncr);
                }
            }
        }
    }

    virtual void _resize (int nNewM, int nNewN) throw (cvmexception)
    {
        if (nNewM != mnM || nNewN != mnN)
        {
            this -> _check_ld();

            if (nNewM < 0) throw cvmexception (CVM_WRONGSIZE, nNewM);
            if (nNewN < 0) throw cvmexception (CVM_WRONGSIZE, nNewN);
            const int nNewSize = nNewM * nNewN;

            if (nNewSize == 0)
            {
                cvmFree<TC>(this -> mpD);
                this -> mnSize = this -> mnIncr = this -> mnM = this -> mnN = this -> mnLD = 0;
            }
            else
            {
                TC* pD = cvmMalloc<TC>(nNewSize);
                if (nNewSize > this -> mnSize) CleanMemory<TC> (pD, nNewSize);
                const int nMinM = _cvm_min<int>(nNewM, this -> mnM);
                const int nMinN = _cvm_min<int>(nNewN, this -> mnN);

                for (int i = 0; i < nMinN; ++i)
                {
                    __copy<TC> (nMinM, this -> mpD + i * this -> mnM, this -> mnIncr, pD + i * nNewM, 1);
                }
                cvmFree<TC>(this -> mpD);
                this -> mpD    = pD;
                this -> mnSize = nNewSize;
                CVM_ASSERT(this -> mpD, this -> mnSize * sizeof(TC))

                this -> mnM    = nNewM;
                this -> mnN    = nNewN;
                this -> mnLD   = nNewM;
                this -> mnIncr = 1;
            }
        }
    }

    virtual int _ld_for_replace () const
    {
        return this -> mnLD;
    }

    virtual int _size_for_replace () const
    {
        return this -> mnSize;
    }

    void _replace (const Matrix& m) throw (cvmexception)
    {
        this -> _check_ld();                                    // submatrix replacement is obviously not possible
        cvmFree<TC>(this -> mpD);
        this -> mnSize = m._size_for_replace();
        this -> mpD = cvmMalloc<TC>(this -> mnSize);
        this -> mnIncr = 1;
        CVM_ASSERT(this -> mpD, (this -> mnSize * sizeof(TC)))
        this -> mnM  = m.mnM;
        this -> mnN  = m.mnN;
        this -> mnLD = m._ld_for_replace();
    }

    void _transp (const Matrix& m)
    {
        int i;
        if (this -> mnM > this -> mnN) for (i = 0; i < this -> mnN; ++i)
        {
            __copy<TC> (m.nsize(), m.get() + i, m.ld(), this -> mpD + i * this -> mnLD, 1);
        }
        else for (i = 0; i < this -> mnM; ++i)
        {
            __copy<TC> (m.msize(), m.get() + i * m.ld(), 1, this -> mpD + i, this -> mnLD);
        }
    }

    virtual type_proxy<TC,TR> _ij_proxy_val (int i, int j)                      // zero based
    {
        CVM_ASSERT(this -> mpD, (this -> mnLD * j + i + 1) * sizeof(TC))
        return type_proxy<TC,TR>(this -> mpD [this -> mnLD * j + i], false);
    }

    virtual TC _ij_val (int i, int j) const                                     // zero based
    {
        CVM_ASSERT(this -> mpD, (this -> mnLD * j + 1) * sizeof(TC))
        return this -> mpD [this -> mnLD * j + i];
    }

    virtual void _swap_rows (int n1, int n2) throw (cvmexception)
    {
        if (n1 <= 0 || n1 > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, n1);
        if (n2 <= 0 || n2 > this -> mnM) throw cvmexception (CVM_OUTOFRANGE2, n2);
        if (n1 != n2)
        {
            __swap<TC> (this -> mnN, this -> mpD + n1 - 1, this -> mnLD, this -> mpD + n2 - 1, this -> mnLD);
        }
    }

    virtual void _swap_cols (int n1, int n2) throw (cvmexception)
    {
        if (n1 <= 0 || n1 > mnN) throw cvmexception (CVM_OUTOFRANGE1, n1);
        if (n2 <= 0 || n2 > mnN) throw cvmexception (CVM_OUTOFRANGE2, n2);
        if (n1 != n2)
        {
            __swap<TC> (this -> mnM, this -> mpD + (n1 - 1) * this -> mnLD, 1, this -> mpD + (n2 - 1) * this -> mnLD, 1);
        }
    }

    virtual const TC* _p (const Matrix& m) const
    {
        return m._pd();
    }

    virtual void _msum (const Matrix& m1, const Matrix& m2)
    {
        if (this -> _continuous() && m1._continuous() && m2._continuous())
        {
            _sum<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, this -> _p(m1), m1.incr(), this -> _p(m2), m2.incr());
        }
        else for (int i = 0; i < mnN; ++i)
        {
            _sum<TR, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, this -> _p(m1) + m1._ldm() * i, m1.incr(), this -> _p(m2) + m2._ldm() * i, m2.incr());
        }
    }

    void _mdiff (const Matrix& m1, const Matrix& m2)
    {
        if (this -> _continuous() && m1._continuous() && m2._continuous())
        {
            _diff<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, this -> _p(m1), m1.incr(), this -> _p(m2), m2.incr());
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            _diff<TR, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, this -> _p(m1) + m1._ldm() * i, m1.incr(), this -> _p(m2) + m2._ldm() * i, m2.incr());
        }
    }

    virtual void _mincr (const Matrix& m)
    {
        if (this -> _continuous() && m._continuous())
        {
            _incr<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, this -> _p(m), m.incr());
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            _incr<TR, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, this -> _p(m) + m._ldm() * i, m.incr());
        }
    }

    virtual void _mdecr (const Matrix& m)
    {
        if (this -> _continuous() && m._continuous())
        {
            _decr<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, this -> _p(m), m.incr());
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            _decr<TR, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, this -> _p(m) + m._ldm() * i, m.incr());
        }
    }

    // matrix cleaning (we ALWAYS have mnIncr = 1 for matrices)
    virtual void _vanish()
    {
        CVM_ASSERT(this -> mpD, this -> mnSize * sizeof(TC))
        if (this -> _continuous())
        {
            memset (this -> mpD, 0, this -> mnSize * sizeof(TC));
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            memset (this -> mpD + this -> mnLD * i, 0, this -> mnM * sizeof(TC));
        }
    }

public:
    friend std::ostream& operator << <> (std::ostream& os, const Matrix<TR,TC>& mOut);
    friend std::istream& operator >> <> (std::istream& os, Matrix<TR,TC>& mIn);
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SqMatrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// generalized square matrix
template <typename TR, typename TC>
class SqMatrix
{
protected:
    SqMatrix()
    {
    }

    virtual ~SqMatrix()
    {
    }

    virtual int _size    () const = 0;
    virtual int _msize   () const = 0;
    virtual int _nsize   () const = 0;
    virtual int _ld      () const = 0;
    virtual const TC* _p () const = 0;
    virtual TC* _p () = 0;

    // it differs from Matrix::_transp because in this case we can do it in-place.
    void _sq_transp()
    {
        const int mnM  = this -> _msize();
        const int mnLD = this -> _ld();
        TC* mpD = this -> _p();
        if (mnM > 1)
        {
            const int nM1 = mnLD + 1, nM1m = mnLD - 1, nM2m = mnM - 1;
            int i = 1, j = 1, m;
            for (;;)
            {
                m = mnM - i;
                __swap<TC> (m, mpD + j, 1, mpD + j + nM1m, mnLD);
                if (i >= nM2m)
                {
                    break;
                }
                ++i;
                j += nM1;
            }
        }
    }

    void _sq_plus_plus()                                        // plus identity
    {
        TC* mpD = this -> _p();
        static const TC one(1.);
        const int mnSize = this -> _size();
        const int nNext = this -> _msize() + 1;
        CVM_ASSERT(mpD, mnSize * sizeof(TC))
        for (int i = 0; i < mnSize; i += nNext)
        {
            mpD[i] += one;
        }
    }

    void _sq_minus_minus()                                      // minus identity
    {
        TC* mpD = this -> _p();
        static const TC one(1.);
        const int mnSize = this -> _size();
        const int nNext = this -> _msize() + 1;
        CVM_ASSERT(mpD, mnSize * sizeof(TC))
        for (int i = 0; i < mnSize; i += nNext)
        {
            mpD[i] -= one;
        }
    }

public:
    void _clean_low_triangle ()
    {
        const int mnM  = this -> _msize();
        const int mnLD = this -> _ld();
        TC* mpD = this -> _p();
        int n = 1;
        static const TR zero(0.);
        for (int i = 1; i < mnM; ++i)
        {
            __scal<TR,TC> (mpD + n, mnM - i, 1, zero);   // column by column
            n += mnLD + 1;
        }
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// rmatrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// matrix of real numbers
template <typename TR>
class basic_rmatrix : public Matrix <TR,TR>
{
    typedef std::complex<TR> TC;
    typedef Array<TR,TR> BaseArray;
    typedef Matrix<TR,TR> BaseMatrix;
    typedef basic_rvector<TR> RVector;

    friend class basic_rvector<TR>; // _multiply

public:
    basic_rmatrix()
    {
    }

    basic_rmatrix (int nM, int nN)
        : BaseMatrix (nM, nN, nM, true)
    {
    }

    basic_rmatrix (TR* pD, int nM, int nN)
        : BaseMatrix (pD, nM, nN, nM, nM * nN)
    {
    }

    basic_rmatrix (const basic_rmatrix& m)
        : BaseMatrix (m.mnM, m.mnN, m.mnM, false)
    {
        this -> _massign(m);
    }

    explicit basic_rmatrix (const RVector& v, bool bBeColumn = true)            // true = column, false = row
        : BaseMatrix (v, bBeColumn)
    {
    }

    // submatrix constructor, 1-based
    basic_rmatrix (basic_rmatrix& m, int nRow, int nCol, int nHeight, int nWidth)
        : BaseMatrix (m._sub_pointer (nRow, nCol, nHeight, nWidth), nHeight, nWidth, m.ld(), nHeight * nWidth)
    {
        m._check_submatrix();
    }

    type_proxy<TR,TR> operator () (int nIm, int nIn) throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_proxy_val (nIm - 1, nIn - 1);
    }

    TR operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    RVector operator () (int nI) throw (cvmexception)                          // returns COLUMN which IS an l-value
    {
        if (nI <= 0 || nI > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _col (nI - 1);
    }

    RVector operator [] (int nI) throw (cvmexception)                          // returns ROW which IS an l-value
    {
        if (nI <= 0 || nI > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _row (nI - 1);
    }

    const RVector operator () (int nI) const throw (cvmexception)              // returns column which IS NOT an l-value
    {
        if (nI <= 0 || nI > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _col (nI - 1);
    }

    const RVector operator [] (int nI) const throw (cvmexception)              // returns ROW which IS NOT an l-value
    {
        if (nI <= 0 || nI > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _row (nI - 1);
    }

    RVector diag (int nDiag) throw (cvmexception)                               // returns diagonal which IS l-value (shares memory) 
    {                                                                           // 0 - main, negative - low, positive - up
        return this -> _diag(nDiag);
    }

    const RVector diag (int nDiag) const throw (cvmexception)                   // returns diagonal which IS NOT l-value
    {                                                                           // 0 - main, negative - low, positive - up
        return this -> _diag(nDiag);
    }

    basic_rmatrix& operator = (const basic_rmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        return *this;
    }

    basic_rmatrix& assign (const RVector& v)                                    // assigns vector
    {
        this -> _assign (v, v.incr());
        return *this;
    }

    basic_rmatrix& assign (const TR* pD)                                        // assigns foregn array (nIncr = 1)
    {
        this -> _assign (pD, 1);
        return *this;
    }

    basic_rmatrix& assign (int nRow, int nCol, const basic_rmatrix& m) throw (cvmexception)    // submatrix assignment
    {
        if (nRow <= 0 || nRow > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nRow);
        if (nCol <= 0 || nCol > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nCol);
        if (m.mnM + nRow - 1 > this -> mnM || m.mnN + nCol - 1 > this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> _sub_pointer_nocheck (nRow, nCol), m._pd(), m.mnM, m.mnN, m.mnLD);
        return *this;
    }

    basic_rmatrix& set (TR d)
    {
        this -> _set (d);
        return *this;
    }

    basic_rmatrix& resize (int nNewM, int nNewN) throw (cvmexception)
    {
        this -> _resize (nNewM, nNewN);
        return *this;
    }

    bool operator == (const basic_rmatrix& m) const
    {
        return this -> mnM == m.mnM && this -> mnN == m.mnN && this -> _mequals (m);
    }

    bool operator != (const basic_rmatrix& m) const
    {
        return !operator == (m);
    }

    basic_rmatrix& operator << (const basic_rmatrix& m) throw (cvmexception)    // matrix replacement
    {
        this -> _replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_rmatrix operator + (const basic_rmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_rmatrix operator - (const basic_rmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    // this = v1 + v2
    basic_rmatrix& sum (const basic_rmatrix& m1, const basic_rmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnN != m1.mnN || this -> mnM != m2.mnM || this -> mnN != m2.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _msum (m1, m2);
        return *this;
    }

    // this = v1 + v2
    basic_rmatrix& diff (const basic_rmatrix& m1, const basic_rmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnN != m1.mnN || this -> mnM != m2.mnM || this -> mnN != m2.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_rmatrix& operator += (const basic_rmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mincr (m);
        return *this;
    }

    basic_rmatrix& operator -= (const basic_rmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdecr (m);
        return *this;
    }

    basic_rmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_rmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    basic_rmatrix operator * (TR dMult) const
    {
        basic_rmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_rmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_rmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_rmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_rmatrix& operator /= (TR dDiv) throw (cvmexception)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_rmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    // transposed Matrix
    basic_rmatrix operator ~ () const throw (cvmexception)
    {
        basic_rmatrix mRes (this -> mnN, this -> mnM);
        mRes._transp (*this);
        return mRes;
    }

    basic_rmatrix& transpose (const basic_rmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnN || this -> mnN != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        if (this -> mpD == m.mpD)
        {
            basic_rmatrix mTmp(m);
            this -> _transp (mTmp);
        }
        else
        {
            this -> _transp (m);
        }
        return *this;
    }

    // in-place transpose Matrix - changes dimensions
    basic_rmatrix& transpose() throw (cvmexception)
    {
        basic_rmatrix mTmp (*this);
        this -> _resize (this -> mnN, this -> mnM);
        this -> _transp (mTmp);
        return *this;
    }

    RVector operator * (const RVector& v) const throw (cvmexception)
    {
        if (this -> mnN != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector vRes (this -> mnM);
        this -> _multiply (vRes, v, false);
        return vRes;
    }

    basic_rmatrix operator * (const basic_rmatrix& m) const throw (cvmexception)
    {
        if (this -> mnN != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rmatrix mRes (this -> mnM, m.mnN);
        mRes.mult (*this, m);
        return mRes;
    }

    // this = m1 * m2
    // overridden in srsmatrix, NOT in srmatrix
    basic_rmatrix& mult (const basic_rmatrix& m1, const basic_rmatrix& m2) throw (cvmexception)
    {
        this -> _mult (m1, m2);
        return *this;
    }

    // this = v_col * v_row (rank-1 update)
    basic_rmatrix& rank1update (const RVector& vCol, const RVector& vRow) throw (cvmexception)
    {
        static const TR one(1.);
        this -> _check_rank1update();
        if (this -> mnM != vCol.size() || this -> mnN != vRow.size()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _vanish();
        __ger<TR, basic_rmatrix, RVector> (*this, vCol, vRow, one);
        return *this;
    }

    basic_rmatrix& swap_rows (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_rows (n1, n2);
        return *this;
    }

    basic_rmatrix& swap_cols (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_cols (n1, n2);
        return *this;
    }

    // linear solvers
    basic_rmatrix& solve (const basic_srmatrix<TR>& mA, const basic_rmatrix& mB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != mB.msize() || mA.msize() != this -> msize() || mB.nsize() != this -> nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (mB, *this, dErr, NULL, NULL);
        return *this;
    }

    basic_rmatrix& solve (const basic_srmatrix<TR>& mA, const basic_rmatrix& mB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (mA, mB, dErr);
    }

    basic_rmatrix& solve_lu (const basic_srmatrix<TR>& mA, const basic_srmatrix<TR>& mLU, const int* pPivots,
                             const basic_rmatrix& mB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != mB.msize()  || mA.msize() != this -> msize() ||
            mA.msize() != mLU.msize() || mB.nsize() != this -> nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (mB, *this, dErr, mLU, pPivots);
        return *this;
    }

    basic_rmatrix& solve_lu (const basic_srmatrix<TR>& mA, const basic_srmatrix<TR>& mLU, const int* pPivots,
                             const basic_rmatrix& mB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mA, mLU, pPivots, mB, dErr);
    }

    // singular value decomposition (values in decreasing order)
    RVector svd() const throw (cvmexception)
    {
        RVector vRes (_cvm_min<int>(this -> mnM, this -> mnN));
        this -> _svd (vRes, NULL, NULL);
        return vRes;
    }

    RVector svd (basic_srmatrix<TR>& mU, basic_srmatrix<TR>& mVH) const throw (cvmexception)
    {
        RVector vRes (_cvm_min<int>(this -> mnM, this -> mnN));
        this -> _svd (vRes, &mU, &mVH);
        return vRes;
    }

    // pseudo (generalized) inversion - const version
    basic_rmatrix pinv (TR threshold = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        basic_rmatrix mAx(this -> mnN, this -> mnM);
        this -> _pinv (mAx, threshold);
        return mAx;
    }

    // pseudo (generalized) inversion - non-const version
    basic_rmatrix& pinv (const basic_rmatrix& mA, TR threshold = basic_cvmMachSp<TR>()) throw (cvmexception)
    {
        if (mA.msize() != this -> nsize() || mA.nsize() != this -> msize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._pinv (*this, threshold);
        return *this;
    }

    int rank (TR eps = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        int nRank = 0;
        RVector vS (_cvm_min<int>(this -> mnM, this -> mnN));
        this -> _svd (vS, NULL, NULL);
        vS.normalize();
        for (; nRank < vS.size(); ++nRank)
        {
            if (vS [nRank * this -> mnIncr + 1] < eps) break;
        }
        return nRank;
    }

    // QR decomposition 
    // Case 1: economy mode, A is (m x n) and Q is (m x n) and R is (n x n)
    void qr (basic_rmatrix<TR>& mQ, basic_srmatrix<TR>& mR) const throw (cvmexception)
    {
        this -> _qr(mQ, mR);
    }
    // Case 2: full mode, A is (m x n) and Q is (m x m) and R is (m x n)
    void qr (basic_srmatrix<TR>& mQ, basic_rmatrix<TR>& mR) const throw (cvmexception)
    {
        this -> _qr(mQ, mR);
    }

    basic_rmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    // this += alpha * v_col * v_row (rank-1 update)
    basic_rmatrix& ger (TR alpha, const RVector& vCol, const RVector& vRow) throw (cvmexception)
    {
        this -> _check_ger();
        if (this -> mnM != vCol.size() || this -> mnN != vRow.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __ger<TR, basic_rmatrix, RVector> (*this, vCol, vRow, alpha);
        return *this;
    }

    // this = alpha * m1 * m2 + beta * this
    basic_rmatrix& gemm (const basic_rmatrix& m1, bool bTrans1, const basic_rmatrix& m2, bool bTrans2, TR dAlpha, TR dBeta) throw (cvmexception)
    {
        this -> _check_gemm();
        if (this -> mnM != (bTrans1 ? m1.mnN : m1.mnM) || 
            this -> mnN != (bTrans2 ? m2.mnM : m2.mnN) || 
            (bTrans1 ? m1.mnM : m1.mnN) != (bTrans2 ? m2.mnN : m2.mnM)) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _gemm (bTrans1, m1, bTrans2, m2, dAlpha, dBeta);
        return *this;
    }

    // this = alpha*a*b + beta*this of this = alpha*b*a + beta*this  where a is symmetric
    basic_rmatrix& symm (bool bLeft, const basic_srsmatrix<TR>& ms, const basic_rmatrix& m, TR dAlpha, TR dBeta) throw (cvmexception)
    {
        this -> _check_symm();
        if (this -> mnM != m.mnM || this -> mnN != m.mnN ||
            ms.msize() != (bLeft ? this -> mnM : this -> mnN)) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _symm (bLeft, ms, m, dAlpha, dBeta);
        return *this;
    }

    // 2-norm (maximum singular value)
    virtual TR norm2() const
    {
        RVector vS (_cvm_min<int>(this -> mnM, this -> mnN));
        this -> _svd (vS, NULL, NULL);
        return vS[1];
    }

    virtual void _svd (RVector& vRes, basic_srmatrix<TR>* pmU, basic_srmatrix<TR>* pmVH) const throw (cvmexception)
    {
        if (pmU != NULL && pmVH != NULL && (this -> mnM != pmU -> msize() || this -> mnN != pmVH -> msize())) throw cvmexception (CVM_SIZESMISMATCH);
        __svd<TR, basic_rmatrix, basic_srmatrix<TR> > (vRes, vRes.size(), vRes.incr(), *this, pmU, pmVH);
    }

    virtual void _pinv (basic_rmatrix& mX, TR threshold) const throw (cvmexception)
    {
        __pinv<TR, basic_rmatrix, basic_rmatrix> (mX, *this, threshold);
    }

    // ?gemm routines perform a matrix-matrix operation with general matrices. The operation is defined as
    // c := alpha*op(a)*op(b) + beta*c,
    // where: op(x) is one of op(x) = x or op(x) = x' or op(x) = conjg(x'),
    void _gemm (bool bTrans1, const basic_rmatrix& m1, bool bTrans2, const basic_rmatrix& m2, TR dAlpha, TR dBeta) throw (cvmexception)     // this = m1 * m2
    {
        basic_rmatrix mTmp1, mTmp2;
        const TR* pD1 = m1.get();
        const TR* pD2 = m2.get();
        if (this -> mpD == pD1) mTmp1 << m1;
        if (this -> mpD == pD2) mTmp2 << m2;
        __gemm<TR, basic_rmatrix> (this -> mpD == pD1 ? mTmp1 : m1, bTrans1, this -> mpD == pD2 ? mTmp2 : m2, bTrans2, dAlpha, *this, dBeta);
    }

    // this = alpha*a*b + beta*this or this = alpha*b*a + beta*this  where a is symmetric
    void _symm (bool bLeft, const basic_srsmatrix<TR>& ms, const basic_rmatrix& m, TR dAlpha, TR dBeta) throw (cvmexception)
    {
        basic_rmatrix mTmp;
        basic_srsmatrix<TR> msTmp;
        const TR* pD1 = ms.get();
        const TR* pD2 = m._pd();
        if (this -> mpD == pD1) msTmp << ms;
        if (this -> mpD == pD2) mTmp << m;
        __symm<TR, basic_srsmatrix<TR>, basic_rmatrix> (bLeft, this -> mpD == pD1 ? msTmp : ms, this -> mpD == pD2 ? mTmp : m, dAlpha, *this, dBeta);
    }

    virtual void _check_submatrix () {}

protected:
    // protected constructors for inherited stuff
    basic_rmatrix (int nM, int nN, int nLD, bool bZeroMemory)
        : BaseMatrix (nM, nN, nLD, bZeroMemory)
    {
    }

    basic_rmatrix (TR* pD, int nM, int nN, int nLD, int nSize)
        : BaseMatrix (pD, nM, nN, nLD, nSize)
    {
    }

    // returns diagonal which IS l-value (shares memory)
    // 0 - main, negative - low, positive - up
    virtual RVector _diag (int nDiag) throw (cvmexception)
    {
        int nShift = 0;
        int nSize = 0;
        if (nDiag >= 0)
        {
            if (nDiag >= this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nDiag);
            nShift = nDiag * this -> mnLD;
            nSize = this -> mnN > this -> mnM ? (nDiag > this -> mnN - this -> mnM ? this -> mnN - nDiag : this -> mnM) : this -> mnN - nDiag;
        }
        else
        {
            nShift = - nDiag;
            if (nShift >= this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nDiag);
            nSize = this -> mnM > this -> mnN ? (nShift > this -> mnM - this -> mnN ? this -> mnM - nShift : this -> mnN) : this -> mnM - nShift;
        }
        return RVector (this -> mpD + nShift, nSize, this -> mnLD + 1);
    }

    // returns diagonal which IS NOT l-value (creates a copy)
    // 0 - main, negative - low, positive - up
    virtual const RVector _diag (int nDiag) const throw (cvmexception)
    {
        RVector vRet (const_cast<basic_rmatrix*>(this) -> _diag(nDiag));
        return vRet;
    }

    // compares matrix elements (equal sizes assumed)
    bool _mequals (const basic_rmatrix& m) const
    {
        return ((*this) - m).norminf() <= basic_cvmMachMin<TR>();
    }

    // ?gemv routines perform a matrix-vector operation defined as
    // vRes = alpha*m*v + beta * vRes or vRes = alpha*v'*m + beta * vRes
    // not virtual since __gemv calls all virtual methods inside
    void _gemv (bool bLeft, TR dAlpha, const RVector& v, TR dBeta, RVector& vRes) const
    {
        RVector vTmp;
        basic_rmatrix mTmp;
        const TR* pDv = v;
        if (vRes.get() == pDv) vTmp << v;
        if (vRes.get() == this -> mpD) mTmp << *this;
        __gemv<TR, basic_rmatrix, RVector> (bLeft, vRes.get() == this -> mpD ? mTmp : *this, dAlpha, 
                                                   vRes.get() == pDv ? vTmp : v, dBeta, vRes);
    }

    // 0-based, returns l-value sharing memory
    virtual RVector _row (int m)
    {
        return RVector (this -> mpD + m, this -> mnN, this -> mnLD);
    }

    // 0-based, returns l-value sharing memory
    virtual RVector _col (int n)
    {
        return RVector (this -> mpD + this -> mnLD * n, this -> mnM);
    }

    // 0-based
    virtual const RVector _row (int m) const
    {
        return RVector (this -> mpD + m, this -> mnN, this -> mnLD);
    }

    // 0-based, returns l-value sharing memory
    virtual const RVector _col (int n) const
    {
        return RVector (this -> mpD + this -> mnLD * n, this -> mnM);
    }

    virtual void _mult (const basic_rmatrix& m1, const basic_rmatrix& m2) throw (cvmexception)
    {
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        if (this -> mnM != m1.mnM || this -> mnN != m2.mnN || m1.mnN != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _gemm (false, m1, false, m2, one, zero);
    }

    virtual void _multiply (RVector& vRes, const RVector& v, bool bLeft) const
    {
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        this -> _gemv (bLeft, one, v, zero, vRes);
    }

    // QR decomposition
    // Case 1: economy mode, A is (m x n) and Q is (m x n) and R is (n x n)
    virtual void _qr(basic_rmatrix<TR>& pmQ, basic_srmatrix<TR>& pmR) const throw (cvmexception)
    {
        if (this -> mnM != pmQ.msize() || this -> mnN != pmQ.nsize() || this -> mnN != pmR.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        __qre<basic_rmatrix, basic_srmatrix<TR> > (*this, pmQ, pmR);
    }

    // Case 2: full mode, A is (m x n) and Q is (m x m) and R is (m x n)
    virtual void _qr(basic_srmatrix<TR>& pmQ, basic_rmatrix<TR>& pmR) const throw (cvmexception)
    {
        if (this -> mnM != pmQ.msize() || this -> mnM != pmR.msize() || this -> mnN != pmR.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        __qrf<basic_rmatrix, basic_srmatrix<TR> > (*this, pmQ, pmR);
    }

    virtual void _check_ger() {}
    virtual void _check_rank1update() {}
    virtual void _check_gemm() {}
    virtual void _check_symm() {}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SrMatrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TR>
class basic_srmatrix : public basic_rmatrix<TR>, public SqMatrix <TR, TR>
{
    typedef std::complex<TR> TC;
    typedef basic_rvector<TR> RVector;
    typedef Array<TR,TR> BaseArray;
    typedef Matrix<TR,TR> BaseMatrix;
    typedef SqMatrix<TR, TR> BaseSqMatrix;
    typedef basic_rmatrix<TR> BaseRMatrix;

public:
    basic_srmatrix()
    {
    }

    explicit basic_srmatrix (int nMN)
        : BaseRMatrix (nMN, nMN)
    {
    }

    basic_srmatrix (TR* pD, int nMN)
        : BaseRMatrix (pD, nMN, nMN)
    {
    }

    basic_srmatrix (const basic_srmatrix& m)
        : BaseRMatrix (m.mnM, m.mnN, m.mnM, false), BaseSqMatrix()
    {
        this -> _massign(m);
    }

    basic_srmatrix (const BaseRMatrix& m)
        : BaseRMatrix (m.msize(), m.nsize(), m.msize(), false)
    {
        if (this -> mnM != this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
    }

    // diagonal square matrix constructor
    explicit basic_srmatrix (const RVector& v)
        : BaseRMatrix (v.size(), v.size(), v.size(), true)
    {
        __copy<TR> (this -> mnM, v, v.incr(), this -> mpD, this -> mnM + 1);
    }

    // submatrix constructor
    // 1-based
    basic_srmatrix (BaseRMatrix& m, int nRow, int nCol, int nSize)
        : BaseRMatrix (m, nRow, nCol, nSize, nSize)
    {
        m._check_submatrix();
    }

    type_proxy<TR,TR> operator () (int nIm, int nIn) throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_proxy_val (nIm - 1, nIn - 1);
    }

    TR operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    // returns COLUMN which CAN be l-value
    RVector operator () (int nFI) throw (cvmexception)
    {
        return this -> _col (nFI - 1);
    }

    // returns column which CAN NOT be l-value
    const RVector operator () (int nFI) const throw (cvmexception)
    {
        return this -> _col (nFI - 1);
    }

    // returns ROW which CAN be l-value
    RVector operator [] (int nFI) throw (cvmexception)
    {
        return this -> _row (nFI - 1);
    }

    // returns ROW which CAN NOT be l-value
    const RVector operator [] (int nFI) const throw (cvmexception)
    {
        return this -> _row (nFI - 1);
    }

    basic_srmatrix& operator = (const basic_srmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        return *this;
    }

    basic_srmatrix& assign (const RVector& v)                                   // assigns vector
    {
        this -> _assign (v, v.incr());
        return *this;
    }

    basic_srmatrix& assign (const TR* pD)                                       // assigns foregn array (nIncr = 1)
    {
        this -> _assign (pD, 1);
        return *this;
    }

    basic_srmatrix& assign (int nRow, int nCol, const BaseRMatrix& m) throw (cvmexception)      // submatrix assignment
    {
        if (nRow <= 0 || nRow > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nRow);
        if (nCol <= 0 || nCol > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nCol);
        if (m.msize() + nRow - 1 > this -> mnM || m.nsize() + nCol - 1 > this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> _sub_pointer_nocheck (nRow, nCol), m._pd(), m.msize(), m.nsize(), m.ld());
        return *this;
    }

    // fills the content
    basic_srmatrix& set (TR d)
    {
        this -> _set (d);
        return *this;
    }

    basic_srmatrix& resize (int nNewMN) throw (cvmexception)
    {
        this -> _resize (nNewMN, nNewMN);
        return *this;
    }

    basic_srmatrix& operator << (const basic_srmatrix& m) throw (cvmexception)
    {
        this -> _replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_srmatrix operator + (const basic_srmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_srmatrix operator - (const basic_srmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    // this = v1 + v2
    basic_srmatrix& sum (const basic_srmatrix& m1, const basic_srmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _msum (m1, m2);
        return *this;
    }

    // this = v1 + v2
    basic_srmatrix& diff (const basic_srmatrix& m1, const basic_srmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_srmatrix& operator += (const basic_srmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mincr (m);
        return *this;
    }

    basic_srmatrix& operator -= (const basic_srmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdecr (m);
        return *this;
    }

    basic_srmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_srmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    // plus identity, prefix
    basic_srmatrix& operator ++ ()
    {
        this -> _plus_plus();
        return *this;
    }

    // plus identity, postfix
    basic_srmatrix& operator ++ (int)
    {
        this -> _plus_plus();
        return *this;
    }

    // minus identity, prefix
    basic_srmatrix& operator -- ()
    {
        this -> _minus_minus();
        return *this;
    }

    // minus identity, postfix
    basic_srmatrix& operator -- (int)
    {
        this -> _minus_minus();
        return *this;
    }

    basic_srmatrix operator * (TR dMult) const
    {
        basic_srmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_srmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_srmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_srmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_srmatrix& operator /= (TR dDiv)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_srmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    basic_srmatrix operator ~ () const throw (cvmexception)
    {
        basic_srmatrix mRes (*this);
        return mRes.transpose();
    }

    basic_srmatrix& transpose (const basic_srmatrix& m) throw (cvmexception)
    {
        (*this) = m;
        return this -> transpose();
    }

    basic_srmatrix& transpose()
    {
        this -> _transp();
        return *this;
    }

    RVector operator * (const RVector& v) const throw (cvmexception)
    {
        return this -> BaseRMatrix::operator * (v);
    }

    // special exclusion since matrix product is not commutative
    BaseRMatrix operator * (const BaseRMatrix& m) const throw (cvmexception)
    {
        return this -> BaseRMatrix::operator * (m);
    }

    basic_srmatrix operator * (const basic_srmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srmatrix mRes (this -> mnM);
        mRes.mult (*this, m);
        return mRes;
    }

    basic_srmatrix& operator *= (const basic_srmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        const basic_srmatrix mTmp (*this);
        this -> mult (mTmp, m);
        return *this;
    }

    basic_srmatrix& swap_rows (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_rows (n1, n2);
        return *this;
    }

    basic_srmatrix& swap_cols (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_cols (n1, n2);
        return *this;
    }

    // linear solvers Ax=b. Also return solution flavor.
    RVector solve (const RVector& vB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != vB.size()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector vX (this -> mnM);
        this -> _solve (vB, vX, dErr, NULL, NULL);
        return vX;
    }

    RVector solve (const RVector& vB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (vB, dErr);
    }

    BaseRMatrix solve (const BaseRMatrix& mB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != mB.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        BaseRMatrix mX (mB.msize(), mB.nsize());
        this -> _solve (mB, mX, dErr, NULL, NULL);
        return mX;
    }

    BaseRMatrix solve (const BaseRMatrix& mB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (mB, dErr);
    }

    RVector solve_lu (const basic_srmatrix& mLU, const int* pPivots, const RVector& vB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != vB.size() || this -> mnM != mLU.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector vX (this -> mnM);
        this -> _solve (vB, vX, dErr, mLU, pPivots);
        return vX;
    }

    RVector solve_lu (const basic_srmatrix& mLU, const int* pPivots, const RVector& vB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mLU, pPivots, vB, dErr);
    }

    BaseRMatrix solve_lu (const basic_srmatrix& mLU, const int* pPivots, const BaseRMatrix& mB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != mB.msize() || this -> mnM != mLU.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        BaseRMatrix mX (mB.msize(), mB.nsize());
        this -> _solve (mB, mX, dErr, mLU, pPivots);
        return mX;
    }

    BaseRMatrix solve_lu (const basic_srmatrix& mLU, const int* pPivots, const BaseRMatrix& mB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mLU, pPivots, mB, dErr);
    }

    // matrix determinant
    TR det() const throw (cvmexception)
    {
        return this -> _det();
    }

    // low-up factorization
    basic_srmatrix& low_up (const basic_srmatrix& m, int* nPivots) throw (cvmexception)
    {
        (*this) = m;
        this -> _low_up (nPivots);
        return *this;
    }

    basic_srmatrix low_up (int* nPivots) const throw (cvmexception)
    {
        basic_srmatrix mRes (*this);
        mRes._low_up (nPivots);
        return mRes;
    }

    // reciprocal of the condition number
    TR cond() const throw (cvmexception)
    {
        TR dCondNum(0.);
        __cond_num<TR, basic_srmatrix>(*this, dCondNum);        // universal method, no need to virtualize
        return dCondNum;
    }

    // matrix inversion
    basic_srmatrix& inv (const basic_srmatrix& mArg) throw (cvmexception)
    {
        __inv<basic_srmatrix>(*this, mArg);                     // overridden in srsmatrix, no need to virtualize
        return *this;
    }

    basic_srmatrix inv() const throw (cvmexception)
    {
        basic_srmatrix mRes (this -> mnM);
        __inv<basic_srmatrix>(mRes, *this);                     // overridden in srsmatrix, no need to virtualize
        return mRes;
    }

    // matrix exponent with given tolerance
    basic_srmatrix& exp (const basic_srmatrix& mArg, TR tol = basic_cvmMachSp<TR>()) throw (cvmexception)
    {
        __exp<basic_srmatrix, TR>(*this, mArg, tol);            // uses universal code inside - no need to virtualize
        return *this;
    }

    basic_srmatrix exp (TR tol = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        basic_srmatrix mRes (this -> mnM);
        __exp<basic_srmatrix, TR>(mRes, *this, tol);
        return mRes;
    }

    // this = v(1)*I + v(2)*m + v(3)*m^2 + ... + v(N)*m^(N-1)
    basic_srmatrix& polynom (const basic_srmatrix& m, const RVector& v) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector v1;
        if (v.incr() > 1) v1 << v;   // to make sure incr = 1
        __polynom<TR, RVector> (this -> mpD, this -> mnLD, this -> mnM, m._pd(), m._ldm(), v.incr() > 1 ? v1 : v);
        return *this;
    }

    // returns v(1)*I + v(2)*this + v(3)*this^2 + ... + v(N)*this^(N-1)
    basic_srmatrix polynom (const RVector& v) const
    {
        basic_srmatrix mRes (this -> mnM);
        RVector v1;
        if (v.incr() > 1) v1 << v;   // to make sure incr = 1
        __polynom<TR, RVector> (mRes.mpD, mRes.mnLD, this -> mnM, this -> mpD, this -> mnLD, v.incr() > 1 ? v1 : v);
        return mRes;
    }

   
    // Cholesky factorization
    basic_srmatrix& cholesky (const basic_srsmatrix<TR>& m) throw (cvmexception)
    {
        this -> _check_cholesky();  // doesn't work for band matrices
        *this = m;
        int nOutInfo = __cholesky<basic_srmatrix> (*this);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDEFINITE, nOutInfo);
        this -> _clean_low_triangle ();
        return *this;
    }

    // Bunch-Kaufman factorization
    basic_srmatrix& bunch_kaufman (const basic_srsmatrix<TR>& m, int* nPivots) throw (cvmexception)
    {
        this -> _check_bunch_kaufman();  // doesn't work for band matrices
        *this = m;
        __bunch_kaufman<basic_srmatrix> (*this, nPivots);
        return *this;
    }

    // QR - always full mode for square matrices
    void qr (basic_srmatrix<TR>& mQ, basic_srmatrix<TR>& mR) const throw (cvmexception)
    {
        this -> _qr(mQ, mR);
    }

    basic_srmatrix& identity()
    {
        this -> _vanish();
        this -> _plus_plus();
        return *this;
    }

    basic_srmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    virtual void _solve (const RVector& vB, RVector& vX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        vX = vB;
        RVector vB1;
        RVector vX1;
        if (vB.incr() > 1) vB1 << vB;   // to make sure incr = 1
        if (vX.incr() > 1) vX1 << vX;
        __solve<TR, TR, basic_srmatrix> (*this, 1, vB.incr() > 1 ? vB1 : vB, vB.size(), vX.incr() > 1 ? vX1 : vX, vX.size(), dErr, pLU, pPivots);
        if (vX.incr() > 1) vX = vX1;
    }

    virtual void _solve (const BaseRMatrix& mB, BaseRMatrix& mX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        mX = mB;
        __solve<TR, TR, basic_srmatrix> (*this, mB.nsize(), mB, mB.ld(), mX, mX.ld(), dErr, pLU, pPivots);
    }

protected:
    // protected constructors for inherited stuff
    basic_srmatrix (int nMN, int nLD, bool bZeroMemory)
        : BaseRMatrix (nMN, nMN, nLD, bZeroMemory)
    {
    }

    basic_srmatrix (TR* pD, int nMN, int nLD, int nSize)
        : BaseRMatrix (pD, nMN, nMN, nLD, nSize)
    {
    }

    virtual int _size   () const {return this -> mnSize;}
    virtual int _msize  () const {return this -> mnM;}
    virtual int _nsize  () const {return this -> mnN;}
    virtual int _ld     () const {return this -> mnLD;}
    virtual const TR* _p() const {return this -> mpD;}
    virtual TR*       _p()       {return this -> mpD;}

    // returns diagonal which IS l-value (shares memory)
    // 0 - main, negative - low, positive - up
    virtual RVector _diag (int nDiag) throw (cvmexception)
    {
        const int nD = abs (nDiag);
        if (nD >= this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nDiag);
        return RVector (this -> mpD + (nDiag > 0 ? nDiag * this -> mnLD : nD), this -> mnM - nD, this -> mnLD + 1);
    }

    // returns diagonal which IS NOT l-value (creates a copy)
    // 0 - main, negative - low, positive - up
    virtual const RVector _diag (int nDiag) const throw (cvmexception)
    {
        RVector vRet (const_cast<basic_srmatrix*>(this) -> _diag(nDiag));
        return vRet;
    }

    // returns main diagonal of low_up factorization
    virtual RVector _low_up_diag (basic_array<int>& naPivots) const throw (cvmexception)
    {
        return this -> low_up (naPivots).diag(0);
    }

    virtual void _transp()
    {
        this -> _sq_transp();
    }

    virtual void _plus_plus()
    {
        this -> _sq_plus_plus();
    }

    virtual void _minus_minus()
    {
        this -> _sq_minus_minus();
    }

    virtual TR _det() const throw (cvmexception)
    {
        TR dDet(0.);
        switch (this -> mnM)
        {
            case 0:
                break;
            case 1:
                dDet = this -> _ij_val (0, 0);
                break;
            case 2:
                dDet = this -> _ij_val (0, 0) * this -> _ij_val (1, 1) - 
                       this -> _ij_val (1, 0) * this -> _ij_val (0, 1);
                break;
            default:
                try
                {
                    static const TR one(1.);
                    basic_array<int> naPivots (this -> mnM);
                    RVector vUpDiag = this -> _low_up_diag (naPivots);

                    dDet = one;
                    for (int i = 1; i <= this -> mnM; ++i)
                    {
                        dDet *= vUpDiag[i];
                        if (i != naPivots[i]) dDet = - dDet;
                    }
                }
                catch (const cvmexception& e)
                {
                    if (e.cause() != CVM_SINGULARMATRIX) throw e;
                }
                break;
        }
        return dDet;
    }

    virtual void _low_up (int* nPivots) throw (cvmexception)
    {
        __low_up<basic_srmatrix>(*this, nPivots);
    }

    // QR - "economy" mode here
    virtual void _qr(basic_srmatrix<TR>& pmQ, basic_srmatrix<TR>& pmR) const throw (cvmexception)
    {
        if (this -> mnM != pmQ.msize() || this -> mnM != pmR.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        __qre<basic_rmatrix<TR>, basic_srmatrix<TR> > (*this, pmQ, pmR);
    }

    virtual void _check_cholesky() {}
    virtual void _check_bunch_kaufman() {}
};

// square band matrix
// way of storage:
/*
L=1, U=2          L=2, U=3        L=0, U=0      L=1, U=0

-                 -               d             d
--                --               d            *d
d**               ---               d            *d
*d**              d***               d            *d
 *d**             *d***               d            *d
  *d**            **d***               d            *d
   *d**            **d**                d            _
    *d*             **d*
     *d              **d
      -               --
                       -
*/
template <typename TR, typename TC>
class BandMatrix
{
protected:
    int mnKL;                                                   // number of non-zero sub-diagonals
    int mnKU;                                                   // number of non-zero super-diagonals

    BandMatrix()
        : mnKL(0), mnKU(0)
    {
    }

    BandMatrix (int nKL, int nKU)
        : mnKL (nKL), mnKU (nKU)
    {
    }

    virtual ~BandMatrix()
    {
    }

    virtual int       _size  () const = 0;
    virtual int       _msize () const = 0;
    virtual int       _nsize () const = 0;
    virtual int       _ld    () const = 0;
    virtual const TC* _p     () const = 0;
    virtual       TC* _p     () = 0;
    virtual void      _set_p (TC*) = 0;
    virtual void      _set   (TC* pD, int nSize, int nM, int nN, int nIncr, int nLD) = 0;

    void _bresize (int nNewM) throw (cvmexception)
    {
        TC* mpD = this -> _p();
        const int mnM = this -> _msize();
        const int mnSize = this -> _size();
        const int mnIncr = 1;
        if (nNewM != mnM)
        {
            if (nNewM < 0) throw cvmexception (CVM_WRONGSIZE, nNewM);
            const int nNewSize = nNewM * (1 + mnKL + mnKU);

            if (nNewSize == 0)
            {
                cvmFree<TC>(mpD);
                this -> _set (NULL, 0, 0, 0, 0, 0);
                mnKL = mnKU = 0;
            }
            else
            {
                TC* pD = cvmMalloc<TC>(nNewSize);
                if (nNewSize > mnSize) CleanMemory<TC> (pD, nNewSize);
                const int nMinSize = _cvm_min<int>(mnSize, nNewSize);
                __copy<TC> (nMinSize, mpD, mnIncr, pD, 1);
                CVM_ASSERT(pD, nNewSize * sizeof(TC))
                cvmFree<TC>(mpD);
                this -> _set (pD, nNewSize, nNewM, nNewM, 1, 1 + mnKL + mnKU);
            }
        }
    }

    void _check_dimensions () const throw (cvmexception)
    {
        if (mnKL < 0) throw cvmexception (CVM_WRONGSIZE, mnKL);
        if (mnKU < 0) throw cvmexception (CVM_WRONGSIZE, mnKU);
    }

    void _check_dimensions (const BandMatrix& m) const throw (cvmexception)
    {
        if (this -> _msize() != m._msize() || mnKL != m.mnKL || mnKU != m.mnKU) throw cvmexception (CVM_SIZESMISMATCH);
    }

    bool _bcontinuous () const
    {
        return this -> _msize() == 0 || this -> _ld() == 1 + mnKL + mnKU;
    }

    void _mbassign (const Matrix<TR,TC>& m)                               // m is a band matrix here
    {
        TC* mpD = this -> _p();
        const int mnSize = this -> _size();
        if (mpD != m.get())
        {
            __copy<TC> (mnSize, m, m.incr(), mpD, 1);
        }
    }

    TR _bnorm1() const                                          // 1-norm
    {
        int i, j;
        int nLen   = 0;
        int nShift = 0;
        TR  rNorm  = TR(0.);
        TR  rSum;
        const int mnN = this -> _nsize();
        basic_array<TC> rv (this -> _msize());

        for (j = 0; j < mnN; ++j)
        {
            rSum = TR(0.);
            this -> _get_col (j, rv, 1, &nLen, &nShift);

            nLen += nShift;
            for (i = nShift + 1; i <= nLen; ++i)
            {
                rSum += _abs (rv[i]);
            }

            if (rSum > rNorm)
            {
                rNorm = rSum;
            }
        }
        return rNorm;
    }

    TR _bnorminf() const                                        // infinity-norm
    {
        int i, j;
        int nLen = 0;
        int nShift = 0;
        TR  rNorm = TR(0.);
        TR  rSum;
        const int mnM = this -> _msize();
        basic_array<TC> rv (this -> _nsize());

        for (i = 0; i < mnM; ++i)
        {
            rSum = TR(0.);
            this -> _get_row (i, rv, 1, &nLen, &nShift);

            nLen += nShift;
            for (j = nShift + 1; j <= nLen; ++j)
            {
                rSum += _abs (rv[j]);
            }

            if (rSum > rNorm)
            {
                rNorm = rSum;
            }
        }
        return rNorm;
    }

    type_proxy<TC,TR> _b_ij_proxy_val (int i, int j)            // zero based
    {
        static const TC zero = TC(0.);
        TC* mpD = this -> _p();
        const int nA = j - mnKU;
        CVM_ASSERT(mpD, (i + j * (1 + mnKL + mnKU) - nA + 1) * sizeof(TC))
        return (nA < 0 || i >= nA) && i <= mnKL + j ? type_proxy<TC,TR>(mpD [i + j * (1 + mnKL + mnKU) - nA], false) : 
                                                      type_proxy<TC,TR>(zero, true);
    }

    TC _b_ij_val (int i, int j) const                           // zero based
    {
        static const TC zero = TC(0.);
        const TC* mpD = this -> _p();
        const int nA = j - mnKU;
        CVM_ASSERT(mpD, (i + j * (1 + mnKL + mnKU) - nA + 1) * sizeof(TC))
        return (nA < 0 || i >= nA) && i <= mnKL + j ? mpD [i + j * (1 + mnKL + mnKU) - nA] : zero;
    }

    void _get_col (int i, TC* pCol, int nIncr, int* pnLen = NULL, int* pnShift = NULL) const
    {
        const TC* mpD = this -> _p();
        const int mnM = this -> _msize();
        const int mnN = this -> _nsize();
        const int nCol = 1 + mnKL + mnKU;
        int nS         = nCol;
        int nShiftSrc  = 0;
        int nShiftDest = 0;

        if (i > mnKU)
        {
            nShiftDest = i - mnKU;
        }
        else
        {
            nShiftSrc = mnKU - i;
            nS -= nShiftSrc;
        }

        if (mnM - i <= mnKL)
        {
            nS -= mnKL + 1 - (mnN - i);
        }

        __copy<TC> (nS,
                   mpD + i * nCol + nShiftSrc,
                   1,
                   pCol + nShiftDest,
                   nIncr);

        if (pnLen) *pnLen = nS;
        if (pnShift) *pnShift = nShiftDest;
    }

    void _get_row (int i, TC* pCol, int nIncr, int* pnLen = NULL, int* pnShift = NULL) const
    {
        const TC* mpD = this -> _p();
        const int mnM = this -> _msize();
        const int mnN = this -> _nsize();
        const int nCol = mnKL + mnKU;
        int nS         = mnN;
        int nShiftSrc  = i + mnKU;
        int nShiftDest = 0;

        if (i > mnKL)
        {
            nShiftDest = i - mnKL;
            nShiftSrc += nShiftDest * nCol;
            nS -= nShiftDest;
        }
        if (mnN - i > mnKU)
        {
            nS -= (mnM - i) - mnKU - 1;
        }

        __copy<TC> (nS,
                   mpD + nShiftSrc,
                   nCol,
                   pCol + nShiftDest,
                   nIncr);

        if (pnLen) *pnLen = nS;
        if (pnShift) *pnShift = nShiftDest;
    }

    void _btransp() throw (cvmexception)
    {
        TC* mpD = this -> _p();
        const int mnN = this -> _nsize();
        if (mnKL > 0 || mnKU > 0)
        {
            const int nLU  = mnKL + mnKU;
            const int nCol = 1 + nLU;
            int nS, nShiftSrc;
            TC* pL;
            TC* pR;
            TC* pD = cvmMalloc<TC>(nCol * mnN);

            for (int i = 0; i < mnN; ++i)
            {
                nS = nCol;
                nShiftSrc = 0;

                if (i < mnKU)
                {
                    nShiftSrc = mnKU - i;
                    nS -= nShiftSrc;
                }
                if (mnN - i <= mnKL)
                {
                    nS -= mnKL + 1 - (mnN - i);
                }

                pL = mpD + i * nCol + nShiftSrc;
                pR = pD + (i > mnKU ? (i - mnKU + 1) * nCol - 1 : mnKL + i);

                __copy<TC> (nS, pL, 1, pR, nLU);
            }

            cvmFree<TC>(mpD);
            this -> _set_p (pD);
            std::swap(mnKL, mnKU);
        }
    }

    void _b_plus_plus()
    {
        TC* mpD = this -> _p();
        static const TC one(1.);
        const int mnN = this -> _nsize();
        const int nNext = 1 + mnKL + mnKU;
        const int nSize = nNext * mnN;
        CVM_ASSERT(mpD, nSize * sizeof(TC))
        for (int i = mnKU; i < nSize; i += nNext)
        {
            mpD[i] += one;
        }
    }

    void _b_minus_minus()
    {
        TC* mpD = this -> _p();
        static const TC one(1.);
        const int mnN = this -> _nsize();
        const int nNext = 1 + mnKL + mnKU;
        const int nSize = nNext * mnN;
        CVM_ASSERT(mpD, nSize * sizeof(TC))
        for (int i = mnKU; i < nSize; i += nNext)
        {
            mpD[i] -= one;
        }
    }

    void _b_replace (const BandMatrix& m) throw (cvmexception)  // matrix replacement, no assignment
    {
        TC* mpD = this -> _p();
        cvmFree<TC>(mpD);
        int mnSize = m._size();
        mpD = cvmMalloc<TC>(mnSize);
        CVM_ASSERT(mpD, (mnSize * sizeof(TC)))
        mnKL = m.mnKL;
        mnKU = m.mnKU;
        this -> _set (mpD, mnSize, m._msize(), m._nsize(), 1, m._ld());
    }

    void _resize_lu (int nNewKL, int nNewKU) throw (cvmexception)
    {
        if (nNewKL != mnKL || nNewKU != mnKU)
        {
            if (mnKL < 0) throw cvmexception (CVM_WRONGSIZE, mnKL);
            if (mnKU < 0) throw cvmexception (CVM_WRONGSIZE, mnKU);
            TC* mpD = this -> _p();
            const int mnM = this -> _msize();
            const int mnN = this -> _nsize();
            const int nOldLD = 1 + mnKL + mnKU;
            const int nNewLD = 1 + nNewKL + nNewKU;
            const int nMinKL = _cvm_min<int>(mnKL, nNewKL);
            const int nMinKU = _cvm_min<int>(mnKU, nNewKU);
            const int nNewSize = mnN * (1 + nNewKL + nNewKU);
            TC* pD = cvmMalloc<TC>(nNewSize);
            CVM_ASSERT(pD, nNewSize * sizeof(TC))
            CleanMemory<TC> (pD, nNewSize);

            for (int i = - nMinKU; i <= nMinKL; ++i)
            {
                __copy<TC> (mnN, mpD + (mnKU + i), nOldLD, pD + (nNewKU + i), nNewLD);
            }

            cvmFree<TC>(mpD);
            mnKL = nNewKL;
            mnKU = nNewKU;
            this -> _set (pD, nNewSize, mnM, mnN, 1, nNewLD);
        }
    }

public:
    int lsize() const
    {
        return mnKL;
    }

    int usize() const
    {
        return mnKU;
    }

    const int* _pl() const
    {
        return &mnKL;
    }

    const int* _pu() const
    {
        return &mnKU;
    }
};


template <typename TR>
class basic_srbmatrix : public basic_srmatrix<TR>, public BandMatrix<TR,TR>
{
    typedef std::complex<TR> TC;
    typedef basic_rvector<TR> RVector;
 
    typedef Array<TR,TR> BaseArray;
    typedef Matrix<TR,TR> BaseMatrix;
    typedef basic_rmatrix<TR> BaseRMatrix;
    typedef basic_srmatrix<TR> BaseSRMatrix;
    typedef BandMatrix<TR,TR> BaseBandMatrix;

  

protected:
    mutable BaseSRMatrix mSM;

    void _bake_SM() const
    {
        mSM.resize(this -> mnM);
        mSM.vanish();
        _copy_b_matrix<TR, TR, BaseSRMatrix, basic_srbmatrix>(mSM, * const_cast<basic_srbmatrix*>(this), false);
    }

public:
    basic_srbmatrix()
    {
    }

    explicit basic_srbmatrix (int nMN)
        : BaseSRMatrix (nMN, 1, true), BaseBandMatrix (0, 0)                    // diagonal matrix
    {
    }

    basic_srbmatrix (int nMN, int nKL, int nKU)
        : BaseSRMatrix (nMN, 1 + nKL + nKU, true), BaseBandMatrix (nKL, nKU)
    {
        this -> _check_dimensions();
    }

    basic_srbmatrix (TR* pD, int nMN, int nKL, int nKU)
        : BaseSRMatrix (pD, nMN, 1 + nKL + nKU, nMN * (1 + nKL + nKU)), BaseBandMatrix (nKL, nKU)
    {
        this -> _check_dimensions();
    }

    basic_srbmatrix (const basic_srbmatrix& m)
        : BaseSRMatrix (m.mnM, 1 + m.mnKL + m.mnKU, false), BaseBandMatrix (m.mnKL, m.mnKU)
    {
        this -> _massign(m);
    }

    basic_srbmatrix (const BaseRMatrix& m, int nKL, int nKU)
        : BaseSRMatrix (m.msize(), 1 + nKL + nKU, false), BaseBandMatrix (nKL, nKU)
    {
        if (this -> mnM != this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _check_dimensions();
        _copy_b_matrix<TR, TR, BaseRMatrix, basic_srbmatrix> (const_cast<BaseRMatrix&>(m), *this, true);
    }

    // diagonal square matrix constructor
    explicit basic_srbmatrix (const RVector& v)
        : BaseSRMatrix (v.size(), 1, false), BaseBandMatrix (0, 0)
    {
        __copy<TR> (this -> mnM, v, v.incr(), this -> mpD, 1);
    }

    type_proxy<TR,TR> operator () (int nIm, int nIn) throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_proxy_val (nIm - 1, nIn - 1);
    }

    TR operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    // returns column which CAN NOT be l-value
    const RVector operator () (int nI) const throw (cvmexception)
    {
        if (nI <= 0 || nI > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _col (nI - 1);
    }

    // returns row which CAN NOT be l-value
    const RVector operator [] (int nI) const throw (cvmexception)
    {
        if (nI <= 0 || nI > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _row (nI - 1);
    }

    basic_srbmatrix& operator = (const basic_srbmatrix& m) throw (cvmexception)
    {
        this -> _check_dimensions(m);
        this -> _massign(m);
        return *this;
    }

    basic_srbmatrix& assign (const RVector& v)                                  // assigns vector
    {
        this -> _assign (v, v.incr());
        return *this;
    }

    basic_srbmatrix& assign (const TR* pD)                                      // assigns foregn array (nIncr = 1)
    {
        this -> _assign (pD, 1);
        return *this;
    }

    // fills the content
    basic_srbmatrix& set (TR d)
    {
        this -> _set (d);
        return *this;
    }

    basic_srbmatrix& resize (int nNewMN) throw (cvmexception)
    {
        this -> _resize (nNewMN, nNewMN);
        return *this;
    }

    basic_srbmatrix& resize_lu (int nNewKL, int nNewKU) throw (cvmexception)
    {
        this -> _resize_lu (nNewKL, nNewKU);
        return *this;
    }

    bool operator == (const basic_srbmatrix& m) const
    {
        return this -> mnM == m.mnM && this -> mnN == m.mnN && this -> mnKL == m.mnKL && this -> mnKU == m.mnKU && this -> _mequals (m);
    }

    bool operator != (const basic_srbmatrix& m) const
    {
        return !(this -> operator == (m));
    }

    basic_srbmatrix& operator << (const basic_srbmatrix& m) throw (cvmexception)
    {
        this -> _check_ld();                                    // submatrix replacement is obviously not possible
        this -> _b_replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_srbmatrix operator + (const basic_srbmatrix& m) const throw (cvmexception)
    {
        this -> _check_dimensions (m);
        basic_srbmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_srbmatrix operator - (const basic_srbmatrix& m) const throw (cvmexception)
    {
        this -> _check_dimensions (m);
        basic_srbmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    basic_srbmatrix& sum (const basic_srbmatrix& m1, const basic_srbmatrix& m2) throw (cvmexception)
    {
        this -> _check_dimensions (m1);
        this -> _check_dimensions (m2);
        this -> _msum (m1, m2);
        return *this;
    }

    basic_srbmatrix& diff (const basic_srbmatrix& m1, const basic_srbmatrix& m2) throw (cvmexception)
    {
        this -> _check_dimensions (m1);
        this -> _check_dimensions (m2);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_srbmatrix& operator += (const basic_srbmatrix& m) throw (cvmexception)
    {
        this -> _check_dimensions (m);
        this -> _mincr (m);
        return *this;
    }

    basic_srbmatrix& operator -= (const basic_srbmatrix& m) throw (cvmexception)
    {
        this -> _check_dimensions (m);
        this -> _mdecr (m);
        return *this;
    }

    basic_srbmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_srbmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    // plus identity, prefix
    basic_srbmatrix& operator ++ ()
    {
        this -> _plus_plus();
        return *this;
    }

    // plus identity, postfix
    basic_srbmatrix& operator ++ (int)
    {
        this -> _plus_plus();
        return *this;
    }

    // minus identity, prefix
    basic_srbmatrix& operator -- ()
    {
        this -> _minus_minus();
        return *this;
    }

    // minus identity, postfix
    basic_srbmatrix& operator -- (int)
    {
        this -> _minus_minus();
        return *this;
    }

    basic_srbmatrix operator * (TR dMult) const
    {
        basic_srbmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_srbmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_srbmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_srbmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_srbmatrix& operator /= (TR dDiv) throw (cvmexception)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_srbmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    // transposed Matrix
    basic_srbmatrix operator ~ () const throw (cvmexception)
    {
        basic_srbmatrix mRes (*this);
        return mRes.transpose();
    }

    // well, not the best possible algorithm, has to be optimized
    basic_srbmatrix& transpose (const basic_srbmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnKL != m.mnKU || this -> mnKU != m.mnKL) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srbmatrix mTmp(m);
        mTmp.transpose();
        __copy<TR> (this -> mnSize, mTmp, mTmp.incr(), this -> mpD, this -> mnIncr);
        return *this;
    }

    basic_srbmatrix& transpose() throw (cvmexception)
    {
        this -> _transp();
        return *this;
    }

    RVector operator * (const RVector& v) const throw (cvmexception)
    {
        if (this -> mnN != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector vRes (this -> mnM);
        this -> _multiply (vRes, v, false);
        return vRes;
    }

    // special exclusion since matrix product is not commutative
    BaseRMatrix operator * (const BaseRMatrix& m) const throw (cvmexception)
    {
        _bake_SM();
        return mSM * m;
    }

    BaseSRMatrix operator * (const BaseSRMatrix& m) const throw (cvmexception)
    {
        _bake_SM();
        return mSM * m;
    }

    basic_srbmatrix operator * (const basic_srbmatrix& m) const throw (cvmexception)
    {
        _bake_SM();
        m._bake_SM();
        return basic_srbmatrix (mSM * m.mSM, 
                                _cvm_min<int>(this -> mnM - 1, this -> mnKL + m.mnKL), 
                                _cvm_min<int>(this -> mnM - 1, this -> mnKU + m.mnKU));
    }

    // low-up factorization
    basic_srbmatrix& low_up (const basic_srbmatrix& m, int* nPivots) throw (cvmexception)
    {
        (*this) = m;
        this -> _low_up (nPivots);
        return *this;
    }

    basic_srbmatrix low_up (int* nPivots) const throw (cvmexception)
    {
        basic_srbmatrix mRes (*this);
        mRes._low_up (nPivots);
        return mRes;
    }

    basic_srbmatrix& identity()
    {
        this -> _vanish();
        this -> _plus_plus();
        return *this;
    }

    basic_srbmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    // ATTENTION!!! This is not a good idea to use the following function. It's provided for
    // overriding of basic_rmatrix<TR>::mult only
    // this = m1 * m2
    basic_srbmatrix& mult (const BaseRMatrix& m1, const BaseRMatrix& m2) throw (cvmexception)
    {
        this -> _mult (m1, m2);
        return *this;
    }

    virtual TR* _pd()                                           // redefinition of Array's function
    {
#ifdef CVM_DEBUG
        assert (false);     // it's abnormal to call this function, this is pointer to copy, not to an object. only const version is allowable
#endif
        _bake_SM();
        return mSM._pd();
    }

    virtual const TR* _pd() const                               // redefinition of Array's function
    {
        _bake_SM();
        return mSM._pd();
    }

    // Euclid norm - band matrices require special care because of tails
    virtual TR norm() const throw (cvmexception)
    {
        TR dNorm = TR(0.), d;
        int i;

        for (i = 0; i <= this -> mnKL; ++i)
        {
            const RVector& dV = const_cast<basic_srbmatrix*>(this) -> diag(-i);
            d = dV.norm();
            dNorm += d * d;
        }
        for (i = 1; i <= this -> mnKU; ++i)
        {
            const RVector& dV = const_cast<basic_srbmatrix*>(this) -> diag(i);
            d = dV.norm();
            dNorm += d * d;
        }

        return _sqrt(dNorm);
    }

    virtual TR norm1() const
    {
        return this -> _bnorm1();
    }

    virtual TR norminf() const
    {
        return this -> _bnorminf();
    }

    // singular values in decreasing order
    virtual void _svd (RVector& vRes, BaseSRMatrix* pmU, BaseSRMatrix* pmVH) const throw (cvmexception)
    {
        if (pmU != NULL && pmVH != NULL && (this -> mnM != pmU -> msize() || this -> mnN != pmVH -> msize())) throw cvmexception (CVM_SIZESMISMATCH);
        __svd<TR, basic_srbmatrix, BaseSRMatrix> (vRes, vRes.size(), vRes.incr(), *this, pmU, pmVH);
    }

    virtual void _pinv (BaseRMatrix& mX, TR threshold) const throw (cvmexception)
    {
        __pinv<TR, basic_srbmatrix, BaseRMatrix> (mX, *this, threshold);
    }

    virtual void _solve (const RVector& vB, RVector& vX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        vX = vB;
        RVector vB1;
        RVector vX1;
        if (vB.incr() > 1) vB1 << vB;       // to make sure incr = 1
        if (vX.incr() > 1) vX1 << vX;
        __solve<TR, TR, basic_srbmatrix> (*this, 1, vB.incr() > 1 ? vB1 : vB, vB.size(), vX.incr() > 1 ? vX1 : vX, vX.size(), dErr, pLU, pPivots);
        if (vX.incr() > 1) vX = vX1;
    }

    virtual void _solve (const BaseRMatrix& mB, BaseRMatrix& mX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        mX = mB;
        __solve<TR, TR, basic_srbmatrix> (*this, mB.nsize(), mB, mB.ld(), mX, mX.ld(), dErr, pLU, pPivots);
    }

    // ?gbmv routines perform a matrix-vector operation defined as
    // vRes = alpha*m*v + beta * vRes or vRes = alpha*v'*m + beta * vRes
    // not virtual since __gbmv calls all virtual methods inside
    void _gbmv (bool bLeft, TR dAlpha, const RVector& v, TR dBeta, RVector& vRes) const
    {
        RVector vTmp;
        basic_srbmatrix mTmp;
        const TR* pDv = v;
        if (vRes.get() == pDv) vTmp << v;
        if (vRes.get() == this -> mpD) mTmp << *this;
        __gbmv<TR, basic_srbmatrix, RVector> (bLeft, vRes.get() == this -> mpD ? mTmp : *this, dAlpha, 
                                                     vRes.get() == pDv ? vTmp : v, dBeta, vRes);
    }

    virtual void _check_submatrix ()  throw (cvmexception) {throw cvmexception (CVM_SUBMATRIXNOTAVAILABLE, "srbmatrix");}

protected:
    virtual int       _size  () const { return this -> mnSize; }
    virtual int       _msize () const { return this -> mnM; }
    virtual int       _nsize () const { return this -> mnN; }
    virtual int       _ld    () const { return this -> mnLD; }
    virtual const TR* _p     () const { return this -> mpD; }
    virtual       TR* _p     ()       { return this -> mpD; }

    virtual void _set_p (TR* pD)
    {
        this -> mpD = pD;
    }

    virtual void _set (TR* pD, int nSize, int nM, int nN, int nIncr, int nLD)
    {
        this -> mpD = pD;
        this -> mnSize = nSize;
        this -> mnM = nM;
        this -> mnN = nN;
        this -> mnIncr = nIncr;
        this -> mnLD = nLD;
    }

    // for _msum _mdiff etc.
    virtual const TR* _p (const BaseMatrix& m) const
    {
        return m.get();
    }

    virtual int _ldm() const
    {
        return this -> mnM;
    }

    virtual const int* _pldm() const
    {
        return &this -> mnM;
    }

    // 0-based
    virtual RVector _row (int m)
    {
        RVector vRet (this -> mnM);
        this -> _get_row (m, vRet, vRet.incr());
        return vRet;
    }

    // 0-based
    virtual RVector _col (int n)
    {
        RVector vRet (this -> mnM);
        this -> _get_col (n, vRet, vRet.incr());
        return vRet;
    }

    // 0-based
    virtual const RVector _row (int m) const
    {
        RVector vRet (this -> mnM);
        this -> _get_row (m, vRet, vRet.incr());
        return vRet;
    }

    // 0-based
    virtual const RVector _col (int n) const
    {
        RVector vRet (this -> mnM);
        this -> _get_col (n, vRet, vRet.incr());
        return vRet;
    }

    virtual RVector _diag (int nDiag) throw (cvmexception)
    {
        const int nD = abs (nDiag);
        if (nDiag < 0 && nD > this -> mnKL || nDiag > 0 && nDiag > this -> mnKU) throw cvmexception (CVM_OUTOFRANGE, nDiag);
        const int nLU = this -> mnKL + this -> mnKU;
        return RVector (this -> mpD + this -> mnKU + (nDiag < 0 ? nD : nD * nLU), this -> mnM - nD, nLU + 1);
    }

    virtual const RVector _diag (int nDiag) const throw (cvmexception)
    {
        RVector vRet (const_cast<basic_srbmatrix*>(this) -> _diag(nDiag));
        return vRet;
    }

    virtual void _swap_rows(int, int)   throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_rows", "srbmatrix");}
    virtual void _swap_cols(int, int)   throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_cols", "srbmatrix");}
    virtual void _check_ger()           throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "ger", "srbmatrix");}
    virtual void _check_rank1update()   throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "rank1update", "srbmatrix");}
    virtual void _check_gemm()          throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "gemm", "srbmatrix");}
    virtual void _check_symm()          throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "symm", "srbmatrix");}
    virtual void _check_cholesky()      throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "cholesky", "srbmatrix");}
    virtual void _check_bunch_kaufman() throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "bunch_kaufman", "srbmatrix");}

    virtual void _resize (int nNewM, int) throw (cvmexception)
    {
        this -> _bresize(nNewM);
    }

    virtual bool _continuous () const
    {
        return this -> _bcontinuous ();
    }

    virtual void _massign (const BaseMatrix& m)
    {
        this -> _mbassign(m);
    }

    // zero based
    virtual type_proxy<TR,TR> _ij_proxy_val (int i, int j)
    {
        return this -> _b_ij_proxy_val (i, j);
    }

    // zero based
    virtual TR _ij_val (int i, int j) const
    {
        return this -> _b_ij_val (i, j);
    }

    virtual void _transp() throw (cvmexception)
    {
        this -> _btransp();
    }

    virtual void _plus_plus()
    {
        this -> _b_plus_plus();
    }

    virtual void _minus_minus()
    {
        this -> _b_minus_minus();
    }

    virtual int _indofmax() const
    {
        this -> _check_ld();
        _bake_SM();
        return mSM.indofmax();
    }

    virtual int _indofmin() const
    {
        this -> _check_ld();
        _bake_SM();
        return mSM.indofmin();
    }

    // returns main diagonal of low_up factorization
    virtual RVector _low_up_diag (basic_array<int>& naPivots) const throw (cvmexception)
    {
        return this -> low_up (naPivots).diag(0);
    }

    virtual void _scal (TR d)
    {
        __scal<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, d);                // zero tails are supposed here
    }

    virtual void _mult (const BaseRMatrix& m1, const BaseRMatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.msize() || this -> mnN != m2.nsize() || m1.nsize() != m2.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        BaseSRMatrix mR (this -> mnM);
        mR.mult (m1, m2);
        this -> _resize_lu (this -> mnM - 1, this -> mnM - 1);
        _copy_b_matrix<TR, TR, BaseSRMatrix, basic_srbmatrix> (const_cast<BaseSRMatrix&>(mR), *this, true);
    }

    virtual void _multiply (RVector& vRes, const RVector& v, bool bLeft) const
    {
        static const TR zero(0.);
        static const TR one(1.);
        this -> _gbmv (bLeft, one, v, zero, vRes);
    }

    virtual void _low_up (int* nPivots) throw (cvmexception)
    {
        __low_up<basic_srbmatrix>(*this, nPivots);
    }

    virtual int _ld_for_replace () const
    {
        return this -> mnM;
    }

    virtual int _size_for_replace () const
    {
        return this -> mnM * this -> mnN;
    }
};

template <typename TR>
class basic_srsmatrix : public basic_srmatrix<TR>
{
    typedef basic_rvector<TR> RVector;
    typedef Array<TR,TR> BaseArray;
    typedef Matrix<TR,TR> BaseMatrix;
    typedef SqMatrix<TR, TR> BaseSqMatrix;
    typedef basic_rmatrix<TR> BaseRMatrix;
    typedef basic_srmatrix<TR> BaseSRMatrix;

public:
    basic_srsmatrix()
    {
    }

    explicit basic_srsmatrix (int nMN)
        : BaseSRMatrix (nMN)
    {
    }

    basic_srsmatrix (TR* pD, int nMN, TR tol = basic_cvmMachSp<TR>())
        : BaseSRMatrix (pD, nMN, nMN, nMN * nMN)
    {
        this -> _check_symmetric(tol);
    }

    basic_srsmatrix (const basic_srsmatrix& m)
        : BaseSRMatrix (m.mnM, m.mnM, false)
    {
        this -> _massign(m);
    }

    explicit basic_srsmatrix (const BaseRMatrix& m, TR tol = basic_cvmMachSp<TR>())
        : BaseSRMatrix (m.msize(), m.nsize(), false)
    {
        if (this -> mnM != this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        this -> _check_symmetric(tol);
    }

    // diagonal square symmetric matrix constructor
    explicit basic_srsmatrix (const RVector& v)
        : BaseSRMatrix (v)
    {
    }

    // submatrix constructor
    basic_srsmatrix (basic_srsmatrix& m, int nRowCol, int nSize)        // 1-based
        : BaseSRMatrix (m._sub_pointer (nRowCol, nRowCol, nSize, nSize), nSize, m.ld(), nSize * nSize)
    {
    }

    TR operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    // returns column which CAN NOT be l-value
    const RVector operator () (int nFI) const throw (cvmexception)
    {
        return this -> _col (nFI - 1);
    }

    // returns row which CAN NOT be l-value
    const RVector operator [] (int nFI) const throw (cvmexception)
    {
        return this -> _row (nFI - 1);
    }

    // returns diagonal which IS NOT l-value (since it could break symmetricity)
    // 0 - main, negative - low, positive - up
    const RVector diag (int nDiag) const throw (cvmexception)
    {
        RVector vRet (this -> _diag(nDiag));
        return vRet;
    }

    basic_srsmatrix& operator = (const basic_srsmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        return *this;
    }

    basic_srsmatrix& assign (const RVector& v, TR tol = basic_cvmMachSp<TR>()) throw (cvmexception)             // assigns a vector
    {
        this -> _assign (v, v.incr());
        this -> _check_symmetric(tol);
        return *this;
    }

    basic_srsmatrix& assign (const TR* pD, TR tol = basic_cvmMachSp<TR>()) throw (cvmexception)                 // assigns a foregn array (nIncr = 1)
    {
        this -> _assign (pD, 1);
        this -> _check_symmetric(tol);
        return *this;
    }

    basic_srsmatrix& assign (int nRowCol, const basic_srsmatrix& m) throw (cvmexception)         // subvector assignment
    {
        if (nRowCol <= 0 || nRowCol > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nRowCol);
        if (m.mnM + nRowCol - 1 > this -> mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> _sub_pointer_nocheck (nRowCol, nRowCol), m._pd(), m.mnM, m.mnN, m.mnLD);
        return *this;
    }

    basic_srsmatrix& set (TR d)
    {
        this -> _set (d);
        return *this;
    }

    // sets both elements to keep matrix symmetric
    basic_srsmatrix& set (int nRow, int nCol, TR d)
    {
        this -> _set_at (nRow - 1, nCol - 1, d);
        return *this;
    }

    // sets both diagonals (keeps matrix symmetric)
    basic_srsmatrix& set_diag (int nDiag, const RVector& vDiag) throw (cvmexception)
    {
        RVector (this -> _diag(nDiag)) = vDiag;
        if (nDiag != 0)
        {
            RVector (this -> _diag(-nDiag)) = vDiag;
        }
        return *this;
    }

    basic_srsmatrix& resize (int nNewMN) throw (cvmexception)
    {
        this -> _resize (nNewMN, nNewMN);
        return *this;
    }

    bool operator == (const basic_srsmatrix& m) const
    {
        return this -> mnM == m.mnM && this -> mnN == m.mnN && this -> _mequals (m);
    }

    bool operator != (const basic_srsmatrix& m) const
    {
        return !(this -> operator == (m));
    }

    basic_srsmatrix& operator << (const basic_srsmatrix& m) throw (cvmexception)
    {
        this -> _replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_srsmatrix operator + (const basic_srsmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srsmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_srsmatrix operator - (const basic_srsmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srsmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    basic_srsmatrix& sum (const basic_srsmatrix& m1, const basic_srsmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _msum (m1, m2);
        return *this;
    }

    basic_srsmatrix& diff (const basic_srsmatrix& m1, const basic_srsmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_srsmatrix& operator += (const basic_srsmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mincr (m);
        return *this;
    }

    basic_srsmatrix& operator -= (const basic_srsmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdecr (m);
        return *this;
    }

    basic_srsmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_srsmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    // plus identity, prefix
    basic_srsmatrix& operator ++ ()
    {
        this -> _plus_plus();
        return *this;
    }

    // plus identity, postfix
    basic_srsmatrix& operator ++ (int)
    {
        this -> _plus_plus();
        return *this;
    }

    // minus identity, prefix
    basic_srsmatrix& operator -- ()
    {
        this -> _minus_minus();
        return *this;
    }

    // minus identity, postfix
    basic_srsmatrix& operator -- (int)
    {
        this -> _minus_minus();
        return *this;
    }

    basic_srsmatrix operator * (TR dMult) const
    {
        basic_srsmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_srsmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_srsmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_srsmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_srsmatrix& operator /= (TR dDiv) throw (cvmexception)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_srsmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    // transposed Matrix - does nothing
    basic_srsmatrix operator ~ () const
    {
        return *this;
    }

    basic_srsmatrix& transpose (const basic_srsmatrix& m) throw (cvmexception)
    {
        (*this) = m;
        return *this;
    }

    basic_srsmatrix& transpose()
    {
        return *this;
    }

    RVector operator * (const RVector& v) const throw (cvmexception)
    {
        if (this -> mnN != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector vRes (this -> mnM);
        this -> _multiply (vRes, v, false);
        return vRes;
    }

    // special exclusion since matrix product is not commutative
    BaseRMatrix operator * (const BaseRMatrix& m) const throw (cvmexception)
    {
        if (this -> mnN != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        BaseRMatrix mRes(m.msize(), m.nsize());
        mRes._symm (true, *this, m, one, zero);
        return mRes;
    }

    BaseSRMatrix operator * (const BaseSRMatrix& m) const throw (cvmexception)
    {
        if (this -> mnN != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        BaseSRMatrix mRes(this -> mnM);
        mRes._symm (true, *this, m, one, zero);
        return mRes;
    }

    // this = alpha*v*v' + beta*this
    basic_srsmatrix& syrk (TR alpha, const RVector& v, TR beta) throw (cvmexception)
    {
        if (this -> mnM != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __syrk<TR, basic_srsmatrix> (false, alpha, 1, v, v.size(), beta, *this);
        this -> _flip();
        return *this;
    }

    // this = alpha*a*a' + beta*this or this = alpha*a'*a + beta*this
    basic_srsmatrix& syrk (bool bTransp, TR alpha, const BaseRMatrix& mA, TR beta) throw (cvmexception)
    {
        if (this -> mnM != (bTransp ? mA.nsize() : mA.msize())) throw cvmexception (CVM_SIZESMISMATCH);
        __syrk<TR, basic_srsmatrix> (bTransp, alpha, bTransp ? mA.msize() : mA.nsize(), mA, mA.ld(), beta, *this);
        this -> _flip();
        return *this;
    }

    // this = alpha*v1*v2' + alpha*v2*v1' + beta*this
    basic_srsmatrix& syr2k (TR alpha, const RVector& v1, const RVector& v2, TR beta) throw (cvmexception)
    {
        if (this -> mnM != v1.size() || this -> mnM != v2.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __syr2k<TR, basic_srsmatrix> (false, alpha, 1, v1, v1.size(), v2, v2.size(), beta, *this);
        this -> _flip();
        return *this;
    }

    // this = alpha*a*b' + alpha*b*a' + beta*this or this = alpha*a'*b + alpha*b'*a + beta*this
    basic_srsmatrix& syr2k (bool bTransp, TR alpha, const BaseRMatrix& mA, const BaseRMatrix& mB, TR beta) throw (cvmexception)
    {
        if (this -> mnM != (bTransp ? mA.nsize() : mA.msize()) ||
            this -> mnM != (bTransp ? mB.nsize() : mB.msize()) ||
            bTransp ? mA.msize() != mB.msize() : mA.nsize() != mB.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        __syr2k<TR, basic_srsmatrix> (bTransp, alpha, bTransp ? mA.msize() : mA.nsize(), mA, mA.ld(), mB, mB.ld(), beta, *this);
        this -> _flip();
        return *this;
    }

    // matrix inversion
    basic_srsmatrix& inv (const basic_srsmatrix& mArg) throw (cvmexception)
    {
        __inv<basic_srsmatrix>(*this, mArg);
        return *this;
    }

    basic_srsmatrix inv() const throw (cvmexception)
    {
        basic_srsmatrix mRes (this -> mnM);
        __inv<basic_srsmatrix>(mRes, *this);
        return mRes;
    }

    // matrix exponent with given tolerance
    basic_srsmatrix& exp (const basic_srsmatrix& mArg, TR tol = basic_cvmMachSp<TR>()) throw (cvmexception)
    {
        __exp_symm<basic_srsmatrix, TR> (*this, mArg, tol);
        return *this;
    }

    basic_srsmatrix exp (TR tol = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        basic_srsmatrix msRes (this -> mnM);
        __exp_symm<basic_srsmatrix, TR> (msRes, *this, tol);
        return msRes;
    }

    // this = v(1)*I + v(2)*m + v(3)*m^2 + ... + v(N)*m^(N-1)
    basic_srsmatrix& polynom (const basic_srsmatrix& m, const RVector& v) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector v1;
        if (v.incr() > 1) v1 << v;   // to make sure incr = 1
        __polynom<TR, RVector> (this -> mpD, this -> mnLD, this -> mnM, m._pd(), m._ldm(), v.incr() > 1 ? v1 : v);
        return *this;
    }

    // returns v(1)*I + v(2)*this + v(3)*this^2 + ... + v(N)*this^(N-1)
    basic_srsmatrix polynom (const RVector& v) const
    {
        basic_srsmatrix msRes (this -> mnM);
        RVector v1;
        if (v.incr() > 1) v1 << v;   // to make sure incr = 1
        __polynom<TR, RVector> (msRes, msRes.mnLD, this -> mnM, this -> mpD, this -> mnLD, v.incr() > 1 ? v1 : v);
        return msRes;
    }

    // eigenvalues
    // we don't use _eig here since this is the special case - symmetric matrix
    RVector eig (BaseSRMatrix& mEigVect) const throw (cvmexception)
    {
        RVector vEig(this -> mnM);
        __eig<RVector, basic_srsmatrix, BaseSRMatrix> (vEig, *this, &mEigVect, true);
        return vEig;
    }

    RVector eig() const throw (cvmexception)
    {
        RVector vEig(this -> mnM);
        __eig<RVector, basic_srsmatrix, BaseSRMatrix> (vEig, *this, NULL, true);
        return vEig;
    }

    // Cholesky factorization
    BaseSRMatrix cholesky () const throw (cvmexception)
    {
        BaseSRMatrix mRes (*this);
        int nOutInfo = __cholesky<BaseSRMatrix> (mRes);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDEFINITE, nOutInfo);
        mRes._clean_low_triangle ();
        return mRes;
    }

    // Bunch-Kaufman factorization
    BaseSRMatrix bunch_kaufman (int* nPivots) const throw (cvmexception)
    {
        BaseSRMatrix mRes (*this);
        __bunch_kaufman<BaseSRMatrix> (mRes, nPivots);
        return mRes;
    }

    basic_srsmatrix& identity()
    {
        this -> _vanish();
        this -> _plus_plus();
        return *this;
    }

    basic_srsmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    // matrix equilibration (useful for further solve and solve_lu calling)
    // returns true if equilibration was needed and performed
    bool equilibrate (RVector& vScalings, RVector& vB) throw (cvmexception)
    {
        if (this -> mnM != vB.size()) throw cvmexception (CVM_SIZESMISMATCH);
        bool bRes = this -> _equilibrate (vScalings);
        if (bRes)
        {
            for (int i = 1; i <= this -> mnM; ++i)
            {
                vB[i] *= vScalings[i];
            }
        }
        return bRes;
    }

    bool equilibrate (RVector& vScalings, BaseRMatrix& mB) throw (cvmexception)
    {
        if (this -> mnM != mB.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        bool bRes = this -> _equilibrate (vScalings);
        if (bRes)
        {
            for (int j = 1; j <= mB.nsize(); ++j)
            {
                for (int i = 1; i <= this -> mnM; ++i)
                {
                    mB(i,j) *= vScalings[i];
                }
            }
        }
        return bRes;
    }

    // special care for symmetric matrices
    basic_srsmatrix& _factorize (const basic_srsmatrix& m, int* nPivots, bool& bPositiveDefinite) throw (cvmexception)
    {
        (*this) = m;
        int nOutInfo = __cholesky<BaseSRMatrix>(*this);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0)
        {
            (*this) = m;
            __bunch_kaufman<BaseSRMatrix>(*this, nPivots);
            bPositiveDefinite = false;
        }
        else {
            bPositiveDefinite = true;
        }
        return *this;
    }

    virtual TR* _pd()                                           // redefinition of Array's function
    {
        return this -> mpD;
    }

    virtual const TR* _pd() const                               // redefinition of Array's function
    {
        return this -> mpD;
    }

    virtual TR norminf() const                                  // infinity-norm - the same as 1-norm for symmetric matrices
    {
        return this -> norm1();
    }

    // makes lower triangle to be equal to upper one
    virtual void _flip()
    {
        if (this -> mnM > 1)
        {
            const int nM1 = this -> mnLD + 1, nM1m = this -> mnLD - 1, nM2m = this -> mnM - 1;
            int i = 1, j = 1, m;
            for (;;)
            {
                m = this -> mnM - i;
                __copy<TR> (m, this -> mpD + j + nM1m, this -> mnLD, this -> mpD + j, 1);
                if (i >= nM2m)
                {
                    break;
                }
                i++;
                j += nM1;
            }
        }
    }

    virtual void _check_submatrix ()  throw (cvmexception) {throw cvmexception (CVM_SUBMATRIXNOTAVAILABLE, "srsmatrix");}

protected:
    virtual const TR* _p (const BaseMatrix& m) const  // for _msum _mdiff etc.
    {
        return m.get();
    }

    virtual void _check_ger()         throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "ger", "srsmatrix");}
    virtual void _check_rank1update() throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "rank1update", "srsmatrix");}
    virtual void _check_gemm()        throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "gemm", "srsmatrix");}
    virtual void _swap_rows(int, int) throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_rows", "srsmatrix");}
    virtual void _swap_cols(int, int) throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_cols", "srsmatrix");}

    // we do nothing here - it's symmetric 
    virtual void _transp()
    {
    }

    // returns main diagonal of low_up factorization
    // this call is useless for symmetric matrices. this call would mean serious CVM internal error
    virtual RVector _low_up_diag (basic_array<int>&) const throw (cvmexception)
    {
        throw cvmexception (CVM_NOTIMPLEMENTED, "_low_up_diag");
    }

    virtual void _scal (TR d)
    {
        __scal<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, d);     // zero tails are supposed here
    }

    virtual void _multiply (RVector& vRes, const RVector& v, bool) const
    {
        RVector vTmp;
        basic_srsmatrix mTmp;
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        const TR* pDm = this -> mpD;
        const TR* pDv = v;
        if (vRes.get() == pDv) vTmp << v;
        if (vRes.get() == pDm) mTmp << *this;
        __symv<TR, basic_srsmatrix, RVector> (vRes.get() == pDm ? mTmp : *this, one, 
                                              vRes.get() == pDv ? vTmp : v, zero, vRes);
    }

    virtual void _solve (const RVector& vB, RVector& vX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        RVector vB1 (vB);    // to make sure incr = 1 and equilibrate
        RVector vScalings (this -> mnM);
        basic_srsmatrix m (*this);
        const bool bEquilibrated = m.equilibrate (vScalings, vB1);
        RVector vX1 (vB1);   // to make sure incr = 1
        __solve<TR, TR, basic_srsmatrix> (m, 1, vB1, vB1.size(), vX1, vX1.size(), dErr, pLU, pPivots);
        if (bEquilibrated)
        {
            for (int i = 1; i <= this -> mnM; ++i)
            {
                vX[i] = vX1[i] * vScalings[i];
            }
        }
        else
        {
            vX = vX1;
        }
    }

    virtual void _solve (const BaseRMatrix& mB, BaseRMatrix& mX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        BaseRMatrix mB1 (mB);    // to equilibrate
        RVector vScalings (this -> mnM);
        basic_srsmatrix m (*this);
        const bool bEquilibrated = m.equilibrate (vScalings, mB1);
        mX = mB1;
        __solve<TR, TR, basic_srsmatrix> (m, mB.nsize(), mB, mB.ld(), mX, mX.ld(), dErr, pLU, pPivots);
        if (bEquilibrated)
        {
            for (int j = 1; j <= mX.nsize(); ++j)
            {
                for (int i = 1; i <= this -> mnM; ++i)
                {
                    mX(i,j) *= vScalings[i];
                }
            }
        }
    }

    // matrix determinant
    virtual TR _det() const throw (cvmexception)
    {
        TR dDet = TR(0.);
        switch (this -> mnM)
        {
            case 0:
                break;
            case 1:
                dDet = this -> _ij_val (0, 0);
                break;
            case 2:
                dDet = this -> _ij_val (0, 0) * this -> _ij_val (1, 1) - 
                       this -> _ij_val (1, 0) * this -> _ij_val (0, 1);
                break;
            default:
                try
                {
                    static const TR one(1.);
                    bool bPositiveDefinite = false;
                    basic_srsmatrix m(this -> mnM);
                    basic_array<int> nPivots (this -> mnM);
                    m._factorize (*this, nPivots, bPositiveDefinite);

                    int i;
                    dDet = one;
                    if (bPositiveDefinite)
                    {
                        const RVector vUpDiag = m.diag(0);
                        for (i = 1; i <= this -> mnM; ++i)
                        {
                            dDet *= vUpDiag[i] * vUpDiag[i];    //here we use Cholesky factorization
                        }
                    }
                    else
                    {
                        const RVector vEig = this -> eig();
                        for (i = 1; i <= this -> mnM; ++i)
                        {
                            dDet *= vEig[i];                    //here we use eigenvalues. probably there is a better way.
                        }
                    }
                }
                catch (const cvmexception& e)
                {
                    if (e.cause() != CVM_WRONGBUNCHKAUFMANFACTOR) throw e;
                }
                break;
        }
        return dDet;
    }

    // matrix equilibration helper
    bool _equilibrate (RVector& vScalings) throw (cvmexception)
    {
        if (this -> mnM != vScalings.size()) throw cvmexception (CVM_SIZESMISMATCH);
        bool bRes = false;
        TR dCond(0.);
        TR dMax  = TR();
        static const TR sp = basic_cvmMachSp<TR>();
        static const TR sp_inv = TR(1.) / sp;

        __poequ<TR, basic_srsmatrix, RVector> (*this, vScalings, dCond, dMax);

        if (dCond < TR(0.1) || _abs(dMax) <= sp || _abs(dMax) >= sp_inv)
        {
            bRes = true;
            for (int i = 0; i < this -> mnM; ++i)
            {
                for (int j = i; j < this -> mnM; ++j)
                {
                    this -> mpD[this -> mnLD * j + i] *= vScalings[i+1] * vScalings[j+1];
                }
            }
        }
        return bRes;
    }

    // sets both elements to keep matrix symmetric, checks ranges
    void _set_at (int nRow, int nCol, TR val) throw (cvmexception)      // zero based
    {
        if (nRow < 0 || nRow >= this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nRow);
        if (nCol < 0 || nCol >= this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nCol);
        this -> mpD[this -> mnLD * nCol + nRow] = val;
        if (nRow != nCol)
        {
            this -> mpD[this -> mnLD * nRow + nCol] = val;
        }
    }

    void _check_symmetric (TR tol) const throw (cvmexception)
    {
        for (int j = 0; j < this -> mnN; ++j)
        {
            for (int i = 0; i < this -> mnM; ++i)
            {
                if (i != j && _abs (this -> mpD[this -> mnLD * j + i] - this -> mpD[this -> mnLD * i + j]) > tol)
                {
                    throw cvmexception (CVM_MATRIXNOTSYMMETRIC);
                }
            }
        }
    }
};

template <typename TR>
inline basic_rvector<TR> operator * (TR d, const basic_rvector<TR>& v)
{
    return v * d;
}
template <typename TR>
inline basic_rmatrix<TR> operator * (TR d, const basic_rmatrix<TR>& m)
{
    return m * d;
}
template <typename TR>
inline basic_srmatrix<TR> operator * (TR d, const basic_srmatrix<TR>& m)
{
    return m * d;
}
template <typename TR>
inline basic_srbmatrix<TR> operator * (TR d, const basic_srbmatrix<TR>& m)
{
    return m * d;
}
template <typename TR>
inline basic_srsmatrix<TR> operator * (TR d, const basic_srsmatrix<TR>& m)
{
    return m * d;
}

template <typename TR>
inline basic_rvector<TR>  operator * (CVM_LONGEST_INT d, const basic_rvector<TR>& v)
{
    return v * static_cast<TR>(d);
}
template <typename TR>
inline basic_rmatrix<TR>  operator * (CVM_LONGEST_INT d, const basic_rmatrix<TR>& m)
{
    return m * static_cast<TR>(d);
}
template <typename TR>
inline basic_srmatrix<TR>  operator * (CVM_LONGEST_INT d, const basic_srmatrix<TR>& m)
{
    return m * static_cast<TR>(d);
}
template <typename TR>
inline basic_srbmatrix<TR>  operator * (CVM_LONGEST_INT d, const basic_srbmatrix<TR>& m)
{
    return m * static_cast<TR>(d);
}
template <typename TR>
inline basic_srsmatrix<TR>  operator * (CVM_LONGEST_INT d, const basic_srsmatrix<TR>& m)
{
    return m * static_cast<TR>(d);
}

#if defined (CVM_FLOAT)
typedef float  treal;
#else
typedef double treal;
#endif

typedef std::complex<treal> tcomplex;

typedef basic_array    <int>             iarray;
typedef basic_rvector  <treal>           rvector;
typedef basic_rmatrix  <treal>           rmatrix;
typedef basic_srmatrix <treal>           srmatrix;
typedef basic_srbmatrix<treal>           srbmatrix;
typedef basic_srsmatrix<treal>           srsmatrix;

// identity matrices creation
template <typename TR>
inline const basic_srmatrix<TR> basic_eye_real (int nM)
{
    basic_srmatrix<TR> mI(nM);
    ++mI;
    return mI;
}

inline const srmatrix eye_real (int nM)
{
    return basic_eye_real<treal>(nM);
}

inline treal cvmMachMin()
{
    return basic_cvmMachMin<treal>();
}
inline treal cvmMachSp()
{
    return basic_cvmMachSp<treal>();
}


CVM_NAMESPACE_END


// BLAS callback error handler
#if !defined (_MSC_VER)
#define XERBLA xerbla_
#endif

extern "C" {
    void   __stdcall XERBLA    (const char* szSubName,
#ifdef CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES
                                const unsigned int nLen,
#endif
                                const int* pnParam) throw (cvm::cvmexception);
}

#if defined (_MSC_VER)
#   pragma warning(pop)
#endif

#endif                  // _CVM_H
