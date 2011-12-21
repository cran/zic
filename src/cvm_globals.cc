//                  CVM Class Library
//                  http://cvmlib.com
//
//          Copyright Sergei Nikolaev 1992-2008
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include "cvm.h"

#if defined (_MSC_VER)
#   pragma warning(disable:4786)
#endif

extern "C" {
    void __stdcall XERBLA (const char* szSubName,
    #ifdef CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES
                        const unsigned int,
    #endif
                        const int* pnParam) throw (cvm::cvmexception)
    {
        throw cvm::cvmexception (CVM_WRONGMKLARG2, *pnParam, szSubName);
    }
}

#if !defined (CVM_STATIC) && (defined (_MSC_VER) || defined (__WATCOMC__))
BOOL APIENTRY DllMain (HANDLE /*hModule*/,
                       DWORD  ul_reason_for_call,
                       LPVOID /*lpReserved*/)
{
    switch (ul_reason_for_call)
    {
        case DLL_PROCESS_ATTACH:
        case DLL_THREAD_ATTACH:
        case DLL_THREAD_DETACH:
        case DLL_PROCESS_DETACH:
            break;
    }
    return TRUE;
}
#endif

//#define CVM_USE_CRITICAL_SECTION_NOT_MUTEX

CVM_NAMESPACE_BEG

// multithreaded synchronizer
class CriticalSection {
#if defined (CVM_MT)
private:
    bool mbOK;

    #if defined (WIN32) || defined (_WIN32)
    #   if defined (CVM_USE_CRITICAL_SECTION_NOT_MUTEX)
            ::CRITICAL_SECTION mCriticalSection;
    #   else
            HANDLE mMutex;
    #   endif
    #else                                                       // POSIX Threads library assumed
        pthread_mutex_t mMutex;
        pthread_mutexattr_t mMutexAttr;
    #endif
#endif

public:
    CriticalSection() 
#if defined (CVM_MT)
        : mbOK (false)
#endif
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
    #   if defined (CVM_USE_CRITICAL_SECTION_NOT_MUTEX)
        if (!::InitializeCriticalSectionAndSpinCount (&mCriticalSection, 0x80000400))
        {
            ::InitializeCriticalSection (&mCriticalSection);
        }
        mbOK = true;
    #   else
            mMutex = ::CreateMutex (NULL, FALSE, NULL);
            mbOK = mMutex != NULL;
    #   endif
    #else
        if (pthread_mutexattr_init (&mMutexAttr) == 0 &&
            pthread_mutexattr_setpshared (&mMutexAttr, PTHREAD_PROCESS_PRIVATE) == 0 &&
            pthread_mutex_init (&mMutex, &mMutexAttr) == 0)
        {
            mbOK = true;
        }
    #endif
#endif
    }

    ~CriticalSection()
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (mbOK)
        {
    #   if defined (CVM_USE_CRITICAL_SECTION_NOT_MUTEX)
            ::DeleteCriticalSection (&mCriticalSection);
    #   else
            ::CloseHandle(mMutex);
    #   endif
        }
    #else
        pthread_mutexattr_destroy (&mMutexAttr);
        pthread_mutex_destroy (&mMutex);
    #endif
#endif
    }

    void enter()
#if defined (CVM_MT)
    throw (cvmexception)
#endif
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (!mbOK)
        {
            throw cvmexception (CVM_SEMAPHOREERROR);
        }
    #   if defined (CVM_USE_CRITICAL_SECTION_NOT_MUTEX)
            ::EnterCriticalSection (&mCriticalSection); 
    #   else
            ::WaitForSingleObject (mMutex, INFINITE);
    #   endif
    #else
        if (!mbOK || pthread_mutex_lock (&mMutex) != 0)
        {
            throw cvmexception (CVM_SEMAPHOREERROR);
        }
    #endif
#endif
    }

    void leave()
#if defined (CVM_MT)
    throw (cvmexception)
#endif
    {
#if defined (CVM_MT)
    #if defined (WIN32) || defined (_WIN32)
        if (!mbOK)
        {
            throw cvmexception (CVM_SEMAPHOREERROR);
        }
    #   if defined (CVM_USE_CRITICAL_SECTION_NOT_MUTEX)
            ::LeaveCriticalSection (&mCriticalSection);
    #   else
            ::ReleaseMutex(mMutex);
    #   endif
    #else
        if (!mbOK || pthread_mutex_unlock (&mMutex) != 0)
        {
            throw cvmexception (CVM_SEMAPHOREERROR);
        }
    #endif
#endif
    }
};


// these species must be global to avoid problems in multithreading environments
CriticalSection gCS;
MemoryPool gPool;

class Lock
{
public:
    Lock ();
    ~Lock ();
};

Lock::Lock ()
{
    gCS.enter();
}
Lock::~Lock ()
{
    gCS.leave();
}


// 5.5.2 - moved out of cvm.h
cvmexception::cvmexception (int nCause, ...)
    : mnCause (nCause)
{
    va_list argList;
    va_start (argList, nCause);
#if defined (CVM_VSNPRINTF_S_DEFINED)
    const int nLength = CVM_VSNPRINTF (mszMsg, sizeof(mszMsg), sizeof(mszMsg) - 1, _get_message(mnCause), argList);
#else
    const int nLength = CVM_VSNPRINTF (mszMsg, sizeof(mszMsg) - 1, _get_message(mnCause), argList);
#endif
    va_end (argList);
    if (nLength >= (int) sizeof(mszMsg))
    {
        mszMsg[sizeof(mszMsg) - 1] = '\0';
    }
}

cvmexception::cvmexception (const cvmexception& e)
    : std::exception(e), mnCause (e.mnCause)
{
#if defined (CVM_STRCPY_S_DEFINED)
    strcpy_s (mszMsg, sizeof(mszMsg), e.mszMsg);
#else
    strcpy (mszMsg, e.mszMsg);
#endif
}


#ifdef CVM_USE_POOL_MANAGER
CVM_API void _cvm_assert (const void* pvBlock, size_t nBytes)
{
    Lock l;
    gPool.Assert (pvBlock, nBytes);
}
#else
CVM_API void _cvm_assert (const void*, size_t)
{
}
#endif  // CVM_USE_POOL_MANAGER

CVM_API tbyte* _cvmMalloc (size_t nBytes) throw (cvmexception)
{
    Lock lock;
    return gPool.Malloc (nBytes);
}

CVM_API tbyte* _cvmAddRef (const tbyte* pD)
{
    Lock lock;
    return gPool.AddRef (pD);
}

CVM_API int _cvmFree (tbyte*& pD)
{
    Lock lock;
    return gPool.FFFree (pD);
}

CVM_API void cvmExit()
{
#ifdef CVM_USE_POOL_MANAGER
    gPool.Clear();
#endif
}

#ifdef CVM_USE_POOL_MANAGER
#define CVM_PAGE_SIZE (0x1000)
#define CVM_HEAP_SIZE ((int) 0x40000000)

size_t _up_value (size_t n)                                     // the least power of 2 multiplied by 2
{
    if (n < CVM_PAGE_SIZE)                                      // let small objects be in one page
    {
        n = CVM_PAGE_SIZE;
    }
    else //if (n < CVM_HEAP_SIZE)
    {
        int i = 0;
        while (n >> i) ++i;
        if (i && n & (1 << (i - 1)) - 1) ++i;                   // obey warning C4554 :)
        n = 1 << i;
    }
    return n;
}
#endif  // CVM_USE_POOL_MANAGER

MemoryPool::MemoryPool()
{
}

MemoryPool::~MemoryPool()
{
#ifdef CVM_USE_POOL_MANAGER
    Clear();
#endif
}

#ifdef CVM_USE_POOL_MANAGER
#ifdef __BORLANDC__
#    pragma warn -8091
#endif

void MemoryPool::Clear()
{
    std::for_each (mOutBlocks.rbegin(), mOutBlocks.rend(), MemoryPool::DeletePtr());
    mOutBlocks.clear();
}

#ifdef __BORLANDC__
#    pragma warn +8091
#endif

#if !defined (CVM_ALLOCATOR)
#   define CVM_ALLOCATOR std::allocator
#endif

template <class T>
inline CVM_ALLOCATOR<T>& AllocatorInstance()
{
    static CVM_ALLOCATOR<T> _A;
    return _A;
}

#endif  // CVM_USE_POOL_MANAGER


tbyte* MemoryPool::Malloc (size_t nBytes) throw (cvmexception)
{
#ifdef CVM_USE_POOL_MANAGER
    if (nBytes >= CVM_HEAP_SIZE) throw cvmexception (CVM_WRONGSIZE, nBytes);
    if (nBytes == 0) return static_cast<tbyte*>(NULL);

    tbyte* pB = mMemoryBlocks.GetFreeBlock (nBytes);

    if (pB == NULL)     // There is no suitable memory block. Let's create a new one.
    {
        const size_t nUpBytes = _up_value (nBytes);
        const size_t nRest    = nUpBytes - nBytes;

        try
        {
            pB = AllocatorInstance<tbyte>().allocate(nUpBytes, NULL);
        }
        catch (const std::bad_alloc&)
        {
        }
        if (pB == NULL)
        {
            throw (cvmexception (CVM_OUTOFMEMORY));
        }

        mOutBlocks.push_back (pB);
        mMemoryBlocks.AddPair (pB, nBytes, nRest);
    }
#else
    tbyte* pB = NULL;
    try
    {
//        pB = AllocatorInstance<tbyte>().allocate(nBytes, NULL);
        pB = new tbyte[nBytes];
    }
    catch (const std::bad_alloc&)
    {
    }
    if (pB == NULL)
    {
        throw (cvmexception (CVM_OUTOFMEMORY));
    }

    mMemoryBlocks.AddNew (pB, nBytes);
#endif
    return pB;
}

tbyte* MemoryPool::AddRef (const tbyte* pD)
{
    return mMemoryBlocks.AddRef (pD);
}

int MemoryPool::FFFree (tbyte*& pToFree) throw (cvmexception)
{
    int nRefCounter = mMemoryBlocks.FreeBlock (pToFree);
    if (!nRefCounter)
    {
        pToFree = NULL;
    }
    return nRefCounter;
}



#ifdef CVM_USE_POOL_MANAGER

void MemoryBlocks::AddBlock (tbyte* pBlock, size_t nBytes, bool bOccupied)
{
    if (!bOccupied)                                             // Add freed block
    {
        itr_FreeIt j;
        itr_Blocks i = mBlocks.upper_bound (pBlock);
        itr_Blocks i_next = i;
                                                                // Is there upper neighboring memory block?
        if (i != mBlocks.end())
        {
            tbyte* pUpperBlock = (*i).first;
            j = mFreeIt.find (pUpperBlock);
            if (j != mFreeIt.end() && pBlock + nBytes == pUpperBlock)           // Yes. It's free and will be concatenated
            {
                nBytes += (*i).second.mnSize;
                ++i_next;
                mBlocks.erase (i);
                i = i_next;
                mFreeBs.erase ((*j).second);
                mFreeIt.erase (j);
            }
        }
                                                                // Is there lower neighboring memory block?
        if (i != mBlocks.begin() && mBlocks.size() > 0)
        {
            --i;
            tbyte* pLowerBlock = (*i).first;
            const size_t nLowerBytes = (*i).second.mnSize;
            j = mFreeIt.find (pLowerBlock);
            if (j != mFreeIt.end() && pLowerBlock + nLowerBytes == pBlock)      // Yes. It's free and will be concatenated
            {
                pBlock = pLowerBlock;
                nBytes += nLowerBytes;
                mBlocks.erase (i);
                mFreeBs.erase ((*j).second);
                mFreeIt.erase (j);
            }
        }
        mFreeIt[pBlock] = mFreeBs.insert (std::pair<int, tbyte*>(nBytes, pBlock));
    }

    mBlocks.insert (std::pair<tbyte*, BlockProperty>(pBlock, BlockProperty(nBytes, 1)));
}

int MemoryBlocks::FreeBlock (tbyte* pBlock)
{
    int nRefCounter = 0;
    itr_Blocks i = mBlocks.find (pBlock);

    if (i != mBlocks.end())
    {
        if (mFreeIt.find (pBlock) == mFreeIt.end())
        {
            nRefCounter = -- (*i).second.mnRefCount;
            if (nRefCounter <= 0)
            {
                const int nBytes = (*i).second.mnSize;
                mBlocks.erase (i);
                AddBlock (pBlock, nBytes, false);                               // return free block to the pool
            }
        }
#ifdef CVM_DEBUG
        else
        {
            assert (mFreeIt.find (pBlock) == mFreeIt.end());
        }
#endif
    }
    else
    {
        nRefCounter = -1;                                                       // foreign array.
    }

    return nRefCounter;
}

tbyte* MemoryBlocks::GetFreeBlock (size_t nBytes)
{
    tbyte* pBlock = NULL;
    if (mFreeBs.size() > 0)
    {
        // Is there a suitable memory block?
        itr_FreeBs i = mFreeBs.lower_bound (nBytes);

        // Yes. Let's use it.
        if (i != mFreeBs.end())
        {
            const size_t nRest = (*i).first - nBytes;
            pBlock = (*i).second;

            mFreeBs.erase (i);
            if (mFreeIt.size() > 0 && mFreeIt.find(pBlock) != mFreeIt.end())
            {
                mFreeIt.erase (pBlock);
            }
            if (mBlocks.size() > 0 && mBlocks.find(pBlock) != mBlocks.end())
            {
                mBlocks.erase (pBlock);
            }

            AddPair (pBlock, nBytes, nRest);
        }
    }
    return pBlock;
}

tbyte* MemoryBlocks::AddRef (const tbyte* pcBlock)
{
    tbyte* pBlock = const_cast<tbyte*>(pcBlock);
    itr_Blocks i = mBlocks.find (pBlock);
    if (i != mBlocks.end())
    {
        ++ (*i).second.mnRefCount;
    }
    else
    {
        // This is a foreign array. Leave it alone.
    }
    return pBlock;
}


#ifdef CVM_DEBUG
void MemoryBlocks::Assert (const void* pvBlock, size_t nBytes)
{
    tbyte* pBlock = (tbyte*) const_cast<void*>(pvBlock);
    itr_Blocks i = mBlocks.find (pBlock);
    if (i != mBlocks.end())
    {
        const size_t nSize = (*i).second.mnSize;
        assert (nSize >= nBytes);
    }
    else
    {
        tbyte* pB;
        size_t nB;
        itr_Blocks end = mBlocks.end();
        for (i = mBlocks.begin(); i != end; ++i)
        {
            pB = (*i).first;
            nB = (*i).second.mnSize;
            if (pBlock >= pB && pBlock < pB + nB)
            {
                tbyte* pBase = pB + nB;
                tbyte* pTest = pBlock + nBytes;
                assert (pTest <= pBase);
            }
        }
    }
}
#endif

void MemoryBlocks::AddPair (tbyte* pBlock, size_t nBytes, size_t nRest)
{
    AddBlock (pBlock, nBytes, true);                            // occupied block...
    if (nRest > 0)
    {
        AddBlock (pBlock + nBytes, nRest, false);               // ...and the rest free block
    }
}

#else // CVM_USE_POOL_MANAGER

tbyte* MemoryBlocks::AddRef (const tbyte* pcBlock)
{
    CVM_PTR_WRAPPER nBlock = reinterpret_cast<CVM_PTR_WRAPPER>(pcBlock);
    itr_Blocks i = mBlocks.find (nBlock);
    if (i != mBlocks.end())
    {
        ++(*i).second.mnRefCount;
    }
    return reinterpret_cast<tbyte*>(nBlock);
}

int MemoryBlocks::FreeBlock (tbyte* pBlock)
{
    int nRefCounter;
    itr_Blocks i = mBlocks.find (reinterpret_cast<CVM_PTR_WRAPPER>(pBlock));
    if (i != mBlocks.end())
    {
        nRefCounter = -- (*i).second.mnRefCount;
        if (nRefCounter <= 0)
        {
//            AllocatorInstance<tbyte>().deallocate(pBlock, (*i).second.mnSize);
            delete[] pBlock;
            mBlocks.erase (i);
        }
    }
    else
    {
        nRefCounter = -1;
    }
    return nRefCounter;
}

void MemoryBlocks::AddNew (tbyte* pBlock, size_t nBytes)
{
    mBlocks.insert (std::pair<CVM_PTR_WRAPPER, BlockProperty>
                   (reinterpret_cast<CVM_PTR_WRAPPER>(pBlock), BlockProperty(nBytes, 1)));
}

#endif // !CVM_USE_POOL_MANAGER

CVM_NAMESPACE_END
