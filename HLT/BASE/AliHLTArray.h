/****************************************************************************
 * This file is property of and copyright by the ALICE HLT Project          *
 * ALICE Experiment at CERN, All rights reserved.                           *
 *                                                                          *
 * Copyright (C) 2009 Matthias Kretz <kretz@kde.org>                        *
 *               for The ALICE HLT Project.                                 *
 *                                                                          *
 * Permission to use, copy, modify and distribute this software and its     *
 * documentation strictly for non-commercial purposes is hereby granted     *
 * without fee, provided that the above copyright notice appears in all     *
 * copies and that both the copyright notice and this permission notice     *
 * appear in the supporting documentation. The authors make no claims       *
 * about the suitability of this software for any purpose. It is            *
 * provided "as is" without express or implied warranty.                    *
 ***************************************************************************/

/**
 * \file AliHLTArray.h
 * \author Matthias Kretz <kretz@kde.org>
 *
 * This file contains the classes AliHLTResizableArray and AliHLTFixedArray with AliHLTArray as base
 * class. It's a drop-in replacement for C-Arrays. It makes it easy to use variable sized arrays on
 * the stack and pass arrays as arguments to other functions with an optional bounds-checking
 * enabled for the whole time.
 */

#ifndef ALIHLTARRAY_H
#define ALIHLTARRAY_H

#ifndef assert
#include <assert.h>
#endif

#if defined(__MMX__) || defined(__SSE__)
#include <mm_malloc.h>
#else
#include <cstdlib>
#endif

namespace AliHLTArrayInternal
{
  template<bool> class STATIC_ASSERT_FAILURE;
  template<> class STATIC_ASSERT_FAILURE<true> {};
}

#define ALIHLTARRAY_STATIC_ASSERT_CONCAT_HELPER(a, b) a##b
#define ALIHLTARRAY_STATIC_ASSERT_CONCAT(a, b) ALIHLTARRAY_STATIC_ASSERT_CONCAT_HELPER(a, b)
#define ALIHLTARRAY_STATIC_ASSERT(cond, msg) \
  typedef AliHLTArrayInternal::STATIC_ASSERT_FAILURE<cond> ALIHLTARRAY_STATIC_ASSERT_CONCAT(_STATIC_ASSERTION_FAILED_##msg, __LINE__); \
  ALIHLTARRAY_STATIC_ASSERT_CONCAT(_STATIC_ASSERTION_FAILED_##msg, __LINE__) Error_##msg; \
  (void) Error_##msg

template<typename T, int Dim> class AliHLTArray;

namespace AliHLTInternal
{
  // XXX
  // The ArrayBoundsCheck and Allocator classes implement a virtual destructor only in order to
  // silence the -Weff-c++ warning. It really is not required for these classes to have a virtual
  // dtor since polymorphism is not used (AliHLTResizableArray and AliHLTFixedArray are allocated on
  // the stack only). The virtual dtor only adds an unnecessary vtable to the code.
#ifndef ENABLE_ARRAY_BOUNDS_CHECKING
  /**
   * no-op implementation that for no-bounds-checking
   */
  class ArrayBoundsCheck
  {
    protected:
      virtual inline ~ArrayBoundsCheck() {}
      inline bool IsInBounds( int ) const { return true; }
      inline void SetBounds( int, int ) {}
      inline void MoveBounds( int ) {}
  };
#define BOUNDS_CHECK(x, y)
#else
  /**
   * implementation for bounds-checking.
   */
  class ArrayBoundsCheck
  {
    protected:
      virtual inline ~ArrayBoundsCheck() {}
      /**
       * checks whether the given offset is valid
       */
      inline bool IsInBounds( int x ) const;
      /**
       * set the start and end offsets that are still valid
       */
      inline void SetBounds( int start, int end ) { fStart = start; fEnd = end; }
      /**
       * move the start and end offsets by the same amount
       */
      inline void MoveBounds( int d ) { fStart += d; fEnd += d; }

    private:
      int fStart;
      int fEnd;
  };
#define BOUNDS_CHECK(x, y) if (AliHLTInternal::ArrayBoundsCheck::IsInBounds(x)) {} else return y
#endif
  template<typename T, int alignment> class Allocator
  {
    protected:
      virtual inline ~Allocator() {}
#if defined(__MMX__) || defined(__SSE__)
      static inline T *Alloc( int s ) { T *p = reinterpret_cast<T *>( _mm_malloc( s * sizeof( T ), alignment ) ); return new( p ) T[s]; }
      static inline void Free( const T *const p ) { /** p->~T(); */ _mm_free( p ); } // XXX: doesn't call dtor because it's an array
#else
      static inline T *Alloc( int s ) { T *p; posix_memalign( &p, alignment, s * sizeof( T ) ); return new( p ) T[s]; }
      static inline void Free( const T *const p ) { std::free( p ); } // XXX: doesn't call dtor because it's an array
#endif
  };
  template<typename T> class Allocator<T, 0>
  {
    protected:
      virtual inline ~Allocator() {}
      static inline T *Alloc( int s ) { return new T[s]; }
      static inline void Free( const T *const p ) { delete[] p; }
  };
  /**
   * Array base class for dimension dependent behavior
   */
  template<typename T, int Dim> class ArrayBase;

  /**
   * 1-dim arrays only have operator[]
   */
  template<typename T>
  class ArrayBase<T, 1> : public ArrayBoundsCheck
  {
      friend class ArrayBase<T, 2>;
    public:
      /**
       * return a reference to the value at the given index
       */
      inline T &operator[]( int x ) { BOUNDS_CHECK( x, fData[0] ); return fData[x]; }
      /**
       * return a const reference to the value at the given index
       */
      inline const T &operator[]( int x ) const { BOUNDS_CHECK( x, fData[0] ); return fData[x]; }

    protected:
      T *fData;
      inline void SetSize( int, int, int ) {}
  };

  /**
   * 2-dim arrays should use operator(int, int)
   * operator[] can be used to return a 1-dim array
   */
  template<typename T>
  class ArrayBase<T, 2> : public ArrayBoundsCheck
  {
      friend class ArrayBase<T, 3>;
    public:
      /**
       * return a reference to the value at the given indexes
       */
      inline T &operator()( int x, int y ) { BOUNDS_CHECK( x * fStride + y, fData[0] ); return fData[x * fStride + y]; }
      /**
       * return a const reference to the value at the given indexes
       */
      inline const T &operator()( int x, int y ) const { BOUNDS_CHECK( x * fStride + y, fData[0] ); return fData[x * fStride + y]; }
      /**
       * return a 1-dim array at the given index. This makes it behave like a 2-dim C-Array.
       */
      inline AliHLTArray<T, 1> operator[]( int x );
      /**
       * return a const 1-dim array at the given index. This makes it behave like a 2-dim C-Array.
       */
      inline const AliHLTArray<T, 1> operator[]( int x ) const;

    protected:
      T *fData;
      int fStride;
      inline void SetSize( int, int y, int ) { fStride = y; }
  };

  /**
   * 3-dim arrays should use operator(int, int, int)
   * operator[] can be used to return a 2-dim array
   */
  template<typename T>
  class ArrayBase<T, 3> : public ArrayBoundsCheck
  {
    public:
      /**
       * return a reference to the value at the given indexes
       */
      inline T &operator()( int x, int y, int z );
      /**
       * return a const reference to the value at the given indexes
       */
      inline const T &operator()( int x, int y, int z ) const;
      /**
       * return a 2-dim array at the given index. This makes it behave like a 3-dim C-Array.
       */
      inline AliHLTArray<T, 2> operator[]( int x );
      /**
       * return a const 2-dim array at the given index. This makes it behave like a 3-dim C-Array.
       */
      inline const AliHLTArray<T, 2> operator[]( int x ) const;

    protected:
      T *fData;
      int fStrideX;
      int fStrideY;
      inline void SetSize( int, int y, int z ) { fStrideX = y * z; fStrideY = z; }
  };

  // XXX AlignedData really is an internal struct, but the RuleChecker doesn't understand that
  template<typename T, unsigned int Size, int alignment> class AlignedData;
  template<typename T, unsigned int Size> class AlignedData<T, Size, 0>
  {
    protected:
      T d[Size];
  };
#ifdef __GNUC__
#define ALIGN(n) __attribute__((aligned(n)))
#else
#define ALIGN(n) __declspec(align(n))
#endif
  template<typename T, unsigned int Size> class AlignedData<T, Size, 4>
  {
    protected:
      ALIGN( 4 ) T d[Size];
  };
  template<typename T, unsigned int Size> class AlignedData<T, Size, 8>
  {
    protected:
      ALIGN( 8 ) T d[Size];
  };
  template<typename T, unsigned int Size> class AlignedData<T, Size, 16>
  {
    protected:
      ALIGN( 16 ) T d[Size];
  };
  template<typename T, unsigned int Size> class AlignedData<T, Size, 32>
  {
    protected:
      ALIGN( 32 ) T d[Size];
  };
  template<typename T, unsigned int Size> class AlignedData<T, Size, 64>
  {
    protected:
      ALIGN( 64 ) T d[Size];
  };
  template<typename T, unsigned int Size> class AlignedData<T, Size, 128>
  {
    protected:
      ALIGN( 128 ) T d[Size];
  };
#undef ALIGN
} // namespace AliHLTInternal

/**
 * C-Array like class with the dimension dependent behavior defined in the ArrayBase class
 */
template < typename T, int Dim = 1 >
class AliHLTArray : public AliHLTInternal::ArrayBase<T, Dim>
{
  public:
    typedef AliHLTInternal::ArrayBase<T, Dim> Parent;
    /**
     * allows you to check for validity of the array by casting to bool
     */
    inline operator bool() const { return Parent::fData != 0; }
    /**
     * allows you to check for validity of the array
     */
    inline bool IsValid() const { return Parent::fData != 0; }

    /**
     * returns a reference to the data at index 0
     */
    inline T &operator*() { BOUNDS_CHECK( 0, Parent::fData[0] ); return *Parent::fData; }
    /**
     * returns a const reference to the data at index 0
     */
    inline const T &operator*() const { BOUNDS_CHECK( 0, Parent::fData[0] ); return *Parent::fData; }

    /**
     * returns a pointer to the data
     * This circumvents bounds checking so it should not be used.
     */
    inline T *Data() { return Parent::fData; }
    /**
     * returns a const pointer to the data
     * This circumvents bounds checking so it should not be used.
     */
    inline const T *Data() const { return Parent::fData; }

    /**
     * moves the array base pointer so that the data that was once at index 0 will then be at index -x
     */
    inline AliHLTArray operator+( int x ) const;
    /**
     * moves the array base pointer so that the data that was once at index 0 will then be at index x
     */
    inline AliHLTArray operator-( int x ) const;
};

/**
 * Owns the data. When it goes out of scope the data is freed.
 *
 * The memory is allocated on the heap.
 *
 * Instantiate this class on the stack. Allocation on the heap is disallowed.
 *
 * \param T type of the entries in the array.
 * \param Dim selects the operator[]/operator() behavior it should have. I.e. makes it behave like a
 * 1-, 2- or 3-dim array. (defaults to 1)
 * \param alignment Defaults to 0 (default alignment). Other valid values are any multiples of 2.
 *                  This is especially useful for aligning data for SIMD vectors.
 *
 * \warning when using alignment the type T may not have a destructor (well it may, but it won't be
 * called)
 *
 * Example:
 * \code
 * void init( AliHLTArray<int> a, int size )
 * {
 *   for ( int i = 0; i < size; ++i ) {
 *     a[i] = i;
 *   }
 * }
 *
 * int size = ...;
 * AliHLTResizableArray<int> foo( size ); // notice that size doesn't have to be a constant like it
 *                                        // has to be for C-Arrays in ISO C++
 * init( foo, size );
 * // now foo[i] == i
 *
 * \endcode
 */
template < typename T, int Dim = 1, int alignment = 0 >
class AliHLTResizableArray : public AliHLTArray<T, Dim>, public AliHLTInternal::Allocator<T, alignment>
{
  public:
    typedef AliHLTInternal::ArrayBase<T, Dim> Parent;
    /**
     * does not allocate any memory
     */
    inline AliHLTResizableArray();
    /**
     * use for 1-dim arrays: allocates x * sizeof(T) bytes for the array
     */
    inline AliHLTResizableArray( int x );
    /**
     * use for 2-dim arrays: allocates x * y * sizeof(T) bytes for the array
     */
    inline AliHLTResizableArray( int x, int y );
    /**
     * use for 3-dim arrays: allocates x * y * z * sizeof(T) bytes for the array
     */
    inline AliHLTResizableArray( int x, int y, int z );

    /**
     * frees the data
     */
    inline ~AliHLTResizableArray() { AliHLTInternal::Allocator<T, alignment>::Free( Parent::fData ); }

    /**
     * use for 1-dim arrays: resizes the memory for the array to x * sizeof(T) bytes.
     *
     * \warning this does not keep your previous data. If you were looking for this you probably
     * want to use std::vector instead.
     */
    inline void Resize( int x );
    /**
     * use for 2-dim arrays: resizes the memory for the array to x * y * sizeof(T) bytes.
     *
     * \warning this does not keep your previous data. If you were looking for this you probably
     * want to use std::vector instead.
     */
    inline void Resize( int x, int y );
    /**
     * use for 3-dim arrays: resizes the memory for the array to x * y * z * sizeof(T) bytes.
     *
     * \warning this does not keep your previous data. If you were looking for this you probably
     * want to use std::vector instead.
     */
    inline void Resize( int x, int y, int z );

  private:
    // disable allocation on the heap
    void *operator new( size_t );

    // disable copy
    AliHLTResizableArray( const AliHLTResizableArray & );
    AliHLTResizableArray &operator=( const AliHLTResizableArray & );
};

/**
 * Owns the data. When it goes out of scope the data is freed.
 *
 * The memory is allocated on the stack.
 *
 * Instantiate this class on the stack.
 *
 * \param T type of the entries in the array.
 * \param Size number of entries in the array.
 * \param Dim selects the operator[]/operator() behavior it should have. I.e. makes it behave like a
 * 1-, 2- or 3-dim array. (defaults to 1)
 */
template < typename T, unsigned int Size, int Dim = 1, int alignment = 0 >
class AliHLTFixedArray : public AliHLTArray<T, Dim>
{
  public:
    typedef AliHLTInternal::ArrayBase<T, Dim> Parent;
    inline AliHLTFixedArray() { Parent::fData = &fDataOnStack.d[0]; Parent::SetBounds( 0, Size - 1 ); }

  private:
    // disable allocation on the heap
    void *operator new( size_t );

    AliHLTInternal::AlignedData<T, Size, alignment> fDataOnStack;

    // disable copy
    AliHLTFixedArray( const AliHLTFixedArray & );
    AliHLTFixedArray &operator=( const AliHLTFixedArray & );
};





////////////////////////
//// implementation ////
////////////////////////




namespace AliHLTInternal
{
#ifdef ENABLE_ARRAY_BOUNDS_CHECKING
  inline bool ArrayBoundsCheck::IsInBounds( int x ) const
  {
    assert( x >= fStart );
    assert( x <= fEnd );
    return ( x >= fStart && x <= fEnd );
  }
#endif

  template<typename T>
  inline AliHLTArray<T, 1> ArrayBase<T, 2>::operator[]( int x )
  {
    x *= fStride;
    typedef AliHLTArray<T, 1> AT1;
    BOUNDS_CHECK( x, AT1() );
    AliHLTArray<T, 1> a;
    a.fData = &fData[x];
    a.ArrayBoundsCheck::operator=( *this );
    a.MoveBounds( -x );
    return a;
  }

  template<typename T>
  inline const AliHLTArray<T, 1> ArrayBase<T, 2>::operator[]( int x ) const
  {
    x *= fStride;
    typedef AliHLTArray<T, 1> AT1;
    BOUNDS_CHECK( x, AT1() );
    AliHLTArray<T, 1> a;
    a.fData = &fData[x];
    a.ArrayBoundsCheck::operator=( *this );
    a.MoveBounds( -x );
    return a;
  }

  template<typename T>
  inline T &ArrayBase<T, 3>::operator()( int x, int y, int z )
  {
    BOUNDS_CHECK( x * fStrideX + y + fStrideY + z, fData[0] );
    return fData[x * fStrideX + y + fStrideY + z];
  }
  template<typename T>
  inline const T &ArrayBase<T, 3>::operator()( int x, int y, int z ) const
  {
    BOUNDS_CHECK( x * fStrideX + y + fStrideY + z, fData[0] );
    return fData[x * fStrideX + y + fStrideY + z];
  }
  template<typename T>
  inline AliHLTArray<T, 2> ArrayBase<T, 3>::operator[]( int x )
  {
    x *= fStrideX;
    typedef AliHLTArray<T, 2> AT2;
    BOUNDS_CHECK( x, AT2() );
    AliHLTArray<T, 2> a;
    a.fData = &fData[x];
    a.fStride = fStrideY;
    a.ArrayBoundsCheck::operator=( *this );
    a.MoveBounds( -x );
    return a;
  }
  template<typename T>
  inline const AliHLTArray<T, 2> ArrayBase<T, 3>::operator[]( int x ) const
  {
    x *= fStrideX;
    typedef AliHLTArray<T, 2> AT2;
    BOUNDS_CHECK( x, AT2() );
    AliHLTArray<T, 2> a;
    a.fData = &fData[x];
    a.fStride = fStrideY;
    a.ArrayBoundsCheck::operator=( *this );
    a.MoveBounds( -x );
    return a;
  }
} // namespace AliHLTInternal


template<typename T, int Dim>
inline AliHLTArray<T, Dim> AliHLTArray<T, Dim>::operator+( int x ) const
{
  AliHLTArray<T, Dim> r( *this );
  r.fData += x;
  r.MoveBounds( -x );
  return r;
}
template<typename T, int Dim>
inline AliHLTArray<T, Dim> AliHLTArray<T, Dim>::operator-( int x ) const
{
  AliHLTArray<T, Dim> r( *this );
  r.fData -= x;
  r.MoveBounds( x );
  return r;
}

template<typename T, int Dim, int alignment>
inline AliHLTResizableArray<T, Dim, alignment>::AliHLTResizableArray()
{
  Parent::fData = 0;
  Parent::SetBounds( 0, -1 );
}
template<typename T, int Dim, int alignment>
inline AliHLTResizableArray<T, Dim, alignment>::AliHLTResizableArray( int x )
{
  ALIHLTARRAY_STATIC_ASSERT( Dim == 1, AliHLTResizableArray1_used_with_incorrect_dimension );
  Parent::fData = AliHLTInternal::Allocator<T, alignment>::Alloc( x );
  Parent::SetBounds( 0, x - 1 );
}
template<typename T, int Dim, int alignment>
inline AliHLTResizableArray<T, Dim, alignment>::AliHLTResizableArray( int x, int y )
{
  ALIHLTARRAY_STATIC_ASSERT( Dim == 2, AliHLTResizableArray2_used_with_incorrect_dimension );
  Parent::fData = AliHLTInternal::Allocator<T, alignment>::Alloc( x * y );
  Parent::SetSize( x, y, 0 );
  Parent::SetBounds( 0, x * y - 1 );
}
template<typename T, int Dim, int alignment>
inline AliHLTResizableArray<T, Dim, alignment>::AliHLTResizableArray( int x, int y, int z )
{
  ALIHLTARRAY_STATIC_ASSERT( Dim == 3, AliHLTResizableArray3_used_with_incorrect_dimension );
  Parent::fData = AliHLTInternal::Allocator<T, alignment>::Alloc( x * y * z );
  Parent::SetSize( x, y, z );
  Parent::SetBounds( 0, x * y * z - 1 );
}
template<typename T, int Dim, int alignment>
inline void AliHLTResizableArray<T, Dim, alignment>::Resize( int x )
{
  ALIHLTARRAY_STATIC_ASSERT( Dim == 1, AliHLTResizableArray1_resize_used_with_incorrect_dimension );
  AliHLTInternal::Allocator<T, alignment>::Free( Parent::fData );
  Parent::fData = ( x == 0 ) ? 0 : AliHLTInternal::Allocator<T, alignment>::Alloc( x );
  Parent::SetBounds( 0, x - 1 );
}
template<typename T, int Dim, int alignment>
inline void AliHLTResizableArray<T, Dim, alignment>::Resize( int x, int y )
{
  ALIHLTARRAY_STATIC_ASSERT( Dim == 2, AliHLTResizableArray2_resize_used_with_incorrect_dimension );
  AliHLTInternal::Allocator<T, alignment>::Free( Parent::fData );
  Parent::fData = ( x == 0 ) ? 0 : AliHLTInternal::Allocator<T, alignment>::Alloc( x * y );
  Parent::SetSize( x, y, 0 );
  Parent::SetBounds( 0, x * y - 1 );
}
template<typename T, int Dim, int alignment>
inline void AliHLTResizableArray<T, Dim, alignment>::Resize( int x, int y, int z )
{
  ALIHLTARRAY_STATIC_ASSERT( Dim == 3, AliHLTResizableArray3_resize_used_with_incorrect_dimension );
  AliHLTInternal::Allocator<T, alignment>::Free( Parent::fData );
  Parent::fData = ( x == 0 ) ? 0 : AliHLTInternal::Allocator<T, alignment>::Alloc( x * y * z );
  Parent::SetSize( x, y, z );
  Parent::SetBounds( 0, x * y * z - 1 );
}

#undef BOUNDS_CHECK

#endif // ALIHLTARRAY_H
