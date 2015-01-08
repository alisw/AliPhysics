/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtMatrix.hh,v 1.1 2009-03-16 16:50:49 robbep Exp $
 *
 * Description:
 *   Class to make simple computations with matrices: assignment, product,
 *      determinant, inverse... Still many functions could be implemented.
 *
 * Modification history:
 *   Jordi Garra Ticó     2008/07/03         File created
 *****************************************************************************/


#ifndef __EVT_MATRIX_HH__
#define __EVT_MATRIX_HH__

#include <vector>
#include <sstream>


template <class T> class EvtMatrix
{
private:
  T** _mat;
  int _range;
public:
  EvtMatrix() : _range( 0 ) {};
  ~EvtMatrix();
  inline void setRange( int range );

  T& operator()( int row, int col ) { return _mat[ row ][ col ]; }
  T* operator[]( int row )          { return _mat[ row ];        }
  T det();
  EvtMatrix* min( int row, int col );
  EvtMatrix* inverse();
  std::string dump();

  template <class M> friend EvtMatrix< M >* operator*( const EvtMatrix< M >& left, const EvtMatrix< M >& right );
};




template <class T> inline void EvtMatrix< T >::setRange( int range )
{
  // If the range is changed, delete any previous matrix stored
  //    and allocate elements with the newly specified range.
  if ( _range != range )
    {
      if ( _range )
	{
	  for ( int row = 0; row < _range; row++ )
	    delete[] _mat[ row ];
	  delete[] _mat;
	}

      _mat = new T*[ range ];
      for ( int row = 0; row < range; row++ )
	_mat[ row ] = new T[ range ];

      // Set the new range.
      _range = range;
    }

  // Since user is willing to change the range, reset the matrix elements.
  for ( int row = 0; row < _range; row++ )
    for ( int col = 0; col < _range; col++ )
      _mat[ row ][ col ] = 0.;
}


template <class T> EvtMatrix< T >::~EvtMatrix()
{
  for( int row = 0; row < _range; row++ )
    delete[] _mat[ row ];
  delete[] _mat;
}


template <class T> std::string EvtMatrix< T >::dump()
{
  std::ostringstream str;

  for ( int row = 0; row < _range; row++ )
    {
      str << "|";
      for ( int col = 0; col < _range; col++ )
	str << "\t" << _mat[ row ][ col ];
      str << "\t|" << std::endl;
    }

  return str.str();
}


template <class T> T EvtMatrix< T >::det()
{
  if ( _range == 1 )
    return _mat[ 0 ][ 0 ];

  // There's no need to define the range 2 determinant manually, but it may
  //    speed up the calculation.
  if ( _range == 2 )
    return _mat[ 0 ][ 0 ] * _mat[ 1 ][ 1 ] - _mat[ 0 ][ 1 ] * _mat[ 1 ][ 0 ];

  T sum = 0.;

  for ( int col = 0; col < _range; col++ )
    {
      EvtMatrix< T >* minor = min( 0, col );
      sum += std::pow( -1., col ) * _mat[ 0 ][ col ] * minor->det();
      delete minor;
    }

  return sum;
}


// Returns the minor at (i, j).
template <class T> EvtMatrix< T >* EvtMatrix< T >::min( int row, int col )
{
  EvtMatrix< T >* minor = new EvtMatrix< T >();
  minor->setRange( _range - 1 );

  int minIndex = 0;

  for ( int r = 0; r < _range; r++ )
    for ( int c = 0; c < _range; c++ )
      if ( ( r != row ) && ( c != col ) )
	{
	  (*minor)( minIndex / ( _range - 1 ), minIndex % ( _range - 1 ) ) = _mat[ r ][ c ];
	  minIndex++;
	}

  return minor;
}


template <class T> EvtMatrix< T >* EvtMatrix< T >::inverse()
{
  EvtMatrix< T >* inv = new EvtMatrix< T >();
  inv->setRange( _range );

  if ( det() == 0 )
    {
      std::cerr << "This matrix has a null determinant and cannot be inverted. Returning zero matrix." << std::endl;
      for ( int row = 0; row < _range; row++ )
	for ( int col = 0; col < _range; col++ )
	  (*inv)( row, col ) = 0.;
      return inv;
    }

  T determinant = det();

  for ( int row = 0; row < _range; row++ )
    for ( int col = 0; col < _range; col++ )
      {
	EvtMatrix< T >* minor = min( row, col );
	inv->_mat[col][row] = std::pow( -1., row + col ) * minor->det() / determinant;
	delete minor;
      }

  return inv;
}


template <class T>
EvtMatrix< T >* operator*( const EvtMatrix< T >& left, const EvtMatrix< T >& right )
{
  // Chech that the matrices have the correct range.
  if ( left._range != right._range )
    {
      std::cerr << "These matrices cannot be multiplied." << std::endl;
      return new EvtMatrix< T >();
    }

  EvtMatrix< T >* mat = new EvtMatrix< T >();
  mat->setRange( left._range );

  // Initialize the elements of the matrix.
  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      (*mat)[ row ][ col ] = 0;

  for ( int row = 0; row < left._range; row++ )
    for ( int col = 0; col < right._range; col++ )
      for ( int line = 0; line < right._range; line++ )
	(*mat)[ row ][ col ] += left._mat[ row ][ line ] * right._mat[ line ][ col ];

  return mat;
}


#endif
