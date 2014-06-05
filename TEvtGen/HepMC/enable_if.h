#ifndef ENABLE_IF_H
#define ENABLE_IF_H
// author: Walter Brown

// ----------------------------------------------------------------------
// prolog

namespace HepMC  {
namespace detail {


// ----------------------------------------------------------------------
// enable_if<>

/// internal - used to decide if a class is arithmetic
template< bool, class >
struct enable_if
{ };

/// internal - use if class T is arithmetic
template< class T >
struct enable_if<true, T>
{
  typedef  T  type;	//!< check type of class T
};


// ----------------------------------------------------------------------
// disable_if<>

/// internal - used by SimpleVector to decide if a class is arithmetic
template< bool, class >
struct disable_if
{ };

/// internal - used by SimpleVector to decide if a class is arithmetic
template< class T >
struct disable_if<false, T>
{
  typedef  T  type;	//!< check type of class T
};


// ----------------------------------------------------------------------
// epilog

}  // namespace detail
}  // namespace HepMC

#endif  // ENABLE_IF_H
