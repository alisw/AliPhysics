#ifndef IS_ARITHMETIC
#define IS_ARITHMETIC
// author: Walter Brown

// ----------------------------------------------------------------------
// prolog

namespace HepMC  {

///
/// \namespace detail
/// internal namespace
///
namespace detail  {


// ----------------------------------------------------------------------
// is_arithmetic<>

/// undefined and therefore non-arithmetic
template< class T >
  struct is_arithmetic
{
  static  bool const  value = false;
};

/// character is arithmetic
template<>
  struct is_arithmetic<char>
{ static  bool const  value = true; };

/// unsigned character is arithmetic
template<>
  struct is_arithmetic<unsigned char>
{ static  bool const  value = true; };

/// signed character is arithmetic
template<>
  struct is_arithmetic<signed char>
{ static  bool const  value = true; };

/// short is arithmetic
template<>
  struct is_arithmetic<short>
{ static  bool const  value = true; };

/// unsigned short is arithmetic
template<>
  struct is_arithmetic<unsigned short>
{ static  bool const  value = true; };

/// int is arithmetic
template<>
  struct is_arithmetic<int>
{ static  bool const  value = true; };

/// unsigned int is arithmetic
template<>
  struct is_arithmetic<unsigned int>
{ static  bool const  value = true; };

/// long is arithmetic
template<>
  struct is_arithmetic<long>
{ static  bool const  value = true; };

/// unsigned long is arithmetic
template<>
  struct is_arithmetic<unsigned long>
{ static  bool const  value = true; };

/// float is arithmetic
template<>
  struct is_arithmetic<float>
{ static  bool const  value = true; };

/// double is arithmetic
template<>
  struct is_arithmetic<double>
{ static  bool const  value = true; };

/// long double is arithmetic
template<>
  struct is_arithmetic<long double>
{ static  bool const  value = true; };


// ----------------------------------------------------------------------
// epilog

}  // namespace detail
}  // namespace HepMC

#endif  // IS_ARITHMETIC
