////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_NEW_HPP
#define dHLT_NEW_HPP

#include "Error.hpp"

/* Overload the global new and delete operators to perform better error checking.
   These new operators throw dHLT::OutOfMemory exceptions if the system is out of
   memory.
 */
void* operator new (size_t size) throw (std::bad_alloc);
void* operator new [] (size_t size) throw (std::bad_alloc);
void operator delete (void* memory) throw ();
void operator delete [] (void* memory) throw ();

#endif // dHLT_NEW_HPP
