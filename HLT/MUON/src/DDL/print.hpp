////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DDL_PRINT_HPP
#define dHLT_DDL_PRINT_HPP

#include <ostream>

namespace dHLT
{
namespace DDL
{


/* Prints user readable bit pattern to the stream.
   bitcount specified the number bits to print of x, starting from the least
   significant bit to the most significant bit.
 */
void PrintBits(std::ostream& os, UInt x, int bitcount);

/* Prints the value x as a hexadecimal string.
 */
void PrintHex(std::ostream& os, UInt x, UChar width);


/* Implementing << operators to be able to print data structures to streams.
   eg.
       DDLHeader h;
       cout << h << endl;
 */
std::ostream& operator << (std::ostream& os, const DDLHeader& h);
std::ostream& operator << (std::ostream& os, const BlockHeader& h);
std::ostream& operator << (std::ostream& os, const DSPHeader& h);
std::ostream& operator << (std::ostream& os, const PatchBusHeader& h);
std::ostream& operator << (std::ostream& os, const ADCData& d);
std::ostream& operator << (std::ostream& os, const LocalData& d);
std::ostream& operator << (std::ostream& os, const RegionalData& d);
std::ostream& operator << (std::ostream& os, const EnhancedHeader& h);
std::ostream& operator << (std::ostream& os, const TriggerData& d);


}; // DDL
}; // dHLT

#endif // dHLT_DDL_PRINT_HPP
