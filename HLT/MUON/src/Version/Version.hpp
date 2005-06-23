////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* This source file defines the routines used to extract the version 
   information for the muon HLT module. 
 */
  
#ifndef dHLT_VERSION_FUNCTIONS_HPP
#define dHLT_VERSION_FUNCTIONS_HPP

#include "BasicTypes.hpp"

namespace dHLT
{


/* Returns a NULL terminated version string of the form: XX.YY.ZZ
   where XX will be the MajorVersion() number, YY the MinorVersion() number
   and ZZ the BuildVersion() build number.
  */
const Char* VersionString();

/* Returns the major version number of the muon HLT module.
   Major version numbers are incremented after major and fundemental changes 
   to the muon HLT module. Fundemental changes include: changing large
   parts of the code, e.g. changing existing class behaviour or completely
   rewriting classes and methods.
 */
UInt MajorVersion();

/* Returns the minor version number of the muon HLT module.
   Minor version numbers are incremented for any small changes made to the 
   module such as adding new classes, fixing methods etc...
 */
UInt MinorVersion();

/* Returns the build number of the muon HLT module.
   The build number is incremented every time the source is built for distribution.
 */
UInt BuildNumber();


} // dHLT

#endif // dHLT_VERSION_FUNCTIONS_HPP
