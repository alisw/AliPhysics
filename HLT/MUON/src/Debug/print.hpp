////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DEBUG_PRINT_ROUTINES_HPP
#define dHLT_DEBUG_PRINT_ROUTINES_HPP

#include <ostream>

#include "EventID.hpp"
#include "Point.hpp"
#include "TriggerRecord.hpp"
#include "RegionOfInterest.hpp"

std::ostream& operator << (std::ostream& os, const AliHLTMUONCoreEventID& id);
std::ostream& operator << (std::ostream& os, const AliHLTMUONCorePoint& p);
std::ostream& operator << (std::ostream& os, const AliHLTMUONCoreParticleSign s);
std::ostream& operator << (std::ostream& os, const AliHLTMUONCoreTriggerRecord& rec);
std::ostream& operator << (std::ostream& os, const AliHLTMUONCoreChamberID chamber);

#endif // dHLT_DEBUG_PRINT_ROUTINES_HPP
