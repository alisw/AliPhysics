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

namespace dHLT
{

std::ostream& operator << (std::ostream& os, const EventID& id);
std::ostream& operator << (std::ostream& os, const Point& p);
std::ostream& operator << (std::ostream& os, const ParticleSign s);
std::ostream& operator << (std::ostream& os, const TriggerRecord& rec);
std::ostream& operator << (std::ostream& os, const ChamberID chamber);

}; // dHLT

#endif // dHLT_DEBUG_PRINT_ROUTINES_HPP
