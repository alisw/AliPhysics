////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_CONVERT_HPP
#define dHLT_ALIROOT_CONVERT_HPP

#include "../Point.hpp"
#include "AliRoot/Point.hpp"
#include "../TriggerRecord.hpp"
#include "AliRoot/TriggerRecord.hpp"
#include "../Track.hpp"
#include "AliRoot/Track.hpp"
#include "../RegionOfInterest.hpp"
#include "AliRoot/Region.hpp"

namespace dHLT
{
namespace AliRoot
{


/* Conversion routines to convert from dHLT structures to AliRoot structures
   and back again.
 */
extern AliMUONHLT::Point         Convert(const dHLT::Point&               point);
extern dHLT::Point               Convert(const AliMUONHLT::Point&         point);
extern AliMUONHLT::TriggerRecord Convert(const dHLT::TriggerRecord&       record, Int_t triggernumber);
extern dHLT::TriggerRecord       Convert(const AliMUONHLT::TriggerRecord& record);
extern AliMUONHLT::Track         Convert(const dHLT::Track&               track);
extern dHLT::Track               Convert(const AliMUONHLT::Track&         track);
extern AliMUONHLT::Region        Convert(const dHLT::ROI                  region);
extern dHLT::ROI                 Convert(const AliMUONHLT::Region&        region, UInt_t chamber);


} // AliRoot
} // dHLT

#endif // dHLT_ALIROOT_CONVERT_HPP
