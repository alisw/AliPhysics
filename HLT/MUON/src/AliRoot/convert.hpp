////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCONVERT_H
#define ALIHLTMUONCONVERT_H

#include "../AliHLTMUONCorePoint.h"
#include "AliRoot/Point.hpp"
#include "../AliHLTMUONCoreTriggerRecord.h"
#include "AliRoot/TriggerRecord.hpp"
#include "../AliHLTMUONCoreTrack.h"
#include "AliRoot/Track.hpp"
#include "../AliHLTMUONCoreRegionOfInterest.h"
#include "AliRoot/Region.hpp"

/* Conversion routines to convert from dHLT structures to AliRoot structures
   and back again.
 */
extern AliHLTMUONPoint             AliHLTMUONConvert(const AliHLTMUONCorePoint&         point);
extern AliHLTMUONCorePoint         AliHLTMUONConvert(const AliHLTMUONPoint&             point);
extern AliHLTMUONTriggerRecord     AliHLTMUONConvert(const AliHLTMUONCoreTriggerRecord& record, Int_t triggernumber);
extern AliHLTMUONCoreTriggerRecord AliHLTMUONConvert(const AliHLTMUONTriggerRecord&     record);
extern AliHLTMUONTrack             AliHLTMUONConvert(const AliHLTMUONCoreTrack&         track);
extern AliHLTMUONCoreTrack         AliHLTMUONConvert(const AliHLTMUONTrack&             track);
extern AliHLTMUONRegion            AliHLTMUONConvert(const AliHLTMUONCoreROI            region);
extern AliHLTMUONCoreROI           AliHLTMUONConvert(const AliHLTMUONRegion&            region, UInt_t chamber);


#endif // ALIHLTMUONCONVERT_H
