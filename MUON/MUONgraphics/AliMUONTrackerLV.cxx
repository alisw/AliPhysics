#include "AliMUONTrackerLV.h"
#include "AliLog.h"

///\cond CLASSIMP
ClassImp(AliMUONTrackerLV)
///\endcond

AliMUONTrackerLV::AliMUONTrackerLV(const char* runlist, const char* ocdbPath)
 : AliMUONTrackerVoltages(runlist,ocdbPath)
{
  fOCDBObjectPath = "MUON/Calib/LV";
}

AliMUONTrackerLV::AliMUONTrackerLV(Int_t runNumber, const char* ocdbPath)
: AliMUONTrackerVoltages(runNumber,ocdbPath)
{
  fOCDBObjectPath = "MUON/Calib/LV";
}

void AliMUONTrackerLV::Scan(Int_t verbose)
{
  AliError("Not implemented (yet ?)");
}
