#ifndef ALIMUONTRACKERLV_H
#define ALIMUONTRACKERLV_H

#ifndef ALIMUONTRACKERVOLTAGES_H
#  include "AliMUONTrackerVoltages.h"
#endif

class AliMUONTrackerLV : public AliMUONTrackerVoltages
{
public:
  AliMUONTrackerLV(const char* runlist, const char* ocdbPath);
  AliMUONTrackerLV(Int_t runNumber, const char* ocdbPath);

  virtual ~AliMUONTrackerLV() {};

  virtual void Scan(Int_t verbose=0);

  void CheckLV(Int_t runNumber, Int_t verbose=0);

private:

  AliMUONTrackerLV(const AliMUONTrackerLV& rhs); // not implemented on purpose
  AliMUONTrackerLV& operator=(const AliMUONTrackerLV& rhs); // not implemented on purpose

  ClassDef(AliMUONTrackerLV,1) // Utility class to inspect MUON Tracker LV values
};

#endif
