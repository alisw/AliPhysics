#ifndef ALI_ITS_ONLINESPDPHYSINFO_H
#define ALI_ITS_ONLINESPDPHYSINFO_H

/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// It holds information needed for a physics run.              //
// This class should only be used through the interface of the //
// AliITSOnlineSPDphys class.                                  //
/////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TArrayI.h>

class AliITSOnlineSPDphysInfo :  public TObject {

 public:
  AliITSOnlineSPDphysInfo();
  virtual ~AliITSOnlineSPDphysInfo();

  virtual void   ClearThis();
  // SET METHODS ***********************************
  void     AddRunNr(UInt_t val);
  void     SetEqNr(UInt_t val) {fEqNr=val;}
  void     SetNrEvents(UInt_t val) {fNrEvents=val;}
  void     AddNrEvents(Int_t val);
  void     IncrementNrEvents() {fNrEvents++;}

  // GET METHODS ***********************************
  UInt_t   GetNrRuns() const {return fNrRuns;}
  UInt_t   GetRunNr(UInt_t posi) const ;
  UInt_t   GetEqNr() const {return fEqNr;}
  UInt_t   GetNrEvents() const {return fNrEvents;}

 protected:
  UInt_t   fNrRuns;               // nr of runs used to fill hitmap
  TArrayI  fRunNrs;               // list of run nrs for the hitmap
  UInt_t   fEqNr;                 // eq nr
  UInt_t   fNrEvents;             // number of events

  ClassDef(AliITSOnlineSPDphysInfo,1)
    };
    
#endif
