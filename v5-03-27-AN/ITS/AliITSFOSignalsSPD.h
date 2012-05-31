#ifndef ALIITS_FOSIGNALSSPD_H
#define ALIITS_FOSIGNALSSPD_H

/////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                          //
//                                                                 //
// This class is used to store information on generated Fast-OR    //
// signals. 1200 bits, one per pixel chip.                         //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TBits.h>

class AliITSFOSignalsSPD : public TObject {

 public:
  AliITSFOSignalsSPD();
  virtual ~AliITSFOSignalsSPD();
  AliITSFOSignalsSPD(const AliITSFOSignalsSPD& fo);
  AliITSFOSignalsSPD& operator=(const AliITSFOSignalsSPD& fo);

  virtual void    ResetSignals() {fSignals.ResetAllBits();}
  virtual void    SetSignal(UInt_t eq, UInt_t hs, UInt_t chip, Bool_t setVal=kTRUE);
  virtual Bool_t  GetSignal(UInt_t eq, UInt_t hs, UInt_t chip) const;

  virtual Bool_t  GetNextSignal(Int_t& eq, Int_t& hs, Int_t& chip) const;
  virtual void    DumpSignals();

 protected:
  TBits fSignals; // FO signals, one bit per chip

  UInt_t  GetChipKey(Int_t eq, Int_t hs, Int_t chip) const;
  void    GetChipFromKey(UInt_t key, Int_t& eq, Int_t& hs, Int_t& chip) const;

  ClassDef(AliITSFOSignalsSPD,1)
};

#endif
