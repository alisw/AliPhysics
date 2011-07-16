#ifndef ALILHCCLOCKPHASE_H
#define ALILHCCLOCKPHASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////////
//                          Class AliLHCClockPhase                          //
//   Container class for storing of the LHC clock phase.                    //
//   The source of the the data are BPTXs - they measure                    //
//   The beam pick-up time w.r.t to the LHC clock distributed by CTP.       //
//   The values stored by DCS are always relative to some fixed reference   //
//   moment.                                                                //
//                                                                          //
//   cvetan.cheshkov@cern.ch 21/07/2010                                     //
//////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TObjArray.h>

class AliDCSValue;

class AliLHCClockPhase : public TObject {

 public:
  AliLHCClockPhase();
  virtual ~AliLHCClockPhase() {}

  Int_t    GetNumberOfPhaseB1DPs()   const {return fPhaseB1.GetEntries();}
  Int_t    GetNumberOfPhaseB2DPs()   const {return fPhaseB2.GetEntries();}
  const AliDCSValue* GetPhaseB1DP(Int_t index) const;
  const AliDCSValue* GetPhaseB2DP(Int_t index) const;

  Float_t  GetMeanPhaseB1() const;
  Float_t  GetMeanPhaseB2() const;
  Float_t  GetMeanPhase()   const;
  Float_t  GetPhaseB1(UInt_t timestamp) const;
  Float_t  GetPhaseB2(UInt_t timestamp) const;
  Float_t  GetPhase(UInt_t timestamp)   const;

  void     AddPhaseB1DP(UInt_t timestamp, Float_t phase);
  void     AddPhaseB2DP(UInt_t timestamp, Float_t phase);

  virtual void Print( const Option_t* opt ="" ) const;

 private:
  AliLHCClockPhase(const AliLHCClockPhase &phase);
  AliLHCClockPhase& operator= (const AliLHCClockPhase& phase);

  TObjArray fPhaseB1;              // Array of AliDCSValue's containing the phase as measure by BPTX on beam1
  TObjArray fPhaseB2;              // Array of AliDCSValue's containing the phase as measure by BPTX on beam2

  ClassDef(AliLHCClockPhase,1)     // LHC-clock phase container class
};

#endif
