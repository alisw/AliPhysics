#ifndef ALIAODVZERO_H
#define ALIAODVZERO_H

//-------------------------------------------------------------------------
//     Container class for AOD VZERO data
//     Author: Cvetan Cheshkov
//     cvetan.cheshkov@cern.ch 2/02/2011
//-------------------------------------------------------------------------

#include "AliVVZERO.h"

class AliAODVZERO : public AliVVZERO 
{
public:
  AliAODVZERO();
  AliAODVZERO(const AliAODVZERO& source);
  AliAODVZERO(const AliVVZERO& source);
  AliAODVZERO &operator=(const AliAODVZERO& source);
  AliAODVZERO &operator=(const AliVVZERO& source);

  virtual ~AliAODVZERO() {};

  void SetMultiplicity(Float_t Multiplicity[64])
    {for(Int_t i=0;i<64;i++) fMultiplicity[i]=Multiplicity[i];}

  // Getters  
  virtual Short_t  GetNbPMV0A() const;
  virtual Short_t  GetNbPMV0C() const;
  virtual Float_t  GetMTotV0A() const;
  virtual Float_t  GetMTotV0C() const; 
  virtual Float_t  GetMRingV0A(Int_t ring) const;
  virtual Float_t  GetMRingV0C(Int_t ring) const;

  virtual Float_t  GetMultiplicity(Int_t i) const;
  virtual Float_t  GetMultiplicityV0A(Int_t i) const;
  virtual Float_t  GetMultiplicityV0C(Int_t i) const;    

  virtual Bool_t   BBTriggerV0A(Int_t i) const;
  virtual Bool_t   BGTriggerV0A(Int_t i) const;
  virtual Bool_t   BBTriggerV0C(Int_t i) const;
  virtual Bool_t   BGTriggerV0C(Int_t i) const;  
  virtual Bool_t   GetBBFlag(Int_t i) const;
  virtual Bool_t   GetBGFlag(Int_t i) const;

  virtual Float_t  GetV0ATime() const { return fV0ATime; }
  virtual Float_t  GetV0CTime() const { return fV0CTime; }

  virtual Decision GetV0ADecision() const { return fV0ADecision; }
  virtual Decision GetV0CDecision() const { return fV0CDecision; }

  virtual UShort_t GetTriggerChargeA() const { return fTriggerChargeA; }
  virtual UShort_t GetTriggerChargeC() const { return fTriggerChargeC; }
  virtual UShort_t GetTriggerBits() const { return fTriggerBits; }
  
protected:

  UInt_t  fBBtriggerV0A;     // bit mask for Beam-Beam trigger in V0A
  UInt_t  fBGtriggerV0A;     // bit mask for Beam-Gas trigger in V0A
  UInt_t  fBBtriggerV0C;     // bit mask for Beam-Beam trigger in V0C
  UInt_t  fBGtriggerV0C;     // bit mask for Beam-Gas trigger in V0C

  Float_t fMultiplicity[64]; //  multiplicity for each channel

  Bool_t  fBBFlag[64];       //  BB Flags from Online V0 Electronics
  Bool_t  fBGFlag[64];       //  BG Flags from Online V0 Electronics

  Float_t fV0ATime;          // Average time in V0A
  Float_t fV0CTime;          // Average time in V0C

  Decision fV0ADecision;     // V0A final decision based on average time of channels
  Decision fV0CDecision;     // V0C final decision based on average time of channels

  UShort_t fTriggerChargeA;  // Sum of the trigger (clock=10) charge on A side
  UShort_t fTriggerChargeC;  // Sum of the trigger (clock=10) charge on C side
  UShort_t fTriggerBits;     // V0 trigger bits as defined in the firmware

  ClassDef(AliAODVZERO,2)
};

#endif
