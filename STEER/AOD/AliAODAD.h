#ifndef ALIAODAD_H
#define ALIAODAD_H

//-------------------------------------------------------------------------
//     Container class for AOD AD data
//     Author: Michal Broz
//     michal.broz@cern.ch 
//-------------------------------------------------------------------------

#include <AliVAD.h>

class AliAODAD : public AliVAD 
{
public:
  AliAODAD();
  AliAODAD(const AliAODAD& source);
  AliAODAD(const AliVAD& source);
  AliAODAD &operator=(const AliAODAD& source);
  AliAODAD &operator=(const AliVAD& source);

  virtual ~AliAODAD() {};
   
  // Getters  
  virtual Short_t  GetNbPMADA() const;
  virtual Short_t  GetNbPMADC() const;
  virtual Float_t  GetMTotADA() const;
  virtual Float_t  GetMTotADC() const; 

  virtual Float_t  GetMultiplicity(Int_t i) const;
  virtual Float_t  GetMultiplicityADA(Int_t i) const;
  virtual Float_t  GetMultiplicityADC(Int_t i) const;    

  virtual Bool_t   BBTriggerADA(Int_t i) const;
  virtual Bool_t   BGTriggerADA(Int_t i) const;
  virtual Bool_t   BBTriggerADC(Int_t i) const;
  virtual Bool_t   BGTriggerADC(Int_t i) const;  
  virtual Bool_t   GetBBFlag(Int_t i) const;
  virtual Bool_t   GetBGFlag(Int_t i) const;

  virtual Float_t  GetADATime() const { return fADATime; }
  virtual Float_t  GetADCTime() const { return fADCTime; }

  virtual Decision GetADADecision() const { return fADADecision; }
  virtual Decision GetADCDecision() const { return fADCDecision; }

  virtual UShort_t GetTriggerChargeA() const { return fTriggerChargeA; }
  virtual UShort_t GetTriggerChargeC() const { return fTriggerChargeC; }
  virtual UShort_t GetTriggerBits() const { return fTriggerBits; }
  
  virtual Bool_t   GetPFBBFlag(Int_t channel, Int_t clock) const { return fIsBB[channel][clock]; } 
  virtual Bool_t   GetPFBGFlag(Int_t channel, Int_t clock) const { return fIsBG[channel][clock]; }  
  
protected:

  UInt_t  fBBtriggerADA;     // bit mask for Beam-Beam trigger in ADA
  UInt_t  fBGtriggerADA;     // bit mask for Beam-Gas trigger in ADA
  UInt_t  fBBtriggerADC;     // bit mask for Beam-Beam trigger in ADC
  UInt_t  fBGtriggerADC;     // bit mask for Beam-Gas trigger in ADC

  Float_t fMultiplicity[16]; //  multiplicity for each channel

  Bool_t  fBBFlag[16];       //  BB Flags from Online AD Electronics
  Bool_t  fBGFlag[16];       //  BG Flags from Online AD Electronics

  Float_t fADATime;          // Average time in ADA
  Float_t fADCTime;          // Average time in ADC

  Decision fADADecision;     // ADA final decision based on average time of channels
  Decision fADCDecision;     // ADC final decision based on average time of channels

  UShort_t fTriggerChargeA;  // Sum of the trigger (clock=10) charge on A side
  UShort_t fTriggerChargeC;  // Sum of the trigger (clock=10) charge on C side
  UShort_t fTriggerBits;     // AD trigger bits as defined in the firmware
  
  Bool_t   fIsBB[16][21];  // BB flag for all channels and 21 clocks
  Bool_t   fIsBG[16][21];  // BG flag for all channels and 21 clocks

  ClassDef(AliAODAD,2)
};

#endif
