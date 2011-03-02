#ifndef ALIZDCRECO_H
#define ALIZDCRECO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Classe for ZDC RecPoints                  //
////////////////////////////////////////////////

#include "TObject.h"

class AliZDCReco : public TObject {

public:
  AliZDCReco();
  AliZDCReco(Float_t* ezn1, Float_t* ezp1, Float_t* ezn2, Float_t* ezp2,  
	     Float_t* ezn1tow, Float_t* ezp1tow, Float_t* ezn2tow, Float_t* ezp2tow, 
	     Float_t* ezem1, Float_t* ezem2, Float_t* ref1, Float_t* ref2, 
	     //	   
	     Int_t detspnSideA,  Int_t detsppSideA, 
	     Int_t detspnSideC, Int_t detsppSideC,  
	     Int_t trsp, Int_t trspSideA, Int_t trspSideC,
	     Int_t npart, Int_t npartSideA, Int_t npartSideC, 
	     Float_t b, Float_t bSideA, Float_t bSideC,
	     UInt_t recoFlag, Bool_t energyFlag, Bool_t scalerOn, 
	     UInt_t* scaler, Int_t tdcData[32][4]);

  AliZDCReco(const AliZDCReco &oldreco);
  virtual ~AliZDCReco() {}

  // Getters 
  virtual Float_t GetZN1HREnergy()   const  {return fZN1Energy[0];}
  virtual Float_t GetZP1HREnergy()   const  {return fZP1Energy[0];}
  virtual Float_t GetZN2HREnergy()   const  {return fZN2Energy[0];}
  virtual Float_t GetZP2HREnergy()   const  {return fZP2Energy[0];}
  //
  virtual Float_t GetZN1LREnergy()   const  {return fZN1Energy[1];}
  virtual Float_t GetZP1LREnergy()   const  {return fZP1Energy[1];}
  virtual Float_t GetZN2LREnergy()   const  {return fZN2Energy[1];}
  virtual Float_t GetZP2LREnergy()   const  {return fZP2Energy[1];}
  //
  virtual Float_t GetZN1HREnTow(Int_t tow)  const {return fZN1EnTow[tow];}
  virtual Float_t GetZP1HREnTow(Int_t tow)  const {return fZP1EnTow[tow];}
  virtual Float_t GetZN2HREnTow(Int_t tow)  const {return fZN2EnTow[tow];}
  virtual Float_t GetZP2HREnTow(Int_t tow)  const {return fZP2EnTow[tow];}
  //
  virtual Float_t GetZN1LREnTow(Int_t tow)  const {return fZN1EnTow[tow+5];}
  virtual Float_t GetZP1LREnTow(Int_t tow)  const {return fZP1EnTow[tow+5];}
  virtual Float_t GetZN2LREnTow(Int_t tow)  const {return fZN2EnTow[tow+5];}
  virtual Float_t GetZP2LREnTow(Int_t tow)  const {return fZP2EnTow[tow+5];}
  //
  virtual Float_t GetZEM1HRsignal()   const  {return fZEM1signal[0];}
  virtual Float_t GetZEM1LRsignal()   const  {return fZEM1signal[1];}
  virtual Float_t GetZEM2HRsignal()   const  {return fZEM2signal[0];}
  virtual Float_t GetZEM2LRsignal()   const  {return fZEM2signal[1];}
  //
  virtual Float_t GetPMRef1HRsignal()   const  {return fZEM1signal[0];}
  virtual Float_t GetPMRef1LRsignal()   const  {return fZEM1signal[1];}
  virtual Float_t GetPMRef2HRsignal()   const  {return fZEM2signal[0];}
  virtual Float_t GetPMRef2LRsignal()   const  {return fZEM2signal[1];}
  //
  virtual Int_t   GetNDetSpecNSideA()  const {return fNDetSpecNSideA;}
  virtual Int_t   GetNDetSpecPSideA()  const {return fNDetSpecPSideA;}
  virtual Int_t   GetNDetSpecNSideC()  const {return fNDetSpecNSideC;}
  virtual Int_t   GetNDetSpecPSideC()  const {return fNDetSpecPSideC;}
  virtual Int_t   GetNTrueSpectators() const {return fNTrueSpectators;}
  virtual Int_t   GetNTrueSpecSideA()  const {return fNTrueSpecSideA;}
  virtual Int_t   GetNTrueSpecSideC()  const {return fNTrueSpecSideC;}
  virtual Int_t   GetNParticipants()   const {return fNParticipants;}
  virtual Int_t   GetNPartSideA()      const {return fNPartSideA;}
  virtual Int_t   GetNPartSideC()      const {return fNPartSideC;}
  virtual Float_t GetImpParameter()    const {return fImpParameter;}
  virtual Float_t GetImpParSideA()     const {return fImpParSideA;}
  virtual Float_t GetImpParSideC()     const {return fImpParSideC;}
  //
  virtual UInt_t  GetRecoFlag()      const {return fRecoFlag;}
  virtual UInt_t  GetZDCPattern()    const {return (fRecoFlag & 0x0000003f);}
  virtual UInt_t  GetChOnFlag()      const {return (fRecoFlag & 0x00000100);}
  virtual UInt_t  GetChOvflwFlag()   const {return (fRecoFlag & 0x00000200);}
  virtual UInt_t  GetChUndflwFlag()  const {return (fRecoFlag & 0x00000400);}
  //
  virtual Bool_t  GetEnergyFlag()    const {return fEnergyFlag;}
  virtual Bool_t  IsScalerOn()          const {return fIsScalerOn;}
  virtual UInt_t  GetZDCScaler(Int_t k) const {return fZDCScaler[k];}
  //
  virtual Int_t   GetZDCTDCData(Int_t j, Int_t k) const {return fZDCTDCData[j][k];}

  // Print method
  virtual void Print(Option_t *) const;

private:
  // Data members
  Float_t fZN1Energy[2]; // Energy detected in ZN1 (sum of 5 tower signals)
  Float_t fZP1Energy[2]; // Energy detected in ZP1 (sum of 5 tower signals)
  Float_t fZN2Energy[2]; // Energy detected in ZN2 (sum of 5 tower signals)
  Float_t fZP2Energy[2]; // Energy detected in ZP2 (sum of 5 tower signals)
  //
  Float_t fZN1EnTow[10]; // Energy in ZN1 towers
  Float_t fZP1EnTow[10]; // Energy in ZP1 towers
  Float_t fZN2EnTow[10]; // Energy in ZN2 towers
  Float_t fZP2EnTow[10]; // Energy in ZP2 towers
  //
  Float_t fZEM1signal[2];// Signal in EM1 ZDC
  Float_t fZEM2signal[2];// Signal in EM2 ZDC
  //
  Float_t fPMRef1[2];	 // Reference PM side C
  Float_t fPMRef2[2];	 // Reference PM side A
  //
  Int_t	  fNDetSpecNSideA; // Number of spectator neutrons detected
  Int_t	  fNDetSpecPSideA; // Number of spectator protons detected
  Int_t	  fNDetSpecNSideC; // Number of spectator neutrons detected
  Int_t	  fNDetSpecPSideC; // Number of spectator protons detected
  Int_t	  fNTrueSpectators;// Estimate of the total number of spectators
  Int_t	  fNTrueSpecSideA; // Estimate of the number of spectators side A
  Int_t	  fNTrueSpecSideC; // Estimate of the number of spectators side C
  Int_t	  fNParticipants;  // Estimate of the total number of participants
  Int_t	  fNPartSideA;	   // Estimate of the number of participants side A
  Int_t	  fNPartSideC;	   // Estimate of the number of participants side C
  Float_t fImpParameter;   // Estimate of the impact parameter
  Float_t fImpParSideA;	   // Estimate of the impact parameter side A
  Float_t fImpParSideC;	   // Estimate of the impact parameter side B
  //
  UInt_t  fRecoFlag;       // Reconstruction flag
  Bool_t  fEnergyFlag;     // Is the reco value in energy?
  Bool_t  fIsScalerOn;     // True if scaler has been read in the event
  UInt_t  fZDCScaler[32];  // Counts from ZDC VME scaler
  //
  Int_t fZDCTDCData[32][4];      // TDC data raw

  ClassDef(AliZDCReco,14)  // RecPoints for the Zero Degree Calorimeters
};
 
#endif
