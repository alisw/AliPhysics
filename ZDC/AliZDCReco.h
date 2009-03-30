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
	     Int_t detspnLeft,  Int_t detsppLeft, Int_t detspnRight, Int_t detsppRight,  
	     Int_t trspnLeft, Int_t trsppLeft, Int_t trspLeft, 
	     Int_t trspnRight, Int_t trsppRight, Int_t trspRight,
	     Int_t partLeft, Int_t partRight, Float_t b);

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
  virtual Int_t   GetNDetSpecNLeft()   const {return fNDetSpecNLeft;}
  virtual Int_t   GetNDetSpecPLeft()   const {return fNDetSpecPLeft;}
  virtual Int_t   GetNDetSpecNRight()  const {return fNDetSpecNRight;}
  virtual Int_t   GetNDetSpecPRight()  const {return fNDetSpecPRight;}
  virtual Int_t   GetNTrueSpecNLeft()  const {return fNTrueSpecNLeft;}
  virtual Int_t   GetNTrueSpecPLeft()  const {return fNTrueSpecPLeft;}
  virtual Int_t   GetNTrueSpecLeft()   const {return fNTrueSpecLeft;}
  virtual Int_t   GetNTrueSpecNRight() const {return fNTrueSpecNRight;}
  virtual Int_t   GetNTrueSpecPRight() const {return fNTrueSpecPRight;}
  virtual Int_t   GetNTrueSpecRight()  const {return fNTrueSpecRight;}
  virtual Int_t   GetNPartLeft()       const {return fNPartLeft;}
  virtual Int_t   GetNPartRight()      const {return fNPartRight;}
  virtual Float_t GetImpPar()          const {return fImpPar;}

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
  Int_t	  fNDetSpecNLeft;  // Number of spectator neutrons detected
  Int_t	  fNDetSpecPLeft;  // Number of spectator protons detected
  Int_t	  fNDetSpecNRight; // Number of spectator neutrons detected
  Int_t	  fNDetSpecPRight; // Number of spectator protons detected
  Int_t	  fNTrueSpecNLeft; // Estimate of the number of spectator neutrons generated
  Int_t	  fNTrueSpecPLeft; // Estimate of the number of spectator protons generated
  Int_t	  fNTrueSpecLeft;  // Estimate of the total number of spectators
  Int_t	  fNTrueSpecNRight;// Estimate of the number of spectator neutrons generated
  Int_t	  fNTrueSpecPRight;// Estimate of the number of spectator protons generated
  Int_t	  fNTrueSpecRight; // Estimate of the total number of spectators
  Int_t	  fNPartLeft;	// Estimate of the number of participants for 1 nucleus
  Int_t	  fNPartRight;	// Estimate of the number of participants for 1 nucleus
  Float_t fImpPar;	// Estimate of the impact parameter


  ClassDef(AliZDCReco,5)  // RecPoints for the Zero Degree Calorimeters
};
 
#endif
