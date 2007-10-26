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
  AliZDCReco(Float_t ezn1, Float_t ezp1, Float_t ezn2, Float_t ezp2,  
	     Float_t* ezn1tow, Float_t* ezp1tow,
	     Float_t* ezn2tow, Float_t* ezp2tow, 
	     Float_t ezem1, Float_t ezem2, 
	     //	   
	     Int_t detspnLeft,  Int_t detsppLeft, Int_t detspnRight, Int_t detsppRight,  
	     Int_t trspnLeft, Int_t trsppLeft, Int_t trspLeft, 
	     Int_t trspnRight, Int_t trsppRight, Int_t trspRight,
	     Int_t partLeft, Int_t partRight,  
	     Float_t b);

  AliZDCReco(const AliZDCReco &oldreco);
  virtual ~AliZDCReco() {}

  // Getters 
  virtual Float_t GetZN1Energy()   const  {return fZN1Energy;}
  virtual Float_t GetZP1Energy()   const  {return fZP1Energy;}
  virtual Float_t GetZN2Energy()   const  {return fZN2Energy;}
  virtual Float_t GetZP2Energy()   const  {return fZP2Energy;}
  //
  virtual Float_t GetZN1EnTow(Int_t tow)  const {return fZN1EnTow[tow];}
  virtual Float_t GetZP1EnTow(Int_t tow)  const {return fZP1EnTow[tow];}
  virtual Float_t GetZN2EnTow(Int_t tow)  const {return fZN2EnTow[tow];}
  virtual Float_t GetZP2EnTow(Int_t tow)  const {return fZP2EnTow[tow];}
  //
  virtual Float_t GetZEM1signal()   const  {return fZEM1signal;}
  virtual Float_t GetZEM2signal()   const  {return fZEM2signal;}
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
  Float_t fZN1Energy;	// Energy detected in ZN1 (sum of 5 tower signals)
  Float_t fZP1Energy;	// Energy detected in ZP1 (sum of 5 tower signals)
  Float_t fZN2Energy;	// Energy detected in ZN2 (sum of 5 tower signals)
  Float_t fZP2Energy;	// Energy detected in ZP2 (sum of 5 tower signals)
  //
  Float_t fZN1EnTow[5];	// Energy in ZN1 towers
  Float_t fZP1EnTow[5]; // Energy in ZP1 towers
  Float_t fZN2EnTow[5];	// Energy in ZN2 towers
  Float_t fZP2EnTow[5]; // Energy in ZP2 towers
  //
  Float_t fZEM1signal;	// Signal in EM1 ZDC
  Float_t fZEM2signal;	// Signal in EM2 ZDC
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


  ClassDef(AliZDCReco,3)  // RecPoints for the Zero Degree Calorimeters
};
 
#endif
