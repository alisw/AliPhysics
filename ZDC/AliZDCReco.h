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
  AliZDCReco(Float_t ezn1, Float_t ezp1, Float_t ezdc1, Float_t ezem, 
  	     Float_t ezn2, Float_t ezp2, Float_t ezdc2, Int_t detspnLeft, 
             Int_t detsppLeft, Int_t detspnRight, Int_t detsppRight, 
	     Int_t trspn, Int_t trspp, Int_t trsp, Int_t part, Float_t b);
  AliZDCReco(const AliZDCReco &oldreco);
  virtual ~AliZDCReco() {}

  // Getters 
  virtual Float_t GetZN1energy()   const  {return fZN1energy;}
  virtual Float_t GetZP1energy()   const  {return fZP1energy;}
  virtual Float_t GetZDC1energy()  const  {return fZDC1energy;}
  virtual Float_t GetZN2energy()   const  {return fZN2energy;}
  virtual Float_t GetZP2energy()   const  {return fZP2energy;}
  virtual Float_t GetZDC2energy()  const  {return fZDC2energy;}
  virtual Float_t GetZEMenergy()   const  {return fZEMenergy;}
  virtual Int_t   GetNDetSpecNLeft()  const  {return fNDetSpecNLeft;}
  virtual Int_t   GetNDetSpecPLeft()  const  {return fNDetSpecPLeft;}
  virtual Int_t   GetNDetSpecNRight() const  {return fNDetSpecNRight;}
  virtual Int_t   GetNDetSpecPRight() const  {return fNDetSpecPRight;}
  virtual Int_t   GetNTrueSpecN()  const  {return fNTrueSpecN;}
  virtual Int_t   GetNTrueSpecP()  const  {return fNTrueSpecP;}
  virtual Int_t   GetNTrueSpec()   const  {return fNTrueSpec;}
  virtual Int_t   GetNPart()       const  {return fNPart;}
  virtual Float_t GetImpPar()      const  {return fImpPar;}

  // Print method
  virtual void Print(Option_t *) const;

private:
  // Data members
  Float_t fZN1energy;	// Energy detected in neutron ZDC
  Float_t fZP1energy;	// Energy detected in proton ZDC
  Float_t fZDC1energy;	// Total hadronic energy detcted in ZDCs
  Float_t fZN2energy;	// Energy detected in neutron ZDC
  Float_t fZP2energy;	// Energy detected in proton ZDC
  Float_t fZDC2energy;	// Total hadronic energy detcted in ZDCs
  Float_t fZEMenergy;	// Energy detected in EM ZDC
  Int_t	  fNDetSpecNLeft;  // Number of spectator neutrons detected
  Int_t	  fNDetSpecPLeft;  // Number of spectator protons detected
  Int_t	  fNDetSpecNRight; // Number of spectator neutrons detected
  Int_t	  fNDetSpecPRight; // Number of spectator protons detected
  Int_t	  fNTrueSpecN;  // Estimate of the number of spectator neutrons generated
  Int_t	  fNTrueSpecP;  // Estimate of the number of spectator protons generated
  Int_t	  fNTrueSpec ;  // Estimate of the total number of spectators
  Int_t	  fNPart;	// Estimate of the number of participants for 1 nucleus
  Float_t fImpPar;	// Estimate of the impact parameter


  ClassDef(AliZDCReco,1)  // RecPoints for the Zero Degree Calorimeters
};
 
#endif
