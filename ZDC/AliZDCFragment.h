#ifndef ALIZDCFRAGMENT_H
#define ALIZDCFRAGMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


////////////////////////////////////////////////////
//                                                //  
//     Generate nuclear fragments parametrizing   //
//       resuslts of SIS and SPS energies	  //
//                                                //
////////////////////////////////////////////////////


#include <TMath.h>

extern int comp(const void *, const void *);
 
class AliZDCFragment : public TNamed {

public:
  AliZDCFragment();
  AliZDCFragment(Float_t b);
  virtual      ~AliZDCFragment() {}
  void GenerateIMF(Int_t* fZZ, Int_t &fNalpha);
  void AttachNeutrons(Int_t* fZZ, Int_t* fNN, Int_t &Ztot, Int_t &Ntot);
  
  // Setting parameters
  virtual void SetImpactParameter(Float_t b) {fB=b;};
  
  // Getting parameters
  Int_t GetFragmentNum() {return fNimf;};
  
 
protected:
  
   Float_t  fB; 	 // Impact parameter
   Float_t  fZbAverage ; // Mean value of Z bound 
   Int_t    fNimf;	 // Number of IMF
   Float_t  fZmax;	 // Mean value of maximum Z of fragment
   Float_t  fTau;	 // Exponent of charge distribution: dN/dZ = Z*exp(-fTau)
   Int_t    fZZ[100];	 // Array of atomic numbers of fragments
   Int_t    fNN[100];	 // Array of number of neutrons of fragments
   Int_t    fNalpha;	 // Number of alpha particles
   Int_t    fZtot;	 // Total number of bound protons
   Int_t    fNtot;	 // Total number of bound neutrons

  
   ClassDef(AliZDCFragment,1)  // Generator for AliZDC fragment class
};

#endif
