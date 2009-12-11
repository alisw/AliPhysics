#ifndef ALITRDCALPIDLQ_H
#define ALITRDCALPIDLQ_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// PID distributions for the LQ method                                    //
//                                                                        //
// Author:                                                                //
// Alex Bercuci <A.Bercuci@gsi.de>                                        //
//                                                                        //          
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDCALPID_H
#include "AliTRDCalPID.h"
#endif
#ifndef ROOT_TMath
#include "TMath.h"
#endif

class AliTRDCalPIDLQ : public AliTRDCalPID
{
public:
  enum ETRDCalPIDLQ {
    kNLength = 4 // No of bins for tracklet length discretization 
   ,kNN2LQtransition = 4 // index of NN slices where first LQ slice ends 
  };
	
  AliTRDCalPIDLQ();
  AliTRDCalPIDLQ(const Text_t *name, const Text_t *title);
  virtual        ~AliTRDCalPIDLQ(){};

  inline static Bool_t  CookdEdx(const Float_t * const dedx, Double_t *x);
  TObject*              GetModel(Int_t ip, Int_t iType, Int_t iPlane) const;
  static Double_t       GetLength(Int_t il) { return (il<0 || il>=kNLength) ? -1. : fgTrackSegLength[il]; }
  inline static Int_t   GetModelID(Int_t mom, Int_t spec, Int_t ii=-1);
  Double_t              GetProbability(Int_t spec, Float_t mom
                               , const Float_t * const dedx
                               , Float_t length, Int_t plane) const;
  Bool_t                LoadReferences(Char_t* refFile);

protected:
  static Float_t  fgTrackSegLength[kNLength]; // Track segment lengths

private:
  AliTRDCalPIDLQ(const AliTRDCalPIDLQ& pd);
  AliTRDCalPIDLQ&   operator=(const AliTRDCalPIDLQ &c);
  void      Init();

  ClassDef(AliTRDCalPIDLQ, 2)                 // LQ PID reference manager
};

//_________________________________________________________________________
inline Int_t AliTRDCalPIDLQ::GetModelID(Int_t mom, Int_t spec, Int_t /*ii*/)
{  
// Returns the ID of the PDF distribution 
// 5 species * 11 momentum ranges

  return spec * AliTRDCalPID::kNMom + mom;
}

//_________________________________________________________________________
inline Bool_t  AliTRDCalPIDLQ::CookdEdx(const Float_t * const dedx, Double_t *x)
{
// Convert NN dEdx slices to the representation used by LQ

  x[0]=0.;x[1]=0.;
  for(Int_t islice=AliTRDCalPID::kNSlicesNN; islice--;){
    Int_t jslice = islice>kNN2LQtransition;
    x[jslice]+=dedx[islice];
  }
  
  // check data integrity
  if(x[0]<1.e-30) return kFALSE;
  if(x[1]<1.e-30) return kFALSE;

  x[0] = TMath::Log(x[0]);
  x[1] = TMath::Log(x[1]);
  return kTRUE;
}

#endif

