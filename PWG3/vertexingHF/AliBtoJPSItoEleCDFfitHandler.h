#ifndef ALIBTOJPSITOELECDFFITHANDLER_H
#define ALIBTOJPSITOELECDFFITHANDLER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//___________________________________________________________________________
//                  Class AliBtoJPSItoEleCDFfitHandler
//            Class to perform unbinned log-likelihood fit
//      
//                         Origin: C. Di Giglio
//     Contact: Carmelo.Digiglio@ba.infn.it; Giuseppe.Bruno@ba.infn.it
//____________________________________________________________________________

#include <TNamed.h>
#include "TMath.h"
#include "AliBtoJPSItoEleCDFfitFCN.h"
//#include "AliLog.h"

class AliBtoJPSItoEleCDFfitHandler : public TNamed {
 public:
  //
  AliBtoJPSItoEleCDFfitHandler();
  AliBtoJPSItoEleCDFfitHandler& operator= (const AliBtoJPSItoEleCDFfitHandler& c);
  AliBtoJPSItoEleCDFfitHandler(const AliBtoJPSItoEleCDFfitHandler& c);
  AliBtoJPSItoEleCDFfitHandler(Double_t* decaytime, Double_t* invariantmass, Int_t ncand);
  ~AliBtoJPSItoEleCDFfitHandler(); 
  Double_t Up() const { return fUp*TMath::Sqrt(fLikely->GetIntegral()); }
  void SetErrorDef(Double_t up) {fUp = up;}

  Double_t operator()(const Double_t* par) const ;
  void CdfFCN(Int_t & /* npar */, Double_t * /* gin */, Double_t &f, Double_t *par, Int_t /* iflag */);

  Double_t* Decaytime() const         { return fX; }
  Double_t* InvariantMass() const     { return fM; }
  AliBtoJPSItoEleCDFfitFCN* LikelihoodPointer() const { return fLikely; }
  Int_t DoMinimization();

 private:
  //
  Double_t fUp;                                      //error definition 
  Double_t* fX; 	                     	     //pseudo-proper decay time X
  Double_t* fM;                                      //invariant mass M
  AliBtoJPSItoEleCDFfitFCN* fLikely;                 //Log likelihood function
  Int_t fNcand;                                      //number of candidates
  //
  ClassDef(AliBtoJPSItoEleCDFfitHandler,0);

}; 
#endif
