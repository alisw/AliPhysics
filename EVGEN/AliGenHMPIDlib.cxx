/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// Library class for particle pt and y distributions used for 
// HMPID simulations.
// To be used with AliGenParam.
// The following particle typed can be simulated:
// phi, lambda, k
//
// Author: Annalisa Mastroserio <Annalisa.Mastroserio@ba.infn.it>
//
//

#include <TPDGCode.h>
#include <TMath.h>
#include <TRandom.h>
#include <TString.h>

#include "AliGenHMPIDlib.h"
#include "AliLog.h"

ClassImp(AliGenHMPIDlib)

  //---------------------------------------
  //    Pi Plus
  //---------------------------------------
Int_t AliGenHMPIDlib::IpPiPlus(TRandom *)
{
//PDG code
  return 211;     
}

Double_t AliGenHMPIDlib::PtPiPlusFlat( Double_t *, Double_t *)
{
  //PiPlus FLAT pt-distribution
  return 1; 
}

Double_t AliGenHMPIDlib::PtPiPlusExp( Double_t *x, Double_t *)
{
  //PiPlus     EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);   
}

Double_t AliGenHMPIDlib::YPiPlusFlat(Double_t *,Double_t *)
{
  //PiPlus            y-distribution 
 return 1;
}


//---------------------------------------
//        Pi Minus
//---------------------------------------
Int_t AliGenHMPIDlib::IpPiMinus(TRandom *)
{
//PDG code
  return -211;     
}

Double_t AliGenHMPIDlib::PtPiMinusFlat( Double_t *, Double_t *)
{
// PiMinus FLAT pt-distribution
  return 1; 
}

Double_t AliGenHMPIDlib::PtPiMinusExp( Double_t *x, Double_t *)
{
//PiMinus     EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);   
}

Double_t AliGenHMPIDlib::YPiMinusFlat(Double_t *,Double_t *)
{
//PiMinus          y-distribution 
  return 1;
}


//--------------------------------------------
//        K Plus
//--------------------------------------------
Int_t AliGenHMPIDlib::IpKPlus(TRandom *)
{
//PDG code
  return 321;
}

Double_t AliGenHMPIDlib::PtKPlusFlat( Double_t *, Double_t *)
{
// K+ FLAT pt-distribution
  return 1;
}

Double_t AliGenHMPIDlib::PtKPlusExp( Double_t *x, Double_t *)
{
// K+   EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);
}

Double_t AliGenHMPIDlib::YKPlusFlat(Double_t *,Double_t *)
{
// K+             y-distribution
  return 1;
}


//-----------------------------------------------
//         K Minus
//-----------------------------------------------
Int_t AliGenHMPIDlib::IpKMinus(TRandom *)
{
//PDG code
  return -321;
}

Double_t AliGenHMPIDlib::PtKMinusFlat( Double_t *, Double_t *)
{
// K- FLAT pt-distribution
  return 1;
}

Double_t AliGenHMPIDlib::PtKMinusExp( Double_t *x, Double_t *)
{
// K-   EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);
}

Double_t AliGenHMPIDlib::YKMinusFlat(Double_t *,Double_t *)
{
// K-             y-distribution
  return 1;
}


//-----------------------------------------------
//       K0 short
//-----------------------------------------------
Int_t AliGenHMPIDlib::IpK0s(TRandom *)
{
//PDG code
  return 310;
}

Double_t AliGenHMPIDlib::PtK0sFlat( Double_t *, Double_t *)
{
// K0s FLAT pt-distribution
  return 1;
}

Double_t AliGenHMPIDlib::PtK0sExp( Double_t *x, Double_t *)
{
// K0s   EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);
}

Double_t AliGenHMPIDlib::YK0sFlat(Double_t *,Double_t *)
{
// K0s             y-distribution
  return 1;
}


//---------------------------------------------
//         Phi(1020)
//---------------------------------------------
Int_t AliGenHMPIDlib::IpPhi(TRandom *)
{
//PDG code
  return 333;     
}

Double_t AliGenHMPIDlib::PtPhiFlat( Double_t *, Double_t *)
{
// Phi FLAT pt-distribution
  return 1; 
}

Double_t AliGenHMPIDlib::PtPhiExp( Double_t *x, Double_t *)
{
//phi     EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);   
}

Double_t AliGenHMPIDlib::YPhiFlat(Double_t *,Double_t *)
{
//phi             y-distribution 
  return 1;
}


//-------------------------------------------------------
//                    PROTON
//-------------------------------------------------------
Int_t AliGenHMPIDlib::IpProton(TRandom *)
{
//PDG code
  return 2122;     
}

Double_t AliGenHMPIDlib::PtProtonFlat( Double_t *, Double_t *)
{
// ProtonFLAT pt-distribution

  return 1; 
}

Double_t AliGenHMPIDlib::PtProtonExp( Double_t *x, Double_t *)
{
//Proton    EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);   
}

Double_t AliGenHMPIDlib::YProtonFlat(Double_t *,Double_t *)
{
  //Proton            y-distribution 
  return 1;
}


//-------------------------------------------------------
//                    PROTON-BAR
//-------------------------------------------------------
Int_t AliGenHMPIDlib::IpProtonBar(TRandom *)
{
//PDG code
  return -2122;     
}

Double_t AliGenHMPIDlib::PtProtonBarFlat( Double_t *, Double_t *)
{
// ProtonBar FLAT pt-distribution

  return 1; 
}

Double_t AliGenHMPIDlib::PtProtonBarExp( Double_t *x, Double_t *)
{
//ProtonBar     EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);   
}

Double_t AliGenHMPIDlib::YProtonBarFlat(Double_t *,Double_t *)
{
  //ProtonBar             y-distribution 
  return 1;
}


//-------------------------------------------------------
//                    LAMBDA
//-------------------------------------------------------
Int_t AliGenHMPIDlib::IpLambda(TRandom *)
{
//PDG code
  return 3122;     
}

Double_t AliGenHMPIDlib::PtLambdaFlat( Double_t *, Double_t *)
{
// Lambda FLAT pt-distribution

  return 1; 
}

Double_t AliGenHMPIDlib::PtLambdaExp( Double_t *x, Double_t *)
{
//Lambda     EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);   
}

Double_t AliGenHMPIDlib::YLambdaFlat(Double_t *,Double_t *)
{
  //Lambda             y-distribution 
  return 1;
}


//-------------------------------------------------------
//                    LAMBDA-BAR
//-------------------------------------------------------
Int_t AliGenHMPIDlib::IpLambdaBar(TRandom *)
{
//PDG code
  return -3122;     
}

Double_t AliGenHMPIDlib::PtLambdaBarFlat( Double_t *, Double_t *)
{
// LambdaBar FLAT pt-distribution

  return 1; 
}

Double_t AliGenHMPIDlib::PtLambdaBarExp( Double_t *x, Double_t *)
{
//LambdaBar     EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);   
}

Double_t AliGenHMPIDlib::YLambdaBarFlat(Double_t *,Double_t *)
{
  //LambdaBar             y-distribution 
  return 1;
}





typedef Double_t (*GenFunc)   (Double_t*,  Double_t*);
typedef Int_t    (*GenFuncIp) (TRandom *);

GenFunc AliGenHMPIDlib::GetPt(Int_t iPID, const char * sForm) const
{
// Return pointer to Pt parameterisation
  AliDebug(1,Form("PID: %i, form: %s",iPID,sForm));	
   TString type(sForm);

   switch(iPID) {

   case kPiPlus:  
      if     (type=="FLAT")                                         return PtPiPlusFlat;
      else if(type=="EXP")                                          return PtPiPlusExp;
      else {
        AliFatal(Form("Unknown Pt distribution form: %s",sForm));   return 0;
       }

   case kPiMinus:  
      if     (type=="FLAT")                                         return PtPiMinusFlat;
      else if(type=="EXP")                                          return PtPiMinusExp;
      else {
        AliFatal(Form("Unknown Pt distribution form: %s",sForm));   return 0;
       }

  case kKPlus:
      if     (type=="FLAT")                                         return PtKPlusFlat;
      else if(type=="EXP")                                          return PtKPlusExp;
      else {
        AliFatal(Form("Unknown Pt distribution form: %s",sForm));   return 0;
       }

   case kKMinus:
      if     (type=="FLAT")                                         return PtKMinusFlat;
      else if(type=="EXP")                                          return PtKMinusExp;
      else {
        AliFatal(Form("Unknown Pt distribution form: %s",sForm));   return 0;
       }


  case kK0Short:  
      if     (type=="FLAT")                                         return PtK0sFlat;
      else if(type=="EXP")                                          return PtK0sExp;
      else {
        AliFatal(Form("Unknown Pt distribution form: %s",sForm));   return 0;
       }


   case kPhi:  
      if     (type=="FLAT")                                         return PtPhiFlat;
      else if(type=="EXP")                                          return PtPhiExp;
      else {
        AliFatal(Form("Unknown Pt distribution form: %s",sForm));   return 0;
       }

   case kProton:  
      if     (type=="FLAT")                                         return PtProtonFlat;
      else if(type=="EXP")                                          return PtProtonExp;
      else {
        AliFatal(Form("Unknown Pt distribution form: %s",sForm));   return 0;
       }  

   case kProtonBar:  
      if     (type=="FLAT")                                         return PtProtonBarFlat;
      else if(type=="EXP")                                          return PtProtonBarExp;
      else {
        AliFatal(Form("Unknown Pt distribution form: %s",sForm));   return 0;
       }  

   case kLambda0:  
      if     (type=="FLAT")                                         return PtLambdaFlat;
      else if(type=="EXP")                                          return PtLambdaExp;
      else {
        AliFatal(Form("Unknown Pt distribution form: %s",sForm));   return 0;
       }  

   case kLambda0Bar:  
      if     (type=="FLAT")                                         return PtLambdaBarFlat;
      else if(type=="EXP")                                          return PtLambdaBarExp;
      else {
        AliFatal(Form("Unknown Pt distribution form: %s",sForm));   return 0;
       }  

   default : AliFatal(Form("Unknown particle type: %i",iPID));      return 0;
   }//switch
}

GenFunc AliGenHMPIDlib::GetY(Int_t iPID, const char *sForm) const
{
  AliDebug(1,Form("PID: %i, form: %s",iPID,sForm));	

   switch (iPID) {

   case kPiPlus:                                                   return YPiPlusFlat;
   case kPiMinus:                                                  return YPiMinusFlat;
   case kKPlus:                                                    return YKPlusFlat;
   case kKMinus:                                                   return YKMinusFlat;
   case kK0Short:                                                  return YK0sFlat;
   case kPhi:                                                      return YPhiFlat;
   case kProton:                                                   return YProtonFlat;
   case kProtonBar:                                                return YProtonBarFlat; 
   case kLambda0:                                                  return YLambdaFlat;
   case kLambda0Bar:                                               return YLambdaBarFlat;

  default  : AliFatal(Form("Unknown particle type: %i",iPID));     return 0;

   }//switch
}

GenFuncIp AliGenHMPIDlib::GetIp(Int_t iPID, const char *sForm) const
{
// Return pointer to particle type parameterisation
  AliDebug(1,Form("PID: %i, form: %s",iPID,sForm));   //////////	

  switch (iPID){

    case kPiPlus:                                                return IpPiPlus;
    case kPiMinus:                                               return IpPiMinus;
    case kKPlus:                                                 return IpKPlus;
    case kKMinus:                                                return IpKMinus;
    case kK0Short:                                               return IpK0s;
    case kPhi:                                                   return IpPhi;
    case kProton:                                                return IpProton; 
    case kProtonBar:                                             return IpProtonBar;
    case kLambda0:                                               return IpLambda;
    case kLambda0Bar:                                            return IpLambdaBar; 

    default  : AliFatal(Form("Unknown particle type: %i",iPID))  return 0;
  }
}

