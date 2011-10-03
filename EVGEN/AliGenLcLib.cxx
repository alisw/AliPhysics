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


// Library class for particle pt and y distributions used for 
// LambdaC simulations.
// To be used with AliGenParam.
//
// Author: Annalisa Mastroserio <Annalisa.Mastroserio@cern.ch>
//

#include <TPDGCode.h>
#include <TMath.h>
#include <TRandom.h>
#include <TString.h>

#include "AliGenLcLib.h"
#include "AliLog.h"

ClassImp(AliGenLcLib)


//---------------------------------------------
//        LambdaC 
//---------------------------------------------
  Int_t AliGenLcLib::IpLcPlus(TRandom *)
{
  //PDG code
  return 4122;
}

Int_t AliGenLcLib::IpLcMinus(TRandom *)
{
  //PDG code
  return -4122;
}

Double_t AliGenLcLib::PtLcFlat( const Double_t *, const Double_t *)
{
  // FLAT pt-distribution
  return 1; 
}

Double_t AliGenLcLib::PtLcExp( const Double_t *x, const Double_t *)
{
  //     EXP  pt-distribution
  return x[0]*TMath::Exp(-x[0]/0.17);   
}

Double_t AliGenLcLib::YLcFlat(const Double_t */*x*/,const Double_t *)
{
  //LambdaC             y-distribution 
  return 5;
}



typedef Double_t (*GenFunc)   (const Double_t*,  const Double_t*);
typedef Int_t    (*GenFuncIp) (TRandom *);

GenFunc AliGenLcLib::GetPt(Int_t iPID, const char * sForm) const
{
  // Return pointer to Pt parameterisation
  printf("PID: %i, form: %s \n",iPID,sForm);	
  TString type(sForm);
  GenFunc func;

  switch(iPID) {

  case kLcPlus:  
    if     (type=="FLAT")                                         {func=PtLcFlat; break;}
    else if(type=="EXP")                                          {func=PtLcExp; break;}
    else {
      AliFatal(Form("Unknown Pt distribution form: %s",sForm));   func=0;
    }

  case kLcMinus:  
    if     (type=="FLAT")                                         {func=PtLcFlat; break;}
    else if(type=="EXP")                                          {func=PtLcExp; break;}
    else {
      AliFatal(Form("Unknown Pt distribution form: %s",sForm));   func=0;
    }

  default : AliFatal(Form("Unknown particle type: %i",iPID));      func=0;
  }//switch

  return func;
}

GenFunc AliGenLcLib::GetY(Int_t iPID, const char *sForm) const
{
  AliDebug(1,Form("PID: %i, form: %s",iPID,sForm));	
  GenFunc func;
  switch (iPID) {

  case kLcPlus:                                                func=YLcFlat; break;
  case kLcMinus:                                               func=YLcFlat; break;

  default  : AliFatal(Form("Unknown particle type: %i",iPID));     func=0; break;

  }//switch
  return func;
}

GenFuncIp AliGenLcLib::GetIp(Int_t iPID, const char *sForm) const
{
  // Return pointer to particle type parameterisation
  AliDebug(1,Form("PID: %i, form: %s",iPID,sForm));   //////////	

  switch (iPID){

  case kLcPlus:                                                  return IpLcPlus;
  case kLcMinus:                                                 return IpLcMinus;

  default  : AliFatal(Form("Unknown particle type: %i",iPID));  return 0;
  }
}
