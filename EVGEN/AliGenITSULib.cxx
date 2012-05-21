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
// ITS Upgrade related signal simulations.
//
// Author: Annalisa Mastroserio <Annalisa.Mastroserio@cern.ch>
//

#include <TPDGCode.h>
#include <TMath.h>
#include <TRandom.h>
#include <TString.h>

#include "AliGenITSULib.h"
#include "AliLog.h"

ClassImp(AliGenITSULib)

 typedef Double_t (*GenFunc)   (const Double_t*,  const Double_t*);
 typedef Int_t    (*GenFuncIp) (TRandom *);


Double_t AliGenITSULib::PtLcDist( const Double_t *x, const Double_t *)
{
 //     pt-distribution
 // same shape as the Lb pt distribution for high pt were overestimated. The Lb fitting funcion has been
 // used from 3.5 GeV/c on. At smaller Pt, the shape is likely as  sqrt(x)*exp(-x) 
 Double_t par[6]={2.781336,1.353902,4.634330,(145479448.743187/202262270.892062),1.789479,-1.329143}; 
 Double_t y=0;
 if(x[0]<3.5) y= par[3]*(TMath::Power(x[0],par[4]))*TMath::Exp(x[0]*par[5]);
 else y= x[0]/TMath::Power((1+TMath::Power(x[0]/par[0],par[1])),par[2]);
 return y;
}


Double_t AliGenITSULib::PtLbDist( const Double_t *x, const Double_t *)
{
 //     pt-distribution (fitted in 0-30 GeV/c)
 Double_t par[4]={3.42500e-03,6.18902,1.76908,3.24823};
 return par[0]*x[0]/TMath::Power((1+TMath::Power(x[0]/par[1],par[2])),par[3]);
}

GenFunc AliGenITSULib::GetPt(Int_t iPID, const char * sForm) const
{
 // Return pointer to Pt parameterisation
 AliDebug(1,Form("PID: %i, form: %s \n",iPID,sForm));	
 TString type(sForm);
 GenFunc func;

 if(type=="FLAT") {
  func=PtFlat;
 } else if(type=="DIST") {

  switch(TMath::Abs(iPID)) {

   case kLb    :    func=PtLbDist; break;  
   case kLc    :    func=PtLcDist; break;  
   case kXi_c  :    func=PtLcDist; break;
   case kB     :    func=PtLbDist; break;
   case kDs    :    func=PtLcDist; break;
   case kDplus :    func=PtLcDist; break;
   default : AliError(Form("Unknown particle type: %i, Pt dist is 0",iPID));      func=0;
  } 
 }else {
  AliError(Form("Unknown Pt distribution %s. Pt distribution is set to 0 ",sForm));
  func=0;
 }

 return func;
}

GenFunc AliGenITSULib::GetY(Int_t iPID, const char *sForm) const
{
 GenFunc func;

 if(TMath::Abs(iPID) != kLc && TMath::Abs(iPID) != kLb && TMath::Abs(iPID) != kXi_c && TMath::Abs(iPID) != kB && TMath::Abs(iPID)!=kDplus && TMath::Abs(iPID)!=kDs) {
  AliError(Form("Unknown PID: %i, form: %s, returning 0",iPID,sForm));   //////////	
  func=0;
 } else { 
  func = YFlat;
 }
 return func;
}

GenFuncIp AliGenITSULib::GetIp(Int_t iPID, const char *sForm) const
{
 AliDebug(1,Form(" %i - %s",iPID,sForm));
 // Return pointer to particle type parameterisation
 GenFuncIp id;

 if(TMath::Abs(iPID) != kLc && TMath::Abs(iPID) != kLb && TMath::Abs(iPID) != kXi_c && TMath::Abs(iPID) != kB && TMath::Abs(iPID)!=kDplus && TMath::Abs(iPID)!=kDs) {
  AliError(Form("Unknown PID: %i, form: %s, return 0",iPID,sForm));   //////////	
  id = 0;
 } else {
  switch (iPID){
   case kLc    :                                  return id=IpLcPlus;
   case -kLc   :                                  return id=IpLcMinus;
   case kLb    :                                  return id=IpLb;
   case -kLb   :                                  return id=IpLbBar;
   case kXi_c  :                                  return id=IpXic;
   case -kXi_c :                                  return id=IpXicBar;
   case kB     :                                  return id=IpBPlus;
   case -kB    :                                  return id=IpBMinus;
   case kDs    :                                  return id=IpDsPlus;
   case -kDs   :                                  return id=IpDsMinus;
   case kDplus :                                  return id=IpDPlus;
   case -kDplus:                                  return id=IpDMinus;
   default  : AliFatal(Form("Unknown particle type: %i",iPID));  id=0;
  }

 }

 return id;
}
