/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: */

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class to store parameters of ITS        //
// response funcions for dE/dx based PID                         //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TFormula.h>
#include <TNamed.h>
#include <TMath.h>
#include "AliITSPidParams.h"
#include "AliPID.h"

ClassImp(AliITSPidParams)

//______________________________________________________________________
AliITSPidParams::AliITSPidParams():
  TNamed("default",""),
  fSDDPionMPV(0),
  fSDDPionLandauWidth(0),
  fSDDPionGaussWidth(0),
  fSSDPionMPV(0),
  fSSDPionLandauWidth(0),
  fSSDPionGaussWidth(0),
  fSDDKaonMPV(0),
  fSDDKaonLandauWidth(0),
  fSDDKaonGaussWidth(0),
  fSSDKaonMPV(0),
  fSSDKaonLandauWidth(0),
  fSSDKaonGaussWidth(0),
  fSDDProtMPV(0),
  fSDDProtLandauWidth(0),
  fSDDProtGaussWidth(0),
  fSSDProtMPV(0),
  fSSDProtLandauWidth(0),
  fSSDProtGaussWidth(0)
{
  // default constructor
  InitDefaults();
}
//______________________________________________________________________
AliITSPidParams::AliITSPidParams(Char_t * name):
  TNamed(name,""),
  fSDDPionMPV(0),
  fSDDPionLandauWidth(0),
  fSDDPionGaussWidth(0),
  fSSDPionMPV(0),
  fSSDPionLandauWidth(0),
  fSSDPionGaussWidth(0),
  fSDDKaonMPV(0),
  fSDDKaonLandauWidth(0),
  fSDDKaonGaussWidth(0),
  fSSDKaonMPV(0),
  fSSDKaonLandauWidth(0),
  fSSDKaonGaussWidth(0),
  fSDDProtMPV(0),
  fSDDProtLandauWidth(0),
  fSDDProtGaussWidth(0),
  fSSDProtMPV(0),
  fSSDProtLandauWidth(0),
  fSSDProtGaussWidth(0)
{
  // standard constructor
  InitDefaults();
}
//______________________________________________________________________
AliITSPidParams::~AliITSPidParams(){
  // 
  if(fSDDPionMPV) delete fSDDPionMPV;
  if(fSDDPionLandauWidth) delete fSDDPionLandauWidth;
  if(fSDDPionGaussWidth) delete fSDDPionGaussWidth;

  if(fSSDPionMPV) delete fSSDPionMPV;
  if(fSSDPionLandauWidth) delete fSSDPionLandauWidth;
  if(fSSDPionGaussWidth) delete fSSDPionGaussWidth;

  if(fSDDKaonMPV) delete fSDDKaonMPV;
  if(fSDDKaonLandauWidth) delete fSDDKaonLandauWidth;
  if(fSDDKaonGaussWidth) delete fSDDKaonGaussWidth;

  if(fSSDKaonMPV) delete fSSDKaonMPV;
  if(fSSDKaonLandauWidth) delete fSSDKaonLandauWidth;
  if(fSSDKaonGaussWidth) delete fSSDKaonGaussWidth;

  if(fSDDProtMPV) delete fSDDProtMPV;
  if(fSDDProtLandauWidth) delete fSDDProtLandauWidth;
  if(fSDDProtGaussWidth) delete fSDDProtGaussWidth;

  if(fSSDProtMPV) delete fSSDProtMPV;
  if(fSSDProtLandauWidth) delete fSSDProtLandauWidth;
  if(fSSDProtGaussWidth) delete fSSDProtGaussWidth;
}

//______________________________________________________________________
void AliITSPidParams::InitDefaults(){
  // initialize TFormulas to default values (=p-p simulations PYTHIA+GEANT)
  // parameter values from Emanuele Biolcati

  // pions
  if(fSDDPionMPV) delete fSDDPionMPV;
  fSDDPionMPV=new TFormula("fSDDPionMPV","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]*TMath::Log(x*x)+[3]");
  fSDDPionMPV->SetParameters(1.4831933,-0.000403,3.9722756,92.710680);

  if(fSDDPionLandauWidth) delete fSDDPionLandauWidth;
  fSDDPionLandauWidth=new TFormula("fSDDPionLandauWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSDDPionLandauWidth->SetParameters(-0.091098,-0.001355,8.3019280);

  if(fSDDPionGaussWidth) delete fSDDPionGaussWidth;
  fSDDPionGaussWidth=new TFormula("fSDDPionGaussWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSDDPionGaussWidth->SetParameters(-0.129570,-0.002686,15.287701);

  if(fSSDPionMPV) delete fSSDPionMPV;
  fSSDPionMPV=new TFormula("fSSDPionMPV","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]*TMath::Log(x*x)+[3]");
  fSSDPionMPV->SetParameters(1.2455667,-0.000743,3.2260119,84.237030);

  if(fSSDPionLandauWidth) delete fSSDPionLandauWidth;
  fSSDPionLandauWidth=new TFormula("fSSDPionLandauWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSSDPionLandauWidth->SetParameters(-0.083918,-0.001246,7.4750478);

  if(fSSDPionGaussWidth) delete fSSDPionGaussWidth;
  fSSDPionGaussWidth=new TFormula("fSSDPionGaussWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSSDPionGaussWidth->SetParameters(-0.140441,-0.002984,10.936747);

  // kaons
  if(fSDDKaonMPV) delete fSDDKaonMPV;
  fSDDKaonMPV=new TFormula("fSDDKaonMPV","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSDDKaonMPV->SetParameters(16.998674,-0.058006,89.660709);

  if(fSDDKaonLandauWidth) delete fSDDKaonLandauWidth;
  fSDDKaonLandauWidth=new TFormula("fSDDKaonLandauWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSDDKaonLandauWidth->SetParameters(1.0014235,-5.61e-15,9.2672327);

  if(fSDDKaonGaussWidth) delete fSDDKaonGaussWidth;
  fSDDKaonGaussWidth=new TFormula("fSDDKaonGaussWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSDDKaonGaussWidth->SetParameters(0.7919025,-0.103579,14.016803);

  if(fSSDKaonMPV) delete fSSDKaonMPV;
  fSSDKaonMPV=new TFormula("fSSDKaonMPV","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSSDKaonMPV->SetParameters(14.090845,-0.087253,81.782088);

  if(fSSDKaonLandauWidth) delete fSSDKaonLandauWidth;
  fSSDKaonLandauWidth=new TFormula("fSSDKaonLandauWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSSDKaonLandauWidth->SetParameters(1.0769127,-9.06e-13,7.5221492);

  if(fSSDKaonGaussWidth) delete fSSDKaonGaussWidth;
  fSSDKaonGaussWidth=new TFormula("fSSDKaonGaussWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSSDKaonGaussWidth->SetParameters(0.5155878,-0.098696,10.771975);

  // protons
  if(fSDDProtMPV) delete fSDDProtMPV;
  fSDDProtMPV=new TFormula("fSDDProtMPV","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSDDProtMPV->SetParameters(70.325146,0.0386808,87.797052);

  if(fSDDProtLandauWidth) delete fSDDProtLandauWidth;
  fSDDProtLandauWidth=new TFormula("fSDDProtLandauWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSDDProtLandauWidth->SetParameters(4.0476840,-3.77e-14,10.294707);

  if(fSDDProtGaussWidth) delete fSDDProtGaussWidth;
  fSDDProtGaussWidth=new TFormula("fSDDProtGaussWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSDDProtGaussWidth->SetParameters(5.9780498,-2.16e-11,13.357744);

  if(fSSDProtMPV) delete fSSDProtMPV;
  fSSDProtMPV=new TFormula("fSSDProtMPV","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSSDProtMPV->SetParameters(58.918831,-0.303855,80.341765);
			     
  if(fSSDProtLandauWidth) delete fSSDProtLandauWidth;
  fSSDProtLandauWidth=new TFormula("fSSDProtLandauWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSSDProtLandauWidth->SetParameters(3.0814273,-1.26e-13,8.8353833);

  if(fSSDProtGaussWidth) delete fSSDProtGaussWidth;
  fSSDProtGaussWidth=new TFormula("fSSDProtGaussWidth","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x*x)+[2]");
  fSSDProtGaussWidth->SetParameters(5.6177229,-1.67e-13,10.538921);

}
//_______________________________________________________________________
Double_t AliITSPidParams::GetLandauGausNormPdgCode(Double_t dedx, Int_t pdgCode, Double_t mom, Int_t lay) const {
  // Computes Landau Gauss convolution for given particle specie and given momentum in a given ITS layer
  if(TMath::Abs(pdgCode)==211) return GetLandauGausNorm(dedx,AliPID::kPion,mom,lay);
  else if(TMath::Abs(pdgCode)==321) return GetLandauGausNorm(dedx,AliPID::kKaon,mom,lay);
  else if(TMath::Abs(pdgCode)==2212) return GetLandauGausNorm(dedx,AliPID::kProton,mom,lay);
  else return 0.;
}
//_______________________________________________________________________
  Double_t AliITSPidParams::GetLandauGausNorm(Double_t dedx, Int_t partType, Double_t mom, Int_t lay) const{
  // Computes Landau Gauss convolution for given particle specie and given momentum in a given ITS layer

  Double_t par[3];
  Bool_t isSet=kFALSE;
  if(partType==AliPID::kPion){
    if(lay==3 || lay==4){
      par[0]=GetSDDPionLandauWidth(mom);
      par[1]=GetSDDPionMPV(mom);
      par[2]=GetSDDPionGaussWidth(mom);
      isSet=kTRUE;
    }
    else if(lay==5 || lay==6){
      par[0]=GetSSDPionLandauWidth(mom);
      par[1]=GetSSDPionMPV(mom);
      par[2]=GetSSDPionGaussWidth(mom);
      isSet=kTRUE;
    }
  }else if(partType==AliPID::kKaon){
    if(lay==3 || lay==4){
      par[0]=GetSDDKaonLandauWidth(mom);
      par[1]=GetSDDKaonMPV(mom);
      par[2]=GetSDDKaonGaussWidth(mom);
      isSet=kTRUE;
    }
    else if(lay==5 || lay==6){
      par[0]=GetSSDKaonLandauWidth(mom);
      par[1]=GetSSDKaonMPV(mom);
      par[2]=GetSSDKaonGaussWidth(mom);
      isSet=kTRUE;
    }
  }else if(partType==AliPID::kProton){
    if(lay==3 || lay==4){
      par[0]=GetSDDProtLandauWidth(mom);
      par[1]=GetSDDProtMPV(mom);
      par[2]=GetSDDProtGaussWidth(mom);
      isSet=kTRUE;
    }
    else if(lay==5 || lay==6){
      par[0]=GetSSDProtLandauWidth(mom);
      par[1]=GetSSDProtMPV(mom);
      par[2]=GetSSDProtGaussWidth(mom);
      isSet=kTRUE;
    }
  }
  if(!isSet) return 0.;
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0];
  // Range of convolution integral
  xlow = dedx - sc * par[2];
  xupp = dedx + sc * par[2];
  if(np!=0) step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
   
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(dedx,xx,par[2]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(dedx,xx,par[2]);
  }
  
  return (step * sum * invsq2pi / par[2]);
}

