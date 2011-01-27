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
  InitMC();
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
  InitMC();
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
void AliITSPidParams::InitMC(){
  // initialize TFormulas to Monte Carlo values (=p-p simulations PYTHIA+GEANT)
  // parameter values from LHC10d1 

  // pions
  if(fSDDPionMPV) delete fSDDPionMPV;
  fSDDPionMPV=new TFormula("fSDDPionMPV","[0]/(x*x)*TMath::Log(x)+[1]/(x*x*x*x)*TMath::Log(x)+[2]*TMath::Log(x)+[3]");
  fSDDPionMPV->SetParameters(-0.690010, 0.002602, 1.185083, 78.454691);

  if(fSDDPionLandauWidth) delete fSDDPionLandauWidth;
  fSDDPionLandauWidth=new TFormula("fSDDPionLandauWidth","[0]/(x*x)+[1]");
  fSDDPionLandauWidth->SetParameters(0.061606, 5.960376);

  if(fSDDPionGaussWidth) delete fSDDPionGaussWidth;
  fSDDPionGaussWidth=new TFormula("fSDDPionGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSDDPionGaussWidth->SetParameters(-0.065307, 7.896339);

  if(fSSDPionMPV) delete fSSDPionMPV;
  fSSDPionMPV=new TFormula("fSSDPionMPV","[0]/(x*x)*TMath::Log(x)+[1]/(x*x*x*x)*TMath::Log(x)+[2]*TMath::Log(x)+[3]");
  fSSDPionMPV->SetParameters(-0.699466, 0.002429, 1.366895, 80.759188);

  if(fSSDPionLandauWidth) delete fSSDPionLandauWidth;
  fSSDPionLandauWidth=new TFormula("fSSDPionLandauWidth","[0]/(x*x)+[1]");
  fSSDPionLandauWidth->SetParameters(0.066319, 5.889438);

  if(fSSDPionGaussWidth) delete fSSDPionGaussWidth;
  fSSDPionGaussWidth=new TFormula("fSSDPionGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSSDPionGaussWidth->SetParameters(-0.077798, 5.903887);

  // kaons
  if(fSDDKaonMPV) delete fSDDKaonMPV;
  fSDDKaonMPV=new TFormula("fSDDKaonMPV","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x)+[2]");
  fSDDKaonMPV->SetParameters(15.924230, 0.085357, 73.528107);

  if(fSDDKaonLandauWidth) delete fSDDKaonLandauWidth;
  fSDDKaonLandauWidth=new TFormula("fSDDKaonLandauWidth","[0]/(x*x)+[1]");
  fSDDKaonLandauWidth->SetParameters(1.121062, 5.925409);

  if(fSDDKaonGaussWidth) delete fSDDKaonGaussWidth;
  fSDDKaonGaussWidth=new TFormula("fSDDKaonGaussWidth","[0]/(x*x)+[1]");
  fSDDKaonGaussWidth->SetParameters(2.010609, 5.973445);

  if(fSSDKaonMPV) delete fSSDKaonMPV;
  fSSDKaonMPV=new TFormula("fSSDKaonMPV","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x)+[2]");
  fSSDKaonMPV->SetParameters(15.197250, 0.016714, 76.446132);

  if(fSSDKaonLandauWidth) delete fSSDKaonLandauWidth;
  fSSDKaonLandauWidth=new TFormula("fSSDKaonLandauWidth","[0]/(x*x)+[1]");
  fSSDKaonLandauWidth->SetParameters(1.036749, 6.106413);

  if(fSSDKaonGaussWidth) delete fSSDKaonGaussWidth;
  fSSDKaonGaussWidth=new TFormula("fSSDKaonGaussWidth","[0]/(x*x)+[1]");
  fSSDKaonGaussWidth->SetParameters(2.426498, 3.383779);

  // protons
  if(fSDDProtMPV) delete fSDDProtMPV;
  fSDDProtMPV=new TFormula("fSDDProtMPV","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x)+[2]");
  fSDDProtMPV->SetParameters(56.888592, 1.115447, 75.416075);

  if(fSDDProtLandauWidth) delete fSDDProtLandauWidth;
  fSDDProtLandauWidth=new TFormula("fSDDProtLandauWidth","[0]/(x*x)+[1]");
  fSDDProtLandauWidth->SetParameters(6.350805, 4.312568);

  if(fSDDProtGaussWidth) delete fSDDProtGaussWidth;
  fSDDProtGaussWidth=new TFormula("fSDDProtGaussWidth","[0]/(x*x)+[1]");
  fSDDProtGaussWidth->SetParameters(6.556759, 5.953683);

  if(fSSDProtMPV) delete fSSDProtMPV;
  fSSDProtMPV=new TFormula("fSSDProtMPV","[0]/(x*x)+[1]/(x*x*x*x)*TMath::Log(x)+[2]");
  fSSDProtMPV->SetParameters(57.385512, 0.884585, 76.138989);

  if(fSSDProtLandauWidth) delete fSSDProtLandauWidth;
  fSSDProtLandauWidth=new TFormula("fSSDProtLandauWidth","[0]/(x*x)+[1]");
  fSSDProtLandauWidth->SetParameters(6.653282, 3.997930);

  if(fSSDProtGaussWidth) delete fSSDProtGaussWidth;
  fSSDProtGaussWidth=new TFormula("fSSDProtGaussWidth","[0]/(x*x)+[1]");
  fSSDProtGaussWidth->SetParameters(8.203296, 1.491822);

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

