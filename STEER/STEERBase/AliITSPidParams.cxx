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
AliITSPidParams::AliITSPidParams(Bool_t isMC):
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
  if (isMC) InitMC();
  else InitData();
}
//______________________________________________________________________
AliITSPidParams::AliITSPidParams(Char_t * name, Bool_t isMC):
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
  if (isMC) InitMC();
  else InitData();
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
  fSDDPionMPV->SetParameters(-0.892291, 0.003630, 1.866484, 78.378179);

  if(fSDDPionLandauWidth) delete fSDDPionLandauWidth;
  fSDDPionLandauWidth=new TFormula("fSDDPionLandauWidth","[0]/(x*x)+[1]");
  fSDDPionLandauWidth->SetParameters(0.080999, 5.917715);

  if(fSDDPionGaussWidth) delete fSDDPionGaussWidth;
  fSDDPionGaussWidth=new TFormula("fSDDPionGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSDDPionGaussWidth->SetParameters(-0.092822, 7.839729);

  if(fSSDPionMPV) delete fSSDPionMPV;
  fSSDPionMPV=new TFormula("fSSDPionMPV","[0]/(x*x)*TMath::Log(x)+[1]/(x*x*x*x)*TMath::Log(x)+[2]*TMath::Log(x)+[3]");
  fSSDPionMPV->SetParameters(-0.896507, 0.003173, 2.017155, 80.682567);

  if(fSSDPionLandauWidth) delete fSSDPionLandauWidth;
  fSSDPionLandauWidth=new TFormula("fSSDPionLandauWidth","[0]/(x*x)+[1]");
  fSSDPionLandauWidth->SetParameters(0.087182, 5.843610);

  if(fSSDPionGaussWidth) delete fSSDPionGaussWidth;
  fSSDPionGaussWidth=new TFormula("fSSDPionGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSSDPionGaussWidth->SetParameters(-0.110444, 5.837737);

  // kaons
  if(fSDDKaonMPV) delete fSDDKaonMPV;
  fSDDKaonMPV=new TFormula("fSDDKaonMPV","[0]/(x*x)+[1]/(x*x*x*x)+[2]");
  fSDDKaonMPV->SetParameters(17.581590, -0.120134, 72.550701);

  if(fSDDKaonLandauWidth) delete fSDDKaonLandauWidth;
  fSDDKaonLandauWidth=new TFormula("fSDDKaonLandauWidth","[0]/(x*x)+[1]");
  fSDDKaonLandauWidth->SetParameters(1.271756, 5.778888);

  if(fSDDKaonGaussWidth) delete fSDDKaonGaussWidth;
  fSDDKaonGaussWidth=new TFormula("fSDDKaonGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSDDKaonGaussWidth->SetParameters(-1.650298, 8.322084);

  if(fSSDKaonMPV) delete fSSDKaonMPV;
  fSSDKaonMPV=new TFormula("fSSDKaonMPV","[0]/(x*x)+[1]/(x*x*x*x)+[2]");
  fSSDKaonMPV->SetParameters(16.238778, 0.039318, 75.863719);

  if(fSSDKaonLandauWidth) delete fSSDKaonLandauWidth;
  fSSDKaonLandauWidth=new TFormula("fSSDKaonLandauWidth","[0]/(x*x)+[1]");
  fSSDKaonLandauWidth->SetParameters(1.179541, 5.961353);

  if(fSSDKaonGaussWidth) delete fSSDKaonGaussWidth;
  fSSDKaonGaussWidth=new TFormula("fSSDKaonGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSSDKaonGaussWidth->SetParameters(-2.019126, 6.155977);

  // protons
  if(fSDDProtMPV) delete fSDDProtMPV;
  fSDDProtMPV=new TFormula("fSDDProtMPV","[0]/(x*x)+[1]/(x*x*x*x)+[2]");
  fSDDProtMPV->SetParameters(64.482762, -1.667823, 71.850731);

  if(fSDDProtLandauWidth) delete fSDDProtLandauWidth;
  fSDDProtLandauWidth=new TFormula("fSDDProtLandauWidth","[0]/(x*x)+[1]");
  fSDDProtLandauWidth->SetParameters(6.948997, 3.928018);

  if(fSDDProtGaussWidth) delete fSDDProtGaussWidth;
  fSDDProtGaussWidth=new TFormula("fSDDProtGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSDDProtGaussWidth->SetParameters(-6.522760, 12.673959);

  if(fSSDProtMPV) delete fSSDProtMPV;
  fSSDProtMPV=new TFormula("fSSDProtMPV","[0]/(x*x)+[1]/(x*x*x*x)+[2]");
  fSSDProtMPV->SetParameters(63.817375, -1.221779, 73.233644);

  if(fSSDProtLandauWidth) delete fSSDProtLandauWidth;
  fSSDProtLandauWidth=new TFormula("fSSDProtLandauWidth","[0]/(x*x)+[1]");
  fSSDProtLandauWidth->SetParameters(7.286942, 3.581451);

  if(fSSDProtGaussWidth) delete fSSDProtGaussWidth;
  fSSDProtGaussWidth=new TFormula("fSSDProtGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSSDProtGaussWidth->SetParameters(-8.327867, 9.723422);

}
//______________________________________________________________________
void AliITSPidParams::InitData(){
  // initialize TFormulas to Real Data values (=p-p simulations PYTHIA+GEANT)
  // parameter values from LHC10b 
  
  // pions
  if(fSDDPionMPV) delete fSDDPionMPV;
  fSDDPionMPV=new TFormula("fSDDPionMPV","[0]/(x*x)*TMath::Log(x)+[1]/(x*x*x*x)*TMath::Log(x)+[2]*TMath::Log(x)+[3]");
  fSDDPionMPV->SetParameters(-0.907609, 0.006521, 3.340347, 81.297942);

  if(fSDDPionLandauWidth) delete fSDDPionLandauWidth;
  fSDDPionLandauWidth=new TFormula("fSDDPionLandauWidth","[0]/(x*x)+[1]");
  fSDDPionLandauWidth->SetParameters(0.077272, 5.478557);

  if(fSDDPionGaussWidth) delete fSDDPionGaussWidth;
  fSDDPionGaussWidth=new TFormula("fSDDPionGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSDDPionGaussWidth->SetParameters(-0.098529, 10.265711);

  if(fSSDPionMPV) delete fSSDPionMPV;
  fSSDPionMPV=new TFormula("fSSDPionMPV","[0]/(x*x)*TMath::Log(x)+[1]/(x*x*x*x)*TMath::Log(x)+[2]*TMath::Log(x)+[3]");
  fSSDPionMPV->SetParameters(-0.920046, 0.006061, 3.428578, 81.401816);

  if(fSSDPionLandauWidth) delete fSSDPionLandauWidth;
  fSSDPionLandauWidth=new TFormula("fSSDPionLandauWidth","[0]/(x*x)+[1]");
  fSSDPionLandauWidth->SetParameters(0.071243, 5.388830);

  if(fSSDPionGaussWidth) delete fSSDPionGaussWidth;
  fSSDPionGaussWidth=new TFormula("fSSDPionGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSSDPionGaussWidth->SetParameters(-0.099189, 7.412309);

  // kaons
  if(fSDDKaonMPV) delete fSDDKaonMPV;
  fSDDKaonMPV=new TFormula("fSDDKaonMPV","[0]/(x*x)+[1]/(x*x*x*x)+[2]");
  fSDDKaonMPV->SetParameters(15.429146, 0.178251, 74.640293);

  if(fSDDKaonLandauWidth) delete fSDDKaonLandauWidth;
  fSDDKaonLandauWidth=new TFormula("fSDDKaonLandauWidth","[0]/(x*x)+[1]");
  fSDDKaonLandauWidth->SetParameters(0.975202, 5.699311);

  if(fSDDKaonGaussWidth) delete fSDDKaonGaussWidth;
  fSDDKaonGaussWidth=new TFormula("fSDDKaonGaussWidth","[0]/(x*x)+[1]");
  fSDDKaonGaussWidth->SetParameters(1.660840, 9.389192);

  if(fSSDKaonMPV) delete fSSDKaonMPV;
  fSSDKaonMPV=new TFormula("fSSDKaonMPV","[0]/(x*x)+[1]/(x*x*x*x)+[2]");
  fSSDKaonMPV->SetParameters(15.170715, 0.181379, 74.951884);

  if(fSSDKaonLandauWidth) delete fSSDKaonLandauWidth;
  fSSDKaonLandauWidth=new TFormula("fSSDKaonLandauWidth","[0]/(x*x)+[1]");
  fSSDKaonLandauWidth->SetParameters(0.756466, 5.818274);

  if(fSSDKaonGaussWidth) delete fSSDKaonGaussWidth;
  fSSDKaonGaussWidth=new TFormula("fSSDKaonGaussWidth","[0]/(x*x)+[1]");
  fSSDKaonGaussWidth->SetParameters(1.546693, 6.389872);

  // protons
  if(fSDDProtMPV) delete fSDDProtMPV;
  fSDDProtMPV=new TFormula("fSDDProtMPV","[0]/(x*x)+[1]/(x*x*x*x)+[2]");
  fSDDProtMPV->SetParameters(61.452534, 0.372908, 71.668352);

  if(fSDDProtLandauWidth) delete fSDDProtLandauWidth;
  fSDDProtLandauWidth=new TFormula("fSDDProtLandauWidth","[0]/(x*x)+[1]");
  fSDDProtLandauWidth->SetParameters(3.667023, 5.430721);

  if(fSDDProtGaussWidth) delete fSDDProtGaussWidth;
  fSDDProtGaussWidth=new TFormula("fSDDProtGaussWidth","[0]/(x*x)+[1]");
  fSDDProtGaussWidth->SetParameters(5.503814, 9.657439);

  if(fSSDProtMPV) delete fSSDProtMPV;
  fSSDProtMPV=new TFormula("fSSDProtMPV","[0]/(x*x)+[1]/(x*x*x*x)+[2]");
  fSSDProtMPV->SetParameters(60.246538, 0.000323, 71.992031);

  if(fSSDProtLandauWidth) delete fSSDProtLandauWidth;
  fSSDProtLandauWidth=new TFormula("fSSDProtLandauWidth","[0]/(x*x)+[1]");
  fSSDProtLandauWidth->SetParameters(2.568323, 5.939774);

  if(fSSDProtGaussWidth) delete fSSDProtGaussWidth;
  fSSDProtGaussWidth=new TFormula("fSSDProtGaussWidth","[0]/(x*x)+[1]");
  fSSDProtGaussWidth->SetParameters(5.050541, 6.290964);

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
  Double_t step = 0.;
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

