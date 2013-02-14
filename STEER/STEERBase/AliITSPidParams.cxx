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

ClassImp(AliITSPidParams)

//______________________________________________________________________
AliITSPidParams::AliITSPidParams(Bool_t isMC):
TNamed("default",""),
  fSDDElecLandauWidth(0),
  fSDDElecGaussWidth(0),
  fSSDElecLandauWidth(0),
  fSSDElecGaussWidth(0),
  fSDDPionLandauWidth(0),
  fSDDPionGaussWidth(0),
  fSSDPionLandauWidth(0),
  fSSDPionGaussWidth(0),
  fSDDKaonLandauWidth(0),
  fSDDKaonGaussWidth(0),
  fSSDKaonLandauWidth(0),
  fSSDKaonGaussWidth(0),
  fSDDProtLandauWidth(0),
  fSDDProtGaussWidth(0),
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
  fSDDElecLandauWidth(0),
  fSDDElecGaussWidth(0),
  fSSDElecLandauWidth(0),
  fSSDElecGaussWidth(0),
  fSDDPionLandauWidth(0),
  fSDDPionGaussWidth(0),
  fSSDPionLandauWidth(0),
  fSSDPionGaussWidth(0),
  fSDDKaonLandauWidth(0),
  fSDDKaonGaussWidth(0),
  fSSDKaonLandauWidth(0),
  fSSDKaonGaussWidth(0),
  fSDDProtLandauWidth(0),
  fSDDProtGaussWidth(0),
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
  if(fSDDElecLandauWidth) delete fSDDElecLandauWidth;
  if(fSDDElecGaussWidth) delete fSDDElecGaussWidth;
  
  if(fSSDElecLandauWidth) delete fSSDElecLandauWidth;
  if(fSSDElecGaussWidth) delete fSSDElecGaussWidth;

  if(fSDDPionLandauWidth) delete fSDDPionLandauWidth;
  if(fSDDPionGaussWidth) delete fSDDPionGaussWidth;

  if(fSSDPionLandauWidth) delete fSSDPionLandauWidth;
  if(fSSDPionGaussWidth) delete fSSDPionGaussWidth;
  
  if(fSDDKaonLandauWidth) delete fSDDKaonLandauWidth;
  if(fSDDKaonGaussWidth) delete fSDDKaonGaussWidth;
  
  if(fSSDKaonLandauWidth) delete fSSDKaonLandauWidth;
  if(fSSDKaonGaussWidth) delete fSSDKaonGaussWidth;

  if(fSDDProtLandauWidth) delete fSDDProtLandauWidth;
  if(fSDDProtGaussWidth) delete fSDDProtGaussWidth;

  if(fSSDProtLandauWidth) delete fSSDProtLandauWidth;
  if(fSSDProtGaussWidth) delete fSSDProtGaussWidth;
}
//______________________________________________________________________
Double_t AliITSPidParams::BetheBloch(Double_t mom, Double_t mass, const Double_t *p) const{
 //
  // This is the empirical parameterization of the Bethe-Bloch formula.
  // It is normalized to 1 at the minimum.
  //
  // bg - beta*gamma
  // 
  Double_t bg = mom/mass;
  Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);
  Double_t aa = TMath::Power(beta,p[3]);
  Double_t bb = TMath::Power(1./bg,p[4]);
  bb=TMath::Log(p[2]+bb);
  return (p[0]-aa-bb)*p[1]/aa;
}
//______________________________________________________________________
Double_t AliITSPidParams::ExtrapolateWidth(Double_t mom, Double_t x1, Double_t y1,Double_t x2, Double_t y2) const{
  //
  // This is a linear extrapolation of Landau width and Gaussian width 
  // for low momentum.
  
 Double_t slope = (y2-y1)/(x2-x1);
 return slope*mom+(y1-slope*x1);
}//______________________________________________________________________
void AliITSPidParams::InitMC(){
  // initialize TFormulas to Monte Carlo values (=p-p simulations PYTHIA+GEANT)
  // parameter values from LHC10d1 
  // MPV BetheBloch parameters;
  
  //sdd MC electrons parameters
  fSDDElecMPVBetheParams[0] = -0.0931934;
  fSDDElecMPVBetheParams[1] = 77.8422;
  fSDDElecMPVBetheParams[2] = -0.889085;
  fSDDElecMPVBetheParams[3] = -154.455;
  fSDDElecMPVBetheParams[4] = -0.000256076;
  
  //ssd MC electrons parameters
  fSSDElecMPVBetheParams[0] = -0.0989358;
  fSSDElecMPVBetheParams[1] = 77.8271;
  fSSDElecMPVBetheParams[2] = -0.900887;
  fSSDElecMPVBetheParams[3] =-1241.14;
  fSSDElecMPVBetheParams[4] = -0.0014204;
  
   // electrons 
  if(fSDDElecLandauWidth) delete fSDDElecLandauWidth;
  fSDDElecLandauWidth=new TFormula("fSDDElecLandauWidth","[0]/(x*x)+[1]");
  fSDDElecLandauWidth->SetParameters(-0.002702, 6.224960);

  if(fSDDElecGaussWidth) delete fSDDElecGaussWidth;
  fSDDElecGaussWidth=new TFormula("fSDDElecGaussWidth","[0]/(x*x)+[1]");
  fSDDElecGaussWidth->SetParameters(0.012402, 7.993106);

  if(fSSDElecLandauWidth) delete fSSDElecLandauWidth;
  fSSDElecLandauWidth=new TFormula("fSSDElecLandauWidth","[0]/(x*x)+[1]");
  fSSDElecLandauWidth->SetParameters(-0.002144, 6.231089);

  if(fSSDElecGaussWidth) delete fSSDElecGaussWidth;
  fSSDElecGaussWidth=new TFormula("fSSDElecGaussWidth","[0]/(x*x)+[1]");
  fSSDElecGaussWidth->SetParameters(0.014530, 6.217153);
  
  //sdd MC hadrons parameters
  fSDDHadronMPVBetheParams[0] = 1.13531;
  fSDDHadronMPVBetheParams[1] = -156.651;
  fSDDHadronMPVBetheParams[2] = 1.87562;
  fSDDHadronMPVBetheParams[3] = 0.45819;
  fSDDHadronMPVBetheParams[4] = 2.26506;
  
  //ssd MC hadrons parameters
  fSSDHadronMPVBetheParams[0] = -0.451908;
  fSSDHadronMPVBetheParams[1] = -55.4368;
  fSSDHadronMPVBetheParams[2] = 0.984636;
  fSSDHadronMPVBetheParams[3] = 0.97078;
  fSSDHadronMPVBetheParams[4] = 2.50883;

  // pions 
  if(fSDDPionLandauWidth) delete fSDDPionLandauWidth;
  fSDDPionLandauWidth=new TFormula("fSDDPionLandauWidth","[0]/(x*x)+[1]");
  fSDDPionLandauWidth->SetParameters(0.08026, 5.87922);

  if(fSDDPionGaussWidth) delete fSDDPionGaussWidth;
  fSDDPionGaussWidth=new TFormula("fSDDPionGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSDDPionGaussWidth->SetParameters(-0.090495, 7.705286);

  if(fSSDPionLandauWidth) delete fSSDPionLandauWidth;
  fSSDPionLandauWidth=new TFormula("fSSDPionLandauWidth","[0]/(x*x)+[1]");
  fSSDPionLandauWidth->SetParameters(0.083882, 5.823419);

  if(fSSDPionGaussWidth) delete fSSDPionGaussWidth;
  fSSDPionGaussWidth=new TFormula("fSSDPionGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSSDPionGaussWidth->SetParameters(-0.105218, 5.650956);

  // kaons
  if(fSDDKaonLandauWidth) delete fSDDKaonLandauWidth;
  fSDDKaonLandauWidth=new TFormula("fSDDKaonLandauWidth","[0]/(x*x)+[1]");
  fSDDKaonLandauWidth->SetParameters(1.430692, 5.581389);

  if(fSDDKaonGaussWidth) delete fSDDKaonGaussWidth;
  fSDDKaonGaussWidth=new TFormula("fSDDKaonGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSDDKaonGaussWidth->SetParameters(-1.6537, 8.071832);

  if(fSSDKaonLandauWidth) delete fSSDKaonLandauWidth;
  fSSDKaonLandauWidth=new TFormula("fSSDKaonLandauWidth","[0]/(x*x)+[1]");
  fSSDKaonLandauWidth->SetParameters(1.368824, 5.639291);

  if(fSSDKaonGaussWidth) delete fSSDKaonGaussWidth;
  fSSDKaonGaussWidth=new TFormula("fSSDKaonGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSSDKaonGaussWidth->SetParameters(-1.901858, 6.0932281);

  // protons
  if(fSDDProtLandauWidth) delete fSDDProtLandauWidth;
  fSDDProtLandauWidth=new TFormula("fSDDProtLandauWidth","[0]/(x*x)+[1]");
  fSDDProtLandauWidth->SetParameters(6.529418, 5.049098);

  if(fSDDProtGaussWidth) delete fSDDProtGaussWidth;
  fSDDProtGaussWidth=new TFormula("fSDDProtGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSDDProtGaussWidth->SetParameters(-9.360599, 11.318026);

  if(fSSDProtLandauWidth) delete fSSDProtLandauWidth;
  fSSDProtLandauWidth=new TFormula("fSSDProtLandauWidth","[0]/(x*x)+[1]");
  fSSDProtLandauWidth->SetParameters(6.419493, 4.925070);

  if(fSSDProtGaussWidth) delete fSSDProtGaussWidth;
  fSSDProtGaussWidth=new TFormula("fSSDProtGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSSDProtGaussWidth->SetParameters(-9.599064, 9.358656);
}
//______________________________________________________________________
void AliITSPidParams::InitData(){
  // initialize TFormulas to Real Data values (=p-p ALICE Experiment)
  // parameter values from LHC10b 
  // MPV BetheBloch parameters;
 
 //sdd data electrons parameters
  fSDDElecMPVBetheParams[0] = -0.130204;
  fSDDElecMPVBetheParams[1] = 75.1267;
  fSDDElecMPVBetheParams[2] = -0.889337;
  fSDDElecMPVBetheParams[3] = 388.372;
  fSDDElecMPVBetheParams[4] = 0.00134649;
  
  //ssd data electrons parameters
  fSSDElecMPVBetheParams[0] = -0.162773;
  fSSDElecMPVBetheParams[1] = 72.9393;
  fSSDElecMPVBetheParams[2] = -0.896944;
  fSSDElecMPVBetheParams[3] = 3233.02;
  fSSDElecMPVBetheParams[4] = 0.00146896;
  
   // electrons 
  if(fSDDElecLandauWidth) delete fSDDElecLandauWidth;
  fSDDElecLandauWidth=new TFormula("fSDDElecLandauWidth","[0]/(x*x)+[1]");
  fSDDElecLandauWidth->SetParameters(-0.002702, 6.224960);

  if(fSDDElecGaussWidth) delete fSDDElecGaussWidth;
  fSDDElecGaussWidth=new TFormula("fSDDElecGaussWidth","[0]/(x*x)+[1]");
  fSDDElecGaussWidth->SetParameters(0.012402, 7.993106);

  if(fSSDElecLandauWidth) delete fSSDElecLandauWidth;
  fSSDElecLandauWidth=new TFormula("fSSDElecLandauWidth","[0]/(x*x)+[1]");
  fSSDElecLandauWidth->SetParameters(-0.002144, 6.231089);

  if(fSSDElecGaussWidth) delete fSSDElecGaussWidth;
  fSSDElecGaussWidth=new TFormula("fSSDElecGaussWidth","[0]/(x*x)+[1]");
  fSSDElecGaussWidth->SetParameters(0.014530, 6.217153);
  
  //sdd data hadrons parameters
  fSDDHadronMPVBetheParams[0] = -18.1867;
  fSDDHadronMPVBetheParams[1] = -3.45806;
  fSDDHadronMPVBetheParams[2] =  2.23635;
  fSDDHadronMPVBetheParams[3] =  2.08328;
  fSDDHadronMPVBetheParams[4] = -1.92331;
 
  //ssd data hadrons parameters
  fSSDHadronMPVBetheParams[0] = -12.6459;
  fSSDHadronMPVBetheParams[1] = -4.84598;
  fSSDHadronMPVBetheParams[2] =  1.50253;
  fSSDHadronMPVBetheParams[3] =  2.08328;
  fSSDHadronMPVBetheParams[4] = -1.3719;
  
  // pions
  if(fSDDPionLandauWidth) delete fSDDPionLandauWidth;
  fSDDPionLandauWidth=new TFormula("fSDDPionLandauWidth","[0]/(x*x)+[1]");
  fSDDPionLandauWidth->SetParameters(0.07694, 5.468728);

  if(fSDDPionGaussWidth) delete fSDDPionGaussWidth;
  fSDDPionGaussWidth=new TFormula("fSDDPionGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSDDPionGaussWidth->SetParameters(-0.098209, 10.409441);

  if(fSSDPionLandauWidth) delete fSSDPionLandauWidth;
  fSSDPionLandauWidth=new TFormula("fSSDPionLandauWidth","[0]/(x*x)+[1]");
  fSSDPionLandauWidth->SetParameters(0.071602, 5.365442);

  if(fSSDPionGaussWidth) delete fSSDPionGaussWidth;
  fSSDPionGaussWidth=new TFormula("fSSDPionGaussWidth","[0]/(x*x)*TMath::Log(x)+[1]");
  fSSDPionGaussWidth->SetParameters(-0.098045, 7.617583);

  // kaons
  if(fSDDKaonLandauWidth) delete fSDDKaonLandauWidth;
  fSDDKaonLandauWidth=new TFormula("fSDDKaonLandauWidth","[0]/(x*x)+[1]");
  fSDDKaonLandauWidth->SetParameters(0.998191, 5.461668);

  if(fSDDKaonGaussWidth) delete fSDDKaonGaussWidth;
  fSDDKaonGaussWidth=new TFormula("fSDDKaonGaussWidth","[0]/(x*x)+[1]");
  fSDDKaonGaussWidth->SetParameters(1.629308, 9.851873);

  if(fSSDKaonLandauWidth) delete fSSDKaonLandauWidth;
  fSSDKaonLandauWidth=new TFormula("fSSDKaonLandauWidth","[0]/(x*x)+[1]");
  fSSDKaonLandauWidth->SetParameters(0.773113, 5.618683);

  if(fSSDKaonGaussWidth) delete fSSDKaonGaussWidth;
  fSSDKaonGaussWidth=new TFormula("fSSDKaonGaussWidth","[0]/(x*x)+[1]");
  fSSDKaonGaussWidth->SetParameters(1.510713, 6.862774);

  // protons
  if(fSDDProtLandauWidth) delete fSDDProtLandauWidth;
  fSDDProtLandauWidth=new TFormula("fSDDProtLandauWidth","[0]/(x*x)+[1]");
  fSDDProtLandauWidth->SetParameters(3.561429, 5.372105);

  if(fSDDProtGaussWidth) delete fSDDProtGaussWidth;
  fSDDProtGaussWidth=new TFormula("fSDDProtGaussWidth","[0]/(x*x)+[1]");
  fSDDProtGaussWidth->SetParameters(5.395926, 10.044613);

  if(fSSDProtLandauWidth) delete fSSDProtLandauWidth;
  fSSDProtLandauWidth=new TFormula("fSSDProtLandauWidth","[0]/(x*x)+[1]");
  fSSDProtLandauWidth->SetParameters(2.647428, 5.678460);

  if(fSSDProtGaussWidth) delete fSSDProtGaussWidth;
  fSSDProtGaussWidth=new TFormula("fSSDProtGaussWidth","[0]/(x*x)+[1]");
  fSSDProtGaussWidth->SetParameters(4.91025, 6.779763);
}
//_______________________________________________________________________
Double_t AliITSPidParams::GetLandauGausNormPdgCode(Double_t dedx, Int_t pdgCode, Double_t mom, Int_t lay) const {
  // Computes Landau Gauss convolution for given particle specie and given momentum in a given ITS layer
  if(TMath::Abs(pdgCode)==11) return GetLandauGausNorm(dedx,AliPID::kElectron,mom,lay);
  else if(TMath::Abs(pdgCode)==211) return GetLandauGausNorm(dedx,AliPID::kPion,mom,lay);
  else if(TMath::Abs(pdgCode)==321) return GetLandauGausNorm(dedx,AliPID::kKaon,mom,lay);
  else if(TMath::Abs(pdgCode)==2212) return GetLandauGausNorm(dedx,AliPID::kProton,mom,lay);
  else return 0.;
}
//_______________________________________________________________________
Double_t AliITSPidParams::GetLandauGausNorm(Double_t dedx, Int_t partType, Double_t mom, Int_t lay) const{
  // Computes Landau Gauss convolution for given particle specie and given momentum in a given ITS layer
  
  Double_t par[4];
  Bool_t isSet=kFALSE;
  if(partType==AliPID::kElectron){
    if(lay==3 || lay==4){
      par[0]=GetSDDElecLandauWidth(mom);
      par[1]=GetSDDElecMPV(mom);
      par[2]=GetSDDElecGaussWidth(mom);
      isSet=kTRUE;
    }
    else if(lay==5 || lay==6){
      par[0]=GetSSDElecLandauWidth(mom);
      par[1]=GetSSDElecMPV(mom);
      par[2]=GetSSDElecGaussWidth(mom);
      isSet=kTRUE;
    }
  }else if(partType==AliPID::kPion){
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
  const Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  const Double_t mpshift  = -0.22278298;       // Landau maximum location
  // Control constants
  const Double_t np = 100.0;      // number of convolution steps
  const Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
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
  step = (xupp-xlow) / np;
  
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

