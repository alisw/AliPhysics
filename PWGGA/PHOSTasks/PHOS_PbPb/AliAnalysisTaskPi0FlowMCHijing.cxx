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

// Extension to Pi0FlowMC, mimicing AliPHOSHijingEfficiency
// by Dmitri Peressounko, 05.02.2013
// Authors: Henrik Qvigstad, Dmitri Peressounko
// Date   : 05.04.2013
/* $Id$ */

#include "AliStack.h"
#include "TParticle.h"
#include "AliCaloPhoton.h"

#include "AliAnalysisTaskPi0FlowMCHijing.h"

ClassImp(AliAnalysisTaskPi0FlowMCHijing);

AliAnalysisTaskPi0FlowMCHijing::AliAnalysisTaskPi0FlowMCHijing(const char* name, AliAnalysisTaskPi0Flow::Period period)
: AliAnalysisTaskPi0FlowMC(name, period)
{
}

AliAnalysisTaskPi0FlowMCHijing::~AliAnalysisTaskPi0FlowMCHijing()
{
}


Double_t AliAnalysisTaskPi0FlowMCHijing::PrimaryWeight(Int_t primary)
{
  //Check who is the primary and introduce weight to correct primary spectrum
  
  if(primary<0 || primary>=fStack->GetNtrack())
    return 1 ;
  //trace primaries up to IP
  TParticle* particle =  fStack->Particle(primary);
  Double_t r=particle->R() ;
  Int_t mother = particle->GetFirstMother() ;
  while(mother>-1){
    if(r<1. && particle->GetPdgCode()==111)
      break ;
    particle =  fStack->Particle(mother);
    mother = particle->GetFirstMother() ;
    r=particle->R() ;
  }

  return TMath::Max(0.,PrimaryParticleWeight(particle)) ;
}


Double_t AliAnalysisTaskPi0FlowMCHijing::PrimaryParticleWeight(TParticle* particle)
{
  Int_t pdg = particle->GetPdgCode() ;
  Int_t type=0 ;
  if(pdg == 111 || TMath::Abs(pdg)==211){
    type =1 ;
  }
  else{
    if(TMath::Abs(pdg)<1000){ //Kaon-like
      type =2 ;    
    }
    else
      type = 3;  //baryons
  }
    
  Double_t pt = particle->Pt() ;
  if(type==1){
   if(fCentBin==0) //0-5
     return (1.662990+1.140890*pt-0.192088*pt*pt)/(1.-0.806630*pt+0.304771*pt*pt)+0.141690*pt ;
   if(fCentBin==1) //5-10
     return (1.474351+0.791492*pt-0.066369*pt*pt)/(1.-0.839338*pt+0.317312*pt*pt)+0.093289*pt ;
   if(fCentBin==2) //10-20
     return (1.174728+0.959681*pt-0.137695*pt*pt)/(1.-0.788873*pt+0.299538*pt*pt)+0.128759*pt ; 
   if(fCentBin==3) //20-40
     return (0.927335+0.475349*pt+0.004364*pt*pt)/(1.-0.817966*pt+0.309787*pt*pt)+0.086899*pt ; 
   if(fCentBin==4) //40-60
     return (0.676878+0.190680*pt+0.077031*pt*pt)/(1.-0.790623*pt+0.305183*pt*pt)+0.064510*pt ; 
   if(fCentBin==5) //60-80
     return (0.684726-0.606262*pt+0.409963*pt*pt)/(1.-1.080061*pt+0.456933*pt*pt)+0.005151*pt ; 
  }
  if(type==2){
   if(fCentBin==0) //0-5
     return (-0.417131+2.253936*pt-0.337731*pt*pt)/(1.-0.909892*pt+0.316820*pt*pt)+0.157312*pt ;
   if(fCentBin==1) //5-10
     return (-0.352275+1.844466*pt-0.248598*pt*pt)/(1.-0.897048*pt+0.316462*pt*pt)+0.132461*pt ; 
   if(fCentBin==2) //10-20
     return (-0.475481+1.975108*pt-0.336013*pt*pt)/(1.-0.801028*pt+0.276705*pt*pt)+0.188164*pt ; 
   if(fCentBin==3) //20-40
     return (-0.198954+1.068789*pt-0.103540*pt*pt)/(1.-0.848354*pt+0.299209*pt*pt)+0.112939*pt ; 
   if(fCentBin==4) //40-60
     return (-0.111052+0.664041*pt-0.019717*pt*pt)/(1.-0.804916*pt+0.300779*pt*pt)+0.095784*pt ;
   if(fCentBin==5) //0-5
     return (0.202788-0.439832*pt+0.564585*pt*pt)/(1.-1.254029*pt+0.679444*pt*pt)+0.016235*pt ;
  }
  if(type==3){
   if(fCentBin==0) //0-5
     return (-1.312732+2.743568*pt-0.375775*pt*pt)/(1.-0.717533*pt+0.164694*pt*pt)+0.164445*pt ;
   if(fCentBin==1) //5-10
     return (-1.229425+2.585889*pt-0.330164*pt*pt)/(1.-0.715892*pt+0.167386*pt*pt)+0.133085*pt ; 
   if(fCentBin==2) //10-20
     return (-1.135677+2.397489*pt-0.320355*pt*pt)/(1.-0.709312*pt+0.164350*pt*pt)+0.146095*pt ; 
   if(fCentBin==3) //20-40
     return (-0.889993+1.928263*pt-0.220785*pt*pt)/(1.-0.715991*pt+0.174729*pt*pt)+0.095098*pt ; 
   if(fCentBin==4) //40-60
     return (-0.539237+1.329118*pt-0.115439*pt*pt)/(1.-0.722906*pt+0.186832*pt*pt)+0.059267*pt ; 
   if(fCentBin==5) //60-80
     return (-0.518126+1.327628*pt-0.130881*pt*pt)/(1.-0.665649*pt+0.184300*pt*pt)+0.081701*pt ;   
  }
  return 1. ;  
}