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

// Extension to Pi0FlowMC, using parametrized weights.
// Authors: Boris Polishchuk
// Date   : 09.07.2013

#include "AliStack.h"
#include "TParticle.h"
#include "AliCaloPhoton.h"
#include "TH1.h"
#include "TF1.h"

#include "AliAnalysisTaskPi0FlowMCParamWeights.h"

ClassImp(AliAnalysisTaskPi0FlowMCParamWeights);

AliAnalysisTaskPi0FlowMCParamWeights::AliAnalysisTaskPi0FlowMCParamWeights(const char* name, AliAnalysisTaskPi0Flow::Period period)
  : AliAnalysisTaskPi0FlowMC(name, period)
{}

AliAnalysisTaskPi0FlowMCParamWeights::~AliAnalysisTaskPi0FlowMCParamWeights()
{}

void AliAnalysisTaskPi0FlowMCParamWeights::UserCreateOutputObjects()
{
  // Do Pi0FlowMC CreateOuputObjects and call Sumw2 for newly created histograms.
  AliAnalysisTaskPi0FlowMC::UserCreateOutputObjects();
  
  TH1 * hist = 0;
  char key[80];

  const int kNPID = 14;
  const char* pidNames[kNPID] = {"All", "Allcore", "Disp", "Disp2", "Dispcore",  "Disp2core", "CPV", "CPVcore", "CPV2", "CPV2core", "Both", "Bothcore", "Both2", "Both2core"};
  
  for(UInt_t iBin=0; iBin<GetNumberOfCentralityBins(); iBin++) {
    for(Int_t ipid=0; ipid < kNPID; ipid++){

      sprintf(key,"hPi0%s_cen%d", pidNames[ipid],iBin);  
      hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key));
      if(hist) { hist->Sumw2(); printf("     ->Sumw2 invoked for %s.\n",key); }

      sprintf(key,"hMiPi0%s_cen%d", pidNames[ipid],iBin);  
      hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key));
      if(hist) { hist->Sumw2(); printf("     ->Sumw2 invoked for %s.\n",key); }
    }
  }
  
}

Double_t AliAnalysisTaskPi0FlowMCParamWeights::PrimaryWeight(Int_t primary){
  //Check who is the primary and introduce weight to correct primary spectrum 
 
  if(primary<0 || primary>=fStack->GetNtrack())
    return 1 ;
  //trace primaries up to IP
  TParticle* particle =  fStack->Particle(primary);

  Double_t r=particle->R() ;
  Int_t mother = particle->GetFirstMother() ;
  while(mother>-1){
    if(r<1.)
      break ;
    particle =  fStack->Particle(mother);
    mother = particle->GetFirstMother() ;
    r=particle->R() ;
  }
  
  return TMath::Max(0.,PrimaryParticleWeight(particle)) ;
}

Double_t AliAnalysisTaskPi0FlowMCParamWeights::PrimaryParticleWeight(TParticle * particle){
  //Weight for particle

    
  // pi0 weight (MC/Data):
  // w = (par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0])/(1.+par[3]*x[0]+par[4]*x[0]*x[0]+par[5]*x[0]*x[0]*x[0]);
    
    Int_t pdg = particle->GetPdgCode() ;
    Double_t x = particle->Pt() ;
    Int_t mother = particle->GetFirstMother() ;
    while(TMath::Abs(pdg)<100){//gamma, electrons, muons
        if(mother>-1){
            TParticle * tmpP=fStack->Particle(mother) ;
            pdg=tmpP->GetPdgCode() ;
            x = tmpP->Pt() ;
            mother = tmpP->GetFirstMother() ;
        }
        else{ //direct photon/electron....
            return 1.;
        }
    }
    

    // //lhc13e7    
    // Double_t w = 15.9695 -26.2752*x + 16.7259*x*x -4.77561*x*x*x +0.602569*x*x*x*x;
    
    // //single pi0: 0-25GeV-flat
    // //Double_t w = 1.48592 -2.78259*x + 2.01851*x*x -0.661467*x*x*x +0.0980217*x*x*x*x;

    // //single pi0: 0-6GeV-flat
    // //Double_t w = 55.6309 -102.175*x + 73.7241*x*x -24.3558*x*x*x +3.74631*x*x*x*x;

    Double_t w = 1.;
    if(fWeights) w = fWeights->Eval(x);

    return 1./w;
    
}



