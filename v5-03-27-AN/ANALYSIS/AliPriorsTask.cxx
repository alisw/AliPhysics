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

//-----------------------------------------------------------------------
// Example of task (running locally) to calculate 
// the a priori concentration within the new class
// AliEventPoolLoop
//-----------------------------------------------------------------------


#include <TParticle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TH1I.h>
#include <TChain.h>
#include <TH2D.h>

#include "AliPriorsTask.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliLog.h"

ClassImp(AliPriorsTask)

//__________________________________________________________________________
AliPriorsTask::AliPriorsTask() :
  fHistEventsProcessed(0x0),
  fCalcPriors(0x0),
  fNiterations(0x0)
{
  //
  //Default ctor
  //
  for(Int_t i=0; i< AliPID::kSPECIES; i++){
  fPriors[i]=0;
  fRecId[i]=0;
  fMCId[i]=0;
  }
}
//___________________________________________________________________________
AliPriorsTask::AliPriorsTask(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fHistEventsProcessed(0x0),
  fCalcPriors(0x0),
  fNiterations(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  
  for(Int_t i=0; i< AliPID::kSPECIES; i++){
  fPriors[i]=0;
  fRecId[i]=0;
  fMCId[i]=0;
  }
  
  
  Info("AliPriorsTask","Calling Constructor");

  /*
    DefineInput(0) and DefineOutput(0)
    are taken care of by AliAnalysisTaskSE constructor
  */
  DefineOutput(1,TH1I::Class());
  DefineOutput(2,TH2D::Class());
}

//___________________________________________________________________________
AliPriorsTask& AliPriorsTask::operator=(const AliPriorsTask& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fHistEventsProcessed = c.fHistEventsProcessed;
    fCalcPriors = c.fCalcPriors;
    for(Int_t i=0; i< AliPID::kSPECIES; i++){
      fPriors[i]=c.fPriors[i];
      fRecId[i]=c.fRecId[i];
      fMCId[i]=c.fMCId[i];
    }
    fNiterations=c.fNiterations;
  }
  return *this;
}

//___________________________________________________________________________
AliPriorsTask::AliPriorsTask(const AliPriorsTask& c) :
  AliAnalysisTaskSE(c),
  fHistEventsProcessed(c.fHistEventsProcessed),
  fCalcPriors(c.fCalcPriors),
  fNiterations(c.fNiterations)
{
  //
  // Copy Constructor
  //
    for(Int_t i=0; i< AliPID::kSPECIES; i++){
      fPriors[i]=c.fPriors[i];
      fRecId[i]=c.fRecId[i];
      fMCId[i]=c.fMCId[i];
    }
}

//___________________________________________________________________________
AliPriorsTask::~AliPriorsTask() {
  //
  // Destructor
  //
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
  if (fCalcPriors) delete fCalcPriors ;
}

//_________________________________________________
void AliPriorsTask::UserExec(Option_t *)
{
  //
  // Main loop function
  //
  AliDebug(1,"UserExec") ;
  
  if(!fInputEvent) {
    Printf("No input event!");
    return;
  }

  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  if (!fESD) {
    Error("UserExec","NO ESD FOUND!");
    return;
  }
  
  if (!fMCEvent) {
    Error("UserExec","NO MC INFO FOUND");
    return;
  }
  
  //loop on tracks
  for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfTracks(); iTrack++) {

    AliESDtrack* track = fESD->GetTrack(iTrack);
    if (!track) continue;

    //some quality/pid requirements
    if (!track->IsOn(AliESDtrack::kITSrefit)) continue;
    if (!track->IsOn(AliESDtrack::kESDpid))   continue;
    if (!track->IsOn(AliESDtrack::kTPCpid))   continue;
    if (!track->IsOn(AliESDtrack::kTOFpid))   continue;
    
    if(track->GetLabel()<0) continue;

    Double_t esdpid[AliPID::kSPECIES]; track->GetESDpid(esdpid);

    AliPID pid(esdpid,kTRUE); 
    pid.SetPriors(fPriors);
    AliPID::EParticleType id = pid.GetMostProbable();
    
    if ((Int_t)id < (Int_t)AliPID::kSPECIES)   fRecId[id]++; 
    
    if(fNiterations == 0 ){
      Int_t pdgcode = TMath::Abs(((fMCEvent->GetTrack(track->GetLabel()))->Particle())->GetPdgCode());
      switch(pdgcode){
      case 11  :         fMCId[AliPID::kElectron]++; break;
      case 13  :         fMCId[AliPID::kMuon]++;     break;
      case 211 :         fMCId[AliPID::kPion]++;     break;
      case 321 :         fMCId[AliPID::kKaon]++;     break;
      case 2212:         fMCId[AliPID::kProton]++;   break;
      default  :                                     break;
      }
    }
  } 
  

  fHistEventsProcessed->Fill(0);

  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCalcPriors) ;
}

//___________________________________________________________________________
void AliPriorsTask::Terminate(Option_t*) 
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  AliAnalysisTaskSE::Terminate();
  
  Double_t MCsum =0; for(Int_t i=0; i< AliPID::kSPECIES; i++ ) MCsum+=fMCId[i];
  if(MCsum>0) Printf(" MC Priors :  %f  %f  %f  %f  %f",fMCId[0]/MCsum,fMCId[1]/MCsum,fMCId[2]/MCsum,fMCId[3]/MCsum,fMCId[4]/MCsum);
  Printf("    Priors :  %f  %f  %f  %f  %f \n",fPriors[0],fPriors[1],fPriors[2],fPriors[3],fPriors[4]);
  
  TH2D *histo= dynamic_cast<TH2D*> (GetOutputData(2));
  for(Int_t ip = 0; ip < AliPID::kSPECIES; ip++) histo->Fill(fNiterMax-1,fPriors[ip]);
  histo->SetMarkerStyle(8);
  histo->DrawCopy();
  
  Int_t ltype[AliPID::kSPECIES]  = {1,2,3,2,1}; 
  Int_t lcolor[AliPID::kSPECIES] = {1,2,8,38,46}; 
  
  TLine *l[AliPID::kSPECIES];
  for(Int_t isp =0; isp < AliPID::kSPECIES; isp++) {
    if(MCsum>0)l[isp]= new TLine(0.,fMCId[isp]/MCsum,fNiterMax,fMCId[isp]/MCsum);
    l[isp]->SetLineWidth(isp+1);
    l[isp]->SetLineColor(lcolor[isp]);
    l[isp]->SetLineStyle(ltype[isp]);  
    l[isp]->Draw("same");
  }
  
  const char *part[AliPID::kSPECIES]={"e","#mu","#pi","K","p"};
  TLegend *mcpriors = new TLegend(0.4,0.6,0.6,0.8);
  for(Int_t ip=0; ip < AliPID::kSPECIES; ip++) mcpriors->AddEntry(l[ip],part[ip],"l");
  mcpriors->Draw("same");
  
}
//___________________________________________________________________________
void AliPriorsTask::UserCreateOutputObjects() 
{
  //
  //HERE ONE CAN CREATE OUTPUT OBJECTS 
  //BEFORE THE EXECUTION OF THE TASK
  //
  OpenFile(1);
  fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;
  
  OpenFile(2);
  fCalcPriors = new TH2D("priors",Form("priors in %i iterations",fNiterMax),fNiterMax*10,-0.5,fNiterMax+0.5,100,0,1);  
  fCalcPriors->SetXTitle("Iteration step");
  fCalcPriors->SetYTitle("concentrations");
}
//________________________________________________________________________

Bool_t AliPriorsTask::NotifyBinChange()
{
  //
  //METHOD CALLED AT THE END OF THE EVENT LOOP. 
  //HERE CONCENTRATIONS ARE CALCULATED, STORED
  //AND USED IN THE NEXT EVENT LOOP
  //
  
  Double_t sum =0;   
  for(Int_t i=0; i< AliPID::kSPECIES; i++) {
    sum+=fRecId[i];
  }
  if(sum == 0) return kFALSE;
  
  Double_t conc[AliPID::kSPECIES];
  for(Int_t i=0; i< AliPID::kSPECIES; i++) {
    conc[i]=fRecId[i]/sum; 
    fCalcPriors->Fill(fNiterations,conc[i]);
    fRecId[i]=0;
  }
  
  SetPriors(conc);
  fNiterations++;
  
  return kTRUE;
}
