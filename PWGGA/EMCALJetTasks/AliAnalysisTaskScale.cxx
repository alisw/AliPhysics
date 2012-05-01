// $Id$
//
// Scale task.
//
//

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TTree.h>
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliVCluster.h"
#include "AliAnalysisTaskScale.h"

ClassImp(AliAnalysisTaskScale)

//________________________________________________________________________
AliAnalysisTaskScale::AliAnalysisTaskScale(const char *name) 
  : AliAnalysisTaskSE(name), 
    fTracksName("tracks"),
    fClustersName("clusters"),
    fESD(0), 
    fOutputList(0), 
    fHistCentrality(0), 
    fHistPtTPCvsCent(0), 
    fHistPtEMCALvsCent(0), 
    fHistEtvsCent(0),  
    fHistScalevsCent(0),  
    fHistDeltaScalevsCent(0), 
    fHistPtTPCvsNtrack(0), 
    fHistPtEMCALvsNtrack(0), 
    fHistEtvsNtrack(0),  
    fHistScalevsNtrack(0),  
    fHistDeltaScalevsNtrack(0)
{
  // Constructor.

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskScale::UserCreateOutputObjects()
{
  // Create my user objects.

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  fHistCentrality = new TH1F("fHistCentrality","Centrality",101,-1,100);

  fHistPtTPCvsCent = new TH2F("fHistPtTPCvsCent","rho vs cent",101,-1,100,500,0,1000);
  fHistPtEMCALvsCent = new TH2F("fHistPtEMCALvsCent","rho vs cent",101,-1,100,500,0,1000);
  fHistEtvsCent = new TH2F("fHistEtvsCent","rho vs cent",101,-1,100,500,0,1000);
  fHistScalevsCent = new TH2F("fHistScalevsCent","rho vs cent",101,-1,100,400,0,4.0);
  fHistDeltaScalevsCent = new TH2F("fHistDeltaScalevsCent","rho vs cent",101,-1,100,400,-2.0,2.0);
  
  fHistPtTPCvsNtrack = new TH2F("fHistPtTPCvsNtrack","rho vs cent",500,0,2500,500,0,1000);
  fHistPtEMCALvsNtrack = new TH2F("fHistPtEMCALvsNtrack","rho vs cent",500,0,2500,500,0,1000);
  fHistEtvsNtrack = new TH2F("fHistEtvsNtrack","rho vs cent",500,0,2500,500,0,1000);
  fHistScalevsNtrack = new TH2F("fHistScalevsNtrack","rho vs cent",500,0,2500,400,0,4.0);
  fHistDeltaScalevsNtrack = new TH2F("fHistDeltaScalevsNtrack","rho vs cent",500,0,2500,400,-2.0,2.0);
    
  fOutputList->Add(fHistCentrality);
  fOutputList->Add(fHistPtTPCvsCent);
  fOutputList->Add(fHistPtEMCALvsCent);
  fOutputList->Add(fHistEtvsCent);
  fOutputList->Add(fHistScalevsCent);
  fOutputList->Add(fHistDeltaScalevsCent);
  fOutputList->Add(fHistPtTPCvsNtrack);
  fOutputList->Add(fHistPtEMCALvsNtrack);
  fOutputList->Add(fHistEtvsNtrack);
  fOutputList->Add(fHistScalevsNtrack);
  fOutputList->Add(fHistDeltaScalevsNtrack);

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskScale::UserExec(Option_t *) 
{
  // Execute on each event.

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  TClonesArray *tracks = 0;
  TClonesArray *clusters = 0;
  AliCentrality *centrality = 0;
  TList *l = InputEvent()->GetList();

  tracks = dynamic_cast<TClonesArray*>(l->FindObject(fTracksName));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
    return;
  }
    
  clusters = dynamic_cast<TClonesArray*>(l->FindObject(fClustersName));
  if (!clusters){
    AliError(Form("Pointer to clusters %s == 0", fClustersName.Data() ));
    return;
  }

  centrality = dynamic_cast<AliCentrality*>(l->FindObject("Centrality"));
  float fCent = centrality->GetCentralityPercentile("V0M");
  fHistCentrality->Fill(fCent);

  float scale = 0.000066*fCent*fCent-0.0015*fCent+1.5; //march 25th fit
  float ptTPC = 0;float ptEMCAL=0;
  const Int_t Ntracks = tracks->GetEntries();
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
    AliVTrack *track = static_cast<AliVTrack*>(tracks->At(iTracks));
    if (!track)
      continue;
    if (fabs(track->Eta())>0.7)
      continue;
    ptTPC+=track->Pt();
    if ((track->Phi()>3.14)||(track->Phi()<1.4))
      continue;
    ptEMCAL+=track->Pt();
  } //track loop 
  
  Double_t vertex[3] = {0, 0, 0};
  
  float Et = 0;
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  const Int_t Nclus = clusters->GetEntries();
  for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
    
    AliVCluster *c = dynamic_cast<AliVCluster*>(clusters->At(iClus));
    if (!c->IsEMCAL())
      continue;
    TLorentzVector nPart;
    c->GetMomentum(nPart, vertex);
    Et += nPart.P();
  }
  
  float scalecalc = ((Et+ptEMCAL)/2.44)*(8.8/ptTPC); //2.44 = EMCAL Acc, 8.8 = TPC Acc
  fHistPtTPCvsCent->Fill(fCent,ptTPC);
  fHistPtEMCALvsCent->Fill(fCent,ptEMCAL);
  fHistEtvsCent->Fill(fCent,Et);
  fHistScalevsCent->Fill(fCent,scalecalc);
  fHistDeltaScalevsCent->Fill(fCent,scalecalc-scale);

  fHistPtTPCvsNtrack->Fill(Ntracks,ptTPC);
  fHistPtEMCALvsNtrack->Fill(Ntracks,ptEMCAL);
  fHistEtvsNtrack->Fill(Ntracks,Et);
  fHistScalevsNtrack->Fill(Ntracks,scalecalc);
  fHistDeltaScalevsNtrack->Fill(Ntracks,scalecalc-scale);

  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskScale::Terminate(Option_t *) 
{
  // Nothing to be done for the moment.
}
