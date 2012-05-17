// $Id$
//
// Scale task.
//
// Author: R.Reed, M.Connors

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliVCluster.h"

#include "AliAnalysisTaskScale.h"

ClassImp(AliAnalysisTaskScale)

//________________________________________________________________________
AliAnalysisTaskScale::AliAnalysisTaskScale(const char *name) :
  AliAnalysisTaskSE(name), 
  fTracksName("Tracks"),
  fClustersName("CaloClusters"),
  fScaleFunction(0),
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

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskScale::UserCreateOutputObjects()
{
  // Create my user objects.

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  fHistCentrality         = new TH1F("fHistCentrality","Centrality",         101, -1, 100);
  fHistPtTPCvsCent        = new TH2F("fHistPtTPCvsCent","rho vs cent",       101, -1, 100, 500,   0, 1000);
  fHistPtEMCALvsCent      = new TH2F("fHistPtEMCALvsCent","rho vs cent",     101, -1, 100, 500,   0, 1000);
  fHistEtvsCent           = new TH2F("fHistEtvsCent","rho vs cent",          101, -1, 100, 500,   0, 1000);
  fHistScalevsCent        = new TH2F("fHistScalevsCent","rho vs cent",       101, -1, 100, 400,   0, 4);
  fHistDeltaScalevsCent   = new TH2F("fHistDeltaScalevsCent","rho vs cent",  101, -1, 100, 400,  -2, 2);
  fHistPtTPCvsNtrack      = new TH2F("fHistPtTPCvsNtrack","rho vs cent",     500,  0, 2500, 500,  0, 1000);
  fHistPtEMCALvsNtrack    = new TH2F("fHistPtEMCALvsNtrack","rho vs cent",   500,  0, 2500, 500,  0, 1000);
  fHistEtvsNtrack         = new TH2F("fHistEtvsNtrack","rho vs cent",        500,  0, 2500, 500,  0, 1000);
  fHistScalevsNtrack      = new TH2F("fHistScalevsNtrack","rho vs cent",     500,  0, 2500, 400,  0, 4);
  fHistDeltaScalevsNtrack = new TH2F("fHistDeltaScalevsNtrack","rho vs cent",500,  0, 2500, 400, -2, 2);

  fNewScaleFunction       = new TF1("sfunc","pol2", -1, 100);
    
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
  //fOutputList->Add(fNewScaleFunction);

  PostData(1, fOutputList);
}

//________________________________________________________________________
Double_t AliAnalysisTaskScale::GetScaleFactor(Double_t cent)
{
  Double_t scale = -1;
  if (fScaleFunction)
    scale = fScaleFunction->Eval(cent);
  return scale;
}

//________________________________________________________________________
void AliAnalysisTaskScale::UserExec(Option_t *) 
{
  // Execute on each event.

  // get input collections
  TClonesArray *tracks = 0;
  TClonesArray *clusters = 0;
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

  // get centrality
  Double_t cent = -1;
  AliCentrality *centrality = InputEvent()->GetCentrality();
  if (centrality)
    cent = centrality->GetCentralityPercentile("V0M");
  if (cent < 0) {
    AliError(Form("Centrality negative: %f", cent));
    return;
  }
  fHistCentrality->Fill(cent);

  // get vertex
  Double_t vertex[3] = {0, 0, 0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // hardcoded geometrical parameters
  const Double_t EmcalMinPhi = 80 * TMath::DegToRad();
  const Double_t EmcalMaxPhi = 180 * TMath::DegToRad();
  const Double_t TpcMinPhi   = 0;
  const Double_t TpcMaxPhi   = 2 * TMath::Pi();
  const Double_t MaxEta      = 0.7;
  const Double_t TpcArea     = (TpcMaxPhi - TpcMinPhi) * (MaxEta * 2);
  const Double_t EmcalArea   = (EmcalMaxPhi - EmcalMinPhi) * (MaxEta * 2);

  Double_t ptTPC   = 0;
  Double_t ptEMCAL = 0;

  const Int_t Ntracks = tracks->GetEntries();
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
    AliVTrack *track = static_cast<AliVTrack*>(tracks->At(iTracks));

    if (!track)
      continue;

    if (TMath::Abs(track->Eta()) > 0.7)   // only accept tracks in the EMCal eta range
      continue;

    ptTPC += track->Pt();
    if ((track->Phi() > EmcalMaxPhi) || (track->Phi() < EmcalMinPhi))
      continue;

    ptEMCAL += track->Pt();
  }
  
  Double_t Et = 0;
  const Int_t Nclus = clusters->GetEntries();
  for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
    AliVCluster *c = static_cast<AliVCluster*>(clusters->At(iClus));
    if (!c)
      continue;
    if (!c->IsEMCAL())
      continue;
    TLorentzVector nPart;
    c->GetMomentum(nPart, vertex);
    Et += nPart.Pt();
  }
  
  const Double_t scalecalc = ((Et + ptEMCAL) / EmcalArea) * (TpcArea / ptTPC);

  const Double_t scale     = GetScaleFactor(cent);

  fHistPtTPCvsCent->Fill(cent, ptTPC);
  fHistPtEMCALvsCent->Fill(cent, ptEMCAL);
  fHistEtvsCent->Fill(cent, Et);
  fHistScalevsCent->Fill(cent, scalecalc);
  fHistDeltaScalevsCent->Fill(cent, scalecalc - scale);
  fHistPtTPCvsNtrack->Fill(Ntracks, ptTPC);
  fHistPtEMCALvsNtrack->Fill(Ntracks, ptEMCAL);
  fHistEtvsNtrack->Fill(Ntracks, Et);
  fHistScalevsNtrack->Fill(Ntracks, scalecalc);
  fHistDeltaScalevsNtrack->Fill(Ntracks, scalecalc - scale);

  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskScale::Terminate(Option_t *) 
{
  /*
  fHistScalevsCent->Fit(fNewScaleFunction, "NO");
  PostData(1, fOutputList);
  */
}
