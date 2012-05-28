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
#include "AliCentrality.h"
#include "AliVEvent.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEMCALGeometry.h"
#include "AliLog.h"

#include "AliAnalysisTaskScale.h"

ClassImp(AliAnalysisTaskScale)

//________________________________________________________________________
AliAnalysisTaskScale::AliAnalysisTaskScale(const char *name) :
  AliAnalysisTaskSE(name), 
  fTracksName("Tracks"),
  fClustersName("CaloClusters"),
  fScaleFunction(0),
  fGeom(0),
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
  fHistDeltaScalevsNtrack(0),
  fHistTrackPtvsCent(0), 
  fHistClusterPtvsCent(0),
  fHistTrackEtaPhi(0), 
  fHistClusterEtaPhi(0),
  fMinTrackPt(0.15),
  fMinClusterPt(0.15)
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

  fHistCentrality         = new TH1F("Centrality","Centrality",              101, -1, 100);
  fHistPtTPCvsCent        = new TH2F("PtTPCvsCent","rho vs cent",            101, -1, 100, 500,   0, 1000);
  fHistPtEMCALvsCent      = new TH2F("PtEMCALvsCent","rho vs cent",          101, -1, 100, 500,   0, 1000);
  fHistEtvsCent           = new TH2F("EtvsCent","rho vs cent",               101, -1, 100, 500,   0, 1000);
  fHistScalevsCent        = new TH2F("ScalevsCent","rho vs cent",            101, -1, 100, 400,   0, 4);
  fHistDeltaScalevsCent   = new TH2F("DeltaScalevsCent","rho vs cent",       101, -1, 100, 400,  -2, 2);
  fHistPtTPCvsNtrack      = new TH2F("PtTPCvsNtrack","rho vs cent",          500,  0, 2500, 500,  0, 1000);
  fHistPtEMCALvsNtrack    = new TH2F("PtEMCALvsNtrack","rho vs cent",        500,  0, 2500, 500,  0, 1000);
  fHistEtvsNtrack         = new TH2F("EtvsNtrack","rho vs cent",             500,  0, 2500, 500,  0, 1000);
  fHistScalevsNtrack      = new TH2F("ScalevsNtrack","rho vs cent",          500,  0, 2500, 400,  0, 4);
  fHistDeltaScalevsNtrack = new TH2F("DeltaScalevsNtrack","rho vs cent",     500,  0, 2500, 400, -2, 2);
  fHistTrackPtvsCent      = new TH2F("TrackPtvsCent","Track pt vs cent",     101, -1, 100,  500,  0, 100);
  fHistClusterPtvsCent    = new TH2F("ClusterPtvsCent","Cluster pt vs cent", 101, -1, 100,  500,  0, 100);
  fHistTrackEtaPhi        = new TH2F("TrackEtaPhi","Track eta phi",          100, -1.0, 1.0, 64,  0, 6.4);
  fHistClusterEtaPhi      = new TH2F("ClusterEtaPhi","Cluster eta phi",      100, -1.0, 1.0, 64, -3.2, 3.2);

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
  fOutputList->Add(fHistTrackPtvsCent);
  fOutputList->Add(fHistClusterPtvsCent);
  fOutputList->Add(fHistTrackEtaPhi);
  fOutputList->Add(fHistClusterEtaPhi);

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

  if (!fGeom)
    fGeom = AliEMCALGeometry::GetInstance();

  if (!fGeom) {
    AliFatal("Can not create geometry");
    return;
  }

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
  
  const Double_t EmcalMinEta = fGeom->GetArm1EtaMin();
  const Double_t EmcalMaxEta = fGeom->GetArm1EtaMax();
  const Double_t EmcalMinPhi = fGeom->GetArm1PhiMin() * TMath::DegToRad();
  const Double_t EmcalMaxPhi = fGeom->GetArm1PhiMax() * TMath::DegToRad();

  const Double_t TpcMinPhi   = 0;
  const Double_t TpcMaxPhi   = 2 * TMath::Pi();

  const Double_t TpcArea     = (TpcMaxPhi - TpcMinPhi) * (EmcalMinEta - EmcalMaxEta);
  const Double_t EmcalArea   = (EmcalMaxPhi - EmcalMinPhi) * (EmcalMinEta - EmcalMaxEta);

  Double_t ptTPC   = 0;
  Double_t ptEMCAL = 0;

  const Int_t Ntracks = tracks->GetEntries();
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
    AliVTrack *track = static_cast<AliVTrack*>(tracks->At(iTracks));

    if (!track)
      continue;

    if (TMath::Abs(track->Eta()) > 0.7)   // only accept tracks in the EMCal eta range
      continue;

    fHistTrackPtvsCent->Fill(cent,track->Pt());
    fHistTrackEtaPhi->Fill(track->Eta(),track->Phi());

    if (track->Pt()< fMinTrackPt)
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

    fHistClusterPtvsCent->Fill(cent,nPart.Pt());
    fHistClusterEtaPhi->Fill(nPart.Eta(),nPart.Phi());

    if (nPart.Pt()< fMinClusterPt)
      continue;

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

}
