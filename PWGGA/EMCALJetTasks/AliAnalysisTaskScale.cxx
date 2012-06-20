// $Id$
//
// Scale task.
//
// Author: R.Reed, M.Connors

#include "AliAnalysisTaskScale.h"

#include <TClonesArray.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVTrack.h"

ClassImp(AliAnalysisTaskScale)

//________________________________________________________________________
AliAnalysisTaskScale::AliAnalysisTaskScale() : 
  AliAnalysisTaskEmcal("AliAnalysisTaskScale", kTRUE), 
  fScaleFunction(0),
  fGeom(0),
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
  fHistClusterEtaPhi(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskScale::AliAnalysisTaskScale(const char *name) :
  AliAnalysisTaskEmcal(name, kTRUE), 
  fScaleFunction(0),
  fGeom(0),
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
  fHistClusterEtaPhi(0)
{
  // Constructor.

}

//________________________________________________________________________
void AliAnalysisTaskScale::UserCreateOutputObjects()
{
  // Create my user objects.

  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();

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

  fOutput->Add(fHistCentrality);
  fOutput->Add(fHistPtTPCvsCent);
  fOutput->Add(fHistPtEMCALvsCent);
  fOutput->Add(fHistEtvsCent);
  fOutput->Add(fHistScalevsCent);
  fOutput->Add(fHistDeltaScalevsCent);
  fOutput->Add(fHistPtTPCvsNtrack);
  fOutput->Add(fHistPtEMCALvsNtrack);
  fOutput->Add(fHistEtvsNtrack);
  fOutput->Add(fHistScalevsNtrack);
  fOutput->Add(fHistDeltaScalevsNtrack);
  fOutput->Add(fHistTrackPtvsCent);
  fOutput->Add(fHistClusterPtvsCent);
  fOutput->Add(fHistTrackEtaPhi);
  fOutput->Add(fHistClusterEtaPhi);

  PostData(1, fOutput);
}

//________________________________________________________________________
Double_t AliAnalysisTaskScale::GetScaleFactor(Double_t cent)
{
  // Get scale function.

  Double_t scale = -1;
  if (fScaleFunction)
    scale = fScaleFunction->Eval(cent);
  return scale;
}

//________________________________________________________________________
void AliAnalysisTaskScale::ExecOnce() 
{
  // Init the analysis.

  fGeom = AliEMCALGeometry::GetInstance();

  if (!fGeom) {
    AliFatal("Can not create geometry");
    return;
  }

  AliAnalysisTaskEmcal::ExecOnce();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskScale::FillHistograms() 
{
  // Execute on each event.

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

  const Int_t Ntracks = fTracks->GetEntries();
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
    AliVTrack *track = static_cast<AliVTrack*>(fTracks->At(iTracks));

    if (!track)
      continue;

    if (!AcceptTrack(track))
      continue;

    if (TMath::Abs(track->Eta()) > 0.7)   // only accept tracks in the EMCal eta range
      continue;

    fHistTrackPtvsCent->Fill(fCent,track->Pt());
    fHistTrackEtaPhi->Fill(track->Eta(),track->Phi());

    ptTPC += track->Pt();
    if ((track->Phi() > EmcalMaxPhi) || (track->Phi() < EmcalMinPhi))
      continue;

    ptEMCAL += track->Pt();
  }

  if (ptTPC == 0) 
    return kFALSE;
  
  Double_t Et = 0;
  const Int_t Nclus = fCaloClusters->GetEntries();
  for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
    AliVCluster *c = static_cast<AliVCluster*>(fCaloClusters->At(iClus));
    if (!c)
      continue;

    if (!AcceptCluster(c))
      continue;

    TLorentzVector nPart;
    c->GetMomentum(nPart, fVertex);

    fHistClusterPtvsCent->Fill(fCent, nPart.Pt());
    fHistClusterEtaPhi->Fill(nPart.Eta(), nPart.Phi());

    Et += nPart.Pt();
  }

  const Double_t scalecalc = ((Et + ptEMCAL) / EmcalArea) * (TpcArea / ptTPC);
  const Double_t scale     = GetScaleFactor(fCent);

  fHistCentrality->Fill(fCent);
  fHistPtTPCvsCent->Fill(fCent, ptTPC);
  fHistPtEMCALvsCent->Fill(fCent, ptEMCAL);
  fHistEtvsCent->Fill(fCent, Et);
  fHistScalevsCent->Fill(fCent, scalecalc);
  fHistDeltaScalevsCent->Fill(fCent, scalecalc - scale);
  fHistPtTPCvsNtrack->Fill(Ntracks, ptTPC);
  fHistPtEMCALvsNtrack->Fill(Ntracks, ptEMCAL);
  fHistEtvsNtrack->Fill(Ntracks, Et);
  fHistScalevsNtrack->Fill(Ntracks, scalecalc);
  fHistDeltaScalevsNtrack->Fill(Ntracks, scalecalc - scale);

  return kTRUE;
}      

//________________________________________________________________________
void AliAnalysisTaskScale::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
