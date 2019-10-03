#include <array>

#include <THistManager.h>
#include <TLorentzVector.h>

#include "AliAnalysisTaskEmcalNeutralJets.h"
#include "AliAnalysisUtils.h"
#include "AliClusterContainer.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVVertex.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalNeutralJets)
/// \endcond

AliAnalysisTaskEmcalNeutralJets::AliAnalysisTaskEmcalNeutralJets() :
  AliAnalysisTaskEmcalJet(),
  fHistos(),
  fTriggerBits(0),
  fTriggerString(""),
  fNameR02jets(""),
  fNameR04jets(""),
  fNameClusters(""),
  fR02jets(nullptr),
  fR04jets(nullptr),
  fClusters(nullptr)
{

}

AliAnalysisTaskEmcalNeutralJets::AliAnalysisTaskEmcalNeutralJets(const char *name) :
  AliAnalysisTaskEmcalJet(name, true),
  fHistos(),
  fTriggerBits(0),
  fTriggerString(""),
  fNameR02jets(""),
  fNameR04jets(""),
  fNameClusters(""),
  fR02jets(nullptr),
  fR04jets(nullptr),
  fClusters(nullptr)
{

}

AliAnalysisTaskEmcalNeutralJets::~AliAnalysisTaskEmcalNeutralJets() {
}

void AliAnalysisTaskEmcalNeutralJets::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fAliAnalysisUtils = new AliAnalysisUtils;

  fHistos = new THistManager("CaloJetHistos");

  fHistos->CreateTH1("hEvents", "Number of events", 1, 0.5, 1.5);
  const std::array<double, 2> kJetRadii = {0.2, 0.4};
  for(auto r : kJetRadii){
    fHistos->CreateTH1(Form("hJetPtRawR%02d", static_cast<Int_t>(10.*r)), Form("Raw jet pt spectrum for R=%.1f jets", r), 20, 0., 200 );
    fHistos->CreateTH1(Form("hJetPtRecR%02d", static_cast<Int_t>(10.*r)), Form("Rec jet pt spectrum for R=%.1f jets", r), 20, 0., 200 );

    for(int i = 10; i < 200; i += 10)
      fHistos->CreateTH1(Form("hR%02dJet%d%dPtRel", static_cast<int>(10.*r), i, i + 10), Form("Pt rel of electron candidates for neutral jets with "), 100, 0., 1.);
  }

  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);
  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalNeutralJets::IsEventSelected(){
  // Data set requires a reduced event selection - no tracking info
  if(!(fInputHandler->IsEventSelected() & fTriggerBits)) return false;
  if(fTriggerString){
    if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerString)) return false;
  }
  if(fAliAnalysisUtils->IsPileUpEvent(fInputEvent)) return false;
  if(TMath::Abs(fVertexSPD[2]) > 10.) return false;
  return true;
}

bool AliAnalysisTaskEmcalNeutralJets::Run(){
  if(!fR02jets) fR02jets = this->GetJetContainer(fNameR02jets.Data());
  if(!fR04jets) fR04jets = this->GetJetContainer(fNameR04jets.Data());
  if(!fClusters) fClusters = this->GetClusterContainer(fNameClusters.Data());
  return true;
}

bool AliAnalysisTaskEmcalNeutralJets::FillHistograms(){
  // Loop over all jets in EMCAL surface
  for(auto j : fR02jets->accepted()){
    if(!IsTrackInEmcalAcceptance(j, 0.2)) continue;   // select only jets fully contained inside EMCAL acceptance
    if(TMath::Abs(j->Eta()) > 0.5) continue;
    if(j->Pt() < 20.) continue;
    fHistos->FillTH1("hJetPtRawR02", j->Pt());

    int jetptmin = static_cast<int>(j->Pt()/10.) * 10;
    TLorentzVector jetvec(j->Px(), j->Py(), j->Pz(), j->E());
    for(int i = 0; j->Nn(); i++){
      AliVCluster *clust = j->ClusterAt(i, fClusters->GetArray());
      if(clust->E() < 1.) continue;
      TLorentzVector clustvec;
      clust->GetMomentum(clustvec, fVertexSPD);
      Double_t ptrel = clustvec.Perp(jetvec.Vect());
      fHistos->FillTH1(Form("hR02Jet%d%dPtRel", jetptmin, jetptmin+10), ptrel);
    }
  }

  for(auto j : fR04jets->accepted()){
    if(!IsTrackInEmcalAcceptance(j, 0.4)) continue;   // select only jets fully contained inside EMCAL acceptance
    if(TMath::Abs(j->Eta())) continue;
    if(j->Pt() < 20) continue;
    fHistos->FillTH1("hJetPtRawR02", j->Pt());

    int jetptmin = static_cast<int>(j->Pt()/10.) * 10;
    TLorentzVector jetvec(j->Px(), j->Py(), j->Pz(), j->E());
    for(int i = 0; j->Nn(); i++){
      AliVCluster *clust = j->ClusterAt(i, fClusters->GetArray());
      if(clust->E() < 1.) continue;
      TLorentzVector clustvec;
      clust->GetMomentum(clustvec, fVertexSPD);
      Double_t ptrel = clustvec.Perp(jetvec.Vect());
      fHistos->FillTH1(Form("hR04Jet%d%dPtRel", jetptmin, jetptmin+10), ptrel);
    }
  }
}
