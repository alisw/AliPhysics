// $Id$
//
// Standalone jet finder
//   CINT-compatible wrapper for AliFJWrapper, can be used in macros and from the ROOT command line.
//   Compiled code can use AliFJWrapper directly
//
// Authors: R.Haake

#include <vector>
#include <TLorentzVector.h>

#include "AliLog.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliJetContainer.h"
#include "AliFJWrapper.h"
#include "AliEmcalJetFinder.h"
#include "TH1.h"

//________________________________________________________________________
AliEmcalJetFinder::AliEmcalJetFinder() :
  TNamed("EmcalJetFinder","EmcalJetFinder"), fFastjetWrapper(0), fInputVectorIndex(0), fJetCount(0), fJetArray(), fGhostArea(0.005), fRadius(0.4), fJetAlgorithm(0), fRecombScheme(-1), fTrackMaxEta(0.9), fJetMaxEta(0.5), fJetMinPt(0), fJetMinArea(0)
{
  // Constructor
  fFastjetWrapper = new AliFJWrapper("FJWrapper", "FJWrapper");
}

//________________________________________________________________________
AliEmcalJetFinder::AliEmcalJetFinder(const char* name) :
  TNamed(name, name), fFastjetWrapper(0), fInputVectorIndex(0), fJetCount(0), fJetArray(), fGhostArea(0.005), fRadius(0.4), fJetAlgorithm(0), fRecombScheme(-1), fTrackMaxEta(0.9), fJetMaxEta(0.5), fJetMinPt(0), fJetMinArea(0)
{
  // Constructor
  fFastjetWrapper = new AliFJWrapper("FJWrapper", "FJWrapper");
}

//________________________________________________________________________
AliEmcalJetFinder::~AliEmcalJetFinder()
{
  if(fFastjetWrapper)
    delete fFastjetWrapper;
}


//________________________________________________________________________
Bool_t AliEmcalJetFinder::FindJets()
{
  // Tidy up and check input
  for (UInt_t i=0; i<fJetArray.size(); i++) {
    delete fJetArray[i];
    fJetArray[i] = 0;
  }
  fJetArray.clear();
  fJetCount = 0;
  if(!fInputVectorIndex)
  {
    AliError("No input vectors added to jet finder!");
    return kFALSE;
  }

  // Pass settings to fastjet
  fFastjetWrapper->SetAreaType(fastjet::active_area_explicit_ghosts);
  fFastjetWrapper->SetGhostArea(fGhostArea);
  fFastjetWrapper->SetR(fRadius);
  if(fJetAlgorithm == 0)
    fFastjetWrapper->SetAlgorithm(fastjet::antikt_algorithm);  
  if(fJetAlgorithm == 1)
    fFastjetWrapper->SetAlgorithm(fastjet::kt_algorithm);  
  if(fRecombScheme>=0)
    fFastjetWrapper->SetRecombScheme(static_cast<fastjet::RecombinationScheme>(fRecombScheme));

  fFastjetWrapper->SetMaxRap(fTrackMaxEta);

  // Run jet finding
  fFastjetWrapper->Run();

  // Save the found jets as light-weight objects
  std::vector<fastjet::PseudoJet> fastjets = fFastjetWrapper->GetInclusiveJets();
  fJetArray.resize(fastjets.size());

  for (UInt_t i=0; i<fastjets.size(); i++)
  {
    // Apply jet cuts
    if (fastjets[i].perp()<fJetMinPt) 
      continue;
    if (fFastjetWrapper->GetJetArea(i)<fJetMinArea)
      continue;
    if (TMath::Abs(fastjets[i].eta())>fJetMaxEta)
      continue;

    AliEmcalJet* jet = new AliEmcalJet(fastjets[i].perp(), fastjets[i].eta(), fastjets[i].phi(), fastjets[i].m());

    // Set the most important properties of the jet
    Int_t nConstituents(fFastjetWrapper->GetJetConstituents(i).size());
    jet->SetArea(fFastjetWrapper->GetJetArea(i));
    jet->SetNumberOfTracks(nConstituents);
    jet->SetNumberOfClusters(nConstituents);
    fJetArray[fJetCount] = jet;
    fJetCount++;
  }

  fJetArray.resize(fJetCount);

  //fastjets.clear(); // will be done by the destructor at the end of the function
  fFastjetWrapper->Clear();
  fInputVectorIndex = 0;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEmcalJetFinder::Filter(AliEmcalJet *pJet, AliJetContainer *pContJets, Double_t dVtx[3])
{
//
//  AliEmcalJetFinder::Filter
//

  fJetCount = 0;
  for (UInt_t j=0; j<fJetArray.size(); j++) {
    if (fJetArray[j]) { delete fJetArray[j]; fJetArray[j] = 0; }
  } fJetArray.clear();
//=============================================================================

  if ((!pJet) || (!pContJets)) return kFALSE;
//=============================================================================

  AliParticleContainer *pContTrks = pContJets->GetParticleContainer();
  if (pContTrks) for (Int_t i=0; i<pJet->GetNumberOfTracks(); i++) {
    AliVParticle *pTrk = pJet->TrackAt(i, pContTrks->GetArray()); if (!pTrk) continue;
    AddInputVector(pTrk->Px(), pTrk->Py(), pTrk->Pz(), pTrk->E(), pJet->TrackAt(i)+100);
  }

  AliClusterContainer *pContClus = pContJets->GetClusterContainer();
  if (pContClus) for (Int_t i=0; i<pJet->GetNumberOfClusters(); i++) {
    AliVCluster *pClu = pJet->ClusterAt(i, pContClus->GetArray()); if (!pClu) continue;

    TLorentzVector vClu; pClu->GetMomentum(vClu, dVtx);
    AddInputVector(vClu.Px(), vClu.Py(), vClu.Pz(), vClu.P(), -1*pJet->ClusterAt(i)-100);
  }
//=============================================================================

  if(!fInputVectorIndex) {
    AliError("No input vectors added to jet finder!");
    return kFALSE;
  }
//=============================================================================

  if (pJet->HasGhost()) {
    const std::vector<TLorentzVector> aGhosts = pJet->GetGhosts();
    for (UInt_t i=0; i<aGhosts.size(); i++) AddInputGhost(aGhosts[i].Px(),
                                                          aGhosts[i].Py(),
                                                          aGhosts[i].Pz(),
                                                          aGhosts[i].E());
  }
//=============================================================================

  fFastjetWrapper->SetR(fRadius);
  fFastjetWrapper->SetMaxRap(fTrackMaxEta);
  fFastjetWrapper->SetGhostArea(fGhostArea);
  if(fJetAlgorithm==0) fFastjetWrapper->SetAlgorithm(fastjet::antikt_algorithm);  
  if(fJetAlgorithm==1) fFastjetWrapper->SetAlgorithm(fastjet::kt_algorithm);  
  if(fRecombScheme>=0) fFastjetWrapper->SetRecombScheme(static_cast<fastjet::RecombinationScheme>(fRecombScheme));
//=============================================================================

  fFastjetWrapper->Filter();
  std::vector<fastjet::PseudoJet> aFilteredJets = fFastjetWrapper->GetFilteredJets();

  fJetArray.resize(aFilteredJets.size());
  for (UInt_t j=0; j<aFilteredJets.size(); j++) {
    if (aFilteredJets[j].perp()<fJetMinPt) continue;
    if (TMath::Abs(aFilteredJets[j].eta())>fJetMaxEta) continue;
    if (fFastjetWrapper->GetDoFilterArea()) if (fFastjetWrapper->GetFilteredJetArea(j)<fJetMinArea) continue;

    AliEmcalJet *piece = new AliEmcalJet(aFilteredJets[j].perp(),
                                         aFilteredJets[j].eta(),
                                         aFilteredJets[j].phi(),
                                         aFilteredJets[j].m());

    piece->SetLabel(j);
    if (fFastjetWrapper->GetDoFilterArea()) {
      fastjet::PseudoJet area(fFastjetWrapper->GetFilteredJetAreaVector(j));
      piece->SetArea(area.perp());
      piece->SetAreaEta(area.eta());
      piece->SetAreaPhi(area.phi());
      piece->SetAreaE(area.E());
    }

    std::vector<fastjet::PseudoJet> aConstis(fFastjetWrapper->GetFilteredJetConstituents(j));

    UInt_t nConstis = aConstis.size();
    piece->SetNumberOfTracks(nConstis);
    piece->SetNumberOfClusters(nConstis);

    Int_t nt = 0;
    Int_t nc = 0;
    for (UInt_t i=0; i<nConstis; i++) {
      Int_t uid = aConstis[i].user_index();
      if (uid>= 100) { piece->AddTrackAt(   1*uid-100, nt); ++nt; }
      if (uid<=-100) { piece->AddClusterAt(-1*uid-100, nc); ++nc; }
    }

    piece->SetNumberOfTracks(nt);
    piece->SetNumberOfClusters(nc);

    fJetArray[fJetCount] = piece;
    fJetCount++;
  }

  fJetArray.resize(fJetCount);
  fFastjetWrapper->Clear();
  fInputVectorIndex = 0;
  return kTRUE;
}

//________________________________________________________________________
void AliEmcalJetFinder::AddInputVector(Double_t px, Double_t py, Double_t pz)
{
  fFastjetWrapper->AddInputVector(px, py, pz, TMath::Sqrt(px*px + py*py + pz*pz), fInputVectorIndex + 100);
  fInputVectorIndex++;
}

//________________________________________________________________________
void AliEmcalJetFinder::AddInputVector(Double_t px, Double_t py, Double_t pz, Double_t E)
{
  fFastjetWrapper->AddInputVector(px, py, pz, E, fInputVectorIndex + 100);
  fInputVectorIndex++;
}

//________________________________________________________________________
void AliEmcalJetFinder::AddInputVector(Double_t px, Double_t py, Double_t pz, Int_t index)
{
  fFastjetWrapper->AddInputVector(px, py, pz, TMath::Sqrt(px*px + py*py + pz*pz), index);
  fInputVectorIndex++;
}

//________________________________________________________________________
void AliEmcalJetFinder::AddInputVector(Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index)
{
  fFastjetWrapper->AddInputVector(px, py, pz, E, index);
  fInputVectorIndex++;
}

//________________________________________________________________________
void AliEmcalJetFinder::AddInputGhost(Double_t px, Double_t py, Double_t pz, Double_t E)
{
  fFastjetWrapper->AddInputGhost(px, py, pz, E, -1);
}

//________________________________________________________________________
void AliEmcalJetFinder::FillPtHistogram(TH1* histogram)
{
  if(!histogram)
    return;
  for (std::size_t i=0; i<fJetArray.size(); i++)
  {
    histogram->Fill(fJetArray[i]->Pt());
  }
}

//________________________________________________________________________
void AliEmcalJetFinder::FillPhiHistogram(TH1* histogram)
{
  if(!histogram)
    return;
  for (std::size_t i=0; i<fJetArray.size(); i++)
  {
    histogram->Fill(fJetArray[i]->Phi());
  }
}

//________________________________________________________________________
void AliEmcalJetFinder::FillEtaHistogram(TH1* histogram)
{
  if(!histogram)
    return;
  for (std::size_t i=0; i<fJetArray.size(); i++)
  {
    histogram->Fill(fJetArray[i]->Eta());
  }
}



//________________________________________________________________________
Double_t AliEmcalJetFinder::Nsubjettiness(AliEmcalJet *pJet, AliJetContainer *pContJets,  Double_t dVtx[3], Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Int_t Option){
  
  fJetCount = 0;
  for (UInt_t j=0; j<fJetArray.size(); j++) {
    if (fJetArray[j]) { delete fJetArray[j]; fJetArray[j] = 0; }
  } fJetArray.clear();
//=============================================================================

  if ((!pJet) || (!pContJets)) return kFALSE;
//=============================================================================

  AliParticleContainer *pContTrks = pContJets->GetParticleContainer();
  if (pContTrks) for (Int_t i=0; i<pJet->GetNumberOfTracks(); i++) {
    AliVParticle *pTrk = pJet->TrackAt(i, pContTrks->GetArray()); if (!pTrk) continue;
    AddInputVector(pTrk->Px(), pTrk->Py(), pTrk->Pz(), pTrk->E(), pJet->TrackAt(i)+100);
  }
    AliClusterContainer *pContClus = pContJets->GetClusterContainer();
  if (pContClus) for (Int_t i=0; i<pJet->GetNumberOfClusters(); i++) {
      AliVCluster *pClu = pJet->ClusterAt(i, pContClus->GetArray()); if (!pClu) continue;

      TLorentzVector vClu; pClu->GetMomentum(vClu, dVtx);
      AddInputVector(vClu.Px(), vClu.Py(), vClu.Pz(), vClu.P(), -1*pJet->ClusterAt(i)-100);
      }
//=============================================================================

  if(!fInputVectorIndex) {
    AliError("No input vectors added to jet finder!");
  }

//=============================================================================

  if (pJet->HasGhost()) {
    const std::vector<TLorentzVector> aGhosts = pJet->GetGhosts();
    for (UInt_t i=0; i<aGhosts.size(); i++) AddInputGhost(aGhosts[i].Px(),
                                                          aGhosts[i].Py(),
                                                          aGhosts[i].Pz(),
                                                          aGhosts[i].E());
  }
  //these are all for the jet not subjets
  fFastjetWrapper->SetR(fRadius);
  fFastjetWrapper->SetMaxRap(fTrackMaxEta);
  fFastjetWrapper->SetGhostArea(fGhostArea);
  fFastjetWrapper->SetMinJetPt(fJetMinPt);
  if(fJetAlgorithm==0) fFastjetWrapper->SetAlgorithm(fastjet::antikt_algorithm);  //this is for the jet clustering not the subjet reclustering. 
  // if(fJetAlgorithm==1) fFastjetWrapper->SetAlgorithm(fastjet::kt_algorithm);  
  if(fRecombScheme>=0) fFastjetWrapper->SetRecombScheme(static_cast<fastjet::RecombinationScheme>(fRecombScheme));
  return fFastjetWrapper->AliFJWrapper::NSubjettiness(N,Algorithm,Radius, Beta, Option);

}





