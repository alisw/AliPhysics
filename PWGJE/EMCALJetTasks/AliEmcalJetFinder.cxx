// $Id$
//
// Standalone jet finder
//
// Authors: R.Haake

#include "AliFJWrapper.h"
#include "AliEmcalJet.h"
#include "AliEmcalJetFinder.h"
#include "AliLog.h"
#include <vector>
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
  if(fRecombScheme>0)
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
    jet->SetNumberOfTracks(nConstituents);
    jet->SetNumberOfClusters(nConstituents);
    fJetArray[fJetCount] = jet;
    fJetCount++;
  }

  fJetArray.resize(fJetCount);

  fastjets.clear();
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
void AliEmcalJetFinder::FillPtHistogram(TH1* histogram)
{
  if(!histogram)
    return;
  for (Int_t i=0; i<fJetArray.size(); i++)
  {
    histogram->Fill(fJetArray[i]->Pt());
  }
}

//________________________________________________________________________
void AliEmcalJetFinder::FillPhiHistogram(TH1* histogram)
{
  if(!histogram)
    return;
  for (Int_t i=0; i<fJetArray.size(); i++)
  {
    histogram->Fill(fJetArray[i]->Phi());
  }
}

//________________________________________________________________________
void AliEmcalJetFinder::FillEtaHistogram(TH1* histogram)
{
  if(!histogram)
    return;
  for (Int_t i=0; i<fJetArray.size(); i++)
  {
    histogram->Fill(fJetArray[i]->Eta());
  }
}
