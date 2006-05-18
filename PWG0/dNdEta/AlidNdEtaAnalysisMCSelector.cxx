#include "AlidNdEtaAnalysisMCSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVector3.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>
#include <AliTracker.h>

#include "../esdTrackCuts/AliESDtrackCuts.h"
#include "dNdEtaCorrection.h"
#include "dNdEtaAnalysis.h"

ClassImp(AlidNdEtaAnalysisMCSelector)

AlidNdEtaAnalysisMCSelector::AlidNdEtaAnalysisMCSelector(TTree * tree) :
  AlidNdEtaAnalysisSelector(tree)
{
  //
  // Constructor. Initialization of pointers
  //
}

AlidNdEtaAnalysisMCSelector::~AlidNdEtaAnalysisMCSelector()
{
  //
  // Destructor
  //
}

Bool_t AlidNdEtaAnalysisMCSelector::Process(Long64_t entry)
{
  //

  if (AliSelector::Process(entry) == kFALSE)
    return kFALSE;

  TTree* particleTree = GetKinematics();
  if (!fHeader || !particleTree)
    return kFALSE;

  // get the MC vertex
  AliGenEventHeader* genHeader = fHeader->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  TParticle* particle = 0;
  particleTree->SetBranchAddress("Particles", &particle);

  Int_t nPrim  = fHeader->GetNprimary();
  Int_t nTotal = fHeader->GetNtrack();

  for (Int_t i_mc = nTotal - nPrim; i_mc < nTotal; ++i_mc)
  {
    particleTree->GetEntry(i_mc);

    if (!particle)
      continue;

    if (IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    fdNdEtaAnalysis->FillTrack(vtxMC[2], particle->Eta(), 1);
  }
  fdNdEtaAnalysis->FillEvent(vtxMC[2]);

  return kTRUE;
}
