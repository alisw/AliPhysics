// 
// Calculate the corrections in the forward regions
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODForwardMult 
// 
// Histograms 
//   
// Corrections used 
// 
#include "AliForwardMCCorrectionsTask.h"
#include "AliForwardCorrectionManager.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliLog.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliESDEvent.h"
#include "AliAODHandler.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAODForwardMult.h"
#include "AliFMDStripIndex.h"
#include "AliFMDCorrSecondaryMap.h"
#include <TH1.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TList.h>
#include <TROOT.h>
#include <TVector3.h>
#include <iostream>

//====================================================================
AliForwardMCCorrectionsTask::AliForwardMCCorrectionsTask()
  : AliBaseMCCorrectionsTask(),
    fTrackDensity(),
    fESDFMD(),
    fSecCorr(0)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
}

//____________________________________________________________________
AliForwardMCCorrectionsTask::AliForwardMCCorrectionsTask(const char* name)
  : AliBaseMCCorrectionsTask(name, &(AliForwardCorrectionManager::Instance())),
    fTrackDensity("trackDensity"),
    fESDFMD(),
    fSecCorr(0)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
}


//____________________________________________________________________
AliBaseMCCorrectionsTask::VtxBin*
AliForwardMCCorrectionsTask::CreateVtxBin(Double_t low, Double_t high)
{
  return new AliForwardMCCorrectionsTask::VtxBin(low,high, fEtaAxis);
}

//____________________________________________________________________
Bool_t
AliForwardMCCorrectionsTask::PreEvent()
{
  // Clear our ESD object 
  fESDFMD.Clear();
  return true;
}

//____________________________________________________________________
Bool_t
AliForwardMCCorrectionsTask::ProcessESD(const AliESDEvent& esd, 
					const AliMCEvent& mc, 
					AliBaseMCCorrectionsTask::VtxBin& bin,
					const TVector3& ip)
{
  AliESDFMD* esdFMD = esd.GetFMDData();
  const Float_t maxMult = 100;
  fTrackDensity.Calculate(*esdFMD, mc, ip, fESDFMD, bin.fPrimary);
  bin.fCounts->Fill(0.5);

  AliForwardMCCorrectionsTask::VtxBin& vb = 
    static_cast<AliForwardMCCorrectionsTask::VtxBin&>(bin);

  // And then bin the data in our vtxbin 
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q=0; q<nr; q++) { 
      Char_t      r = (q == 0 ? 'I' : 'O');
      UShort_t    ns= (q == 0 ?  20 :  40);
      UShort_t    nt= (q == 0 ? 512 : 256);
      TH2D*       h = vb.fHists.Get(d,r);

      for (UShort_t s=0; s<ns; s++) { 
	for (UShort_t t=0; t<nt; t++) {
	  Float_t mult = fESDFMD.Multiplicity(d,r,s,t);
	  
	  if (mult == 0 || mult > maxMult) continue;

	  Float_t phi = fESDFMD.Phi(d,r,s,t) / 180 * TMath::Pi();
	  Float_t eta = fESDFMD.Eta(d,r,s,t);
	  h->Fill(eta,phi,mult);
	} // for t
      } // for s 
    } // for q 
  } // for d
  return true;
}
//____________________________________________________________________
void
AliForwardMCCorrectionsTask::CreateCorrections(TList* results)
{
  fSecCorr = new AliFMDCorrSecondaryMap;
  fSecCorr->SetVertexAxis(fVtxAxis);
  fSecCorr->SetEtaAxis(fEtaAxis);
  results->Add(fSecCorr);
}

//____________________________________________________________________
Bool_t 
AliForwardMCCorrectionsTask::FinalizeVtxBin(AliBaseMCCorrectionsTask::VtxBin* 
					    bin,  UShort_t iVz) 
{
  
  AliForwardMCCorrectionsTask::VtxBin* vb = 
    static_cast<AliForwardMCCorrectionsTask::VtxBin*>(bin);
  vb->Terminate(fList, fResults, iVz, fSecCorr);
  return true;
}


//____________________________________________________________________
void
AliForwardMCCorrectionsTask::Print(Option_t* option) const
{
  AliBaseMCCorrectionsTask::Print(option);
  gROOT->IncreaseDirLevel();
  fTrackDensity.Print(option);
  gROOT->DecreaseDirLevel();
}

//====================================================================
AliForwardMCCorrectionsTask::VtxBin::VtxBin()
  : AliBaseMCCorrectionsTask::VtxBin(),
    fHists()
{
}
//____________________________________________________________________
AliForwardMCCorrectionsTask::VtxBin::VtxBin(Double_t low, 
					    Double_t high, 
					    const TAxis& axis)
  :  AliBaseMCCorrectionsTask::VtxBin(low, high, axis, 40),
    fHists()
{
  fHists.Init(axis);
}


//____________________________________________________________________
TList*
AliForwardMCCorrectionsTask::VtxBin::CreateOutputObjects(TList* l)
{
  TList* d = AliBaseMCCorrectionsTask::VtxBin::CreateOutputObjects(l);

  d->Add(fHists.fFMD1i);
  d->Add(fHists.fFMD2i);
  d->Add(fHists.fFMD2o);
  d->Add(fHists.fFMD3i);
  d->Add(fHists.fFMD3o);

  return d;
}

//____________________________________________________________________
TH2D*
AliForwardMCCorrectionsTask::VtxBin::MakeBg(const TH2D* hits, 
					    const TH2D* primary) const
{
  TH2D* h = static_cast<TH2D*>(hits->Clone());
  h->SetDirectory(0);
  TString n(h->GetName());
  n.ReplaceAll("_cache", "");
  h->SetName(n);
  h->Divide(primary);
  
  return h;
}
  
//____________________________________________________________________
void
AliForwardMCCorrectionsTask::VtxBin::Terminate(const TList* input, 
					    TList* output, 
					    UShort_t iVz,
					    AliFMDCorrSecondaryMap* map)
{
  TList* out = new TList;
  out->SetName(GetName());
  out->SetOwner();
  output->Add(out);

  TList* l = static_cast<TList*>(input->FindObject(GetName()));
  if (!l) { 
    AliError(Form("List %s not found in %s", GetName(), input->GetName()));
    return;
  }

  TH2D*   fmd1i = static_cast<TH2D*>(l->FindObject("FMD1I_cache"));
  TH2D*   fmd2i = static_cast<TH2D*>(l->FindObject("FMD2I_cache"));
  TH2D*   fmd2o = static_cast<TH2D*>(l->FindObject("FMD2O_cache"));
  TH2D*   fmd3i = static_cast<TH2D*>(l->FindObject("FMD3I_cache"));
  TH2D*   fmd3o = static_cast<TH2D*>(l->FindObject("FMD3O_cache"));
  TH2D*   primO = static_cast<TH2D*>(l->FindObject("primary"));
  if (!fmd1i || !fmd2i || !fmd2o || !fmd3i || !fmd3o || !primO) {
    AliError(Form("Missing histogram(s): %p,%p,%p,%p,%p,%p",
		  fmd1i, fmd2i, fmd2o, fmd3i, fmd3o, primO));
    return;
  }

  // Half coverage in phi for inners
  TH2D*   primI = static_cast<TH2D*>(primO->Clone());
  primI->SetDirectory(0);
  primI->RebinY(2); 

  TH2D* bg1i = MakeBg(fmd1i, primI);
  TH2D* bg2i = MakeBg(fmd2i, primI);
  TH2D* bg2o = MakeBg(fmd2o, primO);
  TH2D* bg3i = MakeBg(fmd3i, primI);
  TH2D* bg3o = MakeBg(fmd3o, primO);
  map->SetCorrection(1, 'I', iVz, bg1i);
  map->SetCorrection(2, 'I', iVz, bg2i);
  map->SetCorrection(2, 'O', iVz, bg2o);
  map->SetCorrection(3, 'I', iVz, bg3i);
  map->SetCorrection(3, 'O', iVz, bg3o);
  out->Add(bg1i);
  out->Add(bg2i);
  out->Add(bg2o);
  out->Add(bg3i);
  out->Add(bg3o);
 
}

//
// EOF
//
