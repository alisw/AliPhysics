#include "AliForwardUtil.h"
#include <AliAnalysisManager.h>
#include "AliAODForwardMult.h"
#include <AliLog.h>
#include <AliInputEventHandler.h>
#include <AliESDEvent.h>
#include <AliPhysicsSelection.h>
#include <AliTriggerAnalysis.h>
#include <AliMultiplicity.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TMath.h>

//====================================================================
AliForwardUtil::Histos::~Histos()
{
  if (fFMD1i) delete fFMD1i;
  if (fFMD2i) delete fFMD2i;
  if (fFMD2o) delete fFMD2o;
  if (fFMD3i) delete fFMD3i;
  if (fFMD3o) delete fFMD3o;
}

//____________________________________________________________________
TH2D*
AliForwardUtil::Histos::Make(UShort_t d, Char_t r, 
			     const TAxis& etaAxis) const
{
  Int_t ns = (r == 'I' || r == 'i') ? 20 : 40;
  TH2D* hist = new TH2D(Form("FMD%d%c_cache", d, r), 
			Form("FMD%d%c cache", d, r),
			etaAxis.GetNbins(), etaAxis.GetXmin(), 
			etaAxis.GetXmax(), ns, 0, 2*TMath::Pi());
  hist->SetXTitle("#eta");
  hist->SetYTitle("#phi [radians]");
  hist->SetZTitle("d^{2}N_{ch}/d#etad#phi");
  hist->Sumw2();
  hist->SetDirectory(0);

  return hist;
}
//____________________________________________________________________
void
AliForwardUtil::Histos::Init(const TAxis& etaAxis)
{
  fFMD1i = Make(1, 'I', etaAxis);
  fFMD2i = Make(2, 'I', etaAxis);
  fFMD2o = Make(2, 'O', etaAxis);
  fFMD3i = Make(3, 'I', etaAxis);
  fFMD3o = Make(3, 'O', etaAxis);
}
//____________________________________________________________________
void
AliForwardUtil::Histos::Clear(Option_t* option)
{
  fFMD1i->Reset(option);
  fFMD2i->Reset(option);
  fFMD2o->Reset(option);
  fFMD3i->Reset(option);
  fFMD3o->Reset(option);
}

//____________________________________________________________________
TH2D*
AliForwardUtil::Histos::Get(UShort_t d, Char_t r) const
{
  switch (d) { 
  case 1: return fFMD1i;
  case 2: return (r == 'I' || r == 'i' ? fFMD2i : fFMD2o);
  case 3: return (r == 'I' || r == 'i' ? fFMD3i : fFMD3o);
  }
  return 0;
}
#if 0
//____________________________________________________________________
TH2D*
AliForwardUtil::Histos::Get(UShort_t d, Char_t r)
{
  switch (d) { 
  case 1: return fFMD1i;
  case 2: return (r == 'I' || r == 'i' ? fFMD2i : fFMD2o);
  case 3: return (r == 'I' || r == 'i' ? fFMD3i : fFMD3o);
  }
  default: return 0;
}
#endif
//====================================================================
TList*
AliForwardUtil::RingHistos::DefineOutputList(TList* d) const
{
  if (!d) return 0;
  TList* list = new TList;
  list->SetName(fName.Data());
  d->Add(list);
  return list;
}
//____________________________________________________________________
TList*
AliForwardUtil::RingHistos::GetOutputList(TList* d) const
{
  if (!d) return 0;
  TList* list = static_cast<TList*>(d->FindObject(fName.Data()));
  return list;
}

//____________________________________________________________________
TH1*
AliForwardUtil::RingHistos::GetOutputHist(TList* d, const char* name) const
{
  return static_cast<TH1*>(d->FindObject(name));
}

//====================================================================
Bool_t
AliForwardUtil::ReadTriggers(AliESDEvent* esd, UInt_t& triggers,
			     TH1I* hTriggers)
{
  triggers = 0;

  // Get the analysis manager - should always be there 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  if (!am) { 
    AliWarningGeneral("ReadTriggers","No analysis manager defined!");
    return kFALSE;
  }

  // Get the input handler - should always be there 
  AliInputEventHandler* ih = 
    static_cast<AliInputEventHandler*>(am->GetInputEventHandler());
  if (!ih) { 
    AliWarningGeneral("ReadTriggers","No input handler");
    return kFALSE;
  }
  
  // Get the physics selection - add that by using the macro 
  // AddTaskPhysicsSelection.C 
  AliPhysicsSelection* ps = 
    static_cast<AliPhysicsSelection*>(ih->GetEventSelection());
  if (!ps) { 
    AliWarningGeneral("ReadTriggers","No physics selection");
    return kFALSE;
  }
  
  // Check if this is a collision candidate (INEL)
  Bool_t inel = ps->IsCollisionCandidate(esd);
  if (inel) { 
    triggers |= AliAODForwardMult::kInel;
    hTriggers->Fill(0.5);
  }

  // IF this is inel, see if we have a tracklet 
  if (inel) { 
    const AliMultiplicity* spdmult = esd->GetMultiplicity();
    if (!spdmult) {
      AliWarningGeneral("ReadTriggers","No SPD multiplicity");
    }
    else { 
      Int_t n = spdmult->GetNumberOfTracklets();
      for (Int_t j = 0; j < n; j++) { 
	if(TMath::Abs(spdmult->GetEta(j)) < 1) { 
	  triggers |= AliAODForwardMult::kInelGt0;
	  hTriggers->Fill(1.5);
	  break;
	}
      }
    }
  }

  // Analyse some trigger stuff 
  AliTriggerAnalysis ta;
  if (ta.IsOfflineTriggerFired(esd, AliTriggerAnalysis::kNSD1)) {
    triggers |= AliAODForwardMult::kNSD;
    hTriggers->Fill(2.5);
  }

  // Get trigger stuff 
  TString trigStr = esd->GetFiredTriggerClasses();
  if (trigStr.Contains("CBEAMB-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kEmpty;
    hTriggers->Fill(3.5);
  }

  if (trigStr.Contains("CINT1A-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kA;
    hTriggers->Fill(4.5);
  }

  if (trigStr.Contains("CINT1B-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kB;
    hTriggers->Fill(5.5);
  }


  if (trigStr.Contains("CINT1C-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kC;
    hTriggers->Fill(6.5);
  }

  if (trigStr.Contains("CINT1-E-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kE;
    hTriggers->Fill(7.5);
  }

  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliForwardUtil::ReadVertex(AliESDEvent* esd, Double_t& vz, Double_t maxErr)
{
  // Get the vertex 
  const AliESDVertex* vertex = esd->GetPrimaryVertexSPD();
  if (!vertex) { 
#ifdef VERBOSE
    AliWarningGeneral("ReadVertex","No SPD vertex found in ESD");
#endif
    return kFALSE;
  }

  // Check that enough tracklets contributed 
  if(vertex->GetNContributors() <= 0) {
#ifdef VERBOSE
    AliWarningGeneral("ReadVertex",
		      Form("Number of contributors to vertex is %d<=0",
			   vertex->GetNContributors()));
#endif
    return kFALSE;
  }

  // Check that the uncertainty isn't too large 
  if (vertex->GetZRes() > maxErr) { 
#ifdef VERBOSE
    AliWarningGeneral("ReadVertex",
		      Form("Uncertaintity in Z of vertex is too large %f > %d", 
			   vertex->GetZRes(), maxErr));
#endif
    return kFALSE;
  }

  // Get the z coordiante 
  vz = vertex->GetZ();
	       
  return kTRUE;
}

//
// EOF
//
