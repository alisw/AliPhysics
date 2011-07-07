#ifndef __CINT__
#include "AliForwardCorrectionManager.h"
#include "AliFMDCorrELossFit.h"
#include <TAxis.h>
#include <TFile.h>
#include <TList.h>
#include <TError.h>
#include <TH1.h>
#include <TF1.h>
#include <TClonesArray.h>
#else
class TAxis;
class AliFMDCirrELossFit;
class TH1;

#endif

/**
 * Class to make correction object and save to file 
 *
 * Run like					
 *
 * @verbatim 
 * gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");
 * gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/Compile.C");
 * Compile("$ALICE_ROOT/PWG2/FORWARD/analysis2/MakeELossFit.C"); 
 * MakeELossFit mef(sys, cms, field, mc, "AnalysisResults.root"); 
 * mef.Run();
 * @endverbatim 
 * where @e sys, the collision system, is one of 
 * - 1: pp
 * - 2: PbPb
 * 
 * @e cms is the center of mass energy per nuclean in GeV (e.g., 2760
 * for first PbPb data), @e field is (signed) L3 magnetic in kG, and 
 * @e mc is wether this correction applies to MC data or not. 
 * 
 * The class generates a file like 
 * @verbatim
 * elossfits<sys>_<cms>GeV_<field>kG_<realmc>.root
 * @endverbatim 
 * in the working directory. This file can be moved to 
 * @verbatim
 * $(ALICE_ROOT)/PWG2/FORWARD/corrections/ELossFits
 * @endverbatim 
 * and stored in SVN for later use. 
 *
 * @depcrecated 
 * The class AliFMDELossFitter automatically generates the
 * AliFMDCorrELossFit object.
 *
 * @ingroup pwg2_forward_analysis_scripts_tests
 */

class MakeELossFit 
{
protected: 
public:
  TList* fFitter;
  TAxis* fAxis;
  TClonesArray fFits;
  UShort_t     fSys;
  UShort_t     fCMS;
  Short_t      fField;
  Bool_t       fMC;

  //__________________________________________________________________
  MakeELossFit(UShort_t    sys, 
	       UShort_t    cms, 
	       Short_t     field, 
	       Bool_t      mc, 
	       const char* filename="forward_eloss.root") 
    : fFitter(0), 
      fAxis(0),
      fFits("AliFMDCorrELossFit::ELossFit"),
      fSys(sys), 
      fCMS(cms), 
      fField(field), 
      fMC(mc)
  {
    TFile* file = TFile::Open(filename, "READ");
    if (!file) { 
      Error("MakeELossFit", "Failed to open file %s", filename);
      return;
    }
    TList* forward = static_cast<TList*>(file->Get("Forward"));
    // static_cast<TList*>(file->Get("PWG2forwardDnDeta/Forward"));
    if (!forward) { 
      Error("MakeELossFit", "Couldn't get forward list from %s", filename);
      return;
    }

    fFitter = static_cast<TList*>(forward->FindObject("fmdEnergyFitter"));
    if (!fFitter) { 
      Error("MakeELossFit", "Couldn't get fitter folder");
      return;
    }

    fAxis = static_cast<TAxis*>(fFitter->FindObject("etaAxis"));
    if (!fAxis) { 
      Error("MakeELossFit", "Couldn't get eta axis");
      return;
    }
    file->Close();

#if 0
    AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
    mgr.Init(sys, cms, field, mc, 0, true);
#endif
  }
  //__________________________________________________________________
  TList* GetRing(UShort_t d, Char_t r) const
  {
    TList* rL = static_cast<TList*>(fFitter->FindObject(Form("FMD%d%c",d,r)));
    if (!rL) { 
      Warning("DrawFits", "List FMD%d%c not found", d, r);
      return 0;
    }
    return rL;
  }
  //__________________________________________________________________
  TList* GetEDists(UShort_t d, Char_t r) const
  {
    TList* rList = GetRing(d, r);
    if (!rList) 
      return 0;
    // rList->ls();

    TList* edists = static_cast<TList*>(rList->FindObject("EDists"));
    if (!edists) { 
      Error("PrintFits", "Couldn't get EDists list for FMD%d%c", d,r);
      return 0; 
    }
    return edists;
  }
  //__________________________________________________________________
  TH1* GetEDist(UShort_t d, Char_t r, UShort_t etabin)
  {
    TList* eList = GetEDists(d, r);
    if (!eList) { 
      Warning("GetEDist", "No list for FMD%d%c", d, r);
      return 0;
    }
    // eList->ls();

    TString cmp(Form("FMD%d%c_etabin%03d", d, r, etabin));

    TIter next(eList);
    TObject* o = 0;
    while ((o = next())) { 
      if (!cmp.CompareTo(o->GetName())) {
	return static_cast<TH1*>(o);
      }
    }
    return 0;
  }
#ifndef __CINT__
  //__________________________________________________________________
  AliFMDCorrELossFit::ELossFit* FindBestFit(TH1* dist) 
  {
    TList* funcs = dist->GetListOfFunctions();
    TIter  next(funcs);
    TF1*   func = 0;
    fFits.Clear();
    Int_t  i = 0;
    Info("FindBestFit", "%s", dist->GetName());
    while ((func = static_cast<TF1*>(next()))) { 
      AliFMDCorrELossFit::ELossFit* fit = 
	new(fFits[i++]) AliFMDCorrELossFit::ELossFit(0,*func);
      fit->CalculateQuality(10, .20, 1e-7);
    }
    fFits.Sort(false);
    // fFits.Print();
    return static_cast<AliFMDCorrELossFit::ELossFit*>(fFits.At(0));
  }
#endif
  //__________________________________________________________________
  Bool_t Run()
  {
    if (!fFitter || !fAxis) { 
      Error("Run", "Missing objects");
      return kFALSE;
    }
    AliFMDCorrELossFit* obj = new AliFMDCorrELossFit;
    obj->SetEtaAxis(*fAxis);

    for (UShort_t d = 1; d <= 3; d++) { 
      Info("Run", "detector is FMD%d", d);
      UShort_t nQ = (d == 1 ? 1 : 2);
      for (UShort_t q = 0; q < nQ; q++) { 
	Char_t r = (q == 0 ? 'I' : 'O');
	Int_t nBin = fAxis->GetNbins();
	Info("Run", " ring is FMD%d%c - %d bins", d, r, nBin);
	
	Bool_t oneSeen = kFALSE;
	for (UShort_t b = 1; b <= nBin; b++) { 
	  TH1* h = GetEDist(d, r, b);
	  if (oneSeen && !h) break;
	  if (!h)            continue;
	  if (!oneSeen)      oneSeen = true;

	  AliFMDCorrELossFit::ELossFit* best = FindBestFit(h);
	  best->Print("");
	  best->fDet  = d; 
	  best->fRing = r;
	  best->fBin  = b;
	  // Double_t eta = fAxis->GetBinCenter(b);
	  Info("Run", "    bin=%3d->%6.4f", b, eta);
	  obj->SetFit(d, r, b, new AliFMDCorrELossFit::ELossFit(*best));
	}
      }
    }
    AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
    TString fname(mgr.GetFileName(AliForwardCorrectionManager::kELossFits,
				  fSys, fCMS, fField, fMC));
    TFile* output = TFile::Open(fname.Data(), "RECREATE");
    if (!output) { 
      Warning("Run", "Failed to open output file %s", fname.Data());
      return kFALSE;
    }
    obj->Write(mgr.GetObjectName(AliForwardCorrectionManager::kELossFits));
    output->Write();
    output->Close();
    Info("Run", "File %s created.  It should be copied to %s and stored in SVN",
	 fname.Data(),mgr.GetFileDir(AliForwardCorrectionManager::kELossFits));
    
    return kTRUE;
  }
};
//
// EOF
//
