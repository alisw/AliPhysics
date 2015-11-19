//====================================================================
// 
// Base class for classes that calculate the multiplicity in the
// forward regions event-by-event
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
#include "AliForwardMultiplicityBase.h"
#include "AliForwardCorrectionManager.h"
#include "AliForwardUtil.h"
#include "AliFMDCorrELossFit.h"
#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliFMDEventInspector.h"
// #include "AliMultEventClassifier.h"
#include "AliFMDESDFixer.h"
#include "AliFMDSharingFilter.h"
#include "AliFMDDensityCalculator.h"
#include "AliFMDCorrector.h"
#include "AliFMDHistCollector.h"
#include "AliFMDEventPlaneFinder.h"
#include "AliESDEvent.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TProfile.h>
#include <THStack.h>
#include <iostream>
#include <iomanip>
#define AOD_SLOT 3
#ifdef POST_AOD
# define DEFINE(N) DefineOutput(N,AliAODForwardMult::Class())
# define POST(N)   PostData(N,fAODFMD)
#else
# define DEFINE(N) do { } while(false)
# define POST(N)   do { } while(false)
#endif

//====================================================================
AliForwardMultiplicityBase::AliForwardMultiplicityBase(const char* name) 
  : AliBaseESDTask(name, "AliForwardMultiplicityBase", 
		   &(AliForwardCorrectionManager::Instance())),
    fEnableLowFlux(false), 
    fStorePerRing(false),
    fHData(0),
    fHistos(),
    fAODFMD(false),
    fAODEP(false),
    // fAODRef(),
    fRingSums(),
    fDoTiming(false),
    fHTiming(0),
    fHStatus(0),
    fAddMask(AliAODForwardMult::kInel)
{
  DGUARD(fDebug, 3,"Named CTOR of AliForwardMultiplicityBase %s",name);

  DEFINE(AOD_SLOT);
}


//____________________________________________________________________
void
AliForwardMultiplicityBase::SetDebug(Int_t dbg)
{
  // 
  // Set debug level 
  // 
  // Parameters:
  //    dbg debug level
  //
  AliBaseESDTask::       SetDebug(dbg);
  GetSharingFilter()	.SetDebug(dbg);
  GetDensityCalculator().SetDebug(dbg);
  GetCorrections()	.SetDebug(dbg);
  GetHistCollector()	.SetDebug(dbg);
  GetEventPlaneFinder()	.SetDebug(dbg);
}

//____________________________________________________________________
Bool_t
AliForwardMultiplicityBase::Book()
{
  // 
  // Create output objects 
  // 
  //
  DGUARD(fDebug,1,"Create user ouput");
  UInt_t what = AliForwardCorrectionManager::kAll;
  if (!fEnableLowFlux)
    what ^= AliForwardCorrectionManager::kDoubleHit;
  if (!GetESDFixer().IsUseNoiseCorrection()) 
    what ^= AliForwardCorrectionManager::kNoiseGain;
  if (!GetCorrections().IsUseVertexBias())
    what ^= AliForwardCorrectionManager::kVertexBias;
  if (!GetCorrections().IsUseAcceptance())
    what ^= AliForwardCorrectionManager::kAcceptance;
  if (!GetCorrections().IsUseMergingEfficiency())
    what ^= AliForwardCorrectionManager::kMergingEfficiency;
  fNeededCorrections = what;

  // GetMultEventClassifier().CreateOutputObjects(fList);
  GetESDFixer()           .CreateOutputObjects(fList);
  GetSharingFilter()	  .CreateOutputObjects(fList);
  GetDensityCalculator()  .CreateOutputObjects(fList);
  GetCorrections()	  .CreateOutputObjects(fList);
  GetHistCollector()	  .CreateOutputObjects(fList);
  GetEventPlaneFinder()	  .CreateOutputObjects(fList);

  TAxis tmp(1, 0, 1);
  fHistos.Init(tmp);

  if (fDebug > 1) fDoTiming = true;
  if (fDoTiming) { 
    fHTiming = new TProfile("timing", "Timing of task", 
			    kTimingTotal, 0.5, kTimingTotal+.5);
    fHTiming->SetDirectory(0);
    fHTiming->SetFillColor(kRed+1);
    fHTiming->SetFillStyle(3001);
    fHTiming->SetMarkerStyle(20);
    fHTiming->SetMarkerColor(kBlack);
    fHTiming->SetLineColor(kBlack);
    fHTiming->SetXTitle("Part");
    fHTiming->SetYTitle("#LTt_{part}#GT [CPU]");
    fHTiming->SetStats(0);
    TAxis* xaxis = fHTiming->GetXaxis();
    xaxis->SetBinLabel(kTimingEventInspector,	
		       GetEventInspector()   .GetName());
    xaxis->SetBinLabel(kTimingSharingFilter,	
		       GetSharingFilter()    .GetName());
    xaxis->SetBinLabel(kTimingDensityCalculator,	
		       GetDensityCalculator().GetName());
    xaxis->SetBinLabel(kTimingCorrections,	
		       GetCorrections()      .GetName());
    xaxis->SetBinLabel(kTimingHistCollector,	
		       GetHistCollector()    .GetName());
    xaxis->SetBinLabel(kTimingEventPlaneFinder,	
		       GetEventPlaneFinder() .GetName());
    xaxis->SetBinLabel(kTimingTotal, "Total");
    fList->Add(fHTiming);
  }
  fHStatus = new TH1I("status", "Status of task",16, .5, 16.5);
  fHStatus->SetDirectory(0);
  fHStatus->SetFillColor(kCyan+2);
  fHStatus->SetFillStyle(3001);
  fHStatus->SetLineColor(kBlack);
  fHStatus->SetMarkerStyle(20);
  fHStatus->SetMarkerColor(kBlack);
  fHStatus->SetYTitle("Events");
  TAxis* a = fHStatus->GetXaxis();
  a->SetBinLabel(kStatusNoEvent,	"No event");
  a->SetBinLabel(kStatusNoTrigger,	"No triggers");
  a->SetBinLabel(kStatusNoSPD,	    	"No SPD (not used)");
  a->SetBinLabel(kStatusNoFMD,	    	"No FMD");
  a->SetBinLabel(kStatusNoVertex,	"No Vertex");
  a->SetBinLabel(kStatusPileup,	    	"Pile-up");
  a->SetBinLabel(kStatusSPDOutlier,	"SPD-outlier");
  a->SetBinLabel(kStatusIPzOutOfRange,  "IP_{z} out of range");
  a->SetBinLabel(kStatusFailSharing,	"Merger failed");
  a->SetBinLabel(kStatusFailDensity,	"N_{ch} estimator failed");
  a->SetBinLabel(kStatusFailEventPlane, "#Phi_{R} estimator failed");
  a->SetBinLabel(kStatusOutlier,	"Too many outliers");
  a->SetBinLabel(kStatusFailCorrector,  "Corrector failed");
  a->SetBinLabel(kStatusFailCollector,  "Histogram collector failed");
  a->SetBinLabel(kStatusNotAdded,	"Not added");
  a->SetBinLabel(kStatusAllThrough,     "All through");
  fList->Add(fHStatus);

  POST(AOD_SLOT);
  return true;
}
//____________________________________________________________________
void
AliForwardMultiplicityBase::CreateBranches(AliAODHandler* ah)
{
  TObject* obj   = &fAODFMD; ah->AddBranch("AliAODForwardMult", &obj);
  TObject* epobj = &fAODEP;  ah->AddBranch("AliAODForwardEP", &epobj);
  // TObject* rmobj = &fAODRef; ah->AddBranch("AliAODMultEventClass", &rmobj);

  if (!fStorePerRing) return;
  
  AliWarning("Per-ring histograms in AOD\n"
	     "*********************************************************\n"
	     "* For each event 5 additional 2D histogram are stored   *\n"
	     "* in separate branches of the AODs.  This will increase *\n"
	     "* the size of the AODs - proceed with caution           *\n"
	     "*********************************************************");
  TObject* hists[] = { fHistos.fFMD1i, 
		       fHistos.fFMD2i, fHistos.fFMD2o, 
		       fHistos.fFMD3i, fHistos.fFMD3o };
  for (Int_t i = 0; i < 5; i++) { 
    ah->AddBranch("TH2D", &(hists[i]));
  }
}

//____________________________________________________________________
Bool_t
AliForwardMultiplicityBase::PreData(const TAxis& vertex, const TAxis& eta)
{
  // 
  // Initialise the sub objects and stuff.  Called on first event 
  // 
  //
  DGUARD(fDebug,1,"Initialize sub-algorithms");

  // Force this here so we select the proper quality 
  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();

  UInt_t what = fNeededCorrections;
  // Check that we have the energy loss fits, needed by 
  //   AliFMDSharingFilter 
  //   AliFMDDensityCalculator 
  if (what & AliForwardCorrectionManager::kELossFits && !fcm.GetELossFit()) 
    AliFatal(Form("No energy loss fits"));
  
  // Check that we have the double hit correction - (optionally) used by 
  //  AliFMDDensityCalculator 
  if (what & AliForwardCorrectionManager::kDoubleHit && !fcm.GetDoubleHit()) 
    AliFatal("No double hit corrections"); 

  // Check that we have the secondary maps, needed by 
  //   AliFMDCorrector 
  //   AliFMDHistCollector
  if (what & AliForwardCorrectionManager::kSecondaryMap && 
      !fcm.GetSecondaryMap()) 
    AliFatal("No secondary corrections");

  // Check that we have the vertex bias correction, needed by 
  //   AliFMDCorrector 
  if (what & AliForwardCorrectionManager::kVertexBias && !fcm.GetVertexBias())
    AliFatal("No event vertex bias corrections");

  // Check that we have the merging efficiencies, optionally used by 
  //   AliFMDCorrector 
  if (what & AliForwardCorrectionManager::kMergingEfficiency && 
      !fcm.GetMergingEfficiency()) 
    AliFatal("No merging efficiencies");

  // Check that we have the acceptance correction, needed by 
  //   AliFMDCorrector 
  if (what & AliForwardCorrectionManager::kAcceptance && !fcm.GetAcceptance())
    AliFatal("No acceptance corrections");

  // const AliFMDCorrELossFit* fits = fcm.GetELossFit();
  // fits->CacheBins(GetDensityCalculator().GetMinQuality());

  InitMembers(eta,vertex);
  
  GetDensityCalculator().SetupForData(eta);
  GetSharingFilter()	.SetupForData(eta);
  GetCorrections()	.SetupForData(eta);
  GetHistCollector()	.SetupForData(vertex,eta);

  GetEventPlaneFinder() .SetRunNumber(GetEventInspector().GetRunNumber());
  GetEventPlaneFinder()	.SetupForData(eta);
  
  fAODFMD.SetBit(AliAODForwardMult::kSecondary, 
		 GetCorrections().IsUseSecondaryMap());
  fAODFMD.SetBit(AliAODForwardMult::kVertexBias, 
		 GetCorrections().IsUseVertexBias());
  fAODFMD.SetBit(AliAODForwardMult::kAcceptance, 
		 GetCorrections().IsUseAcceptance());
  fAODFMD.SetBit(AliAODForwardMult::kMergingEfficiency, 
		 GetCorrections().IsUseMergingEfficiency());
  fAODFMD.SetBit(AliAODForwardMult::kSum, 
		 GetHistCollector().GetMergeMethod() == 
		 AliFMDHistCollector::kSum);
  return true;
}

//____________________________________________________________________
void
AliForwardMultiplicityBase::InitMembers(const TAxis& eta, const TAxis& /*pv*/)
{
  fHistos.ReInit(eta);
  fAODFMD.Init(eta);
  fAODEP.Init(eta);
  fRingSums.Init(eta);

  fHData = static_cast<TH2D*>(fAODFMD.GetHistogram().Clone("d2Ndetadphi"));
  fHData->SetStats(0);
  fHData->SetDirectory(0);
  fList->Add(fHData);

  TList* rings = new TList;
  rings->SetName("ringSums");
  rings->SetOwner();
  fList->Add(rings);

  rings->Add(fRingSums.Get(1, 'I'));
  rings->Add(fRingSums.Get(2, 'I'));
  rings->Add(fRingSums.Get(2, 'O'));
  rings->Add(fRingSums.Get(3, 'I'));
  rings->Add(fRingSums.Get(3, 'O'));
  fRingSums.Get(1, 'I')->SetMarkerColor(AliForwardUtil::RingColor(1, 'I'));
  fRingSums.Get(2, 'I')->SetMarkerColor(AliForwardUtil::RingColor(2, 'I'));
  fRingSums.Get(2, 'O')->SetMarkerColor(AliForwardUtil::RingColor(2, 'O'));
  fRingSums.Get(3, 'I')->SetMarkerColor(AliForwardUtil::RingColor(3, 'I'));
  fRingSums.Get(3, 'O')->SetMarkerColor(AliForwardUtil::RingColor(3, 'O'));
}
//____________________________________________________________________
Bool_t
AliForwardMultiplicityBase::PostEvent()
{
  POST(AOD_SLOT);
  return true;
}
//____________________________________________________________________
Bool_t
AliForwardMultiplicityBase::Finalize()
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  DGUARD(fDebug,1,"Processing the merged results");

  TList* list   = fList;
  TList* output = fResults;
  Double_t nTr = 0, nTrVtx = 0, nAcc = 0;
  MakeSimpledNdeta(list, output, nTr, nTrVtx, nAcc);
  AliInfoF("\n"
	   "\t# events w/trigger:                 %f\n"
	   "\t# events w/trigger+vertex:          %f\n"
	   "\t# events accepted by cuts:          %f", 
	   nTr, nTrVtx, nAcc);

  EstimatedNdeta(list, output);

  GetSharingFilter()	.Terminate(list,output,Int_t(nTr));
  GetDensityCalculator().Terminate(list,output,Int_t(nTrVtx));
  GetCorrections()	.Terminate(list,output,Int_t(nTrVtx));

  TProfile* timing  = static_cast<TProfile*>(list->FindObject("timing"));
  Int_t     nTiming = (timing ? timing->GetBinContent(timing->GetNbinsX()) : 0);
  if (timing && nTiming > 0) { 
    TProfile* p = static_cast<TProfile*>(timing->Clone());
    p->SetDirectory(0);
    p->Scale(100. / nTiming);
    p->SetYTitle("#LTt_{part}#GT/#LTt_{total}#GT [%]");
    p->SetTitle("Relative timing of task");
    output->Add(p);
  }
  return true;
}

//____________________________________________________________________
void
AliForwardMultiplicityBase::EstimatedNdeta(const TList* input, 
					   TList*       output) const
{
  MakeRingdNdeta(input, "ringSums", output, "ringResults");
}

//____________________________________________________________________
Bool_t
AliForwardMultiplicityBase::MakeSimpledNdeta(const TList* input, 
					     TList*       output,
					     Double_t&    nTr, 
					     Double_t&    nTrVtx, 
					     Double_t&    nAcc)
{
  // Get our histograms from the container 
  TH1I* hEventsTr    = 0;
  TH1I* hEventsTrVtx = 0;
  TH1I* hEventsAcc   = 0;
  TH1I* hTriggers    = 0;
  TH1*  hStatus      = dynamic_cast<TH1*>(input->FindObject("status"));
  if (!GetEventInspector().FetchHistograms(input, 
					   hEventsTr, 
					   hEventsTrVtx, 
					   hEventsAcc,
					   hTriggers)) { 
    AliError(Form("Didn't get histograms from event selector "
		  "(hEventsTr=%p,hEventsTrVtx=%p,hEventsAcc=%p,hTriggers=%p)", 
		  hEventsTr, hEventsTrVtx, hEventsAcc, hTriggers));
    input->ls();
    return false;
  }
  nTr              = hEventsTr->Integral();
  nTrVtx           = hEventsTrVtx->Integral();
  nAcc             = hEventsAcc->Integral();
  Double_t vtxEff  = nTrVtx / nTr;
  Double_t vtxEff2 = 0;
  if (hStatus) {
    Double_t nTrg    = hStatus->Integral(3,15);
    Double_t nTrgVtx = hStatus->Integral(6,15);
    vtxEff2          = (nTrg > 0 ? nTrgVtx / nTrg : 0);
  }
    
  TH2D*   hData   = static_cast<TH2D*>(input->FindObject("d2Ndetadphi"));
  if (!hData) { 
    AliError(Form("Couldn't get our summed histogram from output "
		  "list %s (d2Ndetadphi=%p)", input->GetName(), hData));
    input->ls();
    return false;
  }

  Int_t nY      = hData->GetNbinsY();
  TH1D* dNdeta  = hData->ProjectionX("dNdeta",  1,     nY, "e");
  TH1D* dNdeta_ = hData->ProjectionX("dNdeta_", 1,     nY, "e");
  TH1D* norm    = hData->ProjectionX("norm",    0,     0,  "");
  TH1D* phi     = hData->ProjectionX("phi",     nY+1,  nY+1,  "");
  dNdeta->SetTitle("dN_{ch}/d#eta in the forward regions");
  dNdeta->SetYTitle("#frac{1}{N}#frac{dN_{ch}}{d#eta}");
  dNdeta->SetMarkerColor(kRed+1);
  dNdeta->SetMarkerStyle(20);
  dNdeta->SetDirectory(0);

  dNdeta_->SetTitle("dN_{ch}/d#eta in the forward regions");
  dNdeta_->SetYTitle("#frac{1}{N}#frac{dN_{ch}}{d#eta}");
  dNdeta_->SetMarkerColor(kMagenta+1);
  dNdeta_->SetMarkerStyle(21);
  dNdeta_->SetDirectory(0);

  norm->SetTitle("Normalization to  #eta coverage");
  norm->SetYTitle("#eta coverage");
  norm->SetLineColor(kBlue+1);
  norm->SetMarkerColor(kBlue+1);
  norm->SetMarkerStyle(21);
  norm->SetFillColor(kBlue+1);
  norm->SetFillStyle(3005);
  norm->SetDirectory(0);

  phi->SetTitle("Normalization to  #phi acceptance");
  phi->SetYTitle("#phi acceptance");
  phi->SetLineColor(kGreen+1);
  phi->SetMarkerColor(kGreen+1);
  phi->SetMarkerStyle(20);
  phi->SetFillColor(kGreen+1);
  phi->SetFillStyle(3004);
  // phi->Scale(1. / nAcc);
  phi->SetDirectory(0);

  // dNdeta->Divide(norm);
  dNdeta->Divide(phi);
  dNdeta->SetStats(0);
  dNdeta->Scale(vtxEff,	"width");

  dNdeta_->Divide(norm);
  dNdeta_->SetStats(0);
  dNdeta_->Scale(vtxEff, "width");

  AliInfoF("All\n"
	   "\tNormalization eta:  %d\n"
	   "\tNormalization phi:  %d\n"
	   "\tVertex efficiency:  %f (%f)",
	   Int_t(norm->GetMaximum()), Int_t(phi->GetMaximum()), 
	   vtxEff, vtxEff2);
  output->Add(dNdeta);
  output->Add(dNdeta_);
  output->Add(norm);
  output->Add(phi);

  return true;
}

					     
//____________________________________________________________________
void
AliForwardMultiplicityBase::MakeRingdNdeta(const TList* input, 
					   const char*  inName,
					   TList*       output,
					   const char*  outName,
					   Int_t        style) const
{
  // Make dN/deta for each ring found in the input list.  
  // 
  // A stack of all the dN/deta is also made for easy drawing. 
  // 
  // Note, that the distributions are normalised to the number of
  // observed events only - they should be corrected for 
  DGUARD(fDebug,3,"Make first-shot ring dN/deta");

  if (!input) return;
  TList* list = static_cast<TList*>(input->FindObject(inName));
  if (!list) { 
    AliWarning(Form("No list %s found in %s", inName, input->GetName()));
    return;
  }
  
  TList* out = new TList;
  out->SetName(outName);
  out->SetOwner();
  output->Add(out);

  THStack*     dndetaRings = new THStack("all", "dN/d#eta per ring");
  const char*  names[]     = { "FMD1I", "FMD2I", "FMD2O", "FMD3I", "FMD3O", 0 };
  const char** ptr         = names;
  
  while (*ptr) { 
    TList* thisList = new TList;
    thisList->SetOwner();
    thisList->SetName(*ptr);
    out->Add(thisList);

    TH2D* h = static_cast<TH2D*>(list->FindObject(Form("%s_cache", *ptr)));
    if (!h) { 
      AliWarning(Form("Didn't find %s_cache in %s", *ptr, list->GetName()));
      ptr++;
      continue;
    }
    TH2D* sumPhi = static_cast<TH2D*>(h->Clone("sum_phi"));
    sumPhi->SetDirectory(0);
    thisList->Add(sumPhi);

    TH2D* sumEta = static_cast<TH2D*>(h->Clone("sum_eta"));
    sumEta->SetDirectory(0);
    thisList->Add(sumEta);
    
    Int_t nY   = sumEta->GetNbinsY();
    TH1D* etaCov =static_cast<TH1D*>(h->ProjectionX("etaCov", 0,    0,    ""));
    TH1D* phiAcc =static_cast<TH1D*>(h->ProjectionX("phiAcc", nY+1, nY+1, ""));

    etaCov->SetTitle("Normalization to #eta coverage");
    etaCov->SetYTitle("#eta coverage");
    etaCov->SetMarkerColor(kBlue+1);
    etaCov->SetFillColor(kBlue+1);
    etaCov->SetFillStyle(3005);
    etaCov->SetDirectory(0);
    
    phiAcc->SetTitle("Normalization to #phi acceptance");
    phiAcc->SetYTitle("#phi acceptance");
    phiAcc->SetMarkerColor(kGreen+1);
    phiAcc->SetFillColor(kGreen+1);
    phiAcc->SetFillStyle(3004);
    // phiAcc->Scale(1. / nAcc);
    phiAcc->SetDirectory(0);

    // Double_t s = (etaCov->GetMaximum() > 0 ? 1. / etaCov->GetMaximum() : 1);
    for (Int_t i = 1; i <= sumEta->GetNbinsX(); i++) { 
      for (Int_t j = 1; j <= nY; j++) { 
	Double_t c = sumEta->GetBinContent(i, j);
	Double_t e = sumEta->GetBinError(i, j);
	Double_t a = etaCov->GetBinContent(i);
	Double_t p = phiAcc->GetBinContent(i);
	// Double_t t = p; // * a
	sumEta->SetBinContent(i, j, a <= 0 ? 0 : c / a);
	sumEta->SetBinError(  i, j, a <= 0 ? 0 : e / a);
	sumPhi->SetBinContent(i, j, p <= 0 ? 0 : c / p);
	sumPhi->SetBinError(  i, j, p <= 0 ? 0 : e / p);
      }
    }
    // etaCov->Scale(s);
    // phiAcc->Scale(s);

    TH1D* resPhi  =static_cast<TH1D*>(sumPhi->ProjectionX("dndeta_phi",
							  1,nY,"e"));
    resPhi->SetMarkerStyle(style);
    resPhi->SetDirectory(0);
    resPhi->Scale(1, "width");

    TH1D* resEta  =static_cast<TH1D*>(sumEta->ProjectionX("dndeta_eta",
							  1,nY,"e"));
    resEta->SetMarkerStyle(style);
    resEta->SetDirectory(0);
    resEta->Scale(1, "width");

    thisList->Add(resEta);
    thisList->Add(etaCov);
    thisList->Add(resPhi);
    thisList->Add(phiAcc);
    dndetaRings->Add(resPhi);

    DMSG(fDebug,1, "%s Normalization eta=%8d phi=%8d",
	 ptr, Int_t(etaCov->GetMaximum()), Int_t(phiAcc->GetMaximum()));

    ptr++;
  }
  out->Add(dndetaRings);
}
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)
//____________________________________________________________________
void
AliForwardMultiplicityBase::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  AliBaseESDTask::Print(option);
  gROOT->IncreaseDirLevel();
  PFB("Enable low flux code", fEnableLowFlux);
  PFB("Store per-ring hists", fStorePerRing);
  PFB("Make timing histogram", fDoTiming);
  PFV("Trigger mask for adding", AliAODForwardMult::GetTriggerString(fAddMask));
  // gROOT->IncreaseDirLevel();
  // GetMultEventClassifier().Print(option);
  GetESDFixer()           .Print(option);        
  GetSharingFilter()      .Print(option);
  GetDensityCalculator()  .Print(option);
  GetCorrections()        .Print(option);
  GetHistCollector()      .Print(option);
  GetEventPlaneFinder()   .Print(option);
  // gROOT->DecreaseDirLevel();
  gROOT->DecreaseDirLevel();
}


//
// EOF
//
