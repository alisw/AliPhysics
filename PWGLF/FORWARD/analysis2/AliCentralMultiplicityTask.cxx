//====================================================================
// 
// Base class for classes that calculate the multiplicity in the
// central region event-by-event
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODCentralMult 
// 
// Histograms 
//   
// Corrections used 
#include "AliCentralMultiplicityTask.h"
#include "AliCentralCorrectionManager.h"
#include "AliCentralCorrAcceptance.h"
#include "AliCentralCorrSecondaryMap.h"
#include "AliAODForwardMult.h"
#include "AliForwardUtil.h"
#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TError.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <iostream>
#include <iomanip>

//====================================================================
AliCentralMultiplicityTask::AliCentralMultiplicityTask(const char* name) 
  : AliAnalysisTaskSE(name),
    fInspector("centralEventInspector"),
    fList(0),
    fAODCentral(kFALSE),
    fUseSecondary(true),
    fUseAcceptance(true),
    fFirstEventSeen(false), 
  fIvz(0),
  fNClusterTracklet(0),
    fClusterPerTracklet(0),
    fNCluster(0),
  fNTracklet(0),
    fVtxList(0),
    fStore(false),
    fHData(0),
    fCorrManager(0)
{
  // 
  // Constructor 
  //   
  DGUARD(fDebug, 3,"Named CTOR of AliCentralMultiplicityTask: %s", name);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

  fCorrManager = &(AliCentralCorrectionManager::Instance());
  fBranchNames = 
    "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,"
    "SPDVertex.,PrimaryVertex.";
}
//____________________________________________________________________
AliCentralMultiplicityTask::AliCentralMultiplicityTask() 
  : AliAnalysisTaskSE(),
    fInspector(),
    fList(0),
    fAODCentral(),
    fUseSecondary(true),
    fUseAcceptance(true),
    fFirstEventSeen(false), 
  fIvz(0),
    fNClusterTracklet(0),
    fClusterPerTracklet(0),
    fNCluster(0),
    fNTracklet(0),
    fVtxList(0),
    fStore(false),
    fHData(0),
    fCorrManager(0)
{
  // 
  // Constructor 
  // 
  DGUARD(fDebug, 3,"Default CTOR of AliCentralMultiplicityTask");
}
//____________________________________________________________________
AliCentralMultiplicityTask::AliCentralMultiplicityTask(const AliCentralMultiplicityTask& o)
  : AliAnalysisTaskSE(o),
    fInspector(o.fInspector),
    fList(o.fList),
    fAODCentral(o.fAODCentral),
    fUseSecondary(o.fUseSecondary),
    fUseAcceptance(o.fUseAcceptance),
    fFirstEventSeen(o.fFirstEventSeen), 
  fIvz(o.fIvz),
    fNClusterTracklet(o.fNClusterTracklet),
    fClusterPerTracklet(o.fClusterPerTracklet),
    fNCluster(o.fNCluster),
    fNTracklet(o.fNTracklet),
    fVtxList(o.fVtxList),
    fStore(o.fStore),
    fHData(o.fHData),
    fCorrManager(o.fCorrManager)
{
  //
  // Copy constructor 
  // 
  DGUARD(fDebug, 3,"COPY CTOR of AliCentralMultiplicityTask");

}
//____________________________________________________________________
AliCentralMultiplicityTask&
AliCentralMultiplicityTask::operator=(const AliCentralMultiplicityTask& o)
{
  // 
  // Assignment operator 
  //
  DGUARD(fDebug,3,"Assignment of AliCentralMultiplicityTask");
  if (&o == this) return *this; 
  fInspector         = o.fInspector;
  fList              = o.fList;
  fAODCentral        = o.fAODCentral;
  fUseSecondary      = o.fUseSecondary;
  fUseAcceptance     = o.fUseAcceptance;
  fFirstEventSeen    = o.fFirstEventSeen;
  fIvz               = o.fIvz;
  fNClusterTracklet  = o.fNClusterTracklet;
  fClusterPerTracklet= o.fClusterPerTracklet;
  fNCluster          = o.fNCluster;
  fNTracklet         = o.fNTracklet;
  fVtxList           = o.fVtxList;
  fCorrManager       = o.fCorrManager;
  fStore             = o.fStore;
  fHData             = o.fHData;
  return *this;
}
//____________________________________________________________________
Bool_t 
AliCentralMultiplicityTask::Configure(const char* macro)
{
  // --- Configure the task ------------------------------------------
  TString macroPath(gROOT->GetMacroPath());
  if (!macroPath.Contains("$(ALICE_ROOT)/PWGLF/FORWARD/analysis2")) { 
    macroPath.Append(":$(ALICE_ROOT)/PWGLF/FORWARD/analysis2");
    gROOT->SetMacroPath(macroPath);
  }
  TString mac(macro);
  if (mac.EqualTo("-default-")) 
    mac = "$(ALICE_ROOT)/PWGLF/FORWARD/analysis2/CentralAODConfig.C";
  
  const char* config = gSystem->Which(gROOT->GetMacroPath(),mac.Data());
  if (!config) {
    AliWarningF("%s not found in %s", macro, gROOT->GetMacroPath());
    return false;
  }

  AliInfoF("Loading configuration of '%s' from %s", ClassName(), config);
  gROOT->Macro(Form("%s((AliCentralMultiplicityTask*)%p)", config, this));
  delete config;
  
  return true;
}

//____________________________________________________________________
void AliCentralMultiplicityTask::UserCreateOutputObjects() 
{
  // 
  // Create output objects 
  // 
  //
  DGUARD(fDebug,1,"Create user output in AliCentralMultiplicityTask");

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (ah) {
    //AliFatal("No AOD output handler set in analysis manager");
    TObject* obj = &fAODCentral;
    ah->AddBranch("AliAODCentralMult", &obj);
  } 
    
  fList = new TList();
  fList->SetOwner();

  fInspector.CreateOutputObjects(fList);

  PostData(1,fList);  
}

//____________________________________________________________________
AliESDEvent*
AliCentralMultiplicityTask::GetESDEvent()
{
  //
  // Get the ESD event. IF this is the first event, initialise
  //
  DGUARD(fDebug,1,"Get ESD event in AliCentralMultiplicityTask");
  if (IsZombie()) return 0;
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliWarning("No ESD event found for input event");
    return 0;
  }

  // Load in the data needed
  LoadBranches();

  // IF we've read the first event already, just return the event 
  if (fFirstEventSeen) return esd;
  
  // Read the details of the rung 
  fInspector.ReadRunDetails(esd);

  // If we weren't initialised before (i.e., in the setup), do so now. 
  AliCentralCorrectionManager& ccm = 
    AliCentralCorrectionManager::Instance();

  if (!ccm.Init(fInspector.GetRunNumber(),
		fInspector.GetCollisionSystem(),
		fInspector.GetEnergy(),
		fInspector.GetField())) {
    AliWarning("Failed  to intialize correction mananger");
  }
  //AliInfo("Manager of corrections in AliCentralMultiplicityTask init");
  Bool_t ok = true;
  if (/*fUseSecondary &&*/ !ccm.GetSecondaryMap()) {
    ok = false;
    AliError("No secondary correction defined!");
  }
  if (/*fUseAcceptance &&*/ !ccm.GetAcceptance()) {
    ok = false;
    AliError("No acceptance correction defined!");
  }
  // If the corrections are not seen, make this a zombie, and prevent
  // further execution of this task.
  if (!ok) { 
    AliError("Missing corrections, make this a zombie");
    SetZombie(true);
    esd = 0;
    fFirstEventSeen = true;
    return esd;
  }

  // Check for existence and get secondary map 
  const AliCentralCorrSecondaryMap* secMap = ccm.GetSecondaryMap(); 
  const TAxis& vaxis = secMap->GetVertexAxis();

  FindEtaLimits();

  fNClusterTracklet = new TH2D("nClusterVsnTracklet", 
			       "Total number of cluster vs number of tracklets",
			       100, 0, 10000, 100, 0, 10000);
  fNClusterTracklet->SetDirectory(0);
  fNClusterTracklet->SetXTitle("# of free clusters");
  fNClusterTracklet->SetYTitle("# of tracklets");
  fNClusterTracklet->SetStats(0);
  fList->Add(fNClusterTracklet);

  Int_t    nEta = 80;
  Double_t lEta = 2;
  fClusterPerTracklet = new TH2D("clusterPerTracklet", 
				 "N_{free cluster}/N_{tracklet} vs. #eta", 
				 nEta,-lEta,lEta, 101, -.05, 10.05);
  fClusterPerTracklet->SetDirectory(0);
  fClusterPerTracklet->SetXTitle("#eta");
  fClusterPerTracklet->SetYTitle("N_{free cluster}/N_{tracklet}");
  fClusterPerTracklet->SetStats(0);
  fList->Add(fClusterPerTracklet);

  // Cache histograms 
  fNCluster = new TH1D("cacheCluster", "", nEta,-lEta,lEta);
  fNCluster->SetDirectory(0);
  fNCluster->Sumw2();
		       
  fNTracklet = new TH1D("cacheTracklet", "", nEta,-lEta,lEta);
  fNTracklet->SetDirectory(0);
  fNTracklet->Sumw2();

  fList->Add(AliForwardUtil::MakeParameter("secondary",  fUseSecondary));
  fList->Add(AliForwardUtil::MakeParameter("acceptance", fUseAcceptance));

  fHData = static_cast<TH2D*>(fAODCentral.GetHistogram().Clone("d2Ndetadphi"));
  fHData->SetStats(0);
  fHData->SetDirectory(0);
  fList->Add(fHData);

  // Initialize the inspecto 
  fInspector.SetupForData(vaxis);
  fFirstEventSeen = kTRUE;

  fAODCentral.SetBit(AliAODCentralMult::kSecondary,  fUseSecondary);
  fAODCentral.SetBit(AliAODCentralMult::kAcceptance, fUseAcceptance);
  
  // Print some information 
  Print("R");

  return esd;
}
//____________________________________________________________________
void
AliCentralMultiplicityTask::MarkEventForStore() const
{
  // Make sure the AOD tree is filled 
  DGUARD(fDebug,1,"Mark AOD event for store in AliCentralMultiplicityTask");
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (ah) {  
    //AliFatal("No AOD output handler set in analysis manager");
    ah->SetFillAOD(kTRUE);
  }
}

//____________________________________________________________________
void AliCentralMultiplicityTask::FindEtaLimits()
{
  // Find our pseudo-rapidity limits 
  // 
  // Uses the secondary map to do so.
  DGUARD(fDebug,1,"Find eta limits in AliCentralMultiplicityTask");
  AliCentralCorrectionManager& ccm = 
    AliCentralCorrectionManager::Instance();
  const AliCentralCorrSecondaryMap* secMap = ccm.GetSecondaryMap();
  const TAxis&                      vaxis  = secMap->GetVertexAxis();

  unsigned short s = 1;
  TH2D* hCoverage = new TH2D("coverage", "#eta coverage per v_{z}", 
			     secMap->GetCorrection(s)->GetXaxis()->GetNbins(),
			     secMap->GetCorrection(s)->GetXaxis()->GetXmin(),
			     secMap->GetCorrection(s)->GetXaxis()->GetXmax(),
			     vaxis.GetNbins(),vaxis.GetXmin(),vaxis.GetXmax());
  hCoverage->SetDirectory(0);
  hCoverage->SetXTitle("#eta");
  hCoverage->SetYTitle("v_{z} [cm]");
  hCoverage->SetZTitle("n_{bins}");
  
  fAODCentral.Init(*(secMap->GetCorrection(s)->GetXaxis()));
  
  UShort_t nVz = vaxis.GetNbins();
  fVtxList     = new TObjArray(nVz, 1);
  fVtxList->SetName("centMultVtxBins");
  fVtxList->SetOwner();
  
  // Bool_t store = false;
  for (Int_t v = 1; v <= nVz; v++) { 
    VtxBin* bin = new VtxBin(v, vaxis.GetBinLowEdge(v), vaxis.GetBinUpEdge(v));
    bin->SetupForData(fList, hCoverage, fStore);
    fVtxList->AddAt(bin, v);
  }
  fList->Add(hCoverage);
}

//____________________________________________________________________
void AliCentralMultiplicityTask::UserExec(Option_t* /*option*/) 
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  
  DGUARD(fDebug,1,"Process event in AliCentralMultiplicityTask");
  fAODCentral.Clear("");

  AliESDEvent* esd = GetESDEvent();
  if (!esd) return;

  fIvz               = 0;
  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip;
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fInspector.Process(esd, triggers, lowFlux, 
					  ivz, ip, cent, nClusters);

  // No event or no trigger 
  if (found & AliFMDEventInspector::kNoEvent)    return;
  if (found & AliFMDEventInspector::kNoTriggers) return;
  
  // Make sure AOD is filled
  MarkEventForStore();

  if (found == AliFMDEventInspector::kNoSPD)      return;
  if (found == AliFMDEventInspector::kNoVertex)   return;
  if (triggers & AliAODForwardMult::kPileUp)      return;
  if (found == AliFMDEventInspector::kBadVertex)  return; // Out of range
  
  //Doing analysis
  const AliMultiplicity* spdmult = esd->GetMultiplicity();

  TH2D& aodHist = fAODCentral.GetHistogram();

  VtxBin* bin = static_cast<VtxBin*>(fVtxList->At(ivz));
  if (!bin) return;

  ProcessESD(aodHist, spdmult);
  bin->Correct(aodHist, fUseSecondary, fUseAcceptance);
  
  if (triggers & AliAODForwardMult::kInel) 
    fHData->Add(&(fAODCentral.GetHistogram()));

  PostData(1,fList);
}
//____________________________________________________________________
void 
AliCentralMultiplicityTask::ProcessESD(TH2D& aodHist, 
				       const AliMultiplicity* spdmult) const
{
  DGUARD(fDebug,1,"Process the ESD in AliCentralMultiplicityTask");
  fNTracklet->Reset();
  fNCluster->Reset();

  //Filling clusters in layer 1 used for tracklets...
  for(Int_t j = 0; j< spdmult->GetNumberOfTracklets();j++) {
    Double_t eta = spdmult->GetEta(j);
    fNTracklet->Fill(eta);
    aodHist.Fill(eta,spdmult->GetPhi(j));
  }

  //...and then the unused ones in layer 1 
  for(Int_t j = 0; j< spdmult->GetNumberOfSingleClusters();j++) {
    Double_t eta = -TMath::Log(TMath::Tan(spdmult->GetThetaSingle(j)/2.));
    fNCluster->Fill(eta);
    aodHist.Fill(eta, spdmult->GetPhiSingle(j));
  }
  fNClusterTracklet->Fill(fNCluster->GetEntries(), 
			  fNTracklet->GetEntries());
  
  fNCluster->Divide(fNTracklet);
  for (Int_t j = 1; j <= fNCluster->GetNbinsX(); j++)  
    fClusterPerTracklet->Fill(fNCluster->GetXaxis()->GetBinCenter(j), 
			      fNCluster->GetBinContent(j));

}

//____________________________________________________________________
void AliCentralMultiplicityTask::Terminate(Option_t* /*option*/) 
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  DGUARD(fDebug,1,"Process merged output in AliCentralMultiplicityTask");

  TList* list = dynamic_cast<TList*>(GetOutputData(1));
  if (!list) {
    AliError(Form("No output list defined (%p)", GetOutputData(1)));
    if (GetOutputData(1)) GetOutputData(1)->Print();
    return;
  }
  TList* output = new TList;
  output->SetName(Form("%sResults", GetName()));
  output->SetOwner();

  Double_t nTr = 0, nTrVtx = 0, nAcc = 0;
  MakeSimpledNdeta(list, output, nTr, nTrVtx, nAcc);
  AliInfoF("\n"
	   "\t# events w/trigger:                 %f\n"
	   "\t# events w/trigger+vertex:          %f\n"
	   "\t# events accepted y cuts:           %f", 
	   nTr, nTrVtx, nAcc);

  PostData(2, output);
}

//____________________________________________________________________
Bool_t
AliCentralMultiplicityTask::MakeSimpledNdeta(const TList* input, 
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
  nTr             = hEventsTr->Integral();
  nTrVtx          = hEventsTrVtx->Integral();
  nAcc            = hEventsAcc->Integral();
  Double_t vtxEff = nTrVtx / nTr;
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

  output->Add(dNdeta);
  output->Add(dNdeta_);
  output->Add(norm);
  output->Add(phi);

  return true;
}

//____________________________________________________________________
void
AliCentralMultiplicityTask::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //

  std::cout << ClassName() << ": " << GetName() << "\n" 
	    << std::boolalpha
	    << "  Use secondary correction:  " << fUseSecondary << '\n'
	    << "  Use acceptance correction: " << fUseAcceptance << '\n' 
	    << "  Off-line trigger mask:  0x" 
	    << std::hex     << std::setfill('0') 
	    << std::setw (8) << fOfflineTriggerMask 
	    << std::dec     << std::setfill (' ') 
	    << std::noboolalpha << std::endl;
  
  AliCentralCorrectionManager& ccm = 
    AliCentralCorrectionManager::Instance();
  if (ccm.IsInit()) {
    const AliCentralCorrSecondaryMap* secMap = ccm.GetSecondaryMap();
    if (secMap) {
      const TAxis& vaxis = secMap->GetVertexAxis();
      fVtxList->ls();
      std::cout << "  Eta ranges:\n"
		<< "     Vertex        | Eta bins\n"
		<< "   bin     range   | \n"
		<< "   ----------------+-----------" << std::endl;
      for (Int_t v = 1; v <= vaxis.GetNbins(); v++) { 
	VtxBin* bin = static_cast<VtxBin*>(fVtxList->At(v));
	if (!bin) continue;
	bin->Print();
      }
    }
  }

  gROOT->IncreaseDirLevel();
  ccm.Print(option);
  fInspector.Print(option);
  gROOT->DecreaseDirLevel();
  
}

//====================================================================
AliCentralMultiplicityTask::VtxBin::VtxBin(Int_t iVz, 
					   Double_t minIpZ, 
					   Double_t maxIpZ) 
  : fId(iVz), 
    fMinIpZ(minIpZ), 
    fMaxIpZ(maxIpZ),
    fEtaMin(999), 
    fEtaMax(0),
    fSec(0),
    fAcc(0),
    fHits(0)
{
}
//____________________________________________________________________
AliCentralMultiplicityTask::VtxBin::VtxBin(const VtxBin& o) 
  : TObject(o),
    fId(o.fId), 
    fMinIpZ(o.fMinIpZ), 
    fMaxIpZ(o.fMaxIpZ),
    fEtaMin(o.fEtaMin), 
    fEtaMax(o.fEtaMax),
    fSec(o.fSec),
    fAcc(o.fAcc),
    fHits(o.fHits)
{
}
//____________________________________________________________________
AliCentralMultiplicityTask::VtxBin&
AliCentralMultiplicityTask::VtxBin::operator=(const VtxBin& o) 
{
  if (&o == this) return *this;
  fId		= o.fId; 
  fMinIpZ	= o.fMinIpZ; 
  fMaxIpZ	= o.fMaxIpZ;
  fEtaMin	= o.fEtaMin; 
  fEtaMax	= o.fEtaMax;
  fSec		= o.fSec;
  fAcc		= o.fAcc;
  fHits		= o.fHits;

  return *this;
}

//____________________________________________________________________
const char*
AliCentralMultiplicityTask::VtxBin::GetName() const
{
  return Form("%c%03d_%c%03d", 
	      (fMinIpZ >= 0 ? 'p' : 'm'), Int_t(TMath::Abs(fMinIpZ)), 
	      (fMaxIpZ >= 0 ? 'p' : 'm'), Int_t(TMath::Abs(fMaxIpZ)));
}

//____________________________________________________________________
void
AliCentralMultiplicityTask::VtxBin::SetupForData(TList* l, 
						 TH2* coverage, 
						 Bool_t store)
{
  TList* out = 0;
  if (store) { 
    out = new TList;
    out->SetName(GetName());
    out->SetOwner();
    l->Add(out);
  }

  AliCentralCorrectionManager& ccm = 
    AliCentralCorrectionManager::Instance();

  // Clean-up 
  if (fSec) { 
    // delete fSec;
    fSec = 0;
  }
  if (fAcc) { 
    // delete fAcc;
    fAcc = 0;
  }
  // Get secondary correction and make a projection onto eta
  TH2* sec = ccm.GetSecondaryMap()->GetCorrection(UShort_t(fId));
  TH1* acc = ccm.GetAcceptance()->GetCorrection(UShort_t(fId));
  fSec = static_cast<TH2*>(sec->Clone());
  fAcc = static_cast<TH1*>(acc->Clone());
  fSec->SetDirectory(0);
  fAcc->SetDirectory(0);

  TH1D* proj = fSec->ProjectionX("secondary");
  proj->SetDirectory(0);
  proj->Scale(1. / fSec->GetNbinsY());

  // Find lower bound on eta 
  fEtaMin = proj->GetNbinsX();
  for (Int_t e = 1; e <= proj->GetNbinsX(); e++) { 
    Double_t c = proj->GetBinContent(e);
    if (c > .5 /*&& TMath::Abs(c - prev) < .1*c*/) {
      fEtaMin = e;
      break;
    }
  }
  // Find upper bound on eta 
  fEtaMax = 1;
  for (Int_t e = proj->GetNbinsX(); e >= 1; e--) { 
    Double_t c = proj->GetBinContent(e);
    if (c > .5 /*&& TMath::Abs(c - prev) < .1*c*/) {
      fEtaMax = e;
      break;
    }
  }
  // Fill our coverage histogram
  for (Int_t nn = fEtaMin; nn<=fEtaMax; nn++) { 
    coverage->SetBinContent(nn,fId,1);
  }
  
  if (!store) {
    // If we're not asked to store anything, clean-up, and get out
    delete proj;
    return;
  }

  // Modify the title of the projection 
  proj->SetTitle(Form("Projection of secondary correction "
		      "for %+5.1f<v_{z}<%+5.1f",fMinIpZ, fMaxIpZ));
  proj->SetYTitle("#LT 2^{nd} correction#GT");
  proj->SetMarkerStyle(20);
  proj->SetMarkerColor(kBlue+1);
  out->Add(proj);

  // Make some histograms to store diagnostics 
  TH2D* obg = static_cast<TH2D*>(fSec->Clone("secondaryMapFiducial"));
  obg->SetTitle(Form("%s - fiducial volume", obg->GetTitle()));
  obg->GetYaxis()->SetTitle("#varphi");
  obg->SetDirectory(0);
  out->Add(obg);
    
  TH1D* after = static_cast<TH1D*>(proj->Clone("secondaryFiducial"));
  after->SetDirectory(0);
  after->GetYaxis()->SetTitle("#LT 2^{nd} correction#GT");
  after->SetTitle(Form("%s - fiducial volume", after->GetTitle()));
  after->SetMarkerColor(kRed+1);
  out->Add(after);

  if (fHits) { 
    // delete fHits;
    fHits = 0;
  }
  fHits = static_cast<TH2D*>(fSec->Clone("hitMap"));
  fHits->SetDirectory(0);
  fHits->SetTitle(Form("d^{2}N/d#eta d#phi for %+5.1f<v_{z}<%+5.1f",
		      fMinIpZ, fMaxIpZ));
  fHits->GetYaxis()->SetTitle("#varphi");
  fHits->GetZaxis()->SetTitle("d^{2}N/d#eta d#varphi");
  fHits->SetMarkerColor(kBlack);
  fHits->SetMarkerStyle(1);
  out->Add(fHits);
    
  // Get the acceptance, and store that 
  TH1D* accClone   = static_cast<TH1D*>(fAcc->Clone("acceptance"));
  accClone->SetTitle(Form("Acceptance for %+5.1f<v_{z}<%+5.1f",
			  fMinIpZ, fMaxIpZ));
  accClone->SetDirectory(0);
  out->Add(accClone);
    
  // Now zero content outside our eta range 
  for (Int_t e = 1; e < fEtaMin; e++) { 
    after->SetBinContent(e, 0);
    after->SetBinError(e, 0);
    for(Int_t nn =1; nn <=obg->GetNbinsY();nn++) 
      obg->SetBinContent(e,nn,0);
  }

  for (Int_t e = fEtaMax+1; e <= proj->GetNbinsX(); e++) { 
    after->SetBinContent(e, 0);
    after->SetBinError(e, 0);
    for(Int_t nn =1; nn <=obg->GetNbinsY();nn++)
      obg->SetBinContent(e,nn,0);
  }
}
//____________________________________________________________________
void
AliCentralMultiplicityTask::VtxBin::Correct(TH2D&  aodHist,
					    Bool_t useSecondary,
					    Bool_t useAcceptance,
					    Bool_t  sum) const
{
  if (useSecondary && fSec) aodHist.Divide(fSec);

  Int_t nY = aodHist.GetNbinsY();
  for(Int_t ix = 1; ix <= aodHist.GetNbinsX(); ix++) {
    Bool_t fiducial = true;
    if (ix < fEtaMin || ix > fEtaMax) fiducial = false;
    //  Bool_t etabinSeen = kFALSE;  

    Double_t eta    = aodHist.GetXaxis()->GetBinCenter(ix);
    Int_t    iax    = fAcc->GetXaxis()->FindBin(eta);
    Float_t  accCor = fAcc->GetBinContent(iax);
    // For test
    // Float_t accErr = fAcc->GetBinError(ix);

    // Loop over phi 
    for(Int_t iy = 1; iy <= nY; iy++) {
      // If outside our fiducial volume, zero content 
      if (!fiducial) { 
	aodHist.SetBinContent(ix, iy, 0);
	aodHist.SetBinError(ix, iy, 0);
	continue;
      }
      // Get currrent value 
      Float_t aodValue = aodHist.GetBinContent(ix,iy);
      Float_t aodErr   = aodHist.GetBinError(ix,iy);

      // Ignore very small values
      if (aodValue < 0.000001) { 
	aodHist.SetBinContent(ix,iy, 0); 
	aodHist.SetBinError(ix,iy, 0); 
	continue; 
      }
      if (!useAcceptance) continue; 

      // Acceptance correction 
      Float_t accTmp = accCor;
      if (accTmp   < 0.000001) accTmp = 1;
      Float_t aodNew   = aodValue / accTmp ;
      aodHist.SetBinContent(ix,iy, aodNew);
      aodHist.SetBinError(ix,iy,aodErr);
      // - Test - 
      // Float_t error    = aodNew*TMath::Sqrt(TMath::Power(aodErr/aodValue,2) +
      // TMath::Power(accErr/accCor,2) );
      // test - aodHist.SetBinError(ix,iy,error);
    } // for (iy)
    //Filling underflow bin if we eta bin is in range
    if (fiducial) {
      aodHist.SetBinContent(ix,0, 1.);
      aodHist.SetBinContent(ix,nY+1, accCor);
    }
  } // for (ix)
  if (sum && fHits) fHits->Add(&aodHist);
}
    
//____________________________________________________________________
void
AliCentralMultiplicityTask::VtxBin::Print(Option_t* /*option*/) const
{
  std::cout << "   " 
	    << std::setw(2) << fId << "  " 
	    << std::setw(5) << fMinIpZ << "-"
	    << std::setw(5) << fMaxIpZ << " | "
	    << std::setw(3) << fEtaMin << "-" 
	    << std::setw(3) << fEtaMax << std::endl;
}

//
// EOF
//
