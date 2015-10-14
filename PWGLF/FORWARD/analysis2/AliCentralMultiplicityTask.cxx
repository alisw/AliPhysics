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
  : AliBaseESDTask(name, "AliCentralMultiplicityTask", 
		   &(AliCentralCorrectionManager::Instance())),
    fInspector("event"),
    fAODCentral(kFALSE),
    fUseSecondary(true),
    fUseAcceptance(true),
    fIvz(0),
    fNClusterTracklet(0),
    fClusterPerTracklet(0),
    fNCluster(0),
    fNTracklet(0),
    fVtxList(0),
    fStore(false),
    fHData(0)
{
  // 
  // Constructor 
  //   
  DGUARD(fDebug, 3,"Named CTOR of AliCentralMultiplicityTask: %s", name);
  fBranchNames = 
    "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,"
    "SPDVertex.,PrimaryVertex.";
}
//____________________________________________________________________
AliCentralMultiplicityTask::AliCentralMultiplicityTask() 
  : AliBaseESDTask(),
    fInspector(),
    fAODCentral(),
    fUseSecondary(true),
    fUseAcceptance(true),
  fIvz(0),
    fNClusterTracklet(0),
    fClusterPerTracklet(0),
    fNCluster(0),
    fNTracklet(0),
    fVtxList(0),
    fStore(false),
    fHData(0)
{
  // 
  // Constructor 
  // 
  DGUARD(fDebug, 3,"Default CTOR of AliCentralMultiplicityTask");
}

//____________________________________________________________________
void
AliCentralMultiplicityTask::CreateBranches(AliAODHandler* ah) 
{
  // 
  // Create output objects 
  // 
  //
  DGUARD(fDebug,1,"Create user output in AliCentralMultiplicityTask");

  if (!ah) 
    //AliFatal("No AOD output handler set in analysis manager");
    return;

  TObject* obj = &fAODCentral;
  ah->AddBranch("AliAODCentralMult", &obj);
}
//____________________________________________________________________
Bool_t
AliCentralMultiplicityTask::Book()
{
  fNeededCorrections = (AliCentralCorrectionManager::kSecondaryMap|
			AliCentralCorrectionManager::kAcceptance);
  return true;
}

//____________________________________________________________________
Bool_t
AliCentralMultiplicityTask::PreData(const TAxis& vertex, const TAxis& eta)
{
  // FindEtaLimits();
  // unsigned short s = 1;
  TH2D* hCoverage = new TH2D("coverage", "#eta coverage per v_{z}", 
			     eta.GetNbins(), eta.GetXmin(), eta.GetXmax(),
			     vertex.GetNbins(),vertex.GetXmin(),
			     vertex.GetXmax());
  hCoverage->SetDirectory(0);
  hCoverage->SetXTitle("#eta");
  hCoverage->SetYTitle("v_{z} [cm]");
  hCoverage->SetZTitle("n_{bins}");
  
  fAODCentral.Init(eta);
  
  UShort_t nVz = vertex.GetNbins();
  fVtxList     = new TObjArray(nVz, 1);
  fVtxList->SetName("centMultVtxBins");
  fVtxList->SetOwner();
  
  // Bool_t store = false;
  for (Int_t v = 1; v <= nVz; v++) { 
    VtxBin* bin = new VtxBin(v, vertex.GetBinLowEdge(v), 
			     vertex.GetBinUpEdge(v));
    bin->SetupForData(fList, hCoverage, fStore);
    fVtxList->AddAt(bin, v);
  }
  fList->Add(hCoverage);

  // Bins are 
  TArrayD bins;
  // N-bins, loweset order, higest order, return array
  AliForwardUtil::MakeLogScale(300, 0, 5, bins);
  fNClusterTracklet = new TH2D("nClusterVsnTracklet", 
			       "Total number of cluster vs number of tracklets",
			       bins.GetSize()-1, bins.GetArray(),
			       bins.GetSize()-1, bins.GetArray());
  fNClusterTracklet->SetDirectory(0);
  fNClusterTracklet->SetXTitle("N_{free cluster}");
  fNClusterTracklet->SetYTitle("N_{tracklet}");
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
  // fInspector.SetupForData(vertex);

  fAODCentral.SetBit(AliAODCentralMult::kSecondary,  fUseSecondary);
  fAODCentral.SetBit(AliAODCentralMult::kAcceptance, fUseAcceptance);

  return true;
}
//____________________________________________________________________
Bool_t 
AliCentralMultiplicityTask::PreEvent()
{
  fAODCentral.Clear("");
  return true;
}
//____________________________________________________________________
Bool_t 
AliCentralMultiplicityTask::Event(AliESDEvent& esd)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  
  DGUARD(fDebug,1,"Process event in AliCentralMultiplicityTask");
  fIvz               = 0;
  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip;
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fInspector.Process(&esd, triggers, lowFlux, 
					  ivz, ip, cent, nClusters);

  // No event or no trigger 
  if (found & AliFMDEventInspector::kNoEvent)    return false;
  if (found & AliFMDEventInspector::kNoTriggers) return false;
  
  // Make sure AOD is filled
  MarkEventForStore();

  if (found    & AliFMDEventInspector::kNoSPD)      return false;
  if (found    & AliFMDEventInspector::kNoVertex)   return false;
  if (triggers & AliAODForwardMult::kPileUp)        return false;
  if (found    & AliFMDEventInspector::kBadVertex)  return false; 
  
  //Doing analysis
  const AliMultiplicity* spdmult = esd.GetMultiplicity();

  TH2D& aodHist = fAODCentral.GetHistogram();

  VtxBin* bin = static_cast<VtxBin*>(fVtxList->At(ivz));
  if (!bin) return false;

  ProcessESD(aodHist, spdmult);
  bin->Correct(aodHist, fUseSecondary, fUseAcceptance);
  
  if (triggers & AliAODForwardMult::kInel) 
    fHData->Add(&(fAODCentral.GetHistogram()));

  return true;
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
Bool_t 
AliCentralMultiplicityTask::Finalize()
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  DGUARD(fDebug,1,"Process merged output in AliCentralMultiplicityTask");

  Double_t nTr = 0, nTrVtx = 0, nAcc = 0;
  MakeSimpledNdeta(fList, fResults, nTr, nTrVtx, nAcc);
  AliInfoF("\n"
	   "\t# events w/trigger:                 %f\n"
	   "\t# events w/trigger+vertex:          %f\n"
	   "\t# events accepted by cuts:          %f", 
	   nTr, nTrVtx, nAcc);
  return true;
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

#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)
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
AliCentralMultiplicityTask::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  AliBaseESDTask::Print(option);
  gROOT->IncreaseDirLevel();
  PFB("Use secondary correction", fUseSecondary);
  PFB("Use acceptance correction", fUseAcceptance);
  
  AliCentralCorrectionManager& ccm = 
    AliCentralCorrectionManager::Instance();
  if (ccm.IsInit()) {
    const AliCentralCorrSecondaryMap* secMap = ccm.GetSecondaryMap();
    if (secMap) {
      const TAxis& vaxis = secMap->GetVertexAxis();
      // fVtxList->ls();
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
