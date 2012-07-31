// This task finds the eventplane
// using the FMD
// 
#include "AliFMDEventPlaneFinder.h"
#include <TList.h>
#include <TMath.h>
#include "AliLog.h"
#include <TH2D.h>
#include <TH3D.h>
#include "TROOT.h"
#include <iostream>
#include <iomanip>
#include <TString.h>
#include <TFile.h>
#include "AliOADBContainer.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliEventplane.h"
#include "AliCentrality.h"
#include "AliAODForwardEP.h"
#include "AliAODForwardMult.h"
#include "AliAODEvent.h"

ClassImp(AliFMDEventPlaneFinder)
#if 0
; // For Emacs
#endif 

//____________________________________________________________________
AliFMDEventPlaneFinder::AliFMDEventPlaneFinder()
  : TNamed(),
    fList(0),
    fEvent(0),
    fQt(),
    fQa(),
    fQc(),
    fQ1(),
    fQ2(),
    fQeta(),
    fHepFMD(0),
    fHepFMDA(0),
    fHepFMDC(0),
    fHepFMDQC1(0),
    fHepFMDQC2(0),
    fHdiffFMDAC(0),
    fHdiffFMDTPC(0),
    fHdiffFMDVZERO(0),
    fHcorrFMDAC(0),
    fHcorrFMDTPC(0),
    fHcorrFMDVZERO(0),
    fHPhi(0),
    fDebug(0),
    fOADBFileName(0),
    fOADBContainer(0),
    fPhiDist(0),
    fRunNumber(0),
    fUsePhiWeights(0)
{
  // 
  // Constructor 
  //
  DGUARD(fDebug,0,"Default CTOR of AliFMDEventPlaneFinder");
   
}

//____________________________________________________________________
AliFMDEventPlaneFinder::AliFMDEventPlaneFinder(const char* title)
  : TNamed("fmdEventPlaneFinder", title), 
    fList(0),
    fEvent(0),
    fQt(),
    fQa(),
    fQc(),
    fQ1(),
    fQ2(),
    fQeta(),
    fHepFMD(0),
    fHepFMDA(0),
    fHepFMDC(0),
    fHepFMDQC1(0),
    fHepFMDQC2(0),
    fHdiffFMDAC(0),
    fHdiffFMDTPC(0),
    fHdiffFMDVZERO(0),
    fHcorrFMDAC(0),
    fHcorrFMDTPC(0),
    fHcorrFMDVZERO(0),
    fHPhi(0),
    fDebug(0),
    fOADBFileName(0),
    fOADBContainer(0),
    fPhiDist(0),
    fRunNumber(0),
    fUsePhiWeights(1)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of object
  //
  DGUARD(fDebug,0,"Named CTOR of AliFMDEventPlaneFinder: %s", title);
}

//____________________________________________________________________
AliFMDEventPlaneFinder::AliFMDEventPlaneFinder(const 
						 AliFMDEventPlaneFinder& o)
  : TNamed(o),
    fList(o.fList),
    fEvent(o.fEvent),
    fQt(o.fQt),
    fQa(o.fQa),
    fQc(o.fQc),
    fQ1(o.fQ1),
    fQ2(o.fQ2),
    fQeta(o.fQeta),
    fHepFMD(o.fHepFMD),
    fHepFMDA(o.fHepFMDA),
    fHepFMDC(o.fHepFMDC),
    fHepFMDQC1(o.fHepFMDQC1),
    fHepFMDQC2(o.fHepFMDQC2),
    fHdiffFMDAC(o.fHdiffFMDAC),
    fHdiffFMDTPC(o.fHdiffFMDTPC),
    fHdiffFMDVZERO(o.fHdiffFMDVZERO),
    fHcorrFMDAC(o.fHcorrFMDAC),
    fHcorrFMDTPC(o.fHcorrFMDTPC),
    fHcorrFMDVZERO(o.fHcorrFMDVZERO),
    fHPhi(o.fHPhi),
    fDebug(o.fDebug),
    fOADBFileName(o.fOADBFileName),
    fOADBContainer(o.fOADBContainer),
    fPhiDist(o.fPhiDist),
    fRunNumber(o.fRunNumber),
    fUsePhiWeights(o.fUsePhiWeights)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  DGUARD(fDebug,0,"Copy CTOR of AliFMDEventPlaneFinder");
}

//____________________________________________________________________
AliFMDEventPlaneFinder::~AliFMDEventPlaneFinder()
{
  // 
  // Destructor 
  //
}

//____________________________________________________________________
AliFMDEventPlaneFinder&
AliFMDEventPlaneFinder::operator=(const AliFMDEventPlaneFinder& o)
{
  // 
  // Assignement operator
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object
  //
  DGUARD(fDebug,3,"Assignment of AliFMDEventPlaneFinder");
  if (&o == this) return *this; 
  TNamed::operator=(o);

  fList                = o.fList;
  fEvent               = o.fEvent;
  fQt                  = o.fQt;
  fQa                  = o.fQa;
  fQc                  = o.fQc;
  fQ1                  = o.fQ1;
  fQ2                  = o.fQ2;
  fQeta                = o.fQeta;
  fHepFMD              = o.fHepFMD;
  fHepFMDA             = o.fHepFMDA;
  fHepFMDC             = o.fHepFMDC;
  fHepFMDQC1           = o.fHepFMDQC1;
  fHepFMDQC2           = o.fHepFMDQC2;
  fHdiffFMDAC          = o.fHdiffFMDAC;
  fHdiffFMDTPC         = o.fHdiffFMDTPC;
  fHdiffFMDVZERO       = o.fHdiffFMDVZERO;
  fHcorrFMDAC          = o.fHcorrFMDAC; 
  fHcorrFMDTPC         = o.fHcorrFMDTPC;
  fHcorrFMDVZERO       = o.fHcorrFMDVZERO;
  fHPhi                = o.fHPhi;
  fDebug               = o.fDebug;
  fOADBFileName        = o.fOADBFileName;
  fOADBContainer       = o.fOADBContainer;
  fPhiDist             = o.fPhiDist;
  fRunNumber           = o.fRunNumber;
  fUsePhiWeights       = o.fUsePhiWeights;

  return *this;
}

//____________________________________________________________________
TH1D*
AliFMDEventPlaneFinder::MakePsiRHist(const char* name, 
				     const char* title,
				     Int_t color)
{
  // Generate a 1D histogram of Psi_R. 
  TH1D* ret = new TH1D(name, Form("#Psi_{R} - %s",title), 100, 0, TMath::Pi());
  ret->SetDirectory(0);
  ret->SetXTitle("#Psi_{R} [radians]");
  ret->SetYTitle("Events");
  ret->SetLineColor(color);
  ret->SetFillColor(color);
  ret->SetFillStyle(3001);
  ret->Sumw2();
  fList->Add(ret);
  return ret;
}
//____________________________________________________________________
TH1D*
AliFMDEventPlaneFinder::MakeDiffHist(const char* name, 
				     const char* first,
				     const char* second,
				     Int_t color)
{
  TH1D* ret = new TH1D(name, Form("#Delta#Psi_{R} - %s minus %s",first,second), 
		       100, -TMath::Pi()/2, +TMath::Pi()/2);
  ret->SetDirectory(0);
  ret->SetXTitle(Form("#Psi_{R,%s}-#Psi_{R,%s} [radians]", first, second));
  ret->SetYTitle("Events");
  ret->SetLineColor(color);
  ret->SetFillColor(color);
  ret->SetFillStyle(3001);
  ret->Sumw2();
  fList->Add(ret);
  return ret;

}
//____________________________________________________________________
TH2F*
AliFMDEventPlaneFinder::MakeCorrHist(const char* name, 
				     const char* first,
				     const char* second)
{
  TH2F* ret = new TH2F(name, Form("#Psi_{R} - %s vs %s", first, second), 
		       100, 0, TMath::Pi(), 100, 0, TMath::Pi());
  ret->SetDirectory(0);
  ret->Sumw2();
  ret->SetXTitle(Form("#Psi_{R,%s} [radians]", first));
  ret->SetYTitle(Form("#Psi_{R,%s} [radians]", second));
  ret->SetZTitle("Events");
  fList->Add(ret);
  return ret;
}
   
//____________________________________________________________________
void
AliFMDEventPlaneFinder::Init(const TAxis& etaAxis)
{
  // Intialize this sub-algorithm 
  //
  // Parameters:
  //   etaAxis   fmd eta axis binning
  //
  DGUARD(fDebug,1,"Initalization of AliFMDEventPlaneFinder");

  fHepFMD    = MakePsiRHist("epFMD", "FMD", kRed+1);
  fHepFMDA   = MakePsiRHist("epFMDA","FMD - A side", kGreen+1);
  fHepFMDC   = MakePsiRHist("epFMDC","FMD - C side", kBlue+1);
  fHepFMDQC1 = MakePsiRHist("epFMDQC1", "FMD QC{1}", kMagenta+1);
  fHepFMDQC2 = MakePsiRHist("epEMDQC2", "FMD QC{2}", kCyan+1);

  fHdiffFMDAC     = MakeDiffHist("diffFMDAC",    "FMD-A", "FMD-C", kRed+1);
  fHdiffFMDTPC    = MakeDiffHist("diffFMDTPC",   "FMD",   "TPC",   kGreen+1);
  fHdiffFMDVZERO  = MakeDiffHist("diffFMDVZERO", "FMD",   "VZERO", kBlue+1);
  fHcorrFMDAC     = MakeCorrHist("corrFMDAC",    "FMD-A", "FMD-C");
  fHcorrFMDTPC    = MakeCorrHist("corrFMDTPC",   "FMD",   "TPC");
  fHcorrFMDVZERO  = MakeCorrHist("corrFMDVZERO", "FMD",   "VZERO");

  //
  fHPhi = new TH2D("hPhi", "Phi distribution in FMD",
		   etaAxis.GetNbins(), etaAxis.GetXmin(), etaAxis.GetXmax(), 
		   20, 0., TMath::TwoPi());
  fHPhi->Sumw2();
  fHPhi->SetXTitle("#eta");
  fHPhi->SetYTitle("#phi [radians]");
  fHPhi->SetZTitle("Events");
  fHPhi->SetDirectory(0);
  fList->Add(fHPhi);

  //

  if (!fUsePhiWeights) return;
  
  if (fOADBFileName.Length()==0)
    fOADBFileName = Form("%s/PWGLF/FORWARD/FMDEVENTPLANE/data/fmdEPoadb.root", 
     			AliAnalysisManager::GetOADBPath());

  TFile foadb(fOADBFileName);
  if(!foadb.IsOpen()) {
    AliError(Form("Cannot open OADB file %s", fOADBFileName.Data()));
    return;
  }
  fOADBContainer = (AliOADBContainer*)foadb.Get("FMDphidist");
  if (!fOADBContainer) 
    AliError(Form("No OADB container found in %s", fOADBFileName.Data()));
  foadb.Close();

  return;
}

//_____________________________________________________________________
void
AliFMDEventPlaneFinder::DefineOutput(TList* dir)
{
  // 
  // Output diagnostic histograms to directory 
  // 
  // Parameters:
  //    dir List to write in
  //  
  DGUARD(fDebug,1,"Define output of AliFMDEventPlaneFinder");
  fList = new TList();
  fList->SetOwner();
  fList->SetName(GetName());
  dir->Add(fList);

  return;
}

//____________________________________________________________________
Bool_t
AliFMDEventPlaneFinder::FindEventplane(AliVEvent* ev,
				       AliAODForwardEP& aodep,
                                       TH2D* h,
                                       AliForwardUtil::Histos* hists)
{
  // 
  // Do the calculations 
  // 
  // Parameters:
  //    hists    Histogram cache
  //    ep       calculated eventplane 
  // 
  // Return:
  //    true on successs 
  DGUARD(fDebug,1,"Find the event plane in AliFMDEventPlaneFinder");

  fQt.Set(0., 0.);
  fQa.Set(0., 0.);
  fQc.Set(0., 0.);
  fQ1.Set(0., 0.);
  fQ2.Set(0., 0.);

  TH1D epEtaHist = aodep.GetHistogram();

  fEvent = ev;
  if (hists) {
    for (UShort_t d=1; d<=3; d++) { 
      UShort_t nr = (d == 1 ? 1 : 2);
      for (UShort_t q=0; q<nr; q++) { 
        Char_t r = (q == 0 ? 'I' : 'O');
        CalcQVectors(hists->Get(d,r), &epEtaHist);
      } // for q
    } // for d
  }
  else if (h) {
    CalcQVectors(h, &epEtaHist);
  }

  aodep.SetEventplane(CalcEventplane(fQt));
  aodep.SetEventplaneA(CalcEventplane(fQa));
  aodep.SetEventplaneC(CalcEventplane(fQc));
  // aodep.SetEventplane1(CalcEventplane(fQ1));
  // aodep.SetEventplane2(CalcEventplane(fQ2));
  
  FillHists(&aodep);

  return kTRUE;
}

//_____________________________________________________________________
void
AliFMDEventPlaneFinder::CalcQVectors(TH2D* h, TH1D* eHist)
{
  //
  // Calculate the Q vectors
  //
  DGUARD(fDebug,2,"Calculate Q-vectors in AliFMDEventPlaneFinder");
  Double_t phi = 0, eta = 0, weight = 0;
  for (Int_t e = 1; e <= h->GetNbinsX(); e++) {
    Double_t qx = 0, qy = 0;
    eta = h->GetXaxis()->GetBinCenter(e);
    for (Int_t p = 1; p <= h->GetNbinsY(); p++) {
      phi = h->GetYaxis()->GetBinCenter(p);
      weight = h->GetBinContent(e, p);
      if (fUsePhiWeights) weight *= GetPhiWeight(e, p);
      // fix for FMD1 hole
      if (e > 168 && p == 14) {
	weight = h->GetBinContent(e, 4);
        if (fUsePhiWeights) weight *= GetPhiWeight(e, 4);
      }
      fHPhi->Fill(eta, phi, weight);
      // for underflowbin total Nch/eta
      fHPhi->Fill(eta, -1., weight);
      
      // increment Q vectors
      qx += weight*TMath::Cos(2.*phi);
      qy += weight*TMath::Sin(2.*phi);
    }
    TVector2 qVec = TVector2(qx, qy);
    fQt += qVec;
    if (eta < 0) fQa += qVec;
    if (eta > 0) fQc += qVec;
    // TODO: Add random ep increments
    fQeta += qVec;
    if (e % 10 == 0) {
      eHist->Fill(eta, CalcEventplane(fQeta));
      fQeta.Set(0., 0.);
    }
  }

  return;
}

//_____________________________________________________________________
Double_t
AliFMDEventPlaneFinder::CalcEventplane(const TVector2& v) const
{
  //
  // Calculate the eventplane
  //
  DGUARD(fDebug,2,"Calculate Event plane in AliFMDEventPlaneFinder");
  Double_t ep = -1;
 
  if (v.Mod() == 0.) return ep;

  ep = v.Phi();
  ep /= 2.;

  if (fDebug > 0)
    AliInfo(Form("Eventplane found to be: %f", ep));
 
  return ep;
}

//_____________________________________________________________________
Double_t 
AliFMDEventPlaneFinder::CalcDifference(Double_t a1, Double_t a2) const
{
  Double_t diff = a1 - a2; 
  if (diff <  TMath::Pi()/2) diff = TMath::Pi() - diff;
  if (diff >= TMath::Pi()/2) diff = TMath::Pi() - diff;
  return diff;
}

//_____________________________________________________________________
void 
AliFMDEventPlaneFinder::FillHists(AliAODForwardEP* fmdEP)
{
  //
  // Fill diagnostics histograms
  //
  DGUARD(fDebug,2,"Fill histograms in AliFMDEventPlaneFinder");
  Double_t fmdEPt = fmdEP->GetEventplane();
  Double_t fmdEPa = fmdEP->GetEventplaneA();
  Double_t fmdEPc = fmdEP->GetEventplaneC();
  // Double_t fmdEP1 = fmdEP->GetEventplane1();
  // Double_t fmdEP2 = fmdEP->GetEventplane2();

  // FMD hists
  fHepFMD->Fill(fmdEPt);
  fHepFMDA->Fill(fmdEPa);
  fHepFMDC->Fill(fmdEPc);

  fHdiffFMDAC->Fill(CalcDifference(fmdEPa,fmdEPc));
  fHcorrFMDAC->Fill(fmdEPa, fmdEPc);
  // fHepFMDQC1->Fill(fmdEP1);
  // fHepFMDQC2->Fill(fmdEP2);

  // TPC comparison
  AliEventplane* ep = fEvent->GetEventplane();
  Double_t tpcEP   = (ep ? ep->GetEventplane("Q") : -1);
  fHdiffFMDTPC->Fill(CalcDifference(fmdEPt, tpcEP));
  fHcorrFMDTPC->Fill(fmdEPt, tpcEP);

  // VZERO comparison
  Double_t vzeroEP = ep ? ep->GetEventplane("V0", fEvent, 2) : -1;
  fHdiffFMDVZERO->Fill(CalcDifference(fmdEPt, vzeroEP));
  fHcorrFMDVZERO->Fill(fmdEPt, vzeroEP);
  
  return;
}

//_____________________________________________________________________
Double_t 
AliFMDEventPlaneFinder::GetPhiWeight(Int_t etaBin, Int_t phiBin) const
{
  //
  // Get phi weight for flattening
  //
  Double_t phiWeight = 1;
  
  if (!fPhiDist) return phiWeight;

  Double_t nParticles   = fPhiDist->GetBinContent(etaBin, 0);
  Double_t nPhiBins     = fPhiDist->GetNbinsY();
  Double_t phiDistValue = fPhiDist->GetBinContent(etaBin, phiBin);
  
  if (phiDistValue > 0) phiWeight = nParticles/nPhiBins/phiDistValue;

  return phiWeight;
}

//____________________________________________________________________
void
AliFMDEventPlaneFinder::SetRunNumber(Int_t run)
{
  // 
  // Set run number
  //
  if (fRunNumber == run) return;
  
  fRunNumber = run;
  if (fUsePhiWeights) GetPhiDist();

  return;
}

//____________________________________________________________________
void
AliFMDEventPlaneFinder::GetPhiDist()
{
  //
  // Get phi dist from OADB
  //
  if (!fOADBContainer) return;

  fPhiDist = static_cast<TH2D*>(fOADBContainer
				->GetObject(fRunNumber, "Default")
				->Clone(Form("fPhiDist_%d", fRunNumber)));
  fList->Add(fPhiDist);

  return;
}

//____________________________________________________________________
void
AliFMDEventPlaneFinder::Print(Option_t* /*option*/) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  char ind[gROOT->GetDirLevel()+3];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << ClassName() << ": " << GetName() << '\n'
	    << std::boolalpha 
	    << ind << "Eventplane finder active!" << '\n'
	    << ind << "Loading OADB object for run number: " 
	    << fRunNumber << '\n'
	    << ind << "Using phi weights: " << fUsePhiWeights  << '\n'
	    << std::noboolalpha
	    << std::endl;
}
//____________________________________________________________________
//
// EOF
//
	  


