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
    fHistFMDEventplane(0),
    fHistFMDEventplaneA(0),
    fHistFMDEventplaneC(0),
    fHistFMDEventplane1(0),
    fHistFMDEventplane2(0),
    fHistPhiDist(0),
    fHistFMDmTPCep(0),
    fHistFMDvsTPCep(0),
    fHistFMDmVZEROep(0),
    fHistFMDvsVZEROep(0),
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
    fHistFMDEventplane(0),
    fHistFMDEventplaneA(0),
    fHistFMDEventplaneC(0),
    fHistFMDEventplane1(0),
    fHistFMDEventplane2(0),
    fHistPhiDist(0),
    fHistFMDmTPCep(0),
    fHistFMDvsTPCep(0),
    fHistFMDmVZEROep(0),
    fHistFMDvsVZEROep(0),
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
    fHistFMDEventplane(o.fHistFMDEventplane),
    fHistFMDEventplaneA(o.fHistFMDEventplaneA),
    fHistFMDEventplaneC(o.fHistFMDEventplaneC),
    fHistFMDEventplane1(o.fHistFMDEventplane1),
    fHistFMDEventplane2(o.fHistFMDEventplane2),
    fHistPhiDist(o.fHistPhiDist),
    fHistFMDmTPCep(o.fHistFMDmTPCep),
    fHistFMDvsTPCep(o.fHistFMDvsTPCep),
    fHistFMDmVZEROep(o.fHistFMDmVZEROep),
    fHistFMDvsVZEROep(o.fHistFMDvsVZEROep),
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

  fList               = o.fList;
  fEvent              = o.fEvent;
  fQt                 = o.fQt;
  fQa                 = o.fQa;
  fQc                 = o.fQc;
  fQ1                 = o.fQ1;
  fQ2                 = o.fQ2;
  fQeta               = o.fQeta;
  fHistFMDEventplane  = o.fHistFMDEventplane;
  fHistFMDEventplaneA = o.fHistFMDEventplaneA;
  fHistFMDEventplaneC = o.fHistFMDEventplaneC;
  fHistFMDEventplane1 = o.fHistFMDEventplane1;
  fHistFMDEventplane2 = o.fHistFMDEventplane2;
  fHistPhiDist        = o.fHistPhiDist;
  fHistFMDmTPCep      = o.fHistFMDmTPCep;
  fHistFMDvsTPCep     = o.fHistFMDvsTPCep;
  fHistFMDmVZEROep    = o.fHistFMDmVZEROep;
  fHistFMDvsVZEROep   = o.fHistFMDvsVZEROep;
  fDebug              = o.fDebug;
  fOADBFileName       = o.fOADBFileName;
  fOADBContainer      = o.fOADBContainer;
  fPhiDist            = o.fPhiDist;
  fRunNumber          = o.fRunNumber;
  fUsePhiWeights      = o.fUsePhiWeights;

  return *this;
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
  fHistFMDEventplane = new TH1D("hFMDEventplane", "FMD eventplane", 
                                100, 0., TMath::Pi());
  fHistFMDEventplane->Sumw2();
  fHistFMDEventplane->SetDirectory(0);
  fList->Add(fHistFMDEventplane);

  //
  fHistFMDEventplaneA = new TH1D("hFMDEventplaneA", "FMD eventplaneA", 
                                100, 0., TMath::Pi());
  fHistFMDEventplaneA->Sumw2();
  fHistFMDEventplaneA->SetDirectory(0);
  fList->Add(fHistFMDEventplaneA);

  //
  fHistFMDEventplaneC = new TH1D("hFMDEventplaneC", "FMD eventplaneC", 
                                100, 0., TMath::Pi());
  fHistFMDEventplaneC->Sumw2();
  fHistFMDEventplaneC->SetDirectory(0);
  fList->Add(fHistFMDEventplaneC);

  //
  fHistFMDEventplane1 = new TH1D("hFMDEventplane1", "FMD eventplane1", 
                                100, 0., TMath::Pi());
  fHistFMDEventplane1->Sumw2();
  fHistFMDEventplane1->SetDirectory(0);
  fList->Add(fHistFMDEventplane1);

  //
  fHistFMDEventplane2 = new TH1D("hFMDEventplane2", "FMD eventplane2", 
                                100, 0., TMath::Pi());
  fHistFMDEventplane2->Sumw2();
  fHistFMDEventplane2->SetDirectory(0);
  fList->Add(fHistFMDEventplane2);

  //
  fHistPhiDist = new TH2D("hPhiDist", "Phi distribution in FMD",
	etaAxis.GetNbins(), etaAxis.GetXmin(), etaAxis.GetXmax(), 
	20, 0., TMath::TwoPi());
  fHistPhiDist->Sumw2();
  fHistPhiDist->SetDirectory(0);
  fList->Add(fHistPhiDist);

  //
  fHistFMDmTPCep = new TH1D("hFMDmTPCep", 
                             "Eventplane from FMD - TPC",
                             100, -TMath::Pi()/2., TMath::Pi()/2.);
  fHistFMDmTPCep->Sumw2();
  fHistFMDmTPCep->SetDirectory(0);
  fList->Add(fHistFMDmTPCep);

  //
  fHistFMDvsTPCep = new TH2F("hFMDvsTPCep", 
                             "Eventplane from FMD vs. TPC",
                             100, 0., TMath::Pi(),
                             100, 0., TMath::Pi());
  fHistFMDvsTPCep->Sumw2();
  fHistFMDvsTPCep->SetDirectory(0);
  fList->Add(fHistFMDvsTPCep);

  //
  fHistFMDmVZEROep = new TH1D("hFMDmVZEROep", 
                             "Eventplane from FMD - VZERO",
                             100, -TMath::Pi()/2., TMath::Pi()/2.);
  fHistFMDmVZEROep->Sumw2();
  fHistFMDmVZEROep->SetDirectory(0);
  fList->Add(fHistFMDmVZEROep);

  //
  fHistFMDvsVZEROep = new TH2F("hFMDvsVZEROep", 
                             "Eventplane from FMD vs. VZERO",
                             100, 0., TMath::Pi(),
                             100, 0., TMath::Pi());
  fHistFMDvsVZEROep->Sumw2();
  fHistFMDvsVZEROep->SetDirectory(0);
  fList->Add(fHistFMDvsVZEROep);

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
  if (!fOADBContainer) AliError(Form("No OADB container found in %s", fOADBFileName.Data()));
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
      fHistPhiDist->Fill(eta, phi, weight);
      // for underflowbin total Nch/eta
      fHistPhiDist->Fill(eta, -1., weight);
      
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
AliFMDEventPlaneFinder::CalcEventplane(TVector2 v) const
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
  fHistFMDEventplane->Fill(fmdEPt);
  fHistFMDEventplaneA->Fill(fmdEPa);
  fHistFMDEventplaneC->Fill(fmdEPc);
//  fHistFMDEventplane1->Fill(fmdEP1);
//  fHistFMDEventplane2->Fill(fmdEP2);

  // TPC comparison
  AliEventplane* ep = fEvent->GetEventplane();
  Double_t tpcEP = -1;
  if (ep) tpcEP = ep->GetEventplane("Q");

  Double_t diffTPC = fmdEPt - tpcEP;

  if (diffTPC <  TMath::Pi()/2.) diffTPC = TMath::Pi() - diffTPC;
  if (diffTPC >= TMath::Pi()/2.) diffTPC = TMath::Pi() - diffTPC;

  fHistFMDmTPCep->Fill(diffTPC);

//  if (fmdEPt >= TMath::Pi()/2. && fmdEPt - tpcEP >= TMath::Pi()/2.) tpcEP = TMath::Pi()-tpcEP;
//  if (fmdEPt <  TMath::Pi()/2. && tpcEP - fmdEPt >= TMath::Pi()/2.) tpcEP = TMath::Pi()-tpcEP;

  fHistFMDvsTPCep->Fill(fmdEPt, tpcEP);

  // VZERO comparison
  Double_t vzeroEP = -1;
  if (ep) vzeroEP = ep->GetEventplane("V0", fEvent, 2);

  Double_t diffVZERO = fmdEPt - vzeroEP;

  if (diffVZERO <  TMath::Pi()/2.) diffVZERO = TMath::Pi() - diffVZERO;
  if (diffVZERO >= TMath::Pi()/2.) diffVZERO = TMath::Pi() - diffVZERO;

  fHistFMDmVZEROep->Fill(diffVZERO);

//  if (fmdEPt >= TMath::Pi()/2. && fmdEPt - vzeroEP >= TMath::Pi()/2.) vzeroEP = TMath::Pi()-vzeroEP;
//  if (fmdEPt <  TMath::Pi()/2. && vzeroEP - fmdEPt >= TMath::Pi()/2.) vzeroEP = TMath::Pi()-vzeroEP;

  fHistFMDvsVZEROep->Fill(fmdEPt, vzeroEP);
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
  
  if (fPhiDist) {
    Double_t nParticles = fPhiDist->GetBinContent(etaBin, 0);
    Double_t nPhiBins = fPhiDist->GetNbinsY();
    Double_t phiDistValue = fPhiDist->GetBinContent(etaBin, phiBin);
    
    if (phiDistValue > 0) phiWeight = nParticles/nPhiBins/phiDistValue;
  }

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

  fPhiDist = (TH2D*)fOADBContainer->GetObject(fRunNumber, "Default")->Clone(Form("fPhiDist_%d", fRunNumber));
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
	    << ind << "Loading OADB object for run number: " << fRunNumber << '\n'
	    << ind << "Using phi weights: " << (fUsePhiWeights ? "true" : "false") << '\n'
	    << std::noboolalpha
	    << std::flush;
  std::cout << std::endl;
}
//____________________________________________________________________
//
// EOF
//
	  


