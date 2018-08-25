#include "AliAnalysisTaskNucleiKineCor.h"
#include <Riostream.h>
#include <AliAnalysisManager.h>
#include <AliCollisionGeometry.h>
#include <AliGenCocktailEventHeader.h>
#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliHeader.h>
#include <AliInputEventHandler.h>
#include <AliLog.h>
#include <AliMCEvent.h>
#include <AliMCEventHandler.h>
#include <AliStack.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TList.h>
#include <TMath.h>
#include <TProfile.h>

ClassImp(AliAnalysisTaskNucleiKineCor);

AliAnalysisTaskNucleiKineCor::AliAnalysisTaskNucleiKineCor(const char* name) :
  AliAnalysisTaskSE{name},
  fPdgCodes{211, -211, 321, -321, 2212, -2212, 2112, -2112, 1000010020, -1000010020,3122,-3122,3312,-3312,3334,-3334},
  fParticleNames{"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}", "n", "#bar{n}", "d", "#bar{d}",
      "#Lambda", "#bar{#Lambda}", "#Xi^{+}", "#Xi^{-}", "#Omega^{+}", "#Omega^{-}"},
  fPt(5),
  fOutputList(0)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskNucleiKineCor::~AliAnalysisTaskNucleiKineCor() 
{
  if (fOutputList) 
    delete fOutputList;
}

void AliAnalysisTaskNucleiKineCor::UserCreateOutputObjects() 
{
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fEvents = new TH1D("fEventCounter",";Number of particles on stack; Events",2500,-0.5,2499.5);
  fOutputList->Add(fEvents);

  fPtHist = new TH1D("fPtHadrons",";p_{T} (GeV); Entries",250,0,25);
  fPtHist->Sumw2();
  fOutputList->Add(fPtHist);
  fPtLead = new TH1D("fPtLeadHadrons",";p_{T} (GeV); Entries",250,0,25);
  fPtLead->Sumw2();
  fOutputList->Add(fPtLead);

  fPtPro = new TH1D("fPtProtons",";p_{T} (GeV); Entries",250,0,25);
  fPtPro->Sumw2();
  fOutputList->Add(fPtPro);
  fPtDeu = new TH1D("fPtDeuterons",";p_{T} (GeV); Entries",250,0,25);
  fPtDeu->Sumw2();
  fOutputList->Add(fPtDeu);

  fHists[0] = new TH1D("hNhadrons",";Number of hadrons at |Y|<1; Events",250,-0.5,249.5);
  fHists[0]->Sumw2();
  fOutputList->Add(fHists[0]);
  fHists[1] = new TH1D("hNprotons",";Number of protons at |Y|<1; Events",50,-0.5,49.5);
  fHists[1]->Sumw2();
  fOutputList->Add(fHists[1]);
  fHists[2] = new TH1D("hNdeuterons",";Number of deuterons at |Y|<1; Events",15,-0.5,14.5);
  fHists[2]->Sumw2();
  fOutputList->Add(fHists[2]);

  fHists[10] = new TH2D("hHadHad","Hadron-hadron correlation;#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,100,0.,10);
  fHists[10]->Sumw2();
  fOutputList->Add(fHists[10]);
  fHists[11] = new TH1D("fHadTrigs",";Number; Entries",25,-0.5,24.5);
  fHists[11]->Sumw2();
  fOutputList->Add(fHists[11]);
  fHists[12] = new TH2D("hHadLeadHad","Leading Hadron-hadron correlation;#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,100,0.,10);
  fHists[12]->Sumw2();
  fOutputList->Add(fHists[12]);
  fHists[13] = new TH1D("fLeadHadTrigs",";Number; Entries",2,-0.5,1.5);
  fHists[13]->Sumw2();
  fOutputList->Add(fHists[13]);

  fHists[14] = new TH2D("hHadProt","Hadron-proton correlation;#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,100,0.,10);
  fHists[14]->Sumw2();
  fOutputList->Add(fHists[14]);
  fHists[15] = new TH2D("hHadLeadProt","Leading Hadron-proton correlation;#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,100,0.,10);
  fHists[15]->Sumw2();
  fOutputList->Add(fHists[15]);

  fHists[16] = new TH2D("hHadDeut","Hadron-deuteron correlation;#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,100,0.,10);
  fHists[16]->Sumw2();
  fOutputList->Add(fHists[16]);
  fHists[17] = new TH2D("hHadLeadDeut","Leading Hadron-deuteron correlation;#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,100,0.,10);
  fHists[17]->Sumw2();
  fOutputList->Add(fHists[17]);

  fHists[97] = new TProfile("hXsec",";Cross section;;",1,-0.5,0.5,"e");
  fOutputList->Add(fHists[97]);

  fHists[98] = new TProfile("hTrial",";Trials;;",1,-0.5,0.5,"e");
  fOutputList->Add(fHists[98]);

  fHists[99] = new TH1D("hStats","Stats;Criterium;Entries",10,-0.5,9.5);
  fHists[99]->Sumw2();
  fOutputList->Add(fHists[99]);

  PostData(1, fOutputList);
}

void AliAnalysisTaskNucleiKineCor::UserExec(Option_t*) 
{
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent)
    AliFatal("Missing MC event!");

  AliStack* stack = mcEvent->Stack();
  if (!stack)
    AliFatal("Missing stack.");

  int nstack = stack->GetNtrack();
  fEvents->Fill(nstack);

  TParticle *leadP = 0;
  Double_t ptl     = 0;
  TObjArray arrh;
  TObjArray arrp;
  TObjArray arrd;
  for (int iTracks = 0; iTracks < nstack; ++iTracks) {
    TParticle* track = stack->Particle(iTracks);
    if (!track) 
      continue;
    if (!stack->IsPhysicalPrimary(iTracks)) 
      continue;
    if (TMath::Abs(track->Y()) <= 1) 
      continue;
    const int pdg = track->GetPdgCode();
    const int apg = TMath::Abs(pdg);
    const double pt = track->Pt();
    if ((apg==211)||(apg==321)||(apg==2212)) {
      fPtHist->Fill(pt);
      if (pt>ptl) {
	leadP = track;
	ptl = pt;
      }
      arrh.Add(track);
    }
    if (apg==2212) {
      fPtPro->Fill(pt);
      arrp.Add(track);
    } 
    if (apg==1000010020) {
      fPtDeu->Fill(pt);
      arrd.Add(track);
    } 
  }
  fPtLead->Fill(leadP->Pt());

  Int_t nh = arrh.GetEntries();
  fHists[0]->Fill(nh);
  Int_t np = arrp.GetEntries();
  fHists[1]->Fill(np);
  Int_t nd = arrd.GetEntries();
  fHists[2]->Fill(nd);

  // hadron - hadron
  TH2 *hh = (TH2*)fHists[10];
  Int_t ntrigs = 0;
  for (Int_t i=0;i<nh;++i) {
    TParticle *part1 = (TParticle*)arrh.At(i);
    const double pt1 = part1->Pt();
    if (pt1<fPt)
      continue;
    const double phi1 = part1->Phi();
    ntrigs++;
    for (Int_t j=0;j<nh;++j) {
      TParticle *part2 = (TParticle*)arrh.At(j);
      if (part2==part1) 
	continue;
      double dphi = DeltaPhi(phi1,part2->Phi());
      hh->Fill(dphi,part2->Pt());
    }
    fHists[11]->Fill(ntrigs);  
  }

  // leading hadron - hadron
  hh = (TH2*)fHists[12];
  if (ptl>fPt) {
    Double_t phi1 = leadP->Phi();
    for (Int_t j=0;j<nh;++j) {
      TParticle *part2 = (TParticle*)arrh.At(j);
      if (part2==leadP) 
	continue;
      double dphi = DeltaPhi(phi1,part2->Phi());
      hh->Fill(dphi,part2->Pt());
    }
    fHists[13]->Fill(1);  
  }

  // hadron - proton
  hh = (TH2*)fHists[14];
  for (Int_t i=0;i<nh;++i) {
    TParticle *part1 = (TParticle*)arrh.At(i);
    const double pt1 = part1->Pt();
    if (pt1<fPt)
      continue;
    const double phi1 = part1->Phi();
    for (Int_t j=0;j<np;++j) {
      TParticle *part2 = (TParticle*)arrp.At(j);
      if (part2==part1) 
	continue;
      double dphi = DeltaPhi(phi1,part2->Phi());
      hh->Fill(dphi,part2->Pt());
    }
  }

  // leading hadron - proton
  hh = (TH2*)fHists[15];
  if (ptl>fPt) {
    Double_t phi1 = leadP->Phi();
    for (Int_t j=0;j<np;++j) {
      TParticle *part2 = (TParticle*)arrp.At(j);
      if (part2==leadP) 
	continue;
      double dphi = DeltaPhi(phi1,part2->Phi());
      hh->Fill(dphi,part2->Pt());
    }
  }

  // hadron - deuteron
  hh = (TH2*)fHists[16];
  for (Int_t i=0;i<nh;++i) {
    TParticle *part1 = (TParticle*)arrh.At(i);
    const double pt1 = part1->Pt();
    if (pt1<fPt)
      continue;
    const double phi1 = part1->Phi();
    for (Int_t j=0;j<nd;++j) {
      TParticle *part2 = (TParticle*)arrd.At(j);
      if (part2==part1) 
	continue;
      double dphi = DeltaPhi(phi1,part2->Phi());
      hh->Fill(dphi,part2->Pt());
    }
  }

  // leading hadron - deuteron
  hh = (TH2*)fHists[17];
  if (ptl>fPt) {
    Double_t phi1 = leadP->Phi();
    for (Int_t j=0;j<nd;++j) {
      TParticle *part2 = (TParticle*)arrd.At(j);
      if (part2==leadP) 
	continue;
      double dphi = DeltaPhi(phi1,part2->Phi());
      hh->Fill(dphi,part2->Pt());
    }
  }

  fHists[99]->Fill(0);  
  if (ntrigs>0)
    fHists[99]->Fill(1);  
  if (nd>0)
    fHists[99]->Fill(2);  
  if ((nd>0)&&(ntrigs>0))
    fHists[99]->Fill(3);  

  AliGenCocktailEventHeader *hd = dynamic_cast<AliGenCocktailEventHeader*>(mcEvent->GenEventHeader());
  if (!hd) 
    AliFatal("Missing cocktail header");
  AliGenPythiaEventHeader *hp = dynamic_cast<AliGenPythiaEventHeader*>(hd->GetHeaders()->At(0));
  Double_t xsec = hp->GetXsection();
  Int_t trials = hp->Trials();
  TProfile* pr = (TProfile*)fHists[97];
  pr->Fill(0.,xsec);
  pr = (TProfile*)fHists[98];
  pr->Fill(0.,trials);

  PostData(1, fOutputList);
}

Double_t AliAnalysisTaskNucleiKineCor::DeltaPhi(Double_t phia, Double_t phib, Double_t rangeMin, Double_t rangeMax) const
{
  Double_t dphi = -999;
  const Double_t pi = TMath::Pi();
  
  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += 2*pi;
  else if (dphi > rangeMax) dphi -= 2*pi;
  
  return dphi;
}
