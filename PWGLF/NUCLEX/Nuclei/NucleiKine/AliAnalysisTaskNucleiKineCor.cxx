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

using std::cout;
using std::endl;

AliAnalysisTaskNucleiKineCor::AliAnalysisTaskNucleiKineCor(const char* name) :
  AliAnalysisTaskSE{name},
  fPdgCodes{211, -211, 321, -321, 2212, -2212, 2112, -2112, 1000010020, -1000010020,3122,-3122,3312,-3312,3334,-3334},
  fParticleNames{"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}", "n", "#bar{n}", "d", "#bar{d}",
      "#Lambda", "#bar{#Lambda}", "#Xi^{+}", "#Xi^{-}", "#Omega^{+}", "#Omega^{-}"},
  fPt(5),
  fP0(0),
  fAcc(1),
  fCutY(0),
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
  fHists[3] = new TH1D("hEtaHadrons",";#eta; Events",100,-fAcc,fAcc);
  fHists[3]->Sumw2();
  fOutputList->Add(fHists[3]);

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

  fHists[18] = new TH3D("hHadHadAcc","Acceptance correction;#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,128,-fAcc,fAcc,100,0.,10);
  fHists[18]->Sumw2();
  fOutputList->Add(fHists[18]);
  fHists[19] = new TH2D("hHadHadwAcc","Hadron-hadron correlation;#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,100,0.,10);
  fHists[19]->Sumw2();
  fOutputList->Add(fHists[19]);

  fHists[20] = new TH2D("hHadEta","Hadron;#eta;Associated p_{T} (GeV/c)",
			128,-fAcc,fAcc,100,0.,10);
  fHists[20]->Sumw2();

  fOutputList->Add(fHists[20]);

  fHists[21] = new TH2D("hProtEta","Proton correlation;#eta;Associated p_{T} (GeV/c)",
			128,-fAcc,fAcc,100,0.,10);
  fHists[21]->Sumw2();
  fOutputList->Add(fHists[21]);

  fHists[22] = new TH2D("hDeutEta","Deuteron correlation;#eta;Associated p_{T} (GeV/c)",
			128,-fAcc,fAcc,100,0.,10);
  fHists[22]->Sumw2();
  fOutputList->Add(fHists[22]);

  fHists[23] = new TH2D("heHadHad","Hadron-hadron correlation (half #eta);#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,100,0.,10);
  fHists[23]->Sumw2();
  fOutputList->Add(fHists[23]);

  fHists[24] = new TH2D("heHadProt","Hadron-proton correlation (half #eta);#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,100,0.,10);
  fHists[24]->Sumw2();
  fOutputList->Add(fHists[24]);

  fHists[25] = new TH2D("heHadDeut","Hadron-deuteron correlation (half #eta);#Delta#varphi;Associated p_{T} (GeV/c)",
			128,-TMath::Pi()/2,3*TMath::Pi()/2,100,0.,10);
  fHists[25]->Sumw2();
  fOutputList->Add(fHists[25]);

  fHists[26] = new TH1D("feHadTrigs",";Number; Entries",25,-0.5,24.5);
  fHists[26]->Sumw2();
  fOutputList->Add(fHists[26]);

  fHists[97] = new TProfile("hXsec",";Cross section;;",1,-0.5,0.5,"e");
  fOutputList->Add(fHists[97]);

  fHists[98] = new TProfile("hTrial",";Trials;;",1,-0.5,0.5,"e");
  fOutputList->Add(fHists[98]);

  fHists[99] = new TH1D("hStats","Stats;Criterium;Entries",10,-0.5,9.5);
  fHists[99]->Sumw2();
  fOutputList->Add(fHists[99]);

  if (fP0>0) {
    fAfterBurner.SetNucleusPdgCode(1000010020);
    fAfterBurner.SetCoalescenceMomentum(fP0);
  }
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

  Int_t nstack = stack->GetNtrack();
  fEvents->Fill(nstack);

  // find nucleons and anti-nucleons
  TObjArray keepParticles;
  if (fP0>0) {
    for (Int_t i=0; i < stack->GetNprimary(); ++i) {
      TParticle* iParticle = stack->Particle(i);
      if (iParticle->GetStatusCode() != 1)
	continue;
      switch (iParticle->GetPdgCode()) {
      case 2212: //kProton:
	keepParticles.Add(iParticle);
	break;
      case 2112: //kNeutron:
	keepParticles.Add(iParticle);
	break;
      case -2212: //kProtonBar:
	keepParticles.Add(iParticle);
	break;
      case -2112: //kNeutronBar:
	keepParticles.Add(iParticle);
	break;
      default:
	break;
      }
    }
    fAfterBurner.SetStack(stack);
    fAfterBurner.Generate();
  }

  TParticle *leadP = 0;
  Double_t ptl     = 0;
  TObjArray arrh;
  TObjArray arrp;
  TObjArray arrd;
  Int_t nstack2 = stack->GetNtrack();

  for (int iTracks = 0; iTracks < nstack2; ++iTracks) {
    TParticle* track = stack->Particle(iTracks);
    if (!track) 
      continue;
    if (!stack->IsPhysicalPrimary(iTracks)) 
      continue;
    if (fCutY) {
      if (TMath::Abs(track->Y()) > fAcc) 
	continue;
    } else {
      if (TMath::Abs(track->Eta()) > fAcc) 
	continue;
    }
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

  if (leadP!=0) 
    fPtLead->Fill(leadP->Pt());

  Int_t nh = arrh.GetEntries();
  fHists[0]->Fill(nh);
  Int_t np = arrp.GetEntries();
  fHists[1]->Fill(np);
  Int_t nd = arrd.GetEntries();
  fHists[2]->Fill(nd);

  // hadron - hadron
  TH2 *hh  = (TH2*)fHists[10];
  TH3 *ac  = (TH3*)fHists[18];
  TH2 *hh2 = (TH2*)fHists[19];
  TH2 *hhe = (TH2*)fHists[20];
  TH2 *hh3 = (TH2*)fHists[23];
  Int_t ntrigs  = 0;
  Int_t ntrigs2 = 0;
  for (Int_t i=0;i<nh;++i) {
    TParticle *part1 = (TParticle*)arrh.At(i);
    const double pt1 = part1->Pt();
    if (pt1<fPt)
      continue;
    const double phi1 = part1->Phi();
    const double eta1 = part1->Eta();
    fHists[3]->Fill(eta1);
    ntrigs++;
    if (TMath::Abs(eta1)<fAcc/2.)
      ntrigs2++;
    for (Int_t j=0;j<nh;++j) {
      TParticle *part2 = (TParticle*)arrh.At(j);
      if (part2==part1) 
	continue;
      const double dphi = DeltaPhi(phi1,part2->Phi());
      const double deta = part2->Eta()-eta1;
      hh->Fill(dphi,part2->Pt());
      hhe->Fill(part2->Eta(),part2->Pt());
      const double acorr = 1-TMath::Abs(deta)/2/fAcc; 
      if (acorr>0) {
	ac->Fill(dphi,deta,part2->Pt());
	hh2->Fill(dphi,part2->Pt(),1./acorr);
      }
      if (TMath::Abs(eta1)<fAcc/2.)
	hh3->Fill(dphi,part2->Pt());
    }
    fHists[11]->Fill(ntrigs);  
    fHists[26]->Fill(ntrigs2);  
  }

  // leading hadron - hadron
  hh = (TH2*)fHists[12];
  hhe = (TH2*)fHists[21];
  if (ptl>fPt) {
    const double phi1 = leadP->Phi();
    // const double eta1 = leadP->Eta();
    for (Int_t j=0;j<nh;++j) {
      TParticle *part2 = (TParticle*)arrh.At(j);
      if (part2==leadP) 
	continue;
      const double dphi = DeltaPhi(phi1,part2->Phi());
      // const double deta = part2->Eta()-eta1;
      hh->Fill(dphi,part2->Pt());
      hhe->Fill(part2->Eta(),part2->Pt());
    }
    fHists[13]->Fill(1);  
  }

  // hadron - proton
  hh = (TH2*)fHists[14];
  hh3 = (TH2*)fHists[24];
  for (Int_t i=0;i<nh;++i) {
    TParticle *part1 = (TParticle*)arrh.At(i);
    const double pt1 = part1->Pt();
    if (pt1<fPt)
      continue;
    const double phi1 = part1->Phi();
    const double eta1 = part1->Eta();
    for (Int_t j=0;j<np;++j) {
      TParticle *part2 = (TParticle*)arrp.At(j);
      if (part2==part1) 
	continue;
      const double dphi = DeltaPhi(phi1,part2->Phi());
      // const double deta = part2->Eta()-eta1;
      hh->Fill(dphi,part2->Pt());
      if (TMath::Abs(eta1)<fAcc/2.)
	hh3->Fill(dphi,part2->Pt());
    }
  }

  // leading hadron - proton
  hh = (TH2*)fHists[15];
  if (ptl>fPt) {
    const double phi1 = leadP->Phi();
    // const double eta1 = leadP->Eta();
    for (Int_t j=0;j<np;++j) {
      TParticle *part2 = (TParticle*)arrp.At(j);
      if (part2==leadP) 
	continue;
      const double dphi = DeltaPhi(phi1,part2->Phi());
      // const double deta = part2->Eta()-eta1;
      hh->Fill(dphi,part2->Pt());
    }
  }

  // hadron - deuteron
  hh = (TH2*)fHists[16];
  hhe = (TH2*)fHists[22];
  hh3 = (TH2*)fHists[25];
  for (Int_t i=0;i<nh;++i) {
    TParticle *part1 = (TParticle*)arrh.At(i);
    const double pt1 = part1->Pt();
    if (pt1<fPt)
      continue;
    const double phi1 = part1->Phi();
    const double eta1 = part1->Eta();
    for (Int_t j=0;j<nd;++j) {
      TParticle *part2 = (TParticle*)arrd.At(j);
      if (part2==part1) 
	continue;
      const double dphi = DeltaPhi(phi1,part2->Phi());
      // const double deta = part2->Eta()-eta1;
      hh->Fill(dphi,part2->Pt());
      hhe->Fill(part2->Eta(),part2->Pt());
      if (TMath::Abs(eta1)<fAcc/2.)
	hh3->Fill(dphi,part2->Pt());
    }
  }

  // leading hadron - deuteron
  hh = (TH2*)fHists[17];
  if (ptl>fPt) {
    const double phi1 = leadP->Phi();
    // const double eta1 = leadP->Eta();
    for (Int_t j=0;j<nd;++j) {
      TParticle *part2 = (TParticle*)arrd.At(j);
      if (part2==leadP) 
	continue;
      const double dphi = DeltaPhi(phi1,part2->Phi());
      // const double deta = part2->Eta()-eta1;
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

  AliGenPythiaEventHeader *hp = 0;
  if (fP0>0) {
    hp = dynamic_cast<AliGenPythiaEventHeader*>(mcEvent->GenEventHeader());
  } else {
    AliGenCocktailEventHeader *hd = dynamic_cast<AliGenCocktailEventHeader*>(mcEvent->GenEventHeader());
    if (!hd) 
      AliFatal("Missing cocktail header");
    hp = dynamic_cast<AliGenPythiaEventHeader*>(hd->GetHeaders()->At(0));
  }
  if (hp) {
    Double_t xsec = hp->GetXsection();
    Int_t trials = hp->Trials();
    TProfile* pr = (TProfile*)fHists[97];
    pr->Fill(0.,xsec);
    pr = (TProfile*)fHists[98];
    pr->Fill(0.,trials);
  }

  // reset nucleons and anti-nucleons
  if (fP0>0) {
    for (Int_t i=0; i<keepParticles.GetEntries(); ++i) {
      TParticle *p = dynamic_cast<TParticle*>(keepParticles.At(i));
      if (!p) 
	continue;
      p->SetStatusCode(1);
    }
    for (Int_t i=nstack;i<nstack2;++i) {
      TParticle* iParticle = stack->Particle(i);
      if (iParticle->GetPdgCode()==1000010020) {
	iParticle->SetPdgCode(0);
	iParticle->SetStatusCode(201);
      }
    }
    stack->SetNtrack(nstack);
    stack->SetHighWaterMark(nstack);
  }

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
