//*-- Author: Aleksei Pavlinov(WSU) 
//            8-feb-2002 new version with AliEMCALJetMicroDst
#include "anaDst.h"
#include "jetDst.h"
#if !defined(__CINT__) || defined(__MAKECINT__)
// These files are needed for compilation only
// (see directory /auto/u/pavlinov/macros on PDSF)
// Using with caution if you are not PAI
#include "macroIncludePai.h"
#include "macroIncludeAlice.h"
extern "C++" {void loadlibs();}
#endif

class TH1F; 
TH1F* hJetEt, *hJetEta, *hJetPhi, *hNtracksInJet, *hPtTracksInJet, *hTrackDR, *hNjet;

class TH2F; TH2F* hPartPt, *hPartEta, *hPartPhi;
TH1F* hPartDiffPt, *hPartDiffEta, * hPartDiffPhi;

TH1F *hJetPartAng, *hDeta, *hDphi, *hDpt; // Comparing between jet and nearest parton

class TCanvas; TCanvas *c1, *c2, *c3;
class TPad;    TPad  *p1;
class TList;   TList   *lPartons, *lJets, *lCompPJ;

class TVector3; 

Int_t nparticles, nJet;
Double_t phiMax = TMath::Pi()*120./180.;
Double_t rJet = 0.5; // this default

TFile *file;
class AliEMCALJetMicroDst; AliEMCALJetMicroDst *jDst=0;
Int_t nevMax = 1000;
TString fname;
// fitting
class TF1 *g1;// gauss

void 
anaDst(Int_t mode) 
{
// Dynamically link some shared libs                    
  if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("../macros/loadlibs.C");
      loadlibs();
  }

//partons after HSC and jet in acceptance.
  TVector3 *vJet;
  TVector3 *vPart1, *vPart2;
  vJet   = new TVector3;
  vPart1 = new TVector3;
  vPart2 = new TVector3;

// for getting information about nFileOut 
//#if !defined(__CINT__)
//  gROOT->ProcessLine(".L jetDst.C+");
//#else
//  gROOT->LoadMacro("jetDst.C");
//#endif
  if(!defineMicroDst(mode)) return;    

  bookHist1();
  bookHistPartonsAfterHardSc();

    Float_t phiT[50], etaT[50], ptT[50];

    Int_t nbytes=0, nb=0, nentries = Int_t(jDst->GetTree()->GetEntries());
    if(nevMax>nentries) nevMax = nentries;
    for (int nev=0; nev< nevMax; nev++) {
	nb = jDst->GetEntry(nev);
        nbytes += nb;
	printf("\n Event .............%d bytes %d\n", nev, nb);
	//        jDst->Print();
        if(nb<0) break; // last event or something wrong
	//        return;
    
        Bool_t inEmc = fillInfoForPartons(vPart1, vPart2);
	//        continue;

	nJet = jDst->GetNjet();
        Double_t et, eta, phi;

	for(Int_t ij=0; ij < nJet; ij++) {

	   jDst->GetJet(ij, 1, (*vJet)); // W
                
	   et  = vJet->Pt();
	   eta = vJet->Eta();
	   phi = vJet->Phi();
		
	   printf("\n Jet:%d E %f phi %f eta %f tracks %d\n", ij, 
		       et, phi, eta, 0);

		//		hNtracksInJet->Fill((Float_t)lJet->NTracks());
		// Safe zone - rJet - 13-feb-2002
	   if(nJet==1 && TMath::Abs(eta)<(0.7-rJet)&& 
               phi>rJet && phi<(phiMax-rJet)) { 
	      hJetEt->Fill(et);
	      hJetEta->Fill(eta);
	      hJetPhi->Fill(phi);
              if(inEmc) compareJetPartons(vJet, vPart1, vPart2);
           }

	} // jet
	cleanUpEvent();
   } // event

    gStyle->SetOptStat(1111111);

    //drawPartons();

#if !defined(__CINT__)
    if(hJetPartAng&&hJetPartAng->GetEntries()>2.){
      _drawListOfHist(lCompPJ, 2, 2);
      fitGauss(hDpt, 0); 
    }
#endif

    gROOT->cd();
    //    printf("\n Variant %1i : Alice Root file => %s \n",
    //mode, nf.Data());
}

Bool_t defineMicroDst(Int_t mode)
{
  TString dir("/auto/u/pavlinov/ALICE/RES/FILE/");
  //  TString  *sTmp = defineInput(mode); // see jetDst.h
  switch(mode){
    //  case -4: 
    //if(sTmp) fname = sTmp->Data();
    //else     return kFALSE;
    //break; 
  case -3: 
    fname = dir + "PythiaMode-3R0.5Seed4.0JEt20.0CEt0.0PtCut0.0.root"; 
    break; 
  case -2: 
    fname = dir + "PythiaMode-2R0.5Seed4.0JEt20.0CEt1.0PtCut1.0.root"; 
    break; 
  case -1: 
    fname = dir + "PythiaMode-1R0.5Seed4.0JEt20.0CEt1.0PtCut2.0.root"; 
    break; 

  case 11: 
    fname = dir + "HijingMode11R0.3Seed4.0JEt40.0CEt1.0PtCut2.0.root"; 
    break;
  default: 
     printf("\n<I> Wrong input mode => %i\n", mode);
     return kFALSE; 
  }
  jDst = new AliEMCALJetMicroDst;
  jDst->Open(fname.Data());
  gROOT->GetListOfBrowsables()->Add(jDst);
  return kTRUE;
}

Bool_t fillInfoForPartons(TVector3 *vec1, TVector3 *vec2) 
{
  // npart=4 for HIJING and 2 for PYTHIA 
  Int_t npart = jDst->GetNpart();
  if(!jDst->GetParton(npart-2, (*vec1))||!jDst->GetParton(npart-1, (*vec2)) ) 
  return kFALSE;

  static Double_t dphi;

  Bool_t inEmc = kFALSE;
  TVector3 *vec = vec1;
  for (Int_t i=0; i<2; i++) {
    Double_t ri = Double_t(i);
    hPartPt->Fill(vec->Pt(), ri);
    hPartEta->Fill(vec->Eta(), ri);
    hPartPhi->Fill(vec->Phi(), ri);
    if(fabs(vec->Eta())<0.7 && vec->Phi()>0 && vec->Phi()< phiMax) inEmc=kTRUE; 
    vec = vec2;
  }
  hPartDiffPt->Fill(vec1->Pt()-vec2->Pt());
  hPartDiffEta->Fill(vec1->Eta()-vec2->Eta());
  dphi = vec1->Phi() - vec2->Phi();
  dphi = fabs(TVector2::Phi_mpi_pi(dphi));
  hPartDiffPhi->Fill(dphi);
  return inEmc;
}

void compareJetPartons(TVector3 *vecJet, TVector3 *vec1, TVector3 *vec2)
{// 30-jan-2002
  TVector3 *v, *vw;
  Double_t ang, angw;

  v    = vec1;
  ang  = vecJet->Angle((*v));
  vw   = vec2;
  angw = vecJet->Angle((*vw));
  if(ang>angw) {
    ang = angw;
    v   = vw;
  }
  hJetPartAng->Fill(ang);
  hDeta->Fill(vecJet->Eta() - v->Eta());
  hDphi->Fill(vecJet->Phi() - v->Phi());
  hDpt->Fill(vecJet->Pt() - v->Pt());
}

void bookHist1()
{
    gROOT->cd();
// Book histos
    hNjet= new TH1F("hNjet","# jets",    11,  -0.5, 10.5);
    hJetEt   = new TH1F("hJetEt","Energy of jet",    150,  0.0, 450.);
    hJetEta  = new TH1F("hJetEta","#eta of jet",     180, -0.9,  0.9);
    hJetPhi = new TH1F("hJetPhi","#phi of jet",      62, -3.1,  3.1);
    hNtracksInJet   = new TH1F("hNtracksInJet","n tracks",   30,  0.5, 29.5);
    hPtTracksInJet  = new TH1F("hPtTracksInJet","p_{T} of tracks in jets cone", 
    100., 0., 100.);
    hTrackDR  = new TH1F("hTrackDR","Track dR", 120., 0.,   6.);    
#if !defined(__CINT__)
    lJets = _moveHistsToList("Hist for jets");
#endif
// QA
    hJetPartAng =  new TH1F("hJetPartAng","angle between jet and nearest parton", 100,0.0,0.2);
    hDeta =  new TH1F("hDeta","#eta (jets) - #eta (partons)", 100,  -.1, +.1);
    hDphi =  new TH1F("hDphi","#phi (jets} - #phi (partons)", 100,  -.1, +.1);
    hDpt  =  new TH1F("hDpt","E_{t} (jets) - E_{t} (partons)",100, -200., +200.);
#if !defined(__CINT__)
    lCompPJ = _moveHistsToList("Hist for jet vs parton");
#endif
}

void bookHistPartonsAfterHardSc()
{
  hPartPt  = new TH2F("hPartPt","Partons p_{t}", 500,  30., 1030., 2, -0.5, 1.5);
  hPartEta = new TH2F("hPartEta","Partons #eta", 100,  -5., 5., 2,-0.5,1.5);
  hPartPhi = new TH2F("hPartPhi","Partons #phi", 120,-0.2,2.2, 2,-0.5,1.5);
  hPartDiffPt  = new TH1F("hPartDiffPt","p_{t} difference for partons",100, -100., 100.);
  hPartDiffEta  = new TH1F("hPartDiffEta","#eta difference for partons",100, -5., 5.);
  hPartDiffPhi  = new TH1F("hPartDiffPhi","#phi angle for partons",320, 0.0, 3.2);
#if !defined(__CINT__)
  lPartons = _moveHistsToList("Hist for partons");
#endif
}

void drawPartons()
{
  gStyle->SetOptStat(1111111);
  c3 = new TCanvas("c3","Partons canvas",400,10,600,700);
  c3->Divide(2,3);
  TH1D *hid;
  
  c3->cd(1);
  hid = hPartPt->ProjectionX("_pt1",1,1); _drawHist(hid, 2);
  hid = hPartPt->ProjectionX("_pt2",2,2); _drawHist(hid, 1, 2, "same");
  c3->cd(2); _drawHist(hPartDiffPt,2);

  c3->cd(3);
  hid = hPartEta->ProjectionX("_eta1",1,1); _drawHist(hid, 2);
  hid = hPartEta->ProjectionX("_eta2",2,2); _drawHist(hid, 1, 2, "same");
  c3->cd(4); _drawHist(hPartDiffEta,2);

  c3->cd(5);
  hid = hPartPhi->ProjectionX("_phi1",1,1); _drawHist(hid, 2);
  hid = hPartPhi->ProjectionX("_phi2",2,2); _drawHist(hid, 1, 2, "same");
  c3->cd(6); _drawHist(hPartDiffPhi,2);
}

void fitGauss(TH1F* hid, Int_t opt)
{ // see /home/pavlinov/cosmic/pythia/geant/lin.C
  Double_t c = 2.5;
  Double_t xmax = hid->GetBinCenter(hid->GetMaximumBin());
  Double_t rms  = hid->GetRMS();
  Double_t xmi=xmax -c*rms, xma=xmax+c*rms;
  _fit(111);
  if(!g1) {
    g1 = new TF1("g1","gaus", xmi, xma);
  }
  g1->SetLineColor(2); // red
  g1->SetLineWidth(3);

  TString stmp = fname(0,fname.Length()-5); // discard .root
  if(!c1) {
    c1 = new TCanvas("c1","Jet resolution",400,10,600,700);
    p1 = _pad(stmp.Data(),c1);
  }
  p1->cd(); p1->Clear();
  hid->Fit("g1","R+");

  if(opt) {
    _drawHist(hid, 2);
    if(opt>=10) {
      stmp += "_res.ps";
      stmp.ReplaceAll("/FILE","");
      p1->Print(stmp.Data());
    }
  }
}

void cleanUpEvent()
{
  /* ???
  if(lPart[0]) {
    lPart[0]->Delete(); 
    lPart[0]=0;
    lPart[1]->Delete();
    lPart[1]=0;
  }
  */
}
