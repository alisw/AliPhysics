/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 
//---------------------------------------------------------------------
// JetAnalysis class 
// Analyse Jets (already found jets)
// Authors: andreas.morsch@cern.ch, jgcn@mail.cern.ch
//          mercedes.lopez.noriega@cern.ch
//---------------------------------------------------------------------
 
#include <Riostream.h>
#include "AliJetAnalysis.h"
ClassImp(AliJetAnalysis)
  
// root
#include <Riostream.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TMath.h>
// aliroot
#include "AliJetProductionDataPDC2004.h"
#include "AliJet.h"
#include "AliUA1JetHeader.h"
#include "AliLeading.h"
#include "AliJetReaderHeader.h"
  
AliJetAnalysis::AliJetAnalysis():
  fReaderHeader(0x0),
  fDirectory(0x0),
  fBkgdDirectory(0x0),
  fFile(0x0),
  fEventMin(0),
  fEventMax(-1),
  fRunMin(0),
  fRunMax(11),
  fminMult(0),
  fPercentage(-1.0),
  fPartPtCut(0.0),
  fdrJt(1.0),
  fdrdNdxi(0.7),
  fdrdEdr(1.0),
  fEfactor(1.0),
  fp0(0.0),
  fPtJ(0.0),
  fEJ(0.0),
  fEtaJ(0.0),
  fPhiJ(0.0),
  fjv3X(0.0),
  fjv3Y(0.0),
  fjv3Z(0.0),
  fPythia(kFALSE),
  fDoPart(kTRUE),
  fDoGenJ(kTRUE),
  fDoRecJ(kTRUE),
  fDoKine(kTRUE),
  fDoCorr(kTRUE),
  fDoCorr50(kFALSE),
  fDoShap(kTRUE),
  fDoFrag(kTRUE),
  fDoTrig(kTRUE),
  fDoJt(kTRUE),
  fDodNdxi(kTRUE),
  fDoBkgd(kTRUE),
  fWeight(1.0),
  fWShapR(0.0),
  fWFragR(0.0),
  fWdEdr(0.0),
  fWJt(0.0),
  fWdNdxi(0.0),
  fPart(0),
  fGenJ(0),
  fRecJ(0),
  fRecB(0),
  fRKineEneH(0),
  fRKinePtH(0),
  fRKinePhiH(0),
  fRKineEtaH(0),
  fGKineEneH(0),
  fGKinePtH(0),
  fGKinePhiH(0),
  fGKineEtaH(0),
  fPKineEneH(0),
  fPKinePtH(0),
  fPKinePhiH(0),
  fPKineEtaH(0),
  fPGCorrEneH(0),
  fPGCorrPtH(0),
  fPGCorrEtaH(0),
  fPGCorrPhiH(0),
  fPRCorrEneH(0),
  fPRCorrPtH(0),
  fPRCorrEtaH(0),
  fPRCorrPhiH(0),
  fRGCorrEneH(0),
  fRGCorrPtH(0),
  fRGCorrEtaH(0),
  fRGCorrPhiH(0),
  fPRCorr50EneH(0),
  fPRCorr50PtH(0),
  fPRCorr50EtaH(0),
  fPRCorr50PhiH(0),
  fRGCorr50EneH(0),
  fRGCorr50PtH(0),
  fRGCorr50EtaH(0),
  fRGCorr50PhiH(0),
  fRFragSelH(0),
  fRFragRejH(0),
  fRFragAllH(0),
  fRShapSelH(0),
  fRShapRejH(0),
  fRShapAllH(0),
  fGTriggerEneH(0),
  fRTriggerEneH(0),
  fGPTriggerEneH(0),
  fPTriggerEneH(0),
  fdEdrH(0),
  fdEdrB(0),
  fPtEneH2(0),
  fdEdrW(0),
  fJtH(0),
  fJtB(0),
  fJtW(0),
  fdNdxiH(0),
  fdNdxiB(0),
  fdNdxiW(0),
  fPtEneH(0)
{
  // Default constructor
  fFile          = "anaJets.root";   
  
  // initialize weight for dE/dr histo
  SetdEdrWeight();
}

AliJetAnalysis::~AliJetAnalysis()
{
  // Destructor
}


////////////////////////////////////////////////////////////////////////
// define histogrames 

void AliJetAnalysis::DefineHistograms()
{
  // Define the histograms to be filled
  if (fDoKine) DefineKineH();
  if (fDoCorr) DefineCorrH();
  if (fDoCorr50) DefineCorr50H();
  if (fDoShap) DefineShapH();
  if (fDoFrag) DefineFragH();
  if (fDoTrig) DefineTrigH();
  if (fDoJt) DefineJtH();
  if (fDodNdxi) DefinedNdxiH();
}

void AliJetAnalysis::DefineKineH()
{
  // leading particle    
  if (fDoPart) {
    fPKineEneH = new TH1F("PKineEne","Energy of leading particle",50,0,200);
    SetProperties(fPKineEneH,"Energy (GeV)","Entries");
    fPKinePtH = new TH1F("PKinePt","Pt of leading particle",50,0,200);
    SetProperties(fPKinePtH,"P_{T} (GeV)","Entries");
    fPKinePhiH = new TH1F("PKinePhiH","Azimuthal angle of leading particle",
			  90,0.,2.0*TMath::Pi());
    SetProperties(fPKinePhiH,"#phi","Entries");
    fPKineEtaH = new TH1F("PKineEtaH","Pseudorapidity of leading particle",
			  40,-1.0,1.0);
    SetProperties(fPKineEtaH,"#eta","Entries");
  }
  // leading generated jet
  if (fDoGenJ) {
    fGKineEneH = new TH1F("GKineEne","Energy of generated jet",50,0,200);
    SetProperties(fGKineEneH,"Energy (GeV)","Entries");
    fGKinePtH = new TH1F("GKinePt","Pt of generated jet",50,0,200);
    SetProperties(fGKinePtH,"P_{T} (GeV)","Entries");
    fGKinePhiH = new TH1F("GKinePhiH","Azimuthal angle of generated jet",
			  90,0.,2.0*TMath::Pi());
    SetProperties(fGKinePhiH,"#phi","Entries");
    fGKineEtaH = new TH1F("GKineEtaH","Pseudorapidity of generated jet",
			  40,-1.0,1.0);
    SetProperties(fGKineEtaH,"#eta","Entries");
  }
  // leading reconstructed jet
  if (fDoRecJ) {
    fRKineEneH = new TH1F("RKineEne","Energy of reconstructed jet",50,0,200);
    SetProperties(fRKineEneH,"Energy (GeV)","Entries");
    fRKinePtH = new TH1F("RKinePt","Pt of reconstructed jet",50,0,200);
    SetProperties(fRKinePtH,"P_{T} (GeV)","Entries");
    fRKinePhiH = new TH1F("RKinePhiH","Azimuthal angle of reconstructed jet",
			  90,0.,2.0*TMath::Pi());
    SetProperties(fRKinePhiH,"#phi","Entries");
    fRKineEtaH = new TH1F("RKineEtaH","Pseudorapidity of reconstructed jet",
			  40,-1.0,1.0);
    SetProperties(fRKineEtaH,"#eta","Entries");
  } 
}

void AliJetAnalysis::DefineCorrH()
{
  // correlation 
  if (fDoPart && fDoGenJ) {
    fPGCorrEneH = new TH2F("PGCorrEne","Energy correlation part-gen jet",
			   40,0,200,40,0,200);
    SetProperties(fPGCorrEneH,"Part Energy (GeV)","Gen Jet Energy (GeV)");
    fPGCorrPtH = new TH2F("PGCorrPt","Pt correlation part-gen jet",
			  40,0,200,40,0,200);
    SetProperties(fPGCorrPtH,"Part P_{T} (GeV)","Gen Jet P_{T} (GeV)");
    fPGCorrEtaH = new TH2F("PGCorrEta","Pseudorapidity correlation part-gen jet",
			   40,-1.0,1.0,40,-1.0,1.0);
    SetProperties(fPGCorrEtaH,"Part #eta","Gen Jet #eta");
    fPGCorrPhiH = new TH2F("PGCorrPhi","Azimuthal angle correlation part-gen jet",
			   90,0.,2.0*TMath::Pi(),90,0.,2.0*TMath::Pi());
    SetProperties(fPGCorrPhiH,"Part #phi","Gen Jet #phi");
  }  
  if (fDoPart && fDoRecJ) {
    fPRCorrEneH = new TH2F("PRCorrEne","Energy correlation part-rec jet",
			   40,0,200,40,0,200);
    SetProperties(fPRCorrEneH,"Part Jet Energy (GeV)","Rec Jet Energy (GeV)");
    fPRCorrPtH = new TH2F("PRCorrPt","Pt correlation part-rec jet",
			  40,0,200,40,0,200);
    SetProperties(fPRCorrPtH,"Part Jet P_{T} (GeV)","Rec Jet P_{T} (GeV)");
    fPRCorrEtaH = new TH2F("PRCorrEta","Pseudorapidity correlation part-rec jet",
			   40,-1.0,1.0,40,-1.0,1.0);
    SetProperties(fPRCorrEtaH,"part #eta","Rec Jet #eta");
    fPRCorrPhiH = new TH2F("PRCorrPhi","Azimuthal angle correlation part-rec jet",
			   90,0.,2.0*TMath::Pi(),90,0.,2.0*TMath::Pi());
    SetProperties(fPRCorrPhiH,"Part #phi","Rec Jet #phi");
  }  
  if (fDoGenJ && fDoRecJ) {
    fRGCorrEneH = new TH2F("RGCorrEne","Energy correlation rec jet-gen jet",
			   40,0,200,40,0,200);
    SetProperties(fRGCorrEneH,"Rec Jet Energy (GeV)","Gen Jet Energy (GeV)");
    fRGCorrPtH = new TH2F("RGCorrPt","Pt correlation rec jet-gen jet",
			  40,0,200,40,0,200);
    SetProperties(fRGCorrPtH,"Rec Jet P_{T} (GeV)","Gen Jet P_{T} (GeV)");
    fRGCorrEtaH = new TH2F("RGCorrEta","Pseudorapidity correlation rec jet-gen jet",
			   40,-1.0,1.0,40,-1.0,1.0);
    SetProperties(fRGCorrEtaH,"Rec Jet #eta","Gen Jet #eta");
    fRGCorrPhiH = new TH2F("RGCorrPhi","Azimuthal angle correlation rec jet-gen jet",
			   90,0.,2.0*TMath::Pi(),90,0.,2.0*TMath::Pi());
    SetProperties(fRGCorrPhiH,"Rec Jet #phi","Gen Jet #phi");
  }
}

void AliJetAnalysis::DefineCorr50H()
{
  // correlation 
  if (fDoPart && fDoRecJ) {
    fPRCorr50EneH = new TH2F("PRCorr50Ene","Energy correlation part-rec jet",
			     40,0,200,40,0,200);
    SetProperties(fPRCorr50EneH,"Part Jet Energy (GeV)","Rec Jet Energy (GeV)");
    fPRCorr50PtH = new TH2F("PRCorr50Pt","Pt correlation part-rec jet",
			    40,0,200,40,0,200);
    SetProperties(fPRCorr50PtH,"Part Jet P_{T} (GeV)","Rec Jet P_{T} (GeV)");
    fPRCorr50EtaH = new TH2F("PRCorr50Eta","Pseudorapidity correlation part-rec jet",
			     40,-1.0,1.0,40,-1.0,1.0);
    SetProperties(fPRCorr50EtaH,"part #eta","Rec Jet #eta");
    fPRCorr50PhiH = new TH2F("PRCorr50Phi","Azimuthal angle correlation part-rec jet",
			     90,0.,2.0*TMath::Pi(),90,0.,2.0*TMath::Pi());
    SetProperties(fPRCorr50PhiH,"Part #phi","Rec Jet #phi");
  }
  
  if (fDoGenJ && fDoRecJ) {
    fRGCorr50EneH = new TH2F("RGCorr50Ene","Energy correlation rec jet-gen jet",
			     40,0,200,40,0,200);
    SetProperties(fRGCorr50EneH,"Rec Jet Energy (GeV)","Gen Jet Energy (GeV)");
    fRGCorr50PtH = new TH2F("RGCorr50Pt","Pt correlation rec jet-gen jet",
			    40,0,200,40,0,200);
    SetProperties(fRGCorr50PtH,"Rec Jet P_{T} (GeV)","Gen Jet P_{T} (GeV)");
    fRGCorr50EtaH = new TH2F("RGCorr50Eta","Pseudorapidity correlation rec jet-gen jet",
			     40,-1.0,1.0,40,-1.0,1.0);
    SetProperties(fRGCorr50EtaH,"Rec Jet #eta","Gen Jet #eta");
    fRGCorr50PhiH = new TH2F("RGCorr50Phi","Azimuthal angle correlation rec jet-gen jet",
			     90,0.,2.0*TMath::Pi(),90,0.,2.0*TMath::Pi());
    SetProperties(fRGCorr50PhiH,"Rec Jet #phi","Gen Jet #phi");
  }
}

void AliJetAnalysis::DefineShapH()
{
  // leading reconstructed jet
  if (fDoRecJ) {
    fdEdrH = new TH2F("fdEdrH","dE/dr histo",20,0,1,40,0,200);
    SetProperties(fdEdrH,"r","Rec Jet P_{T}");
    fdEdrB = new TH2F("fdEdrB","dE/dr bkgdhisto",20,0,1,40,0,200);
    SetProperties(fdEdrB,"r","Rec P_{T}");
    fPtEneH2 = new TH2F("fPtEneH2","fPtEneH2",40,0,200,40,0,200);
    SetProperties(fPtEneH2,"Rec Jet E","Rec Jet P_{T}");
    fdEdrW = new TH1F("fdEdrW","weights for dE/dr",40,0,200);
    SetProperties(fdEdrW,"Rec Jet P_{T}","weights");
    
    fRShapSelH = new TH1F("RShapSel","Shape of generated jets (sel part)",20,0.,1.);
    SetProperties(fRShapSelH,"r","1/N_{JET}#Sigma P_{T}(0,r)/P_{T}(0,R_{JET})");
    fRShapRejH = new TH1F("RShapRej","Shape of generated jets (rej part)",20,0.,1.);
    SetProperties(fRShapRejH,"r","1/N_{JET}#Sigma P_{T}(0,r)/P_{T}(0,R_{JET})");
    fRShapAllH = new TH1F("RShapAll","Shape of generated jets (all part)",20,0.,1.);
    SetProperties(fRShapAllH,"r","1/N_{JET}#Sigma P_{T}(0,r)/P_{T}(0,R_{JET})");
  }
}

void AliJetAnalysis::DefineJtH()
{
  // Define the histogram for J_T
  if (fDoRecJ) {
    fJtH = new TH2F("fJtH","J_{T} histogram",80,0.,4.,40,0.,200.);
    SetProperties(fJtH,"J_{T}","Rec Jet P_{T}");
    fJtB = new TH2F("fJtB","J_{T} bkgd histogram",80,0.,4.,40,0.,200.);
    SetProperties(fJtB,"J_{T}","Rec P_{T}");
    fJtW = new TH1F("fJtW","J_{T} weight",40,0,200);
    SetProperties(fJtW,"J_{T}W","weight");
  }
}

void AliJetAnalysis::DefinedNdxiH()
{
  // Define the histogram for dN/dxi
  if (fDoRecJ) {
    fdNdxiH = new TH2F("fdNdxiH","dN/d#xi histo",200,0,10,40,0,200);
    SetProperties(fdNdxiH,"xi","Rec Jet P_{T}");
    fdNdxiB = new TH2F("fdNdxiB","dN/d#xi bkgd histo",200,0,10,40,0,200);
    SetProperties(fdNdxiB,"xi","Rec P_{T}");
    fdNdxiW = new TH1F("fdNdxiW","dN/d#xi histo",40,0,200);
    SetProperties(fdNdxiW,"Rec Jet P_{T}","weights");
    fPtEneH = new TH2F("fPtEneH","fPtEneH",40,0,200,40,0,200);
    SetProperties(fPtEneH,"Rec Jet E","Rec Jet P_{T}");
  }
}

void AliJetAnalysis::DefineFragH()
{
  // leading reconstructed jet
  if (fDoRecJ) {
    fRFragSelH = new TH1F("RFragSel","Frag Fun of reconstructed jets (sel part)",20,0.,10.);
    SetProperties(fRFragSelH,"#xi=ln(E_{JET}/E_{i})","1/N_{JET}dN_{ch}/d#xi");
    fRFragRejH = new TH1F("RFragRej","Frag Fun of reconstructed jets (rej part)",20,0.,10.);
    SetProperties(fRFragRejH,"#xi=ln(E_{JET}/E_{i})","1/N_{JET}dN_{ch}/d#xi");
    fRFragAllH = new TH1F("RFragAll","Frag Fun of reconstructed jets (all part)",20,0.,10.);
    SetProperties(fRFragAllH,"#xi=ln(E_{JET}/E_{i})","1/N_{JET}dN_{ch}/d#xi");
  }
}

void AliJetAnalysis::DefineTrigH()
{
  // generated energy
  fGTriggerEneH = new TProfile("GTriggerEne","Generated energy (trigger bias)", 
			       20, 0., 200., 0., 1., "S");    
  fGTriggerEneH->SetXTitle("E_{gen}");
  fGTriggerEneH->SetYTitle("E_{rec}/E_{gen}");
  // reconstructed energy
  fRTriggerEneH = new TProfile("RTriggerEne","Reconstructed energy (trigger bias)", 
			       20, 0., 200., 0., 1., "S");  
  fRTriggerEneH->SetXTitle("E_{rec}");
  fRTriggerEneH->SetYTitle("E_{rec}/E_{gen}");
  // generated energy vs generated/leading
  fGPTriggerEneH = new TProfile("GPTriggerEne","Generated energy (trigger bias)", 
				20, 0., 200., 0., 1., "S");    
  fGPTriggerEneH->SetXTitle("E_{gen}");
  fGPTriggerEneH->SetYTitle("E_{L}/E_{gen}");
  // leading particle energy
  fPTriggerEneH = new TProfile("PTriggerEne","Leading particle energy (trigger bias)", 
			       20, 0., 200., 0., 1., "S");  
  fPTriggerEneH->SetXTitle("E_{L}/0.2");
  fPTriggerEneH->SetYTitle("E_{L}/E_{gen}");

}

void AliJetAnalysis::SetProperties(TH1* h,const char* x, const char* y) const
{
  //Set properties of histos (x and y title and error propagation)
  h->SetXTitle(x);
  h->SetYTitle(y);
  h->Sumw2();
}
 
////////////////////////////////////////////////////////////////////////
// compute weights for dE/dr histogram

void AliJetAnalysis::SetdEdrWeight()
{
  // Due to the limited acceptance, each bin in the dE/dr has a different
  // weight given by the ratio of the area of a disk dR to the area
  // within the eta acceptance. The weight depends on the eta position 
  // of the jet. Here a look up table for the weights is computed. It
  // assumes |etaJet|<0.5 and |eta_lego|<0.9. It makes bins of 0.05
  // units in eta. Note that this is fixed by the bin width chosen for
  // the histogram --> this weights are tailored for the specific 
  // histogram definition used here!
  
  // two dimensional table: first index, bin in eta of jet, second
  // index bin of dE/dr histo

  Int_t nBins = 20;
  Float_t xMin = 0.0;
  Float_t xMax = 1.0;
  Float_t binSize = (xMax-xMin)/nBins;

  Float_t da,ds,r1,r2,h1,h2,a1,a2,theta1,theta2,area1,area2;
  Int_t ji;
  for (Int_t i=0;i<(nBins/2);i++) {
    // index of first histo bin needing a scale factor
    ji=(nBins-2)-i;
    for(Int_t j=0;j<nBins;j++) {
      // area of ring.
      da = TMath::Pi()*(binSize*binSize)*(1.0+2.0*j);
      ds = 1.0;
      if (j>=ji) { // compute missing area using segments of circle
	r1=j*binSize;
	r2=(j+1)*binSize;
	h1=(j-ji)*binSize;
	h2=(j+1-ji)*binSize;
	a1=2.0*TMath::Sqrt(2.0*h1*r1-h1*h1);
	a2=2.0*TMath::Sqrt(2.0*h2*r2-h2*h2);
	theta1=2*TMath::ACos((r1-h1)/r1);
	theta2=2*TMath::ACos((r2-h2)/r2);
	area1=binSize*(r1*r1*theta1-a1*(r1-h1));
	area2=binSize*(r2*r2*theta2-a2*(r2-h2));
	ds = (da-(area2-area1))/da;
      }
      fWeightdEdr[i][j]=ds/da;
    }
  }
}

// get weight for dE/dr histogram
Float_t AliJetAnalysis::GetdEdrWeight(Float_t etaJet, Float_t r)
{
  // Return the correponding weight for the dE/dr plot
  Int_t nBins = 20;
  Float_t xMin = 0.0;
  Float_t xMax = 1.0;
  Float_t binSize = (xMax-xMin)/nBins;

  Float_t eta = TMath::Abs(etaJet);
  if ((eta > 0.5) || (r > fdrdEdr)) return 0.0;
  Int_t binJet = (Int_t) TMath::Floor(eta/binSize);
  Int_t binR = (Int_t) TMath::Floor(r/binSize);
  Float_t w = fWeightdEdr[binJet][binR];
  return w;
}


////////////////////////////////////////////////////////////////////////
// fill histograms 

void AliJetAnalysis::FillHistograms()
{
  // fill histograms 

  // Run data 
  AliJetProductionDataPDC2004* runData = new AliJetProductionDataPDC2004();
  
  // Loop over runs
  TFile* jFile = 0x0;
  
  for (Int_t iRun = fRunMin; iRun <= fRunMax; iRun++) {
    // Open file
    char fn[256];
    sprintf(fn, "%s/%s.root", fDirectory, (runData->GetRunTitle(iRun)).Data());
    jFile = new TFile(fn);
    printf("  Analyzing run: %d %s\n", iRun,fn);	
    
    // Get reader header and events to be looped over
    AliJetReaderHeader *jReaderH = (AliJetReaderHeader*)(jFile->Get(fReaderHeader));
    if (fEventMin == -1) fEventMin =  jReaderH->GetFirstEvent();
    if (fEventMax == -1) {
      fEventMax =  jReaderH->GetLastEvent();
    } else {
      fEventMax = TMath::Min(fEventMax, jReaderH->GetLastEvent());
    }
    
    AliUA1JetHeader *jH = (AliUA1JetHeader *) (jFile->Get("AliUA1JetHeader"));
    
    // Calculate weight
    fWeight = runData->GetWeight(iRun)/ ( (Float_t) (fEventMax - fEventMin + 1));

    // Loop over events
    for (Int_t i = fEventMin; i < fEventMax; i++) {
      if (i%100 == 0) printf("  Analyzing event %d / %d \n",i,fEventMax);
      
      // Get next tree with AliJet
      char nameT[100];
      sprintf(nameT, "TreeJ%d",i);
      TTree *jetT =(TTree *)(jFile->Get(nameT));
      if (fDoRecJ) jetT->SetBranchAddress("FoundJet",    &fRecJ);
      if (fDoGenJ) jetT->SetBranchAddress("GenJet",      &fGenJ);
      if (fDoPart) jetT->SetBranchAddress("LeadingPart", &fPart);
      jetT->GetEntry(0);
      
      int nJets = fRecJ->GetNJets();  
      
      TArrayI inJet = fRecJ->GetInJet();
      if(inJet.GetSize()>fminMult){     // removes events with low multiplicity
	if (fDoKine) FillKineH();
	if (fDoCorr) FillCorrH();
	if (fDoCorr50) FillCorr50H();
	if (fDoShap) FillShapH(jH->GetRadius());
	if (fDoFrag) FillFragH();    
	if (fDoTrig) FillTrigH();
	if (fDoJt) FillJtH();
	if (fDodNdxi) FilldNdxiH();
      }
      delete jetT;                       // jet should be deleted before creating a new one
      if(inJet.GetSize()>fminMult){      // removes events with low multiplicity
	if (fDoBkgd && nJets>0) FillBkgd(i,iRun);
      }      
    } // end loop over events in one file
    if (jFile) jFile->Close();
    delete jFile;
  } // end loop over files
}

void AliJetAnalysis::FillKineH()
{
  // leading particle 
  if (fDoPart && fPart->LeadingFound()){
    fPKineEneH->Fill(fPart->GetE(),fWeight);
    fPKinePtH->Fill(fPart->GetPt(),fWeight);
    fPKinePhiH->Fill(fPart->GetPhi(),fWeight);
    fPKineEtaH->Fill(fPart->GetEta(),fWeight);
  } 
  // leading generated jet
  if (fDoGenJ && fGenJ->GetNJets()> 0){
    fGKineEneH->Fill(fGenJ->GetE(0),fWeight);
    fGKinePtH->Fill(fGenJ->GetPt(0),fWeight);
    fGKinePhiH->Fill(fGenJ->GetPhi(0),fWeight);
    fGKineEtaH->Fill(fGenJ->GetEta(0),fWeight);
  } 
  // leading reconstructed jet
  if (fDoRecJ && fRecJ->GetNJets()> 0) {
    TArrayF p=fRecJ->GetPtFromSignal();
    if (p[0]>fPercentage) {
      fRKineEneH->Fill(fRecJ->GetE(0)/fEfactor,fWeight);
      fRKinePtH->Fill(fRecJ->GetPt(0),fWeight);
      fRKinePhiH->Fill(fRecJ->GetPhi(0),fWeight);
      fRKineEtaH->Fill(fRecJ->GetEta(0),fWeight);
    }
  }
}

void AliJetAnalysis::FillCorrH()
{
  // Fill correlation histograms
  if (fDoPart && fPart->LeadingFound() && fDoGenJ && fGenJ->GetNJets()> 0) 
    Correlation(fPart->GetLeading(),fGenJ->GetJet(0), 
		fPGCorrEneH,fPGCorrPtH,fPGCorrEtaH,fPGCorrPhiH);
  if (fDoPart && fPart->LeadingFound() && fDoRecJ && fRecJ->GetNJets()> 0) {
    TArrayF p=fRecJ->GetPtFromSignal();
    if (p[0]>fPercentage) 
      Correlation(fPart->GetLeading(),fRecJ->GetJet(0), 
		  fPRCorrEneH,fPRCorrPtH,fPRCorrEtaH,fPRCorrPhiH);
  }
  if (fDoGenJ && fGenJ->GetNJets()> 0  && fDoRecJ && fRecJ->GetNJets()> 0) {
    TArrayF p=fRecJ->GetPtFromSignal();
    if (p[0]>fPercentage) 
      Correlation(fRecJ->GetJet(0), fGenJ->GetJet(0),
		  fRGCorrEneH,fRGCorrPtH,fRGCorrEtaH,fRGCorrPhiH);
  }
}

void AliJetAnalysis::Correlation(TLorentzVector *lv1,TLorentzVector *lv2,
				 TH2F *h1, TH2F *h2, TH2F *h3, TH2F *h4)
{
  // Correlation histograms
  h1->Fill(lv1->E(),lv2->E(),fWeight);
  h2->Fill(lv1->Pt(),lv2->Pt(),fWeight);
  h3->Fill(lv1->Eta(),lv2->Eta(),fWeight);
  Float_t p1, p2;
  p1 = ((lv1->Phi() < 0) ? (lv1->Phi()) + 2. * TMath::Pi() : lv1->Phi());
  p2 = ((lv2->Phi() < 0) ? (lv2->Phi()) + 2. * TMath::Pi() : lv2->Phi());
  h4->Fill(p1,p2,fWeight);
}

void AliJetAnalysis::FillCorr50H()
{
  // Fill correlation histograms when one particle has > 50% of jet energy
  if (fDoRecJ && fRecJ->GetNJets()> 0) {
    TArrayF p = fRecJ->GetPtFromSignal();
    if (p[0]>fPercentage) {
      if (fDoPart && fPart->LeadingFound()) 
	Correlation50(fRecJ, fPart->GetLeading(),fRecJ->GetJet(0), 
		      fPRCorr50EneH,fPRCorr50PtH,fPRCorr50EtaH,fPRCorr50PhiH);
      if (fDoGenJ && fGenJ->GetNJets()> 0) 
	Correlation50(fRecJ, fRecJ->GetJet(0), fGenJ->GetJet(0),
		      fRGCorr50EneH,fRGCorr50PtH,fRGCorr50EtaH,fRGCorr50PhiH);
    }
  }
}

void AliJetAnalysis::Correlation50(AliJet *j,TLorentzVector *lv1,TLorentzVector *lv2,
				   TH2F *h1, TH2F *h2, TH2F *h3, TH2F *h4)
{
  // Correlation histograms when one particle has > 50% of jet energy
    TArrayF ptin = j->GetPtIn();
    TArrayF etain = j->GetEtaIn();
    TArrayI inJet = j->GetInJet();
    
    Int_t flag = 0;
    for(Int_t i=0;i<(inJet.GetSize());i++) {
      if (inJet[i] == 1) {
	Float_t t1 = TMath::Tan(2.0*TMath::ATan(TMath::Exp(-etain[i])));
	Float_t x = (ptin[i]*TMath::Sqrt(1.+1./(t1*t1)))/(j->GetE(0));
	if (x>0.5) flag++;
      }
    }
    if (flag>1) cout << " Flag = " << flag << endl;
    if (flag == 1) {
      h1->Fill(lv1->E(),lv2->E(),fWeight);
      h2->Fill(lv1->Pt(),lv2->Pt(),fWeight);
      h3->Fill(lv1->Eta(),lv2->Eta(),fWeight);
      Float_t p1, p2;
      p1 = ((lv1->Phi() < 0) ? 
	    (lv1->Phi()) + 2. * TMath::Pi() : lv1->Phi());
      p2 = ((lv2->Phi() < 0) ? 
	    (lv2->Phi()) + 2. * TMath::Pi() : lv2->Phi());
      h4->Fill(p1,p2,fWeight);
    }
}

void AliJetAnalysis::FillJtH()
{
  // Fill J_T histogram
  if (fRecJ->GetNJets()> 0) {
    fjv3X = 0.0; fjv3Y = 0.0; fjv3Z = 0.0;
    TArrayF p=fRecJ->GetPtFromSignal();
    if (p[0]>fPercentage) {
      // initialize
      const TVector3 kJv3 = fRecJ->GetJet(0)->Vect();
      fjv3X =  kJv3.X(); fjv3Y =  kJv3.Y(); fjv3Z =  kJv3.Z();
      TVector3 trk;
      TArrayI inJet = fRecJ->GetInJet();
      TArrayF etain = fRecJ->GetEtaIn();
      TArrayF ptin = fRecJ->GetPtIn();
      TArrayF phiin = fRecJ->GetPhiIn();
      Float_t deta, dphi,jt, dr;
      for(Int_t i=0;i<inJet.GetSize();i++) {
	deta = etain[i] - fRecJ->GetEta(0);
	dphi = phiin[i] - fRecJ->GetPhi(0);
	if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	dr = TMath::Sqrt(deta * deta + dphi * dphi);
	if ((dr<fdrJt) && (ptin[i] > fPartPtCut)) { 
	  trk.SetPtEtaPhi(ptin[i],etain[i],phiin[i]);
	  jt = TMath::Sin(trk.Angle(kJv3))*trk.Mag();
	  fJtH->Fill(jt,fRecJ->GetPt(0),fWeight);
	}
      }
      fJtW->Fill(fRecJ->GetPt(0),fWeight);
      fWJt+=fWeight;
    }
  }
}

void AliJetAnalysis::FilldNdxiH()
{
  // Fill dN/d#xi histograms
  if (fRecJ->GetNJets()> 0) {
    TArrayF p=fRecJ->GetPtFromSignal();
    if (p[0]>fPercentage) {
      fp0 = p[0];
      TArrayI inJet = fRecJ->GetInJet();
      TArrayF etain = fRecJ->GetEtaIn();
      TArrayF ptin = fRecJ->GetPtIn();
      TArrayF phiin = fRecJ->GetPhiIn();
      Float_t xi,t1,ene,dr,deta,dphi;
      for(Int_t i=0;i<inJet.GetSize();i++) {
	deta = etain[i] - fRecJ->GetEta(0);
	dphi = phiin[i] - fRecJ->GetPhi(0);
	if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
	if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
	dr = TMath::Sqrt(deta * deta + dphi * dphi);
	if ((dr<fdrdNdxi) && (ptin[i] > fPartPtCut)) { 
	  t1 = TMath::Tan(2.0*TMath::ATan(TMath::Exp(-etain[i])));
	  ene = ptin[i]*TMath::Sqrt(1.+1./(t1*t1));
	  xi = (Float_t) TMath::Log((fRecJ->GetE(0)/fEfactor)/ene);
	  fdNdxiH->Fill(xi,fRecJ->GetPt(0),fWeight);
	}  
      } 
      fdNdxiW->Fill(fRecJ->GetPt(0),fWeight);
      fPtEneH->Fill(fRecJ->GetE(0)/fEfactor,fRecJ->GetPt(0),fWeight);
      fWdNdxi+=fWeight;
    }
  }
}

void AliJetAnalysis::FillBkgd(Int_t eventN, Int_t runN)
{
  // Background calculated with hijing events (no pythia)  
  if (fp0>fPercentage) {
    fPtJ=0.,fEJ=0.,fEtaJ=0.,fPhiJ=0.;
    // Background calculated with hijing events (no pythia)
    AliJetProductionDataPDC2004* runDataB = new AliJetProductionDataPDC2004();
    TFile* jFileB =0x0;;
    char fnB[256];
    
    sprintf(fnB, "%s/%s.root", fBkgdDirectory, (runDataB->GetRunTitle(runN)).Data());
    jFileB = new TFile(fnB);

    char nameB[100];
    sprintf(nameB, "TreeJ%d",eventN);
    TTree *jetB =(TTree *)(jFileB->Get(nameB));
    jetB->SetBranchAddress("FoundJet",    &fRecB);
    jetB->GetEntry(0);
    
    TArrayI inJetB = fRecB->GetInJet();
    TArrayF etainB = fRecB->GetEtaIn();
    TArrayF ptinB = fRecB->GetPtIn();
    TArrayF phiinB = fRecB->GetPhiIn();
    fPtJ = fRecJ->GetPt(0);
    fEJ = fRecJ->GetE(0);
    fEtaJ = fRecJ->GetEta(0);
    fPhiJ = fRecJ->GetPhi(0);
    Float_t t1,ene,xi,detaB,dphiB,drB,jt,wB;
    TVector3 trkB;
    TVector3 jv3B;
    jv3B.SetX(fjv3X); jv3B.SetY(fjv3Y); jv3B.SetZ(fjv3Z);
    
    for(Int_t k=0;k<inJetB.GetSize();k++){
      if(ptinB[k] > fPartPtCut){
	detaB = etainB[k] - fEtaJ;
	dphiB = phiinB[k] - fPhiJ;
	if (dphiB < -TMath::Pi()) dphiB= -dphiB - 2.0 * TMath::Pi();
	if (dphiB > TMath::Pi()) dphiB = 2.0 * TMath::Pi() - dphiB;
	drB = TMath::Sqrt(detaB * detaB + dphiB * dphiB);
	t1 = TMath::Tan(2.0*TMath::ATan(TMath::Exp(-etainB[k])));
	ene = ptinB[k]*TMath::Sqrt(1.+1./(t1*t1));
	trkB.SetPtEtaPhi(ptinB[k],etainB[k],phiinB[k]);
	// --- dN/dxi
	if (drB<fdrdNdxi) {
	  xi = (Float_t) TMath::Log(fEJ/ene);
	  fdNdxiB->Fill(xi,fPtJ,fWeight);
	}
	// --- Jt
	if (drB<fdrJt) {
	  jt = TMath::Sin(trkB.Angle(jv3B))*(trkB.Mag());
	  fJtB->Fill(jt,fPtJ,fWeight);
	}
	// --- dE/dr
	if (drB<fdrdEdr) {
	  wB = GetdEdrWeight(fEtaJ,drB)*fWeight*ene;
	  fdEdrB->Fill(drB,fPtJ,wB);
	} 
      }
    }
    delete jetB;
    if (jFileB) jFileB->Close();
    delete jFileB; 
  }
}

void AliJetAnalysis::FillShapH(Float_t r)
{
  // Fill jet shape histograms
  if (fDoRecJ && fRecJ->GetNJets()> 0) {
    TArrayF p=fRecJ->GetPtFromSignal();
    if (p[0]>fPercentage) {
      Shape(fRecJ,fRShapSelH,fRShapRejH,fRShapAllH,fdEdrH,fPtEneH2,fdEdrW,r);
      fWdEdr+=fWeight;
      fWShapR+=fWeight;
    }
  }
}

void AliJetAnalysis::Shape(AliJet *j,TH1F* hs, TH1F* hr, TH1F* ha, 
			   TH2F *hdedr, TH2F *hptene, TH1F* wdedr, Float_t radius)
{
  // initialize
  TArrayI inJet = j->GetInJet();
  TArrayF etain = j->GetEtaIn();
  TArrayF ptin = j->GetPtIn();
  TArrayF phiin = j->GetPhiIn();
  
  // first compute dE/dr histo
  Float_t etaj = j->GetEta(0);
  Float_t ene,w,deta,dphi,dr,t1;
  for(Int_t i=0;i<inJet.GetSize();i++) {
    deta = etain[i] - j->GetEta(0);
    dphi = phiin[i] - j->GetPhi(0);
    if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
    if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
    dr = TMath::Sqrt(deta * deta + dphi * dphi);
    if ((dr<fdrdEdr) && (ptin[i] > fPartPtCut)) {
      t1 = TMath::Tan(2.0*TMath::ATan(TMath::Exp(-etain[i])));
      ene = ptin[i]*TMath::Sqrt(1.+1./(t1*t1));
      w = GetdEdrWeight(etaj,dr)*fWeight*ene;
      hdedr->Fill(dr,j->GetPt(0),w);
    }
  }
  hptene->Fill(fRecJ->GetE(0),fRecJ->GetPt(0),fWeight);
  wdedr->Fill(j->GetPt(0),fWeight);
  
  // now compute shape histos
  Int_t nBins = ha->GetNbinsX();
  Float_t xMin = ha->GetXaxis()->GetXmin();
  Float_t xMax = ha->GetXaxis()->GetXmax();
  Float_t binSize,halfBin;
  binSize = (xMax-xMin)/nBins;
  halfBin = binSize/2;
  Float_t rptS[20], rptR[20], rptA[20];
  for(Int_t i=0;i<nBins;i++) rptS[i]=rptR[i]=rptA[i]=0.0;
  // fill bins in r for leading jet
  for(Int_t i=0;i<inJet.GetSize();i++) {
    if (inJet[i] == 1 || (inJet[i] == -1 && fPythia)) {
      deta = etain[i] - j->GetEta(0);
      dphi = phiin[i] - j->GetPhi(0);
      if (dphi < -TMath::Pi()) dphi= -dphi - 2.0 * TMath::Pi();
      if (dphi > TMath::Pi()) dphi = 2.0 * TMath::Pi() - dphi;
      dr = TMath::Sqrt(deta * deta + dphi * dphi);
      Float_t rR = dr/radius;
      Int_t bin = (Int_t) TMath::Floor(rR/binSize);
      rptA[bin]+=ptin[i]/(j->GetPt(0));
      if (inJet[i] == 1) rptS[bin]+=ptin[i]/(j->GetPt(0));
      if (fPythia && inJet[i] == -1) 
	rptR[bin]+=ptin[i]/(j->GetPt(0));
    }
  }
  
  // compute shape and fill histogram
  Float_t ptS,ptR,ptA,r;
  ptS=ptR=ptA=0.0;
  for (Int_t i=0;i<nBins;i++) {
    ptS+=rptS[i];
    ptR+=rptR[i];
    ptA+=rptA[i];
    r=(i+1)*binSize-halfBin;
    hs->Fill(r,ptS*fWeight);
    if(fPythia) {
      hr->Fill(r,ptR*fWeight);
      ha->Fill(r,ptA*fWeight);
    }
  }
}

void AliJetAnalysis::FillFragH()
{
  // Fill fragmentation histogram
  if (fDoRecJ && fRecJ->GetNJets()> 0) {
    TArrayF p=fRecJ->GetPtFromSignal();
    if (p[0]>fPercentage) {
      FragFun(fRecJ,fRFragSelH,fRFragRejH,fRFragAllH);
      fWFragR+=fWeight;
    }
  }
}

void AliJetAnalysis::FragFun(AliJet *j,TH1F* hs, TH1F* hr, TH1F* ha)
{
  // Calculate de fragmentation function
  TArrayI inJet = j->GetInJet();
  TArrayF etain = j->GetEtaIn();
  TArrayF ptin = j->GetPtIn();
  
  Float_t xi,t1,ene;
  
  for(Int_t i=0;i<(inJet.GetSize());i++) {
    if (inJet[i] == 1 || (inJet[i] == -1 && fPythia)) {
      t1 = TMath::Tan(2.0*TMath::ATan(TMath::Exp(-etain[i])));
      ene = ptin[i]*TMath::Sqrt(1.+1./(t1*t1));
      xi = (Float_t) TMath::Log((j->GetE(0))/ene);
      if (fPythia) ha->Fill(xi,fWeight);
      if (inJet[i] == 1) hs->Fill(xi,fWeight);
      if (fPythia && inJet[i] == -1) hr->Fill(xi,fWeight);
    }
  }
}

void AliJetAnalysis:: FillTrigH()
{
  // Fill trigger bias histograms
  if(fDoTrig && fGenJ->GetNJets()>0 && fRecJ->GetNJets()>0){	
    TArrayF p=fRecJ->GetPtFromSignal();
    if (p[0]>fPercentage) {
      float genEne = fGenJ->GetE(0);
      float recEne = fRecJ->GetE(0);
      float eneRatio = (recEne/fEfactor)/genEne;
      
      fGTriggerEneH->Fill(genEne, eneRatio, fWeight);
      fRTriggerEneH->Fill(recEne/fEfactor, eneRatio, fWeight);
      
      if (fPart->LeadingFound()){
	float leaEne = fPart->GetE();
	float eneRatio2 = leaEne/genEne;
	fGPTriggerEneH->Fill(genEne, eneRatio2, fWeight);
	fPTriggerEneH->Fill(leaEne/0.2, eneRatio2, fWeight);
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Normalize histogrames 

void AliJetAnalysis::NormHistograms()
{
  // normalize shape histograms
  if (fDoShap) {
    if (fDoRecJ && fWShapR>0.) {  // leading reconstructed jet
      fRShapSelH->Scale(1.0/fWShapR);
      fRShapRejH->Scale(1.0/fWShapR);
      fRShapAllH->Scale(1.0/fWShapR);
    }
  }
  
  // normalize fragmentation function histograms
  if (fDoFrag) {
    if (fDoRecJ && fWFragR>0.) {  // leading reconstructed jet
      fRFragSelH->Scale(2.0/fWFragR);
      fRFragRejH->Scale(2.0/fWFragR);
      fRFragAllH->Scale(2.0/fWFragR);
    }
  }
}

////////////////////////////////////////////////////////////////////////

void AliJetAnalysis::PlotHistograms()
{
  // Plot histogramas (to be done...)
  if (fDoKine) PlotKineH();
  if (fDoCorr) PlotCorrH();
  if (fDoCorr50) PlotCorr50H();
  if (fDoShap) PlotShapH();
  if (fDoFrag) PlotFragH();
  if (fDoTrig) PlotTrigH();
}
 
void AliJetAnalysis::PlotKineH() const
{
    // missing    
    if (fDoPart) ;
    if (fDoGenJ) ;
    if (fDoRecJ) ;
}

void AliJetAnalysis::PlotCorrH() const
{
    // missing    
    if (fDoPart && fDoGenJ) ;
    if (fDoPart && fDoRecJ) ; 
    if (fDoGenJ && fDoRecJ) ; 
}
void AliJetAnalysis::PlotCorr50H() const
{
    // missing    
    if (fDoPart && fDoGenJ) ;
    if (fDoPart && fDoRecJ) ; 
    if (fDoGenJ && fDoRecJ) ; 
}

void AliJetAnalysis::PlotShapH() const
{
    // missing    
    if (fDoGenJ) ;
    if (fDoRecJ) ;
}

void AliJetAnalysis::PlotFragH() const
{
    // missing    
    if (fDoGenJ) ;
    if (fDoRecJ) ;
}

void AliJetAnalysis::PlotTrigH()
{
    // missing    

}

////////////////////////////////////////////////////////////////////////

void AliJetAnalysis::SaveHistograms()
{
  // Save histograms
  TFile *fOut = new TFile(fFile,"recreate");
  fOut->cd();
  if (fDoKine) SaveKineH();
  if (fDoCorr) SaveCorrH();
  if (fDoCorr50) SaveCorr50H();
  if (fDoShap) SaveShapH();
  if (fDoFrag) SaveFragH();
  if (fDoTrig) SaveTrigH();
  if (fDoJt) SaveJtH();
  if (fDodNdxi) SavedNdxiH();
  fOut->Close();
}

void AliJetAnalysis::SaveKineH()
{
  // Save kinematic histograms
  if (fDoPart) {
    fPKineEneH->Write();
    fPKinePtH->Write();
    fPKinePhiH->Write();
    fPKineEtaH->Write();
  }
  
  if (fDoGenJ) {
    fGKineEneH->Write();
    fGKinePtH->Write();
    fGKinePhiH->Write();
    fGKineEtaH->Write();
  }
  
  if (fDoRecJ) {
    fRKineEneH->Write();
    fRKinePtH->Write();
    fRKinePhiH->Write();
    fRKineEtaH->Write();
  }
}

void AliJetAnalysis::SaveCorrH()
{
  // Save correlation histograms
  if (fDoPart && fDoGenJ) {
    fPGCorrEneH->Write();
    fPGCorrPtH->Write();
    fPGCorrEtaH->Write();
    fPGCorrPhiH->Write();
  }
  
  if (fDoPart && fDoRecJ) {
    fPRCorrEneH->Write();
    fPRCorrPtH->Write();
    fPRCorrEtaH->Write();
    fPRCorrPhiH->Write();
  }
  
  if (fDoGenJ && fDoRecJ) {
    fRGCorrEneH->Write();
    fRGCorrPtH->Write();
    fRGCorrEtaH->Write();
    fRGCorrPhiH->Write();
  } 
}

void AliJetAnalysis::SaveCorr50H()
{
  // Save correlation histograms (special case)
  if (fDoPart && fDoRecJ) {
    fPRCorr50EneH->Write();
    fPRCorr50PtH->Write();
    fPRCorr50EtaH->Write();
    fPRCorr50PhiH->Write();
  }
  if (fDoGenJ && fDoRecJ) {
    fRGCorr50EneH->Write();
    fRGCorr50PtH->Write();
    fRGCorr50EtaH->Write();
    fRGCorr50PhiH->Write();
  } 
}

void AliJetAnalysis::SaveShapH()
{
  // Save jet shape histograms
  if (fDoRecJ) {
    fRShapSelH->Write();
    fdEdrH->Write();
    if(fDoBkgd) fdEdrB->Write();
    fPtEneH2->Write();
    fdEdrW->Write();
    if (fPythia){
      fRShapRejH->Write();
      fRShapAllH->Write();  
    }
  }    
}

void AliJetAnalysis::SaveJtH()
{
  // Save J_T histograms
  if (fDoRecJ) {
    fJtH->Write();
    if(fDoBkgd) fJtB->Write();
    fJtW->Write();
  }
}

void AliJetAnalysis::SavedNdxiH()
{
  // Save dN/d#xi histograms
  if (fDoRecJ) {
    fdNdxiH->Write();
    if(fDoBkgd) fdNdxiB->Write();
    fPtEneH->Write();
    fdNdxiW->Write();
  }
}

void AliJetAnalysis::SaveFragH()
{
  // Save fragmentation histograms
  if (fDoRecJ) {
    fRFragSelH->Write();
    if(fPythia){
      fRFragRejH->Write();
      fRFragAllH->Write();  
    }
  }
}

void AliJetAnalysis::SaveTrigH()
{
  // Save trigger bias histograms
  if(fDoTrig){
    fGTriggerEneH->Write();
    fRTriggerEneH->Write();
    fGPTriggerEneH->Write();
    fPTriggerEneH->Write();
  }
}

////////////////////////////////////////////////////////////////////////
// main Analysis function

void AliJetAnalysis::Analyze()
    
{
    // Kinematics for
    //   leading particle
    //   leading generated jet
    //   leading reconstructed jet

    // Correlations amd resolutions
    //    a) correlations in energy, pt, phi, eta
    //    b) resolutions in energy, pt, phi, eta, r 
    //   leading particle and leading generated jet
    //   leading particle and leading reconstructed jet
    //   leading generated jet and leading reconstructed jet

    // Fragmentation functions and Shapes
    //    a) integrated over all pt
    //    b) in 3 flavors:
    //       b.1) only for user selected particles in jet
    //       b.2) only for user rejected particles in jet
    //       b.3) for all particles in jet

    DefineHistograms();
    FillHistograms();
    NormHistograms();
    PlotHistograms();
    SaveHistograms();
}

