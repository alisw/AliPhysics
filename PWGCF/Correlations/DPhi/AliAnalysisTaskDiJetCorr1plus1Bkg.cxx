#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "THnSparse.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h" 
#include "AliCentrality.h"
#include "AliMultiplicity.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"

#include "AliEventPoolManager.h"

#include "AliAnalysisTaskDiJetCorr1plus1Bkg.h"

ClassImp(AliAnalysisTaskDiJetCorr1plus1Bkg)

//________________________________________________________________________
AliAnalysisTaskDiJetCorr1plus1Bkg::AliAnalysisTaskDiJetCorr1plus1Bkg(const char *name) 
  : AliAnalysisTaskSE(name)
, fAOD(0)
, centrality(0)
, fOutputList(0)
, fAodTracks(0)
, fAodTracksT(0)
, fAodTracksA(0)
, vertex(0)
, vtxSPD(0)
, fCentrOrMult(0)
, fBit(0)
, fTriggerpTLowThr(0)
, fTriggerpTHighThr(0)
, fHistTrigDPhi(0)
, fHistTrigDPhiM(0)
, fHistTrigPEtaS(0)
, fHistTrigSEtaS(0) 
, fHistTrigPEtaM(0)
, fHistTrigSEtaM(0)
, fHistDeltaPhiT1T2BC(0)
, fHistDeltaPhiT1T2AC(0)
, fHistDeltaPhiT1T2C1(0)
, fHistDeltaEtaT1T2BC(0)
, fHistDeltaEtaT1T2AC(0)

//, fHistEff(0)
, fEventCounter(0)
, fHistCent(0)
, f3DEffCor(0)
, fControlConvResT1(0)
, fControlConvResT2(0)
, fControlConvResMT1(0)
, fControlConvResMT2(0)
, fTHnCentZvtxDEtaDPhi1SE(0)
//, fTHnCentZvtxDEtaDPhi2SE(0)
, fTHnTrigCentZvtxpTtrig1SE(0)
//, fTHnTrigCentZvtxpTtrig2SE(0)
, fTHnCentZvtxDEtaDPhi1ME(0)
//, fTHnCentZvtxDEtaDPhi2ME(0)
, fTHnTrigCentZvtxpTtrig1ME(0)
, fThnEff(0)
, fEffCheck(0)
, fNoMixedEvents(0)
, fMixStatCentorMult(0)
, fMixStatZvtx(0)
//, fTHnTrigCentZvtxpTtrig2ME(0)
, fTrackArray(0)
, fPoolMgr(0x0)
, fFillMixed(kTRUE)
, fMixingTracks(15000)//i ran on 25000 for without resonance correction
, useVarBins(kTRUE)
, fCutConversions(kTRUE)
, fCutResonances(kTRUE)
, twoTrackEfficiencyCut(kTRUE)
, fEtaOrdering(kTRUE)
, fSetSystemValue(kTRUE)


{
  // Constructor
  for ( Int_t i = 0; i < 9; i++) { 
    fHistQA[i] = NULL;
  }
  // Define input and output slots here
  // Input slot #0 works with a TChain.
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________

AliAnalysisTaskDiJetCorr1plus1Bkg::~AliAnalysisTaskDiJetCorr1plus1Bkg() 

{
  delete centrality;
  delete fOutputList;
  delete vtxSPD;
  delete vertex;
  if(fThnEff) {delete fThnEff; fThnEff = 0;}
   

}

//________________________________________________________________________
void AliAnalysisTaskDiJetCorr1plus1Bkg::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputList = new TList();
  fOutputList->SetOwner();
  
  fEventCounter = new TH1F("fEventCounter","EventCounter",10, 0, 10);
  // fEventCounter->GetXaxis()->SetBinLabel(1,"Total Events Accessed");
  //fEventCounter->GetXaxis()->SetBinLabel(4,"After vertex cut");
  //fEventCounter->GetXaxis()->SetBinLabel(7,"Events Analysed");
  fEventCounter->Sumw2();
  fOutputList->Add(fEventCounter);
    
    
    
  if (fSetSystemValue) fHistCent = new TH1F("fHistCent", "centrality distribution", 100, 0, 100);
    
  if(!fSetSystemValue) fHistCent = new TH1F("fHistCent", "centrality distribution", 100, 0, 250);
  fHistCent->Sumw2();
  fOutputList->Add(fHistCent);

  
  //QA histograms
  fHistQA[0] = new TH1F("fHistVtx", "Z vertex distribution", 100, -11., 11.);
  fHistQA[1] = new TH1F("fHistnTpcCluster", "n TPC Cluster", 100, 0., 200.);
  fHistQA[2] = new TH1F("fHistnTpcClusterF", "n TPC Cluster findable", 100, 0., 200.);
  fHistQA[3] = new TH1F("fHistDcaXY", "dca x y", 10, -5., 5.);
  fHistQA[4] = new TH1F("fHistDcaZ", "dca Z", 10, -5., 5.);
  fHistQA[5] = new TH1F("fHistChi2Tpc", "Chi2 TPC", 100, 0., 10.);
  fHistQA[6] = new TH1F("fHistpT", "pT distribution",100 ,0.,20.);
  fHistQA[7] = new TH1F("fHistPhi", "Phi distribution" , 36, -0.5, 2*TMath::Pi()+0.5);
  fHistQA[8] = new TH1F("fHistEta", "Eta distribution" , 18, -1.8, 1.8);

  for( Int_t i = 0; i < 9; i++)
    {
      fHistQA[i]->Sumw2();
      fOutputList->Add(fHistQA[i]);
    }

    fEffCheck = new TH1F("fEffCheck","eff values: for check",10, 0, 10.0);
    fOutputList->Add(fEffCheck);
    
    fNoMixedEvents = new TH1F("fNoMixedEvents","Mixed event stat",1, 0, 1) ;
    fOutputList->Add(fNoMixedEvents);
    
    fMixStatCentorMult = new TH2F("fMixStatCentorMult","no of events in pool  vs Centrality;Nevent in pool;centOrMul",50,0,200,100,0,200);
    fOutputList->Add(fMixStatCentorMult);
    
    fMixStatZvtx = new TH2F("fMixStatZvtx","no of events in pool  vs zvtx;Nevents in pool;zvtx",50,0,200,10,-10,10);
    fOutputList->Add(fMixStatZvtx);

  
  //Histos for same events.


  //fHistTrigDPhi = new TH1F("fHistTrigDPhi", " Trig Phi Difference Same",100, 0, 0.5);
  //fHistTrigDPhiM = new TH1F("fHistTrigDPhiM", " Trig Phi Difference Mixed ",100, 0, 0.5);
 // fHistTrigPEtaS = new TH1F("fHistTrigPEtaS", " Prim trig - same events",100, -1, 1);
 // fHistTrigSEtaS = new TH1F("fHistTrigSEtaS", " Sec trig - same events ",100, -1, 1);
 // fHistTrigPEtaM = new TH1F("fHistTrigPEtaM", " Prim trig - mixed events",100, -1, 1);
 // fHistTrigSEtaM = new TH1F("fHistTrigSEtaM", " Sec trig - mixed events ",100, -1, 1);
  //fHistDeltaPhiT1T2BC = new TH3F("fHistDeltaPhiT1T2BC", " TrgPhi1-TrgPhi2 ",36, -1*(0.5)*TMath::Pi(), 3*(0.5)*TMath::Pi(), 11, 0., 100., 30, 0., 30.);
//  fHistDeltaPhiT1T2C1 = new TH2F("fHistDeltaPhiT1T2C1", " TrgPhi1-TrgPhi2 ",36, -1*(0.5)*TMath::Pi(), 3*(0.5)*TMath::Pi(), 11, 0., 100.);

//  fHistDeltaPhiT1T2AC = new TH2F("fHistDeltaPhiT1T2AC", " TrgPhi1-TrgPhi2 ",36, -1*(0.5)*TMath::Pi(), 3*(0.5)*TMath::Pi(), 11, 0., 100.);
 
  //fHistDeltaEtaT1T2BC = new TH3F("fHistDeltaEtaT1T2BC", " TrgEta1-TrgEta2 ",36, -1.8, 1.8, 4, 0., 20., 30, 0., 30.);
 // fHistDeltaEtaT1T2AC = new TH2F("fHistDeltaEtaT1T2AC", " TrgEta1-TrgEta2 ",36, -1.8, 1.8, 10, 0., 20.);
    
  /*fOutputList->Add(fHistTrigPEtaS);
  fOutputList->Add(fHistTrigSEtaS);
  fOutputList->Add(fHistTrigPEtaM);
  fOutputList->Add(fHistTrigSEtaM);*/
  //fOutputList->Add(fHistDeltaPhiT1T2BC);
  //fOutputList->Add(fHistDeltaPhiT1T2AC);
  //fOutputList->Add(fHistDeltaPhiT1T2C1);

  //fOutputList->Add(fHistDeltaEtaT1T2BC);
 // fOutputList->Add(fHistDeltaEtaT1T2AC);
    

  //fHistEff = new TH1D("fHistEff", " Eff histo ",1000, 0, 100);

  //fOutputList->Add(fHistEff);
 // fOutputList->Add(fHistTrigDPhi);
  //if(fFillMixed)fOutputList->Add(fHistTrigDPhiM);

  fControlConvResT1 = new TH2F("fControlConvResT1", ";id;delta mass;T1", 3, -0.5, 2.5, 100, -0.1, 0.1);
  fControlConvResT2 = new TH2F("fControlConvResT2", ";id;delta mass;T2", 3, -0.5, 2.5, 100, -0.1, 0.1);
  if(fFillMixed)
{
  fControlConvResMT1 = new TH2F("fControlConvResMT1", ";id;delta mass;MT1", 3, -0.5, 2.5, 100, -0.1, 0.1);
  fControlConvResMT2 = new TH2F("fControlConvResMT2", ";id;delta mass;MT2", 3, -0.5, 2.5, 100, -0.1, 0.1);
 }

  fOutputList->Add(fControlConvResT1);
  fOutputList->Add(fControlConvResT2);
  fOutputList->Add(fControlConvResMT1);
  fOutputList->Add(fControlConvResMT2);
  /*********************************************************/

   /*******************************************************/

    Int_t nBinsCentorMult = 0; Double_t fMinCentorMult = 0.0; Double_t fMaxCentorMult = 0.0;
    
    if(fSetSystemValue){
        nBinsCentorMult = 12; fMinCentorMult = 0.0, fMaxCentorMult = 100.0;}
    
    if(!fSetSystemValue){
         nBinsCentorMult = 2;  fMinCentorMult = 0.0;  fMaxCentorMult = 250.0;}
    

    const Int_t pTbinTrigger1plus1 = Int_t(fTriggerpTHighThr - fTriggerpTLowThr);
    Int_t   fBinsTrg1plus1[3]   = {nBinsCentorMult,       10,   pTbinTrigger1plus1};
    Double_t fMinTrg1plus1[3]   = {fMinCentorMult,   -10.0,   fTriggerpTLowThr};
    Double_t fMaxTrg1plus1[3]   = {fMaxCentorMult,  10.0,   fTriggerpTHighThr};
    
    
    fTHnTrigCentZvtxpTtrig1SE = new THnSparseD("fTHnTrigCentZvtxpTtrig1SE","Cent-Zvtx-pTtr1",3, fBinsTrg1plus1, fMinTrg1plus1, fMaxTrg1plus1);
   
    
    
    Int_t   fBins121plus1[6] = {nBinsCentorMult,     10,   18,                36,  pTbinTrigger1plus1,  10};
    Double_t  fMin121plus1[6] = {fMinCentorMult,   -10., -1.8, -0.5*TMath::Pi(), fTriggerpTLowThr,  0.5};
    Double_t  fMax121plus1[6] = {fMaxCentorMult,  10.,  1.8,  1.5*TMath::Pi(), fTriggerpTHighThr, 10};
    
    fTHnCentZvtxDEtaDPhi1SE = new THnSparseD("fTHnCentZvtxDEtaDPhi1SE","Cent-zVtx-DEta1-DPhi1-pTT-pTA",6, fBins121plus1, fMin121plus1, fMax121plus1);
    
  

  
  
  
  if(fFillMixed){

    fTHnCentZvtxDEtaDPhi1ME = new THnSparseD("fTHnCentZvtxDEtaDPhi1ME","Cent-zVtx-DEta1-DPhi1-pTTrig1-pTAsso1-MixEvt",6, fBins121plus1, fMin121plus1, fMax121plus1);

  

    fTHnTrigCentZvtxpTtrig1ME = new THnSparseD("fTHnTrigCentZvtxpTtrig1ME","Cent-zvtx-pTtrig1-MixEvt",3, fBinsTrg1plus1, fMinTrg1plus1, fMaxTrg1plus1);

    
}

 //Set variable bins
    if(fSetSystemValue){

  if(useVarBins){
   
      const Int_t nvarBinsCent = 12;
      Double_t varBinsCent[nvarBinsCent+1] = {0., 1., 2., 3., 4., 5., 7.5, 10., 20., 30., 40., 50., 100.1};
      
     // fHistDeltaPhiT1T2BC->GetYaxis()->Set(nvarBinsCent, varBinsCent);
     // fHistDeltaPhiT1T2C1->GetYaxis()->Set(nvarBinsCent, varBinsCent);
     // fHistDeltaPhiT1T2AC->GetYaxis()->Set(nvarBinsCent, varBinsCent);

    fTHnCentZvtxDEtaDPhi1SE->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
    fTHnTrigCentZvtxpTtrig1SE->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
      

    if(fFillMixed) {
      fTHnCentZvtxDEtaDPhi1ME->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
      fTHnTrigCentZvtxpTtrig1ME->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
        

    }   
  }
    }
    
    
    
    
    if(!fSetSystemValue){
        
        if(useVarBins){
            
            const Int_t nvarBinsCent = 2;
            Double_t varBinsCent[nvarBinsCent+1] = {0., 35., 250.};
            
            fTHnCentZvtxDEtaDPhi1SE->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
            fTHnTrigCentZvtxpTtrig1SE->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
            
            
            if(fFillMixed) {
                fTHnCentZvtxDEtaDPhi1ME->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
                fTHnTrigCentZvtxpTtrig1ME->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
                
                
            }
        }
    }
    
    
    
    const Int_t nvarBinspT = 10;
    Double_t varBinspT[nvarBinspT+1] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
    
      fTHnCentZvtxDEtaDPhi1SE->GetAxis(5)->Set(nvarBinspT, varBinspT);
    
      
      if(fFillMixed)
      fTHnCentZvtxDEtaDPhi1ME->GetAxis(5)->Set(nvarBinspT, varBinspT);
    
  fTHnCentZvtxDEtaDPhi1SE->Sumw2();
  fTHnTrigCentZvtxpTtrig1SE->Sumw2();
    
  fOutputList->Add(fTHnCentZvtxDEtaDPhi1SE);
  fOutputList->Add(fTHnTrigCentZvtxpTtrig1SE);
  
  if(fFillMixed){
      
  fTHnCentZvtxDEtaDPhi1ME->Sumw2();
  fTHnTrigCentZvtxpTtrig1ME->Sumw2();
  fOutputList->Add(fTHnCentZvtxDEtaDPhi1ME);
  fOutputList->Add(fTHnTrigCentZvtxpTtrig1ME);
  }
 // pool for event mixing
  
 // Int_t trackDepth = 50000;
    
    if(fFillMixed){
        if(fSetSystemValue) {Bool_t DefPool = DefineMixedEventPoolPbPb();
            if(!DefPool){
                AliInfo("UserCreateOutput: Pool is not define properly");
                return;
            }}
        
        
        if(!fSetSystemValue) {Bool_t DefPool = DefineMixedEventPoolpp();
            if(!DefPool){
                AliInfo("UserCreateOutput: Pool is not define properly");
                return;
            }}
    }

    
    
  // fPoolMgr->SetTargetValues(fMixingTracks, 0.1, 5);
  // fPoolMgr->SetTargetValues(fMixingTracks, 0.1, 5);
}

//________________________________________________________________________
void  AliAnalysisTaskDiJetCorr1plus1Bkg::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Post output data.

  Int_t counter = 0;
  TObjArray* fTrackArray = new TObjArray;
  fTrackArray->SetOwner();

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;
  }
 
 fEventCounter->Fill(1);
 Float_t bSign = 0;
 if(fAOD)
   {
     bSign = (fAOD->GetMagneticField() > 0) ? 1 : -1;
   }
  
    
    
    
    AliCentrality *centralityObj = 0x0;
    if(fSetSystemValue){ // pPb, PbPb
        centralityObj = ((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP();
        fCentrOrMult = centralityObj->GetCentralityPercentileUnchecked("V0M");
        if(centralityObj->GetQuality()!=0) return;
        if((abs(fCentrOrMult)) < 0. || (abs(fCentrOrMult)) > 100.1)return;
    }
    else if(!fSetSystemValue){ // pp, pPb
        Double_t count = -1, mineta = -1.0, maxeta = 1.0;
        AliAODTracklets* tracklets = fAOD->GetTracklets();
        Int_t nTr=tracklets->GetNumberOfTracklets();
        for(Int_t iTr=0; iTr<nTr; iTr++){
            Double_t theta=tracklets->GetTheta(iTr);
            Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
            if(eta>mineta && eta<maxeta) count++;
        }
        fCentrOrMult = count;
    }

    
    
  fHistCent->Fill(fCentrOrMult);
  
  const AliAODVertex* vertex = fAOD->GetPrimaryVertex();
  if(!vertex || vertex->GetNContributors()<=0) return;
  TString vtxTtl = vertex->GetTitle();
  if(!vtxTtl.Contains("VertexerTracks")) return;
  Float_t zVertex = vertex->GetZ();
  const AliAODVertex* vtxSPD = fAOD->GetPrimaryVertexSPD();
  if(vtxSPD->GetNContributors()<=0) return;
  TString vtxTyp = vtxSPD->GetTitle();
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if(vtxTyp.Contains("vertexer:Z") && (zRes > 0.25)) return;
  if(TMath::Abs(vtxSPD->GetZ() - vertex->GetZ()) > 0.5) return;
  if (TMath::Abs(zVertex) > 10) return;
  //cout << "zvertex"<< zVertex << endl;
  fHistQA[0]->Fill(zVertex);
  fEventCounter->Fill(4);



  
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) 
    
     {
      
      AliAODTrack* fAodTracks = (AliAODTrack*)fAOD->GetTrack(iTracks);
      
      if (!fAodTracks)continue;
      
      if(!fAodTracks->TestFilterBit(fBit)) continue;
      
      
      if(fAodTracks->Eta() < -0.9 || fAodTracks->Eta() > 0.9){

	continue; }
      
    /*  if (fAodTracks->Pt() < 0.5 || fAodTracks->Pt() > 20.){
	
     continue;}*/
         
         
         
        if (fAodTracks->Pt() < 0.5 || fAodTracks->Pt() > 20.){
	
          continue;}
      
      
     
      if (!fAodTracks)continue;
      
      fHistQA[1]->Fill(fAodTracks->GetTPCClusterInfo(2,1)); 
      
      fHistQA[3]->Fill(fAodTracks->DCA());
      fHistQA[4]->Fill(fAodTracks->ZAtDCA());
      Double_t chi2Tpc = fAodTracks->Chi2perNDF();
     
      fHistQA[5]->Fill(chi2Tpc);
      fHistQA[6]->Fill(fAodTracks->Pt());
      fHistQA[7]->Fill(fAodTracks->Phi());
      fHistQA[8]->Fill(fAodTracks->Eta());  
      
      fTrackArray->Add(fAodTracks);
    }
    

    for(int k=0; k < fTrackArray->GetEntries(); k++)
    
    {
		
      TObject* obj1 = fTrackArray->At(k);
      
        AliAODTrack* fAodTracksT = (AliAODTrack*) obj1;
        
     //   if (fAodTracksT->Pt() <= 20. && fAodTracksT->Pt() >= 4.){
            
            if (fAodTracksT->Pt() <= fTriggerpTHighThr && fAodTracksT->Pt() >= fTriggerpTLowThr){
                
            Double_t effvalueT1 = GetTrackWeight( fAodTracksT->Eta(), fAodTracksT->Pt(),fCentrOrMult, zVertex);

            
            Double_t fCentZvtxpT[3] = {fCentrOrMult, zVertex, fAodTracksT->Pt()};
            fTHnTrigCentZvtxpTtrig1SE->Fill(fCentZvtxpT,effvalueT1);
        for(int m=0;  m< fTrackArray->GetEntries(); m++){
            
         TObject* obj2 = fTrackArray->At(m);
      
            AliAODTrack* fAodTracksAS = (AliAODTrack*) obj2;
    
        
            if (fAodTracksAS->Pt() >= fAodTracksT->Pt()) continue;
	  
	//conversions
	     if (fCutConversions && (fAodTracksAS->Charge() * fAodTracksT->Charge()) < 0)
	      {
              Float_t mass = GetInvMassSquaredCheap(fAodTracksT->Pt(), fAodTracksT->Eta(), fAodTracksT->Phi(), fAodTracksAS->Pt(), fAodTracksAS->Eta(), fAodTracksAS->Phi(), 0.510e-3, 0.510e-3);
	         if (mass < 0.1)
	          {
                  mass = GetInvMassSquared(fAodTracksT->Pt(), fAodTracksT->Eta(), fAodTracksT->Phi(), fAodTracksAS->Pt(), fAodTracksAS->Eta(), fAodTracksAS->Phi(), 0.510e-3, 0.510e-3);
		          if(mass < 0.04 * 0.04)
		          continue;
	           }
	      
	       }

	//K0s
	if (fCutResonances && (fAodTracksAS->Charge() * fAodTracksT->Charge()) < 0)
	  {
	    Float_t mass = GetInvMassSquaredCheap(fAodTracksT->Pt(), fAodTracksT->Eta(), fAodTracksT->Phi(), fAodTracksAS->Pt(), fAodTracksAS->Eta(), fAodTracksAS->Phi(), 0.1396, 0.1396);
	    
	    const Float_t kK0smass = 0.4976;
	    if (TMath::Abs(mass -kK0smass * kK0smass)  < 0.1)
	      {
		mass = GetInvMassSquared(fAodTracksT->Pt(), fAodTracksT->Eta(), fAodTracksT->Phi(), fAodTracksAS->Pt(), fAodTracksAS->Eta(), fAodTracksAS->Phi(), 0.1396, 0.1396);
		fControlConvResT1->Fill(1, mass -kK0smass * kK0smass); 
		if(mass > (kK0smass-0.02)*(kK0smass-0.02) && mass < (kK0smass+0.02)*(kK0smass+0.02))
		  continue;
	      }
	  }
	
	//lambda

	if (fCutResonances && (fAodTracksAS->Charge() * fAodTracksT->Charge()) < 0)
	  {
	    Float_t mass1 = GetInvMassSquaredCheap(fAodTracksT->Pt(), fAodTracksT->Eta(), fAodTracksT->Phi(), fAodTracksAS->Pt(), fAodTracksAS->Eta(), fAodTracksAS->Phi(), 0.1396, 0.9383);
	    Float_t mass2 = GetInvMassSquaredCheap(fAodTracksT->Pt(), fAodTracksT->Eta(), fAodTracksT->Phi(), fAodTracksAS->Pt(), fAodTracksAS->Eta(), fAodTracksAS->Phi(), 0.9383, 0.1396);
	    
	    const Float_t kLambdaMass = 1.115;
	    if (TMath::Abs(mass1 -kLambdaMass * kLambdaMass)  < 0.1)
	      {
		mass1 = GetInvMassSquared(fAodTracksT->Pt(), fAodTracksT->Eta(), fAodTracksT->Phi(), fAodTracksAS->Pt(), fAodTracksAS->Eta(), fAodTracksAS->Phi(), 0.1396, 0.9383);
		fControlConvResT1->Fill(2, mass1 - kLambdaMass * kLambdaMass);	
		if(mass1 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass1 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
		  continue;
	      }

	    if (TMath::Abs(mass2 -kLambdaMass * kLambdaMass)  < 0.1)
	      {
		mass2 = GetInvMassSquared(fAodTracksT->Pt(), fAodTracksT->Eta(), fAodTracksT->Phi(), fAodTracksAS->Pt(), fAodTracksAS->Eta(), fAodTracksAS->Phi(), 0.1396, 0.9383);
		fControlConvResT1->Fill(2, mass2 - kLambdaMass * kLambdaMass);	
		if(mass2 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass2 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
		  continue;
	      }
	    
	    
	  }
	///////////
	if (twoTrackEfficiencyCut)
	  {
	    Float_t phi1 = fAodTracksT->Phi();
	    Float_t pt1 = fAodTracksT->Pt();
	    Float_t charge1 = fAodTracksT->Charge();

	    Float_t phi2 = fAodTracksAS->Phi();
	    Float_t pt2 = fAodTracksAS->Pt();
	    Float_t charge2 = fAodTracksAS->Charge();

	    Float_t deta = fAodTracksT->Eta() - fAodTracksAS->Eta();

	    if (TMath::Abs(deta) < 0.02 * 2.5 * 3)
	    {
	    // check first boundaries to see if is worth to loop and find the minimum
	    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
	    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);
	  
	    const Float_t kLimit = 0.02 * 3;

	    Float_t dphistarminabs = 1e5;
	    Float_t dphistarmin = 1e5;
	    if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0)
	    {
	      for (Double_t rad=0.8; rad<2.51; rad+=0.01) 
	      {
		Float_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSign);

		Float_t dphistarabs = TMath::Abs(dphistar);
		
		if (dphistarabs < dphistarminabs)
		{
		  dphistarmin = dphistar;
		  dphistarminabs = dphistarabs;
		}
	      }
	      
	    
	      if (dphistarminabs < 0.02 && TMath::Abs(deta) < 0.02)
	      {
		//Printf("Trigger1 Removed track pair with %f %f %f %f %f %f %f %f %f", deta, dphistarminabs, phi1, pt1, charge1, phi2, pt2, charge2, bSign);
		continue;
	      }

    	      
	    }
	  }

	  }

	///////////
           
      Double_t deltaPhi1 = fAodTracksT->Phi() - fAodTracksAS->Phi();
      if (deltaPhi1 > 1.5 * TMath::Pi()) 
	  deltaPhi1 -= TMath::TwoPi();
      
      if (deltaPhi1 < -0.5 * TMath::Pi())
	  deltaPhi1 += TMath::TwoPi();
            

      //fHistDphi1->Fill(deltaPhi1);
     Double_t deltaEta1 = fAodTracksT->Eta() - fAodTracksAS->Eta();
            
         // Double_t ptLim_Sparse1 = 0.0 ;
     Double_t ptTrack1 = fAodTracksAS->Pt();
            
            
     Double_t effvalueAS = GetTrackWeight(fAodTracksAS->Eta(), ptTrack1, fCentrOrMult, zVertex);
    
     Double_t effvalue = effvalueT1*effvalueAS;
            
         //   cout<<"effvalue"<<effvalue<<endl;
            
     fEffCheck->Fill(effvalue);


        Double_t fCentzVtxDEtaDPhiSE1[6] = {fCentrOrMult, zVertex, deltaEta1, deltaPhi1, fAodTracksT->Pt(), ptTrack1};
       // Double_t effvalue = 1.0;//GetTrackbyTrackEffValue(fAodTracksAS, fCentrOrMult);
      // cout<<"Efficiency value = "<<effvalue<<endl;
            
           // cout<<"deltaeta SE:"<<" "<<deltaEta1<<" "<<"deltaphi3 SE:"<<" "<<deltaPhi1<<endl;
   // cout<<"   centrality:"<<fCentrOrMult<<"   zvertex:"<<zVertex<<" "<<"  deltaEta1:"<<deltaEta1<<"  deltaPhi1"<<deltaPhi1<<"  trig pT S:"<<fAodTracksT->Pt()<<"  ass pT S:"<<ptTrack1<<"  effvalue:"<<effvalue<<endl;
      
     fTHnCentZvtxDEtaDPhi1SE->Fill(fCentzVtxDEtaDPhiSE1, effvalue);
        }
      
      
    }
    }
  
  if (fFillMixed)
    
  {
      
      
      AliEventPool* pool = fPoolMgr->GetEventPool(fCentrOrMult, zVertex);
      
    
      if (!pool){
          
          //	printf("No pool found for centrality = %f, zVertex = %f\n", fCentrOrMult, zVertex);
          
          return;
      }
      
      if (pool->IsReady())
          
      {
          //cout<<"pool found"<<endl;
          if (!fTrackArray)return;
          
          Int_t nMix = pool->GetCurrentNEvents();
          fNoMixedEvents->Fill(0);
          fMixStatCentorMult->Fill(nMix,fCentrOrMult);
          fMixStatZvtx->Fill(nMix,zVertex);

          
         // if(nMix < 5) return;
          
          for (Int_t jMix=0; jMix < nMix; jMix++) {
              
              TObjArray* bgTracks = pool->GetEvent(jMix);
              
              if(!bgTracks)return;
              
             
              
              for(int i=0; i<fTrackArray->GetEntries(); i++)
                  
              {
                  
                  TObject* obj3 = fTrackArray->At(i);
                  
                  
                  AliAODTrack* fAodTracksT3 = (AliAODTrack*) obj3;
                  
                //  if (fAodTracksT3->Pt() <= 20. && fAodTracksT3->Pt() >= 4.){
                      
            if (fAodTracksT3->Pt() <= fTriggerpTHighThr && fAodTracksT3->Pt() >= fTriggerpTLowThr){
                      
              
              for(int l=0; l < bgTracks->GetEntries(); l++)
              {
                  TObject* obj4 = bgTracks->At(l);
                  
                  AliAODTrack* fAodTracksASM = (AliAODTrack*) obj4;
                  
                  if(!fAodTracksASM) continue;
                  
                  if (fAodTracksASM->Pt() >= fAodTracksT3->Pt()) continue;
                  
                  
                      //conversions
                      if (fCutConversions && (fAodTracksASM->Charge() * fAodTracksT3->Charge()) < 0)
                      {
                          Float_t mass = GetInvMassSquaredCheap(fAodTracksT3->Pt(), fAodTracksT3->Eta(), fAodTracksT3->Phi(), fAodTracksASM->Pt(), fAodTracksASM->Eta(), fAodTracksASM->Phi(), 0.510e-3, 0.510e-3);
                          if (mass < 0.1)
                          {
                              mass = GetInvMassSquared(fAodTracksT3->Pt(), fAodTracksT3->Eta(), fAodTracksT3->Phi(), fAodTracksASM->Pt(), fAodTracksASM->Eta(), fAodTracksASM->Phi(), 0.510e-3, 0.510e-3);
                              
                              fControlConvResMT1->Fill(0.0, mass);
                              
                              if(mass < 0.04 * 0.04)
                                  continue;
                          }
                      }
                      
                      //K0s
                     if (fCutResonances && (fAodTracksASM->Charge() * fAodTracksT3->Charge()) < 0)
                      {
                          Float_t mass = GetInvMassSquaredCheap(fAodTracksT3->Pt(), fAodTracksT3->Eta(), fAodTracksT3->Phi(), fAodTracksASM->Pt(), fAodTracksASM->Eta(), fAodTracksASM->Phi(), 0.1396, 0.1396);
                          
                          const Float_t kK0smass = 0.4976;
                          if (TMath::Abs(mass -kK0smass * kK0smass)  < 0.1)
                          {
                              mass = GetInvMassSquared(fAodTracksT3->Pt(), fAodTracksT3->Eta(), fAodTracksT3->Phi(), fAodTracksASM->Pt(), fAodTracksASM->Eta(), fAodTracksASM->Phi(), 0.1396, 0.1396);
                              fControlConvResMT1->Fill(1, mass -kK0smass * kK0smass);
                              if(mass > (kK0smass-0.02)*(kK0smass-0.02) && mass < (kK0smass+0.02)*(kK0smass+0.02))
                                  continue;
                          }
                      }
                      
                      //lambda
                  
                      if (fCutResonances && (fAodTracksASM->Charge() * fAodTracksT3->Charge()) < 0)
                      {
                          Float_t mass1 = GetInvMassSquaredCheap(fAodTracksT3->Pt(), fAodTracksT3->Eta(), fAodTracksT3->Phi(), fAodTracksASM->Pt(), fAodTracksASM->Eta(), fAodTracksASM->Phi(), 0.1396, 0.9383);
                          Float_t mass2 = GetInvMassSquaredCheap(fAodTracksT3->Pt(), fAodTracksT3->Eta(), fAodTracksT3->Phi(), fAodTracksASM->Pt(), fAodTracksASM->Eta(), fAodTracksASM->Phi(), 0.9383, 0.1396);
                          
                          const Float_t kLambdaMass = 1.115;
                          if (TMath::Abs(mass1 -kLambdaMass * kLambdaMass)  < 0.1)
                          {
                              mass1 = GetInvMassSquared(fAodTracksT3->Pt(), fAodTracksT3->Eta(), fAodTracksT3->Phi(), fAodTracksASM->Pt(), fAodTracksASM->Eta(), fAodTracksASM->Phi(), 0.1396, 0.9383);
                              fControlConvResMT1->Fill(2, mass1 - kLambdaMass * kLambdaMass);
                              if(mass1 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass1 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
                                  continue;
                          }
                          
                          if (TMath::Abs(mass2 -kLambdaMass * kLambdaMass)  < 0.1)
                          {
                              mass2 = GetInvMassSquared(fAodTracksT3->Pt(), fAodTracksT3->Eta(), fAodTracksT3->Phi(), fAodTracksASM->Pt(), fAodTracksASM->Eta(), fAodTracksASM->Phi(), 0.1396, 0.9383);
                              fControlConvResMT1->Fill(2, mass1 - kLambdaMass * kLambdaMass);
                              if(mass2 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass2 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
                                  continue;
                          }
                          
                          
                      }
                 
                  if (twoTrackEfficiencyCut)
                      {
                          Float_t phi1 = fAodTracksT3->Phi();
                          Float_t pt1 = fAodTracksT3->Pt();
                          Float_t charge1 = fAodTracksT3->Charge();
                          
                          Float_t phi2 = fAodTracksASM->Phi();
                          Float_t pt2 = fAodTracksASM->Pt();
                          Float_t charge2 = fAodTracksASM->Charge();
                          
                          Float_t deta = fAodTracksT3->Eta() - fAodTracksASM->Eta();
                          
                          if (TMath::Abs(deta) < 0.02 * 2.5 * 3)
                          {
                              // check first boundaries to see if is worth to loop and find the minimum
                              Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
                              Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);
                              
                              const Float_t kLimit = 0.02 * 3;
                              
                              Float_t dphistarminabs = 1e5;
                              Float_t dphistarmin = 1e5;
                              if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0)
                              {
                                  for (Double_t rad=0.8; rad<2.51; rad+=0.01)
                                  {
                                      Float_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSign);
                                      
                                      Float_t dphistarabs = TMath::Abs(dphistar);
                                      
                                      if (dphistarabs < dphistarminabs)
                                      {
                                          dphistarmin = dphistar;
                                          dphistarminabs = dphistarabs;
                                      }
                                  }
                                  
                                  
                                  if (dphistarminabs < 0.02 && TMath::Abs(deta) < 0.02)
                                  {
                                      //  Printf("Removed track pair with %f %f %f %f %f %f %f %f %f", deta, dphistarminabs, phi1, pt1, charge1, phi2, pt2, charge2, bSign);
                                      continue;
                                  }
                                  
                                  
                              }
                          }
                          
                      }
                      
                   
                   
                      Double_t deltaPhi3 = fAodTracksT3->Phi() - fAodTracksASM->Phi();
                      if (deltaPhi3 > 1.5 * TMath::Pi())
                          deltaPhi3 -= TMath::TwoPi();
                      
                      if (deltaPhi3 < -0.5 * TMath::Pi())
                          deltaPhi3 += TMath::TwoPi();
                      
                      //fHistDphiM1->Fill(deltaPhi3);
                      
                      
                      Double_t deltaEta3 = fAodTracksT3->Eta() - fAodTracksASM->Eta();
                      
                     // Double_t ptLim_Sparse1ME = 0.0 ;
                      Double_t ptTrack1ME = fAodTracksASM->Pt();
                 
                  Double_t effvalueT1 = GetTrackWeight( fAodTracksT3->Eta(), fAodTracksT3->Pt(),fCentrOrMult, zVertex);
                  
                  Double_t effvalueAS = GetTrackWeight( fAodTracksASM->Eta(), ptTrack1ME, fCentrOrMult, zVertex);
                  
                  Double_t effvalue = effvalueT1*effvalueAS;
                  
                  fEffCheck->Fill(effvalue);
                  
                      Double_t fCentzVtxDEtaDPhiME1[6] = {fCentrOrMult, zVertex, deltaEta3, deltaPhi3, fAodTracksT3->Pt(), ptTrack1ME};
                  
                 // cout<<"  centrality:"<<fCentrOrMult<<"  zvertex:"<<zVertex<<"  deltaEta3:"<<deltaEta3<<"  deltaPhi3"<<deltaPhi3<<"  trig pT S:"<<fAodTracksT3->Pt()<<"  ass pT S:"<<ptTrack1ME<<"  effvalue:"<<effvalue<<endl;
                      //Double_t effvalue = 1.0;//GetTrackbyTrackEffValue(fAodTracksASM, fCentrOrMult);
                      fTHnCentZvtxDEtaDPhi1ME->Fill(fCentzVtxDEtaDPhiME1, effvalue);
                      
                      
                  
                  
                      // fHistDetaDphiM1->Fill(deltaEta3,deltaPhi3);
                  
                  
              }
				  
                  }
          }	 
              }
      }
      
      // clone and give ownership to event pool
      TObjArray* tracksClone = (TObjArray*) fTrackArray->Clone();
      tracksClone->SetOwner(kTRUE);
      
      pool->UpdatePool(tracksClone);
      // delete tracksClone;
      
      //pool->PrintInfo();
      
  }
    
    if(!fFillMixed){
        
        delete fTrackArray;
        fTrackArray = NULL;
    }
  
  fEventCounter->Fill(7);
  PostData(1, fOutputList);
  //delete f3DEffCor;
  
    } 
  



  //________________________________________________________________________
  void AliAnalysisTaskDiJetCorr1plus1Bkg::Terminate(Option_t *) 
  {
    // Draw result to the screen
    // Called once at the end of the query
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));

    if (!fOutputList) {
      printf("ERROR: Output list not available\n");
      return;
    }
    
    
  }


//______________________________|  Track Efficiency


Double_t AliAnalysisTaskDiJetCorr1plus1Bkg::GetTrackWeight(Double_t eta, Double_t pt, Double_t CentrOrMult, Double_t zVertex){
    
    // cout<<" histo.....returning............"<<endl;
    Double_t efficiency = 0;
    if(!fThnEff){
        
        //cout<<"No histo.....returning............"<<endl;
        
        return 1;
        
    }
    
    //  TAxis *x = fThnEff.GetAxis(0);
    
    Int_t bin[4];
    bin[0] = fThnEff->GetAxis(0)->FindBin(eta);
    bin[1] = fThnEff->GetAxis(1)->FindBin(pt);
    bin[2] = fThnEff->GetAxis(2)->FindBin(CentrOrMult);
    bin[3] = fThnEff->GetAxis(3)->FindBin(zVertex);
    
    
    // Int_t bin=fThnEff.FindBin(CentrOrMult,eta,pt);
    
    
    efficiency = fThnEff->GetBinContent(bin);
    
    if(efficiency == 0) efficiency = 1;
    
    return efficiency;
    
}

  
