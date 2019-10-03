#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDPmdTrack.h"
#include "AliESDVertex.h"
#include "AliAnalysisTaskPMDSim.h"
#include <TParticle.h>
#include <AliLog.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"

// AnalysisTask For PMD
// Authors: Sudipan De, Subhash Singha

ClassImp(AliAnalysisTaskPMDSim)

//________________________________________________________________________
AliAnalysisTaskPMDSim::AliAnalysisTaskPMDSim(const char *name) 
  : AliAnalysisTaskSE(name), 
  fESD(0), 
  fOutputList(0), 
  fHistTotEvent(0),
  fHistTotEventAfterPhySel(0),
  fHistTotEventAfterVtx(0),
  fVtxZ(0),
  fHistXYPre(0),
  fHistEtaPhM(0),
  fHistEtaPhM1(0),
  fHistEtaT(0),
  fMultMeasured(0),
  fMultMeasured1(0),
  fMultTrue(0),
  fMultCorr(0),
  fMultCorr1(0)
{
  for(Int_t i=0; i<10; i++){
    fHistMultMeasEtaBinA[i] = 0;
    fHistMultMeasEtaBinA1[i] = 0;
    fHistMultTrueEtaBinA[i] = 0;
    fHistMultCorrEtaBinA[i] = 0;
    fHistMultCorrEtaBinA1[i] = 0;
  }
  
  // Constructor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskPMDSim::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  fOutputList = new TList();
  
  Int_t kNbinsMult = 100;  Float_t XminMult = -0.5; 
  Float_t XmaxMult = 99.5;//mult
  Int_t kNbinsMultA = 50;  Float_t XminMultA = -0.5; 
  Float_t XmaxMultA = 49.5;//mult
  Int_t kNbinEtaPMDCov = 10; Float_t XminEtaPMDCov = 2.1; Float_t XmaxEtaPMDCov = 4.1;//etaPMD
  Int_t kNBinsEvent = 10; Float_t XminEvent = 0; Float_t XmaxEvent = 10;
  Int_t kNbinsXY = 200; Float_t XminXY = -100.0; Float_t XmaxXY  = 100.0;
  fHistTotEvent = new TH1F("TotEvent","TotEvent",
			   kNBinsEvent,XminEvent,XmaxEvent); 
  fOutputList->Add(fHistTotEvent);
  fHistTotEventAfterPhySel = new TH1F("TotEventAfterPhySel","TotEventAfterPhySel",
				      kNBinsEvent,XminEvent,XmaxEvent); 
  fOutputList->Add(fHistTotEventAfterPhySel);
  fHistTotEventAfterVtx = new TH1F("TotEventAfterVtx","TotEventAfterVtx",
				  kNBinsEvent,XminEvent,XmaxEvent); 
  fOutputList->Add(fHistTotEventAfterVtx);
  fVtxZ     = new TH1F("VtxZ","VtxZ",200,-50.0,50.0);
  fOutputList->Add(fVtxZ);
  fHistXYPre = new TH2F("XYPre","XYPre",kNbinsXY,XminXY,XmaxXY,
                        kNbinsXY,XminXY,XmaxXY);
  fOutputList->Add(fHistXYPre);
  fHistEtaPhM = new TH1F("fHistEtaPhM","fHistEtaPhM",kNbinEtaPMDCov,XminEtaPMDCov,XmaxEtaPMDCov);
  fOutputList->Add(fHistEtaPhM);
  fHistEtaPhM1 = new TH1F("fHistEtaPhM1","fHistEtaPhM1",kNbinEtaPMDCov,XminEtaPMDCov,XmaxEtaPMDCov);
  fOutputList->Add(fHistEtaPhM1);
  fHistEtaT = new TH1F("fHistEtaT","fHistEtaT",kNbinEtaPMDCov,XminEtaPMDCov,XmaxEtaPMDCov);
  fOutputList->Add(fHistEtaT);
  fMultMeasured = new TH1F("MultM","MultM",kNbinsMult,XminMult,XmaxMult);
  fOutputList->Add(fMultMeasured);
  fMultMeasured1 = new TH1F("MultM1","MultM1",kNbinsMult,XminMult,XmaxMult);
  fOutputList->Add(fMultMeasured1);
  fMultTrue = new TH1F("MultT","MultT",kNbinsMult,XminMult,XmaxMult);
  fOutputList->Add(fMultTrue);
  fMultCorr = new TH2F("MultCorr","MultCorr",kNbinsMult,XminMult,XmaxMult,
  		       kNbinsMult,XminMult,XmaxMult);
  fOutputList->Add(fMultCorr);
  fMultCorr1 = new TH2F("MultCorr1","MultCorr1",kNbinsMult,XminMult,XmaxMult,
  		       kNbinsMult,XminMult,XmaxMult);
  fOutputList->Add(fMultCorr1);
   
  
  Char_t nameT[256], nameM[256], nameCorr[256],nameM1[256], nameCorr1[256];
  for(Int_t i=0; i<10; i++){
    sprintf(nameM,"MultM_EtaBin%d",i+1);
    sprintf(nameM1,"MultM1_EtaBin%d",i+1);
    sprintf(nameT,"MultT_EtaBin%d",i+1);
    sprintf(nameCorr,"MultCorr_EtaBin%d",i+1);
    sprintf(nameCorr1,"MultCorr1_EtaBin%d",i+1);
    fHistMultMeasEtaBinA[i] = 
      new TH1F(nameM,nameM,kNbinsMultA,XminMultA,XmaxMultA);
    fHistMultMeasEtaBinA1[i] = 
      new TH1F(nameM1,nameM1,kNbinsMultA,XminMultA,XmaxMultA);
    fHistMultTrueEtaBinA[i] = 
      new TH1F(nameT,nameT,kNbinsMultA,XminMultA,XmaxMultA);
    fHistMultCorrEtaBinA[i] = 
      new TH2F(nameCorr,nameCorr,kNbinsMultA,XminMultA,XmaxMultA,
	       kNbinsMultA,XminMultA,XmaxMultA);
    fHistMultCorrEtaBinA1[i] = 
      new TH2F(nameCorr1,nameCorr1,kNbinsMultA,XminMultA,XmaxMultA,
	       kNbinsMultA,XminMultA,XmaxMultA);
    fOutputList->Add(fHistMultMeasEtaBinA[i]);
    fOutputList->Add(fHistMultMeasEtaBinA1[i]);
    fOutputList->Add(fHistMultTrueEtaBinA[i]);
    fOutputList->Add(fHistMultCorrEtaBinA[i]);
    fOutputList->Add(fHistMultCorrEtaBinA1[i]);
  }//i loop
  //fOutputList->Add(fHistPt);
}

//________________________________________________________________________
void AliAnalysisTaskPMDSim::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  Float_t MipCut1 = 432;//MPV=432=> 6*MPV=432
  Float_t MipCut2 = 648;//MPV=432=> 6*MPV=648 
  Float_t etacls=0., theta=0., rdist=0.;
  Int_t PhotonCls = 0;
  Int_t PhotonCls1 = 0;
  Int_t PhotonClsAEtaBin[10] = {0};//#of photon measured with diff Eta bin
  Int_t PhotonClsAEtaBin1[10] = {0};//#of photon measured with diff Eta bin
  Int_t PhotonTrueAEtaBin[10] = {0};//#of photon Incident with diff Eta bin
  // Post output data.
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  fHistTotEvent->Fill(5);
  //event selection
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  if (! isSelected) return;
  fHistTotEventAfterPhySel->Fill(5);
  //Vertex selection
  const AliESDVertex *vertex = fESD->GetPrimaryVertex();
  Float_t Vz = vertex->GetZ();    
  // Float_t Vx = vertex->GetX();    
  //Float_t Vy = vertex->GetY();
  Bool_t zVerStatus = vertex->GetStatus();
  if(zVerStatus){
    fVtxZ->Fill(Vz);
    if(TMath::Abs(Vz)<10){   
      fHistTotEventAfterVtx->Fill(5);   
      Int_t ptracks = fESD->GetNumberOfPmdTracks();
      for(Int_t kk=0;kk<ptracks;kk++){
	AliESDPmdTrack *pmdtr = fESD->GetPmdTrack(kk);
	Int_t   det   = pmdtr->GetDetector();
	Float_t clsX  = pmdtr->GetClusterX();
	Float_t clsY  = pmdtr->GetClusterY();
	Float_t clsZ  = pmdtr->GetClusterZ();
	clsZ -= Vz;
	Float_t ncell = pmdtr->GetClusterCells();
	Float_t adc   = pmdtr->GetClusterADC();
	//Float_t pid   = pmdtr->GetClusterPID();
	//Float_t isotag = pmdtr->GetClusterSigmaX();
	//Int_t trno = pmdtr->GetClusterTrackNo();
	//calculation of #eta
	rdist = TMath::Sqrt(clsX*clsX + clsY*clsY);
	if(clsZ!=0) theta = TMath::ATan2(rdist,clsZ);
	etacls   = -TMath::Log(TMath::Tan(0.5*theta));
	
	if(det==0 && adc>0)fHistXYPre->Fill(clsX,clsY);
	if(det==0 && adc>MipCut1 && ncell>2){
	  fHistEtaPhM->Fill(etacls);
	  if(etacls > 2.3 && etacls <= 3.9)  PhotonCls++;
	  if(etacls > 2.1 && etacls <= 2.3)  PhotonClsAEtaBin[0]++;
	  if(etacls > 2.3 && etacls <= 2.5)  PhotonClsAEtaBin[1]++;
	  if(etacls > 2.5 && etacls <= 2.7)  PhotonClsAEtaBin[2]++;
	  if(etacls > 2.7 && etacls <= 2.9)  PhotonClsAEtaBin[3]++;
	  if(etacls > 2.9 && etacls <= 3.1)  PhotonClsAEtaBin[4]++;
	  if(etacls > 3.1 && etacls <= 3.3)  PhotonClsAEtaBin[5]++;
	  if(etacls > 3.3 && etacls <= 3.5)  PhotonClsAEtaBin[6]++;
	  if(etacls > 3.5 && etacls <= 3.7)  PhotonClsAEtaBin[7]++;
	  if(etacls > 3.7 && etacls <= 3.9)  PhotonClsAEtaBin[8]++;
	  if(etacls > 3.9 && etacls <= 4.1)  PhotonClsAEtaBin[9]++;
	}//if MipCut1
	if( det==0 && adc>MipCut2 && ncell>2){
	  fHistEtaPhM1->Fill(etacls);
	  if(etacls > 2.3 && etacls <= 3.9)  PhotonCls1++;
	  if(etacls > 2.1 && etacls <= 2.3)  PhotonClsAEtaBin1[0]++;
	  if(etacls > 2.3 && etacls <= 2.5)  PhotonClsAEtaBin1[1]++;
	  if(etacls > 2.5 && etacls <= 2.7)  PhotonClsAEtaBin1[2]++;
	  if(etacls > 2.7 && etacls <= 2.9)  PhotonClsAEtaBin1[3]++;
	  if(etacls > 2.9 && etacls <= 3.1)  PhotonClsAEtaBin1[4]++;
	  if(etacls > 3.1 && etacls <= 3.3)  PhotonClsAEtaBin1[5]++;
	  if(etacls > 3.3 && etacls <= 3.5)  PhotonClsAEtaBin1[6]++;
	  if(etacls > 3.5 && etacls <= 3.7)  PhotonClsAEtaBin1[7]++;
	  if(etacls > 3.7 && etacls <= 3.9)  PhotonClsAEtaBin1[8]++;
	  if(etacls > 3.9 && etacls <= 4.1)  PhotonClsAEtaBin1[9]++;
	}//if MipCut2
      } //track loop 
      //reading MC info.
      AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>
	(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!eventHandler) {
	Printf("ERROR: Could not retrieve MC event handler");
	return;
      }
      AliMCEvent* mcEvent = eventHandler->MCEvent();
      if (!mcEvent) {
	Printf("ERROR: Could not retrieve MC event");
	return;
      }
      AliStack* stack = mcEvent->Stack();
      if (!stack)
	{
	  AliDebug(AliLog::kError, "Stack not available");
	  return;
	}
      if(stack){
	Int_t nPrim  = stack->GetNprimary();
	Int_t PhotonTrue = 0;
	for (Int_t iMc = 0; iMc < nPrim; ++iMc)
	  {
	    TParticle *MPart = stack->Particle(iMc);
	    Int_t mpart  = MPart->GetPdgCode();
	    Float_t eta   = MPart->Eta();
	    if(mpart == 22){  
	      fHistEtaT->Fill(eta);
	      if(eta > 2.3 && eta <= 3.9)PhotonTrue++;
	      if(eta > 2.1 && eta <= 2.3)PhotonTrueAEtaBin[0]++;
	      if(eta > 2.3 && eta <= 2.5)PhotonTrueAEtaBin[1]++;
	      if(eta > 2.5 && eta <= 2.7)PhotonTrueAEtaBin[2]++;
	      if(eta > 2.7 && eta <= 2.9)PhotonTrueAEtaBin[3]++;
	      if(eta > 2.9 && eta <= 3.1)PhotonTrueAEtaBin[4]++;
	      if(eta > 3.1 && eta <= 3.3)PhotonTrueAEtaBin[5]++;
	      if(eta > 3.3 && eta <= 3.5)PhotonTrueAEtaBin[6]++;
	      if(eta > 3.5 && eta <= 3.7)PhotonTrueAEtaBin[7]++;
	      if(eta > 3.7 && eta <= 3.9)PhotonTrueAEtaBin[8]++;
	      if(eta > 3.9 && eta <= 4.1)PhotonTrueAEtaBin[9]++;
	    }//mpart
	  }//track loop
	for(Int_t i=0; i<10; i++){
	  fHistMultMeasEtaBinA[i]->Fill(PhotonClsAEtaBin[i]);
	  fHistMultMeasEtaBinA1[i]->Fill(PhotonClsAEtaBin1[i]);
	  fHistMultTrueEtaBinA[i]->Fill(PhotonTrueAEtaBin[i]);
	  fHistMultCorrEtaBinA[i]->Fill(PhotonTrueAEtaBin[i],PhotonClsAEtaBin[i]);
	  fHistMultCorrEtaBinA1[i]->Fill(PhotonTrueAEtaBin[i],PhotonClsAEtaBin1[i]);
	}//i loop
	fMultMeasured->Fill(PhotonCls);
	fMultMeasured1->Fill(PhotonCls1);
	fMultTrue->Fill(PhotonTrue);
	fMultCorr->Fill(PhotonTrue,PhotonCls);
	fMultCorr1->Fill(PhotonTrue,PhotonCls1);
      }//if stack
    }//if Vz!= 0
  }//Vz loop
  PostData(1, fOutputList);
}//UserExec()     

//_______________________________________________________________________
void AliAnalysisTaskPMDSim::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}//Terminate
