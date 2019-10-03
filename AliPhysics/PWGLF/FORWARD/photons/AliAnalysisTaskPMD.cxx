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
#include "AliESDVZERO.h"
#include "AliAnalysisTaskPMD.h"

// AnalysisTask For PMD
// Authors: Sudipan De, Subhash Singha

ClassImp(AliAnalysisTaskPMD)

//________________________________________________________________________
AliAnalysisTaskPMD::AliAnalysisTaskPMD(const char *name) 
: AliAnalysisTaskSE(name), 
  fESD(0), 
  fOutputList(0), 
  fHistTotEvent(0),
  fHistTotEventAfterPhySel(0),
  fHistTotEventAfterVtx(0),
  fHistVtxZ(0),
  fHistXYPre(0),
  fHistEta(0),
  fHistEta1(0),
  fHistMultMeas(0),
  fHistMultMeas1(0)
{
  for(Int_t i=0; i<10; i++){
    fHistMultMeasEtaBinA[i] = 0;
    fHistMultMeasEtaBinA1[i] = 0;
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
void AliAnalysisTaskPMD::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  fOutputList = new TList();
  Int_t kNbinsMultA = 50;  Float_t XminMultA = -0.5; Float_t XmaxMultA = 49.5;
  Int_t kNbinsXY = 200; Float_t XminXY = -100.0; Float_t XmaxXY  = 100.0;
  Int_t kNBinsEvent = 10; Float_t XminEvent = 0; Float_t XmaxEvent = 10;
  Int_t kNBinsEta = 10; Float_t XminEta = 2.1; Float_t XmaxEta = 4.1;
 
  fHistTotEvent = new TH1F("TotEvent","TotEvent",
			   kNBinsEvent,XminEvent,XmaxEvent); 
  fOutputList->Add(fHistTotEvent);
  fHistTotEventAfterPhySel = new TH1F("TotEventAfterPhySel","TotEventAfterPhySel",
				      kNBinsEvent,XminEvent,XmaxEvent); 
  fOutputList->Add(fHistTotEventAfterPhySel);
  fHistTotEventAfterVtx = new TH1F("TotEventAfterVtx","TotEventAfterVtx",
				   kNBinsEvent,XminEvent,XmaxEvent); 
  fOutputList->Add(fHistTotEventAfterVtx);
  fHistVtxZ = new TH1F("VtxZ","VtxZ",100,-50,50);
  fOutputList->Add(fHistVtxZ);
  fHistXYPre = new TH2F("XYPre","XYPre",kNbinsXY,XminXY,XmaxXY,
			kNbinsXY,XminXY,XmaxXY);
  fOutputList->Add(fHistXYPre);
  fHistEta = new TH1F ("Eta","Eta",kNBinsEta,XminEta,XmaxEta);
  fOutputList->Add(fHistEta);
  fHistEta1 = new TH1F ("Eta1","Eta1",kNBinsEta,XminEta,XmaxEta);
  fOutputList->Add(fHistEta1);
  fHistMultMeas = new TH1F("MultM","MultM",100,-0.5,99.5);
  fOutputList->Add(fHistMultMeas);
  fHistMultMeas1 = new TH1F("MultM1","MultM1",100,-0.5,99.5);
  fOutputList->Add(fHistMultMeas1);
  
  Char_t nameM[256], nameM1[256];
  for(Int_t i=0; i<10; i++){
    sprintf(nameM,"MultM_EtaBin%d",i+1);
    fHistMultMeasEtaBinA[i] =
      new TH1F(nameM,nameM,kNbinsMultA,
	       XminMultA,XmaxMultA);
    fOutputList->Add(fHistMultMeasEtaBinA[i]);
    sprintf(nameM1,"MultM1_EtaBin%d",i+1);
    fHistMultMeasEtaBinA1[i] =
      new TH1F(nameM1,nameM1,kNbinsMultA,
	       XminMultA,XmaxMultA);
    fOutputList->Add(fHistMultMeasEtaBinA1[i]);
  }//i loop
  
}

//________________________________________________________________________
void AliAnalysisTaskPMD::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  Float_t MipCut1 = 432;//MPV=72=> 6*MPV=432
  Float_t MipCut2 = 648;//MPV=72=> 9*MPV=648
  Float_t etacls=0., theta=0., rdist=0.;
  Int_t PhotonCls = 0;//# of photon measured within 2.3 to 3.9 #eta
  Int_t PhotonCls1 = 0;//# of photon measured within 2.3 to 3.9 #eta
  Int_t PhotonClsAEtaBin[10] = {0};//# of photon measured with diff Eta bin
  Int_t PhotonClsAEtaBin1[10] = {0};//# of photon measured with diff Eta bin
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
  //printf("There are %d tracks in this event\n", fESD->GetNumberOfPmdTracks());
  //vertex selection
  const AliESDVertex *vertex = fESD->GetPrimaryVertex();
  Float_t Vz = vertex->GetZ();    
  Bool_t zVerStatus = vertex->GetStatus();
  if(zVerStatus){
    fHistVtxZ->Fill(Vz);
    if(TMath::Abs(Vz)<10){
      fHistTotEventAfterVtx->Fill(5);
      Int_t ptracks = fESD->GetNumberOfPmdTracks();
      //PMDtrack Loop
      for(Int_t kk=0;kk<ptracks;kk++){
	AliESDPmdTrack *pmdtr = fESD->GetPmdTrack(kk);
	Int_t   det   = pmdtr->GetDetector();
	Float_t clsX  = pmdtr->GetClusterX();
	Float_t clsY  = pmdtr->GetClusterY();
	Float_t clsZ  = pmdtr->GetClusterZ();
	clsZ -= Vz;
	Float_t ncell = pmdtr->GetClusterCells();
	Float_t adc   = pmdtr->GetClusterADC();
	//Int_t smn = pmdtr->GetSmn();
	//Calculation of #eta
	rdist = TMath::Sqrt(clsX*clsX + clsY*clsY);
	if(clsZ!=0) theta = TMath::ATan2(rdist,clsZ);
	etacls  = -TMath::Log(TMath::Tan(0.5*theta));
	
	if(det==0 && adc>72)fHistXYPre->Fill(clsX,clsY);
	if(det==0 && adc>MipCut1 && ncell>2){
	  if(etacls > 2.1 && etacls<= 2.3)PhotonClsAEtaBin[0]++;
	  if(etacls > 2.3 && etacls<= 2.5)PhotonClsAEtaBin[1]++;
	  if(etacls > 2.5 && etacls<= 2.7)PhotonClsAEtaBin[2]++;
	  if(etacls > 2.7 && etacls<= 2.9)PhotonClsAEtaBin[3]++;
	  if(etacls > 2.9 && etacls<= 3.1)PhotonClsAEtaBin[4]++;
	  if(etacls > 3.1 && etacls<= 3.3)PhotonClsAEtaBin[5]++;
	  if(etacls > 3.3 && etacls<= 3.5)PhotonClsAEtaBin[6]++;
	  if(etacls > 3.5 && etacls<= 3.7)PhotonClsAEtaBin[7]++;
	  if(etacls > 3.7 && etacls<= 3.9)PhotonClsAEtaBin[8]++;
	  if(etacls > 3.9 && etacls<= 4.1)PhotonClsAEtaBin[9]++;
	  if(etacls > 2.3 && etacls<= 3.9)PhotonCls++;
	  fHistEta->Fill(etacls);
	}//if MipCut1
	if(det==0 && adc>MipCut2 && ncell>2){
	  if(etacls > 2.1 && etacls<= 2.3)PhotonClsAEtaBin1[0]++;
	  if(etacls > 2.3 && etacls<= 2.5)PhotonClsAEtaBin1[1]++;
	  if(etacls > 2.5 && etacls<= 2.7)PhotonClsAEtaBin1[2]++;
	  if(etacls > 2.7 && etacls<= 2.9)PhotonClsAEtaBin1[3]++;
	  if(etacls > 2.9 && etacls<= 3.1)PhotonClsAEtaBin1[4]++;
	  if(etacls > 3.1 && etacls<= 3.3)PhotonClsAEtaBin1[5]++;
	  if(etacls > 3.3 && etacls<= 3.5)PhotonClsAEtaBin1[6]++;
	  if(etacls > 3.5 && etacls<= 3.7)PhotonClsAEtaBin1[7]++;
	  if(etacls > 3.7 && etacls<= 3.9)PhotonClsAEtaBin1[8]++;
	  if(etacls > 3.9 && etacls<= 4.1)PhotonClsAEtaBin1[9]++;
	  if(etacls > 2.3 && etacls<= 3.9)PhotonCls1++;
	  fHistEta1->Fill(etacls);
	}//if Mipcut2
      } //PMDtrack loop 
      for(Int_t i=0; i<10; i++){
	fHistMultMeasEtaBinA[i]->Fill(PhotonClsAEtaBin[i]);
	fHistMultMeasEtaBinA1[i]->Fill(PhotonClsAEtaBin1[i]);
      }//i loop
      fHistMultMeas->Fill(PhotonCls);
      fHistMultMeas1->Fill(PhotonCls1);
    }// Vz loop 
  }// Vzcut loop
  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskPMD::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}
