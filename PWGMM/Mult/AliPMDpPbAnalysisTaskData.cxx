#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TArrayF.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDPmdTrack.h"
#include "AliHeader.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDVZERO.h"
#include "AliPMDpPbAnalysisTaskData.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliPPVsMultUtils.h"
#include "AliAnalysisUtils.h"
#include "AliVEvent.h"

//PMD pA data analysis
//Author: Abhi Modak (abhi.modak@cern.ch)
//Co-Authors: Sudipan De (Sudipan.De@cern.ch) and Sidharth Kumar Prasad (sidharth.kumar.prasad@cern.ch)
//Taken the help from the code: alisw/AliPhysics/PWGLF/FORWARD/photons/AliAnalysisTaskPMD.cxx

ClassImp(AliPMDpPbAnalysisTaskData)

//.......................................................//

AliPMDpPbAnalysisTaskData::AliPMDpPbAnalysisTaskData(const char *name)
: AliAnalysisTaskSE(name),
  fESD(0), 
  fEsdV0(0),
  fOutputList(0),
  fCentEstimator("V0A"),
  fTrigSel("kINT7"),
  fHistClsXYPre(0),
  fHistClsXYCpv(0),
  fHistVtxZ(0),
  fHistTotEvent(0),
  fHistEtaCut1(0),
  fHistEtaCut2(0),
  fHistNcellCut1(0),
  fHistNcellCut2(0),
  ntMeas1(0),
  ntMeas2(0),
  ntCent(0),
  ntCorr(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliPMDpPbAnalysisTaskData::UserCreateOutputObjects() {

  fOutputList = new TList();
  fOutputList->SetOwner();

  Int_t NbinsEvt = 10, XminEvt  = 0, XmaxEvt  = 0, NbinsClsXY = 200, NbinsVtxZ = 100, NbinsEta = 10, NbinsNcell = 50;
  Float_t XYmin = -100.0, XYmax = 100.0, XminEta = 2.1, XmaxEta = 4.1, XminNcell = 0.0, XmaxNcell = 50.0, VtxZmin = -50.0, VtxZmax = 50.0;
  
  fHistTotEvent = new TH1F("fHistTotEvent","fHistTotEvent",5,0,5);
  fHistTotEvent->GetXaxis()->SetBinLabel(1, "Processed events");
  fHistTotEvent->GetXaxis()->SetBinLabel(2, (fTrigSel+" trigger").Data());
  fHistTotEvent->GetXaxis()->SetBinLabel(3, "Pileup Multiple Vertices");                                                                      
  fHistTotEvent->GetXaxis()->SetBinLabel(4, "NContributors > 0");
  fHistTotEvent->GetXaxis()->SetBinLabel(5, "|VtxZ|<10:Selected for Analysis");
  fHistClsXYPre         = new TH2F("fHistClsXYPre","fHistClsXYPre",NbinsClsXY,XYmin,XYmax,NbinsClsXY,XYmin,XYmax);
  fHistClsXYCpv         = new TH2F("fHistClsXYCpv","fHistClsXYCpv",NbinsClsXY,XYmin,XYmax,NbinsClsXY,XYmin,XYmax);
  fHistVtxZ             = new TH1F("fHistVtxZ","fHistVtxZ",NbinsVtxZ,VtxZmin,VtxZmax);
  fHistEtaCut1          = new TH1F("EtaCut1","EtaCut1",NbinsEta,XminEta,XmaxEta);
  fHistEtaCut2          = new TH1F("EtaCut2","EtaCut2",NbinsEta,XminEta,XmaxEta);
  fHistNcellCut1        = new TH1F("NcellCut1","NcellCut1",NbinsNcell,XminNcell,XmaxNcell);
  fHistNcellCut2        = new TH1F("NcellCut2","NcellCut2",NbinsNcell,XminNcell,XmaxNcell);

  ntMeas1 = new TNtuple("ntMeas1","MeasPhotonCls1","PhotonCls1:PhotonCls1EtaBin1:PhotonCls1EtaBin2:PhotonCls1EtaBin3:PhotonCls1EtaBin4:PhotonCls1EtaBin5:PhotonCls1EtaBin6:PhotonCls1EtaBin7:PhotonCls1EtaBin8");
  ntMeas2 = new TNtuple("ntMeas2","MeasPhotonCls2","PhotonCls2:PhotonCls2EtaBin1:PhotonCls2EtaBin2:PhotonCls2EtaBin3:PhotonCls2EtaBin4:PhotonCls2EtaBin5:PhotonCls2EtaBin6:PhotonCls2EtaBin7:PhotonCls2EtaBin8");
  ntCent  = new TNtuple((fCentEstimator+"Centrality").Data(),"Centrality","nCentrality");
  ntCorr  = new TNtuple("ntCorr","Correlation","fV0mult:fV0Amult:fV0Cmult:GlobalTracks");

  fOutputList->Add(fHistTotEvent);
  fOutputList->Add(fHistClsXYPre);
  fOutputList->Add(fHistClsXYCpv);
  fOutputList->Add(fHistVtxZ);
  fOutputList->Add(fHistEtaCut1);
  fOutputList->Add(fHistEtaCut2);
  fOutputList->Add(fHistNcellCut1);
  fOutputList->Add(fHistNcellCut2);
  fOutputList->Add(ntMeas1);
  fOutputList->Add(ntMeas2);
  fOutputList->Add(ntCent);
  fOutputList->Add(ntCorr);

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliPMDpPbAnalysisTaskData::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  Int_t MipCut1 = 432, MipCut2 = 648;
  Float_t etacls, theta, rdist, phicls;

  Int_t PhClsEtaBinCut1[10], PhClsEtaBinCut2[10];

  for(Int_t i = 0; i < 10; i++){
    PhClsEtaBinCut1[i] = 0;
    PhClsEtaBinCut2[i] = 0;    
  }

  Int_t PhotonCls1 = 0, PhotonCls1EtaBin1 = 0, PhotonCls1EtaBin2 = 0, PhotonCls1EtaBin3 = 0, PhotonCls1EtaBin4 = 0, PhotonCls1EtaBin5 = 0, PhotonCls1EtaBin6 = 0, PhotonCls1EtaBin7 = 0, PhotonCls1EtaBin8 = 0;

  Int_t PhotonCls2 = 0, PhotonCls2EtaBin1 = 0, PhotonCls2EtaBin2 = 0, PhotonCls2EtaBin3 = 0, PhotonCls2EtaBin4 = 0, PhotonCls2EtaBin5 = 0, PhotonCls2EtaBin6 = 0, PhotonCls2EtaBin7 = 0, PhotonCls2EtaBin8 = 0;

  Float_t fV0mult = 0.0, fV0Amult = 0.0, fV0Cmult = 0.0, GlobalTracks = 0.0, nCentrality = 300.0, VtxZ;

  AliMultSelection *MultSelection = 0x0; 

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  fHistTotEvent->Fill(0.5);
  //Printf("Event accepted");

  //event selection
  Bool_t isSelected = 0;
  if(fTrigSel == "kINT7")
    isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  if(fTrigSel == "kINT1")
    isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT1);
  
  if (!isSelected) return;

  fHistTotEvent->Fill(1.5);

  AliAnalysisUtils util;                                                                                                                      
  util.SetMinPlpContribMV(5);
  util.SetMaxPlpChi2MV(5);
  util.SetMinWDistMV(15);
  util.SetCheckPlpFromDifferentBCMV(kFALSE);
  Bool_t IsPileUpMV = util.IsPileUpMV(fESD);                                                                                                
  if(IsPileUpMV)return;                                                                                                                       
  fHistTotEvent->Fill(2.5);

  const AliESDVertex *vertex = fESD->GetPrimaryVertex();
  VtxZ = vertex->GetZ();
  if(vertex){
    if(vertex->GetNContributors() > 0){
      fHistTotEvent->Fill(3.5);
      //if(vertex->GetZRes() != 0){
      if(TMath::Abs(VtxZ) < 10){
	fHistTotEvent->Fill(4.5);
	fHistVtxZ->Fill(VtxZ);
	
	MultSelection = (AliMultSelection *) fESD->FindListObject("MultSelection");
	if( !MultSelection) {
	  //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
	  AliWarning("AliMultSelection object not found!");
	}//if loop
	else{
	  nCentrality = MultSelection->GetMultiplicityPercentile(fCentEstimator);
	}//else loop

	fEsdV0   = fESD->GetVZEROData();
	fV0Amult = fEsdV0->GetMTotV0A(); //returns total multiplicity in V0A 
	fV0Cmult = fEsdV0->GetMTotV0C(); //returns total multiplicity in V0C
	fV0mult  = fV0Amult + fV0Cmult;   //returns total multiplicity in V0A+V0C
	GlobalTracks  = fESD->GetNumberOfTracks();  
	
	Int_t PMDTrackks = fESD->GetNumberOfPmdTracks();
	for(Int_t trk = 0; trk < PMDTrackks; trk++){
	  AliESDPmdTrack *pmdtr = fESD->GetPmdTrack(trk);
	  Int_t   det           = pmdtr->GetDetector();
	  Float_t clsX          = pmdtr->GetClusterX();
	  Float_t clsY          = pmdtr->GetClusterY();
	  Float_t clsZ          = pmdtr->GetClusterZ();
	  clsZ -= VtxZ;
	  Float_t ncell         = pmdtr->GetClusterCells();
	  Float_t adc           = pmdtr->GetClusterADC();
	  Float_t pid           = pmdtr->GetClusterPID();
	  Float_t isotag        = pmdtr->GetClusterSigmaX();
	  Int_t trno            = pmdtr->GetClusterTrackNo();
	  Float_t trpid         = pmdtr->GetClusterTrackPid();
	  Int_t   smn           = pmdtr->GetSmn();
	  
	  rdist                 = TMath::Sqrt(clsX*clsX + clsY*clsY);
	  if(clsZ!=0) theta     = TMath::ATan2(rdist,clsZ);
	  etacls                = -TMath::Log(TMath::Tan(0.5*theta));
	  
	  if(det == 0 && adc > 0)fHistClsXYPre->Fill(clsX,clsY);
	  if(det == 1 && adc > 0)fHistClsXYCpv->Fill(clsX,clsY);
	  
	  if(det == 0 && adc > MipCut1 && ncell > 2){
	    if(etacls > 2.3 && etacls <= 3.9){
	      fHistEtaCut1->Fill(etacls);
	      fHistNcellCut1->Fill(ncell);
	      PhotonCls1++;
	    }//etacls loop
	    for(Int_t i=0;i<10;i++){
	      Float_t etamin = 2.1 + i*0.2;
	      Float_t etamax = etamin + 0.2;
	      if(etacls > etamin && etacls <= etamax)PhClsEtaBinCut1[i]++;
	    }//for loop
	  }//MipCut1 loop
	  
	  if(det == 0 && adc > MipCut2 && ncell > 2){
	    if(etacls > 2.3 && etacls <= 3.9){
	      fHistEtaCut2->Fill(etacls);
	      fHistNcellCut2->Fill(ncell);
	      PhotonCls2++;
	    }//etacls loop
	    for(Int_t i=0;i<10;i++){
	      Float_t etamin = 2.1 + i*0.2;
	      Float_t etamax = etamin + 0.2;
	      if(etacls > etamin && etacls <= etamax)PhClsEtaBinCut2[i]++;
	    }//for loop
	  }//MipCut2 loop
	}//pmd track loop
	
	PhotonCls1EtaBin1 = PhClsEtaBinCut1[1];
	PhotonCls1EtaBin2 = PhClsEtaBinCut1[2];
	PhotonCls1EtaBin3 = PhClsEtaBinCut1[3];
	PhotonCls1EtaBin4 = PhClsEtaBinCut1[4];
	PhotonCls1EtaBin5 = PhClsEtaBinCut1[5];
	PhotonCls1EtaBin6 = PhClsEtaBinCut1[6];
	PhotonCls1EtaBin7 = PhClsEtaBinCut1[7];
	PhotonCls1EtaBin8 = PhClsEtaBinCut1[8];
	
	PhotonCls2EtaBin1 = PhClsEtaBinCut2[1];
	PhotonCls2EtaBin2 = PhClsEtaBinCut2[2];
	PhotonCls2EtaBin3 = PhClsEtaBinCut2[3];
	PhotonCls2EtaBin4 = PhClsEtaBinCut2[4];
	PhotonCls2EtaBin5 = PhClsEtaBinCut2[5];
	PhotonCls2EtaBin6 = PhClsEtaBinCut2[6];
	PhotonCls2EtaBin7 = PhClsEtaBinCut2[7];
	PhotonCls2EtaBin8 = PhClsEtaBinCut2[8];
	
	ntMeas1->Fill(PhotonCls1,PhotonCls1EtaBin1,PhotonCls1EtaBin2,PhotonCls1EtaBin3,PhotonCls1EtaBin4,PhotonCls1EtaBin5,PhotonCls1EtaBin6,PhotonCls1EtaBin7,PhotonCls1EtaBin8);
	ntMeas2->Fill(PhotonCls2,PhotonCls2EtaBin1,PhotonCls2EtaBin2,PhotonCls2EtaBin3,PhotonCls2EtaBin4,PhotonCls2EtaBin5,PhotonCls2EtaBin6,PhotonCls2EtaBin7,PhotonCls2EtaBin8);
	ntCorr->Fill(fV0mult,fV0Amult,fV0Cmult,GlobalTracks);
	ntCent->Fill(nCentrality);
      }//VtxZ<10
      //}//GetZRes
    }//Contributors
  }//vertex
  
}

//________________________________________________________________________
void AliPMDpPbAnalysisTaskData::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}


