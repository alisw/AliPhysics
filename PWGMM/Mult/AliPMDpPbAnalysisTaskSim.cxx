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
#include "AliPMDpPbAnalysisTaskSim.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliGenEventHeader.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisUtils.h"
#include "AliMCEvent.h"

//PMD pA data analysis
//Author: Abhi Modak (abhi.modak@cern.ch)
//Co-Authors: Sudipan De (Sudipan.De@cern.ch) and Sidharth Kumar Prasad (sidharth.kumar.prasad@cern.ch)
//Take the help from the code: alisw/AliPhysics/PWGLF/FORWARD/photons/AliAnalysisTaskPMD.cxx

ClassImp(AliPMDpPbAnalysisTaskSim)

//.......................................................//

AliPMDpPbAnalysisTaskSim::AliPMDpPbAnalysisTaskSim(const char *name)
: AliAnalysisTaskSE(name),
  fESD(0), 
  fEsdV0(0),
  fOutputList(0),
  fCentEstimator("V0A"),
  fTrigSel("kINT7"),
  fHistClsXYPre(0),
  fHistClsXYCpv(0),
  fHistADCPre(0),
  fHistADCCpv(0),
  fHistVtxZ(0),
  fHistTotEvent(0),
  fHistEtaCut1(0),
  fHistEtaCut2(0),
  fHistNcellCut1(0),
  fHistNcellCut2(0),
  fHistEtaChTrue(0),
  fHistEtaPhTrue(0),
  ntMeas1(0),
  ntMeas2(0),
  ntCent(0),
  ntCorr(0),
  ntTrue(0),
  ntDet1(0),
  ntDet2(0)
{

  for(Int_t i=0; i<7; i++){
    fHistEtaChCentClass[i] = 0;
    fHistEtaGammaCentClass[i] = 0;
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
void AliPMDpPbAnalysisTaskSim::UserCreateOutputObjects() {

  fOutputList = new TList();
  fOutputList->SetOwner();

  Int_t NbinsEvt = 10, XminEvt  = 0, XmaxEvt  = 10, NbinsClsXY = 200, NbinsVtxZ = 100, NbinsEta = 10, NbinsNcell = 50, NbinsADC = 2500;
  Float_t XYmin = -100.0, XYmax = 100.0, XminEta = 2.1, XmaxEta = 4.1, XminNcell = 0.0, XmaxNcell = 50.0, VtxZmin = -50.0, VtxZmax = 50.0, XminADC = 0.0, XmaxADC = 5000.0;

  fHistTotEvent = new TH1F("fHistTotEvent","fHistTotEvent",6,0,6);
  fHistTotEvent->GetXaxis()->SetBinLabel(1, "Processed events");
  fHistTotEvent->GetXaxis()->SetBinLabel(2, (fTrigSel+" trigger").Data());
  fHistTotEvent->GetXaxis()->SetBinLabel(3, "Pileup Multiple Vertices");
  fHistTotEvent->GetXaxis()->SetBinLabel(4, "NContributors > 0");
  fHistTotEvent->GetXaxis()->SetBinLabel(5, "GetZRes != 0");
  fHistTotEvent->GetXaxis()->SetBinLabel(6, "|VtxZ|<10:Selected for Analysis");
  fHistClsXYPre         = new TH2F("fHistClsXYPre","fHistClsXYPre",NbinsClsXY,XYmin,XYmax,NbinsClsXY,XYmin,XYmax);
  fHistClsXYCpv         = new TH2F("fHistClsXYCpv","fHistClsXYCpv",NbinsClsXY,XYmin,XYmax,NbinsClsXY,XYmin,XYmax);
  fHistADCPre           = new TH1F("fHistADCPre","fHistADCPre",NbinsADC,XminADC,XmaxADC);
  fHistADCCpv           = new TH1F("fHistADCCpv","fHistADCCpv",NbinsADC,XminADC,XmaxADC);
  fHistVtxZ             = new TH1F("fHistVtxZ","fHistVtxZ",NbinsVtxZ,VtxZmin,VtxZmax);
  fHistEtaCut1          = new TH1F("EtaCut1","EtaCut1",NbinsEta,XminEta,XmaxEta);
  fHistEtaCut2          = new TH1F("EtaCut2","EtaCut2",NbinsEta,XminEta,XmaxEta);
  fHistNcellCut1        = new TH1F("NcellCut1","NcellCut1",NbinsNcell,XminNcell,XmaxNcell);
  fHistNcellCut2        = new TH1F("NcellCut2","NcellCut2",NbinsNcell,XminNcell,XmaxNcell);
  fHistEtaChTrue        = new TH1F("fHistEtaChTrue","fHistEtaChTrue",200,-10.0,10.0);
  fHistEtaPhTrue        = new TH1F("fHistEtaPhTrue","fHistEtaPhTrue",200,-10.0,10.0);

  ntMeas1 = new TNtuple("ntMeas1","MeasPhotonCls1","PhotonCls1:PhotonCls1EtaBin1:PhotonCls1EtaBin2:PhotonCls1EtaBin3:PhotonCls1EtaBin4:PhotonCls1EtaBin5:PhotonCls1EtaBin6:PhotonCls1EtaBin7:PhotonCls1EtaBin8");
  ntMeas2 = new TNtuple("ntMeas2","MeasPhotonCls2","PhotonCls2:PhotonCls2EtaBin1:PhotonCls2EtaBin2:PhotonCls2EtaBin3:PhotonCls2EtaBin4:PhotonCls2EtaBin5:PhotonCls2EtaBin6:PhotonCls2EtaBin7:PhotonCls2EtaBin8");
  ntCent  = new TNtuple((fCentEstimator+"Centrality").Data(),"Centrality","nCentrality");
  ntCorr  = new TNtuple("ntCorr","Correlation","fV0mult:fV0Amult:fV0Cmult:GlobalTracks");
  ntTrue  = new TNtuple("ntTrue","TruePhoton","PhotonTrue:PhotonTrueEtaBin1:PhotonTrueEtaBin2:PhotonTrueEtaBin3:PhotonTrueEtaBin4:PhotonTrueEtaBin5:PhotonTrueEtaBin6:PhotonTrueEtaBin7:PhotonTrueEtaBin8");
  ntDet1  = new TNtuple("ntDet1","DetectedPhoton","PhotonDet1:PhotonDet1EtaBin1:PhotonDet1EtaBin2:PhotonDet1EtaBin3:PhotonDet1EtaBin4:PhotonDet1EtaBin5:PhotonDet1EtaBin6:PhotonDet1EtaBin7:PhotonDet1EtaBin8");
  ntDet2  = new TNtuple("ntDet2","DetectedPhoton","PhotonDet2:PhotonDet2EtaBin1:PhotonDet2EtaBin2:PhotonDet2EtaBin3:PhotonDet2EtaBin4:PhotonDet2EtaBin5:PhotonDet2EtaBin6:PhotonDet2EtaBin7:PhotonDet2EtaBin8");

  fOutputList->Add(fHistTotEvent);
  fOutputList->Add(fHistClsXYPre);
  fOutputList->Add(fHistClsXYCpv);
  fOutputList->Add(fHistADCPre);
  fOutputList->Add(fHistADCCpv);
  fOutputList->Add(fHistVtxZ);
  fOutputList->Add(fHistEtaCut1);
  fOutputList->Add(fHistEtaCut2);
  fOutputList->Add(fHistNcellCut1);
  fOutputList->Add(fHistNcellCut2);
  fOutputList->Add(ntMeas1);
  fOutputList->Add(ntMeas2);
  fOutputList->Add(ntCent);
  fOutputList->Add(ntCorr);
  fOutputList->Add(ntTrue);
  fOutputList->Add(ntDet1);
  fOutputList->Add(ntDet2);
  fOutputList->Add(fHistEtaChTrue);
  fOutputList->Add(fHistEtaPhTrue);

  Char_t nameCh[256], nameGamma[256];
  for(Int_t i=0; i<7; i++){
    sprintf(nameCh,"EtaDistChCentClass%d",i+1);
    sprintf(nameGamma,"EtaDistGammaCentClass%d",i+1);
    fHistEtaChCentClass[i] = new TH1F(nameCh,nameCh,200,-10.0,10.0);
    fHistEtaGammaCentClass[i] = new TH1F(nameGamma,nameGamma,200,-10.0,10.0);
    fOutputList->Add(fHistEtaChCentClass[i]);
    fOutputList->Add(fHistEtaGammaCentClass[i]);
  }

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliPMDpPbAnalysisTaskSim::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  Int_t MipCut1 = 432, MipCut2 = 648;
  const Int_t maxTrackNo = 100000;
  Float_t etacls, theta, rdist, phicls;

  Int_t PhClsEtaBinCut1[10], PhClsEtaBinCut2[10], PhDetEtaBinCut1[10], PhDetEtaBinCut2[10], PhTrueEtaBin[10];

  for(Int_t i = 0; i < 10; i++){
    PhClsEtaBinCut1[i] = 0;
    PhClsEtaBinCut2[i] = 0;
    PhDetEtaBinCut1[i] = 0;
    PhDetEtaBinCut2[i] = 0;
    PhTrueEtaBin[i]    = 0;
  }

  Int_t PhotonCls1 = 0, PhotonCls1EtaBin1 = 0, PhotonCls1EtaBin2 = 0, PhotonCls1EtaBin3 = 0, PhotonCls1EtaBin4 = 0, PhotonCls1EtaBin5 = 0, PhotonCls1EtaBin6 = 0, PhotonCls1EtaBin7 = 0, PhotonCls1EtaBin8 = 0;

  Int_t PhotonCls2 = 0, PhotonCls2EtaBin1 = 0, PhotonCls2EtaBin2 = 0, PhotonCls2EtaBin3 = 0, PhotonCls2EtaBin4 = 0, PhotonCls2EtaBin5 = 0, PhotonCls2EtaBin6 = 0, PhotonCls2EtaBin7 = 0, PhotonCls2EtaBin8 = 0;

  Int_t PhotonTrue = 0, PhotonTrueEtaBin1 = 0, PhotonTrueEtaBin2 = 0, PhotonTrueEtaBin3 = 0, PhotonTrueEtaBin4 = 0, PhotonTrueEtaBin5 = 0, PhotonTrueEtaBin6 = 0, PhotonTrueEtaBin7 = 0, PhotonTrueEtaBin8 = 0;

  Int_t PhotonDet1 = 0, PhotonDet1EtaBin1 = 0, PhotonDet1EtaBin2 = 0, PhotonDet1EtaBin3 = 0, PhotonDet1EtaBin4 = 0, PhotonDet1EtaBin5 = 0, PhotonDet1EtaBin6 = 0, PhotonDet1EtaBin7 = 0, PhotonDet1EtaBin8 = 0;

  Int_t PhotonDet2 = 0, PhotonDet2EtaBin1 = 0, PhotonDet2EtaBin2 = 0, PhotonDet2EtaBin3 = 0, PhotonDet2EtaBin4 = 0, PhotonDet2EtaBin5 = 0, PhotonDet2EtaBin6 = 0, PhotonDet2EtaBin7 = 0, PhotonDet2EtaBin8 = 0;

  Float_t fV0mult = 0.0, fV0Amult = 0.0, fV0Cmult = 0.0, GlobalTracks = 0.0;

  Float_t nCentrality = 300.0;
  Float_t CentValue[] = {0., 5., 10., 20., 40., 60., 80., 100.};
  AliMultSelection *MultSelection = 0x0;

  Int_t TrStatusEtabinCut1[maxTrackNo], TrStatusEtabinCut2[maxTrackNo], TrStatusCut1[maxTrackNo], TrStatusCut2[maxTrackNo];

  for(Int_t tr = 0; tr < maxTrackNo; tr++){
    TrStatusEtabinCut1[tr] = -1;
    TrStatusEtabinCut2[tr] = -1;
    TrStatusCut1[tr]       = -1;
    TrStatusCut2[tr]       = -1;
  }

  Float_t VtxZ, VtxY, VtxX;

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  fHistTotEvent->Fill(0.5);

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
  VtxX = vertex->GetX();
  VtxY = vertex->GetY();
  
  if(vertex){
    if(vertex->GetNContributors() > 0){
      fHistTotEvent->Fill(3.5);
      //if(vertex->GetZRes() != 0){
      fHistTotEvent->Fill(4.5);
      if(TMath::Abs(VtxZ) < 10){
	fHistTotEvent->Fill(5.5);
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
	  
	  rdist                 = TMath::Sqrt(clsX*clsX + clsY*clsY);
	  if(clsZ!=0) theta     = TMath::ATan2(rdist,clsZ);
	  etacls                = -TMath::Log(TMath::Tan(0.5*theta));
	  
	  if(det == 0 && adc > 0)fHistClsXYPre->Fill(clsX,clsY);
	  if(det == 1 && adc > 0)fHistClsXYCpv->Fill(clsX,clsY);
	  
	  if(det == 0 && ncell == 1)fHistADCPre->Fill(adc);
	  if(det == 1 && ncell == 1)fHistADCCpv->Fill(adc);
	  
	  if(det == 0 && adc > MipCut1 && ncell > 2){
	    if(etacls > 2.3 && etacls <= 3.9){
	      fHistEtaCut1->Fill(etacls);
	      fHistNcellCut1->Fill(ncell);
	      PhotonCls1++;
	      if(trpid == 22){
		if(TrStatusCut1[trno] == -1){
		  PhotonDet1++;
		  TrStatusCut1[trno] = 1;
		}//trstatus
	      }//trpid
	    }//etacls loop
	    for(Int_t i=0;i<10;i++){
	      Float_t etamin = 2.1 + i*0.2;
	      Float_t etamax = etamin + 0.2;
	      if(etacls > etamin && etacls <= etamax){
		PhClsEtaBinCut1[i]++;
		if(trpid == 22){
		  if(TrStatusEtabinCut1[trno] == -1){
		    PhDetEtaBinCut1[i]++;
		    TrStatusEtabinCut1[trno] = 1;
		  }//trstatus
		}//trpid
	      }//etacls
	    }//for loop
	  }//MipCut1 loop
	  
	  if(det == 0 && adc > MipCut2 && ncell > 2){
	    if(etacls > 2.3 && etacls <= 3.9){
	      fHistEtaCut2->Fill(etacls);
	      fHistNcellCut2->Fill(ncell);
	      PhotonCls2++;
	      if(trpid == 22){
		if(TrStatusCut2[trno] == -1){
		  PhotonDet2++;
		  TrStatusCut2[trno] = 1;
		}//trstatus
	      }//trpid
	    }//etacls loop
	    for(Int_t i=0;i<10;i++){
	      Float_t etamin = 2.1 + i*0.2;
	      Float_t etamax = etamin + 0.2;
	      if(etacls > etamin && etacls <= etamax){
		PhClsEtaBinCut2[i]++;
		if(trpid == 22){
		  if(TrStatusEtabinCut2[trno] == -1){
		    PhDetEtaBinCut2[i]++;
		    TrStatusEtabinCut2[trno] = 1;
		  }//trstatus
		}//trpid
	      }//etacls
	    }//for loop
	  }//MipCut2 loop
	}//pmd track loop
	
	//True level information	    
	AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>
	  (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
	
	if (!eventHandler) {
	  Printf("ERROR: Could not retrieve MC event handler");
	  return;
	}//eventHandler
	
	AliMCEvent* mcEvent = eventHandler->MCEvent();
	if (!mcEvent) {
	  Printf("ERROR: Could not retrieve MC event");
	  return;
	}//mcEvent
	
	for(Int_t mcTrack = 0; mcTrack < mcEvent->GetNumberOfTracks(); mcTrack++) {
	  AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(mcTrack);
	  if (!track) {
	    Printf("ERROR: Could not receive track %d", mcTrack);
	    continue;
	  }//!track
	  if(!mcEvent->IsPhysicalPrimary(mcTrack)) continue;
	  
	  Int_t mpart  = track->PdgCode();
	  Float_t eta  = track->Eta();
	  Double_t pt  = track->Pt();
	  Int_t ch     = track->Charge();
	  Float_t phi  = track->Phi();

	  if(ch!=0){
	    for(Int_t i=0; i<7; i++){
	      if(nCentrality>CentValue[i] && nCentrality<=CentValue[i+1]) fHistEtaChCentClass[i]->Fill(eta);
	    }
	    fHistEtaChTrue->Fill(eta);
	  }
	  if(mpart == 22) {
	    fHistEtaPhTrue->Fill(eta);
	    for(Int_t i=0; i<7; i++){
	      if(nCentrality>CentValue[i] && nCentrality<=CentValue[i+1]) fHistEtaGammaCentClass[i]->Fill(eta);
	    }
	    if(eta > 2.3 && eta <= 3.9) PhotonTrue++;
	    for(Int_t i=0;i<10;i++) {
	      Float_t etamin = 2.1 + i*0.2;
	      Float_t etamax = etamin + 0.2;
	      if(eta > etamin && eta <= etamax){
		PhTrueEtaBin[i]++;
	      }//eta loop
	    }//for loop
	  }//mpart loop	       
	}//mcTrack loop
	
	PhotonTrueEtaBin1 = PhTrueEtaBin[1];
	PhotonTrueEtaBin2 = PhTrueEtaBin[2];
	PhotonTrueEtaBin3 = PhTrueEtaBin[3];
	PhotonTrueEtaBin4 = PhTrueEtaBin[4];
	PhotonTrueEtaBin5 = PhTrueEtaBin[5];
	PhotonTrueEtaBin6 = PhTrueEtaBin[6];
	PhotonTrueEtaBin7 = PhTrueEtaBin[7];
	PhotonTrueEtaBin8 = PhTrueEtaBin[8];
	
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
	
	PhotonDet1EtaBin1 = PhDetEtaBinCut1[1];
	PhotonDet1EtaBin2 = PhDetEtaBinCut1[2];
	PhotonDet1EtaBin3 = PhDetEtaBinCut1[3];
	PhotonDet1EtaBin4 = PhDetEtaBinCut1[4];
	PhotonDet1EtaBin5 = PhDetEtaBinCut1[5];
	PhotonDet1EtaBin6 = PhDetEtaBinCut1[6];
	PhotonDet1EtaBin7 = PhDetEtaBinCut1[7];
	PhotonDet1EtaBin8 = PhDetEtaBinCut1[8];
	
	PhotonDet2EtaBin1 = PhDetEtaBinCut2[1];
	PhotonDet2EtaBin2 = PhDetEtaBinCut2[2];
	PhotonDet2EtaBin3 = PhDetEtaBinCut2[3];
	PhotonDet2EtaBin4 = PhDetEtaBinCut2[4];
	PhotonDet2EtaBin5 = PhDetEtaBinCut2[5];
	PhotonDet2EtaBin6 = PhDetEtaBinCut2[6];
	PhotonDet2EtaBin7 = PhDetEtaBinCut2[7];
	PhotonDet2EtaBin8 = PhDetEtaBinCut2[8];
	
	ntTrue->Fill(PhotonTrue,PhotonTrueEtaBin1,PhotonTrueEtaBin2,PhotonTrueEtaBin3,PhotonTrueEtaBin4,PhotonTrueEtaBin5,PhotonTrueEtaBin6,PhotonTrueEtaBin7,PhotonTrueEtaBin8);
	ntMeas1->Fill(PhotonCls1,PhotonCls1EtaBin1,PhotonCls1EtaBin2,PhotonCls1EtaBin3,PhotonCls1EtaBin4,PhotonCls1EtaBin5,PhotonCls1EtaBin6,PhotonCls1EtaBin7,PhotonCls1EtaBin8);
	ntMeas2->Fill(PhotonCls2,PhotonCls2EtaBin1,PhotonCls2EtaBin2,PhotonCls2EtaBin3,PhotonCls2EtaBin4,PhotonCls2EtaBin5,PhotonCls2EtaBin6,PhotonCls2EtaBin7,PhotonCls2EtaBin8);
	ntCorr->Fill(fV0mult,fV0Amult,fV0Cmult,GlobalTracks);
	ntDet1->Fill(PhotonDet1,PhotonDet1EtaBin1,PhotonDet1EtaBin2,PhotonDet1EtaBin3,PhotonDet1EtaBin4,PhotonDet1EtaBin5,PhotonDet1EtaBin6,PhotonDet1EtaBin7,PhotonDet1EtaBin8);
	ntDet2->Fill(PhotonDet2,PhotonDet2EtaBin1,PhotonDet2EtaBin2,PhotonDet2EtaBin3,PhotonDet2EtaBin4,PhotonDet2EtaBin5,PhotonDet2EtaBin6,PhotonDet2EtaBin7,PhotonDet2EtaBin8);
	ntCent->Fill(nCentrality);
      }//VtxZ<10
      //}//GetZRes
    }//Contributors
  }//vertex
  
}

//________________________________________________________________________
void AliPMDpPbAnalysisTaskSim::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}


