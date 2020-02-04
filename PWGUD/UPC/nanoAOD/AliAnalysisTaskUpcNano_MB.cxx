/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TROOT.h"
#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TColor.h"
#include "TRandom.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliVVZERO.h"
#include "AliAODZDC.h"
#include "AliESDZDC.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliAODVertex.h"
#include "AliTOFTriggerMask.h"
#include "TObjArray.h"
#include "AliDataFile.h"
#include "TString.h"
#include "AliTimeRangeCut.h"
#include "AliOADBContainer.h"
#include "AliVVertex.h"

#include "AliESDEvent.h" 
#include "AliESDtrack.h" 
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"


// my headers
#include "AliAnalysisTaskUpcNano_MB.h"

ClassImp(AliAnalysisTaskUpcNano_MB);

using std::cout;
using std::endl;

//analysis skeleton of UPC nano AODs,

//_____________________________________________________________________________
AliAnalysisTaskUpcNano_MB::AliAnalysisTaskUpcNano_MB() 
  : AliAnalysisTaskSE(),
    fPIDResponse(0),
    fTrackCutsBit0(0),
    fTrackCutsBit1(0),
    fTrackCutsBit4(0),
    fTrackCutsBit5(0),
    fTrackCutsMatching(0),
    isMC(kFALSE),
    isESD(kFALSE),
    cutEta(0.9),
    storeRho(kFALSE),
    fOutputList(0),
    fHistEvents(0),
    fHistMCTriggers(0),
    fTreeJPsi(0),
    fTreePsi2s(0),
    fTreeGen(0),
    hTPCPIDMuonCorr(0),
    hTPCPIDMuon(0),
    hTPCPIDElectronCorr(0),
    hTPCPIDElectron(0),
    hTPCPIDPion(0),
    hTPCPIDPionCorr(0),
    hTOFPIDProton(0),
    hTOFPIDProtonCorr(0),
    hITSPIDKaon(0),
    hITSPIDKaonCorr(0),
    hTPCdEdxCorr(0),
    hTriggerCounter(0),
    hADdecision(0),
    hV0decision(0),
    hVertexZ(0),
    hVertexContrib(0),
    hTPCdEdxCorrMuon(0),
    hTPCdEdxCorrElectron(0),
    hChannelPIDCorr(0),
    hChannelPIDCorrJpsi(0),
    fPt(0),
    fY(0),
    fM(0),
    fPhi(0),
    fPtGen(0),
    fYGen(0),
    fMGen(0),
    fPhiGen(0),
    fDiLeptonM(0),
    fDiLeptonPt(0),
    fZNAenergy(0),
    fZNCenergy(0),
    fChannel(0),
    fSign(0),
    fRunNumber(0),
    fADAdecision(0),
    fADCdecision(0),
    fV0Adecision(0),
    fV0Cdecision(0),
    fNGoodTracksITS(0),
    fNGoodTracksLoose(0),
    fNGoodTracksDCA(0),
    fInEtaGen(0),
    fInEtaRec(0),
    fMatchITS(0), 
    fMatchTOF(0),
    fPIDsigma(0),
    fSPDfile(0),
    fTOFfile(0),
    fLoadedRun(-1),
    hTOFeff(0),
    hSPDeff(0),
    fTOFmask(0)
{

//Dummy constructor

}//AliAnalysisTaskUpcNano_MB


//_____________________________________________________________________________
AliAnalysisTaskUpcNano_MB::AliAnalysisTaskUpcNano_MB(const char *name) 
  : AliAnalysisTaskSE(name),
    fPIDResponse(0),
    fTrackCutsBit0(0),
    fTrackCutsBit1(0),
    fTrackCutsBit4(0),
    fTrackCutsBit5(0),
    fTrackCutsMatching(0),
    isMC(kFALSE),
    isESD(kFALSE),
    cutEta(0.9),
    storeRho(kFALSE),
    fOutputList(0),
    fHistEvents(0),
    fHistMCTriggers(0),
    fTreeJPsi(0),
    fTreePsi2s(0),
    fTreeGen(0),
    hTPCPIDMuonCorr(0),
    hTPCPIDMuon(0),
    hTPCPIDElectronCorr(0),
    hTPCPIDElectron(0),
    hTPCPIDPion(0),
    hTPCPIDPionCorr(0),
    hTOFPIDProton(0),
    hTOFPIDProtonCorr(0),
    hITSPIDKaon(0),
    hITSPIDKaonCorr(0),
    hTPCdEdxCorr(0),		
    hTriggerCounter(0),
    hADdecision(0),
    hV0decision(0),
    hVertexZ(0),
    hVertexContrib(0),
    hTPCdEdxCorrMuon(0),
    hTPCdEdxCorrElectron(0),
    hChannelPIDCorr(0),
    hChannelPIDCorrJpsi(0),
    fPt(0),
    fY(0),
    fM(0),
    fPhi(0),
    fPtGen(0),
    fYGen(0),
    fMGen(0),
    fPhiGen(0),
    fDiLeptonM(0),
    fDiLeptonPt(0),
    fZNAenergy(0),
    fZNCenergy(0),
    fChannel(0),
    fSign(0),
    fRunNumber(0),
    fADAdecision(0),
    fADCdecision(0),
    fV0Adecision(0),
    fV0Cdecision(0),
    fNGoodTracksITS(0),
    fNGoodTracksLoose(0),
    fNGoodTracksDCA(0),
    fInEtaGen(0),
    fInEtaRec(0),
    fMatchITS(0), 
    fMatchTOF(0), 
    fPIDsigma(0),
    fSPDfile(0),
    fTOFfile(0),
    fLoadedRun(-1),
    hTOFeff(0),
    hSPDeff(0),
    fTOFmask(0)
{
  for(Int_t i = 0; i<NTRIGGERINPUTS; i++) fTriggerInputsMC[i] = kFALSE;
  for(Int_t i = 0; i<NTRIGGERS; i++) fTriggers[i] = kFALSE;
  for(Int_t i = 0; i<2;  i++){ 
    fPtGenDaughter[i] = -1;
    fPtDaughter[i] = -1;
    fSignDaughter[i] = 0;
    }
  for(Int_t i = 0; i<3;  i++)fTriggerClass[i] = kFALSE;
  DefineOutput(1, TList::Class());

}//AliAnalysisTaskUpcNano_MB

//_____________________________________________________________________________
AliAnalysisTaskUpcNano_MB::~AliAnalysisTaskUpcNano_MB() 
{
  // Destructor
  
  // Destructor
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
     delete fOutputList;
     fOutputList = 0x0;
  }

}//~AliAnalysisTaskUpcNano_MB


//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::UserCreateOutputObjects()
{
  
AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  
  fTrackCutsBit0 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fTrackCutsBit0->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  fTrackCutsBit5 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
  fTrackCutsBit4 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  
  fTrackCutsBit1 = AliESDtrackCuts::GetStandardITSSATrackCuts2010(kFALSE,kTRUE);
  
  fTrackCutsMatching = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE,0);
  fTrackCutsMatching->SetEtaRange(-0.8, 0.8);
  fTrackCutsMatching->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  
  fOutputList = new TList();
  fOutputList ->SetOwner();

  TString gCutName[9] = {"Analyzed","2 TPC"," 2TPC+2ITS","4TPC","3TPC+1ITS","4TPC+ITS","2TPC+1ITS","2ITS","4ITS"};		  
  fHistEvents = new TH1D("fHistEvents","NEvents",9,0.5,9.5);
  for (Int_t i = 0; i<9; i++) fHistEvents->GetXaxis()->SetBinLabel(i+1,gCutName[i].Data());
  fOutputList->Add(fHistEvents);
  
  TString gTriggerName[NTRIGGERINPUTS+1] = {"0VBA","0VBC","0UBA","0UBC","0OMU","0OM2","0SMB","0SM2","0STP","0SH1","0STG","Offline OMU"};
  fHistMCTriggers = new TH1D("fHistMCTriggers","Fired MC trigger inputs",NTRIGGERINPUTS+1,0.5,NTRIGGERINPUTS+1.5);
  for (Int_t i = 0; i<NTRIGGERINPUTS+1; i++) fHistMCTriggers->GetXaxis()->SetBinLabel(i+1,gTriggerName[i].Data());
  fOutputList->Add(fHistMCTriggers);
  
  fTreeJPsi = new TTree("fTreeJPsi", "fTreeJPsi");
  fTreeJPsi ->Branch("fPt", &fPt, "fPt/F");
  fTreeJPsi ->Branch("fPtDaughter", &fPtDaughter[0], "fPtDaughter[2]/F");
  //fTreeJPsi ->Branch("fVectDaughter0", &fVectDaughter[0]);
  //fTreeJPsi ->Branch("fVectDaughter1", &fVectDaughter[1]);
  //fTreeJPsi ->Branch("fSignDaughter", &fSignDaughter[0], "fSignDaughter[2]/I");
  fTreeJPsi ->Branch("fY", &fY, "fY/F");
  fTreeJPsi ->Branch("fM", &fM, "fM/F");
  fTreeJPsi ->Branch("fPhi", &fPhi, "fPhi/F");
  fTreeJPsi ->Branch("fChannel", &fChannel, "fChannel/I");
  fTreeJPsi ->Branch("fSign", &fSign, "fSign/I");
  fTreeJPsi ->Branch("fZNAenergy", &fZNAenergy,"fZNAenergy/F");
  fTreeJPsi ->Branch("fZNCenergy", &fZNCenergy,"fZNCenergy/F");
  fTreeJPsi ->Branch("fZNAtime", &fZNAtime[0],"fZNAtime[4]/F");
  fTreeJPsi ->Branch("fZNCtime", &fZNCtime[0],"fZNCtime[4]/F");
  fTreeJPsi ->Branch("fPIDsigma", &fPIDsigma,"fPIDsigma/F");
  fTreeJPsi ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
  fTreeJPsi ->Branch("fInEtaRec", &fInEtaRec, "fInEtaRec/O");
  fTreeJPsi ->Branch("fTriggers", &fTriggers, Form("fTriggers[%i]/O",NTRIGGERS));
  fTreeJPsi ->Branch("fADAdecision", &fADAdecision, "fADAdecision/I");
  fTreeJPsi ->Branch("fADCdecision", &fADCdecision, "fADCdecision/I");
  fTreeJPsi ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fTreeJPsi ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");
  fTreeJPsi ->Branch("fNGoodTracksITS", &fNGoodTracksITS, "fNGoodTracksITS/I");
  fTreeJPsi ->Branch("fNGoodTracksDCA", &fNGoodTracksDCA, "fNGoodTracksDCA/I");
  fTreeJPsi ->Branch("fNGoodTracksLoose", &fNGoodTracksLoose, "fNGoodTracksLoose/I");
  //fTreeJPsi ->Branch("fTrackLenght", &fTrackLenght[0],"fTrackLenght[6]/F");
  //fTreeJPsi ->Branch("fTrackPhiPos", &fTrackPhiPos[0],"fTrackPhiPos[2]/F");
  if(isMC){
	fTreeJPsi ->Branch("fPtGen", &fPtGen, "fPtGen/F");
	fTreeJPsi ->Branch("fPtGenDaughter", &fPtGenDaughter[0], "fPtGenDaughter[2]/F");
	fTreeJPsi ->Branch("fYGen", &fYGen, "fYGen/F");
	fTreeJPsi ->Branch("fMGen", &fMGen, "fMGen/F");
	fTreeJPsi ->Branch("fPhiGen", &fPhiGen, "fPhiGen/F");
	fTreeJPsi ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], Form("fTriggerInputsMC[%i]/O",NTRIGGERINPUTS));
	}
  fOutputList->Add(fTreeJPsi);

  fTreePsi2s = new TTree("fTreePsi2s", "fTreePsi2s");
  fTreePsi2s ->Branch("fPt", &fPt, "fPt/F");
  fTreePsi2s ->Branch("fY", &fY, "fY/F");
  fTreePsi2s ->Branch("fM", &fM, "fM/F");
  fTreePsi2s ->Branch("fPhi", &fPhi, "fPhi/F");
  fTreePsi2s ->Branch("fChannel", &fChannel, "fChannel/I");
  fTreePsi2s ->Branch("fSign", &fSign, "fSign/I");
  fTreePsi2s ->Branch("fDiLeptonM", &fDiLeptonM, "fDiLeptonM/F");
  fTreePsi2s ->Branch("fDiLeptonPt", &fDiLeptonPt, "fDiLeptonPt/F");
  fTreePsi2s ->Branch("fZNAenergy", &fZNAenergy,"fZNAenergy/F");
  fTreePsi2s ->Branch("fZNCenergy", &fZNCenergy,"fZNCenergy/F");
  fTreePsi2s ->Branch("fZNAtime", &fZNAtime[0],"fZNAtime[4]/F");
  fTreePsi2s ->Branch("fZNCtime", &fZNCtime[0],"fZNCtime[4]/F");
  fTreePsi2s ->Branch("fPIDsigma", &fPIDsigma,"fPIDsigma/F");
  fTreePsi2s ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
  fTreePsi2s ->Branch("fInEtaRec", &fInEtaRec, "fInEtaRec/O");
  fTreePsi2s ->Branch("fTriggers", &fTriggers, Form("fTriggers[%i]/O",NTRIGGERS));
  fTreePsi2s ->Branch("fADAdecision", &fADAdecision, "fADAdecision/I");
  fTreePsi2s ->Branch("fADCdecision", &fADCdecision, "fADCdecision/I");
  fTreePsi2s ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  fTreePsi2s ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");
  fTreePsi2s ->Branch("fNGoodTracksITS", &fNGoodTracksITS, "fNGoodTracksITS/I");
  fTreePsi2s ->Branch("fNGoodTracksDCA", &fNGoodTracksDCA, "fNGoodTracksDCA/I");
  fTreePsi2s ->Branch("fNGoodTracksLoose", &fNGoodTracksLoose, "fNGoodTracksLoose/I");
  //fTreePsi2s ->Branch("fTrackLenght", &fTrackLenght[0],"fTrackLenght[4]/F");
  if(isMC){ 
    fTreePsi2s ->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], Form("fTriggerInputsMC[%i]/O",NTRIGGERINPUTS));
    }
  fOutputList->Add(fTreePsi2s);
    
  fTreeGen = new TTree("fTreeGen", "fTreeGen");
  fTreeGen ->Branch("fPt", &fPt, "fPt/F");
  fTreeGen ->Branch("fY", &fY, "fY/F");
  fTreeGen ->Branch("fM", &fM, "fM/F");
  fTreeGen ->Branch("fPhi", &fPhi, "fPhi/F");
  fTreeGen ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
  fTreeGen ->Branch("fInEtaGen", &fInEtaGen, "fInEtaGen/O");
  if(isMC) fOutputList->Add(fTreeGen);
  
  fTreeMatch = new TTree("fTreeMatch", "fTreeMatch");
  fTreeMatch ->Branch("fPtDaughter", &fPtDaughter, "fPtDaughter/F");
  fTreeMatch ->Branch("fNGoodTracksLoose", &fNGoodTracksLoose, "fNGoodTracksLoose/I");
  fTreeMatch ->Branch("fMatchITS", &fMatchITS, "fMatchITS/O");
  fTreeMatch ->Branch("fMatchTOF", &fMatchTOF, "fMatchTOF/O");
  fTreeMatch ->Branch("fTriggerClass", &fTriggerClass[0],"fTriggerClass[3]/O");
  fOutputList->Add(fTreeMatch);
   
  hTPCPIDMuonCorr = new TH2D("hTPCPIDMuonCorr"," ",100,-10.0,10.0,100,-10.0,10.0);
  hTPCPIDMuonCorr->GetXaxis()->SetTitle("TPC PID N #sigma (#mu)");
  hTPCPIDMuonCorr->GetYaxis()->SetTitle("TPC PID N #sigma (#mu)");
  fOutputList->Add(hTPCPIDMuonCorr);
  
  hTPCPIDElectronCorr = new TH2D("hTPCPIDElectronCorr"," ",100,-10.0,10.0,100,-10.0,10.0);
  hTPCPIDElectronCorr->GetXaxis()->SetTitle("TPC PID N #sigma (e)");
  hTPCPIDElectronCorr->GetYaxis()->SetTitle("TPC PID N #sigma (e)");
  fOutputList->Add(hTPCPIDElectronCorr);

  hTPCPIDPionCorr = new TH2D("hTPCPIDPionCorr"," ",100,-10.0,10.0,100,-10.0,10.0);
  hTPCPIDPionCorr->GetXaxis()->SetTitle("TPC PID N #sigma (#pi)");
  hTPCPIDPionCorr->GetYaxis()->SetTitle("TPC PID N #sigma (#pi)");
  fOutputList->Add(hTPCPIDPionCorr);
  
  hTOFPIDProtonCorr = new TH2D("hTOFPIDProtonCorr"," ",100,-10.0,10.0,100,-10.0,10.0);
  hTOFPIDProtonCorr->GetXaxis()->SetTitle("TOF PID N #sigma (p)");
  hTOFPIDProtonCorr->GetYaxis()->SetTitle("TOF PID N #sigma (p)");
  fOutputList->Add(hTOFPIDProtonCorr);
  
  hITSPIDKaonCorr = new TH2D("hITSPIDKaonCorr"," ",100,-10.0,10.0,100,-10.0,10.0);
  hITSPIDKaonCorr->GetXaxis()->SetTitle("ITS PID N #sigma (K)");
  hITSPIDKaonCorr->GetYaxis()->SetTitle("ITS PID N #sigma (K)");
  fOutputList->Add(hITSPIDKaonCorr);
  
  hTPCPIDMuon = new TH1D("hTPCPIDMuon"," ",200,-10,10);
  hTPCPIDMuon->GetXaxis()->SetTitle("TPC PID N #sigma (#mu)");
  fOutputList->Add(hTPCPIDMuon);
  
  hTPCPIDElectron = new TH1D("hTPCPIDElectron"," ",200,-10,10);
  hTPCPIDElectron->GetXaxis()->SetTitle("TPC PID N #sigma (#e)");
  fOutputList->Add(hTPCPIDElectron);
  
  hTPCPIDPion = new TH1D("hTPCPIDPion"," ",200,-10,10);
  hTPCPIDPion->GetXaxis()->SetTitle("TPC PID N #sigma (#pi)");
  fOutputList->Add(hTPCPIDPion);
  
  hTOFPIDProton = new TH1D("hTOFPIDProton"," ",200,-10,10);
  hTOFPIDProton->GetXaxis()->SetTitle("TOF PID N #sigma (p)");
  fOutputList->Add(hTOFPIDProton);
  
  hITSPIDKaon = new TH1D("hITSPIDKaon"," ",200,-10,10);
  hITSPIDKaon->GetXaxis()->SetTitle("ITS PID N #sigma (K)");
  fOutputList->Add(hITSPIDKaon);
  
  hTPCdEdxCorr = new TH2D("hTPCdEdxCorr"," ",200,0,200,200,0,200);
  hTPCdEdxCorr->GetXaxis()->SetTitle("dE/dx^{TPC} (a.u.)");
  hTPCdEdxCorr->GetYaxis()->SetTitle("dE/dx^{TPC} (a.u.)");
  fOutputList->Add(hTPCdEdxCorr);
  
  hTriggerCounter = new TH2I("hTriggerCounter","Number of analyzed UPC triggers per run",3,1,4,3000,295000,298000);
  fOutputList->Add(hTriggerCounter);
  
  hADdecision = new TH2I("hADdecision","hADdecision",7,-2,5,7,-2,5);
  fOutputList->Add(hADdecision);
  hV0decision = new TH2I("hV0decision","hV0decision",7,-2,5,7,-2,5);
  fOutputList->Add(hV0decision);
  
  hVertexZ = new TH1D("hVertexZ"," ",600,-30,30);
  fOutputList->Add(hVertexZ);
  hVertexContrib = new TH1D("hVertexContrib"," ",103,-2,100);
  fOutputList->Add(hVertexContrib);
    
  hChannelPIDCorr = new TH2D("hChannelPIDCorr"," ",100,0.0,20.0,100,0.0,20.0);
  hChannelPIDCorr->GetXaxis()->SetTitle("TPC PID N #sigma (e)");
  hChannelPIDCorr->GetYaxis()->SetTitle("TPC PID N #sigma (#mu)");
  fOutputList->Add(hChannelPIDCorr);
  
  hChannelPIDCorrJpsi = new TH2D("hChannelPIDCorrJpsi"," ",100,0.0,20.0,100,0.0,20.0);
  hChannelPIDCorrJpsi->GetXaxis()->SetTitle("TPC PID N #sigma (e)");
  hChannelPIDCorrJpsi->GetYaxis()->SetTitle("TPC PID N #sigma (#mu)");
  fOutputList->Add(hChannelPIDCorrJpsi);	
  
  hTPCdEdxCorrMuon = new TH2D("hTPCdEdxCorrMuon"," ",200,0,200,200,0,200);
  hTPCdEdxCorrMuon->GetXaxis()->SetTitle("dE/dx^{TPC} (a.u.)");
  hTPCdEdxCorrMuon->GetYaxis()->SetTitle("dE/dx^{TPC} (a.u.)");
  fOutputList->Add(hTPCdEdxCorrMuon);
  
  hTPCdEdxCorrElectron = new TH2D("hTPCdEdxCorrElectron"," ",200,0,200,200,0,200);
  hTPCdEdxCorrElectron->GetXaxis()->SetTitle("dE/dx^{TPC} (a.u.)");
  hTPCdEdxCorrElectron->GetYaxis()->SetTitle("dE/dx^{TPC} (a.u.)");
  fOutputList->Add(hTPCdEdxCorrElectron);	
     
  PostData(1, fOutputList);

}//UserCreateOutputObjects


//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::UserExec(Option_t *) 
{

  AliVEvent *fEvent = InputEvent();
  if(!fEvent) return;

  fRunNumber = fEvent->GetRunNumber();
  
  fTimeRangeCut.InitFromEvent(InputEvent());
  if(fTimeRangeCut.CutEvent(InputEvent()))return;
  
  TString trigger = fEvent->GetFiredTriggerClasses();
  
  if(!isMC){
    if(fRunNumber>=295881 && fRunNumber<296594)
      if(!trigger.Contains("CCUP29-B-SPD2-CENTNOTRD") && !trigger.Contains("CCUP30-B-SPD2-CENTNOTRD") && !trigger.Contains("CCUP31-B-SPD2-CENTNOTRD"))return;
    if(fRunNumber>=296594)
      if(!trigger.Contains("CCUP29-U-SPD2-CENTNOTRD") && !trigger.Contains("CCUP30-B-SPD2-CENTNOTRD") && !trigger.Contains("CCUP31-B-SPD2-CENTNOTRD"))return;
    if(fRunNumber<295881)
      if(!trigger.Contains("CCUP29-B-NOPF-CENTNOTRD") && !trigger.Contains("CCUP30-B-NOPF-CENTNOTRD") && !trigger.Contains("CCUP31-B-NOPF-CENTNOTRD"))return;
    }

  UInt_t fL0inputs = fEvent->GetHeader()->GetL0TriggerInputs();
  fTriggers[0] = trigger.Contains("CCUP29-B-SPD2-CENTNOTRD");
  fTriggers[1] = trigger.Contains("CCUP29-B-NOPF-CENTNOTRD");
  fTriggers[2] = trigger.Contains("CCUP29-U-SPD2-CENTNOTRD");
  fTriggers[3] = trigger.Contains("CCUP30-B-NOPF-CENTNOTRD");
  fTriggers[4] = trigger.Contains("CCUP30-B-SPD2-CENTNOTRD");
  fTriggers[5] = trigger.Contains("CCUP31-B-NOPF-CENTNOTRD");
  fTriggers[6] = trigger.Contains("CCUP31-B-SPD2-CENTNOTRD");
  fTriggers[7] =  fL0inputs & (1 << 11); //OM2
  fTriggers[8] =  fL0inputs & (1 << 12); //OMU
  
  if(trigger.Contains("CCUP29-B") || trigger.Contains("CCUP29-U"))hTriggerCounter->Fill(1,fRunNumber); 
  if(trigger.Contains("CCUP30-B"))hTriggerCounter->Fill(2,fRunNumber); 
  if(trigger.Contains("CCUP31-B"))hTriggerCounter->Fill(3,fRunNumber);
  
  AliVVZERO *fV0data = fEvent->GetVZEROData();
  AliVAD *fADdata = fEvent->GetADData();
    
  if(isMC){
    if(fRunNumber != fLoadedRun){
      if(!fSPDfile)fSPDfile = AliDataFile::OpenOADB("PWGUD/UPC/SPDEfficiency18qr.root");
      if(!fTOFfile)fTOFfile = AliDataFile::OpenOADB("PWGUD/TOFTriggerEfficiency.root");
      AliOADBContainer* fTOFcont = (AliOADBContainer*)fTOFfile->Get("TOFTriggerEfficiency");
      hTOFeff  = (TH2F*)fTOFcont->GetObject(fRunNumber,"Default");
      
      if(fSPDfile->Get(Form("eff%i",fRunNumber))) hSPDeff  = (TH1D*) fSPDfile->Get(Form("eff%i",fRunNumber));
      Int_t tempRun = fRunNumber;
      while(!hSPDeff){
        tempRun--;
        hSPDeff  = (TH1D*) fSPDfile->Get(Form("eff%i",tempRun));
      }
      fLoadedRun = fRunNumber;
      }
    RunMC(fEvent);
    }
  
  for(Int_t i = 0; i<NTRIGGERINPUTS; i++)if(fTriggerInputsMC[i])fHistMCTriggers->Fill(i+1);
    
  if(isESD){
  	AliESDZDC *fZDCdata = (AliESDZDC*)fEvent->GetZDCData();
  	fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  	fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];

	Int_t detChZNA  = fZDCdata->GetZNATDCChannel();
        Int_t detChZNC  = fZDCdata->GetZNCTDCChannel();
	if (fEvent->GetRunNumber()>=245726 && fEvent->GetRunNumber()<=245793) detChZNA = 10;
  	for (Int_t i=0;i<4;i++){ 
  		fZNAtime[i] = fZDCdata->GetZDCTDCCorrected(detChZNA,i);
  		fZNCtime[i] = fZDCdata->GetZDCTDCCorrected(detChZNC,i);
		}
  }
  else{
  	AliAODZDC *fZDCdata = (AliAODZDC*)fEvent->GetZDCData();
  	fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  	fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];

  	for (Int_t i=0;i<4;i++){ 
  		fZNAtime[i] = fZDCdata->GetZNATDCm(i);
  		fZNCtime[i] = fZDCdata->GetZNCTDCm(i);
		}
  }
  
  const AliVVertex *fVertex = fEvent->GetPrimaryVertex();
  hVertexZ->Fill(fVertex->GetZ());
  hVertexContrib->Fill(fVertex->GetNContributors());
  
  if(fVertex->GetNContributors()<2)return;
  if(TMath::Abs(fVertex->GetZ())>15)return;
 
 
  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  hV0decision->Fill(fV0Adecision,fV0Cdecision);
  
  fADAdecision = fADdata->GetADADecision();
  fADCdecision = fADdata->GetADCDecision();
  hADdecision->Fill(fADAdecision,fADCdecision);
  
  fHistEvents->Fill(1);
  //cout<<"Event, tracks = "<<fEvent ->GetNumberOfTracks()<<endl; 
  
  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  TParticlePDG *partKaon = pdgdat->GetParticle( 321 );
  Float_t kaonMass = partKaon->Mass();
  
  TParticlePDG *partMuon = pdgdat->GetParticle( 13 );
  Float_t muonMass = partMuon->Mass();
  
  TParticlePDG *partElectron = pdgdat->GetParticle( 11 );
  Float_t electronMass = partElectron->Mass();
  
  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  Float_t pionMass = partPion->Mass();
  
  TParticlePDG *partProton = pdgdat->GetParticle( 2212 );
  Float_t protonMass = partProton->Mass();
  
  Short_t qTrack[5];
  TLorentzVector vMuon[5],vElectron[5],vProton[5], vJPsiCandidate, vJPsiGenerated;
  TLorentzVector vKaon[5], vPhiCandidate;
  TLorentzVector vPion[5], vRhoCandidate;

  Float_t nSigmaMuon[5], nSigmaElectron[5], nSigmaPion[5], nSigmaProton[5], dEdx[5];
  Short_t qPion[5];
  TLorentzVector vLepton[5], vDilepton, vPsi2sCandidate;
  Short_t qLepton[5];
  UInt_t nPion = 0, nLepton = 0, nHighPt = 0;
  UInt_t nGoodTracksTPC=0;
  fNGoodTracksDCA = 0;
  fNGoodTracksITS = 0; 
  UInt_t nGoodTracksSPD=0;
  UInt_t nGoodTracksTOF=0;
  fNGoodTracksLoose=0;
  Int_t TrackIndexTPC[5] = {-1,-1,-1,-1,-1};
  Int_t TrackIndexITS[5] = {-1,-1,-1,-1,-1};
  Int_t TrackIndexALL[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  Double_t TrackPtTPC[5]={0,0,0,0,0};
  Double_t TrackPtALL[10]={0,0,0,0,0,0,0,0,0,0};
  Double_t MeanPt = -1;
  //Track loop
  for(Int_t iTrack=0; iTrack < fEvent->GetNumberOfTracks(); iTrack++) {
  Bool_t goodTPCTrack = kTRUE;
  Bool_t goodITSTrack = kTRUE;
    if(isESD){ 
    	AliESDtrack *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(iTrack));
	if( !trk ) continue;
	if(fTrackCutsBit0->AcceptTrack(trk) && (trk->HasPointOnITSLayer(0) || trk->HasPointOnITSLayer(1))) fNGoodTracksLoose++;
	
    	if(!fTrackCutsBit4->AcceptTrack(trk)) goodTPCTrack = kFALSE;
	else{ 
	  if(trk->HasPointOnITSLayer(0) && trk->HasPointOnITSLayer(1))nGoodTracksSPD++;
	  if(fTrackCutsBit5->AcceptTrack(trk))fNGoodTracksDCA++;
	  }

    	if(!fTrackCutsBit1->AcceptTrack(trk)) goodITSTrack = kFALSE;
	
	}
    else{ 
    	AliAODTrack *trk = dynamic_cast<AliAODTrack*>(fEvent->GetTrack(iTrack));
    	if( !trk ) continue;
    
    	if(trk->TestFilterBit(1<<0) && (trk->HasPointOnITSLayer(0) || trk->HasPointOnITSLayer(1)))fNGoodTracksLoose++;
    	if(!(trk->TestFilterBit(1<<5)))goodTPCTrack = kFALSE;
    	else{
    		if(trk->HasPointOnITSLayer(0) && trk->HasPointOnITSLayer(1))nGoodTracksSPD++;
		if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trk) == AliPIDResponse::kDetPidOk)nGoodTracksTOF++;
    		}
    
    	if(!(trk->TestFilterBit(1<<1))) goodITSTrack = kFALSE;
	}
    AliVTrack *trk = dynamic_cast<AliVTrack*>(fEvent->GetTrack(iTrack));

    if(goodTPCTrack){
    	TrackIndexTPC[nGoodTracksTPC] = iTrack;
    	TrackPtTPC[nGoodTracksTPC] = trk->Pt();
	TrackIndexALL[nGoodTracksTPC+fNGoodTracksITS] = iTrack;
    	TrackPtALL[nGoodTracksTPC+fNGoodTracksITS] = trk->Pt();
    	nGoodTracksTPC++;
	}
     	
    if(goodITSTrack){
    	TrackIndexITS[fNGoodTracksITS] = iTrack;
	TrackIndexALL[nGoodTracksTPC+fNGoodTracksITS] = iTrack;
    	TrackPtALL[nGoodTracksTPC+fNGoodTracksITS] = trk->Pt();
    	fNGoodTracksITS++;
    	}
     
     
    if(nGoodTracksTPC > 4 || fNGoodTracksITS > 4) break;
    }
    
  //{"Analyzed","2 TPC"," 2TPC+2ITS","4TPC","3TPC+1ITS","4TPC+ITS","2TPC+1ITS","2ITS","4ITS"};
  if(nGoodTracksTPC == 2 && fNGoodTracksITS == 0)fHistEvents->Fill(2);
  if(nGoodTracksTPC == 2 && fNGoodTracksITS == 2)fHistEvents->Fill(3);
  if(nGoodTracksTPC == 4 && fNGoodTracksITS == 0)fHistEvents->Fill(4);
  if(nGoodTracksTPC == 3 && fNGoodTracksITS == 1)fHistEvents->Fill(5);
  if(nGoodTracksTPC == 4 && fNGoodTracksITS != 0)fHistEvents->Fill(6);
  if(nGoodTracksTPC == 2 && fNGoodTracksITS != 0 && fNGoodTracksITS != 2)fHistEvents->Fill(7);
  if(nGoodTracksTPC == 0 && fNGoodTracksITS == 2)fHistEvents->Fill(8);
  if(nGoodTracksTPC == 0 && fNGoodTracksITS == 4)fHistEvents->Fill(9);
  
  Int_t crossedFO[4];
  TBits fFOCrossedChips(1200); 
  const AliVMultiplicity *mult = fEvent->GetMultiplicity();
  TBits fFOFiredChips = mult->GetFastOrFiredChips();
  
  fInEtaRec = kTRUE;   
  if(nGoodTracksTPC+fNGoodTracksITS == 4 && nGoodTracksTPC > 1 && nGoodTracksSPD > 1){
  	fFOCrossedChips.ResetAllBits(kFALSE);
    	MeanPt = GetMedian(TrackPtALL);
  	for(Int_t iTrack=0; iTrack<4; iTrack++) {

	if(isESD){ 
    	  AliESDtrack *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(TrackIndexALL[iTrack]));
	  if( !trk ) continue;
	  crossedFO[0] = trk->GetITSModuleIndex(0);
	  crossedFO[1] = trk->GetITSModuleIndex(1);
	  crossedFO[2] = trk->GetITSModuleIndex(6);
	  crossedFO[3] = trk->GetITSModuleIndex(7);
	  SetCrossed(crossedFO, fFOCrossedChips);
	  if(trk->GetInnerParam())fTrackLenght[iTrack] = trk->GetLengthInActiveZone(1, 2.0, 220, fEvent->GetMagneticField());
	  else fTrackLenght[iTrack] = trk->GetLengthInActiveZone(0, 2.0, 220, fEvent->GetMagneticField());
	  }
	
	AliVTrack *trk = dynamic_cast<AliVTrack*>(fEvent->GetTrack(TrackIndexALL[iTrack]));
	if(TMath::Abs(trk->Eta())>cutEta) fInEtaRec = kFALSE;
	
	if(trk->Pt() > 1.0) nHighPt++;
      		if(trk->Pt() > MeanPt){   
      			qLepton[nLepton] = trk->Charge();
			Float_t fPIDTPCMuon = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
    			Float_t fPIDTPCElectron = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
			Float_t fPIDTOFProton = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);

			nSigmaMuon[nLepton] = fPIDTPCMuon;
    			nSigmaElectron[nLepton] = fPIDTPCElectron;
			
			nSigmaProton[nLepton] = fPIDTOFProton;
			vProton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), protonMass);
			
			dEdx[nLepton] = trk->GetTPCsignal();
			
			if(TMath::Abs(fPIDTPCMuon) < TMath::Abs(fPIDTPCElectron)){	
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
				nLepton++;
				}
      			if(TMath::Abs(fPIDTPCMuon) > TMath::Abs(fPIDTPCElectron)){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
				nLepton++;
				}
			}
		else{
			qPion[nPion] = trk->Charge();
			vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
			nPion++;
			}
    		}
		
	fFOCrossFiredChips = fFOCrossedChips & fFOFiredChips;
        fTriggers[9] = IsSTGFired(fFOCrossFiredChips,fRunNumber >= 295753 ? 9 : 3);
	
	if((nLepton == 2) && (nPion == 2)){
		vPsi2sCandidate = vLepton[0]+vLepton[1]+vPion[0]+vPion[1];
  		vDilepton = vLepton[0]+vLepton[1];
		fDiLeptonM = vDilepton.M();
  		fDiLeptonPt = vDilepton.Pt();
		if((qLepton[0]*qLepton[1] < 0) && (qPion[0]*qPion[1] < 0)) fSign = -1;
		else fSign = 1;
		Float_t nSigmaDistMuon = TMath::Sqrt(TMath::Power(nSigmaMuon[0],2) + TMath::Power(nSigmaMuon[1],2));
		Float_t nSigmaDistElectron = TMath::Sqrt(TMath::Power(nSigmaElectron[0],2) + TMath::Power(nSigmaElectron[1],2));
		Float_t nSigmaDistProton = TMath::Sqrt(TMath::Power(nSigmaProton[0],2) + TMath::Power(nSigmaProton[1],2));
		fChannel = 0;
		if((nSigmaDistMuon < nSigmaDistElectron)){
			fPIDsigma = nSigmaDistMuon; 
			if(fNGoodTracksITS == 0)fChannel = 1;
			else fChannel = 10;
			}
		if((nSigmaDistMuon > nSigmaDistElectron)){
			fPIDsigma = nSigmaDistElectron;
			if(fNGoodTracksITS == 0)fChannel = -1;
			else fChannel = -10;
			}
		if(nSigmaDistProton < 3){ 
			vJPsiCandidate = vProton[0]+vProton[1]+vPion[0]+vPion[1];
			fPIDsigma = nSigmaDistProton;
			fChannel = 3;
			FillTree(fTreeJPsi,vJPsiCandidate);
			}
		
		if(nHighPt > 0 )FillTree(fTreePsi2s,vPsi2sCandidate);
		}
  	}
	
  //Two track loop
  
  TLorentzVector vMC, vLabelPart;
  Int_t nTOFPID = 0;
  fInEtaRec = kTRUE;
  if(nGoodTracksTPC == 2 && nGoodTracksSPD == 2){
  	fFOCrossedChips.ResetAllBits(kFALSE);
  	for(Int_t iTrack=0; iTrack<2; iTrack++) {

	if(isESD){ 
    	  AliESDtrack *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(TrackIndexTPC[iTrack]));
	  if( !trk ) continue;
	  crossedFO[0] = trk->GetITSModuleIndex(0);
	  crossedFO[1] = trk->GetITSModuleIndex(1);
	  crossedFO[2] = trk->GetITSModuleIndex(6);
	  crossedFO[3] = trk->GetITSModuleIndex(7);
	  SetCrossed(crossedFO, fFOCrossedChips);
	  fTrackLenght[iTrack] = trk->GetLengthInActiveZone(1, 1.5, 220, fEvent->GetMagneticField());
	  fTrackLenght[iTrack+2] = trk->GetLengthInActiveZone(1, 2.0, 220, fEvent->GetMagneticField());
	  fTrackLenght[iTrack+4] = trk->GetLengthInActiveZone(1, 2.5, 220, fEvent->GetMagneticField());
	  fTrackPhiPos[iTrack] = trk->GetInnerParam()->PhiPos();
	  }
	  
    	AliVTrack *trk = dynamic_cast<AliVTrack*>(fEvent->GetTrack(TrackIndexTPC[iTrack]));
	
	if(isMC && trk->GetLabel() >= 0){
	  AliMCEvent *mc = MCEvent();
          if(!mc) return;
          AliMCParticle *mcPart = (AliMCParticle*) mc->GetTrack(trk->GetLabel());
	  
          fPtGenDaughter[iTrack] = mcPart->Pt();
          TParticlePDG *partGen = pdgdat->GetParticle(mcPart->PdgCode());
          vLabelPart.SetXYZM(mcPart->Px(),mcPart->Py(), mcPart->Pz(),partGen->Mass());
          vMC += vLabelPart;
	  }
	fPtDaughter[iTrack] = trk->Pt();
	fSignDaughter[iTrack] = trk->Charge();
	
	if(TMath::Abs(trk->Eta())>cutEta) fInEtaRec = kFALSE;
	
	Float_t fPIDTPCMuon = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
    	Float_t fPIDTPCElectron = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
	Float_t fPIDTPCPion = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
	
	Float_t fPIDTPCProton = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
	Float_t fPIDTOFProton = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);
	
    	qTrack[iTrack] = trk->Charge();
	
	if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trk) == AliPIDResponse::kDetPidOk) nTOFPID++;
	
	vElectron[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
    	vMuon[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
    	nSigmaMuon[iTrack] = fPIDTPCMuon;
    	nSigmaElectron[iTrack] = fPIDTPCElectron;
    
    	vPion[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
    	nSigmaPion[iTrack] = fPIDTPCPion;
	
    	vProton[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), protonMass);
	if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trk) == AliPIDResponse::kDetPidOk) nSigmaProton[iTrack] = fPIDTOFProton;
	else nSigmaProton[iTrack] = fPIDTPCProton;
 	
	dEdx[iTrack] = trk->GetTPCsignal();
    	}
	
  
  fFOCrossFiredChips = fFOCrossedChips & fFOFiredChips;
  fTriggers[9] = IsSTGFired(fFOCrossFiredChips,fRunNumber >= 295753 ? 9 : 3);

  Float_t nSigmaDistMuon = TMath::Sqrt(TMath::Power(nSigmaMuon[0],2) + TMath::Power(nSigmaMuon[1],2));
  Float_t nSigmaDistPion = TMath::Sqrt(TMath::Power(nSigmaPion[0],2) + TMath::Power(nSigmaPion[1],2));
  Float_t nSigmaDistElectron = TMath::Sqrt(TMath::Power(nSigmaElectron[0],2) + TMath::Power(nSigmaElectron[1],2));
  Float_t nSigmaDistProton = TMath::Sqrt(TMath::Power(nSigmaProton[0],2) + TMath::Power(nSigmaProton[1],2));
  
  hChannelPIDCorr->Fill(nSigmaDistElectron,nSigmaDistMuon);
  
  if(qTrack[0]*qTrack[1]<0)fSign = -1;
  if(qTrack[0]*qTrack[1]>0)fSign = 1;

  if(nSigmaDistProton < 4){ 
  	  fPIDsigma = nSigmaDistProton;
	  if(nTOFPID == 0) fPIDsigma = 666;
  	  vJPsiCandidate = vProton[0]+vProton[1];
	  fVectDaughter[0] = vProton[0];
	  fVectDaughter[1] = vProton[1];
  	  fChannel = 2;
  	  if(!isMC) FillTree(fTreeJPsi,vJPsiCandidate);
	  else FillTree(fTreeJPsi,vJPsiCandidate,vMC);
  	  }
  if(nSigmaDistMuon < nSigmaDistElectron){
  	  fPIDsigma = nSigmaDistMuon; 
  	  vJPsiCandidate = vMuon[0]+vMuon[1];
	  fVectDaughter[0] = vMuon[0];
	  fVectDaughter[1] = vMuon[1];
  	  fChannel = 1;
  	  if(!isMC) FillTree(fTreeJPsi,vJPsiCandidate);
	  else FillTree(fTreeJPsi,vJPsiCandidate,vMC);
  	  }
  if(nSigmaDistPion < nSigmaDistElectron){
  	  fPIDsigma = nSigmaDistPion; 
	  fChannel = 0;
  	  vRhoCandidate = vPion[0]+vPion[1];
	  fVectDaughter[0] = vPion[0];
	  fVectDaughter[1] = vPion[1];
  	  if(storeRho)FillTree(fTreeJPsi,vRhoCandidate);
  	  }  
  if(nSigmaDistMuon > nSigmaDistElectron){ 
  	  fPIDsigma = nSigmaDistElectron;
  	  vJPsiCandidate = vElectron[0]+vElectron[1];
	  fVectDaughter[0] = vElectron[0];
	  fVectDaughter[1] = vElectron[1];
  	  fChannel = -1;
	  if(!isMC) FillTree(fTreeJPsi,vJPsiCandidate);
	  else FillTree(fTreeJPsi,vJPsiCandidate,vMC);
  	  }

  hTPCPIDPionCorr->Fill(nSigmaPion[0],nSigmaPion[1]);
  hTPCPIDPion->Fill(nSigmaPion[0]);hTPCPIDPion->Fill(nSigmaPion[1]);

  if(vJPsiCandidate.M()>2.2 && vJPsiCandidate.M()<5.0 && fSign == -1){
          hChannelPIDCorrJpsi->Fill(nSigmaDistElectron,nSigmaDistMuon); 
	  if(nSigmaDistMuon < nSigmaDistElectron){
	    hTPCdEdxCorrMuon->Fill(dEdx[0],dEdx[1]);
	    hTPCdEdxCorrMuon->Fill(dEdx[1],dEdx[0]);
	    }
	  if(nSigmaDistMuon > nSigmaDistElectron){
	    hTPCdEdxCorrElectron->Fill(dEdx[0],dEdx[1]);
  	    hTPCdEdxCorrElectron->Fill(dEdx[1],dEdx[0]);
	    } 
  	  hTPCdEdxCorr->Fill(dEdx[0],dEdx[1]);
  	  hTPCdEdxCorr->Fill(dEdx[1],dEdx[0]);
  	  hTPCPIDMuonCorr->Fill(nSigmaMuon[0],nSigmaMuon[1]);
  	  hTPCPIDMuonCorr->Fill(nSigmaMuon[1],nSigmaMuon[0]);
  	  hTPCPIDMuon->Fill(nSigmaMuon[0]);hTPCPIDMuon->Fill(nSigmaMuon[1]);
  	  hTPCPIDElectronCorr->Fill(nSigmaElectron[0],nSigmaElectron[1]);
  	  hTPCPIDElectronCorr->Fill(nSigmaElectron[1],nSigmaElectron[0]);
  	  hTPCPIDElectron->Fill(nSigmaElectron[0]);hTPCPIDElectron->Fill(nSigmaElectron[1]); 
  	  hTOFPIDProtonCorr->Fill(nSigmaProton[0],nSigmaProton[1]);
  	  hTOFPIDProtonCorr->Fill(nSigmaProton[1],nSigmaProton[0]);
  	  hTOFPIDProton->Fill(nSigmaProton[0]);hTOFPIDProton->Fill(nSigmaProton[1]);
  	  }
  }
  
  //
  if(isESD){
  
  Int_t TrackIndexSyst[5] = {-1,-1,-1,-1,-1};
  fNGoodTracksLoose=0;
  
  Bool_t refit    = fTrackCutsMatching->GetRequireITSRefit();
  Float_t chi2tpc = fTrackCutsMatching->GetMaxChi2TPCConstrainedGlobal();
  Float_t chi2its = fTrackCutsMatching->GetMaxChi2PerClusterITS();
  
  Bool_t isCUP29, isCUP30, isCUP31;
  if(isMC){
    isCUP29 = !fTriggerInputsMC[0] && !fTriggerInputsMC[1] && !fTriggerInputsMC[2] && !fTriggerInputsMC[3] && fTriggerInputsMC[10];
    isCUP30 = !fTriggerInputsMC[0] && !fTriggerInputsMC[1] && !fTriggerInputsMC[2] && !fTriggerInputsMC[3] && fTriggerInputsMC[10] && fTriggerInputsMC[5];
    isCUP31 = !fTriggerInputsMC[0] && !fTriggerInputsMC[1] && !fTriggerInputsMC[2] && !fTriggerInputsMC[3] && fTriggerInputsMC[10] && fTriggerInputsMC[4];
    }
  else{
    isCUP29 = fTriggers[0] || fTriggers[1] || fTriggers[2];
    isCUP30 = fTriggers[3] || fTriggers[4];
    isCUP31 = fTriggers[5] || fTriggers[6];
    }
  fTriggerClass[0] = isCUP29; fTriggerClass[1] = isCUP30; fTriggerClass[2] = isCUP31;
  
  // remove all ITS requirements
  fTrackCutsMatching->SetRequireITSRefit(kFALSE);
  fTrackCutsMatching->SetMaxChi2TPCConstrainedGlobal(99999.);
  fTrackCutsMatching->SetMaxChi2PerClusterITS(999999.);
  fTrackCutsMatching->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
  //DCA
  fTrackCutsMatching->SetMaxDCAToVertexXY(2.4);
  fTrackCutsMatching->SetMaxDCAToVertexZ(3.2);
  
  for(Int_t iTrack=0; iTrack < fEvent->GetNumberOfTracks(); iTrack++) {
    AliESDtrack *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(iTrack));
    if( !trk ) continue;
    if(fTrackCutsMatching->AcceptTrack(trk)){
      TrackIndexSyst[fNGoodTracksLoose] = iTrack; 
      fNGoodTracksLoose++;
      }
    if(fNGoodTracksLoose>4)break; 
    }
  //set the ITS requirements
  fTrackCutsMatching->SetRequireITSRefit(refit);
  fTrackCutsMatching->SetMaxChi2TPCConstrainedGlobal(chi2tpc);
  fTrackCutsMatching->SetMaxChi2PerClusterITS(chi2its);
  //set the SPD cluster requirement
  fTrackCutsMatching->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
  
  if(fNGoodTracksLoose == 2 || fNGoodTracksLoose == 4){
    for(Int_t iTrack=0; iTrack < fNGoodTracksLoose; iTrack++) {
      AliESDtrack *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(TrackIndexSyst[iTrack]));
      if( !trk ) continue;
      fMatchITS = kFALSE;
      fMatchTOF = kFALSE;
      fPtDaughter[iTrack] = trk->Pt();      
      if(fTrackCutsMatching->AcceptTrack(trk)){
        fMatchITS = kTRUE; 
        if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trk) == AliPIDResponse::kDetPidOk)fMatchTOF = kTRUE;
        }
      fTreeMatch->Fill();
      }
    }
  }
  
  PostData(1, fOutputList);

}//UserExec

//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::SetCrossed(Int_t spd[4], TBits &crossed){

  Int_t chipId2;
  for(Int_t iLayer = 0; iLayer<4 ;iLayer++)
    if(spd[iLayer]>0) { crossed.SetBitNumber(GetChipId(spd[iLayer],chipId2)); crossed.SetBitNumber(chipId2); }
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskUpcNano_MB::GetChipId(Int_t index, Int_t &chipId2, Bool_t debug){
  Int_t status   = (index%1000000)/100000;
  Int_t iModule  = index/1000000;           // 0 - 239
  Int_t iPhi     = iModule/4;               // 0-19 - inner, 20-59 outer
  Int_t iModuleZ = iModule%4;               // 0-3
  Int_t iSign    = (index%100000)/10000;    // 1-4
  Int_t signZ    = iPhi<20 ? (iSign%2==1 ? 1 : -1) : (iSign%2==0 ? 1 : -1); // 1 or -1
  Int_t iX       = (index%10000)/100;       // ??
  Int_t iZ       = index%100;               // 0-36 [mm]
  Int_t signZiZ  = (36-signZ*iZ);
  Int_t chipId   = iModule*5+signZiZ*5/72;
  if (chipId<0) return 1200;
  if (chipId>=1200) return 1201;
  if (signZiZ<0) return 1202;
  if (signZiZ>72) return 1203;
  if (signZiZ==72 && chipId%20==0 && chipId>=400) return 1204;
  chipId2=chipId;
  
  if (signZiZ==0  && chipId%20!=0)  chipId2=chipId-1;
  if (signZiZ==72 && chipId%20!=19) chipId2=chipId+1;
  if (signZiZ==13)  chipId2=chipId+1;
  if (signZiZ==14)  chipId2=chipId+1;
  if (signZiZ==15)  chipId2=chipId-1;
  if (signZiZ==16)  chipId2=chipId-1;
  if (signZiZ==27)  chipId2=chipId+1;
  if (signZiZ==28)  chipId2=chipId+1;
  if (signZiZ==29)  chipId2=chipId-1;
  if (signZiZ==30)  chipId2=chipId-1;
  if (signZiZ==42)  chipId2=chipId+1;
  if (signZiZ==43)  chipId2=chipId+1;
  if (signZiZ==44)  chipId2=chipId-1;
  if (signZiZ==45)  chipId2=chipId-1;
  if (signZiZ==56)  chipId2=chipId+1;
  if (signZiZ==57)  chipId2=chipId+1;
  if (signZiZ==58)  chipId2=chipId-1;
  if (signZiZ==59)  chipId2=chipId-1;
  if (debug) printf("%4i %4i %3i %3i %3i\n",chipId,chipId2,iX,signZiZ,iSign);
  return chipId;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUpcNano_MB::IsSTGFired(TBits bits, Int_t dphiMin, Int_t dphiMax, Bool_t tolerance){
  Int_t n1 = bits.CountBits(400);
  Int_t n0 = bits.CountBits()-n1;
  //cout<<n0<<" "<<n1<<endl;
  if (n0<1 || n1<1) return 0;
  Bool_t stg = 0;
  Bool_t l0[20]={0};
  Bool_t l1[40]={0};
  Bool_t phi[20]={0};
  for (Int_t i=0;   i< 400; ++i) if (bits.TestBitNumber(i)) l0[      i/20] = 1;
  for (Int_t i=400; i<1200; ++i) if (bits.TestBitNumber(i)) l1[(i-400)/20] = 1;
  for (Int_t i=0; i<20; ++i) {
    if (tolerance) phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40] | l1[(2*i+2)%40] | l1[(2*i+39)%40]);
    else           phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40]);
  }
  for (Int_t dphi=dphiMin;dphi<=dphiMax;dphi++)
    for (Int_t i=0; i<20; ++i) stg |= phi[i] & phi[(i+dphi)%20];
  return stg;
}


//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::RunMC(AliVEvent *fEvent)
{

  // **************************************************************************
  //    TRIGGER INPUTS
  // **************************************************************************
  for(Int_t i=0; i<NTRIGGERINPUTS; i++) fTriggerInputsMC[i] = kFALSE;
  UShort_t fTriggerAD = fEvent->GetADData()->GetTriggerBits();
  UShort_t fTriggerVZERO = fEvent->GetVZEROData()->GetTriggerBits();
  UInt_t fL0inputs = fEvent->GetHeader()->GetL0TriggerInputs();

  fTriggerInputsMC[0] = fTriggerVZERO & (1 << 12); //0VBA VZERO A
  fTriggerInputsMC[1] = fTriggerVZERO & (1 << 13); //0VBC VZERO C
  fTriggerInputsMC[2] = fTriggerAD & (1 << 12);   //0UBA ADA
  fTriggerInputsMC[3] = fTriggerAD & (1 << 13);   //0UBC ADC

  // ------------------------
  // TOF trigger input
  // ------------------------
  const AliTOFHeader *tofH = fEvent->GetTOFHeader();
  fTOFmask = tofH->GetTriggerMask();
  
  Bool_t firedMaxiPhi[36] = {0};
  Int_t NfiredMaxiPads = 0;
 
  for(Int_t ltm=0;ltm<72;ltm++){
    Int_t ip = ltm%36;
    for(Int_t cttm=0;cttm<23;cttm++){
      if(fTOFmask->IsON(ltm,cttm) && gRandom->Rndm(1.0)<hTOFeff->GetBinContent(ltm+1,cttm+1)){
    	firedMaxiPhi[ip] = kTRUE;
    	NfiredMaxiPads++;
    	}
      }
    }
  Bool_t offlineOMU = kFALSE;
  for(Int_t ip=0;ip<36;ip++){
    if(!firedMaxiPhi[ip])continue;
    for(Int_t jp=ip+1;jp<36;jp++){
      if(!firedMaxiPhi[jp])continue;
      Int_t DeSlots = jp-ip;
      Int_t AntiDeSlots = 36 - DeSlots;
      if(DeSlots >= 15 && DeSlots <= 18)offlineOMU = kTRUE;
      else if(AntiDeSlots >= 15 && AntiDeSlots <= 18)offlineOMU = kTRUE;
      }
    }
  if(NfiredMaxiPads>6)offlineOMU = kFALSE;

  if(offlineOMU)fTriggerInputsMC[4] = kTRUE;  //0OMU TOF two hits with topology
  if(NfiredMaxiPads >= 2)fTriggerInputsMC[5] = kTRUE;	//0OM2 TOF two hits	

  // ------------------------
  // STG trigger input
  // ------------------------
  //SPD inputs
  const AliVMultiplicity *mult = fEvent->GetMultiplicity();
  Bool_t vPhiInner[20]; for (Int_t i=0; i<20; ++i) vPhiInner[i]=kFALSE;
  Bool_t vPhiOuter[40]; for (Int_t i=0; i<40; ++i) vPhiOuter[i]=kFALSE;

  Int_t nInner(0), nOuter(0);
  for (Int_t i(0); i<1200; ++i) {
    Bool_t isFired = (mult->TestFastOrFiredChips(i)) && (gRandom->Rndm(1.0) < hSPDeff->GetBinContent(i+1));
    if (i<400) {
      if(isFired)vPhiInner[i/20] = kTRUE;
      nInner += isFired;
    } else {
      if(isFired)vPhiOuter[(i-400)/20] = kTRUE;
      nOuter += isFired;
    }
  }

  // STG input
  Int_t dphiMax=10; 
  Int_t dphiMin=9;
  Bool_t tolerance = 1;
  Bool_t firedSTG = 0;
  Bool_t phi[20]={0};
  if(fRunNumber<295753){ dphiMin = 3; }

  for (Int_t i=0; i<20; ++i) {
    if (tolerance) phi[i] = vPhiInner[i] & (vPhiOuter[(2*i)%40] | vPhiOuter[(2*i+1)%40] | vPhiOuter[(2*i+2)%40] | vPhiOuter[(2*i+39)%40]);
    else           phi[i] = vPhiInner[i] & (vPhiOuter[(2*i)%40] | vPhiOuter[(2*i+1)%40]);
  }
  for (Int_t dphi=dphiMin;dphi<=dphiMax;dphi++)
    for (Int_t i=0; i<20; ++i) firedSTG |= phi[i] & phi[(i+dphi)%20];

  // STP input
  Int_t fired = 0;
  for (Int_t i(0); i<10; ++i) {
    for (Int_t j(0); j<2; ++j) {
      const Int_t k(2*i+j);
      fired += ((   vPhiOuter[k]    || vPhiOuter[k+1]       ||
                    vPhiOuter[k+2]      )
                && (vPhiOuter[k+20] || vPhiOuter[(k+21)%40] ||
                    vPhiOuter[(k+22)%40])
                && (vPhiInner[i]    || vPhiInner[i+1]       )
                && (vPhiInner[i+10] || vPhiInner[(i+11)%20]));
    }
  }
  //0SMB - At least one hit in SPD
  if (nOuter > 0 || nInner > 0) fTriggerInputsMC[6] = kTRUE;
  //0SM2 - Two hits on outer layer
  if (nOuter > 1) fTriggerInputsMC[7] = kTRUE;
  //0STP - Topological SPD trigger (two pairs)
  if (fired != 0) fTriggerInputsMC[8] = kTRUE;
  //0SH1 - More then 6 hits on outer layer
  if (nOuter >= 7) fTriggerInputsMC[9] = kTRUE;
  if(firedSTG) fTriggerInputsMC[10] = kTRUE;

  // **************************************************************************
  //    LOCAL VARIABLES DEFINITION
  // **************************************************************************

  TLorentzVector vGenerated, vDecayProduct;
  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  vGenerated.SetPtEtaPhiM(0.,0.,0.,0.);
  fInEtaGen = kTRUE;

  // **************************************************************************
  //    MC PARTICLES LOOP
  // **************************************************************************
  
  AliMCEvent *mc = MCEvent();
  if(!mc) return;

  for(Int_t imc=0; imc<mc->GetNumberOfTracks(); imc++) {
    AliMCParticle *mcPart = (AliMCParticle*) mc->GetTrack(imc);
    if(!mcPart) continue;
    
    if(TMath::Abs(mcPart->PdgCode()) == 11 || 
       TMath::Abs(mcPart->PdgCode()) == 13 ||
       TMath::Abs(mcPart->PdgCode()) == 2212 ||
       TMath::Abs(mcPart->PdgCode()) == 211||
       TMath::Abs(mcPart->PdgCode()) == 22){
       
       if(mcPart->GetMother() == -1){
         if(TMath::Abs(mcPart->PdgCode()) != 22 && TMath::Abs(mcPart->Eta())>cutEta) fInEtaGen = kFALSE;
    
         TParticlePDG *partGen = pdgdat->GetParticle(mcPart->PdgCode());
         vDecayProduct.SetXYZM(mcPart->Px(),mcPart->Py(), mcPart->Pz(),partGen->Mass());
         vGenerated += vDecayProduct;
         }
       else{
         AliMCParticle *mcMother = (AliMCParticle*) mc->GetTrack(mcPart->GetMother());
         if(TMath::Abs(mcMother->PdgCode()) != 443 && TMath::Abs(mcMother->PdgCode()) != 100443)continue;
         if(TMath::Abs(mcPart->PdgCode()) != 22 && TMath::Abs(mcPart->Eta())>cutEta) fInEtaGen = kFALSE;
    
         TParticlePDG *partGen = pdgdat->GetParticle(mcPart->PdgCode());
         vDecayProduct.SetXYZM(mcPart->Px(),mcPart->Py(), mcPart->Pz(),partGen->Mass());
         vGenerated += vDecayProduct;
       	 }
    }
  }//loop over mc particles 

  FillTree(fTreeGen,vGenerated);

}//RunMC



//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::Terminate(Option_t *) 
{
  cout<<"Analysis complete."<<endl;
}//Terminate

//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::FillTree(TTree *t, TLorentzVector v) {

  fPt      = v.Pt();
  if(v.E() != v.Pz())fY = v.Rapidity();
  else fY = -999;
  fM       = v.M();
  fPhi     = v.Phi();

  t->Fill();

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcNano_MB::FillTree(TTree *t, TLorentzVector v, TLorentzVector vgen) {

  fPt      = v.Pt();
  fPtGen   = vgen.Pt();
  fY       = v.Rapidity();
  fYGen    = vgen.Rapidity();
  fM       = v.M();
  fMGen    = vgen.M();
  fPhi     = v.Phi();
  fPhiGen  = vgen.Phi();

  t->Fill();

}

//_____________________________________________________________________________
Double_t AliAnalysisTaskUpcNano_MB::GetMedian(Double_t *daArray) {
    // Allocate an array of the same size and sort it.
    Double_t dpSorted[4];
    for (Int_t i = 0; i < 4; ++i) {
        dpSorted[i] = daArray[i];
    }
    for (Int_t i = 3; i > 0; --i) {
        for (Int_t j = 0; j < i; ++j) {
            if (dpSorted[j] > dpSorted[j+1]) {
                Double_t dTemp = dpSorted[j];
                dpSorted[j] = dpSorted[j+1];
                dpSorted[j+1] = dTemp;
            }
        }
    }

    // Middle or average of middle values in the sorted array.
    Double_t dMedian = 0.0;
    dMedian = (dpSorted[2] + dpSorted[1])/2.0;
    
    return dMedian;
}
