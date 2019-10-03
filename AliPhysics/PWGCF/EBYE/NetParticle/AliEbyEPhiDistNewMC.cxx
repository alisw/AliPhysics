/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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


//=========================================================================//
//             AliEbyE Analysis for Net-Particle study                     //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                            Surya Prakash Pathak                         //
//                       surya.prakash.pathak@cern.ch                      //
//                         (Last Modified 2019/03/06)                      //
//                 Dealing with Wide pT Window Modified to ESDs            //
//Some parts of the code are taken from J. Thaeder/ M. Weber NetParticle analysis code//
//=========================================================================//

#include <Riostream.h>
#include "TList.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODHeader.h"
#include "AliAODpidUtil.h"
#include "AliEventCuts.h"
#include "AliMultSelection.h"
//#include "AliHelperPID.h"
#include "AliAnalysisTaskSE.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliEbyEPhiDistNewMC.h"

using std::endl;
using std::cout;

ClassImp(AliEbyEPhiDistNewMC)

//-----------------------------------------------------------------------
AliEbyEPhiDistNewMC::AliEbyEPhiDistNewMC() 
: AliAnalysisTaskSE(), 
  fThnList(NULL),
  fArrayMC(0x0),
  fESDtrackCuts(0x0),
  fMCEvent(0x0),
  fMCStack(0x0),
  fEventCuts(NULL),
  fanaUtils(0x0),
  fRun("LHC10h"),
  fCentralityEstimator("V0M"),
  
  fAODtrackCutBit(128),
  fVxMax(3.), 
  fVyMax(3.), 
  fVzMax(10.), 
  fPtMin(0.15),   
  fPtMax(3.15),
  fPhiMin(0.0),
  fPhiMax(6.28),
  fEtaMin(-1.),
  fEtaMax(1.),
  fNptBins(22),
  fNphiBins(18),
  fDcaXy(10.),
  fDcaZ(10.),  
  fNcrossRows(80),
  fChi2NDF(4),

  
  fIsMC(kFALSE),
  fIsAOD(kFALSE),
  fIsQA(kFALSE),
  fIsTrig(kFALSE),
  fIsThn(kFALSE),

  fNTracks(0),
  fCentrality(-1),
  fEventCounter(NULL), 
  fHistCent(0x0), 

  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fPidType(1),
  fMcPid(211),
  fPidStrategy(0),
  fNSigmaMaxITS(3.),
  fNSigmaMaxTPC(3.),
  fNSigmaMaxTPClow(3.),
  fNSigmaMaxTOF(3.),
  fMinPtForTOFRequired(0.69),
  fMaxPtForTPClow(0.69),
  fParticleSpecies(AliPID::kPion),

  fHistTPC(NULL), 
  fHistTOF(NULL),
  fHistITS(NULL),

  fHistTPCc(NULL),
  fHistTOFc(NULL),
  fHistITSc(NULL),

  fHistTPCTOF(NULL),
  fHistTPCTOFc(NULL),
  

  fHistNsTPC(NULL), 
  fHistNsTOF(NULL), 
  fHistNsITS(NULL), 

  fHistNsTPCc(NULL), 
  fHistNsTOFc(NULL), 
  fHistNsITSc(NULL),
  
  fPhiBinNplusNminusCh(NULL),
  fPhiBinNplusNminusChTruth(NULL),
  fTHnCentNplusNminusCh(NULL)

{ 
  
  for (Int_t i = 0; i < 2; i++) {
    fHitCentRec[i] = NULL;
    fHitCentGen[i] = NULL;
    nPidRec[i] = 0.;
    nPidRecP[i] = 0.;
    nPidRecMid[i] = 0.;
    nPidRecSec[i] = 0.;
    nPidRecWD[i] = 0.;
    nPidWoPID[i] = 0.;
    for (Int_t j = 0; j < 3; j++) {
      fHistERec[i][j]     = NULL;
      fHistERecPri[i][j]  = NULL;      
      fHistEGen[i][j]     = NULL;      
      
      fHistCSec[i][j] = NULL;      
      fHistCMat[i][j] = NULL;      
      fHistCMisId[i][j] = NULL;      
      //  fHistCLepton[i][j] = NULL;      
    }
  }
}

//-----------------------------------------------------------------------
AliEbyEPhiDistNewMC::AliEbyEPhiDistNewMC( const char *name ) 
  : AliAnalysisTaskSE(name), 
    fThnList(NULL), 
    fArrayMC(NULL),
    fESDtrackCuts(NULL),
    fMCEvent(NULL),
    fMCStack(NULL),
    fEventCuts(NULL),
    fanaUtils(NULL),
    fRun("LHC10h"),
    fCentralityEstimator("V0M"),
    
    fAODtrackCutBit(128),
    fVxMax(3.), 
    fVyMax(3.), 
    fVzMax(10.), 
    fPtMin(0.15),   
    fPtMax(3.15),
    fPhiMin(0.0),
    fPhiMax(6.28),
    fEtaMin(-1.),
    fEtaMax(1.),
    fNptBins(22),
    fNphiBins(18),
    fDcaXy(10.),
    fDcaZ(10.),
    fNcrossRows(80),
    fChi2NDF(4),
    
    fIsMC(kFALSE),
    fIsAOD(kFALSE),
    fIsQA(kFALSE),
    fIsTrig(kFALSE),
    fIsThn(kFALSE),
    
    fNTracks(0),
    
    fCentrality(-1),
    fEventCounter(NULL), 
    fHistCent(NULL), 
   
    fPIDResponse(0x0),
    fPIDCombined(0x0),
    fPidType(1),
    fMcPid(211),
    fPidStrategy(0),
    fNSigmaMaxITS(4.),
    fNSigmaMaxTPC(4.),
    fNSigmaMaxTPClow(3.),
    fNSigmaMaxTOF(4.),
    fMinPtForTOFRequired(0.69),
    fMaxPtForTPClow(0.69),
    fParticleSpecies(AliPID::kPion),
    fHistTPC(NULL), 
    fHistTOF(NULL),
    fHistITS(NULL),
    
    fHistTPCc(NULL),
    fHistTOFc(NULL),
    fHistITSc(NULL),
    
    fHistTPCTOF(NULL),
    fHistTPCTOFc(NULL),
    
    
    fHistNsTPC(NULL), 
    fHistNsTOF(NULL), 
    fHistNsITS(NULL), 
    
    fHistNsTPCc(NULL), 
    fHistNsTOFc(NULL), 
    fHistNsITSc(NULL),

    fPhiBinNplusNminusCh(NULL),
    fPhiBinNplusNminusChTruth(NULL),
    fTHnCentNplusNminusCh(NULL)
{  
  for (Int_t i = 0; i < 2; i++) {
    fHitCentRec[i] = NULL;
    fHitCentGen[i] = NULL;
    nPidRec[i] = 0.;
    nPidRecP[i] = 0.;
    nPidRecMid[i] = 0.;
    nPidRecSec[i] = 0.;
    nPidRecWD[i] = 0.;
    nPidWoPID[i] = 0.;
    for (Int_t j = 0; j < 3; j++) {
      fHistERec[i][j]     = NULL;
      fHistERecPri[i][j]  = NULL;      
      fHistEGen[i][j]     = NULL;      
      
      fHistCSec[i][j] = NULL;      
      fHistCMat[i][j] = NULL;      
      fHistCMisId[i][j] = NULL;      
      // fHistCLepton[i][j] = NULL;      
    }
  }
  DefineOutput(1, TList::Class()); 
}

//---------------------------------------------------------------------------------
AliEbyEPhiDistNewMC::~AliEbyEPhiDistNewMC() {
  
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fThnList;
    //delete [] fPtArray;
  }
}

//---------------------------------------------------------------------------------
void AliEbyEPhiDistNewMC::UserCreateOutputObjects(){
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  inputHandler->SetNeedField();
  fPIDResponse = inputHandler->GetPIDResponse();

  fanaUtils = new AliAnalysisUtils();
  
  fThnList = new TList();
  fThnList->SetOwner(kTRUE);
  
  if (!fIsAOD) {
    
    if(!fESDtrackCuts) {

      fESDtrackCuts = new AliESDtrackCuts;
      fESDtrackCuts->SetName(Form("NetPesdCut_%d",fPidType));
      fESDtrackCuts->SetEtaRange(fEtaMin, fEtaMax);
      fESDtrackCuts->SetPtRange(0.1, 1e10);
      
      // TPC track
      fESDtrackCuts->SetMinNCrossedRowsTPC( fNcrossRows );// default = 80, sys: 60 and 100 
      fESDtrackCuts->SetMaxChi2PerClusterTPC( fChi2NDF );//default = 4, sys: 3.5, 4.5
      fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);//
      fESDtrackCuts->SetRequireTPCRefit(kTRUE);//
      fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
      
      if (fPidType == 0) {
	fESDtrackCuts->SetMaxDCAToVertexXY( fDcaXy ); //default=2. sys: 1.5, 2.5
      } else {
	if( fDcaXy < 2) fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0156+0.03/pt^1.01"); //tight cut
	else if( fDcaXy > 2) fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0208+0.04/pt^1.01");//losse cut
	else if( fDcaXy == 2) fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.035/pt^1.01");//default
      }
      
      
      fESDtrackCuts->SetRequireITSRefit(kTRUE);//
      fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); // Reason for the structure in phi
      fESDtrackCuts->SetMaxChi2PerClusterITS(36);
      
      fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);//
      fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);//    
      
      // Vertex Related
      fESDtrackCuts->SetMaxDCAToVertexZ( fDcaZ );//
      fESDtrackCuts->SetDCAToVertex2D(kFALSE);
      fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
    } else 
      Printf(" >>>>  User Track Cuts <<<< ");
    
}
  
  CreatePhiHistMC();
  
  fHistCent = new TH1F("fHistCentPid","Centrality", 100, -0.5, 99.5);			 
  fThnList->Add(fHistCent);
  
  fEventCounter = new TH1D("fEventCounter","EventCounter", 20, -0.5, 19.5);
  fThnList->Add(fEventCounter);
  
  PostData(1, fThnList);
  
}

//----------------------------------------------------------------------------------
void AliEbyEPhiDistNewMC::CreatePhiHistMC() {

  if(fIsQA && fPidType != 0){
    
    const Char_t *pidname[4] = {" ", "Pion", "Kaon", "Proton"};
    
    fHistTOF = new TH2F("fHistPID_TOF","TOF Signal;#it{p_{T}} (GeV/#it{c});", 50, 0, 5, 440, -1.1, 1.1);
    fHistTPC = new TH2F("fHistPID_TPC","TPC Signal;#it{p_{T}} (GeV/#it{c});", 50, 0, 5, 1600, -800, 800);
    fHistITS = new TH2F("fHistPID_ITS","ITS Signal;#it{p_{T}} (GeV/#it{c});", 50, 0, 5, 1600, -800, 800);
    
    fThnList->Add(fHistTOF);
    fThnList->Add(fHistTPC);
    fThnList->Add(fHistITS);
    
    fHistTOFc  = new TH2F("fHistPID_TOFcut",Form("TOF Signal (Selected %s);#it{p_{T}} (GeV/#it{c});",pidname[fPidType]),50, 0, 5,  440, -1.1,1.1);
    fHistTPCc  = new TH2F("fHistPID_TPCcut",Form("TPC Signal (Selected %s);#it{p_{T}} (GeV/#it{c});",pidname[fPidType]),50, 0, 5, 1600, -800,800);
    fHistITSc  = new TH2F("fHistPID_ITScut",Form("ITS Signal (Selected %s);#it{p_{T}} (GeV/#it{c});",pidname[fPidType]),50, 0, 5, 1600, -800,800);
    fThnList->Add(fHistTOFc);
    fThnList->Add(fHistTPCc);
    fThnList->Add(fHistITSc);
    
    fHistTPCTOF  = new TH2F("fHistPID_TPCTOF","TPC vs TOF Signal;TPC;TOF",1600, -800,800, 440, -1.1,1.1);
    fHistTPCTOFc = new TH2F("fHistPID_TPCTOFcut",Form("TPC vs TOF Signal (Selected %s);TPC;TOF",pidname[fPidType]),1600, -800,800, 440, -1.1,1.1);
    
    fThnList->Add(fHistTPCTOF);
    fThnList->Add(fHistTPCTOFc);
    
    fHistNsTPC  = new TH2F("fHistNsPID_TPC","TPS N_{#sigma} after Cut;#it{p_{T}} (GeV/#it{c});N_{#sigma}",50,0,5, 400, -40,40);
    fHistNsTOF  = new TH2F("fHistNsPID_TOF","TOF N_{#sigma} after Cut;#it{p_{T}} (GeV/#it{c});N_{#sigma}",50,0,5, 400, -40,40);
    fHistNsITS  = new TH2F("fHistNsPID_ITS","ITS N_{#sigma} after Cut;#it{p_{T}} (GeV/#it{c});N_{#sigma}",50,0,5, 400, -40,40);
    
    fHistNsTPCc = new TH2F("fHistNsPID_TPCcut",Form("TPS N_{#sigma} (Selected %s);#it{p_{T}} (GeV/#it{c});N_{#sigma}",pidname[fPidType]),50,0,5, 400, -40,40);
    fHistNsTOFc = new TH2F("fHistNsPID_TOFcut",Form("TOF N_{#sigma} (Selected %s);#it{p_{T}} (GeV/#it{c});N_{#sigma}",pidname[fPidType]),50,0,5, 400, -40,40);
    fHistNsITSc = new TH2F("fHistNsPID_ITScut",Form("ITS N_{#sigma} (Selected %s);#it{p_{T}} (GeV/#it{c});N_{#sigma}",pidname[fPidType]),50,0,5, 400, -40,40);
    
    fThnList->Add(fHistNsTPC);
    fThnList->Add(fHistNsTOF);
    fThnList->Add(fHistNsITS);
    
    fThnList->Add(fHistNsTPCc);
    fThnList->Add(fHistNsTOFc);
    fThnList->Add(fHistNsITSc);
    
  }
  
  
  const Char_t *fgkHistName[4] = {"Nch", "Npi", "Nka", "Npr"};
  const Char_t *fgkHistLat[2][4]  = {{"N^{-}", "#pi^{-}", "K^{-}", "P^{-}"},
				     {"N^{+}", "#pi^{+}", "K^{+}", "P^{+}"}
  };
  const Char_t *fgkHistCharge[2] = {"Minus", "Plus"};
  
  Int_t ybins[3] = {30, 32, 18};
  Double_t yaxis[2][3] = {{0.15,-0.8,0.0}, {3.15,0.8,6.6}};
  
  const Int_t xNbins = 100;
  Double_t xBinEdge[xNbins+1];
  for(Int_t iBin = 0 ; iBin <= xNbins; iBin++){
    xBinEdge[iBin] = iBin - 0.5;
  }

  //fPtArray = new Double_t[fNptBins+1];
  
    Double_t chPtBins[17] = { 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 1.0, 1.3, 1.6, 2.0, 2.1 };
    Double_t PtArray[25] = { 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55 };
    
    //fphiarray
    
    Double_t pidPhiBins[19] = {0.0, 0.348, 0.697, 1.04, 1.39, 1.74, 2.09, 2.44, 2.79, 3.14, 3.48, 3.83, 4.18, 4.53, 4.88, 5.23, 5.58, 5.93,  6.28};
  
  //if( fPidType == 0) fPtArray = chPtBins;
  //else fPtArray = pidPtBins;
  
  const Char_t *gstName[3] = {"Pt","Eta","Phi"};
  const Char_t *gstLat[3]  = {"p_{T}","#eta","#phi"};
  
  const Char_t *PidCut[4] = {
    "ChHad",
    "TPC", // 0 only TPC
    "TPC",
    "TPC"
  };
  
  
  
  for(Int_t k = 0; k < 2; k++) {//loop over charge type; k = 0 -ve charge, k = 1, +ve ch.
    
    Int_t i = fPidType;
    TString PIDtype = "Charge";
    if (i != 0) PIDtype = Form("%s",PidCut[fPidType]);
    
    fHitCentRec[k] = new TH2F(Form("fHist%s%sRecCent", fgkHistCharge[k], fgkHistName[i]),Form(" Stg %s : Rec%s", PIDtype.Data(), fgkHistName[i]), 100, -0.5, 99.5, 1900, -0.5, 1899.5);
    fHitCentRec[k]->Sumw2();
    fThnList->Add(fHitCentRec[k]);
    
    if(fIsMC){
      fHitCentGen[k] = new TH2F(Form("fHist%s%sGenCent", fgkHistCharge[k], fgkHistName[i]),Form(" Stg %s : Gen%s", PIDtype.Data(), fgkHistName[i]), 100, -0.5, 99.5, 1900, -0.5, 1899.5);
      fHitCentGen[k]->Sumw2();
      fThnList->Add(fHitCentGen[k]);
    }
    
    for(Int_t j = 0; j < 3; j++) {
      
      if( j == 0){ //for pt only--(unequal bin width )
	
	fHistERec[k][j] = new TH2F(Form("fHist%s%s%sRec",gstName[j],fgkHistCharge[k], fgkHistName[i]),Form(" Stg %s : Rec%s : %s ;#it{Bin};%s",PIDtype.Data(), gstLat[j],fgkHistLat[k][i], gstLat[j]), xNbins, xBinEdge, fNptBins, PtArray);
	fHistERec[k][j]->Sumw2();
	fThnList->Add(fHistERec[k][j]);
	
	if(fIsMC){ 
	  fHistERecPri[k][j] = new TH2F(Form("fHist%s%s%sRecPri",gstName[j], fgkHistCharge[k], fgkHistName[i]),Form("Stg- %s : Rec.Primary %s :  %s ;#it{Bin};%s",PIDtype.Data(),gstLat[j],fgkHistLat[k][i], gstLat[j]), xNbins, xBinEdge, fNptBins, PtArray);
	  fHistERecPri[k][j]->Sumw2();
	  fThnList->Add(fHistERecPri[k][j]);
	  
	  fHistEGen[k][j] = new TH2F(Form("fHist%s%s%sGen",gstName[j], fgkHistCharge[k], fgkHistName[i]),Form("Stg- %s : %s : Generated %s ;#it{Bin};%s",PIDtype.Data(),gstLat[j],fgkHistLat[k][i], gstLat[j]), xNbins, xBinEdge, fNptBins, PtArray);
	  fHistEGen[k][j]->Sumw2();
	  fThnList->Add(fHistEGen[k][j]);
	  
	  fHistCSec[k][j] = new TH2F(Form("fHist%s%s%sContSec",gstName[j], fgkHistCharge[k], fgkHistName[i]),Form("Stg- %s : %s : Secondary %s ;#it{Bin};%s",PIDtype.Data(),gstLat[j],fgkHistLat[k][i], gstLat[j]), xNbins, xBinEdge, fNptBins, PtArray);
	  fHistCSec[k][j]->Sumw2();
	  fThnList->Add(fHistCSec[k][j]);
	  
	  fHistCMat[k][j] = new TH2F(Form("fHist%s%s%sContMat",gstName[j], fgkHistCharge[k], fgkHistName[i]),Form("Stg- %s : %s : Material %s ;#it{Bin};%s",PIDtype.Data(),gstLat[j],fgkHistLat[k][i], gstLat[j]), xNbins, xBinEdge, fNptBins, PtArray);
	  fHistCMat[k][j]->Sumw2();
	  fThnList->Add(fHistCMat[k][j]);
	  
	  fHistCMisId[k][j] = new TH2F(Form("fHist%s%s%sContMisId",gstName[j], fgkHistCharge[k], fgkHistName[i]),Form(" Stg- %s : %s : Mis.Id %s ;#it{Bin};%s",PIDtype.Data(),gstLat[j],fgkHistLat[k][i], gstLat[j]), xNbins, xBinEdge, fNptBins, PtArray);
	  fHistCMisId[k][j]->Sumw2();
	  fThnList->Add(fHistCMisId[k][j]);
	  
	}
	
      }//j == 0; pT only
      
      else {
	//*/
	fHistERec[k][j] = new TH2F(Form("fHist%s%s%sRec",gstName[j],fgkHistCharge[k], fgkHistName[i]),Form(" Stg %s : Rec%s : %s ;#it{Bin};%s",PIDtype.Data(), gstLat[j],fgkHistLat[k][i], gstLat[j]), 100, -0.5, 99.5,ybins[j],yaxis[0][j],yaxis[1][j]);
	fHistERec[k][j]->Sumw2();
	fThnList->Add(fHistERec[k][j]);
	
	if(fIsMC){ 
	  fHistERecPri[k][j] = new TH2F(Form("fHist%s%s%sRecPri",gstName[j], fgkHistCharge[k], fgkHistName[i]),Form("Stg- %s : Rec.Primary %s :  %s ;#it{Bin};%s",PIDtype.Data(),gstLat[j],fgkHistLat[k][i], gstLat[j]), 100, -0.5, 99.5,ybins[j],yaxis[0][j],yaxis[1][j]);
	  fHistERecPri[k][j]->Sumw2();
	  fThnList->Add(fHistERecPri[k][j]);
	  
	  fHistEGen[k][j] = new TH2F(Form("fHist%s%s%sGen",gstName[j], fgkHistCharge[k], fgkHistName[i]),Form("Stg- %s : %s : Generated %s ;#it{Bin};%s",PIDtype.Data(),gstLat[j],fgkHistLat[k][i], gstLat[j]), 100, -0.5, 99.5,ybins[j],yaxis[0][j],yaxis[1][j]);
	  fHistEGen[k][j]->Sumw2();
	  fThnList->Add(fHistEGen[k][j]);
	  
	  fHistCSec[k][j] = new TH2F(Form("fHist%s%s%sContSec",gstName[j], fgkHistCharge[k], fgkHistName[i]),Form("Stg- %s : %s : Secondary %s ;#it{Bin};%s",PIDtype.Data(),gstLat[j],fgkHistLat[k][i], gstLat[j]), 100, -0.5, 99.5,ybins[j],yaxis[0][j],yaxis[1][j]);
	  fHistCSec[k][j]->Sumw2();
	  fThnList->Add(fHistCSec[k][j]);
	  
	  fHistCMat[k][j] = new TH2F(Form("fHist%s%s%sContMat",gstName[j], fgkHistCharge[k], fgkHistName[i]),Form("Stg- %s : %s : Material %s ;#it{Bin};%s",PIDtype.Data(),gstLat[j],fgkHistLat[k][i], gstLat[j]), 100, -0.5, 99.5,ybins[j],yaxis[0][j],yaxis[1][j]);
	  fHistCMat[k][j]->Sumw2();
	  fThnList->Add(fHistCMat[k][j]);
	  
	  fHistCMisId[k][j] = new TH2F(Form("fHist%s%s%sContMisId",gstName[j], fgkHistCharge[k], fgkHistName[i]),Form(" Stg- %s : %s : Mis.Id %s ;#it{Bin};%s",PIDtype.Data(),gstLat[j],fgkHistLat[k][i], gstLat[j]), 100, -0.5, 99.5,ybins[j],yaxis[0][j],yaxis[1][j]);
	  fHistCMisId[k][j]->Sumw2();
	  fThnList->Add(fHistCMisId[k][j]);
	}
	
      }//else
      
    }//j ---kinematics
  }//k -----charge
  
  
    if( fPidType == 0){//for charged hadron
        const Int_t dim = 37; //1 centrality bin +  19 phi bins) * 2
        Int_t bin[dim]    = { 100,
            500, 500, 500,
            500, 500, 500, 500, 500, 500, 500, 500,
            200, 200, 200, 200, 200, 200, 200, 200,
            500, 500, 500,
            500, 500, 500, 500, 500, 500, 500, 500,
            200, 200, 200, 200, 200, 200};
        
        Double_t min[dim] = { -0.5,
            -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
            -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
            -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
            -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5};
        
        Double_t max[dim] = { 99.5,
            499.5, 499.5, 499.5,
            499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5,
            199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5,
            499.5, 499.5, 499.5,
            499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5,
            199.5, 199.5, 199.5, 199.5, 199.5, 199.5 };
    
    fPhiBinNplusNminusCh = new THnSparseI("fPhiBinNplusNminusCh","cent-nplus-nminus", dim, bin, min, max);
    fThnList->Add(fPhiBinNplusNminusCh);
    
    if( fIsMC ){
      fPhiBinNplusNminusChTruth = new THnSparseI("fPhiBinNplusNminusChTruth","cent-nplus-nminus", dim, bin, min, max);
      fThnList->Add(fPhiBinNplusNminusChTruth);
    }
    
  }
  else {//for PID
      const Int_t dim = 37; //1 centrality bin +  19 phi bins) * 2
      Int_t bin[dim]    = { 100,
          500, 500, 500,
          500, 500, 500, 500, 500, 500, 500, 500,
          200, 200, 200, 200, 200, 200, 200, 200,
          500, 500, 500,
          500, 500, 500, 500, 500, 500, 500, 500,
          200, 200, 200, 200, 200, 200};
      
      Double_t min[dim] = { -0.5,
          -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
          -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
          -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
          -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5};
      
      Double_t max[dim] = { 99.5,
          499.5, 499.5, 499.5,
          499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5,
          199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5,
          499.5, 499.5, 499.5,
          499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5,
          199.5, 199.5, 199.5, 199.5, 199.5, 199.5 };

    
    fPhiBinNplusNminusCh = new THnSparseI("fPhiBinNplusNminusCh","cent-nplus-nminus", dim, bin, min, max);
    fThnList->Add(fPhiBinNplusNminusCh);
    
    if( fIsMC ){
      fPhiBinNplusNminusChTruth = new THnSparseI("fPhiBinNplusNminusChTruth","cent-nplus-nminus", dim, bin, min, max);
      fThnList->Add(fPhiBinNplusNminusChTruth);
    }

   
    
  } //else PID
  
  if(fIsQA){
    const Int_t nDim = 3;
    Int_t fBinsCh[nDim] = {100, 1900, 1900};
    Double_t fMinCh[nDim] = { -0.5, -0.5, -0.5 };
    Double_t fMaxCh[nDim] = { 99.5, 1899.5, 1899.5};
    fTHnCentNplusNminusCh = new THnSparseI("fTHnCentNplusNminusCh","Cent-NplusChrg-NminusChrg", nDim, fBinsCh, fMinCh, fMaxCh); 
    fTHnCentNplusNminusCh->GetAxis(0)->SetTitle("Centrality");
    fTHnCentNplusNminusCh->GetAxis(1)->SetTitle("Nplus");
    fTHnCentNplusNminusCh->GetAxis(2)->SetTitle("Nminus");
    fThnList->Add(fTHnCentNplusNminusCh);
  }
  //cout <<" Hisotgrams booked " << endl;
  
}


//----------------------------------------------------------------------------------
void AliEbyEPhiDistNewMC::LocalPost() {
  PostData(1, fThnList);

}

//----------------------------------------------------------------------------------
void AliEbyEPhiDistNewMC::UserExec( Option_t * ){

  const Int_t dim = fNphiBins*2;
  Int_t phiCh[dim];
  Int_t phiChMC[dim];
  for(Int_t idx = 0; idx < dim; idx++){
    phiCh[idx] = 0.;
    phiChMC[idx] = 0;
  }
  
  fEventCounter->Fill(1);
  AliVEvent *event = InputEvent(); 
  if (!event) { LocalPost(); return; }

  AliInputEventHandler* fInputEventHandler = static_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fInputEventHandler) return;
 
  //---Check pileup
  if(!fanaUtils){
    cout << "No pileup check ! " << endl;
    return;
  }
  fanaUtils->SetUseMVPlpSelection(kTRUE);
  fanaUtils->SetUseOutOfBunchPileUp(kTRUE);
  
  if(fanaUtils->IsPileUpMV(event)) return;
  if(fanaUtils->IsOutOfBunchPileUp(event)) return;
  if(fanaUtils->IsSPDClusterVsTrackletBG(event)) return;

  const AliVVertex *vertex = event->GetPrimaryVertex();
  if(!vertex) { LocalPost(); return; }

  Bool_t vtest = kFALSE;
  Double32_t fCov[6];
  vertex->GetCovarianceMatrix(fCov);
  if(vertex->GetNContributors() > 0) {
    if(fCov[5] != 0) {
      vtest = kTRUE;
    }
  }
  if(!vtest) { LocalPost(); return; }
  
  if(TMath::Abs(vertex->GetX()) > fVxMax) { LocalPost(); return; }
  if(TMath::Abs(vertex->GetY()) > fVyMax) { LocalPost(); return; }
  if(TMath::Abs(vertex->GetZ()) > fVzMax) { LocalPost(); return; }

  if( fRun == "LHC10h" || fRun == "LHC11h" ){
    AliCentrality *centrality = event->GetCentrality();
    if(!centrality) return;
    if (centrality->GetQuality() != 0) { LocalPost(); return; }
    
    fCentrality = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
  }
  
  else if (fRun == "LHC15o" ){
    AliMultSelection *fMultSelection= (AliMultSelection *) event->FindListObject("MultSelection");
    
    if(!fMultSelection ){
      cout << "AliMultSelection object not found!" << endl;
      return;
    }
    else fCentrality = fMultSelection->GetMultiplicityPercentile(fCentralityEstimator.Data(), false);

  }
 
  if( fCentrality < 0 || fCentrality >=10 ) return;

  fHistCent->Fill(fCentrality);

  fEventCounter->Fill(2);
  
  //---------- Initiate MC
  if (fIsMC){
    fMCEvent = NULL;
    fEventCounter->Fill(8);
      cout << " I am in fill 8" << endl;
    if(fIsAOD) {
      fArrayMC = NULL;
      fArrayMC = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fArrayMC)
	AliFatal("No array of MC particles found !!!");
    }
    else {
      fMCEvent = MCEvent();
      if (!fMCEvent) {
	Printf("ERROR: Could not retrieve MC event");
	LocalPost(); return; 
      }
      
      fMCStack = fMCEvent->Stack();
      if (!fMCStack) {
	Printf("ERROR: Could not retrieve MC stack");
	LocalPost(); return; 
      }
    }
  }//if---fIsMC-----

  //----------
  
  fEventCounter->Fill(3);
  fNTracks  = event->GetNumberOfTracks();

  Int_t iTracks = 0;
  Double_t nRec[2] = {0., 0.};
  Double_t nGen[2] = {0., 0.};
  
  nPidRec[0] = 0.; nPidRec[1] =  0.;
  nPidRecP[0] = 0.; nPidRecP[1] = 0.;
  nPidRecMid[0] = 0.; nPidRecMid[1] = 0.;
  nPidRecSec[0] = 0.; nPidRecSec[1] = 0.;
  nPidRecWD[0] = 0.;  nPidRecWD[1] = 0.;
  nPidWoPID[0] = 0.; nPidWoPID[1] = 0.;
  
  //track loop
  
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    
    AliVTrack *track = static_cast<AliVTrack*>(event->GetTrack(idxTrack)); 
    if(!AcceptTrackL(track)) continue;
    
    Int_t icharge = track->Charge() < 0 ? 0 : 1;    
    Float_t lPt  = (Float_t)track->Pt();
    Float_t lEta = (Float_t)track->Eta();
    Float_t lPhi = (Float_t)track->Phi();
    
    Int_t iphibin = GetPhiBin(lPhi);
    if( iphibin < 0 || iphibin > fNphiBins-1 ) continue;
    Bool_t isPid = kFALSE;
    
    if (fPidType == 0 ) {
      fHistERec[icharge][0]->Fill(fCentrality,lPt); 
      fHistERec[icharge][1]->Fill(fCentrality,lEta); 
      fHistERec[icharge][2]->Fill(fCentrality,lPhi);
      nRec[icharge] += 1.;
      
      if(icharge == 1){
	phiCh[iphibin] += 1;
      }
      if(icharge == 0){
	phiCh[iphibin+fNphiBins] += 1;
      }
      
    }//fPidType == 0
    
    if (fPidType != 0 ) {
      isPid = IsPidPassed(track);
      if(isPid){
	nPidRec[icharge] += 1.;
	fHistERec[icharge][0]->Fill(fCentrality,lPt); 
	fHistERec[icharge][1]->Fill(fCentrality,lEta); 
	fHistERec[icharge][2]->Fill(fCentrality,lPhi);
	nRec[icharge] += 1.;
	
	if(icharge == 1){
	  phiCh[iphibin] += 1;
	}
	if(icharge == 0){
	  phiCh[iphibin+fNphiBins] += 1;
	}
      }
      
    }//PID----

    //---------------------------------
    // - MC Loop for Physical Primary -
    //--------------------------------- 
    if (fIsMC) {
      
      Int_t label  = TMath::Abs(track->GetLabel()); 
      
      Bool_t isPhysicalPrimary        = 0;
      Bool_t isSecondaryFromWeakDecay = 0;
      Bool_t isSecondaryFromMaterial  = 0;
      AliVParticle* particle = NULL;
      
      if (track->InheritsFrom("AliESDtrack")) {
	particle = static_cast<AliVParticle*>(fMCEvent->GetTrack(label));
	if (!particle) return;
	isPhysicalPrimary        = fMCStack->IsPhysicalPrimary(label);
	isSecondaryFromWeakDecay = fMCStack->IsSecondaryFromWeakDecay(label);
	isSecondaryFromMaterial  = fMCStack->IsSecondaryFromMaterial(label);
      }
      else{
	particle                 =  static_cast<AliVParticle*>(fArrayMC->At(label));
	isPhysicalPrimary        =  (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();
	isSecondaryFromWeakDecay =  (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromWeakDecay();
	isSecondaryFromMaterial  =  (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromMaterial();
      }
      

      Float_t fpTRec   = particle->Pt();
      Float_t fEtaRec  = particle->Eta();
      Float_t fPhiRec  = particle->Phi();
      Int_t pdg        = TMath::Abs(particle->PdgCode());
      
      Bool_t isLep = kFALSE;
      if ((pdg == 11) ||(pdg == 13))  isLep = kTRUE; // Is Leptons ?
      
      if (isPhysicalPrimary) {
	if (fPidType == 0) { // For Charge
	  if (!isLep) {
	    nPidRecP[icharge] += 1.;
	    fHistERecPri[icharge][0]->Fill(fCentrality,fpTRec); 
	    fHistERecPri[icharge][1]->Fill(fCentrality,fEtaRec); 
	    fHistERecPri[icharge][2]->Fill(fCentrality,fPhiRec); 
	  } else {
	    nPidRecMid[icharge] += 1.;
	    fHistCMisId[icharge][0]->Fill(fCentrality,fpTRec); 
	    fHistCMisId[icharge][1]->Fill(fCentrality,fEtaRec); 
	    fHistCMisId[icharge][2]->Fill(fCentrality,fPhiRec); 
	  }//mis-Id
	}//fPidType == 0---ch. particle
      }//Phys. Primary

      //PID------
      if(isPid){
	if(isPhysicalPrimary){
	  if ( pdg == fMcPid ) {// For PID
	    nPidRecP[icharge] += 1.;
	    fHistERecPri[icharge][0]->Fill(fCentrality,fpTRec); 
	    fHistERecPri[icharge][1]->Fill(fCentrality,fEtaRec); 
	    fHistERecPri[icharge][2]->Fill(fCentrality,fPhiRec); 
	  }
	  else {
	    nPidRecMid[icharge] += 1.;
	    fHistCMisId[icharge][0]->Fill(fCentrality,fpTRec); 
	    fHistCMisId[icharge][1]->Fill(fCentrality,fEtaRec); 
	    fHistCMisId[icharge][2]->Fill(fCentrality,fPhiRec); 
	  }
	}//------------isPhysics Primary -------------
	else if (isSecondaryFromWeakDecay) {
	  nPidRecWD[icharge] += 1.;
	  fHistCSec[icharge][0]->Fill(fCentrality,fpTRec); 
	  fHistCSec[icharge][1]->Fill(fCentrality,fEtaRec); 
	  fHistCSec[icharge][2]->Fill(fCentrality,fPhiRec); 
	}//------isSecondaryFromWeakDecay--------------
	else if (isSecondaryFromMaterial) {
	  nPidRecSec[icharge] += 1.;
	  fHistCMat[icharge][0]->Fill(fCentrality,fpTRec); 
	  fHistCMat[icharge][1]->Fill(fCentrality,fEtaRec); 
	  fHistCMat[icharge][2]->Fill(fCentrality,fPhiRec); 
	}//------isSecondaryFromMaterial--------------------
	
      }// if(IsPid)---PID finish
      
    }//IsMC
    
    iTracks++; //reco track numbers
  }//rec track loop --

  fHitCentRec[0]->Fill(fCentrality,nRec[0]);
  fHitCentRec[1]->Fill(fCentrality,nRec[1]);

  
  //cout << " Rec positve " << nRec[1] <<" and -ve particle " << nRec[0] << endl;
  
  
  Double_t phiContainer[dim+1];
  phiContainer[0] = (Double_t)fCentrality;
  for(Int_t i = 1; i <= dim; i++){
    phiContainer[i] = phiCh[i-1];
  }
  
  fPhiBinNplusNminusCh->Fill(phiContainer);
  
    
  fEventCounter->Fill(7);
  //---- - -- - - - - -   -  -- - - - ---- - - - ---
  if (fIsMC) {
    fEventCounter->Fill(8);
    if (fIsAOD) {
      for (Int_t idxMC = 0; idxMC < fArrayMC->GetEntries(); idxMC++) {
	AliAODMCParticle *particle = static_cast<AliAODMCParticle*>(fArrayMC->At(idxMC));
	if (!particle) 
	  continue;

	if (!particle->IsPhysicalPrimary()) continue;
	if (!AcceptTrackLMC((AliVParticle*)particle)) continue;
	Int_t icharge = (particle->PdgCode() < 0) ? 0 : 1;

	Float_t fpTGen   = particle->Pt();
	Float_t fEtaGen  = particle->Eta();
	Float_t fPhiGen  = particle->Phi();

	Int_t pdg  = TMath::Abs(particle->PdgCode());
	
	if (fPidType != 0) { // if not charge and if not pdg = pidtype
	  if (pdg != fMcPid) continue;
	} else { // if charge and if lepton 
	  if (pdg == 11 || pdg == 13) continue;
	}
	
	fHistEGen[icharge][0]->Fill(fCentrality,fpTGen); 
	fHistEGen[icharge][1]->Fill(fCentrality,fEtaGen); 
	fHistEGen[icharge][2]->Fill(fCentrality,fPhiGen);
	nGen[icharge] += 1.;

	Int_t iphibinMC = GetPhiBin(fPhiGen);
	if( iphibinMC < 0 || iphibinMC > fNphiBins-1 ) continue;
	
	if(icharge == 1){
	  phiChMC[iphibinMC] += 1;
	}
	if(icharge == 0){
	  phiChMC[iphibinMC+fNphiBins] += 1;
	}
	
      } //AOD track
      fEventCounter->Fill(9);
    }//if(AOD)----
    else{
      fEventCounter->Fill(10);
      for (Int_t idxMC = 0; idxMC < fMCStack->GetNprimary(); ++idxMC) {
	AliVParticle* particle = fMCEvent->GetTrack(idxMC);
	if (!particle) 
	  continue;
	if(!fMCStack->IsPhysicalPrimary(idxMC))  continue;
	
	if (!AcceptTrackLMC(particle)) continue;
	Int_t icharge = (particle->PdgCode() < 0) ? 0 : 1;
	
	Float_t fpTGen   = particle->Pt();
	Float_t fEtaGen  = particle->Eta();
	Float_t fPhiGen  = particle->Phi();
	Int_t pdg = TMath::Abs(particle->PdgCode());
	
	if (fPidType != 0) {
	  if (pdg != fMcPid) continue;
	} else {
	  if (pdg == 11 || pdg == 13) continue;
	}
	
	fHistEGen[icharge][0]->Fill(fCentrality,fpTGen); 
	fHistEGen[icharge][1]->Fill(fCentrality,fEtaGen); 
	fHistEGen[icharge][2]->Fill(fCentrality,fPhiGen); 
	nGen[icharge] += 1.;
	
	Int_t iphibinMC = GetPhiBin(fPhiGen);
	if( iphibinMC < 0 || iphibinMC > fNphiBins-1 ) continue;
	
	if(icharge == 1){
	  phiChMC[iphibinMC] += 1;
	}
	if(icharge == 0){
	  phiChMC[iphibinMC+fNphiBins] += 1;
	}	
      }
      
      
      fEventCounter->Fill(11);
    }//else (ESD)
    
    fHitCentGen[0]->Fill(fCentrality, nGen[0]);
    fHitCentGen[1]->Fill(fCentrality, nGen[1]);

    
    
    //cout << " Gen positve " << nGen[1] <<" and -ve particle " << nGen[0] << endl;
   
    if(fIsQA){
      Double_t fContainerCh[3] = { (double)fCentrality, nGen[1], nGen[0] };
      fTHnCentNplusNminusCh->Fill(fContainerCh);
    }
    
    Double_t phiContainerMC[dim+1];
    phiContainerMC[0] = (Double_t) fCentrality;
    for(Int_t i = 1; i <= dim; i++){
      phiContainerMC[i] = phiChMC[i-1];
    }
    
    fPhiBinNplusNminusChTruth->Fill(phiContainerMC);
    
  }
  
  fEventCounter->Fill(12);
  PostData(1, fThnList);
}

//___________________________________________________________
Bool_t AliEbyEPhiDistNewMC::AcceptTrackL(AliVTrack *track) const {
  
  if (!track) return kFALSE; 
  if (track->Charge() == 0) return kFALSE; 
    
  if (fIsAOD) {  // AOD
    AliAODTrack * trackAOD = dynamic_cast<AliAODTrack*>(track);
    if (!trackAOD) {
      AliError("Pointer to dynamic_cast<AliAODTrack*>(track) = ZERO");
      return kFALSE; 
    }
    if (!trackAOD->TestFilterBit(fAODtrackCutBit))
      return kFALSE;
  } else {      // ESDs
    if(!fESDtrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))  return kFALSE;
  }
  
    Double_t ptot = track->P();
  
  if(track->Pt() < fPtMin || track->Pt() > fPtMax )  return kFALSE;
    if (track->Phi() <fPhiMin || track->Phi() >fPhiMax) return kFALSE;


  Double_t partMass = AliPID::ParticleMass(fParticleSpecies);
  Double_t pz = track->Pz();
  Double_t en = TMath::Sqrt( ptot*ptot + partMass*partMass );
  Double_t rap = -999.;
  if( en != TMath::Abs(pz) ){
    rap = 0.5*TMath::Log( (en + pz)/(en - pz) );
  }
  else rap = -999.;
  
  if( TMath::Abs(rap) > 0.5 ) return kFALSE;//rapidity cut
  
  if (TMath::Abs(track->Eta()) > fEtaMax) return kFALSE; 

  return kTRUE;
}


//___________________________________________________________
Bool_t AliEbyEPhiDistNewMC::AcceptTrackLMC(AliVParticle *particle) const {
  if(!particle) return kFALSE;
  if (particle->Charge() == 0.0) return kFALSE;
    
    Double_t ptotMC = particle->P();
  
  if (particle->Pt() < fPtMin || particle->Pt() > fPtMax) return kFALSE;
    if (particle->Phi() < fPhiMin || particle->Phi() > fPhiMax) return kFALSE;
    
  //rapidity cut
  Double_t partMass = AliPID::ParticleMass(fParticleSpecies);
  Double_t pz = particle->Pz();
  Double_t en = TMath::Sqrt( ptotMC*ptotMC + partMass*partMass );
  Double_t rap;
  if( en != TMath::Abs(pz) ){
    rap = 0.5*TMath::Log( (en + pz)/(en - pz) );
  }
  else rap = -999;
  
  if( TMath::Abs(rap) > 0.5 ) return kFALSE;//rapidity cut
  if (TMath::Abs(particle->Eta()) > fEtaMax) return kFALSE;
  
  return kTRUE;
}
//---------------------------------------------------------------
Int_t AliEbyEPhiDistNewMC::GetPtBin(Double_t pt){
  
  Int_t bin = -1;
  
  if( fPidType == 0){
    Double_t chPtBins[17] = { 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 1.0, 1.3, 1.6, 2.0, 2.1 };
    for(Int_t iBin = 0; iBin < fNptBins; iBin++){
      
      if( iBin == fNptBins-1){
	if( pt >= chPtBins[iBin] && pt <= chPtBins[iBin+1]){
	  bin = iBin;
	  break;
	}
      }
      else{
	if( pt >= chPtBins[iBin] && pt < chPtBins[iBin+1]){
	  bin = iBin;
	  break;
	  
	}
      }
    }//for 
  }
  else {
      Double_t pidPtBins[25] = { 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55 };

    for(Int_t iBin = 0; iBin < fNptBins; iBin++){
      
      if( iBin == fNptBins-1){
	if( pt >= pidPtBins[iBin] && pt <= pidPtBins[iBin+1]){
	  bin = iBin;
	  break;
	}
      }
      else{
	if( pt >= pidPtBins[iBin] && pt < pidPtBins[iBin+1]){
	  bin = iBin;
	  break;
	  
	}
      }
    }//for 
  }
  return bin;
  
}

//-----------------------------------------GetPhiBin--------------------------------

Int_t AliEbyEPhiDistNewMC::GetPhiBin(Double_t Phi){
    
    Int_t phibin = -1;
    
    Double_t pidPhiBins[19] = {0.0, 0.348, 0.697, 1.04, 1.39, 1.74, 2.09, 2.44, 2.79, 3.14, 3.48, 3.83, 4.18, 4.53, 4.88, 5.23, 5.58, 5.93,  6.28};
    
    for (Int_t pBin = 0; pBin < fNphiBins; pBin++){
        
        if (pBin == fNphiBins - 1){
            if (Phi >= pidPhiBins[pBin] && Phi <= pidPhiBins[pBin+1]){
                phibin = pBin ;
                break;
            }
        }
        else {
            if( Phi >= pidPhiBins[pBin] && Phi < pidPhiBins[pBin +1]){
                phibin = pBin;
                break;
            }
        }
        
    }
    return phibin;
}

//___________________________________________________________
//___________________________________________________________
void AliEbyEPhiDistNewMC::Terminate( Option_t * ){
  //
  //
  Printf(" ------------------------------------------\n"
	 "             Teminating the task           \n"
	 "--------------------------------------------\n");
  
  //
}

void AliEbyEPhiDistNewMC::SetPidType(Int_t i) {
  fPidType = i; fMcPid = -9999;
  if (fPidType == 1) { fParticleSpecies = AliPID::kPion; fMcPid = 211; }
  if (fPidType == 2) { fParticleSpecies = AliPID::kKaon; fMcPid = 321;}
  if (fPidType == 3) { fParticleSpecies = AliPID::kProton; fMcPid = 2212;}
 }

Double_t AliEbyEPhiDistNewMC::TOFBetaCalc(AliVTrack *track) const{
  //TOF beta calculation
  Double_t tofTime=track->GetTOFsignal();
  
  Double_t c=TMath::C()*1.E-9;// m/ns
  Float_t startTime = fPIDResponse->GetTOFResponse().GetStartTime(((AliVTrack*)track)->P());//in ps
  Double_t length= fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron)*1E-3*c;
  tofTime -= startTime;      // subtract startTime to the signal
  Double_t tof= tofTime*1E-3; // ns, average T0 fill subtracted, no info from T0detector 	 
  tof=tof*c;
  return length/tof;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Double_t AliEbyEPhiDistNewMC::GetMass(AliPID::EParticleType id) const{
  //return Mass according to AliHelperParticleSpecies_t. If undefined return -999.
  Double_t mass=-999.;
  if (id == AliPID::kProton) { mass=9.38271999999999995e-01; }
  if (id == AliPID::kKaon)   { mass=4.93676999999999977e-01; }
  if (id == AliPID::kPion)   { mass=1.39570000000000000e-01; }
  return mass;
}


//___________________________________________________________
Bool_t AliEbyEPhiDistNewMC::IsPidPassed(AliVTrack * track) {
  //return  kTRUE; - PID strategy is from Jochen, Hans and Deepika

  Bool_t isAcceptedITS    = kFALSE;
  Bool_t isAcceptedTPC    = kFALSE;
  Bool_t isAcceptedTPClow = kFALSE;
  Bool_t isAcceptedTOF    = kFALSE;
  Bool_t isAccepted       = kFALSE;
  Bool_t hasPIDTOF = kFALSE;
  
  Double_t *pid = new Double_t[3];
  pid[0] = 10.;
  pid[1] = 10.;
  pid[2] = 10.;
  
  Double_t pt = track->Pt();
  
  //---------------------------| el, mu,  pi,  k,    p   | Pt cut offs from spectra
    //ITS--------------
    Double_t ptLowITS[5]       = { 0., 0., 0.2,  0.2,  0.2  };
    Double_t ptHighITS[5]      = { 0., 0., 0.2,  0.2,  0.2  };
    //TPC---------------
    Double_t ptLowTPC[5]       = { 0., 0., 0.3,  0.3, 0.3  };
    Double_t ptHighTPC[5]      = { 0., 0., 1.55,  1.55,   1.55  };
    //TOF----
    Double_t ptLowTOF[5]       = { 0., 0., 0.2,  0.2,  0.2  };
    Double_t ptHighTOF[5]      = { 0., 0., 0.2,  0.2,    0.2  };
    //TPCTOF----------
    Double_t ptLowTPCTOF[5]    = { 0., 0., 0.65, 0.69,   0.69  };
    Double_t ptHighTPCTOF[5]   = { 0., 0., 0.65,  0.69,   0.69  };
    
  

  //--------------------------TPC PID----------------
  if (fPIDResponse->CheckPIDStatus((AliPIDResponse::EDetector)AliPIDResponse::kTPC, track) == AliPIDResponse::kDetPidOk) {
    pid[1] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)AliPIDResponse::kTPC, track, fParticleSpecies);
    
      if (TMath::Abs(pid[1]) < fNSigmaMaxTPC)  // Anywhere withing Max Nsigma TPC
	isAcceptedTPC = kTRUE;
      
      if (TMath::Abs(pid[1]) < fNSigmaMaxTPClow)  // Anywhere withing Mim Nsigma TPC
	isAcceptedTPClow = kTRUE;                 // extra identifier
      
      if (track->Pt() < fMaxPtForTPClow)         // if less than a pt when low nsigma should be applied
	isAcceptedTPC = isAcceptedTPClow;        // --------nslow-----|ptlow|----------------nshigh------
      
    Double_t nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron));
    if (TMath::Abs(pid[1]) > nSigma) isAcceptedTPC = kFALSE;
  }
  
  //--------------------------------ITS PID--------------------------
  if (fPIDResponse->CheckPIDStatus((AliPIDResponse::EDetector)AliPIDResponse::kITS, track) == AliPIDResponse::kDetPidOk) {
    pid[0] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)AliPIDResponse::kITS, track, fParticleSpecies);
    if (TMath::Abs(pid[0]) < fNSigmaMaxITS) 
      isAcceptedITS = kTRUE;
    
    Double_t nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kElectron));
    if (TMath::Abs(pid[0]) > nSigma)
      isAcceptedITS = kFALSE;
  }
  
  
  //-------------------------------TOF--------------------------
  if ( fPIDResponse->CheckPIDStatus((AliPIDResponse::EDetector)AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk) {
    pid[2] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)AliPIDResponse::kTOF, track, fParticleSpecies);
    hasPIDTOF = kTRUE;
    if (TMath::Abs(pid[2]) < fNSigmaMaxTOF) 
      isAcceptedTOF = kTRUE;
    
    Double_t nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kElectron));
    if (TMath::Abs(pid[2]) > nSigma)
      isAcceptedTOF = kFALSE;
  }

  
  if (fIsMC && isAcceptedTOF) {
    Int_t tofLabel[3];  
    if (track->InheritsFrom("AliESDtrack")) {
      (dynamic_cast<AliESDtrack*>(track))->GetTOFLabel(tofLabel);
    } else if (track->InheritsFrom("AliAODTrack")) {
      (dynamic_cast<AliAODTrack*>(track))->GetTOFLabel(tofLabel);
    }
    
    Bool_t hasMatchTOF = kTRUE;
    if (TMath::Abs(track->GetLabel()) != TMath::Abs(tofLabel[0]) || tofLabel[1] > 0) {
      hasMatchTOF = kFALSE;
    }
    
    if(fIsAOD) {
      // AliAODMCParticle  *matchedTrack = dynamic_cast<AliAODMCParticle*>(fArrayMC->At(TMath::Abs(tofLabel[0])));
      // if (TMath::Abs(matchedTrack->GetMother()) == TMath::Abs(track->GetLabel())) 
      //	hasMatchTOF = kTRUE;
      // if (TMath::Abs(track->GetLabel()) ==  TMath::Abs(matchedTrack->GetMother()))
      //	cout << mcMother->GetMother() << "  " << tofLabel[0]  << "  " << track->GetLabel() << endl;
    } else {
      TParticle *matchedTrack = fMCStack->Particle(TMath::Abs(tofLabel[0]));
      if (TMath::Abs(matchedTrack->GetFirstMother()) == TMath::Abs(track->GetLabel())) 
	hasMatchTOF = kTRUE;
    }
    isAcceptedTOF = hasMatchTOF;
  }
  
  
  Short_t c    = 0;
  Double_t p   = 0;
  Double_t tpc = 0;
  Double_t tof = 0;
  Double_t its = 0;
  
  if (fIsQA && (fCentrality == 0 || fCentrality == 1)) {
    c   = track->Charge();
    p   = track->P();
    tpc = track->GetTPCsignal();
    tof = TOFBetaCalc(track);    // => GetTOFsignal();
    its = track->GetITSsignal();
    
    fHistTOF->Fill(p,c*tof);
    fHistTPC->Fill(p,c*tpc);
    fHistITS->Fill(p,c*its);
    fHistTPCTOF->Fill(c*tpc,c*tof);
    
    fHistNsITS->Fill(p,pid[0]);
    fHistNsTPC->Fill(p,pid[1]);
    fHistNsTOF->Fill(p,pid[2]);
  }
  
  
  if (fParticleSpecies == 2) {//for Pion: TPC
    isAccepted = isAcceptedTPC;;
  }
  
  if( fParticleSpecies == 3){  //for kaon: TPC
     isAccepted = isAcceptedTPC;
  }
  
  if( fParticleSpecies == 4){//for proton
    
    isAccepted = isAcceptedTPC;
  }//for proton
  
  
  if (fIsQA && isAccepted && (fCentrality == 0 || fCentrality == 1)) {
    fHistTOFc->Fill(p,c*tof);
    fHistTPCc->Fill(p,c*tpc);
    fHistITSc->Fill(p,c*its);
    
    fHistTPCTOFc->Fill(c*tpc,c*tof);
    
    fHistNsITSc->Fill(p,pid[0]);
    fHistNsTPCc->Fill(p,pid[1]);
    fHistNsTOFc->Fill(p,pid[2]);
  }
  
  delete [] pid;
  return isAccepted;
  
  PostData(1, fThnList);
  
}//IsPidPassed


//----------------------------------------------------------------------------------

