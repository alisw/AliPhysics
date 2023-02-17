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

 /* $Id: AliAnalysisTaskRecEff.cxx  Rihan Haque, 18/09/2019 (ver1) $ */

 //-- general include---
 #include "TChain.h"

 #include "TTree.h"

 #include "TGrid.h"

 #include "TROOT.h"

 #include "TArrayI.h"

 #include "TObjArray.h"

 #include "TMatrixDSym.h"

 #include "TParticle.h"

 #include "TMath.h"

 #include "stdio.h"

 #include "Riostream.h"


 //---- manager and handler---
 #include "AliAnalysisManager.h"

 #include "AliInputEventHandler.h"

 //---V0 and ZDC info---
 #include "AliAODZDC.h"

 #include "AliAODVZERO.h"

 #include "AliAODVertex.h"

 //---AOD,ESD event--
 #include "AliESDEvent.h"

 #include "AliAODHeader.h"

 #include "AliAODEvent.h"

 #include "AliAODTrack.h"

 //----- For MC event------
 #include "AliMCEvent.h"

 #include "AliStack.h"

 #include "AliAODMCParticle.h"

 #include "AliAODMCHeader.h"

 //----for PID-----
 #include "AliPIDResponse.h"

 #include "AliPIDCombined.h"

 //----- Vevent and tracks
 #include "AliVEventHandler.h"

 #include "AliVEvent.h"

 #include "AliVTrack.h"

 #include "AliVParticle.h"

 #include "AliCentrality.h"

 //----- must include-------
 #include "AliMultSelection.h"

 #include "AliAnalysisUtils.h"

 #include "AliPhysicsSelection.h"

 #include "AliAnalysisTaskRecEff.h"

 #include "AliAnalysisTaskSE.h"

 #include "AliAnalysisTaskRecEff.h"

 using std::cout;
 using std::endl;
 using std::vector;

 ClassImp(AliAnalysisTaskRecEff)

 AliAnalysisTaskRecEff::AliAnalysisTaskRecEff(const char * name): AliAnalysisTaskSE(name),
 fVevent(NULL),
 fESD(NULL),
 fAOD(NULL),
 fPIDResponse(NULL),
 fMultSelection(NULL),
 fAnalysisUtil(NULL),
 fListHist(NULL),
 fCentralityMin(0),
 fCentralityMax(90),
 fFilterBit(1),
 fTPCclustMin(70),
 bUseKinkTracks(kFALSE),
 fNSigmaCut(2.0),
 fMinPtCut(0.2),
 fMaxPtCut(5.0),
 fDCAxyMax(2.4),
 fDCAzMax(3.2),
 fMinEtaCut(-0.8),
 fMaxEtaCut(0.8),
 fTrkChi2Min(0.1),
 fdEdxMin(10.0),
 fMinVzCut(-10.0),
 fMaxVzCut(10.0),
 sCentrEstimator("V0M"),
 fCentDistBeforCut(NULL),
 fCentDistAfterCut(NULL),
 fHistEtaPtBeforCut(NULL),
 fHistEtaPhiBeforCut(NULL),
 fHistEtaPhiAfterCut(NULL),
 fHistPtDistMBwithCut(NULL),

 fHistPtRecChrgPos(NULL),
 fHistPtGenChrgPos(NULL),
 fHistPtRecPionPos(NULL),
 fHistPtGenPionPos(NULL),
 fHistPtRecKaonPos(NULL),
 fHistPtGenKaonPos(NULL),
 fHistPtRecProtPos(NULL),
 fHistPtGenProtPos(NULL),
 fHistPtRecChrgNeg(NULL),
 fHistPtGenChrgNeg(NULL),
 fHistPtRecPionNeg(NULL),
 fHistPtGenPionNeg(NULL),
 fHistPtRecKaonNeg(NULL),
 fHistPtGenKaonNeg(NULL),
 fHistPtRecProtNeg(NULL),
 fHistPtGenProtNeg(NULL),

 fHistPtRecPionPosimp(NULL),
 fHistPtRecKaonPosimp(NULL),
 fHistPtRecProtPosimp(NULL),
 fHistPtRecPionNegimp(NULL),
 fHistPtRecKaonNegimp(NULL),
 fHistPtRecProtNegimp(NULL),

 fHistnSigmaTPCPionCent1(NULL),
 fHistnSigmaTPCKaonCent1(NULL),
 fHistnSigmaTPCProtCent1(NULL),
 fHistnSigmaTOFPionCent1(NULL),
 fHistnSigmaTOFKaonCent1(NULL),
 fHistnSigmaTOFProtCent1(NULL),

 fHistPtContaminationPrimariesChrgPos(NULL),
 fHistPtContaminationPrimariesChrgNeg(NULL),
 fHistPtContaminationSecondariesChrgPos(NULL),
 fHistPtContaminationSecondariesChrgNeg(NULL),
 fHistPtContaminationSecondariesMaterialChrgPos(NULL),
 fHistPtContaminationSecondariesMaterialChrgNeg(NULL),
 fHistPtContaminationSecondariesWeakDecChrgPos(NULL),
 fHistPtContaminationSecondariesWeakDecChrgNeg(NULL),
 
 fMCVzDistBeforCut(NULL),
 fParticlePdgCode(0),
 fHistAChrggenVsCent(NULL),
 fHistAChrgrecVsCent(NULL),

 fHistTPConlyVsCL1Before(NULL),
 fHistTPConlyVsV0MBefore(NULL),
 fHistCL0VsV0MBefore(NULL),
 fHistTPConlyVsCL1After(NULL),
 fHistTPConlyVsV0MAfter(NULL),
 fHistCL0VsV0MAfter(NULL),
 fHistGlobalVsV0MBefore(NULL),
 fHistGlobalVsV0MAfter(NULL),
 fTPCvsGlobalTrkBefore(NULL),
 fTPCvsGlobalTrkAfter(NULL),
 fHistTPCVsESDTrkBefore(NULL),
 fHistTPCVsESDTrkAfter(NULL),
 fHistEventCount(NULL) 
 {

	 for (int i = 0; i < 10; i++) {
		 fHistAch[i] = NULL;

	 }

	 //Must be here:
	 DefineInput(0, TChain::Class());
	 DefineOutput(1, TList::Class());
 }

 //_______________________empty constructor_______________________
AliAnalysisTaskRecEff::AliAnalysisTaskRecEff():
 fVevent(NULL),
 fESD(NULL),
 fAOD(NULL),
 fPIDResponse(NULL),
 fMultSelection(NULL),
 fAnalysisUtil(NULL),
 fListHist(NULL),
 fCentralityMin(0),
 fCentralityMax(90),
 fFilterBit(1),
 fTPCclustMin(70),
 bUseKinkTracks(kFALSE),
 fNSigmaCut(2.0),
 fMinPtCut(0.2),
 fMaxPtCut(5.0),
 fDCAxyMax(2.4),
 fDCAzMax(3.2),
 fMinEtaCut(-0.8),
 fMaxEtaCut(0.8),
 fTrkChi2Min(0.1),
 fdEdxMin(10.0),
 fMinVzCut(-10.0),
 fMaxVzCut(10.0),
 sCentrEstimator("V0M"),
 fCentDistBeforCut(NULL),
 fCentDistAfterCut(NULL),
 fHistEtaPtBeforCut(NULL),
 fHistEtaPhiBeforCut(NULL),
 fHistEtaPhiAfterCut(NULL),
 fHistPtDistMBwithCut(NULL),

 fHistPtRecChrgPos(NULL),
 fHistPtGenChrgPos(NULL),
 fHistPtRecPionPos(NULL),
 fHistPtGenPionPos(NULL),
 fHistPtRecKaonPos(NULL),
 fHistPtGenKaonPos(NULL),
 fHistPtRecProtPos(NULL),
 fHistPtGenProtPos(NULL),
 fHistPtRecChrgNeg(NULL),
 fHistPtGenChrgNeg(NULL),
 fHistPtRecPionNeg(NULL),
 fHistPtGenPionNeg(NULL),
 fHistPtRecKaonNeg(NULL),
 fHistPtGenKaonNeg(NULL),
 fHistPtRecProtNeg(NULL),
 fHistPtGenProtNeg(NULL),

 fHistPtRecPionPosimp(NULL),
 fHistPtRecKaonPosimp(NULL),
 fHistPtRecProtPosimp(NULL),
 fHistPtRecPionNegimp(NULL),
 fHistPtRecKaonNegimp(NULL),
 fHistPtRecProtNegimp(NULL),

 fHistnSigmaTPCPionCent1(NULL),
 fHistnSigmaTPCKaonCent1(NULL),
 fHistnSigmaTPCProtCent1(NULL),
 fHistnSigmaTOFPionCent1(NULL),
 fHistnSigmaTOFKaonCent1(NULL),
 fHistnSigmaTOFProtCent1(NULL),

 fHistPtContaminationPrimariesChrgPos(NULL),
 fHistPtContaminationPrimariesChrgNeg(NULL),
 fHistPtContaminationSecondariesChrgPos(NULL),
 fHistPtContaminationSecondariesChrgNeg(NULL),
 fHistPtContaminationSecondariesMaterialChrgPos(NULL),
 fHistPtContaminationSecondariesMaterialChrgNeg(NULL),
 fHistPtContaminationSecondariesWeakDecChrgPos(NULL),
 fHistPtContaminationSecondariesWeakDecChrgNeg(NULL),
  
 fMCVzDistBeforCut(NULL),
 fParticlePdgCode(0),
 fHistAChrggenVsCent(NULL),
 fHistAChrgrecVsCent(NULL),

 fHistTPConlyVsCL1Before(NULL),
 fHistTPConlyVsV0MBefore(NULL),
 fHistCL0VsV0MBefore(NULL),
 fHistTPConlyVsCL1After(NULL),
 fHistTPConlyVsV0MAfter(NULL),
 fHistCL0VsV0MAfter(NULL),
 fHistGlobalVsV0MBefore(NULL),
 fHistGlobalVsV0MAfter(NULL),
 fTPCvsGlobalTrkBefore(NULL),
 fTPCvsGlobalTrkAfter(NULL),
 fHistTPCVsESDTrkBefore(NULL),
 fHistTPCVsESDTrkAfter(NULL),
 fHistEventCount(NULL) 
 {

	 for (int i = 0; i < 10; i++) {
		 fHistAch[i] = NULL;
	 }

 }

 //__________________ destructor ___________________
 AliAnalysisTaskRecEff::~AliAnalysisTaskRecEff() {
     if (fListHist) delete fListHist;
     if (fAnalysisUtil) delete fAnalysisUtil; // because its 'new' !!
 }

 //________________ Define Histograms _______________
 void AliAnalysisTaskRecEff::UserCreateOutputObjects() {
     //std::cout<<"\n UserCreateOutputObject: function begins...\n"<<endl; 
     //Get The Input Hander:
     AliAnalysisManager * mgr = AliAnalysisManager::GetAnalysisManager();
     AliInputEventHandler * inputHandler = dynamic_cast < AliInputEventHandler * > (mgr->GetInputEventHandler());

     if (!inputHandler) {
         printf("\n***** ERROR *****\n Input handler missing, Status:QUIT!\n");
         exit(1);
     }

     //PileUp Multi-Vertex
     fAnalysisUtil = new AliAnalysisUtils();
     fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
     fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);

     //// obtain the PID response object if needed:
     fPIDResponse = inputHandler->GetPIDResponse();
     if (!fPIDResponse) {
         printf("\n***** ERROR *****\n fPIDResponse missing, Status:QUIT!\n");
         exit(1);
     }

     fListHist = new TList();
     fListHist->SetOwner(kTRUE);

     SetupEventAndTaskConfigInfo();

     fCentDistBeforCut = new TH1F("fCentDistBeforCut", "Cent w/o any Cuts; Cent (%); no.Events ", 100, 0, 100);
     fListHist->Add(fCentDistBeforCut);

     fCentDistAfterCut = new TH1F("fCentDistAfterCut", "Cent with all Cut; Cent (%); no.Events ", 100, 0, 100);
     fListHist->Add(fCentDistAfterCut);

     fHistEtaPtBeforCut = new TH2F("fHistEtaPtBeforCut", "#eta vs p_{T} (wFB, w/o  cut)", 24, -1.2, 1.2, 100, 0, 10);
     fListHist->Add(fHistEtaPtBeforCut);

     fHistEtaPhiBeforCut = new TH2F("fHistPhiEtaBeforCut", "#phi vs #eta (wFB, w/o cut)", 50, 0, 6.2835, 16, -0.8, 0.8);
     fListHist->Add(fHistEtaPhiBeforCut);

     fHistEtaPhiAfterCut = new TH2F("fHistEtaPhiAfterCut", "#phi vs #eta (wFB, w/o cut)", 50, 0, 6.2835, 16, -0.8, 0.8);
     fListHist->Add(fHistEtaPhiAfterCut);

     //Analysis Histograms:
     fHistPtDistMBwithCut = new TH1F("fHistPtDistMBwithCut", "p_{T} with Cut; p_{T} (GeV/c); entries ", 200, 0, 10);
     fListHist->Add(fHistPtDistMBwithCut);

     Int_t gMaxGlobalmult = 4000;
     Int_t gMaxTPCcorrmult = 5000;
     Int_t gMaxESDtracks = 20000;

     fHistGlobalVsV0MBefore = new TH2F("fHistGlobalVsV0MBefore", "Before;Cent(V0M);Global", 100, 0, 100, 500, 0, gMaxGlobalmult);
     fListHist->Add(fHistGlobalVsV0MBefore);
     fHistGlobalVsV0MAfter = new TH2F("fHistGlobalVsV0MAfter", " After; Cent(V0M);Global", 100, 0, 100, 500, 0, gMaxGlobalmult);
     fListHist->Add(fHistGlobalVsV0MAfter);

     fHistTPConlyVsCL1Before = new TH2F("fHistTPConlyVsCL1Before", "Before;Cent(CL1); TPC(FB128)", 100, 0, 100, 250, 0, gMaxTPCcorrmult);
     fListHist->Add(fHistTPConlyVsCL1Before);
     fHistTPConlyVsCL1After = new TH2F("fHistTPConlyVsCL1After", "After; Cent(CL1); TPC(FB128) ", 100, 0, 100, 250, 0, gMaxTPCcorrmult);
     fListHist->Add(fHistTPConlyVsCL1After);

     fHistTPConlyVsV0MBefore = new TH2F("fHistTPConlyVsV0MBefore", "Before;Cent(V0M); TPC(FB128)", 100, 0, 100, 250, 0, gMaxTPCcorrmult);
     fListHist->Add(fHistTPConlyVsV0MBefore);
     fHistTPConlyVsV0MAfter = new TH2F("fHistTPConlyVsV0MAfter", "After; Cent(V0M); TPC(FB128) ", 100, 0, 100, 250, 0, gMaxTPCcorrmult);
     fListHist->Add(fHistTPConlyVsV0MAfter);

     fHistCL0VsV0MBefore = new TH2F("fHistCL0VsV0MBefore", "Before;Cent(V0M); Cent(CL0)", 100, 0, 100, 100, 0, 100);
     fListHist->Add(fHistCL0VsV0MBefore);
     fHistCL0VsV0MAfter = new TH2F("fHistCL0VsV0MAfter", "After; Cent(V0M); Cent(CL0) ", 100, 0, 100, 100, 0, 100);
     fListHist->Add(fHistCL0VsV0MAfter);

     fTPCvsGlobalTrkBefore = new TH2F("fTPCvsGlobalTrkBefore", "Global(FB32) vs TPC(FB128)", 250, 0, 5000, 250, 0, 5000);
     fListHist->Add(fTPCvsGlobalTrkBefore);
     fTPCvsGlobalTrkAfter = new TH2F("fTPCvsGlobalTrkAfter", "Global(FB32) vs TPC(FB128)", 250, 0, 5000, 250, 0, 5000);
     fListHist->Add(fTPCvsGlobalTrkAfter);

     fHistTPCVsESDTrkBefore = new TH2F("fHistTPCVsESDTrkBefore", "Before; TPC1; ESD trk", 500, 0, gMaxTPCcorrmult, 200, 0, gMaxESDtracks);
     fListHist->Add(fHistTPCVsESDTrkBefore);
     fHistTPCVsESDTrkAfter = new TH2F("fHistTPCVsESDTrkAfter", " After;  TPC1; ESD trk", 500, 0, gMaxTPCcorrmult, 200, 0, gMaxESDtracks);
     fListHist->Add(fHistTPCVsESDTrkAfter);

     // Prottay's histograms:

     Char_t name[1000];
     Char_t title[1000];
     Double_t fCentArray[11] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};

     ////-------- Particles ------------
     fHistPtRecChrgPos = new TH2F("fHistPtRecChrgPos", "RecoMatch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtGenChrgPos = new TH2F("fHistPtGenChrgPos", "Generated;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecChrgPos);
     fListHist->Add(fHistPtGenChrgPos);
     fHistPtRecPionPos = new TH2F("fHistPtRecPionPos", "RecoMatch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtGenPionPos = new TH2F("fHistPtGenPionPos", "Generated;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecPionPos);
     fListHist->Add(fHistPtGenPionPos);

     fHistPtRecPionPosimp = new TH2F("fHistPtRecPionPosimp", "RecoMatchimp;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecPionPosimp);

     fHistPtRecKaonPos = new TH2F("fHistPtRecKaonPos", "RecoMatch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtGenKaonPos = new TH2F("fHistPtGenKaonPos", "Generated;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecKaonPos);
     fListHist->Add(fHistPtGenKaonPos);

     fHistPtRecKaonPosimp = new TH2F("fHistPtRecKaonPosimp", "RecoMatchimp;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecKaonPosimp);

     fHistPtRecProtPos = new TH2F("fHistPtRecProtPos", "RecoMatch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtGenProtPos = new TH2F("fHistPtGenProtPos", "Generated;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecProtPos);
     fListHist->Add(fHistPtGenProtPos);

     fHistPtRecProtPosimp = new TH2F("fHistPtRecProtPosimp", "RecoMatchimp;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecProtPosimp);

     ////------ Anti-particles -----------
     fHistPtRecChrgNeg = new TH2F("fHistPtRecChrgNeg", "RecoMatch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtGenChrgNeg = new TH2F("fHistPtGenChrgNeg", "Generated;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecChrgNeg);
     fListHist->Add(fHistPtGenChrgNeg);
     fHistPtRecPionNeg = new TH2F("fHistPtRecPionNeg", "RecoMatch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtGenPionNeg = new TH2F("fHistPtGenPionNeg", "Generated;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecPionNeg);
     fListHist->Add(fHistPtGenPionNeg);

     fHistPtRecPionNegimp = new TH2F("fHistPtRecPionNegimp", "RecoMatchimp;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecPionNegimp);

     fHistPtRecKaonNeg = new TH2F("fHistPtRecKaonNeg", "RecoMatch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtGenKaonNeg = new TH2F("fHistPtGenKaonNeg", "Generated;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecKaonNeg);
     fListHist->Add(fHistPtGenKaonNeg);

     fHistPtRecKaonNegimp = new TH2F("fHistPtRecKaonNegimp", "RecoMatchimp;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecKaonNegimp);

     fHistPtRecProtNeg = new TH2F("fHistPtRecProtNeg", "RecoMatch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtGenProtNeg = new TH2F("fHistPtGenProtNeg", "Generated;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecProtNeg);
     fListHist->Add(fHistPtGenProtNeg);

     fHistPtRecProtNegimp = new TH2F("fHistPtRecProtNegimp", "RecoMatchimp;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtRecProtNegimp);

     ///Change centrality binning.
     fHistnSigmaTPCPionCent1 = new TH2F("fHistnSigmaTPCPionCent1", "Reco,dEdx; p(mom), GeV/c ;n#sigma_{TPC};", 100, 0.1, 10.1, 300, -15, 15);
     fListHist->Add(fHistnSigmaTPCPionCent1);
     fHistnSigmaTPCKaonCent1 = new TH2F("fHistnSigmaTPCKaonCent1", "Reco,dEdx; p(mom), GeV/c ;n#sigma_{TPC};", 100, 0.1, 10.1, 300, -15, 15);
     fListHist->Add(fHistnSigmaTPCKaonCent1);
     fHistnSigmaTPCProtCent1 = new TH2F("fHistnSigmaTPCProtCent1", "Reco,dEdx; p(mom), GeV/c ;n#sigma_{TPC};", 100, 0.1, 10.1, 300, -15, 15);
     fListHist->Add(fHistnSigmaTPCProtCent1);
     fHistnSigmaTOFPionCent1 = new TH2F("fHistnSigmaTOFPionCent1", "Reco,dEdx; p(mom), GeV/c ;n#sigma_{TOF};", 100, 0.1, 10.1, 300, -15, 15);
     fListHist->Add(fHistnSigmaTOFPionCent1);
     fHistnSigmaTOFKaonCent1 = new TH2F("fHistnSigmaTOFKaonCent1", "Reco,dEdx; p(mom), GeV/c ;n#sigma_{TOF};", 100, 0.1, 10.1, 300, -15, 15);
     fListHist->Add(fHistnSigmaTOFKaonCent1);
     fHistnSigmaTOFProtCent1 = new TH2F("fHistnSigmaTOFProtCent1", "Reco,dEdx; p(mom), GeV/c ;n#sigma_{TOF};", 100, 0.1, 10.1, 300, -15, 15);
     fListHist->Add(fHistnSigmaTOFProtCent1);


     // contamination primary and secondary
     fHistPtContaminationPrimariesChrgPos = new TH2F("fHistPtContaminationPrimariesChrgPos", "Primaries Pos Ch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtContaminationPrimariesChrgNeg = new TH2F("fHistPtContaminationPrimariesChrgNeg", "Primaries Neg Ch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtContaminationPrimariesChrgPos);
     fListHist->Add(fHistPtContaminationPrimariesChrgNeg);
     
     fHistPtContaminationSecondariesChrgPos = new TH2F("fHistPtContaminationSecondariesChrgPos", "Secondaries Pos Ch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtContaminationSecondariesChrgNeg = new TH2F("fHistPtContaminationSecondariesChrgNeg", "Secondaries Neg Ch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtContaminationSecondariesChrgPos);
     fListHist->Add(fHistPtContaminationSecondariesChrgNeg);
     
     fHistPtContaminationSecondariesMaterialChrgPos = new TH2F("fHistPtContaminationSecondariesMaterialChrgPos", "Secondaries Material Pos Ch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtContaminationSecondariesMaterialChrgNeg = new TH2F("fHistPtContaminationSecondariesMaterialChrgNeg", "Secondaries Material Neg Ch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtContaminationSecondariesMaterialChrgPos);
     fListHist->Add(fHistPtContaminationSecondariesMaterialChrgNeg);
     
     fHistPtContaminationSecondariesWeakDecChrgPos = new TH2F("fHistPtContaminationSecondariesWeakDecChrgPos", "Secondaries Weak Dec Pos Ch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fHistPtContaminationSecondariesWeakDecChrgNeg = new TH2F("fHistPtContaminationSecondariesWeakDecChrgNeg", "Secondaries Weak Dec Neg Ch;Centrality;p_{T}(GeV/c);", 10, fCentArray, 200, 0.1, 10.1);
     fListHist->Add(fHistPtContaminationSecondariesWeakDecChrgPos);
     fListHist->Add(fHistPtContaminationSecondariesWeakDecChrgNeg);

  
     fMCVzDistBeforCut = new TH2F("fMCVzDistBeforCut", ";Centrality; MC_Trk Vz(cm);", 10, fCentArray, 80, -20, 20);
     fListHist->Add(fMCVzDistBeforCut);

     for (int i = 0; i < 10; i++) {
         ////Charge:
         sprintf(name, "fHistAch_Cent%d", i);
         sprintf(title, "Cent %2.0f-%2.0f; A_{ch}; v_{2}", i, i + 1);
         fHistAch[i] = new TProfile(name, title, 100, -1, 1, "");
         fHistAch[i]->Sumw2();
         fListHist->Add(fHistAch[i]);

     }

     fHistAChrggenVsCent = new TH2F("fHistAChrggenVsCent", "GenCent:Ach", 18, 0, 90, 1000, -1.0, 1.0);
     fListHist->Add(fHistAChrggenVsCent);
     fHistAChrgrecVsCent = new TH2F("fHistAChrgrecVsCent", "RecCent:Ach", 18, 0, 90, 1000, -1.0, 1.0);
     fListHist->Add(fHistAChrgrecVsCent);

     TString PIDlabel[26] = {"#gamma", "#pi^{0}", "#rho^{0}", "K_{L}^{0}", "#pi^{#pm}", "#rho^{#pm}", "#eta", "#omega", "K_{S}^{0}", "K^{0}", "K^{*0}", "K^{+}", "K^{*+}", "#eta'", "#Delta^{-}", "#Delta^{0}", "n", "p", "#Delta^{+}", "#Delta^{++}", "#Sigma^{-}", "#Lambda", "3124", "#Sigma^{0}", "#Sigma^{*0}", "#Sigma^{+}"};

     fHistPIDDecayMotherDaughterPhysicalPrimary = new TH2D("fHistPIDDecayMotherDaughterPhysicalPrimary", "fHistPIDDecayMotherDaughterPhysicalPrimary; Mother; Daughter", 26, 0, 26, 26, 0, 26);
     fListHist->Add(fHistPIDDecayMotherDaughterPhysicalPrimary);

     fHistPIDDecayMotherDaughterNotPhysicalPrimary = new TH2D("fHistPIDDecayMotherDaughterNotPhysicalPrimary", "fHistPIDDecayMotherDaughterNotPhysicalPrimary; Mother; Daughter", 26, 0, 26, 26, 0, 26);
     fListHist->Add(fHistPIDDecayMotherDaughterNotPhysicalPrimary);

     fHistPIDDecayDaughterNoMotherPhysicalPrimary = new TH1D("fHistPIDDecayDaughterNoMotherPhysicalPrimary", "fHistPIDDecayDaughterNoMotherPhysicalPrimary; Species; Entries", 26, 0, 26);
     fListHist->Add(fHistPIDDecayDaughterNoMotherPhysicalPrimary);

     fHistPIDDecayDaughterNoMotherNotPhysicalPrimary = new TH1D("fHistPIDDecayDaughterNoMotherNotPhysicalPrimary", "fHistPIDDecayDaughterNoMotherNotPhysicalPrimary; Species; Entries", 26, 0, 26);
     fListHist->Add(fHistPIDDecayDaughterNoMotherNotPhysicalPrimary);

     for (Int_t i = 0; i < 26; i++) {
         fHistPIDDecayMotherDaughterPhysicalPrimary->GetXaxis()->SetBinLabel(i + 1, PIDlabel[i].Data());
         fHistPIDDecayMotherDaughterPhysicalPrimary->GetYaxis()->SetBinLabel(i + 1, PIDlabel[i].Data());
         fHistPIDDecayMotherDaughterNotPhysicalPrimary->GetXaxis()->SetBinLabel(i + 1, PIDlabel[i].Data());
         fHistPIDDecayMotherDaughterNotPhysicalPrimary->GetYaxis()->SetBinLabel(i + 1, PIDlabel[i].Data());
         fHistPIDDecayDaughterNoMotherPhysicalPrimary->GetXaxis()->SetBinLabel(i + 1, PIDlabel[i].Data());
         fHistPIDDecayDaughterNoMotherNotPhysicalPrimary->GetXaxis()->SetBinLabel(i + 1, PIDlabel[i].Data());
     }

     TString GeneEventCounterLabel[8] = {"Total MC", "Pass VtxZ", "Pass #eta", "Pass hasGenerator", "Pass NoMother", "Pass PhysPrim", "Pass PileUp", "Charged"};
     fQAGeneEventCounter = new TH1D("fQAGeneEventCounter", "fQAGeneEventCounter; ; Entries", 8, 0, 8);
     fListHist->Add(fQAGeneEventCounter);
     for (Int_t i = 0; i < 8; i++) {
         fQAGeneEventCounter->GetXaxis()->SetBinLabel(i + 1, GeneEventCounterLabel[i].Data());
     }

     TString RecoEventCounterLabel[7] = {"Total", "Pass FB", "Pass MatchMC", "Pass hasGenerator", "Pass NoMother", "Pass PhysPrim", "Pass KineCut chged"};
     fQARecoEventCounter = new TH1D("fQARecoEventCounter", "fQARecoEventCounter; ; Entries", 7, 0, 7);
     fListHist->Add(fQARecoEventCounter);
     for (Int_t i = 0; i < 7; i++) {
         fQARecoEventCounter->GetXaxis()->SetBinLabel(i + 1, RecoEventCounterLabel[i].Data());
     }

     PostData(1, fListHist);
     //std::cout<<"\n UserCreateOutputObject: function Ends...\n"<<endl;
 }

 //____________________________ Call Event by Event ___________________________________
 void AliAnalysisTaskRecEff::UserExec(Option_t * ) {

     //cout<<"\n Info:UserExec() called ..!!!\n";
     //watch.Start(kTRUE);

     Float_t stepCount = 0.5;

     fHistEventCount->Fill(stepCount); //1
     stepCount++;

     fAOD = dynamic_cast < AliAODEvent * > (InputEvent());
     fESD = dynamic_cast < AliESDEvent * > (InputEvent());
     if (!(fESD || fAOD)) {
         printf("ERROR: fESD & fAOD not available\n");
         return;
     }

     fVevent = dynamic_cast < AliVEvent * > (InputEvent());
     if (!fVevent) {
         printf("ERROR: fVevent not available\n");
         return;
     }

     AliMCEvent * mcEvent = dynamic_cast < AliMCEvent * > (MCEvent());
     //AliMCEvent *mcEvent = dynamic_cast<AliMCEvent *>(InputEvent());
     if (!mcEvent) {
         Printf("ERROR: Could not retrieve MC event. Exit");
         return;
     }

     fHistEventCount->Fill(stepCount);
     stepCount++;

     ///  get stack ///
     //AliStack *mcStack = mcEvent->Stack();
     //if (!mcStack) {
     //Printf("ERROR: Could not retrieve MC Stack. Exit");
     //return;
     //}

     //fHistEventCount->Fill(stepCount); //2
     //stepCount++;

     //-------------- Vtx cuts ---------------
     const AliVVertex * pVtx = fVevent->GetPrimaryVertex();
     Double_t pVtxZ = -999;
     pVtxZ = pVtx->GetZ();

     if (pVtxZ < fMinVzCut || pVtxZ > fMaxVzCut) {
         return;
     }

     fHistEventCount->Fill(stepCount); //3
     stepCount++;

     Float_t centrality = -99.0;
     Float_t centrV0M = -99.0;
     Float_t centrCL1 = -99.0;

     fMultSelection = (AliMultSelection * ) InputEvent()->FindListObject("MultSelection"); // Must never comment this
     if (!fMultSelection) {
         printf("\n...**ERROR**...\n UserExec() AliMultSelection object not found\n Status:Quit!! \n");
         exit(1);
     }

     centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
     //centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");

     centrality = centrV0M; // This Is Always Default, changes below for other options:
     /*  
     if(sCentrEstimator=="CL0"){
       centrality = fMultSelection->GetMultiplicityPercentile("CL0");;
     }
     else if(sCentrEstimator=="CL1"){
       centrality = fMultSelection->GetMultiplicityPercentile("CL1");;
     }
     else if(sCentrEstimator=="V0C"){
       centrality = fMultSelection->GetMultiplicityPercentile("V0C");
     }
     else if(sCentrEstimator=="V0A"){
       centrality = fMultSelection->GetMultiplicityPercentile("V0A");
     }
     else if(sCentrEstimator=="TRK"){
       centrality = fMultSelection->GetMultiplicityPercentile("TRK");
     }
     else{
       centrality = centrV0M;
     }
     */

     fCentDistBeforCut->Fill(centrality);

     Bool_t kPileupEvent = kFALSE;

     //kPileupEvent = CheckEventIsPileUp2018(fAOD);
     //if (!bSkipPileUpCut)
     //if(kPileupEvent)  return;

     if (centrality < fCentralityMin || centrality > fCentralityMax) {
         return;
     }

     ///----> Get 0-10 index of centrality:
     Int_t cent10bin = -99, iCent = -99;

     if (centrality < 5.0) {
         cent10bin = 0;
     } else if (centrality >= 5.0 && centrality < 10) {
         cent10bin = 1;
     } else if (centrality >= 10.0) {
         cent10bin = abs(centrality / 10.0) + 1;
     }

     iCent = cent10bin;

     fHistEventCount->Fill(stepCount); //4
     stepCount++;

     Int_t ntracks = fAOD->GetNumberOfTracks();
     if (ntracks < 4) return; // Minimum 4 tracks per event. 

     fHistEventCount->Fill(stepCount); //5
     stepCount++;

     //////----> Get Magnetic field and RunNo.---------
     // Float_t fMagField = fAOD->GetMagneticField();
     // const Int_t QAindex = (fMagField > 0) ? 1 : 0;
     // Int_t runNumber = fAOD->GetRunNumber();
     //------------------------------------------------

     Float_t fMultTPCFull = 0; // TPC mult estimate
     Float_t fMultGlobal = 0; // global track multiplicity
     Float_t fMultwRawFB = 0; // Uncorrected Multiplicity
     Float_t fMultCorrFB = 0; // Corrected Multiplicity

     Int_t gMultEtaNeg = 0;
     Int_t gMultEtaPos = 0;
     Int_t gMultEtaAll = 0;

     Float_t trkPt = 0.0, trkPhi = 0.0, trkEta = 0.0, trkRap = 0.0;
     Float_t trkP = 0.0, trkChi2 = 0.0, trkdEdx = 0.0, trkDCAxy = 0.0, trkDCAz = 0.0;
     Int_t trkChrg = 0, trkTpcNC = 0, trkTPCsharedNC = 0;
     Float_t NPosgen = 0., NNeggen = 0.;
     Float_t NPosrec = 0., NNegrec = 0.;
     Float_t Achgen = 0.0, Achrec = 0.0;

     //Prottay's variable:
     Double_t Q2xEtaPos = 0., Q2yEtaPos = 0.;
     Double_t Q2xEtaNeg = 0., Q2yEtaNeg = 0.;

     ////PID related variables
     Double_t nSigTOFpion, nSigTPCpion;
     Double_t nSigTOFkaon, nSigTPCkaon;
     Double_t nSigTOFprot, nSigTPCprot;

     Double_t nSigRMSTPCTOFpion;
     Double_t nSigRMSTPCTOFkaon;
     Double_t nSigRMSTPCTOFprot;

     Bool_t bTOFmatch = kFALSE;
     Bool_t isItPion = kFALSE, isItKaon = kFALSE, isItProt = kFALSE;
     Float_t cent = centrality;

     fHistEventCount->Fill(stepCount);
     stepCount++;

     ///////// Rihan: MC Code. Credit to Davide Caffari.

     AliGenEventHeader * eventHeader = 0;
     //For AOD only
     AliAODMCHeader * header = (AliAODMCHeader * ) fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
     if (!header) {
         AliFatal("fInjectedSignals set but no MC header found");
         return;
     }

     eventHeader = header->GetCocktailHeader(0);

     if (!eventHeader) {
         // We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
         // (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
         AliError("First event header not found. Skipping this event.");
         return;
     }

     fHistEventCount->Fill(stepCount);
     stepCount++;
     //skipParticlesAbove = eventHeader->NProduced();
     //AliInfo(Form("Injected signals in this event (%d headers). Keeping particles/tracks of %s. Will skip particles/tracks above %d.", headers, eventHeader->ClassName(), skipParticlesAbove)); 

     TClonesArray * arrayMC = (TClonesArray * ) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
     if (!arrayMC) {
         ::Fatal("AliAnalysisTaskSimle::UserExec", " dMC particles branch not found.");
     }

     Int_t PIDarray[26] = {22, 111, 113, 130, 211, 213, 221, 223, 310, 311, 313, 321, 323, 331, 1114, 2114, 2112, 2212, 2214, 2224, 3112, 3122, 3124, 3212, 3214, 3222};

     //++++++++++++++++++ MC tracks+++++++++++++++++++++//
     Int_t nMCLabelCounter = 0;
     Int_t nMCParticles = mcEvent->GetNumberOfTracks();
     Int_t pdgIndex = 0;
     Int_t LabelOfMother = 0;
     Int_t pdgcode = 0;
     Double_t vzMC = 0;
     TArrayI labelMCArray(nMCParticles);

     for (Int_t iTracks = 0; iTracks < nMCParticles; iTracks++) {
         fQAGeneEventCounter->Fill(0.5, 1);
         AliAODMCParticle * mcTrack = (AliAODMCParticle * ) mcEvent->GetTrack(iTracks);

         if (!mcTrack) {
             AliError(Form("ERROR: Could not receive track %d (mc loop)", iTracks));
             continue;
         }

         //exclude particles generated out of the acceptance
         vzMC = mcTrack->Zv();
         fMCVzDistBeforCut->Fill(cent, vzMC);

         //if(TMath::Abs(vzMC) > 10.0) continue;
         if (TMath::Abs(vzMC) > fMaxVzCut) continue;
         fQAGeneEventCounter->Fill(1.5, 1);
         trkEta = mcTrack->Eta();
         //acceptance
         if (fabs(trkEta) > 0.8) continue;
         fQAGeneEventCounter->Fill(2.5, 1);
         Double_t y = mcTrack->Y();
         //if(TMath::Abs(y) < 0.5) continue;  //  why do I need this cut?
         //if(TMath::Abs(y) > 0.8) continue;  //  why do I need this cut?

         trkPt = mcTrack->Pt();
         //if(mcTrack->Pt() > 5.0 || mcTrack->Pt() < 0.1)   continue;

         TString generatorName;
         Bool_t hasGenerator = mcEvent->GetCocktailGenerator(iTracks, generatorName);
         //cout<<"NAME OF GENERATOR---------------------------------"<<generatorName<<endl;
         //if((!hasGenerator) || (!generatorName.Contains("AMPT")))  continue;
         if ((!hasGenerator)) continue;
         fQAGeneEventCounter->Fill(3.5, 1);
         //if (nMCPassNoMother<200) 
         //cout<<"===> before mother cut: mcTrack->GetPdgCode() = "<<mcTrack->GetPdgCode()<<", mcTrack->GetMother() = "<<mcTrack->GetMother()<<", mcTrack->IsPhysicalPrimary() = "<<mcTrack->IsPhysicalPrimary()<<endl;

         //if (nMCPassNoMother < 200 && mcTrack->GetMother() != -1) {
         //AliAODMCParticle *mcTrackFFF = (AliAODMCParticle*) mcEvent->GetTrack(mcTrack->GetMother()); 
         //cout<<"===> mother PID = "<<mcTrackFFF->GetPdgCode()<<endl;

         //}
         AliAODMCParticle * mcTrackMother = (AliAODMCParticle * ) mcEvent->GetTrack(mcTrack->GetMother());

         Bool_t findMother = kFALSE;
         Bool_t findDaughter = kFALSE;
         Int_t MotherIndex = -1;
         Int_t DaughterIndex = -1;

         for (Int_t i = 0; i < 26; i++) {
             if (!(mcTrack->GetMother() == -1)) {
                 if (mcTrackMother->GetPdgCode() == PIDarray[i]) {
                     findMother = kTRUE;
                     MotherIndex = i;
                 }
             }
             if (mcTrack->GetPdgCode() == PIDarray[i]) {
                 findDaughter = kTRUE;
                 DaughterIndex = i;
             }

         }

         if (mcTrack->IsPhysicalPrimary()) {
             if (findDaughter) {
                 if (mcTrack->GetMother() == -1) {
                     fHistPIDDecayDaughterNoMotherPhysicalPrimary->Fill(DaughterIndex + 0.5, 1);
                 } else if (findMother) {
                     fHistPIDDecayMotherDaughterPhysicalPrimary->Fill(MotherIndex + 0.5, DaughterIndex + 0.5, 1);
                 }
             }
         } else {
             if (findDaughter) {
                 if (mcTrack->GetMother() == -1) {
                     fHistPIDDecayDaughterNoMotherNotPhysicalPrimary->Fill(DaughterIndex + 0.5, 1);
                 } else if (findMother) {
                     fHistPIDDecayMotherDaughterNotPhysicalPrimary->Fill(MotherIndex + 0.5, DaughterIndex + 0.5, 1);
                 }
             }
         }

         //LabelOfMother = mcTrack->GetMother();
         //if(LabelOfMother>0)               continue; //only using primary track
         fQAGeneEventCounter->Fill(4.5, 1);
         if (!mcTrack->IsPhysicalPrimary()) continue; //only using primary track
         fQAGeneEventCounter->Fill(5.5, 1);
         if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTracks, header, arrayMC))
             continue;
         fQAGeneEventCounter->Fill(6.5, 1);

         //--------- Now load other mcTrack parameters ------------
         trkChrg = mcTrack->Charge();
         if (!trkChrg) continue;
         fQAGeneEventCounter->Fill(7.5, 1);
         trkPhi = mcTrack->Phi();

         pdgcode = TMath::Abs(mcTrack->GetPdgCode());
         //if(pdgcode==0) continue;

         if (trkChrg > 0) {
             trkChrg = 1;
             if (trkPt > 0.2 && trkPt < 10.0)
                 NPosgen += 1;

         } else if (trkChrg < 0) {
             trkChrg = -1;
             if (trkPt > 0.2 && trkPt < 10.0)
                 NNeggen += 1;
         }

         if (trkChrg > 0) {
             fHistPtGenChrgPos->Fill(cent, trkPt);

             if (pdgcode == 211) {
                 fHistPtGenPionPos->Fill(cent, trkPt);
             } else if (pdgcode == 321) {
                 fHistPtGenKaonPos->Fill(cent, trkPt);
             } else if (pdgcode == 2212) {
                 fHistPtGenProtPos->Fill(cent, trkPt);
             }
         } else if (trkChrg < 0) {
             fHistPtGenChrgNeg->Fill(cent, trkPt);
             ////negetive MC tracks  
             if (pdgcode == 211) {
                 fHistPtGenPionNeg->Fill(cent, trkPt);
             } else if (pdgcode == 321) {
                 fHistPtGenKaonNeg->Fill(cent, trkPt);
             } else if (pdgcode == 2212) {
                 fHistPtGenProtNeg->Fill(cent, trkPt);
             }
         }

         fHistEtaPhiBeforCut->Fill(trkPhi, trkEta);

         labelMCArray.AddAt(iTracks, nMCLabelCounter);

         if (nMCLabelCounter > 50000) break;

         nMCLabelCounter++;

     } //mc track loop

     if ((NPosgen + NNeggen) != 0)
         Achgen = (NPosgen - NNeggen) / (NPosgen + NNeggen);

     //cout<<"value of achgen------------------"<<NPosgen<<" "<<NNeggen<<" "<<Achgen<<endl;

     fHistEventCount->Fill(stepCount);
     stepCount++;

     Int_t nMatchedRecoTrk = 0;

     //------------- reco tracks--------------------------
     Int_t nGoodTracks = fAOD->GetNumberOfTracks();
     TArrayI labelArray(nGoodTracks);
     Int_t labelCounter = 0;

     for (Int_t itrack = 0; itrack < nGoodTracks; itrack++) {
         fQARecoEventCounter->Fill(0.5, 1);
         AliAODTrack * AODtrack = dynamic_cast < AliAODTrack * > (fVevent->GetTrack(itrack));
         if (!AODtrack) continue;

         if (!AODtrack->TestFilterBit(fFilterBit)) continue;
         fQARecoEventCounter->Fill(1.5, 1);
         Int_t label = TMath::Abs(AODtrack->GetLabel());
         //if(itrack%10==0) cout<<"Trk "<<itrack<<" Gene = "<<generatorName.Data()<<"\t dPt1 = "<<dPt1<<"\t dEta1 = "<<dEta1<<endl;

         //-------- Get back matched MC particle ----------------
         Int_t mcGoods = nMCLabelCounter;

         for (Int_t k = 0; k < mcGoods; k++) {

             Int_t mcLabel = labelMCArray.At(k);
             //if (mcLabel != label) continue;

             if (mcLabel == label) {
                 AliAODMCParticle * mcTracMatchedWithReco = (AliAODMCParticle * ) mcEvent->GetTrack(label);

                 if (!mcTracMatchedWithReco) {
                     std::cout << " ERROR: Could not receive MC generated track for #" << label << std::endl;
                     break;
                 }
                 fQARecoEventCounter->Fill(2.5, 1);
                 pdgcode = TMath::Abs(mcTracMatchedWithReco->GetPdgCode());

                 TString generatorName;
                 Bool_t hasGenerator = mcEvent->GetCocktailGenerator(label, generatorName);
                 //if((!hasGenerator) || (!generatorName.Contains("Hijing")))
                 if ((!hasGenerator))
                     break;
                 fQARecoEventCounter->Fill(3.5, 1);

                 //LabelOfMother = mcTracMatchedWithReco->GetMother();
                 //if(LabelOfMother>0) break;
                 fQARecoEventCounter->Fill(4.5, 1);
                 //if(itrack%10==0) cout<<"MC track "<<k<<" pdgcode = "<<pdgcode<<"\t dPt1 = "<<dPt1<<"\t dEta1 = "<<dEta1<<endl;

                 if (!mcTracMatchedWithReco->IsPhysicalPrimary()) break;
                 fQARecoEventCounter->Fill(5.5, 1);
                 //Double_t y = mcTracMatchedWithReco->Y();
                 //if(TMath::Abs(y) > 1.0) break;  // dont need in AOD tracks as I apply eta cut.

                 //	if(fAnalysisUtil->IsParticleFromOutOfBunchPileupCollision(label,mcEvent)) ///remove PileUp tracks      
                 //continue; 

                 //--------- Now load AOD track parameters ------------

                 trkP = AODtrack->P();
                 trkPt = AODtrack->Pt();
                 trkPhi = AODtrack->Phi();
                 trkEta = AODtrack->Eta();
                 trkRap = TMath::Abs(AODtrack->Y());
                 trkChrg = AODtrack->Charge();
                 trkChi2 = AODtrack->Chi2perNDF();
                 trkTpcNC = AODtrack->GetTPCNcls();
                 trkTPCsharedNC = AODtrack->GetTPCnclsS();
                 trkDCAxy = AODtrack->DCA();
                 trkDCAz = AODtrack->ZAtDCA();

                 /*

	Double_t dTrackXYZ[3] = {0};
	Double_t dVertexXYZ[3] = {0.};
	Double_t dDCAXYZ[3] = {0.};
	//const AliAODVertex* vertex = fAOD->GetPrimaryVertex();
	//if(!vertex) return kFALSE; //event does not have a PV
        
	AODtrack->GetXYZ(dTrackXYZ);
	pVtx->GetXYZ(dVertexXYZ);
        
	for(Short_t i(0); i < 3; i++)
	  dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i];
    
	trkDCAxy=TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]);
        trkDCAz=dDCAXYZ[2];
	*/
                 /// This Next function is called After FilterBit is validated!! Otherwise code breaks.!
                 trkdEdx = AODtrack->GetDetPid()->GetTPCsignal();

                 //Apply track cuts here:
                 if ((trkPt <= fMaxPtCut) && (trkPt >= fMinPtCut) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkTPCsharedNC <= fTPCclustMin * 0.8) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= 4.0) && TMath::Abs(trkChrg)) {
                     //if((trkPt <= fMaxPtCut) && (trkPt >= fMinPtCut) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkRap <= 0.8) && (trkRap >= 0.5) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= 5.0) && TMath::Abs(trkChrg)) {
                     ///if((AODtrack->DCA() <= fDCAxyMax) && (AODtrack->ZAtDCA() <= fDCAzMax)) /// DCA cut not needed if I am using FB.
                     fQARecoEventCounter->Fill(6.5, 1);
                     isItPion = kFALSE;
                     isItKaon = kFALSE;
                     isItProt = kFALSE;

                     ///=========> Get TPC/TOF nSigma for PID
                     nSigTPCpion = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kPion);
                     nSigTPCkaon = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kKaon);
                     nSigTPCprot = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kProton);

                     nSigTOFpion = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kPion);
                     nSigTOFkaon = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kKaon);
                     nSigTOFprot = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kProton);

                     nSigRMSTPCTOFpion = TMath::Sqrt(nSigTPCpion * nSigTPCpion + nSigTOFpion * nSigTOFpion);
                     nSigRMSTPCTOFkaon = TMath::Sqrt(nSigTPCkaon * nSigTPCkaon + nSigTOFkaon * nSigTOFkaon);
                     nSigRMSTPCTOFprot = TMath::Sqrt(nSigTPCprot * nSigTPCprot + nSigTOFprot * nSigTOFprot);

                     //----- Pion
                     if (trkPt <= 0.5 && TMath::Abs(nSigTPCpion) <= 3.0) {
                         isItPion = kTRUE;
                     }
                     // else if(trkPt>0.6 && TMath::Abs(nSigTPCpion)<=3.0 && TMath::Abs(nSigTOFpion)<=3.0 ){
                     else if (trkPt > 0.5 && TMath::Abs(nSigRMSTPCTOFpion) <= 3.0) {
                         isItPion = kTRUE;
                     }
                     //----- Kaon
                     if (trkPt <= 0.45 && TMath::Abs(nSigTPCkaon) <= 3.0) {
                         isItKaon = kTRUE;
                     } else if (trkPt > 0.45 && TMath::Abs(nSigRMSTPCTOFkaon) <= 2.5) {
                         //else if(trkPt>0.45 && TMath::Abs(nSigTPCkaon)<=3.0 && TMath::Abs(nSigTOFkaon)<=2.5 ){
                         isItKaon = kTRUE;
                     }
                     //----- Proton 
                     if (trkPt <= 0.6 && TMath::Abs(nSigTPCprot) <= 3.0) {
                         isItProt = kTRUE;
                         if (trkChrg > 0 && trkPt < 0.4) isItProt = kFALSE; //Proton below 0.4 GeV has beam Pipe Contamination
                     }
                     //else if(trkPt>0.8 && TMath::Abs(nSigTPCprot)<=3.0 && TMath::Abs(nSigTOFprot)<=3.0 ){
                     else if (trkPt > 0.6 && TMath::Abs(nSigRMSTPCTOFprot) <= 3.0) {
                         isItProt = kTRUE;
                     }

                     ///Fill nSigma for purity estimation
                     /*
                     fHistnSigmaTPCPionCent1->Fill(trkP,nSigTPCpion);
                     fHistnSigmaTOFPionCent1->Fill(trkP,nSigTOFpion);	  

                     fHistnSigmaTPCKaonCent1->Fill(trkP,nSigTPCkaon);
                     fHistnSigmaTOFKaonCent1->Fill(trkP,nSigTOFkaon);

                     fHistnSigmaTPCProtCent1->Fill(trkP,nSigTPCprot);
                     fHistnSigmaTOFProtCent1->Fill(trkP,nSigTOFprot);
                     */

                     //Get MC Generated pT, as we did not take care of Energyloss:
                     trkPt = mcTracMatchedWithReco->Pt();
                     trkPhi = mcTracMatchedWithReco->Phi();

                     //Reco tracks:
                     if (trkChrg > 0) {
                         if (trkPt > 0.2 && trkPt < 10.0)
                             NPosrec += 1;
                         fHistPtRecChrgPos->Fill(cent, trkPt);
                         if (isItPion && pdgcode == 211) {
                             fHistPtRecPionPos->Fill(cent, trkPt);
                         }
                         if (isItPion)
                             fHistPtRecPionPosimp->Fill(cent, trkPt);
                         if (isItKaon && pdgcode == 321) {
                             fHistPtRecKaonPos->Fill(cent, trkPt);
                         }
                         if (isItKaon)
                             fHistPtRecKaonPosimp->Fill(cent, trkPt);
                         if (isItProt && pdgcode == 2212) {
                             fHistPtRecProtPos->Fill(cent, trkPt);
                         }
                         if (isItProt)
                             fHistPtRecProtPosimp->Fill(cent, trkPt);
                     } ///// Pos charge      
                     else if (trkChrg < 0) {
                         if (trkPt > 0.2 && trkPt < 10.0)
                             NNegrec += 1;
                         fHistPtRecChrgNeg->Fill(cent, trkPt);
                         if (isItPion && pdgcode == 211) {
                             fHistPtRecPionNeg->Fill(cent, trkPt);
                         }
                         if (isItPion)
                             fHistPtRecPionNegimp->Fill(cent, trkPt);
                         if (isItKaon && pdgcode == 321) {
                             fHistPtRecKaonNeg->Fill(cent, trkPt);
                         }
                         if (isItKaon)
                             fHistPtRecKaonNegimp->Fill(cent, trkPt);
                         if (isItProt && pdgcode == 2212) {
                             fHistPtRecProtNeg->Fill(cent, trkPt);
                         }
                         if (isItProt)
                             fHistPtRecProtNegimp->Fill(cent, trkPt);

                     } ///// Neg charge 

                     fHistEtaPhiAfterCut->Fill(trkPhi, trkEta);

                     nMatchedRecoTrk++;
                 }
             } //(if AOD track == MC track)

         } //MC track loop

     } //===================== AOD track loop ends ============================

     labelMCArray.Reset();
	 
	 
	 for (Int_t itrack = 0; itrack < nGoodTracks; itrack++) {
         AliAODTrack * AODtrack = dynamic_cast < AliAODTrack * > (fVevent->GetTrack(itrack));
         if (!AODtrack) continue;

         if (!AODtrack->TestFilterBit(fFilterBit)) continue;
         Int_t label = TMath::Abs(AODtrack->GetLabel());
         if (label > nMCParticles) continue;
         
         AliAODMCParticle * AODmcTrack = (AliAODMCParticle * ) mcEvent -> GetTrack(label);
         //Short_t gAODmcCharge = AODmcTrack -> Charge(); ////
		 
         trkP = AODtrack->P();
		 trkPt = AODtrack->Pt();
		 trkPhi = AODtrack->Phi();
		 trkEta = AODtrack->Eta();
		 trkRap = TMath::Abs(AODtrack->Y());
		 trkChrg = AODtrack->Charge();
		 trkChi2 = AODtrack->Chi2perNDF();
		 trkTpcNC = AODtrack->GetTPCNcls();
		 trkTPCsharedNC = AODtrack->GetTPCnclsS();
		 trkDCAxy = AODtrack->DCA();
		 trkDCAz = AODtrack->ZAtDCA();
         trkdEdx = AODtrack->GetDetPid()->GetTPCsignal();
         
		 //Apply track cuts here:
		 if ((trkPt <= fMaxPtCut) && (trkPt >= fMinPtCut) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkTPCsharedNC <= fTPCclustMin * 0.8) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= 4.0) && TMath::Abs(trkChrg)) {
			 if (AODmcTrack -> IsPhysicalPrimary()) { 
			     if (trkChrg > 0) {
					 fHistPtContaminationPrimariesChrgPos->Fill(cent, trkPt);
			     }
			     if (trkChrg < 0) {
					 fHistPtContaminationPrimariesChrgNeg->Fill(cent, trkPt);
			     }
			 } else {
				 Bool_t isFromMaterial = kFALSE;
				 if (AODmcTrack -> IsSecondaryFromMaterial()) isFromMaterial = kTRUE;
				 
				 if (trkChrg > 0) {
					 fHistPtContaminationSecondariesChrgPos->Fill(cent, trkPt);
					 if (isFromMaterial) {
						 fHistPtContaminationSecondariesMaterialChrgPos->Fill(cent, trkPt);
					 } else {
						 fHistPtContaminationSecondariesWeakDecChrgPos->Fill(cent, trkPt);
					 }
			     }
			     if (trkChrg < 0) {
					 fHistPtContaminationSecondariesChrgNeg->Fill(cent, trkPt);
					 if (isFromMaterial) {
						 fHistPtContaminationSecondariesMaterialChrgNeg->Fill(cent, trkPt);
					 } else {
						 fHistPtContaminationSecondariesWeakDecChrgNeg->Fill(cent, trkPt);
					 }
			     }
			 }
		 }
	 }
	 
     if ((NPosrec + NNegrec) != 0)
         Achrec = (NPosrec - NNegrec) / (NPosrec + NNegrec);

     //cout<<"value of achrec------------------"<<NPosrec<<" "<<NNegrec<<" "<<Achrec<<endl;

     fHistAch[iCent]->Fill(Achrec, Achgen, 1);

     fHistAChrggenVsCent->Fill(cent, Achgen);
     fHistAChrgrecVsCent->Fill(cent, Achrec);

     //Last lines in Event loop
     fCentDistAfterCut->Fill(centrality);
     fHistEventCount->Fill(14.5); //final Event count.
     PostData(1, fListHist);
 } //---------------- UserExec ----------------------

 void AliAnalysisTaskRecEff::SetupEventAndTaskConfigInfo() {
     fHistEventCount = new TH1F("fHistEventCount", "Event counts", 15, 0, 15);
     fHistEventCount->GetXaxis()->SetBinLabel(1, "Called UserExec()");
     fHistEventCount->GetXaxis()->SetBinLabel(2, "Called Exec()");
     fHistEventCount->GetXaxis()->SetBinLabel(3, "AOD Exist");
     fHistEventCount->GetXaxis()->SetBinLabel(4, "Vz < 10");
     fHistEventCount->GetXaxis()->SetBinLabel(5, Form("%2.0f<Cent<%2.0f", fCentralityMin, fCentralityMax));
     fHistEventCount->GetXaxis()->SetBinLabel(6, "noAODtrack > 2 ");
     fHistEventCount->GetXaxis()->SetBinLabel(7, "TPC vs Global");
     fHistEventCount->GetXaxis()->SetBinLabel(8, "TPC128 vs ESD");
     fHistEventCount->GetXaxis()->SetBinLabel(9, "Cent vs TPC");
     fHistEventCount->GetXaxis()->SetBinLabel(10, "mult eta+/-> 2");
     fHistEventCount->GetXaxis()->SetBinLabel(11, "centCL1 < 90");
     fHistEventCount->GetXaxis()->SetBinLabel(15, "Survived Events");
     fListHist->Add(fHistEventCount);
     //fHistEventCount->Fill(1);
 } //----------SetupEventAndTaskConfigInfo-----------

 Int_t AliAnalysisTaskRecEff::GetCentralityScaled0to10(Double_t fCent) {

     Int_t cIndex = 0;

     if (fCent < 5.0) {
         cIndex = 0;
     } else if (fCent >= 5.0 && fCent < 10) {
         cIndex = 1;
     } else if (fCent >= 10.0) {
         cIndex = abs(fCent / 10.0) + 1;
     }
     return cIndex;

 } //------------GetCentralityScaled0to10------------

 Bool_t AliAnalysisTaskRecEff::CheckEventIsPileUp2018(AliAODEvent * faod) {

     /*
  TF1 *fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);
    

  Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  TF1 *fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);
    

  Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  TF1 *fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  TF1 *fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);
    

  Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  TF1 *fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  fMultCutPU->SetParameters(parFB32);




  Bool_t BisPileup=kFALSE;


  Double_t centrV0M=-99.0;
  Double_t centrCL1=-99.0;
  Double_t centrCL0=-99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
 
  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
  }


  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");



  Int_t nITSClsLy0 = faod->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = faod->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;
    

  AliAODTracklets* aodTrkl = (AliAODTracklets*)faod->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();


  const Int_t nTracks = faod->GetNumberOfTracks();
    

  Int_t multTrk = 0;
    

  for (Int_t it = 0; it < nTracks; it++) {
        

    AliAODTrack* aodTrk = (AliAODTrack*)faod->GetTrack(it);
        

    if (!aodTrk){
      delete aodTrk;
      continue;
    }


        

    if (aodTrk->TestFilterBit(32)){
      if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
	multTrk++;
            
    }     

  }




    

  AliAODVZERO* aodV0 = faod->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;

  Int_t tpcClsTot = faod->GetNumberOfTPCClusters();
  Float_t nclsDif = Float_t(tpcClsTot) - (60932.9 + 69.2897*multV0Tot - 0.000217837*multV0Tot*multV0Tot);


  if (centrCL0 < fCenCutLowPU->Eval(centrV0M))
    BisPileup=kTRUE;

  if (centrCL0 > fCenCutHighPU->Eval(centrV0M))
    BisPileup=kTRUE;        


  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls))
    BisPileup=kTRUE;
        

  if (multV0On < fV0CutPU->Eval(multV0Tot))
    BisPileup=kTRUE;


  if (Float_t(multTrk) < fMultCutPU->Eval(centrV0M))
    BisPileup=kTRUE;

        
  if (((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0)
    BisPileup=kTRUE;


  if (faod->IsIncompleteDAQ())
    BisPileup=kTRUE;
        
  //if (nclsDif > 200000)//can be increased to 200000
  // BisPileup=kTRUE;



  //fHistTPConlyVsV0MBefore->Fill(centrV0M,multTrk);
  //fHistCL0VsV0MBefore->Fill(centrV0M,centrCL0);


  if (!BisPileup)
    {
      //fHistTPConlyVsV0MAfter->Fill(centrV0M,multTrk);
      //fHistCL0VsV0MAfter->Fill(centrV0M,centrCL0);
    }

  */

     TF1 * fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);
     Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
     TF1 * fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
     fV0CutPU->SetParameters(parV0);
     Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
     TF1 * fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
     fCenCutLowPU->SetParameters(parV0CL0);
     TF1 * fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
     fCenCutHighPU->SetParameters(parV0CL0);
     Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
     TF1 * fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
     fMultCutPU->SetParameters(parFB32);

     Bool_t BisPileup = kFALSE;

     Double_t centrV0M = -99.0;
     Double_t centrCL1 = -99.0;
     Double_t centrCL0 = -99.0;

     fMultSelection = (AliMultSelection * ) InputEvent()->FindListObject("MultSelection");

     if (!fMultSelection) {
         printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
         exit(1);
     }

     centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
     centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
     centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

     Int_t nITSClsLy0 = faod->GetNumberOfITSClusters(0);
     Int_t nITSClsLy1 = faod->GetNumberOfITSClusters(1);
     Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

     AliAODTracklets * aodTrkl = (AliAODTracklets * ) faod->GetTracklets();
     Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

     const Int_t nTracks = faod->GetNumberOfTracks();

     Int_t multTrk = 0;

     for (Int_t it = 0; it < nTracks; it++) {

         AliAODTrack * aodTrk = (AliAODTrack * ) faod->GetTrack(it);

         if (!aodTrk) {
             delete aodTrk;
             continue;
         }

         if (aodTrk->TestFilterBit(32)) {
             if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
                 multTrk++;

         }

     }

     AliAODVZERO * aodV0 = faod->GetVZEROData();
     Float_t multV0a = aodV0->GetMTotV0A();
     Float_t multV0c = aodV0->GetMTotV0C();
     Float_t multV0Tot = multV0a + multV0c;
     UShort_t multV0aOn = aodV0->GetTriggerChargeA();
     UShort_t multV0cOn = aodV0->GetTriggerChargeC();
     UShort_t multV0On = multV0aOn + multV0cOn;

     Int_t tpcClsTot = faod->GetNumberOfTPCClusters();
     Float_t nclsDif = Float_t(tpcClsTot) - (60932.9 + 69.2897 * multV0Tot - 0.000217837 * multV0Tot * multV0Tot);

     if (centrCL0 < fCenCutLowPU->Eval(centrV0M)) {
         //cout<<"*****************hi i am in 1st**************************:"<<endl;
         BisPileup = kTRUE;
     }

     if (centrCL0 > fCenCutHighPU->Eval(centrV0M)) {
         //cout<<"*****************hi i am in 2nd**************************:"<<endl;
         BisPileup = kTRUE;
     }

     if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) {
         //cout<<"*****************hi i am in 3rd**************************:"<<endl;
         BisPileup = kTRUE;
     }

     if (multV0On < fV0CutPU->Eval(multV0Tot)) {
         //cout<<"*****************hi i am in 4th**************************:"<<endl;
         BisPileup = kTRUE;
     }

     if (Float_t(multTrk) < fMultCutPU->Eval(centrV0M)) {
         //cout<<"*****************hi i am in 5th**************************:"<<endl;
         BisPileup = kTRUE;
     }

     if (((AliAODHeader * ) faod->GetHeader())->GetRefMultiplicityComb08() < 0)
         BisPileup = kTRUE;

     if (faod->IsIncompleteDAQ())
         BisPileup = kTRUE;

     //if (nclsDif > 200000)//can be increased to 200000
     // BisPileup=kTRUE;

     //fHistTPConlyVsV0MBefore->Fill(centrV0M,multTrk);
     //fHistCL0VsV0MBefore->Fill(centrV0M,centrCL0);

     if (!BisPileup) {
         //fHistTPConlyVsV0MAfter->Fill(centrV0M,multTrk);
         //fHistCL0VsV0MAfter->Fill(centrV0M,centrCL0);
     }

     return BisPileup;

 }
