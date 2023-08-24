
/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                        Basic Analysis Task                            //
//                      for PID        Analysis                          //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TList.h>

#include <AliVEvent.h>
#include <AliInputEventHandler.h>
#include <AliAnalysisManager.h>

#include "AliAnalysisTaskSE.h"
#include "AliStack.h"
#include "AliMultSelection.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h" 
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

#include <AliLog.h>
#include <AliPID.h>
#include <AliPIDCombined.h>
#include <AliPIDResponse.h>
#include "AliAODTrack.h"
#include "AliAnalysisTaskPIDPerformCombIDPtDep.h"


using std::cout;
using std::endl;


const char *AliAnalysisTaskPIDPerformCombIDPtDep::fgkBinMomDesc[AliAnalysisTaskPIDPerformCombIDPtDep::kPtBins] = {
" 0 <= p < 0.5 GeV/c",
  " 0.5 <= p < 0.7 GeV/c",
  " 0.7 <= p < 1.0 GeV/c",
  " 1.0 <= p < 1.5 GeV/c",
  " 1.5 <= p < 2.0 GeV/c",
  " 2.0 <= p < 2.5 GeV/c",
  " 2.5 <= p < 3.0 GeV/c",
  " 3.0 <= p < 4.0 GeV/c",
  " p >= 4.0 GeV/c"
  };

ClassImp(AliAnalysisTaskPIDPerformCombIDPtDep)

//_________________________________________________________________________________
AliAnalysisTaskPIDPerformCombIDPtDep::AliAnalysisTaskPIDPerformCombIDPtDep() :
AliAnalysisTaskSE(),
  fHistList(0),
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fAOD(0x0),
  fUseCentrality(kFALSE),
  fCentralityEstimator("V0M"), 
  fCentralityPercentileMin(0.0), 
  fCentralityPercentileMax(5.0),
  fHistCentrality(0),
  hTrue(),
  hTrueInAccTPC(),
  hTrueInAccTOF(),
  hTrueInAccTPCTOF(),
  hTrueInAccTPCTOFBayes(),
  hIdTPConly2s(),
  hIdTPConly3s(),
  hIDnSigmaComb1(),
  hIDnSigmaComb2(),
  hIDnSigmaComb3(),
  hEffPlotsTPConly2s(),
  hEffPlotsTPConly3s(),
  hIdTOFonly2s(),
  hIdTOFonly3s(),
  hIDBayes1(),
  hIDBayes2(),
  hIDBayes3(),
  hEffPlotsTOFonly2s(),
  hEffPlotsTOFonly3s(),
  hEffPlotsnSigmaComb1(),
  hEffPlotsnSigmaComb2(),
  hEffPlotsnSigmaComb3(),
  hEffPlotsBayes1(),
  hEffPlotsBayes2(),
  hEffPlotsBayes3(),
  fPriorsUsed(),
  fpartOfInterest(),
  fFB(),
  fRejectCheckGenName(kFALSE),
  fGenToBeKept("Hijing"),
  fPIDMomCut(0.6),
  fbayesth1(0.75),
  fbayesth2(0.80),
  fbayesth3(0.9),
  fTPCchi2Cut(-1),
  fNClustersTPCCut(-1),
  fMinTPCCrossedRows(-1),
  fMinTPCRowsOverFindableCls(-1)
{
  
  
}

//_________________________________________________________________________________
AliAnalysisTaskPIDPerformCombIDPtDep::AliAnalysisTaskPIDPerformCombIDPtDep(const char *name) :
  AliAnalysisTaskSE(name),
  fHistList(0),
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fAOD(0x0),
  fUseCentrality(kFALSE),
  fCentralityEstimator("V0M"), 
  fCentralityPercentileMin(0.0), 
  fCentralityPercentileMax(5.0),
  fHistCentrality(0),
  hTrue(),
  hTrueInAccTPC(),
  hTrueInAccTOF(),
  hTrueInAccTPCTOF(),
  hTrueInAccTPCTOFBayes(),
  hIdTPConly2s(),
  hIdTPConly3s(),
  hEffPlotsTPConly2s(),
  hEffPlotsTPConly3s(),
  hIdTOFonly2s(),
  hIdTOFonly3s(),
  hIDnSigmaComb1(),
  hIDnSigmaComb2(),
  hIDnSigmaComb3(),
  hIDBayes1(),
  hIDBayes2(),
  hIDBayes3(),
  hEffPlotsTOFonly2s(),
  hEffPlotsTOFonly3s(),
  hEffPlotsnSigmaComb1(),
  hEffPlotsnSigmaComb2(),
  hEffPlotsnSigmaComb3(),
  hEffPlotsBayes1(),
  hEffPlotsBayes2(),
  hEffPlotsBayes3(),
  fPriorsUsed(),
  fpartOfInterest(),
  fFB(),
  fRejectCheckGenName(kFALSE),
  fGenToBeKept("Hijing"),
  fPIDMomCut(0.6),
  fbayesth1(0.75),
  fbayesth2(0.80),
  fbayesth3(0.9),
  fTPCchi2Cut(-1),
  fNClustersTPCCut(-1),
  fMinTPCCrossedRows(-1),
  fMinTPCRowsOverFindableCls(-1)
{
  
  
  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
}

//_______________________________________________________________________________________
void AliAnalysisTaskPIDPerformCombIDPtDep::UserCreateOutputObjects()
{
  //
  // Initialise the framework objects
  //

  fHistList=new TList();
  fHistList->SetName("fHistList");
  fHistList->SetOwner();
  
  // ------- setup PIDCombined
  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);

   //AOD analysis
  fHistCentrality = new TH1F("fHistCentrality",";Centrality bin;Events",
			     1001,-0.5,100.5);
  fHistList->Add(fHistCentrality);

  Int_t ptBin = 36;
  Double_t nArrayPt[37]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0};
  
  hIdTPConly2s[0] = new TH1F("hIdTPConly2sPos", "hIdTPConly2sPos", ptBin, nArrayPt);
  fHistList->Add(hIdTPConly2s[0]);
  hIdTPConly2s[1] = new TH1F("hIdTPConly2sNeg", "hIdTPConly2sNeg",ptBin, nArrayPt);
  fHistList->Add(hIdTPConly2s[1]);
  
  hIdTPConly3s[0] = new TH1F("hIdTPConly3sPos", "hIdTPConly3sPos", ptBin, nArrayPt);
  fHistList->Add(hIdTPConly3s[0]);
  hIdTPConly3s[1] = new TH1F("hIdTPConly3sNeg", "hIdTPConly3sNeg", ptBin, nArrayPt);
  fHistList->Add(hIdTPConly3s[1]);
  
  hIdTOFonly2s[0] = new TH1F("hIdTOFonly2sPos", "hIdTOFonly2sPos", ptBin, nArrayPt);
  fHistList->Add(hIdTOFonly2s[0]);
  hIdTOFonly2s[1] = new TH1F("hIdTOFonly2sNeg", "hIdTOFonly2sNeg", ptBin, nArrayPt);
  fHistList->Add(hIdTOFonly2s[1]);
  
  hIdTOFonly3s[0] = new TH1F("hIdTOFonly3sPos", "hIdTOFonly3sPos", ptBin, nArrayPt);
  fHistList->Add(hIdTOFonly3s[0]);
  hIdTOFonly3s[1] = new TH1F("hIdTOFonly3sNeg", "hIdTOFonly3sNeg",ptBin, nArrayPt);
  fHistList->Add(hIdTOFonly3s[1]);

  hIDnSigmaComb1[0] = new TH1F("hIDnSigmaComb1Pos", "hIDnSigmaComb1Pos", ptBin, nArrayPt);
  fHistList->Add(hIDnSigmaComb1[0]);
  hIDnSigmaComb1[1] = new TH1F("hIDnSigmaComb1Neg", "hIDnSigmaComb1Neg", ptBin, nArrayPt);
  fHistList->Add(hIDnSigmaComb1[1]);

  hIDnSigmaComb2[0] = new TH1F("hIDnSigmaComb2Pos", "hIDnSigmaComb2Pos", ptBin, nArrayPt);
  fHistList->Add(hIDnSigmaComb2[0]);
  hIDnSigmaComb2[1] = new TH1F("hIDnSigmaComb2Neg", "hIDnSigmaComb2Neg",ptBin, nArrayPt);
  fHistList->Add(hIDnSigmaComb2[1]);

  hIDnSigmaComb3[0] = new TH1F("hIDnSigmaComb3Pos", "hIDnSigmaComb3Pos", ptBin, nArrayPt);
  fHistList->Add(hIDnSigmaComb3[0]);
  hIDnSigmaComb3[1] = new TH1F("hIDnSigmaComb3Neg", "hIDnSigmaComb3Neg", ptBin, nArrayPt);
  fHistList->Add(hIDnSigmaComb3[1]);

  hIDBayes1[0]= new TH1F("hIDBayes1Pos", "hIDBayes1Pos", ptBin, nArrayPt);
  fHistList->Add(hIDBayes1[0]);
  hIDBayes1[1]= new TH1F("hIDBayes1Neg", "hIDBayes1Neg", ptBin, nArrayPt);
  fHistList->Add(hIDBayes1[1]);

  hIDBayes2[0]= new TH1F("hIDBayes2Pos", "hIDBayes2Pos", ptBin, nArrayPt);
  fHistList->Add(hIDBayes2[0]);
  hIDBayes2[1]= new TH1F("hIDBayes2Neg", "hIDBayes2Neg", ptBin, nArrayPt);
  fHistList->Add(hIDBayes2[1]);

  hIDBayes3[0]= new TH1F("hIDBayes3Pos", "hIDBayes3Pos", ptBin, nArrayPt);
  fHistList->Add(hIDBayes3[0]);
  hIDBayes3[1]= new TH1F("hIDBayes3Neg", "hIDBayes3Neg", ptBin, nArrayPt);
  fHistList->Add(hIDBayes3[1]);
   
  
  for (Int_t ispec1=0; ispec1<5; ispec1++){

    hTrue[ispec1][0] = new TH1F(Form("hTruePosSpec%d",ispec1), Form("hTruePosSpec%d",ispec1),ptBin, nArrayPt);
    fHistList->Add(hTrue[ispec1][0]);
    
    hTrue[ispec1][1] = new TH1F(Form("hTrueNegSpec%d",ispec1), Form("hTrueNegSpec%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hTrue[ispec1][1]);

    hTrueInAccTPC[ispec1][0] = new TH1F(Form("hTrueInAccTPCPosSpec%d",ispec1), Form("hTrueInAccTPCPosSpec%d",ispec1),ptBin, nArrayPt);
    fHistList->Add(hTrueInAccTPC[ispec1][0]);
    
    hTrueInAccTPC[ispec1][1] = new TH1F(Form("hTrueInAccTPCNegSpec%d",ispec1), Form("hTrueInAccTPCNegSpec%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hTrueInAccTPC[ispec1][1]);

    hTrueInAccTOF[ispec1][0] = new TH1F(Form("hTrueInAccTOFPosSpec%d",ispec1), Form("hTrueInAccTOFPosSpec%d",ispec1),ptBin, nArrayPt);
    fHistList->Add(hTrueInAccTOF[ispec1][0]);
    
    hTrueInAccTOF[ispec1][1] = new TH1F(Form("hTrueInAccTOFNegSpec%d",ispec1), Form("hTrueInAccTOFNegSpec%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hTrueInAccTOF[ispec1][1]);
    
    hTrueInAccTPCTOF[ispec1][0] = new TH1F(Form("hTrueInAccTPCTOFPosSpec%d",ispec1), Form("hTrueInAccTPCTOFPosSpec%d",ispec1),ptBin, nArrayPt);
    fHistList->Add(hTrueInAccTPCTOF[ispec1][0]);

    hTrueInAccTPCTOF[ispec1][1] = new TH1F(Form("hTrueInAccTPCTOFNegSpec%d",ispec1), Form("hTrueInAccTPCTOFNegSpec%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hTrueInAccTPCTOF[ispec1][1]);

    hTrueInAccTPCTOFBayes[ispec1][0] = new TH1F(Form("hTrueInAccTPCTOFBayesPosSpec%d",ispec1), Form("hTrueInAccTPCTOFBayesPosSpec%d",ispec1),ptBin, nArrayPt);
    fHistList->Add(hTrueInAccTPCTOFBayes[ispec1][0]);

    hTrueInAccTPCTOFBayes[ispec1][1] = new TH1F(Form("hTrueInAccTPCTOFBayesNegSpec%d",ispec1), Form("hTrueInAccTPCTOFNegBayesSpec%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hTrueInAccTPCTOFBayes[ispec1][1]);
    
    fPriorsUsed[ispec1] = new TH2D(Form("%s_priorsUsed",AliPID::ParticleName(ispec1)),
				   Form("%s priors vs transverse momentum;p_{t} (GeV/c);priors",AliPID::ParticleName(ispec1)),
				   100,0.,20.,101,0,1.01);      
    fHistList->Add(fPriorsUsed[ispec1]);
    
    
    hEffPlotsTPConly2s[ispec1][0] = new TH1F (Form("hEffPlotsTPConly2sPosTrue%d",ispec1), Form("hEffPlotsTPConly2sPosTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsTPConly2s[ispec1][0]);
    
    hEffPlotsTPConly2s[ispec1][1] = new TH1F (Form("hEffPlotsTPConly2sNegTrue%d",ispec1), Form("hEffPlotsTPConly2sNegTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsTPConly2s[ispec1][1]);

    hEffPlotsTPConly3s[ispec1][0] = new TH1F (Form("hEffPlotsTPConly3sPosTrue%d",ispec1), Form("hEffPlotsTPConly3sPosTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsTPConly3s[ispec1][0]);

    hEffPlotsTPConly3s[ispec1][1] = new TH1F (Form("hEffPlotsTPConly3sNegTrue%d",ispec1), Form("hEffPlotsTPConly3sNegTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsTPConly3s[ispec1][1]);
    
    hEffPlotsTOFonly2s[ispec1][0] = new TH1F (Form("hEffPlotsTOFonly2sPosTrue%d",ispec1), Form("hEffPlotsTOFonly2sPosTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsTOFonly2s[ispec1][0]);
    
    hEffPlotsTOFonly2s[ispec1][1] = new TH1F (Form("hEffPlotsTOFonly2sNegTrue%d",ispec1), Form("hEffPlotsTOFonly2sNegTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add( hEffPlotsTOFonly2s[ispec1][1]);

    hEffPlotsTOFonly3s[ispec1][0] = new TH1F (Form("hEffPlotsTOFonly3sPosTrue%d",ispec1), Form("hEffPlotsTOFonly3sPosTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsTOFonly3s[ispec1][0]);

    hEffPlotsTOFonly3s[ispec1][1] = new TH1F (Form("hEffPlotsTOFonly3sNegrue%d",ispec1), Form("hEffPlotsTOFonly3sNegTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsTOFonly3s[ispec1][1]);
    
    hEffPlotsnSigmaComb1[ispec1][0] = new TH1F (Form("hEffPlotsnSigmaComb1PosTrue%d",ispec1), Form("hEffPlotsnSigmaComb1PosTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsnSigmaComb1[ispec1][0]);
    hEffPlotsnSigmaComb1[ispec1][1] = new TH1F (Form("hEffPlotsnSigmaComb1NegTrue%d",ispec1), Form("hEffPlotsnSigmaComb1NegTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsnSigmaComb1[ispec1][1]);

    hEffPlotsnSigmaComb2[ispec1][0] = new TH1F (Form("hEffPlotsnSigmaComb2PosTrue%d",ispec1), Form("hEffPlotsnSigmaComb2PosTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsnSigmaComb2[ispec1][0]);
    hEffPlotsnSigmaComb2[ispec1][1] = new TH1F (Form("hEffPlotsnSigmaComb2NegTrue%d",ispec1), Form("hEffPlotsnSigmaComb2NegTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsnSigmaComb2[ispec1][1]);

    hEffPlotsnSigmaComb3[ispec1][0] = new TH1F (Form("hEffPlotsnSigmaComb3PosTrue%d",ispec1), Form("hEffPlotsnSigmaComb3PosTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsnSigmaComb3[ispec1][0]);
    hEffPlotsnSigmaComb3[ispec1][1] = new TH1F (Form("hEffPlotsnSigmaComb3NegTrue%d",ispec1), Form("hEffPlotsnSigmaComb3NegTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsnSigmaComb3[ispec1][1]);

    hEffPlotsBayes1[ispec1][0]= new TH1F (Form("hEffPlotsBayes1PosTrue%d",ispec1), Form("hEffPlotsBayes1PosTrue%d",ispec1), ptBin, nArrayPt);
    hEffPlotsBayes1[ispec1][1]= new TH1F (Form("hEffPlotsBayes1NegTrue%d",ispec1), Form("hEffPlotsBayes1NegTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsBayes1[ispec1][0]);
    fHistList->Add(hEffPlotsBayes1[ispec1][1]);

    hEffPlotsBayes2[ispec1][0]= new TH1F (Form("hEffPlotsBayes2PosTrue%d",ispec1), Form("hEffPlotsBayes2PosTrue%d",ispec1), ptBin, nArrayPt);
    hEffPlotsBayes2[ispec1][1]= new TH1F (Form("hEffPlotsBayes2NegTrue%d",ispec1), Form("hEffPlotsBayes2NegTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsBayes2[ispec1][0]);
    fHistList->Add(hEffPlotsBayes2[ispec1][1]);
    
    hEffPlotsBayes3[ispec1][0]= new TH1F (Form("hEffPlotsBayes3PosTrue%d",ispec1), Form("hEffPlotsBayes3PosTrue%d",ispec1), ptBin, nArrayPt);
    hEffPlotsBayes3[ispec1][1]= new TH1F (Form("hEffPlotsBayes3NegTrue%d",ispec1), Form("hEffPlotsBayes3NegTrue%d",ispec1), ptBin, nArrayPt);
    fHistList->Add(hEffPlotsBayes3[ispec1][0]);
    fHistList->Add(hEffPlotsBayes3[ispec1][1]);
    
  }

  
  fHistList->SetOwner();
  PostData(1,fHistList);
  
}


//_________________________________________________________________________________
void AliAnalysisTaskPIDPerformCombIDPtDep::UserExec(Option_t *)
{
  //
  // Main loop. Called for every event
  //

  //Printf("main loop starts");
  
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;
  }
  
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    AliError("ERROR: Could not retrieve MC event");
    return;
  }

  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");

   AliAODHeader *headerAOD = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
    if (!headerAOD){
      AliFatal("AOD header found");
      return;
    }

    AliMultSelection *multSelection;

    //Centrality stuff
    Double_t nCentrality = 0;
    if(fUseCentrality){
      if (fAOD->GetRunNumber()<244824) {
	AliCentrality *centrality = headerAOD->GetCentralityP();
	nCentrality =centrality->GetCentralityPercentile(fCentralityEstimator.Data());
	if(!centrality->IsEventInCentralityClass(fCentralityPercentileMin,
						 fCentralityPercentileMax,
						 fCentralityEstimator.Data()))
	  return;
      }
      
      else {
	multSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
          
    if(!multSelection) {
            AliWarning("AliMultSelection object not found!");
    }
          
	nCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator, kTRUE);

	if ((nCentrality < fCentralityPercentileMin) || (nCentrality >= fCentralityPercentileMax)) return;
	fHistCentrality->Fill(nCentrality);
	
      }

    }
    
  
  Double_t probTPC[AliPID::kSPECIES]={0.};
  Double_t probTOF[AliPID::kSPECIES]={0.};
  Double_t probTPCTOF[AliPID::kSPECIES]={0.};
  
  /*Double_t nSigmaTPCTPConly[AliPID::kSPECIES]={0.};
  Double_t nSigmaTOFTOFonly[AliPID::kSPECIES]={0.};
  Double_t nSigmaTPCTPCTOFcomb[AliPID::kSPECIES]={0.};
  Double_t nSigmaTOFTPCTOFcomb[AliPID::kSPECIES]={0.};
  */
  
  Double_t nSigmaProtonsTPCOnly=0., nSigmaProtonsTOFOnly=0., nSigmaProtonsTPCNsigcomb=0, nSigmaProtonsTOFNsigcomb=0.;
  Double_t combSquaredSigma=0;
  
  AliAODTrack *track=0x0;
  Int_t ntracks=fAOD->GetNumberOfTracks();

  Int_t tofLabel[3] = {0,0,0};
  
  // Printf("ntracks =%d", ntracks);
   //loop over all tracks
  for (Int_t itrack=0; itrack<ntracks; ++itrack){

    track=(AliAODTrack*)fAOD->GetTrack(itrack);
    //Printf("itrack = %d, ntracks=%d,  track=%p",itrack, ntracks, track);
	
    if (!track) {
	AliError(Form("Could not receive track %d", itrack));
	continue;
    }
    
    if(track->TestFilterBit(fFB)){

      if( fTPCchi2Cut != -1 && track->Chi2perNDF() > fTPCchi2Cut){
	continue;
      }
      if( fNClustersTPCCut != -1 && track->GetTPCNcls() < fNClustersTPCCut){
	continue;
      }
      
      if(fMinTPCCrossedRows != -1){
	if ((Float_t)track->GetTPCNCrossedRows() < (120 - (5/(Float_t)track->Pt())) ){
	  continue;
	}
      }
      
      if (fMinTPCRowsOverFindableCls != -1){
	Float_t nTPCCrossedRowsOverFindCls = (((Float_t)track->GetTPCNCrossedRows())/((Float_t)track->GetTPCNclsF()));
	if (nTPCCrossedRowsOverFindCls < fMinTPCRowsOverFindableCls){
	  continue;
	}
      }
      
      //Printf("fb test quiiiiii");
            
      Double_t mom=track->GetTPCmomentum();
      Double_t pt=track->Pt();
      Int_t ibin=GetMomBin(mom);
      Short_t charge = track->Charge();
      Int_t chargeBin; 
      if (charge > 0) chargeBin=0;
      else chargeBin = 1;
      //  Printf("track =%p, mom =%f, ibin =%d", track, track->GetTPCmomentum(), GetMomBin(track->GetTPCmomentum()));

      
     
      AliAODPid* pidObj = track->GetDetPid();
      
      //MC part
      Int_t label = track->GetLabel();
      if((label >  mcEvent->GetNumberOfTracks()) || (label < 0))  continue;
      AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) mcEvent->GetTrack(label); 
      Short_t gAODmcCharge = AODmcTrack->Charge();////
      Int_t pdgcode = AODmcTrack->GetPdgCode();
      //Printf("AODmcTrack=%p, pdgcode =%d, gAODmcCharge=%d", AODmcTrack, pdgcode, gAODmcCharge);
      Int_t pdgIndex =-1;
     
        
      if (fRejectCheckGenName){
      TString generatorName;
      Bool_t hasGenerator = mcEvent->GetCocktailGenerator(label,generatorName);
      if((!hasGenerator) || (!generatorName.Contains(fGenToBeKept.Data())))
        continue;
      }
        
      if (!(AODmcTrack->IsPhysicalPrimary())) continue;

      Double_t y = AODmcTrack->Y();
      if (TMath::Abs(y) > 0.5) continue;
      
      //Printf("mother =%d, generatorName=%s", label, generatorName.Data()); 
      
      
      // if (pdgcode*charge < 0) Printf(" ********** charge mismatch Reco / Gen ********");
      //check PID strategy
      switch (TMath::Abs(pdgcode)) {
      case 11:
	pdgIndex = 0; //e
	hTrue[pdgIndex][chargeBin]->Fill(pt);
	break;
      case 13: 
	pdgIndex = 1; //mu 
	hTrue[pdgIndex][chargeBin]->Fill(pt);
	break;
      case 211:
	pdgIndex = 2; //pi 
	hTrue[pdgIndex][chargeBin]->Fill(pt);
	break;	
      case 321:
	pdgIndex = 3; //K
	hTrue[pdgIndex][chargeBin]->Fill(pt);
	break;	
      case 2212:
	pdgIndex = 4; //p
	hTrue[pdgIndex][chargeBin]->Fill(pt);
	break;

      }

      if (pdgIndex ==-1) {
	Printf("none of the considered species");
	continue; 
      }
      
      // Printf("pdgIndex=%d", pdgIndex);

      
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
      UInt_t detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
      
      Double_t nSigmaTPCPions = 0.;
      Double_t nSigmaTPCKaons = 0.;
      Double_t nSigmaTPCProtons = 0.;
      Double_t nSigmaTOFPions = 0.;
      Double_t nSigmaTOFKaons = 0.;
      Double_t nSigmaTOFProtons = 0.;
      Double_t nSigmaTPCTOFPions = 0.;
      Double_t nSigmaTPCTOFKaons = 0.;
      Double_t nSigmaTPCTOFProtons = 0.;
      //Printf("detUsed=%d", detUsed);
      
      if (detUsed  == (UInt_t)fPIDCombined->GetDetectorMask() ) {  // TPC is available
	if (pdgIndex > 4) continue;
        
	hTrueInAccTPC[pdgIndex][chargeBin]->Fill(pt);
        
	nSigmaProtonsTPCOnly = fPIDResponse->NumberOfSigmasTPC(track,fpartOfInterest);
	nSigmaTPCPions = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
	nSigmaTPCKaons = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
	nSigmaTPCProtons = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));
        
	if(TMath::Abs(nSigmaProtonsTPCOnly)<2) {
	  hIdTPConly2s[chargeBin]->Fill(pt);
	  hEffPlotsTPConly2s[pdgIndex][chargeBin]->Fill(pt);
	}
        
	if (TMath::Abs(nSigmaProtonsTPCOnly) < 3)  {
	  // Printf("inside 3 sigma");
	  hIdTPConly3s[chargeBin]->Fill(pt);
	  hEffPlotsTPConly3s[pdgIndex][chargeBin]->Fill(pt);
	}
        
	//if(mom < fPIDMomCut){
	if(pt < fPIDMomCut){
	  
	  hTrueInAccTPCTOF[pdgIndex][chargeBin]->Fill(pt);
                
	  if (fpartOfInterest == (AliPID::kPion)){
	    
	    if (TMath::Abs(nSigmaProtonsTPCOnly)<2 && !(TMath::Abs(nSigmaTPCKaons)<3.) && !(TMath::Abs(nSigmaTPCProtons)<3.)) {
	      hIDnSigmaComb1[chargeBin]->Fill(pt);
	      hEffPlotsnSigmaComb1[pdgIndex][chargeBin]->Fill(pt);
              
	    }
            
	    if (TMath::Abs(nSigmaProtonsTPCOnly)<3 && !(TMath::Abs(nSigmaTPCKaons)<3.) && !(TMath::Abs(nSigmaTPCProtons)<3.)) {
	      hIDnSigmaComb2[chargeBin]->Fill(pt);
	      hEffPlotsnSigmaComb2[pdgIndex][chargeBin]->Fill(pt);
	    }
	    
	    if (TMath::Abs(nSigmaProtonsTPCOnly)<2) {
	      hIDnSigmaComb3[chargeBin]->Fill(pt);
	      hEffPlotsnSigmaComb3[pdgIndex][chargeBin]->Fill(pt);
	    }
            
            
	  }//end of pions
          
	  if (fpartOfInterest == (AliPID::kKaon)){
	    
	    if (TMath::Abs(nSigmaProtonsTPCOnly)<2 && !(TMath::Abs(nSigmaTPCPions)<3.) && !(TMath::Abs(nSigmaTPCProtons)<3.)) {
	      
	      hIDnSigmaComb1[chargeBin]->Fill(pt);
	      hEffPlotsnSigmaComb1[pdgIndex][chargeBin]->Fill(pt);
	    }
            
            
	    if (TMath::Abs(nSigmaProtonsTPCOnly)<3 && !(TMath::Abs(nSigmaTPCPions)<3.) && !(TMath::Abs(nSigmaTPCProtons)<3.)) {
	      
	      hIDnSigmaComb2[chargeBin]->Fill(pt);
	      hEffPlotsnSigmaComb2[pdgIndex][chargeBin]->Fill(pt);
	    }
            
	    if (TMath::Abs(nSigmaProtonsTPCOnly)<2) {
	      
	      hIDnSigmaComb3[chargeBin]->Fill(pt);
	      hEffPlotsnSigmaComb3[pdgIndex][chargeBin]->Fill(pt);
	    }
            
	  }//end of kaons
          
	  if (fpartOfInterest == (AliPID::kProton)){
	    
	    if (TMath::Abs(nSigmaProtonsTPCOnly)<2 && !(TMath::Abs(nSigmaTPCPions)<3.) && !(TMath::Abs(nSigmaTPCKaons)<3.)) {
	      hIDnSigmaComb1[chargeBin]->Fill(pt);
	      hEffPlotsnSigmaComb1[pdgIndex][chargeBin]->Fill(pt);
	    }
            
	    if (TMath::Abs(nSigmaProtonsTPCOnly)<3 && !(TMath::Abs(nSigmaTPCPions)<3.) && !(TMath::Abs(nSigmaTPCKaons)<3.)) {
	      
	      hIDnSigmaComb2[chargeBin]->Fill(pt);
	      hEffPlotsnSigmaComb2[pdgIndex][chargeBin]->Fill(pt);
	    }
            
	    if (TMath::Abs(nSigmaProtonsTPCOnly)<2) {
	      
	      hIDnSigmaComb3[chargeBin]->Fill(pt);
	      hEffPlotsnSigmaComb3[pdgIndex][chargeBin]->Fill(pt);
	    }
            
	  }//end of protons
	  
          
	}
        
        
	//Printf("after TPC");
        
	//TOF ONLY
	fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
	detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
        
	Double_t priors[5];     // check priors used for TOF
	fPIDCombined->GetPriors(track,priors,fPIDResponse,detUsed);
	for(Int_t ispec=0;ispec<5;ispec++) fPriorsUsed[ispec]->Fill(TMath::Abs(track->Pt()),priors[ispec]);
        
	if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()) {  // TOF is available
	  if (pdgIndex > 4) continue;
          
	  Int_t detStatus = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,track);
	  // Printf("CheckPIDStatus(kDetTPC) =%d, CheckPIDStatus(kDetTOF)=%d", fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track),fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track) );
          
          
	  if(pidObj && pidObj->GetTOFsignal() < 99999){
	    // track->GetTOFLabel(tofLabel);
	    //printf("%i compared to (%i,%i,%i) = %i\n",track->GetLabel(),tofLabel[0],tofLabel[1],tofLabel[2],TMath::Abs(track->GetLabel()) == TMath::Abs(tofLabel[0]));
	    //if (track->GetLabel() != tofLabel[0]) continue;
	    // if ((fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track)) !=  (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track))) continue;
            
	    hTrueInAccTOF[pdgIndex][chargeBin]->Fill(pt);
            
	    nSigmaProtonsTOFOnly = fPIDResponse->NumberOfSigmasTOF(track,fpartOfInterest);
	    
	    if (TMath::Abs(nSigmaProtonsTOFOnly)< 2){
	      hIdTOFonly2s[chargeBin]->Fill(pt);
	      hEffPlotsTOFonly2s[pdgIndex][chargeBin]->Fill(pt);
	    }
	    if (TMath::Abs(nSigmaProtonsTOFOnly) < 3){
	      hIdTOFonly3s[chargeBin]->Fill(pt);
	      hEffPlotsTOFonly3s[pdgIndex][chargeBin]->Fill(pt);
	    }
	  }
	}
	//Printf("after TOF");
        
        
	fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
	detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPCTOF);
        
	if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()){
	  if (pdgIndex > 4) continue;
          
	  if ((detUsed >= AliPIDResponse::kDetTOF) && (pidObj && pidObj->GetTOFsignal() < 99999)){
	    // track->GetTOFLabel(tofLabel);
	    //if (track->GetLabel() != tofLabel[0]) continue;
            
            
	    nSigmaProtonsTPCNsigcomb = fPIDResponse->NumberOfSigmasTPC(track,fpartOfInterest);
	    nSigmaProtonsTOFNsigcomb = fPIDResponse->NumberOfSigmasTOF(track,fpartOfInterest);
	    combSquaredSigma = TMath::Sqrt((nSigmaProtonsTPCNsigcomb*nSigmaProtonsTPCNsigcomb) + (nSigmaProtonsTOFNsigcomb*nSigmaProtonsTOFNsigcomb));
            
	    nSigmaTOFPions = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
	    nSigmaTOFKaons = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
	    nSigmaTOFProtons = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);
            
	    nSigmaTPCTOFPions = TMath::Sqrt(nSigmaTPCPions*nSigmaTPCPions + nSigmaTOFPions*nSigmaTOFPions);
	    nSigmaTPCTOFKaons = TMath::Sqrt(nSigmaTPCKaons*nSigmaTPCKaons + nSigmaTOFKaons*nSigmaTOFKaons);
	    nSigmaTPCTOFProtons = TMath::Sqrt(nSigmaTPCProtons*nSigmaTPCProtons + nSigmaTOFProtons*nSigmaTOFProtons);
            
	    //if (mom > fPIDMomCut){
	    if (pt > fPIDMomCut){
	      
	      hTrueInAccTPCTOF[pdgIndex][chargeBin]->Fill(pt);
              
	      if (fpartOfInterest == (AliPID::kPion)){
		//first nsigma strategy
                
		if (TMath::Abs(combSquaredSigma)<2. && !(TMath::Abs(nSigmaTPCTOFKaons)<3.) && !(TMath::Abs(nSigmaTPCTOFProtons)<3.)) {
		  
		  hIDnSigmaComb1[chargeBin]->Fill(pt);
		  hEffPlotsnSigmaComb1[pdgIndex][chargeBin]->Fill(pt);
		}
                
		if (TMath::Abs(combSquaredSigma)<3. && !(TMath::Abs(nSigmaTPCTOFKaons)<3.) && !(TMath::Abs(nSigmaTPCTOFProtons)<3.)) {
		  
		  hIDnSigmaComb2[chargeBin]->Fill(pt);
		  hEffPlotsnSigmaComb2[pdgIndex][chargeBin]->Fill(pt);
		}
                
		if (TMath::Abs(combSquaredSigma)<2.) {
		  
		  hIDnSigmaComb3[chargeBin]->Fill(pt);
		  hEffPlotsnSigmaComb3[pdgIndex][chargeBin]->Fill(pt);
		}
                
	      }//end of pions
              
	      if (fpartOfInterest == (AliPID::kKaon)){
		//first nsigma strategy
                
		if (TMath::Abs(combSquaredSigma)<2. && !(TMath::Abs(nSigmaTPCTOFPions)<3.) && !(TMath::Abs(nSigmaTPCTOFProtons)<3.)) {
		  
		  hIDnSigmaComb1[chargeBin]->Fill(pt);
		  hEffPlotsnSigmaComb1[pdgIndex][chargeBin]->Fill(pt);
		}
                
		if (TMath::Abs(combSquaredSigma)<2.) {
		  
		  hIDnSigmaComb3[chargeBin]->Fill(pt);
		  hEffPlotsnSigmaComb3[pdgIndex][chargeBin]->Fill(pt);
		}
                
		if (pt<=2){
		  
		  if (TMath::Abs(combSquaredSigma)<2) {
		    
		    hIDnSigmaComb2[chargeBin]->Fill(pt);
		    hEffPlotsnSigmaComb2[pdgIndex][chargeBin]->Fill(pt);
                    
		  }
                  
		}
                
		else if (pt>2){
		  
		  if (TMath::Abs(combSquaredSigma)<1.5) {
		    
		    hIDnSigmaComb2[chargeBin]->Fill(pt);
		    hEffPlotsnSigmaComb2[pdgIndex][chargeBin]->Fill(pt);
		  }
		}
	      }//end of kaons
              
	      if (fpartOfInterest == (AliPID::kProton)){
		//first nsigma strategy
                
		if (TMath::Abs(combSquaredSigma)<2. && !(TMath::Abs(nSigmaTPCTOFPions)<3.) && !(TMath::Abs(nSigmaTPCTOFKaons)<3.)) {
		  hIDnSigmaComb1[chargeBin]->Fill(pt);
		  hEffPlotsnSigmaComb1[pdgIndex][chargeBin]->Fill(pt);
		}
                
		if (TMath::Abs(combSquaredSigma)<2.) {
		  hIDnSigmaComb3[chargeBin]->Fill(pt);
		  hEffPlotsnSigmaComb3[pdgIndex][chargeBin]->Fill(pt);
		}
                
		if (pt<=2){
		  
		  if (TMath::Abs(nSigmaProtonsTOFNsigcomb)<3) {
		    
		    hIDnSigmaComb2[chargeBin]->Fill(pt);
		    hEffPlotsnSigmaComb2[pdgIndex][chargeBin]->Fill(pt);
                    
		  }
		}
                
		else if (pt>2){
		  
		  if (TMath::Abs(combSquaredSigma)<1.5) {
		    
		    hIDnSigmaComb2[chargeBin]->Fill(pt);
		    hEffPlotsnSigmaComb2[pdgIndex][chargeBin]->Fill(pt);
		  }
		}
                
                
	      }//end of protons
              
	    }
            
	  }
          
	}
      }
      
      
    }
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskPIDPerformCombIDPtDep::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  
  fHistList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fHistList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//_________________________________________________________________________________
Int_t AliAnalysisTaskPIDPerformCombIDPtDep::GetMomBin(Float_t mom)
{
  //
  // Given momentum return histogram to be filled
  //
  if (mom>0. && mom < 0.5) return 0;
  if (mom>=0.5 && mom < 0.7) return 1;
  if (mom>=0.7 && mom < 1.0) return 2;
  if (mom>=1.0 && mom < 1.5) return 3;
  if (mom>=1.5 && mom < 2.0) return 4;
  if (mom>=2.0 && mom < 2.5) return 5;
  if (mom>=2.5 && mom < 3.0) return 6;
  if (mom>=3.0 && mom < 4.0) return 7;
  return kPtBins-1;
}

