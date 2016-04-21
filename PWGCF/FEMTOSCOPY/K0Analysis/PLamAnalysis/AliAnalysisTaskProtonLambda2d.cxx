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
#if !defined(__CINT__) || defined(__MAKECINT__)
// This only for low level memory profiling,
// needs Root.ObjectStat in .rootrc
// #include <TObjectTable.h>

#include <TH1F.h>
// #include <TH2F.h>
#include <TH3F.h>

#include <TAxis.h>
// #include "TObjArray.h"

#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>

#include <AliAODEvent.h>
#include <AliAODVertex.h>
#include <AliAODv0.h>
#include <AliAODInputHandler.h>

#include "AliAnalysisTaskProtonLambda2d.h"
#include <AliCentrality.h>
//#include "AliAODpid.h"
#include <AliPID.h>
#include <AliPIDResponse.h>
// #include <../STEER/STEER/AliV0.h>
#include <AliExternalTrackParam.h>
//#include <AliAODTrack.h>
//#include <AliESDtrack.h>

//#include "EventCollection.h"

// Task to study femtoscopic proton-lambda correlations
// Author: Hans Beck, haarigerhans at gmail dot com

ClassImp(AliAnalysisTaskProtonLambda2d)
//ClassImp(AliAnalysisTaskProtonLambda2d::GlobalTrackInfo)
//ClassImp(AliAnalysisTaskProtonLambda2d::GTIContainer)

#endif  
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::AliAnalysisTaskProtonLambda2d() 
  : AliAnalysisTaskSE(),
  fkDoStdHists(kFALSE),
  fkDoDCAHists(kFALSE),
  fkDoBgLamALam(kFALSE),
  fkUseOnTheFly(kTRUE),
  fkAbsZvertexCut(10.0),
  fkCentCutLo(0),
  fkCentCutHi(10.0),
  fkMinvCut(.004),
  fkCentEst("V0M"),
//  fdoLamOnly(kFALSE),
  fdEtaSCut(.147),
  fdPhiSCut(.0283),
  fkLamMass(1.115683),
  fkProMass(0.9382720),
  fkPioMass(0.13957018),
  fkDCAFlag(kStd), // kStd=0,kWide=1,kStrict=2
  fPIDResponse(0), 
  fTpcResponse(0),
  fFemtoBuffer(0),
  fAOD(0), fPrimaryVtx(0), fOutputList(0), fOutputPrimaries(0),
  fOutput2Part(0),
  fGTI(0),fTrackBuffSize(19000),
  fLamPur(),fALamPur(), 
  fLamFdPur010(),fLamFdPur1020(),fLamFdPur2040(),fLamFdPur4060(),
  fALamFdPur010(),fALamFdPur1020(),fALamFdPur2040(),fALamFdPur4060(),
  fProFdPur(),fAProFdPur(),
// Histograms to be filled
  fHistGoodEvent(0),
  fHistCentComp(0),
  fHistSideBandOffLam(0), fHistSideBandOffALam(0), 
  fHistMassLambdaOff(0), fHistMassAntiLambdaOff(0),        
  fHistYPtMassLamOff(0), fHistYPtMassALamOff(0),
  fHistSideBandOnLam(0), fHistSideBandOnALam(0),
  fHistMassLambdaOn(0),fHistMassAntiLambdaOn(0),
  fHistYPtMassLamOn(0),fHistYPtMassALamOn(0),
  fPriHistShare(0),
  fPriHistTOFsignalPosVsP(0),      
  fPriHistTOFsignalNegVsP(0),
  fPriHistHybridTOFsigPosWoTPC(0), 
  fPriHistHybridTOFsigPosTPCok(0),fPriHistHybridTOFsigNegWoTPC(0),fPriHistHybridTOFsigNegTPCok(0),
  fPriHistTPCsignalPos(0),
  fPriHistTPCsignalLowPPos(0),fPriHistTPCsignalMedPPos(0),fPriHistTPCsignalHigPPos(0),
  fPriHistTPCsignalNeg(0),
  fPriHistTPCsignalLowPNeg(0),fPriHistTPCsignalMedPNeg(0),fPriHistTPCsignalHigPNeg(0),
  fPriHistDCAxyzYPtPro(0),fPriHistDCAxyzYPtAPro(0),

// Monitor needed buffer size
  f2HistNLamBefClean(0),f2HistNProBefClean(0),
  f2HistNALamBefClean(0),f2HistNAProBefClean(0),

// 2d angular distance at R=1.2
  f2HistLamProAngDistSft2dAtR12Real(0),
  f2HistALamAProAngDistSft2dAtR12Real(0),
  f2HistBgLamProAngDistSft2dAtR12Real(0),
  f2HistBgALamAProAngDistSft2dAtR12Real(0),
  
  f2HistLamProAngDistSft2dAtR12Mixed(0),
  f2HistALamAProAngDistSft2dAtR12Mixed(0),
  f2HistBgLamProAngDistSft2dAtR12Mixed(0),
  f2HistBgALamAProAngDistSft2dAtR12Mixed(0),
  f2HistMtLamProReal(0), 
  f2HistMtALamAProReal(0),
  f2HistMtLowQLamProReal(0), 
  f2HistMtLowQALamAProReal(0),
  f2HistMtVsPhiSLowQ(0),
  fLamProReal(0),fALamAProReal(0),
// purities vs qinv
  f2HistLamProWLamPurReal(0),f2HistLamProWLamFdPurReal(0),
  f2HistLamProWProFdPurReal(0),f2HistLamProWoPurReal(0),
  f2HistALamAProWALamPurReal(0),f2HistALamAProWALamFdPurReal(0),
  f2HistALamAProWAProFdPurReal(0),f2HistALamAProWoPurReal(0),
// purities vs qinv vs mt
  f2HistLamProWLamPurmTReal(0),f2HistLamProWLamFdPurmTReal(0),
  f2HistLamProWProFdPurmTReal(0),f2HistLamProWoPurmTReal(0),
  f2HistALamAProWALamPurmTReal(0),f2HistALamAProWALamFdPurmTReal(0),
  f2HistALamAProWAProFdPurmTReal(0),f2HistALamAProWoPurmTReal(0),
// purity fluctuations
  f2HistPurFlucLamPro(0),f2HistPurFlucALamAPro(0),
  f2HistLamYPtPair(0), f2HistALamYPtPair(0),
  f2HistProYPtPair(0), f2HistAProYPtPair(0),
  f2HistYKtPair(0), f2HistYKtAPair(0),
  f2HistLamYPtHiPair(0), f2HistALamYPtHiPair(0),
  fBgLamProReal(0),fBgALamAProReal(0),
  fLamProMixed(0),fALamAProMixed(0),
  fBgLamProMixed(0),fBgALamAProMixed(0)
{
  // Dummy constructor
  fPrimaryVtxPosition[0]=0;
  fPrimaryVtxPosition[1]=0;
  fPrimaryVtxPosition[2]=0;
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::AliAnalysisTaskProtonLambda2d(const char *name
							     ,const Float_t centCutLo
							     ,const Float_t centCutHi
							     ,const Float_t minvCut
							     ,const Bool_t doBgLam
							     ,const Bool_t useOnTheFly
							     ,const TString centEst // V0M, TRK
							     ,const Int_t dcaFlag)  // kStd=0,kWide=1,kStrict=2
  : AliAnalysisTaskSE(name),
    fkDoStdHists(kFALSE),
    fkDoDCAHists(kFALSE),
    fkDoBgLamALam(doBgLam),
    fkUseOnTheFly(useOnTheFly),
    fkAbsZvertexCut(10.0),
    fkCentCutLo(centCutLo),
    fkCentCutHi(centCutHi),
    fkMinvCut(minvCut),
    fkCentEst(centEst.Data()),
    //  fdoLamOnly(kFALSE),
  fdEtaSCut(.147),
  fdPhiSCut(.0283),
  fkLamMass(1.115683),
  fkProMass(0.9382720),
  fkPioMass(0.13957018),
    fkDCAFlag(dcaFlag),
  fPIDResponse(0), 
  fTpcResponse(0),
  fFemtoBuffer(0),
  fAOD(0), fPrimaryVtx(0), fOutputList(0), fOutputPrimaries(0),
  fOutput2Part(0),
  fGTI(0),fTrackBuffSize(19000),
    fLamPur(),fALamPur(),
    fLamFdPur010(),fLamFdPur1020(),fLamFdPur2040(),fLamFdPur4060(),
    fALamFdPur010(),fALamFdPur1020(),fALamFdPur2040(),fALamFdPur4060(),
    fProFdPur(),fAProFdPur(),
  fHistGoodEvent(0),
  fHistCentComp(0),
  fHistSideBandOffLam(0), fHistSideBandOffALam(0),
  fHistMassLambdaOff(0), fHistMassAntiLambdaOff(0),        
  fHistYPtMassLamOff(0), fHistYPtMassALamOff(0),
  fHistSideBandOnLam(0), fHistSideBandOnALam(0),
  fHistMassLambdaOn(0),fHistMassAntiLambdaOn(0),
  fHistYPtMassLamOn(0),fHistYPtMassALamOn(0),
  fPriHistShare(0),
    fPriHistTOFsignalPosVsP(0),      
    fPriHistTOFsignalNegVsP(0),
    fPriHistHybridTOFsigPosWoTPC(0), 
  fPriHistHybridTOFsigPosTPCok(0),fPriHistHybridTOFsigNegWoTPC(0),fPriHistHybridTOFsigNegTPCok(0),
  fPriHistTPCsignalPos(0),
  fPriHistTPCsignalLowPPos(0),fPriHistTPCsignalMedPPos(0),fPriHistTPCsignalHigPPos(0),
  fPriHistTPCsignalNeg(0),
  fPriHistTPCsignalLowPNeg(0),fPriHistTPCsignalMedPNeg(0),fPriHistTPCsignalHigPNeg(0),
  fPriHistDCAxyzYPtPro(0),fPriHistDCAxyzYPtAPro(0),

  // Monitor needed buffer size
  f2HistNLamBefClean(0),f2HistNProBefClean(0),
  f2HistNALamBefClean(0),f2HistNAProBefClean(0),
    
    // 2d angular distance at r=1.2m with shifted vertices 
  f2HistLamProAngDistSft2dAtR12Real(0),
  f2HistALamAProAngDistSft2dAtR12Real(0),
  f2HistBgLamProAngDistSft2dAtR12Real(0),
  f2HistBgALamAProAngDistSft2dAtR12Real(0),
  
  f2HistLamProAngDistSft2dAtR12Mixed(0),
  f2HistALamAProAngDistSft2dAtR12Mixed(0),
  f2HistBgLamProAngDistSft2dAtR12Mixed(0),
  f2HistBgALamAProAngDistSft2dAtR12Mixed(0),

  f2HistMtLamProReal(0), 
  f2HistMtALamAProReal(0),
  f2HistMtLowQLamProReal(0), 
  f2HistMtLowQALamAProReal(0),
  f2HistMtVsPhiSLowQ(0),
    // C2: mt vs eta* vs phi* qinv
  fLamProReal(0),fALamAProReal(0),
    // purities vs qinv
  f2HistLamProWLamPurReal(0),f2HistLamProWLamFdPurReal(0),
  f2HistLamProWProFdPurReal(0),f2HistLamProWoPurReal(0),
  f2HistALamAProWALamPurReal(0),f2HistALamAProWALamFdPurReal(0),
    f2HistALamAProWAProFdPurReal(0),f2HistALamAProWoPurReal(0),
    // purities vs qinv vs mt
    f2HistLamProWLamPurmTReal(0),f2HistLamProWLamFdPurmTReal(0),
    f2HistLamProWProFdPurmTReal(0),f2HistLamProWoPurmTReal(0),
    f2HistALamAProWALamPurmTReal(0),f2HistALamAProWALamFdPurmTReal(0),
    f2HistALamAProWAProFdPurmTReal(0),f2HistALamAProWoPurmTReal(0),
    // purity fluc
  f2HistPurFlucLamPro(0),f2HistPurFlucALamAPro(0),
    // phase space
  f2HistLamYPtPair(0), f2HistALamYPtPair(0),
  f2HistProYPtPair(0), f2HistAProYPtPair(0),
  f2HistYKtPair(0), f2HistYKtAPair(0),
  f2HistLamYPtHiPair(0), f2HistALamYPtHiPair(0),
  fBgLamProReal(0),fBgALamAProReal(0),
  fLamProMixed(0),fALamAProMixed(0),
  fBgLamProMixed(0),fBgALamAProMixed(0)
{
  // Constructor
  fPrimaryVtxPosition[0]=0;
  fPrimaryVtxPosition[1]=0;
  fPrimaryVtxPosition[2]=0;

  // Define output slots only here
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::~AliAnalysisTaskProtonLambda2d() {
  // Destructor, go through the data member and delete them

  // fPIDResponse is just a pointer to the pid response task,
  // we don't create it so we don't delete it. It comes from 
  // the AliInputEventHandler

  if (fTpcResponse){
    delete fTpcResponse;
    fTpcResponse=0;
  }
  if(fFemtoBuffer){
    delete fFemtoBuffer;
    fFemtoBuffer=0;
  }
  // fAOD also just comes from a function of the AliAnalysisTaskSE
  // fPrimaryVtx comes from the fAOD

  // The lists containing the histograms
  if (fOutputList){
    fOutputList->Delete();
    delete fOutputList;
    fOutputList=0;
  }
  if (fOutputPrimaries){
    fOutputPrimaries->Delete();
    delete fOutputPrimaries;
    fOutputPrimaries=0;
  }
  if (fOutput2Part){
    fOutput2Part->Delete();
    delete fOutput2Part;
    fOutput2Part=0;
  }

  // Arrays of pointers to DCAxyz hists
  // First we have to check that the arrays were created
  if(fPriHistDCAxyzYPtPro){
    for(UChar_t iY(0);iY<6;iY++){
      if(fPriHistDCAxyzYPtPro[iY])
	delete[]fPriHistDCAxyzYPtPro[iY];
      fPriHistDCAxyzYPtPro[iY]=0;
    }
    delete[]fPriHistDCAxyzYPtPro;
    fPriHistDCAxyzYPtPro=0;
  }
  if(fPriHistDCAxyzYPtAPro){
    for(UChar_t iY(0);iY<6;iY++){
      if(fPriHistDCAxyzYPtAPro[iY])
	delete[]fPriHistDCAxyzYPtAPro[iY];
      fPriHistDCAxyzYPtAPro[iY]=0;
    }
    delete[]fPriHistDCAxyzYPtAPro;
    fPriHistDCAxyzYPtAPro=0;
  }
  // Array, note the [] with the delete
  if (fGTI)
    delete[] fGTI;
  fGTI=0;
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::UserCreateOutputObjects()
{
  // Create histograms and other objects and variables
  // Called once

  // Get the PID response object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(!man){AliError("Couldn't get the analysis manager!");}
  else{
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    if(!inputHandler){AliError("Couldn't get the input handler!");}
    else{
      fPIDResponse = inputHandler->GetPIDResponse();
      if(!fPIDResponse){AliError("Couldn't get the PID response task!");}
    }
  }
  // Create dE/dx spectra cut. use it by calling
  // fTpcResponse->GetExpectedSignal(mom, AliPID::kProton)
  // Use a factor of 0.95 for 11h data!
  fTpcResponse = new AliTPCPIDResponse();
  // Set parameters of the AlephParametrization.
  // They are only valid for data, see $ALICE_ROOT/PWG2/SPECTRA/AliProtonAnalysisBase.cxx
  // for monte carlo parameters
  fTpcResponse->SetBetheBlochParameters(0.0283086,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00);
  
  // Create the buffer for event mixing
  // Standard values are
  //  fkZvertexBins(10),
  //  fkCentBins(10),
  //  fkMixBuff(5),
  //  fkPriTrackLim(100),
  //  fkV0Lim(50),
  //  fFemtoBuffer = new FemtoBuffer(10,10,5,100,50,fkAbsZvertexCut,fkCentCutHi);
  // const UChar_t ZvertexBins,const UChar_t CentBins,const UChar_t MixBuff,const UShort_t PriTrackLim,const UShort_t V0Lim, const Float_t AbsZvertexCut,const Float_t CentCutLo,const Float_t CentCutHi

  // Do mixing in 10% centrality bins, calculate from  
  // fkCentCutLo and fkCentCutHi
  // Do the buffer only if we need it
  //  if(!fdoLamOnly){
  const UChar_t nCentBins = (UChar_t) (((fkCentCutHi-fkCentCutLo) + .1) /10.);
  const UChar_t priTrackLim(40),V0Lim(25); // Checked no need for more than 25 for lam but more than 30 for protons!
  fFemtoBuffer = new FemtoBuffer(4,nCentBins,8,priTrackLim,V0Lim,fkAbsZvertexCut,fkCentCutLo,fkCentCutHi,fkDoBgLamALam);
    //  }

  // In AODs, TPC only tracks don't have the pid information stored.
  // Also, the TPC only tracks don't have any resolution in the DCAxy
  // to distinguish between primaries and secondaries so we need the
  // corresponding global track for every TPC only track. The way to do 
  // this is to just store the pointer to the global track for every id.
  fGTI = new const AliAODTrack *[fTrackBuffSize]; // Array of pointers 

  // Create the output list
  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputPrimaries = new TList();
  fOutputPrimaries->SetOwner();
  fOutput2Part = new TList();
  fOutput2Part->SetOwner();

  // Binning for lambdas in minv, y and pt
  const UChar_t nMinvBins = 180,
    nRapBins = 30;
  const Float_t minvLoEdge=1.070683, minvHiEdge=1.160683,
    rapLoEdge = -1.5, rapHiEdge = 1.5;
  // Non-uniformly wide pt bins
  const Int_t nPtBins=18;
  const Float_t ptBins[nPtBins+1]={0.,.25,.5,.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,4.,5.,6.,8.,10.,15.};
  // Do the rapBins
  Float_t rapBins[nRapBins+1];
  for(Int_t iRapBin=0;iRapBin<=nRapBins;iRapBin++){
    rapBins[iRapBin]=rapLoEdge + (rapHiEdge - rapLoEdge)/nRapBins * iRapBin;
  }
  // Do the minvBins
  Float_t minvBins[nMinvBins+1];
  for(Int_t iMinvBin=0;iMinvBin<=nMinvBins;iMinvBin++){
    minvBins[iMinvBin]=minvLoEdge + (minvHiEdge - minvLoEdge)/nMinvBins * iMinvBin;
  }
  

  //Try to keep memory low
  if(fkDoStdHists){
    // Control hist for event cuts
    fHistGoodEvent = new TH1F("h1GoodEvent","No of events passing the cuts.",10,-.5,9.5);
    fOutputList->Add(fHistGoodEvent);

    fHistCentComp = new TH2F ("CentComp"
			      ,Form("Comparision of cent estimators;%s;TRK",fkCentEst.Data())
			      ,101,0.,101
			      ,101,0.,101);
    fOutputList->Add(fHistCentComp);
    
    //
    // V0 offline distributons
    //
    
    // Only create them if we use offline V0s
    if(!fkUseOnTheFly){
      // Invariant mass distribution for the side band background
      fHistSideBandOffLam = new TH1F ("h1SideBandOffLam"
				      ,"m_{inv}(#Lambda) w/o any cuts;m_{inv}(#Lambda)"
				      ,nMinvBins,minvLoEdge,minvHiEdge);
      fOutputList->Add(fHistSideBandOffLam);
      fHistSideBandOffALam  = new TH1F ("h1SideBandOffALam"
					,"m_{inv}(#bar{#Lambda}) w/o any cuts;m_{inv}(#bar{#Lambda})"
					,nMinvBins,minvLoEdge,minvHiEdge);
      fOutputList->Add(fHistSideBandOffALam);
    
      // Invariant mass, invariant mass vs pt and y-pt
      fHistMassLambdaOff = new TH1F("h1MassLambdaOff"
				    , "#Lambda Offline candidates;M(p#pi^{-}) (GeV/c^{2});Counts"
				    ,nMinvBins,minvLoEdge,minvHiEdge);
      fOutputList->Add(fHistMassLambdaOff);
      fHistMassAntiLambdaOff = new TH1F("h1MassAntiLambdaOff"
					,"#bar{#Lambda} Offline candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts"
					,nMinvBins,minvLoEdge,minvHiEdge);
      fOutputList->Add(fHistMassAntiLambdaOff);
      
      // 3d y pt mass
      fHistYPtMassLamOff = new TH3F ("h3YPtMassLamOff","m_{inv}(#Lambda) vs y and pt;y;pt;mass"
				     ,nRapBins,rapBins,nPtBins,ptBins,nMinvBins,minvBins);
      fOutputList->Add(fHistYPtMassLamOff);
      fHistYPtMassALamOff = new TH3F ("h3YPtMassALamOff","m_{inv}(#bar{#Lambda}) vs y and pt;y;pt;mass"
				    ,nRapBins,rapBins,nPtBins,ptBins,nMinvBins,minvBins);
      fOutputList->Add(fHistYPtMassALamOff);
      
    }// End of only if we use offline V0s
    else{
    
      //
      // V0 on-the-fly distributons
      //
    
      // Invariant mass distribution for the side band background
      fHistSideBandOnLam = new TH1F ("h1SideBandOnLam","m_{inv}(#Lambda) w/o any cuts;m_{inv}(#Lambda)"
				     ,nMinvBins,minvLoEdge,minvHiEdge);
      fOutputList->Add(fHistSideBandOnLam);
      fHistSideBandOnALam  = new TH1F ("h1SideBandOnALam","m_{inv}(#bar{#Lambda}) w/o any cuts;m_{inv}(#bar{#Lambda})"
				       ,nMinvBins,minvLoEdge,minvHiEdge);
      fOutputList->Add(fHistSideBandOnALam);
    
      // Invariant mass, invariant mass vs pt and y-pt
      fHistMassLambdaOn  = new TH1F("h1MassLambdaOn"
				    ,"#Lambda on-the-fly candidates;M(p#pi^{-}) (GeV/c^{2});Counts"
				    ,nMinvBins,minvLoEdge,minvHiEdge);
      fOutputList->Add(fHistMassLambdaOn);
      fHistMassAntiLambdaOn = new TH1F("h1MassAntiLambdaOn"
				       ,"#bar{#Lambda} on-the-fly candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts"
				       ,nMinvBins,minvLoEdge,minvHiEdge);
      fOutputList->Add(fHistMassAntiLambdaOn);
    
      // 3d y pt mass
      fHistYPtMassLamOn = new TH3F ("h3YPtMassLamOn","m_{inv}(#Lambda) vs y and pt;y;pt;mass"
				    ,nRapBins,rapBins,nPtBins,ptBins,nMinvBins,minvBins);
      fOutputList->Add(fHistYPtMassLamOn);
      fHistYPtMassALamOn = new TH3F ("h3YPtMassALamOn","m_{inv}(#bar{#Lambda}) vs y and pt;y;pt;mass"
				     ,nRapBins,rapBins,nPtBins,ptBins,nMinvBins,minvBins);
      fOutputList->Add(fHistYPtMassALamOn);
    }

    // Do these histograms only if we need 'em
    //  if(!fdoLamOnly){
    //
    // Distributions for the primaries 
    //
    // Shared clusters
    fPriHistShare = new TH1F ("h1PriShare","Shared clusters, primaries;#shared clusters;counts",
			      160,0,160);
    fOutputPrimaries->Add(fPriHistShare);
    
    // TOF signal distribution when forcing TOF
    const UChar_t nBinsSigTOF(120),nBinsTofMom(28);
    const Float_t lowSigTOF(-5000.),hiSigTOF(1000.)
      ,lowTofMom(1.),hiTofMom(8.);
    fPriHistTOFsignalPosVsP = new TH2F ("h2TOFsignalPosVsP"
					,"tof signal vs p (positives);p [GeV/c];t_{meas} - t_{0} - t_{expected} [ps]"
					,nBinsTofMom,lowTofMom,hiTofMom,nBinsSigTOF,lowSigTOF,hiSigTOF);
    fOutputPrimaries->Add(fPriHistTOFsignalPosVsP);

    fPriHistTOFsignalNegVsP = new TH2F ("h2TOFsignalNegVsP"
					,"tof signal vs p (negatives);p [GeV/c];t_{meas} - t_{0} - t_{expected} [ps]"
					,nBinsTofMom,lowTofMom,hiTofMom,nBinsSigTOF,lowSigTOF,hiSigTOF);
    fOutputPrimaries->Add(fPriHistTOFsignalNegVsP);

    // Hybrid analysis
    fPriHistHybridTOFsigPosWoTPC = new TH1F ("h1HybridTOFsigPosWoTPC"
					     ,"tof signal pos (p=.75-1.0GeV) w/o dedx sel.;t_{meas} - t_{0} - t_{expected} [ps]"
					     ,100,-8000.0,2000.0);
    fOutputPrimaries->Add(fPriHistHybridTOFsigPosWoTPC);
    fPriHistHybridTOFsigPosTPCok = new TH1F ("h1HybridTOFsigPosTPCok"
					     ,"tof signal pos (p=.75-1.0GeV) with dedx sel.;t_{meas} - t_{0} - t_{expected} [ps]"
					     ,100,-8000.0,2000.0);
    fOutputPrimaries->Add(fPriHistHybridTOFsigPosTPCok);
    fPriHistHybridTOFsigNegWoTPC = new TH1F ("h1HybridTOFsigNegWoTPC"
					     ,"tof signal neg (p=.75-1.0GeV) w/o dedx sel.;t_{meas} - t_{0} - t_{expected} [ps]"
					     ,100,-8000.0,2000.0);
    fOutputPrimaries->Add(fPriHistHybridTOFsigNegWoTPC);
    fPriHistHybridTOFsigNegTPCok = new TH1F ("h1HybridTOFsigNegTPCok"
					     ,"tof signal neg (p=.75-1.0GeV) with dedx sel.;t_{meas} - t_{0} - t_{expected} [ps]"
					     ,100,-8000.0,2000.0);
    fOutputPrimaries->Add(fPriHistHybridTOFsigNegTPCok);
    // dEdx analysis
    fPriHistTPCsignalPos = new TH2F ("h2TPCsignalPos"
				     ,"TPC signal for positives;p_{tot};dEdx"
				     ,40,0,4,100,0,400);
    fOutputPrimaries->Add(fPriHistTPCsignalPos);

    fPriHistTPCsignalNeg = new TH2F ("h2TPCsignalNeg"
				     ,"TPC signal for negatives;p_{tot};dEdx"
				     ,40,0.0,4.0,100,0.0,400.0);
    fOutputPrimaries->Add(fPriHistTPCsignalNeg);
    
    fPriHistTPCsignalLowPPos = new TH2F ("h2TPCsignalLowPPos"
					 ,"dEdx for low momenta, positives"
					 ,20,0.1,0.3,300,0,3000);
    fOutputPrimaries->Add(fPriHistTPCsignalLowPPos);
    fPriHistTPCsignalMedPPos = new TH2F ("h2TPCsignalMedPPos"
					 ,"dEdx for medium momenta, positives"
					 ,60,0.3,0.9,500,0,500);
    fOutputPrimaries->Add(fPriHistTPCsignalMedPPos);
    fPriHistTPCsignalHigPPos = new TH2F ("h2TPCsignalHigPPos"
					 ,"dEdx for high momenta, positives"
					 ,100,0.9,1.9,120,0,120);
    fOutputPrimaries->Add(fPriHistTPCsignalHigPPos);
    fPriHistTPCsignalLowPNeg = new TH2F ("h2TPCsignalLowPNeg"
					 ,"dEdx for low momenta, negatives"
					 ,20,0.1,0.3,300,0,3000);
    fOutputPrimaries->Add(fPriHistTPCsignalLowPNeg);
    fPriHistTPCsignalMedPNeg = new TH2F ("h2TPCsignalMedPNeg"
					 ,"dEdx for medium momenta, negatives"
					 ,60,0.3,0.9,500,0,500);
    fOutputPrimaries->Add(fPriHistTPCsignalMedPNeg);
    fPriHistTPCsignalHigPNeg = new TH2F ("h2TPCsignalHigPNeg"
					 ,"dEdx for high momenta, negatives"
					 ,100,0.9,1.9,120,0,120);
    fOutputPrimaries->Add(fPriHistTPCsignalHigPNeg);
  
    //  Common for all protons

    // DCA xy and z distribution to determine primaries, secondaries from weak decay and secondaries from material
    // We create an array of hists. One for each phase space bin. Each hist is 2d in DCAxy vs DCAz. Since Sept 28
    // Only create them (heavy in memory) if requested
    if(fkDoDCAHists){
      Float_t DCAxyBins[47];
      UChar_t iEdge(0);
      for(Float_t dcaEdge=-2.4;dcaEdge<2.41;iEdge++){
	//printf("iEdge %u, dca %.2f\n",iEdge,dcaEdge);
	DCAxyBins[iEdge]=dcaEdge;
	if(dcaEdge<-1.1){dcaEdge+=.2;}
	else if(dcaEdge<-.35){dcaEdge+=.1;}
	else if(dcaEdge<-.125){dcaEdge+=.05;}
	else if(dcaEdge<.095){dcaEdge+=.02;}
	else if(dcaEdge<.275){dcaEdge+=.05;}
	else if(dcaEdge<.95){dcaEdge+=.1;}
	else dcaEdge+=.2;
      }
      Float_t DCAzBins[55];
      iEdge=0;
      for(Float_t dcaEdge=-3.2;dcaEdge<3.21;iEdge++){
	//printf("iEdge %u, dca %.2f\n",iEdge,dcaEdge);
	DCAzBins[iEdge]=dcaEdge;
	if(dcaEdge<-1.1){dcaEdge+=.2;}
	else if(dcaEdge<-.35){dcaEdge+=.1;}
	else if(dcaEdge<-.125){dcaEdge+=.05;}
	else if(dcaEdge<.095){dcaEdge+=.02;}
	else if(dcaEdge<.275){dcaEdge+=.05;}
	else if(dcaEdge<.95){dcaEdge+=.1;}
	else dcaEdge+=.2;
      }
      // We have 36 bins in (y,pt), 6 in each dim
      const Float_t DCAxyzRapBins[]={-1.,-.5,-.25,0.,.25,.5,1.};
      const Float_t DCAxyzPtBins[]={0.,.5,1.,1.5,2.,3.,5.};
      // Create the array of hists in phase space
      fPriHistDCAxyzYPtPro = new TH2F**[6];
      fPriHistDCAxyzYPtAPro = new TH2F**[6];
      for(UChar_t iY(0);iY<6;iY++){
	fPriHistDCAxyzYPtPro[iY] = new TH2F*[6];
	fPriHistDCAxyzYPtAPro[iY] = new TH2F*[6];
	for(UChar_t iPt(0);iPt<6;iPt++){
	  // Create hist for protons
	  fPriHistDCAxyzYPtPro[iY][iPt] = new TH2F (Form("h2DCAxyzY%uPt%uPro",iY,iPt)
						    ,Form("DCA_{xy} vs DCA_{z} p (y%.2f-%.2f,pt%.1f-%.1f)"
							  ,DCAxyzRapBins[iY],DCAxyzRapBins[iY+1]
							  ,DCAxyzPtBins[iPt],DCAxyzPtBins[iPt+1])
						    ,sizeof(DCAxyBins)/sizeof(Float_t)-1
						    ,DCAxyBins
						    ,sizeof(DCAzBins)/sizeof(Float_t)-1
						    ,DCAzBins);
	  fOutputPrimaries->Add(fPriHistDCAxyzYPtPro[iY][iPt]);
	  // .. and anti-protons
	  fPriHistDCAxyzYPtAPro[iY][iPt] = new TH2F (Form("h2DCAxyzY%uPt%uAPro",iY,iPt)
						     ,Form("DCA_{xy} vs DCA_{z} #bar{p} (y%.2f-%.2f,pt%.1f-%.1f)"
							   ,DCAxyzRapBins[iY],DCAxyzRapBins[iY+1]
							   ,DCAxyzPtBins[iPt],DCAxyzPtBins[iPt+1])
						     ,sizeof(DCAxyBins)/sizeof(Float_t)-1
						     ,DCAxyBins
						     ,sizeof(DCAzBins)/sizeof(Float_t)-1
						     ,DCAzBins);
	  fOutputPrimaries->Add(fPriHistDCAxyzYPtAPro[iY][iPt]);
	}// Loop over iPt
      }// Loop over iY
    }// Only create DCA histograms if requested

    //
    //  2 particle histograms fOutput2Part
    //

    // Binning for 2d ang dist at R=1.2m
    const UChar_t nPhiSBins(40),nEtaBins(60);
    const Float_t phiSMin(0.),phiSMax(.1),etaMin(-.3),etaMax(.3);
    // 2d Angular distance at R=1.2m with shifted vertices
    // ..real
    f2HistLamProAngDistSft2dAtR12Real = new TH2F("h2LamProAngDistSft2dAtR12Real"
						 ,"2d ang. dist.(R=1.2m) p(#Lambda)-p"
						 ";#Delta#eta;#Delta#phi*"
						 ,nEtaBins,etaMin,etaMax
						 ,nPhiSBins,phiSMin,phiSMax);
    fOutput2Part->Add(f2HistLamProAngDistSft2dAtR12Real);
    f2HistALamAProAngDistSft2dAtR12Real = new TH2F("h2ALamAProAngDistSft2dAtR12Real"
						   ,"2d ang. dist.(R=1.2m) #bar{p}(#bar{#Lambda})-#bar{p}"
						   ";#Delta#eta;#Delta#phi*"
						   ,nEtaBins,etaMin,etaMax
						   ,nPhiSBins,phiSMin,phiSMax);
    fOutput2Part->Add(f2HistALamAProAngDistSft2dAtR12Real);
    // Bg lambdas only if enabled
    if(fkDoBgLamALam){
      f2HistBgLamProAngDistSft2dAtR12Real = new TH2F("h2BgLamProAngDistSft2dAtR12Real"
						     ,"2d ang. dist.(R=1.2m) p(Bg#Lambda)-p"
						     ";#Delta#eta;#Delta#phi*"
						     ,nEtaBins,etaMin,etaMax
						     ,nPhiSBins,phiSMin,phiSMax);
      fOutput2Part->Add(f2HistBgLamProAngDistSft2dAtR12Real);
      f2HistBgALamAProAngDistSft2dAtR12Real = new TH2F("h2BgALamAProAngDistSft2dAtR12Real"
						       ,"2d ang. dist.(R=1.2m) #bar{p}(#bar{Bg#Lambda})-#bar{p}"
						       ";#Delta#eta;#Delta#phi*"
						       ,nEtaBins,etaMin,etaMax
						       ,nPhiSBins,phiSMin,phiSMax);
      fOutput2Part->Add(f2HistBgALamAProAngDistSft2dAtR12Real);
    }
    // ..mixed
    f2HistLamProAngDistSft2dAtR12Mixed = new TH2F("h2LamProAngDistSft2dAtR12Mixed"
						  ,"2d ang. dist.(R=1.2m) p(#Lambda)-p"
						  ";#Delta#eta;#Delta#phi*"
						  ,nEtaBins,etaMin,etaMax
						  ,nPhiSBins,phiSMin,phiSMax);
    fOutput2Part->Add(f2HistLamProAngDistSft2dAtR12Mixed);
    f2HistALamAProAngDistSft2dAtR12Mixed = new TH2F("h2ALamAProAngDistSft2dAtR12Mixed"
						    ,"2d ang. dist.(R=1.2m) #bar{p}(#bar{#Lambda})-#bar{p}"
						    ";#Delta#eta;#Delta#phi*"
						    ,nEtaBins,etaMin,etaMax
						    ,nPhiSBins,phiSMin,phiSMax);
    fOutput2Part->Add(f2HistALamAProAngDistSft2dAtR12Mixed);
    // Bg lambdas only if enabled
    if(fkDoBgLamALam){
      f2HistBgLamProAngDistSft2dAtR12Mixed = new TH2F("h2BgLamProAngDistSft2dAtR12Mixed"
						      ,"2d ang. dist.(R=1.2m) p(Bg#Lambda)-p"
						      ";#Delta#eta;#Delta#phi*"
						      ,nEtaBins,etaMin,etaMax
						      ,nPhiSBins,phiSMin,phiSMax);
      fOutput2Part->Add(f2HistBgLamProAngDistSft2dAtR12Mixed);
      f2HistBgALamAProAngDistSft2dAtR12Mixed = new TH2F("h2BgALamAProAngDistSft2dAtR12Mixed"
							,"2d ang. dist.(R=1.2m) #bar{p}(#bar{Bg#Lambda})-#bar{p}"
							";#Delta#eta;#Delta#phi*"
							,nEtaBins,etaMin,etaMax
							,nPhiSBins,phiSMin,phiSMax);
      fOutput2Part->Add(f2HistBgALamAProAngDistSft2dAtR12Mixed);
    }
    // Mt of the pairs
    //if(fDoMtHists) {
    const UChar_t nMtBins=25;
    const Float_t mtLow=1.0,mtHig=3.5;
    f2HistMtLamProReal = new TH1F("h1MtLamProReal"
				  ,"m_{t}(p #Lambda);m_{T} [GeV];counts"
				  ,nMtBins,mtLow,mtHig);
    fOutput2Part->Add(f2HistMtLamProReal);
    f2HistMtALamAProReal  = new TH1F("h1MtALamAProReal"
				     ,"m_{t}(#bar{p} #bar{#Lambda});m_{T} [GeV];counts"
				     ,nMtBins,mtLow,mtHig);
    fOutput2Part->Add(f2HistMtALamAProReal);
    // The same only filling for low q pairs
    f2HistMtLowQLamProReal = new TH1F("h1MtLowQLamProReal"
				      ,"m_{t}(p #Lambda);m_{T} [GeV];counts"
				      ,nMtBins,mtLow,mtHig);
    fOutput2Part->Add(f2HistMtLowQLamProReal);
    f2HistMtLowQALamAProReal  = new TH1F("h1MtLowQALamAProReal"
					 ,"m_{t}(#bar{p} #bar{#Lambda});m_{T} [GeV];counts"
					 ,nMtBins,mtLow,mtHig);
    fOutput2Part->Add(f2HistMtLowQALamAProReal);
    
    // Mt vs phi*
    f2HistMtVsPhiSLowQ = new TH2F("h2MtVsPhiSLowQ","m_{T} vs #phi* p#Lambda&#bar{p}#bar{#Lambda} q < 0.2"
				  ,50,0.,.5
				  ,50,1.,3.5);
    fOutput2Part->Add(f2HistMtVsPhiSLowQ);
    //}

    // Monitor needed buffer size
    f2HistNLamBefClean = new TH1F("h1NLamBefClean"
				  ,"Number of #Lambda before cleaning proc."
				  "; N(#Lambda)"
				  ,31,-.5,30.5);
    fOutput2Part->Add(f2HistNLamBefClean);
    f2HistNProBefClean = new TH1F("h1NProBefClean"
				  ,"Number of protons before cleaning proc."
				  "; N(p)"
				  ,31,-.5,30.5);
    fOutput2Part->Add(f2HistNProBefClean);
    f2HistNALamBefClean = new TH1F("h1NALamBefClean"
				   ,"Number of #bar{#Lambda} before cleaning proc."
				   "; N(#bar{#Lambda})"
				   ,31,-.5,30.5);
    fOutput2Part->Add(f2HistNALamBefClean);
    f2HistNAProBefClean = new TH1F("h1NAProBefClean"
				   ,"Number of anti-protons before cleaning proc."
				   "; N(#bar{p})"
				   ,31,-.5,30.5);
    fOutput2Part->Add(f2HistNAProBefClean);

  } // End of if(fkDoStdHists){ 
  
  
  // Common qinv binning
  const UChar_t nQinvBins = 50; // also for minv
  const Float_t QinvLow = 0.0;
  const Float_t QinvHig = .5;
  // Binning for the c2 in ang dist at r=1.2m
  const UChar_t nPhiSBinsC2(20),nEtaBinsC2(10);
  const Float_t phiSMinC2(0.),phiSMaxC2(.05),etaMinC2(0.),etaMaxC2(.2);
  // Mt binning of c2
  const UChar_t nMtBinsC2 = 15;
  const Float_t mtLowC2 = 1.;
  const Float_t mtHigC2 = 2.5;
  
  // Sept'12 Use a THn for (Bg)(A)Lam(A)Pro with 4 dimensions:
  // qinv (ppri,ppri), mean dist (ppri,ppri), min dist(ppri,ppri)
  // qinv (lam pro)
  const Int_t HnBins[4]={nMtBinsC2,nEtaBinsC2,nPhiSBinsC2,nQinvBins};
  const Double_t HnMin[4]={mtLowC2,etaMinC2,phiSMinC2,QinvLow};
  const Double_t HnMax[4]={mtHigC2,etaMaxC2,phiSMaxC2,QinvHig};
  fLamProReal = new THnF("Hn4LamProReal","lamProRealMtdEtaSdPhiSQlamp"
			 ,4,HnBins,HnMin,HnMax);
  fOutput2Part->Add(fLamProReal);
  fALamAProReal = new THnF("Hn4ALamAProReal","alamAProRealMtdEtaSdPhiSQlamp"
			   ,4,HnBins,HnMin,HnMax);
  fOutput2Part->Add(fALamAProReal);
  
  // Pair purity as fct of qinv,
  // really unneccessary if we have it as a fct of qinv and mT below
  if(fkDoStdHists){
    f2HistLamProWLamPurReal = new TH1F ("h1LamProWLamPurReal"
					,"p#Lambda pairs with #Lambda purity;q_{inv} (GeV/c);pairs"
					,nQinvBins,QinvLow,QinvHig);
    fOutput2Part->Add(f2HistLamProWLamPurReal);
    f2HistLamProWLamFdPurReal = new TH1F("h1LamProWLamFdPurReal"
					 ,"p#Lambda pairs with #Lambda FD purity;q_{inv} (GeV/c);pairs"
					 ,nQinvBins,QinvLow,QinvHig);
    fOutput2Part->Add(f2HistLamProWLamFdPurReal);
    f2HistLamProWProFdPurReal = new TH1F("h1LamProWProFdPurReal"
					 ,"p#Lambda pairs with p FD purity;q_{inv} (GeV/c);pairs"
					 ,nQinvBins,QinvLow,QinvHig);
    fOutput2Part->Add(f2HistLamProWProFdPurReal);
    f2HistLamProWoPurReal = new TH1F ("h1LamProWoPurReal"
				      ,"p#Lambda pairs without purity;q_{inv} (GeV/c);pairs"
				      ,nQinvBins,QinvLow,QinvHig);
    fOutput2Part->Add(f2HistLamProWoPurReal);
    f2HistALamAProWALamPurReal = new TH1F ("h1ALamAProWALamPurReal"
					   ,"#bar{p}#bar{#Lambda} pairs with #bar{#Lambda} purity"
					   ";q_{inv} (GeV/c);pairs"
					   ,nQinvBins,QinvLow,QinvHig);
    fOutput2Part->Add(f2HistALamAProWALamPurReal);
    f2HistALamAProWALamFdPurReal = new TH1F("h1ALamAProWALamFdPurReal"
					    ,"#bar{p}#bar{#Lambda} pairs with #bar{#Lambda} FD purity"
					    ";q_{inv} (GeV/c);pairs"
					    ,nQinvBins,QinvLow,QinvHig);
    fOutput2Part->Add(f2HistALamAProWALamFdPurReal);
    f2HistALamAProWAProFdPurReal = new TH1F("h1ALamAProWAProFdPurReal"
					    ,"#bar{p}#bar{#Lambda} pairs with #bar{p} FD purity"
					    ";q_{inv} (GeV/c);pairs"
					    ,nQinvBins,QinvLow,QinvHig);
    fOutput2Part->Add(f2HistALamAProWAProFdPurReal);
    f2HistALamAProWoPurReal = new TH1F ("h1ALamAProWoPurReal"
					,"#bar{p}#bar{#Lambda} pairs with #bar{p} FD purity"
					";q_{inv} (GeV/c);pairs"
					,nQinvBins,QinvLow,QinvHig);
    fOutput2Part->Add(f2HistALamAProWoPurReal);
  }// End of if(fkDoStdHists){
  
  // Pair purity as fct of qinv and mT
  f2HistLamProWLamPurmTReal = new TH2F ("h2LamProWLamPurmTReal"
					,"p#Lambda pairs with #Lambda purity;q_{inv} (GeV/c);m_{T} (GeV/c);pairs"
					,nQinvBins,QinvLow,QinvHig,nMtBinsC2,mtLowC2,mtHigC2);
  fOutput2Part->Add(f2HistLamProWLamPurmTReal);
  f2HistLamProWLamFdPurmTReal = new TH2F("h2LamProWLamFdPurmTReal"
					 ,"p#Lambda pairs with #Lambda FD purity;q_{inv} (GeV/c);m_{T} (GeV/c);pairs"
					 ,nQinvBins,QinvLow,QinvHig,nMtBinsC2,mtLowC2,mtHigC2);
  fOutput2Part->Add(f2HistLamProWLamFdPurmTReal);
  f2HistLamProWProFdPurmTReal = new TH2F("h2LamProWProFdPurmTReal"
					 ,"p#Lambda pairs with p FD purity;q_{inv} (GeV/c);m_{T} (GeV/c);pairs"
					 ,nQinvBins,QinvLow,QinvHig,nMtBinsC2,mtLowC2,mtHigC2);
  fOutput2Part->Add(f2HistLamProWProFdPurmTReal);
  f2HistLamProWoPurmTReal = new TH2F ("h2LamProWoPurmTReal"
				      ,"p#Lambda pairs without purity;q_{inv} (GeV/c);m_{T} (GeV/c);pairs"
				      ,nQinvBins,QinvLow,QinvHig,nMtBinsC2,mtLowC2,mtHigC2);
  fOutput2Part->Add(f2HistLamProWoPurmTReal);
  f2HistALamAProWALamPurmTReal = new TH2F ("h2ALamAProWALamPurmTReal"
					   ,"#bar{p}#bar{#Lambda} pairs with #bar{#Lambda} purity"
					   ";q_{inv} (GeV/c);m_{T} (GeV/c);pairs"
					   ,nQinvBins,QinvLow,QinvHig,nMtBinsC2,mtLowC2,mtHigC2);
  fOutput2Part->Add(f2HistALamAProWALamPurmTReal);
  f2HistALamAProWALamFdPurmTReal = new TH2F("h2ALamAProWALamFdPurmTReal"
					    ,"#bar{p}#bar{#Lambda} pairs with #bar{#Lambda} FD purity"
					    ";q_{inv} (GeV/c);m_{T} (GeV/c);pairs"
					    ,nQinvBins,QinvLow,QinvHig,nMtBinsC2,mtLowC2,mtHigC2);
  fOutput2Part->Add(f2HistALamAProWALamFdPurmTReal);
  f2HistALamAProWAProFdPurmTReal = new TH2F("h2ALamAProWAProFdPurmTReal"
					    ,"#bar{p}#bar{#Lambda} pairs with #bar{p} FD purity"
					    ";q_{inv} (GeV/c);m_{T} (GeV/c);pairs"
					    ,nQinvBins,QinvLow,QinvHig,nMtBinsC2,mtLowC2,mtHigC2);
  fOutput2Part->Add(f2HistALamAProWAProFdPurmTReal);
  f2HistALamAProWoPurmTReal = new TH2F ("h2ALamAProWoPurmTReal"
					,"#bar{p}#bar{#Lambda} pairs with #bar{p} FD purity"
					";q_{inv} (GeV/c);m_{T} (GeV/c);pairs"
					,nQinvBins,QinvLow,QinvHig,nMtBinsC2,mtLowC2,mtHigC2);
  fOutput2Part->Add(f2HistALamAProWoPurmTReal);
  
  // Fluctuation of the purity
  f2HistPurFlucLamPro = new TH1F("f2HistPurFlucLamPro","f2HistPurFlucLamPro qinv={.0-.2};pur;counts",100,.7,1.);
  fOutput2Part->Add(f2HistPurFlucLamPro);
  f2HistPurFlucALamAPro= new TH1F("f2HistPurFlucALamAPro","f2HistPurFlucALamAPro qinv={.0-.2};pur;counts",100,.7,1.);
  fOutput2Part->Add(f2HistPurFlucALamAPro);

  if(fkDoStdHists){
    //
    // Origin in phase space of the pairs..
    // 
    // ..lambda (y,pt) for each pair
    f2HistLamYPtPair = new TH2F("f2HistLamYPtPair","#Lambda y,pt for each pair w/ q_{inv}< .2GeV/c;y;p_{T} (GeV/c)"
				//,nRapBins,rapLoEdge,rapHiEdge,nPtBins,ptLoEdge,ptHiEdge);
				,nRapBins,rapBins,nPtBins,ptBins);
    fOutput2Part->Add(f2HistLamYPtPair);
    f2HistALamYPtPair = new TH2F("f2HistALamYPtPair","#bar{#Lambda} y,pt for each pair w/ q_{inv}< .2GeV/c;y;p_{T} (GeV/c)"
				 //,nRapBins,rapLoEdge,rapHiEdge,nPtBins,ptLoEdge,ptHiEdge);
				 ,nRapBins,rapBins,nPtBins,ptBins);
    fOutput2Part->Add(f2HistALamYPtPair);
    // ..proton (y,pt) for each pair
    f2HistProYPtPair = new TH2F("f2HistProYPtPair","proton y,pt for each pair w/ q_{inv}< .2GeV/c;y;p_{T} (GeV/c)"
				//,nRapBins,rapLoEdge,rapHiEdge,nPtBins,ptLoEdge,ptHiEdge);
				,nRapBins,rapBins,nPtBins,ptBins);
    fOutput2Part->Add(f2HistProYPtPair);
    f2HistAProYPtPair = new TH2F("f2HistAProYPtPair","#bar{proton} y,pt for each pair w/ q_{inv}< .2GeV/c;y;p_{T} (GeV/c)"
				 //,nRapBins,rapLoEdge,rapHiEdge,nPtBins,ptLoEdge,ptHiEdge);
				 ,nRapBins,rapBins,nPtBins,ptBins);
    fOutput2Part->Add(f2HistAProYPtPair);
    // ..pair (y,kt) for each pair, low qinv
    f2HistYKtPair = new TH2F("f2HistYKtPair","y_{p#Lambda},kt for each pair w/ q_{inv}< .2GeV/c;y_{p#Lambda};k_{T} (GeV/c)"
			     //,nRapBins,rapLoEdge,rapHiEdge,nPtBins,ptLoEdge,ptHiEdge);
			     ,nRapBins,rapBins,nPtBins,ptBins);
    fOutput2Part->Add(f2HistYKtPair);
    f2HistYKtAPair = new TH2F("f2HistYKtAPair","y_{#bar{p#Lambda}},kt for each anti particle pair w/ q_{inv}< .2GeV/c;y_{#bar{p#Lambda}};k_{T} (GeV/c)"
			      //,nRapBins,rapLoEdge,rapHiEdge,nPtBins,ptLoEdge,ptHiEdge);
			      ,nRapBins,rapBins,nPtBins,ptBins);
    fOutput2Part->Add(f2HistYKtAPair);

    // ..lamdba (y,pt) for each pair w/ high qinv; 
    // checks the origin of the rising purity with higher qinv
    f2HistLamYPtHiPair = new TH2F("f2HistLamYPtHiPair"
				  ,"y_{#Lambda},pt for each pair w/ .4 #leq q_{inv}< .5GeV/c"
				  ";y_{#Lambda};p_{T} (GeV/c)"
				  ,nRapBins,rapBins,nPtBins,ptBins);
    fOutput2Part->Add(f2HistLamYPtHiPair);
    f2HistALamYPtHiPair  = new TH2F("f2HistALamYPtHiPair"
				    ,"y_{#bar{#Lambda}},pt for each pair w/ .4 #leq q_{inv}< .5GeV/c"
				    ";y_{#bar{#Lambda}};p_{T} (GeV/c)"
				    ,nRapBins,rapBins,nPtBins,ptBins);
    fOutput2Part->Add(f2HistALamYPtHiPair);
  }

  // 
  // Mixed events
  //

  // Sept'12 Use a THn for (Bg)(A)Lam(A)Pro with 4 dimensions:
  // qinv (ppri,ppri), mean dist (ppri,ppri), min dist(ppri,ppri)
  // qinv (lam pro)
  fLamProMixed = new THnF("Hn4LamProMixed","lamProMixedQppdEtaSdPhiSQlamp"
			  ,4,HnBins,HnMin,HnMax);
  fOutput2Part->Add(fLamProMixed);
  fALamAProMixed = new THnF("Hn4ALamAProMixed","alamAProMixedQppdEtaSdPhiSQlamp"
			    ,4,HnBins,HnMin,HnMax);
  fOutput2Part->Add(fALamAProMixed);

  // Sept'12 Use a THn for (Bg)(A)Lam(A)Pro with 4 dimensions:
  // qinv (ppri,ppri), mean dist (ppri,ppri), min dist(ppri,ppri)
  // qinv (lam pro)
        
  // Bg lambdas only if enabled
  if(fkDoBgLamALam){
    fBgLamProReal = new THnF("Hn4BgLamProReal","lamProRealQppdEtaSdPhiSQlamp"
			     ,4,HnBins,HnMin,HnMax);
    fOutput2Part->Add(fBgLamProReal);
    fBgALamAProReal = new THnF("Hn4BgALamAProReal","alamAProRealQppdEtaSdPhiSQlamp"
			       ,4,HnBins,HnMin,HnMax);
    fOutput2Part->Add(fBgALamAProReal);
    
    // Sept'12 Use a THn for (Bg)(A)Lam(A)Pro with 4 dimensions:
    // qinv (ppri,ppri), mean dist (ppri,ppri), min dist(ppri,ppri)
    // qinv (lam pro)
    fBgLamProMixed = new THnF("Hn4BgLamProMixed","lamProMixedQppdEtaSdPhiSQlamp"
			      ,4,HnBins,HnMin,HnMax);
    fOutput2Part->Add(fBgLamProMixed);
    fBgALamAProMixed = new THnF("Hn4BgALamAProMixed","alamAProMixedQppdEtaSdPhiSQlamp"
				,4,HnBins,HnMin,HnMax);
    fOutput2Part->Add(fBgALamAProMixed);
  }
  
    //  } // End of only do the hists if we need 'em

  // Post the data
  PostData(1, fOutputList);
  PostData(2, fOutputPrimaries);
  PostData(3, fOutput2Part);

}

//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Fill a control histogram
  if(fkDoStdHists)
    fHistGoodEvent->Fill(0.0);

  // Get the event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;
  }

  // Fill a control histogram
  if(fkDoStdHists)
    fHistGoodEvent->Fill(1.0);  

  // Get the centrality selection
  const AliCentrality *centrality = fAOD->GetCentrality();
  if (!centrality) {
    printf ("ERROR: couldn't get the AliCentrality\n");
    return;
  }
  
  // Fill a control histogram
  if(fkDoStdHists)
    fHistGoodEvent->Fill(2.0);  

  // Check the fQuality flag of the centrality task
  // for details see
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies#How_we_determine_centrality
  if (centrality->GetQuality()){
    return;
  }

  // Fill a control histogram
  if(fkDoStdHists)
    fHistGoodEvent->Fill(3.0);  

  // Analyze only 20% most central events using multiplicity in V0 detector (standard)
  const Float_t centPerc = centrality->GetCentralityPercentileUnchecked(fkCentEst.Data());
  if ( (centPerc < fkCentCutLo) ||
       (centPerc >= fkCentCutHi) ) {
    return;
  }

  // Fill a control histogram
  if(fkDoStdHists)
    fHistGoodEvent->Fill(4.0);  

  // Primary vertex, GetPrimaryVertex() returns the "best" reconstructed vertex
  fPrimaryVtx = fAOD->GetPrimaryVertex();
  if (!fPrimaryVtx){
    printf ("ERROR: no primary vertex\n");
    return;
  }

  // Fill a control histogram
  if(fkDoStdHists)
    fHistGoodEvent->Fill(5.0);  
  fPrimaryVtx->GetXYZ(fPrimaryVtxPosition);
  // fHistPrimaryVertexPosXY->Fill(fPrimaryVtxPosition[0],fPrimaryVtxPosition[1]);
  // fHistPrimaryVertexPosZ->Fill(fPrimaryVtxPosition[2]);
  
  // Zvertex cut, probably done anyhow in centrality task
  if (TMath::Abs(fPrimaryVtxPosition[2]) > fkAbsZvertexCut)
    return;
  
  // Fill a control histogram
  if(fkDoStdHists)
    fHistGoodEvent->Fill(6.0);

  // Multiplicity
  if (!(fAOD->GetNumberOfTracks())) {
    return;
  }
  
  // Fill a control histogram
  if(fkDoStdHists)
    fHistGoodEvent->Fill(7.0);
  
  // Simple memory profiling
#ifdef ROOT_TObjectTable
  static ULong_t nEvt=0;
  if(nEvt==1||nEvt==10){
    printf(" Object table after %lu events\n",nEvt);
    printf("------------------------------\n");
    gObjectTable->Print();
    printf("------------------------------\n");
  }
  nEvt++;
#endif

  // Compare diff cent estimators
  if(fkDoStdHists){
    const Float_t centTRK = centrality->GetCentralityPercentileUnchecked("TRK");
    fHistCentComp->Fill(centPerc,centTRK);
  }

  // fHistTrackMultiplicity->Fill(fAOD->GetNumberOfTracks());

  // Set up the event buffer to store this event
  //  if(!fdoLamOnly) //only if we need it
  fFemtoBuffer->ShiftAndAdd(fAOD,fkCentEst);

  // // Debugging: print number of stored tracks in the event
  // for(UChar_t i=0;i<fFemtoBuffer->GetMixBuffSize();i++)
  //   printf("iMix: %u, NPro %u, NAPro %u, NLam %u, NALam %u"
  // 	   "NBgLam %u, NBgALam %u\n"
  // 	   ,i
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNPro()
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNAPro()
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNLam()
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNALam()
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNBgLam()
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNBgALam()
  // 	   );
  // printf("\n");

  // Reset the reference array to the global tracks..
  ResetGlobalTrackReference();
  // ..and set it
  for (Int_t iTrack=0;iTrack<fAOD->GetNumberOfTracks();iTrack++){
    // cast needed since the event now returns AliVTrack instead of AliAODTrack
    const AliAODTrack *track = static_cast<const AliAODTrack *>(fAOD->GetTrack(iTrack));
    if (!track) continue;
    
    // Store the reference of the global tracks
    StoreGlobalTrackReference(track);
  }
  
  // V0 loop
  for (Int_t iV0 = 0; iV0 < fAOD->GetNumberOfV0s(); iV0++) {
    const AliAODv0 *v0 = fAOD->GetV0(iV0);

    // Skip if V0 is not there
    if((!v0))
      continue;

    // Check that the array fGTI isn't too small
    // for the track ids
    if(v0->GetPosID() >= fTrackBuffSize||
       v0->GetNegID() >= fTrackBuffSize)
      continue;

    // This is AODs: find the track for given id:
    const AliAODTrack *pTrack=fGTI[v0->GetPosID()];
    const AliAODTrack *nTrack=fGTI[v0->GetNegID()];
	
    // Skip if one of the daughter is not there
    if ((!pTrack) || (!nTrack)) continue;

    // Famous crossed rows / findable clusters cut,
    // rejects split tracks very well
    // (Don't do it for the V0s as we require 80 clusters 
    // and reject shared clusters)
    //    if( (!acceptTrack(pTrack)) || (!acceptTrack(nTrack)) )
    //      continue;

    // Reject tracks with shared clusters
    if(!GoodTPCFitMapSharedMap(pTrack,nTrack))
      continue;

    // Analysis done seperately for offline and on-the-fly
    if (!(v0->GetOnFlyStatus())){
      // Only do the offline V0s if we use them
      if(!fkUseOnTheFly){
	ProcessOffline(v0, pTrack, nTrack);
      }
    }
    else
      ProcessOnTheFly(v0, pTrack, nTrack);

    // V0s get added to the mixed events in the 'Process..' fcts
    
  } // End of V0 loop

  // Stop here if we only want the lambdas
  // if(fdoLamOnly)
  //   return;

  // Loop over primary tracks
  for (Int_t iTrack=0;iTrack<fAOD->GetNumberOfTracks();iTrack++){
        // cast needed since the event now returns AliVTrack instead of AliAODTrack
    const AliAODTrack *track = static_cast<const AliAODTrack *>(fAOD->GetTrack(iTrack));
    if (!track) continue;
    
    if(!track->TestFilterBit(128))
      continue;

    // Famous crossed rows / findable clusters cut,
    // rejects split tracks very well
    if(!acceptTrack(track))
      continue;

    // Reject tracks with shared clusters
    if(!GoodTPCFitMapSharedMap(track))
      continue;

    // Check that the array fGTI isn't too small
    // for the track id
    if(-track->GetID()-1 >= fTrackBuffSize)
      continue;

    // Without a corresponding global track it's useless
    if(!fGTI[-track->GetID()-1]){
      printf ("No global info! iTrack %d, ID %d\n",iTrack,track->GetID());
      continue;
    }

    // Visualization of TPC dE/dx
    FillDedxHist(track);

    // Depending on momentum choose pid method
    if (track->P() < 0.75){
       ProcessTPC(track);
    }
    else if (track->P() < 1.0){
       ProcessHybrid(track);
    }
    else {
      // TOF spectrum will get filled for tracks up to 8 GeV (histogram bounds),
      // used are still only tracks up to p_tot = 5.0 GeV/c
      ProcessTOF(track);
    }

    
    // Tracks get added to the mixed events in the 'Process..' fcts

  } // End of loop over primary tracks

  // Track cuts do not allow for split tracks

  // Cleaning procedure for lambdas & lambdas, lambdas & protons,
  // anti-lambdas & anti-lambdas, anti-lambdas & protons + (anti-)lambda background
  CleaningProcedure();

  // Process real events
  ProcessReal();
  if(fkDoBgLamALam)
    ProcessRealBackground();
  
  // Process mixed events
  ProcessMixed();
  if(fkDoBgLamALam)
    ProcessMixedBackground();

  // Post output data.
  PostData(1, fOutputList);
  PostData(2, fOutputPrimaries);
  PostData(3, fOutput2Part);

}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::ProcessOffline(const AliAODv0 *v0,const AliAODTrack *pTrack,const AliAODTrack *nTrack) 
{

  // For clarity in code: Fill some hists with on-the-fly status
  //  const Float_t kOnTheFlyStat = 0.0;

  // All cuts are checked with invariant mass histograms
  //  v0->ChangeMassHypothesis(3122);
  Float_t minvLam = v0->MassLambda();
  //  v0->ChangeMassHypothesis(-3122);
  Float_t minvALam = v0->MassAntiLambda();
  // Cosine as local variable as this is some computation
  const Float_t lCosPoint = v0->CosPointingAngle(fPrimaryVtxPosition);

  // Also calculate a V0 momentum with TPC only daughters
  //  Double_t TPConlyV0Mom[3], TPConlyV0MinvLam=0, TPConlyV0MinvALam=0;
  //  getTPConlyV0Info(pTrack, nTrack,
  //		   TPConlyV0Mom, TPConlyV0MinvLam, TPConlyV0MinvALam);

  // Fill a minv hist w/o any cuts. Select background from the sideband
  if(fkDoStdHists){
    fHistSideBandOffLam->Fill(minvLam);
    fHistSideBandOffALam->Fill(minvALam);
  }
  // Fill the event buffer w/ background if enabled
  if(fkDoBgLamALam){
    if (!fkUseOnTheFly){
      if ( TMath::Abs(minvLam - fkLamMass) > 0.015 &&
	   TMath::Abs(minvLam - fkLamMass) < 0.035 ){
	fFemtoBuffer->GetEvt(0)->AddBgLam(v0, pTrack, nTrack);
      }
      if ( TMath::Abs(minvALam - fkLamMass) > 0.015 &&
	   TMath::Abs(minvALam - fkLamMass) < 0.035 ){
	fFemtoBuffer->GetEvt(0)->AddBgALam(v0, pTrack, nTrack);
      }
    }
  }

  // Control histogram: fill all v0s
  // fHistGoodV0->Fill(0.0,kOnTheFlyStat);
  // fHistGoodV0->Fill(1.0,kOnTheFlyStat);

  // Require 80 TPC clusters for both pos and neg daughter
  // fHistTPCNclsPosOffLam->Fill(pTrack->GetTPCNcls(),minvLam);
  // fHistTPCNclsNegOffLam->Fill(nTrack->GetTPCNcls(),minvLam);
  // fHistTPCNclsPosOffALam->Fill(pTrack->GetTPCNcls(),minvALam);
  // fHistTPCNclsNegOffALam->Fill(nTrack->GetTPCNcls(),minvALam);

  if ( ( (pTrack->GetTPCNcls()) < 80 ) || ( (nTrack->GetTPCNcls()) < 80 ) ) 
    return;
  //  fHistGoodV0->Fill(2.0,kOnTheFlyStat);

  // Require a maximum dca of the daughters of 0.6cm
  // fHistDcaV0DaughtersOffLam->Fill(v0->DcaV0Daughters(),minvLam);
  // fHistDcaV0DaughtersOffALam->Fill(v0->DcaV0Daughters(),minvALam);
  // fHistDcaV0Daughters->Fill(v0->DcaV0Daughters(),kOnTheFlyStat);
  if (v0->DcaV0Daughters() > 0.6)
    return;
  //  fHistGoodV0->Fill(3.0,kOnTheFlyStat);
  
  // Force TPC PID to be present
  if (!(pTrack->GetStatus() & AliAODTrack::kTPCpid) ||
      !(nTrack->GetStatus() & AliAODTrack::kTPCpid))
    return;
  //  fHistGoodV0->Fill(4.0,kOnTheFlyStat);

  // Perform cut on TPC dE/dx
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton)) > 3.4)
    minvLam=0.0;
  // else 
  //   fHistPosTpcAfterCut->Fill(pTrack->P(),pTrack->GetTPCsignal());
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)) > 4.4)
    minvLam=0.0;
  // else
  //   fHistNegTpcAfterCut->Fill(nTrack->P(),nTrack->GetTPCsignal());
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)) > 4.2)
    minvALam=0.0;
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton)) > 3.4)
    minvALam=0.0;

  // Don't use a tof cut for pions

  // Check whether to use a 5sigma tof cut or none for protons
  // if (pTrack->GetStatus() & AliAODTrack::kTOFpid){
  //   if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kProton)) > 5.0)
  //     fHistUseTofOffLam->Fill(1.0,minvLam);
  //   else
  //     fHistUseTofOffLam->Fill(0.0,minvLam);
  // }
  // else
  //   fHistUseTofOffLam->Fill(0.0,minvLam);
  // Check whether to use a 5sigma tof cut or none for anti-protons
  // if (nTrack->GetStatus() & AliAODTrack::kTOFpid){
  //   if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kProton)) > 5.0)
  //     fHistUseTofOffALam->Fill(1.0,minvALam);
  //   else
  //     fHistUseTofOffALam->Fill(0.0,minvALam);
  // }
  // else
  //   fHistUseTofOffALam->Fill(0.0,minvALam);

  // Don't use a TOF cut for offline
  
  // Don't need to check for sign of pairs as this is always
  // correct for offline finder

  // Don't need to check for TPC refit as it is required
  // by the offline finder itself

  //
  // Require a minimum distance between daughters and primary vertex
  //
  // Fill histograms with the distributions before cutting
  // fHistDcaPosOffLam->Fill(v0->DcaPosToPrimVertex(),minvLam);
  // fHistDcaPosOffALam->Fill(v0->DcaPosToPrimVertex(),minvALam);
  // fHistDcaNegOffLam->Fill(v0->DcaNegToPrimVertex(),minvLam);
  // fHistDcaNegOffALam->Fill(v0->DcaNegToPrimVertex(),minvALam);
  
  // fHistDcaPosToPrimVertex->Fill(v0->DcaPosToPrimVertex(),kOnTheFlyStat);
  // fHistDcaNegToPrimVertex->Fill(v0->DcaNegToPrimVertex(),kOnTheFlyStat);
  // fHistDcaPosToPrimVertexZoom->Fill(v0->DcaPosToPrimVertex(),kOnTheFlyStat);
  // fHistDcaNegToPrimVertexZoom->Fill(v0->DcaNegToPrimVertex(),kOnTheFlyStat);
  
  // Do the cut
  if (v0->DcaPosToPrimVertex() < 0.1)
    minvLam=0.0;
  if (v0->DcaPosToPrimVertex() < 0.3)
    minvALam=0.0;
  if (v0->DcaNegToPrimVertex() < 0.1)
    minvALam=0.0;
  if (v0->DcaNegToPrimVertex() < 0.3)
    minvLam=0.0;

  // Cosine of pointing angle. Computed at the beginning.
  // Fill historgrams before cutting
  // fHistCosPointLamOff->Fill(lCosPoint,minvLam);
  // fHistCosPointALamOff->Fill(lCosPoint,minvALam);
  // fHistCosPointLamZoomOff->Fill(lCosPoint,minvLam);
  // fHistCosPointALamZoomOff->Fill(lCosPoint,minvALam);
  
  // fHistCosPointAngle->Fill(lCosPoint,kOnTheFlyStat);
  // fHistCosPointAngleZoom->Fill(lCosPoint,kOnTheFlyStat);

  // Do the cut in cos (pointing angle) 
  // (note the difference 0.9996 for offline and 0.9999 for on-the-fly)
  if (lCosPoint < 0.9996)
    return;
  
  // fHistGoodV0->Fill(7.0,kOnTheFlyStat);
  
  // Fill some histograms with cut variables
  // fHistChi2->Fill(v0->Chi2V0(),kOnTheFlyStat);
  
  // Idea to cut on the radius
  // fHistRadiusV0->Fill(v0->RadiusV0(),kOnTheFlyStat);
  // fHistV0RadiusLamOff->Fill(v0->RadiusV0(),minvLam);
  // fHistV0RadiusALamOff->Fill(v0->RadiusV0(),minvALam);

  // Idea to cut on the decay length
  // fHistDecayLengthV0->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),kOnTheFlyStat);
  // fHistV0DecayLengthLamOff->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),minvLam);
  // fHistV0DecayLengthALamOff->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),minvALam);
 
  // Idea to cut on DCA of V0 and primay vertex
  // fHistDcaV0PriVertexLamOff->Fill(v0->DcaV0ToPrimVertex(),minvLam);
  // fHistDcaV0PriVertexALamOff->Fill(v0->DcaV0ToPrimVertex(),minvALam);
     
  // Fill some invariant mass distributions
  if(fkDoStdHists){
    fHistMassLambdaOff->Fill(minvLam);
    fHistMassAntiLambdaOff->Fill(minvALam);
  }
  // fHistPtVsMassLambdaOff->Fill(v0->Pt(),minvLam);
  // fHistPtVsMassAntiLambdaOff->Fill(v0->Pt(),minvALam);

  // 3d histogram: rapidity, pt and mass
  if(fkDoStdHists){
    fHistYPtMassLamOff->Fill(v0->Y(3122),v0->Pt(),minvLam);
    fHistYPtMassALamOff->Fill(v0->Y(-3122),v0->Pt(),minvALam);
  }
  // Invariant mass cut lambda :: fill a y-pt hist
  // if ( TMath::Abs(minvLam - fkLamMass) < 0.01 ){
  //   fHistPtVsYLambdaOff->Fill(v0->Pt(),v0->Y(3122));
  // }
  // // Invariant mass cut anti-lambda :: fill a y-pt hist
  // if ( TMath::Abs(minvALam - fkLamMass) < 0.01 ){
  //   fHistPtVsYAntiLambdaOff->Fill(v0->Pt(),v0->Y(-3122));
  // }

  // Fill the mixed events when offline V0 finder is used,
  // and we want to do more than just getting the lambdas
  //  if (!fkUseOnTheFly && !fdoLamOnly){
  if (!fkUseOnTheFly){
    // Highest significance for minv +/- 4 MeV (std value for cut)
    if ( TMath::Abs(minvLam - fkLamMass) < fkMinvCut ){
      
      // We introduce a pt cut - only on lambdas! - of 0.5 GeV/c.
      // The reason for this is that below 0.5 GeV/c we have a 
      // huge background which will only ruin our whole lambda
      // sample. There's also the highest track density, so ex-
      // cluding this region will in addition probably help us 
      // fighting two-track effects. The signal is tiny, so we
      // don't loose much. The reason for the background is pro-
      // bably the protons from material?!
      if(v0->Pt() >= 0.5){
	fFemtoBuffer->GetEvt(0)->AddLam(v0, pTrack, nTrack);
      }
    }
    if ( TMath::Abs(minvALam - fkLamMass) < fkMinvCut ){
      fFemtoBuffer->GetEvt(0)->AddALam(v0, pTrack, nTrack);
    }
  }
}   
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::ProcessOnTheFly(const AliAODv0 *v0,const AliAODTrack *pTrack,const AliAODTrack *nTrack) 
{
  // For clarity in code: Fill some hists with on-the-fly status
  //  const Float_t kOnTheFlyStat = 1.0;

  // All cuts are checked with invariant mass histograms
  Float_t minvLam = v0->MassLambda();
  Float_t minvALam = v0->MassAntiLambda();
  const Float_t lCosPoint = v0->CosPointingAngle(fPrimaryVtxPosition);

  // Control histogram: fill all v0s
  //  fHistGoodV0->Fill(0.0,kOnTheFlyStat);
  // Control hist: after require two daughter tracks
  //  fHistGoodV0->Fill(1.0,kOnTheFlyStat);
  // Check the right sign of the tracks (mainly on-the-fly)
  if (pTrack->Charge() > 0 && nTrack->Charge() < 0){
    // Correct assignment
    // fHistCorrectSigns->Fill(0.0,kOnTheFlyStat);

    // fHistLikeSignOnLam->Fill(0.0,minvLam);
    // fHistLikeSignOnALam->Fill(0.0,minvALam);    
  }
  else if (pTrack->Charge() < 0 && nTrack->Charge() > 0){
    // Swapped sings
    //    fHistCorrectSigns->Fill(1.0,kOnTheFlyStat);

    pTrack = fGTI[v0->GetNegID()];
    nTrack = fGTI[v0->GetPosID()];
    
    // See http://savannah.cern.ch/bugs/?90749
    // For AODs it depends on with which root version 
    // the AODs got produced.

    // See above: swapping mass assignment
    minvLam = v0->MassAntiLambda();
    minvALam = v0->MassLambda();

    //    fHistLikeSignOnLam->Fill(1.0,minvLam);
    //    fHistLikeSignOnALam->Fill(1.0,minvALam);    
  }
  else {
    // Like sign pairs
    //    fHistCorrectSigns->Fill(2.0,kOnTheFlyStat);
    
    //    fHistLikeSignOnLam->Fill(2.0,minvLam);
    //    fHistLikeSignOnALam->Fill(2.0,minvALam);    

    // Don't use like sign-pairs
    return;
  }
  //  fHistGoodV0->Fill(2.0,kOnTheFlyStat);

  // V0 momentum
  Double_t V0Mom[3];
  v0->PxPyPz(V0Mom);
  // Also calculate a V0 momentum with TPC only daughters
  //  Double_t TPConlyV0Mom[3], TPConlyV0MinvLam=0, TPConlyV0MinvALam=0;
  //  getTPConlyV0Info(pTrack, nTrack,
  //		   TPConlyV0Mom, TPConlyV0MinvLam, TPConlyV0MinvALam);

  // Fill a minv hist w/o any cuts. Select background from the sideband
  if(fkDoStdHists){
    fHistSideBandOnLam->Fill(minvLam);
    fHistSideBandOnALam->Fill(minvALam);
  }
  // Fill the event buffer w/ background if enabled
  if(fkDoBgLamALam){
    if (fkUseOnTheFly){
      // Select side band aka background lambdas
      if (TMath::Abs(minvLam - fkLamMass) > 0.015 &&
	  TMath::Abs(minvLam - fkLamMass) < 0.035 ){
	
	//      if(!fdoLamOnly)
	fFemtoBuffer->GetEvt(0)->AddBgLam(v0, pTrack, nTrack);
	// Momentum difference of standard V0 / TPC only V0
	//      fHistMomDiffBgLam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
	//			    V0Mom[1] - TPConlyV0Mom[1],
	//			    V0Mom[2] - TPConlyV0Mom[2]);
	// Same excluding V0s with daughters with SPD hits
	//      if( !(pTrack->HasPointOnITSLayer(0) || pTrack->HasPointOnITSLayer(1) ||
	//	    nTrack->HasPointOnITSLayer(0) || nTrack->HasPointOnITSLayer(1) )){
	// No SPD hits
	//	fHistMomDiffWoSPDBgLam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
	//				   V0Mom[1] - TPConlyV0Mom[1],
	//				   V0Mom[2] - TPConlyV0Mom[2]);
	
	//    }
      } // End of background lambdas
      // Select side band aka background anti-lambdas
      if ( TMath::Abs(minvALam - fkLamMass) > 0.015 &&
	   TMath::Abs(minvALam - fkLamMass) < 0.035 ){
	//      if(!fdoLamOnly)
	fFemtoBuffer->GetEvt(0)->AddBgALam(v0, pTrack, nTrack);
	// Momentum difference of standard V0 / TPC only V0
	//      fHistMomDiffBgALam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
	//			       V0Mom[1] - TPConlyV0Mom[1],
	//			       V0Mom[2] - TPConlyV0Mom[2]);
	// Same excluding V0s with daughters with SPD hits
	//      if( !(pTrack->HasPointOnITSLayer(0) || pTrack->HasPointOnITSLayer(1) ||
	//	    nTrack->HasPointOnITSLayer(0) || nTrack->HasPointOnITSLayer(1) )){
	// No SPD hits
	//	fHistMomDiffWoSPDBgALam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
	//			      V0Mom[1] - TPConlyV0Mom[1],
	//				      V0Mom[2] - TPConlyV0Mom[2]);
	//      } // No SPD hits
      } // End of background anti-lambda
    } // End of if use on-the-fly finder
  } // End of if bg lam enabled

  //
  // Require 80 TPC clusters for both daughters
  //
  // There's a lambda signal for 0-9 clusters of the proton 
  // as it's for 110-120?!
  // There was a bug in the finding of the global track, since 
  // fixing it, offline is fine (and the problem looks less 
  // severe for on-the-fly). Still there is a problem here. 
  // There are tracks with 0 clusters. This is not the case
  // for the offline finder. The speculation would be that 
  // 1-9 clusters are treated correctly also here, it's just
  // the 0 cluster tracks. Should be a filter issue: on-the-fly
  // finds a V0, stores the daughter but info doesn't get written.

  // SOLUTION, see http://savannah.cern.ch/bugs/?97544:
  // The cluster map was not always properly updated. Needs a
  // reconstrunction from raw data
  
  // if(pTrack->GetTPCNcls()){
  //   // More than zero clusters
  //   fHistTPCNclsPosOnLam->Fill(pTrack->GetTPCNcls(),minvLam);
  //   fHistTPCNclsPosOnALam->Fill(pTrack->GetTPCNcls(),minvALam);
  // }
  // else {
  //   // Zero clusters, fill the underflow to distinguish
  //   fHistTPCNclsPosOnLam->Fill(-1,minvLam);
  //   fHistTPCNclsPosOnALam->Fill(-1,minvALam);
  // }
  // if(nTrack->GetTPCNcls()){
  //   // More than zero clusters
  //   fHistTPCNclsNegOnLam->Fill(nTrack->GetTPCNcls(),minvLam);
  //   fHistTPCNclsNegOnALam->Fill(nTrack->GetTPCNcls(),minvALam);
  // }
  // else {
  //   // Zero clusters, fill the underflow to distinguish
  //   fHistTPCNclsNegOnLam->Fill(-1,minvLam);
  //   fHistTPCNclsNegOnALam->Fill(-1,minvALam);
  // }
  
  // Do the cut on the TPC clusters, 0 OR at least 80
  if ( ( pTrack->GetTPCNcls() < 80 && pTrack->GetTPCNcls() ) ||
       ( nTrack->GetTPCNcls() < 80 && nTrack->GetTPCNcls() ) ) 
    return;
  //  fHistGoodV0->Fill(3.0,kOnTheFlyStat);

  // Require a maximum dca of the daughters of 0.2cm
  // fHistDcaV0DaughtersOnLam->Fill(v0->DcaV0Daughters(),minvLam);
  // fHistDcaV0DaughtersOnALam->Fill(v0->DcaV0Daughters(),minvALam);
  // fHistDcaV0Daughters->Fill(v0->DcaV0Daughters(),kOnTheFlyStat);
  if (v0->DcaV0Daughters() > 0.2)
    return;
  //  fHistGoodV0->Fill(4.0,kOnTheFlyStat);
  
  // Require cosine of pointing angle bigger than 0.9999
  // fHistCosPointAngle->Fill(lCosPoint,kOnTheFlyStat);
  // fHistCosPointAngleZoom->Fill(lCosPoint,kOnTheFlyStat);
  // fHistCosPointLamOn->Fill(lCosPoint,minvLam);
  // fHistCosPointALamOn->Fill(lCosPoint,minvALam);
  // fHistCosPointLamZoomOn->Fill(lCosPoint,minvLam);
  // fHistCosPointALamZoomOn->Fill(lCosPoint,minvALam);
  if (lCosPoint<0.9999)
    return;
  //  fHistGoodV0->Fill(5.0,kOnTheFlyStat);
  // Force TPC PID to be present
  if (!(pTrack->GetStatus() & AliAODTrack::kTPCpid) ||
      !(nTrack->GetStatus() & AliAODTrack::kTPCpid)) {
    // No TPC pid present for this track
    return;
  }
  //  fHistGoodV0->Fill(6.0,kOnTheFlyStat);
  // Visualize TPC signal before performing selection
  // fHistPosTpcBeforeCut->Fill(pTrack->P(),pTrack->GetTPCsignal());
  // fHistNegTpcBeforeCut->Fill(nTrack->P(),nTrack->GetTPCsignal());
  // // The Nsigma distribution for TPC dE/dx
  // fHistPosNsigmaTpcOnLam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton)),minvLam);
  // fHistPosNsigmaTpcOnALam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)),minvALam);
  // fHistNegNsigmaTpcOnLam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)),minvLam);
  // fHistNegNsigmaTpcOnALam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton)),minvALam);

  // Perform cut on TPC dE/dx
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton)) > 3.7)
    minvLam=0.0;
  // else 
  //   fHistPosTpcAfterCut->Fill(pTrack->P(),pTrack->GetTPCsignal());
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)) > 3.8)
    minvLam=0.0;
  // else
  //   fHistNegTpcAfterCut->Fill(nTrack->P(),nTrack->GetTPCsignal());
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)) > 4.2)
    minvALam=0.0;
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton)) > 3.9)
    minvALam=0.0;

  // Don't use a tof cut for pions

  // Check whether to use a 5sigma tof cut or none for protons
  // if (pTrack->GetStatus() & AliAODTrack::kTOFpid){
  //   if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kProton)) > 5.0)
  //     fHistUseTofOnLam->Fill(1.0,minvLam);
  //   else
  //     fHistUseTofOnLam->Fill(0.0,minvLam);
  // }
  // else
  //   fHistUseTofOnLam->Fill(0.0,minvLam);
  // // Check whether to use a 5sigma tof cut or none for anti-protons
  // if (nTrack->GetStatus() & AliAODTrack::kTOFpid){
  //   if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kProton)) > 5.0)
  //     fHistUseTofOnALam->Fill(1.0,minvALam);
  //   else
  //     fHistUseTofOnALam->Fill(0.0,minvALam);
  // }
  // else
  //   fHistUseTofOnALam->Fill(0.0,minvALam);

  // Reject (anti-)protons with more than 5sigma TOF
  if (nTrack->GetStatus() & AliAODTrack::kTOFpid){
    if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kProton)) > 5.0)
      minvALam=0.0;
  }
  if (pTrack->GetStatus() & AliAODTrack::kTOFpid){
    if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kProton)) > 5.0)
      minvLam=0.0;
  }
  
  // Don't require TPC refit. You would kill nearly your whole signal  

  // Distance between daughters and primary vertex
  // fHistDcaPosToPrimVertex->Fill(v0->DcaPosToPrimVertex(),kOnTheFlyStat);
  // fHistDcaNegToPrimVertex->Fill(v0->DcaNegToPrimVertex(),kOnTheFlyStat);
  // fHistDcaPosToPrimVertexZoom->Fill(v0->DcaPosToPrimVertex(),kOnTheFlyStat);
  // fHistDcaNegToPrimVertexZoom->Fill(v0->DcaNegToPrimVertex(),kOnTheFlyStat);
  // fHistDcaPosOnLam->Fill(v0->DcaPosToPrimVertex(),minvLam);
  // fHistDcaPosOnALam->Fill(v0->DcaPosToPrimVertex(),minvALam);
  // fHistDcaNegOnLam->Fill(v0->DcaNegToPrimVertex(),minvLam);
  // fHistDcaNegOnALam->Fill(v0->DcaNegToPrimVertex(),minvALam);
  // Require at least 0.02 cm distance from the primary vertex for the (anti-)protons
  if (v0->DcaPosToPrimVertex() < 0.02)
    minvLam=0.0;
  if (v0->DcaNegToPrimVertex() < 0.02)
    minvALam=0.0;
  // Require at least 0.05 cm distance from the primary vertex for the pions
  if (v0->DcaPosToPrimVertex() < 0.05)
    minvALam=0.0;
  if (v0->DcaNegToPrimVertex() < 0.05)
    minvLam=0.0;
  
  // Fill some histograms with cut variables
  //  fHistChi2->Fill(v0->Chi2V0(),kOnTheFlyStat);
  
  
  // Idea to cut on the radius
  // fHistRadiusV0->Fill(v0->RadiusV0(),kOnTheFlyStat);
  // fHistV0RadiusLamOn->Fill(v0->RadiusV0(),minvLam);
  // fHistV0RadiusALamOn->Fill(v0->RadiusV0(),minvALam);
  
  // Idea to cut on the decay length
  // fHistDecayLengthV0->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),kOnTheFlyStat);
  // fHistV0DecayLengthLamOn->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),minvLam);
  // fHistV0DecayLengthALamOn->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),minvALam);

  // Idea to cut on DCA of V0 and primay vertex
  // fHistDcaV0PriVertexLamOn->Fill(v0->DcaV0ToPrimVertex(),minvLam);
  // fHistDcaV0PriVertexALamOn->Fill(v0->DcaV0ToPrimVertex(),minvALam);

  // TPC Chi2 / number of degrees of freedom
  // A cut on at least 80 clusters is already done before,
  // no concern to divide by zero
  // fHistChi2TPCPosLamOn->Fill(pTrack->Chi2perNDF(),minvLam);
  // fHistChi2TPCPosALamOn->Fill(pTrack->Chi2perNDF(),minvALam);
  // fHistChi2TPCNegLamOn->Fill(nTrack->Chi2perNDF(),minvLam);
  // fHistChi2TPCNegALamOn->Fill(nTrack->Chi2perNDF(),minvALam);
  // Don't cut like Chi2/ndf < 4! One might throw away the tracks
  // with Chi2/ndf roughly one as they are good primaries

  // Fill some invariant mass distributions
  if(fkDoStdHists){
    fHistMassLambdaOn->Fill(minvLam);
    fHistMassAntiLambdaOn->Fill(minvALam);
  }
  // fHistPtVsMassLambdaOn->Fill(v0->Pt(),minvLam);
  // fHistPtVsMassAntiLambdaOn->Fill(v0->Pt(),minvALam);

  // TPC only invariant mass distributions
  //  if(minvLam > .1){
    // Lambda is good
    //    fHistMinvTPConlyLamOn->Fill(TPConlyV0MinvLam);
  //  }
  //  if (minvALam > .1){
    // Anti-lambda is good
    //    fHistMinvTPConlyALamOn->Fill(TPConlyV0MinvALam);
  //  }
  
  // 3d histogram: rapidity, pt and mass
  if(fkDoStdHists){
    fHistYPtMassLamOn->Fill(v0->Y(3122),v0->Pt(),minvLam);
    fHistYPtMassALamOn->Fill(v0->Y(-3122),v0->Pt(),minvALam);
  }
  // // Invariant mass cut lambda :: fill a y-pt hists
  // if ( TMath::Abs(minvLam - fkLamMass) < 0.01 ){
  //   fHistPtVsYLambdaOn->Fill(v0->Pt(),v0->Y(3122));
  // }
  // // Invariant mass cut anti-lambda :: fill a y-pt hists
  // if ( TMath::Abs(minvALam - fkLamMass) < 0.01 ){
  //   fHistPtVsYAntiLambdaOn->Fill(v0->Pt(),v0->Y(-3122));
  // }
  
  // Fill the mixed events when on-the-fly V0 finder is used
  if (fkUseOnTheFly){

    // Highest significance for minv +/- 4 MeV
    if ( TMath::Abs(minvLam - fkLamMass) < fkMinvCut ){

      // We introduce a pt cut - only on lambdas! - of 0.5 GeV/c.
      // The reason for this is that below 0.5 GeV/c we have a 
      // huge background which will only ruin our whole lambda
      // sample. There's also the highest track density, so ex-
      // cluding this region will in addition probably help us 
      // fighting two-track effects. The signal is tiny, so we
      // don't loose much. The reason for the background is pro-
      // bably the protons from material?!
      if(v0->Pt() >= 0.5){
	fFemtoBuffer->GetEvt(0)->AddLam(v0, pTrack, nTrack);
      }

      // Momentum difference of standard V0 / TPC only V0
      //      fHistMomDiffLam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //		    V0Mom[1] - TPConlyV0Mom[1],
      //		    V0Mom[2] - TPConlyV0Mom[2]);
      // Same excluding V0s with daughters with SPD hits
      //      if( !(pTrack->HasPointOnITSLayer(0) || pTrack->HasPointOnITSLayer(1) ||
      //	    nTrack->HasPointOnITSLayer(0) || nTrack->HasPointOnITSLayer(1) )){
      //	// No SPD hits
      //	fHistMomDiffWoSPDLam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //				   V0Mom[1] - TPConlyV0Mom[1],
      //				   V0Mom[2] - TPConlyV0Mom[2]);
      // } // No SPD hits
    } // Good lambda
    if ( TMath::Abs(minvALam - fkLamMass) < fkMinvCut ) {

      fFemtoBuffer->GetEvt(0)->AddALam(v0, pTrack, nTrack);
      // Momentum difference of standard V0 / TPC only V0
      //      fHistMomDiffALam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //			    V0Mom[1] - TPConlyV0Mom[1],
      //			    V0Mom[2] - TPConlyV0Mom[2]);
      // Same excluding V0s with daughters with SPD hits
      //      if( !(pTrack->HasPointOnITSLayer(0) || pTrack->HasPointOnITSLayer(1) ||
      //	    nTrack->HasPointOnITSLayer(0) || nTrack->HasPointOnITSLayer(1) )){
	// No SPD hits
      //	fHistMomDiffWoSPDALam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //				   V0Mom[1] - TPConlyV0Mom[1],
      //				   V0Mom[2] - TPConlyV0Mom[2]);
      //    } // No SPD hits
    } // Good anti-lambda
  } // Use on-the-fly finder for Femto analysis
} // ProcessOnTheFly
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::ProcessTOF(const AliAODTrack* track) 
{
  // Request the kTOFpid bit. There are tracks with kTOFout and wihthout kTOFpid,
  // but these tracks have a bad TOF signal.
  if(!((fGTI[-track->GetID()-1])->GetStatus() & AliAODTrack::kTOFpid))
    return;

  // TOF signal corrected for expected time and (if neccessary) for start time
  const Float_t corrTOFsig = GetCorrectedTOFSignal(track);

  // We do a mild dEdx cut for momenta < 1.25 as there's still a slight mis-
  // match problem.
  if(track->P()<1.25){
    // Check for TPC pid
    if((fGTI[-track->GetID()-1])->GetStatus() & AliAODTrack::kTPCpid){
      // Compare signal with bb parametrization
      if ((fGTI[-track->GetID()-1])->GetTPCsignal() < 
	  fTpcResponse->GetExpectedSignal((fGTI[-track->GetID()-1])->GetTPCmomentum(),
					  AliPID::kProton)
	  - 10. -35.*((fGTI[-track->GetID()-1])->GetTPCmomentum()-1.) ) {
	// Track fails dEdx cut. (To view the cut see findProCutsTPC2013.C )
	// We get less strict with the cut at higher momenta, i.e. ~1.25 GeV/c as we
	// really loose any analyzing power with the dEdx and can't see what we're doing,
	// so we cut less strict. Also, the most yield for the 1.0 - 1.25 bin sits at 
	// 1.0 GeV. In addition, the mismatch is a low pt problem. So if we cut a lot
	// around 1.0 GeV, we're done, as that's the mismatch we want to get rid of. :)
	return;    
      } // End of fails dEdx cut
    } // End of TPC pid
  } // End of p < 1.25
  
  // Distinguish between charges
  if (track->Charge() > 0){
    // Simple Nsigma TOF distribution
    //    fPriHistPosNsigmaTof->Fill(fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    // Nsigma TOF in bins of total momentum
    // fPriHistPosNsigmaTofVsP->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    // Nsigma TOF in bins of transverse momentum
    //    fPriHistPosNsigmaTofVsPt->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    
    // Try the tof signal instead of nsigma
    if(fkDoStdHists)
      fPriHistTOFsignalPosVsP->Fill(track->P(), corrTOFsig);
    //    fPriHistTOFsignalPosVsPt->Fill(track->Pt(), corrTOFsig);
    
  }
  else if (track->Charge() < 0){
    // Simple Nsigma TOF distribution
    //    fPriHistNegNsigmaTof->Fill(fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    // Nsigma TOF in bins of total momentum
    // fPriHistNegNsigmaTofVsP->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    // Nsigma TOF in bins of transverse momentum
    //    fPriHistNegNsigmaTofVsPt->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    
    // Try the tof signal instead of nsigma
    if(fkDoStdHists)
      fPriHistTOFsignalNegVsP->Fill(track->P(), corrTOFsig);
    //    fPriHistTOFsignalNegVsPt->Fill(track->Pt(), corrTOFsig);
  }
  

  // We use an upper cut on the TOF time. 
  // It's mainly to have a well defined 
  // region to determine the purity of
  // the proton sample.
  if(corrTOFsig >= 500.){
    return;
  }
  
  // Here comes the TOF cut. It was determined
  // via fits to the TOF spectrum and the cuts
  // result in a proton purity >= 99.0%
  if(track->P()<2. && corrTOFsig < -500.)
    return;
  else if(track->P()<2.25 && corrTOFsig < -350.)
    return;
  else if(track->P()<2.50 && corrTOFsig < -300.)
    return;
  else if(track->P()<2.75 && corrTOFsig < -250.)
    return;
  else if(track->P()<3.00 && corrTOFsig < -200.)
    return;
  else if(track->P()<3.75 && corrTOFsig < -150.)
    return;
  else if(track->P()<4.25 && corrTOFsig < -100.)
    return;
  else if(track->P()<4.50 && corrTOFsig < -50.)
    return;
  else if(track->P()<5.00 && corrTOFsig < 0.)
    return;

  // We don't go beyond 5.0 GeV/c momentum yet
  if(track->P()>5.00)
    return;

  if (track->Charge()>0){
    // Cut on DCAxy and fill a histogram
    if(goodDCA(track)){
      
      // We reject protons with 
      // pt<.5 && |y|<.5 as in this region
      // we have purities of about 40%.
      // The significance grows if one 
      // rejects this region. See 
      // ~/alice/wikigrid/analyzeOutput/Fit2dDCA.C
      if(track->Pt()<.5 &&
	 TMath::Abs(RapidityProton(track)<.5)){
	return;
      }

      // Add to the femto event
      fFemtoBuffer->GetEvt(0)->AddPro(track);
    }
  }
  else{
    // Cut on DCAxy and fill a histogram
    if(goodDCA(track)){
      // We reject anti-protons at large 
      // pseudorapidity. Probably the material
      // is not well described in this region
      // plus the fact that geant3 doesn't 
      // describe well the anti-proton cross-
      // sections. See, e.g., Marco van Leeuwen's
      // talk at the "LHC Detector Simulation 
      // workshop", 6-7 October 2011:
      // http://indico.cern.ch/getFile.py/access?contribId=4&sessionId=0&resId=0&materialId=slides&confId=144956
      // The result is that the DCA_xy and DCA_z 
      // distributions are not well described. 
      // Additionally the statistics in MC are 
      // not the nicest in this region. See also
      // ~/alice/wikigrid/analyzeOutput/Fit2dDCA.C
      if(track->Pt()<0.5 && 
	 TMath::Abs(RapidityProton(track)>.5)){
	return;
      }

      // Add to the femto event
      fFemtoBuffer->GetEvt(0)->AddAPro(track);
    }
  }
} // End of void ProcessTOF
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::ProcessTPC(const AliAODTrack* track){

  // Require the TPCpid bit
  if (!((fGTI[-track->GetID()-1])->GetStatus() & AliAODTrack::kTPCpid))
    return;
    
  // In contrast to ESDs one doesn't check for AliESDtrack::kTOFpid
  // but for AliAODTrack::kTOFout?? 
  // Check how many particles have TOFout bit
  // if (track->Charge() > 0){
  //   if ((fGTI[-track->GetID()-1])->GetStatus() & AliAODTrack::kTOFpid)
  //     fPriHistHasTofPos->Fill(1.0);
  //   else 
  //     fPriHistHasTofPos->Fill(0.0);
  // }
  // else{
  //   if ((fGTI[-track->GetID()-1])->GetStatus() & AliAODTrack::kTOFpid)
  //     fPriHistHasTofNeg->Fill(1.0);
  //   else 
  //     fPriHistHasTofNeg->Fill(0.0);
  // }

  // For all plots <dE/dx> vs p one should use
  // the momentum at the inner wall of the TPC.

  // Use a TOF cut and fill the same dE/dx histograms
  // Bool_t acceptedTOF=kFALSE;
  // if ((fGTI[-track->GetID()-1])->GetStatus() & AliAODTrack::kTOFpid){
  //   if (fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton) > -10.0)
  //     acceptedTOF=kTRUE;
  // }
  // if (acceptedTOF){
  //     if (track->Charge() > 0){
  // 	fPriHistTPCsignalTOFcutPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
  // 					 (fGTI[-track->GetID()-1])->GetTPCsignal());
  // 	fPriHistNsigmaTPCTOFcutPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
  // 					 fPIDResponse->NumberOfSigmasTPC((fGTI[-track->GetID()-1]), AliPID::kProton));
  //     }
  //     else{
  // 	fPriHistTPCsignalTOFcutNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
  // 					 (fGTI[-track->GetID()-1])->GetTPCsignal());
  // 	fPriHistNsigmaTPCTOFcutNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
  // 					 fPIDResponse->NumberOfSigmasTPC((fGTI[-track->GetID()-1]), AliPID::kProton));
  //     }
  // }	
    
  // A first idea of a cut: use the spectra cut.
  // (should perhaps change for momenta ~ 0.75 GeV)
  if ( ((fGTI[-track->GetID()-1])->GetTPCsignal() > 
	0.95 * fTpcResponse->GetExpectedSignal((fGTI[-track->GetID()-1])->GetTPCmomentum(),
					AliPID::kProton))
      // New since Sept 10th 2012: Also use a cut to reject deuterons.
      // I checked: The cut is good!
       && ((fGTI[-track->GetID()-1])->GetTPCsignal() <
	 0.95 * 2.0*fTpcResponse->GetExpectedSignal((fGTI[-track->GetID()-1])->GetTPCmomentum(),
				      AliPID::kProton))
                        ) {
    // Distinguish between charges
    if (track->Charge()>0){
      // Cut on DCAxy and fill a histogram
      if(goodDCA(track)){
	// Add to the femto event
	fFemtoBuffer->GetEvt(0)->AddPro(track);
      }
    }
    else{
      // Cut on DCAxy and fill a histogram
      if(goodDCA(track)){
	// Add to the femto event
	fFemtoBuffer->GetEvt(0)->AddAPro(track);
      }
    }
  }
} // End of void ProcessTPC
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::ProcessHybrid(const AliAODTrack *track){
  
  // Intermediate momentum: use dEdx for a pre-selection
  // and do the pid with tof
  
  // Boolean for extra! tpc pid cuts
  Bool_t acceptTPC = kTRUE;

  // Require the TPCpid bit
  if (!((fGTI[-track->GetID()-1])->GetStatus() & AliAODTrack::kTPCpid))
    acceptTPC = kFALSE;
 
  // Pre-selection cut with TPC, don't return immediately to be able
  // to visualize the effect
  if (acceptTPC){
    // Do a mild dEdx cut
    if ((fGTI[-track->GetID()-1])->GetTPCsignal() < 
	fTpcResponse->GetExpectedSignal((fGTI[-track->GetID()-1])->GetTPCmomentum(),
					AliPID::kElectron) *
	(((fGTI[-track->GetID()-1])->GetTPCmomentum() > 0.85)?0.9:1.0)) // from 0.85 to 1.0 GeV the cut was a little too strict
      acceptTPC = kFALSE;
  }
    
  // Ask for TOF pid flag and fill
  if (!((fGTI[-track->GetID()-1])->GetStatus() & AliAODTrack::kTOFpid))
    return;
  
  // The corrected TOF signal
  const Float_t corrTOFsig = GetCorrectedTOFSignal(track);
  
  // Distinguish between charges
  if (track->Charge() > 0) {
    // Fill the tof signal w/o dedx pre-selection
    if(fkDoStdHists)
      fPriHistHybridTOFsigPosWoTPC->Fill(corrTOFsig);
    // Do the pre-selection
    if (acceptTPC){
      if(fkDoStdHists)
	fPriHistHybridTOFsigPosTPCok->Fill(corrTOFsig);

      // Do the tof cut
      // Sept '12: also include an upper cut
      if ( (corrTOFsig > -1000.0) && (corrTOFsig < 1250.) ){
	// Create additional TPC only constrained to pri. vtx track parameters
	//	constrainTrack(track);
	// Cut on DCAxy and fill a histogram
	if(goodDCA(track)){
	  // Add to the femto event
	  fFemtoBuffer->GetEvt(0)->AddPro(track);
	}
      }
    }
  }
  else {
    // Fill the tof signal w/o dedx pre-selection
    if(fkDoStdHists)
      fPriHistHybridTOFsigNegWoTPC->Fill(corrTOFsig);
    // Do the pre-selection
    if (acceptTPC){
      if(fkDoStdHists)
	fPriHistHybridTOFsigNegTPCok->Fill(corrTOFsig);
      
      // Do the tof cut
      // Sept '12: also include an upper cut
      if ( (corrTOFsig > -1000.0) && (corrTOFsig < 1250.) ){
	// Create additional TPC only constrained to pri. vtx track parameters
	//	constrainTrack(track);
	// Cut on DCAxy and fill a histogram
	if(goodDCA(track)){
	  // add to the femto event
	  fFemtoBuffer->GetEvt(0)->AddAPro(track);
	}
      }
    }
  }
} // End of ProcessHybrid
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::CleaningProcedure() {
  // fFemtoBuffer->GetEvt(0) pointer must be set
  // Checks that no tracks are shared between Lam & Lam, Lam & Pro, ALam & ALam, ALam & APro

  // printf ("Cleaning procedure. Lam: %d, ALam: %d, Pro: %d, APro:%d\n"
  // 	  ,fFemtoBuffer->GetEvt(0)->GetNLam(),fFemtoBuffer->GetEvt(0)->GetNALam(),fFemtoBuffer->GetEvt(0)->GetNPro(),fFemtoBuffer->GetEvt(0)->GetNAPro());

  // Monitor the number of particles in the buffer, 
  // so the buffer size can be adjusted
  if(fkDoStdHists){ 
    f2HistNLamBefClean->Fill(fFemtoBuffer->GetEvt(0)->GetNLam());
    f2HistNProBefClean->Fill(fFemtoBuffer->GetEvt(0)->GetNPro());
    f2HistNALamBefClean->Fill(fFemtoBuffer->GetEvt(0)->GetNALam());
    f2HistNAProBefClean->Fill(fFemtoBuffer->GetEvt(0)->GetNAPro());
  }

  //
  // Check for lambdas..
  //
  for (Int_t i=0;i<fFemtoBuffer->GetEvt(0)->GetNLam();i++) {
    if (!fFemtoBuffer->GetEvt(0)->fLamTracks[i].UseIt())
      continue;
    // Unique track ids for first V0
    const UShort_t posId1 = fFemtoBuffer->GetEvt(0)->fLamTracks[i].fPosDaughter.fID;
    const UShort_t negId1 = fFemtoBuffer->GetEvt(0)->fLamTracks[i].fNegDaughter.fID;

    // .. & lambdas
    for (Int_t j=i+1;j<fFemtoBuffer->GetEvt(0)->GetNLam();j++){
      if (!fFemtoBuffer->GetEvt(0)->fLamTracks[j].UseIt())
  	continue;
      // Unique track ids for second V0
      const UShort_t posId2 = fFemtoBuffer->GetEvt(0)->fLamTracks[j].fPosDaughter.fID;
      const UShort_t negId2 = fFemtoBuffer->GetEvt(0)->fLamTracks[j].fPosDaughter.fID;
      
      // If V0s share a track remove one
      if (posId1 == posId2 || negId1 == negId2){

	// printf ("shared track lamlam! posId1: %d, posId2: %d, negId1: %d, negId2: %d\n",
  	// 	posId1, posId2, negId1, negId2);
	
  	// Use a criterion to select best V0
  	if (fFemtoBuffer->GetEvt(0)->fLamTracks[i].fCosPoint > fFemtoBuffer->GetEvt(0)->fLamTracks[j].fCosPoint){
  	  fFemtoBuffer->GetEvt(0)->fLamTracks[j].SetBadFlag();
  	}
  	else{
  	  fFemtoBuffer->GetEvt(0)->fLamTracks[i].SetBadFlag();
  	}
      }
      
    } // Scnd V0 loop

    if (!fFemtoBuffer->GetEvt(0)->fLamTracks[i].UseIt())
      continue;

    // .. & protons
    for (Int_t j=0;j<fFemtoBuffer->GetEvt(0)->GetNPro();j++){
      if (!fFemtoBuffer->GetEvt(0)->fProTracks[j].UseIt())
  	continue;
      // Unique track ids for proton
      const UShort_t posId2 = fFemtoBuffer->GetEvt(0)->fProTracks[j].fID;
      
      // If V0 and proton share a track
      if (posId1 == posId2){
	//  	printf ("shared track lam p! id:%d\n",posId1);
  	
	// Remove the proton
  	fFemtoBuffer->GetEvt(0)->fProTracks[j].SetBadFlag();
      }
      
    } // Proton loop

  } // First V0 loop

  //
  // Check for anti-lambdas..
  //
  for (Int_t i=0;i<fFemtoBuffer->GetEvt(0)->GetNALam();i++){
    if (!fFemtoBuffer->GetEvt(0)->fALamTracks[i].UseIt())
      continue;
    // Unique track ids for first V0
    const UShort_t posId1 = fFemtoBuffer->GetEvt(0)->fALamTracks[i].fPosDaughter.fID;
    const UShort_t negId1 = fFemtoBuffer->GetEvt(0)->fALamTracks[i].fNegDaughter.fID;

    // .. & anti-lambdas
    for (Int_t j=i+1;j<fFemtoBuffer->GetEvt(0)->GetNALam();j++){
      if (!fFemtoBuffer->GetEvt(0)->fALamTracks[j].UseIt())
  	continue;
      // Unique track ids for second V0
      const UShort_t posId2 = fFemtoBuffer->GetEvt(0)->fALamTracks[j].fPosDaughter.fID;
      const UShort_t negId2 = fFemtoBuffer->GetEvt(0)->fALamTracks[j].fNegDaughter.fID;
      
      // If V0s share a track remove one
      if (posId1 == posId2 || negId1 == negId2){
	
	// printf ("shared track ALamALam! posId1: %d, posId2: %d, negId1: %d, negId2: %d\n",
  	// 	posId1, posId2, negId1, negId2);

  	// Use a criterion to select best V0
  	if (fFemtoBuffer->GetEvt(0)->fALamTracks[i].fCosPoint > fFemtoBuffer->GetEvt(0)->fALamTracks[j].fCosPoint){
  	  fFemtoBuffer->GetEvt(0)->fALamTracks[j].SetBadFlag();
  	}
  	else{
  	  fFemtoBuffer->GetEvt(0)->fALamTracks[i].SetBadFlag();
  	}
      }
      
    } // Scnd anti-V0 loop

    if (!fFemtoBuffer->GetEvt(0)->fALamTracks[i].UseIt())
      continue;
    
    // .. & anti-protons
    for (Int_t j=0;j<fFemtoBuffer->GetEvt(0)->GetNAPro();j++){
      if (!fFemtoBuffer->GetEvt(0)->fAProTracks[j].UseIt())
  	continue;
      // Unique track ids for anti-proton
      const UShort_t negId2 = fFemtoBuffer->GetEvt(0)->fAProTracks[j].fID;
      
      // If V0 and proton share a track
      if (negId1 == negId2){
	//  	printf ("shared track alam ap! id:%d\n",posId1);

  	// Remove the proton
  	fFemtoBuffer->GetEvt(0)->fAProTracks[j].SetBadFlag();
      }
      
    } // Anti-proton loop

  } // First anti-V0 loop

  //
  // Do the same with the side band background.
  // Discard background when sharing track with primary proton.
  //
  if(fkDoBgLamALam){
    for (Int_t i=0;i<fFemtoBuffer->GetEvt(0)->GetNBgLam();i++){
      if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].UseIt())
	continue;
      // Unique track id's for first V0
      const UShort_t posId1 = fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].fPosDaughter.fID;
      const UShort_t negId1 = fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].fNegDaughter.fID;

      // .. & lambdas
      for (Int_t j=i+1;j<fFemtoBuffer->GetEvt(0)->GetNBgLam();j++){
	if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[j].UseIt())
	  continue;
	// Unique track id's for second V0
	const UShort_t posId2 = fFemtoBuffer->GetEvt(0)->fBgLamTracks[j].fPosDaughter.fID;
	const UShort_t negId2 = fFemtoBuffer->GetEvt(0)->fBgLamTracks[j].fNegDaughter.fID;
      
	// If V0s share a track remove one
	if (posId1 == posId2 || negId1 == negId2){

	  // printf ("shared track bglambglam! posId1: %d, posId2: %d, negId1: %d, negId2: %d\n",
	  // 	posId1, posId2, negId1, negId2);
	
	  // Use a criterion to select best V0
	  if (fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].fCosPoint > fFemtoBuffer->GetEvt(0)->fBgLamTracks[j].fCosPoint){
	    fFemtoBuffer->GetEvt(0)->fBgLamTracks[j].SetBadFlag();
	  }
	  else{
	    fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].SetBadFlag();
	  }
	}
      
      } // Scnd V0 loop

      if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].UseIt())
	continue;

      // .. & protons
      for (Int_t j=0;j<fFemtoBuffer->GetEvt(0)->GetNPro();j++) {
	if (!fFemtoBuffer->GetEvt(0)->fProTracks[j].UseIt())
	  continue;
	if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].UseIt())
	  continue;

	// Unique track id's for second V0
	const UShort_t posId2 = fFemtoBuffer->GetEvt(0)->fProTracks[j].fID;
      
	// If V0 and proton share a track
	if (posId1 == posId2){
	  //  	printf ("shared track bglam p! id:%d\n",posId1);
	  // Remove the background lambda
	  fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].SetBadFlag();
	}
      
      } // Proton loop

    } // First V0 loop

    //
    // Check for anti-lambdas..
    //
    for (Int_t i=0;i<fFemtoBuffer->GetEvt(0)->GetNBgALam();i++){
      if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].UseIt())
	continue;
      // Unique track id's for first V0
      const UShort_t posId1 = fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].fPosDaughter.fID;
      const UShort_t negId1 = fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].fNegDaughter.fID;

      // .. & anti-lambdas
      for (Int_t j=i+1;j<fFemtoBuffer->GetEvt(0)->GetNBgALam();j++){
	if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[j].UseIt())
	  continue;
	// Unique track id's for second V0
	const UShort_t posId2 = fFemtoBuffer->GetEvt(0)->fBgALamTracks[j].fPosDaughter.fID;
	const UShort_t negId2 = fFemtoBuffer->GetEvt(0)->fBgALamTracks[j].fNegDaughter.fID;
      
	// If V0s share a track remove one
	if (posId1 == posId2 || negId1 == negId2){
	
	  // printf ("shared track BgALamBgALam! posId1: %d, posId2: %d, negId1: %d, negId2: %d\n",
	  // 	posId1, posId2, negId1, negId2);

	  // Use a criterion to select best V0
	  if (fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].fCosPoint > fFemtoBuffer->GetEvt(0)->fBgALamTracks[j].fCosPoint){
	    fFemtoBuffer->GetEvt(0)->fBgALamTracks[j].SetBadFlag();
	  }
	  else{
	    fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].SetBadFlag();
	  }
	}
      
      } // Scnd anti-V0 loop

      if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].UseIt())
	continue;
    
      // .. & anti-protons
      for (Int_t j=0;j<fFemtoBuffer->GetEvt(0)->GetNAPro();j++){
	if (!fFemtoBuffer->GetEvt(0)->fAProTracks[j].UseIt())
	  continue;
      
	// Unique track id's for second V0
	const UShort_t negId2 = fFemtoBuffer->GetEvt(0)->fAProTracks[j].fID;
      
	// If V0 and proton share a track
	if (negId1 == negId2){
	  //  	printf ("shared track bgalam ap! id:%d\n",posId1);
	  // Remove the background anti-lambda
	  fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].SetBadFlag();
	}
      
      } // Anti-proton loop

    } // First anti-V0 loop
  }// End of if DoBgLamALam 
  
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::ProcessReal() {
  // Process real events
  
  // Declare numbers
  UChar_t iLam,iPro,iALam,iAPro;
  Double_t x[4];
  // printf("nLam %u nPro %u nALam %u nAPro %u\n"
  // 	 ,fFemtoBuffer->GetEvt(0)->GetNLam()
  // 	 ,fFemtoBuffer->GetEvt(0)->GetNPro()
  // 	 ,fFemtoBuffer->GetEvt(0)->GetNALam()
  // 	 ,fFemtoBuffer->GetEvt(0)->GetNAPro());

  // Event centrality
  const Float_t evtCent(fAOD->GetCentrality()->GetCentralityPercentileUnchecked(fkCentEst.Data()));
  
  //
  // Lambda loop
  //
  for (iLam = 0; iLam < fFemtoBuffer->GetEvt(0)->GetNLam(); iLam++) {
    
    // Skip if track flagged by cleaning procedure
    if(!fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].UseIt())
      continue;
    // Skip if track wasn't well propagated to R=1.2(5)m
    if(!fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter.GoodPropR12())
      continue;

    // Lambda properties:
    // .. phase space
    const Float_t rapLam(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].RapLam())
      ,pTLam(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].Pt());
    // .. purity
    const Float_t lamPur(GetLamPur(rapLam,pTLam))
      ,lamFdPur(GetLamFdPur(rapLam,pTLam,evtCent));      
    
    //
    // Proton loop
    //
    for (iPro=0;iPro<fFemtoBuffer->GetEvt(0)->GetNPro();iPro++){

      // Skip if track flagged by cleaning procedure
      if (!fFemtoBuffer->GetEvt(0)->fProTracks[iPro].UseIt())
  	continue;
      // Skip if track wasn't well propagated to R=1.2(5)m
      if (!fFemtoBuffer->GetEvt(0)->fProTracks[iPro].GoodPropR12())
  	continue;
      
      // Calculate the pair variables
      x[0]=mt(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
	      fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
      x[1]=dEtaS(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter
		,fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
      x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter
		      ,fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
      x[3]=Qinv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
		fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);


      if(fkDoStdHists){
	// 2d distance at r=1.2
	f2HistLamProAngDistSft2dAtR12Real->Fill(x[1],x[2]);
	// Mt of the pair
	f2HistMtLamProReal->Fill(x[0]);
      }

      // THn with qinvpropro, |dEtaS|, dPhiSR12, qinv lampro
      x[1]=TMath::Abs(x[1]);
      fLamProReal->Fill(x);

      // Since April 4th 2013 do the 99% TTR cut here.
      // x[1] is |dEtaS|, x[2] is dPhiS
      if(TMath::Power(x[1]/fdEtaSCut,2) + TMath::Power(x[2]/fdPhiSCut,2) < 1.){
	// Too close, reject
	continue;
      }
      
      // Proton properties:
      // .. phase space
      const Float_t rapPro(fFemtoBuffer->GetEvt(0)->fProTracks[iPro].RapPro())
	,pTPro(fFemtoBuffer->GetEvt(0)->fProTracks[iPro].Pt());
      // .. purity
      const Float_t proFdPur(GetProFdPur(rapPro,pTPro));

      // Fill four histograms: three with purity as a weight,
      // one without, vs qinv
      if(fkDoStdHists){
	f2HistLamProWLamPurReal->Fill(x[3],lamPur);
	f2HistLamProWLamFdPurReal->Fill(x[3],lamFdPur);
	f2HistLamProWProFdPurReal->Fill(x[3],proFdPur);
	f2HistLamProWoPurReal->Fill(x[3]);
      }
      // Also (vs qinv) vs mT
      f2HistLamProWLamPurmTReal->Fill(x[3],x[0],lamPur);
      f2HistLamProWLamFdPurmTReal->Fill(x[3],x[0],lamFdPur);
      f2HistLamProWProFdPurmTReal->Fill(x[3],x[0],proFdPur);
      f2HistLamProWoPurmTReal->Fill(x[3],x[0]);


      // Get purity fluctuations for qinv .0-.2
      if(x[3] <.2){
	f2HistPurFlucLamPro->Fill(lamPur);
      }

      // Origin in phase space (y,pt)
      if(fkDoStdHists){
	// For low q pairs, qinv < .2 GeV/c
	if(x[3] <.2){
	  // .. for lambdas
	  f2HistLamYPtPair->Fill(rapLam,pTLam);
	  // .. protons
	  f2HistProYPtPair->Fill(rapPro,pTPro);
	  // ..pairs
	  f2HistYKtPair->Fill(yPair(fFemtoBuffer->GetEvt(0)->fProTracks[iPro],
				    fFemtoBuffer->GetEvt(0)->fLamTracks[iLam]),
			      TMath::Sqrt(ktSquared(fFemtoBuffer->GetEvt(0)->fProTracks[iPro],
						    fFemtoBuffer->GetEvt(0)->fLamTracks[iLam])));
	  // mT of the low q pairs
	  f2HistMtLowQLamProReal->Fill(x[0]);
	  // Also  mT vs phi* (x[2]) vs for lowq pairs
	  f2HistMtVsPhiSLowQ->Fill(x[2],x[0]);
	} // End of low q < .2
	// Origin in phase space (y,pt) for each pair with .4 < qinv < .5 GeV/c
	else if(x[3] < .5 && x[3] >= .4){
	  // just for lambdas
	  f2HistLamYPtHiPair->Fill(rapLam,pTLam);
	}
      } // End of origin in phase space

    }// Proton loop
  }// Lambda loop

  //
  // Anti-lambda loop
  //
  for (iALam = 0; iALam < fFemtoBuffer->GetEvt(0)->GetNALam(); iALam++){

    // Skip if track was flagged by cleaning procedure
    if (!fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].UseIt())
      continue;
    // Skip if track was not well propagated to R=1.2(5)m
    if (!fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter.GoodPropR12())
      continue;
    
    // Anti-lambda properties..
    // ..phase space
    const Float_t rapALam(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].RapLam())
      ,pTALam(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].Pt());
    // ..purity
    const Float_t alamPur(GetALamPur(rapALam,pTALam))
      ,alamFdPur(GetALamFdPur(rapALam,pTALam,evtCent));
    
    //
    // AProton loop
    //
    for (iAPro=0;iAPro<fFemtoBuffer->GetEvt(0)->GetNAPro();iAPro++){

      // Skip if track was flagged by cleaning procedure
      if (!fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].UseIt())
  	continue;
      // Skip if track was not well propagated to R=1.2(5)m
      if (!fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].GoodPropR12())
  	continue;

      // Calculate pair variables
      x[0]=mt(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
	      fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
      x[1]=dEtaS(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter
		,fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
      x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter
		      ,fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
      x[3]=Qinv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
		fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);

      if(fkDoStdHists){
	// 2d distance at R=1.2
	f2HistALamAProAngDistSft2dAtR12Real->Fill(x[1],x[2]);
	// Mt of the pair
	f2HistMtALamAProReal->Fill(x[0]);
      }

      // Use THn since sept '12, with |dEtaS| instead of dEtaS since 2014
      x[1] = TMath::Abs(x[1]);
      fALamAProReal->Fill(x);

      // Since April 22nd 2013 do the 99% TTR cut also here.
      // x[1] is |dEtaS|, x[2] is dPhiS
      if(TMath::Power(x[1]/fdEtaSCut,2) + TMath::Power(x[2]/fdPhiSCut,2) < 1.){
	// Too close, reject
	continue;
      }

      // Anti-proton properties..
      // ..phase space
      const Float_t rapAPro(fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].RapPro())
	,pTAPro(fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].Pt());
      // .. purity
      const Float_t aproFdPur(GetAProFdPur(rapAPro,pTAPro));

      // Four histograms with pairs as a fct of qinv:
      // three with purity as weight one without
      if(fkDoStdHists){
	f2HistALamAProWALamPurReal->Fill(x[3],alamPur);
	f2HistALamAProWALamFdPurReal->Fill(x[3],alamFdPur);
	f2HistALamAProWAProFdPurReal->Fill(x[3],aproFdPur);
	f2HistALamAProWoPurReal->Fill(x[3]);
      }
      // Also (vs qinv) vs mT
      f2HistALamAProWALamPurmTReal->Fill(x[3],x[0],alamPur);
      f2HistALamAProWALamFdPurmTReal->Fill(x[3],x[0],alamFdPur);
      f2HistALamAProWAProFdPurmTReal->Fill(x[3],x[0],aproFdPur);
      f2HistALamAProWoPurmTReal->Fill(x[3],x[0]);

      // Get purity fluctuations for low q
      if(x[3] <.2){
	f2HistPurFlucALamAPro->Fill(GetALamPur(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].RapLam(),
					      fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].Pt()));
      }

      // Origin in phase space (y,pt) for each pair
      if(fkDoStdHists){
	// For pairs with low q, i.e. qinv < .2 GeV/c
	if(x[3] <.2){
	  // .. for lambdas
	  f2HistALamYPtPair->Fill(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].RapLam(),
				fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].Pt());
	  // .. protons
	  f2HistAProYPtPair->Fill(fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].RapPro(),
				  fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].Pt());
	  // ..pairs
	  f2HistYKtAPair->Fill(yPair(fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro],
				     fFemtoBuffer->GetEvt(0)->fALamTracks[iALam]),
			       TMath::Sqrt(ktSquared(fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro],
						     fFemtoBuffer->GetEvt(0)->fALamTracks[iALam])));
	  // Mt of the pair for low q pairs only
	  f2HistMtLowQALamAProReal->Fill(x[0]);
	  // Also  mT vs phi* (x[2])
	  f2HistMtVsPhiSLowQ->Fill(x[2],x[0]);
	}
	// Origin in phase space (y,pt) for each pair with .4 < qinv < .5 GeV/c
	else if(x[3] < .5 && x[3] >= .4){
	  // just for lambdas
	  f2HistALamYPtHiPair->Fill(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].RapLam(),
				    fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].Pt());
	}
      }	// End of origin in phase space

      
    }// AProton loop
  }// ALambda loop

}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::ProcessMixed() {
  // Process mixed events

  // Declare numbers
  UChar_t iLam,iPro,iALam,iAPro;
  Double_t x[4];
  
  // Loop over the event buffer
  for (UChar_t iMix = 1;iMix<fFemtoBuffer->GetMixBuffSize();iMix++){
    
    // Lambda loop
    for (iLam = 0; iLam < fFemtoBuffer->GetEvt(0)->GetNLam(); iLam++){
      
      // Skip if track was flagged by the cleaning procedure
      if (!fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].UseIt())
  	continue;
      // Skip if track was not well propated
      if (!fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter.GoodPropR12())
  	continue;

      
      // Proton loop
      for (iPro=0;iPro<(fFemtoBuffer->GetEvt(iMix))->GetNPro();iPro++){
	
	// Skip if track was flagged by the cleaning procedure
  	if (!(fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro].UseIt())
  	  continue;
	// Skip if track was not well propated
  	if (!(fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro].GoodPropR12())
  	  continue;

	
	// Calculate pair variables
	x[0]=mt(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
		 fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]);
	x[1]=dEtaS(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter
		  ,fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]);
	x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter
			,fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]);
	x[3]=Qinv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
		  fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]);

	// 2d distance at R=1.2
	if(fkDoStdHists){
	  f2HistLamProAngDistSft2dAtR12Mixed->Fill(x[1],x[2]);
	}

	// Use THn since sept '12, since 2014 use |dEtaS| instead of dEtaS
	x[1] = TMath::Abs(x[1]);
	fLamProMixed->Fill(x);

      }// Proton loop
    }// Lambda loop

    // The other way around, keep the proton from the current event
    // and mix with the lambdas from the buffered event
      
    // Proton loop
    for (iPro=0;iPro<(fFemtoBuffer->GetEvt(0))->GetNPro();iPro++){
      
      // Skip if track was flagged by the cleaning procedure
      if (!(fFemtoBuffer->GetEvt(0))->fProTracks[iPro].UseIt())
	continue;
      // Skip if track was not well propated
      if (!(fFemtoBuffer->GetEvt(0))->fProTracks[iPro].GoodPropR12())
	continue;
      
      // Lambda loop
      for (iLam = 0; iLam < fFemtoBuffer->GetEvt(iMix)->GetNLam(); iLam++){
	
	// Skip if track was flagged by the cleaning procedure
	if (!fFemtoBuffer->GetEvt(iMix)->fLamTracks[iLam].UseIt())
	  continue;
	// Skip if track was not well propated
	if (!fFemtoBuffer->GetEvt(iMix)->fLamTracks[iLam].fPosDaughter.GoodPropR12())
	  continue;


	// Calculate pair variables
	x[0]=mt(fFemtoBuffer->GetEvt(iMix)->fLamTracks[iLam],
		fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
	x[1]=dEtaS(fFemtoBuffer->GetEvt(iMix)->fLamTracks[iLam].fPosDaughter
		  ,fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
	x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(iMix)->fLamTracks[iLam].fPosDaughter
			,fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
	x[3]=Qinv(fFemtoBuffer->GetEvt(iMix)->fLamTracks[iLam],
		  fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);

	// 2d distance at R=1.2
	if(fkDoStdHists){
	  f2HistLamProAngDistSft2dAtR12Mixed->Fill(x[1],x[2]);
	}

	// Use THn since sept '12, since 2014 use |dEtaS| instead of dEtaS
	x[1] = TMath::Abs(x[1]);
	fLamProMixed->Fill(x);

      }// Proton loop
    }// Lambda loop
    
    
    // Anti-lambda loop
    for (iALam = 0; iALam < fFemtoBuffer->GetEvt(0)->GetNALam(); iALam++){
      
      // Skip if track was flagged by cleaning procedure
      if (!fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].UseIt())
  	continue;
      // Skip if track was not well propagated
      if (!fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter.GoodPropR12())
  	continue;

      // AProton loop
      for (iAPro=0;iAPro<(fFemtoBuffer->GetEvt(iMix))->GetNAPro();iAPro++){
	
  	// Skip if track was flagged by cleaning procedure
  	if (!(fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro].UseIt())
  	  continue;
  	// Skip if track was not well propagated
  	if (!(fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro].GoodPropR12())
  	  continue;

	// Calculate pair variables
	x[0]=mt(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
		fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]);
	x[1]=dEtaS(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter
		  ,fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]);
	x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter
			,fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]);
	x[3]=Qinv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
		  fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]);
	// 2d distance at R=1.2
	if(fkDoStdHists){
	  f2HistALamAProAngDistSft2dAtR12Mixed->Fill(x[1],x[2]);
	}

	// Use THn since sept '12, since 2014 use |dEtaS| instead of dEtaS
	x[1] = TMath::Abs(x[1]);
	fALamAProMixed->Fill(x);
	
      }// AProton loop
    }// ALambda loop


    // The other way around: Keep the anti-proton from the current event
    // and mix it with lambdas from the buffered events

    // AProton loop
    for (iAPro=0;iAPro<(fFemtoBuffer->GetEvt(0))->GetNAPro();iAPro++){
	
      // Skip if track was flagged by cleaning procedure
      if (!(fFemtoBuffer->GetEvt(0))->fAProTracks[iAPro].UseIt())
	continue;
      // Skip if track was not well propagated
      if (!(fFemtoBuffer->GetEvt(0))->fAProTracks[iAPro].GoodPropR12())
	continue;
      
      // Anti-lambda loop
      for (iALam = 0; iALam < fFemtoBuffer->GetEvt(iMix)->GetNALam(); iALam++){
	
	// Skip if track was flagged by cleaning procedure
	if (!fFemtoBuffer->GetEvt(iMix)->fALamTracks[iALam].UseIt())
	  continue;
	// Skip if track was not well propagated
	if (!fFemtoBuffer->GetEvt(iMix)->fALamTracks[iALam].fNegDaughter.GoodPropR12())
	  continue;

	// Calculate pair variables
	x[0]=mt(fFemtoBuffer->GetEvt(iMix)->fALamTracks[iALam],
		fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
	x[1]=dEtaS(fFemtoBuffer->GetEvt(iMix)->fALamTracks[iALam].fNegDaughter
		  ,fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
	x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(iMix)->fALamTracks[iALam].fNegDaughter
			,fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
	x[3]=Qinv(fFemtoBuffer->GetEvt(iMix)->fALamTracks[iALam],
		  fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
	// 2d distance at R=1.2
	if(fkDoStdHists){
	  f2HistALamAProAngDistSft2dAtR12Mixed->Fill(x[1],x[2]);
	}

	// Use THn since sept '12, since 2014 use |dEtaS| instead of dEtaS
	x[1] = TMath::Abs(x[1]);
	fALamAProMixed->Fill(x);
	
      }// AProton loop
    }// ALambda loop


  }// Event buffer loop

}// End of void ProcessMixed 
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::ProcessRealBackground() {
  // Process real events with background lambdas
  
  // Declare numbers
  UChar_t iBgLam,iPro,iBgALam,iAPro;
  Double_t x[4];

  // BgLambda loop
  for (iBgLam = 0; iBgLam < fFemtoBuffer->GetEvt(0)->GetNBgLam(); iBgLam++){

    // Skip if flagged by cleaning procedure
    if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].UseIt())
      continue;
    // Skip if not well propagated
    if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter.GoodPropR12())
      continue;
    
    // Proton loop
    for (iPro=0;iPro<fFemtoBuffer->GetEvt(0)->GetNPro();iPro++){

      // Skip if flagged by cleaning procedure
      if (!fFemtoBuffer->GetEvt(0)->fProTracks[iPro].UseIt())
  	continue;
      // Skip if not well propagated
      if (!fFemtoBuffer->GetEvt(0)->fProTracks[iPro].GoodPropR12())
  	continue;

      // Calculate pair variables
      x[0]=mt(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam],
	      fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
      x[1]=dEtaS(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter
		,fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
      x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter
		      ,fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
      x[3]=Qinv(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam],
		fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
      // 2d distance at R=1.2
      if(fkDoStdHists){
	f2HistBgLamProAngDistSft2dAtR12Real->Fill(x[1],x[2]);
      }

      // Since sept '12 do THn, since 2014 use |dEtaS| instead of dEtaS
      x[1] = TMath::Abs(x[1]);
      fBgLamProReal->Fill(x);

    }// Proton loop
  }// BgLambda loop


  // Anti-lambda loop
  for (iBgALam = 0; iBgALam < fFemtoBuffer->GetEvt(0)->GetNBgALam(); iBgALam++){

    // Skip if flagged by cleaning procedure
    if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].UseIt())
      continue;
    // Skip if not well propagated
    if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter.GoodPropR12())
      continue;
    
    // AProton loop
    for (iAPro=0;iAPro<fFemtoBuffer->GetEvt(0)->GetNAPro();iAPro++){

      // Skip if flagged by cleaning procedure
      if (!fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].UseIt())
  	continue;
      // Skip if not well propagated
      if (!fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].GoodPropR12())
  	continue;
      
      // Calculate pair variables
      x[0]=mt(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam],
	      fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
      x[1]=dEtaS(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter
		,fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
      x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter
		      ,fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
      x[3]=Qinv(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam], 
		fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
      
      // 2d distance at R=1.2
      if(fkDoStdHists){
	f2HistBgALamAProAngDistSft2dAtR12Real->Fill(x[1],x[2]);
      }

      // Since sept '12 do THn, since 2014 use |dEtaS| instead of dEtaS
      x[1] = TMath::Abs(x[1]);
      fBgALamAProReal->Fill(x);

    }// AProton loop
  }// BgALambda loop
} // End of void ProcessRealBackground
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::ProcessMixedBackground() {
  // Process mixed events

  // Declare numbers
  UChar_t iBgLam,iPro,iBgALam,iAPro;
  Double_t x[4];
  
  // Loop over the event buffer
  for (UChar_t iMix = 1;iMix<fFemtoBuffer->GetMixBuffSize();iMix++){

    // BgLambda loop
    for (iBgLam = 0; iBgLam < fFemtoBuffer->GetEvt(0)->GetNBgLam(); iBgLam++){

      // Skip if flagged by cleaning procedure
      if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].UseIt())
  	continue;
      // Skip if not well propagated
      if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter.GoodPropR12())
  	continue;
      
      // Proton loop
      for (iPro=0;iPro<(fFemtoBuffer->GetEvt(iMix))->GetNPro();iPro++){
	
  	// Skip if flagged by cleaning procedure
  	if (!(fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro].UseIt())
  	  continue;
  	// Skip if not well propagated
  	if (!(fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro].GoodPropR12())
  	  continue;
	
	// Calculate pair variables
	x[0]=mt(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam],
		fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]);
	x[1]=dEtaS(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter
		  ,fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]);
	x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter
			,fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]);
	x[3]=Qinv(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam],
		  fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]);
						   
	// 2d distance at r=1.2
	if(fkDoStdHists){
	  f2HistBgLamProAngDistSft2dAtR12Mixed->Fill(x[1],x[2]);
	}

	// Since sept '12 do THn, since 2014 use |dEtaS| instead of dEtaS
	x[1] = TMath::Abs(x[1]);
	fBgLamProMixed->Fill(x);

      }// Proton loop
    }// BgLambda loop


    // The other way around, keep the proton from the current 
    // event and mix with bg lambdas from the buffered events
      
    // Proton loop
    for (iPro=0;iPro<(fFemtoBuffer->GetEvt(0))->GetNPro();iPro++){
      // Skip if flagged by cleaning procedure
      if (!(fFemtoBuffer->GetEvt(0))->fProTracks[iPro].UseIt())
	continue;
      // Skip if not well propagated
      if (!(fFemtoBuffer->GetEvt(0))->fProTracks[iPro].GoodPropR12())
	continue;

      // BgLambda loop
      for (iBgLam = 0; iBgLam < fFemtoBuffer->GetEvt(iMix)->GetNBgLam(); iBgLam++){
	// Skip if flagged by cleaning procedure
	if (!fFemtoBuffer->GetEvt(iMix)->fBgLamTracks[iBgLam].UseIt())
	  continue;
	// Skip if not well propagated
	if (!fFemtoBuffer->GetEvt(iMix)->fBgLamTracks[iBgLam].fPosDaughter.GoodPropR12())
	  continue;
	
	// Calculate pair variables
	x[0]=mt(fFemtoBuffer->GetEvt(iMix)->fBgLamTracks[iBgLam],
		fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
	x[1]=dEtaS(fFemtoBuffer->GetEvt(iMix)->fBgLamTracks[iBgLam].fPosDaughter
		  ,fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
	x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(iMix)->fBgLamTracks[iBgLam].fPosDaughter
			,fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
	x[3]=Qinv(fFemtoBuffer->GetEvt(iMix)->fBgLamTracks[iBgLam],
		  fFemtoBuffer->GetEvt(0)->fProTracks[iPro]);
						   
	// 2d distance at r=1.2
	if(fkDoStdHists){
	  f2HistBgLamProAngDistSft2dAtR12Mixed->Fill(x[1],x[2]);
	}

	// Since sept '12 do THn, since 2014 use |dEtaS| instead of dEtaS
	x[1] = TMath::Abs(x[1]);
	fBgLamProMixed->Fill(x);

      }// Proton loop
    }// BgLambda loop
    
    
    // Anti-lambda loop
    for (iBgALam = 0; iBgALam < fFemtoBuffer->GetEvt(0)->GetNBgALam(); iBgALam++){
      // Skip if flagged by the cleaning procedure
      if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].UseIt())
  	continue;
      // Skip if not well propagated
      if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter.GoodPropR12())
  	continue;

      // AProton loop
      for (iAPro=0;iAPro<(fFemtoBuffer->GetEvt(iMix))->GetNAPro();iAPro++){
  	// Skip if flagged by the cleaning procedure
  	if (!(fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro].UseIt())
  	  continue;
  	// Skip if not well propagated
  	if (!(fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro].GoodPropR12())
  	  continue;
	
	// Calculate pair variables
	x[0]=mt(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam],
		fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]);
	x[1]=dEtaS(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter
		  ,fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]);
	x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter
			,fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]);
	x[3]=Qinv(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam], 
		  fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]);

	// 2d distance at r=1.2
	if(fkDoStdHists){
	  f2HistBgALamAProAngDistSft2dAtR12Mixed->Fill(x[1],x[2]);
	}

	// Use THn since Sept '12, since 2014 use |dEtaS| instead of dEtaS
	x[1] = TMath::Abs(x[1]);
	fBgALamAProMixed->Fill(x);
	
	}// AProton loop
    }// BgALambda loop


    // The other way around: Keep the proton from the current
    // event and mix with a bg lambda from a buffered event

    // AProton loop
    for (iAPro=0;iAPro<(fFemtoBuffer->GetEvt(0))->GetNAPro();iAPro++){
      // Skip if flagged by the cleaning procedure
      if (!(fFemtoBuffer->GetEvt(0))->fAProTracks[iAPro].UseIt())
	continue;
      // Skip if not well propagated
      if (!(fFemtoBuffer->GetEvt(0))->fAProTracks[iAPro].GoodPropR12())
	continue;

      // Anti-lambda loop
      for (iBgALam = 0; iBgALam < fFemtoBuffer->GetEvt(iMix)->GetNBgALam(); iBgALam++){
	// Skip if flagged by the cleaning procedure
	if (!fFemtoBuffer->GetEvt(iMix)->fBgALamTracks[iBgALam].UseIt())
	  continue;
	// Skip if not well propagated
	if (!fFemtoBuffer->GetEvt(iMix)->fBgALamTracks[iBgALam].fNegDaughter.GoodPropR12())
	  continue;
	
	// Calculate pair variables
	x[0]=mt(fFemtoBuffer->GetEvt(iMix)->fBgALamTracks[iBgALam],
		fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
	x[1]=dEtaS(fFemtoBuffer->GetEvt(iMix)->fBgALamTracks[iBgALam].fNegDaughter
		  ,fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
	x[2]=dPhiSAtR12(fFemtoBuffer->GetEvt(iMix)->fBgALamTracks[iBgALam].fNegDaughter
			,fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);
	x[3]=Qinv(fFemtoBuffer->GetEvt(iMix)->fBgALamTracks[iBgALam], 
		  fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]);

	// 2d distance at r=1.2
	if(fkDoStdHists){
	  f2HistBgALamAProAngDistSft2dAtR12Mixed->Fill(x[1],x[2]);
	}

	// Use THn since Sept '12, since 2014 use |dEtaS| instead of dEtaS
	x[1] = TMath::Abs(x[1]);
	fBgALamAProMixed->Fill(x);
	
	}// AProton loop
    }// BgALambda loop

  }// Event buffer loop

}// End of void ProcessMixedBackground
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::dEtaS(const FemtoBufferTrack &track1,
					   const FemtoBufferTrack &track2){
  // Returns the pseudorapidity star difference

  // It is important to keep the calculations easy and separated.
  // The calculation of EtaS is straight forward, one just has to
  // do it step by step to not get confused.
  return track1.EtaS() - track2.EtaS();
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::dThetaS(const FemtoBufferTrack &track1,
					       const FemtoBufferTrack &track2){
  // Angular distance with shifted vertex

  // It is important to keep the calculations easy and separated.
  // The calculation of ThetaS is straight forward, one just has to
  // do it step by step to not get confused.

  // Contrary to phi, theta -90 degree is not at all close to theta 90 degree
  return track1.ThetaS() - track2.ThetaS();
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::dPhiSAtR12(const FemtoBufferTrack &track1,
						  const FemtoBufferTrack &track2){
  // returns delta phi star at R=1.2m
  // position at R=1.2m is stored as second radius
  // const Float_t distSft= TMath::Sqrt(TMath::Power(track1.fXshifted[2][0] - track2.fXshifted[2][0],2)
  // 				     +TMath::Power(track1.fXshifted[2][1] - track2.fXshifted[2][1],2));
  const Float_t distSft= TMath::Sqrt(TMath::Power(track1.fXSftR125[0] - track2.fXSftR125[0],2)
				     +TMath::Power(track1.fXSftR125[1] - track2.fXSftR125[1],2));
  return 2.0 * TMath::ATan(distSft/2./(125.));
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::Qinv(const FemtoBufferV0 &v01, const FemtoBufferV0 &v02){
  // Copied from NA49. See http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Qinv
  // Always using lambda mass (no mass difference found yet for lam <-> alam (see PDG))

  // printf("v01 px %3.2f py %3.2f pz %3.2f"
  // 	 "v02 px %3.2f py %3.2f pz %3.2f"
  // 	 "\n"
  // 	 ,v01.fP[0],v01.fP[1],v01.fP[2]
  // 	 ,v02.fP[0],v02.fP[1],v02.fP[2]
  // 	 );

  //Double_t e1 = t1->GetE(mPart1);
  const Double_t e1 = TMath::Sqrt(fkLamMass*fkLamMass + v01.fP[0]*v01.fP[0]+v01.fP[1]*v01.fP[1]+v01.fP[2]*v01.fP[2]);
  //Double_t e2 = t2->GetE(mPart2);
  const Double_t e2 = TMath::Sqrt(fkLamMass*fkLamMass + v02.fP[0]*v02.fP[0]+v02.fP[1]*v02.fP[1]+v02.fP[2]*v02.fP[2]);
  
  // First calculate -Qinv^2  as usual : 
  //qinvL = (e1-e2) * (e1-e2) - Q(t1,t2) * Q(t1,t2);
  const Double_t qinvL = (e1-e2) * (e1-e2) - ( (v01.fP[0]-v02.fP[0])*(v01.fP[0]-v02.fP[0]) + (v01.fP[1]-v02.fP[1])*(v01.fP[1]-v02.fP[1]) + (v01.fP[2]-v02.fP[2])*(v01.fP[2]-v02.fP[2]) );  

  //Qx(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()-t2->GetPx()); };
  //static Double_t Qy(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPy()-t2->GetPy()); };
  //static Double_t Qz(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPz()-t2->GetPz()); };
  //static Double_t  Q(T49ParticleRoot* t1,T49ParticleRoot* t2)
  //               { return  TMath::Sqrt(Qx(t1,t2)*Qx(t1,t2)+Qy(t1,t2)*Qy(t1,t2)+Qz(t1,t2)*Qz(t1,t2)); };
  


  //static Double_t Px(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()+t2->GetPx()); };
  //qP    = (e1-e2)   * (e1+e2)
  //      - Qx(t1,t2) * Px(t1,t2)
  //      - Qy(t1,t2) * Py(t1,t2)
  //      - Qz(t1,t2) * Pz(t1,t2);
  const Double_t qP = (e1-e2)               * (e1+e2)
                    - (v01.fP[0]-v02.fP[0]) * (v01.fP[0]+v02.fP[0])
  	            - (v01.fP[1]-v02.fP[1]) * (v01.fP[1]+v02.fP[1])
                    - (v01.fP[2]-v02.fP[2]) * (v01.fP[2]+v02.fP[2]);

  //pinv  = (e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2);
  const Double_t pinv  = (e1+e2) * (e1+e2) - ( (v01.fP[0]+v02.fP[0])*(v01.fP[0]+v02.fP[0])
					      +(v01.fP[1]+v02.fP[1])*(v01.fP[1]+v02.fP[1])
					      +(v01.fP[2]+v02.fP[2])*(v01.fP[2]+v02.fP[2]));

  return TMath::Sqrt(qP*qP/pinv - qinvL);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::Qinv(const FemtoBufferV0 &v0, const FemtoBufferTrack &track) {
  // Copied from NA49. See http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Qinv
  //  Always using lambda mass (no mass difference found yet for lam <-> alam (see PDG))
  
  //  Double_t e1 = t1->GetE(mPart1);
  const Double_t e1 = TMath::Sqrt(fkLamMass*fkLamMass + v0.fP[0]*v0.fP[0]+v0.fP[1]*v0.fP[1]+v0.fP[2]*v0.fP[2]);
  //  Double_t e2 = t2->GetE(mPart2);
  const Double_t e2 = TMath::Sqrt(fkProMass*fkProMass + track.fP[0]*track.fP[0]+track.fP[1]*track.fP[1]+track.fP[2]*track.fP[2]);

  // First calculate -Qinv^2  as usual : 
  //qinvL = (e1-e2) * (e1-e2) - Q(t1,t2) * Q(t1,t2);
  const Double_t qinvL = (e1-e2) * (e1-e2) - ( (v0.fP[0]-track.fP[0])*(v0.fP[0]-track.fP[0]) + (v0.fP[1]-track.fP[1])*(v0.fP[1]-track.fP[1]) + (v0.fP[2]-track.fP[2])*(v0.fP[2]-track.fP[2]) );  

  //Qx(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()-t2->GetPx()); };
  //static Double_t Qy(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPy()-t2->GetPy()); };
  //static Double_t Qz(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPz()-t2->GetPz()); };
  //static Double_t  Q(T49ParticleRoot* t1,T49ParticleRoot* t2)
  //               { return  TMath::Sqrt(Qx(t1,t2)*Qx(t1,t2)+Qy(t1,t2)*Qy(t1,t2)+Qz(t1,t2)*Qz(t1,t2)); };
  


  //static Double_t Px(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()+t2->GetPx()); };
  //qP    = (e1-e2)   * (e1+e2)
  //      - Qx(t1,t2) * Px(t1,t2)
  //      - Qy(t1,t2) * Py(t1,t2)
  //      - Qz(t1,t2) * Pz(t1,t2);
  const Double_t qP = (e1-e2)                * (e1+e2)
                    - (v0.fP[0]-track.fP[0]) * (v0.fP[0]+track.fP[0])
  	            - (v0.fP[1]-track.fP[1]) * (v0.fP[1]+track.fP[1])
  	            - (v0.fP[2]-track.fP[2]) * (v0.fP[2]+track.fP[2]);

  //pinv  = (e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2);
  const Double_t pinv  = (e1+e2) * (e1+e2) - ( (v0.fP[0]+track.fP[0])*(v0.fP[0]+track.fP[0])
					      +(v0.fP[1]+track.fP[1])*(v0.fP[1]+track.fP[1])
					      +(v0.fP[2]+track.fP[2])*(v0.fP[2]+track.fP[2]));

  return TMath::Sqrt(qP*qP/pinv - qinvL);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::Qinv(const FemtoBufferTrack &track, const FemtoBufferV0 &v0){
  return Qinv(v0, track);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::QinvProPro(const FemtoBufferTrack &proTrack1, const FemtoBufferTrack &proTrack2) {
  // Same as above, with different masses for the tracks,
  // here both tracks are protons

  // Copied from NA49. See http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Qinv
  
  //  PDG_t e1 = t1->GetE(mPart1);
  const Double_t e1 = TMath::Sqrt(fkProMass*fkProMass + proTrack1.fP[0]*proTrack1.fP[0]+proTrack1.fP[1]*proTrack1.fP[1]+proTrack1.fP[2]*proTrack1.fP[2]);
  //  Double_t e2 = t2->GetE(mPart2);
  const Double_t e2 = TMath::Sqrt(fkProMass*fkProMass + proTrack2.fP[0]*proTrack2.fP[0]+proTrack2.fP[1]*proTrack2.fP[1]+proTrack2.fP[2]*proTrack2.fP[2]);
  
  // First calculate -Qinv^2  as usual : 
  //qinvL = (e1-e2) * (e1-e2) - Q(t1,t2) * Q(t1,t2);
  const Double_t qinvL = (e1-e2) * (e1-e2) - ( (proTrack1.fP[0]-proTrack2.fP[0])*(proTrack1.fP[0]-proTrack2.fP[0]) + (proTrack1.fP[1]-proTrack2.fP[1])*(proTrack1.fP[1]-proTrack2.fP[1]) + (proTrack1.fP[2]-proTrack2.fP[2])*(proTrack1.fP[2]-proTrack2.fP[2]) );  

  //Qx(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()-t2->GetPx()); };
  //static Double_t Qy(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPy()-t2->GetPy()); };
  //static Double_t Qz(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPz()-t2->GetPz()); };
  //static Double_t  Q(T49ParticleRoot* t1,T49ParticleRoot* t2)
  //               { return  TMath::Sqrt(Qx(t1,t2)*Qx(t1,t2)+Qy(t1,t2)*Qy(t1,t2)+Qz(t1,t2)*Qz(t1,t2)); };
  


  //static Double_t Px(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()+t2->GetPx()); };
  //qP    = (e1-e2)   * (e1+e2)
  //      - Qx(t1,t2) * Px(t1,t2)
  //      - Qy(t1,t2) * Py(t1,t2)
  //      - Qz(t1,t2) * Pz(t1,t2);
  const Double_t qP = (e1-e2)                           * (e1+e2)
                    - (proTrack1.fP[0]-proTrack2.fP[0]) * (proTrack1.fP[0]+proTrack2.fP[0])
                    - (proTrack1.fP[1]-proTrack2.fP[1]) * (proTrack1.fP[1]+proTrack2.fP[1])
                    - (proTrack1.fP[2]-proTrack2.fP[2]) * (proTrack1.fP[2]+proTrack2.fP[2]);

  //pinv  = (e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2);
  const Double_t pinv = (e1+e2) * (e1+e2) - ( (proTrack1.fP[0]+proTrack2.fP[0])*(proTrack1.fP[0]+proTrack2.fP[0])
					     +(proTrack1.fP[1]+proTrack2.fP[1])*(proTrack1.fP[1]+proTrack2.fP[1])
					     +(proTrack1.fP[2]+proTrack2.fP[2])*(proTrack1.fP[2]+proTrack2.fP[2]));

  return TMath::Sqrt(qP*qP/pinv - qinvL);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::QinvPioPro(const FemtoBufferTrack &pioTrack, const FemtoBufferTrack &proTrack) {
  // Same as above, with different masses for the tracks,
  // here both tracks are protons

  // Copied from NA49. See http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Qinv
  
  //  PDG_t e1 = t1->GetE(mPart1);
  const Double_t e1 = TMath::Sqrt(fkPioMass*fkPioMass + pioTrack.fP[0]*pioTrack.fP[0]+pioTrack.fP[1]*pioTrack.fP[1]+pioTrack.fP[2]*pioTrack.fP[2]);
  //  Double_t e2 = t2->GetE(mPart2);
  const Double_t e2 = TMath::Sqrt(fkProMass*fkProMass + proTrack.fP[0]*proTrack.fP[0]+proTrack.fP[1]*proTrack.fP[1]+proTrack.fP[2]*proTrack.fP[2]);
  
  // First calculate -Qinv^2  as usual : 
  //qinvL = (e1-e2) * (e1-e2) - Q(t1,t2) * Q(t1,t2);
  const Double_t qinvL = (e1-e2) * (e1-e2) - ( (pioTrack.fP[0]-proTrack.fP[0])*(pioTrack.fP[0]-proTrack.fP[0]) + (pioTrack.fP[1]-proTrack.fP[1])*(pioTrack.fP[1]-proTrack.fP[1]) + (pioTrack.fP[2]-proTrack.fP[2])*(pioTrack.fP[2]-proTrack.fP[2]) );  

  //Qx(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()-t2->GetPx()); };
  //static Double_t Qy(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPy()-t2->GetPy()); };
  //static Double_t Qz(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPz()-t2->GetPz()); };
  //static Double_t  Q(T49ParticleRoot* t1,T49ParticleRoot* t2)
  //               { return  TMath::Sqrt(Qx(t1,t2)*Qx(t1,t2)+Qy(t1,t2)*Qy(t1,t2)+Qz(t1,t2)*Qz(t1,t2)); };
  


  //static Double_t Px(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()+t2->GetPx()); };
  //qP    = (e1-e2)   * (e1+e2)
  //      - Qx(t1,t2) * Px(t1,t2)
  //      - Qy(t1,t2) * Py(t1,t2)
  //      - Qz(t1,t2) * Pz(t1,t2);
  const Double_t qP = (e1-e2)                         * (e1+e2)
                    - (pioTrack.fP[0]-proTrack.fP[0]) * (pioTrack.fP[0]+proTrack.fP[0])
                    - (pioTrack.fP[1]-proTrack.fP[1]) * (pioTrack.fP[1]+proTrack.fP[1])
                    - (pioTrack.fP[2]-proTrack.fP[2]) * (pioTrack.fP[2]+proTrack.fP[2]);

  //pinv  = (e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2);
  const Double_t pinv  = (e1+e2) * (e1+e2) - ( (pioTrack.fP[0]+proTrack.fP[0])*(pioTrack.fP[0]+proTrack.fP[0])
					      +(pioTrack.fP[1]+proTrack.fP[1])*(pioTrack.fP[1]+proTrack.fP[1])
					      +(pioTrack.fP[2]+proTrack.fP[2])*(pioTrack.fP[2]+proTrack.fP[2]));

  return TMath::Sqrt(qP*qP/pinv - qinvL);
}
//________________________________________________________________________
// Float_t AliAnalysisTaskProtonLambda2d::QinvConstr(const FemtoBufferV0 &v0, const FemtoBufferTrack &track) {
//   // Same as Qinv(v0,track) but with constrained momentum for the track

//   // Check whether constrained momentum is there
//   if ((track.fPconstr[0]<0.00001)&&(track.fPconstr[1]<0.00001)&&(track.fPconstr[2]<0.00001))
//     return Qinv(v0,track);

//   // Standard Qinv(v0, track), just with constrained momentum instead of TPC only momentum
//   Double_t e1 = TMath::Sqrt(fkLamMass*fkLamMass + v0.fP[0]*v0.fP[0]+v0.fP[1]*v0.fP[1]+v0.fP[2]*v0.fP[2]);
//   Double_t e2 = TMath::Sqrt(fkProMass*fkProMass + track.fPconstr[0]*track.fPconstr[0]+track.fPconstr[1]*track.fPconstr[1]+track.fPconstr[2]*track.fPconstr[2]);
//   Double_t qinvL;
//   Double_t qP;
//   Double_t pinv;
//   qinvL = (e1-e2) * (e1-e2) - ( (v0.fP[0]-track.fPconstr[0])*(v0.fP[0]-track.fPconstr[0]) + (v0.fP[1]-track.fPconstr[1])*(v0.fP[1]-track.fPconstr[1]) + (v0.fP[2]-track.fPconstr[2])*(v0.fP[2]-track.fPconstr[2]) );  
//   qP    = (e1-e2)   * (e1+e2)
//         - (v0.fP[0]-track.fPconstr[0]) * (v0.fP[0]+track.fPconstr[0])
//   	- (v0.fP[1]-track.fPconstr[1]) * (v0.fP[1]+track.fPconstr[1])
//   	- (v0.fP[2]-track.fPconstr[2]) * (v0.fP[2]+track.fPconstr[2]);
//   pinv  = (e1+e2) * (e1+e2) - ( (v0.fP[0]+track.fPconstr[0])*(v0.fP[0]+track.fPconstr[0])
//   			       +(v0.fP[1]+track.fPconstr[1])*(v0.fP[1]+track.fPconstr[1])
//   			       +(v0.fP[2]+track.fPconstr[2])*(v0.fP[2]+track.fPconstr[2]));

//   return TMath::Sqrt(qP*qP/pinv - qinvL);
// }
// //________________________________________________________________________
// Float_t AliAnalysisTaskProtonLambda2d::QinvConstr(const FemtoBufferTrack &track, const FemtoBufferV0 &v0){
//   return QinvConstr(v0, track);
// }
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::Minv(const FemtoBufferV0 &v01, const FemtoBufferV0 &v02){
  // Taken from NA49. See 
  // http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Minv
  
  //  Double_t e1 = t1->GetE(mPart1);
  //  Double_t e2 = t2->GetE(mPart2);
  //  GetE(Float_t mass)  { return sqrt(GetP()*GetP()+mass*mass); }  
  const Float_t e1 = TMath::Sqrt(v01.fP[0]*v01.fP[0] + v01.fP[1]*v01.fP[1] + v01.fP[2]*v01.fP[2] 
				 + fkLamMass*fkLamMass);
  const Float_t e2 = TMath::Sqrt(v02.fP[0]*v02.fP[0] + v02.fP[1]*v02.fP[1] + v02.fP[2]*v02.fP[2] 
				 + fkLamMass*fkLamMass);
  
  // return TMath::Sqrt((e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2));
  return TMath::Sqrt((e1+e2) * (e1+e2) - (  (v01.fP[0]+v02.fP[0])*(v01.fP[0]+v02.fP[0])
  					   +(v01.fP[1]+v02.fP[1])*(v01.fP[1]+v02.fP[1])
  					   +(v01.fP[2]+v02.fP[2])*(v01.fP[2]+v02.fP[2])));

}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::Minv(const FemtoBufferV0 &v0, const FemtoBufferTrack &track){
  // Taken from NA49. See 
  // http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Minv
  
  //  Double_t e1 = t1->GetE(mPart1);
  //  Double_t e2 = t2->GetE(mPart2);
  //  GetE(Float_t mass)  { return sqrt(GetP()*GetP()+mass*mass); }  
  const Float_t e1 = TMath::Sqrt(v0.fP[0]*v0.fP[0] + v0.fP[1]*v0.fP[1] + v0.fP[2]*v0.fP[2] 
  			   + fkLamMass*fkLamMass);
  const Float_t e2 = TMath::Sqrt(track.fP[0]*track.fP[0] + track.fP[1]*track.fP[1] + track.fP[2]*track.fP[2] 
  			   + fkProMass*fkProMass);
  
  // return TMath::Sqrt((e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2));
  return TMath::Sqrt((e1+e2) * (e1+e2) - (  (v0.fP[0]+track.fP[0])*(v0.fP[0]+track.fP[0])
  					   +(v0.fP[1]+track.fP[1])*(v0.fP[1]+track.fP[1])
  					   +(v0.fP[2]+track.fP[2])*(v0.fP[2]+track.fP[2])));

}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::Minv(const FemtoBufferTrack &track, const FemtoBufferV0 &v0){
    return Minv(v0, track);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskProtonLambda2d::goodDCA(const AliAODTrack *track) {
  // Get the DCAxy and DCAz. There also exists a TPC only 
  // impact parameter, but this has not enough resolution 
  // to discriminate between primaries, secondaries and material

  // Find the phase space bin
  // Rejects |y| > 1. or pT > 5.0 GeV/c
  UChar_t rapBin(0),ptBin(0);
  if(!getRapAndPtBin(RapidityProton(track),track->Pt(),rapBin,ptBin))
    return kFALSE;

  // Get the DCAxy and DCAz global
  Double_t dca[2]={-9999.,-9999.};
  if(!DCAxyz(fGTI[-track->GetID()-1], fAOD,dca))
    return kFALSE;

  // Fill the DCAxy histograms, only if requested
  if(fkDoDCAHists){
    if (track->Charge() > 0){
      fPriHistDCAxyzYPtPro[rapBin][ptBin]->Fill(dca[0],dca[1]);
    }
    else{
      fPriHistDCAxyzYPtAPro[rapBin][ptBin]->Fill(dca[0],dca[1]);
    }
  }

  // Do a cut. See ~/alice/wikigrid/analyzeOutput/Fit2dDCA.C
  // protons have tighter cuts because of their higher
  // contamination

  // Allow for variation of cut with fkDCAFlag (0=kStd,1=kWide,2=kStrict)
  // The cut values as std, wide, strict
  static const Double_t proXY[3] = {.10,.20,.08};   //kStd,kWide,kStrict
  static const Double_t proZ[3]  = {.15,.30,.15};
  static const Double_t aproXY[3]= {.15,.30,.15};
  static const Double_t aproZ[3] = {.20,.40,.15};
  // New, cut adjustable with flag
  if(track->Charge()>0) {
    // Protons
    if (TMath::Abs(dca[0])>proXY[fkDCAFlag] ||
	TMath::Abs(dca[1])>proZ[fkDCAFlag]  ) {
      // Bad, reject
      return kFALSE;
    }
  }
  else{
    // Anti-protons
    if (TMath::Abs(dca[0])>aproXY[fkDCAFlag] ||
	TMath::Abs(dca[1])>aproZ[fkDCAFlag]  ) {
      // Bad, reject
      return kFALSE;
    }
  }

  // old, hard coded cuts
  // if(track->Charge()>0) {
  //   // Protons
  //   // std
  //   if (TMath::Abs(dca[0])>.10 || TMath::Abs(dca[1])>.15) {
  //     // wide
  //   // if (TMath::Abs(dca[0])>.20 || TMath::Abs(dca[1])>.30) {
  //   // strict
  //   // if (TMath::Abs(dca[0])>.08 || TMath::Abs(dca[1])>.15) {
  //     // Bad, reject
  //     return kFALSE;
  //   }
  // }
  // else{
  //   // Anti-protons
  //   // std
  //   if (TMath::Abs(dca[0])>.15 || TMath::Abs(dca[1])>.20) {
  //   // wide
  //   // if (TMath::Abs(dca[0])>.30 || TMath::Abs(dca[1])>.40) {
  //   // strict
  //   // if (TMath::Abs(dca[0])>.15 || TMath::Abs(dca[1])>.15) {
  //     // Bad, reject
  //     return kFALSE;
  //   }
  // }
  
  // Good track
  return kTRUE;
}
//_______________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::RapidityProton(const AliAODTrack *track){
  // Can't find how to set the assumed mass for the AliAODTrack.
  // Same stuff as in AliAODTrack::Y() just with proton mass
  const Double_t e = TMath::Sqrt(track->P()*track->P() + fkProMass*fkProMass);
  const Double_t pz = track->Pz();
  if (e != TMath::Abs(pz)) { // energy was not equal to pz
    return 0.5*TMath::Log((e+pz)/(e-pz));
  } else { // energy was equal to pz
    return -999.;
  }
}
//________________________________________________________________________
Bool_t AliAnalysisTaskProtonLambda2d::DCAxyz(const AliAODTrack *track, 
					     const AliVEvent *evt,
					     Double_t dca[2]){
  // Note that AliAODTrack::PropagateToDCA() changes the track. 
  // Don't know whether this is what one wants?
  if(!track){
    printf("Pointer to track is zero!\n");
    return kFALSE;
  }

  // See https://savannah.cern.ch/bugs/?102721
  // which is one of the two 11h re-filtering follow-ups:
  // Andrea Dainese now first does the beam pipe
  // check and then copies from the vtrack (was the other
  // way around) to avoid the crash in the etp::Set()
  Double_t xyz[3]={0.,0.,0.};
  track->GetPosition(xyz);
  if(xyz[0]*xyz[0]+xyz[1]*xyz[1] > 3.*3.){
    // printf("This method can be used only for propagation inside the beam pipe\n");
    // printf("  id: %d, filtermap: %d\n",track->GetID(),track->GetFilterMap());
    return kFALSE; 
  }

  // Create an external parameter from the AODtrack
  AliExternalTrackParam etp; etp.CopyFromVTrack(track);

  // Do the propagation (xyz is abused as unused covariance)
  if(!etp.PropagateToDCA(evt->GetPrimaryVertex(),evt->GetMagneticField(),10.,dca,xyz))
    return kFALSE;
  
  // We did it, good dca
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskProtonLambda2d::getRapAndPtBin(const Float_t rap,const Float_t pt,
						   UChar_t &rapBin,UChar_t &ptBin){
  // The DCAxyz histograms are in an array
  // Float_t DCAxyzRapBins[]={-1.,-.5,-.25,0.,.25,.5,1.};
  // Float_t DCAxyzPtBins[]={0.,.5,1.,1.5,2.,3.,5.};

  // Get the rapidity bin
  if(rap<-1.0)return kFALSE;
  else if(rap<-.5)rapBin=0;
  else if(rap<-.25)rapBin=1;
  else if(rap<0.)rapBin=2;
  else if(rap<.25)rapBin=3;
  else if(rap<.5)rapBin=4;
  else if(rap<1.0)rapBin=5;
  else return kFALSE;
  // Get the pt bin
  if(pt<0.)return kFALSE;
  else if(pt<.5)ptBin=0;
  else if(pt<1.)ptBin=1;
  else if(pt<1.5)ptBin=2;
  else if(pt<2.)ptBin=3;
  else if(pt<3.)ptBin=4;
  else if(pt<5.)ptBin=5;
  else return kFALSE;
  // Done
  return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FillDedxHist(const AliAODTrack *track){
  // This is for visualization. Fill the the dE/dx histograms
  // for all tracks, not only for those, where only the TPC
  // is used for PID. Thus avoiding the sharp cut off at a 
  // momentum of 0.75 GeV/c.

  // Only do hists if requested
  if(!fkDoStdHists){
    return;
  }

  if(!(fGTI[-track->GetID()-1])){
    printf("Warning: No global track info there!\n");
    return;
  }

  // TPC signal and Nsigma. See STEER/STEERBase/AliPIDResponse.h for how 
  // NSigmaTPC works (and refrain from banging your head against the wall
  // when you see it).
  // Positive tracks
  if (track->Charge() > 0){
    fPriHistTPCsignalPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
			       (fGTI[-track->GetID()-1])->GetTPCsignal());
    // Fill histograms in three momentum ranges
    fPriHistTPCsignalLowPPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
				   (fGTI[-track->GetID()-1])->GetTPCsignal());
    fPriHistTPCsignalMedPPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
				   (fGTI[-track->GetID()-1])->GetTPCsignal());
    fPriHistTPCsignalHigPPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
				   (fGTI[-track->GetID()-1])->GetTPCsignal());  
    
  }
  // Negative tracks
  else{ 
    fPriHistTPCsignalNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
			       (fGTI[-track->GetID()-1])->GetTPCsignal());
    // Fill histograms in three momentum ranges
    fPriHistTPCsignalLowPNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
				   (fGTI[-track->GetID()-1])->GetTPCsignal());
    fPriHistTPCsignalMedPNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
				   (fGTI[-track->GetID()-1])->GetTPCsignal());
    fPriHistTPCsignalHigPNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
				   (fGTI[-track->GetID()-1])->GetTPCsignal());  
  }
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::GetLamPur(Float_t y, Float_t pt){
  // Returns the lambda purity for given y,pt 

  // Using common fct, we only select the hist
  return GetPur(y,pt,&fLamPur);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::GetALamPur(Float_t y, Float_t pt){
  // Returns the anti lambda purity for given y,pt 

  // Using common fct, we only select the hist
  return GetPur(y,pt,&fALamPur);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::GetLamFdPur(Float_t y, Float_t pt,const Float_t cent){
  // Returns the lambda purity for given y,pt and event centrality
  
  // We choose the correct histogram according to the centrality
  if(cent < 10.){return GetPur(y,pt,&fLamFdPur010);}
  else if(cent < 20.){return GetPur(y,pt,&fLamFdPur1020);}
  else if(cent < 40.){return GetPur(y,pt,&fLamFdPur2040);}
  else if(cent < 60.){return GetPur(y,pt,&fLamFdPur4060);}
  
  AliFatal(Form("Out of cent range %2.1f",cent));
  return 0.;
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::GetALamFdPur(Float_t y, Float_t pt,const Float_t cent){
  // Returns the lambda purity for given y,pt and event centrality
  
  // We choose the correct histogram according to the centrality
  if(cent < 10.){return GetPur(y,pt,&fALamFdPur010);}
  else if(cent < 20.){return GetPur(y,pt,&fALamFdPur1020);}
  else if(cent < 40.){return GetPur(y,pt,&fALamFdPur2040);}
  else if(cent < 60.){return GetPur(y,pt,&fALamFdPur4060);}
  
  AliFatal(Form("Out of cent range %2.1f",cent));
  return 0.;
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::GetProFdPur(Float_t y, Float_t pt){
  // Returns the anti lambda purity for given y,pt 

  // Using common fct, we only select the hist
  return GetPur(y,pt,&fProFdPur);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::GetAProFdPur(Float_t y, Float_t pt){
  // Returns the anti lambda purity for given y,pt 

  // Using common fct, we only select the hist
  return GetPur(y,pt,&fAProFdPur);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::GetPur(Float_t y, Float_t pt,const TH2 *h){
  // Just returns the bin content of h at y and pt
  
  // Check we're inside the histogram range, with safety margin
  if(y < h->GetXaxis()->GetBinCenter(1)) 
    y = h->GetXaxis()->GetBinCenter(1);
  else if(y > h->GetXaxis()->GetBinCenter(h->GetXaxis()->GetNbins()))
    y = h->GetXaxis()->GetBinCenter(h->GetXaxis()->GetNbins());
  if(pt < h->GetYaxis()->GetBinCenter(1))
    pt = h->GetYaxis()->GetBinCenter(1);
  else if(pt > h->GetYaxis()->GetBinCenter(h->GetYaxis()->GetNbins()))
    pt = h->GetYaxis()->GetBinCenter(h->GetYaxis()->GetNbins());
  
  // Get the purity
  return h->GetBinContent(h->GetXaxis()->FindBin(y)
  			  ,h->GetYaxis()->FindBin(pt));
  
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::StoreGlobalTrackReference(const AliAODTrack *track){
  // Stores the pointer to the global track

  // This was AOD073
  // // Don't use the filter bits 2 (ITS standalone) and 128 TPC only
  // // Remove this return statement and you'll see they don't have
  // // any TPC signal
  // if(track->TestFilterBit(128) || track->TestFilterBit(2))
  //   return;
  // This is AOD086
  // Another set of tracks was introduced: Global constrained.
  // We only want filter bit 1 <-- NO! we also want no 
  // filter bit at all, which are the v0 tracks
  //  if(!track->TestFilterBit(1))
  //    return;

  // There are also tracks without any filter bit, i.e. filter map 0,
  // at the beginning of the event: they have ~id 1 to 5, 1 to 12
  // This are tracks that didn't survive the primary track filter but
  // got written cause they are V0 daughters

  // Check whether the track has some info
  // I don't know: there are tracks with filter bit 0
  // and no TPC signal. ITS standalone V0 daughters?
  // if(!track->GetTPCsignal()){
  //   printf("Warning: track has no TPC signal, "
  // 	   //	   "not adding it's info! "
  // 	   "ID: %d FilterMap: %d\n"
  // 	   ,track->GetID(),track->GetFilterMap());
  //   //    return;
  // }
  
  // Check that the id is positive
  if(track->GetID()<0){
    //    printf("Warning: track has negative ID: %d\n",track->GetID());
    return;
  }

  // Check id is not too big for buffer
  if(track->GetID()>=fTrackBuffSize){
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
	   ,track->GetID(),fTrackBuffSize);
    return;
  }

  // Warn if we overwrite a track
  if(fGTI[track->GetID()]){
    // Seems like there are FilterMap 0 tracks
    // that have zero TPCNcls, don't store these!
    if( (!track->GetFilterMap()) &&
	(!track->GetTPCNcls())   )
      return;

    // Imagine the other way around,
    // the zero map zero clusters track
    // is stored and the good one wants 
    // to be added. We ommit the warning
    // and just overwrite the 'bad' track
    if( fGTI[track->GetID()]->GetFilterMap() ||
	fGTI[track->GetID()]->GetTPCNcls()   ){
      // If we come here, there's a problem
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
	     (fGTI[track->GetID()])->GetTPCNcls(),track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
	     (fGTI[track->GetID()])->GetFilterMap(),track->GetFilterMap());
    }
  } // Two tracks same id

  // // There are tracks with filter bit 0,
  // // do they have TPCNcls stored?
  // if(!track->GetFilterMap()){
  //   printf("Filter map is zero, TPCNcls: %u\n"
  // 	   ,track->GetTPCNcls());
  // }

  // Assign the pointer
  (fGTI[track->GetID()]) = track;
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::ResetGlobalTrackReference(){
  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for(UShort_t i=0;i<fTrackBuffSize;i++){
    fGTI[i]=0;
  }
}
//________________________________________________________________________
Bool_t AliAnalysisTaskProtonLambda2d::acceptTrack(const AliAODTrack *track){
  // Apply additional track cuts

  // In the documents
  // https://alisoft.cern.ch/AliRoot/trunk/TPC/doc/Definitions/Definitions.pdf
  // TPC people describe the cut strategy for the TPC. It is explicitly
  // stated that a cut on the number of crossed rows and a cut on the
  // number of crossed rows over findable clusters is recommended to 
  // remove fakes. In the pdf a cut value of .83 on the ratio 
  // is stated, no value for the number of crossed rows. Looking at the 
  // AliESDtrackCuts.cxx one sees that exactly this cut is used with
  // 0.8 on the ratio and 70 on the crossed rows.

  // Checked the filter task and AliAODTrack and AliESDtrack and
  // AliESDtrackCuts and the Definitions.pdf:
  // The function to get the findable clusters is GetTPCNclsF()
  
  // For the number fo crossed rows for ESD tracks, the function
  // GetTPCCrossedRows() usually is used. Looking at the AliESDtrack.cxx
  // one sees that it's just an alias (with additional caching) for
  // GetTPCClusterInfo(2, 1); The identical function exists in the
  // AliAODTrack.cxx

  // I checked: for AOD073 both, the number of crossed rows and the
  // number of findable clusters, are there.

  // WARNING: in LHC10h pass2 the cluster map is wrong for 
  // sector 0 / 18. It's used in the calculation of
  // the number of crossed rows!

  const Float_t nCrossed = track->GetTPCClusterInfo(2, 1);
  if(nCrossed<70)
    return kFALSE;
  if(!track->GetTPCNclsF())
    return kFALSE; // Note that the AliESDtrackCuts would here return kTRUE
  if((nCrossed/track->GetTPCNclsF()) < .8)
    return kFALSE;
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskProtonLambda2d::GoodTPCFitMapSharedMap(const AliAODTrack *pTrack,
							     const AliAODTrack *nTrack){
  // Rejects tracks with shared clusters after filling a control histogram
  // This overload is used for positive and negative daughters from V0s

  // Get the shared maps
  const TBits posSharedMap = pTrack->GetTPCSharedMap();
  const TBits negSharedMap = nTrack->GetTPCSharedMap();
  // Reject shared clusters
  if( ((posSharedMap.CountBits()) >= 1) ||
      ((negSharedMap.CountBits()) >= 1)){
    // Bad tracks, have too many shared clusters!
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskProtonLambda2d::GoodTPCFitMapSharedMap(const AliAODTrack *track){
  // Rejects tracks with shared clusters after filling a control histogram
  // This overload is used for primaries

  // Get the shared maps
  const TBits sharedMap = track->GetTPCSharedMap();
  // Fill a control histogram
  if(fkDoStdHists){
    fPriHistShare->Fill(sharedMap.CountBits());
  }
  // Reject shared clusters
  if((sharedMap.CountBits()) >= 1){
    // Bad track, has too many shared clusters!
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::GetCorrectedTOFSignal(const AliAODTrack *track){
  // Return the corrected TOF signal, see https://twiki.cern.ch/twiki/bin/viewauth/ALICE/TOF
  if(!track){
    printf("E-GetCorrectedTOFSignal::Nullpointer!\n");
    return -9999;
  }

  // Get the global track if ID is negative
  const AliAODTrack *globalTrack(0);
  if(track->GetID()<0){
    // TPC only track or similar, get the global track from fGTI,
    // check for boundaries was done before
    globalTrack = fGTI[-track->GetID()-1];
    // Check it worked
    if(!globalTrack){
      printf("E-GetCorrectedTOFSignal::No global track!\n");
      return -9999;
    }
  }
  else{
    // Track is already global track
    globalTrack=track;
  }

  // Request the TOFpid bit
  if(! (globalTrack->GetStatus()&AliAODTrack::kTOFpid))
    return -9999.;

  // Changed to global track on June 28th 2014
  const Float_t protonTime = GetProtonTime(globalTrack);

  // Check for TOF header
  if(fAOD->GetTOFHeader()){
    // New AODs without start time subtraction
    return (GetTOFSignal(globalTrack)
	    - protonTime
	    - fPIDResponse->GetTOFResponse().GetStartTime(globalTrack->P()));
  }

  // Old AODs with start time already subtracted
  return (GetTOFSignal(globalTrack)
	  - protonTime);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::GetTOFSignal(const AliAODTrack *track){
  // Wrapper for AliAODTrack::GetTOFSignal to 
  // make the offending line more obvious
  return track->GetTOFsignal();
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::GetProtonTime(const AliAODTrack *track){
  // Wrapper for AliAODTrack::GetTOFSignal to 
  // make the offending line more obvious

  // The expected time
  Double_t expectedTimes[AliPID::kSPECIES];
  // Should give number of particles
  //track->GetIntegratedTimes(expectedTimes);
  track->GetIntegratedTimes(expectedTimes,AliPID::kSPECIES);
  

  return expectedTimes[AliPID::kProton];
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::mt(const FemtoBufferTrack &track, const FemtoBufferV0 &v0) {
  // Overloaded function
  return mt(v0,track);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::mt(const FemtoBufferV0 &v0, const FemtoBufferTrack &track){
  // Returns the transverse mass of the pair assuming 
  // proton mass for track and lambda mass for v0

  // Following Phys Rev C 83, 054906
  return TMath::Sqrt(ktSquared(v0,track) +
		     TMath::Power((0.5*(fkLamMass + fkProMass)),2));
}

//________________________________________________________________________
 Float_t AliAnalysisTaskProtonLambda2d::mt(const FemtoBufferV0 &v01, const FemtoBufferV0 &v02){
  // Returns the transverse mass of the pair assuming 
  // lambda mass for both v0

  // Following Phys Rev C 83, 054906
  return TMath::Sqrt(ktSquared(v01,v02) +
		     TMath::Power(fkLamMass,2));
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::ktSquared(const FemtoBufferV0 &v01, const FemtoBufferV0 &v02){
  // Returns the kt squared
  // kt = 1/2 * | (vector{pt1} + vector{pt2}) |
  // kt = 1/2 * | ({px1+px2}, {py1+py2}) |
  // kt2 = 1/2*1/2 * ( (px1+px2)*(px1+px2) + (py1+py2)*(py1+py2) )
  return .5*.5*(  (v01.fP[0] + v02.fP[0])*(v01.fP[0] + v02.fP[0])
		+ (v01.fP[1] + v02.fP[1])*(v01.fP[1] + v02.fP[1]));
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::ktSquared(const FemtoBufferTrack &track, const FemtoBufferV0 &v0){
  // Overloaded function
  return ktSquared(v0,track);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::ktSquared(const FemtoBufferV0 &v0, const FemtoBufferTrack &track){
  // Returns the kt squared
  // kt = 1/2 * | (vector{pt1} + vector{pt2}) |
  // kt = 1/2 * | ({px1+px2}, {py1+py2}) |
  // kt2 = 1/2*1/2 * ( (px1+px2)*(px1+px2) + (py1+py2)*(py1+py2) )
  return .5*.5*(  (v0.fP[0] + track.fP[0])*(v0.fP[0] + track.fP[0])
		+ (v0.fP[1] + track.fP[1])*(v0.fP[1] + track.fP[1]));
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::yPair(const FemtoBufferTrack &track, const FemtoBufferV0 &v0){
  // Over loaded
 return yPair(v0, track);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda2d::yPair(const FemtoBufferV0 &v0, const FemtoBufferTrack &track){
  // Returns the pair rapidity (see Phys Rev C83, 054906) 
  // assuming proton and lambda mass
  // y_pL = 1/2 * ln ( (Ep + EL + pzp + pzL) / (Ep + EL - pzp - pzL) )
  
  // Calc Ep and EL
  const Float_t Ep = TMath::Sqrt(track.fP[0]*track.fP[0] +
				 track.fP[1]*track.fP[1] +
				 track.fP[2]*track.fP[2] +
				 0.9382720*0.9382720);
  const Float_t EL =  TMath::Sqrt(v0.fP[0]*v0.fP[0] +
				  v0.fP[1]*v0.fP[1] +
				  v0.fP[2]*v0.fP[2] +
				  1.115683*1.115683);
  
  // Calc y_pair
  return .5* TMath::Log( (Ep + EL + track.fP[2] + v0.fP[2]) /
			 (Ep + EL - track.fP[2] - v0.fP[2]));
}
// //________________________________________________________________________
// AliAnalysisTaskProtonLambda2d::AliAnalysisTaskProtonLambda2d(const AliAnalysisTaskProtonLambda2d& atpl)
//   // Not implemented, only initializing the const data member as the compiler complains.
//   // Implementation is straight forward, though.
//   : AliAnalysisTaskSE(atpl),
//     fkDoStdHists(atpl.fkDoStdHists),
//     fkDoDCAHists(atpl.fkDoDCAHists),
//     fkDoBgLamALam(atpl.fkDoBgLamALam),
//     fkUseOnTheFly(atpl.fkUseOnTheFly),
//     fkAbsZvertexCut(atpl.fkAbsZvertexCut),
//     fkCentCutLo(atpl.fkCentCutLo),
//     fkCentCutHi(atpl.fkCentCutHi),
//     fkMinvCut(atpl.fkMinvCut),
//     fkCentEst(atpl.fkCentEst),
//     fdEtaSCut(atpl.fdEtaSCut),
//     fdPhiSCut(atpl.fdPhiSCut),
//     fkLamMass(atpl.fkLamMass),
//     fkProMass(atpl.fkProMass),
//     fkPioMass(atpl.fkPioMass),
    
//     fPIDResponse(0), 
//     fTpcResponse(0),
//     fFemtoBuffer(0),
//     fAOD(0), fPrimaryVtx(0), fOutputList(0), fOutputPrimaries(0),
//     fOutput2Part(0),
//     fGTI(0),    
//     fTrackBuffSize(atpl.fTrackBuffSize),
//     fLamPur(atpl.fLamPur),
//     fALamPur(atpl.fALamPur),
//     fHistGoodEvent(0),
//     fHistSideBandOffLam(0), fHistSideBandOffALam(0),
//     fHistMassLambdaOff(0), fHistMassAntiLambdaOff(0),        
//     fHistYPtMassLamOff(0), fHistYPtMassALamOff(0),
//     fHistSideBandOnLam(0), fHistSideBandOnALam(0),
//     fHistMassLambdaOn(0),fHistMassAntiLambdaOn(0),
//     fHistYPtMassLamOn(0),fHistYPtMassALamOn(0),
//     fPriHistShare(0),
//     fPriHistTOFsignalPosVsP(0),      
//     fPriHistTOFsignalNegVsP(0),
//     fPriHistHybridTOFsigPosWoTPC(0), 
//     fPriHistHybridTOFsigPosTPCok(0),fPriHistHybridTOFsigNegWoTPC(0),fPriHistHybridTOFsigNegTPCok(0),
//     fPriHistTPCsignalPos(0),
//     fPriHistTPCsignalLowPPos(0),fPriHistTPCsignalMedPPos(0),fPriHistTPCsignalHigPPos(0),
//     fPriHistTPCsignalNeg(0),
//     fPriHistTPCsignalLowPNeg(0),fPriHistTPCsignalMedPNeg(0),fPriHistTPCsignalHigPNeg(0),
//     fPriHistDCAxyzYPtPro(0),fPriHistDCAxyzYPtAPro(0),
//     f2HistMtLamProReal(0), 
//     f2HistMtALamAProReal(0),
//     f2HistMtLowQLamProReal(0), 
//     f2HistMtLowQALamAProReal(0),
//     fLamProReal(0),fALamAProReal(0),
//     f2HistLamProWLamPurReal(0),f2HistLamProWoPurReal(0),
//     fBgLamProReal(0),fBgALamAProReal(0),
//     fLamProMixed(0),fALamAProMixed(0),
//     fBgLamProMixed(0),fBgALamAProMixed(0)
// {
//   // Copy constructor
//   printf("Copy constructor not implemented\n");
// }
// //________________________________________________________________________
// AliAnalysisTaskProtonLambda2d& AliAnalysisTaskProtonLambda2d::operator=(const AliAnalysisTaskProtonLambda2d& atpl)
// {
//   if(this!=&atpl){
//   // One operation with the atpl to get rid of the warning unused parameter
//   fPrimaryVtxPosition[0]=atpl.fPrimaryVtxPosition[0];
//   printf("Assignment operator not implemented\n");
//   }
//   return *this;
// }
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  // Simply memory profiling
#ifdef ROOT_TObjectTable
  printf(" Object table in Terminate() \n");
  printf("------------------------------\n");
  gObjectTable->Print();
  printf("------------------------------\n");
#endif
}
//________________________________________________________________________
//
//
//     Classes in the class AliAnalysisTaskProtonLambda2d
//         FemtoBuffer, FemtoBufferEvent, FemtoBufferV0 and FemtoBufferTrack
//
//________________________________________________________________________
//
//                        FemtoBufferTrack
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::FemtoBufferTrack():
  fID(65535)
  // ,fP{-9999.,-9999.,-9999.}    //C++11
  // ,fXR125{-9999.,-9999.,-9999.} //C++11
{
  // Standard constructor, initialize everything with values indicating 
  // a track that should not be used
  
  // Initialize arrays, old C++03 way
  for (UChar_t i=0;i<3;i++){
    fP[i]=-9999.;
    fXSftR125[i]=-9999.;
    // for (UChar_t j=0;j<9;j++){
    //   fXglobal[j][i]=-9999.;
    //   fXshifted[j][i]=-9999.;
    // }
  }
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::FemtoBufferTrack(const AliAODTrack *track,const Float_t bfield,const Float_t priVtx[3]):
  fID(65535)  
{
  // Constructor

  // Use the function to have the code in one place
  Set(track,bfield,priVtx);
}
// //________________________________________________________________________
// void AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::GetGlobalPositionAtGlobalRadii(const AliAODTrack *track, const Float_t bfield){
//   // // Function not used, do dummy operations to get rid of warnings
//   // Float_t a=bfield;
//   // a=track->P();
  
//   // Gets the global position of the track at nine different radii in the TPC
//   // track is the track you want to propagate
//   // bfield is the magnetic field of your event
//   // globalPositionsAtRadii is the array of global positions in the radii and xyz
  
//   // We have two versions of the two track resolution plots in our proton-lambda task:
//   // a) with all events shifted to (0,0,0), b) without shift.
//   // For a) we should compare the tracks at shifted radii,
//   // for b) we should still use the global radii. This function here is for b).

//   // Initialize the array to something indicating there was no propagation
//   for(Int_t i=0;i<9;i++){
//     for(Int_t j=0;j<3;j++){
//       fXglobal[i][j]=-9999.;
//     }
//   }

//    // Make a copy of the track to not change parameters of the track
//   AliExternalTrackParam etp; etp.CopyFromVTrack(track);
//   //  printf("\nAfter CopyFromVTrack\n");
//   //  etp.Print();
 
//   // The global position of the the track
//   Double_t xyz[3]={-9999.,-9999.,-9999.};  

//   // Counter for which radius we want
//   Int_t iR=0; 
//   // The radii at which we get the global positions
//   // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
//   // Compare squared radii for faster code
//   Float_t RSquaredWanted[9]={85.*85.,105.*105.,125.*125.,145.*145.,165.*165.,
//   			     185.*185.,205.*205.,225.*225.,245.*245.}; 
//   // The global radius we are at, squared. Compare squared radii for faster code
//   Float_t globalRadiusSquared=0;

//   // Propagation is done in local x of the track
//   for (Float_t x = 58.;x<247.;x+=1.){
//     // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
//     // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
//     // the track is straight, i.e. has inifinite pt and doesn't get bent. 
//     // If the track's momentum is smaller than infinite, it will develop a y-component,
//     // which adds to the global radius

//     // Stop if the propagation was not succesful. This can happen for low pt tracks
//     // that don't reach outer radii
//     if(!etp.PropagateTo(x,bfield))break;
//     etp.GetXYZ(xyz); // GetXYZ returns global coordinates

//     // No shifting for global radii
//     globalRadiusSquared = (xyz[0])*(xyz[0])
//                         + (xyz[1])*(xyz[1]);

//     // Roughly reached the radius we want
//     if(globalRadiusSquared > RSquaredWanted[iR]){
      
//       // Bigger loop has bad precision, we're nearly one centimeter too far, 
//       // go back in small steps.
//       while (globalRadiusSquared>RSquaredWanted[iR]){
//   	x-=.1;
//   	//	printf("propagating to x %5.2f\n",x);
//   	if(!etp.PropagateTo(x,bfield))break;
//   	etp.GetXYZ(xyz); // GetXYZ returns global coordinates

//   	// No shifting for global radii
//   	globalRadiusSquared = (xyz[0])*(xyz[0])
//   	                    + (xyz[1])*(xyz[1]);
//       }
//       //      printf("At Radius:%05.2f (local x %5.2f). Setting position to x %4.1f y %4.1f z %4.1f\n",TMath::Sqrt(globalRadiusSquared),x,xyz[0],xyz[1],xyz[2]);
//       fXglobal[iR][0]=xyz[0];
//       fXglobal[iR][1]=xyz[1];
//       fXglobal[iR][2]=xyz[2];
//       // Indicate we want the next radius    
//       iR+=1;
//     }
//     if(iR>=8){
//       // TPC edge reached
//       return;
//     }
//   }
// }
// //________________________________________________________________________
// void AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::GetShiftedPositionAtShiftedRadii(const AliAODTrack *track, const Float_t bfield, const Float_t priVtx[3]){
//   // Gets the global position of the track at nine different radii in the TPC
//   // track is the track you want to propagate
//   // bfield is the magnetic field of your event
//   // globalPositionsAtRadii is the array of global positions in the radii and xyz
  
//   // Initialize the array to something indicating there was no propagation
//   for(Int_t i=0;i<9;i++){
//     for(Int_t j=0;j<3;j++){
//       fXshifted[i][j]=-9999.;
//     }
//   }

//    // Make a copy of the track to not change parameters of the track
//   AliExternalTrackParam etp; etp.CopyFromVTrack(track);
//   //  printf("\nAfter CopyFromVTrack\n");
//   //  etp.Print();
 
//   // The global position of the the track
//   Double_t xyz[3]={-9999.,-9999.,-9999.};  

//   // Counter for which radius we want
//   Int_t iR=0; 
//   // The radii at which we get the global positions
//   // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
//   // Compare squared radii for faster code
//   Float_t RSquaredWanted[9]={85.*85.,105.*105.,125.*125.,145.*145.,165.*165.,
// 			     185.*185.,205.*205.,225.*225.,245.*245.}; 
//   // The shifted radius we are at, squared. Compare squared radii for faster code
//   Float_t shiftedRadiusSquared=0;

//   // Propagation is done in local x of the track
//   for (Float_t x = 58.;x<247.;x+=1.){
//     // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
//     // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
//     // the track is straight, i.e. has inifinite pt and doesn't get bent. 
//     // If the track's momentum is smaller than infinite, it will develop a y-component,
//     // which adds to the global radius

//     // Stop if the propagation was not succesful. This can happen for low pt tracks
//     // that don't reach outer radii
//     if(!etp.PropagateTo(x,bfield))break;
//     etp.GetXYZ(xyz); // GetXYZ returns global coordinates

//     // Without shifting the primary vertex to (0.,0.,0.) the next line would just be
//     // WRONG: globalRadiusSquared = xyz[0]*xyz[0]+xyz[1]*xyz[1];
//     // but as we shift the primary vertex we want to compare positions at shifted radii.
//     // I can't draw in ASCII but please take a piece of paper and just visualize it once.

//     // Changing plus to minus on July10th2012
//     shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
//                          + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

//     // Roughly reached the radius we want
//     if(shiftedRadiusSquared > RSquaredWanted[iR]){
      
//       // Bigger loop has bad precision, we're nearly one centimeter too far, 
//       // go back in small steps.
//       while (shiftedRadiusSquared>RSquaredWanted[iR]){
// 	x-=.1;
// 	//	printf("propagating to x %5.2f\n",x);
// 	if(!etp.PropagateTo(x,bfield))break;
// 	etp.GetXYZ(xyz); // GetXYZ returns global coordinates
// 	// Added the shifting also here on July11th2012
// 	shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
// 	                     + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);
//       }
//       //      printf("At Radius:%05.2f (local x %5.2f). Setting position to x %4.1f y %4.1f z %4.1f\n",TMath::Sqrt(globalRadiusSquared),x,xyz[0],xyz[1],xyz[2]);
//       fXshifted[iR][0]=xyz[0]-priVtx[0];
//       fXshifted[iR][1]=xyz[1]-priVtx[1];
//       fXshifted[iR][2]=xyz[2]-priVtx[2];
//       // Indicate we want the next radius    
//       iR+=1;
//     }
//     if(iR>=8){
//       // TPC edge reached
//       return;
//     }
//   }
// }
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::SetSftPosR125(const AliAODTrack *track
								    ,const Float_t bfield
								    ,const Float_t priVtx[3]){
  // Sets the spatial position of the track at the radius R=1.25m in the shifted coordinate system
  
  // Initialize the array to something indicating there was no propagation
  fXSftR125[0]=-9999.;
  fXSftR125[1]=-9999.;
  fXSftR125[2]=-9999.;
   // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);
  
  // The global position of the the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};  

  // The radius we want to propagate to, squared
  const Float_t RSquaredWanted(125.*125.);


  // Propagation is done in local x of the track
  for (Float_t x = 58.;x<247.;x+=1.){
    // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
    // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
    // the track is straight, i.e. has inifinite pt and doesn't get bent. 
    // If the track's momentum is smaller than infinite, it will develop a y-component,
    // which adds to the global radius
    // We don't change the propagation steps to not mess up things!

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
    if(!etp.PropagateTo(x,bfield))break;
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates

    // Calculate the shifted radius we are at, squared. 
    // Compare squared radii for faster code
    Float_t shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
                                 + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

    // Roughly reached the radius we want
    if(shiftedRadiusSquared > RSquaredWanted){
      
      // Bigger loop has bad precision, we're nearly one centimeter too far, 
      // go back in small steps.
      while (shiftedRadiusSquared>RSquaredWanted){
	// Propagate a mm inwards
	x-=.1;
	if(!etp.PropagateTo(x,bfield)){
	  // Propagation failed but we're already with a
	  // cm precision at R=1.25m so we only break the 
	  // inner loop
	  break;
	}
	// Get the global position
	etp.GetXYZ(xyz);
	// Calculate shifted radius, squared
	shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
	                     + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);
      }
      // We reached R=1.25m with a precission of a cm to a mm,
      // set the spatial position
      fXSftR125[0]=xyz[0]-priVtx[0];
      fXSftR125[1]=xyz[1]-priVtx[1];
      fXSftR125[2]=xyz[2]-priVtx[2];
      // Done
      return;
    } // End of if roughly reached radius
  } // End of coarse propagation loop
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::Set(const AliAODTrack *track,const Float_t bfield,const Double_t priVtx[3]){
  // Overloaded function, using static casts because of narrowing
  const Float_t priVtxPos[3]={static_cast<Float_t>(priVtx[0]),
			      static_cast<Float_t>(priVtx[1]),
			      static_cast<Float_t>(priVtx[2])};
  Set(track,bfield,priVtxPos);
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::Set(const AliAODTrack *track,const Float_t bfield,const Float_t priVtx[3]){
  // Set the properties of this to the AliAODtrack
  //
  //    UShort_t fID;               //! Unique track id (->AliAODTrack.h), UShort_t goes to 65000
  //    Double_t fP[3];             //! Momentum of track
  //    Float_t  fXglobal[9][3];    //! Global positions at different global radii
  //    Float_t  fXshifted[9][3];   //! Shifted positions at different shifted radii


  // Set the ID, a good ID also indicates to use the track
  if(track->GetID() >=0){
    // global tracks, e.g. v0 daughters
    fID = track->GetID();
  }
  else {
    // e.g. tpc only tracks, i.e. primary protons
    fID = -track->GetID()-1;
  }
  // Set the momentum
  track->PxPyPz(fP);
  // Set spatial positions in the TPC
  SetSftPosR125(track,bfield,priVtx);
  // GetGlobalPositionAtGlobalRadii(track,bfield);
  // GetShiftedPositionAtShiftedRadii(track,bfield,priVtx);
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::FemtoBufferTrack(const FemtoBufferTrack& fbt):
  fID(fbt.fID)
 {
  // Copy constructor

  for (UChar_t i=0;i<3;i++){
    fP[i]=fbt.fP[i];
    fXSftR125[i]=fbt.fXSftR125[i];
    // for (UChar_t j=0;j<9;j++){
    //   fXglobal[j][i]=fbt.fXglobal[j][i];
    //   fXshifted[j][i]=fbt.fXshifted[j][i];
    // }
  }
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferTrack& AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::operator=(const FemtoBufferTrack& fbt){
  // Assignment operator, from wikipedia :)
  
  // Protect against self-assignment
  if(this != &fbt){
    fID = fbt.fID;
    for (UChar_t i=0;i<3;i++){
      fP[i]=fbt.fP[i];
      fXSftR125[i]=fbt.fXSftR125[i];
      // for (UChar_t j=0;j<9;j++){
      // 	fXglobal[j][i]=fbt.fXglobal[j][i];
      // 	fXshifted[j][i]=fbt.fXshifted[j][i];
      // }
    }
  }
  // By convention, always return *this (Could it be the convention is called c++?)
  return *this;
}
//_______________________________________________________________
Double_t AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::RapPro() const{
  // Calculate the rapidity assuming proton mass
  const Double_t e = TMath::Sqrt(fP[0]*fP[0] + fP[1]*fP[1]
				 +fP[2]*fP[2] + 0.9382720*0.9382720);
  if (e != TMath::Abs(fP[2])) { // energy was not equal to pz
    return 0.5*TMath::Log((e+fP[2])/(e-fP[2]));
  } else { // energy was equal to pz
    return -999.;
  }
}
//_______________________________________________________________
Double_t AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::Pt()const{
  // Calcs pt = sqrt (px^2 + py^2)
  return TMath::Sqrt(fP[0]*fP[0]+fP[1]*fP[1]);
}
//_______________________________________________________________
Double_t AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::ThetaS()const{
  // Returns the longitudinal angle of the particles propagated
  // position at R=1.25m. See
  // https://edms.cern.ch/file/406391/2/ALICE-INT-2003-038.pdf
  // for the ALICE coordinate system. Theta is zero at positive z,
  // pi/2 at z = 0 aka the xy plane and pi at negative z 

  // R^    ^  
  //  |   /
  //  |θ'/
  //  | / θ
  //  |/____>z
  // 
  // Let's compute θ' and θ = π/2 - θ'
  // where θ' can even be and should 
  // sometimes be negative
  // tan(θ') = z/R
  // θ' = arctan(z/R)
  // θ = π/2 - θ'
  //   = π/2 - arctan(z/R)
  // Note that in the doc above theta
  // is calculated as arccos(z/sqrt(x^2+y^2+z^2))

  // Array of positions is 85,105,125,..cm,
  // we take the z position at R=1.25m
  // return TMath::Pi()/2. - TMath::ATan(fXshifted[2][2]/125.);
  return TMath::Pi()/2. - TMath::ATan(fXSftR125[2]/125.);
}
//_______________________________________________________________
Double_t AliAnalysisTaskProtonLambda2d::FemtoBufferTrack::EtaS()const{
  // Returns the corresponding eta of a pri. part. 
  // with this particles pos at R=1.25m

  // http://en.wikipedia.org/wiki/Pseudorapidity
  // η = -ln[ tan(θ/2)]
  // printf("z: %+04.0f, thetaS %+03.2f etaS %+1.2f\n"
  // 	 ,fXshifted[2][2],ThetaS(),-TMath::Log( TMath::Tan(ThetaS()/2.) ));
  return -TMath::Log( TMath::Tan(ThetaS()/2.) );
}
//________________________________________________________________________
//
//                        FemtoBufferV0
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferV0::FemtoBufferV0():
  fCosPoint(-9999.),
  fPosDaughter(),
  fNegDaughter()
  // ,fP{-9999.,-9999.,-9999.} //C++11
{
  // Dummy constructor, set everything so it
  // indicates a V0 which should not be used
  fP[0]=-9999.; //C++03
  fP[1]=-9999.;
  fP[2]=-9999.;
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferV0::FemtoBufferV0(const AliAODv0 *v0, const AliAODTrack *posDaughter, const AliAODTrack *negDaughter, const Double_t bfield, Double_t priVtxPos[3]):
  fCosPoint(-9999.),
  fPosDaughter(),
  fNegDaughter()
{
  // Constructor, set the properties of this to these of the AliAODv0

  // Use Set function to keep code in one place. Only constant data member
  // would require the FemtoBuff() : fbla(), fblup() {} method
  Set(v0,posDaughter,negDaughter,bfield,priVtxPos);
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferV0::Set(const AliAODv0 *v0, const AliAODTrack *posDaughter, const AliAODTrack *negDaughter, const Double_t bfield, Double_t priVtxPos[3])
{
  // Set the properties of this to these of the AliAODv0
  fCosPoint=v0->CosPointingAngle(priVtxPos);
  v0->PxPyPz(fP);
  // printf("Set px %3.2f, py %3.2f, pz %3.2f\n"
  // 	 ,fP[0],fP[1],fP[2]
  // 	 );
  // The daughters
  fPosDaughter.Set(posDaughter,bfield,priVtxPos);
  fNegDaughter.Set(negDaughter,bfield,priVtxPos);
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferV0::FemtoBufferV0(const FemtoBufferV0 &fbv):
  fCosPoint(fbv.fCosPoint),
  fPosDaughter(fbv.fPosDaughter),
  fNegDaughter(fbv.fNegDaughter)
  //,fP{fbv.fP[0],fbv.fP[1],fbv.fP[2]} // C++11
{
  // Copy constructor
  fP[0] = fbv.fP[0]; // C++03
  fP[1] = fbv.fP[1];
  fP[2] = fbv.fP[2];
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferV0& AliAnalysisTaskProtonLambda2d::FemtoBufferV0::operator=(const FemtoBufferV0 &fbv){
  // Assignment operator

  // Protect against self-assignment
  if(this != &fbv){
    fCosPoint=fbv.fCosPoint;
    fP[0]=fbv.fP[0];
    fP[1]=fbv.fP[1];
    fP[2]=fbv.fP[2];
    fPosDaughter=fbv.fPosDaughter;
    fNegDaughter=fbv.fNegDaughter;
  }
  return *this;
}
//_______________________________________________________________
Double_t AliAnalysisTaskProtonLambda2d::FemtoBufferV0::RapLam() const{
  // Calculate the rapidity assuming lambda mass
  Double_t e = TMath::Sqrt(fP[0]*fP[0] + fP[1]*fP[1]
			   +fP[2]*fP[2] + 1.115683*1.115683);
  if (e != TMath::Abs(fP[2])) { // energy was not equal to pz
    return 0.5*TMath::Log((e+fP[2])/(e-fP[2]));
  } else { // energy was equal to pz
    return -999.;
  }
}
//_______________________________________________________________
Double_t AliAnalysisTaskProtonLambda2d::FemtoBufferV0::Pt()const{
  // Calcs pt = sqrt (px^2 + py^2)
  return TMath::Sqrt(fP[0]*fP[0]+fP[1]*fP[1]);
}
//________________________________________________________________________
//
//                        FemtoBufferEvent
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::FemtoBufferEvent():
fDoBgLamALam(kFALSE)
  ,fPriTrackLim(0),fV0Lim(0)
  ,fProTracks(0),fAProTracks(0)
  ,fLamTracks(0),fALamTracks(0)
  ,fBgLamTracks(0),fBgALamTracks(0)
  ,fNProTracks(0),fNAProTracks(0),fNLamTracks(0),fNALamTracks(0)
  ,fNBgLamTracks(0),fNBgALamTracks(0)
  ,fBfield(-9999.)
{
  // Standard constructor, all pointer to zero
  fPriVtxPos[0]=-9999.;
  fPriVtxPos[1]=-9999.;
  fPriVtxPos[2]=-9999.;

  printf("This constructor has zero size in the arrays!\n");
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::FemtoBufferEvent(const UShort_t priTrackBuff,const UShort_t V0Buff,const Bool_t DoBgLamALam,const Double_t bfield,const Double_t priVtxPos[3]):
fDoBgLamALam(DoBgLamALam)
  ,fPriTrackLim(priTrackBuff),fV0Lim(V0Buff)
  ,fProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fAProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fLamTracks (new FemtoBufferV0[fV0Lim])
  ,fALamTracks(new FemtoBufferV0[fV0Lim])
  ,fBgLamTracks(fDoBgLamALam?new FemtoBufferV0[fV0Lim]:0)
  ,fBgALamTracks(fDoBgLamALam?new FemtoBufferV0[fV0Lim]:0)
  ,fNProTracks(0),fNAProTracks(0),fNLamTracks(0),fNALamTracks(0)
  ,fNBgLamTracks(0),fNBgALamTracks(0)
  ,fBfield(-bfield)
  //  ,fPriVtxPos{priVtxPos[0],priVtxPos[1],priVtxPos[2]} // This is C++11
{
  // Constructor.
  fPriVtxPos[0] = priVtxPos[0]; // This is some old C++
  fPriVtxPos[1] = priVtxPos[1];
  fPriVtxPos[2] = priVtxPos[2];
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::FemtoBufferEvent(const UShort_t priTrackBuff,const UShort_t V0Buff
								  ,const Bool_t DoBgLamALam):
fDoBgLamALam(DoBgLamALam)
  ,fPriTrackLim(priTrackBuff),fV0Lim(V0Buff)
  ,fProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fAProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fLamTracks (new FemtoBufferV0[fV0Lim])
  ,fALamTracks(new FemtoBufferV0[fV0Lim])
  ,fBgLamTracks(fDoBgLamALam?new FemtoBufferV0[fV0Lim]:0)
  ,fBgALamTracks(fDoBgLamALam?new FemtoBufferV0[fV0Lim]:0)
  ,fNProTracks(0),fNAProTracks(0),fNLamTracks(0),fNALamTracks(0)
  ,fNBgLamTracks(0),fNBgALamTracks(0)
  ,fBfield(-9999.)
  //  ,fPriVtxPos{-9999.,-9999.,-9999.} // This is C++11
{  
  // Constructor. fBfield and fPriVtxPos not needed yet, can be set later.
  fPriVtxPos[0] = -9999.; // This is C++03
  fPriVtxPos[1] = -9999.;
  fPriVtxPos[2] = -9999.;

  //  printf("constructed eventwith NBgLam: %u\n",fNBgLamTracks);
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::FemtoBufferEvent(const FemtoBufferEvent &fbe):
fDoBgLamALam(fbe.GetDoBgLamALam())
  ,fPriTrackLim(fbe.GetPriTrackLim())
  ,fV0Lim(fbe.GetV0Lim())
  ,fProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fAProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fLamTracks (new FemtoBufferV0[fV0Lim])
  ,fALamTracks(new FemtoBufferV0[fV0Lim])
  ,fBgLamTracks(fDoBgLamALam?new FemtoBufferV0[fV0Lim]:0)
  ,fBgALamTracks(fDoBgLamALam?new FemtoBufferV0[fV0Lim]:0)
  ,fNProTracks(fbe.GetNPro()),fNAProTracks(fbe.GetNAPro())
  ,fNLamTracks(fbe.GetNLam()),fNALamTracks(fbe.GetNALam())
  ,fNBgLamTracks(fbe.GetNBgLam()),fNBgALamTracks(fbe.GetNBgALam())
  ,fBfield(fbe.GetBfield())
{
  // Copy constructor
  fbe.GetVtxPos(fPriVtxPos);
  // Avoid to much creation and deletion of objects
  UShort_t i;
  // Copy the primary tracks
  for (i=0;i<fPriTrackLim;i++){
    fProTracks[i]=fbe.fProTracks[i];
    fAProTracks[i]=fbe.fAProTracks[i];
  }
  // Copy the V0s
  for (i=0;i<fV0Lim;i++){
    fLamTracks[i]=fbe.fLamTracks[i];
    fALamTracks[i]=fbe.fALamTracks[i];
  }
  // Copy the bg V0s
  if(fDoBgLamALam){
    for (i=0;i<fV0Lim;i++){
      fBgLamTracks[i]=fbe.fBgLamTracks[i];
      fBgALamTracks[i]=fbe.fBgALamTracks[i];
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferEvent& AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::operator=(const FemtoBufferEvent &fbe){
  // Assignment operator

  // Protect against self-assignment
  if(this!=&fbe){
    // Well, we use arrays of a constant size to avoid
    // excessive memory allocation and won't give this up.
    // So we'll only copy as much as fits on the left side
    // from the right side.
    // DON'T COPY THE ARRAY SIZES fV0Lim AND fPriTrackLim !!!
    if(fPriTrackLim < fbe.GetPriTrackLim() 
       || fV0Lim < fbe.GetV0Lim()){
      // AliWarning(Form("Trying to assign too big event (buffer %d/%d) to"
      // 		    " this (buffer %d/%d). Only partially copying.",
      // 		    fbe.GetPriTrackLim(),fbe.GetV0Lim(),
      // 		    fPriTrackLim,fV0Lim));
      printf("Trying to assign too big event (buffer %d/%d) to"
    		    " this (buffer %d/%d). Only partially copying.\n",
	     fbe.GetPriTrackLim(),fbe.GetV0Lim(),
	     fPriTrackLim,fV0Lim);
    }
    if(!fDoBgLamALam&&fbe.GetDoBgLamALam()){
      printf("Trying to copy Bg lam/alam although this buffer has no storage for that"
	     "this->GetDoBgLamALam() %d\n",GetDoBgLamALam());
    }
   
    // Always start with the easy stuff :)
    fbe.GetVtxPos(fPriVtxPos);
    fBfield = fbe.GetBfield();
    // Number of tracks is minimum of array size of 'this'
    // and the number of tracks from the right side
    fNProTracks = TMath::Min(fPriTrackLim,fbe.GetNPro());
    fNAProTracks = TMath::Min(fPriTrackLim,fbe.GetNAPro());
    fNLamTracks = TMath::Min(fV0Lim,fbe.GetNLam());
    fNALamTracks = TMath::Min(fV0Lim,fbe.GetNALam());
    fNBgLamTracks = fDoBgLamALam?TMath::Min(fV0Lim,fbe.GetNBgLam()):0;
    fNBgALamTracks = fDoBgLamALam?TMath::Min(fV0Lim,fbe.GetNBgALam()):0;
    
    // Avoid creation and deletion of 'i' for every loop
    UShort_t i;
    // Copy primary tracks. No need to set a 'bad track'
    // flag for the entries above GetNPro() (...) as
    // above everything is bad by definition.
    // Protons
    for (i=0;i<GetNPro();i++)
      fProTracks[i]=fbe.fProTracks[i];
    // Anti-protons
    for (i=0;i<GetNAPro();i++)
      fAProTracks[i]=fbe.fAProTracks[i];
    // Copy the V0s 
    // Lambdas
    for (i=0;i<GetNLam();i++){
      fLamTracks[i]=fbe.fLamTracks[i];
    }
    // Anti-lambdas
    for (i=0;i<GetNALam();i++){
      fALamTracks[i]=fbe.fALamTracks[i];
    }
    // Background lambdas
    for (i=0;i<GetNBgLam();i++){
      fBgLamTracks[i]=fbe.fBgLamTracks[i];
    }
    // Background anti-lambdas
    for (i=0;i<GetNBgALam();i++){
      fBgALamTracks[i]=fbe.fBgALamTracks[i];
    }  
  }
  return *this;
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::~FemtoBufferEvent(){
  // Destructor

  // Delete the arrays of tracks,
  // note the [] with the delete
  if(fProTracks){
    delete[] fProTracks;
    fProTracks=0;
  }
  if(fAProTracks){
    delete[] fAProTracks;
    fAProTracks=0;
  }
  if(fLamTracks){
    delete[] fLamTracks;
    fLamTracks=0;
  }
  if(fALamTracks){
    delete[] fALamTracks;
    fALamTracks=0;
  }
  if(fBgLamTracks){
    delete[] fBgLamTracks;
    fBgLamTracks=0;
  }
  if(fBgALamTracks){
    delete[] fBgALamTracks;
    fBgALamTracks=0;
  }
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::Reset(const Double_t bfield, const Double_t priVtxPos[3]){
  // Reset the old event, i.e., make clear 'here is no info'
  // by setting the 'number of stored ...' to zero
  fNProTracks=0;
  fNAProTracks=0;
  fNLamTracks=0;
  fNALamTracks=0;
  fNBgLamTracks=0;
  fNBgALamTracks=0;
  
  // And set the new event properties 
  fBfield = bfield;
  fPriVtxPos[0]=priVtxPos[0];
  fPriVtxPos[1]=priVtxPos[1];
  fPriVtxPos[2]=priVtxPos[2];
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::AddPro(const AliAODTrack *track){
  // Add a proton to this event

  // Check whether there is still space in the array
  if(fNProTracks > fPriTrackLim-1){
    // AliWarning(Form("Cannot add proton, array size (%d) too small"
    // 		    ,fPriTrackLim));
    printf("Cannot add proton, array size (%d) too small\n"
    		    ,fPriTrackLim);
    return;
  }
  // Add the V0 at the end of the array
  fProTracks[fNProTracks].Set(track,fBfield,fPriVtxPos);
  fNProTracks++;
  //  printf("Added proton %d/%d\n",fNProTracks,fPriTrackLim);

}  
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::AddAPro(const AliAODTrack *track){
  // Add a anti-proton to this event

  // Check whether there is still space in the array
  if(fNAProTracks > fPriTrackLim-1){
    // AliWarning(Form("Cannot add anti-proton, array size (%d) too small"
    // 		    ,fPriTrackLim));
    printf("Cannot add anti-proton, array size (%d) too small\n"
		    ,fPriTrackLim);
    return;
  }
  // Add the V0 at the end of the array
  fAProTracks[fNAProTracks].Set(track,fBfield,fPriVtxPos);
  fNAProTracks++;
}  
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::AddLam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter){
  // Adds a lambda with it's daughters to the event

  // Check whether there is still space in the array
  if(fNLamTracks > fV0Lim-1){
    // AliWarning(Form("Cannot add lambda, array size (%d) too small"
    // 		    ,fV0Lim));
    printf("Cannot add lambda, array size (%d) too small"
		    ,fV0Lim);
    return;
  }
 
  // Add the V0 at the end of the array
  fLamTracks[fNLamTracks].Set(v0,posDaughter,negDaughter,
			      fBfield,fPriVtxPos);
  fNLamTracks++;
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::AddALam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter){
  // Adds a lambda with it's daughters to the event

  // Check whether there is still space in the array
  if(fNALamTracks > fV0Lim-1){
    // AliWarning(Form("Cannot add anti-lambda, array size (%d) too small"
    // 		    ,fV0Lim));
    printf("Cannot add anti-lambda, array size (%d) too small\n"
    		    ,fV0Lim);
    return;
  }
 
  // Add the V0 at the end of the array
  fALamTracks[fNALamTracks].Set(v0,posDaughter,negDaughter,
				fBfield,fPriVtxPos);
  fNALamTracks++;
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::AddBgLam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter){
  // Adds a bg-lambda with it's daughters to the event

  // Check whether we want to store bg lam/alam
  if(!GetDoBgLamALam()){
    printf("This buffer is not intended to store bg lam/alam\n");
    return;
  }
  // As with the lambdas, we reject bg-lambdas
  // with pt < 0.5 GeV/c
  if(v0->Pt() < 0.5)
    return;

  // Check whether there is still space in the array
  if(fNBgLamTracks > fV0Lim-1){
    // AliWarning(Form("Cannot add background lambda,"
    // 		    " array size (%d) too small"
    // 		    ,fV0Lim));
    // printf("Cannot add background lambda,"
    // 	   "already stored %d" 
    // 	   " array size (%d) too small\n"
    // 	   ,fNBgLamTracks
    // 	   ,fV0Lim);
    return;
  }
 
  // Add the V0 at the end of the array
  fBgLamTracks[fNALamTracks].Set(v0,posDaughter,negDaughter,
				fBfield,fPriVtxPos);
  fNBgLamTracks++;
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBufferEvent::AddBgALam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter){
  // Adds a lambda with it's daughters to the event

  // Check whether we want to store bg lam/alam
  if(!GetDoBgLamALam()){
    printf("This buffer is not intended to store bg lam/alam\n");
    return;
  }

  // Check whether there is still space in the array
  if(fNBgALamTracks > fV0Lim-1){
    // AliWarning(Form("Cannot add background anti-lambda,"
    // 		    " array size (%d) too small"
    // 		    ,fV0Lim));
    //    printf("Cannot add background anti-lambda,"
    //		    " array size (%d) too small\n"
    //		    ,fV0Lim);
    return;
  }
 
  // Add the V0 at the end of the array
  fBgALamTracks[fNALamTracks].Set(v0,posDaughter,negDaughter,
				fBfield,fPriVtxPos);
  fNBgALamTracks++;
}
//________________________________________________________________________
//
//                        FemtoBuffer
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBuffer::FemtoBuffer() :
  fkZvertexBins(0),
  fkCentBins(0),
  fkMixBuffSize(0),
  fDoBgLamALam(kFALSE),
  fkPriTrackLim(0),
  fkV0Lim(0),
  fZvertexAxis(0),
  fCentAxis(0),
  fCurEvt(0),
  fEC(0)
{
  // Dummy constructor, create arrays with zero size
  // Note that some data member are constant, you
  // won't be able to create the FemtoBuffer first with this
  // constructor and then set the appropiate size.
  
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBuffer::FemtoBuffer(const UChar_t ZvertexBins,const UChar_t CentBins,const UChar_t MixBuff,const UShort_t PriTrackLim,const UShort_t V0Lim, const Float_t AbsZvertexCut,const Float_t CentCutLo,const Float_t CentCutHi,const Bool_t DoBgLamALam) :
  fkZvertexBins(ZvertexBins),
  fkCentBins(CentBins),
  fkMixBuffSize(MixBuff),
  fDoBgLamALam(DoBgLamALam),
  fkPriTrackLim(PriTrackLim),
  fkV0Lim(V0Lim),
  fZvertexAxis(new TAxis(fkZvertexBins,-AbsZvertexCut,AbsZvertexCut)),
  fCentAxis(new TAxis (fkCentBins,CentCutLo,CentCutHi)),
  fCurEvt(new FemtoBufferEvent *[fkMixBuffSize]),
  fEC(new FemtoBufferEvent ***[fkZvertexBins])
{
  // Constructor, creates at once all events with all tracks
  //  printf ("Creating with pritracklim %d and v0lim %d\n",fkPriTrackLim,fkV0Lim);

  // Create the array step by step
  // Bins in z of the primary vertex position. Do this as
  // the detector looks different from a different z coordinate
  for (UChar_t iZBin=0;iZBin<fkZvertexBins;iZBin++){
    fEC[iZBin] = new FemtoBufferEvent **[fkCentBins];
    // Bins in centrality
    for (UChar_t iCentBin=0;iCentBin<fkCentBins;iCentBin++){
      fEC[iZBin][iCentBin] = new FemtoBufferEvent *[fkMixBuffSize];
      // The number of events to keep for one mixing class
      for(UChar_t iMixBuff=0;iMixBuff<fkMixBuffSize;iMixBuff++){
	// Create an event to hold the info for mixing
	fEC[iZBin][iCentBin][iMixBuff] = new FemtoBufferEvent(fkPriTrackLim,fkV0Lim,fDoBgLamALam);
      }
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBuffer::FemtoBuffer(const AliAnalysisTaskProtonLambda2d::FemtoBuffer &fb) :
  fkZvertexBins(fb.fkZvertexBins),
  fkCentBins(fb.fkCentBins),
  fkMixBuffSize(fb.fkMixBuffSize),
  fDoBgLamALam(fb.fDoBgLamALam),
  fkPriTrackLim(fb.fkPriTrackLim),
  fkV0Lim(fb.fkV0Lim),
  fZvertexAxis(new TAxis(*(fb.fZvertexAxis))),
  fCentAxis(new TAxis (*(fb.fCentAxis))),
  fCurEvt(new FemtoBufferEvent *[fkMixBuffSize]),
  fEC(new FemtoBufferEvent ***[fkZvertexBins])
{
  // Copy constructor. Linux complains not having this and 
  // compiling this task with aliroot

  printf("FemtoBuffer ctor not tested yet, be cautious\n");
  
  // Create the array step by step
  // Bins in z of the primary vertex position. Do this as
  // the detector looks different from a different z coordinate
  for (UChar_t iZBin=0;iZBin<fkZvertexBins;iZBin++){
    fEC[iZBin] = new FemtoBufferEvent **[fkCentBins];
    // Bins in centrality
    for (UChar_t iCentBin=0;iCentBin<fkCentBins;iCentBin++){
      fEC[iZBin][iCentBin] = new FemtoBufferEvent *[fkMixBuffSize];
      // The number of events to keep for one mixing class
      for(UChar_t iMixBuff=0;iMixBuff<fkMixBuffSize;iMixBuff++){
	// Create an event to hold the info for mixing
	fEC[iZBin][iCentBin][iMixBuff] = new FemtoBufferEvent(*(fb.fEC[iZBin][iCentBin][iMixBuff]));
      }
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBuffer& AliAnalysisTaskProtonLambda2d::FemtoBuffer::operator=(const AliAnalysisTaskProtonLambda2d::FemtoBuffer& fb){
  //Assignment operator
  if(this!=&fb){
    printf("FemtoBuffer assignment operator not implemented\n");
  }
  return *this;
  
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda2d::FemtoBuffer::~FemtoBuffer(){
  // Destructor
  // The axes to fin the correct bins
  if(fZvertexAxis){
    delete fZvertexAxis;
    fZvertexAxis=0;
  }
  if(fCentAxis){
    delete fCentAxis;
    fCentAxis=0;
  }
  // fCurEvt is an array of pointer
  if(fCurEvt){
    delete[] fCurEvt;
    fCurEvt=0;
  }
  // Delete all the events and the pointer to them
  for (UChar_t iZBin=0;iZBin<fkZvertexBins;iZBin++){
    for (UChar_t iCentBin=0;iCentBin<fkCentBins;iCentBin++){
      for(UChar_t iMixBuff=0;iMixBuff<fkMixBuffSize;iMixBuff++){
	if(fEC[iZBin][iCentBin][iMixBuff]){
	  delete fEC[iZBin][iCentBin][iMixBuff];
	  fEC[iZBin][iCentBin][iMixBuff]=0;
	}
      }
      if(fEC[iZBin][iCentBin]){
	delete[] fEC[iZBin][iCentBin];
	fEC[iZBin][iCentBin]=0;
      }
    }
    if(fEC[iZBin]){
      delete[] fEC[iZBin];
      fEC[iZBin]=0;
    }
  }
  if(fEC){
    delete[] fEC;
    fEC=0;
  }
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBuffer::ShiftAndAdd(AliAODEvent *evt,
							     const TString centEst){
  // Shift the events in the appropiate centrality / zvertex bin and set the 
  // current event pointer correctly
  Double_t priVtxPos[3];
  evt->GetPrimaryVertex()->GetXYZ(priVtxPos);
  //  printf("Mag field: %f\n",evt->GetMagneticField());
  ShiftAndAdd(evt->GetMagneticField(),
	      priVtxPos,
	      evt->GetCentrality()->GetCentralityPercentileUnchecked(centEst.Data()));
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda2d::FemtoBuffer::ShiftAndAdd(const Double_t bfield,const Double_t priVtxPos[3],const Float_t centrality){
  // Shift the events in the appropiate centrality / zvertex bin and set the 
  // current event pointer correctly

  // Find the correct centrality/zvertex bin 
  const UChar_t ZvertexBin = fZvertexAxis->FindFixBin(priVtxPos[2]) - 1; // -1 for array starting at 0
  const UChar_t CentBin = fCentAxis->FindFixBin(centrality) - 1;// -1 for array starting at 0

  // The new current event is the old last event
  fCurEvt[0] = fEC[ZvertexBin][CentBin][fkMixBuffSize-1];

  // Shift the pointer, starting from the back
  UChar_t iMix;
  for(iMix=fkMixBuffSize-1;iMix>0;iMix--){
    fEC[ZvertexBin][CentBin][iMix] = fEC[ZvertexBin][CentBin][iMix-1];
  }
  // And reset the zero'th one
  fEC[ZvertexBin][CentBin][0] = fCurEvt[0];
  fEC[ZvertexBin][CentBin][0]->Reset(bfield,priVtxPos);
  // Also set the pointer to the other events..
  for (iMix=1;iMix<fkMixBuffSize;iMix++){
    fCurEvt[iMix] = fEC[ZvertexBin][CentBin][iMix];
  }
}
