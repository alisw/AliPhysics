/****************************************************************************
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
 *                                                                        *
 * Original task: AliAnalysisTaskFlatenicity.cxx                          *
 * Author: Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)                *
 * Anatask to compute flatenicity (arXiv:2204.13733)                      *
 * Modified by: Gyula Bencedi (bencedi.gyula@wigner.hu)                   *
 **************************************************************************/

class TTree;

class AliPPVsMultUtils;
class AliESDtrackCuts;

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDUtils.h"
#include "AliESDVZERO.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliPPVsMultUtils.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TParticle.h"
#include "TProfile.h"
#include "TVector3.h"
#include <AliAnalysisFilter.h>
#include <AliESDVertex.h>
#include <AliHeader.h>
#include <AliMultiplicity.h>
#include <Riostream.h>
#include <TBits.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>
using std::cout;
using std::endl;

#include "AliAnalysisTaskSpectraFlatenicity.h"

const Int_t nPtbinsFlatSpecFlatSpec = 60;
Double_t PtbinsFlatSpec[nPtbinsFlatSpecFlatSpec+1] = {
    0.0,  0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 
    0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 
    1.1 , 1.2, 1.3,  1.4, 1.5 , 1.6, 1.7 , 1.8, 1.9 , 2.0, 
    2.2 , 2.4, 2.6,  2.8, 3.0 , 3.2, 3.4 , 3.6, 3.8 , 4.0, 
    4.5 , 5.0, 5.5,  6.0, 6.5 , 7.0, 8.0 , 9.0, 10.0,11.0,
    12.0,13.0,14.0, 15.0, 16.0, 18.0,20.0, 22.0,24.0,26.0,
    30.0};     
     
const Int_t nCent = 9;
Double_t centClassFlatSpec[nCent + 1] = {0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

using namespace std; // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskSpectraFlatenicity) // classimp: necessary for root

    AliAnalysisTaskSpectraFlatenicity::AliAnalysisTaskSpectraFlatenicity()
    : AliAnalysisTaskSE(), 
    fESD(0), 
    fEventCuts(0x0), 
    fMCStack(0), 
    fMC(0),
    fUseMC(kFALSE), 
    fV0Mindex(-1), 
    fmultV0A(-1), 
    fmultV0C(-1), 
    fmultTPC(-1), 
    fmultV0Amc(-1), 
    fmultV0Cmc(-1), 
    fmultTPCmc(-1), 
    fRemoveTrivialScaling(kFALSE),
    fSysVarTrkCuts(0),
    fnGen(-1), 
    fnRec(-1), 
    fnRecWoDCA(-1), 
    fPIDResponse(0x0), 
    fTrackFilter(0x0),
    fTrackFilterwoDCA(0x0),
    fOutputList(0), 
    fEtaCut(0.8), 
    fPtMin(0.5), 
    fv0mpercentile(0), 
    fdcaxy(-999), 
    fdcaz(-999),
    fFlat(-1), 
    fFlatMC(-1), 
    fMultSelection(0x0), 
    hFlatenicity(0), 
    hFlatenicityMC(0), 
    hFlatResponse(0), 
    hFlatVsPt(0), 
    hFlatVsPtMC(0), 
    hFlatVsNch(0), 
    hFlatVsNchMC(0), 
    hNchV0M(0), 
    hNchV0MMC(0), 
    hNchMidRap(0), 
    hNchMidRapMC(0), 
    hNchV0a(0), 
    hNchV0c(0), 
    hNchV0aMC(0), 
    hNchV0cMC(0), 
    hPtPrimIn(0), 
    hPtOutRec(0), 
    hPtInPrim(0),
    hPtInPrimLambda(0),
    hPtInPrimPion(0),
    hPtInPrimKaon(0),
    hPtInPrimProton(0),
    hPtInPrimSigmap(0),
    hPtInPrimSigmam(0),
    hPtInPrimOmega(0),
    hPtInPrimXi(0),
    hPtInPrimRest(0),
    hPtOut(0),
    hPtOutPrim(0),
    hPtOutSec(0),
    hPtOutPrimLambda(0),
    hPtOutPrimPion(0),
    hPtOutPrimKaon(0),
    hPtOutPrimProton(0),
    hPtOutPrimSigmap(0),
    hPtOutPrimSigmam(0),
    hPtOutPrimOmega(0),
    hPtOutPrimXi(0),
    hPtOutPrimRest(0),
    hPtVsDCAData(0), 
    hPtVsDCAPrim(0), 
    hPtVsDCADec(0), 
    hPtVsDCAMat(0), 
    hPtVsDCAAll(0), 
    hFlatVsV0M(0), 
    hEta(0), 
    hEtamc(0), 
    hCounter(0),
    hV0MBad(0)
{
    for (Int_t i_c = 0; i_c < nCent; ++i_c) {
        hFlatVsPtV0M[i_c] = 0;
        hFlatVsNchTPCV0M[i_c] = 0;
    }
    for (Int_t i_c = 0; i_c < nCent; ++i_c) {
        hFlatVsPtV0MMC[i_c] = 0;
        hFlatVsNchTPCV0MMC[i_c] = 0;
    }
}

//_____________________________________________________________________________
AliAnalysisTaskSpectraFlatenicity::AliAnalysisTaskSpectraFlatenicity(const char *name)
    : AliAnalysisTaskSE(name), 
    fESD(0), 
    fEventCuts(0x0), 
    fMCStack(0), 
    fMC(0),
    fUseMC(kFALSE), 
    fV0Mindex(-1), 
    fmultV0A(-1), 
    fmultV0C(-1), 
    fmultTPC(-1), 
    fmultV0Amc(-1), 
    fmultV0Cmc(-1), 
    fmultTPCmc(-1), 
    fRemoveTrivialScaling(kFALSE),
    fSysVarTrkCuts(0),
    fnGen(-1), 
    fnRec(-1), 
    fnRecWoDCA(-1), 
    fPIDResponse(0x0), 
    fTrackFilter(0x0),
    fTrackFilterwoDCA(0x0),
    fOutputList(0), 
    fEtaCut(0.8), 
    fPtMin(0.5), 
    fv0mpercentile(0), 
    fdcaxy(-999), 
    fdcaz(-999),
    fFlat(-1), 
    fFlatMC(-1), 
    fMultSelection(0x0), 
    hFlatenicity(0), 
    hFlatenicityMC(0), 
    hFlatResponse(0), 
    hFlatVsPt(0), 
    hFlatVsPtMC(0), 
    hFlatVsNch(0), 
    hFlatVsNchMC(0), 
    hNchV0M(0), 
    hNchV0MMC(0), 
    hNchMidRap(0), 
    hNchMidRapMC(0), 
    hNchV0a(0), 
    hNchV0c(0), 
    hNchV0aMC(0), 
    hNchV0cMC(0), 
    hPtPrimIn(0), 
    hPtOutRec(0), 
    hPtInPrim(0),
    hPtInPrimLambda(0),
    hPtInPrimPion(0),
    hPtInPrimKaon(0),
    hPtInPrimProton(0),
    hPtInPrimSigmap(0),
    hPtInPrimSigmam(0),
    hPtInPrimOmega(0),
    hPtInPrimXi(0),
    hPtInPrimRest(0),
    hPtOut(0),
    hPtOutPrim(0),
    hPtOutSec(0),
    hPtOutPrimLambda(0),
    hPtOutPrimPion(0),
    hPtOutPrimKaon(0),
    hPtOutPrimProton(0),
    hPtOutPrimSigmap(0),
    hPtOutPrimSigmam(0),
    hPtOutPrimOmega(0),
    hPtOutPrimXi(0),
    hPtOutPrimRest(0),
    hPtVsDCAData(0), 
    hPtVsDCAPrim(0), 
    hPtVsDCADec(0), 
    hPtVsDCAMat(0), 
    hPtVsDCAAll(0), 
    hFlatVsV0M(0), 
    hEta(0), 
    hEtamc(0), 
    hCounter(0),
    hV0MBad(0)
{
    for (Int_t i_c = 0; i_c < nCent; ++i_c) {
        hFlatVsPtV0M[i_c] = 0;
        hFlatVsNchTPCV0M[i_c] = 0;
    }
    for (Int_t i_c = 0; i_c < nCent; ++i_c) {
        hFlatVsPtV0MMC[i_c] = 0;
        hFlatVsNchTPCV0MMC[i_c] = 0;
    }

  DefineInput(0, TChain::Class()); // define the input of the analysis: in this
                                   // case you take a 'chain' of events
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this
                                   // case it's a list of histograms
}

//_____________________________________________________________________________
AliAnalysisTaskSpectraFlatenicity::~AliAnalysisTaskSpectraFlatenicity() {
  // destructor
  if (fOutputList) {
    delete fOutputList; // at the end of your task, it is deleted from memory by
                        // calling this function
    fOutputList = 0x0;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskSpectraFlatenicity::UserCreateOutputObjects() {

  // create track filters
  fTrackFilter = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts *fCuts = new AliESDtrackCuts();
  
  fCuts->SetAcceptKinkDaughters(kFALSE);
  fCuts->SetRequireTPCRefit(kTRUE);
  fCuts->SetRequireITSRefit(kTRUE);
  fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fCuts->SetDCAToVertex2D(kFALSE);
  fCuts->SetRequireSigmaToVertex(kFALSE);
  fCuts->SetEtaRange(-0.8, 0.8);
  fCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
  
  if(fSysVarTrkCuts==1){ //! Lower: SetMinNCrossedRowsTPC(60)
      fCuts->SetMinNCrossedRowsTPC(60);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      fCuts->SetMaxChi2PerClusterTPC(4);
      fCuts->SetMaxDCAToVertexZ(2);
      fCuts->SetMaxChi2PerClusterITS(36);
  }
  else if(fSysVarTrkCuts==2){ //! Higher: SetMinNCrossedRowsTPC(100)
      fCuts->SetMinNCrossedRowsTPC(100);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      fCuts->SetMaxChi2PerClusterTPC(4);
      fCuts->SetMaxDCAToVertexZ(2);
      fCuts->SetMaxChi2PerClusterITS(36);
  }
  else if(fSysVarTrkCuts==3){ //! Lower: SetMinRatioCrossedRowsOverFindableClustersTPC(0.7)
      fCuts->SetMinNCrossedRowsTPC(70);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
      fCuts->SetMaxChi2PerClusterTPC(4);
      fCuts->SetMaxDCAToVertexZ(2);
      fCuts->SetMaxChi2PerClusterITS(36);
  }
  else if(fSysVarTrkCuts==4){ //! Higher: SetMinRatioCrossedRowsOverFindableClustersTPC(0.9)
      fCuts->SetMinNCrossedRowsTPC(70);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
      fCuts->SetMaxChi2PerClusterTPC(4);
      fCuts->SetMaxDCAToVertexZ(2);
      fCuts->SetMaxChi2PerClusterITS(36);
  }
  else if(fSysVarTrkCuts==5){ //! Lower: SetMaxChi2PerClusterTPC(3)
      fCuts->SetMinNCrossedRowsTPC(70);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      fCuts->SetMaxChi2PerClusterTPC(3);
      fCuts->SetMaxDCAToVertexZ(2);
      fCuts->SetMaxChi2PerClusterITS(36);
  }
  else if(fSysVarTrkCuts==6){ //! Higher: SetMaxChi2PerClusterTPC(5)
      fCuts->SetMinNCrossedRowsTPC(70);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      fCuts->SetMaxChi2PerClusterTPC(5);
      fCuts->SetMaxDCAToVertexZ(2);
      fCuts->SetMaxChi2PerClusterITS(36);
  }
  else if(fSysVarTrkCuts==7){ //! Lower: SetMaxChi2PerClusterITS(25)
      fCuts->SetMinNCrossedRowsTPC(70);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      fCuts->SetMaxChi2PerClusterTPC(4);
      fCuts->SetMaxDCAToVertexZ(2);
      fCuts->SetMaxChi2PerClusterITS(25);
  }
  else if(fSysVarTrkCuts==8){ //! Higher: SetMaxChi2PerClusterITS(49)
      fCuts->SetMinNCrossedRowsTPC(70);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      fCuts->SetMaxChi2PerClusterTPC(4);
      fCuts->SetMaxDCAToVertexZ(2);
      fCuts->SetMaxChi2PerClusterITS(49);
  }
  else if(fSysVarTrkCuts==9){ //! Lower: SetMaxDCAToVertexZ(1)
      fCuts->SetMinNCrossedRowsTPC(70);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      fCuts->SetMaxChi2PerClusterTPC(4);
      fCuts->SetMaxDCAToVertexZ(1);
      fCuts->SetMaxChi2PerClusterITS(36);
  }
  else if(fSysVarTrkCuts==10){ //! Lower: SetMaxDCAToVertexZ(5)
      fCuts->SetMinNCrossedRowsTPC(70);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      fCuts->SetMaxChi2PerClusterTPC(4);
      fCuts->SetMaxDCAToVertexZ(5);
      fCuts->SetMaxChi2PerClusterITS(36);
  }
  else{ //! Default values
      fCuts->SetMinNCrossedRowsTPC(70);
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);	
      fCuts->SetMaxChi2PerClusterTPC(4);
      fCuts->SetMaxDCAToVertexZ(2);
      fCuts->SetMaxChi2PerClusterITS(36);
  }

  fTrackFilter->AddCuts(fCuts);

  // wo DCA cut  
  fTrackFilterwoDCA = new AliAnalysisFilter("trackFilterwoDCA");
  AliESDtrackCuts *fCutswoDCA = new AliESDtrackCuts();
  SetCutsFilterWoDCA(fCutswoDCA); 
  
  // create output objects

  OpenFile(1);
  fOutputList = new TList(); // this is a list which will contain all of your histograms
  fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all

  hFlatenicity = new TH1D("hFlatenicity", "counter", 1020, -0.01, 1.01);
  fOutputList->Add(hFlatenicity);

  hEta = new TH1D("hEta", "Eta rec; #eta; counts", 200, -1.0, 1.0); hEta->Sumw2();
  fOutputList->Add(hEta);      
  
  hFlatVsPt = new TH2D("hFlatVsPt", "Measured; Flatenicity; #it{p}_{T} (GeV/#it{c})", 1020, -0.01, 1.01, nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
  fOutputList->Add(hFlatVsPt);

  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    hFlatVsPtV0M[i_c] = new TH2D(Form("hFlatVsPtV0M_c%d", i_c), Form("Measured %1.0f-%1.0f%%V0M; Flatenicity; #it{p}_{T} (GeV/#it{c})",centClassFlatSpec[i_c], centClassFlatSpec[i_c + 1]),1020, -0.01, 1.01, nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hFlatVsPtV0M[i_c]);
    hFlatVsNchTPCV0M[i_c] = new TH2D(Form("hFlatVsNchTPCV0M_c%d", i_c), Form("Measured %1.0f-%1.0f%%V0M; Flatenicity; #it{N}_{ch}",centClassFlatSpec[i_c], centClassFlatSpec[i_c + 1]),1020, -0.01, 1.01, 100, -0.5, 99.5);
    fOutputList->Add(hFlatVsNchTPCV0M[i_c]);    
  }

  hNchMidRap = new TH1D("hNchMidRap", ";Nch; counts", 400, -0.5, 399.5);
  hNchMidRap->Sumw2();
  fOutputList->Add(hNchMidRap);
  
  if (fUseMC) {
      
    hEtamc = new TH1D("hEtamc", "Eta mc.; #eta; counts", 200, -1.0, 1.0); hEtamc->Sumw2();
    fOutputList->Add(hEtamc);      
    
    for (Int_t i_c = 0; i_c < nCent; ++i_c) {
        hFlatVsPtV0MMC[i_c] = new TH2D( Form("hFlatVsPtV0MMC_c%d", i_c), Form("Measured %1.0f-%1.0f%%V0M; Flatenicity; #it{p}_{T} (GeV/#it{c})",centClassFlatSpec[i_c], centClassFlatSpec[i_c + 1]), 1020, -0.01, 1.01, nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
        fOutputList->Add(hFlatVsPtV0MMC[i_c]);
        hFlatVsNchTPCV0MMC[i_c] = new TH2D( Form("hFlatVsNchTPCV0MMC_c%d", i_c), Form("Measured %1.0f-%1.0f%%V0M; Flatenicity; #it{N}_{ch}",centClassFlatSpec[i_c], centClassFlatSpec[i_c + 1]), 1020, -0.01, 1.01, 100, -0.5, 99.5);
        fOutputList->Add(hFlatVsNchTPCV0MMC[i_c]);        
    }

    hPtPrimIn = new TH1D("hPtPrimIn", "Prim In; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtPrimIn);

    hPtOutRec = new TH1D("hPtOutRec", "all Out; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutRec);
    
    hPtInPrim = new TH1D("hPtInPrim", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtInPrim);
    
    hPtInPrimLambda = new TH1D("hPtInPrimLambda", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtInPrimLambda);

    hPtInPrimPion = new TH1D("hPtInPrimPion", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtInPrimPion);

    hPtInPrimKaon = new TH1D("hPtInPrimKaon", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtInPrimKaon);

    hPtInPrimProton = new TH1D("hPtInPrimProton", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtInPrimProton);

    hPtInPrimSigmap = new TH1D("hPtInPrimSigmap", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtInPrimSigmap);

    hPtInPrimSigmam = new TH1D("hPtInPrimSigmam", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtInPrimSigmam);

    hPtInPrimOmega = new TH1D("hPtInPrimOmega", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtInPrimOmega);

    hPtInPrimXi = new TH1D("hPtInPrimXi", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtInPrimXi);

    hPtInPrimRest = new TH1D("hPtInPrimRest", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtInPrimRest);

    hPtOut = new TH1D("hPtOut", "pT all rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOut);
    
    hPtOutPrim = new TH1D("hPtOutPrim", "pT prim rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutPrim);

    hPtOutSec = new TH1D("hPtOutSec", "pT sec rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutSec);
    
    hPtOutPrimLambda = new TH1D("hPtOutPrimLambda", "pT prim true; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutPrimLambda);

    hPtOutPrimPion = new TH1D("hPtOutPrimPion", "pT prim rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutPrimPion);

    hPtOutPrimKaon = new TH1D("hPtOutPrimKaon", "pT prim rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutPrimKaon);

    hPtOutPrimProton = new TH1D("hPtOutPrimProton", "pT prim rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutPrimProton);

    hPtOutPrimSigmap = new TH1D("hPtOutPrimSigmap", "pT prim rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutPrimSigmap);

    hPtOutPrimSigmam = new TH1D("hPtOutPrimSigmam", "pT prim rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutPrimSigmam);

    hPtOutPrimOmega = new TH1D("hPtOutPrimOmega", "pT prim rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutPrimOmega);

    hPtOutPrimXi = new TH1D("hPtOutPrimXi", "pT prim rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutPrimXi);

    hPtOutPrimRest = new TH1D("hPtOutPrimRest", "pT prim rec; #it{p}_{T} (GeV/#it{c}); counts", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hPtOutPrimRest);
    
    hFlatenicityMC = new TH1D("hFlatenicityMC", "counter", 1020, -0.01, 1.01);
    fOutputList->Add(hFlatenicityMC);

    hFlatResponse = new TH2D("hFlatResponse", "; true flat; measured flat", 1020, -0.01, 1.01, 1020, -0.01, 1.01);
    fOutputList->Add(hFlatResponse);

    hFlatVsPtMC = new TH2D("hFlatVsPtMC", "MC true; Flatenicity; #it{p}_{T} (GeV/#it{c})", 1020, -0.01, 1.01, nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec);
    fOutputList->Add(hFlatVsPtMC);

    hFlatVsNchMC = new TH2D("hFlatVsNchMC", "; true flat; true Nch", 1020, -0.01, 1.01, 100, -0.5, 99.5);
    fOutputList->Add(hFlatVsNchMC);
    
    /// Added V0M multiplicity distribtion 
    hNchV0MMC = new TH1D("hNchV0MMC", ";true Nch; counts", 400, -0.5, 399.5);
    hNchV0MMC->Sumw2();
    fOutputList->Add(hNchV0MMC);

    hNchMidRapMC = new TH1D("hNchMidRapMC", ";true Nch; counts", 400, -0.5, 399.5);
    hNchMidRapMC->Sumw2();
    fOutputList->Add(hNchMidRapMC);

    hNchV0aMC = new TH1D("hNchV0aMC", ";rec Nch; counts", 400, -0.5, 399.5);
    hNchV0aMC->Sumw2();
    fOutputList->Add(hNchV0aMC);

    hNchV0cMC = new TH1D("hNchV0cMC", ";rec Nch; counts", 400, -0.5, 399.5);
    hNchV0cMC->Sumw2();
    fOutputList->Add(hNchV0cMC);  
    
  }

  hPtVsDCAData = new TH2D("hPtVsDCAData","; #it{p}_{T} (GeV/#it{c}); DCA_{xy} data",nPtbinsFlatSpecFlatSpec,PtbinsFlatSpec, 560, -3.5, 3.5);
  fOutputList->Add(hPtVsDCAData);    
  
  hPtVsDCAPrim = new TH2D("hPtVsDCAPrim", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} primaries", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec, 560, -3.5, 3.5);
  fOutputList->Add(hPtVsDCAPrim);

  hPtVsDCADec = new TH2D("hPtVsDCADec", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} decays", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec, 560, -3.5, 3.5);
  fOutputList->Add(hPtVsDCADec);

  hPtVsDCAMat = new TH2D("hPtVsDCAMat", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} material", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec, 560, -3.5, 3.5);
  fOutputList->Add(hPtVsDCAMat);

  hPtVsDCAAll = new TH2D("hPtVsDCAAll", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} all", nPtbinsFlatSpecFlatSpec, PtbinsFlatSpec, 560, -3.5, 3.5);
  fOutputList->Add(hPtVsDCAAll);
  
  hFlatVsNch = new TH2D("hFlatVsNch", "; rec flat; rec Nch", 1020, -0.01, 1.01, 100, -0.5, 99.5);
  fOutputList->Add(hFlatVsNch);
  
  hFlatVsV0M = new TH2D("hFlatVsV0M", "", nCent, centClassFlatSpec, 1020, -0.01, 1.01);
  fOutputList->Add(hFlatVsV0M);

  hNchV0M = new TH1D("hNchV0M", ";rec Nch; counts", 400, -0.5, 399.5);
  hNchV0M->Sumw2();
  fOutputList->Add(hNchV0M);

  hNchV0a = new TH1D("hNchV0a", ";rec Nch; counts", 400, -0.5, 399.5);
  hNchV0a->Sumw2();
  fOutputList->Add(hNchV0a);

  hNchV0c = new TH1D("hNchV0c", ";rec Nch; counts", 400, -0.5, 399.5);
  hNchV0c->Sumw2();
  fOutputList->Add(hNchV0c);  
  
  hCounter = new TH1D("hCounter", "counter", 10, -0.5, 9.5);
  fOutputList->Add(hCounter);
  
  hV0MBad = new TH1D("hV0MBad", "hV0MBad", 1000, -0.5, 999.5);
  fOutputList->Add(hV0MBad);  

  fEventCuts.AddQAplotsToList(fOutputList);
  PostData(1, fOutputList); // postdata will notify the analysis manager of
                            // changes / updates to the
}

//_____________________________________________________________________________
void AliAnalysisTaskSpectraFlatenicity::UserExec(Option_t *) {

  AliVEvent *event = InputEvent();
  if (!event) { Error("UserExec", "Could not retrieve event"); return; }

  fESD = dynamic_cast<AliESDEvent *>(event);
  if (!fESD) { Printf("%s:%d ESDEvent not found in Input Manager", (char *)__FILE__,__LINE__); this->Dump(); return; }

  if (fUseMC) {
    //      E S D
    fMC = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMC) { Printf("%s:%d MCEvent not found in Input Manager", (char *)__FILE__,__LINE__); this->Dump(); return; }
    fMCStack = fMC->Stack();
  }

  AliHeader *headerMC;
  Bool_t isGoodVtxPosMC = kFALSE;

  vector<Float_t> ptMc;
  vector<Int_t> idMc;

  if (fUseMC) {
    headerMC = fMC->Header();
    AliGenEventHeader *genHeader = headerMC->GenEventHeader();
    TArrayF vtxMC(3); // primary vertex  MC
    vtxMC[0] = 9999;
    vtxMC[1] = 9999;
    vtxMC[2] = 9999; // initialize with dummy
    if (genHeader) {
      genHeader->PrimaryVertex(vtxMC);
    }
    if (TMath::Abs(vtxMC[2]) <= 15.)
      isGoodVtxPosMC = kTRUE;

    fnGen = FillArrayMC(ptMc, idMc);
  }
  
  // w/ DCA cut
  vector<Float_t> ptWDCA;
  vector<Float_t> dcaxyWDCA;
  vector<Int_t> isprimWDCA;
  vector<Int_t> idWDCA;

  fnRec = FillArray(ptWDCA, dcaxyWDCA, isprimWDCA, idWDCA, kTRUE);
  

  // w/o DCA cut
  vector<Float_t> ptWoDCA;
  vector<Float_t> dcaxyWoDCA;
  vector<Int_t> isprimWoDCA;
  vector<Int_t> idWoDCA;

  fnRecWoDCA = FillArray(ptWoDCA, dcaxyWoDCA, isprimWoDCA, idWoDCA, kFALSE);

  // Trigger selection
  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;
  if (!isINT7selected)
    return;

  // Good events
  if (!fEventCuts.AcceptEvent(event)) {
    PostData(1, fOutputList);
    return;
  }

  // Good vertex
  Bool_t hasRecVertex = kFALSE;
  hasRecVertex = HasRecVertex();
  if (!hasRecVertex)
    return;

  // Multiplicity Estimation
  fv0mpercentile = -999;
  fMultSelection = (AliMultSelection *)fESD->FindListObject("MultSelection");
  if (!fMultSelection)
    cout << "------- No AliMultSelection Object Found --------" << fMultSelection << endl;
  fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  hCounter->Fill(1);
  float v0mult = fMultSelection->GetEstimator("V0M")->GetValue();

  for (Int_t i_c = 0; i_c < nCent; ++i_c) {
    if (fv0mpercentile >= centClassFlatSpec[i_c] && fv0mpercentile < centClassFlatSpec[i_c + 1]) {
      fV0Mindex = i_c;
    } else {
      continue;
    }
  }
/*
  // INEL>0 selection
	if( AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0) < 1. )
        return;
*/

  CheckMultiplicities(ptWoDCA, dcaxyWoDCA, fnRecWoDCA);
  
  fFlat = GetFlatenicity();
  
  // DATA
  if (fFlat >= 0) {
    hFlatenicity->Fill(fFlat);
    if (fV0Mindex >= 0) {
        if ((fV0Mindex == nCent - 1) && (v0mult > 400.)) { // to reject 70-100% multiplicity class events that have large V0M amplitude
            hV0MBad->Fill(v0mult);
        } else {
            hFlatVsV0M->Fill(fv0mpercentile, fFlat);
            MakeDataanalysis();
        }
    }
  }

  // MC
  if (fUseMC) {
    fFlatMC = GetFlatenicityMC();
    if (fFlatMC >= 0) {
        hFlatenicityMC->Fill(fFlatMC);
        hFlatResponse->Fill(fFlatMC, fFlat);
        MakeMCanalysis();
    }
    
    if (isGoodVtxPosMC){
        GetCorrections(fnGen, fnRec, ptMc, ptWDCA, idMc, idWDCA, isprimWDCA);
// //         GetCorrections(fnGen, fnRecWoDCA, ptMc, ptWoDCA, idMc, idWoDCA, isprimWoDCA);
        CheckMultiplicitiesMC(ptWoDCA, dcaxyWoDCA, isprimWoDCA, fnRecWoDCA);
    }
  } // MC

  PostData(1, fOutputList); // stream the result of this event to the output
                            // manager which will write it to a file
}

//______________________________________________________________________________
void AliAnalysisTaskSpectraFlatenicity::Terminate(Option_t *) {
}

//______________________________________________________________________________
void AliAnalysisTaskSpectraFlatenicity::MakeDataanalysis() {

  // rec
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t iT = 0; iT < nTracks; ++iT) {

    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(fESD->GetTrack(iT)); // get a track (type AliesdTrack)
    if (!esdtrack)
      continue;
    if (!fTrackFilter->IsSelected(esdtrack))
      continue;
    if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
      continue;
    if (esdtrack->Pt() < fPtMin)
      continue;
    
    hFlatVsPt->Fill(fFlat, esdtrack->Pt());
    hFlatVsPtV0M[fV0Mindex]->Fill(fFlat, esdtrack->Pt());
    hEta->Fill(esdtrack->Eta());
    hFlatVsNchTPCV0M[fV0Mindex]->Fill(fFlat, fmultTPC);
  }
}

//______________________________________________________________________________
void AliAnalysisTaskSpectraFlatenicity::MakeMCanalysis() {

  for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (TMath::Abs(particle->Eta()) > fEtaCut)
      continue;
    if (particle->Pt() < fPtMin)
      continue;
    if (TMath::Abs(particle->Charge()) < 0.1)
      continue;
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
      continue;    
    
    hFlatVsPtMC->Fill(fFlatMC, particle->Pt());
    hFlatVsPtV0MMC[fV0Mindex]->Fill(fFlatMC, particle->Pt());
    hPtPrimIn->Fill(particle->Pt());
    hEtamc->Fill(particle->Eta());
    hFlatVsNchTPCV0MMC[fV0Mindex]->Fill(fFlatMC, fmultTPCmc);
  }
  // rec
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t iT = 0; iT < nTracks; ++iT) {

    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(fESD->GetTrack(iT)); // get a track (type AliesdTrack)
    if (!esdtrack)
      continue;
    if (!fTrackFilter->IsSelected(esdtrack))
      continue;
    if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
      continue;
    if (esdtrack->Pt() < fPtMin)
      continue;
    
    hPtOutRec->Fill(esdtrack->Pt());
  }
}

//______________________________________________________________________________
void AliAnalysisTaskSpectraFlatenicity::CheckMultiplicitiesMC(const vector<Float_t> &ptRecWoDCA, const vector<Float_t> &dcaxyRecWoDCA, const vector<Int_t> &isprimRecWoDCA, Int_t multRecWoDCA) {

  fmultV0Amc = 0;
  fmultV0Cmc = 0;
  fmultTPCmc = 0;

  for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (particle->Pt() <= 0.0)
      continue;
    if (TMath::Abs(particle->Charge()) < 0.1)
      continue;
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
      continue;    

    Double_t eta_a = particle->Eta();
    if (eta_a >= 2.8 && eta_a < 4.5) { // v0a acceptance (excluding first ring)
      fmultV0Amc++;
    }
    if (eta_a >= -3.7 && eta_a < -1.7) { // v0c
      fmultV0Cmc++;
    }
    if (TMath::Abs(eta_a) < 0.8) { // adc
      fmultTPCmc++;
    }
  }
  
  hNchMidRapMC->Fill(fmultTPCmc);
  hNchV0aMC->Fill(fmultV0Amc);
  hNchV0cMC->Fill(fmultV0Cmc);

  for (Int_t i = 0; i < multRecWoDCA; ++i) {

    hPtVsDCAAll->Fill(ptRecWoDCA[i], dcaxyRecWoDCA[i]);

    if (fUseMC) {
      if (isprimRecWoDCA[i] == 0) {
        hPtVsDCAPrim->Fill(ptRecWoDCA[i], dcaxyRecWoDCA[i]);
      } else if (isprimRecWoDCA[i] == 1) {
        hPtVsDCADec->Fill(ptRecWoDCA[i], dcaxyRecWoDCA[i]);
      } else if (isprimRecWoDCA[i] == 2) {
        hPtVsDCAMat->Fill(ptRecWoDCA[i], dcaxyRecWoDCA[i]);
      }
    }
  }
}

//______________________________________________________________________________
void AliAnalysisTaskSpectraFlatenicity::CheckMultiplicities(const vector<Float_t> &ptRecWoDCA, const vector<Float_t> &dcaxyRecWoDCA, Int_t multRecWoDCA) {

  fmultTPC = 0;
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t iT = 0; iT < nTracks; ++iT) {
    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(fESD->GetTrack(iT)); // get a track (type AliesdTrack)
    if (!esdtrack)
      continue;
    if (!fTrackFilter->IsSelected(esdtrack))
      continue;
    if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
      continue;
    if (esdtrack->Pt() < fPtMin)
      continue;
    fmultTPC++;
  }

  AliVVZERO *lVV0 = 0x0;
  AliVEvent *lVevent = 0x0;
  lVevent = dynamic_cast<AliVEvent *>(InputEvent());
  if (!lVevent) {
    AliWarning("ERROR: ESD / AOD event not available \n");
    return;
  }
  // Get VZERO Information for multiplicity later
  lVV0 = lVevent->GetVZEROData();
  if (!lVV0) {
    AliError("AliVVZERO not available");
    return;
  }

  const Int_t nChannels = 64;
  fmultV0C = 0;
  fmultV0A = 0;
  for (Int_t iCh = 0; iCh < nChannels; iCh++) {
    Float_t mult = lVV0->GetMultiplicity(iCh);
    if (iCh < 32) { // V0C
      fmultV0C += mult;
    } else if (iCh >= 32 &&
               iCh < 40) { // exclude first ring to avoid overlap with ADA
      continue;
    } else { // V0A
      fmultV0A += mult;
    }
  }
  
  hNchMidRap->Fill(fmultTPC);
  hNchV0a->Fill(fmultV0A);
  hNchV0c->Fill(fmultV0C);
  
  // Auxiliar distribution to calculate the contamination from secondary
  // particles, it runs over tracks wo DCA cut
  for (Int_t i = 0; i < multRecWoDCA; ++i) {
    hPtVsDCAData->Fill(ptRecWoDCA[i], dcaxyRecWoDCA[i]);
  }
}

//______________________________________________________________________________
Double_t AliAnalysisTaskSpectraFlatenicity::GetFlatenicity() {

  AliVVZERO *lVV0 = 0x0;
  AliVEvent *lVevent = 0x0;
  lVevent = dynamic_cast<AliVEvent *>(InputEvent());
  if (!lVevent) {
    AliWarning("ERROR: ESD / AOD event not available \n");
    return -1;
  }
  // Get VZERO Information for multiplicity later
  lVV0 = lVevent->GetVZEROData();
  if (!lVV0) {
    AliError("AliVVZERO not available");
    return -1;
  }
  // Flatenicity calculation
  const Int_t nRings = 4;
  const Int_t nSectors = 8;
  Float_t minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
  Float_t maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
  Float_t maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
  Float_t minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};
  
  // Grid
  const Int_t nCells = nRings * 2 * nSectors;
  Float_t RhoLattice[nCells];
  Float_t multLattice[nCells];
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    RhoLattice[iCh] = 0.0;
    multLattice[iCh] = 0.0;
  }

  Int_t nringA = 0;
  Int_t nringC = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    Float_t detaV0 = -1;
    Float_t mult = lVV0->GetMultiplicity(iCh);
    if (iCh < 32) { // V0C
      if (iCh < 8) {
        nringC = 0;
      } else if (iCh >= 8 && iCh < 16) {
        nringC = 1;
      } else if (iCh >= 16 && iCh < 24) {
        nringC = 2;
      } else {
        nringC = 3;
      }
      detaV0 = maxEtaV0C[nringC] - minEtaV0C[nringC];
    } else { // V0A
      if (iCh < 40) {
        nringA = 0;
      } else if (iCh >= 40 && iCh < 48) {
        nringA = 1;
      } else if (iCh >= 48 && iCh < 56) {
        nringA = 2;
      } else {
        nringA = 3;
      }
      detaV0 = maxEtaV0A[nringA] - minEtaV0A[nringA];
    }
    RhoLattice[iCh] = mult / detaV0; // needed to consider the different eta coverage
    multLattice[iCh] = mult;
  }

  Float_t mRho = 0;
//   Float_t multRho = 0;
  Float_t flatenicity = -1;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    mRho    += RhoLattice[iCh];
//     multRho += multLattice[iCh];
  }
  Float_t multV0Mdeta = mRho;
//   Float_t multV0M = multRho;

  // average activity per cell
  mRho /= (1.0 * nCells);
  
  // get sigma
  Double_t sRho_tmp = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  Float_t sRho = TMath::Sqrt(sRho_tmp);
  if (mRho > 0) {
    if (fRemoveTrivialScaling) {
      flatenicity = TMath::Sqrt(multV0Mdeta) * sRho / mRho; // scaling by absolute tot mult
    } else {
      flatenicity = sRho / mRho;
    }
  } else {
    flatenicity = -1;
  }
  
  hFlatVsNch->Fill(flatenicity, multV0Mdeta);
  hNchV0M->Fill(multV0Mdeta);

  return flatenicity;
}

//______________________________________________________________________________
Double_t AliAnalysisTaskSpectraFlatenicity::GetFlatenicityMC() {

  // Flatenicity calculation
  const Int_t nRings = 8;
  Float_t maxEta[nRings] = {-3.2, -2.7, -2.2, -1.7, 5.1, 4.5, 3.9, 3.4};
  Float_t minEta[nRings] = {-3.7, -3.2, -2.7, -2.2, 4.5, 3.9, 3.4, 2.8};

  const Int_t nSectors = 8;
  Float_t PhiBins[nSectors + 1];
  Float_t deltaPhi = (2.0 * TMath::Pi()) / (1.0 * nSectors);
  for (int i_phi = 0; i_phi < nSectors + 1; ++i_phi) {
    PhiBins[i_phi] = 0;
    if (i_phi < nSectors) {
      PhiBins[i_phi] = i_phi * deltaPhi;
    } else {
      PhiBins[i_phi] = 2.0 * TMath::Pi();
    }
  }

  // Grid
  const Int_t nCells = nRings * nSectors;
  Float_t RhoLattice[nCells];
  Float_t multLattice[nCells];
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    RhoLattice[iCh] = 0.0;
    multLattice[iCh] = 0.0;
  }

  Int_t nMult = 0;
  for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (particle->Pt() <= 0.0)
      continue;
    if (TMath::Abs(particle->Charge()) < 0.1)
      continue;
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
      continue;    

    Double_t phi = particle->Phi();
    Double_t eta = particle->Eta();

    Int_t i_segment = 0;
    for (int i_eta = 0; i_eta < nRings; ++i_eta) {

      for (int i_phi = 0; i_phi < nSectors; ++i_phi) {

        if (eta >= minEta[i_eta] && eta < maxEta[i_eta] &&
            phi >= PhiBins[i_phi] && phi < PhiBins[i_phi + 1]) {
          nMult++;
          RhoLattice[i_segment] += 1.0;
          multLattice[i_segment] += 1.0;
        }
        i_segment++;
      }
    }
  }

  Int_t i_seg = 0;
  for (int i_eta = 0; i_eta < nRings; ++i_eta) {
    for (int i_phi = 0; i_phi < nSectors; ++i_phi) {
      Float_t deltaEta = TMath::Abs(maxEta[i_eta] - minEta[i_eta]);
      RhoLattice[i_seg] /= deltaEta;
      i_seg++;
    }
  }

  Float_t mRho = 0;
  Float_t multRho = 0;
  Float_t flatenicity = -1;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    mRho    += RhoLattice[iCh];
    multRho += multLattice[iCh];
  }
  Float_t multiplicityV0M = mRho;
  Float_t multV0M = multRho;

  // average activity per cell
  mRho /= (1.0 * nCells);

  // get sigma
  Float_t sRho_tmp = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  Float_t sRho = TMath::Sqrt(sRho_tmp);
  if (mRho > 0) {
    if (fRemoveTrivialScaling) {
      flatenicity = TMath::Sqrt(1.0 * multV0M) * sRho / mRho;
    } else {
      flatenicity = sRho / mRho;
    }
  } else {
    sRho = -1;
  }
  
  hFlatVsNchMC->Fill(flatenicity, nMult);
  hNchV0MMC->Fill(nMult);
  
  return flatenicity;
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskSpectraFlatenicity::HasRecVertex() {

  float fMaxDeltaSpdTrackAbsolute = 0.5f;
  float fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  float fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  float fMaxResolutionSPDvertex = 0.25f;
  float fMaxDispersionSPDvertex = 1.e14f;

  Bool_t fRequireTrackVertex = true;
  unsigned long fFlag;
  fFlag = BIT(AliEventCuts::kNoCuts);

  const AliVVertex *vtTrc = fESD->GetPrimaryVertex();
  bool isTrackV = true;
  if (vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ())
    isTrackV = false;
  const AliVVertex *vtSPD = fESD->GetPrimaryVertexSPD();

  if (vtSPD->GetNContributors() > 0)
    fFlag |= BIT(AliEventCuts::kVertexSPD);

  if (vtTrc->GetNContributors() > 1 && isTrackV)
    fFlag |= BIT(AliEventCuts::kVertexTracks);

  if (((fFlag & BIT(AliEventCuts::kVertexTracks)) || !fRequireTrackVertex) &&
      (fFlag & BIT(AliEventCuts::kVertexSPD)))
    fFlag |= BIT(AliEventCuts::kVertex);

  const AliVVertex *&vtx =
      bool(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
  AliVVertex *fPrimaryVertex = const_cast<AliVVertex *>(vtx);
  if (!fPrimaryVertex)
    return kFALSE;

  /// Vertex quality cuts
  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = bool(fFlag & AliEventCuts::kVertexSPD) &&
                      bool(fFlag & AliEventCuts::kVertexTracks)
                  ? vtTrc->GetZ() - vtSPD->GetZ()
                  : 0.; /// If one of the two vertices is not available this cut
                        /// is always passed.
  double errTot = TMath::Sqrt(covTrc[5] + covSPD[5]);
  double errTrc =
      bool(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
  double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
  /// vertex dispersion for run1, only for ESD, AOD code to be added here
  const AliESDVertex *vtSPDESD = dynamic_cast<const AliESDVertex *>(vtSPD);
  double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
  if ((TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute &&
       nsigTot <= fMaxDeltaSpdTrackNsigmaSPD &&
       nsigTrc <=
           fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
      (!vtSPD->IsFromVertexerZ() ||
       TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
      (!vtSPD->IsFromVertexerZ() ||
       vtSPDdispersion <= fMaxDispersionSPDvertex) /// vertex dispersion cut for
                                                   /// run1, only for ESD
      ) // quality cut on vertexer SPD z
    fFlag |= BIT(AliEventCuts::kVertexQuality);

  Bool_t hasVtx = (TESTBIT(fFlag, AliEventCuts::kVertex)) &&
                  (TESTBIT(fFlag, AliEventCuts::kVertexQuality));

  return hasVtx;
}

//______________________________________________________________________________
Int_t AliAnalysisTaskSpectraFlatenicity::FillArrayMC(vector<Float_t> &ptArray, vector<Int_t> &idArray)
{
//      id 0: pion, 1: kaon, 2: proton, 3: sigma plus, 4: sigma minus, 5: Omega, 6: Xi, 7: other charged

  ptArray.clear();
  idArray.clear();
  Int_t nNchGen = 0;
  
  for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (TMath::Abs(particle->Eta()) > fEtaCut)
      continue;
    if (particle->Pt() < fPtMin)
      continue;
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
      continue;    

    Int_t idPart = -1;
    Int_t partPDG = TMath::Abs(particle->PdgCode());
    if (partPDG == 3122)
      idPart = 0; // lambda
    if (particle->Charge() != 0) {
      if (partPDG == 211)
        idPart = 1; // pions
      else if (partPDG == 321)
        idPart = 2; // kaons
      else if (partPDG == 2212)
        idPart = 3; // protons
      else if (partPDG == 3222)
        idPart = 4; // sigma plus
      else if (partPDG == 3112)
        idPart = 5; // sigma minus
      else if (partPDG == 3334)
        idPart = 6; // Omega
      else if (partPDG == 3312)
        idPart = 7; // Xi
      else
        idPart = 8; // rest of the charged particles
    }
    ptArray.push_back(particle->Pt());
    idArray.push_back(idPart);
    nNchGen++;
  }
  return nNchGen;
}

//______________________________________________________________________________
Int_t AliAnalysisTaskSpectraFlatenicity::FillArray(vector<Float_t> &ptArray, vector<Float_t> &dcaxyArray, vector<Int_t> &isprimArray, vector<Int_t> &idArray, const Bool_t wDcaCut)
{
//      id 0: pion, 1: kaon, 2: proton, 3: sigma plus, 4: sigma minus, 5: Omega, 6: Xi, 7: other charged
  ptArray.clear();
  dcaxyArray.clear();
  isprimArray.clear();
  idArray.clear();

  Int_t nTracks = fESD->GetNumberOfTracks();
  Int_t nNchRec = 0;

  if (wDcaCut) 
  { // with DCA cut
    for (Int_t iT = 0; iT < nTracks; ++iT) 
    {
      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(fESD->GetTrack(iT));
      if (!esdtrack)
        continue;

      fdcaxy = -999;
      fdcaz = -999;

      if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
        continue;
      if (esdtrack->Pt() < fPtMin)
        continue;

      AliESDtrack *newTrack = 0x0;
      Int_t isPrim = -1;
      Int_t idTrack = -1;
      Int_t mcLabel = -1;
      
      if (fUseMC) 
      {
        // get label: 0: prim, 1: weak decays, 2: material
        mcLabel = TMath::Abs(esdtrack->GetLabel());
        TParticle *mcParticle = fMC->GetTrack(mcLabel)->Particle();
        if (!mcParticle) {
          printf("----ERROR: mcParticle not available------------------\n");
          continue;
        }
        
        Int_t partPDG_rec = TMath::Abs(mcParticle->GetPdgCode());
        if (partPDG_rec == 3122)
          idTrack = 0; // lambdas
        else if (partPDG_rec == 211)
          idTrack = 1; // pions
        else if (partPDG_rec == 321)
          idTrack = 2; // kaons
        else if (partPDG_rec == 2212)
          idTrack = 3; // protons
        else if (partPDG_rec == 3222)
          idTrack = 4; // sigma plus
        else if (partPDG_rec == 3112)
          idTrack = 5; // sigma minus
        else if (partPDG_rec == 3334)
          idTrack = 6; // Omega
        else if (partPDG_rec == 3312)
          idTrack = 7; // Xi
        else
          idTrack = 8; // rest of the charged particles
      } // MC 
      
      if (!fTrackFilter->IsSelected(esdtrack))
          continue;
      
      newTrack = new AliESDtrack(*esdtrack);
      newTrack->GetImpactParameters(fdcaxy, fdcaz);

      ptArray.push_back(newTrack->Pt());
      dcaxyArray.push_back(fdcaxy);
        
      if (fUseMC) {
          // get label: 0: prim, 1: weak decays, 2: material
          if (fMC->IsPhysicalPrimary(mcLabel))
              isPrim = 0;
          else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
              isPrim = 1;
          else if (fMC->IsSecondaryFromMaterial(mcLabel))
              isPrim = 2;
          else
              continue;
      }

      isprimArray.push_back(isPrim);
      idArray.push_back(idTrack);
      nNchRec++;
      delete newTrack;
    }
  } else { // w/o DCA cut
    for (Int_t iT = 0; iT < nTracks; ++iT)
    {
      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(fESD->GetTrack(iT)); 
      if (!esdtrack)
        continue;

      fdcaxy = -999;
      fdcaz = -999;

      if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
        continue;
      if (esdtrack->Pt() < fPtMin)
        continue;

      AliESDtrack *newTrack = 0x0;
      Int_t isPrim = -1;
      Int_t idTrack = -1;
      Int_t mcLabel = -1;
      TParticle *mcParticle = 0;

      if (fUseMC) { // get label: 0: prim, 1: weak decays, 2: material
        mcLabel = TMath::Abs(esdtrack->GetLabel());
        mcParticle = fMC->GetTrack(mcLabel)->Particle();
        if (!mcParticle) {
          printf("----ERROR: mcParticle not available------------------\n");
          continue;
        }

        Int_t partPDG_rec = TMath::Abs(mcParticle->GetPdgCode());
        if (partPDG_rec == 3122)
          idTrack = 0; // lambdas
        else if (partPDG_rec == 211)
          idTrack = 1; // pions
        else if (partPDG_rec == 321)
          idTrack = 2; // kaons
        else if (partPDG_rec == 2212)
          idTrack = 3; // protons
        else if (partPDG_rec == 3222)
          idTrack = 4; // sigma plus
        else if (partPDG_rec == 3112)
          idTrack = 5; // sigma minus
        else if (partPDG_rec == 3334)
          idTrack = 6; // Omega
        else if (partPDG_rec == 3312)
          idTrack = 7; // Xi
        else
          idTrack = 8; // rest of the charged particles
      }
      
      if (!fTrackFilterwoDCA->IsSelected(esdtrack))
          continue;
      
      newTrack = new AliESDtrack(*esdtrack);
      newTrack->GetImpactParameters(fdcaxy, fdcaz);

      ptArray.push_back(newTrack->Pt());
      dcaxyArray.push_back(fdcaxy);

      if (fUseMC) { // get label: 0: prim, 1: weak decays, 2: material
          if (fMC->IsPhysicalPrimary(mcLabel))
              isPrim = 0;
          else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
              isPrim = 1;
          else if (fMC->IsSecondaryFromMaterial(mcLabel))
              isPrim = 2;
          else
              continue;
      }

      isprimArray.push_back(isPrim);
      idArray.push_back(idTrack);
      nNchRec++;
      delete newTrack;
    }
  }
  return nNchRec;
}

//______________________________________________________________________________
void AliAnalysisTaskSpectraFlatenicity::SetCutsFilterWoDCA(AliESDtrackCuts *cFilt) 
{
  cFilt->SetMinNCrossedRowsTPC(70);
  cFilt->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  cFilt->SetMaxChi2PerClusterTPC(4);
  cFilt->SetMaxDCAToVertexZ(2);
  cFilt->SetAcceptKinkDaughters(kFALSE);
  cFilt->SetRequireTPCRefit(kTRUE);
  cFilt->SetRequireITSRefit(kTRUE);
  cFilt->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  cFilt->SetDCAToVertex2D(kFALSE);
  cFilt->SetRequireSigmaToVertex(kFALSE);
  cFilt->SetEtaRange(-0.8, 0.8);
  fTrackFilterwoDCA->AddCuts(cFilt);
}

//_____________________________________________________________________________
void AliAnalysisTaskSpectraFlatenicity::GetCorrections(Int_t multGen, Int_t multRec, const vector<Float_t> &ptGen, const vector<Float_t> &ptRec, const vector<Int_t> &idGen, const vector<Int_t> &idRec, const vector<Int_t> &isprimRec)
{
  // Histos for efficiencyxacceptance
    
  for (Int_t i = 0; i < multGen; ++i) {
    if (idGen[i] >= 0 && idGen[i] <= 8) {
      if (idGen[i] > 0) {
        hPtInPrim->Fill(ptGen[i]); // inital pT distribution (MC gen only charged particles)
      }
      if (idGen[i] == 0)
          hPtInPrimLambda->Fill(ptGen[i]); // lambdas
      else if (idGen[i] == 1)
          hPtInPrimPion->Fill(ptGen[i]); // pions
      else if (idGen[i] == 2)
          hPtInPrimKaon->Fill(ptGen[i]); // kaons
      else if (idGen[i] == 3)
          hPtInPrimProton->Fill(ptGen[i]); // protons
      else if (idGen[i] == 4)
          hPtInPrimSigmap->Fill(ptGen[i]); // sigma plus
      else if (idGen[i] == 5)
          hPtInPrimSigmam->Fill(ptGen[i]); // sigma minus
      else if (idGen[i] == 6)
          hPtInPrimOmega->Fill(ptGen[i]); // Omega
      else if (idGen[i] == 7)
          hPtInPrimXi->Fill(ptGen[i]); // Xi
      else
          hPtInPrimRest->Fill(ptGen[i]); // rest of the charged particles
    }
  }

  for (Int_t i = 0; i < multRec; ++i) { // loop over all these tracks
    hPtOut->Fill(ptRec[i]);
    if (idRec[i] >= 0 && idRec[i] <= 8) {
      if (isprimRec[i] == 0) {
        hPtOutPrim->Fill(ptRec[i]);

        if (idRec[i] == 0)
            hPtOutPrimLambda->Fill(ptRec[i]); // lambdas
        else if (idRec[i] == 1)
            hPtOutPrimPion->Fill(ptRec[i]); // pions
        else if (idRec[i] == 2)
            hPtOutPrimKaon->Fill(ptRec[i]); // kaons
        else if (idRec[i] == 3)
            hPtOutPrimProton->Fill(ptRec[i]); // protons
        else if (idRec[i] == 4)
            hPtOutPrimSigmap->Fill(ptRec[i]); // sigma plus
        else if (idRec[i] == 5)
            hPtOutPrimSigmam->Fill(ptRec[i]); // sigma minus
        else if (idRec[i] == 6)
            hPtOutPrimOmega->Fill(ptRec[i]); // Omega
        else if (idRec[i] == 7)
            hPtOutPrimXi->Fill(ptRec[i]); // Xi
        else
            hPtOutPrimRest->Fill(ptRec[i]); // rest of the charged particles
      }
    }
    if (isprimRec[i] == 1 || isprimRec[i] == 2) {
      hPtOutSec->Fill(ptRec[i]);
    }
  }
}

