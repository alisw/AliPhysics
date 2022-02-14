/**************************************************************************
 * Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////
// Author: Jianhui Zhu (1,2), Jeremy Wilkinson (2), Annalena Kalteyer (2)
// (1) Central China Normal University
// (2) GSI Helmholtz Centre for Heavy Ion Research
// E-mail: zjh@mail.ccnu.edu.cn, jeremy.wilkinson@cern.ch, annalena.sophie.kalteyer@cern.ch
/////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <TDatabasePDG.h>
#include <vector>
#include <TVector3.h>
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TLine.h"
#include "TList.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliCentrality.h"
#include "AliNormalizationCounter.h"
#include "AliAnalysisTaskSELc2pKs0fromKFP.h"
#include "AliPIDResponse.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisVertexingHF.h"

#include "AliAODMCParticle.h"

// includes added to play with KFParticle
#ifndef HomogeneousField
#define HomogeneousField
#endif

class AliAnalysisTaskSELc2pKs0fromKFP;    // your analysis class

ClassImp(AliAnalysisTaskSELc2pKs0fromKFP) // classimp: necessary for root

AliAnalysisTaskSELc2pKs0fromKFP::AliAnalysisTaskSELc2pKs0fromKFP() :
  AliAnalysisTaskSE(),
  fPID(0),
  fPIDCombined(0),
  fAnaCuts(0),
  fpVtx(0),
  fpVtxOff(0),
  fMCEvent(0),
  fBzkG(0),
  fCentrality(0),
  fAodTrackInd(0),
  fOutputList(0),
  fOutputWeight(0),
  fTree_Event(0),
  fVar_Event(0),
  fTree_Lc(0),
  fVar_Lc(0),
  fTree_Lc_QA(0),
  fVar_Lc_QA(0),
  fTree_LcMCGen(0),
  fVar_LcMCGen(0),
  fListCuts(0),
  fIsMC(kFALSE),
  fUseWeights(kFALSE),
  fKeepOnlyMCSignal(kTRUE),
  fKeepAllVariables(kFALSE),
  fIsAnaLc2Lpi(kFALSE),
  fCounter(0),
  fHistEvents(0),
  fHTrigger(0),
  fWriteLcMCGenTree(kFALSE),
  fWriteLcTree(kFALSE),
  fWriteLcQATree(kFALSE),
  fWeight(0),
  fHistMCGen_LcPt_weight(0),
  f2DHistMCRec_LcPt_weight(0),
  fFuncWeightPythia(0),
  fFuncWeightFONLL5overLHC13d3(0),
  fFuncWeightFONLL5overLHC13d3Lc(0),
  fUseMult(kFALSE),
  fRefMult(0),
  fAnalysisType(kpPb2016),
  fUseOnTheFlyV0(kFALSE)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
    for (Int_t i=0; i<4; i++) fMultEstimatorAvg[i] = 0;

}
//_____________________________________________________________________________
AliAnalysisTaskSELc2pKs0fromKFP::AliAnalysisTaskSELc2pKs0fromKFP(const char* name, AliRDHFCutsKFP* cuts) :
  AliAnalysisTaskSE(name),
  fPID(0),
  fPIDCombined(0),
  fAnaCuts(cuts),
  fpVtx(0),
  fpVtxOff(0),
  fMCEvent(0),
  fBzkG(0),
  fCentrality(0),
  fAodTrackInd(0),
  fOutputList(0),
  fOutputWeight(0),
  fTree_Event(0),
  fVar_Event(0),
  fTree_Lc(0),
  fVar_Lc(0),
  fTree_Lc_QA(0),
  fVar_Lc_QA(0),
  fTree_LcMCGen(0),
  fVar_LcMCGen(0),
  fListCuts(0),
  fIsMC(kFALSE),
  fUseWeights(kFALSE),
  fKeepOnlyMCSignal(kTRUE),
  fKeepAllVariables(kFALSE),
  fIsAnaLc2Lpi(kFALSE),
  fCounter(0),
  fHistEvents(0),
  fHTrigger(0),
  fWriteLcMCGenTree(kFALSE),
  fWriteLcTree(kFALSE),
  fWriteLcQATree(kFALSE),
  fWeight(0),
  fHistMCGen_LcPt_weight(0),
  f2DHistMCRec_LcPt_weight(0),
  fFuncWeightPythia(0),
  fFuncWeightFONLL5overLHC13d3(0),
  fFuncWeightFONLL5overLHC13d3Lc(0),
  fUseMult(kFALSE),
  fRefMult(0),
  fAnalysisType(kpPb2016),
  fUseOnTheFlyV0(kFALSE)
{
    // constructor
    for (Int_t i=0; i<4; i++) fMultEstimatorAvg[i] = 0;
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it,
                                        // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
  DefineOutput(2, AliNormalizationCounter::Class());
  DefineOutput(3, TTree::Class()); // event
  DefineOutput(4, TTree::Class()); // Lc
  DefineOutput(5, TTree::Class()); // Lc MCGen
  DefineOutput(6, TList::Class()); // Lc weight of MC pt shape
  DefineOutput(7, TTree::Class()); // Lc QA

}
//_____________________________________________________________________________
AliAnalysisTaskSELc2pKs0fromKFP::~AliAnalysisTaskSELc2pKs0fromKFP()
{
    // destructor
    if (fOutputList) {
      delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
      fOutputList = 0;
    }

    if (fOutputWeight) {
      delete fOutputWeight;     // at the end of your task, it is deleted from memory by calling this function
      fOutputWeight = 0;
    }
    if (fPID) {
      delete fPID;
      fPID = 0;
    }

    if (fPIDCombined) {
       delete fPIDCombined;
       fPIDCombined = 0;
    }


    if (fListCuts) {
      delete fListCuts;
      fListCuts = 0;
    }

    if (fAnaCuts) {
      delete fAnaCuts;
      fAnaCuts = 0;
    }

    if (fTree_Event) {
      delete fTree_Event;
      fTree_Event = 0;
    }

    if (fVar_Event) {
      delete fVar_Event;
      fVar_Event = 0;
    }

    if (fTree_Lc) {
      delete fTree_Lc;
      fTree_Lc = 0;
    }

    if (fVar_Lc) {
      delete fVar_Lc;
      fVar_Lc = 0;
    }

    if (fTree_Lc_QA) {
      delete fTree_Lc_QA;
      fTree_Lc_QA = 0;
    }

    if (fVar_Lc_QA) {
      delete fVar_Lc_QA;
      fVar_Lc_QA = 0;
    }

    if (fTree_LcMCGen) {
      delete fTree_LcMCGen;
      fTree_LcMCGen = 0;
    }

    if (fVar_LcMCGen) {
      delete fVar_LcMCGen;
      fVar_LcMCGen = 0;
    }

    if (fCounter) {
      delete fCounter;
      fCounter = 0;
    }

    if (fFuncWeightPythia) {
      delete fFuncWeightPythia;
      fFuncWeightPythia = 0;
    }
    if (fFuncWeightFONLL5overLHC13d3) {
      delete fFuncWeightFONLL5overLHC13d3;
      fFuncWeightFONLL5overLHC13d3 = 0;
    }
    if (fFuncWeightFONLL5overLHC13d3Lc) {
      delete fFuncWeightFONLL5overLHC13d3Lc;
      fFuncWeightFONLL5overLHC13d3Lc = 0;
    }
    for(Int_t i=0; i<4; i++) { if(fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i]; }

}
//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::Init()
{
  // Initialization

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsKFP(*fAnaCuts));
  PostData(1, fListCuts);

  return;

}
//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
  fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)

  DefineAnaHist(); // define analysis histograms

    // example of a histogram
  fHistEvents = new TH1F("fHistEvents", "fHistEvents", 18, 0.5, 18.5);
  fHistEvents->GetXaxis()->SetBinLabel(1,"Analyzed events");
  fHistEvents->GetXaxis()->SetBinLabel(2,"AliAODVertex exists");
  fHistEvents->GetXaxis()->SetBinLabel(3,"TriggerOK");
  fHistEvents->GetXaxis()->SetBinLabel(4,"IsEventSelected");
  fHistEvents->GetXaxis()->SetBinLabel(5,"V0 exists");
  fHistEvents->GetXaxis()->SetBinLabel(6,"Cascade exists");

  fHistEvents->GetXaxis()->SetBinLabel(7,"MCarray exists");
  fHistEvents->GetXaxis()->SetBinLabel(8,"MCheader exists");
  fHistEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
  fHistEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
  fHistEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
  fHistEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
  fHistEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
  fHistEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
  fHistEvents->GetXaxis()->SetBinLabel(15,Form("zVtx<=%2.0fcm", fAnaCuts->GetMaxVtxZ()));
  fHistEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
  fHistEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  fHistEvents->GetXaxis()->SetBinLabel(18,Form("zVtxMC<=%2.0fcm",fAnaCuts->GetMaxVtxZ()));
  fHistEvents->GetYaxis()->SetTitle("counts");

  fHTrigger = new TH1F("fHTrigger", "counter", 18, -0.5, 17.5);
  fHTrigger->SetStats(kTRUE);
  fHTrigger->GetXaxis()->SetBinLabel(1,"X1");
  fHTrigger->GetXaxis()->SetBinLabel(2,"kMB");
  fHTrigger->GetXaxis()->SetBinLabel(3,"kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(4,"kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(5,"kINT7");
  fHTrigger->GetXaxis()->SetBinLabel(6,"kEMC7");
  //fHTrigger->GetXaxis()->SetBinLabel(7,"Space");
  fHTrigger->GetXaxis()->SetBinLabel(8,"kMB|kSemiCentral|kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(9,"kINT7|kEMC7");
  fHTrigger->GetXaxis()->SetBinLabel(11,"kMB&kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(12,"kMB&kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(13,"kINT7&kEMC7");

  fOutputList->Add(fHistEvents); // don't forget to add it to the list! the list will be written to file, so if you want
  fOutputList->Add(fHTrigger);

  // Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(2)->GetContainer();
  if(cont) normName = (TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();
  PostData(2, fCounter);
  DefineEvent();
  PostData(3, fTree_Event);  // postdata will notify the analysis manager of changes / updates to the

  DefineTreeLc_Rec();
  PostData(4, fTree_Lc);

  DefineTreeLc_Gen();
  PostData(5, fTree_LcMCGen);

  fOutputWeight = new TList();
  fOutputWeight->SetOwner(kTRUE);
  fHistMCGen_LcPt_weight = new TH1D("fHistMCGen_LcPt_weight", "", 11, 1., 12.);
  f2DHistMCRec_LcPt_weight = new TH2D("f2DHistMCRec_LcPt_weight", "", 11, 1., 12., 495, 0.5, 50);
  fOutputWeight->Add(fHistMCGen_LcPt_weight);
  fOutputWeight->Add(f2DHistMCRec_LcPt_weight);
  PostData(6, fOutputWeight);

  DefineTreeLc_Rec_QA();
  PostData(7, fTree_Lc_QA);


  //initialise AliPIDCombined object for Bayesian PID
  fPIDCombined = new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);


  // weight function from ratio of flat value (1/30) to pythia
  // use to normalise to flat distribution (should lead to flat pT distribution)
  fFuncWeightPythia = new TF1("funcWeightPythia","1./(30. *[0]*x/TMath::Power(1.+(TMath::Power((x/[1]),[3])),[2]))",0.15,30);
  fFuncWeightPythia->SetParameters(0.36609,1.94635,1.40463,2.5);

  // weight function from the ratio of the LHC13d3 MC
  // and FONLL calculations for pp data
  fFuncWeightFONLL5overLHC13d3 = new TF1("funcWeightFONLL5overLHC13d3","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,30.);
  fFuncWeightFONLL5overLHC13d3->SetParameters(2.94999e+00,3.47032e+00,2.81278e+00,2.5,1.93370e-02,3.86865e+00,-1.54113e-01,8.86944e-02,2.56267e-02);

  fFuncWeightFONLL5overLHC13d3Lc = new TF1("funcWeightFONLL5overLHC13d3Lc","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,20.);
  fFuncWeightFONLL5overLHC13d3Lc->SetParameters(5.94428e+01,1.63585e+01,9.65555e+00,6.71944e+00,8.88338e-02,2.40477e+00,-4.88649e-02,-6.78599e-01,-2.10951e-01);




  return;
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::UserExec(Option_t *)
{
  // user exec
  // this function is called once for each event
  // the manager will take care of reading the events from file, and with the static function InputEvent() you
  // have access to the current event.
  // once you return from the UserExec function, the manager will retrieve the next event from the chain

  if (!fInputEvent) { // if the event is empty (getting it failed) skip this event
    AliError("NO EVENT FOUND!");
    return;
  }
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);    // get an event (called AODEvent) from the input file
                                                        // there's another event format (ESD) which works in a similar way
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with aod's

  fHistEvents->Fill(1);


  TClonesArray *arrayLc2pKs0orLpi = NULL;
  if (!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler *aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if (aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayLc2pKs0orLpi = (TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
    }
  } else {
    arrayLc2pKs0orLpi = (TClonesArray*)aodEvent->GetList()->FindObject("CascadesHF");
  }


  //--------------------------------------------------------------
  // First check if the event has magnetic field and proper vertex
  //--------------------------------------------------------------

  fBzkG = (Double_t)aodEvent->GetMagneticField();
  if (TMath::Abs(fBzkG)<0.001) return;
  // Setting magnetic field for KF vertexing
  KFParticle::SetField(fBzkG);

  fpVtx = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!fpVtx) return;
  fHistEvents->Fill(2);

  fCounter->StoreEvent(aodEvent,fAnaCuts,fIsMC);

  /// Recalculate PV with diamond constraint off
  AliVertexerTracks *vertexer = new AliVertexerTracks(aodEvent->GetMagneticField());
  vertexer->SetConstraintOff();
  fpVtxOff = vertexer->FindPrimaryVertex(aodEvent);


  //------------------------------------------------
  // MC analysis setting
  //------------------------------------------------

  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if (fIsMC) {
    fMCEvent = MCEvent(); // get the corresponding MC event fMCEvent
    if (!fMCEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if ( !mcArray ) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    fHistEvents->Fill(7); // number of MC array exist

    // load MC header
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if ( !mcHeader ) {
      AliError("AliAnalysisTaskSELc2pKs0fromKFP::UserExec: MC header branch not found!\n");
      return;
    }
    fHistEvents->Fill(8); // number of MC header exist

    Double_t zMCvtx = mcHeader->GetVtxZ();
    if ( TMath::Abs(zMCvtx) > fAnaCuts->GetMaxVtxZ() ) {
      AliDebug(2,Form("Event rejected: fabs(zVtxMC)=%f > fAnaCuts->GetMaxVtxZ()=%f", zMCvtx, fAnaCuts->GetMaxVtxZ()));
      return;
    } else {
      fHistEvents->Fill(18);
    }
    if ((TMath::Abs(zMCvtx) < fAnaCuts->GetMaxVtxZ()) && (!fAnaCuts->IsEventRejectedDuePhysicsSelection()) && (!fAnaCuts->IsEventRejectedDueToTrigger())) {
      Bool_t selevt = MakeMCAnalysis(mcArray, mcHeader, aodEvent);
      if(!selevt) return;
    }
  }


  //------------------------------------------------
  // Event selection
  //------------------------------------------------
  Bool_t IsTriggerNotOK = fAnaCuts->IsEventRejectedDueToTrigger();
  Bool_t IsPhysSelNotOK = fAnaCuts->IsEventRejectedDuePhysicsSelection();
  Bool_t IsNoVertex = fAnaCuts->IsEventRejectedDueToNotRecoVertex();
  if( !IsTriggerNotOK && !IsPhysSelNotOK && !IsNoVertex && fabs(fpVtx->GetZ())<fAnaCuts->GetMaxVtxZ() ) fHistEvents->Fill(3);

  Bool_t IsEventSelected = fAnaCuts->IsEventSelected(aodEvent);
  if(!IsEventSelected) {
//    cout<<"Why: "<<fAnaCuts->GetWhyRejection()<<endl;
    return;
  }
  fHistEvents->Fill(4);


  Bool_t IsMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB)==(AliVEvent::kMB);
  Bool_t IsSemi = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kSemiCentral)==(AliVEvent::kSemiCentral);
  Bool_t IsCent = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kCentral)==(AliVEvent::kCentral);
  Bool_t IsINT7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7)==(AliVEvent::kINT7);
  Bool_t IsEMC7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMC7)==(AliVEvent::kEMC7);
  if(IsMB) fHTrigger->Fill(1);
  if(IsSemi) fHTrigger->Fill(2);
  if(IsCent) fHTrigger->Fill(3);
  if(IsINT7) fHTrigger->Fill(4);
  if(IsEMC7) fHTrigger->Fill(5);
  if(IsMB||IsSemi||IsCent) fHTrigger->Fill(7);
  if(IsINT7||IsEMC7) fHTrigger->Fill(8);
  if(IsMB&&IsSemi) fHTrigger->Fill(10);
  if(IsMB&&IsCent) fHTrigger->Fill(11);
  if(IsINT7&&IsEMC7) fHTrigger->Fill(12);

//  AliCentrality *cent = aodEvent->GetCentrality();
//  Float_t Centrality = cent->GetCentralityPercentile("V0M");

  // set primary vertex
  KFPVertex pVertex;
  Double_t pos[3],cov[6];
  fpVtx->GetXYZ(pos);
  if ( fabs(pos[2])>10. ) return; // vertex cut on z-axis direction
  fpVtx->GetCovarianceMatrix(cov);

//  if ( !AliVertexingHFUtils::CheckAODvertexCov(fpVtx) ) cout << "Vertex Cov. is wrong!!!" << endl;
  pVertex.SetXYZ((Float_t)pos[0], (Float_t)pos[1], (Float_t)pos[2]);
  Float_t covF[6];
  for (Int_t i=0; i<6; i++) { covF[i] = (Float_t)cov[i]; }
  pVertex.SetCovarianceMatrix(covF);
  pVertex.SetChi2(fpVtx->GetChi2());
  pVertex.SetNDF(fpVtx->GetNDF());
  pVertex.SetNContributors(fpVtx->GetNContributors());

  KFParticle PV(pVertex);

  if(!fAnaCuts) return;

  FillEventROOTObjects(aodEvent);

//------------------------------------------------
// Main analysis done in this function
//------------------------------------------------

  fPID = fInputHandler->GetPIDResponse();
  MakeAnaLcFromCascadeHF(arrayLc2pKs0orLpi, aodEvent, mcArray, PV);

  PostData(2, fCounter);
  PostData(3, fTree_Event);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
  PostData(4, fTree_Lc);
  PostData(5, fTree_LcMCGen);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
/*
    TCanvas *c1 = new TCanvas();
    fHistDecayLLambda->Draw();
    TLine *lLPDG = new TLine(7.89, 0, 7.89, 1e10);
    lLPDG->SetLineColor(2);
    lLPDG->Draw();
    TCanvas *c2 = new TCanvas();
    fHistDecayLXiMinus->Draw();
    TLine *lXiPDG = new TLine(4.91, 0, 4.91, 1e10);
    lXiPDG->SetLineColor(2);
    lXiPDG->Draw();
*/
    return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSELc2pKs0fromKFP::MakeMCAnalysis(TClonesArray *mcArray, AliAODMCHeader *mcHeader, AliAODEvent *aodEvent)
{

  // method to fill MC histo: how many Lc --> Ks0 + p are there at MC level
  for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) {
    AliAODMCParticle *mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
    if (!mcPart) {
      AliError("Failed casting particle from MC array!, Skipping particle");
      continue;
    }
    Int_t pdg = mcPart->GetPdgCode();
    if ( TMath::Abs(pdg)!=4122 ) {
      AliDebug(2, Form("MC particle %d is not a Lc: its pdg code is %d", iPart, pdg));
      continue;
    }
    AliDebug(2, Form("Step 0 ok: MC particle %d is a Lc: its pdg code is %d", iPart, pdg));
    Int_t labeldaugh0 = mcPart->GetDaughterLabel(0);
    Int_t labeldaugh1 = mcPart->GetDaughterLabel(1);
    if (labeldaugh0 <= 0 || labeldaugh1 <= 0){
      AliDebug(2, Form("The MC particle doesn't have correct daughters, skipping!!"));
      continue;
    }
    else if ( (labeldaugh1 - labeldaugh0) == 1 ) {
      AliDebug(2, Form("Step 1 ok: The MC particle has correct daughters!!"));
      AliAODMCParticle* daugh0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labeldaugh0));
      AliAODMCParticle* daugh1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labeldaugh1));
      if ( !daugh0 || !daugh1 ) {
        AliDebug(2, "Particle daughters not properly retrieved!");
        continue;
      }
      Int_t pdgCodeDaugh0 = daugh0->GetPdgCode();
      Int_t pdgCodeDaugh1 = daugh1->GetPdgCode();
      AliAODMCParticle* bachelorMC = daugh0;
      AliAODMCParticle* v0MC = daugh1;
      if (!fIsAnaLc2Lpi) {
        if ( TMath::Abs(pdgCodeDaugh0)==2212 || TMath::Abs(pdgCodeDaugh1)==2212 ) {
          AliDebug(1, Form("pdgCodeDaugh0 = %d, pdgCodeDaugh1 = %d", pdgCodeDaugh0, pdgCodeDaugh1));
        }
        if ( (TMath::Abs(pdgCodeDaugh0)==311 && TMath::Abs(pdgCodeDaugh1)==2212) || (TMath::Abs(pdgCodeDaugh0)==2212 && TMath::Abs(pdgCodeDaugh1)==311) ) {
          if ( TMath::Abs(pdgCodeDaugh0)==311 ) {
            bachelorMC = daugh1;
            v0MC = daugh0;
          }
          if ( v0MC->GetNDaughters()!=1 ) {
            AliDebug(2, "The K0 does not decay in 1 body only! Impossible... Continuing...");
            continue;
          }
          AliDebug(2, "Step 2 ok: The K0 does decay in 1 body only! ");
          Int_t labelK0daugh = v0MC->GetDaughterLabel(0);
          AliAODMCParticle* partK0S = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelK0daugh));
          if ( !partK0S ) {
            AliError("Error while casting particle! returning a NULL array");
            continue;
          }
//        cout << "=== Lc" << pdg << endl;
//        cout << "=== pdgCodeDaugh0 " << pdgCodeDaugh0 << endl;
//        cout << "=== pdgCodeDaugh1 " << pdgCodeDaugh1 << endl;
//        cout << "=== K0 " << partK0S->GetPdgCode() << endl;
          if ( TMath::Abs(partK0S->GetPdgCode()) != 310 ) {
            AliDebug(2, "The K0 daughter is not a K0S");
            continue;
          }
          //=== check K0S decay ===
          if ( partK0S->GetNDaughters() != 2 ) {
            AliDebug(2, "The K0S does not decay in 2 bodies");
            continue;
          }
          AliDebug(2, "Step 3 ok: The K0 daughter is a K0S and does decay in 2 bodies");
          Int_t labelK0Sdaugh0 = partK0S->GetDaughterLabel(0);
          Int_t labelK0Sdaugh1 = partK0S->GetDaughterLabel(1);
          AliAODMCParticle* daughK0S0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelK0Sdaugh0));
          AliAODMCParticle* daughK0S1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelK0Sdaugh1));
          if ( !daughK0S0 || !daughK0S1 ) {
            AliDebug(2, "Could not access K0S daughters, continuing...");
            continue;
          }
          AliDebug(2, "Step 4 ok: Could access K0S daughters, continuing...");
          Int_t pdgK0Sdaugh0 = daughK0S0->GetPdgCode();
          Int_t pdgK0Sdaugh1 = daughK0S1->GetPdgCode();
          if ( TMath::Abs(pdgK0Sdaugh0)!=211 || TMath::Abs(pdgK0Sdaugh1)!=211 ) {
            AliDebug(2, "The K0S does not decay in pi+pi-, continuing");
            continue;
          }
          if ( TMath::Abs(mcPart->Y()) < 0.8 ) {
            Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
            FillTreeGenLc(mcArray, mcPart, CheckOrigin, mcHeader, aodEvent);
          }
        }
      }
      if (fIsAnaLc2Lpi) {
        if ( TMath::Abs(pdgCodeDaugh0)==3122 || TMath::Abs(pdgCodeDaugh1)==3122 ) {
          AliDebug(1, Form("pdgCodeDaugh0 = %d, pdgCodeDaugh1 = %d", pdgCodeDaugh0, pdgCodeDaugh1));
        }
        if ( (TMath::Abs(pdgCodeDaugh0)==211 && TMath::Abs(pdgCodeDaugh1)==3122) || (TMath::Abs(pdgCodeDaugh0)==3122 && TMath::Abs(pdgCodeDaugh1)==211) ) {
          if ( TMath::Abs(pdgCodeDaugh0)==3122 ) {
            bachelorMC = daugh1;
            v0MC = daugh0;
          }
          if ( v0MC->GetNDaughters()!=2 ) {
            AliDebug(2, "The Lambda does not decay in 2 bodies");
            continue;
          }
          AliDebug(2, "Step 2 ok: The Lambda does decay in 2 bodies");
          Int_t labelLamDaugh0 = v0MC->GetDaughterLabel(0);
          Int_t labelLamDaugh1 = v0MC->GetDaughterLabel(1);
          AliAODMCParticle* daughLam0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelLamDaugh0));
          AliAODMCParticle* daughLam1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labelLamDaugh1));
          if ( !daughLam0 || !daughLam1 ) {
            AliDebug(2, "Could not access Lambda daughters, continuing...");
            continue;
          }
          AliDebug(2, "Step 3 ok: Could access Lambda daughters, continuing...");
          Int_t pdgLamDaugh0 = daughLam0->GetPdgCode();
          Int_t pdgLamDaugh1 = daughLam1->GetPdgCode();
          if ( (TMath::Abs(pdgLamDaugh0)==2212 && TMath::Abs(pdgLamDaugh1)==211) || (TMath::Abs(pdgLamDaugh0)==211 && TMath::Abs(pdgLamDaugh1)==2212) ) {
            if ( TMath::Abs(mcPart->Y()) < 0.8 ) {
              Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
              FillTreeGenLc(mcArray, mcPart, CheckOrigin, mcHeader, aodEvent);
            }
          }
        }
      }
    }
  }

  return kTRUE;

}

//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::FillTreeGenLc(TClonesArray *mcArray,AliAODMCParticle *mcpart, Int_t CheckOrigin, AliAODMCHeader *mcHeader, AliAODEvent *aodEvent)
{
  // Fill histograms or tree depending

  if(!mcpart) return;

  for(Int_t i=0;i<13;i++){
    fVar_LcMCGen[i] = -9999.;
  }

  fVar_LcMCGen[ 0] = fCentrality;
//  if (mcpart->IsPrimary() && (!mcpart->IsPhysicalPrimary())) fVar_LcMCGen[1] = 1;
//  if (mcpart->IsPhysicalPrimary()) fVar_LcMCGen[1] = 2;
//  if (mcpart->IsSecondaryFromWeakDecay()) fVar_LcMCGen[1] = 3;
//  if (mcpart->IsSecondaryFromMaterial()) fVar_LcMCGen[1] = 4;
//  if (mcpart->IsFromSubsidiaryEvent()) fVar_LcMCGen[1] = 5;
  fVar_LcMCGen[ 1] = mcpart->Y();
  fVar_LcMCGen[ 2] = mcpart->Pt();
  fVar_LcMCGen[ 3] = AliVertexingHFUtils::GetBeautyMotherPt(mcArray, mcpart);
  fVar_LcMCGen[ 4] = CheckOrigin;
  fVar_LcMCGen[5] = mcpart->Xv();
  fVar_LcMCGen[6] = mcpart->Yv();
  fVar_LcMCGen[7] = mcpart->Zv();
  if (fUseWeights && CheckOrigin>=0) { //add branches for MC pT weights
    fVar_LcMCGen[8] = fFuncWeightPythia->Eval(mcpart->Pt()); // weight pT flat
    fVar_LcMCGen[9] = fFuncWeightFONLL5overLHC13d3->Eval(mcpart->Pt()); // weight pT flat
    fVar_LcMCGen[10] = fFuncWeightFONLL5overLHC13d3Lc->Eval(mcpart->Pt()); // weight pT flat

  }
  if (fUseMult) { // add multiplicity branches for MC gen
    Double_t zPrimVertex = mcHeader->GetVtxZ();
    TProfile *estimatorAvg = GetEstimatorHistogram(aodEvent);
    Double_t nTrackletsEta10 = static_cast<Double_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.));
    Double_t nTrackletsEta10Corr = static_cast<Double_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,nTrackletsEta10,zPrimVertex,fRefMult));

    fVar_LcMCGen[10] = nTrackletsEta10;
    fVar_LcMCGen[11] = nTrackletsEta10Corr;

  }

  if (fWriteLcMCGenTree) fTree_LcMCGen->Fill();

}

//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::MakeAnaLcFromCascadeHF(TClonesArray *arrayLc2pKs0orLpi, AliAODEvent *aodEvent, TClonesArray *mcArray, KFParticle PV)
{
  // Main analysis called from "UserExec"

//  std::cout.setf(std::ios::fixed);
//  std::cout.setf(std::ios::showpoint);
//  std::cout.precision(3);

  UInt_t nCasc = arrayLc2pKs0orLpi->GetEntriesFast();
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  Double_t covP[21], covN[21], covB[21];
  const Int_t NDaughters = 2;
  const Float_t massKs0_PDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  //const Float_t massLc_PDG  = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  const Float_t massLam_PDG = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  for (UInt_t iCasc=0; iCasc<nCasc; iCasc++) {
    // Lc candidates
    AliAODRecoCascadeHF *Lc2pKs0orLpi = dynamic_cast<AliAODRecoCascadeHF*>(arrayLc2pKs0orLpi->At(iCasc));
    if (!Lc2pKs0orLpi) {
      AliDebug(2,Form("Cascade %d doens't exist, skipping",iCasc));
      continue;
    }

    // --- Needed for newer samples ---
    if (!vHF->FillRecoCasc(aodEvent, Lc2pKs0orLpi, kFALSE)) { // Fill data members of candidate if not done
      AliDebug(2,Form("Cascade %d not refilled, skipping",iCasc));
      continue;
    }
    // --------------------------------

    if (!(Lc2pKs0orLpi->CheckCascadeFlags())) {
      AliDebug(2,Form("Cascade %d is not flagged as Lc candidate",iCasc));
      continue;
    }

    if (Lc2pKs0orLpi->GetNDaughters()!=2) {
      AliDebug(2,Form("Cascade %d does not have 2 daughters (nDaughters=%d)",iCasc,Lc2pKs0orLpi->GetNDaughters()));
      continue;
    }

    AliAODv0 *v0part = dynamic_cast<AliAODv0*>(Lc2pKs0orLpi->Getv0());
    AliAODTrack *bachPart = dynamic_cast<AliAODTrack*>(Lc2pKs0orLpi->GetBachelor());
    if (!v0part || !bachPart) {
      AliDebug(2,Form("Cascade %d has no V0 or no bachelor object",iCasc));
      continue;
    }

    // primary track cuts
//    if ( !fAnaCuts->PassedTrackQualityCuts_Primary(bachPart) ) continue;

    if (!v0part->GetSecondaryVtx()) {
      AliDebug(2,Form("No secondary vertex for V0 by cascade %d",iCasc));
      continue;
    }

    if (v0part->GetNDaughters()!=2) {
      AliDebug(2,Form("current V0 does not have 2 daughters (onTheFly=%d, nDaughters=%d)",v0part->GetOnFlyStatus(),v0part->GetNDaughters()));
      continue;
    }

    if (v0part->GetOnFlyStatus() && !fUseOnTheFlyV0) {
      AliDebug(2,Form("V0 for cascade %d is on-the-fly but only offline requested",iCasc));
      continue;
    }



    AliAODTrack * v0Pos = dynamic_cast<AliAODTrack*>(Lc2pKs0orLpi->Getv0PositiveTrack());
    AliAODTrack * v0Neg = dynamic_cast<AliAODTrack*>(Lc2pKs0orLpi->Getv0NegativeTrack());
    if (!v0Pos || !v0Neg) {
      AliDebug(2,Form("V0 by cascade %d has no V0positive or V0negative object",iCasc));
      continue;
    }

    // check charge of v0
    Int_t v0charge = v0Pos->Charge() + v0Neg->Charge();
    if (v0charge!=0) {
      AliDebug(2,Form("V0 by cascade %d has charge: IMPOSSIBLE!",iCasc));
      continue;
    }

    // check charge of the first daughter, if negative, define it as the second one
    if (v0Pos->Charge()<0) {
      v0Pos = (AliAODTrack*) (Lc2pKs0orLpi->Getv0NegativeTrack());
      v0Neg = (AliAODTrack*) (Lc2pKs0orLpi->Getv0PositiveTrack());
    }

    // === pre-selection for Lc ===
    if (!fIsAnaLc2Lpi) {
      if ( !fAnaCuts->PreSelForLc2pKs0(Lc2pKs0orLpi) ) continue;
      if ( fabs(v0part->MassLambda()-massLam_PDG) <= fAnaCuts->GetProdMassRejLambda() ) continue;
      if ( fabs(v0part->MassAntiLambda()-massLam_PDG) <= fAnaCuts->GetProdMassRejLambda() ) continue;
    }
    if (fIsAnaLc2Lpi) {
      if ( !fAnaCuts->PreSelForLc2Lpi(Lc2pKs0orLpi) ) continue;
      if ( fabs(v0part->MassK0Short()-massKs0_PDG) <= fAnaCuts->GetProdMassRejKs0() ) continue;
    }
//    if ( v0part->InvMass2Prongs(0,1,11,11) <= fAnaCuts->GetMassCutEplusEminus() ) continue; // InvMass of e+e-
    if ( v0part->InvMass2Prongs(0,1,11,11) <= 0.1 ) continue; // InvMass of e+e-
    // ============================

    if ( !bachPart->GetCovarianceXYZPxPyPz(covB) || !v0Pos->GetCovarianceXYZPxPyPz(covP) || !v0Neg->GetCovarianceXYZPxPyPz(covN) ) continue;
    if ( !AliVertexingHFUtils::CheckAODtrackCov(bachPart) || !AliVertexingHFUtils::CheckAODtrackCov(v0Pos) || !AliVertexingHFUtils::CheckAODtrackCov(v0Neg) ) continue;


    //// Recalculate primary vertex without daughter tracks, if requested
    /// IMPORTANT: Own primary vertex must be unset before continue/return, else memory leak
    bool recVtx = false;
    AliAODVertex *origOwnVtx = 0x0;
    AliAODVertex *ownPVtx = 0x0; // recalculated primary vertex for candidate
    if (fAnaCuts->GetIsPrimaryWithoutDaughters())
    {
      if (Lc2pKs0orLpi->GetOwnPrimaryVtx() )
               origOwnVtx = new AliAODVertex(*Lc2pKs0orLpi->GetOwnPrimaryVtx());
      if (fAnaCuts->RecalcOwnPrimaryVtx(Lc2pKs0orLpi,aodEvent)) {
         ownPVtx = Lc2pKs0orLpi->GetOwnPrimaryVtx();
         recVtx = true;
         KFPVertex pVertex;
         Double_t pos[3],cov[6];
         ownPVtx->GetXYZ(pos);
         if ( fabs(pos[2])>10. ) {Lc2pKs0orLpi->UnsetOwnPrimaryVtx(); fAnaCuts->CleanOwnPrimaryVtx(Lc2pKs0orLpi,aodEvent,origOwnVtx); continue;} // vertex cut on z-axis direction
         ownPVtx->GetCovarianceMatrix(cov);
         //  if ( !AliVertexingHFUtils::CheckAODvertexCov(fpVtx) ) cout << "Vertex Cov. is wrong!!!" << endl;
         pVertex.SetXYZ((Float_t)pos[0], (Float_t)pos[1], (Float_t)pos[2]);
         Float_t covF[6];
         for (Int_t i=0; i<6; i++) { covF[i] = (Float_t)cov[i]; }
         pVertex.SetCovarianceMatrix(covF);
         pVertex.SetChi2(ownPVtx->GetChi2());
         pVertex.SetNDF(ownPVtx->GetNDF());
         pVertex.SetNContributors(ownPVtx->GetNContributors());
         PV = KFParticle(pVertex);
      }  else {
         Lc2pKs0orLpi->UnsetOwnPrimaryVtx();
         fAnaCuts->CleanOwnPrimaryVtx(Lc2pKs0orLpi,aodEvent,origOwnVtx);
         continue;
      }
    }

    KFParticle kfpBach;
    Bool_t isRej = kFALSE;

    if (!fIsAnaLc2Lpi) {
      if (bachPart->Charge()>0) kfpBach = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart, 2212); // proton
      if (bachPart->Charge()<0) kfpBach = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart, -2212); // anti-proton
      KFParticle kfpPionPlus   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, 211);
      KFParticle kfpPionMinus  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, -211);
      KFParticle kfpKs0;
      const KFParticle *Ks0Daughters[2] = {&kfpPionPlus, &kfpPionMinus};
      kfpKs0.Construct(Ks0Daughters, NDaughters);
      Float_t massKs0_rec=0., err_massKs0_rec=0.;
      kfpKs0.GetMass(massKs0_rec, err_massKs0_rec);
      // check rapidity of Ks0
      if ( TMath::Abs(kfpKs0.GetE())<=TMath::Abs(kfpKs0.GetPz()) ) isRej = kTRUE;

      // chi2>0 && NDF>0 for selecting Ks0
      if ( (kfpKs0.GetNDF()<=1.e-10 || kfpKs0.GetChi2()<=1.e-10) ) isRej = kTRUE;

      // check cov. of Ks0
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpKs0) ) isRej = kTRUE;

      // err_mass(Ks0) > 0
      if ( err_massKs0_rec<=1.e-10 ) isRej = kTRUE;

      // Chi2geo cut of Ks0
      if ( (kfpKs0.GetChi2()/kfpKs0.GetNDF()) >= fAnaCuts->GetKFPKs0_Chi2geoMax() ) isRej = kTRUE;


      // Calculate l/Δl for Ks0
      Double_t ldl_Ks0 = AliVertexingHFUtils::ldlFromKF(kfpKs0, PV);

      // l/Deltal cut of Ks0
      if ( ldl_Ks0 <= fAnaCuts->GetKFPKs0_lDeltalMin() ) isRej = kTRUE;

      // mass window cut of Ks0
      if ( TMath::Abs(massKs0_rec-massKs0_PDG) > (fAnaCuts->GetProdMassTolKs0()) ) isRej = kTRUE;

      // mass constraint for Ks0
      KFParticle kfpKs0_massConstraint = kfpKs0;
      kfpKs0_massConstraint.SetNonlinearMassConstraint(massKs0_PDG);

      // QA check after mass constraint
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpKs0_massConstraint) || TMath::Abs(kfpKs0_massConstraint.GetE()) <= TMath::Abs(kfpKs0_massConstraint.GetPz()) ) isRej = kTRUE;

      // reconstruct Lc with mass constraint
      KFParticle kfpLc;
      const KFParticle *LcDaughters[2] = {&kfpBach, &kfpKs0_massConstraint};
      kfpLc.Construct(LcDaughters, NDaughters);

      // reconstruct Lc without mass constraint
      KFParticle kfpLc_woKs0MassConst;
      const KFParticle *LcDaughters_woKs0MassConst[2] = {&kfpBach, &kfpKs0};
      kfpLc_woKs0MassConst.Construct(LcDaughters_woKs0MassConst, NDaughters);

      // === for Lc with mass constraint ===
      // check rapidity of Lc
      if ( TMath::Abs(kfpLc.GetE())<=TMath::Abs(kfpLc.GetPz()) ) isRej = kTRUE;

      // chi2>0 && NDF>0
      if ( kfpLc.GetNDF()<=1.e-10 || kfpLc.GetChi2()<=1.e-10 ) isRej = kTRUE;

      // check covariance matrix
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLc) ) isRej = kTRUE;

      // Prefilter
      if ( kfpLc.GetChi2()/kfpLc.GetNDF() >= fAnaCuts->GetKFPLc_Chi2geoMax() ) isRej = kTRUE;
      if ( kfpLc.GetPt() < fAnaCuts->GetPtMinLc() ) isRej = kTRUE;

      // err_mass(Lc) > 0
      Float_t massLc_rec=0., err_massLc_rec=0.;
      kfpLc.GetMass(massLc_rec, err_massLc_rec);
      if (err_massLc_rec <= 1.e-10 ) isRej = kTRUE;
      // ===================================

      // === for Lc without mass constraint ===
      // check rapidity of Lc
      if ( TMath::Abs(kfpLc_woKs0MassConst.GetE())<=TMath::Abs(kfpLc_woKs0MassConst.GetPz()) ) isRej = kTRUE;

      // chi2>0 && NDF>0
      if ( kfpLc_woKs0MassConst.GetNDF()<=1.e-10 || kfpLc_woKs0MassConst.GetChi2()<=1.e-10 ) isRej = kTRUE;

      // check covariance matrix
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLc_woKs0MassConst) ) isRej = kTRUE;

      // Prefilter
      if ( kfpLc_woKs0MassConst.GetChi2()/kfpLc_woKs0MassConst.GetNDF() >= fAnaCuts->GetKFPLc_Chi2geoMax() ) isRej = kTRUE;
      if ( kfpLc_woKs0MassConst.GetPt() < fAnaCuts->GetPtMinLc() ) isRej = kTRUE;

      // err_mass(Lc) > 0
      Float_t massLc_woKs0MassConst_rec=0., err_massLc_woKs0MassConst_rec=0.;
      kfpLc_woKs0MassConst.GetMass(massLc_woKs0MassConst_rec, err_massLc_woKs0MassConst_rec);
      if (err_massLc_woKs0MassConst_rec <= 1.e-10 ) isRej = kTRUE;
      // ===================================
      if (isRej) {
         if (recVtx) {
            Lc2pKs0orLpi->UnsetOwnPrimaryVtx();
            fAnaCuts->CleanOwnPrimaryVtx(Lc2pKs0orLpi,aodEvent,origOwnVtx);
         }
         continue;
      }
      if (fWriteLcTree) {
        Int_t lab_Lc  = -9999;
        Int_t lab_Ks0 = -9999;
        if ( fIsMC ) {
          lab_Ks0 = MatchToMCKs0(v0Pos, v0Neg, mcArray);
          lab_Lc  = MatchToMCLc2pKs0(v0Pos, v0Neg, bachPart, mcArray);
        }
        FillTreeRecLcFromCascadeHF(Lc2pKs0orLpi, kfpLc, bachPart, kfpBach, kfpKs0, kfpKs0_massConstraint, v0Pos, v0Neg, PV, mcArray, lab_Ks0, lab_Lc, kfpLc_woKs0MassConst, aodEvent, ownPVtx);
      }
      kfpLc.Clear();
      kfpKs0_massConstraint.Clear();
      kfpKs0.Clear();
      kfpPionMinus.Clear();
      kfpPionPlus.Clear();
    }

    if (fIsAnaLc2Lpi) {
      if (bachPart->Charge()>0) {
        kfpBach = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart, 211); // pion+
        KFParticle kfpProton    = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, 2212);
        KFParticle kfpPionMinus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, -211);
        KFParticle kfpLam;
        const KFParticle *LamDaughters[2] = {&kfpProton, &kfpPionMinus};
        kfpLam.Construct(LamDaughters, 2);
        Float_t massLam_rec=0., err_massLam_rec=0.;
        kfpLam.GetMass(massLam_rec, err_massLam_rec);

        // check rapidity of Lam
        if ( TMath::Abs(kfpLam.GetE())<=TMath::Abs(kfpLam.GetPz()) ) isRej = kTRUE;

        // chi2>0 && NDF>0 for selecting Lam
        if ( (kfpLam.GetNDF()<=1.e-10 || kfpLam.GetChi2()<=1.e-10) ) isRej = kTRUE;

        // check cov. of Lam
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLam) ) isRej = kTRUE;

        // err_mass(Lam) > 0
        if ( err_massLam_rec<=1.e-10 ) isRej = kTRUE;

        // Chi2geo cut of Lam
        if ( (kfpLam.GetChi2()/kfpLam.GetNDF()) >= fAnaCuts->GetKFPLam_Chi2geoMax() ) isRej = kTRUE;

        // Calculate l/Δl for Lam
        Double_t ldl_Lam = AliVertexingHFUtils::ldlFromKF(kfpLam, PV);

        // l/Deltal cut of Lam
        if ( ldl_Lam <= fAnaCuts->GetKFPLam_lDeltalMin() ) isRej = kTRUE;

        // mass window cut of Lam
        if ( TMath::Abs(massLam_rec-massLam_PDG) > (fAnaCuts->GetProdMassTolLambda()) ) isRej = kTRUE;

        // mass constraint for Lam
        KFParticle kfpLam_massConstraint = kfpLam;
        kfpLam_massConstraint.SetNonlinearMassConstraint(massLam_PDG);

        // QA check after mass constraint
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLam_massConstraint) || TMath::Abs(kfpLam_massConstraint.GetE()) <= TMath::Abs(kfpLam_massConstraint.GetPz()) ) isRej = kTRUE;

        // reconstruct Lc with mass constraint
        KFParticle kfpLc;
        const KFParticle *LcDaughters[2] = {&kfpBach, &kfpLam_massConstraint};
        kfpLc.Construct(LcDaughters, 2);

        // reconstruct Lc without mass constraint
        KFParticle kfpLc_woLamMassConst;
        const KFParticle *LcDaughters_woLamMassConst[2] = {&kfpBach, &kfpLam};
        kfpLc_woLamMassConst.Construct(LcDaughters_woLamMassConst, 2);

        // === for Lc with mass constraint ===
        // check rapidity of Lc
        if ( TMath::Abs(kfpLc.GetE())<=TMath::Abs(kfpLc.GetPz()) ) isRej = kTRUE;

        // chi2>0 && NDF>0
        if ( kfpLc.GetNDF()<=1.e-10 || kfpLc.GetChi2()<=1.e-10 ) isRej = kTRUE;

        // check covariance matrix
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLc) ) isRej = kTRUE;

        // Prefilter
        if ( kfpLc.GetChi2()/kfpLc.GetNDF() >= fAnaCuts->GetKFPLc_Chi2geoMax() ) isRej = kTRUE;
        if ( kfpLc.GetPt() < fAnaCuts->GetPtMinLc() ) isRej = kTRUE;

        // err_mass(Lc) > 0
        Float_t massLc_rec=0., err_massLc_rec=0.;
        kfpLc.GetMass(massLc_rec, err_massLc_rec);
        if (err_massLc_rec <= 1.e-10 ) isRej = kTRUE;
        // ===================================

        // === for Lc without mass constraint ===
        // check rapidity of Lc
        if ( TMath::Abs(kfpLc_woLamMassConst.GetE())<=TMath::Abs(kfpLc_woLamMassConst.GetPz()) ) isRej = kTRUE;

        // chi2>0 && NDF>0
        if ( kfpLc_woLamMassConst.GetNDF()<=1.e-10 || kfpLc_woLamMassConst.GetChi2()<=1.e-10 ) isRej = kTRUE;

        // check covariance matrix
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLc_woLamMassConst) ) isRej = kTRUE;

        // Prefilter
        if ( kfpLc_woLamMassConst.GetChi2()/kfpLc_woLamMassConst.GetNDF() >= fAnaCuts->GetKFPLc_Chi2geoMax() ) isRej = kTRUE;
        if ( kfpLc_woLamMassConst.GetPt() < fAnaCuts->GetPtMinLc() ) isRej = kTRUE;

        // err_mass(Lc) > 0
        Float_t massLc_woLamMassConst_rec=0., err_massLc_woLamMassConst_rec=0.;
        kfpLc_woLamMassConst.GetMass(massLc_woLamMassConst_rec, err_massLc_woLamMassConst_rec);
        if (err_massLc_woLamMassConst_rec <= 1.e-10 ) isRej = kTRUE;
        // ===================================
        if (isRej) {
           if (recVtx) {
              Lc2pKs0orLpi->UnsetOwnPrimaryVtx();
              fAnaCuts->CleanOwnPrimaryVtx(Lc2pKs0orLpi,aodEvent,origOwnVtx);
           }
           continue;
        }

        if (fWriteLcTree) {
          Int_t lab_Lc  = -9999;
          Int_t lab_Lam = -9999;
          if ( fIsMC ) {
            lab_Lam = MatchToMCLam(v0Pos, v0Neg, mcArray, kTRUE);
            lab_Lc  = MatchToMCLc2Lpi(v0Pos, v0Neg, bachPart, mcArray, kTRUE);
          }
          FillTreeRecLcFromCascadeHF(Lc2pKs0orLpi, kfpLc, bachPart, kfpBach, kfpLam, kfpLam_massConstraint, v0Pos, v0Neg, PV, mcArray, lab_Lam, lab_Lc, kfpLc_woLamMassConst, aodEvent, ownPVtx);
        }
        kfpLc_woLamMassConst.Clear();
        kfpLc.Clear();
        kfpLam_massConstraint.Clear();
        kfpLam.Clear();
        kfpPionMinus.Clear();
        kfpProton.Clear();
      }

      if (bachPart->Charge()<0) {
        kfpBach = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart, -211); // pion-
        KFParticle kfpAntiProton = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, -2212);
        KFParticle kfpPionPlus   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, 211);
        KFParticle kfpAntiLam;
        const KFParticle *AntiLamDaughters[2] = {&kfpAntiProton, &kfpPionPlus};
        kfpAntiLam.Construct(AntiLamDaughters, 2);
        Float_t massAntiLam_rec=0., err_massAntiLam_rec=0.;
        kfpAntiLam.GetMass(massAntiLam_rec, err_massAntiLam_rec);

        // check rapidity of AntiLam
        if ( TMath::Abs(kfpAntiLam.GetE())<=TMath::Abs(kfpAntiLam.GetPz()) ) isRej = kTRUE;

        // chi2>0 && NDF>0 for selecting AntiLam
        if ( (kfpAntiLam.GetNDF()<=1.e-10 || kfpAntiLam.GetChi2()<=1.e-10) ) isRej = kTRUE;

        // check cov. of AntiLam
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLam) ) isRej = kTRUE;

        // err_mass(AntiLam) > 0
        if ( err_massAntiLam_rec<=1.e-10 ) isRej = kTRUE;

        // Chi2geo cut of AntiLam
        if ( (kfpAntiLam.GetChi2()/kfpAntiLam.GetNDF()) >= fAnaCuts->GetKFPLam_Chi2geoMax() ) isRej = kTRUE;

        // Calculate l/Δl for AntiLam
        Double_t ldl_AntiLam = AliVertexingHFUtils::ldlFromKF(kfpAntiLam, PV);

        // l/Deltal cut of AntiLam
        if ( ldl_AntiLam <= fAnaCuts->GetKFPLam_lDeltalMin() ) isRej = kTRUE;

        // mass window cut of AntiLam
        if ( TMath::Abs(massAntiLam_rec-massLam_PDG) > (fAnaCuts->GetProdMassTolLambda()) ) isRej = kTRUE;

        // mass constraint for AntiLam
        KFParticle kfpAntiLam_massConstraint = kfpAntiLam;
        kfpAntiLam_massConstraint.SetNonlinearMassConstraint(massLam_PDG);

        // QA check after mass constraint
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLam_massConstraint) || TMath::Abs(kfpAntiLam_massConstraint.GetE()) <= TMath::Abs(kfpAntiLam_massConstraint.GetPz()) ) isRej = kTRUE;

        // reconstruct AntiLc with mass constraint
        KFParticle kfpAntiLc;
        const KFParticle *AntiLcDaughters[2] = {&kfpBach, &kfpAntiLam_massConstraint};
        kfpAntiLc.Construct(AntiLcDaughters, 2);

        // reconstruct AntiLc without mass constraint
        KFParticle kfpAntiLc_woAntiLamMassConst;
        const KFParticle *AntiLcDaughters_woAntiLamMassConst[2] = {&kfpBach, &kfpAntiLam};
        kfpAntiLc_woAntiLamMassConst.Construct(AntiLcDaughters_woAntiLamMassConst, 2);

        // === for AntiLc with mass constraint ===
        // check rapidity of AntiLc
        if ( TMath::Abs(kfpAntiLc.GetE())<=TMath::Abs(kfpAntiLc.GetPz()) ) isRej = kTRUE;

        // chi2>0 && NDF>0
        if ( kfpAntiLc.GetNDF()<=1.e-10 || kfpAntiLc.GetChi2()<=1.e-10 ) isRej = kTRUE;

        // check covariance matrix
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLc) ) isRej = kTRUE;

        // Prefilter
        if ( kfpAntiLc.GetChi2()/kfpAntiLc.GetNDF() >= fAnaCuts->GetKFPLc_Chi2geoMax() ) isRej = kTRUE;
        if ( kfpAntiLc.GetPt() < fAnaCuts->GetPtMinLc() ) isRej = kTRUE;

        // err_mass(AntiLc) > 0
        Float_t massAntiLc_rec=0., err_massAntiLc_rec=0.;
        kfpAntiLc.GetMass(massAntiLc_rec, err_massAntiLc_rec);
        if (err_massAntiLc_rec <= 1.e-10 ) isRej = kTRUE;
        // ===================================

        // === for AntiLc without mass constraint ===
        // check rapidity of AntiLc
        if ( TMath::Abs(kfpAntiLc_woAntiLamMassConst.GetE())<=TMath::Abs(kfpAntiLc_woAntiLamMassConst.GetPz()) ) isRej = kTRUE;

        // chi2>0 && NDF>0
        if ( kfpAntiLc_woAntiLamMassConst.GetNDF()<=1.e-10 || kfpAntiLc_woAntiLamMassConst.GetChi2()<=1.e-10 ) isRej = kTRUE;

        // check covariance matrix
        if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLc_woAntiLamMassConst) ) isRej = kTRUE;

        // Prefilter
        if ( kfpAntiLc_woAntiLamMassConst.GetChi2()/kfpAntiLc_woAntiLamMassConst.GetNDF() >= fAnaCuts->GetKFPLc_Chi2geoMax() ) isRej = kTRUE;
        if ( kfpAntiLc_woAntiLamMassConst.GetPt() < fAnaCuts->GetPtMinLc() ) isRej = kTRUE;

        // err_mass(AntiLc) > 0
        Float_t massAntiLc_woAntiLamMassConst_rec=0., err_massAntiLc_woAntiLamMassConst_rec=0.;
        kfpAntiLc_woAntiLamMassConst.GetMass(massAntiLc_woAntiLamMassConst_rec, err_massAntiLc_woAntiLamMassConst_rec);
        if (err_massAntiLc_woAntiLamMassConst_rec <= 1.e-10 ) isRej = kTRUE;
        // ===================================

        if (isRej) {
           if (recVtx) {
              Lc2pKs0orLpi->UnsetOwnPrimaryVtx();
              fAnaCuts->CleanOwnPrimaryVtx(Lc2pKs0orLpi,aodEvent,origOwnVtx);
           }
           continue;
        }
        if (fWriteLcTree) {
          Int_t lab_AntiLc  = -9999;
          Int_t lab_AntiLam = -9999;
          if ( fIsMC ) {
            lab_AntiLam = MatchToMCLam(v0Pos, v0Neg, mcArray, kFALSE);
            lab_AntiLc  = MatchToMCLc2Lpi(v0Pos, v0Neg, bachPart, mcArray, kFALSE);
          }
          FillTreeRecLcFromCascadeHF(Lc2pKs0orLpi, kfpAntiLc, bachPart, kfpBach, kfpAntiLam, kfpAntiLam_massConstraint, v0Pos, v0Neg, PV, mcArray, lab_AntiLam, lab_AntiLc, kfpAntiLc_woAntiLamMassConst, aodEvent, ownPVtx);
        }
        kfpAntiLc_woAntiLamMassConst.Clear();
        kfpAntiLc.Clear();
        kfpAntiLam_massConstraint.Clear();
        kfpAntiLam.Clear();
        kfpPionPlus.Clear();
        kfpAntiProton.Clear();
      }
    }
    kfpBach.Clear();
    if (recVtx) {
       Lc2pKs0orLpi->UnsetOwnPrimaryVtx();
       fAnaCuts->CleanOwnPrimaryVtx(Lc2pKs0orLpi,aodEvent,origOwnVtx);
    }
  }

  delete vHF;
  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::SelectTrack(AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks, Bool_t *seleFlags)
{
  // Select good tracks using fAnaCuts (AliRDHFCuts object)
  if(trkEntries==0) return;

  nSeleTrks=0;
  for(Int_t i=0; i<trkEntries; i++) {
    seleFlags[i] = kFALSE;

    AliVTrack *track;
    track = (AliVTrack*)event->GetTrack(i);

//    if(track->GetID()<0) continue;
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;

//    AliAODTrack *aodt = (AliAODTrack*)track;

/*
    if(!fAnaCuts) continue;
    if(fAnaCuts->SingleTrkCuts(aodt)){
      seleFlags[i]=kTRUE;
      nSeleTrks++;
//      fHistoPiPtRef->Fill(aodt->Pt());
    }
*/
  } // end loop on tracks
}

//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::DefineEvent()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(3)->GetContainer()->GetName();
  fTree_Event = new TTree(nameoutput, "Event");
  Int_t nVar = 13;
  fVar_Event = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  fVarNames[0]  = "centrality";
  fVarNames[1]  = "x_vtx_reco";
  fVarNames[2]  = "y_vtx_reco";
  fVarNames[3]  = "z_vtx_reco";
  fVarNames[4]  = "n_vtx_contributors";
  fVarNames[5]  = "n_tracks";
  fVarNames[6]  = "is_ev_sel";
  fVarNames[7]  = "run_number";
  fVarNames[8]  = "ev_id";
  fVarNames[9]  = "x_vtx_reco_constOff";
  fVarNames[10]  = "y_vtx_reco_constOff";
  fVarNames[11]  = "z_vtx_reco_constOff";
  fVarNames[12]  = "n_vtx_contributors_constOff";


  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_Event->Branch(fVarNames[ivar].Data(), &fVar_Event[ivar], Form("%s/F", fVarNames[ivar].Data()));
  }

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::DefineTreeLc_Rec()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fTree_Lc = new TTree(nameoutput, "Lc variables tree");
  Int_t nVar = 0;
  TString *fVarNames;
  if (!fIsAnaLc2Lpi){
    if (!fKeepAllVariables){
      if (fIsMC){
        nVar = 37;
        fVar_Lc = new Float_t[nVar];
        fVarNames = new TString[nVar];
      }
      else{
        nVar = 32;
        fVar_Lc = new Float_t[nVar];
        fVarNames = new TString[nVar];
      }
    }
    if (fKeepAllVariables){
      nVar = 50;
      fVar_Lc = new Float_t[nVar];
      fVarNames = new TString[nVar];
    }
  }
  if (fIsAnaLc2Lpi) {
    nVar = 41;
    fVar_Lc = new Float_t[nVar];
    fVarNames = new TString[nVar];
  }


  if (!fIsAnaLc2Lpi) {
    fVarNames[0]  = "nSigmaTPC_PiPlus"; //TPC nsigma for pion+ coming from K0s
    fVarNames[1]  = "nSigmaTPC_PiMinus"; //TPC nsigma for pion- coming from K0s
    fVarNames[2]  = "nSigmaTPC_Pr"; //TPC nsigma for proton
    fVarNames[3]  = "nSigmaTOF_PiPlus"; //TOF nsigma for pion+ coming from K0s
    fVarNames[4]  = "nSigmaTOF_PiMinus"; //TOF nsigma for pion coming from K0s
    fVarNames[5]  = "nSigmaTOF_Pr"; //TOF nsigma for proton
    fVarNames[6]  = "chi2topo_Ks0_PV"; //chi2_topological of K0s (with mass constraint) to PV
    fVarNames[7]  = "ldl_Ks0"; //l/dl of K0s
    fVarNames[8]  = "chi2topo_Lc"; //chi2_topological of Lc (with mass const. of Ks0) to PV
    fVarNames[9]  = "ldl_Lc"; //l/dl of Lc (with mass const. of Ks0)
    fVarNames[10] = "ldlxy_Lc"; //l/dl of Lc (with mass const. of Ks0)
    fVarNames[11] = "DecayLxy_Ks0"; //decay length of Ks0 in x-y plane
    fVarNames[12] = "ct_Ks0"; // life time of Ks0
    fVarNames[13] = "DCA_Ks0"; // DCA of Ks0 to PV
    fVarNames[14] = "DecayLxy_Lc"; //decay length of Lc in x-y plane
    fVarNames[15] = "DecayL_Lc"; //decay length of Lc
    fVarNames[16] = "DCA_Lc"; // DCA of Lc to PV
    fVarNames[17] = "PA_Ks0"; //pointing angle of Ks0 (pointing back to Lc)
    fVarNames[18] = "PA_Lc"; //pointing angle of Lc (pointing back to PV)
    fVarNames[19] = "cosPAxy_Lc"; //pointing angle of Lc in xy plane (pointing back to PV)
    fVarNames[20] = "d0xy_LcToPV"; // impact parameter of Lc w.r.t PV
    fVarNames[21] = "pt_Lc"; //pt of Lc (with mass const. of Ks0 and with PV const.)
    fVarNames[22] = "rap_Lc"; //rapidity of Lc (with mass const. of Ks0 and with PV const.)
    fVarNames[23] = "mass_Lc"; //mass of Lc (with mass const. of Ks0 and with PV const.)
    fVarNames[24] = "pt_Pr"; //pt of proton
    fVarNames[25] = "d0_PrToPV"; //rphi impact params of proton w.r.t. Primary Vtx [cm]
    fVarNames[26] = "signd0_p"; // signed d0 of proton
    fVarNames[27] = "nTrackletsCorr"; // corrected Ntrk
    fVarNames[28] = "CombinedPIDProb_Pr"; // Bayesian PID probability of proton for bachelor track
    fVarNames[29] = "CombinedPIDProb_Pr_TPCOnly"; // Bayesian PID probability of proton for bachelor track
    fVarNames[30] = "nSigmaCombined_Pr"; // nSigma-combined for proton
    fVarNames[31] = "nSigmaCombined_Pi_bach"; // nSigma-combined for proton from pions (for exclusion)
    // fVarNames[25] = "AODVertex_X"; //Primary vertex from AOD X position
    // fVarNames[26] = "AODVertex_Y"; //Primary vertex from AOD Y position
    // fVarNames[27] = "AODVertex_Z"; //Primary vertex from AOD Z position
    // fVarNames[28] = "AODVertex_X_pRemoved"; //Primary vertex from AOD X position after proton removal
    // fVarNames[29] = "AODVertex_Y_pRemoved"; //Primary vertex from AOD Y position after proton removal
    // fVarNames[30] = "AODVertex_Z_pRemoved"; //Primary vertex from AOD Z position after proton removal
    if (fIsMC && !fKeepAllVariables) {
      ///Only needed in MC
      fVarNames[32] = "Source_Lc"; //flag for Lc MC truth (“>=0” signal, “<0” background)
      fVarNames[33] = "pt_B"; //pt of B hadron
      fVarNames[34] = "weightPtFlat"; // flat pT weight for MC
      fVarNames[35] = "weightFONLL5overLHC13d3"; // FONLL / LHC13d3 weight (default D meson)
      fVarNames[36] = "weightFONLL5overLHC13d3Lc"; // FONLL/LHC13d3 weight (modified for baryon)
    }
    if (fKeepAllVariables){
      fVarNames[32] = "Source_Lc"; //flag for Lc MC truth (“>=0” signal, “<0” background)
      fVarNames[33] = "pt_B"; //pt of B hadron
      fVarNames[34] = "weightPtFlat"; // flat pT weight for MC
      fVarNames[35] = "weightFONLL5overLHC13d3"; // FONLL / LHC13d3 weight (default D meson)
      fVarNames[36] = "weightFONLL5overLHC13d3Lc"; // FONLL/LHC13d3 weight (modified for baryon)
      /// Additional variables
      fVarNames[37]  = "DCA_Ks0Dau"; //Distance between pions coming from K0s (calculated from AOD v0)
      fVarNames[38]  = "chi2geo_Ks0"; //chi2_geometry of K0s (without mass constraint)
      fVarNames[39] = "chi2geo_Lc"; //chi2_geometry of Lc (with mass const. of Ks0)
      fVarNames[40] = "pt_Ks0"; //pt of Ks0 (without mass const.)
      fVarNames[41] = "mass_Ks0"; //mass of Ks0 (without mass const.)
      fVarNames[42] = "pt_PiPlus"; //pt of pion+
      fVarNames[43] = "pt_PiMinus"; //pt of pion-
      fVarNames[44] = "d0_Ks0ToPV"; //rphi impact params of Ks0 w.r.t. Primary Vtx [cm]
      fVarNames[45] = "cosThetaStar"; //cos-thetastar of decay
      fVarNames[46] = "armenteros_K0s"; // armenteros qT/|alpha| for cascade
      fVarNames[47] = "cos_p_K0s";   // cos pointing angle of V0 from RecoCascadeHF
      fVarNames[48] = "d_len_K0s";    // decay length of V0 from RecoCascadeHF
      fVarNames[49] = "nTrackletsRaw"; // raw Ntrk

    }


  }
  if (fIsAnaLc2Lpi) {
    fVarNames[0]  = "nSigmaTPC_V0Pr"; //TPC nsigma for proton coming from Lam
    fVarNames[1]  = "nSigmaTPC_V0Pi"; //TPC nsigma for pion coming from Lam
    fVarNames[2]  = "nSigmaTPC_Bach"; //TPC nsigma for bachlor
    fVarNames[3]  = "nSigmaTOF_V0Pr"; //TOF nsigma for proton coming from Lam
    fVarNames[4]  = "nSigmaTOF_V0Pi"; //TOF nsigma for pion coming from Lam
    fVarNames[5]  = "nSigmaTOF_Bach"; //TOF nsigma for bachlor

    fVarNames[6]  = "DCA_LamDau"; //Distance between daughters coming from Lam (calculated from AOD v0)

    fVarNames[7]  = "chi2geo_Lam"; //chi2_geometry of Lam (without mass constraint)
    fVarNames[8]  = "chi2topo_Lam_PV"; //chi2_topological of Lam (with mass constraint) to PV
    fVarNames[9]  = "ldl_Lam"; //l/dl of Lam

    fVarNames[10] = "chi2geo_Lc"; //chi2_geometry of Lc (with mass const. of Lam)
    fVarNames[11] = "chi2topo_Lc"; //chi2_topological of Lc (with mass const. of Lam) to PV
    fVarNames[12] = "ldl_Lc"; //l/dl of Lc (with mass const. of Lam)

    fVarNames[13] = "DecayLxy_Lam"; //decay length of Lam in x-y plane
    fVarNames[14] = "ct_Lam"; // life time of Lam
    fVarNames[15] = "DecayLxy_Lc"; //decay length of Lc in x-y plane
    fVarNames[16] = "PA_Lam"; //pointing angle of Lam (pointing back to Lc)
    fVarNames[17] = "PA_Lc"; //pointing angle of Lc (pointing back to PV)

    fVarNames[18] = "pt_Lam"; //pt of Lam (without mass const.)
    fVarNames[19] = "mass_Lam"; //mass of Lam (without mass const.)
    fVarNames[20] = "pt_Lc"; //pt of Lc (with mass const. of Lam and with PV const.)
    fVarNames[21] = "rap_Lc"; //rapidity of Lc (with mass const. of Lam and with PV const.)
    fVarNames[22] = "mass_Lc"; //mass of Lc (with mass const. of Lam and with PV const.)
    fVarNames[23] = "pt_Bach"; //pt of bachlor
    fVarNames[24] = "pt_V0Pr"; //pt of proton
    fVarNames[25] = "pt_V0Pi"; //pt of pion
    fVarNames[26] = "d0_BachToPV"; //rphi impact params of bachlor w.r.t. Primary Vtx [cm]
    fVarNames[27] = "d0_LamToPV"; //rphi impact params of Lam w.r.t. Primary Vtx [cm]
    fVarNames[28] = "Source_Lc"; //flag for Lc MC truth (“>=0” signal, “<0” background)
    fVarNames[29] = "cosThetaStar"; //cos-thetastar of decay
    fVarNames[30] = "CombinedPIDProb_V0Pr"; // Bayesian PID probability of proton from Lam decay
    fVarNames[31] = "armenteros_Lam"; // armenteros qT/|alpha| for cascade
    fVarNames[32] = "nSigmaCombined_V0Pr"; // nSigma-combined for proton from Lam decay
    fVarNames[33] = "nSigmaCombined_Pi_V0Pr"; // nSigma-combined for pion from Lam decay (for exclusion)
    fVarNames[34] = "cos_p_Lam"; // cosine pointing angle of cascade
    fVarNames[35] = "d_len_Lam"; // dlen of cascade
    fVarNames[36] = "weightPtFlat"; // flat pT weight for MC
    fVarNames[37] = "weightFONLL5overLHC13d3"; // FONLL / LHC13d3 weight (default D meson)
    fVarNames[38] = "weightFONLL5overLHC13d3Lc"; // FONLL/LHC13d3 weight (modified for baryon)
    fVarNames[39] = "nTrackletsRaw"; // raw Ntrk
    fVarNames[40] = "nTrackletsCorr"; // corrected Ntrk


  }

//  fVarNames[]  = "chi2geo_Ks0_wMassConst"; //chi2_geometry of K0s (with mass constraint)
//  fVarNames[] = "DecayL_Ks0"; //decay length of K0s in 3D
//  fVarNames[] = "DecayL_Lc"; //decay length of Lc in 3D
//  fVarNames[] = "Source_Ks0"; //flag for Ks0 MC truth (“>=0” signal, “<0” background)
//  fVarNames[] = "Source_Lc"; //flag for Lc MC truth (“>=0” signal, “<0” background)
//  fVarNames[] = "DCA_PrToPV_KF"; //DCA of proton to PV from KF in 3D
//  fVarNames[] = "ct_Lc"; // life time of Lc

//  fVarNames[] = "mass_Lam"; //mass of Lambda
//  fVarNames[] = "mass_AntiLam"; //mass of Anti-Lambda
//  fVarNames[] = "mass_Gamma"; //mass of e+e-

//  fVarNames[] = "NtrkCorr"; //number of tracks after correction

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_Lc->Branch(fVarNames[ivar].Data(), &fVar_Lc[ivar], Form("%s/F", fVarNames[ivar].Data()));
  }

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::DefineTreeLc_Rec_QA()
{

  const char* nameoutput = GetOutputSlot(7)->GetContainer()->GetName();
  fTree_Lc_QA = new TTree(nameoutput, "QA of Lc variables tree");
  Int_t nVar = 26;
  fVar_Lc_QA = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  if (!fIsAnaLc2Lpi) {
    // check pt, rapidity and mass with and without mass constraint
    // pt and mass without mass constraint in analysis tree
    fVarNames[0] = "rap_Ks0_woMassConst"; //rapidity of Ks0 (without mass const.)
    fVarNames[1] = "pt_Ks0_wMassConst"; //pt of Ks0 (with mass const.)
    fVarNames[2] = "rap_Ks0_wMassConst"; //rapidity of Ks0 (with mass const.)
    fVarNames[3] = "mass_Ks0_wMassConst"; //mass of Ks0 (with mass const.)
    fVarNames[4] = "chi2mass_Ks0"; //chi2_MassConst. of Ks0

    // check chi2_topo of Ks0 pointing to PV without Ks0 mass const.
    fVarNames[5] = "chi2topo_Ks0woMassConst_PV"; //chi2_topo of Ks0 (without mass const.) to PV

    // check Lc reconstruction without Ks0 mass const.
    fVarNames[6] = "chi2geo_Lc_woKs0MassConst"; //chi2_geometry of Lc (without mass constraint of Ks0)
    fVarNames[7] = "pt_Lc_woKs0MassConst"; //pt of Lc reconstructed without Ks0 mass const.
    fVarNames[8] = "rap_Lc_woKs0MassConst"; //rap of Lc reconstructed without Ks0 mass const.
    fVarNames[9] = "mass_Lc_woKs0MassConst"; //mass of Lc reconstructed without Ks0 mass const. (case 1 requested by Silvia)

    // check Lc reconstruction with Ks0 mass const. and without PV const.
    fVarNames[10] = "pt_Lc_woPV"; //pt of Lc (with mass const. of Ks0 and without PV const.)
    fVarNames[11] = "rap_Lc_woPV"; //rap of Lc (with mass const. of Ks0 and without PV const.)
    fVarNames[12] = "mass_Lc_woPV"; //mass of Lc (with mass const. of Ks0 and without PV const.) (case 2 requested by Silvia)

    // --- reconstruction logic 2 (set production vertex of Ks0 to PV) ---
    fVarNames[13] = "pt_Ks0_PV"; //pt of Ks0 (with mass and PV const.)
    fVarNames[14] = "rap_Ks0_PV"; //rapidity of Ks0 (with mass and PV const.)
    fVarNames[15] = "mass_Ks0_PV"; //mass of Ks0 (with mass and PV const.)

    fVarNames[16] = "chi2geo_Lc_wKs0MassPVConst"; //chi2_geo of Lc (with mass and PV const. of Ks0)
    fVarNames[17] = "pt_Lc_wKs0MassPVConst"; //pt of Lc (with mass and PV const. of Ks0)
    fVarNames[18] = "rap_Lc_wKs0MassPVConst"; //rapidity of Lc (with mass and PV const. of Ks0)
    fVarNames[19] = "mass_Lc_wKs0MassPVConst"; //mass of Lc (with mass and PV const. of Ks0)
    fVarNames[20] = "ldl_Lc_wKs0MassPVConst"; //l/dl of Lc (with mass and PV const. of Ks0)

    fVarNames[21] = "chi2topo_Lc_wKs0MassPVConst"; //chi2_topo of Lc (with mass and PV const. of Ks0)
    // -------------------------------------------------------------------

    // flags for signal and background
    fVarNames[22] = "Source_Ks0"; //flag for Ks0 MC truth (“>=0” signal, “<0” background)
    fVarNames[23] = "Source_Lc"; //flag for Lc MC truth (“>=0” signal, “<0” background)
    fVarNames[24] = "mass_Lc_AODRecoCascadeHF"; //mass of Lc (from AliAODRecoCascadeHF)
    fVarNames[25] = "pt_Lc_AODRecoCascadeHF"; //pt of Lc (from AliAODRecoCascadeHF)
  }
  if (fIsAnaLc2Lpi) {
    // check pt, rapidity and mass with and without mass constraint
    // pt and mass without mass constraint in analysis tree
    fVarNames[0] = "rap_Lam_woMassConst"; //rapidity of Lam (without mass const.)
    fVarNames[1] = "pt_Lam_wMassConst"; //pt of Lam (with mass const.)
    fVarNames[2] = "rap_Lam_wMassConst"; //rapidity of Lam (with mass const.)
    fVarNames[3] = "mass_Lam_wMassConst"; //mass of Lam (with mass const.)
    fVarNames[4] = "chi2mass_Lam"; //chi2_MassConst. of Lam

    // check chi2_topo of Lam pointing to PV without Lam mass const.
    fVarNames[5] = "chi2topo_LamwoMassConst_PV"; //chi2_topo of Lam (without mass const.) to PV

    // check Lc reconstruction without Lam mass const.
    fVarNames[6] = "chi2geo_Lc_woLamMassConst"; //chi2_geometry of Lc (without mass constraint of Lam)
    fVarNames[7] = "pt_Lc_woLamMassConst"; //pt of Lc reconstructed without Lam mass const.
    fVarNames[8] = "rap_Lc_woLamMassConst"; //rap of Lc reconstructed without Lam mass const.
    fVarNames[9] = "mass_Lc_woLamMassConst"; //mass of Lc reconstructed without Lam mass const. (case 1 requested by Silvia)

    // check Lc reconstruction with Lam mass const. and without PV const.
    fVarNames[10] = "pt_Lc_woPV"; //pt of Lc (with mass const. of Lam and without PV const.)
    fVarNames[11] = "rap_Lc_woPV"; //rap of Lc (with mass const. of Lam and without PV const.)
    fVarNames[12] = "mass_Lc_woPV"; //mass of Lc (with mass const. of Lam and without PV const.) (case 2 requested by Silvia)

    // --- reconstruction logic 2 (set production vertex of Lam to PV) ---
    fVarNames[13] = "pt_Lam_PV"; //pt of Lam (with mass and PV const.)
    fVarNames[14] = "rap_Lam_PV"; //rapidity of Lam (with mass and PV const.)
    fVarNames[15] = "mass_Lam_PV"; //mass of Lam (with mass and PV const.)

    fVarNames[16] = "chi2geo_Lc_wLamMassPVConst"; //chi2_geo of Lc (with mass and PV const. of Lam)
    fVarNames[17] = "pt_Lc_wLamMassPVConst"; //pt of Lc (with mass and PV const. of Lam)
    fVarNames[18] = "rap_Lc_wLamMassPVConst"; //rapidity of Lc (with mass and PV const. of Lam)
    fVarNames[19] = "mass_Lc_wLamMassPVConst"; //mass of Lc (with mass and PV const. of Lam)
    fVarNames[20] = "ldl_Lc_wLamMassPVConst"; //l/dl of Lc (with mass and PV const. of Lam)

    fVarNames[21] = "chi2topo_Lc_wLamMassPVConst"; //chi2_topo of Lc (with mass and PV const. of Lam)
    // -------------------------------------------------------------------

    // flags for signal and background
    fVarNames[22] = "Source_Lam"; //flag for Lam MC truth (“>=0” signal, “<0” background)
    fVarNames[23] = "Source_Lc"; //flag for Lc MC truth (“>=0” signal, “<0” background)
    fVarNames[24] = "mass_Lc_AODRecoCascadeHF"; //mass of Lc (from AliAODRecoCascadeHF)
    fVarNames[25] = "pt_Lc_AODRecoCascadeHF"; //pt of Lc (from AliAODRecoCascadeHF)
  }

//  fVarNames[] = "CosThetaStar_Pr"; //cosine angle between the proton momentum in the Lc rest frame and the boost direction

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_Lc_QA->Branch(fVarNames[ivar].Data(), &fVar_Lc_QA[ivar], Form("%s/F", fVarNames[ivar].Data()));
  }

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::DefineTreeLc_Gen()
{
  const char* nameoutput = GetOutputSlot(5)->GetContainer()->GetName();
  fTree_LcMCGen = new TTree(nameoutput,"Lc MC variables tree");
  Int_t nVar = 13;
  fVar_LcMCGen = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  fVarNames[ 0] = "Centrality";
  fVarNames[ 1] = "LcY";
  fVarNames[ 2] = "LcPt";
  fVarNames[ 3] = "BPt";
  fVarNames[ 4] = "LcSource";
  fVarNames[ 5] = "Vertex_X"; //Primary vertex X position
  fVarNames[ 6] = "Vertex_Y"; //Primary vertex Y position
  fVarNames[ 7] = "Vertex_Z"; //Primary vertex Z position
  fVarNames[ 8] = "weightPtFlat"; // flat pT weight for MC
  fVarNames[ 9] = "weightFONLL5overLHC13d3"; // FONLL / LHC13d3 weight (default D meson)
  fVarNames[ 10] = "weightFONLL5overLHC13d3Lc"; // FONLL/LHC13d3 weight (modified for baryon)
  fVarNames[11] = "nTrackletsRaw"; // raw Ntrk
  fVarNames[12] = "nTrackletsCorr"; // corrected Ntrk

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_LcMCGen->Branch(fVarNames[ivar].Data(),&fVar_LcMCGen[ivar],Form("%s/F",fVarNames[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskSELc2pKs0fromKFP::InvMassV0atPV(AliAODTrack *trk1, AliAODTrack *trk2, Int_t pdg1, Int_t pdg2)
{

  Double_t mass1 = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
  Double_t mass2 = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();
  Double_t E1 = TMath::Sqrt(mass1*mass1 + trk1->P()*trk1->P());
  Double_t E2 = TMath::Sqrt(mass2*mass2 + trk2->P()*trk2->P());
  Double_t mass = TMath::Sqrt( (E1+E2)*(E1+E2) - (trk1->Px()+trk2->Px())*(trk1->Px()+trk2->Px()) - (trk1->Py()+trk2->Py())*(trk1->Py()+trk2->Py()) - (trk1->Pz()+trk2->Pz())*(trk1->Pz()+trk2->Pz()) );

  return mass;
}

//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::DefineAnaHist()
{
  // Define analysis histograms
}

//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::FillEventROOTObjects(AliAODEvent* aodEvent)
{

  for (Int_t i=0; i<9; i++) {
    fVar_Event[i] = 0.;
  }

  Double_t pos[3];
  fpVtx->GetXYZ(pos);


  fVar_Event[1] = pos[0];
  fVar_Event[2] = pos[1];
  fVar_Event[3] = pos[2];
  fVar_Event[4] = fpVtx->GetNContributors();
  fVar_Event[5] = aodEvent->GetNumberOfTracks();
  fVar_Event[7] = aodEvent->GetRunNumber();
  AliAODHeader *header = dynamic_cast<AliAODHeader*>(aodEvent->GetHeader());
  ULong64_t eventId = header->GetEventIdAsLong();
  fVar_Event[8] = eventId;
  fVar_Event[9] = fpVtxOff->GetX();
  fVar_Event[10] = fpVtxOff->GetY();
  fVar_Event[11] = fpVtxOff->GetZ();
  fVar_Event[12] = fpVtxOff->GetNContributors();

  fTree_Event->Fill();

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSELc2pKs0fromKFP::FillTreeRecLcFromCascadeHF(AliAODRecoCascadeHF *Lc2pKs0orLpi, KFParticle kfpLc, AliAODTrack *trackBach, KFParticle kfpBach, KFParticle kfpV0, KFParticle kfpV0_massConstraint, AliAODTrack *v0Pos, AliAODTrack *v0Neg, KFParticle PV, TClonesArray *mcArray, Int_t lab_V0, Int_t lab_Lc, KFParticle kfpLc_woV0MassConst, AliAODEvent *aodEvent, AliAODVertex *ownPVtx)
{

  if (!fIsAnaLc2Lpi){
    if (!fKeepAllVariables){
      if (fIsMC){
        for (Int_t i=0; i<37; i++) {
          fVar_Lc[i] = -9999.;
        }
      }
      else{
        for (Int_t i=0; i<32; i++) {
          fVar_Lc[i] = -9999.;
        }
      }
    }
    if (fKeepAllVariables){
      for (Int_t i=0; i<50; i++) {
        fVar_Lc[i] = -9999.;
      }
    }
  }
  if (fIsAnaLc2Lpi) {
    for (Int_t i=0; i<41; i++) {
      fVar_Lc[i] = -9999.;
    }
  }

  for (Int_t i=0; i<26; i++) {
    fVar_Lc_QA[i] = -9999.;
  }
  /// KF particle Lambda_c
  KFParticle kfpLc_PV = kfpLc;
  kfpLc_PV.SetProductionVertex(PV);

  /// pt cut for Lc
  Float_t pT_Lc=0.;
  pT_Lc = kfpLc_PV.GetPt();
  if ( pT_Lc <= fAnaCuts->GetPtMinLc() ) return;

  /// mass window cut for Lc
  const Float_t massLc_PDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  Float_t massLc_rec=0., err_massLc_rec=0.;
  kfpLc_PV.GetMass(massLc_rec, err_massLc_rec);
  if ( fabs(massLc_rec-massLc_PDG) > fAnaCuts->GetProdMassTolLc() ) return;

  /// Partice Identification
  Float_t nSigmaTPC_v0Pos = 0.;
  Float_t nSigmaTPC_v0Neg = 0.;
  Float_t nSigmaTPC_v0Pos_excl = 0.;
  Float_t nSigmaTPC_v0Neg_excl = 0.;
  Float_t nSigmaTPC_bach  = 0.;
  Float_t nSigmaTPC_bach_pi  = 0.;
  Float_t nSigmaTOF_v0Pos = 0.;
  Float_t nSigmaTOF_v0Neg = 0.;
  Float_t nSigmaTOF_v0Pos_excl = 0.;
  Float_t nSigmaTOF_v0Neg_excl = 0.;
  Float_t nSigmaTOF_bach  = 0.;
  Float_t nSigmaTOF_bach_pi  = 0.;

  if (!fIsAnaLc2Lpi) {
    nSigmaTPC_v0Pos = fPID->NumberOfSigmasTPC(v0Pos, AliPID::kPion);
    nSigmaTPC_v0Neg = fPID->NumberOfSigmasTPC(v0Neg, AliPID::kPion);
    nSigmaTPC_bach  = fPID->NumberOfSigmasTPC(trackBach, AliPID::kProton);
    nSigmaTPC_bach_pi = fPID->NumberOfSigmasTPC(trackBach, AliPID::kPion);

    nSigmaTOF_v0Pos = fPID->NumberOfSigmasTOF(v0Pos, AliPID::kPion);
    nSigmaTOF_v0Neg = fPID->NumberOfSigmasTOF(v0Neg, AliPID::kPion);
    nSigmaTOF_bach  = fPID->NumberOfSigmasTOF(trackBach, AliPID::kProton);
    nSigmaTOF_bach_pi  = fPID->NumberOfSigmasTOF(trackBach, AliPID::kPion);
  }
  if (fIsAnaLc2Lpi) {
    nSigmaTPC_bach  = fPID->NumberOfSigmasTPC(trackBach, AliPID::kPion);
    nSigmaTOF_bach  = fPID->NumberOfSigmasTOF(trackBach, AliPID::kPion);
    if (trackBach->Charge()>0) {
      nSigmaTPC_v0Pos = fPID->NumberOfSigmasTPC(v0Pos, AliPID::kProton);
      nSigmaTPC_v0Neg = fPID->NumberOfSigmasTPC(v0Neg, AliPID::kPion);
      nSigmaTPC_v0Pos_excl = fPID->NumberOfSigmasTPC(v0Pos, AliPID::kPion);
      nSigmaTPC_v0Neg_excl = fPID->NumberOfSigmasTPC(v0Neg, AliPID::kProton);

      nSigmaTOF_v0Pos = fPID->NumberOfSigmasTOF(v0Pos, AliPID::kProton);
      nSigmaTOF_v0Neg = fPID->NumberOfSigmasTOF(v0Neg, AliPID::kPion);
      nSigmaTOF_v0Pos_excl = fPID->NumberOfSigmasTOF(v0Pos, AliPID::kPion);
      nSigmaTOF_v0Neg_excl = fPID->NumberOfSigmasTOF(v0Neg, AliPID::kProton);
    }
    if (trackBach->Charge()<0) {
      nSigmaTPC_v0Pos = fPID->NumberOfSigmasTPC(v0Pos, AliPID::kPion);
      nSigmaTPC_v0Neg = fPID->NumberOfSigmasTPC(v0Neg, AliPID::kProton);
      nSigmaTPC_v0Pos_excl = fPID->NumberOfSigmasTPC(v0Pos, AliPID::kProton);
      nSigmaTPC_v0Neg_excl = fPID->NumberOfSigmasTPC(v0Neg, AliPID::kPion);

      nSigmaTOF_v0Pos = fPID->NumberOfSigmasTOF(v0Pos, AliPID::kPion);
      nSigmaTOF_v0Neg = fPID->NumberOfSigmasTOF(v0Neg, AliPID::kProton);
      nSigmaTOF_v0Pos_excl = fPID->NumberOfSigmasTOF(v0Pos, AliPID::kProton);
      nSigmaTOF_v0Neg_excl = fPID->NumberOfSigmasTOF(v0Neg, AliPID::kPion);

    }
  }
  /// apply 4 sigma cut for TPC
  if ( fabs(nSigmaTPC_v0Pos)>=4. || fabs(nSigmaTPC_v0Neg)>=4. || fabs(nSigmaTPC_bach)>=4. ) return;

  /// KF particle V0 constrained with Lc vertex as production vertex
  KFParticle kfpV0_Lc = kfpV0_massConstraint;
  kfpV0_Lc.SetProductionVertex(kfpLc);
//  if ( kfpV0_Lc.GetChi2()/kfpV0_Lc.GetNDF() >= fAnaCuts->GetKFPKs0_Chi2topoMax() ) return;

  /// KF particle V0 constrained with primary vertex as production vertex
  KFParticle kfpV0_PV = kfpV0_massConstraint;
  kfpV0_PV.SetProductionVertex(PV);
//  if ( kfpV0_PV.GetChi2()/kfpV0_PV.GetNDF() >= fAnaCuts->GetKFPKs0_Chi2topoMax() ) return;

  /*
  // === mass of Lam ===
  KFParticle kfpProton     = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, 2212);
  KFParticle kfpPionMinus  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, -211);
  KFParticle kfpLam;
  const KFParticle *vLamDaughters[2] = {&kfpProton, &kfpPionMinus};
  kfpLam.Construct(vLamDaughters, 2);
  Float_t massLam_rec, err_massLam;
  kfpLam.GetMass(massLam_rec, err_massLam);
  // === mass of AntiLam ===
  KFParticle kfpAntiProton = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, -2212);
  KFParticle kfpPionPlus   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, 211);
  KFParticle kfpAntiLam;
  const KFParticle *vAntiLamDaughters[2] = {&kfpAntiProton, &kfpPionPlus};
  kfpAntiLam.Construct(vAntiLamDaughters, 2);
  Float_t massAntiLam_rec, err_massAntiLam;
  kfpAntiLam.GetMass(massAntiLam_rec, err_massAntiLam);
  // === mass of Gamma ===
  KFParticle kfpElePlus    = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, -11);
  KFParticle kfpEleMinus   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, 11);
  KFParticle kfpGamma;
  const KFParticle *vGammaDaughters[2] = {&kfpElePlus, &kfpEleMinus};
  kfpGamma.Construct(vGammaDaughters, 2);
  Float_t massGamma_rec, err_massGamma;
  kfpGamma.GetMass(massGamma_rec, err_massGamma);
  */

  /// KF particle bachelor
  KFParticle kfpBach_Lc = kfpBach;
  kfpBach_Lc.SetProductionVertex(kfpLc);

  /// calculate CosPointingAngle
  Double_t cosPA_V0 = AliVertexingHFUtils::CosPointingAngleFromKF(kfpV0_massConstraint, kfpLc);
  Double_t cosPA_Lc  = AliVertexingHFUtils::CosPointingAngleFromKF(kfpLc, PV);

  AliAODv0 *v0 = dynamic_cast<AliAODv0*>(Lc2pKs0orLpi->Getv0());

  /// Decay length, lifetime V0 candidate
  Float_t DecayLxy_V0=0., err_DecayLxy_V0=0.;
  kfpV0_Lc.GetDecayLengthXY(DecayLxy_V0, err_DecayLxy_V0);
  Float_t ct_V0=0., err_ct_V0=0.;
  kfpV0_Lc.GetLifeTime(ct_V0, err_ct_V0);

  /// Decay length, lifetime Lc candidate
  Float_t DecayLxy_Lc=0., err_DecayLxy_Lc=0.;
  kfpLc_PV.GetDecayLengthXY(DecayLxy_Lc, err_DecayLxy_Lc);
  Float_t DecayL_Lc=0., err_DecayL_Lc=0.;
  kfpLc_PV.GetDecayLength(DecayL_Lc, err_DecayL_Lc);

  
  Double_t pos[3], cov[6];
  fpVtx->GetXYZ(pos);
  fpVtx->GetCovarianceMatrix(cov);
  Float_t posF[3], covF[6];
  for(int iEl = 0; iEl < 3; iEl++)
    posF[iEl] = (float)pos[iEl];
  for(int iEl = 0; iEl < 6; iEl++)
    covF[iEl] = (float)cov[iEl];

  Float_t dcaPointV0[8], dcaPointV0Cov[36];
  kfpV0_Lc.GetParametersAtPoint(posF, covF, dcaPointV0, dcaPointV0Cov);
  Float_t dcaV02 = 0;
  for(int i = 0; i < 3; i++){
      dcaV02 += (dcaPointV0[i] - pos[i]) * (dcaPointV0[i] - pos[i]);
  }
  Float_t DCA_V0 = TMath::Sqrt(dcaV02);

  Float_t dcaPoint[8], dcaPointCov[36];
  kfpLc_PV.GetParametersAtPoint(posF, covF, dcaPoint, dcaPointCov);
  Float_t dca2 = 0;
  for(int i = 0; i < 3; i++){
    dca2 += (dcaPoint[i] - pos[i]) * (dcaPoint[i] - pos[i]);
  }
  Float_t DCA_Lc = TMath::Sqrt(dca2);

  Float_t mass_V0_rec=0., err_mass_V0_rec=0.;
  kfpV0.GetMass(mass_V0_rec, err_mass_V0_rec);
  //  fVar_Lc[] = kfpPr.GetDistanceFromVertex(PV); //DCA of proton to PV from KF in 3D

  // Sign of d0 proton (different from regular d0)
  double d0z0bach[2], covd0z0bach[3];
  trackBach->PropagateToDCA(fpVtx, fBzkG, kVeryBig, d0z0bach, covd0z0bach);
  double tx[3];
  trackBach->GetXYZ(tx);
  tx[0] -= fpVtx->GetX();
  tx[1] -= fpVtx->GetY();
  tx[2] -= fpVtx->GetZ();
  double innerpro = tx[0]*kfpLc_PV.Px()+tx[1]*kfpLc_PV.Py();
  double signd0 = 1.;
  if(innerpro<0.) signd0 = -1.;
  signd0 = signd0*TMath::Abs(d0z0bach[0]);

  if ( TMath::Abs(kfpLc_PV.GetE())<TMath::Abs(kfpLc_PV.GetPz()) ) return;

  /// Combined PID response (Bayesian probability)
  Double_t probProton = -1.;
  Double_t probProtonTPC = -1.;
  if (!fIsAnaLc2Lpi) {
    Double_t probTPCTOF[AliPID::kSPECIES] = {-1.};
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
    UInt_t detUsed = fPIDCombined->ComputeProbabilities(trackBach, fPID, probTPCTOF);

    if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()) { //TPC+TOF both present
      probProton = probTPCTOF[AliPID::kProton];
    }
    else {   // if TOF information not available, try only TPC
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
      detUsed = fPIDCombined->ComputeProbabilities(trackBach, fPID, probTPCTOF);
      if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()) { // Check that TPC-only worked. If not, then return -1 as probability
         probProton = probTPCTOF[AliPID::kProton];
      }
      //Reset detector mask for PIDCombined object to TPC+TOF
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
    }
  }
  if (fIsAnaLc2Lpi) {
    // Combined PID response (Bayesian probability) [30]
    Double_t probTPCTOF[AliPID::kSPECIES] = {-1.};
    UInt_t detUsed = 0;
    if (trackBach->Charge()>0) fPIDCombined->ComputeProbabilities(v0Pos, fPID, probTPCTOF);
    if (trackBach->Charge()<0) fPIDCombined->ComputeProbabilities(v0Neg, fPID, probTPCTOF);
    if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()) { //TPC+TOF both present
      probProton = probTPCTOF[AliPID::kProton];
    }
    else {   // if TOF information not available, try only TPC
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
      if (trackBach->Charge()>0) detUsed = fPIDCombined->ComputeProbabilities(v0Pos, fPID, probTPCTOF);
      if (trackBach->Charge()<0) detUsed = fPIDCombined->ComputeProbabilities(v0Neg, fPID, probTPCTOF);
      if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()) { // Check that TPC-only worked. If not, then return -1 as probability
         probProton = probTPCTOF[AliPID::kProton];
      }
      //Reset detector mask for PIDCombined object to TPC+TOF
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
    }
  }

  if (!fIsAnaLc2Lpi) {
    // Combined PID response (Bayesian probability) using only TPC
    Double_t probTPC[AliPID::kSPECIES] = {-1.};
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
    UInt_t detUsed = fPIDCombined->ComputeProbabilities(trackBach, fPID, probTPC);
    if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()) {// Check that TPC-only worked. If not, then return -1 as probability
      probProtonTPC = probTPC[AliPID::kProton];
    }
    //Reset detector mask for PIDCombined object to TPC+TOF
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
  }

  AliAODMCParticle *mcProton;
  AliAODMCParticle *mcLc;
  if (fIsMC && fUseWeights && lab_Lc >= 0) { //add branches for MC pT weights
    Int_t labelProton = fabs(trackBach->GetLabel());
    mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
    Int_t IndexLc = mcProton->GetMother();
    mcLc = static_cast<AliAODMCParticle*>(mcArray->At(IndexLc));
  }

  Double_t nTrackletsEta10 = 0.;
  Double_t nTrackletsEta10Corr = 0.;
  if (fUseMult) {
    AliAODVertex *aodVtx = (AliAODVertex*)aodEvent->GetPrimaryVertex();
    Double_t zPrimVertex = aodVtx->GetZ();
    TProfile *estimatorAvg = GetEstimatorHistogram(aodEvent);
    nTrackletsEta10 = static_cast<Double_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.));
    nTrackletsEta10Corr = static_cast<Double_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,nTrackletsEta10,zPrimVertex,fRefMult));
  }

  if (!fIsAnaLc2Lpi){
    fVar_Lc[0]  = nSigmaTPC_v0Pos;
    fVar_Lc[1]  = nSigmaTPC_v0Neg;
    fVar_Lc[2]  = nSigmaTPC_bach;
    fVar_Lc[3]  = nSigmaTOF_v0Pos;
    fVar_Lc[4]  = nSigmaTOF_v0Neg;
    fVar_Lc[5]  = nSigmaTOF_bach;
    fVar_Lc[6]  = kfpV0_PV.GetChi2()/kfpV0_PV.GetNDF(); //chi2_topological of V0 (with mass constraint) to PV
    fVar_Lc[7]  = AliVertexingHFUtils::ldlFromKF(kfpV0, PV); // ldl_V0
    fVar_Lc[8]  = kfpLc_PV.GetChi2()/kfpLc_PV.GetNDF(); // chi2topo_Lc
    fVar_Lc[9]  = AliVertexingHFUtils::ldlFromKF(kfpLc, PV); // ldl_Lc
    fVar_Lc[10]  = AliVertexingHFUtils::ldlXYFromKF(kfpLc, PV); // ldl_Lc in xy plane
    fVar_Lc[11] = DecayLxy_V0;
    fVar_Lc[12] = ct_V0;
    fVar_Lc[13]   = DCA_V0;
    fVar_Lc[14] = DecayLxy_Lc;
    fVar_Lc[15] = DecayL_Lc;
    fVar_Lc[16]   = DCA_Lc; // DCA of Lc
    fVar_Lc[17] = TMath::ACos(cosPA_V0); // PA_V0
    fVar_Lc[18] = TMath::ACos(cosPA_Lc);  // PA_Lc
    fVar_Lc[19] = AliVertexingHFUtils::CosPointingAngleXYFromKF(kfpLc, PV); // cos PA Lc in xy plane
    fVar_Lc[20] = kfpLc_PV.GetDistanceFromVertexXY(PV); //impact parameter of Lc w.r.t PV
    fVar_Lc[21] = pT_Lc; //pt of Lc (with mass const. of Ks0 and with PV const.)
    fVar_Lc[22] = kfpLc_PV.GetRapidity();
    fVar_Lc[23] = massLc_rec; //mass of Lc (with mass const. of Ks0 and with PV const.)
    fVar_Lc[24] = kfpBach.GetPt();
    fVar_Lc[25] = Lc2pKs0orLpi->Getd0Prong(0); //rphi impact params of bachlor w.r.t. Primary Vtx [cm]
    fVar_Lc[26]   = signd0; // signed d0 proton
    if (fUseMult){
      fVar_Lc[27] = nTrackletsEta10Corr;
    }
    fVar_Lc[28] = probProton;
    fVar_Lc[29] = probProtonTPC;
    fVar_Lc[30] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nSigmaTPC_bach,nSigmaTOF_bach); // nsigma_combined for proton bachelor from Lc
    fVar_Lc[31] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nSigmaTPC_bach_pi,nSigmaTOF_bach_pi);
    // fVar_Lc[25] = fpVtx->GetX();
    // fVar_Lc[26] = fpVtx->GetY();
    // fVar_Lc[27] = fpVtx->GetZ();
    // fVar_Lc[28] = ownPVtx->GetX();
    // fVar_Lc[29] = ownPVtx->GetY();
    // fVar_Lc[30] = ownPVtx->GetZ();
    if (fIsMC && lab_Lc >= 0){
      fVar_Lc[32] = lab_Lc;
      fVar_Lc[33] = AliVertexingHFUtils::GetBeautyMotherPt(mcArray, mcLc);
      if (fUseWeights) {
        fVar_Lc[34] = fFuncWeightPythia->Eval(mcLc->Pt()); // weight pT flat
        fVar_Lc[35] = fFuncWeightFONLL5overLHC13d3->Eval(mcLc->Pt()); // weight pT flat
        fVar_Lc[36] = fFuncWeightFONLL5overLHC13d3Lc->Eval(mcLc->Pt()); // weight pT flat
      }
    }
    if (fKeepAllVariables) {
      fVar_Lc[37]  = v0->GetDCA(); // DCA_V0Dau
      fVar_Lc[38]  = kfpV0.GetChi2()/kfpV0.GetNDF(); //chi2_geometry of V0 (without mass constraint)
      fVar_Lc[39] = kfpLc.GetChi2()/kfpLc.GetNDF(); // chi2geo_Lc
      fVar_Lc[40] = kfpV0.GetPt();
      fVar_Lc[41] = mass_V0_rec;
      fVar_Lc[42] = v0Pos->Pt(); // pion+
      fVar_Lc[43] = v0Neg->Pt(); // pion-
      fVar_Lc[44] = Lc2pKs0orLpi->Getd0Prong(1); ////rphi impact params of V0 w.r.t. Primary Vtx [cm]
      fVar_Lc[45] = AliVertexingHFUtils::CosThetaStarFromKF(0, 4122, 2212, 310, kfpLc, kfpBach_Lc, kfpV0_Lc);  ///cos theta-star
      fVar_Lc[46] = v0->PtArmV0() / TMath::Abs(v0->AlphaV0()); //armenteros qT/|alpha|
      fVar_Lc[47] = cosPA_V0;
      fVar_Lc[48] = AliVertexingHFUtils::DecayLengthFromKF(kfpV0,PV) ;   //d_len_K0s;
      if (fUseMult){
        fVar_Lc[49] = nTrackletsEta10;
      }
    }
  }
  if (fIsAnaLc2Lpi){
    if (trackBach->Charge()>0) {
      fVar_Lc[0]  = nSigmaTPC_v0Pos;
      fVar_Lc[1]  = nSigmaTPC_v0Neg;
      fVar_Lc[2]  = nSigmaTPC_bach;
      fVar_Lc[3]  = nSigmaTOF_v0Pos;
      fVar_Lc[4]  = nSigmaTOF_v0Neg;
      fVar_Lc[5]  = nSigmaTOF_bach;
    }
    else {
      fVar_Lc[0]  = nSigmaTPC_v0Neg; // Anti-Proton
      fVar_Lc[1]  = nSigmaTPC_v0Pos; // Pion+
      fVar_Lc[2]  = nSigmaTPC_bach;
      fVar_Lc[3]  = nSigmaTOF_v0Neg; // Anti-Proton
      fVar_Lc[4]  = nSigmaTOF_v0Pos; // Pion+
      fVar_Lc[5]  = nSigmaTOF_bach;
    }
    fVar_Lc[6]  = v0->GetDCA(); // DCA_V0Dau
    fVar_Lc[7]  = kfpV0.GetChi2()/kfpV0.GetNDF(); //chi2_geometry of V0 (without mass constraint)
    fVar_Lc[8]  = kfpV0_PV.GetChi2()/kfpV0_PV.GetNDF(); //chi2_topological of V0 (with mass constraint) to PV
    fVar_Lc[9]  = AliVertexingHFUtils::ldlFromKF(kfpV0, PV); // ldl_V0
    fVar_Lc[10] = kfpLc.GetChi2()/kfpLc.GetNDF(); // chi2geo_Lc
    fVar_Lc[11] = kfpLc_PV.GetChi2()/kfpLc_PV.GetNDF(); // chi2topo_Lc
    fVar_Lc[12] = AliVertexingHFUtils::ldlFromKF(kfpLc, PV); // ldl_Lc
    fVar_Lc[13] = DecayLxy_V0;
    fVar_Lc[14] = ct_V0;
    fVar_Lc[15] = DecayLxy_Lc;
    fVar_Lc[16] = TMath::ACos(cosPA_V0); // PA_V0
    fVar_Lc[17] = TMath::ACos(cosPA_Lc);  // PA_Lc
    fVar_Lc[18] = kfpV0.GetPt();
    fVar_Lc[19] = mass_V0_rec;
    fVar_Lc[20] = pT_Lc; //pt of Lc (with mass const. of Ks0 and with PV const.)
    fVar_Lc[21] = kfpLc_PV.GetRapidity();
    fVar_Lc[22] = massLc_rec; //mass of Lc (with mass const. of Ks0 and with PV const.)
    fVar_Lc[23] = kfpBach.GetPt();
    if (trackBach->Charge()>0){
      fVar_Lc[24] = v0Pos->Pt(); // proton
      fVar_Lc[25] = v0Neg->Pt(); // pion-
    }
    else {
      fVar_Lc[24] = v0Neg->Pt(); // anti-proton
      fVar_Lc[25] = v0Pos->Pt(); // pion+
    }
    fVar_Lc[26] = Lc2pKs0orLpi->Getd0Prong(0); //rphi impact params of bachlor w.r.t. Primary Vtx [cm]
    fVar_Lc[27] = Lc2pKs0orLpi->Getd0Prong(1); ////rphi impact params of V0 w.r.t. Primary Vtx [cm]
    fVar_Lc[28] = lab_Lc;
    fVar_Lc[29] = AliVertexingHFUtils::CosThetaStarFromKF(0, 4122, 211, 3122, kfpLc, kfpBach_Lc, kfpV0_Lc);  ///cos theta-star
    fVar_Lc[30] = probProton;
    fVar_Lc[31] = v0->PtArmV0() / TMath::Abs(v0->AlphaV0()); //armenteros qT/|alpha|
    if( trackBach->Charge()>0) {
      fVar_Lc[32] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nSigmaTPC_v0Pos, nSigmaTOF_v0Pos);
      fVar_Lc[33] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nSigmaTPC_v0Pos_excl, nSigmaTOF_v0Pos_excl);
    }
    else {
      fVar_Lc[32] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nSigmaTPC_v0Neg, nSigmaTOF_v0Neg);
      fVar_Lc[33] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nSigmaTPC_v0Neg_excl, nSigmaTOF_v0Neg_excl);
    }
    fVar_Lc[34] = cosPA_V0;
    fVar_Lc[35] = AliVertexingHFUtils::DecayLengthFromKF(kfpV0,PV) ;   //d_len_K0s;
    if (fIsMC && fUseWeights && lab_Lc >= 0){
      fVar_Lc[36] = fFuncWeightPythia->Eval(mcLc->Pt()); // weight pT flat
      fVar_Lc[37] = fFuncWeightFONLL5overLHC13d3->Eval(mcLc->Pt()); // weight pT flat
      fVar_Lc[38] = fFuncWeightFONLL5overLHC13d3Lc->Eval(mcLc->Pt()); // weight pT flat
    }
    if (fUseMult) {
      fVar_Lc[39] = nTrackletsEta10;
      fVar_Lc[40] = nTrackletsEta10Corr;
    }
  }


  // === QA tree ===
  fVar_Lc_QA[0]  = kfpV0.GetRapidity(); //rapidity of v0 (without mass const.)
  fVar_Lc_QA[1]  = kfpV0_massConstraint.GetPt(); //pt of V0 (with mass const.)
  fVar_Lc_QA[2]  = kfpV0_massConstraint.GetRapidity(); //rapidity of V0 (with mass const.)
  Float_t mass_V0_massConst_rec=0., err_mass_V0_massConst_rec=0.;
  kfpV0_massConstraint.GetMass(mass_V0_massConst_rec, err_mass_V0_massConst_rec);
  fVar_Lc_QA[3]  = mass_V0_massConst_rec; //mass of V0 (with mass const.)
  fVar_Lc_QA[4]  = kfpV0_massConstraint.GetChi2()/kfpV0_massConstraint.GetNDF(); //chi2_MassConst. of V0

  KFParticle kfpV0_woMassConst_PV = kfpV0;
  kfpV0_woMassConst_PV.SetProductionVertex(PV);
  fVar_Lc_QA[5]  = kfpV0_woMassConst_PV.GetChi2()/kfpV0_woMassConst_PV.GetNDF();

  fVar_Lc_QA[6]  = kfpLc_woV0MassConst.GetChi2()/kfpLc_woV0MassConst.GetNDF();
  fVar_Lc_QA[7]  = kfpLc_woV0MassConst.GetPt();
  fVar_Lc_QA[8]  = kfpLc_woV0MassConst.GetRapidity();
  Float_t mass_Lc_woV0MassConst_rec=0., err_mass_Lc_woV0MassConst_rec=0.;
  kfpLc_woV0MassConst.GetMass(mass_Lc_woV0MassConst_rec, err_mass_Lc_woV0MassConst_rec);
  fVar_Lc_QA[9]  = mass_Lc_woV0MassConst_rec;

  fVar_Lc_QA[10] = kfpLc.GetPt();
  fVar_Lc_QA[11] = kfpLc.GetRapidity();
  Float_t mass_Lc_woPV_rec=0., err_mass_Lc_woPV_rec=0.;
  kfpLc.GetMass(mass_Lc_woPV_rec, err_mass_Lc_woPV_rec);
  fVar_Lc_QA[12] = mass_Lc_woPV_rec;

  // --- reconstruction logic 2 (set production vertex of V0 to PV) ---
  //
  // reconstruct Lc (with mass and PV const. of V0)
  fVar_Lc_QA[13] = kfpV0_PV.GetPt();
  fVar_Lc_QA[14] = kfpV0_PV.GetRapidity();
  Float_t mass_V0_wMassPVConst_rec=0., err_mass_V0_wMassPVConst_rec=0.;
  kfpV0_PV.GetMass(mass_V0_wMassPVConst_rec, err_mass_V0_wMassPVConst_rec);
  fVar_Lc_QA[15] = mass_V0_wMassPVConst_rec;

  KFParticle kfpLc_wV0PV;
  const KFParticle *LcDaughters_wV0PV[2] = {&kfpBach, &kfpV0_PV};
  kfpLc_wV0PV.Construct(LcDaughters_wV0PV, 2);
  fVar_Lc_QA[16] = kfpLc_wV0PV.GetChi2()/kfpLc_wV0PV.GetNDF();
  fVar_Lc_QA[17] = kfpLc_wV0PV.GetPt();
  fVar_Lc_QA[18] = kfpLc_wV0PV.GetRapidity();
  Float_t mass_Lc_wV0PV_rec=0., err_mass_Lc_wV0PV_rec=0.;
  kfpLc_wV0PV.GetMass(mass_Lc_wV0PV_rec, err_mass_Lc_wV0PV_rec);
  fVar_Lc_QA[19] = mass_Lc_wV0PV_rec;
  fVar_Lc_QA[20] = AliVertexingHFUtils::ldlFromKF(kfpLc_wV0PV, PV);

  KFParticle kfpLc_PV_wV0PV = kfpLc_wV0PV;
  kfpLc_PV_wV0PV.SetProductionVertex(PV);
  fVar_Lc_QA[21] = kfpLc_PV_wV0PV.GetChi2()/kfpLc_PV_wV0PV.GetNDF();

  fVar_Lc_QA[22] = lab_V0;
  fVar_Lc_QA[23] = lab_Lc;

  fVar_Lc_QA[24] = Lc2pKs0orLpi->InvMassLctoLambdaPi();
  fVar_Lc_QA[25] = Lc2pKs0orLpi->Pt();

  /*
  //Method to get tracklet multiplicity from event
  Int_t countTreta1corr = 0;
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TProfile *estimatorAvg = GetEstimatorHistogram(aodEvent);
  if (estimatorAvg) {
    countTreta1corr = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,fNTracklets,fVtx1->GetZ(),fRefMult));
  }
  fVar_Lc[43] = countTreta1corr; // NtrkCorr
  */

  if (fIsMC && !(lab_Lc < 0 && fKeepOnlyMCSignal)) fTree_Lc->Fill();
  if (!fIsMC) fTree_Lc->Fill();

  if (fWriteLcQATree) {fTree_Lc_QA->Fill();}

  /*
  cout << "==========" << endl;
  cout << "kfpV0: " << kfpV0 << endl;
  cout << "kfpV0_massConstraint: " << kfpV0_massConstraint << endl;
  cout << "kfpV0_PV: " << kfpV0_PV << endl;
  cout << "PV: " << PV << endl;
  cout << "==========" << endl;
  */

  return;

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSELc2pKs0fromKFP::MatchToMCLc2pKs0(AliAODTrack *v0Pos, AliAODTrack *v0Neg, AliAODTrack *bachPart, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelPionPlus = fabs(v0Pos->GetLabel());
  if (labelPionPlus<=0) return -1;
  AliAODMCParticle *mcPionPlus = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus));

  Int_t labelPionMinus = fabs(v0Neg->GetLabel());
  if (labelPionMinus<=0) return -1;
  AliAODMCParticle *mcPionMinus = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus));

  Int_t labelProton = fabs(bachPart->GetLabel());
  if (labelProton<=0) return -1;
  AliAODMCParticle *mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));

  if ( mcPionPlus->GetPdgCode()!=211 || mcPionMinus->GetPdgCode()!=-211 || fabs(mcProton->GetPdgCode())!=2212 ) return -1;

  Int_t IndexMother[2] = {-9999, -9999};

  IndexMother[0] = mcPionPlus->GetMother();
  IndexMother[1] = mcPionMinus->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( fabs(mcMother->GetPdgCode()) != 310 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Ks0 and have two daughters

  IndexMother[0] = mcMother->GetMother(); // mother of Ks0
  if ( IndexMother[0]<0 ) return -1; // check mother exist
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( fabs(mcMother->GetPdgCode()) != 311 || mcMother->GetNDaughters()!=1 ) return -1; // check mother is K0 and have only one daughter

  IndexMother[0] = mcMother->GetMother(); // mother of K0
  IndexMother[1] = mcProton->GetMother(); // mother of proton
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( fabs(mcMother->GetPdgCode()) != 4122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Lc and have two daughters

//  if ( mcMother->IsPrimary() ) return 1;
//  if ( mcMother->IsPhysicalPrimary() ) return 2;
//  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
//  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
//  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

  Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcMother,kTRUE);
  return CheckOrigin;

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSELc2pKs0fromKFP::MatchToMCKs0(AliAODTrack *v0Pos, AliAODTrack *v0Neg, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelPionPlus = fabs(v0Pos->GetLabel());
  if (labelPionPlus<=0) return -1;
  AliAODMCParticle *mcPionPlus = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus));

  Int_t labelPionMinus = fabs(v0Neg->GetLabel());
  if (labelPionMinus<=0) return -1;
  AliAODMCParticle *mcPionMinus = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus));

  if ( mcPionPlus->GetPdgCode()!=211 || mcPionMinus->GetPdgCode()!=-211 ) return -1;

  Int_t IndexMother[2] = {-9999, -9999};

  IndexMother[0] = mcPionPlus->GetMother();
  IndexMother[1] = mcPionMinus->GetMother();

  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( fabs(mcMother->GetPdgCode()) != 310 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Ks0 and have two daughters

  return 1;

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSELc2pKs0fromKFP::MatchToMCLam(AliAODTrack *v0Pos, AliAODTrack *v0Neg, TClonesArray *mcArray, Bool_t IsParticle)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelPos = fabs(v0Pos->GetLabel());
  if (labelPos<=0) return -1;
  AliAODMCParticle *mcPos = static_cast<AliAODMCParticle*>(mcArray->At(labelPos));

  Int_t labelNeg = fabs(v0Neg->GetLabel());
  if (labelNeg<=0) return -1;
  AliAODMCParticle *mcNeg = static_cast<AliAODMCParticle*>(mcArray->At(labelNeg));

  if (IsParticle) {
    if ( mcPos->GetPdgCode()!=2212 || mcNeg->GetPdgCode()!=-211 ) return -1;
  }
  if (!IsParticle) {
    if ( mcPos->GetPdgCode()!=211 || mcNeg->GetPdgCode()!=-2212 ) return -1;
  }

  Int_t IndexMother[2] = {-9999, -9999};

  IndexMother[0] = mcPos->GetMother();
  IndexMother[1] = mcNeg->GetMother();

  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( fabs(mcMother->GetPdgCode()) != 3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Ks0 and have two daughters

  return 1;

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSELc2pKs0fromKFP::MatchToMCLc2Lpi(AliAODTrack *v0Pos, AliAODTrack *v0Neg, AliAODTrack *bachPart, TClonesArray *mcArray, Bool_t IsParticle)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelPos = fabs(v0Pos->GetLabel());
  if (labelPos<=0) return -1;
  AliAODMCParticle *mcPos = static_cast<AliAODMCParticle*>(mcArray->At(labelPos));

  Int_t labelNeg = fabs(v0Neg->GetLabel());
  if (labelNeg<=0) return -1;
  AliAODMCParticle *mcNeg = static_cast<AliAODMCParticle*>(mcArray->At(labelNeg));

  Int_t labelBach = fabs(bachPart->GetLabel());
  if (labelBach<=0) return -1;
  AliAODMCParticle *mcBach = static_cast<AliAODMCParticle*>(mcArray->At(labelBach));

  if (IsParticle) {
    if ( mcPos->GetPdgCode()!=2212 || mcNeg->GetPdgCode()!=-211 || mcBach->GetPdgCode()!=211 ) return -1;
  }
  if (!IsParticle) {
    if ( mcPos->GetPdgCode()!=211 || mcNeg->GetPdgCode()!=-2212 || mcBach->GetPdgCode()!=-211 ) return -1;
  }

  Int_t IndexMother[2] = {-9999, -9999};

  IndexMother[0] = mcPos->GetMother();
  IndexMother[1] = mcNeg->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( fabs(mcMother->GetPdgCode()) != 3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is K0 and have only one daughter

  IndexMother[0] = mcMother->GetMother(); // mother of Lam
  IndexMother[1] = mcBach->GetMother(); // mother of bachlor
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( fabs(mcMother->GetPdgCode()) != 4122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Lc and have two daughters

  Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcMother,kTRUE);
  return CheckOrigin;

}


TProfile* AliAnalysisTaskSELc2pKs0fromKFP::GetEstimatorHistogram(const AliVEvent* event)  {


  Int_t runNo = event->GetRunNumber();
  Int_t period = -1;
  switch (fAnalysisType) {    // flag to set which system and year is being used
      case kpPb2016: //0 = LHC16q, 265499 -- 265525 || 265309 -- 265387, 1 = LHC16q, 265435, 2 = LHC16q, 265388 -- 265427, 3 = LHC16t, 267163 -- 267166
         if ((runNo >=265499 && runNo <=265525) || (runNo >= 265309 && runNo <= 265387)) period = 0;
         else if (runNo == 265435) period = 1;
         else if (runNo >= 265388 && runNo <= 265427) period = 2;
         else if (runNo >=267163 && runNo <=276166) period = 3;
         if (period < 0 || period > 3) { AliInfo(Form("Run number %d not found for LHC16 pPb!",runNo)); return 0;}
         break;
      case kpp2016: //0 = LHC16j, 1 = LHC16k, 2 = LHC16l
         if (runNo >= 256219 && runNo <= 256418) period = 0;
         else if (runNo >= 256504 && runNo <= 258537) period = 1;
         else if (runNo >= 258883 && runNo <= 260187) period = 2;
         if (period < 0 || period > 2) {AliInfo(Form("Run number %d not found for LHC16 pp!",runNo)); return 0;}
         break;
      default:       //no valid switch
         return 0;
      break;
  }

  return fMultEstimatorAvg[period];

}