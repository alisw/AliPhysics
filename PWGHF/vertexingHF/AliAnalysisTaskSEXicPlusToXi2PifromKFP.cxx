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
// Author: Jianhui Zhu (1,2)
// (1) Central China Normal University
// (2) GSI Helmholtz Centre for Heavy Ion Research
// E-mail: zjh@ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <TObjString.h>
#include <iostream>
#include <iomanip>
#include <TDatabasePDG.h>
#include <vector>
#include <TVector3.h>
#include <TFile.h>
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLine.h"
#include "TList.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliNormalizationCounter.h"
#include "AliAnalysisTaskSEXicPlusToXi2PifromKFP.h"
#include "AliPIDResponse.h"

#include "AliAODMCParticle.h"

// includes added to play with KFParticle
#ifndef HomogeneousField
#define HomogeneousField 
#endif

using std::cout;
using std::endl;

class AliAnalysisTaskSEXicPlusToXi2PifromKFP;    // your analysis class

ClassImp(AliAnalysisTaskSEXicPlusToXi2PifromKFP) // classimp: necessary for root

AliAnalysisTaskSEXicPlusToXi2PifromKFP::AliAnalysisTaskSEXicPlusToXi2PifromKFP() :
  AliAnalysisTaskSE(),
  fIsMC(kFALSE),
  fPID(0),
  fAnaCuts(0),
  fpVtx(0),
  fMCEvent(0),
  fBzkG(0),
  fCentrality(0),
  fAodTrackInd(0),
  fOutputList(0),
  fListCuts(0),
  fTree_Event(0),
  fVar_Event(0),
  fTree_XicPlus(0),
  fVar_XicPlus(0),
  fTree_XicPlus_QA(0),
  fVar_XicPlus_QA(0),
  fTree_XicPlus_QA_woMassConstForLamAndXi(0),
  fVar_XicPlus_QA_woMassConstForLamAndXi(0),
  fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi(0),
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi(0),
  fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi(0),
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi(0),
  fTree_XicPlus_QA_wMassConstForLamAndXi(0),
  fVar_XicPlus_QA_wMassConstForLamAndXi(0),
  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi(0),
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi(0),
  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi(0),
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi(0),
  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic(0),
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic(0),
  fVar_XicPlus_EvtID(0),
  fTree_XicPlusMCGen(0),
  fVar_XicPlusMCGen(0),
  fCounter(0),
  fHistEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fHCountUsedForPrimVtxFit(0),
  fHNumberOfCasc(0),
  fHPrimVtx_woDau_x(0),
  fHPrimVtx_woDau_y(0),
  fHPrimVtx_woDau_z(0),
  fHPrimVtx_woDau_err_x(0),
  fHPrimVtx_woDau_err_y(0),
  fHPrimVtx_woDau_err_z(0),
  fHNumOfCandidatePerEvent_In3sigma(0),
  fCount_NumOfCandidatePerEvent_In3Sigma(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_x(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_y(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_z(0),
  fHPrimVtx_recalPV_Minus_PVrec_x(0),
  fHPrimVtx_recalPV_Minus_PVrec_y(0),
  fHPrimVtx_recalPV_Minus_PVrec_z(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_x(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_y(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_z(0),
  fHPrimVtx_recalPV_KF_Minus_PV_x(0),
  fHPrimVtx_recalPV_KF_Minus_PV_y(0),
  fHPrimVtx_recalPV_KF_Minus_PV_z(0),
  fHPrimVtx_PVrec_Minus_PVgen_x(0),
  fHPrimVtx_PVrec_Minus_PVgen_y(0),
  fHPrimVtx_PVrec_Minus_PVgen_z(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_x(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_y(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_z(0),
  fHPrimVtx_recalPV_Minus_PVgen_x(0),
  fHPrimVtx_recalPV_Minus_PVgen_y(0),
  fHPrimVtx_recalPV_Minus_PVgen_z(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_x(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_y(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_z(0),
  fHPrimVtx_PV_PULL_x(0),
  fHPrimVtx_PV_PULL_y(0),
  fHPrimVtx_PV_PULL_z(0),
  fHPrimVtx_PV_KF_PULL_x(0),
  fHPrimVtx_PV_KF_PULL_y(0),
  fHPrimVtx_PV_KF_PULL_z(0),
  fHPrimVtx_PV_KF_Refit_PULL_x(0),
  fHPrimVtx_PV_KF_Refit_PULL_y(0),
  fHPrimVtx_PV_KF_Refit_PULL_z(0),
  fHPrimVtx_recalPV_PULL_x(0),
  fHPrimVtx_recalPV_PULL_y(0),
  fHPrimVtx_recalPV_PULL_z(0),
  fHPrimVtx_recalPV_KF_Refit_PULL_x(0),
  fHPrimVtx_recalPV_KF_Refit_PULL_y(0),
  fHPrimVtx_recalPV_KF_Refit_PULL_z(0),
  fHPrimVtx_recalPV_KF_Minus_PVgen_x(0),
  fHPrimVtx_recalPV_KF_Minus_PVgen_y(0),
  fHPrimVtx_recalPV_KF_Minus_PVgen_z(0),
  fHPrimVtx_recalPV_KF_PULL_x(0),
  fHPrimVtx_recalPV_KF_PULL_y(0),
  fHPrimVtx_recalPV_KF_PULL_z(0),
  fFileName(""),
  fEventNumber(0),
  fDirNumber(0),
  fWriteXicPlusTree(kFALSE),
  fWriteXicPlusQATree(kFALSE),
  fWriteXicPlusMCGenTree(kFALSE)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskSEXicPlusToXi2PifromKFP::AliAnalysisTaskSEXicPlusToXi2PifromKFP(const char* name, AliRDHFCutsKFP* cuts) :
  AliAnalysisTaskSE(name),
  fIsMC(kFALSE),
  fPID(0),
  fAnaCuts(cuts),
  fpVtx(0),
  fMCEvent(0),
  fBzkG(0),
  fCentrality(0),
  fAodTrackInd(0),
  fOutputList(0),
  fListCuts(0),
  fTree_Event(0),
  fVar_Event(0),
  fTree_XicPlus(0),
  fVar_XicPlus(0),
  fTree_XicPlus_QA(0),
  fVar_XicPlus_QA(0),
  fTree_XicPlus_QA_woMassConstForLamAndXi(0),
  fVar_XicPlus_QA_woMassConstForLamAndXi(0),
  fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi(0),
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi(0),
  fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi(0),
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi(0),
  fTree_XicPlus_QA_wMassConstForLamAndXi(0),
  fVar_XicPlus_QA_wMassConstForLamAndXi(0),
  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi(0),
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi(0),
  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi(0),
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi(0),
  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic(0),
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic(0),
  fVar_XicPlus_EvtID(0),
  fTree_XicPlusMCGen(0),
  fVar_XicPlusMCGen(0),
  fCounter(0),
  fHistEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fHCountUsedForPrimVtxFit(0),
  fHNumberOfCasc(0),
  fHPrimVtx_woDau_x(0),
  fHPrimVtx_woDau_y(0),
  fHPrimVtx_woDau_z(0),
  fHPrimVtx_woDau_err_x(0),
  fHPrimVtx_woDau_err_y(0),
  fHPrimVtx_woDau_err_z(0),
  fHNumOfCandidatePerEvent_In3sigma(0),
  fCount_NumOfCandidatePerEvent_In3Sigma(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_x(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_y(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_z(0),
  fHPrimVtx_recalPV_Minus_PVrec_x(0),
  fHPrimVtx_recalPV_Minus_PVrec_y(0),
  fHPrimVtx_recalPV_Minus_PVrec_z(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_x(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_y(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_z(0),
  fHPrimVtx_recalPV_KF_Minus_PV_x(0),
  fHPrimVtx_recalPV_KF_Minus_PV_y(0),
  fHPrimVtx_recalPV_KF_Minus_PV_z(0),
  fHPrimVtx_PVrec_Minus_PVgen_x(0),
  fHPrimVtx_PVrec_Minus_PVgen_y(0),
  fHPrimVtx_PVrec_Minus_PVgen_z(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_x(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_y(0),
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_z(0),
  fHPrimVtx_recalPV_Minus_PVgen_x(0),
  fHPrimVtx_recalPV_Minus_PVgen_y(0),
  fHPrimVtx_recalPV_Minus_PVgen_z(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_x(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_y(0),
  fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_z(0),
  fHPrimVtx_PV_PULL_x(0),
  fHPrimVtx_PV_PULL_y(0),
  fHPrimVtx_PV_PULL_z(0),
  fHPrimVtx_PV_KF_PULL_x(0),
  fHPrimVtx_PV_KF_PULL_y(0),
  fHPrimVtx_PV_KF_PULL_z(0),
  fHPrimVtx_PV_KF_Refit_PULL_x(0),
  fHPrimVtx_PV_KF_Refit_PULL_y(0),
  fHPrimVtx_PV_KF_Refit_PULL_z(0),
  fHPrimVtx_recalPV_PULL_x(0),
  fHPrimVtx_recalPV_PULL_y(0),
  fHPrimVtx_recalPV_PULL_z(0),
  fHPrimVtx_recalPV_KF_Refit_PULL_x(0),
  fHPrimVtx_recalPV_KF_Refit_PULL_y(0),
  fHPrimVtx_recalPV_KF_Refit_PULL_z(0),
  fHPrimVtx_recalPV_KF_Minus_PVgen_x(0),
  fHPrimVtx_recalPV_KF_Minus_PVgen_y(0),
  fHPrimVtx_recalPV_KF_Minus_PVgen_z(0),
  fHPrimVtx_recalPV_KF_PULL_x(0),
  fHPrimVtx_recalPV_KF_PULL_y(0),
  fHPrimVtx_recalPV_KF_PULL_z(0),
  fFileName(""),
  fEventNumber(0),
  fDirNumber(0),
  fWriteXicPlusTree(kFALSE),
  fWriteXicPlusQATree(kFALSE),
  fWriteXicPlusMCGenTree(kFALSE)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
  DefineOutput(2, AliNormalizationCounter::Class());
  DefineOutput(3, TTree::Class()); // event
  DefineOutput(4, TTree::Class()); // XicPlus
  DefineOutput(5, TTree::Class()); // XicPlus MCGen
  DefineOutput(6, TList::Class()); // XicPlus event & trigger
  DefineOutput(7, TTree::Class()); // XicPlus QA
  DefineOutput(8, TTree::Class()); // XicPlus QA woMassConstForLamAndXi
  DefineOutput(9, TTree::Class()); // XicPlus QA wMassConstForLam_woMassConstForXi
  DefineOutput(10, TTree::Class()); // XicPlus QA woMassConstForLam_wMassConstForXi
  DefineOutput(11, TTree::Class()); // XicPlus QA wMassConstForLamAndXi
  DefineOutput(12, TTree::Class()); // XicPlus QA wMassAndTopoConstForLam_wMassConstForXi
  DefineOutput(13, TTree::Class()); // XicPlus QA wMassAndTopoConstForLam_wMassAndTopoConstForXi
  DefineOutput(14, TTree::Class()); // XicPlus QA wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic

}
//_____________________________________________________________________________
AliAnalysisTaskSEXicPlusToXi2PifromKFP::~AliAnalysisTaskSEXicPlusToXi2PifromKFP()
{
    // destructor
    if (fOutputList) {
      delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
      fOutputList = 0;
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

    if (fTree_XicPlus) {
      delete fTree_XicPlus;
      fTree_XicPlus = 0;
    }

    if (fVar_XicPlus) {
      delete fVar_XicPlus;
      fVar_XicPlus = 0;
    }

    if (fTree_XicPlus_QA) {
      delete fTree_XicPlus_QA;
      fTree_XicPlus_QA = 0;
    }

    if (fVar_XicPlus_QA) {
      delete fVar_XicPlus_QA;
      fVar_XicPlus_QA = 0;
    }

    if (fTree_XicPlus_QA_woMassConstForLamAndXi) {
      delete fTree_XicPlus_QA_woMassConstForLamAndXi;
      fTree_XicPlus_QA_woMassConstForLamAndXi = 0;
    }

    if (fVar_XicPlus_QA_woMassConstForLamAndXi) {
      delete fVar_XicPlus_QA_woMassConstForLamAndXi;
      fVar_XicPlus_QA_woMassConstForLamAndXi = 0;
    }

    if (fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi) {
      delete fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi;
      fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi = 0;
    }

    if (fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi) {
      delete fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi;
      fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi = 0;
    }

    if (fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi) {
      delete fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi;
      fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi = 0;
    }

    if (fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi) {
      delete fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi;
      fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi = 0;
    }

    if (fTree_XicPlus_QA_wMassConstForLamAndXi) {
      delete fTree_XicPlus_QA_wMassConstForLamAndXi;
      fTree_XicPlus_QA_wMassConstForLamAndXi = 0;
    }

    if (fVar_XicPlus_QA_wMassConstForLamAndXi) {
      delete fVar_XicPlus_QA_wMassConstForLamAndXi;
      fVar_XicPlus_QA_wMassConstForLamAndXi = 0;
    }

    if (fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi) {
      delete fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi;
      fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi = 0;
    }

    if (fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi) {
      delete fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi;
      fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi = 0;
    }

    if (fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi) {
      delete fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi;
      fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi = 0;
    }

    if (fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi) {
      delete fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi;
      fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi = 0;
    }

    if (fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic) {
      delete fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic;
      fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic = 0;
    }

    if (fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic) {
      delete fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic;
      fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic = 0;
    }

    if (fTree_XicPlusMCGen) {
      delete fTree_XicPlusMCGen;
      fTree_XicPlusMCGen = 0;
    }

    if (fVar_XicPlusMCGen) {
      delete fVar_XicPlusMCGen;
      fVar_XicPlusMCGen = 0;
    }

    if (fCounter) {
      delete fCounter;
      fCounter = 0;
    }


}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::Init()
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
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::UserCreateOutputObjects()
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

  fHCentrality = new TH1F("fHCentrality", "counter", 100, 0., 100.);

  fHCountUsedForPrimVtxFit = new TH1F("fHCountUsedForPrimVtxFit", "counter", 6, -0.5, 5.5);
  fHCountUsedForPrimVtxFit->GetXaxis()->SetBinLabel(1,"#pi_{0}^{+}");
  fHCountUsedForPrimVtxFit->GetXaxis()->SetBinLabel(2,"#pi_{1}^{+}");
  fHCountUsedForPrimVtxFit->GetXaxis()->SetBinLabel(3,"#pi_{2}^{#minus}");
  fHCountUsedForPrimVtxFit->GetXaxis()->SetBinLabel(4,"#pi_{3}^{#minus}");
  fHCountUsedForPrimVtxFit->GetXaxis()->SetBinLabel(5,"p");
  fHCountUsedForPrimVtxFit->GetXaxis()->SetBinLabel(6,"#Xi_{c}^{+}");

  fHNumberOfCasc = new TH1F("fHNumberOfCasc", "counter", 21, -0.5, 20.5);

  fHPrimVtx_woDau_x = new TH1F("fHPrimVtx_woDau_x", "PV_woDau_x", 20000, -10, 10);
  fHPrimVtx_woDau_y = new TH1F("fHPrimVtx_woDau_y", "PV_woDau_y", 20000, -10, 10);
  fHPrimVtx_woDau_z = new TH1F("fHPrimVtx_woDau_z", "PV_woDau_z", 20000, -10, 10);
  fHPrimVtx_woDau_err_x = new TH1F("fHPrimVtx_woDau_err_x", "PV_woDau_err_x", 10000, 0, 10);
  fHPrimVtx_woDau_err_y = new TH1F("fHPrimVtx_woDau_err_y", "PV_woDau_err_y", 10000, 0, 10);
  fHPrimVtx_woDau_err_z = new TH1F("fHPrimVtx_woDau_err_z", "PV_woDau_err_z", 10000, 0, 10);
  fHNumOfCandidatePerEvent_In3sigma = new TH1F("fHNumOfCandidatePerEvent_In3sigma", "NumOfCandidatePerEvent_In3sigma", 101, -0.5, 100.5);
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_x = new TH1F("fHPrimVtx_PV_KF_Refit_Minus_PVrec_x", "(PV_KF_Refit - PVrec) x", 2000, -0.1, 0.1);
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_y = new TH1F("fHPrimVtx_PV_KF_Refit_Minus_PVrec_y", "(PV_KF_Refit - PVrec) y", 2000, -0.1, 0.1);
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_z = new TH1F("fHPrimVtx_PV_KF_Refit_Minus_PVrec_z", "(PV_KF_Refit - PVrec) z", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_Minus_PVrec_x = new TH1F("fHPrimVtx_recalPV_Minus_PVrec_x", "(recalPV - PVrec) x", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_Minus_PVrec_y = new TH1F("fHPrimVtx_recalPV_Minus_PVrec_y", "(recalPV - PVrec) y", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_Minus_PVrec_z = new TH1F("fHPrimVtx_recalPV_Minus_PVrec_z", "(recalPV - PVrec) z", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_x = new TH1F("fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_x", "(recalPV_KF_Refit - PV_KF_Refit) x", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_y = new TH1F("fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_y", "(recalPV_KF_Refit - PV_KF_Refit) y", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_z = new TH1F("fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_z", "(recalPV_KF_Refit - PV_KF_Refit) z", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Minus_PV_x = new TH1F("fHPrimVtx_recalPV_KF_Minus_PV_x", "(recalPV_KF - PV) x", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Minus_PV_y = new TH1F("fHPrimVtx_recalPV_KF_Minus_PV_y", "(recalPV_KF - PV) y", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Minus_PV_z = new TH1F("fHPrimVtx_recalPV_KF_Minus_PV_z", "(recalPV_KF - PV) z", 2000, -0.1, 0.1);
  fHPrimVtx_PVrec_Minus_PVgen_x = new TH1F("fHPrimVtx_PVrec_Minus_PVgen_x", "(PVrec - PVgen) x", 2000, -0.1, 0.1);
  fHPrimVtx_PVrec_Minus_PVgen_y = new TH1F("fHPrimVtx_PVrec_Minus_PVgen_y", "(PVrec - PVgen) y", 2000, -0.1, 0.1);
  fHPrimVtx_PVrec_Minus_PVgen_z = new TH1F("fHPrimVtx_PVrec_Minus_PVgen_z", "(PVrec - PVgen) z", 2000, -0.1, 0.1);
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_x = new TH1F("fHPrimVtx_PV_KF_Refit_Minus_PVgen_x", "(PV_KF_Refit - PVgen) x", 2000, -0.1, 0.1);
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_y = new TH1F("fHPrimVtx_PV_KF_Refit_Minus_PVgen_y", "(PV_KF_Refit - PVgen) y", 2000, -0.1, 0.1);
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_z = new TH1F("fHPrimVtx_PV_KF_Refit_Minus_PVgen_z", "(PV_KF_Refit - PVgen) z", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_Minus_PVgen_x = new TH1F("fHPrimVtx_recalPV_Minus_PVgen_x", "(recalPV - PVgen) x", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_Minus_PVgen_y = new TH1F("fHPrimVtx_recalPV_Minus_PVgen_y", "(recalPV - PVgen) y", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_Minus_PVgen_z = new TH1F("fHPrimVtx_recalPV_Minus_PVgen_z", "(recalPV - PVgen) z", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_x = new TH1F("fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_x", "(recalPV_KF_Refit - PVgen) x", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_y = new TH1F("fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_y", "(recalPV_KF_Refit - PVgen) y", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_z = new TH1F("fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_z", "(recalPV_KF_Refit - PVgen) z", 2000, -0.1, 0.1);
  fHPrimVtx_PV_PULL_x = new TH1F("fHPrimVtx_PV_PULL_x", "PV_PULL_x", 200, -10, 10);
  fHPrimVtx_PV_PULL_y = new TH1F("fHPrimVtx_PV_PULL_y", "PV_PULL_y", 200, -10, 10);
  fHPrimVtx_PV_PULL_z = new TH1F("fHPrimVtx_PV_PULL_z", "PV_PULL_z", 200, -10, 10);
  fHPrimVtx_PV_KF_PULL_x = new TH1F("fHPrimVtx_PV_KF_PULL_x", "PV_KF_PULL_x", 200, -10, 10);
  fHPrimVtx_PV_KF_PULL_y = new TH1F("fHPrimVtx_PV_KF_PULL_y", "PV_KF_PULL_y", 200, -10, 10);
  fHPrimVtx_PV_KF_PULL_z = new TH1F("fHPrimVtx_PV_KF_PULL_z", "PV_KF_PULL_z", 200, -10, 10);
  fHPrimVtx_PV_KF_Refit_PULL_x = new TH1F("fHPrimVtx_PV_KF_Refit_PULL_x", "PV_KF_Refit_PULL_x", 200, -10, 10);
  fHPrimVtx_PV_KF_Refit_PULL_y = new TH1F("fHPrimVtx_PV_KF_Refit_PULL_y", "PV_KF_Refit_PULL_y", 200, -10, 10);
  fHPrimVtx_PV_KF_Refit_PULL_z = new TH1F("fHPrimVtx_PV_KF_Refit_PULL_z", "PV_KF_Refit_PULL_z", 200, -10, 10);
  fHPrimVtx_recalPV_PULL_x = new TH1F("fHPrimVtx_recalPV_PULL_x", "recalPV_PULL_x", 200, -10, 10);
  fHPrimVtx_recalPV_PULL_y = new TH1F("fHPrimVtx_recalPV_PULL_y", "recalPV_PULL_y", 200, -10, 10);
  fHPrimVtx_recalPV_PULL_z = new TH1F("fHPrimVtx_recalPV_PULL_z", "recalPV_PULL_z", 200, -10, 10);
  fHPrimVtx_recalPV_KF_Refit_PULL_x = new TH1F("fHPrimVtx_recalPV_KF_Refit_PULL_x", "recalPV_KF_Refit_PULL_x", 200, -10, 10);
  fHPrimVtx_recalPV_KF_Refit_PULL_y = new TH1F("fHPrimVtx_recalPV_KF_Refit_PULL_y", "recalPV_KF_Refit_PULL_y", 200, -10, 10);
  fHPrimVtx_recalPV_KF_Refit_PULL_z = new TH1F("fHPrimVtx_recalPV_KF_Refit_PULL_z", "recalPV_KF_Refit_PULL_z", 200, -10, 10);
  fHPrimVtx_recalPV_KF_Minus_PVgen_x = new TH1F("fHPrimVtx_recalPV_KF_Minus_PVgen_x", "(recalPV_KF - PVgen) x", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Minus_PVgen_y = new TH1F("fHPrimVtx_recalPV_KF_Minus_PVgen_y", "(recalPV_KF - PVgen) y", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_Minus_PVgen_z = new TH1F("fHPrimVtx_recalPV_KF_Minus_PVgen_z", "(recalPV_KF - PVgen) z", 2000, -0.1, 0.1);
  fHPrimVtx_recalPV_KF_PULL_x = new TH1F("fHPrimVtx_recalPV_KF_PULL_x", "recalPV_KF_PULL_x", 200, -10, 10);
  fHPrimVtx_recalPV_KF_PULL_y = new TH1F("fHPrimVtx_recalPV_KF_PULL_y", "recalPV_KF_PULL_y", 200, -10, 10);
  fHPrimVtx_recalPV_KF_PULL_z = new TH1F("fHPrimVtx_recalPV_KF_PULL_z", "recalPV_KF_PULL_z", 200, -10, 10);

  fOutputList->Add(fHistEvents); // don't forget to add it to the list! the list will be written to file, so if you want
  fOutputList->Add(fHTrigger);
  fOutputList->Add(fHCountUsedForPrimVtxFit);
  fOutputList->Add(fHNumberOfCasc);
  fOutputList->Add(fHPrimVtx_woDau_x);
  fOutputList->Add(fHPrimVtx_woDau_y);
  fOutputList->Add(fHPrimVtx_woDau_z);
  fOutputList->Add(fHPrimVtx_woDau_err_x);
  fOutputList->Add(fHPrimVtx_woDau_err_y);
  fOutputList->Add(fHPrimVtx_woDau_err_z);
  fOutputList->Add(fHNumOfCandidatePerEvent_In3sigma);
  fOutputList->Add(fHPrimVtx_PV_KF_Refit_Minus_PVrec_x);
  fOutputList->Add(fHPrimVtx_PV_KF_Refit_Minus_PVrec_y);
  fOutputList->Add(fHPrimVtx_PV_KF_Refit_Minus_PVrec_z);
  fOutputList->Add(fHPrimVtx_recalPV_Minus_PVrec_x);
  fOutputList->Add(fHPrimVtx_recalPV_Minus_PVrec_y);
  fOutputList->Add(fHPrimVtx_recalPV_Minus_PVrec_z);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_x);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_y);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_z);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Minus_PV_x);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Minus_PV_y);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Minus_PV_z);
  fOutputList->Add(fHPrimVtx_PVrec_Minus_PVgen_x);
  fOutputList->Add(fHPrimVtx_PVrec_Minus_PVgen_y);
  fOutputList->Add(fHPrimVtx_PVrec_Minus_PVgen_z);
  fOutputList->Add(fHPrimVtx_PV_KF_Refit_Minus_PVgen_x);
  fOutputList->Add(fHPrimVtx_PV_KF_Refit_Minus_PVgen_y);
  fOutputList->Add(fHPrimVtx_PV_KF_Refit_Minus_PVgen_z);
  fOutputList->Add(fHPrimVtx_recalPV_Minus_PVrec_x);
  fOutputList->Add(fHPrimVtx_recalPV_Minus_PVrec_y);
  fOutputList->Add(fHPrimVtx_recalPV_Minus_PVrec_z);
  fOutputList->Add(fHPrimVtx_recalPV_Minus_PVgen_x);
  fOutputList->Add(fHPrimVtx_recalPV_Minus_PVgen_y);
  fOutputList->Add(fHPrimVtx_recalPV_Minus_PVgen_z);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_x);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_y);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_z);
  fOutputList->Add(fHPrimVtx_PV_PULL_x);
  fOutputList->Add(fHPrimVtx_PV_PULL_y);
  fOutputList->Add(fHPrimVtx_PV_PULL_z);
  fOutputList->Add(fHPrimVtx_PV_KF_PULL_x);
  fOutputList->Add(fHPrimVtx_PV_KF_PULL_y);
  fOutputList->Add(fHPrimVtx_PV_KF_PULL_z);
  fOutputList->Add(fHPrimVtx_PV_KF_Refit_PULL_x);
  fOutputList->Add(fHPrimVtx_PV_KF_Refit_PULL_y);
  fOutputList->Add(fHPrimVtx_PV_KF_Refit_PULL_z);
  fOutputList->Add(fHPrimVtx_recalPV_PULL_x);
  fOutputList->Add(fHPrimVtx_recalPV_PULL_y);
  fOutputList->Add(fHPrimVtx_recalPV_PULL_z);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Refit_PULL_x);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Refit_PULL_y);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Refit_PULL_z);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Minus_PVgen_x);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Minus_PVgen_y);
  fOutputList->Add(fHPrimVtx_recalPV_KF_Minus_PVgen_z);
  fOutputList->Add(fHPrimVtx_recalPV_KF_PULL_x);
  fOutputList->Add(fHPrimVtx_recalPV_KF_PULL_y);
  fOutputList->Add(fHPrimVtx_recalPV_KF_PULL_z);

  // Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(2)->GetContainer();
  if(cont) normName = (TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();
  PostData(2, fCounter);
  DefineEvent();
  PostData(3, fTree_Event);  // postdata will notify the analysis manager of changes / updates to the 

  DefineTreeRecXicPlus();
  PostData(4, fTree_XicPlus);

  DefineTreeGenXicPlus();
  PostData(5, fTree_XicPlusMCGen);

  PostData(6, fOutputList);

  DefineTreeQAXicPlus();
  PostData(7, fTree_XicPlus_QA);

  DefineTreeQAXicPlus_woMassConstForLamAndXi();
  PostData(8, fTree_XicPlus_QA_woMassConstForLamAndXi);

  DefineTreeQAXicPlus_wMassConstForLam_woMassConstForXi();
  PostData(9, fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi);

  DefineTreeQAXicPlus_woMassConstForLam_wMassConstForXi();
  PostData(10, fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi);

  DefineTreeQAXicPlus_wMassConstForLamAndXi();
  PostData(11, fTree_XicPlus_QA_wMassConstForLamAndXi);

  DefineTreeQAXicPlus_wMassAndTopoConstForLam_wMassConstForXi();
  PostData(12, fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi);

  DefineTreeQAXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi();
  PostData(13, fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi);

  DefineTreeQAXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic();
  PostData(14, fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic);

  return;
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::UserExec(Option_t *)
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
  AliAODEvent* AODEvent = dynamic_cast<AliAODEvent*>(fInputEvent);    // get an event (called AODEvent) from the input file
                                                        // there's another event format (ESD) which works in a similar way
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with aod's

  fHistEvents->Fill(1);

  //--------------------------------------------------------------
  // First check if the event has magnetic field and proper vertex
  //--------------------------------------------------------------

  fBzkG = (Double_t)AODEvent->GetMagneticField();
  if (TMath::Abs(fBzkG)<0.001) return;
  KFParticle::SetField(fBzkG);

  fpVtx = (AliAODVertex*)AODEvent->GetPrimaryVertex();
  if (!fpVtx) return;
  fHistEvents->Fill(2);

  //cout << "Title: " << fpVtx->GetTitle() << endl;

  fCounter->StoreEvent(AODEvent,fAnaCuts,fIsMC);

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
    mcArray = dynamic_cast<TClonesArray*>(AODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if ( !mcArray ) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    fHistEvents->Fill(7); // number of MC array exist

    // load MC header
    mcHeader = (AliAODMCHeader*)AODEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if ( !mcHeader ) {
      AliError("AliAnalysisTaskSEXicPlusToXi2PifromKFP::UserExec: MC header branch not found!\n");
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
      Bool_t selevt = MakeMCAnalysis(mcArray);
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

  Bool_t IsEventSelected = fAnaCuts->IsEventSelected(AODEvent);
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

//  AliCentrality *cent = AODEvent->GetCentrality();
//  Float_t Centrality = cent->GetCentralityPercentile("V0M");
//  fHCentrality->Fill(Centrality);

  //------------------------------------------------
  // Check if the event has v0 candidate
  //------------------------------------------------
  Int_t num_v0 = AODEvent->GetNumberOfV0s();
  if (num_v0>0) fHistEvents->Fill(5);

  //------------------------------------------------
  // Check if the event has cascade candidate
  //------------------------------------------------
  Int_t num_casc = AODEvent->GetNumberOfCascades();
  if (num_casc<=0) return;
  fHistEvents->Fill(6);
  fHNumberOfCasc->Fill(num_casc);

  // set primary vertex
  KFPVertex pVertex;
  Double_t pos[3],cov[6];
  fpVtx->GetXYZ(pos);
  if ( fabs(pos[2])>10. ) return; // vertex cut on z-axis direction
  fpVtx->GetCovarianceMatrix(cov);
  pVertex.SetXYZ((Float_t)pos[0], (Float_t)pos[1], (Float_t)pos[2]);
  Float_t covF[6];
  for (Int_t i=0; i<6; i++) { covF[i] = (Float_t)cov[i]; }
  pVertex.SetCovarianceMatrix(covF);
  pVertex.SetChi2(fpVtx->GetChi2());
  pVertex.SetNDF(fpVtx->GetNDF());
  pVertex.SetNContributors(fpVtx->GetNContributors());

  KFParticle PV(pVertex);

  if(!fAnaCuts) return;

  FillEventROOTObjects();

//------------------------------------------------
// Main analysis done in this function
//------------------------------------------------
  
  fPID = fInputHandler->GetPIDResponse();

  fCount_NumOfCandidatePerEvent_In3Sigma = 0;
  MakeAnaXicPlusFromCasc(AODEvent, mcArray, PV);
  fHNumOfCandidatePerEvent_In3sigma->Fill(fCount_NumOfCandidatePerEvent_In3Sigma);

  PostData(2, fCounter);
  PostData(3, fTree_Event);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
  PostData(4, fTree_XicPlus);
  PostData(5, fTree_XicPlusMCGen);
  PostData(7, fTree_XicPlus_QA);
  PostData(8, fTree_XicPlus_QA_woMassConstForLamAndXi);
  PostData(9, fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi);
  PostData(10, fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi);
  PostData(11, fTree_XicPlus_QA_wMassConstForLamAndXi);
  PostData(12, fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi);
  PostData(13, fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi);
  PostData(14, fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
/*
    TCanvas *c1 = new TCanvas();
    TLine *lLPDG = new TLine(7.89, 0, 7.89, 1e10);
    lLPDG->SetLineColor(2);
    lLPDG->Draw();

    TCanvas *c2 = new TCanvas();
    TLine *lXiPDG = new TLine(4.91, 0, 4.91, 1e10);
    lXiPDG->SetLineColor(2);
    lXiPDG->Draw();
*/
    return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MakeMCAnalysis(TClonesArray *mcArray)
{
  // Analyse AliAODMCParticle
  
  for(Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) {
    AliAODMCParticle *mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
    if (!mcPart) {
      AliError("Failed casting particle from MC array!, Skipping particle");
      continue;
    }
    Int_t pdg = mcPart->GetPdgCode();
    if ( TMath::Abs(pdg)!=4232 ) {
      AliDebug(2, Form("MC particle %d is not a Xic+: its pdg code is %d", iPart, pdg));
      continue;
    }
    AliDebug(2, Form("Step 0 ok: MC particle %d is a Xic+: its pdg code is %d", iPart, pdg));

    // 3 daughters
    if (mcPart->GetNDaughters()==3) {

    Int_t index_FirstDau = mcPart->GetDaughterFirst();
    if (index_FirstDau<0) continue;
    AliAODMCParticle* mcXicPlusDau_0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(index_FirstDau));
    AliAODMCParticle* mcXicPlusDau_1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(index_FirstDau+1));
    AliAODMCParticle* mcXicPlusDau_2 = dynamic_cast<AliAODMCParticle*>(mcArray->At(index_FirstDau+2));
    if ( TMath::Abs(mcXicPlusDau_0->GetPdgCode())!=211 && TMath::Abs(mcXicPlusDau_0->GetPdgCode())!=3312 ) continue;
    if ( TMath::Abs(mcXicPlusDau_1->GetPdgCode())!=211 && TMath::Abs(mcXicPlusDau_1->GetPdgCode())!=3312 ) continue;
    if ( TMath::Abs(mcXicPlusDau_2->GetPdgCode())!=211 && TMath::Abs(mcXicPlusDau_2->GetPdgCode())!=3312 ) continue;

    // -------------------------------- First daughter is Xi --------------------------------
    if ( TMath::Abs(mcXicPlusDau_0->GetPdgCode())==3312 && TMath::Abs(mcXicPlusDau_1->GetPdgCode())==211 && TMath::Abs(mcXicPlusDau_2->GetPdgCode())==211 ) { // First is Xi
      if (mcXicPlusDau_1->GetPdgCode()!=mcXicPlusDau_2->GetPdgCode()) continue; // Check pion charge is same
      if ( (mcXicPlusDau_0->GetPdgCode()>0&&mcXicPlusDau_1->GetPdgCode()<0) || (mcXicPlusDau_0->GetPdgCode()<0&&mcXicPlusDau_1->GetPdgCode()>0) ) continue; // Check charge of pion and Xi is different
      // Start to check Xi decay
      if ( mcXicPlusDau_0->GetNDaughters()!=2 ) continue;
      Bool_t pifromXi_flag = kFALSE;
      Bool_t v0_flag = kFALSE;
      Bool_t pifromLam_flag = kFALSE;
      Bool_t prfromLam_flag = kFALSE;
      AliAODMCParticle *mcv0part = NULL;
      for(Int_t idau=mcXicPlusDau_0->GetDaughterFirst();idau<=mcXicPlusDau_0->GetDaughterLast();idau++) {
        if(idau<0) break;
        AliAODMCParticle *mcdau = dynamic_cast<AliAODMCParticle*>(mcArray->At(idau));
        if ( TMath::Abs(mcdau->GetPdgCode())==211 && (idau!=(index_FirstDau+1)) && (idau!=(index_FirstDau+2)) ) { // 211: pion
          pifromXi_flag = kTRUE;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==3122 && mcdau->GetNDaughters()==2) { // 3122: Lambda
          v0_flag = kTRUE;
          mcv0part = mcdau;
          for(Int_t jdau=mcv0part->GetDaughterFirst();jdau<=mcv0part->GetDaughterLast();jdau++) {
            if (jdau<0) break;
            AliAODMCParticle *mcDau_Lam = (AliAODMCParticle*) mcArray->At(jdau);
            if ( TMath::Abs(mcDau_Lam->GetPdgCode())==211 && (jdau!=idau) && (jdau!=(index_FirstDau+1)) && (jdau!=(index_FirstDau+2)) ) pifromLam_flag = kTRUE;
            if(TMath::Abs(mcDau_Lam->GetPdgCode())==2212) prfromLam_flag = kTRUE;
          }
        }
      }
      if ( pifromXi_flag && v0_flag && pifromLam_flag && prfromLam_flag ) {
        if ( TMath::Abs(mcPart->Y()) < 0.8 ) { // Acceptance
          Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
          AliAODMCParticle *mcDau_0 = (AliAODMCParticle*) mcArray->At(mcPart->GetDaughterFirst());
          Double_t MLoverP = sqrt( pow(mcPart->Xv()-mcDau_0->Xv(),2.)+pow(mcPart->Yv()-mcDau_0->Yv(),2.)+pow(mcPart->Zv()-mcDau_0->Zv(),2.) ) * mcPart->M() / mcPart->P()*1.e4; // c*(proper lifetime) in um
          FillTreeGenXicPlus(mcPart, CheckOrigin, MLoverP);
        }
      }

      /*
      // Start to check Xi decay (method 2)
      if ( mcXicPlusDau_0->GetNDaughters()!=2 ) continue;
      Int_t label_XiDau_0 = mcXicPlusDau_0->GetDaughterLabel(0);
      Int_t label_XiDau_1 = mcXicPlusDau_0->GetDaughterLabel(1);
      AliAODMCParticle* mcXiDau_0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_XiDau_0));
      AliAODMCParticle* mcXiDau_1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_XiDau_1));
      if ( !mcXiDau_0 || !mcXiDau_1 ) {
        AliDebug(2, "Could not access Xi daughters, continuing...");
        continue;
      }
      Int_t pdgXiDau_0 = mcXiDau_0->GetPdgCode();
      Int_t pdgXiDau_1 = mcXiDau_1->GetPdgCode();
      if ( TMath::Abs(pdgXiDau_0)!=3122 || TMath::Abs(pdgXiDau_0)!=211 ) continue;
      if ( TMath::Abs(pdgXiDau_1)!=3122 || TMath::Abs(pdgXiDau_1)!=211 ) continue;
      if ( TMath::Abs(pdgXiDau_0)==3122 && TMath::Abs(pdgXiDau_1)==211) { // First is Lam
        if ( mcXiDau_0->GetNDaughters()!=2 ) continue;
        Int_t label_LamDau_0 = mcXiDau_0->GetDaughterLabel(0);
        Int_t label_LamDau_1 = mcXiDau_0->GetDaughterLabel(1);
        AliAODMCParticle* mcLamDau_0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_LamDau_0));
        AliAODMCParticle* mcLamDau_1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_LamDau_1));
        if ( !mcLamDau_0 || !mcLamDau_1 ) {
          AliDebug(2, "Could not access Lambda daughters, continuing...");
          continue;
        }
        Int_t pdgLamDau_0 = mcLamDau_0->GetPdgCode();
        Int_t pdgLamDau_1 = mcLamDau_1->GetPdgCode();
        if ( (TMath::Abs(pdgLamDau_0)==2212 && TMath::Abs(pdgLamDau_1)==211) || (TMath::Abs(pdgLamDau_1)==2212 && TMath::Abs(pdgLamDau_0)==211) ) {
          if ( TMath::Abs(mcPart->Y()) < 0.8 ) {
            Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
            FillTreeGenXicPlus(mcPart, CheckOrigin);
          }
        }
      } // First is Lam
      if ( TMath::Abs(pdgXiDau_1)==3122 && TMath::Abs(pdgXiDau_0)==211) { // Second is Lam
        if ( mcXiDau_1->GetNDaughters()!=2 ) continue;
        Int_t label_LamDau_0 = mcXiDau_1->GetDaughterLabel(0);
        Int_t label_LamDau_1 = mcXiDau_1->GetDaughterLabel(1);
        AliAODMCParticle* mcLamDau_0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_LamDau_0));
        AliAODMCParticle* mcLamDau_1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(label_LamDau_1));
        if ( !mcLamDau_0 || !mcLamDau_1 ) {
          AliDebug(2, "Could not access Lambda daughters, continuing...");
          continue;
        }
        Int_t pdgLamDau_0 = mcLamDau_0->GetPdgCode();
        Int_t pdgLamDau_1 = mcLamDau_1->GetPdgCode();
        if ( (TMath::Abs(pdgLamDau_0)==2212 && TMath::Abs(pdgLamDau_1)==211) || (TMath::Abs(pdgLamDau_1)==2212 && TMath::Abs(pdgLamDau_0)==211) ) {
          if ( TMath::Abs(mcPart->Y()) < 0.8 ) {
            Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
            FillTreeGenXicPlus(mcPart, CheckOrigin);
          }
        }
      } // Second is Lam
      */
    } // First is Xi

    // -------------------------------- Second daughter is Xi --------------------------------
    if ( TMath::Abs(mcXicPlusDau_1->GetPdgCode())==3312 && TMath::Abs(mcXicPlusDau_0->GetPdgCode())==211 && TMath::Abs(mcXicPlusDau_2->GetPdgCode())==211 ) { // Second is Xi
      if (mcXicPlusDau_0->GetPdgCode()!=mcXicPlusDau_2->GetPdgCode()) continue; // Check pion charge is same
      if ( (mcXicPlusDau_1->GetPdgCode()>0&&mcXicPlusDau_0->GetPdgCode()<0) || (mcXicPlusDau_1->GetPdgCode()<0&&mcXicPlusDau_0->GetPdgCode()>0) ) continue; // Check charge of pion and Xi is different
      // Start to check Xi decay
      if ( mcXicPlusDau_1->GetNDaughters()!=2 ) continue;
      Bool_t pifromXi_flag = kFALSE;
      Bool_t v0_flag = kFALSE;
      Bool_t pifromLam_flag = kFALSE;
      Bool_t prfromLam_flag = kFALSE;
      AliAODMCParticle *mcv0part = NULL;
      for(Int_t idau=mcXicPlusDau_1->GetDaughterFirst();idau<=mcXicPlusDau_1->GetDaughterLast();idau++) {
        if(idau<0) break;
        AliAODMCParticle *mcdau = dynamic_cast<AliAODMCParticle*>(mcArray->At(idau));
        if ( TMath::Abs(mcdau->GetPdgCode())==211 && (idau!=index_FirstDau) && (idau!=(index_FirstDau+2)) ) { // 211: pion
          pifromXi_flag = kTRUE;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==3122 && mcdau->GetNDaughters()==2) { // 3122: Lambda
          v0_flag = kTRUE;
          mcv0part = mcdau;
          for(Int_t jdau=mcv0part->GetDaughterFirst();jdau<=mcv0part->GetDaughterLast();jdau++) {
            if (jdau<0) break;
            AliAODMCParticle *mcDau_Lam = (AliAODMCParticle*) mcArray->At(jdau);
            if ( TMath::Abs(mcDau_Lam->GetPdgCode())==211 && (jdau!=idau) && (jdau!=index_FirstDau) && (jdau!=(index_FirstDau+2)) ) pifromLam_flag = kTRUE;
            if(TMath::Abs(mcDau_Lam->GetPdgCode())==2212) prfromLam_flag = kTRUE;
          }
        }
      }
      if ( pifromXi_flag && v0_flag && pifromLam_flag && prfromLam_flag ) {
        if ( TMath::Abs(mcPart->Y()) < 0.8 ) { // Acceptance
          Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
          AliAODMCParticle *mcDau_0 = (AliAODMCParticle*) mcArray->At(mcPart->GetDaughterFirst());
          Double_t MLoverP = sqrt( pow(mcPart->Xv()-mcDau_0->Xv(),2.)+pow(mcPart->Yv()-mcDau_0->Yv(),2.)+pow(mcPart->Zv()-mcDau_0->Zv(),2.) ) * mcPart->M() / mcPart->P()*1.e4; // c*(proper lifetime) in um
          FillTreeGenXicPlus(mcPart, CheckOrigin, MLoverP);
        }
      }
    } // Second is Xi

    // -------------------------------- Third daughter is Xi --------------------------------
    if ( TMath::Abs(mcXicPlusDau_2->GetPdgCode())==3312 && TMath::Abs(mcXicPlusDau_0->GetPdgCode())==211 && TMath::Abs(mcXicPlusDau_1->GetPdgCode())==211 ) { // Third is Xi
      if (mcXicPlusDau_0->GetPdgCode()!=mcXicPlusDau_1->GetPdgCode()) continue; // Check pion charge is same
      if ( (mcXicPlusDau_2->GetPdgCode()>0&&mcXicPlusDau_0->GetPdgCode()<0) || (mcXicPlusDau_2->GetPdgCode()<0&&mcXicPlusDau_0->GetPdgCode()>0) ) continue; // Check charge of pion and Xi is different
      // Start to check Xi decay
      if ( mcXicPlusDau_2->GetNDaughters()!=2 ) continue;
      Bool_t pifromXi_flag = kFALSE;
      Bool_t v0_flag = kFALSE;
      Bool_t pifromLam_flag = kFALSE;
      Bool_t prfromLam_flag = kFALSE;
      AliAODMCParticle *mcv0part = NULL;
      for(Int_t idau=mcXicPlusDau_2->GetDaughterFirst();idau<=mcXicPlusDau_2->GetDaughterLast();idau++) {
        if(idau<0) break;
        AliAODMCParticle *mcdau = dynamic_cast<AliAODMCParticle*>(mcArray->At(idau));
        if ( TMath::Abs(mcdau->GetPdgCode())==211 && (idau!=index_FirstDau) && (idau!=(index_FirstDau+1)) ) { // 211: pion
          pifromXi_flag = kTRUE;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==3122 && mcdau->GetNDaughters()==2) { // 3122: Lambda
          v0_flag = kTRUE;
          mcv0part = mcdau;
          for(Int_t jdau=mcv0part->GetDaughterFirst();jdau<=mcv0part->GetDaughterLast();jdau++) {
            if (jdau<0) break;
            AliAODMCParticle *mcDau_Lam = (AliAODMCParticle*) mcArray->At(jdau);
            if ( TMath::Abs(mcDau_Lam->GetPdgCode())==211 && (jdau!=idau) && (jdau!=index_FirstDau) && (jdau!=(index_FirstDau+1)) ) pifromLam_flag = kTRUE;
            if(TMath::Abs(mcDau_Lam->GetPdgCode())==2212) prfromLam_flag = kTRUE;
          }
        }
      }
      if ( pifromXi_flag && v0_flag && pifromLam_flag && prfromLam_flag ) {
        if ( TMath::Abs(mcPart->Y()) < 0.8 ) { // Acceptance
          Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
          AliAODMCParticle *mcDau_0 = (AliAODMCParticle*) mcArray->At(mcPart->GetDaughterFirst());
          Double_t MLoverP = sqrt( pow(mcPart->Xv()-mcDau_0->Xv(),2.)+pow(mcPart->Yv()-mcDau_0->Yv(),2.)+pow(mcPart->Zv()-mcDau_0->Zv(),2.) ) * mcPart->M() / mcPart->P()*1.e4; // c*(proper lifetime) in um
          FillTreeGenXicPlus(mcPart, CheckOrigin, MLoverP);
        }
      }
    } // Third is Xi
    } // 3 daughters
    if (mcPart->GetNDaughters()==2) {
      Bool_t PifromXicPlus_flag = kFALSE;
      Bool_t XiStar0_flag = kFALSE;
      Bool_t PifromXiStar0_flag = kFALSE;
      Bool_t Xi_flag = kFALSE;
      Bool_t PifromXi_flag = kFALSE;
      Bool_t Lam_flag = kFALSE;
      Bool_t PifromLam_flag = kFALSE;
      Bool_t PrfromLam_flag = kFALSE;
      for(Int_t idau=mcPart->GetDaughterFirst();idau<=mcPart->GetDaughterLast();idau++) {
        if (idau<0) break;
        AliAODMCParticle *mcDau_XicPlus = dynamic_cast<AliAODMCParticle*>(mcArray->At(idau));
        if (TMath::Abs(mcDau_XicPlus->GetPdgCode())==211) PifromXicPlus_flag = kTRUE;
        if (TMath::Abs(mcDau_XicPlus->GetPdgCode())==3324 && mcDau_XicPlus->GetNDaughters()==2) { // Xi*0
          XiStar0_flag = kTRUE;
          for (Int_t jdau=mcDau_XicPlus->GetDaughterFirst();jdau<=mcDau_XicPlus->GetDaughterLast();jdau++) { // Xi*0 daughters
            if (jdau<0) break;
            AliAODMCParticle *mcDau_XiStar0 = (AliAODMCParticle*) mcArray->At(jdau);
            if (TMath::Abs(mcDau_XiStar0->GetPdgCode())==211) PifromXiStar0_flag = kTRUE;
            if (TMath::Abs(mcDau_XiStar0->GetPdgCode())==3312 && mcDau_XiStar0->GetNDaughters()==2) { // Xi-
              Xi_flag = kTRUE;
              for (Int_t kdau=mcDau_XiStar0->GetDaughterFirst();kdau<=mcDau_XiStar0->GetDaughterLast();kdau++) { // Xi- daughters
                if (kdau<0) break;
                AliAODMCParticle *mcDau_Xi = (AliAODMCParticle*) mcArray->At(kdau);
                if (TMath::Abs(mcDau_Xi->GetPdgCode())==211) PifromXi_flag = kTRUE;
                if (TMath::Abs(mcDau_Xi->GetPdgCode())==3122 && mcDau_Xi->GetNDaughters()==2) { // Lambda
                  Lam_flag = kTRUE;
                  for(Int_t ldau=mcDau_Xi->GetDaughterFirst();ldau<=mcDau_Xi->GetDaughterLast();ldau++) {
                    if (ldau<0) break;
                    AliAODMCParticle *mcDau_Lam = (AliAODMCParticle*) mcArray->At(ldau);
                    if(TMath::Abs(mcDau_Lam->GetPdgCode())==211) PifromLam_flag = kTRUE;
                    if(TMath::Abs(mcDau_Lam->GetPdgCode())==2212) PrfromLam_flag = kTRUE;
                  }
                }
              }
            }
          }
        }
      }
      if ( PifromXicPlus_flag && XiStar0_flag && PifromXiStar0_flag && Xi_flag && PifromXi_flag && Lam_flag && PifromLam_flag && PrfromLam_flag ) {
        if ( TMath::Abs(mcPart->Y()) < 0.8 ) { // Acceptance
          Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcPart,kTRUE);
          if (CheckOrigin==4) CheckOrigin=-4;
          if (CheckOrigin==5) CheckOrigin=-5;
          AliAODMCParticle *mcDau_0 = (AliAODMCParticle*) mcArray->At(mcPart->GetDaughterFirst());
          Double_t MLoverP = sqrt( pow(mcPart->Xv()-mcDau_0->Xv(),2.)+pow(mcPart->Yv()-mcDau_0->Yv(),2.)+pow(mcPart->Zv()-mcDau_0->Zv(),2.) ) * mcPart->M() / mcPart->P()*1.e4; // c*(proper lifetime) in um
          if (fWriteXicPlusMCGenTree) FillTreeGenXicPlus(mcPart, CheckOrigin, MLoverP);
        }
      }
    } // 2 daughters
  }

  return kTRUE;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillTreeGenXicPlus(AliAODMCParticle *mcpart, Int_t CheckOrigin, Double_t MLoverP)
{
  // Fill histograms or tree depending

  for(Int_t i=0;i<5;i++){
    fVar_XicPlusMCGen[i] = -9999.;
  }

  fVar_XicPlusMCGen[0] = mcpart->Y();
  fVar_XicPlusMCGen[1] = mcpart->Pt();
  fVar_XicPlusMCGen[2] = CheckOrigin;
  fVar_XicPlusMCGen[3] = mcpart->GetPdgCode();
  fVar_XicPlusMCGen[4] = MLoverP;

  if (fVar_XicPlusMCGen[1]>0.9999) fTree_XicPlusMCGen->Fill();

//  fVar_XicPlusMCGen[ 0] = fCentrality;
//  fVar_XicPlusMCGen[ 1] = decaytype;
//  if (mcpart->IsPrimary() && (!mcpart->IsPhysicalPrimary())) fVar_XicPlusMCGen[2] = 1;
//  if (mcpart->IsPhysicalPrimary()) fVar_XicPlusMCGen[2] = 2;
//  if (mcpart->IsSecondaryFromWeakDecay()) fVar_XicPlusMCGen[2] = 3;
//  if (mcpart->IsSecondaryFromMaterial()) fVar_XicPlusMCGen[2] = 4;
//  if (mcpart->IsFromSubsidiaryEvent()) fVar_XicPlusMCGen[2] = 5;
//  fVar_XicPlusMCGen[ 3] = mcpart->Eta();
//  fVar_XicPlusMCGen[ 4] = mcpart->Y();
//  fVar_XicPlusMCGen[ 5] = mcpart->Px();
//  fVar_XicPlusMCGen[ 6] = mcpart->Py();
//  fVar_XicPlusMCGen[ 7] = mcpart->Pz();
//  fVar_XicPlusMCGen[ 8] = mcpipart->Px();
//  fVar_XicPlusMCGen[ 9] = mcpipart->Py();
//  fVar_XicPlusMCGen[10] = mcpipart->Pz();
//  fVar_XicPlusMCGen[11] = mccascpart->Px();
//  fVar_XicPlusMCGen[12] = mccascpart->Py();
//  fVar_XicPlusMCGen[13] = mccascpart->Pz();
//  fVar_XicPlusMCGen[14] = mcpart->GetPdgCode();
//  fVar_XicPlusMCGen[15] = mcpipart->GetPdgCode();
//  fVar_XicPlusMCGen[16] = mccascpart->GetPdgCode();
//  fVar_XicPlusMCGen[17] = fRunNumber;
//  fVar_XicPlusMCGen[18] = fEvNumberCounter;
//  fVar_XicPlusMCGen[19] = CheckOrigin;

  /*
  const Double_t massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t massXi   = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  Double_t pipx = mcpipart->Px();
  Double_t pipy = mcpipart->Py();
  Double_t pipz = mcpipart->Pz();
//  Double_t piE  = sqrt(pipx*pipx+pipy*pipy+pipz*pipz+0.000511*0.000511);
  Double_t piE  = sqrt(pipx*pipx+pipy*pipy+pipz*pipz+massPion*massPion);
  Double_t cascpx = mccascpart->Px();
  Double_t cascpy = mccascpart->Py();
  Double_t cascpz = mccascpart->Pz();
  Double_t cascE  = sqrt(cascpx*cascpx+cascpy*cascpy+cascpz*cascpz+massXi*massXi);

  Double_t InvMassPiXi = sqrt(pow(piE+cascE,2)-pow(pipx+cascpx,2)-pow(pipy+cascpy,2)-pow(pipz+cascpz,2));

  Double_t contXicPlusMC[3];
  contXicPlusMC[0] = mcpart->Pt();
  contXicPlusMC[1] = mcpart->Y();
  contXicPlusMC[2] = fCentrality;

  Double_t contPionMC[3];
  contPionMC[0] = mcpipart->Pt();
  contPionMC[1] = mcpipart->Eta();
  contPionMC[2] = fCentrality;

  Double_t contXiMC[3];
  contXiMC[0] = mccascpart->Pt();
  contXiMC[1] = mccascpart->Y();
  contXiMC[2] = fCentrality;

  Double_t contPiXiMassMCGen[3];
  contPiXiMassMCGen[0] = InvMassPiXi;
  contPiXiMassMCGen[1] = mcpart->Pt();
  contPiXiMassMCGen[2] = fCentrality;

  Double_t contPiXiMassvsPiPtMCGen[3];
  contPiXiMassvsPiPtMCGen[0] = InvMassPiXi;
  contPiXiMassvsPiPtMCGen[1] = mcpipart->Pt();
  contPiXiMassvsPiPtMCGen[2] = fCentrality;

  if (decaytype==0) {
    if (fabs(mcpipart->Eta())<fAnaCuts->GetProdTrackEtaRange()) {
    }
  }
  */
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::MakeAnaXicPlusFromCasc(AliAODEvent *AODEvent, TClonesArray *mcArray, KFParticle PV)
{
  // Main analysis called from "UserExec"

//  std::cout.setf(std::ios::fixed);
//  std::cout.setf(std::ios::showpoint);
//  std::cout.precision(3);

  // set the magnetic field
  KFParticle::SetField(fBzkG);

  const UInt_t nCasc = AODEvent->GetNumberOfCascades();

  Double_t xyzP[3], xyzN[3];
  Double_t xvyvzvP[3], xvyvzvN[3];
  Double_t covP[21], covN[21], covB[21];
  const Float_t massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  const Float_t massXi     = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  const Float_t massK0S    = TDatabasePDG::Instance()->GetParticle(310)->Mass();

  // select good candidates for pion
  const UInt_t nTracks = AODEvent->GetNumberOfTracks();
  AliAODTrack *trackP[nTracks], *trackN[nTracks];
  Int_t flag_trkP = 0, flag_trkN = 0;

  Int_t count_ForPV = 0;
  Int_t count_ForPV_AfterQAcheck = 0;
  Int_t count_ForPV_wPosID = 0;
  Int_t Index_AODtrkUsedForPVfit[999] = {-999};

  for (UInt_t itrk=0; itrk<nTracks; itrk++) {
    AliAODTrack *trk = static_cast<AliAODTrack*>(AODEvent->GetTrack(itrk));

    if (trk->GetUsedForPrimVtxFit()) count_ForPV++;

    if (trk->GetID()>=-0.00001 && trk->GetUsedForPrimVtxFit()) count_ForPV_wPosID++;

    Double_t covtest[21];
    if ( !trk || trk->GetID()<0 || !trk->GetCovarianceXYZPxPyPz(covtest) || !AliVertexingHFUtils::CheckAODtrackCov(trk) ) continue;

    if  (trk->Charge() > 0 ) {
      trackP[flag_trkP] = trk;
      flag_trkP++;
    }
    if  (trk->Charge() < 0 ) {
      trackN[flag_trkN] = trk;
      flag_trkN++;
    }
    if (trk->GetUsedForPrimVtxFit()) {
      Index_AODtrkUsedForPVfit[count_ForPV_AfterQAcheck] = itrk;
      if (fIsMC) {
        AliAODMCParticle* mcTrk = static_cast<AliAODMCParticle*>(mcArray->At(fabs(trk->GetLabel())));
        //cout << "PDG: " << mcTrk->GetPdgCode() << endl;
        AliAODMCParticle* mcTrk_Mother = static_cast<AliAODMCParticle*>(mcArray->At(mcTrk->GetMother()));
        //cout << "Mother PDG: " << mcTrk_Mother->GetPdgCode() << endl;
      }
      count_ForPV_AfterQAcheck++;
    }
  }

  // Rebuild PV with KF
  KFVertex PV_KF_Refit;
  if ( strstr(fpVtx->GetTitle(),"VertexerTracksWithConstraint") || strstr(fpVtx->GetTitle(),"VertexerTracksMVWithConstraint") ) {
    Float_t diamondcovxy[3], sigma2DiamondZ;
    AODEvent->GetDiamondCovXY(diamondcovxy);
    sigma2DiamondZ = AODEvent->GetSigma2DiamondZ();
    Double_t pos[3]={AODEvent->GetDiamondX(),AODEvent->GetDiamondY(),AODEvent->GetDiamondZ()};
  	Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,sigma2DiamondZ};
    PV_KF_Refit.SetBeamConstraint(pos[0], pos[1], pos[2], sqrt(cov[0]), sqrt(cov[2]), sqrt(cov[5]));
  } else {
    PV_KF_Refit.SetBeamConstraintOff();
  }
  Bool_t vtxFlag[count_ForPV_AfterQAcheck];
  KFParticle PVdau_tmp[count_ForPV_AfterQAcheck];
  const KFParticle **PV_dau = new const KFParticle*[count_ForPV_AfterQAcheck+1];
  for (UInt_t i=0; i<count_ForPV_AfterQAcheck; i++) {
    AliAODTrack *trk = static_cast<AliAODTrack*>(AODEvent->GetTrack(Index_AODtrkUsedForPVfit[i]));
    PVdau_tmp[i] = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trk, 0);
    PV_dau[i] = &PVdau_tmp[i];
  }
  PV_KF_Refit.ConstructPrimaryVertex(PV_dau, count_ForPV_AfterQAcheck, vtxFlag, 1.e9);
  delete [] PV_dau;

  /*
  cout << "===========================" << endl;
  cout << "PV_KF_Refit.GetNContributors: " << PV_KF_Refit.GetNContributors()+2 << endl; // bug in KF, should be added by 2
  cout << "count_ForPV_AfterQAcheck: " << count_ForPV_AfterQAcheck << endl;
  cout << "===========================" << endl;
  */

  /*
  cout << "===========================" << endl;
  cout << "count_ForPV_AfterQAcheck: " << count_ForPV_AfterQAcheck << endl;
  cout << "PV_KF_Refit(GetNContributors): " << PV_KF_Refit.GetNContributors() << endl;
  cout << "PV_KF_Refit(X): " << PV_KF_Refit.GetX() << endl;
  cout << "PV(X): " << PV.GetX() << endl;
  cout << "---------------------------" << endl;
  cout << "PV (GetNContributors): " << fpVtx->GetNContributors() << endl;
  cout << "PV (GetNDaughters): " << fpVtx->GetNDaughters() << endl;
  cout << "PV (CountRealContributors): " << fpVtx->CountRealContributors() << endl;
  cout << "Number of AOD track used for PV fit: " << count_ForPV << endl;
  cout << "Number of AOD track used for PV fit (w/ positive ID): " << count_ForPV_wPosID << endl;
  cout << "Number of AOD track used for PV fit (w/ positive ID and good covariance): " << count_ForPV_AfterQAcheck << endl;
  cout << "===========================" << endl;
  */

  for (UInt_t iCasc=0; iCasc<nCasc; iCasc++) {
    AliAODcascade *casc = AODEvent->GetCascade(iCasc);

    // cascade cut
    if ( !fAnaCuts->SingleCascCuts(casc, kFALSE) ) continue;

    AliAODTrack *ptrack = (AliAODTrack*) (casc->GetDaughter(0));
    AliAODTrack *ntrack = (AliAODTrack*) (casc->GetDaughter(1));
    AliAODTrack *btrack = (AliAODTrack*) (casc->GetDecayVertexXi()->GetDaughter(0));

    if ( !ptrack||!ntrack||!btrack ) continue;

    /*
    if (nCasc>1) {
      cout << "~~~~~~~~~~~~" << endl;
      AliAODMCParticle* mcCasc_Dau = static_cast<AliAODMCParticle*>(mcArray->At(fabs(btrack->GetLabel())));
      AliAODMCParticle* mcCasc = static_cast<AliAODMCParticle*>(mcArray->At(mcCasc_Dau->GetMother()));
      AliAODMCParticle* mcCasc_Mother = static_cast<AliAODMCParticle*>(mcArray->At(mcCasc->GetMother()));
      AliAODMCParticle* mcCasc_MoMother = static_cast<AliAODMCParticle*>(mcArray->At(mcCasc_Mother->GetMother()));
      AliAODMCParticle* mc_ptrack = static_cast<AliAODMCParticle*>(mcArray->At(fabs(ptrack->GetLabel())));
      AliAODMCParticle* mc_ntrack = static_cast<AliAODMCParticle*>(mcArray->At(fabs(ntrack->GetLabel())));
      if ( fabs(mcCasc_Mother->GetPdgCode())==4232 ) {
      cout << "Casc_MoMother PDG: " << mcCasc_MoMother->GetPdgCode() << ", Index: " << mcCasc_Mother->GetMother() << endl;
      cout << "Casc_Mother PDG: " << mcCasc_Mother->GetPdgCode() << ", Index: " << mcCasc->GetMother() << endl;
      cout << "btrack PDG: " << mcCasc_Dau->GetPdgCode() << ", ID: " << btrack->GetID() << endl;
      cout << "ptrack PDG: " << mc_ptrack->GetPdgCode() << ", ID: " << ptrack->GetID() << endl;
      cout << "ntrack PDG: " << mc_ntrack->GetPdgCode() << ", ID: " << ntrack->GetID() << endl;
      }
    }
    */

    // check charge of the first daughter, if negative, define it as the second one
    if ( ptrack->Charge()<0 ) {
      ptrack = (AliAODTrack*) (casc->GetDaughter(1));
      ntrack = (AliAODTrack*) (casc->GetDaughter(0));
    }

    if ( !ptrack->GetCovarianceXYZPxPyPz(covP) || !ntrack->GetCovarianceXYZPxPyPz(covN) || !btrack->GetCovarianceXYZPxPyPz(covB) ) continue;

    if ( !AliVertexingHFUtils::CheckAODtrackCov(ptrack) || !AliVertexingHFUtils::CheckAODtrackCov(ntrack) || !AliVertexingHFUtils::CheckAODtrackCov(btrack) ) continue;

    KFParticle kfpProton     = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ptrack, 2212);
    KFParticle kfpPionMinus  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ntrack, -211);
    KFParticle kfpAntiProton = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ntrack, -2212);
    KFParticle kfpPionPlus   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ptrack, 211);

    KFParticle kfpElePlus    = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ptrack, -11);
    KFParticle kfpEleMinus   = AliVertexingHFUtils::CreateKFParticleFromAODtrack(ntrack, 11);

    // === K0S ===
    KFParticle kfpK0Short;
    const KFParticle *vk0sDaughters[2]  = {&kfpPionPlus, &kfpPionMinus};
    kfpK0Short.Construct(vk0sDaughters, 2);
    // ============
    // === Gamma ===
    KFParticle kfpGamma;
    const KFParticle *vGammaDaughters[2]  = {&kfpElePlus, &kfpEleMinus};
    kfpGamma.Construct(vGammaDaughters, 2);
    // =============

    if ( btrack->Charge()<0 ) { // Xi^-

      const KFParticle *vDaughters[2] = {&kfpProton, &kfpPionMinus};

      KFParticle kfpLambda;
      kfpLambda.Construct(vDaughters, 2);
      Float_t massLambda_Rec, err_massLambda_Rec;
      kfpLambda.GetMass(massLambda_Rec, err_massLambda_Rec);

      // check rapidity of lambda
      if ( TMath::Abs(kfpLambda.GetE())<=TMath::Abs(kfpLambda.GetPz()) ) continue;

      // chi2>0 && NDF>0 for selecting Lambda
      if ( (kfpLambda.GetNDF()<=0 || kfpLambda.GetChi2()<=0) ) continue;

      // check cov. of Lambda
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLambda) ) continue;

      // err_mass_Rec>0 of Lambda
      if ( err_massLambda_Rec<=0 ) continue;

      // Chi2geo cut of Lambda
      if ( (kfpLambda.GetChi2()/kfpLambda.GetNDF()) >= fAnaCuts->GetKFPLam_Chi2geoMax() ) continue; 

      //************************** calculate l/Δl for Lambda *************************************
      Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
      Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
      Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
      Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
      Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
      if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
      dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
      Double_t nErr_l_Lambda = l_Lambda/dl_Lambda;
      //***************************************************************************************

      // l/Deltal cut of Lambda
      if ( nErr_l_Lambda <= fAnaCuts->GetKFPLam_lDeltalMin() ) continue;

      // mass window cut of Lambda
      if ( TMath::Abs(massLambda_Rec-massLambda) > (fAnaCuts->GetProdMassTolLambda()) ) continue;

      // ================================ QA study (2022.09.07) ================================
      if (fWriteXicPlusQATree) {
        KFParticle kfpPion_ForXi;
        kfpPion_ForXi = AliVertexingHFUtils::CreateKFParticleFromAODtrack(btrack, -211); // pion-

        // reconstruct Xi-
        KFParticle kfpXiMinus_woMassConstForLamAndXi;
        const KFParticle *vXiDs_woMassConstForLamAndXi[2] = {&kfpPion_ForXi, &kfpLambda};
        kfpXiMinus_woMassConstForLamAndXi.Construct(vXiDs_woMassConstForLamAndXi, 2);
        Float_t massXiMinus_Rec_woMassConstForLamAndXi, err_massXiMinus_Rec_woMassConstForLamAndXi;
        kfpXiMinus_woMassConstForLamAndXi.GetMass(massXiMinus_Rec_woMassConstForLamAndXi, err_massXiMinus_Rec_woMassConstForLamAndXi);



        ////// Without mass constraint for Lambda
        if ( TMath::Abs(kfpXiMinus_woMassConstForLamAndXi.GetE())>TMath::Abs(kfpXiMinus_woMassConstForLamAndXi.GetPz()) && // check rapidity of Xi-
             err_massXiMinus_Rec_woMassConstForLamAndXi>0 && // err_massXi > 0
             kfpXiMinus_woMassConstForLamAndXi.GetNDF()>0 && kfpXiMinus_woMassConstForLamAndXi.GetChi2()>0 && // chi2>0 && NDF>0
             AliVertexingHFUtils::CheckKFParticleCov(kfpXiMinus_woMassConstForLamAndXi) && // check covariance matrix
             kfpXiMinus_woMassConstForLamAndXi.GetChi2()/kfpXiMinus_woMassConstForLamAndXi.GetNDF() < fAnaCuts->GetKFPXi_Chi2geoMax() && // Prefilter
             TMath::Abs(massXiMinus_Rec_woMassConstForLamAndXi-massXi) <= (fAnaCuts->GetProdMassTolXi()) // mass window cut of Xi-
           ) {
          ////// Without mass constraint for Xi
          // Loop for bachelor pions
          for (Int_t itrkBP_trk1=0; itrkBP_trk1<(flag_trkP-1); itrkBP_trk1++) { // Loop for first bachelor pion+
            for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkP; itrkBP_trk2++) { // Loop for second bachelor pion+
              if ( trackP[itrkBP_trk1]->GetID() == ptrack->GetID() || trackP[itrkBP_trk2]->GetID() == ptrack->GetID() || trackP[itrkBP_trk1]->GetID() == trackP[itrkBP_trk2]->GetID() ) continue;
              if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk1]) ) continue;
              if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk2]) ) continue;

              // === pt(pion_0)<pt(pion_1) ===
              AliAODTrack *trackPiFromXicPlus_LowPt = NULL, *trackPiFromXicPlus_HighPt = NULL;
              trackPiFromXicPlus_LowPt  = trackP[itrkBP_trk1];
              trackPiFromXicPlus_HighPt = trackP[itrkBP_trk2];
              if (trackP[itrkBP_trk1]->Pt() > trackP[itrkBP_trk2]->Pt()) {
                trackPiFromXicPlus_LowPt  = trackP[itrkBP_trk2];
                trackPiFromXicPlus_HighPt = trackP[itrkBP_trk1];
              }
              // =============================
              KFParticle kfpBP_LowPt  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicPlus_LowPt, 211);
              KFParticle kfpBP_HighPt = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicPlus_HighPt, 211);

              // reconstruct Xic+
              KFParticle kfpXicPlus_woMassConstForLamAndXi;
              const KFParticle *vXicPlusDs_woMassConstForLamAndXi[3] = {&kfpXiMinus_woMassConstForLamAndXi, &kfpBP_LowPt, &kfpBP_HighPt};
              kfpXicPlus_woMassConstForLamAndXi.Construct(vXicPlusDs_woMassConstForLamAndXi, 3);

              Float_t massXicPlus_Rec_woMassConstForLamAndXi, err_massXicPlus_Rec_woMassConstForLamAndXi;
              kfpXicPlus_woMassConstForLamAndXi.GetMass(massXicPlus_Rec_woMassConstForLamAndXi, err_massXicPlus_Rec_woMassConstForLamAndXi);

              if ( kfpXicPlus_woMassConstForLamAndXi.GetNDF()>0 && kfpXicPlus_woMassConstForLamAndXi.GetChi2()>0 && // chi2>0 && NDF>0
                   TMath::Abs(kfpXicPlus_woMassConstForLamAndXi.GetE())>TMath::Abs(kfpXicPlus_woMassConstForLamAndXi.GetPz()) && // check rapidity of XicPlus
                   AliVertexingHFUtils::CheckKFParticleCov(kfpXicPlus_woMassConstForLamAndXi) && // check covariance matrix
                   err_massXicPlus_Rec_woMassConstForLamAndXi>0 && // err_massXicPlus > 0
                   kfpXicPlus_woMassConstForLamAndXi.GetChi2()/kfpXicPlus_woMassConstForLamAndXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                   kfpXicPlus_woMassConstForLamAndXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                 ) {
                Int_t lab_XicPlus = -9999;
                if (fIsMC) lab_XicPlus = MatchToMCXicPlus(ptrack, ntrack, btrack, trackP[itrkBP_trk1], trackP[itrkBP_trk2], mcArray);
                FillQATreeXicPlusFromCasc_woMassConstForLamAndXi(kfpLambda, kfpXiMinus_woMassConstForLamAndXi, kfpXicPlus_woMassConstForLamAndXi, PV, trackPiFromXicPlus_HighPt, mcArray, lab_XicPlus, AODEvent);
              }
              kfpXicPlus_woMassConstForLamAndXi.Clear();
              kfpBP_HighPt.Clear();
              kfpBP_LowPt.Clear();
            }
          }

          ////// (With mass Constraint for Xi) + (without mass constraint for Lambda)
          KFParticle kfpXiMinus_woMassConstForLam_wMassConstForXi = kfpXiMinus_woMassConstForLamAndXi;
          kfpXiMinus_woMassConstForLam_wMassConstForXi.SetNonlinearMassConstraint(massXi);
          if ( AliVertexingHFUtils::CheckKFParticleCov(kfpXiMinus_woMassConstForLam_wMassConstForXi) && TMath::Abs(kfpXiMinus_woMassConstForLam_wMassConstForXi.GetE()) > TMath::Abs(kfpXiMinus_woMassConstForLam_wMassConstForXi.GetPz()) ) {
            // Loop for bachelor pions
            for (Int_t itrkBP_trk1=0; itrkBP_trk1<(flag_trkP-1); itrkBP_trk1++) { // Loop for first bachelor pion+
              for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkP; itrkBP_trk2++) { // Loop for second bachelor pion+
                if ( trackP[itrkBP_trk1]->GetID() == ptrack->GetID() || trackP[itrkBP_trk2]->GetID() == ptrack->GetID() || trackP[itrkBP_trk1]->GetID() == trackP[itrkBP_trk2]->GetID() ) continue;

                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk1]) ) continue;
                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk2]) ) continue;

                // === pt(pion_0)<pt(pion_1) ===
                AliAODTrack *trackPiFromXicPlus_LowPt = NULL, *trackPiFromXicPlus_HighPt = NULL;
                trackPiFromXicPlus_LowPt  = trackP[itrkBP_trk1];
                trackPiFromXicPlus_HighPt = trackP[itrkBP_trk2];
                if (trackP[itrkBP_trk1]->Pt() > trackP[itrkBP_trk2]->Pt()) {
                  trackPiFromXicPlus_LowPt  = trackP[itrkBP_trk2];
                  trackPiFromXicPlus_HighPt = trackP[itrkBP_trk1];
                }
                // =============================
                KFParticle kfpBP_LowPt  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicPlus_LowPt, 211);
                KFParticle kfpBP_HighPt = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicPlus_HighPt, 211);

                // reconstruct Xic+
                KFParticle kfpXicPlus_woMassConstForLam_wMassConstForXi;
                const KFParticle *vXicPlusDs_woMassConstForLam_wMassConstForXi[3] = {&kfpXiMinus_woMassConstForLam_wMassConstForXi, &kfpBP_LowPt, &kfpBP_HighPt};
                kfpXicPlus_woMassConstForLam_wMassConstForXi.Construct(vXicPlusDs_woMassConstForLam_wMassConstForXi, 3);

                Float_t massXicPlus_Rec_woMassConstForLam_wMassConstForXi, err_massXicPlus_Rec_woMassConstForLam_wMassConstForXi;
                kfpXicPlus_woMassConstForLam_wMassConstForXi.GetMass(massXicPlus_Rec_woMassConstForLam_wMassConstForXi, err_massXicPlus_Rec_woMassConstForLam_wMassConstForXi);

                if ( kfpXicPlus_woMassConstForLam_wMassConstForXi.GetNDF()>0 && kfpXicPlus_woMassConstForLam_wMassConstForXi.GetChi2()>0 && // chi2>0 && NDF>0
                     TMath::Abs(kfpXicPlus_woMassConstForLam_wMassConstForXi.GetE())>TMath::Abs(kfpXicPlus_woMassConstForLam_wMassConstForXi.GetPz()) && // check rapidity of XicPlus
                     AliVertexingHFUtils::CheckKFParticleCov(kfpXicPlus_woMassConstForLam_wMassConstForXi) && // check covariance matrix
                     err_massXicPlus_Rec_woMassConstForLam_wMassConstForXi>0 && // err_massXicPlus > 0
                     kfpXicPlus_woMassConstForLam_wMassConstForXi.GetChi2()/kfpXicPlus_woMassConstForLam_wMassConstForXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                     kfpXicPlus_woMassConstForLam_wMassConstForXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                   ) {
                  Int_t lab_XicPlus = -9999;
                  if (fIsMC) lab_XicPlus = MatchToMCXicPlus(ptrack, ntrack, btrack, trackP[itrkBP_trk1], trackP[itrkBP_trk2], mcArray);
                  FillQATreeXicPlusFromCasc_woMassConstForLam_wMassConstForXi(kfpLambda, kfpXiMinus_woMassConstForLam_wMassConstForXi, kfpXicPlus_woMassConstForLam_wMassConstForXi, PV, trackPiFromXicPlus_HighPt, mcArray, lab_XicPlus, AODEvent);
                }
                kfpXicPlus_woMassConstForLam_wMassConstForXi.Clear();
                kfpBP_HighPt.Clear();
                kfpBP_LowPt.Clear();
              }
            }
          }
          kfpXiMinus_woMassConstForLam_wMassConstForXi.Clear();
          kfpXiMinus_woMassConstForLamAndXi.Clear();
        }



        ////// With mass constraint for Lambda
        // Mass Constraint for Lambda
        KFParticle kfpLambda_wMassConst = kfpLambda;
        kfpLambda_wMassConst.SetNonlinearMassConstraint(massLambda);
        if ( AliVertexingHFUtils::CheckKFParticleCov(kfpLambda_wMassConst) && TMath::Abs(kfpLambda_wMassConst.GetE()) > TMath::Abs(kfpLambda_wMassConst.GetPz()) ) {
          KFParticle kfpXiMinus_wMassConstForLam_woMassConstForXi;
          const KFParticle *vXiDs_wMassConstForLam_woMassConstForXi[2] = {&kfpPion_ForXi, &kfpLambda_wMassConst};
          kfpXiMinus_wMassConstForLam_woMassConstForXi.Construct(vXiDs_wMassConstForLam_woMassConstForXi, 2);

          Float_t massXiMinus_Rec_wMassConstForLam_woMassConstForXi, err_massXiMinus_Rec_wMassConstForLam_woMassConstForXi;
          kfpXiMinus_wMassConstForLam_woMassConstForXi.GetMass(massXiMinus_Rec_wMassConstForLam_woMassConstForXi, err_massXiMinus_Rec_wMassConstForLam_woMassConstForXi);

          if ( TMath::Abs(kfpXiMinus_wMassConstForLam_woMassConstForXi.GetE())>TMath::Abs(kfpXiMinus_wMassConstForLam_woMassConstForXi.GetPz()) && // check rapidity of Xi-
               err_massXiMinus_Rec_wMassConstForLam_woMassConstForXi>0 && // err_massXi > 0
               kfpXiMinus_wMassConstForLam_woMassConstForXi.GetNDF()>0 && kfpXiMinus_wMassConstForLam_woMassConstForXi.GetChi2()>0 && // chi2>0 && NDF>0
               AliVertexingHFUtils::CheckKFParticleCov(kfpXiMinus_wMassConstForLam_woMassConstForXi) && // check covariance matrix
               kfpXiMinus_wMassConstForLam_woMassConstForXi.GetChi2()/kfpXiMinus_wMassConstForLam_woMassConstForXi.GetNDF() < fAnaCuts->GetKFPXi_Chi2geoMax() && // Prefilter
               TMath::Abs(massXiMinus_Rec_wMassConstForLam_woMassConstForXi-massXi) <= (fAnaCuts->GetProdMassTolXi()) // mass window cut of Xi-
             ) {
            // Loop for bachelor pions
            for (Int_t itrkBP_trk1=0; itrkBP_trk1<(flag_trkP-1); itrkBP_trk1++) { // Loop for first bachelor pion+
              for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkP; itrkBP_trk2++) { // Loop for second bachelor pion+
                if ( trackP[itrkBP_trk1]->GetID() == ptrack->GetID() || trackP[itrkBP_trk2]->GetID() == ptrack->GetID() || trackP[itrkBP_trk1]->GetID() == trackP[itrkBP_trk2]->GetID() ) continue;
                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk1]) ) continue;
                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk2]) ) continue;

                // === pt(pion_0)<pt(pion_1) ===
                AliAODTrack *trackPiFromXicPlus_LowPt = NULL, *trackPiFromXicPlus_HighPt = NULL;
                trackPiFromXicPlus_LowPt  = trackP[itrkBP_trk1];
                trackPiFromXicPlus_HighPt = trackP[itrkBP_trk2];

                if (trackP[itrkBP_trk1]->Pt() > trackP[itrkBP_trk2]->Pt()) {
                  trackPiFromXicPlus_LowPt  = trackP[itrkBP_trk2];
                  trackPiFromXicPlus_HighPt = trackP[itrkBP_trk1];
                }
                // =============================
                KFParticle kfpBP_LowPt  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicPlus_LowPt, 211);
                KFParticle kfpBP_HighPt = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicPlus_HighPt, 211);

                // reconstruct Xic+
                KFParticle kfpXicPlus_wMassConstForLam_woMassConstForXi;
                const KFParticle *vXicPlusDs_wMassConstForLam_woMassConstForXi[3] = {&kfpXiMinus_wMassConstForLam_woMassConstForXi, &kfpBP_LowPt, &kfpBP_HighPt};
                kfpXicPlus_wMassConstForLam_woMassConstForXi.Construct(vXicPlusDs_wMassConstForLam_woMassConstForXi, 3);

                Float_t massXicPlus_Rec_wMassConstForLam_woMassConstForXi, err_massXicPlus_Rec_wMassConstForLam_woMassConstForXi;
                kfpXicPlus_wMassConstForLam_woMassConstForXi.GetMass(massXicPlus_Rec_wMassConstForLam_woMassConstForXi, err_massXicPlus_Rec_wMassConstForLam_woMassConstForXi);

                if ( kfpXicPlus_wMassConstForLam_woMassConstForXi.GetNDF()>0 && kfpXicPlus_wMassConstForLam_woMassConstForXi.GetChi2()>0 && // chi2>0 && NDF>0
                     TMath::Abs(kfpXicPlus_wMassConstForLam_woMassConstForXi.GetE())>TMath::Abs(kfpXicPlus_wMassConstForLam_woMassConstForXi.GetPz()) && // check rapidity of XicPlus
                     AliVertexingHFUtils::CheckKFParticleCov(kfpXicPlus_wMassConstForLam_woMassConstForXi) && // check covariance matrix
                     err_massXicPlus_Rec_wMassConstForLam_woMassConstForXi>0 && // err_massXi > 0
                     kfpXicPlus_wMassConstForLam_woMassConstForXi.GetChi2()/kfpXicPlus_wMassConstForLam_woMassConstForXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                     kfpXicPlus_wMassConstForLam_woMassConstForXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                   ) {
                  Int_t lab_XicPlus = -9999;
                  if (fIsMC) lab_XicPlus = MatchToMCXicPlus(ptrack, ntrack, btrack, trackP[itrkBP_trk1], trackP[itrkBP_trk2], mcArray);
                  FillQATreeXicPlusFromCasc_wMassConstForLam_woMassConstForXi(kfpLambda_wMassConst, kfpXiMinus_wMassConstForLam_woMassConstForXi, kfpXicPlus_wMassConstForLam_woMassConstForXi, PV, trackPiFromXicPlus_HighPt, mcArray, lab_XicPlus, AODEvent);
                }
                kfpXicPlus_wMassConstForLam_woMassConstForXi.Clear();
                kfpBP_HighPt.Clear();
                kfpBP_LowPt.Clear();
              }
            }
          }

          ////// (With mass constraint for Xi) + (with mass constraint for Lambda)
          // Mass Constraint for Xi
          KFParticle kfpXiMinus_wMassConstForLamAndXi = kfpXiMinus_wMassConstForLam_woMassConstForXi;
          kfpXiMinus_wMassConstForLamAndXi.SetNonlinearMassConstraint(massXi);
          if ( AliVertexingHFUtils::CheckKFParticleCov(kfpXiMinus_wMassConstForLamAndXi) && TMath::Abs(kfpXiMinus_wMassConstForLamAndXi.GetE()) > TMath::Abs(kfpXiMinus_wMassConstForLamAndXi.GetPz()) ) {
            // Lambda with mass and topo constraint to Xi
            KFParticle kfpLambda_wMassConst_To_Xi_wMassConst = kfpLambda_wMassConst;
            kfpLambda_wMassConst_To_Xi_wMassConst.SetProductionVertex(kfpXiMinus_wMassConstForLamAndXi);
            // Reconstruct Xi with Lambda (both mass and topo constraint)
            KFParticle kfpXiMinus_wMassAndTopoConstForLam_woMassConstForXi;
            const KFParticle *vXiDs_wMassAndTopoConstForLam_woMassConstForXi[2] = {&kfpPion_ForXi, &kfpLambda_wMassConst_To_Xi_wMassConst};
            kfpXiMinus_wMassAndTopoConstForLam_woMassConstForXi.Construct(vXiDs_wMassAndTopoConstForLam_woMassConstForXi, 2);

            // Xi with mass constraint
            KFParticle kfpXiMinus_wMassAndTopoConstForLam_wMassConstForXi = kfpXiMinus_wMassAndTopoConstForLam_woMassConstForXi;
            kfpXiMinus_wMassAndTopoConstForLam_wMassConstForXi.SetNonlinearMassConstraint(massXi);

            // Xi with mass and topo constraint
            KFParticle kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi = kfpXiMinus_wMassAndTopoConstForLam_wMassConstForXi;
            kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.SetProductionVertex(PV);

            // Loop for bachelor pions
            for (Int_t itrkBP_trk1=0; itrkBP_trk1<(flag_trkP-1); itrkBP_trk1++) { // Loop for first bachelor pion+
              for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkP; itrkBP_trk2++) { // Loop for second bachelor pion+
                if ( trackP[itrkBP_trk1]->GetID() == ptrack->GetID() || trackP[itrkBP_trk2]->GetID() == ptrack->GetID() || trackP[itrkBP_trk1]->GetID() == trackP[itrkBP_trk2]->GetID() ) continue;

                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk1]) ) continue;
                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk2]) ) continue;

                // === pt(pion_0)<pt(pion_1) ===
                AliAODTrack *trackPiFromXicPlus_LowPt = NULL, *trackPiFromXicPlus_HighPt = NULL;
                trackPiFromXicPlus_LowPt  = trackP[itrkBP_trk1];
                trackPiFromXicPlus_HighPt = trackP[itrkBP_trk2];
                if (trackP[itrkBP_trk1]->Pt() > trackP[itrkBP_trk2]->Pt()) {
                  trackPiFromXicPlus_LowPt  = trackP[itrkBP_trk2];
                  trackPiFromXicPlus_HighPt = trackP[itrkBP_trk1];
                }
                // =============================
                KFParticle kfpBP_LowPt  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicPlus_LowPt, 211);
                KFParticle kfpBP_HighPt = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicPlus_HighPt, 211);

                // reconstruct Xic+
                KFParticle kfpXicPlus_wMassConstForLamAndXi;
                const KFParticle *vXicPlusDs_wMassConstForLamAndXi[3] = {&kfpXiMinus_wMassConstForLamAndXi, &kfpBP_LowPt, &kfpBP_HighPt};
                kfpXicPlus_wMassConstForLamAndXi.Construct(vXicPlusDs_wMassConstForLamAndXi, 3);

                Float_t massXicPlus_Rec_wMassConstForLamAndXi, err_massXicPlus_Rec_wMassConstForLamAndXi;
                kfpXicPlus_wMassConstForLamAndXi.GetMass(massXicPlus_Rec_wMassConstForLamAndXi, err_massXicPlus_Rec_wMassConstForLamAndXi);

                // reconstruct Xic+ (wMassAndTopoConstForLam_wMassConstForXi)
                KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi;
                const KFParticle *vXicPlusDs_wMassAndTopoConstForLam_wMassConstForXi[3] = {&kfpXiMinus_wMassAndTopoConstForLam_wMassConstForXi, &kfpBP_LowPt, &kfpBP_HighPt};
                kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.Construct(vXicPlusDs_wMassAndTopoConstForLam_wMassConstForXi, 3);

                Float_t massXicPlus_Rec_wMassAndTopoConstForLam_wMassConstForXi, err_massXicPlus_Rec_wMassAndTopoConstForLam_wMassConstForXi;
                kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetMass(massXicPlus_Rec_wMassAndTopoConstForLam_wMassConstForXi, err_massXicPlus_Rec_wMassAndTopoConstForLam_wMassConstForXi);

                if ( kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetNDF()>0 && kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetChi2()>0 && // chi2>0 && NDF>0
                     TMath::Abs(kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetE())>TMath::Abs(kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetPz()) && // check rapidity of XicPlus
                     AliVertexingHFUtils::CheckKFParticleCov(kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi) && // check covariance matrix
                     err_massXicPlus_Rec_wMassAndTopoConstForLam_wMassConstForXi>0 && // err_massXicPlus > 0
                     kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetChi2()/kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                     kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                   ) {
                  Int_t lab_XicPlus = -9999;
                  if (fIsMC) lab_XicPlus = MatchToMCXicPlus(ptrack, ntrack, btrack, trackP[itrkBP_trk1], trackP[itrkBP_trk2], mcArray);
                  FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassConstForXi(kfpLambda_wMassConst_To_Xi_wMassConst, kfpXiMinus_wMassAndTopoConstForLam_wMassConstForXi, kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi, PV, trackPiFromXicPlus_HighPt, mcArray, lab_XicPlus, AODEvent);
                }
                kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.Clear();

                // reconstruct Xic+ (wMassAndTopoConstForLam_wMassAndTopoConstForXi)
                KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi;
                const KFParticle *vXicPlusDs_wMassAndTopoConstForLam_wMassAndTopoConstForXi[3] = {&kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi, &kfpBP_LowPt, &kfpBP_HighPt};
                kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.Construct(vXicPlusDs_wMassAndTopoConstForLam_wMassAndTopoConstForXi, 3);

                Float_t massXicPlus_Rec_wMassAndTopoConstForLam_wMassAndTopoConstForXi, err_massXicPlus_Rec_wMassAndTopoConstForLam_wMassAndTopoConstForXi;
                kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetMass(massXicPlus_Rec_wMassAndTopoConstForLam_wMassAndTopoConstForXi, err_massXicPlus_Rec_wMassAndTopoConstForLam_wMassAndTopoConstForXi);

                if ( kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetNDF()>0 && kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetChi2()>0 && // chi2>0 && NDF>0
                     TMath::Abs(kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetE())>TMath::Abs(kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetPz()) && // check rapidity of XicPlus
                     AliVertexingHFUtils::CheckKFParticleCov(kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi) && // check covariance matrix
                     err_massXicPlus_Rec_wMassAndTopoConstForLam_wMassAndTopoConstForXi>0 && // err_massXicPlus > 0
                     kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetChi2()/kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                     kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                   ) {
                  Int_t lab_XicPlus = -9999;
                  if (fIsMC) lab_XicPlus = MatchToMCXicPlus(ptrack, ntrack, btrack, trackP[itrkBP_trk1], trackP[itrkBP_trk2], mcArray);
                  FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassAndTopoConstForXi(kfpLambda_wMassConst_To_Xi_wMassConst, kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi, kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi, PV, trackPiFromXicPlus_HighPt, mcArray, lab_XicPlus, AODEvent);

                  // reconstruct Xic+ (wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic)
                  KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi;
                  kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.SetProductionVertex(PV);
                  Float_t massXicPlus_Rec_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic, err_massXicPlus_Rec_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic;
                  kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetMass(massXicPlus_Rec_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic, err_massXicPlus_Rec_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic);

                  if ( kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetNDF()>0 && kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetChi2()>0 && // chi2>0 && NDF>0
                       TMath::Abs(kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetE())>TMath::Abs(kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetPz()) && // check rapidity of XicPlus
                       AliVertexingHFUtils::CheckKFParticleCov(kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic) && // check covariance matrix
                       err_massXicPlus_Rec_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic>0 && // err_massXicPlus > 0
                       kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetChi2()/kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                       kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                     ) {
                    FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic(kfpLambda_wMassConst_To_Xi_wMassConst, kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi, kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic, PV, trackPiFromXicPlus_HighPt, trackPiFromXicPlus_LowPt, btrack, ptrack, ntrack, kfpBP_HighPt, kfpBP_LowPt, kfpPion_ForXi, kfpProton, kfpPionMinus, mcArray, lab_XicPlus, AODEvent);
                  }
                  kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.Clear();
                }
                kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.Clear();

                if ( kfpXicPlus_wMassConstForLamAndXi.GetNDF()>0 && kfpXicPlus_wMassConstForLamAndXi.GetChi2()>0 && // chi2>0 && NDF>0
                     TMath::Abs(kfpXicPlus_wMassConstForLamAndXi.GetE())>TMath::Abs(kfpXicPlus_wMassConstForLamAndXi.GetPz()) && // check rapidity of XicPlus
                     AliVertexingHFUtils::CheckKFParticleCov(kfpXicPlus_wMassConstForLamAndXi) && // check covariance matrix
                     err_massXicPlus_Rec_wMassConstForLamAndXi>0 && // err_massXicPlus > 0
                     kfpXicPlus_wMassConstForLamAndXi.GetChi2()/kfpXicPlus_wMassConstForLamAndXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                     kfpXicPlus_wMassConstForLamAndXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                   ) {
                  Int_t lab_XicPlus = -9999;
                  if (fIsMC) lab_XicPlus = MatchToMCXicPlus(ptrack, ntrack, btrack, trackP[itrkBP_trk1], trackP[itrkBP_trk2], mcArray);
                  FillQATreeXicPlusFromCasc_wMassConstForLamAndXi(kfpLambda_wMassConst, kfpXiMinus_wMassConstForLamAndXi, kfpXicPlus_wMassConstForLamAndXi, PV, trackPiFromXicPlus_HighPt, mcArray, lab_XicPlus, AODEvent);
                }
                kfpXicPlus_wMassConstForLamAndXi.Clear();
                kfpBP_HighPt.Clear();
                kfpBP_LowPt.Clear();
              }
            }
            kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.Clear();
            kfpXiMinus_wMassAndTopoConstForLam_wMassConstForXi.Clear();
            kfpXiMinus_wMassAndTopoConstForLam_woMassConstForXi.Clear();
            kfpLambda_wMassConst_To_Xi_wMassConst.Clear();
          }
          kfpXiMinus_wMassConstForLamAndXi.Clear();
          kfpXiMinus_wMassConstForLam_woMassConstForXi.Clear();
        }
        kfpLambda_wMassConst.Clear();
        kfpPion_ForXi.Clear();

      }
      // ================================ (2022.09.07) ================================

      KFParticle kfpLambda_m = kfpLambda;
      kfpLambda_m.SetNonlinearMassConstraint(massLambda);

      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpLambda_m) || TMath::Abs(kfpLambda_m.GetE()) <= TMath::Abs(kfpLambda_m.GetPz()) ) continue;

      KFParticle kfpPionOrKaon;
      kfpPionOrKaon = AliVertexingHFUtils::CreateKFParticleFromAODtrack(btrack, -211); // pion-
      KFParticle kfpXiMinus;
      const KFParticle *vXiDs[2] = {&kfpPionOrKaon, &kfpLambda_m};
      kfpXiMinus.Construct(vXiDs, 2);

      // check rapidity of Xi-
      if ( TMath::Abs(kfpXiMinus.GetE())<=TMath::Abs(kfpXiMinus.GetPz()) ) continue;

      // err_massXi > 0
      Float_t massXiMinus_Rec, err_massXiMinus_Rec;
      kfpXiMinus.GetMass(massXiMinus_Rec, err_massXiMinus_Rec);
      if ( err_massXiMinus_Rec<=0 ) continue;

      // chi2>0 && NDF>0
      if ( kfpXiMinus.GetNDF()<=0 || kfpXiMinus.GetChi2()<=0 ) continue;

      // Prefilter
      if ( kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF() >= fAnaCuts->GetKFPXi_Chi2geoMax() ) continue;

      // check covariance matrix
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXiMinus) ) continue;

      // mass window cut of Xi-
      if ( TMath::Abs(massXiMinus_Rec-massXi) > (fAnaCuts->GetProdMassTolXi()) ) continue;

      KFParticle kfpXiMinus_m = kfpXiMinus;
      kfpXiMinus_m.SetNonlinearMassConstraint(massXi);

      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXiMinus_m) || TMath::Abs(kfpXiMinus_m.GetE()) <= TMath::Abs(kfpXiMinus_m.GetPz()) ) continue;

      for (Int_t itrkBP_trk1=0; itrkBP_trk1<(flag_trkP-1); itrkBP_trk1++) { // Loop for first bachelor pion+
        for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkP; itrkBP_trk2++) { // Loop for second bachelor pion+

          if ( trackP[itrkBP_trk1]->GetID() == ptrack->GetID() || trackP[itrkBP_trk2]->GetID() == ptrack->GetID() || trackP[itrkBP_trk1]->GetID() == trackP[itrkBP_trk2]->GetID() ) continue;

          if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk1]) ) continue;
          if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP_trk2]) ) continue;

          // === pt(pion_0)<pt(pion_1) ===
          AliAODTrack *trackPiFromXicPlus_LowPt = NULL, *trackPiFromXicPlus_HighPt = NULL;
          trackPiFromXicPlus_LowPt  = trackP[itrkBP_trk1];
          trackPiFromXicPlus_HighPt = trackP[itrkBP_trk2];

          if (trackP[itrkBP_trk1]->Pt() > trackP[itrkBP_trk2]->Pt()) {
            trackPiFromXicPlus_LowPt  = trackP[itrkBP_trk2];
            trackPiFromXicPlus_HighPt = trackP[itrkBP_trk1];
          }
          // =============================
          
          // Check if the track is used for primary vertex fit
//          if (trackPiFromXicPlus_LowPt->GetUsedForPrimVtxFit()) cout << "!!!!! Pi0 is used for PV fit !!!!!" << endl;
//          if (trackPiFromXicPlus_HighPt->GetUsedForPrimVtxFit()) cout << "!!!!! Pi1 is used for PV fit !!!!!" << endl;

          KFParticle kfpBP_LowPt  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicPlus_LowPt, 211);
          KFParticle kfpBP_HighPt = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicPlus_HighPt, 211);

          // reconstruct XicPlus
          KFParticle kfpXicPlus;
          const KFParticle *vXicPlusDs[3] = {&kfpXiMinus_m, &kfpBP_LowPt, &kfpBP_HighPt};
          kfpXicPlus.Construct(vXicPlusDs, 3);

          // chi2>0 && NDF>0
          if ( kfpXicPlus.GetNDF()<=0 || kfpXicPlus.GetChi2()<=0 ) continue;

          // Prefilter
          if ( kfpXicPlus.GetChi2()/kfpXicPlus.GetNDF() >= fAnaCuts->GetKFPXicPlus_Chi2geoMax() ) continue;
          if ( kfpXicPlus.GetPt() < fAnaCuts->GetPtMinXicPlus() ) continue;

          // check rapidity of XicPlus
          if ( TMath::Abs(kfpXicPlus.GetE())<=TMath::Abs(kfpXicPlus.GetPz()) ) continue;

          // check covariance matrix
          if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXicPlus) ) continue;

          // err_massXicPlus > 0
          Float_t massXicPlus_Rec, err_massXicPlus_Rec;
          kfpXicPlus.GetMass(massXicPlus_Rec, err_massXicPlus_Rec);
          if ( err_massXicPlus_Rec<=0 ) continue;

          if (fWriteXicPlusTree) {
            fHCountUsedForPrimVtxFit->Fill(5);
            if (btrack->GetUsedForPrimVtxFit()) {
              // cout << "!!!!! bachelor is used for PV fit !!!!!" << endl;
              fHCountUsedForPrimVtxFit->Fill(2);}
            if (ptrack->GetUsedForPrimVtxFit()) {
              //cout << "!!!!! V0 positive is used for PV fit !!!!!" << endl;
              fHCountUsedForPrimVtxFit->Fill(4);}
            if (ntrack->GetUsedForPrimVtxFit()) {
              //cout << "!!!!! V0 negative is used for PV fit !!!!!" << endl;
              fHCountUsedForPrimVtxFit->Fill(3);}
            if (trackPiFromXicPlus_LowPt->GetUsedForPrimVtxFit()) {
              //cout << "!!!!! Pi0 is used for PV fit !!!!!" << endl;
              fHCountUsedForPrimVtxFit->Fill(0);}
            if (trackPiFromXicPlus_HighPt->GetUsedForPrimVtxFit()) {
              //cout << "!!!!! Pi1 is used for PV fit !!!!!" << endl;
              fHCountUsedForPrimVtxFit->Fill(1);}

            Int_t lab_XicPlus = -9999;
            if (fIsMC) {
              lab_XicPlus = MatchToMCXicPlus(ptrack, ntrack, btrack, trackP[itrkBP_trk1], trackP[itrkBP_trk2], mcArray);
              //if (lab_XicPlus>=0) FillTreeRecXicPlusFromCasc(AODEvent, casc, kfpXicPlus, trackPiFromXicPlus_LowPt, kfpBP_LowPt, kfpXiMinus, kfpXiMinus_m, kfpPionOrKaon, btrack, kfpK0Short, kfpGamma, kfpLambda, kfpLambda_m, ptrack, ntrack, trackPiFromXicPlus_HighPt, kfpBP_HighPt, kfpProton, kfpPionMinus, PV, mcArray, lab_XicPlus);
              FillTreeRecXicPlusFromCasc(AODEvent, casc, kfpXicPlus, trackPiFromXicPlus_LowPt, kfpBP_LowPt, kfpXiMinus, kfpXiMinus_m, kfpPionOrKaon, btrack, kfpK0Short, kfpGamma, kfpLambda, kfpLambda_m, ptrack, ntrack, trackPiFromXicPlus_HighPt, kfpBP_HighPt, kfpProton, kfpPionMinus, PV, PV_KF_Refit, mcArray, lab_XicPlus);
            }
            if (!fIsMC) {
             FillTreeRecXicPlusFromCasc(AODEvent, casc, kfpXicPlus, trackPiFromXicPlus_LowPt, kfpBP_LowPt, kfpXiMinus, kfpXiMinus_m, kfpPionOrKaon, btrack, kfpK0Short, kfpGamma, kfpLambda, kfpLambda_m, ptrack, ntrack, trackPiFromXicPlus_HighPt, kfpBP_HighPt, kfpProton, kfpPionMinus, PV, PV_KF_Refit, mcArray, lab_XicPlus);
            }
          }
          kfpXicPlus.Clear();
          kfpBP_HighPt.Clear();
          kfpBP_LowPt.Clear();
        } // Loop for second bachelor pion+
      } // Loop for first bachelor pion+
      kfpXiMinus_m.Clear();
      kfpXiMinus.Clear();
      kfpPionOrKaon.Clear();
      kfpLambda_m.Clear();
      kfpLambda.Clear();
    }

    if ( btrack->Charge()>0 ) { // Xi^+

      const KFParticle *vAntiDaughters[2] = {&kfpPionPlus, &kfpAntiProton};

      KFParticle kfpAntiLambda;
      kfpAntiLambda.Construct(vAntiDaughters, 2);
      Float_t massAntiLambda_Rec, err_massAntiLambda_Rec;
      kfpAntiLambda.GetMass(massAntiLambda_Rec, err_massAntiLambda_Rec);

      // check rapidity of Anti-Lambda
      if ( TMath::Abs(kfpAntiLambda.GetE())<=TMath::Abs(kfpAntiLambda.GetPz()) ) continue;

      // chi2>0 && NDF>0 for selecting Anti-Lambda
      if ( kfpAntiLambda.GetNDF()<=0 || kfpAntiLambda.GetChi2()<=0 ) continue;

      // check cov. of Anti-Lambda
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLambda) ) continue;

      // err_mass>0 of Anti-Lambda
      if ( err_massAntiLambda_Rec<=0 ) continue;

      // Chi2geo cut of Anti-Lambda
      if ( (kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF()) >= fAnaCuts->GetKFPLam_Chi2geoMax() ) continue;

      //************************** calculate l/Δl for Anti-Lambda *************************************
      Double_t dx_AntiLambda = PV.GetX()-kfpAntiLambda.GetX();
      Double_t dy_AntiLambda = PV.GetY()-kfpAntiLambda.GetY();
      Double_t dz_AntiLambda = PV.GetZ()-kfpAntiLambda.GetZ();
      Double_t l_AntiLambda = TMath::Sqrt(dx_AntiLambda*dx_AntiLambda + dy_AntiLambda*dy_AntiLambda + dz_AntiLambda*dz_AntiLambda);
      Double_t dl_AntiLambda = (PV.GetCovariance(0)+kfpAntiLambda.GetCovariance(0))*dx_AntiLambda*dx_AntiLambda + (PV.GetCovariance(2)+kfpAntiLambda.GetCovariance(2))*dy_AntiLambda*dy_AntiLambda + (PV.GetCovariance(5)+kfpAntiLambda.GetCovariance(5))*dz_AntiLambda*dz_AntiLambda + 2*( (PV.GetCovariance(1)+kfpAntiLambda.GetCovariance(1))*dx_AntiLambda*dy_AntiLambda + (PV.GetCovariance(3)+kfpAntiLambda.GetCovariance(3))*dx_AntiLambda*dz_AntiLambda + (PV.GetCovariance(4)+kfpAntiLambda.GetCovariance(4))*dy_AntiLambda*dz_AntiLambda );
      if ( fabs(l_AntiLambda)<1.e-8f ) l_AntiLambda = 1.e-8f;
      dl_AntiLambda = dl_AntiLambda<0. ? 1.e8f : sqrt(dl_AntiLambda)/l_AntiLambda;
      Double_t nErr_l_AntiLambda = l_AntiLambda/dl_AntiLambda;
      //***************************************************************************************

      // l/Deltal cut of Anti-Lambda
      if ( nErr_l_AntiLambda <= fAnaCuts->GetKFPLam_lDeltalMin() ) continue;

      // mass window cut of Anti-Lambda
      if ( TMath::Abs(massAntiLambda_Rec-massLambda) > (fAnaCuts->GetProdMassTolLambda()) ) continue;

      // ================================ QA study (2022.09.07) ================================
      if (fWriteXicPlusQATree) {
        KFParticle kfpPion_ForXiPlus;
        kfpPion_ForXiPlus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(btrack, 211); // pion+

        // reconstruct Xi+
        KFParticle kfpXiPlus_woMassConstForLamAndXi;
        const KFParticle *vXiDs_woMassConstForLamAndXi[2] = {&kfpPion_ForXiPlus, &kfpAntiLambda};
        kfpXiPlus_woMassConstForLamAndXi.Construct(vXiDs_woMassConstForLamAndXi, 2);
        Float_t massXiPlus_Rec_woMassConstForLamAndXi, err_massXiPlus_Rec_woMassConstForLamAndXi;
        kfpXiPlus_woMassConstForLamAndXi.GetMass(massXiPlus_Rec_woMassConstForLamAndXi, err_massXiPlus_Rec_woMassConstForLamAndXi);



        ////// Without mass constraint for Anti-Lambda
        if ( TMath::Abs(kfpXiPlus_woMassConstForLamAndXi.GetE())>TMath::Abs(kfpXiPlus_woMassConstForLamAndXi.GetPz()) && // check rapidity of Xi+
             err_massXiPlus_Rec_woMassConstForLamAndXi>0 && // err_massXi > 0
             kfpXiPlus_woMassConstForLamAndXi.GetNDF()>0 && kfpXiPlus_woMassConstForLamAndXi.GetChi2()>0 && // chi2>0 && NDF>0
             AliVertexingHFUtils::CheckKFParticleCov(kfpXiPlus_woMassConstForLamAndXi) && // check covariance matrix
             kfpXiPlus_woMassConstForLamAndXi.GetChi2()/kfpXiPlus_woMassConstForLamAndXi.GetNDF() < fAnaCuts->GetKFPXi_Chi2geoMax() && // Prefilter
             TMath::Abs(massXiPlus_Rec_woMassConstForLamAndXi-massXi) <= (fAnaCuts->GetProdMassTolXi()) // mass window cut of Xi+
           ) {
          ////// Without mass constraint for Xi
          // Loop for bachelor pions
          for (Int_t itrkBP_trk1=0; itrkBP_trk1<(flag_trkN-1); itrkBP_trk1++) { // Loop for first bachelor pion-
            for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkN; itrkBP_trk2++) { // Loop for second bachelor pion-
              if ( trackN[itrkBP_trk1]->GetID() == ntrack->GetID() || trackN[itrkBP_trk2]->GetID() == ntrack->GetID() || trackN[itrkBP_trk1]->GetID() == trackN[itrkBP_trk2]->GetID() ) continue;
              if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk1]) ) continue;
              if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk2]) ) continue;

              // === pt(pion_0)<pt(pion_1) ===
              AliAODTrack *trackPiFromXicMinus_LowPt = NULL, *trackPiFromXicMinus_HighPt = NULL;
              trackPiFromXicMinus_LowPt  = trackN[itrkBP_trk1];
              trackPiFromXicMinus_HighPt = trackN[itrkBP_trk2];
              if (trackN[itrkBP_trk1]->Pt() > trackN[itrkBP_trk2]->Pt()) {
                trackPiFromXicMinus_LowPt  = trackN[itrkBP_trk2];
                trackPiFromXicMinus_HighPt = trackN[itrkBP_trk1];
              }
              // =============================
              
              KFParticle kfpBP_LowPt  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicMinus_LowPt, -211);
              KFParticle kfpBP_HighPt = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicMinus_HighPt, -211);

              // reconstruct Xic-
              KFParticle kfpXicMinus_woMassConstForLamAndXi;
              const KFParticle *vXicMinusDs_woMassConstForLamAndXi[3] = {&kfpXiPlus_woMassConstForLamAndXi, &kfpBP_LowPt, &kfpBP_HighPt};
              kfpXicMinus_woMassConstForLamAndXi.Construct(vXicMinusDs_woMassConstForLamAndXi, 3);

              Float_t massXicMinus_Rec_woMassConstForLamAndXi, err_massXicMinus_Rec_woMassConstForLamAndXi;
              kfpXicMinus_woMassConstForLamAndXi.GetMass(massXicMinus_Rec_woMassConstForLamAndXi, err_massXicMinus_Rec_woMassConstForLamAndXi);

              if ( kfpXicMinus_woMassConstForLamAndXi.GetNDF()>0 && kfpXicMinus_woMassConstForLamAndXi.GetChi2()>0 && // chi2>0 && NDF>0
                   TMath::Abs(kfpXicMinus_woMassConstForLamAndXi.GetE())>TMath::Abs(kfpXicMinus_woMassConstForLamAndXi.GetPz()) && // check rapidity of Xic-
                   AliVertexingHFUtils::CheckKFParticleCov(kfpXicMinus_woMassConstForLamAndXi) && // check covariance matrix
                   err_massXicMinus_Rec_woMassConstForLamAndXi>0 && // err_massXicMinus > 0
                   kfpXicMinus_woMassConstForLamAndXi.GetChi2()/kfpXicMinus_woMassConstForLamAndXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                   kfpXicMinus_woMassConstForLamAndXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                 ) {
                Int_t lab_XicMinus = -9999;
                if (fIsMC) lab_XicMinus = MatchToMCAntiXicPlus(ntrack, ptrack, btrack, trackN[itrkBP_trk1], trackN[itrkBP_trk2], mcArray);
                FillQATreeXicPlusFromCasc_woMassConstForLamAndXi(kfpAntiLambda, kfpXiPlus_woMassConstForLamAndXi, kfpXicMinus_woMassConstForLamAndXi, PV, trackPiFromXicMinus_HighPt, mcArray, lab_XicMinus, AODEvent);
              }
              kfpXicMinus_woMassConstForLamAndXi.Clear();
              kfpBP_HighPt.Clear();
              kfpBP_LowPt.Clear();
            }
          }

          ////// (With mass Constraint for Xi) + (without mass constraint for Anti-Lambda)
          KFParticle kfpXiPlus_woMassConstForLam_wMassConstForXi = kfpXiPlus_woMassConstForLamAndXi;
          kfpXiPlus_woMassConstForLam_wMassConstForXi.SetNonlinearMassConstraint(massXi);
          if ( AliVertexingHFUtils::CheckKFParticleCov(kfpXiPlus_woMassConstForLam_wMassConstForXi) && TMath::Abs(kfpXiPlus_woMassConstForLam_wMassConstForXi.GetE()) > TMath::Abs(kfpXiPlus_woMassConstForLam_wMassConstForXi.GetPz()) ) {
            // Loop for bachelor pions
            for (Int_t itrkBP_trk1=0; itrkBP_trk1<(flag_trkN-1); itrkBP_trk1++) { // Loop for first bachelor pion-
              for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkN; itrkBP_trk2++) { // Loop for second bachelor pion-
                if ( trackN[itrkBP_trk1]->GetID() == ntrack->GetID() || trackN[itrkBP_trk2]->GetID() == ntrack->GetID() || trackN[itrkBP_trk1]->GetID() == trackN[itrkBP_trk2]->GetID() ) continue;

                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk1]) ) continue;
                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk2]) ) continue;

                // === pt(pion_0)<pt(pion_1) ===
                AliAODTrack *trackPiFromXicMinus_LowPt = NULL, *trackPiFromXicMinus_HighPt = NULL;
                trackPiFromXicMinus_LowPt  = trackN[itrkBP_trk1];
                trackPiFromXicMinus_HighPt = trackN[itrkBP_trk2];
                if (trackN[itrkBP_trk1]->Pt() > trackN[itrkBP_trk2]->Pt()) {
                  trackPiFromXicMinus_LowPt  = trackN[itrkBP_trk2];
                  trackPiFromXicMinus_HighPt = trackN[itrkBP_trk1];
                }
                // =============================
                KFParticle kfpBP_LowPt  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicMinus_LowPt, -211);
                KFParticle kfpBP_HighPt = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicMinus_HighPt, -211);

                // reconstruct Xic-
                KFParticle kfpXicMinus_woMassConstForLam_wMassConstForXi;
                const KFParticle *vXicMinusDs_woMassConstForLam_wMassConstForXi[3] = {&kfpXiPlus_woMassConstForLam_wMassConstForXi, &kfpBP_LowPt, &kfpBP_HighPt};
                kfpXicMinus_woMassConstForLam_wMassConstForXi.Construct(vXicMinusDs_woMassConstForLam_wMassConstForXi, 3);

                Float_t massXicMinus_Rec_woMassConstForLam_wMassConstForXi, err_massXicMinus_Rec_woMassConstForLam_wMassConstForXi;
                kfpXicMinus_woMassConstForLam_wMassConstForXi.GetMass(massXicMinus_Rec_woMassConstForLam_wMassConstForXi, err_massXicMinus_Rec_woMassConstForLam_wMassConstForXi);

                if ( kfpXicMinus_woMassConstForLam_wMassConstForXi.GetNDF()>0 && kfpXicMinus_woMassConstForLam_wMassConstForXi.GetChi2()>0 && // chi2>0 && NDF>0
                     TMath::Abs(kfpXicMinus_woMassConstForLam_wMassConstForXi.GetE())>TMath::Abs(kfpXicMinus_woMassConstForLam_wMassConstForXi.GetPz()) && // check rapidity of Xic-
                     AliVertexingHFUtils::CheckKFParticleCov(kfpXicMinus_woMassConstForLam_wMassConstForXi) && // check covariance matrix
                     err_massXicMinus_Rec_woMassConstForLam_wMassConstForXi>0 && // err_massXicMinus > 0
                     kfpXicMinus_woMassConstForLam_wMassConstForXi.GetChi2()/kfpXicMinus_woMassConstForLam_wMassConstForXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                     kfpXicMinus_woMassConstForLam_wMassConstForXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                   ) {
                  Int_t lab_XicMinus = -9999;
                  if (fIsMC) lab_XicMinus = MatchToMCAntiXicPlus(ntrack, ptrack, btrack, trackN[itrkBP_trk1], trackN[itrkBP_trk2], mcArray);
                  FillQATreeXicPlusFromCasc_woMassConstForLam_wMassConstForXi(kfpAntiLambda, kfpXiPlus_woMassConstForLam_wMassConstForXi, kfpXicMinus_woMassConstForLam_wMassConstForXi, PV, trackPiFromXicMinus_HighPt, mcArray, lab_XicMinus, AODEvent);
                }
                kfpXicMinus_woMassConstForLam_wMassConstForXi.Clear();
                kfpBP_HighPt.Clear();
                kfpBP_LowPt.Clear();
              }
            }
          }
          kfpXiPlus_woMassConstForLam_wMassConstForXi.Clear();
          kfpXiPlus_woMassConstForLamAndXi.Clear();
        }



        ////// With mass constraint for Anti-Lambda
        // Mass Constraint for Anti-Lambda
        KFParticle kfpAntiLambda_wMassConst = kfpAntiLambda;
        kfpAntiLambda_wMassConst.SetNonlinearMassConstraint(massLambda);
        if ( AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLambda_wMassConst) && TMath::Abs(kfpAntiLambda_wMassConst.GetE()) > TMath::Abs(kfpAntiLambda_wMassConst.GetPz()) ) {
          KFParticle kfpXiPlus_wMassConstForLam_woMassConstForXi;
          const KFParticle *vXiDs_wMassConstForAntiLam_woMassConstForXi[2] = {&kfpPion_ForXiPlus, &kfpAntiLambda_wMassConst};
          kfpXiPlus_wMassConstForLam_woMassConstForXi.Construct(vXiDs_wMassConstForAntiLam_woMassConstForXi, 2);

          Float_t massXiPlus_Rec_wMassConstForLam_woMassConstForXi, err_massXiPlus_Rec_wMassConstForLam_woMassConstForXi;
          kfpXiPlus_wMassConstForLam_woMassConstForXi.GetMass(massXiPlus_Rec_wMassConstForLam_woMassConstForXi, err_massXiPlus_Rec_wMassConstForLam_woMassConstForXi);

          if ( TMath::Abs(kfpXiPlus_wMassConstForLam_woMassConstForXi.GetE())>TMath::Abs(kfpXiPlus_wMassConstForLam_woMassConstForXi.GetPz()) && // check rapidity of Xi+
               err_massXiPlus_Rec_wMassConstForLam_woMassConstForXi>0 && // err_massXi > 0
               kfpXiPlus_wMassConstForLam_woMassConstForXi.GetNDF()>0 && kfpXiPlus_wMassConstForLam_woMassConstForXi.GetChi2()>0 && // chi2>0 && NDF>0
               AliVertexingHFUtils::CheckKFParticleCov(kfpXiPlus_wMassConstForLam_woMassConstForXi) && // check covariance matrix
               kfpXiPlus_wMassConstForLam_woMassConstForXi.GetChi2()/kfpXiPlus_wMassConstForLam_woMassConstForXi.GetNDF() < fAnaCuts->GetKFPXi_Chi2geoMax() && // Prefilter
               TMath::Abs(massXiPlus_Rec_wMassConstForLam_woMassConstForXi-massXi) <= (fAnaCuts->GetProdMassTolXi()) // mass window cut of Xi-
             ) {
            // Loop for bachelor pions
            for (Int_t itrkBP_trk1=0; itrkBP_trk1<(flag_trkN-1); itrkBP_trk1++) { // Loop for first bachelor pion-
              for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkN; itrkBP_trk2++) { // Loop for second bachelor pion-
                if ( trackN[itrkBP_trk1]->GetID() == ntrack->GetID() || trackN[itrkBP_trk2]->GetID() == ntrack->GetID() || trackN[itrkBP_trk1]->GetID() == trackN[itrkBP_trk2]->GetID() ) continue;
                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk1]) ) continue;
                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk2]) ) continue;

                // === pt(pion_0)<pt(pion_1) ===
                AliAODTrack *trackPiFromXicMinus_LowPt = NULL, *trackPiFromXicMinus_HighPt = NULL;
                trackPiFromXicMinus_LowPt  = trackN[itrkBP_trk1];
                trackPiFromXicMinus_HighPt = trackN[itrkBP_trk2];

                if (trackN[itrkBP_trk1]->Pt() > trackN[itrkBP_trk2]->Pt()) {
                  trackPiFromXicMinus_LowPt  = trackN[itrkBP_trk2];
                  trackPiFromXicMinus_HighPt = trackN[itrkBP_trk1];
                }
                // =============================
                KFParticle kfpBP_LowPt  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicMinus_LowPt, -211);
                KFParticle kfpBP_HighPt = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicMinus_HighPt, -211);

                // reconstruct Xic-
                KFParticle kfpXicMinus_wMassConstForLam_woMassConstForXi;
                const KFParticle *vXicMinusDs_wMassConstForLam_woMassConstForXi[3] = {&kfpXiPlus_wMassConstForLam_woMassConstForXi, &kfpBP_LowPt, &kfpBP_HighPt};
                kfpXicMinus_wMassConstForLam_woMassConstForXi.Construct(vXicMinusDs_wMassConstForLam_woMassConstForXi, 3);

                Float_t massXicMinus_Rec_wMassConstForLam_woMassConstForXi, err_massXicMinus_Rec_wMassConstForLam_woMassConstForXi;
                kfpXicMinus_wMassConstForLam_woMassConstForXi.GetMass(massXicMinus_Rec_wMassConstForLam_woMassConstForXi, err_massXicMinus_Rec_wMassConstForLam_woMassConstForXi);

                if ( kfpXicMinus_wMassConstForLam_woMassConstForXi.GetNDF()>0 && kfpXicMinus_wMassConstForLam_woMassConstForXi.GetChi2()>0 && // chi2>0 && NDF>0
                     TMath::Abs(kfpXicMinus_wMassConstForLam_woMassConstForXi.GetE())>TMath::Abs(kfpXicMinus_wMassConstForLam_woMassConstForXi.GetPz()) && // check rapidity of Xic-
                     AliVertexingHFUtils::CheckKFParticleCov(kfpXicMinus_wMassConstForLam_woMassConstForXi) && // check covariance matrix
                     err_massXicMinus_Rec_wMassConstForLam_woMassConstForXi>0 && // err_massXicMinus > 0
                     kfpXicMinus_wMassConstForLam_woMassConstForXi.GetChi2()/kfpXicMinus_wMassConstForLam_woMassConstForXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                     kfpXicMinus_wMassConstForLam_woMassConstForXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                   ) {
                  Int_t lab_XicMinus = -9999;
                  if (fIsMC) lab_XicMinus = MatchToMCAntiXicPlus(ntrack, ptrack, btrack, trackN[itrkBP_trk1], trackN[itrkBP_trk2], mcArray);
                  FillQATreeXicPlusFromCasc_wMassConstForLam_woMassConstForXi(kfpAntiLambda_wMassConst, kfpXiPlus_wMassConstForLam_woMassConstForXi, kfpXicMinus_wMassConstForLam_woMassConstForXi, PV, trackPiFromXicMinus_HighPt, mcArray, lab_XicMinus, AODEvent);
                }
                kfpXicMinus_wMassConstForLam_woMassConstForXi.Clear();
                kfpBP_HighPt.Clear();
                kfpBP_LowPt.Clear();
              }
            }
          }

          ////// (With mass constraint for Xi) + (with mass constraint for Anti-Lambda)
          // Mass Constraint for Xi
          KFParticle kfpXiPlus_wMassConstForLamAndXi = kfpXiPlus_wMassConstForLam_woMassConstForXi;
          kfpXiPlus_wMassConstForLamAndXi.SetNonlinearMassConstraint(massXi);
          if ( AliVertexingHFUtils::CheckKFParticleCov(kfpXiPlus_wMassConstForLamAndXi) && TMath::Abs(kfpXiPlus_wMassConstForLamAndXi.GetE()) > TMath::Abs(kfpXiPlus_wMassConstForLamAndXi.GetPz()) ) {
            // Anti-Lambda with mass and topo constraint to Xi
            KFParticle kfpAntiLambda_wMassConst_To_XiPlus_wMassConst = kfpAntiLambda_wMassConst;
            kfpAntiLambda_wMassConst_To_XiPlus_wMassConst.SetProductionVertex(kfpXiPlus_wMassConstForLamAndXi);
            // Reconstruct Xi with Anti-Lambda (both mass and topo constraint)
            KFParticle kfpXiPlus_wMassAndTopoConstForAntiLam_woMassConstForXi;
            const KFParticle *vXiDs_wMassAndTopoConstForAntiLam_woMassConstForXi[2] = {&kfpPion_ForXiPlus, &kfpAntiLambda_wMassConst_To_XiPlus_wMassConst};
            kfpXiPlus_wMassAndTopoConstForAntiLam_woMassConstForXi.Construct(vXiDs_wMassAndTopoConstForAntiLam_woMassConstForXi, 2);

            // XicPlus with mass constraint
            KFParticle kfpXiPlus_wMassAndTopoConstForAntiLam_wMassConstForXi = kfpXiPlus_wMassAndTopoConstForAntiLam_woMassConstForXi;
            kfpXiPlus_wMassAndTopoConstForAntiLam_wMassConstForXi.SetNonlinearMassConstraint(massXi);

            // Xi with mass and topo constraint
            KFParticle kfpXiPlus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi = kfpXiPlus_wMassAndTopoConstForAntiLam_wMassConstForXi;
            kfpXiPlus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.SetProductionVertex(PV);

            // Loop for bachelor pions
            for (Int_t itrkBP_trk1=0; itrkBP_trk1<(flag_trkN-1); itrkBP_trk1++) { // Loop for first bachelor pion-
              for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkN; itrkBP_trk2++) { // Loop for second bachelor pion-
                if ( trackN[itrkBP_trk1]->GetID() == ntrack->GetID() || trackN[itrkBP_trk2]->GetID() == ntrack->GetID() || trackN[itrkBP_trk1]->GetID() == trackN[itrkBP_trk2]->GetID() ) continue;

                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk1]) ) continue;
                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk2]) ) continue;

                // === pt(pion_0)<pt(pion_1) ===
                AliAODTrack *trackPiFromXicMinus_LowPt = NULL, *trackPiFromXicMinus_HighPt = NULL;
                trackPiFromXicMinus_LowPt  = trackN[itrkBP_trk1];
                trackPiFromXicMinus_HighPt = trackN[itrkBP_trk2];
                if (trackN[itrkBP_trk1]->Pt() > trackN[itrkBP_trk2]->Pt()) {
                  trackPiFromXicMinus_LowPt  = trackN[itrkBP_trk2];
                  trackPiFromXicMinus_HighPt = trackN[itrkBP_trk1];
                }
                // =============================
                KFParticle kfpBP_LowPt  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicMinus_LowPt, -211);
                KFParticle kfpBP_HighPt = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXicMinus_HighPt, -211);

                // reconstruct Xic-
                KFParticle kfpXicMinus_wMassConstForLamAndXi;
                const KFParticle *vXicMinusDs_wMassConstForLamAndXi[3] = {&kfpXiPlus_wMassConstForLamAndXi, &kfpBP_LowPt, &kfpBP_HighPt};
                kfpXicMinus_wMassConstForLamAndXi.Construct(vXicMinusDs_wMassConstForLamAndXi, 3);

                Float_t massXicMinus_Rec_wMassConstForLamAndXi, err_massXicMinus_Rec_wMassConstForLamAndXi;
                kfpXicMinus_wMassConstForLamAndXi.GetMass(massXicMinus_Rec_wMassConstForLamAndXi, err_massXicMinus_Rec_wMassConstForLamAndXi);

                // reconstruct Xic- (wMassAndTopoConstForAntiLam_wMassConstForXi)
                KFParticle kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi;
                const KFParticle *vXicMinusDs_wMassAndTopoConstForAntiLam_wMassConstForXi[3] = {&kfpXiPlus_wMassAndTopoConstForAntiLam_wMassConstForXi, &kfpBP_LowPt, &kfpBP_HighPt};
                kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi.Construct(vXicMinusDs_wMassAndTopoConstForAntiLam_wMassConstForXi, 3);

                Float_t massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassConstForXi, err_massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassConstForXi;
                kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi.GetMass(massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassConstForXi, err_massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassConstForXi);

                if ( kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi.GetNDF()>0 && kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi.GetChi2()>0 && // chi2>0 && NDF>0
                     TMath::Abs(kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi.GetE())>TMath::Abs(kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi.GetPz()) && // check rapidity of XicPlus
                     AliVertexingHFUtils::CheckKFParticleCov(kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi) && // check covariance matrix
                     err_massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassConstForXi>0 && // err_massXicPlus > 0
                     kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi.GetChi2()/kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                     kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                   ) {
                  Int_t lab_XicMinus = -9999;
                  if (fIsMC) lab_XicMinus = MatchToMCAntiXicPlus(ntrack, ptrack, btrack, trackN[itrkBP_trk1], trackN[itrkBP_trk2], mcArray);
                  FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassConstForXi(kfpAntiLambda_wMassConst_To_XiPlus_wMassConst, kfpXiPlus_wMassAndTopoConstForAntiLam_wMassConstForXi, kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi, PV, trackPiFromXicMinus_HighPt, mcArray, lab_XicMinus, AODEvent);
                }
                kfpXicMinus_wMassAndTopoConstForAntiLam_wMassConstForXi.Clear();

                // reconstruct Xic- (wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi)
                KFParticle kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi;
                const KFParticle *vXicMinusDs_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi[3] = {&kfpXiPlus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi, &kfpBP_LowPt, &kfpBP_HighPt};
                kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.Construct(vXicMinusDs_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi, 3);

                Float_t massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi, err_massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi;
                kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.GetMass(massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi, err_massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi);

                if ( kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.GetNDF()>0 && kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.GetChi2()>0 && // chi2>0 && NDF>0
                     TMath::Abs(kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.GetE())>TMath::Abs(kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.GetPz()) && // check rapidity of XicPlus
                     AliVertexingHFUtils::CheckKFParticleCov(kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi) && // check covariance matrix
                     err_massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi>0 && // err_massXicPlus > 0
                     kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.GetChi2()/kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                     kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                   ) {
                  Int_t lab_XicMinus = -9999;
                  if (fIsMC) lab_XicMinus = MatchToMCAntiXicPlus(ntrack, ptrack, btrack, trackN[itrkBP_trk1], trackN[itrkBP_trk2], mcArray);
                  FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassAndTopoConstForXi(kfpAntiLambda_wMassConst_To_XiPlus_wMassConst, kfpXiPlus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi, kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi, PV, trackPiFromXicMinus_HighPt, mcArray, lab_XicMinus, AODEvent);

                  // reconstruct Xic- (wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic)
                  KFParticle kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic = kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi;
                  kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic.SetProductionVertex(PV);
                  Float_t massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic, err_massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic;
                  kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic.GetMass(massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic, err_massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic);

                  if ( kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic.GetNDF()>0 && kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic.GetChi2()>0 && // chi2>0 && NDF>0
                       TMath::Abs(kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic.GetE())>TMath::Abs(kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic.GetPz()) && // check rapidity of XicPlus
                       AliVertexingHFUtils::CheckKFParticleCov(kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic) && // check covariance matrix
                       err_massXicMinus_Rec_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic>0 && // err_massXicPlus > 0
                       kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic.GetChi2()/kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                       kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                     ) {
                    FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic(kfpAntiLambda_wMassConst_To_XiPlus_wMassConst, kfpXiPlus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi, kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic, PV, trackPiFromXicMinus_HighPt, trackPiFromXicMinus_LowPt, btrack, ntrack, ptrack, kfpBP_HighPt, kfpBP_LowPt, kfpPion_ForXiPlus, kfpAntiProton, kfpPionPlus, mcArray, lab_XicMinus, AODEvent);
                  }
                  kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi_wTopoConstForXic.Clear();
                }
                kfpXicMinus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.Clear();

                if ( kfpXicMinus_wMassConstForLamAndXi.GetNDF()>0 && kfpXicMinus_wMassConstForLamAndXi.GetChi2()>0 && // chi2>0 && NDF>0
                     TMath::Abs(kfpXicMinus_wMassConstForLamAndXi.GetE())>TMath::Abs(kfpXicMinus_wMassConstForLamAndXi.GetPz()) && // check rapidity of Xic-
                     AliVertexingHFUtils::CheckKFParticleCov(kfpXicMinus_wMassConstForLamAndXi) && // check covariance matrix
                     err_massXicMinus_Rec_wMassConstForLamAndXi>0 && // err_massXicMinus > 0
                     kfpXicMinus_wMassConstForLamAndXi.GetChi2()/kfpXicMinus_wMassConstForLamAndXi.GetNDF() < fAnaCuts->GetKFPXicPlus_Chi2geoMax() && // Prefilter
                     kfpXicMinus_wMassConstForLamAndXi.GetPt() >= fAnaCuts->GetPtMinXicPlus() // Prefilter
                   ) {
                  Int_t lab_XicMinus = -9999;
                  if (fIsMC) lab_XicMinus = MatchToMCAntiXicPlus(ntrack, ptrack, btrack, trackN[itrkBP_trk1], trackN[itrkBP_trk2], mcArray);
                  FillQATreeXicPlusFromCasc_wMassConstForLamAndXi(kfpAntiLambda_wMassConst, kfpXiPlus_wMassConstForLamAndXi, kfpXicMinus_wMassConstForLamAndXi, PV, trackPiFromXicMinus_HighPt, mcArray, lab_XicMinus, AODEvent);
                }
                kfpXicMinus_wMassConstForLamAndXi.Clear();
                kfpBP_HighPt.Clear();
                kfpBP_LowPt.Clear();
              }
            }
            kfpXiPlus_wMassAndTopoConstForAntiLam_wMassAndTopoConstForXi.Clear();
            kfpXiPlus_wMassAndTopoConstForAntiLam_wMassConstForXi.Clear();
            kfpXiPlus_wMassAndTopoConstForAntiLam_woMassConstForXi.Clear();
            kfpAntiLambda_wMassConst_To_XiPlus_wMassConst.Clear();
          }
          kfpXiPlus_wMassConstForLamAndXi.Clear();
          kfpXiPlus_wMassConstForLam_woMassConstForXi.Clear();
        }
        kfpAntiLambda_wMassConst.Clear();
        kfpPion_ForXiPlus.Clear();

      }
      // ================================ (2022.09.07) ================================

      KFParticle kfpAntiLambda_m = kfpAntiLambda;
      kfpAntiLambda_m.SetNonlinearMassConstraint(massLambda);

      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiLambda_m) || TMath::Abs(kfpAntiLambda_m.GetE()) <= TMath::Abs(kfpAntiLambda_m.GetPz()) ) continue;

      KFParticle kfpPionOrKaon;
      kfpPionOrKaon = AliVertexingHFUtils::CreateKFParticleFromAODtrack(btrack, 211); // pion+
      KFParticle kfpXiPlus;
      const KFParticle *vXiDs[2] = {&kfpPionOrKaon, &kfpAntiLambda_m};
      kfpXiPlus.Construct(vXiDs, 2);

      // check rapidity of Xi+
      if ( TMath::Abs(kfpXiPlus.GetE())<=TMath::Abs(kfpXiPlus.GetPz()) ) continue;

      // err_massXi > 0
      Float_t massXiPlus_Rec, err_massXiPlus_Rec;
      kfpXiPlus.GetMass(massXiPlus_Rec, err_massXiPlus_Rec);
      if ( err_massXiPlus_Rec<=0 ) continue;

      // chi2>0 && NDF>0
      if ( kfpXiPlus.GetNDF()<=0 || kfpXiPlus.GetChi2()<=0 ) continue;

      // Prefilter
      if ( kfpXiPlus.GetChi2()/kfpXiPlus.GetNDF() >= fAnaCuts->GetKFPXi_Chi2geoMax() ) continue;

      // check covariance matrix
      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXiPlus) ) continue;

      // mass window cut of Xi+
      if ( TMath::Abs(massXiPlus_Rec-massXi) > (fAnaCuts->GetProdMassTolXi()) ) continue;

      KFParticle kfpXiPlus_m = kfpXiPlus;
      kfpXiPlus_m.SetNonlinearMassConstraint(massXi);

      if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpXiPlus_m) || TMath::Abs(kfpXiPlus_m.GetE()) <= TMath::Abs(kfpXiPlus_m.GetPz()) ) continue;

      for (Int_t itrkBP_trk1=0; itrkBP_trk1<(flag_trkN-1); itrkBP_trk1++) { // Loop for first bachelor pion-
        for (Int_t itrkBP_trk2=itrkBP_trk1+1; itrkBP_trk2<flag_trkN; itrkBP_trk2++) { // Loop for second bachelor pion-

          if ( trackN[itrkBP_trk1]->GetID() == ntrack->GetID() || trackN[itrkBP_trk2]->GetID() == ntrack->GetID() || trackN[itrkBP_trk1]->GetID() == trackN[itrkBP_trk2]->GetID() ) continue;

          if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk1]) ) continue;
          if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP_trk2]) ) continue;

          // === pt(pion_0)<pt(pion_1) ===
          AliAODTrack *trackPiFromAntiXicPlus_LowPt = NULL, *trackPiFromAntiXicPlus_HighPt = NULL;
          trackPiFromAntiXicPlus_LowPt  = trackN[itrkBP_trk1];
          trackPiFromAntiXicPlus_HighPt = trackN[itrkBP_trk2];

          if (trackN[itrkBP_trk1]->Pt() > trackN[itrkBP_trk2]->Pt()) {
            trackPiFromAntiXicPlus_LowPt  = trackN[itrkBP_trk2];
            trackPiFromAntiXicPlus_HighPt = trackN[itrkBP_trk1];
          }
          // =============================

          KFParticle kfpBP_LowPt  = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromAntiXicPlus_LowPt, -211);
          KFParticle kfpBP_HighPt = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromAntiXicPlus_HighPt, -211);

          // reconstruct Anti-XicPlus
          KFParticle kfpAntiXicPlus;
          const KFParticle *vXicPlusDs[3] = {&kfpXiPlus_m, &kfpBP_LowPt, &kfpBP_HighPt};
          kfpAntiXicPlus.Construct(vXicPlusDs, 3);

          // chi2>0 && NDF>0
          if ( kfpAntiXicPlus.GetNDF()<=0 || kfpAntiXicPlus.GetChi2()<=0 ) continue;

          // Prefilter
          if ( kfpAntiXicPlus.GetChi2()/kfpAntiXicPlus.GetNDF() >= fAnaCuts->GetKFPXicPlus_Chi2geoMax() ) continue;
          if ( kfpAntiXicPlus.GetPt() < fAnaCuts->GetPtMinXicPlus() ) continue;

          // check rapidity of Anti-XicPlus
          if ( TMath::Abs(kfpAntiXicPlus.GetE())<=TMath::Abs(kfpAntiXicPlus.GetPz()) ) continue;

          // check covariance matrix
          if ( !AliVertexingHFUtils::CheckKFParticleCov(kfpAntiXicPlus) ) continue;

          // err_massAntiXicPlus > 0
          Float_t massAntiXicPlus_Rec, err_massAntiXicPlus_Rec;
          kfpAntiXicPlus.GetMass(massAntiXicPlus_Rec, err_massAntiXicPlus_Rec);
          if ( err_massAntiXicPlus_Rec<=0 ) continue;

          if (fWriteXicPlusTree) {
            fHCountUsedForPrimVtxFit->Fill(5);
            if (btrack->GetUsedForPrimVtxFit()) {
              //cout << "!!!!! bachelor is used for PV fit !!!!!" << endl;
              fHCountUsedForPrimVtxFit->Fill(2);}
            if (ptrack->GetUsedForPrimVtxFit()) {
              //cout << "!!!!! V0 positive is used for PV fit !!!!!" << endl;
              fHCountUsedForPrimVtxFit->Fill(3);}
            if (ntrack->GetUsedForPrimVtxFit()) {
              //cout << "!!!!! V0 negative is used for PV fit !!!!!" << endl;
              fHCountUsedForPrimVtxFit->Fill(4);}
            if (trackPiFromAntiXicPlus_LowPt->GetUsedForPrimVtxFit()) {
              //cout << "!!!!! Pi0 is used for PV fit !!!!!" << endl;
              fHCountUsedForPrimVtxFit->Fill(0);}
            if (trackPiFromAntiXicPlus_HighPt->GetUsedForPrimVtxFit()) {
              //cout << "!!!!! Pi1 is used for PV fit !!!!!" << endl;
              fHCountUsedForPrimVtxFit->Fill(1);}

            Int_t lab_AntiXicPlus = -9999.;
            if (fIsMC) {
              lab_AntiXicPlus = MatchToMCAntiXicPlus(ntrack, ptrack, btrack, trackN[itrkBP_trk1], trackN[itrkBP_trk2], mcArray);
              //if (lab_AntiXicPlus>=0) FillTreeRecXicPlusFromCasc(AODEvent, casc, kfpAntiXicPlus, trackPiFromAntiXicPlus_LowPt, kfpBP_LowPt, kfpXiPlus, kfpXiPlus_m, kfpPionOrKaon, btrack, kfpK0Short, kfpGamma, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, trackPiFromAntiXicPlus_HighPt, kfpBP_HighPt, kfpAntiProton, kfpPionPlus, PV, mcArray, lab_AntiXicPlus);
              FillTreeRecXicPlusFromCasc(AODEvent, casc, kfpAntiXicPlus, trackPiFromAntiXicPlus_LowPt, kfpBP_LowPt, kfpXiPlus, kfpXiPlus_m, kfpPionOrKaon, btrack, kfpK0Short, kfpGamma, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, trackPiFromAntiXicPlus_HighPt, kfpBP_HighPt, kfpAntiProton, kfpPionPlus, PV, PV_KF_Refit, mcArray, lab_AntiXicPlus);
            }
            if (!fIsMC) {
              FillTreeRecXicPlusFromCasc(AODEvent, casc, kfpAntiXicPlus, trackPiFromAntiXicPlus_LowPt, kfpBP_LowPt, kfpXiPlus, kfpXiPlus_m, kfpPionOrKaon, btrack, kfpK0Short, kfpGamma, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, trackPiFromAntiXicPlus_HighPt, kfpBP_HighPt, kfpAntiProton, kfpPionPlus, PV, PV_KF_Refit, mcArray, lab_AntiXicPlus);
            }
          }
          kfpAntiXicPlus.Clear();
          kfpBP_HighPt.Clear();
          kfpBP_LowPt.Clear();
        } // Loop for second bachelor pion-
      } // Loop for first bachelor pion-
      kfpXiPlus_m.Clear();
      kfpXiPlus.Clear();
      kfpPionOrKaon.Clear();
      kfpAntiLambda_m.Clear();
      kfpAntiLambda.Clear();
    }

    kfpK0Short.Clear();
    kfpPionPlus.Clear();
    kfpAntiProton.Clear();
    kfpPionMinus.Clear();
    kfpProton.Clear();

  }

  return;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCXicPlus(AliAODTrack *trackProton, AliAODTrack *trackPionMinus_2, AliAODTrack *trackPionMinus_1, AliAODTrack *trackPionPlus_trk1, AliAODTrack *trackPionPlus_trk2, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle
  
  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<0) return -99;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));

  Int_t labelPionMinus_2  = fabs(trackPionMinus_2->GetLabel());
  if (labelPionMinus_2<0) return -99;
  AliAODMCParticle* mcPionMinus_2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus_2));

  Int_t labelPionMinus_1  = fabs(trackPionMinus_1->GetLabel());
  if (labelPionMinus_1<0) return -99;
  AliAODMCParticle* mcPionMinus_1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus_1));

  Int_t labelPionPlus_trk1  = fabs(trackPionPlus_trk1->GetLabel());
  if (labelPionPlus_trk1<0) return -99;
  AliAODMCParticle* mcPionPlus_trk1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus_trk1));

  Int_t labelPionPlus_trk2  = fabs(trackPionPlus_trk2->GetLabel());
  if (labelPionPlus_trk2<0) return -99;
  AliAODMCParticle* mcPionPlus_trk2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus_trk2));

  if ( mcProton->GetPdgCode() != 2212 || mcPionMinus_2->GetPdgCode() != -211 || mcPionMinus_1->GetPdgCode() != -211 || mcPionPlus_trk1->GetPdgCode() != 211 || mcPionPlus_trk2->GetPdgCode() != 211 ) return -90; // check pdg

  // === check Lambda ===
  Int_t IndexMother[2];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPionMinus_2->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -80; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -70; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 || mcMother->GetNDaughters()!=2 ) return -60; // check mother is lambda
  // ====================

  // === check Xi- ===
  IndexMother[0] = mcMother->GetMother(); // mother of lambda
  IndexMother[1] = mcPionMinus_1->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -50; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -40; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3312 || mcMother->GetNDaughters()!=2 ) return -30; // check mother is Xi-
  // =================

  // === check XicPlus ===
  Int_t IndexMother_XicPlusDau[3];
  IndexMother_XicPlusDau[0] = mcMother->GetMother(); // mother of Xi-
  IndexMother_XicPlusDau[1] = mcPionPlus_trk1->GetMother();
  IndexMother_XicPlusDau[2] = mcPionPlus_trk2->GetMother();
  if ( IndexMother_XicPlusDau[0]<0 || IndexMother_XicPlusDau[1]<0 || IndexMother_XicPlusDau[2]<0 ) return -20; // check mother exist
  if ( (IndexMother_XicPlusDau[0] != IndexMother_XicPlusDau[1]) || (IndexMother_XicPlusDau[0] != IndexMother_XicPlusDau[2]) ) { // come from different mother
    if ( (IndexMother_XicPlusDau[0] == IndexMother_XicPlusDau[1]) || (IndexMother_XicPlusDau[0] == IndexMother_XicPlusDau[2]) ) { // Xi and Pi from same mother
      AliAODMCParticle *mcMother_Xi  = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother_XicPlusDau[0]));
      AliAODMCParticle *mcMother_Pi0 = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother_XicPlusDau[1]));
      AliAODMCParticle *mcMother_Pi1 = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother_XicPlusDau[2]));
      if ( mcMother_Xi->GetPdgCode() == 3324 && ( (mcMother_Pi0->GetPdgCode()==4232) || (mcMother_Pi1->GetPdgCode()==4232) ) ) {
        Int_t Index_XiStar0_Mother = mcMother_Xi->GetMother();
        if (Index_XiStar0_Mother<0) return -15;
        AliAODMCParticle *mcMother_XiStar0 = static_cast<AliAODMCParticle*>(mcArray->At(Index_XiStar0_Mother));
        if ( mcMother_XiStar0->GetPdgCode()==4232 && mcMother_XiStar0->GetNDaughters()==2 && (Index_XiStar0_Mother==IndexMother_XicPlusDau[0] || Index_XiStar0_Mother==IndexMother_XicPlusDau[1]) ) {
          if (AliVertexingHFUtils::CheckOrigin(mcArray,mcMother_XiStar0,kTRUE)==4) return -4;
          if (AliVertexingHFUtils::CheckOrigin(mcArray,mcMother_XiStar0,kTRUE)==5) return -5;
        }
      }
    }
    return -10; // check the same mother
  }
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother_XicPlusDau[0]));
  if ( mcMother->GetPdgCode() != 4232 || mcMother->GetNDaughters()!=3 ) {
    return -8; // check mother is XicPlus
  }
  // ==================

//  if ( mcMother->IsPrimary() ) return 1;
//  if ( mcMother->IsPhysicalPrimary() ) return 2;
//  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
//  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
//  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

  Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcMother,kTRUE);
  return CheckOrigin;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCAntiXicPlus(AliAODTrack *trackAntiProton, AliAODTrack *trackPionPlus_2, AliAODTrack *trackPionPlus_1, AliAODTrack *trackPionMinus_trk1, AliAODTrack *trackPionMinus_trk2, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle
  
  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<0) return -99;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));

  Int_t labelPionPlus_2  = fabs(trackPionPlus_2->GetLabel());
  if (labelPionPlus_2<0) return -99;
  AliAODMCParticle* mcPionPlus_2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus_2));

  Int_t labelPionPlus_1  = fabs(trackPionPlus_1->GetLabel());
  if (labelPionPlus_1<0) return -99;
  AliAODMCParticle* mcPionPlus_1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionPlus_1));

  Int_t labelPionMinus_trk1  = fabs(trackPionMinus_trk1->GetLabel());
  if (labelPionMinus_trk1<0) return -99;
  AliAODMCParticle* mcPionMinus_trk1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus_trk1));

  Int_t labelPionMinus_trk2  = fabs(trackPionMinus_trk2->GetLabel());
  if (labelPionMinus_trk2<0) return -99;
  AliAODMCParticle* mcPionMinus_trk2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPionMinus_trk2));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcPionPlus_2->GetPdgCode() != 211 || mcPionPlus_1->GetPdgCode() != 211 || mcPionMinus_trk1->GetPdgCode() != -211 || mcPionMinus_trk2->GetPdgCode() != -211 ) return -90; // check pdg

  // === check Anti-Lambda ===
  Int_t IndexMother[2];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcPionPlus_2->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -80; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -70; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 || mcMother->GetNDaughters()!=2 ) return -60; // check mother is Anti-Lambda

  // === check Xi+ ===
  IndexMother[0] = mcMother->GetMother(); // mother of Anti-Lambda
  IndexMother[1] = mcPionPlus_1->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -50; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -40; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3312 || mcMother->GetNDaughters()!=2 ) return -30; // check mother is Xi+
  // =================

  // === check Anti-XicPlus ===
  Int_t IndexMother_XicPlusDau[3];
  IndexMother_XicPlusDau[0] = mcMother->GetMother(); // mother of Xi+
  IndexMother_XicPlusDau[1] = mcPionMinus_trk1->GetMother();
  IndexMother_XicPlusDau[2] = mcPionMinus_trk2->GetMother();
  // =======================
  if ( IndexMother_XicPlusDau[0]<0 || IndexMother_XicPlusDau[1]<0 || IndexMother_XicPlusDau[2]<0 ) return -20; // check mother exist
  if ( (IndexMother_XicPlusDau[0] != IndexMother_XicPlusDau[1]) || (IndexMother_XicPlusDau[0] != IndexMother_XicPlusDau[2]) ) { // come from different mother
    if ( (IndexMother_XicPlusDau[0] == IndexMother_XicPlusDau[1]) || (IndexMother_XicPlusDau[0] == IndexMother_XicPlusDau[2]) ) { // Xi and Pi from same mother
      AliAODMCParticle *mcMother_Xi  = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother_XicPlusDau[0]));
      AliAODMCParticle *mcMother_Pi0 = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother_XicPlusDau[1]));
      AliAODMCParticle *mcMother_Pi1 = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother_XicPlusDau[2]));
      if ( mcMother_Xi->GetPdgCode() == -3324 && ( (mcMother_Pi0->GetPdgCode()==-4232) || (mcMother_Pi1->GetPdgCode()==-4232) ) ) {
        Int_t Index_XiStar0_Mother = mcMother_Xi->GetMother();
        if (Index_XiStar0_Mother<0) return -15;
        AliAODMCParticle *mcMother_XiStar0 = static_cast<AliAODMCParticle*>(mcArray->At(Index_XiStar0_Mother));
        if ( mcMother_XiStar0->GetPdgCode()==-4232 && mcMother_XiStar0->GetNDaughters()==2 && (Index_XiStar0_Mother==IndexMother_XicPlusDau[0] || Index_XiStar0_Mother==IndexMother_XicPlusDau[1]) ) {
          if (AliVertexingHFUtils::CheckOrigin(mcArray,mcMother_Xi,kTRUE)==4) return -4;
          if (AliVertexingHFUtils::CheckOrigin(mcArray,mcMother_Xi,kTRUE)==5) return -5;
        }
      }
    }
    return -10; // check the same mother
  }
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother_XicPlusDau[0]));
  if ( mcMother->GetPdgCode() != -4232 || mcMother->GetNDaughters()!=3 ) {
    return -8; // check mother is Anti-XicPlus
  }

//  if ( mcMother->IsPrimary() ) return 1;
//  if ( mcMother->IsPhysicalPrimary() ) return 2;
//  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
//  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
//  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

  Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcMother,kTRUE);
  return CheckOrigin;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::SelectTrack(AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks, Bool_t *seleFlags)
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
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCXiMinus(AliAODTrack *trackProton, AliAODTrack *trackPion3, AliAODTrack *trackPion2, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
  Int_t labelPion3  = fabs(trackPion3->GetLabel());
  if (labelPion3<0) return -1;
  AliAODMCParticle* mcPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion3));
  Int_t labelPion2  = fabs(trackPion2->GetLabel());
  if (labelPion2<0) return -1;
  AliAODMCParticle* mcPion2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion2));

  if ( mcProton->GetPdgCode() != 2212 || mcPion3->GetPdgCode() != -211 || mcPion2->GetPdgCode() != -211 ) return -1; // check pdg

  Int_t IndexMother[4];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPion3->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is lambda

  IndexMother[2] = mcMother->GetMother(); // mother of lambda
  IndexMother[3] = mcPion2->GetMother();
  if ( IndexMother[2]<0 || IndexMother[3]<0 ) return -1; // check mother exist
  if ( IndexMother[2] != IndexMother[3] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[2]));
  if ( mcMother->GetPdgCode() != 3312 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xi-

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[2];
  return 1;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCXiPlus(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, AliAODTrack *trackAntiPion2, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
  Int_t labelAntiPion3  = fabs(trackAntiPion3->GetLabel());
  if (labelAntiPion3<0) return -1;
  AliAODMCParticle* mcAntiPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion3));
  Int_t labelAntiPion2  = fabs(trackAntiPion2->GetLabel());
  if (labelAntiPion2<0) return -1;
  AliAODMCParticle* mcAntiPion2 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion2));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcAntiPion3->GetPdgCode() != 211 || mcAntiPion2->GetPdgCode() != 211 ) return -1; // check pdg

  Int_t IndexMother[4];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcAntiPion3->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Anti-lambda

  IndexMother[2] = mcMother->GetMother(); // mother of lambda
  IndexMother[3] = mcAntiPion2->GetMother();
  if ( IndexMother[2]<0 || IndexMother[3]<0 ) return -1; // check mother exist
  if ( IndexMother[2] != IndexMother[3] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[2]));
  if ( mcMother->GetPdgCode() != -3312 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xi+

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[2];
  return 1;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCLambda(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<=0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
  Int_t labelPion3  = fabs(trackPion3->GetLabel());
  if (labelPion3<=0) return -1;
  AliAODMCParticle* mcPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion3));

  if ( mcProton->GetPdgCode() != 2212 || mcPion3->GetPdgCode() != -211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is lambda and only have two daughters

//  AliAODMCParticle* mcMother2 = static_cast<AliAODMCParticle*>(mcArray->At(fabs(mcMother->GetMother())));

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[0];
  return 1;

}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCLambdaFromXi(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<=0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
  Int_t labelPion3  = fabs(trackPion3->GetLabel());
  if (labelPion3<=0) return -1;
  AliAODMCParticle* mcPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion3));

  if ( mcProton->GetPdgCode() != 2212 || mcPion3->GetPdgCode() != -211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 ) return -1; // check mother is lambda
  
  Int_t Index_Mother_Xi = mcMother->GetMother();
  if ( Index_Mother_Xi<=0 ) return -1;
  AliAODMCParticle* mcMother_Xi = static_cast<AliAODMCParticle*>(mcArray->At(Index_Mother_Xi));
  if ( mcMother_Xi->GetPdgCode() != 3312 ) return -1;

//  if ( !mcMother_Xi->IsPhysicalPrimary() ) return -1; // check IsPhysicalPrimary()

  return IndexMother[0];

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCAntiLambda(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<=0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
  Int_t labelAntiPion3  = fabs(trackAntiPion3->GetLabel());
  if (labelAntiPion3<=0) return -1;
  AliAODMCParticle* mcAntiPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion3));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcAntiPion3->GetPdgCode() != 211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcAntiPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Anti-lambda and only have two daughters

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[0];
  return 1;

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCAntiLambdaFromXi(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<=0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
  Int_t labelAntiPion3  = fabs(trackAntiPion3->GetLabel());
  if (labelAntiPion3<=0) return -1;
  AliAODMCParticle* mcAntiPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion3));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcAntiPion3->GetPdgCode() != 211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcAntiPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 ) return -1; // check mother is Anti-lambda

  Int_t Index_Mother_Xi = mcMother->GetMother();
  if ( Index_Mother_Xi<=0 ) return -1;
  AliAODMCParticle* mcMother_Xi = static_cast<AliAODMCParticle*>(mcArray->At(Index_Mother_Xi));
  if ( mcMother_Xi->GetPdgCode() != -3312 ) return -1;

//  if ( !mcMother_Xi->IsPhysicalPrimary() ) return -1; // check IsPhysicalPrimary()

  return IndexMother[0];

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToMCPion(AliAODTrack *track, TClonesArray *mcArray)
{
  Int_t labelPion = fabs(track->GetLabel());
  if (labelPion<=0) return -1;
  AliAODMCParticle* mcPion = static_cast<AliAODMCParticle*>(mcArray->At(labelPion));
  if ( TMath::Abs(mcPion->GetPdgCode()) != 211 ) return -1;

  return labelPion;
}

/*
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::MatchToXicPlusMC(TClonesArray *mcArray, Int_t PDGXicPlus, const Int_t nDaughters, const Int_t *daughterIndex, const Int_t *daughterPDG)
{
  ///
  /// Check if this candidate is matched to a MC signal XicPlus
  /// If yes, return Index (>=0) of the AliAODMCParticle
  /// If no, return -1
  ///

  Int_t IndexMom[10] = {0};
  Int_t IndexMother=-1;
  Bool_t pdgUsed[10] = {0};

  // loop on daughter Index
  for(Int_t i=0; i<nDaughters; i++) {
    IndexMom[i]=-1;
    Int_t Index = daughterIndex[i];
    if(Index<0) {
      printf("daughter with negative index %d\n", Index);
      return -1;
    }
    AliAODMCParticle *part = (AliAODMCParticle*)mcArray->At(Index);
    if(!part) { 
      printf("no MC particle\n");
      return -1;
    }

    // check the PDG of the daughter
    Int_t pdgPart = part->GetPdgCode();
    for(Int_t j=0; j<nDaughters; j++) {
      if(!pdgUsed[j] && pdgPart==daughterPDG[j]) {
        pdgUsed[j]=kTRUE;
        break;
      }
    }

    AliAODMCParticle *mother = part;
    while ( mother->GetMother()>=0 ) {
      IndexMother = mother->GetMother();
      mother = (AliAODMCParticle*)mcArray->At(IndexMother);
      if (!mother) {
        printf("no MC mother particle\n");
        break;
      }
      Int_t pdgMother = mother->GetPdgCode();
      if ( pdgMother==PDGXicPlus ) { // check mother is XicPlus
        IndexMom[i]=IndexMother;
        break;
      } else if( pdgMother>PDGXicPlus || pdgMother<10 ) {
        break;
      }
    }

    if( IndexMom[i]==-1 ) return -1; // mother PDG not ok for this daughter

  } // end loop on daughters

  IndexMother=IndexMom[0];
  for(Int_t i=0; i<nDaughters; i++) {
    // all Index have to be the same and !=-1
    if(IndexMom[i]==-1)          return -1;
    if(IndexMom[i]!=IndexMother) return -1;
    // check that all daughter PDGs are matched
    if(pdgUsed[i]==kFALSE)       return -1;
  }

  return IndexMother;
}
*/

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineEvent()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(3)->GetContainer()->GetName();
  fTree_Event = new TTree(nameoutput, "Event");
  Int_t nVar = 12;
  fVar_Event = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  fVarNames[0]  = "centrality";
  fVarNames[1]  = "x_vtx_reco";
  fVarNames[2]  = "y_vtx_reco";
  fVarNames[3]  = "z_vtx_reco";
  fVarNames[4]  = "err_x_vtx_reco";
  fVarNames[5]  = "err_y_vtx_reco";
  fVarNames[6]  = "err_z_vtx_reco";
  fVarNames[7]  = "n_vtx_contributors";
  fVarNames[8]  = "n_tracks";
  fVarNames[9]  = "is_ev_rej";
  fVarNames[10] = "run_number";
  fVarNames[11] = "ev_id";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_Event->Branch(fVarNames[ivar].Data(), &fVar_Event[ivar], Form("%s/F", fVarNames[ivar].Data()));
  }

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeQAXicPlus()
{
  const char* nameoutput = GetOutputSlot(7)->GetContainer()->GetName();
  fTree_XicPlus_QA = new TTree(nameoutput, "XicPlus variables QA tree");
  Int_t nVar = 62;
  fVar_XicPlus_QA = new Float_t[nVar-1];
  TString *fVarNames_QA = new TString[nVar];

  fVarNames_QA[0]  = "PV_X_MC"; // MC
  fVarNames_QA[1]  = "PV_Y_MC"; // MC
  fVarNames_QA[2]  = "PV_Z_MC"; // MC
  fVarNames_QA[3]  = "PV_X_rec";
  fVarNames_QA[4]  = "PV_Y_rec";
  fVarNames_QA[5]  = "PV_Z_rec";
  fVarNames_QA[6]  = "PV_sigma_X_rec";
  fVarNames_QA[7]  = "PV_sigma_Y_rec";
  fVarNames_QA[8]  = "PV_sigma_Z_rec";
  fVarNames_QA[9]  = "recalPV_X_rec";
  fVarNames_QA[10] = "recalPV_Y_rec";
  fVarNames_QA[11] = "recalPV_Z_rec";
  fVarNames_QA[12] = "recalPV_sigma_X_rec";
  fVarNames_QA[13] = "recalPV_sigma_Y_rec";
  fVarNames_QA[14] = "recalPV_sigma_Z_rec";
  fVarNames_QA[15] = "PV_NContributors";
  fVarNames_QA[16] = "recalPV_NContributors";
  fVarNames_QA[17] = "SV_X_MC";
  fVarNames_QA[18] = "SV_Y_MC";
  fVarNames_QA[19] = "SV_Z_MC";
  fVarNames_QA[20] = "SV_X_rec";
  fVarNames_QA[21] = "SV_Y_rec";
  fVarNames_QA[22] = "SV_Z_rec";
  fVarNames_QA[23] = "SV_sigma_X_rec";
  fVarNames_QA[24] = "SV_sigma_Y_rec";
  fVarNames_QA[25] = "SV_sigma_Z_rec";
  fVarNames_QA[26] = "SV_X_rec_wTopoConst";
  fVarNames_QA[27] = "SV_Y_rec_wTopoConst";
  fVarNames_QA[28] = "SV_Z_rec_wTopoConst";
  fVarNames_QA[29] = "SV_sigma_X_rec_wTopoConst";
  fVarNames_QA[30] = "SV_sigma_Y_rec_wTopoConst";
  fVarNames_QA[31] = "SV_sigma_Z_rec_wTopoConst";
  fVarNames_QA[32] = "pt_XicPlus";
  fVarNames_QA[33] = "pt_XicPlus_wTopoConst";
  fVarNames_QA[34] = "SV_X_rec_wTopoConst_recalPV";
  fVarNames_QA[35] = "SV_Y_rec_wTopoConst_recalPV";
  fVarNames_QA[36] = "SV_Z_rec_wTopoConst_recalPV";
  fVarNames_QA[37] = "SV_sigma_X_rec_wTopoConst_recalPV";
  fVarNames_QA[38] = "SV_sigma_Y_rec_wTopoConst_recalPV";
  fVarNames_QA[39] = "SV_sigma_Z_rec_wTopoConst_recalPV";
  fVarNames_QA[40] = "pt_XicPlus_wTopoConst_recalPV";
  fVarNames_QA[41] = "PV_CountRealContributors";
  fVarNames_QA[42] = "PV_KF_Refit_X";
  fVarNames_QA[43] = "PV_KF_Refit_Y";
  fVarNames_QA[44] = "PV_KF_Refit_Z";
  fVarNames_QA[45] = "PV_KF_Refit_sigma_X";
  fVarNames_QA[46] = "PV_KF_Refit_sigma_Y";
  fVarNames_QA[47] = "PV_KF_Refit_sigma_Z";
  fVarNames_QA[48] = "recalPV_KF_Refit_X";
  fVarNames_QA[49] = "recalPV_KF_Refit_Y";
  fVarNames_QA[50] = "recalPV_KF_Refit_Z";
  fVarNames_QA[51] = "recalPV_KF_Refit_sigma_X";
  fVarNames_QA[52] = "recalPV_KF_Refit_sigma_Y";
  fVarNames_QA[53] = "recalPV_KF_Refit_sigma_Z";
  fVarNames_QA[54] = "recalPV_KF_X";
  fVarNames_QA[55] = "recalPV_KF_Y";
  fVarNames_QA[56] = "recalPV_KF_Z";
  fVarNames_QA[57] = "recalPV_KF_sigma_X";
  fVarNames_QA[58] = "recalPV_KF_sigma_Y";
  fVarNames_QA[59] = "recalPV_KF_sigma_Z";
  fVarNames_QA[60] = "Source_XicPlus";
  fVarNames_QA[61] = "event_ID";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    if (ivar<(nVar-1))  fTree_XicPlus_QA->Branch(fVarNames_QA[ivar].Data(), &fVar_XicPlus_QA[ivar], Form("%s/F", fVarNames_QA[ivar].Data()));
    if (ivar==(nVar-1)) fTree_XicPlus_QA->Branch(fVarNames_QA[ivar].Data(), &fVar_XicPlus_EvtID, Form("%s/l", fVarNames_QA[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeQAXicPlus_woMassConstForLamAndXi()
{
  const char* nameoutput = GetOutputSlot(8)->GetContainer()->GetName();
  fTree_XicPlus_QA_woMassConstForLamAndXi = new TTree(nameoutput, "XicPlus variables QA tree woMassConstForLamAndXi");
  Int_t nVar = 17;
  fVar_XicPlus_QA_woMassConstForLamAndXi = new Float_t[nVar-1];
  TString *fVarNames_QA_woMassConstForLamAndXi = new TString[nVar];

  fVarNames_QA_woMassConstForLamAndXi[0]  = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi
  fVarNames_QA_woMassConstForLamAndXi[1]  = "chi2topo_XiToXicPlus"; // chi2_topo of Xi to Xic+
  fVarNames_QA_woMassConstForLamAndXi[2]  = "chi2topo_XicPlusToPV"; // chi2_topo of Xic+ to PV
  fVarNames_QA_woMassConstForLamAndXi[3]  = "pt_XicPlus";
  fVarNames_QA_woMassConstForLamAndXi[4]  = "mass_XicPlus";
  fVarNames_QA_woMassConstForLamAndXi[5]  = "SV_X_rec";
  fVarNames_QA_woMassConstForLamAndXi[6]  = "SV_Y_rec";
  fVarNames_QA_woMassConstForLamAndXi[7]  = "SV_Z_rec";
  fVarNames_QA_woMassConstForLamAndXi[8]  = "SV_sigma_X_rec";
  fVarNames_QA_woMassConstForLamAndXi[9]  = "SV_sigma_Y_rec";
  fVarNames_QA_woMassConstForLamAndXi[10] = "SV_sigma_Z_rec";
  fVarNames_QA_woMassConstForLamAndXi[11] = "SV_X_MC";
  fVarNames_QA_woMassConstForLamAndXi[12] = "SV_Y_MC";
  fVarNames_QA_woMassConstForLamAndXi[13] = "SV_Z_MC";
  fVarNames_QA_woMassConstForLamAndXi[14] = "pt_XicPlus_MC";
  fVarNames_QA_woMassConstForLamAndXi[15] = "Source_XicPlus";
  fVarNames_QA_woMassConstForLamAndXi[16] = "event_ID";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    if (ivar<(nVar-1))  fTree_XicPlus_QA_woMassConstForLamAndXi->Branch(fVarNames_QA_woMassConstForLamAndXi[ivar].Data(), &fVar_XicPlus_QA_woMassConstForLamAndXi[ivar], Form("%s/F", fVarNames_QA_woMassConstForLamAndXi[ivar].Data()));
    if (ivar==(nVar-1)) fTree_XicPlus_QA_woMassConstForLamAndXi->Branch(fVarNames_QA_woMassConstForLamAndXi[ivar].Data(), &fVar_XicPlus_EvtID, Form("%s/l", fVarNames_QA_woMassConstForLamAndXi[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeQAXicPlus_wMassConstForLam_woMassConstForXi()
{
  const char* nameoutput = GetOutputSlot(9)->GetContainer()->GetName();
  fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi = new TTree(nameoutput, "XicPlus variables QA tree wMassConstForLam_woMassConstForXi");
  Int_t nVar = 17;
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi = new Float_t[nVar-1];
  TString *fVarNames_QA_wMassConstForLam_woMassConstForXi = new TString[nVar];

  fVarNames_QA_wMassConstForLam_woMassConstForXi[0]  = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi
  fVarNames_QA_wMassConstForLam_woMassConstForXi[1]  = "chi2topo_XiToXicPlus"; // chi2_topo of Xi to Xic+
  fVarNames_QA_wMassConstForLam_woMassConstForXi[2]  = "chi2topo_XicPlusToPV"; // chi2_topo of Xic+ to PV
  fVarNames_QA_wMassConstForLam_woMassConstForXi[3]  = "pt_XicPlus";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[4]  = "mass_XicPlus";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[5]  = "SV_X_rec";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[6]  = "SV_Y_rec";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[7]  = "SV_Z_rec";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[8]  = "SV_sigma_X_rec";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[9]  = "SV_sigma_Y_rec";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[10] = "SV_sigma_Z_rec";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[11] = "SV_X_MC";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[12] = "SV_Y_MC";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[13] = "SV_Z_MC";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[14] = "pt_XicPlus_MC";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[15] = "Source_XicPlus";
  fVarNames_QA_wMassConstForLam_woMassConstForXi[16] = "event_ID";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    if (ivar<(nVar-1))  fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi->Branch(fVarNames_QA_wMassConstForLam_woMassConstForXi[ivar].Data(), &fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[ivar], Form("%s/F", fVarNames_QA_wMassConstForLam_woMassConstForXi[ivar].Data()));
    if (ivar==(nVar-1)) fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi->Branch(fVarNames_QA_wMassConstForLam_woMassConstForXi[ivar].Data(), &fVar_XicPlus_EvtID, Form("%s/l", fVarNames_QA_wMassConstForLam_woMassConstForXi[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeQAXicPlus_woMassConstForLam_wMassConstForXi()
{
  const char* nameoutput = GetOutputSlot(10)->GetContainer()->GetName();
  fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi = new TTree(nameoutput, "XicPlus variables QA tree woMassConstForLam_wMassConstForXi");
  Int_t nVar = 17;
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi = new Float_t[nVar-1];
  TString *fVarNames_QA_woMassConstForLam_wMassConstForXi = new TString[nVar];

  fVarNames_QA_woMassConstForLam_wMassConstForXi[0]  = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi
  fVarNames_QA_woMassConstForLam_wMassConstForXi[1]  = "chi2topo_XiToXicPlus"; // chi2_topo of Xi to Xic+
  fVarNames_QA_woMassConstForLam_wMassConstForXi[2]  = "chi2topo_XicPlusToPV"; // chi2_topo of Xic+ to PV
  fVarNames_QA_woMassConstForLam_wMassConstForXi[3]  = "pt_XicPlus";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[4]  = "mass_XicPlus";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[5]  = "SV_X_rec";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[6]  = "SV_Y_rec";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[7]  = "SV_Z_rec";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[8]  = "SV_sigma_X_rec";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[9]  = "SV_sigma_Y_rec";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[10] = "SV_sigma_Z_rec";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[11] = "SV_X_MC";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[12] = "SV_Y_MC";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[13] = "SV_Z_MC";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[14] = "pt_XicPlus_MC";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[15] = "Source_XicPlus";
  fVarNames_QA_woMassConstForLam_wMassConstForXi[16] = "event_ID";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    if (ivar<(nVar-1))  fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi->Branch(fVarNames_QA_woMassConstForLam_wMassConstForXi[ivar].Data(), &fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[ivar], Form("%s/F", fVarNames_QA_woMassConstForLam_wMassConstForXi[ivar].Data()));
    if (ivar==(nVar-1)) fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi->Branch(fVarNames_QA_woMassConstForLam_wMassConstForXi[ivar].Data(), &fVar_XicPlus_EvtID, Form("%s/l", fVarNames_QA_woMassConstForLam_wMassConstForXi[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeQAXicPlus_wMassConstForLamAndXi()
{
  const char* nameoutput = GetOutputSlot(11)->GetContainer()->GetName();
  fTree_XicPlus_QA_wMassConstForLamAndXi = new TTree(nameoutput, "XicPlus variables QA tree wMassConstForLamAndXi");
  Int_t nVar = 17;
  fVar_XicPlus_QA_wMassConstForLamAndXi = new Float_t[nVar-1];
  TString *fVarNames_QA_wMassConstForLamAndXi = new TString[nVar];

  fVarNames_QA_wMassConstForLamAndXi[0]  = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi
  fVarNames_QA_wMassConstForLamAndXi[1]  = "chi2topo_XiToXicPlus"; // chi2_topo of Xi to Xic+
  fVarNames_QA_wMassConstForLamAndXi[2]  = "chi2topo_XicPlusToPV"; // chi2_topo of Xic+ to PV
  fVarNames_QA_wMassConstForLamAndXi[3]  = "pt_XicPlus";
  fVarNames_QA_wMassConstForLamAndXi[4]  = "mass_XicPlus";
  fVarNames_QA_wMassConstForLamAndXi[5]  = "SV_X_rec";
  fVarNames_QA_wMassConstForLamAndXi[6]  = "SV_Y_rec";
  fVarNames_QA_wMassConstForLamAndXi[7]  = "SV_Z_rec";
  fVarNames_QA_wMassConstForLamAndXi[8]  = "SV_sigma_X_rec";
  fVarNames_QA_wMassConstForLamAndXi[9]  = "SV_sigma_Y_rec";
  fVarNames_QA_wMassConstForLamAndXi[10] = "SV_sigma_Z_rec";
  fVarNames_QA_wMassConstForLamAndXi[11] = "SV_X_MC";
  fVarNames_QA_wMassConstForLamAndXi[12] = "SV_Y_MC";
  fVarNames_QA_wMassConstForLamAndXi[13] = "SV_Z_MC";
  fVarNames_QA_wMassConstForLamAndXi[14] = "pt_XicPlus_MC";
  fVarNames_QA_wMassConstForLamAndXi[15] = "Source_XicPlus";
  fVarNames_QA_wMassConstForLamAndXi[16] = "event_ID";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    if (ivar<(nVar-1))  fTree_XicPlus_QA_wMassConstForLamAndXi->Branch(fVarNames_QA_wMassConstForLamAndXi[ivar].Data(), &fVar_XicPlus_QA_wMassConstForLamAndXi[ivar], Form("%s/F", fVarNames_QA_wMassConstForLamAndXi[ivar].Data()));
    if (ivar==(nVar-1)) fTree_XicPlus_QA_wMassConstForLamAndXi->Branch(fVarNames_QA_wMassConstForLamAndXi[ivar].Data(), &fVar_XicPlus_EvtID, Form("%s/l", fVarNames_QA_wMassConstForLamAndXi[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeQAXicPlus_wMassAndTopoConstForLam_wMassConstForXi()
{
  const char* nameoutput = GetOutputSlot(12)->GetContainer()->GetName();
  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi = new TTree(nameoutput, "XicPlus variables QA tree wMassAndTopoConstForLam_wMassConstForXi");
  Int_t nVar = 17;
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi = new Float_t[nVar-1];
  TString *fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi = new TString[nVar];

  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[0]  = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[1]  = "chi2mass_Xi"; // chi2_mass of Xi
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[2]  = "chi2topo_XicPlusToPV"; // chi2_topo of Xic+ to PV
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[3]  = "pt_XicPlus";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[4]  = "mass_XicPlus";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[5]  = "SV_X_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[6]  = "SV_Y_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[7]  = "SV_Z_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[8]  = "SV_sigma_X_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[9]  = "SV_sigma_Y_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[10] = "SV_sigma_Z_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[11] = "SV_X_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[12] = "SV_Y_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[13] = "SV_Z_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[14] = "pt_XicPlus_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[15] = "Source_XicPlus";
  fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[16] = "event_ID";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    if (ivar<(nVar-1))  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi->Branch(fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[ivar].Data(), &fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[ivar], Form("%s/F", fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[ivar].Data()));
    if (ivar==(nVar-1)) fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi->Branch(fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[ivar].Data(), &fVar_XicPlus_EvtID, Form("%s/l", fVarNames_QA_wMassAndTopoConstForLam_wMassConstForXi[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeQAXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi()
{
  const char* nameoutput = GetOutputSlot(13)->GetContainer()->GetName();
  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi = new TTree(nameoutput, "XicPlus variables QA tree wMassAndTopoConstForLam_wMassAndTopoConstForXi");
  Int_t nVar = 17;
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi = new Float_t[nVar-1];
  TString *fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi = new TString[nVar];

  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[0]  = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[1]  = "chi2topo_XiToPV"; // chi2_mass of Xi
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[2]  = "chi2topo_XicPlusToPV"; // chi2_topo of Xic+ to PV
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[3]  = "pt_XicPlus";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[4]  = "mass_XicPlus";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[5]  = "SV_X_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[6]  = "SV_Y_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[7]  = "SV_Z_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[8]  = "SV_sigma_X_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[9]  = "SV_sigma_Y_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[10] = "SV_sigma_Z_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[11] = "SV_X_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[12] = "SV_Y_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[13] = "SV_Z_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[14] = "pt_XicPlus_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[15] = "Source_XicPlus";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[16] = "event_ID";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    if (ivar<(nVar-1))  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi->Branch(fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[ivar].Data(), &fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[ivar], Form("%s/F", fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[ivar].Data()));
    if (ivar==(nVar-1)) fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi->Branch(fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[ivar].Data(), &fVar_XicPlus_EvtID, Form("%s/l", fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeQAXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic()
{
  const char* nameoutput = GetOutputSlot(14)->GetContainer()->GetName();
  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic = new TTree(nameoutput, "XicPlus variables QA tree wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic");
  Int_t nVar = 32;
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic = new Float_t[nVar-1];
  TString *fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic = new TString[nVar];

  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[0]  = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[1]  = "chi2topo_XiToPV"; // chi2_mass of Xi
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[2]  = "chi2topo_XicPlusToPV"; // chi2_topo of Xic+ to PV
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[3]  = "pt_XicPlus";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[4]  = "mass_XicPlus";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[5]  = "SV_X_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[6]  = "SV_Y_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[7]  = "SV_Z_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[8]  = "SV_sigma_X_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[9]  = "SV_sigma_Y_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[10] = "SV_sigma_Z_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[11] = "SV_X_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[12] = "SV_Y_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[13] = "SV_Z_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[14] = "PV_X_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[15] = "PV_Y_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[16] = "PV_Z_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[17] = "PV_X_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[18] = "PV_Y_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[19] = "PV_Z_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[20] = "PV_sigma_X_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[21] = "PV_sigma_Y_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[22] = "PV_sigma_Z_rec";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[23] = "recalPV_woKFrefit_X";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[24] = "recalPV_woKFrefit_Y";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[25] = "recalPV_woKFrefit_Z";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[26] = "recalPV_woKFrefit_sigma_X";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[27] = "recalPV_woKFrefit_sigma_Y";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[28] = "recalPV_woKFrefit_sigma_Z";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[29] = "pt_XicPlus_MC";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[30] = "Source_XicPlus";
  fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[31] = "event_ID";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    if (ivar<(nVar-1))  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic->Branch(fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[ivar].Data(), &fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[ivar], Form("%s/F", fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[ivar].Data()));
    if (ivar==(nVar-1)) fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic->Branch(fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[ivar].Data(), &fVar_XicPlus_EvtID, Form("%s/l", fVarNames_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeRecXicPlus()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fTree_XicPlus = new TTree(nameoutput, "XicPlus variables tree");
  Int_t nVar = 66;
  fVar_XicPlus = new Float_t[nVar-1];
  TString *fVarNames = new TString[nVar];

  fVarNames[0]  = "nSigmaTPC_Pi0FromXicPlus"; // TPC nsigma for pion_0 coming from XicPlus
  fVarNames[1]  = "nSigmaTOF_Pi0FromXicPlus"; // TOF nsigma for pion_0 coming from XicPlus
  fVarNames[2]  = "nSigmaTPC_Pi1FromXicPlus"; // TPC nsigma for pion_1 coming from XicPlus
  fVarNames[3]  = "nSigmaTOF_Pi1FromXicPlus"; // TOF nsigma for pion_1 coming from XicPlus
  fVarNames[4]  = "nSigmaTPC_PiFromXi"; // TPC nsigma for pion coming from Xi
  fVarNames[5]  = "nSigmaTOF_PiFromXi"; // TOF nsigma for pion coming from Xi
  fVarNames[6]  = "nSigmaTPC_PiFromLam"; // TPC nsigma for pion coming from Lambda
  fVarNames[7]  = "nSigmaTPC_PrFromLam"; // TPC nsigma for proton coming from Lambda

  fVarNames[8]  = "chi2geo_Lam"; // chi2_geometry of Lambda (without mass constraint)
  fVarNames[9]  = "ldl_Lam"; // l/dl of Lambda
  fVarNames[10] = "chi2topo_LamToPV"; // chi2_topo of Lambda (with mass constraint) to PV

  fVarNames[11] = "chi2geo_Xi"; // chi2_geometry of Xi (with Lambda mass const.)
  fVarNames[12] = "ldl_Xi"; // l/dl of Xi (with Lambda mass const.)
  fVarNames[13] = "chi2topo_XiToPV"; // chi2_topo of Xi (with mass constraint) to PV
  fVarNames[14] = "chi2MassConst_Xi"; // chi2_MassConst of Xi

  fVarNames[15] = "DecayLxy_Lam"; // decay length of Lambda in x-y plane
  fVarNames[16] = "ct_Lam"; // life time of Lambda
  fVarNames[17] = "DecayLxy_Xi"; // decay length of Xi in x-y plane
  fVarNames[18] = "ct_Xi"; // life time of Xi

  fVarNames[19] = "PA_LamToXi"; // pointing angle of Lmabda (pointing back to Xi)
  fVarNames[20] = "PA_LamToPV"; // pointing angle of Lambda (pointing back to PV)
  fVarNames[21] = "PA_XiToPV"; // pointing angle of Xi (pointing back to PV)

  fVarNames[22] = "mass_Lam"; // mass of Lambda (without mass const.)
  fVarNames[23] = "mass_Xi"; // mass of Xi (without mass const.)

  fVarNames[24] = "pt_Pi0FromXicPlus"; // pt of pion_0 coming from XicPlus
  fVarNames[25] = "pt_Pi1FromXicPlus"; // pt of pion_1 coming from XicPlus

  fVarNames[26] = "pt_XicPlus"; // pt of XicPlus
  fVarNames[27] = "rap_XicPlus"; // rapidity of XicPlus
  fVarNames[28] = "mass_XicPlus"; // mass of XicPlus
  fVarNames[29] = "chi2geo_XicPlus"; // chi2_geometry of XicPlus

  fVarNames[30] = "chi2prim_Pi0FromXicPlus"; // chi2_topo of pion_0 to PV
  fVarNames[31] = "chi2prim_Pi1FromXicPlus"; // chi2_topo of pion_1 to PV
  fVarNames[32] = "DCAxy_Pi0FromXicPlusToPV_KF"; // DCA of pion_0 coming from XicPlus in x-y plane
  fVarNames[33] = "DCAxy_Pi1FromXicPlusToPV_KF"; // DCA of pion_1 coming from XicPlus in x-y plane
  fVarNames[34] = "mass_K0S"; // mass of Ks0
  fVarNames[35] = "mass_Gamma"; // mass of e+e-

  fVarNames[36] = "DCAxy_LamDau"; // DCA of Lam's daughters
  fVarNames[37] = "DCAxy_XiDau"; // DCA of Xi's daughters (calculated from KF after Lambda mass constraint)
  fVarNames[38] = "DCAxy_XiToPV"; // DCA of Xi to PV in x-y plane (calculated from KF after Xi mass constraint)
  fVarNames[39] = "DCAxy_PiToPi"; // DCA of pi to pi in x-y plane
  fVarNames[40] = "DCA_PiToPi"; // DCA of pi to pi
  fVarNames[41] = "DCAxy_Pi0ToXi"; // DCA of pion_0 to Xi in x-y plane
  fVarNames[42] = "DCAxy_Pi1ToXi"; // DCA of pion_1 to Xi in x-y plane
  fVarNames[43] = "PA_XicPlusToPV"; // pointing angle of XicPlus (pointing back to PV)
  fVarNames[44] = "DecayLxy_XicPlus"; // decay length of XicPlus in x-y plane
  fVarNames[45] = "chi2topo_XicPlus"; // chi2_topo of XicPlus to PV
  fVarNames[46] = "ldl_XicPlus"; // l/dl of XicPlus
  fVarNames[47] = "ct_XicPlus"; // lifetime of XicPlus
  fVarNames[48] = "PA_XiToXicPlus"; // pointing angle of Xi (pointing back to XicPlus)
  fVarNames[49] = "pt_PiFromXi"; // pt of pion coming from Xi
  fVarNames[50] = "mass_Omega"; // mass of Omega
  fVarNames[51] = "PA_XicPlusToRecalPVfromKF_Refit"; // pointing angle of XicPlus (pointing back to recalPV from KF) (Refit)
  fVarNames[52] = "chi2topo_XicPlusToRecalPVfromKF_Refit"; // chi2_topo of XicPlus to recalPV from KF (Refit)
  fVarNames[53] = "PA_XicPlusToRecalPVfromKF"; // pointing angle of XicPlus (pointing back to recalPV from KF)
  fVarNames[54] = "chi2topo_XicPlusToRecalPVfromKF"; // chi2_topo of XicPlus to recalPV from KF
  fVarNames[55] = "PA_XicPlusToRecalPV"; // pointing angle of XicPlus (pointing back to recalPV)
  fVarNames[56] = "chi2topo_XicPlusToRecalPV"; // chi2_topo of XicPlus to recalPV
  fVarNames[57] = "PAXY_XicPlusToPV"; // pointing angle (x-Y) of XicPlus (pointing back to PV)
  fVarNames[58] = "PAXY_XicPlusToRecalPVKF_Refit"; // pointing angle (X-Y) of XicPlus (pointing back to recalPV from KF) (Refit)
  fVarNames[59] = "PAXY_XicPlusToRecalPVfromKF"; // pointing angle (X-Y) of XicPlus (pointing back to recalPV from KF)
  fVarNames[60] = "PAXY_XicPlusToRecalPV"; // pointing angle (X-Y) of XicPlus (pointing back to recalPV)
  fVarNames[61] = "PA_XicPlusToRecalPVfromKF_Refit_woAddMother"; // pointing angle of XicPlus without adding XicPlus into PV
  fVarNames[62] = "chi2topo_XicPlusToRecalPVfromKF_Refit_woAddMother";
  fVarNames[63] = "PV_NContributors"; // number of tracks used for PV fit + 1
  fVarNames[64] = "Source_XicPlus"; // flag for XicPlus MC truth (“4” prompt, "5" feed-down, “<0” background)
  fVarNames[65] = "event_ID"; // event ID

//  fVarNames[26] = "CosThetaStar_PiFromXicPlus"; // CosThetaStar of pion coming from XicPlus
//  fVarNames[27] = "CosThetaStar_Xi"; // CosThetaStar of Xi coming from XicPlus
//  fVarNames[41] = "DCA_XiDau_Cascade"; // DCA of Xi's daughters (calculated from AOD cascade)
//  fVarNames[33] = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi
//  fVarNames[34] = "chi2topo_XiToXicPlus"; // chi2_topo of Xi to XicPlus
//  fVarNames[35] = "DecayLxy_XicPlus"; // decay length of XicPlus in x-y plane
//  fVarNames[36] = "ct_XicPlus"; // life time of XicPlus
//  fVarNames[38] = "DCA_LamDau"; // DCA of Lambda's daughters (calculated from AOD cascade)
//  fVarNames[40] = "DCA_XicPlusDau_KF"; // DCA of XicPlus's daughters (calculated from KF after Xi mass constraint)

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    if (ivar<(nVar-1))  fTree_XicPlus->Branch(fVarNames[ivar].Data(), &fVar_XicPlus[ivar], Form("%s/F", fVarNames[ivar].Data()));
    if (ivar==(nVar-1)) fTree_XicPlus->Branch(fVarNames[ivar].Data(), &fVar_XicPlus_EvtID, Form("%s/l", fVarNames[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineTreeGenXicPlus()
{
  const char* nameoutput = GetOutputSlot(5)->GetContainer()->GetName();
  fTree_XicPlusMCGen = new TTree(nameoutput,"XicPlus MC variables tree");
  Int_t nVar = 5;
  fVar_XicPlusMCGen = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];
  
  fVarNames[0] = "rap_XicPlus";
  fVarNames[1] = "pt_XicPlus";
  fVarNames[2] = "Source_XicPlus";
  fVarNames[3] = "PDG_XicPlus";
  fVarNames[4] = "MLoverP"; // c*(proper lifetime)

  /*
  fVarNames[ 0]="Centrality";
  fVarNames[ 1]="DecayType";
  fVarNames[ 2]="XicSource";
  fVarNames[ 3]="XicEta";
  fVarNames[ 4]="XicY";
  fVarNames[ 5]="XicPx";
  fVarNames[ 6]="XicPy";
  fVarNames[ 7]="XicPz";
  fVarNames[ 8]="PiPx";
  fVarNames[ 9]="PiPy";
  fVarNames[10]="PiPz";
  fVarNames[11]="CascPx";
  fVarNames[12]="CascPy";
  fVarNames[13]="CascPz";
  fVarNames[14]="XicPdgCode";
  fVarNames[15]="PiPdgCode";
  fVarNames[16]="CascPdgCode";
  fVarNames[17]="RunNumber";
  fVarNames[18]="EvNumber";
  fVarNames[19]="IsPrompt";
  */

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_XicPlusMCGen->Branch(fVarNames[ivar].Data(),&fVar_XicPlusMCGen[ivar],Form("%s/F",fVarNames[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::InvMassV0atPV(AliAODTrack *trk1, AliAODTrack *trk2, Int_t pdg1, Int_t pdg2)
{
  
  Double_t mass1 = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
  Double_t mass2 = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();
  Double_t E1 = TMath::Sqrt(mass1*mass1 + trk1->P()*trk1->P());
  Double_t E2 = TMath::Sqrt(mass2*mass2 + trk2->P()*trk2->P());
  Double_t mass = TMath::Sqrt( (E1+E2)*(E1+E2) - (trk1->Px()+trk2->Px())*(trk1->Px()+trk2->Px()) - (trk1->Py()+trk2->Py())*(trk1->Py()+trk2->Py()) - (trk1->Pz()+trk2->Pz())*(trk1->Pz()+trk2->Pz()) );

  return mass;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::DefineAnaHist()
{
  // Define analysis histograms
  
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillEventROOTObjects()
{

  for (Int_t i=0; i<7; i++) {
    fVar_Event[i] = -9999.;
  }

  Double_t pos[3];
  fpVtx->GetXYZ(pos);
  Double_t cov[6];
  //fpVtx->GetSigmaXYZ(sigma);
  fpVtx->GetCovarianceMatrix(cov);

  fVar_Event[1] = pos[0];
  fVar_Event[2] = pos[1];
  fVar_Event[3] = pos[2];
  if (cov[0]>=0 && cov[2]>=0 && cov[5]>=0) {
    fVar_Event[4] = sqrt(cov[0]);
    fVar_Event[5] = sqrt(cov[2]);
    fVar_Event[6] = sqrt(cov[5]);
  }
  fVar_Event[7] = fpVtx->GetNContributors();

  fTree_Event->Fill();

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillTreeRecXicPlusFromCasc(AliAODEvent *AODEvent, AliAODcascade *casc, KFParticle kfpXicPlus, AliAODTrack *trackPiFromXicPlus_trk1, KFParticle kfpBP_trk1, KFParticle kfpXiMinus, KFParticle kfpXiMinus_m, KFParticle kfpPionOrKaon, AliAODTrack *trackPiFromXiOrKaonFromOmega, KFParticle kfpK0Short, KFParticle kfpGamma, KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkProton, AliAODTrack *trkPion, AliAODTrack *trackPiFromXicPlus_trk2, KFParticle kfpBP_trk2, KFParticle kfpProtonFromLam, KFParticle kfpPionFromLam, KFParticle PV, KFParticle PV_KF_Refit, TClonesArray *mcArray, Int_t lab_XicPlus)
{

  for (Int_t i=0; i<(66-1); i++) {
    fVar_XicPlus[i] = -9999.;
  }

  Double_t nSigmaTOF_PiFromXicPlus_trk1 = -999., nSigmaTOF_PiFromXicPlus_trk2 = -999., nSigmaTOF_PiFromXi = -999.;
  AliAODPidHF *Pid_HF = fAnaCuts->GetPidHF();
  if (Pid_HF) {
    Pid_HF->GetnSigmaTOF(trackPiFromXicPlus_trk1, AliPID::kPion, nSigmaTOF_PiFromXicPlus_trk1);
    Pid_HF->GetnSigmaTOF(trackPiFromXicPlus_trk2, AliPID::kPion, nSigmaTOF_PiFromXicPlus_trk2);
    Pid_HF->GetnSigmaTOF(trackPiFromXiOrKaonFromOmega, AliPID::kPion, nSigmaTOF_PiFromXi);
  }

  Float_t nSigmaTPC_PiFromXicPlus_trk1 = fPID->NumberOfSigmasTPC(trackPiFromXicPlus_trk1,AliPID::kPion);
//  Float_t nSigmaTOF_PiFromXicPlus_trk1 = fPID->NumberOfSigmasTOF(trackPiFromXicPlus_trk1,AliPID::kPion);
  Float_t nSigmaTPC_PiFromXicPlus_trk2 = fPID->NumberOfSigmasTPC(trackPiFromXicPlus_trk2,AliPID::kPion);
//  Float_t nSigmaTOF_PiFromXicPlus_trk2 = fPID->NumberOfSigmasTOF(trackPiFromXicPlus_trk2,AliPID::kPion);

  Float_t nSigmaTPC_PrFromLam  = fPID->NumberOfSigmasTPC(trkProton,AliPID::kProton);
  Float_t nSigmaTPC_PiFromLam  = fPID->NumberOfSigmasTPC(trkPion,AliPID::kPion);

  Float_t nSigmaTPC_PiFromXi = -9999.;
  nSigmaTPC_PiFromXi = fPID->NumberOfSigmasTPC(trackPiFromXiOrKaonFromOmega,AliPID::kPion);
//  nSigmaTOF_PiFromXi = fPID->NumberOfSigmasTOF(trackPiFromXiOrKaonFromOmega,AliPID::kPion);

  if ( fabs(nSigmaTPC_PiFromXicPlus_trk1)>=4. || fabs(nSigmaTPC_PiFromXicPlus_trk2)>=4. || fabs(nSigmaTPC_PiFromXi)>=4. || fabs(nSigmaTPC_PrFromLam)>=4. || fabs(nSigmaTPC_PiFromLam)>=4. ) return;

//  AliAODTrack *trk0 = (AliAODTrack*) (v0->GetDaughter(0));

//  Double_t alpha_FirstDaugPos = (v0->MomPosAlongV0() - v0->MomNegAlongV0())/(v0->MomPosAlongV0() + v0->MomNegAlongV0());

//  if (trk0->Charge()<0) {
//    alpha_FirstDaugPos = (v0->MomNegAlongV0() - v0->MomPosAlongV0())/(v0->MomPosAlongV0() + v0->MomNegAlongV0());
//  }

  KFParticle kfpXicPlus_PV = kfpXicPlus;
  kfpXicPlus_PV.SetProductionVertex(PV);

  KFParticle kfpXiMinus_XicPlus = kfpXiMinus_m;
  kfpXiMinus_XicPlus.SetProductionVertex(kfpXicPlus);
  KFParticle kfpXiMinus_PV = kfpXiMinus_m;
  kfpXiMinus_PV.SetProductionVertex(PV);

  KFParticle kfpLambda_Xi = kfpLambda_m;
  kfpLambda_Xi.SetProductionVertex(kfpXiMinus);
  KFParticle kfpLambda_PV = kfpLambda_m;
  kfpLambda_PV.SetProductionVertex(PV);
//  KFParticle kfpXiMinus_pv = kfpXiMinus;
//  kfpXiMinus_pv.SetProductionVertex(PV);

//  KFParticle kfpBP_XicPlus = kfpBP_trk1;
//  kfpBP_XicPlus.SetProductionVertex(kfpXicPlus);

  kfpBP_trk1.GetDistanceFromVertexXY(PV);
  kfpBP_trk2.GetDistanceFromVertexXY(PV);

  //************************** calculate l/Δl for Lambda *************************************
  Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
  Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
  Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
  Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
  Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
  if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
  dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
  if ( dl_Lambda<=0 ) return;
  //***************************************************************************************
  //************************** calculate l/Δl for Xi- *************************************
  Double_t dx_Xi = PV.GetX()-kfpXiMinus.GetX();
  Double_t dy_Xi = PV.GetY()-kfpXiMinus.GetY();
  Double_t dz_Xi = PV.GetZ()-kfpXiMinus.GetZ();
  Double_t l_Xi = TMath::Sqrt(dx_Xi*dx_Xi + dy_Xi*dy_Xi + dz_Xi*dz_Xi);
  Double_t dl_Xi = (PV.GetCovariance(0)+kfpXiMinus.GetCovariance(0))*dx_Xi*dx_Xi + (PV.GetCovariance(2)+kfpXiMinus.GetCovariance(2))*dy_Xi*dy_Xi + (PV.GetCovariance(5)+kfpXiMinus.GetCovariance(5))*dz_Xi*dz_Xi + 2*( (PV.GetCovariance(1)+kfpXiMinus.GetCovariance(1))*dx_Xi*dy_Xi + (PV.GetCovariance(3)+kfpXiMinus.GetCovariance(3))*dx_Xi*dz_Xi + (PV.GetCovariance(4)+kfpXiMinus.GetCovariance(4))*dy_Xi*dz_Xi );
  if ( fabs(l_Xi)<1.e-8f ) l_Xi = 1.e-8f;
  dl_Xi = dl_Xi<0. ? 1.e8f : sqrt(dl_Xi)/l_Xi;
  if ( dl_Xi<=0 ) return;
  //***************************************************************************************
  //************************** calculate l/Δl for Xic+ *************************************
  Double_t dx_XicPlus = PV.GetX()-kfpXicPlus.GetX();
  Double_t dy_XicPlus = PV.GetY()-kfpXicPlus.GetY();
  Double_t dz_XicPlus = PV.GetZ()-kfpXicPlus.GetZ();
  Double_t l_XicPlus = TMath::Sqrt(dx_XicPlus*dx_XicPlus + dy_XicPlus*dy_XicPlus + dz_XicPlus*dz_XicPlus);
  Double_t dl_XicPlus = (PV.GetCovariance(0)+kfpXicPlus.GetCovariance(0))*dx_XicPlus*dx_XicPlus + (PV.GetCovariance(2)+kfpXicPlus.GetCovariance(2))*dy_XicPlus*dy_XicPlus + (PV.GetCovariance(5)+kfpXicPlus.GetCovariance(5))*dz_XicPlus*dz_XicPlus + 2*( (PV.GetCovariance(1)+kfpXicPlus.GetCovariance(1))*dx_XicPlus*dy_XicPlus + (PV.GetCovariance(3)+kfpXicPlus.GetCovariance(3))*dx_XicPlus*dz_XicPlus + (PV.GetCovariance(4)+kfpXicPlus.GetCovariance(4))*dy_XicPlus*dz_XicPlus );
  if ( fabs(l_XicPlus)<1.e-8f ) l_XicPlus = 1.e-8f;
  dl_XicPlus = dl_XicPlus<0. ? 1.e8f : sqrt(dl_XicPlus)/l_XicPlus;
  if ( dl_XicPlus<=0 ) return;
  //***************************************************************************************

  if ( kfpLambda_PV.GetChi2()/kfpLambda_PV.GetNDF() <= fAnaCuts->GetKFPLam_Chi2topoMin() ) return;
  if ( kfpXiMinus_PV.GetChi2()/kfpXiMinus_PV.GetNDF() >= fAnaCuts->GetKFPXi_Chi2topoMax() ) return;
  if ( l_Xi/dl_Xi <= fAnaCuts->GetKFPXi_lDeltalMin() ) return;

  const Float_t PDGmassXicPlus = TDatabasePDG::Instance()->GetParticle(4232)->Mass();
  Float_t mass_XicPlus_PV, err_mass_XicPlus_PV;
  kfpXicPlus_PV.GetMass(mass_XicPlus_PV, err_mass_XicPlus_PV);
  fVar_XicPlus[28] = mass_XicPlus_PV; // mass of XicPlus

  if ( fabs(mass_XicPlus_PV-PDGmassXicPlus) > fAnaCuts->GetProdMassTolXicPlus() ) return;


  fVar_XicPlus[0]  = nSigmaTPC_PiFromXicPlus_trk1; // TPC nsigma for pion_0 coming from XicPlus
  fVar_XicPlus[1]  = nSigmaTOF_PiFromXicPlus_trk1; // TOF nsigma for pion_0 coming from XicPlus
  fVar_XicPlus[2]  = nSigmaTPC_PiFromXicPlus_trk2; // TPC nsigma for pion_1 coming from XicPlus
  fVar_XicPlus[3]  = nSigmaTOF_PiFromXicPlus_trk2; // TOF nsigma for pion_1 coming from XicPlus
  fVar_XicPlus[4]  = nSigmaTPC_PiFromXi; // TPC nsigma for pion coming from Xi
  fVar_XicPlus[5]  = nSigmaTOF_PiFromXi; // TOF nsigma for pion coming from Xi
  fVar_XicPlus[6]  = nSigmaTPC_PiFromLam; // TPC nsigma for pion coming from Lambda
  fVar_XicPlus[7]  = nSigmaTPC_PrFromLam; // TPC nsigma for proton coming from Lambda

  fVar_XicPlus[8] = kfpLambda.GetChi2()/kfpLambda.GetNDF(); // chi2geo_Lam
  fVar_XicPlus[9] = l_Lambda/dl_Lambda; // l/dl of Lambda
  fVar_XicPlus[10] = kfpLambda_PV.GetChi2()/kfpLambda_PV.GetNDF(); // chi2_topo of Lambda (with mass constraint) to PV

  fVar_XicPlus[11] = kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF(); // chi2_geometry of Xi (with Lambda mass const.)
  fVar_XicPlus[12] = l_Xi/dl_Xi; // l/dl of Xi (with Lambda mass const.)
  fVar_XicPlus[13] = kfpXiMinus_PV.GetChi2()/kfpXiMinus_PV.GetNDF(); // chi2_topo of Xi (with mass constraint) to PV
  fVar_XicPlus[14] = kfpXiMinus_m.GetChi2()/kfpXiMinus_m.GetNDF(); // chi2_MassConst of Xi

  Float_t DecayLxy_Lam, err_DecayLxy_Lam;
  kfpLambda_Xi.GetDecayLengthXY(DecayLxy_Lam, err_DecayLxy_Lam);
  fVar_XicPlus[15] = DecayLxy_Lam; // decay length of Lambda in x-y plane
  Float_t ct_Lam=0., err_ct_Lam=0.;
  kfpLambda_Xi.GetLifeTime(ct_Lam, err_ct_Lam);
  fVar_XicPlus[16] = ct_Lam; // life time of Lambda

  Float_t DecayLxy_Xi, err_DecayLxy_Xi;
  kfpXiMinus_PV.GetDecayLengthXY(DecayLxy_Xi, err_DecayLxy_Xi);
  fVar_XicPlus[17] = DecayLxy_Xi; // decay length of Xi in x-y plane
  Float_t ct_Xi=0., err_ct_Xi=0.;
  kfpXiMinus_PV.GetLifeTime(ct_Xi, err_ct_Xi);
  fVar_XicPlus[18] = ct_Xi; // life time of Xi

  // calculate CosPointingAngle
  Double_t cosPA_v0toXi      = AliVertexingHFUtils::CosPointingAngleFromKF(kfpLambda_m, kfpXiMinus);
  Double_t cosPA_v0toPV      = AliVertexingHFUtils::CosPointingAngleFromKF(kfpLambda_m, PV);
  Double_t cosPA_XiToPV      = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXiMinus_m, PV);
  Double_t cosPA_XiToXicPlus = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXiMinus_m, kfpXicPlus);
  fVar_XicPlus[19] = TMath::ACos(cosPA_v0toXi); // pointing angle of Lmabda (pointing back to Xi)
  fVar_XicPlus[20] = TMath::ACos(cosPA_v0toPV); // pointing angle of Lambda (pointing back to PV)
  fVar_XicPlus[21] = TMath::ACos(cosPA_XiToPV); // pointing angle of Xi (pointing back to PV)

  Float_t mass_Lam_Rec, err_mass_Lam_Rec;
  kfpLambda.GetMass(mass_Lam_Rec, err_mass_Lam_Rec);
  fVar_XicPlus[22] = mass_Lam_Rec; // mass of Lambda (without mass const.)

  Float_t mass_Xi_Rec, err_mass_Xi_Rec;
  kfpXiMinus.GetMass(mass_Xi_Rec, err_mass_Xi_Rec);
  fVar_XicPlus[23] = mass_Xi_Rec; // mass of Xi (without mass const.)

  fVar_XicPlus[24] = trackPiFromXicPlus_trk1->Pt(); // pt of pion_0 coming from XicPlus
  fVar_XicPlus[25] = trackPiFromXicPlus_trk2->Pt(); // pt of pion_1 coming from XicPlus

  fVar_XicPlus[26] = kfpXicPlus_PV.GetPt(); // pt of XicPlus
  if ( TMath::Abs(kfpXicPlus_PV.GetE())>TMath::Abs(kfpXicPlus_PV.GetPz()) ) {
    fVar_XicPlus[27] = kfpXicPlus_PV.GetRapidity(); // rapidity of XicPlus
  }

  fVar_XicPlus[29] = kfpXicPlus.GetChi2()/kfpXicPlus.GetNDF(); // chi2_geometry of XicPlus

  // --- chi2_prim of Pion to PV ---
  KFParticle kfpBP_trk1_PV = kfpBP_trk1;
  kfpBP_trk1_PV.SetProductionVertex(PV);
  fVar_XicPlus[30] = kfpBP_trk1_PV.GetChi2()/kfpBP_trk1_PV.GetNDF(); // chi2_topo of pion_0 to PV
  KFParticle kfpBP_trk2_PV = kfpBP_trk2;
  kfpBP_trk2_PV.SetProductionVertex(PV);
  fVar_XicPlus[31] = kfpBP_trk2_PV.GetChi2()/kfpBP_trk2_PV.GetNDF(); // chi2_topo of pion_1 to PV
  // -------------------------------
  
  // --- DCA of Pion to PV ---
  fVar_XicPlus[32] = kfpBP_trk1.GetDistanceFromVertexXY(PV); // DCA of pion_0 coming from XicPlus in x-y plane
  fVar_XicPlus[33] = kfpBP_trk2.GetDistanceFromVertexXY(PV); // DCA of pion_1 coming from XicPlus in x-y plane
  // -------------------------

  Float_t massK0S_Rec, err_massK0S_Rec;
  kfpK0Short.GetMass(massK0S_Rec, err_massK0S_Rec);
  fVar_XicPlus[34] = massK0S_Rec; // mass of Ks0
  Float_t massGamma_Rec, err_massGamma_Rec;
  kfpGamma.GetMass(massGamma_Rec, err_massGamma_Rec);
  fVar_XicPlus[35] = massGamma_Rec; // mass of e+e-

  fVar_XicPlus[36] = kfpPionFromLam.GetDistanceFromParticleXY(kfpProtonFromLam); // DCA of Lam's daughters
  fVar_XicPlus[37] = kfpPionOrKaon.GetDistanceFromParticleXY(kfpLambda_m); // DCA of Xi's daughters (calculated from KF after Lambda mass constraint)
  fVar_XicPlus[38] = kfpXiMinus_m.GetDistanceFromVertexXY(PV); // DCA of Xi to PV in x-y plane (calculated from KF after Xi mass constraint)
  fVar_XicPlus[39] = kfpBP_trk1.GetDistanceFromParticleXY(kfpBP_trk2); // DCA of pi to pi in x-y plane
  fVar_XicPlus[40] = kfpBP_trk1.GetDistanceFromParticle(kfpBP_trk2); // DCA of pi to pi
  fVar_XicPlus[41] = kfpBP_trk1.GetDistanceFromParticleXY(kfpXiMinus_m); // DCA of pi_trk1 to Xi in x-y plane
  fVar_XicPlus[42] = kfpBP_trk2.GetDistanceFromParticleXY(kfpXiMinus_m); // DCA of pi_trk2 to Xi in x-y plane

  // --- CosPointingAngle ---
  Double_t cosPA_XicPlusToPV = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXicPlus, PV);
  fVar_XicPlus[43] = TMath::ACos(cosPA_XicPlusToPV); // pointing angle of XicPlus (pointing back to PV)
  // -------------------------
  Float_t DecayLxy_XicPlus, err_DecayLxy_XicPlus;
  kfpXicPlus_PV.GetDecayLengthXY(DecayLxy_XicPlus, err_DecayLxy_XicPlus);
  fVar_XicPlus[44] = DecayLxy_XicPlus; // decay length of XicPlus in x-y plane
  fVar_XicPlus[45] = kfpXicPlus_PV.GetChi2()/kfpXicPlus_PV.GetNDF(); // chi2_topo of XicPlus to PV
  fVar_XicPlus[46] = l_XicPlus/dl_XicPlus; // l/dl of XicPlus

//  fVar_XicPlus[26] = AliVertexingHFUtils::CosThetaStarFromKF(0, 4132, 211, 3312, kfpXicPlus, kfpBP_XicPlus, kfpXiMinus_XicPlus);
//  fVar_XicPlus[27] = AliVertexingHFUtils::CosThetaStarFromKF(1, 4132, 211, 3312, kfpXicPlus, kfpBP_XicPlus, kfpXiMinus_XicPlus);


//  fVar_XicPlus[33] = kfpLambda_Xi.GetChi2()/kfpLambda_Xi.GetNDF();
//  fVar_XicPlus[34] = kfpXiMinus_XicPlus.GetChi2()/kfpXiMinus_XicPlus.GetNDF();

  Float_t ct_XicPlus=0., err_ct_XicPlus=0.;
  kfpXicPlus_PV.GetLifeTime(ct_XicPlus, err_ct_XicPlus);
  fVar_XicPlus[47] = ct_XicPlus;
//  fVar_XicPlus[38] = casc->DcaV0Daughters(); // DCA_LamDau
//  fVar_XicPlus[40] = kfpBP_trk1.GetDistanceFromParticle(kfpXiMinus_m); // DCA_XicPlusDau_KF

  fVar_XicPlus[48] = TMath::ACos(cosPA_XiToXicPlus); // pointing angle of Xi (pointing back to XicPlus)

  fVar_XicPlus[63] = fpVtx->GetNContributors();

  if (fIsMC) {
    fVar_XicPlus_EvtID = GetMCEventID(); // Event ID for MC
    fVar_XicPlus[64] = lab_XicPlus;
    // === weight ===
    /*
    if (lab_XicPlus>0) {
      Int_t labelPion_trk1 = fabs(trackPiFromXicPlus_trk1->GetLabel());
      AliAODMCParticle* mcPion_trk1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion_trk1));
      AliAODMCParticle* mcXicPlus  = static_cast<AliAODMCParticle*>(mcArray->At(mcPion_trk1->GetMother()));
    }
    */
  }
  if (!fIsMC) {
    fVar_XicPlus_EvtID = GetEventIdAsLong(AODEvent->GetHeader()); // Event ID for Data
    //fVar_XicPlus_QA[76] = AODEvent->GetHeader()->GetEventIdAsLong(); // Event ID for Data
  }

  KFParticle kfpPionOrKaon_Rej;
  kfpPionOrKaon_Rej = AliVertexingHFUtils::CreateKFParticleFromAODtrack(trackPiFromXiOrKaonFromOmega, -321); // -321 kaon- for Xi analysis
  KFParticle kfpCasc_Rej;
  const KFParticle *vCasc_Rej_Ds[2] = {&kfpPionOrKaon_Rej, &kfpLambda_m};
  kfpCasc_Rej.Construct(vCasc_Rej_Ds, 2);
  Float_t massCasc_Rej, err_massCasc_Rej;
  kfpCasc_Rej.GetMass(massCasc_Rej, err_massCasc_Rej);
  fVar_XicPlus[50] = massCasc_Rej;

  fVar_XicPlus[49] = trackPiFromXiOrKaonFromOmega->Pt(); // pt of pion coming from Xi

  //======= Fill QA tree =======
  for (Int_t i=0; i<(62-1); i++) {
    fVar_XicPlus_QA[i] = -9999.;
  }

  if (fIsMC) {
    Double_t PV_Gen[3]={0.};
    AliAODMCHeader *mcHeader = (AliAODMCHeader*)AODEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    mcHeader->GetVertex(PV_Gen);
    fVar_XicPlus_QA[0] = PV_Gen[0]; // PV: X_{MC}
    fVar_XicPlus_QA[1] = PV_Gen[1]; // PV: Y_{MC}
    fVar_XicPlus_QA[2] = PV_Gen[2]; // PV: Z_{MC}
  }

  Double_t pos_PV[3]={0.}, cov_PV[6]={0.};
  fpVtx->GetXYZ(pos_PV);
  fpVtx->GetCovarianceMatrix(cov_PV);
  fVar_XicPlus_QA[3] = pos_PV[0]; // PV: X_{rec}
  fVar_XicPlus_QA[4] = pos_PV[1]; // PV: Y_{rec}
  fVar_XicPlus_QA[5] = pos_PV[2]; // PV: Z_{rec}
  fVar_XicPlus_QA[6] = sqrt(cov_PV[0]); // PV: sigma_X^{rec}
  fVar_XicPlus_QA[7] = sqrt(cov_PV[2]); // PV: sigma_Y^{rec}
  fVar_XicPlus_QA[8] = sqrt(cov_PV[5]); // PV: sigma_Z^{rec}

  //// Recalculate primary vertex without daughter tracks, if requested
  //// IMPORTANT: Own primary vertex must be unset before continue/return, else memory leak
  if (fAnaCuts->GetIsPrimaryWithoutDaughters()) {
    AliAODVertex *PV_woDau = CallPrimaryVertex(casc, trackPiFromXicPlus_trk1, trackPiFromXicPlus_trk2, AODEvent);
    if (!PV_woDau) PV_woDau=fpVtx;
    Double_t pos_recalPV[3]={0.}, cov_recalPV[6]={0.};
    PV_woDau->GetXYZ(pos_recalPV);
    PV_woDau->GetCovarianceMatrix(cov_recalPV);
    //PV_woDau->GetSigmaXYZ(sigma);
    fVar_XicPlus_QA[9]  = pos_recalPV[0]; // recal_PV: X_{rec}
    fVar_XicPlus_QA[10] = pos_recalPV[1]; // recal_PV: Y_{rec}
    fVar_XicPlus_QA[11] = pos_recalPV[2]; // recal_PV: Z_{rec}
    fHPrimVtx_woDau_x->Fill(pos_recalPV[0]);
    fHPrimVtx_woDau_y->Fill(pos_recalPV[1]);
    fHPrimVtx_woDau_z->Fill(pos_recalPV[2]);
    if (cov_recalPV[0]<0 || cov_recalPV[2]<0 || cov_recalPV[5]<0) {
      fHPrimVtx_woDau_err_x->Fill(9999.);
      fHPrimVtx_woDau_err_y->Fill(9999.);
      fHPrimVtx_woDau_err_z->Fill(9999.);
    }
    if (cov_recalPV[0]>=0 && cov_recalPV[2]>=0 && cov_recalPV[5]>=0) {
      fVar_XicPlus_QA[12] = sqrt(cov_recalPV[0]); // recal_PV: sigma_X^{rec}
      fVar_XicPlus_QA[13] = sqrt(cov_recalPV[2]); // recal_PV: sigma_Y^{rec}
      fVar_XicPlus_QA[14] = sqrt(cov_recalPV[5]); // recal_PV: sigma_Z^{rec}
      fHPrimVtx_woDau_err_x->Fill(sqrt(cov_recalPV[0]));
      fHPrimVtx_woDau_err_y->Fill(sqrt(cov_recalPV[2]));
      fHPrimVtx_woDau_err_z->Fill(sqrt(cov_recalPV[5]));
    }
    fVar_XicPlus_QA[15] = fpVtx->GetNContributors(); // PV: NContributors
    fVar_XicPlus_QA[16] = PV_woDau->GetNContributors(); // recal_PV: NContributors
    //cout << "PV (NContributors): " << fpVtx->GetNContributors() << endl;
    //cout << "PV (NDaughters): " << fpVtx->GetNDaughters() << endl;
    //cout << "PV (CountRealContributors): " << fpVtx->CountRealContributors() << endl;
    //cout << "recal_PV (NContributors): " << PV_woDau->GetNContributors() << endl;
    //cout << "recal_PV (NDaughters): " << PV_woDau->GetNDaughters() << endl;
    //cout << "recal_PV (CountRealContributors): " << PV_woDau->CountRealContributors() << endl;

    // set recalculated primary vertex to KF
    KFPVertex pVertex_recal;
    pVertex_recal.SetXYZ((Float_t)pos_recalPV[0], (Float_t)pos_recalPV[1], (Float_t)pos_recalPV[2]);
    Float_t covF[6];
    for (Int_t i=0; i<6; i++) { covF[i] = (Float_t)cov_recalPV[i]; }
    pVertex_recal.SetCovarianceMatrix(covF);
    pVertex_recal.SetChi2(PV_woDau->GetChi2());
    pVertex_recal.SetNDF(PV_woDau->GetNDF());
    pVertex_recal.SetNContributors(PV_woDau->GetNContributors());
    KFParticle PV_recal(pVertex_recal);

    KFParticle kfpXicPlus_recalPV = kfpXicPlus;
    kfpXicPlus_recalPV.SetProductionVertex(PV_recal);
    KFParticle kfpXicPlus_recalPV_DecayVtx = kfpXicPlus_recalPV;
    kfpXicPlus_recalPV_DecayVtx.TransportToDecayVertex();
    fVar_XicPlus_QA[34] = kfpXicPlus_recalPV_DecayVtx.GetX(); // SV: X_{rec} (w/ topo. constraint recalPV)
    fVar_XicPlus_QA[35] = kfpXicPlus_recalPV_DecayVtx.GetY(); // SV: Y_{rec} (w/ topo. constraint recalPV)
    fVar_XicPlus_QA[36] = kfpXicPlus_recalPV_DecayVtx.GetZ(); // SV: Z_{rec} (w/ topo. constraint recalPV)
    fVar_XicPlus_QA[37] = kfpXicPlus_recalPV_DecayVtx.GetErrX(); // SV: sigma_X^{rec} (w/ topo. constraint recalPV)
    fVar_XicPlus_QA[38] = kfpXicPlus_recalPV_DecayVtx.GetErrY(); // SV: sigma_Y^{rec} (w/ topo. constraint recalPV)
    fVar_XicPlus_QA[39] = kfpXicPlus_recalPV_DecayVtx.GetErrZ(); // SV: sigma_Z^{rec} (w/ topo. constraint recalPV)
    fVar_XicPlus_QA[40] = kfpXicPlus_recalPV.GetPt(); // SV: pt of Xic+ (w/ topo. constraint recalPV)

    // === recalPV_KF =========
    // --- CosPointingAngle ---
    Double_t cosPA_XicPlusToRecalPV = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXicPlus, PV_recal);
    fVar_XicPlus[55] = TMath::ACos(cosPA_XicPlusToRecalPV); // pointing angle of XicPlus (pointing back to recalPV)
    // --- chi2_topo ---
    fVar_XicPlus[56] = kfpXicPlus_recalPV.GetChi2()/kfpXicPlus_recalPV.GetNDF(); // chi2_topo of XicPlus to recalPV
    // ========================

  // === Recalculate PV after removing Xic+ daughters and adding Xic+ (Refit) ===
  KFParticle recalPV_KF_Refit = PV_KF_Refit;
  if (trackPiFromXicPlus_trk1->GetUsedForPrimVtxFit()) kfpBP_trk1.SubtractFromVertex(recalPV_KF_Refit);
  if (trackPiFromXicPlus_trk2->GetUsedForPrimVtxFit()) kfpBP_trk2.SubtractFromVertex(recalPV_KF_Refit);
  if (trackPiFromXiOrKaonFromOmega->GetUsedForPrimVtxFit()) kfpPionOrKaon.SubtractFromVertex(recalPV_KF_Refit);
  if (trkProton->GetUsedForPrimVtxFit()) kfpProtonFromLam.SubtractFromVertex(recalPV_KF_Refit);
  if (trkPion->GetUsedForPrimVtxFit()) kfpPionFromLam.SubtractFromVertex(recalPV_KF_Refit);
  // --- without adding mother particle ---
  KFParticle recalPV_KF_Refit_woAddMother = recalPV_KF_Refit;
  Double_t cosPA_XicPlusToRecalPVKF_Refit_woAddMother = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXicPlus, recalPV_KF_Refit_woAddMother);
  fVar_XicPlus[61] = TMath::ACos(cosPA_XicPlusToRecalPVKF_Refit_woAddMother);
  KFParticle kfpXicPlus_recalPVKF_Refit_woAddMother = kfpXicPlus;
  kfpXicPlus_recalPVKF_Refit_woAddMother.SetProductionVertex(recalPV_KF_Refit_woAddMother);
  fVar_XicPlus[62] = kfpXicPlus_recalPVKF_Refit_woAddMother.GetChi2()/kfpXicPlus_recalPVKF_Refit_woAddMother.GetNDF();
  // --------------------------------------

  recalPV_KF_Refit.AddDaughter(kfpXicPlus);

  // --- CosPointingAngle ---
  Double_t cosPA_XicPlusToRecalPVKF_Refit = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXicPlus, recalPV_KF_Refit);
  fVar_XicPlus[51] = TMath::ACos(cosPA_XicPlusToRecalPVKF_Refit); // pointing angle of XicPlus (pointing back to recalPV_KF_Refit)
  // --- chi2_topo ---
  KFParticle kfpXicPlus_recalPVKF_Refit = kfpXicPlus;
  kfpXicPlus_recalPVKF_Refit.SetProductionVertex(recalPV_KF_Refit);
  fVar_XicPlus[52] = kfpXicPlus_recalPVKF_Refit.GetChi2()/kfpXicPlus_recalPVKF_Refit.GetNDF(); // chi2_topo of XicPlus to recalPV_KF_Refit
  // --- CosPointingAngle (TransportToDecayVertex) ---
  /*
  KFParticle kfpXiMinus_m_FullyFitted = kfpXiMinus_m;
  kfpXiMinus_m_FullyFitted.SetProductionVertex(kfpXicPlus_recalPVKF_Refit);
  KFParticle kfpBP_trk1_FullyFitted = kfpBP_trk1;
  kfpBP_trk1_FullyFitted.SetProductionVertex(kfpXicPlus_recalPVKF_Refit);
  KFParticle kfpBP_trk2_FullyFitted = kfpBP_trk2;
  kfpBP_trk2_FullyFitted.SetProductionVertex(kfpXicPlus_recalPVKF_Refit);
  cout << "-----------------------------" << endl;
  cout << "kfpXiMinus_m: " << kfpXiMinus_m.GetX() << endl;
  cout << "kfpXiMinus_m_FullyFitted: " <<  kfpXiMinus_m_FullyFitted.GetX() << endl;
  cout << "kfpBP_trk1: " << kfpBP_trk1.GetX() << endl;
  cout << "kfpBP_trk1_FullyFitted: " << kfpBP_trk1_FullyFitted.GetX() << endl;
  cout << "kfpBP_trk2: " << kfpBP_trk2.GetX() << endl;
  cout << "kfpBP_trk2_FullyFitted: " << kfpBP_trk2_FullyFitted.GetX() << endl;
  cout << "kfpXicPlus_recalPVKF_Refit: " << kfpXicPlus_recalPVKF_Refit.GetX() << endl;
  */
  /*
  cout << "kfpXicPlus_recalPVKF_Refit (TransportToDecayVertex): " << kfpXicPlus_recalPVKF_Refit.GetX() << endl;
  cout << "kfpXicPlus: " << kfpXicPlus.GetX() << endl;
  KFParticle kfpBP_trk1_FullyFitted_XicPlusDecayVertex = kfpBP_trk1;
  kfpBP_trk1_FullyFitted_XicPlusDecayVertex.SetProductionVertex(kfpXicPlus_recalPVKF_Refit);
  KFParticle kfpBP_trk2_FullyFitted_XicPlusDecayVertex = kfpBP_trk2;
  kfpBP_trk2_FullyFitted_XicPlusDecayVertex.SetProductionVertex(kfpXicPlus_recalPVKF_Refit);
  KFParticle kfpXiMinus_m_FullyFitted_XicPlusDecayVertex = kfpXiMinus_m;
  kfpXiMinus_m_FullyFitted_XicPlusDecayVertex.SetProductionVertex(kfpXicPlus_recalPVKF_Refit);
  cout << "kfpBP_trk1_FullyFitted_XicPlusDecayVertex: " << kfpBP_trk1_FullyFitted_XicPlusDecayVertex.GetX() << endl;
  cout << "kfpBP_trk2_FullyFitted_XicPlusDecayVertex: " << kfpBP_trk2_FullyFitted_XicPlusDecayVertex.GetX() << endl;
  cout << "kfpXiMinus_m_FullyFitted_XicPlusDecayVertex: " << kfpXiMinus_m_FullyFitted_XicPlusDecayVertex.GetX() << endl;
  */
  // ============================================================================

  // === Recalculate PV after removing Xic+ daughters and adding Xic+ (w/o Refit) ===
  KFParticle recalPV_KF = PV;
  if (trackPiFromXicPlus_trk1->GetUsedForPrimVtxFit()) kfpBP_trk1.SubtractFromVertex(recalPV_KF);
  if (trackPiFromXicPlus_trk2->GetUsedForPrimVtxFit()) kfpBP_trk2.SubtractFromVertex(recalPV_KF);
  if (trackPiFromXiOrKaonFromOmega->GetUsedForPrimVtxFit()) kfpPionOrKaon.SubtractFromVertex(recalPV_KF);
  if (trkProton->GetUsedForPrimVtxFit()) kfpProtonFromLam.SubtractFromVertex(recalPV_KF);
  if (trkPion->GetUsedForPrimVtxFit()) kfpPionFromLam.SubtractFromVertex(recalPV_KF);
  recalPV_KF.AddDaughter(kfpXicPlus);

  // --- CosPointingAngle ---
  Double_t cosPA_XicPlusToRecalPVKF = AliVertexingHFUtils::CosPointingAngleFromKF(kfpXicPlus, recalPV_KF);
  fVar_XicPlus[53] = TMath::ACos(cosPA_XicPlusToRecalPVKF); // pointing angle of XicPlus (pointing back to recalPV_KF)
  // --- chi2_topo ---
  KFParticle kfpXicPlus_recalPVKF = kfpXicPlus;
  kfpXicPlus_recalPVKF.SetProductionVertex(recalPV_KF);
  fVar_XicPlus[54] = kfpXicPlus_recalPVKF.GetChi2()/kfpXicPlus_recalPVKF.GetNDF(); // chi2_topo of XicPlus to recalPV_KF
  // ============================================================================

  // --- pointing angle (X-Y) ---
  Double_t cosPAXY_XicPlusToRecalPVKF_Refit = AliVertexingHFUtils::CosPointingAngleXYFromKF(kfpXicPlus, recalPV_KF_Refit);
  fVar_XicPlus[58] = TMath::ACos(cosPAXY_XicPlusToRecalPVKF_Refit); // pointing angle (X-Y) of XicPlus (pointing back to recalPV_KF_Refit)
  Double_t cosPAXY_XicPlusToRecalPVKF = AliVertexingHFUtils::CosPointingAngleXYFromKF(kfpXicPlus, recalPV_KF);
  fVar_XicPlus[59] = TMath::ACos(cosPAXY_XicPlusToRecalPVKF); // pointing angle (X-Y) of XicPlus (pointing back to recalPV_KF)
  Double_t cosPAXY_XicPlusToRecalPV = AliVertexingHFUtils::CosPointingAngleXYFromKF(kfpXicPlus, PV_recal);
  fVar_XicPlus[60] = TMath::ACos(cosPAXY_XicPlusToRecalPV); // pointing angle (X-Y) of XicPlus (pointing back to recalPV)
  //-----------------------------

  fVar_XicPlus_QA[42] = PV_KF_Refit.GetX(); // PV_KF_Refit: X_{rec}
  fVar_XicPlus_QA[43] = PV_KF_Refit.GetY(); // PV_KF_Refit: Y_{rec}
  fVar_XicPlus_QA[44] = PV_KF_Refit.GetZ(); // PV_KF_Refit: Z_{rec}
  fVar_XicPlus_QA[45] = PV_KF_Refit.GetErrX(); // PV_KF_Refit: sigma_X^{rec}
  fVar_XicPlus_QA[46] = PV_KF_Refit.GetErrY(); // PV_KF_Refit: sigma_Y^{rec}
  fVar_XicPlus_QA[47] = PV_KF_Refit.GetErrZ(); // PV_KF_Refit: sigma_Z^{rec}

  fVar_XicPlus_QA[48] = recalPV_KF_Refit.GetX(); // recalPV_KF_Refit: X_{rec}
  fVar_XicPlus_QA[49] = recalPV_KF_Refit.GetY(); // recalPV_KF_Refit: Y_{rec}
  fVar_XicPlus_QA[50] = recalPV_KF_Refit.GetZ(); // recalPV_KF_Refit: Z_{rec}
  fVar_XicPlus_QA[51] = recalPV_KF_Refit.GetErrX(); // recalPV_KF_Refit: sigma_X^{rec}
  fVar_XicPlus_QA[52] = recalPV_KF_Refit.GetErrY(); // recalPV_KF_Refit: sigma_Y^{rec}
  fVar_XicPlus_QA[53] = recalPV_KF_Refit.GetErrZ(); // recalPV_KF_Refit: sigma_Z^{rec}

  fVar_XicPlus_QA[54] = recalPV_KF.GetX(); // recalPV_KF: X_{rec}
  fVar_XicPlus_QA[55] = recalPV_KF.GetY(); // recalPV_KF: Y_{rec}
  fVar_XicPlus_QA[56] = recalPV_KF.GetZ(); // recalPV_KF: Z_{rec}
  fVar_XicPlus_QA[57] = recalPV_KF.GetErrX(); // recalPV_KF: sigma_X^{rec}
  fVar_XicPlus_QA[58] = recalPV_KF.GetErrY(); // recalPV_KF: sigma_Y^{rec}
  fVar_XicPlus_QA[59] = recalPV_KF.GetErrZ(); // recalPV_KF: sigma_Z^{rec}

  fHPrimVtx_PV_KF_Refit_Minus_PVrec_x->Fill(PV_KF_Refit.GetX()-pos_PV[0]);
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_y->Fill(PV_KF_Refit.GetY()-pos_PV[1]);
  fHPrimVtx_PV_KF_Refit_Minus_PVrec_z->Fill(PV_KF_Refit.GetZ()-pos_PV[2]);
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_x->Fill(recalPV_KF.GetX()-PV_KF_Refit.GetX());
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_y->Fill(recalPV_KF.GetY()-PV_KF_Refit.GetY());
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_z->Fill(recalPV_KF.GetZ()-PV_KF_Refit.GetZ());
  fHPrimVtx_recalPV_Minus_PVrec_x->Fill(fVar_XicPlus_QA[9]-pos_PV[0]);
  fHPrimVtx_recalPV_Minus_PVrec_y->Fill(fVar_XicPlus_QA[10]-pos_PV[1]);
  fHPrimVtx_recalPV_Minus_PVrec_z->Fill(fVar_XicPlus_QA[11]-pos_PV[2]);
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_x->Fill(recalPV_KF_Refit.GetX()-PV_KF_Refit.GetX());
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_y->Fill(recalPV_KF_Refit.GetY()-PV_KF_Refit.GetY());
  fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_z->Fill(recalPV_KF_Refit.GetZ()-PV_KF_Refit.GetZ());

  fHPrimVtx_PV_KF_Refit_Minus_PVgen_x->Fill(fVar_XicPlus_QA[76]-fVar_XicPlus_QA[0]);
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_y->Fill(fVar_XicPlus_QA[77]-fVar_XicPlus_QA[1]);
  fHPrimVtx_PV_KF_Refit_Minus_PVgen_z->Fill(fVar_XicPlus_QA[78]-fVar_XicPlus_QA[2]);
  fHPrimVtx_recalPV_Minus_PVgen_x->Fill(fVar_XicPlus_QA[9]-fVar_XicPlus_QA[0]);
  fHPrimVtx_recalPV_Minus_PVgen_y->Fill(fVar_XicPlus_QA[10]-fVar_XicPlus_QA[1]);
  fHPrimVtx_recalPV_Minus_PVgen_z->Fill(fVar_XicPlus_QA[11]-fVar_XicPlus_QA[2]);
  fHPrimVtx_PV_PULL_x->Fill((pos_PV[0]-fVar_XicPlus_QA[0])/fVar_XicPlus_QA[6]);
  fHPrimVtx_PV_PULL_y->Fill((pos_PV[1]-fVar_XicPlus_QA[1])/fVar_XicPlus_QA[7]);
  fHPrimVtx_PV_PULL_z->Fill((pos_PV[2]-fVar_XicPlus_QA[2])/fVar_XicPlus_QA[8]);
  fHPrimVtx_PV_KF_PULL_x->Fill((PV.GetX()-fVar_XicPlus_QA[0])/PV.GetErrX());
  fHPrimVtx_PV_KF_PULL_y->Fill((PV.GetY()-fVar_XicPlus_QA[1])/PV.GetErrY());
  fHPrimVtx_PV_KF_PULL_z->Fill((PV.GetZ()-fVar_XicPlus_QA[2])/PV.GetErrZ());
  fHPrimVtx_PV_KF_Refit_PULL_x->Fill((PV_KF_Refit.GetX()-fVar_XicPlus_QA[0])/PV_KF_Refit.GetErrX());
  fHPrimVtx_PV_KF_Refit_PULL_y->Fill((PV_KF_Refit.GetY()-fVar_XicPlus_QA[1])/PV_KF_Refit.GetErrY());
  fHPrimVtx_PV_KF_Refit_PULL_z->Fill((PV_KF_Refit.GetZ()-fVar_XicPlus_QA[2])/PV_KF_Refit.GetErrZ());
  fHPrimVtx_recalPV_PULL_x->Fill((fVar_XicPlus_QA[9]-fVar_XicPlus_QA[0])/fVar_XicPlus_QA[12]);
  fHPrimVtx_recalPV_PULL_y->Fill((fVar_XicPlus_QA[10]-fVar_XicPlus_QA[1])/fVar_XicPlus_QA[13]);
  fHPrimVtx_recalPV_PULL_z->Fill((fVar_XicPlus_QA[11]-fVar_XicPlus_QA[2])/fVar_XicPlus_QA[14]);

  fHPrimVtx_recalPV_KF_Minus_PV_x->Fill(recalPV_KF.GetX()-PV.GetX());
  fHPrimVtx_recalPV_KF_Minus_PV_y->Fill(recalPV_KF.GetY()-PV.GetY());
  fHPrimVtx_recalPV_KF_Minus_PV_z->Fill(recalPV_KF.GetZ()-PV.GetZ());

  if (fIsMC) {
    fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_x->Fill(recalPV_KF_Refit.GetX()-fVar_XicPlus_QA[0]);
    fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_y->Fill(recalPV_KF_Refit.GetY()-fVar_XicPlus_QA[1]);
    fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_z->Fill(recalPV_KF_Refit.GetZ()-fVar_XicPlus_QA[2]);
    fHPrimVtx_recalPV_KF_Refit_PULL_x->Fill((recalPV_KF_Refit.GetX()-fVar_XicPlus_QA[0])/recalPV_KF_Refit.GetErrX());
    fHPrimVtx_recalPV_KF_Refit_PULL_y->Fill((recalPV_KF_Refit.GetY()-fVar_XicPlus_QA[1])/recalPV_KF_Refit.GetErrY());
    fHPrimVtx_recalPV_KF_Refit_PULL_z->Fill((recalPV_KF_Refit.GetZ()-fVar_XicPlus_QA[2])/recalPV_KF_Refit.GetErrZ());
    fHPrimVtx_recalPV_KF_Minus_PVgen_x->Fill(recalPV_KF.GetX()-fVar_XicPlus_QA[0]);
    fHPrimVtx_recalPV_KF_Minus_PVgen_y->Fill(recalPV_KF.GetY()-fVar_XicPlus_QA[1]);
    fHPrimVtx_recalPV_KF_Minus_PVgen_z->Fill(recalPV_KF.GetZ()-fVar_XicPlus_QA[2]);
    fHPrimVtx_recalPV_KF_PULL_x->Fill((recalPV_KF.GetX()-fVar_XicPlus_QA[0])/recalPV_KF.GetErrX());
    fHPrimVtx_recalPV_KF_PULL_y->Fill((recalPV_KF.GetY()-fVar_XicPlus_QA[1])/recalPV_KF.GetErrY());
    fHPrimVtx_recalPV_KF_PULL_z->Fill((recalPV_KF.GetZ()-fVar_XicPlus_QA[2])/recalPV_KF.GetErrZ());
  }
  }

  // --- pointing angle (X-Y) ---
  Double_t cosPAXY_XicPlusToPV = AliVertexingHFUtils::CosPointingAngleXYFromKF(kfpXicPlus, PV);
  fVar_XicPlus[57] = TMath::ACos(cosPAXY_XicPlusToPV); // pointing angle (X-Y) of XicPlus (pointing back to PV)

  fVar_XicPlus_QA[20] = kfpXicPlus.GetX(); // SV: X_{rec}
  fVar_XicPlus_QA[21] = kfpXicPlus.GetY(); // SV: Y_{rec}
  fVar_XicPlus_QA[22] = kfpXicPlus.GetZ(); // SV: Z_{rec}
  fVar_XicPlus_QA[23] = kfpXicPlus.GetErrX(); // SV: sigma_X^{rec}
  fVar_XicPlus_QA[24] = kfpXicPlus.GetErrY(); // SV: sigma_Y^{rec}
  fVar_XicPlus_QA[25] = kfpXicPlus.GetErrZ(); // SV: sigma_Z^{rec}
  KFParticle kfpXicPlus_PV_DecayVtx = kfpXicPlus_PV;
  kfpXicPlus_PV_DecayVtx.TransportToDecayVertex();
  fVar_XicPlus_QA[26] = kfpXicPlus_PV_DecayVtx.GetX(); // SV: X_{rec} (w/ topo. constraint)
  fVar_XicPlus_QA[27] = kfpXicPlus_PV_DecayVtx.GetY(); // SV: Y_{rec} (w/ topo. constraint)
  fVar_XicPlus_QA[28] = kfpXicPlus_PV_DecayVtx.GetZ(); // SV: Z_{rec} (w/ topo. constraint)
  fVar_XicPlus_QA[29] = kfpXicPlus_PV_DecayVtx.GetErrX(); // SV: sigma_X^{rec} (w/ topo. constraint)
  fVar_XicPlus_QA[30] = kfpXicPlus_PV_DecayVtx.GetErrY(); // SV: sigma_Y^{rec} (w/ topo. constraint)
  fVar_XicPlus_QA[31] = kfpXicPlus_PV_DecayVtx.GetErrZ(); // SV: sigma_Z^{rec} (w/ topo. constraint)
  fVar_XicPlus_QA[32] = kfpXicPlus.GetPt(); // SV: pt of Xic+
  fVar_XicPlus_QA[33] = kfpXicPlus_PV.GetPt(); // SV: pt of Xic+ (w/ topo. constraint)

  fVar_XicPlus_QA[41] = fpVtx->CountRealContributors(); // PV: count daughter primary tracks

  if (fIsMC) {
  // SV
  Int_t labelPiFromXicPlus_trk1 = fabs(trackPiFromXicPlus_trk1->GetLabel());
  AliAODMCParticle* mcPiFromXicPlus_trk1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPiFromXicPlus_trk1));
  fVar_XicPlus_QA[17] = mcPiFromXicPlus_trk1->Xv(); // SV: X_{MC}
  fVar_XicPlus_QA[18] = mcPiFromXicPlus_trk1->Yv(); // SV: Y_{MC}
  fVar_XicPlus_QA[19] = mcPiFromXicPlus_trk1->Zv(); // SV: Z_{MC}

  fVar_XicPlus_QA[60] = lab_XicPlus;
  fHPrimVtx_PVrec_Minus_PVgen_x->Fill(pos_PV[0]-fVar_XicPlus_QA[0]);
  fHPrimVtx_PVrec_Minus_PVgen_y->Fill(pos_PV[1]-fVar_XicPlus_QA[1]);
  fHPrimVtx_PVrec_Minus_PVgen_z->Fill(pos_PV[2]-fVar_XicPlus_QA[2]);
  }

  // pt(XicPlus)>=1
  //if (fVar_XicPlus[26]>0.9999) fTree_XicPlus->Fill();
  fTree_XicPlus->Fill();

  if (fWriteXicPlusQATree) fTree_XicPlus_QA->Fill();
  //============================

  if ( fabs(fVar_XicPlus[28]-PDGmassXicPlus) < 0.03 ) fCount_NumOfCandidatePerEvent_In3Sigma++;

//  fVar_XicPlus[10] = casc->DcaXiDaughters(); // DCA_XiDau
//  fVar_XicPlus[32] = AliVertexingHFUtils::CosThetaStarFromKF(0, 4132, 211, 3312, kfpXicPlus, kfpBP_trk1, kfpXiMinus_m);
//  fVar_XicPlus[33] = AliVertexingHFUtils::CosThetaStarFromKF(1, 4132, 211, 3312, kfpXicPlus, kfpBP_trk1, kfpXiMinus_m);

//  fVar_XicPlus[3]  = TMath::Prob(kfpLambda.GetChi2(), kfpLambda.GetNDF());
//  fVar_XicPlus[7]  = kfpLambda.GetPz();
//  fVar_XicPlus[8]  = kfpLambda.GetX();
//  fVar_XicPlus[9]  = kfpLambda.GetY();
//  fVar_XicPlus[10] = kfpLambda.GetZ();
//  fVar_XicPlus[8] = kfpLambda_pv.GetChi2()/kfpLambda_pv.GetNDF();
//  fVar_XicPlus[9] = TMath::Prob(kfpLambda_pv.GetChi2(), kfpLambda_pv.GetNDF());
//  fVar_XicPlus[15] = alpha_FirstDaugPos;
//  fVar_XicPlus[16] = v0->PtArmV0();

//  fVar_XicPlus[17] = trackPiFromXi->Charge();
//  fVar_XicPlus[20] = TMath::Prob(kfpXiMinus.GetChi2(), kfpXiMinus.GetNDF());
//  fVar_XicPlus[21] = kfpXiMinus_pv.GetChi2()/kfpXiMinus_pv.GetNDF();
//  fVar_XicPlus[22] = TMath::Prob(kfpXiMinus_pv.GetChi2(), kfpXiMinus_pv.GetNDF());
//  Float_t mass_Xi_PV, err_mass_Xi_PV;
//  kfpXiMinus_pv.GetMass(mass_Xi_PV, err_mass_Xi_PV);
//  fVar_XicPlus[23] = mass_Xi_PV;
//  if ( TMath::Abs(kfpXiMinus_pv.GetE())>TMath::Abs(kfpXiMinus_pv.GetPz()) ) {
//    fVar_XicPlus[24] = kfpXiMinus_pv.GetRapidity();
//  }
//  fVar_XicPlus[25] = kfpXiMinus_pv.GetPt();
//  fVar_XicPlus[26] = kfpXiMinus_pv.GetPz();
//  fVar_XicPlus[27] = kfpXiMinus_pv.GetX();
//  fVar_XicPlus[28] = kfpXiMinus_pv.GetY();
//  fVar_XicPlus[29] = kfpXiMinus_pv.GetZ();
//  fVar_XicPlus[33] = kfpXiMinus.GetPz();
//  fVar_XicPlus[34] = kfpXiMinus.GetX();
//  fVar_XicPlus[35] = kfpXiMinus.GetY();
//  fVar_XicPlus[36] = kfpXiMinus.GetZ();

//  fVar_XicPlus[39] = trackPiFromXicPlus_trk1->Charge();
//  fVar_XicPlus[41] = trackPiFromXicPlus_trk1->Pz();

//  fVar_XicPlus[20] = kfpXicPlus.GetChi2()/kfpXicPlus.GetNDF();
//  fVar_XicPlus[45] = TMath::Prob(kfpXicPlus.GetChi2(), kfpXicPlus.GetNDF());
//  fVar_XicPlus[21] = kfpXicPlus_PV.GetChi2()/kfpXicPlus_PV.GetNDF();
//  fVar_XicPlus[47] = TMath::Prob(kfpXicPlus_PV.GetChi2(), kfpXicPlus_PV.GetNDF());
//  fVar_XicPlus[51] = kfpXicPlus_PV.GetPz();
//  fVar_XicPlus[52] = kfpXicPlus_PV.GetX();
//  fVar_XicPlus[53] = kfpXicPlus_PV.GetY();
//  fVar_XicPlus[54] = kfpXicPlus_PV.GetZ();
//  Float_t mass_Xic, err_mass_Xic;
//  kfpXicPlus.GetMass(mass_Xic, err_mass_Xic);
//  fVar_XicPlus[55] = mass_Xic;
//  fVar_XicPlus[56] = kfpXicPlus.GetRapidity();
//  fVar_XicPlus[57] = kfpXicPlus.GetPt();
//  fVar_XicPlus[58] = kfpXicPlus.GetPz();
//  fVar_XicPlus[59] = kfpXicPlus.GetX();
//  fVar_XicPlus[60] = kfpXicPlus.GetY();
//  fVar_XicPlus[61] = kfpXicPlus.GetZ();

//  fVar_XicPlus[27] = cosPA_XiToXic;


//  Float_t CT_Lam, err_CT_Lam;
//  Float_t CT_Xi, err_CT_Xi;
//  Float_t CT_XicPlus, err_CT_XicPlus;
//  kfpLambda_Xi.GetLifeTime(CT_Lam, err_CT_Lam);
//  kfpXiMinus_pv.GetLifeTime(CT_Xi, err_CT_Xi);
//  kfpXicPlus_PV.GetLifeTime(CT_XicPlus, err_CT_XicPlus);
//  fVar_XicPlus[38] = CT_Lam;
//  fVar_XicPlus[39] = CT_Xi;
//  fVar_XicPlus[40] = CT_XicPlus;

//  if (fIsMC) {
//    fVar_XicPlus[78] = lab_XicPlus;
//  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillQATreeXicPlusFromCasc_woMassConstForLamAndXi(KFParticle kfpLambda, KFParticle kfpXiMinus_woMassConstForLamAndXi, KFParticle kfpXicPlus_woMassConstForLamAndXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent)
{
  for (Int_t i=0;i<16;i++) {
    fVar_XicPlus_QA_woMassConstForLamAndXi[i] = -9999.;
  }

  // chi2topo_LamToXi
  KFParticle kfpLambda_Xi = kfpLambda;
  kfpLambda_Xi.SetProductionVertex(kfpXiMinus_woMassConstForLamAndXi);
  fVar_XicPlus_QA_woMassConstForLamAndXi[0] = kfpLambda_Xi.GetChi2()/kfpLambda_Xi.GetNDF();

  // chi2topo_XiToXicPlus
  KFParticle kfpXi_XicPlus = kfpXiMinus_woMassConstForLamAndXi;
  kfpXi_XicPlus.SetProductionVertex(kfpXicPlus_woMassConstForLamAndXi);
  fVar_XicPlus_QA_woMassConstForLamAndXi[1] = kfpXi_XicPlus.GetChi2()/kfpXi_XicPlus.GetNDF();

  // chi2topo_XicPlusToPV
  KFParticle kfpXicPlus_PV = kfpXicPlus_woMassConstForLamAndXi;
  kfpXicPlus_PV.SetProductionVertex(PV);
  fVar_XicPlus_QA_woMassConstForLamAndXi[2] = kfpXicPlus_PV.GetChi2()/kfpXicPlus_PV.GetNDF();

  // pt_XicPlus
  fVar_XicPlus_QA_woMassConstForLamAndXi[3] = kfpXicPlus_woMassConstForLamAndXi.GetPt();

  // mass_XicPlus
  Float_t massXicPlus, err_massXicPlus;
  kfpXicPlus_woMassConstForLamAndXi.GetMass(massXicPlus, err_massXicPlus);
  fVar_XicPlus_QA_woMassConstForLamAndXi[4] = massXicPlus;
  
  // SV_X_rec
  fVar_XicPlus_QA_woMassConstForLamAndXi[5] = kfpXicPlus_woMassConstForLamAndXi.GetX();
  // SV_Y_rec
  fVar_XicPlus_QA_woMassConstForLamAndXi[6] = kfpXicPlus_woMassConstForLamAndXi.GetY();
  // SV_Z_rec
  fVar_XicPlus_QA_woMassConstForLamAndXi[7] = kfpXicPlus_woMassConstForLamAndXi.GetZ();
  // SV_sigma_X_rec
  fVar_XicPlus_QA_woMassConstForLamAndXi[8] = kfpXicPlus_woMassConstForLamAndXi.GetErrX();
  // SV_sigma_Y_rec
  fVar_XicPlus_QA_woMassConstForLamAndXi[9] = kfpXicPlus_woMassConstForLamAndXi.GetErrY();
  // SV_sigma_Z_rec
  fVar_XicPlus_QA_woMassConstForLamAndXi[10] = kfpXicPlus_woMassConstForLamAndXi.GetErrZ();

  if (fIsMC) {
    // SV_X_MC
    Int_t labelPiFromXicPlus_HighPt = fabs(trackPiFromXicPlus_HighPt->GetLabel());
    AliAODMCParticle* mcPiFromXicPlus_HighPt = static_cast<AliAODMCParticle*>(mcArray->At(labelPiFromXicPlus_HighPt));
    fVar_XicPlus_QA_woMassConstForLamAndXi[11] = mcPiFromXicPlus_HighPt->Xv();
    // SV_Y_MC
    fVar_XicPlus_QA_woMassConstForLamAndXi[12] = mcPiFromXicPlus_HighPt->Yv();
    // SV_Z_MC
    fVar_XicPlus_QA_woMassConstForLamAndXi[13] = mcPiFromXicPlus_HighPt->Zv();

    AliAODMCParticle* mcXicPlusCand = static_cast<AliAODMCParticle*>(mcArray->At(mcPiFromXicPlus_HighPt->GetMother()));
    if (abs(mcXicPlusCand->GetPdgCode()) == 4232) fVar_XicPlus_QA_woMassConstForLamAndXi[14] = mcXicPlusCand->Pt();
    if (abs(mcXicPlusCand->GetPdgCode()) == 3324) {
      AliAODMCParticle* mcXicPlusCand_DecayToResonance = static_cast<AliAODMCParticle*>(mcArray->At(mcXicPlusCand->GetMother()));
      if (abs(mcXicPlusCand_DecayToResonance->GetPdgCode()) == 4232) fVar_XicPlus_QA_woMassConstForLamAndXi[14] = mcXicPlusCand_DecayToResonance->Pt();
    }
  }

  // Source_XicPlus
  fVar_XicPlus_QA_woMassConstForLamAndXi[15] = lab_XicPlus;

  // event_ID
  if (fIsMC) fVar_XicPlus_EvtID = GetMCEventID();
  if (!fIsMC) fVar_XicPlus_EvtID = GetEventIdAsLong(AODEvent->GetHeader());

  fTree_XicPlus_QA_woMassConstForLamAndXi->Fill();

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillQATreeXicPlusFromCasc_woMassConstForLam_wMassConstForXi(KFParticle kfpLambda, KFParticle kfpXiMinus_woMassConstForLam_wMassConstForXi, KFParticle kfpXicPlus_woMassConstForLam_wMassConstForXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent)
{
  for (Int_t i=0;i<15;i++) {
    fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[i] = -9999.;
  }

  // chi2topo_LamToXi
  KFParticle kfpLambda_Xi = kfpLambda;
  kfpLambda_Xi.SetProductionVertex(kfpXiMinus_woMassConstForLam_wMassConstForXi);
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[0] = kfpLambda_Xi.GetChi2()/kfpLambda_Xi.GetNDF();

  // chi2topo_XiToXicPlus
  KFParticle kfpXi_XicPlus = kfpXiMinus_woMassConstForLam_wMassConstForXi;
  kfpXi_XicPlus.SetProductionVertex(kfpXicPlus_woMassConstForLam_wMassConstForXi);
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[1] = kfpXi_XicPlus.GetChi2()/kfpXi_XicPlus.GetNDF();

  // chi2topo_XicPlusToPV
  KFParticle kfpXicPlus_PV = kfpXicPlus_woMassConstForLam_wMassConstForXi;
  kfpXicPlus_PV.SetProductionVertex(PV);
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[2] = kfpXicPlus_PV.GetChi2()/kfpXicPlus_PV.GetNDF();

  // pt_XicPlus
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[3] = kfpXicPlus_woMassConstForLam_wMassConstForXi.GetPt();

  // mass_XicPlus
  Float_t massXicPlus, err_massXicPlus;
  kfpXicPlus_woMassConstForLam_wMassConstForXi.GetMass(massXicPlus, err_massXicPlus);
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[4] = massXicPlus;

  // SV_X_rec
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[5] = kfpXicPlus_woMassConstForLam_wMassConstForXi.GetX();
  // SV_Y_rec
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[6] = kfpXicPlus_woMassConstForLam_wMassConstForXi.GetY();
  // SV_Z_rec
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[7] = kfpXicPlus_woMassConstForLam_wMassConstForXi.GetZ();
  // SV_sigma_X_rec
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[8] = kfpXicPlus_woMassConstForLam_wMassConstForXi.GetErrX();
  // SV_sigma_Y_rec
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[9] = kfpXicPlus_woMassConstForLam_wMassConstForXi.GetErrY();
  // SV_sigma_Z_rec
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[10] = kfpXicPlus_woMassConstForLam_wMassConstForXi.GetErrZ();

  if (fIsMC) {
    // SV_X_MC
    Int_t labelPiFromXicPlus_HighPt = fabs(trackPiFromXicPlus_HighPt->GetLabel());
    AliAODMCParticle* mcPiFromXicPlus_HighPt = static_cast<AliAODMCParticle*>(mcArray->At(labelPiFromXicPlus_HighPt));
    fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[11] = mcPiFromXicPlus_HighPt->Xv();
    // SV_Y_MC
    fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[12] = mcPiFromXicPlus_HighPt->Yv();
    // SV_Z_MC
    fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[13] = mcPiFromXicPlus_HighPt->Zv();

    AliAODMCParticle* mcXicPlusCand = static_cast<AliAODMCParticle*>(mcArray->At(mcPiFromXicPlus_HighPt->GetMother()));
    if (abs(mcXicPlusCand->GetPdgCode()) == 4232) fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[14] = mcXicPlusCand->Pt();
    if (abs(mcXicPlusCand->GetPdgCode()) == 3324) {
      AliAODMCParticle* mcXicPlusCand_DecayToResonance = static_cast<AliAODMCParticle*>(mcArray->At(mcXicPlusCand->GetMother()));
      if (abs(mcXicPlusCand_DecayToResonance->GetPdgCode()) == 4232) fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[14] = mcXicPlusCand_DecayToResonance->Pt();
    }
  }

  // Source_XicPlus
  fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi[15] = lab_XicPlus;

  // event_ID
  if (fIsMC) fVar_XicPlus_EvtID = GetMCEventID();
  if (!fIsMC) fVar_XicPlus_EvtID = GetEventIdAsLong(AODEvent->GetHeader());

  fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi->Fill();

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillQATreeXicPlusFromCasc_wMassConstForLam_woMassConstForXi(KFParticle kfpLambda_wMassConst, KFParticle kfpXiMinus_wMassConstForLam_woMassConstForXi, KFParticle kfpXicPlus_wMassConstForLam_woMassConstForXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent)
{
  for (Int_t i=0;i<16;i++) {
    fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[i] = -9999.;
  }

  // chi2topo_LamToXi
  KFParticle kfpLambda_Xi = kfpLambda_wMassConst;
  kfpLambda_Xi.SetProductionVertex(kfpXiMinus_wMassConstForLam_woMassConstForXi);
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[0] = kfpLambda_Xi.GetChi2()/kfpLambda_Xi.GetNDF();

  // chi2topo_XiToXicPlus
  KFParticle kfpXi_XicPlus = kfpXiMinus_wMassConstForLam_woMassConstForXi;
  kfpXi_XicPlus.SetProductionVertex(kfpXicPlus_wMassConstForLam_woMassConstForXi);
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[1] = kfpXi_XicPlus.GetChi2()/kfpXi_XicPlus.GetNDF();

  // chi2topo_XicPlusToPV
  KFParticle kfpXicPlus_PV = kfpXicPlus_wMassConstForLam_woMassConstForXi;
  kfpXicPlus_PV.SetProductionVertex(PV);
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[2] = kfpXicPlus_PV.GetChi2()/kfpXicPlus_PV.GetNDF();

  // pt_XicPlus
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[3] = kfpXicPlus_wMassConstForLam_woMassConstForXi.GetPt();

  // mass_XicPlus
  Float_t massXicPlus, err_massXicPlus;
  kfpXicPlus_wMassConstForLam_woMassConstForXi.GetMass(massXicPlus, err_massXicPlus);
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[4] = massXicPlus;

  // SV_X_rec
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[5] = kfpXicPlus_wMassConstForLam_woMassConstForXi.GetX();
  // SV_Y_rec
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[6] = kfpXicPlus_wMassConstForLam_woMassConstForXi.GetY();
  // SV_Z_rec
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[7] = kfpXicPlus_wMassConstForLam_woMassConstForXi.GetZ();
  // SV_sigma_X_rec
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[8] = kfpXicPlus_wMassConstForLam_woMassConstForXi.GetErrX();
  // SV_sigma_Y_rec
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[9] = kfpXicPlus_wMassConstForLam_woMassConstForXi.GetErrY();
  // SV_sigma_Z_rec
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[10] = kfpXicPlus_wMassConstForLam_woMassConstForXi.GetErrZ();

  if (fIsMC) {
    // SV_X_MC
    Int_t labelPiFromXicPlus_HighPt = fabs(trackPiFromXicPlus_HighPt->GetLabel());
    AliAODMCParticle* mcPiFromXicPlus_HighPt = static_cast<AliAODMCParticle*>(mcArray->At(labelPiFromXicPlus_HighPt));
    fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[11] = mcPiFromXicPlus_HighPt->Xv();
    // SV_Y_MC
    fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[12] = mcPiFromXicPlus_HighPt->Yv();
    // SV_Z_MC
    fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[13] = mcPiFromXicPlus_HighPt->Zv();

    AliAODMCParticle* mcXicPlusCand = static_cast<AliAODMCParticle*>(mcArray->At(mcPiFromXicPlus_HighPt->GetMother()));
    if (abs(mcXicPlusCand->GetPdgCode()) == 4232) fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[14] = mcXicPlusCand->Pt();
    if (abs(mcXicPlusCand->GetPdgCode()) == 3324) {
      AliAODMCParticle* mcXicPlusCand_DecayToResonance = static_cast<AliAODMCParticle*>(mcArray->At(mcXicPlusCand->GetMother()));
      if (abs(mcXicPlusCand_DecayToResonance->GetPdgCode()) == 4232) fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[14] = mcXicPlusCand_DecayToResonance->Pt();
    }
  }

  // Source_XicPlus
  fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi[15] = lab_XicPlus;

  // event_ID
  if (fIsMC) fVar_XicPlus_EvtID = GetMCEventID();
  if (!fIsMC) fVar_XicPlus_EvtID = GetEventIdAsLong(AODEvent->GetHeader());

  fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi->Fill();

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillQATreeXicPlusFromCasc_wMassConstForLamAndXi(KFParticle kfpLambda_wMassConst, KFParticle kfpXiMinus_wMassConstForLamAndXi, KFParticle kfpXicPlus_wMassConstForLamAndXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent)
{
  for (Int_t i=0;i<16;i++) {
    fVar_XicPlus_QA_wMassConstForLamAndXi[i] = -9999.;
  }

  // chi2topo_LamToXi
  KFParticle kfpLambda_Xi = kfpLambda_wMassConst;
  kfpLambda_Xi.SetProductionVertex(kfpXiMinus_wMassConstForLamAndXi);
  fVar_XicPlus_QA_wMassConstForLamAndXi[0] = kfpLambda_Xi.GetChi2()/kfpLambda_Xi.GetNDF();

  // chi2topo_XiToXicPlus
  KFParticle kfpXi_XicPlus = kfpXiMinus_wMassConstForLamAndXi;
  kfpXi_XicPlus.SetProductionVertex(kfpXicPlus_wMassConstForLamAndXi);
  fVar_XicPlus_QA_wMassConstForLamAndXi[1] = kfpXi_XicPlus.GetChi2()/kfpXi_XicPlus.GetNDF();

  // chi2topo_XicPlusToPV
  KFParticle kfpXicPlus_PV = kfpXicPlus_wMassConstForLamAndXi;
  kfpXicPlus_PV.SetProductionVertex(PV);
  fVar_XicPlus_QA_wMassConstForLamAndXi[2] = kfpXicPlus_PV.GetChi2()/kfpXicPlus_PV.GetNDF();

  // pt_XicPlus
  fVar_XicPlus_QA_wMassConstForLamAndXi[3] = kfpXicPlus_wMassConstForLamAndXi.GetPt();

  // mass_XicPlus
  Float_t massXicPlus, err_massXicPlus;
  kfpXicPlus_wMassConstForLamAndXi.GetMass(massXicPlus, err_massXicPlus);
  fVar_XicPlus_QA_wMassConstForLamAndXi[4] = massXicPlus;
  
  // SV_X_rec
  fVar_XicPlus_QA_wMassConstForLamAndXi[5] = kfpXicPlus_wMassConstForLamAndXi.GetX();
  // SV_Y_rec
  fVar_XicPlus_QA_wMassConstForLamAndXi[6] = kfpXicPlus_wMassConstForLamAndXi.GetY();
  // SV_Z_rec
  fVar_XicPlus_QA_wMassConstForLamAndXi[7] = kfpXicPlus_wMassConstForLamAndXi.GetZ();
  // SV_sigma_X_rec
  fVar_XicPlus_QA_wMassConstForLamAndXi[8] = kfpXicPlus_wMassConstForLamAndXi.GetErrX();
  // SV_sigma_Y_rec
  fVar_XicPlus_QA_wMassConstForLamAndXi[9] = kfpXicPlus_wMassConstForLamAndXi.GetErrY();
  // SV_sigma_Z_rec
  fVar_XicPlus_QA_wMassConstForLamAndXi[10] = kfpXicPlus_wMassConstForLamAndXi.GetErrZ();

  if (fIsMC) {
    // SV_X_MC
    Int_t labelPiFromXicPlus_HighPt = fabs(trackPiFromXicPlus_HighPt->GetLabel());
    AliAODMCParticle* mcPiFromXicPlus_HighPt = static_cast<AliAODMCParticle*>(mcArray->At(labelPiFromXicPlus_HighPt));
    fVar_XicPlus_QA_wMassConstForLamAndXi[11] = mcPiFromXicPlus_HighPt->Xv();
    // SV_Y_MC
    fVar_XicPlus_QA_wMassConstForLamAndXi[12] = mcPiFromXicPlus_HighPt->Yv();
    // SV_Z_MC
    fVar_XicPlus_QA_wMassConstForLamAndXi[13] = mcPiFromXicPlus_HighPt->Zv();

    AliAODMCParticle* mcXicPlusCand = static_cast<AliAODMCParticle*>(mcArray->At(mcPiFromXicPlus_HighPt->GetMother()));
    if (abs(mcXicPlusCand->GetPdgCode()) == 4232) fVar_XicPlus_QA_wMassConstForLamAndXi[14] = mcXicPlusCand->Pt();
    if (abs(mcXicPlusCand->GetPdgCode()) == 3324) {
      AliAODMCParticle* mcXicPlusCand_DecayToResonance = static_cast<AliAODMCParticle*>(mcArray->At(mcXicPlusCand->GetMother()));
      if (abs(mcXicPlusCand_DecayToResonance->GetPdgCode()) == 4232) fVar_XicPlus_QA_wMassConstForLamAndXi[14] = mcXicPlusCand_DecayToResonance->Pt();
    }
  }

  // Source_XicPlus
  fVar_XicPlus_QA_wMassConstForLamAndXi[15] = lab_XicPlus;

  // event_ID
  if (fIsMC) fVar_XicPlus_EvtID = GetMCEventID();
  if (!fIsMC) fVar_XicPlus_EvtID = GetEventIdAsLong(AODEvent->GetHeader());

  fTree_XicPlus_QA_wMassConstForLamAndXi->Fill();

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassConstForXi(KFParticle kfpLambda_wMassConst_To_Xi_wMassConst, KFParticle kfpXiMinus_wMassAndTopoConstForLam_wMassConstForXi, KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent)
{
  for (Int_t i=0;i<16;i++) {
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[i] = -9999.;
  }

  // chi2topo_LamToXi
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[0] = kfpLambda_wMassConst_To_Xi_wMassConst.GetChi2()/kfpLambda_wMassConst_To_Xi_wMassConst.GetNDF();

  // chi2mass_Xi
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[1] = kfpXiMinus_wMassAndTopoConstForLam_wMassConstForXi.GetChi2()/kfpXiMinus_wMassAndTopoConstForLam_wMassConstForXi.GetNDF();

  // chi2topo_XicPlusToPV
  KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi_PV = kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi;
  kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi_PV.SetProductionVertex(PV);
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[2] = kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi_PV.GetChi2()/kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi_PV.GetNDF();

  // pt_XicPlus
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[3] = kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetPt();

  // mass_XicPlus
  Float_t massXicPlus, err_massXicPlus;
  kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetMass(massXicPlus, err_massXicPlus);
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[4] = massXicPlus;
  
  // SV_X_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[5] = kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetX();
  // SV_Y_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[6] = kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetY();
  // SV_Z_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[7] = kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetZ();
  // SV_sigma_X_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[8] = kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetErrX();
  // SV_sigma_Y_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[9] = kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetErrY();
  // SV_sigma_Z_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[10] = kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi.GetErrZ();

  if (fIsMC) {
    // SV_X_MC
    Int_t labelPiFromXicPlus_HighPt = fabs(trackPiFromXicPlus_HighPt->GetLabel());
    AliAODMCParticle* mcPiFromXicPlus_HighPt = static_cast<AliAODMCParticle*>(mcArray->At(labelPiFromXicPlus_HighPt));
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[11] = mcPiFromXicPlus_HighPt->Xv();
    // SV_Y_MC
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[12] = mcPiFromXicPlus_HighPt->Yv();
    // SV_Z_MC
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[13] = mcPiFromXicPlus_HighPt->Zv();

    AliAODMCParticle* mcXicPlusCand = static_cast<AliAODMCParticle*>(mcArray->At(mcPiFromXicPlus_HighPt->GetMother()));
    if (abs(mcXicPlusCand->GetPdgCode()) == 4232) fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[14] = mcXicPlusCand->Pt();
    if (abs(mcXicPlusCand->GetPdgCode()) == 3324) {
      AliAODMCParticle* mcXicPlusCand_DecayToResonance = static_cast<AliAODMCParticle*>(mcArray->At(mcXicPlusCand->GetMother()));
      if (abs(mcXicPlusCand_DecayToResonance->GetPdgCode()) == 4232) fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[14] = mcXicPlusCand_DecayToResonance->Pt();
    }
  }

  // Source_XicPlus
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi[15] = lab_XicPlus;

  // event_ID
  if (fIsMC) fVar_XicPlus_EvtID = GetMCEventID();
  if (!fIsMC) fVar_XicPlus_EvtID = GetEventIdAsLong(AODEvent->GetHeader());

  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi->Fill();

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassAndTopoConstForXi(KFParticle kfpLambda_wMassConst_To_Xi_wMassConst, KFParticle kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi, KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent)
{
  for (Int_t i=0;i<16;i++) {
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[i] = -9999.;
  }

  // chi2topo_LamToXi
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[0] = kfpLambda_wMassConst_To_Xi_wMassConst.GetChi2()/kfpLambda_wMassConst_To_Xi_wMassConst.GetNDF();

  // chi2topo_XiToPV
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[1] = kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetChi2()/kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetNDF();

  // chi2topo_XicPlusToPV
  KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_PV = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi;
  kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_PV.SetProductionVertex(PV);
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[2] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_PV.GetChi2()/kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_PV.GetNDF();

  // pt_XicPlus
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[3] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_PV.GetPt();

  // mass_XicPlus
  Float_t massXicPlus, err_massXicPlus;
  kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_PV.GetMass(massXicPlus, err_massXicPlus);
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[4] = massXicPlus;
  
  // SV_X_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[5] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetX();
  // SV_Y_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[6] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetY();
  // SV_Z_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[7] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetZ();
  // SV_sigma_X_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[8] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetErrX();
  // SV_sigma_Y_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[9] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetErrY();
  // SV_sigma_Z_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[10] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetErrZ();

  if (fIsMC) {
    // SV_X_MC
    Int_t labelPiFromXicPlus_HighPt = fabs(trackPiFromXicPlus_HighPt->GetLabel());
    AliAODMCParticle* mcPiFromXicPlus_HighPt = static_cast<AliAODMCParticle*>(mcArray->At(labelPiFromXicPlus_HighPt));
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[11] = mcPiFromXicPlus_HighPt->Xv();
    // SV_Y_MC
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[12] = mcPiFromXicPlus_HighPt->Yv();
    // SV_Z_MC
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[13] = mcPiFromXicPlus_HighPt->Zv();

    AliAODMCParticle* mcXicPlusCand = static_cast<AliAODMCParticle*>(mcArray->At(mcPiFromXicPlus_HighPt->GetMother()));
    if (abs(mcXicPlusCand->GetPdgCode()) == 4232) fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[14] = mcXicPlusCand->Pt();
    if (abs(mcXicPlusCand->GetPdgCode()) == 3324) {
      AliAODMCParticle* mcXicPlusCand_DecayToResonance = static_cast<AliAODMCParticle*>(mcArray->At(mcXicPlusCand->GetMother()));
      if (abs(mcXicPlusCand_DecayToResonance->GetPdgCode()) == 4232) fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[14] = mcXicPlusCand_DecayToResonance->Pt();
    }
  }

  // Source_XicPlus
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi[15] = lab_XicPlus;

  // event_ID
  if (fIsMC) fVar_XicPlus_EvtID = GetMCEventID();
  if (!fIsMC) fVar_XicPlus_EvtID = GetEventIdAsLong(AODEvent->GetHeader());

  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi->Fill();

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicPlusToXi2PifromKFP::FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic(KFParticle kfpLambda_wMassConst_To_Xi_wMassConst, KFParticle kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi, KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, AliAODTrack *trackPiFromXicPlus_LowPt, AliAODTrack *trackPiFromXiOrKaonFromOmega, AliAODTrack *trackPrFromLam, AliAODTrack *trackPiFromLam, KFParticle kfpBP_HighPt, KFParticle kfpBP_LowPt, KFParticle kfpPion_ForXi, KFParticle kfpProton_ForLam, KFParticle kfpPion_ForLam, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent)
{
  for (Int_t i=0;i<31;i++) {
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[i] = -9999.;
  }

  // chi2topo_LamToXi
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[0] = kfpLambda_wMassConst_To_Xi_wMassConst.GetChi2()/kfpLambda_wMassConst_To_Xi_wMassConst.GetNDF();

  // chi2topo_XiToPV
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[1] = kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetChi2()/kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi.GetNDF();

  // chi2topo_XicPlusToPV
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[2] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetChi2()/kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetNDF();

  // pt_XicPlus
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[3] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetPt();

  // mass_XicPlus
  Float_t massXicPlus, err_massXicPlus;
  kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic.GetMass(massXicPlus, err_massXicPlus);
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[4] = massXicPlus;
  
  KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic_TransportToDecayVertex = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic;
  kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic_TransportToDecayVertex.TransportToDecayVertex();

  // SV_X_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[5] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic_TransportToDecayVertex.GetX();
  // SV_Y_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[6] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic_TransportToDecayVertex.GetY();
  // SV_Z_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[7] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic_TransportToDecayVertex.GetZ();
  // SV_sigma_X_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[8] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic_TransportToDecayVertex.GetErrX();
  // SV_sigma_Y_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[9] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic_TransportToDecayVertex.GetErrY();
  // SV_sigma_Z_rec
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[10] = kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic_TransportToDecayVertex.GetErrZ();

  if (fIsMC) {
    // SV_X_MC
    Int_t labelPiFromXicPlus_HighPt = fabs(trackPiFromXicPlus_HighPt->GetLabel());
    AliAODMCParticle* mcPiFromXicPlus_HighPt = static_cast<AliAODMCParticle*>(mcArray->At(labelPiFromXicPlus_HighPt));
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[11] = mcPiFromXicPlus_HighPt->Xv();
    // SV_Y_MC
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[12] = mcPiFromXicPlus_HighPt->Yv();
    // SV_Z_MC
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[13] = mcPiFromXicPlus_HighPt->Zv();
    // PV_X_MC
    Double_t PV_Gen[3]={0.};
    AliAODMCHeader *mcHeader = (AliAODMCHeader*)AODEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    mcHeader->GetVertex(PV_Gen);
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[14] = PV_Gen[0];
    // PV_Y_MC
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[15] = PV_Gen[1];
    // PV_Z_MC
    fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[16] = PV_Gen[2];

    AliAODMCParticle* mcXicPlusCand = static_cast<AliAODMCParticle*>(mcArray->At(mcPiFromXicPlus_HighPt->GetMother()));
    if (abs(mcXicPlusCand->GetPdgCode()) == 4232) fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[29] = mcXicPlusCand->Pt();
    if (abs(mcXicPlusCand->GetPdgCode()) == 3324) {
      AliAODMCParticle* mcXicPlusCand_DecayToResonance = static_cast<AliAODMCParticle*>(mcArray->At(mcXicPlusCand->GetMother()));
      if (abs(mcXicPlusCand_DecayToResonance->GetPdgCode()) == 4232) fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[29] = mcXicPlusCand_DecayToResonance->Pt();
    }
  }

  // PV_rec
  Double_t pos_PV[3]={0.}, cov_PV[6]={0.};
  fpVtx->GetXYZ(pos_PV);
  fpVtx->GetCovarianceMatrix(cov_PV);
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[17] = pos_PV[0]; // PV: X_{rec}
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[18] = pos_PV[1]; // PV: Y_{rec}
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[19] = pos_PV[2]; // PV: Z_{rec}
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[20] = sqrt(cov_PV[0]); // PV: sigma_X^{rec}
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[21] = sqrt(cov_PV[2]); // PV: sigma_Y^{rec}
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[22] = sqrt(cov_PV[5]); // PV: sigma_Z^{rec}

  // === Recalculate PV after removing Xic+ daughters and adding Xic+ (w/o KF Refit) ===
  KFParticle recalPV_woKFrefit = PV;
  if (trackPiFromXicPlus_HighPt->GetUsedForPrimVtxFit()) kfpBP_HighPt.SubtractFromVertex(recalPV_woKFrefit);
  if (trackPiFromXicPlus_LowPt->GetUsedForPrimVtxFit()) kfpBP_LowPt.SubtractFromVertex(recalPV_woKFrefit);
  if (trackPiFromXiOrKaonFromOmega->GetUsedForPrimVtxFit()) kfpPion_ForXi.SubtractFromVertex(recalPV_woKFrefit);
  if (trackPrFromLam->GetUsedForPrimVtxFit()) kfpProton_ForLam.SubtractFromVertex(recalPV_woKFrefit);
  if (trackPiFromLam->GetUsedForPrimVtxFit()) kfpPion_ForLam.SubtractFromVertex(recalPV_woKFrefit);
  recalPV_woKFrefit.AddDaughter(kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic);

  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[23] = recalPV_woKFrefit.GetX(); // recalPV_woKFrefit: X
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[24] = recalPV_woKFrefit.GetY(); // recalPV_woKFrefit: Y
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[25] = recalPV_woKFrefit.GetZ(); // recalPV_woKFrefit: Z
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[26] = recalPV_woKFrefit.GetErrX(); // recalPV_woKFrefit: sigma_X
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[27] = recalPV_woKFrefit.GetErrY(); // recalPV_woKFrefit: sigma_Y
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[28] = recalPV_woKFrefit.GetErrZ(); // recalPV_woKFrefit: sigma_Z

  recalPV_woKFrefit.Clear();

  // Source_XicPlus
  fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic[30] = lab_XicPlus;

  // event_ID
  if (fIsMC) fVar_XicPlus_EvtID = GetMCEventID();
  if (!fIsMC) fVar_XicPlus_EvtID = GetEventIdAsLong(AODEvent->GetHeader());

  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic->Fill();

  return;
}

//_____________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEXicPlusToXi2PifromKFP::PrimaryVertex(const TObjArray *trkArray, AliVEvent *event)
{
  //
  //Used only for pp
  //copied from AliAnalysisVertexingHF (except for the following 3 lines)
  //

  Bool_t fRecoPrimVtxSkippingTrks = kTRUE;
  Bool_t fRmTrksFromPrimVtx = kFALSE;

  AliESDVertex *vertexESD = NULL;
  AliAODVertex *vertexAOD = NULL;
  
  Double_t pos_PV[3]={0.}, cov_PV[6]={0.};
  fpVtx->GetXYZ(pos_PV);
  fpVtx->GetCovarianceMatrix(cov_PV);

  AliESDVertex *fV1 = new AliESDVertex(pos_PV,cov_PV,100.,100,fpVtx->GetName());
  
  Int_t nTrks = trkArray->GetEntriesFast();
  AliVertexerTracks *vertexer = new AliVertexerTracks(event->GetMagneticField());
    
  if(fRecoPrimVtxSkippingTrks) {
    // recalculating the vertex
      
    if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraint")) {
	    Float_t diamondcovxy[3];
	    event->GetDiamondCovXY(diamondcovxy);
      Double_t pos[3]={event->GetDiamondX(),event->GetDiamondY(),0.};
	    Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
      AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
      vertexer->SetVtxStart(diamond);
      delete diamond; diamond=NULL;
      if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraintOnlyFitter")) 
      vertexer->SetOnlyFitter();
    }
    Int_t skipped[1000];
    Int_t nTrksToSkip=0,id;
    AliExternalTrackParam *t = 0;
    for(Int_t i=0; i<nTrks; i++) {
      t = (AliExternalTrackParam*)trkArray->UncheckedAt(i);
      id = (Int_t)t->GetID();
      if(id<0) continue;
      skipped[nTrksToSkip++] = id;
    }
    // TEMPORARY FIX
    // For AOD, skip also tracks without covariance matrix
    Double_t covtest[21];
    for(Int_t j=0; j<event->GetNumberOfTracks(); j++) {
	    AliVTrack *vtrack = (AliVTrack*)event->GetTrack(j);
	    if(!vtrack->GetCovarianceXYZPxPyPz(covtest)) {
	      id = (Int_t)vtrack->GetID();
	      if(id<0) continue;
	      skipped[nTrksToSkip++] = id;
    	}
    }
    for(Int_t ijk=nTrksToSkip; ijk<1000; ijk++) skipped[ijk]=-1;
    //
    vertexer->SetSkipTracks(nTrksToSkip,skipped);
    vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(event); 

  } else if(fRmTrksFromPrimVtx && nTrks>0) { 
    // removing the prongs tracks

    TObjArray rmArray(nTrks);
    UShort_t *rmId = new UShort_t[nTrks];
    AliESDtrack *esdTrack = 0;
    AliESDtrack *t = 0;
    for(Int_t i=0; i<nTrks; i++) {
      t = (AliESDtrack*)trkArray->UncheckedAt(i);
      esdTrack = new AliESDtrack(*t);
      rmArray.AddLast(esdTrack);
      if(esdTrack->GetID()>=0) {
        rmId[i]=(UShort_t)esdTrack->GetID();
      } else {
        rmId[i]=9999;
      }
    }
    Float_t diamondxy[2]={static_cast<Float_t>(event->GetDiamondX()),static_cast<Float_t>(event->GetDiamondY())};
    vertexESD = vertexer->RemoveTracksFromVertex(fV1,&rmArray,rmId,diamondxy);
    delete [] rmId; rmId=NULL;
    rmArray.Delete();

  }
    
  delete vertexer; vertexer=NULL;
  if(!vertexESD) return vertexAOD;
  if(vertexESD->GetNContributors()<=0) { 
    //AliDebug(2,"vertexing failed"); 
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }
  
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();

  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF);
  vertexAOD->SetNContributors(vertexESD->GetNContributors());
  
  delete vertexESD; vertexESD=NULL;
  delete fV1; fV1=NULL;
  return vertexAOD;
}

//_____________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEXicPlusToXi2PifromKFP::CallPrimaryVertex(AliAODcascade *casc, AliAODTrack *trk1, AliAODTrack *trk2, AliAODEvent *aodEvent)
{
  //
  // Make an array of tracks which should not be used in primary vertex calculation and 
  // Call PrimaryVertex function
  //
  
  TObjArray *twoTrackArrayPlusXi = new TObjArray(5);

  AliESDtrack *cptrk1 = new AliESDtrack((AliVTrack*)trk1);
  twoTrackArrayPlusXi->AddAt(cptrk1,0);
  AliESDtrack *cptrk2 = new AliESDtrack((AliVTrack*)trk2);
  twoTrackArrayPlusXi->AddAt(cptrk2,1);

  AliESDtrack *cascptrack = new AliESDtrack((AliVTrack*)casc->GetDaughter(0));
  twoTrackArrayPlusXi->AddAt(cascptrack,2);
  AliESDtrack *cascntrack = new AliESDtrack((AliVTrack*)casc->GetDaughter(1));
  twoTrackArrayPlusXi->AddAt(cascntrack,3);
  AliESDtrack *cascbtrack = new AliESDtrack((AliVTrack*)casc->GetDecayVertexXi()->GetDaughter(0));
  twoTrackArrayPlusXi->AddAt(cascbtrack,4);

  AliAODVertex *newVert  = PrimaryVertex(twoTrackArrayPlusXi,aodEvent);

  for(Int_t i=0;i<5;i++)
  {
    AliESDtrack *tesd = (AliESDtrack*)twoTrackArrayPlusXi->UncheckedAt(i);
    delete tesd;
  }
  twoTrackArrayPlusXi->Clear();
  delete twoTrackArrayPlusXi;
  
  return newVert;
}

//_____________________________________________________________________________
ULong64_t AliAnalysisTaskSEXicPlusToXi2PifromKFP::GetEventIdAsLong(AliVHeader* header)
{
	// To have a unique id for each event in a run!
	// Modified from AliRawReader.h
	return ((ULong64_t)header->GetBunchCrossNumber()+
			(ULong64_t)header->GetOrbitNumber()*3564+
			(ULong64_t)header->GetPeriodNumber()*16777215*3564);
}

//_____________________________________________________________________________
unsigned int AliAnalysisTaskSEXicPlusToXi2PifromKFP::GetMCEventID()
{
  TString currentfilename = ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetCurrentFile()))->GetName();
  if(!fFileName.EqualTo(currentfilename)) {
    fEventNumber = 0;
    fFileName = currentfilename;
    TObjArray *path = fFileName.Tokenize("/");
    TString s = ((TObjString*)path->At( ((path->GetLast())-1) ))->GetString();
    fDirNumber = (unsigned int)s.Atoi();
    delete path;
  }
  Long64_t ev_number = Entry();
  if(fIsMC){
    ev_number = fEventNumber;
  }
  unsigned int evID = (unsigned int)ev_number + (unsigned int)(fDirNumber<<17);
  fEventNumber++;
  return evID;
}
