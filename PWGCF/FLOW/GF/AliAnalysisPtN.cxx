/* AliAnaysisTaskPtN
 * Maintainer: Yifan Zhang
 * calculating the mean Pt and fluctuations with respect to Nch
 */

#include "TChain.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliAODEvent.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisPtN.h"
#include <cmath>
#include <TComplex.h>
# include <TRandom3.h>
#include <iostream>
//#include "CorrelationCalculator.h"

class AliAnalysisPtN;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisPtN) // classimp: necessary for root

AliAnalysisPtN::AliAnalysisPtN() : AliAnalysisTaskSE(), 
    fAOD(nullptr), 
    fOutputList(nullptr), 
    fBstList(nullptr),
    fWeightsListNUE(nullptr),
    fWeightNUE(nullptr),
    fMinPt(0.2),
    fMaxPt(3.0),
    fVz(10.0),
    fDCAz(2.0),
    fDCAxy(7.0),
    fTPCcls(70),
    filterBit(96),
    fTrigger(1),
    ctrType("V0M"),
    fPeriod("LHC15o"),
    fNUE("LHC20e3a"),
    fSysflg(0),
    //correlator(),
    fTestNonWeight(nullptr),
    fNchDistri(nullptr),
    fPtQA(nullptr),
    fDCAxyQA(nullptr),
    fDCAzQA(nullptr),
    fetaQA(nullptr),
    fTPCclsQA(nullptr),
    fDCAxyPt(nullptr),
    fDCAzPt(nullptr),
    fPtNch(nullptr),
    dPtNch(nullptr),
    dPt2Nch(nullptr),
    dPt3Nch(nullptr),
    dPt4Nch(nullptr),
    dPtNchNB(nullptr),
    dPtNch_0(nullptr),
    dPtNch_1(nullptr),
    dPtNch_2(nullptr),
    dPtNch_3(nullptr),
    dPtNch_4(nullptr),
    dPtNch_5(nullptr),
    dPtNch_6(nullptr),
    dPtNch_7(nullptr),
    dPtNch_8(nullptr),
    dPtNch_9(nullptr),
    dPt2Nch_0(nullptr),
    dPt2Nch_1(nullptr),
    dPt2Nch_2(nullptr),
    dPt2Nch_3(nullptr),
    dPt2Nch_4(nullptr),
    dPt2Nch_5(nullptr),
    dPt2Nch_6(nullptr),
    dPt2Nch_7(nullptr),
    dPt2Nch_8(nullptr),
    dPt2Nch_9(nullptr),
    dPt3Nch_0(nullptr),
    dPt3Nch_1(nullptr),
    dPt3Nch_2(nullptr),
    dPt3Nch_3(nullptr),
    dPt3Nch_4(nullptr),
    dPt3Nch_5(nullptr),
    dPt3Nch_6(nullptr),
    dPt3Nch_7(nullptr),
    dPt3Nch_8(nullptr),
    dPt3Nch_9(nullptr),
    dPt4Nch_0(nullptr),
    dPt4Nch_1(nullptr),
    dPt4Nch_2(nullptr),
    dPt4Nch_3(nullptr),
    dPt4Nch_4(nullptr),
    dPt4Nch_5(nullptr),
    dPt4Nch_6(nullptr),
    dPt4Nch_7(nullptr),
    dPt4Nch_8(nullptr),
    dPt4Nch_9(nullptr),
    fPtNchUCC(nullptr),
    TestPtCtr(nullptr),
    TestPt2Ctr(nullptr),
    TestPt3Ctr(nullptr),
    TestPt4Ctr(nullptr),
    TestNchCtr(nullptr),
    TestNchSelectedCtr(nullptr),
    TestPtCtr_0(nullptr),
    TestPtCtr_1(nullptr),
    TestPtCtr_2(nullptr),
    TestPtCtr_3(nullptr),
    TestPtCtr_4(nullptr),
    TestPtCtr_5(nullptr),
    TestPtCtr_6(nullptr),
    TestPtCtr_7(nullptr),
    TestPtCtr_8(nullptr),
    TestPtCtr_9(nullptr),
    TestPt2Ctr_0(nullptr),
    TestPt2Ctr_1(nullptr),
    TestPt2Ctr_2(nullptr),
    TestPt2Ctr_3(nullptr),
    TestPt2Ctr_4(nullptr),
    TestPt2Ctr_5(nullptr),
    TestPt2Ctr_6(nullptr),
    TestPt2Ctr_7(nullptr),
    TestPt2Ctr_8(nullptr),
    TestPt2Ctr_9(nullptr),
    TestPt3Ctr_0(nullptr),
    TestPt3Ctr_1(nullptr),
    TestPt3Ctr_2(nullptr),
    TestPt3Ctr_3(nullptr),
    TestPt3Ctr_4(nullptr),
    TestPt3Ctr_5(nullptr),
    TestPt3Ctr_6(nullptr),
    TestPt3Ctr_7(nullptr),
    TestPt3Ctr_8(nullptr),
    TestPt3Ctr_9(nullptr),
    TestPt4Ctr_0(nullptr),
    TestPt4Ctr_1(nullptr),
    TestPt4Ctr_2(nullptr),
    TestPt4Ctr_3(nullptr),
    TestPt4Ctr_4(nullptr),
    TestPt4Ctr_5(nullptr),
    TestPt4Ctr_6(nullptr),
    TestPt4Ctr_7(nullptr),
    TestPt4Ctr_8(nullptr),
    TestPt4Ctr_9(nullptr),
    TestNchCtr_0(nullptr),
    TestNchCtr_1(nullptr),
    TestNchCtr_2(nullptr),
    TestNchCtr_3(nullptr),
    TestNchCtr_4(nullptr),
    TestNchCtr_5(nullptr),
    TestNchCtr_6(nullptr),
    TestNchCtr_7(nullptr),
    TestNchCtr_8(nullptr),
    TestNchCtr_9(nullptr),
    TestNchSelectedCtr_0(nullptr),
    TestNchSelectedCtr_1(nullptr),
    TestNchSelectedCtr_2(nullptr),
    TestNchSelectedCtr_3(nullptr),
    TestNchSelectedCtr_4(nullptr),
    TestNchSelectedCtr_5(nullptr),
    TestNchSelectedCtr_6(nullptr),
    TestNchSelectedCtr_7(nullptr),
    TestNchSelectedCtr_8(nullptr),
    TestNchSelectedCtr_9(nullptr),
    fEventCuts(0),
    multSelection(nullptr),
    radm(32213)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisPtN::AliAnalysisPtN(const char* name) : AliAnalysisTaskSE(name),
    fAOD(nullptr), 
    fOutputList(nullptr), 
    fBstList(nullptr),
    fWeightsListNUE(nullptr),
    fWeightNUE(nullptr),
    fMinPt(0.2),
    fMaxPt(3.0),
    fVz(10.0),
    fDCAz(2.0),
    fDCAxy(7.0),
    fTPCcls(70),
    filterBit(96),
    fTrigger(1),
    ctrType("V0M"),
    fPeriod("LHC15o"),
    fNUE("LHC20e3a"),
    fSysflg(0),
    //correlator(),
    fTestNonWeight(nullptr),
    fNchDistri(nullptr),
    fPtQA(nullptr),
    fDCAxyQA(nullptr),
    fDCAzQA(nullptr),
    fetaQA(nullptr),
    fTPCclsQA(nullptr),
    fDCAxyPt(nullptr),
    fDCAzPt(nullptr),
    fPtNch(nullptr),
    dPtNch(nullptr),
    dPt2Nch(nullptr),
    dPt3Nch(nullptr),
    dPt4Nch(nullptr),
    dPtNchNB(nullptr),
    dPtNch_0(nullptr),
    dPtNch_1(nullptr),
    dPtNch_2(nullptr),
    dPtNch_3(nullptr),
    dPtNch_4(nullptr),
    dPtNch_5(nullptr),
    dPtNch_6(nullptr),
    dPtNch_7(nullptr),
    dPtNch_8(nullptr),
    dPtNch_9(nullptr),
    dPt2Nch_0(nullptr),
    dPt2Nch_1(nullptr),
    dPt2Nch_2(nullptr),
    dPt2Nch_3(nullptr),
    dPt2Nch_4(nullptr),
    dPt2Nch_5(nullptr),
    dPt2Nch_6(nullptr),
    dPt2Nch_7(nullptr),
    dPt2Nch_8(nullptr),
    dPt2Nch_9(nullptr),
    dPt3Nch_0(nullptr),
    dPt3Nch_1(nullptr),
    dPt3Nch_2(nullptr),
    dPt3Nch_3(nullptr),
    dPt3Nch_4(nullptr),
    dPt3Nch_5(nullptr),
    dPt3Nch_6(nullptr),
    dPt3Nch_7(nullptr),
    dPt3Nch_8(nullptr),
    dPt3Nch_9(nullptr),
    dPt4Nch_0(nullptr),
    dPt4Nch_1(nullptr),
    dPt4Nch_2(nullptr),
    dPt4Nch_3(nullptr),
    dPt4Nch_4(nullptr),
    dPt4Nch_5(nullptr),
    dPt4Nch_6(nullptr),
    dPt4Nch_7(nullptr),
    dPt4Nch_8(nullptr),
    dPt4Nch_9(nullptr),
    fPtNchUCC(nullptr),
    TestPtCtr(nullptr),
    TestPt2Ctr(nullptr),
    TestPt3Ctr(nullptr),
    TestPt4Ctr(nullptr),
    TestNchCtr(nullptr),
    TestNchSelectedCtr(nullptr),
    TestPtCtr_0(nullptr),
    TestPtCtr_1(nullptr),
    TestPtCtr_2(nullptr),
    TestPtCtr_3(nullptr),
    TestPtCtr_4(nullptr),
    TestPtCtr_5(nullptr),
    TestPtCtr_6(nullptr),
    TestPtCtr_7(nullptr),
    TestPtCtr_8(nullptr),
    TestPtCtr_9(nullptr),
    TestPt2Ctr_0(nullptr),
    TestPt2Ctr_1(nullptr),
    TestPt2Ctr_2(nullptr),
    TestPt2Ctr_3(nullptr),
    TestPt2Ctr_4(nullptr),
    TestPt2Ctr_5(nullptr),
    TestPt2Ctr_6(nullptr),
    TestPt2Ctr_7(nullptr),
    TestPt2Ctr_8(nullptr),
    TestPt2Ctr_9(nullptr),
    TestPt3Ctr_0(nullptr),
    TestPt3Ctr_1(nullptr),
    TestPt3Ctr_2(nullptr),
    TestPt3Ctr_3(nullptr),
    TestPt3Ctr_4(nullptr),
    TestPt3Ctr_5(nullptr),
    TestPt3Ctr_6(nullptr),
    TestPt3Ctr_7(nullptr),
    TestPt3Ctr_8(nullptr),
    TestPt3Ctr_9(nullptr),
    TestPt4Ctr_0(nullptr),
    TestPt4Ctr_1(nullptr),
    TestPt4Ctr_2(nullptr),
    TestPt4Ctr_3(nullptr),
    TestPt4Ctr_4(nullptr),
    TestPt4Ctr_5(nullptr),
    TestPt4Ctr_6(nullptr),
    TestPt4Ctr_7(nullptr),
    TestPt4Ctr_8(nullptr),
    TestPt4Ctr_9(nullptr),
    TestNchCtr_0(nullptr),
    TestNchCtr_1(nullptr),
    TestNchCtr_2(nullptr),
    TestNchCtr_3(nullptr),
    TestNchCtr_4(nullptr),
    TestNchCtr_5(nullptr),
    TestNchCtr_6(nullptr),
    TestNchCtr_7(nullptr),
    TestNchCtr_8(nullptr),
    TestNchCtr_9(nullptr),
    TestNchSelectedCtr_0(nullptr),
    TestNchSelectedCtr_1(nullptr),
    TestNchSelectedCtr_2(nullptr),
    TestNchSelectedCtr_3(nullptr),
    TestNchSelectedCtr_4(nullptr),
    TestNchSelectedCtr_5(nullptr),
    TestNchSelectedCtr_6(nullptr),
    TestNchSelectedCtr_7(nullptr),
    TestNchSelectedCtr_8(nullptr),
    TestNchSelectedCtr_9(nullptr),
    fEventCuts(0),
    multSelection(nullptr),
    radm(32213)
{
    // constructor
    
    DefineInput(0, TChain::Class()); 
    DefineInput(1, TList::Class()); 
      // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisPtN::~AliAnalysisPtN()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisPtN::UserCreateOutputObjects()
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
    fEventCuts.AddQAplotsToList(fOutputList);
    
    fTestNonWeight = new TH1F("fHisPhi", "fHisPhi", 100, 0, 6.28);
    fNchDistri = new TH1F("fNchDistri", "fNchDistri", 90, 0, 4500);
    fPtQA = new TH1F("fPtQA", "fPtQA", 100, 0, 5);
    fDCAxyQA = new TH1F("fDCAxyQA", "fDCAxyQA", 100, 0, 0.5);
    fDCAzQA = new TH1F("fDCAzQA", "fDCAzQA", 200, -2.5, 2.5);
    fetaQA = new TH1F("fetaQA", "fetaQA", 100, -1, 1);
    fTPCclsQA = new TH1F("fTPCclsQA", "fTPCclsQA", 200, 0, 200);
    fDCAxyPt = new TH2F("fDCAxyPt", "fDCAxyPt", 500, 0, 0.5, 500, 0.2, 3);
    fDCAzPt = new TH2F("fDCAzPt", "fDCAzPt", 1000, -2, 2, 500, 0.2, 3);
    fPtNch = new TH2F("fPtNch", "fPtNch", 900, 0, 4500, 500, 0, 3);
    dPtNch = new TProfile("dPtNch", "dPtNch", 90, 0, 4500);
    dPt2Nch = new TProfile("dPt2Nch", "dPt2Nch", 90, 0, 4500);
    dPt3Nch = new TProfile("dPt3Nch", "dPt3Nch", 90, 0, 4500);
    dPt4Nch = new TProfile("dPt4Nch", "dPt4Nch", 90, 0, 4500);
    dPtNchNB = new TProfile("dPtNchNB", "dPtNchNB", 90, 0, 4500);
    fPtNchUCC = new TH2F("fPtNchUCC", "fPtNchUCC", 700, 1000, 4500, 500, 0.2, 5.0);
    TestPtCtr = new TProfile("TestPtCtr", "TestPtCtr", 200, 0, 100);
    TestPt2Ctr = new TProfile("TestPt2Ctr", "TestPt2Ctr", 200, 0, 100);
    TestPt3Ctr = new TProfile("TestPt3Ctr", "TestPt3Ctr", 200, 0, 100);
    TestPt4Ctr = new TProfile("TestPt4Ctr", "TestPt4Ctr", 200, 0, 100);
    TestNchCtr = new TProfile("TestNchCtr", "TestNchCtr", 200, 0, 100);
    TestNchSelectedCtr = new TProfile("TestNchSelectedCtr", "TestNchSelectedCtr", 200, 0, 100);

    fBstList = new TList();
    dPtNch_0 = new TProfile("dPtNch_0", "dPtNch_0", 90, 0, 4500);
    dPtNch_1 = new TProfile("dPtNch_1", "dPtNch_1", 90, 0, 4500);
    dPtNch_2 = new TProfile("dPtNch_2", "dPtNch_2", 90, 0, 4500);
    dPtNch_3 = new TProfile("dPtNch_3", "dPtNch_3", 90, 0, 4500);
    dPtNch_4 = new TProfile("dPtNch_4", "dPtNch_4", 90, 0, 4500);
    dPtNch_5 = new TProfile("dPtNch_5", "dPtNch_5", 90, 0, 4500);
    dPtNch_6 = new TProfile("dPtNch_6", "dPtNch_6", 90, 0, 4500);
    dPtNch_7 = new TProfile("dPtNch_7", "dPtNch_7", 90, 0, 4500);
    dPtNch_8 = new TProfile("dPtNch_8", "dPtNch_8", 90, 0, 4500);
    dPtNch_9 = new TProfile("dPtNch_9", "dPtNch_9", 90, 0, 4500);
    
    dPt2Nch_0 = new TProfile("dPt2Nch_0", "dPt2Nch_0", 90, 0, 4500);
    dPt2Nch_1 = new TProfile("dPt2Nch_1", "dPt2Nch_1", 90, 0, 4500);
    dPt2Nch_2 = new TProfile("dPt2Nch_2", "dPt2Nch_2", 90, 0, 4500);
    dPt2Nch_3 = new TProfile("dPt2Nch_3", "dPt2Nch_3", 90, 0, 4500);
    dPt2Nch_4 = new TProfile("dPt2Nch_4", "dPt2Nch_4", 90, 0, 4500);
    dPt2Nch_5 = new TProfile("dPt2Nch_5", "dPt2Nch_5", 90, 0, 4500);
    dPt2Nch_6 = new TProfile("dPt2Nch_6", "dPt2Nch_6", 90, 0, 4500);
    dPt2Nch_7 = new TProfile("dPt2Nch_7", "dPt2Nch_7", 90, 0, 4500);
    dPt2Nch_8 = new TProfile("dPt2Nch_8", "dPt2Nch_8", 90, 0, 4500);
    dPt2Nch_9 = new TProfile("dPt2Nch_9", "dPt2Nch_9", 90, 0, 4500);
  
    dPt3Nch_0 = new TProfile("dPt3Nch_0", "dPt3Nch_0", 90, 0, 4500);
    dPt3Nch_1 = new TProfile("dPt3Nch_1", "dPt3Nch_1", 90, 0, 4500);
    dPt3Nch_2 = new TProfile("dPt3Nch_2", "dPt3Nch_2", 90, 0, 4500);
    dPt3Nch_3 = new TProfile("dPt3Nch_3", "dPt3Nch_3", 90, 0, 4500);
    dPt3Nch_4 = new TProfile("dPt3Nch_4", "dPt3Nch_4", 90, 0, 4500);
    dPt3Nch_5 = new TProfile("dPt3Nch_5", "dPt3Nch_5", 90, 0, 4500);
    dPt3Nch_6 = new TProfile("dPt3Nch_6", "dPt3Nch_6", 90, 0, 4500);
    dPt3Nch_7 = new TProfile("dPt3Nch_7", "dPt3Nch_7", 90, 0, 4500);
    dPt3Nch_8 = new TProfile("dPt3Nch_8", "dPt3Nch_8", 90, 0, 4500);
    dPt3Nch_9 = new TProfile("dPt3Nch_9", "dPt3Nch_9", 90, 0, 4500);

    dPt4Nch_0 = new TProfile("dPt4Nch_0", "dPt4Nch_0", 90, 0, 4500);
    dPt4Nch_1 = new TProfile("dPt4Nch_1", "dPt4Nch_1", 90, 0, 4500);
    dPt4Nch_2 = new TProfile("dPt4Nch_2", "dPt4Nch_2", 90, 0, 4500);
    dPt4Nch_3 = new TProfile("dPt4Nch_3", "dPt4Nch_3", 90, 0, 4500);
    dPt4Nch_4 = new TProfile("dPt4Nch_4", "dPt4Nch_4", 90, 0, 4500);
    dPt4Nch_5 = new TProfile("dPt4Nch_5", "dPt4Nch_5", 90, 0, 4500);
    dPt4Nch_6 = new TProfile("dPt4Nch_6", "dPt4Nch_6", 90, 0, 4500);
    dPt4Nch_7 = new TProfile("dPt4Nch_7", "dPt4Nch_7", 90, 0, 4500);
    dPt4Nch_8 = new TProfile("dPt4Nch_8", "dPt4Nch_8", 90, 0, 4500);
    dPt4Nch_9 = new TProfile("dPt4Nch_9", "dPt4Nch_9", 90, 0, 4500);

    TestPtCtr_0 = new TProfile("TestPtCtr_0", "TestPtCtr_0", 200, 0, 100);
    TestPtCtr_1 = new TProfile("TestPtCtr_1", "TestPtCtr_1", 200, 0, 100);
    TestPtCtr_2 = new TProfile("TestPtCtr_2", "TestPtCtr_2", 200, 0, 100);
    TestPtCtr_3 = new TProfile("TestPtCtr_3", "TestPtCtr_3", 200, 0, 100);
    TestPtCtr_4 = new TProfile("TestPtCtr_4", "TestPtCtr_4", 200, 0, 100);
    TestPtCtr_5 = new TProfile("TestPtCtr_5", "TestPtCtr_5", 200, 0, 100);
    TestPtCtr_6 = new TProfile("TestPtCtr_6", "TestPtCtr_6", 200, 0, 100);
    TestPtCtr_7 = new TProfile("TestPtCtr_7", "TestPtCtr_7", 200, 0, 100);
    TestPtCtr_8 = new TProfile("TestPtCtr_8", "TestPtCtr_8", 200, 0, 100);
    TestPtCtr_9 = new TProfile("TestPtCtr_9", "TestPtCtr_9", 200, 0, 100);

    TestPt2Ctr_0 = new TProfile("TestPt2Ctr_0", "TestPt2Ctr_0", 200, 0, 100);
    TestPt2Ctr_1 = new TProfile("TestPt2Ctr_1", "TestPt2Ctr_1", 200, 0, 100);
    TestPt2Ctr_2 = new TProfile("TestPt2Ctr_2", "TestPt2Ctr_2", 200, 0, 100);
    TestPt2Ctr_3 = new TProfile("TestPt2Ctr_3", "TestPt2Ctr_3", 200, 0, 100);
    TestPt2Ctr_4 = new TProfile("TestPt2Ctr_4", "TestPt2Ctr_4", 200, 0, 100);
    TestPt2Ctr_5 = new TProfile("TestPt2Ctr_5", "TestPt2Ctr_5", 200, 0, 100);
    TestPt2Ctr_6 = new TProfile("TestPt2Ctr_6", "TestPt2Ctr_6", 200, 0, 100);
    TestPt2Ctr_7 = new TProfile("TestPt2Ctr_7", "TestPt2Ctr_7", 200, 0, 100);
    TestPt2Ctr_8 = new TProfile("TestPt2Ctr_8", "TestPt2Ctr_8", 200, 0, 100);
    TestPt2Ctr_9 = new TProfile("TestPt2Ctr_9", "TestPt2Ctr_9", 200, 0, 100);

    TestPt3Ctr_0 = new TProfile("TestPt3Ctr_0", "TestPt3Ctr_0", 200, 0, 100);
    TestPt3Ctr_1 = new TProfile("TestPt3Ctr_1", "TestPt3Ctr_1", 200, 0, 100);
    TestPt3Ctr_2 = new TProfile("TestPt3Ctr_2", "TestPt3Ctr_2", 200, 0, 100);
    TestPt3Ctr_3 = new TProfile("TestPt3Ctr_3", "TestPt3Ctr_3", 200, 0, 100);
    TestPt3Ctr_4 = new TProfile("TestPt3Ctr_4", "TestPt3Ctr_4", 200, 0, 100);
    TestPt3Ctr_5 = new TProfile("TestPt3Ctr_5", "TestPt3Ctr_5", 200, 0, 100);
    TestPt3Ctr_6 = new TProfile("TestPt3Ctr_6", "TestPt3Ctr_6", 200, 0, 100);
    TestPt3Ctr_7 = new TProfile("TestPt3Ctr_7", "TestPt3Ctr_7", 200, 0, 100);
    TestPt3Ctr_8 = new TProfile("TestPt3Ctr_8", "TestPt3Ctr_8", 200, 0, 100);
    TestPt3Ctr_9 = new TProfile("TestPt3Ctr_9", "TestPt3Ctr_9", 200, 0, 100);

    TestPt4Ctr_0 = new TProfile("TestPt4Ctr_0", "TestPt4Ctr_0", 200, 0, 100);
    TestPt4Ctr_1 = new TProfile("TestPt4Ctr_1", "TestPt4Ctr_1", 200, 0, 100);
    TestPt4Ctr_2 = new TProfile("TestPt4Ctr_2", "TestPt4Ctr_2", 200, 0, 100);
    TestPt4Ctr_3 = new TProfile("TestPt4Ctr_3", "TestPt4Ctr_3", 200, 0, 100);
    TestPt4Ctr_4 = new TProfile("TestPt4Ctr_4", "TestPt4Ctr_4", 200, 0, 100);
    TestPt4Ctr_5 = new TProfile("TestPt4Ctr_5", "TestPt4Ctr_5", 200, 0, 100);
    TestPt4Ctr_6 = new TProfile("TestPt4Ctr_6", "TestPt4Ctr_6", 200, 0, 100);
    TestPt4Ctr_7 = new TProfile("TestPt4Ctr_7", "TestPt4Ctr_7", 200, 0, 100);
    TestPt4Ctr_8 = new TProfile("TestPt4Ctr_8", "TestPt4Ctr_8", 200, 0, 100);
    TestPt4Ctr_9 = new TProfile("TestPt4Ctr_9", "TestPt4Ctr_9", 200, 0, 100);

    TestNchCtr_0 = new TProfile("TestNchCtr_0", "TestNchCtr_0", 200, 0, 100);
    TestNchCtr_1 = new TProfile("TestNchCtr_1", "TestNchCtr_1", 200, 0, 100);
    TestNchCtr_2 = new TProfile("TestNchCtr_2", "TestNchCtr_2", 200, 0, 100);
    TestNchCtr_3 = new TProfile("TestNchCtr_3", "TestNchCtr_3", 200, 0, 100);
    TestNchCtr_4 = new TProfile("TestNchCtr_4", "TestNchCtr_4", 200, 0, 100);
    TestNchCtr_5 = new TProfile("TestNchCtr_5", "TestNchCtr_5", 200, 0, 100);
    TestNchCtr_6 = new TProfile("TestNchCtr_6", "TestNchCtr_6", 200, 0, 100);
    TestNchCtr_7 = new TProfile("TestNchCtr_7", "TestNchCtr_7", 200, 0, 100);
    TestNchCtr_8 = new TProfile("TestNchCtr_8", "TestNchCtr_8", 200, 0, 100);
    TestNchCtr_9 = new TProfile("TestNchCtr_9", "TestNchCtr_9", 200, 0, 100);

    TestNchSelectedCtr_0 = new TProfile("TestNchSelectedCtr_0", "TestNchSelectedCtr_0", 200, 0, 100);
    TestNchSelectedCtr_1 = new TProfile("TestNchSelectedCtr_1", "TestNchSelectedCtr_1", 200, 0, 100);
    TestNchSelectedCtr_2 = new TProfile("TestNchSelectedCtr_2", "TestNchSelectedCtr_2", 200, 0, 100);
    TestNchSelectedCtr_3 = new TProfile("TestNchSelectedCtr_3", "TestNchSelectedCtr_3", 200, 0, 100);
    TestNchSelectedCtr_4 = new TProfile("TestNchSelectedCtr_4", "TestNchSelectedCtr_4", 200, 0, 100);
    TestNchSelectedCtr_5 = new TProfile("TestNchSelectedCtr_5", "TestNchSelectedCtr_5", 200, 0, 100);
    TestNchSelectedCtr_6 = new TProfile("TestNchSelectedCtr_6", "TestNchSelectedCtr_6", 200, 0, 100);
    TestNchSelectedCtr_7 = new TProfile("TestNchSelectedCtr_7", "TestNchSelectedCtr_7", 200, 0, 100);
    TestNchSelectedCtr_8 = new TProfile("TestNchSelectedCtr_8", "TestNchSelectedCtr_8", 200, 0, 100);
    TestNchSelectedCtr_9 = new TProfile("TestNchSelectedCtr_9", "TestNchSelectedCtr_9", 200, 0, 100);
     // create your histogram
    fBstList->Add(dPtNch_0);
    fBstList->Add(dPtNch_1);
    fBstList->Add(dPtNch_2);
    fBstList->Add(dPtNch_3);
    fBstList->Add(dPtNch_4);
    fBstList->Add(dPtNch_5);
    fBstList->Add(dPtNch_6);
    fBstList->Add(dPtNch_7);
    fBstList->Add(dPtNch_8);
    fBstList->Add(dPtNch_9);
    fBstList->Add(dPt2Nch_0);
    fBstList->Add(dPt2Nch_1);
    fBstList->Add(dPt2Nch_2);
    fBstList->Add(dPt2Nch_3);
    fBstList->Add(dPt2Nch_4);
    fBstList->Add(dPt2Nch_5);
    fBstList->Add(dPt2Nch_6);
    fBstList->Add(dPt2Nch_7);
    fBstList->Add(dPt2Nch_8);
    fBstList->Add(dPt2Nch_9);
    fBstList->Add(dPt3Nch_0);
    fBstList->Add(dPt3Nch_1);
    fBstList->Add(dPt3Nch_2);
    fBstList->Add(dPt3Nch_3);
    fBstList->Add(dPt3Nch_4);
    fBstList->Add(dPt3Nch_5);
    fBstList->Add(dPt3Nch_6);
    fBstList->Add(dPt3Nch_7);
    fBstList->Add(dPt3Nch_8);
    fBstList->Add(dPt3Nch_9);
    fBstList->Add(dPt4Nch_0);
    fBstList->Add(dPt4Nch_1);
    fBstList->Add(dPt4Nch_2);
    fBstList->Add(dPt4Nch_3);
    fBstList->Add(dPt4Nch_4);
    fBstList->Add(dPt4Nch_5);
    fBstList->Add(dPt4Nch_6);
    fBstList->Add(dPt4Nch_7);
    fBstList->Add(dPt4Nch_8);
    fBstList->Add(dPt4Nch_9);
    fBstList->Add(TestPtCtr_0);
    fBstList->Add(TestPtCtr_1);
    fBstList->Add(TestPtCtr_2);
    fBstList->Add(TestPtCtr_3);
    fBstList->Add(TestPtCtr_4);
    fBstList->Add(TestPtCtr_5);
    fBstList->Add(TestPtCtr_6);
    fBstList->Add(TestPtCtr_7);
    fBstList->Add(TestPtCtr_8);
    fBstList->Add(TestPtCtr_9);
    fBstList->Add(TestPt2Ctr_0);
    fBstList->Add(TestPt2Ctr_1);
    fBstList->Add(TestPt2Ctr_2);
    fBstList->Add(TestPt2Ctr_3);
    fBstList->Add(TestPt2Ctr_4);
    fBstList->Add(TestPt2Ctr_5);
    fBstList->Add(TestPt2Ctr_6);
    fBstList->Add(TestPt2Ctr_7);
    fBstList->Add(TestPt2Ctr_8);
    fBstList->Add(TestPt2Ctr_9);
    fBstList->Add(TestPt3Ctr_0);
    fBstList->Add(TestPt3Ctr_1);
    fBstList->Add(TestPt3Ctr_2);
    fBstList->Add(TestPt3Ctr_3);
    fBstList->Add(TestPt3Ctr_4);
    fBstList->Add(TestPt3Ctr_5);
    fBstList->Add(TestPt3Ctr_6);
    fBstList->Add(TestPt3Ctr_7);
    fBstList->Add(TestPt3Ctr_8);
    fBstList->Add(TestPt3Ctr_9);
    fBstList->Add(TestPt4Ctr_0);
    fBstList->Add(TestPt4Ctr_1);
    fBstList->Add(TestPt4Ctr_2);
    fBstList->Add(TestPt4Ctr_3);
    fBstList->Add(TestPt4Ctr_4);
    fBstList->Add(TestPt4Ctr_5);
    fBstList->Add(TestPt4Ctr_6);
    fBstList->Add(TestPt4Ctr_7);
    fBstList->Add(TestPt4Ctr_8);
    fBstList->Add(TestPt4Ctr_9);
    fBstList->Add(TestNchCtr_0);
    fBstList->Add(TestNchCtr_1);
    fBstList->Add(TestNchCtr_2);
    fBstList->Add(TestNchCtr_3);
    fBstList->Add(TestNchCtr_4);
    fBstList->Add(TestNchCtr_5);
    fBstList->Add(TestNchCtr_6);
    fBstList->Add(TestNchCtr_7);
    fBstList->Add(TestNchCtr_8);
    fBstList->Add(TestNchCtr_9);
    fBstList->Add(TestNchSelectedCtr_0);
    fBstList->Add(TestNchSelectedCtr_1);
    fBstList->Add(TestNchSelectedCtr_2);
    fBstList->Add(TestNchSelectedCtr_3);
    fBstList->Add(TestNchSelectedCtr_4);
    fBstList->Add(TestNchSelectedCtr_5);
    fBstList->Add(TestNchSelectedCtr_6);
    fBstList->Add(TestNchSelectedCtr_7);
    fBstList->Add(TestNchSelectedCtr_8);
    fBstList->Add(TestNchSelectedCtr_9);


    fOutputList->Add(fBstList);
    fOutputList->Add(fTestNonWeight);
    fOutputList->Add(fNchDistri);
    fOutputList->Add(fPtQA);
    fOutputList->Add(fDCAxyQA);
    fOutputList->Add(fDCAzQA);
    fOutputList->Add(fetaQA);
    fOutputList->Add(fTPCclsQA);
    fOutputList->Add(fDCAxyPt);
    fOutputList->Add(fDCAzPt);
    fOutputList->Add(fPtNch);
    fOutputList->Add(dPtNch);
    fOutputList->Add(dPt2Nch);
    fOutputList->Add(dPt3Nch);
    fOutputList->Add(dPt4Nch);
    fOutputList->Add(dPtNchNB);
    fOutputList->Add(fPtNchUCC);
    fOutputList->Add(TestPtCtr);
    fOutputList->Add(TestPt2Ctr);
    fOutputList->Add(TestPt3Ctr);
    fOutputList->Add(TestPt4Ctr);
    fOutputList->Add(TestNchCtr);
    fOutputList->Add(TestNchSelectedCtr);
         // don't forget to add it to the list! the list will be written to file, so if you want
                                        // your histogram in the output file, add it to the list!
    
    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisPtN::UserExec(Option_t *)
{
    
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
                                                        // there's another event format (ESD) which works in a similar way
    // cout<< "event inputted"<<endl;                                                    // but is more cpu/memory unfriendly. for now, we'll stick with aod's
    if(!fAOD) return;
    // cout<< "event is AOD"<<endl;
    

    UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isTrigselected = false;
    if(fTrigger==0){
      isTrigselected = fSelectMask&(AliVEvent::kINT7+AliVEvent::kMB);
    } else if(fTrigger==1){
      isTrigselected = fSelectMask&(AliVEvent::kINT7+AliVEvent::kMB+AliVEvent::kCentral);
    }
    if(isTrigselected == false) return;
    
    if(!fEventCuts.AcceptEvent(fAOD)) return;                                  // if the pointer to the event is empty (getting it failed) skip this event
        // example part: i'll show how to loop over the tracks in an event 
        // and extract some information from them which we'll store in a histogram
    // cout<< "event selected"<<endl;
    Float_t vertexZ = fAOD->GetPrimaryVertex()->GetZ();
    if(vertexZ<=(-fVz) || vertexZ>=fVz) return;             //filter Vz range

    Float_t centrality(0);
    AliMultSelection* multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if(multSelection) centrality = multSelection->GetMultiplicityPercentile(ctrType);
    //if(centrality>90) return;

    //initilize the weight list
    fWeightsListNUE = (TList*)GetInputData(1);
    if (fSysflg==0) {
      fWeightNUE = (TH1D*)fWeightsListNUE->FindObject(Form("EffRescaled_Cent0"));
    } else {
      fWeightNUE = (TH1D*)fWeightsListNUE->FindObject(Form("EffRescaled_Cent0_SystFlag%i_", fSysflg));
    }
    

    Float_t wtE;

    Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
    Double_t sump = 0, sumw = 0, sump2 = 0, sumw2 = 0, sump3 = 0, sumw3 = 0, sump4 = 0, sumw4 = 0;
    Int_t nTrackSelected = 0;
                                            // how many tracks are selected
    for(Int_t i(0); i < iTracks; i++) {
                                                            // loop ove rall these tracks
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track || !track->TestFilterBit(filterBit)) continue; // filterbit if we failed, skip this track
        if(track->Pt()<=fMinPt || track->Pt()>=fMaxPt) continue; //filter Pt range                           
        if(track->Eta()<=(-0.8) || track->Eta()>=0.8) continue; //filter Eta range
        if(track->GetTPCNcls()<=fTPCcls) continue; //filter TPC cluster numbers
        Double_t vtxXYZ[3], trXYZ[3];
        track->GetXYZ(trXYZ);
        fAOD->GetPrimaryVertex()->GetXYZ(vtxXYZ);
        trXYZ[2] -= vtxXYZ[2];
        if(TMath::Abs(trXYZ[2]) > fDCAz) continue;
        trXYZ[0] -= vtxXYZ[0];
        trXYZ[1] -= vtxXYZ[1];
        Double_t trDcaxy = TMath::Sqrt(trXYZ[0]*trXYZ[0]+trXYZ[1]*trXYZ[1]);
        Double_t cutDcaxy = 0.0026+0.005/TMath::Power(track->Pt(),1.01);
        if(trDcaxy > fDCAxy*cutDcaxy) continue;

        //all the selection are before here
        ++nTrackSelected;
        wtE = GetWeightNUE(track->Pt());
        fTestNonWeight->Fill(track->Phi());
        fPtQA->Fill(track->Pt());
        fDCAxyQA->Fill(trDcaxy);
        fDCAzQA->Fill(trXYZ[2]);
        fetaQA->Fill(track->Eta());
        fTPCclsQA->Fill(track->GetTPCNcls());
        fDCAxyPt->Fill(trDcaxy, track->Pt());
        fDCAzPt->Fill(trXYZ[2], track->Pt());


          sump += wtE*track->Pt();
          sump2 += wtE*wtE*track->Pt()*track->Pt();
          sump3 += wtE*wtE*wtE*track->Pt()*track->Pt()*track->Pt();
          sump4 += wtE*wtE*wtE*wtE*track->Pt()*track->Pt()*track->Pt()*track->Pt();

          sumw += wtE;
          sumw2 += wtE*wtE;
          sumw3 += wtE*wtE*wtE;
          sumw4 += wtE*wtE*wtE*wtE;
    
    }
    
    Int_t rd = int(floor(radm.Rndm()*10));

    //if(M<100) return; //we need enough tracks in one event
    if(sumw!=0 && (sumw*sumw-sumw2)!=0 && (sumw*sumw*sumw-3*sumw2*sumw+2*sumw3)!=0 && (pow(sumw,4)-6*sumw2*sumw*sumw+3*sumw2*sumw2+8*sumw3*sumw-6*sumw4)!=0) {
      Float_t pt = sump/sumw;
      fPtNch->Fill(nTrackSelected,pt);
      dPtNch->Fill(nTrackSelected,pt);

      Float_t pt2 = (sump*sump-sump2)/(sumw*sumw-sumw2);
      dPt2Nch->Fill(nTrackSelected,pt2);

      Float_t pt3 = (sump*sump*sump-3*sump2*sump+2*sump3)/(sumw*sumw*sumw-3*sumw2*sumw+2*sumw3);
      dPt3Nch->Fill(nTrackSelected,pt3);

      Float_t pt4 = (pow(sump,4)-6*sump2*sump*sump+3*sump2*sump2+8*sump3*sump-6*sump4)/(pow(sumw,4)-6*sumw2*sumw*sumw+3*sumw2*sumw2+8*sumw3*sumw-6*sumw4);
      dPt4Nch->Fill(nTrackSelected,pt4);

      if(centrality<1) fPtNchUCC->Fill(nTrackSelected,pt);
      TestPtCtr->Fill(centrality,pt);
      TestPt2Ctr->Fill(centrality,pt2);
      TestPt3Ctr->Fill(centrality,pt3);
      TestPt4Ctr->Fill(centrality,pt4);

    //here we fill the boostrap profiles
      
      switch (rd){
        case 0:
          dPtNch_0->Fill(nTrackSelected,pt);
          dPt2Nch_0->Fill(nTrackSelected,pt2);
          dPt3Nch_0->Fill(nTrackSelected,pt3);
          dPt4Nch_0->Fill(nTrackSelected,pt4);
          TestPtCtr_0->Fill(centrality,pt);
          TestPt2Ctr_0->Fill(centrality,pt2);
          TestPt3Ctr_0->Fill(centrality,pt3);
          TestPt4Ctr_0->Fill(centrality,pt4);
          break;
        case 1:
          dPtNch_1->Fill(nTrackSelected,pt);
          dPt2Nch_1->Fill(nTrackSelected,pt2);
          dPt3Nch_1->Fill(nTrackSelected,pt3);
          dPt4Nch_1->Fill(nTrackSelected,pt4);
          TestPtCtr_1->Fill(centrality,pt);
          TestPt2Ctr_1->Fill(centrality,pt2);
          TestPt3Ctr_1->Fill(centrality,pt3);
          TestPt4Ctr_1->Fill(centrality,pt4);
          break;
        case 2:
          dPtNch_2->Fill(nTrackSelected,pt);
          dPt2Nch_2->Fill(nTrackSelected,pt2);
          dPt3Nch_2->Fill(nTrackSelected,pt3);
          dPt4Nch_2->Fill(nTrackSelected,pt4);
          TestPtCtr_2->Fill(centrality,pt);
          TestPt2Ctr_2->Fill(centrality,pt2);
          TestPt3Ctr_2->Fill(centrality,pt3);
          TestPt4Ctr_2->Fill(centrality,pt4);
          break;
        case 3:
          dPtNch_3->Fill(nTrackSelected,pt);
          dPt2Nch_3->Fill(nTrackSelected,pt2);
          dPt3Nch_3->Fill(nTrackSelected,pt3);
          dPt4Nch_3->Fill(nTrackSelected,pt4);
          TestPtCtr_3->Fill(centrality,pt);
          TestPt2Ctr_3->Fill(centrality,pt2);
          TestPt3Ctr_3->Fill(centrality,pt3);
          TestPt4Ctr_3->Fill(centrality,pt4);
          break;
        case 4:
          dPtNch_4->Fill(nTrackSelected,pt);
          dPt2Nch_4->Fill(nTrackSelected,pt2);
          dPt3Nch_4->Fill(nTrackSelected,pt3);
          dPt4Nch_4->Fill(nTrackSelected,pt4);
          TestPtCtr_4->Fill(centrality,pt);
          TestPt2Ctr_4->Fill(centrality,pt2);
          TestPt3Ctr_4->Fill(centrality,pt3);
          TestPt4Ctr_4->Fill(centrality,pt4);
          break;
        case 5:
          dPtNch_5->Fill(nTrackSelected,pt);
          dPt2Nch_5->Fill(nTrackSelected,pt2);
          dPt3Nch_5->Fill(nTrackSelected,pt3);
          dPt4Nch_5->Fill(nTrackSelected,pt4);
          TestPtCtr_5->Fill(centrality,pt);
          TestPt2Ctr_5->Fill(centrality,pt2);
          TestPt3Ctr_5->Fill(centrality,pt3);
          TestPt4Ctr_5->Fill(centrality,pt4);
          break;
        case 6:
          dPtNch_6->Fill(nTrackSelected,pt);
          dPt2Nch_6->Fill(nTrackSelected,pt2);
          dPt3Nch_6->Fill(nTrackSelected,pt3);
          dPt4Nch_6->Fill(nTrackSelected,pt4);
          TestPtCtr_6->Fill(centrality,pt);
          TestPt2Ctr_6->Fill(centrality,pt2);
          TestPt3Ctr_6->Fill(centrality,pt3);
          TestPt4Ctr_6->Fill(centrality,pt4);
          break;
        case 7:
          dPtNch_7->Fill(nTrackSelected,pt);
          dPt2Nch_7->Fill(nTrackSelected,pt2);
          dPt3Nch_7->Fill(nTrackSelected,pt3);
          dPt4Nch_7->Fill(nTrackSelected,pt4);
          TestPtCtr_7->Fill(centrality,pt);
          TestPt2Ctr_7->Fill(centrality,pt2);
          TestPt3Ctr_7->Fill(centrality,pt3);
          TestPt4Ctr_7->Fill(centrality,pt4);
          break;
        case 8:
          dPtNch_8->Fill(nTrackSelected,pt);
          dPt2Nch_8->Fill(nTrackSelected,pt2);
          dPt3Nch_8->Fill(nTrackSelected,pt3);
          dPt4Nch_8->Fill(nTrackSelected,pt4);
          TestPtCtr_8->Fill(centrality,pt);
          TestPt2Ctr_8->Fill(centrality,pt2);
          TestPt3Ctr_8->Fill(centrality,pt3);
          TestPt4Ctr_8->Fill(centrality,pt4);
          break;
        case 9:
          dPtNch_9->Fill(nTrackSelected,pt);
          dPt2Nch_9->Fill(nTrackSelected,pt2);
          dPt3Nch_9->Fill(nTrackSelected,pt3);
          dPt4Nch_9->Fill(nTrackSelected,pt4);
          TestPtCtr_9->Fill(centrality,pt);
          TestPt2Ctr_9->Fill(centrality,pt2);
          TestPt3Ctr_9->Fill(centrality,pt3);
          TestPt4Ctr_9->Fill(centrality,pt4);
          break;
      }
      TestNchCtr->Fill(centrality,iTracks);
      TestNchSelectedCtr->Fill(centrality,nTrackSelected);

      switch (rd){
        case 0:
          TestNchCtr_0->Fill(centrality,iTracks);
          TestNchSelectedCtr_0->Fill(centrality,nTrackSelected);
          break;
        case 1:
          TestNchCtr_1->Fill(centrality,iTracks);
          TestNchSelectedCtr_1->Fill(centrality,nTrackSelected);
          break;
        case 2:
          TestNchCtr_2->Fill(centrality,iTracks);
          TestNchSelectedCtr_2->Fill(centrality,nTrackSelected);
          break;
        case 3:
          TestNchCtr_3->Fill(centrality,iTracks);
          TestNchSelectedCtr_3->Fill(centrality,nTrackSelected);
          break;
        case 4:
          TestNchCtr_4->Fill(centrality,iTracks);
          TestNchSelectedCtr_4->Fill(centrality,nTrackSelected);
          break;
        case 5:
          TestNchCtr_5->Fill(centrality,iTracks);
          TestNchSelectedCtr_5->Fill(centrality,nTrackSelected);
          break;
        case 6:
          TestNchCtr_6->Fill(centrality,iTracks);
          TestNchSelectedCtr_6->Fill(centrality,nTrackSelected);
          break;
        case 7:
          TestNchCtr_7->Fill(centrality,iTracks);
          TestNchSelectedCtr_7->Fill(centrality,nTrackSelected);
          break;
        case 8:
          TestNchCtr_8->Fill(centrality,iTracks);
          TestNchSelectedCtr_8->Fill(centrality,nTrackSelected);
          break;
        case 9:
          TestNchCtr_9->Fill(centrality,iTracks);
          TestNchSelectedCtr_9->Fill(centrality,nTrackSelected);
          break;
      }

    }
    
    if(sumw!=0) {
       Float_t ptnb = sump/sumw;
       dPtNchNB->Fill(nTrackSelected,ptnb);
    }

    fNchDistri->Fill(nTrackSelected);
    
    

    
                                                        // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}
//____________________________________________________________________________
// double AliAnalysisPtN::GetWeightNUA(double phi, double eta, double vz) {
//   double weight = fWeightNUA->GetBinContent(fWeightNUA->GetXaxis()->FindBin(phi),
//       fWeightNUA->GetYaxis()->FindBin(eta),
//       fWeightNUA->GetZaxis()->FindBin(vz));
//   return weight;
// }
//_____________________________________________________________________________
double AliAnalysisPtN::GetWeightNUE(double pt)
{
  double binPt = fWeightNUE->GetXaxis()->FindBin(pt);
  double eff = fWeightNUE->GetBinContent(binPt);
  double error = fWeightNUE->GetBinError(binPt);
  double weight = 1;
  //..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error


  if((eff < 0.03) || ((error/eff) > 0.1)) error = 0.00001;
  if((eff < 0.03)) return 1;

  TRandom3 r(0);
  double efficiency = 0;
  efficiency = r.Gaus(eff, error);
  weight = 1./efficiency; 
  return weight;
}
//_____________________________________________________________________________
void AliAnalysisPtN::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
