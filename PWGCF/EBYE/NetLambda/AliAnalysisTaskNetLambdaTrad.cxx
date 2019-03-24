
// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Naomi Umaka Apr 2018
// Updated Mar20


#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliEventCuts.h"
#include "AliExternalTrackParam.h"
#include "AliAnalysisFilter.h"
#include "AliVMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliStack.h"
#include "AliPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliV0vertexer.h"
#include "AliESDv0Cuts.h"
#include "AliMultSelection.h"
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH3F.h"
#include "THnSparse.h"


#include "AliAnalysisTaskNetLambdaTrad.h"

ClassImp(AliAnalysisTaskNetLambdaTrad)

//-----------------------------------------------------------------------------
AliAnalysisTaskNetLambdaTrad::AliAnalysisTaskNetLambdaTrad(const char* name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fPIDResponse(0x0),
fEventCuts(0),
fListHist(0x0),
fHistEventCounter(0x0),
fHistCentrality(0x0),

f2fHistGenCentVsPtLambda(0x0),
f2fHistGenCentVsPtAntiLambda(0x0),
f2fHistXiPlus(0x0),
f2fHistXiMinus(0x0),

f2fHistRecPrimariesCentVsPtLambdaFour(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaFour(0x0),
f2fHistRecPrimariesCentVsPtLambdaThree(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaThree(0x0),
f2fHistRecPrimariesCentVsPtLambdaTwo(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaTwo(0x0),
f2fHistRecPrimariesCentVsPtLambdaOne(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaOne(0x0),

f2fHistLambdaMisIdFour(0x0),
f2fHistAntiLambdaMisIdFour(0x0),
f2fHistLambdaMisIdThree(0x0),
f2fHistAntiLambdaMisIdThree(0x0),
f2fHistLambdaMisIdTwo(0x0),
f2fHistAntiLambdaMisIdTwo(0x0),
f2fHistLambdaMisIdOne(0x0),
f2fHistAntiLambdaMisIdOne(0x0),

f3fHistCentInvMassVsPtLambdaRecFour(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecFour(0x0),
f3fHistCentInvMassVsPtLambdaRecThree(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecThree(0x0),
f3fHistCentInvMassVsPtLambdaRecTwo(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecTwo(0x0),
f3fHistCentInvMassVsPtLambdaRecOne(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecOne(0x0),

f3fHistLambdafromXiFour(0x0),
f3fHistAntiLambdafromXiFour(0x0),
f3fHistLambdafromXiThree(0x0),
f3fHistAntiLambdafromXiThree(0x0),
f3fHistLambdafromXiTwo(0x0),
f3fHistAntiLambdafromXiTwo(0x0),
f3fHistLambdafromXiOne(0x0),
f3fHistAntiLambdafromXiOne(0x0),

//SIG3
f2fHistRecPrimariesCentVsPtLambdaFourSigthree(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree(0x0),
f2fHistRecPrimariesCentVsPtLambdaThreeSigthree(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaThreeSigthree(0x0),
f2fHistRecPrimariesCentVsPtLambdaTwoSigthree(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaTwoSigthree(0x0),
f2fHistRecPrimariesCentVsPtLambdaOneSigthree(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaOneSigthree(0x0),

f2fHistLambdaMisIdFourSigthree(0x0),
f2fHistAntiLambdaMisIdFourSigthree(0x0),
f2fHistLambdaMisIdThreeSigthree(0x0),
f2fHistAntiLambdaMisIdThreeSigthree(0x0),
f2fHistLambdaMisIdTwoSigthree(0x0),
f2fHistAntiLambdaMisIdTwoSigthree(0x0),
f2fHistLambdaMisIdOneSigthree(0x0),
f2fHistAntiLambdaMisIdOneSigthree(0x0),

f3fHistCentInvMassVsPtLambdaRecFourSigthree(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree(0x0),
f3fHistCentInvMassVsPtLambdaRecThreeSigthree(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecThreeSigthree(0x0),
f3fHistCentInvMassVsPtLambdaRecTwoSigthree(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecTwoSigthree(0x0),
f3fHistCentInvMassVsPtLambdaRecOneSigthree(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecOneSigthree(0x0),

f3fHistLambdafromXiFourSigthree(0x0),
f3fHistAntiLambdafromXiFourSigthree(0x0),
f3fHistLambdafromXiThreeSigthree(0x0),
f3fHistAntiLambdafromXiThreeSigthree(0x0),
f3fHistLambdafromXiTwoSigthree(0x0),
f3fHistAntiLambdafromXiTwoSigthree(0x0),
f3fHistLambdafromXiOneSigthree(0x0),
f3fHistAntiLambdafromXiOneSigthree(0x0),

f2fHistLRecstat(0x0),
f2fHistARecstat(0x0),
f2fHistLGenstat(0x0),
f2fHistAGenstat(0x0),

fCentrality(-1),
fTreeVariablePID(-1),

fNptBins(31),
fIsMC(kTRUE),
fEvSel(AliVEvent::kINT7),


fPtBinNplusNminusChTagFour(NULL),
fPtBinNplusNminusChTagThree(NULL),
fPtBinNplusNminusChTagTwo(NULL),
fPtBinNplusNminusChTagOne(NULL),

fPtBinNplusNminusChTagFourSigThree(NULL),
fPtBinNplusNminusChTagThreeSigThree(NULL),
fPtBinNplusNminusChTagTwoSigThree(NULL),
fPtBinNplusNminusChTagOneSigThree(NULL),

fPtBinNplusNminusChTruth(NULL)

{
    Info("AliAnalysisTaskNetLambdaTrad","Calling Constructor");
    
    DefineInput(0,TChain::Class());
    DefineOutput(1,TList::Class());
    
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void AliAnalysisTaskNetLambdaTrad::UserCreateOutputObjects()
{
    fListHist = new TList();
    fListHist->SetOwner();
    
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",2,0,2);
    fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
    fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
    fListHist->Add(fHistEventCounter);
    
    fHistCentrality = new TH1D( "fHistCentrality", "fHistCentrality",100,0,100);
    fListHist->Add(fHistCentrality);
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    const Int_t CentbinNum = 81;
    Double_t CentBins[CentbinNum+1];
    for(Int_t ic = 0; ic <= CentbinNum; ic++) CentBins[ic] = ic - 0.5;
    Double_t LambdaPtBins[32] = {1.0,1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1};
    Double_t xibinlimits[24] = {0.00,  0.20,  0.40,  0.60,  0.80,  0.90,1.00,  1.10,  1.20,  1.30,  1.40,  1.50, 1.70,  1.90,  2.20,  2.60,  3.10,  3.90,4.90,  6.00,  7.20,  8.50 ,10.00, 12.00};
    Long_t xibinnumb = sizeof(xibinlimits)/sizeof(Double_t) - 1;
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //THNSPARSE BINNING
    const Int_t dim = 63; //31 pt bins*2 + 1 cent bin
    Int_t bin[dim];
    bin[0] = 81;
    for(Int_t ibin = 1; ibin < dim; ibin++) bin[ibin] = 300;
    Double_t min[dim];
    for(Int_t jbin = 0; jbin < dim; jbin++) min[jbin] =  -0.5;
    Double_t max[dim];
    max[0] = 80.5;
    for(Int_t jbin = 1; jbin < dim; jbin++) max[jbin] = 299.5;

    
    if(fIsMC)
    {
        //--------------------------------------------------------THNSPARSE----------------------------------------------------------------------------------------------------------------

        fPtBinNplusNminusChTruth = new THnSparseI("fPtBinNplusNminusChTruth","fPtBinNplusNminusChTruth", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTruth); //gen
    
        fPtBinNplusNminusChTagFour = new THnSparseI("fPtBinNplusNminusChTagFour","fPtBinNplusNminusChTagFour", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTagFour);
        
        fPtBinNplusNminusChTagThree = new THnSparseI("fPtBinNplusNminusChTagThree","fPtBinNplusNminusChTagThree", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTagThree);
        
        fPtBinNplusNminusChTagTwo = new THnSparseI("fPtBinNplusNminusChTagTwo","fPtBinNplusNminusChTagTwo", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTagTwo);
        
        fPtBinNplusNminusChTagOne = new THnSparseI("fPtBinNplusNminusChTagOne","fPtBinNplusNminusChTagOne", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTagOne);
        
        fPtBinNplusNminusChTagFourSigThree = new THnSparseI("fPtBinNplusNminusChTagFourSigThree","fPtBinNplusNminusChTagFourSigThree", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTagFourSigThree);
        
        fPtBinNplusNminusChTagThreeSigThree = new THnSparseI("fPtBinNplusNminusChTagThreeSigThree","fPtBinNplusNminusChTagThreeSigThree", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTagThreeSigThree);
        
        fPtBinNplusNminusChTagTwoSigThree = new THnSparseI("fPtBinNplusNminusChTagTwoSigThree","fPtBinNplusNminusChTagTwoSigThree", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTagTwoSigThree);
        
        fPtBinNplusNminusChTagOneSigThree = new THnSparseI("fPtBinNplusNminusChTagOneSigThree","fPtBinNplusNminusChTagOneSigThree", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTagOneSigThree);
        //--------------------------------------------------------------STATS---------------------------------------------------------------------------------------------------------

        
        f2fHistLRecstat = new TH2F("f2fHistLRecstat","f2fHistLRecstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
        fListHist->Add(f2fHistLRecstat);
        f2fHistARecstat = new TH2F("f2fHistARecstat","f2fHistARecstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
        fListHist->Add(f2fHistARecstat);
        f2fHistLGenstat = new TH2F("f2fHistLGenstat","f2fHistLGenstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
        fListHist->Add(f2fHistLGenstat);
        f2fHistAGenstat = new TH2F("f2fHistAGenstat","f2fHistAGenstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
        fListHist->Add(f2fHistAGenstat);
        
        //-------------------------------------------------------------------GEN-----------------------------------------------------------------------------------------------

        f2fHistGenCentVsPtLambda = new TH2F( "f2fHistGenCentVsPtLambda", "Centrality Vs #Lambda Gen Pt",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistGenCentVsPtLambda);
        
        f2fHistGenCentVsPtAntiLambda = new TH2F( "f2fHistGenCentVsPtAntiLambda", "Centrality Vs #bar{#Lambda} Gen Pt",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistGenCentVsPtAntiLambda);
        
        f2fHistXiPlus = new TH2F("f2fHistXiPlus","f2fHistXiPlus",81,-0.5,80.5, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistXiPlus);
        
        f2fHistXiMinus = new TH2F("f2fHistXiMinus","f2fHistXiMinus",81,-0.5,80.5, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistXiMinus);
        
        //-------------------------------------------------------------------MC REC-----------------------------------------------------------------------------------------------
        //CONT
        f2fHistLambdaMisIdFour = new TH2F("f2fHistLambdaMisIdFour","#Lambda Mis Identified (deltaEta 1.6) ",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistLambdaMisIdFour);
        
        f2fHistAntiLambdaMisIdFour = new TH2F("f2fHistAntiLambdaMisIdFour","#bar{#Lambda}  Mis Identified (deltaEta 1.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistAntiLambdaMisIdFour);
        
        f2fHistLambdaMisIdThree = new TH2F("f2fHistLambdaMisIdThree","#Lambda Mis Identified (deltaEta 1.0) ",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistLambdaMisIdThree);
        
        f2fHistAntiLambdaMisIdThree = new TH2F("f2fHistAntiLambdaMisIdThree","#bar{#Lambda}  Mis Identified (deltaEta 1.0)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistAntiLambdaMisIdThree);
        
        f2fHistLambdaMisIdTwo = new TH2F("f2fHistLambdaMisIdTwo","#Lambda Mis Identified (deltaEta 0.6) ",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistLambdaMisIdTwo);
        
        f2fHistAntiLambdaMisIdTwo = new TH2F("f2fHistAntiLambdaMisIdTwo","#bar{#Lambda}  Mis Identified (deltaEta 0.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistAntiLambdaMisIdTwo);
        
        f2fHistLambdaMisIdOne = new TH2F("f2fHistLambdaMisIdOne","#Lambda Mis Identified (deltaEta 0.2) ",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistLambdaMisIdOne);
        
        f2fHistAntiLambdaMisIdOne = new TH2F("f2fHistAntiLambdaMisIdOne","#bar{#Lambda}  Mis Identified (deltaEta 0.2)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistAntiLambdaMisIdOne);
        
        //---
        f2fHistLambdaMisIdFourSigthree = new TH2F("f2fHistLambdaMisIdFourSigthree","#Lambda Mis Identified (deltaEta 1.6) ",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistLambdaMisIdFourSigthree);
        
        f2fHistAntiLambdaMisIdFourSigthree = new TH2F("f2fHistAntiLambdaMisIdFourSigthree","#bar{#Lambda}  Mis Identified (deltaEta 1.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistAntiLambdaMisIdFourSigthree);
        
        f2fHistLambdaMisIdThreeSigthree = new TH2F("f2fHistLambdaMisIdThreeSigthree","#Lambda Mis Identified (deltaEta 1.0) ",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistLambdaMisIdThreeSigthree);
        
        f2fHistAntiLambdaMisIdThreeSigthree = new TH2F("f2fHistAntiLambdaMisIdThreeSigthree","#bar{#Lambda}  Mis Identified (deltaEta 1.0)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistAntiLambdaMisIdThreeSigthree);
        
        f2fHistLambdaMisIdTwoSigthree = new TH2F("f2fHistLambdaMisIdTwoSigthree","#Lambda Mis Identified (deltaEta 0.6) ",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistLambdaMisIdTwoSigthree);
        
        f2fHistAntiLambdaMisIdTwoSigthree = new TH2F("f2fHistAntiLambdaMisIdTwoSigthree","#bar{#Lambda}  Mis Identified (deltaEta 0.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistAntiLambdaMisIdTwoSigthree);
        
        f2fHistLambdaMisIdOneSigthree = new TH2F("f2fHistLambdaMisIdOneSigthree","#Lambda Mis Identified (deltaEta 0.2) ",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistLambdaMisIdOneSigthree);
        
        f2fHistAntiLambdaMisIdOneSigthree = new TH2F("f2fHistAntiLambdaMisIdOneSigthree","#bar{#Lambda}  Mis Identified (deltaEta 0.2)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistAntiLambdaMisIdOneSigthree);
        
         //REC PRI
        f2fHistRecPrimariesCentVsPtLambdaFour = new TH2F("f2fHistRecPrimariesCentVsPtLambdaFour","#Lambda primaries (deltaEta 1.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaFour);
        
        f2fHistRecPrimariesCentVsPtAntiLambdaFour = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaFour","#bar{#Lambda} primaries (deltaEta 1.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaFour);
        
        f2fHistRecPrimariesCentVsPtLambdaThree = new TH2F("f2fHistRecPrimariesCentVsPtLambdaThree","#Lambda primaries (deltaEta 1.0)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaThree);
        
        f2fHistRecPrimariesCentVsPtAntiLambdaThree = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaThree","#bar{#Lambda} primaries (deltaEta 1.0)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaThree);
        
        f2fHistRecPrimariesCentVsPtLambdaTwo = new TH2F("f2fHistRecPrimariesCentVsPtLambdaTwo","#Lambda primaries (deltaEta 0.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaTwo);
        
        f2fHistRecPrimariesCentVsPtAntiLambdaTwo = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaTwo","#bar{#Lambda} primaries (deltaEta 0.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaTwo);
        
        f2fHistRecPrimariesCentVsPtLambdaOne = new TH2F("f2fHistRecPrimariesCentVsPtLambdaOne","#Lambda primaries (deltaEta 0.2)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaOne);
        
        f2fHistRecPrimariesCentVsPtAntiLambdaOne = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaOne","#bar{#Lambda} primaries (deltaEta 0.2)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaOne);
        
        //---
        f2fHistRecPrimariesCentVsPtLambdaFourSigthree = new TH2F("f2fHistRecPrimariesCentVsPtLambdaFourSigthree","#Lambda primaries (deltaEta 1.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaFourSigthree);
        
        f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree","#bar{#Lambda} primaries (deltaEta 1.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree);
        
        f2fHistRecPrimariesCentVsPtLambdaThreeSigthree = new TH2F("f2fHistRecPrimariesCentVsPtLambdaThreeSigthree","#Lambda primaries (deltaEta 1.0)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaThreeSigthree);
        
        f2fHistRecPrimariesCentVsPtAntiLambdaThreeSigthree = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaThreeSigthree","#bar{#Lambda} primaries (deltaEta 1.0)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaThreeSigthree);
        
        f2fHistRecPrimariesCentVsPtLambdaTwoSigthree = new TH2F("f2fHistRecPrimariesCentVsPtLambdaTwoSigthree","#Lambda primaries (deltaEta 0.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaTwoSigthree);
        
        f2fHistRecPrimariesCentVsPtAntiLambdaTwoSigthree = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaTwoSigthree","#bar{#Lambda} primaries (deltaEta 0.6)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaTwoSigthree);
        
        f2fHistRecPrimariesCentVsPtLambdaOneSigthree = new TH2F("f2fHistRecPrimariesCentVsPtLambdaOneSigthree","#Lambda primaries (deltaEta 0.2)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaOneSigthree);
        
        f2fHistRecPrimariesCentVsPtAntiLambdaOneSigthree = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaOneSigthree","#bar{#Lambda} primaries (deltaEta 0.2)",81,-0.5,80.5,40,1,5);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaOneSigthree);
    
        //FD
        f3fHistLambdafromXiFour = new TH3F("f3fHistLambdafromXiFour","f3fHistLambdafromXiFour (deltaEta 1.6)", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiFour);
        
        f3fHistAntiLambdafromXiFour = new TH3F("f3fHistAntiLambdafromXiFour","f3fHistAntiLambdafromXiFour (deltaEta 1.6)",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiFour);
        
        f3fHistLambdafromXiThree = new TH3F("f3fHistLambdafromXiThree","f3fHistLambdafromXiThree (deltaEta 1.0)", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiThree);
        
        f3fHistAntiLambdafromXiThree = new TH3F("f3fHistAntiLambdafromXiThree","f3fHistAntiLambdafromXiThree (deltaEta 1.0)",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiThree);
        
        f3fHistLambdafromXiTwo = new TH3F("f3fHistLambdafromXiTwo","f3fHistLambdafromXiTwo (deltaEta 0.6)", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiTwo);
        
        f3fHistAntiLambdafromXiTwo = new TH3F("f3fHistAntiLambdafromXiTwo","f3fHistAntiLambdafromXiTwo (deltaEta 0.6)",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiTwo);
        
        f3fHistLambdafromXiOne = new TH3F("f3fHistLambdafromXiOne","f3fHistLambdafromXiOne (deltaEta 0.2)", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiOne);
        
        f3fHistAntiLambdafromXiOne = new TH3F("f3fHistAntiLambdafromXiOne","f3fHistAntiLambdafromXiOne (deltaEta 0.2)",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiOne);
        
        //---
        f3fHistLambdafromXiFourSigthree = new TH3F("f3fHistLambdafromXiFourSigthree","f3fHistLambdafromXiFourSigthree (deltaEta 1.6)", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiFourSigthree);
        
        f3fHistAntiLambdafromXiFourSigthree = new TH3F("f3fHistAntiLambdafromXiFourSigthree","f3fHistAntiLambdafromXiFourSigthree (deltaEta 1.6)",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiFourSigthree);
        
        f3fHistLambdafromXiThreeSigthree = new TH3F("f3fHistLambdafromXiThreeSigthree","f3fHistLambdafromXiThreeSigthree (deltaEta 1.0)", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiThreeSigthree);
        
        f3fHistAntiLambdafromXiThreeSigthree = new TH3F("f3fHistAntiLambdafromXiThreeSigthree","f3fHistAntiLambdafromXiThreeSigthree (deltaEta 1.0)",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiThreeSigthree);
        
        f3fHistLambdafromXiTwoSigthree = new TH3F("f3fHistLambdafromXiTwoSigthree","f3fHistLambdafromXiTwoSigthree (deltaEta 0.6)", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiTwoSigthree);
        
        f3fHistAntiLambdafromXiTwoSigthree = new TH3F("f3fHistAntiLambdafromXiTwoSigthree","f3fHistAntiLambdafromXiTwoSigthree (deltaEta 0.6)",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiTwoSigthree);
        
        f3fHistLambdafromXiOneSigthree = new TH3F("f3fHistLambdafromXiOneSigthree","f3fHistLambdafromXiOneSigthree (deltaEta 0.2)", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiOneSigthree);
        
        f3fHistAntiLambdafromXiOneSigthree = new TH3F("f3fHistAntiLambdafromXiOneSigthree","f3fHistAntiLambdafromXiOneSigthree (deltaEta 0.2)",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiOneSigthree);
        
        //REC
        f3fHistCentInvMassVsPtLambdaRecFour = new TH3F("f3fHistCentInvMassVsPtLambdaRecFour","f3fHistCentInvMassVsPtLambdaRecFour (deltaEta 1.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecFour);
        
        f3fHistCentInvMassVsPtAntiLambdaRecFour = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecFour","f3fHistCentInvMassVsPtAntiLambdaRecFour (deltaEta 1.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecFour);
        
        f3fHistCentInvMassVsPtLambdaRecThree = new TH3F("f3fHistCentInvMassVsPtLambdaRecThree","f3fHistCentInvMassVsPtLambdaRecThree (deltaEta 1.0)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecThree);
        
        f3fHistCentInvMassVsPtAntiLambdaRecThree = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecThree","f3fHistCentInvMassVsPtAntiLambdaRecThree (deltaEta 1.0)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecThree);
        
        f3fHistCentInvMassVsPtLambdaRecTwo = new TH3F("f3fHistCentInvMassVsPtLambdaRecTwo","f3fHistCentInvMassVsPtLambdaRecTwo (deltaEta 0.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecTwo);
        
        f3fHistCentInvMassVsPtAntiLambdaRecTwo = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecTwo","f3fHistCentInvMassVsPtAntiLambdaRecTwo (deltaEta 0.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecTwo);
        
        f3fHistCentInvMassVsPtLambdaRecOne = new TH3F("f3fHistCentInvMassVsPtLambdaRecOne","f3fHistCentInvMassVsPtLambdaRecOne (deltaEta 0.2)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecOne);
        
        f3fHistCentInvMassVsPtAntiLambdaRecOne = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecOne","f3fHistCentInvMassVsPtAntiLambdaRecOne (deltaEta 0.2)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecOne);
        
        //---
        f3fHistCentInvMassVsPtLambdaRecFourSigthree = new TH3F("f3fHistCentInvMassVsPtLambdaRecFourSigthree","f3fHistCentInvMassVsPtLambdaRecFourSigthree (deltaEta 1.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecFourSigthree);
        
        f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree","f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree (deltaEta 1.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree);
        
        f3fHistCentInvMassVsPtLambdaRecThreeSigthree = new TH3F("f3fHistCentInvMassVsPtLambdaRecThreeSigthree","f3fHistCentInvMassVsPtLambdaRecThreeSigthree (deltaEta 1.0)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecThreeSigthree);
        
        f3fHistCentInvMassVsPtAntiLambdaRecThreeSigthree = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecThreeSigthree","f3fHistCentInvMassVsPtAntiLambdaRecThreeSigthree (deltaEta 1.0)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecThreeSigthree);
        
        f3fHistCentInvMassVsPtLambdaRecTwoSigthree = new TH3F("f3fHistCentInvMassVsPtLambdaRecTwoSigthree","f3fHistCentInvMassVsPtLambdaRecTwoSigthree (deltaEta 0.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecTwoSigthree);
        
        f3fHistCentInvMassVsPtAntiLambdaRecTwoSigthree = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecTwoSigthree","f3fHistCentInvMassVsPtAntiLambdaRecTwoSigthree (deltaEta 0.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecTwoSigthree);
        
        f3fHistCentInvMassVsPtLambdaRecOneSigthree = new TH3F("f3fHistCentInvMassVsPtLambdaRecOneSigthree","f3fHistCentInvMassVsPtLambdaRecOneSigthree (deltaEta 0.2)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecOneSigthree);
        
        f3fHistCentInvMassVsPtAntiLambdaRecOneSigthree = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecOneSigthree","f3fHistCentInvMassVsPtAntiLambdaRecOneSigthree (deltaEta 0.2)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecOneSigthree);
    }
    
    PostData(1,fListHist);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNetLambdaTrad::UserExec(Option_t *)
{
    
    const Int_t dim = fNptBins*2;
    
    Int_t ptChMC[dim];
    Int_t ptChTagOne[dim];
    Int_t ptChTagTwo[dim];
    Int_t ptChTagThree[dim];
    Int_t ptChTagFour[dim];
    
    Int_t ptChTagOneSigthree[dim];
    Int_t ptChTagTwoSigthree[dim];
    Int_t ptChTagThreeSigthree[dim];
    Int_t ptChTagFourSigthree[dim];
    
    for(Int_t idx = 0; idx < dim; idx++)
    {
        ptChMC[idx] = 0;
        ptChTagOne[idx] = 0.;
        ptChTagTwo[idx] = 0.;
        ptChTagThree[idx] = 0.;
        ptChTagFour[idx] = 0.;
        
        ptChTagOneSigthree[idx] = 0.;
        ptChTagTwoSigthree[idx] = 0.;
        ptChTagThreeSigthree[idx] = 0.;
        ptChTagFourSigthree[idx] = 0.;
    }
    
    
    if (!fInputEvent) return;
    
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) return;
 
    fPIDResponse = fInputHandler->GetPIDResponse();
    if(!fPIDResponse) return;
    
    AliStack *stack = 0x0;
    if(fIsMC)
    {
            AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
            if(!mcH) return;
            fMCEvent=mcH->MCEvent();
        
            if(!fMCEvent) return;
        
            stack = fMCEvent->Stack();
            if(!stack) return;
    }
    Double_t lMagneticField = -10;
    lMagneticField = fESD->GetMagneticField();
    
    AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent->FindListObject("MultSelection");
    if(!MultSelection) return;
    
    if(!(fInputHandler->IsEventSelected() & fEvSel)) return;
    
    Double_t vVtx[3];
    
    AliVVertex *vvertex = (AliVVertex*)fInputEvent->GetPrimaryVertex();
    if (!vvertex) return;
    vVtx[0] = vvertex->GetX();
    vVtx[1] = vvertex->GetY();
    vVtx[2] = vvertex->GetZ();
    
    if(vVtx[2] < -10. || vVtx[2] > 10.) return;
    
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    
    if( fCentrality < 0 || fCentrality >=80 ) return;
    if (!fEventCuts.AcceptEvent(fInputEvent)) return;//pileup cut
    
    fHistEventCounter->Fill(1.5);
    fHistCentrality->Fill(fCentrality);
    
    const AliESDVertex *lPrimaryBestESDVtx     = fESD->GetPrimaryVertex();
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    
    
    Int_t nGen = 0;
    Double_t lRap = 0.0;
    Double_t nRecL = 0.0;
    Double_t nRecA = 0.0;
    Double_t nGenL = 0.0;
    Double_t nGenA = 0.0;
    
    if(fIsMC)
    {
        
        nGen = stack->GetNtrack();
        for(Int_t iGen = 0; iGen < nGen; iGen++)
        {
            Int_t genpid = -1;
            Float_t gpt = 0.0, eta = 0.0,abseta =0.0;
   
                TParticle* mctrack = stack->Particle(iGen);
                if(!mctrack) continue;
                if(!(stack->IsPhysicalPrimary(iGen))) continue;
                genpid = mctrack->GetPdgCode();
                gpt = mctrack->Pt();
                eta = mctrack->Eta();
    
            abseta = TMath::Abs(eta);
            if(abseta > 0.8) continue;
            
            Int_t iptbinMC = GetPtBin(gpt);
            
            if(genpid == 3122)
            {
                f2fHistGenCentVsPtLambda->Fill(fCentrality, gpt);
                nGenL += 1.;
                ptChMC[iptbinMC] += 1;
            }
            
            if(genpid == -3122)
            {
                f2fHistGenCentVsPtAntiLambda->Fill(fCentrality,gpt);
                nGenA += 1.;
                ptChMC[iptbinMC+fNptBins] += 1;
            }
            
            if(genpid == -3312)
            {
                f2fHistXiPlus->Fill(fCentrality,gpt);
            }
            
            if(genpid == 3312)
            {
                f2fHistXiMinus->Fill(fCentrality,gpt);
            }
            
            
        } // end loop over generated particles
        //------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        f2fHistLGenstat->Fill(fCentrality, nGenL);
        f2fHistAGenstat->Fill(fCentrality, nGenA);
        
        //------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        Double_t ptContainerMC[dim+1];
        ptContainerMC[0] = (Double_t) fCentrality;
        for(Int_t i = 1; i <= dim; i++)
        {
            ptContainerMC[i] = ptChMC[i-1];
        }
        fPtBinNplusNminusChTruth->Fill(ptContainerMC);
        //------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    }//MC condition on gen
    
    Int_t nV0 = 0;
    nV0 = fESD->GetNumberOfV0s();
    AliESDv0 *esdv0 = 0x0;
    AliESDtrack *esdpTrack = 0x0;
    AliESDtrack *esdnTrack = 0x0;
 
    Double_t fMinV0Pt = 0.5;
    Double_t fMaxV0Pt = 4.5;
    
    for(Int_t iV0 = 0; iV0 < nV0; iV0++)
    {
        esdv0 = 0x0;
        esdpTrack = 0x0;
        esdnTrack = 0x0;
        
        Double_t  vertx[3];
        
        Float_t invMassLambda = -999, invMassAntiLambda = -999;
        Float_t V0pt = -999, eta = -999, pmom = -999;
        Float_t ppt = -999,  peta = -999, posprnsg = -999, pospion =-999, v0Radius =-999, lRapLambda=-999;
        Float_t npt = -999,  neta = -999, negprnsg = -999, negpion =-999;
        Bool_t  ontheflystat = kFALSE;
        Float_t dcaPosToVertex = -999, dcaNegToVertex = -999, dcaDaughters = -999, dcaV0ToVertex = -999, cosPointingAngle = -999;
        
        esdv0 = fESD->GetV0(iV0);
        if(!esdv0) continue;
        esdpTrack =  (AliESDtrack*)fESD->GetTrack(TMath::Abs(esdv0->GetPindex()));
        if(!esdpTrack) continue;
        esdnTrack = (AliESDtrack*)fESD->GetTrack(TMath::Abs(esdv0->GetNindex()));
        if(!esdnTrack) continue;
        
        if(esdpTrack->Charge() == esdnTrack->Charge())
        {
            continue;
        }
        esdv0->GetXYZ(vertx[0], vertx[1], vertx[2]); //decay vertex
        v0Radius = TMath::Sqrt(vertx[0]*vertx[0]+vertx[1]*vertx[1]);
        
        lRapLambda  = esdv0->RapLambda();
        V0pt = esdv0->Pt();
        if ((V0pt<fMinV0Pt)||(fMaxV0Pt<V0pt)) continue;
 //--------------------------------------------------------------------Track selection-------------------------------------------------------------------------
            Float_t lPosTrackCrossedRows = esdpTrack->GetTPCClusterInfo(2,1);
            Float_t lNegTrackCrossedRows = esdnTrack->GetTPCClusterInfo(2,1);
            
            fTreeVariableLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
            if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
                fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
            
            if( !(esdpTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
            if( !(esdnTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
            
            if ( ( ( esdpTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( esdnTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;
            
            if( esdpTrack->GetKinkIndex(0)>0 || esdnTrack->GetKinkIndex(0)>0 ) continue;
            
            if( esdpTrack->GetTPCNclsF()<=0 || esdnTrack->GetTPCNclsF()<=0 ) continue;
            
            Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(esdpTrack->GetTPCNclsF()));
            Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(esdnTrack->GetTPCNclsF()));
            
            fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
            if( lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable )
                fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
            
            if ( fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) continue;
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            pmom = esdv0->P();
            eta = esdv0->Eta();
            ppt = esdpTrack->Pt();
            peta = esdpTrack->Eta();
            npt = esdnTrack->Pt();
            neta = esdnTrack->Eta();
            ontheflystat = esdv0->GetOnFlyStatus();
            
            dcaPosToVertex = TMath::Abs(esdpTrack->GetD(lBestPrimaryVtxPos[0],
                                                        lBestPrimaryVtxPos[1],
                                                        lMagneticField) );
            
            dcaNegToVertex = TMath::Abs(esdnTrack->GetD(lBestPrimaryVtxPos[0],
                                                        lBestPrimaryVtxPos[1],
                                                        lMagneticField) );
            
            cosPointingAngle = esdv0->GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
            dcaDaughters = esdv0->GetDcaV0Daughters();
            dcaV0ToVertex = esdv0->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
            
            
            esdv0->ChangeMassHypothesis(3122);
            invMassLambda = esdv0->GetEffMass();
            esdv0->ChangeMassHypothesis(-3122);
            invMassAntiLambda = esdv0->GetEffMass();
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        posprnsg = fPIDResponse->NumberOfSigmasTPC(esdpTrack, AliPID::kProton);
        negprnsg = fPIDResponse->NumberOfSigmasTPC(esdnTrack, AliPID::kProton);
        pospion  = fPIDResponse->NumberOfSigmasTPC( esdpTrack, AliPID::kPion );
        negpion  = fPIDResponse->NumberOfSigmasTPC( esdnTrack, AliPID::kPion );
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if(TMath::Abs(peta) > 1) continue;
        if(TMath::Abs(neta) > 1) continue;
        if(cosPointingAngle < 0.99) continue;
        if(dcaDaughters > 0.8) continue;
        if(v0Radius < 5.0) continue;
        if(v0Radius > 200.) continue;
    
        if( ontheflystat == 0 )
        {
            if(fIsMC)
            {
                    fTreeVariablePID = -999, fTreeVariablePIDParent = -999;
                    Float_t mcpt = -999, mceta = -999, fTreeVariablePtParent = -999;
                    Bool_t isPrim = kFALSE, isSecFromMaterial = kFALSE, isSecFromWeakDecay = kFALSE, isPrimParent =kFALSE;
                    fTreeVariablePIDPositive = -999;
                    fTreeVariablePIDNegative = -999;
        
                    if(TMath::Abs(esdpTrack->GetLabel()) >= nGen || TMath::Abs(esdnTrack->GetLabel()) >= nGen) continue;
                    Int_t lblPosV0Dghter = (Int_t) TMath::Abs( esdpTrack->GetLabel() );
                    Int_t lblNegV0Dghter = (Int_t) TMath::Abs( esdnTrack->GetLabel() );
                    
                    TParticle* esdGenTrackPos = stack->Particle( lblPosV0Dghter );
                    TParticle* esdGenTrackNeg = stack->Particle( lblNegV0Dghter );
                    
                    Int_t lPIDPositive = esdGenTrackPos -> GetPdgCode();
                    Int_t lPIDNegative = esdGenTrackNeg -> GetPdgCode();
                    
                    fTreeVariablePIDPositive = lPIDPositive;
                    fTreeVariablePIDNegative = lPIDNegative;
                    
                    Int_t posTparticle = esdGenTrackPos->GetFirstMother();
                    Int_t negTparticle = esdGenTrackNeg->GetFirstMother();
                    
                    if( posTparticle == negTparticle && posTparticle > 0 )
                    {
                        TParticle *esdlthisV0 = stack->Particle(posTparticle);
                        if(!esdlthisV0) continue;
                        fTreeVariablePID = esdlthisV0->GetPdgCode();
                        
                        mcpt = esdlthisV0->Pt();
                        mceta = esdlthisV0->Eta();
                        isSecFromMaterial = stack->IsSecondaryFromMaterial(posTparticle);
                        isSecFromWeakDecay = stack->IsSecondaryFromWeakDecay(posTparticle);
                        isPrim = stack->IsPhysicalPrimary(posTparticle);
                        
                        Int_t esdlthisV0parent = esdlthisV0->GetFirstMother();
                        if (esdlthisV0parent > 0)
                        {
                            TParticle *lbV0parent = stack->Particle(esdlthisV0parent);
                            if(!lbV0parent) continue;
                            fTreeVariablePIDParent = lbV0parent->GetPdgCode();
                            fTreeVariablePtParent = lbV0parent->Pt();
                            isPrimParent =  stack->IsPhysicalPrimary(esdlthisV0parent);
                        }
                    }
                    
                    Int_t iptbinRecPRI = GetPtBin(mcpt);
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex >  0.1  && TMath::Abs(posprnsg)  <= 2. && TMath::Abs(negpion)  <= 2.)
                    {
                    if(TMath::Abs(eta) < 0.8) //MC-TEST
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtLambdaFour->Fill(fCentrality,mcpt);}
                            else{f2fHistLambdaMisIdFour->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == 3312)
                                    {
                                        f3fHistLambdafromXiFour->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            f3fHistCentInvMassVsPtLambdaRecFour->Fill(fCentrality,invMassLambda,mcpt);
                            nRecL += 1.;
                            ptChTagFour[iptbinRecPRI] += 1;
                        }
                    }
                        
                        if(TMath::Abs(eta) < 0.5)
                        {
                            if(fTreeVariablePID == 3122)
                            {
                                if(isPrim){f2fHistRecPrimariesCentVsPtLambdaThree->Fill(fCentrality,mcpt);}
                                else{f2fHistLambdaMisIdThree->Fill(fCentrality,mcpt);}
                                if(isSecFromWeakDecay)
                                {
                                    if(isPrimParent)
                                    {
                                        if(fTreeVariablePIDParent == 3312)
                                        {
                                            f3fHistLambdafromXiThree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                        }
                                    }
                                }
                                f3fHistCentInvMassVsPtLambdaRecThree->Fill(fCentrality,invMassLambda,mcpt);
                                ptChTagThree[iptbinRecPRI] += 1;
                            }
                        }
                        
                        if(TMath::Abs(eta) < 0.3)
                        {
                            if(fTreeVariablePID == 3122)
                            {
                                if(isPrim){f2fHistRecPrimariesCentVsPtLambdaTwo->Fill(fCentrality,mcpt);}
                                else{f2fHistLambdaMisIdTwo->Fill(fCentrality,mcpt);}
                                if(isSecFromWeakDecay)
                                {
                                    if(isPrimParent)
                                    {
                                        if(fTreeVariablePIDParent == 3312)
                                        {
                                            f3fHistLambdafromXiTwo->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                        }
                                    }
                                }
                                f3fHistCentInvMassVsPtLambdaRecTwo->Fill(fCentrality,invMassLambda,mcpt);
                                ptChTagTwo[iptbinRecPRI] += 1;
                            }
                        }
                        
                        if(TMath::Abs(eta) < 0.1)
                        {
                            if(fTreeVariablePID == 3122)
                            {
                                if(isPrim){f2fHistRecPrimariesCentVsPtLambdaOne->Fill(fCentrality,mcpt);}
                                else{f2fHistLambdaMisIdOne->Fill(fCentrality,mcpt);}
                                if(isSecFromWeakDecay)
                                {
                                    if(isPrimParent)
                                    {
                                        if(fTreeVariablePIDParent == 3312)
                                        {
                                            f3fHistLambdafromXiOne->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                        }
                                    }
                                }
                                f3fHistCentInvMassVsPtLambdaRecOne->Fill(fCentrality,invMassLambda,mcpt);
                                ptChTagOne[iptbinRecPRI] += 1;
                            }
                        }
                    }
                Int_t iptbinRec = GetPtBin(mcpt);
                if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex >  0.1  && TMath::Abs(posprnsg)  <= 3. && TMath::Abs(negpion)  <= 3.)
                {
                    if(TMath::Abs(eta) < 0.8) //MC-TEST
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtLambdaFourSigthree->Fill(fCentrality,mcpt);}
                            else{f2fHistLambdaMisIdFourSigthree->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == 3312)
                                    {
                                        f3fHistLambdafromXiFourSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            f3fHistCentInvMassVsPtLambdaRecFourSigthree->Fill(fCentrality,invMassLambda,mcpt);
                            ptChTagFourSigthree[iptbinRec] += 1;
                        }
                    }
                    
                    if(TMath::Abs(eta) < 0.5)
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtLambdaThreeSigthree->Fill(fCentrality,mcpt);}
                            else{f2fHistLambdaMisIdThreeSigthree->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == 3312)
                                    {
                                        f3fHistLambdafromXiThreeSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            f3fHistCentInvMassVsPtLambdaRecThreeSigthree->Fill(fCentrality,invMassLambda,mcpt);
                            ptChTagThreeSigthree[iptbinRec] += 1;
                        }
                    }
                    
                    if(TMath::Abs(eta) < 0.3)
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtLambdaTwoSigthree->Fill(fCentrality,mcpt);}
                            else{f2fHistLambdaMisIdTwoSigthree->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == 3312)
                                    {
                                        f3fHistLambdafromXiTwoSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            f3fHistCentInvMassVsPtLambdaRecTwoSigthree->Fill(fCentrality,invMassLambda,mcpt);
                            ptChTagTwoSigthree[iptbinRec] += 1;
                        }
                    }
                    
                    if(TMath::Abs(eta) < 0.1)
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtLambdaOneSigthree->Fill(fCentrality,mcpt);}
                            else{f2fHistLambdaMisIdOneSigthree->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == 3312)
                                    {
                                        f3fHistLambdafromXiOneSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            f3fHistCentInvMassVsPtLambdaRecOneSigthree->Fill(fCentrality,invMassLambda,mcpt);
                            ptChTagOneSigthree[iptbinRec] += 1;
                        }
                    }
                }
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.25 && TMath::Abs(negprnsg)  <= 2. && TMath::Abs(pospion)  <= 2.)
                    {
                    if(TMath::Abs(eta) < 0.8)//MC-TEST
                    {
                       if(fTreeVariablePID == -3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambdaFour->Fill(fCentrality,mcpt);}
                            else{f2fHistAntiLambdaMisIdFour->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == -3312)
                                    {
                                        f3fHistAntiLambdafromXiFour->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            f3fHistCentInvMassVsPtAntiLambdaRecFour->Fill(fCentrality,invMassAntiLambda,mcpt);
                            nRecA += 1.;
                            ptChTagFour[iptbinRecPRI+fNptBins] += 1;
                        }
                    }
                        if(TMath::Abs(eta) < 0.5)
                        {
                            if(fTreeVariablePID == -3122)
                            {
                                if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambdaThree->Fill(fCentrality,mcpt);}
                                else{f2fHistAntiLambdaMisIdThree->Fill(fCentrality,mcpt);}
                                if(isSecFromWeakDecay)
                                {
                                    if(isPrimParent)
                                    {
                                        if(fTreeVariablePIDParent == -3312)
                                        {
                                            f3fHistAntiLambdafromXiThree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                        }
                                    }
                                }
                                f3fHistCentInvMassVsPtAntiLambdaRecThree->Fill(fCentrality,invMassAntiLambda,mcpt);
                                ptChTagThree[iptbinRecPRI+fNptBins] += 1;
                            }
                        }
                        if(TMath::Abs(eta) < 0.3)
                        {
                            if(fTreeVariablePID == -3122)
                            {
                                if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambdaTwo->Fill(fCentrality,mcpt);}
                                else{f2fHistAntiLambdaMisIdTwo->Fill(fCentrality,mcpt);}
                                if(isSecFromWeakDecay)
                                {
                                    if(isPrimParent)
                                    {
                                        if(fTreeVariablePIDParent == -3312)
                                        {
                                            f3fHistAntiLambdafromXiTwo->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                        }
                                    }
                                }
                                f3fHistCentInvMassVsPtAntiLambdaRecTwo->Fill(fCentrality,invMassAntiLambda,mcpt);
                                ptChTagTwo[iptbinRecPRI+fNptBins] += 1;
                            }
                        }
                        if(TMath::Abs(eta) < 0.1)
                        {
                            if(fTreeVariablePID == -3122)
                            {
                                if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambdaOne->Fill(fCentrality,mcpt);}
                                else{f2fHistAntiLambdaMisIdOne->Fill(fCentrality,mcpt);}
                                if(isSecFromWeakDecay)
                                {
                                    if(isPrimParent)
                                    {
                                        if(fTreeVariablePIDParent == -3312)
                                        {
                                            f3fHistAntiLambdafromXiOne->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                        }
                                    }
                                }
                                f3fHistCentInvMassVsPtAntiLambdaRecOne->Fill(fCentrality,invMassAntiLambda,mcpt);
                                ptChTagOne[iptbinRecPRI+fNptBins] += 1;
                            }
                        }
                    }
                if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.25 && TMath::Abs(negprnsg)  <= 3. && TMath::Abs(pospion)  <= 3.)
                {
                    if(TMath::Abs(eta) < 0.8)//MC-TEST
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree->Fill(fCentrality,mcpt);}
                            else{f2fHistAntiLambdaMisIdFourSigthree->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == -3312)
                                    {
                                        f3fHistAntiLambdafromXiFourSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree->Fill(fCentrality,invMassAntiLambda,mcpt);
                            ptChTagFourSigthree[iptbinRec+fNptBins] += 1;
                        }
                    }
                    if(TMath::Abs(eta) < 0.5)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambdaThreeSigthree->Fill(fCentrality,mcpt);}
                            else{f2fHistAntiLambdaMisIdThreeSigthree->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == -3312)
                                    {
                                        f3fHistAntiLambdafromXiThreeSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            f3fHistCentInvMassVsPtAntiLambdaRecThreeSigthree->Fill(fCentrality,invMassAntiLambda,mcpt);
                            ptChTagThreeSigthree[iptbinRec+fNptBins] += 1;
                        }
                    }
                    if(TMath::Abs(eta) < 0.3)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambdaTwoSigthree->Fill(fCentrality,mcpt);}
                            else{f2fHistAntiLambdaMisIdTwoSigthree->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == -3312)
                                    {
                                        f3fHistAntiLambdafromXiTwoSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            f3fHistCentInvMassVsPtAntiLambdaRecTwoSigthree->Fill(fCentrality,invMassAntiLambda,mcpt);
                            ptChTagTwoSigthree[iptbinRec+fNptBins] += 1;
                        }
                    }
                    if(TMath::Abs(eta) < 0.1)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambdaOneSigthree->Fill(fCentrality,mcpt);}
                            else{f2fHistAntiLambdaMisIdOneSigthree->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == -3312)
                                    {
                                        f3fHistAntiLambdafromXiOneSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            f3fHistCentInvMassVsPtAntiLambdaRecOneSigthree->Fill(fCentrality,invMassAntiLambda,mcpt);
                            ptChTagOneSigthree[iptbinRec+fNptBins] += 1;
                        }
                    }
                }
            } //MC condition
        }// zero onfly V0
    }// end of V0 loop
    
    
    f2fHistLRecstat->Fill(fCentrality, nRecL);
    f2fHistARecstat->Fill(fCentrality, nRecA);
    //-------------------------------------------------
    Double_t ptContainerFour[dim+1];
    ptContainerFour[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerFour[i] = ptChTagFour[i-1];
    }
    fPtBinNplusNminusChTagFour->Fill(ptContainerFour);
    //-------------------------------------------------
    Double_t ptContainerThree[dim+1];
    ptContainerThree[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerThree[i] = ptChTagThree[i-1];
    }
    fPtBinNplusNminusChTagThree->Fill(ptContainerThree);
    //-------------------------------------------------
    Double_t ptContainerTwo[dim+1];
    ptContainerTwo[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerTwo[i] = ptChTagTwo[i-1];
    }
    fPtBinNplusNminusChTagTwo->Fill(ptContainerTwo);
    //-------------------------------------------------
    Double_t ptContainerOne[dim+1];
    ptContainerOne[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerOne[i] = ptChTagOne[i-1];
    }
    fPtBinNplusNminusChTagOne->Fill(ptContainerOne);
    //-------------------------------------------------
    
    //SIGTHREE
    //-------------------------------------------------
    Double_t ptContainerFourSigThree[dim+1];
    ptContainerFourSigThree[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerFourSigThree[i] = ptChTagFourSigthree[i-1];
    }
    fPtBinNplusNminusChTagFourSigThree->Fill(ptContainerFourSigThree);
    //-------------------------------------------------
    Double_t ptContainerThreeSigThree[dim+1];
    ptContainerThreeSigThree[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerThreeSigThree[i] = ptChTagThreeSigthree[i-1];
    }
    fPtBinNplusNminusChTagThreeSigThree->Fill(ptContainerThreeSigThree);
    //-------------------------------------------------
    Double_t ptContainerTwoSigThree[dim+1];
    ptContainerTwoSigThree[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerTwoSigThree[i] = ptChTagTwoSigthree[i-1];
    }
    fPtBinNplusNminusChTagTwoSigThree->Fill(ptContainerTwoSigThree);
    //-------------------------------------------------
    Double_t ptContainerOneSigThree[dim+1];
    ptContainerOneSigThree[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerOneSigThree[i] = ptChTagOneSigthree[i-1];
    }
    fPtBinNplusNminusChTagOneSigThree->Fill(ptContainerOneSigThree);
    //-------------------------------------------------
   
    PostData(1,fListHist);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskNetLambdaTrad::GetPtBin(Double_t pt)
{
    Int_t bin = -1;
    
    Double_t LambdaPtBins[32] = {1.0,1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1};

    for(Int_t iBin = 0; iBin < fNptBins; iBin++)
    {
        
        if( iBin == fNptBins-1){
            if( pt >= LambdaPtBins[iBin] && pt <= LambdaPtBins[iBin+1]){
                bin = iBin;
                break;
            }
        }
        else{
            if( pt >= LambdaPtBins[iBin] && pt < LambdaPtBins[iBin+1]){
                bin = iBin;
                break;
                
            }
        }
    }
    
    return bin;
    
}



