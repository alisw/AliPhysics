
// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Naomi Umaka Apr 2018
// Updated May 6


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



f3fHistCentVsInvMassLambda1point0(0x0),
f3fHistCentVsInvMassLambda1point0Masscut(0x0),
f3fHistCentVsInvMassLambda1point0Sigtwo(0x0),
f3fHistCentVsInvMassLambda1point0SigtwoMasscut(0x0),
f3fHistCentVsInvMassLambda1point0Sigfour(0x0),
f3fHistCentVsInvMassLambda1point0SigfourMasscut(0x0),
f3fHistCentVsInvMassLambda1point0postight(0x0),

f3fHistCentVsInvMassLambda1point0postightMasscut(0x0),
f3fHistCentVsInvMassLambda1point0posloose(0x0),
f3fHistCentVsInvMassLambda1point0poslooseMasscut(0x0),
f3fHistCentVsInvMassLambda1point0negtight(0x0),
f3fHistCentVsInvMassLambda1point0negtightMasscut(0x0),
f3fHistCentVsInvMassLambda1point0negloose(0x0),
f3fHistCentVsInvMassLambda1point0neglooseMasscut(0x0),

f3fHistCentVsInvMassAntiLambda1point0(0x0),
f3fHistCentVsInvMassAntiLambda1point0Masscut(0x0),
f3fHistCentVsInvMassAntiLambda1point0Sigtwo(0x0),
f3fHistCentVsInvMassAntiLambda1point0SigtwoMasscut(0x0),
f3fHistCentVsInvMassAntiLambda1point0Sigfour(0x0),
f3fHistCentVsInvMassAntiLambda1point0SigfourMasscut(0x0),
f3fHistCentVsInvMassAntiLambda1point0postight(0x0),

f3fHistCentVsInvMassAntiLambda1point0postightMasscut(0x0),
f3fHistCentVsInvMassAntiLambda1point0posloose(0x0),
f3fHistCentVsInvMassAntiLambda1point0poslooseMasscut(0x0),
f3fHistCentVsInvMassAntiLambda1point0negtight(0x0),
f3fHistCentVsInvMassAntiLambda1point0negtightMasscut(0x0),
f3fHistCentVsInvMassAntiLambda1point0negloose(0x0),
f3fHistCentVsInvMassAntiLambda1point0neglooseMasscut(0x0),

fCentrality(-1),
fNptBins(23),

fEvSel(AliVEvent::kINT7),

fPtBinNplusNminusChnegtightBKG(NULL),
fPtBinNplusNminusChnegtight(NULL),
fPtBinNplusNminusChneglooseBKG(NULL),
fPtBinNplusNminusChnegloose(NULL),
fPtBinNplusNminusChposlooseBKG(NULL),
fPtBinNplusNminusChposloose(NULL),
fPtBinNplusNminusChpostightBKG(NULL),
fPtBinNplusNminusChpostight(NULL),
fPtBinNplusNminusChSigfourBKG(NULL),
fPtBinNplusNminusChSigfour(NULL),
fPtBinNplusNminusChSigtwoBKG(NULL),
fPtBinNplusNminusChSigtwo(NULL),
fPtBinNplusNminusChBKG(NULL),
fPtBinNplusNminusCh(NULL)


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
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    const Int_t CentbinNum = 81;
    Double_t CentBins[CentbinNum+1];
    for(Int_t ic = 0; ic <= CentbinNum; ic++) CentBins[ic] = ic - 0.5;
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //pt binning
    Double_t LambdaPtBins[24] = {0.8,0.9,1.0,1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,4.2};
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Double_t MassBins[103]
    = {1.0788,1.0796,1.0804,1.0812,1.082,1.0828,1.0836,1.0844,1.0852,1.086,1.0868,1.0876,1.0884,1.0892,1.09,1.0908,1.0916,1.0924,1.0932,1.094,1.0948,1.0956,1.0964,1.0972,1.098,1.0988,1.0996,1.1004,1.1012,1.102,
        1.1028,1.1036,1.1044,1.1052,1.106,1.1068,1.1076,1.1084,1.1092,1.11,1.1108,1.1116,1.1124,1.1132,1.114,1.1148,1.1156,1.1164,1.1172,1.118,1.1188,1.1196,1.1204,1.1212,1.122,1.1228,1.1236,1.1244,
        1.1252,1.126,1.1268,1.1276,1.1284,1.1292,1.13,1.1308,1.1316,1.1324,1.1332,1.134,1.1348,1.1356,1.1364,1.1372,1.138,1.1388,1.1396,1.1404,1.1412,1.142,1.1428,1.1436,1.1444,1.1452,1.146,1.1468,1.1476,
        1.1484,1.1492,1.15,1.1508,1.1516,1.1524,1.1532,1.154,1.1548,1.1556,1.1564,1.1572,1.158,1.1588,1.1596,1.1604};
    Long_t Massbinnumb = sizeof(MassBins)/sizeof(Double_t) - 1;
    
    //V0 hists//
    
    
    f3fHistCentVsInvMassLambda1point0 = new TH3F("f3fHistCentVsInvMassLambda1point0","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0);
    
    f3fHistCentVsInvMassLambda1point0Masscut = new TH3F("f3fHistCentVsInvMassLambda1point0Masscut","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0Masscut);
    
    f3fHistCentVsInvMassLambda1point0Sigtwo = new TH3F("f3fHistCentVsInvMassLambda1point0Sigtwo","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0Sigtwo);
    
    f3fHistCentVsInvMassLambda1point0SigtwoMasscut = new TH3F("f3fHistCentVsInvMassLambda1point0SigtwoMasscut","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0SigtwoMasscut);
    
    
    f3fHistCentVsInvMassLambda1point0Sigfour = new TH3F("f3fHistCentVsInvMassLambda1point0Sigfour","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0Sigfour);
    
    f3fHistCentVsInvMassLambda1point0SigfourMasscut = new TH3F("f3fHistCentVsInvMassLambda1point0SigfourMasscut","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0SigfourMasscut);
    
    
    f3fHistCentVsInvMassLambda1point0postight = new TH3F("f3fHistCentVsInvMassLambda1point0postight","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0postight);
    
    f3fHistCentVsInvMassLambda1point0postightMasscut = new TH3F("f3fHistCentVsInvMassLambda1point0postightMasscut","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0postightMasscut);
    
    
    f3fHistCentVsInvMassLambda1point0posloose = new TH3F("f3fHistCentVsInvMassLambda1point0posloose","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0posloose);
    
    
    f3fHistCentVsInvMassLambda1point0poslooseMasscut = new TH3F("f3fHistCentVsInvMassLambda1point0poslooseMasscut","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0poslooseMasscut);
    
    f3fHistCentVsInvMassLambda1point0negtight = new TH3F("f3fHistCentVsInvMassLambda1point0negtight","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0negtight);
    
    f3fHistCentVsInvMassLambda1point0negtightMasscut = new TH3F("f3fHistCentVsInvMassLambda1point0negtightMasscut","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0negtightMasscut);
    
    f3fHistCentVsInvMassLambda1point0negloose = new TH3F("f3fHistCentVsInvMassLambda1point0negloose","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0negloose);
    
    f3fHistCentVsInvMassLambda1point0neglooseMasscut = new TH3F("f3fHistCentVsInvMassLambda1point0neglooseMasscut","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0neglooseMasscut);
    
    //////////////////
    f3fHistCentVsInvMassAntiLambda1point0 = new TH3F("f3fHistCentVsInvMassAntiLambda1point0","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0);
    
    f3fHistCentVsInvMassAntiLambda1point0Masscut = new TH3F("f3fHistCentVsInvMassAntiLambda1point0Masscut","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0Masscut);
    
    f3fHistCentVsInvMassAntiLambda1point0Sigtwo = new TH3F("f3fHistCentVsInvMassAntiLambda1point0Sigtwo","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0Sigtwo);
    
    
    f3fHistCentVsInvMassAntiLambda1point0SigtwoMasscut = new TH3F("f3fHistCentVsInvMassAntiLambda1point0SigtwoMasscut","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0SigtwoMasscut);
    
    
    f3fHistCentVsInvMassAntiLambda1point0Sigfour = new TH3F("f3fHistCentVsInvMassAntiLambda1point0Sigfour","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0Sigfour);
    
    f3fHistCentVsInvMassAntiLambda1point0SigfourMasscut = new TH3F("f3fHistCentVsInvMassAntiLambda1point0SigfourMasscut","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0SigfourMasscut);
    
    f3fHistCentVsInvMassAntiLambda1point0postight = new TH3F("f3fHistCentVsInvMassAntiLambda1point0postight","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0postight);
    
    f3fHistCentVsInvMassAntiLambda1point0postightMasscut = new TH3F("f3fHistCentVsInvMassAntiLambda1point0postightMasscut","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0postightMasscut);
    
    f3fHistCentVsInvMassAntiLambda1point0posloose = new TH3F("f3fHistCentVsInvMassAntiLambda1point0posloose","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0posloose);
    
    
    f3fHistCentVsInvMassAntiLambda1point0poslooseMasscut = new TH3F("f3fHistCentVsInvMassAntiLambda1point0poslooseMasscut","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0poslooseMasscut);
    
    f3fHistCentVsInvMassAntiLambda1point0negtight = new TH3F("f3fHistCentVsInvMassAntiLambda1point0negtight","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0negtight);
    
    
    f3fHistCentVsInvMassAntiLambda1point0negtightMasscut = new TH3F("f3fHistCentVsInvMassAntiLambda1point0negtightMasscut","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1.)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0negtightMasscut);
    
    f3fHistCentVsInvMassAntiLambda1point0negloose = new TH3F("f3fHistCentVsInvMassAntiLambda1point0negloose","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0negloose);
    
    f3fHistCentVsInvMassAntiLambda1point0neglooseMasscut = new TH3F("f3fHistCentVsInvMassAntiLambda1point0neglooseMasscut","Cent vs. #bar{#Lambda}Inv Mass vs. pT(deltaEta 1)",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0neglooseMasscut);
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //THNSPARSE BINNING
    
    const Int_t dim = 47; //23 pt bins*2 + 1 cent bin
    Int_t bin[dim];
    bin[0] = 81;
    for(Int_t ibin = 1; ibin < dim; ibin++) bin[ibin] = 500;
    Double_t min[dim];
    for(Int_t jbin = 0; jbin < dim; jbin++) min[jbin] =  -0.5;
    Double_t max[dim];
    max[0] = 80.5;
    for(Int_t jbin = 1; jbin < dim; jbin++) max[jbin] = 499.5;
    
    fPtBinNplusNminusChnegtightBKG = new THnSparseI("fPtBinNplusNminusChnegtightBKG","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChnegtightBKG);
    
    fPtBinNplusNminusChnegtight = new THnSparseI("fPtBinNplusNminusChnegtight","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChnegtight);
    
    
    fPtBinNplusNminusChneglooseBKG = new THnSparseI("fPtBinNplusNminusChneglooseBKG","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChneglooseBKG);
    
    
    fPtBinNplusNminusChnegloose = new THnSparseI("fPtBinNplusNminusChnegloose","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChnegloose);
    
    
    fPtBinNplusNminusChposlooseBKG = new THnSparseI("fPtBinNplusNminusChposlooseBKG","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChposlooseBKG);
    
    
    fPtBinNplusNminusChposloose = new THnSparseI("fPtBinNplusNminusChposloose","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChposloose);
    
    
    
    fPtBinNplusNminusChpostightBKG = new THnSparseI("fPtBinNplusNminusChpostightBKG","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChpostightBKG);
    
    
    fPtBinNplusNminusChpostight = new THnSparseI("fPtBinNplusNminusChpostight","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChpostight);
    
    
    
    fPtBinNplusNminusChSigfourBKG = new THnSparseI("fPtBinNplusNminusChSigfourBKG","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChSigfourBKG);
    
    
    
    fPtBinNplusNminusChSigfour = new THnSparseI("fPtBinNplusNminusChSigfour","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChSigfour);
    
    fPtBinNplusNminusChSigtwoBKG = new THnSparseI("fPtBinNplusNminusChSigtwoBKG","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChSigtwoBKG);
    
    
    fPtBinNplusNminusChSigtwo = new THnSparseI("fPtBinNplusNminusChSigtwo","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChSigtwo);
    
    fPtBinNplusNminusChBKG = new THnSparseI("fPtBinNplusNminusChBKG","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChBKG);
    
    
    fPtBinNplusNminusCh = new THnSparseI("fPtBinNplusNminusCh","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusCh);
    
    
    
    
    
    PostData(1,fListHist);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNetLambdaTrad::UserExec(Option_t *)
{
    
    const Int_t dim = fNptBins*2;
    
    Int_t ptChEta1point0[dim];
    Int_t ptChEta1point0LF[dim];
    Int_t ptChEta1point0Sigtwo[dim];
    Int_t ptChEta1point0LFSigtwo[dim];
    Int_t ptChEta1point0Sigfour[dim];
    Int_t ptChEta1point0LFSigfour[dim];
    Int_t ptChEta1point0postight[dim];
    Int_t ptChEta1point0LFpostight[dim];
    Int_t ptChEta1point0posloose[dim];
    Int_t ptChEta1point0LFposloose[dim];
    Int_t ptChEta1point0negtight[dim];
    Int_t ptChEta1point0LFnegtight[dim];
    Int_t ptChEta1point0negloose[dim];
    Int_t ptChEta1point0LFnegloose[dim];
    
    
    for(Int_t idx = 0; idx < dim; idx++)
    {
        ptChEta1point0[idx] = 0.;
        ptChEta1point0LF[idx] = 0.;
        ptChEta1point0Sigtwo[idx] = 0.;
        ptChEta1point0LFSigtwo[idx] = 0.;
        ptChEta1point0Sigfour[idx] = 0.;
        ptChEta1point0LFSigfour[idx] = 0.;
        ptChEta1point0postight[idx] = 0.;
        ptChEta1point0LFpostight[idx] = 0.;
        ptChEta1point0posloose[idx] = 0.;
        ptChEta1point0LFposloose[idx] = 0.;
        ptChEta1point0negtight[idx] = 0.;
        ptChEta1point0LFnegtight[idx] = 0.;
        ptChEta1point0negloose[idx] = 0.;
        ptChEta1point0LFnegloose[idx] = 0.;
        
    }
    
    if (!fInputEvent) return;
    
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) return;
    
    
    fPIDResponse = fInputHandler->GetPIDResponse();
    if(!fPIDResponse) return;
    
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
    
    Int_t nV0 = 0;
    Double_t fMinV0Pt = 0.5;
    Double_t fMaxV0Pt = 4.5;
    nV0 = fESD->GetNumberOfV0s();
    AliESDv0 *esdv0 = 0x0;
    AliESDtrack *esdpTrack = 0x0;
    AliESDtrack *esdnTrack = 0x0;
    
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
//        if(TMath::Abs(lRapLambda)> 0.5 ) continue;
        
        
        //--------------------------------------------------------------------Track selection-------------------------------------------------------------------------------------------------------------------------------------
        
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
        
        posprnsg = fPIDResponse->NumberOfSigmasTPC(esdpTrack, AliPID::kProton);
        negprnsg = fPIDResponse->NumberOfSigmasTPC(esdnTrack, AliPID::kProton);
        pospion  = fPIDResponse->NumberOfSigmasTPC( esdpTrack, AliPID::kPion );
        negpion  = fPIDResponse->NumberOfSigmasTPC( esdnTrack, AliPID::kPion );
        
        if(TMath::Abs(peta) > 0.8) continue;
        if(TMath::Abs(neta) > 0.8) continue;
        if(cosPointingAngle < 0.99) continue;
        if(dcaDaughters > 0.8) continue;
        if(v0Radius < 5.0) continue;
        if(v0Radius > 200.) continue;
        
        
        
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Int_t iptbin = GetPtBin(V0pt);
        
        if( ontheflystat == 0 )
        {
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex >  0.1  && TMath::Abs(posprnsg)  <= 3. && TMath::Abs(negpion)  <= 3.) //Default
            {
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassLambda1point0->Fill(fCentrality,invMassLambda,V0pt);
                    if(invMassLambda > 1.11 && invMassLambda < 1.122)
                    {
                        ptChEta1point0[iptbin] += 1;
                        f3fHistCentVsInvMassLambda1point0Masscut->Fill(fCentrality,invMassLambda,V0pt);
                    }
                    if(invMassLambda > 1.094 && invMassLambda < 1.107)
                    {
                        ptChEta1point0LF[iptbin] += 1;
                    }
                }
            }
            
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex >  0.1  && TMath::Abs(posprnsg)  <= 2.5 && TMath::Abs(negpion)  <= 2.5) //tight
            {
                
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassLambda1point0Sigtwo->Fill(fCentrality,invMassLambda,V0pt);
                    if(invMassLambda > 1.11 && invMassLambda < 1.122)
                    {
                        ptChEta1point0Sigtwo[iptbin] += 1;
                        f3fHistCentVsInvMassLambda1point0SigtwoMasscut->Fill(fCentrality,invMassLambda,V0pt);
                    }
                    if(invMassLambda > 1.094 && invMassLambda < 1.107)
                    {
                        ptChEta1point0LFSigtwo[iptbin] += 1;
                    }
                }
            }
            
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex >  0.1  && TMath::Abs(posprnsg)  <= 4 && TMath::Abs(negpion)  <= 4)//loose
            {
                
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassLambda1point0Sigfour->Fill(fCentrality,invMassLambda,V0pt);
                    if(invMassLambda > 1.11 && invMassLambda < 1.122)
                    {
                        ptChEta1point0Sigfour[iptbin] += 1;
                        f3fHistCentVsInvMassLambda1point0SigfourMasscut->Fill(fCentrality,invMassLambda,V0pt);
                    }
                    if(invMassLambda > 1.094 && invMassLambda < 1.106)
                    {
                        ptChEta1point0LFSigfour[iptbin] += 1;
                    }
                }
            }
            
            //L Pos to PV
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex >  0.13  && TMath::Abs(posprnsg)  <= 3 && TMath::Abs(negpion)  <= 3) //tight
            {
                
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassLambda1point0postight->Fill(fCentrality,invMassLambda,V0pt);
                    if(invMassLambda > 1.11 && invMassLambda < 1.122)
                    {
                        ptChEta1point0postight[iptbin] += 1;
                        f3fHistCentVsInvMassLambda1point0postightMasscut->Fill(fCentrality,invMassLambda,V0pt);
                    }
                    if(invMassLambda > 1.095 && invMassLambda < 1.108)
                    {
                        ptChEta1point0LFpostight[iptbin] += 1;
                    }
                }
            }
            
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex >  0.08  && TMath::Abs(posprnsg)  <= 3 && TMath::Abs(negpion)  <= 3) //loose
            {
                
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassLambda1point0posloose->Fill(fCentrality,invMassLambda,V0pt);
                    if(invMassLambda > 1.11 && invMassLambda < 1.122)
                    {
                        ptChEta1point0posloose[iptbin] += 1;
                        f3fHistCentVsInvMassLambda1point0poslooseMasscut->Fill(fCentrality,invMassLambda,V0pt);
                    }
                    if(invMassLambda > 1.096 && invMassLambda < 1.108)
                    {
                        ptChEta1point0LFposloose[iptbin] += 1;
                    }
                }
            }
            
            
            //L Neg to PV
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.3 && dcaPosToVertex >  0.1  && TMath::Abs(posprnsg)  <= 3 && TMath::Abs(negpion)  <= 3) //tight
            {
                
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassLambda1point0negtight->Fill(fCentrality,invMassLambda,V0pt);
                    if(invMassLambda > 1.11 && invMassLambda < 1.122)
                    {
                        ptChEta1point0negtight[iptbin] += 1;
                        f3fHistCentVsInvMassLambda1point0negtightMasscut->Fill(fCentrality,invMassLambda,V0pt);
                    }
                    if(invMassLambda > 1.095 && invMassLambda < 1.108)
                    {
                        ptChEta1point0LFnegtight[iptbin] += 1;
                    }
                }
            }
            
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.2 && dcaPosToVertex >  0.1  && TMath::Abs(posprnsg)  <= 3 && TMath::Abs(negpion)  <= 3) //loose
            {
                
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassLambda1point0negloose->Fill(fCentrality,invMassLambda,V0pt);
                    if(invMassLambda > 1.11 && invMassLambda < 1.122)
                    {
                        ptChEta1point0negloose[iptbin] += 1;
                        f3fHistCentVsInvMassLambda1point0neglooseMasscut->Fill(fCentrality,invMassLambda,V0pt);
                    }
                    if(invMassLambda > 1.096 && invMassLambda < 1.108)
                    {
                        ptChEta1point0LFnegloose[iptbin] += 1;
                    }
                }
            }
            
            // Bar-L
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.25 && TMath::Abs(negprnsg)  <= 3. && TMath::Abs(pospion)  <= 3.) //default
            {
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassAntiLambda1point0->Fill(fCentrality,invMassAntiLambda,V0pt);
                    if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.122)
                    {
                        ptChEta1point0[iptbin+fNptBins] += 1;
                        f3fHistCentVsInvMassAntiLambda1point0Masscut->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    if(invMassAntiLambda > 1.094 && invMassAntiLambda < 1.1075)
                    {
                        ptChEta1point0LF[iptbin+fNptBins] += 1;
                    }
                }
            }
            
            
            //bar-L nsig
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.25 && TMath::Abs(negprnsg)  <= 2.5 && TMath::Abs(pospion)  <= 2.5) //tight
            {
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassAntiLambda1point0Sigtwo->Fill(fCentrality,invMassAntiLambda,V0pt);
                    if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.122)
                    {
                        ptChEta1point0Sigtwo[iptbin+fNptBins] += 1;
                        f3fHistCentVsInvMassAntiLambda1point0SigtwoMasscut->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    if(invMassAntiLambda > 1.094 && invMassAntiLambda < 1.108)
                    {
                        ptChEta1point0LFSigtwo[iptbin+fNptBins] += 1;
                    }
                }
            }
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.25 && TMath::Abs(negprnsg)  <= 4 && TMath::Abs(pospion)  <= 4) //loose
            {
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassAntiLambda1point0Sigfour->Fill(fCentrality,invMassAntiLambda,V0pt);
                    if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.122)
                    {
                        ptChEta1point0Sigfour[iptbin+fNptBins] += 1;
                        f3fHistCentVsInvMassAntiLambda1point0SigfourMasscut->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    if(invMassAntiLambda > 1.094 && invMassAntiLambda < 1.1062)
                    {
                        ptChEta1point0LFSigfour[iptbin+fNptBins] += 1;
                    }
                }
            }
            
            
            //Bar-L Pos to PV
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.3 && TMath::Abs(negprnsg)  <= 3. && TMath::Abs(pospion)  <= 3.) //tight
            {
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassAntiLambda1point0postight->Fill(fCentrality,invMassAntiLambda,V0pt);
                    if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.122)
                    {
                        ptChEta1point0postight[iptbin+fNptBins] += 1;
                        f3fHistCentVsInvMassAntiLambda1point0postightMasscut->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    if(invMassAntiLambda > 1.095 && invMassAntiLambda < 1.108)
                    {
                        ptChEta1point0LFpostight[iptbin+fNptBins] += 1;
                    }
                }
            }
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.2 && TMath::Abs(negprnsg)  <= 3. && TMath::Abs(pospion)  <= 3.) //loose
            {
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassAntiLambda1point0posloose->Fill(fCentrality,invMassAntiLambda,V0pt);
                    if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.122)
                    {
                        ptChEta1point0posloose[iptbin+fNptBins] += 1;
                        f3fHistCentVsInvMassAntiLambda1point0poslooseMasscut->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    if(invMassAntiLambda > 1.096 && invMassAntiLambda < 1.108)
                    {
                        ptChEta1point0LFposloose[iptbin+fNptBins] += 1;
                    }
                }
            }
            //Bar-L Neg to PV
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.13 && dcaPosToVertex >  0.25 && TMath::Abs(negprnsg)  <= 3. && TMath::Abs(pospion)  <= 3.) //tight
            {
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassAntiLambda1point0negtight->Fill(fCentrality,invMassAntiLambda,V0pt);
                    if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.122)
                    {
                        ptChEta1point0negtight[iptbin+fNptBins] += 1;
                        f3fHistCentVsInvMassAntiLambda1point0negtightMasscut->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    if(invMassAntiLambda > 1.096 && invMassAntiLambda < 1.1088)
                    {
                        ptChEta1point0LFnegtight[iptbin+fNptBins] += 1;
                    }
                }
            }
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.08 && dcaPosToVertex >  0.25 && TMath::Abs(negprnsg)  <= 3. && TMath::Abs(pospion)  <= 3.) //loose
            {
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassAntiLambda1point0negloose->Fill(fCentrality,invMassAntiLambda,V0pt);
                    if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.122)
                    {
                        ptChEta1point0negloose[iptbin+fNptBins] += 1;
                        f3fHistCentVsInvMassAntiLambda1point0neglooseMasscut->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    if(invMassAntiLambda > 1.096 && invMassAntiLambda < 1.1088)
                    {
                        ptChEta1point0LFnegloose[iptbin+fNptBins] += 1;
                    }
                }
            }
            
        }// zero onfly V0
    }// end of V0 loop
    ///////////////////
    
    Double_t ptContainer[dim+1];
    ptContainer[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainer[i] = ptChEta1point0[i-1];
    }
    fPtBinNplusNminusCh->Fill(ptContainer);
    Double_t ptContainerBKG[dim+1];
    ptContainerBKG[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerBKG[i] = ptChEta1point0LF[i-1];
    }
    fPtBinNplusNminusChBKG->Fill(ptContainerBKG);
    
    
    /////////
    Double_t ptContainerSigtwo[dim+1];
    ptContainerSigtwo[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerSigtwo[i] = ptChEta1point0Sigtwo[i-1];
    }
    fPtBinNplusNminusChSigtwo->Fill(ptContainerSigtwo);
    Double_t ptContainerSigtwoBKG[dim+1];
    ptContainerSigtwoBKG[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerSigtwoBKG[i] = ptChEta1point0LFSigtwo[i-1];
    }
    fPtBinNplusNminusChSigtwoBKG->Fill(ptContainerSigtwoBKG);
    
    ////////
    Double_t ptContainerSigfour[dim+1];
    ptContainerSigfour[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerSigfour[i] = ptChEta1point0Sigfour[i-1];
    }
    fPtBinNplusNminusChSigfour->Fill(ptContainerSigfour);
    Double_t ptContainerSigfourBKG[dim+1];
    ptContainerSigfourBKG[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerSigfourBKG[i] = ptChEta1point0LFSigfour[i-1];
    }
    fPtBinNplusNminusChSigfourBKG->Fill(ptContainerSigfourBKG);
    
    /////////////
    Double_t ptContainerpostight[dim+1];
    ptContainerpostight[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerpostight[i] = ptChEta1point0postight[i-1];
    }
    fPtBinNplusNminusChpostight->Fill(ptContainerpostight);
    Double_t ptContainerpostightBKG[dim+1];
    ptContainerpostightBKG[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerpostightBKG[i] = ptChEta1point0LFpostight[i-1];
    }
    fPtBinNplusNminusChpostightBKG->Fill(ptContainerpostightBKG);
    
    ///////
    Double_t ptContainerposloose[dim+1];
    ptContainerposloose[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerposloose[i] = ptChEta1point0posloose[i-1];
    }
    fPtBinNplusNminusChposloose->Fill(ptContainerposloose);
    Double_t ptContainerposlooseBKG[dim+1];
    ptContainerposlooseBKG[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerposlooseBKG[i] = ptChEta1point0LFposloose[i-1];
    }
    fPtBinNplusNminusChposlooseBKG->Fill(ptContainerposlooseBKG);
    
    //////
    
    Double_t ptContainernegloose[dim+1];
    ptContainernegloose[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainernegloose[i] = ptChEta1point0negloose[i-1];
    }
    fPtBinNplusNminusChnegloose->Fill(ptContainernegloose);
    ///////////
    Double_t ptContainerneglooseBKG[dim+1];
    ptContainerneglooseBKG[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerneglooseBKG[i] = ptChEta1point0LFnegloose[i-1];
    }
    fPtBinNplusNminusChneglooseBKG->Fill(ptContainerneglooseBKG);
    
    ///////////////
    
    Double_t ptContainernegtight[dim+1];
    ptContainernegtight[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainernegtight[i] = ptChEta1point0negtight[i-1];
    }
    fPtBinNplusNminusChnegtight->Fill(ptContainernegtight);
    
    Double_t ptContainernegtightBKG[dim+1];
    ptContainernegtightBKG[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainernegtightBKG[i] = ptChEta1point0LFnegtight[i-1];
    }
    fPtBinNplusNminusChnegtightBKG->Fill(ptContainernegtightBKG);
    
    /////////////
    
    PostData(1,fListHist);
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskNetLambdaTrad::GetPtBin(Double_t pt)
{
    Int_t bin = -1;
    
    Double_t LambdaPtBins[24] = {0.8,0.9,1.0,1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,4.2};
    
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
