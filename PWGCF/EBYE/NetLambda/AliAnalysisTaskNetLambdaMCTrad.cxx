// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Naomi Umaka Apr 2018
// email: ejiro.naomi.umaka@cern.ch
//Apr 2 syst.


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


#include "AliAnalysisTaskNetLambdaMCTrad.h"

ClassImp(AliAnalysisTaskNetLambdaMCTrad)

//-----------------------------------------------------------------------------
AliAnalysisTaskNetLambdaMCTrad::AliAnalysisTaskNetLambdaMCTrad(const char* name) :
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


f2fHistRecPrimariesCentVsPtLambdaFourSigthree(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree(0x0),
f3fHistLambdafromXiFourSigthree(0x0),
f3fHistAntiLambdafromXiFourSigthree(0x0),

f2fHistRecPrimariesCentVsPtLambdaFourSigthreensigtight(0x0),
f3fHistLambdafromXiFourSigthreensigtight(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreensigtight(0x0),
f3fHistAntiLambdafromXiFourSigthreensigtight(0x0),

f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegloose(0x0),
f3fHistLambdafromXiFourSigthreenegloose(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegloose(0x0),
f3fHistAntiLambdafromXiFourSigthreenegloose(0x0),

f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegtight(0x0),
f3fHistLambdafromXiFourSigthreenegtight(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegtight(0x0),
f3fHistAntiLambdafromXiFourSigthreenegtight(0x0),

f2fHistRecPrimariesCentVsPtLambdaFourSigthreeposloose(0x0),
f3fHistLambdafromXiFourSigthreeposloose(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreeposloose(0x0),
f3fHistAntiLambdafromXiFourSigthreeposloose(0x0),

f2fHistRecPrimariesCentVsPtLambdaFourSigthreepostight(0x0),
f3fHistLambdafromXiFourSigthreepostight(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreepostight(0x0),
f3fHistAntiLambdafromXiFourSigthreepostight(0x0),

fCentrality(-1),
fTreeVariablePID(-1),

fNptBins(23),
fIsMC(kTRUE),
fEvSel(AliVEvent::kINT7),


fPtBinNplusNminusChnsigtight(NULL),
fPtBinNplusNminusChnegloose(NULL),
fPtBinNplusNminusChnegtight(NULL),
fPtBinNplusNminusChposloose(NULL),
fPtBinNplusNminusChpostight(NULL),
fPtBinNplusNminusCh(NULL),
fPtBinNplusNminusChTruth(NULL)

{
    Info("AliAnalysisTaskNetLambdaMCTrad","Calling Constructor");
    
    DefineInput(0,TChain::Class());
    DefineOutput(1,TList::Class());
    
}
//---------------------------------------------------------------------------------------------------------------------------------------

AliAnalysisTaskNetLambdaMCTrad::~AliAnalysisTaskNetLambdaMCTrad()
{
    // Default destructor
    if( fListHist ) delete fListHist;
}
//----------------------------------------------------------------------------------------------------

void AliAnalysisTaskNetLambdaMCTrad::UserCreateOutputObjects()
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
    Double_t CentBins[101] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100};
    Long_t CentbinNum = sizeof(CentBins)/sizeof(Double_t) - 1;
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     Double_t LambdaPtBins[24] = {0.9,1.0,1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,4.2, 4.4};
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Double_t xibinlimits[26] = {0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.5,6.0,7.0,8.0};
    Long_t xibinnumb = sizeof(xibinlimits)/sizeof(Double_t) - 1;
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Double_t MassBins[103]
    = {1.0788,1.0796,1.0804,1.0812,1.082,1.0828,1.0836,1.0844,1.0852,1.086,1.0868,1.0876,1.0884,1.0892,1.09,1.0908,1.0916,1.0924,1.0932,1.094,1.0948,1.0956,1.0964,1.0972,1.098,1.0988,1.0996,1.1004,1.1012,1.102,
        1.1028,1.1036,1.1044,1.1052,1.106,1.1068,1.1076,1.1084,1.1092,1.11,1.1108,1.1116,1.1124,1.1132,1.114,1.1148,1.1156,1.1164,1.1172,1.118,1.1188,1.1196,1.1204,1.1212,1.122,1.1228,1.1236,1.1244,
        1.1252,1.126,1.1268,1.1276,1.1284,1.1292,1.13,1.1308,1.1316,1.1324,1.1332,1.134,1.1348,1.1356,1.1364,1.1372,1.138,1.1388,1.1396,1.1404,1.1412,1.142,1.1428,1.1436,1.1444,1.1452,1.146,1.1468,1.1476,
        1.1484,1.1492,1.15,1.1508,1.1516,1.1524,1.1532,1.154,1.1548,1.1556,1.1564,1.1572,1.158,1.1588,1.1596,1.1604};
    Long_t Massbinnumb = sizeof(MassBins)/sizeof(Double_t) - 1;
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
    
    
    if(fIsMC)
    {
        //--------------------------------------------------------THNSPARSE----------------------------------------------------------------------------------------------------------------
        
        fPtBinNplusNminusChTruth = new THnSparseI("fPtBinNplusNminusChTruth","fPtBinNplusNminusChTruth", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTruth); //gen
        
        fPtBinNplusNminusCh = new THnSparseI("fPtBinNplusNminusCh","fPtBinNplusNminusCh", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusCh);
        
        fPtBinNplusNminusChnsigtight = new THnSparseI("fPtBinNplusNminusChnsigtight","fPtBinNplusNminusChnsigtight", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChnsigtight);
        
        fPtBinNplusNminusChnegloose = new THnSparseI("fPtBinNplusNminusChnegloose","fPtBinNplusNminusChnegloose", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChnegloose);
        
        fPtBinNplusNminusChnegtight = new THnSparseI("fPtBinNplusNminusChnegtight","fPtBinNplusNminusChnegtight", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChnegtight);
        
        fPtBinNplusNminusChposloose = new THnSparseI("fPtBinNplusNminusChposloose","fPtBinNplusNminusChposloose", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChposloose);
        
        fPtBinNplusNminusChpostight = new THnSparseI("fPtBinNplusNminusChpostight","fPtBinNplusNminusChpostight", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChpostight);
        
        
        
        //Gen
        
        f2fHistGenCentVsPtLambda = new TH2F( "f2fHistGenCentVsPtLambda", "Centrality Vs #Lambda Gen Pt",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistGenCentVsPtLambda);
        
        f2fHistGenCentVsPtAntiLambda = new TH2F( "f2fHistGenCentVsPtAntiLambda", "Centrality Vs #bar{#Lambda} Gen Pt",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistGenCentVsPtAntiLambda);
        
        f2fHistXiPlus = new TH2F("f2fHistXiPlus","f2fHistXiPlus",CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistXiPlus);
        
        f2fHistXiMinus = new TH2F("f2fHistXiMinus","f2fHistXiMinus",CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistXiMinus);
        ///
        f2fHistRecPrimariesCentVsPtLambdaFourSigthree = new TH2F("f2fHistRecPrimariesCentVsPtLambdaFourSigthree","f2fHistRecPrimariesCentVsPtLambdaFourSigthree",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaFourSigthree);
        
        f3fHistLambdafromXiFourSigthree = new TH3F("f3fHistLambdafromXiFourSigthree","f3fHistLambdafromXiFourSigthree ", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiFourSigthree);
        
        f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree","f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree);
        
        f3fHistAntiLambdafromXiFourSigthree = new TH3F("f3fHistAntiLambdafromXiFourSigthree","f3fHistAntiLambdafromXiFourSigthree ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiFourSigthree);
        ///
        f2fHistRecPrimariesCentVsPtLambdaFourSigthreensigtight = new TH2F("f2fHistRecPrimariesCentVsPtLambdaFourSigthreensigtight","f2fHistRecPrimariesCentVsPtLambdaFourSigthreensigtight",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaFourSigthreensigtight);
        
        f3fHistLambdafromXiFourSigthreensigtight = new TH3F("f3fHistLambdafromXiFourSigthreensigtight","f3fHistLambdafromXiFourSigthreensigtight ", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiFourSigthreensigtight);
        ///
        f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreensigtight = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreensigtight","f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreensigtight",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreensigtight);
        
        f3fHistAntiLambdafromXiFourSigthreensigtight = new TH3F("f3fHistAntiLambdafromXiFourSigthreensigtight","f3fHistAntiLambdafromXiFourSigthreensigtight ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiFourSigthreensigtight);
        ///
        f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegloose = new TH2F("f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegloose","f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegloose",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegloose);
        
        f3fHistLambdafromXiFourSigthreenegloose = new TH3F("f3fHistLambdafromXiFourSigthreenegloose","f3fHistLambdafromXiFourSigthreenegloose ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiFourSigthreenegloose);
        ///
        f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegloose = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegloose","f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegloose",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegloose);
        
        f3fHistAntiLambdafromXiFourSigthreenegloose = new TH3F("f3fHistAntiLambdafromXiFourSigthreenegloose","f3fHistAntiLambdafromXiFourSigthreenegloose ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiFourSigthreenegloose);
        ///
        f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegtight = new TH2F("f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegtight","f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegtight",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegtight);
        
        f3fHistLambdafromXiFourSigthreenegtight = new TH3F("f3fHistLambdafromXiFourSigthreenegtight","f3fHistLambdafromXiFourSigthreenegtight ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiFourSigthreenegtight);
        ///
        f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegtight = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegtight","f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegtight",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegtight);
        
        f3fHistAntiLambdafromXiFourSigthreenegtight = new TH3F("f3fHistAntiLambdafromXiFourSigthreenegtight","f3fHistAntiLambdafromXiFourSigthreenegtight ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiFourSigthreenegtight);
      ///
        f2fHistRecPrimariesCentVsPtLambdaFourSigthreeposloose = new TH2F("f2fHistRecPrimariesCentVsPtLambdaFourSigthreeposloose","f2fHistRecPrimariesCentVsPtLambdaFourSigthreeposloose",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaFourSigthreeposloose);
        
        f3fHistLambdafromXiFourSigthreeposloose = new TH3F("f3fHistLambdafromXiFourSigthreeposloose","f3fHistLambdafromXiFourSigthreeposloose ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiFourSigthreeposloose);
        ///
        f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreeposloose = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreeposloose","f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreeposloose",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreeposloose);
        
        f3fHistAntiLambdafromXiFourSigthreeposloose = new TH3F("f3fHistAntiLambdafromXiFourSigthreeposloose","f3fHistAntiLambdafromXiFourSigthreeposloose ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiFourSigthreeposloose);
        ///
        f2fHistRecPrimariesCentVsPtLambdaFourSigthreepostight = new TH2F("f2fHistRecPrimariesCentVsPtLambdaFourSigthreepostight","f2fHistRecPrimariesCentVsPtLambdaFourSigthreepostight",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaFourSigthreepostight);
        
        f3fHistLambdafromXiFourSigthreepostight = new TH3F("f3fHistLambdafromXiFourSigthreepostight","f3fHistLambdafromXiFourSigthreepostight ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiFourSigthreepostight);
        ///
        f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreepostight = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreepostight","f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreepostight",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreepostight);
        
        f3fHistAntiLambdafromXiFourSigthreepostight = new TH3F("f3fHistAntiLambdafromXiFourSigthreepostight","f3fHistAntiLambdafromXiFourSigthreepostight ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiFourSigthreepostight);
       

        
    }
    
    PostData(1,fListHist);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNetLambdaMCTrad::UserExec(Option_t *)
{
    
    const Int_t dim = fNptBins*2;
    
    Int_t ptChMC[dim];
    Int_t ptChEta1point0[dim];
    Int_t ptChEta1point0nsigtight[dim];
    Int_t ptChEta1point0negloose[dim];
    Int_t ptChEta1point0negtight[dim];
    Int_t ptChEta1point0posloose[dim];
    Int_t ptChEta1point0postight[dim];

    
    for(Int_t idx = 0; idx < dim; idx++)
    {
        ptChMC[idx] = 0;
        ptChEta1point0[idx] = 0.;
        ptChEta1point0nsigtight[idx] = 0.;
        ptChEta1point0negloose[idx] = 0.;
        ptChEta1point0negtight[idx] = 0.;
        ptChEta1point0posloose[idx] = 0.;
        ptChEta1point0postight[idx] = 0.;
   
    }
    
    
    if (!fInputEvent) return;
    
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) return;
    
    fPIDResponse = fInputHandler->GetPIDResponse();
    if(!fPIDResponse) return;
    
    if(!(fInputHandler->IsEventSelected() & fEvSel)) return;
    
    
    
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    if(!mcH) return;
    fMCEvent=mcH->MCEvent();
    if(!fMCEvent) return;
    
    
    
    Double_t lMagneticField = -10;
    lMagneticField = fESD->GetMagneticField();
    
    AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent->FindListObject("MultSelection");
    if(!MultSelection) return;
    
    
    Double_t vVtx[3];
    
    AliVVertex *vvertex = (AliVVertex*)fInputEvent->GetPrimaryVertex();
    if (!vvertex) return;
    vVtx[0] = vvertex->GetX();
    vVtx[1] = vvertex->GetY();
    vVtx[2] = vvertex->GetZ();
    
    if(vVtx[2] < -10. || vVtx[2] > 10.) return;
    
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    
    if( fCentrality < 0 || fCentrality > 80 ) return;
    if (!fEventCuts.AcceptEvent(fInputEvent)) return;//pileup cut
    
    fHistEventCounter->Fill(1.5);
    fHistCentrality->Fill(fCentrality);
    
    const AliESDVertex *lPrimaryBestESDVtx     = fESD->GetPrimaryVertex();
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    
    
    Int_t nGen = 0;
    if(fIsMC)
    {
        
        nGen = fMCEvent->GetNumberOfTracks();
        
        for(Int_t iGen = 0; iGen < nGen; iGen++)
        {
            Int_t genpid = -1;
            Double_t lThisRap  = 0;
            Float_t gpt = 0.0, eta = 0.0,abseta =0.0;
            
            AliMCParticle* mctrack = (AliMCParticle*)fMCEvent->GetTrack(iGen);
            if(!mctrack) continue;
            if(!(fMCEvent->IsPhysicalPrimary(iGen))) continue;
            
            TParticle *part = mctrack->Particle();
            genpid = part->GetPdgCode();
            gpt = part->Pt();
            eta = part->Eta();
            abseta = TMath::Abs(eta);
            lThisRap   = MyRapidity(part->Energy(),part->Pz());
            
            
            Int_t iptbinMC = GetPtBin(gpt);
            if( iptbinMC < 0 || iptbinMC > fNptBins-1 ) continue;

            
            if(abseta < 0.5)
            {
                if(genpid == 3122)
                {
                    f2fHistGenCentVsPtLambda->Fill(fCentrality,gpt);
                    ptChMC[iptbinMC] += 1;
                }
                
                if(genpid == -3122)
                {
                    f2fHistGenCentVsPtAntiLambda->Fill(fCentrality,gpt);
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
            }
            
        } // end loop over generated particles
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
        Float_t ppt = -999,  peta = -999, posprnsg = -999, pospion =-999, v0Radius =-999, v0DecayLength =-999, v0dlength = -999,proLT =-999, proLTbar =-999, lRapLambda=-999;
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
        v0DecayLength = TMath::Sqrt(TMath::Power(vertx[0] - vVtx[0],2) +
                                    TMath::Power(vertx[1] - vVtx[1],2) +
                                    TMath::Power(vertx[2] - vVtx[2],2 ));
        v0dlength= TMath::Sqrt(TMath::Power(vertx[0] - vVtx[0],2) +
                               TMath::Power(vertx[1] - vVtx[1],2) +
                               TMath::Power(vertx[2] - vVtx[2],2 ));
        
        lRapLambda  = esdv0->RapLambda();
        V0pt = esdv0->Pt();
        if ((V0pt<fMinV0Pt)||(fMaxV0Pt<V0pt)) continue;
        
        Double_t tV0mom[3];
        esdv0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] );
        Double_t lV0TotalMomentum = TMath::Sqrt(
                                                tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );
        v0DecayLength /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure
        //--------------------------------------------------------------------Track selection-------------------------------------------------------------------------
        Float_t lPosTrackCrossedRows = esdpTrack->GetTPCClusterInfo(2,1);
        Float_t lNegTrackCrossedRows = esdnTrack->GetTPCClusterInfo(2,1);
        
        fTreeVariableLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
        if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
            fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
        
        if( !(esdpTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        if( !(esdnTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        
        //Extra track quality: min track length
        Float_t lSmallestTrackLength = 1000;
        Float_t lPosTrackLength = -1;
        Float_t lNegTrackLength = -1;
        
        if (esdpTrack->GetInnerParam()) lPosTrackLength = esdpTrack->GetLengthInActiveZone(1, 2.0, 220.0, fESD->GetMagneticField());
        if (esdnTrack->GetInnerParam()) lNegTrackLength = esdnTrack->GetLengthInActiveZone(1, 2.0, 220.0, fESD->GetMagneticField());
        
        if ( lPosTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lPosTrackLength;
        if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
        
        if ( ( ( ( esdpTrack->GetTPCClusterInfo(2,1) ) < 80 ) || ( ( esdnTrack->GetTPCClusterInfo(2,1) ) < 80 ) ) && lSmallestTrackLength < 90 ) continue;
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
        //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        proLT = v0DecayLength*invMassLambda;
        proLTbar = v0DecayLength*invMassAntiLambda;
        
        
        if(TMath::Abs(peta) > 0.8) continue;
        if(TMath::Abs(neta) > 0.8) continue;
        if(dcaDaughters > 0.8) continue;
        if(v0Radius < 5.0) continue;
        if(cosPointingAngle < 0.99) continue;
        
        if( ontheflystat == 0)
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
                
                
                AliMCParticle *esdGenTrackPos = (AliMCParticle*)fMCEvent->GetTrack(lblPosV0Dghter);
                if(!esdGenTrackPos) continue;
                AliMCParticle *esdGenTrackNeg = (AliMCParticle*)fMCEvent->GetTrack(lblNegV0Dghter);
                if(!esdGenTrackNeg) continue;
                
                
                Int_t posTparticle = esdGenTrackPos->GetMother();
                Int_t negTparticle = esdGenTrackNeg->GetMother();
                
                if( posTparticle == negTparticle && posTparticle > 0 )
                {
                    AliMCParticle *esdlthisV0 = (AliMCParticle*)fMCEvent->GetTrack(posTparticle);
                    if(!esdlthisV0) continue;
                    TParticle *partRecMom = esdlthisV0->Particle();
                    fTreeVariablePID = partRecMom->GetPdgCode();
                    
                    mcpt = partRecMom->Pt();
                    mceta = partRecMom->Eta();
                    
                    isSecFromMaterial = fMCEvent->IsSecondaryFromMaterial(posTparticle);
                    isSecFromWeakDecay = fMCEvent->IsSecondaryFromWeakDecay(posTparticle);
                    isPrim = fMCEvent->IsPhysicalPrimary(posTparticle);
                    
                    
                    Int_t esdlthisV0parent = esdlthisV0->GetMother();
                    
                    if (esdlthisV0parent > 0)
                    {
                        AliMCParticle *lbV0parent = (AliMCParticle*)fMCEvent->GetTrack(esdlthisV0parent);
                        
                        if(!lbV0parent) continue;
                        TParticle *partRecGMom = lbV0parent->Particle();
                        fTreeVariablePIDParent = partRecGMom->GetPdgCode();
                        fTreeVariablePtParent = partRecGMom->Pt();
                        isPrimParent =  fMCEvent->IsPhysicalPrimary(esdlthisV0parent);
                        
                    }
                }
                
                Int_t iptbin = GetPtBin(mcpt);
                if( iptbin < 0 || iptbin > fNptBins-1 ) continue;
                
                if(TMath::Abs(eta) < 0.5)
                {
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex > 0.1 && TMath::Abs(posprnsg) <= 4 && TMath::Abs(negpion) <= 4)
                    {
                      if(fTreeVariablePID == 3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtLambdaFourSigthree->Fill(fCentrality,mcpt);
                                ptChEta1point0[iptbin] += 1;
                            }
                            
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == 3312) || (fTreeVariablePIDParent == 3322))
                                    {
                                        f3fHistLambdafromXiFourSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                    
                   if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex > 0.25 && TMath::Abs(negprnsg) <= 4 && TMath::Abs(pospion) <= 4)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree->Fill(fCentrality,mcpt);
                                 ptChEta1point0[iptbin+fNptBins] += 1;
                            }
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == -3312) || (fTreeVariablePIDParent == -3322))
                                        
                                    {
                                        f3fHistAntiLambdafromXiFourSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex > 0.1 && TMath::Abs(posprnsg) <= 2.5 && TMath::Abs(negpion) <= 2.5)
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtLambdaFourSigthreensigtight->Fill(fCentrality,mcpt);
                                ptChEta1point0nsigtight[iptbin] += 1;
                            }
                            
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == 3312) || (fTreeVariablePIDParent == 3322))
                                    {
                                        f3fHistLambdafromXiFourSigthreensigtight->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                  if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex > 0.25 && TMath::Abs(negprnsg) <= 2.5 && TMath::Abs(pospion) <= 2.5)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreensigtight->Fill(fCentrality,mcpt);
                                ptChEta1point0nsigtight[iptbin+fNptBins] += 1;
                            }
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == -3312) || (fTreeVariablePIDParent == -3322))
                                        
                                    {
                                        f3fHistAntiLambdafromXiFourSigthreensigtight->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.2 && dcaPosToVertex > 0.1 && TMath::Abs(posprnsg) <= 3 && TMath::Abs(negpion) <= 3)
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegloose->Fill(fCentrality,mcpt);
                                ptChEta1point0negloose[iptbin] += 1;
                            }
                            
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == 3312) || (fTreeVariablePIDParent == 3322))
                                    {
                                        f3fHistLambdafromXiFourSigthreenegloose->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.08 && dcaPosToVertex > 0.25 && TMath::Abs(negprnsg) <= 3 && TMath::Abs(pospion) <= 3)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegloose->Fill(fCentrality,mcpt);
                                ptChEta1point0negloose[iptbin+fNptBins] += 1;
                            }
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == -3312) || (fTreeVariablePIDParent == -3322))
                                        
                                    {
                                        f3fHistAntiLambdafromXiFourSigthreenegloose->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                  if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.3 && dcaPosToVertex > 0.1 && TMath::Abs(posprnsg) <= 3 && TMath::Abs(negpion) <= 3)
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegtight->Fill(fCentrality,mcpt);
                                ptChEta1point0negtight[iptbin] += 1;
                            }
                            
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == 3312) || (fTreeVariablePIDParent == 3322))
                                    {
                                        f3fHistLambdafromXiFourSigthreenegtight->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.13 && dcaPosToVertex > 0.25 && TMath::Abs(negprnsg) <= 3 && TMath::Abs(pospion) <= 3)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegtight->Fill(fCentrality,mcpt);
                                ptChEta1point0negtight[iptbin+fNptBins] += 1;
                            }
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == -3312) || (fTreeVariablePIDParent == -3322))
                                        
                                    {
                                        f3fHistAntiLambdafromXiFourSigthreenegtight->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex > 0.08 && TMath::Abs(posprnsg) <= 3 && TMath::Abs(negpion) <= 3)
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtLambdaFourSigthreeposloose->Fill(fCentrality,mcpt);
                                ptChEta1point0posloose[iptbin] += 1;
                            }
                            
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == 3312) || (fTreeVariablePIDParent == 3322))
                                    {
                                        f3fHistLambdafromXiFourSigthreeposloose->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                   if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex > 0.2 && TMath::Abs(negprnsg) <= 3 && TMath::Abs(pospion) <= 3)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreeposloose->Fill(fCentrality,mcpt);
                                ptChEta1point0posloose[iptbin+fNptBins] += 1;
                            }
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == -3312) || (fTreeVariablePIDParent == -3322))
                                        
                                    {
                                        f3fHistAntiLambdafromXiFourSigthreeposloose->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                  if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex > 0.13 && TMath::Abs(posprnsg) <= 3 && TMath::Abs(negpion) <= 3)
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtLambdaFourSigthreepostight->Fill(fCentrality,mcpt);
                                ptChEta1point0postight[iptbin] += 1;
                            }
                            
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == 3312) || (fTreeVariablePIDParent == 3322))
                                    {
                                        f3fHistLambdafromXiFourSigthreepostight->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                   if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex > 0.3 && TMath::Abs(negprnsg) <= 3 && TMath::Abs(pospion) <= 3)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreepostight->Fill(fCentrality,mcpt);
                                ptChEta1point0postight[iptbin+fNptBins] += 1;
                            }
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == -3312) || (fTreeVariablePIDParent == -3322))
                                        
                                    {
                                        f3fHistAntiLambdafromXiFourSigthreepostight->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                    
                    
                }//|eta|<0.5
            } //MC condition
        }// zero onfly V0
    }// end of V0 loop
    
    
    //-------------------------------------------------
    ///
    Double_t ptContainer[dim+1];
    ptContainer[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainer[i] = ptChEta1point0[i-1];
    }
    fPtBinNplusNminusCh->Fill(ptContainer); //nsig4
    ///
    Double_t ptContainer2[dim+1];
    ptContainer2[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainer2[i] = ptChEta1point0nsigtight[i-1];
    }
    fPtBinNplusNminusChnsigtight->Fill(ptContainer2); //nsig2.5
    ///
    Double_t ptContainer3[dim+1];
    ptContainer3[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainer3[i] = ptChEta1point0negloose[i-1];
    }
    fPtBinNplusNminusChnegloose->Fill(ptContainer3); //negloose
    ///
    Double_t ptContainer4[dim+1];
    ptContainer4[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainer4[i] = ptChEta1point0negtight[i-1];
    }
    fPtBinNplusNminusChnegtight->Fill(ptContainer4); //negtight
    ///
    Double_t ptContainer5[dim+1];
    ptContainer5[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainer5[i] = ptChEta1point0posloose[i-1];
    }
    fPtBinNplusNminusChposloose->Fill(ptContainer5); //posloose
    ///
    Double_t ptContainer6[dim+1];
    ptContainer6[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainer6[i] = ptChEta1point0postight[i-1];
    }
    fPtBinNplusNminusChpostight->Fill(ptContainer6); //postight

    
    
    PostData(1,fListHist);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskNetLambdaMCTrad::GetPtBin(Double_t pt)
{
    Int_t bin = -1;
    
    Double_t LambdaPtBins[24] = {0.9,1.0,1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,4.2,4.4};
    
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

Double_t AliAnalysisTaskNetLambdaMCTrad::MyRapidity(Double_t rE, Double_t rPz) const
{
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}
