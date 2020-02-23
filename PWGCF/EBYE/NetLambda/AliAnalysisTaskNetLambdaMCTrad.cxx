// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Naomi Umaka Apr 2018
// email: ejiro.naomi.umaka@cern.ch
// Updated Feb 22


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

f2fHistGenCentVsPtLambdacopy(0x0),
f2fHistGenCentVsPtAntiLambdacopy(0x0),
f2fHistXiPluscopy(0x0),
f2fHistXiMinuscopy(0x0),


f2fHistRecSecCentVsPtLambdaFourSigthree(0x0),
f2fHistRecSecCentVsPtAntiLambdaFourSigthree(0x0),
f2fHistRecMatCentVsPtLambdaFourSigthree(0x0),
f2fHistRecMatCentVsPtAntiLambdaFourSigthree(0x0),


f2fHistRecPrimariesCentVsPtLambdaFourSigthree(0x0),
f2fHistRecPrimariesCentVsPtLambdaFourSigthreetight(0x0),

f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree(0x0),
f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreetight(0x0),

f3fHistLambdafromXiFourSigthree(0x0),
f3fHistAntiLambdafromXiFourSigthree(0x0),
f3fHistLambdafromXiFourSigthreetight(0x0),
f3fHistAntiLambdafromXiFourSigthreetight(0x0),

f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntag(0x0),
f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntagCut(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntag(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntagCut(0x0),

f3fHistCentInvMassVsPtLambdaRecFourSigthree(0x0),
f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree(0x0),

fCentrality(-1),
fTreeVariablePID(-1),

fNptBins(23),
fIsMC(kTRUE),
fEvSel(AliVEvent::kINT7),


fPtBinNplusNminusChRec(NULL),
fPtBinNplusNminusChRecTag(NULL),
fPtBinNplusNminusChTruth(NULL)

{
    Info("AliAnalysisTaskNetLambdaMCTrad","Calling Constructor");
    
    DefineInput(0,TChain::Class());
    DefineOutput(1,TList::Class());
    
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
        
        fPtBinNplusNminusChRec = new THnSparseI("fPtBinNplusNminusChRec","fPtBinNplusNminusChRec", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChRec); //
        
        fPtBinNplusNminusChRecTag = new THnSparseI("fPtBinNplusNminusChRecTag","fPtBinNplusNminusChRecTag", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChRecTag); //
        
        
        
        //-------------------------------------------------------------------GEN-----------------------------------------------------------------------------------------------
        
        f2fHistGenCentVsPtLambda = new TH2F( "f2fHistGenCentVsPtLambda", "Centrality Vs #Lambda Gen Pt",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistGenCentVsPtLambda);
        
        f2fHistGenCentVsPtAntiLambda = new TH2F( "f2fHistGenCentVsPtAntiLambda", "Centrality Vs #bar{#Lambda} Gen Pt",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistGenCentVsPtAntiLambda);
        
        f2fHistXiPlus = new TH2F("f2fHistXiPlus","f2fHistXiPlus",CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistXiPlus);
        
        f2fHistXiMinus = new TH2F("f2fHistXiMinus","f2fHistXiMinus",CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistXiMinus);
        
        
        f2fHistGenCentVsPtLambdacopy = new TH2F( "f2fHistGenCentVsPtLambdacopy", "f2fHistGenCentVsPtLambdacopy",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistGenCentVsPtLambdacopy);
        
        f2fHistGenCentVsPtAntiLambdacopy = new TH2F( "f2fHistGenCentVsPtAntiLambdacopy", "f2fHistGenCentVsPtAntiLambdacopy",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistGenCentVsPtAntiLambdacopy);
        
        f2fHistXiPluscopy = new TH2F("f2fHistXiPluscopy","f2fHistXiPluscopy",CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistXiPluscopy);
        
        f2fHistXiMinuscopy = new TH2F("f2fHistXiMinuscopy","f2fHistXiMinuscopy",CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistXiMinuscopy);
        
        //-------------------------------------------------------------------MC REC-----------------------------------------------------------------------------------------------
        
        //Sec
        f2fHistRecSecCentVsPtLambdaFourSigthree = new TH2F("f2fHistRecSecCentVsPtLambdaFourSigthree","#Lambda SEC  ",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecSecCentVsPtLambdaFourSigthree);
        
        f2fHistRecSecCentVsPtAntiLambdaFourSigthree = new TH2F("f2fHistRecSecCentVsPtAntiLambdaFourSigthree","#bar{#Lambda}  SEC ",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecSecCentVsPtAntiLambdaFourSigthree);
        
        //Mat
        f2fHistRecMatCentVsPtLambdaFourSigthree = new TH2F("f2fHistRecMatCentVsPtLambdaFourSigthree","#Lambda Mat  ",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecMatCentVsPtLambdaFourSigthree);
        
        f2fHistRecMatCentVsPtAntiLambdaFourSigthree = new TH2F("f2fHistRecMatCentVsPtAntiLambdaFourSigthree","#bar{#Lambda}  Mat ",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecMatCentVsPtAntiLambdaFourSigthree);
        
        
        //REC PRI
        
        f2fHistRecPrimariesCentVsPtLambdaFourSigthree = new TH2F("f2fHistRecPrimariesCentVsPtLambdaFourSigthree","#Lambda primaries",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaFourSigthree);
        
        //---
        f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree","#bar{#Lambda} primaries",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree);
        
        f2fHistRecPrimariesCentVsPtLambdaFourSigthreetight = new TH2F("f2fHistRecPrimariesCentVsPtLambdaFourSigthreetight","#Lambda primaries",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambdaFourSigthreetight);
        
        //---
        f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreetight = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreetight","#bar{#Lambda} primaries",CentbinNum, CentBins,fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreetight);
        
        //FD
        
        f3fHistLambdafromXiFourSigthree = new TH3F("f3fHistLambdafromXiFourSigthree","f3fHistLambdafromXiFourSigthree ", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiFourSigthree);
        
        //-
        f3fHistAntiLambdafromXiFourSigthree = new TH3F("f3fHistAntiLambdafromXiFourSigthree","f3fHistAntiLambdafromXiFourSigthree ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiFourSigthree);
        
        f3fHistLambdafromXiFourSigthreetight = new TH3F("f3fHistLambdafromXiFourSigthreetight","f3fHistLambdafromXiFourSigthreetight ", fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistLambdafromXiFourSigthreetight);
        
        //-
        f3fHistAntiLambdafromXiFourSigthreetight = new TH3F("f3fHistAntiLambdafromXiFourSigthreetight","f3fHistAntiLambdafromXiFourSigthreetight ",fNptBins, LambdaPtBins,CentbinNum, CentBins, xibinnumb, xibinlimits);
        fListHist->Add(f3fHistAntiLambdafromXiFourSigthreetight);
        
        //REC
        
        f3fHistCentInvMassVsPtLambdaRecFourSigthree = new TH3F("f3fHistCentInvMassVsPtLambdaRecFourSigthree","f3fHistCentInvMassVsPtLambdaRecFourSigthree",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecFourSigthree);
        
        f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree","f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree);
        
        f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntag = new TH3F("f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntag","f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntag",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntag);
        
        f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntagCut = new TH3F("f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntagCut","f3fHistCentInvMassVsPtLambdaRecFourUntagCuttight",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
        fListHist->Add(f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntagCut);
        
        /////
        f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntag = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntag","f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntag",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntag);
        
        f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntagCut = new TH3F("f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntagCut","f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntagCut",CentbinNum, CentBins, Massbinnumb,MassBins,fNptBins, LambdaPtBins);
        fListHist->Add(f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntagCut);
        
    }
    
    PostData(1,fListHist);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNetLambdaMCTrad::UserExec(Option_t *)
{
    
    const Int_t dim = fNptBins*2;
    
    Int_t ptChMC[dim];
    Int_t ptChRec[dim];
    Int_t ptChRecTag[dim];
    
    for(Int_t idx = 0; idx < dim; idx++)
    {
        ptChMC[idx] = 0;
        ptChRec[idx] = 0;
        ptChRecTag[idx] = 0;
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
            
            ////
            if(abseta < 0.1)
            {
                if(genpid == 3122)
                {
                    f2fHistGenCentVsPtLambdacopy->Fill(fCentrality,gpt);
                }
                if(genpid == -3122)
                {
                    f2fHistGenCentVsPtAntiLambdacopy->Fill(fCentrality,gpt);
                }
                
                if(genpid == -3312)
                {
                    f2fHistXiPluscopy->Fill(fCentrality,gpt);
                }
                if(genpid == 3312)
                {
                    f2fHistXiMinuscopy->Fill(fCentrality,gpt);
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
                Int_t iptbinRecTag = GetPtBin(mcpt);
                Int_t iptbinRec = GetPtBin(mcpt);
                if( iptbinRec < 0 || iptbinRec > fNptBins-1 ) continue;
                if( iptbinRecTag < 0 || iptbinRecTag > fNptBins-1 ) continue;
                
                
                //L
                if(TMath::Abs(eta) < 0.5)
                {
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex > 0.1 && TMath::Abs(posprnsg)  <= 3. && TMath::Abs(negpion)  <= 3.)
                    {
                        f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntag->Fill(fCentrality,invMassLambda,mcpt);
                        if(invMassLambda > 1.11 && invMassLambda < 1.122)
                        {
                            ptChRec[iptbinRec] += 1;
                            f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntagCut->Fill(fCentrality,invMassLambda,mcpt);
                        }
                        
                        if(fTreeVariablePID == 3122)
                        {
                            f3fHistCentInvMassVsPtLambdaRecFourSigthree->Fill(fCentrality,invMassLambda,mcpt);
                            if(isPrim){f2fHistRecPrimariesCentVsPtLambdaFourSigthree->Fill(fCentrality,mcpt); ptChRecTag[iptbinRecTag] += 1;}
                            if(isSecFromWeakDecay)
                            {
                                f2fHistRecSecCentVsPtLambdaFourSigthree->Fill(fCentrality,mcpt);
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == 3312) || (fTreeVariablePIDParent == 3322))
                                    {
                                        f3fHistLambdafromXiFourSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            else if(isSecFromMaterial) {f2fHistRecMatCentVsPtLambdaFourSigthree->Fill(fCentrality,mcpt);}
                            
                        }
                    }
                    
                    //L-BAR
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex > 0.25 && TMath::Abs(negprnsg)  <= 3 && TMath::Abs(pospion)  <= 3)
                    {
                        f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntag->Fill(fCentrality,invMassAntiLambda,mcpt);
                        if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.122)
                        {
                            ptChRec[iptbinRec+fNptBins] += 1;
                            f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntagCut->Fill(fCentrality,invMassAntiLambda,mcpt);
                        }
                        
                        if(fTreeVariablePID == -3122)
                        {
                            f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree->Fill(fCentrality,invMassAntiLambda,mcpt);
                            if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree->Fill(fCentrality,mcpt); ptChRecTag[iptbinRecTag+fNptBins] += 1;}
                            if(isSecFromWeakDecay)
                            {
                                f2fHistRecSecCentVsPtAntiLambdaFourSigthree->Fill(fCentrality,mcpt);
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == -3312) || (fTreeVariablePIDParent == -3322))
                                        
                                    {
                                        f3fHistAntiLambdafromXiFourSigthree->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            else if(isSecFromMaterial) {f2fHistRecMatCentVsPtAntiLambdaFourSigthree->Fill(fCentrality,mcpt);}
                        }
                    }
                }// |eta| < 0.8
                
                if(TMath::Abs(eta) < 0.1)
                {
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex > 0.1 && TMath::Abs(posprnsg)  <= 3. && TMath::Abs(negpion)  <= 3.)
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtLambdaFourSigthreetight->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == 3312) || (fTreeVariablePIDParent == 3322))
                                    {
                                        f3fHistLambdafromXiFourSigthreetight->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                    
                    //L-BAR
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex > 0.25 && TMath::Abs(negprnsg)  <= 3 && TMath::Abs(pospion)  <= 3)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreetight->Fill(fCentrality,mcpt); }
                            if(isSecFromWeakDecay)
                            {
                                if(isPrimParent)
                                {
                                    if ((fTreeVariablePIDParent == -3312) || (fTreeVariablePIDParent == -3322))
                                        
                                    {
                                        f3fHistAntiLambdafromXiFourSigthreetight->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                        }
                    }
                }// |eta| < 0.3
            } //MC condition
        }// zero onfly V0
    }// end of V0 loop
    
    
    //-------------------------------------------------
    Double_t ptContainerRec[dim+1];
    ptContainerRec[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerRec[i] = ptChRec[i-1];
    }
    fPtBinNplusNminusChRec->Fill(ptContainerRec);
    
    
    Double_t ptContainerRecTag[dim+1];
    ptContainerRecTag[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerRecTag[i] =  ptChRecTag[i-1];
    }
    fPtBinNplusNminusChRecTag->Fill(ptContainerRecTag);
    
    
    
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




