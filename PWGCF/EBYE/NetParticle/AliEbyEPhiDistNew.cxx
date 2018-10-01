/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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


//=========================================================================//
//                 AliEbyE Analysis for Net-Particle Study                 //
//                           Surya Prakash Pathak                          //
//                      surya.prakash.pathak@cern.ch                       //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                                                                         //
//                        (Last Modified 2018/08/27)                       //
//                 Dealing with Wide pT Window Modified to ESDs            //
//Some parts of the code are taken from J. Thaeder/ M. Weber NetParticle   //
//analysis task.                                                           //
//=========================================================================//

//ver: 2018/09/11

#include <Riostream.h>
#include "TList.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODHeader.h"
#include "AliAODpidUtil.h"
#include "AliEventCuts.h"
#include "AliMultSelection.h"
//#include "AliHelperPID.h"
#include "AliAnalysisTaskSE.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliEbyEPhiDistNew.h"

class AliEbyEPhiDistNew;

using namespace std;

ClassImp(AliEbyEPhiDistNew)

//------------------------------------------------------------------------
AliEbyEPhiDistNew::AliEbyEPhiDistNew():AliAnalysisTaskSE(),
 fThnList(NULL), fInputHandler(NULL), fMCEventHandler(NULL), fVevent(NULL), fArrayMC(NULL), fESDtrackCuts(NULL), fMCEvent(NULL), fMCStack(NULL), fEventCuts(NULL), fRun("LHC10h") , fCentralityEstimator("V0M"),

 fAODtrackCutBit(272), fVxMax(3.), fVyMax(3.), fVzMax(10.), fPtMin(0.15), fPtMax(1.5), fPhiMin(0.0), fPhiMax(6.28),fEtaMin(-1), fEtaMax(1), fNptBins(19), fNphiBins(6), fDcaXy(10.), fDcaZ(10.), fNcrossRows(80), fChi2NDF(4),

 fIsMC(kFALSE), fIsAOD(kFALSE), fIsRapCut(kFALSE),

 fNTracks(0), fCentrality(-1), fEventCounter(NULL), fHistCent(0x0),

fPIDResponse(0x0), fPIDCombined(0x0), fPidType(1), fMcPid(211), fPidStrategy(0), fNSigmaMaxITS(3.), fNSigmaMaxTPC(3.), fNSigmaMaxTOF(3.), fParticleSpecies(AliPID::kPion),

fPtBinNplusNminus(NULL), fPtBinNplusNminusTruth(NULL)
{
    for (Int_t i = 0; i < 2 ; i++){
        fHistCentRec[i] = NULL;
        fHistCentGen[i] = NULL;
        fCentPtEtaPhiThnGen[i] = NULL;
        fCentPtEtaPhiThnRec[i] = NULL;
        fCentPtEtaPhiThnRecPrim[i] = NULL;
        fCentPtEtaPhiThnSec[i] = NULL;
        fCentPtEtaPhiThnMat[i] = NULL;
        fCentPtEtaPhiThnMisId[i] = NULL;
        for (Int_t j = 0 ; j < 3; j++){
            fHistERec[i][j] = NULL;
            fHistERecPri[i][j] = NULL;
            fHistEGen[i][j] = NULL;
        }
    }
}

//---------------------------------------------------------------------------------------------

AliEbyEPhiDistNew::AliEbyEPhiDistNew(const char *name): AliAnalysisTaskSE(name),
fThnList(NULL), fInputHandler(NULL), fMCEventHandler(NULL), fVevent(NULL), fArrayMC(NULL), fESDtrackCuts(NULL), fMCEvent(NULL), fMCStack(NULL), fEventCuts(NULL), fRun("LHC10h") , fCentralityEstimator("V0M"),

fAODtrackCutBit(272), fVxMax(3.), fVyMax(3.), fVzMax(10.), fPtMin(0.15), fPtMax(1.5), fPhiMin(0.0), fPhiMax(6.28),fEtaMin(-1), fEtaMax(1), fNptBins(19), fNphiBins(6), fDcaXy(10.), fDcaZ(10.), fNcrossRows(80), fChi2NDF(4),

fIsMC(kFALSE), fIsAOD(kFALSE), fIsRapCut(kFALSE),

fNTracks(0), fCentrality(-1), fEventCounter(NULL), fHistCent(0x0),

fPIDResponse(0x0), fPIDCombined(0x0), fPidType(1), fMcPid(211), fPidStrategy(0), fNSigmaMaxITS(3.), fNSigmaMaxTPC(3.), fNSigmaMaxTOF(3.), fParticleSpecies(AliPID::kPion),

fPtBinNplusNminus(NULL), fPtBinNplusNminusTruth(NULL)
{
    for (Int_t i = 0; i < 2 ; i++){
        fHistCentRec[i] = NULL;
        fHistCentGen[i] = NULL;
        fCentPtEtaPhiThnGen[i] = NULL;
        fCentPtEtaPhiThnRec[i] = NULL;
        fCentPtEtaPhiThnRecPrim[i] = NULL;
        fCentPtEtaPhiThnSec[i] = NULL;
        fCentPtEtaPhiThnMat[i] = NULL;
        fCentPtEtaPhiThnMisId[i] = NULL;
        for (Int_t j = 0 ; j < 3; j++){
            fHistERec[i][j] = NULL;
            fHistERecPri[i][j] = NULL;
            fHistEGen[i][j] = NULL;
        }
    }
    DefineOutput(1, TList::Class()); // define the output of the analysis: in this case it's a list of histogram
}

//--------------------------------------------------------------------------------------------

AliEbyEPhiDistNew::~AliEbyEPhiDistNew(){
    //default destructor
    if (fThnList) delete fThnList;
    if (fEventCuts) delete fEventCuts;
    if (fESDtrackCuts) delete fESDtrackCuts;
}


//---------------------Create output objects--------------------------------------------------

void AliEbyEPhiDistNew::UserCreateOutputObjects(){
    
    fInputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!fInputHandler){
        AliError("No PID response task found !!");
    }
    fPIDResponse = dynamic_cast<AliInputEventHandler *>(fInputHandler)->GetPIDResponse();
    if (!fPIDResponse){
        AliError("No PID response task found !!");
    }
    if (fRun == "LHC15o"){
        fEventCuts = new AliEventCuts();
    }
    
    fThnList = new TList;
    fThnList->SetOwner(kTRUE);
    
    if (!fIsAOD) {
        if (!fESDtrackCuts) {
            fESDtrackCuts = new AliESDtrackCuts;
            fESDtrackCuts->SetName(Form("NetPesdCut_%d",fPidType));
            fESDtrackCuts->SetEtaRange(fEtaMin,fEtaMax);
            fESDtrackCuts->SetPtRange(0.1, 1e10);
            //fESDtrackCuts->SetPhiRange(fPhiMin,fPhiMax); // I added this
            
            //TPC track
            fESDtrackCuts->SetMinNCrossedRowsTPC(fNcrossRows); // default 80, sys 60 and 100
            fESDtrackCuts->SetMaxChi2PerClusterTPC(fChi2NDF); // 4 default,
            fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
            fESDtrackCuts->SetRequireTPCRefit(kTRUE);
            fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
            
            if (fPidType == 0){
                fESDtrackCuts->SetMaxDCAToVertexXY(fDcaXy); // default 2 ,
            } else {
                if (fDcaXy == 2) fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.035/pt^1.01"); // default
            }
            
            fESDtrackCuts->SetRequireITSRefit(kTRUE);
            fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
            fESDtrackCuts->SetMaxChi2PerClusterITS(36);
            fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
            fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);
            
            //Vertex related
            fESDtrackCuts->SetMaxDCAToVertexZ(fDcaZ);
            fESDtrackCuts->SetDCAToVertex2D(kFALSE);
            fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
        }
        else cout << ">>>> User Track Cuts <<<<" << endl;
    }
    CreatePhiHist();
    
    fHistCent = new TH1F("fHistCentPID", "Centrality", 100, -0.5, 99.5);
    fThnList->Add(fHistCent);
    
    fEventCounter =  new TH1D("fEventCounter","EventCounter", 20, -0.5, 19.5);
    fThnList->Add(fEventCounter);
    
    PostData(1,fThnList);

}

//--------------------------------------------------------------------------------------

void AliEbyEPhiDistNew::CreatePhiHist() {
    const Char_t *fgkHistName[4] = {"Nch","Npi","Npi", "Nka"};
    const Char_t *fgkHistLat[2][4] = {{"N^{-}","#pi^{-}","K^{-}","P^{-}"},{"N^{+}","#pi^{+}","K^{+}","P^{+}"}};
    
    const Char_t *fgkHistCharge[2] = {"Minus","Plus"};
    
    const Int_t cbin = 81;
    Double_t CentBins[cbin+1] ;
    for (Int_t ic = 0; ic <=cbin ; ic++) CentBins[ic] = ic - 0.5;
    
    const Int_t ebin = 8;
    Double_t EtaBins[ebin+1];
    for (Int_t ie = 0; ie <= ebin; ie++) EtaBins[ie] = ie - 0.5;
    
    const Int_t ptBins = 19;
    Double_t pidPtBins[ptBins+1] = {0.35, 0.4, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.55};
    
    const Int_t phiBins = 6;
    Double_t pidPhiBins[phiBins+1] = {0.0, 1.04, 2.09, 3.14, 4.18, 5.23, 6.28};
    
    
    const Char_t *gstName[3] = {"pT", "Eta", "Phi"};
    const Char_t *gstLat[3] = {"p_{T}","#eta", "#phi"};
    const Char_t *PidCut[3] = {"TPC", "TPC+TOF", "ITS+TPC+TOF"};
    
    for (Int_t k = 0; k<2; k++) {
        
        Int_t pidtype = fPidType;
        
        fHistCentRec[k] = new TH2F(Form("fHist%s%sRecCent",fgkHistCharge[k], fgkHistName[pidtype]),Form(" Rec%s",fgkHistName[pidtype]), cbin, CentBins, phiBins, pidPhiBins);
        fHistCentRec[k]->Sumw2();
        fThnList->Add(fHistCentRec[k]);
        
        //Rec
        fCentPtEtaPhiThnRec[k] = new TH3F(Form("fCentPtPhiThnRec%s%s",fgkHistCharge[k],fgkHistName[pidtype]),"Rec:cent-pT-Phi",cbin,CentBins,ptBins,pidPtBins,phiBins,pidPhiBins);
        fCentPtEtaPhiThnRec[k]->GetXaxis()->SetTitle("Centrality");
        fCentPtEtaPhiThnRec[k]->GetYaxis()->SetTitle("p_{T}");
        fCentPtEtaPhiThnRec[k]->GetZaxis()->SetTitle("#phi");
        fCentPtEtaPhiThnRec[k]->Sumw2();
        fThnList->Add(fCentPtEtaPhiThnRec[k]);
        
        if (fIsMC){
            
            fHistCentRec[k] = new TH2F(Form("fHist%s%sRecCent",fgkHistCharge[k], fgkHistName[pidtype]),Form(" Rec%s",fgkHistName[pidtype]), cbin, CentBins, phiBins, pidPhiBins);
            fHistCentRec[k]->Sumw2();
            fThnList->Add(fHistCentRec[k]);
            
            fHistCentGen[k] = new TH2F(Form("fHist%s%sGenCent",fgkHistCharge[k], fgkHistName[pidtype]),Form(" Gen%s",fgkHistName[pidtype]), cbin, CentBins, phiBins, pidPhiBins);
            fHistCentGen[k]->Sumw2();
            fThnList->Add(fHistCentGen[k]);
            
            
            //Gen
            fCentPtEtaPhiThnGen[k]=  new TH3F(Form("fCentPtEtsThnGen%s%s",fgkHistCharge[k],fgkHistName[pidtype]),"Gen:cent-pT-Phi",cbin,CentBins,ptBins,pidPtBins,phiBins,pidPhiBins);
            fCentPtEtaPhiThnGen[k]->GetXaxis()->SetTitle("Centrality");
            fCentPtEtaPhiThnGen[k]->GetYaxis()->SetTitle("p_{T}");
            fCentPtEtaPhiThnGen[k]->GetZaxis()->SetTitle("#phi");
            fCentPtEtaPhiThnGen[k]->Sumw2();
            fThnList->Add(fCentPtEtaPhiThnGen[k]);
            
            //Rec--
            fCentPtEtaPhiThnRec[k] = new TH3F(Form("fCentPtPhiThnRec%s%s",fgkHistCharge[k],fgkHistName[pidtype]),"Rec:cent-pT-Phi",cbin,CentBins,ptBins,pidPtBins,phiBins,pidPhiBins);
            fCentPtEtaPhiThnRec[k]->GetXaxis()->SetTitle("Centrality");
            fCentPtEtaPhiThnRec[k]->GetYaxis()->SetTitle("p_{T}");
            fCentPtEtaPhiThnRec[k]->GetZaxis()->SetTitle("#phi");
            fCentPtEtaPhiThnRec[k]->Sumw2();
            fThnList->Add(fCentPtEtaPhiThnRec[k]);
            
            //Rec Primary
            fCentPtEtaPhiThnRecPrim[k]=new TH3F (Form("fCentPtEtaThnRecPrim%s%s",fgkHistCharge[k],fgkHistName[pidtype]),"Rec-Prim-cent-pT-Phi",cbin,CentBins,ptBins,pidPtBins,phiBins,pidPhiBins);
            fCentPtEtaPhiThnRecPrim[k]->GetXaxis()->SetTitle("Centrality");
            fCentPtEtaPhiThnRecPrim[k]->GetYaxis()->SetTitle("p_{T}");
            fCentPtEtaPhiThnRecPrim[k]->GetZaxis()->SetTitle("#phi");
            fCentPtEtaPhiThnRecPrim[k]->Sumw2();
            fThnList->Add(fCentPtEtaPhiThnRecPrim[k]);
            
            //Rec. Secondary
            fCentPtEtaPhiThnSec[k]=new TH3F (Form("fCentPtEtaThnSec%s%s",fgkHistCharge[k],fgkHistName[pidtype]),"Rec-Sec-cent-pT-Phi",cbin,CentBins,ptBins,pidPtBins,phiBins,pidPhiBins);
            fCentPtEtaPhiThnSec[k]->GetXaxis()->SetTitle("Centrality");
            fCentPtEtaPhiThnSec[k]->GetYaxis()->SetTitle("p_{T}");
            fCentPtEtaPhiThnSec[k]->GetZaxis()->SetTitle("#phi");
            fCentPtEtaPhiThnSec[k]->Sumw2();
            fThnList->Add(fCentPtEtaPhiThnSec[k]);
            
            //Rec. Material
            fCentPtEtaPhiThnMat[k]=new TH3F (Form("fCentPtEtaThnMat%s%s",fgkHistCharge[k],fgkHistName[pidtype]),"Rec-Mat-cent-pT-Phi",cbin,CentBins,ptBins,pidPtBins,phiBins,pidPhiBins);
            fCentPtEtaPhiThnMat[k]->GetXaxis()->SetTitle("Centrality");
            fCentPtEtaPhiThnMat[k]->GetYaxis()->SetTitle("p_{T}");
            fCentPtEtaPhiThnMat[k]->GetZaxis()->SetTitle("#phi");
            fCentPtEtaPhiThnMat[k]->Sumw2();
            fThnList->Add(fCentPtEtaPhiThnMat[k]);
            
            //Rec. Mis ID
            fCentPtEtaPhiThnMisId[k]=new TH3F (Form("fCentPtEtaThnMisId%s%s",fgkHistCharge[k],fgkHistName[pidtype]),"Rec-MisId-cent-pT-Phi",cbin,CentBins,ptBins,pidPtBins,phiBins,pidPhiBins);
            fCentPtEtaPhiThnMisId[k]->GetXaxis()->SetTitle("Centrality");
            fCentPtEtaPhiThnMisId[k]->GetYaxis()->SetTitle("p_{T}");
            fCentPtEtaPhiThnMisId[k]->GetZaxis()->SetTitle("#phi");
            fCentPtEtaPhiThnMisId[k]->Sumw2();
            fThnList->Add(fCentPtEtaPhiThnMisId[k]);
            
        }
    }
    
    //-----For Thn Sparse----------------------
    const Int_t dim = 229; // 1 centrality bin + ( 19 pt bins )*2 *6
        Int_t bin[dim];
        bin[0] = 81;
        for (Int_t ibin = 1; ibin<dim ; ibin++) bin[ibin] = 500;
        
        Double_t min[dim];
        for (Int_t jbin = 0; jbin<dim ; jbin++) min[jbin] = -0.5;
        
        Double_t max[dim];
        max[0] = 80.5 ;
        for(Int_t jbin = 1; jbin<dim; jbin++) max[jbin] = 499.5;
        
        fPtBinNplusNminus = new THnSparseI("fPtBinNplusNminus","cent-nplus-nminus",dim, bin, min, max);
        fThnList->Add(fPtBinNplusNminus);
        
        if(fIsMC){
            fPtBinNplusNminusTruth = new THnSparseI("fPtBinNplusNminusTruth","cent-nplus-nminus",dim, bin, min, max);
            fThnList->Add(fPtBinNplusNminusTruth);
        }
}

    //------------------------------------------------------
    
    void AliEbyEPhiDistNew::LocalPost(){
        PostData(1,fThnList);
    }

    void AliEbyEPhiDistNew::UserExec (Option_t *){
        
        if (fPidType < 1 || fPidType >3){
            AliError("PID are not supported");
            return;
        }
        
        const Int_t dim = 38; // number of pT bins * 2
        const Int_t kPhi = 6; // number of Phi Bins
        Int_t pTPhi[dim][kPhi];
        Int_t pTPhiMC[dim][kPhi];
        
        for (Int_t idx= 0; idx < dim ; idx++){
            for (Int_t it = 0; it <kPhi ; it++){
                pTPhi[idx][it] = 0.;
                pTPhiMC[idx][it] = 0.;
            }
        }
        
        fEventCounter->Fill(1);
        if (!fInputHandler)
            fInputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
        if (!fInputHandler){
            AliError("No InputHandler");
            return;
        }
        
        fVevent = dynamic_cast<AliVEvent*>(fInputHandler->GetEvent());
        if (!fVevent){
            cout << "Error : fVEvent not available\n" << endl;
            LocalPost(); return;
            
        }
        
        //--------------Initiate MC-----------
        if (fIsMC){
            fMCEvent = NULL;
            fEventCounter->Fill(8);
            if(fIsAOD){
                fArrayMC = NULL;
                fArrayMC = dynamic_cast<TClonesArray*>(fVevent->FindListObject(AliAODMCParticle::StdBranchName()));
                if(!fArrayMC)
                    AliFatal("No array of MC particles found !!");
            }
            else{
                fMCEventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
                if(!fMCEventHandler){
                    AliError("No MC event Handler available");
                    LocalPost(); return;
                }
                fMCEvent = fMCEventHandler->MCEvent();
                if (!fMCEvent){
                    cout << "Error: Could not retrieve MC event" << endl;
                    LocalPost();return;
                }
                fMCStack = fMCEvent->Stack();
                if(!fMCStack){
                    cout << "ERROR: Could not retrive MC stack" << endl;
                    LocalPost(); return;
                }
            }
        } // if -------MC--
        
        
        //Plie up cout for Run2
        
        if(fRun == "LHC15o"){
            if(!fEventCuts->AcceptEvent(fVevent)){
                LocalPost(); return;
            }
        }
        
        const AliVVertex *vertex = fVevent->GetPrimaryVertex();
        if(!vertex) { LocalPost(); return;}
        
        Bool_t vtest = kFALSE;
        Double32_t fCov[6];
        vertex->GetCovarianceMatrix(fCov);
        if(vertex->GetNContributors() > 0) {
            if(fCov[5] != 0){
                vtest = kTRUE;
            }
        }
        if(!vtest) {LocalPost(); return;}
        
        if(TMath::Abs(vertex->GetX()) > fVxMax) {LocalPost(); return; }
        if(TMath::Abs(vertex->GetY()) > fVyMax) {LocalPost(); return; }
        if(TMath::Abs(vertex->GetZ()) > fVzMax) {LocalPost(); return; }
        
        //--------Centrality task----------------------------------------
        
        if( fRun == "LHC10h" || fRun == "LHC11h") {
            AliCentrality *centrality = fVevent -> GetCentrality();
            if(!centrality) return;
            if (centrality->GetQuality() != 0) {LocalPost(); return; }
            
            fCentrality = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
        }
        else if (fRun == "LHC15o") {
            AliMultSelection *fMultSelection = (AliMultSelection *) fVevent->FindListObject("MultSelection");
            if(!fMultSelection){
                cout << "AliMultSelection object not found! " << endl;
                return;
            }
            else fCentrality = fMultSelection->GetMultiplicityPercentile(fCentralityEstimator.Data(),false);
        }
        else if (fRun == "LHC10hAMPT"){
            //centrality range from impact parameter
            AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
            if(!genHeader) return;
            
            Double_t impactPar = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
            
            if( impactPar >= 0. && impactPar < 3.51 ) fCentrality = 1.; //0-5
            if( impactPar >= 3.51 && impactPar < 4.96 ) fCentrality = 6.; //5-10
            if( impactPar >= 4.96 && impactPar < 6.08 ) fCentrality = 11.;//10-15
            if( impactPar >= 6.08 && impactPar < 7.01 ) fCentrality = 16.;//15-20
            if( impactPar >= 7.01 && impactPar < 7.84 ) fCentrality = 21.;//20-25
            if( impactPar >= 7.84 && impactPar < 8.59 ) fCentrality = 26.;//25-30
            if( impactPar >= 8.59 && impactPar < 9.27 ) fCentrality = 31.;//30-35
            if( impactPar >= 9.27 && impactPar < 9.92 ) fCentrality = 36.;//35-40
            if( impactPar >= 9.92 && impactPar < 10.5 ) fCentrality = 41.;//40-45
            if( impactPar >= 10.5 && impactPar < 11.1 ) fCentrality = 46.;//45-50
            if( impactPar >= 11.1 && impactPar < 11.6 ) fCentrality = 51.;//50-55
            if( impactPar >= 11.6 && impactPar < 12.1 ) fCentrality = 56.;//55-60
            if( impactPar >= 12.1 && impactPar < 12.6 ) fCentrality = 61.;//60-65
            if( impactPar >= 12.6 && impactPar < 13.1 ) fCentrality = 66.;//65-70
            if( impactPar >= 13.1 && impactPar < 13.6 ) fCentrality = 71.;//70-75
            if( impactPar >= 13.6 && impactPar < 14.0 ) fCentrality = 76.;//75-80
            
        }
        else {
            cout << "Wrong run period for centrality" << endl;
            LocalPost(); return;
        }
        
        if (fCentrality < 0 || fCentrality >=80) return;
        fHistCent->Fill(fCentrality);
        fEventCounter->Fill(2);
        
        fEventCounter->Fill(3);
        //-------tracks--------------
        fNTracks = fVevent->GetNumberOfTracks();
        
        Double_t nRec[2]= {0.,0.};
        Double_t nGen[2]= {0.,0.};
        
        //-----------------track loop---------------
        for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack){
            AliVTrack *track = static_cast<AliVTrack*>(fVevent->GetTrack(idxTrack));
            if(!AcceptTrackL(track))continue;
            
            Int_t icharge = track->Charge() < 0 ? 0 : 1 ;
            Float_t lPt = (Float_t)track->Pt();
            Float_t lPz = track->Pz();
            Float_t lP = 0.;
            if(fTotP) lP = track->GetInnerParam()->GetP();
            Float_t lEta = (Float_t)track->Eta();
            Float_t lPhi = (Float_t)track->Phi();
            
            //Get the pT (p) bin
            Int_t iptbin = -1;
            if(fTotP) iptbin = GetPtBin(lP);
            else iptbin = GetPtBin(lPt);
            if(iptbin < 0 || iptbin > fNptBins - 1) continue;
            
            //Get the Phi bins
            Int_t iphibin = -1 ;
            iphibin = GetPhiBin(lPhi);
            if (iphibin < 0 || iphibin > fNphiBins-1) continue;
            
            //Get the Eta Bin--
            Int_t etabin = -1;
            etabin = GetEtaBin(TMath::Abs(lEta));
            if ( etabin < 0 || etabin >7) continue;
            
            //--Fill Histos-------------------
            
            Float_t RecContainer[3];
            RecContainer[0] = fCentrality;
            if(fTotP) RecContainer[1] = lP;
            else RecContainer[1] = lPt;
            RecContainer[2] = lPhi;
            
            Bool_t isPid = kFALSE;
            
            if (fPidType != 0){
                isPid = IsPidPassed(track);
                
                if(isPid){
                    if(fIsMC){
                        cout << "yes it is inside " << endl;
                        fCentPtEtaPhiThnRec[icharge]->Fill(RecContainer[0],RecContainer[1], RecContainer[2]);
                    }
                    else {
                        fHistCentRec[icharge]->Fill(RecContainer[0],RecContainer[2]);
                        fCentPtEtaPhiThnRec[icharge]->Fill(RecContainer[0],RecContainer[1], RecContainer[2]);
                    }
                    nRec[icharge] += 1;
                    
                    if (icharge == 1){
                        pTPhi[iptbin][iphibin] += 1;
                        
                    }
                    if(icharge == 0){
                        pTPhi[iptbin+fNphiBins][iphibin] += 1;
                    }
                }
            }
            
            //MC loop for Physical primary
            
            if (fIsMC){
                Int_t label = TMath::Abs(track->GetLabel());
                
                Bool_t isPhysicalPrimary = 0;
                Bool_t isSecondaryFromWeakDecay = 0;
                Bool_t isSecondaryFromMaterial = 0;
                AliVParticle* particle = NULL;
                
                if(track->InheritsFrom("AliESDtrack")) {
                    particle = static_cast<AliVParticle*>(fMCEvent->GetTrack(label));
                    if(!particle) return;
                    isPhysicalPrimary = fMCStack->IsPhysicalPrimary(label);
                    isSecondaryFromWeakDecay = fMCStack->IsSecondaryFromWeakDecay(label);
                    isSecondaryFromMaterial = fMCStack->IsSecondaryFromMaterial(label);
                }
                else {
                    particle = static_cast<AliVParticle*>(fArrayMC->At(label));
                    isPhysicalPrimary = (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();
                    isSecondaryFromWeakDecay = (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromWeakDecay();
                    isSecondaryFromMaterial = (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromMaterial();
                }
                Float_t fpTRec = particle->Pt();
                Float_t fpZRec = particle->Pz();
                Float_t fPRec = particle ->P();
                Float_t fEtaRec = particle->Eta();
                Float_t fPhiRec = particle->Phi();
                Int_t pdg = TMath::Abs(particle->PdgCode());
                
                //Get the eta (Rap) bin---
                
                Int_t etabinRecMC = -1;
                etabinRecMC = GetEtaBin(TMath::Abs(fEtaRec));
                if(etabinRecMC < 0 || etabinRecMC > 7) continue;
                
                Double_t RecMCContainer[3];
                RecMCContainer[0] = fCentrality;
                if(fTotP) RecMCContainer[1] = fPRec;
                else RecMCContainer[1] = fpTRec;
                RecMCContainer[2] = fPhiRec;
                
                //PID--------
                if (isPid) {
                    if(isPhysicalPrimary){
                        if (pdg == fMcPid){
                        fCentPtEtaPhiThnRecPrim[icharge]->Fill (RecMCContainer[0], RecMCContainer[1], RecMCContainer[2]);
                    }
                    else {
                        fCentPtEtaPhiThnMisId[icharge]->Fill(RecMCContainer[0], RecMCContainer[1], RecMCContainer[2]);
                    }
                }
                    else if (isSecondaryFromWeakDecay) {
                        fCentPtEtaPhiThnSec[icharge] -> Fill(RecMCContainer[0], RecMCContainer[1], RecMCContainer[2]);
                    }
                    else if (isSecondaryFromMaterial){
                        fCentPtEtaPhiThnMat[icharge] -> Fill(RecMCContainer[0], RecMCContainer[1], RecMCContainer[2]);
                    }
            }
        
        }
        }// Reconstructed track loop
            
            if (fIsMC){
                fHistCentRec[0] -> Fill(fCentrality, nRec[0]);
                fHistCentRec[1] -> Fill(fCentrality, nRec[1]);
            }
            
            const Int_t thndim = dim * kPhi;
            Double_t ptContainer[thndim+1];
            
            ptContainer[0] = (Double_t)fCentrality;
            
            for(Int_t ipt = 0; ipt <  dim ; ipt++) {
                for (Int_t jphi = 0 ; jphi < kPhi ; jphi++){
                    Int_t k = (ipt * kPhi) + jphi;
                    ptContainer[k+1] = pTPhi[ipt][jphi];
                    if(pTPhi[ipt][jphi] > 4) fEventCounter->Fill(17);
                }
            }
            fPtBinNplusNminus -> Fill(ptContainer);
            
            fEventCounter -> Fill(7);
            //---------------------------------------------------------------------
            
            if (fIsMC){
                fEventCounter -> Fill(8);
                if (fIsAOD){
                    for(Int_t idxMC = 0; idxMC < fArrayMC-> GetEntries(); idxMC++) {
                        AliAODMCParticle *particle = static_cast<AliAODMCParticle*> (fArrayMC->At(idxMC));
                        if(!particle) continue;
                        
                        if (!particle->IsPhysicalPrimary()) continue;
                        if (!AcceptTrackLMC((AliVParticle*)particle)) continue;
                        Int_t icharge = (particle->PdgCode() < 0) ? 0 :1;
                        
                        Float_t fpTGen = particle->Pt();
                        Float_t fpzGen = particle->Pz();
                        Float_t fpGen  = particle->P();
                        Float_t fEtaGen = particle-> Eta();
                        Float_t fPhiGen = particle-> Phi();
                        
                        Int_t pdg = TMath::Abs(particle->PdgCode());
                        if (pdg != fMcPid) continue;
                        
                        Int_t iptbinMC = -1;
                        if (fTotP) iptbinMC = GetPtBin(fpGen); // total p bin
                        else iptbinMC = GetPtBin(fpTGen);
                        if (iptbinMC < 0 || iptbinMC > fNptBins -1) continue;
                        
                        Int_t iphibinMC = -1 ;
                        iphibinMC = GetPhiBin(fPhiGen);
                        if (iphibinMC < 0 || iphibinMC > fNphiBins-1) continue;
                        
                        Int_t etabinMC = -1;
                        etabinMC = GetEtaBin(TMath :: Abs(fEtaGen));
                        if (etabinMC < 0 || etabinMC > 7) continue;
                        
                        Double_t GenContainer[3];
                        GenContainer[0] = fCentrality;
                        if(fTotP) GenContainer[1] = fpGen;
                        else GenContainer[1] = fpTGen;
                        GenContainer[2] = fPhiGen;
                    
                        
                        fCentPtEtaPhiThnGen[icharge]->Fill(GenContainer[0], GenContainer[1], GenContainer[2]);
                        
                        nGen[icharge] += 1;
                        if (icharge == 1){
                            pTPhiMC[iptbinMC][iphibinMC] += 1;
                        }
                        if (icharge == 0){
                            pTPhiMC[iphibinMC][iphibinMC] += 1;
                        }
                    }
                    fEventCounter -> Fill(9);
                }
                else {
                    fEventCounter->Fill(10);
                    for(Int_t idxMC = 0; idxMC < fMCStack->GetNprimary(); ++idxMC) {
                        AliVParticle* particle = fMCEvent->GetTrack(idxMC);
                        if(!particle)
                            continue;
                        if (!fMCStack->IsPhysicalPrimary(idxMC)) continue;
                        
                        if (!AcceptTrackLMC(particle)) continue;
                        Int_t icharge = (particle->PdgCode() < 0) ? 0 : 1;
                        
                        Float_t fpTGen   = particle->Pt();
                        Float_t fpzGen   = particle->Pz();
                        Float_t fpGen    = particle->P();
                        Float_t fEtaGen  = particle->Eta();
                        Float_t fPhiGen  = particle->Phi();
                        
                        Int_t pdg = TMath::Abs(particle->PdgCode());
                        if (pdg != fMcPid) continue;
                        
                        //pt (p) bin-------
                        
                        Int_t iptbinMC = -1;
                        if (fTotP) iptbinMC = GetPtBin(fpGen); // p bin
                        else iptbinMC = GetPtBin(fpTGen); // pt bin
                        if (iptbinMC < 0 || iptbinMC > fNptBins - 1)continue;
                        
                        //phibins
                        Int_t iphibinMC = -1 ;
                        iphibinMC = GetPhiBin(fPhiGen);
                        if (iphibinMC < 0 || iphibinMC > fNphiBins-1) continue;
                        
                        //Eta Bin
                        Int_t etabinMC = -1;
                        etabinMC = GetEtaBin(TMath :: Abs(fEtaGen));
                        if (etabinMC < 0 || etabinMC > 7) continue;
                        
                        Double_t GenContainer[3];
                        GenContainer[0] = fCentrality;
                        if (fTotP) GenContainer[1] = fpGen;
                        else GenContainer[1] = fpTGen;
                        GenContainer[2] = fPhiGen;
                        
                        fCentPtEtaPhiThnGen[icharge] -> Fill(GenContainer[0], GenContainer[1], GenContainer[2]);
                        
                        nGen[icharge] += 1;
                        
                        //-------------------------
                        
                        if(icharge == 1){
                            pTPhiMC [iptbinMC][iphibinMC] += 1;
                        }
                        if (icharge == 0) {
                            pTPhiMC [iptbinMC+fNptBins][iphibinMC] += 1;
                        }
                    }
                    
                    fEventCounter->Fill(11) ;
                } // else is for ESD
                
                fHistCentGen[0] -> Fill(fCentrality, nGen[0]);
                fHistCentGen[1] -> Fill(fCentrality, nGen[1]);
                
                Double_t ptContainerMC[thndim+1];
                ptContainerMC[0] = (Double_t) fCentrality;
                
                for (Int_t ipt = 0; ipt < dim ; ipt++){
                    for (Int_t jphi = 0 ; jphi < kPhi ; jphi++){
                        Int_t k = (ipt*kPhi) + jphi;
                        ptContainerMC[k+1] = pTPhiMC[ipt][jphi];
                        if (pTPhiMC[ipt][jphi] > 4) fEventCounter-> Fill(18);
                    }
                }
                fPtBinNplusNminusTruth->Fill(ptContainerMC);
            }
            
            fEventCounter-> Fill(12);
            PostData(1,fThnList);
        }

//--------------------------------------

Bool_t AliEbyEPhiDistNew::AcceptTrackL(AliVTrack *track) const {
    
    if (!track) return kFALSE;
    if (track->Charge() == 0) return kFALSE;
    
    if (fIsAOD) {
        AliAODTrack * trackAOD = dynamic_cast<AliAODTrack*>(track);
        if(!track){
            AliError("Pointer to dynamic_cast<AliAODTrack*>(track) = ZERO");
            return kFALSE;
        }
        if(!trackAOD->TestFilterBit(fAODtrackCutBit))
            return kFALSE;
    } else { // ESD
        if (!fESDtrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track))) return kFALSE;
    }
    
    Double_t pt = track->Pt();
    Double_t pz = track->Pz();
    
    if (fTotP){
        if(!track->GetInnerParam() ) return kFALSE;
        Double_t ptotTPC = track->GetInnerParam()->GetP(); // total momentum
        if (ptotTPC < fPtMin || ptotTPC >fPtMax) return kFALSE; // cut on momentum
    }
    else {
        if (pt < fPtMin || pt > fPtMax) return kFALSE; // pT cut
    }
    
    if (fIsRapCut) {
        Double_t rap = GetRapidity(pt, pz);
        if (TMath::Abs(rap) > 0.5) return kFALSE; // rapidity cut
        if (TMath::Abs(track->Eta()) > fEtaMax) return kFALSE;
    }
    else {
        if (TMath::Abs(track->Eta()) > fEtaMax) return kFALSE;
    }
    
    return kTRUE;
    
}
//---------------------------------------------------------------------------
Bool_t AliEbyEPhiDistNew::AcceptTrackLMC(AliVParticle *particle) const {
    
    if (!particle) return kFALSE;
    if (particle->Charge() == 0.0) return kFALSE;
    
    Double_t ptotMC = particle->P();
    Double_t ptMC = particle->Pt();
    Double_t pzMC = particle->Pz();
    
    if (fTotP){
        if (ptotMC < fPtMin || ptotMC >fPtMax) return kFALSE; // cut on momentum
    }
    else {
        if (ptMC < fPtMin || ptMC > fPtMax) return kFALSE;
    }
    
    //rapidity cut
    if (fIsRapCut) {
        Double_t rapMC = GetRapidity(ptMC,pzMC);
        if (TMath::Abs(rapMC) > 0.5) return kFALSE; // rapidity cut
        if (TMath::Abs(particle->Eta()) > fEtaMax) return kFALSE;
    }
    else {
        if (TMath::Abs(particle->Eta()) > fEtaMax ) return kFALSE;
    }
    return kTRUE;
    
}

//----------------------------------------------------------------

Double_t AliEbyEPhiDistNew::GetRapidity(Float_t pt, Float_t pz) const{
    Double_t partMass = AliPID::ParticleMass(fParticleSpecies);
    Double_t en = TMath::Sqrt(pt*pt + pz*pz + partMass * partMass);
    Double_t rap = -999.;
    if (en != TMath::Abs(pz)) {
        rap = 0.5 * TMath::Log( (en + pz)/(en - pz) );
    }
    else rap = -999.;
    
    return rap;
}

//--------------------------------------------------------------------

Int_t AliEbyEPhiDistNew::GetPtBin(Double_t pt){
    Int_t bin = -1;
    Double_t pidPtBins[20] = {0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.55};
    for (Int_t iBin = 0 ; iBin < fNptBins; iBin++){
        
        if (iBin == fNptBins-1){
            if (pt >= pidPtBins[iBin] && pt <= pidPtBins[iBin + 1]){
                bin = iBin;
                break;
            }
        }
        else {
            if (pt >= pidPtBins[iBin] && pt < pidPtBins[iBin + 1]){
                bin = iBin;
                break;
            }
        }
    }
    return bin;
}

//---------------------Eta------------------------
Int_t AliEbyEPhiDistNew::GetEtaBin(Float_t eta){
    
    Int_t etabin = -1;
    Float_t EtaRange[9] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
    
    for(Int_t iBin = 0; iBin < 8 ; iBin++){
        if (iBin == 7){
            
            if(eta >= EtaRange[iBin] && eta <= EtaRange[iBin + 1]){
                etabin = iBin;
                break;
            }
        }
        else {
            if (eta >= EtaRange[iBin] && eta < EtaRange[iBin + 1]){
                etabin = iBin;
                break;
            }
        }
    }
    return etabin;
}

//-----------------------Phi----------------------------------------

Int_t AliEbyEPhiDistNew::GetPhiBin(Double_t Phi){
    
    Int_t phibin = -1;
    
    Double_t pidPhiBins[7] = {0.0, 1.04, 2.09, 3.14, 4.18, 5.23, 6.28};
    
    for (Int_t pBin = 0; pBin < fNphiBins; pBin++){
        
        if (pBin == fNphiBins - 1){
            if (Phi >= pidPhiBins[pBin] && Phi <= pidPhiBins[pBin+1]){
                phibin = pBin ;
                break;
            }
        }
        else {
            if( Phi >= pidPhiBins[pBin] && Phi < pidPhiBins[pBin +1]){
                phibin = pBin;
                break;
            }
        }
        
    }
    return phibin;
}

//------------------------------------------------------------------

void AliEbyEPhiDistNew::Terminate(Option_t *) {
    cout << "-------------------------------------\n"
            "         Terminating the task        \n"
            "-------------------------------------\n";
}


void AliEbyEPhiDistNew::SetPidType(Int_t i) {
    fPidType = i, fMcPid = -9999;
    
    if (fPidType == 1) {fParticleSpecies = AliPID::kPion; fMcPid = 211; }
    if (fPidType == 2) {fParticleSpecies = AliPID::kKaon; fMcPid = 321; }
    if (fPidType == 3) {fParticleSpecies = AliPID::kProton; fMcPid = 2212;}
}

Double_t AliEbyEPhiDistNew::TOFBetaCalc(AliVTrack *track) const {
    //TOF beta calculation
    
    Double_t tofTime = track->GetTOFsignal();
    
    Double_t c = TMath::C()*1.E-9; // m/ns
    Float_t startTime = fPIDResponse->GetTOFResponse().GetStartTime(((AliVTrack*)track)->P()); // in ps
    Double_t length = fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron)*1E-3*c ;
    tofTime -= startTime; // substract starttime to the signal
    Double_t tof = tofTime*1E-3 ; // ns, average T0 fill substracted, no info from T0 detector
    tof = tof * c;
    return length/tof;
}

//-------------------------------------------------------------------------

Double_t AliEbyEPhiDistNew::GetMass(AliPID::EParticleType id) const{
    //return mass according to AliHelperParticleSpecies_t. If undefined return -999.
    
    Double_t mass = -999;
    
    if (id == AliPID::kProton)  {mass = 9.38271999999999995e-01;}
    if (id == AliPID::kKaon)   { mass=4.93676999999999977e-01; }
    if (id == AliPID::kPion)   { mass=1.39570000000000000e-01; }
    return mass;
}



//------------------------------------------------------------------------------

Bool_t AliEbyEPhiDistNew::IsPidPassed(AliVTrack * track) {
    //return kTRUE ; PID strategy is from Jochen, Hans and deepika
    
    Bool_t isAcceptedITS = kFALSE;
    Bool_t isAcceptedTPC = kFALSE;
    Bool_t isAcceptedTOF = kFALSE;
    Bool_t isAccepted    = kFALSE;
    Bool_t hasPIDTOF     = kFALSE;
    
    Double_t *pid = new Double_t[3];
    pid[0] = 10.;
    pid[1] = 10.;
    pid[2] = 10.;
    
    Double_t pt = track->Pt();
    Double_t ptot = track->GetInnerParam()->GetP(); // total momentum
    
    //---------------------------| el, mu,  pi,  k,    p   | Pt cut offs from spectra
    //ITS--------------
    Double_t ptLowITS[5]       = { 0., 0., 0.2,  0.2,  0.3  };
    Double_t ptHighITS[5]      = { 0., 0., 0.6,  0.6,  1.1  };
    //TPC---------------
    Double_t ptLowTPC[5]       = { 0., 0., 0.2,  0.325, 0.3  };
    Double_t ptHighTPC[5]      = { 0., 0., 2.0,  2.0,   2.0  };
    //TOF----
    Double_t ptLowTOF[5]       = { 0., 0., 0.2,  0.625,  1.1  };
    Double_t ptHighTOF[5]      = { 0., 0., 2.0,  2.0,    2.0  };
    //TPCTOF----------
    Double_t ptLowTPCTOF[5]    = { 0., 0., 0.65, 0.69,   0.8  };
    Double_t ptHighTPCTOF[5]   = { 0., 0., 2.0,  2.00,   2.0  };
    
    
    //----------------------------------ITS PID---------------------------------------
    if (fPIDResponse->CheckPIDStatus((AliPIDResponse::EDetector)AliPIDResponse::kITS, track) == AliPIDResponse::kDetPidOk)
    {
        pid[0] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)AliPIDResponse::kITS, track, fParticleSpecies);
        
        if (TMath::Abs(pid[0]) < fNSigmaMaxITS) isAcceptedITS = kTRUE;
        Double_t nSigmaPion = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kPion));
        Double_t nSigmaKaon = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kKaon));
        Double_t nSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kElectron));
        
        if (TMath::Abs(pid[0]) > nSigmaEl) isAcceptedITS = kFALSE;
        
        if( fPidStrategy == 2){
            if( TMath::Abs(pid[0]) > nSigmaPion) isAcceptedITS = kFALSE;
            if( TMath::Abs(pid[0]) > nSigmaKaon) isAcceptedITS = kFALSE;
        }
        
    }
    
    //----------------------------------TPC PID------------------------------------------------
    if (fPIDResponse->CheckPIDStatus((AliPIDResponse::EDetector)AliPIDResponse::kTPC, track ) == AliPIDResponse::kDetPidOk){
        
        pid[1] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)AliPIDResponse::kTPC, track, fParticleSpecies);
        
        if (fParticleSpecies == 2){ //kaon
            if (track->Pt() > 0.525 && track->Pt() < 1.5) {
                if (TMath::Abs(pid[1]) < 2.) // cut on nsigma
                    isAcceptedTPC = kTRUE;
            }
        }// Pion

        if (fParticleSpecies == 3){ //kaon
            if (track->Pt() > 0.525 && track->Pt() < 1.5) {
                if (TMath::Abs(pid[1]) < 2.) // cut on nsigma
                    isAcceptedTPC = kTRUE;
            }
    }// kaon
        else {
            if (TMath::Abs(pid[1]) < fNSigmaMaxTPC) isAcceptedTPC = kTRUE;
        }
        
        Double_t nSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron));
        Double_t nSigmaPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kPion));
        Double_t nSigmaKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kKaon));
        
        if (TMath::Abs(pid[1]) > nSigmaEl) isAcceptedTPC = kFALSE;
        
        if (fPidStrategy == 2){
            if (pt > 0.9) {
                if (TMath::Abs(pid[1]) > nSigmaPion) isAcceptedTPC = kFALSE;
                if (TMath::Abs(pid[1]) > nSigmaKaon) isAcceptedTPC = kFALSE;
            }
        }
    }
    
    //-------------------------------TOF-------------------------------------------
    
    if (fPIDResponse->CheckPIDStatus((AliPIDResponse::EDetector)AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk) {
        pid[2] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)AliPIDResponse::kTOF, track, fParticleSpecies);
        hasPIDTOF = kTRUE;
        if (TMath::Abs(pid[2]) < fNSigmaMaxTOF) isAcceptedTOF = kTRUE;
        
        Double_t nSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kElectron));
        Double_t nSigmaPion = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kPion));
        Double_t nSigmaKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kKaon));
        
        if (TMath::Abs(pid[2]) > nSigmaEl) isAcceptedTOF = kFALSE;
        
        if( fPidStrategy == 2){
            if( pt > 0.9 ){
                if (TMath::Abs(pid[2]) > nSigmaPion) isAcceptedTOF = kFALSE;
                if (TMath::Abs(pid[2]) > nSigmaKaon) isAcceptedTOF = kFALSE;
            }
        }
        
    }
    
    
    if (fIsMC && isAcceptedTOF) {
        Int_t tofLabel[3];
        if (track->InheritsFrom("AliESDtrack")) {
            (dynamic_cast<AliESDtrack*>(track))->GetTOFLabel(tofLabel);
        } else if (track->InheritsFrom("AliAODTrack")) {
            (dynamic_cast<AliAODTrack*>(track))->GetTOFLabel(tofLabel);
        }
        
        Bool_t hasMatchTOF = kTRUE;
        if (TMath::Abs(track->GetLabel()) != TMath::Abs(tofLabel[0]) || tofLabel[1] > 0) {
            hasMatchTOF = kFALSE;
        }
        
        if(fIsAOD) {
            //------
        } else {
            TParticle *matchedTrack = fMCStack->Particle(TMath::Abs(tofLabel[0]));
            if (TMath::Abs(matchedTrack->GetFirstMother()) == TMath::Abs(track->GetLabel()))
                hasMatchTOF = kTRUE;
        }
        isAcceptedTOF = hasMatchTOF;
    }
    
    //----------------COmbined PID-----------------------------
    if (fParticleSpecies == 2){//for Pion: TPC+TOF---
        
        if(fPidStrategy == 0){
            isAccepted = isAcceptedTPC && isAcceptedTOF;
        }
        else if( fPidStrategy == 1){
            Double_t nsigCombined = TMath::Sqrt( pid[1]*pid[1] +  pid[2]*pid[2] );
            if( nsigCombined < fNSigmaMaxTOF ) isAccepted = kTRUE;
        }
        else if( fPidStrategy == 2){
            isAccepted = isAcceptedTOF;
        }
    }
    else if( fParticleSpecies == 3){//for kaon: TPC and/or TOF
        
        if ( pt > ptLowTPC[fParticleSpecies] && pt < ptHighTPC[fParticleSpecies] ) isAccepted = isAcceptedTPC;
        else isAccepted =  isAcceptedTOF;
    }
    
    else if( fParticleSpecies == 4){//for proton
        
        if(fPidStrategy == 0){
            //ITS+TPC and TPC+TOF
            if( pt >= ptLowITS[fParticleSpecies] && pt <= ptHighITS[fParticleSpecies] ) isAccepted = isAcceptedITS && isAcceptedTPC;
            else isAccepted = isAcceptedTPC && isAcceptedTOF;
        }
        else if( fPidStrategy == 1){
            
            if( pt >= ptLowITS[fParticleSpecies] && pt <= ptHighITS[fParticleSpecies] ) {
                Double_t nsigCompITSTPC = TMath::Sqrt( pid[1]*pid[1] +  pid[0]*pid[0] );
                if( nsigCompITSTPC < fNSigmaMaxTPC ) isAccepted = kTRUE;
            }
            else {
                Double_t nsigCompTPCTOF = TMath::Sqrt( pid[1]*pid[1] +  pid[2]*pid[2] );
                if( nsigCompTPCTOF < fNSigmaMaxTPC ) isAccepted = kTRUE;
            }
        }
        else if( fPidStrategy == 2){
            isAccepted = isAcceptedTOF && isAcceptedTPC;
        }
        else if( fPidStrategy == 3){
            if( pt >= 0.3 && pt <= 0.575 ) isAccepted = isAcceptedITS && isAcceptedTPC;
            else if( pt >= 0.825 && pt < 2.0 ) isAccepted = isAcceptedTPC && isAcceptedTOF;
            else isAccepted =  isAcceptedTPC;
            
        }
        else if( fPidStrategy == 4){
            if( pt >= 0.825 && pt < 2.0 ) isAccepted = isAcceptedTPC && isAcceptedTOF;
            else isAccepted =  isAcceptedTPC;
        }
        
    }
    
    
    delete [] pid;
    return isAccepted;
    
    
} // end

//-----------------------------------------------------------------------------------------------
