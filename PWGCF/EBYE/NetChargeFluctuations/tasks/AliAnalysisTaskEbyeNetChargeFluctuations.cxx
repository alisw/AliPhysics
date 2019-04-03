#include "Riostream.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "THn.h"
#include "TTree.h"
#include "TList.h"
#include "TObjArray.h"
#include "TTreeStream.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliVEvent.h"
#include "AliAODHeader.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliHeader.h"
#include "AliMultSelection.h"
#include "AliMultSelectionBase.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliVEventHandler.h"

#include "AliAnalysisTaskEbyeNetChargeFluctuations.h"

class AliAnalysisTaskEbyeNetChargeFluctuations;

#include "TFile.h"
#include <iostream>
#include <fstream>

using namespace std;

ClassImp(AliAnalysisTaskEbyeNetChargeFluctuations)

//================================================
//--------------Constructor-----------------------
//================================================

AliAnalysisTaskEbyeNetChargeFluctuations::AliAnalysisTaskEbyeNetChargeFluctuations()
: AliAnalysisTaskSE(),
fAOD(0x0),
fOutputList(0),
fTreeSRedirector(0x0),
fEtaDown(-0.8),
fEtaUp(0.8),
fHistPosEffMatrixRec(0),
fHistNegEffMatrixRec(0),
fHistPosEffMatrixGen(0),
fHistNegEffMatrixGen(0),
fHistCentralityMultSelection(0),
fHistVertexNconributors(0),
fHistVertexStats(0),
fHistZVertexCent(0),
fEventStatistics(0),
fHistClustersTPC(0),
fHistChi2perNDF(0),
fHistVx(0),
fHistVy(0),
fHistVz(0),
fRunNumber(-1),
fArrayMC(0),
hTrackPt(0),
hTrackPhi(0),
hTrackEta(0),
hTrackPtallrec(0),
hTrackPhiallrec(0),
hTrackEtaallrec(0),
hGenPt(0),
hGenPhi(0),
hGenEta(0),
fEventCuts(0),
EventNumber(0),
//fnCentbinsData(10),
//fcentDownArr(0),
//fcentUpArr(0),
fEta(0),
fpT(0),
fPhi(0),
fCharge(0),
fEtaMCall(0),
fpTMCall(0),
fPhiMCall(0),
fChargeMCall(0),
fCentrality(0),
fhCent(0),
fPos(0),
fNeg(0),
fMomentsCross(0),
fMomentsPos(0),
fMomentsNeg(0),
fEtaMCgen(0),
fpTMCgen(0),
fPhiMCgen(0),
fChargegen(0),
fEtaMCallgen(0),
fpTMCallgen(0),
fPhiMCallgen(0),
fChargeMCallgen(0),
//fEtaMC(0),
//fpTMC(0),
genPos(0),
genNeg(0),
allgenPos(0),
allgenNeg(0),
Nch(0),
genNch(0),
recNch(0),
allrecNch(0),
allgenNch(0),
genMomentsCross(0),
genMomentsPos(0),
genMomentsNeg(0),
allgenMomentsCross(0),
allgenMomentsPos(0),
allgenMomentsNeg(0),
recPos(0),
recNeg(0),
recMomentsCross(0),
recMomentsPos(0),
recMomentsNeg(0),
allrecPos(0),
allrecNeg(0),
allrecMomentsCross(0),
allrecMomentsPos(0),
allrecMomentsNeg(0),
fAnalysisType(0),
fTree(0x0),
fTreeTrackCuts(0x0),
fTreeMCTrackCuts(0x0),
fTreeMCrec(0x0),
fTreeMCallrec(0x0),
fTreeMCallgen(0x0),
fTreeMCgen(0x0)
{
    // default constructor
    
}

//-------------------------------------------------

AliAnalysisTaskEbyeNetChargeFluctuations::AliAnalysisTaskEbyeNetChargeFluctuations(const char* name)
:AliAnalysisTaskSE(name),
fAOD(0x0),
fOutputList(0),
fTreeSRedirector(0x0),
fEtaDown(-0.8),
fEtaUp(0.8),
fHistPosEffMatrixRec(0),
fHistNegEffMatrixRec(0),
fHistPosEffMatrixGen(0),
fHistNegEffMatrixGen(0),
fHistCentralityMultSelection(0),
fHistVertexNconributors(0),
fHistVertexStats(0),
fHistZVertexCent(0),
fEventStatistics(0),
fHistClustersTPC(0),
fHistChi2perNDF(0),
fHistVx(0),
fHistVy(0),
fHistVz(0),
fRunNumber(-1),
fArrayMC(0),
hTrackPt(0),
hTrackPhi(0),
hTrackEta(0),
hTrackPtallrec(0),
hTrackPhiallrec(0),
hTrackEtaallrec(0),
hGenPt(0),
hGenPhi(0),
hGenEta(0),
fEventCuts(0),
EventNumber(0),
//fnCentbinsData(10),
//fcentDownArr(0),
//fcentUpArr(0),
fEta(0),
fpT(0),
fPhi(0),
fCharge(0),
fEtaMCall(0),
fpTMCall(0),
fPhiMCall(0),
fChargeMCall(0),
fCentrality(0),
fhCent(0),
fPos(0),
fNeg(0),
fMomentsCross(0),
fMomentsPos(0),
fMomentsNeg(0),
fEtaMCgen(0),
fpTMCgen(0),
fPhiMCgen(0),
fChargegen(0),
fEtaMCallgen(0),
fpTMCallgen(0),
fPhiMCallgen(0),
fChargeMCallgen(0),
//fEtaMC(0),
//fpTMC(0),
genPos(0),
genNeg(0),
allgenPos(0),
allgenNeg(0),
Nch(0),
genNch(0),
recNch(0),
allrecNch(0),
allgenNch(0),
genMomentsCross(0),
genMomentsPos(0),
genMomentsNeg(0),
allgenMomentsCross(0),
allgenMomentsPos(0),
allgenMomentsNeg(0),
recPos(0),
recNeg(0),
recMomentsCross(0),
recMomentsPos(0),
recMomentsNeg(0),
allrecPos(0),
allrecNeg(0),
allrecMomentsCross(0),
allrecMomentsPos(0),
allrecMomentsNeg(0),
fAnalysisType(0),
fTree(0x0),
fTreeTrackCuts(0x0),
fTreeMCTrackCuts(0x0),
fTreeMCrec(0x0),
fTreeMCallrec(0x0),
fTreeMCallgen(0x0),
fTreeMCgen(0x0)
{
    
    cout << "============================================================" << endl;
    cout << "=================End of Constructor====================="  << endl;
    cout << "============================================================" << endl;
    
    //constructor
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TTree::Class());
    DefineOutput(4, TTree::Class());
    DefineOutput(5, TTree::Class());
    DefineOutput(6, TTree::Class());
    
}
//==========================================================
//--------------Destructor----------------------------------
//==========================================================
AliAnalysisTaskEbyeNetChargeFluctuations::~AliAnalysisTaskEbyeNetChargeFluctuations()
{
    //destructor
    if(fOutputList) 			        delete fOutputList;
    if (fHistPosEffMatrixRec) 		    delete fHistPosEffMatrixRec;
    if (fHistNegEffMatrixRec) 		    delete fHistNegEffMatrixRec;
    if (fHistPosEffMatrixGen) 		    delete fHistPosEffMatrixGen;
    if (fHistNegEffMatrixGen) 		    delete fHistNegEffMatrixGen;
    if (fhCent)               		    delete fhCent;
    if (fHistCentralityMultSelection)   delete fHistCentralityMultSelection;
    if (fEventStatistics)               delete fEventStatistics;
}
//============================================================
//--------------UserCreateOutputObjects-----------------------
//============================================================
void AliAnalysisTaskEbyeNetChargeFluctuations::UserCreateOutputObjects()
{
    
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    // Event statistics
    fEventStatistics = new TH1I("fEventStatistics","",10,0,10);
    fOutputList->Add(fEventStatistics);
    
    // ****************** Efficiency matrix histograms ************************
    // 0 --> particle type: 0 positive, 2, negative
    // 1 --> Centrality
    // 2 --> eta
    // 3 --> momentum
    // 4 --> phi
    const Int_t ndim=5;
    
    //						   0, 1,   2,   3,   4
    Int_t nbins0[ndim]      = {2, 8,   48,  16 , 40  };
    Double_t xmin0[ndim]    = {0, 0.0, 0.2,-0.8, 0.  };
    Double_t xmax0[ndim]    = {2, 80.0,5.0, 0.8, 6.25};
    fHistPosEffMatrixRec    = new THnF("fHistPosEffMatrixRec","fHistPosEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixRec    = new THnF("fHistNegEffMatrixRec","fHistNegEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistPosEffMatrixGen    = new THnF("fHistPosEffMatrixGen","fHistPosEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixGen    = new THnF("fHistNegEffMatrixGen","fHistNegEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    
    TString axisNameEff[ndim]  = {"particle"      ,"Centrality"     ,"momentum"      ,"eta"  ,"phi"};
    TString axisTitleEff[ndim] = {"particle type" ,"Centrality (%)" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta" ,"#phi"};
    for (Int_t iEff=0; iEff<ndim;iEff++){
        fHistPosEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
        fHistNegEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
        fHistPosEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
        fHistNegEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
    }
    fOutputList->Add(fHistPosEffMatrixRec);
    fOutputList->Add(fHistNegEffMatrixRec);
    fOutputList->Add(fHistPosEffMatrixGen);
    fOutputList->Add(fHistNegEffMatrixGen);
    
    // single-track QA plots
    hTrackPt = new TH1F("hTrackPt","track p_{T};p_{T} (GeV/c);",120,0.0,6.0);
    hTrackPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hTrackPt->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hTrackPt);
    hTrackPhi = new TH1F("hTrackPhi","track #varphi;#varphi;",160,0,TMath::TwoPi());
    hTrackPhi->GetXaxis()->SetTitle("#Phi ");
    hTrackPhi->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hTrackPhi);
    hTrackEta = new TH1F("hTrackEta","track #eta;#eta;",32, -0.8, 0.8);
    hTrackEta->GetXaxis()->SetTitle("#eta ");
    hTrackEta->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hTrackEta);
    
    hTrackPtallrec = new TH1F("hTrackPtallrec","track p_{T};p_{T} (GeV/c);",200,0.0,10.0);
    hTrackPtallrec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hTrackPtallrec->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hTrackPtallrec);
    hTrackPhiallrec = new TH1F("hTrackPhiallrec","track #varphi;#varphi;",160,0,TMath::TwoPi());
    hTrackPhiallrec->GetXaxis()->SetTitle("#Phi ");
    hTrackPhiallrec->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hTrackPhiallrec);
    hTrackEtaallrec = new TH1F("hTrackEtaallrec","track #eta;#eta;",32, -0.8, 0.8);
    hTrackEtaallrec->GetXaxis()->SetTitle("#eta ");
    hTrackEtaallrec->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hTrackEtaallrec);
    
    hGenPt = new TH1F("hGenPt","generated p_{T};p_{T} (GeV/c);",200,0.0,10.0);
    hGenPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hGenPt->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hGenPt);
    hGenPhi = new TH1F("hGenPhi","generated #varphi;#varphi;",160,0,TMath::TwoPi());
    hGenPhi->GetXaxis()->SetTitle("#Phi ");
    hGenPhi->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hGenPhi);
    hGenEta = new TH1F("hGenEta","generated #eta;#eta;",32, -0.8, 0.8);
    hGenEta->GetXaxis()->SetTitle("#eta ");
    hGenEta->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hGenEta);
    
    // AliVEvent QA Cuts
    fEventCuts.AddQAplotsToList(fOutputList);
    
    //Vertex distributions
    fHistVx = new TH1D("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm)",400,-20, 20);
    fOutputList->Add(fHistVx);
    fHistVy = new TH1D("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm)",400,-20, 20);
    fOutputList->Add(fHistVy);
    fHistVz = new TH1D("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm)",400,-20,20);
    fOutputList->Add(fHistVz);
    
    fHistCentralityMultSelection = new TH1D("fHistCentralityMultSelection","Centrality Percentile ;Centrality;Entries",80,0.,80.);  //new sk
    fHistCentralityMultSelection->GetXaxis()->SetTitle("Centrality (%)");
    fOutputList->Add(fHistCentralityMultSelection);
    
    fHistZVertexCent = new TH2F("fHistZVertexCent"," Vz of primary vertex Vs Centrality; V_{Z} {cm} ; Centrality {%}", 60, -15, 15,100,0,100);
    fOutputList->Add(fHistZVertexCent);
    
    fHistClustersTPC = new TH1D("fHistClustersTPC","N Clusters TPC;N_{TPC clusters};Entries",181,-0.5,180.5);
    fOutputList->Add(fHistClustersTPC);
    
    fHistChi2perNDF = new TH1D("fHistChi2perNDF",";#chi^{2} / ndf;n tracks", 700, -1, 6);
    fOutputList->Add(fHistChi2perNDF);
    
    // ************************************************************************
    
    // Tree for pt eta and centrality checks
    fTreeSRedirector    = new TTreeSRedirector();
    fTree               = ((*fTreeSRedirector)<<"Realdata").GetTree();
    fTreeMCrec          = ((*fTreeSRedirector)<<"MCrec").GetTree();
    fTreeMCallrec       = ((*fTreeSRedirector)<<"MCallrec").GetTree();
    fTreeMCgen          = ((*fTreeSRedirector)<<"MCgen").GetTree();
    fTreeMCallgen       = ((*fTreeSRedirector)<<"MCallgen").GetTree();
    
    cout << "============================================================" << endl;
    cout << "=================End of User Create Object=================="  << endl;
    cout << "============================================================" << endl;
    
    //Send data to the container
    PostData(1, fOutputList);
    PostData(2, fTree);
    PostData(3, fTreeMCrec);
    PostData(4, fTreeMCallrec);
    PostData(5, fTreeMCgen);
    PostData(6, fTreeMCallgen);
    
}
//============================================================
//----------------------UserExec------------------------------
//============================================================
void AliAnalysisTaskEbyeNetChargeFluctuations::UserExec(Option_t *){
    
    if(fAnalysisType == 0) {
        
        doAODEvent();
        
    } //====================Read data AOD-analysis
    
    else if(fAnalysisType == 1) {
        
        doMCAODEvent();
        FillMCEffMatrix();
        
    } //====================Monte Carlo AOD-analysis
    
    else return;
    
}
//============================================================
//--------------------- FUNCTIONS-----------------------------
//============================================================

void AliAnalysisTaskEbyeNetChargeFluctuations::doAODEvent(){
    
    fEventStatistics->Fill("before cuts",1);
    
    //  AliVEvent *event = fInputEvent();
    if (!fInputEvent) return;
    
    fEventStatistics->Fill("after event check",1);
    
    fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
    
    if(!fAOD) { Printf("ERROR: fAOD not available"); return; }
    fEventStatistics->Fill("after aod check",1);
    
    //Physics selection
    UInt_t fSelectMask= fInputHandler->IsEventSelected();
    Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7; // from https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsEventProp
    if(!isINT7selected)return ;
    
    fEventStatistics->Fill("physics selection",1);
    
    AliAODVertex *fPrimaryVtx = (AliAODVertex*)fInputEvent->GetPrimaryVertex();
    if (!fPrimaryVtx){
        printf ("ERROR: no primary vertex\n");
        return;
    }
    
    Double_t xv=fPrimaryVtx->GetX();
    Double_t yv=fPrimaryVtx->GetY();
    Double_t zv=fPrimaryVtx->GetZ();
    
    fHistVx->Fill(xv);
    fHistVy->Fill(yv);
    fHistVz->Fill(zv);
    
    fEventStatistics->Fill("found primary vertex",1);
    
    //===========================Centrality calculation =====================
    fCentrality = -2;
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    fEventStatistics->Fill("found MultSelection object",1);
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
    
    if (TMath::Abs(zv) > 10.0) return;
    fEventStatistics->Fill("vz cut",1);
    
    if(fCentrality < 0 || fCentrality >= 80) return;
    fHistCentralityMultSelection->Fill(fCentrality);
    fEventStatistics->Fill("centrality selection",1);
    
    fHistZVertexCent->Fill(zv, fCentrality);
    
    if (!fEventCuts.AcceptEvent(fInputEvent)) return;
    fEventStatistics->Fill("AliEventCuts",1);
    
    // =========Different Acceptance cuts for Eta momentum and Centrality=============
    Double_t momDownArray[4] = {0.2, 0.2, 0.6, 0.6};
    Double_t momUpArray[4]   = {2.0, 5.0, 2.0, 5.0};
    Double_t etaDownArray[8] = {-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8};
    Double_t etaUpArray[8]   = { 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8};
    //case1:
    Double_t centDownArray[9]= {0., 5.,  10., 20., 30., 40., 50., 60., 70.};
    Double_t centUpArray[9]  = {5., 10., 20., 30., 40., 50., 60., 70., 80.};
    
    //case 3:
    //    Double_t centDownArray[18] =  {0,2.5,5.0,7.5,10,15,20,25,30,35,40,45,50,55,60,65,70,75};
    //    Double_t centUpArray[18] = {2.5,5.0,7.5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80};
    
    for (Int_t imom=0; imom<4; imom++){
        for (Int_t ieta=0; ieta<8; ieta++){
            for (Int_t icent=0; icent<9; icent++){
                
                Double_t centBin = (centDownArray[icent]+centUpArray[icent])/2.;
                
                Double_t etaBin  = (TMath::Abs(etaDownArray[ieta])+etaUpArray[ieta]);
                
                // Initialize the positive and negative particles
                fPos = 0;
                fNeg = 0;
                Int_t subsample1 = 0 ;
                Int_t trCount=0, TotalCharge = 0;
                
                subsample1 = EventNumber%30;
                //=============================== track loop ==============================
                Int_t iTracks(fAOD->GetNumberOfTracks());
                for(Int_t i(0); i < iTracks; i++) {
                    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
                    
                    if(!track) {
                        AliError(Form("ERROR: Could not retrieve AODtrack %d",i));
                        continue;
                    }
                    
                    if(!AcceptTrack(track)) continue;
                
                    if ((track->Eta()<etaDownArray[ieta]) || (track->Eta()>etaUpArray[ieta])) continue;  // eta Cut
                    
                    fpT     = track->Pt();
                    fEta    = track->Eta();
                    fPhi    = track->Phi();
                    fCharge = track->Charge();
                    
                    hTrackPt->Fill(fpT);
                    hTrackPhi->Fill(fPhi);
                    hTrackEta->Fill(fEta);
                    
                    //===================apply cuts on pt eta and centrality=====================
                   if ((fpT>=momDownArray[imom])
                        &&(fpT<momUpArray[imom])
                        &&(fCentrality>=centDownArray[icent])
                        &&(fCentrality<centUpArray[icent])){

                        // calculate first moments
                        if (fCharge < 0 || fCharge > 0) Nch++;
                        if(fCharge > 0) fPos++;
                        if(fCharge < 0) fNeg++;
                        
                        trCount++;
                    }
                    
                } //============end of track loop
                
                fRunNumber = fInputEvent->GetRunNumber();
                
                fEventStatistics->Fill("after event loop",1);
                
                // calculate second moments
                fMomentsCross =  fPos * fNeg;
                fMomentsPos    = fPos * fPos;
                fMomentsNeg    = fNeg * fNeg;
                
                // Fill the pt eta and centrality check in the tree
               if ( trCount>0 ){
                    if(!fTreeSRedirector) return;
                    (*fTreeSRedirector)<<"Realdata"<<
                    // upper edge of momentum bin
                    "pos="          << fPos                <<       // positive particles
                    "neg="          << fNeg                <<       // negative particles
                    "isample1="     << subsample1          <<       // sample id for subsample method
                    "momentscross=" << fMomentsCross       <<       // second cross moments
                    "momentspos="   << fMomentsPos         <<       // second pos moments
                    "momentsneg="   << fMomentsNeg         <<       // second neg moments
                    "etaDown="      << etaDownArray[ieta]  <<       // lower edge of eta bin
                    "etaUp="        << etaUpArray[ieta]    <<       // upper edge of eta bin
                    "centDown="     << centDownArray[icent]<<       // lower edge of cent bin
                    "centUp="       << centUpArray[icent]  <<       // upper edge of cent bin
                    "momDown="      << momDownArray[imom]  <<       // lower edge of mom bin
                    "momUp="        << momUpArray[imom]    <<       // upper edge of mom bin
                    "ieta="         << ieta                <<
                    "icent="        << icent               <<
                    "imom="         << imom                <<
                    "centBin="      << centBin             <<       // centrality bining
                    "etabin="       << etaBin              <<       // eta bin
                    "fRunNumber="   << fRunNumber          <<       // Run Number
                    "fPt="          << fpT                 <<
                    "fEta="         << fEta                <<
                    "fPhi="         << fPhi                <<
                    "\n";
                } //===========tree filling========
                
            } //===========end of momentum loop
        } //===========end of eta loop
    } //===========end of centralityy loop
    
    EventNumber++;
    PostData(1, fOutputList);
    PostData(2, fTree);
    PostData(3, fTreeMCrec);
    PostData(4, fTreeMCallrec);
    PostData(5, fTreeMCgen);
    PostData(6, fTreeMCallgen);
    
}

void AliAnalysisTaskEbyeNetChargeFluctuations::doMCAODEvent(){
    
    fEventStatistics->Fill("before cuts",1);
    
    AliVEvent *event = InputEvent();
    if (!event) { Printf("UserExec: NO EVENT FOUND!"); return; }
    
    fEventStatistics->Fill("after event check",1);
    
    AliAODEvent* fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) { Printf("ERROR: fAOD not available"); return; }
    
    AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    fEventStatistics->Fill("after aod check",1);
    
    //Physics selection
    UInt_t fSelectMask= eventHandler->IsEventSelected();
    Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7; // from https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsEventProp
    if(!isINT7selected)return ;
    
    fEventStatistics->Fill("physics selection",1);
    
    AliAODVertex *fPrimaryVtx = (AliAODVertex*)InputEvent()->GetPrimaryVertex();
    if (!fPrimaryVtx){
        printf ("ERROR: no primary vertex\n");
        return;
    }
    
    Double_t xv=fPrimaryVtx->GetX();
    Double_t yv=fPrimaryVtx->GetY();
    Double_t zv=fPrimaryVtx->GetZ();
    
    fHistVx->Fill(xv);
    fHistVy->Fill(yv);
    fHistVz->Fill(zv);
    
    fEventStatistics->Fill("found primary vertex",1);
    
    //===========================Centrality calculation =====================
    fCentrality = -2;
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    fEventStatistics->Fill("found MultSelection object",1);
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
    
    if (TMath::Abs(zv) > 10.0) return;
    fEventStatistics->Fill("vz cut",1);
    
    if(fCentrality < 0 || fCentrality >= 80) return;
    fHistCentralityMultSelection->Fill(fCentrality);
    fEventStatistics->Fill("centrality selection",1);
    
    fHistZVertexCent->Fill(zv, fCentrality);
        
    fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fArrayMC)
        AliFatal("No array of MC particles found !!!");
    
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
        AliError("ERROR: Could not retrieve MC event");
        return;
    }
    
    //    fEventStatistics->Fill("AliMCEventCuts",1);
    
    // =========Different Acceptance cuts for Eta momentum and Centrality=============
    Double_t momDownArray[4] = {0.2, 0.2, 0.6, 0.6};
    Double_t momUpArray[4]   = {2.0, 5.0, 2.0, 5.0};
    Double_t etaDownArray[8] = {-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8};
    Double_t etaUpArray[8]   = { 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8};
    //case1:
    Double_t centDownArray[9]= {0., 5.,  10., 20., 30., 40., 50., 60., 70.};
    Double_t centUpArray[9]  = {5., 10., 20., 30., 40., 50., 60., 70., 80.};
    
    //case 2:
    //  Double_t centDownArray[18] =  {0,2.5,5.0,7.5,10,15,20,25,30,35,40,45,50,55,60,65,70,75};
    //  Double_t centUpArray[18] = {2.5,5.0,7.5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80};
    
    
    for (Int_t imom=0; imom<4; imom++){
        for (Int_t ieta=0; ieta<8; ieta++){
            for (Int_t icent=0; icent<9; icent++){
                
                // -----------------------------------------------------------------------------------------
                // cout<< " ------------------reconstructed MC particles------------------------" << endl;
                // -----------------------------------------------------------------------------------------
                
                Double_t centBin = (centDownArray[icent]+centUpArray[icent])/2.;
                Double_t etaBin  = (TMath::Abs(etaDownArray[ieta])+etaUpArray[ieta]);
                
                // Initialize
                Int_t trCountMCrec=0, subsample1 = 0;
                Int_t recPos = 0, recNeg = 0;
                Int_t recTotalCharge = 0, primPhysicalCount = 0;
                
                subsample1 = EventNumber%30;
                
                // Loop over the reconstructed tracks
                Int_t iTracks(fAOD->GetNumberOfTracks());
                for(Int_t i(0); i < iTracks; i++) { // track loop
                    
                    // Initialize the variables
                    pdg = 0, pdgMom=0, pdgMomPhysicalPrim = 0, pdgMomPrim = 0, pdgPhysicalPrim = 0, pdgPrim = 0;
                    
                    // get the AOD Tracks
                    AliAODTrack* trackReal = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
                    if(!trackReal) {
                        AliError(Form("ERROR: Could not retrieve AODtrack %d",i));
                        continue;
                    }
                    
                    // Track cuts from detector
                    if(!AcceptTrack(trackReal)) continue;
                    
                    fpT     = trackReal->Pt();
                    fEta    = trackReal->Eta();
                    fPhi    = trackReal->Phi();
                    fCharge = trackReal->Charge();
                    
                    // MC track cuts
                    if ((fEta < etaDownArray[ieta]) || (fEta > etaUpArray[ieta])) continue;  // eta cut
                    
                    AliAODMCParticle* trackMCrec = (AliAODMCParticle*)fArrayMC->At(TMath::Abs(trackReal->GetLabel()));
                    if(!trackMCrec) continue;
                    
                    if (!trackMCrec->IsPhysicalPrimary()) continue;
                    
                    Int_t pdg = trackMCrec->GetPdgCode();
                    
                    // get the mother of the daughter particle
                    Int_t momlab  = trackMCrec->GetMother();
                    AliVParticle *aodMother = mcEvent->GetTrack(momlab);
                    if(aodMother) pdgMom  = aodMother->PdgCode();   // pdg of mother particle
                    
                    //===================apply cuts on pt eta and centrality=====================
                    if ((fpT >=momDownArray[imom])
                        &&(fpT <momUpArray[imom])
                        &&(fCentrality>=centDownArray[icent])
                        &&(fCentrality<centUpArray[icent])){
                        
                        // Select Physical primary particles and fill it on the tree
                        if (trackMCrec->IsPhysicalPrimary()) {
                            pdgMomPhysicalPrim = pdgMom;
                            primPhysicalCount++;
                            // calculate first moments
                            if (fCharge < 0 || fCharge > 0) recNch++;
                            if(fCharge > 0) recPos++;
                            if(fCharge < 0) recNeg++;
                            
                            trCountMCrec++;
                        }
                        
                    } // end loop of apply cuts
                    
                    hTrackPt->Fill(fpT);
                    hTrackPhi->Fill(fPhi);
                    hTrackEta->Fill(fEta);
                    
                } //============end of track loop====================
                
                fRunNumber = fInputEvent->GetRunNumber();
                
                fEventStatistics->Fill("after event loop",1);
                
                // calculate second moments
                recMomentsCross  = recPos * recNeg;
                recMomentsPos    = recPos * recPos;
                recMomentsNeg    = recNeg * recNeg;
                
                // Tree for all the cuts variables
                if ( trCountMCrec>0 ){
                    if(!fTreeSRedirector) return;
                    (*fTreeSRedirector)<<"MCrec"<<
                    "pos="                 << recPos                  <<    // MC gen positive particles
                    "neg="                 << recNeg                  <<    // MC gen negative particles
                    "isample1="      	   << subsample1              <<    // sample id for subsample method
                    "momentscross="        << recMomentsCross         <<    // MC gen second moments of pos and neg
                    "momentspos="          << recMomentsPos           <<    // MC gen second moments of pos
                    "momentsneg="          << recMomentsNeg           <<    // MC gen second moments of neg
                    "etaDown="             << etaDownArray[ieta]      <<    // lower edge of eta bin
                    "etaUp="               << etaUpArray[ieta]        <<    // upper edge of eta bin
                    "centDown="            << centDownArray[icent]    <<    // lower edge of cent bin
                    "centUp="              << centUpArray[icent]      <<    // upper edge of cent bin
                    "momDown="             << momDownArray[imom]      <<    // lower edge of mom bin
                    "momUp="               << momUpArray[imom]        <<    // upper edge of mom bin
                    "ieta="                << ieta                    <<    // eta index loop
                    "icent="               << icent                   <<    // cent index loop
                    "imom="                << imom                    <<    // mom index loop
                    "centBin="             << centBin                 <<    // centrality bining
                    "etabin="              << etaBin                  <<    // eta bin
                    "fRunNumber="          << fRunNumber              <<    // Run Number
                    "\n";
                    
                } // trackcuts tree fill========
                
                // -----------------------------------------------------------------------------------------
                // cout<< " -----------All reconstructed MC particles------------------------" << endl;
                // -----------------------------------------------------------------------------------------
                
                
                // Initialize
                Int_t trCountMCallrec=0, sampleNo = 0;;
                Int_t allrecPos = 0, allrecNeg = 0;
                Int_t allrecTotalCharge = 0, allprimPhysicalCount = 0;
                
                // Loop over the reconstructed tracks
                Int_t nTracks(fAOD->GetNumberOfTracks());
                for(Int_t irec = 0; irec < nTracks; irec++) { // track loop
                    
                    // Initialize the variables
                    allpdg = 0, allpdgMom=0, allpdgMomPhysicalPrim = 0, allpdgMomPrim = 0, allpdgPhysicalPrim = 0, allpdgPrim = 0;
                    
                    // get the AOD Tracks
                    AliAODTrack* trackall = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(irec));
                    if(!trackall) {
                        AliError(Form("ERROR: Could not retrieve AODtrack %d",irec));
                        continue;
                    }
                    
                    // Track cuts from detector
                    if(!AcceptTrack(trackall)) continue;
                    
                    fpTMCall      = trackall->Pt();
                    fEtaMCall     = trackall->Eta();
                    fPhiMCall     = trackall->Phi();
                    fChargeMCall  = trackall->Charge();
                    
                    // MC track cuts
                    if ((fEtaMCall < etaDownArray[ieta]) || fEtaMCall > etaUpArray[ieta]) continue;  // eta cut
                    
                    AliAODMCParticle* trackMCallrec = (AliAODMCParticle*)fArrayMC->At(TMath::Abs(trackall->GetLabel()));
                    if(!trackMCallrec) continue;
                    
                    Int_t allpdg = trackMCallrec->GetPdgCode();
                    
                    // get the mother of the daughter particle
                    Int_t allmomlab  = trackMCallrec->GetMother();
                    AliVParticle *aodMother = mcEvent->GetTrack(allmomlab);
                    if(aodMother) allpdgMom  = aodMother->PdgCode(); // pdg of mother particle
                    
                    //===================apply cuts on pt eta and centrality=====================
                    if ((fpTMCall >= momDownArray[imom])
                        &&(fpTMCall < momUpArray[imom])
                        &&(fCentrality >= centDownArray[icent])
                        &&(fCentrality < centUpArray[icent])){
                        
                        // calculate first moments
                        if (fChargeMCall < 0 || fChargeMCall > 0) allrecNch++;
                        if(fChargeMCall > 0) allrecPos++;
                        if(fChargeMCall < 0) allrecNeg++;
                        trCountMCallrec++;
                       
                    } // end loop of apply cuts
                    
                    hTrackPtallrec->Fill(fpTMCall);
                    hTrackPhiallrec->Fill(fPhiMCall);
                    hTrackEtaallrec->Fill(fEtaMCall);
                    
                } //============end of track loop====================
                
                fRunNumber = fInputEvent->GetRunNumber();
                
                // calculate second moments
                allrecMomentsCross  = allrecPos * allrecNeg;
                allrecMomentsPos    = allrecPos * allrecPos;
                allrecMomentsNeg    = allrecNeg * allrecNeg;
                
                // Tree for all the cuts variables
                if ( trCountMCallrec>0 ){
                    if(!fTreeSRedirector) return;
                    (*fTreeSRedirector)<<"MCallrec"<<
                    "pos="                 << allrecPos               <<    // MC gen positive particles
                    "neg="                 << allrecNeg               <<    // MC gen negative particles
                    "isample1="      	   << subsample1              <<    // sample id for subsample method
                    "momentscross="        << allrecMomentsCross      <<    // MC gen second moments of pos and neg
                    "momentspos="          << allrecMomentsPos        <<    // MC gen second moments of pos
                    "momentsneg="          << allrecMomentsNeg        <<    // MC gen second moments of neg
                    "etaDown="             << etaDownArray[ieta]      <<    // lower edge of eta bin
                    "etaUp="               << etaUpArray[ieta]        <<    // upper edge of eta bin
                    "centDown="            << centDownArray[icent]    <<    // lower edge of cent bin
                    "centUp="              << centUpArray[icent]      <<    // upper edge of cent bin
                    "momDown="             << momDownArray[imom]      <<    // lower edge of mom bin
                    "momUp="               << momUpArray[imom]        <<    // upper edge of mom bin
                    "ieta="                << ieta                    <<    // eta index loop
                    "icent="               << icent                   <<    // cent index loop
                    "imom="                << imom                    <<    // mom index loop
                    "centBin="             << centBin                 <<    // centrality bining
                    "etabin="              << etaBin                  <<    // eta bin
                    "fRunNumber="          << fRunNumber              <<    // Run Number
                    "\n";
                    
                } // trackcuts tree fill========
                
                // -----------------------------------------------------------------------------------------
                // cout<< " ------------------Generated MC particles------------------------" << endl;
                // -----------------------------------------------------------------------------------------
                // Initialize the positive and negative particles
                genPos = 0, genNeg = 0;
                Int_t trCountMCgen = 0;
                Int_t genTotalCharge = 0, primPhysicalCountgen = 0;
                
                //=============================== track loop ==============================
                
                // Loop over the MC gen tracks
                Int_t nMCgenTracks = mcEvent->GetNumberOfTracks();
                
                for(Int_t iPart = 0; iPart < nMCgenTracks; iPart++) {
                    
                    // Initialize the variables
                    pdgMomPhysicalPrimgen = 0, pdgMomPrimgen = 0,pdgPhysicalPrimgen = 0, pdgPrimgen = 0;
                    
                    AliAODMCParticle *trackMCgen  = (AliAODMCParticle*)mcEvent->GetTrack(iPart);
                    if(!trackMCgen) continue;
                    
                    fpTMCgen     = trackMCgen->Pt();
                    fPhiMCgen    = trackMCgen->Phi();
                    fChargegen   = trackMCgen->Charge();
                    fEtaMCgen    = trackMCgen->Eta();
                    
                    if (fChargegen == 0) continue;
                    if (fpTMCgen < 0.2) continue;
                    
                    // MC eta cut
                    if ((fEtaMCgen < etaDownArray[ieta]) || (fEtaMCgen > etaUpArray[ieta])) continue;
                    
                    if (!trackMCgen->IsPhysicalPrimary()) continue;
                    Int_t pdggen=trackMCgen->GetPdgCode();
                    
                    // get the mother of the daughter particle
                    Int_t momlabgen  = trackMCgen->GetMother();
                    AliVParticle *aodMothergen = mcEvent->GetTrack(momlabgen);
                    if(aodMothergen) pdgMomgen = aodMothergen->PdgCode(); // pdg of mother particle
                    
                    //===================apply cuts on pt eta and centrality=====================
                    if ((fpTMCgen>=momDownArray[imom])
                        &&(fpTMCgen<momUpArray[imom])
                        &&(fCentrality>=centDownArray[icent])
                        &&(fCentrality<centUpArray[icent])){
                        
                        // Select Physical primary particles and fill it on the tree
                        if (trackMCgen->IsPhysicalPrimary()) {
                            pdgMomPhysicalPrimgen = pdgMomgen;
                            primPhysicalCountgen++;
                            // calculate first moments
                            if (fChargegen < 0 || fChargegen > 0) genNch++;
                            if(fChargegen > 0) genPos++;
                            if(fChargegen < 0) genNeg++;
                            trCountMCgen++;
                            
                        }
                    } // end loop of apply cuts
                    
                    hGenPt->Fill(fpTMCgen);
                    hGenPhi->Fill(fPhiMCgen);
                    hGenEta->Fill(fEtaMCgen);
                    
                } //============end of track loop====================
                
                fRunNumber = fInputEvent->GetRunNumber();
                
                // calculate second moments
                genMomentsCross  = genPos * genNeg;
                genMomentsPos    = genPos * genPos;
                genMomentsNeg    = genNeg * genNeg;
                
                
                // Fill the pt eta and centrality check in the tree
                if ( trCountMCgen>0 ){
                    if(!fTreeSRedirector) return;
                    (*fTreeSRedirector)<<"MCgen"<<
                    "pos="                 << genPos                    <<    // MC gen positive particles
                    "neg="                 << genNeg                    <<    // MC gen negative particles
                    "isample1="      	   << subsample1                <<    // sample id for subsample method
                    "momentscross="        << genMomentsCross           <<    // MC gen second moments of pos and neg
                    "momentspos="          << genMomentsPos             <<    // MC gen second moments of pos
                    "momentsneg="          << genMomentsNeg             <<    // MC gen second moments of neg
                    "etaDown="             << etaDownArray[ieta]        <<    // lower edge of eta bin
                    "etaUp="               << etaUpArray[ieta]          <<    // upper edge of eta bin
                    "centDown="            << centDownArray[icent]      <<    // lower edge of cent bin
                    "centUp="              << centUpArray[icent]        <<    // upper edge of cent bin
                    "momDown="             << momDownArray[imom]        <<    // lower edge of mom bin
                    "momUp="               << momUpArray[imom]          <<    // upper edge of mom bin
                    "ieta="                << ieta                      <<    // eta index loop
                    "icent="               << icent                     <<    // cent index loop
                    "imom="                << imom                      <<    // mom index loop
                    "centBin="             << centBin                   <<    // centrality bin
                    "etabin="              << etaBin                    <<    // eta bin
                    "fRunNumber="          << fRunNumber                <<    // Run Number
                    "\n";
                    
                } //===========tree filling========
                
                // -----------------------------------------------------------------------------------------
                // cout<< " ------------------ALL Generated MC particles------------------------" << endl;
                // -----------------------------------------------------------------------------------------
                // Initialize the positive and negative particles
                allgenPos = 0, allgenNeg = 0, allgenNch = 0;
                Int_t alltrCountMCgen = 0;
                
                //=============================== track loop ==============================
                
                // Loop over the MC gen tracks
                Int_t nMCallgenTracks = mcEvent->GetNumberOfTracks();
                
                for(Int_t iMC = 0; iMC < nMCallgenTracks; iMC++) {
                    
                    AliAODMCParticle *trackMCallgen  = (AliAODMCParticle*)mcEvent->GetTrack(iMC);
                    if(!trackMCallgen) continue;
                    
                    fpTMCallgen         = trackMCallgen->Pt();
                    fPhiMCallgen        = trackMCallgen->Phi();
                    fChargeMCallgen     = trackMCallgen->Charge();
                    fEtaMCallgen        = trackMCallgen->Eta();
                    
                    if (fChargeMCallgen == 0) continue;
                    if (fpTMCallgen < 0.2) continue;
                    
                    // MC eta cut
                    if ((fEtaMCallgen < etaDownArray[ieta]) || (fEtaMCallgen > etaUpArray[ieta])) continue;
                    
                    //     if (!trackMCgen->IsPhysicalPrimary()) continue;
                    Int_t allpdggen     =   trackMCallgen->GetPdgCode();
                    
                    // get the mother of the daughter particle
                    Int_t allmomlabgen                  = trackMCallgen->GetMother();
                    AliVParticle *allaodMothergen       = mcEvent->GetTrack(allmomlabgen);
                    if(allaodMothergen) allpdgMomgen    = allaodMothergen->PdgCode(); // pdg of mother particle
                    
                    //===================apply cuts on pt eta and centrality=====================
                    if ((fpTMCallgen>=momDownArray[imom])
                        &&(fpTMCallgen<momUpArray[imom])
                        &&(fCentrality>=centDownArray[icent])
                        &&(fCentrality<centUpArray[icent])){
                        
                        // calculate first moments
                        if (fChargeMCallgen < 0 || fChargeMCallgen > 0) allgenNch++;
                        if(fChargeMCallgen > 0) allgenPos++;
                        if(fChargeMCallgen < 0) allgenNeg++;
                        alltrCountMCgen++;
                    
                    } // end loop of apply cuts
                    
                } //============end of track loop====================
                
                fRunNumber = fInputEvent->GetRunNumber();
                
                // calculate second moments
                allgenMomentsCross  = allgenPos * allgenNeg;
                allgenMomentsPos    = allgenPos * allgenPos;
                allgenMomentsNeg    = allgenNeg * allgenNeg;
                
                // Fill the pt eta and centrality check in the tree
                if (alltrCountMCgen>0 ){
                    if(!fTreeSRedirector) return;
                    (*fTreeSRedirector)<<"MCallgen"<<
                    "pos="                 << allgenPos                    <<    // MC gen positive particles
                    "neg="                 << allgenNeg                    <<    // MC gen negative particles
                    "isample1="      	   << subsample1                   <<    // sample id for subsample method
                    "momentscross="        << allgenMomentsCross           <<    // MC gen second moments of pos and neg
                    "momentspos="          << allgenMomentsPos             <<    // MC gen second moments of pos
                    "momentsneg="          << allgenMomentsNeg             <<    // MC gen second moments of neg
                    "etaDown="             << etaDownArray[ieta]           <<    // lower edge of eta bin
                    "etaUp="               << etaUpArray[ieta]             <<    // upper edge of eta bin
                    "centDown="            << centDownArray[icent]         <<    // lower edge of cent bin
                    "centUp="              << centUpArray[icent]           <<    // upper edge of cent bin
                    "momDown="             << momDownArray[imom]           <<    // lower edge of mom bin
                    "momUp="               << momUpArray[imom]             <<    // upper edge of mom bin
                    "ieta="                << ieta                         <<    // eta index loop
                    "icent="               << icent                        <<    // cent index loop
                    "imom="                << imom                         <<    // mom index loop
                    "centBin="             << centBin                      <<    // centrality bin
                    "etabin="              << etaBin                       <<    // eta bin
                    "fRunNumber="          << fRunNumber                   <<    // Run Number
                    "\n";
                    
                } //===========tree filling========
                
            } //===========end of momentum loop
        } //===========end of eta loop
    } //===========end of centrality loop
    
    EventNumber++;
    
    PostData(1, fOutputList);
    PostData(2, fTree);
    PostData(3, fTreeMCrec);
    PostData(4, fTreeMCallrec);
    PostData(5, fTreeMCgen);
    PostData(6, fTreeMCallgen);
    
}
//--------------------------------------------------------------------------------

void AliAnalysisTaskEbyeNetChargeFluctuations::FillMCEffMatrix(){
    
    AliVEvent *event = InputEvent();
    if (!event) { Printf("UserExec: NO EVENT FOUND!"); return; }
    
    AliAODEvent* fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) { Printf("ERROR: fAOD not available"); return; }
    
    AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    //Physics selection
    UInt_t fSelectMask= eventHandler->IsEventSelected();
    Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7; // from https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsEventProp
    if(!isINT7selected)return ;
    
    AliAODVertex *fPrimaryVtx = (AliAODVertex*)InputEvent()->GetPrimaryVertex();
    if (!fPrimaryVtx){
        printf ("ERROR: no primary vertex\n");
        return;
    }
    
    Double_t xv=fPrimaryVtx->GetX();
    Double_t yv=fPrimaryVtx->GetY();
    Double_t zv=fPrimaryVtx->GetZ();

    //===========================Centrality calculation =====================
    fCentrality = -2;
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    fEventStatistics->Fill("found MultSelection object",1);
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
    
    if (TMath::Abs(zv) > 10.0) return;
    
    if(fCentrality < 0 || fCentrality >= 80) return;
    
    fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fArrayMC)
        AliFatal("No array of MC particles found !!!");
    
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
        AliError("ERROR: Could not retrieve MC event");
        return;
    }
    
    //-----------------------------------------------------------------------------------------
    // ----------------------------   reconstructed MC particles  ------------------------------
    // -----------------------------------------------------------------------------------------
    
    Int_t primPhysicalCount = 0;
    
    // Loop over tracks
    for(Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {   // track loop
        
        // get the AOD Tracks
        AliAODTrack *trackReal = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!trackReal) continue;
        
        // Initialize the variables
        pdgMom=0, pdgMomPhysicalPrim = 0, pdgMomPrim = 0, pdgPhysicalPrim = 0, pdgPrim = 0;
        
        if (!AcceptTrack(trackReal)) continue;

        // get track info
        Float_t fpTRec     = trackReal->Pt();
        Float_t fEtaRec    = trackReal->Eta();
        Float_t fPhiRec    = trackReal->Phi();
        Short_t fChargeRec = trackReal->Charge();
        Int_t   fPartID    = -10;
        
        // Eta cut
        if ((trackReal->Eta()<fEtaDown) || (trackReal->Eta()>fEtaUp)) continue;
        
        //MC tracks
        AliAODMCParticle* trackMCrec = (AliAODMCParticle*)fArrayMC->At(TMath::Abs(trackReal->GetLabel()));
        if(!trackMCrec) continue;
        
        if (!trackMCrec->IsPhysicalPrimary()) continue;
        
        Int_t fCentRec = fHistCentralityMultSelection->FindBin(fCentrality)-1; //new sk
    
        Int_t pdg = trackMCrec->GetPdgCode();
        
        // get the mother of the daughter particle
        Int_t momlab  = trackMCrec->GetMother();
        AliVParticle *aodMother = mcEvent->GetTrack(momlab);
        if(aodMother) pdgMom  = aodMother->PdgCode(); // pdg of mother particle
    
        // Select Physical primary particles and fill it on the tree
        if (trackMCrec->IsPhysicalPrimary()) {
            pdgMomPhysicalPrim = pdgMom;
            primPhysicalCount++;
            // calculate first moments
            if(fChargeRec > 0)  fPartID= 0;  // select positive
            if(fChargeRec < 0)  fPartID= 1;  // select negative
            
        }
        
        if (fPartID == -10) continue;
        
        Double_t xxxRec[5]={Float_t(fPartID),Float_t(fCentrality),fpTRec,fEtaRec,fPhiRec};
        
        if (fChargeRec>0){ fHistPosEffMatrixRec->Fill(xxxRec);}
        if (fChargeRec<0) {fHistNegEffMatrixRec->Fill(xxxRec);}
        
    } // ======= end of track loop =======
    
    // -----------------------------------------------------------------------------------------
    // ----------------------------   MC generated pure MC particles  --------------------------
    // -----------------------------------------------------------------------------------------
    
    Int_t primPhysicalCountgen = 0;
    // Loop over the MC gen tracks
    Int_t nMCgenTracks = mcEvent->GetNumberOfTracks();
    
    for(Int_t iPart = 0; iPart < nMCgenTracks; iPart++) {
        // Initialize the variables
        pdgMomPhysicalPrimgen = 0, pdgMomPrimgen = 0,pdgPhysicalPrimgen = 0, pdgPrimgen = 0;
        
        AliAODMCParticle *trackMCgen  = (AliAODMCParticle*)mcEvent->GetTrack(iPart);
        
        if(!trackMCgen) continue;
        
        // get track info
        Float_t fpTGen    = trackMCgen->Pt();
        Float_t fEtaGen   = trackMCgen->Eta();
        Float_t fPhiGen   = trackMCgen->Phi();
        Short_t fChargeGen= trackMCgen->Charge();
        Int_t fPartID  = -10;
        
        if (trackMCgen->Charge() == 0) continue;
        if (trackMCgen->Pt() < 0.2) continue;
        if ((trackMCgen->Eta()<fEtaDown) || (trackMCgen->Eta()>fEtaUp)) continue;
        
        if (!trackMCgen->IsPhysicalPrimary()) continue;
        
        Int_t fCentGen = fHistCentralityMultSelection->FindBin(fCentrality)-1; //new sk
        
        Int_t pdggen = trackMCgen->GetPdgCode();
        
        // get the mother of the daughter particle
        Int_t momlabgen  = trackMCgen->GetMother();
        AliVParticle *aodMothergen = mcEvent->GetTrack(momlabgen);
        if(aodMothergen) pdgMomgen = aodMothergen->PdgCode(); // pdg of mother particle
        
        // Select Physical primary particles and fill it on the tree
        if (trackMCgen->IsPhysicalPrimary()) {
            pdgMomPhysicalPrimgen = pdgMomgen;
            primPhysicalCountgen++;
            // calculate first moments
            
            if(fChargeGen > 0) fPartID= 0;  // select positive
            if(fChargeGen < 0) fPartID= 1;  // select negative;
        }
        
        if (fPartID == -10) continue;
        
        Double_t xxxGen[5]={Float_t(fPartID),Float_t(fCentrality),fpTGen,fEtaGen,fPhiGen};
        if (fChargeGen>0) fHistPosEffMatrixGen->Fill(xxxGen);
        if (fChargeGen<0) fHistNegEffMatrixGen->Fill(xxxGen);
        
    } // ======= end of track loop =======
    
}
//--------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskEbyeNetChargeFluctuations::AcceptTrack(AliAODTrack* aodtrack) const{
    
    if(!aodtrack) return kFALSE;
    Double_t pt = aodtrack->Pt();
    
    if(pt< 0.2) return kFALSE;
    if( aodtrack->Charge() == 0 ) return kFALSE;
    if(!aodtrack->TestFilterBit(768)) return kFALSE;   // for hybrid tracks
   // if(aodtrack->GetTPCNcls() < 80) return kFALSE;
   // if(aodtrack->Chi2perNDF() > 4.0) return kFALSE;
    
    fHistClustersTPC->Fill(aodtrack->GetTPCNcls());
    fHistChi2perNDF->Fill(aodtrack->Chi2perNDF());
    
    return kTRUE;
}

//------------------------------------------------------------------------------------------
void AliAnalysisTaskEbyeNetChargeFluctuations::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
    Info("AliAnalysisEbyeNetChargeFluctuations"," Task Successfully finished");
    AliInfo(Form("Found  %d MC events",EventNumber));
    
    
}
//____________________________________________________________________________________________
