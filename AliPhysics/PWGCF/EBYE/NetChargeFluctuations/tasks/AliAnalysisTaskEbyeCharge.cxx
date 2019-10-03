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
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliVEvent.h"
#include "AliAODHeader.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
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
#include "AliGenHijingEventHeader.h"
#include "AliGenEventHeader.h"
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliESDtrackCuts.h>

#include "AliAnalysisTaskEbyeCharge.h"

class AliAnalysisTaskEbyeCharge;

#include "TFile.h"
#include <iostream>
#include <fstream>

using namespace std;

ClassImp(AliAnalysisTaskEbyeCharge)

//================================================
//--------------Constructor-----------------------
//================================================

AliAnalysisTaskEbyeCharge::AliAnalysisTaskEbyeCharge()
: AliAnalysisTaskSE(),
fAOD(0x0),
fVxMax(3.),
fVyMax(3.),
fVzMax(10.),
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
fHistVx(0),
fHistVy(0),
fHistVz(0),
hTrackPt(0),
hTrackPhi(0),
hTrackEta(0),
hGenPt(0),
hGenPhi(0),
hGenEta(0),
fEventCuts(0),
fEta(0),
fpT(0),
fCharge(0),
fCentrality(0),
fhCent(0),
fPos(0),
fNeg(0),
fMomentsCross(0),
fMomentsPos(0),
fMomentsNeg(0),
fEtaMCgen(0),
fpTMCgen(0),
fEtaMC(0),
fpTMC(0),
genPos(0),
genNeg(0),
genMomentsCross(0),
genMomentsPos(0),
genMomentsNeg(0),
recPos(0),
recNeg(0),
recMomentsCross(0),
recMomentsPos(0),
recMomentsNeg(0),
fAnalysisType(0),
fTree(0x0),
fTreeTrackCuts(0x0),
fTreeMCTrackCuts(0x0),
fTreeMCrec(0x0),
fTreeMCgen(0x0)
{
    // default constructor
    
}

//-------------------------------------------------

AliAnalysisTaskEbyeCharge::AliAnalysisTaskEbyeCharge(const char* name)
:AliAnalysisTaskSE(name),
fAOD(0x0),
fVxMax(3.),
fVyMax(3.),
fVzMax(10.),
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
fHistVx(0),
fHistVy(0),
fHistVz(0),
hTrackPt(0),
hTrackPhi(0),
hTrackEta(0),
hGenPt(0),
hGenPhi(0),
hGenEta(0),
fEventCuts(0),
fEta(0),
fpT(0),
fCharge(0),
fCentrality(0),
fhCent(0),
fPos(0),
fNeg(0),
fMomentsCross(0),
fMomentsPos(0),
fMomentsNeg(0),
fEtaMCgen(0),
fpTMCgen(0),
fEtaMC(0),
fpTMC(0),
genPos(0),
genNeg(0),
genMomentsCross(0),
genMomentsPos(0),
genMomentsNeg(0),
recPos(0),
recNeg(0),
recMomentsCross(0),
recMomentsPos(0),
recMomentsNeg(0),
fAnalysisType(0),
fTree(0x0),
fTreeTrackCuts(0x0),
fTreeMCTrackCuts(0x0),
fTreeMCrec(0x0),
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
    
    
}
//==========================================================
//--------------Destructor----------------------------------
//==========================================================
AliAnalysisTaskEbyeCharge::~AliAnalysisTaskEbyeCharge()
{
    //destructor
    if(fOutputList) 			        delete fOutputList;
    if (fHistPosEffMatrixRec) 		    delete fHistPosEffMatrixRec;
    if (fHistNegEffMatrixRec) 		    delete fHistNegEffMatrixRec;
    if (fHistPosEffMatrixGen) 		    delete fHistPosEffMatrixGen;
    if (fHistNegEffMatrixGen) 		    delete fHistNegEffMatrixGen;
    if (fhCent)               		    delete fhCent;
    if (fHistCentralityMultSelection)   delete fHistCentralityMultSelection;
}
//============================================================
//--------------UserCreateOutputObjects-----------------------
//============================================================
void AliAnalysisTaskEbyeCharge::UserCreateOutputObjects()
{
    
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    // Event statistics
    fEventStatistics = new TH1I("fEventStatistics","",10,0,10);
   // fEventStatistics->SetBit(TH1::kCanRebin);
    fOutputList->Add(fEventStatistics);
    
    // ****************** Efficiency matrix histograms ************************
    
    const Int_t ndim=5;
    // const Int_t nEtaBins=(fEtaUp-fEtaDown)*20;
    Int_t nbins0[ndim]  ={2, 9, 50, 16 , 40  };
    Double_t xmin0[ndim]={0, 0, 0.2, -0.8, 0.  };
    Double_t xmax0[ndim]={2, 9, 3.2,  0.8, 6.25};
    fHistPosEffMatrixRec  =new THnF("fHistPosEffMatrixRec","fHistPosEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixRec  =new THnF("fHistNegEffMatrixRec","fHistNegEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistPosEffMatrixGen  =new THnF("fHistPosEffMatrixGen","fHistPosEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixGen  =new THnF("fHistNegEffMatrixGen","fHistNegEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    TString axisNameEff[ndim]  = {"charge"      ,"Centrality"     ,"momentum"      ,"eta"  ,"phi"};
    TString axisTitleEff[ndim] = {"charge type" ,"Centrality (%)" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta" ,"#phi"};
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
    hTrackPt = new TH1F("hTrackPt","track p_{T};p_{T} (GeV/c);",100,0,10);
    fOutputList->Add(hTrackPt);
    hTrackPhi = new TH1F("hTrackPhi","track #varphi;#varphi;",100,0,2*TMath::Pi());
    fOutputList->Add(hTrackPhi);
    hTrackEta = new TH1F("hTrackEta","track #eta;#eta;",100,-0.8,0.8);
    fOutputList->Add(hTrackEta);
    
    hGenPt = new TH1F("hGenPt","generated p_{T};p_{T} (GeV/c);",100,0,10);
    fOutputList->Add(hGenPt);
    hGenPhi = new TH1F("hGenPhi","generated #varphi;#varphi;",100,0,2*TMath::Pi());
    fOutputList->Add(hGenPhi);
    hGenEta = new TH1F("hGenEta","generated #eta;#eta;",100,-0.8,0.8);
    fOutputList->Add(hGenEta);
    
    // AliVEvent QA Cuts
    fEventCuts.AddQAplotsToList(fOutputList);
    
    fHistCentralityMultSelection = new TH1D("fHistCentralityMultSelection","Centrality Percentile ;Centrality;Entries",10,0.,100.);  //new sk
    fHistCentralityMultSelection->GetXaxis()->SetTitle("Centrality (%)");
    fOutputList->Add(fHistCentralityMultSelection);
    
    fHistZVertexCent = new TH2F("fHistZVertexCent"," Vz of primary vertex Vs Centrality; V_{Z} {cm} ; Centrality {%}", 60, -15, 15,100,0,100);
    fOutputList->Add(fHistZVertexCent);
    // ************************************************************************
    
    // Tree for pt eta and centrality checks
    fTreeSRedirector = new TTreeSRedirector();
    fTree       = ((*fTreeSRedirector)<<"Realdata").GetTree();
    fTreeMCrec  = ((*fTreeSRedirector)<<"MCrec").GetTree();
    fTreeMCgen  = ((*fTreeSRedirector)<<"MCgen").GetTree();
    
    cout << "============================================================" << endl;
    cout << "=================End of User Create Object=================="  << endl;
    cout << "============================================================" << endl;
    
    //Send data to the container
    PostData(1, fOutputList);
    PostData(2, fTree);
    PostData(3, fTreeMCrec);
    PostData(4, fTreeMCgen);
    
}
//============================================================
//----------------------UserExec------------------------------
//============================================================
void AliAnalysisTaskEbyeCharge::UserExec(Option_t *){
    
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

void AliAnalysisTaskEbyeCharge::doAODEvent(){
    
    fEventStatistics->Fill("before cuts",1);
    
    //  AliVEvent *event = fInputEvent();
    if (!fInputEvent) return;
    
    fEventStatistics->Fill("after event check",1);
    
    fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
    
    if(!fAOD) { Printf("ERROR: fAOD not available"); return; }
    fEventStatistics->Fill("after aod check",1);
    
    //===========================Centrality calculation =====================
    fCentrality = -2;
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    fEventStatistics->Fill("found MultSelection object",1);
    
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
    
    fEventStatistics->Fill("found primary vertex",1);
    
    if (TMath::Abs(zv) > 10.0) return;
    
    fEventStatistics->Fill("vz cut",1);
    
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
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
    Double_t centDownArray[9]= {0., 5.,  10., 20., 30., 40., 50., 60., 70.};
    Double_t centUpArray[9]  = {5., 10., 20., 30., 40., 50., 60., 70., 80.};
    for (Int_t imom=0; imom<4; imom++){
        for (Int_t ieta=0; ieta<8; ieta++){
            for (Int_t icent=0; icent<9; icent++){
                
                Double_t centBin = (centDownArray[icent]+centUpArray[icent])/2.;
                
                // Initialize the positive and negative particles
                fPos = 0;
                fNeg = 0;
                Int_t trCount=0;
                
                //=============================== track loop ==============================
                Int_t iTracks(fAOD->GetNumberOfTracks());
                for(Int_t i(0); i < iTracks; i++) {
                    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
                    
                    if(!track) {
                        AliError(Form("ERROR: Could not retrieve AODtrack %d",i));
                        continue;
                    }
                    
                    if(!AcceptTrack(track)) continue;
                    
                    fpT     = track->Pt();
                    fEta    = track->Eta();
                    fPhi    = track->Phi();
                    fCharge = track->Charge();
                    
                    hTrackPt->Fill(fpT);
                    hTrackPhi->Fill(fPhi);
                    hTrackEta->Fill(fEta);
                    
                    if ((fEta<etaDownArray[ieta]) || (fEta>etaUpArray[ieta])) continue;  // eta Cut
                    
                    //===================apply cuts on pt eta and centrality=====================
                    if ((fpT>=momDownArray[imom])
                        &&(fpT<momUpArray[imom])
                        &&(fCentrality>=centDownArray[icent])
                        &&(fCentrality<centUpArray[icent])){
                        
                        // calculate first moments
                        if(fCharge > 0) fPos++;
                        if(fCharge < 0) fNeg++;
                        trCount++;
                        
                    }
                } //============end of track loop
                
                // calculate second moments
                fMomentsCross =  fPos * fNeg;
                fMomentsPos    = fPos * fPos;
                fMomentsNeg    = fNeg * fNeg;
                
                // Fill the pt eta and centrality check in the tree
                if ( trCount>0 ){
                    if(!fTreeSRedirector) return;
                    (*fTreeSRedirector)<<"Realdata"<<
                    // upper edge of momentum bin
                    "etaDown="      << etaDownArray[ieta]  <<       // lower edge of eta bin
                    "etaUp="        << etaUpArray[ieta]    <<       // upper edge of eta bin
                    "centDown="     << centDownArray[icent]<<       // lower edge of cent bin
                    "centUp="       << centUpArray[icent]  <<       // upper edge of cent bin
                    "pos="          << fPos                <<       // positive particles
                    "neg="          << fNeg                <<       // negative particles
                    "momentscross=" << fMomentsCross       <<       // second moments
                    "momentspos="   << fMomentsPos         <<
                    "momentsneg="   << fMomentsNeg         <<
                    "ieta="         << ieta                <<
                    "icent="        << icent               <<
                    "imom="         << imom                <<
                    "centBin="      << centBin             <<      // centrality bining
                    "\n";
                } //===========tree filling========
                
            } //===========end of momentum loop
        } //===========end of eta loop
    } //===========end of centralityy loop
    
    PostData(1, fOutputList);
    PostData(2, fTree);
    PostData(3, fTreeMCrec);
    PostData(4, fTreeMCgen);
    
}

void AliAnalysisTaskEbyeCharge::doMCAODEvent(){
    
    fEventStatistics->Fill("before cuts",1);
    
    //  AliVEvent *event = fInputEvent();
    if (!fInputEvent) return;
    
    fEventStatistics->Fill("after event check",1);
    
    fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
    
    if(!fAOD) { Printf("ERROR: fAOD not available"); return; }
    fEventStatistics->Fill("after aod check",1);
    
    //===========================Centrality calculation =====================
    fCentrality = -2;
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    fEventStatistics->Fill("found MultSelection object",1);
    
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
    
    fEventStatistics->Fill("found primary vertex",1);
    
    if (TMath::Abs(zv) > 10.0) return;
    
    fEventStatistics->Fill("vz cut",1);
    
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    if(fCentrality < 0 || fCentrality >= 80) return;
    fHistCentralityMultSelection->Fill(fCentrality);
    fEventStatistics->Fill("centrality selection",1);
    
    fHistZVertexCent->Fill(zv, fCentrality);
    
    if (!fEventCuts.AcceptEvent(fInputEvent)) return;
    fEventStatistics->Fill("AliEventCuts",1);
    
    // =========Different Acceptance cuts for Eta momentum and Centrality=============
    Double_t momDownArray[4]  = {0.2, 0.2, 0.6, 0.6};
    Double_t momUpArray[4]   = {2.0, 5.0, 2.0, 5.0};
    Double_t etaDownArray[8] = {-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8};
    Double_t etaUpArray[8]   = { 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8};
    Double_t centDownArray[9]= {0., 5.,  10., 20., 30., 40., 50., 60., 70.};
    Double_t centUpArray[9]  = {5., 10., 20., 30., 40., 50., 60., 70., 80.};
    
    
    for (Int_t imom=0; imom<4; imom++){
        for (Int_t ieta=0; ieta<8; ieta++){
            for (Int_t icent=0; icent<9; icent++){
                
                // -----------------------------------------------------------------------------------------
                // cout<< " ------------------reconstructed MC particles------------------------" << endl;
                // -----------------------------------------------------------------------------------------
                
                Double_t centBin = (centDownArray[icent]+centUpArray[icent])/2.;
                
                //cout << "Info::Centbin================" << centBin << endl;
                // Initialize
                Int_t trCountMCrec=0;
                Int_t recPos = 0, recNeg = 0;
                Int_t momlab = 0;
                Int_t primCount = 0, primPhysicalCount = 0;
                
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
                    
                    AliAODMCParticle* trackMCrec = (AliAODMCParticle*)fMCEvent->GetTrack(i);
                    if(!trackMCrec) continue;
                    if (!trackMCrec->IsPhysicalPrimary()) continue;
                    
                    pdg = trackMCrec->GetPdgCode(); // pdg of the primary particle
                    
                    // get the mother of the daughter particle
                    momlab  = trackMCrec->GetMother();
                    AliVParticle *aodMother = fMCEvent->GetTrack(momlab);
                    if(aodMother) pdgMom  = aodMother->PdgCode(); // pdg of mother particl
                    
                    fpT     = trackReal->Pt();
                    fEta    = trackReal->Eta();
                    fPhi    = trackReal->Phi();
                    fCharge = trackReal->Charge();
                    
                    hTrackPt->Fill(fpT);
                    hTrackPhi->Fill(fPhi);
                    hTrackEta->Fill(fEta);
                    
                    // MC track cuts
                    if ((fEta<etaDownArray[ieta]) || (fEta>etaUpArray[ieta])) continue;  // eta Cut
                    
                    //===================apply cuts on pt eta and centrality=====================
                    if ((fpT>=momDownArray[imom])
                        &&(fpT<momUpArray[imom])
                        &&(fCentrality>=centDownArray[icent])
                        &&(fCentrality<centUpArray[icent])){
                        
                        // Select Physical primary particles and fill it on the tree
                        if (trackMCrec->IsPhysicalPrimary()) {
                            pdgMomPhysicalPrim = pdgMom;
                            primPhysicalCount++;
                            // calculate first moments
                            if(fCharge > 0) recPos++;
                            if(fCharge < 0) recNeg++;
                            trCountMCrec++;
                            
                        }
                        // Select Primary particles and fill it on the tree
                        if (trackMCrec->IsPrimary()) {
                            pdgMomPrim = pdgMom;
                            primCount++;
                        }
                    } // end loop of apply cuts
                    
                } //============end of track loop====================
                
                // calculate second moments
                recMomentsCross  = recPos * recNeg;
                recMomentsPos    = recPos * recPos;
                recMomentsNeg    = recNeg * recNeg;
                
                // Tree for all the cuts variables
                if ( trCountMCrec>0 ){
                    if(!fTreeSRedirector) return;
                    (*fTreeSRedirector)<<"MCrec"<<
                    "Physcount="           << primPhysicalCount       <<    // no. of Isphysical primary particles
                    "Prmcount="            << primCount               <<    // no. of Isprimary particles
                    "pos="                 << recPos                  <<    // MC gen positive particles
                    "neg="                 << recNeg                  <<    // MC gen negative particles
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
                    "\n";
                    
                } // trackcuts tree fill========
                
                // -----------------------------------------------------------------------------------------
                // cout<< " ------------------Generated MC particles------------------------" << endl;
                // -----------------------------------------------------------------------------------------
                // Initialize the positive and negative particles
                genPos = 0, genNeg = 0;
                Int_t trCountMCgen = 0;
                Int_t momlabgen = 0;
                Int_t primCountgen = 0, primPhysicalCountgen = 0;
                
                //=============================== track loop ==============================
                
                // Loop over the MC gen tracks
                Int_t nMCgenStackTracks = fMCEvent->GetNumberOfTracks();
                for(int iPart = 0; iPart < nMCgenStackTracks; iPart++) {
                    
                    // Initialize the variables
                    pdgMomgen = 0, pdgMomPhysicalPrimgen = 0, pdgMomPrimgen = 0,pdgPhysicalPrimgen = 0, pdgPrimgen = 0;
                    
                    AliAODMCParticle *trackMCgen  = (AliAODMCParticle*)fMCEvent->GetTrack(iPart);
                    
                    if(!trackMCgen) continue;
                    if (!trackMCgen->IsPhysicalPrimary()) continue;
                    pdggen = trackMCgen->GetPdgCode(); // pdg of the primary particle
                    
                    // get the mother of the daughter particle
                    momlabgen  = trackMCgen->GetMother();
                    AliVParticle *aodMothergen = fMCEvent->GetTrack(momlabgen);
                    if(aodMothergen) pdgMomgen = aodMothergen->PdgCode(); // pdg of mother particle
                    
                    fpTMCgen  = trackMCgen->Pt();
                    fPhiMCgen = trackMCgen->Pt();
                    fCharge   = trackMCgen->Charge();
                    fEtaMCgen = trackMCgen->Eta();
                    
                    hGenPt->Fill(fpTMCgen);
                    hGenPhi->Fill(fPhiMCgen);
                    hGenEta->Fill(fEtaMCgen);
                    
                    // apply primary vertex and eta cut
                    if ((fEtaMCgen < etaDownArray[ieta]) || (fEtaMCgen > etaUpArray[ieta])) continue;
                    
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
                            if(fCharge > 0) genPos++;
                            if(fCharge < 0) genNeg++;
                            trCountMCgen++;
                            
                        }
                        // Select Primary particles and fill it on the tree
                        if (trackMCgen->IsPrimary()) {
                            pdgMomPrimgen = pdgMomgen;
                            primCountgen++;
                        }
                    } // end loop of apply cuts
                    
                    
                } //============end of track loop====================
                
                fEventStatistics->Fill("after event loop",1);
                
                // calculate second moments
                genMomentsCross  = genPos * genNeg;
                genMomentsPos    = genPos * genPos;
                genMomentsNeg    = genNeg * genNeg;
                
                
                // Fill the pt eta and centrality check in the tree
                if ( trCountMCgen>0 ){
                    if(!fTreeSRedirector) return;
                    (*fTreeSRedirector)<<"MCgen"<<
                    "Physcount="           << primPhysicalCountgen      <<    // no. of Isphysical primary particles
                    "Prmcount="            << primCountgen              <<    // no. of Isprimary particles
                    "pos="                 << genPos                    <<    // MC gen positive particles
                    "neg="                 << genNeg                    <<    // MC gen negative particles
                    "momentscross="        << genMomentsCross           <<    // MC gen second moments of pos and neg
                    "momentspos="          << genMomentsPos             <<    // MC gen second moments of pos
                    "momentsneg="          << genMomentsNeg             <<    // MC gen second moments of neg
                    "etaDown="             << etaDownArray[ieta]       <<    // lower edge of eta bin
                    "etaUp="               << etaUpArray[ieta]         <<    // upper edge of eta bin
                    "centDown="            << centDownArray[icent]     <<    // lower edge of cent bin
                    "centUp="              << centUpArray[icent]       <<    // upper edge of cent bin
                    "momDown="             << momDownArray[imom]       <<    // lower edge of mom bin
                    "momUp="               << momUpArray[imom]         <<    // upper edge of mom bin
                    "ieta="                << ieta                     <<    // eta index loop
                    "icent="               << icent                    <<    // cent index loop
                    "imom="                << imom                     <<    // mom index loop
                    "\n";
                    
                } //===========tree filling========
                
            } //===========end of momentum loop
        } //===========end of eta loop
    } //===========end of centrality loop
    
    
    PostData(1, fOutputList);
    PostData(2, fTree);
    PostData(3, fTreeMCrec);
    PostData(4, fTreeMCgen);
    
}
//--------------------------------------------------------------------------------

void AliAnalysisTaskEbyeCharge::FillMCEffMatrix(){
    
    //  AliVEvent *event = fInputEvent();
    if (!fInputEvent) return;
    
    fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
    
    if(!fAOD) { Printf("ERROR: fAOD not available"); return; }
    
    //===========================Centrality calculation =====================
    fCentrality = -2;
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    
    //Physics selection
    UInt_t fSelectMask= fInputHandler->IsEventSelected();
    Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7; // from https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsEventProp
    if(!isINT7selected)return ;
    
    AliAODVertex *fPrimaryVtx = (AliAODVertex*)fInputEvent->GetPrimaryVertex();
    if (!fPrimaryVtx){
        printf ("ERROR: no primary vertex\n");
        return;
    }
    
    Double_t xv=fPrimaryVtx->GetX();
    Double_t yv=fPrimaryVtx->GetY();
    Double_t zv=fPrimaryVtx->GetZ();
    
    if (TMath::Abs(zv) > 10.0) return;
    
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    if(fCentrality < 0 || fCentrality >= 80) return;
    
    if (!fEventCuts.AcceptEvent(fInputEvent)) return;
    
    //-----------------------------------------------------------------------------------------
    // ----------------------------   reconstructed MC particles  ------------------------------
    // -----------------------------------------------------------------------------------------
    
    // Loop over tracks
    for(Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {   // track loop
        
        // get the AOD Tracks
        AliAODTrack *trackReal = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!trackReal) continue;
        //    fEtaMC    = trackReal->Eta();
        if ((trackReal->Eta()<fEtaDown) || (trackReal->Eta()>fEtaUp)) continue;
        if(!AcceptTrack(trackReal)) continue;            // real track cuts
        //MC tracks
        AliAODMCParticle* trackMCrec = (AliAODMCParticle*)fMCEvent->GetTrack(i);
        if(!trackMCrec) continue;
        if (!trackMCrec->IsPhysicalPrimary()) continue;  // MC primary check
        
        // get track info
        Float_t fpTRec     = trackReal->Pt();
        Float_t fYRec      = trackReal->Y();
        Float_t fEtaRec    = trackReal->Eta();
        Float_t fPhiRec    = trackReal->Phi();
        Short_t fChargeRec = trackReal->Charge();
        
        //Int_t fCentRec = fhCent->FindBin(fCentrality)-1;
        Int_t fCentRec = fHistCentralityMultSelection->FindBin(fCentrality)-1; //new sk
        //Int_t fCentRec = fCentrality;
        Int_t fPartID  = -10;
        
        if(fChargeRec > 0) recPos++;
        if(fChargeRec < 0) recNeg++;
        
        if(fChargeRec > 0)  fPartID=0;  // select positive
        if(fChargeRec < 0)  fPartID=1;  // select negative
        if (fPartID == -10) continue;
        
        Double_t xxxRec[5]={Float_t(fPartID),Float_t(fCentRec),fpTRec,fEtaRec,fPhiRec};
        
        if (fChargeRec>0){ fHistPosEffMatrixRec->Fill(xxxRec);}
        if (fChargeRec<0) {fHistNegEffMatrixRec->Fill(xxxRec);}
        
    } // ======= end of track loop =======
    
    // -----------------------------------------------------------------------------------------
    // ----------------------------   MC generated pure MC particles  --------------------------
    // -----------------------------------------------------------------------------------------
    
    // Loop over the MC gen tracks
    for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++){ // track loop
        
        AliAODMCParticle *trackMCgen  = (AliAODMCParticle*)fMCEvent->GetTrack(iTrack);
        if(!trackMCgen) continue;
        if ((trackMCgen->Eta()<fEtaDown) || (trackMCgen->Eta()>fEtaUp)) continue;
        
        if (!trackMCgen->IsPhysicalPrimary()) continue;
        
        // get track info
        Float_t fpTGen    = trackMCgen->Pt();
        Float_t fYGen     = trackMCgen->Y();
        Float_t fEtaGen   = trackMCgen->Eta();
        Float_t fPhiGen   = trackMCgen->Phi();
        Short_t fChargeGen= trackMCgen->Charge();
        
        //Int_t fCentGen = fhCent->FindBin(fCentrality)-1;
        Int_t fCentGen = fHistCentralityMultSelection->FindBin(fCentrality)-1; //new sk
        // Int_t fCentGen = fCentrality;
        Int_t fPartID  = -10;
        
        // Efficiency matices for positive and negative particles
        if(fChargeGen > 0) genPos++;
        if(fChargeGen < 0) genNeg++;
        
        if(fChargeGen > 0)  fPartID=0;  // select positive
        if(fChargeGen < 0)  fPartID=1;  // select negative
        if (fPartID == -10) continue;
        
        Double_t xxxGen[5]={Float_t(fPartID),Float_t(fCentGen),fpTGen,fEtaGen,fPhiGen};
        if (fChargeGen>0) fHistPosEffMatrixGen->Fill(xxxGen);
        if (fChargeGen<0) fHistNegEffMatrixGen->Fill(xxxGen);
        
    } // ======= end of track loop =======
    
}
//---------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskEbyeCharge::AcceptTrack(AliAODTrack* aodtrack) const{
    
    if(!aodtrack) return kFALSE;
    if( aodtrack->Charge() == 0 ) return kFALSE;
    
    if(!aodtrack->TestFilterBit(768)) return kFALSE;
    
    return kTRUE;
}
//------------------------------------------------------------------------------------------
void AliAnalysisTaskEbyeCharge::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
    Info("AliAnalysisEbyECharge"," Task Successfully finished");
    
}
//____________________________________________________________________________________________
