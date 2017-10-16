//
//  AliAnalysisTaskTPCCalBeauty.cxx
//  
//
//  Created by Erin Gauger
//
//

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskTPCCalBeauty.h"
#include "AliKFParticle.h"
#include "AliAODMCParticle.h"
#include "AliGenHijingEventHeader.h"

//#include "AliCentralitySelectionTask.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"

class AliAnalysisTaskTPCCalBeauty;

using namespace std;

ClassImp(AliAnalysisTaskTPCCalBeauty)

AliAnalysisTaskTPCCalBeauty::AliAnalysisTaskTPCCalBeauty() :
    AliAnalysisTaskSE(),
    fAOD(0),
    fMCHeader(0),
    fMCarray(0),
    fMCparticle(0),
    fpidResponse(0),
    fOutputList(0),
    fMultSelection(0),
    fCentrality(-1),
    fCentralityMin(0),
    fCentralityMax(100),
    fEMCEG1(kFALSE),
    fDCalDG1(kFALSE),
    fFlagClsTypeEMC(kTRUE),
    fFlagClsTypeDCAL(kTRUE),
    fUseTender(kTRUE),
    fFlagULS(kFALSE),
    fFlagLS(kFALSE),
    fNevents(0),
    fVtX(0),
    fVtY(0),
    fVtZ(0),
    fTrkPtB4TC(0),
    fTrkPt(0),
    fTrkP(0),
    fTrkClsPhi(0),
    fTrkClsEta(0),
    fClsPhi(0),
    fClsEta(0),
    fClsEamDCal(0),
    fClsEamEMCal(0),
    fTrkPhi(0),
    fTrkEta(0),
    fdEdx(0),
    fCentCheck(0),
    fTrigCheck(0),
    fEMCTrkMatch(0),
    fInvmassLS(0),
    fInvmassULS(0),
    fPhotonicElecYield(0),
    fULSdcaBelow(0),
    fLSdcaBelow(0),
    fPhotonicDCA(0),
    fDalitzDCA(0),
    fInclElecDCA(0),
    fInclElecEoP(0),
    fHadronEoP(0),
    fHadronDCA(0),
    fPi0Weight(0),
    fEtaWeight(0),
    //fPi0EtaWeight(0),
    fWeight(1),
    fPi0DCA(0),
    fEtaDCA(0),
    //fPi0EtaDCA(0),
    fPhotonMCTrueDCA(0),
    fPi0MCTrueWeightEnhDCA(0),
    fEtaMCTrueWeightEnhDCA(0),
    //fPi0EtaMCTrueWeightEnhDCA(0),
    fPi0MCTrueHijingDCA(0),
    fEtaMCTrueHijingDCA(0),
    //fPi0EtaMCTrueHijingDCA(0),
    //fEtaPt_Hijing(0),
    //fPionPt_Hijing(0),
    //fEtaPt_Enh(0),
    //fPionPt_Enh(0),
    fNtotMCpart(0),
    fNpureMC(0),
    fNembMCpi0(0),
    fNembMCeta(0),
    fSprsPi0EtaWeightCal(0),
    fSprsTemplates(0)

    //fElectronSprs(0)
{
    //Root IO constructor, don't allocate memory here
}
//_____________________________________________________________________
AliAnalysisTaskTPCCalBeauty::AliAnalysisTaskTPCCalBeauty(const char *name) :
    AliAnalysisTaskSE(name),
    fAOD(0),
    fMCHeader(0),
    fMCarray(0),
    fMCparticle(0),
    fpidResponse(0),
    fOutputList(0),
    fMultSelection(0),
    fCentrality(-1),
    fCentralityMin(0),
    fCentralityMax(100),
    fEMCEG1(kFALSE),
    fDCalDG1(kFALSE),
    fFlagClsTypeEMC(kTRUE),
    fFlagClsTypeDCAL(kTRUE),
    fUseTender(kTRUE),
    fFlagULS(kFALSE),
    fFlagLS(kFALSE),
    fNevents(0),
    fVtX(0),
    fVtY(0),
    fVtZ(0),
    fTrkPtB4TC(0),
    fTrkPt(0),
    fTrkP(0),
    fTrkClsPhi(0),
    fTrkClsEta(0),
    fClsPhi(0),
    fClsEta(0),
    fClsEamDCal(0),
    fClsEamEMCal(0),
    fTrkPhi(0),
    fTrkEta(0),
    fdEdx(0),
    fCentCheck(0),
    fTrigCheck(0),
    fEMCTrkMatch(0),
    fInvmassLS(0),
    fInvmassULS(0),
    fPhotonicElecYield(0),
    fULSdcaBelow(0),
    fLSdcaBelow(0),
    fPhotonicDCA(0),
    fDalitzDCA(0),
    fInclElecDCA(0),
    fInclElecEoP(0),
    fHadronEoP(0),
    fHadronDCA(0),
    fPi0Weight(0),
    fEtaWeight(0),
    //fPi0EtaWeight(0),
    fPi0DCA(0),
    fEtaDCA(0),
    //fPi0EtaDCA(0),
    fWeight(1),
    fPhotonMCTrueDCA(0),
    fPi0MCTrueWeightEnhDCA(0),
    fEtaMCTrueWeightEnhDCA(0),
    //fPi0EtaMCTrueWeightEnhDCA(0),
    fPi0MCTrueHijingDCA(0),
    fEtaMCTrueHijingDCA(0),
    //fPi0EtaMCTrueHijingDCA(0),
    //fEtaPt_Hijing(0),
    //fPionPt_Hijing(0),
    //fEtaPt_Enh(0),
    //fPionPt_Enh(0),
    fNtotMCpart(0),
    fNpureMC(0),
    fNembMCpi0(0),
    fNembMCeta(0),
    fSprsPi0EtaWeightCal(0),
    fSprsTemplates(0)
    //fElectronSprs(0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskTPCCalBeauty::~AliAnalysisTaskTPCCalBeauty()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;
    }
}
//_____________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::UserCreateOutputObjects()
{
    //Weights for pho reco?
    
    /////////////////
    // Output List //
    /////////////////
    
    //create a new TList that owns its objects
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    //create our histos and add them to the list
    fNevents = new TH1F("fNevents", "No. of Events; Counts", 3,-0.5,2.5);
    fOutputList->Add(fNevents);
    fNevents->GetXaxis()->SetBinLabel(1,"All");
    fNevents->GetXaxis()->SetBinLabel(2,">2 Trks");
    fNevents->GetXaxis()->SetBinLabel(3,">2 Trks, Vtx_{z}<10cm");
    
    fVtX = new TH1F("fVtX","X Vertex Position;Vtx_{X};Counts",50,-5,5);
    fOutputList->Add(fVtX);
    
    fVtY = new TH1F("fVtY","Y Vertex Position;Vtx_{Y};Counts",50,-5,5);
    fOutputList->Add(fVtY);
    
    fVtZ = new TH1F("fVtZ","Z Vertex Position;Vtx_{Z};Counts",100,-10,10);
    fOutputList->Add(fVtZ);
    
    fTrkPtB4TC = new TH1F("fTrkPtB4TC","Track p_{T} Distribution before track cuts;p_{T} (GeV/c);Counts",100,0,50);
    fOutputList->Add(fTrkPtB4TC);
    
    fTrkPt = new TH1F("fTrkPt","Track p_{T} Distribution;p_{T} (GeV/c);Counts",100,0,50);
    fOutputList->Add(fTrkPt);
    
    fTrkP = new TH1F("fTrkP","Track p Distribution;p (GeV/c);Counts",100,0,50);
    fOutputList->Add(fTrkP);
    
    fTrkClsPhi = new TH1F("fTrkClsPhi","Track and Cluster #Delta #phi Distribution;#Delta #phi;Counts",100,0,6.3);
    fOutputList->Add(fTrkClsPhi);
    
    fTrkClsEta = new TH1F("fTrkClsEta","Track and Cluster #Delta #eta Distribution;#Delta #eta;Counts",100,-1.5,1.5);
    fOutputList->Add(fTrkClsEta);
    
    fClsPhi = new TH1F("fClsPhi","Cluster #phi Distribution;#phi;Counts",100,0,6.3);
    fOutputList->Add(fClsPhi);
    
    fClsEta = new TH1F("fClsEta","Cluster #eta Distribution;#eta;Counts",100,-1.5,1.5);
    fOutputList->Add(fClsEta);
    
    fClsEamDCal = new TH1F("fClsEamDCal","Cluster Energy after track matching to DCal;Cluster E;Counts",250,0.,50);
    fOutputList->Add(fClsEamDCal);

    fClsEamEMCal = new TH1F("fClsEamEMCal","Cluster Energy after track matching to EMCal;Cluster E;Counts",250,0.,50.);
    fOutputList->Add(fClsEamEMCal);
    
    fTrkPhi = new TH1F("fTrkPhi","Track #phi Distribution after matching;#phi;Counts",100,0,6.3);
    fOutputList->Add(fTrkPhi);
    
    fTrkEta = new TH1F("fTrkEta","Track #eta Distribution after matching;#eta;Counts",100,-1.5,1.5);
    fOutputList->Add(fTrkEta);
    
    fdEdx = new TH1F("fdEdx","Track dE/dx Distribution;dE/dx;Counts",160,0,160);
    fOutputList->Add(fdEdx);
    
    fCentCheck = new TH1F("fCentCheck","Event Centrality Distribution;Centrality;Counts",100,0,100);
    fOutputList->Add(fCentCheck);
    
    fTrigCheck = new TH1F("fTrigCheck", "No. of Events; Counts",3,-0.5,2.5);
    fOutputList->Add(fTrigCheck);
    fTrigCheck->GetXaxis()->SetBinLabel(1,"INT7");
    fTrigCheck->GetXaxis()->SetBinLabel(2,"EG1");
    fTrigCheck->GetXaxis()->SetBinLabel(3,"DGl");
    
    fEMCTrkMatch = new TH2F("fEMCTrkMatch","EMCal cluster distance from closest track", 100, -0.3, 0.3,100,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch);
    
    fInvmassLS = new TH1F("fInvmassLS", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassLS);
    
    fInvmassULS = new TH1F("fInvmassULS", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassULS);
    
    fPhotonicElecYield = new TH1F("fPhotonicElecYield","Photonic Elec Yield; p_{T}(GeV/c); yield;", 60,0,30.);
    fOutputList->Add(fPhotonicElecYield);
    
    fULSdcaBelow = new TH2F("fULSdcaBelow","ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fULSdcaBelow);
    
    fLSdcaBelow = new TH2F("fLSdcaBelow","LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fLSdcaBelow);
    
    fPhotonicDCA = new TH2F("fPhotonicDCA","Photonic DCA using MC PID; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPhotonicDCA);
    
    fDalitzDCA = new TH2F("fDalitzDCA","Pion and eta DCA using MC PID; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fDalitzDCA);
    
    fInclElecDCA = new TH2F("fInclElecDCA","Incl Elec DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fInclElecDCA);
    
    fInclElecEoP = new TH2F("fInclElecEoP","Incl Elec E/p; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fInclElecEoP);
    
    fHadronEoP = new TH2F("fHadronEoP","Unscaled Hadron E/p; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fHadronEoP);
    
    fHadronDCA = new TH2F("fHadronDCA","Unscaled Hadron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fHadronDCA);
    
    fPi0Weight = new TF1("fPi0Weight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    fEtaWeight = new TF1("fEtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    //fPi0EtaWeight = new TF1("fPi0EtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    fPi0Weight->SetParameters(3.10605e+04,-2.88504e-01,4.90349e-03,4.82424e-01,3.88146e+00);
    fEtaWeight->SetParameters(1.48495e+03,-2.20500e-01,3.35819e-03,1.14265e+00,4.06488e+00);
    //fPi0EtaWeight->SetParameters(3.68627e+03,-2.15364e-01,3.26753e-03,1.05213e+00,4.20391e+00);
    
    fPi0DCA = new TH2F("fPi0DCA","Pi0 DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPi0DCA);
    
    fEtaDCA = new TH2F("fEtaDCA","Eta DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fEtaDCA);
    
    //fPi0EtaDCA = new TH2F("fPi0EtaDCA","Pi0+Eta DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    //fOutputList->Add(fPi0EtaDCA);
    
    fPhotonMCTrueDCA = new TH2F("fPhotonMCTrueDCA","MC True gamma->e DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPhotonMCTrueDCA);
    
    fPi0MCTrueWeightEnhDCA = new TH2F("fPi0MCTrueWeightEnhDCA","MC True Weighted pi0->e DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPi0MCTrueWeightEnhDCA);
    
    fEtaMCTrueWeightEnhDCA = new TH2F("fEtaMCTrueWeightEnhDCA","MC True Weighted eta->e DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fEtaMCTrueWeightEnhDCA);
    
    //fPi0EtaMCTrueWeightEnhDCA = new TH2F("fPi0EtaMCTrueWeightEnhDCA","MC True Weighted pi0,eta->e DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    //fOutputList->Add(fPi0EtaMCTrueWeightEnhDCA);
    
    fPi0MCTrueHijingDCA = new TH2F("fPi0MCTrueHijingDCA","MC True Hijing pi0->e DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPi0MCTrueHijingDCA);
    
    fEtaMCTrueHijingDCA = new TH2F("fEtaMCTrueHijingDCA","MC True Hijing eta->e DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fEtaMCTrueHijingDCA);
    
    //fPi0EtaMCTrueHijingDCA = new TH2F("fPi0EtaMCTrueHijingDCA","MC True Hijing pi0+eta->e DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    //fOutputList->Add(fPi0EtaMCTrueHijingDCA);
    
    Int_t bin[4] = {60,3,2,2}; //pT, PDG, HijingOrNot, MotherOrNot
    Double_t xmin[4] = {0,0,0,0};
    Double_t xmax[4] = {30,3,2,2};
    fSprsPi0EtaWeightCal = new THnSparseD("fSprsPi0EtaWeightCal","Sparse to calculate #pi^{0} and #eta weight;p_{T};PDG ID;HijingOrNot;MotherOrNot;",4,bin,xmin,xmax);
    fSprsPi0EtaWeightCal->GetAxis(1)->SetBinLabel(1,"Pi0");
    fSprsPi0EtaWeightCal->GetAxis(1)->SetBinLabel(2,"Eta");
    fSprsPi0EtaWeightCal->GetAxis(1)->SetBinLabel(3,"Photon");
    fSprsPi0EtaWeightCal->GetAxis(2)->SetBinLabel(1,"Hijing MB");
    fSprsPi0EtaWeightCal->GetAxis(2)->SetBinLabel(2,"Enhanced");
    fSprsPi0EtaWeightCal->GetAxis(3)->SetBinLabel(1,"hasMom");
    fSprsPi0EtaWeightCal->GetAxis(3)->SetBinLabel(2,"noMom");
    fOutputList->Add(fSprsPi0EtaWeightCal);
    
    Int_t binTemp[4] = {60,200,13,3}; //pT, DCA, Mom PID
    Double_t xminTemp[4] = {0.,-0.2,-0.5,-0.5};
    Double_t xmaxTemp[4] = {30.,0.2,12.5,2.5};
    fSprsTemplates = new THnSparseD("fSprsTemplates","Sparse for Templates;p_{T};DCA;MomPID;",4,binTemp,xminTemp,xmaxTemp);
    fOutputList->Add(fSprsTemplates);
    
    /*fEtaPt_Hijing = new TH1F("fEtaPtHijing","Hijing Eta p_{T} Distribution;p_{T} (GeV/c);Counts",100,0,50);
    fOutputList->Add(fEtaPt_Hijing);
    
    fPionPt_Hijing = new TH1F("fPionPtHijing","Hijing Pion p_{T} Distribution;p_{T} (GeV/c);Counts",100,0,50);
    fOutputList->Add(fPionPt_Hijing);
    
    fEtaPt_Enh = new TH1F("fEtaPtEnh","Enhanced Eta p_{T} Distribution;p_{T} (GeV/c);Counts",100,0,50);
    fOutputList->Add(fEtaPt_Enh);
    
    fPionPt_Enh = new TH1F("fPionPtEnh","Enhanced Pion p_{T} Distribution;p_{T} (GeV/c);Counts",100,0,50);
    fOutputList->Add(fPionPt_Enh);
     */
    
    //Electron THnSparse
    //Int_t bins1[5]=  {/*280*/60,  160, 100, 100,  200}; // pT;nSigma;eop;m20;DCA
    //Double_t xmin1[5]={ /*2*/0,   -8,   0,   0, -0.2};
    //Double_t xmax1[5]={30,    8,   2,   1,  0.2};
    //fElectronSprs = new THnSparseD("Electron","Electron;pT;nSigma;eop;m20;DCA;",5,bins1,xmin1,xmax1);
    //fOutputList->Add(fElectronSprs);

    //add the list to our output file
    PostData(1, fOutputList);
}
//_____________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::UserExec(Option_t*)
{
    // Get AOD event from the analysis manager
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;
    
    //PID initialized
    fpidResponse = fInputHandler->GetPIDResponse();
    
    ////////////////
    // Centrality //
    ////////////////
    Bool_t pass = kFALSE;
    if(fCentralityMin > -0.5){
        fCentrality = CheckCentrality(fAOD,pass);
        if(!pass)return;
    }
    fCentCheck->Fill(fCentrality);
    
    ////////////////////
    // Get MC Headers //
    ////////////////////
    fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    /*if(!fMCarray){
        AliError("Array of MC particles not found");
        return;
    }*/
    fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    /*if (!fMCHeader) {
        AliError("Could not find MC Header in AOD");
        //return;
    }*/
    
    //Get NParticles from the generators
    if (fMCarray && fMCHeader) {
        GetNMCPartProduced();
    }
    
    // Make Pi0 and Eta Weight Sparse
    if (fMCarray) {
        GetPi0EtaWeight(fSprsPi0EtaWeightCal);
    }
    
    ///////////////////
    // Trigger Check //
    ///////////////////
    TString firedTrigger;
    TString TriggerEG1("EG1");
    TString TriggerDG1("DG1");
    if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();
    if(fEMCEG1){
        if(!firedTrigger.Contains(TriggerEG1))return;
        fTrigCheck->Fill(1);
    }else if(fDCalDG1){
        if(!firedTrigger.Contains(TriggerDG1))return;
        fTrigCheck->Fill(2);
    }else{fTrigCheck->Fill(0);}
    
    ////////////////
    // Mag. field //
    ////////////////
    
    Int_t MagSign = 1;
    if(fAOD->GetMagneticField()<0) MagSign = -1;
    
    //////////////////////////////
    // Event Vertex & Selection //
    //////////////////////////////
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t NcontV = pVtx->GetNContributors();
    
    //Making n track and vertex cut
    fNevents->Fill(0);
    if(NcontV<2)return;
    fNevents->Fill(1);
    if(TMath::Abs(pVtx->GetZ())>10.0) return;
    fNevents->Fill(2);
    
    fVtX->Fill(pVtx->GetX());
    fVtY->Fill(pVtx->GetY());
    fVtZ->Fill(pVtx->GetZ());
    
    //////////////////////
    // Find Kink Mother //
    //////////////////////
    Int_t numberofvertices = 100;
    if(fAOD) numberofvertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[numberofvertices];
    Int_t numberofmotherkink = 0;
    
    for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
        AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
        if(!aodvertex) continue;
        if(aodvertex->GetType()==AliAODVertex::kKink) {
            AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
            if(!mother) continue;
            Int_t idmother = mother->GetID();
            listofmotherkink[numberofmotherkink] = idmother;
            numberofmotherkink++;
        }
    }

    ////////////////
    // track loop //
    ////////////////
    Int_t nTracks(fAOD->GetNumberOfTracks());
    for(Int_t i=0; i<nTracks; i++){
        
        AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track) continue;
        
        fTrkPtB4TC->Fill(track->Pt());
        
        //////////////////////
        // Apply track cuts //
        //////////////////////
        if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //global cuts with loose DCA cut
        
        Bool_t kinkmotherpass = kTRUE;
        for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
            if(track->GetID() == listofmotherkink[kinkmother]) {
                kinkmotherpass = kFALSE;
                continue;
            }
        }
        if(!kinkmotherpass) continue; //kink rejection
        
        Double_t d0z0[2]={-999,-999}, cov[3];
        Double_t DCAxyCut = 2.4, DCAzCut = 3.2;
        if(track->GetTPCNcls() < 80) continue;
        if(track->GetITSNcls() < 3) continue;
        if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) continue;
        
        double phiMatchIts = track->Phi();
        
        if(track->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov)){
            if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
        }
        Double_t DCA = d0z0[0]*track->Charge()*MagSign;
        
        //fill the track histograms
        fTrkPt->Fill(track->Pt());
        fTrkP->Fill(track->P());
        fTrkPhi->Fill(track->Phi());
        fTrkEta->Fill(track->Eta());
        fdEdx->Fill(track->GetTPCsignal());
        
        
        ///////////////////////////
        // Match tracks to EMCal //
        ///////////////////////////
        Int_t EMCalIndex = -1;
        EMCalIndex = track->GetEMCALcluster();
        if(EMCalIndex < 0) continue;
        
        AliVCluster *clustMatch=0x0;
        clustMatch = (AliAODCaloCluster*)fAOD->GetCaloCluster(EMCalIndex);
        Double_t emcphi = -999, emceta=-999;
        Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
        
        if(clustMatch && clustMatch->IsEMCAL())
        {
            Double_t fPhiDiff = -999, fEtaDiff = -999;
            GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
            
            if(TMath::Abs(fPhiDiff) > 0.05 || TMath::Abs(fEtaDiff)> 0.05) continue;
            
            /////////////////////////////////
            //Select EMCAL or DCAL clusters//
            /////////////////////////////////
            Float_t  emcx[3]; // cluster pos
            clustMatch->GetPosition(emcx);
            TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
            emcphi = clustpos.Phi();
            emceta = clustpos.Eta();
            if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
            if(emcphi > 1.39 && emcphi < 3.265) {
                fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
                fClsEamEMCal->Fill(clustMatch->E());
            }
            if(emcphi > 4.53 && emcphi < 5.708) {
                fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327
                fClsEamDCal->Fill(clustMatch->E());
            }
            
            //----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
            if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
                if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
            
            if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
                if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
            
            if(fClsTypeEMC) fClsEamEMCal->Fill(clustMatch->E());
            if(fClsTypeDCAL) fClsEamDCal->Fill(clustMatch->E());
            
            fEMCTrkMatch->Fill(fPhiDiff,fEtaDiff);
            fTrkClsPhi->Fill(fPhiDiff);
            fTrkClsEta->Fill(fEtaDiff);
            fClsPhi->Fill(emcphi);
            fClsEta->Fill(emceta);
            
            /////////////////////////
            // Get Mother PID info //
            /////////////////////////
            Int_t fpidSort = -99;
            Bool_t kEmbEta = kFALSE;
            Bool_t kEmbPi0 = kFALSE;
            Bool_t kHijing = kFALSE;
            Int_t fMomGen = 99;
            Double_t momPt = -99;
            Int_t pdg = -99;
            Int_t ilabel = -99;
            ilabel = TMath::Abs(track->GetLabel()); //get MC label of track
            
            //if MC--------------------------------
            if(ilabel>0 && fMCarray)
            {
                //cout<<"TESTING1"<<endl;
                fMCparticle = (AliAODMCParticle*) fMCarray->At(ilabel);
                pdg = fMCparticle->GetPdgCode(); //get pid of track
                
                //if electron--------------------------------
                if(TMath::Abs(pdg)==11){
                    
                    //cout<<"TESTING2"<<endl;
                    FindMother(fMCparticle, fpidSort, kEmbEta, kEmbPi0, kHijing, momPt); //get its mom
                    
                    if (kHijing) fMomGen = 0;
                    if (kEmbPi0) fMomGen = 1;
                    if (kEmbEta) fMomGen = 2;
                    
                    //Fill template sparse
                    Double_t tempValue[4] = {-999,-999,-999,-999};
                    tempValue[0] = track->Pt();
                    tempValue[1] = DCA;
                    tempValue[2] = fpidSort;
                    tempValue[3] = fMomGen;
                    
                    fSprsTemplates->Fill(tempValue);
                    
                }
                
            }
            
            ///////////////////////////
            // Hadron Contam. Histos //
            ///////////////////////////
            Double_t nsigma = -999;
            nsigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            Double_t EovP = (clustMatch->E())/(track->P());
            Double_t M20 = clustMatch->GetM20();
            Double_t M02 = clustMatch->GetM02();
            if(nsigma<-4.) {
                if(M20>0.01 && M20<0.35) {
                    fHadronEoP->Fill(track->Pt(),EovP);
                    if(EovP>0.9 && EovP<1.2){
                        fHadronDCA->Fill(track->Pt(),DCA);
                    }
                }
            }
            
            ///////////////////
            // Electron Cuts //
            ///////////////////
            if((nsigma<-1) || (nsigma>3)) continue;
            if((M20<0.01) || (M20>0.35)) continue;
            fInclElecEoP->Fill(track->Pt(),EovP);
            
            if((EovP<0.9) || (EovP>1.2)) continue;
            
            
            /////////////////////////
            // Plot Reco Electrons //
            /////////////////////////
            fInclElecDCA->Fill(track->Pt(),DCA);
            
            //////////////////////
            // Get MC True DCAs //
            //////////////////////
            //if MC--------------------------------
            if(ilabel>0 && fMCarray)
            {
                //cout<<"TESTING1"<<endl;
                
                //if electron--------------------------------
                if(TMath::Abs(pdg)==11){
                    
                    //cout<<"TESTING2"<<endl;
                    
                    //if mom is Pi0--------------------------------
                    if(fpidSort==3) {
                        fPi0DCA->Fill(track->Pt(),DCA);
                        if(!kEmbEta && !kEmbEta && kHijing) {
                            fPi0MCTrueHijingDCA->Fill(track->Pt(),DCA);
                        }
                        if(kEmbPi0) {
                            fWeight = fPi0Weight->Eval(momPt);
                            fPi0MCTrueWeightEnhDCA->Fill(track->Pt(),DCA,fWeight);
                        }
                        if(kEmbEta) {
                            fWeight = fEtaWeight->Eval(momPt);
                            fEtaMCTrueWeightEnhDCA->Fill(track->Pt(),DCA,fWeight);
                        }
                        
                    }
                    //if Eta--------------------------------
                    if(fpidSort==4){
                        fEtaDCA->Fill(track->Pt(),DCA);
                        if(!kEmbEta && !kEmbPi0 && kHijing) fEtaMCTrueHijingDCA->Fill(track->Pt(),DCA);
                        if(kEmbEta) {
                            fWeight = fEtaWeight->Eval(momPt);
                            fEtaMCTrueWeightEnhDCA->Fill(track->Pt(),DCA,fWeight);
                        }
                    }
                    /*if(fpidSort==3 || fpidSort==4){
                        fPi0EtaDCA->Fill(track->Pt(),DCA);
                        if(!kEmbEta && !kEmbPi0) fPi0EtaMCTrueHijingDCA->Fill(track->Pt(),DCA);
                        if(kEmbEta || kEmbPi0) {
                            fWeight = fPi0EtaWeight->Eval(momPt);
                            fPi0EtaMCTrueWeightEnhDCA->Fill(track->Pt(),DCA,fWeight);
                        }
                    }*/
                    //if photon--------------------------------
                    if(fpidSort==5){
                        if(!kEmbEta && !kEmbPi0 && kHijing) fPhotonMCTrueDCA->Fill(track->Pt(),DCA);
                        if(kEmbPi0) {
                            fWeight = fPi0Weight->Eval(momPt);
                            fPi0MCTrueWeightEnhDCA->Fill(track->Pt(),DCA,fWeight);
                            //fWeight = fPi0EtaWeight->Eval(momPt);
                            //fPi0EtaMCTrueWeightEnhDCA->Fill(track->Pt(),DCA,fWeight);
                        }
                        if(kEmbEta) {
                            fWeight = fEtaWeight->Eval(momPt);
                            fEtaMCTrueWeightEnhDCA->Fill(track->Pt(),DCA,fWeight);
                            //fWeight = fPi0EtaWeight->Eval(momPt);
                            //fPi0EtaMCTrueWeightEnhDCA->Fill(track->Pt(),DCA,fWeight);
                        }
                    }

                }
                
            }
        
            InvMassCheck(i, track, d0z0, MagSign);
            
            //Make incl electron and photonic electron plots
            /*if(nsigma>-1 && nsigma<3) {
                if(M20>0.01 && M20<0.35) {
                    fInclElecEoP->Fill(track->Pt(),EovP);
                    if(EovP>0.9 && EovP<1.2){
                        InvMassCheck(i, track, d0z0, MagSign);
                        fInclElecDCA->Fill(track->Pt(),DCA);
                    }
                }
            }*/
            
            ////////////////////////////////////////////
            // Label Electrons - Mother and Generator //
            ////////////////////////////////////////////
            /*if(fMCarray){
                
                Int_t fpidSort = 3;
                Int_t pdg = -999;
                Int_t pidM = -99;
                Int_t pid_ele = -99;
                Int_t ilabelM = -1;
                Int_t ilabel = TMath::Abs(track->GetLabel()); //get MC label of track
                if(ilabel>0 && fMCarray)
                {
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(ilabel);
                    pdg = fMCparticle->GetPdgCode(); //get pid of track
                    
                    if(TMath::Abs(pdg)==11)pid_ele = 1.0; //if electron...
                    if(pid_ele==1.0){
                        FindMother(fMCparticle, ilabelM, pidM, fpidSort);//get its mom
                        if(fpidSort==3) fDalitzDCA->Fill(track->Pt(),DCA);
                        if(fpidSort==4) fPhotonicDCA->Fill(track->Pt(),DCA);
                    }
                    
                }
            }

            
            /////////////////////
            // Electron sparse //
            /////////////////////
            
            
            Double_t fvalueElectron[5] = {-999,-999,-999,-999,-999};
            fvalueElectron[0] = track->Pt();
            fvalueElectron[1] = nsigma;
            fvalueElectron[2] = EovP;
            fvalueElectron[3] = M20;
            fvalueElectron[4] = DCA;
            
            fvalueElectron[5] = -999;
            if (fClsTypeEMC){
                fvalueElectron[5] = 0; //0=EMCal, 1=DCal
            }
            if (fClsTypeDCAL){
                fvalueElectron[5] = 1; //0=EMCal, 1=DCal
            }*/
            //fElectronSprs->Fill(fvalueElectron);
            
        }
        
    }
    
    //save the data gathered in this iteration
    PostData(1,fOutputList);
}
//___________________________________________
Double_t AliAnalysisTaskTPCCalBeauty::CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass)
{
    //check centrality, Run 2
    if(fAOD)fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
    if(!fMultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    }else{
        fCentrality = fMultSelection->GetMultiplicityPercentile("V0M", false);
    }
    
    //AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
    //if(!header) AliFatal("Not a standard AOD");
    //fMultiplicity = header->GetRefMultiplicity();

    if ((fCentrality <= fCentralityMin) || (fCentrality > fCentralityMax)){
        //fCentralityNoPass->Fill(fCentrality);
        //fMultiplicityNoPass->Fill(fMultiplicity);
        centralitypass = kFALSE;
    }else{
        //fCentralityPass->Fill(fCentrality);
        //fMultiplicityPass->Fill(fMultiplicity);
        centralitypass = kTRUE;
    }
    
    return fCentrality;
}
//_________________________________________________________________
Bool_t AliAnalysisTaskTPCCalBeauty::GetNMCPartProduced()
{
    //Get number of MC particles produced by generators.
    
    //list of headers
    TList *lh = fMCHeader->GetCocktailHeaders();
    fNtotMCpart = 0;
    fNembMCpi0 = 0;
    fNembMCeta = 0;
    fNpureMC = 0;
    TString MCgen;
    TString embpi0("pi");
    TString embeta("eta");
    
    if(!lh){
        AliError("no MC header");
        return (0);
    }
    //loop through headers
    for(int igene=0; igene<lh->GetEntries(); igene++)
    {
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
        if(!gh) continue;
        
        MCgen =  gh->GetName();
        //cout << "Gen name, N produced = " << gh->GetName() << ", " << gh->NProduced() << endl;
        
        //the first header is pure MC
        if(igene==0) fNpureMC = gh->NProduced();  // generated by HIJING
        
        //if(MCgen.Contains(embpi0))cout << MCgen << endl;
        //if(MCgen.Contains(embeta))cout << MCgen << endl;
        
        //if header has "pi" or "eta", note number in stack where it starts
        if(MCgen.Contains(embpi0))fNembMCpi0 = fNtotMCpart;
        if(MCgen.Contains(embeta))fNembMCeta = fNtotMCpart;
        fNtotMCpart += gh->NProduced();
    }
    //cout << "fNpureMC, fNembMCpi0, fNembMCeta, fNtotMCpart : " <<fNpureMC << ", " << fNembMCpi0 << ", " << fNembMCeta << ", " << fNtotMCpart << endl;
    
    return kTRUE;
}
//_________________________________________
void AliAnalysisTaskTPCCalBeauty::GetPi0EtaWeight(THnSparse *SparseWeight)
{
    //Get pi0 and eta information for weight calculation
    
    for(int imc=0; imc< fNtotMCpart; imc++)
    {
        AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(imc);
        if(TMath::Abs(AODMCtrack->Eta()) > 0.9) continue;
        
        //-------Get PDG
        Int_t TrackPDG = TMath::Abs(AODMCtrack->GetPdgCode());
        if((TrackPDG != 111) && (TrackPDG != 221) && (TrackPDG != 22)) continue;
        
        Double_t fPartPDGid = -999;
        if (TrackPDG == 111) fPartPDGid = 0.2; //pi0
        if (TrackPDG == 221) fPartPDGid = 1.2; //eta
        if (TrackPDG == 22) fPartPDGid = 2.2; //photon
        
        //-------Check if the particle is from hijing or not
        Double_t fFromHijing = 0.2; //Hijing MB
        if(imc >= fNpureMC)fFromHijing = 1.2; //enhanced
        
        //------Get type of the particle
        Double_t fMother = 0.2; //has mother
        Int_t motherlabel = AODMCtrack->GetMother();
        if(motherlabel<0) fMother = 1.2; //No mom
        
        Double_t fvalue[4] = {-999,-999,-999,-999};
        fvalue[0] = AODMCtrack->Pt();
        fvalue[1] = fPartPDGid;
        fvalue[2] = fFromHijing;
        fvalue[3] = fMother;
        
        SparseWeight->Fill(fvalue);
    }
}
//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
{
    // Calculate phi and eta difference between a track and a cluster. The position of the track is obtained on the EMCAL surface

    phidiff = 999;
    etadiff = 999;
    
    if (!t||!v) return;

    Double_t veta = t->GetTrackEtaOnEMCal();
    Double_t vphi = t->GetTrackPhiOnEMCal();

    Float_t pos[3] = {0};
    v->GetPosition(pos);
    TVector3 cpos(pos);
    Double_t ceta     = cpos.Eta();
    Double_t cphi     = cpos.Phi();
    etadiff=veta-ceta;
    phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
}
//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::FindMother(AliAODMCParticle* part, Int_t &fpidSort, Bool_t &kEmbEta, Bool_t &kEmbPi0, Bool_t &kHijing, Double_t &momPt)
{
    //gets the pid of mother track
    
    //Get the mother, grandma, and great grandma pid
    Int_t pdg = -999;
    Int_t pidM = -99;
    Int_t pid_ele = -99;
    Int_t ilabelM = -1;
    Int_t ilabelGM = -1;
    Int_t ilabelGGM = -1;
    
    //cout<<"TESTING3"<<endl;
    
    //Get the mother, grandma, and great grandma pid
    ilabelM = part->GetMother(); //get MC label for Mom
    if(ilabelM>0){
        AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(ilabelM); //get mom particle
        pidM = TMath::Abs(partM->GetPdgCode()); //ask for the Mom's pid
        momPt = partM->Pt();
        
        if(ilabelM<fNpureMC) kHijing = kTRUE; //mark whether mom is from Hijing
        
        ilabelGM = partM->GetMother();//get MC for Grandma
        
        //sort according to mother
        if(pidM>500 && pidM<599){
            fpidSort = 1; //Mom is B
        }
        else if(pidM>400 && pidM<499){
            fpidSort = 2; //Mom is D
        }
        else if(pidM==111){
            fpidSort = 3; //Mom is pi0
            if(ilabelM >= fNembMCpi0 && ilabelM < fNembMCeta  && ilabelGM<0) kEmbPi0 = kTRUE;
        }
        else if(pidM==221){
            fpidSort = 4; //Mom is eta
            if(ilabelM >= fNembMCeta && ilabelM < fNtotMCpart && ilabelGM<0) kEmbEta = kTRUE;
        }
        else if(pidM==22){
            fpidSort = 5; //Mom is gamma
        }
        else if(pidM==443){
            fpidSort = 6; //Mom is J/psi
        }
        else if(pidM==23){
            fpidSort = 7; //Mom is Z
        }
        else if(pidM==24){
            fpidSort = 8; //Mom is W
        }
        else if(pidM>4000 && pidM<4999){
            fpidSort = 9; //Mom is c Baryon
        }
        else if(pidM>5000 && pidM<5999){
            fpidSort = 10; //Mom is b Baryon
        }
        else{
            //cout<<"TESTING4"<<endl;
            fpidSort = 11; //Mom is something else
        }
        
        if(ilabelGM>0){
            AliAODMCParticle *partGM = (AliAODMCParticle*)fMCarray->At(ilabelGM); // get GMa particle
            Int_t pidGM = TMath::Abs(partGM->GetPdgCode()); //ask for grandma's pid
            ilabelGGM = partGM->GetMother();//get MC for Great Grandma
            
            //check if pi0 grandma is eta
            if(pidM==111){
                if(pidGM==221){
                    fpidSort = 4; //label as eta electron
                    if(ilabelGM >= fNembMCeta && ilabelGM < fNtotMCpart && ilabelGGM<0) kEmbEta = kTRUE; //GMa is eta
                }
            }
            //check if gamma grandma is eta/pion
            if(pidM==22){
                if(pidGM==221){
                    if(ilabelGM >= fNembMCeta && ilabelGM<fNtotMCpart && ilabelGGM<0) kEmbEta = kTRUE; //GMa is eta
                }
                if(pidGM==111){
                    if(ilabelM >= fNembMCpi0 && ilabelM<fNembMCeta && ilabelGGM<0) kEmbPi0 = kTRUE; //GMa is pi0
                }
            }
            
            
            //check if D grandma is B
            if(pidM>400 && pidM<499){
                if(pidGM>500 && pidGM<599){
                    fpidSort = 1; //GMa is B
                }
            }
            if(ilabelGGM>0 && fpidSort!=2){
                AliAODMCParticle *partGGM = (AliAODMCParticle*)fMCarray->At(ilabelGGM); // get GGMa particle
                Int_t pidGGM = TMath::Abs(partGGM->GetPdgCode()); //ask for ggma's pid
                
                //check if D great grandma is B
                if(pidM>400 && pidM<499){
                    if(pidGGM>500 && pidGGM<599){
                        fpidSort = 1; //GGMa is B
                    }
                }
            }
        }
    }else{
        fpidSort = 12; //No mother
    }
}

//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::InvMassCheck(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign)
{
    // Flags photonic electrons with inv mass cut
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t d0z0Asso[2]={-999,-999}, covAsso[3];
    Double_t DCAxyCut = 0.25, DCAzCut = 1;
    Int_t fPDGe1 = 11, fPDGe2 = 11;
    
    Double_t ptAsso=-999., nsigmaAsso=-999.;
    Int_t chargeAsso=0;
    Int_t charge=track->Charge();
    Double_t mass=-999., width = -999;
    Int_t MassCorrect;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    Int_t Nuls=0, Nls=0;
    
    Int_t ntracks = fAOD->GetNumberOfTracks();
    for (int jtrack=0; jtrack<ntracks; jtrack++) {
        if (jtrack==itrack) {continue;} //asso track != selected track
        
        fFlagLS=kFALSE;
        fFlagULS=kFALSE;
        
        AliAODTrack *trackAsso = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(jtrack));
        if(!trackAsso) continue;
        if(!trackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
        if(trackAsso->GetTPCNcls() < 80) continue;
        
        nsigmaAsso = fpidResponse->NumberOfSigmasTPC(trackAsso, AliPID::kElectron);
        ptAsso = trackAsso->Pt();
        chargeAsso = trackAsso->Charge();
        
        //Some cuts on the associated track
        if(ptAsso < 0.3) continue;
        if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
        if(nsigmaAsso < -3 || nsigmaAsso > 3) continue;
        
        if(trackAsso->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0Asso, covAsso))
            if(TMath::Abs(d0z0Asso[0]) > DCAxyCut || TMath::Abs(d0z0Asso[1]) > DCAzCut) continue;
        
        if(charge>0) fPDGe1 = -11; //-11 in PDG is for positron, just to be confusing
        if(chargeAsso>0) fPDGe2 = -11;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        AliKFParticle::SetField(fAOD->GetMagneticField());
        
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*trackAsso, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        MassCorrect = recg.GetMass(mass,width); //returns 1, not the mass
        if(fFlagLS && track->Pt()>1) fInvmassLS->Fill(mass);
        if(fFlagULS && track->Pt()>1) fInvmassULS->Fill(mass);
        
        if(fFlagLS && mass<0.1) Nls++;
        if(fFlagULS && mass<0.1) Nuls++;
        
        if (fFlagULS && mass<0.1 && track->Pt()>1) {
            fULSdcaBelow->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign);
            
        }else if(fFlagLS && mass<0.1 && track->Pt()>1){
            fLSdcaBelow->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign);
        }
        
    }
    
    fPhotonicElecYield->Fill(track->Pt(),Nuls-Nls);
}
//_____________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::Terminate(Option_t *)
{
    // terminate
}
//_____________________________________________________________________
