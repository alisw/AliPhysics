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
fClsEAll(0),
fClsEamElecEMC(0),
fClsEamElecDC(0),
fTrkPhi(0),
fTrkEta(0),
fdEdx(0),
fCentCheck(0),
fTrigCheck(0),
fEMCTrkMatch(0),
fInvmassLS(0),
fInvmassULS(0),
fInvmassLSWeightEnhEta(0),
fInvmassULSWeightEnhEta(0),
fInvmassLSWeightEnhPi0(0),
fInvmassULSWeightEnhPi0(0),
fInvmassLSHijingEta(0),
fInvmassULSHijingEta(0),
fInvmassLSHijingPi0(0),
fInvmassULSHijingPi0(0),
fInvmassLSHijingPhoton(0),
fInvmassULSHijingPhoton(0),
fInvmassLSEnhPhoton(0),
fInvmassULSEnhPhoton(0),
fULSdcaBelow(0),
fLSdcaBelow(0),
fLSWeightEnhEta(0),
fULSWeightEnhEta(0),
fLSWeightEnhPi0(0),
fULSWeightEnhPi0(0),
fLSHijingEta(0),
fULSHijingEta(0),
fLSHijingPi0(0),
fULSHijingPi0(0),
fLSHijingPhoton(0),
fULSHijingPhoton(0),
fLSEnhPhoton(0),
fULSEnhPhoton(0),
fPhotonicDCA(0),
fInclElecDCA(0),
fInclElecDCAnoSign(0),
fInclElecEoP(0),
fHadronEoP(0),
fHadronDCA(0),
fHadronCamDCA(0),
fPi0Weight(0),
fEtaWeight(0),
fDWeight(0),
fEnhEtaDCA(0),
fEnhEtaWeightedPt(0),
fEnhPi0DCA(0),
fEnhPi0WeightedPt(0),
fEtaHijingDCA(0),
fEtaHijingPt(0),
fPi0HijingDCA(0),
fPi0HijingPt(0),
fPhotonHijingDCA(0),
fPhotonHijingPt(0),
fEnhPhotonDCA(0),
fEnhPhotonWeightedPt(0),
ComboNumWeight(0),
ComboNumNoWeight(0),
ComboDenomWeight(0),
ComboDenomNoWeight(0),
DMesonPDG(0),
fD0MesonPt(0),
fD0MesonFromDStarPt(0),
fDPlusMesonPt(0),
fDsMesonPt(0),
fDStarMesonPt(0),
fAllDMesonPt(0),
fLambdaCPt(0),
fEtaCPt(0),
fCBaryonPt(0),
fPromptD0DCAWeight(0),
fD0FromDStarDCAWeight(0),
fPromptD0DCANoWeight(0),
fD0FromDStarDCANoWeight(0),
fNtotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),
fSprsPi0EtaWeightCal(0),
fSprsTemplatesNoWeight(0),
fSprsTemplatesWeight(0),
fDTemplateWeight(0),
fDTemplateNoWeight(0)
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
fClsEAll(0),
fClsEamElecEMC(0),
fClsEamElecDC(0),
    fTrkPhi(0),
    fTrkEta(0),
    fdEdx(0),
    fCentCheck(0),
    fTrigCheck(0),
    fEMCTrkMatch(0),
    fInvmassLS(0),
    fInvmassULS(0),
fInvmassLSWeightEnhEta(0),
fInvmassULSWeightEnhEta(0),
fInvmassLSWeightEnhPi0(0),
fInvmassULSWeightEnhPi0(0),
fInvmassLSHijingEta(0),
fInvmassULSHijingEta(0),
fInvmassLSHijingPi0(0),
fInvmassULSHijingPi0(0),
fInvmassLSHijingPhoton(0),
fInvmassULSHijingPhoton(0),
fInvmassLSEnhPhoton(0),
fInvmassULSEnhPhoton(0),
    fULSdcaBelow(0),
    fLSdcaBelow(0),
fLSWeightEnhEta(0),
fULSWeightEnhEta(0),
fLSWeightEnhPi0(0),
fULSWeightEnhPi0(0),
fLSHijingEta(0),
fULSHijingEta(0),
fLSHijingPi0(0),
fULSHijingPi0(0),
fLSHijingPhoton(0),
fULSHijingPhoton(0),
fLSEnhPhoton(0),
fULSEnhPhoton(0),
    fPhotonicDCA(0),
    fInclElecDCA(0),
fInclElecDCAnoSign(0),
fInclElecEoP(0),
fHadronEoP(0),
fHadronDCA(0),
fHadronCamDCA(0),
    fPi0Weight(0),
    fEtaWeight(0),
fDWeight(0),
fEnhEtaDCA(0),
fEnhEtaWeightedPt(0),
fEnhPi0DCA(0),
fEnhPi0WeightedPt(0),
fEtaHijingDCA(0),
fEtaHijingPt(0),
fPi0HijingDCA(0),
fPi0HijingPt(0),
fPhotonHijingDCA(0),
fPhotonHijingPt(0),
fEnhPhotonDCA(0),
fEnhPhotonWeightedPt(0),
ComboNumWeight(0),
ComboNumNoWeight(0),
ComboDenomWeight(0),
ComboDenomNoWeight(0),
DMesonPDG(0),
fD0MesonPt(0),
fD0MesonFromDStarPt(0),
fDPlusMesonPt(0),
fDsMesonPt(0),
fDStarMesonPt(0),
fAllDMesonPt(0),
fLambdaCPt(0),
fEtaCPt(0),
fCBaryonPt(0),
fPromptD0DCAWeight(0),
fD0FromDStarDCAWeight(0),
fPromptD0DCANoWeight(0),
fD0FromDStarDCANoWeight(0),
    fNtotMCpart(0),
    fNpureMC(0),
    fNembMCpi0(0),
    fNembMCeta(0),
    fSprsPi0EtaWeightCal(0),
    fSprsTemplatesNoWeight(0),
fSprsTemplatesWeight(0),
fDTemplateWeight(0),
fDTemplateNoWeight(0)
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
    
    fClsEAll = new TH1F("fClsEAll","Cluster Energy, All Clusters;Cluster E;Counts",250,0.,50);
    fOutputList->Add(fClsEAll);
    fClsEamElecEMC = new TH1F("fClsEamElecEMC","Cluster Energy of e- after track matching to DCal;Cluster E;Counts",250,0.,50);
    fOutputList->Add(fClsEamElecEMC);
    fClsEamElecDC = new TH1F("fClsEamElecDC","Cluster Energy of e- after track matching to DCal;Cluster E;Counts",250,0.,50);
    fOutputList->Add(fClsEamElecDC);
    
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
    
    fInvmassLSWeightEnhEta = new TH1F("fInvmassLSWeightEnhEta", "Inv mass of Weighted Enh Eta LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fInvmassLSWeightEnhEta->Sumw2();
    fOutputList->Add(fInvmassLSWeightEnhEta);
    fInvmassULSWeightEnhEta = new TH1F("fInvmassULSWeightEnhEta", "Inv mass of Weighted Enh Eta ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fInvmassULSWeightEnhEta->Sumw2();
    fOutputList->Add(fInvmassULSWeightEnhEta);
    
    fInvmassLSWeightEnhPi0 = new TH1F("fInvmassLSWeightEnhPi0", "Inv mass of Weighted Enh Pi0 LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fInvmassLSWeightEnhPi0->Sumw2();
    fOutputList->Add(fInvmassLSWeightEnhPi0);
    fInvmassULSWeightEnhPi0 = new TH1F("fInvmassULSWeightEnhPi0", "Inv mass of Weighted Enh Pi0 ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassULSWeightEnhPi0);
    fInvmassULSWeightEnhPi0->Sumw2();
    
    fInvmassLSHijingEta = new TH1F("fInvmassLSHijingEta", "Inv mass of Hijing Eta LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassLSHijingEta);
    fInvmassULSHijingEta = new TH1F("fInvmassULSHijingEta", "Inv mass of Hijing Eta ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassULSHijingEta);

    fInvmassLSHijingPi0 = new TH1F("fInvmassLSHijingPi0", "Inv mass of Hijing Pi0 LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassLSHijingPi0);
    fInvmassULSHijingPi0 = new TH1F("fInvmassULSHijingPi0", "Inv mass of Hijing Pi0 ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassULSHijingPi0);
    
    fInvmassLSHijingPhoton = new TH1F("fInvmassLSHijingPhoton", "Inv mass of Hijing Photon LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassLSHijingPhoton);
    fInvmassULSHijingPhoton = new TH1F("fInvmassULSHijingPhoton", "Inv mass of Hijing Photon ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassULSHijingPhoton);
    
    fInvmassLSEnhPhoton = new TH1F("fInvmassLSEnhPhoton", "Inv mass of Photon LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassLSEnhPhoton);
    fInvmassLSEnhPhoton->Sumw2();
    fInvmassULSEnhPhoton = new TH1F("fInvmassULSEnhPhoton", "Inv mass of Photon ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fInvmassULSEnhPhoton->Sumw2();
    fOutputList->Add(fInvmassULSEnhPhoton);
    
    fULSdcaBelow = new TH2F("fULSdcaBelow","ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fULSdcaBelow);
    
    fLSdcaBelow = new TH2F("fLSdcaBelow","LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fLSdcaBelow);
    
    fLSWeightEnhEta = new TH1F("fLSWeightEnhEta","Weighted Enh Eta LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fLSWeightEnhEta->Sumw2();
    fOutputList->Add(fLSWeightEnhEta);
    
    fULSWeightEnhEta = new TH1F("fULSWeightEnhEta","Weighted Enh Eta ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fULSWeightEnhEta->Sumw2();
    fOutputList->Add(fULSWeightEnhEta);
    
    fLSWeightEnhPi0 = new TH1F("fLSWeightEnhPi0","Weighted Enh Pi0 LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fLSWeightEnhPi0->Sumw2();
    fOutputList->Add(fLSWeightEnhPi0);
    
    fULSWeightEnhPi0 = new TH1F("fULSWeightEnhPi0","Weighted Enh Pi0 ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fULSWeightEnhPi0->Sumw2();
    fOutputList->Add(fULSWeightEnhPi0);
    
    fLSHijingEta = new TH1F("fLSHijingEta","Hijing Eta LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fOutputList->Add(fLSHijingEta);
    
    fULSHijingEta = new TH1F("fULSHijingEta","Hijing Eta ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fOutputList->Add(fULSHijingEta);
    
    fLSHijingPi0 = new TH1F("fLSHijingPi0","Hijing Pi0 LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30.);
    fOutputList->Add(fLSHijingPi0);
    
    fULSHijingPi0 = new TH1F("fULSHijingPi0","Hijing Pi0 ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fOutputList->Add(fULSHijingPi0);
    
    fLSHijingPhoton = new TH1F("fLSHijingPhoton","Hijing Photon LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fOutputList->Add(fLSHijingPhoton);
    
    fULSHijingPhoton = new TH1F("fULSHijingPhoton","Hijing Photon ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fOutputList->Add(fULSHijingPhoton);
    
    fLSEnhPhoton = new TH1F("fLSEnhPhoton","All Photon LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fLSEnhPhoton->Sumw2();
    fOutputList->Add(fLSEnhPhoton);
    
    fULSEnhPhoton = new TH1F("fULSEnhPhoton","All Photon ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); counts;", 60,0,30.);
    fULSEnhPhoton->Sumw2();
    fOutputList->Add(fULSEnhPhoton);
    
    fPhotonicDCA = new TH2F("fPhotonicDCA","Photonic DCA using MC PID; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPhotonicDCA);
    
    fInclElecDCA = new TH2F("fInclElecDCA","Incl Elec DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fInclElecDCA);
    
    fInclElecDCAnoSign = new TH2F("fInclElecDCAnoSign","Incl Elec DCA (no Charge); p_{T}(GeV/c); DCAxMagField; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fInclElecDCAnoSign);
    
    fInclElecEoP = new TH2F("fInclElecEoP","Incl Elec E/p; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fInclElecEoP);
    
    fHadronEoP = new TH2F("fHadronEoP","Unscaled Hadron E/p; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fHadronEoP);
    
    fHadronDCA = new TH2F("fHadronDCA","Unscaled Hadron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fHadronDCA);
    
    fHadronCamDCA = new TH2F("fHadronCamDCA","Unscaled Hadron DCA, no E/p cut; p_{T}(GeV/c); DCAxMagField; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fHadronCamDCA);
    
    
    fPi0Weight = new TF1("fPi0Weight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    fEtaWeight = new TF1("fEtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    //fPi0EtaWeight = new TF1("fPi0EtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    //fPi0Weight->SetParameters(3.10605e+04,-2.88504e-01,4.90349e-03,4.82424e-01,3.88146e+00);
    //fEtaWeight->SetParameters(1.48495e+03,-2.20500e-01,3.35819e-03,1.14265e+00,4.06488e+00);
    //fPi0EtaWeight->SetParameters(3.68627e+03,-2.15364e-01,3.26753e-03,1.05213e+00,4.20391e+00);
    
    //Deepa's Weight
    fPi0Weight->SetParameters(8.96715e+02,-1.77016e-01,2.47560e-03,1.50783e+00,4.41819e+00);
    fEtaWeight->SetParameters(4.18393e+02,-5.20750e-02,-3.11517e-04,2.04739e+00,5.30788e+00);
    
    //D Meson pt weighting
    Int_t nbins = 13;
    Double_t xbins[14] = {1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,16.,24.,36.,50.};
    //Double_t err[13] = {};
    fDWeight = new TH1F("fDWeight","D^{0}_data/AllD_MC;p_{T} (GeV/c);Weight;",nbins,xbins);
    //Double_t ratio[13] = {2.03552,1.0201,0.45925,0.211574,0.11987,0.0898116,0.0631282,0.0546798,0.0477205,0.0410021,0.0307936,0.0398483,0.0175335};
    //Double_t err[13] = {0.541651,0.146443,0.0498454,0.024907,0.01438,0.0107908,0.00848616,0.0061723,0.00587082,0.00566712,0.00597994,0.00811015,0.00693105};
    Double_t ratio[13] = {0.106888,0.0650239,0.0343858,0.0172579,0.00957876,0.00640323,0.00399907,0.00269269,0.00163078,0.000942387,0.000441093,0.000353811,0.000143011};
    Double_t err[13] = {0.0284416,0.0093333,0.00373075,0.00203067,0.00114824,0.000768388,0.00053676,0.000303334,0.000199878,0.000129785,8.53822e-05,7.13313e-05,5.61316e-05};
    for (int idata=1; idata<14; idata++) {
        fDWeight->SetBinContent(idata,ratio[idata-1]);
        fDWeight->SetBinError(idata,err[idata-1]);
    }
    fDWeight->Sumw2();
    fOutputList->Add(fDWeight);
    
    fPi0DCA = new TH2F("fPi0DCA","Pi0 DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPi0DCA);
    
    fEtaDCA = new TH2F("fEtaDCA","Eta DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fEtaDCA);
    
    fEnhEtaDCA = new TH2F("fEnhEtaDCA","Enh Eta DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fEnhEtaDCA->Sumw2();
    fOutputList->Add(fEnhEtaDCA);
    fEnhEtaWeightedPt = new TH1F("fEnhEtaWeightedPt","Enh Eta Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
    fEnhEtaWeightedPt->Sumw2();
    fOutputList->Add(fEnhEtaWeightedPt);
    fEnhPi0DCA = new TH2F("fEnhPi0DCA","Enh Pi0 DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fEnhPi0DCA->Sumw2();
    fOutputList->Add(fEnhPi0DCA);
    fEnhPi0WeightedPt = new TH1F("fEnhPi0WeightedPt","Enh Pi0 Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
    fEnhPi0WeightedPt->Sumw2();
    fOutputList->Add(fEnhPi0WeightedPt);
    fEtaHijingDCA = new TH2F("fEtaHijingDCA","Hijing Eta DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fEtaHijingDCA);
    fEtaHijingPt = new TH1F("fEtaHijingPt","Hijing Eta Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
    fOutputList->Add(fEtaHijingPt);
    fPi0HijingDCA = new TH2F("fPi0HijingDCA","Hijing Pi0 DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPi0HijingDCA);
    fPi0HijingPt = new TH1F("fPi0HijingPt","Hijing Pi0 Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
    fOutputList->Add(fPi0HijingPt);
    fPhotonHijingDCA = new TH2F("fPhotonHijingDCA","Hijing Photon DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPhotonHijingDCA);
    fPhotonHijingPt = new TH1F("fPhotonHijingPt","Hijing Photon Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
    fOutputList->Add(fPhotonHijingPt);
    fEnhPhotonDCA = new TH2F("fEnhPhotonDCA","Enh Photon DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fEnhPhotonDCA->Sumw2();
    fOutputList->Add(fEnhPhotonDCA);
    fEnhPhotonWeightedPt = new TH1F("fEnhPhotonWeightedPt","Enh Eta Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
    fEnhPhotonWeightedPt->Sumw2();
    fOutputList->Add(fEnhPhotonWeightedPt);
    
    ComboNumWeight = new TH1F("ComboNumWeight","Eff Num with Weight; p_{T}(GeV/c); counts;", 60,0,30.);
    ComboNumWeight->Sumw2();
    fOutputList->Add(ComboNumWeight);
    ComboNumNoWeight = new TH1F("ComboNumNoWeight","Eff Num Without Weight; p_{T}(GeV/c); counts;", 60,0,30.);
    ComboNumNoWeight->Sumw2();
    fOutputList->Add(ComboNumNoWeight);
    ComboDenomWeight = new TH1F("ComboDenomWeight","Eff Denom with Weight; p_{T}(GeV/c); counts;", 60,0,30.);
    ComboDenomWeight->Sumw2();
    fOutputList->Add(ComboDenomWeight);
    ComboDenomNoWeight = new TH1F("ComboDenomNoWeight","Eff Denom without Weight; p_{T}(GeV/c); counts;", 60,0,30.);
    ComboDenomNoWeight->Sumw2();
    fOutputList->Add(ComboDenomNoWeight);
    
    DMesonPDG = new TH1F("DMesonPDG","D Meson PDG values; abs(PDG); counts;", 51,399.5,450.5);
    //DMesonPDG->Sumw2();
    fOutputList->Add(DMesonPDG);
    
    fD0MesonPt = new TH1F("fD0MesonPt","D^{0} Meson Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
    fD0MesonPt->Sumw2();
    fOutputList->Add(fD0MesonPt);
    
    fD0MesonFromDStarPt = new TH1F("fD0MesonFromDStarPt","D^{0} from D* Meson Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
    fD0MesonFromDStarPt->Sumw2();
    fOutputList->Add(fD0MesonFromDStarPt);
    
    fDPlusMesonPt = new TH1F("fDPlusMesonPt","D+ Meson Spectrum after track cuts; p_{T}(GeV/c); counts;",100,0,50.);
    fDPlusMesonPt->Sumw2();
    fOutputList->Add(fDPlusMesonPt);
    
    fDsMesonPt = new TH1F("fDsMesonPt","Ds Meson Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
    fDsMesonPt->Sumw2();
    fOutputList->Add(fDsMesonPt);

    fDStarMesonPt = new TH1F("fDStarMesonPt","D* Meson Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
    fDStarMesonPt->Sumw2();
    fOutputList->Add(fDStarMesonPt);
    
    fAllDMesonPt = new TH1F("fAllDMesonPt","All D Meson Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
    fAllDMesonPt->Sumw2();
    fOutputList->Add(fAllDMesonPt);
    
    fLambdaCPt = new TH1F("fLambdaCPt","Lambda_c Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
    fLambdaCPt->Sumw2();
    fOutputList->Add(fLambdaCPt);
    
    fEtaCPt = new TH1F("fEtaCPt","Eta_c Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
    fEtaCPt->Sumw2();
    fOutputList->Add(fEtaCPt);
    
    fCBaryonPt = new TH1F("fCBaryonPt","Charm Baryon Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
    fCBaryonPt->Sumw2();
    fOutputList->Add(fCBaryonPt);
    
    fPromptD0DCAWeight = new TH2F("fPromptD0DCAWeight","Prompt D0 DCA with Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPromptD0DCAWeight);

    fD0FromDStarDCAWeight = new TH2F("fD0FromDStarDCAWeight","Prompt D0 DCA with Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fD0FromDStarDCAWeight);
    
    fPromptD0DCANoWeight = new TH2F("fPromptD0DCANoWeight","Prompt D0 DCA without Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fPromptD0DCANoWeight);
    
    fD0FromDStarDCANoWeight = new TH2F("fD0FromDStarDCANoWeight","Prompt D0 DCA without Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., 200,-0.2,0.2);
    fOutputList->Add(fD0FromDStarDCANoWeight);
    
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
    
    Int_t binTemp[5] = {60,200,19,3,50}; //pT, DCA, Mom PID, Mom Gen, mompT
    Double_t xminTemp[5] = {0.,-0.2,0.5,-0.5,0.,};
    Double_t xmaxTemp[5] = {30.,0.2,19.5,2.5,50.};
    fSprsTemplatesWeight = new THnSparseD("fSprsTemplatesWeight","Sparse for Templates, D meson weight applied;p_{T};DCA;MomPID;MomGen;",5,binTemp,xminTemp,xmaxTemp);
    fOutputList->Add(fSprsTemplatesWeight);
    fSprsTemplatesWeight->Sumw2();
    fSprsTemplatesNoWeight = new THnSparseD("fSprsTemplatesNoWeight","Sparse for Templates, No weight applied;p_{T};DCA;MomPID;MomGen;",5,binTemp,xminTemp,xmaxTemp);
    fOutputList->Add(fSprsTemplatesNoWeight);
    fSprsTemplatesNoWeight->Sumw2();
    
    fDTemplateWeight = new TH2F("fDTemplateWeight","D Meson DCA template", 100,0,50., 200,-0.2,0.2);
    fOutputList->Add(fDTemplateWeight);
    
    fDTemplateNoWeight = new TH2F("fDTemplateNoWeight","D Meson DCA template w/Weight", 100,0,50., 200,-0.2,0.2);
    fOutputList->Add(fDTemplateNoWeight);
    
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
        //return;
    }
    fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!fMCHeader) {
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
    
        //Fill Charm species pT histos
        for(int i=0; i<fNtotMCpart; i++)
        {
            //if (i<fNpureMC) continue; //reject plain Hijing
            
            Bool_t fromDStar = kFALSE;
            
            AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(i);
            if(TMath::Abs(AODMCtrack->Eta()) > 0.9) continue;
        
            //-------Get PDG
            Int_t TrackPDG = TMath::Abs(AODMCtrack->GetPdgCode());
        
            Int_t ilabelM = AODMCtrack->GetMother();
            if (ilabelM>0) {
                AliAODMCParticle *momPart = (AliAODMCParticle*)fMCarray->At(ilabelM); //get mom particle
                Int_t pidM = TMath::Abs(momPart->GetPdgCode());
                if (pidM==413) fromDStar = kTRUE;
                if (pidM>500 && pidM<599) {
                    continue; //reject beauty feed down
                }
                if (pidM>5000 && pidM<5999) {
                    continue; //reject beauty feed down
                }
                Int_t ilabelGM = momPart->GetMother();
                if (ilabelGM>0) {
                    AliAODMCParticle *gmomPart = (AliAODMCParticle*)fMCarray->At(ilabelGM);//get grandma particle
                    Int_t pidGM = TMath::Abs(momPart->GetPdgCode());
                    if (pidGM>500 && pidGM<599) {
                        continue; //reject beauty feed down
                    }
                    if (pidGM>5000 && pidGM<5999) {
                        continue; //reject beauty feed down
                    }
                }
            
                if(TrackPDG>4000 && TrackPDG<4999) fCBaryonPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==4122)fLambdaCPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==4132)fEtaCPt->Fill(AODMCtrack->Pt());
        
                if (TrackPDG<400 || TrackPDG>499) {
                    continue; //reject stuff that's not in charm meson range
                }
                DMesonPDG->Fill(TrackPDG);
                fAllDMesonPt->Fill(AODMCtrack->Pt());
        
                if(TrackPDG==421) fD0MesonPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==421 && fromDStar==kTRUE) fD0MesonFromDStarPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==411) fDPlusMesonPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==431) fDsMesonPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==413) fDStarMesonPt->Fill(AODMCtrack->Pt());
            }
        
        }
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
    //Cluster loop//
    ////////////////
    Int_t nclus = -999;
    nclus = fAOD->GetNumberOfCaloClusters;
    for (Int_t icl = 0; icl < nclus; icl++) {
        //ESD and AOD CaloCells carries the same information
        AliVCluster* clus = (AliAODCaloCluster*)fAOD->GetCaloCluster(icl);
        if(clus && clus->IsEMCAL()){
            fClsEAll->Fill(clus->E()); //E of all clusters
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
        
        Double_t nsigma = -999;
        nsigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        if(nsigma>-5.&&nsigma<-3.) {
            fHadronCamDCA->Fill(track->Pt(),d0z0[0]);
        }
        
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
                //fClsEamEMCal->Fill(clustMatch->E());
            }
            if(emcphi > 4.53 && emcphi < 5.708) {
                fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327
                //fClsEamDCal->Fill(clustMatch->E());
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
            Double_t dWeight = -99;
            Bool_t kEmbEta = kFALSE;
            Bool_t kEmbPi0 = kFALSE;
            Bool_t kHijing = kFALSE;
            Bool_t kFlagReco = kFALSE;
            Int_t fMomGen = 99;
            Double_t momPt = -99;
            Int_t pdg = -99;
            Int_t pidGM = -99;
            Int_t ilabel = -99;
            Int_t ilabelM = -99;
            Int_t ilabelGM = -99;
            ilabel = TMath::Abs(track->GetLabel()); //get MC label of track
            if(ilabel == 0) continue;
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
                    Double_t tempValue[5] = {-999,-999,-999,-999,-999};
                    tempValue[0] = track->Pt();
                    tempValue[1] = DCA;
                    tempValue[2] = fpidSort;
                    tempValue[3] = fMomGen;
                    tempValue[4] = momPt;
                    
                    fSprsTemplatesNoWeight->Fill(tempValue);
                    
                    //Took out Lambda_c (fpidSort==17) to the weighting
                    if (fpidSort==2||fpidSort==11||fpidSort==12||fpidSort==14||fpidSort==15||fpidSort==16) {
                        if (momPt>1 && momPt<50.) {
                            dWeight = fDWeight->GetBinContent(fDWeight->FindBin(momPt));
                            fDTemplateWeight->Fill(track->Pt(), DCA, dWeight);
                            fDTemplateNoWeight->Fill(track->Pt(), DCA);
                            fSprsTemplatesWeight->Fill(tempValue,dWeight);
                        }
                    }
                    else{
                        fSprsTemplatesWeight->Fill(tempValue);
                    }
        
                    //if electron from D0
                    if(fpidSort==12){
                        ilabelM = fMCparticle->GetMother(); //get MC label for e Mom
                        AliAODMCParticle *MCpartM = (AliAODMCParticle*)fMCarray->At(ilabelM); //get e mom particle
                        ilabelGM = MCpartM->GetMother(); //get MC label for grandma
                        //if no grandma.............
                        if (ilabelGM<0) {
                            fPromptD0DCANoWeight->Fill(track->Pt(),DCA);
                            if (momPt>1 && momPt<50.) {
                                dWeight = fDWeight->GetBinContent(fDWeight->FindBin(momPt));
                                fPromptD0DCAWeight->Fill(track->Pt(),DCA,dWeight);
                            }
                        }
                        //if grandma.................
                        if (ilabelGM>0) {
                            AliAODMCParticle *MCpartGM = (AliAODMCParticle*)fMCarray->At(ilabelGM); //get grandma particle
                            pidGM = TMath::Abs(MCpartGM->GetPdgCode()); //ask for grandma's pid
                            //if grandma is D*+..............
                            if (pidGM==413) {
                                fD0FromDStarDCANoWeight->Fill(track->Pt(),DCA);
                                if (momPt>1 && momPt<50.) {
                                    dWeight = fDWeight->GetBinContent(fDWeight->FindBin(momPt));
                                    fD0FromDStarDCAWeight->Fill(track->Pt(),DCA,dWeight);
                                }
                            }
                        }
                    }
                    
                }else{continue;}
                
            }
            
            //////////////////////
            // Get MC True DCAs //
            //////////////////////
            //if MC--------------------------------
            if(ilabel>0 && fMCarray)
            {
                //cout<<"TESTING1"<<endl;
                
                //if electron--------------------------------
                if(TMath::Abs(pdg)==11){
                    
                    //if mom is Pi0--------------------------------
                    if(fpidSort==3) {
                        fPi0DCA->Fill(track->Pt(),DCA);
                        if(!kEmbEta && !kEmbEta && kHijing) {
                            fPi0HijingDCA->Fill(track->Pt(),DCA);
                            fPi0HijingPt->Fill(track->Pt());
                            //ComboDenomWeight->Fill(track->Pt());
                            //ComboDenomNoWeight->Fill(track->Pt());
                        }
                        if(kEmbPi0) {
                            fWeight = fPi0Weight->Eval(momPt);
                            fEnhPi0DCA->Fill(track->Pt(),DCA);
                            fEnhPi0WeightedPt->Fill(track->Pt(),fWeight);
                            ComboDenomWeight->Fill(track->Pt(),fWeight);
                            ComboDenomNoWeight->Fill(track->Pt());
                        }
                        if(kEmbEta) {
                            fWeight = fEtaWeight->Eval(momPt);
                            fEnhEtaDCA->Fill(track->Pt(),DCA);
                            fEnhEtaWeightedPt->Fill(track->Pt(),fWeight);
                            ComboDenomWeight->Fill(track->Pt(),fWeight);
                            ComboDenomNoWeight->Fill(track->Pt());
                        }
                        
                    }
                    //if Eta--------------------------------
                    if(fpidSort==4){
                        fEtaDCA->Fill(track->Pt(),DCA);
                        if(!kEmbEta && !kEmbPi0 && kHijing) {
                            fEtaHijingDCA->Fill(track->Pt(),DCA);
                            fEtaHijingPt->Fill(track->Pt());
                            //ComboDenomWeight->Fill(track->Pt());
                            //ComboDenomNoWeight->Fill(track->Pt());
                        }
                        if(kEmbEta) {
                            fWeight = fEtaWeight->Eval(momPt);
                            fEnhEtaDCA->Fill(track->Pt(),DCA);
                            fEnhEtaWeightedPt->Fill(track->Pt(),fWeight);
                            ComboDenomWeight->Fill(track->Pt(),fWeight);
                            ComboDenomNoWeight->Fill(track->Pt());
                        }
                    }
                    //if photon--------------------------------
                    if(fpidSort==5){
                        if(!kEmbEta && !kEmbPi0 && kHijing) {
                            fPhotonHijingDCA->Fill(track->Pt(),DCA);
                            fPhotonHijingPt->Fill(track->Pt());
                            //ComboDenomWeight->Fill(track->Pt());
                            //ComboDenomNoWeight->Fill(track->Pt());
                        }
                        if(kEmbPi0) {
                            fWeight = fPi0Weight->Eval(momPt);
                            fEnhPhotonDCA->Fill(track->Pt(),DCA);
                            fEnhPhotonWeightedPt->Fill(track->Pt(),fWeight);
                            ComboDenomWeight->Fill(track->Pt(),fWeight);
                            ComboDenomNoWeight->Fill(track->Pt());
                        }
                        if(kEmbEta) {
                            fWeight = fEtaWeight->Eval(momPt);
                            fEnhPhotonDCA->Fill(track->Pt(),DCA);
                            fEnhPhotonWeightedPt->Fill(track->Pt(),fWeight);
                            ComboDenomWeight->Fill(track->Pt(),fWeight);
                            ComboDenomNoWeight->Fill(track->Pt());
                        }
                    }

                    
                    InvMassCheckMC(i, track, d0z0, MagSign, kHijing, kEmbEta, kEmbPi0, kFlagReco,fWeight, fpidSort);
                
                    //cout<<"TESTING2"<<endl;
                    if(kFlagReco){
                        
                        //if mom is Pi0--------------------------------
                        if(fpidSort==3) {
                            fPi0DCA->Fill(track->Pt(),DCA);
                            if(!kEmbEta && !kEmbEta && kHijing) {
                                fULSHijingPi0->Fill(track->Pt()); //pi0 mama
                                //ComboNumWeight->Fill(track->Pt());
                                //ComboNumNoWeight->Fill(track->Pt());
                            }
                            if(kEmbPi0) {
                                fWeight = fPi0Weight->Eval(momPt);
                                fULSWeightEnhPi0->Fill(track->Pt(),fWeight); //pi0 mama
                                ComboNumWeight->Fill(track->Pt(),fWeight);
                                ComboNumNoWeight->Fill(track->Pt());
                            }
                            if(kEmbEta) {
                                fWeight = fEtaWeight->Eval(momPt);
                                fULSWeightEnhEta->Fill(track->Pt(),fWeight); //eta mama
                                ComboNumWeight->Fill(track->Pt(),fWeight);
                                ComboNumNoWeight->Fill(track->Pt());
                            }
                        
                        }
                        //if Eta--------------------------------
                        if(fpidSort==4){
                            fEtaDCA->Fill(track->Pt(),DCA);
                            if(!kEmbEta && !kEmbPi0 && kHijing) {
                                fULSHijingEta->Fill(track->Pt()); //eta mama
                                //ComboNumWeight->Fill(track->Pt());
                                //ComboNumNoWeight->Fill(track->Pt());
                            }
                            if(kEmbEta) {
                                fWeight = fEtaWeight->Eval(momPt);
                                fULSWeightEnhEta->Fill(track->Pt(),fWeight); //eta mama
                                ComboNumWeight->Fill(track->Pt(),fWeight);
                                ComboNumNoWeight->Fill(track->Pt());
                            }
                        }
                        //if photon--------------------------------
                        if(fpidSort==5){
                            if(!kEmbEta && !kEmbPi0 && kHijing) {
                                fULSHijingPhoton->Fill(track->Pt()); //photon mama
                                //ComboNumWeight->Fill(track->Pt());
                                //ComboNumNoWeight->Fill(track->Pt());
                            }
                            if(kEmbPi0) {
                                fWeight = fPi0Weight->Eval(momPt);
                                fULSEnhPhoton->Fill(track->Pt(),fWeight); //photon mama
                                ComboNumWeight->Fill(track->Pt(),fWeight);
                                ComboNumNoWeight->Fill(track->Pt());
                            }
                            if(kEmbEta) {
                                fWeight = fEtaWeight->Eval(momPt);
                                fULSEnhPhoton->Fill(track->Pt(),fWeight); //photon mama
                                ComboNumWeight->Fill(track->Pt(),fWeight);
                                ComboNumNoWeight->Fill(track->Pt());
                            }
                        }
                    }
                }
                
            }
            
            if(!fMCarray){
                InvMassCheckData(i, track, d0z0, MagSign);
            }
            
            ///////////////////////////
            // Hadron Contam. Histos //
            ///////////////////////////
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
            //if(nsigma>-5.&&nsigma<-3.) {
            //    fHadronCamDCA->Fill(track->Pt(),d0z0[0]*MagSign);
            //}
            
            ///////////////////
            // Electron Cuts //
            ///////////////////
            if((nsigma<-1) || (nsigma>3)) continue;
            if((M20<0.01) || (M20>0.35)) continue;
            fInclElecEoP->Fill(track->Pt(),EovP);
            
            if((EovP<0.9) || (EovP>1.2)) continue;
            
            if(fClsTypeDCAL) fClsEamElecDC->Fill(clustMatch->E());
            if(fClsTypeEMC) fClsEamElecEMC->Fill(clustMatch->E());
            
            /////////////////////////
            // Plot Reco Electrons //
            /////////////////////////
            fInclElecDCA->Fill(track->Pt(),DCA);
            fInclElecDCAnoSign->Fill(track->Pt(),d0z0[0]);
            
            
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
    Double_t fvalue[4] = {-999,-999,-999,-999};
    
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
    Int_t ilabelGGGM = -1;
    
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
            fpidSort = 18; //Mom is something else
        }
        //looking for specific particles in the ranges
        if(pidM==411){
            fpidSort = 11; //Mom is D+
        }
        else if(pidM==421){
            fpidSort = 12; //Mom is D0
        }
        else if(pidM==413){
            fpidSort = 14; //Mom is D*+
        }
        else if(pidM==431){
            fpidSort = 15; //Ds
        }
        else if(pidM>430 && pidM<436){
            fpidSort = 16; //other Ds
        }else if(pidM==4122){
            fpidSort = 17; //Lambda c
        }else if(pidM==443){
            fpidSort = 6; //Mom is J/psi
        }
        
        if(ilabelGM>0){
            AliAODMCParticle *partGM = (AliAODMCParticle*)fMCarray->At(ilabelGM); // get GMa particle
            Int_t pidGM = TMath::Abs(partGM->GetPdgCode()); //ask for grandma's pid
            ilabelGGM = partGM->GetMother();//get MC for Great Grandma
            
            //check if pi0 grandma is eta
            if(pidM==111){
                if(pidGM==221){
                    fpidSort = 4; //label as eta electron
                    if(ilabelGM >= fNembMCeta && ilabelGM < fNtotMCpart && ilabelGGM<0) { //GMa is eta
                        kEmbPi0 = kFALSE;
                        kEmbEta = kTRUE;
                        momPt = partGM->Pt(); //make eta pt mompt for weighting
                    }
                }
            }
            //check if gamma grandma is eta/pion
            if(pidM==22){
                if(pidGM==221){
                    if(ilabelGM >= fNembMCeta && ilabelGM<fNtotMCpart && ilabelGGM<0) {
                        kEmbPi0 = kFALSE;
                        kEmbEta = kTRUE; //GMa is enh eta
                        momPt = partGM->Pt(); //make eta pt mompt for weighting
                    }
                }
                if(pidGM==111){
                    if(ilabelM >= fNembMCpi0 && ilabelM<fNembMCeta && ilabelGGM<0) { //GMa is emb pi0
                        kEmbEta = kFALSE;
                        kEmbPi0 = kTRUE;
                        momPt = partGM->Pt(); //make pi0 pt mompt for weighting
                    }//GMa is pi0
                }
            }
            
            
            //check if D grandma is B
            if(pidM>400 && pidM<499){
                if(pidGM>500 && pidGM<599){
                    fpidSort = 1; //GMa is B
                }
                if(pidGM>5000 && pidGM<5999){
                    fpidSort = 10; //GMa is b baryon
                }
            }
            //check if charm baryon grandma is B
            if(pidM>4000 && pidM<4999){
                if(pidGM>500 && pidGM<599){
                    fpidSort = 1; //GMa is B
                }
                if(pidGM>5000 && pidGM<5999){
                    fpidSort = 10; //GMa is b baryon
                }
            }
            if(ilabelGGM>0){
                AliAODMCParticle *partGGM = (AliAODMCParticle*)fMCarray->At(ilabelGGM); // get GGMa particle
                Int_t pidGGM = TMath::Abs(partGGM->GetPdgCode()); //ask for ggma's pid
                ilabelGGGM = partGGM->GetMother();//get MC for Great Grandma
                
                //check if D great grandma is B
                if(pidM>400 && pidM<499){
                    if(pidGGM>500 && pidGGM<599){
                        fpidSort = 1; //GGMa is B
                    }
                    if(pidGGM>5000 && pidGGM<5999){
                        fpidSort = 10; //GGMa is b baryon
                    }
                }
                //check if charm baryon great grandma is B
                if(pidM>4000 && pidM<4999){
                    if(pidGGM>500 && pidGGM<599){
                        fpidSort = 1; //GGMa is B
                    }
                    if(pidGGM>5000 && pidGGM<5999){
                        fpidSort = 10; //GGMa is b baryon
                    }
                }
                
                //check if gamma great grandma is eta
                if(pidM==22){
                    if(pidGM==111){ //grandma is pion
                        if(pidGGM==221){ //great grandma is eta
                            if(ilabelGGM >= fNembMCeta && ilabelGGM<fNtotMCpart && ilabelGGGM<0) {
                                kEmbEta = kTRUE; //GMa is enh eta
                                kEmbPi0 = kFALSE;
                                momPt = partGGM->Pt(); //make eta pt mompt for weighting
                            }
                        }
                    }
                }
            }
        }
    }else{
        fpidSort = 19; //No mother
    }
}
//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::InvMassCheckData(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign)
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
    
    //fPhotonicElecYield->Fill(track->Pt(),Nuls-Nls);
}
//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::InvMassCheckMC(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign, Bool_t kHijing, Bool_t kEmbEta, Bool_t kEmbPi0, Bool_t &kFlagReco, Double_t fWeight, Int_t fpidSort)
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
        if(fFlagLS && track->Pt()>1){
            fInvmassLS->Fill(mass);
        }
        if(fFlagULS && track->Pt()>1){
            fInvmassULS->Fill(mass);
        }
        
        if(fFlagLS && mass<0.1) Nls++;
        if(fFlagULS && mass<0.1) Nuls++;
        
        //CHANGED FROM pt>1
        if (fFlagULS && mass<0.1) {
            kFlagReco = kTRUE;
            /*fULSdcaBelow->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign);
            if (kHijing) {
                if (fpidSort==5){
                    fULSHijingPhoton->Fill(track->Pt()); //photon mama
                    ComboNumWeight->Fill(track->Pt());
                    ComboNumNoWeight->Fill(track->Pt());
                }
                if (fpidSort==4){
                    fULSHijingEta->Fill(track->Pt()); //eta mama
                    ComboNumWeight->Fill(track->Pt());
                    ComboNumNoWeight->Fill(track->Pt());
                }
                if (fpidSort==3) {
                    fULSHijingPi0->Fill(track->Pt()); //pi0 mama
                    ComboNumWeight->Fill(track->Pt());
                    ComboNumNoWeight->Fill(track->Pt());
                }
            }
            if (kEmbEta) {
                if (fpidSort==4) {
                    fULSWeightEnhEta->Fill(track->Pt(),fWeight); //eta mama
                    ComboNumWeight->Fill(track->Pt(),fWeight);
                    ComboNumNoWeight->Fill(track->Pt());
                }
                if (fpidSort==5) {
                    fULSEnhPhoton->Fill(track->Pt(),fWeight); //photon mama
                    ComboNumWeight->Fill(track->Pt(),fWeight);
                    ComboNumNoWeight->Fill(track->Pt());
                }
            }
            if (kEmbPi0) {
                if (fpidSort==3) {
                    fULSWeightEnhPi0->Fill(track->Pt(),fWeight); //pi0 mama
                    ComboNumWeight->Fill(track->Pt(),fWeight);
                    ComboNumNoWeight->Fill(track->Pt());
                }
                if (fpidSort==5) {
                    fULSEnhPhoton->Fill(track->Pt(),fWeight); //photon mama
                    ComboNumWeight->Fill(track->Pt(),fWeight);
                    ComboNumNoWeight->Fill(track->Pt());
                }
            }*/
            
        }else if(fFlagLS && mass<0.1){
            kFlagReco = kFALSE;
            /*fLSdcaBelow->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign);
            if (kHijing) {
                if (fpidSort==5) {
                    fLSHijingPhoton->Fill(track->Pt()); //photon mama
                    //ComboNumWeight->Fill(track->Pt());
                    //ComboNumNoWeight->Fill(track->Pt());
                }
                if (fpidSort==4) {
                    fLSHijingEta->Fill(track->Pt()); //eta mama
                    //ComboNumWeight->Fill(track->Pt());
                    //ComboNumNoWeight->Fill(track->Pt());
                }
                if (fpidSort==3) {
                    fLSHijingPi0->Fill(track->Pt()); //pi0 mama
                    //ComboNumWeight->Fill(track->Pt());
                    //ComboNumNoWeight->Fill(track->Pt());
                }
            }
            if (kEmbEta) {
                if (fpidSort==4) {
                    fLSWeightEnhEta->Fill(track->Pt(),fWeight); //eta mama
                    //ComboNumWeight->Fill(track->Pt(),fWeight);
                    //ComboNumNoWeight->Fill(track->Pt());
                }
                if (fpidSort==5) {
                    fLSEnhPhoton->Fill(track->Pt(),fWeight); //photon mama
                    //ComboNumWeight->Fill(track->Pt(),fWeight);
                    //ComboNumNoWeight->Fill(track->Pt());
                }
            }
            if (kEmbPi0) {
                if (fpidSort==3) {
                    fLSWeightEnhPi0->Fill(track->Pt(),fWeight); //pi0 mama
                    //ComboNumWeight->Fill(track->Pt(),fWeight);
                    //ComboNumNoWeight->Fill(track->Pt());
                }
                if (fpidSort==5) {
                    fLSEnhPhoton->Fill(track->Pt(),fWeight); //photon mama
                    //ComboNumWeight->Fill(track->Pt(),fWeight);
                    //ComboNumNoWeight->Fill(track->Pt());
                }
            }*/
        }
        
    }
    
}

//_____________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::Terminate(Option_t *)
{
    // terminate
}
//_____________________________________________________________________
