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
fMaxM20Cut(0),
fMinEoPCut(0),
fMinNSigCut(0),
fMinNSigAssoCut(0),
fMinPtAssoCut(0),
fDCABinSize(0),
fApplyCentrality(kTRUE),
fFlagFillSprs(kFALSE),
fFlagFillMCHistos(kFALSE),
fFlagRunStackLoop(kFALSE),
fNclusTPC(80),
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
fTrkMatch(0),
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
fnSigma(0),
fnSigmaAftTrkMatch(0),
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
fnSigaftEoPCut(0),
fnSigaftSysEoPCut(0),
fnSigaftM20EoPCut(0),
fnSigaftSysM20EoPCut(0),
fInclElecDCAnoSign(0),
fElecEoPnoSig(0),
fInclElecEoP(0),
fTPCElecEoP(0),
fHadronEoP(0),
fHadronDCA(0),
fHadronCamDCAHij(0),
fHadronCamDCA(0),
fPi0Weight(0),
fEtaWeight(0),
fDWeight(0),
fDWeightNew(0),
fDWeightVar1(0),
fDWeightVar2(0),
fBWeight(0),
fBWeightNew(0),
fBWeightVar1(0),
fBWeightVar2(0),
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
fBMesonPt(0),
fBMesonPtATLAS(0),
fBPlusPtATLAS(0),
fBMesonPtCMS(0),
fBPlusPtCMS(0),
fBMesonPtLHCb(0),
fBPlusPtLHCb(0),
fBBaryonPt(0),
fBMesonElecPt(0),
fBBaryonElecPt(0),
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
fSprsTemplatesWeightVar1(0),
fSprsTemplatesWeightVar2(0),
fDTemplateWeight(0),
fDTemplateNoWeight(0),
fDTemplateWeightNew(0),
fDTemplateWeightVar1(0),
fDTemplateWeightVar2(0),

fBTemplateWeight(0),
fBTemplateNoWeight(0),
fBTemplateWeightNew(0),
fBTemplateWeightVar1(0),
fBTemplateWeightVar2(0),

fAllElecStack(0),
fHFElecStack(0),
fBElecStack(0),

fElecTPCTrk(0),
fHFElecTPCTrk(0),
fBElecTPCTrk(0),

fElecAftTrkCuts(0),
fHFElecAftTrkCuts(0),
fBElecAftTrkCuts(0),

fElecAftTrkMatch(0),
fHFElecAftTrkMatch(0),
fBElecAftTrkMatch(0),

fElecAftTPCeID(0),
fHFElecAftTPCeID(0),
fBElecAftTPCeID(0),

fElecAftEMCeID(0),
fHFElecAftEMCeID(0),
fBElecAftEMCeID(0),

fElectronSprs(0)
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
fMaxM20Cut(0),
fMinEoPCut(0),
fMinNSigCut(0),
fMinNSigAssoCut(0),
fMinPtAssoCut(0),
fDCABinSize(0),
fApplyCentrality(kTRUE),
fFlagFillSprs(kFALSE),
fFlagFillMCHistos(kFALSE),
fFlagRunStackLoop(kFALSE),
fNclusTPC(80),
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
fTrkMatch(0),
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
fnSigma(0),
fnSigmaAftTrkMatch(0),
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
fnSigaftEoPCut(0),
fnSigaftSysEoPCut(0),
fnSigaftM20EoPCut(0),
fnSigaftSysM20EoPCut(0),
fInclElecDCAnoSign(0),
fElecEoPnoSig(0),
fInclElecEoP(0),
fTPCElecEoP(0),
fHadronEoP(0),
fHadronDCA(0),
fHadronCamDCAHij(0),
fHadronCamDCA(0),
fPi0Weight(0),
fEtaWeight(0),
fDWeight(0),
fDWeightNew(0),
fDWeightVar1(0),
fDWeightVar2(0),
fBWeight(0),
fBWeightNew(0),
fBWeightVar1(0),
fBWeightVar2(0),
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
fBMesonPt(0),
fBMesonPtATLAS(0),
fBPlusPtATLAS(0),
fBMesonPtCMS(0),
fBPlusPtCMS(0),
fBMesonPtLHCb(0),
fBPlusPtLHCb(0),
fBBaryonPt(0),
fBMesonElecPt(0),
fBBaryonElecPt(0),
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
fSprsTemplatesWeightVar1(0),
fSprsTemplatesWeightVar2(0),
fDTemplateWeight(0),
fDTemplateNoWeight(0),
fDTemplateWeightNew(0),
fDTemplateWeightVar1(0),
fDTemplateWeightVar2(0),
fBTemplateWeight(0),
fBTemplateNoWeight(0),
fBTemplateWeightNew(0),
fBTemplateWeightVar1(0),
fBTemplateWeightVar2(0),
fAllElecStack(0),
fHFElecStack(0),
fBElecStack(0),

fElecTPCTrk(0),
fHFElecTPCTrk(0),
fBElecTPCTrk(0),

fElecAftTrkCuts(0),
fHFElecAftTrkCuts(0),
fBElecAftTrkCuts(0),

fElecAftTrkMatch(0),
fHFElecAftTrkMatch(0),
fBElecAftTrkMatch(0),

fElecAftTPCeID(0),
fHFElecAftTPCeID(0),
fBElecAftTPCeID(0),

fElecAftEMCeID(0),
fHFElecAftEMCeID(0),
fBElecAftEMCeID(0),
fElectronSprs(0)
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
    /////////////////
    // Output List //
    /////////////////
    
    Int_t nDCAbins = 0.4/fDCABinSize;
    
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
    
    if(fFlagClsTypeDCAL && !fFlagClsTypeEMC){
        fClsEamDCal = new TH1F("fClsEamDCal","Cluster Energy after track matching to DCal;Cluster E;Counts",250,0.,50);
        fOutputList->Add(fClsEamDCal);
    }
    
    if(fFlagClsTypeEMC && !fFlagClsTypeDCAL){
        fClsEamEMCal = new TH1F("fClsEamEMCal","Cluster Energy after track matching to EMCal;Cluster E;Counts",250,0.,50.);
        fOutputList->Add(fClsEamEMCal);
    }
        
    fClsEAll = new TH1F("fClsEAll","Cluster Energy, All Clusters;Cluster E;Counts",250,0.,50);
    fOutputList->Add(fClsEAll);
    if(fFlagClsTypeEMC && !fFlagClsTypeDCAL){
        fClsEamElecEMC = new TH1F("fClsEamElecEMC","Cluster Energy of e- after track matching to DCal;Cluster E;Counts",250,0.,50);
        fOutputList->Add(fClsEamElecEMC);
    }
    if(fFlagClsTypeDCAL && !fFlagClsTypeEMC){
        fClsEamElecDC = new TH1F("fClsEamElecDC","Cluster Energy of e- after track matching to DCal;Cluster E;Counts",250,0.,50);
        fOutputList->Add(fClsEamElecDC);
    }
        
    fTrkPhi = new TH1F("fTrkPhi","Track #phi Distribution after matching;#phi;Counts",100,0,6.3);
    fOutputList->Add(fTrkPhi);
    
    fTrkEta = new TH1F("fTrkEta","Track #eta Distribution after matching;#eta;Counts",100,-1.5,1.5);
    fOutputList->Add(fTrkEta);
    
    fdEdx = new TH1F("fdEdx","Track dE/dx Distribution;dE/dx;Counts",160,0,160);
    fOutputList->Add(fdEdx);
    
    fnSigma = new TH2F("fnSigma","Track fnSigma Distribution;pT;fnSigma",30,0,30,100,-10,10);
    fOutputList->Add(fnSigma);
    
    fnSigmaAftTrkMatch = new TH2F("fnSigmaAftTrkMatch","Track fnSigma Distribution after track matching to cal;pT;fnSigma",30,0,30,100,-10,10);
    fOutputList->Add(fnSigmaAftTrkMatch);
    
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
    
    if (fFlagFillMCHistos) {
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
    }
    
    fULSdcaBelow = new TH2F("fULSdcaBelow","ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
    fOutputList->Add(fULSdcaBelow);
    
    fLSdcaBelow = new TH2F("fLSdcaBelow","LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
    fOutputList->Add(fLSdcaBelow);
    
    if (fFlagFillMCHistos) {
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
    
        fPhotonicDCA = new TH2F("fPhotonicDCA","Photonic DCA using MC PID; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fPhotonicDCA);
    }
    fInclElecDCA = new TH2F("fInclElecDCA","Incl Elec DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
    fOutputList->Add(fInclElecDCA);
    
    fnSigaftEoPCut = new TH2F("fnSigaftEoPCut","nSig with 0.9<E/p<1.2; p_{T}(GeV/c); nSigma; counts;", 60,0,30., 160,-8,8);
    fOutputList->Add(fnSigaftEoPCut);
    
    fnSigaftSysEoPCut = new TH2F("fnSigaftSysEoPCut","nSig with Sys E/p cut; p_{T}(GeV/c); nSigma; counts;", 60,0,30., 160,-8,8);
    fOutputList->Add(fnSigaftSysEoPCut);
    
    fnSigaftM20EoPCut = new TH2F("fnSigaftM20EoPCut","nSig with 0.9<E/p<1.2, 0.01<M20<0.35; p_{T}(GeV/c); nSigma; counts;", 60,0,30., 160,-8,8);
    fOutputList->Add(fnSigaftM20EoPCut);
    
    fnSigaftSysM20EoPCut = new TH2F("fnSigaftSysM20EoPCut","nSig with systematic E/p and M20 cut; p_{T}(GeV/c); nSigma; counts;", 60,0,30., 160,-8,8);
    fOutputList->Add(fnSigaftSysM20EoPCut);
    
    fInclElecDCAnoSign = new TH2F("fInclElecDCAnoSign","Incl Elec DCA (no Charge); p_{T}(GeV/c); DCAxMagField; counts;", 60,0,30., nDCAbins,-0.2,0.2);
    fOutputList->Add(fInclElecDCAnoSign);
    
    fElecEoPnoSig = new TH2F("fElecEoPnoSig","Elec E/p, no nSig cut; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fElecEoPnoSig);
    
    fInclElecEoP = new TH2F("fInclElecEoP","Incl Elec E/p; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fInclElecEoP);
    
    fTPCElecEoP = new TH2F("fTPCElecEoP","Elec E/p, -0.1<nsig<3; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fTPCElecEoP);
    
    fHadronEoP = new TH2F("fHadronEoP","Unscaled Hadron E/p; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fHadronEoP);
    
    fHadronDCA = new TH2F("fHadronDCA","Unscaled Hadron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
    fOutputList->Add(fHadronDCA);
    
    if (fFlagFillMCHistos) {
        fHadronCamDCAHij = new TH2F("fHadronCamDCAHij","Unscaled Hadron DCA, no E/p cut, Hijing; p_{T}(GeV/c); DCAxMagField; counts;", 60,0,30., 400,-0.2,0.2);
        fOutputList->Add(fHadronCamDCAHij);
    }
    
    fHadronCamDCA = new TH2F("fHadronCamDCA","Unscaled Hadron DCA, no E/p cut, Enh+Hij; p_{T}(GeV/c); DCAxMagField; counts;", 60,0,30., 400,-0.2,0.2);
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
    fDWeightNew = new TH1F("fDWeightNew","D^{0}_data/AllD_MCNew;p_{T} (GeV/c);Weight;",nbins,xbins);
    fDWeightVar1 = new TH1F("fDWeightVar1","D^{0}_data/AllD_MC;p_{T} (GeV/c);Weight Var1;",nbins,xbins);
    fDWeightVar2 = new TH1F("fDWeightVar2","D^{0}_data/AllD_MC;p_{T} (GeV/c);Weight Var2;",nbins,xbins);
    //Double_t ratio[13] = {2.03552,1.0201,0.45925,0.211574,0.11987,0.0898116,0.0631282,0.0546798,0.0477205,0.0410021,0.0307936,0.0398483,0.0175335};
    //Double_t err[13] = {0.541651,0.146443,0.0498454,0.024907,0.01438,0.0107908,0.00848616,0.0061723,0.00587082,0.00566712,0.00597994,0.00811015,0.00693105};
    if (fCentralityMin==0 && fCentralityMax==10 && fApplyCentrality) {
        Double_t ratio[13] = {0.106888,0.0650239,0.0343858,0.0172579,0.00957876,0.00640323,0.00399907,0.00269269,0.00163078,0.000942387,0.000441093,0.000353811,0.000143011};
        Double_t err[13] = {0.0284416,0.0093333,0.00373075,0.00203067,0.00114824,0.000768388,0.00053676,0.000303334,0.000199878,0.000129785,8.53822e-05,7.13313e-05,5.61316e-05};
        Double_t ratioNew[13] = {0.197449,0.118714,0.0627949,0.0321233,0.0182153,0.0124903,0.00801369,0.00553768,0.00340667,0.00193131,0.00089526,0.000678224,0.00026223};
        Double_t ratioVar1[13] = {0.249988,0.132914,0.0673368,0.0340132,0.0189431,0.01274,0.00801369,0.00543372,0.00326751,0.00179833,0.000779735,0.000564287,0.000159311};
        Double_t ratioVar2[13] = {0.144911,0.104514,0.058253,0.0302335,0.0174875,0.0122405,0.00801369,0.00564164,0.00354583,0.00206429,0.00101078,0.00079216,0.000365149};
        for (int idata=1; idata<14; idata++) {
            fDWeight->SetBinContent(idata,ratio[idata-1]);
            fDWeight->SetBinError(idata,err[idata-1]);
            
            fDWeightNew->SetBinContent(idata,ratioNew[idata-1]);
            fDWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
            fDWeightVar2->SetBinContent(idata,ratioVar2[idata-1]);
        }
    }else if (fCentralityMin==30 && fCentralityMax==50 && fApplyCentrality) {
        Double_t ratio[13] = {0.472316,0.277093,0.143805,0.0715906,0.0394568,0.0262827,0.0164828,0.0112023,0.00680627,0.00387925,0.00180795,0.0014053,0.000562486};
        Double_t err[13] = {0.125678,0.0397729,0.0156024,0.00842381,0.00472987,0.00315395,0.00221236,0.00126196,0.000834225,0.000534252,0.000349965,0.00028332,0.000220775};
        Double_t ratioNew[13] = {0.472316,0.277093,0.143805,0.0715906,0.0394568,0.0262827,0.0164828,0.0112023,0.00680627,0.00387925,0.00180795,0.0014053,0.000562486};
        Double_t ratioVar1[13] = {0.597993,0.310236,0.154206,0.0758023,0.0410333,0.0268083,0.0164828,0.010992,0.00652823,0.00361215,0.00157465,0.00116922,0.000341724};
        Double_t ratioVar2[13] = {0.346638,0.243949,0.133404,0.0673789,0.0378803,0.0257571,0.0164828,0.0114126,0.0070843,0.00414635,0.00204124,0.00164137,0.000783247};
        for (int idata=1; idata<14; idata++) {
            fDWeight->SetBinContent(idata,ratio[idata-1]);
            fDWeight->SetBinError(idata,err[idata-1]);
            
            fDWeightNew->SetBinContent(idata,ratioNew[idata-1]);
            fDWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
            fDWeightVar2->SetBinContent(idata,ratioVar2[idata-1]);
        }
    }else{
        Double_t ratio[13];
        Double_t err[13];
        Double_t ratioNew[13];
        Double_t ratioVar1[13];
        Double_t ratioVar2[13];
    
        for (int idata=1; idata<14; idata++) {
            fDWeight->SetBinContent(idata,ratio[idata-1]);
            fDWeight->SetBinError(idata,err[idata-1]);
            
            fDWeightNew->SetBinContent(idata,ratioNew[idata-1]);
            fDWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
            fDWeightVar2->SetBinContent(idata,ratioVar2[idata-1]);
        }
    }
    
    fDWeight->Sumw2();
    fOutputList->Add(fDWeight);
    fOutputList->Add(fDWeightNew);
    fOutputList->Add(fDWeightVar1);
    fOutputList->Add(fDWeightVar2);
    
    //B Meson pt weighting
    Int_t nbinsB = 250;
    Double_t xbinsB[251] = {0.,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8,8,8.2,8.4,8.6,8.8,9,9.2,9.4,9.6,9.8,10,10.2,10.4,10.6,10.8,11,11.2,11.4,11.6,11.8,12,12.2,12.4,12.6,12.8,13,13.2,13.4,13.6,13.8,14,14.2,14.4,14.6,14.8,15,15.2,15.4,15.6,15.8,16,16.2,16.4,16.6,16.8,17,17.2,17.4,17.6,17.8,18,18.2,18.4,18.6,18.8,19,19.2,19.4,19.6,19.8,20,20.2,20.4,20.6,20.8,21,21.2,21.4,21.6,21.8,22,22.2,22.4,22.6,22.8,23,23.2,23.4,23.6,23.8,24,24.2,24.4,24.6,24.8,25,25.2,25.4,25.6,25.8,26,26.2,26.4,26.6,26.8,27,27.2,27.4,27.6,27.8,28,28.2,28.4,28.6,28.8,29,29.2,29.4,29.6,29.8,30,30.2,30.4,30.6,30.8,31,31.2,31.4,31.6,31.8,32,32.2,32.4,32.6,32.8,33,33.2,33.4,33.6,33.8,34,34.2,34.4,34.6,34.8,35,35.2,35.4,35.6,35.8,36,36.2,36.4,36.6,36.8,37,37.2,37.4,37.6,37.8,38,38.2,38.4,38.6,38.8,39,39.2,39.4,39.6,39.8,40,40.2,40.4,40.6,40.8,41,41.2,41.4,41.6,41.8,42,42.2,42.4,42.6,42.8,43,43.2,43.4,43.6,43.8,44,44.2,44.4,44.6,44.8,45,45.2,45.4,45.6,45.8,46,46.2,46.4,46.6,46.8,47,47.2,47.4,47.6,47.8,48,48.2,48.4,48.6,48.8,49,49.2,49.4,49.6,49.8,50.};
    fBWeight = new TH1F("fBWeight","TAMU RAA x FONLL/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    fBWeightNew = new TH1F("fBWeightNew","TAMU RAA x FONLL(New)/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    fBWeightVar1 = new TH1F("fBWeightVar1","TAMU RAA(Max) x FONLL(New)/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    fBWeightVar2 = new TH1F("fBWeightVar2","TAMU RAA(Min) x FONLL(New)/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    
    if (fCentralityMin==0 && fCentralityMax==10 && fApplyCentrality) {
        Double_t ratioB[250] = {0.498558,0.62782,0.672398,0.718151,0.731677,0.734297,0.770069,0.776389,0.809054,0.825417,0.859501,0.874636,0.898427,0.923358,0.92293,0.926125,0.926343,0.916683,0.917184,0.911809,0.886661,0.868603,0.852465,0.837679,0.834915,0.813126,0.776888,0.762151,0.745226,0.726064,0.689477,0.675392,0.657326,0.639663,0.610833,0.599204,0.574208,0.551114,0.534405,0.508084,0.503657,0.468409,0.459998,0.447575,0.429903,0.413056,0.396865,0.385802,0.36568,0.358106,0.351127,0.342732,0.324592,0.325566,0.317839,0.306017,0.292589,0.292114,0.28069,0.274267,0.268707,0.260392,0.261063,0.248518,0.245039,0.237828,0.226948,0.221301,0.216985,0.214584,0.210782,0.197669,0.195025,0.190053,0.186186,0.179784,0.171765,0.169136,0.159454,0.161479,0.157118,0.153996,0.14812,0.141091,0.138405,0.137488,0.134021,0.127508,0.126735,0.119841,0.117923,0.11417,0.109935,0.110609,0.106261,0.102496,0.0988583,0.0988566,0.0962846,0.0971962,0.0912081,0.0906153,0.0909837,0.089602,0.0851854,0.0829322,0.0854181,0.0793186,0.0828916,0.079228,0.0765657,0.0755108,0.0746448,0.0749831,0.0726684,0.0741421,0.0725614,0.0725931,0.0708997,0.0691765,0.0675163,0.0638974,0.0636498,0.0696325,0.0621972,0.0654069,0.0613243,0.0604615,0.0588438,0.0602627,0.0602133,0.0562468,0.058719,0.0581855,0.0582384,0.0562759,0.053968,0.053529,0.0549016,0.0524293,0.0523908,0.0529672,0.0515757,0.0493422,0.0481627,0.0482736,0.0459697,0.0458283,0.0453141,0.0433415,0.045522,0.043535,0.0411514,0.0456961,0.0441323,0.0430352,0.0428453,0.0426922,0.0436595,0.040605,0.0381453,0.03843,0.0410627,0.0374946,0.0391381,0.0379844,0.0367569,0.036893,0.0372399,0.0355532,0.0336221,0.0344674,0.0330855,0.0332478,0.0308431,0.0318583,0.0313297,0.03175,0.0314691,0.0320693,0.0312073,0.0297707,0.0285189,0.0289505,0.0291229,0.0291655,0.0305036,0.0281665,0.0276302,0.0291551,0.026976,0.027417,0.027251,0.0245767,0.0246722,0.024308,0.0247152,0.0273125,0.0242685,0.0232144,0.0220058,0.023307,0.0224514,0.0219948,0.0208518,0.0216804,0.0215576,0.0195631,0.0195621,0.0211566,0.0185032,0.0196948,0.0190276,0.0193687,0.0191364,0.0195232,0.0190942,0.0177003,0.0186524,0.0180672,0.017813,0.0173329,0.0163168,0.0163818,0.0158271,0.0169163,0.0164239,0.0157851,0.0169332,0.0180851,0.0152819,0.0149112,0.0141721,0.015732,0.0163534,0.014337,0.0143059,0.014786,0.0140978,0.0145358,0.0145749,0.0134334,0.0130966,0.0138116,0.013262,0.0143179,0.0132258,0.0135375,0.0132521,0.0132662};
        Double_t ratioBNew[250] = {0.527854,0.676468,0.697991,0.711774,0.722505,0.744184,0.771669,0.791643,0.815676,0.842682,0.869755,0.886323,0.900165,0.921397,0.927047,0.93348,0.936192,0.926124,0.917184,0.909077,0.889903,0.876369,0.858462,0.838199,0.821877,0.805409,0.785763,0.762809,0.741829,0.72182,0.704218,0.683978,0.660144,0.639746,0.617747,0.596334,0.573373,0.553575,0.533483,0.516636,0.49255,0.475004,0.460341,0.441742,0.4269,0.412142,0.40238,0.387698,0.37121,0.360204,0.352943,0.341443,0.332603,0.323606,0.315956,0.305826,0.295788,0.291485,0.28239,0.275661,0.268983,0.260012,0.25416,0.24767,0.24152,0.234847,0.229889,0.2255,0.216294,0.21002,0.204699,0.200694,0.194504,0.189242,0.185419,0.179002,0.173614,0.169286,0.162982,0.158871,0.15605,0.149979,0.146925,0.143822,0.138876,0.135658,0.131925,0.128432,0.125932,0.121191,0.11921,0.116162,0.111859,0.11006,0.106977,0.105162,0.103074,0.100297,0.0993633,0.0962979,0.0940254,0.0921904,0.0906611,0.0891579,0.0880537,0.0863541,0.0855137,0.0817595,0.0811291,0.0799165,0.0793575,0.0781609,0.0766887,0.0750305,0.0738189,0.073269,0.072447,0.0706597,0.0707376,0.0681811,0.0680168,0.0670508,0.065782,0.0663052,0.0645379,0.0638748,0.0625307,0.0614706,0.0613937,0.0600465,0.0599109,0.0578851,0.0589849,0.0577687,0.0577089,0.0564293,0.0548461,0.0541342,0.0537008,0.0534945,0.0522089,0.0527365,0.0517209,0.0505692,0.0501077,0.0491331,0.0493484,0.0485754,0.0477572,0.0471697,0.0461767,0.0455698,0.0463311,0.0454191,0.0441808,0.0445879,0.0432137,0.043078,0.0429534,0.041554,0.0404588,0.040246,0.0399505,0.0400038,0.0393302,0.0378221,0.0377077,0.0375886,0.0364395,0.0361215,0.0359106,0.0354513,0.035101,0.0343332,0.0340229,0.033936,0.0329639,0.0324529,0.0320335,0.0311993,0.0318253,0.0314773,0.0295871,0.0297416,0.0294255,0.0293431,0.0287639,0.028382,0.0275853,0.0278151,0.0270327,0.026975,0.0267576,0.0258747,0.0254334,0.0253068,0.0250752,0.0246926,0.0236031,0.0242901,0.0232739,0.0232222,0.0229169,0.0232431,0.0216698,0.0220481,0.0212689,0.0212692,0.0212703,0.0207071,0.0207659,0.0199056,0.019779,0.0194439,0.0190936,0.0189077,0.0188995,0.017954,0.0178651,0.0179176,0.0177795,0.0169843,0.0170956,0.0171893,0.0168782,0.0165971,0.017014,0.0161936,0.0161965,0.0158903,0.0157885,0.0152633,0.0151995,0.0154153,0.0155522,0.0153605,0.0147381,0.0147705,0.0144276,0.0145398,0.0142471,0.0144064,0.0138917,0.0140977,0.014126,0.0138707,0.0134346,0.0137999,0.013464,0.0130542};
        Double_t ratioBVar1[250] = {0.54281,0.695606,0.717742,0.731959,0.743084,0.765527,0.794016,0.814866,0.839998,0.86832,0.896867,0.914756,0.930023,0.953155,0.96042,0.968762,0.97354,0.965327,0.958592,0.953066,0.936273,0.925755,0.910985,0.894068,0.881733,0.869652,0.854532,0.83615,0.820235,0.805694,0.794136,0.779845,0.761552,0.747232,0.730978,0.715231,0.697306,0.682818,0.667487,0.655676,0.633963,0.619851,0.608774,0.591691,0.578793,0.565203,0.557723,0.542693,0.524328,0.512981,0.506383,0.493153,0.483236,0.472629,0.463576,0.450508,0.437227,0.432141,0.419709,0.410571,0.401326,0.388495,0.380185,0.370806,0.36184,0.352007,0.344679,0.338147,0.324347,0.314906,0.306865,0.300774,0.291389,0.283384,0.277522,0.267772,0.259559,0.252931,0.243351,0.237051,0.232676,0.22346,0.218746,0.213963,0.206445,0.201502,0.195801,0.190463,0.186603,0.179431,0.176352,0.171701,0.165204,0.162412,0.15773,0.154925,0.151721,0.147509,0.146013,0.14139,0.137937,0.135132,0.132778,0.130467,0.128743,0.126153,0.12482,0.119241,0.118222,0.116358,0.115448,0.113612,0.11138,0.108881,0.107034,0.106149,0.104871,0.1022,0.102228,0.0984524,0.0981344,0.0966614,0.0947546,0.0954302,0.0928108,0.0917825,0.0897781,0.0881844,0.0880027,0.0860019,0.0857383,0.0827723,0.0842771,0.0824731,0.0823216,0.0804317,0.0781125,0.0770371,0.0763594,0.0760055,0.07412,0.0748097,0.0733108,0.0716217,0.070912,0.0694779,0.0697275,0.0685814,0.0673733,0.0664924,0.0650417,0.0641368,0.0651575,0.0638252,0.0620368,0.0625599,0.060585,0.060348,0.0601271,0.0581234,0.0565479,0.0562073,0.0557518,0.0557834,0.0548022,0.0526606,0.0524614,0.0522559,0.05062,0.0501403,0.0498098,0.0491355,0.0486133,0.0475141,0.0470493,0.0468939,0.0455165,0.0447775,0.0441656,0.0429835,0.0438132,0.043302,0.0406715,0.0408535,0.0403895,0.0402465,0.0394231,0.038871,0.037752,0.0380386,0.0369416,0.0368358,0.0365122,0.0352817,0.0346547,0.034457,0.034117,0.0335721,0.0320675,0.0329771,0.0315746,0.0314818,0.0310455,0.0314648,0.029314,0.0298044,0.0287306,0.0287104,0.0286915,0.0279118,0.0279712,0.0267935,0.0266041,0.026135,0.025646,0.0253783,0.0253495,0.0240644,0.0239284,0.023982,0.0237805,0.0227009,0.0228338,0.022943,0.0225121,0.0221218,0.0226618,0.0215541,0.0215431,0.0211211,0.0209714,0.0202598,0.0201613,0.0204336,0.0206009,0.0203331,0.0194958,0.0195254,0.0190591,0.0191943,0.0187951,0.0189924,0.0183014,0.0185603,0.018585,0.0182369,0.0176516,0.0181193,0.0176665,0.0171173};
        Double_t ratioBVar2[250] = {0.512898,0.65733,0.67824,0.69159,0.701926,0.722841,0.749322,0.768421,0.791354,0.817044,0.842643,0.85789,0.870306,0.889639,0.893675,0.898198,0.898843,0.886921,0.875776,0.865088,0.843533,0.826984,0.80594,0.782329,0.762022,0.741167,0.716994,0.689468,0.663422,0.637946,0.6143,0.588111,0.558735,0.53226,0.504516,0.477438,0.449439,0.424332,0.399479,0.377596,0.351136,0.330158,0.311908,0.291794,0.275008,0.259082,0.247036,0.232703,0.218092,0.207428,0.199503,0.189734,0.181971,0.174584,0.168336,0.161143,0.154349,0.150828,0.145071,0.14075,0.13664,0.131529,0.128135,0.124534,0.121201,0.117687,0.1151,0.112853,0.108241,0.105134,0.102532,0.100614,0.0976182,0.0951008,0.0933164,0.0902324,0.0876688,0.0856413,0.0826121,0.080691,0.0794234,0.0764967,0.0751029,0.0736805,0.0713075,0.0698141,0.0680495,0.0664012,0.0652605,0.0629505,0.0620669,0.0606222,0.0585142,0.0577086,0.0562239,0.0554,0.0544271,0.0530847,0.0527133,0.051206,0.0501136,0.0492492,0.0485438,0.0478486,0.0473641,0.0465557,0.0462072,0.0442785,0.0440359,0.0434749,0.0432671,0.0427093,0.0419975,0.0411799,0.0406037,0.0403892,0.0400228,0.0391198,0.0392472,0.0379098,0.0378991,0.0374402,0.0368093,0.0371801,0.0362649,0.0359672,0.0352834,0.0347569,0.0347848,0.0340911,0.0340835,0.0329978,0.0336927,0.0330643,0.0330962,0.0324268,0.0315796,0.0312313,0.0310422,0.0309835,0.0302978,0.0306634,0.030131,0.0295167,0.0293034,0.0287882,0.0289693,0.0285695,0.0281411,0.027847,0.0273117,0.0270028,0.0275047,0.027013,0.0263247,0.0266158,0.0258425,0.025808,0.0257798,0.0249847,0.0243697,0.0242847,0.0241492,0.0242241,0.0238581,0.0229835,0.022954,0.0229212,0.022259,0.0221028,0.0220115,0.021767,0.0215887,0.0211522,0.0209965,0.0209781,0.0204113,0.0201284,0.0199013,0.0194151,0.0198373,0.0196527,0.0185027,0.0186296,0.0184616,0.0184396,0.0181047,0.017893,0.0174185,0.0175915,0.0171238,0.0171143,0.017003,0.0164678,0.0162122,0.0161565,0.0160335,0.0158132,0.0151387,0.0156032,0.0149732,0.0149626,0.0147882,0.0150214,0.0140256,0.0142919,0.0138073,0.013828,0.0138492,0.0135023,0.0135606,0.0130178,0.0129538,0.0127529,0.0125413,0.012437,0.0124495,0.0118436,0.0118018,0.0118533,0.0117786,0.0112676,0.0113574,0.0114356,0.0112443,0.0110724,0.0113663,0.0108331,0.01085,0.0106594,0.0106056,0.0102667,0.0102377,0.0103971,0.0105035,0.010388,0.00998035,0.0100156,0.00979607,0.0098853,0.00969908,0.0098204,0.00948191,0.00963511,0.00966696,0.00950457,0.0092176,0.0094804,0.00926154,0.00899106};
        for (int idata=1; idata<251; idata++) {
            fBWeight->SetBinContent(idata,ratioB[idata-1]);
            fBWeightNew->SetBinContent(idata,ratioBNew[idata-1]);
            fBWeightVar1->SetBinContent(idata,ratioBVar1[idata-1]);
            fBWeightVar2->SetBinContent(idata,ratioBVar2[idata-1]);
        }
    }else if (fCentralityMin==30 && fCentralityMax==50 && fApplyCentrality) {
        Double_t ratioB[250] = {15.776,20.1256,20.8146,21.5019,21.3878,22.3697,22.8105,23.6098,24.2267,24.9541,25.6441,26.3561,26.9386,27.3885,27.6019,27.7579,27.7871,27.5959,27.2793,27.0913,26.5028,25.9281,25.4498,24.8388,24.4589,23.9615,23.4909,22.769,22.123,21.4844,20.8086,20.2888,19.5455,19.0463,18.2068,17.5922,17.0383,16.3843,15.6861,15.2584,14.798,13.9779,13.4579,13.1398,12.6415,12.2117,11.8548,11.3553,11.0465,10.7081,10.3911,10.1305,9.92061,9.56648,9.25112,9.09984,8.91987,8.44315,8.34427,8.01732,7.9745,7.66226,7.56115,7.3246,7.0889,6.87951,6.71628,6.49561,6.2974,6.19828,6.01367,5.90582,5.68053,5.51326,5.3825,5.23421,5.18718,4.99365,4.80072,4.65244,4.49083,4.3721,4.29031,4.16016,3.99989,3.97594,3.86849,3.7549,3.6126,3.60325,3.44631,3.38378,3.33835,3.24526,3.17999,3.05991,3.00964,2.90524,2.86489,2.79226,2.7782,2.66336,2.59954,2.55306,2.52766,2.47114,2.43289,2.36425,2.36679,2.33415,2.33202,2.24345,2.20614,2.19795,2.15301,2.1011,2.09064,2.09047,2.08161,2.0386,2.00781,1.95784,1.90601,1.89256,1.83627,1.77392,1.806,1.78043,1.76157,1.7597,1.68611,1.69052,1.6739,1.65687,1.61229,1.58701,1.56778,1.5686,1.52242,1.52299,1.50742,1.50004,1.46443,1.46295,1.47865,1.43253,1.42391,1.37535,1.37515,1.36354,1.32573,1.29165,1.28718,1.26891,1.23962,1.25188,1.24054,1.21876,1.23364,1.21827,1.19191,1.18391,1.12391,1.14973,1.07696,1.08079,1.10159,1.06866,1.05734,1.04826,1.04148,1.0283,0.977329,0.97373,0.970759,0.944777,0.92652,0.920294,0.906087,0.879841,0.880005,0.86038,0.85176,0.887516,0.842141,0.841606,0.783892,0.800277,0.777027,0.787906,0.743009,0.759993,0.749196,0.721056,0.722033,0.692535,0.697764,0.676892,0.68242,0.65176,0.653058,0.663174,0.652346,0.619086,0.623192,0.601236,0.619329,0.600425,0.603537,0.596844,0.591831,0.557848,0.56453,0.54059,0.538191,0.513374,0.538588,0.498087,0.49361,0.526396,0.480338,0.481123,0.475169,0.473394,0.456706,0.47354,0.466587,0.483567,0.463179,0.445174,0.449093,0.420737,0.42934,0.445405,0.425301,0.423512,0.432645,0.424743,0.407295,0.415857,0.401262,0.401077,0.391331,0.378896,0.376927,0.372601,0.373813,0.369854,0.392131,0.363538};
        Double_t ratioBNew[250] = {15.776,20.1256,20.8146,21.5019,21.3878,22.3697,22.8105,23.6098,24.2267,24.9541,25.6441,26.3561,26.9386,27.3885,27.6019,27.7579,27.7871,27.5959,27.2793,27.0913,26.5028,25.9281,25.4498,24.8388,24.4589,23.9615,23.4909,22.769,22.123,21.4844,20.8086,20.2888,19.5455,19.0463,18.2068,17.5922,17.0383,16.3843,15.6861,15.2584,14.798,13.9779,13.4579,13.1398,12.6415,12.2117,11.8548,11.3553,11.0465,10.7081,10.3911,10.1305,9.92061,9.56648,9.25112,9.09984,8.91987,8.44315,8.34427,8.01732,7.9745,7.66226,7.56115,7.3246,7.0889,6.87951,6.71628,6.49561,6.2974,6.19828,6.01367,5.90582,5.68053,5.51326,5.3825,5.23421,5.18718,4.99365,4.80072,4.65244,4.49083,4.3721,4.29031,4.16016,3.99989,3.97594,3.86849,3.7549,3.6126,3.60325,3.44631,3.38378,3.33835,3.24526,3.17999,3.05991,3.00964,2.90524,2.86489,2.79226,2.7782,2.66336,2.59954,2.55306,2.52766,2.47114,2.43289,2.36425,2.36679,2.33415,2.33202,2.24345,2.20614,2.19795,2.15301,2.1011,2.09064,2.09047,2.08161,2.0386,2.00781,1.95784,1.90601,1.89256,1.83627,1.77392,1.806,1.78043,1.76157,1.7597,1.68611,1.69052,1.6739,1.65687,1.61229,1.58701,1.56778,1.5686,1.52242,1.52299,1.50742,1.50004,1.46443,1.46295,1.47865,1.43253,1.42391,1.37535,1.37515,1.36354,1.32573,1.29165,1.28718,1.26891,1.23962,1.25188,1.24054,1.21876,1.23364,1.21827,1.19191,1.18391,1.12391,1.14973,1.07696,1.08079,1.10159,1.06866,1.05734,1.04826,1.04148,1.0283,0.977329,0.97373,0.970759,0.944777,0.92652,0.920294,0.906087,0.879841,0.880005,0.86038,0.85176,0.887516,0.842141,0.841606,0.783892,0.800277,0.777027,0.787906,0.743009,0.759993,0.749196,0.721056,0.722033,0.692535,0.697764,0.676892,0.68242,0.65176,0.653058,0.663174,0.652346,0.619086,0.623192,0.601236,0.619329,0.600425,0.603537,0.596844,0.591831,0.557848,0.56453,0.54059,0.538191,0.513374,0.538588,0.498087,0.49361,0.526396,0.480338,0.481123,0.475169,0.473394,0.456706,0.47354,0.466587,0.483567,0.463179,0.445174,0.449093,0.420737,0.42934,0.445405,0.425301,0.423512,0.432645,0.424743,0.407295,0.415857,0.401262,0.401077,0.391331,0.378896,0.376927,0.372601,0.373813,0.369854,0.392131,0.363538};
        Double_t ratioBVar1[250] = {16.223,20.695,21.4036,22.1116,21.997,23.0112,23.4711,24.3024,24.9491,25.7134,26.4434,27.2016,27.8322,28.3325,28.5955,28.8071,28.8957,28.764,28.5108,28.4022,27.8838,27.3892,27.0069,26.4944,26.2402,25.8727,25.5468,24.9582,24.4612,23.9809,23.4655,23.1325,22.548,22.2463,21.5441,21.0997,20.7211,20.2096,19.6262,19.3648,19.0466,18.2403,17.7972,17.6001,17.1394,16.7468,16.4316,15.8949,15.603,15.2499,14.9085,14.6317,14.4136,13.9719,13.5734,13.4049,13.1851,12.5174,12.4019,11.9411,11.8981,11.4485,11.3103,10.9662,10.6204,10.3116,10.0699,9.74045,9.44335,9.29376,9.01513,8.85086,8.51009,8.25591,8.05614,7.82993,7.75502,7.46102,7.16805,6.94189,6.696,6.5142,6.38756,6.18906,5.946,5.90573,5.74155,5.56846,5.35308,5.33485,5.09829,5.00164,4.93039,4.78891,4.68867,4.50785,4.43007,4.27281,4.20992,4.09975,4.07568,3.90392,3.80717,3.73596,3.69569,3.61003,3.55117,3.44809,3.44891,3.39851,3.39258,3.26101,3.20411,3.18957,3.12177,3.04398,3.02632,3.02358,3.00829,2.9437,2.89686,2.82245,2.74549,2.72388,2.6407,2.54897,2.59295,2.55416,2.52505,2.52033,2.41299,2.41735,2.39165,2.36542,2.29992,2.26205,2.23286,2.23224,2.16479,2.16389,2.14006,2.12789,2.07573,2.07199,2.09258,2.0257,2.01194,1.94179,1.93999,1.9221,1.86734,1.81792,1.81022,1.78314,1.74062,1.75648,1.73922,1.70736,1.72687,1.70404,1.66589,1.65344,1.56844,1.60325,1.50063,1.50481,1.53261,1.48566,1.46881,1.45508,1.44458,1.42522,1.35356,1.34756,1.34244,1.30553,1.27934,1.26979,1.24925,1.21216,1.21149,1.18359,1.17086,1.21911,1.15592,1.15434,1.07438,1.09603,1.06341,1.0775,1.01536,1.03781,1.02232,0.983202,0.983816,0.942938,0.949368,0.920302,0.927146,0.884851,0.885974,0.89905,0.883734,0.838074,0.843028,0.812745,0.836603,0.810489,0.814109,0.804508,0.797184,0.750877,0.759332,0.726617,0.722882,0.689063,0.722397,0.667603,0.661139,0.704558,0.642462,0.643063,0.634661,0.63185,0.609153,0.631167,0.621469,0.64364,0.616076,0.591719,0.596518,0.558468,0.569497,0.5904,0.563366,0.560612,0.572312,0.561475,0.538044,0.548981,0.529354,0.528752,0.515555,0.498835,0.495908,0.489886,0.49115,0.485622,0.514526,0.476689};
        Double_t ratioBVar2[250] = {15.329,19.5563,20.2256,20.8922,20.7786,21.7281,22.1499,22.9172,23.5043,24.1949,24.8447,25.5106,26.0451,26.4445,26.6083,26.7088,26.6786,26.4278,26.0477,25.7804,25.1218,24.467,23.8927,23.1832,22.6776,22.0502,21.435,20.5799,19.7847,18.988,18.1516,17.4451,16.543,15.8462,14.8696,14.0847,13.3555,12.5591,11.7459,11.152,10.5494,9.71551,9.11848,8.67954,8.14362,7.67653,7.27814,6.81564,6.49001,6.1664,5.87361,5.62934,5.42766,5.16106,4.92883,4.79482,4.65459,4.36889,4.28666,4.09357,4.05095,3.87602,3.81197,3.68297,3.55739,3.44747,3.36267,3.25077,3.15145,3.10279,3.01222,2.96077,2.85096,2.7706,2.70887,2.63849,2.61934,2.52627,2.43338,2.36299,2.28566,2.22999,2.19306,2.13127,2.05379,2.04614,1.99544,1.94134,1.87212,1.87165,1.79433,1.76592,1.74631,1.70161,1.67131,1.61198,1.58921,1.53767,1.51985,1.48477,1.48072,1.4228,1.3919,1.37015,1.35963,1.33225,1.31461,1.2804,1.28466,1.26979,1.27146,1.22588,1.20816,1.20633,1.18425,1.15822,1.15496,1.15736,1.15494,1.13349,1.11876,1.09323,1.06654,1.06124,1.03183,0.998874,1.01905,1.00669,0.998078,0.999059,0.959233,0.963695,0.956146,0.948324,0.924652,0.911968,0.902708,0.904963,0.880048,0.882103,0.874787,0.872191,0.853132,0.853909,0.864728,0.839351,0.835888,0.808904,0.810313,0.804977,0.784115,0.765379,0.764144,0.754685,0.738616,0.747284,0.741864,0.730156,0.740405,0.732492,0.717926,0.714379,0.679377,0.696217,0.653297,0.656768,0.670578,0.65166,0.645874,0.641428,0.638373,0.631374,0.6011,0.599902,0.599082,0.584029,0.573703,0.570798,0.562921,0.54752,0.548524,0.537172,0.532661,0.555925,0.528359,0.528877,0.493402,0.504522,0.490647,0.498308,0.470658,0.482177,0.476074,0.458911,0.460249,0.442133,0.446161,0.433482,0.437694,0.418669,0.420143,0.427299,0.420959,0.400098,0.403357,0.389728,0.402054,0.390361,0.392965,0.38918,0.386478,0.364819,0.369727,0.354562,0.3535,0.337686,0.35478,0.32857,0.326081,0.348233,0.318214,0.319183,0.315676,0.314937,0.304259,0.315912,0.311705,0.323494,0.310281,0.298628,0.301669,0.283005,0.289184,0.300409,0.287236,0.286411,0.292979,0.288011,0.276546,0.282733,0.273169,0.273402,0.267107,0.258957,0.257946,0.255315,0.256477,0.254087,0.269736,0.250387};
        for (int idata=1; idata<251; idata++) {
            fBWeight->SetBinContent(idata,ratioB[idata-1]);
            fBWeightNew->SetBinContent(idata,ratioBNew[idata-1]);
            fBWeightVar1->SetBinContent(idata,ratioBVar1[idata-1]);
            fBWeightVar2->SetBinContent(idata,ratioBVar2[idata-1]);
        }
    }else{
        Double_t ratioB[250];
        Double_t ratioBNew[250];
        Double_t ratioBVar1[250];
        Double_t ratioBVar2[250];
        for (int idata=1; idata<251; idata++) {
            fBWeight->SetBinContent(idata,ratioB[idata-1]);
            fBWeightNew->SetBinContent(idata,ratioBNew[idata-1]);
            fBWeightVar1->SetBinContent(idata,ratioBVar1[idata-1]);
            fBWeightVar2->SetBinContent(idata,ratioBVar2[idata-1]);
        }
    }

    fOutputList->Add(fBWeight);
    fOutputList->Add(fBWeightNew);
    fOutputList->Add(fBWeightVar1);
    fOutputList->Add(fBWeightVar2);
    
    if (fFlagFillMCHistos) {
        fPi0DCA = new TH2F("fPi0DCA","Pi0 DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fPi0DCA);
    
        fEtaDCA = new TH2F("fEtaDCA","Eta DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fEtaDCA);
    
        fEnhEtaDCA = new TH2F("fEnhEtaDCA","Enh Eta DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fEnhEtaDCA->Sumw2();
        fOutputList->Add(fEnhEtaDCA);
        fEnhEtaWeightedPt = new TH1F("fEnhEtaWeightedPt","Enh Eta Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
        fEnhEtaWeightedPt->Sumw2();
        fOutputList->Add(fEnhEtaWeightedPt);
        fEnhPi0DCA = new TH2F("fEnhPi0DCA","Enh Pi0 DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fEnhPi0DCA->Sumw2();
        fOutputList->Add(fEnhPi0DCA);
        fEnhPi0WeightedPt = new TH1F("fEnhPi0WeightedPt","Enh Pi0 Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
        fEnhPi0WeightedPt->Sumw2();
        fOutputList->Add(fEnhPi0WeightedPt);
        fEtaHijingDCA = new TH2F("fEtaHijingDCA","Hijing Eta DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fEtaHijingDCA);
        fEtaHijingPt = new TH1F("fEtaHijingPt","Hijing Eta Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
        fOutputList->Add(fEtaHijingPt);
        fPi0HijingDCA = new TH2F("fPi0HijingDCA","Hijing Pi0 DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fPi0HijingDCA);
        fPi0HijingPt = new TH1F("fPi0HijingPt","Hijing Pi0 Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
        fOutputList->Add(fPi0HijingPt);
        fPhotonHijingDCA = new TH2F("fPhotonHijingDCA","Hijing Photon DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fPhotonHijingDCA);
        fPhotonHijingPt = new TH1F("fPhotonHijingPt","Hijing Photon Weighted pT; p_{T}(GeV/c); counts;", 60,0,30.);
        fOutputList->Add(fPhotonHijingPt);
        fEnhPhotonDCA = new TH2F("fEnhPhotonDCA","Enh Photon DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
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
    
        fBMesonPt = new TH1F("fBMesonPt","B Meson Spectrum; p_{T}(GeV/c); counts;",250,0,50.);
        fBMesonPt->Sumw2();
        fOutputList->Add(fBMesonPt);
    
        fBMesonPtATLAS = new TH1F("fBMesonPtATLAS","ATLAS B Meson Spectrum; p_{T}(GeV/c); counts;",240,0,120.);
        fBMesonPtATLAS->Sumw2();
        fOutputList->Add(fBMesonPtATLAS);
        
        fBPlusPtATLAS = new TH1F("fBPlusPtATLAS","ATLAS B Plus Spectrum; p_{T}(GeV/c); counts;",240,0,120.);
        fBPlusPtATLAS->Sumw2();
        fOutputList->Add(fBPlusPtATLAS);
        
        fBMesonPtCMS = new TH1F("fBMesonPtCMS","CMS B Meson Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
        fBMesonPtCMS->Sumw2();
        fOutputList->Add(fBMesonPtCMS);
        
        fBPlusPtCMS = new TH1F("fBPlusPtCMS","CMS B Plus Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
        fBPlusPtCMS->Sumw2();
        fOutputList->Add(fBPlusPtCMS);
        
        fBMesonPtLHCb = new TH1F("fBMesonPtLHCb","LHCb B Meson Spectrum; p_{T}(GeV/c); counts;",400,0,40.);
        fBMesonPtLHCb->Sumw2();
        fOutputList->Add(fBMesonPtLHCb);
        
        fBPlusPtLHCb = new TH1F("fBPlusPtLHCb","LHCb B Plus Spectrum; p_{T}(GeV/c); counts;",400,0,40.);
        fBPlusPtLHCb->Sumw2();
        fOutputList->Add(fBPlusPtLHCb);
        
        fBBaryonPt = new TH1F("fBBaryonPt","Beauty Baryon Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
        fBBaryonPt->Sumw2();
        fOutputList->Add(fBBaryonPt);
    
        fBMesonElecPt = new TH1F("fBMesonElecPt","Beauty Meson->e Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
        fBMesonElecPt->Sumw2();
        fOutputList->Add(fBMesonElecPt);
        
        fBBaryonElecPt = new TH1F("fBBaryonElecPt","Beauty Baryon->e Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
        fBBaryonElecPt->Sumw2();
        fOutputList->Add(fBBaryonElecPt);
        
        fPromptD0DCAWeight = new TH2F("fPromptD0DCAWeight","Prompt D0 DCA with Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fPromptD0DCAWeight);
    
        fD0FromDStarDCAWeight = new TH2F("fD0FromDStarDCAWeight","Prompt D0 DCA with Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fD0FromDStarDCAWeight);
    
        fPromptD0DCANoWeight = new TH2F("fPromptD0DCANoWeight","Prompt D0 DCA without Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fPromptD0DCANoWeight);
    
        fD0FromDStarDCANoWeight = new TH2F("fD0FromDStarDCANoWeight","Prompt D0 DCA without Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
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
    
        Int_t binTemp[5] = {60,nDCAbins,19,3,50}; //pT, DCA, Mom PID, Mom Gen, mompT
        Double_t xminTemp[5] = {0.,-0.2,0.5,-0.5,0.};
        Double_t xmaxTemp[5] = {30.,0.2,19.5,2.5,50.};
        fSprsTemplatesNoWeight = new THnSparseD("fSprsTemplatesNoWeight","Sparse for Templates, No weight applied;p_{T};DCA;MomPID;MomGen;",5,binTemp,xminTemp,xmaxTemp);
        fOutputList->Add(fSprsTemplatesNoWeight);
        fSprsTemplatesNoWeight->Sumw2();
        fSprsTemplatesWeight = new THnSparseD("fSprsTemplatesWeight","Sparse for Templates, D meson weight applied;p_{T};DCA;MomPID;MomGen;",5,binTemp,xminTemp,xmaxTemp);
        fOutputList->Add(fSprsTemplatesWeight);
        fSprsTemplatesWeight->Sumw2();
        fSprsTemplatesWeightVar1 = new THnSparseD("fSprsTemplatesWeightVar1","Sparse for Templates, D meson WeightVar1 applied;p_{T};DCA;MomPID;MomGen;",5,binTemp,xminTemp,xmaxTemp);
        fOutputList->Add(fSprsTemplatesWeightVar1);
        fSprsTemplatesWeightVar1->Sumw2();
        fSprsTemplatesWeightVar2 = new THnSparseD("fSprsTemplatesWeightVar2","Sparse for Templates, D meson WeightVar2 applied;p_{T};DCA;MomPID;MomGen;",5,binTemp,xminTemp,xmaxTemp);
        fOutputList->Add(fSprsTemplatesWeightVar2);
        fSprsTemplatesWeightVar2->Sumw2();
    
        fDTemplateWeight = new TH2F("fDTemplateWeight","D Meson DCA template", 100,0,50., nDCAbins,-0.2,0.2);
        fOutputList->Add(fDTemplateWeight);
    
        fDTemplateNoWeight = new TH2F("fDTemplateNoWeight","D Meson DCA template w/o Weight", 100,0,50., nDCAbins,-0.2,0.2);
        fOutputList->Add(fDTemplateNoWeight);
    
        fDTemplateWeightNew = new TH2F("fDTemplateWeightNew","New D Meson DCA template", 100,0,50., nDCAbins,-0.2,0.2);
        fOutputList->Add(fDTemplateWeightNew);
        fDTemplateWeightVar1 = new TH2F("fDTemplateWeightVar1","D Meson DCA template Var1", 100,0,50., nDCAbins,-0.2,0.2);
        fOutputList->Add(fDTemplateWeightVar1);
        fDTemplateWeightVar2 = new TH2F("fDTemplateWeightVar2","D Meson DCA template Var2", 100,0,50., nDCAbins,-0.2,0.2);
        fOutputList->Add(fDTemplateWeightVar2);
        
        fBTemplateWeight = new TH2F("fBTemplateWeight","B Meson DCA template w/Weight", 100,0,50., nDCAbins,-0.2,0.2);
        fOutputList->Add(fBTemplateWeight);
    
        fBTemplateNoWeight = new TH2F("fBTemplateNoWeight","B Meson DCA template", 100,0,50., nDCAbins,-0.2,0.2);
        fOutputList->Add(fBTemplateNoWeight);
    
        fBTemplateWeightNew = new TH2F("fBTemplateWeightNew","B Meson DCA template w/Weight New", 100,0,50., nDCAbins,-0.2,0.2);
        fOutputList->Add(fBTemplateWeightNew);
        fBTemplateWeightVar1 = new TH2F("fBTemplateWeightVar1","B Meson DCA template w/Weight Var1", 100,0,50., nDCAbins,-0.2,0.2);
        fOutputList->Add(fBTemplateWeightVar1);
        fBTemplateWeightVar2 = new TH2F("fBTemplateWeightVar2","B Meson DCA template w/Weight Var2", 100,0,50., nDCAbins,-0.2,0.2);
        fOutputList->Add(fBTemplateWeightVar2);
        
        fAllElecStack = new TH1F("fAllElecStack","All Elec from Stack; p_{T}(GeV/c); counts;",100,0,50.);
        fAllElecStack->Sumw2();
        fOutputList->Add(fAllElecStack);
    
        fHFElecStack = new TH1F("fHFElecStack","HF Elec from Stack; p_{T}(GeV/c); counts;",100,0,50.);
        fHFElecStack->Sumw2();
        fOutputList->Add(fHFElecStack);
    
        fBElecStack = new TH1F("fBElecStack","B Elec from Stack; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecStack->Sumw2();
        fOutputList->Add(fBElecStack);
    
        fElecTPCTrk = new TH1F("fElecTPCTrk","Elec TPC tracks; p_{T}(GeV/c); counts;",100,0,50.);
        fElecTPCTrk->Sumw2();
        fOutputList->Add(fElecTPCTrk);
    
        fHFElecTPCTrk = new TH1F("fHFElecTPCTrk","HF Elec TPC tracks; p_{T}(GeV/c); counts;",100,0,50.);
        fHFElecTPCTrk->Sumw2();
        fOutputList->Add(fHFElecTPCTrk);
    
        fBElecTPCTrk = new TH1F("fBElecTPCTrk","B Elec TPC tracks; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecTPCTrk->Sumw2();
        fOutputList->Add(fBElecTPCTrk);
    
        fElecAftTrkCuts = new TH1F("fElecAftTrkCuts","Elec after trk cuts; p_{T}(GeV/c); counts;",100,0,50.);
        fElecAftTrkCuts->Sumw2();
        fOutputList->Add(fElecAftTrkCuts);
    
        fHFElecAftTrkCuts = new TH1F("fHFElecAftTrkCuts","HF Elec after trk cuts; p_{T}(GeV/c); counts;",100,0,50.);
        fHFElecAftTrkCuts->Sumw2();
        fOutputList->Add(fHFElecAftTrkCuts);
    
        fBElecAftTrkCuts = new TH1F("fBElecAftTrkCuts","B Elec after trk cuts; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecAftTrkCuts->Sumw2();
        fOutputList->Add(fBElecAftTrkCuts);
    
        fElecAftTrkMatch = new TH1F("fElecAftTrkMatch","Elec after trk Match; p_{T}(GeV/c); counts;",100,0,50.);
        fElecAftTrkMatch->Sumw2();
        fOutputList->Add(fElecAftTrkMatch);
    
        fHFElecAftTrkMatch = new TH1F("fHFElecAftTrkMatch","HF Elec after trk Match; p_{T}(GeV/c); counts;",100,0,50.);
        fHFElecAftTrkMatch->Sumw2();
        fOutputList->Add(fHFElecAftTrkMatch);
    
        fBElecAftTrkMatch = new TH1F("fBElecAftTrkMatch","B Elec after trk Match; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecAftTrkMatch->Sumw2();
        fOutputList->Add(fBElecAftTrkMatch);
    
        fElecAftTPCeID = new TH1F("fElecAftTPCeID","Elec after TPC eID; p_{T}(GeV/c); counts;",100,0,50.);
        fElecAftTPCeID->Sumw2();
        fOutputList->Add(fElecAftTPCeID);
    
        fHFElecAftTPCeID = new TH1F("fHFElecAftTPCeID","HF Elec after TPC eID; p_{T}(GeV/c); counts;",100,0,50.);
        fHFElecAftTPCeID->Sumw2();
        fOutputList->Add(fHFElecAftTPCeID);
    
        fBElecAftTPCeID = new TH1F("fBElecAftTPCeID","B Elec after TPC eID; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecAftTPCeID->Sumw2();
        fOutputList->Add(fBElecAftTPCeID);
    
        fElecAftEMCeID = new TH1F("fElecAftEMCeID","Elec after EMC eID; p_{T}(GeV/c); counts;",100,0,50.);
        fElecAftEMCeID->Sumw2();
        fOutputList->Add(fElecAftEMCeID);
    
        fHFElecAftEMCeID = new TH1F("fHFElecAftEMCeID","HF Elec after EMC eID; p_{T}(GeV/c); counts;",100,0,50.);
        fHFElecAftEMCeID->Sumw2();
        fOutputList->Add(fHFElecAftEMCeID);
    
        fBElecAftEMCeID = new TH1F("fBElecAftEMCeID","B Elec after EMC eID; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecAftEMCeID->Sumw2();
        fOutputList->Add(fBElecAftEMCeID);
    }
    
    if (fFlagFillSprs && fFlagFillMCHistos) {
        Int_t bins1[6]=  {/*280*/60,  160, 100, 100,  nDCAbins, 4}; // pT;nSigma;eop;M20;DCA;MCTruth
        Double_t xmin1[6]={ /*2*/0,   -8,   0,   0, -0.2, -0.5};
        Double_t xmax1[6]={30,    8,   2,   1,  0.2, 3.5};
        fElectronSprs = new THnSparseD("Electron","Electron;pT;nSigma;eop;DCA;MCTruth;",6,bins1,xmin1,xmax1);
        fOutputList->Add(fElectronSprs);
    }
    if (fFlagFillSprs && !fFlagFillMCHistos) {
        Int_t bins1[5]=  {60,  160, 100, 100,  nDCAbins}; // pT;nSigma;eop;M20;DCA;MCTruth
        Double_t xmin1[5]={ /*2*/0,   -8,   0,   0, -0.2};
        Double_t xmax1[5]={30,    8,   2,   1,  0.2};
        fElectronSprs = new THnSparseD("Electron","Electron;pT;nSigma;eop;DCA;MCTruth;",5,bins1,xmin1,xmax1);
        fOutputList->Add(fElectronSprs);
    }
    
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
    if(fApplyCentrality){
        Bool_t pass = kFALSE;
        if(fCentralityMin > -0.5){
            fCentrality = CheckCentrality(fAOD,pass);
            if(!pass)return;
        }
        fCentCheck->Fill(fCentrality);
    }
    
    ////////////////////
    // Get MC Headers //
    ////////////////////
    fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    
    //Get NParticles from the generators
    if (fFlagFillMCHistos) {
        if (fMCarray && fMCHeader) {
            //cout<<"Test111..................................."<<endl;
            GetNMCPartProduced();
            // cout<<"Total Number of Particles = "<<fNtotMCpart<< "," << fMCarray->GetEntries() <<endl;
        }
    }
    
    ///////////////////
    //Loop over Stack//
    ///////////////////
    /*if (fFlagFillMCHistos && !fFlagRunStackLoop) {
        Int_t TrackPDG = -999;
        Int_t ilabelM = -99;
        for(int i=0; i<(fMCarray->GetEntries()); i++)
        {
            AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(i);
            if(TMath::Abs(AODMCtrack->Eta()) > 0.6) continue;
            
            //-------Get PDG
            TrackPDG = TMath::Abs(AODMCtrack->GetPdgCode());
            ilabelM = AODMCtrack->GetMother();
            
            //Electrons only
            if(TrackPDG != 11 || !AODMCtrack->IsPhysicalPrimary()) continue;
            fAllElecStack->Fill(AODMCtrack->Pt());
            
            //Fill Charm species pT histos
            if (ilabelM>0) {
                //cout<<"Test4..................................."<<endl;
                AliAODMCParticle *momPart = (AliAODMCParticle*)fMCarray->At(ilabelM); //get mom particle
                Int_t pidM = TMath::Abs(momPart->GetPdgCode());
                
                if(pidM>400 && pidM<600 && AODMCtrack->IsPhysicalPrimary()) {
                    fHFElecStack->Fill(AODMCtrack->Pt());
                    
                }
                if(pidM>500 && pidM<600 && AODMCtrack->IsPhysicalPrimary()) {
                    fBElecStack->Fill(AODMCtrack->Pt());
                }
            }
        }
    }*/
    
    
    if (fFlagFillMCHistos && fFlagRunStackLoop) {
        Int_t eleinStack=0;
        if (fMCarray) {
            //cout<<"Test2..................................."<<endl;
            // Make Pi0 and Eta Weight Sparse
            GetPi0EtaWeight(fSprsPi0EtaWeightCal);
            Int_t TrackPDG = -999;
            Int_t ilabelM = -99;
            Int_t ilabelGM = -99;
            Int_t ilabelGGM = -99;
            Int_t pidM, pidGM, pidGGM;
            
            for(int i=0; i<(fMCarray->GetEntries()); i++)
            {
                //cout<<"Test "<<fMCarray->GetEntries()<<" ..................................."<<i<<endl;
                //if (i<fNpureMC) continue; //reject plain Hijing
            
                Bool_t fromDStar = kFALSE;
            
                AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(i);
            
                //-------Get PDG
                TrackPDG = TMath::Abs(AODMCtrack->GetPdgCode());
                ilabelM = AODMCtrack->GetMother();
            
                if(TrackPDG == 11 && AODMCtrack->IsPhysicalPrimary()) {
                    // cout<<"TESTINGGGGGGGGGGGG"<<endl;
                    fAllElecStack->Fill(AODMCtrack->Pt());
                    eleinStack++;
                }
                //
                if(TMath::Abs(AODMCtrack->Y()) < 2.25) {
                    if (TrackPDG>500 && TrackPDG<599) fBMesonPtATLAS->Fill(AODMCtrack->Pt());
                    if (TrackPDG == 521) fBPlusPtATLAS->Fill(AODMCtrack->Pt());
                }
                if(TMath::Abs(AODMCtrack->Y()) < 2.4) {
                    if (TrackPDG>500 && TrackPDG<599) fBMesonPtCMS->Fill(AODMCtrack->Pt());
                    if (TrackPDG == 521) fBPlusPtCMS->Fill(AODMCtrack->Pt());
                    /*if(TrackPDG==11 && ilabelM>0){
                        //cout<<"TESTING2"<<endl;
                        AliAODMCParticle *momPart = (AliAODMCParticle*)fMCarray->At(ilabelM); //get mom particle
                        pidM = TMath::Abs(momPart->GetPdgCode());
                        
                        if(pidM>500 && pidM<599) fBMesonElecPt->Fill(AODMCtrack->Pt());
                        if(pidM>5000 && pidM<5999) fBBaryonElecPt->Fill(AODMCtrack->Pt());
                        
                        ilabelGM = momPart->GetMother();
                        if(ilabelGM>0){
                            pidGM = TMath::Abs(momPart->GetPdgCode());
                            if(pidGM>500 && pidGM<599) fBMesonElecPt->Fill(AODMCtrack->Pt());
                            AliAODMCParticle *gMomPart = (AliAODMCParticle*)fMCarray->At(ilabelGM); //get mom particle
                            ilabelGGM = gMomPart->GetMother();
                            if(ilabelGGM>0){
                                pidGGM = TMath::Abs(gMomPart->GetPdgCode());
                                if(pidGGM>500 && pidGGM<599) fBMesonElecPt->Fill(AODMCtrack->Pt());
                            }
                        }
                    }*/
                }
                if(TMath::Abs(AODMCtrack->Y()) > 2.0 && TMath::Abs(AODMCtrack->Y()) < 4.5) {
                    if (TrackPDG>500 && TrackPDG<599) fBMesonPtLHCb->Fill(AODMCtrack->Pt());
                    if (TrackPDG == 521) fBPlusPtLHCb->Fill(AODMCtrack->Pt());
                }
                
                if(TMath::Abs(AODMCtrack->Eta()) > 0.6) continue;
                if (TrackPDG>500 && TrackPDG<599) {
                    fBMesonPt->Fill(AODMCtrack->Pt());
                }
                if (TrackPDG>5000 && TrackPDG<5999) {
                    fBBaryonPt->Fill(AODMCtrack->Pt());
                }
            
                //Fill Charm species pT histos
                if (ilabelM>0) {
                    //cout<<"Test4..................................."<<endl;
                
                    if(TrackPDG == 11 && pidM>400 && pidM<600 && AODMCtrack->IsPhysicalPrimary()) {
                        fHFElecStack->Fill(AODMCtrack->Pt());
                    
                    }
                    if(TrackPDG == 11 && pidM>500 && pidM<600 && AODMCtrack->IsPhysicalPrimary()) {
                        fBElecStack->Fill(AODMCtrack->Pt());
                        fBMesonElecPt->Fill(AODMCtrack->Pt());
                    }
                    if(TrackPDG == 11 && pidM>5000 && pidM<5999 && AODMCtrack->IsPhysicalPrimary()) {
                        fBBaryonElecPt->Fill(AODMCtrack->Pt());
                    }
                
                    if (pidM==413) fromDStar = kTRUE;
                    if (pidM>500 && pidM<599) {
                        continue; //reject beauty feed down
                    }
                    if (pidM>5000 && pidM<5999) {
                        continue; //reject beauty feed down
                    }
                    AliAODMCParticle *momPart = (AliAODMCParticle*)fMCarray->At(ilabelM); //get mom particle
                    ilabelGM = momPart->GetMother();
                    if (ilabelGM>0) {
                        AliAODMCParticle *gmomPart = (AliAODMCParticle*)fMCarray->At(ilabelGM);//get grandma particle
                        pidGM = TMath::Abs(gmomPart->GetPdgCode());
                        if (pidGM>500 && pidGM<599) {
                            if(TrackPDG == 11 && AODMCtrack->IsPhysicalPrimary()) fBMesonElecPt->Fill(AODMCtrack->Pt());
                            continue; //reject beauty feed down
                        }
                        if (pidGM>5000 && pidGM<5999) {
                            continue; //reject beauty feed down
                        }
                        ilabelGGM = gmomPart->GetMother();
                        if (ilabelGGM>0) {
                            AliAODMCParticle *ggmomPart = (AliAODMCParticle*)fMCarray->At(ilabelGGM); //get great grandma particle
                            pidGGM = TMath::Abs(ggmomPart->GetPdgCode());
                            if (pidGGM>500 && pidGGM<599) {
                                if(TrackPDG == 11 && AODMCtrack->IsPhysicalPrimary()) fBMesonElecPt->Fill(AODMCtrack->Pt());
                            }
                        }
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
                // cout<<"Total Number of Particles = "<<fNtotMCpart<<endl;
            }
        }
    }
    //cout << "Electron in stack -------------- " << eleinStack <<endl;
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
    Double_t clsphi = -999, clseta=-999;
    nclus = fAOD->GetNumberOfCaloClusters();
    /*for (Int_t icl = 0; icl < nclus; icl++) {
        //ESD and AOD CaloCells carries the same information
        AliVCluster* clus = (AliAODCaloCluster*)fAOD->GetCaloCluster(icl);
        if(clus && clus->IsEMCAL()){
            fClsEAll->Fill(clus->E()); //E of all clusters
        }
    }*/
    for (Int_t icl = 0; icl < nclus; icl++) {
        //ESD and AOD CaloCells carries the same information
        AliVCluster* clus = (AliAODCaloCluster*)fAOD->GetCaloCluster(icl);
        if(clus && clus->IsEMCAL()){
            Float_t  emcPos[3]; // cluster pos
            clus->GetPosition(emcPos);
            TVector3 clustVec(emcPos[0],emcPos[1],emcPos[2]);
            clsphi = clustVec.Phi();
            clseta = clustVec.Eta();
            if(clsphi < 0) clsphi = clsphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
            if(fFlagClsTypeEMC && !fFlagClsTypeDCAL){ //if we want EMCal
                if(clsphi > 4.53 && clsphi < 5.708) { //DCAL  : 260 < phi < 327 but it's DCal
                    continue;//leave out DCal
                }
            }
            if(fFlagClsTypeDCAL && !fFlagClsTypeEMC){//if we want DCal
                if(clsphi > 1.39 && clsphi < 3.265) {//EMCAL : 80 < phi < 187 but it's EMCal
                    continue; //leave out EMCal
                }
            }
            
            fClsEAll->Fill(clus->E()); //E of all clusters
        }
    }
    
    Int_t eleinTrkLoop=0;
    ////////////////
    // track loop //
    ////////////////
    Int_t nTracks(fAOD->GetNumberOfTracks());
    for(Int_t i=0; i<nTracks; i++){
        
        AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track) continue;
        
        if(TMath::Abs(track->Eta()) > 0.6) continue;
        
        fTrkPtB4TC->Fill(track->Pt());
        
        //See if true electron
        Bool_t kTruElec = kFALSE;
        Bool_t kTruHFElec = kFALSE;
        Bool_t kTruBElec = kFALSE;
        Int_t pdg = -99;
        Int_t pidM = -99;
        Int_t ilabelM = -99;
        Int_t ilabel = -99;
        if (fFlagFillMCHistos) {
            if(fMCarray)
            {
                ilabel = TMath::Abs(track->GetLabel()); //get MC label of track
                if(ilabel == 0) continue;
                // ilabel = track->GetLabel();
                // if(ilabel < -1) continue;
                //cout <<"ilabel = "<<ilabel<<"************************"<<endl;
                fMCparticle = (AliAODMCParticle*) fMCarray->At(ilabel);
                pdg = TMath::Abs(fMCparticle->GetPdgCode()); //get pid of track
            
                if(TMath::Abs(fMCparticle->Eta()) > 0.6) continue;
                //cout<<"TESTING1234"<<endl;
            
                //if electron--------------------------------
                if(pdg==11 && fMCparticle->IsPhysicalPrimary()){
                    eleinTrkLoop++;
                    kTruElec = kTRUE;
                    //  cout<<"TESTING12345674"<<endl;
                    ilabelM = fMCparticle->GetMother();
                    if(ilabelM>0){
                        AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(ilabelM); //get mom particle
                        pidM = TMath::Abs(partM->GetPdgCode()); //ask for the Mom's pid
                        //    cout << "Test pidM = "<<pidM<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                        if(pidM>400 && pidM<600) kTruHFElec = kTRUE;
                        if(pidM>500 && pidM<600) kTruBElec = kTRUE;
                    }
                }
            }
            if(kTruElec == kTRUE) fElecTPCTrk->Fill(track->Pt());
            if(kTruHFElec == kTRUE) fHFElecTPCTrk->Fill(track->Pt());
            if(kTruBElec == kTRUE) fBElecTPCTrk->Fill(track->Pt());
        }
        
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
        if(track->GetTPCNcls() < fNclusTPC) continue;
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
        
        if(kTruElec == kTRUE) fElecAftTrkCuts->Fill(track->Pt());
        if(kTruHFElec == kTRUE) fHFElecAftTrkCuts->Fill(track->Pt());
        if(kTruBElec == kTRUE) fBElecAftTrkCuts->Fill(track->Pt());
        
        Double_t nsigma = -999;
        nsigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        
        fnSigma->Fill(track->Pt(),nsigma);
        if(nsigma>-5.&&nsigma<-3.) {
            fHadronCamDCA->Fill(track->Pt(),d0z0[0]);
        }
        if (fFlagFillMCHistos) {
            if(fMCarray)
            {
                if(TMath::Abs(fMCparticle->Eta()) > 0.6) continue;
                //Fill Hijing Hadron DCA
                if(ilabel<fNpureMC && nsigma>-5. && nsigma<-3.) {
                fHadronCamDCAHij->Fill(track->Pt(),d0z0[0]);
                }
            }
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
            
            if(TMath::Abs(fPhiDiff) > fTrkMatch || TMath::Abs(fEtaDiff)> fTrkMatch) continue;
            
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
            
            if(kTruElec == kTRUE) fElecAftTrkMatch->Fill(track->Pt());
            if(kTruHFElec == kTRUE) fHFElecAftTrkMatch->Fill(track->Pt());
            if(kTruBElec == kTRUE) fBElecAftTrkMatch->Fill(track->Pt());
            
            fnSigmaAftTrkMatch->Fill(track->Pt(),nsigma);
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
            Double_t bWeight = -99;
            Bool_t kEmbEta = kFALSE;
            Bool_t kEmbPi0 = kFALSE;
            Bool_t kHijing = kFALSE;
            Bool_t kFlagReco = kFALSE;
            Int_t fMomGen = 99;
            Double_t momPt = -99;
            Int_t pidGM = -99;
            //Int_t ilabel = -99;
            Int_t ilabelM = -99;
            Int_t ilabelGM = -99;
            
            //cout<<"TESTING0"<<endl;
            //if MC--------------------------------
            if (fFlagFillMCHistos) {
                if(ilabel>0 && fMCarray)
                {
                    //cout<<"TESTING1"<<endl;
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(ilabel);
                    pdg = fMCparticle->GetPdgCode(); //get pid of track
                
                    //cout<<"TESTING1"<<endl;
                
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
                        if (fpidSort==2||fpidSort==11||fpidSort==12||fpidSort==14||fpidSort==15||fpidSort==16) { //if from D meson
                            //cout<<"TESTING3"<<endl;
                            if (momPt>1 && momPt<50.) { //in proper pt range
                                //cout<<"TESTING4"<<endl;
                                dWeight = fDWeight->GetBinContent(fDWeight->FindBin(momPt));
                                fDTemplateWeight->Fill(track->Pt(), DCA, dWeight);
                                fDTemplateNoWeight->Fill(track->Pt(), DCA);
                                
                                dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(momPt));
                                fSprsTemplatesWeight->Fill(tempValue,dWeight);
                                fDTemplateWeightNew->Fill(track->Pt(), DCA, dWeight);
                                
                                dWeight = fDWeightVar1->GetBinContent(fDWeightVar1->FindBin(momPt));
                                fDTemplateWeightVar1->Fill(track->Pt(), DCA, dWeight);
                                fSprsTemplatesWeightVar1->Fill(tempValue,dWeight);
                                
                                dWeight = fDWeightVar2->GetBinContent(fDWeightVar2->FindBin(momPt));
                                fDTemplateWeightVar2->Fill(track->Pt(), DCA, dWeight);
                                fSprsTemplatesWeightVar2->Fill(tempValue,dWeight);
                            }
                        }
                        else if (fpidSort==1) {//if from B meson
                            //cout<<"TESTING5"<<endl;
                            if (momPt>0. && momPt<50.) { //in proper pt range
                                //cout<<"TESTING6"<<endl;
                                bWeight = fBWeight->GetBinContent(fBWeight->FindBin(momPt));
                                fBTemplateWeight->Fill(track->Pt(), DCA, bWeight);
                                fBTemplateNoWeight->Fill(track->Pt(), DCA);
                                
                                bWeight = fBWeightNew->GetBinContent(fBWeightNew->FindBin(momPt));
                                fBTemplateWeightNew->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeight->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar1->GetBinContent(fBWeightVar1->FindBin(momPt));
                                fBTemplateWeightVar1->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeightVar1->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar2->GetBinContent(fBWeightVar2->FindBin(momPt));
                                fBTemplateWeightVar2->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeightVar2->Fill(tempValue,bWeight);
                            }
                        }
                        else{
                            //cout<<"TESTING7"<<endl;
                            fSprsTemplatesWeight->Fill(tempValue);
                            fSprsTemplatesWeightVar1->Fill(tempValue);
                            fSprsTemplatesWeightVar2->Fill(tempValue);
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
            }
            //////////////////////
            // Get MC True DCAs //
            //////////////////////
            //if MC--------------------------------
            if (fFlagFillMCHistos) {
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
            }
            /////////////////////
            // Electron sparse //
            /////////////////////
            Double_t EovP = (clustMatch->E())/(track->P());
            Double_t M20 = clustMatch->GetM20();
            Double_t M02 = clustMatch->GetM02();
            
            if(fFlagFillMCHistos && fFlagFillSprs) {
                Double_t fvalueElectron[6] = {-999,-999,-999,-999,-999,-999};
                fvalueElectron[0] = track->Pt();
                fvalueElectron[1] = nsigma;
                fvalueElectron[2] = EovP;
                fvalueElectron[3] = M20;
                fvalueElectron[4] = DCA;
                fvalueElectron[5] = 0;
                if(kTruElec) fvalueElectron[5] = 1;
                if(kTruHFElec) fvalueElectron[5] = 2;
                if(kTruBElec) fvalueElectron[5] = 3;
                fElectronSprs->Fill(fvalueElectron);
            }
            if(!fFlagFillMCHistos && fFlagFillSprs){
                Double_t fvalueElectron[5] = {-999,-999,-999,-999,-999};
                fvalueElectron[0] = track->Pt();
                fvalueElectron[1] = nsigma;
                fvalueElectron[2] = EovP;
                fvalueElectron[3] = M20;
                fvalueElectron[4] = DCA;
                fElectronSprs->Fill(fvalueElectron);
            }
            
            ///////////////////////////
            // Hadron Contam. Histos //
            ///////////////////////////
            if(nsigma<-4.) {
                if(M20>0.01 && M20<fMaxM20Cut) {
                    fHadronEoP->Fill(track->Pt(),EovP);
                    if(EovP>fMinEoPCut && EovP<1.2){
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
            
            if((EovP>0.9) && (EovP<1.2)) {
                fnSigaftEoPCut->Fill(track->Pt(),nsigma);
                if (M20>0.01 && M20<0.35) {
                    fnSigaftM20EoPCut->Fill(track->Pt(),nsigma);
                }
            }
            
            if((EovP>fMinEoPCut) && (EovP<1.2)) {
                fTPCElecEoP->Fill(track->Pt(),EovP);
                fnSigaftSysEoPCut->Fill(track->Pt(),nsigma);
            }
            
            //Apply M20 cut for electrons
            if((M20<0.01) || (M20>fMaxM20Cut)) continue;
            fElecEoPnoSig->Fill(track->Pt(),EovP);
            
            if((nsigma>fMinNSigCut) && (nsigma<3)) fInclElecEoP->Fill(track->Pt(),EovP);
            
            //Apply E/p Cut for electrons
            if((EovP<fMinEoPCut) || (EovP>1.2)) continue;
            fnSigaftSysM20EoPCut->Fill(track->Pt(),nsigma);
            
            if(kTruElec == kTRUE) fElecAftEMCeID->Fill(track->Pt());
            if(kTruHFElec == kTRUE) fHFElecAftEMCeID->Fill(track->Pt());
            if(kTruBElec == kTRUE) fBElecAftEMCeID->Fill(track->Pt());
            
            //Apply TPC nSigma cut for electrons
            if((nsigma<fMinNSigCut) || (nsigma>3)) continue;
            
            if(kTruElec == kTRUE) fElecAftTPCeID->Fill(track->Pt());
            if(kTruHFElec == kTRUE) fHFElecAftTPCeID->Fill(track->Pt());
            if(kTruBElec == kTRUE) fBElecAftTPCeID->Fill(track->Pt());
            
            if(fClsTypeDCAL) fClsEamElecDC->Fill(clustMatch->E());
            if(fClsTypeEMC) fClsEamElecEMC->Fill(clustMatch->E());
            
            /////////////////////////
            // Plot Reco Electrons //
            /////////////////////////
            fInclElecDCA->Fill(track->Pt(),DCA);
            fInclElecDCAnoSign->Fill(track->Pt(),d0z0[0]);
            
            if(!fFlagFillMCHistos){
                InvMassCheckData(i, track, d0z0, MagSign);
            }
            //Make incl electron and photonic electron plots
            /*if(nsigma>fMinNSigCut && nsigma<3) {
             if(M20>0.01 && M20<fMaxM20Cut) {
             fInclElecEoP->Fill(track->Pt(),EovP);
             if(EovP>fMinEoPCut && EovP<1.2){
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
        //   cout<<"Testiessssssss2"<<endl;
        AliError("no MC header");
        return (0);
    }
    //loop through headers
    for(int igene=0; igene<lh->GetEntries(); igene++)
    {
        //   cout<<"Testiessssssss"<<endl;
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
    // cout<<"TestB..................................."<<endl;
    for(int imc=0; imc < fNtotMCpart; imc++)
    {
        AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(imc);
        if(TMath::Abs(AODMCtrack->Eta()) > 0.9) continue;
        
        //cout<<"TestA..................................."<<endl;
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
                    momPt = partGM->Pt();
                }
                if(pidGM>5000 && pidGM<5999){
                    fpidSort = 10; //GMa is b baryon
                }
            }
            //check if charm baryon grandma is B
            if(pidM>4000 && pidM<4999){
                if(pidGM>500 && pidGM<599){
                    fpidSort = 1; //GMa is B
                    momPt = partGM->Pt();
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
                        momPt = partGGM->Pt();
                    }
                    if(pidGGM>5000 && pidGGM<5999){
                        fpidSort = 10; //GGMa is b baryon
                    }
                }
                //check if charm baryon great grandma is B
                if(pidM>4000 && pidM<4999){
                    if(pidGGM>500 && pidGGM<599){
                        fpidSort = 1; //GGMa is B
                        momPt = partGGM->Pt();
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
        //if(ptAsso < 0.3) continue;
        if(ptAsso < fMinPtAssoCut) continue;
        if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
        //if(nsigmaAsso < -3 || nsigmaAsso > 3) continue;
        if(nsigmaAsso < fMinNSigAssoCut || nsigmaAsso > 3) continue;
        
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
        if(ptAsso < fMinPtAssoCut) continue;
        if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
        if(nsigmaAsso < fMinNSigAssoCut || nsigmaAsso > 3) continue;
        
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
