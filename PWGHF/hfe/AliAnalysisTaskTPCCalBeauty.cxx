//
//  AliAnalysisTaskTPCCalBeauty.cxx
//
//
//  Created by Erin Gauger
//
//

#include "TChain.h"
#include "TH1F.h"
#include "TH3F.h"
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
#include "AliESDtrack.h"

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
fApplyM02Cut(kFALSE),
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
fDCAzCut(3.2),
fMinMass(0.),
fMaxMass(0.1),
fMinEta(-0.6),
fMaxEta(0.6),
fAssoDCAxy(0.25),
fAssoDCAz(1.),
fAssoTPCnCls(80),
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
fTrkMatch(0),
fUseTender(kTRUE),
fApplyHadEoPCut(kTRUE),
fVtxZCut(0.),
fDCAxyCut(2.4),
fFlagULS(kFALSE),
fFlagLS(kFALSE),
fNevents(0),
fVtX(0),
fVtY(0),
fVtZ(0),
fTrkPtB4TC(0),
fDCAxyz(0),
fTrkPt(0),
fTrkP(0),
fTrkClsPhi(0),
fTrkClsEta(0),
fClsPhi(0),
fClsEta(0),
fClsE(0),
fClsEnoTimeCut(0),
//fClsEamDCal(0),
//fClsEamEMCal(0),
//fClsEAll(0),
//fClsEamElecEMC(0),
//fClsEamElecDC(0),
fTrkPhi(0),
fTrkEta(0),
fdEdx(0),
fnSigma(0),
fnSigmaAftTrkMatch(0),
fCentCheck(0),
fTrigCheck(0),
fITSLayerCheck(0),
fEMCTrkMatch(0),
fBasicSprs(0),
fvalueBasic(0),
//fElecDcaB4TrkMatch(0),
//fDcaB4TrkMatch(0),
//fElecDcaTrkMatch(0),
//fDcaTrkMatch(0),
fInvmassLS(0),
fInvmassULS(0),
//fInvmassLSWeightEnhEta(0),
//fInvmassULSWeightEnhEta(0),
//fInvmassLSWeightEnhPi0(0),
//fInvmassULSWeightEnhPi0(0),
//fInvmassLSHijingEta(0),
//fInvmassULSHijingEta(0),
//fInvmassLSHijingPi0(0),
//fInvmassULSHijingPi0(0),
//fInvmassLSHijingPhoton(0),
//fInvmassULSHijingPhoton(0),
//fInvmassLSEnhPhoton(0),
//fInvmassULSEnhPhoton(0),
fULSdcaBelow(0),
fLSdcaBelow(0),
fULSdcaBelowWeight(0),
fLSdcaBelowWeight(0),
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
//fPhotonicDCA(0),
fInclElecDCA(0),
fnSigaftEoPCut(0),
//fnSigaftSysEoPCut(0),
fnSigaftM20EoPCut(0),
//fnSigaftSysM20EoPCut(0),
//fInclElecDCAnoSign(0),
//fElecEoPnoSig(0),
//fInclElecEoPnoShift(0),
//fHadronEoPnoShift(0),
fInclElecEoP(0),
fInclElecEoPNoM20(0),
//fTPCElecEoP(0),
fHadronEoP(0),
fHadronEoPNoM20(0),
fHadronDCA(0),
//fHadronCamDCAHij(0),
//fHadronCamDCA(0),
fPi0Weight(0),
fEtaWeight(0),
//fDWeight(0),
fDWeightNew(0),
fDWeightVar1(0),
fDWeightVar2(0),
fDPlusWeightVar1(0),
fDsWeightVar1(0),
fLcWeightVar1(0),
fLcWeightVar2(0),
//fBWeight(0),
fBWeightNew(0),
fBWeightVar1(0),
fBWeightVar2(0),
fBPlusTauWeight(0),
fB0TauWeight(0),
fBsTauWeight(0),
fDPlusTauWeight(0),
fD0TauWeight(0),
fDsTauWeight(0),

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

fPhotonHijingTagDCA(0),
fEnhPhotonTagDCA(0),

fComboNumWeight(0),
fComboNumNoWeight(0),
fComboDenomWeight(0),
fComboDenomNoWeight(0),
//fDMesonPDG(0),
fD0MesonPt(0),
fD0MesonFromDStarPt(0),
fDPlusMesonPt(0),
fDsMesonPt(0),
fDStarMesonPt(0),
fAllDMesonPt(0),
fLambdaCPt(0),
fD0MesonPtWeight(0),
fLambdaCPtWeight(0),
fDPlusMesonPtWeight(0),
fDsMesonPtWeight(0),
fEtaCPt(0),
fCBaryonPt(0),
fBMesonPt(0),
//fBMesonPtATLAS(0),
//fBPlusPtATLAS(0),
//fBMesonPtCMS(0),
//fBPlusPtCMS(0),
//fBMesonPtLHCb(0),
//fBPlusPtLHCb(0),
fBBaryonPt(0),
fBMesonElecPt(0),
fBBaryonElecPt(0),
//fPromptD0DCAWeight(0),
//fD0FromDStarDCAWeight(0),
//fPromptD0DCANoWeight(0),
//fD0FromDStarDCANoWeight(0),
fNtotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),
fSprsPi0EtaWeightCal(0),
fSprsTemplatesNoWeight(0),
fSprsTemplatesWeight(0),
fSprsTemplatesWeightVar1(0),
fSprsTemplatesWeightVar2(0),
fSprsClosureTest(0),
fSprsClosureTestWeight(0),
fSprsULSdca(0),
fSprsULSdcaWeight(0),
fSprsLSdca(0),
fSprsLSdcaWeight(0),

//fDTemplateWeight(0),
//fDTemplateNoWeight(0),
//fDTemplateWeightNew(0),
//fDTemplateWeightVar1(0),
//fDTemplateWeightVar2(0),

//fBTemplateWeight(0),
//fBTemplateNoWeight(0),
//fBTemplateWeightNew(0),
//fBTemplateWeightVar1(0),
//fBTemplateWeightVar2(0),

fAllElecStack(0),
fHFElecStack(0),
fBElecStack(0),

fAllElecStackDiffPID(0),
fDElecStackDiffPID(0),
fBElecStackDiffPID(0),

fElecTPCTrk(0),
fHFElecTPCTrk(0),
fBElecTPCTrk(0),

fElecAftTrkCuts(0),
fHFElecAftTrkCuts(0),
fBElecAftTrkCuts(0),

fElecAftLooseTrkCuts(0),
fHFElecAftLooseTrkCuts(0),
fBElecAftLooseTrkCuts(0),

//fElecAftLooseTrkCutsDiffPID(0),
//fDElecAftLooseTrkCutsDiffPID(0),
//fBElecAftLooseTrkCutsDiffPID(0),

fElecAftTrkMatch(0),
fHFElecAftTrkMatch(0),
fBElecAftTrkMatch(0),

fElecAftTPCeID(0),
fHFElecAftTPCeID(0),
fBElecAftTPCeID(0),

fElecAftEMCeID(0),
fHFElecAftEMCeID(0),
fBElecAftEMCeID(0),

fElecAftEoP(0),
fHFElecAftEoP(0),
fBElecAftEoP(0),

fElectronSprs(0)
//fvalueElectron(0)
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
fApplyM02Cut(kFALSE),
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
fDCAzCut(3.2),
fMinMass(0.),
fMaxMass(0.1),
fMinEta(-0.6),
fMaxEta(0.6),
fAssoDCAxy(0.25),
fAssoDCAz(1.),
fAssoTPCnCls(80),
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
fTrkMatch(0),
fUseTender(kTRUE),
fApplyHadEoPCut(kTRUE),
fVtxZCut(0.),
fDCAxyCut(2.4),
fFlagULS(kFALSE),
fFlagLS(kFALSE),
fNevents(0),
fVtX(0),
fVtY(0),
fVtZ(0),
fTrkPtB4TC(0),
fDCAxyz(0),
fTrkPt(0),
fTrkP(0),
fTrkClsPhi(0),
fTrkClsEta(0),
fClsPhi(0),
fClsEta(0),
fClsE(0),
fClsEnoTimeCut(0),
//fClsEamDCal(0),
//fClsEamEMCal(0),
//fClsEAll(0),
//fClsEamElecEMC(0),
//fClsEamElecDC(0),
fTrkPhi(0),
fTrkEta(0),
fdEdx(0),
fnSigma(0),
fnSigmaAftTrkMatch(0),
fCentCheck(0),
fTrigCheck(0),
fITSLayerCheck(0),
fEMCTrkMatch(0),
fBasicSprs(0),
fvalueBasic(0),
//fElecDcaB4TrkMatch(0),
//fDcaB4TrkMatch(0),
//fElecDcaTrkMatch(0),
//fDcaTrkMatch(0),
fInvmassLS(0),
fInvmassULS(0),
//fInvmassLSWeightEnhEta(0),
//fInvmassULSWeightEnhEta(0),
//fInvmassLSWeightEnhPi0(0),
//fInvmassULSWeightEnhPi0(0),
//fInvmassLSHijingEta(0),
//fInvmassULSHijingEta(0),
//fInvmassLSHijingPi0(0),
//fInvmassULSHijingPi0(0),
//fInvmassLSHijingPhoton(0),
//fInvmassULSHijingPhoton(0),
//fInvmassLSEnhPhoton(0),
//fInvmassULSEnhPhoton(0),
fULSdcaBelow(0),
fLSdcaBelow(0),
fULSdcaBelowWeight(0),
fLSdcaBelowWeight(0),
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
//fPhotonicDCA(0),
fInclElecDCA(0),
fnSigaftEoPCut(0),
//fnSigaftSysEoPCut(0),
fnSigaftM20EoPCut(0),
//fnSigaftSysM20EoPCut(0),
//fInclElecDCAnoSign(0),
//fElecEoPnoSig(0),
//fInclElecEoPnoShift(0),
//fHadronEoPnoShift(0),
fInclElecEoP(0),
fInclElecEoPNoM20(0),
//fTPCElecEoP(0),
fHadronEoP(0),
fHadronEoPNoM20(0),
fHadronDCA(0),
//fHadronCamDCAHij(0),
//fHadronCamDCA(0),
fPi0Weight(0),
fEtaWeight(0),
//fDWeight(0),
fDWeightNew(0),
fDWeightVar1(0),
fDWeightVar2(0),
fDPlusWeightVar1(0),
fDsWeightVar1(0),
fLcWeightVar1(0),
fLcWeightVar2(0),
//fBWeight(0),
fBWeightNew(0),
fBWeightVar1(0),
fBWeightVar2(0),
fBPlusTauWeight(0),
fB0TauWeight(0),
fBsTauWeight(0),
fDPlusTauWeight(0),
fD0TauWeight(0),
fDsTauWeight(0),

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

fPhotonHijingTagDCA(0),
fEnhPhotonTagDCA(0),

fComboNumWeight(0),
fComboNumNoWeight(0),
fComboDenomWeight(0),
fComboDenomNoWeight(0),
//fDMesonPDG(0),
fD0MesonPt(0),
fD0MesonFromDStarPt(0),
fDPlusMesonPt(0),
fDsMesonPt(0),
fDStarMesonPt(0),
fAllDMesonPt(0),
fLambdaCPt(0),
fD0MesonPtWeight(0),
fLambdaCPtWeight(0),
fDPlusMesonPtWeight(0),
fDsMesonPtWeight(0),
fEtaCPt(0),
fCBaryonPt(0),
fBMesonPt(0),
//fBMesonPtATLAS(0),
//fBPlusPtATLAS(0),
//fBMesonPtCMS(0),
//fBPlusPtCMS(0),
//fBMesonPtLHCb(0),
//fBPlusPtLHCb(0),
fBBaryonPt(0),
fBMesonElecPt(0),
fBBaryonElecPt(0),
//fPromptD0DCAWeight(0),
//fD0FromDStarDCAWeight(0),
//fPromptD0DCANoWeight(0),
//fD0FromDStarDCANoWeight(0),
fNtotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),
fSprsPi0EtaWeightCal(0),
fSprsTemplatesNoWeight(0),
fSprsTemplatesWeight(0),
fSprsTemplatesWeightVar1(0),
fSprsTemplatesWeightVar2(0),
fSprsClosureTest(0),
fSprsClosureTestWeight(0),
fSprsULSdca(0),
fSprsULSdcaWeight(0),
fSprsLSdca(0),
fSprsLSdcaWeight(0),

//fDTemplateWeight(0),
//fDTemplateNoWeight(0),
//fDTemplateWeightNew(0),
//fDTemplateWeightVar1(0),
//fDTemplateWeightVar2(0),
//fBTemplateWeight(0),
//fBTemplateNoWeight(0),
//fBTemplateWeightNew(0),
//fBTemplateWeightVar1(0),
//fBTemplateWeightVar2(0),
fAllElecStack(0),
fHFElecStack(0),
fBElecStack(0),

fAllElecStackDiffPID(0),
fDElecStackDiffPID(0),
fBElecStackDiffPID(0),

fElecTPCTrk(0),
fHFElecTPCTrk(0),
fBElecTPCTrk(0),

fElecAftTrkCuts(0),
fHFElecAftTrkCuts(0),
fBElecAftTrkCuts(0),

fElecAftLooseTrkCuts(0),
fHFElecAftLooseTrkCuts(0),
fBElecAftLooseTrkCuts(0),

//fElecAftLooseTrkCutsDiffPID(0),
//fDElecAftLooseTrkCutsDiffPID(0),
//fBElecAftLooseTrkCutsDiffPID(0),

fElecAftTrkMatch(0),
fHFElecAftTrkMatch(0),
fBElecAftTrkMatch(0),

fElecAftTPCeID(0),
fHFElecAftTPCeID(0),
fBElecAftTPCeID(0),

fElecAftEMCeID(0),
fHFElecAftEMCeID(0),
fBElecAftEMCeID(0),

fElecAftEoP(0),
fHFElecAftEoP(0),
fBElecAftEoP(0),
fElectronSprs(0)
//fvalueElectron(0)
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
        delete []fvalueBasic;
        /*if(fBPlusTauWeight)    {delete fBPlusTauWeight; fBPlusTauWeight=0;}
        if(fB0TauWeight)    {delete fB0TauWeight; fB0TauWeight=0;}
        if(fBsTauWeight)    {delete fBsTauWeight; fBsTauWeight=0;}*/
        //delete []fvalueElectron;
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
    fNevents = new TH1F("fNevents", "No. of Events; Counts", 6,-0.5,5.5);
    fOutputList->Add(fNevents);
    fNevents->GetXaxis()->SetBinLabel(1,"All");
    fNevents->GetXaxis()->SetBinLabel(2,">2 Trks");
    fNevents->GetXaxis()->SetBinLabel(3,">2 Trks, Vtx_{z}<10cm");
    fNevents->GetXaxis()->SetBinLabel(4,"Vtx_{z}<10cm");
    fNevents->GetXaxis()->SetBinLabel(5,"Vtx_{z}<10cm, Pile-up cuts");
    fNevents->GetXaxis()->SetBinLabel(6,">2 Trks, Vtx_{z}<10cm, Pile-Up cuts");
    
    fVtX = new TH1F("fVtX","X Vertex Position;Vtx_{X};Counts",50,-5,5);
    fOutputList->Add(fVtX);
    
    fVtY = new TH1F("fVtY","Y Vertex Position;Vtx_{Y};Counts",50,-5,5);
    fOutputList->Add(fVtY);
    
    fVtZ = new TH1F("fVtZ","Z Vertex Position;Vtx_{Z};Counts",100,-10,10);
    fOutputList->Add(fVtZ);
    
    fTrkPtB4TC = new TH1F("fTrkPtB4TC","Track p_{T} Distribution before track cuts;p_{T} (GeV/c);Counts",100,0,50);
    fOutputList->Add(fTrkPtB4TC);
    
    fDCAxyz = new TH2F("fDCAxyz","Track DCAz vs. DCAxy;DCAxy;DCAz",40,0,4,40,0,4);
    fOutputList->Add(fDCAxyz);
    
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
    
    fClsE = new TH1F("fClsE","Cluster Energy after matching;Cluster E;Counts",250,0.,50);
    fOutputList->Add(fClsE);
    
    fClsEnoTimeCut = new TH1F("fClsEnoTimeCut","Cluster Energy after matching, no cluster time cut;Cluster E;Counts",250,0.,50);
    fOutputList->Add(fClsEnoTimeCut);
    
    /*if(fFlagClsTypeDCAL && !fFlagClsTypeEMC){
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
    }*/
        
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
    
    fITSLayerCheck = new TH1F("fITSLayerCheck", "No. of tracks; Counts",3,-0.5,2.5);
    fOutputList->Add(fITSLayerCheck);
    fITSLayerCheck->GetXaxis()->SetBinLabel(1,"kAny");
    fITSLayerCheck->GetXaxis()->SetBinLabel(2,"kFirst");
    fITSLayerCheck->GetXaxis()->SetBinLabel(3,"kBoth");
    
    fEMCTrkMatch = new TH2F("fEMCTrkMatch","EMCal cluster distance from closest track", 100, -0.3, 0.3,100,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch);
    
    Int_t binBasic[7] = {60,nDCAbins,100,100,2,2,2}; //pT, DCA, phi, eta, charge, elec bool, trk match bool
    Double_t xminBasic[7] = {0,-0.2,0,-1.5,-2,-1,-1};
    Double_t xmaxBasic[7] = {30,0.2,6.3,1.5,2,2,2};
    fBasicSprs = new THnSparseD("fBasicSprs","Sprs with eta/phi info;p_{T};DCA;#phi;#eta;Charge;Elec Bool;Trk Match Bool;",7,binBasic,xminBasic,xmaxBasic);
    fBasicSprs->GetAxis(5)->SetBinLabel(1,"Not elec");
    fBasicSprs->GetAxis(5)->SetBinLabel(2,"-1<n#sigma<3");
    fBasicSprs->GetAxis(6)->SetBinLabel(1,"Not matched");
    fBasicSprs->GetAxis(6)->SetBinLabel(2,"Matched");
    fOutputList->Add(fBasicSprs);
        
    fvalueBasic = new Double_t[7];
    
    /*fElecDcaB4TrkMatch = new TH3F("fElecDcaB4TrkMatch","DCA of e- before trk match;p_{T}(GeV/c);DCA;#phi;", 60,0,30., nDCAbins,-0.2,0.2, 100,0,6.3);
    fOutputList->Add(fElecDcaB4TrkMatch);
    
    fDcaB4TrkMatch = new TH3F("fDcaB4TrkMatch","DCA before trk match;p_{T}(GeV/c);DCA;#phi;", 60,0,30., nDCAbins,-0.2,0.2, 100,0,6.3);
    fOutputList->Add(fDcaB4TrkMatch);
    
    fElecDcaTrkMatch = new TH3F("fElecDcaTrkMatch","DCA of e- after trk match;p_{T}(GeV/c);DCA;#phi;", 60,0,30., nDCAbins,-0.2,0.2, 100,0,6.3);
    fOutputList->Add(fElecDcaTrkMatch);
    
    fDcaTrkMatch = new TH3F("fDcaTrkMatch","DCA after trk match;p_{T}(GeV/c);DCA;#phi;", 60,0,30., nDCAbins,-0.2,0.2, 100,0,6.3);
    fOutputList->Add(fDcaTrkMatch);*/
    
    fInvmassLS = new TH1F("fInvmassLS", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassLS);
    fInvmassULS = new TH1F("fInvmassULS", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
    fOutputList->Add(fInvmassULS);
    
    if (fFlagFillMCHistos) {
        /*fInvmassLSWeightEnhEta = new TH1F("fInvmassLSWeightEnhEta", "Inv mass of Weighted Enh Eta LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 100,0,1.0);
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
        fOutputList->Add(fInvmassULSEnhPhoton);*/
    }
    
    fULSdcaBelow = new TH2F("fULSdcaBelow","ULS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
    fULSdcaBelow->Sumw2();
    fOutputList->Add(fULSdcaBelow);
    
    fLSdcaBelow = new TH2F("fLSdcaBelow","LS Elec DCA m<0.1GeV/c^{2}; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
    fLSdcaBelow->Sumw2();
    fOutputList->Add(fLSdcaBelow);
    
    //if (fFlagFillMCHistos) {
        fULSdcaBelowWeight = new TH2F("fULSdcaBelowWeight","ULS Elec DCA m<0.1GeV/c^{2} w/ pi0+eta weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fULSdcaBelowWeight->Sumw2();
        fOutputList->Add(fULSdcaBelowWeight);
        
        fLSdcaBelowWeight = new TH2F("fLSdcaBelowWeight","LS Elec DCA m<0.1GeV/c^{2} w/ pi0+eta weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fLSdcaBelowWeight->Sumw2();
        fOutputList->Add(fLSdcaBelowWeight);
    //}
    
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
    
        //fPhotonicDCA = new TH2F("fPhotonicDCA","Photonic DCA using MC PID; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        //fOutputList->Add(fPhotonicDCA);
    }
    fInclElecDCA = new TH2F("fInclElecDCA","Incl Elec DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
    fOutputList->Add(fInclElecDCA);
    
    fnSigaftEoPCut = new TH2F("fnSigaftEoPCut","nSig with 0.9<E/p<1.2; p_{T}(GeV/c); nSigma; counts;", 60,0,30., 160,-8,8);
    fOutputList->Add(fnSigaftEoPCut);
    
    //fnSigaftSysEoPCut = new TH2F("fnSigaftSysEoPCut","nSig with Sys E/p cut; p_{T}(GeV/c); nSigma; counts;", 60,0,30., 160,-8,8);
    //fOutputList->Add(fnSigaftSysEoPCut);
    
    fnSigaftM20EoPCut = new TH2F("fnSigaftM20EoPCut","nSig with 0.9<E/p<1.2, 0.01<M20<0.35; p_{T}(GeV/c); nSigma; counts;", 60,0,30., 160,-8,8);
    fOutputList->Add(fnSigaftM20EoPCut);
    
    //fnSigaftSysM20EoPCut = new TH2F("fnSigaftSysM20EoPCut","nSig with systematic E/p and M20 cut; p_{T}(GeV/c); nSigma; counts;", 60,0,30., 160,-8,8);
    //fOutputList->Add(fnSigaftSysM20EoPCut);
    
    //fInclElecDCAnoSign = new TH2F("fInclElecDCAnoSign","Incl Elec DCA (no Charge); p_{T}(GeV/c); DCAxMagField; counts;", 60,0,30., nDCAbins,-0.2,0.2);
    //fOutputList->Add(fInclElecDCAnoSign);
    
    //fElecEoPnoSig = new TH2F("fElecEoPnoSig","Elec E/p, no nSig cut; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    //fOutputList->Add(fElecEoPnoSig);
    
    //fInclElecEoPnoShift = new TH2F("fInclElecEoPnoShift","Incl Elec E/p no Shift; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    //fOutputList->Add(fInclElecEoPnoShift);
    
    //fHadronEoPnoShift = new TH2F("fHadronEoPnoShift","Incl Elec E/p no Shift; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    //fOutputList->Add(fHadronEoPnoShift);
    
    fInclElecEoP = new TH2F("fInclElecEoP","Incl Elec E/p; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fInclElecEoP);
    
    fInclElecEoPNoM20 = new TH2F("fInclElecEoPNoM20","Incl Elec E/p, No M20 cut; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fInclElecEoPNoM20);
    
    //fTPCElecEoP = new TH2F("fTPCElecEoP","Elec E/p, -0.1<nsig<3; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    //fOutputList->Add(fTPCElecEoP);
    
    fHadronEoP = new TH2F("fHadronEoP","Unscaled Hadron E/p; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fHadronEoP);
    
    fHadronEoPNoM20 = new TH2F("fHadronEoPNoM20","Unscaled Hadron E/p, No M20 Cut; p_{T}(GeV/c); E/p; counts;", 60,0,30., 100,0.,2.);
    fOutputList->Add(fHadronEoPNoM20);
    
    fHadronDCA = new TH2F("fHadronDCA","Unscaled Hadron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
    fOutputList->Add(fHadronDCA);
    
    /*if (fFlagFillMCHistos) {
        fHadronCamDCAHij = new TH2F("fHadronCamDCAHij","Unscaled Hadron DCA, no E/p cut, Hijing; p_{T}(GeV/c); DCAxMagField; counts;", 60,0,30., 400,-0.2,0.2);
        fOutputList->Add(fHadronCamDCAHij);
    }
    
    fHadronCamDCA = new TH2F("fHadronCamDCA","Unscaled Hadron DCA, no E/p cut, Enh+Hij; p_{T}(GeV/c); DCAxMagField; counts;", 60,0,30., 400,-0.2,0.2);
    fOutputList->Add(fHadronCamDCA);*/
    
    
    fPi0Weight = new TF1("fPi0Weight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    fEtaWeight = new TF1("fEtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    //fPi0EtaWeight = new TF1("fPi0EtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    
    //Deepa's Weight
    //fPi0Weight->SetParameters(8.96715e+02,-1.77016e-01,2.47560e-03,1.50783e+00,4.41819e+00);
    //fEtaWeight->SetParameters(4.18393e+02,-5.20750e-02,-3.11517e-04,2.04739e+00,5.30788e+00);
    
    //Default 0
    fPi0Weight->SetParameters(0,1,1,1,1);
    fEtaWeight->SetParameters(0,1,1,1,1);
    
    //0-10%
    if (fCentralityMin==0 && fCentralityMax==10 && fApplyCentrality){
        fPi0Weight->SetParameters(3.10605e+04,-2.88504e-01,4.90349e-03,4.82424e-01,3.88146e+00);
        fEtaWeight->SetParameters(1.48495e+03,-2.20500e-01,3.35819e-03,1.14265e+00,4.06488e+00);
        //fPi0EtaWeight->SetParameters(3.68627e+03,-2.15364e-01,3.26753e-03,1.05213e+00,4.20391e+00);
    }
    //30-50%
    if (fCentralityMin==30 && fCentralityMax==50){
        fPi0Weight->SetParameters(21889.2, -0.296368, 0.00520277, 0.449904, 3.92565);
        fEtaWeight->SetParameters(1554.85, -0.243678, 0.00379987, 0.827357, 3.9573);
    }
    fOutputList->Add(fPi0Weight);
    fOutputList->Add(fEtaWeight);
        
    //D Meson pt weighting
    Int_t nbins = 13;
    Double_t xbins[14] = {1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,16.,24.,36.,50.};
    //Double_t err[13] = {};
    //fDWeight = new TH1F("fDWeight","D^{0}_data/AllD_MC;p_{T} (GeV/c);Weight;",nbins,xbins);
    fDWeightNew = new TH1F("fDWeightNew","D^{0}_data/AllD_MCNew;p_{T} (GeV/c);Weight;",nbins,xbins);
    fDWeightVar1 = new TH1F("fDWeightVar1","D^{0}_data/AllD_MC;p_{T} (GeV/c);Weight Var1;",nbins,xbins);
    fDWeightVar2 = new TH1F("fDWeightVar2","D^{0}_data/AllD_MC;p_{T} (GeV/c);Weight Var2;",nbins,xbins);
    
    fDPlusWeightVar1 = new TH1F("fDPlusWeightVar1","(D^{+}/D^{0})_data*(D^{0}_data/D^{+}_{MC});p_{T} (GeV/c);Weight Var1;",nbins,xbins);
    fDsWeightVar1 = new TH1F("fDsWeightVar1","(D^{s}/D^{0})_data*(D^{0}_data/D^{s}_{MC});p_{T} (GeV/c);Weight Var1;",nbins,xbins);
    
    fLcWeightVar1 = new TH1F("fLcWeightVar1","Lc weight, Lc/D0 from model;p_{T} (GeV/c);Weight Var1;",nbins,xbins);
    fLcWeightVar2 = new TH1F("fLcWeightVar2","Lc weight, Lc/D0 from data;p_{T} (GeV/c);Weight Var2;",nbins,xbins);
    //Double_t ratio[13] = {2.03552,1.0201,0.45925,0.211574,0.11987,0.0898116,0.0631282,0.0546798,0.0477205,0.0410021,0.0307936,0.0398483,0.0175335};
    //Double_t err[13] = {0.541651,0.146443,0.0498454,0.024907,0.01438,0.0107908,0.00848616,0.0061723,0.00587082,0.00566712,0.00597994,0.00811015,0.00693105};
    
    if (fCentralityMin==0 && fCentralityMax==10 && fApplyCentrality) {
        Double_t ratio[13] = {0.106888,0.0650239,0.0343858,0.0172579,0.00957876,0.00640323,0.00399907,0.00269269,0.00163078,0.000942387,0.000441093,0.000353811,0.000143011};
        Double_t err[13] = {0.0284416,0.0093333,0.00373075,0.00203067,0.00114824,0.000768388,0.00053676,0.000303334,0.000199878,0.000129785,8.53822e-05,7.13313e-05,5.61316e-05};
        //Double_t ratioNew[13] = {0.197449,0.118714,0.0627949,0.0321233,0.0182153,0.0124903,0.00801369,0.00553768,0.00340667,0.00193131,0.00089526,0.000678224,0.00026223};
        Double_t ratioNew[13] = {0.382299,0.223983,0.116971,0.0573891,0.0292338,0.0192381,0.011778,0.00863768,0.00534462,0.00301279,0.00159646,0.00131757,0.000523965};
        Double_t ratioVar1[13] = {0.249988,0.132914,0.0673368,0.0340132,0.0189431,0.01274,0.00801369,0.00543372,0.00326751,0.00179833,0.000779735,0.000564287,0.000159311};
        Double_t ratioVar2[13] = {0.144911,0.104514,0.058253,0.0302335,0.0174875,0.0122405,0.00801369,0.00564164,0.00354583,0.00206429,0.00101078,0.00079216,0.000365149};
        Double_t wLcVar1[13] = {1.57532,1.46238,0.948202,0.497431,0.250927,0.151636,0.0827256,0.0487993,0.0218917,0.0076259,0.00180605,0.00055039,8.79344e-05};
        Double_t wLcVar2[13] = {3.47398,1.72259,0.839403,0.514328,0.261815,0.145702,0.0895206,0.0380316,0.0234082,0.0051148,0.00266179,0.00049151,0};
        for (int idata=1; idata<14; idata++) {
            //fDWeight->SetBinContent(idata,ratio[idata-1]);
            //fDWeight->SetBinError(idata,err[idata-1]);
            
            fDWeightNew->SetBinContent(idata,ratioNew[idata-1]);
            fDWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
            fDWeightVar2->SetBinContent(idata,ratioVar2[idata-1]);
            fLcWeightVar1->SetBinContent(idata,wLcVar1[idata-1]);
            fLcWeightVar2->SetBinContent(idata,wLcVar2[idata-1]);
            fDPlusWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
            fDsWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
        }
    }else if (fCentralityMin==30 && fCentralityMax==50 && fApplyCentrality) {
        Double_t ratio[13] = {0.079428,0.0402934,0.0258836,0.0165168,0.0117076,0.00807683,0.00545914,0.00413535,0.00218055,0.00147282,0.000578039,0.000286482,0.000286482};
        Double_t err[13] = {0.0112672,0.00212923,0.000929744,0.000665396,0.000505178,0.00041285,0.000325958,0.000210603,0.000152176,0.000110271,6.60032e-05,5.6821e-05,0};
        //Double_t ratioNew[13] = {0.0611667,0.0299809,0.0186145,0.0114545,0.00789462,0.0053161,0.00355211,0.002672,0.00141544,0.000970815,0.000391017,0.000201985,0.000107469};
        Double_t ratioNew[13] = {0.171312,0.0872016,0.0553106,0.034332,0.0236477,0.0159165,0.0106375,0.00797653,0.00422176,0.00292325,0.00118981,0.000624289,0.0003415285};
        Double_t ratioVar1[13] = {0.0705071,0.031752,0.019312,0.0118197,0.00807804,0.00538886,0.00355211,0.0026429,0.00137442,0.000931529,0.000360442,0.000168603,0.0000726835};
        Double_t ratioVar2[13] = {0.0518263,0.0282098,0.0179171,0.0110893,0.00771121,0.00524334,0.00355211,0.00270111,0.00145645,0.0010101,0.000421592,0.000235367,0.0001422545};
        Double_t wLcVar1[13] = {2.02089,0.889089,0.490217,0.266978,0.158602,0.0931538,0.0536931,0.0323656,0.0131682,0.00593832,0.00123755,0.000273619,0.};
        Double_t wLcVar2[13] = {0.621365,0.207712,0.0824139,0.0377918,0.0179345,0.0140968,0.00861207,0.00516433,0.00249769,0.0017025,0.000517206,6.8216e-05,0.};
        Double_t wDPlusVar1[13] = {0.125876,0.0642961,0.0414463,0.026278,0.01833,0.0124752,0.00830149,0.00628302,0.00333208,0.00227467,0.000910711,0.000469747,0.000028783};
        Double_t wDsVar1[13] = {0.317507,0.147254,0.0890879,0.0533869,0.0365636,0.0245759,0.0163403,0.0121869,0.00646731,0.00444424,0.00180416,0.000929189,0.000054218};
        /*Double_t ratio[13] = {0.566977,0.233989,0.0909109,0.0346338,0.0155742,0.00734675,0.00362088,0.00189595,0.000670163,0.000284606,5.13181e-05,8.84883e-06,8.84883e-06};
        Double_t err[13] = {0.0804276,0.0123637,0.00326455,0.0013947,0.000671652,0.000375317,0.000216074,9.65006e-05,4.67481e-05,2.13019e-05,5.85888e-06,1.75491e-06,0};
        Double_t ratioNew[13] = {0.566977,0.233989,0.0909109,0.0346338,0.0155742,0.00734675,0.00362088,0.00189595,0.000670163,0.000284606,5.13181e-05,8.84883e-06,8.84883e-06};
        Double_t ratioVar1[13] = {0.566977,0.233989,0.0909109,0.0346338,0.0155742,0.00734675,0.00362088,0.00189595,0.000670163,0.000284606,5.13181e-05,8.84883e-06,8.84883e-06};
        Double_t ratioVar2[13] = {0.566977,0.233989,0.0909109,0.0346338,0.0155742,0.00734675,0.00362088,0.00189595,0.000670163,0.000284606,5.13181e-05,8.84883e-06,8.84883e-06};*/
        for (int idata=1; idata<14; idata++) {
            //fDWeight->SetBinContent(idata,ratio[idata-1]);
            //fDWeight->SetBinError(idata,err[idata-1]);
            fDWeightNew->SetBinContent(idata,ratioNew[idata-1]);
            fDWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
            fDWeightVar2->SetBinContent(idata,ratioVar2[idata-1]);
            fLcWeightVar1->SetBinContent(idata,wLcVar1[idata-1]);
            fLcWeightVar2->SetBinContent(idata,wLcVar2[idata-1]);
            fDPlusWeightVar1->SetBinContent(idata,wDPlusVar1[idata-1]);
            fDsWeightVar1->SetBinContent(idata,wDsVar1[idata-1]);
        }
    }else{
        for (int idata=1; idata<14; idata++) {
            //fDWeight->SetBinContent(idata,ratio[idata-1]);
            //fDWeight->SetBinError(idata,err[idata-1]);
            fDWeightNew->SetBinContent(idata,1);
            fDWeightVar1->SetBinContent(idata,1);
            fDWeightVar2->SetBinContent(idata,1);
            fLcWeightVar1->SetBinContent(idata,1);
            fLcWeightVar2->SetBinContent(idata,1);
            fDPlusWeightVar1->SetBinContent(idata,1);
            fDsWeightVar1->SetBinContent(idata,1);
        }
    }
    
    //fDWeight->Sumw2();
    //fOutputList->Add(fDWeight);
    fOutputList->Add(fDWeightNew);
    fOutputList->Add(fDWeightVar1);
    fOutputList->Add(fDWeightVar2);
    fOutputList->Add(fLcWeightVar1);
    fOutputList->Add(fLcWeightVar2);
    fOutputList->Add(fDPlusWeightVar1);
    fOutputList->Add(fDsWeightVar1);
    
    
    //B Meson pt weighting
    Int_t nbinsB = 250;
    Double_t xbinsB[251] = {0.,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8,8,8.2,8.4,8.6,8.8,9,9.2,9.4,9.6,9.8,10,10.2,10.4,10.6,10.8,11,11.2,11.4,11.6,11.8,12,12.2,12.4,12.6,12.8,13,13.2,13.4,13.6,13.8,14,14.2,14.4,14.6,14.8,15,15.2,15.4,15.6,15.8,16,16.2,16.4,16.6,16.8,17,17.2,17.4,17.6,17.8,18,18.2,18.4,18.6,18.8,19,19.2,19.4,19.6,19.8,20,20.2,20.4,20.6,20.8,21,21.2,21.4,21.6,21.8,22,22.2,22.4,22.6,22.8,23,23.2,23.4,23.6,23.8,24,24.2,24.4,24.6,24.8,25,25.2,25.4,25.6,25.8,26,26.2,26.4,26.6,26.8,27,27.2,27.4,27.6,27.8,28,28.2,28.4,28.6,28.8,29,29.2,29.4,29.6,29.8,30,30.2,30.4,30.6,30.8,31,31.2,31.4,31.6,31.8,32,32.2,32.4,32.6,32.8,33,33.2,33.4,33.6,33.8,34,34.2,34.4,34.6,34.8,35,35.2,35.4,35.6,35.8,36,36.2,36.4,36.6,36.8,37,37.2,37.4,37.6,37.8,38,38.2,38.4,38.6,38.8,39,39.2,39.4,39.6,39.8,40,40.2,40.4,40.6,40.8,41,41.2,41.4,41.6,41.8,42,42.2,42.4,42.6,42.8,43,43.2,43.4,43.6,43.8,44,44.2,44.4,44.6,44.8,45,45.2,45.4,45.6,45.8,46,46.2,46.4,46.6,46.8,47,47.2,47.4,47.6,47.8,48,48.2,48.4,48.6,48.8,49,49.2,49.4,49.6,49.8,50.};
    //fBWeight = new TH1F("fBWeight","TAMU RAA x FONLL/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    fBWeightNew = new TH1F("fBWeightNew","TAMU RAA x FONLL(New)/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    fBWeightVar1 = new TH1F("fBWeightVar1","TAMU RAA(Max) x FONLL(New)/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    fBWeightVar2 = new TH1F("fBWeightVar2","TAMU RAA(Min) x FONLL(New)/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    
    if (fCentralityMin==0 && fCentralityMax==10 && fApplyCentrality) {
        Double_t ratioB[250] = {0.498558,0.62782,0.672398,0.718151,0.731677,0.734297,0.770069,0.776389,0.809054,0.825417,0.859501,0.874636,0.898427,0.923358,0.92293,0.926125,0.926343,0.916683,0.917184,0.911809,0.886661,0.868603,0.852465,0.837679,0.834915,0.813126,0.776888,0.762151,0.745226,0.726064,0.689477,0.675392,0.657326,0.639663,0.610833,0.599204,0.574208,0.551114,0.534405,0.508084,0.503657,0.468409,0.459998,0.447575,0.429903,0.413056,0.396865,0.385802,0.36568,0.358106,0.351127,0.342732,0.324592,0.325566,0.317839,0.306017,0.292589,0.292114,0.28069,0.274267,0.268707,0.260392,0.261063,0.248518,0.245039,0.237828,0.226948,0.221301,0.216985,0.214584,0.210782,0.197669,0.195025,0.190053,0.186186,0.179784,0.171765,0.169136,0.159454,0.161479,0.157118,0.153996,0.14812,0.141091,0.138405,0.137488,0.134021,0.127508,0.126735,0.119841,0.117923,0.11417,0.109935,0.110609,0.106261,0.102496,0.0988583,0.0988566,0.0962846,0.0971962,0.0912081,0.0906153,0.0909837,0.089602,0.0851854,0.0829322,0.0854181,0.0793186,0.0828916,0.079228,0.0765657,0.0755108,0.0746448,0.0749831,0.0726684,0.0741421,0.0725614,0.0725931,0.0708997,0.0691765,0.0675163,0.0638974,0.0636498,0.0696325,0.0621972,0.0654069,0.0613243,0.0604615,0.0588438,0.0602627,0.0602133,0.0562468,0.058719,0.0581855,0.0582384,0.0562759,0.053968,0.053529,0.0549016,0.0524293,0.0523908,0.0529672,0.0515757,0.0493422,0.0481627,0.0482736,0.0459697,0.0458283,0.0453141,0.0433415,0.045522,0.043535,0.0411514,0.0456961,0.0441323,0.0430352,0.0428453,0.0426922,0.0436595,0.040605,0.0381453,0.03843,0.0410627,0.0374946,0.0391381,0.0379844,0.0367569,0.036893,0.0372399,0.0355532,0.0336221,0.0344674,0.0330855,0.0332478,0.0308431,0.0318583,0.0313297,0.03175,0.0314691,0.0320693,0.0312073,0.0297707,0.0285189,0.0289505,0.0291229,0.0291655,0.0305036,0.0281665,0.0276302,0.0291551,0.026976,0.027417,0.027251,0.0245767,0.0246722,0.024308,0.0247152,0.0273125,0.0242685,0.0232144,0.0220058,0.023307,0.0224514,0.0219948,0.0208518,0.0216804,0.0215576,0.0195631,0.0195621,0.0211566,0.0185032,0.0196948,0.0190276,0.0193687,0.0191364,0.0195232,0.0190942,0.0177003,0.0186524,0.0180672,0.017813,0.0173329,0.0163168,0.0163818,0.0158271,0.0169163,0.0164239,0.0157851,0.0169332,0.0180851,0.0152819,0.0149112,0.0141721,0.015732,0.0163534,0.014337,0.0143059,0.014786,0.0140978,0.0145358,0.0145749,0.0134334,0.0130966,0.0138116,0.013262,0.0143179,0.0132258,0.0135375,0.0132521,0.0132662};
        Double_t ratioBNew[250] = {0.527854,0.676468,0.697991,0.711774,0.722505,0.744184,0.771669,0.791643,0.815676,0.842682,0.869755,0.886323,0.900165,0.921397,0.927047,0.93348,0.936192,0.926124,0.917184,0.909077,0.889903,0.876369,0.858462,0.838199,0.821877,0.805409,0.785763,0.762809,0.741829,0.72182,0.704218,0.683978,0.660144,0.639746,0.617747,0.596334,0.573373,0.553575,0.533483,0.516636,0.49255,0.475004,0.460341,0.441742,0.4269,0.412142,0.40238,0.387698,0.37121,0.360204,0.352943,0.341443,0.332603,0.323606,0.315956,0.305826,0.295788,0.291485,0.28239,0.275661,0.268983,0.260012,0.25416,0.24767,0.24152,0.234847,0.229889,0.2255,0.216294,0.21002,0.204699,0.200694,0.194504,0.189242,0.185419,0.179002,0.173614,0.169286,0.162982,0.158871,0.15605,0.149979,0.146925,0.143822,0.138876,0.135658,0.131925,0.128432,0.125932,0.121191,0.11921,0.116162,0.111859,0.11006,0.106977,0.105162,0.103074,0.100297,0.0993633,0.0962979,0.0940254,0.0921904,0.0906611,0.0891579,0.0880537,0.0863541,0.0855137,0.0817595,0.0811291,0.0799165,0.0793575,0.0781609,0.0766887,0.0750305,0.0738189,0.073269,0.072447,0.0706597,0.0707376,0.0681811,0.0680168,0.0670508,0.065782,0.0663052,0.0645379,0.0638748,0.0625307,0.0614706,0.0613937,0.0600465,0.0599109,0.0578851,0.0589849,0.0577687,0.0577089,0.0564293,0.0548461,0.0541342,0.0537008,0.0534945,0.0522089,0.0527365,0.0517209,0.0505692,0.0501077,0.0491331,0.0493484,0.0485754,0.0477572,0.0471697,0.0461767,0.0455698,0.0463311,0.0454191,0.0441808,0.0445879,0.0432137,0.043078,0.0429534,0.041554,0.0404588,0.040246,0.0399505,0.0400038,0.0393302,0.0378221,0.0377077,0.0375886,0.0364395,0.0361215,0.0359106,0.0354513,0.035101,0.0343332,0.0340229,0.033936,0.0329639,0.0324529,0.0320335,0.0311993,0.0318253,0.0314773,0.0295871,0.0297416,0.0294255,0.0293431,0.0287639,0.028382,0.0275853,0.0278151,0.0270327,0.026975,0.0267576,0.0258747,0.0254334,0.0253068,0.0250752,0.0246926,0.0236031,0.0242901,0.0232739,0.0232222,0.0229169,0.0232431,0.0216698,0.0220481,0.0212689,0.0212692,0.0212703,0.0207071,0.0207659,0.0199056,0.019779,0.0194439,0.0190936,0.0189077,0.0188995,0.017954,0.0178651,0.0179176,0.0177795,0.0169843,0.0170956,0.0171893,0.0168782,0.0165971,0.017014,0.0161936,0.0161965,0.0158903,0.0157885,0.0152633,0.0151995,0.0154153,0.0155522,0.0153605,0.0147381,0.0147705,0.0144276,0.0145398,0.0142471,0.0144064,0.0138917,0.0140977,0.014126,0.0138707,0.0134346,0.0137999,0.013464,0.0130542};
        Double_t ratioBVar1[250] = {0.54281,0.695606,0.717742,0.731959,0.743084,0.765527,0.794016,0.814866,0.839998,0.86832,0.896867,0.914756,0.930023,0.953155,0.96042,0.968762,0.97354,0.965327,0.958592,0.953066,0.936273,0.925755,0.910985,0.894068,0.881733,0.869652,0.854532,0.83615,0.820235,0.805694,0.794136,0.779845,0.761552,0.747232,0.730978,0.715231,0.697306,0.682818,0.667487,0.655676,0.633963,0.619851,0.608774,0.591691,0.578793,0.565203,0.557723,0.542693,0.524328,0.512981,0.506383,0.493153,0.483236,0.472629,0.463576,0.450508,0.437227,0.432141,0.419709,0.410571,0.401326,0.388495,0.380185,0.370806,0.36184,0.352007,0.344679,0.338147,0.324347,0.314906,0.306865,0.300774,0.291389,0.283384,0.277522,0.267772,0.259559,0.252931,0.243351,0.237051,0.232676,0.22346,0.218746,0.213963,0.206445,0.201502,0.195801,0.190463,0.186603,0.179431,0.176352,0.171701,0.165204,0.162412,0.15773,0.154925,0.151721,0.147509,0.146013,0.14139,0.137937,0.135132,0.132778,0.130467,0.128743,0.126153,0.12482,0.119241,0.118222,0.116358,0.115448,0.113612,0.11138,0.108881,0.107034,0.106149,0.104871,0.1022,0.102228,0.0984524,0.0981344,0.0966614,0.0947546,0.0954302,0.0928108,0.0917825,0.0897781,0.0881844,0.0880027,0.0860019,0.0857383,0.0827723,0.0842771,0.0824731,0.0823216,0.0804317,0.0781125,0.0770371,0.0763594,0.0760055,0.07412,0.0748097,0.0733108,0.0716217,0.070912,0.0694779,0.0697275,0.0685814,0.0673733,0.0664924,0.0650417,0.0641368,0.0651575,0.0638252,0.0620368,0.0625599,0.060585,0.060348,0.0601271,0.0581234,0.0565479,0.0562073,0.0557518,0.0557834,0.0548022,0.0526606,0.0524614,0.0522559,0.05062,0.0501403,0.0498098,0.0491355,0.0486133,0.0475141,0.0470493,0.0468939,0.0455165,0.0447775,0.0441656,0.0429835,0.0438132,0.043302,0.0406715,0.0408535,0.0403895,0.0402465,0.0394231,0.038871,0.037752,0.0380386,0.0369416,0.0368358,0.0365122,0.0352817,0.0346547,0.034457,0.034117,0.0335721,0.0320675,0.0329771,0.0315746,0.0314818,0.0310455,0.0314648,0.029314,0.0298044,0.0287306,0.0287104,0.0286915,0.0279118,0.0279712,0.0267935,0.0266041,0.026135,0.025646,0.0253783,0.0253495,0.0240644,0.0239284,0.023982,0.0237805,0.0227009,0.0228338,0.022943,0.0225121,0.0221218,0.0226618,0.0215541,0.0215431,0.0211211,0.0209714,0.0202598,0.0201613,0.0204336,0.0206009,0.0203331,0.0194958,0.0195254,0.0190591,0.0191943,0.0187951,0.0189924,0.0183014,0.0185603,0.018585,0.0182369,0.0176516,0.0181193,0.0176665,0.0171173};
        Double_t ratioBVar2[250] = {0.512898,0.65733,0.67824,0.69159,0.701926,0.722841,0.749322,0.768421,0.791354,0.817044,0.842643,0.85789,0.870306,0.889639,0.893675,0.898198,0.898843,0.886921,0.875776,0.865088,0.843533,0.826984,0.80594,0.782329,0.762022,0.741167,0.716994,0.689468,0.663422,0.637946,0.6143,0.588111,0.558735,0.53226,0.504516,0.477438,0.449439,0.424332,0.399479,0.377596,0.351136,0.330158,0.311908,0.291794,0.275008,0.259082,0.247036,0.232703,0.218092,0.207428,0.199503,0.189734,0.181971,0.174584,0.168336,0.161143,0.154349,0.150828,0.145071,0.14075,0.13664,0.131529,0.128135,0.124534,0.121201,0.117687,0.1151,0.112853,0.108241,0.105134,0.102532,0.100614,0.0976182,0.0951008,0.0933164,0.0902324,0.0876688,0.0856413,0.0826121,0.080691,0.0794234,0.0764967,0.0751029,0.0736805,0.0713075,0.0698141,0.0680495,0.0664012,0.0652605,0.0629505,0.0620669,0.0606222,0.0585142,0.0577086,0.0562239,0.0554,0.0544271,0.0530847,0.0527133,0.051206,0.0501136,0.0492492,0.0485438,0.0478486,0.0473641,0.0465557,0.0462072,0.0442785,0.0440359,0.0434749,0.0432671,0.0427093,0.0419975,0.0411799,0.0406037,0.0403892,0.0400228,0.0391198,0.0392472,0.0379098,0.0378991,0.0374402,0.0368093,0.0371801,0.0362649,0.0359672,0.0352834,0.0347569,0.0347848,0.0340911,0.0340835,0.0329978,0.0336927,0.0330643,0.0330962,0.0324268,0.0315796,0.0312313,0.0310422,0.0309835,0.0302978,0.0306634,0.030131,0.0295167,0.0293034,0.0287882,0.0289693,0.0285695,0.0281411,0.027847,0.0273117,0.0270028,0.0275047,0.027013,0.0263247,0.0266158,0.0258425,0.025808,0.0257798,0.0249847,0.0243697,0.0242847,0.0241492,0.0242241,0.0238581,0.0229835,0.022954,0.0229212,0.022259,0.0221028,0.0220115,0.021767,0.0215887,0.0211522,0.0209965,0.0209781,0.0204113,0.0201284,0.0199013,0.0194151,0.0198373,0.0196527,0.0185027,0.0186296,0.0184616,0.0184396,0.0181047,0.017893,0.0174185,0.0175915,0.0171238,0.0171143,0.017003,0.0164678,0.0162122,0.0161565,0.0160335,0.0158132,0.0151387,0.0156032,0.0149732,0.0149626,0.0147882,0.0150214,0.0140256,0.0142919,0.0138073,0.013828,0.0138492,0.0135023,0.0135606,0.0130178,0.0129538,0.0127529,0.0125413,0.012437,0.0124495,0.0118436,0.0118018,0.0118533,0.0117786,0.0112676,0.0113574,0.0114356,0.0112443,0.0110724,0.0113663,0.0108331,0.01085,0.0106594,0.0106056,0.0102667,0.0102377,0.0103971,0.0105035,0.010388,0.00998035,0.0100156,0.00979607,0.0098853,0.00969908,0.0098204,0.00948191,0.00963511,0.00966696,0.00950457,0.0092176,0.0094804,0.00926154,0.00899106};
        for (int idata=1; idata<251; idata++) {
            //fBWeight->SetBinContent(idata,ratioB[idata-1]);
            fBWeightNew->SetBinContent(idata,ratioBNew[idata-1]);
            fBWeightVar1->SetBinContent(idata,ratioBVar1[idata-1]);
            fBWeightVar2->SetBinContent(idata,ratioBVar2[idata-1]);
        }
    }else if (fCentralityMin==30 && fCentralityMax==50 && fApplyCentrality) {
        Double_t ratioB[250] = {0.613763,0.784834,0.810916,0.836203,0.830989,0.868295,0.883505,0.913221,0.934821,0.961457,0.985556,1.01031,1.03038,1.04399,1.04894,1.05124,1.04879,1.0374,1.02127,1.00945,0.982878,0.956372,0.934072,0.906477,0.887155,0.864324,0.842351,0.811834,0.784108,0.757401,0.730328,0.709669,0.681074,0.661868,0.631406,0.60966,0.590193,0.567588,0.54466,0.531031,0.516403,0.489309,0.472915,0.463391,0.447918,0.434322,0.422943,0.406763,0.397189,0.386448,0.376006,0.367809,0.360888,0.348902,0.337959,0.332709,0.326533,0.309368,0.305947,0.29407,0.29256,0.281315,0.277646,0.268647,0.25997,0.252176,0.245948,0.237849,0.230473,0.226727,0.219848,0.215782,0.207203,0.201134,0.196046,0.190369,0.188318,0.181247,0.17419,0.168473,0.162579,0.157984,0.15497,0.150132,0.144187,0.143221,0.139225,0.134864,0.129645,0.129193,0.123249,0.12106,0.119186,0.115747,0.113364,0.108896,0.107027,0.103195,0.101673,0.0990433,0.09838,0.0941805,0.0918415,0.0901931,0.0891644,0.0870442,0.0855439,0.0830778,0.0831012,0.0818048,0.0817035,0.0784107,0.0771535,0.0767316,0.0750625,0.0731386,0.0727318,0.0727746,0.0723,0.0707661,0.0694893,0.067761,0.0659144,0.065337,0.063405,0.0611122,0.062151,0.0613304,0.0606145,0.0604442,0.0577981,0.057869,0.0572931,0.0565873,0.0549912,0.0541402,0.0534318,0.0533835,0.0517918,0.0516942,0.051115,0.0508591,0.0495408,0.0494513,0.0499725,0.0483911,0.0479715,0.0464031,0.0463092,0.0458376,0.0445392,0.0433594,0.0431873,0.0424705,0.0415192,0.04182,0.0414201,0.0406436,0.041191,0.0404933,0.0395873,0.0392351,0.0372448,0.0381429,0.035662,0.0358194,0.0363743,0.0352793,0.0348819,0.0344932,0.0342752,0.0337997,0.0320732,0.0319788,0.031834,0.0309246,0.0303152,0.0300479,0.0296263,0.0287396,0.0286671,0.0280178,0.0277545,0.0288423,0.027348,0.0273515,0.0254216,0.0258547,0.0251372,0.0254501,0.0239668,0.0245532,0.0241397,0.0232339,0.0232454,0.0222495,0.0223577,0.021668,0.021833,0.0208466,0.0208539,0.0211989,0.0208092,0.0197581,0.0198031,0.0191019,0.0196992,0.0190488,0.019162,0.0189245,0.0187646,0.0176105,0.0178578,0.0170691,0.0170069,0.0161859,0.0169536,0.0156688,0.0155009,0.0164949,0.0150437,0.0150612,0.0148613,0.0148058,0.0142766,0.0148055,0.0145336,0.0150874,0.0144156,0.0138823,0.0139659,0.0130704,0.0133007,0.0138011,0.0131655,0.0131122,0.0133701,0.0131427,0.0125634,0.012817,0.0123759,0.0123669,0.0120754,0.0116675,0.0115908,0.0114222,0.0114497,0.0113251,0.012022,0.011139};
        Double_t ratioBNew[250] = {0.629879,0.788083,0.813433,0.837649,0.83473,0.871466,0.88658,0.917419,0.935438,0.961492,0.987841,1.01214,1.03055,1.04675,1.05152,1.05398,1.05146,1.03976,1.02127,1.01217,0.985053,0.956285,0.936065,0.908407,0.889645,0.863525,0.844254,0.814193,0.788496,0.759985,0.73241,0.710624,0.68301,0.663172,0.632443,0.610222,0.59217,0.569207,0.545748,0.531126,0.516459,0.489473,0.472985,0.465141,0.448287,0.435402,0.42446,0.408086,0.397874,0.386006,0.376618,0.368262,0.36165,0.349324,0.339087,0.332359,0.328113,0.309907,0.305498,0.294595,0.294602,0.282511,0.27876,0.270066,0.260252,0.253447,0.247503,0.238828,0.230165,0.226219,0.219926,0.217022,0.208488,0.201821,0.196422,0.190846,0.1887,0.181658,0.174466,0.169161,0.162862,0.158847,0.154782,0.150732,0.1443,0.143006,0.139508,0.135351,0.129897,0.129112,0.123742,0.121913,0.11981,0.116241,0.113353,0.108844,0.107325,0.10318,0.102095,0.099406,0.0986112,0.0947324,0.0920699,0.0910207,0.0894683,0.087799,0.0859247,0.0831379,0.083139,0.0824289,0.0816908,0.07866,0.0770417,0.0769627,0.0754568,0.0732459,0.0730796,0.0729844,0.0725179,0.0702804,0.0694963,0.0676671,0.0658698,0.0653475,0.0636846,0.0609512,0.0621307,0.061395,0.0607265,0.0608395,0.0579561,0.0583977,0.0574282,0.0567415,0.0553473,0.0545031,0.0534561,0.0533514,0.0518442,0.0519239,0.0513135,0.0509472,0.0498151,0.049784,0.0502312,0.0484379,0.0479072,0.0464586,0.046161,0.0456644,0.0446926,0.0436307,0.0437529,0.0425312,0.0417728,0.0419363,0.0415331,0.0405329,0.0409199,0.0405974,0.0398935,0.0395114,0.0375411,0.0380874,0.0360252,0.0357682,0.0366069,0.0356826,0.0350621,0.0344909,0.0342995,0.0336487,0.0321502,0.0318939,0.0315562,0.0312153,0.0307212,0.0299212,0.0297284,0.0287068,0.0285927,0.0281519,0.0276306,0.0287455,0.027473,0.0273274,0.0257916,0.0256825,0.0252547,0.0255688,0.0242197,0.0247519,0.0242898,0.0233756,0.0233464,0.0223576,0.0222939,0.0216614,0.0219721,0.0209188,0.0210538,0.0211145,0.0208338,0.0198707,0.0201319,0.0192302,0.0196697,0.0191208,0.0193517,0.0191076,0.0189315,0.0179213,0.0181558,0.0172182,0.0169685,0.0161772,0.0170894,0.0156943,0.0156453,0.0167149,0.015151,0.0151612,0.0150915,0.0149835,0.0143214,0.0148453,0.0144279,0.01508,0.0145104,0.0138127,0.0140673,0.0131155,0.0132893,0.0139538,0.0131358,0.0130515,0.0133015,0.0131268,0.0125072,0.0127295,0.0124395,0.0124556,0.0118749,0.0117609,0.0116434,0.0114951,0.0116084,0.011298,0.0121,0.0111589};
        Double_t ratioBVar1[250] = {0.573077,0.717766,0.741751,0.764901,0.763463,0.798537,0.814112,0.84448,0.863463,0.890334,0.918046,0.944505,0.966183,0.98657,0.997,1.00608,1.01129,1.00854,1,1.00155,0.986125,0.969671,0.962596,0.948575,0.944539,0.933351,0.930149,0.915439,0.905731,0.892737,0.880535,0.874948,0.861608,0.857325,0.837869,0.828285,0.82316,0.809796,0.793977,0.789405,0.783337,0.756718,0.744393,0.744272,0.728337,0.717369,0.708315,0.688908,0.67871,0.664664,0.653959,0.644245,0.636892,0.618813,0.603806,0.594538,0.589303,0.558567,0.552316,0.534034,0.535292,0.514361,0.508419,0.493305,0.475994,0.464062,0.453608,0.438059,0.422455,0.415447,0.404079,0.398897,0.383331,0.371165,0.361304,0.351096,0.347183,0.334247,0.321023,0.311261,0.299663,0.292261,0.284764,0.277289,0.265432,0.263023,0.256558,0.248881,0.238819,0.237341,0.227436,0.224038,0.220139,0.213545,0.208205,0.199887,0.197062,0.189418,0.187392,0.182424,0.180933,0.173784,0.168868,0.166912,0.164035,0.160944,0.157478,0.152342,0.152315,0.150985,0.149604,0.144026,0.141036,0.140865,0.138082,0.13401,0.13368,0.13348,0.132602,0.128485,0.127027,0.12366,0.120352,0.119375,0.116315,0.111301,0.113433,0.112068,0.110827,0.111011,0.10573,0.106515,0.104726,0.103454,0.100893,0.0993345,0.0974076,0.097198,0.094434,0.0945608,0.0934312,0.0927463,0.0906679,0.0905938,0.0913901,0.0881104,0.0871283,0.0844774,0.0839202,0.0830013,0.0812194,0.0792744,0.079481,0.0772469,0.0758548,0.0761371,0.0753906,0.0735609,0.0742491,0.0736497,0.0723589,0.071652,0.0680658,0.069043,0.0652924,0.0648141,0.0663211,0.0646342,0.063498,0.0624515,0.0620931,0.0609033,0.0581799,0.057705,0.0570831,0.0564556,0.0555514,0.0540944,0.0537356,0.051879,0.051663,0.0508568,0.0499054,0.0519092,0.0496018,0.0493295,0.0465484,0.0463426,0.045562,0.0461198,0.043678,0.0446292,0.0437878,0.0421317,0.0420709,0.0402814,0.0401591,0.0390122,0.0395643,0.0376604,0.0378963,0.0379983,0.0374861,0.0357463,0.0362092,0.0345808,0.0353644,0.034371,0.0347796,0.0343342,0.0340113,0.0321904,0.0326054,0.0309157,0.0304616,0.0290355,0.0306669,0.0281581,0.0280648,0.0299779,0.0271678,0.0271809,0.0270508,0.0268521,0.0256608,0.0265943,0.0258418,0.0270047,0.0259797,0.0247258,0.0251768,0.0234689,0.0237754,0.0249596,0.0234919,0.0233367,0.0237791,0.0234624,0.0223507,0.0227437,0.0222214,0.0222459,0.0212047,0.0209973,0.0207835,0.0205149,0.0207132,0.0201556,0.0215822,0.0198998};
        Double_t ratioBVar2[250] = {0.360924,0.463659,0.490858,0.517781,0.527756,0.562621,0.583404,0.614127,0.635702,0.661931,0.687458,0.71048,0.728122,0.742842,0.748039,0.750202,0.747529,0.737199,0.721136,0.711003,0.687762,0.663216,0.644631,0.621141,0.604108,0.582581,0.566292,0.543469,0.524326,0.504077,0.485197,0.470848,0.453267,0.441402,0.42274,0.41011,0.400569,0.387892,0.374939,0.368077,0.361179,0.345513,0.337034,0.33457,0.325438,0.318942,0.313639,0.304058,0.298802,0.292062,0.286966,0.282445,0.279074,0.271095,0.264535,0.260545,0.258367,0.245036,0.242464,0.234625,0.235382,0.226387,0.223988,0.217547,0.210128,0.205075,0.200668,0.193996,0.187287,0.18438,0.179529,0.17742,0.170682,0.165445,0.161225,0.156841,0.155261,0.149639,0.143875,0.139651,0.134594,0.131412,0.128179,0.12495,0.119736,0.118777,0.115983,0.112633,0.108196,0.107642,0.103261,0.101827,0.100162,0.097266,0.094935,0.0912398,0.0900464,0.0866456,0.0858103,0.0836242,0.0830286,0.0798326,0.0776566,0.0768387,0.0755938,0.0742479,0.0727258,0.0704279,0.0704896,0.0699478,0.069381,0.0668643,0.0655448,0.0655336,0.0643063,0.0624754,0.0623867,0.0623584,0.0620125,0.0601501,0.0595294,0.0580116,0.0565186,0.0561177,0.0547358,0.0524305,0.0534901,0.0529011,0.052369,0.0525104,0.0500636,0.0504872,0.0496905,0.0491372,0.0479698,0.0472773,0.0464077,0.0463552,0.045083,0.0451896,0.0446953,0.0444128,0.0434617,0.0434703,0.043897,0.0423645,0.0419348,0.0407001,0.0404724,0.0400697,0.0392491,0.0383478,0.0384864,0.0374422,0.0368044,0.0369785,0.0366527,0.0357989,0.03617,0.0359139,0.0353197,0.0350095,0.0332905,0.0338021,0.0319976,0.0317948,0.0325664,0.0317695,0.031242,0.0307575,0.0306113,0.0300544,0.0287388,0.0285323,0.0282526,0.0279695,0.0275486,0.0268524,0.0267005,0.0258032,0.025721,0.0253444,0.0248946,0.0259194,0.0247914,0.0246794,0.0233107,0.0232302,0.0228611,0.0231634,0.0219583,0.0224583,0.0220562,0.0212425,0.0212324,0.0203489,0.0203066,0.0197457,0.0200444,0.0190981,0.0192363,0.0193066,0.0190645,0.0181971,0.0184505,0.0176375,0.0180544,0.017564,0.0177897,0.0175787,0.0174299,0.0165124,0.0167412,0.0158887,0.0156701,0.0149506,0.0158056,0.0145263,0.0144918,0.0154943,0.0140551,0.0140751,0.0140209,0.013931,0.0133255,0.0138232,0.0134447,0.0140628,0.0135417,0.0129002,0.0131477,0.0122673,0.0124391,0.0130708,0.0123136,0.0122436,0.0124873,0.0123325,0.011759,0.0119768,0.0117126,0.0117363,0.0111974,0.0110981,0.0109952,0.0108631,0.0109782,0.0106925,0.0114598,0.0105762};
        /*Double_t ratioB[250] = {0.528468,0.676528,0.69985,0.722598,0.719079,0.75247,0.766871,0.794031,0.81433,0.839234,0.862171,0.885957,0.905921,0.920518,0.92776,0.932943,0.93418,0.927708,0.917184,0.910725,0.891085,0.871534,0.855825,0.835218,0.822144,0.805686,0.789807,0.765574,0.74352,0.721928,0.69942,0.682462,0.657241,0.640435,0.612102,0.591613,0.572795,0.55046,0.527428,0.513092,0.497552,0.469891,0.452484,0.441642,0.425178,0.410609,0.398269,0.381579,0.371263,0.360024,0.349237,0.340697,0.333486,0.32174,0.311097,0.305812,0.299775,0.283748,0.280412,0.269394,0.267931,0.257601,0.254249,0.246051,0.238174,0.231127,0.225531,0.218232,0.211603,0.208313,0.202149,0.198573,0.190843,0.18542,0.180898,0.175828,0.174104,0.167734,0.161368,0.156232,0.150923,0.146812,0.144165,0.139813,0.134422,0.133665,0.130076,0.126139,0.12139,0.121099,0.115653,0.113723,0.112086,0.108971,0.106844,0.102746,0.101093,0.0975802,0.0962457,0.0938589,0.0933322,0.0894457,0.0873193,0.0858455,0.0849587,0.0830286,0.081686,0.0794172,0.0795254,0.0783695,0.078357,0.0752801,0.0741529,0.0738267,0.0722983,0.0705208,0.0702036,0.0703201,0.069936,0.0685253,0.0673606,0.065755,0.064031,0.0635373,0.0617238,0.0595547,0.060631,0.0598935,0.0592567,0.0591523,0.0566221,0.056751,0.0562451,0.0556102,0.054098,0.0533164,0.0526735,0.0526806,0.0511629,0.0511194,0.0505989,0.0503977,0.0491419,0.0491037,0.0496723,0.0481498,0.0477813,0.0462665,0.0462201,0.0457961,0.0445443,0.0434086,0.0432802,0.0426051,0.041693,0.0420377,0.0416778,0.0409378,0.041531,0.0408687,0.0399945,0.0396785,0.0377035,0.0386513,0.0361735,0.0363695,0.0369697,0.0358926,0.0355235,0.0351626,0.0349751,0.034524,0.0327929,0.0327288,0.0326127,0.0317123,0.031118,0.0308738,0.0304706,0.0295875,0.0295417,0.0289009,0.0286572,0.0298095,0.0282926,0.0283237,0.0263508,0.0268257,0.0261065,0.0264569,0.0249391,0.0255739,0.0251674,0.0242464,0.0242816,0.0232636,0.0233992,0.0226991,0.0228937,0.0218802,0.0219088,0.0222925,0.0219034,0.0208167,0.020884,0.0201636,0.0208137,0.0201455,0.0202843,0.0200517,0.019901,0.0186946,0.0189748,0.0181538,0.0181046,0.0172467,0.0180816,0.0167268,0.016563,0.0176415,0.0161043,0.016138,0.0159385,0.0158937,0.0153398,0.0159227,0.0156447,0.0162558,0.0155462,0.0149848,0.0150888,0.0141342,0.0143964,0.0149517,0.0142761,0.0142313,0.0145244,0.0142902,0.0136728,0.0139614,0.0134931,0.0134955,0.0131893,0.0127552,0.0126828,0.0125095,0.0125509,0.0124254,0.0132018,0.0122431};
        Double_t ratioBNew[250] = {0.528468,0.676528,0.69985,0.722598,0.719079,0.75247,0.766871,0.794031,0.81433,0.839234,0.862171,0.885957,0.905921,0.920518,0.92776,0.932943,0.93418,0.927708,0.917184,0.910725,0.891085,0.871534,0.855825,0.835218,0.822144,0.805686,0.789807,0.765574,0.74352,0.721928,0.69942,0.682462,0.657241,0.640435,0.612102,0.591613,0.572795,0.55046,0.527428,0.513092,0.497552,0.469891,0.452484,0.441642,0.425178,0.410609,0.398269,0.381579,0.371263,0.360024,0.349237,0.340697,0.333486,0.32174,0.311097,0.305812,0.299775,0.283748,0.280412,0.269394,0.267931,0.257601,0.254249,0.246051,0.238174,0.231127,0.225531,0.218232,0.211603,0.208313,0.202149,0.198573,0.190843,0.18542,0.180898,0.175828,0.174104,0.167734,0.161368,0.156232,0.150923,0.146812,0.144165,0.139813,0.134422,0.133665,0.130076,0.126139,0.12139,0.121099,0.115653,0.113723,0.112086,0.108971,0.106844,0.102746,0.101093,0.0975802,0.0962457,0.0938589,0.0933322,0.0894457,0.0873193,0.0858455,0.0849587,0.0830286,0.081686,0.0794172,0.0795254,0.0783695,0.078357,0.0752801,0.0741529,0.0738267,0.0722983,0.0705208,0.0702036,0.0703201,0.069936,0.0685253,0.0673606,0.065755,0.064031,0.0635373,0.0617238,0.0595547,0.060631,0.0598935,0.0592567,0.0591523,0.0566221,0.056751,0.0562451,0.0556102,0.054098,0.0533164,0.0526735,0.0526806,0.0511629,0.0511194,0.0505989,0.0503977,0.0491419,0.0491037,0.0496723,0.0481498,0.0477813,0.0462665,0.0462201,0.0457961,0.0445443,0.0434086,0.0432802,0.0426051,0.041693,0.0420377,0.0416778,0.0409378,0.041531,0.0408687,0.0399945,0.0396785,0.0377035,0.0386513,0.0361735,0.0363695,0.0369697,0.0358926,0.0355235,0.0351626,0.0349751,0.034524,0.0327929,0.0327288,0.0326127,0.0317123,0.031118,0.0308738,0.0304706,0.0295875,0.0295417,0.0289009,0.0286572,0.0298095,0.0282926,0.0283237,0.0263508,0.0268257,0.0261065,0.0264569,0.0249391,0.0255739,0.0251674,0.0242464,0.0242816,0.0232636,0.0233992,0.0226991,0.0228937,0.0218802,0.0219088,0.0222925,0.0219034,0.0208167,0.020884,0.0201636,0.0208137,0.0201455,0.0202843,0.0200517,0.019901,0.0186946,0.0189748,0.0181538,0.0181046,0.0172467,0.0180816,0.0167268,0.016563,0.0176415,0.0161043,0.016138,0.0159385,0.0158937,0.0153398,0.0159227,0.0156447,0.0162558,0.0155462,0.0149848,0.0150888,0.0141342,0.0143964,0.0149517,0.0142761,0.0142313,0.0145244,0.0142902,0.0136728,0.0139614,0.0134931,0.0134955,0.0131893,0.0127552,0.0126828,0.0125095,0.0125509,0.0124254,0.0132018,0.0122431};
        Double_t ratioBVar1[250] = {0.543441,0.695667,0.719653,0.743089,0.73956,0.774051,0.78908,0.817324,0.838612,0.864768,0.889046,0.914378,0.93597,0.952245,0.961158,0.968204,0.971448,0.966978,0.958592,0.954794,0.937517,0.920647,0.908186,0.890889,0.882019,0.86995,0.858929,0.839181,0.822106,0.805814,0.788725,0.778117,0.758203,0.748037,0.724299,0.709568,0.696604,0.678976,0.659911,0.651178,0.640401,0.613179,0.598384,0.591556,0.576458,0.563099,0.552026,0.534127,0.524402,0.512725,0.501066,0.492075,0.484518,0.469903,0.456447,0.450488,0.443121,0.420671,0.416769,0.401237,0.399756,0.384892,0.380318,0.368382,0.356827,0.346432,0.338145,0.327248,0.317312,0.312347,0.303043,0.297595,0.285905,0.27766,0.270755,0.263024,0.260292,0.250612,0.240941,0.233114,0.225032,0.218743,0.214638,0.208,0.199823,0.198542,0.193057,0.187062,0.179873,0.179295,0.17109,0.168097,0.165539,0.160804,0.157534,0.151365,0.148804,0.143514,0.141432,0.137809,0.13692,0.131108,0.127884,0.12562,0.124218,0.121294,0.119233,0.115824,0.115885,0.114106,0.113992,0.109425,0.107697,0.107134,0.104829,0.102167,0.101624,0.101708,0.10107,0.0989493,0.0971877,0.0947934,0.0922324,0.0914466,0.088764,0.0855749,0.0870506,0.0859218,0.0849395,0.0847212,0.0810318,0.0811507,0.0803625,0.0793915,0.0771707,0.0759947,0.0750183,0.0749686,0.0727507,0.072631,0.0718343,0.0714919,0.0696553,0.069546,0.0702958,0.0680875,0.0675132,0.0653215,0.0652048,0.0645562,0.0627425,0.061095,0.0608669,0.0598709,0.0585437,0.0589818,0.0584316,0.0573498,0.0581359,0.0571648,0.055899,0.0554147,0.0526161,0.0538975,0.0504038,0.0506381,0.0514347,0.0498981,0.0493476,0.0488092,0.0485122,0.0478504,0.0454168,0.0452938,0.0450992,0.0438211,0.0429676,0.0425987,0.0420108,0.0407629,0.0406695,0.0397577,0.0393932,0.0409468,0.0388344,0.0388484,0.0361157,0.0367395,0.0357283,0.0361813,0.0340805,0.0349224,0.0343423,0.0330614,0.0330853,0.0316751,0.0318366,0.0308616,0.0311037,0.0297053,0.0297226,0.0302213,0.0296726,0.0281802,0.0282509,0.0272569,0.0281156,0.0271936,0.0273614,0.0270285,0.0268063,0.0251633,0.0255225,0.0244009,0.0243176,0.0231489,0.0242524,0.0224195,0.0221844,0.0236123,0.0215398,0.0215699,0.0212883,0.0212137,0.0204601,0.0212229,0.0208379,0.0216368,0.0206781,0.0199176,0.020042,0.0187612,0.0190961,0.019819,0.0189105,0.0188383,0.0192131,0.0188905,0.018062,0.0184308,0.0178004,0.0177915,0.0173761,0.0167928,0.0166862,0.0164472,0.0164905,0.0163147,0.0173225,0.0160538};
        Double_t ratioBVar2[250] = {0.513494,0.657388,0.680046,0.702107,0.698598,0.73089,0.744663,0.770738,0.790048,0.813701,0.835295,0.857536,0.875872,0.888791,0.894362,0.897681,0.896912,0.888437,0.875776,0.866656,0.844654,0.822421,0.803464,0.779547,0.762269,0.741421,0.720684,0.691967,0.664935,0.638041,0.610115,0.586807,0.556278,0.532833,0.499906,0.473658,0.448987,0.421944,0.394945,0.375005,0.354702,0.326604,0.306584,0.291727,0.273898,0.258118,0.244512,0.22903,0.218123,0.207324,0.197408,0.189319,0.182453,0.173577,0.165747,0.161136,0.156429,0.146825,0.144055,0.13755,0.136105,0.130309,0.12818,0.12372,0.119522,0.115823,0.112918,0.109216,0.105894,0.104279,0.101255,0.099551,0.095781,0.0931801,0.0910409,0.0886326,0.0879164,0.0848564,0.0817939,0.0793508,0.0768141,0.0748818,0.0736922,0.071627,0.0690201,0.0687883,0.0670956,0.0652155,0.0629067,0.0629027,0.060215,0.0593496,0.0586328,0.0571375,0.0561542,0.054127,0.0533809,0.0516468,0.0510594,0.0499091,0.0497441,0.0477829,0.0467545,0.0460709,0.0456993,0.0447628,0.0441389,0.0430099,0.0431655,0.0426333,0.0427216,0.0411352,0.0406088,0.0405192,0.0397673,0.0388742,0.0387834,0.0389317,0.0388025,0.0381012,0.0375335,0.0367166,0.0358295,0.0356281,0.0346837,0.0335346,0.0342115,0.0338651,0.033574,0.0335835,0.0322125,0.0323513,0.0321277,0.0318289,0.0310254,0.030638,0.0303287,0.0303927,0.0295751,0.0296079,0.0293635,0.0293035,0.0286286,0.0286613,0.0290488,0.0282121,0.0280493,0.0272114,0.0272353,0.0270361,0.0263462,0.0257221,0.0256935,0.0253394,0.0248424,0.0250935,0.024924,0.0245258,0.0249261,0.0245726,0.0240901,0.0239423,0.0227909,0.0234051,0.0219433,0.0221008,0.0225047,0.021887,0.0216995,0.021516,0.021438,0.0211977,0.0201691,0.0201638,0.0201262,0.0196034,0.0192683,0.019149,0.0189303,0.0184121,0.018414,0.0180441,0.0179212,0.0186722,0.0177508,0.017799,0.0165858,0.0169118,0.0164847,0.0167326,0.0157976,0.0162253,0.0159925,0.0154314,0.015478,0.0148521,0.0149618,0.0145365,0.0146837,0.0140551,0.0140949,0.0143636,0.0141343,0.0134533,0.013517,0.0130702,0.0135118,0.0130974,0.0132072,0.013075,0.0129958,0.0122258,0.0124272,0.0119067,0.0118916,0.0113445,0.0119107,0.0110341,0.0109416,0.0116706,0.0106688,0.0107062,0.0105887,0.0105737,0.0102194,0.0106225,0.0104515,0.0108747,0.0104143,0.010052,0.0101356,0.00950728,0.00969676,0.0100843,0.00964163,0.00962429,0.00983562,0.00968996,0.00928355,0.00949209,0.00918577,0.00919945,0.0090025,0.00871755,0.00867932,0.00857183,0.00861128,0.00853617,0.00908116,0.00843244};*/
        for (int idata=1; idata<251; idata++) {
            //fBWeight->SetBinContent(idata,ratioB[idata-1]);
            fBWeightNew->SetBinContent(idata,ratioBNew[idata-1]);
            fBWeightVar1->SetBinContent(idata,ratioBVar1[idata-1]);
            fBWeightVar2->SetBinContent(idata,ratioBVar2[idata-1]);
        }
    }else{
        for (int idata=1; idata<251; idata++) {
            //fBWeight->SetBinContent(idata,ratioB[idata-1]);
            fBWeightNew->SetBinContent(idata,1);
            fBWeightVar1->SetBinContent(idata,1);
            fBWeightVar2->SetBinContent(idata,1);
        }
    }

    //fOutputList->Add(fBWeight);
    fOutputList->Add(fBWeightNew);
    fOutputList->Add(fBWeightVar1);
    fOutputList->Add(fBWeightVar2);
    
    if (fFlagFillMCHistos) {
        //fPi0DCA = new TH2F("fPi0DCA","Pi0 DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        //fOutputList->Add(fPi0DCA);
    
        //fEtaDCA = new TH2F("fEtaDCA","Eta DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        //fOutputList->Add(fEtaDCA);
        
        //Decay lengths and weights
        /*Double_t tau_B0Py6 = 468;
        Double_t tau_BPlusPy6 = 462;
        Double_t tau_BsPy6 = 483;*/
        Double_t tau_B0Py6 = 444;
        Double_t tau_BPlusPy6 = 446;
        Double_t tau_BsPy6 = 434;
        Double_t tau_D0Py6 = 124.4;
        Double_t tau_DPlusPy6 = 317;
        Double_t tau_DsPy6 = 140;
        
        /*Double_t tau_B0Py8 = 458.7;
        Double_t tau_BPlusPy8 = 491.1;
        Double_t tau_BsPy8 = 439;*/
        Double_t tau_B0Py8 = 454;
        Double_t tau_BPlusPy8 = 471;
        Double_t tau_BsPy8 = 362;
        Double_t tau_D0Py8 = 122.9;
        Double_t tau_DPlusPy8 = 311.8;
        Double_t tau_DsPy8 = 149.9;
        
        if (fUseTauWeight) {
            fBPlusTauWeight = new TF1("fBPlusTauWeight","exp(-x/[0])/exp(-x/[1])",0,10000);
            fBPlusTauWeight->FixParameter(0,tau_BPlusPy8);
            fBPlusTauWeight->FixParameter(1,tau_BPlusPy6);
            fOutputList->Add(fBPlusTauWeight);
            
            fB0TauWeight = new TF1("fB0TauWeight","exp(-x/[0])/exp(-x/[1])",0,10000);
            fB0TauWeight->FixParameter(0,tau_B0Py8);
            fB0TauWeight->FixParameter(1,tau_B0Py6);
            fOutputList->Add(fB0TauWeight);
            
            fBsTauWeight = new TF1("fBsTauWeight","exp(-x/[0])/exp(-x/[1])",0,10000);
            fBsTauWeight->FixParameter(0,tau_BsPy8);
            fBsTauWeight->FixParameter(1,tau_BsPy6);
            fOutputList->Add(fBsTauWeight);
            
            fDPlusTauWeight = new TF1("fDPlusTauWeight","exp(-x/[0])/exp(-x/[1])",0,10000);
            fDPlusTauWeight->FixParameter(0,tau_DPlusPy8);
            fDPlusTauWeight->FixParameter(1,tau_DPlusPy6);
            fOutputList->Add(fDPlusTauWeight);
            
            fD0TauWeight = new TF1("fD0TauWeight","exp(-x/[0])/exp(-x/[1])",0,10000);
            fD0TauWeight->FixParameter(0,tau_D0Py8);
            fD0TauWeight->FixParameter(1,tau_D0Py6);
            fOutputList->Add(fD0TauWeight);
            
            fDsTauWeight = new TF1("fDsTauWeight","exp(-x/[0])/exp(-x/[1])",0,10000);
            fDsTauWeight->FixParameter(0,tau_DsPy8);
            fDsTauWeight->FixParameter(1,tau_DsPy6);
            fOutputList->Add(fDsTauWeight);
        }else{
            fBPlusTauWeight = new TF1("fBPlusTauWeight","1",0,10000);
            fOutputList->Add(fBPlusTauWeight);
            
            fB0TauWeight = new TF1("fB0TauWeight","1",0,10000);
            fOutputList->Add(fB0TauWeight);
            
            fBsTauWeight = new TF1("fBsTauWeight","1",0,10000);
            fOutputList->Add(fBsTauWeight);
            
            fDPlusTauWeight = new TF1("fDPlusTauWeight","1",0,10000);
            fOutputList->Add(fDPlusTauWeight);
            
            fD0TauWeight = new TF1("fD0TauWeight","1",0,10000);
            fOutputList->Add(fD0TauWeight);
            
            fDsTauWeight = new TF1("fDsTauWeight","1",0,10000);
            fOutputList->Add(fDsTauWeight);
        }
        
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
    
        fPhotonHijingTagDCA = new TH2F("fPhotonHijingTagDCA","Tagged Hijing Photon DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fPhotonHijingTagDCA->Sumw2();
        fOutputList->Add(fPhotonHijingTagDCA);
        fEnhPhotonTagDCA = new TH2F("fEnhPhotonTagDCA","Tagged Enh Photon DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fEnhPhotonTagDCA->Sumw2();
        fOutputList->Add(fEnhPhotonTagDCA);
        
        fComboNumWeight = new TH3F("fComboNumWeight","Eff Num with Weight; p_{T}(GeV/c); prod. radius; mother pid", 60,0,30.,250,0.,50.,3,2.5,5.5);
        fComboNumWeight->Sumw2();
        fOutputList->Add(fComboNumWeight);
        fComboNumNoWeight = new TH3F("fComboNumNoWeight","Eff Num Without Weight; p_{T}(GeV/c); prod. radius; mother pid", 60,0,30.,250,0.,50.,3,2.5,5.5);
        fComboNumNoWeight->Sumw2();
        fOutputList->Add(fComboNumNoWeight);
        fComboDenomWeight = new TH3F("fComboDenomWeight","Eff Denom with Weight; p_{T}(GeV/c); prod. radius; mother pid", 60,0,30.,250,0.,50.,3,2.5,5.5);
        fComboDenomWeight->Sumw2();
        fOutputList->Add(fComboDenomWeight);
        fComboDenomNoWeight = new TH3F("fComboDenomNoWeight","Eff Denom without Weight; p_{T}(GeV/c); prod. radius; mother pid", 60,0,30.,250,0.,50.,3,2.5,5.5);
        fComboDenomNoWeight->Sumw2();
        fOutputList->Add(fComboDenomNoWeight);
    
        //fDMesonPDG = new TH1F("fDMesonPDG","D Meson PDG values; abs(PDG); counts;", 51,399.5,450.5);
        //DMesonPDG->Sumw2();
        //fOutputList->Add(fDMesonPDG);
    
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
    
        fD0MesonPtWeight = new TH1F("fD0MesonPtWeight","D0 Spectrum w/weight; p_{T}(GeV/c); counts;",100,0,50.);
        fD0MesonPtWeight->Sumw2();
        fOutputList->Add(fD0MesonPtWeight);
        
        fLambdaCPtWeight = new TH1F("fLambdaCPtWeight","Lc Spectrum w/weight; p_{T}(GeV/c); counts;",100,0,50.);
        fLambdaCPtWeight->Sumw2();
        fOutputList->Add(fLambdaCPtWeight);
        
        fDPlusMesonPtWeight = new TH1F("fDPlusMesonPtWeight","D+ Spectrum w/D+ weight; p_{T}(GeV/c); counts;",100,0,50.);
        fDPlusMesonPtWeight->Sumw2();
        fOutputList->Add(fDPlusMesonPtWeight);
        
        fDsMesonPtWeight = new TH1F("fDsMesonPtWeight","Ds Spectrum w/Ds weight; p_{T}(GeV/c); counts;",100,0,50.);
        fDsMesonPtWeight->Sumw2();
        fOutputList->Add(fDsMesonPtWeight);
        
        fEtaCPt = new TH1F("fEtaCPt","Eta_c Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
        fEtaCPt->Sumw2();
        fOutputList->Add(fEtaCPt);
    
        fCBaryonPt = new TH1F("fCBaryonPt","Charm Baryon Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
        fCBaryonPt->Sumw2();
        fOutputList->Add(fCBaryonPt);
    
        fBMesonPt = new TH1F("fBMesonPt","B Meson Spectrum; p_{T}(GeV/c); counts;",250,0,50.);
        fBMesonPt->Sumw2();
        fOutputList->Add(fBMesonPt);
    
        /*fBMesonPtATLAS = new TH1F("fBMesonPtATLAS","ATLAS B Meson Spectrum; p_{T}(GeV/c); counts;",240,0,120.);
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
        fOutputList->Add(fBPlusPtLHCb);*/
        
        fBBaryonPt = new TH1F("fBBaryonPt","Beauty Baryon Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
        fBBaryonPt->Sumw2();
        fOutputList->Add(fBBaryonPt);
    
        fBMesonElecPt = new TH1F("fBMesonElecPt","Beauty Meson->e Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
        fBMesonElecPt->Sumw2();
        fOutputList->Add(fBMesonElecPt);
        
        fBBaryonElecPt = new TH1F("fBBaryonElecPt","Beauty Baryon->e Spectrum; p_{T}(GeV/c); counts;",100,0,50.);
        fBBaryonElecPt->Sumw2();
        fOutputList->Add(fBBaryonElecPt);
        
        /*fPromptD0DCAWeight = new TH2F("fPromptD0DCAWeight","Prompt D0 DCA with Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fPromptD0DCAWeight);
    
        fD0FromDStarDCAWeight = new TH2F("fD0FromDStarDCAWeight","Prompt D0 DCA with Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fD0FromDStarDCAWeight);
    
        fPromptD0DCANoWeight = new TH2F("fPromptD0DCANoWeight","Prompt D0 DCA without Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fPromptD0DCANoWeight);
    
        fD0FromDStarDCANoWeight = new TH2F("fD0FromDStarDCANoWeight","Prompt D0 DCA without Weight; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., nDCAbins,-0.2,0.2);
        fOutputList->Add(fD0FromDStarDCANoWeight);*/
    
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
    
        Int_t binTemp[5] = {60,nDCAbins,22,20,200}; //pT, DCA, Mom PID, momGamma, momTime
        Double_t xminTemp[5] = {0.,-0.2,0.5,1.,0.};
        Double_t xmaxTemp[5] = {30.,0.2,22.5,21.,2000.};
        fSprsTemplatesNoWeight = new THnSparseD("fSprsTemplatesNoWeight","Sparse for Templates, No weight applied;p_{T};DCA;MomPID;Mom #gamma;Mom ct;",5,binTemp,xminTemp,xmaxTemp);
        fOutputList->Add(fSprsTemplatesNoWeight);
        fSprsTemplatesNoWeight->Sumw2();
        fSprsTemplatesWeight = new THnSparseD("fSprsTemplatesWeight","Sparse for Templates, D meson weight applied;p_{T};DCA;MomPID;Mom #gamma;Mom ct;",5,binTemp,xminTemp,xmaxTemp);
        fOutputList->Add(fSprsTemplatesWeight);
        fSprsTemplatesWeight->Sumw2();
        fSprsTemplatesWeightVar1 = new THnSparseD("fSprsTemplatesWeightVar1","Sparse for Templates, D meson WeightVar1 applied;p_{T};DCA;MomPID;Mom #gamma;Mom ct;",5,binTemp,xminTemp,xmaxTemp);
        fOutputList->Add(fSprsTemplatesWeightVar1);
        fSprsTemplatesWeightVar1->Sumw2();
        fSprsTemplatesWeightVar2 = new THnSparseD("fSprsTemplatesWeightVar2","Sparse for Templates, D meson WeightVar2 applied;p_{T};DCA;MomPID;Mom #gamma;Mom ct;",5,binTemp,xminTemp,xmaxTemp);
        fOutputList->Add(fSprsTemplatesWeightVar2);
        fSprsTemplatesWeightVar2->Sumw2();
        
        Int_t binClos[4] = {60,nDCAbins,22,250}; //pT, DCA, Mom PID, prod. radius
        Double_t xminClos[4] = {0.,-0.2,0.5,0.};
        Double_t xmaxClos[4] = {30.,0.2,22.5,50.};
        fSprsClosureTest = new THnSparseD("fSprsClosureTest","Sparse for Closure Test;p_{T};DCA;MomPID;prod.R",4,binClos,xminClos,xmaxClos);
        fOutputList->Add(fSprsClosureTest);
        fSprsClosureTest->Sumw2();
        
        fSprsClosureTestWeight = new THnSparseD("fSprsClosureTestWeight","Sparse for Closure Test w/ pi0+eta weight;p_{T};DCA;MomPID;prod.R",4,binClos,xminClos,xmaxClos);
        fOutputList->Add(fSprsClosureTestWeight);
        fSprsClosureTestWeight->Sumw2();
        
        Int_t binULS[4] = {60,nDCAbins,3,250}; //pT, DCA, Mom PID, prod. radius
        Double_t xminULS[4] = {0.,-0.2,2.5,0.};
        Double_t xmaxULS[4] = {30.,0.2,5.5,50.};
        fSprsULSdca = new THnSparseD("fSprsULSdca","ULS Sparse for Closure Test;p_{T};DCA;MomPID;prod.R",4,binULS,xminULS,xmaxULS);
        fOutputList->Add(fSprsULSdca);
        fSprsULSdca->Sumw2();
        
        fSprsULSdcaWeight = new THnSparseD("fSprsULSdcaWeight","ULS Sparse for Closure Test w/weight;p_{T};DCA;MomPID;prod.R",4,binULS,xminULS,xmaxULS);
        fOutputList->Add(fSprsULSdcaWeight);
        fSprsULSdcaWeight->Sumw2();
        
        fSprsLSdca = new THnSparseD("fSprsLSdca","LS Sparse for Closure Test;p_{T};DCA;MomPID;prod.R",4,binULS,xminULS,xmaxULS);
        fOutputList->Add(fSprsLSdca);
        fSprsLSdca->Sumw2();
        
        fSprsLSdcaWeight = new THnSparseD("fSprsLSdcaWeight","LS Sparse for Closure Test w/weight;p_{T};DCA;MomPID;prod.R",4,binULS,xminULS,xmaxULS);
        fOutputList->Add(fSprsLSdcaWeight);
        fSprsLSdcaWeight->Sumw2();
    
        /*fDTemplateWeight = new TH2F("fDTemplateWeight","D Meson DCA template", 100,0,50., nDCAbins,-0.2,0.2);
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
        fOutputList->Add(fBTemplateWeightVar2);*/
        
        fAllElecStack = new TH1F("fAllElecStack","All Elec from Stack; p_{T}(GeV/c); counts;",100,0,50.);
        fAllElecStack->Sumw2();
        fOutputList->Add(fAllElecStack);
    
        fHFElecStack = new TH1F("fHFElecStack","HF Elec from Stack; p_{T}(GeV/c); counts;",100,0,50.);
        fHFElecStack->Sumw2();
        fOutputList->Add(fHFElecStack);
    
        fBElecStack = new TH1F("fBElecStack","B Elec from Stack; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecStack->Sumw2();
        fOutputList->Add(fBElecStack);
        
        fAllElecStackDiffPID = new TH1F("fAllElecStackDiffPID","All Elec from Stack, Shingo mother PID; p_{T}(GeV/c); counts;",100,0,50.);
        fAllElecStackDiffPID->Sumw2();
        fOutputList->Add(fAllElecStackDiffPID);
        
        fDElecStackDiffPID = new TH1F("fDElecStackDiffPID","D Elec from Stack, Shingo mother PID; p_{T}(GeV/c); counts;",100,0,50.);
        fDElecStackDiffPID->Sumw2();
        fOutputList->Add(fDElecStackDiffPID);
        
        fBElecStackDiffPID = new TH1F("fBElecStackDiffPID","B Elec from Stack, Shingo mother PID; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecStackDiffPID->Sumw2();
        fOutputList->Add(fBElecStackDiffPID);
    
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
        
        fElecAftLooseTrkCuts = new TH1F("fElecAftLooseTrkCuts","Elec after loose trk cuts; p_{T}(GeV/c); counts;",100,0,50.);
        fElecAftLooseTrkCuts->Sumw2();
        fOutputList->Add(fElecAftLooseTrkCuts);
        
        fHFElecAftLooseTrkCuts = new TH1F("fHFElecAftLooseTrkCuts","HF Elec after loose trk cuts; p_{T}(GeV/c); counts;",100,0,50.);
        fHFElecAftLooseTrkCuts->Sumw2();
        fOutputList->Add(fHFElecAftLooseTrkCuts);
        
        fBElecAftLooseTrkCuts = new TH1F("fBElecAftLooseTrkCuts","B Elec after loose trk cuts; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecAftLooseTrkCuts->Sumw2();
        fOutputList->Add(fBElecAftLooseTrkCuts);
        
        /*fElecAftLooseTrkCutsDiffPID = new TH1F("fElecAftLooseTrkCutsDiffPID","Elec after loose trk cuts, Shingo's mother PID; p_{T}(GeV/c); counts;",100,0,50.);
        fElecAftLooseTrkCutsDiffPID->Sumw2();
        fOutputList->Add(fElecAftLooseTrkCutsDiffPID);
        
         fDElecAftLooseTrkCutsDiffPID = new TH1F("fDElecAftLooseTrkCutsDiffPID","D Elec after loose trk cuts, Shingo's mother PID; p_{T}(GeV/c); counts;",100,0,50.);
        fDElecAftLooseTrkCutsDiffPID->Sumw2();
        fOutputList->Add(fDElecAftLooseTrkCutsDiffPID);
        
        fBElecAftLooseTrkCutsDiffPID = new TH1F("fBElecAftLooseTrkCutsDiffPID","B Elec after loose trk cuts, Shingo's mother PID; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecAftLooseTrkCutsDiffPID->Sumw2();
        fOutputList->Add(fBElecAftLooseTrkCutsDiffPID);*/
    
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
        
        fElecAftEoP = new TH1F("fElecAftEoP","Elec after EoP cut; p_{T}(GeV/c); counts;",100,0,50.);
        fElecAftEoP->Sumw2();
        fOutputList->Add(fElecAftEoP);
        
        fHFElecAftEoP = new TH1F("fHFElecAftEoP","HF Elec after EoP cut; p_{T}(GeV/c); counts;",100,0,50.);
        fHFElecAftEoP->Sumw2();
        fOutputList->Add(fHFElecAftEoP);
        
        fBElecAftEoP = new TH1F("fBElecAftEoP","B Elec after EoP cut; p_{T}(GeV/c); counts;",100,0,50.);
        fBElecAftEoP->Sumw2();
        fOutputList->Add(fBElecAftEoP);
    }
    
    if (fFlagFillSprs && fFlagFillMCHistos) {
        Int_t bins1[6]=  {/*280*/60,  160, 100, 100,  nDCAbins, 4}; // pT;nSigma;eop;M20/M02;DCA;MCTruth
        Double_t xmin1[6]={ /*2*/0,   -8,   0,   0, -0.2, -0.5};
        Double_t xmax1[6]={30,    8,   2,   1,  0.2, 3.5};
        fElectronSprs = new THnSparseD("Electron","Electron;pT;nSigma;eop;M20/M02;DCA;MCTruth;",6,bins1,xmin1,xmax1);
        fOutputList->Add(fElectronSprs);
    }
    if (fFlagFillSprs && !fFlagFillMCHistos) {
        Int_t bins1[5]=  {60,  160, 100, 100,  nDCAbins}; // pT;nSigma;eop;M20/M02;DCA
        Double_t xmin1[5]={ /*2*/0,   -8,   0,   0, -0.2};
        Double_t xmax1[5]={30,    8,   2,   1,  0.2};
        fElectronSprs = new THnSparseD("Electron","Electron;pT;nSigma;eop;M20/M02;DCA;",5,bins1,xmin1,xmax1);
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
    
    
    /*if (fFlagFillMCHistos && fFlagRunStackLoop) {
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
            
                if(TrackPDG == 11 && AODMCtrack->IsPhysicalPrimary() && TMath::Abs(AODMCtrack->Eta()) <= 0.6) {
                    // cout<<"TESTINGGGGGGGGGGGG"<<endl;
                    fAllElecStack->Fill(AODMCtrack->Pt());
                    eleinStack++;
                }
                //
                *//*if(TMath::Abs(AODMCtrack->Y()) < 2.25) {
                    if (TrackPDG>500 && TrackPDG<599) fBMesonPtATLAS->Fill(AODMCtrack->Pt());
                    if (TrackPDG == 521) fBPlusPtATLAS->Fill(AODMCtrack->Pt());
                }
                if(TMath::Abs(AODMCtrack->Y()) < 2.4) {
                    if (TrackPDG>500 && TrackPDG<599) fBMesonPtCMS->Fill(AODMCtrack->Pt());
                    if (TrackPDG == 521) fBPlusPtCMS->Fill(AODMCtrack->Pt());
                
                    if(TrackPDG==11 && ilabelM>0){
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
                    }
                }
                if(TMath::Abs(AODMCtrack->Y()) > 2.0 && TMath::Abs(AODMCtrack->Y()) < 4.5) {
                    if (TrackPDG>500 && TrackPDG<599) fBMesonPtLHCb->Fill(AODMCtrack->Pt());
                    if (TrackPDG == 521) fBPlusPtLHCb->Fill(AODMCtrack->Pt());
                }*/
                
                /*if(TMath::Abs(AODMCtrack->Eta()) > 0.6) continue;
                if (TrackPDG>500 && TrackPDG<599) {
                    fBMesonPt->Fill(AODMCtrack->Pt());
                }
                if (TrackPDG>5000 && TrackPDG<5999) {
                    fBBaryonPt->Fill(AODMCtrack->Pt());
                }
            
                //Fill Charm species pT histos
                if (ilabelM>0 && TrackPDG == 11) {
                    //cout<<"Test4..................................."<<endl;
                
                    AliAODMCParticle *momPart = (AliAODMCParticle*)fMCarray->At(ilabelM); //get mom particle
                    pidM = TMath::Abs(momPart->GetPdgCode());
                    
                    fAllElecStackDiffPID->Fill(AODMCtrack->Pt());
                    if(pidM==411 || pidM==421 || pidM==413 || pidM==423 || pidM==431 || pidM==433 || pidM==4122){
                        fDElecStackDiffPID->Fill(AODMCtrack->Pt());
                    }
                    if(pidM==511 || pidM==521 || pidM==513 || pidM==523 || pidM==531 || pidM==533){
                        fBElecStackDiffPID->Fill(AODMCtrack->Pt());
                    }
                    
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
                fDMesonPDG->Fill(TrackPDG);
                fAllDMesonPt->Fill(AODMCtrack->Pt());
            
                if(TrackPDG==421) fD0MesonPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==421 && fromDStar==kTRUE) fD0MesonFromDStarPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==411) fDPlusMesonPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==431) fDsMesonPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==413) fDStarMesonPt->Fill(AODMCtrack->Pt());
                // cout<<"Total Number of Particles = "<<fNtotMCpart<<endl;
            }
        }
    }*/
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
    
    //Making event cuts
    fNevents->Fill(0); //all events
    if(NcontV>=2) fNevents->Fill(1); //>2 Trks
    if(NcontV>=2 && TMath::Abs(pVtx->GetZ())<=fVtxZCut) fNevents->Fill(2); //>2 Trks, with Vtx_Z cut
    //if(TMath::Abs(pVtx->GetZ())<=10.0) fNevents->Fill(3); //with Vtx_Z cut
    
    //make cut in Vtx_Z
    if(TMath::Abs(pVtx->GetZ())>fVtxZCut) return; //make cut in Vtx_Z
    fNevents->Fill(3);
    
    //Pile-up cuts
    // remove event 1
    const AliVVertex *spdVtx = fAOD->GetPrimaryVertexSPD();
    double covTrc[6],covSPD[6];
    pVtx->GetCovarianceMatrix(covTrc);
    spdVtx->GetCovarianceMatrix(covSPD);
    double dz = pVtx->GetZ()-spdVtx->GetZ();
    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
    //cout << TMath::Abs(dz) << " ; " << nsigTot << " ; " << nsigTrc << endl;
    if (fEnablePileupCut1){
        if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20)return;
    }
        
    // remove event2
    Int_t nTPCout=0;
    Float_t mTotV0=0;
    AliAODVZERO* v0data=(AliAODVZERO*) fAOD->GetVZEROData();
    Float_t mTotV0A=v0data->GetMTotV0A();
    Float_t mTotV0C=v0data->GetMTotV0C();
    mTotV0=mTotV0A+mTotV0C;
    Int_t ntracksEv = fAOD->GetNumberOfTracks();
    for(Int_t itrack=0; itrack<ntracksEv; itrack++) { // loop on tacks
        AliAODTrack * atrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
        if(!atrack) {AliFatal("Not a standard AOD");}
        if(atrack->GetID()<0)continue;
        if((atrack->GetFlags())&(AliVTrack::kTPCout)) nTPCout++;
        else continue;
    }
    Float_t mV0Cut=-2200.+(2.5*nTPCout)+(0.000012*nTPCout*nTPCout); //function to apply to pile-up rejection
    if(fEnablePileupRejVZEROTPCout)
    {
        if(mTotV0<mV0Cut) return;
    }
    
    fNevents->Fill(4); //with Vtx_Z cut + pile-up cuts
    
    //make cut in Ncontributors
    if(NcontV<2)return;
    
    fNevents->Fill(5); //with Vtx_Z cut + >2 Trks + pile-up cuts
    
    fVtX->Fill(pVtx->GetX());
    fVtY->Fill(pVtx->GetY());
    fVtZ->Fill(pVtx->GetZ());
    
    ///////////////////
    //Loop over Stack//
    ///////////////////
    
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
                
                if(TrackPDG == 11 && AODMCtrack->IsPhysicalPrimary() && AODMCtrack->Eta() >= fMinEta && AODMCtrack->Eta() <= fMaxEta) {
                    // cout<<"TESTINGGGGGGGGGGGG"<<endl;
                    fAllElecStack->Fill(AODMCtrack->Pt());
                    eleinStack++;
                }
                
                if(AODMCtrack->Eta() < fMinEta || AODMCtrack->Eta() > fMaxEta) continue;
                if (TrackPDG>500 && TrackPDG<599) {
                    fBMesonPt->Fill(AODMCtrack->Pt());
                }
                if (TrackPDG>5000 && TrackPDG<5999) {
                    fBBaryonPt->Fill(AODMCtrack->Pt());
                }
                
                //Fill Charm species pT histos
                if (ilabelM>0 && TrackPDG == 11) {
                    //cout<<"Test4..................................."<<endl;
                    
                    AliAODMCParticle *momPart = (AliAODMCParticle*)fMCarray->At(ilabelM); //get mom particle
                    pidM = TMath::Abs(momPart->GetPdgCode());
                    
                    fAllElecStackDiffPID->Fill(AODMCtrack->Pt());
                    if(pidM==411 || pidM==421 || pidM==413 || pidM==423 || pidM==431 || pidM==433 || pidM==4122){
                        fDElecStackDiffPID->Fill(AODMCtrack->Pt());
                    }
                    if(pidM==511 || pidM==521 || pidM==513 || pidM==523 || pidM==531 || pidM==533){
                        fBElecStackDiffPID->Fill(AODMCtrack->Pt());
                    }
                    
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
                if(TrackPDG==4122){
                    fLambdaCPt->Fill(AODMCtrack->Pt());
                    fLambdaCPtWeight->Fill(AODMCtrack->Pt(),fLcWeightVar1->GetBinContent(fLcWeightVar1->FindBin(AODMCtrack->Pt())));
                }
                if(TrackPDG==4132)fEtaCPt->Fill(AODMCtrack->Pt());
                
                if (TrackPDG<400 || TrackPDG>499) {
                    continue; //reject stuff that's not in charm meson range
                }
                //fDMesonPDG->Fill(TrackPDG);
                fAllDMesonPt->Fill(AODMCtrack->Pt());
                
                if(TrackPDG==421) {
                    fD0MesonPt->Fill(AODMCtrack->Pt());
                    fD0MesonPtWeight->Fill(AODMCtrack->Pt(),fDWeightNew->GetBinContent(fDWeightNew->FindBin(AODMCtrack->Pt())));
                }
                if(TrackPDG==421 && fromDStar==kTRUE) fD0MesonFromDStarPt->Fill(AODMCtrack->Pt());
                if(TrackPDG==411) {
                    fDPlusMesonPt->Fill(AODMCtrack->Pt());
                    fDPlusMesonPtWeight->Fill(AODMCtrack->Pt(),fDPlusWeightVar1->GetBinContent(fDPlusWeightVar1->FindBin(AODMCtrack->Pt())));
                }
                if(TrackPDG==431) {
                    fDsMesonPt->Fill(AODMCtrack->Pt());
                    fDsMesonPtWeight->Fill(AODMCtrack->Pt(),fDsWeightVar1->GetBinContent(fDsWeightVar1->FindBin(AODMCtrack->Pt())));
                }
                if(TrackPDG==413) fDStarMesonPt->Fill(AODMCtrack->Pt());
                // cout<<"Total Number of Particles = "<<fNtotMCpart<<endl;
            }
        }
    }
    
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
            
            if (fApplyTimeCut) {
                Float_t tof1 = clus->GetTOF()*1e+9; // ns
                if(!fMCarray && fUseTender && fFlagFillMCHistos)
                {
                    if(tof1<-30 || tof1>30)continue; // timing cut on data
                }
            }
            
            //fClsEAll->Fill(clus->E()); //E of all clusters
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
        
        if(track->Eta() < fMinEta || track->Eta() > fMaxEta) continue;
        
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
            
                if(fMCparticle->Eta() < fMinEta || fMCparticle->Eta() > fMaxEta) continue;
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
        
        if (fApplyITSLayer == 0) { //kAny (also included in kTrkGlobalNoDCA)
            if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) continue;
        }
        if (fApplyITSLayer == 1) { //kFirst
            if(!(track->HasPointOnITSLayer(0))) continue;
        }
        if (fApplyITSLayer == 2) { //kBoth
            if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) continue;
        }
        fITSLayerCheck->Fill(fApplyITSLayer);
        
        Bool_t kinkmotherpass = kTRUE;
        for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
            if(track->GetID() == listofmotherkink[kinkmother]) {
                kinkmotherpass = kFALSE;
                continue;
            }
        }
        if(!kinkmotherpass) continue; //kink rejection
        
        Double_t d0z0[2]={-999,-999}, cov[3];
        Double_t DCAxyCut = fDCAxyCut, DCAzCut = fDCAzCut;
        
        if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) continue;
        
        double phiMatchIts = track->Phi();
        
        if(track->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov)){
            fDCAxyz->Fill(d0z0[0],d0z0[1]);
            if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
        }
        Double_t DCA = d0z0[0]*track->Charge()*MagSign;
        
        //Refit
        if((!(track->GetStatus()&AliESDtrack::kITSrefit)|| (!(track->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        
        //Looser cuts to fill histo
        if(track->GetTPCNcls()>=60 && track->GetITSNcls()>=2 && track->GetTPCNCrossedRows()>=80){
            if(kTruElec == kTRUE) fElecAftLooseTrkCuts->Fill(track->Pt());
            if(kTruHFElec == kTRUE) fHFElecAftLooseTrkCuts->Fill(track->Pt());
            if(kTruBElec == kTRUE) fBElecAftLooseTrkCuts->Fill(track->Pt());
        }
        
        //Looser cuts and Shingo's method of getting the mother PID
        /*if(track->GetTPCNcls()>=60 && track->GetITSNcls()>=2 && track->GetTPCNCrossedRows()>=80){
            if(pdg==11 && ilabelM>0){
                fElecAftLooseTrkCutsDiffPID->Fill(track->Pt());
                if(pidM==411 || pidM==421 || pidM==413 || pidM==423 || pidM==431 || pidM==433 || pidM==4122){
                    fDElecAftLooseTrkCutsDiffPID->Fill(track->Pt());
                }
                if(pidM==511 || pidM==521 || pidM==513 || pidM==523 || pidM==531 || pidM==533){
                    fBElecAftLooseTrkCutsDiffPID->Fill(track->Pt());
                }
            }
        }*/
            
        //Tighter track cuts
        if(track->GetTPCNcls() < fNclusTPC) continue;
        if(track->GetITSNcls() < 3) continue;
        if(track->GetTPCNCrossedRows() < fNCrossRows) continue;
        if (fItsChi2>=0) {
            if(track->GetITSchi2() > fItsChi2) continue;
        }
        
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
        
        /*fDcaB4TrkMatch->Fill(track->Pt(),DCA,track->Phi());
        if (nsigma>-1.0 && nsigma<3.) {
            fElecDcaB4TrkMatch->Fill(track->Pt(),DCA,track->Phi());
        }*/
        
        fvalueBasic[0] = track->Pt();
        fvalueBasic[1] = DCA;
        fvalueBasic[2] = track->Phi();
        fvalueBasic[3] = track->Eta();
        fvalueBasic[4] = track->Charge();
        if (nsigma>-1.0 && nsigma<3.) {
            fvalueBasic[5] = 1;
        }else{
            fvalueBasic[5] = 0;
        }
        fvalueBasic[6] = 0;
        
        /*if(nsigma>-5.&&nsigma<-3.) {
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
        }*/
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
        
        //Fill sparse
        
        Double_t fPhiDiff = -999, fEtaDiff = -999;
        if(clustMatch && clustMatch->IsEMCAL())
        {
            //Double_t fPhiDiff = -999, fEtaDiff = -999;
            GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
            
            if(TMath::Abs(fPhiDiff) < fTrkMatch && TMath::Abs(fEtaDiff) < fTrkMatch) {
                fvalueBasic[6] = 1;
            }
        }
        fBasicSprs->Fill(fvalueBasic);
                
        if(clustMatch && clustMatch->IsEMCAL())
        {
            //Double_t fPhiDiff = -999, fEtaDiff = -999;
            //GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
            
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
            
            fClsEnoTimeCut->Fill(clustMatch->E());
            
            if (fApplyTimeCut) {
                Float_t tof = clustMatch->GetTOF()*1e+9; // ns
                if(!fMCarray && fUseTender && fFlagFillMCHistos)
                {
                    if(tof<-30 || tof>30)continue; // timing cut on data
                }
            }
            
            fClsE->Fill(clustMatch->E());
            
            //if(fClsTypeEMC) fClsEamEMCal->Fill(clustMatch->E());
            //if(fClsTypeDCAL) fClsEamDCal->Fill(clustMatch->E());
            
            if(kTruElec == kTRUE) fElecAftTrkMatch->Fill(track->Pt());
            if(kTruHFElec == kTRUE) fHFElecAftTrkMatch->Fill(track->Pt());
            if(kTruBElec == kTRUE) fBElecAftTrkMatch->Fill(track->Pt());
            
            fnSigmaAftTrkMatch->Fill(track->Pt(),nsigma);
            fEMCTrkMatch->Fill(fPhiDiff,fEtaDiff);
            fTrkClsPhi->Fill(fPhiDiff);
            fTrkClsEta->Fill(fEtaDiff);
            fClsPhi->Fill(emcphi);
            fClsEta->Fill(emceta);
            
            /*fDcaTrkMatch->Fill(track->Pt(),DCA,track->Phi());
            if (nsigma>-1.0 && nsigma<3.) {
                fElecDcaTrkMatch->Fill(track->Pt(),DCA,track->Phi());
            }*/
            
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
            Bool_t kFlagULS = kFALSE;
            Bool_t kFlagLS = kFALSE;
            Int_t fMomGen = 99;
            Double_t momPt = -99;
            Double_t momGamma = -99;
            Double_t momTime = -99;
            Int_t pidGM = -99;
            //Int_t ilabel = -99;
            Int_t ilabelM = -99;
            Int_t ilabelGM = -99;
            Double_t prodR = -99;
            
            //cout<<"TESTING0"<<endl;
            //if MC--------------------------------
            if (fFlagFillMCHistos) {
                if(ilabel>0 && fMCarray)
                {
                    //cout<<"TESTING1"<<endl;
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(ilabel);
                    pdg = fMCparticle->GetPdgCode(); //get pid of track
                    prodR = TMath::Sqrt(fMCparticle->Xv()*fMCparticle->Xv()+fMCparticle->Yv()*fMCparticle->Yv());
                
                    //cout<<"TESTING1"<<endl;
                
                    //if electron--------------------------------
                    if(TMath::Abs(pdg)==11){
                    
                        //cout<<"TESTING2"<<endl;
                        FindMother(fMCparticle, fpidSort, kEmbEta, kEmbPi0, kHijing, momPt, momGamma, momTime); //get its mom
                    
                        if (kHijing) fMomGen = 0;
                        if (kEmbPi0) fMomGen = 1;
                        if (kEmbEta) fMomGen = 2;
                    
                        //Fill template sparse
                        Double_t tempValue[5] = {-999,-999,-999,-999,-999};
                        tempValue[0] = track->Pt();
                        tempValue[1] = DCA;
                        tempValue[2] = fpidSort;
                        //tempValue[3] = fMomGen;
                        //tempValue[4] = momPt;
                        tempValue[3] = momGamma;
                        tempValue[4] = momTime;
                        fSprsTemplatesNoWeight->Fill(tempValue);
                    
                        //Took out Lambda_c (fpidSort==17) to the weighting
                        if (fpidSort==2||fpidSort==14||fpidSort==16) { //if from D meson
                            //cout<<"TESTING3"<<endl;
                            if (momPt>1 && momPt<50.) { //in proper pt range
                                //cout<<"TESTING4"<<endl;
                                //dWeight = fDWeight->GetBinContent(fDWeight->FindBin(momPt));
                                //fDTemplateWeight->Fill(track->Pt(), DCA, dWeight);
                                //fDTemplateNoWeight->Fill(track->Pt(), DCA);
                                
                                dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(momPt));
                                fSprsTemplatesWeight->Fill(tempValue,dWeight);
                                //fDTemplateWeightNew->Fill(track->Pt(), DCA, dWeight);
                                
                                dWeight = fDWeightVar1->GetBinContent(fDWeightVar1->FindBin(momPt));
                                //fDTemplateWeightVar1->Fill(track->Pt(), DCA, dWeight);
                                fSprsTemplatesWeightVar1->Fill(tempValue,dWeight);
                                
                                dWeight = fDWeightVar2->GetBinContent(fDWeightVar2->FindBin(momPt));
                                //fDTemplateWeightVar2->Fill(track->Pt(), DCA, dWeight);
                                fSprsTemplatesWeightVar2->Fill(tempValue,dWeight);
                            }
                        }else if (fpidSort==12) { //if from D0 meson
                            //cout<<"TESTING3"<<endl;
                            if (momPt>1 && momPt<50.) { //in proper pt range
                                //cout<<"TESTING4"<<endl;
                                //dWeight = fDWeight->GetBinContent(fDWeight->FindBin(momPt));
                                //fDTemplateWeight->Fill(track->Pt(), DCA, dWeight);
                                //fDTemplateNoWeight->Fill(track->Pt(), DCA);
                                
                                dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(momPt))*fD0TauWeight->Eval(momTime);
                                fSprsTemplatesWeight->Fill(tempValue,dWeight);
                                //fDTemplateWeightNew->Fill(track->Pt(), DCA, dWeight);
                                
                                dWeight = fDWeightVar1->GetBinContent(fDWeightVar1->FindBin(momPt));
                                //fDTemplateWeightVar1->Fill(track->Pt(), DCA, dWeight);
                                fSprsTemplatesWeightVar1->Fill(tempValue,dWeight);
                                
                                dWeight = fDWeightVar2->GetBinContent(fDWeightVar2->FindBin(momPt));
                                //fDTemplateWeightVar2->Fill(track->Pt(), DCA, dWeight);
                                fSprsTemplatesWeightVar2->Fill(tempValue,dWeight);
                            }
                        }
                        else if (fpidSort==11) { //if from D+ meson
                            //cout<<"TESTING3"<<endl;
                            if (momPt>1 && momPt<50.) { //in proper pt range
                                //cout<<"TESTING4"<<endl;
                                //dWeight = fDWeight->GetBinContent(fDWeight->FindBin(momPt));
                                //fDTemplateWeight->Fill(track->Pt(), DCA, dWeight);
                                //fDTemplateNoWeight->Fill(track->Pt(), DCA);
                                
                                dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(momPt))*fDPlusTauWeight->Eval(momTime);
                                fSprsTemplatesWeight->Fill(tempValue,dWeight);
                                //fDTemplateWeightNew->Fill(track->Pt(), DCA, dWeight);
                                
                                dWeight = fDPlusWeightVar1->GetBinContent(fDPlusWeightVar1->FindBin(momPt));
                                //fDTemplateWeightVar1->Fill(track->Pt(), DCA, dWeight);
                                fSprsTemplatesWeightVar1->Fill(tempValue,dWeight);
                                
                                dWeight = fDWeightVar2->GetBinContent(fDWeightVar2->FindBin(momPt));
                                //fDTemplateWeightVar2->Fill(track->Pt(), DCA, dWeight);
                                fSprsTemplatesWeightVar2->Fill(tempValue,dWeight);
                            }
                        }else if (fpidSort==15) { //if from Ds meson
                            //cout<<"TESTING3"<<endl;
                            if (momPt>1 && momPt<50.) { //in proper pt range
                                //cout<<"TESTING4"<<endl;
                                //dWeight = fDWeight->GetBinContent(fDWeight->FindBin(momPt));
                                //fDTemplateWeight->Fill(track->Pt(), DCA, dWeight);
                                //fDTemplateNoWeight->Fill(track->Pt(), DCA);
                                
                                dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(momPt))*fDsTauWeight->Eval(momTime);
                                fSprsTemplatesWeight->Fill(tempValue,dWeight);
                                //fDTemplateWeightNew->Fill(track->Pt(), DCA, dWeight);
                                
                                dWeight = fDsWeightVar1->GetBinContent(fDsWeightVar1->FindBin(momPt));
                                //fDTemplateWeightVar1->Fill(track->Pt(), DCA, dWeight);
                                fSprsTemplatesWeightVar1->Fill(tempValue,dWeight);
                                
                                dWeight = fDWeightVar2->GetBinContent(fDWeightVar2->FindBin(momPt));
                                //fDTemplateWeightVar2->Fill(track->Pt(), DCA, dWeight);
                                fSprsTemplatesWeightVar2->Fill(tempValue,dWeight);
                            }
                        }
                        else if (fpidSort==17) { //if from Lc
                            if (momPt>1 && momPt<50.) { //in proper pt range
                                dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(momPt));
                                fSprsTemplatesWeight->Fill(tempValue,dWeight);
                                
                                dWeight = fLcWeightVar1->GetBinContent(fLcWeightVar1->FindBin(momPt));
                                fSprsTemplatesWeightVar1->Fill(tempValue,dWeight);
                                
                                dWeight = fLcWeightVar2->GetBinContent(fLcWeightVar2->FindBin(momPt));
                                fSprsTemplatesWeightVar2->Fill(tempValue,dWeight);
                            }
                        }
                        else if (fpidSort==1) {//if from B meson
                            //cout<<"TESTING5"<<endl;
                            if (momPt>0. && momPt<50.) { //in proper pt range
                                //cout<<"TESTING6"<<endl;
                                //bWeight = fBWeight->GetBinContent(fBWeight->FindBin(momPt));
                                //fBTemplateWeight->Fill(track->Pt(), DCA, bWeight);
                                //fBTemplateNoWeight->Fill(track->Pt(), DCA);
                                
                                bWeight = fBWeightNew->GetBinContent(fBWeightNew->FindBin(momPt));
                                //fBTemplateWeightNew->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeight->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar1->GetBinContent(fBWeightVar1->FindBin(momPt));
                                //fBTemplateWeightVar1->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeightVar1->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar2->GetBinContent(fBWeightVar2->FindBin(momPt));
                                //fBTemplateWeightVar2->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeightVar2->Fill(tempValue,bWeight);
                            }
                        }else if (fpidSort==20) {//if from B+ meson
                            //cout<<"TESTING5"<<endl;
                            if (momPt>0. && momPt<50.) { //in proper pt range
                                bWeight = fBWeightNew->GetBinContent(fBWeightNew->FindBin(momPt))*fBPlusTauWeight->Eval(momTime);
                                //fBTemplateWeightNew->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeight->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar1->GetBinContent(fBWeightVar1->FindBin(momPt));
                                //fBTemplateWeightVar1->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeightVar1->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar2->GetBinContent(fBWeightVar2->FindBin(momPt));
                                //fBTemplateWeightVar2->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeightVar2->Fill(tempValue,bWeight);
                            }
                        }
                        else if (fpidSort==21) {//if from B0 meson
                            //cout<<"TESTING5"<<endl;
                            if (momPt>0. && momPt<50.) { //in proper pt range
                                bWeight = fBWeightNew->GetBinContent(fBWeightNew->FindBin(momPt))*fB0TauWeight->Eval(momTime);
                                //fBTemplateWeightNew->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeight->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar1->GetBinContent(fBWeightVar1->FindBin(momPt));
                                //fBTemplateWeightVar1->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeightVar1->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar2->GetBinContent(fBWeightVar2->FindBin(momPt));
                                //fBTemplateWeightVar2->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeightVar2->Fill(tempValue,bWeight);
                            }
                        }
                        else if (fpidSort==22) {//if from Bs meson
                            //cout<<"TESTING5"<<endl;
                            if (momPt>0. && momPt<50.) { //in proper pt range
                                bWeight = fBWeightNew->GetBinContent(fBWeightNew->FindBin(momPt))*fBsTauWeight->Eval(momTime);
                                //fBTemplateWeightNew->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeight->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar1->GetBinContent(fBWeightVar1->FindBin(momPt));
                                //fBTemplateWeightVar1->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeightVar1->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar2->GetBinContent(fBWeightVar2->FindBin(momPt));
                                //fBTemplateWeightVar2->Fill(track->Pt(), DCA, bWeight);
                                fSprsTemplatesWeightVar2->Fill(tempValue,bWeight);
                            }
                        }
                        else if (fpidSort==10) {//if from B baryon
                            if (momPt>0. && momPt<50.) { //in proper pt range
                                
                                bWeight = fBWeightNew->GetBinContent(fBWeightNew->FindBin(momPt));
                                fSprsTemplatesWeight->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar1->GetBinContent(fBWeightVar1->FindBin(momPt));
                                fSprsTemplatesWeightVar1->Fill(tempValue,bWeight);
                                
                                bWeight = fBWeightVar2->GetBinContent(fBWeightVar2->FindBin(momPt));
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
                        /*if(fpidSort==12){
                            ilabelM = fMCparticle->GetMother(); //get MC label for e Mom
                            AliAODMCParticle *MCpartM = (AliAODMCParticle*)fMCarray->At(ilabelM); //get e mom particle
                            ilabelGM = MCpartM->GetMother(); //get MC label for grandma
                            //if no grandma.............
                            if (ilabelGM<0) {
                                fPromptD0DCANoWeight->Fill(track->Pt(),DCA);
                                if (momPt>1 && momPt<50.) {
                                    dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(momPt));
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
                                        dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(momPt));
                                        fD0FromDStarDCAWeight->Fill(track->Pt(),DCA,dWeight);
                                    }
                                }
                            }
                        }*/
                    
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
                        
                        Double_t prodR = TMath::Sqrt(fMCparticle->Xv()*fMCparticle->Xv()+fMCparticle->Yv()*fMCparticle->Yv());
                        
                        //if mom is Pi0--------------------------------
                        if(fpidSort==3) {
                            //fPi0DCA->Fill(track->Pt(),DCA);
                            if(!kEmbEta && !kEmbPi0 && kHijing) {
                                fPi0HijingDCA->Fill(track->Pt(),DCA);
                                fPi0HijingPt->Fill(track->Pt());
                                //ComboDenomWeight->Fill(track->Pt());
                                //ComboDenomNoWeight->Fill(track->Pt());
                            }
                            if(kEmbPi0) {
                                fWeight = fPi0Weight->Eval(momPt);
                                fEnhPi0DCA->Fill(track->Pt(),DCA);
                                fEnhPi0WeightedPt->Fill(track->Pt(),fWeight);
                                fComboDenomWeight->Fill(track->Pt(),prodR,fpidSort,fWeight);
                                fComboDenomNoWeight->Fill(track->Pt(),prodR,fpidSort);
                            }
                            if(kEmbEta) {
                                fWeight = fEtaWeight->Eval(momPt);
                                fEnhEtaDCA->Fill(track->Pt(),DCA);
                                fEnhEtaWeightedPt->Fill(track->Pt(),fWeight);
                                fComboDenomWeight->Fill(track->Pt(),prodR,fpidSort,fWeight);
                                fComboDenomNoWeight->Fill(track->Pt(),prodR,fpidSort);
                            }
                        
                        }
                        //if Eta--------------------------------
                        if(fpidSort==4){
                            //fEtaDCA->Fill(track->Pt(),DCA);
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
                                fComboDenomWeight->Fill(track->Pt(),prodR,fpidSort,fWeight);
                                fComboDenomNoWeight->Fill(track->Pt(),prodR,fpidSort);
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
                                fComboDenomWeight->Fill(track->Pt(),prodR,fpidSort,fWeight);
                                fComboDenomNoWeight->Fill(track->Pt(),prodR,fpidSort);
                            }
                            if(kEmbEta) {
                                fWeight = fEtaWeight->Eval(momPt);
                                fEnhPhotonDCA->Fill(track->Pt(),DCA);
                                fEnhPhotonWeightedPt->Fill(track->Pt(),fWeight);
                                fComboDenomWeight->Fill(track->Pt(),prodR,fpidSort,fWeight);
                                fComboDenomNoWeight->Fill(track->Pt(),prodR,fpidSort);
                            }
                        }
                    
                    
                        InvMassCheckMC(i, track, d0z0, MagSign, kHijing, kEmbEta, kEmbPi0, kFlagReco,fWeight, fpidSort, prodR);
                    
                        //cout<<"TESTING2"<<endl;
                        if(kFlagReco){
                        
                            //if mom is Pi0--------------------------------
                            if(fpidSort==3) {
                                //fPi0DCA->Fill(track->Pt(),DCA);
                                if(!kEmbEta && !kEmbPi0 && kHijing) {
                                    fULSHijingPi0->Fill(track->Pt()); //pi0 mama
                                    //ComboNumWeight->Fill(track->Pt());
                                    //ComboNumNoWeight->Fill(track->Pt());
                                }
                                if(kEmbPi0) {
                                    fWeight = fPi0Weight->Eval(momPt);
                                    fULSWeightEnhPi0->Fill(track->Pt(),fWeight); //pi0 mama
                                    fComboNumWeight->Fill(track->Pt(),prodR,fpidSort,fWeight);
                                    fComboNumNoWeight->Fill(track->Pt(),prodR,fpidSort);
                                }
                                if(kEmbEta) {
                                    fWeight = fEtaWeight->Eval(momPt);
                                    fULSWeightEnhEta->Fill(track->Pt(),fWeight); //eta mama
                                    fComboNumWeight->Fill(track->Pt(),prodR,fpidSort,fWeight);
                                    fComboNumNoWeight->Fill(track->Pt(),prodR,fpidSort);
                                }
                            
                            }
                            //if Eta--------------------------------
                            if(fpidSort==4){
                                //fEtaDCA->Fill(track->Pt(),DCA);
                                if(!kEmbEta && !kEmbPi0 && kHijing) {
                                    fULSHijingEta->Fill(track->Pt()); //eta mama
                                    //ComboNumWeight->Fill(track->Pt());
                                    //ComboNumNoWeight->Fill(track->Pt());
                                }
                                if(kEmbEta) {
                                    fWeight = fEtaWeight->Eval(momPt);
                                    fULSWeightEnhEta->Fill(track->Pt(),fWeight); //eta mama
                                    fComboNumWeight->Fill(track->Pt(),prodR,fpidSort,fWeight);
                                    fComboNumNoWeight->Fill(track->Pt(),prodR,fpidSort);
                                }
                            }
                            //if photon--------------------------------
                            if(fpidSort==5){
                                if(!kEmbEta && !kEmbPi0 && kHijing) {
                                    fULSHijingPhoton->Fill(track->Pt()); //photon mama
                                    fPhotonHijingTagDCA->Fill(track->Pt(),DCA);
                                    //ComboNumWeight->Fill(track->Pt());
                                    //ComboNumNoWeight->Fill(track->Pt());
                                }
                                if(kEmbPi0) {
                                    fWeight = fPi0Weight->Eval(momPt);
                                    fEnhPhotonTagDCA->Fill(track->Pt(),DCA);
                                    fULSEnhPhoton->Fill(track->Pt(),fWeight); //photon mama
                                    fComboNumWeight->Fill(track->Pt(),prodR,fpidSort,fWeight);
                                    fComboNumNoWeight->Fill(track->Pt(),prodR,fpidSort);
                                }
                                if(kEmbEta) {
                                    fWeight = fEtaWeight->Eval(momPt);
                                    fEnhPhotonTagDCA->Fill(track->Pt(),DCA);
                                    fULSEnhPhoton->Fill(track->Pt(),fWeight); //photon mama
                                    fComboNumWeight->Fill(track->Pt(),prodR,fpidSort,fWeight);
                                    fComboNumNoWeight->Fill(track->Pt(),prodR,fpidSort);
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
            
            /*if((nsigma>fMinNSigCut) && (nsigma<3)) {
                if ((M20>0.01) && (M20<fMaxM20Cut)) {
                    fInclElecEoPnoShift->Fill(track->Pt(),EovP);
                }
            }
            if(nsigma<-4.) {
                if(M20>0.01 && M20<fMaxM20Cut) {
                    fHadronEoPnoShift->Fill(track->Pt(),EovP);
                }
            }*/
            
            if(fShiftEoP && fMCarray && fFlagFillMCHistos)  // E/p MC mean shift correction
            {
                if(fCentralityMin==30 && fCentralityMax==50)
                {
                    EovP += 0.04; //30-50%
                }
                else if(fCentralityMin==60 && fCentralityMax==80)
                {
                    EovP += 0.045; //60-80% (tuned up to 18 GeV/c)
                }
                else
                {
                    EovP += 0.0;
                }
            }
            
            if(fFlagFillMCHistos && fFlagFillSprs) {
                Double_t fvalueElectron[6] = {-999,-999,-999,-999,-999,-999};
                fvalueElectron[0] = track->Pt();
                fvalueElectron[1] = nsigma;
                fvalueElectron[2] = EovP;
                if (fApplyM02Cut) {
                    fvalueElectron[3] = M02;
                }else{
                    fvalueElectron[3] = M20;
                }
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
                if (fApplyM02Cut) {
                    fvalueElectron[3] = M02;
                }else{
                    fvalueElectron[3] = M20;
                }
                fvalueElectron[4] = DCA;
                fElectronSprs->Fill(fvalueElectron);
            }
            
            ///////////////////////////
            // Hadron Contam. Histos //
            ///////////////////////////
            if (fApplyM02Cut) {
                if(nsigma<-4.) {
                    if(M02>0.01 && M02<0.7) {
                        fHadronEoP->Fill(track->Pt(),EovP);
                        if (!fApplyHadEoPCut) {
                            fHadronDCA->Fill(track->Pt(),DCA);
                        }
                        if(fApplyHadEoPCut && EovP>fMinEoPCut && EovP<fMaxEoPCut){
                            fHadronDCA->Fill(track->Pt(),DCA);
                        }
                    }
                }
            } else {
                if(nsigma<-4.) {
                    if(M20>0.01 && M20<fMaxM20Cut) {
                        fHadronEoP->Fill(track->Pt(),EovP);
                        if (!fApplyHadEoPCut) {
                            fHadronDCA->Fill(track->Pt(),DCA);
                        }
                        if(fApplyHadEoPCut && EovP>fMinEoPCut && EovP<fMaxEoPCut){
                            fHadronDCA->Fill(track->Pt(),DCA);
                        }
                    }
                }
            }
            
            if(nsigma<-4.) {
                    fHadronEoPNoM20->Fill(track->Pt(),EovP);
            }
            //if(nsigma>-5.&&nsigma<-3.) {
            //    fHadronCamDCA->Fill(track->Pt(),d0z0[0]*MagSign);
            //}
            
            ///////////////////
            // Electron Cuts //
            ///////////////////
            if (fApplyM02Cut) {
                if((EovP>fMinEoPCut) && (EovP<fMaxEoPCut)) {
                    fnSigaftEoPCut->Fill(track->Pt(),nsigma);
                    if (M02>0.01 && M02<0.7) {
                        fnSigaftM20EoPCut->Fill(track->Pt(),nsigma);
                    }
                }
            } else {
                if((EovP>fMinEoPCut) && (EovP<fMaxEoPCut)) {
                    fnSigaftEoPCut->Fill(track->Pt(),nsigma);
                    if (M20>0.01 && M20<0.35) {
                        fnSigaftM20EoPCut->Fill(track->Pt(),nsigma);
                    }
                }
            }
            
            /*if((EovP>fMinEoPCut) && (EovP<fMaxEoPCut)) {
                fTPCElecEoP->Fill(track->Pt(),EovP);
                fnSigaftSysEoPCut->Fill(track->Pt(),nsigma);
            }*/
            
            if((EovP>fMinEoPCut) && (EovP<fMaxEoPCut)){
                if(kTruElec == kTRUE) fElecAftEoP->Fill(track->Pt());
                if(kTruHFElec == kTRUE) fHFElecAftEoP->Fill(track->Pt());
                if(kTruBElec == kTRUE) fBElecAftEoP->Fill(track->Pt());
            }
            if((nsigma>fMinNSigCut) && (nsigma<3)) fInclElecEoPNoM20->Fill(track->Pt(),EovP);
            
            //Apply M20 cut for electrons
            if (fApplyM02Cut) {
                if((M02<0.01) || (M02>0.7)) continue;
            } else {
                if((M20<0.01) || (M20>fMaxM20Cut)) continue;
            }
            
            //fElecEoPnoSig->Fill(track->Pt(),EovP);
            
            if((nsigma>fMinNSigCut) && (nsigma<3)) fInclElecEoP->Fill(track->Pt(),EovP);
            
            //Apply E/p Cut for electrons
            if((EovP<fMinEoPCut) || (EovP>fMaxEoPCut)) continue;
            //fnSigaftSysM20EoPCut->Fill(track->Pt(),nsigma);
            
            if(kTruElec == kTRUE) fElecAftEMCeID->Fill(track->Pt());
            if(kTruHFElec == kTRUE) fHFElecAftEMCeID->Fill(track->Pt());
            if(kTruBElec == kTRUE) fBElecAftEMCeID->Fill(track->Pt());
            
            //Apply TPC nSigma cut for electrons
            if((nsigma<fMinNSigCut) || (nsigma>3)) continue;
            
            if(kTruElec == kTRUE) fElecAftTPCeID->Fill(track->Pt());
            if(kTruHFElec == kTRUE) fHFElecAftTPCeID->Fill(track->Pt());
            if(kTruBElec == kTRUE) fBElecAftTPCeID->Fill(track->Pt());
            
            //if(fClsTypeDCAL) fClsEamElecDC->Fill(clustMatch->E());
            //if(fClsTypeEMC) fClsEamElecEMC->Fill(clustMatch->E());
            
            /////////////////////////
            // Plot Reco Electrons //
            /////////////////////////
            fInclElecDCA->Fill(track->Pt(),DCA);
            //fInclElecDCAnoSign->Fill(track->Pt(),d0z0[0]);
            
            //if(!fFlagFillMCHistos){
            //fWeight = 1;
            //InvMassCheckData(i, track, d0z0, MagSign, fWeight);
            //}
            
            //Apply weights to inv mass histos for closure test
            if(fFlagFillMCHistos){
                if(ilabel>0 && fMCarray)
                {
                    //cout<<"BIG TEST__________________________________________________________"<<endl;
                    //if mom is Pi0--------------------------------
                    if(fpidSort==3) {
                        //cout<<"Test 2"<<endl;
                        if(kEmbPi0) {
                            //cout<<"Test 3"<<endl;
                            fWeight = fPi0Weight->Eval(momPt);
                            //InvMassCheckData(i, track, d0z0, MagSign, fWeight);
                            FillULSSparse(i, track, d0z0, MagSign, fWeight, prodR, fpidSort);
                        }
                        if(kEmbEta) {
                            fWeight = fEtaWeight->Eval(momPt);
                            //InvMassCheckData(i, track, d0z0, MagSign, fWeight);
                            FillULSSparse(i, track, d0z0, MagSign, fWeight, prodR, fpidSort);
                        }
                        
                    }
                    //if Eta--------------------------------
                    if(fpidSort==4){
                        //cout<<"Test 2a"<<endl;
                        if(kEmbEta) {
                            //cout<<"Test 3a"<<endl;
                            fWeight = fEtaWeight->Eval(momPt);
                            //InvMassCheckData(i, track, d0z0, MagSign, fWeight);
                            FillULSSparse(i, track, d0z0, MagSign, fWeight, prodR, fpidSort);
                        }
                    }
                    //if photon--------------------------------
                    if(fpidSort==5){
                        //cout<<"Test 2b"<<endl;
                        if(kEmbPi0) {
                            //cout<<"Test 3b"<<endl;
                            fWeight = fPi0Weight->Eval(momPt);
                            //InvMassCheckData(i, track, d0z0, MagSign, fWeight);
                            FillULSSparse(i, track, d0z0, MagSign, fWeight, prodR, fpidSort);
                        }
                        if(kEmbEta) {
                            fWeight = fEtaWeight->Eval(momPt);
                            //InvMassCheckData(i, track, d0z0, MagSign, fWeight);
                            FillULSSparse(i, track, d0z0, MagSign, fWeight, prodR, fpidSort);
                        }
                    }
                }
            }else{
                fWeight = 1;
                InvMassCheckData(i, track, d0z0, MagSign, fWeight);
            }
            
            //Fill DCA true for closure test
            if (fFlagFillMCHistos) {
                if(ilabel>0 && fMCarray)
                {
                    if(TMath::Abs(pdg)==11){
                        //Fill closure test sparse
                        Double_t closValue[4] = {-999,-999,-999,-999};
                        closValue[0] = track->Pt();
                        closValue[1] = DCA;
                        closValue[2] = fpidSort;
                        closValue[3] = prodR;
                        fSprsClosureTest->Fill(closValue);
                    
                        //if mom is Pi0--------------------------------
                        if(fpidSort==3) {
                            if(kEmbPi0) {
                                fWeight = fPi0Weight->Eval(momPt);
                                fSprsClosureTestWeight->Fill(closValue,fWeight);
                            }
                            if(kEmbEta) {
                                fWeight = fEtaWeight->Eval(momPt);
                                fSprsClosureTestWeight->Fill(closValue,fWeight);
                            }
                        
                        }
                        //if Eta--------------------------------
                        if(fpidSort==4){
                            if(kEmbEta) {
                                fWeight = fEtaWeight->Eval(momPt);
                                fSprsClosureTestWeight->Fill(closValue,fWeight);
                            }
                        }
                        //if photon--------------------------------
                        if(fpidSort==5){
                            if(kEmbPi0) {
                                fWeight = fPi0Weight->Eval(momPt);
                                fSprsClosureTestWeight->Fill(closValue,fWeight);
                            }
                            if(kEmbEta) {
                                fWeight = fEtaWeight->Eval(momPt);
                                fSprsClosureTestWeight->Fill(closValue,fWeight);
                            }
                        }
                    }
                }
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
void AliAnalysisTaskTPCCalBeauty::FindMother(AliAODMCParticle* part, Int_t &fpidSort, Bool_t &kEmbEta, Bool_t &kEmbPi0, Bool_t &kHijing, Double_t &momPt, Double_t &momGamma,Double_t &momTime)
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
    Double_t decayL = -99;
    
    //cout<<"TESTING3"<<endl;
    
    //Get the mother, grandma, and great grandma pid
    ilabelM = part->GetMother(); //get MC label for Mom
    if(ilabelM>0){
        AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(ilabelM); //get mom particle
        pidM = TMath::Abs(partM->GetPdgCode()); //ask for the Mom's pid
        momPt = partM->Pt();
        if (partM->M()>0) {
            momGamma = partM->E()/partM->M();
        }
        
        decayL = TMath::Sqrt(TMath::Power(partM->Xv()-part->Xv(),2)+TMath::Power(partM->Yv()-part->Yv(),2)+TMath::Power(partM->Zv()-part->Zv(),2));
        momTime = (10000*decayL*partM->M())/partM->P();
        
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
        }else if(pidM==521){
            fpidSort = 20; //Mom is B+
        }else if(pidM==511){
            fpidSort = 21; //Mom is B0
        }else if(pidM==531){
            fpidSort = 22; //Mom is Bs
        }
        
        //Using Jonghan's method to find beauty feeddown for the D mesons
        if((int(pidM/100.)%10) == 4 || (int(pidM/1000.)%10) == 4) {
            
            // iterate until you find B hadron as a mother or become top ancestor
            AliAODMCParticle *dummyPart; //dummy particle for iteration
            AliAODMCParticle *dummyPartDaughter; //2nd dummy particle for iteration
            
            int grandMaPDG;
            
            for (int i=1; i<100; i++){
                int jLabel = partM->GetMother();
                if (jLabel == -1) {
                    break;
                }
                if ((jLabel<0)){
                    AliDebug(1, "Stack label is negative, return\n");
                    break;
                }
                
                // if there is an ancestor
                if(!(dummyPart = dynamic_cast<AliAODMCParticle *>(fMCarray->At(TMath::Abs(jLabel))))) {
                    break;
                }
                //cout<<"Before GetDaughterFirst"<<endl;
                dummyPartDaughter = dynamic_cast<AliAODMCParticle *>(fMCarray->At(dummyPart->GetDaughterFirst()));
                //cout<<"After GetDaughterFirst"<<endl;
                
                grandMaPDG = TMath::Abs(dummyPart->GetPdgCode());
                if (grandMaPDG>500 && grandMaPDG<599){
                    fpidSort = 1; //B mother feeddown
                    if (grandMaPDG==521) {
                        fpidSort = 20; //B+
                    }else if(grandMaPDG==511){
                        fpidSort = 21; //Mom is B0
                    }else if(grandMaPDG==531){
                        fpidSort = 22; //Mom is Bs
                    }
                    momPt = dummyPart->Pt();
                    if (partM->M()>0) {
                        momGamma = dummyPart->E()/dummyPart->M();
                    }
                    decayL = TMath::Sqrt(TMath::Power(dummyPart->Xv()-dummyPartDaughter->Xv(),2)+TMath::Power(dummyPart->Yv()-dummyPartDaughter->Yv(),2)+TMath::Power(dummyPart->Zv()-dummyPartDaughter->Zv(),2));
                    momTime = (10000*decayL*dummyPart->M())/dummyPart->P();
                    break;
                }
                if (grandMaPDG>5000 && grandMaPDG<5999){
                    fpidSort = 10; //b baryon mother feeddown
                    momPt = dummyPart->Pt();
                    if (partM->M()>0) {
                        momGamma = dummyPart->E()/dummyPart->M();
                    }
                    decayL = TMath::Sqrt(TMath::Power(dummyPart->Xv()-dummyPartDaughter->Xv(),2)+TMath::Power(dummyPart->Yv()-dummyPartDaughter->Yv(),2)+TMath::Power(dummyPart->Zv()-dummyPartDaughter->Zv(),2));
                    momTime = (10000*decayL*dummyPart->M())/dummyPart->P();
                    break;
                }
                partM = dummyPart;
            } // end of iteration
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
                        if (partM->M()>0) {
                            momGamma = partGM->E()/partGM->M();
                        }
                        decayL = TMath::Sqrt(TMath::Power(partGM->Xv()-partM->Xv(),2)+TMath::Power(partGM->Yv()-partM->Yv(),2)+TMath::Power(partGM->Zv()-partM->Zv(),2));
                        momTime = (10000*decayL*partGM->M())/partGM->P();
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
                        if (partM->M()>0) {
                            momGamma = partGM->E()/partGM->M();
                        }
                        decayL = TMath::Sqrt(TMath::Power(partGM->Xv()-partM->Xv(),2)+TMath::Power(partGM->Yv()-partM->Yv(),2)+TMath::Power(partGM->Zv()-partM->Zv(),2));
                        momTime = (10000*decayL*partGM->M())/partGM->P();
                    }
                }
                if(pidGM==111){
                    if(ilabelM >= fNembMCpi0 && ilabelM<fNembMCeta && ilabelGGM<0) { //GMa is emb pi0
                        kEmbEta = kFALSE;
                        kEmbPi0 = kTRUE;
                        momPt = partGM->Pt(); //make pi0 pt mompt for weighting
                        if (partM->M()>0) {
                            momGamma = partGM->E()/partGM->M();
                        }
                        decayL = TMath::Sqrt(TMath::Power(partGM->Xv()-partM->Xv(),2)+TMath::Power(partGM->Yv()-partM->Yv(),2)+TMath::Power(partGM->Zv()-partM->Zv(),2));
                        momTime = (10000*decayL*partGM->M())/partGM->P();
                    }//GMa is pi0
                }
            }
            
            
            //check if D grandma is B
            /*if(pidM>400 && pidM<499){
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
            }*/
            if(ilabelGGM>0){
                AliAODMCParticle *partGGM = (AliAODMCParticle*)fMCarray->At(ilabelGGM); // get GGMa particle
                Int_t pidGGM = TMath::Abs(partGGM->GetPdgCode()); //ask for ggma's pid
                ilabelGGGM = partGGM->GetMother();//get MC for Great Grandma
                
                //check if D great grandma is B
                /*if(pidM>400 && pidM<499){
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
                }*/
                
                //check if gamma great grandma is eta
                if(pidM==22){
                    if(pidGM==111){ //grandma is pion
                        if(pidGGM==221){ //great grandma is eta
                            if(ilabelGGM >= fNembMCeta && ilabelGGM<fNtotMCpart && ilabelGGGM<0) {
                                kEmbEta = kTRUE; //GMa is enh eta
                                kEmbPi0 = kFALSE;
                                momPt = partGGM->Pt(); //make eta pt mompt for weighting
                                if (partM->M()>0) {
                                    momGamma = partGGM->E()/partGGM->M();
                                }
                                decayL = TMath::Sqrt(TMath::Power(partGGM->Xv()-partGM->Xv(),2)+TMath::Power(partGGM->Yv()-partGM->Yv(),2)+TMath::Power(partGGM->Zv()-partGM->Zv(),2));
                                momTime = (10000*decayL*partGGM->M())/partGGM->P();
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
void AliAnalysisTaskTPCCalBeauty::InvMassCheckData(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign, Double_t fWeight)
{
    // Flags photonic electrons with inv mass cut
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t d0z0Asso[2]={-999,-999}, covAsso[3];
    //Double_t DCAxyCut = 0.25, DCAzCut = 1;
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
        if(trackAsso->GetTPCNcls() < fAssoTPCnCls) continue;
        
        //Refit
        if((!(trackAsso->GetStatus()&AliESDtrack::kITSrefit)|| (!(trackAsso->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        
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
            if(TMath::Abs(d0z0Asso[0]) > fAssoDCAxy || TMath::Abs(d0z0Asso[1]) > fAssoDCAz) continue;
        
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
        
        if(fFlagLS && mass>fMinMass && mass<fMaxMass) Nls++;
        if(fFlagULS && mass>fMinMass && mass<fMaxMass) Nuls++;
        
        if (fFlagULS && mass>fMinMass && mass<fMaxMass && track->Pt()>1) {
            fULSdcaBelow->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign);
            fULSdcaBelowWeight->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign,fWeight);
                
        }else if(fFlagLS && mass>fMinMass && mass<fMaxMass && track->Pt()>1){
            fLSdcaBelow->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign);
            fLSdcaBelowWeight->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign,fWeight);
        }
        
    }
    
    //fPhotonicElecYield->Fill(track->Pt(),Nuls-Nls);
}
//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::InvMassCheckMC(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign, Bool_t kHijing, Bool_t kEmbEta, Bool_t kEmbPi0, Bool_t &kFlagReco, Double_t fWeight, Int_t fpidSort, Double_t prodRadius)
{
    // Flags photonic electrons with inv mass cut
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t d0z0Asso[2]={-999,-999}, covAsso[3];
    //Double_t DCAxyCut = 0.25, DCAzCut = 1;
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
        if(trackAsso->GetTPCNcls() < fAssoTPCnCls) continue;
        
        //Refit
        if((!(trackAsso->GetStatus()&AliESDtrack::kITSrefit)|| (!(trackAsso->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        
        nsigmaAsso = fpidResponse->NumberOfSigmasTPC(trackAsso, AliPID::kElectron);
        ptAsso = trackAsso->Pt();
        chargeAsso = trackAsso->Charge();
        
        //Some cuts on the associated track
        if(ptAsso < fMinPtAssoCut) continue;
        if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
        if(nsigmaAsso < fMinNSigAssoCut || nsigmaAsso > 3) continue;
        
        if(trackAsso->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0Asso, covAsso))
            if(TMath::Abs(d0z0Asso[0]) > fAssoDCAxy || TMath::Abs(d0z0Asso[1]) > fAssoDCAz) continue;
        
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
        
        if(fFlagLS && mass>fMinMass && mass<fMaxMass) Nls++;
        if(fFlagULS && mass>fMinMass && mass<fMaxMass) Nuls++;
        
        //CHANGED FROM pt>1
        if (fFlagULS && mass>fMinMass && mass<fMaxMass) {
            kFlagReco = kTRUE;
        }else if(fFlagLS && mass>fMinMass && mass<fMaxMass){
            kFlagReco = kFALSE;
        }
        
        /*if (fFlagULS && mass>fMinMass && mass<fMaxMass && track->Pt()>1) {
            fULSdcaBelow->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign,prodRadius);
            
        }else if(fFlagLS && mass>fMinMass && mass<fMaxMass && track->Pt()>1){
            fLSdcaBelow->Fill(track->Pt(),d0z0[0]*track->Charge()*MagSign,prodRadius);
        }*/
        
    }
    
}
//________________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::FillULSSparse(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign, Double_t fWeight,Double_t prodRadius,Int_t pidSort)
{
    // Flags photonic electrons with inv mass cut
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t d0z0Asso[2]={-999,-999}, covAsso[3];
    //Double_t DCAxyCut = 0.25, DCAzCut = 1;
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
        if(trackAsso->GetTPCNcls() < fAssoTPCnCls) continue;
        
        //Refit
        if((!(trackAsso->GetStatus()&AliESDtrack::kITSrefit)|| (!(trackAsso->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        
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
            if(TMath::Abs(d0z0Asso[0]) > fAssoDCAxy || TMath::Abs(d0z0Asso[1]) > fAssoDCAz) continue;
        
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
        
        if(fFlagLS && mass>fMinMass && mass<fMaxMass) Nls++;
        if(fFlagULS && mass>fMinMass && mass<fMaxMass) Nuls++;
        
        Double_t sprsFillValue[4] = {-999,-999,-999,-999};
        sprsFillValue[0] = track->Pt();
        sprsFillValue[1] = d0z0[0]*track->Charge()*MagSign;
        sprsFillValue[2] = pidSort;
        sprsFillValue[3] = prodRadius;
        
        if (fFlagULS && mass>fMinMass && mass<fMaxMass && track->Pt()>1) {
            fSprsULSdca->Fill(sprsFillValue);
            fSprsULSdcaWeight->Fill(sprsFillValue,fWeight);
                
        }else if(fFlagLS && mass>fMinMass && mass<fMaxMass && track->Pt()>1){
            fSprsLSdca->Fill(sprsFillValue);
            fSprsLSdcaWeight->Fill(sprsFillValue,fWeight);
        }
        
    }
    
    //fPhotonicElecYield->Fill(track->Pt(),Nuls-Nls);
}
//________________________________________________________________________
//________________________________________________________________________
/*void AliAnalysisTaskTPCCalBeauty::SetBmesonTauWeight(TF2 *BPlus, TF2 *B0, TF2 *Bs)
{
    fBPlusTauWeight = (TF2 *)BPlus->Clone();
    fB0TauWeight = (TF2 *)B0->Clone();
    fBsTauWeight = (TF2 *)Bs->Clone();
}*/
//_____________________________________________________________________
void AliAnalysisTaskTPCCalBeauty::Terminate(Option_t *)
{
    // terminate
}
