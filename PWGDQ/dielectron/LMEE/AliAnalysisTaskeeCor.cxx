#include "TChain.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TList.h"
#include "TString.h"
#include "TRandom.h"

#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliMCEvent.h>
#include <AliInputEventHandler.h>
#include <AliAODInputHandler.h>
#include <AliMCEventHandler.h>
#include "AliStack.h"
#include "AliVVertex.h"
#include "AliMCParticle.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliGenEventHeader.h"

#include "AliAODEvent.h"
#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "AliAODTrack.h"

#include "AliAnalysisTaskeeCor.h"

class AliAnalysisTaskeeCor;
using namespace std;

ClassImp(AliAnalysisTaskeeCor)

//________________________________________________________________________
AliAnalysisTaskeeCor::AliAnalysisTaskeeCor()
: AliAnalysisTaskSE()
,mcEvent(0)
,aodEvent(0)
,fOutputList(0)
,fOutputListGen(0)
,fOutputListRec(0)
,fOutputListGenCC(0)
,fOutputListRecCC(0)
,fOutputListGenBB(0)
,fOutputListRecBB(0)
,fEventStat(0)
,fMotherPDGCC(0)
,fHistEleGenPt(0)
,fHistPosGenPt(0)
,fHistEleGenEtaPhi(0)
,fHistPosGenEtaPhi(0)
,fHistEleGenPhi(0)
,fHistPosGenPhi(0)
,fHistEleGenDCAxy(0)
,fHistPosGenDCAxy(0)
,fHistEleGenPtCC(0)
,fHistPosGenPtCC(0)
,fHistEleGenEtaPhiCC(0)
,fHistPosGenEtaPhiCC(0)
,fHistPairsGen(0)
,fHistPairsCCgen(0)
,fHistPairsCCgen2(0)
,fHistEleRecPt(0)
,fHistPosRecPt(0)
,fHistEleRecEtaPhi(0)
,fHistPosRecEtaPhi(0)
,fHistEleRecPhi(0)
,fHistPosRecPhi(0)
,fHistEleRecDCAxy(0)
,fHistPosRecDCAxy(0)
,fHistEleRecDCAxySig(0)
,fHistPosRecDCAxySig(0)
,fHistEleRecPtDCAxySig(0)
,fHistPosRecPtDCAxySig(0)
,fHistEleRecPtCC(0)
,fHistPosRecPtCC(0)
,fHistEleRecEtaPhiCC(0)
,fHistPosRecEtaPhiCC(0)
,fHistNclsTPC(0)
,fHistNclsSTPC(0)
,fHistNcrTPC(0)
,fHistNclsITS(0)
,fHistDCAxy(0)
,fHistDCAz(0)
,fHistChi2perNDF(0)
,fHistTPCnSigmaEle(0)
,fHistTPCnSigmaPio(0)
,fHistTOFnSigmaEle(0)
,fHistPairsRec(0)
,fHistPairsCCrec(0)
,fHistPairsCCrec2(0)
,fMotherPDGBB(0)
,fHistEleGenPtBB(0)
,fHistPosGenPtBB(0)
,fHistEleGenEtaPhiBB(0)
,fHistPosGenEtaPhiBB(0)
,fHistPairsBBgen(0)
,fHistPairsBBgenLS(0)
,fHistPairsBBgen2(0)
,fHistPairsBBgenLS2(0)
,fHistEleRecPtBB(0)
,fHistPosRecPtBB(0)
,fHistEleRecEtaPhiBB(0)
,fHistPosRecEtaPhiBB(0)
,fHistPairsBBrec(0)
,fHistPairsBBrecLS(0)
,fHistPairsBBrec2(0)
,fHistPairsBBrecLS2(0)
,fHistPairsGenPt(0)
,fHistPairsGenMass(0)
,fHistPairsGenDphi(0)
,fHistPairsGenMassPt(0)
,fHistPairsGenPhiPt(0)
,fHistPairsGenPtMasPhi(0)
,fHistPairsRecPt(0)
,fHistPairsRecMass(0)
,fHistPairsRecDphi(0)
,fHistPairsRecMassPt(0)
,fHistPairsRecPhiPt(0)
,fHistPairsRecPtMasPhi(0)
,fHistPairsGenMassPt2(0)
,fHistPairsGenPhiPt2(0)
,fHistPairsRecMassPt2(0)
,fHistPairsRecPhiPt2(0)
,fHistPairsGenPtCC(0)
,fHistPairsGenMassCC(0)
,fHistPairsGenDphiCC(0)
,fHistPairsGenMassPtCC(0)
,fHistPairsGenPhiPtCC(0)
,fHistPairsGenPtMasPhiCC(0)
,fHistPairsRecPtCC(0)
,fHistPairsRecMassCC(0)
,fHistPairsRecDphiCC(0)
,fHistPairsRecMassPtCC(0)
,fHistPairsRecPhiPtCC(0)
,fHistPairsRecPtMasPhiCC(0)
,fHistPairsGenMassPt2CC(0)
,fHistPairsGenPhiPt2CC(0)
,fHistPairsRecMassPt2CC(0)
,fHistPairsRecPhiPt2CC(0)
,fHistPairsGenPtBB(0)
,fHistPairsGenMassBB(0)
,fHistPairsGenDphiBB(0)
,fHistPairsGenMassPtBB(0)
,fHistPairsGenPhiPtBB(0)
,fHistPairsGenMassPt2BB(0)
,fHistPairsGenPhiPt2BB(0)
,fHistPairsGenPtMasPhiBB(0)
,fHistPairsGenPtLSBB(0)
,fHistPairsGenMassLSBB(0)
,fHistPairsGenDphiLSBB(0)
,fHistPairsGenMassPtLSBB(0)
,fHistPairsGenPhiPtLSBB(0)
,fHistPairsGenMassPt2LSBB(0)
,fHistPairsGenPhiPt2LSBB(0)
,fHistPairsGenPtMasPhiLSBB(0)
,fHistPairsRecPtBB(0)
,fHistPairsRecMassBB(0)
,fHistPairsRecDphiBB(0)
,fHistPairsRecMassPtBB(0)
,fHistPairsRecPhiPtBB(0)
,fHistPairsRecMassPt2BB(0)
,fHistPairsRecPhiPt2BB(0)
,fHistPairsRecPtMasPhiBB(0)
,fHistPairsRecPtLSBB(0)
,fHistPairsRecMassLSBB(0)
,fHistPairsRecDphiLSBB(0)
,fHistPairsRecMassPtLSBB(0)
,fHistPairsRecPhiPtLSBB(0)
,fHistPairsRecMassPt2LSBB(0)
,fHistPairsRecPhiPt2LSBB(0)
,fHistPairsRecPtMasPhiLSBB(0)
,fHistPairsDCAxy(0)
,fHistPairsDCAxySig(0)
,fHistPairsPtDCAxySig(0)
,fHistPairsDCAxyCC(0)
,fHistPairsDCAxySigCC(0)
,fHistPairsPtDCAxySigCC(0)
,fHistPairsDCAxyBB(0)
,fHistPairsDCAxySigBB(0)
,fHistPairsPtDCAxySigBB(0)
,fHistPairsDCAxyLSBB(0)
,fHistPairsDCAxySigLSBB(0)
,fHistPairsPtDCAxySigLSBB(0)
,fHistPairsGenDCAxy(0)
,fHistPairsGenDCAxySig(0)
,fHistPairsGenPtDCAxySig(0)
,fHistPairsGenDCAxyCC(0)
,fHistPairsGenDCAxySigCC(0)
,fHistPairsGenPtDCAxySigCC(0)
,fHistPairsGenDCAxyBB(0)
,fHistPairsGenDCAxySigBB(0)
,fHistPairsGenPtDCAxySigBB(0)
,fHistPairsGenDCAxyLSBB(0)
,fHistPairsGenDCAxySigLSBB(0)
,fHistPairsGenPtDCAxySigLSBB(0)
,fHistEleRecDCAxyCC(0)
,fHistPosRecDCAxyCC(0)
,fHistEleRecDCAxyBB(0)
,fHistPosRecDCAxyBB(0)
,fPtRec_Gen(0)
,fPtRecOverGen(0)
,fPRec_Gen(0)
,fEtaRec_Gen(0)
,fElePhiRec_Gen(0)
,fPosPhiRec_Gen(0)
,fPairDCARec_Gen(0)
,fEleDCARes(0)
,fPosDCARes(0)
,fEleDCARes2(0)
,fPosDCARes2(0)
,fEleDCARes3(0)
,fPosDCARes3(0)
,fPIDResponse(0)
,fTrackCuts(0)
,fPtMinCut(0.2)
,fRecabPID(kFALSE)
,fTPCmean(0x0)
,fTPCwidth(0x0)
,fTOFmean(0x0)
,fTOFwidth(0x0)
,fUseSmearing(kFALSE)
,fDCASmearingByMath(kFALSE)
,fDCASmearingByMaps(kFALSE)
,fDCASmearingByPars(kFALSE)
,fDCAparSmr(kFALSE)
,fPtSmr(0x0)
,fEtaSmr(0x0)
,fPhiEleSmr(0x0)
,fPhiPosSmr(0x0)
,fDCASmrMap0(0x0)
,fDCASmrMap1(0x0)
,fDCASmrCen(0x0)
,fDCASmrSig(0x0)
,fDCASmrMax(0x0)
,fDCASmrFromMC(kFALSE)
,fTPCmin(-3.)
,fTPCmax(3.)
,fTOFmin(-3.)
,fTOFmax(3.)
,fTPCminPio(-100.)
,fTPCmaxPio(4.)
,fTPCminPro(-4.)
,fTPCmaxPro(4.)
,fTPCminKao(-4.)
,fTPCmaxKao(4.)
,fFillPIDrecMaps(kTRUE)
,fFillSmrMaps(kTRUE)
,TPCnSigmaEle_Eta_P_lin(0)
,TOFnSigmaEle_Eta_P_lin(0)
,TPCnSigmaEle_Eta_P_lin2(0)
,TOFnSigmaEle_Eta_P_lin2(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//________________________________________________________________________
AliAnalysisTaskeeCor::AliAnalysisTaskeeCor(const char* name)
: AliAnalysisTaskSE(name)
,mcEvent(0)
,aodEvent(0)
,fOutputList(0)
,fOutputListGen(0)
,fOutputListRec(0)
,fOutputListGenCC(0)
,fOutputListRecCC(0)
,fOutputListGenBB(0)
,fOutputListRecBB(0)
,fEventStat(0)
,fMotherPDGCC(0)
,fHistEleGenPt(0)
,fHistPosGenPt(0)
,fHistEleGenEtaPhi(0)
,fHistPosGenEtaPhi(0)
,fHistEleGenPhi(0)
,fHistPosGenPhi(0)
,fHistEleGenDCAxy(0)
,fHistPosGenDCAxy(0)
,fHistEleGenPtCC(0)
,fHistPosGenPtCC(0)
,fHistEleGenEtaPhiCC(0)
,fHistPosGenEtaPhiCC(0)
,fHistPairsGen(0)
,fHistPairsCCgen(0)
,fHistPairsCCgen2(0)
,fHistEleRecPt(0)
,fHistPosRecPt(0)
,fHistEleRecEtaPhi(0)
,fHistPosRecEtaPhi(0)
,fHistEleRecPhi(0)
,fHistPosRecPhi(0)
,fHistEleRecDCAxy(0)
,fHistPosRecDCAxy(0)
,fHistEleRecDCAxySig(0)
,fHistPosRecDCAxySig(0)
,fHistEleRecPtDCAxySig(0)
,fHistPosRecPtDCAxySig(0)
,fHistEleRecPtCC(0)
,fHistPosRecPtCC(0)
,fHistEleRecEtaPhiCC(0)
,fHistPosRecEtaPhiCC(0)
,fHistNclsTPC(0)
,fHistNclsSTPC(0)
,fHistNcrTPC(0)
,fHistNclsITS(0)
,fHistDCAxy(0)
,fHistDCAz(0)
,fHistChi2perNDF(0)
,fHistTPCnSigmaEle(0)
,fHistTPCnSigmaPio(0)
,fHistTOFnSigmaEle(0)
,fHistPairsRec(0)
,fHistPairsCCrec(0)
,fHistPairsCCrec2(0)
,fMotherPDGBB(0)
,fHistEleGenPtBB(0)
,fHistPosGenPtBB(0)
,fHistEleGenEtaPhiBB(0)
,fHistPosGenEtaPhiBB(0)
,fHistPairsBBgen(0)
,fHistPairsBBgenLS(0)
,fHistPairsBBgen2(0)
,fHistPairsBBgenLS2(0)
,fHistEleRecPtBB(0)
,fHistPosRecPtBB(0)
,fHistEleRecEtaPhiBB(0)
,fHistPosRecEtaPhiBB(0)
,fHistPairsBBrec(0)
,fHistPairsBBrecLS(0)
,fHistPairsBBrec2(0)
,fHistPairsBBrecLS2(0)
,fHistPairsGenPt(0)
,fHistPairsGenMass(0)
,fHistPairsGenDphi(0)
,fHistPairsGenMassPt(0)
,fHistPairsGenPhiPt(0)
,fHistPairsGenPtMasPhi(0)
,fHistPairsRecPt(0)
,fHistPairsRecMass(0)
,fHistPairsRecDphi(0)
,fHistPairsRecMassPt(0)
,fHistPairsRecPhiPt(0)
,fHistPairsRecPtMasPhi(0)
,fHistPairsGenMassPt2(0)
,fHistPairsGenPhiPt2(0)
,fHistPairsRecMassPt2(0)
,fHistPairsRecPhiPt2(0)
,fHistPairsGenPtCC(0)
,fHistPairsGenMassCC(0)
,fHistPairsGenDphiCC(0)
,fHistPairsGenMassPtCC(0)
,fHistPairsGenPhiPtCC(0)
,fHistPairsGenPtMasPhiCC(0)
,fHistPairsRecPtCC(0)
,fHistPairsRecMassCC(0)
,fHistPairsRecDphiCC(0)
,fHistPairsRecMassPtCC(0)
,fHistPairsRecPhiPtCC(0)
,fHistPairsRecPtMasPhiCC(0)
,fHistPairsGenMassPt2CC(0)
,fHistPairsGenPhiPt2CC(0)
,fHistPairsRecMassPt2CC(0)
,fHistPairsRecPhiPt2CC(0)
,fHistPairsGenPtBB(0)
,fHistPairsGenMassBB(0)
,fHistPairsGenDphiBB(0)
,fHistPairsGenMassPtBB(0)
,fHistPairsGenPhiPtBB(0)
,fHistPairsGenMassPt2BB(0)
,fHistPairsGenPhiPt2BB(0)
,fHistPairsGenPtMasPhiBB(0)
,fHistPairsGenPtLSBB(0)
,fHistPairsGenMassLSBB(0)
,fHistPairsGenDphiLSBB(0)
,fHistPairsGenMassPtLSBB(0)
,fHistPairsGenPhiPtLSBB(0)
,fHistPairsGenMassPt2LSBB(0)
,fHistPairsGenPhiPt2LSBB(0)
,fHistPairsGenPtMasPhiLSBB(0)
,fHistPairsRecPtBB(0)
,fHistPairsRecMassBB(0)
,fHistPairsRecDphiBB(0)
,fHistPairsRecMassPtBB(0)
,fHistPairsRecPhiPtBB(0)
,fHistPairsRecMassPt2BB(0)
,fHistPairsRecPhiPt2BB(0)
,fHistPairsRecPtMasPhiBB(0)
,fHistPairsRecPtLSBB(0)
,fHistPairsRecMassLSBB(0)
,fHistPairsRecDphiLSBB(0)
,fHistPairsRecMassPtLSBB(0)
,fHistPairsRecPhiPtLSBB(0)
,fHistPairsRecMassPt2LSBB(0)
,fHistPairsRecPhiPt2LSBB(0)
,fHistPairsRecPtMasPhiLSBB(0)
,fHistPairsDCAxy(0)
,fHistPairsDCAxySig(0)
,fHistPairsPtDCAxySig(0)
,fHistPairsDCAxyCC(0)
,fHistPairsDCAxySigCC(0)
,fHistPairsPtDCAxySigCC(0)
,fHistPairsDCAxyBB(0)
,fHistPairsDCAxySigBB(0)
,fHistPairsPtDCAxySigBB(0)
,fHistPairsDCAxyLSBB(0)
,fHistPairsDCAxySigLSBB(0)
,fHistPairsPtDCAxySigLSBB(0)
,fHistPairsGenDCAxy(0)
,fHistPairsGenDCAxySig(0)
,fHistPairsGenPtDCAxySig(0)
,fHistPairsGenDCAxyCC(0)
,fHistPairsGenDCAxySigCC(0)
,fHistPairsGenPtDCAxySigCC(0)
,fHistPairsGenDCAxyBB(0)
,fHistPairsGenDCAxySigBB(0)
,fHistPairsGenPtDCAxySigBB(0)
,fHistPairsGenDCAxyLSBB(0)
,fHistPairsGenDCAxySigLSBB(0)
,fHistPairsGenPtDCAxySigLSBB(0)
,fHistEleRecDCAxyCC(0)
,fHistPosRecDCAxyCC(0)
,fHistEleRecDCAxyBB(0)
,fHistPosRecDCAxyBB(0)
,fPtRec_Gen(0)
,fPtRecOverGen(0)
,fPRec_Gen(0)
,fEtaRec_Gen(0)
,fElePhiRec_Gen(0)
,fPosPhiRec_Gen(0)
,fPairDCARec_Gen(0)
,fEleDCARes(0)
,fPosDCARes(0)
,fEleDCARes2(0)
,fPosDCARes2(0)
,fEleDCARes3(0)
,fPosDCARes3(0)
,fPIDResponse(0)
,fTrackCuts(0)
,fPtMinCut(0.2)
,fRecabPID(kFALSE)
,fTPCmean(0x0)
,fTPCwidth(0x0)
,fTOFmean(0x0)
,fTOFwidth(0x0)
,fUseSmearing(kFALSE)
,fDCASmearingByMath(kFALSE)
,fDCASmearingByMaps(kFALSE)
,fDCASmearingByPars(kFALSE)
,fDCAparSmr(kFALSE)
,fPtSmr(0x0)
,fEtaSmr(0x0)
,fPhiEleSmr(0x0)
,fPhiPosSmr(0x0)
,fDCASmrMap0(0x0)
,fDCASmrMap1(0x0)
,fDCASmrCen(0x0)
,fDCASmrSig(0x0)
,fDCASmrMax(0x0)
,fDCASmrFromMC(kFALSE)
,fTPCmin(-3.)
,fTPCmax(3.)
,fTOFmin(-3.)
,fTOFmax(3.)
,fTPCminPio(-100.)
,fTPCmaxPio(4.)
,fTPCminPro(-4.)
,fTPCmaxPro(4.)
,fTPCminKao(-4.)
,fTPCmaxKao(4.)
,TPCnSigmaEle_Eta_P_lin(0)
,TOFnSigmaEle_Eta_P_lin(0)
,TPCnSigmaEle_Eta_P_lin2(0)
,TOFnSigmaEle_Eta_P_lin2(0)
,fFillPIDrecMaps(kTRUE)
,fFillSmrMaps(kTRUE)
{
    // constructor
    DefineInput(0, TChain::Class());    // a 'chain' of events created by the analysis manager automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
    // one can add more output objects by calling DefineOutput(2, classname::Class())
    // make sure to connect them properly in AddTask and to call PostData() for all of them!

}

//________________________________________________________________________
AliAnalysisTaskeeCor::~AliAnalysisTaskeeCor()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;
    }
}

//________________________________________________________________________
void AliAnalysisTaskeeCor::UserCreateOutputObjects()
{

    // retrieve PID object from the input handler
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliAODInputHandler* aodHandler=(AliAODInputHandler*)man->GetInputEventHandler();
    fPIDResponse = aodHandler->GetPIDResponse();

    // Create histograms
    // Called once at the start of the analysis

    fOutputList = new TList();
    fOutputList->SetName("Output");
    fOutputList->SetOwner(kTRUE); // the list is owner of all objects it contains and will delete them

    fOutputListGen = new TList();
    fOutputListGen->SetName("OutputGEN");
    fOutputListGen->SetOwner(kTRUE); // the list is owner of all objects it contains and will delete them

    fOutputListRec = new TList();
    fOutputListRec->SetName("OutputRec");
    fOutputListRec->SetOwner(kTRUE); // the list is owner of all objects

    fOutputListGenCC = new TList();
    fOutputListGenCC->SetName("OutputGENCC");
    fOutputListGenCC->SetOwner(kTRUE); // the list is owner of all objects it contains and will delete them

    fOutputListRecCC = new TList();
    fOutputListRecCC->SetName("OutputRecCC");
    fOutputListRecCC->SetOwner(kTRUE); // the list is owner of all objects it contains and will delete them

    fOutputListGenBB = new TList();
    fOutputListGenBB->SetName("OutputGENBB");
    fOutputListGenBB->SetOwner(kTRUE); // the list is owner of all objects it contains and will delete them

    fOutputListRecBB = new TList();
    fOutputListRecBB->SetName("OutputRecBB");
    fOutputListRecBB->SetOwner(kTRUE); // the list is owner of all objects it contains and will delete them

    fEventStat = new TH1I("fEventStat", "Event gen. statistics", kNbinsEvent, 0, kNbinsEvent);
    fEventStat->GetXaxis()->SetBinLabel(1,"Pythia CC");
    fEventStat->GetXaxis()->SetBinLabel(2,"Pythia BB");
    fEventStat->GetXaxis()->SetBinLabel(3,"Pythia B");
    fEventStat->GetXaxis()->SetBinLabel(4,"Other");
    fOutputList->Add(fEventStat);

    fMotherPDGCC = new TH1I("fMotherPDGCC", "Mother PDG code", 5, 0, 5);
    fMotherPDGCC->GetXaxis()->SetBinLabel(1,"411");
    fMotherPDGCC->GetXaxis()->SetBinLabel(2,"421");
    fMotherPDGCC->GetXaxis()->SetBinLabel(3,"431");
    fMotherPDGCC->GetXaxis()->SetBinLabel(4,"4122");
    fMotherPDGCC->GetXaxis()->SetBinLabel(5,"Other");
    fOutputListGenCC->Add(fMotherPDGCC);

    fMotherPDGBB = new TH1I("fMotherPDGBB", "Mother PDG code", 5, 0, 5);
    fMotherPDGBB->GetXaxis()->SetBinLabel(1,"mesons");
    fMotherPDGBB->GetXaxis()->SetBinLabel(2,"barions");
    fMotherPDGBB->GetXaxis()->SetBinLabel(3,"D-meson");
    fMotherPDGBB->GetXaxis()->SetBinLabel(4,"Other");
    fOutputListGenBB->Add(fMotherPDGBB);

    // gen histos
    fHistEleGenPt = new TH1F("fHistEleGenPt","Gen electrons pt",100,0.0,10.0);
    fHistEleGenPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListGen->Add(fHistEleGenPt);

    fHistPosGenPt = new TH1F("fHistPosGenPt","Gen positrons pt",100,0.0,10.0);
    fHistPosGenPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListGen->Add(fHistPosGenPt);

    fHistEleGenEtaPhi = new TH2F("fHistEleGenEtaPhi","Gen electrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistEleGenEtaPhi->GetXaxis()->SetTitle("#eta");
    fHistEleGenEtaPhi->GetYaxis()->SetTitle("#varphi");
    fOutputListGen->Add(fHistEleGenEtaPhi);

    fHistPosGenEtaPhi = new TH2F("fHistPosGenEtaPhi","Gen positrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistPosGenEtaPhi->GetXaxis()->SetTitle("#eta");
    fHistPosGenEtaPhi->GetYaxis()->SetTitle("#varphi");
    fOutputListGen->Add(fHistPosGenEtaPhi);

    fHistEleGenPhi = new TH1F("fHistEleGenPhi","Gen electrons phi",100,0.0,2*TMath::Pi());
    fHistEleGenPhi->GetXaxis()->SetTitle("#varphi");
    fOutputListGen->Add(fHistEleGenPhi);

    fHistPosGenPhi = new TH1F("fHistPosGenPhi","Gen positrons phi",100,0.0,2*TMath::Pi());
    fHistPosGenPhi->GetXaxis()->SetTitle("#varphi");
    fOutputListGen->Add(fHistPosGenPhi);

    fHistEleGenDCAxy = new TH1F("fHistEleGenDCAxy","Gen electrons DCAxy",800,-1.0,1.0);
    fHistEleGenDCAxy->GetXaxis()->SetTitle("DCAxy");
    fOutputListGen->Add(fHistEleGenDCAxy);

    fHistPosGenDCAxy = new TH1F("fHistPosGenDCAxy","Gen positrons DCAxy",800,-1.0,1.0);
    fHistPosGenDCAxy->GetXaxis()->SetTitle("DCAxy");
    fOutputListGen->Add(fHistPosGenDCAxy);

    fHistEleRecDCAxyCC = new TH1F("fHistEleRecDCAxyCC","Gen electrons DCAxy",800,-1.0,1.0);
    fHistEleRecDCAxyCC->GetXaxis()->SetTitle("DCAxy");
    fOutputListRecCC->Add(fHistEleRecDCAxyCC);

    fHistPosRecDCAxyCC = new TH1F("fHistPosRecDCAxyCC","Gen positrons DCAxy",800,-1.0,1.0);
    fHistPosRecDCAxyCC->GetXaxis()->SetTitle("DCAxy");
    fOutputListRecCC->Add(fHistPosRecDCAxyCC);

    fHistEleRecDCAxyBB = new TH1F("fHistEleRecDCAxyBB","Gen electrons DCAxy",800,-1.0,1.0);
    fHistEleRecDCAxyBB->GetXaxis()->SetTitle("DCAxy");
    fOutputListRecBB->Add(fHistEleRecDCAxyBB);

    fHistPosRecDCAxyBB = new TH1F("fHistPosRecDCAxyBB","Gen positrons DCAxy",800,-1.0,1.0);
    fHistPosRecDCAxyBB->GetXaxis()->SetTitle("DCAxy");
    fOutputListRecBB->Add(fHistPosRecDCAxyBB);

    fHistEleGenPtCC = new TH1F("fHistEleGenPtCC","Gen electrons pt",100,0.0,10.0);
    fHistEleGenPtCC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListGenCC->Add(fHistEleGenPtCC);

    fHistPosGenPtCC = new TH1F("fHistPosGenPtCC","Gen positrons pt",100,0.0,10.0);
    fHistPosGenPtCC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListGenCC->Add(fHistPosGenPtCC);

    fHistEleGenEtaPhiCC = new TH2F("fHistEleGenEtaPhiCC","Gen electrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistEleGenEtaPhiCC->GetXaxis()->SetTitle("#eta");
    fHistEleGenEtaPhiCC->GetYaxis()->SetTitle("#varphi");
    fOutputListGenCC->Add(fHistEleGenEtaPhiCC);

    fHistPosGenEtaPhiCC = new TH2F("fHistPosGenEtaPhiCC","Gen positrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistPosGenEtaPhiCC->GetXaxis()->SetTitle("#eta");
    fHistPosGenEtaPhiCC->GetYaxis()->SetTitle("#varphi");
    fOutputListGenCC->Add(fHistPosGenEtaPhiCC);

    fHistPairsGen = new TH2F("fHistPairsGen", "Gen e+e- from all sources", 50, 0.,5.0,10,0.,10.);
    fHistPairsGen->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsGen->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListGen->Add(fHistPairsGen);

    fHistPairsCCgen = new TH2F("fHistPairsCCgen", "Gen e+e- from cc", 500, 0.,5.0,100,0.,10.);
    fHistPairsCCgen->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsCCgen->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListGenCC->Add(fHistPairsCCgen);

    fHistPairsCCgen2 = new TH2F("fHistPairsCCgen2", "Gen e+e- from cc", 50, 0.,5.0,10,0.,10.);
    fHistPairsCCgen2->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsCCgen2->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListGenCC->Add(fHistPairsCCgen2);

    fHistEleGenPtBB = new TH1F("fHistEleGenPtBB","Gen electrons pt",100,0.0,10.0);
    fHistEleGenPtBB->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListGenBB->Add(fHistEleGenPtBB);

    fHistPosGenPtBB = new TH1F("fHistPosGenPtBB","Gen positrons pt",100,0.0,10.0);
    fHistPosGenPtBB->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListGenBB->Add(fHistPosGenPtBB);

    fHistEleGenEtaPhiBB = new TH2F("fHistEleGenEtaPhiBB","Gen electrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistEleGenEtaPhiBB->GetXaxis()->SetTitle("#eta");
    fHistEleGenEtaPhiBB->GetYaxis()->SetTitle("#varphi");
    fOutputListGenBB->Add(fHistEleGenEtaPhiBB);

    fHistPosGenEtaPhiBB = new TH2F("fHistPosGenEtaPhiBB","Gen positrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistPosGenEtaPhiBB->GetXaxis()->SetTitle("#eta");
    fHistPosGenEtaPhiBB->GetYaxis()->SetTitle("#varphi");
    fOutputListGenBB->Add(fHistPosGenEtaPhiBB);

    fHistPairsBBgen = new TH2F("fHistPairsBBgen", "Gen e+e- from bb", 500, 0.,5.0,100,0.,10.);
    fHistPairsBBgen->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsBBgen->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListGenBB->Add(fHistPairsBBgen);

    fHistPairsBBgenLS = new TH2F("fHistPairsBBgenLS", "Gen e+e- from bb", 500, 0.,5.0,100,0.,10.);
    fHistPairsBBgenLS->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsBBgenLS->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListGenBB->Add(fHistPairsBBgenLS);

    fHistPairsBBgen2 = new TH2F("fHistPairsBBgen2", "Gen e+e- from bb", 50, 0.,5.0,10,0.,10.);
    fHistPairsBBgen2->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsBBgen2->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListGenBB->Add(fHistPairsBBgen2);

    fHistPairsBBgenLS2 = new TH2F("fHistPairsBBgenLS2", "Gen e+e- from bb", 50, 0.,5.0,10,0.,10.);
    fHistPairsBBgenLS2->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsBBgenLS2->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListGenBB->Add(fHistPairsBBgenLS2);

    // rec histos
    fHistEleRecPt = new TH1F("fHistEleRecPt","Rec electrons pt",100,0.0,10.0);
    fHistEleRecPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListRec->Add(fHistEleRecPt);

    fHistPosRecPt = new TH1F("fHistPosRecPt","Rec positrons pt",100,0.0,10.0);
    fHistPosRecPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListRec->Add(fHistPosRecPt);

    fHistEleRecEtaPhi = new TH2F("fHistEleRecEtaPhi","Rec electrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistEleRecEtaPhi->GetXaxis()->SetTitle("#eta");
    fHistEleRecEtaPhi->GetXaxis()->SetTitle("#varphi");
    fOutputListRec->Add(fHistEleRecEtaPhi);

    fHistPosRecEtaPhi = new TH2F("fHistPosRecEtaPhi","Rec positrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistPosRecEtaPhi->GetXaxis()->SetTitle("#eta");
    fHistPosRecEtaPhi->GetXaxis()->SetTitle("#varphi");
    fOutputListRec->Add(fHistPosRecEtaPhi);

    fHistEleRecPhi = new TH1F("fHistEleRecPhi","Rec electrons phi",100,0.0,2*TMath::Pi());
    fHistEleRecPhi->GetXaxis()->SetTitle("#varphi");
    fOutputListRec->Add(fHistEleRecPhi);

    fHistPosRecPhi = new TH1F("fHistPosRecPhi","Rec positrons phi",100,0.0,2*TMath::Pi());
    fHistPosRecPhi->GetXaxis()->SetTitle("#varphi");
    fOutputListRec->Add(fHistPosRecPhi);

    fHistEleRecDCAxy = new TH1F("fHistEleRecDCAxy","Rec electrons DCAxy",800,-1.0,1.0);
    fHistEleRecDCAxy->GetXaxis()->SetTitle("DCAxy");
    fOutputListRec->Add(fHistEleRecDCAxy);

    fHistPosRecDCAxy = new TH1F("fHistPosRecDCAxy","Rec positrons DCAxy",800,-1.0,1.0);
    fHistPosRecDCAxy->GetXaxis()->SetTitle("DCAxy");
    fOutputListRec->Add(fHistPosRecDCAxy);
    
    fHistEleRecDCAxySig = new TH1F("fHistEleRecDCAxySig","Rec electrons DCAxySig",800,-20.0,20.0);
    fHistEleRecDCAxySig->GetXaxis()->SetTitle("DCAxy");
    fOutputListRec->Add(fHistEleRecDCAxySig);

    fHistPosRecDCAxySig = new TH1F("fHistPosRecDCAxySig","Rec positrons DCAxySig",800,-20.0,20.0);
    fHistPosRecDCAxySig->GetXaxis()->SetTitle("DCAxy");
    fOutputListRec->Add(fHistPosRecDCAxySig);
    
    /*fHistEleRecPtDCAxySig = new TH1F("fHistEleRecPtDCAxySig","Rec electrons DCAxySig",800,-20.0,20.0);
    fHistEleRecPtDCAxySig->GetXaxis()->SetTitle("DCAxy");
    fOutputListRec->Add(fHistEleRecPtDCAxySig);

    fHistPosRecPtDCAxySig = new TH1F("fHistPosRecPtDCAxySig","Rec positrons DCAxySig",800,-20.0,20.0);
    fHistPosRecPtDCAxySig->GetXaxis()->SetTitle("DCAxy");
    fOutputListRec->Add(fHistPosRecPtDCAxySig);*/

    fHistEleRecPtCC = new TH1F("fHistEleRecPtCC","Rec electrons pt",100,0.0,10.0);
    fHistEleRecPtCC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListRecCC->Add(fHistEleRecPtCC);

    fHistPosRecPtCC = new TH1F("fHistPosRecPtCC","Rec positrons pt",100,0.0,10.0);
    fHistPosRecPtCC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListRecCC->Add(fHistPosRecPtCC);

    fHistEleRecEtaPhiCC = new TH2F("fHistEleRecEtaPhiCC","Rec electrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistEleRecEtaPhiCC->GetXaxis()->SetTitle("#eta");
    fHistEleRecEtaPhiCC->GetXaxis()->SetTitle("#varphi");
    fOutputListRecCC->Add(fHistEleRecEtaPhiCC);

    fHistPosRecEtaPhiCC = new TH2F("fHistPosRecEtaPhiCC","Rec positrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistPosRecEtaPhiCC->GetXaxis()->SetTitle("#eta");
    fHistPosRecEtaPhiCC->GetXaxis()->SetTitle("#varphi");
    fOutputListRecCC->Add(fHistPosRecEtaPhiCC);

    fHistEleRecPtBB = new TH1F("fHistEleRecPtBB","Rec electrons pt",100,0.0,10.0);
    fHistEleRecPtBB->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListRecBB->Add(fHistEleRecPtBB);

    fHistPosRecPtBB = new TH1F("fHistPosRecPtBB","Rec positrons pt",100,0.0,10.0);
    fHistPosRecPtBB->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fOutputListRecBB->Add(fHistPosRecPtBB);

    fHistEleRecEtaPhiBB = new TH2F("fHistEleRecEtaPhiBB","Rec electrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistEleRecEtaPhiBB->GetXaxis()->SetTitle("#eta");
    fHistEleRecEtaPhiBB->GetXaxis()->SetTitle("#varphi");
    fOutputListRecBB->Add(fHistEleRecEtaPhiBB);

    fHistPosRecEtaPhiBB = new TH2F("fHistPosRecEtaPhiBB","Rec positrons eta phi",100,-1.0,1.0,100,0.0,2*TMath::Pi());
    fHistPosRecEtaPhiBB->GetXaxis()->SetTitle("#eta");
    fHistPosRecEtaPhiBB->GetXaxis()->SetTitle("#varphi");
    fOutputListRecBB->Add(fHistPosRecEtaPhiBB);

    // tracking info

    fHistNclsTPC = new TH1I("fHistNclsTPC","N clusters TPC",200,0,200);
    fHistNclsTPC->GetXaxis()->SetTitle("N clusters TPC");
    fOutputListRec->Add(fHistNclsTPC);

    fHistNclsSTPC = new TH1I("fHistNclsSTPC","N shared clusters TPC",200,0,200);
    fHistNclsSTPC->GetXaxis()->SetTitle("N shared clusters TPC");
    fOutputListRec->Add(fHistNclsSTPC);

    fHistNcrTPC = new TH1I("fHistNcrTPC","N crossed rows TPC",200,0,200);
    fHistNcrTPC->GetXaxis()->SetTitle("N crossed rows TPC");
    fOutputListRec->Add(fHistNcrTPC);

    fHistNclsITS = new TH1I("fHistNclsITS","N clusters ITS",10,0,10);
    fHistNclsITS->GetXaxis()->SetTitle("N clustersITS");
    fOutputListRec->Add(fHistNclsITS);

    fHistDCAxy = new TH1F("fHistDCAxy","DCAxy",800,-2.0,2.0);
    fHistDCAxy->GetXaxis()->SetTitle("DCAxy");
    fOutputListRec->Add(fHistDCAxy);

    fHistDCAz = new TH1F("fHistDCAz","DCAz",800,-5.0,5.0);
    fHistDCAz->GetXaxis()->SetTitle("DCAz");
    fOutputListRec->Add(fHistDCAz);

    fHistChi2perNDF = new TH1F("fHistChi2perNDF","Chi2 per NDF",100,0.0,10.0);
    fHistChi2perNDF->GetXaxis()->SetTitle("#chi^{2}/NDF");
    fOutputListRec->Add(fHistChi2perNDF);

    fHistTPCnSigmaEle = new TH2F("fHistTPCnSigmaEle", "TPC nSigma ele", 100,0.,10.,100,-5.,5.);
    fHistTPCnSigmaEle->GetXaxis()->SetTitle("#it{p}, GeV/#it{c}");
    fHistTPCnSigmaEle->GetYaxis()->SetTitle("TPCn#sigma_{ele}");
    fOutputListRec->Add(fHistTPCnSigmaEle);

    fHistTPCnSigmaPio= new TH2F("fHistTPCnSigmaPio", "TPC nSigma pion", 100,0.,10.,150,0.,15.);
    fHistTPCnSigmaPio->GetXaxis()->SetTitle("#it{p}, GeV/#it{c}");
    fHistTPCnSigmaPio->GetYaxis()->SetTitle("TPCn#sigma_{pio}");
    fOutputListRec->Add(fHistTPCnSigmaPio);

    fHistTOFnSigmaEle = new TH2F("fHistTOFnSigmaEle", "TOF nSigma ele", 100,0.,10.,100,-5.,5.);
    fHistTOFnSigmaEle->GetXaxis()->SetTitle("#it{p}, GeV/#it{c}");
    fHistTOFnSigmaEle->GetYaxis()->SetTitle("TOFn#sigma_{ele}");
    fOutputListRec->Add(fHistTOFnSigmaEle);

    fHistPairsRec = new TH2F("fHistPairsRec", "Rec e+e- from all sources", 50, 0.,5.0,10,0.,10.);
    fHistPairsRec->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsRec->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListRec->Add(fHistPairsRec);

    fHistPairsCCrec = new TH2F("fHistPairsCCrec", "Rec e+e- from cc", 500, 0.,5.0,100,0.,10.);
    fHistPairsCCrec->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsCCrec->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListRecCC->Add(fHistPairsCCrec);

    fHistPairsCCrec2 = new TH2F("fHistPairsCCrec2", "Rec e+e- from cc", 50, 0.,5.0,10,0.,10.);
    fHistPairsCCrec2->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsCCrec2->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListRecCC->Add(fHistPairsCCrec2);

    fHistPairsBBrec = new TH2F("fHistPairsBBrec", "Rec e+e- from cc", 500, 0.,5.0,100,0.,10.);
    fHistPairsBBrec->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsBBrec->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListRecBB->Add(fHistPairsBBrec);

    fHistPairsBBrecLS = new TH2F("fHistPairsBBrecLS", "Rec e+e- from cc", 500, 0.,5.0,100,0.,10.);
    fHistPairsBBrecLS->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsBBrecLS->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListRecBB->Add(fHistPairsBBrecLS);

    fHistPairsBBrec2 = new TH2F("fHistPairsBBrec2", "Rec e+e- from cc", 50, 0.,5.0,10,0.,10.);
    fHistPairsBBrec2->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsBBrec2->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListRecBB->Add(fHistPairsBBrec2);

    fHistPairsBBrecLS2 = new TH2F("fHistPairsBBrecLS2", "Rec e+e- from cc", 50, 0.,5.0,10,0.,10.);
    fHistPairsBBrecLS2->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistPairsBBrecLS2->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fOutputListRecBB->Add(fHistPairsBBrecLS2);

    TPCnSigmaEle_Eta_P_lin = new TH3F("TPCnSigmaEle_Eta_P_lin","TPC number of sigmas Electrons vs Eta and P", 200, -1., 1., 200, -10., 10., 100, 0., 5.);
    TOFnSigmaEle_Eta_P_lin = new TH3F("TOFnSigmaEle_Eta_P_lin","TOF number of sigmas Electrons vs Eta and P", 200, -1., 1., 200, -10., 10., 100, 0., 5.);
    fOutputListRec->Add(TPCnSigmaEle_Eta_P_lin);
    fOutputListRec->Add(TOFnSigmaEle_Eta_P_lin);

    TPCnSigmaEle_Eta_P_lin2 = new TH2F("TPCnSigmaEle_Eta_P_lin2", "TPC nSigma ele hadron rej", 100,0.,10.,100,-5.,5.);
    TPCnSigmaEle_Eta_P_lin2->GetXaxis()->SetTitle("#it{p}, GeV/#it{c}");
    TPCnSigmaEle_Eta_P_lin2->GetYaxis()->SetTitle("TPCn#sigma_{ele}");
    fOutputListRec->Add(TPCnSigmaEle_Eta_P_lin2);

    TOFnSigmaEle_Eta_P_lin2 = new TH2F("TOFnSigmaEle_Eta_P_lin2", "TOF nSigma ele hadron rej", 100,0.,10.,100,-5.,5.);
    TOFnSigmaEle_Eta_P_lin2->GetXaxis()->SetTitle("#it{p}, GeV/#it{c}");
    TOFnSigmaEle_Eta_P_lin2->GetYaxis()->SetTitle("TOFn#sigma_{ele}");
    fOutputListRec->Add(TOFnSigmaEle_Eta_P_lin2);

	fHistPairsGenPt = new TH1F("fHistPairsGenPt", "Gen e+e- Pt", 100,0.,10.);
	fHistPairsGenPt->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGen->Add(fHistPairsGenPt);

	fHistPairsRecPt = new TH1F("fHistPairsRecPt", "Rec e+e- Pt", 100,0.,10.);
	fHistPairsRecPt->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRec->Add(fHistPairsRecPt);

	fHistPairsGenMass = new TH1F("fHistPairsGenMass", "Gen e+e- Mass", 500,0.,5.);
	fHistPairsGenMass->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fOutputListGen->Add(fHistPairsGenMass);

	fHistPairsRecMass = new TH1F("fHistPairsRecMass", "Rec e+e- Mass", 500,0.,5.);
	fHistPairsRecMass->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fOutputListRec->Add(fHistPairsRecMass);

	fHistPairsGenDphi = new TH1F("fHistPairsGenDphi", "Gen e+e- DeltaPhi", 50,0.,TMath::Pi());
	fHistPairsGenDphi->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListGen->Add(fHistPairsGenDphi);

	fHistPairsRecDphi = new TH1F("fHistPairsRecDphi", "Rec e+e- DeltaPhi", 50,0.,TMath::Pi());
	fHistPairsRecDphi->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListRec->Add(fHistPairsRecDphi);

	fHistPairsGenMassPt = new TH2F("fHistPairsGenMassPt", "Gen e+e- Mass Pt", 500, 0.,5.0,100,0.,10.);
	fHistPairsGenMassPt->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsGenMassPt->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGen->Add(fHistPairsGenMassPt);

	fHistPairsRecMassPt = new TH2F("fHistPairsRecMassPt", "Rec e+e- Mass Pt", 500, 0.,5.0,100,0.,10.);
	fHistPairsRecMassPt->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsRecMassPt->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRec->Add(fHistPairsRecMassPt);

	fHistPairsGenPhiPt = new TH2F("fHistPairsGenPhiPt", "Gen e+e- dPhi Pt", 50,0.,TMath::Pi(),100,0.,10.);
	fHistPairsGenPhiPt->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsGenPhiPt->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGen->Add(fHistPairsGenPhiPt);

	fHistPairsRecPhiPt = new TH2F("fHistPairsRecPhiPt", "Rec e+e- dPhi Pt", 50,0.,TMath::Pi(),100,0.,10.);
	fHistPairsRecPhiPt->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsRecPhiPt->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRec->Add(fHistPairsRecPhiPt);

	fHistPairsGenMassPt2 = new TH2F("fHistPairsGenMassPt2", "Gen e+e- Mass Pt", 50, 0.,5.0,10,0.,10.);
	fHistPairsGenMassPt2->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsGenMassPt2->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGen->Add(fHistPairsGenMassPt2);

	fHistPairsRecMassPt2 = new TH2F("fHistPairsRecMassPt2", "Rec e+e- Mass Pt", 50, 0.,5.0,10,0.,10.);
	fHistPairsRecMassPt2->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsRecMassPt2->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRec->Add(fHistPairsRecMassPt2);

	fHistPairsGenPhiPt2 = new TH2F("fHistPairsGenPhiPt2", "Gen e+e- dPhi Pt", 9,0.,TMath::Pi(),10,0.,10.);
	fHistPairsGenPhiPt2->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsGenPhiPt2->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGen->Add(fHistPairsGenPhiPt2);

	fHistPairsRecPhiPt2 = new TH2F("fHistPairsRecPhiPt2", "Rec e+e- dPhi Pt", 9,0.,TMath::Pi(),10,0.,10.);
	fHistPairsRecPhiPt2->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsRecPhiPt2->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRec->Add(fHistPairsRecPhiPt2);

	fHistPairsGenPtMasPhi = new TH3F("fHistPairsGenPtMasPhi", "Gen e+e- Pt Mass Phi", 10,0.,10., 50, 0.,5.0, 9,0.,TMath::Pi());
	fHistPairsGenPtMasPhi->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fHistPairsGenPtMasPhi->GetYaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");;
	fHistPairsGenPtMasPhi->GetZaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListGen->Add(fHistPairsGenPtMasPhi);

	fHistPairsRecPtMasPhi = new TH3F("fHistPairsRecPtMasPhi", "Rec e+e- Pt Mass Phi", 10,0.,10., 50, 0.,5.0, 9,0.,TMath::Pi());
	fHistPairsRecPtMasPhi->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fHistPairsRecPtMasPhi->GetYaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");;
	fHistPairsRecPtMasPhi->GetZaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListRec->Add(fHistPairsRecPtMasPhi);

	fHistPairsGenPtCC = new TH1F("fHistPairsGenPtCC", "Gen e+e- Pt", 100,0.,10.);
	fHistPairsGenPtCC->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenCC->Add(fHistPairsGenPtCC);

	fHistPairsRecPtCC = new TH1F("fHistPairsRecPtCC", "Rec e+e- Pt", 100,0.,10.);
	fHistPairsRecPtCC->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecCC->Add(fHistPairsRecPtCC);

	fHistPairsGenMassCC = new TH1F("fHistPairsGenMassCC", "Gen e+e- Mass", 500,0.,5.);
	fHistPairsGenMassCC->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fOutputListGenCC->Add(fHistPairsGenMassCC);

	fHistPairsRecMassCC = new TH1F("fHistPairsRecMassCC", "Rec e+e- Mass", 500,0.,5.);
	fHistPairsRecMassCC->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fOutputListRecCC->Add(fHistPairsRecMassCC);

	fHistPairsGenDphiCC = new TH1F("fHistPairsGenDphiCC", "Gen e+e- DeltaPhi", 50,0.,TMath::Pi());
	fHistPairsGenDphiCC->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListGenCC->Add(fHistPairsGenDphiCC);

	fHistPairsRecDphiCC = new TH1F("fHistPairsRecDphiCC", "Rec e+e- DeltaPhi", 50,0.,TMath::Pi());
	fHistPairsRecDphiCC->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListRecCC->Add(fHistPairsRecDphiCC);

	fHistPairsGenMassPtCC = new TH2F("fHistPairsGenMassPtCC", "Gen e+e- Mass Pt", 500, 0.,5.0,100,0.,10.);
	fHistPairsGenMassPtCC->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsGenMassPtCC->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenCC->Add(fHistPairsGenMassPtCC);

	fHistPairsRecMassPtCC = new TH2F("fHistPairsRecMassPtCC", "Rec e+e- Mass Pt", 500, 0.,5.0,100,0.,10.);
	fHistPairsRecMassPtCC->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsRecMassPtCC->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecCC->Add(fHistPairsRecMassPtCC);

	fHistPairsGenPhiPtCC = new TH2F("fHistPairsGenPhiPtCC", "Gen e+e- dPhi Pt", 50,0.,TMath::Pi(),100,0.,10.);
	fHistPairsGenPhiPtCC->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsGenPhiPtCC->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenCC->Add(fHistPairsGenPhiPtCC);

	fHistPairsRecPhiPtCC = new TH2F("fHistPairsRecPhiPtCC", "Rec e+e- dPhi Pt", 50,0.,TMath::Pi(),100,0.,10.);
	fHistPairsRecPhiPtCC->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsRecPhiPtCC->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecCC->Add(fHistPairsRecPhiPtCC);

	fHistPairsGenMassPt2CC = new TH2F("fHistPairsGenMassPt2CC", "Gen e+e- Mass Pt", 50, 0.,5.0,10,0.,10.);
	fHistPairsGenMassPt2CC->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsGenMassPt2CC->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenCC->Add(fHistPairsGenMassPt2CC);

	fHistPairsRecMassPt2CC = new TH2F("fHistPairsRecMassPt2CC", "Rec e+e- Mass Pt", 50, 0.,5.0,10,0.,10.);
	fHistPairsRecMassPt2CC->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsRecMassPt2CC->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecCC->Add(fHistPairsRecMassPt2CC);

	fHistPairsGenPhiPt2CC = new TH2F("fHistPairsGenPhiPt2CC", "Gen e+e- dPhi Pt", 9,0.,TMath::Pi(),10,0.,10.);
	fHistPairsGenPhiPt2CC->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsGenPhiPt2CC->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenCC->Add(fHistPairsGenPhiPt2CC);

	fHistPairsRecPhiPt2CC = new TH2F("fHistPairsRecPhiPt2CC", "Rec e+e- dPhi Pt", 9,0.,TMath::Pi(),10,0.,10.);
	fHistPairsRecPhiPt2CC->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsRecPhiPt2CC->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecCC->Add(fHistPairsRecPhiPt2CC);

	fHistPairsGenPtMasPhiCC = new TH3F("fHistPairsGenPtMasPhiCC", "Gen e+e- Pt Mass Phi", 10,0.,10., 50, 0.,5.0, 9,0.,TMath::Pi());
	fHistPairsGenPtMasPhiCC->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fHistPairsGenPtMasPhiCC->GetYaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");;
	fHistPairsGenPtMasPhiCC->GetZaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListGenCC->Add(fHistPairsGenPtMasPhiCC);

	fHistPairsRecPtMasPhiCC = new TH3F("fHistPairsRecPtMasPhiCC", "Rec e+e- Pt Mass Phi", 10,0.,10., 50, 0.,5.0, 9,0.,TMath::Pi());
	fHistPairsRecPtMasPhiCC->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fHistPairsRecPtMasPhiCC->GetYaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");;
	fHistPairsRecPtMasPhiCC->GetZaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListRecCC->Add(fHistPairsRecPtMasPhiCC);

	fHistPairsGenPtBB = new TH1F("fHistPairsGenPtBB", "Gen e+e- Pt", 100,0.,10.);
	fHistPairsGenPtBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenBB->Add(fHistPairsGenPtBB);

	fHistPairsGenPtLSBB = new TH1F("fHistPairsGenPtLSBB", "Gen e+e- Pt", 100,0.,10.);
	fHistPairsGenPtLSBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenBB->Add(fHistPairsGenPtLSBB);

	fHistPairsRecPtBB = new TH1F("fHistPairsRecPtBB", "Rec e+e- Pt", 100,0.,10.);
	fHistPairsRecPtBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecBB->Add(fHistPairsRecPtBB);

	fHistPairsRecPtLSBB = new TH1F("fHistPairsRecPtLSBB", "Rec e+e- Pt", 100,0.,10.);
	fHistPairsRecPtLSBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecBB->Add(fHistPairsRecPtLSBB);

	fHistPairsGenMassBB = new TH1F("fHistPairsGenMassBB", "Gen e+e- Mass", 500,0.,5.);
	fHistPairsGenMassBB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fOutputListGenBB->Add(fHistPairsGenMassBB);

	fHistPairsGenMassLSBB = new TH1F("fHistPairsGenMassLSBB", "Gen e+e- Mass", 500,0.,5.);
	fHistPairsGenMassLSBB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fOutputListGenBB->Add(fHistPairsGenMassLSBB);

	fHistPairsRecMassBB = new TH1F("fHistPairsRecMassBB", "Rec e+e- Mass", 500,0.,5.);
	fHistPairsRecMassBB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fOutputListRecBB->Add(fHistPairsRecMassBB);

	fHistPairsRecMassLSBB = new TH1F("fHistPairsRecMassLSBB", "Rec e+e- Mass", 500,0.,5.);
	fHistPairsRecMassLSBB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fOutputListRecBB->Add(fHistPairsRecMassLSBB);

	fHistPairsGenDphiBB = new TH1F("fHistPairsGenDphiBB", "Gen e+e- DeltaPhi", 50,0.,TMath::Pi());
	fHistPairsGenDphiBB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListGenBB->Add(fHistPairsGenDphiBB);

	fHistPairsGenDphiLSBB = new TH1F("fHistPairsGenDphiLSBB", "Gen e+e- DeltaPhi", 50,0.,TMath::Pi());
	fHistPairsGenDphiLSBB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListGenBB->Add(fHistPairsGenDphiLSBB);

	fHistPairsRecDphiBB = new TH1F("fHistPairsRecDphiBB", "Rec e+e- DeltaPhi", 50,0.,TMath::Pi());
	fHistPairsRecDphiBB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListRecBB->Add(fHistPairsRecDphiBB);

	fHistPairsRecDphiLSBB = new TH1F("fHistPairsRecDphiLSBB", "Rec e+e- DeltaPhi", 50,0.,TMath::Pi());
	fHistPairsRecDphiLSBB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListRecBB->Add(fHistPairsRecDphiLSBB);

	fHistPairsGenMassPtBB = new TH2F("fHistPairsGenMassPtBB", "Gen e+e- Mass Pt", 500, 0.,5.0,100,0.,10.);
	fHistPairsGenMassPtBB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsGenMassPtBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenBB->Add(fHistPairsGenMassPtBB);

	fHistPairsGenMassPtLSBB = new TH2F("fHistPairsGenMassPtLSBB", "Gen e+e- Mass Pt", 500, 0.,5.0,100,0.,10.);
	fHistPairsGenMassPtLSBB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsGenMassPtLSBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenBB->Add(fHistPairsGenMassPtLSBB);

	fHistPairsRecMassPtBB = new TH2F("fHistPairsRecMassPtBB", "Rec e+e- Mass Pt", 500, 0.,5.0,100,0.,10.);
	fHistPairsRecMassPtBB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsRecMassPtBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecBB->Add(fHistPairsRecMassPtBB);

	fHistPairsRecMassPtLSBB = new TH2F("fHistPairsRecMassPtLSBB", "Rec e+e- Mass Pt", 500, 0.,5.0,100,0.,10.);
	fHistPairsRecMassPtLSBB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsRecMassPtLSBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecBB->Add(fHistPairsRecMassPtLSBB);

	fHistPairsGenPhiPtBB = new TH2F("fHistPairsGenPhiPtBB", "Gen e+e- dPhi Pt", 50,0.,TMath::Pi(),100,0.,10.);
	fHistPairsGenPhiPtBB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsGenPhiPtBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenBB->Add(fHistPairsGenPhiPtBB);

	fHistPairsGenPhiPtLSBB = new TH2F("fHistPairsGenPhiPtLSBB", "Gen e+e- dPhi Pt", 50,0.,TMath::Pi(),100,0.,10.);
	fHistPairsGenPhiPtLSBB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsGenPhiPtLSBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenBB->Add(fHistPairsGenPhiPtLSBB);

	fHistPairsRecPhiPtBB = new TH2F("fHistPairsRecPhiPtBB", "Rec e+e- dPhi Pt", 50,0.,TMath::Pi(),100,0.,10.);
	fHistPairsRecPhiPtBB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsRecPhiPtBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecBB->Add(fHistPairsRecPhiPtBB);

	fHistPairsRecPhiPtLSBB = new TH2F("fHistPairsRecPhiPtLSBB", "Rec e+e- dPhi Pt", 50,0.,TMath::Pi(),100,0.,10.);
	fHistPairsRecPhiPtLSBB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsRecPhiPtLSBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecBB->Add(fHistPairsRecPhiPtLSBB);

	fHistPairsGenMassPt2BB = new TH2F("fHistPairsGenMassPt2BB", "Gen e+e- Mass Pt", 50, 0.,5.0,10,0.,10.);
	fHistPairsGenMassPt2BB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsGenMassPt2BB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenBB->Add(fHistPairsGenMassPt2BB);

	fHistPairsGenMassPt2LSBB = new TH2F("fHistPairsGenMassPt2LSBB", "Gen e+e- Mass Pt", 50, 0.,5.0,10,0.,10.);
	fHistPairsGenMassPt2LSBB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsGenMassPt2LSBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenBB->Add(fHistPairsGenMassPt2LSBB);

	fHistPairsRecMassPt2BB = new TH2F("fHistPairsRecMassPt2BB", "Rec e+e- Mass Pt", 50, 0.,5.0,10,0.,10.);
	fHistPairsRecMassPt2BB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsRecMassPt2BB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecBB->Add(fHistPairsRecMassPt2BB);

	fHistPairsRecMassPt2LSBB = new TH2F("fHistPairsRecMassPt2LSBB", "Rec e+e- Mass Pt", 50, 0.,5.0,10,0.,10.);
	fHistPairsRecMassPt2LSBB->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
	fHistPairsRecMassPt2LSBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecBB->Add(fHistPairsRecMassPt2LSBB);

	fHistPairsGenPhiPt2BB = new TH2F("fHistPairsGenPhiPt2BB", "Gen e+e- dPhi Pt", 9,0.,TMath::Pi(),10,0.,10.);
	fHistPairsGenPhiPt2BB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsGenPhiPt2BB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenBB->Add(fHistPairsGenPhiPt2BB);

	fHistPairsGenPhiPt2LSBB = new TH2F("fHistPairsGenPhiPt2LSBB", "Gen e+e- dPhi Pt", 9,0.,TMath::Pi(),10,0.,10.);
	fHistPairsGenPhiPt2LSBB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsGenPhiPt2LSBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListGenBB->Add(fHistPairsGenPhiPt2LSBB);

	fHistPairsRecPhiPt2BB = new TH2F("fHistPairsRecPhiPt2BB", "Rec e+e- dPhi Pt", 9,0.,TMath::Pi(),10,0.,10.);
	fHistPairsRecPhiPt2BB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsRecPhiPt2BB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecBB->Add(fHistPairsRecPhiPt2BB);

	fHistPairsRecPhiPt2LSBB = new TH2F("fHistPairsRecPhiPt2LSBB", "Rec e+e- dPhi Pt", 9,0.,TMath::Pi(),10,0.,10.);
	fHistPairsRecPhiPt2LSBB->GetXaxis()->SetTitle("#Delta#varphi, rad");
	fHistPairsRecPhiPt2LSBB->GetYaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fOutputListRecBB->Add(fHistPairsRecPhiPt2LSBB);

	fHistPairsGenPtMasPhiBB = new TH3F("fHistPairsGenPtMasPhiBB", "Gen e+e- Pt Mass Phi", 10,0.,10., 50, 0.,5.0, 9,0.,TMath::Pi());
	fHistPairsGenPtMasPhiBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fHistPairsGenPtMasPhiBB->GetYaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");;
	fHistPairsGenPtMasPhiBB->GetZaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListGenBB->Add(fHistPairsGenPtMasPhiBB);

	fHistPairsGenPtMasPhiLSBB = new TH3F("fHistPairsGenPtMasPhiLSBB", "Gen e+e- Pt Mass Phi", 10,0.,10., 50, 0.,5.0, 9,0.,TMath::Pi());
	fHistPairsGenPtMasPhiLSBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fHistPairsGenPtMasPhiLSBB->GetYaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");;
	fHistPairsGenPtMasPhiLSBB->GetZaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListGenBB->Add(fHistPairsGenPtMasPhiLSBB);

	fHistPairsRecPtMasPhiBB = new TH3F("fHistPairsRecPtMasPhiBB", "Rec e+e- Pt Mass Phi", 10,0.,10., 50, 0.,5.0, 9,0.,TMath::Pi());
	fHistPairsRecPtMasPhiBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fHistPairsRecPtMasPhiBB->GetYaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");;
	fHistPairsRecPtMasPhiBB->GetZaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListRecBB->Add(fHistPairsRecPtMasPhiBB);

	fHistPairsRecPtMasPhiLSBB = new TH3F("fHistPairsRecPtMasPhiLSBB", "Rec e+e- Pt Mass Phi", 10,0.,10., 50, 0.,5.0, 9,0.,TMath::Pi());
	fHistPairsRecPtMasPhiLSBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
	fHistPairsRecPtMasPhiLSBB->GetYaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");;
	fHistPairsRecPtMasPhiLSBB->GetZaxis()->SetTitle("#Delta#varphi, rad");
	fOutputListRecBB->Add(fHistPairsRecPtMasPhiLSBB);

	fHistPairsGenDCAxy = new TH1F("fHistPairsGenDCAxy","Rec dielectrons DCAxy",400,0.0,2.0);
    fHistPairsGenDCAxy->GetXaxis()->SetTitle("DCAxy");
    fOutputListGen->Add(fHistPairsGenDCAxy);
    
    fHistPairsGenDCAxySig = new TH1F("fHistPairsGenDCAxySig","Rec dielectrons DCAxy",400,0.0,20.0);
    fHistPairsGenDCAxySig->GetXaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListGen->Add(fHistPairsGenDCAxySig);
    
    fHistPairsGenPtDCAxySig = new TH2F("fHistPairsGenPtDCAxySig","Rec dielectrons DCAxy",10,0.,10.,400,0.0,20.0);
	fHistPairsGenPtDCAxySig->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fHistPairsGenPtDCAxySig->GetYaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListGen->Add(fHistPairsGenPtDCAxySig);

    fHistPairsGenDCAxyCC = new TH1F("fHistPairsGenDCAxyCC","Rec dielectrons DCAxy",400,0.0,1.0);
    fHistPairsGenDCAxyCC->GetXaxis()->SetTitle("DCAxy");
    fOutputListGenCC->Add(fHistPairsGenDCAxyCC);
    
    fHistPairsGenDCAxySigCC = new TH1F("fHistPairsGenDCAxySigCC","Rec dielectrons DCAxy",400,0.0,20.0);
    fHistPairsGenDCAxySigCC->GetXaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListGenCC->Add(fHistPairsGenDCAxySigCC);
    
    fHistPairsGenPtDCAxySigCC = new TH2F("fHistPairsGenPtDCAxySigCC","Rec dielectrons DCAxy",10,0.,10.,400,0.0,20.0);
	fHistPairsGenPtDCAxySigCC->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fHistPairsGenPtDCAxySigCC->GetYaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListGenCC->Add(fHistPairsGenPtDCAxySigCC);

    fHistPairsGenDCAxyBB = new TH1F("fHistPairsGenDCAxyBB","Rec dielectrons DCAxy",400,0.0,1.0);
    fHistPairsGenDCAxyBB->GetXaxis()->SetTitle("DCAxy");
    fOutputListGenBB->Add(fHistPairsGenDCAxyBB);
    
    fHistPairsGenDCAxySigBB = new TH1F("fHistPairsGenDCAxySigBB","Rec dielectrons DCAxy",400,0.0,20.0);
    fHistPairsGenDCAxySigBB->GetXaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListGenBB->Add(fHistPairsGenDCAxySigBB);
    
    fHistPairsGenPtDCAxySigBB = new TH2F("fHistPairsGenPtDCAxySigBB","Rec dielectrons DCAxy",10,0.,10.,400,0.0,20.0);
	fHistPairsGenPtDCAxySigBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fHistPairsGenPtDCAxySigBB->GetYaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListGenBB->Add(fHistPairsGenPtDCAxySigBB);

    fHistPairsGenDCAxyLSBB = new TH1F("fHistPairsGenDCAxyLSBB","Rec dielectrons DCAxy",400,0.0,1.0);
    fHistPairsGenDCAxyLSBB->GetXaxis()->SetTitle("DCAxy");
    fOutputListGenBB->Add(fHistPairsGenDCAxyLSBB);
    
    fHistPairsGenDCAxySigLSBB = new TH1F("fHistPairsGenDCAxySigLSBB","Rec dielectrons DCAxy",400,0.0,20.0);
    fHistPairsGenDCAxySigLSBB->GetXaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListGenBB->Add(fHistPairsGenDCAxySigLSBB);
    
    fHistPairsGenPtDCAxySigLSBB = new TH2F("fHistPairsGenPtDCAxySigLSBB","Rec dielectrons DCAxy",10,0.,10.,400,0.0,20.0);
	fHistPairsGenPtDCAxySigLSBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fHistPairsGenPtDCAxySigLSBB->GetYaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListGenBB->Add(fHistPairsGenPtDCAxySigLSBB);

	fHistPairsDCAxy = new TH1F("fHistPairsDCAxy","Rec dielectrons DCAxy",400,0.0,2.0);
    fHistPairsDCAxy->GetXaxis()->SetTitle("DCAxy");
    fOutputListRec->Add(fHistPairsDCAxy);
    
    fHistPairsDCAxySig = new TH1F("fHistPairsDCAxySig","Rec dielectrons DCAxy",400,0.0,20.0);
    fHistPairsDCAxySig->GetXaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListRec->Add(fHistPairsDCAxySig);
    
    fHistPairsPtDCAxySig = new TH2F("fHistPairsPtDCAxySig","Rec dielectrons DCAxy",10,0.,10.,400,0.0,20.0);
	fHistPairsPtDCAxySig->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fHistPairsPtDCAxySig->GetYaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListRec->Add(fHistPairsPtDCAxySig);

    fHistPairsDCAxyCC = new TH1F("fHistPairsDCAxyCC","Rec dielectrons DCAxy",400,0.0,1.0);
    fHistPairsDCAxyCC->GetXaxis()->SetTitle("DCAxy");
    fOutputListRecCC->Add(fHistPairsDCAxyCC);
    
    fHistPairsDCAxySigCC = new TH1F("fHistPairsDCAxySigCC","Rec dielectrons DCAxy",400,0.0,20.0);
    fHistPairsDCAxySigCC->GetXaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListRecCC->Add(fHistPairsDCAxySigCC);
    
    fHistPairsPtDCAxySigCC = new TH2F("fHistPairsPtDCAxySigCC","Rec dielectrons DCAxy",10,0.,10.,400,0.0,20.0);
	fHistPairsPtDCAxySigCC->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fHistPairsPtDCAxySigCC->GetYaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListRecCC->Add(fHistPairsPtDCAxySigCC);

    fHistPairsDCAxyBB = new TH1F("fHistPairsDCAxyBB","Rec dielectrons DCAxy",400,0.0,1.0);
    fHistPairsDCAxyBB->GetXaxis()->SetTitle("DCAxy");
    fOutputListRecBB->Add(fHistPairsDCAxyBB);
    
    fHistPairsDCAxySigBB = new TH1F("fHistPairsDCAxySigBB","Rec dielectrons DCAxy",400,0.0,20.0);
    fHistPairsDCAxySigBB->GetXaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListRecBB->Add(fHistPairsDCAxySigBB);
    
    fHistPairsPtDCAxySigBB = new TH2F("fHistPairsPtDCAxySigBB","Rec dielectrons DCAxy",10,0.,10.,400,0.0,20.0);
	fHistPairsPtDCAxySigBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fHistPairsPtDCAxySigBB->GetYaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListRecBB->Add(fHistPairsPtDCAxySigBB);

    fHistPairsDCAxyLSBB = new TH1F("fHistPairsDCAxyLSBB","Rec dielectrons DCAxy",400,0.0,1.0);
    fHistPairsDCAxyLSBB->GetXaxis()->SetTitle("DCAxy");
    fOutputListRecBB->Add(fHistPairsDCAxyLSBB);
    
    fHistPairsDCAxySigLSBB = new TH1F("fHistPairsDCAxySigLSBB","Rec dielectrons DCAxy",400,0.0,20.0);
    fHistPairsDCAxySigLSBB->GetXaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListRecBB->Add(fHistPairsDCAxySigLSBB);
    
    fHistPairsPtDCAxySigLSBB = new TH2F("fHistPairsPtDCAxySigLSBB","Rec dielectrons DCAxy",10,0.,10.,400,0.0,20.0);
	fHistPairsPtDCAxySigLSBB->GetXaxis()->SetTitle("#it{p}_{T,ee}, GeV/#it{c}");
    fHistPairsPtDCAxySigLSBB->GetYaxis()->SetTitle("DCAxy (#sigma)");
    fOutputListRecBB->Add(fHistPairsPtDCAxySigLSBB);

    fPRec_Gen = new TH2F("fPRec_Gen", "Rec - Gen e+- P", 1500,0.,15.,1000,-5,5);
  	fPRec_Gen->GetXaxis()->SetTitle("#it{p}^{gen} (GeV/#it{c})");
  	fPRec_Gen->GetYaxis()->SetTitle("#it{p}^{rec} - #it{p}^{gen} (GeV/#it{c})");
  	fOutputListRec->Add(fPRec_Gen);
    fPtRec_Gen = new TH2F("fPtRec_Gen", "Rec - Gen e+- Pt", 1500,0.,15.,1000,-5.,5.);
  	fPtRec_Gen->GetXaxis()->SetTitle("#it{p}_{T,ee}^{gen} (GeV/#it{c})");
  	fPtRec_Gen->GetYaxis()->SetTitle("#it{p}_{T,e}^{rec} - #it{p}_{T,e}^{gen} (GeV/#it{c})");
  	fOutputListRec->Add(fPtRec_Gen);
    fPtRecOverGen = new TH2F("fPtRecOverGen", "Rec / Gen e+- Pt", 1500,0.,15.,400,0.,1.2);
  	fPtRecOverGen->GetXaxis()->SetTitle("#it{p}_{T,ee}^{gen} (GeV/#it{c})");
  	fPtRecOverGen->GetYaxis()->SetTitle("#it{p}_{T,e}^{rec} - #it{p}_{T,e}^{gen} (GeV/#it{c})");
  	fOutputListRec->Add(fPtRecOverGen);
    fEtaRec_Gen = new TH2F("fEtaRec_Gen", "Rec - Gen e+- Eta", 1500,0.,15.,600,-0.3,0.3);
  	fEtaRec_Gen->GetXaxis()->SetTitle("#eta^{gen}");
  	fEtaRec_Gen->GetYaxis()->SetTitle("#eta^{rec} - #eta^{gen}");
  	fOutputListRec->Add(fEtaRec_Gen);
    fElePhiRec_Gen = new TH2F("fElePhiRec_Gen", "Rec - Gen e- Eta", 1500,0.,15.,400,-0.2,0.2);
  	fElePhiRec_Gen->GetXaxis()->SetTitle("#eta^{gen}");
  	fElePhiRec_Gen->GetYaxis()->SetTitle("#eta^{rec} - #eta^{gen}");
  	fOutputListRec->Add(fElePhiRec_Gen);
    fPosPhiRec_Gen = new TH2F("fPosPhiRec_Gen", "Rec - Gen e+ Eta", 1500,0.,15.,400,-0.2,0.2);
  	fPosPhiRec_Gen->GetXaxis()->SetTitle("#eta^{gen}");
  	fPosPhiRec_Gen->GetYaxis()->SetTitle("#eta^{rec} - #eta^{gen}");
  	fOutputListRec->Add(fPosPhiRec_Gen);
  	
  	fPairDCARec_Gen = new TH2F("fPairDCARec_Gen", "Rec - Gen e+ DCA", 1500,0.,15.,800,-40.,40.);
  	fPairDCARec_Gen->GetXaxis()->SetTitle("#it{p}_{T,ee}^{gen} (GeV/#it{c})");
  	fPairDCARec_Gen->GetYaxis()->SetTitle("DCA_{xy}^{rec} - DCA_{xy}^{gen}");
  	fOutputListRec->Add(fPairDCARec_Gen);
  	fEleDCARes = new TH2F("fEleDCARes", "e- DCA res", 500,0.,15.,5000,0.,0.001);
  	fEleDCARes->GetXaxis()->SetTitle("#it{p}_{T,ee}^{gen} (GeV/#it{c})");
  	fEleDCARes->GetYaxis()->SetTitle("#sigma{DCA_{xy}}");
  	fOutputListRec->Add(fEleDCARes);
  	fPosDCARes = new TH2F("fPosDCARes", "e+ DCA res", 500,0.,15.,5000,0.,0.001);
  	fPosDCARes->GetXaxis()->SetTitle("#it{p}_{T,ee}^{gen} (GeV/#it{c})");
  	fPosDCARes->GetYaxis()->SetTitle("#sigma{DCA_{xy}}");
  	fOutputListRec->Add(fPosDCARes);
  	fEleDCARes2 = new TH2F("fEleDCARes2", "e- DCA res", 500,0.,15.,1000,0.,0.0004);
  	fEleDCARes2->GetXaxis()->SetTitle("#it{p}_{T,ee}^{gen} (GeV/#it{c})");
  	fEleDCARes2->GetYaxis()->SetTitle("#sigma{DCA_{xy}}");
  	fOutputListRec->Add(fEleDCARes2);
  	fPosDCARes2 = new TH2F("fPosDCARes2", "e+ DCA res", 500,0.,15.,1000,0.,0.0004);
  	fPosDCARes2->GetXaxis()->SetTitle("#it{p}_{T,ee}^{gen} (GeV/#it{c})");
  	fPosDCARes2->GetYaxis()->SetTitle("#sigma{DCA_{xy}}");
  	fOutputListRec->Add(fPosDCARes2);
  	fEleDCARes3 = new TH2F("fEleDCARes3", "e- DCA res", 500,0.,15.,1000,0.,0.00004);
  	fEleDCARes3->GetXaxis()->SetTitle("#it{p}_{T,ee}^{gen} (GeV/#it{c})");
  	fEleDCARes3->GetYaxis()->SetTitle("#sigma{DCA_{xy}}");
  	fOutputListRec->Add(fEleDCARes3);
  	fPosDCARes3 = new TH2F("fPosDCARes3", "e+ DCA res", 500,0.,15.,1000,0.,0.00004);
  	fPosDCARes3->GetXaxis()->SetTitle("#it{p}_{T,ee}^{gen} (GeV/#it{c})");
  	fPosDCARes3->GetYaxis()->SetTitle("#sigma{DCA_{xy}}");
  	fOutputListRec->Add(fPosDCARes3);

    fOutputList->Add(fOutputListGen);
    fOutputList->Add(fOutputListRec);
    fOutputList->Add(fOutputListGenCC);
    fOutputList->Add(fOutputListRecCC);
    fOutputList->Add(fOutputListGenBB);
    fOutputList->Add(fOutputListRecBB);

//;Eta;n#sigma_{ele}^{TPC};#it{p}(GeV/#it{c})
    PostData(1, fOutputList); // Postdata will notify the analysis manager of changes/updates to the
    // fOutputList object. the manager will in the end take care of writing your output to file
    // so it needs to know what's in the output
}

//________________________________________________________________________
void AliAnalysisTaskeeCor::UserExec(Option_t *)
{
    // Main loop
    // Called for each event

    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    if (!man) {
        // Printf("no analysis manager!\n");
        return;
    }

    AliMCEventHandler *mcHandler = (AliMCEventHandler*)man->GetMCtruthEventHandler();
    if (!mcHandler) {
        // Printf("no MC handler!\n");
        return;
    }

    AliAODInputHandler* aodHandler=(AliAODInputHandler*)man->GetInputEventHandler();
    if (!aodHandler) {
        // Printf("ERROR: AOD handler not available\n");
        return;
    }

    mcEvent = aodHandler->MCEvent();
    if (!mcEvent) {
        // Printf("ERROR: mcEvent not available\n");
        return;
    }

    aodEvent = aodHandler->GetEvent();
    if (!aodEvent){
        Printf("ERROR: aodEvent not available\n");
        return;
    }

    //get generator header
    TList *l = (TList*)mcEvent->GetCocktailList();
    // Printf("GetCocktailList entries = %i", l->GetEntries());
    AliGenEventHeader* gh=(AliGenEventHeader*)l->At(0);
    TString genname=gh->GetName();

    // Printf("Genname: %s\n",genname.Data());
    Bool_t pythiaCCfile = kFALSE;
    Bool_t pythiaBBfile = kFALSE;
    Bool_t pythiaBfile = kFALSE;

    if(genname.Contains("Pythia CC")){
        fEventStat->Fill(kPythiaCC);
        pythiaCCfile = kTRUE;
    } else if(genname.Contains("Pythia BB")){
        fEventStat->Fill(kPythiaBB);
        pythiaBBfile = kTRUE;
    } else if(genname.Contains("Pythia B_")){
        fEventStat->Fill(kPythiaB);
        pythiaBfile = kTRUE;
    } else {
		fEventStat->Fill(kOther);
		PostData(1, fOutputList);
		return;
	}

    //if (!pythiaCCfile) return;
    // check the number of primary charmed hadrons produced in event
    Int_t number_of_cquarks = 0;
    Int_t number_of_cbarquarks = 0;

    Int_t number_of_charmedH = 0;
    Int_t number_of_anticharmedH = 0;

    Int_t number_of_bquarks = 0;
    Int_t number_of_bbarquarks = 0;

    Int_t nMCtracks = mcEvent->GetNumberOfTracks();
    for(Int_t iMCtrack = 0; iMCtrack < nMCtracks; iMCtrack++){
        if(mcEvent->GetTrack(iMCtrack)->PdgCode()== 4) number_of_cquarks++;
        if(mcEvent->GetTrack(iMCtrack)->PdgCode()==-4) number_of_cbarquarks++;

        if(mcEvent->GetTrack(iMCtrack)->PdgCode()== 5) number_of_bquarks++;
        if(mcEvent->GetTrack(iMCtrack)->PdgCode()==-5) number_of_bbarquarks++;
    }

    // remove events with bbbar quarks for charm file
    if (pythiaCCfile){
		if (number_of_bquarks != 0 || number_of_bbarquarks != 0){
			PostData(1, fOutputList);
			return;
		}
		// analyse events only with 1 charmed hadron-antihadron pair to have clear environment
		if (number_of_cquarks != 1 || number_of_cbarquarks != 1){
			PostData(1, fOutputList);
			return;
		}
	}// remove events with ccbar quarks for beauty file
	else { // should be pythiaBBfile || pythiaBfile?
		if (number_of_cquarks != 0 || number_of_cbarquarks != 0) {
			PostData(1, fOutputList);
			return;
		}
		// analyse events only with 1 beauty-antibeauty pair to have clear environment
		if (pythiaBBfile && (number_of_bquarks != 1 || number_of_bbarquarks != 1)){
			PostData(1, fOutputList);
			return;
		}
	}
	
	const AliVVertex *pVtxMC = mcEvent->GetPrimaryVertex();
	Double_t fPVxMC = pVtxMC->GetX();//cm
	Double_t fPVyMC = pVtxMC->GetY();//cm
	Double_t fPVzMC = pVtxMC->GetZ();//cm

    // check single tracks
    gRandom->SetSeed(0);
    std::vector<LMEEparticle> LMEEelectrons,LMEEpositrons,LMEEelectrons2,LMEEpositrons2;

    for(Int_t iMC = 0; iMC < nMCtracks; iMC++){
        AliMCParticle *part = (AliMCParticle*)(mcEvent->GetTrack(iMC));
        if(!part) continue;

		Double_t minPtCutBefSmr = fPtMinCut - 0.1;
		Double_t minEtaCutBefSmr = 1.1;
        if(part->Pt() < minPtCutBefSmr || TMath::Abs(part->Eta()) > minEtaCutBefSmr) continue;
		if(!part->IsPrimary()) continue;
        if(TMath::Abs(part->PdgCode()) != 11) continue;

        Int_t mLab = part->GetMother();
        AliMCParticle *mother(0x0);
        if(mLab >= 0) mother = (AliMCParticle*)(mcEvent->GetTrack(mLab));
        if(!mother) continue;
        Int_t mPDG = mother->PdgCode();
        Int_t grmLab = mother->GetMother();
        AliMCParticle *grmother(0x0);
        if(grmLab >= 0) grmother = (AliMCParticle*)(mcEvent->GetTrack(grmLab));
        Int_t grmPDG(0);
        if(grmother) grmPDG = grmother->PdgCode();

        LMEEparticle lmeeLeg;
        lmeeLeg.genP     = part->P();
        lmeeLeg.genPt    = part->Pt();
        lmeeLeg.genTheta = part->Theta();
        lmeeLeg.genEta   = part->Eta();
        lmeeLeg.genPhi   = part->Phi();
        lmeeLeg.charge   = part->Charge();
        lmeeLeg.label    = iMC;
        lmeeLeg.mlabel   = mLab;
        lmeeLeg.mPDG     = mPDG;
        lmeeLeg.grmlabel = grmLab;
        lmeeLeg.grmPDG   = grmPDG;

        if(fUseSmearing){
          if (fPtSmr)  lmeeLeg.genPt  = lmeeLeg.genPt*GetPtSmr(part->Pt());
          if (fEtaSmr) lmeeLeg.genEta = lmeeLeg.genEta - GetEtaSmr(part->Pt());
          if (fPhiEleSmr && fPhiPosSmr) lmeeLeg.genPhi = lmeeLeg.genPhi - GetPhiSmr(part->Pt(),lmeeLeg.charge);
        }
        if(lmeeLeg.genPt < fPtMinCut || TMath::Abs(lmeeLeg.genEta) > 0.8) continue; //Real cuts only after smearing

        lmeeLeg.isCharm  = IsCharmedEle(iMC);
        lmeeLeg.isBeauty = IsBeautyEle(iMC);
        
        Double_t partVx = part->Xv();//cm
        Double_t partVy = part->Yv();//cm
        Double_t q = lmeeLeg.charge/3;
        if (lmeeLeg.charge > 0) q = 1.;
        else q = -1.;
        Double_t Radius = lmeeLeg.genPt/(0.3*0.5*q);
        
        Double_t x = partVx + 100*Radius*TMath::Sin(lmeeLeg.genPhi);
        Double_t y = partVy - 100*Radius*TMath::Cos(lmeeLeg.genPhi);
        Double_t DCAxGen = x - fPVxMC;
        Double_t DCAyGen = y - fPVyMC;
        Double_t DCAxyGen = TMath::Sqrt(DCAxGen*DCAxGen + DCAyGen*DCAyGen) - 100*TMath::Abs(Radius);
        
        Double_t DCAxySigGen = DCAxyGen;
        Double_t pointRes = -1.0;
        Double_t rndm;
        if (fDCASmearingByMath){
			//Double_t theta = part->Theta();
			Double_t theta = (180./TMath::Pi())*2.*atan(exp(-lmeeLeg.genEta)); //EtaToTheta -> The smearing is applied in Eta, need to be propagated to theta
			Double_t sinTheta = TMath::Sin(theta);
			Double_t absSinTheta = TMath::Abs(sinTheta);
			Double_t sqrtSinTheta = TMath::Sqrt(absSinTheta);
			Double_t A = 10.;
			Double_t B = 53./(lmeeLeg.genPt*sqrtSinTheta);
			Double_t trackRes = (A + B)/10000; //cm
			Double_t C = 6.;
			Double_t D = 63./TMath::Sqrt(lmeeLeg.genPt);
			Double_t E = -8./lmeeLeg.genPt;
			Double_t vtxRes = (C + D + E)/10000; //cm OLD
			pointRes = trackRes + vtxRes;
		}
		else if (fDCASmearingByMaps){
			pointRes = GetDCASmr(lmeeLeg.genPt); //cm
			if (fDCASmrFromMC) pointRes*100;
		}
		else if (fDCASmearingByPars){
			pointRes = GetDCASmrByPar(lmeeLeg.genPt); //cm
			if (fDCASmrFromMC) pointRes*100;
		}
		
		if (pointRes > 0.) rndm = gRandom->Gaus(DCAxyGen,pointRes);
		else rndm = -999.;
		DCAxyGen = rndm;
		DCAxySigGen = DCAxyGen/pointRes; //cm Correct
        
        lmeeLeg.genDCAxy = DCAxyGen;
        lmeeLeg.genDCAxySig = DCAxySigGen;

        /*if (pythiaCCfile && !IsCharmedEle(iMC)) continue;
        else if (pythiaBBfile && !IsBeautyEle(iMC)) continue;
        else if (pythiaBfile) continue; //TO DO pythiaBfile*/

        lmeeLeg.MakeGenLV();

        //rec track loop
        for (Int_t iTracks = 0; iTracks < aodEvent->GetNumberOfTracks(); iTracks++){
            AliAODTrack* track = (AliAODTrack*)aodEvent->GetTrack(iTracks);
            if (!track) { Printf("ERROR: Could not receive track %d", iTracks); continue; }
            Int_t label = track->GetLabel();
            Int_t abslabel = TMath::Abs( track->GetLabel() );
            if(abslabel != iMC) continue;
            if(part->Charge()/3 != track->Charge()) continue;

            // track cuts: use AliESDtrackCuts functionality
            if (!(fTrackCuts->IsSelected(track))) continue;

            Double_t trackPt  = track->Pt();
            Double_t trackP   = track->P();
            Double_t trackEta = track->Eta();
            Double_t trackPhi = track->Phi();
            
            lmeeLeg.recP   = trackP;
            lmeeLeg.recPt  = trackPt;
            lmeeLeg.recEta = trackEta;
            lmeeLeg.recPhi = trackPhi;
            lmeeLeg.MakeRecLV();

            if (fFillSmrMaps){
				fPRec_Gen->Fill(lmeeLeg.genP,lmeeLeg.recP-lmeeLeg.genP);
				fPtRec_Gen->Fill(lmeeLeg.genPt,lmeeLeg.recPt-lmeeLeg.genPt);
				fPtRecOverGen->Fill(lmeeLeg.genPt,lmeeLeg.recPt/lmeeLeg.genPt);
				fEtaRec_Gen->Fill(lmeeLeg.genPt,lmeeLeg.recEta-lmeeLeg.genEta);
				if (track->Charge() > 0){
					fPosPhiRec_Gen->Fill(lmeeLeg.genPt,lmeeLeg.recPhi-lmeeLeg.genPhi);
					//fPosDCARes->Fill(lmeeLeg.genPt,lmeeLeg.recDCAxySig-lmeeLeg.genDCAxySig);
				}
				else{
					fElePhiRec_Gen->Fill(lmeeLeg.genPt,lmeeLeg.recPhi-lmeeLeg.genPhi);
					//fEleDCARes->Fill(lmeeLeg.genPt,lmeeLeg.recDCAxySig-lmeeLeg.genDCAxySig);
				}
			}

            Float_t nSigmaTPCele = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            Float_t nSigmaTPCpio = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
            Float_t nSigmaTPCpro = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
            Float_t nSigmaTPCkao = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);

            Bool_t isTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
            Float_t nSigmaTOFele = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
            Bool_t TOFpid = (isTOF && TMath::Abs(nSigmaTOFele) < 3);

			if (fRecabPID && trackP < 1.5){
				Double_t TPCm, TPCw, TOFm, TOFw;
				GetRecalibrationPID(trackP,trackEta,&TPCm,&TPCw,&TOFm,&TOFw);
				nSigmaTPCele = (nSigmaTPCele-TPCm)/TPCw;
				nSigmaTOFele = (nSigmaTOFele-TOFm)/TOFw;
			}
			
			if (nSigmaTPCele < fTPCmin || nSigmaTPCele > fTPCmax) continue;				//  -3 < TPCele < 3 pass
			if (nSigmaTPCpio > fTPCminPio && nSigmaTPCpio < fTPCmaxPio) continue;		//-inf < TPCpio < 4 reject

			Bool_t passedA = kFALSE;
			Bool_t passedB = kFALSE;

			if (trackP < 0.4) passedA = kTRUE;
			
			if (isTOF && trackP >= 0.4){
				if (nSigmaTOFele >= fTOFmin && nSigmaTOFele <= fTOFmax) passedA = kTRUE;
			}
			else if ((nSigmaTPCpro < fTPCminPro || nSigmaTPCpro > fTPCmaxPro) && (nSigmaTPCkao < fTPCminKao || nSigmaTPCkao > fTPCmaxKao)) passedB = kTRUE;

			if (!passedA && !passedB) continue;
			lmeeLeg.passedCuts = kTRUE;

			if (passedA && fFillPIDrecMaps){
				TPCnSigmaEle_Eta_P_lin->Fill(trackEta,nSigmaTPCele,trackP);
				if (trackPt >= 0.4) TOFnSigmaEle_Eta_P_lin->Fill(trackEta,nSigmaTOFele,trackP);
			}
			if (passedB && fFillPIDrecMaps){
				TPCnSigmaEle_Eta_P_lin2->Fill(trackP, nSigmaTPCele);
				if (isTOF && trackPt >= 0.4) TOFnSigmaEle_Eta_P_lin2->Fill(trackP, nSigmaTOFele);
			}

            // fill tracking info
            fHistNclsTPC    ->Fill(track->GetTPCNcls());
            fHistNclsSTPC   ->Fill(track->GetTPCnclsS());
            fHistNclsITS    ->Fill(track->GetITSNcls());
            fHistChi2perNDF ->Fill(track->Chi2perNDF());

			/*Float_t DCAxy, DCAz;
			track->GetImpactParameters(&DCAxy,&DCAz);
    		fHistDCAxy  ->Fill(DCAxy);
    		fHistDCAz   ->Fill(DCAz);
    		lmeeLeg.recDCAxy = DCAxy;*/
    				
    		Double_t DCA[2]       = {-999.,-999.};      // xy,z absolute values
			Double_t DCASig[2]    = {-999.,-999.};      // xy,z sigma values
			Double_t DCARes[3]    = {-999.,-999.,-999.};// Covariance matrix
					
			AliVParticle* vPart = (AliVParticle*)aodEvent->GetTrack(iTracks);
			GetDCA(static_cast<AliAODTrack*>(vPart), DCA, DCARes);
			if(DCARes[0] > 0.) DCASig[0] = DCA[0]/TMath::Sqrt(DCARes[0]);
			if(DCARes[2] > 0.) DCASig[1] = DCA[1]/TMath::Sqrt(DCARes[2]);
					
			lmeeLeg.recDCAxy = DCA[0];
			lmeeLeg.recDCAxySig = DCASig[0];

    		fHistTPCnSigmaEle->Fill(trackP, nSigmaTPCele);
    		fHistTPCnSigmaPio->Fill(trackP, nSigmaTPCpio);
    		if (isTOF) fHistTOFnSigmaEle->Fill(trackP, nSigmaTOFele);

            if (track->Charge() > 0){
				fHistPosRecPt    ->Fill(trackPt);
				fHistPosRecEtaPhi->Fill(trackEta,trackPhi);
				fHistPosRecPhi   ->Fill(trackPhi);
				fHistPosRecDCAxy ->Fill(DCA[0]);
				if(DCARes[0] > 0.){
					fPosDCARes->Fill(lmeeLeg.recPt,DCARes[0]);
					fPosDCARes2->Fill(lmeeLeg.recPt,DCARes[0]);
					fPosDCARes3->Fill(lmeeLeg.recPt,DCARes[0]);
				}
				if (pythiaCCfile && lmeeLeg.isCharm){
					fHistPosRecPtCC    ->Fill(trackPt);
					fHistPosRecEtaPhiCC->Fill(trackEta,trackPhi);
					fHistPosRecDCAxyCC ->Fill(DCA[0]);
      			}
      			else if (pythiaBBfile && lmeeLeg.isBeauty){
					fHistPosRecPtBB    ->Fill(trackPt);
					fHistPosRecEtaPhiBB->Fill(trackEta,trackPhi);
					fHistPosRecDCAxyBB ->Fill(DCA[0]);
      			}
            }
            else {
				fHistEleRecPt    ->Fill(trackPt);
				fHistEleRecEtaPhi->Fill(trackEta,trackPhi);
				fHistEleRecPhi   ->Fill(trackPhi);
				fHistEleRecDCAxy ->Fill(DCA[0]);
				if(DCARes[0] > 0.){
					fEleDCARes->Fill(lmeeLeg.recPt,DCARes[0]);
					fEleDCARes2->Fill(lmeeLeg.recPt,DCARes[0]);
					fEleDCARes3->Fill(lmeeLeg.recPt,DCARes[0]);
				}
				if (pythiaCCfile && lmeeLeg.isCharm){
					fHistEleRecPtCC    ->Fill(trackPt);
					fHistEleRecEtaPhiCC->Fill(trackEta,trackPhi);
					fHistEleRecDCAxyCC ->Fill(DCA[0]);
				}
				else if (pythiaBBfile && lmeeLeg.isBeauty){
					fHistEleRecPtBB    ->Fill(trackPt);
					fHistEleRecEtaPhiBB->Fill(trackEta,trackPhi);
					fHistEleRecDCAxyBB ->Fill(DCA[0]);
				}
            }
        } // rec loop

        if(part->PdgCode() == +11){
			LMEEpositrons.push_back(lmeeLeg);
			LMEEpositrons2.push_back(lmeeLeg);
      		fHistPosGenPt    ->Fill(lmeeLeg.genPt);
      		fHistPosGenEtaPhi->Fill(lmeeLeg.genEta,lmeeLeg.genPhi);
      		fHistEleGenPhi   ->Fill(lmeeLeg.genPhi);
      		if (pythiaCCfile && lmeeLeg.isCharm){
				fHistPosGenPtCC    ->Fill(lmeeLeg.genPt);
      			fHistPosGenEtaPhiCC->Fill(lmeeLeg.genEta,lmeeLeg.genPhi);
    		}
    		else if (pythiaBBfile && lmeeLeg.isBeauty){
				fHistPosGenPtBB    ->Fill(lmeeLeg.genPt);
				fHistPosGenEtaPhiBB->Fill(lmeeLeg.genEta,lmeeLeg.genPhi);
			}
        }
        else if(part->PdgCode() == -11){
            LMEEelectrons.push_back(lmeeLeg);
            LMEEelectrons2.push_back(lmeeLeg);
            fHistEleGenPt    ->Fill(lmeeLeg.genPt);
            fHistEleGenEtaPhi->Fill(lmeeLeg.genEta,lmeeLeg.genPhi);
            fHistPosGenPhi   ->Fill(lmeeLeg.genPhi);
            if (pythiaCCfile){
				fHistEleGenPtCC    ->Fill(lmeeLeg.genPt);
				fHistEleGenEtaPhiCC->Fill(lmeeLeg.genEta,lmeeLeg.genPhi);
			}
			else if (pythiaBBfile){
				fHistEleGenPtBB    ->Fill(lmeeLeg.genPt);
				fHistEleGenEtaPhiBB->Fill(lmeeLeg.genEta,lmeeLeg.genPhi);
			}
        }

        // mother PDG
        if (pythiaCCfile){
			switch (TMath::Abs(lmeeLeg.mPDG)) {
				case 411:
					fMotherPDGCC->Fill(0);
					break;
    			case 421:
    				fMotherPDGCC->Fill(1);
    				break;
    			case 431:
    				fMotherPDGCC->Fill(2);
    				break;
    			case 4122:
    				fMotherPDGCC->Fill(3);
    				break;
    			default:
					fMotherPDGCC->Fill(4);
					break;
			}
		}
        else if (pythiaBBfile){
			if (TMath::Abs(lmeeLeg.mPDG) > 500 && TMath::Abs(lmeeLeg.mPDG) < 600) fMotherPDGBB->Fill(0);
			else if (TMath::Abs(lmeeLeg.mPDG) > 5000 && TMath::Abs(lmeeLeg.mPDG) < 6000) fMotherPDGBB->Fill(1);
			else if ((TMath::Abs(lmeeLeg.mPDG) > 410 && TMath::Abs(lmeeLeg.mPDG) < 432) || (TMath::Abs(lmeeLeg.mPDG) == 4122)) fMotherPDGBB->Fill(2);
    		else fMotherPDGBB->Fill(3);
    	}
    } // MC track loop

    // pairing ULS
    // in the case of 1 charmed hadron-antihadron pair there are no like-sign pairs
    if(LMEEelectrons.size() > 0 && LMEEpositrons.size() > 0 && pythiaCCfile){
		for(std::vector<LMEEparticle>::iterator it1 = LMEEelectrons.begin(); it1 != LMEEelectrons.end(); ++it1){
            Bool_t charm1 = it1->isCharm;

            for(std::vector<LMEEparticle>::iterator it2 = LMEEpositrons.begin(); it2 != LMEEpositrons.end(); ++it2){
                Bool_t charm2 = it2->isCharm;

                Double_t mee_gen  = (it1->genLv + it2->genLv).M();
                Double_t ptee_gen = (it1->genLv + it2->genLv).Pt();
                Double_t mee_rec  = (it1->recLv + it2->recLv).M();
                Double_t ptee_rec = (it1->recLv + it2->recLv).Pt();
                Double_t dPhiee_gen = it1->genPhi - it2->genPhi;
                if (dPhiee_gen > TMath::Pi()) dPhiee_gen = dPhiee_gen - 2*TMath::Pi();
                if (dPhiee_gen < -TMath::Pi()) dPhiee_gen = dPhiee_gen + 2*TMath::Pi();
                Double_t dPhiee_rec = it1->recPhi - it2->recPhi;
                if (dPhiee_rec > TMath::Pi()) dPhiee_rec = dPhiee_rec - 2*TMath::Pi();
                if (dPhiee_rec < -TMath::Pi()) dPhiee_rec = dPhiee_rec + 2*TMath::Pi();
                Double_t DCAxyee_gen = TMath::Sqrt(((it1->genDCAxy)*(it1->genDCAxy) + (it2->genDCAxy)*(it2->genDCAxy))/2);
                Double_t DCAxySigee_gen = TMath::Sqrt(((it1->genDCAxySig)*(it1->genDCAxySig) + (it2->genDCAxySig)*(it2->genDCAxySig))/2);
                Double_t DCAxyee = TMath::Sqrt(((it1->recDCAxy)*(it1->recDCAxy) + (it2->recDCAxy)*(it2->recDCAxy))/2);
                Double_t DCAxySigee = TMath::Sqrt(((it1->recDCAxySig)*(it1->recDCAxySig) + (it2->recDCAxySig)*(it2->recDCAxySig))/2);

                fHistPairsGen->Fill(mee_gen,ptee_gen);
                fHistPairsGenPt->Fill(ptee_gen);
                fHistPairsGenMass->Fill(mee_gen);
                fHistPairsGenMassPt->Fill(mee_gen,ptee_gen);
                fHistPairsGenMassPt2->Fill(mee_gen,ptee_gen);
                fHistPairsGenDphi->Fill(dPhiee_gen);
                fHistPairsGenPhiPt->Fill(dPhiee_gen,ptee_gen);
                fHistPairsGenPhiPt2->Fill(dPhiee_gen,ptee_gen);
                fHistPairsGenPtMasPhi->Fill(ptee_gen,mee_gen,dPhiee_gen);
                
                fHistPairsGenDCAxy->Fill(DCAxyee_gen);
				fHistPairsGenDCAxySig->Fill(DCAxySigee_gen);
				fHistPairsGenPtDCAxySig->Fill(ptee_rec,DCAxySigee_gen);

                if (it1->passedCuts && it2->passedCuts){
                    fHistPairsRec->Fill(mee_rec,ptee_rec);
					fHistPairsRecPt->Fill(ptee_rec);
					fHistPairsRecMass->Fill(mee_rec);
					fHistPairsRecMassPt->Fill(mee_rec,ptee_rec);
					fHistPairsRecMassPt2->Fill(mee_rec,ptee_rec);
					fHistPairsRecDphi->Fill(dPhiee_rec);
					fHistPairsRecPhiPt->Fill(dPhiee_rec,ptee_rec);
					fHistPairsRecPhiPt2->Fill(dPhiee_rec,ptee_rec);
					fHistPairsRecPtMasPhi->Fill(ptee_rec,mee_rec,dPhiee_rec);

					fHistPairsDCAxy->Fill(DCAxyee);
					fHistPairsDCAxySig->Fill(DCAxySigee);
					fHistPairsPtDCAxySig->Fill(ptee_rec,DCAxySigee);
					fPairDCARec_Gen->Fill(ptee_rec,DCAxySigee-DCAxySigee_gen);
                }

                if(charm1 && charm2){
                    if(it1->mlabel == it2->mlabel){
                        // resonance decay! should not be the case
                    }
                    else{
                        fHistPairsCCgen->Fill(mee_gen,ptee_gen);
                        fHistPairsCCgen2->Fill(mee_gen,ptee_gen);
                        fHistPairsGenPtCC->Fill(ptee_gen);
						fHistPairsGenMassCC->Fill(mee_gen);
						fHistPairsGenMassPtCC->Fill(mee_gen,ptee_gen);
						fHistPairsGenMassPt2CC->Fill(mee_gen,ptee_gen);
						fHistPairsGenDphiCC->Fill(dPhiee_gen);
						fHistPairsGenPhiPtCC->Fill(dPhiee_gen,ptee_gen);
						fHistPairsGenPhiPt2CC->Fill(dPhiee_gen,ptee_gen);
						fHistPairsGenPtMasPhiCC->Fill(ptee_gen,mee_gen,dPhiee_gen);
						
						fHistPairsGenDCAxyCC->Fill(DCAxyee_gen);
						fHistPairsGenDCAxySigCC->Fill(DCAxySigee_gen);
						fHistPairsGenPtDCAxySigCC->Fill(ptee_rec,DCAxySigee_gen);

                        if (it1->passedCuts && it2->passedCuts){
                            fHistPairsCCrec->Fill(mee_rec,ptee_rec);
                            fHistPairsCCrec2->Fill(mee_rec,ptee_rec);
                            fHistPairsRecPtCC->Fill(ptee_rec);
							fHistPairsRecMassCC->Fill(mee_rec);
							fHistPairsRecMassPtCC->Fill(mee_rec,ptee_rec);
							fHistPairsRecMassPt2CC->Fill(mee_rec,ptee_rec);
							fHistPairsRecDphiCC->Fill(dPhiee_rec);
							fHistPairsRecPhiPtCC->Fill(dPhiee_rec,ptee_rec);
							fHistPairsRecPhiPt2CC->Fill(dPhiee_rec,ptee_rec);
							fHistPairsRecPtMasPhiCC->Fill(ptee_rec,mee_rec,dPhiee_rec);
							fHistPairsDCAxyCC->Fill(DCAxyee);
							fHistPairsDCAxySigCC->Fill(DCAxySigee);
							fHistPairsPtDCAxySigCC->Fill(ptee_rec,DCAxySigee);
                        }
                    }
                }
            } // positron iteration
        } // electron iteration
    } // pairing ULS
    else if (pythiaBBfile){
        for(std::vector<LMEEparticle>::iterator it1 = LMEEelectrons.begin(); it1 != LMEEelectrons.end(); ++it1){

            for(std::vector<LMEEparticle>::iterator it2 = LMEEpositrons.begin(); it2 != LMEEpositrons.end(); ++it2){

                Double_t mee_gen  = (it1->genLv + it2->genLv).M();
                Double_t ptee_gen = (it1->genLv + it2->genLv).Pt();
                Double_t mee_rec  = (it1->recLv + it2->recLv).M();
                Double_t ptee_rec = (it1->recLv + it2->recLv).Pt();
                Double_t dPhiee_gen = it1->genPhi - it2->genPhi;
                if (dPhiee_gen > TMath::Pi()) dPhiee_gen = dPhiee_gen - 2*TMath::Pi();
                if (dPhiee_gen < -TMath::Pi()) dPhiee_gen = dPhiee_gen + 2*TMath::Pi();
                Double_t dPhiee_rec = it1->recPhi - it2->recPhi;
                if (dPhiee_rec > TMath::Pi()) dPhiee_rec = dPhiee_rec - 2*TMath::Pi();
                if (dPhiee_rec < -TMath::Pi()) dPhiee_rec = dPhiee_rec + 2*TMath::Pi();
                Double_t DCAxyee_gen = TMath::Sqrt(((it1->genDCAxy)*(it1->genDCAxy) + (it2->genDCAxy)*(it2->genDCAxy))/2);
                Double_t DCAxySigee_gen = TMath::Sqrt(((it1->genDCAxySig)*(it1->genDCAxySig) + (it2->genDCAxySig)*(it2->genDCAxySig))/2);
				Double_t DCAxyee = TMath::Sqrt(((it1->recDCAxy)*(it1->recDCAxy) + (it2->recDCAxy)*(it2->recDCAxy))/2);
                Double_t DCAxySigee = TMath::Sqrt(((it1->recDCAxySig)*(it1->recDCAxySig) + (it2->recDCAxySig)*(it2->recDCAxySig))/2);

				fHistPairsGen->Fill(mee_gen,ptee_gen);
                fHistPairsGenPt->Fill(ptee_gen);
                fHistPairsGenMass->Fill(mee_gen);
                fHistPairsGenMassPt->Fill(mee_gen,ptee_gen);
                fHistPairsGenMassPt2->Fill(mee_gen,ptee_gen);
                fHistPairsGenDphi->Fill(dPhiee_gen);
                fHistPairsGenPhiPt->Fill(dPhiee_gen,ptee_gen);
                fHistPairsGenPhiPt2->Fill(dPhiee_gen,ptee_gen);
                fHistPairsGenPtMasPhi->Fill(ptee_gen,mee_gen,dPhiee_gen);
                
                fHistPairsGenDCAxy->Fill(DCAxyee_gen);
				fHistPairsGenDCAxySig->Fill(DCAxySigee_gen);
				fHistPairsGenPtDCAxySig->Fill(ptee_rec,DCAxySigee_gen);
				fPairDCARec_Gen->Fill(ptee_rec,DCAxySigee-DCAxySigee_gen);

                if (it1->passedCuts && it2->passedCuts){
                    fHistPairsRec->Fill(mee_rec,ptee_rec);
					fHistPairsRecPt->Fill(ptee_rec);
					fHistPairsRecMass->Fill(mee_rec);
					fHistPairsRecMassPt->Fill(mee_rec,ptee_rec);
					fHistPairsRecMassPt2->Fill(mee_rec,ptee_rec);
					fHistPairsRecDphi->Fill(dPhiee_rec);
					fHistPairsRecPhiPt->Fill(dPhiee_rec,ptee_rec);
					fHistPairsRecPhiPt2->Fill(dPhiee_rec,ptee_rec);
					fHistPairsRecPtMasPhi->Fill(ptee_rec,mee_rec,dPhiee_rec);

					fHistPairsDCAxy->Fill(DCAxyee);
					fHistPairsDCAxySig->Fill(DCAxySigee);
					fHistPairsPtDCAxySig->Fill(ptee_rec,DCAxySigee);
                }

                if (it1->isBeauty && it2->isBeauty){
					if(it1->mlabel == it2->mlabel){
						// resonance decay! should not be the case
					}
					else{
						fHistPairsBBgen->Fill(mee_gen,ptee_gen);
          				fHistPairsBBgen2->Fill(mee_gen,ptee_gen);
          				fHistPairsGenPtBB->Fill(ptee_gen);
						fHistPairsGenMassBB->Fill(mee_gen);
						fHistPairsGenMassPtBB->Fill(mee_gen,ptee_gen);
						fHistPairsGenMassPt2BB->Fill(mee_gen,ptee_gen);
						fHistPairsGenDphiBB->Fill(dPhiee_gen);
						fHistPairsGenPhiPtBB->Fill(dPhiee_gen,ptee_gen);
						fHistPairsGenPhiPt2BB->Fill(dPhiee_gen,ptee_gen);
						fHistPairsGenPtMasPhiBB->Fill(ptee_gen,mee_gen,dPhiee_gen);
							
						fHistPairsGenDCAxyBB->Fill(DCAxyee_gen);
						fHistPairsGenDCAxySigBB->Fill(DCAxySigee_gen);
						fHistPairsGenPtDCAxySigBB->Fill(ptee_rec,DCAxySigee_gen);
							
          				if (it1->passedCuts && it2->passedCuts){
          					fHistPairsBBrec->Fill(mee_rec,ptee_rec);
          					fHistPairsBBrec2->Fill(mee_rec,ptee_rec);
							fHistPairsRecPtBB->Fill(ptee_rec);
							fHistPairsRecMassBB->Fill(mee_rec);
							fHistPairsRecMassPtBB->Fill(mee_rec,ptee_rec);
							fHistPairsRecMassPt2BB->Fill(mee_rec,ptee_rec);
							fHistPairsRecDphiBB->Fill(dPhiee_rec);
							fHistPairsRecPhiPtBB->Fill(dPhiee_rec,ptee_rec);
							fHistPairsRecPhiPt2BB->Fill(dPhiee_rec,ptee_rec);
							fHistPairsRecPtMasPhiBB->Fill(ptee_rec,mee_rec,dPhiee_rec);
							fHistPairsDCAxyBB->Fill(DCAxyee);
							fHistPairsDCAxySigBB->Fill(DCAxySigee);
							fHistPairsPtDCAxySigBB->Fill(ptee_rec,DCAxySigee);
          				}
          			}
          		}
            } // positron iteration
        } // electron iteration
        for(std::vector<LMEEparticle>::iterator it1 = LMEEelectrons.begin(); it1 != LMEEelectrons.end(); ++it1){
			if (!it1->isBeauty) continue;
            for(std::vector<LMEEparticle>::iterator it2 = LMEEelectrons2.begin(); it2 != LMEEelectrons2.end(); ++it2){
                if (it1->label <= it2->label) continue;
                if (!it2->isBeauty) continue;

                Double_t mee_gen  = (it1->genLv + it2->genLv).M();
                Double_t ptee_gen = (it1->genLv + it2->genLv).Pt();
                Double_t mee_rec  = (it1->recLv + it2->recLv).M();
                Double_t ptee_rec = (it1->recLv + it2->recLv).Pt();
                Double_t dPhiee_gen = it1->genPhi - it2->genPhi;
                if (dPhiee_gen > TMath::Pi()) dPhiee_gen = dPhiee_gen - 2*TMath::Pi();
                if (dPhiee_gen < -TMath::Pi()) dPhiee_gen = dPhiee_gen + 2*TMath::Pi();
                Double_t dPhiee_rec = it1->recPhi - it2->recPhi;
                if (dPhiee_rec > TMath::Pi()) dPhiee_rec = dPhiee_rec - 2*TMath::Pi();
                if (dPhiee_rec < -TMath::Pi()) dPhiee_rec = dPhiee_rec + 2*TMath::Pi();
                Double_t DCAxyee_gen = TMath::Sqrt(((it1->genDCAxy)*(it1->genDCAxy) + (it2->genDCAxy)*(it2->genDCAxy))/2);
                Double_t DCAxySigee_gen = TMath::Sqrt(((it1->genDCAxySig)*(it1->genDCAxySig) + (it2->genDCAxySig)*(it2->genDCAxySig))/2);
				Double_t DCAxyee = TMath::Sqrt(((it1->recDCAxy)*(it1->recDCAxy) + (it2->recDCAxy)*(it2->recDCAxy))/2);
                Double_t DCAxySigee = TMath::Sqrt(((it1->recDCAxySig)*(it1->recDCAxySig) + (it2->recDCAxySig)*(it2->recDCAxySig))/2);

                if(it1->mlabel == it2->mlabel){
					// resonance decay! should not be the case
                }
                else{
					fHistPairsBBgenLS->Fill(mee_gen,ptee_gen);
					fHistPairsBBgenLS2->Fill(mee_gen,ptee_gen);
					fHistPairsGenPtLSBB->Fill(ptee_gen);
					fHistPairsGenMassLSBB->Fill(mee_gen);
					fHistPairsGenMassPtLSBB->Fill(mee_gen,ptee_gen);
					fHistPairsGenMassPt2LSBB->Fill(mee_gen,ptee_gen);
					fHistPairsGenDphiLSBB->Fill(dPhiee_gen);
					fHistPairsGenPhiPtLSBB->Fill(dPhiee_gen,ptee_gen);
					fHistPairsGenPhiPt2LSBB->Fill(dPhiee_gen,ptee_gen);
					fHistPairsGenPtMasPhiLSBB->Fill(ptee_gen,mee_gen,dPhiee_gen);
					
					fHistPairsGenDCAxyLSBB->Fill(DCAxyee_gen);
					fHistPairsGenDCAxySigLSBB->Fill(DCAxySigee_gen);
					fHistPairsGenPtDCAxySigLSBB->Fill(ptee_rec,DCAxySigee_gen);
					
          			if (it1->passedCuts && it2->passedCuts){
          				fHistPairsBBrecLS->Fill(mee_rec,ptee_rec);
          				fHistPairsBBrecLS2->Fill(mee_rec,ptee_rec);
						fHistPairsRecPtLSBB->Fill(ptee_rec);
						fHistPairsRecMassLSBB->Fill(mee_rec);
						fHistPairsRecMassPtLSBB->Fill(mee_rec,ptee_rec);
						fHistPairsRecMassPt2LSBB->Fill(mee_rec,ptee_rec);
						fHistPairsRecDphiLSBB->Fill(dPhiee_rec);
						fHistPairsRecPhiPtLSBB->Fill(dPhiee_rec,ptee_rec);
						fHistPairsRecPhiPt2LSBB->Fill(dPhiee_rec,ptee_rec);
						fHistPairsRecPtMasPhiLSBB->Fill(ptee_rec,mee_rec,dPhiee_rec);

						fHistPairsDCAxyLSBB->Fill(DCAxyee);
						fHistPairsDCAxySigLSBB->Fill(DCAxySigee);
						fHistPairsPtDCAxySigLSBB->Fill(ptee_rec,DCAxySigee);
          			}
                }
            }
        }
        for(std::vector<LMEEparticle>::iterator it1 = LMEEpositrons.begin(); it1 != LMEEpositrons.end(); ++it1){
			if (!it1->isBeauty) continue;
			
            for(std::vector<LMEEparticle>::iterator it2 = LMEEpositrons2.begin(); it2 != LMEEpositrons2.end(); ++it2){
                if (it1->label <= it2->label) continue;
                if (!it2->isBeauty) continue;

                Double_t mee_gen  = (it1->genLv + it2->genLv).M();
                Double_t ptee_gen = (it1->genLv + it2->genLv).Pt();
                Double_t mee_rec  = (it1->recLv + it2->recLv).M();
                Double_t ptee_rec = (it1->recLv + it2->recLv).Pt();
                Double_t dPhiee_gen = it1->genPhi - it2->genPhi;
                if (dPhiee_gen > TMath::Pi()) dPhiee_gen = dPhiee_gen - 2*TMath::Pi();
                if (dPhiee_gen < -TMath::Pi()) dPhiee_gen = dPhiee_gen + 2*TMath::Pi();
                Double_t dPhiee_rec = it1->recPhi - it2->recPhi;
                if (dPhiee_rec > TMath::Pi()) dPhiee_rec = dPhiee_rec - 2*TMath::Pi();
                if (dPhiee_rec < -TMath::Pi()) dPhiee_rec = dPhiee_rec + 2*TMath::Pi();
                Double_t DCAxyee_gen = TMath::Sqrt(((it1->genDCAxy)*(it1->genDCAxy) + (it2->genDCAxy)*(it2->genDCAxy))/2);
                Double_t DCAxySigee_gen = TMath::Sqrt(((it1->genDCAxySig)*(it1->genDCAxySig) + (it2->genDCAxySig)*(it2->genDCAxySig))/2);
				Double_t DCAxyee = TMath::Sqrt(((it1->recDCAxy)*(it1->recDCAxy) + (it2->recDCAxy)*(it2->recDCAxy))/2);
                Double_t DCAxySigee = TMath::Sqrt(((it1->recDCAxySig)*(it1->recDCAxySig) + (it2->recDCAxySig)*(it2->recDCAxySig))/2);

                if(it1->mlabel == it2->mlabel){
					// resonance decay! should not be the case
                }
                else{
					fHistPairsBBgenLS->Fill(mee_gen,ptee_gen);
					fHistPairsBBgenLS2->Fill(mee_gen,ptee_gen);
					fHistPairsGenPtLSBB->Fill(ptee_gen);
					fHistPairsGenMassLSBB->Fill(mee_gen);
					fHistPairsGenMassPtLSBB->Fill(mee_gen,ptee_gen);
					fHistPairsGenMassPt2LSBB->Fill(mee_gen,ptee_gen);
					fHistPairsGenDphiLSBB->Fill(dPhiee_gen);
					fHistPairsGenPhiPtLSBB->Fill(dPhiee_gen,ptee_gen);
					fHistPairsGenPhiPt2LSBB->Fill(dPhiee_gen,ptee_gen);
					fHistPairsGenPtMasPhiLSBB->Fill(ptee_gen,mee_gen,dPhiee_gen);
					
					fHistPairsGenDCAxyLSBB->Fill(DCAxyee_gen);
					fHistPairsGenDCAxySigLSBB->Fill(DCAxySigee_gen);
					fHistPairsGenPtDCAxySigLSBB->Fill(ptee_rec,DCAxySigee_gen);
					
          			if (it1->passedCuts && it2->passedCuts){
          				fHistPairsBBrecLS->Fill(mee_rec,ptee_rec);
          				fHistPairsBBrecLS2->Fill(mee_rec,ptee_rec);
						fHistPairsRecPtLSBB->Fill(ptee_rec);
						fHistPairsRecMassLSBB->Fill(mee_rec);
						fHistPairsRecMassPtLSBB->Fill(mee_rec,ptee_rec);
						fHistPairsRecMassPt2LSBB->Fill(mee_rec,ptee_rec);
						fHistPairsRecDphiLSBB->Fill(dPhiee_rec);
						fHistPairsRecPhiPtLSBB->Fill(dPhiee_rec,ptee_rec);
						fHistPairsRecPhiPt2LSBB->Fill(dPhiee_rec,ptee_rec);
						fHistPairsRecPtMasPhiLSBB->Fill(ptee_rec,mee_rec,dPhiee_rec);

						fHistPairsDCAxyLSBB->Fill(DCAxyee);
						fHistPairsDCAxySigLSBB->Fill(DCAxySigee);
						fHistPairsPtDCAxySigLSBB->Fill(ptee_rec,DCAxySigee);
          			}
                }
            } // positron iteration
        } // electron iteration
    }

    PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskeeCor::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}

Bool_t AliAnalysisTaskeeCor::IsCharmedEle(Int_t label){
	AliMCParticle *part = (AliMCParticle*)(mcEvent->GetTrack(label));
	if (!part)  return kFALSE;
	Int_t pdg = TMath::Abs(part->PdgCode());
	if (pdg != 11) return kFALSE;

	AliMCParticle *mother = (AliMCParticle*)(mcEvent->GetTrack(part->GetMother()));
	if (!mother) return kFALSE;
	Int_t mpdg = TMath::Abs(mother->PdgCode());
	if (!((mpdg > 410 && mpdg < 436) || (mpdg > 4113 && mpdg < 4445))) return kFALSE;

	AliMCParticle *gmother = (AliMCParticle*)(mcEvent->GetTrack(mother->GetMother()));
	if (!gmother) return kFALSE;
	Int_t gmpdg = TMath::Abs(gmother->PdgCode());
	if (gmpdg != 4) return kFALSE;

	return kTRUE;
}

Bool_t AliAnalysisTaskeeCor::IsBeautyEle(Int_t label){
	AliMCParticle *part = (AliMCParticle*)(mcEvent->GetTrack(label));
	if (!part)  return kFALSE;
	Int_t pdg = TMath::Abs(part->PdgCode());
	if (pdg != 11) return kFALSE;

	Bool_t motherD = 0;
	Bool_t motherB = 0;
	Bool_t gmotherB = 0;
	Bool_t ggmotherB = 0;
	Bool_t gggmotherB = 0;

	AliMCParticle *mother = (AliMCParticle*)(mcEvent->GetTrack(part->GetMother()));
	if (!mother) return kFALSE;
	Int_t mpdg = TMath::Abs(mother->PdgCode());
	if ((mpdg > 410 && mpdg < 436) || (mpdg > 4113 && mpdg < 4445)) motherD = 1;
	else if ((mpdg > 500 && mpdg < 600) || (mpdg > 5000 && mpdg < 6000)) motherB = 1;
	else return kFALSE;

	AliMCParticle *gmother = (AliMCParticle*)(mcEvent->GetTrack(mother->GetMother()));
	if (!gmother) return kFALSE;
	Int_t gmpdg = TMath::Abs(gmother->PdgCode());
	if ((gmpdg > 500 && gmpdg < 600) || (gmpdg > 5000 && gmpdg < 6000)) gmotherB = 1;
	else if (gmpdg == 5) return kTRUE;
	else return kFALSE;

	AliMCParticle *ggmother = (AliMCParticle*)(mcEvent->GetTrack(gmother->GetMother()));
	if (!ggmother) return kFALSE;
	Int_t ggmpdg = TMath::Abs(ggmother->PdgCode());
	if ((ggmpdg > 500 && ggmpdg < 600) || (ggmpdg > 5000 && ggmpdg < 6000)) ggmotherB = 1;
	else if (ggmpdg == 5) return kTRUE;
	else return kFALSE;

	AliMCParticle *gggmother = (AliMCParticle*)(mcEvent->GetTrack(ggmother->GetMother()));
	if (!gggmother) return kFALSE;
	Int_t gggmpdg = TMath::Abs(gggmother->PdgCode());
	if (gggmpdg != 5) return kFALSE;

	return kTRUE;
}

void AliAnalysisTaskeeCor::GetRecalibrationPID(Double_t mom, Double_t eta, Double_t *meanTPC, Double_t *widthTPC, Double_t *meanTOF, Double_t *widthTOF){
	*meanTPC = 0.;
	*widthTPC = 1.;
	*meanTOF = 0.;
	*widthTOF = 1.;

	if (fTPCmean){
		Double_t binX = fTPCmean->GetXaxis()->FindBin(mom);
		Double_t binY = fTPCmean->GetYaxis()->FindBin(eta);
		if (binX > 0 && binY > 0 && binX <= fTPCmean->GetNbinsX() && binY <= fTPCmean->GetNbinsY()) *meanTPC = fTPCmean->GetBinContent(binX,binY);
	}
	if (fTPCwidth){
		Double_t binX = fTPCwidth->GetXaxis()->FindBin(mom);
		Double_t binY = fTPCwidth->GetYaxis()->FindBin(eta);
		if (binX > 0 && binY > 0 && binX <= fTPCwidth->GetNbinsX() && binY <= fTPCwidth->GetNbinsY()) *widthTPC = fTPCwidth->GetBinContent(binX,binY);
	}
	if (fTOFmean){
		Double_t binX = fTOFmean->GetXaxis()->FindBin(mom);
		Double_t binY = fTOFmean->GetYaxis()->FindBin(eta);
		if (binX > 0 && binY > 0 && binX <= fTOFmean->GetNbinsX() && binY <= fTOFmean->GetNbinsY()) *meanTOF = fTOFmean->GetBinContent(binX,binY);
	}
	if (fTOFwidth){
		Double_t binX = fTOFwidth->GetXaxis()->FindBin(mom);
		Double_t binY = fTOFwidth->GetYaxis()->FindBin(eta);
		if (binX > 0 && binY > 0 && binX <= fTOFwidth->GetNbinsX() && binY <= fTOFwidth->GetNbinsY()) *widthTOF = fTOFwidth->GetBinContent(binX,binY);
	}
}

Double_t AliAnalysisTaskeeCor::GetPtSmr(Double_t pt){
	Int_t binPt = fPtSmr->GetXaxis()->FindBin(pt);
	if (binPt <= 0) return 1.;					//underflow protection
	if (binPt > fPtSmr->GetNbinsX()) return 1.; //overflow protection
	fPtSmr->GetXaxis()->SetRange(binPt,binPt);
	TH1F *projPt = (TH1F*)fPtSmr->ProjectionY();
	return projPt->GetRandom();
}

Double_t AliAnalysisTaskeeCor::GetEtaSmr(Double_t pt){
	Int_t binPt = fEtaSmr->GetXaxis()->FindBin(pt);
	if (binPt <= 0) return 0.;					//underflow protection
	if (binPt > fEtaSmr->GetNbinsX()) return 0.;//overflow protection
	fEtaSmr->GetXaxis()->SetRange(binPt,binPt);
	TH1F *projEta = (TH1F*)fEtaSmr->ProjectionY();
	return projEta->GetRandom();
}

Double_t AliAnalysisTaskeeCor::GetPhiSmr(Double_t pt, Double_t q){
	if (q > 0.){
		Int_t binPt = fPhiPosSmr->GetXaxis()->FindBin(pt);
		if (binPt <= 0) return 0.;						//underflow protection
		if (binPt > fPhiPosSmr->GetNbinsX()) return 0.;	//overflow protection
		fPhiPosSmr->GetXaxis()->SetRange(binPt,binPt);
		TH1F *projPhi = (TH1F*)fPhiPosSmr->ProjectionY();
		return projPhi->GetRandom();
	}
	else{
		Int_t binPt = fPhiEleSmr->GetXaxis()->FindBin(pt);
		if (binPt <= 0) return 0.;						//underflow protection
		if (binPt > fPhiEleSmr->GetNbinsX()) return 0.;	//overflow protection
		fPhiEleSmr->GetXaxis()->SetRange(binPt,binPt);
		TH1F *projPhi = (TH1F*)fPhiEleSmr->ProjectionY();
		return projPhi->GetRandom();
	}
}

Double_t AliAnalysisTaskeeCor::GetDCASmr(Double_t pt){
	if (pt < 0.8){
		Int_t binPt = fDCASmrMap0->GetXaxis()->FindBin(pt);
		if (binPt <= 0) return -1.;						//underflow protection
		if (binPt > fDCASmrMap0->GetNbinsX()) binPt = fDCASmrCen->GetNbinsX(); //overflow correction
		fDCASmrMap0->GetXaxis()->SetRange(binPt,binPt);
		TH1F *projDCA = (TH1F*)fDCASmrMap0->ProjectionY();
		return projDCA->GetRandom();
	}
	else{
		Int_t binPt = fDCASmrMap1->GetXaxis()->FindBin(pt);
		if (binPt <= 0) return -1.;						//underflow protection
		if (binPt > fDCASmrMap1->GetNbinsX()) binPt = fDCASmrCen->GetNbinsX(); //overflow correction
		fDCASmrMap1->GetXaxis()->SetRange(binPt,binPt);
		TH1F *projDCA = (TH1F*)fDCASmrMap1->ProjectionY();
		return projDCA->GetRandom();
	}
}

Double_t AliAnalysisTaskeeCor::GetDCASmrByPar(Double_t pt){

    Int_t binPt = fDCASmrCen->GetXaxis()->FindBin(pt);
    if (binPt <= 0) return -1.;						//underflow protection
	if (binPt > fDCASmrCen->GetNbinsX()) binPt = fDCASmrCen->GetNbinsX(); //overflow correction
    Double_t center = fDCASmrCen->GetBinContent(binPt);
    Double_t sigma = fDCASmrSig->GetBinContent(binPt);
    Int_t binPt2 = fDCASmrMax->GetXaxis()->FindBin(pt);
    Double_t maximum = fDCASmrMax->GetBinContent(binPt2);
    if (fDCAparSmr){
		Double_t rndm;
		Double_t minResAcc = 0.00038; //not good, take too long and we can't rise it to something more real like a minimum of 4.8um
		if (fDCASmrFromMC) minResAcc/100;
		do {
			rndm = gRandom->Gaus(center,sigma);
		} while (rndm <= minResAcc);
	}
    else{
		//return center;
		if (maximum > 0) return maximum;
		else return -1.;
	}
}

Bool_t AliAnalysisTaskeeCor::GetDCA(const AliAODTrack *track, Double_t* d0z0, Double_t* covd0z0){
	if (track->TestBit(AliAODTrack::kIsDCA)){
		d0z0[0] = track->DCA();
		d0z0[1] = track->ZAtDCA();
		// the covariance matrix is not stored in case of AliAODTrack::kIsDCA
		return kTRUE;
	}
	
	Bool_t ok = kFALSE;
	if (aodEvent){
		AliExternalTrackParam etp; etp.CopyFromVTrack(track);
		
		Float_t xstart = etp.GetX();
		if(xstart > 3.) {
			d0z0[0] = -999.;
			d0z0[1] = -999.;
			return kFALSE;
		}
		
		AliAODVertex *vtx = (AliAODVertex*)(aodEvent->GetPrimaryVertex());
		Double_t fBzkG = aodEvent->GetMagneticField(); // z componenent of field in kG
		ok = etp.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
	}
	if(!ok){
		d0z0[0] = -999.;
		d0z0[1] = -999.;
	}
	return ok;
}
