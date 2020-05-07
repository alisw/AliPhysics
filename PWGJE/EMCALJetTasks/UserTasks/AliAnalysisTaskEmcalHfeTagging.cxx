//
// HFE tagged jet shape analysis task.
//
// Authors: D. Caffarri, L. Cunqueiro (jet); D. Godoy (HFE)
//
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "TVector2.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalPythiaInfo.h"
#include "TRandom3.h"
#include "AliAODInputHandler.h"
#include "AliPIDResponse.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliMultiEventInputHandler.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"

#include "AliAODEvent.h"
#include "AliAnalysisTaskEmcalHfeTagging.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalHfeTagging)

//________________________________________________________________________
AliAnalysisTaskEmcalHfeTagging::AliAnalysisTaskEmcalHfeTagging() :
AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalHfeTagging", kTRUE),
fAOD(0),
fVevent(0),
fpidResponse(0),
fTracksTender(0),
pVtx(0),
spdVtx(0),
fMC(0),
fStack(0),
fMCparticle(0),
fMCarray(0),
fMCheader(0),
fContainer(0),
fMinFractionShared(0),
fJetShapeType(kData),
fJetShapeSub(kNoSub),
fJetSelection(kInclusive),
fPtThreshold(-9999.),
fRMatching(0.2),
fSelectedShapes(0),
fminpTTrig(20.),
fmaxpTTrig(50.),
fangWindowRecoil(0.6),
fSemigoodCorrect(0),
fHolePos(0),
fHoleWidth(0),
fCentSelectOn(kTRUE),
fCentMin(0),
fCentMax(10),
fOneConstSelectOn(kFALSE),
fDerivSubtrOrder(0),
fMCweight(0),
fRunNumber(0),
fAssPtCut(0.1),
fITSncut(3),
fAssTPCnCut(60),
fTPCnCut(100),
fAssITSrefitCut(kTRUE),
fSigmaTOFcut(3.),
fSigmaTPCcutLowPt(-1.),
fSigmaTPCcutHighPt(-1.5),
fSigmTPCcutExcElec(3.5),
fDcaXYcut(1.),
fDcaZcut(2.),
fIMcut(0.1),
fEtaCut(0.7),
fMinEoPcut(0.9),
fMaxEoPcut(1.3),
fM20cut(0.35),
fMinPtTPC(0.5),
fMaxPtTPC(4.),
fMinPtEMCal(4.),
fMaxPtEMCal(25.),
fNeventV0(0),
fNeventT0(0),
fh2ResponseUW(0x0),
fh2ResponseW(0x0),
fPhiJetCorr6(0x0),
fPhiJetCorr7(0x0),
fEtaJetCorr6(0x0),
fEtaJetCorr7(0x0),
fPtJetCorr(0x0),
fPtJet(0x0),
fPtGenJet(0),
fPhiJet(0x0),
fEtaJet(0x0),
fEtaPhiJet(0x0),
fAreaJet(0x0),
fJetProbDensityDetPart(0x0),
fJetProbDensityPartDet(0x0),
fNbOfConstvspT(0x0),
fnTPCnTOFnocut(0),
fnTPCnocutP(0),
fnTOFnocutP(0),
fnTPCcutP(0),
fnTPCcutPt(0),
fnULSmLSpairsPerElectron(0),
fnPartPerJet(0),
fnElecOverPartPerJet(0),
fnInclElecPerJet(0),
fnPhotElecPerJet(0),
fnIncSubPhotElecPerJet(0),
fnTrueElecPerJet(0),
fnTrueHFElecPerJet(0),
fnTruePElecPerJet(0),
fPi0PtGen(0),
fPi0PtEnh(0),
fEtaPtGen(0),
fEtaPtEnh(0),
fGenHfePt(0),
fGenPePt(0),
fPtP(0),
fptJetIE(0),
fptJetPE(0),
fptJetHFE(0),
fptRecPE(0),
fptTruePE(0),
fPtTrack(0),
fPhiTrack(0x0),
fEtaTrack(0x0),
fEtaPhiTrack(0x0),
fPhiRecElecTPC(0x0),
fEtaRecElecTPC(0x0),
fEtaPhiRecElecTPC(0x0),
fPhiRecElecEMCal(0x0),
fEtaRecElecEMCal(0x0),
fEtaPhiRecElecEMCal(0x0),
fPhiTrueElec(0x0),
fEtaTrueElec(0x0),
fEtaPhiTrueElec(0x0),
fnEovPelecNoTPCcut(0x0),
fnEovPelecTPCcut(0x0),
fnEovPelecTPCEMCalcut(0x0),
fnEovPbackg(0x0),
fnClsE(0x0),
fnM20(0x0),
fnM02(0x0),
fnClsTime(0x0),
fAngULS(0x0),
fAngLS(0x0),
fAngChargPart(0x0),
fAngHadron(0x0),
fAngIncElec(0x0),
fAngPhotElec(0x0),
fAngElecFromD(0x0),
fAngElecFromB(0x0),
fAngElecFromDFromB(0x0),
fAngD(0x0),
fAngB(0x0),
fAngCharm(0x0),
fAngBeauty(0x0),
fAngQuark(0x0),
fAngGluon(0x0),
fDispULS(0x0),
fDispLS(0x0),
fDispChargPart(0x0),
fDispHadron(0x0),
fDispIncElec(0x0),
fDispPhotElec(0x0),
fDispElecFromD(0x0),
fDispElecFromB(0x0),
fDispElecFromDFromB(0x0),
fDispD(0x0),
fDispB(0x0),
fDispCharm(0x0),
fDispBeauty(0x0),
fDispQuark(0x0),
fDispGluon(0x0),
fTreeObservableTagging(0)

{
    for(Int_t i=0;i<26;i++){
        fShapesVar[i]=0;
    }
    
    for(Int_t i=0;i<5;i++){
        fptTrueHFEeffTPCTOF[i] = NULL;
        fptTrueHFEeffEMCal[i] = NULL;
        fnEovPelecTPCsscut[i] = NULL;
    }
    
    for(Int_t i=0;i<2;i++){
        fptTrueHFEeffTPCTOFang[i] = NULL;
        fptTrueHFEeffEMCalang[i] = NULL;
        fptTrueHFEeffTPCTOFdisp[i] = NULL;
        fptTrueHFEeffEMCaldisp[i] = NULL;

    }
    
    
    for(Int_t i=0;i<5;i++){
        fInvmassLS[i] = NULL;
        fInvmassULS[i] = NULL;
        for(Int_t j=0;j<18;j++){
            fnTPCSigma[i][j] = NULL;
        }
    }
    
    for(Int_t i=0;i<5;i++){
        for(Int_t j=0;j<5;j++){
            fRecPEAng[i][j] = NULL;
            fTotPEAng[i][j] = NULL;
            fULSptAng[i][j] = NULL;
            fLSptAng[i][j] = NULL;
        }
        for(Int_t j=0;j<6;j++){
            fRecPEDisp[i][j] = NULL;
            fTotPEDisp[i][j] = NULL;
            fULSptDisp[i][j] = NULL;
            fLSptDisp[i][j] = NULL;
        }
    }
    
    
    SetMakeGeneralHistograms(kTRUE);
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalHfeTagging::AliAnalysisTaskEmcalHfeTagging(const char *name) :
AliAnalysisTaskEmcalJet(name, kTRUE),
fAOD(0),
fVevent(0),
fpidResponse(0),
fTracksTender(0),
pVtx(0),
spdVtx(0),
fMC(0),
fStack(0),
fMCparticle(0),
fMCarray(0),
fMCheader(0),
fContainer(0),
fMinFractionShared(0),
fJetShapeType(kData),
fJetShapeSub(kNoSub),
fJetSelection(kInclusive),
fPtThreshold(-9999.),
fRMatching(0.2),
fSelectedShapes(0),
fminpTTrig(20.),
fmaxpTTrig(50.),
fangWindowRecoil(0.6),
fSemigoodCorrect(0),
fHolePos(0),
fHoleWidth(0),
fCentSelectOn(kTRUE),
fCentMin(0),
fCentMax(10),
fOneConstSelectOn(kFALSE),
fDerivSubtrOrder(0),
fMCweight(0),
fRunNumber(0),
fAssPtCut(0.1),
fITSncut(3),
fAssTPCnCut(60),
fTPCnCut(100),
fAssITSrefitCut(kTRUE),
fSigmaTOFcut(3.),
fSigmaTPCcutLowPt(-1.),
fSigmaTPCcutHighPt(-1.5),
fSigmTPCcutExcElec(3.5),
fDcaXYcut(1.),
fDcaZcut(2.),
fIMcut(0.1),
fEtaCut(0.7),
fMinEoPcut(0.9),
fMaxEoPcut(1.3),
fM20cut(0.35),
fMinPtTPC(0.5),
fMaxPtTPC(4.),
fMinPtEMCal(4.),
fMaxPtEMCal(25.),
fNeventV0(0),
fNeventT0(0),
fh2ResponseUW(0x0),
fh2ResponseW(0x0),
fPhiJetCorr6(0x0),
fPhiJetCorr7(0x0),
fEtaJetCorr6(0x0),
fEtaJetCorr7(0x0),
fPtJetCorr(0x0),
fPtJet(0x0),
fPtGenJet(0),
fPhiJet(0x0),
fEtaJet(0x0),
fEtaPhiJet(0x0),
fAreaJet(0x0),
fJetProbDensityDetPart(0x0),
fJetProbDensityPartDet(0x0),
fNbOfConstvspT(0x0),
fnTPCnTOFnocut(0),
fnTPCnocutP(0),
fnTOFnocutP(0),
fnTPCcutP(0),
fnTPCcutPt(0),
fnULSmLSpairsPerElectron(0),
fnPartPerJet(0),
fnElecOverPartPerJet(0),
fnInclElecPerJet(0),
fnPhotElecPerJet(0),
fnIncSubPhotElecPerJet(0),
fnTrueElecPerJet(0),
fnTrueHFElecPerJet(0),
fnTruePElecPerJet(0),
fPi0PtGen(0),
fPi0PtEnh(0),
fEtaPtGen(0),
fEtaPtEnh(0),
fGenHfePt(0),
fGenPePt(0),
fPtP(0),
fptJetIE(0),
fptJetPE(0),
fptJetHFE(0),
fptRecPE(0),
fptTruePE(0),
fPtTrack(0),
fPhiTrack(0x0),
fEtaTrack(0x0),
fEtaPhiTrack(0x0),
fPhiRecElecTPC(0x0),
fEtaRecElecTPC(0x0),
fEtaPhiRecElecTPC(0x0),
fPhiRecElecEMCal(0x0),
fEtaRecElecEMCal(0x0),
fEtaPhiRecElecEMCal(0x0),
fPhiTrueElec(0x0),
fEtaTrueElec(0x0),
fEtaPhiTrueElec(0x0),
fnEovPelecNoTPCcut(0x0),
fnEovPelecTPCcut(0x0),
fnEovPelecTPCEMCalcut(0x0),
fnEovPbackg(0x0),
fnClsE(0x0),
fnM20(0x0),
fnM02(0x0),
fnClsTime(0x0),
fAngULS(0x0),
fAngLS(0x0),
fAngChargPart(0x0),
fAngHadron(0x0),
fAngIncElec(0x0),
fAngPhotElec(0x0),
fAngElecFromD(0x0),
fAngElecFromB(0x0),
fAngElecFromDFromB(0x0),
fAngD(0x0),
fAngB(0x0),
fAngCharm(0x0),
fAngBeauty(0x0),
fAngQuark(0x0),
fAngGluon(0x0),
fDispULS(0x0),
fDispLS(0x0),
fDispChargPart(0x0),
fDispHadron(0x0),
fDispIncElec(0x0),
fDispPhotElec(0x0),
fDispElecFromD(0x0),
fDispElecFromB(0x0),
fDispElecFromDFromB(0x0),
fDispD(0x0),
fDispB(0x0),
fDispCharm(0x0),
fDispBeauty(0x0),
fDispQuark(0x0),
fDispGluon(0x0),
fTreeObservableTagging(0)
{
    // Standard constructor.
    for(Int_t i=0;i<26;i++){
        fShapesVar[i]=0;
    }
    
    for(Int_t i=0;i<5;i++){
        fptTrueHFEeffTPCTOF[i] = NULL;
        fptTrueHFEeffEMCal[i] = NULL;
        fnEovPelecTPCsscut[i] = NULL;
    }
    
    for(Int_t i=0;i<2;i++){
        fptTrueHFEeffTPCTOFang[i] = NULL;
        fptTrueHFEeffEMCalang[i] = NULL;
        fptTrueHFEeffTPCTOFdisp[i] = NULL;
        fptTrueHFEeffEMCaldisp[i] = NULL;
    }
    
    for(Int_t i=0;i<5;i++){
        fInvmassLS[i] = NULL;
        fInvmassULS[i] = NULL;
        for(Int_t j=0;j<18;j++){
            fnTPCSigma[i][j] = NULL;
        }
    }
    
    for(Int_t i=0;i<5;i++){
        for(Int_t j=0;j<5;j++){
            fRecPEAng[i][j] = NULL;
            fTotPEAng[i][j] = NULL;
            fULSptAng[i][j] = NULL;
            fLSptAng[i][j] = NULL;
        }
        for(Int_t j=0;j<6;j++){
            fRecPEDisp[i][j] = NULL;
            fTotPEDisp[i][j] = NULL;
            fULSptDisp[i][j] = NULL;
            fLSptDisp[i][j] = NULL;
        }
    }
    
    SetMakeGeneralHistograms(kTRUE);
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalHfeTagging::~AliAnalysisTaskEmcalHfeTagging()
{
    // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHfeTagging::UserCreateOutputObjects()
{
    // Create user output.
    
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    
    Double_t ptRange[34] = {0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4,
        1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 5.,
        6., 8., 10., 12., 14., 16., 19., 22., 25., 30.,
        35., 40., 45., 50.};
    
    int nbin = 59;
    double xbins[60] =  {0.01,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,
        0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,
        0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,
        1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,
        3.6,3.8,4,4.5,5,5.5,6,6.5,7,8,
        9,10,11,12,13,14,15,16,18,20};
    
    Double_t bin_JetPt[6] = {5.,20.,40.,60.,80.,120.};
    Double_t bin_g[9] = {0.,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.12};//g
    Double_t bin_ptd[7] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};//pTD
    
    Double_t angRange[6] = {0.,0.02,0.04,0.06,0.08,0.12};//for syst.
    
    fNeventV0 = new TH1F("fNeventV0","Number of Events (V0)",5,-0.5,4.5);
    fOutput->Add(fNeventV0);
    
    fNeventT0 = new TH1F("fNeventT0","Number of Events (T0)",5,-0.5,4.5);
    fOutput->Add(fNeventT0);
    
    fh2ResponseUW= new TH2F("fh2ResponseUW", "fh2ResponseUW", 100, 0, 200,  100, 0, 200);
    fOutput->Add(fh2ResponseUW);
    
    fh2ResponseW= new TH2F("fh2ResponseW", "fh2ResponseW", 100, 0, 200,  100, 0, 200);
    fOutput->Add(fh2ResponseW);
    
    fPhiJetCorr6= new TH2F("fPhiJetCorr6", "fPhiJetCorr6", 50, 0, 2*TMath::Pi(), 50, 0, 2*TMath::Pi());
    fOutput->Add(fPhiJetCorr6);
    
    fEtaJetCorr6= new TH2F("fEtaJetCorr6", "fEtaJetCorr6", 50, -1.5, 1.5, 50, -1.5, 1.5);
    fOutput->Add(fEtaJetCorr6);
    
    fPhiJetCorr7= new TH2F("fPhiJetCorr7", "fPhiJetCorr7", 50, 0, 2*TMath::Pi(), 50, 0, 2*TMath::Pi());
    fOutput->Add(fPhiJetCorr7);
    
    fEtaJetCorr7= new TH2F("fEtaJetCorr7", "fEtaJetCorr7", 50, -1.5, 1.5, 50, -1.5, 1.5);
    fOutput->Add(fEtaJetCorr7);
    
    fPtJetCorr= new TH2F("fPtJetCorr", "fPtJetCorr", 100, 0, 200,  100, 0, 200);
    fOutput->Add(fPtJetCorr);
    
    fPtJet= new TH1F("fPtJet", "fPtJet", 100, 0, 200);
    fOutput->Add(fPtJet);
    
    fPtGenJet= new TH1F("fPtGenJet", "fPtGenJet", 100, 0, 200);
    fOutput->Add(fPtGenJet);
    
    fPhiJet= new TH2F("fPhiJet", "fPhiJet", 100, 0, 200, 100, 0, TMath::TwoPi());
    fOutput->Add(fPhiJet);
    
    fEtaJet= new TH2F("fEtaJet", "fEtaJet", 100, 0, 200, 100, -1,1);
    fOutput->Add(fEtaJet);
    
    fEtaPhiJet= new TH2F("fEtaPhiJet", "fEtaPhiJet", 100, 0, TMath::TwoPi(), 100, -1,1);
    fOutput->Add(fEtaPhiJet);
    
    fAreaJet= new TH2F("fAreaJet", "fAreaJet", 100, 0, 200, 100, 0,1.5);
    fOutput->Add(fAreaJet);
    
    fJetProbDensityDetPart=new TH2F("fJetProbDensityDetPart", "fJetProbDensityDetPart",100,-2,2,100, 0,200);
    fOutput->Add(fJetProbDensityDetPart);
    
    fJetProbDensityPartDet=new TH2F("fJetProbDensityPartDet", "fJetProbDensityPartDet",100,-2,2,100, 0,200);
    fOutput->Add(fJetProbDensityPartDet);
    
    fNbOfConstvspT=new TH2F("fNbOfConstvspT", "fNbOfConstvspT", 100, 0, 100, 200, 0, 200);
    fOutput->Add(fNbOfConstvspT);
    
    fnTPCnTOFnocut=new TH2F("fnTPCnTOFnocut", "fnTPCnTOFnocut", 100, -10, 10, 100, -10, 10);
    fOutput->Add(fnTPCnTOFnocut);
    
    fnTPCnocutP=new TH2F("fnTPCnocutP", "fnTPCnocutP", 50, 0, 5, 200, -10, 10);
    fOutput->Add(fnTPCnocutP);
    
    fnTOFnocutP=new TH2F("fnTOFnocutP", "fnTOFnocutP", 50, 0, 5, 200, -10, 10);
    fOutput->Add(fnTOFnocutP);
    
    fnTPCcutP=new TH2F("fnTPCcutP", "fnTPCcutP", 50, 0, 5, 200, -10, 10);
    fOutput->Add(fnTPCcutP);
    
    fnTPCcutPt=new TH2F("fnTPCcutPt", "fnTPCcutPt", 50, 0, 5, 200, -10, 10);
    fOutput->Add(fnTPCcutPt);
    
    for(Int_t i=0;i<5;i++){
        for(Int_t j=0;j<18;j++){
            fnTPCSigma[i][j] = new TH1F(Form("fnTPCSigma%d%d",i,j),Form("fnTPCSigma%d%d",i,j), 100,-15,15);
            fOutput->Add(fnTPCSigma[i][j]);
        }
    }
    
    fnULSmLSpairsPerElectron =new TH2F("fnULSmLSpairsPerElectron", "fnULSmLSpairsPerElectron", 50, 0, 5, 100, 0, 100);
    fOutput->Add(fnULSmLSpairsPerElectron);
    
    for(Int_t i=0;i<5;i++){
        fInvmassLS[i] = new TH2F(Form("fInvmassLS%d",i), Form("fInvmassLS%d",i), 33,ptRange, 100, 0, 0.5);
        fOutput->Add(fInvmassLS[i]);
    }
    
    for(Int_t i=0;i<5;i++){
        fInvmassULS[i] = new TH2F(Form("fInvmassULS%d",i), Form("fInvmassULS%d",i), 33,ptRange, 100, 0, 0.5);
        fOutput->Add(fInvmassULS[i]);
    }
    
    fnPartPerJet=new TH1F("fnPartPerJet", "fnPartPerJet", 50,0,50);
    fOutput->Add(fnPartPerJet);
    
    fnElecOverPartPerJet=new TH1F("fnElecOverPartPerJet", "fnElecOverPartPerJet", 50,0,0.1);
    fOutput->Add(fnElecOverPartPerJet);
    
    fnInclElecPerJet=new TH1F("fnInclElecPerJet", "fnInclElecPerJet", 50,0,50);
    fOutput->Add(fnInclElecPerJet);
    
    fnPhotElecPerJet=new TH1F("fnPhotElecPerJet", "fnPhotElecPerJet", 50,0,50);
    fOutput->Add(fnPhotElecPerJet);
    
    fnIncSubPhotElecPerJet=new TH1F("fnIncSubPhotElecPerJet", "fnIncSubPhotElecPerJet", 51,-1,50);
    fOutput->Add(fnIncSubPhotElecPerJet);
    
    fnTrueElecPerJet=new TH1F("fnTrueElecPerJet", "fnTrueElecPerJet", 50,0,50);
    fOutput->Add(fnTrueElecPerJet);
    
    fnTrueHFElecPerJet=new TH1F("fnTrueHFElecPerJet", "fnTrueHFElecPerJet", 50,0,50);
    fOutput->Add(fnTrueHFElecPerJet);
    
    fnTruePElecPerJet=new TH1F("fnTruePElecPerJet", "fnTruePElecPerJet", 50,0,50);
    fOutput->Add(fnTruePElecPerJet);
    
    fPi0PtGen = new TH1F("fPi0PtGen","fPi0PtGen",nbin,xbins);
    fOutput->Add(fPi0PtGen);
    
    fPi0PtEnh = new TH1F("fPi0PtEnh","fPi0PtEnh",nbin,xbins);
    fOutput->Add(fPi0PtEnh);
    
    fEtaPtGen = new TH1F("fEtaPtGen", "fEtaPtGen",nbin,xbins);
    fOutput->Add(fEtaPtGen);
    
    fEtaPtEnh = new TH1F("fEtaPtEnh", "fEtaPtEnh",nbin,xbins);
    fOutput->Add(fEtaPtEnh);
    
    fGenHfePt = new TH1F("fGenHfePt","fGenHfePt",nbin,xbins);
    fOutput->Add(fGenHfePt);
    
    fGenPePt = new TH1F("fGenPePt", "fGenPePt",nbin,xbins);
    fOutput->Add(fGenPePt);
    
    for(Int_t i=0;i<5;i++){
        
        for(Int_t j=0;j<5;j++){
            fRecPEAng[i][j] = new TH1F(Form("fRecPEAng%d%d",i,j), Form("fRecPEAng%d%d",i,j), 33,ptRange);
            fOutput->Add(fRecPEAng[i][j]);

            fTotPEAng[i][j] = new TH1F(Form("fTotPEAng%d%d",i,j), Form("fTotPEAng%d%d",i,j), 33,ptRange);
            fOutput->Add(fTotPEAng[i][j]);
            
            fULSptAng[i][j] = new TH1F(Form("fULSptAng%d%d",i,j), Form("fULSptAng%d%d",i,j), 33,ptRange);
            fOutput->Add(fULSptAng[i][j]);

            fLSptAng[i][j] = new TH1F(Form("fLSptAng%d%d",i,j), Form("fLSptAng%d%d",i,j), 33,ptRange);
            fOutput->Add(fLSptAng[i][j]);
        }
        
        for(Int_t j=0;j<6;j++){
            fRecPEDisp[i][j] = new TH1F(Form("fRecPEDisp%d%d",i,j), Form("fRecPEDisp%d%d",i,j), 33,ptRange);
            fOutput->Add(fRecPEDisp[i][j]);
            
            fTotPEDisp[i][j] = new TH1F(Form("fTotPEDisp%d%d",i,j), Form("fTotPEDisp%d%d",i,j), 33,ptRange);
            fOutput->Add(fTotPEDisp[i][j]);
            
            fULSptDisp[i][j] = new TH1F(Form("fULSptDisp%d%d",i,j), Form("fULSptDisp%d%d",i,j), 33,ptRange);
            fOutput->Add(fULSptDisp[i][j]);
            
            fLSptDisp[i][j] = new TH1F(Form("fLSptDisp%d%d",i,j), Form("fLSptDisp%d%d",i,j), 33,ptRange);
            fOutput->Add(fLSptDisp[i][j]);
        }
    }
    

    
    fPtP=new TH2F("fPtP", "fPtP", 33,ptRange,33,ptRange);
    fOutput->Add(fPtP);
    
    fptJetIE= new TH1F("fptJetIE", "fptJetIE", 100, 0, 200);
    fOutput->Add(fptJetIE);
    
    fptJetPE= new TH1F("fptJetPE", "fptJetPE", 100, 0, 200);
    fOutput->Add(fptJetPE);
    
    fptJetHFE= new TH1F("fptJetHFE", "fptJetHFE", 100, 0, 200);
    fOutput->Add(fptJetHFE);
    
    fptRecPE= new TH1F("fptRecPE", "fptRecPE", 33,ptRange);
    fOutput->Add(fptRecPE);
    
    fptTruePE= new TH1F("fptTruePE", "fptTruePE", 33,ptRange);
    fOutput->Add(fptTruePE);
    
    for(Int_t i=0;i<5;i++){
        fptTrueHFEeffTPCTOF[i] = new TH1F(Form("fptTrueHFEeffTPCTOF%d",i), Form("fptTrueHFEeffTPCTOF%d",i), 33,ptRange);
        fOutput->Add(fptTrueHFEeffTPCTOF[i]);
    }
    
    for(Int_t i=0;i<2;i++){
        fptTrueHFEeffTPCTOFang[i] = new TH3F(Form("fptTrueHFEeffTPCTOFang%d",i),
                                             Form("fptTrueHFEeffTPCTOFang%d",i), 33,ptRange,5, bin_JetPt, 5, angRange);
        fOutput->Add(fptTrueHFEeffTPCTOFang[i]);
    }
    
    for(Int_t i=0;i<2;i++){
        fptTrueHFEeffTPCTOFdisp[i] = new TH3F(Form("fptTrueHFEeffTPCTOFdisp%d",i),
                                              Form("fptTrueHFEeffTPCTOFdisp%d",i), 33,ptRange,5, bin_JetPt,6, bin_ptd);
        fOutput->Add(fptTrueHFEeffTPCTOFdisp[i]);
    }
    
    
    for(Int_t i=0;i<5;i++){
        fptTrueHFEeffEMCal[i] = new TH1F(Form("fptTrueHFEeffEMCal%d",i), Form("fptTrueHFEeffEMCal%d",i), 33,ptRange);
        fOutput->Add(fptTrueHFEeffEMCal[i]);
    }
    
    for(Int_t i=0;i<2;i++){
        fptTrueHFEeffEMCalang[i] = new TH3F(Form("fptTrueHFEeffEMCalang%d",i),
                                            Form("fptTrueHFEeffEMCalang%d",i), 33,ptRange,5, bin_JetPt, 5, angRange);
        fOutput->Add(fptTrueHFEeffEMCalang[i]);
    }
    
    for(Int_t i=0;i<2;i++){
        fptTrueHFEeffEMCaldisp[i] = new TH3F(Form("fptTrueHFEeffEMCaldisp%d",i),
                                             Form("fptTrueHFEeffEMCaldisp%d",i), 33,ptRange,5, bin_JetPt,6, bin_ptd);
        fOutput->Add(fptTrueHFEeffEMCaldisp[i]);
    }
    
    fPtTrack= new TH1F("fPtTrack", "fPtTrack", 100, 0, 200);
    fOutput->Add(fPtTrack);
    
    fPhiTrack= new TH2F("fPhiTrack", "fPhiTrack", 100, 0, 200, 100, 0, TMath::TwoPi());
    fOutput->Add(fPhiTrack);
    
    fEtaTrack= new TH2F("fEtaTrack", "fEtaTrack", 100, 0, 200, 100, -1,1);
    fOutput->Add(fEtaTrack);
    
    fEtaPhiTrack= new TH2F("fEtaPhiTrack", "fEtaPhiTrack", 100, 0, TMath::TwoPi(), 100, -1,1);
    fOutput->Add(fEtaPhiTrack);
    
    fPhiRecElecTPC= new TH2F("fPhiRecElecTPC", "fPhiRecElecTPC", 100, 0, 50, 100, 0, TMath::TwoPi());
    fOutput->Add(fPhiRecElecTPC);
    
    fEtaRecElecTPC= new TH2F("fEtaRecElecTPC", "fEtaRecElecTPC", 100, 0, 50, 100, -1,1);
    fOutput->Add(fEtaRecElecTPC);
    
    fEtaPhiRecElecTPC= new TH2F("fEtaPhiRecElecTPC", "fEtaPhiRecElecTPC", 100, 0, TMath::TwoPi(), 100, -1,1);
    fOutput->Add(fEtaPhiRecElecTPC);
    
    fPhiRecElecEMCal= new TH2F("fPhiRecElecEMCal", "fPhiRecElecEMCal", 100, 0, 50, 100, 0, TMath::TwoPi());
    fOutput->Add(fPhiRecElecEMCal);
    
    fEtaRecElecEMCal= new TH2F("fEtaRecElecEMCal", "fEtaRecElecEMCal", 100, 0, 50, 100, -1,1);
    fOutput->Add(fEtaRecElecEMCal);
    
    fEtaPhiRecElecEMCal= new TH2F("fEtaPhiRecElecEMCal", "fEtaPhiRecElecEMCal", 100, 0, TMath::TwoPi(), 100, -1,1);
    fOutput->Add(fEtaPhiRecElecEMCal);
    
    fPhiTrueElec= new TH2F("fPhiTrueElec", "fPhiTrueElec", 100, 0, 50, 100, 0, TMath::TwoPi());
    fOutput->Add(fPhiTrueElec);
    
    fEtaTrueElec= new TH2F("fEtaTrueElec", "fEtaTrueElec", 100, 0, 50, 100, -1,1);
    fOutput->Add(fEtaTrueElec);
    
    fEtaPhiTrueElec= new TH2F("fEtaPhiTrueElec", "fEtaPhiTrueElec", 100, 0, TMath::TwoPi(), 100, -1,1);
    fOutput->Add(fEtaPhiTrueElec);
    
    fnEovPelecNoTPCcut = new TH2F("fnEovPelecNoTPCcut", "fnEovPelecNoTPCcut", 100, 0, 100, 40,0,2);
    fOutput->Add(fnEovPelecNoTPCcut);
    
    fnEovPelecTPCcut = new TH2F("fnEovPelecTPCcut", "fnEovPelecTPCcut", 100, 0, 100, 40,0,2);
    fOutput->Add(fnEovPelecTPCcut);
    
    fnEovPelecTPCEMCalcut = new TH2F("fnEovPelecTPCEMCalcut", "fnEovPelecTPCEMCalcut", 100, 0, 100,40,0,2);
    fOutput->Add(fnEovPelecTPCEMCalcut);
    
    for(Int_t i=0;i<5;i++){
        fnEovPelecTPCsscut[i] = new TH2F(Form("fnEovPelecTPCsscut%d",i), Form("fnEovPelecTPCsscut%d",i), 33,ptRange,40,0,2);
        fOutput->Add(fnEovPelecTPCsscut[i]);
    }
    
    fnEovPbackg = new TH2F("fnEovPbackg", "fnEovPbackg", 100, 0, 100,40,0,2);
    fOutput->Add(fnEovPbackg);
    
    fnClsE = new TH2F("fnClsE", "fnClsE", 100, 0, 100, 100, 0,100);
    fOutput->Add(fnClsE);
    
    fnM20 = new TH2F("fnM20", "fnM20", 100, 0, 100, 100, 0,2);
    fOutput->Add(fnM20);
    
    fnM02 = new TH2F("fnM02", "fnM02", 100, 0, 100, 100, 0,2);
    fOutput->Add(fnM02);
    
    fnClsTime = new TH2F("fnClsTime", "fnClsTime", 100, 0, 100, 100, -200,200);
    fOutput->Add(fnClsTime);
    
    fAngULS = new TH2F("fAngULS", "fAngULS", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngULS);
    
    fAngLS = new TH2F("fAngLS", "fAngLS", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngLS);
    
    fAngChargPart = new TH2F("fAngChargPart", "fAngChargPart", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngChargPart);
    
    fAngHadron = new TH2F("fAngHadron", "fAngHadron", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngHadron);
    
    fAngIncElec = new TH2F("fAngIncElec", "fAngIncElec", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngIncElec);
    
    fAngPhotElec = new TH2F("fAngPhotElec", "fAngPhotElec", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngPhotElec);
    
    fAngElecFromD = new TH2F("fAngElecFromD", "fAngElecFromD", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngElecFromD);
    
    fAngElecFromB = new TH2F("fAngElecFromB", "fAngElecFromB", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngElecFromB);
    
    fAngElecFromDFromB = new TH2F("fAngElecFromDFromB", "fAngElecFromDFromB", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngElecFromDFromB);
    
    fAngD = new TH2F("fAngD", "fAngD", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngD);
    
    fAngB = new TH2F("fAngB", "fAngB", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngB);
    
    fAngCharm = new TH2F("fAngCharm", "fAngCharm", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngCharm);
    
    fAngBeauty = new TH2F("fAngBeauty", "fAngBeauty", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngBeauty);
    
    fAngQuark = new TH2F("fAngQuark", "fAngQuark", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngQuark);
    
    fAngGluon = new TH2F("fAngGluon", "fAngGluon", 5, bin_JetPt, 8, bin_g);
    fOutput->Add(fAngGluon);
    
    fDispULS = new TH2F("fDispULS", "fDispULS", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispULS);
    
    fDispLS = new TH2F("fDispLS", "fDispLS", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispLS);
    
    fDispChargPart = new TH2F("fDispChargPart", "fDispChargPart", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispChargPart);
    
    fDispHadron = new TH2F("fDispHadron", "fDispHadron", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispHadron);
    
    fDispIncElec = new TH2F("fDispIncElec", "fDispIncElec", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispIncElec);
    
    fDispPhotElec = new TH2F("fDispPhotElec", "fDispPhotElec", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispPhotElec);
    
    fDispElecFromD = new TH2F("fDispElecFromD", "fDispElecFromD", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispElecFromD);
    
    fDispElecFromB = new TH2F("fDispElecFromB", "fDispElecFromB", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispElecFromB);
    
    fDispElecFromDFromB = new TH2F("fDispElecFromDFromB", "fDispElecFromDFromB", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispElecFromDFromB);
    
    fDispD = new TH2F("fDispD", "fDispD", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispD);
    
    fDispB = new TH2F("fDispB", "fDispB", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispB);
    
    fDispCharm = new TH2F("fDispCharm", "fDispCharm", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispCharm);
    
    fDispBeauty = new TH2F("fDispBeauty", "fDispBeauty", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispBeauty);
    
    fDispQuark = new TH2F("fDispQuark", "fDispQuark", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispQuark);
    
    fDispGluon = new TH2F("fDispGluon", "fDispGluon", 5, bin_JetPt, 6, bin_ptd);
    fOutput->Add(fDispGluon);
    
    // =========== Switch on Sumw2 for all histos ===========
    for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
        TH1 *h = dynamic_cast<TH1*>(fOutput->At(i));
        if (h){
            h->Sumw2();
            continue;
        }
        THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
        if(hn)hn->Sumw2();
    }
    
    
    TH1::AddDirectory(oldStatus);
    const Int_t nVar = 26;
    
    fTreeObservableTagging = new TTree(GetOutputSlot(2)->GetContainer()->GetName(), GetOutputSlot(2)->GetContainer()->GetName());
    
    TString *fShapesVarNames = new TString [nVar];
    
    fShapesVarNames[0] = "partonCode";
    fShapesVarNames[1] = "ptJet";
    fShapesVarNames[2] = "ptDJet";
    fShapesVarNames[3] = "mJet";
    // fShapesVarNames[4] = "nbOfConst";
    fShapesVarNames[4] = "angularity";
    fShapesVarNames[5] = "circularity";
    fShapesVarNames[6] = "lesub";
    fShapesVarNames[7] = "coronna";
    
    fShapesVarNames[8] = "ptJetMatch";
    fShapesVarNames[9] = "ptDJetMatch";
    fShapesVarNames[10] = "mJetMatch";
    // fShapesVarNames[12] = "nbOfConstMatch";
    fShapesVarNames[11] = "angularityMatch";
    fShapesVarNames[12] = "circularityMatch";
    fShapesVarNames[13] = "lesubMatch";
    fShapesVarNames[14] = "coronnaMatch";
    fShapesVarNames[15] = "weightPythia";
    //fShapesVarNames[14]="ntrksEvt";
    fShapesVarNames[16] = "rhoVal";
    fShapesVarNames[17] = "rhoMassVal";
    fShapesVarNames[18] = "ptUnsub";
    fShapesVarNames[19] = "ptTrueHFE";
    fShapesVarNames[20] = "ptElec";
    fShapesVarNames[21] = "nInclElec";
    fShapesVarNames[22] = "nPhotElec";
    fShapesVarNames[23] = "hasElec";
    fShapesVarNames[24] = "nTrueElectronsMC";
    fShapesVarNames[25] = "nTrueHFElecMC";
    
    
    for(Int_t ivar=0; ivar < nVar; ivar++){
        //cout<<"looping over variables"<<endl;
        fTreeObservableTagging->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));
    }
    
    PostData(1,fOutput);
    PostData(2,fTreeObservableTagging);
    
    delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::Run()
{
    //  // Run analysis code here, if needed. It will be executed before FillHistograms().
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD){
        printf("ERROR: fAOD not available\n");
        return kFALSE;
    }
    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if(!fVevent){
        printf("ERROR: fVEvent not available\n");
        return kFALSE;
    }
    
    fRunNumber = fVevent->GetRunNumber();
    
    //remaining event selection
    pVtx = fVevent->GetPrimaryVertex();
    spdVtx = fVevent->GetPrimaryVertexSPD();
    
    fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    
    fNeventV0->Fill(0);
    if (!fMCheader){
        TString firedTriggerClasses = static_cast<const AliAODEvent*>(InputEvent())->GetFiredTriggerClasses();
        if ((firedTriggerClasses.Contains("C0TVX-B-NOPF-CENT"))){
            fNeventT0->Fill(0);
        }
    }
    
    Int_t NpureMC = -1, NpureMCproc = -1;
    
    if (fMCheader){
        TList *lh=fMCheader->GetCocktailHeaders();
        if(lh){
            for(int igene=0; igene<lh->GetEntries(); igene++){
                AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
                if(gh){
                    //  cout << "<------- imc = "<< gh->GetName() << endl;
                    if(igene==0)NpureMC = gh->NProduced();  // generate by PYTHIA or HIJING
                    NpureMCproc += gh->NProduced();
                }
            }
        }
    }
    // pileup rejection with SPD multiple vertices and track multivertexer (using default configuration)
    
    Bool_t isPileupfromSPDmulbins = kFALSE, isPileupFromMV = kFALSE;
    isPileupfromSPDmulbins = fAOD->IsPileupFromSPDInMultBins();
    if (isPileupfromSPDmulbins) return kFALSE;
    fNeventV0->Fill(1);
    
    AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(5);
    utils.SetMaxPlpChi2MV(5.);
    utils.SetMinWDistMV(15);
    utils.SetCheckPlpFromDifferentBCMV(kFALSE);
    isPileupFromMV = utils.IsPileUpMV(fAOD);
    if (isPileupFromMV) return kFALSE;
    fNeventV0->Fill(2);
    
    // Selection of pi0 and eta in MC to compute the weight
    
    Bool_t isPrimary = kFALSE, isFromLMdecay = kTRUE, isFromHFdecay=kTRUE, hasMother = kTRUE;
    
    Double_t MCweight = 1.;
    Int_t iDecay = 0;
    
    if(fMCarray){
        //Int_t nParticles = fMCarray->GetEntries();
        for (Int_t iParticle = 0; iParticle < NpureMCproc; iParticle++) {
            AliAODMCParticle* particle = (AliAODMCParticle*) fMCarray->At(iParticle);
            int fPDG = particle->GetPdgCode();
            Double_t pTMC = particle->Pt();
            Double_t phiMC = particle->Phi();
            Double_t etaMC = particle->Eta();
            
            Bool_t iEnhance = kFALSE;
            if(iParticle>=NpureMC)iEnhance = kTRUE;
            
            //if (TMath::Abs(etaMC)>1.2)continue;
            
            isPrimary = IsPrimary(particle);
            isFromLMdecay = IsFromLMdecay(particle);
            isFromHFdecay = IsFromHFdecay(particle);
            hasMother = HasMother(particle);
            
            MCweight = 1.;
            iDecay = 0;
            
            GetWeightAndDecay(particle,iDecay,MCweight);
            
            
            if (TMath::Abs(etaMC)<0.8 && TMath::Abs(fPDG)==11){
                fPhiTrueElec->Fill(pTMC,phiMC);
                fEtaTrueElec->Fill(pTMC,etaMC);
                if (pTMC > 0.5) fEtaPhiTrueElec->Fill(phiMC,etaMC);
                
                if (isFromHFdecay) fGenHfePt->Fill(pTMC);
                if (iDecay>0 && iDecay<7) fGenPePt->Fill(pTMC);
            }
            
            
            Double_t yMC = particle->Y();
            if (TMath::Abs(yMC)>1.0)continue;
            
            if (isPrimary){
                if (!hasMother || (!isFromLMdecay && !isFromHFdecay)){
                    if(fPDG==111 && iEnhance==kFALSE) fPi0PtGen->Fill(pTMC);
                    if(fPDG==111 && iEnhance==kTRUE)  fPi0PtEnh->Fill(pTMC);
                    
                    if(fPDG==221 && iEnhance==kFALSE) fEtaPtGen->Fill(pTMC);
                    if(fPDG==221 && iEnhance==kTRUE)  fEtaPtEnh->Fill(pTMC);
                }
            }
        }
    }//MC
    return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::FillHistograms()
{
    // Fill histograms.
    //cout<<"base container"<<endl;
    AliEmcalJet* jet1 = NULL;
    AliJetContainer *jetCont = GetJetContainer(0);
    
    Float_t kWeight=1;
    if (fCentSelectOn)
        if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
    
    // list of kink mothers
    Int_t nVertices = 1;
    nVertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[nVertices];
    Int_t nMotherKink = 0;
    for(Int_t ivertex=0; ivertex < nVertices; ivertex++) {
        AliAODVertex *vertex = fAOD->GetVertex(ivertex);
        if(!vertex) continue;
        if(vertex->GetType()==AliAODVertex::kKink) {
            AliAODTrack *mother = (AliAODTrack *) vertex->GetParent();
            if(!mother) continue;
            Int_t idmother = mother->GetID();
            listofmotherkink[nMotherKink] = idmother;
            nMotherKink++;
        }
    }
    
    Float_t rhoVal=0, rhoMassVal = 0.;
    
    if(jetCont) {
        jetCont->ResetCurrentID();
        if ((fJetShapeSub==kConstSub) || (fJetShapeSub==kDerivSub)){
            //rho
            AliRhoParameter* rhoParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoSparseR020"));
            if (!rhoParam) {
                Printf("%s: Could not retrieve rho %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoName().Data());
            } else rhoVal = rhoParam->GetVal();
            //rhom
            AliRhoParameter* rhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoMassSparseR020"));
            if (!rhomParam) {
                Printf("%s: Could not retrieve rho_m %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoMassName().Data());
            } else rhoMassVal = rhomParam->GetVal();
        }
        
        while((jet1 = jetCont->GetNextAcceptJet())) {
            if (!jet1) continue;
            AliEmcalJet* jet2 = 0x0;
            AliEmcalJet* jet3 = 0x0;
            fPtJet->Fill(jet1->Pt());
            fPhiJet->Fill(jet1->Pt(),jet1->Phi());
            fEtaJet->Fill(jet1->Pt(),jet1->Eta());
            if (jet1->Pt() > 5.) fEtaPhiJet->Fill(jet1->Phi(),jet1->Eta());
            fAreaJet->Fill(jet1->Pt(),jet1->Area());
            AliEmcalJet *jetUS = NULL;
            Int_t ifound=0, ilab=-1;
            
            fShapesVar[0] = 0.;
            if(fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kDetEmbPartPythia){
                //AliJetContainer *jetContTrue = GetJetContainer(1);
                AliJetContainer *jetContUS = GetJetContainer(2);
                
                if(fJetShapeSub == AliAnalysisTaskEmcalHfeTagging::kConstSub){
                    for(Int_t i = 0; i<jetContUS->GetNJets(); i++) {
                        jetUS = jetContUS->GetJet(i);
                        if(jetUS->GetLabel()==jet1->GetLabel()) {
                            ifound++;
                            if(ifound==1) ilab = i;
                        }
                    }
                    if(ilab==-1) continue;
                    jetUS=jetContUS->GetJet(ilab);
                    jet2=jetUS->ClosestJet();
                }
                
                if(!(fJetShapeSub == AliAnalysisTaskEmcalHfeTagging::kConstSub)) jet2 = jet1->ClosestJet();
                if (!jet2) {
                    Printf("jet2 does not exist, returning");
                    continue;
                }
                
                //AliJetContainer *jetContPart=GetJetContainer(3);
                jet3=jet2->ClosestJet();
                
                if(!jet3){
                    Printf("jet3 does not exist, returning");
                    continue;
                }
                //cout<<"jet 3 exists"<<jet3->Pt()<<endl;
                
                
                fh2ResponseUW->Fill(jet1->Pt(),jet2->Pt());
                
                Double_t fraction=0;
                if(!(fJetShapeSub==kConstSub))  fraction = jetCont->GetFractionSharedPt(jet1);
                if(fJetShapeSub==kConstSub) fraction = jetContUS->GetFractionSharedPt(jetUS);
                //if (fraction > 0.1) cout<<"***** hey a jet matched with fraction"<<fraction<<"  "<<jet1->Pt()<<" "<<jet2->Pt()<<" "<<fCent<<endl;
                
                if(fraction<fMinFractionShared) continue;
                //InputEvent()->Print();
                
            }
            
            if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kPythiaDef){
                jet3 = jet1->ClosestJet();
                if (!jet3) continue;
                
                fh2ResponseUW->Fill(jet1->Pt(),jet3->Pt());
                
                Double_t probDensDetPart = -999., probDensPartDet = -999.;
                
                if (jet1->Pt()>0) probDensPartDet = (jet3->Pt()-jet1->Pt())/jet1->Pt();
                if (jet3->Pt()>0) probDensDetPart = (jet1->Pt()-jet3->Pt())/jet3->Pt();
                
                fJetProbDensityDetPart->Fill(probDensDetPart,jet3->Pt());
                fJetProbDensityPartDet->Fill(probDensPartDet,jet1->Pt());
                
            }
            
            if (fJetShapeType == kGenOnTheFly){
                const AliEmcalPythiaInfo *partonsInfo = 0x0;
                partonsInfo = GetPythiaInfo();
                Double_t jp1=RelativePhi(jet1->Phi(),partonsInfo->GetPartonPhi6());
                Double_t detap1=(jet1->Eta())-(partonsInfo->GetPartonEta6());
                kWeight=partonsInfo->GetPythiaEventWeight();
                fh2ResponseW->Fill(jet1->Pt(),jet1->Pt(),kWeight);
                
                Float_t dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
                fEtaJetCorr6->Fill(jet1->Eta(), partonsInfo->GetPartonEta6());
                fPhiJetCorr6->Fill(jet1->Phi(), partonsInfo->GetPartonPhi6());
                if(dRp1 < fRMatching) {
                    fShapesVar[0] = partonsInfo->GetPartonFlag6();
                    fPtJetCorr ->Fill(partonsInfo->GetPartonPt6(), jet1->Pt());
                }
                else {
                    jp1=RelativePhi(jet1->Phi(),partonsInfo->GetPartonPhi7());
                    detap1=(jet1->Eta())-(partonsInfo->GetPartonEta7());
                    dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
                    fEtaJetCorr7->Fill(jet1->Eta(), partonsInfo->GetPartonEta7());
                    fPhiJetCorr7->Fill(jet1->Phi(), partonsInfo->GetPartonPhi7());
                    if(dRp1 < fRMatching) {
                        fShapesVar[0] = partonsInfo->GetPartonFlag7();
                        fPtJetCorr->Fill(partonsInfo->GetPartonPt7(), jet1->Pt());
                    }
                    else fShapesVar[0]=0;
                }
            }
            
            Double_t ptSubtracted = 0;
            if (fJetShapeSub==kConstSub) ptSubtracted= jet1->Pt();
            
            else if (fJetShapeSub==kDerivSub)  {
                ptSubtracted=jet1->Pt()-GetRhoVal(0)*jet1->Area();
            }
            
            else if (fJetShapeSub==kNoSub) {
                if ((fJetShapeType==kData) || (fJetShapeType==kDetEmbPartPythia)) ptSubtracted=jet1->Pt()-GetRhoVal(0)*jet1->Area();
                else if ((fJetShapeType==kPythiaDef) || (fJetShapeType==kMCTrue) || (fJetShapeType==kGenOnTheFly)) ptSubtracted= jet1->Pt();
            }
            
            //Printf("ptSubtracted=%f,fPtThreshold =%f ", ptSubtracted, fPtThreshold);
            if (ptSubtracted < fPtThreshold) continue;
            
            if (fOneConstSelectOn == kTRUE) fNbOfConstvspT->Fill(GetJetNumberOfConstituents(jet1,0), ptSubtracted);
            if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1)) continue;
            
            fShapesVar[1] = ptSubtracted;
            fShapesVar[2] = GetJetpTD(jet1,0);
            fShapesVar[3] = GetJetMass(jet1,0);
            fShapesVar[4] = GetJetAngularity(jet1,0);
            fShapesVar[5] = GetJetCircularity(jet1,0);
            fShapesVar[6] = GetJetLeSub(jet1,0);
            fShapesVar[7] = GetJetCoronna(jet1,0);
            
            Float_t ptMatch=0., ptDMatch=0., massMatch=0., angulMatch=0.,circMatch=0., lesubMatch=0., coronnaMatch=0;
            //Float constMatch=0., sigma2Match=0.;
            Int_t kMatched = 0;
            
            if (fJetShapeType==kPythiaDef) {
                kMatched =1;
                if(fJetShapeSub==kConstSub) kMatched = 3;
                
                if (!jet3) {
                    Printf("jet3 does not exist, returning");
                    continue;
                }
                
                ptMatch=jet3->Pt();
                ptDMatch=GetJetpTD(jet3, kMatched);
                massMatch=GetJetMass(jet3,kMatched);
                //constMatch=1.*GetJetNumberOfConstituents(jet2,kMatched);
                angulMatch=GetJetAngularity(jet3, kMatched);
                circMatch=GetJetCircularity(jet3, kMatched);
                lesubMatch=GetJetLeSub(jet3, kMatched);
                coronnaMatch=GetJetCoronna(jet3,kMatched);
                //sigma2Match = GetSigma2(jet2, kMatched);
            }
            
            if (fJetShapeType==kDetEmbPartPythia) {
                if(fJetShapeSub==kConstSub) kMatched = 3;
                if(fJetShapeSub==kDerivSub) kMatched = 2;
                ptMatch=jet3->Pt();
                ptDMatch=GetJetpTD(jet3, kMatched);
                massMatch=GetJetMass(jet3,kMatched);
                // constMatch=1.*GetJetNumberOfConstituents(jet3,kMatched);
                angulMatch=GetJetAngularity(jet3, kMatched);
                circMatch=GetJetCircularity(jet3, kMatched);
                lesubMatch=GetJetLeSub(jet3, kMatched);
                coronnaMatch = GetJetCoronna(jet3, kMatched);
            }
            
            
            
            if (fJetShapeType == kMCTrue || fJetShapeType == kData || fJetShapeType == kGenOnTheFly) {
                kMatched = 0;
                ptMatch=0.;
                ptDMatch=0.;
                massMatch=0.;
                //constMatch=0.;
                angulMatch=0.;
                circMatch=0.;
                lesubMatch=0.;
                coronnaMatch =0.;
                
            }
            
            fShapesVar[8] = ptMatch;
            fShapesVar[9] = ptDMatch;
            fShapesVar[10] = massMatch;
            fShapesVar[11] = angulMatch;
            fShapesVar[12] = circMatch;
            fShapesVar[13] = lesubMatch;
            fShapesVar[14] = coronnaMatch;
            fShapesVar[15] = kWeight;
            fShapesVar[16] = rhoVal;
            fShapesVar[17] = rhoMassVal;
            fShapesVar[18] = jet1->Pt();
            
            Int_t nInclusiveElectrons = 0, nPhotonicElectrons = 0, nTrueElectronsMC= 0, nTrueHFElecMC= 0;
            Double_t pElec = 0., ptElec = 0.;
            Bool_t hasElectrons = kFALSE;
            
            GetNumberOfElectrons(jet1, 0,nMotherKink,listofmotherkink,nInclusiveElectrons,nPhotonicElectrons,pElec,ptElec,hasElectrons);
            
            fAngChargPart->Fill(jet1->Pt(),GetJetAngularity(jet1,0));
            fDispChargPart->Fill(jet1->Pt(),GetJetpTD(jet1,0));
            
            if(nInclusiveElectrons==1){
                fAngIncElec->Fill(jet1->Pt(),GetJetAngularity(jet1,0));
                fDispIncElec->Fill(jet1->Pt(),GetJetpTD(jet1,0));
            }
            
            if(!hasElectrons){
                fAngHadron->Fill(jet1->Pt(),GetJetAngularity(jet1,0));
                fDispHadron->Fill(jet1->Pt(),GetJetpTD(jet1,0));
            }
            
            
            
            // generated HFE jets
            
            AliVParticle *vp1 = 0x0;
            Int_t pdgCode = 0;
            Int_t elecCounter = 0;
            Double_t ptTrueHFE = -1.;
            
            if(fMCarray){
                
                GetNumberOfTrueElectrons(jet1, 0,nMotherKink,listofmotherkink,nTrueElectronsMC,nTrueHFElecMC, ptTrueHFE);
                
                for(UInt_t i = 0; i < jet1->GetNumberOfTracks(); i++) {
                    vp1 = static_cast<AliVParticle*>(jet1->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
                    
                    if (!vp1){
                        Printf("AliVParticle associated to constituent not found");
                        continue;
                    }
                    
                    pdgCode = vp1->PdgCode();
                    if (TMath::Abs(pdgCode)==11) elecCounter++;
                }
                if (elecCounter==1) fPtGenJet->Fill(jet1->Pt());
            }
            
            fShapesVar[19] = ptTrueHFE;
            fShapesVar[20] = ptElec;
            fShapesVar[21] = nInclusiveElectrons;
            fShapesVar[22] = nPhotonicElectrons;
            fShapesVar[23] = hasElectrons;
            fShapesVar[24] = nTrueElectronsMC;
            fShapesVar[25] = nTrueHFElecMC;
            
            fTreeObservableTagging->Fill();
            
            
        }//jet loop
    }
    return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetMass(AliEmcalJet *jet,Int_t jetContNb=0){
    //calc subtracted jet mass
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtracted();
        else return jet->GetShapeProperties()->GetSecondOrderSubtracted();
        else
            return jet->M();
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHfeTagging::GetNumberOfElectrons(AliEmcalJet *jet, Int_t jetContNb , Int_t nMother, Double_t listMother[] ,  Int_t &nIncElec,  Int_t &nPhotElec, Double_t &pElec, Double_t &ptElec, Bool_t &hasElec){
    // count the number of inclusive electrons per jet and per event
    
    AliVParticle *vp1 = 0x0;
    Int_t nIE=0, nPairs=0, nPE=0, sub = -1;
    Float_t ratioElec = -1.;
    Double_t pte=0., pe=0.;
    Bool_t hasElecCand = kFALSE;
    
    Double_t jetPt = jet->Pt();
    
    Double_t ptRange[34] = {0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4,
        1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 5.,
        6., 8., 10., 12., 14., 16., 19., 22., 25., 30.,
        35., 40., 45., 50.};
    Double_t ptJetRange[6] = {5,20,40,60,80,120};
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (jet->GetNumberOfTracks()){
        
        fnPartPerJet->Fill(jet->GetNumberOfTracks());
        fpidResponse = fInputHandler->GetPIDResponse();
        if(!fpidResponse){
            printf("ERROR: pid response not available\n");
        }
        
        for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
            vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
            
            if (!vp1){
                Printf("AliVParticle associated to constituent not found");
                continue;
            }
            
            AliVTrack *vtrack = dynamic_cast<AliVTrack*>(vp1);
            if (!vtrack) {
                printf("ERROR: Could not receive track%d\n", i);
                continue;
            }
            
            AliAODTrack *track = dynamic_cast<AliAODTrack*>(vtrack);
            if (!track) continue;
            
            // track cuts for electron identification
            Bool_t passTrackCut = kFALSE;
            passTrackCut = InclElecTrackCuts(pVtx,track,nMother,listMother);
            if (!passTrackCut) continue;
            
            Double_t p=-9., pt=-9., fTPCnSigma=-99., fTOFnSigma=-99., phi = -9., eta = -99.;
            p = track->P();
            pt = track->Pt();
            phi = track->Phi();
            eta = track->Eta();
            
            fPtTrack->Fill(pt);
            fPhiTrack->Fill(pt,phi);
            fEtaTrack->Fill(pt,eta);
            if (pt > 0.5) fEtaPhiTrack->Fill(phi,eta);
            
            nPairs=0;
            
            fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            fTOFnSigma = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
            
            fnTPCnTOFnocut->Fill(fTPCnSigma,fTOFnSigma);
            fnTPCnocutP->Fill(p,fTPCnSigma);
            fnTOFnocutP->Fill(p,fTOFnSigma);
            if (TMath::Abs(fTOFnSigma)<fSigmaTOFcut){
                fnTPCcutP->Fill(p,fTPCnSigma);
                fnTPCcutPt->Fill(pt,fTPCnSigma);
                
                for (Int_t l=0;l<5;l++){// pt jet range
                    for(Int_t k=0;k<18;k++){// pt electron range (up to 4 GeV/c)
                        if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1] && p>=ptRange[k] && p<ptRange[k+1]) fnTPCSigma[l][k]->Fill(fTPCnSigma);
                    }
                }
            }
            
            
            if(TMath::Abs(fTPCnSigma)<fSigmTPCcutExcElec) hasElecCand = kTRUE;
            
            if (TMath::Abs(fTOFnSigma)<fSigmaTOFcut) fPtP->Fill(pt,p);
            
            if (fTPCnSigma>fSigmaTPCcutLowPt  && fTPCnSigma<3 && TMath::Abs(fTOFnSigma)<fSigmaTOFcut && pt>=fMinPtTPC && pt<fMaxPtTPC){
                fPhiRecElecTPC->Fill(pt,phi);
                fEtaRecElecTPC->Fill(pt,eta);
                fEtaPhiRecElecTPC->Fill(phi,eta);
                
                nIE++;
                pte=pt;
                pe=p;
                nPairs = GetNumberOfPairs(jet,track,pVtx,nMother,listMother,0,1);
                if (nPairs>0) nPE++;
            }
            
            // Electron ID with EMCal
            
            Double_t clsE = -9., m02 = -9., m20 = -9., clsTime = -9, EovP = -9.;
            
            Double_t emcphimim = 1.39;
            Double_t emcphimax = 3.265;
            
            Int_t clsId = track->GetEMCALcluster();
            AliVCluster* cluster=0x0;
            
            AliClusterContainer* clusterCont = GetClusterContainer(0);
            
            if (clsId>=0 && clusterCont){
                
                cluster = clusterCont->GetCluster(clsId);
                
                if(cluster && cluster->IsEMCAL()){
                    
                    Float_t  emcx[3]; // cluster pos
                    cluster->GetPosition(emcx);
                    TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
                    Double_t emcphi = clustpos.Phi();
                    
                    if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi());
                    
                    if (emcphi>emcphimim && emcphi<emcphimax){
                        clsE = cluster->GetNonLinCorrEnergy();
                        m20 = cluster->GetM20();
                        m02 = cluster->GetM02();
                        clsTime = cluster->GetTOF()*1e+9; // ns
                        
                        EovP = clsE/p;
                    }
                }
            }
            
            fnEovPelecNoTPCcut->Fill(pt,EovP);
            
            
            if (fTPCnSigma>fSigmaTPCcutHighPt  && fTPCnSigma<3) fnEovPelecTPCcut->Fill(pt,EovP);
            
            if (fTPCnSigma>fSigmaTPCcutHighPt  && fTPCnSigma<3 && m20 > 0.01 && m20 < fM20cut){
                fnEovPelecTPCEMCalcut->Fill(pt,EovP);
                
                for (Int_t l=0;l<5;l++){// pt jet range
                    if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1]) fnEovPelecTPCsscut[l]->Fill(pt,EovP);
                }
            }
            
            if (TMath::Abs(fTPCnSigma) > fSigmTPCcutExcElec) fnEovPbackg->Fill(pt,EovP);
            fnClsE->Fill(pt,clsE);
            fnM20->Fill(pt,m20);
            fnM02->Fill(pt,m02);
            fnClsTime->Fill(pt,clsTime);
            
            if (fTPCnSigma>fSigmaTPCcutHighPt  && fTPCnSigma<3 && EovP>fMinEoPcut && EovP<fMaxEoPcut && m20 > 0.01 && m20 < fM20cut && pt>=fMinPtEMCal && pt<fMaxPtEMCal){
                fPhiRecElecEMCal->Fill(pt,phi);
                fEtaRecElecEMCal->Fill(pt,eta);
                fEtaPhiRecElecEMCal->Fill(phi,eta);
                
                nIE++;
                pte=pt;
                pe=p;
                nPairs = GetNumberOfPairs(jet,track,pVtx,nMother,listMother,0,1);
                if (nPairs>0) nPE++;
            }
        }//tagged track
    }
    
    sub = nIE - nPE;
    if (jet->GetNumberOfTracks()>0) ratioElec = nIE/jet->GetNumberOfTracks();
    fnInclElecPerJet->Fill(nIE);
    fnPhotElecPerJet->Fill(nPE);
    fnIncSubPhotElecPerJet->Fill(sub);
    fnElecOverPartPerJet->Fill(ratioElec);
    
    nIncElec = nIE;
    nPhotElec = nPE;
    pElec = pe;//to be used for jets with only one IE
    ptElec = pte;//to be used for jets with only one IE
    hasElec = hasElecCand;
}
//________________________________________________________________________
void AliAnalysisTaskEmcalHfeTagging::GetNumberOfTrueElectrons(AliEmcalJet *jet, Int_t jetContNb ,  Int_t nMother, Double_t listMother[] ,  Int_t &nTrueElec,  Int_t &nTrueHFElec, Double_t &ptTrueHFElec){
    // count the number of inclusive and HF electrons per jet and per event (true MC)
    
    AliVParticle *vp1 = 0x0;
    Int_t nIE=0, nRecIE=0, nHFE=0, nPE=0, nPairs=0, iDecay = 0, nDmeson = 0, nBmeson = 0, nElecFromB = 0, nElecFromD = 0, nElecFromDfromB = 0, nQuark = 0, nGluon = 0, nBeauty = 0, nCharm = 0;
    Double_t p=-9., pt=-9., fTPCnSigma=-99., fTOFnSigma=-99., MCweight = 1., eta = -99., phi = -99., pte=0.;
    
    Double_t ptJetRange[6] = {5.,20.,40.,60.,80.,120.};
    Double_t angRange[6] = {0.,0.02,0.04,0.06,0.08,0.12};
    Double_t dispRange[7] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
    
    Bool_t isFromHFdecay=kFALSE;
    Bool_t isFromLMdecay=kFALSE;
    Bool_t passTrackCut = kFALSE;
    
    Double_t jetPt = jet->Pt();
    Double_t ang = GetJetAngularity(jet,0);
    Double_t disp = GetJetpTD(jet,0);
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (jet->GetNumberOfTracks()){
        
        fpidResponse = fInputHandler->GetPIDResponse();
        if(!fpidResponse){
            printf("ERROR: pid response not available\n");
        }
        
        for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
            vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
            
            if (!vp1){
                Printf("AliVParticle associated to constituent not found");
                continue;
            }
            
            AliVTrack *vtrack = dynamic_cast<AliVTrack*>(vp1);
            if (!vtrack) {
                printf("ERROR: Could not receive track%d\n", i);
                continue;
            }
            
            AliAODTrack *track = dynamic_cast<AliAODTrack*>(vtrack);
            if (!track) continue;
            
            passTrackCut = kFALSE;
            isFromHFdecay=kFALSE;
            isFromLMdecay=kFALSE;
            
            p = track->P();
            pt = track->Pt();
            eta = track->Eta();
            phi = track->Phi();
            
            fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            fTOFnSigma = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
            
            nPairs=0;
            MCweight = 1.;
            iDecay = 0;
            
            // Electron ID with EMCal
            
            Double_t EovP = -9., clsE = -9., m20 = -9.;
            
            Double_t emcphimim = 1.39;
            Double_t emcphimax = 3.265;
            
            Int_t clsId = track->GetEMCALcluster();
            
            AliVCluster* cluster=0x0;
            
            AliClusterContainer* clusterCont = GetClusterContainer(0);
            
            if (clsId>=0 && clusterCont){

                cluster = clusterCont->GetCluster(clsId);
                
                if(cluster && cluster->IsEMCAL() && phi > emcphimim && phi < emcphimax){
                    clsE = cluster->GetNonLinCorrEnergy();
                    m20 = cluster->GetM20();
                }
            }
            
            EovP = clsE/p;
            
            if(fMCarray){
                Int_t label = track->GetLabel();
                if(label!=0){
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                    if(fMCparticle){
                        Int_t partPDG = TMath::Abs(fMCparticle->GetPdgCode());
                        
                        if (((partPDG/100)%10) == 4) nDmeson++;
                        if (((partPDG/100)%10) == 5) nBmeson++;
                        
                        
                        Int_t idMother = fMCparticle->GetMother();
                        
                        if (idMother>0){
                            AliAODMCParticle* mother = (AliAODMCParticle*) fMCarray->At(idMother);
                            Int_t motherPDG = TMath::Abs(mother->GetPdgCode());
                            
                            if ((partPDG==11) && ((motherPDG/100)%10) == 5) nElecFromB++;
                            
                            if (((motherPDG/100)%10) == 4) nDmeson++;
                            if (((motherPDG/100)%10) == 5) nBmeson++;
                            
                            
                            Int_t idSecondMother = mother->GetMother();
                            if (idSecondMother>0){
                                AliAODMCParticle* secondMother = (AliAODMCParticle*) fMCarray->At(idSecondMother);
                                Int_t secondMotherPDG = TMath::Abs(secondMother->GetPdgCode());
                                
                                if ((partPDG==11) &&  (((secondMotherPDG/100)%10) != 5) && (((motherPDG/100)%10) == 4))
                                    nElecFromD++;
                                if ((partPDG==11) &&  (((secondMotherPDG/100)%10) == 5) && (((motherPDG/100)%10) == 4))
                                    nElecFromDfromB++;
                                
                                if (((secondMotherPDG/100)%10) == 4) nDmeson++;
                                if (((secondMotherPDG/100)%10) == 5) nBmeson++;
                                
                                
                                Int_t idThirdMother = secondMother->GetMother();
                                if (idThirdMother>0){
                                    AliAODMCParticle* thirdMother = (AliAODMCParticle*) fMCarray->At(idThirdMother);
                                    Int_t thirdMotherPDG = TMath::Abs(thirdMother->GetPdgCode());
                                    
                                    if ((partPDG==11) &&  (((motherPDG/100)%10) == 5) && secondMotherPDG == 5){
                                        if (thirdMotherPDG < 9) nQuark++;
                                        if (thirdMotherPDG == 21) nGluon++;
                                    }
                                    
                                    if ((partPDG==11) &&  (((motherPDG/100)%10) == 4) && secondMotherPDG == 4){
                                        if (thirdMotherPDG < 9) nQuark++;
                                        if (thirdMotherPDG == 21) nGluon++;
                                    }
                                    
                                    if (((thirdMotherPDG/100)%10) == 4) nDmeson++;
                                    if (((thirdMotherPDG/100)%10) == 5) nBmeson++;
                                    
                                    
                                    Int_t idForthMother = thirdMother->GetMother();
                                    if (idForthMother>0){
                                        AliAODMCParticle* forthMother = (AliAODMCParticle*) fMCarray->At(idForthMother);
                                        Int_t forthMotherPDG = TMath::Abs(forthMother->GetPdgCode());
                                        
                                        if ((partPDG==11) &&  (((secondMotherPDG/100)%10) == 5) && (((motherPDG/100)%10) == 4) && thirdMotherPDG == 5) {
                                            if (forthMotherPDG < 9) nQuark++;
                                            if (forthMotherPDG == 21) nGluon++;
                                        }
                                        
                                        if (((forthMotherPDG/100)%10) == 4) nDmeson++;
                                        if (((forthMotherPDG/100)%10) == 5) nBmeson++;
                                        
                                        if (thirdMotherPDG == 4 || forthMotherPDG == 4) nCharm++;
                                        if (thirdMotherPDG == 5 || forthMotherPDG == 5) nBeauty++;
                                        
                                    }
                                }//3rd mother
                            }//2nd mother
                        }//1st mother
                        
                        GetWeightAndDecay(fMCparticle,iDecay,MCweight);
                        isFromHFdecay = IsFromHFdecay(fMCparticle);
                        isFromLMdecay = IsFromLMdecay(fMCparticle);
                        
                        
                        if (partPDG == 11){
                            nIE++;
                            pte=pt;
                            
                            if (isFromHFdecay) nHFE++;
                            if (iDecay>0 && iDecay<7) nPE++;
                        }
                        
                        if ((partPDG==11) && isFromHFdecay && nHFE<2){
                            fptTrueHFEeffTPCTOF[0]->Fill(pt);
                            fptTrueHFEeffEMCal[0]->Fill(pt);
                            
                            fptTrueHFEeffTPCTOFang[0]->Fill(pt,jet->Pt(),GetJetAngularity(jet,0));
                            fptTrueHFEeffEMCalang[0]->Fill(pt,jet->Pt(),GetJetAngularity(jet,0));
                            
                            fptTrueHFEeffTPCTOFdisp[0]->Fill(pt,jet->Pt(),GetJetpTD(jet,0));
                            fptTrueHFEeffEMCaldisp[0]->Fill(pt,jet->Pt(),GetJetpTD(jet,0));
                        }
                        
                        // track cuts
                        passTrackCut = InclElecTrackCuts(pVtx,track,nMother,listMother);
                        if (!passTrackCut) continue;
                        
                        if ((partPDG==11) && isFromHFdecay  && nHFE<2){
                            fptTrueHFEeffTPCTOF[1]->Fill(pt);
                            fptTrueHFEeffEMCal[1]->Fill(pt);
                            
                            if (fTPCnSigma>fSigmaTPCcutLowPt && fTPCnSigma<3)
                                fptTrueHFEeffTPCTOF[2]->Fill(pt);
                            
                            if (TMath::Abs(fTOFnSigma)<fSigmaTOFcut)
                                fptTrueHFEeffTPCTOF[3]->Fill(pt);
                            
                            if (fTPCnSigma>fSigmaTPCcutLowPt && fTPCnSigma<3 && TMath::Abs(fTOFnSigma)<fSigmaTOFcut){
                                fptTrueHFEeffTPCTOF[4]->Fill(pt);
                                fptTrueHFEeffTPCTOFang[1]->Fill(pt,jet->Pt(),GetJetAngularity(jet,0));
                                fptTrueHFEeffTPCTOFdisp[1]->Fill(pt,jet->Pt(),GetJetpTD(jet,0));
                            }
                            
                            
                            if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && m20 > 0.01 && m20 < fM20cut)
                                fptTrueHFEeffEMCal[2]->Fill(pt);
                            
                            if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && EovP>fMinEoPcut && EovP<fMaxEoPcut)
                                fptTrueHFEeffEMCal[3]->Fill(pt);
                            
                            if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && EovP>fMinEoPcut && EovP<fMaxEoPcut && m20 > 0.01 && m20 < fM20cut){
                                
                                fptTrueHFEeffEMCal[4]->Fill(pt);
                                fptTrueHFEeffEMCalang[1]->Fill(pt,jet->Pt(),GetJetAngularity(jet,0));
                                fptTrueHFEeffEMCaldisp[1]->Fill(pt,jet->Pt(),GetJetpTD(jet,0));
                            }
                        }
                        
                        if (fTPCnSigma>fSigmaTPCcutLowPt && fTPCnSigma<3 && TMath::Abs(fTOFnSigma)<fSigmaTOFcut && pt>=fMinPtTPC && pt<fMaxPtTPC) nRecIE++;
                        if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && EovP>fMinEoPcut && EovP<fMaxEoPcut &&
                            m20 > 0.01 && m20 < fM20cut && pt>=fMinPtEMCal && pt<fMaxPtEMCal) nRecIE++;
                        
                        // TPC-TOF
                        if (fTPCnSigma>fSigmaTPCcutLowPt && fTPCnSigma<3 && TMath::Abs(fTOFnSigma)<fSigmaTOFcut && pt>=fMinPtTPC && pt<fMaxPtTPC & nRecIE <2){
                            nPairs = GetNumberOfPairs(jet,track,pVtx,nMother,listMother,iDecay,MCweight);
                            if (nPairs>0 && iDecay>0 && iDecay<7) fptRecPE->Fill(pt,MCweight);
                            if (iDecay>0 && iDecay<7) fptTruePE->Fill(pt,MCweight);
                            
                            for (Int_t l=0;l<5;l++){// pt jet range
                                if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1]){
                                    
                                    for (Int_t m=0;m<5;m++){// angularity
                                        if (ang>=angRange[m] && ang<angRange[m+1]){
                                            if (iDecay>0 && iDecay<7) fTotPEAng[l][m]->Fill(pt,MCweight);
                                            if (nPairs>0 && iDecay>0 && iDecay<7) fRecPEAng[l][m]->Fill(pt,MCweight);
                                        }
                                    }
                                    
                                    for (Int_t m=0;m<6;m++){// dispersion
                                        if (disp>=dispRange[m] && disp<dispRange[m+1]){
                                            if (iDecay>0 && iDecay<7) fTotPEDisp[l][m]->Fill(pt,MCweight);
                                            if (nPairs>0 && iDecay>0 && iDecay<7) fRecPEDisp[l][m]->Fill(pt,MCweight);
                                        }
                                    }
                                }
                            }
                        }// PID cuts
                        
                        //EMCal
                        if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && EovP>fMinEoPcut && EovP<fMaxEoPcut &&
                            m20 > 0.01 && m20 < fM20cut && pt>=fMinPtEMCal && pt<fMaxPtEMCal  & nRecIE <2){
                            nPairs = GetNumberOfPairs(jet,track,pVtx,nMother,listMother,iDecay,MCweight);
                            if (nPairs>0 && iDecay>0 && iDecay<7) fptRecPE->Fill(pt,MCweight);
                            if (iDecay>0 && iDecay<7) fptTruePE->Fill(pt,MCweight);
                            
                            for (Int_t l=0;l<5;l++){// pt jet range
                                if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1]){
                                    
                                    for (Int_t m=0;m<5;m++){// angularity
                                        if (ang>=angRange[m] && ang<angRange[m+1]){
                                            if (iDecay>0 && iDecay<7) fTotPEAng[l][m]->Fill(pt,MCweight);
                                            if (nPairs>0 && iDecay>0 && iDecay<7) fRecPEAng[l][m]->Fill(pt,MCweight);
                                        }
                                    }
                                    
                                    for (Int_t m=0;m<6;m++){// dispersion
                                        if (disp>=dispRange[m] && disp<dispRange[m+1]){
                                            if (iDecay>0 && iDecay<7) fTotPEDisp[l][m]->Fill(pt,MCweight);
                                            if (nPairs>0 && iDecay>0 && iDecay<7) fRecPEDisp[l][m]->Fill(pt,MCweight);
                                        }
                                    }
                                }
                            }
                        }// PID cuts
                    }
                }
            }
        }//track loop
    }
    
    if (nIE==1) fptJetIE->Fill(jet->Pt());
    if (nIE==1 && nPE==1) fptJetPE->Fill(jet->Pt());
    if (nIE==1 && nHFE==1) fptJetHFE->Fill(jet->Pt());
    
    fnTrueElecPerJet->Fill(nIE);
    fnTrueHFElecPerJet->Fill(nHFE);
    fnTruePElecPerJet->Fill(nPE);
    
    nTrueElec = nIE;
    nTrueHFElec = nHFE;
    
    if (nIE==1 && nHFE==1) ptTrueHFElec = pte;
    
    if (nCharm > 1){
        fAngCharm->Fill(jet->Pt(),GetJetAngularity(jet,0));
        fDispCharm->Fill(jet->Pt(),GetJetpTD(jet,0));
    }
    
    if (nBeauty > 1){
        fAngBeauty->Fill(jet->Pt(),GetJetAngularity(jet,0));
        fDispBeauty->Fill(jet->Pt(),GetJetpTD(jet,0));
    }
    
    if (nDmeson > 1){
        fAngD->Fill(jet->Pt(),GetJetAngularity(jet,0));
        fDispD->Fill(jet->Pt(),GetJetpTD(jet,0));
    }
    
    if (nBmeson > 1){
        fAngB->Fill(jet->Pt(),GetJetAngularity(jet,0));
        fDispB->Fill(jet->Pt(),GetJetpTD(jet,0));
    }
    
    if (nElecFromB == 1){
        fAngElecFromB->Fill(jet->Pt(),GetJetAngularity(jet,0));
        fDispElecFromB->Fill(jet->Pt(),GetJetpTD(jet,0));
    }
    
    if (nElecFromD == 1){
        fAngElecFromD->Fill(jet->Pt(),GetJetAngularity(jet,0));
        fDispElecFromD->Fill(jet->Pt(),GetJetpTD(jet,0));
    }
    
    if (nElecFromDfromB == 1){
        fAngElecFromDFromB->Fill(jet->Pt(),GetJetAngularity(jet,0));
        fDispElecFromDFromB->Fill(jet->Pt(),GetJetpTD(jet,0));
    }
    
    if (nQuark == 1){
        fAngQuark->Fill(jet->Pt(),GetJetAngularity(jet,0));
        fDispQuark->Fill(jet->Pt(),GetJetpTD(jet,0));
    }
    if (nGluon == 1){
        fAngGluon->Fill(jet->Pt(),GetJetAngularity(jet,0));
        fDispGluon->Fill(jet->Pt(),GetJetpTD(jet,0));
    }
}

//_________________________________________
Int_t AliAnalysisTaskEmcalHfeTagging::GetNumberOfPairs(AliEmcalJet *jet, AliAODTrack *track,const AliVVertex *pVtx,
                                                       Int_t nMother, Double_t listMother[],Int_t decay, Double_t weight)
{
    //Get number of ULS and LS pairs per event
    
    Int_t ntracks = 0;
    ntracks = fVevent->GetNumberOfTracks();
    
    Int_t nULSpairs = 0;
    Int_t nLSpairs = 0;
    Int_t sub = -1;
    
    Double_t ptJetRange[6] = {5,20,40,60,80,120};
    Double_t angRange[6] = {0.,0.02,0.04,0.06,0.08,0.12};
    Double_t dispRange[7] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
    
    Double_t jetPt = jet->Pt();
    Double_t  ang = GetJetAngularity(jet,0);
    Double_t disp = GetJetpTD(jet,0);
    
    for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
        
        AliVParticle* Vassotrack = fVevent->GetTrack(jTracks);
        
        if (!Vassotrack) {
            printf("ERROR: Could not receive associated track %d\n", jTracks);
            continue;
        }
        AliAODTrack *trackAsso = dynamic_cast<AliAODTrack*>(Vassotrack);
        if(!trackAsso) continue;
        
        if((track->GetID())==(trackAsso->GetID())) continue;
        
        // track cuts
        Bool_t passAssoTrackCut = kFALSE;
        passAssoTrackCut = PhotElecTrackCuts(pVtx,trackAsso,nMother,listMother);
        if (!passAssoTrackCut) continue;
        
        // tagged particle
        Double_t p=-9.,pt=-9.,eta =-9.,phi=-9.;
        Int_t charge = 0;
        
        pt = track->Pt();
        p = track->P();
        phi = track->Phi();
        eta = track->Eta();
        charge = track->Charge();
        
        
        // associated particle variables
        Double_t pAsso=-9.,ptAsso=-9.,etaAsso =-9.,fITSnSigmaAsso=-9.,fTPCnSigmaAsso=-9.,phiAsso=-9.;
        Int_t chargeAsso = 0;
        
        ptAsso = trackAsso->Pt();
        pAsso = trackAsso->P();
        phiAsso = trackAsso->Phi();
        etaAsso = trackAsso->Eta();
        chargeAsso = trackAsso->Charge();
        
        
        // looser PID cuts
        fITSnSigmaAsso = fpidResponse->NumberOfSigmasITS(trackAsso, AliPID::kElectron);
        fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(trackAsso, AliPID::kElectron);
        
        if(TMath::Abs(fTPCnSigmaAsso)>3.) continue;
        
        // invariant mass
        Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
        Double_t openingAngle = -999., mass=999., width = -999;
        
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        AliKFParticle::SetField(fVevent->GetMagneticField());
        AliKFParticle ge1(*track, fPDGe1);
        AliKFParticle ge2(*trackAsso, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        openingAngle = ge1.GetAngle(ge2);
        //if(openingAngle > fOpeningAngleCut) continue;
        
        Int_t MassCorrect=-9;
        MassCorrect = recg.GetMass(mass,width);
        
        for (Int_t l=0;l<5;l++){// pt jet range
            if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1] && fFlagULS) fInvmassULS[l]->Fill(pt,mass);
            if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1] && fFlagLS)  fInvmassLS[l]->Fill(pt,mass);
        }
        
        if(mass<fIMcut && fFlagULS){
            nULSpairs++;
            fAngULS->Fill(jetPt,ang);
            fDispULS->Fill(jetPt,disp);
        }
        
        if(mass<fIMcut && fFlagLS){
            nLSpairs++;
            fAngLS->Fill(jetPt,ang);
            fDispLS->Fill(jetPt,disp);
        }
        
        
        for (Int_t l=0;l<5;l++){// pt jet range
            if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1]){
                
                for (Int_t m=0;m<5;m++){// angularity
                    if (ang>=angRange[m] && ang<angRange[m+1] && decay>0 && decay<7 && mass<fIMcut){
                        if (fFlagULS) fULSptAng[l][m]->Fill(pt,weight);
                        if (fFlagLS) fLSptAng[l][m]->Fill(pt,weight);
                    }
                }
                
                for (Int_t m=0;m<6;m++){// dispersion
                    if (disp>=dispRange[m] && disp<dispRange[m+1]  && decay>0 && decay<7 && mass<fIMcut){
                        if (fFlagULS) fULSptDisp[l][m]->Fill(pt,weight);
                        if (fFlagLS) fLSptDisp[l][m]->Fill(pt,weight);
                    }
                }
            }
        }
        
        
        
    }//track loop
    
    sub = nULSpairs-nLSpairs;
    fnULSmLSpairsPerElectron->Fill(track->Pt(),sub);
    
    if (sub>0){
        fAngPhotElec->Fill(jet->Pt(),GetJetAngularity(jet,0));
        fDispPhotElec->Fill(jet->Pt(),GetJetpTD(jet,0));
    }
    
    return sub;
}
//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::IsFromHFdecay(AliAODMCParticle *particle)
{
    // Check if the mother comes from heavy-flavour decays
    Bool_t isHFdecay = kFALSE;
    
    Int_t idMother = particle->GetMother();
    if (idMother>0){
        AliAODMCParticle* mother = (AliAODMCParticle*) fMCarray->At(idMother);
        Int_t motherPDG = TMath::Abs(mother->GetPdgCode());
        
        // c decays
        if( (((motherPDG/1000)%10) == 4) || (((motherPDG/100)%10) == 4) ) isHFdecay = kTRUE;
        
        // b decays
        if( (((motherPDG/1000)%10) == 5) || (((motherPDG/100)%10) == 5) ) isHFdecay = kTRUE;
    }
    
    return isHFdecay;
}
//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::IsFromLMdecay(AliAODMCParticle *particle)
{
    // Check if  mother comes from light-meson decays
    Bool_t isLMdecay = kFALSE;
    
    Int_t idMother = particle->GetMother();
    if (idMother>0){
        AliAODMCParticle* mother = (AliAODMCParticle*) fMCarray->At(idMother);
        Int_t motherPDG = TMath::Abs(mother->GetPdgCode());
        
        if(motherPDG == 111 || motherPDG == 221 || motherPDG==223 || motherPDG==333 || motherPDG==331 ||
           motherPDG==113 || motherPDG==213 || motherPDG==313 || motherPDG==323) isLMdecay = kTRUE;
    }
    
    return isLMdecay;
}
//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::IsPrimary(AliAODMCParticle *particle)
{
    // Check if  particle is primary
    Bool_t isprimary = kFALSE;
    if (particle->IsPrimary()) isprimary = kTRUE;
    
    return isprimary;
}
//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::HasMother(AliAODMCParticle *particle)
{
    // Check if the particle has mother
    Bool_t hasMother = kTRUE;
    
    Int_t idMother = particle->GetMother();
    if (idMother < 0) hasMother = kFALSE;
    
    return hasMother;
}
//_________________________________________
Double_t AliAnalysisTaskEmcalHfeTagging::GetPi0weight(Double_t mcPi0pT) const
{
    //Get Pi0 weight
    double weight = 1.;
    double parPi0_enh[5] = {1.12,0.9988,0.227,7.493,-0.006416};
    
//    if (fRunNumber >= 252235 && fRunNumber <= 264347){//LHC17h8b (o/p)
//        parPi0_enh[0] = 11.4919;
//        parPi0_enh[1] = 0.021972;
//        parPi0_enh[2] = 0.133756;
//        parPi0_enh[3] = 1.67114;
//        parPi0_enh[4] = 5.80327;
//    }
//
//    if (fRunNumber >= 256504 && fRunNumber <= 259888){//LHC18f4b (k/l)
//        parPi0_enh[0] = 17.2349;
//        parPi0_enh[1] = 0.101117;
//        parPi0_enh[2] = 0.165525;
//        parPi0_enh[3] = 1.52289;
//        parPi0_enh[4] = 5.72523;
//    }
    
    if (fMCweight==1) weight = parPi0_enh[0] / TMath::Power((exp(parPi0_enh[1]*mcPi0pT - parPi0_enh[2]*mcPi0pT*mcPi0pT) + (mcPi0pT/parPi0_enh[3])) , parPi0_enh[4]);
    
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskEmcalHfeTagging::GetEtaweight(Double_t mcEtapT) const
{
    //Get Pi0 weight
    double weight = 1.;
    double parEta_enh[5] = {1.042,-1.171,0.1734,1.594,0.06675};
    
//    if (fRunNumber >= 252235 && fRunNumber <= 264347){//LHC17h8b (o/p)
//        parEta_enh[0] = 3.89483;
//        parEta_enh[1] = 2.8197;
//        parEta_enh[2] = 0.0662309;
//        parEta_enh[3] = 0.00563749;
//        parEta_enh[4] = 0.435908;
//    }
//
//    if (fRunNumber >= 256504 && fRunNumber <= 259888){//LHC18f4b (k/l)
//        parEta_enh[0] = 0.328435;
//        parEta_enh[1] = -0.450847;
//        parEta_enh[2] = 0.00972497;
//        parEta_enh[3] = 3.13337;
//        parEta_enh[4] = 5.91296;
//    }
    
    if (fMCweight==1) weight = parEta_enh[0] / TMath::Power((exp(parEta_enh[1]*mcEtapT - parEta_enh[2]*mcEtapT*mcEtapT) + (mcEtapT/parEta_enh[3])), parEta_enh[4]);
    
    return weight;
}
//________________________________________________________________________
void AliAnalysisTaskEmcalHfeTagging::GetWeightAndDecay(AliAODMCParticle *particle, Int_t &decay, Double_t &weight)
{
    //Get decay channel and weight for pi0 and eta (MC with enhanced signal)
    Double_t w = 1.;
    Int_t d = 0;
    Int_t partPDG = particle->GetPdgCode();
    
    if (TMath::Abs(partPDG)==11){
        Int_t idMother = particle->GetMother();
        
        if (idMother>0){
            AliAODMCParticle* mother = (AliAODMCParticle*) fMCarray->At(idMother);
            Int_t motherPDG = mother->GetPdgCode();
            Double_t motherPt = mother->Pt();
            
            Bool_t isMotherPrimary = IsPrimary(mother);
            Bool_t isMotherFromHF = IsFromHFdecay(mother);
            Bool_t isMotherFromLM = IsFromLMdecay(mother);
            Bool_t hasSecondMother = HasMother(mother);
            
            if (motherPDG==111 && isMotherPrimary){// pi0 -> e
                if (!hasSecondMother || (!isMotherFromHF && !isMotherFromLM)){
                    d = 1;
                    w = GetPi0weight(motherPt);
                }
            }
            
            if (motherPDG==221 && isMotherPrimary){// eta -> e
                if (!hasSecondMother || (!isMotherFromHF && !isMotherFromLM)){
                    d = 2;
                    w = GetEtaweight(motherPt);
                }
            }
            
            Int_t idSecondMother = mother->GetMother();
            
            if (idSecondMother>0){
                AliAODMCParticle* secondMother = (AliAODMCParticle*) fMCarray->At(idSecondMother);
                Int_t secondMotherPDG = secondMother->GetPdgCode();
                Double_t secondMotherPt = secondMother->Pt();
                
                Bool_t isSecondMotherPrimary = IsPrimary(secondMother);
                Bool_t isSecondMotherFromHF = IsFromHFdecay(secondMother);
                Bool_t isSecondMotherFromLM = IsFromLMdecay(secondMother);
                Bool_t hasThirdMother = HasMother(secondMother);
                
                if (motherPDG==22 && secondMotherPDG==111 && isSecondMotherPrimary){ //pi0 -> g -> e
                    if (!hasThirdMother || (!isSecondMotherFromHF && !isSecondMotherFromLM)){
                        d = 3;
                        w = GetPi0weight(secondMotherPt);
                    }
                }
                
                if (motherPDG==22 && secondMotherPDG==221  && isSecondMotherPrimary){ //eta -> g -> e
                    if (!hasThirdMother || (!isSecondMotherFromHF && !isSecondMotherFromLM)){
                        d = 4;
                        w = GetEtaweight(secondMotherPt);
                    }
                }
                
                if (motherPDG==111 && secondMotherPDG==221  && isSecondMotherPrimary){ //eta -> pi0 -> e
                    if (!hasThirdMother || (!isSecondMotherFromHF && !isSecondMotherFromLM)){
                        d = 5;
                        w = GetEtaweight(secondMotherPt);
                    }
                }
                
                Int_t idThirdMother = secondMother->GetMother();
                
                if (idThirdMother>0){
                    AliAODMCParticle* thirdMother = (AliAODMCParticle*) fMCarray->At(idThirdMother);
                    Int_t thirdMotherPDG = thirdMother->GetPdgCode();
                    Double_t thirdMotherPt = thirdMother->Pt();
                    
                    Bool_t isThirdMotherPrimary = IsPrimary(thirdMother);
                    Bool_t isThirdMotherFromHF = IsFromHFdecay(thirdMother);
                    Bool_t isThirdMotherFromLM = IsFromLMdecay(thirdMother);
                    Bool_t hasFourthMother = HasMother(thirdMother);
                    
                    if (motherPDG==22 && secondMotherPDG==111 && thirdMotherPDG==221 && isThirdMotherPrimary){//eta->pi0->g-> e
                        if (!hasFourthMother || (!isThirdMotherFromHF && !isThirdMotherFromLM)){
                            d = 6;
                            w = GetEtaweight(thirdMotherPt);
                        }
                    }
                }//third mother
            }//second mother
        }//mother
    }// if electron
    decay = d;
    weight = w;
}
//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::Angularity(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        Double_t dphi = RelativePhi(vp1->Phi(),jet->Phi());
        Double_t dr2 = (vp1->Eta()-jet->Eta())*(vp1->Eta()-jet->Eta()) + dphi*dphi;
        Double_t dr = TMath::Sqrt(dr2);
        num=num+vp1->Pt()*dr;
        den=den+vp1->Pt();
    }
    return num/den;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    if((fJetShapeSub==kDerivSub) && (jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity();
        else
            return Angularity(jet, jetContNb);
    
}

//____________________________________________________________________________

Float_t AliAnalysisTaskEmcalHfeTagging::Coronna(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    AliTrackContainer *PartCont = NULL;
    AliParticleContainer *PartContMC = NULL;
    
    
    if (fJetShapeSub==kConstSub){
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) PartContMC = GetParticleContainer(1);
        else PartCont = GetTrackContainer(1);
    }
    else{
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) PartContMC = GetParticleContainer(0);
        else PartCont = GetTrackContainer(0);
    }
    
    TClonesArray *TracksArray = NULL;
    TClonesArray *TracksArrayMC = NULL;
    
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) TracksArrayMC = PartContMC->GetArray();
    else TracksArray = PartCont->GetArray();
    
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly){
        if(!PartContMC || !TracksArrayMC) return -2;
    }
    else {
        if(!PartCont || !TracksArray) return -2;
    }
    
    
    AliAODTrack *Track = 0x0;
    Float_t sumpt=0;
    Int_t NTracks=0;
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) NTracks = TracksArrayMC->GetEntriesFast();
    else NTracks = TracksArray->GetEntriesFast();
    
    for(Int_t i=0; i < NTracks; i++){
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly){
            if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
                if (!Track) continue;
                if(TMath::Abs(Track->Eta())>0.9) continue;
                Double_t dphi = RelativePhi(Track->Phi(),jet->Phi());
                Double_t dr2 = (Track->Eta()-jet->Eta())*(Track->Eta()-jet->Eta()) + dphi*dphi;
                Double_t dr = TMath::Sqrt(dr2);
                if((dr>=0.8) && (dr<1)) sumpt=sumpt+Track->Pt();
            }
        }
        else{
            if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
                if (!Track) continue;
                if(TMath::Abs(Track->Eta())>0.9) continue;
                if (Track->Pt()<0.15) continue;
                Double_t dphi = RelativePhi(Track->Phi(),jet->Phi());
                Double_t dr2 = (Track->Eta()-jet->Eta())*(Track->Eta()-jet->Eta()) + dphi*dphi;
                Double_t dr = TMath::Sqrt(dr2);
                if((dr>=0.8) && (dr<1)) sumpt=sumpt+Track->Pt();
                
            }
        }
    }
    return sumpt;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetCoronna(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    if((fJetShapeSub==kDerivSub) && (jetContNb==0)) return -2;
    else
        return Coronna(jet, jetContNb);
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::PTD(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        num=num+vp1->Pt()*vp1->Pt();
        den=den+vp1->Pt();
    }
    return TMath::Sqrt(num)/den;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetpTD(AliEmcalJet *jet, Int_t jetContNb = 0){
    //calc subtracted jet mass
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedpTD();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedpTD();
        else
            return PTD(jet, jetContNb);
    
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::Circularity(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t mxx    = 0.;
    Double_t myy    = 0.;
    Double_t mxy    = 0.;
    int  nc     = 0;
    Double_t sump2  = 0.;
    Double_t pxjet=jet->Px();
    Double_t pyjet=jet->Py();
    Double_t pzjet=jet->Pz();
    
    
    //2 general normalized vectors perpendicular to the jet
    TVector3  ppJ1(pxjet, pyjet, pzjet);
    TVector3  ppJ3(- pxjet* pzjet, - pyjet * pzjet, pxjet * pxjet + pyjet * pyjet);
    ppJ3.SetMag(1.);
    TVector3  ppJ2(-pyjet, pxjet, 0);
    ppJ2.SetMag(1.);
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        TVector3 pp(vp1->Px(), vp1->Py(), vp1->Pz());
        
        //local frame
        TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
        TVector3 pPerp = pp - pLong;
        //projection onto the two perpendicular vectors defined above
        
        Float_t ppjX = pPerp.Dot(ppJ2);
        Float_t ppjY = pPerp.Dot(ppJ3);
        Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);
        if(ppjT<=0) return 0;
        
        mxx += (ppjX * ppjX / ppjT);
        myy += (ppjY * ppjY / ppjT);
        mxy += (ppjX * ppjY / ppjT);
        nc++;
        sump2 += ppjT;}
    
    if(nc<2) return 0;
    if(sump2==0) return 0;
    // Sphericity Matrix
    Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};
    TMatrixDSym m0(2,ele);
    
    // Find eigenvectors
    TMatrixDSymEigen m(m0);
    TVectorD eval(2);
    TMatrixD evecm = m.GetEigenVectors();
    eval  = m.GetEigenValues();
    // Largest eigenvector
    int jev = 0;
    //  cout<<eval[0]<<" "<<eval[1]<<endl;
    if (eval[0] < eval[1]) jev = 1;
    TVectorD evec0(2);
    // Principle axis
    evec0 = TMatrixDColumn(evecm, jev);
    Double_t compx=evec0[0];
    Double_t compy=evec0[1];
    TVector2 evec(compx, compy);
    Double_t circ=0;
    if(jev==1) circ=2*eval[0];
    if(jev==0) circ=2*eval[1];
    
    return circ;
    
    
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb =0 ){
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedCircularity();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedCircularity();
        else
            return Circularity(jet, jetContNb);
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::LeSub(AliEmcalJet *jet, Int_t jetContNb =0 ){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    AliVParticle *vp2 = 0x0;
    std::vector<int> ordindex;
    ordindex=jet->GetPtSortedTrackConstituentIndexes(jetCont->GetParticleContainer()->GetArray());
    //Printf("Nbof const = %d", jet->GetNumberOfTracks());
    //Printf("ordindex[0] = %d, ordindex[1] = %d", ordindex[0], ordindex[1]);
    
    if(ordindex.size()<2) return -1;
    
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[0], jetCont->GetParticleContainer()->GetArray()));
    if (!vp1){
        Printf("AliVParticle associated to Leading constituent not found");
        return -1;
    }
    
    vp2 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[1], jetCont->GetParticleContainer()->GetArray()));
    if (!vp2){
        Printf("AliVParticle associated to Subleading constituent not found");
        return -1;
    }
    
    
    num=vp1->Pt();
    den=vp2->Pt();
    //Printf("vp1->Pt() =%f, vp2->Pt() =%f", vp1->Pt(), vp2->Pt());
    
    return num-den;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb =0) {
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedLeSub();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub();
        else
            return LeSub(jet, jetContNb);
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb=0){
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedConstituent();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedConstituent();
        else
            return jet->GetNumberOfTracks();
    
}


//______________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::Sigma2(AliEmcalJet *jet, Int_t jetContNb=0){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t mxx    = 0.;
    Double_t myy    = 0.;
    Double_t mxy    = 0.;
    int  nc     = 0;
    Double_t sump2  = 0.;
    
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        Double_t ppt=vp1->Pt();
        Double_t dphi = RelativePhi(vp1->Phi(),jet->Phi());
        
        Double_t deta = vp1->Eta()-jet->Eta();
        mxx += ppt*ppt*deta*deta;
        myy += ppt*ppt*dphi*dphi;
        mxy -= ppt*ppt*deta*TMath::Abs(dphi);
        nc++;
        sump2 += ppt*ppt;
        
    }
    if(nc<2) return 0;
    if(sump2==0) return 0;
    // Sphericity Matrix
    Double_t ele[4] = {mxx , mxy , mxy , myy };
    TMatrixDSym m0(2,ele);
    
    // Find eigenvectors
    TMatrixDSymEigen m(m0);
    TVectorD eval(2);
    TMatrixD evecm = m.GetEigenVectors();
    eval  = m.GetEigenValues();
    // Largest eigenvector
    int jev = 0;
    //  cout<<eval[0]<<" "<<eval[1]<<endl;
    if (eval[0] < eval[1]) jev = 1;
    TVectorD evec0(2);
    // Principle axis
    evec0 = TMatrixDColumn(evecm, jev);
    Double_t compx=evec0[0];
    Double_t compy=evec0[1];
    TVector2 evec(compx, compy);
    Double_t sig=0;
    if(jev==1) sig=TMath::Sqrt(TMath::Abs(eval[0])/sump2);
    if(jev==0) sig=TMath::Sqrt(TMath::Abs(eval[1])/sump2);
    
    return sig;
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetSigma2(AliEmcalJet *jet, Int_t jetContNb=0){
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedSigma2();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedSigma2();
        else
            return Sigma2(jet, jetContNb);
    
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalHfeTagging::SelectTrigger(Float_t minpT, Float_t maxpT){
    
    AliTrackContainer *PartCont = NULL;
    AliParticleContainer *PartContMC = NULL;
    
    
    if (fJetShapeSub==kConstSub){
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) PartContMC = GetParticleContainer(1);
        else PartCont = GetTrackContainer(1);
    }
    else{
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) PartContMC = GetParticleContainer(0);
        else PartCont = GetTrackContainer(0);
    }
    
    TClonesArray *TracksArray = NULL;
    TClonesArray *TracksArrayMC = NULL;
    
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) TracksArrayMC = PartContMC->GetArray();
    else TracksArray = PartCont->GetArray();
    
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly){
        if(!PartContMC || !TracksArrayMC) return -99999;
    }
    else {
        if(!PartCont || !TracksArray) return -99999;
    }
    
    
    AliAODTrack *Track = 0x0;
    
    
    
    //TList *trackList = new TList();
    Int_t triggers[100];
    for (Int_t iTrigger=0; iTrigger<100; iTrigger++) triggers[iTrigger] = 0;
    Int_t iTT = 0;
    Int_t NTracks=0;
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) NTracks = TracksArrayMC->GetEntriesFast();
    else NTracks = TracksArray->GetEntriesFast();
    
    for(Int_t i=0; i < NTracks; i++){
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly){
            if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
                if (!Track) continue;
                if(TMath::Abs(Track->Eta())>0.9) continue;
                if (Track->Pt()<0.15) continue;
                if ((Track->Pt() >= minpT) && (Track->Pt()< maxpT)) {
                    triggers[iTT] = i;
                    iTT++;
                }
            }
        }
        else{
            if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
                if (!Track) continue;
                if(TMath::Abs(Track->Eta())>0.9) continue;
                if (Track->Pt()<0.15) continue;
                if ((Track->Pt() >= minpT) && (Track->Pt()< maxpT)) {
                    triggers[iTT] = i;
                    iTT++;
                }
            }
        }
    }
    
    
    if (iTT == 0) return -99999;
    Int_t nbRn = 0, index = 0 ;
    TRandom3* random = new TRandom3(0);
    nbRn = random->Integer(iTT);
    index = triggers[nbRn];
    //Printf("iTT Total= %d, nbRn = %d, Index = %d",iTT, nbRn, index );
    return index;
    
}

//__________________________________________________________________________________
Double_t AliAnalysisTaskEmcalHfeTagging::RelativePhi(Double_t mphi,Double_t vphi){
    
    if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
    else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
    if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
    else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
    double dphi = mphi-vphi;
    if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
    else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
    return dphi;//dphi in [-Pi, Pi]
}


//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::RetrieveEventObjects() {
    //
    // retrieve event objects
    //
    if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
        return kFALSE;
    
    return kTRUE;
}

//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::InclElecTrackCuts(const AliVVertex *pVietx,AliAODTrack *ietrack,Int_t nMother, Double_t listMother[])
{
    // track cuts for inclusive electrons
    
    if(!ietrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
    if(TMath::Abs(ietrack->Eta())>fEtaCut) return kFALSE;
    if(ietrack->GetTPCNcls() < fTPCnCut) return kFALSE;
    if (ietrack->GetITSNcls() < fITSncut) return kFALSE;
    if(!ietrack->IsOn(AliAODTrack::kITSrefit)) return kFALSE;
    if(!ietrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
    if(!(ietrack->HasPointOnITSLayer(0) && ietrack->HasPointOnITSLayer(1))) return kFALSE;
    if(ietrack->GetTPCFoundFraction() < 0.6) return kFALSE;
    if(ietrack->Chi2perNDF() > 4) return kFALSE;
    
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < nMother; kinkmother++) {
        if(ietrack->GetID() == listMother[kinkmother]) {
            kinkmotherpass = kFALSE;
            continue;
        }
    }
    if(!kinkmotherpass) return kFALSE;
    
    Double_t d0z0[2]={-999,-999}, cov[3];
    
    if(ietrack->PropagateToDCA(pVietx, fVevent->GetMagneticField(), 20., d0z0, cov))
        if(TMath::Abs(d0z0[0]) > fDcaXYcut || TMath::Abs(d0z0[1]) > fDcaZcut) return kFALSE;
    
    return kTRUE;
    
}
//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::PhotElecTrackCuts(const AliVVertex *pVaetx,AliAODTrack *aetrack,Int_t nMother, Double_t listMother[])
{
    // track cuts for associate tracks of photonic electrons
    
    if(!aetrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) return kFALSE;
    if(aetrack->Pt() < fAssPtCut) return kFALSE;
    if(TMath::Abs(aetrack->Eta())>0.9) return kFALSE;
    if(aetrack->GetTPCNcls() < fAssTPCnCut) return kFALSE;
    if (fAssITSrefitCut && !(aetrack->IsOn(AliAODTrack::kITSrefit))) return kFALSE;
    if(!aetrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
    
    return kTRUE;
    
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalHfeTagging::Terminate(Option_t *)
{
    // Called once at the end of the analysis.
    
    // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
    // if (!fTreeObservableTagging){
    //   Printf("ERROR: fTreeObservableTagging not available");
    //   return;
    // }
    
}

