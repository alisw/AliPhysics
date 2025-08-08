//
// HFE tagged jet shape analysis task.
//
// Authors: D. Caffarri, L. Cunqueiro (jet); D. Godoy (HFE), L. Barreto
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
fMinM02cut(0.01),
fMaxM02cut(2.00),
fMinPtTPC(0.5),
fMaxPtTPC(4.),
fMinPtEMCal(4.),
fMaxPtEMCal(25.),
fMinPtSemiInclusive(5.),
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
fnTrueElecPerJetPt(0),
fnTrueHFElecPerJetPt(0),
fnTruePElecPerJetPt(0),
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
fptJetHadron(0),
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
fnEovPelecEMCalcut(0x0),
fnEovPelecTPCEMCalcut(0x0),
fnM20backg(0x0),
fnM02backg(0x0),
fnEovPbackg(0x0),
fnEovPbackgEMCalcut(0x0),
fnClsE(0x0),
fnM20(0x0),
fnM02(0x0),
fnClsTime(0x0),
fnM20TrueElec(0x0),
fnM02TrueElec(0x0),
fnShowerShapeTrueElec(0x0),
fnEovPTrueElecnocut(0x0),
fnEovPTrueElecTPCcut(0x0),
fnEovPTrueElecEMCalcut(0x0),
fnEovPTrueElecTPCEMCalcut(0x0),
fnM20TrueBkg(0x0),
fnM02TrueBkg(0x0),
fnShowerShapeTrueBkg(0x0),
fnEovPTrueBkgnocut(0x0),
fnEovPTrueBkgTPCcut(0x0),
fnEovPTrueBkgEMCalcut(0x0),
fnEovPTrueBkgTPCEMCalcut(0x0),
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
fptJetNoElectrons(0),
fAngJetNoElectrons(0),
fDispJetNoElectrons(0),
fTreeObservableTagging(0)
{
    for(Int_t i = 0; i < nbranches; i++){
        fShapesVar[i]=0;
    }
    
	for(Int_t i = 0; i < ncutseff; i++){
        fptTrueHFEeffTPCTOF[i] = NULL;
        fptTrueHFEeffEMCal[i] = NULL;
    }
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        fnEovPelecTPCsscut[i] = NULL;
    }
    
    for(Int_t i = 0; i < ncutseffsubs; i++){
        fptTrueHFEeffTPCTOFang[i] = NULL;
        fptTrueHFEeffEMCalang[i] = NULL;
        fptTrueHFEeffTPCTOFdisp[i] = NULL;
        fptTrueHFEeffEMCaldisp[i] = NULL;

    }
    
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        fInvmassLS[i] = NULL;
        fInvmassULS[i] = NULL;
        for(Int_t j = 0; j < nbins_eptTPC; j++){
            fnTPCSigma[i][j] = NULL;
        }
    }
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        for(Int_t j = 0; j < nbins_gsys; j++){
            fRecPEAng[i][j] = NULL;
            fTotPEAng[i][j] = NULL;
            fULSptAng[i][j] = NULL;
            fLSptAng[i][j] = NULL;
        }
        for(Int_t j = 0; j < nbins_ptdsys; j++){
            fRecPEDisp[i][j] = NULL;
            fTotPEDisp[i][j] = NULL;
            fULSptDisp[i][j] = NULL;
            fLSptDisp[i][j] = NULL;
        }
    }
    
	for(Int_t i = 0; i < nbins_jetpt; i++){
        fPtSemiInclJet[i] = NULL;
        fAngSemiInclJet[i] = NULL;
		fDispSemiInclJet[i] = NULL;
		fRMPtSemiInclJet[i] = NULL; 
		fRMAngSemiInclJet[i] = NULL; 
		fRMDispSemiInclJet[i] = NULL; 
		fRMUWPtSemiInclJet[i] = NULL; 
		fRMUWAngSemiInclJet[i] = NULL; 
		fRMUWDispSemiInclJet[i] = NULL; 
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
fMinM02cut(0.01),
fMaxM02cut(2.00),
fMinPtTPC(0.5),
fMaxPtTPC(4.),
fMinPtEMCal(4.),
fMaxPtEMCal(25.),
fMinPtSemiInclusive(5.),
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
fnTrueElecPerJetPt(0),
fnTrueHFElecPerJetPt(0),
fnTruePElecPerJetPt(0),
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
fptJetHadron(0),
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
fnEovPelecEMCalcut(0x0),
fnEovPelecTPCEMCalcut(0x0),
fnM20backg(0x0),
fnM02backg(0x0),
fnEovPbackg(0x0),
fnEovPbackgEMCalcut(0x0),
fnClsE(0x0),
fnM20(0x0),
fnM02(0x0),
fnClsTime(0x0),
fnM20TrueElec(0x0),
fnM02TrueElec(0x0),
fnShowerShapeTrueElec(0x0),
fnEovPTrueElecnocut(0x0),
fnEovPTrueElecTPCcut(0x0),
fnEovPTrueElecEMCalcut(0x0),
fnEovPTrueElecTPCEMCalcut(0x0),
fnM20TrueBkg(0x0),
fnM02TrueBkg(0x0),
fnShowerShapeTrueBkg(0x0),
fnEovPTrueBkgnocut(0x0),
fnEovPTrueBkgTPCcut(0x0),
fnEovPTrueBkgEMCalcut(0x0),
fnEovPTrueBkgTPCEMCalcut(0x0),
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
fptJetNoElectrons(0),
fAngJetNoElectrons(0),
fDispJetNoElectrons(0),
fTreeObservableTagging(0)
{
    // Standard constructor.
    for(Int_t i = 0; i < nbranches; i++){
        fShapesVar[i]=0;
    }
    
	for(Int_t i = 0; i < ncutseff; i++){
        fptTrueHFEeffTPCTOF[i] = NULL;
        fptTrueHFEeffEMCal[i] = NULL;
    }
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        fnEovPelecTPCsscut[i] = NULL;
    }
    
    for(Int_t i = 0; i < ncutseffsubs; i++){
        fptTrueHFEeffTPCTOFang[i] = NULL;
        fptTrueHFEeffEMCalang[i] = NULL;
        fptTrueHFEeffTPCTOFdisp[i] = NULL;
        fptTrueHFEeffEMCaldisp[i] = NULL;
    }
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        fInvmassLS[i] = NULL;
        fInvmassULS[i] = NULL;
        for(Int_t j=0;j<18;j++){
            fnTPCSigma[i][j] = NULL;
        }
    }
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        for(Int_t j = 0; j < nbins_gsys; j++){
            fRecPEAng[i][j] = NULL;
            fTotPEAng[i][j] = NULL;
            fULSptAng[i][j] = NULL;
            fLSptAng[i][j] = NULL;
        }
        for(Int_t j = 0; j < nbins_ptdsys; j++){
            fRecPEDisp[i][j] = NULL;
            fTotPEDisp[i][j] = NULL;
            fULSptDisp[i][j] = NULL;
            fLSptDisp[i][j] = NULL;
        }
    }
	
	for(Int_t i = 0; i < nbins_jetpt; i++){
        fPtSemiInclJet[i] = NULL;
        fAngSemiInclJet[i] = NULL;
		fDispSemiInclJet[i] = NULL;
		fRMPtSemiInclJet[i] = NULL; 
		fRMAngSemiInclJet[i] = NULL; 
		fRMDispSemiInclJet[i] = NULL; 
		fRMUWPtSemiInclJet[i] = NULL; 
		fRMUWAngSemiInclJet[i] = NULL; 
		fRMUWDispSemiInclJet[i] = NULL; 
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
    
	/* const int nbins_ept = 33, nbins_jetpt = 6, nbins_g = 14, nbins_ptd = 14; */
	/* const int nbins_gsys = 6, nbins_ptdsys = 7, nbins_ptauxiliary = 59; */

    Double_t bin_ept[nbins_ept + 1] = {0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4,
        1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 5.,
        6., 8., 10., 12., 14., 16., 19., 22., 25., 30.,
        35., 40., 45., 50.};
    
    Double_t bin_jetpt[nbins_jetpt + 1] = {5., 10., 20., 40., 60., 80., 120.};
    Double_t bin_g[nbins_g + 1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
								   0.5, 0.55, 0.6, 0.65, 0.7};
    Double_t bin_ptd[nbins_ptd + 1] = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 
									   0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    
    Double_t bin_ptauxiliary[nbins_ptauxiliary + 1] = {0.01,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,
        0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,
        0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,
        1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,
        3.6,3.8,4,4.5,5,5.5,6,6.5,7,8,
        9,10,11,12,13,14,15,16,18,20};
    
    /* Double_t bin_gsys[nbins_gsys + 1] = {0., 0.1, 0.2, 0.3, 0.4, 0.6, 0.7}; */
    /* Double_t bin_ptdsys[nbins_gsys + 1] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; */
    
	Int_t nbin_jetptRM[2] = {24, 24};
	Double_t max_jetptRM[2] = {120., 120};
   	Double_t min_jetptRM[2] = {0., 0.};

	Int_t nbin_gRM[4] = {24, 24, 20, 20};
	Double_t max_gRM[4] = {120., 120., 1., 1.};
	Double_t min_gRM[4] = {0., 0., 0., 0.};
	
	Int_t nbin_ptdRM[4] = {24, 24, 20, 20};
	Double_t max_ptdRM[4] = {120., 120., 1., 1.};
	Double_t min_ptdRM[4] = {0., 0., 0., 0.};


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
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        for(Int_t j=0; j < nbins_eptTPC; j++){
            fnTPCSigma[i][j] = new TH1F(Form("fnTPCSigma%d%d",i,j),Form("fnTPCSigma%d%d",i,j), 100,-15,15);
            fOutput->Add(fnTPCSigma[i][j]);
        }
    }
    
    fnULSmLSpairsPerElectron =new TH2F("fnULSmLSpairsPerElectron", "fnULSmLSpairsPerElectron", 50, 0, 5, 100, 0, 100);
    fOutput->Add(fnULSmLSpairsPerElectron);
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        fInvmassLS[i] = new TH2F(Form("fInvmassLS%d",i), Form("fInvmassLS%d",i), nbins_ept,bin_ept, 100, 0, 0.5);
        fOutput->Add(fInvmassLS[i]);
    }
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        fInvmassULS[i] = new TH2F(Form("fInvmassULS%d",i), Form("fInvmassULS%d",i), nbins_ept,bin_ept, 100, 0, 0.5);
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
    
	fnTrueElecPerJetPt = new TH2F("fnTrueElecPerJetPt", "fnTrueElecPerJetPt", nbins_jetpt, bin_jetpt, 10, 0, 10);
    fOutput->Add(fnTrueElecPerJetPt);
    
    fnTrueHFElecPerJetPt = new TH2F("fnTrueHFElecPerJetPt", "fnTrueHFElecPerJetPt", nbins_jetpt, bin_jetpt, 10, 0, 10);
    fOutput->Add(fnTrueHFElecPerJetPt);
    
	fnTruePElecPerJetPt = new TH2F("fnTruePElecPerJetPt", "fnTruePElecPerJetPt", nbins_jetpt, bin_jetpt, 10, 0, 10);
    fOutput->Add(fnTruePElecPerJetPt);
    
    fPi0PtGen = new TH1F("fPi0PtGen","fPi0PtGen",nbins_ptauxiliary,bin_ptauxiliary);
    fOutput->Add(fPi0PtGen);
    
    fPi0PtEnh = new TH1F("fPi0PtEnh","fPi0PtEnh",nbins_ptauxiliary,bin_ptauxiliary);
    fOutput->Add(fPi0PtEnh);
    
    fEtaPtGen = new TH1F("fEtaPtGen", "fEtaPtGen",nbins_ptauxiliary,bin_ptauxiliary);
    fOutput->Add(fEtaPtGen);
    
    fEtaPtEnh = new TH1F("fEtaPtEnh", "fEtaPtEnh",nbins_ptauxiliary,bin_ptauxiliary);
    fOutput->Add(fEtaPtEnh);
    
    fGenHfePt = new TH1F("fGenHfePt","fGenHfePt",nbins_ptauxiliary,bin_ptauxiliary);
    fOutput->Add(fGenHfePt);
    
    fGenPePt = new TH1F("fGenPePt", "fGenPePt",nbins_ptauxiliary,bin_ptauxiliary);
    fOutput->Add(fGenPePt);
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        for(Int_t j = 0; j < nbins_gsys; j++){
            fRecPEAng[i][j] = new TH1F(Form("fRecPEAng%d%d",i,j), Form("fRecPEAng%d%d",i,j), nbins_ept,bin_ept);
            fOutput->Add(fRecPEAng[i][j]);

            fTotPEAng[i][j] = new TH1F(Form("fTotPEAng%d%d",i,j), Form("fTotPEAng%d%d",i,j), nbins_ept,bin_ept);
            fOutput->Add(fTotPEAng[i][j]);
            
            fULSptAng[i][j] = new TH1F(Form("fULSptAng%d%d",i,j), Form("fULSptAng%d%d",i,j), nbins_ept,bin_ept);
            fOutput->Add(fULSptAng[i][j]);

            fLSptAng[i][j] = new TH1F(Form("fLSptAng%d%d",i,j), Form("fLSptAng%d%d",i,j), nbins_ept,bin_ept);
            fOutput->Add(fLSptAng[i][j]);
        }
        
        for(Int_t j = 0; j < nbins_ptdsys; j++){
            fRecPEDisp[i][j] = new TH1F(Form("fRecPEDisp%d%d",i,j), Form("fRecPEDisp%d%d",i,j), nbins_ept,bin_ept);
            fOutput->Add(fRecPEDisp[i][j]);
            
            fTotPEDisp[i][j] = new TH1F(Form("fTotPEDisp%d%d",i,j), Form("fTotPEDisp%d%d",i,j), nbins_ept,bin_ept);
            fOutput->Add(fTotPEDisp[i][j]);
            
            fULSptDisp[i][j] = new TH1F(Form("fULSptDisp%d%d",i,j), Form("fULSptDisp%d%d",i,j), nbins_ept,bin_ept);
            fOutput->Add(fULSptDisp[i][j]);
            
            fLSptDisp[i][j] = new TH1F(Form("fLSptDisp%d%d",i,j), Form("fLSptDisp%d%d",i,j), nbins_ept,bin_ept);
            fOutput->Add(fLSptDisp[i][j]);
        }
    }
    
    fPtP=new TH2F("fPtP", "fPtP", nbins_ept,bin_ept,nbins_ept,bin_ept);
    fOutput->Add(fPtP);
    
    fptJetIE= new TH1F("fptJetIE", "fptJetIE", 100, 0, 200);
    fOutput->Add(fptJetIE);
    
    fptJetPE= new TH1F("fptJetPE", "fptJetPE", 100, 0, 200);
    fOutput->Add(fptJetPE);
    
    fptJetHFE= new TH1F("fptJetHFE", "fptJetHFE", 100, 0, 200);
    fOutput->Add(fptJetHFE);
	
	fptJetHadron = new TH1F("fptJetHadron", "fptJetHadron", 24, 0., 120.);
    fOutput -> Add(fptJetHadron);
    
    fptRecPE= new TH1F("fptRecPE", "fptRecPE", nbins_ept,bin_ept);
    fOutput->Add(fptRecPE);
    
    fptTruePE= new TH1F("fptTruePE", "fptTruePE", nbins_ept,bin_ept);
    fOutput->Add(fptTruePE);
    
    for(Int_t i = 0; i < ncutseff; i++){
        fptTrueHFEeffTPCTOF[i] = new TH1F(Form("fptTrueHFEeffTPCTOF%d",i), Form("fptTrueHFEeffTPCTOF%d",i), nbins_ept,bin_ept);
        fOutput->Add(fptTrueHFEeffTPCTOF[i]);
    }
    
    for(Int_t i = 0; i < ncutseffsubs; i++){
        fptTrueHFEeffTPCTOFang[i] = new TH3F(Form("fptTrueHFEeffTPCTOFang%d",i),
                                             Form("fptTrueHFEeffTPCTOFang%d",i), nbins_ept,bin_ept,nbins_jetpt, bin_jetpt, nbins_g, bin_g);
        fOutput->Add(fptTrueHFEeffTPCTOFang[i]);
    }
    
    for(Int_t i = 0; i < ncutseffsubs; i++){
        fptTrueHFEeffTPCTOFdisp[i] = new TH3F(Form("fptTrueHFEeffTPCTOFdisp%d",i),
                                              Form("fptTrueHFEeffTPCTOFdisp%d",i), nbins_ept,bin_ept,nbins_jetpt, bin_jetpt,nbins_ptd, bin_ptd);
        fOutput->Add(fptTrueHFEeffTPCTOFdisp[i]);
    }
    
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        fptTrueHFEeffEMCal[i] = new TH1F(Form("fptTrueHFEeffEMCal%d",i), Form("fptTrueHFEeffEMCal%d",i), nbins_ept,bin_ept);
        fOutput->Add(fptTrueHFEeffEMCal[i]);
    }
    
    for(Int_t i = 0; i < ncutseffsubs; i++){
        fptTrueHFEeffEMCalang[i] = new TH3F(Form("fptTrueHFEeffEMCalang%d",i),
                                            Form("fptTrueHFEeffEMCalang%d",i), nbins_ept,bin_ept,nbins_jetpt, bin_jetpt, nbins_g, bin_g);
        fOutput->Add(fptTrueHFEeffEMCalang[i]);
    }
    
    for(Int_t i = 0; i < ncutseffsubs; i++){
        fptTrueHFEeffEMCaldisp[i] = new TH3F(Form("fptTrueHFEeffEMCaldisp%d",i),
                                             Form("fptTrueHFEeffEMCaldisp%d",i), nbins_ept,bin_ept,nbins_jetpt, bin_jetpt,nbins_ptd, bin_ptd);
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
    
	fnEovPelecEMCalcut = new TH2F("fnEovPelecEMCalcut", "fnEovPelecEMCalcut", 100, 0, 100, 40, 0, 2);
    fOutput->Add(fnEovPelecEMCalcut);
    
    fnEovPelecTPCEMCalcut = new TH2F("fnEovPelecTPCEMCalcut", "fnEovPelecTPCEMCalcut", 100, 0, 100,40,0,2);
    fOutput->Add(fnEovPelecTPCEMCalcut);
    
    for(Int_t i = 0; i < nbins_jetpt; i++){
        fnEovPelecTPCsscut[i] = new TH2F(Form("fnEovPelecTPCsscut%d",i), Form("fnEovPelecTPCsscut%d",i), nbins_ept,bin_ept,40,0,2);
        fOutput->Add(fnEovPelecTPCsscut[i]);
    }
    
	fnM20backg = new TH2F("fnM20backg", "fnM20backg", 50, 0.5, 50, 100, 0, 2);
    fOutput->Add(fnM20backg);
	
	fnM02backg = new TH2F("fnM02backg", "fnM02backg", 50, 0.5, 50, 100, 0, 2);
    fOutput->Add(fnM02backg);
    
    fnEovPbackg = new TH2F("fnEovPbackg", "fnEovPbackg", 100, 0, 100,40,0,2);
    fOutput->Add(fnEovPbackg);
    
	fnEovPbackgEMCalcut = new TH2F("fnEovPbackgEMCalcut", "fnEovPbackgEMCalcut", 100, 0, 100, 40, 0, 2);
    fOutput->Add(fnEovPbackgEMCalcut);
    
    fnClsE = new TH2F("fnClsE", "fnClsE", 100, 0, 100, 100, 0,100);
    fOutput->Add(fnClsE);
    
    fnM20 = new TH2F("fnM20", "fnM20", 100, 0, 100, 100, 0,2);
    fOutput->Add(fnM20);
    
    fnM02 = new TH2F("fnM02", "fnM02", 100, 0, 100, 100, 0,2);
    fOutput->Add(fnM02);
    
    fnClsTime = new TH2F("fnClsTime", "fnClsTime", 100, 0, 100, 100, -200,200);
    fOutput->Add(fnClsTime);
    
	fnM20TrueElec = new TH2F("fnM20TrueElec", "fnM20TrueElec", 50, 0., 50., 100, 0.01, 2.01);
    fOutput->Add(fnM20TrueElec);
	
	fnM02TrueElec = new TH2F("fnM02TrueElec", "fnM02TrueElec", 50, 0., 50., 100, 0.01, 2.01);
    fOutput->Add(fnM02TrueElec);
	
	fnShowerShapeTrueElec = new TH3F("fnShowerShapeTrueElec", "fnShowerShapeTrueElec", 50, 0., 50., 50, 0.01, 2.01, 50, 0.01, 2.01);
    fOutput->Add(fnShowerShapeTrueElec);
	
	fnEovPTrueElecnocut = new TH2F("fnEovPTrueElecnocut", "fnEovPTrueElecnocut", 50, 0., 50., 40, 0, 2);
    fOutput->Add(fnEovPTrueElecnocut);
	
	fnEovPTrueElecTPCcut = new TH2F("fnEovPTrueElecTPCcut", "fnEovPTrueElecTPCcut", 50, 0., 50., 40, 0, 2);
    fOutput->Add(fnEovPTrueElecTPCcut);
	
	fnEovPTrueElecEMCalcut = new TH2F("fnEovPTrueElecEMCalcut", "fnEovPTrueElecEMCalcut", 50, 0., 50., 40, 0, 2);
    fOutput->Add(fnEovPTrueElecEMCalcut);
    
	fnEovPTrueElecTPCEMCalcut = new TH2F("fnEovPTrueElecTPCEMCalcut", "fnEovPTrueElecTPCEMCalcut", 50, 0., 50., 40, 0, 2);
    fOutput->Add(fnEovPTrueElecTPCEMCalcut);
	
	fnM20TrueBkg = new TH2F("fnM20TrueBkg", "fnM20TrueBkg", 50, 0., 50., 100, 0.01, 2.01);
    fOutput->Add(fnM20TrueBkg);
	
	fnM02TrueBkg = new TH2F("fnM02TrueBkg", "fnM02TrueBkg", 50, 0., 50., 100, 0.01, 2.01);
    fOutput->Add(fnM02TrueBkg);
	
	fnShowerShapeTrueBkg = new TH3F("fnShowerShapeTrueBkg", "fnShowerShapeTrueBkg", 50, 0., 50., 50, 0.01, 2.01, 50, 0.01, 2.01);
    fOutput->Add(fnShowerShapeTrueBkg);
	
	fnEovPTrueBkgnocut = new TH2F("fnEovPTrueBkgnocut", "fnEovPTrueBkgnocut", 50, 0., 50., 40, 0, 2);
    fOutput->Add(fnEovPTrueBkgnocut);

	fnEovPTrueBkgTPCcut = new TH2F("fnEovPTrueBkgTPCcut", "fnEovPTrueBkgTPCcut", 50, 0., 50., 40, 0, 2);
    fOutput->Add(fnEovPTrueBkgTPCcut);
	
	fnEovPTrueBkgEMCalcut = new TH2F("fnEovPTrueBkgEMCalcut", "fnEovPTrueBkgEMCalcut", 50, 0., 50., 40, 0, 2);
    fOutput->Add(fnEovPTrueBkgEMCalcut);

	fnEovPTrueBkgTPCEMCalcut = new TH2F("fnEovPTrueBkgTPCEMCalcut", "fnEovPTrueBkgTPCEMCalcut", 50, 0., 50., 40, 0, 2);
    fOutput->Add(fnEovPTrueBkgTPCEMCalcut);
    
    fAngULS = new TH2F("fAngULS", "fAngULS", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngULS);
    
    fAngLS = new TH2F("fAngLS", "fAngLS", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngLS);
    
    fAngChargPart = new TH2F("fAngChargPart", "fAngChargPart", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngChargPart);
    
    fAngHadron = new TH2F("fAngHadron", "fAngHadron", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngHadron);
    
    fAngIncElec = new TH2F("fAngIncElec", "fAngIncElec", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngIncElec);
    
    fAngPhotElec = new TH2F("fAngPhotElec", "fAngPhotElec", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngPhotElec);
    
    fAngElecFromD = new TH2F("fAngElecFromD", "fAngElecFromD", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngElecFromD);
    
    fAngElecFromB = new TH2F("fAngElecFromB", "fAngElecFromB", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngElecFromB);
    
    fAngElecFromDFromB = new TH2F("fAngElecFromDFromB", "fAngElecFromDFromB", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngElecFromDFromB);
    
    fAngD = new TH2F("fAngD", "fAngD", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngD);
    
    fAngB = new TH2F("fAngB", "fAngB", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngB);
    
    fAngCharm = new TH2F("fAngCharm", "fAngCharm", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngCharm);
    
    fAngBeauty = new TH2F("fAngBeauty", "fAngBeauty", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngBeauty);
    
    fAngQuark = new TH2F("fAngQuark", "fAngQuark", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngQuark);
    
    fAngGluon = new TH2F("fAngGluon", "fAngGluon", nbins_jetpt, bin_jetpt, nbins_g, bin_g);
    fOutput->Add(fAngGluon);
    
    fDispULS = new TH2F("fDispULS", "fDispULS", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispULS);
    
    fDispLS = new TH2F("fDispLS", "fDispLS", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispLS);
    
    fDispChargPart = new TH2F("fDispChargPart", "fDispChargPart", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispChargPart);
    
    fDispHadron = new TH2F("fDispHadron", "fDispHadron", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispHadron);
    
    fDispIncElec = new TH2F("fDispIncElec", "fDispIncElec", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispIncElec);
    
    fDispPhotElec = new TH2F("fDispPhotElec", "fDispPhotElec", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispPhotElec);
    
    fDispElecFromD = new TH2F("fDispElecFromD", "fDispElecFromD", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispElecFromD);
    
    fDispElecFromB = new TH2F("fDispElecFromB", "fDispElecFromB", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispElecFromB);
    
    fDispElecFromDFromB = new TH2F("fDispElecFromDFromB", "fDispElecFromDFromB", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispElecFromDFromB);
    
    fDispD = new TH2F("fDispD", "fDispD", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispD);
    
    fDispB = new TH2F("fDispB", "fDispB", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispB);
    
    fDispCharm = new TH2F("fDispCharm", "fDispCharm", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispCharm);
    
    fDispBeauty = new TH2F("fDispBeauty", "fDispBeauty", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispBeauty);
    
    fDispQuark = new TH2F("fDispQuark", "fDispQuark", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispQuark);
    
    fDispGluon = new TH2F("fDispGluon", "fDispGluon", nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
    fOutput->Add(fDispGluon);

	// Inclusive and semi-inclusive histograms
    for(Int_t i = 0; i < nbins_jetpt; i++){
        fPtSemiInclJet[i] = new TH1F(Form("fPtSemiInclJet%d",i), Form("fPtSemiInclJet%d",i), 24, 0., 120.);
		fRMPtSemiInclJet[i] = new THnSparseF(Form("fRMPtSemiInclJet%d", i), Form("fRMPtSemiInclJet%d", i), 2, 
						                    nbin_jetptRM, min_jetptRM, max_jetptRM);
		fRMUWPtSemiInclJet[i] = new THnSparseF(Form("fRMUWPtSemiInclJet%d", i), Form("fRMUWPtSemiInclJet%d", i), 2, 
						                    nbin_jetptRM, min_jetptRM, max_jetptRM);
        fOutput->Add(fPtSemiInclJet[i]);
        fOutput->Add(fRMPtSemiInclJet[i]);
        fOutput->Add(fRMUWPtSemiInclJet[i]);

        fAngSemiInclJet[i] = new TH2F(Form("fAngSemiInclJet%d",i), Form("fAngSemiInclJet%d",i), nbins_jetpt, bin_jetpt, nbins_g, bin_g);
		fRMAngSemiInclJet[i] = new THnSparseF(Form("fRMAngSemiInclJet%d", i), Form("fRMAngSemiInclJet%d", i), 4, 
						                    nbin_gRM, min_gRM, max_gRM);
		fRMUWAngSemiInclJet[i] = new THnSparseF(Form("fRMUWAngSemiInclJet%d", i), Form("fRMUWAngSemiInclJet%d", i), 4, 
						                    nbin_gRM, min_gRM, max_gRM);
        fOutput->Add(fAngSemiInclJet[i]);
        fOutput->Add(fRMAngSemiInclJet[i]);
        fOutput->Add(fRMUWAngSemiInclJet[i]);
        
		fDispSemiInclJet[i] = new TH2F(Form("fDispSemiInclJet%d",i), Form("fDispSemiInclJet%d",i), nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
		fRMDispSemiInclJet[i] = new THnSparseF(Form("fRMDispSemiInclJet%d", i), Form("fRMDispSemiInclJet%d", i), 4, 
						                    nbin_ptdRM, min_ptdRM, max_ptdRM);
		fRMUWDispSemiInclJet[i] = new THnSparseF(Form("fRMUWDispSemiInclJet%d", i), Form("fRMUWDispSemiInclJet%d", i), 4, 
						                    nbin_ptdRM, min_ptdRM, max_ptdRM);
        fOutput->Add(fDispSemiInclJet[i]);
        fOutput->Add(fRMDispSemiInclJet[i]);
        fOutput->Add(fRMUWDispSemiInclJet[i]);
    }

	// No electrons 3D methodology
    fptJetNoElectrons = new TH2F("fptJetNoElectrons", "fptJetNoElectrons", nbins_ept, bin_ept, nbins_jetpt, bin_jetpt);
	fOutput -> Add(fptJetNoElectrons);

    fAngJetNoElectrons = new TH3F("fAngJetNoElectrons", "fAngJetNoElectrons", nbins_ept, bin_ept, nbins_jetpt, bin_jetpt, nbins_g, bin_g);
	fOutput -> Add(fAngJetNoElectrons);

    fDispJetNoElectrons = new TH3F("fDispJetNoElectrons", "fDispJetNoElectrons", nbins_ept, bin_ept, nbins_jetpt, bin_jetpt, nbins_ptd, bin_ptd);
	fOutput -> Add(fDispJetNoElectrons);
    

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
    
    fTreeObservableTagging = new TTree(GetOutputSlot(2)->GetContainer()->GetName(), GetOutputSlot(2)->GetContainer()->GetName());
    
    TString *fShapesVarNames = new TString [nbranches];
    
	// Jet info
    fShapesVarNames[0] = "partonCode";
    fShapesVarNames[1] = "ptJet";
    fShapesVarNames[2] = "ptDJet";
    fShapesVarNames[3] = "mJet";
    fShapesVarNames[4] = "angularity";
    fShapesVarNames[5] = "ptJetMatch";
    fShapesVarNames[6] = "ptDJetMatch";
    fShapesVarNames[7] = "mJetMatch";
    fShapesVarNames[8] = "angularityMatch";
    fShapesVarNames[9] = "weightPythia";
    fShapesVarNames[10] = "rhoVal";
    fShapesVarNames[11] = "rhoMassVal";
    fShapesVarNames[12] = "ptUnsub";

	// Electron info
    fShapesVarNames[13] = "ptTrueHFE";
    fShapesVarNames[14] = "ptElec";
    fShapesVarNames[15] = "nInclElec";
    fShapesVarNames[16] = "nPhotElec";
    fShapesVarNames[17] = "hasElec";
    fShapesVarNames[18] = "nTrueElectronsMC";
    fShapesVarNames[19] = "nTrueHFElecMC";
    fShapesVarNames[20] = "pdgOnlyOneTrueElectron";
    
    
    for(Int_t ivar=0; ivar < nbranches; ivar++){
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

	// Fix user-defined minimum pT cut for semi-inclusive if needed as calculations need
	// its array to be ordered in MaxPtBinForSemiInclusiveJet
	if (fMinPtSemiInclusive <= 2.5) {
		cout << "User-defined MinPtSemiInclusive <= 2.5 GeV, setting new value to 5 GeV" << endl;
		fMinPtSemiInclusive = 5.;
	}
    
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

	// Always load all possible 4 containers
	// Order will change depending on JetShapeType and Sub, but main is always first
    AliJetContainer* jetCont = GetJetContainer(0);
    if (!jetCont) return 0;
	
    AliJetContainer* jetContExtra1 = GetJetContainer(1);
    AliJetContainer* jetContExtra2 = GetJetContainer(2);
    AliJetContainer* jetContExtra3 = GetJetContainer(3);
    
	AliJetContainer* jetContTrue = nullptr;
	AliJetContainer* jetContUS = nullptr;
	AliJetContainer* jetContPartLevel = nullptr;
		
	// Determine extra jet containers 
	// No need to do anything if kMCTrue, kData or kGenOnTheFly as they only use first
	if (fJetShapeType == kDetEmbPartPythia || fJetShapeType == kPythiaDef) {
		jetContTrue = jetContExtra1;
		if (fJetShapeSub == kConstSub) {
			jetContUS = jetContExtra2;			
			jetContPartLevel = jetContExtra3;			
		} else {
			jetContPartLevel = jetContExtra2;			
		}
	}
	
	// Determine container for jet matching
	Int_t kMatched = 0;
    if (fJetShapeType == kPythiaDef) {
        if (fJetShapeSub == kConstSub) {
			kMatched = 3;
		} else {
    		kMatched = 1;
		}
	} else if (fJetShapeType == kDetEmbPartPythia) {
        if (fJetShapeSub == kConstSub) {
			kMatched = 3;
		/* } else if (fJetShapeSub == kDerivSub) { */
		} else {
			kMatched = 2;
		}
	}


    AliEmcalJet* jetbase = nullptr;
    
    Float_t kWeight=1;
    if (fCentSelectOn) {
        if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
    }

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
    
    
    if (jetCont) {
        jetCont -> ResetCurrentID();
        
		// Determine Rho if subtracted container
    	Float_t rhoVal = 0., rhoMassVal = 0.;
    	if ((fJetShapeSub==kConstSub) || (fJetShapeSub==kDerivSub)){
    	    //rho
    	    AliRhoParameter* rhoParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoSparseR020"));
    	    if (!rhoParam) {
    	        Printf("%s: Could not retrieve rho %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoName().Data());

			} else {
				rhoVal = rhoParam->GetVal();
			}

    	    //rhom
    	    AliRhoParameter* rhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoMassSparseR020"));
    	    
			if (!rhomParam) {
    	        Printf("%s: Could not retrieve rho_m %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoMassName().Data());

		   	} else {
					rhoMassVal = rhomParam->GetVal();
			}
    	}

        while ((jetbase = jetCont -> GetNextAcceptJet())) {
            if (!jetbase) continue;

            AliEmcalJet* jettrue = nullptr;  // Matched unsubtracted jet
            AliEmcalJet* jetpartlevel = nullptr;  // Matched particle-level jet
            AliEmcalJet* jetUS = nullptr;
									  
			// Fill kinematics histograms
            fPtJet->Fill(jetbase->Pt());
            fPhiJet->Fill(jetbase->Pt(),jetbase->Phi());
            fEtaJet->Fill(jetbase->Pt(),jetbase->Eta());
            if (jetbase->Pt() > 5.) fEtaPhiJet->Fill(jetbase->Phi(),jetbase->Eta());
            fAreaJet->Fill(jetbase->Pt(),jetbase->Area());

			// Jet matching
			// If kConstSub, get match to equivalent unsubtracted jet jetUS by jettrue
			// Particle level jet will be described by jetpartlevel
            if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kDetEmbPartPythia){
                Double_t fractionpt = 0;  // For min pt shared cut

                if (fJetShapeSub == AliAnalysisTaskEmcalHfeTagging::kConstSub){
					// Original logic
                    for(Int_t i = 0; i < jetContUS -> GetNJets(); i++) {
                        if(jetContUS -> GetJet(i) -> GetLabel() == jetbase -> GetLabel()) {
                    		jetUS = jetContUS -> GetJet(i);
							break;
                        }
                    }

					if (!jetUS) {
                    	/* Printf("Equivalent unsubtracted jet (jetUS) not found, returning"); */
						continue;
					}

            		jettrue = GetClosestOnOtherJetContainer(jetUS, jetContTrue);
					fractionpt = GetFractionSharedPtBetweenJets(jetUS, jettrue);

                } else {
            		jettrue = GetClosestOnOtherJetContainer(jetbase, jetContTrue);
					fractionpt = GetFractionSharedPtBetweenJets(jetbase, jettrue);
				}
                
				// Min pt shared cut
				if (fractionpt < fMinFractionShared) continue;
                
                if (!jettrue) {
                    Printf("jettrue does not exist, returning");
                    continue;
                }
                    	    
                jetpartlevel = GetClosestOnOtherJetContainer(jettrue, jetContPartLevel);
                
                if(!jetpartlevel){
                    Printf("jetpartlevel does not exist, returning");
                    continue;
                }
					
                fh2ResponseUW->Fill(jetbase->Pt(),jettrue->Pt());
            }
            
            if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kPythiaDef){
                jetpartlevel = GetClosestOnOtherJetContainer(jetbase, jetContPartLevel);
                if(!jetpartlevel){
                    Printf("jetpartlevel does not exist, returning");
                    continue;
                }
                
                fh2ResponseUW->Fill(jetbase->Pt(),jetpartlevel->Pt());
                
                Double_t probDensDetPart = -999., probDensPartDet = -999.;
                
                if (jetbase->Pt()>0) probDensPartDet = (jetpartlevel->Pt()-jetbase->Pt())/jetbase->Pt();
                if (jetpartlevel->Pt()>0) probDensDetPart = (jetbase->Pt()-jetpartlevel->Pt())/jetpartlevel->Pt();
                
                fJetProbDensityDetPart->Fill(probDensDetPart,jetpartlevel->Pt());
                fJetProbDensityPartDet->Fill(probDensPartDet,jetbase->Pt());
                
            }
            
            fShapesVar[0] = 0.;
            if (fJetShapeType == kGenOnTheFly){
                const AliEmcalPythiaInfo *partonsInfo = 0x0;
                partonsInfo = GetPythiaInfo();
                Double_t jp1=RelativePhi(jetbase->Phi(),partonsInfo->GetPartonPhi6());
                Double_t detap1=(jetbase->Eta())-(partonsInfo->GetPartonEta6());
                kWeight=partonsInfo->GetPythiaEventWeight();
                fh2ResponseW->Fill(jetbase->Pt(),jetbase->Pt(),kWeight);
                
                Float_t dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
                fEtaJetCorr6->Fill(jetbase->Eta(), partonsInfo->GetPartonEta6());
                fPhiJetCorr6->Fill(jetbase->Phi(), partonsInfo->GetPartonPhi6());
                if(dRp1 < fRMatching) {
                    fShapesVar[0] = partonsInfo->GetPartonFlag6();
                    fPtJetCorr ->Fill(partonsInfo->GetPartonPt6(), jetbase->Pt());
                }
                else {
                    jp1=RelativePhi(jetbase->Phi(),partonsInfo->GetPartonPhi7());
                    detap1=(jetbase->Eta())-(partonsInfo->GetPartonEta7());
                    dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
                    fEtaJetCorr7->Fill(jetbase->Eta(), partonsInfo->GetPartonEta7());
                    fPhiJetCorr7->Fill(jetbase->Phi(), partonsInfo->GetPartonPhi7());
                    if(dRp1 < fRMatching) {
                        fShapesVar[0] = partonsInfo->GetPartonFlag7();
                        fPtJetCorr->Fill(partonsInfo->GetPartonPt7(), jetbase->Pt());
                    }
                    else fShapesVar[0]=0;
                }
            }
            
			// Pt subtraction procedure
            Double_t ptSubtracted = 0., jetbaseang = -1., jetbaseptd = -1.;
            if (fJetShapeSub == kConstSub) { 
				ptSubtracted = jetbase -> Pt();

			} else if (fJetShapeSub == kDerivSub)  {
                ptSubtracted = jetbase -> Pt() -GetRhoVal(0) * jetbase -> Area();

            } else if (fJetShapeSub==kNoSub) {
                if ((fJetShapeType == kData) || (fJetShapeType == kDetEmbPartPythia)) ptSubtracted = jetbase -> Pt() - GetRhoVal(0) * jetbase -> Area();

                else if ((fJetShapeType == kPythiaDef) || (fJetShapeType == kMCTrue) || (fJetShapeType == kGenOnTheFly)) ptSubtracted = jetbase -> Pt();
            }
            
            if (ptSubtracted < fPtThreshold) continue;
            
            if (fOneConstSelectOn == kTRUE) fNbOfConstvspT->Fill(GetJetNumberOfConstituents(jetbase,0), ptSubtracted);
            if ((fCentSelectOn == kFALSE) && (jetbase->GetNumberOfTracks() <= 1)) continue;
            
			jetbaseptd = GetJetpTD(jetbase, 0);
			jetbaseang = GetJetAngularity(jetbase, 0); 

            fShapesVar[1] = ptSubtracted;
            fShapesVar[2] = jetbaseptd;
            fShapesVar[3] = GetJetMass(jetbase, 0);
            fShapesVar[4] = jetbaseang;
            
            Float_t ptMatch=0., ptDMatch=0., massMatch=0., angulMatch=0.;
            
            if (kMatched != 0) {
                ptMatch = jetpartlevel -> Pt();
                ptDMatch = GetJetpTD(jetpartlevel, kMatched);
				massMatch = GetJetMass(jetpartlevel, kMatched);
                angulMatch = GetJetAngularity(jetpartlevel, kMatched);
            }
            
            fShapesVar[5] = ptMatch;
            fShapesVar[6] = ptDMatch;
            fShapesVar[7] = massMatch;
            fShapesVar[8] = angulMatch;
            fShapesVar[9] = kWeight;
            fShapesVar[10] = rhoVal;
            fShapesVar[11] = rhoMassVal;
            fShapesVar[12] = jetbase->Pt();
            
            Int_t nInclusiveElectrons = 0, nPhotonicElectrons = 0, nTrueElectronsMC= 0, nTrueHFElecMC= 0;
            Double_t pElec = 0., ptElec = 0.;
            Bool_t hasElectrons = kFALSE;
            
            GetNumberOfElectrons(jetbase, 0,nMotherKink,listofmotherkink,nInclusiveElectrons,nPhotonicElectrons,pElec,ptElec,hasElectrons);
            
            fAngChargPart->Fill(ptSubtracted, jetbaseang);
            fDispChargPart->Fill(ptSubtracted, jetbaseptd);

			Double_t ptRMvalues[2] = {ptSubtracted, ptMatch};
			Double_t angRMvalues[4] = {ptSubtracted, ptMatch, jetbaseang, angulMatch};
			Double_t ptdRMvalues[4] = {ptSubtracted, ptMatch, jetbaseptd, ptDMatch};

			for (int sibin = MaxPtBinForSemiInclusiveJet(jetbase, 0); sibin >= 0; sibin--) {
				fPtSemiInclJet[sibin] -> Fill(ptSubtracted);
				fAngSemiInclJet[sibin] -> Fill(ptSubtracted, jetbaseang);
				fDispSemiInclJet[sibin] -> Fill(ptSubtracted, jetbaseptd);

				fRMPtSemiInclJet[sibin] -> Fill(ptRMvalues, kWeight);
				fRMUWPtSemiInclJet[sibin] -> Fill(ptRMvalues, 1.);
				
				fRMAngSemiInclJet[sibin] -> Fill(angRMvalues, kWeight);
				fRMUWAngSemiInclJet[sibin] -> Fill(angRMvalues, 1.);
				
				fRMDispSemiInclJet[sibin] -> Fill(ptdRMvalues, kWeight);
				fRMUWDispSemiInclJet[sibin] -> Fill(ptdRMvalues, 1.);
			}	
            
            if (nInclusiveElectrons == 1) {
                fAngIncElec -> Fill(ptSubtracted, jetbaseang);
                fDispIncElec -> Fill(ptSubtracted, jetbaseptd);
            }
            
            if (!hasElectrons) {
                fptJetHadron -> Fill(ptSubtracted);
                fAngHadron -> Fill(ptSubtracted, jetbaseang);
                fDispHadron -> Fill(ptSubtracted, jetbaseptd);

				// Choose one constituent to act as a random tagger
				int randomtaggernumber = gRandom -> Integer(jetbase -> GetNumberOfTracks());
				AliVParticle* randomtagger = jetbase -> Track(randomtaggernumber);
				if (randomtagger) {
					fptJetNoElectrons -> Fill(randomtagger -> Pt(), ptSubtracted);
					fAngJetNoElectrons -> Fill(randomtagger -> Pt(), ptSubtracted, jetbaseang);
					fDispJetNoElectrons -> Fill(randomtagger -> Pt(), ptSubtracted, jetbaseptd);
				}
            }
            
            
            
            // generated HFE jets
            
            AliVParticle* vp1 = nullptr;
			AliAODMCParticle* mcparticle = nullptr, *mcmother = nullptr;
            Int_t trackLabel = 0, pdgCode = 0, pdgMother = 0, motherIndex = -1;
            Int_t elecCounter = 0;
            Double_t ptTrueHFE = -1.;
            
            if(fMCarray){
                GetNumberOfTrueElectrons(jetbase, 0, nMotherKink, listofmotherkink, nTrueElectronsMC, nTrueHFElecMC, ptTrueHFE);
               
                for(UInt_t i = 0; i < jetbase -> GetNumberOfTracks(); i++) {
                    vp1 = jetbase -> Track(i);
                    
                    if (!vp1){
                        Printf("AliVParticle associated to constituent not found");
                        continue;
                    }

					trackLabel = vp1 -> GetLabel();
					if (trackLabel < 0) continue;
                    	
					mcparticle = static_cast<AliAODMCParticle*>(fMCarray->At(trackLabel)); 
                    if (!mcparticle) continue;

				    pdgCode = mcparticle -> PdgCode();
                    
                    if (TMath::Abs(pdgCode)==11) {
				    	elecCounter++;

				    	if (nTrueElectronsMC == 1) {
				    		motherIndex = mcparticle -> GetMother();
				    		mcmother = static_cast<AliAODMCParticle*>(fMCarray -> At(motherIndex));

				    		if (mcmother) pdgMother = mcmother -> PdgCode();
				    	}
				    }
				}

                if (elecCounter==1) fPtGenJet->Fill(jetbase->Pt());
				
            }
            
            fShapesVar[13] = ptTrueHFE;
            fShapesVar[14] = ptElec;
            fShapesVar[15] = nInclusiveElectrons;
            fShapesVar[16] = nPhotonicElectrons;
            fShapesVar[17] = hasElectrons;
            fShapesVar[18] = nTrueElectronsMC;
            fShapesVar[19] = nTrueHFElecMC;
            fShapesVar[20] = pdgMother;
            
			// Only fill tree if electron candidate or MC electron is found, ignore soft jets
			if ((nInclusiveElectrons > 0 || nTrueElectronsMC > 0) && ptSubtracted > 5.) {
            	fTreeObservableTagging->Fill();
			}
            
            
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
    
    Double_t bin_ept[nbins_ept + 1] = {0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4,
        1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 5.,
        6., 8., 10., 12., 14., 16., 19., 22., 25., 30.,
        35., 40., 45., 50.};
    
    Double_t bin_jetpt[nbins_jetpt + 1] = {5., 10., 20., 40., 60., 80., 120.};
    
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
                
                for (Int_t l = 0; l < nbins_jetpt; l++){// pt jet range
                    for(Int_t k = 0; k < nbins_eptTPC; k++){// pt electron range (up to 4 GeV/c)
                        if (jetPt>=bin_jetpt[l] && jetPt<bin_jetpt[l+1] && p>=bin_ept[k] && p<bin_ept[k+1]) fnTPCSigma[l][k]->Fill(fTPCnSigma);
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
            
            
            if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3) fnEovPelecTPCcut->Fill(pt,EovP);
            
			if (m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut) fnEovPelecEMCalcut -> Fill(pt, EovP);
            
            if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut){
                fnEovPelecTPCEMCalcut->Fill(pt,EovP);
                
                for (Int_t l = 0; l < nbins_jetpt; l++){// pt jet range
                    if (jetPt>=bin_jetpt[l] && jetPt<bin_jetpt[l+1]) fnEovPelecTPCsscut[l]->Fill(pt,EovP);
                }
            }
            
			// Not electrons
            if (TMath::Abs(fTPCnSigma) > fSigmTPCcutExcElec) { 
				fnM20backg -> Fill(pt, m20);
				fnM02backg -> Fill(pt, m02);
				fnEovPbackg -> Fill(pt, EovP);

				// Also consider not electrons with same shower shape cut as electron selection
				if (m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut) fnEovPbackgEMCalcut -> Fill(pt, EovP);
			}

            fnClsE->Fill(pt,clsE);
            fnM20->Fill(pt,m20);
            fnM02->Fill(pt,m02);
            fnClsTime->Fill(pt,clsTime);
            
            if (fTPCnSigma>fSigmaTPCcutHighPt  && fTPCnSigma<3 && EovP>fMinEoPcut && EovP<fMaxEoPcut && m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut && pt>=fMinPtEMCal && pt<fMaxPtEMCal){
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
	Int_t nIEptmin = 0, nHFEptmin = 0, nPEptmin = 0;
    Double_t p=-9., pt=-9., fTPCnSigma=-99., fTOFnSigma=-99., MCweight = 1., phi = -99., pte=0.;
    
    Double_t bin_jetpt[nbins_jetpt + 1] = {5., 10., 20., 40., 60., 80., 120.};
    Double_t bin_g[nbins_gsys + 1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
								   0.5, 0.55, 0.6, 0.65, 0.7};
    Double_t bin_ptd[nbins_ptdsys + 1] = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 
									   0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    
    Bool_t isFromHFdecay=kFALSE;
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
            
            p = track->P();
            pt = track->Pt();
            phi = track->Phi();
            
            fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            fTOFnSigma = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
            
            nPairs=0;
            MCweight = 1.;
            iDecay = 0;
            
            // Electron ID with EMCal
            
            Double_t EovP = -9., clsE = -9., m20 = -9., m02 = -9.;
            
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
                    m02 = cluster->GetM02();
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
                        
                        if (partPDG == 11){
                            nIE++;
                            pte=pt;
                            
                            if (isFromHFdecay) nHFE++;
                            if (iDecay>0 && iDecay<7) nPE++;

							if (pte > 0.5) {
								nIEptmin++;
                            	if (isFromHFdecay) nHFEptmin++;
                            	if (iDecay>0 && iDecay<7) nPEptmin++;
							}
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

						// EMCal shower shape studies
						if (m20 > 0.001 and m02 > 0.001) { // Ignore cluster with close to 0 axis lengths
							if (partPDG == 11) {
								fnM20TrueElec -> Fill(pt, m20);
								fnM02TrueElec -> Fill(pt, m02);
								fnShowerShapeTrueElec -> Fill(pt, m20, m02);
								fnEovPTrueElecnocut -> Fill(pt, EovP); 
                            	
								if (fTPCnSigma > fSigmaTPCcutHighPt && fTPCnSigma < 3) {
									fnEovPTrueElecTPCcut -> Fill(pt, EovP); 
								}
								
								if (m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut) {
									fnEovPTrueElecEMCalcut -> Fill(pt, EovP); 
								}

								if (fTPCnSigma > fSigmaTPCcutHighPt && fTPCnSigma < 3 && m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut) {
									fnEovPTrueElecTPCEMCalcut -> Fill(pt, EovP); 
								}
                        	} else {
								fnM20TrueBkg -> Fill(pt, m20);
								fnM02TrueBkg -> Fill(pt, m02);
								fnShowerShapeTrueBkg -> Fill(pt, m20, m02);
								fnEovPTrueBkgnocut -> Fill(pt, EovP); 
								
								if (fTPCnSigma > fSigmaTPCcutHighPt && fTPCnSigma < 3) {
									fnEovPTrueBkgTPCcut -> Fill(pt, EovP); 
								}
								
								if (m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut) {
									fnEovPTrueBkgEMCalcut -> Fill(pt, EovP); 
								}

								if (fTPCnSigma > fSigmaTPCcutHighPt && fTPCnSigma < 3 && m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut) {
									fnEovPTrueBkgTPCEMCalcut -> Fill(pt, EovP); 
								}
							}
						}
                        
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
                            
                            
                            if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut)
                                fptTrueHFEeffEMCal[2]->Fill(pt);
                            
                            if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && EovP>fMinEoPcut && EovP<fMaxEoPcut)
                                fptTrueHFEeffEMCal[3]->Fill(pt);
                            
                            if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && EovP>fMinEoPcut && EovP<fMaxEoPcut && m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut){
                                
                                fptTrueHFEeffEMCal[4]->Fill(pt);
                                fptTrueHFEeffEMCalang[1]->Fill(pt,jet->Pt(),GetJetAngularity(jet,0));
                                fptTrueHFEeffEMCaldisp[1]->Fill(pt,jet->Pt(),GetJetpTD(jet,0));
                            }
                        }
                        
                        if (fTPCnSigma>fSigmaTPCcutLowPt && fTPCnSigma<3 && TMath::Abs(fTOFnSigma)<fSigmaTOFcut && pt>=fMinPtTPC && pt<fMaxPtTPC) nRecIE++;
                        if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && EovP>fMinEoPcut && EovP<fMaxEoPcut &&
                            m20 > 0.01 && m20 < fM20cut && pt>=fMinPtEMCal && pt<fMaxPtEMCal) nRecIE++;
                        
                        // TPC-TOF
                        if (fTPCnSigma>fSigmaTPCcutLowPt && fTPCnSigma<3 && TMath::Abs(fTOFnSigma)<fSigmaTOFcut && pt>=fMinPtTPC && pt<fMaxPtTPC && nRecIE <2){
                            nPairs = GetNumberOfPairs(jet,track,pVtx,nMother,listMother,iDecay,MCweight);
                            if (nPairs>0 && iDecay>0 && iDecay<7) fptRecPE->Fill(pt,MCweight);
                            if (iDecay>0 && iDecay<7) fptTruePE->Fill(pt,MCweight);
                            
                            for (Int_t l = 0; l < nbins_jetpt; l++){// pt jet range
                                if (jetPt>=bin_jetpt[l] && jetPt<bin_jetpt[l+1]){
                                    
                                    for (Int_t m = 0; m < nbins_gsys; m++){// angularity
                                        if (ang>=bin_g[m] && ang<bin_g[m+1]){
                                            if (iDecay>0 && iDecay<7) fTotPEAng[l][m]->Fill(pt,MCweight);
                                            if (nPairs>0 && iDecay>0 && iDecay<7) fRecPEAng[l][m]->Fill(pt,MCweight);
                                        }
                                    }
                                    
                                    for (Int_t m = 0; m < nbins_ptdsys; m++){// dispersion
                                        if (disp>=bin_ptd[m] && disp<bin_ptd[m+1]){
                                            if (iDecay>0 && iDecay<7) fTotPEDisp[l][m]->Fill(pt,MCweight);
                                            if (nPairs>0 && iDecay>0 && iDecay<7) fRecPEDisp[l][m]->Fill(pt,MCweight);
                                        }
                                    }
                                }
                            }
                        }// PID cuts
                        
                        //EMCal
                        if (fTPCnSigma>fSigmaTPCcutHighPt && fTPCnSigma<3 && EovP>fMinEoPcut && EovP<fMaxEoPcut &&
                            m20 > 0.01 && m20 < fM20cut && m02 > fMinM02cut && m02 < fMaxM02cut && pt>=fMinPtEMCal && pt<fMaxPtEMCal && nRecIE <2){
                            nPairs = GetNumberOfPairs(jet,track,pVtx,nMother,listMother,iDecay,MCweight);
                            if (nPairs>0 && iDecay>0 && iDecay<7) fptRecPE->Fill(pt,MCweight);
                            if (iDecay>0 && iDecay<7) fptTruePE->Fill(pt,MCweight);
                            
                            for (Int_t l = 0; l < nbins_jetpt; l++){// pt jet range
                                if (jetPt>=bin_jetpt[l] && jetPt<bin_jetpt[l+1]){
                                    
                                    for (Int_t m = 0; m < nbins_gsys; m++){// angularity
                                        if (ang>=bin_g[m] && ang<bin_g[m+1]){
                                            if (iDecay>0 && iDecay<7) fTotPEAng[l][m]->Fill(pt,MCweight);
                                            if (nPairs>0 && iDecay>0 && iDecay<7) fRecPEAng[l][m]->Fill(pt,MCweight);
                                        }
                                    }
                                    
                                    for (Int_t m = 0; m < nbins_ptdsys; m++){// dispersion
                                        if (disp>=bin_ptd[m] && disp<bin_ptd[m+1]){
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
    
	if (nIEptmin > 0) fnTrueElecPerJetPt -> Fill(jet -> Pt(), nIEptmin);
    if (nHFEptmin > 0) fnTrueHFElecPerJetPt -> Fill(jet -> Pt(), nHFEptmin);
    if (nPEptmin > 0) fnTruePElecPerJetPt -> Fill(jet -> Pt(), nPEptmin);
    
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
    
    Double_t bin_jetpt[nbins_jetpt + 1] = {5., 10., 20., 40., 60., 80., 120.};
    Double_t bin_g[nbins_gsys + 1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
								   0.5, 0.55, 0.6, 0.65, 0.7};
    Double_t bin_ptd[nbins_ptdsys + 1] = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 
									   0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    
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
        Double_t pt = -9;
        Int_t charge = 0;
        
        pt = track->Pt();
        charge = track->Charge();
        
        
        // associated particle variables
        Double_t fTPCnSigmaAsso=-9.;
        Int_t chargeAsso = 0;
        
        chargeAsso = trackAsso->Charge();
        
        // looser PID cuts
        fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(trackAsso, AliPID::kElectron);
        
        if(TMath::Abs(fTPCnSigmaAsso)>3.) continue;
        
        // invariant mass
        Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
		Double_t mass=999., width = -999;
        
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        AliKFParticle::SetField(fVevent->GetMagneticField());
        AliKFParticle ge1(*track, fPDGe1);
        AliKFParticle ge2(*trackAsso, fPDGe2);
        AliKFParticle recg(ge1, ge2);

		Int_t MassCorrect = recg.GetMass(mass, width);

        if(recg.GetNDF()<1) continue;

        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
		if (TMath::Abs(chi2recg) > 3.) continue;
        
        for (Int_t l = 0; l < nbins_jetpt; l++){// pt jet range
            if (jetPt>=bin_jetpt[l] && jetPt<bin_jetpt[l+1] && fFlagULS) fInvmassULS[l]->Fill(pt,mass);
            if (jetPt>=bin_jetpt[l] && jetPt<bin_jetpt[l+1] && fFlagLS)  fInvmassLS[l]->Fill(pt,mass);
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
        
        
        for (Int_t l = 0; l < nbins_jetpt; l++){// pt jet range
            if (jetPt>=bin_jetpt[l] && jetPt<bin_jetpt[l+1]){
                
                for (Int_t m = 0; m < nbins_gsys; m++){// angularity
                    if (ang>=bin_g[m] && ang<bin_g[m+1] && decay>0 && decay<7 && mass<fIMcut){
                        if (fFlagULS) fULSptAng[l][m]->Fill(pt,weight);
                        if (fFlagLS) fLSptAng[l][m]->Fill(pt,weight);
                    }
                }
                
                for (Int_t m = 0; m < nbins_ptdsys; m++){// dispersion
                    if (disp>=bin_ptd[m] && disp<bin_ptd[m+1]  && decay>0 && decay<7 && mass<fIMcut){
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
	Double_t jetradius = jetCont -> GetJetRadius();

    if (!jet->GetNumberOfTracks()) return 0;

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

    return num / (den * jetradius);
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb = 0){
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
	Double_t jetradius = jetCont -> GetJetRadius();
    
    if((fJetShapeSub==kDerivSub) && (jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity() / jetradius;
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity() / jetradius;
        else
            return Angularity(jet, jetContNb);
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::PTD(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks()) return 0;

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

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb=0){
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedConstituent();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedConstituent();
        else
            return jet->GetNumberOfTracks();
    
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalHfeTagging::MaxPtBinForSemiInclusiveJet(AliEmcalJet *jet, Int_t jetContNb=0){
    // Give the bin with max value for minimum pT bin for semi-inclusive observables due to the 
	// min(pT) = [0.5, 1.0, 2.5, user defined] GeV
	// No need for transverse mass rescaling since lowest pT cut >> mass of electron
	
	// This array must be ordered
	double minpT[5] = {0, 0.5, 1.0, 2.5, fMinPtSemiInclusive};

    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks()) return -1;

	// Get largest constituent pt
	double largestpt = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        if (vp1 -> Pt() > largestpt) largestpt = vp1 -> Pt();
    }

	// If there is a particle with pt > min pt[j], the jet should populate all the bins previous
	// to j (inclusive) 
	for (int j = 4; j > 0; j--) {
    	if (largestpt > minpT[j]) {
        	return j;
        }
	}
	
	return 0;
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

//_________________________________________
Double_t AliAnalysisTaskEmcalHfeTagging::AngularDifference(AliVParticle* jet1, AliVParticle* jet2) {
  Double_t phi1 = jet1 -> Phi(), eta1 = jet1 -> Eta();
  Double_t phi2 = jet2 -> Phi(), eta2 = jet2 -> Eta();

  Double_t deltaeta = eta1 - eta2, deltaphi = RelativePhi(phi1, phi2);

  if (deltaphi > TMath::Pi()) deltaphi -= 2 * TMath::Pi();
  if (deltaphi < -TMath::Pi()) deltaphi += 2 * TMath::Pi();

  return TMath::Sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
}

//_________________________________________
AliEmcalJet* AliAnalysisTaskEmcalHfeTagging::GetClosestOnOtherJetContainer(AliEmcalJet* jet1, AliJetContainer* othercontainer) {
    AliEmcalJet* closestjet = nullptr;
	AliEmcalJet* jet2 = nullptr;
	// If distance is always higher than pi, return null;
	Double_t distmin = TMath::Pi();

	for (Int_t i = 0; i < othercontainer -> GetNJets(); i++) {
		jet2 = othercontainer -> GetJet(i);
		Double_t dist = AngularDifference(jet1, jet2);
		/* cout << "dist: " << dist << endl; */

		if (dist < distmin) {
			distmin = dist;
			closestjet = jet2;
		}
	}

	return closestjet;
}

//_________________________________________
Double_t AliAnalysisTaskEmcalHfeTagging::GetFractionSharedPtBetweenJets(AliEmcalJet* jet1, AliEmcalJet* jetmatched) {
	// Similar to AliJetContainer::GetFractionSharedPt() but using any two jets instead of closest and geometric matching
	// The angular distance threshold is set as 0.02 and 5% for pT matching
	
	if (!jet1 || !jetmatched) return -1;
	
	Double_t jetmatchedpt = jetmatched -> Pt(), sumpt = 0.;
	if (jetmatchedpt < 0. || jet1 -> Pt() < 0) return -1;

	AliVParticle* p1 = nullptr, *p2 = nullptr;
	
	for (Int_t ic2 = 0; ic2 < jetmatched -> GetNumberOfTracks(); ic2++) {
		p2 = jetmatched -> Track(ic2);
		if (!p2) continue;
		
		for (Int_t ic1 = 0; ic1 < jet1 -> GetNumberOfTracks(); ic1++) {
			p1 = jet1 -> Track(ic1);
			if (!p1) continue;

			if (AngularDifference(p1, p2) < 0.02 && abs((p1 -> Pt() - p2 -> Pt()) / p2 -> Pt()) < 0.05) {
				sumpt += p2 -> Pt();
				break;
			}
		}
	}

	return sumpt / jetmatchedpt;
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

