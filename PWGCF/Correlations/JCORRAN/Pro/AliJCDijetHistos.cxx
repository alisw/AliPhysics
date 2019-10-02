/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

// Container class for histograms needed in the analysis.

#include "AliJCDijetHistos.h"
#include <TGrid.h>
#include <TPRegexp.h>

//Double_t AliJCDijetHistos::pttJacek[74+16] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500};
//UInt_t AliJCDijetHistos::NpttJacek = sizeof(AliJCDijetHistos::pttJacek)/sizeof(AliJCDijetHistos::pttJacek[0])-1;
//int const nALICEBins = 8;
//double ptBinsALICE[nALICEBins+1] = { 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0 };
vector<double> AliJCDijetHistos::CentBin;
int AliJCDijetHistos::fNCentBin;

//______________________________________________________________________________
AliJCDijetHistos::AliJCDijetHistos() :
    fHMG(NULL),
    fHistCentBin(),
    fJetBin(),
    fh_events(),
    fh_centrality(),
    fh_zvtx(),
    fh_pt(),
    fh_eta(),
    fh_phi(),
    fh_rho(),
    fh_rhom(),
    fh_jetPt(),
    fh_jetPt_ALICE(),
    fh_jetEta(),
    fh_jetPhi(),
    fh_jetEtaPhi(),
    fh_jetArea(),
    fh_jetAreaRho(),
    fh_dijetInvM(),
    fh_dijetPtPair(),
    fh_dijetDeltaPhi(),
    fh_dijetPtPairDeltaPhiCut(),
    fh_dijetInvMDeltaPhiCut(),
    fh_dijetDeltaPhiWithCut(),
    fh_responseInfo(),
    fh_jetResponseDeltaR(),
    fh_jetResponseDeltaRClosest(),
    fh_jetResponseDeltaPt(),
    fh_jetDeltaRMin(),
    fh_jetBGSubtrDeltaR(),
    fh_jetResponse(),
    fh_dijetResponse(),
    fh_dijetResponseDeltaPhiCut()
{

}

//______________________________________________________________________________
AliJCDijetHistos::AliJCDijetHistos(const AliJCDijetHistos& obj) :
    fHMG(obj.fHMG),
    fHistCentBin(obj.fHistCentBin),
    fJetBin(obj.fJetBin),
    fh_events(obj.fh_events),
    fh_centrality(obj.fh_centrality),
    fh_zvtx(obj.fh_zvtx),
    fh_pt(obj.fh_pt),
    fh_eta(obj.fh_eta),
    fh_phi(obj.fh_phi),
    fh_rho(obj.fh_rho),
    fh_rhom(obj.fh_rhom),
    fh_jetPt(obj.fh_jetPt),
    fh_jetPt_ALICE(obj.fh_jetPt_ALICE),
    fh_jetEta(obj.fh_jetEta),
    fh_jetPhi(obj.fh_jetPhi),
    fh_jetEtaPhi(obj.fh_jetEtaPhi),
    fh_jetArea(obj.fh_jetArea),
    fh_jetAreaRho(obj.fh_jetAreaRho),
    fh_dijetInvM(obj.fh_dijetInvM),
    fh_dijetPtPair(obj.fh_dijetPtPair),
    fh_dijetDeltaPhi(obj.fh_dijetDeltaPhi),
    fh_dijetPtPairDeltaPhiCut(obj.fh_dijetPtPairDeltaPhiCut),
    fh_dijetInvMDeltaPhiCut(obj.fh_dijetInvMDeltaPhiCut),
    fh_dijetDeltaPhiWithCut(obj.fh_dijetDeltaPhiWithCut),
    fh_responseInfo(obj.fh_responseInfo),
    fh_jetResponseDeltaR(obj.fh_jetResponseDeltaR),
    fh_jetResponseDeltaRClosest(obj.fh_jetResponseDeltaRClosest),
    fh_jetResponseDeltaPt(obj.fh_jetResponseDeltaPt),
    fh_jetDeltaRMin(obj.fh_jetDeltaRMin),
    fh_jetBGSubtrDeltaR(obj.fh_jetBGSubtrDeltaR),
    fh_jetResponse(obj.fh_jetResponse),
    fh_jetResponse_ALICE(obj.fh_jetResponse_ALICE),
    fh_dijetResponse(obj.fh_dijetResponse),
    fh_dijetResponseDeltaPhiCut(obj.fh_dijetResponseDeltaPhiCut)
{
    // copy constructor
}

//______________________________________________________________________________
AliJCDijetHistos& AliJCDijetHistos::operator=(const AliJCDijetHistos& obj){
    // copy constructor
    return *this;
}

//______________________________________________________________________________
AliJCDijetHistos::~AliJCDijetHistos() {
    // destructor
    delete fHMG;
}


//______________________________________________________________________________
void AliJCDijetHistos::CreateEventTrackHistos(){
    // Create basic event histograms
    fHMG = new AliJHistManager(Form("AliJCDijetHistManager%s",sMngrName.Data()),sMngrName.Data());
    // set AliJBin here //
    fHistCentBin.Set("CentBin","CentBin","Cent:",AliJBin::kSingle).SetBin(fNCentBin);
    fJetBin.Set("JetBin","JetBin","Jet bin:",AliJBin::kSingle).SetBin(5);

    // fh_events counts several things:
    // 0:  Number of events
    // 1:  Number of ch. particles
    // 2:  Number of accepted ch. particles
    // 3:  Number of events with no rho calculations
    // 4:  Number of events with proper rho calculations
    // 5:  Number of jets
    // 6:  Number of accepted jets
    // 7:  Number of accepted jets after const. cut
    // 8:  Number of accepted bg subtracted jets
    // 9:  Number of accepted bg subtracted const. cut jets
    // 10: Number of kt-jets
    // 11: Number of accepted kt-jets
    // 12: Number of jets that drop under leading pt cut after bg subtraction
    // 13: Number of jets that drop under subleading pt cut after bg subtraction
    // 14: Number of raw dijets
    // 15: Number of raw dijets after leading pt cut
    // 16: Number of accepted raw dijets
    // 17: Number of accepted raw dijets with delta phi cut
    // 14: Number of bg subtr. dijets
    // 15: Number of bg subtr. dijets after leading pt cut
    // 16: Number of accepted bg subtr. dijets
    // 17: Number of accepted bg subtr. dijets with delta phi cut
    // 18: Number of bg subtr. const. cut dijets
    // 19: Number of bg subtr. const. cut dijets after leading pt cut
    // 20: Number of accepted bg subtr. const. cut dijets
    // 21: Number of accepted bg subtr. const. cut dijets with delta phi cut
    // 22: Number of const. cut dijets
    // 23: Number of const. cut dijets after leading pt cut
    // 24: Number of accepted const. cut dijets
    // 25: Number of accepted const. cut dijets with delta phi cut
    // 26: Number of kt-dijets
    // 27: Number of kt-dijets after leading pt cut
    // 28: Number of accepted kt-dijets
    // 29: Number of accepted kt-dijets with delta phi cut
    fh_events
        << TH1D("h_events", "h_events", 40, 0.0, 40.0 )
        << fHistCentBin
        << "END" ;

    fh_info
        << TH1D("h_info", "h_info", 40, 0.0, 40.0 )
        << "END" ;

    fh_centrality
        << TH1D("h_centrality", "h_centrality", 100, 0.0, 100.0 )
        << "END" ;

    fh_zvtx
        << TH1D("h_zvtx", "h_zvtx", 40, -20.0, 20.0 )
        << "END" ;

    int NBINSJet=150;
    double LogBinsXJet[NBINSJet+1], LimLJet=0.1, LimHJet=500;
    double logBWJet = (log(LimHJet)-log(LimLJet))/NBINSJet;
    for(int ijetBin=0;ijetBin<=NBINSJet;ijetBin++) LogBinsXJet[ijetBin]=LimLJet*exp(ijetBin*logBWJet);

    // ============= CHARGED PARTICLE HISTOS ============= 
    fh_pt
        //<< TH1D("h_pt", "h_pt", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
        << TH1D("h_pt","h_pt",NBINSJet, LogBinsXJet )
        << fHistCentBin
        << "END" ;

    fh_eta
        << TH1D("h_eta", "h_eta", 100, -1.0, 1.0 )
        << fHistCentBin
        << "END" ;

    fh_phi
        << TH1D("h_phi", "h_phi", 100, -TMath::Pi(), TMath::Pi())
        << fHistCentBin
        << "END" ;

    fh_etaPhi
        << TH2D("h_etaPhi", "h_etaPhi", 100, -1.0, 1.0, 100, -TMath::Pi(), TMath::Pi())
        << fHistCentBin
        << "END" ;

    // ============= JET HISTOS ============= 
    fh_rho
        << TH1D("h_rho", "h_rho", NBINSJet, LogBinsXJet)
        << fHistCentBin
        << "END" ;

    fh_rhom
        << TH1D("h_rhom", "h_rhom", NBINSJet, LogBinsXJet)
        << fHistCentBin
        << "END" ;

    fh_jetPt
        //<< TH1D("h_jetPt", "h_jetPt", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
        << TH1D("h_jetPt","h_jetPt",NBINSJet, LogBinsXJet )
        << fHistCentBin << fJetBin
        << "END" ;

    fh_jetPt_ALICE
        << TH1D("h_jetPt_ALICE","h_jetPt_ALICE", 200, 10, 210 )
        << fHistCentBin << fJetBin
        << "END" ;

    fh_jetEta
        << TH1D("h_jetEta", "h_jetEta", 100, -1.0, 1.0)
        << fHistCentBin << fJetBin
        << "END" ;

    fh_jetPhi
        << TH1D("h_jetPhi", "h_jetPhi", 100, -TMath::Pi(), TMath::Pi())
        << fHistCentBin << fJetBin
        << "END" ;

    fh_jetEtaPhi
        << TH2D("h_jetEtaPhi", "h_jetEtaPhi", 100, -1.0, 1.0, 100, -TMath::Pi(), TMath::Pi())
        << fHistCentBin << fJetBin
        << "END" ;

    fh_jetArea
        //<< TH1D("h_jetArea", "h_jetArea", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
        << TH1D("h_jetArea", "h_jetArea", NBINSJet, LogBinsXJet )
        << fHistCentBin << fJetBin
        << "END" ;

    fh_jetAreaRho
        //<< TH1D("h_jetAreaRho", "h_jetAreaRho", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek)
        << TH1D("h_jetAreaRho", "h_jetAreaRho", NBINSJet, LogBinsXJet )
        << fHistCentBin << fJetBin
        << "END" ;

    int NBINSDijet=170;
    double logBinsXDijet[NBINSDijet+1], LimLDijet=0.1, LimHDijet=1000;
    double logBWDijet = (log(LimHDijet)-log(LimLDijet))/NBINSDijet;
    for(int iDijet=0;iDijet<=NBINSDijet;iDijet++) logBinsXDijet[iDijet]=LimLDijet*exp(iDijet*logBWDijet);

    // ============= DIJET HISTOS ============= 
    fh_dijetInvM
        << TH1D("h_dijetInvM", "h_dijetInvM", NBINSDijet, logBinsXDijet)
        << fHistCentBin << fJetBin
        << "END" ;

    fh_dijetPtPair
        //<< TH1D("h_dijetPtPair", "h_dijetPtPair", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek )
        << TH1D("h_dijetPtPair", "h_dijetPtPair", NBINSDijet, logBinsXDijet )
        << fHistCentBin << fJetBin
        << "END" ;

    fh_dijetDeltaPhi
        << TH1D("h_dijetDeltaPhi", "h_dijetDeltaPhi", 100, 0, 10)
        << fHistCentBin << fJetBin
        << "END" ;

    fh_dijetPtPairDeltaPhiCut
        //<< TH1D("h_dijetPtPair", "h_dijetPtPair", AliJCDijetHistos::NpttJacek, AliJCDijetHistos::pttJacek )
        << TH1D("h_dijetPtPairDeltaPhiCut", "h_dijetPtPairDeltaPhiCut", NBINSDijet, logBinsXDijet )
        << fHistCentBin << fJetBin
        << "END" ;

    fh_dijetInvMDeltaPhiCut
        << TH1D("h_dijetInvMDeltaPhiCut", "h_dijetInvMDeltaPhiCut", NBINSDijet, logBinsXDijet)
        << fHistCentBin << fJetBin
        << "END" ;

    fh_dijetDeltaPhiWithCut
        << TH1D("h_dijetDeltaPhiWithCut", "h_dijetDeltaPhiWithCut", 100, 0, 10)
        << fHistCentBin << fJetBin
        << "END" ;


    // ============ Response histograms ===========
    fh_responseInfo
        << TH1D("h_responseInfo", "h_responseInfo", 40, 0.0, 40.0 )
        << "END" ;

    fh_jetResponseDeltaR
        << TH1D("h_jetResponseDeltaR", "h_jetResponseDeltaR", 100, 0.0, 1.0)
        << "END" ;

    fh_jetResponseDeltaRClosest
        << TH1D("h_jetResponseDeltaRClosest", "h_jetResponseDeltaRClosest", 100, 0.0, 1.0)
        << "END" ;

    fh_jetResponseDeltaPt
        << TH1D("h_jetResponseDeltaPt", "h_jetResponseDeltaPt", 200, -2.0, 1.0)
        << "END" ;

    fh_jetDeltaRMin
        << TH1D("h_jetDeltaRMin", "h_jetDeltaRMin", 400, 0.0, 4.0)
        << fHistCentBin << fJetBin
        << "END" ;

    fh_jetBGSubtrDeltaR
        << TH1D("h_jetBGSubtrDeltaR", "h_jetBGSubtrDeltaR", 400, 0.0, 4.0)
        << fHistCentBin
        << "END" ;

    fh_jetResponse
        << TH2D("h_jetResponse", "h_jetResponse", NBINSJet, LogBinsXJet, NBINSJet, LogBinsXJet )
        << "END" ;

    fh_jetResponse_ALICE
        << TH2D("h_jetResponse_ALICE", "h_jetResponse_ALICE", 200, 10, 210, 200, 10, 210 )
        << "END" ;

    fh_dijetResponse
        << TH2D("h_dijetResponse", "h_dijetResponse", NBINSDijet, logBinsXDijet, NBINSDijet, logBinsXDijet )
        << "END" ;

    fh_dijetResponseDeltaPhiCut
        << TH2D("h_dijetResponseDeltaPhiCut", "h_dijetResponseDeltaPhiCut", NBINSDijet, logBinsXDijet, NBINSDijet, logBinsXDijet )
        << "END" ;
}

int AliJCDijetHistos::GetCentralityClass(Double_t fCent){
    for(int iCbin = 0; iCbin < fNCentBin; iCbin++){
        if(fCent > CentBin[iCbin] && fCent < CentBin[iCbin+1])
            return iCbin;
    }
    return -1;
}

