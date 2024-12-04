//Task for ENC in pbpb collisions
//Code inherited and modified from AliAnalysisTaskLundPlane.cxx authored by Leticia Cunqueiro and Laura Havener
//Authors: Ananya Rai, Laura Havener
#include "AliAODMCHeader.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliEmcalPythiaInfo.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetContainer.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliMCParticleContainer.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliTLorentzVector.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include <AliAnalysisDataContainer.h>
#include <AliAnalysisDataSlot.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TMath.h>
#include "AliAODEvent.h"
#include "AliEmcalContainerUtils.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliAnalysisTaskJetsEECpbpb.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetsEECpbpb)

//________________________________________________________________________
AliAnalysisTaskJetsEECpbpb::AliAnalysisTaskJetsEECpbpb(): AliAnalysisTaskEmcalJet("AliAnalysisTaskJetsEECpbpb", kTRUE),
fContainer(0), fMinFractionShared(0), fJetShapeType(kData),
fJetShapeSub(kNoSub), fJetSelection(kInclusive), fPtThreshold(-9999.), fMinENCtrackPt(1.0), fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE), fCheckResolution(kFALSE),
fMinPtConst(1), fHardCutoff(0), fDoTwoTrack(kFALSE), fCutDoubleCounts(kTRUE),
fPowerAlgo(1), fPhiCutValue(0.02),
fEtaCutValue(0.9), fDerivSubtrOrder(0),
fStoreDetLevelJets(0), fDoFillEncMC(kFALSE), fMatchR(0.2), fjetMinPtSub(20), fjetMinArea(0), fStoreTrig(kFALSE), fpTcorr(0), fpaircut(0), 
fpairfastsim(0), fUnfolding(0),fMatchJetTrack(1),fMissJetTrack(1),fFakeJetTrack(1), fMaxPtTrack(0), fJetPtMin(0),
fDoEmbedding(kFALSE), fIsEmbeddedEvent(kFALSE),fDoPerpCone(kFALSE),fDoRandCone(kFALSE), fCout(kTRUE), fGeneratorLevelName(), fDetectorLevelName(), fGeneratorLevel(0), fDetectorLevel(0), 
fJet_truCont(0),jet_pt_hist(0), EEC_hist(0), EEC_pt_hist(0), E3C_hist(0), E3C_pt_hist(0),  EEC_det_pt_hist_3d(0), 
EEC_tru_pt_hist_3d(0), E3C_det_pt_hist_3d(0), E3C_tru_pt_hist_3d(0), N2_det_pt_hist_3d(0), N2_tru_pt_hist_3d(0), 
N3_det_pt_hist_3d(0), N3_tru_pt_hist_3d(0), EEC_det_match_pt_det(0), EEC_tru_match_pt_tru(0), 
E3C_det_match_pt_det(0), E3C_tru_match_pt_tru(0), pt_tru(0), pt_tru_match(0), pt_det(0), 
pt_det_match(0), test_hist(0), R_matrix(0), JES(0), JES_scaled(0), JER(0), pair_det_EEC(0), 
pair_tru_EEC(0), pair_det_E3C(0), pair_tru_E3C(0), qpt_tru(0), qpt_det(0), track_pt_tru(0), 
track_pt_det(0), track_pt_matched(0), R_match_eec(0), wt_match_eec(0), R_match_e3c(0), 
wt_match_e3c(0), qpt_tru1(0), qpt_tru2(0),pt_tru1(0), pt_tru2(0), eec_Mm(0), eec_mm(0), 
e3c_MMm(0), e3c_Mmm(0), e3c_mmm(0), eec_Mf(0), eec_ff(0), e3c_MMf(0), e3c_Mff(0), e3c_fff(0), 
eec_matched_det(0), eec_matched_tru(0), e3c_matched_det(0), e3c_matched_tru(0), wt_res_eec(0), 
wt_res_e3c(0), R_res_eec(0), R_res_e3c(0), wtnojet_match_eec(0), wtnojet_match_e3c(0), 
wtnojet_res_eec(0), wtnojet_res_e3c(0), R_match_eec_tru(0), wt_match_eec_tru(0), R_match_e3c_tru(0), 
wt_match_e3c_tru(0), wt_res_eec_tru(0), wt_res_e3c_tru(0), R_res_eec_tru(0), R_res_e3c_tru(0), 
wtnojet_match_eec_tru(0), wtnojet_match_e3c_tru(0), wtnojet_res_eec_tru(0), wtnojet_res_e3c_tru(0), 
constituentId(0),R_res_eec_tru_debug(0),track_pt_res_debug(0),track_eta_debug(0), track_rap_debug(0), 
track_phi_debug(0), track_R_debug(0), track_R_debug_rap(0), track_eta_res_debug(0), track_rap_res_debug(0), 
track_phi_res_debug(0), track_R_res_debug(0), track_R_rap_res_debug(0), R_match_eec_tru_rap(0), 
track_pt_response_debug(0), track_pt_wt_res_debug(0), track_pt_wt_response_debug(0), jetpt_res_w_R(0), jetpt_res_w_wt(0), 
h_MJ(0),h_MJ0(0),h_MJ1(0),h_MJ2(0),h_MJ_tru(0),h_MJ0_tru(0),h_MJ1_tru(0),h_MJ2_tru(0), h_MB1(0), h_JMB(0), h_SMB(0), h_BMB(0),
h_MB1MB2(0), h_MB1_tru(0), h_JMB_tru(0), h_SMB_tru(0), h_BMB_tru(0), h_MB1MB2_tru(0),h_MJ_e3c(0),h_MJ0_e3c(0),h_MJ1_e3c(0),h_MJ2_e3c(0),h_MJ3_e3c(0),
h_MJ_e3c_tru(0),h_MJ0_e3c_tru(0),h_MJ1_e3c_tru(0),h_MJ2_e3c_tru(0),h_MJ3_e3c_tru(0),h_MB1MB1MB1(0),h_MB1MB1MB1_tru(0),
h_JJMB(0), h_JJMB_tru(0), h_JMBMB(0), h_JMBMB_tru(0), h_MB1MB1MB2(0), h_MB1MB1MB2_tru(0), h_MB1MB2MB2(0),
h_MB1MB2MB2_tru(0), h_JMB1MB2(0), h_JMB1MB2_tru(0), h_MB1MB2MB3(0), h_MB1MB2MB3_tru(0),h_BMBMB(0), h_SMBMB(0),
h_BBMB(0), h_BMB1MB2(0), h_SMB1MB2(0), h_SBMB(0), h_SSMB(0), h_BBMB_tru(0), h_SBMB_tru(0), h_SSMB_tru(0),
h_BMBMB_tru(0), h_SMBMB_tru(0), h_BMB1MB2_tru(0), h_SMB1MB2_tru(0), fTreeMatchTracks(0), fTreeData(0), fMCParticleArrayName("mcparticles"), fMCParticleArray(0),
fRandom(0x0) 
{
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}


//________________________________________________________________________
AliAnalysisTaskJetsEECpbpb::AliAnalysisTaskJetsEECpbpb(const char *name): AliAnalysisTaskEmcalJet(name, kTRUE),
fContainer(0), fMinFractionShared(0), fJetShapeType(kData),
fJetShapeSub(kNoSub), fJetSelection(kInclusive), fPtThreshold(-9999.), fMinENCtrackPt(1.0), fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE), fCheckResolution(kFALSE),
fMinPtConst(1), fHardCutoff(0), fDoTwoTrack(kFALSE), fCutDoubleCounts(kTRUE),
fPowerAlgo(1), fPhiCutValue(0.02),
fEtaCutValue(0.02), fDerivSubtrOrder(0),
fStoreDetLevelJets(0), fDoFillEncMC(kFALSE), fMatchR(0.2), fjetMinPtSub(20), fjetMinArea(0), fStoreTrig(kFALSE), fpTcorr(0), fpaircut(0), 
fpairfastsim(0), fUnfolding(0),fMatchJetTrack(1),fMissJetTrack(1),fFakeJetTrack(1), fMaxPtTrack(0), fJetPtMin(0),
fDoEmbedding(kFALSE), fIsEmbeddedEvent(kFALSE),fDoPerpCone(kFALSE),fDoRandCone(kFALSE), fCout(kTRUE), fGeneratorLevelName(), fDetectorLevelName(), fGeneratorLevel(0), fDetectorLevel(0), 
fJet_truCont(0),jet_pt_hist(0), EEC_hist(0), EEC_pt_hist(0), E3C_hist(0), E3C_pt_hist(0),  EEC_det_pt_hist_3d(0), 
EEC_tru_pt_hist_3d(0), E3C_det_pt_hist_3d(0), E3C_tru_pt_hist_3d(0), N2_det_pt_hist_3d(0), N2_tru_pt_hist_3d(0), 
N3_det_pt_hist_3d(0), N3_tru_pt_hist_3d(0), EEC_det_match_pt_det(0), EEC_tru_match_pt_tru(0), 
E3C_det_match_pt_det(0), E3C_tru_match_pt_tru(0), pt_tru(0), pt_tru_match(0), pt_det(0), 
pt_det_match(0), test_hist(0), R_matrix(0), JES(0), JES_scaled(0), JER(0), pair_det_EEC(0), 
pair_tru_EEC(0), pair_det_E3C(0), pair_tru_E3C(0), qpt_tru(0), qpt_det(0), track_pt_tru(0), 
track_pt_det(0), track_pt_matched(0), R_match_eec(0), wt_match_eec(0), R_match_e3c(0), 
wt_match_e3c(0), qpt_tru1(0), qpt_tru2(0),pt_tru1(0), pt_tru2(0), eec_Mm(0), eec_mm(0), 
e3c_MMm(0), e3c_Mmm(0), e3c_mmm(0), eec_Mf(0), eec_ff(0), e3c_MMf(0), e3c_Mff(0), e3c_fff(0), 
eec_matched_det(0), eec_matched_tru(0), e3c_matched_det(0), e3c_matched_tru(0), wt_res_eec(0), 
wt_res_e3c(0), R_res_eec(0), R_res_e3c(0), wtnojet_match_eec(0), wtnojet_match_e3c(0), 
wtnojet_res_eec(0), wtnojet_res_e3c(0), R_match_eec_tru(0), wt_match_eec_tru(0), R_match_e3c_tru(0), 
wt_match_e3c_tru(0), wt_res_eec_tru(0), wt_res_e3c_tru(0), R_res_eec_tru(0), R_res_e3c_tru(0), 
wtnojet_match_eec_tru(0), wtnojet_match_e3c_tru(0), wtnojet_res_eec_tru(0), wtnojet_res_e3c_tru(0), 
constituentId(0),R_res_eec_tru_debug(0),track_pt_res_debug(0),track_eta_debug(0), track_rap_debug(0), 
track_phi_debug(0), track_R_debug(0), track_R_debug_rap(0), track_eta_res_debug(0), track_rap_res_debug(0), 
track_phi_res_debug(0), track_R_res_debug(0), track_R_rap_res_debug(0), R_match_eec_tru_rap(0), 
track_pt_response_debug(0), track_pt_wt_res_debug(0), track_pt_wt_response_debug(0), jetpt_res_w_R(0), jetpt_res_w_wt(0), 
h_MJ(0),h_MJ0(0),h_MJ1(0),h_MJ2(0),h_MJ_tru(0),h_MJ0_tru(0),h_MJ1_tru(0),h_MJ2_tru(0), h_MB1(0), h_JMB(0), h_SMB(0), h_BMB(0),
h_MB1MB2(0), h_MB1_tru(0), h_JMB_tru(0), h_SMB_tru(0), h_BMB_tru(0), h_MB1MB2_tru(0),h_MJ_e3c(0),h_MJ0_e3c(0),h_MJ1_e3c(0),h_MJ2_e3c(0),h_MJ3_e3c(0),
h_MJ_e3c_tru(0),h_MJ0_e3c_tru(0),h_MJ1_e3c_tru(0),h_MJ2_e3c_tru(0),h_MJ3_e3c_tru(0),h_MB1MB1MB1(0),h_MB1MB1MB1_tru(0),
h_JJMB(0), h_JJMB_tru(0), h_JMBMB(0), h_JMBMB_tru(0), h_MB1MB1MB2(0), h_MB1MB1MB2_tru(0), h_MB1MB2MB2(0),
h_MB1MB2MB2_tru(0), h_JMB1MB2(0), h_JMB1MB2_tru(0), h_MB1MB2MB3(0), h_MB1MB2MB3_tru(0),h_BMBMB(0), h_SMBMB(0),
h_BBMB(0), h_BMB1MB2(0), h_SMB1MB2(0), h_SBMB(0), h_SSMB(0), h_BBMB_tru(0), h_SBMB_tru(0), h_SSMB_tru(0),
h_BMBMB_tru(0), h_SMBMB_tru(0), h_BMB1MB2_tru(0), h_SMB1MB2_tru(0), fTreeMatchTracks(0), fTreeData(0), fMCParticleArrayName("mcparticles"), fMCParticleArray(0),
fRandom(0x0)
{
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}


//________________________________________________________________________
AliAnalysisTaskJetsEECpbpb::~AliAnalysisTaskJetsEECpbpb() {
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::UserCreateOutputObjects() {

if(fCout)  std::cout << " Info::anrai: ===== In the UserCreateOutputObjects ===== " << std::endl;

// Create user output.
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    TH1::AddDirectory(oldStatus);


  TRandom3* fRandom = new TRandom3(0);

  for(Int_t iCont=0; iCont<fParticleCollArray.GetEntriesFast(); iCont++){
    if (GetParticleContainer(iCont)->GetIsEmbedding()){
      fIsEmbeddedEvent = kTRUE;
    }
  }

  if (fDoEmbedding)
  {
    if(!fIsEmbeddedEvent){
      fMCParticleArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCParticleArrayName.Data()));
    }
    else
    {
      // In case of embedding, the MC particle array needs to be fetched differently
      AliVEvent* event = AliEmcalContainerUtils::GetEvent(InputEvent(), kTRUE);
      fMCParticleArray = dynamic_cast<TClonesArray*>(event->FindListObject(fMCParticleArrayName.Data()));
    }
  }

   if(fCout) cout<<"histos being initialized"<<endl;
    //RL bins
    Double_t from = -4;
    Double_t to = 0;
    Int_t bins = 100;
    Double_t width = (to-from)/bins;
    Double_t new_bins[101] = {};
    for (Int_t i = 0; i <= bins; i++)
    {
        new_bins[i] = TMath::Power(10.0, from + i * width);
    }

    // jet pT bins
    Double_t from_const = 10;
    Double_t to_const = 120;
    Int_t bins_const = 22;
    Double_t width_const = (to_const-from_const)/bins_const;
    Double_t new_bins_const[23] = {};
    for (Int_t i = 0; i <= bins_const; i++)
    {
    new_bins_const[i] = (from_const + i * width_const);
    }
    
    //x and Y bins for wt with jet pt included --eec_wt and e3c_wt
    Double_t wt_from = 0;
    Double_t wt_to = 1;
    Int_t wt_bins = 200;
    Double_t wt_width = (wt_to-wt_from)/wt_bins;
    Double_t wt_new_bins[201] = {};
    for (Int_t i = 0; i <= wt_bins; i++)
    {
        wt_new_bins[i] = (wt_from + i * wt_width);
    }
    
    //match bins for wt with jet pt included for e3c and Xbins for wt_res_e3c
    Double_t wt_from_e3c = -6;
    Double_t wt_to_e3c = 0;
    Int_t wt_bins_e3c = 400;
    Double_t wt_width_e3c = (wt_to_e3c-wt_from_e3c)/wt_bins_e3c;
    Double_t wt_new_bins_e3c[401] = {};
    for (Int_t i = 0; i <= wt_bins_e3c; i++)
    {
        wt_new_bins_e3c[i] = TMath::Power(10.0, wt_from_e3c + i * wt_width_e3c);
    }

    //Y bins for resolution plot should go from -5 to 5 for eec_wt & e3c_wt
    Double_t res_from = -5;
    Double_t res_to = 5;
    Int_t res_bins = 200;
    Double_t res_width = (res_to-res_from)/res_bins;
    Double_t res_new_bins[201] = {};
    for (Int_t i = 0; i <= res_bins; i++)
    {
        res_new_bins[i] = (res_from + i * res_width);
    }
    
    //Y bins for resolution plot should go from -2 to 2 (fine bins)
    Double_t R_res_from = -1;
    Double_t R_res_to = 1;
    Int_t R_res_bins = 2000;
    Double_t R_res_width = (R_res_to-R_res_from)/R_res_bins;
    Double_t R_res_new_bins[2001] = {};
    for (Int_t i = 0; i <= R_res_bins; i++)
    {
        R_res_new_bins[i] = (R_res_from + i * R_res_width);
    }
  
    //Bins for THIS IS BOTH X AND Y FOR THE wtnojet_match_eec
    Double_t nj_eec_from = -1;
    Double_t nj_eec_to = 4;
    Int_t nj_eec_bins = 200;
    Double_t nj_eec_width = (nj_eec_to-nj_eec_from)/nj_eec_bins;
    Double_t nj_eec_new_bins[201] = {};
    for (Int_t i = 0; i <= nj_eec_bins; i++)
    {
        nj_eec_new_bins[i] = TMath::Power(10.0, nj_eec_from + i * nj_eec_width);
    }
    
    //Bins for THIS IS BOTH X AND Y FOR THE wtnojet_match_e3c
    Double_t nj_e3c_from = -1;
    Double_t nj_e3c_to = 6;
    Int_t nj_e3c_bins = 400;
    Double_t nj_e3c_width = (nj_e3c_to-nj_e3c_from)/nj_e3c_bins;
    Double_t nj_e3c_new_bins[401] = {};
    for (Int_t i = 0; i <= nj_e3c_bins; i++)
    {
        nj_e3c_new_bins[i] = TMath::Power(10.0, nj_e3c_from + i * nj_e3c_width);
    }
    

    
    cout<<"HERE BINS ARE DECLARED!!!!!!!!!!!!!"<<endl;
    jet_pt_hist = new TH1D("jet_pt_hist", "Jet Pt", 22, 10, 120);
    fOutput->Add(jet_pt_hist);
    
    //EEC (data or det level)
    EEC_hist = new TH1D("EEC_hist","EEC", 100, new_bins);
    EEC_hist->RebinX(2);
    fOutput->Add(EEC_hist);
    
    EEC_pt_hist = new TH2D("EEC_pt_hist", "EEC and jet_pt 2D", 100, new_bins, 22, 10, 120);
    EEC_pt_hist->RebinX(2);
    fOutput->Add(EEC_pt_hist);
    
    //EEEC histograms (data or det level)
    E3C_hist = new TH1D("E3C_hist","E3C", 100, new_bins);
    E3C_hist->RebinX(2);
    fOutput->Add(E3C_hist);
    
    E3C_pt_hist = new TH2D("E3C_pt_hist", "EEEC and jet_pt 2D", 100, new_bins, 22, 10, 120);
    E3C_pt_hist->RebinX(2);
    fOutput->Add(E3C_pt_hist);
    
    //MC true and det level histograms 3d (det level fills through the data loop)
    EEC_det_pt_hist_3d = new TH3D("EEC_det_pt_hist", "EEC det and jet_pt 3D", 100, new_bins, 22, new_bins_const, 22, new_bins_const);
    EEC_det_pt_hist_3d->RebinX(2);
    fOutput->Add(EEC_det_pt_hist_3d);
    
    EEC_tru_pt_hist_3d = new TH3D("EEC_tru_pt_hist", "EEC tru and jet_pt 3d", 100, new_bins, 22, new_bins_const, 22, new_bins_const);
    EEC_tru_pt_hist_3d->RebinX(2);
    fOutput->Add(EEC_tru_pt_hist_3d);
    
    E3C_det_pt_hist_3d = new TH3D("E3C_det_pt_hist", "E3C det and jet_pt 3D", 100, new_bins, 22, new_bins_const, 22, new_bins_const);
    E3C_det_pt_hist_3d->RebinX(2);
    fOutput->Add(E3C_det_pt_hist_3d);
    
    E3C_tru_pt_hist_3d = new TH3D("E3C_tru_pt_hist", "E3C tru and jet_pt 3d)", 100, new_bins, 22, new_bins_const, 22, new_bins_const);
    E3C_tru_pt_hist_3d->RebinX(2);
    fOutput->Add(E3C_tru_pt_hist_3d);

    // //////////////////////////////////////////////////////////////////////////////////////////
    if(fCout) cout<<"subtraction histograms for EEC"<<endl;
      //For MIN BIAS SUBTRACTION
    //matched subtracted jet
    h_MJ = new TH3D("h_MJ", "h_MJ", 22, new_bins_const,22, new_bins_const,100, new_bins);//all tracks fake + real
    fOutput->Add(h_MJ);
    h_MJ0 = new TH3D("h_MJ0", "h_MJ0",22, new_bins_const,22, new_bins_const,100, new_bins);//0 real tracks
    fOutput->Add(h_MJ0);
    h_MJ1 = new TH3D("h_MJ1", "h_MJ1", 22, new_bins_const,22, new_bins_const,100, new_bins);//1 real track
    fOutput->Add(h_MJ1);
    h_MJ2 = new TH3D("h_MJ2", "h_MJ2", 22, new_bins_const,22, new_bins_const,100, new_bins);//2 real tracks
    fOutput->Add(h_MJ2);
    
    //matched subtracted jet with the pythia jet pt
    h_MJ_tru = new TH3D("h_MJ_tru", "h_MJ_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);//all tracks fake + real
    fOutput->Add(h_MJ_tru);
    h_MJ0_tru = new TH3D("h_MJ0_tru", "h_MJ0_tru",22, new_bins_const,22, new_bins_const,100, new_bins);//0 real tracks
    fOutput->Add(h_MJ0_tru);
    h_MJ1_tru = new TH3D("h_MJ1_tru", "h_MJ1_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);//1 real track
    fOutput->Add(h_MJ1_tru);
    h_MJ2_tru = new TH3D("h_MJ2_tru", "h_MJ2_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);//2 real tracks
    fOutput->Add(h_MJ2_tru);
    
    //Min bias 1
    h_MB1 = new TH3D("h_MB1", "h_MB1", 22, new_bins_const,22, new_bins_const,100, new_bins);//min bias
    fOutput->Add(h_MB1);
    
    //Jet and min bias and contributions
    h_JMB = new TH3D("h_JMB", "h_JMB", 22, new_bins_const,22, new_bins_const,100, new_bins);//Jet*MINBIAS
    fOutput->Add(h_JMB);
    h_SMB = new TH3D("h_SMB", "h_SMB", 22, new_bins_const,22, new_bins_const,100, new_bins);//JetSignal*MINBIAS
    fOutput->Add(h_SMB);
    h_BMB = new TH3D("h_BMB", "h_BMB", 22, new_bins_const,22, new_bins_const,100, new_bins);//JetBkg*MINBIAS
    fOutput->Add(h_BMB);
    
    //Min bias1 and min bias 2
    h_MB1MB2 = new TH3D("h_MB1MB2", "h_MB1MB2", 22, new_bins_const,22, new_bins_const,100, new_bins);//MinBias1*MinBias2
    fOutput->Add(h_MB1MB2);
    
    h_MB1_tru = new TH3D("h_MB1_tru", "h_MB1_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);//min bias
    fOutput->Add(h_MB1_tru);
    
    //Jet and min bias and contributions
    h_JMB_tru = new TH3D("h_JMB_tru", "h_JMB_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);//Jet*MINBIAS
    fOutput->Add(h_JMB_tru);
    h_SMB_tru = new TH3D("h_SMB_tru", "h_SMB_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);//JetSignal*MINBIAS
    fOutput->Add(h_SMB_tru);
    h_BMB_tru = new TH3D("h_BMB_tru", "h_BMB_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);//JetBkg*MINBIAS
    fOutput->Add(h_BMB_tru);
    
    //Min bias1 and min bias 2
    h_MB1MB2_tru = new TH3D("h_MB1MB2_tru", "h_MB1MB2_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);//MinBias1*MinBias2
    fOutput->Add(h_MB1MB2_tru);
    //////////////////////////////////////////////////
    
    // /////////////////////////////////////////////////
    if(fCout) cout<<"subtraction histograms for E3C"<<endl;
//     /////FOR E3C MIN BIAS SUBTRCTION////////
     h_MJ_e3c = new TH3D("h_MJ_e3c", "h_MJ_e3c", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MJ_e3c);
     h_MJ0_e3c = new TH3D("h_MJ0_e3c", "h_MJ0_e3c", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MJ0_e3c);
     h_MJ1_e3c = new TH3D("h_MJ1_e3c", "h_MJ1_e3c",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MJ1_e3c);
     h_MJ2_e3c = new TH3D("h_MJ2_e3c", "h_MJ2_e3c",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MJ2_e3c);
     h_MJ3_e3c = new TH3D("h_MJ3_e3c", "h_MJ3_e3c", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MJ3_e3c);
     h_MJ_e3c_tru = new TH3D("h_MJ_e3c_tru", "h_MJ_e3c_tru",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MJ_e3c_tru);
     h_MJ0_e3c_tru = new TH3D("h_MJ0_e3c_tru", "h_MJ0_e3c_tru",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MJ0_e3c_tru);
     h_MJ1_e3c_tru = new TH3D("h_MJ1_e3c_tru", "h_MJ1_e3c_tru",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MJ1_e3c_tru);
     h_MJ2_e3c_tru = new TH3D("h_MJ2_e3c_tru", "h_MJ2_e3c_tru",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MJ2_e3c_tru);
     h_MJ3_e3c_tru = new TH3D("h_MJ3_e3c_tru", "h_MJ3_e3c_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MJ3_e3c_tru);

    if(fCout) cout<<"subtraction histograms for E3C 1"<<endl; 

     h_MB1MB1MB1 = new TH3D("h_MB1MB1MB1", "h_MB1MB1MB1",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MB1MB1MB1);
     h_MB1MB1MB1_tru = new TH3D("h_MB1MB1MB1_tru", "h_MB1MB1MB1_tru",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MB1MB1MB1_tru);

    if(fCout) cout<<"subtraction histograms for E3C 2"<<endl; 

     h_JJMB = new TH3D("h_JJMB", "h_JJMB", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_JJMB);
     h_JJMB_tru = new TH3D("h_JJMB_tru", "h_JJMB_tru",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_JJMB_tru);
     h_JMBMB = new TH3D("h_JMBMB", "h_JMBMB",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_JMBMB);
     h_JMBMB_tru = new TH3D("h_JMBMB_tru", "h_JMBMB_tru",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_JMBMB_tru);
     h_MB1MB1MB2 = new TH3D("h_MB1MB1MB2", "h_MB1MB1MB2",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MB1MB1MB2);
     h_MB1MB1MB2_tru = new TH3D("h_MB1MB1MB2_tru", "h_MB1MB1MB2_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MB1MB1MB2_tru);
     h_MB1MB2MB2 = new TH3D("h_MB1MB2MB2", "h_MB1MB2MB2",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MB1MB2MB2);
     h_MB1MB2MB2_tru = new TH3D("h_MB1MB2MB2_tru", "h_MB1MB2MB2_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MB1MB2MB2_tru);

    if(fCout) cout<<"subtraction histograms for E3C 3"<<endl; 

     h_JMB1MB2 = new TH3D("h_JMB1MB2", "h_JMB1MB2", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_JMB1MB2);
     h_JMB1MB2_tru = new TH3D("h_JMB1MB2_tru", "h_JMB1MB2_tru",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_JMB1MB2_tru);
     h_MB1MB2MB3 = new TH3D("h_MB1MB2MB3", "h_MB1MB2MB3",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MB1MB2MB3);
     h_MB1MB2MB3_tru = new TH3D("h_MB1MB2MB3_tru", "h_MB1MB2MB3_tru",22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_MB1MB2MB3_tru);

    if(fCout) cout<<"subtraction histograms for E3C 4"<<endl; 

     h_BMBMB = new TH3D("h_BMBMB", "h_BMBMB", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_BMBMB);
     h_SMBMB = new TH3D("h_SMBMB", "h_SMBMB", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_SMBMB);
     h_BMB1MB2 = new TH3D("h_BMB1MB2", "h_BMB1MB2", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_BMB1MB2);
     h_SMB1MB2 = new TH3D("h_SMB1MB2", "h_SMB1MB2", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_SMB1MB2);
     h_BBMB = new TH3D("h_BBMB", "h_BBMB", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_BBMB);
     h_SBMB = new TH3D("h_SBMB", "h_SBMB", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_SBMB);
     h_SSMB = new TH3D("h_SSMB", "h_SSMB", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_SSMB);

if(fCout) cout<<"subtraction histograms for E3C 5"<<endl;

     h_BMBMB_tru = new TH3D("h_BMBMB_tru", "h_BMBMB_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_BMBMB_tru);
     h_SMBMB_tru = new TH3D("h_SMBMB_tru", "h_SMBMB_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_SMBMB_tru);
     h_BMB1MB2_tru = new TH3D("h_BMB1MB2_tru", "h_BMB1MB2_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_BMB1MB2_tru);
     h_BBMB_tru = new TH3D("h_BBMB_tru", "h_BBMB_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_BBMB_tru);
     h_SBMB_tru = new TH3D("h_SBMB_tru", "h_SBMB_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_SBMB_tru);
     h_SSMB_tru = new TH3D("h_SSMB_tru", "h_SSMB_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_SSMB_tru);
     h_SMB1MB2_tru = new TH3D("h_SMB1MB2_tru", "h_SMB1MB2_tru", 22, new_bins_const,22, new_bins_const,100, new_bins);
     fOutput->Add(h_SMB1MB2_tru);
    

//   // /////////////////////////////////////////////////
    if(fCout) cout<<"num pair histos"<<endl;

    //Num of pairs at det and true level with det pt and tru pt
    N2_det_pt_hist_3d = new TH3D("N2_det_pt_hist", "Num pairs det and jet_pt 3d",100, new_bins, 22, new_bins_const, 22, new_bins_const);
    N2_det_pt_hist_3d->RebinX(2);
    fOutput->Add(N2_det_pt_hist_3d);
    
    N2_tru_pt_hist_3d = new TH3D("N2_tru_pt_hist", "Num pairs tru and jet_pt 3d",100, new_bins, 22, new_bins_const, 22, new_bins_const);
    N2_tru_pt_hist_3d->RebinX(2);
    fOutput->Add(N2_tru_pt_hist_3d);
    
    N3_det_pt_hist_3d = new TH3D("N3_det_pt_hist", "Num pairs det and jet_pt 3d",100, new_bins, 22, new_bins_const, 22, new_bins_const);
    N3_det_pt_hist_3d->RebinX(2);
    fOutput->Add(N3_det_pt_hist_3d);
    
    N3_tru_pt_hist_3d = new TH3D("N3_tru_pt_hist", "Num pairs tru and jet_pt 3d",100, new_bins, 22, new_bins_const, 22, new_bins_const);
    N3_tru_pt_hist_3d->RebinX(2);
    fOutput->Add(N3_tru_pt_hist_3d);
    
    //MC true level histograms (detector level goes through the data loop so don't need det level histograms)
    EEC_tru_match_pt_tru = new TH2D("EEC_pt_tru_match_hist", "EEC and pT for tru matched jets", 100, new_bins, 22, 10, 120);
    EEC_tru_match_pt_tru->RebinX(2);
    fOutput->Add(EEC_tru_match_pt_tru);
    
    E3C_tru_match_pt_tru = new TH2D("E3C_pt_tru_match_hist", "E3C and pT for tru matched jets", 100, new_bins, 22, 10, 120);
    E3C_tru_match_pt_tru->RebinX(2);
    fOutput->Add(E3C_tru_match_pt_tru);

    pt_tru = new TH1D("jet_pt_tru_hist", "Jet Pt", 22, 10, 120);
    fOutput->Add(pt_tru);
    
    R_matrix = new TH2D("R_matrix_jetpt", "Jet pT response matrix", 22, 10, 120, 22, 10, 120);
    fOutput->Add(R_matrix);
    
    JES = new TH2D("JES","Jet energy scale", 22, 10, 120, 200, -2, 8);
    fOutput->Add(JES);
    
    JES_scaled = new TH2D("JES scaled","Jet energy scale for scaled det pT", 22, 10, 120, 200, -2, 8);
    fOutput->Add(JES_scaled);
    
    qpt_det = new TH2D("qpt_det","q/pt vs R det", 100, new_bins, 100, 0, 2);
    fOutput->Add(qpt_det);
    
    qpt_tru = new TH2D("qpt_tru","q/pt vs R tru", 100, new_bins, 100, 0, 2);
    fOutput->Add(qpt_tru);
    
    track_pt_tru = new TH1D("track_pt_tru_hist", "Track Pt", 100, 0, 200);
    fOutput->Add(track_pt_tru);
    
    track_pt_det = new TH1D("track_pt_det_hist", "Track Pt Det", 100, 0, 200);
    fOutput->Add(track_pt_det);
    
    track_pt_matched = new TH2D("track_pt_match_hist", "Match Track Pt", 100, 0, 200, 100, 0, 200);
    fOutput->Add(track_pt_matched);

    // Set track container pointers
    if(fCout) cout<<"Set track container pointers"<<endl;
    fDetectorLevel  = GetTrackContainer(fDetectorLevelName);
    fGeneratorLevel = GetMCParticleContainer(fGeneratorLevelName);
    
    qpt_tru1 = new TH1D("qpt_tru1","q/pt tru1",100, 0, 2);
    fOutput->Add(qpt_tru1);
    
    qpt_tru2 = new TH1D("qpt_tru2","q/pt tru2",100, 0, 2);
    fOutput->Add(qpt_tru2);
    
    if(fCout) cout<<"track matching histograms"<<endl;
//     //Track matching ENC histograms for missed and fake tracks
    if(fMissJetTrack==1)
    {
       if(fCout) cout<<"miss jet track"<<endl;
        eec_Mm = new TH2D("EEC_Mm", "EEC match miss tracks", 100, new_bins, 22, 10, 120);
        eec_Mm->RebinX(2);
        fOutput->Add(eec_Mm);
        
        eec_mm = new TH2D("EEC_mm", "EEC miss miss tracks", 100, new_bins, 22, 10, 120);
        eec_mm->RebinX(2);
        fOutput->Add(eec_mm);
        
        e3c_MMm = new TH2D("E3C_MMm", "E3C match match miss tracks", 100, new_bins, 22, 10, 120);
        e3c_MMm->RebinX(2);
        fOutput->Add(e3c_MMm);
        
        e3c_Mmm = new TH2D("E3C_Mmm", "E3C match miss miss tracks", 100, new_bins, 22, 10, 120);
        e3c_Mmm->RebinX(2);
        fOutput->Add(e3c_Mmm);
        
        e3c_mmm = new TH2D("E3C_mmm", "E3C miss miss miss tracks", 100, new_bins, 22, 10, 120);
        e3c_mmm->RebinX(2);
        fOutput->Add(e3c_mmm);
    }

    if(fFakeJetTrack==1)
    {
       if(fCout) cout<<"fake jet track"<<endl;

        eec_Mf = new TH2D("EEC_Mf", "EEC match fake tracks", 100, new_bins, 22, 10, 120);
        eec_Mf->RebinX(2);
        fOutput->Add(eec_Mf);
        
        eec_ff = new TH2D("EEC_ff", "EEC fake fake tracks", 100, new_bins, 22, 10, 120);
        eec_ff->RebinX(2);
        fOutput->Add(eec_ff);
        
        e3c_MMf = new TH2D("E3C_MMf", "E3C match match fake tracks", 100, new_bins, 22, 10, 120);
        e3c_MMf->RebinX(2);
        fOutput->Add(e3c_MMf);
        
        e3c_Mff = new TH2D("E3C_Mff", "E3C match fake fake tracks", 100, new_bins, 22, 10, 120);
        e3c_Mff->RebinX(2);
        fOutput->Add(e3c_Mff);
        
        e3c_fff = new TH2D("E3C_fff", "E3C fake fake fake tracks", 100, new_bins, 22, 10, 120);
        e3c_fff->RebinX(2);
        fOutput->Add(e3c_fff);
    }
    
    if(fMatchJetTrack==1)
    {
       if(fCout) cout<<"match jet track"<<endl;

       eec_matched_det = new TH2D("EEC_MM_det", "EEC match match det tracks", 100, new_bins, 22, 10, 120);
        eec_matched_det->RebinX(2);
        fOutput->Add(eec_matched_det);

        eec_matched_tru = new TH2D("EEC_MM_tru", "EEC match match tru tracks", 100, new_bins, 22, 10, 120);
        eec_matched_tru->RebinX(2);
        fOutput->Add(eec_matched_tru);

        e3c_matched_det = new TH2D("E3C_MMM_det", "E3C match match match det tracks", 100, new_bins, 22, 10, 120);
        e3c_matched_det->RebinX(2);
        fOutput->Add(e3c_matched_det);

        e3c_matched_tru = new TH2D("E3C_MMM_tru", "E3C match match match tru tracks", 100, new_bins, 22, 10, 120);
        e3c_matched_tru->RebinX(2);
        fOutput->Add(e3c_matched_tru);
        
        R_match_eec_tru = new TH3D("R_match_eec_tru", "Matched Track R eec", 100, new_bins, 100, new_bins, 22, new_bins_const);
        R_match_eec_tru->RebinX(2);
        R_match_eec_tru->RebinY(2);
        fOutput->Add(R_match_eec_tru);

        R_match_e3c_tru = new TH3D("R_match_e3c_tru", "Matched Track R e3c", 100, new_bins, 100, new_bins, 22, new_bins_const);
        R_match_e3c_tru->RebinX(2);
        R_match_e3c_tru->RebinY(2);
        fOutput->Add(R_match_e3c_tru);
        
        R_res_eec_tru = new TH3D("R_res_eec_tru", "R resolution scale EEC", 100, new_bins, 2000, R_res_new_bins, 22, new_bins_const);//R_tru bins, diff bins, jet pT bins
        R_res_eec_tru->RebinX(2);
        fOutput->Add(R_res_eec_tru);

        R_res_e3c_tru = new TH3D("R_res_e3c_tru", "R resolution scale E3C", 100, new_bins, 2000, R_res_new_bins, 22, new_bins_const);
        R_res_e3c_tru->RebinX(2);
        fOutput->Add(R_res_e3c_tru);

        wt_match_eec_tru = new TH3D("wt_match_eec_tru", "Matched Track wt eec", 200, wt_new_bins, 200, wt_new_bins, 22, new_bins_const);
        fOutput->Add(wt_match_eec_tru);

        wt_match_e3c_tru = new TH3D("wt_match_e3c_tru", "Matched Track Wt e3c",400, wt_new_bins_e3c, 400, wt_new_bins_e3c, 22, new_bins_const);
        fOutput->Add(wt_match_e3c_tru);

        wt_res_eec_tru = new TH3D("wt_res_eec_tru", "Weight resolution scale EEC", 200, wt_new_bins, 200, res_new_bins,  22, new_bins_const);
        fOutput->Add(wt_res_eec_tru);

        wt_res_e3c_tru = new TH3D("wt_res_e3c_tru", "Weight resolution scale E3C", 400, wt_new_bins_e3c, 200, res_new_bins,  22, new_bins_const);
        fOutput->Add(wt_res_e3c_tru);
        //-------
        wtnojet_match_eec_tru = new TH3D("wtnojet_match_eec_tru", "Matched Track Wt excluding jet pT EEC", 200, nj_eec_new_bins, 200, nj_eec_new_bins, 22, new_bins_const);
        fOutput->Add(wtnojet_match_eec_tru);

        wtnojet_match_e3c_tru = new TH3D("wtnojet_match_e3c_tru", "Matched Track Wt excluding jet pT E3C", 400, nj_e3c_new_bins, 400, nj_e3c_new_bins, 22, new_bins_const);
//        wtnojet_match_e3c_tru->RebinX(2);
//        wtnojet_match_e3c_tru->RebinY(2);
        fOutput->Add(wtnojet_match_e3c_tru);

        wtnojet_res_eec_tru = new TH3D("wtnojet_res_eec_tru", "Weight resolution (excluding jet pT) scale EEC", 200, nj_eec_new_bins, 2000, R_res_new_bins, 22, new_bins_const);
        fOutput->Add(wtnojet_res_eec_tru);

        wtnojet_res_e3c_tru = new TH3D("wtnojet_res_e3c_tru", "Weight resolution (excluding jet pT) scale E3C", 400, nj_e3c_new_bins, 2000, R_res_new_bins, 22, new_bins_const);
        fOutput->Add(wtnojet_res_e3c_tru);
        //-------
        jetpt_res_w_R = new TH3D("jetpt_res_w_R", "Jet pT resolution with R", 100, new_bins, 200, res_new_bins, 22, new_bins_const);
        fOutput->Add(jetpt_res_w_R);
        
        jetpt_res_w_wt = new TH3D("jetpt_res_w_wt", "Jet pT resolution with wt", 200,  nj_eec_new_bins, 200, res_new_bins, 22, new_bins_const);
        fOutput->Add(jetpt_res_w_wt);
//

    }

    if(fUnfolding==1)
     {
        if(fCout) cout<<"Unfolding tree"<<endl;

         fTreeMatchTracks = new TTree("MatchTracksTree", "MatchTracksTree");
         
         fTreeMatchTracks->Branch("fJet_pt_det", &fJet_pt_det, "fJet_pt_det/D");
         fTreeMatchTracks->Branch("fJet_pt_tru", &fJet_pt_tru, "fJet_pt_tru/D");
         fTreeMatchTracks->Branch("fTrack_pt_det", &fTrack_pt_det, "fTrack_pt_det/D");
         fTreeMatchTracks->Branch("fTrack_pt_tru", &fTrack_pt_tru, "fTrack_pt_tru/D");
         fTreeMatchTracks->Branch("fTrack_eta_det", &fTrack_eta_det, "fTrack_eta_det/D");
         fTreeMatchTracks->Branch("fTrack_eta_tru", &fTrack_eta_tru, "fTrack_eta_tru/D");
         fTreeMatchTracks->Branch("fTrack_phi_det", &fTrack_phi_det, "fTrack_phi_det/D");
         fTreeMatchTracks->Branch("fTrack_phi_tru", &fTrack_phi_tru, "fTrack_phi_tru/D");

         fTreeMatchTracks->Branch("fTrack_pt_miss", &fTrack_pt_miss, "fTrack_pt_miss/D");
         fTreeMatchTracks->Branch("fTrack_eta_miss", &fTrack_eta_miss, "fTrack_eta_miss/D");
         fTreeMatchTracks->Branch("fTrack_phi_miss", &fTrack_phi_miss, "fTrack_phi_miss/D");

         fTreeMatchTracks->Branch("fTrack_pt_fake", &fTrack_pt_fake, "fTrack_pt_fake/D");
         fTreeMatchTracks->Branch("fTrack_eta_fake", &fTrack_eta_fake, "fTrack_eta_fake/D");
         fTreeMatchTracks->Branch("fTrack_phi_fake", &fTrack_phi_fake, "fTrack_phi_fake/D");
    
     }
    
    if(fCout) cout<<"histos and/or trees initialized"<<endl;

    PostData(1, fOutput);
    if (fUnfolding==1) 
    { 
        PostData(2, fTreeMatchTracks);
    }
        
    
}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEECpbpb::Run() {
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEECpbpb::FillHistograms()
{
//    if(fpairfastsim == 1){ComputeDelqpt(); MatchTracks();}//compute track matching and pair eff at track level
    
AliEmcalJet *jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(0);
  // container zero is always the base containe: the data container, the
  // embedded subtracted in the case of embedding or the detector level in case
  // of pythia


  if (fCentSelectOn)
    if ((fCent > fCentMax) || (fCent < fCentMin))
      return 0;

  Float_t rhoVal = 0, rhoMassVal = 0.;
  if (jetCont) {
    jetCont->ResetCurrentID();
    if ((fJetShapeSub == kConstSub) || (fJetShapeSub == kDerivSub)) {
      // rho
      AliRhoParameter *rhoParam = dynamic_cast<AliRhoParameter *>(
								  InputEvent()->FindListObject("RhoSparseR020"));
      if (!rhoParam) {
        Printf("%s: Could not retrieve rho %s (some histograms will be filled "
               "with zero)!",
               GetName(), jetCont->GetRhoName().Data());
      } else
        rhoVal = rhoParam->GetVal();
      // rhom
      AliRhoParameter *rhomParam = dynamic_cast<AliRhoParameter *>(
								   InputEvent()->FindListObject("RhoMassSparseR020"));
      if (!rhomParam) {
        Printf("%s: Could not retrieve rho_m %s (some histograms will be "
               "filled with zero)!",
               GetName(), jetCont->GetRhoMassName().Data());
      } else
        rhoMassVal = rhomParam->GetVal();
    }

    while ((jet1 = jetCont->GetNextAcceptJet())) {
      if (!jet1)
        continue;
      AliEmcalJet *jet2 = 0x0;
      AliEmcalJet *jet3 = 0x0;
      if (jet1->Pt() < fJetPtMin) {continue;} //Cuts on jet_pt
      AliEmcalJet *jetUS = NULL;
      Int_t ifound = 0, jfound = 0;
      Int_t ilab = -1, jlab = -1;

      int sub1 = -1;
      int sub2 = -1;

      //  if (fJetShapeType == kData)
      //       {
      //         if(fCout) cout<<"data loop"<<endl;
      //           ComputeENC(jet1, jetCont);//Computing ENC on raw data
      //           if(fCout) cout<<"computed ENC in data loop"<<endl;
      //       }

      // The embedding mode
      // the matching is done between unsubtracted embedded jets and detector
      // level jets unsubtracted and subtracted jets share the label. Once we
      // identify the corresponding unsubtracted jet, jetUS, then we fetch jet2,
      // which is the matched detector level jet. In the case we are not
      // considering constituent subtraction, then the detector-level matched jet
      // is the one that was directly matched to the base jet1. Then, the
      // particle-level jet jet3 is obtained as the matched one to jet2. In short,
      // there are 2 consecutive matchinges, between particle-level (jet3) and
      // detector-level (jet2) pythia jets and between jet2 and the embedding
      // unsubtracted jet. Note that the matching obtained via ClosestJet is
      // purely geometrical. So below we require a fraction of the probe momentum
      // to be reconstructed in the embedded jet.
      if (fJetShapeType == kDetEmbPartPythia) {

        AliJetContainer *jetContUS = GetJetContainer(2);

        if (fJetShapeSub == kConstSub) {

          for (Int_t i = 0; i < jetContUS->GetNJets(); i++) {
            jetUS = jetContUS->GetJet(i);
            if (jetUS->GetLabel() == jet1->GetLabel()) {
              ifound++;
              if (ifound == 1)
                ilab = i;
            }
          }
          if (ilab == -1)
            continue;
          jetUS = jetContUS->GetJet(ilab);
          jet2 = jetUS->ClosestJet();
        }

	if(fJetShapeSub==kEventSub){
	  
	  jetUS = jet1->ClosestJet();
	  if (!jetUS) continue;
	  jet2 = jetUS->ClosestJet();
	}
	
        if (!(fJetShapeSub == kConstSub) && !(fJetShapeSub == kEventSub)) jet2 = jet1->ClosestJet();
        
	
        if (!jet2) {
	  Printf("jet2 does not exist, returning");
          continue;
        }

        // AliJetContainer *jetContPart=GetJetContainer(3);
        jet3 = jet2->ClosestJet(); //particle level matched jet

        if (!jet3) {
	  Printf("jet3 does not exist, returning");
          continue;
        }

	AliJetContainer *jetContTrue = GetJetContainer(1);
	AliJetContainer *jetContPart = GetJetContainer(3);


    //require a fraction of the probe momentum
    // to be reconstructed in the embedded jet.
  Double_t fraction = 0;
        if (!(fJetShapeSub == kConstSub) && !(fJetShapeSub==kEventSub))
          fraction = jetCont->GetFractionSharedPt(jet1);
	  if ((fJetShapeSub == kConstSub) || (fJetShapeSub == kEventSub))
          fraction = jetContUS->GetFractionSharedPt(jetUS);

        if (fraction < fMinFractionShared)
          continue;
      }
    
    //need this for pp analysis. Mode to run over pythia and match jets 
      if (fJetShapeType == kPythiaDef) {

        AliJetContainer *jetContTrue = GetJetContainer(1);
        AliJetContainer *jetContUS = GetJetContainer(2);
        AliJetContainer *jetContPart = GetJetContainer(3);
	
        if (fJetShapeSub == kConstSub) {
	  
          for (Int_t i = 0; i < jetContUS->GetNJets(); i++) {
            jetUS = jetContUS->GetJet(i);
            if (jetUS->GetLabel() == jet1->GetLabel()) {
              ifound++;
              if (ifound == 1)
                ilab = i;
            }
          }
          if (ilab == -1)
            continue;
          jetUS = jetContUS->GetJet(ilab);
          jet2 = jetUS->ClosestJet();

          if (!jet2) {
            Printf("jet2 does not exist, returning");
            continue;
          }

          for (Int_t j = 0; j < jetContPart->GetNJets(); j++) {
	    
            jet3 = jetContPart->GetJet(j);
            if (!jet3)
              continue;
            if (jet3->GetLabel() == jet2->GetLabel()) {
              jfound++;
              if (jfound == 1)
                jlab = j;
            }
          }
          if (jlab == -1)
            continue;
          jet3 = jetContPart->GetJet(jlab);
          if (!jet3) {
            Printf("jet3 does not exist, returning");
            continue;
          }
        }
        if (!(fJetShapeSub == kConstSub))
          jet3 = jet1->ClosestJet();
        if (!jet3) {
          Printf("jet3 does not exist, returning");
          continue;
        }

        Double_t ptMatch=0;
        Int_t kMatched = 0;
            if (fJetShapeType == kPythiaDef)
            {
                kMatched = 1;
                if (fJetShapeSub == kConstSub)
                    kMatched = 3;
                ptMatch = jet3->Pt();
                //                cout<<"the matched jet "<<jet3->Pt()<<" "<<kMatched<<endl;
                // ComputeEncMC(jet1, jetCont, jet3, kMatched);
            }
    }
 
      Double_t ptSubtracted = 0;
      if (fJetShapeSub == kConstSub || fJetShapeSub == kEventSub) {  
        ptSubtracted = jet1->Pt();
     }
      else if (fJetShapeSub == kDerivSub) {
        ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
      }
    //we will almost exclusively use kNoSub
      else if (fJetShapeSub == kNoSub) {
        if ((fJetShapeType == kData) || (fJetShapeType == kDetEmbPartPythia))
          {ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
          if(fCout) cout<<"data loop"<<endl;
          ComputeENC(jet1, ptSubtracted, jetCont);//Computing ENC on raw data
          if(fCout) cout<<"computed ENC in data loop"<<endl;
          }
        else if ((fJetShapeType == kPythiaDef) || (fJetShapeType == kMCTrue) ||
                 (fJetShapeType == kGenOnTheFly))
          ptSubtracted = jet1->Pt();
      }


    if(ptSubtracted < fPtThreshold) continue;
    if(jet1->Area() < fjetMinArea) continue;

    ///if doEmbedding is true and the pT of the subtracted jet is greater than min pT
     if (fDoEmbedding){
          if(fCout) cout<<"performing embedding corrections"<<endl;
        std::vector<fastjet::PseudoJet> jetConstituents; 
                jetConstituents.clear();
                fastjet::PseudoJet PseudoJetTracks; //Creating a pseudojet object called PseduoTracks
                unsigned int JetconstituentIndex = 0;
                for (auto part: jet1->GetParticleConstituents())
                {
                PseudoJetTracks.reset(part.Px(), part.Py(), part.Pz(), part.E()); 
                const AliVParticle* part2 = part.GetParticle(); //"hack", leave this in , to get the index of the jet from AliPhysics
                //does getconstituentID get indices of background mc particle -- need to check this 
                //if its not a backgroun particle, set ID 
                if(GetConstituentID(JetconstituentIndex, part2, jet1)!= -1){
                    PseudoJetTracks.set_user_index(GetConstituentID(JetconstituentIndex, part2, jet1));
                    JetconstituentIndex++;
                }
                
                if (PseudoJetTracks.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
                jetConstituents.push_back(PseudoJetTracks);
                
                 }

         if(fCout) cout<<"creating all combinations of histograms by running FillEmbJets"<<endl;

        // FillEmbJetsEEC(jetConstituents, jetConstituents, ptSubtracted, jet1->Pt(),true, "sameJet", 1, 1, -1, -1);
        // FillEmbJetsE3C(jetConstituents, jetConstituents, jetConstituents, ptSubtracted, jet1->Pt(),"all", "sameJet", 1, 1,-1, -1, -1);


       
        if(fDoPerpCone){
       if(fCout) cout<<"perp cone"<<endl;

        }
        if(fDoRandCone){
      if(fCout) cout<<"rand cone"<<endl;

        }
     }






    }
  } 
    
    return kTRUE;
}

// //_______________________________________________________________________________________________
double AliAnalysisTaskJetsEECpbpb::delR(const fastjet::PseudoJet& ps1,const fastjet::PseudoJet& ps2)
    {
        double dphi = abs(ps1.phi()-ps2.phi());
        if(dphi>TMath::Pi()){dphi = (2.*TMath::Pi() - dphi);}
        double deta = ps1.eta()-ps2.eta();
        double dR = sqrt(dphi*dphi + deta*deta);
        return dR;
    }

// //_______________________________________________________________________________________________
double AliAnalysisTaskJetsEECpbpb::Calculate_pX( double pT, double eta, double phi)
{
	return(pT*TMath::Cos(phi));
}

double AliAnalysisTaskJetsEECpbpb::Calculate_pY( double pT, double eta, double phi)
{
	return(pT*TMath::Sin(phi));
}

double AliAnalysisTaskJetsEECpbpb::Calculate_pZ( double pT, double eta, double phi)
{
	return( pT*TMath::SinH(eta) );
}

double AliAnalysisTaskJetsEECpbpb::Calculate_E( double pT, double eta, double phi)
{
	double pZ = Calculate_pZ(pT,eta, phi);

	return( TMath::Sqrt(pT*pT + pZ*pZ) );
}

// //________________________________________________________________________
int AliAnalysisTaskJetsEECpbpb::GetConstituentID(int constituentIndex, const AliVParticle* part, AliEmcalJet * jet)
{
  // NOTE: Usually, we would use the global offset defined for the general subtracter extraction task. But we don't want to
  //       depend on that task, so we just define it here locally.
  //Get the id of the particle. If the label is not equal to -1, get the label and assign it to id. If it is -1, get the value of of track at id + 20000
  int id = part->GetLabel() != -1 ? part->GetLabel() : (jet->TrackAt(constituentIndex) + 2000000);
  return id;
}


//_______________________________________________________________________
Double_t AliAnalysisTaskJetsEECpbpb::GetDownscaleWeight(string trigString)
{
  Double_t weight = 1.;
  TString triggerclass;
  if(trigString == "INT7") triggerclass = "CINT7-B-NOPF-CENT";
  else if(trigString == "EJ1") triggerclass = "CEMC7EJ1-B-NOPF-CENTNOTRD";
  else if(trigString == "EJ2") triggerclass = "CEMC7EJ2-B-NOPF-CENT";
  if(triggerclass.Length()) weight = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(triggerclass);
  return weight;
}


//______________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::MatchTracks()
{
    AliVTrack* track;
    for (auto trackIterator : fDetectorLevel->accepted_momentum() )
    {
        track = trackIterator.second;
        Byte_t type = fDetectorLevel->GetTrackType(track);
        if (type <= 2) {
            Double_t sigma = 0;

            if (fIsEsd) {
                continue;
            } else { // AOD

                AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
                if(!aodtrack) AliFatal("Not a standard AOD");

                Int_t label = TMath::Abs(track->GetLabel());

                //            FillDetectorLevelTHnSparse(fCent, track->Eta(), track->Phi(), track->Pt(), sigma, type);

                if (fGeneratorLevel && label > 0) {
                    AliAODMCParticle *part =  fGeneratorLevel->GetAcceptMCParticleWithLabel(label); //label is used to match the generator to detector level
                    if (part) {
                        if (part->GetGeneratorIndex() == 0) {
                            Int_t pdg = TMath::Abs(part->PdgCode());
                            // select charged pions, protons, kaons, electrons, muons
//                            if (pdg == 211 || pdg == 2212 || pdg == 321 || pdg == 11 || pdg == 13) {
                                track_pt_matched->Fill(part->Pt(),track->Pt());
//                            }
                        }
                    }
                }
            }
        }
        else {
            AliError(Form("Track %d has type %d not recognized!", fDetectorLevel->GetCurrentID(), type));
        }
    }
}

//________________________________________________________________________
//bkgindex = -1 for jet background particle, -2 for a thermal cone, -3 for second thermal cone
//pt is true (gen) level jet pT, jetpt is the pT of the embedded subtracted jet
void AliAnalysisTaskJetsEECpbpb::FillEmbJetsEEC(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2, double jetpt, float pt, bool typeSame, std::string type, float MCweight, double n, double bkgIndex1, double bkgIndex2)
{
double w_eec = 0.0;
double w_eec_tru = 0.0;
    if(typeSame){
        int mult = particles.size();
        for (int i = 0; i < mult; i++)
        {
             if(particles.at(i).pt()<fMinENCtrackPt) continue;
            for (int j = i+1; j < mult; j++)
            {
                if(particles.at(j).pt()<fMinENCtrackPt) continue;
                if(n==1)
                {
                    w_eec = (1./(jetpt*jetpt))*(2*particles.at(i).pt()*particles.at(j).pt()*MCweight);
                    w_eec_tru = (1./(pt*pt))*(2*particles.at(i).pt()*particles.at(j).pt()*MCweight);
                }
                
                double dR = delR(particles.at(i), particles.at(j));
                int i1 = particles.at(i).user_index();
                int i2 = particles.at(j).user_index();
                
                if(type == "sameJet")
                {
                    //ss, sb, bb
                    h_MJ->Fill(jetpt, pt, dR, w_eec);
                    if ((i1 == bkgIndex1) && (i2 == bkgIndex1)) h_MJ0->Fill(jetpt, pt, dR, w_eec);//both background
                    else if ((i1 == bkgIndex1) || (i2 == bkgIndex1)) h_MJ1->Fill(jetpt, pt, dR, w_eec);// one is pythia
                    else h_MJ2->Fill(jetpt, pt, dR, w_eec);//both are pythia
                    
                    h_MJ_tru->Fill(jetpt, pt, dR, w_eec_tru); //with true jet pT
                    if ((i1 == bkgIndex1) && (i2 == bkgIndex1)) h_MJ0_tru->Fill(jetpt, pt, dR, w_eec_tru);//both background
                    else if ((i1 == bkgIndex1) || (i2 == bkgIndex1)) h_MJ1_tru->Fill(jetpt, pt, dR, w_eec_tru);// one is pythia
                    else h_MJ2_tru->Fill(jetpt, pt, dR, w_eec_tru);//both are pythia
                }
                if(type == "sameMB"){
                    h_MB1->Fill(jetpt, pt, dR, w_eec);
                    
                    h_MB1_tru->Fill(jetpt, pt, dR, w_eec_tru);
                }
            }
        }
    }//else we are combining two different backgrounds or a jet + background
    else{
        int mult = particles.size();
        int mult2 = particles2.size();
        for (int i = 0; i < mult; i++)
        {
            if(particles.at(i).pt()<fMinENCtrackPt) continue;
            for (int j = 0; j < mult2; j++)
            {
                 if(particles2.at(j).pt()<fMinENCtrackPt) continue;
                if(n==1)
                {
                    
                    w_eec = (1./(jetpt*jetpt))*(particles.at(i).pt()*particles2.at(j).pt()*MCweight);
                    w_eec_tru = (1./(pt*pt))*(particles.at(i).pt()*particles2.at(j).pt()*MCweight);
                }
                double dR = delR(particles.at(i), particles2.at(j));
                
                int i1 = particles.at(i).user_index();
                int i2 = particles2.at(j).user_index();
                
                if(type == "jetMB"){
                    h_JMB->Fill(jetpt, pt, dR, w_eec);
                    if (((i1 == bkgIndex1) && (i2 == bkgIndex2)) || ((i1 == bkgIndex2) && (i2 == bkgIndex1))) h_BMB->Fill(jetpt, pt, dR, w_eec); //background + MB
                    else{h_SMB->Fill(jetpt, pt, dR, w_eec);// pythia + MB
                    }
                    h_JMB_tru->Fill(jetpt, pt, dR, w_eec_tru);
                    if(((i1 == bkgIndex1) && (i2 == bkgIndex2)) || ((i1 == bkgIndex2) && (i2 == bkgIndex1))) h_BMB_tru->Fill(jetpt, pt, dR, w_eec_tru); //background + MB
                    else{h_SMB_tru->Fill(jetpt, pt, dR, w_eec_tru);// pythia + MB
                    }
                }
                if(type == "diffMB"){
                    if ((i1 ==  bkgIndex1) && (i2 ==  bkgIndex2)){
                        h_MB1MB2->Fill(jetpt, pt, dR, w_eec);
                        
                        h_MB1MB2_tru->Fill(jetpt, pt, dR, w_eec_tru);
                    }
                }
            }
        }
    }
}

// //________________________________________________________________________
// //bkgindex = -1 for jet background particle, -2 for first thermal cone, -3 for second thermal cone, -4 for third thermal cone
// //pt is true (gen) level jet pT, jetpt is the pT of the embedded subtracted jet
// void AliAnalysisTaskJetsEECpbpb::FillEmbJetsE3C(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2, std::vector<fastjet::PseudoJet> particles3, double jetpt, float pt,string typeSame, std::string type, float MCweight, double n, double bkgIndex1, double bkgIndex2,  double bkgIndex3)
// {
//     int mult = particles.size(); 
//     int mult2 = particles2.size();
//     int mult3 = particles3.size();
//     double w = -1;
//     double w_tru = -1;
//     double w_iij = -1;
//     double w_tru_iij = -1;
//     double w_ijs = -1;
//     double w_ijs_tru = -1;
//     double w_ijk_tru = -1;
//     double w_ijk = -1;


    
//     if(typeSame == "all"){
//         for (int i = 0; i < mult; i++)
//         {
//             if(particles.at(i).pt()<fMinENCtrackPt) continue;
//             for (int j = i+1; j < mult; j++)
//             {
//                 if(particles.at(j).pt()<fMinENCtrackPt) continue;
                
//                 if(n==1)
//                 {
//                     w = (1./(jetpt*jetpt*jetpt))*(3*particles.at(i).pt()*particles.at(j).pt()*particles.at(j).pt()*MCweight);
//                     w_tru = (1./(pt*pt*pt))*(3*particles.at(i).pt()*particles.at(j).pt()*particles.at(j).pt()*MCweight);
                    
//                     w_iij = (1./(jetpt*jetpt*jetpt))*(3*particles.at(i).pt()*particles.at(i).pt()*particles.at(j).pt()*MCweight);
//                     w_tru_iij = (1./(pt*pt*pt))*(3*particles.at(i).pt()*particles.at(i).pt()*particles.at(j).pt()*MCweight);
//                 }
                
                
                
//                 double dR = delR(particles.at(i), particles.at(j));
//                 int i1 = particles.at(i).user_index();
//                 int i2 = particles.at(j).user_index();
                
//                 if(type =="sameJet")
//                 {
//                     h_MJ_e3c->Fill(jetpt, pt, dR, w);h_MJ_e3c->Fill(jetpt, pt, dR, w_iij);
//                     if ((i1 == -1) && (i2 == -1)) {h_MJ0_e3c->Fill(jetpt, pt, dR, w);h_MJ0_e3c->Fill(jetpt, pt, dR, w_iij);}//fake fake fake
//                     else if (i2 == -1) {h_MJ1_e3c->Fill(jetpt, pt, dR, w);h_MJ2_e3c->Fill(jetpt, pt, dR, w_iij);}//real fake fake for ijj & real real fake for iij
//                     else if (i1 == -1) {h_MJ2_e3c->Fill(jetpt, pt, dR, w);h_MJ1_e3c->Fill(jetpt, pt, dR, w_iij);}//real real fake
//                     else {h_MJ3_e3c->Fill(jetpt, pt, dR, w);h_MJ3_e3c->Fill(jetpt, pt, dR, w_iij);}//both are pythia
                    
//                     h_MJ_e3c_tru->Fill(jetpt, pt, dR, w_tru);h_MJ_e3c_tru->Fill(jetpt, pt, dR, w_tru_iij);
//                     if ((i1 == -1) && (i2 == -1)) {h_MJ0_e3c_tru->Fill(jetpt, pt, dR, w_tru);h_MJ0_e3c_tru->Fill(jetpt, pt, dR, w_tru_iij);}//fake fake fake
//                     else if (i2 == -1) {h_MJ1_e3c_tru->Fill(jetpt, pt, dR, w_tru);h_MJ2_e3c_tru->Fill(jetpt, pt, dR, w_tru_iij);}//real fake fake
//                     else if (i1 == -1) {h_MJ2_e3c_tru->Fill(jetpt, pt, dR, w_tru);h_MJ1_e3c_tru->Fill(jetpt, pt, dR, w_tru_iij);}//real real fake
//                     else {h_MJ3_e3c_tru->Fill(jetpt, pt, dR, w_tru);h_MJ3_e3c_tru->Fill(jetpt, pt, dR, w_tru_iij);}//both are pythia
                    
//                 }
                
                
//                 if(type=="sameMB"){
//                     if(i1 == -2 && i2 == -2){
//                         h_MB1MB1MB1->Fill(jetpt, pt, dR, w);
//                         h_MB1MB1MB1->Fill(jetpt, pt, dR, w_iij);
                        
//                         h_MB1MB1MB1_tru->Fill(jetpt, pt, dR, w_tru);
//                         h_MB1MB1MB1_tru->Fill(jetpt, pt, dR, w_tru_iij);
                        
                        
//                     }          
//                 }       
//                 for (int s = j+1; s < mult; s++)
//                 {
//                     if(particles.at(s).pt()<fMinENCtrackPt) continue;
                    
//                     if(n==1)
//                     {
//                         w_ijs = (1./(jetpt*jetpt*jetpt))*(6*particles.at(i).pt()*particles.at(j).pt()*particles.at(s).pt()*MCweight);
//                         w_ijs_tru = (1./(pt*pt*pt))*(6*particles.at(i).pt()*particles.at(j).pt()*particles.at(s).pt()*MCweight);
//                     }
                    
//                     if(n==2)
//                     {
//                         w_ijs = (6*MCweight*(pow(particles.at(i).pt(),2))*(pow(particles.at(j).pt(),2))*(pow(particles.at(s).pt(),2)))/(pow(jetpt,6));
//                         w_ijs_tru = (6*MCweight*(pow(particles.at(i).pt(),2))*(pow(particles.at(j).pt(),2))*(pow(particles.at(s).pt(),2)))/(pow(pt,6));
//                     }
                    
//                     double dR_ij = delR(particles.at(i), particles.at(j));
//                     double dR_js = delR(particles.at(j), particles.at(s));
//                     double dR_is = delR(particles.at(s), particles.at(i));
                    
//                     double R_L = -1; 
//                     if(dR_ij>dR_js && dR_ij>dR_is){R_L = dR_ij;}
//                     else if(dR_js>dR_ij && dR_js>dR_is){R_L = dR_js;}
//                     else{R_L = dR_is;}
                    
//                     int i1 = particles.at(i).user_index();
//                     int i2 = particles.at(j).user_index();
//                     int i3 = particles.at(s).user_index();
                    
//                     if(type=="sameJet"){
//                         h_MJ_e3c->Fill(jetpt, pt, R_L, w_ijs);
//                         if ((i1 == -1) && (i2 == -1) && (i3 ==-1)) h_MJ0_e3c->Fill(jetpt, pt, R_L, w_ijs);//fake fake fake
//                         else if (((i1 == -1) && (i2 == -1) && (i3 != -1)) ||((i1 == -1) && (i2 != -1) && (i3 == -1)) ||((i1 != -1) && (i2 == -1) && (i3 == -1))) h_MJ1_e3c->Fill(jetpt, pt, R_L, w_ijs);//fake fake real
//                         else if (((i1 == -1) && (i2 != -1) && (i3 != -1)) ||((i1 != -1) && (i2 != -1) && (i3 == -1)) ||((i1 != -1) && (i2 == -1) && (i3 != -1))) h_MJ2_e3c->Fill(jetpt, pt, R_L, w_ijs);//real real fake
//                         else h_MJ3_e3c->Fill(jetpt, pt, R_L, w_ijs);//real real real
                      
                        
//                         h_MJ_e3c_tru->Fill(jetpt, pt, R_L, w_ijs_tru);
//                         if ((i1 == -1) && (i2 == -1) && (i3 ==-1)) h_MJ0_e3c_tru->Fill(jetpt, pt, R_L, w_ijs_tru);//fake fake fake
//                         else if (((i1 == -1) && (i2 == -1) && (i3 != -1)) ||((i1 == -1) && (i2 != -1) && (i3 == -1)) ||((i1 != -1) && (i2 == -1) && (i3 == -1))) h_MJ1_e3c_tru->Fill(jetpt, pt, R_L, w_ijs_tru);//fake fake real
//                         else if (((i1 == -1) && (i2 != -1) && (i3 != -1)) ||((i1 != -1) && (i2 != -1) && (i3 == -1)) ||((i1 != -1) && (i2 == -1) && (i3 != -1))) h_MJ2_e3c_tru->Fill(jetpt, pt, R_L, w_ijs_tru);//real real fake
//                         else h_MJ3_e3c_tru->Fill(jetpt, pt, R_L, w_ijs_tru);//real real real
//                     }
                    
//                     if(type=="sameMB"){
//                         if(i1 == -2 && i2 == -2 && i3 == -2){
//                             h_MB1MB1MB1->Fill(jetpt, pt, dR, w_ijs);
//                             h_MB1MB1MB1->Fill(jetpt, pt, dR, w_ijs);
                            
//                             h_MB1MB1MB1_tru->Fill(jetpt, pt, dR, w_ijs_tru);
//                             h_MB1MB1MB1_tru->Fill(jetpt, pt, dR, w_ijs_tru);
//                         }          
//                     }    
//                 }
//             }
//         }
//     }
//     else if (typeSame == "two"){
//         //particles2 and particles3 is the same list
        
//         for (int i = 0; i < mult; i++)
//         {
//             if(particles.at(i).pt()<fMinENCtrackPt) continue;
//             int i1 = particles.at(i).user_index();
            
//             for (int  j = 0; j < mult2; j++)
//             {
                
//                 int i2 = particles2.at(j).user_index();
//                 if(particles2.at(j).pt()<fMinENCtrackPt) continue;
                
//                 double w_twosame = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles2.at(j).pt()*particles2.at(j).pt()*MCweight);
//                 double w_twosame_tru = (1./(pt*pt*pt))*(particles.at(i).pt()*particles2.at(j).pt()*particles2.at(j).pt()*MCweight);
                
//                 double R_L = delR(particles.at(i), particles2.at(j));
                
                
//                 if(i1==-2 && (i2 > 0 || i2 == -1)) h_JJMB->Fill(jetpt, pt, R_L,w_twosame);
//                 if(i1==-2 && i2 == -1 && i2 == -1) h_BBMB->Fill(jetpt, pt, R_L,w_twosame);
//                 if(i1==-2 && i2 > 0 && i2 > 0) h_SSMB->Fill(jetpt, pt, R_L,w_twosame);
                
//                 if((i1== -1 || i1 > 0) && (i2 == -2)) h_JMBMB->Fill(jetpt, pt, R_L,w_twosame);
//                 if((i1== -1)  && (i2 == -2) && (i2==-2)) h_BMBMB->Fill(jetpt, pt, R_L,w_twosame);
//                 if((i1 > 0) && (i2 == -2) && (i2==-2)) h_SMBMB->Fill(jetpt, pt, R_L,w_twosame);
                
//                 if(i1 == -3 && i2 == -2) h_MB1MB1MB2->Fill(jetpt, pt, R_L, w_twosame);
//                 if(i1 == -2 && i2 == -3) h_MB1MB2MB2->Fill(jetpt, pt, R_L, w_twosame);
                
//                 if(i1==-2 && (i2 > 0 || i2 == -1)) h_JJMB_tru->Fill(jetpt, pt, R_L,w_twosame_tru);
//                 if((i1== -1 || i1 > 0) && (i2 == -2)) h_JMBMB_tru->Fill(jetpt, pt, R_L,w_twosame_tru);
//                 if(i1 == -3 && i2 == -2) h_MB1MB1MB2_tru->Fill(jetpt, pt, R_L, w_twosame_tru);
//                 if(i1 == -2 && i2 == -3) h_MB1MB2MB2_tru->Fill(jetpt, pt, R_L, w_twosame_tru);
                
//                 for (int k = j+1; k < mult2; k++)
//                 {
//                     if(particles2.at(k).pt()<fMinENCtrackPt) continue;
//                     int i3 = particles2.at(k).user_index();
                    
//                     if(n==1)
//                     {
//                         w_ijk = (1./(jetpt*jetpt*jetpt))*(2*particles.at(i).pt()*particles2.at(j).pt()*particles2.at(k).pt()*MCweight);
//                         w_ijk_tru = (1./(pt*pt*pt))*(2*particles.at(i).pt()*particles2.at(j).pt()*particles2.at(k).pt()*MCweight);
//                     }
                    
//                     if(n==2)
//                     {
//                         w_ijs = (MCweight*(pow(particles.at(i).pt(),2))*(pow(particles2.at(j).pt(),2))*(pow(particles2.at(k).pt(),2)))/(pow(jetpt,6));
//                         w_ijs_tru = (MCweight*(pow(particles.at(i).pt(),2))*(pow(particles2.at(j).pt(),2))*(pow(particles2.at(k).pt(),2)))/(pow(pt,6));
//                     }
                    
//                     double dR_ij = delR(particles.at(i), particles2.at(j));
//                     double dR_js = delR(particles2.at(j), particles2.at(k));
//                     double dR_is = delR(particles2.at(k), particles.at(i));
                    
//                     double R_L = -1; 
//                     if(dR_ij>dR_js && dR_ij>dR_is){R_L = dR_ij;}
//                     else if(dR_js>dR_ij && dR_js>dR_is){R_L = dR_js;}
//                     else{R_L = dR_is;}
                    
                    
//                     if(i1== -2 && (i2 > 0 || i2 == -1) && (i3 > 0 || i3 == -1)) h_JJMB->Fill(jetpt, pt, R_L,w_ijk);
//                     if(i1== -2 && i2 == -1 && i3 == -1) h_BBMB->Fill(jetpt, pt, R_L,w_ijk);
//                     if(i1== -2 && i2 > 0 && i3 > 0) h_SSMB->Fill(jetpt, pt, R_L,w_ijk);
//                     if(i1== -2 && ((i2 == -1 && i3 > 0)||(i2 > 0 && i3 == -1))) h_SBMB->Fill(jetpt, pt, R_L,w_ijk);
                    
//                     if((i1== -1 || i1 > 0) && (i2 == -2) && (i3== -2)) h_JMBMB->Fill(jetpt, pt, R_L,w_ijk);
//                     if((i1== -1)  && (i2 == -2) && (i3== -2)) h_BMBMB->Fill(jetpt, pt, R_L,w_ijk);
//                     if((i1 > 0) && (i2 == -2) && (i3== -2)) h_SMBMB->Fill(jetpt, pt, R_L,w_ijk);
                    
//                     if(i1 == -3 && i2 == -2 && i3 == -2) h_MB1MB1MB2->Fill(jetpt, pt, R_L, w_ijk);
//                     if(i1 == -2 && i2 == -3 && i3 == -3) h_MB1MB2MB2->Fill(jetpt, pt, R_L, w_ijk);
                    
                    
//                     if(i1== -2 && (i2 > 0 || i2 == -1) && (i3 > 0 || i3 == -1)) h_JJMB_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
//                     if(i1== -2 && i2 == -1 && i3 == -1) h_BBMB_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
//                     if(i1== -2 && i2 > 0 && i3 > 0) h_SSMB_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
//                     if(i1== -2 && ((i2 == -1 && i3 > 0)||(i2 > 0 && i3 == -1))) h_SBMB_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
                    
//                     if((i1== -1 || i1 > 0) && (i2 == -2) && (i3== -2)) h_JMBMB_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
//                     if((i1== -1)  && (i2 == -2) && (i3== -2)) h_BMBMB_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
//                     if((i1 > 0) && (i2 == -2) && (i3== -2)) h_SMBMB_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
                    
//                     if(i1 == -3 && i2 == -2 && i3 == -2) h_MB1MB1MB2_tru->Fill(jetpt, pt, R_L, w_ijk_tru);
//                     if(i1 == -2 && i2 == -3 && i3 == -3) h_MB1MB2MB2_tru->Fill(jetpt, pt, R_L, w_ijk_tru);
                    
//                 }
//             }
//         }
//     }
//     else{
        
//         for (int i = 0; i < mult; i++)
//         {
//             if(particles.at(i).pt()<fMinENCtrackPt) continue;
//             int i1 = particles.at(i).user_index();
            
//             for (int  j = 0; j < mult2; j++)
//             {
//                 if(particles2.at(j).pt()<fMinENCtrackPt) continue;
//                 int i2 = particles2.at(j).user_index();
                
//                 for (int k = 0; k < mult3; k++)
//                 {
//                     if(particles3.at(k).pt()<fMinENCtrackPt) continue;
//                     int i3 = particles3.at(k).user_index();
                    
//                     if(n==1)
//                     {
//                         w_ijk = (1./(jetpt*jetpt*jetpt))*(particles.at(i).pt()*particles2.at(j).pt()*particles3.at(k).pt()*MCweight);
//                         w_ijk_tru = (1./(pt*pt*pt))*(particles.at(i).pt()*particles2.at(j).pt()*particles3.at(k).pt()*MCweight);
//                     }
                    
//                     if(n==2)
//                     {
//                         w_ijs = (MCweight*(pow(particles.at(i).pt(),2))*(pow(particles2.at(j).pt(),2))*(pow(particles3.at(k).pt(),2)))/(pow(jetpt,6));
//                         w_ijs_tru = (MCweight*(pow(particles.at(i).pt(),2))*(pow(particles2.at(j).pt(),2))*(pow(particles3.at(k).pt(),2)))/(pow(pt,6));
//                     }
                    
//                     double dR_ij = delR(particles.at(i), particles2.at(j));
//                     double dR_js = delR(particles2.at(j), particles3.at(k));
//                     double dR_is = delR(particles.at(i), particles3.at(k));
                    
//                     double R_L = -1; 
//                     if(dR_ij>dR_js && dR_ij>dR_is){R_L = dR_ij;}
//                     else if(dR_js>dR_ij && dR_js>dR_is){R_L = dR_js;}
//                     else{R_L = dR_is;}
                    
//                     if((i1 > 0 || i1 == -1) && i2== -2 && i3 == -3) h_JMB1MB2->Fill(jetpt, pt, R_L,w_ijk);
//                     if(i1 == -1 && i2== -2 && i3 == -3) h_BMB1MB2->Fill(jetpt, pt, R_L,w_ijk);
//                     if(i1 > 0 && i2== -2 && i3 == -3) h_SMB1MB2->Fill(jetpt, pt, R_L,w_ijk);
                    
//                     if(i1 == -2 && i2== -3 && i3 == -4) h_MB1MB2MB3->Fill(jetpt, pt, R_L,w_ijk);
                    
//                     if((i1 > 0 || i1 == -1) && i2== -2 && i3 == -3) h_JMB1MB2_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
//                     if(i1 == -1 && i2== -2 && i3 == -3) h_BMB1MB2_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
//                     if(i1 > 0 && i2== -2 && i3 == -3) h_SMB1MB2_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
                    
//                     if(i1 == -2 && i2== -3 && i3 == -4) h_MB1MB2MB3_tru->Fill(jetpt, pt, R_L,w_ijk_tru);
                    
                    
//                 }
//             }
//         }
//     }
// }

// //______________________________________________________________________
// //for embedded, det level and truth jets 
// std::vector<fastjet::PseudoJet> AliAnalysisTaskJetsEECpbpb::FindThermalCone(AliEmcalJet *fJetEmb, AliJetContainer *fJetContEmb, double ptSub, AliEmcalJet *fJet, AliJetContainer *fJet_detCont, AliEmcalJet *fJet_tru, AliJetContainer *fJet_truCont, std::string axisType)
// {
//     Float_t jet_embpt = fJetEmb->Pt();
//     Float_t jet_embphi = fJetEmb->Phi();
//     Float_t jet_embeta = fJetEmb->Eta();
//     Float_t jet_embArea = fJetEmb->AreaPt();
//     Float_t jet_embptsub = ptSub;


//     Double_t dPhi1 = 999.;
//     Double_t dPhi2 = 999.;
//     Double_t dPhi = 999.;
//     Double_t dEta = 999.;
//     Double_t Axis1 = 999, Axis2 = 999;
//     Double_t fEtaMC = 999;
//     fastjet::PseudoJet PseudoTracks; //Creating a pseudojet object called PseduoTracks

//     if (fDoRandCone){
//       double OffsetRandom = (1 + fRandom->Rndm()) * TMath::Pi()/3; // creating random Phi in the range of PI/3 and 2PI/3
//       Axis1 = jet_embphi + OffsetRandom;  // adding the random Phi to the Phi of leading jet, to get the Phi of random cone
//       Axis2 = jet_embphi - OffsetRandom; // mirroring it, so you get Phi of the 2nd random cone
//       if(Axis1 > TMath::TwoPi()) Axis1 -= TMath::TwoPi(); // checking if it's larger than 2PI (leading jet phi is distributed in range of 0 and 2PI
//       if(Axis2 < 0) Axis2 += TMath::TwoPi(); // checking if 2nd cone is smaller than 0
//     }

//     if (fDoPerpCone){
//       Axis1 = ((jet_embphi + (TMath::Pi() / 2.)) > TMath::TwoPi()) ? jet_embphi - ((3. / 2.) * TMath::Pi()) : jet_embphi + (TMath::Pi() / 2.);
//       Axis2 = ((jet_embphi - (TMath::Pi() / 2.)) < 0) ? jet_embphi + ((3. / 2.) * TMath::Pi()) : jet_embphi - (TMath::Pi() / 2.);
//     }

//     AliParticleContainer * partCont = 0;
//     TIter nextPartCont(&fParticleCollArray);
//     while ((partCont = static_cast<AliParticleContainer*>(nextPartCont()))) {
//       AliParticleIterableMomentumContainer itcont = partCont->accepted_momentum();
//       for (AliParticleIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
//         AliVTrack * particle = static_cast<AliVTrack*>(it->second); //partIter.second);
//         AliAODTrack* trackReal = (AliAODTrack*)(particle);
        
//         if (trackReal==NULL) {
//         //   if (fUseCouts) std::cout << "Didn't have a reco track" << std::endl;
//           continue;
//         }

//         fEtaMC = trackReal->Eta();
//         if (TMath::Abs(fEtaMC) > fEtaCutValue) {
//           continue;
//         }

//         if (trackReal->Pt() < fMinENCtrackPt) {
//           continue;
//         }
                      
//         Float_t mod_track_phi = trackReal->Phi() + TMath::Pi();
//         //Check if the track is within the R=0.4 cone in eta, phi
//         dPhi1 = TMath::Abs(mod_track_phi - Axis1);
//         dPhi1 = (dPhi1 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi1 : dPhi1;
//         dPhi2 = TMath::Abs(mod_track_phi - Axis2);
//         dPhi2 = (dPhi2 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi2 : dPhi2;
//         dEta = jet_embeta - trackReal->Eta();

//         if(axisType=="plus"){
//           dPhi = dPhi1;
//         }
//         else{
//           dPhi = dPhi2;
//         }
//         // if ((TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) > 0.4) && (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) > 0.4)) continue; //scale the yields by 1/(2Ncones*piR^2) offline
        
//         if ((TMath::Sqrt(dPhi * dPhi + dEta * dEta) > 0.4)) continue; //scale the yields by 1/(2Ncones*piR^2) offline
//         // PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E());
//         PseudoTracks.reset(trackReal->Px(), trackReal->Py(), trackReal->Pz(), trackReal->E());
//       }
//     }
    
// }
// //______________________________________________________________________
// void AliAnalysisTaskJetsEECpbpb::ComputeEncMC(AliEmcalJet *fJet, AliJetContainer *fJetCont, AliEmcalJet *fJet_tru, Int_t km)
//     //(jet, jet container, vector of jets, vector of constituents within each jet?)
// {
//     //Det level
//     std::vector<fastjet::PseudoJet> fConstituents; //Is a pseudojet object with constituents of the jet
//     fConstituents.clear();
//     //This snippet of code is getting particles within a single jet (fjet) and turning them into pseudojet objects so that fastjet capabilities can be used
//     fastjet::PseudoJet PseudoTracks; //Creating a pseudojet object called PseduoTracks
//     unsigned int constituentIndex = 0;
//     int constituentCharge_det = 0;
//     //The line below gets constituent particles within fjet. C++ syntax[ for (auto elem : container)    // capture elements by value ]
//     for (auto part: fJet->GetParticleConstituents())
//     {
//         PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E()); //part is the constituent at that point in the loop, part keeps getting redefined in each step.
//         const AliVParticle* part2 = part.GetParticle(); //"hack", leave this in , to get the index of the jet from AliPhysics
//         PseudoTracks.set_user_index(GetConstituentID(constituentIndex, part2, fJet)); //leave this in for the same reason as above
// //        PseudoTracks.set_user_index(GetConstituentCharge(constituentCharge_det, part2, fJet_tru));
//         if (PseudoTracks.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
//         fConstituents.push_back(PseudoTracks);
//         constituentIndex++;
//     }

//     //Truth level
//     AliJetContainer *jetCont = GetJetContainer(km); //get the container for the matched true level jet
//     std::vector<fastjet::PseudoJet> fConstituents_tru; //Is a pseudojet object with constituents of the jet
//     fConstituents_tru.clear();
//     fastjet::PseudoJet PseudoTracks_tru; //Creating a pseudojet object called PseduoTracks
//     unsigned int constituentIndex_tru = 0;
//     int constituentCharge_tru = 0;
//     for (auto part_tru: fJet_tru->GetParticleConstituents())
//     {
//         PseudoTracks_tru.reset(part_tru.Px(), part_tru.Py(), part_tru.Pz(), part_tru.E()); //part is the constituent at that point in the loop, part keeps getting redefined in each step.
//         const AliVParticle* part_tru2 = part_tru.GetParticle(); //"hack", leave this in , to get the index of the jet from AliPhysics
//         PseudoTracks_tru.set_user_index(GetConstituentID(constituentIndex_tru, part_tru2, fJet_tru)); //leave this in for the same reason as above
// //        PseudoTracks_tru.set_user_index(GetConstituentCharge(constituentCharge_tru, part_tru2, fJet_tru)); //leave this in for the same reason as above
//         if (PseudoTracks_tru.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
//         fConstituents_tru.push_back(PseudoTracks_tru);
//         constituentIndex_tru++;
//     }

//     double jet_pt = fJet->Pt();
//     double jet_pt_tru = fJet_tru->Pt();
//     pt_tru->Fill(jet_pt_tru); //filling histogram with momentum of jets
    
    
//     std::vector<Double_t> R_dist_tru, R_dist_det;
//     std::vector<int> det_index,tru_index;
//     std::vector<fastjet::PseudoJet> matchtracks_det, matchtracks_tru, misstracks ,faketracks,matchtracks_tru_cop;
    

//     for(int j=0; j<int(fConstituents.size()); j++)
//     {
//         det_index.push_back(fConstituents[j].user_index());
//     }
//     for(int i=0; i<int(fConstituents_tru.size()); i++)
//     {
//         tru_index.push_back(fConstituents_tru[i].user_index());
//     }
    
//     if(fFakeJetTrack==1)
//     {
//         for(int j = 0;j<int(fConstituents.size()); j++)
//         {
//             int valueToCheck = fConstituents[j].user_index();
//             auto it = std::find(tru_index.begin(), tru_index.end(), valueToCheck);
//             if (it != tru_index.end()){matchtracks_det.push_back(fConstituents[j]);}
//             else {faketracks.push_back(fConstituents[j]);}
//         }
//     }
    
//     if(fMissJetTrack==1)
//     {
//         for(int i=0; i<int(fConstituents_tru.size()); i++)
//         {
//             int valueToCheck = fConstituents_tru[i].user_index();
//             auto it = std::find(det_index.begin(), det_index.end(), valueToCheck);
//             if (it != det_index.end()){matchtracks_tru.push_back(fConstituents_tru[i]);matchtracks_tru_cop.push_back(fConstituents_tru[i]);}
//             else {misstracks.push_back(fConstituents_tru[i]);}
//         }
//     }
    
//     //sort the tracks according to their indices
//     struct CustomComparator {
//         bool operator()(const fastjet::PseudoJet& track1, const fastjet::PseudoJet& track2) const {
//             return track1.user_index() < track2.user_index();
//         }
//     };
    
    
//     std::sort(matchtracks_tru.begin(), matchtracks_tru.end(), CustomComparator());
//     std::sort(matchtracks_det.begin(), matchtracks_det.end(), CustomComparator());
    


//     if(fMatchJetTrack==1)
//     {
        
//         for(int j = 0; j < int(matchtracks_det.size()); j++) // match
//         {
//             for(int s = j+1; s < int(matchtracks_det.size()); s++)
//             {
                
//                 double ee_det = matchtracks_det[j].pt()*matchtracks_det[s].pt()/(pow(jet_pt,2));
//                 double ee_tru = matchtracks_tru[j].pt()*matchtracks_tru[s].pt()/(pow(jet_pt_tru,2));
                
//                 double delRdet = delR(matchtracks_det[j],matchtracks_det[s]);
//                 double delRtru = delR(matchtracks_tru[j],matchtracks_tru[s]);
                
//                 double R_det = matchtracks_det[j].delta_R(matchtracks_det[s]); //for checks
//                 double R_tru = matchtracks_tru[j].delta_R(matchtracks_tru[s]); //for checks
                
//                 eec_matched_det->Fill(delRdet,jet_pt,2*ee_det);
//                 eec_matched_tru->Fill(delRtru,jet_pt_tru,2*ee_tru);
                
//                 R_match_eec_tru->Fill(delRtru,delRdet,jet_pt_tru);
//                 R_res_eec_tru->Fill(delRtru,(delRtru-delRdet)/(delRtru),jet_pt_tru);
                
//                 R_match_e3c_tru->Fill(delRtru,delRdet,jet_pt_tru);
//                 R_res_e3c_tru->Fill(delRtru,(delRtru-delRdet)/(delRtru),jet_pt_tru);
                
//                 R_match_eec_tru_rap->Fill(R_tru,R_det,jet_pt_tru);
//                 //                R_res_eec_tru_debug->Fill(R_tru,(R_tru-R_det)/(R_tru),jet_pt_tru);
                
//                 wt_match_eec_tru->Fill(ee_tru,ee_det,jet_pt_tru);
//                 wt_res_eec_tru->Fill(ee_tru,(ee_tru-ee_det)/ee_tru,jet_pt_tru);
                
//                 wtnojet_match_eec_tru->Fill(2*matchtracks_tru[j].pt()*matchtracks_tru[s].pt(),2*matchtracks_det[j].pt()*matchtracks_det[s].pt(),jet_pt_tru);
//                 double wt_diff_eec_nojet = ((matchtracks_tru[j].pt()*matchtracks_tru[s].pt())-(matchtracks_det[j].pt()*matchtracks_det[s].pt()))/(matchtracks_tru[j].pt()*matchtracks_tru[s].pt());
//                 wtnojet_res_eec_tru->Fill(matchtracks_tru[j].pt()*matchtracks_tru[s].pt(),wt_diff_eec_nojet,jet_pt_tru);
                

//                 double eee_det_jjs = (matchtracks_det[j].pt()*matchtracks_det[j].pt()*matchtracks_det[s].pt())/(pow(jet_pt,3));
//                 double eee_tru_jjs = (matchtracks_tru[j].pt()*matchtracks_tru[j].pt()*matchtracks_tru[s].pt())/(pow(jet_pt_tru,3));
                
//                 double eee_det = matchtracks_det[j].pt()*matchtracks_det[s].pt()*matchtracks_det[s].pt()/(pow(jet_pt,3));
//                 double eee_tru = matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[s].pt()/(pow(jet_pt_tru,3));

                
                  
//                 R_match_eec_tru->Fill(delRtru,delRdet,jet_pt_tru);
//                 R_res_eec_tru->Fill(delRtru,(delRtru-delRdet)/(delRtru),jet_pt_tru);

//                 wt_match_eec_tru->Fill(ee_tru,ee_det,jet_pt_tru);
//                 wt_res_eec_tru->Fill(ee_tru,(ee_tru-ee_det)/ee_tru,jet_pt_tru);
                
                
//                 jetpt_res_w_R->Fill(delRtru,(jet_pt_tru-jet_pt)/jet_pt_tru,jet_pt_tru);
//                 jetpt_res_w_wt->Fill(matchtracks_tru[j].pt()*matchtracks_tru[s].pt(),(jet_pt_tru-jet_pt)/jet_pt_tru,jet_pt_tru);


//                 R_match_e3c_tru->Fill(delRtru,delRdet,jet_pt_tru);
//                 R_res_e3c_tru->Fill(delRtru,(delRtru-delRdet)/delRtru,jet_pt_tru);
                
//                 e3c_matched_tru->Fill(3*eee_tru,delRtru,jet_pt_tru);
//                 e3c_matched_det->Fill(3*eee_det,delRdet,jet_pt);

//                 e3c_matched_tru->Fill(3*eee_tru_jjs,delRtru,jet_pt_tru);
//                 e3c_matched_det->Fill(3*eee_det_jjs,delRdet,jet_pt);
                 
//                 wt_match_e3c_tru->Fill(eee_tru,eee_det,jet_pt_tru);
//                 wt_match_e3c_tru->Fill(eee_tru_jjs,eee_det_jjs,jet_pt_tru);

//                 wt_res_e3c_tru->Fill(eee_tru,(eee_tru-eee_det)/eee_tru,jet_pt_tru);
//                 wt_res_e3c_tru->Fill(eee_tru_jjs,(eee_tru_jjs-eee_det_jjs)/eee_tru_jjs,jet_pt_tru);
                
//                 wtnojet_match_e3c_tru->Fill(3*matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[s].pt(),3*matchtracks_det[j].pt()*matchtracks_det[s].pt()*matchtracks_det[s].pt(),jet_pt_tru);
//                 double wt_diff_e3c_nojet = ((matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[s].pt())-(matchtracks_det[j].pt()*matchtracks_det[s].pt()*matchtracks_det[s].pt()))/(matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[s].pt());
//                 wtnojet_res_e3c_tru->Fill(3*matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[s].pt(),wt_diff_e3c_nojet,jet_pt_tru);
                
//                 for(int m = s+1; m < int(matchtracks_det.size()) ; m++) // match
//                 {
                   
//                     double eee_jsm_tru = matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[m].pt()/(pow(jet_pt_tru,3));
//                     double eee_jsm = matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[m].pt()/(pow(jet_pt_tru,3));
                    
//                     double Rjs = delR(matchtracks_det[j], matchtracks_det[s]);
//                     double Rjm = delR(matchtracks_det[j], matchtracks_det[m]);
//                     double Rsm = delR(matchtracks_det[s], matchtracks_det[m]);
                    
//                     std::vector<double> R_dist= {Rjs, Rjm, Rsm};
//                     int max_R_index = std::distance(R_dist.begin(), std::max_element(R_dist.begin(), R_dist.end()));
                    
//                     double Rjs_tru = delR(matchtracks_tru[j], matchtracks_tru[s]);
//                     double Rjm_tru = delR(matchtracks_tru[j], matchtracks_tru[m]);
//                     double Rsm_tru = delR(matchtracks_tru[s], matchtracks_tru[m]);
                    
//                     std::vector<double> R_dist_tru = {Rjs_tru, Rjm_tru, Rsm_tru};
//                     int max_R_index_tru = std::distance(R_dist_tru.begin(), std::max_element(R_dist_tru.begin(), R_dist_tru.end()));
                    
//                     R_match_e3c_tru->Fill(R_dist_tru[max_R_index_tru],R_dist[max_R_index],jet_pt_tru);
//                     R_res_e3c_tru->Fill(R_dist_tru[max_R_index_tru],(R_dist_tru[max_R_index_tru]-R_dist[max_R_index])/(R_dist_tru[max_R_index_tru]),jet_pt_tru);
                    
//                     e3c_matched_tru->Fill(R_dist[max_R_index_tru],jet_pt_tru,6*eee_jsm_tru);
//                     e3c_matched_det->Fill(R_dist[max_R_index],jet_pt,6*eee_jsm);
                    
//                     wt_match_e3c_tru->Fill(eee_jsm_tru,eee_jsm,jet_pt_tru);
//                     wt_res_e3c_tru->Fill(eee_jsm_tru,(eee_jsm_tru-eee_jsm)/eee_jsm_tru,jet_pt_tru);
                    
//                     double wt_diff_e3c_nojet_jsm = ((matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[m].pt())-(matchtracks_det[j].pt()*matchtracks_det[s].pt()*matchtracks_det[m].pt()))/(matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[m].pt());//diff/true
                    
//                     wtnojet_match_e3c_tru->Fill(matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[m].pt(),matchtracks_det[j].pt()*matchtracks_det[s].pt()*matchtracks_det[m].pt(),jet_pt_tru);
//                     wtnojet_res_e3c_tru->Fill(matchtracks_tru[j].pt()*matchtracks_tru[s].pt()*matchtracks_tru[m].pt(),wt_diff_e3c_nojet_jsm,jet_pt_tru);
                    
//                 }
//             }
//         }
//     }

//     matchtracks_tru.clear();
//     matchtracks_det.clear();
//     faketracks.clear();
//     misstracks.clear();
//     tru_index.clear();
//     det_index.clear();

//     if(fpTcorr == 1)
//     {   jet_pt = jet_pt/(0.85); //applying JES correction to jet pT spectra to study spectra shape dependence of ENC
//         jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
        
//         double diff_scaled = (jet_pt - jet_pt_tru)/jet_pt_tru;
//         JES_scaled->Fill(jet_pt_tru,diff_scaled);
        
//         R_matrix->Fill(jet_pt_tru,jet_pt); //Filling the response matrix}
//     }
    
//     else
//     {
//         jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
//         R_matrix->Fill(jet_pt_tru,jet_pt); //Filling the response matrix
        
//         double diff = (jet_pt - jet_pt_tru)/jet_pt_tru;
//         JES->Fill(jet_pt_tru,diff);
//     }
    
//     //Det level
//     //Initializing objects
//     std::vector<Double_t> R_dist;
    
//     //Truth level
//     //Initializing objects
//     std::vector<Double_t> R_dist_part;

    
//     //Looping over the det jet
//      for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
//         {
//             for(int s=j+1; s<int(fConstituents.size()) ; s++)
//             {
//                 double eee_jss =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
//                 double eee_jjs =((3*fConstituents[j].pt()*fConstituents[j].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
//                 double deltaR = delR(fConstituents[j],fConstituents[s]);

//                 double ee_js = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));

//                 if(fpaircut == 1)
//                 {
//                     double j_eta = fConstituents[j].eta();
//                     double s_eta = fConstituents[s].eta();
//                     double del_js_eta = abs(j_eta-s_eta);
//                     if (del_js_eta < 0.008) continue;
//                     else
//                     {

//                         E3C_hist->Fill(deltaR,eee_jss);
//                         E3C_pt_hist->Fill(deltaR,jet_pt,eee_jss);

//                         N3_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru);
//                         E3C_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru, eee_jss);

//                         E3C_hist->Fill(deltaR,eee_jjs);
//                         E3C_pt_hist->Fill(deltaR,jet_pt,eee_jjs);

//                         N3_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru);
//                         E3C_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru, eee_jjs);


//                         EEC_hist->Fill(deltaR,ee_js);
//                         EEC_pt_hist->Fill(deltaR,jet_pt,ee_js);

//                         N2_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru);
//                         EEC_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru, ee_js);

//                     }
//                 }
//                 else
//                 {
//                      E3C_hist->Fill(deltaR,eee_jss);
//                         E3C_pt_hist->Fill(deltaR,jet_pt,eee_jss);

//                         N3_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru);
//                         E3C_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru, eee_jss);

//                         E3C_hist->Fill(deltaR,eee_jjs);
//                         E3C_pt_hist->Fill(deltaR,jet_pt,eee_jjs);

//                         N3_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru);
//                         E3C_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru, eee_jjs);


//                         EEC_hist->Fill(deltaR,ee_js);
//                         EEC_pt_hist->Fill(deltaR,jet_pt,ee_js);

//                         N2_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru);
//                         EEC_det_pt_hist_3d->Fill(deltaR, jet_pt, jet_pt_tru, ee_js);
//                 }
//                 for(int k=s+1; k<int(fConstituents.size()) ; k++)
//                 { 
//                 double eee_jsk =((6*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[k].pt())/(pow(jet_pt,3)));
              
//                 double dR_js = delR(fConstituents[j],fConstituents[s]);
//                 double dR_sk = delR(fConstituents[s],fConstituents[k]);
//                 double dR_kj = delR(fConstituents[j],fConstituents[k]);

//                 double RL = -1;

//                 if(dR_js>dR_sk && dR_js>dR_sk){RL = dR_js;}
//                 else if(dR_sk>dR_js && dR_sk>dR_kj){RL = dR_sk;}
//                 else{RL = dR_kj;}

//                  if(fpaircut == 1)
//                 {
//                     double s_eta = fConstituents[s].eta();
//                     double k_eta = fConstituents[k].eta();
//                     double j_eta = fConstituents[j].eta();
//                     double del_js_eta = abs(j_eta-s_eta);
//                     double del_sk_eta = abs(s_eta-k_eta);
//                     double del_kj_eta = abs(k_eta-j_eta);
//                     if (del_sk_eta < 0.008 || del_js_eta < 0.008 || del_kj_eta < 0.008) continue;
//                     else
//                     {
//                         E3C_hist->Fill(RL,eee_jsk);
//                         E3C_pt_hist->Fill(RL,jet_pt,eee_jsk);

//                         N3_det_pt_hist_3d->Fill(RL, jet_pt, jet_pt_tru);
//                         E3C_det_pt_hist_3d->Fill(RL, jet_pt, jet_pt_tru, eee_jsk);

//                     }
//                 }
//                 else
//                 {
//                         E3C_hist->Fill(RL,eee_jsk);
//                         E3C_pt_hist->Fill(RL,jet_pt,eee_jsk);

//                         N3_det_pt_hist_3d->Fill(RL, jet_pt, jet_pt_tru);
//                         E3C_det_pt_hist_3d->Fill(RL, jet_pt, jet_pt_tru, eee_jsk);

//                 }
//             }
//         }
//     }

   
//     //Looping over truth level jet
//   for(int j=0; j<int(fConstituents_tru.size()); j++)  //looping over constituents of the fConstituents object
//         {
//             for(int s=j+1; s<int(fConstituents_tru.size()) ; s++)
//             {
//                 double eee_jss =((3*fConstituents_tru[j].pt()*fConstituents_tru[s].pt()*fConstituents_tru[s].pt())/(pow(jet_pt_tru,3)));
//                 double eee_jjs =((3*fConstituents_tru[j].pt()*fConstituents_tru[j].pt()*fConstituents_tru[s].pt())/(pow(jet_pt_tru,3)));
//                 double deltaR = delR(fConstituents_tru[j],fConstituents_tru[s]);

//                 double ee_js = (2*fConstituents_tru[j].pt()*fConstituents_tru[s].pt())/(pow((jet_pt_tru),2));
              
//                 //Pair cut
//                 if(fpaircut == 1)
//                 {
//                     double j_eta = fConstituents_tru[j].eta();
//                     double s_eta = fConstituents_tru[s].eta();
//                     double del_js_eta = abs(j_eta-s_eta);
//                     if (del_js_eta < 0.008) continue;
//                     else
//                     {

//                         E3C_tru_match_pt_tru->Fill(deltaR, jet_pt_tru, eee_jss);
//                         E3C_tru_match_pt_tru->Fill(deltaR, jet_pt_tru, eee_jjs);

//                         N3_tru_pt_hist_3d->Fill(deltaR, jet_pt_tru, jet_pt);
//                         E3C_tru_pt_hist_3d->Fill(deltaR, jet_pt_tru, jet_pt, eee_jss);
//                         E3C_tru_pt_hist_3d->Fill(deltaR, jet_pt_tru, jet_pt, eee_jjs);

//                         EEC_tru_match_pt_tru->Fill(deltaR, jet_pt_tru, ee_js);

//                         N2_tru_pt_hist_3d->Fill(deltaR, jet_pt_tru, jet_pt);
//                         EEC_tru_pt_hist_3d->Fill(deltaR, jet_pt_tru, jet_pt, ee_js);
//                     }
//                 }
//                 else
//                 {
//                         E3C_tru_match_pt_tru->Fill(deltaR, jet_pt_tru, eee_jss);
//                         E3C_tru_match_pt_tru->Fill(deltaR, jet_pt_tru, eee_jjs);

//                         N3_tru_pt_hist_3d->Fill(deltaR, jet_pt_tru, jet_pt);
//                         E3C_tru_pt_hist_3d->Fill(deltaR, jet_pt_tru, jet_pt, eee_jss);
//                         E3C_tru_pt_hist_3d->Fill(deltaR, jet_pt_tru, jet_pt, eee_jjs);

//                         EEC_tru_match_pt_tru->Fill(deltaR, jet_pt_tru, ee_js);

//                         N2_tru_pt_hist_3d->Fill(deltaR, jet_pt_tru, jet_pt);
//                         EEC_tru_pt_hist_3d->Fill(deltaR, jet_pt_tru, jet_pt, ee_js);
//                 }
//               for(int k=s+1; k<int(fConstituents_tru.size()) ; k++)
//                 { 
//                 double eee_jsk =((6*fConstituents_tru[j].pt()*fConstituents_tru[s].pt()*fConstituents_tru[k].pt())/(pow(jet_pt_tru,3)));
              
//                 double dR_js = delR(fConstituents_tru[j],fConstituents_tru[s]);
//                 double dR_sk = delR(fConstituents_tru[s],fConstituents_tru[k]);
//                 double dR_kj = delR(fConstituents_tru[j],fConstituents_tru[k]);

//                 double RL = -1;

//                 if(dR_js>dR_sk && dR_js>dR_sk){RL = dR_js;}
//                 else if(dR_sk>dR_js && dR_sk>dR_kj){RL = dR_sk;}
//                 else{RL = dR_kj;}

//                  if(fpaircut == 1)
//                 {
//                     double s_eta = fConstituents_tru[s].eta();
//                     double k_eta = fConstituents_tru[k].eta();
//                     double j_eta = fConstituents_tru[j].eta();
//                     double del_js_eta = abs(j_eta-s_eta);
//                     double del_sk_eta = abs(s_eta-k_eta);
//                     double del_kj_eta = abs(k_eta-j_eta);
//                     if (del_sk_eta < 0.008 || del_js_eta < 0.008 || del_kj_eta < 0.008) continue;
//                     else
//                     {
//                         E3C_tru_match_pt_tru->Fill(RL, jet_pt_tru, eee_jsk);

//                         N3_tru_pt_hist_3d->Fill(RL, jet_pt_tru, jet_pt);
//                         E3C_tru_pt_hist_3d->Fill(RL, jet_pt_tru, jet_pt, eee_jsk);

//                     }
//                 }
//                 else
//                 {
//                         E3C_tru_match_pt_tru->Fill(RL, jet_pt_tru, eee_jsk);

//                         N3_tru_pt_hist_3d->Fill(RL, jet_pt_tru, jet_pt);
//                         E3C_tru_pt_hist_3d->Fill(RL, jet_pt_tru, jet_pt, eee_jsk);
                        

//                 }
//             }
//         }
//     }

    
//     //catch (fastjet::Error)
//     // {
//     //    AliError(" [w] FJ Exception caught.");
//     //    // return -1;
//     // } //end error message
//     return;
// }

    
//ENC computation-------------------------------------------------------
void AliAnalysisTaskJetsEECpbpb::ComputeENC(AliEmcalJet *fJet, double ptSub, AliJetContainer *fJetCont)
{
    //General ENC computation: Need to loop over jets, then loop within a single jet to compute EECs.
    //NOTE: Already in the jet loop defined in the Fill Histogram function. This puts you in the jet loop. The event loop stuff is being taken care of by AliAnalysis manager
    //fjet is my single jet. AliEmCalJet is pointing to the fJet object.
    
    std::vector<fastjet::PseudoJet> fConstituents; //Is a pseudojet object with constituents of the jet
    fConstituents.clear();
    //This snippet of code is getting particles within a single jet (fjet) and turning them into pseudojet objects so that fastjet capabilities can be used
    fastjet::PseudoJet PseudoTracks; //Creating a pseudojet object called PseduoTracks
    unsigned int constituentIndex = 0;
    //The line below gets constituent particles within fjet. C++ syntax[ for (auto elem : container)    // capture elements by value ]
    for (auto part: fJet->GetParticleConstituents())
    {
        PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E()); //part is the constituent at that point in the loop, part keeps getting redefined in each step.
        const AliVParticle* part2 = part.GetParticle(); //"hack", leave this in , to get the index of the jet from AliPhysics
        PseudoTracks.set_user_index(GetConstituentID(constituentIndex, part2, fJet)); //leave this in for the same reason as above
        if (PseudoTracks.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
        fConstituents.push_back(PseudoTracks);
        constituentIndex++;
    }
    
        //Initializing objects for det level
        std::vector<Double_t> R_dist;
       
        //Looping over the jet
        // double jet_pt = fJet->Pt();
        double jet_pt = ptSub;
        
        if(fpTcorr == 1)
         {   jet_pt = jet_pt/(0.85); //applying JES correction to jet pT spectra to study spectra shape dependence of ENC
             jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
             
         }
         
         else
         {
             jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
         
         }
        
          //Looping over the det jet
     for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
        {
            for(int s=j+1; s<int(fConstituents.size()) ; s++)
            {
                double eee_jss =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                double eee_jjs =((3*fConstituents[j].pt()*fConstituents[j].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                double deltaR = delR(fConstituents[j],fConstituents[s]);

                double ee_js = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));

                if(fpaircut == 1)
                {
                    double j_eta = fConstituents[j].eta();
                    double s_eta = fConstituents[s].eta();
                    double del_js_eta = abs(j_eta-s_eta);
                    if (del_js_eta < 0.008) continue;
                    else
                    {

                        E3C_hist->Fill(deltaR,eee_jss);
                        E3C_pt_hist->Fill(deltaR,jet_pt,eee_jss);

                        E3C_hist->Fill(deltaR,eee_jjs);
                        E3C_pt_hist->Fill(deltaR,jet_pt,eee_jjs);


                        EEC_hist->Fill(deltaR,ee_js);
                        EEC_pt_hist->Fill(deltaR,jet_pt,ee_js);


                    }
                }
                else
                {
                        E3C_hist->Fill(deltaR,eee_jss);
                        E3C_pt_hist->Fill(deltaR,jet_pt,eee_jss);

                        E3C_hist->Fill(deltaR,eee_jjs);
                        E3C_pt_hist->Fill(deltaR,jet_pt,eee_jjs);

                        EEC_hist->Fill(deltaR,ee_js);
                        EEC_pt_hist->Fill(deltaR,jet_pt,ee_js);

                      
                }
            
           for(int k=s+1; k<int(fConstituents.size()) ; k++)
            { 
                double eee_jsk =((6*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[k].pt())/(pow(jet_pt,3)));
              
                double dR_js = delR(fConstituents[j],fConstituents[s]);
                double dR_sk = delR(fConstituents[s],fConstituents[k]);
                double dR_kj = delR(fConstituents[j],fConstituents[k]);

                double RL = -1;

                if(dR_js>dR_sk && dR_js>dR_sk){RL = dR_js;}
                else if(dR_sk>dR_js && dR_sk>dR_kj){RL = dR_sk;}
                else{RL = dR_kj;}

                 if(fpaircut == 1)
                {
                    double s_eta = fConstituents[s].eta();
                    double k_eta = fConstituents[k].eta();
                    double j_eta = fConstituents[j].eta();
                    double del_js_eta = abs(j_eta-s_eta);
                    double del_sk_eta = abs(s_eta-k_eta);
                    double del_kj_eta = abs(k_eta-j_eta);
                    if (del_sk_eta < 0.008 || del_js_eta < 0.008 || del_kj_eta < 0.008) continue;
                    else
                    {
                        E3C_hist->Fill(RL,eee_jsk);
                        E3C_pt_hist->Fill(RL,jet_pt,eee_jsk);

                    }
                }
                else
                {
                        E3C_hist->Fill(RL,eee_jsk);
                        E3C_pt_hist->Fill(RL,jet_pt,eee_jsk);

                }
             }
            
        }
    }


    
    //catch (fastjet::Error)
    // {
    //    AliError(" [w] FJ Exception caught.");
    //    // return -1;
    // } //end error message
    return;
}


//_________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::RunChanged(Int_t newrun)
{
  if(fStoreTrig)
  {
    auto downscalehandler = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance();
    if(downscalehandler->GetCurrentRun() != newrun)
    {
      downscalehandler->SetRun(newrun);
    }
  }
}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEECpbpb::RetrieveEventObjects()
{
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;
  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskJetsEECpbpb::Terminate(Option_t *)
{
  // Called once at the end of the analysis.
  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available");
  //   return;
  // }
}

