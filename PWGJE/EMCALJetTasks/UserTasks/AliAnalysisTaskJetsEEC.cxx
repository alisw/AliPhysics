//Task for ENC in pp collisions
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
#include "AliAnalysisTaskJetsEEC.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetsEEC)

//________________________________________________________________________
AliAnalysisTaskJetsEEC::AliAnalysisTaskJetsEEC(): AliAnalysisTaskEmcalJet("AliAnalysisTaskJetsEEC", kTRUE),
fContainer(0), fMinFractionShared(0), fJetShapeType(kData),
fJetShapeSub(kNoSub), fJetSelection(kInclusive), fPtThreshold(-9999.), fMinENCtrackPt(1.0), fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE), fCheckResolution(kFALSE),
fMinPtConst(1), fHardCutoff(0), fDoTwoTrack(kFALSE), fCutDoubleCounts(kTRUE),
fPowerAlgo(1), fPhiCutValue(0.02),
fEtaCutValue(0.02), fDerivSubtrOrder(0),
fStoreDetLevelJets(0), fDoFillEncMC(kFALSE), fStoreTrig(kFALSE), fpTcorr(0), fpaircut(0), fpairfastsim(1), fMatchJetTrack(1),fMissJetTrack(1),fFakeJetTrack(1), fMaxPtTrack(0), fGeneratorLevelName(), fDetectorLevelName(), fGeneratorLevel(0), fDetectorLevel(0), fJet_truCont(0), jet_pt_hist(0), EEC_hist(0), EEC_pt_hist(0), E3C_hist(0), E3C_pt_hist(0),  EEC_det_pt_hist_3d(0), EEC_tru_pt_hist_3d(0), E3C_det_pt_hist_3d(0), E3C_tru_pt_hist_3d(0), N2_det_pt_hist_3d(0), N2_tru_pt_hist_3d(0), N3_det_pt_hist_3d(0), N3_tru_pt_hist_3d(0), EEC_det_match_pt_det(0), EEC_tru_match_pt_tru(0), E3C_det_match_pt_det(0), E3C_tru_match_pt_tru(0), pt_tru(0), pt_tru_match(0), pt_det(0), pt_det_match(0), test_hist(0), R_matrix(0), JES(0), JES_scaled(0), JER(0), pair_det_EEC(0), pair_tru_EEC(0), pair_det_E3C(0), pair_tru_E3C(0), qpt_tru(0), qpt_det(0), track_pt_tru(0), track_pt_det(0), track_pt_matched(0), R_match_eec(0), wt_match_eec(0), R_match_e3c(0), wt_match_e3c(0), qpt_tru1(0), qpt_tru2(0),pt_tru1(0), pt_tru2(0), eec_Mm(0), eec_mm(0), e3c_MMm(0), e3c_Mmm(0), e3c_mmm(0), eec_Mf(0), eec_ff(0), e3c_MMf(0), e3c_Mff(0), e3c_fff(0), eec_matched_det(0), eec_matched_tru(0), e3c_matched_det(0), e3c_matched_tru(0), wt_res_eec(0), wt_res_e3c(0), R_res_eec(0), R_res_e3c(0), wtnojet_match_eec(0), wtnojet_match_e3c(0), wtnojet_res_eec(0), wtnojet_res_e3c(0)
{
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}


//________________________________________________________________________
AliAnalysisTaskJetsEEC::AliAnalysisTaskJetsEEC(const char *name): AliAnalysisTaskEmcalJet(name, kTRUE), fContainer(0),
fMinFractionShared(0), fJetShapeType(kData), fJetShapeSub(kNoSub),
fJetSelection(kInclusive), fPtThreshold(-9999.), fMinENCtrackPt(1.0), fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE), fCheckResolution(kFALSE),
fMinPtConst(1), fHardCutoff(0), fDoTwoTrack(kFALSE), fCutDoubleCounts(kTRUE),
fPowerAlgo(1), fPhiCutValue(0.02),
fEtaCutValue(0.02), fDerivSubtrOrder(0),
fStoreDetLevelJets(0), fDoFillEncMC(kFALSE), fStoreTrig(kFALSE), fpTcorr(0), fpaircut(0), fpairfastsim(1), fMatchJetTrack(1),fMissJetTrack(1),fFakeJetTrack(1), fMaxPtTrack(0), fGeneratorLevelName(), fDetectorLevelName(), fGeneratorLevel(0), fDetectorLevel(0), fJet_truCont(0), jet_pt_hist(0), EEC_hist(0), EEC_pt_hist(0), E3C_hist(0), E3C_pt_hist(0),  EEC_det_pt_hist_3d(0), EEC_tru_pt_hist_3d(0), E3C_det_pt_hist_3d(0), E3C_tru_pt_hist_3d(0), N2_det_pt_hist_3d(0), N2_tru_pt_hist_3d(0), N3_det_pt_hist_3d(0), N3_tru_pt_hist_3d(0), EEC_det_match_pt_det(0), EEC_tru_match_pt_tru(0), E3C_det_match_pt_det(0), E3C_tru_match_pt_tru(0), pt_tru(0), pt_tru_match(0), pt_det(0), pt_det_match(0), test_hist(0), R_matrix(0), JES(0), JES_scaled(0), JER(0), pair_det_EEC(0), pair_tru_EEC(0), pair_det_E3C(0), pair_tru_E3C(0), qpt_tru(0), qpt_det(0), track_pt_tru(0), track_pt_det(0), track_pt_matched(0), R_match_eec(0), wt_match_eec(0), R_match_e3c(0), wt_match_e3c(0), qpt_tru1(0), qpt_tru2(0),pt_tru1(0), pt_tru2(0), eec_Mm(0), eec_mm(0), e3c_MMm(0), e3c_Mmm(0), e3c_mmm(0), eec_Mf(0), eec_ff(0), e3c_MMf(0), e3c_Mff(0), e3c_fff(0), eec_matched_det(0), eec_matched_tru(0), e3c_matched_det(0), e3c_matched_tru(0), wt_res_eec(0), wt_res_e3c(0), R_res_eec(0), R_res_e3c(0), wtnojet_match_eec(0), wtnojet_match_e3c(0), wtnojet_res_eec(0), wtnojet_res_e3c(0)

{
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}


//________________________________________________________________________
AliAnalysisTaskJetsEEC::~AliAnalysisTaskJetsEEC() {
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetsEEC::UserCreateOutputObjects() {

// Create user output.
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    TH1::AddDirectory(oldStatus);
    
    const char *nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
    
    Double_t from = -4;
    Double_t to = 0;
    Int_t bins = 100;
    Double_t width = (to-from)/bins;
    Double_t new_bins[101] = {};
    for (Int_t i = 0; i <= bins; i++)
    {
        new_bins[i] = TMath::Power(10.0, from + i * width);
    }

    Double_t from_const = 15;
    Double_t to_const = 120;
    Int_t bins_const = 21;
    Double_t width_const = (to_const-from_const)/bins_const;
    Double_t new_bins_const[22] = {};
    for (Int_t i = 0; i <= bins_const; i++)
    {
    new_bins_const[i] = (from_const + i * width_const);
    }

    
    jet_pt_hist = new TH1D("jet_pt_hist", "Jet Pt", 21, 15, 120);
    fOutput->Add(jet_pt_hist);
    
    //EEC (data or det level)
    EEC_hist = new TH1D("EEC_hist","EEC", 100, new_bins);
    EEC_hist->RebinX(2);
    fOutput->Add(EEC_hist);
    
    EEC_pt_hist = new TH2D("EEC_pt_hist", "EEC and jet_pt 2D", 100, new_bins, 21, 15, 120);
    EEC_pt_hist->RebinX(2);
    fOutput->Add(EEC_pt_hist);
    
    //EEEC histograms (data or det level)
    E3C_hist = new TH1D("E3C_hist","E3C", 100, new_bins);
    E3C_hist->RebinX(2);
    fOutput->Add(E3C_hist);
    
    E3C_pt_hist = new TH2D("E3C_pt_hist", "EEEC and jet_pt 2D", 100, new_bins, 21, 15, 120);
    E3C_pt_hist->RebinX(2);
    fOutput->Add(E3C_pt_hist);
    
    //MC true and det level histograms 3d (det level fills through the data loop)
    EEC_det_pt_hist_3d = new TH3D("EEC_det_pt_hist", "EEC det and jet_pt 3D", 100, new_bins, 21, new_bins_const, 21, new_bins_const);
    EEC_det_pt_hist_3d->RebinX(2);
    fOutput->Add(EEC_det_pt_hist_3d);
    
    EEC_tru_pt_hist_3d = new TH3D("EEC_tru_pt_hist", "EEC tru and jet_pt 3d", 100, new_bins, 21, new_bins_const, 21, new_bins_const);
    EEC_tru_pt_hist_3d->RebinX(2);
    fOutput->Add(EEC_tru_pt_hist_3d);
    
    E3C_det_pt_hist_3d = new TH3D("E3C_det_pt_hist", "E3C det and jet_pt 3D", 100, new_bins, 21, new_bins_const, 21, new_bins_const);
    E3C_det_pt_hist_3d->RebinX(2);
    fOutput->Add(E3C_det_pt_hist_3d);
    
    E3C_tru_pt_hist_3d = new TH3D("E3C_tru_pt_hist", "E3C tru and jet_pt 3d)", 100, new_bins, 21, new_bins_const, 21, new_bins_const);
    E3C_tru_pt_hist_3d->RebinX(2);
    fOutput->Add(E3C_tru_pt_hist_3d);
    
//
    //Num of pairs at det and true level with det pt and tru pt
    N2_det_pt_hist_3d = new TH3D("N2_det_pt_hist", "Num pairs det and jet_pt 3d",100, new_bins, 21, new_bins_const, 21, new_bins_const);
    N2_det_pt_hist_3d->RebinX(2);
    fOutput->Add(N2_det_pt_hist_3d);
    
    N2_tru_pt_hist_3d = new TH3D("N2_tru_pt_hist", "Num pairs tru and jet_pt 3d",100, new_bins, 21, new_bins_const, 21, new_bins_const);
    N2_tru_pt_hist_3d->RebinX(2);
    fOutput->Add(N2_tru_pt_hist_3d);
    
    N3_det_pt_hist_3d = new TH3D("N3_det_pt_hist", "Num pairs det and jet_pt 3d",100, new_bins, 21, new_bins_const, 21, new_bins_const);
    N3_det_pt_hist_3d->RebinX(2);
    fOutput->Add(N3_det_pt_hist_3d);
    
    N3_tru_pt_hist_3d = new TH3D("N3_tru_pt_hist", "Num pairs tru and jet_pt 3d",100, new_bins, 21, new_bins_const, 21, new_bins_const);
    N3_tru_pt_hist_3d->RebinX(2);
    fOutput->Add(N3_tru_pt_hist_3d);
    
    //MC true level histograms (detector level goes through the data loop so don't need det level histograms)
    EEC_tru_match_pt_tru = new TH2D("EEC_pt_tru_match_hist", "EEC and pT for tru matched jets", 100, new_bins, 21, 15, 120);
    EEC_tru_match_pt_tru->RebinX(2);
    fOutput->Add(EEC_tru_match_pt_tru);
    
    E3C_tru_match_pt_tru = new TH2D("E3C_pt_tru_match_hist", "E3C and pT for tru matched jets", 100, new_bins, 21, 15, 120);
    E3C_tru_match_pt_tru->RebinX(2);
    fOutput->Add(E3C_tru_match_pt_tru);

    pt_tru = new TH1D("jet_pt_tru_hist", "Jet Pt", 21, 15, 120);
    fOutput->Add(pt_tru);
    
    test_hist = new TH1D("test_hist_id", "test constituent id", 3000001, -1, 3000000);
    fOutput->Add(test_hist);
    
    R_matrix = new TH2D("R_matrix_jetpt", "Jet pT response matrix", 21, 15, 120, 21, 15, 120);
    fOutput->Add(R_matrix);
    
    JES = new TH2D("JES","Jet energy scale", 21, 15, 120, 200, -2, 8);
    fOutput->Add(JES);
    
    JES_scaled = new TH2D("JES scaled","Jet energy scale for scaled det pT", 21, 15, 120, 200, -2, 8);
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
    fDetectorLevel  = GetTrackContainer(fDetectorLevelName);
    fGeneratorLevel = GetMCParticleContainer(fGeneratorLevelName);
    
    R_match_eec = new TH2D("R_match_eec", "Matched Track R", 100, new_bins, 100, new_bins);
    R_match_eec->RebinX(2);
    fOutput->Add(R_match_eec);
    
    wt_match_eec = new TH2D("wt_match_eec", "Matched Track Wt", 200, 0 ,1 , 200, 0, 1);
    fOutput->Add(wt_match_eec);
    
    R_match_e3c = new TH2D("R_match_e3c", "Matched Track R e3c", 100, new_bins, 100, new_bins);
    R_match_e3c->RebinX(2);
    fOutput->Add(R_match_e3c);
    
    wt_match_e3c = new TH2D("wt_match_e3c", "Matched Track Wt e3c", 200, 0 ,1 , 200, 0, 1);
    fOutput->Add(wt_match_e3c);
    
    qpt_tru1 = new TH1D("qpt_tru1","q/pt tru1",100, 0, 2);
    fOutput->Add(qpt_tru1);
    
    qpt_tru2 = new TH1D("qpt_tru2","q/pt tru2",100, 0, 2);
    fOutput->Add(qpt_tru2);
    
    pt_tru1 = new TH1D("pt_tru1","pt tru1",1000, 0, 200);
    fOutput->Add(pt_tru1);
    
    pt_tru2 = new TH1D("pt_tru2","pt tru2",1000, 0, 200);
    fOutput->Add(pt_tru2);
    
    //Track matching ENC histograms
    eec_Mm = new TH2D("EEC_Mm", "EEC match miss tracks", 100, new_bins, 21, 15, 120);
    eec_Mm->RebinX(2);
    fOutput->Add(eec_Mm);
    
    eec_mm = new TH2D("EEC_mm", "EEC miss miss tracks", 100, new_bins, 21, 15, 120);
    eec_mm->RebinX(2);
    fOutput->Add(eec_mm);
    
    e3c_MMm = new TH2D("E3C_MMm", "E3C match match miss tracks", 100, new_bins, 21, 15, 120);
    e3c_MMm->RebinX(2);
    fOutput->Add(e3c_MMm);
    
    e3c_Mmm = new TH2D("E3C_Mmm", "E3C match miss miss tracks", 100, new_bins, 21, 15, 120);
    e3c_Mmm->RebinX(2);
    fOutput->Add(e3c_Mmm);
    
    e3c_mmm = new TH2D("E3C_mmm", "E3C miss miss miss tracks", 100, new_bins, 21, 15, 120);
    e3c_mmm->RebinX(2);
    fOutput->Add(e3c_mmm);
  
    eec_Mf = new TH2D("EEC_Mf", "EEC match fake tracks", 100, new_bins, 21, 15, 120);
    eec_Mf->RebinX(2);
    fOutput->Add(eec_Mf);
    
    eec_ff = new TH2D("EEC_ff", "EEC fake fake tracks", 100, new_bins, 21, 15, 120);
    eec_ff->RebinX(2);
    fOutput->Add(eec_ff);
    
    e3c_MMf = new TH2D("E3C_MMf", "E3C match match fake tracks", 100, new_bins, 21, 15, 120);
    e3c_MMf->RebinX(2);
    fOutput->Add(e3c_MMf);
    
    e3c_Mff = new TH2D("E3C_Mff", "E3C match fake fake tracks", 100, new_bins, 21, 15, 120);
    e3c_Mff->RebinX(2);
    fOutput->Add(e3c_Mff);
    
    e3c_fff = new TH2D("E3C_fff", "E3C fake fake fake tracks", 100, new_bins, 21, 15, 120);
    e3c_fff->RebinX(2);
    fOutput->Add(e3c_fff);
    
    eec_matched_det = new TH2D("EEC_MM_det", "EEC match match det tracks", 100, new_bins, 21, 15, 120);
    eec_matched_det->RebinX(2);
    fOutput->Add(eec_matched_det);
  
    eec_matched_tru = new TH2D("EEC_MM_tru", "EEC match match tru tracks", 100, new_bins, 21, 15, 120);
    eec_matched_tru->RebinX(2);
    fOutput->Add(eec_matched_tru);
    
    e3c_matched_det = new TH2D("E3C_MMM_det", "E3C match match match det tracks", 100, new_bins, 21, 15, 120);
    e3c_matched_det->RebinX(2);
    fOutput->Add(e3c_matched_det);
  
    e3c_matched_tru = new TH2D("E3C_MMM_tru", "E3C match match match tru tracks", 100, new_bins, 21, 15, 120);
    e3c_matched_tru->RebinX(2);
    fOutput->Add(e3c_matched_tru);
  
    wt_res_eec = new TH2D("wt_res_eec", "Weight resolution scale EEC", 40, -1, 1, 200, 0, 1);
    fOutput->Add(wt_res_eec);
    
    wt_res_e3c = new TH2D("wt_res_e3c", "Weight resolution scale E3C", 40, -1, 1, 200, 0, 1);
    fOutput->Add(wt_res_e3c);
    
    R_res_eec = new TH2D("R_res_eec", "Weight resolution scale E3C", 40, -1, 1, 200, 0, 1);
    fOutput->Add(R_res_eec);
    
    R_res_e3c = new TH2D("R_res_e3c", "Weight resolution scale E3C", 40, -1, 1, 200, 0, 1);
    fOutput->Add(R_res_e3c);
    
    wtnojet_match_eec = new TH2D("wtnojet_match_eec", "Matched Track Wt excluding jet pT EEC", 200, 0 ,10000 , 200, 0, 10000);
    fOutput->Add(wtnojet_match_eec);
    
    wtnojet_match_e3c = new TH2D("wtnojet_match_e3c", "Matched Track Wt excluding jet pT E3C", 200, 0 ,100000 , 200, 0, 100000);
    fOutput->Add(wtnojet_match_e3c);
    
    wtnojet_res_eec = new TH2D("wtnojet_res_eec", "Weight resolution (excluding jet pT) scale EEC", 40, -1, 1, 200, 0, 10000);
    fOutput->Add(wtnojet_res_eec);
    
    wtnojet_res_e3c = new TH2D("wtnojet_res_e3c", "Weight resolution (excluding jet pT) scale E3C", 40, -1, 1, 200, 0, 10000);
    fOutput->Add(wtnojet_res_e3c);
    
    PostData(1, fOutput);


}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEEC::Run() {
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEEC::FillHistograms()
{
    if(fpairfastsim == 1){ComputeDelqpt(); MatchTracks();}//compute track matching and pair eff at track level
    
    //fill histogram goes through event loops
    AliEmcalJet *jet1 = NULL;
    //container zero is always the base container: the data container, the embedded subtracted in the case of embedding or the detector level in case of pythia.
    AliJetContainer *jetCont = GetJetContainer(0);// Get jet container.
    
    if (fCentSelectOn)
        if ((fCent > fCentMax) || (fCent < fCentMin))
            return 0;
    
    Float_t rhoVal = 0, rhoMassVal = 0.;
    if (jetCont)
    {
        jetCont->ResetCurrentID();
        if ((fJetShapeSub == kConstSub) || (fJetShapeSub == kDerivSub))
        {
            // rho
            AliRhoParameter *rhoParam = dynamic_cast<AliRhoParameter *>(
                                                                        InputEvent()->FindListObject("RhoSparseR020"));
            if (!rhoParam)
            {
                Printf("%s: Could not retrieve rho %s (some histograms will be filled "
                       "with zero)!",
                       GetName(), jetCont->GetRhoName().Data());
            }
            else
                rhoVal = rhoParam->GetVal();
            // rhom
            AliRhoParameter *rhomParam = dynamic_cast<AliRhoParameter *>(
                                                                         InputEvent()->FindListObject("RhoMassSparseR020"));
            if (!rhomParam)
            {
                Printf("%s: Could not retrieve rho_m %s (some histograms will be "
                       "filled with zero)!",
                       GetName(), jetCont->GetRhoMassName().Data());
            }
            else
                rhoMassVal = rhomParam->GetVal();
        }
        
        //Keeping triggers for full jet analysis (if desired)
        //     Int_t mytrig=-1;
        //     Bool_t mytrigmb=kFALSE;
        //     Bool_t mytrigej1=kFALSE;
        //     Bool_t mytrigej2=kFALSE;
        //     Double_t weightmb=-1;
        //     Double_t weightej1=-1;
        //     Double_t weightej2=-1;
        //
        //    if(fStoreTrig==kTRUE)
        //    {
        //    if(fInputHandler->IsEventSelected() & AliVEvent::kINT7 && fInputEvent->GetFiredTriggerClasses().Contains("INT7")) mytrigmb=kTRUE;
        //    if(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE && fInputEvent->GetFiredTriggerClasses().Contains("EJ1")) mytrigej1=kTRUE;
        //    if(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE && fInputEvent->GetFiredTriggerClasses().Contains("EJ2")) mytrigej2=kTRUE;
        //    if(mytrigmb==kTRUE && mytrigej1==kTRUE && mytrigej2==kTRUE) mytrig=3;
        //    if(mytrigmb==kTRUE && mytrigej1==kFALSE && mytrigej2==kFALSE) mytrig=0;
        //    if(mytrigmb==kFALSE && mytrigej1==kTRUE && mytrigej2==kFALSE) mytrig=1;
        //    if(mytrigmb==kFALSE && mytrigej1==kFALSE && mytrigej2==kTRUE) mytrig=2;
        //    if(mytrigmb==kTRUE && mytrigej1==kTRUE && mytrigej2==kFALSE) mytrig=4;
        //    if(mytrigmb==kTRUE && mytrigej1==kFALSE && mytrigej2==kTRUE) mytrig=5;
        //    if(mytrigmb==kFALSE && mytrigej1==kTRUE && mytrigej2==kTRUE) mytrig=6;
        //    if(mytrig==-1) return 0;
        //    }
        //
        //  UInt_t newrun=InputEvent()->GetRunNumber();
        //
        //if(fStoreTrig==kTRUE)
        //  {
        //    RunChanged(newrun);
        //      if(mytrig==3)
        //    {
        //      weightmb = 1./GetDownscaleWeight("INT7");
        //        weightej1 = 1./GetDownscaleWeight("EJ1");
        //      weightej2 = 1./GetDownscaleWeight("EJ2");
        //    }
        //      if(mytrig==0)
        //    {
        //   weightmb = 1./GetDownscaleWeight("INT7");
        //    }
        //      if(mytrig==1)
        //      {
        //          weightej1 = 1./GetDownscaleWeight("EJ1");
        //      }
        //      if(mytrig==2)
        //      {
        //          weightej2 = 1./GetDownscaleWeight("EJ2");
        //      }
        //      if(mytrig==4)
        //      {
        //      weightmb = 1./GetDownscaleWeight("INT7");
        //      weightej1 = 1./GetDownscaleWeight("EJ1");
        //      }
        //      if(mytrig==5)
        //      {
        //    weightmb = 1./GetDownscaleWeight("INT7");
        //    weightej2 = 1./GetDownscaleWeight("EJ2");
        //      }
        //
        //      if(mytrig==6)
        //      {
        //      weightej1 = 1./GetDownscaleWeight("EJ1");
        //      weightej2 = 1./GetDownscaleWeight("EJ2");
        //      }
        //  }
        
        //Jet Loop
        while ((jet1 = jetCont->GetNextAcceptJet()))
        {
            if (!jet1) continue;
            if (jet1->Pt() < 15) {continue;} //Cuts on jet_pt
            if (fJetShapeType == kData)
            {
                ComputeEEC(jet1, jetCont);//Computing the eec on the jet object
            }
            
            AliEmcalJet *jet2 = 0x0;
            AliEmcalJet *jet3 = 0x0;
            
            AliEmcalJet *jetUS = NULL;
            Int_t ifound = 0, jfound = 0;
            Int_t ilab = -1, jlab = -1;
            
            // All for MC. This is the mode to run over pythia to produce a det-part response.
            //Here we have also added the constituent-subtraction case, but we don't use it normally in pp the matching is purely geometrical
            AliJetContainer *jetContTrue;
            AliJetContainer *jetContUS;
            AliJetContainer *jetContPart;
            
            if (fJetShapeType == kPythiaDef)
            {
                jetContTrue = GetJetContainer(1); //we never use this in pp, usually an empty one
                jetContUS = GetJetContainer(2); //unsubtracted one, don't need for pp
                jetContPart = GetJetContainer(3); //this is pythia particle/true level
                
                if (fJetShapeSub == kConstSub)
                {
                    for (Int_t i = 0; i < jetContUS->GetNJets(); i++)
                    {
                        jetUS = jetContUS->GetJet(i);
                        if (jetUS->GetLabel() == jet1->GetLabel())
                        {
                            ifound++;
                            if (ifound == 1) ilab = i;
                        }
                        
                    }
                    if (ilab == -1) continue;
                    jetUS = jetContUS->GetJet(ilab);
                    jet2 = jetUS->ClosestJet();
                    if (!jet2)
                    {
                        Printf("jet2 does not exist, returning");
                        continue;
                        
                    }
                    for (Int_t j = 0; j < jetContPart->GetNJets(); j++)
                    {
                        jet3 = jetContPart->GetJet(j);
                        if (!jet3) continue;
                        if (jet3->GetLabel() == jet2->GetLabel())
                        {
                            jfound++;
                            if (jfound == 1) jlab = j;
                        }
                        
                    }
                    if (jlab == -1) continue;
                    jet3 = jetContPart->GetJet(jlab);
                    if (!jet3)
                    {
                        Printf("jet3 does not exist, returning");
                        continue;
                    }
                    
                }
                if (!(fJetShapeSub == kConstSub))
                    jet3 = jet1->ClosestJet();
                if (!jet3)
                {
                    //   Printf("jet3 does not exist, returning");
                    continue;
                }
            }
            
            //These are options for pb-pb analysis
            Double_t ptSubtracted = 0;
            if (fJetShapeSub == kConstSub || fJetShapeSub == kEventSub) ptSubtracted = jet1->Pt();
            else if (fJetShapeSub == kDerivSub)
            {
                ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
            }
            else if (fJetShapeSub == kNoSub) ptSubtracted = jet1->Pt();
            if (ptSubtracted < fPtThreshold) continue;
            if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1)) continue;
            
            
            Double_t ptMatch=0;
            Int_t kMatched = 0;
            if (fJetShapeType == kPythiaDef)
            {
                kMatched = 1;
                if (fJetShapeSub == kConstSub)
                    kMatched = 3;
                ptMatch = jet3->Pt();
                //                cout<<"the matched jet "<<jet3->Pt()<<" "<<kMatched<<endl;
                ComputeEncMC(jet1, jetCont, jet3, kMatched);
            }
        } //close while loop
    } //close the jet cont loop
    
    return kTRUE;
}



//________________________________________________________________________
int AliAnalysisTaskJetsEEC::GetConstituentID(int constituentIndex, const AliVParticle* part, AliEmcalJet * jet)
{
  // NOTE: Usually, we would use the global offset defined for the general subtracter extraction task. But we don't want to
  //       depend on that task, so we just define it here locally.
  //Get the id of the particle. If the label is not equal to -1, get the label and assign it to id. If it is -1, get the value of of track at id + 20000
  int id = part->GetLabel() != -1 ? part->GetLabel() : (jet->TrackAt(constituentIndex) + 2000000);
  return id;
}


//_______________________________________________________________________
Double_t AliAnalysisTaskJetsEEC::GetDownscaleWeight(string trigString)
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
void AliAnalysisTaskJetsEEC::MatchTracks()
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


//______________________________________________________________________
void AliAnalysisTaskJetsEEC::ComputeDelqpt()
{
    AliVTrack* track1;
    AliVTrack* track2;
    AliAODTrack *aodtrack1;
    AliAODTrack *aodtrack2;
    Short_t charge_det1;
    Short_t charge_det2;
    double delqpt_det;
    double dphi_det;
    double dphi_det_sq;
    double deta_det_sq;
    double delR_det;
    Double_t pt_det1;
    Double_t pt_det2;
    Double_t phi_det1;
    Double_t phi_det2;
    Double_t eta_det1;
    Double_t eta_det2;
    Byte_t type1;
    Byte_t type2;
 
    auto iterable = fDetectorLevel->accepted_momentum(); //creates an iterable container interface over accepted objects in the EMcal cont
    for (auto trackIterator = iterable.begin(); trackIterator != iterable.end(); trackIterator++) {
        track1 = trackIterator->second; //why second? Because trackIterator is of the form std::pair<AliTLorentzVector, AliVTrack *>'
        type1 = fDetectorLevel->GetTrackType(track1);
        if (type1 <= 2) {
            if (fIsEsd) {
                continue;
            } else { // AOD
                aodtrack1 = dynamic_cast<AliAODTrack*>(track1);
                if(!track1) AliFatal("Not a standard AOD");
                charge_det1 = track1->Charge();
                pt_det1 = track1->Pt();
                phi_det1 = track1->Phi();
                eta_det1 = track1->Eta();
                track_pt_det->Fill(pt_det1);
            }
            for (auto trackIterator = iterable.begin(); trackIterator != iterable.end(); trackIterator++) {
                track2 = trackIterator->second;
                type2 = fDetectorLevel->GetTrackType(track2);
                if (type2 <= 2) {
                    if (fIsEsd) {
                        continue;
                    } else { // AOD
                        aodtrack2 = dynamic_cast<AliAODTrack*>(track2);
                        if(!track2) AliFatal("Not a standard AOD");
                        charge_det2 = track2->Charge();
                        pt_det2 = track2->Pt();
                        phi_det2 = track2->Phi();
                        eta_det2 = track2->Eta();
                        delqpt_det = TMath::Abs(charge_det2/pt_det2 - charge_det1/pt_det1);
                        dphi_det = (phi_det2-phi_det1);
                        if (dphi_det>TMath::Pi()) {dphi_det = dphi_det - 2.*TMath::Pi();}
                        if (dphi_det<(-1.*TMath::Pi())) {dphi_det = dphi_det + 2.*TMath::Pi();}
                        dphi_det_sq = dphi_det*dphi_det;
//                        deta_det_sq = (eta_det2-eta_det1)*(phi_det2-phi_det1); // This was wrong 
                        deta_det_sq = (eta_det2-eta_det1)*(eta_det2-eta_det1);
                        delR_det = TMath::Sqrt(dphi_det_sq + deta_det_sq);
                        qpt_det->Fill(delR_det,delqpt_det);
                    }
                }
            }
        }
    }
    
    AliAODMCParticle* part1;
    AliAODMCParticle* part2;
    Short_t charge_part1;
    Short_t charge_part2;
    double delqpt_part;
    double dphi_part;
    double dphi_part_sq;
    double deta_part_sq;
    double delR_part;
    Double_t pt_part1;
    Double_t phi_part1;
    Double_t eta_part1;
    Double_t pt_part2;
    Double_t phi_part2;
    Double_t eta_part2;
    
    if (fGeneratorLevel) {
        auto iterablepart = fGeneratorLevel->accepted_momentum(); //creates an iterable container interface over accepted objects in the EMcal cont
        for (auto partIterator = iterablepart.begin(); partIterator != iterablepart.end(); partIterator++) {
            part1 = partIterator->second;
            Short_t charge_part1 = part1->Charge();
            pt_part1 = part1->Pt();
//            if(pt_part1 < 0.150) continue;
            phi_part1 = part1->Phi();
            eta_part1 = part1->Eta();
            track_pt_tru->Fill(pt_part1);
            for (auto partIterator = iterablepart.begin(); partIterator != iterablepart.end(); partIterator++) {
                part2 = partIterator->second;
                Short_t charge_part2 = part2->Charge();
                pt_part2 = part2->Pt();
//                if(pt_part2 < 0.150) continue;
                phi_part2 = part2->Phi();
                eta_part2 = part2->Eta();
                dphi_part = (phi_part2-phi_part1);
                if (dphi_part>TMath::Pi()) {dphi_part =dphi_part - 2.*TMath::Pi();}
                if (dphi_part<(-1.*TMath::Pi())) {dphi_part =dphi_part + 2.*TMath::Pi();}
                dphi_part_sq = dphi_part*dphi_part;
                deta_part_sq = (eta_part2-eta_part1)*(eta_part2-eta_part1);
                delqpt_part = TMath::Abs(charge_part2/pt_part2 - charge_part1/pt_part1);
                delR_part = TMath::Sqrt(dphi_part_sq + deta_part_sq);
//                pt_tru1->Fill(pt_part1);
//                pt_tru2->Fill(pt_part2);
//                qpt_tru2->Fill(charge_part2/pt_part2);
//                qpt_tru1->Fill(charge_part1/pt_part1);
                qpt_tru->Fill(delR_part,delqpt_part);
            }
        }
    }
}


//______________________________________________________________________
void AliAnalysisTaskJetsEEC::ComputeEncMC(AliEmcalJet *fJet, AliJetContainer *fJetCont, AliEmcalJet *fJet_tru, Int_t km)
    //(jet, jet container, vector of jets, vector of constituents within each jet?)
{
    //Det level
    std::vector<fastjet::PseudoJet> fConstituents; //Is a pseudojet object with constituents of the jet
    fConstituents.clear();
    //This snippet of code is getting particles within a single jet (fjet) and turning them into pseudojet objects so that fastjet capabilities can be used
    fastjet::PseudoJet PseudoTracks; //Creating a pseudojet object called PseduoTracks
    unsigned int constituentIndex = 0;
    int constituentCharge_det = 0;
    //The line below gets constituent particles within fjet. C++ syntax[ for (auto elem : container)    // capture elements by value ]
    for (auto part: fJet->GetParticleConstituents())
    {
        PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E()); //part is the constituent at that point in the loop, part keeps getting redefined in each step.
        const AliVParticle* part2 = part.GetParticle(); //"hack", leave this in , to get the index of the jet from AliPhysics
        PseudoTracks.set_user_index(GetConstituentID(constituentIndex, part2, fJet)); //leave this in for the same reason as above
//        PseudoTracks.set_user_index(GetConstituentCharge(constituentCharge_det, part2, fJet_tru));
        if (PseudoTracks.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
        fConstituents.push_back(PseudoTracks);
        constituentIndex++;
    }

    //Truth level
    AliJetContainer *jetCont = GetJetContainer(km); //get the container for the matched true level jet
    std::vector<fastjet::PseudoJet> fConstituents_tru; //Is a pseudojet object with constituents of the jet
    fConstituents_tru.clear();
    //This snippet of code is getting particles within a single jet (fjet) and turning them into pseudojet objects so that fastjet capabilities can be used
    fastjet::PseudoJet PseudoTracks_tru; //Creating a pseudojet object called PseduoTracks
    unsigned int constituentIndex_tru = 0;
    int constituentCharge_tru = 0;
    for (auto part_tru: fJet_tru->GetParticleConstituents())
    {
        PseudoTracks_tru.reset(part_tru.Px(), part_tru.Py(), part_tru.Pz(), part_tru.E()); //part is the constituent at that point in the loop, part keeps getting redefined in each step.
        const AliVParticle* part_tru2 = part_tru.GetParticle(); //"hack", leave this in , to get the index of the jet from AliPhysics
        PseudoTracks_tru.set_user_index(GetConstituentID(constituentIndex_tru, part_tru2, fJet_tru)); //leave this in for the same reason as above
//        test_hist->Fill(PseudoTracks_tru.user_index());
//        PseudoTracks_tru.set_user_index(GetConstituentCharge(constituentCharge_tru, part_tru2, fJet_tru)); //leave this in for the same reason as above
        if (PseudoTracks_tru.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
        fConstituents_tru.push_back(PseudoTracks_tru);
        constituentIndex_tru++;
    }

    double jet_pt = fJet->Pt();
    double jet_pt_tru = fJet_tru->Pt();
    pt_tru->Fill(jet_pt_tru); //filling histogram with momentum of jets
    
    
    std::vector<Double_t> R_dist_tru, R_dist_det;
    
    if(fMatchJetTrack==1)
    {
        for(int j=0; j<int(fConstituents.size()); j++)//det level loop
        {
            int j_index = fConstituents[j].user_index();
            
            for(int i=0; i<int(fConstituents_tru.size()); i++)//tru level loop
            {
                int i_index = fConstituents_tru[i].user_index();
                
                if(i_index == j_index) //does the det particle have a matching truth particle
                {
                    for(int s=0; s<j; s++)//det level loop
                    {
                        int s_index = fConstituents[s].user_index();
                        
                        for(int l=0; l<i; l++)//tru level loop
                        {
                            int l_index = fConstituents_tru[l].user_index();
                            if(l_index == s_index && l_index != j_index) //does the det particle have a matching truth particle & is it a unique match
                            {
                                double ee_det = 2*fConstituents[j].pt()*fConstituents[s].pt()/(pow(jet_pt,2));
                                double ee_tru = 2*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()/(pow(jet_pt_tru,2));
                                double R_det = fConstituents[j].delta_R(fConstituents[s]);
                                double R_tru = fConstituents_tru[i].delta_R(fConstituents_tru[l]);
                                
                                double ee_jss_det = 3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt()/(pow(jet_pt,3));
                                double ee_jss_tru = 3*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()*fConstituents_tru[l].pt()/(pow(jet_pt_tru,3));
                                double R_jss_det = fConstituents[j].delta_R(fConstituents[s]);
                                double R_jss_tru = fConstituents_tru[i].delta_R(fConstituents_tru[l]);
                                
                                eec_matched_det->Fill(R_det,jet_pt,ee_det);
                                eec_matched_tru->Fill(R_tru,jet_pt_tru,ee_tru);
                                
                                R_match_eec->Fill(R_tru,R_det);
                                R_res_eec->Fill(R_tru,(R_tru-R_det)/(R_tru));
                                wt_match_eec->Fill(ee_tru,ee_det);
                                wtnojet_match_eec->Fill(2*fConstituents_tru[i].pt()*fConstituents_tru[l].pt(),2*fConstituents[j].pt()*fConstituents[s].pt());
                                
                                double wt_diff_eec = (ee_tru-ee_det)/ee_tru;
                                wt_res_eec->Fill(ee_tru,wt_diff_eec);
                                
                                double wt_diff_eec_nojet = ((2*fConstituents_tru[i].pt()*fConstituents_tru[l].pt())-(2*fConstituents[j].pt()*fConstituents[s].pt()))/(2*fConstituents_tru[i].pt()*fConstituents_tru[l].pt());
                                wtnojet_res_eec->Fill(wt_diff_eec_nojet,2*fConstituents_tru[i].pt()*fConstituents_tru[l].pt());
                                
                                
                                e3c_matched_det->Fill(R_jss_det,jet_pt,ee_jss_det);
                                e3c_matched_tru->Fill(R_jss_tru,jet_pt_tru,ee_jss_tru);
                                
                                R_match_e3c->Fill(R_jss_tru,R_jss_det);
                                R_res_e3c->Fill(R_jss_tru,(R_jss_tru-R_jss_det)/(R_jss_tru));
                                wt_match_e3c->Fill(ee_jss_tru,ee_jss_det);
                                
                                wtnojet_match_e3c->Fill(3*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()*fConstituents_tru[l].pt(),3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt());
                                
                                double wt_diff_e3c_jss = (ee_jss_tru-ee_jss_det)/ee_jss_tru;
                                wt_res_e3c->Fill(wt_diff_e3c_jss,ee_jss_tru);
                                
                                double wt_diff_e3c_nojet_jss = ((3*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()*fConstituents_tru[l].pt())-(3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt()))/(3*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()*fConstituents_tru[l].pt());//diff/true
                                wtnojet_res_e3c->Fill(wt_diff_e3c_nojet_jss, 3*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()*fConstituents_tru[l].pt());
                                
                                for(int m=0; m!=j && m!=s; m++)
                                {
                                    if(s>j) continue;
                                    int m_index = fConstituents[m].user_index();
                                    
                                    for(int k=0; k!=l && k!=i; k++)
                                    {
                                        int k_index = fConstituents_tru[k].user_index();
                                        if(k_index == m_index && k_index != s_index && k_index != j_index )
                                        {
                                            double ee_jsm_det = 6*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[m].pt()/(pow(jet_pt,3));
                                            double ee_jsm_tru = 6*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()*fConstituents_tru[k].pt()/(pow(jet_pt_tru,3));
                                            
                                            double deltaR_js_det = fConstituents[j].delta_R(fConstituents[s]);
                                            double deltaR_jm_det = fConstituents[j].delta_R(fConstituents[m]);
                                            double deltaR_sm_det = fConstituents[s].delta_R(fConstituents[m]);
                                            R_dist_det.push_back(deltaR_js_det);
                                            R_dist_det.push_back(deltaR_jm_det);
                                            R_dist_det.push_back(deltaR_sm_det);
                                            int max_R_det = distance(R_dist_det.begin(), max_element(R_dist_det.begin(), R_dist_det.end()));
                                            
                                            
                                            double deltaR_js_tru = fConstituents_tru[i].delta_R(fConstituents_tru[l]);
                                            double deltaR_jm_tru = fConstituents_tru[i].delta_R(fConstituents_tru[k]);
                                            double deltaR_sm_tru = fConstituents_tru[l].delta_R(fConstituents_tru[k]);
                                            R_dist_tru.push_back(deltaR_js_tru);
                                            R_dist_tru.push_back(deltaR_jm_tru);
                                            R_dist_tru.push_back(deltaR_sm_tru);
                                            int max_R_tru = distance(R_dist_tru.begin(), max_element(R_dist_tru.begin(), R_dist_tru.end()));
                                            
                                            e3c_matched_det->Fill(R_dist_det[max_R_det],jet_pt,ee_jsm_det);
                                            e3c_matched_tru->Fill(R_dist_tru[max_R_tru],jet_pt_tru,ee_jsm_tru);
                                            
                                            R_match_e3c->Fill(R_dist_tru[max_R_tru],R_dist_det[max_R_det]);
                                            R_res_e3c->Fill(R_dist_tru[max_R_tru],(R_dist_tru[max_R_tru]-R_dist_det[max_R_det])/(R_dist_tru[max_R_tru]));
                                            wt_match_e3c->Fill(ee_jsm_tru,ee_jsm_det);
                                            
                                            wtnojet_match_e3c->Fill(6*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()*fConstituents_tru[k].pt(),6*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[m].pt());
                                            
                                            double wt_diff_e3c_jsm = (ee_jsm_tru-ee_jsm_det)/ee_jsm_tru;
                                            wt_res_e3c->Fill(wt_diff_e3c_jsm,ee_jsm_tru);
                                            
                                            double wt_diff_e3c_nojet_jsm = ((6*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()*fConstituents_tru[k].pt())-(6*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[m].pt()))/(6*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()*fConstituents_tru[k].pt());//diff/true
                                            wtnojet_res_e3c->Fill(wt_diff_e3c_nojet_jsm, 6*fConstituents_tru[i].pt()*fConstituents_tru[l].pt()*fConstituents_tru[k].pt());
                                            
                                            R_dist_det.clear();
                                            R_dist_tru.clear();
                                            
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    std::vector<int> det_index,tru_index;
    std::vector<fastjet::PseudoJet> matchtracks_det, matchtracks_tru, misstracks ,faketracks;
    
    for(int j=0; j<int(fConstituents.size()); j++)
    {
        det_index.push_back(fConstituents[j].user_index());
    }
    for(int i=0; i<int(fConstituents_tru.size()); i++)
    {
        tru_index.push_back(fConstituents_tru[i].user_index());
    }
    
    for(int j = 0;j<int(fConstituents.size()); j++)
    {
        int valueToCheck = fConstituents[j].user_index();
        auto it = std::find(tru_index.begin(), tru_index.end(), valueToCheck);
        if (it != tru_index.end()){matchtracks_det.push_back(fConstituents[j]);}
        else {faketracks.push_back(fConstituents[j]);}
    }
    
    for(int i=0; i<int(fConstituents_tru.size()); i++)
    {
        int valueToCheck = fConstituents_tru[i].user_index();
        auto it = std::find(det_index.begin(), det_index.end(), valueToCheck);
        if (it != det_index.end()){matchtracks_tru.push_back(fConstituents_tru[i]);}
        else {misstracks.push_back(fConstituents_tru[i]);}
    }
   
    if(fMissJetTrack==1) //compute (match,miss) & (miss,miss) for eec and (match,match,miss), (match,miss,miss) & (miss,miss,miss) for the E3C
    {
        for (int j = 0; j < int(matchtracks_tru.size()); j++) // match
        {
            for (int s = 0; s < int(misstracks.size()); s++) // miss
            {
                double ee_js_tru = matchtracks_tru[j].pt() * misstracks[s].pt() / (pow(jet_pt_tru, 2));
                double R_js_tru = matchtracks_tru[j].delta_R(misstracks[s]);
                
                // Fill histograms for match, miss
                eec_Mm->Fill(R_js_tru, jet_pt_tru, ee_js_tru);
                
                for (int m = 0; m < int(misstracks.size()); m++) // miss
                {
                    if (j != 0)
                        break;
                    
                    if (m == s)
                        continue;
                    
                    double ee_jm_tru = matchtracks_tru[j].pt() * misstracks[m].pt() / (pow(jet_pt_tru, 2));
                    double R_jm_tru = matchtracks_tru[j].delta_R(misstracks[m]);
                    
                    // Fill histograms for miss, miss
                    eec_mm->Fill(R_jm_tru, jet_pt_tru, ee_jm_tru);
                    
                    double ee_jsm_tru = matchtracks_tru[j].pt() * misstracks[s].pt() * misstracks[m].pt() / (pow(jet_pt_tru, 3));
                    double deltaR_js_tru = matchtracks_tru[j].delta_R(misstracks[s]);
                    double deltaR_jm_tru = matchtracks_tru[j].delta_R(misstracks[m]);
                    double deltaR_sm_tru = misstracks[s].delta_R(misstracks[m]);
                    
                    std::vector<double> R_dist_tru_Mmm = {deltaR_js_tru, deltaR_jm_tru, deltaR_sm_tru};
                    int max_R_index_tru = std::distance(R_dist_tru_Mmm.begin(), std::max_element(R_dist_tru_Mmm.begin(), R_dist_tru_Mmm.end()));
                    
                    // Fill histograms for match miss miss
                    e3c_Mmm->Fill(R_dist_tru_Mmm[max_R_index_tru], jet_pt_tru, ee_jsm_tru);
                    
                    for (int n = 0; n < int(misstracks.size()); n++) // miss
                    {
                        if (j != 0)
                            break;
                        
                        if (n == m || n == s)
                            continue;
                        
                        double ee_snm_tru = misstracks[n].pt() * misstracks[s].pt() * misstracks[m].pt() / (pow(jet_pt_tru, 3));
                        double deltaR_nm_tru = misstracks[n].delta_R(misstracks[m]);
                        double deltaR_sn_tru = misstracks[s].delta_R(misstracks[n]);
                        double deltaR_sm_tru = misstracks[s].delta_R(misstracks[m]);
                        
                        std::vector<double> R_dist_tru_mmm = {deltaR_nm_tru, deltaR_sn_tru, deltaR_sm_tru};
                        int max_R_index_tru = std::distance(R_dist_tru_mmm.begin(), std::max_element(R_dist_tru_mmm.begin(), R_dist_tru_mmm.end()));
                        
                        // Fill histograms for miss miss miss
                        e3c_mmm->Fill(R_dist_tru_mmm[max_R_index_tru], jet_pt_tru, ee_snm_tru);
                    }
                }
            }
            
            for(int k=0; k<int(matchtracks_tru.size()); k++) //match match miss
            {
                if(k==j) continue;
                for(int l=0; l<int(misstracks.size()); l++) //miss
                {
                    double ee_jkl_tru = matchtracks_tru[j].pt()*matchtracks_tru[k].pt()*misstracks[l].pt()/(pow(jet_pt_tru,3));
                    double deltaR_jk_tru = matchtracks_tru[j].delta_R(matchtracks_tru[k]);
                    double deltaR_jl_tru = matchtracks_tru[j].delta_R(misstracks[l]);
                    double deltaR_kl_tru = matchtracks_tru[k].delta_R(misstracks[l]);
                    std::vector<double> R_dist_tru_MMm = {deltaR_jk_tru,deltaR_jl_tru,deltaR_kl_tru};
                    
                    int max_R_tru = distance(R_dist_tru_MMm.begin(), max_element(R_dist_tru_MMm.begin(), R_dist_tru_MMm.end()));
                    
                    e3c_MMm->Fill(R_dist_tru_MMm[max_R_tru],jet_pt_tru,ee_jkl_tru);
                    
                    R_dist_tru_MMm.clear();
                }
            }
        }
    }
    
    if(fFakeJetTrack==1)
    {
        for (int j = 0; j < int(matchtracks_tru.size()); j++) // match
        {
            for (int s = 0; s < int(faketracks.size()); s++) // fake
            {
                double ee_js_tru = matchtracks_tru[j].pt() * faketracks[s].pt() / (pow(jet_pt_tru, 2));
                double R_js_tru = matchtracks_tru[j].delta_R(faketracks[s]);
                
                // Fill histograms for match, fake
                eec_Mf->Fill(R_js_tru, jet_pt_tru, ee_js_tru);
                
                for (int m = 0; m < int(faketracks.size()); m++) // fake
                {
                    if (j != 0)
                        break;
                    
                    if (m == s)
                        continue;
                    
                    double ee_jm_tru = matchtracks_tru[j].pt() * faketracks[m].pt() / (pow(jet_pt_tru, 2));
                    double R_jm_tru = matchtracks_tru[j].delta_R(faketracks[m]);
                    
                    // Fill histograms for fake fake
                    eec_ff->Fill(R_jm_tru, jet_pt_tru, ee_jm_tru);
                    
                    double ee_jsm_tru = matchtracks_tru[j].pt() * faketracks[s].pt() * faketracks[m].pt() / (pow(jet_pt_tru, 3));
                    double deltaR_js_tru = matchtracks_tru[j].delta_R(faketracks[s]);
                    double deltaR_jm_tru = matchtracks_tru[j].delta_R(faketracks[m]);
                    double deltaR_sm_tru = faketracks[s].delta_R(faketracks[m]);
                    
                    std::vector<double> R_dist_tru_Mff = {deltaR_js_tru, deltaR_jm_tru, deltaR_sm_tru};
                    int max_R_index_tru = std::distance(R_dist_tru_Mff.begin(), std::max_element(R_dist_tru_Mff.begin(), R_dist_tru_Mff.end()));
                    
                    // Fill histograms for match fake fake
                    e3c_Mff->Fill(R_dist_tru_Mff[max_R_index_tru], jet_pt_tru, ee_jsm_tru);
                    
                    for (int n = 0; n < int(faketracks.size()); n++) // fake
                    {
                        if (j != 0)
                            break;
                        
                        if (n == m || n == s)
                            continue;
                        
                        double ee_snm_tru = faketracks[n].pt() * faketracks[s].pt() * faketracks[m].pt() / (pow(jet_pt_tru, 3));
                        double deltaR_nm_tru = faketracks[n].delta_R(faketracks[m]);
                        double deltaR_sn_tru = faketracks[s].delta_R(faketracks[n]);
                        double deltaR_sm_tru = faketracks[s].delta_R(faketracks[m]);
                        
                        std::vector<double> R_dist_tru_fff = {deltaR_nm_tru, deltaR_sn_tru, deltaR_sm_tru};
                        int max_R_index_tru = std::distance(R_dist_tru_fff.begin(), std::max_element(R_dist_tru_fff.begin(), R_dist_tru_fff.end()));
                        
                        // Fill histograms for fake fake fake
                        e3c_fff->Fill(R_dist_tru_fff[max_R_index_tru], jet_pt_tru, ee_snm_tru);
                    }
                }
            }
            
            for(int k=0; k<int(matchtracks_tru.size()); k++) //match match fake
            {
                if(k==j) continue;
                for(int l=0; l<int(faketracks.size()); l++) //fake
                {
                    double ee_jkl_tru = matchtracks_tru[j].pt()*matchtracks_tru[k].pt()*faketracks[l].pt()/(pow(jet_pt_tru,3));
                    double deltaR_jk_tru = matchtracks_tru[j].delta_R(matchtracks_tru[k]);
                    double deltaR_jl_tru = matchtracks_tru[j].delta_R(faketracks[l]);
                    double deltaR_kl_tru = matchtracks_tru[k].delta_R(faketracks[l]);
                    std::vector<double> R_dist_tru_MMf = {deltaR_jk_tru,deltaR_jl_tru,deltaR_kl_tru};
                    
                    int max_R_tru = distance(R_dist_tru_MMf.begin(), max_element(R_dist_tru_MMf.begin(), R_dist_tru_MMf.end()));
                    
                    e3c_MMf->Fill(R_dist_tru_MMf[max_R_tru],jet_pt_tru,ee_jkl_tru);
                    
                    R_dist_tru_MMf.clear();
                }
            }
        }
    }
    
    matchtracks_tru.clear();
    matchtracks_det.clear();
    faketracks.clear();
    misstracks.clear();

    if(fpTcorr == 1)
    {   jet_pt = jet_pt/(0.85); //applying JES correction to jet pT spectra to study spectra shape dependence of ENC
        jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
        
        double diff_scaled = (jet_pt - jet_pt_tru)/jet_pt_tru;
        JES_scaled->Fill(jet_pt_tru,diff_scaled);
        
        R_matrix->Fill(jet_pt_tru,jet_pt); //Filling the response matrix}
    }
    
    else
    {
        jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
        R_matrix->Fill(jet_pt_tru,jet_pt); //Filling the response matrix
        
        double diff = (jet_pt - jet_pt_tru)/jet_pt_tru;
        JES->Fill(jet_pt_tru,diff);
    }
    
    //Det level
    //Initializing objects
    std::vector<Double_t> R_dist;
    
    //Truth level
    //Initializing objects
    std::vector<Double_t> R_dist_part;

    
    //Looping over the det jet
    //For jets with 2 constituents
    if(int(fConstituents.size()) == 2)
    {
        for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
        {
            //For 3 point correlator
            for(int s=0; s<int(fConstituents.size()) ; s++)
            {
                if(s==j) continue; //if s=j this would be 0
                double eee_jss_2 =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                double deltaR_jss_2 = fConstituents[j].delta_R(fConstituents[s]);
                
                if(fpaircut == 1)
                {
                    double j_eta_3pt_2_det = fConstituents[j].eta();
                    double s_eta_3pt_2_det = fConstituents[s].eta();
                    double del_js_eta_3pt_2_det = abs(j_eta_3pt_2_det-s_eta_3pt_2_det);
                    if (del_js_eta_3pt_2_det < 0.008) continue;
                    else
                    {
                        
                        E3C_hist->Fill(deltaR_jss_2,eee_jss_2);
                        E3C_pt_hist->Fill(deltaR_jss_2,jet_pt,eee_jss_2);
                       
                        N3_det_pt_hist_3d->Fill(deltaR_jss_2, jet_pt, jet_pt_tru);
                        E3C_det_pt_hist_3d->Fill(deltaR_jss_2, jet_pt, jet_pt_tru, eee_jss_2);
                       
                    }
                }
                else
                {
                    E3C_hist->Fill(deltaR_jss_2,eee_jss_2);
                    E3C_pt_hist->Fill(deltaR_jss_2,jet_pt,eee_jss_2);
                  
                    N3_det_pt_hist_3d->Fill(deltaR_jss_2, jet_pt, jet_pt_tru);
                    E3C_det_pt_hist_3d->Fill(deltaR_jss_2, jet_pt, jet_pt_tru, eee_jss_2);
                }
                
            }//close s loop for the 3 point correlator
            //For 2 point correlator
            for(int s=0; s<j ; s++)
            {
                double delta_R_js_2 = fConstituents[j].delta_R(fConstituents[s]);
                double ee_js_2 = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));
                
                //Pair cut
                if(fpaircut == 1)
                {
                    double j_eta_2pt_2_det = fConstituents[j].eta();
                    double s_eta_2pt_2_det = fConstituents[s].eta();
                    double del_js_eta_2pt_2_det = abs(j_eta_2pt_2_det-s_eta_2pt_2_det);
                    if (del_js_eta_2pt_2_det < 0.008) continue;
                    else
                    {
                        
                        EEC_hist->Fill(delta_R_js_2,ee_js_2);
                        EEC_pt_hist->Fill(delta_R_js_2,jet_pt,ee_js_2);
                      
                        N2_det_pt_hist_3d->Fill(delta_R_js_2, jet_pt, jet_pt_tru);
                        EEC_det_pt_hist_3d->Fill(delta_R_js_2, jet_pt, jet_pt_tru, ee_js_2);
               
                    }
                }
                else
                {
                    
                    EEC_hist->Fill(delta_R_js_2,ee_js_2);
                    EEC_pt_hist->Fill(delta_R_js_2,jet_pt,ee_js_2);
                    
                    N2_det_pt_hist_3d->Fill(delta_R_js_2, jet_pt, jet_pt_tru);
                    EEC_det_pt_hist_3d->Fill(delta_R_js_2, jet_pt, jet_pt_tru, ee_js_2);
                }
            }//close s loop for eec
        }//close j loop
    }//close if loop
    
    //For jets with more than 2 constituents
    else
    {
        for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
        {
            
            for(int s=0; s<int(fConstituents.size()) ; s++)
            {
                if(s==j) continue; //This ensures I don't get stuff like (000) for (jss)
                
                double eee_jss =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                double deltaR_jss = fConstituents[j].delta_R(fConstituents[s]);
                
                //Pair cut
                if(fpaircut == 1)
                {
                    double j_eta_3pt_det = fConstituents[j].eta();
                    double s_eta_3pt_det = fConstituents[s].eta();
                    double del_js_eta_3pt_det = abs(j_eta_3pt_det-s_eta_3pt_det);
                    if (del_js_eta_3pt_det < 0.008) continue;
                    else
                    {
                        E3C_hist->Fill(deltaR_jss,eee_jss);
                        E3C_pt_hist->Fill(deltaR_jss,jet_pt,eee_jss);
                        
                        N3_det_pt_hist_3d->Fill(deltaR_jss, jet_pt, jet_pt_tru);
                        E3C_det_pt_hist_3d->Fill(deltaR_jss, jet_pt, jet_pt_tru, eee_jss);
                    }
                }
                else
                {
                    E3C_hist->Fill(deltaR_jss,eee_jss);
                    E3C_pt_hist->Fill(deltaR_jss,jet_pt,eee_jss);
                    
                    N3_det_pt_hist_3d->Fill(deltaR_jss, jet_pt, jet_pt_tru);
                    E3C_det_pt_hist_3d->Fill(deltaR_jss, jet_pt, jet_pt_tru, eee_jss);
                }
                //For 3 point correlator
                for( int m=0; m!=j && m!=s; m++)
                {
                    if(s>j) continue;
                    
                    double eee_jsm = ((6*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[m].pt())/(pow(jet_pt,3)));
                    double deltaR_js = fConstituents[j].delta_R(fConstituents[s]);
                   
                    double deltaR_jm = fConstituents[j].delta_R(fConstituents[m]);
                  
                    double deltaR_sm = fConstituents[s].delta_R(fConstituents[m]);
                  
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_3pt_det = fConstituents[j].eta();
                        double m_eta_3pt_det = fConstituents[m].eta();
                        double s_eta_3pt_det = fConstituents[s].eta();
                        double del_jm_eta_3pt_det = abs(j_eta_3pt_det-m_eta_3pt_det);
                        double del_sm_eta_3pt_det = abs(s_eta_3pt_det-m_eta_3pt_det);
                        double del_js_eta_3pt_det = abs(j_eta_3pt_det-s_eta_3pt_det);
                        if (del_jm_eta_3pt_det < 0.008 || del_sm_eta_3pt_det < 0.008 || del_js_eta_3pt_det < 0.008) continue;
                        else
                        {
                        
                            R_dist.push_back(deltaR_js);
                            R_dist.push_back(deltaR_jm);
                            R_dist.push_back(deltaR_sm);
                            int max_R = distance(R_dist.begin(), max_element(R_dist.begin(), R_dist.end()));//pick the longest side to compute the correlators with

                            E3C_hist->Fill(R_dist[max_R],eee_jsm);
                            E3C_pt_hist->Fill(R_dist[max_R],jet_pt,eee_jsm);
                            
                            N3_det_pt_hist_3d->Fill(R_dist[max_R], jet_pt, jet_pt_tru);
                            E3C_det_pt_hist_3d->Fill(R_dist[max_R], jet_pt, jet_pt_tru, eee_jsm);
                           
                            R_dist.clear();
                          
                        }
                    }
                    else
                    {
                        R_dist.push_back(deltaR_js);
                        R_dist.push_back(deltaR_jm);
                        R_dist.push_back(deltaR_sm);
                        int max_R = distance(R_dist.begin(), max_element(R_dist.begin(), R_dist.end()));//pick the longest side to compute the correlators with
                       
                        E3C_hist->Fill(R_dist[max_R],eee_jsm);
                        E3C_pt_hist->Fill(R_dist[max_R],jet_pt,eee_jsm);
                        
                        N3_det_pt_hist_3d->Fill(R_dist[max_R], jet_pt, jet_pt_tru);
                        E3C_det_pt_hist_3d->Fill(R_dist[max_R], jet_pt, jet_pt_tru, eee_jsm);
                    
                        R_dist.clear();
                    }
                }//close m loop
            }//close s loop for the 3 point correlator
            //For loop for EEC
            for(int s=0; s<j ; s++)
            {
                
                double delta_R_js = fConstituents[j].delta_R(fConstituents[s]);
                double ee_js = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));
                
                //Pair cut
                if(fpaircut == 1)
                {
                    double j_eta_2pt_det = fConstituents[j].eta();
                    double s_eta_2pt_det = fConstituents[s].eta();
                    double del_js_eta_2pt_det = abs(j_eta_2pt_det-s_eta_2pt_det);
                    if (del_js_eta_2pt_det < 0.008) continue;
                    else
                    {
  
                        EEC_hist->Fill(delta_R_js,ee_js);
                        EEC_pt_hist->Fill(delta_R_js,jet_pt,ee_js);
                        
                        N2_det_pt_hist_3d->Fill(delta_R_js, jet_pt, jet_pt_tru);
                        EEC_det_pt_hist_3d->Fill(delta_R_js, jet_pt, jet_pt_tru, ee_js);
                    }
                }
                else
                {
                    
                    EEC_hist->Fill(delta_R_js,ee_js);
                    EEC_pt_hist->Fill(delta_R_js,jet_pt,ee_js);
                    
                    N2_det_pt_hist_3d->Fill(delta_R_js, jet_pt, jet_pt_tru);
                    EEC_det_pt_hist_3d->Fill(delta_R_js, jet_pt, jet_pt_tru, ee_js);
                }
                
            }//close s loop for the 2 point correlator
        } //close j loop
    }//close else loop
    
    
    //Looping over truth level jet
    //For jets with 2 constituents
    if(int(fConstituents_tru.size()) == 2)
    {
        for(int j=0; j<int(fConstituents_tru.size()); j++)  //looping over constituents of the fConstituents_tru object
        {
            //For 3 point correlator
            for(int s=0; s<int(fConstituents_tru.size()); s++)
            {
                
                if(s==j) continue; //if s=j this would be 0
                
                double eee_jss_2_tru =((3*fConstituents_tru[j].pt()*fConstituents_tru[s].pt()*fConstituents_tru[s].pt())/(pow(jet_pt_tru,3)));
                double deltaR_jss_2_tru = fConstituents_tru[j].delta_R(fConstituents_tru[s]);
                
                //Pair cut
                if(fpaircut == 1)
                {
                    double j_eta_3pt_2_tru = fConstituents[j].eta();
                    double s_eta_3pt_2_tru = fConstituents[s].eta();
                    double del_js_eta_3pt_2_tru = abs(j_eta_3pt_2_tru-s_eta_3pt_2_tru);
                    if (del_js_eta_3pt_2_tru < 0.008) continue;
                    else
                    {
                       
                        E3C_tru_match_pt_tru->Fill(deltaR_jss_2_tru, jet_pt_tru, eee_jss_2_tru);
                     
                        N3_tru_pt_hist_3d->Fill(deltaR_jss_2_tru, jet_pt_tru, jet_pt);
                        E3C_tru_pt_hist_3d->Fill(deltaR_jss_2_tru, jet_pt_tru, jet_pt, eee_jss_2_tru);
                    }
                }
                else
                {
                    E3C_tru_match_pt_tru->Fill(deltaR_jss_2_tru, jet_pt_tru, eee_jss_2_tru);
                  
                    N3_tru_pt_hist_3d->Fill(deltaR_jss_2_tru, jet_pt_tru, jet_pt);
                    E3C_tru_pt_hist_3d->Fill(deltaR_jss_2_tru, jet_pt_tru, jet_pt, eee_jss_2_tru);
                }
            }//close s loop for the 3 point correlator
            //For 2 point correlator
            for(int s=0; s<j ; s++)
            {
                
                double delta_R_js_2_tru = fConstituents_tru[j].delta_R(fConstituents_tru[s]);
                double ee_js_2_tru = (2*fConstituents_tru[j].pt()*fConstituents_tru[s].pt())/(pow((jet_pt_tru),2));
                
                //Pair cut
                if(fpaircut == 1)
                {
                    double j_eta_2pt_2_tru = fConstituents[j].eta();
                    double s_eta_2pt_2_tru = fConstituents[s].eta();
                    double del_js_eta_2pt_2_tru = abs(j_eta_2pt_2_tru-s_eta_2pt_2_tru);
                    if (del_js_eta_2pt_2_tru < 0.008) continue;
                    else
                    {
                        EEC_tru_match_pt_tru->Fill(delta_R_js_2_tru, jet_pt_tru, ee_js_2_tru);
                       
                        N2_tru_pt_hist_3d->Fill(delta_R_js_2_tru, jet_pt_tru, jet_pt);
                        EEC_tru_pt_hist_3d->Fill(delta_R_js_2_tru, jet_pt_tru, jet_pt, ee_js_2_tru);
                    }
                }
                else
                {
                    EEC_tru_match_pt_tru->Fill(delta_R_js_2_tru, jet_pt_tru, ee_js_2_tru);
                    
                    N2_tru_pt_hist_3d->Fill(delta_R_js_2_tru, jet_pt_tru, jet_pt);
                    EEC_tru_pt_hist_3d->Fill(delta_R_js_2_tru, jet_pt_tru, jet_pt, ee_js_2_tru);
                }
            }//close s loop for eec
        }//close j loop
    }//close if loop
    
    //For jets with more than 2 constituents
    else
    {
        for(int j=0; j<int(fConstituents_tru.size()); j++)  //looping over constituents of the fConstituents_tru object
        {
            for(int s=0; s<int(fConstituents_tru.size()) ; s++)
            {
                if(s==j) continue; //This ensures I don't get stuff like (000) for (jss)
                
                double eee_jss_tru =((3*fConstituents_tru[j].pt()*fConstituents_tru[s].pt()*fConstituents_tru[s].pt())/(pow(jet_pt_tru,3)));
                double deltaR_jss_tru = fConstituents_tru[j].delta_R(fConstituents_tru[s]);
             
                //Pair cut
                if(fpaircut == 1)
                {
                    double j_eta_3pt_tru = fConstituents[j].eta();
                    double s_eta_3pt_tru = fConstituents[s].eta();
                    double del_js_eta_3pt_tru = abs(j_eta_3pt_tru-s_eta_3pt_tru);
                    if (del_js_eta_3pt_tru < 0.008) continue;
                    else
                    {
                        E3C_tru_match_pt_tru->Fill(deltaR_jss_tru,jet_pt_tru, eee_jss_tru);

                        N3_tru_pt_hist_3d->Fill(deltaR_jss_tru, jet_pt_tru, jet_pt);
                        E3C_tru_pt_hist_3d->Fill(deltaR_jss_tru, jet_pt_tru, jet_pt, eee_jss_tru);
                    }
                }
                else
                {
                    E3C_tru_match_pt_tru->Fill(deltaR_jss_tru,jet_pt_tru, eee_jss_tru);
                    
                    N3_tru_pt_hist_3d->Fill(deltaR_jss_tru, jet_pt_tru, jet_pt);
                    E3C_tru_pt_hist_3d->Fill(deltaR_jss_tru, jet_pt_tru, jet_pt, eee_jss_tru);
                }
                //For 3 point correlator
                for( int m=0; m!=j && m!=s; m++)
                {
                    if(s>j) continue;
                    
                    double eee_jsm_tru = ((6*fConstituents_tru[j].pt()*fConstituents_tru[s].pt()*fConstituents_tru[m].pt())/(pow(jet_pt_tru,3)));
                    double deltaR_js_tru = fConstituents_tru[j].delta_R(fConstituents_tru[s]);
                   
                    double deltaR_jm_tru = fConstituents_tru[j].delta_R(fConstituents_tru[m]);
                   
                    double deltaR_sm_tru = fConstituents_tru[s].delta_R(fConstituents_tru[m]);
                    
                    R_dist_part.push_back(deltaR_js_tru);
                    R_dist_part.push_back(deltaR_jm_tru);
                    R_dist_part.push_back(deltaR_sm_tru);
                    
                    int max_R_tru = distance(R_dist_part.begin(), max_element(R_dist_part.begin(), R_dist_part.end()));//pick the longest side to compute the correlators with
                    
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_3pt_tru = fConstituents[j].eta();
                        double m_eta_3pt_tru = fConstituents[m].eta();
                        double s_eta_3pt_tru = fConstituents[s].eta();
                        double del_jm_eta_3pt_tru = abs(j_eta_3pt_tru-m_eta_3pt_tru);
                        double del_sm_eta_3pt_tru = abs(s_eta_3pt_tru-m_eta_3pt_tru);
                        double del_js_eta_3pt_tru = abs(j_eta_3pt_tru-s_eta_3pt_tru);
                        if (del_jm_eta_3pt_tru < 0.008 || del_sm_eta_3pt_tru < 0.008 || del_js_eta_3pt_tru < 0.008) continue;
                        else
                        {
                            E3C_tru_match_pt_tru->Fill(R_dist_part[max_R_tru], jet_pt_tru, eee_jsm_tru);
                           
                            N3_tru_pt_hist_3d->Fill(R_dist_part[max_R_tru], jet_pt_tru, jet_pt);
                            E3C_tru_pt_hist_3d->Fill(R_dist_part[max_R_tru], jet_pt_tru, jet_pt, eee_jsm_tru);
                          
                            R_dist_part.clear();
                        }
                    }
                    else
                    {
                        E3C_tru_match_pt_tru->Fill(R_dist_part[max_R_tru], jet_pt_tru, eee_jsm_tru);
                       
                        N3_tru_pt_hist_3d->Fill(R_dist_part[max_R_tru], jet_pt_tru, jet_pt);
                        E3C_tru_pt_hist_3d->Fill(R_dist_part[max_R_tru], jet_pt_tru, jet_pt, eee_jsm_tru);
                       
                        R_dist_part.clear();
                    }
                }//close m loop
            }//close s loop for the 3 point correlator
            //For loop for EEC
            for(int s=0; s<j ; s++)
            {
                double delta_R_js_tru = fConstituents_tru[j].delta_R(fConstituents_tru[s]);
                double ee_js_tru = (2*fConstituents_tru[j].pt()*fConstituents_tru[s].pt())/(pow((jet_pt_tru),2));
                
                //Pair cut
                if(fpaircut == 1)
                {
                    double j_eta_2pt_tru = fConstituents[j].eta();
                    double s_eta_2pt_tru = fConstituents[s].eta();
                    double del_js_eta_2pt_tru = abs(j_eta_2pt_tru-s_eta_2pt_tru);
                    if (del_js_eta_2pt_tru < 0.008) continue;
                    else{
                        
                        EEC_tru_match_pt_tru->Fill(delta_R_js_tru,jet_pt_tru,ee_js_tru);
                     
                        N2_tru_pt_hist_3d->Fill(delta_R_js_tru, jet_pt_tru, jet_pt);
                        EEC_tru_pt_hist_3d->Fill(delta_R_js_tru, jet_pt_tru, jet_pt, ee_js_tru);
                    }
                }
                else
                {
                    EEC_tru_match_pt_tru->Fill(delta_R_js_tru,jet_pt_tru,ee_js_tru);
                   
                    N2_tru_pt_hist_3d->Fill(delta_R_js_tru, jet_pt_tru, jet_pt);
                    EEC_tru_pt_hist_3d->Fill(delta_R_js_tru, jet_pt_tru, jet_pt, ee_js_tru);
                }
            }//close s loop for the 2 point correlator
        } //close j loop
    }//close else loop
    
    
    //catch (fastjet::Error)
    // {
    //    AliError(" [w] FJ Exception caught.");
    //    // return -1;
    // } //end error message
    return;
}

    
//EEC computation-------------------------------------------------------
void AliAnalysisTaskJetsEEC::ComputeEEC(AliEmcalJet *fJet, AliJetContainer *fJetCont)
    //(jet, jet container, vector of jets, vector of constituents within each jet?)
{
    //General EEC computation: Need to loop over jets, then loop within a single jet to compute EECs.
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
        double jet_pt = fJet->Pt();
        
        if(fpTcorr == 1)
         {   jet_pt = jet_pt/(0.85); //applying JES correction to jet pT spectra to study spectra shape dependence of ENC
             jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
             
         }
         
         else
         {
             jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
         
         }
        
        //For jets with 2 constituents
        if(int(fConstituents.size()) == 2)
        {
            for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
            {
                //For 3 point correlator
                for(int s=0; s<j ; s++)
                {
                    if(s==j) continue; //this ensure i dont get 000
                    
                    double eee_jss_2 =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                    double deltaR_jss_2 = fConstituents[j].delta_R(fConstituents[s]);
                   
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_3pt_2 = fConstituents[j].eta();
                        double s_eta_3pt_2 = fConstituents[s].eta();
                        double del_js_eta_3pt_2 = abs(j_eta_3pt_2-s_eta_3pt_2);
                        if (del_js_eta_3pt_2 < 0.008) continue;
                        else
                        {
                            E3C_hist->Fill(deltaR_jss_2,eee_jss_2);
                            E3C_pt_hist->Fill(deltaR_jss_2,jet_pt,eee_jss_2);
                        }
                    }
                    else
                    {
                        E3C_hist->Fill(deltaR_jss_2,eee_jss_2);
                        E3C_pt_hist->Fill(deltaR_jss_2,jet_pt,eee_jss_2);
                    }
                }//close s loop for the 3 point correlator
                //For 2 point correlator
                for(int s=0; s<j ; s++)
                {
                    
                    double delta_R_js_2 = fConstituents[j].delta_R(fConstituents[s]);
                    double ee_js_2 = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_2pt_2 = fConstituents[j].eta();
                        double s_eta_2pt_2 = fConstituents[s].eta();
                        double del_js_eta_2pt_2 = abs(j_eta_2pt_2-s_eta_2pt_2);
                        if (del_js_eta_2pt_2 < 0.008) continue;
                        else
                        {
                            EEC_hist->Fill(delta_R_js_2,ee_js_2);
                            EEC_pt_hist->Fill(delta_R_js_2,jet_pt,ee_js_2);
                        }
                    }
                    else
                    {
                        EEC_hist->Fill(delta_R_js_2,ee_js_2);
                        EEC_pt_hist->Fill(delta_R_js_2,jet_pt,ee_js_2);
                    }
                    
                }//close s loop for eec
            }//close j loop
        }//close if loop
        
        //For jets with more than 2 constituents
        else
        {
            for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
            {
                
                for(int s=0; s<int(fConstituents.size()) ; s++)
                {
                    if(s==j) continue; //This ensures I don't get stuff like (000) for (jss)
                    
                    double eee_jss =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                    double deltaR_jss = fConstituents[j].delta_R(fConstituents[s]);
                  
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_3pt = fConstituents[j].eta();
                        double s_eta_3pt = fConstituents[s].eta();
                        double del_js_eta_3pt = abs(j_eta_3pt-s_eta_3pt);
                        if (del_js_eta_3pt < 0.008) continue;
                        else
                        {
                            E3C_hist->Fill(deltaR_jss,eee_jss);
                            E3C_pt_hist->Fill(deltaR_jss,jet_pt,eee_jss);
                        }
                    }
                    else
                    {

                        E3C_hist->Fill(deltaR_jss,eee_jss);
                        E3C_pt_hist->Fill(deltaR_jss,jet_pt,eee_jss);
                    }
                    //For 3 point correlator
                    for( int m=0; m!=j && m!=s; m++)
                    {
                        if(s>j) continue;
                        
                            double eee_jsm = ((6*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[m].pt())/(pow(jet_pt,3)));
                            double deltaR_js = fConstituents[j].delta_R(fConstituents[s]);
                           
                            double deltaR_jm = fConstituents[j].delta_R(fConstituents[m]);
                         
                            double deltaR_sm = fConstituents[s].delta_R(fConstituents[m]);
                         
                            R_dist.push_back(deltaR_js);
                            R_dist.push_back(deltaR_jm);
                            R_dist.push_back(deltaR_sm);
                            
                            int max_R = distance(R_dist.begin(), max_element(R_dist.begin(), R_dist.end()));//pick the longest side to compute the correlators with
                            
                            if(fpaircut==1)
                            {
                                //Pair cut
                                double j_eta_3pt = fConstituents[j].eta();
                                double m_eta_3pt = fConstituents[m].eta();
                                double s_eta_3pt = fConstituents[s].eta();
                                double del_jm_eta_3pt = abs(j_eta_3pt-m_eta_3pt);
                                double del_sm_eta_3pt = abs(s_eta_3pt-m_eta_3pt);
                                double del_js_eta_3pt = abs(j_eta_3pt-s_eta_3pt);
                                if (del_jm_eta_3pt < 0.008 || del_sm_eta_3pt < 0.008 || del_js_eta_3pt < 0.008) continue;
                                else
                                {
                                
                                    E3C_hist->Fill(R_dist[max_R],eee_jsm);
                                    E3C_pt_hist->Fill(R_dist[max_R],jet_pt,eee_jsm);
//
                                    R_dist.clear();
                                }
                            }
                            else
                            {
                                E3C_hist->Fill(R_dist[max_R],eee_jsm);
                                E3C_pt_hist->Fill(R_dist[max_R],jet_pt,eee_jsm);
//
                                R_dist.clear();
                            }

                    }//close m loop
                }//close s loop for the 3 point correlator
                //For loop for EEC
                for(int s=0; s<j ; s++)
                {
                    
                    double delta_R_js = fConstituents[j].delta_R(fConstituents[s]);
                    double ee_js = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));
                    
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_2pt = fConstituents[j].eta();
                        double s_eta_2pt = fConstituents[s].eta();
                        double del_js_eta_2pt = abs(j_eta_2pt-s_eta_2pt);
                        if (del_js_eta_2pt < 0.008) continue;
                        else
                        {
                            EEC_hist->Fill(delta_R_js,ee_js);
                            EEC_pt_hist->Fill(delta_R_js,jet_pt,ee_js);
                        }
                    }
                    else
                    {
                        EEC_hist->Fill(delta_R_js,ee_js);
                        EEC_pt_hist->Fill(delta_R_js,jet_pt,ee_js);
                    }
                }//close s loop for the 2 point correlator
            } //close j loop
        }//close else loop

    
    //catch (fastjet::Error)
    // {
    //    AliError(" [w] FJ Exception caught.");
    //    // return -1;
    // } //end error message
    return;
}


//_________________________________________________________________
void AliAnalysisTaskJetsEEC::RunChanged(Int_t newrun)
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
Bool_t AliAnalysisTaskJetsEEC::RetrieveEventObjects()
{
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;
  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskJetsEEC::Terminate(Option_t *)
{
  // Called once at the end of the analysis.
  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available");
  //   return;
  // }
}

