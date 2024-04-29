#ifndef ALIANALYSISTASKJetsEEC_H
#define ALIANALYSISTASKJetsEEC_H

class TH1;
class TH2;
class TH3;
class TH3F;
class TTree;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisManager;
class AliJetContainer;
class AliEmcalJetFinder;
class AliFJWrapper;
#include "AliAnalysisTaskEmcalJet.h"
#include "AliFJWrapper.h"
#include <vector>
#include <TObjArray.h>
#include "AliReducedMatchedTrack.h"

class AliAnalysisTaskJetsEEC : public AliAnalysisTaskEmcalJet {
public:
    enum JetShapeType {
        kMCTrue = 0,  // generated jets only
        kTrueDet = 1, // detector and generated jets
        kData = 2,    // raw data
        kDetEmb = 3,  // detector embedded jets
        kDetEmbPart = 4,
        kPythiaDef = 5,
        kDetEmbPartPythia = 6,
        kGenOnTheFly = 7
    };
    enum JetShapeSub { kNoSub = 0, kConstSub = 1, kDerivSub = 2, kEventSub=3 };
    enum JetSelectionType { kInclusive = 0, kRecoil = 1 };
    
    enum DerivSubtrOrder { kSecondOrder = 0, kFirstOrder = 1 };
    
    AliAnalysisTaskJetsEEC();
    AliAnalysisTaskJetsEEC(const char *name);
    virtual ~AliAnalysisTaskJetsEEC();
    
    void UserCreateOutputObjects();
    void Terminate(Option_t *option);
    
    // Setters
    void SetJetContainer(Int_t c) { fContainer = c; }
    void SetJetShapeType(JetShapeType t) { fJetShapeType = t; }
    void SetJetShapeSub(JetShapeSub t) { fJetShapeSub = t; }
    void SetJetSelection(JetSelectionType t) { fJetSelection = t; }
    void SetJetPtThreshold(Float_t f) { fPtThreshold = f; }
    //  void SetRMatching(Float_t f) { fRMatching = f; } //MIGHT NEED
    void SetCentralitySelectionOn(Bool_t t) { fCentSelectOn = t; }
    void SetOneConstSelectionOn(Bool_t t) { fOneConstSelectOn = t; }
    void SetCheckTracksOn(Bool_t t) { fTrackCheckPlots = t; }
    void SetCheckResolution(Bool_t t) { fCheckResolution = t; }
    void SetMinPtConst(Float_t t) { fMinPtConst = t;}
    void SetMinTrackPtEncs(Float_t t) { fMinENCtrackPt = t;}
    void SetHardCutoff(Float_t t) { fHardCutoff = t; }
    void SetDoTwoTrack(Bool_t t) { fDoTwoTrack = t; }
    void SetCutDoubleCounts(Bool_t t) {fCutDoubleCounts = t;}
    //  void SetDoAreaIterative(Bool_t t) { fDoAreaIterative = t; }
    void SetPowerAlgorithm(Float_t t) { fPowerAlgo = t; }
    void SetMinCentrality(Float_t t) { fCentMin = t; }
    void SetMaxCentrality(Float_t t) { fCentMax = t; }
    void SetDerivativeSubtractionOrder(Int_t c) { fDerivSubtrOrder = c; }
    void SetDetLevelJetsOn(Bool_t t) { fStoreDetLevelJets = t; }
    //  void SetMatchRadius(Float_t t) { fMatchR = t; } //MIGHT NEED
    void SetStoreTrig(Bool_t t) {fStoreTrig = t;}
    
    void SetFillEncMC(Bool_t t) { fDoFillEncMC = t; } //For MC
    void SetpTcorrect(Int_t c) {fpTcorr = c;}
    void SetPairCut(Int_t c) {fpaircut = c;}
    
    //Adding for tracks at the particle level (not jet level)
    void SetFastSimPairEff(Int_t c) {fpairfastsim = c;}
    void SetUnfolding(Int_t c){fUnfolding = c;} //For unfolding 
    void SetMaxTrackPtCut(Float_t c) {fMaxPtTrack = c;} //set max value of track pT
    void SetGeneratorLevelName(const char* name) { fGeneratorLevelName = name; }
    void SetDetectorLevelName(const char* name) { fDetectorLevelName  = name; }
    
protected:
    Bool_t RetrieveEventObjects();
    Bool_t Run();
    Bool_t FillHistograms();
    
    Float_t GetJetMass(AliEmcalJet *jet, Int_t jetContNb);
    Float_t Angularity(AliEmcalJet *jet, Int_t jetContNb); //MIGHT NOT NEED
    Float_t GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb); //MIGHT NOT NEED
    
    int GetConstituentID(int constituentIndex, const AliVParticle* part, AliEmcalJet * jet);
    double delR(const fastjet::PseudoJet& ps1, const fastjet::PseudoJet& ps2);
    int GetCharge(int detCharge, const AliVParticle* part);
    int GetPartCharge(int partCharge, const AliVParticle* part_gen);
    
    Double_t GetDownscaleWeight(string tstring);
    void ComputeEEC(AliEmcalJet *fJet, AliJetContainer *fJetCont);
    
    void ComputeEncMC(AliEmcalJet *fJet, AliJetContainer *fJetCont, AliEmcalJet *fJet_tru, Int_t km); //MC and det correlations
    
    void ComputeJetTrackMatch(AliEmcalJet *fJet, AliJetContainer *fJetCont, AliEmcalJet *fJet_tru, AliJetContainer *fJet_truCont, Int_t km);
    
//    void FillMatchedParticlesTHnSparse(Double_t partEta, Double_t partPhi, Double_t partPt, Double_t trackEta, Double_t trackPhi, Double_t trackPt, Byte_t trackType); //for tracks matched between detector and truth level
    
    void MatchTracks();//track matching
    void ComputeDelqpt(); //pair eff as a function of q/pt
    
    void RunChanged(Int_t nr);
    Int_t fContainer; ///< jets to be analyzed 0 for Base, 1 for subtracted.
    Float_t fMinFractionShared; ///< only fill histos for jets if shared fraction larger than x
    JetShapeType fJetShapeType; ///< jet type to be used
    JetShapeSub fJetShapeSub;   ///< jet subtraction to be used //MIGHT NOT NEED THIS
    JetSelectionType fJetSelection; ///< Jet selection: inclusive/recoil jet
    Float_t fPtThreshold; ///<
    Float_t fMinENCtrackPt; ///< min track pT for the EECs
    
    
    Bool_t fCentSelectOn;      ///< switch on/off centrality selection
    Float_t fCentMin;          ///< min centrality value
    Float_t fCentMax;          ///< max centrality value
    Bool_t fOneConstSelectOn;  ///< switch on/off one constituent selection
    Bool_t fTrackCheckPlots;   ///< switch on qa plots
    
    Bool_t fCheckResolution;   ///< check subjet energy resolution
    
    Float_t fMinPtConst;       ///< constituent pt cutoff
    Float_t fHardCutoff;       ///< hard cutoff in the iterative declustering
    Bool_t fDoTwoTrack;        ///< switch to consider 2 track effects
    Bool_t fCutDoubleCounts;   ///< turn off to avoid true-hybrid cuts to suppress double counting
    //Bool_t fDoAreaIterative;   ///<  subtract the area in the declustering
    Float_t fPowerAlgo;        ///< power of the generickt algorithm
    Float_t fPhiCutValue;      ///< cuts from HBT
    Float_t fEtaCutValue;      ///< cuts from HBT
    Int_t fDerivSubtrOrder; //What does this do?????????
    Bool_t fStoreDetLevelJets; ///< store the detector level jet quantities
    
    
    ///MIGHT NEED MIGHT NEED MIGHT NEED
    Bool_t fDoFillEncMC;      ///< to fill the matched mc plane
    //  Float_t fRMatching; ///<
    
    
    //MIGHT NEED MIGHT NEED
    //  Float_t fMatchR; ///<the matching radius
    
    Bool_t fStoreTrig; ///<storing the trigger class
    //  Bool_t fMatch; ///< do the matching in the task
    
    Int_t fpTcorr; ///<flag for pT migration checks
    Int_t fpaircut; ///<flag for pT pair cut
    Int_t fpairfastsim; ///<flag for pair eff for fast sim
    Int_t fUnfolding; ///<flag for unfolding 
    Int_t fMatchJetTrack; ///<flag for pair eff for fast sim
    Int_t fMissJetTrack; ///<flag for looking at missed tracks in matched jets
    Int_t fFakeJetTrack; ///<flag for looking at missed tracks in matched jets
    Float_t fMaxPtTrack; ///< max track pt cutoff

    Double_t fJet_pt_det; ///< jet pt det 
    Double_t fJet_pt_tru; ///< jet pt tru
    Double_t fTrack_pt_det; ///< track pt det 
    Double_t fTrack_pt_tru; ///< track pt tru
    Double_t fTrack_eta_det; ///< track eta det 
    Double_t fTrack_eta_tru; ///< track eta tru 
    Double_t fTrack_phi_det; ///< track phi det 
    Double_t fTrack_phi_tru; ///< track phi tru

    //Adding for tracks at the particle level (not jet level)
    TString fGeneratorLevelName;
    TString fDetectorLevelName;
    
    AliMCParticleContainer* fGeneratorLevel;
    AliTrackContainer* fDetectorLevel;
    AliJetContainer* fJet_truCont;
    
    //Histograms
    TH1D *jet_pt_hist; //!<! initializing histogram with jet pt
    
    TH1D *EEC_hist; //!<! initializing histogram with correlators
    TH2D *EEC_pt_hist; //!<! initializing 2D histogram with correlators and jet pt
 
    TH1D *E3C_hist; //!<! initializing histogram with correlators
    TH2D *E3C_pt_hist; //!<! initializing 2D histogram with correlators and jet pt
    
    //For corrections
    TH3D *EEC_det_pt_hist_3d; //!<! initializing 3d histogram with det level EEC and det and true jet pt
    TH3D *EEC_tru_pt_hist_3d; //!<! initializing 3d histogram with tru level EEC and det and true jet pt
    
    TH3D *E3C_det_pt_hist_3d; //!<! initializing 3d histogram with det level E3C and det and true jet pt
    TH3D *E3C_tru_pt_hist_3d; //!<! initializing 3d histogram with tru level E3C and det and true jet pt
    
    TH3D *N2_det_pt_hist_3d; //!<! initializing 3d histogram with det level EEC num pairs and det and true jet pt
    TH3D *N2_tru_pt_hist_3d; //!<! initializing 3d histogram with tru level EEC num pairs and det and true jet pt
    
    TH3D *N3_det_pt_hist_3d;  //!<! initializing 3d histogram with det level E3C num pairs and det and true jet pt
    TH3D *N3_tru_pt_hist_3d; //!<! initializing 3d histogram with tru level E3C num pairs and det and true jet pt
    
    TH2D *EEC_det_match_pt_det; //!<! initializing 2D histogram of det level EEC and det jet pt
    TH2D *EEC_tru_match_pt_tru; //!<! initializing 2D histogram of tru level EEC and tru jet pt
    
    TH2D *E3C_det_match_pt_det; //!<! initializing 2D histogram of det level E3C and det jet pt
    TH2D *E3C_tru_match_pt_tru; //!<! initializing 2D histogram of tru level E3C and tru jet pt
    
    TH1D *pt_tru; //!<! all tru level
    TH1D *pt_tru_match; //!<! tru level matched
    TH1D *pt_det; //!<! all det level
    TH1D *pt_det_match; //!<! det level matched
    
    TH1D *test_hist; //!<! test histogram
    
    TH2D *R_matrix; //!<! response matrix for jet pt
    
    TH2D *JES; //!<! jet energy scale
    TH2D *JES_scaled; //!<! jet energy scale for scaled pt_det
    TH2D *JER; //!<! jet energy resolution
    
    TH3D *pair_det_EEC;//!<! histogram for computing pair efficiency at det level
    TH3D *pair_tru_EEC;//!<! histogram for computing pair efficiency at truth level
    
    TH3D *pair_det_E3C;//!<! histogram for computing pair efficiency at det level
    TH3D *pair_tru_E3C;//!<! histogram for computing pair efficiency at truth level
    
    TH2D *qpt_tru;//!<! q/pt vs R histogram for input to fast sim
    TH2D *qpt_det;//!<! q/pt vs R histogram for input to fast sim
    TH1D *track_pt_tru; //!<! all tru level tracks
    TH1D *track_pt_det; //!<! all det level tracks
    TH2D *track_pt_matched; //!<! all det level tracks matched with truth level
    
    TH3D *R_match_eec; //!<! pair distance of all det level tracks in det jet matched with truth level tracks in truth jet eec as a function of jet pT
    TH3D *wt_match_eec; //!<! pair wts of all det level tracks in det jet matched with truth level tracks in truth jet eec as a function of jet pT
    TH3D *R_match_e3c; //!<! pair distance of all det level tracks in det jet matched with truth level tracks in truth jet e3c as a function of jet pT
    TH3D *wt_match_e3c; //!<! pair wts of all det level tracks in det jet matched with truth level tracks in truth jet e3c as a function of jet pT
    
    TH1D *qpt_tru1;//!<! q/pt vs R histogram for input to fast sim
    TH1D *qpt_tru2;//!<! q/pt vs R histogram for input to fast sim
    //!
    TH1D *pt_tru1;//!<! q/pt vs R histogram for input to fast sim
    TH1D *pt_tru2;//!<! q/pt vs R histogram for input to fast sim
    
    TH2D *eec_Mm; //!<! initializing 2D histogram of EEC matched and missed
    TH2D *eec_mm; //!<! initializing 2D histogram of tru level EEC missed and missed
    TH2D *e3c_MMm; //!<! initializing 2D histogram of det level E3C matched matched missed
    TH2D *e3c_Mmm; //!<! initializing 2D histogram of tru level E3C  matched missed missed
    TH2D *e3c_mmm; //!<! initializing 2D histogram of det level E3C missed missed missed
    
    TH2D *eec_Mf; //!<! initializing 2D histogram of det level EEC matched and fake
    TH2D *eec_ff; //!<! initializing 2D histogram of tru level EEC fake and fake
    TH2D *e3c_MMf; //!<! initializing 2D histogram of det level E3C matched matched fake
    TH2D *e3c_Mff; //!<! initializing 2D histogram of tru level E3C mateched fake fake
    TH2D *e3c_fff; //!<! initializing 2D histogram of det level E3C fake fake fake

    TH2D *eec_matched_det; //!<! initializing 2D histogram of det level EEC and det jet pt
    TH2D *eec_matched_tru; //!<! initializing 2D histogram of tru level EEC and tru jet pt
    TH2D *e3c_matched_det; //!<! initializing 2D histogram of det level E3C and det jet pt
    TH2D *e3c_matched_tru; //!<! initializing 2D histogram of tru level E3C and tru jet pt
    
    TH3D *wt_res_eec; //!<! initializing 2D histogram for wt resolution EEC
    TH3D *wt_res_e3c; //!<! initializing 2D histogram for wt resolution E3C
    TH3D *R_res_eec; //!<! initializing 2D histogram for R resolution EEC
    TH3D *R_res_e3c; //!<! initializing 2D histogram for R resolution E3C
    TH3D *wtnojet_match_eec;//!<! initializing 2D histogram for EEC wt without jet pt
    TH3D *wtnojet_match_e3c;//!<! initializing 2D histogram for E3C wt without jet pt
    TH3D *wtnojet_res_eec; //!<! initializing 2D histogram for wt resolution with jet pt excluded EEC
    TH3D *wtnojet_res_e3c; //!<! initializing 2D histogram for wt resolution with jet pt excluded E3C
    
    TH3D *R_match_eec_tru; //!<! pair distance of all det level tracks in det jet matched with truth level tracks in truth jet eec as a function of jet pT tru
    TH3D *wt_match_eec_tru; //!<! pair wts of all det level tracks in det jet matched with truth level tracks in truth jet eec as a function of jet pT tru
    TH3D *R_match_e3c_tru; //!<! pair distance of all det level tracks in det jet matched with truth level tracks in truth jet e3c as a function of jet pT tru
    TH3D *wt_match_e3c_tru; //!<! pair wts of all det level tracks in det jet matched with truth level tracks in truth jet e3c as a function of jet pT tru

    TH3D *wt_res_eec_tru; //!<! initializing 3D histogram for wt resolution EEC as a function of jet pT tru
    TH3D *wt_res_e3c_tru; //!<! initializing 3D histogram for wt resolution E3C as a function of jet pT tru
    TH3D *R_res_eec_tru; //!<! initializing 3D histogram for R resolution EEC as a function of jet pT tru
    TH3D *R_res_e3c_tru; //!<! initializing 3D histogram for R resolution E3C as a function of jet pT tru
    TH3D *wtnojet_match_eec_tru;//!<! initializing 3D histogram for EEC wt without jet pt as a function of jet pT tru
    TH3D *wtnojet_match_e3c_tru;//!<! initializing 3D histogram for E3C wt without jet pt as a function of jet pT tru
    TH3D *wtnojet_res_eec_tru; //!<! initializing 3D histogram for wt resolution with jet pt excluded EEC as a function of jet pT tru
    TH3D *wtnojet_res_e3c_tru; //!<! initializing 3D histogram for wt resolution with jet pt excluded E3C as a function of jet pT tru

    TH1D *constituentId; //!<! initializing 1D histogram for debugging
    TH3D *R_res_eec_tru_debug; //!<! initializing 3D histogram for R resolution EEC as a function of jet pT tru debug
    TH2D *track_pt_res_debug; //!<! initializing 3D histogram for pt resolution debug
    TH3D *track_eta_debug; //!<! initializing 3D histogram for eta match debug
    TH3D *track_rap_debug; //!<! initializing 3D histogram for rap match debug
    TH2D *track_phi_debug; //!<! initializing 3D histogram phi match debug
    TH2D *track_R_debug; //!<! initializing 3D histogram for R match debug
    TH2D *track_R_debug_rap; //!<! initializing 3D histogram for R-rap debug
    
    TH3D *track_eta_res_debug; //!<! initializing 3D histogram for eta resolution of single matched track as a function of jet pT tru debug
    TH3D *track_rap_res_debug; //!<! initializing 3D histogram for rap resolution of single matched track as a function of jet pT tru debug
    TH3D *track_phi_res_debug; //!<! initializing 3D histogram for phi resolution of single matched track as a function of jet pT tru debug
    TH3D *track_R_res_debug; //!<! initializing 3D histogram for R resolution of single matched track as a function of jet pT tru debug
    TH3D *track_R_rap_res_debug; //!<! initializing 3D histogram for R-rap resolution of single matched track as a function of jet pT tru debug

    TH3D *R_match_eec_tru_rap; //!<! initializing 3D histogram for R resolution of single matched track as a function of jet pT tru debug

    TH2D *track_pt_response_debug; //!<! initializing 2D histogram for debug pt
    TH2D *track_pt_wt_res_debug; //!<! initializing 2D histogram for debug wt
    TH2D *track_pt_wt_response_debug; //!<! initializing 2D histogram for debug wt_response
    
    TH3D* jetpt_res_w_R; //!<! initializing 3D histogram for R resolution of single matched track as a function of jet pT tru debug
    TH3D* jetpt_res_w_wt; //!<! initializing 3D histogram for weight resolution of single matched track as a function of jet pT tru debug
    
    TH3D* OptUn_eec; //!<! initializing 3D histogram for optimizing binning eec
    TH3D* OptUn_e3c; //!<! initializing 3D histogram for optimizing binning e3c

    TTree *fTreeMatchTracks; ///< Tree with matched tracks from MC
private:
    AliAnalysisTaskJetsEEC(
                           const AliAnalysisTaskJetsEEC &); // not implemented
    AliAnalysisTaskJetsEEC &
    operator=(const AliAnalysisTaskJetsEEC &); // not implemented
    
    ClassDef(AliAnalysisTaskJetsEEC, 70) //change this to 71 if you add something new
};
#endif

