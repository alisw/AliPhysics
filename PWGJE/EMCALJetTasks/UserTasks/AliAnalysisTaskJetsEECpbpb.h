#ifndef ALIANALYSISTASKJetsEECPBPB_H
#define ALIANALYSISTASKJetsEECPBPB_H

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
class AliAnalysisTaskRho;
#include "AliAnalysisTaskEmcalJet.h"
#include "AliFJWrapper.h"
#include <vector>
#include <TObjArray.h>
#include "AliReducedMatchedTrack.h"
#include "TClonesArray.h"

class AliAnalysisTaskJetsEECpbpb : public AliAnalysisTaskEmcalJet {
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
    
    AliAnalysisTaskJetsEECpbpb();
    AliAnalysisTaskJetsEECpbpb(const char *name);
    virtual ~AliAnalysisTaskJetsEECpbpb();

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
    void SetMatchRadius(Float_t t) { fMatchR = t; } //MIGHT NEED
    void SetMinJetPtSub(Float_t t) { fjetMinPtSub = t; } //MIGHT NEED
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
    void SetJetPtMin(Float_t t) { fJetPtMin = t; }
    void SetJetMinArea(Float_t t) { fjetMinArea = t; }

    void SetDoEmbedding(Bool_t t) { fDoEmbedding = t; }
    void SetDoPerpCone(Bool_t t) { fDoPerpCone = t; }
    void SetDoRandCone(Bool_t t) {fDoRandCone = t;}
    void SetCoutStatements(Bool_t t) {fCout = t;}

    
    
protected:
    Bool_t RetrieveEventObjects();
    Bool_t Run();
    Bool_t FillHistograms();
    // Bool_t FillEmbJetsEEC();
    
    // Float_t GetJetMass(AliEmcalJet *jet, Int_t jetContNb);
    // Float_t Angularity(AliEmcalJet *jet, Int_t jetContNb); //MIGHT NOT NEED
    // Float_t GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb); //MIGHT NOT NEED
    
    int GetConstituentID(int constituentIndex, const AliVParticle* part, AliEmcalJet * jet);
    double delR(const fastjet::PseudoJet& ps1, const fastjet::PseudoJet& ps2);
    // int GetCharge(int detCharge, const AliVParticle* part);
    // int GetPartCharge(int partCharge, const AliVParticle* part_gen);
    
    Double_t GetDownscaleWeight(string tstring);
    void ComputeENC(AliEmcalJet *fJet, double ptSub, AliJetContainer *fJetCont);
    
    void ComputeEncMC(AliEmcalJet *fJet, AliJetContainer *fJetCont, AliEmcalJet *fJet_tru, Int_t km); //MC and det correlations
    
    // void ComputeJetTrackMatch(AliEmcalJet *fJet, AliJetContainer *fJetCont, AliEmcalJet *fJet_tru, AliJetContainer *fJet_truCont, Int_t km);
    
    void FillEmbJetsEEC(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2, double jetpt, float pt, bool typeSame, std::string type, float MCweight, double n, double bkgIndex1, double bkgIndex2);

    void FillEmbJetsE3C(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2, std::vector<fastjet::PseudoJet> particles3, double jetpt, float pt,string typeSame, std::string type, float MCweight, double n, double bkgIndex1, double bkgIndex2,  double bkgIndex3);

    std::vector<fastjet::PseudoJet> FindThermalCone(AliEmcalJet *fJetEmb, AliJetContainer *fJetContEmb, double ptSub, AliEmcalJet *fJet, AliJetContainer *fJet_detCont, AliEmcalJet *fJet_tru, AliJetContainer *fJet_truCont, std::string axisType);

    double Calculate_pX(double pT, double eta, double phi);
    double Calculate_pY(double pT, double eta, double phi);
    double Calculate_pZ(double pT, double eta, double phi);
    double Calculate_E(double pT, double eta, double phi);
    
    void MatchTracks();//track matching
    // void ComputeDelqpt(); //pair eff as a function of q/pt

    void RunChanged(Int_t nr);
    Int_t fContainer; ///< jets to be analyzed 0 for Base, 1 for subtracted.
    Float_t fMinFractionShared; ///< only fill histos for jets if shared fraction larger than x
    JetShapeType fJetShapeType; ///< jet type to be used
    JetShapeSub fJetShapeSub;   ///< jet subtraction to be used 
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
    Float_t fMatchR; ///<the matching radius
    Float_t fjetMinPtSub; ///<min subtracted jet pT 
    Float_t fjetMinArea; ///<min area of subtracted jet 
    
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
    Float_t fJetPtMin; ///<flag for min jet pT allowed

    Bool_t fDoEmbedding; ///<Do embedding or not
    Bool_t fIsEmbeddedEvent; ///<embedded event or not 
    Bool_t fDoPerpCone; ///<create perp cone
    Bool_t fDoRandCone; ///<create random cone
    Bool_t fCout; ///<print cout statements 

    Double_t fJet_pt_det; ///< jet pt det 
    Double_t fJet_pt_tru; ///< jet pt tru
    Double_t fTrack_pt_det; ///< track pt det 
    Double_t fTrack_pt_tru; ///< track pt tru
    Double_t fTrack_eta_det; ///< track eta det 
    Double_t fTrack_eta_tru; ///< track eta tru 
    Double_t fTrack_phi_det; ///< track phi det 
    Double_t fTrack_phi_tru; ///< track phi tru
    Double_t fTrack_pt_miss; ///< miss track pt 
    Double_t fTrack_pt_fake; ///< fake track pt 
    Double_t fTrack_eta_miss; ///< miss track eta 
    Double_t fTrack_eta_fake; ///< fake track eta  
    Double_t fTrack_phi_miss; ///< miss track phi  
    Double_t fTrack_phi_fake; ///< fake track phi

    Double_t fJet_pt_dat; ///< jet pt dat
    Double_t fTrack_pt_dat; ///< track pt dat
    Double_t fTrack_eta_dat; ///< track eta dat
    Double_t fTrack_phi_dat; ///< track phi dat
    
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
    TH1F* NumJetEvent; //!<! number of jet events

    TH3D* h_MJ;//!<!matchedjet all tracks fake + real
    TH3D* h_MJ0;//!<!matchedjet 0 real tracks
    TH3D* h_MJ1;//!<!matchedjet 1 real track
    TH3D* h_MJ2;//!<!matchedjet 2 real tracks

    TH3D* h_MJ_tru;//!<!matchedjet all tracks fake + real
    TH3D* h_MJ0_tru;//!<!matchedjet 0 real tracks
    TH3D* h_MJ1_tru;//!<!matchedjet 1 real track
    TH3D* h_MJ2_tru;//!<!matchedjet 2 real tracks

    TH3D* h_MB1;//!<!bkg-bkg correlations
    TH3D* h_JMB;//!<!jet-bkg correlations
    TH3D* h_SMB;//!<!signal-bkg correlations
    TH3D* h_BMB;//!<!jetbkg-bkg correlations
    TH3D* h_MB1MB2;//!<!bkg1-bkg2 correlations

    TH3D* h_MB1_tru;//!<!bkg-bkg correlations true
    TH3D* h_JMB_tru;//!<!jet-bkg correlations true
    TH3D* h_SMB_tru;//!<!signal-bkg correlations true
    TH3D* h_BMB_tru;//!<!jetbkg-bkg correlations true
    TH3D* h_MB1MB2_tru;//!<!bkg1-bkg2 correlations true

    TH3D* h_MJ_e3c;//!<!matchedjet all tracks fake + real
    TH3D* h_MJ0_e3c;//!<!matchedjet 0 real tracks
    TH3D* h_MJ1_e3c;//!<!matchedjet 1 real track
    TH3D* h_MJ2_e3c;//!<!matchedjet 2 real tracks
    TH3D* h_MJ3_e3c;//!<!matchedjet all real tracks

    TH3D* h_MJ_e3c_tru;//!<!matchedjet all tracks fake + real true
    TH3D* h_MJ0_e3c_tru;//!<!matchedjet 0 real tracks true
    TH3D* h_MJ1_e3c_tru;//!<!matchedjet 1 real track true
    TH3D* h_MJ2_e3c_tru;//!<!matchedjet 2 real tracks true
    TH3D* h_MJ3_e3c_tru;//!<!matchedjet all real tracks true

    TH3D* h_MB1MB1MB1;//!<!bkg-bkg-bkg correlations
    TH3D* h_MB1MB1MB1_tru;//!<!bkg-bkg-bkg correlations true

    TH3D* h_JJMB;//!<!jet-jet-bkg1 correlations
    TH3D* h_JJMB_tru;//!<!jet-jet-bkg1 correlations true
    TH3D* h_JMBMB;//!<!jet-bkg1-bkg1 correlations
    TH3D* h_JMBMB_tru;//!<!jet-bkg1-bkg1 correlations true
    TH3D* h_MB1MB1MB2;//!<!bkg1-bkg1-bkg2 correlations
    TH3D* h_MB1MB1MB2_tru;//!<!bkg1-bkg1-bkg2 correlations true
    TH3D* h_MB1MB2MB2;//!<!bkg1-bkg2-bkg2 correlations
    TH3D* h_MB1MB2MB2_tru;//!<!bkg1-bkg2-bkg2 correlations true

    TH3D* h_JMB1MB2;//!<!jet-bkg1-bkg2 correlations
    TH3D* h_JMB1MB2_tru;//!<!jet-bkg1-bkg2 correlations true
    TH3D* h_MB1MB2MB3;//!<!bkg1-bkg2-bk3 correlations
    TH3D* h_MB1MB2MB3_tru;//!<!bkg1-bkg2-bkg3 correlations true


    TH3D* h_BMBMB;//!<!bkg-bkg1-bkg1 correlations
    TH3D* h_SMBMB;//!<!sig-bkg-bkg1 correlations

    TH3D* h_BBMB;//!<!bkg-bkg-bkg1 correlations
    TH3D* h_BMB1MB2;//!<!bkg-bkg1-bkg2 correlations
    TH3D* h_SMB1MB2;//!<!sig-bkg1-bkg2 correlations
    TH3D* h_SBMB;//!<!sig-bkg-bkg1 correlations
    TH3D* h_SSMB;//!<!sig-sig-bkg1 correlations

    TH3D* h_BBMB_tru;//!<!bkg-bkg-bkg1 correlations true
    TH3D* h_SBMB_tru;//!<!sig-bkg-bkg1 correlations true
    TH3D* h_SSMB_tru;//!<!sig-sig-bkg1 correlations true
    TH3D* h_BMBMB_tru;//!<!bkg-bkg1-bkg1 correlations true
    TH3D* h_SMBMB_tru;//!<!sig-bkg1-bkg1 correlations true
    TH3D* h_BMB1MB2_tru;//!<!bkg-bkg1-bkg2 correlations true
    TH3D* h_SMB1MB2_tru;//!<!sig-bkg1-bkg2 correlations true


    TTree *fTreeMatchTracks; ///< Tree with matched tracks from MC
    TTree *fTreeData; ///< Tree with tracks from Data

    TString             fMCParticleArrayName; ///< Array name of MC particles in event (mcparticles)
    TClonesArray*       fMCParticleArray;         //!<! Array of MC particles in event (usually mcparticles)
    TRandom3* fRandom; //!<! random number generator
    
private:
    AliAnalysisTaskJetsEECpbpb(
                           const AliAnalysisTaskJetsEECpbpb &); // not implemented
    AliAnalysisTaskJetsEECpbpb &
    operator=(const AliAnalysisTaskJetsEECpbpb &); // not implemented
    
    ClassDef(AliAnalysisTaskJetsEECpbpb, 16) //change this to 15 if you add something new
};
#endif

