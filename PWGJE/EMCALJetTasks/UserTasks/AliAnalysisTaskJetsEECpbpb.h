#ifndef ALIANALYSISTASKJetsEECPBPB_H
#define ALIANALYSISTASKJetsEECPBPB_H

class TH1;
class TH2;
class TH3;
class TH3F;
class TTree;
class THnSparse;
class TClonesArray;
class TObjArray;
class TArrayI;
class AliAnalysisManager;
class AliJetContainer;
class AliAODTrack;
class AliAODEvent;
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
    void SetCentralitySelectionOn(Bool_t t) { fCentSelectOn = t; }
    void SetOneConstSelectionOn(Bool_t t) { fOneConstSelectOn = t; }
    void SetCheckTracksOn(Bool_t t) { fTrackCheckPlots = t; }
    void SetCheckResolution(Bool_t t) { fCheckResolution = t; }
    void SetMinPtConst(Float_t t) { fMinPtConst = t;}
    void SetMinTrackPtEncs(Float_t t) { fMinENCtrackPt = t;}
    void SetHardCutoff(Float_t t) { fHardCutoff = t; }
    void SetDoTwoTrack(Bool_t t) { fDoTwoTrack = t; }
    void SetCutDoubleCounts(Bool_t t) {fCutDoubleCounts = t;}
    void SetPowerAlgorithm(Float_t t) { fPowerAlgo = t; }
    void SetMinCentrality(Float_t t) { fCentMin = t; }
    void SetMaxCentrality(Float_t t) { fCentMax = t; }
    void SetDerivativeSubtractionOrder(Int_t c) { fDerivSubtrOrder = c; }
    void SetDetLevelJetsOn(Bool_t t) { fStoreDetLevelJets = t; }
    void SetMatchRadius(Float_t t) { fMatchR = t; } 
    void SetMinJetPtSub(Float_t t) { fjetMinPtSub = t; } 
    void SetStoreTrig(Bool_t t) {fStoreTrig = t;}
    void SetConeRadius(Float_t t) { fConeR = t; } 
    
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
    void SetEmbTruthJetPtMin(Float_t t) {fEmbTuthJetPtMin = t;}

    void SetDoEmbedding(Bool_t t) { fDoEmbedding = t; }
    void SetDoPerpCone(Bool_t t) { fDoPerpCone = t; }
    void SetDoRandCone(Bool_t t) {fDoRandCone = t;}
    void SetCoutStatements(Bool_t t) {fCout = t;}

    void SetTruthLabelRange(Int_t min, Int_t max) { fTruthMinLabel = min; fTruthMaxLabel = max; }
    void SetJetMatchingSharedPtFraction(Double_t val) { fJetMatchingSharedPtFraction = val; }
    void SetEEC(Bool_t t) { ifeec = t; }
    void SetE3C(Bool_t t) { ife3c = t; }
    void SetMinPtHist(Bool_t t) { ifMinPtHist = t; }
    void SetcFactorHist(Bool_t t) { ifcFactorHist = t; }
    void SetAddEventCuts(Bool_t t) { fAddEventCuts = t; }
    void SetHighPtTrackCutEvent(Float_t t) { fHighPtTrackCutEvent = t; }
    void SetAxisShift(Float_t t) { fDeltaAxisShift = t; }

    
    
protected:
    Bool_t RetrieveEventObjects();
    Bool_t Run();
    Bool_t FillHistograms();
    Bool_t CountEmptyEvents();  // Just count if there is empty events
    
    int GetConstituentID(int constituentIndex, const AliVParticle* part, AliEmcalJet * jet);
    double delR(const fastjet::PseudoJet& ps1, const fastjet::PseudoJet& ps2);
  
    Double_t GetDownscaleWeight(string tstring);
    void ComputeENC(AliEmcalJet *fJet, double ptSub, AliJetContainer *fJetCont);
    
    void ComputeEncMC(AliEmcalJet *fJet, AliJetContainer *fJetCont, AliEmcalJet *fJet_tru, Int_t km); //MC and det correlations
    
    // void ComputeJetTrackMatch(AliEmcalJet *fJet, AliJetContainer *fJetCont, AliEmcalJet *fJet_tru, AliJetContainer *fJet_truCont, Int_t km);
    
    void FillEmbJetsEEC(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2, double jetpt, double pt, bool typeSame, std::string type, double bkgIndex1, double bkgIndex2, bool cfactor, std::string tag);

    void FillEmbJetsE3C(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2, std::vector<fastjet::PseudoJet> particles3, double jetpt, double pt, std::string typeSame, std::string type, bool ifMatchedJet);

    std::vector<fastjet::PseudoJet> FindThermalConeEEC(AliEmcalJet *fJetEmb, std::string axisType, int index);

    std::tuple<std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>> FindThermalConeE3C(AliEmcalJet *fJetEmb, int index1, int index2, int index3);

    void FillBkgSubJetsDataEEC(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2, double jetpt, bool typeSame, std::string type);
    void FillBkgSubJetsDataE3C(std::vector<fastjet::PseudoJet> particles, std::vector<fastjet::PseudoJet> particles2, std::vector<fastjet::PseudoJet> particles3, double jetpt, std::string typeSame, std::string type);

    std::vector<fastjet::PseudoJet> FindConesDataEEC(AliEmcalJet *fJet,std::string axisType);
    std::tuple<std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>> FindConesDataE3C(AliEmcalJet *fJet);

    bool PerformGeometricalJetMatching(AliJetContainer& contBase, AliJetContainer& contTag, double maxDist);
    void DoJetMatching();
    void GetTrueJetPtFraction(AliEmcalJet* jet, Double_t& truePtFraction, Double_t& truePtFraction_mcparticles);
    void GetMatchedJetObservables(AliEmcalJet* jet, Double_t& detJetPt, Double_t& partJetPt, Double_t& detJetPhi, Double_t& detJetEta, Double_t& partJetPhi, Double_t& partJetEta, Double_t& detJetDistance, Double_t& partJetDistance, Float_t& JetEmbPtSub);
    Long64_t GetUniqueEventID(AliVEvent* inputEvent);
    Bool_t CheckHighPtTrackInEvent(AliVEvent* inputEvent);

    Double_t  GetDistance(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2)
    {
    Double_t deltaPhi = TMath::Min(TMath::Abs(phi1-phi2),TMath::TwoPi() - TMath::Abs(phi1-phi2));
    return TMath::Sqrt((eta1-eta2)*(eta1-eta2) + deltaPhi*deltaPhi);
    }

    void FillMatchedTrackTree(AliEmcalJet* fJet_Hyb, AliEmcalJet* fJet_Tru, float_t& JetEmbPtSub);
    void FillEmbHistograms(Long64_t EvCounter);

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
    Float_t fjetMaxPtSub; ///<max subtracted jet pT 
    Float_t fjetMinArea; ///<min area of subtracted jet 
    Float_t fEmbTuthJetPtMin; ///<min pT of truth embedded jet
    
    Bool_t fStoreTrig; ///<storing the trigger class
    Float_t fConeR; ///<cone radius for perp or random cones
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
    Bool_t fisGoodIncEvent; ///<check if event is good 

    Float_t fJet_pt_det; ///< jet pt det 
    Float_t fJet_pt_tru; ///< jet pt tru
    Float_t fTrack_pt_det; ///< track pt det 
    Float_t fTrack_pt_tru; ///< track pt tru
    Float_t fTrack_eta_det; ///< track eta det 
    Float_t fTrack_eta_tru; ///< track eta tru 
    Float_t fTrack_phi_det; ///< track phi det 
    Float_t fTrack_phi_tru; ///< track phi tru
    Float_t fTrack_pt_miss; ///< miss track pt 
    Float_t fTrack_pt_fake; ///< fake track pt 
    Float_t fTrack_eta_miss; ///< miss track eta 
    Float_t fTrack_eta_fake; ///< fake track eta  
    Float_t fTrack_phi_miss; ///< miss track phi  
    Float_t fTrack_phi_fake; ///< fake track phi

    Float_t fJet_pt_dat; ///< jet pt dat
    Float_t fTrack_pt_dat; ///< track pt dat
    Float_t fTrack_eta_dat; ///< track eta dat
    Float_t fTrack_phi_dat; ///< track phi dat

    //From jet extractor task
    Bool_t              fDoPartLevelMatching;
    Bool_t              fDoDetLevelMatching;
    Int_t               fTruthMinLabel;           ///< min track label to consider it as true particle
    Int_t               fTruthMaxLabel;           ///< max track label to consider it as true particle
    Bool_t              fSaveMCInformation;       ///< save MC information
    Double_t            fJetMatchingSharedPtFraction; ///< Shared pT fraction required in matching
    Float_t             fjetRhoVal;
    Float_t             fjetRecoRhoVal;
    Float_t             fjetEmbRhoVal;
    Float_t             fjetGenRhoVal;
    Bool_t              fhasAcceptedEMCjet;
    
    AliJetContainer*   fJetContainer;   //!<! signal jet container
    AliJetContainer*   fbgJetContainer;   //!<! background jet container
    AliJetContainer*   fEmbJetContainer;   //!<! signal embedded jet container
    AliJetContainer*   fRecoJetContainer;   //!<! signal reco jet container
    AliJetContainer*   fGenJetContainer;   //!<! signal gen jet container
    
    //Adding for tracks at the particle level (not jet level)
    TString fGeneratorLevelName;
    TString fDetectorLevelName;
    
    AliMCParticleContainer* fGeneratorLevel;
    AliTrackContainer* fDetectorLevel;
    AliJetContainer* fJet_truCont;

    AliAODEvent* fAOD;//!<! AOD object

    //Histograms
    TH1D *jet_pt_hist; //!<! initializing histogram with jet pt
    
    TH1D *EEC_hist; //!<! initializing histogram with correlators
    TH2F *EEC_pt_hist; //!<! initializing 2D histogram with correlators and jet pt
 
    TH1D *E3C_hist; //!<! initializing histogram with correlators
    TH2F *E3C_pt_hist; //!<! initializing 2D histogram with correlators and jet pt
    
    //For corrections
    TH3F *EEC_det_pt_hist_3d; //!<! initializing 3d histogram with det level EEC and det and true jet pt
    TH3F *EEC_tru_pt_hist_3d; //!<! initializing 3d histogram with tru level EEC and det and true jet pt
    
    TH3F *E3C_det_pt_hist_3d; //!<! initializing 3d histogram with det level E3C and det and true jet pt
    TH3F *E3C_tru_pt_hist_3d; //!<! initializing 3d histogram with tru level E3C and det and true jet pt
    
    TH3F *N2_det_pt_hist_3d; //!<! initializing 3d histogram with det level EEC num pairs and det and true jet pt
    TH3F *N2_tru_pt_hist_3d; //!<! initializing 3d histogram with tru level EEC num pairs and det and true jet pt
    
    TH3F *N3_det_pt_hist_3d;  //!<! initializing 3d histogram with det level E3C num pairs and det and true jet pt
    TH3F *N3_tru_pt_hist_3d; //!<! initializing 3d histogram with tru level E3C num pairs and det and true jet pt
    
    TH2F *EEC_det_match_pt_det; //!<! initializing 2D histogram of det level EEC and det jet pt
    TH2F *EEC_tru_match_pt_tru; //!<! initializing 2D histogram of tru level EEC and tru jet pt
    
    TH2F *E3C_det_match_pt_det; //!<! initializing 2D histogram of det level E3C and det jet pt
    TH2F *E3C_tru_match_pt_tru; //!<! initializing 2D histogram of tru level E3C and tru jet pt
    
    TH1D *pt_tru; //!<! all tru level
    TH1D *pt_tru_match; //!<! tru level matched
    TH1D *pt_det; //!<! all det level
    TH1D *pt_det_match; //!<! det level matched
    
    TH1D *test_hist; //!<! test histogram
    
    TH2F *R_matrix; //!<! response matrix for jet pt
    
    TH2F *JES; //!<! jet energy scale
    TH2F *JES_scaled; //!<! jet energy scale for scaled pt_det
    TH2F *JER; //!<! jet energy resolution
    
    TH3F *pair_det_EEC;//!<! histogram for computing pair efficiency at det level
    TH3F *pair_tru_EEC;//!<! histogram for computing pair efficiency at truth level
    
    TH3F *pair_det_E3C;//!<! histogram for computing pair efficiency at det level
    TH3F *pair_tru_E3C;//!<! histogram for computing pair efficiency at truth level
    
    TH2F *qpt_tru;//!<! q/pt vs R histogram for input to fast sim
    TH2F *qpt_det;//!<! q/pt vs R histogram for input to fast sim
    TH1D *track_pt_tru; //!<! all tru level tracks
    TH1D *track_pt_det; //!<! all det level tracks
    TH2F *track_pt_matched; //!<! all det level tracks matched with truth level
    
    TH3F *R_match_eec; //!<! pair distance of all det level tracks in det jet matched with truth level tracks in truth jet eec as a function of jet pT
    TH3F *wt_match_eec; //!<! pair wts of all det level tracks in det jet matched with truth level tracks in truth jet eec as a function of jet pT
    TH3F *R_match_e3c; //!<! pair distance of all det level tracks in det jet matched with truth level tracks in truth jet e3c as a function of jet pT
    TH3F *wt_match_e3c; //!<! pair wts of all det level tracks in det jet matched with truth level tracks in truth jet e3c as a function of jet pT
    
    TH1D *qpt_tru1;//!<! q/pt vs R histogram for input to fast sim
    TH1D *qpt_tru2;//!<! q/pt vs R histogram for input to fast sim
    //!
    TH1D *pt_tru1;//!<! q/pt vs R histogram for input to fast sim
    TH1D *pt_tru2;//!<! q/pt vs R histogram for input to fast sim
    
    TH2F *eec_Mm; //!<! initializing 2D histogram of EEC matched and missed
    TH2F *eec_mm; //!<! initializing 2D histogram of tru level EEC missed and missed
    TH2F *e3c_MMm; //!<! initializing 2D histogram of det level E3C matched matched missed
    TH2F *e3c_Mmm; //!<! initializing 2D histogram of tru level E3C  matched missed missed
    TH2F *e3c_mmm; //!<! initializing 2D histogram of det level E3C missed missed missed
    
    TH2F *eec_Mf; //!<! initializing 2D histogram of det level EEC matched and fake
    TH2F *eec_ff; //!<! initializing 2D histogram of tru level EEC fake and fake
    TH2F *e3c_MMf; //!<! initializing 2D histogram of det level E3C matched matched fake
    TH2F *e3c_Mff; //!<! initializing 2D histogram of tru level E3C mateched fake fake
    TH2F *e3c_fff; //!<! initializing 2D histogram of det level E3C fake fake fake

    TH2F *eec_matched_det; //!<! initializing 2D histogram of det level EEC and det jet pt
    TH2F *eec_matched_tru; //!<! initializing 2D histogram of tru level EEC and tru jet pt
    TH2F *e3c_matched_det; //!<! initializing 2D histogram of det level E3C and det jet pt
    TH2F *e3c_matched_tru; //!<! initializing 2D histogram of tru level E3C and tru jet pt
    
    TH3F *wt_res_eec; //!<! initializing 2D histogram for wt resolution EEC
    TH3F *wt_res_e3c; //!<! initializing 2D histogram for wt resolution E3C
    TH3F *R_res_eec; //!<! initializing 2D histogram for R resolution EEC
    TH3F *R_res_e3c; //!<! initializing 2D histogram for R resolution E3C
    TH3F *wtnojet_match_eec;//!<! initializing 2D histogram for EEC wt without jet pt
    TH3F *wtnojet_match_e3c;//!<! initializing 2D histogram for E3C wt without jet pt
    TH3F *wtnojet_res_eec; //!<! initializing 2D histogram for wt resolution with jet pt excluded EEC
    TH3F *wtnojet_res_e3c; //!<! initializing 2D histogram for wt resolution with jet pt excluded E3C
    
    TH3F *R_match_eec_tru; //!<! pair distance of all det level tracks in det jet matched with truth level tracks in truth jet eec as a function of jet pT tru
    TH3F *wt_match_eec_tru; //!<! pair wts of all det level tracks in det jet matched with truth level tracks in truth jet eec as a function of jet pT tru
    TH3F *R_match_e3c_tru; //!<! pair distance of all det level tracks in det jet matched with truth level tracks in truth jet e3c as a function of jet pT tru
    TH3F *wt_match_e3c_tru; //!<! pair wts of all det level tracks in det jet matched with truth level tracks in truth jet e3c as a function of jet pT tru

    TH3F *wt_res_eec_tru; //!<! initializing 3D histogram for wt resolution EEC as a function of jet pT tru
    TH3F *wt_res_e3c_tru; //!<! initializing 3D histogram for wt resolution E3C as a function of jet pT tru
    TH3F *R_res_eec_tru; //!<! initializing 3D histogram for R resolution EEC as a function of jet pT tru
    TH3F *R_res_e3c_tru; //!<! initializing 3D histogram for R resolution E3C as a function of jet pT tru
    TH3F *wtnojet_match_eec_tru;//!<! initializing 3D histogram for EEC wt without jet pt as a function of jet pT tru
    TH3F *wtnojet_match_e3c_tru;//!<! initializing 3D histogram for E3C wt without jet pt as a function of jet pT tru
    TH3F *wtnojet_res_eec_tru; //!<! initializing 3D histogram for wt resolution with jet pt excluded EEC as a function of jet pT tru
    TH3F *wtnojet_res_e3c_tru; //!<! initializing 3D histogram for wt resolution with jet pt excluded E3C as a function of jet pT tru

    TH1D *constituentId; //!<! initializing 1D histogram for debugging
    TH3F *R_res_eec_tru_debug; //!<! initializing 3D histogram for R resolution EEC as a function of jet pT tru debug
    TH2F *track_pt_res_debug; //!<! initializing 3D histogram for pt resolution debug
    TH3F *track_eta_debug; //!<! initializing 3D histogram for eta match debug
    TH3F *track_rap_debug; //!<! initializing 3D histogram for rap match debug
    TH2F *track_phi_debug; //!<! initializing 3D histogram phi match debug
    TH2F *track_R_debug; //!<! initializing 3D histogram for R match debug
    TH2F *track_R_debug_rap; //!<! initializing 3D histogram for R-rap debug
    
    TH3F *track_eta_res_debug; //!<! initializing 3D histogram for eta resolution of single matched track as a function of jet pT tru debug
    TH3F *track_rap_res_debug; //!<! initializing 3D histogram for rap resolution of single matched track as a function of jet pT tru debug
    TH3F *track_phi_res_debug; //!<! initializing 3D histogram for phi resolution of single matched track as a function of jet pT tru debug
    TH3F *track_R_res_debug; //!<! initializing 3D histogram for R resolution of single matched track as a function of jet pT tru debug
    TH3F *track_R_rap_res_debug; //!<! initializing 3D histogram for R-rap resolution of single matched track as a function of jet pT tru debug

    TH3F *R_match_eec_tru_rap; //!<! initializing 3D histogram for R resolution of single matched track as a function of jet pT tru debug

    TH2F *track_pt_response_debug; //!<! initializing 2D histogram for debug pt
    TH2F *track_pt_wt_res_debug; //!<! initializing 2D histogram for debug wt
    TH2F *track_pt_wt_response_debug; //!<! initializing 2D histogram for debug wt_response
    
    TH3F* jetpt_res_w_R; //!<! initializing 3D histogram for R resolution of single matched track as a function of jet pT tru debug
    TH3F* jetpt_res_w_wt; //!<! initializing 3D histogram for weight resolution of single matched track as a function of jet pT tru debug
    
    TH1F* NumJetEvent; //!<! number of jet events

//For subtraction in data EEC 
//data cases have "_dat" in them 
 TH2F* h_JMB_dat;//!<!jet-bkg correlations data 2Dhist
 TH2F* h_MB1_dat;//!<!bkg-bkg correlations data 2Dhist
 TH2F* h_MB1MB2_dat;//!<!bkg-bkg correlations data 2Dhist

 TH3F* h3Jet_deltaR_JMB_dat;//!<!jet-bkg correlations data
 TH3F* h3Jet_deltaR_MB1_dat;//!<!bkg1-bkg1 correlations data
 TH3F* h3Jet_deltaR_MB1MB2_dat;//!<!bkg1-bkg2 correlations data

//Create histograms for all matched ("_m") and unmatched ("_um") jets with ("_c") and without c-factor 
//True histograms are only for matched jets "_tru_m"
//For the minpt case there is no "_m" because it is understood that truth min pT cut can only be applied for matched jets

   
    TH3F* h_MJ;//!<!all jets all tracks fake + real
    TH3F* h_MJ0;//!<!all jets 0 real tracks
    TH3F* h_MJ1;//!<!all jets 1 real track
    TH3F* h_MJ2;//!<!all jets 2 real tracks

    TH3F* h_MJ_m;//!<!matchedjet all tracks fake + real
    TH3F* h_MJ0_m;//!<!matchedjet 0 real tracks
    TH3F* h_MJ1_m;//!<!matchedjet 1 real track
    TH3F* h_MJ2_m;//!<!matchedjet 2 real tracks

    TH3F* h_MJ_tru_m;//!<!matchedjet all tracks fake + real tru pt
    TH3F* h_MJ0_tru_m;//!<!matchedjet 0 real tracks tru pt
    TH3F* h_MJ1_tru_m;//!<!matchedjet 1 real track tru pt
    TH3F* h_MJ2_tru_m;//!<!matchedjet 2 real tracks tru pt

    TH3F* h_MJ_um;//!<!unmatchedjet all tracks fake + real
    TH3F* h_MJ0_um;//!<!unmatchedjet 0 real tracks
    TH3F* h_MJ1_um;//!<!unmatchedjet 1 real track
    TH3F* h_MJ2_um;//!<!unmatchedjet 2 real tracks

    TH3F* h_MJ_c;//!<!all jets all tracks fake + real c-factor 
    TH3F* h_MJ0_c;//!<!all jets all tracks fake + real c-factor 
    TH3F* h_MJ1_c;//!<!all jets all tracks fake + real c-factor 
    TH3F* h_MJ2_c;//!<!all jets all tracks fake + real c-factor 

    TH3F* h_MJ_c_m;//!<!matchedjet all tracks fake + real c-factor
    TH3F* h_MJ0_c_m;//!<!matchedjet 0 real tracks c-factor 
    TH3F* h_MJ1_c_m;//!<!matchedjet 1 real track c-factor 
    TH3F* h_MJ2_c_m;//!<!matchedjet 2 real tracks c-factor 
   
    TH3F* h_MJ_tru_c_m;//!<!matchedjet all tracks fake + real c-factor tru pt
    TH3F* h_MJ0_tru_c_m;//!<!matchedjet 0 real tracks c-factor tru pt
    TH3F* h_MJ1_tru_c_m;//!<!matchedjet 1 real track c-factor tru pt
    TH3F* h_MJ2_tru_c_m;//!<!matchedjet 2 real tracks c-factor tru pt

    TH3F* h_MJ_c_um;//!<!unmatchedjet all tracks fake + real c-factor
    TH3F* h_MJ0_c_um;//!<!unmatchedjet 0 real tracks c-factor 
    TH3F* h_MJ1_c_um;//!<!unmatchedjet 1 real track c-factor 
    TH3F* h_MJ2_c_um;//!<!unmatchedjet 2 real tracks c-factor 

    TH3F* h_MB1;//!<!all bkg-bkg correlations
    TH3F* h_MB1_m;//!<!matchedjet bkg-bkg correlations
    TH3F* h_MB1_tru_m;//!<!matchedjet bkg-bkg correlations tru pt
    TH3F* h_MB1_um;//!<!unmatchedjet bkg-bkg correlations
    TH3F* h_MB1_c;//!<!alljet bkg-bkg correlations c-factor
    TH3F* h_MB1_c_m;//!<!matchedjet bkg-bkg correlations c-factor
    TH3F* h_MB1_tru_c_m;//!<!matchedjet bkg-bkg correlations c-factor
    TH3F* h_MB1_c_um;//!<!unmatchedjet bkg-bkg correlations c-factor

    TH3F* h_JMB;//!<!all jet-bkg correlations
    TH3F* h_JMB_m;//!<!matchedjet jet-bkg correlations
    TH3F* h_JMB_tru_m;//!<!matchedjet jet-bkg correlations tru pt
    TH3F* h_JMB_um;//!<!unmatchedjet jet-bkg correlations
    TH3F* h_JMB_c;//!<!alljet jet-bkg correlations c-factor
    TH3F* h_JMB_c_m;//!<!matchedjet jet-bkg correlations c-factor
    TH3F* h_JMB_tru_c_m;//!<!matchedjet jet-bkg correlations c-factor
    TH3F* h_JMB_c_um;//!<!unmatchedjet jet-bkg correlations c-factor

    TH3F* h_SMB;//!<!all signal-bkg correlations
    TH3F* h_SMB_m;//!<!matchedjet signal-bkg correlations
    TH3F* h_SMB_tru_m;//!<!matchedjet signal-bkg correlations tru pt
    TH3F* h_SMB_um;//!<!unmatchedjet signal-bkg correlations
    TH3F* h_SMB_c;//!<!alljet signal-bkg correlations c-factor
    TH3F* h_SMB_c_m;//!<!matchedjet signal-bkg correlations c-factor
    TH3F* h_SMB_tru_c_m;//!<!matchedjet signal-bkg correlations c-factor
    TH3F* h_SMB_c_um;//!<!unmatchedjet signal-bkg correlations c-factor

    TH3F* h_BMB;//!<!all jetbkg-bkg correlations
    TH3F* h_BMB_m;//!<!matchedjet jetbkg-bkg correlations
    TH3F* h_BMB_tru_m;//!<!matchedjet jetbkg-bkg correlations tru pt
    TH3F* h_BMB_um;//!<!unmatchedjet jetbkg-bkg correlations
    TH3F* h_BMB_c;//!<!alljet jetbkg-bkg correlations c-factor
    TH3F* h_BMB_c_m;//!<!matchedjet jetbkg-bkg correlations c-factor
    TH3F* h_BMB_tru_c_m;//!<!matchedjet jetbkg-bkg correlations c-factor
    TH3F* h_BMB_c_um;//!<!unmatchedjet jetbkg-bkg correlations c-factor

    TH3F* h_MB1MB2;//!<!all bkg1-bkg2 correlations
    TH3F* h_MB1MB2_m;//!<!matchedjet bkg1-bkg2 correlations
    TH3F* h_MB1MB2_tru_m;//!<!matchedjet bkg1-bkg2 correlations tru pt
    TH3F* h_MB1MB2_um;//!<!unmatchedjet bkg1-bkg2 correlations
    TH3F* h_MB1MB2_c;//!<!alljet bkg1-bkg2 correlations c-factor
    TH3F* h_MB1MB2_c_m;//!<!matchedjet bkg1-bkg2 correlations c-factor
    TH3F* h_MB1MB2_tru_c_m;//!<!matchedjet bkg1-bkg2 correlations c-factor
    TH3F* h_MB1MB2_c_um;//!<!unmatchedjet bkg1-bkg2 correlations c-factor

   //Declarations for E3C 
    TH3F* h_MJ_e3c;//!<!matchedjet all tracks fake + real
    TH3F* h_MJ0_e3c;//!<!matchedjet 0 real tracks
    TH3F* h_MJ1_e3c;//!<!matchedjet 1 real track
    TH3F* h_MJ2_e3c;//!<!matchedjet 2 real tracks
    TH3F* h_MJ3_e3c;//!<!matchedjet all real tracks

    TH3F* h_MJ_e3c_tru;//!<!matchedjet all tracks fake + real true
    TH3F* h_MJ0_e3c_tru;//!<!matchedjet 0 real tracks true
    TH3F* h_MJ1_e3c_tru;//!<!matchedjet 1 real track true
    TH3F* h_MJ2_e3c_tru;//!<!matchedjet 2 real tracks true
    TH3F* h_MJ3_e3c_tru;//!<!matchedjet all real tracks true

    TH3F* h_MB1MB1MB1;//!<!bkg-bkg-bkg correlations
    TH2F* h_MB1MB1MB1_dat;//!<!bkg-bkg-bkg correlations data
    TH3F* h_MB1MB1MB1_tru;//!<!bkg-bkg-bkg correlations true

    TH3F* h_JJMB;//!<!jet-jet-bkg1 correlations
    TH2F* h_JJMB_dat;//!<!jet-jet-bkg1 correlations data
    TH3F* h_JJMB_tru;//!<!jet-jet-bkg1 correlations true
    TH3F* h_JMBMB;//!<!jet-bkg1-bkg1 correlations
    TH2F* h_JMBMB_dat;//!<!jet-bkg1-bkg1 correlations data
    TH3F* h_JMBMB_tru;//!<!jet-bkg1-bkg1 correlations true
    TH3F* h_MB1MB1MB2;//!<!bkg1-bkg1-bkg2 correlations
    TH2F* h_MB1MB1MB2_dat;//!<!bkg1-bkg1-bkg2 correlations data
    TH3F* h_MB1MB1MB2_tru;//!<!bkg1-bkg1-bkg2 correlations true
    TH3F* h_MB1MB2MB2;//!<!bkg1-bkg2-bkg2 correlations
    TH2F* h_MB1MB2MB2_dat;//!<!bkg1-bkg1-bkg2 correlations data
    TH3F* h_MB1MB2MB2_tru;//!<!bkg1-bkg2-bkg2 correlations true

    TH3F* h_JMB1MB2;//!<!jet-bkg1-bkg2 correlations
    TH2F* h_JMB1MB2_dat;//!<!jet-bkg1-bkg2 correlations data
    TH3F* h_JMB1MB2_tru;//!<!jet-bkg1-bkg2 correlations true
    TH3F* h_MB1MB2MB3;//!<!bkg1-bkg2-bk3 correlations
    TH2F* h_MB1MB2MB3_dat;//!<!bkg1-bkg2-bk3 correlations data
    TH3F* h_MB1MB2MB3_tru;//!<!bkg1-bkg2-bkg3 correlations true


    TH3F* h_BMBMB;//!<!bkg-bkg1-bkg1 correlations
    TH3F* h_SMBMB;//!<!sig-bkg-bkg1 correlations

    TH3F* h_BBMB;//!<!bkg-bkg-bkg1 correlations
    TH3F* h_BMB1MB2;//!<!bkg-bkg1-bkg2 correlations
    TH3F* h_SMB1MB2;//!<!sig-bkg1-bkg2 correlations
    TH3F* h_SBMB;//!<!sig-bkg-bkg1 correlations
    TH3F* h_SSMB;//!<!sig-sig-bkg1 correlations

    TH3F* h_BBMB_tru;//!<!bkg-bkg-bkg1 correlations true
    TH3F* h_SBMB_tru;//!<!sig-bkg-bkg1 correlations true
    TH3F* h_SSMB_tru;//!<!sig-sig-bkg1 correlations true
    TH3F* h_BMBMB_tru;//!<!bkg-bkg1-bkg1 correlations true
    TH3F* h_SMBMB_tru;//!<!sig-bkg1-bkg1 correlations true
    TH3F* h_BMB1MB2_tru;//!<!bkg-bkg1-bkg2 correlations true
    TH3F* h_SMB1MB2_tru;//!<!sig-bkg1-bkg2 correlations true

//Histograms that will be used for subtracting and unfolding (jetpT, RL, weight)
    TH3F* h3Jet_deltaR_MJ;//!<!all jets all jet correlations 
    TH3F* h3Jet_deltaR_MJ0;//!<!all jets all fake tracks
    TH3F* h3Jet_deltaR_MJ1;//!<!all jets one fake track
    TH3F* h3Jet_deltaR_MJ2;//!<!all jets all matched tracks

    TH3F* h3Jet_deltaR_MJ_m;//!<!matchedjets all jet correlations 
    TH3F* h3Jet_deltaR_MJ0_m;//!<!matchedjets all fake tracks
    TH3F* h3Jet_deltaR_MJ1_m;//!<!matchedjets one fake track
    TH3F* h3Jet_deltaR_MJ2_m;//!<!matchedjets all matched tracks

    TH3F* h3Jet_deltaR_MJ_tru_m;//!<!matchedjets all jet correlations tru pt
    TH3F* h3Jet_deltaR_MJ0_tru_m;//!<!matchedjets all fake tracks tru pt
    TH3F* h3Jet_deltaR_MJ1_tru_m;//!<!matchedjets one fake track tru pt 
    TH3F* h3Jet_deltaR_MJ2_tru_m;//!<!matchedjets all matched tracks tru pt 

    TH3F* h3Jet_deltaR_MJ_minpt;//!<!matchedjets all jet correlations with truth pT min cut
    TH3F* h3Jet_deltaR_MJ0_minpt;//!<!matchedjets all fake tracks with truth pT min cut
    TH3F* h3Jet_deltaR_MJ1_minpt;//!<!matchedjets one fake track with truth pT min cut
    TH3F* h3Jet_deltaR_MJ2_minpt;//!<!matchedjets all matched tracks with truth pT min cut

    TH3F* h3Jet_deltaR_MJ_tru_minpt;//!<!matchedjets all correlations with true jet pt with truth pT min cut
    TH3F* h3Jet_deltaR_MJ0_tru_minpt;//!<!matchedjets fake tracks with true jet pt with truth pT min cut
    TH3F* h3Jet_deltaR_MJ1_tru_minpt;//!<!matchedjets one fake track with true jet pt with truth pT min cut
    TH3F* h3Jet_deltaR_MJ2_tru_minpt;//!<!matchedjets all matched tracks with true jet pt with truth pT min cut

    TH3F* h3Jet_deltaR_MJ_um;//!<!unmatchedjets all jet correlations 
    TH3F* h3Jet_deltaR_MJ0_um;//!<!unmatchedjets all fake tracks
    TH3F* h3Jet_deltaR_MJ1_um;//!<!unmatchedjets one fake track
    TH3F* h3Jet_deltaR_MJ2_um;//!<!unmatchedjets all matched tracks

    TH3F* h3Jet_deltaR_MJ_c;//!<!all jet correlations (for c-factor)
    TH3F* h3Jet_deltaR_MJ0_c;//!<!all fake tracks  (for c-factor)
    TH3F* h3Jet_deltaR_MJ1_c;//!<!one fake track  (for c-factor)
    TH3F* h3Jet_deltaR_MJ2_c;//!<!all matched tracks  (for c-factor)

    TH3F* h3Jet_deltaR_MJ_c_m;//!<!matchedjets all jet correlations (for c-factor)
    TH3F* h3Jet_deltaR_MJ0_c_m;//!<!matchedjets all fake tracks (for c-factor)
    TH3F* h3Jet_deltaR_MJ1_c_m;//!<!matchedjets one fake track (for c-factor)
    TH3F* h3Jet_deltaR_MJ2_c_m;//!<!matchedjets all matched tracks (for c-factor)

    TH3F* h3Jet_deltaR_MJ_tru_c_m;//!<!matchedjets all jet correlations with true jet pt (for c-factor)
    TH3F* h3Jet_deltaR_MJ0_tru_c_m;//!<!matchedjets all fake tracks with true jet pt (for c-factor)
    TH3F* h3Jet_deltaR_MJ1_tru_c_m;//!<!matchedjets one fake track with true jet pt (for c-factor)
    TH3F* h3Jet_deltaR_MJ2_tru_c_m;//!<!matchedjets all matched tracks with true jet pt (for c-factor)

    TH3F* h3Jet_deltaR_MJ_minpt_c;//!<!matchedjets all jet correlations with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_MJ0_minpt_c;//!<!matchedjets all fake tracks with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_MJ1_minpt_c;//!<!matchedjets one fake track with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_MJ2_minpt_c;//!<!matchedjets all matched tracks with truth pT min cut (for c-factor)

    TH3F* h3Jet_deltaR_MJ_tru_minpt_c;//!<!matchedjets all jet correlations with true jet pt with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_MJ0_tru_minpt_c;//!<!amatchedjets ll fake tracks with true jet pt with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_MJ1_tru_minpt_c;//!<!matchedjets one fake track with true jet pt with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_MJ2_tru_minpt_c;//!<!matchedjets all matched tracks with true jet pt with truth pT min cut (for c-factor)

    TH3F* h3Jet_deltaR_MJ_c_um;//!<!unmatchedjets all jet correlations (for c-factor)
    TH3F* h3Jet_deltaR_MJ0_c_um;//!<!unmatchedjets all fake tracks (for c-factor)
    TH3F* h3Jet_deltaR_MJ1_c_um;//!<!unmatchedjets one fake track (for c-factor)
    TH3F* h3Jet_deltaR_MJ2_c_um;//!<!unmatchedjets all matched tracks (for c-factor)

    TH3F* h3Jet_deltaR_MB1;//!<!all jets all bkg bkg correlations in one cone
    TH3F* h3Jet_deltaR_MB1_m;//!<!matchedjets all bkg bkg correlations in one cone
    TH3F* h3Jet_deltaR_MB1_tru_m;//!<!matchedjets all bkg bkg correlations in one cone tru pt
    TH3F* h3Jet_deltaR_MB1_minpt;//!<!matchedjets all bkg bkg correlations in one cone with truth pT min cut
    TH3F* h3Jet_deltaR_MB1_tru_minpt;//!<!matchedjets all bkg bkg correlations in one cone with true jet pt with truth pT min cut
    TH3F* h3Jet_deltaR_MB1_um;//!<!unmatchedjets all bkg bkg correlations in one cone
    TH3F* h3Jet_deltaR_MB1_c;//!<!all jets all bkg bkg correlations in one cones  (for c-factor)
    TH3F* h3Jet_deltaR_MB1_c_m;//!<!matchedjets all bkg bkg correlations in one cone (for c-factor)
    TH3F* h3Jet_deltaR_MB1_tru_c_m;//!<!matchedjets all bkg bkg correlations in one cone with true jet pt (for c-factor)
    TH3F* h3Jet_deltaR_MB1_minpt_c;//!<!matchedjets all bkg bkg correlations in one cone with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_MB1_tru_minpt_c;//!<!matchedjets all bkg bkg correlations in one cone true jet pt with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_MB1_c_um;//!<!unmatchedjets all bkg bkg correlations in one cone (for c-factor)

    TH3F* h3Jet_deltaR_JMB;//!<!all jets all jet-bkg correlations in one cone
    TH3F* h3Jet_deltaR_JMB_m;//!<!matchedjets all jet-bkg correlations in one cone
    TH3F* h3Jet_deltaR_JMB_tru_m;//!<!matchedjets all jet-bkg correlations in one cone tru pt
    TH3F* h3Jet_deltaR_JMB_minpt;//!<!matchedjets all jet-bkg correlations in one cone with truth pT min cut
    TH3F* h3Jet_deltaR_JMB_tru_minpt;//!<!matchedjets all jet-bkg correlations in one cone with true jet pt with truth pT min cut
    TH3F* h3Jet_deltaR_JMB_um;//!<!unmatchedjets all jet-bkg correlations in one cone
    TH3F* h3Jet_deltaR_JMB_c;//!<!all jets all jet-bkg correlations in one cones  (for c-factor)
    TH3F* h3Jet_deltaR_JMB_c_m;//!<!matchedjets all jet-bkg correlations in one cone (for c-factor)
    TH3F* h3Jet_deltaR_JMB_tru_c_m;//!<!matchedjets all jet-bkg correlations in one cone with true jet pt (for c-factor)
    TH3F* h3Jet_deltaR_JMB_minpt_c;//!<!matchedjets all jet-bkg correlations in one cone with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_JMB_tru_minpt_c;//!<!matchedjets all jet-bkg correlations in one cone true jet pt with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_JMB_c_um;//!<!unmatchedjets all jet-bkg correlations in one cone (for c-factor)

    TH3F* h3Jet_deltaR_SMB;//!<!all jets all signal-bkg correlations in one cone
    TH3F* h3Jet_deltaR_SMB_m;//!<!matchedjets all signal-bkg correlations in one cone
    TH3F* h3Jet_deltaR_SMB_tru_m;//!<!matchedjets all signal-bkg correlations in one cone tru pt
    TH3F* h3Jet_deltaR_SMB_minpt;//!<!matchedjets all signal-bkg correlations in one cone with truth pT min cut
    TH3F* h3Jet_deltaR_SMB_tru_minpt;//!<!matchedjets all signal-bkg correlations in one cone with true jet pt with truth pT min cut
    TH3F* h3Jet_deltaR_SMB_um;//!<!unmatchedjets all signal-bkg correlations in one cone
    TH3F* h3Jet_deltaR_SMB_c;//!<!all jets all signal-bkg correlations in one cones  (for c-factor)
    TH3F* h3Jet_deltaR_SMB_c_m;//!<!matchedjets all signal-bkg correlations in one cone (for c-factor)
    TH3F* h3Jet_deltaR_SMB_tru_c_m;//!<!matchedjets all signal-bkg correlations in one cone with true jet pt (for c-factor)
    TH3F* h3Jet_deltaR_SMB_minpt_c;//!<!matchedjets all signal-bkg correlations in one cone with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_SMB_tru_minpt_c;//!<!matchedjets all signal-bkg correlations in one cone true jet pt with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_SMB_c_um;//!<!unmatchedjets all signal-bkg correlations in one cone (for c-factor)

    TH3F* h3Jet_deltaR_BMB;//!<!all jets all bkg-bkg1 correlations in one cone
    TH3F* h3Jet_deltaR_BMB_m;//!<!matchedjets all bkg-bkg1 correlations in one cone
    TH3F* h3Jet_deltaR_BMB_tru_m;//!<!matchedjets all bkg-bkg1 correlations in one cone tru pt
    TH3F* h3Jet_deltaR_BMB_minpt;//!<!matchedjets all bkg-bkg1 correlations in one cone with truth pT min cut
    TH3F* h3Jet_deltaR_BMB_tru_minpt;//!<!matchedjets all bkg-bkg1 correlations in one cone with true jet pt with truth pT min cut
    TH3F* h3Jet_deltaR_BMB_um;//!<!unmatchedjets all bkg-bkg1 correlations in one cone
    TH3F* h3Jet_deltaR_BMB_c;//!<!all jets all bkg-bkg1 correlations in one cones  (for c-factor)
    TH3F* h3Jet_deltaR_BMB_c_m;//!<!matchedjets all bkg-bkg1 correlations in one cone (for c-factor)
    TH3F* h3Jet_deltaR_BMB_tru_c_m;//!<!matchedjets all bkg-bkg1 correlations in one cone with true jet pt (for c-factor)
    TH3F* h3Jet_deltaR_BMB_minpt_c;//!<!matchedjets all bkg-bkg1 correlations in one cone with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_BMB_tru_minpt_c;//!<!matchedjets all bkg-bkg1 correlations in one cone true jet pt with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_BMB_c_um;//!<!unmatchedjets all bkg-bkg1 correlations in one cone (for c-factor)

    TH3F* h3Jet_deltaR_MB1MB2;//!<!all jets all bkg1-bkg2 correlations in one cone
    TH3F* h3Jet_deltaR_MB1MB2_m;//!<!matchedjets all bkg1-bkg2 correlations in one cone
    TH3F* h3Jet_deltaR_MB1MB2_tru_m;//!<!matchedjets all bkg1-bkg2 correlations in one cone tru pt
    TH3F* h3Jet_deltaR_MB1MB2_minpt;//!<!matchedjets all bkg1-bkg2 correlations in one cone with truth pT min cut
    TH3F* h3Jet_deltaR_MB1MB2_tru_minpt;//!<!matchedjets all bkg1-bkg2 correlations in one cone with true jet pt with truth pT min cut
    TH3F* h3Jet_deltaR_MB1MB2_um;//!<!unmatchedjets all bkg1-bkg2 correlations in one cone
    TH3F* h3Jet_deltaR_MB1MB2_c;//!<!all jets all bkg1-bkg2 correlations in one cones  (for c-factor)
    TH3F* h3Jet_deltaR_MB1MB2_c_m;//!<!matchedjets all bkg1-bkg2 correlations in one cone (for c-factor)
    TH3F* h3Jet_deltaR_MB1MB2_tru_c_m;//!<!matchedjets all bkg1-bkg2 correlations in one cone with true jet pt (for c-factor)
    TH3F* h3Jet_deltaR_MB1MB2_minpt_c;//!<!matchedjets all bkg1-bkg2 correlations in one cone with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_MB1MB2_tru_minpt_c;//!<!matchedjets all bkg1-bkg2 correlations in one cone true jet pt with truth pT min cut (for c-factor)
    TH3F* h3Jet_deltaR_MB1MB2_c_um;//!<!unmatchedjets all bkg1-bkg2 correlations in one cone (for c-factor)

    TH3F* OptUn_eec;//!<!Histograms for optimizing unfolding bins eec
    TH3F* OptUn_e3c;//!<!Histograms for optimizing unfolding bins e3c

    TH1F* hJet_num_um;//!<!Histograms num of unmatched jets 
    TH1F* hJet_num_m;//!<!Histograms num of matched jets 
    TH1F* h_jetpt_m;//!<!Histograms matched jet pT
    TH1F* h_jetpt_um;//!<!Histograms unmatched jet pT
    TH2F* h_jetpt_m_tru;//!<!Histograms matched jet pT hyb vs true

    TH3F* h3_MB1MB1MB1_dat;//!<!Histograms for B'B'B' data 3D
    TH3F* h3_JJMB_dat;//!<!Histograms for JJB' data 3D
    TH3F* h3_JMB1MB2_dat;//!<!Histograms for JB'B'' data 3D
    TH3F* h3_JMBMB_dat;//!<!Histograms for JB'B' data 3D
    TH3F* h3_MB1MB1MB2_dat;//!<!Histograms for B'B'B'' data 3D
    TH3F* h3_MB1MB2MB2_dat;//!<!Histograms for B'B''B'' data 3D
    TH3F* h3_MB1MB2MB3_dat;//!<!Histograms for B'B''B''' data 3D

    TH3F* hJet_deltaR_MJ_e3c;//!<!Histograms for embedding 2D
    TH3F* hJet_deltaR_MJ0_e3c;//!<!Histograms for embedding 2D all fake tracks
    TH3F* hJet_deltaR_MJ1_e3c;//!<!Histograms for embedding 2D one real track
    TH3F* hJet_deltaR_MJ2_e3c;//!<!Histograms for embedding 2D two real tracks
    TH3F* hJet_deltaR_MJ3_e3c;//!<!Histograms for embedding 2D all real tracks
    

    TH3F* hJet_deltaR_MJ_e3c_m;//!<!Histograms for embedding 2D matchedjet
    TH3F* hJet_deltaR_MJ0_e3c_m;//!<!Histograms for embedding 2D all fake tracks mj
    TH3F* hJet_deltaR_MJ1_e3c_m;//!<!Histograms for embedding 2D one real tracks mj
    TH3F* hJet_deltaR_MJ2_e3c_m;//!<!Histograms for embedding 2D two real tracks mj
    TH3F* hJet_deltaR_MJ3_e3c_m;//!<!Histograms for embedding 2D all real tracks mj

    TH3F* h_MB1MB1MB1_m;//!<!Histograms for embedding 2D mbmbmb mj
    TH3F* h_JJMB_m;//!<!Histograms for embedding 2D jjmb mj
    TH3F* h_JMBMB_m;//!<!Histograms for embedding 2D jmbmb mj
    TH3F* h_MB1MB1MB2_m;//!<!Histograms for embedding 2D mb1mb1mb2 mj
    TH3F* h_MB1MB2MB2_m;//!<!Histograms for embedding 2D mb1mb2mb2 mj
    TH3F* h_JMB1MB2_m;//!<!Histograms for embedding 2D jmb1mb2 mj
    TH3F* h_MB1MB2MB3_m;//!<!Histograms for embedding 2D mb1mb2mb3 mj
    TH3F* h_BMBMB_m;//!<!Histograms for embedding 2D bmbmb mj
    TH3F* h_SMBMB_m;//!<!Histograms for embedding 2D smbmb mj
    TH3F* h_BMB1MB2_m;//!<!Histograms for embedding 2D bmb1mb2 mj
    TH3F* h_SMB1MB2_m;//!<!Histograms for embedding 2D smb1mb2 mj
    TH3F* h_BBMB_m;//!<!Histograms for embedding 2D bbmb mj
    TH3F* h_SBMB_m;//!<!Histograms for embedding 2D sbmb mj
    TH3F* h_SSMB_m;//!<!Histograms for embedding 2D ssmb mj

    TH3F* hJet_deltaR_MJ_e3c_tru_m;//!<!Histograms for embedding 2D mbmbmb mj tru
    TH3F* hJet_deltaR_MJ0_e3c_tru_m;//!<!Histograms for embedding 2D all fake tracks mj tru
    TH3F* hJet_deltaR_MJ1_e3c_tru_m;//!<!Histograms for embedding 2D one true tracks mj tru
    TH3F* hJet_deltaR_MJ2_e3c_tru_m;//!<!Histograms for embedding 2D two true tracks mj tru
    TH3F* hJet_deltaR_MJ3_e3c_tru_m;//!<!Histograms for embedding 2D all true tracks mj tru

    TH3F* h_MB1MB1MB1_tru_m;//!<!Histograms for embedding 2D mbmbmb mj tru
    TH3F* h_JJMB_tru_m;//!<!Histograms for embedding 2D jjmb mj tru
    TH3F* h_JMBMB_tru_m;//!<!Histograms for embedding 2D jmbmb mj tru
    TH3F* h_MB1MB1MB2_tru_m;//!<!Histograms for embedding 2D mb1mb1mb2 mj tru
    TH3F* h_MB1MB2MB2_tru_m;//!<!Histograms for embedding 2D mb1mb2mb2 mj tru
    TH3F* h_JMB1MB2_tru_m;//!<!Histograms for embedding 2D jmb1mb2 mj tru
    TH3F* h_MB1MB2MB3_tru_m;//!<!Histograms for embedding 2D mb1mb2mb3 mj tru
    TH3F* h_BMBMB_tru_m;//!<!Histograms for embedding 2D bmbmb mj tru
    TH3F* h_SMBMB_tru_m;//!<!Histograms for embedding 2D smbmb mj tru
    TH3F* h_BMB1MB2_tru_m;//!<!Histograms for embedding 2D bmb1mb2 mj tru
    TH3F* h_BBMB_tru_m;//!<!Histograms for embedding 2D bbmb mj tru
    TH3F* h_SBMB_tru_m;//!<!Histograms for embedding 2D sbmb tru
    TH3F* h_SSMB_tru_m;//!<!Histograms for embedding 2D ssmb tru
    TH3F* h_SMB1MB2_tru_m;//!<!Histograms for embedding 2D smb1mb2 mj tru

    TH3F* hJet_deltaR_MJ_e3c_um;//!<!Histograms for embedding 2D unmatched jet 
    TH3F* hJet_deltaR_MJ0_e3c_um;//!<!Histograms for embedding 2D unmatchedjet all bkg
    TH3F* hJet_deltaR_MJ1_e3c_um;//!<!Histograms for embedding 2D unmatchedjet 1 pythia
    TH3F* hJet_deltaR_MJ2_e3c_um;//!<!Histograms for embedding 2D unmatchedjet 2 pythia
    TH3F* hJet_deltaR_MJ3_e3c_um;//!<!Histograms for embedding 2D unmatchedjet 3 pythia

    TH3F* h_MB1MB1MB1_um;//!<!Histograms for embedding 2D mbmbmb umj
    TH3F* h_JJMB_um;//!<!Histograms for embedding 2D jjmb umj
    TH3F* h_JMBMB_um;//!<!Histograms for embedding 2D jmbmb umj
    TH3F* h_MB1MB1MB2_um;//!<!Histograms for embedding 2D mb1mb1mb2 umj
    TH3F* h_MB1MB2MB2_um;//!<!Histograms for embedding 2D mb1mb2mb2 umj
    TH3F* h_JMB1MB2_um;//!<!Histograms for embedding 2D jmb1mb2 umj
    TH3F* h_MB1MB2MB3_um;//!<!Histograms for embedding 2D mb1mb2mb3 umj
    TH3F* h_BMBMB_um;//!<!Histograms for embedding 2D bmbmb umj
    TH3F* h_SMBMB_um;//!<!Histograms for embedding 2D smbmb umj
    TH3F* h_BMB1MB2_um;//!<!Histograms for embedding 2D bmb1mb2 umj
    TH3F* h_SMB1MB2_um;//!<!Histograms for embedding 2D smb1mb2 umj
    TH3F* h_BBMB_um;//!<!Histograms for embedding 2D bbmb umj
    TH3F* h_SBMB_um;//!<!Histograms for embedding 2D sbmb umj
    TH3F* h_SSMB_um;//!<!Histograms for embedding 2D ssmb umj
    

    TH3F* h3Jet_deltaR_MJ_e3c;//!<!Histograms for embedding 3D
    TH3F* h3Jet_deltaR_MJ0_e3c;//!<!Histograms for embedding 3D all fake tracks
    TH3F* h3Jet_deltaR_MJ1_e3c;//!<!Histograms for embedding 3D one real track
    TH3F* h3Jet_deltaR_MJ2_e3c;//!<!Histograms for embedding 3D two real tracks
    TH3F* h3Jet_deltaR_MJ3_e3c;//!<!Histograms for embedding 3D all real tracks

    TH3F* h3_MB1MB1MB1;//!<!Histograms for embedding 3D mbmbmb
    TH3F* h3_JJMB;//!<!Histograms for embedding 3D jjmb
    TH3F* h3_JMBMB;//!<!Histograms for embedding 3D jmbmb
    TH3F* h3_MB1MB1MB2;//!<!Histograms for embedding 3D mb1mb1mb2
    TH3F* h3_MB1MB2MB2;//!<!Histograms for embedding 3D mb1mb2mb2
    TH3F* h3_JMB1MB2;//!<!Histograms for embedding 3D jmb1mb2
    TH3F* h3_MB1MB2MB3;//!<!Histograms for embedding 3D mb1mb2mb3
    TH3F* h3_BMBMB;//!<!Histograms for embedding 3D bmbmb
    TH3F* h3_SMBMB;//!<!Histograms for embedding 3D smbmb
    TH3F* h3_BMB1MB2;//!<!Histograms for embedding 3D bmb1mb2
    TH3F* h3_SMB1MB2;//!<!Histograms for embedding 3D smb1mb2
    TH3F* h3_BBMB;//!<!Histograms for embedding 3D bbmb
    TH3F* h3_SBMB;//!<!Histograms for embedding 3D sbmb
    TH3F* h3_SSMB;//!<!Histograms for embedding 3D ssmb

    TH3F* h3Jet_deltaR_MJ_e3c_m;//!<!Histograms for embedding 3D matchedjet
    TH3F* h3Jet_deltaR_MJ0_e3c_m;//!<!Histograms for embedding 3D all fake tracks mj
    TH3F* h3Jet_deltaR_MJ1_e3c_m;//!<!Histograms for embedding 3D one real tracks mj
    TH3F* h3Jet_deltaR_MJ2_e3c_m;//!<!Histograms for embedding 3D two real tracks mj
    TH3F* h3Jet_deltaR_MJ3_e3c_m;//!<!Histograms for embedding 3D all real tracks mj

    TH3F* h3_MB1MB1MB1_m;//!<!Histograms for embedding 3D mbmbmb mj
    TH3F* h3_JJMB_m;//!<!Histograms for embedding 3D jjmb mj
    TH3F* h3_JMBMB_m;//!<!Histograms for embedding 3D jmbmb mj
    TH3F* h3_MB1MB1MB2_m;//!<!Histograms for embedding 3D mb1mb1mb2 mj
    TH3F* h3_MB1MB2MB2_m;//!<!Histograms for embedding 3D mb1mb2mb2 mj
    TH3F* h3_JMB1MB2_m;//!<!Histograms for embedding 3D jmb1mb2 mj
    TH3F* h3_MB1MB2MB3_m;//!<!Histograms for embedding 3D mb1mb2mb3 mj
    TH3F* h3_BMBMB_m;//!<!Histograms for embedding 3D bmbmb mj
    TH3F* h3_SMBMB_m;//!<!Histograms for embedding 3D smbmb mj
    TH3F* h3_BMB1MB2_m;//!<!Histograms for embedding 3D bmb1mb2 mj
    TH3F* h3_SMB1MB2_m;//!<!Histograms for embedding 3D smb1mb2 mj
    TH3F* h3_BBMB_m;//!<!Histograms for embedding 3D bbmb mj
    TH3F* h3_SBMB_m;//!<!Histograms for embedding 3D sbmb mj
    TH3F* h3_SSMB_m;//!<!Histograms for embedding 3D ssmb mj

    TH3F* h3Jet_deltaR_MJ_e3c_tru_m;//!<!Histograms for embedding 3D mbmbmb mj tru
    TH3F* h3Jet_deltaR_MJ0_e3c_tru_m;//!<!Histograms for embedding 3D all fake tracks mj tru
    TH3F* h3Jet_deltaR_MJ1_e3c_tru_m;//!<!Histograms for embedding 3D one true tracks mj tru
    TH3F* h3Jet_deltaR_MJ2_e3c_tru_m;//!<!Histograms for embedding 3D two true tracks mj tru
    TH3F* h3Jet_deltaR_MJ3_e3c_tru_m;//!<!Histograms for embedding 3D all true tracks mj tru

    TH3F* h3_MB1MB1MB1_tru_m;//!<!Histograms for embedding 3D mbmbmb mj tru
    TH3F* h3_JJMB_tru_m;//!<!Histograms for embedding 3D jjmb mj tru
    TH3F* h3_JMBMB_tru_m;//!<!Histograms for embedding 3D jmbmb mj tru
    TH3F* h3_MB1MB1MB2_tru_m;//!<!Histograms for embedding 3D mb1mb1mb2 mj tru
    TH3F* h3_MB1MB2MB2_tru_m;//!<!Histograms for embedding 3D mb1mb2mb2 mj tru
    TH3F* h3_JMB1MB2_tru_m;//!<!Histograms for embedding 3D jmb1mb2 mj tru
    TH3F* h3_MB1MB2MB3_tru_m;//!<!Histograms for embedding 3D mb1mb2mb3 mj tru
    TH3F* h3_BMBMB_tru_m;//!<!Histograms for embedding 3D bmbmb mj tru
    TH3F* h3_SMBMB_tru_m;//!<!Histograms for embedding 3D smbmb mj tru
    TH3F* h3_BMB1MB2_tru_m;//!<!Histograms for embedding 3D bmb1mb2 mj tru
    TH3F* h3_BBMB_tru_m;//!<!Histograms for embedding 3D bbmb mj tru
    TH3F* h3_SBMB_tru_m;//!<!Histograms for embedding 3D sbmb tru
    TH3F* h3_SSMB_tru_m;//!<!Histograms for embedding 3D ssmb tru
    TH3F* h3_SMB1MB2_tru_m;//!<!Histograms for embedding 3D smb1mb2 mj tru

    TH3F* h3Jet_deltaR_MJ_e3c_um;//!<!Histograms for embedding 3D unmatched jet 
    TH3F* h3Jet_deltaR_MJ0_e3c_um;//!<!Histograms for embedding 3D unmatchedjet all bkg
    TH3F* h3Jet_deltaR_MJ1_e3c_um;//!<!Histograms for embedding 3D unmatchedjet 1 pythia
    TH3F* h3Jet_deltaR_MJ2_e3c_um;//!<!Histograms for embedding 3D unmatchedjet 2 pythia
    TH3F* h3Jet_deltaR_MJ3_e3c_um;//!<!Histograms for embedding 3D unmatchedjet 3 pythia

    TH3F* h3_MB1MB1MB1_um;//!<!Histograms for embedding 3D mbmbmb umj
    TH3F* h3_JJMB_um;//!<!Histograms for embedding 3D jjmb umj
    TH3F* h3_JMBMB_um;//!<!Histograms for embedding 3D jmbmb umj
    TH3F* h3_MB1MB1MB2_um;//!<!Histograms for embedding 3D mb1mb1mb2 umj
    TH3F* h3_MB1MB2MB2_um;//!<!Histograms for embedding 3D mb1mb2mb2 umj
    TH3F* h3_JMB1MB2_um;//!<!Histograms for embedding 3D jmb1mb2 umj
    TH3F* h3_MB1MB2MB3_um;//!<!Histograms for embedding 3D mb1mb2mb3 umj
    TH3F* h3_BMBMB_um;//!<!Histograms for embedding 3D bmbmb umj
    TH3F* h3_SMBMB_um;//!<!Histograms for embedding 3D smbmb umj
    TH3F* h3_BMB1MB2_um;//!<!Histograms for embedding 3D bmb1mb2 umj
    TH3F* h3_SMB1MB2_um;//!<!Histograms for embedding 3D smb1mb2 umj
    TH3F* h3_BBMB_um;//!<!Histograms for embedding 3D bbmb umj
    TH3F* h3_SBMB_um;//!<!Histograms for embedding 3D sbmb umj
    TH3F* h3_SSMB_um;//!<!Histograms for embedding 3D ssmb umj

    TTree *fTreeMatchTracks; ///< Tree with matched tracks from MC
    TTree *fTreeData; ///< Tree with tracks from Data

    TString             fMCParticleArrayName; ///< Array name of MC particles in event (mcparticles)
    TClonesArray*       fMCParticleArray;         //!<! Array of MC particles in event (usually mcparticles)
    TRandom3* fRandom; //!<! random number generator

    Bool_t ifeec;///< initialize eec histograms and run eec code
    Bool_t ife3c;///< initialize e3c histograms and run e3c code
    Bool_t ifMinPtHist;///< initialize minpthist and fill minpthistograms 
    Bool_t ifcFactorHist;///< initialize c_factor_hist and fill c_factor_histograms
    Bool_t fAddEventCuts;///< ignore embedding events with tracks with pT greater than fHighPtTrackCutEvent
    Float_t fHighPtTrackCutEvent;///< ignore embedding events with tracks with pT greater than fHighPtTrackCutEvent
    Float_t fDeltaAxisShift;///< shift cone axis in phi wrt jet 

    TH1F* h_dpt;//!<!Histograms for checking delta_pt 
    
private:
    AliAnalysisTaskJetsEECpbpb(
                           const AliAnalysisTaskJetsEECpbpb &); // not implemented
    AliAnalysisTaskJetsEECpbpb &
    operator=(const AliAnalysisTaskJetsEECpbpb &); // not implemented
    
    ClassDef(AliAnalysisTaskJetsEECpbpb, 34) //change this to 35 if you add something new
};
#endif

