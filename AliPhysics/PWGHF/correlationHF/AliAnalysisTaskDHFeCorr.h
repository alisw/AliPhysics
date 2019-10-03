/* *
 * \class AliAnalysisTaskDHFeCorr
 *
 * \brief D meson - heavy-flavour electron correlations
 *
 * This class performs the data filtering for the D meson - HFe correlations.
 * The task provides the same interface to analyse D0 mesons and D+. The D* is planned to be implemented.
 * The electrons are identified and only events with electrons are saved. The non-HF identification is also performed
 * in this task.
 *
 * The configuration is done using a yaml file. The following parameters can be adjusted:
 *
 * -# Event selection: corresponds to the yaml parameter 'event'. All the selected events use the trigger kINT7 for
 * the moment.
 *   - Vertex Z: corresponds to the yaml parameter 'zvxt_range'. The value should of the format [min,max].
 *   - Multiplicity: corresponds to the yaml parameters 'estimator' and 'multiplicity_range'. The value of 'estimator'
 *   should be also defined added to AliAnalysisTaskDHFeCorr::fgkMultiplicityEstimators. The values allowed are set
 * using 'multiplicity_range' in the format [min,max].
 *
 * -# D-meson selection: corresponds to the yaml parameter 'D_meson' and the sub parameter 'selection'. The cuts are
 * implemented in the standard D-meson selection framework, with the AliRDHFCuts as the base class.
 *   - Root file containing the cuts: should be specified in the field 'cut_file'. The name of the object that has the
 *   cuts should be provided in 'cut_name'.
 *   - Minimum transverse momentum: the task rejects low momentum D mesons to reduce the size of the output. Please set
 *   it with 'min_pt'.
 *   - Maximum transverse momentum used for particle identification: particle identification (PID) is used to reduce
 *   the size of the output. The 'max_pt_pid' value can be used to set the maximum value which PID will be used. You can
 *   set the PID selection as usual in the cut file and it will be ignored for pT > max_pt_pid.
 *
 * -# Electron selection is divided in 'main' electrons, the electrons which will be used to the correlation, and partner
 * electrons, which will are the ones used to identify non-HF electrons. The partner electrons usually have looser cuts
 * to achieve a high non-HF tagging efficiency. They are defined in the yaml file as 'main_electron' and
 * 'partner_electron'. The selection for both of them share the same parameters. They are:
 *   - Tracking, under 'track' sub parameter
 *     - Transverse momentum range: set 'pt_range' with [min, max]
 *     - Eta: set 'eta_range' with [min, max]
 *     - AOD filter bit: set with 'filter_bit'. the supported type should cover any declared with AliAODTrack and
 *     it should be also declared one of the types declared in fgkAODFilterBitMap
 *     - Minimum number of cluster in the TPC: set with 'min_TPC_cls'
 *     - Minimum number of cluster in the TPC used for the calculation of the dE/dx: set with 'min_TPC_cls_dedx'
 *     - Minimum number of clusters in the ITS: set with 'min_ITS_hits'
 *     - Requirement in the ITS SPD: set with 'pixel_req'. Possible values can be found at
 *     AliAnalysisTaskDHFeCorr::ITSPixel_t
 *     - The DCA maximum values can be set using 'dca_z' and 'dca_xy'
 *   - Particle identification, under 'PID' sub parameter
 *     - TPC N sigma selection: set the maximun and minimum with 'TPC_selection'
 *     - TOF N sigma selection: set the maximun and minimum with 'TOF_selection'. It can be turned on and off by
 *     'require_TOF'.
 */
#ifndef AliAnalysisTaskDHFeCorr_H
#define AliAnalysisTaskDHFeCorr_H

//ROOT includes
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

//AliROOT includes
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliPID.h"

//AliPhysics includes
#include "AliYAMLConfiguration.h"
#include "AliRDHFCuts.h"
#include "AliAODRecoDecayHF.h"

class TClonesArray;
#include <vector>
#include <memory>
#include <string>


namespace AliDHFeCorr {
    //structure to hold the D meson or electron during the selection process in the task
    typedef struct AliEvent{
        UInt_t fGridPID{0};
        UInt_t fEventNumber{0};
        Float_t fCentrality{-1.};
        Float_t fVtxZ{-999.};
    } AliEvent;
    
    typedef struct AliDMeson {
        AliAODRecoDecayHF *fRecoObj{nullptr};
        Int_t fGridPID{0}; ///<PID of the grid job used to create the tree
        Int_t fEventNumber{0}; ///< Number of the event
        UInt_t fID{0}; ///< D meson id in the event
        Bool_t fIsParticleCandidate{kFALSE}; ///< Particle hypotheses at reconstruction level
        
        //Basic Information
        Float_t fPt{-999.};
        Float_t fEta{-999.};
        Float_t fPhi{-999.};
        Float_t fY{-999.};
        Float_t fInvMass{0.};
        Float_t fReducedChi2{-999.};
        
        //Topologic information
        Float_t fDecayLength{-999.};
        Float_t fDecayLengthXY{-999.};
        
        Float_t fNormDecayLength{-999.};
        Float_t fNormDecayLengthXY{-999.};
        
        Float_t fCosP{-999.};
        Float_t fCosPXY{-999.};
        
        Float_t fImpParXY{-999.};
        Float_t fDCA{-999.};
        
        //D0 and D+, D*+ shared information
        Float_t fNormd0MeasMinusExp{-999.};
        UChar_t fSelectionStatusDefaultPID{99}; ///< Default PID selection status
        
        //D0 information (also used by the D*, since it has a D0)
        Float_t fCosTs{-999.};
        
        //D+ information
        Float_t fSigmaVertex{-999.};
        
        //D* information
        Float_t fAngleD0dkpPisoft{-999.};
        
        //Single-track information
        std::vector<Float_t> fPtDaughters;
        std::vector<Float_t> fD0Daughters;
        std::vector<UInt_t> fIDDaughters; ///< ID obtained using GetID()
        
        std::vector<Float_t> fNSigmaTPCDaughters; ///< The PID TPC response (n sigma)
        std::vector<Float_t> fNSigmaTOFDaughters; ///< The PID TOF response (n sigma)

        //MC Level information
        Float_t fPtMC{-999.};///< Transverse momentum (MC information)
        Bool_t fIsD{kFALSE}; // Is it signal or background? (MC information)
        Bool_t fIsParticle{kFALSE}; //is it a D0 or D0bar? (MC information)
        Bool_t fIsPrompt{kFALSE}; // Does it comes from charm? (MC information)
    } AliDMeson;

    typedef struct AliElectron {
    public:
        AliAODTrack *fTrack{nullptr};
        
        UInt_t fGridPID{0};
        UInt_t fEventNumber{0};
        
        UInt_t fID{0};
        Char_t fCharge{0};
        Float_t fPt{-999.};
        Float_t fP{-999.};
        Float_t fEta{-999.};
        Float_t fPhi{-999.};
        
        UShort_t fNClsTPC{999};
        UShort_t fNClsTPCDeDx{999};
        UChar_t fNITSCls{99};
        
        Bool_t fITSHitFirstLayer{kFALSE};
        Bool_t fITSHitSecondLayer{kFALSE};
        Float_t fDCAxy{-999.};
        Float_t fDCAz{-999.};
        
        Float_t fTPCNSigma{-999.};
        Float_t fTOFNSigma{-999.};
        
        std::vector<Float_t> fInvMassPartnersULS; //mass of the ULS partners
        std::vector<Float_t> fInvMassPartnersLS; //mass of the LS partners
        std::vector<UInt_t> fPartnersULSID; //unique ID of the ULS partners
        std::vector<UInt_t> fPartnersLSID; //unique ID of the LS partners
        
        //MC information
        Float_t fPtMC{-999.};
        Float_t fPhiMC{-999.};
        Float_t fEtaMC{-999.};
        Char_t fOrigin{0}; // track origin from AliVertexingHFUtils::CheckOrigin
        Int_t fPDG{0};
        
        Int_t fFirstMotherPDG{0};
        Float_t fFirstMotherPt{-999.};
        
        Int_t fSecondMotherPDG{0};
        Float_t fSecondMotherPt{-999.};
        
    } AliElectron;
    
    //Structs to hold the selection values ("cuts")
    typedef struct AliEventSelection {
        Float_t fVertexZMin{-10.};
        Float_t fVertexZMax{10.};

        Float_t fMultMin{0.};
        Float_t fMultMax{101.};
        std::string fMultiEstimator = "VOM";
    } AliEventSelection;

    typedef struct AliElectronSelection {
        //Track selection
        Float_t fPtMin{-999.};
        Float_t fPtMax{999.};
        Float_t fEtaMin{-1.};
        Float_t fEtaMax{1.};

        AliAODTrack::AODTrkFilterBits_t fFilterBit{AliAODTrack::kTrkGlobalNoDCA};

        Int_t fTPCClsMin{0};
        Int_t fTPCClsDeDxMin{0};
        Int_t fITSClsMin{0};

        Int_t fITSPixel{0};

        Bool_t fRecalculateDCA{kTRUE};
        Float_t fDCAz{999.};
        Float_t fDCAxy{999.};

        //PID selection
        Bool_t fRequireTOF{kFALSE};
        Float_t fTOFNSigmaMin{-3.};
        Float_t fTOFNSigmaMax{3.};

        Float_t fTPCNSigmaMin{-1};
        Float_t fTPCNSigmaMax{3.};
    } AliElectronSelection;

    typedef struct AliPhotonSelection {
        Float_t fInvMassMax{0.3};
        Float_t fReducedChi2Max{3.};
    } AliPhotonSelection;

    typedef struct AliDMesonSelection {
        std::unique_ptr<AliRDHFCuts> fDMesonCuts;
        Float_t fPtMin{-999.};
        Float_t fPtMaxPID{999.};
    } AliDMesonSelection;

    //Structs to hold QA plots
    typedef struct AliEventQAHistograms {
        std::unique_ptr<TH1F> fNEvents;
        std::unique_ptr<TH1F> fVertexZ;
        std::unique_ptr<TH1F> fCentrality;

    } AliEventQAHistograms;

    typedef struct AliElectronQAHistograms {
        //Tracking
        std::unique_ptr<TH3F> fPtEtaPhi;
        std::unique_ptr<TH1F> fTPCNCls;
        std::unique_ptr<TH1F> fTPCNClsDeDx;
        std::unique_ptr<TH1F> fITSCls;
        std::unique_ptr<TH1F> fDCAz;
        std::unique_ptr<TH1F> fDCAxy;

        //PID
        std::unique_ptr<TH2F> fTPCNsigmaPt;
        std::unique_ptr<TH2F> fTPCNsigmaP;
        std::unique_ptr<TH2F> fTOFNsigmaPt;
        std::unique_ptr<TH2F> fTOFNsigmaP;

        std::unique_ptr<TH2F> fTPCTOFNSigma;
    } AliElectronQAHistograms;

    typedef struct AliDMesonQAHistos {
        std::unique_ptr<TH3F> fPtEtaPhi;
        std::unique_ptr<TH2F> fPtInvMass;
    } AliDMesonQAHistos;

    //Structs to hold configurations for output
    typedef struct AliConfigureElectronOpt {
        std::vector<Float_t> fPtBins;
        std::vector<Float_t> fPBins;

        Int_t fNBinsEta{20};
        Int_t fNBinsPhi{40};
        Int_t fNBinsTPCCls{160};
        Int_t fNBinsNsigma{100};
    } AliConfigureElectronOpt;

    typedef struct AliConfigureDMesonOpt {
        std::vector<Float_t> fPtBins;
        Float_t fInvMassMin{1.3};
        Float_t fInvMassMax{2.5};
        Int_t fNBinsInvMass{100};

        Int_t fNBinsPhi{40};
        Int_t fNBinsEta{16};
    } AliConfigureDMesonOpt;

    namespace Utils {
        template<typename T>
        bool GetPropertyRange(const PWG::Tools::AliYAMLConfiguration &yaml_config,
                              const std::vector<std::string> property_path, T &property_min, T &property_max,
                              const bool required_property) {

            std::vector<T> range;

            if (!yaml_config.GetProperty(property_path, range, required_property)) {
                std::cout << "Check Config file. It is not possible to set " << property_path[property_path.size() - 1]
                          << std::endl;
                return false;
            }

            try {
                property_min = range[0];
                property_max = range[1];
            }
            catch (std::exception &exp) {
                std::cout << "Problem to read the " << property_path[property_path.size() - 1]
                          << "variable. The following error was raised: " << exp.what() << std::endl;
                return false;
            }
            return true;
        }
    }
}

class AliAnalysisTaskDHFeCorr : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskDHFeCorr();

    //D-meson species
    typedef enum {
        kD0, //D0
        kDplus, // D+
        kDstar //D*
    } DMeson_t; ///< Basic type used to identify the D mesons
    
    typedef enum {
        kNotSelected = 0,
        kParticle = 1,
        kAntiParticle = 2,
        kAmbiguous = 3 //Selected for both
    } SelectionStatus_t;

    //Electron ITS pixels
    typedef enum {
        kFirst = 0, //at least one on first
        kSecond = 1, // at least one on second
        kBoth = 2,
        kNone = 3,
        kAny = 4, //have at least one hit on layer 1 or 2
        kExclusiveSecond = 5, // only on layer 2
        kExclusiveFirst = 6 // only on layer 1
    } ITSPixel_t;
    
    explicit AliAnalysisTaskDHFeCorr(const char *name);

    virtual ~AliAnalysisTaskDHFeCorr() = default;

    static AliAnalysisTaskDHFeCorr *AddTask(const std::string &name = "DHFeCorrelation",
                                            const std::string &config_file = "$ALICE_PHYSICS/PWGHF/correlationHF/macros/default_config_d_hfe.yaml");

    virtual void UserCreateOutputObjects();

    virtual void UserExec(Option_t *option);

    virtual void Terminate(Option_t *option);

    //D-meson information
    const std::map<std::string, DMeson_t> fgkDMesonNames{
            {"D0",    AliAnalysisTaskDHFeCorr::kD0},
            {"Dplus", AliAnalysisTaskDHFeCorr::kDplus},
            {"Dstar", AliAnalysisTaskDHFeCorr::kDstar}
    }; ///< D meson names

    const std::map<DMeson_t, std::string> fgkDMesonListName{
            {AliAnalysisTaskDHFeCorr::kD0,    "D0toKpi"},
            {AliAnalysisTaskDHFeCorr::kDplus, "Charm3Prong"},
            {AliAnalysisTaskDHFeCorr::kDstar, "Dstar"}
    }; ///< Name of the list (TClonesArray*) containing the pre-filtered mesons

    const std::map<DMeson_t, Int_t> fgkDMesonPDG{
            {AliAnalysisTaskDHFeCorr::kD0,    421},
            {AliAnalysisTaskDHFeCorr::kDplus, 411},
            {AliAnalysisTaskDHFeCorr::kDstar, 413}
    };///< PDG number of each meson

    const std::map<DMeson_t, std::vector<Int_t>> fgkDMesonDaughterPDG{
            {AliAnalysisTaskDHFeCorr::kD0,    {211, 321}},
            {AliAnalysisTaskDHFeCorr::kDplus, {321, 211, 211}},
            {AliAnalysisTaskDHFeCorr::kDstar, {321, 321, 211}}
    }; ////< PDG of the D-meson daughters

    //Electron information
    const std::map<std::string, AliAODTrack::AODTrkFilterBits_t> fgkAODFilterBitMap = {
            {"kTrkTPCOnly",            AliAODTrack::kTrkTPCOnly},
            {"kTrkITSsa",              AliAODTrack::kTrkITSsa},
            {"kTrkITSConstrained",     AliAODTrack::kTrkITSConstrained},
            {"kTrkElectronsPID",       AliAODTrack::kTrkElectronsPID},
            {"kTrkGlobalNoDCA",        AliAODTrack::kTrkGlobalNoDCA},
            {"kTrkGlobal",             AliAODTrack::kTrkGlobal},
            {"kTrkGlobalSDD",          AliAODTrack::kTrkGlobalSDD},
            {"kTrkTPCOnlyConstrained", AliAODTrack::kTrkTPCOnlyConstrained}
    }; ///< AOD filter bits supported

    const std::map<std::string, ITSPixel_t> fgkITSPixelMap = {
            {"kFirst",           AliAnalysisTaskDHFeCorr::kFirst},
            {"kSecond",          AliAnalysisTaskDHFeCorr::kSecond},
            {"kBoth",            AliAnalysisTaskDHFeCorr::kBoth},
            {"kNone",            AliAnalysisTaskDHFeCorr::kNone},
            {"kAny",             AliAnalysisTaskDHFeCorr::kAny},
            {"kExclusiveSecond", AliAnalysisTaskDHFeCorr::kExclusiveSecond},
            {"kExclusiveFirst",  AliAnalysisTaskDHFeCorr::kExclusiveFirst}
    }; ///< ITS pixel 

    const std::map<DMeson_t, std::vector<AliPID::EParticleType> > fgkDMesonDaughterAliPID = {
            {AliAnalysisTaskDHFeCorr::kD0,    {AliPID::kPion, AliPID::kKaon}},
            {AliAnalysisTaskDHFeCorr::kDplus, {AliPID::kPion, AliPID::kKaon, AliPID::kPion}},
            {AliAnalysisTaskDHFeCorr::kDstar, {AliPID::kPion, AliPID::kKaon, AliPID::kPion}}
    };

    //Event information
    const std::vector<std::string> fgkMultiplicityEstimators = {"ZNA", "V0A"}; ////< Multiplicity estimator supported

    bool Configure(std::string config_file, std::string config_name = "custom_config",
                   std::string default_file = "$ALICE_PHYSICS/PWGHF/correlationHF/macros/default_config_d_hfe.yaml");

private:
    TList fOptEvent; ///< List with histograms from the event QA
    TList fOptElectron; ///< List with with histograms from the electron QA
    TList fOptDMeson; ///< List with with histograms from the D meson QA
    
    std::unique_ptr<TTree> fEventTree; ///< Tree with the event information
    std::unique_ptr<TTree> fElectronTree; ///< Tree with the electron information
    std::unique_ptr<TTree> fDmesonTree; ///< Tree with the D meson information

    UInt_t fGridPID{0}; ///< ID of the process on GRID
    UInt_t fEventNumber{0}; ///< Unique number for each event
    Float_t fVtxZ{-999.}; ///< Vertex Z
    Float_t fCentrality{-1.}; ///< Centrality

    AliDHFeCorr::AliElectron fElectron; ///< Electron information that will be used in the fElectronTree
    AliDHFeCorr::AliDMeson fDmeson;///< D meson information that will be used in the fDmesonTree

    bool fIsMC{false}; ///< Flag for MC analysis
    bool fIsEffMode{false}; ///< Fills only the trees with true D mesons and electrons. No correlation analysis.
    bool fKeepAllCandidates{true}; ///< Keep all candidates (even if no electron is present)
    DMeson_t fDmesonSpecies{kD0}; ///< The D meson species (D0, D+ or D*)

    //YAML configuration
    PWG::Tools::AliYAMLConfiguration fYAMLConfig; //<< YAML config object, used to set the parameters of the task

    bool ConfigureEventSelection(const std::string &name, AliDHFeCorr::AliEventSelection &event_selection);
    bool ConfigureElectrons(const std::string &name, AliDHFeCorr::AliElectronSelection &electron_selection);
    bool ConfigureElectronOpt(const std::string &name, AliDHFeCorr::AliConfigureElectronOpt &opt_config,
                              const std::string &base = "output");
    bool ConfigureDMesonOpt(const std::string &name, AliDHFeCorr::AliConfigureDMesonOpt &opt_config,
                            const std::string &base = "output");
    bool ConfigureDMesons(const std::string &name, AliDHFeCorr::AliDMesonSelection &opt_config);
    
    void CheckConfiguration() const;
    
    void AddEventVariables(std::unique_ptr<TTree> &tree);
    void AddElectronVariables(std::unique_ptr<TTree> &tree);
    void AddDMesonVariables(std::unique_ptr<TTree> &tree, DMeson_t meson_species);
    void AddElectronMCVariables(std::unique_ptr<TTree> &tree);
    void AddDMesonMCVariables(std::unique_ptr<TTree> &tree);

    AliDHFeCorr::AliEventSelection fEventSelection; ///<< Event selection configuration
    
    AliDHFeCorr::AliElectronSelection fElectronRequirements; ///<< Electron selection configuration
    AliDHFeCorr::AliElectronSelection fPartnerElectronRequirements; ///<< Partner selection configuration
    AliDHFeCorr::AliConfigureElectronOpt fElectronOptConfig; ///<< Electron output (histograms) configuration

    AliDHFeCorr::AliPhotonSelection fPhotonSelection; //<< Photon (electron pair) selection configuration

    AliDHFeCorr::AliDMesonSelection fDMesonRequirements; ///<< Electron selection configuration 
    AliDHFeCorr::AliConfigureDMesonOpt fDMesonOptConfig; ///<< Electron output (histograms) configuration

    //QA histograms
    AliDHFeCorr::AliEventQAHistograms fEventQABeforeCuts; //! QA histograms for events before the selection
    AliDHFeCorr::AliEventQAHistograms fEventQAAfterCuts; //! QA histograms for events after the selection

    AliDHFeCorr::AliElectronQAHistograms fElectronQABeforeCuts; //! QA histograms for electrons before the selection
    AliDHFeCorr::AliElectronQAHistograms fElectronQAAfterTrackCuts; //! QA histograms for electrons before the selection
    AliDHFeCorr::AliElectronQAHistograms fElectronQAAfterCuts;  //! QA histograms for electrons after the selection

    AliDHFeCorr::AliDMesonQAHistos fDMesonQABeforeCuts; //! QA histograms for D mesons before the selection
    AliDHFeCorr::AliDMesonQAHistos fDMesonQAAfterCuts;  //! QA histograms for D mesons before the selection

    //Functions for the QA
    static AliDHFeCorr::AliEventQAHistograms CreateQAEvents(const std::string &stage, TList &opt_list);

    static AliDHFeCorr::AliElectronQAHistograms CreateQAElectrons(AliDHFeCorr::AliConfigureElectronOpt config,
                                                                  const std::string &type_particle,
                                                                  const std::string &stage, TList &opt_list);

    static AliDHFeCorr::AliDMesonQAHistos
    CreateQADMeson(AliDHFeCorr::AliConfigureDMesonOpt config, const std::string &type_particle,
                   const std::string &stage, TList &opt_list);

    void FillEventQA(AliDHFeCorr::AliEventQAHistograms &histograms);

    void FillElectronQA(const std::vector<AliDHFeCorr::AliElectron> &tracks,
                        AliDHFeCorr::AliElectronQAHistograms &histograms);

    static void FillDmesonQA(const std::vector<AliDHFeCorr::AliDMeson> &d_mesons,
                             AliDHFeCorr::AliDMesonQAHistos &histogram, DMeson_t meson_species);

    void SetGridPID();
    void PostOutput(); //< Function to post the data

    static float GetMaxd0MeasMinusExp(AliAODRecoDecayHF *candidate, float b_field);
    static std::vector<Float_t> MakeBins(Float_t start, Float_t end, Int_t n_bins);
    
    //Analysis steps
    std::vector<AliDHFeCorr::AliElectron> ElectronAnalysis();
    std::vector<AliDHFeCorr::AliDMeson> DMesonAnalysis();

    //Select Events
    bool IsSelectedEvent(AliAODEvent *event, const AliDHFeCorr::AliEventSelection &selection);

    //Filter Electrons, Partner Electrons
    static bool FulfilPixelSelection(const AliDHFeCorr::AliElectron &track, Int_t requirement);
    void FillElectronInformation(std::vector<AliDHFeCorr::AliElectron> &electrons);
    
    std::vector<AliDHFeCorr::AliElectron>
    FilterElectronsTracking(AliDHFeCorr::AliElectronSelection electronSelection,
                            const std::vector<AliDHFeCorr::AliElectron> &electrons);

    std::vector<AliDHFeCorr::AliElectron> FilterElectronsPID(AliDHFeCorr::AliElectronSelection electronSelection,
                                                             const std::vector<AliDHFeCorr::AliElectron> &electrons);
    void FillElectronMCInfo(std::vector<AliDHFeCorr::AliElectron> &electrons);

    //Filter D mesons
    std::vector<AliDHFeCorr::AliDMeson>
    FillDMesonInfo(const TClonesArray *dmeson_candidates, DMeson_t meson_species,
                        AliDHFeCorr::AliDMesonSelection &dmeson_selection) const;

    void FillDmesonMCInfo(std::vector<AliDHFeCorr::AliDMeson> &d_mesons, AliAnalysisTaskDHFeCorr::DMeson_t meson_species) const;

    static std::vector<AliDHFeCorr::AliDMeson> FilterTrueDMesons(const std::vector<AliDHFeCorr::AliDMeson> &d_mesons);

    static std::vector<AliDHFeCorr::AliDMeson>
    FilterDmesons(const std::vector<AliDHFeCorr::AliDMeson> &dmeson_candidates,
                  const AliDHFeCorr::AliDMesonSelection &selectionDMeson, AliAODEvent *aod_event, int pdg_dmeson);

    //Select ULS and LS pairs
    void FindNonHFe(AliDHFeCorr::AliElectron &main_electron, const std::vector<AliDHFeCorr::AliElectron> &partners) const;

    //Save to tree
    void FillElectronTree(const std::vector<AliDHFeCorr::AliElectron> &tracks);
    void FillDMesonTree(const std::vector<AliDHFeCorr::AliDMeson> &d_mesons);
    
ClassDef(AliAnalysisTaskDHFeCorr, 2);
    
};

#endif
