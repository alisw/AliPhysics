#include <climits>

// This is purely to suppress warning inside ROOT once instanciating
// std::vector<bool>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include <TCollectionProxyInfo.h>
#pragma GCC diagnostic pop

#include <Compression.h>
#include <TFile.h>
#include <TSystem.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <AliVEvent.h>
#include <AliAODEvent.h>
#pragma GCC diagnostic pop
#include <AliESDEvent.h>
#include <AliESDHeader.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <AliVCaloCells.h>
#include <AliOADBContainer.h>
#include <AliAnalysisManager.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>
#include <AliAODMCParticle.h>
#include <AliESDMuonTrack.h>

#include <AliMultSelection.h>
#include <AliEventplane.h>

#include <AliMagF.h>

#include <AliAODTrack.h>
#include <AliAODVertex.h>
#include <AliAODTracklets.h>
#include <AliAODMCHeader.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include <AliAnalysisTaskEmcalEmbeddingHelper.h>
#pragma GCC diagnostic pop

#include "AliAnalysisTaskNTGJ.h"

#ifndef __CINT__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wdeprecated-register"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>

#include "eLut.cpp"
#include "half.cpp"

#include "keras_model.cc"

// #include "tiny_dnn/tiny_dnn.h"

#pragma GCC diagnostic pop
#include "special_function.h"
#include "mc_truth.h"
#include "emcal.h"
#include "jet.h"
#include "track_cuts.h"
#ifdef WITH_EFP7
#include "einstein_sum.h"
#include "efp7.cc"
#else // WITH_EFP7
#define FILL_EFP7 {}
#endif // WITH_EFP7
#include "isolation.h"
#endif // __CINT__

/////////////////////////////////////////////////////////////////////
// Replacement dlfcn.h for CINT

#ifndef RTLD_NOW
#define RTLD_NOW      0x00002
#endif // RTLD_NOW
#ifndef RTLD_GLOBAL
#define RTLD_GLOBAL   0x00100
#endif // RTLD_GLOBAL
#ifndef RTLD_NODELETE
#define RTLD_NODELETE 0x01000
#endif // RTLD_NODELETE

extern "C" {
    extern void *dlopen(const char *, int);
    extern void *dlsym(void *, const char *);
    extern char *dlerror(void);
}

/////////////////////////////////////////////////////////////////////

// FIXME: This needs to be moved to somewhere else, later on

namespace {

class alice_jec_t {
protected:
    const char *_version;
public:
    alice_jec_t(void)
        : _version("lhc16c2-2017-06-pre1")
    {
    }
    TLorentzVector jec_p(TLorentzVector p) const
    {
        const double pt = p.Pt();

        // Dummy value
        const double pt_calib = 2.0 * pt;

        const double eta = p.Eta();

        // Mean value for jet pT raw within 2.5-18 GeV
        static const double par_eta[] = { 0.015, 7.7, 0.65 };

        const double eta_calib = par_eta[0] *
                                 (TMath::Erf(par_eta[1] * (eta - par_eta[2])) +
                                  TMath::Erf(par_eta[1] * (eta + par_eta[2])));

        TLorentzVector calib;

        calib.SetPtEtaPhiM(pt_calib, eta_calib, p.Phi(), p.M());

        return calib;
    }
    const char *version(void) const
    {
        return _version;
    }
};

}

ClassImp(AliAnalysisTaskNTGJ);

#define EMCAL_GEOMETRY_NAME "EMCAL_COMPLETE12SMV1_DCAL_8SM"

#define ntrigger_class NTRIGGER_CLASS_MAX

#define B CHAR_MIN
#define b 0
#define S SHRT_MIN
#define s 0
#define I INT_MIN
#define i 0
#define F NAN
#define D NAN
#define L LONG_MIN
#define l 0
#define O false
#define BRANCH(b, t)                            \
    _branch_ ## b((t)),
#define BRANCH_ARRAY(b, d, t)
#define BRANCH_ARRAY2(b, d, e, t)
#define BRANCH_STR(b)
#define BRANCH_STR_ARRAY(b, d)                          \
    _branch_ ## b(std::vector<std::string>()),

#define CLASS_INITIALIZATION                                \
    _emcal_geometry_name(EMCAL_GEOMETRY_NAME),              \
_tree_event(NULL),                                      \
MEMBER_BRANCH                                           \
_run_number_current(INT_MIN),                           \
_f1_ncluster_tpc_linear_pt_dep(NULL),                   \
_track_cut(std::vector<AliESDtrackCuts>()),             \
_reco_util(new AliEMCALRecoUtils),                      \
_emcal_geometry(NULL),                                  \
_muon_track_cut(new AliMuonTrackCuts),                  \
_ncell(EMCAL_NCELL),                                    \
_emcal_geometry_filename(""),                           \
_emcal_local2master_filename(""),                       \
_force_ue_subtraction(false),                           \
_skim_cluster_min_e(-INFINITY),                         \
_skim_track_min_pt(-INFINITY),                          \
_skim_muon_track_min_pt(-INFINITY),                     \
_skim_jet_min_pt(std::vector<double>(3, -INFINITY)),    \
_skim_jet_average_pt(-INFINITY),                        \
_skim_multiplicity_tracklet_min_n(INT_MIN),             \
_skim_sum_eg_ntrial(0),                                 \
_stored_track_min_pt(-INFINITY),                        \
_stored_jet_min_pt_raw(-INFINITY),                      \
_nrandom_isolation(0),                                  \
_is_embed(false),                                       \
_emcal_mask(std::vector<bool>()),                       \
_emcal_cell_position(NULL),                             \
_keras_model_photon_discrimination(NULL),               \
_alien_plugin(NULL),                                    \
_metadata_filled(false)

AliAnalysisTaskNTGJ::AliAnalysisTaskNTGJ(const char *name)
    : AliAnalysisTaskEmcal(name), CLASS_INITIALIZATION
{
    DefineOutput(1, TTree::Class());
}

AliAnalysisTaskNTGJ::AliAnalysisTaskNTGJ(
    const AliAnalysisTaskNTGJ & x)
    : AliAnalysisTaskEmcal(), CLASS_INITIALIZATION
{
    _emcal_geometry_filename = x._emcal_geometry_filename;
    _emcal_local2master_filename = x._emcal_local2master_filename;

    _force_ue_subtraction = x._force_ue_subtraction;
    _skim_cluster_min_e = x._skim_cluster_min_e;
    _skim_track_min_pt = x._skim_track_min_pt;
    _skim_muon_track_min_pt = x._skim_muon_track_min_pt;
    _skim_jet_min_pt = x._skim_jet_min_pt;
    _skim_jet_average_pt = x._skim_jet_average_pt;
    _skim_multiplicity_tracklet_min_n =
    x._skim_multiplicity_tracklet_min_n;
    _stored_track_min_pt = x._stored_track_min_pt;
    _stored_jet_min_pt_raw = x._stored_jet_min_pt_raw;
    _nrandom_isolation = x._nrandom_isolation;
    _is_embed = x._is_embed;
}

#undef ntrigger_class

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
#undef BRANCH_STR
#undef BRANCH_STR_ARRAY
#undef B
#undef b
#undef S
#undef s
#undef I
#undef i
#undef F
#undef D
#undef L
#undef l
#undef O

AliAnalysisTaskNTGJ & AliAnalysisTaskNTGJ::operator=(
    const AliAnalysisTaskNTGJ & x)
{
    if (this == &x) {
        return *this;
    }

    _emcal_geometry_filename = x._emcal_geometry_filename;
    _emcal_local2master_filename = x._emcal_local2master_filename;

    _force_ue_subtraction = x._force_ue_subtraction;
    _skim_cluster_min_e = x._skim_cluster_min_e;
    _skim_track_min_pt = x._skim_track_min_pt;
    _skim_muon_track_min_pt = x._skim_muon_track_min_pt;
    _skim_jet_min_pt = x._skim_jet_min_pt;
    _skim_jet_average_pt = x._skim_jet_average_pt;
    _skim_multiplicity_tracklet_min_n =
        x._skim_multiplicity_tracklet_min_n;
    _stored_track_min_pt = x._stored_track_min_pt;
    _stored_jet_min_pt_raw = x._stored_jet_min_pt_raw;
    _nrandom_isolation = x._nrandom_isolation;
    _is_embed = x._is_embed;

    return *this;
}

AliAnalysisTaskNTGJ::~AliAnalysisTaskNTGJ(void)
{
    if (_tree_event != NULL) {
        delete _tree_event;
    }
    if (_emcal_cell_position != NULL) {
        delete reinterpret_cast<std::vector<point_2d_t> *>
        (_emcal_cell_position);
    }
    if (_keras_model_photon_discrimination != NULL) {
        delete reinterpret_cast<KerasModel *>
        (_keras_model_photon_discrimination);
    }

    delete _reco_util;
}

void AliAnalysisTaskNTGJ::UserCreateOutputObjects(void)
{
    AliAnalysisTaskEmcal::UserCreateOutputObjects();
    // fAliEventCuts.fGreenLight = true;

    TFile *file = OpenFile(1);

    if (file != NULL) {
        file->SetCompressionSettings(ROOT::CompressionSettings(
                                         ROOT::kZLIB, 9));
    }

    /////////////////////////////////////////////////////////////////

    _tree_event = new TTree("_tree_event", "");

#define BRANCH(b, t)                    \
    _tree_event->Branch(                \
#b, &_branch_ ## b, #b "/" #t);
#define BRANCH_ARRAY(b, d, t)                   \
    _tree_event->Branch(                        \
#b, _branch_ ## b, #b "[" #d "]/" #t);
#define BRANCH_ARRAY2(b, d, e, t)                       \
    _tree_event->Branch(                                \
#b, _branch_ ## b, #b "[" #d "][" #e "]/" #t);
#define BRANCH_STR(b)                           \
    _tree_event->Branch(                        \
#b, _branch_ ## b, #b "/C");
#define BRANCH_STR_ARRAY(b, d)                  \
    _tree_event->Branch(                        \
#b, &_branch_ ## b);

    MEMBER_BRANCH;

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
#undef BRANCH_STR
#undef BRANCH_STR_ARRAY

    PostData(1, _tree_event);

    /////////////////////////////////////////////////////////////////
}

#undef MEMBER_BRANCH

// lines 1138-1160; not sure where else to put this
#define SAFE_MC_TRUTH_INDEX_TO_USHRT(mc_truth_index)\
    !(mc_truth_index >= 0 &&                        \
      static_cast<size_t>(mc_truth_index) <         \
      stored_mc_truth_index.size()) ? USHRT_MAX :   \
        stored_mc_truth_index[mc_truth_index] ==    \
        ULONG_MAX ?                                 \
        USHRT_MAX :                                 \
        std::min(static_cast<size_t>(USHRT_MAX),    \
                 std::max(static_cast<size_t>(0),   \
                          stored_mc_truth_index     \
                          [mc_truth_index]));

Bool_t AliAnalysisTaskNTGJ::Run()
{
    AliVEvent *event = NULL;
    AliESDEvent *esd_event = NULL;
    AliAODEvent *aod_event = NULL;
    const AliVVertex *primary_vertex = NULL;

    AliVCaloCells *emcal_cell = NULL;
    AliClusterContainer *cluster_container = NULL;
    std::vector<AliTrackContainer*> track_containers;
    AliMCParticleContainer *mc_container = NULL;

    // for some reason, moving loadEmcalGeometry later causes a change in cell_position
    loadEmcalGeometry();

    // event-level initialization and branches
    if (!getEvent(event, esd_event, aod_event)) { // getEvent returns false if the event is null
        return false;
    }

    getContainers(cluster_container, track_containers, mc_container, emcal_cell, event);
    if (mc_container) {
        loadMC(aod_event);
    }

    getMultiplicityCentralityEventPlane(event);
    getBeamProperties(event, esd_event, aod_event, primary_vertex);
    getMetadata(esd_event, aod_event);

    if (!skimMultiplicityTracklet(event) || // skimMultiplicityTracklet returns true if we should keep the event
            !skimClusterE(cluster_container)) { // skimClusterE returns true if we should keep the event
        return false;
    }

    // MC particles
    std::vector<size_t> stored_mc_truth_index;
    std::vector<Int_t> reverse_stored_mc_truth_index;
    std::vector<Int_t> reverse_stored_parton_algorithmic_index;
    if (mc_container) {
        getPrimaryMCParticles(mc_container,
                              stored_mc_truth_index,
                              reverse_stored_mc_truth_index,
                              reverse_stored_parton_algorithmic_index);

        doMCParticleLoop(mc_container,
                         esd_event,
                         reverse_stored_mc_truth_index);
    }

    // tracks
    setTrackCuts(); // currently empty
    doTrackLoop(event, aod_event, track_containers, mc_container, stored_mc_truth_index, primary_vertex);

    // clusters and cells
    std::vector<bool> cell_pass_basic_quality;
    loadPhotonNNModel();
    doClusterLoop(event, emcal_cell, cluster_container, stored_mc_truth_index, cell_pass_basic_quality);
    fillCellBranches(emcal_cell, stored_mc_truth_index);

    // voronoi area, UE, jets, isolation
    getUEJetsIsolation(esd_event,
                       mc_container,
                       cluster_container,
                       track_containers,
                       stored_mc_truth_index,
                       reverse_stored_mc_truth_index,
                       reverse_stored_parton_algorithmic_index,
                       cell_pass_basic_quality);

    if (!skimJets()) { // skimJets returns true if we should keep the event
        return false;
    }
    fillMuonBranches(); // currently empty
    fillEgNtrial();

    _tree_event->Fill();
    return true;
}

void AliAnalysisTaskNTGJ::loadEmcalGeometry()
{
    // lines 341-404
    if (!_emcal_geometry_filename.empty() &&
            !_emcal_local2master_filename.empty() &&
            _emcal_geometry == NULL) {
        if (!gSystem->
                AccessPathName(gSystem->ExpandPathName(
                                   _emcal_geometry_filename.c_str())))
        {
            TGeoManager::Import(_emcal_geometry_filename.c_str());
            _emcal_geometry = AliEMCALGeometry::
                              GetInstance(_emcal_geometry_name);
        }

        AliOADBContainer emcal_geometry_container("emcal");

        if (!gSystem->AccessPathName(gSystem->ExpandPathName(
                                         _emcal_local2master_filename.c_str())))
        {
            emcal_geometry_container.
            InitFromFile(_emcal_local2master_filename.c_str(), "AliEMCALgeo");
        }

        TObjArray *geometry_matrix = dynamic_cast<TObjArray *>(
                                         emcal_geometry_container.GetObject(
                                             _branch_run_number, "EmcalMatrices"));

        if (_emcal_geometry != NULL && geometry_matrix != NULL) {
            const Int_t nsm = _emcal_geometry->GetEMCGeometry()->
                              GetNumberOfSuperModules();

            for (Int_t sm = 0; sm < nsm; sm++) {

                _emcal_geometry->SetMisalMatrix(
                    dynamic_cast<TGeoHMatrix *>(
                        geometry_matrix->At(sm)), sm);
            }
            _branch_has_misalignment_matrix = true;
        }
    }

    if (_emcal_mask.size() != EMCAL_NCELL) {
        _emcal_mask.resize(EMCAL_NCELL);
#if 0 // Keep = 1 to for an actual EMCAL mask (and not all channels
        // turned on)
        for (unsigned int i = 0; i < EMCAL_NCELL; i++) {
            _emcal_mask[i] = inside_edge(i, 1);
        }
        // #include "bad_channel.h"
        for (unsigned int i = 0; bad_channel_emcal[i] != -1; i++) {
            if (inside_edge(bad_channel_emcal[i], 1)) {
                unsigned int bad_cell_3_3[9];

                cell_3_3(bad_cell_3_3, bad_channel_emcal[i]);
                for (size_t j = 0; j < 9; j++) {
                    _emcal_mask[bad_cell_3_3[j]] = false;
                }
            }
        }
#else
        for (unsigned int i = 0; i < EMCAL_NCELL; i++) {
            _emcal_mask[i] = true;
        }
#endif
    }
}

bool AliAnalysisTaskNTGJ::getEvent(AliVEvent *&event,
                                   AliESDEvent *&esd_event,
                                   AliAODEvent *&aod_event)
{
    // lines 406-432, 436-449
    event = InputEvent();

    if (event == NULL) {
        return false;
    }

    if (event->GetRunNumber() != _run_number_current) {
        _metadata_filled = false;
        _run_number_current = event->GetRunNumber();
        if (_muon_track_cut != NULL) {
            _muon_track_cut->SetAllowDefaultParams(kTRUE);
            _muon_track_cut->SetRun(
                (AliInputEventHandler *)
                ((AliAnalysisManager::GetAnalysisManager())->
                 GetInputEventHandler()));
        }
    }

    esd_event = dynamic_cast<AliESDEvent *>(event);

    if (esd_event == NULL) {
        aod_event = dynamic_cast<AliAODEvent *>(event);
        if (aod_event == NULL) {
            return false;
        }
    }

    _branch_run_number = event->GetRunNumber();
    _branch_period_number = event->GetPeriodNumber();
    _branch_orbit_number = event->GetOrbitNumber();
    _branch_bunch_crossing_number = event->GetBunchCrossNumber();
    _branch_time_stamp = esd_event != NULL ?
                         esd_event->GetTimeStamp() :
                         aod_event->GetTimeStamp();
    _branch_trigger_mask[0] = esd_event != NULL ?
                              esd_event->GetTriggerMask() :
                              aod_event->GetTriggerMask();
    _branch_trigger_mask[1] = esd_event != NULL ?
                              esd_event->GetTriggerMaskNext50() :
                              aod_event->GetTriggerMaskNext50();
    _branch_has_misalignment_matrix = false;

    return true;
}

void AliAnalysisTaskNTGJ::getContainers(AliClusterContainer *&cluster_container,
                                        std::vector<AliTrackContainer*> &track_containers,
                                        AliMCParticleContainer *&mc_container,
                                        AliVCaloCells *&emcal_cell,
                                        AliVEvent *event)
{
    // set cluster_container, track_containers, and mc_container
    // by taking in references to pointers
    cluster_container = GetClusterContainer(0);
    track_containers.push_back(GetTrackContainer(0));
    if (_is_embed) {
        track_containers.push_back(GetTrackContainer(1));
    }
    mc_container = GetMCParticleContainer("mcparticles");

    AliDebugStream(2) << "Event has " << cluster_container->GetNAcceptEntries() << " clusters" << std::endl;
    int ntrack = 0;
    for (auto track_container : track_containers) {
        ntrack += track_container->GetNAcceptEntries();
    }
    AliDebugStream(2) << "Event has " << ntrack << " tracks" << std::endl;
    if (mc_container) {
        AliDebugStream(2) << "Event has " << mc_container->GetNParticles() << " MC particles" << std::endl;
    }

    if (_is_embed) {
        emcal_cell = (AliVCaloCells*) event->FindListObject("emcalCellsCombined");
    } else {
        emcal_cell = event->GetEMCALCells();
    }
}

void AliAnalysisTaskNTGJ::setTrackCuts()
{

}

void AliAnalysisTaskNTGJ::getMultiplicityCentralityEventPlane(AliVEvent *event)
{
    // lines 547-604
    AliVVZERO *v0 = event->GetVZEROData();

    if (v0 != NULL) {
        for (size_t i = 0; i < 64; i++) {
            _branch_multiplicity_v0[i] = v0->GetMultiplicity(i);
        }
    }

    static const char *centrality_method[9] = {
        "V0M", "CL0", "CL1",
        "V0Mplus05", "V0Mplus10", "V0Mminus05", "V0Mminus10",
        "SPDClustersCorr", "SPDTracklets"
    };

    AliMultSelection *mult_selection = static_cast<AliMultSelection *>
                                       (event->FindListObject("MultSelection"));

    std::fill(_branch_centrality,
              _branch_centrality + sizeof(_branch_centrality) /
              sizeof(*_branch_centrality), NAN);
    if (mult_selection != NULL) {
        for (size_t i = 0; i < 9; i++) {
            _branch_centrality[i] =
                mult_selection->GetMultiplicityPercentile(
                    centrality_method[i]);
        }
    }
    else {
        AliCentrality *centrality = event->GetCentrality();

        if (centrality != NULL) {
            for (size_t i = 0; i < 9; i++) {
                _branch_centrality[i] =
                    centrality->GetCentralityPercentile(
                        centrality_method[i]);
            }
        }
    }
    // Copy for easy cutting (where duplicated global variable are not
    // very wasteful)
    _branch_centrality_v0m = _branch_centrality[0];

    std::fill(_branch_event_plane_psi_v0,
              _branch_event_plane_psi_v0 +
              sizeof(_branch_event_plane_psi_v0) /
              sizeof(*_branch_event_plane_psi_v0), NAN);

    AliEventplane *event_plane = event->GetEventplane();

    if (event_plane != NULL) {
        for (size_t i = 1; i < 4; i++) {
            _branch_event_plane_psi_v0[i - 1] =
                event->GetEventplane()->CalculateVZEROEventPlane(
                    event, 10, i,
                    _branch_event_plane_q_v0[i - 1][0],
                    _branch_event_plane_q_v0[i - 1][1]);
        }
    }
}

void AliAnalysisTaskNTGJ::loadPhotonNNModel()
{
    // lines 606-611
    if (_keras_model_photon_discrimination == NULL) {
        _keras_model_photon_discrimination = new KerasModel;
        reinterpret_cast<KerasModel *>(
            _keras_model_photon_discrimination)->
        LoadModel("photon_discr.model");
    }
}

void AliAnalysisTaskNTGJ::loadMC(AliAODEvent *aod_event)
{
    // lines 616-683
    _branch_eg_signal_process_id = INT_MIN;
    _branch_eg_mpi = INT_MIN;
    _branch_eg_pt_hat = NAN;
    _branch_eg_cross_section = NAN;
    _branch_eg_weight = NAN;
    // This should be default to zero to avoid counting mishap
    _branch_eg_ntrial = 0;

    _branch_eg_scale_pdf = NAN;
    _branch_eg_alpha_qcd = NAN;
    _branch_eg_alpha_qed = NAN;
    std::fill(_branch_eg_pdf_id,
              _branch_eg_pdf_id + sizeof(_branch_eg_pdf_id) /
              sizeof(*_branch_eg_pdf_id), INT_MIN);
    std::fill(_branch_eg_pdf_x,
              _branch_eg_pdf_x + sizeof(_branch_eg_pdf_x) /
              sizeof(*_branch_eg_pdf_x), NAN);
    std::fill(_branch_eg_pdf_x_pdf,
              _branch_eg_pdf_x_pdf + sizeof(_branch_eg_pdf_x_pdf) /
              sizeof(*_branch_eg_pdf_x_pdf), NAN);

    AliGenPythiaEventHeader *mc_truth_pythia_header = NULL;

    // embedding MC header is done separately for now
    if (_is_embed) {
        mc_truth_pythia_header = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetPythiaHeader();
    }
    else {
        AliAODMCHeader *aod_mc_header = dynamic_cast<AliAODMCHeader*>
                                        (aod_event->GetList()->FindObject(AliAODMCHeader::StdBranchName()));

        Int_t nGenerators = aod_mc_header->GetNCocktailHeaders();
        for (Int_t igen = 0; igen < nGenerators; igen++) {
            AliGenEventHeader *eventHeaderGen = aod_mc_header->GetCocktailHeader(igen);
            TString name = eventHeaderGen->GetName();

            if (name.CompareTo("Pythia") == 0) {
                mc_truth_pythia_header = dynamic_cast<AliGenPythiaEventHeader*>(eventHeaderGen);
            }
        }
    }

    if (mc_truth_pythia_header != NULL) {
        _branch_eg_weight = mc_truth_pythia_header->EventWeight();

        TArrayF eg_primary_vertex(3);

        mc_truth_pythia_header->PrimaryVertex(eg_primary_vertex);
        for (Int_t i = 0; i < 3; i++) {
            _branch_eg_primary_vertex[i] = eg_primary_vertex.At(i);
        }

        _branch_eg_signal_process_id = mc_truth_pythia_header->ProcessType();
        _branch_eg_mpi = mc_truth_pythia_header->GetNMPI();
        _branch_eg_pt_hat = mc_truth_pythia_header->GetPtHard();

        _branch_eg_cross_section =  mc_truth_pythia_header->GetXsection();

        // Count ntrial, because the event might get skimmed away
        // by the ntuplizer
        _skim_sum_eg_ntrial += mc_truth_pythia_header->Trials();
    }
}

void AliAnalysisTaskNTGJ::getBeamProperties(AliVEvent *event,
        AliESDEvent *esd_event,
        AliAODEvent *aod_event,
        const AliVVertex *&primary_vertex)
{
    // lines 685-792
    if (esd_event != NULL) {
        esd_event->InitMagneticField();
    }
    else if (aod_event != NULL) {
        aod_event->InitMagneticField();
    }

    getPrimaryVertex(event, esd_event, primary_vertex);
    getPrimaryVertexSPD(event, esd_event, aod_event);
    getPileup(esd_event, aod_event);

    _branch_event_selected = fInputHandler->IsEventSelected();
}

void AliAnalysisTaskNTGJ::getPrimaryVertex(AliVEvent *event,
        AliESDEvent *esd_event,
        const AliVVertex *&primary_vertex)
{
    // lines 692-717
    std::fill(_branch_primary_vertex,
              _branch_primary_vertex +
              sizeof(_branch_primary_vertex) /
              sizeof(*_branch_primary_vertex), NAN);
    std::fill(_branch_primary_vertex_sigma,
              _branch_primary_vertex_sigma +
              sizeof(_branch_primary_vertex_sigma) /
              sizeof(*_branch_primary_vertex_sigma), NAN);
    _branch_primary_vertex_ncontributor = INT_MIN;

    primary_vertex = event->GetPrimaryVertex();

    if (primary_vertex != NULL) {
        primary_vertex->GetXYZ(_branch_primary_vertex);
    }
    if (esd_event != NULL) {
        const AliESDVertex *esd_primary_vertex =
            esd_event->GetPrimaryVertex();

        if (esd_primary_vertex != NULL) {
            esd_primary_vertex->GetSigmaXYZ(
                _branch_primary_vertex_sigma);
            _branch_primary_vertex_ncontributor =
                esd_primary_vertex->GetNContributors();
        }
    }
}

void AliAnalysisTaskNTGJ::getPrimaryVertexSPD(AliVEvent *event,
        AliESDEvent *esd_event,
        AliAODEvent *aod_event)
{
    // lines 719-745
    std::fill(_branch_primary_vertex_spd,
              _branch_primary_vertex_spd +
              sizeof(_branch_primary_vertex_spd) /
              sizeof(*_branch_primary_vertex_spd), NAN);
    std::fill(_branch_primary_vertex_spd_sigma,
              _branch_primary_vertex_spd_sigma +
              sizeof(_branch_primary_vertex_spd_sigma) /
              sizeof(*_branch_primary_vertex_spd_sigma), NAN);
    _branch_primary_vertex_spd_ncontributor = INT_MIN;

    const AliVVertex *primary_vertex_spd =
        event->GetPrimaryVertexSPD();

    if (primary_vertex_spd != NULL) {
        primary_vertex_spd->GetXYZ(_branch_primary_vertex_spd);
    }
    if (esd_event != NULL) {
        const AliESDVertex *esd_primary_vertex_spd =
            esd_event->GetPrimaryVertexSPD();

        if (esd_primary_vertex_spd != NULL) {
            esd_primary_vertex_spd->GetSigmaXYZ(
                _branch_primary_vertex_spd_sigma);
            _branch_primary_vertex_spd_ncontributor =
                esd_primary_vertex_spd->GetNContributors();
        }
    }
    else if (aod_event != NULL) {
        // lines 762-770
        const AliAODVertex *aod_primary_vertex_spd =
            aod_event->GetPrimaryVertexSPD();

        /*if (aod_primary_vertex_spd != NULL) {
          aod_primary_vertex_spd->GetSigmaXYZ(
          _branch_primary_vertex_spd_sigma);
          _branch_primary_vertex_spd_ncontributor =
          aod_primary_vertex_spd->GetNContributors();
          }//*/
    }
}

void AliAnalysisTaskNTGJ::getPileup(AliESDEvent *esd_event,
                                    AliAODEvent *aod_event)
{
    // lines 747-761, 772-780
    _branch_is_pileup_from_spd_3_08 = false;
    _branch_is_pileup_from_spd_5_08 = false;
    _branch_npileup_vertex_spd = INT_MIN;
    _branch_ncluster_tpc = INT_MIN;
    if (esd_event != NULL) {
        _branch_is_pileup_from_spd_3_08 =
            esd_event->IsPileupFromSPD(3, 0.8);
        _branch_is_pileup_from_spd_5_08 =
            esd_event->IsPileupFromSPD(5, 0.8);
        _branch_npileup_vertex_spd =
            esd_event->GetNumberOfPileupVerticesSPD();
        _branch_ncluster_tpc =
            esd_event->GetNumberOfTPCClusters();
    }
    else if (aod_event != NULL) {
        _branch_is_pileup_from_spd_3_08 =
            aod_event->IsPileupFromSPD(3, 0.8);
        _branch_is_pileup_from_spd_5_08 =
            aod_event->IsPileupFromSPD(5, 0.8);
        _branch_npileup_vertex_spd =
            aod_event->GetNumberOfPileupVerticesSPD();
        _branch_ncluster_tpc =
            aod_event->GetNumberOfTPCClusters();
    }
}

bool AliAnalysisTaskNTGJ::skimMultiplicityTracklet(AliVEvent *event)
{
    // lines 784-792
    // returns true if we should KEEP the event
    if (_skim_multiplicity_tracklet_min_n > INT_MIN) {
        AliVMultiplicity *multiplicity = event->GetMultiplicity();

        if (multiplicity == NULL ||
                !(multiplicity->GetNumberOfTracklets() >=
                  _skim_multiplicity_tracklet_min_n)) {
            return false;
        }
    }

    return true;
}

bool AliAnalysisTaskNTGJ::skimClusterE(AliClusterContainer *cluster_container)
{
    // lines 811-829
    if (_skim_cluster_min_e > -INFINITY) {
        double cluster_e_max = -INFINITY;

        for (auto c : cluster_container->accepted()) {
            if (!(c->GetNCells() > 1) ||
                    cell_masked(c, _emcal_mask)) {
                continue;
            }

            TLorentzVector p;

            c->GetMomentum(p, _branch_primary_vertex);
            cluster_e_max = std::max(cluster_e_max, p.E());
        }
        if (!(cluster_e_max >= _skim_cluster_min_e)) {
            // Discard this event
            return false;
        }
    }

    return true; //if no min_e set, return true to keep event
}

void AliAnalysisTaskNTGJ::getMetadata(AliESDEvent *esd_event,
                                      AliAODEvent *aod_event)
{
    // line 434, moved here because it doesn't seem to be used anywhere else
    alice_jec_t jec;

    // lines 831-889
    if (!_metadata_filled) {
        // Use gitattributes ident mechanism to track the blob object
        // name
        strncpy(_branch_id_git, "$Id$", BUFSIZ);
        _branch_id_git[BUFSIZ - 1] = '\0';
        strncpy(_branch_version_jec, jec.version(), BUFSIZ);
        _branch_version_jec[BUFSIZ - 1] = '\0';
        if (esd_event != NULL) {
            _branch_beam_energy = esd_event->GetBeamEnergy();
            for (size_t i = 0; i < 2; i++) {
                _branch_beam_particle[i] =
                    esd_event->GetBeamParticle(i);
            }

            const AliESDRun *esd_run = esd_event->GetESDRun();

            if (esd_run != NULL) {
                _branch_trigger_class.clear();
                for (size_t i = 0; i < NTRIGGER_CLASS_MAX; i++) {
                    _branch_trigger_class.push_back(
                        std::string(esd_run->GetTriggerClass(i)));
                }
            }
        }
        else if (aod_event != NULL) {
            // FIXME: AOD not handled
            _branch_beam_energy = NAN;
            std::fill(_branch_beam_particle,
                      _branch_beam_particle +
                      sizeof(_branch_beam_particle) /
                      sizeof(*_branch_beam_particle), 0);
        }

#ifdef WITH_EFP7
        strncpy(_branch_debug_blas_version,
                cblas.version_str().c_str(), BUFSIZ);
#else // WITH_EFP7
        _branch_debug_blas_version[0] = '\0';
#endif // WITH_EFP7

        _metadata_filled = true;
    }
    else {
        // To make sure no space is wasted, set metadata in all
        // subsequent events to empty
        _branch_id_git[0] = '\0';
        _branch_version_aliroot[0] = '\0';
        _branch_version_aliphysics[0] = '\0';
        _branch_version_jec[0] = '\0';
        _branch_grid_data_dir[0] = '\0';
        _branch_grid_data_pattern[0] = '\0';
        std::fill(_branch_beam_particle,
                  _branch_beam_particle +
                  sizeof(_branch_beam_particle) /
                  sizeof(*_branch_beam_particle), 0);
        _branch_trigger_class.clear();

        _branch_debug_blas_version[0] = '\0';
    }
}

// =================================================================================================================
void AliAnalysisTaskNTGJ::getPrimaryMCParticles(AliMCParticleContainer *mc_container,
        std::vector<size_t> &stored_mc_truth_index,
        std::vector<Int_t> &reverse_stored_mc_truth_index,
        std::vector<Int_t> &reverse_stored_parton_algorithmic_index)
{
    // lines 948-998
    AliDebugStream(3) << "loops 1 and 2 through MC container" << std::endl;

    stored_mc_truth_index.resize(mc_container->GetNParticles(), ULONG_MAX);
    size_t nmc_truth = 0;
    for (Int_t i = 0; i < mc_container->GetNParticles(); i++) {
        // Bookkeeping for primary final state particles
        if (final_state_primary(mc_container, i)) {
            stored_mc_truth_index.at(i) = nmc_truth;
            reverse_stored_mc_truth_index.push_back(i);
            nmc_truth++;
        }


        // Bookkeeping for partons
        if (parton_cms_algorithmic(mc_container, i)) {
            reverse_stored_parton_algorithmic_index.push_back(i);
        }
    }

    // Assign secondaries to the primaries
    for (Int_t i = 0;
            i < mc_container->GetNParticles(); i++) {
        // Skip primaries
        if (final_state_primary(mc_container, i)) {
            continue;
        }

        Int_t j = i;
        bool has_physical_primary_ancestor = false;

        // Ensure termination even if there is a loop
        for (size_t k = 0; k < 1000; k++) {
            const AliAODMCParticle *p = mc_container->GetMCParticle(j);
            if (p == NULL) {
                break;
            }
            j = p->GetMother();
            if (!(j >= 0 &&
                    j < mc_container->GetNParticles())) {
                break;
            }
            if (final_state_primary(mc_container, j)) {
                has_physical_primary_ancestor = true;
                break;
            }
        }
        if (has_physical_primary_ancestor) {
            stored_mc_truth_index.at(i) = stored_mc_truth_index.at(j);
        }
    }
}
// =================================================================================================================
void AliAnalysisTaskNTGJ::doTrackLoop(AliVEvent *event,
                                      AliAODEvent *aod_event,
                                      std::vector<AliTrackContainer*> track_containers,
                                      AliMCParticleContainer *mc_container,
                                      std::vector<size_t> stored_mc_truth_index,
                                      const AliVVertex *primary_vertex)
{
    // lines 1004-1009, 1016-1136, 1150-1160
    std::fill(_branch_met_tpc,
              _branch_met_tpc + sizeof(_branch_met_tpc) /
              sizeof(*_branch_met_tpc), 0);
    std::fill(_branch_met_its,
              _branch_met_its + sizeof(_branch_met_its) /
              sizeof(*_branch_met_its), 0);

    double met_tpc_kahan_error[2] = { 0, 0 };
    double met_its_kahan_error[2] = { 0, 0 };

    _branch_ntrack = 0;
    if (aod_event != NULL) {
        Int_t itrack = 0;
        for (auto track_container : track_containers) {
            for (auto track : track_container->accepted()) {
                if (track == NULL)
                    continue;

                AliAODTrack * t = static_cast<AliAODTrack*>(track);

                /*bits 1 and 2 are Standard hybrid track cuts for pp and pPb
                  bits 3 and 4 are PbPb cuts
                  bit 5 is ITS only cuts
                */
                UInt_t _local_track_cut_bits = get_local_track_cut_bits(t, event);
                if (_local_track_cut_bits == 0)
                    continue;

                if ((_local_track_cut_bits & 15) != 0) {
                    kahan_sum(_branch_met_tpc[0], met_tpc_kahan_error[0], t->Px());
                    kahan_sum(_branch_met_tpc[1], met_tpc_kahan_error[1], t->Py());
                }

                if ((_local_track_cut_bits & 16) != 0) {
                    kahan_sum(_branch_met_its[0], met_its_kahan_error[0], t->Px());
                    kahan_sum(_branch_met_its[1], met_its_kahan_error[1], t->Py());
                }

                _branch_track_e[_branch_ntrack] = half(t->E());
                _branch_track_pt[_branch_ntrack] = half(t->Pt());
                _branch_track_eta[_branch_ntrack] = half(t->Eta());
                _branch_track_phi[_branch_ntrack] =
                    half(angular_range_reduce(t->Phi()));
                if (gGeoManager != NULL) {
                    _branch_track_eta_emcal[_branch_ntrack] =
                        half(t->GetTrackEtaOnEMCal());
                    _branch_track_phi_emcal[_branch_ntrack] =
                        half(angular_range_reduce(t->GetTrackPhiOnEMCal()));
                }
                else {
                    _branch_track_eta_emcal[_branch_ntrack] = NAN;
                    _branch_track_phi_emcal[_branch_ntrack] = NAN;
                }
                _branch_track_charge[_branch_ntrack] =
                    std::min(static_cast<Short_t>(CHAR_MAX),
                             std::max(static_cast<Short_t>(CHAR_MIN), t->Charge()));

                _branch_track_quality[_branch_ntrack] = _local_track_cut_bits;

                _branch_track_tpc_dedx[_branch_ntrack] = half(t->GetTPCsignal());

                static const Int_t mode_inner_wall = 1;
                static const Double_t dead_zone_width_cm = 2;
                static const Double_t max_z_cm = 220;

                _branch_track_tpc_length_active_zone[_branch_ntrack] = NAN;

                /*if (t->GetInnerParam() != NULL) {
                  _branch_track_tpc_length_active_zone
                  [_branch_ntrack] =
                  half(t->GetLengthInActiveZone
                  (mode_inner_wall, dead_zone_width_cm,
                  max_z_cm, event->GetMagneticField()));
                  }//*/

                _branch_track_tpc_xrow[_branch_ntrack] =
                    std::min(static_cast<Float_t>(UCHAR_MAX),
                             std::max(0.0F, t->GetTPCCrossedRows()));
                _branch_track_tpc_ncluster[_branch_ntrack] = std::min(static_cast<UShort_t>(UCHAR_MAX), std::max(static_cast<UShort_t>(0), t->GetTPCNcls()));
                // _branch_track_tpc_ncluster[_branch_ntrack] =
                //     std::min(UCHAR_MAX,
                //              std::max(0, t->GetTPCNcls()));
                _branch_track_tpc_ncluster_dedx[_branch_ntrack] =
                    std::min(static_cast<UShort_t>(UCHAR_MAX),
                             std::max(static_cast<UShort_t>(0),
                                      t->GetTPCsignalN()));
                _branch_track_tpc_ncluster_findable[_branch_ntrack] =
                    std::min(static_cast<UShort_t>(UCHAR_MAX),
                             std::max(static_cast<UShort_t>(0),
                                      t->GetTPCNclsF()));
                _branch_track_its_ncluster[_branch_ntrack] =
                    t->GetITSNcls();
                _branch_track_its_chi_square[_branch_ntrack] =
                    half(t->GetITSchi2());

                Double_t dz[2] = { NAN, NAN };
                Double_t cov[3] = { NAN, NAN, NAN };

                if (t->PropagateToDCA
                        (primary_vertex, event->GetMagneticField(),
                         kVeryBig, dz, cov) == kTRUE) {
                    _branch_track_dca_xy[_branch_ntrack] =
                        half(dz[0]);
                    _branch_track_dca_z[_branch_ntrack] =
                        half(dz[1]);
                }

                const Int_t mc_truth_index = t->GetLabel();

                _branch_track_mc_truth_index[_branch_ntrack] =
                    SAFE_MC_TRUTH_INDEX_TO_USHRT(mc_truth_index);
                _branch_track_voronoi_area[_branch_ntrack] = 0;

                _branch_ntrack++;
                if (_branch_ntrack >= NTRACK_MAX) {
                    break;
                }

                itrack++;
            }
        }
    }
}

// =================================================================================================================
void AliAnalysisTaskNTGJ::doMCParticleLoop(AliMCParticleContainer *mc_container,
        AliESDEvent *esd_event,
        std::vector<Int_t> reverse_stored_mc_truth_index)
{
    // lines 1336-1363, 1412-1474
    AliDebugStream(3) << "loop 3 through MC container" << std::endl;

    // may need to move this outside of this function
    _branch_nmc_truth = 0;

    for (std::vector<Int_t>::const_iterator iterator =
                reverse_stored_mc_truth_index.begin();
            iterator != reverse_stored_mc_truth_index.end();
            iterator++) {
        const AliAODMCParticle *p = mc_container->GetMCParticle(*iterator);

        if (p == NULL) {
            // Keep consistent indexing, though this should never
            // happen
            _branch_mc_truth_e[_branch_nmc_truth] = NAN;
            _branch_mc_truth_pt[_branch_nmc_truth] = NAN;
            _branch_mc_truth_eta[_branch_nmc_truth] = NAN;
            _branch_mc_truth_phi[_branch_nmc_truth] = NAN;
            _branch_mc_truth_charge[_branch_nmc_truth] =
                CHAR_MIN;
            _branch_mc_truth_pdg_code[_branch_nmc_truth] =
                SHRT_MIN;
            _branch_mc_truth_status[_branch_nmc_truth] =
                UCHAR_MAX;
            _branch_mc_truth_generator_index[_branch_nmc_truth] =
                UCHAR_MAX;
            _branch_nmc_truth++;
            continue;
        }

        _branch_mc_truth_e[_branch_nmc_truth] = half(p->E());
        _branch_mc_truth_pt[_branch_nmc_truth] = half(p->Pt());
        _branch_mc_truth_eta[_branch_nmc_truth] = half(p->Eta());
        _branch_mc_truth_phi[_branch_nmc_truth] =
            half(angular_range_reduce(p->Phi()));
        _branch_mc_truth_charge[_branch_nmc_truth] =
            std::min(static_cast<Short_t>(CHAR_MAX),
                     std::max(static_cast<Short_t>(CHAR_MIN),
                              p->Charge()));
        _branch_mc_truth_pdg_code[_branch_nmc_truth] =
            std::min(SHRT_MAX,
                     std::max(SHRT_MIN, p->PdgCode()));
        _branch_mc_truth_status[_branch_nmc_truth] =
            std::min(static_cast<UInt_t>(UCHAR_MAX),
                     std::max(static_cast<UInt_t>(0),
                              mc_container->
                              GetMCParticle(*iterator)->
                              MCStatusCode()));
        // _branch_mc_truth_final_state_primary[_branch_nmc_truth] =
        //     final_state_primary(mc_container, *iterator);
        // _branch_mc_truth_first_parent[_branch_nmc_truth] =
        //     p->GetMother();
        // _branch_mc_truth_first_child[_branch_nmc_truth] =
        //     p->GetDaughterFirst();
        // _branch_mc_truth_first_child[_branch_nmc_truth] =
        //     p->GetDaughterLast();
        _branch_mc_truth_generator_index[_branch_nmc_truth] =
            std::min(static_cast<Short_t>(UCHAR_MAX),
                     std::max(static_cast<Short_t>(0),
                              p->GetGeneratorIndex()));

        if (p->GetMother() >= 0 &&
                p->GetMother() <
                mc_container->GetNParticles()) {
            const AliAODMCParticle *pp =
                mc_container->GetMCParticle(
                    p->GetMother());

            _branch_mc_truth_first_parent_pdg_code
            [_branch_nmc_truth] = pp->PdgCode();
            _branch_mc_truth_first_parent_e
            [_branch_nmc_truth] = half(pp->E());
            _branch_mc_truth_first_parent_pt
            [_branch_nmc_truth] = half(pp->Pt());
            _branch_mc_truth_first_parent_eta
            [_branch_nmc_truth] = half(pp->Eta());
            _branch_mc_truth_first_parent_phi
            [_branch_nmc_truth] = half(pp->Phi());
            _branch_mc_truth_sibling_index[_branch_nmc_truth] =
                *iterator ==
                pp->GetDaughterFirst() ?
                pp->GetDaughterLast() :
                *iterator ==
                pp->GetDaughterLast() ?
                pp->GetDaughterFirst() :
                USHRT_MAX;
        }

        _branch_nmc_truth++;
        if (_branch_nmc_truth >= NMC_TRUTH_MAX) {
            break;
        }
    }

}
// =================================================================================================================
void AliAnalysisTaskNTGJ::doClusterLoop(AliVEvent *event,
                                        AliVCaloCells *emcal_cell,
                                        AliClusterContainer *cluster_container,
                                        std::vector<size_t> stored_mc_truth_index,
                                        std::vector<bool> &cell_pass_basic_quality)
{
    AliVVZERO *v0 = event->GetVZEROData();

    // lines 1666-2319, broken up
    std::fill(_branch_cell_cluster_index,
              _branch_cell_cluster_index +
              sizeof(_branch_cell_cluster_index) /
              sizeof(*_branch_cell_cluster_index), USHRT_MAX);
    _branch_ncluster = 0;

    // FIXME: Turn this into a switch
    static const bool fill_cluster = true;

    const Int_t ncalo_cluster =
        fill_cluster && _nrandom_isolation > 0 ?
        _nrandom_isolation : cluster_container->GetNAcceptEntries();
    AliESDCaloCluster dummy_cluster;

    cell_pass_basic_quality.assign(cluster_container->GetNAcceptEntries(), false);

    AliDebugStream(3) << "loop 3 through clusters; inner loops through cells, tracks, clusters, MC container" << std::endl;
    Int_t i = 0;
    for (auto cluster : cluster_container->accepted()) {
        AliVCluster *c = _nrandom_isolation > 0 ? &dummy_cluster : cluster;

        // cell basic quality -- moved from lines 1193-1202 since this is now the first cluster loop
        Int_t cell_id_max = -1;
        Double_t cell_energy_max = -INFINITY;
        Double_t cell_cross = NAN;

        cell_max_cross(cell_id_max, cell_energy_max, cell_cross,
                       c, emcal_cell);
        if (c->GetNCells() > 1 &&
                cell_cross / cell_energy_max > 0.05 &&
                !cell_masked(c, _emcal_mask)) {
            cell_pass_basic_quality[i] = true;
        }

        if (i >= ncalo_cluster) {
            break;
        }

        // Omitting #if 0 clause (lines 1685-1697)

        // ----------------
        // Cluster Branches
        // ----------------
        // lines 1701-1830
        // set arrays with USHRT_Max. Fill with data in loop
        fillClusterBranches(emcal_cell, c, i, stored_mc_truth_index);

        //-------------------
        // Neural Net Branches
        //-------------------
        // lines 2297-2313
        fillPhotonNNBranches(c, stored_mc_truth_index, emcal_cell, v0);

        // lines 2315-2319
        _branch_ncluster++;
        if (_branch_ncluster >= NCLUSTER_MAX) {
            break;
        }
        i++;
    }
}
void AliAnalysisTaskNTGJ::fillClusterBranches(AliVCaloCells *emcal_cell,
        AliVCluster *c,
        Int_t i,
        std::vector<size_t> stored_mc_truth_index)
{
    // lines 1701-1830
    //FIXME: change away from single letter variables...
    TLorentzVector p;
    c->GetMomentum(p, _branch_primary_vertex);

    _branch_cluster_e[_branch_ncluster] = half(p.E());
    _branch_cluster_pt[_branch_ncluster] = half(p.Pt());
    _branch_cluster_eta[_branch_ncluster] = half(p.Eta());
    _branch_cluster_phi[_branch_ncluster] =
        half(angular_range_reduce(p.Phi()));

    _branch_cluster_lambda_square[_branch_ncluster][0] =
        half(c->GetM02());
    _branch_cluster_lambda_square[_branch_ncluster][1] =
        half(c->GetM20());

    std::fill(_branch_cluster_lambda_square_angle
              [_branch_ncluster],
              _branch_cluster_lambda_square_angle
              [_branch_ncluster] + 2, NAN);
    if (_emcal_geometry != NULL) {
        Float_t l0      = 0;
        Float_t l1      = 0;
        Float_t dispp   = 0;
        Float_t dEta    = 0;
        Float_t dPhi    = 0;
        Float_t sEta    = 0;
        Float_t sPhi    = 0;
        Float_t sEtaPhi = 0;

        _reco_util->SetShowerShapeCellLocationType(1);
        _reco_util->RecalculateClusterShowerShapeParameters(
            _emcal_geometry, emcal_cell, c, l0, l1, dispp,
            dEta, dPhi, sEta, sPhi, sEtaPhi);
        _branch_cluster_lambda_square_angle[_branch_ncluster][0]
            = half(l0);
        _branch_cluster_lambda_square_angle[_branch_ncluster][1]
            = half(l1);
    }

    _branch_cluster_tof[_branch_ncluster] =
        half(c->GetTOF() * 1e+9);
    _branch_cluster_ncell[_branch_ncluster] = c->GetNCells();
    _branch_cluster_nlocal_maxima[_branch_ncluster] =
        c->GetNExMax();
    _branch_cluster_distance_to_bad_channel[_branch_ncluster] =
        half(c->GetDistanceToBadChannel());

    const UShort_t *cell_index = c->GetCellsAbsId();

    if (cell_index != NULL) {
        for (Int_t j = 0; j < c->GetNCells(); j++) {
            if (_branch_cell_cluster_index[cell_index[j]] ==
                    USHRT_MAX) {
                _branch_cell_cluster_index[cell_index[j]] = i;
            }
            else {
                // USHRT_MAX - 1 if there is multiple association
                _branch_cell_cluster_index[cell_index[j]] =
                    USHRT_MAX - 1;
            }
        }
    }

    Int_t cell_id_max = -1;
    Double_t cell_energy_max = -INFINITY;
    Double_t energy_cross = NAN;

    cell_max_cross(cell_id_max, cell_energy_max,
                   energy_cross, c, emcal_cell);

    _branch_cluster_cell_id_max[_branch_ncluster] =
        cell_id_max;

    _branch_cluster_e_max[_branch_ncluster] =
        half(cell_energy_max);

    _branch_cluster_e_cross[_branch_ncluster] =
        half(energy_cross);

    _branch_cluster_nmc_truth[_branch_ncluster] =
        c->GetNLabels();

    std::vector<std::pair<float, unsigned short> > mc_truth_energy_index;
    std::set<Int_t> cluster_mc_truth_index;

    for (UInt_t j = 0; j < c->GetNLabels(); j++) {
        const Int_t mc_truth_index = c->GetLabelAt(j);

        if (mc_truth_index >= 0 &&
                static_cast<size_t>(mc_truth_index) <
                stored_mc_truth_index.size() &&
                stored_mc_truth_index[mc_truth_index] !=
                ULONG_MAX &&
                cluster_mc_truth_index.find(
                    stored_mc_truth_index[mc_truth_index]) ==
                cluster_mc_truth_index.end()) {
            const unsigned short u =
                SAFE_MC_TRUTH_INDEX_TO_USHRT(mc_truth_index);

            mc_truth_energy_index.push_back(
                std::pair<float, unsigned short>(
                    _branch_mc_truth_e[stored_mc_truth_index[
                                           mc_truth_index]],
                    u));
        }
        // Test (above) and insertion must be done via
        // stored_mc_truth_index[], since mc_truth_index are
        // potentially secondaries
        cluster_mc_truth_index.insert(
            stored_mc_truth_index[mc_truth_index]);
    }

    std::sort(mc_truth_energy_index.begin(),
              mc_truth_energy_index.end());
    std::fill(_branch_cluster_mc_truth_index[_branch_ncluster],
              _branch_cluster_mc_truth_index[_branch_ncluster] +
              CLUSTER_NMC_TRUTH_MAX, USHRT_MAX);

    size_t cluster_nmctruth_index = 0;

    for (std::vector<std::pair<float, unsigned short> >::
            const_reverse_iterator iterator =
                mc_truth_energy_index.rbegin();
            iterator != mc_truth_energy_index.rend(); iterator++) {
        if (cluster_nmctruth_index >= CLUSTER_NMC_TRUTH_MAX) {
            break;
        }
        _branch_cluster_mc_truth_index[_branch_ncluster]
        [cluster_nmctruth_index] = iterator->second;
        cluster_nmctruth_index++;
    }
}

void AliAnalysisTaskNTGJ::fillPhotonNNBranches(AliVCluster *c,
        std::vector<size_t> stored_mc_truth_index,
        AliVCaloCells *emcal_cell,
        AliVVZERO *v0)
{
    // lines 2297-2313
    std::fill(_branch_cluster_s_nphoton[_branch_ncluster],
              _branch_cluster_s_nphoton[_branch_ncluster] + 4,
              NAN);
    std::fill(_branch_cluster_s_ncharged_hadron[_branch_ncluster],
              _branch_cluster_s_ncharged_hadron
              [_branch_ncluster] + 4, NAN);

    std::vector<float> output_tensor =
        cluster_cell_keras_inference(
            c, emcal_cell, _branch_primary_vertex, v0,
            *reinterpret_cast<KerasModel *>(
                _keras_model_photon_discrimination));

    if (output_tensor.size() == 2) {
        std::copy(output_tensor.begin(), output_tensor.end(),
                  _branch_cluster_s_nphoton[_branch_ncluster] + 1);
    }
}


void AliAnalysisTaskNTGJ::getUEJetsIsolation(AliESDEvent *esd_event,
        AliMCParticleContainer *mc_container,
        AliClusterContainer *cluster_container,
        std::vector<AliTrackContainer*> track_containers,
        std::vector<size_t> stored_mc_truth_index,
        std::vector<Int_t> reverse_stored_mc_truth_index,
        std::vector<Int_t> reverse_stored_parton_algorithmic_index,
        std::vector<bool> cell_pass_basic_quality)
{
    // A value of 2^(-30) < 10^(-9) would map a 10 TeV particle to
    // less than 10 MeV, sufficient to remove any significant momentum
    // bias while not being too greedy to limit the exponent range.
    // Ghost scaling factor is chosen as power of two to maintain the
    // multiplication/division being numerically exact (provided px,
    // py, pz >= 2^(-1022+30) which is of the order 10^(-290) eV).
    static const double scale_ghost = pow(2.0, -30.0);

    std::vector<fastjet::PseudoJet> jet_truth_ak04;                 // truth anti-kT jets
    std::vector<fastjet::PseudoJet> jet_charged_truth_ak04;         // truth anti-kT charged jets
    fastjet::ClusterSequenceArea *cluster_sequence_truth = NULL;    // truth anti-kT jets

    getTruthJetsAndIsolation(esd_event,
                             mc_container,
                             cluster_container,
                             stored_mc_truth_index,
                             reverse_stored_mc_truth_index,
                             jet_truth_ak04,
                             jet_charged_truth_ak04,
                             cluster_sequence_truth);

    getTpcUEJetsIsolation(track_containers,
                          cluster_container,
                          cell_pass_basic_quality,
                          mc_container,
                          reverse_stored_parton_algorithmic_index,
                          jet_truth_ak04,
                          jet_charged_truth_ak04,
                          cluster_sequence_truth);

    getItsUEJetsIsolation(track_containers,
                          cluster_container,
                          cell_pass_basic_quality,
                          mc_container,
                          reverse_stored_parton_algorithmic_index,
                          jet_truth_ak04,
                          jet_charged_truth_ak04,
                          cluster_sequence_truth);

    getClusterUEIsolation(cluster_container,
                          cell_pass_basic_quality);

    if (cluster_sequence_truth != NULL) {
        delete cluster_sequence_truth;
    }
}

void AliAnalysisTaskNTGJ::getTruthJetsAndIsolation(
    AliESDEvent *esd_event,
    AliMCParticleContainer *mc_container,
    AliClusterContainer *cluster_container,
    std::vector<size_t> stored_mc_truth_index,
    std::vector<Int_t> reverse_stored_mc_truth_index,
    std::vector<fastjet::PseudoJet> &jet_truth_ak04,
    std::vector<fastjet::PseudoJet> &jet_charged_truth_ak04,
    fastjet::ClusterSequenceArea *&cluster_sequence_truth)
{
    // lines 1283-1292
    enum {
        BEAM_PARTICLE_P = 1001
    };

    const bool subtract_ue =
        _force_ue_subtraction ? true :
        esd_event != NULL ?
        !(esd_event->GetBeamParticle(0) == BEAM_PARTICLE_P &&
          esd_event->GetBeamParticle(1) == BEAM_PARTICLE_P) :
        false;

    getTruthJets(subtract_ue,
                 mc_container,
                 reverse_stored_mc_truth_index,
                 jet_truth_ak04,
                 jet_charged_truth_ak04,
                 cluster_sequence_truth);
    getTruthIsolation(cluster_container, mc_container, subtract_ue, stored_mc_truth_index);
}

void AliAnalysisTaskNTGJ::getTruthJets(
    bool subtract_ue,
    AliMCParticleContainer *mc_container,
    std::vector<Int_t> reverse_stored_mc_truth_index,
    std::vector<fastjet::PseudoJet> &jet_truth_ak04,
    std::vector<fastjet::PseudoJet> &jet_charged_truth_ak04,
    fastjet::ClusterSequenceArea *&cluster_sequence_truth)
{
    // lines 1294-1320, 1365-1411, 1476-1490, 1521, 1527-1528, 1540-1545
    static const double jet_antikt_d_04 = 0.4;

    std::vector<fastjet::PseudoJet> particle_truth;
    std::vector<fastjet::PseudoJet> particle_charged_truth;
    fastjet::ClusterSequenceArea *cluster_sequence_charged_truth = NULL;

    // PDG Monte Carlo number scheme

    enum {
        PDG_CODE_ELECTRON_MINUS             = 11,
        PDG_CODE_ELECTRON_NEUTRINO          = 12,
        PDG_CODE_MUON_MINUS                 = 13,
        PDG_CODE_MUON_NEUTRINO              = 14,
        PDG_CODE_TAU_MINUS                  = 15,
        PDG_CODE_TAU_NEUTRINO               = 16,
        PDG_CODE_PHOTON                     = 22,
    };

    std::fill(_branch_met_truth,
              _branch_met_truth + sizeof(_branch_met_truth) /
              sizeof(*_branch_met_truth), 0);

    if (mc_container) {
        double met_truth_kahan_error[2] = { 0, 0 };

        for (std::vector<Int_t>::const_iterator iterator =
                    reverse_stored_mc_truth_index.begin();
                iterator != reverse_stored_mc_truth_index.end();
                iterator++) {
            const AliAODMCParticle *p = mc_container->GetMCParticle(*iterator);

            if (p->GetGeneratorIndex() != 0 || !subtract_ue) {
                const unsigned int abs_pdg_code =
                    std::abs(p->PdgCode());

                switch (abs_pdg_code) {
                case PDG_CODE_ELECTRON_NEUTRINO:
                case PDG_CODE_MUON_NEUTRINO:
                case PDG_CODE_TAU_NEUTRINO:
                    // Remove all (stable) neutrinos from the truth
                    // jet definition
                    break;
                default:
                    particle_truth.push_back(fastjet::PseudoJet(
                                                 p->Px(), p->Py(), p->Pz(), p->P()));
                    switch (abs_pdg_code) {
                    case PDG_CODE_ELECTRON_MINUS:
                    case PDG_CODE_PHOTON:
                        particle_truth.back().
                        set_user_index(USER_INDEX_EM);
                        break;
                    case PDG_CODE_MUON_MINUS:
                        particle_truth.back().
                        set_user_index(USER_INDEX_MUON);
                        break;
                    }
                    // Fill again for charged particle
                    if (p->Charge() != 0) {
                        particle_charged_truth.push_back(fastjet::PseudoJet(
                                                             p->Px(), p->Py(), p->Pz(), p->P()));
                        switch (abs_pdg_code) {
                        case PDG_CODE_ELECTRON_MINUS:
                        case PDG_CODE_PHOTON:
                            particle_charged_truth.back().
                            set_user_index(USER_INDEX_EM);
                            break;
                        case PDG_CODE_MUON_MINUS:
                            particle_charged_truth.back().
                            set_user_index(USER_INDEX_MUON);
                            break;
                        }
                    }
                    kahan_sum(_branch_met_truth[0],
                              met_truth_kahan_error[0], p->Px());
                    kahan_sum(_branch_met_truth[1],
                              met_truth_kahan_error[1], p->Py());
                }
            }
        }

        // run anti-kT clustering algorithm
        cluster_sequence_truth = new fastjet::ClusterSequenceArea(
            particle_truth,
            fastjet::JetDefinition(fastjet::JetDefinition(
                                       fastjet::antikt_algorithm, jet_antikt_d_04)),
            fastjet::VoronoiAreaSpec());
        jet_truth_ak04 = cluster_sequence_truth->inclusive_jets(0);

        cluster_sequence_charged_truth =
            new fastjet::ClusterSequenceArea(
            particle_charged_truth,
            fastjet::JetDefinition(fastjet::JetDefinition(
                                       fastjet::antikt_algorithm, jet_antikt_d_04)),
            fastjet::VoronoiAreaSpec());
        jet_charged_truth_ak04 =
            cluster_sequence_charged_truth->inclusive_jets(0);
    }

    FILL_BRANCH_JET_TRUTH(truth, ak04, jet_truth_ak04);
    FILL_BRANCH_JET_TRUTH(charged_truth, ak04,
                          jet_charged_truth_ak04);

    if (cluster_sequence_charged_truth != NULL) {
        delete cluster_sequence_charged_truth;
    }
}

void AliAnalysisTaskNTGJ::getTruthIsolation(AliClusterContainer *cluster_container,
        AliMCParticleContainer *mc_container,
        bool subtract_ue,
        std::vector<size_t> stored_mc_truth_index)
{
    UInt_t icluster = 0;
    for (auto cluster : cluster_container->accepted()) {
        AliVCluster *c = cluster;
        TLorentzVector p;
        c->GetMomentum(p, _branch_primary_vertex);

        // remake cluster_mc_truth_index here for each cluster
        // rather than pass around a vector of sets (one set per cluster)
        // see fillClusterBranches for more notes
        std::set<Int_t> cluster_mc_truth_index;

        for (UInt_t j = 0; j < c->GetNLabels(); j++) {
            const Int_t mc_truth_index = c->GetLabelAt(j);

            cluster_mc_truth_index.insert(
                stored_mc_truth_index[mc_truth_index]);
        }

        if (mc_container) {
            double cluster_iso_01_truth = 0;
            double cluster_iso_02_truth = 0;
            double cluster_iso_03_truth = 0;
            double cluster_iso_04_truth = 0;

            std::vector<std::pair<double, double> > delta_vs_iso;

            for (Int_t j = 0;
                    j < mc_container->GetNParticles(); j++) {
                if (final_state_primary(mc_container, j) &&
                        cluster_mc_truth_index.find(j) ==
                        cluster_mc_truth_index.end()) {
                    const AliAODMCParticle *t = mc_container->GetMCParticle(j);

                    // charged particles only in jet truth
                    if (t->Charge() != 0) {
                        if (t->GetGeneratorIndex() != 0 || !subtract_ue) {
                            const double dpseudorapidity =
                                t->Eta() - p.Eta();
                            const double dazimuth = angular_range_reduce(
                                                        angular_range_reduce(t->Phi()) -
                                                        angular_range_reduce(p.Phi()));
                            const double dr_2 =
                                std::pow(dpseudorapidity, 2) +
                                std::pow(dazimuth, 2);
                            if (dr_2 < 0.1 * 0.1) {
                                cluster_iso_01_truth += t->Pt();
                            }
                            if (dr_2 < 0.2 * 0.2) {
                                cluster_iso_02_truth += t->Pt();
                            }
                            if (dr_2 < 0.3 * 0.3) {
                                cluster_iso_03_truth += t->Pt();
                            }
                            if (dr_2 < 0.4 * 0.4) {
                                cluster_iso_04_truth += t->Pt();
                                delta_vs_iso.push_back(
                                    std::pair<double, double>(
                                        sqrt(dr_2), t->Pt()));
                            }
                        }
                    }
                }
            }
            _branch_cluster_iso_01_truth[icluster] =
                half(cluster_iso_01_truth);
            _branch_cluster_iso_02_truth[icluster] =
                half(cluster_iso_02_truth);
            _branch_cluster_iso_03_truth[icluster] =
                half(cluster_iso_03_truth);
            _branch_cluster_iso_04_truth[icluster] =
                half(cluster_iso_04_truth);

            _branch_cluster_frixione_04_02_truth[icluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso,
                                              0.4, 0.2));
            _branch_cluster_frixione_04_05_truth[icluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso,
                                              0.4, 0.5));
            _branch_cluster_frixione_04_10_truth[icluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso,
                                              0.4, 1.0));

            _branch_cluster_anti_frixione_04_02_truth
            [icluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso,
                                                   0.4, 0.2));
            _branch_cluster_anti_frixione_04_05_truth
            [icluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso,
                                                   0.4, 0.5));
            _branch_cluster_anti_frixione_04_10_truth
            [icluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso,
                                                   0.4, 1.0));
        } else {
            _branch_cluster_iso_01_truth[icluster] = NAN;
            _branch_cluster_iso_02_truth[icluster] = NAN;
            _branch_cluster_iso_03_truth[icluster] = NAN;
            _branch_cluster_iso_04_truth[icluster] = NAN;

            _branch_cluster_frixione_04_02_truth[icluster] =
                NAN;
            _branch_cluster_frixione_04_05_truth[icluster] =
                NAN;
            _branch_cluster_frixione_04_10_truth[icluster] =
                NAN;

            _branch_cluster_anti_frixione_04_02_truth
            [icluster] = NAN;
            _branch_cluster_anti_frixione_04_05_truth
            [icluster] = NAN;
            _branch_cluster_anti_frixione_04_10_truth
            [icluster] = NAN;
        }

        icluster++;
        if (icluster >= NCLUSTER_MAX) {
            break;
        }
    }
}

void AliAnalysisTaskNTGJ::getTpcUEJetsIsolation(
    std::vector<AliTrackContainer*> track_containers,
    AliClusterContainer *cluster_container,
    std::vector<bool> cell_pass_basic_quality,
    AliMCParticleContainer *mc_container,
    std::vector<Int_t> reverse_stored_parton_algorithmic_index,
    std::vector<fastjet::PseudoJet> jet_truth_ak04,
    std::vector<fastjet::PseudoJet> jet_charged_truth_ak04,
    fastjet::ClusterSequenceArea *cluster_sequence_truth)
{
    std::vector<fastjet::PseudoJet> particle_reco_tpc;  // track pseudojets
    std::map<size_t, size_t> track_reco_index_tpc;      // index for isolation calculation
    std::vector<double> particle_reco_area_tpc;         // track areas
    std::pair<std::pair<std::vector<double>, std::vector<double> >,
        double> ue_estimate_tpc;                        // object containing UE information

    getVoronoiAreaTpc(track_containers,
                      particle_reco_tpc,
                      track_reco_index_tpc,
                      particle_reco_area_tpc);

    // this needs to happen in this order (i.e. copy first, then set_user_index)
    // otherwise some of the jet distributions are wrong
    std::vector<fastjet::PseudoJet> particle_reco_tagged_ak04tpc = particle_reco_tpc;
    for (size_t i = 0; i < particle_reco_tpc.size(); i++) {
        particle_reco_tpc[i].set_user_index(static_cast<int>(i));
    }

    getUETpc(particle_reco_tpc,
             particle_reco_area_tpc,
             ue_estimate_tpc);

    getIsolationTpc(cluster_container,
                    track_containers,
                    track_reco_index_tpc,
                    particle_reco_area_tpc,
                    ue_estimate_tpc);

    getJetsTpc(particle_reco_tpc,
               particle_reco_tagged_ak04tpc,
               jet_truth_ak04,
               jet_charged_truth_ak04,
               cluster_sequence_truth,
               particle_reco_area_tpc,
               ue_estimate_tpc,
               cluster_container,
               cell_pass_basic_quality,
               mc_container,
               reverse_stored_parton_algorithmic_index);
}

void AliAnalysisTaskNTGJ::getItsUEJetsIsolation(
    std::vector<AliTrackContainer*> track_containers,
    AliClusterContainer *cluster_container,
    std::vector<bool> cell_pass_basic_quality,
    AliMCParticleContainer *mc_container,
    std::vector<Int_t> reverse_stored_parton_algorithmic_index,
    std::vector<fastjet::PseudoJet> jet_truth_ak04,
    std::vector<fastjet::PseudoJet> jet_charged_truth_ak04,
    fastjet::ClusterSequenceArea *cluster_sequence_truth)
{
    std::vector<fastjet::PseudoJet> particle_reco_its;  // track pseudojets
    std::map<size_t, size_t> track_reco_index_its;      // index for isolation calculation
    std::vector<double> particle_reco_area_its;         // track areas
    std::pair<std::pair<std::vector<double>, std::vector<double> >,
        double> ue_estimate_its;                        // object containing UE information

    getVoronoiAreaIts(track_containers,
                      particle_reco_its,
                      track_reco_index_its,
                      particle_reco_area_its);

    // this needs to happen in this order (i.e. copy first, then set_user_index)
    // otherwise some of the jet distributions are wrong
    std::vector<fastjet::PseudoJet> particle_reco_tagged_ak04its = particle_reco_its;
    for (size_t i = 0; i < particle_reco_its.size(); i++) {
        particle_reco_its[i].set_user_index(static_cast<int>(i));
    }

    getUEIts(particle_reco_its,
             particle_reco_area_its,
             ue_estimate_its);

    getIsolationIts(cluster_container,
                    track_containers,
                    track_reco_index_its,
                    particle_reco_area_its,
                    ue_estimate_its);

    getJetsIts(particle_reco_its,
               particle_reco_tagged_ak04its,
               jet_truth_ak04,
               jet_charged_truth_ak04,
               cluster_sequence_truth,
               particle_reco_area_its,
               ue_estimate_its,
               cluster_container,
               cell_pass_basic_quality,
               mc_container,
               reverse_stored_parton_algorithmic_index);
}

void AliAnalysisTaskNTGJ::getClusterUEIsolation(
    AliClusterContainer *cluster_container,
    std::vector<bool> cell_pass_basic_quality)
{
    // initialize vectors
    std::vector<fastjet::PseudoJet> particle_reco_cluster;  // cluster pseudojets
    std::map<size_t, size_t> cluster_reco_index;            // index for isolation calculation
    std::vector<double> particle_reco_area_cluster;         // cluster areas
    std::pair<std::pair<std::vector<double>, std::vector<double> >,
        double> ue_estimate_cluster;                        // object containing UE information

    getVoronoiAreaCluster(cluster_container,
                          cell_pass_basic_quality,
                          particle_reco_cluster,
                          cluster_reco_index,
                          particle_reco_area_cluster);

    getUECluster(particle_reco_cluster,
                 particle_reco_area_cluster,
                 ue_estimate_cluster);

    getIsolationCluster(cluster_container,
                        cell_pass_basic_quality,
                        cluster_reco_index,
                        particle_reco_area_cluster,
                        ue_estimate_cluster);
}

void AliAnalysisTaskNTGJ::getVoronoiAreaTpc(
    std::vector<AliTrackContainer*> track_containers,
    std::vector<fastjet::PseudoJet> &particle_reco_tpc,
    std::map<size_t, size_t> &track_reco_index_tpc,
    std::vector<double> &particle_reco_area_tpc)
{
    // lines 1041-1044, 1163-1171, 1243-1254
    std::vector<size_t> reco_stored_track_index_tpc;            // index for filling voronoi area branch
    std::vector<point_2d_t> particle_reco_area_estimation_tpc;  // used for area calculation
    std::vector<std::set<size_t> > particle_reco_incident_tpc;  // unused variable as part of area calculation

    // fill pseudojet and index vectors
    Int_t itrack = 0;
    for (auto track_container : track_containers) {
        for (auto track : track_container->accepted()) {
            if (track == NULL) {
                continue;
            }

            AliAODTrack * t = static_cast<AliAODTrack*>(track);

            if (trackPassesTPCCuts(t)) {
                track_reco_index_tpc[itrack] = particle_reco_tpc.size();
                reco_stored_track_index_tpc.push_back(itrack);
                particle_reco_tpc.push_back(fastjet::PseudoJet(
                                                t->Px(), t->Py(), t->Pz(), t->P()));
            }

            itrack++;
        }
    }

    // calculate and fill voronoi area branch
    for (std::vector<fastjet::PseudoJet>::const_iterator iterator =
                particle_reco_tpc.begin();
            iterator != particle_reco_tpc.end(); iterator++) {
        particle_reco_area_estimation_tpc.push_back(
            point_2d_t(iterator->pseudorapidity(),
                       iterator->phi_std()));
    }

    voronoi_area_incident(particle_reco_area_tpc,
                          particle_reco_incident_tpc,
                          particle_reco_area_estimation_tpc);

    for (size_t i = 0; i < particle_reco_area_tpc.size(); i++) {
        if (i < reco_stored_track_index_tpc.size() &&
                reco_stored_track_index_tpc[i] < _branch_ntrack) {
            _branch_track_voronoi_area
            [reco_stored_track_index_tpc[i]] =
                half(particle_reco_area_tpc[i]);
        }
    }
}

void AliAnalysisTaskNTGJ::getVoronoiAreaIts(
    std::vector<AliTrackContainer*> track_containers,
    std::vector<fastjet::PseudoJet> &particle_reco_its,
    std::map<size_t, size_t> &track_reco_index_its,
    std::vector<double> &particle_reco_area_its)
{
    // lines 1051-1055, 1173-1181, 1256-1270
    std::vector<size_t> reco_stored_track_index_its;            // index for filling voronoi area branch
    std::vector<point_2d_t> particle_reco_area_estimation_its;  // used for area calculation
    std::vector<std::set<size_t> > particle_reco_incident_its;  // unused variable as part of area calculation

    // fill pseudojet and index vectors
    Int_t itrack = 0;
    for (auto track_container : track_containers) {
        for (auto track : track_container->accepted()) {
            if (track == NULL) {
                continue;
            }

            AliAODTrack * t = static_cast<AliAODTrack*>(track);

            // eventually switch this back to trackPassesCut4 when that doesn't need the event itself anymore
            if ((_branch_track_quality[itrack] & 16) != 0) {
                track_reco_index_its[itrack] = particle_reco_its.size();
                reco_stored_track_index_its.push_back(itrack);
                particle_reco_its.push_back(fastjet::PseudoJet(
                                                t->Px(), t->Py(), t->Pz(), t->P()));
            }

            itrack++;
        }
    }

    // calculate and fill voronoi area branch
    for (std::vector<fastjet::PseudoJet>::const_iterator iterator =
                particle_reco_its.begin();
            iterator != particle_reco_its.end(); iterator++) {
        particle_reco_area_estimation_its.push_back(
            point_2d_t(iterator->pseudorapidity(),
                       iterator->phi_std()));
    }

    voronoi_area_incident(particle_reco_area_its,
                          particle_reco_incident_its,
                          particle_reco_area_estimation_its);

    for (size_t i = 0; i < particle_reco_area_its.size(); i++) {
        if (i < reco_stored_track_index_its.size() &&
                reco_stored_track_index_its[i] < _branch_ntrack) {
            _branch_track_voronoi_area
            [reco_stored_track_index_its[i]] =
                half(particle_reco_area_its[i]);
        }
    }
}

void AliAnalysisTaskNTGJ::getVoronoiAreaCluster(
    AliClusterContainer *cluster_container,
    std::vector<bool> cell_pass_basic_quality,
    std::vector<fastjet::PseudoJet> &particle_reco_cluster,
    std::map<size_t, size_t> &cluster_reco_index,
    std::vector<double> &particle_reco_area_cluster)
{
    // lines 1184-1238, 1272-1279, 1591-1593
    std::vector<size_t> reco_stored_cluster_index;                  // unused after filling
    std::vector<point_2d_t> particle_reco_area_estimation_cluster;  // used for area calculation
    std::vector<std::set<size_t> > particle_reco_incident_cluster;  // unused variable as part of area calculation

    // fill pseudojet and index vectors
    Int_t icluster = 0;
    for (auto c : cluster_container->accepted()) {
        if (cell_pass_basic_quality[icluster]) {
            TLorentzVector p;
            c->GetMomentum(p, _branch_primary_vertex);

            const fastjet::PseudoJet
            pj(p.Px(), p.Py(), p.Pz(), p.P());

            // FIXME: This needs to be moved somewhere in emcal.h
            static const double pseudorapidity_limit = 0.661;
            static const double azimuth_limit_0 = 1.415;
            static const double azimuth_limit_1 = 3.123;

            // Fill only when inside the EMCAL (not DCAL)
            if (fabs(pj.pseudorapidity()) < pseudorapidity_limit &&
                    pj.phi_std() >= azimuth_limit_0 &&
                    pj.phi_std() < azimuth_limit_1) {
                cluster_reco_index[icluster] = particle_reco_cluster.size();
                // Note that all clusters are stored, at the moment,
                // so this stored index is identical to i.
                reco_stored_cluster_index.push_back(icluster);
                particle_reco_cluster.push_back(pj);
            }
        }
        icluster++;
    }

    // calculate voronoi area
    for (std::vector<fastjet::PseudoJet>::const_iterator iterator =
                particle_reco_cluster.begin();
            iterator != particle_reco_cluster.end(); iterator++) {
        particle_reco_area_estimation_cluster.push_back(
            point_2d_t(iterator->pseudorapidity(),
                       iterator->phi_std()));
    }

    voronoi_area_incident(particle_reco_area_cluster,
                          particle_reco_incident_cluster,
                          particle_reco_area_estimation_cluster,
                          *reinterpret_cast<std::vector<point_2d_t> *>
                          (_emcal_cell_position));

    for (size_t i = 0; i < particle_reco_cluster.size(); i++) {
        particle_reco_cluster[i].set_user_index(static_cast<int>(i));
    }
}

void AliAnalysisTaskNTGJ::getUETpc(
    std::vector<fastjet::PseudoJet> particle_reco_tpc,
    std::vector<double> particle_reco_area_tpc,
    std::pair<std::pair<std::vector<double>, std::vector<double>>, double> &ue_estimate_tpc)
{
    // lines 1596-1603, 1623-1646, 1656-1658

    // kT clustering
    static const double jet_kt_d_ue_estimation = 0.3;
    const fastjet::ClusterSequenceArea
    cluster_sequence_ue_estimation_tpc(
        particle_reco_tpc,
        fastjet::JetDefinition(fastjet::JetDefinition(
                                   fastjet::kt_algorithm, jet_kt_d_ue_estimation)),
        fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_ue_estimation_tpc =
        cluster_sequence_ue_estimation_tpc.inclusive_jets(0);

    // fill "debug" branches
    _branch_debug_njet_ue_estimation = 0;
    for (std::vector<fastjet::PseudoJet>::const_iterator
            iterator_jet = jet_ue_estimation_tpc.begin();
            iterator_jet != jet_ue_estimation_tpc.end();
            iterator_jet++) {
        _branch_debug_jet_ue_estimation_pt_raw
        [_branch_debug_njet_ue_estimation] =
            iterator_jet->perp();
        _branch_debug_jet_ue_estimation_eta_raw
        [_branch_debug_njet_ue_estimation] =
            iterator_jet->pseudorapidity();
        _branch_debug_jet_ue_estimation_phi_raw
        [_branch_debug_njet_ue_estimation] =
            iterator_jet->phi_std();
        _branch_debug_jet_ue_estimation_area_raw
        [_branch_debug_njet_ue_estimation] =
            iterator_jet->area();
        _branch_debug_njet_ue_estimation++;
    }

    // UE estimation
    ue_estimate_tpc = ue_estimation_median(cluster_sequence_ue_estimation_tpc,
                                           particle_reco_area_tpc);
    _branch_ue_estimate_tpc_const = evaluate_ue_constant(ue_estimate_tpc.first);
    _branch_ue_estimate_tpc_const_se = ue_estimate_tpc.second;
}

void AliAnalysisTaskNTGJ::getUEIts(
    std::vector<fastjet::PseudoJet> particle_reco_its,
    std::vector<double> particle_reco_area_its,
    std::pair<std::pair<std::vector<double>, std::vector<double>>, double> &ue_estimate_its)
{
    // lines 1604-1611, 1647-1650, 1659-1661

    // kT clustering
    static const double jet_kt_d_ue_estimation = 0.3;
    const fastjet::ClusterSequenceArea
    cluster_sequence_ue_estimation_its(
        particle_reco_its,
        fastjet::JetDefinition(fastjet::JetDefinition(
                                   fastjet::kt_algorithm, jet_kt_d_ue_estimation)),
        fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_ue_estimation_its =
        cluster_sequence_ue_estimation_its.inclusive_jets(0);

    // UE estimation
    ue_estimate_its = ue_estimation_median(cluster_sequence_ue_estimation_its,
                                           particle_reco_area_its);
    _branch_ue_estimate_its_const = evaluate_ue_constant(ue_estimate_its.first);
    _branch_ue_estimate_its_const_se = ue_estimate_its.second;
}

void AliAnalysisTaskNTGJ::getUECluster(
    std::vector<fastjet::PseudoJet> particle_reco_cluster,
    std::vector<double> particle_reco_area_cluster,
    std::pair<std::pair<std::vector<double>, std::vector<double>>, double> &ue_estimate_cluster)
{
    // lines 1612-1619, 1651-1654, 1662-1664

    // kT clustering
    static const double jet_kt_d_ue_estimation = 0.3;
    const fastjet::ClusterSequenceArea
    cluster_sequence_ue_estimation_cluster(
        particle_reco_cluster,
        fastjet::JetDefinition(fastjet::JetDefinition(
                                   fastjet::kt_algorithm, jet_kt_d_ue_estimation)),
        fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_ue_estimation_cluster =
        cluster_sequence_ue_estimation_cluster.inclusive_jets(0);

    // UE estimation
    ue_estimate_cluster = ue_estimation_median(cluster_sequence_ue_estimation_cluster,
                          particle_reco_area_cluster);

    _branch_ue_estimate_cluster_const = evaluate_ue_constant(ue_estimate_cluster.first);
    _branch_ue_estimate_cluster_const_se = ue_estimate_cluster.second;
}

void AliAnalysisTaskNTGJ::getIsolationTpc(
    AliClusterContainer *cluster_container,
    std::vector<AliTrackContainer*> track_containers,
    std::map<size_t, size_t> track_reco_index_tpc,
    std::vector<double> particle_reco_area_tpc,
    std::pair<std::pair<std::vector<double>, std::vector<double>>, double> ue_estimate_tpc)
{
    UInt_t icluster = 0;
    for (auto cluster : cluster_container->accepted()) {
        AliVCluster *c = cluster;
        TLorentzVector p;
        c->GetMomentum(p, _branch_primary_vertex);

        double cluster_iso_tpc_01 = 0;
        double cluster_iso_tpc_02 = 0;
        double cluster_iso_tpc_03 = 0;
        double cluster_iso_tpc_04 = 0;
        double cluster_iso_tpc_01_ue = 0;
        double cluster_iso_tpc_02_ue = 0;
        double cluster_iso_tpc_03_ue = 0;
        double cluster_iso_tpc_04_ue = 0;

        std::vector<std::pair<double, double> > delta_vs_iso_tpc;
        std::vector<std::pair<double, double> > delta_vs_iso_tpc_with_ue;

        Int_t itrack = 0;
        for (auto track_container : track_containers) {
            for (auto track : track_container->accepted()) {
                if (track == NULL) {
                    continue;
                }

                AliAODTrack * t = static_cast<AliAODTrack*>(track);

                // Apply PWG-JE cuts (track cuts 0 and 1)
                if (trackPassesCut0(t) || trackPassesCut1(t)) {
                    const double dpseudorapidity = t->Eta() - p.Eta();
                    const double dazimuth = angular_range_reduce(
                                                angular_range_reduce(t->Phi()) -
                                                angular_range_reduce(p.Phi()));
                    const double dr_2 =
                        std::pow(dpseudorapidity, 2) +
                        std::pow(dazimuth, 2);
                    const double ue =
                        dr_2 < 0.4 * 0.4 ?
                        track_reco_index_tpc.find(itrack) !=
                        track_reco_index_tpc.end() ?
                        evaluate_ue(ue_estimate_tpc.first, t->Eta(),
                                    t->Phi()) *
                        particle_reco_area_tpc
                        [track_reco_index_tpc[itrack]] :
                        0 : NAN;
                    const double track_pt_minus_ue =
                        dr_2 < 0.4 * 0.4 ?
                        track_reco_index_tpc.find(itrack) !=
                        track_reco_index_tpc.end() ?
                        t->Pt() - ue :
                        0 : NAN;

                    if (dr_2 < 0.1 * 0.1) {
                        cluster_iso_tpc_01 += track_pt_minus_ue;
                        cluster_iso_tpc_01_ue += ue;
                    }
                    if (dr_2 < 0.2 * 0.2) {
                        cluster_iso_tpc_02 += track_pt_minus_ue;
                        cluster_iso_tpc_02_ue += ue;
                    }
                    if (dr_2 < 0.3 * 0.3) {
                        cluster_iso_tpc_03 += track_pt_minus_ue;
                        cluster_iso_tpc_03_ue += ue;
                    }
                    if (dr_2 < 0.4 * 0.4) {
                        cluster_iso_tpc_04 += track_pt_minus_ue;
                        cluster_iso_tpc_04_ue += ue;
                        delta_vs_iso_tpc.push_back(
                            std::pair<double, double>(
                                sqrt(dr_2), track_pt_minus_ue));
                        delta_vs_iso_tpc_with_ue.push_back(
                            std::pair<double, double>(
                                sqrt(dr_2), track_pt_minus_ue + ue));
                    }
                }

                itrack++;
            }
        }

        _branch_cluster_iso_tpc_01[icluster] =
            half(cluster_iso_tpc_01);
        _branch_cluster_iso_tpc_02[icluster] =
            half(cluster_iso_tpc_02);
        _branch_cluster_iso_tpc_03[icluster] =
            half(cluster_iso_tpc_03);
        _branch_cluster_iso_tpc_04[icluster] =
            half(cluster_iso_tpc_04);
        _branch_cluster_iso_tpc_01_ue[icluster] =
            half(cluster_iso_tpc_01_ue);
        _branch_cluster_iso_tpc_02_ue[icluster] =
            half(cluster_iso_tpc_02_ue);
        _branch_cluster_iso_tpc_03_ue[icluster] =
            half(cluster_iso_tpc_03_ue);
        _branch_cluster_iso_tpc_04_ue[icluster] =
            half(cluster_iso_tpc_04_ue);

        _branch_cluster_frixione_tpc_04_02[icluster] =
            half(frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                          0.4, 0.2));
        _branch_cluster_frixione_tpc_04_05[icluster] =
            half(frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                          0.4, 0.5));
        _branch_cluster_frixione_tpc_04_10[icluster] =
            half(frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                          0.4, 1.0));
        _branch_cluster_frixione_tpc_04_02_with_ue[icluster] =
            half(frixione_iso_max_x_e_eps(
                     delta_vs_iso_tpc_with_ue, 0.4, 0.2));
        _branch_cluster_frixione_tpc_04_05_with_ue[icluster] =
            half(frixione_iso_max_x_e_eps(
                     delta_vs_iso_tpc_with_ue, 0.4, 0.5));
        _branch_cluster_frixione_tpc_04_10_with_ue[icluster] =
            half(frixione_iso_max_x_e_eps(
                     delta_vs_iso_tpc_with_ue, 0.4, 1.0));
        _branch_cluster_anti_frixione_tpc_04_02[icluster] =
            half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                               0.4, 0.2));
        _branch_cluster_anti_frixione_tpc_04_05[icluster] =
            half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                               0.4, 0.5));
        _branch_cluster_anti_frixione_tpc_04_10[icluster] =
            half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                               0.4, 1.0));

        icluster++;
        if (icluster >= NCLUSTER_MAX) {
            break;
        }
    }
}

void AliAnalysisTaskNTGJ::getIsolationIts(
    AliClusterContainer *cluster_container,
    std::vector<AliTrackContainer*> track_containers,
    std::map<size_t, size_t> track_reco_index_its,
    std::vector<double> particle_reco_area_its,
    std::pair<std::pair<std::vector<double>, std::vector<double>>, double> ue_estimate_its)
{
    UInt_t icluster = 0;
    for (auto cluster : cluster_container->accepted()) {
        AliVCluster *c = cluster;
        TLorentzVector p;
        c->GetMomentum(p, _branch_primary_vertex);

        double cluster_iso_its_01 = 0;
        double cluster_iso_its_02 = 0;
        double cluster_iso_its_03 = 0;
        double cluster_iso_its_04 = 0;
        double cluster_iso_its_01_ue = 0;
        double cluster_iso_its_02_ue = 0;
        double cluster_iso_its_03_ue = 0;
        double cluster_iso_its_04_ue = 0;

        std::vector<std::pair<double, double> > delta_vs_iso_its;
        std::vector<std::pair<double, double> > delta_vs_iso_its_with_ue;

        Int_t itrack = 0;
        for (auto track_container : track_containers) {
            for (auto track : track_container->accepted()) {
                if (track == NULL) {
                    continue;
                }

                AliAODTrack * t = static_cast<AliAODTrack*>(track);

                // eventually switch this back to trackPassesCut4 when that doesn't need the event itself anymore
                if ((_branch_track_quality[itrack] & 16) != 0) {
                    const double dpseudorapidity = t->Eta() - p.Eta();
                    const double dazimuth = angular_range_reduce(
                                                angular_range_reduce(t->Phi()) -
                                                angular_range_reduce(p.Phi()));
                    const double dr_2 =
                        std::pow(dpseudorapidity, 2) +
                        std::pow(dazimuth, 2);
                    const double ue =
                        dr_2 < 0.4 * 0.4 ?
                        track_reco_index_its.find(itrack) !=
                        track_reco_index_its.end() ?
                        evaluate_ue(ue_estimate_its.first, t->Eta(),
                                    t->Phi()) *
                        particle_reco_area_its
                        [track_reco_index_its[itrack]] :
                        0 : NAN;
                    const double track_pt_minus_ue =
                        dr_2 < 0.4 * 0.4 ?
                        track_reco_index_its.find(itrack) !=
                        track_reco_index_its.end() ?
                        t->Pt() - ue :
                        0 : NAN;

                    if (dr_2 < 0.1 * 0.1) {
                        cluster_iso_its_01 += track_pt_minus_ue;
                        cluster_iso_its_01_ue += ue;
                    }
                    if (dr_2 < 0.2 * 0.2) {
                        cluster_iso_its_02 += track_pt_minus_ue;
                        cluster_iso_its_02_ue += ue;
                    }
                    if (dr_2 < 0.3 * 0.3) {
                        cluster_iso_its_03 += track_pt_minus_ue;
                        cluster_iso_its_03_ue += ue;
                    }
                    if (dr_2 < 0.4 * 0.4) {
                        cluster_iso_its_04 += track_pt_minus_ue;
                        cluster_iso_its_04_ue += ue;
                        delta_vs_iso_its.push_back(
                            std::pair<double, double>(
                                sqrt(dr_2), track_pt_minus_ue));
                        delta_vs_iso_its_with_ue.push_back(
                            std::pair<double, double>(
                                sqrt(dr_2), track_pt_minus_ue + ue));
                    }
                }

                itrack++;
            }
        }

        _branch_cluster_iso_its_01[icluster] =
            half(cluster_iso_its_01);
        _branch_cluster_iso_its_02[icluster] =
            half(cluster_iso_its_02);
        _branch_cluster_iso_its_03[icluster] =
            half(cluster_iso_its_03);
        _branch_cluster_iso_its_04[icluster] =
            half(cluster_iso_its_04);
        _branch_cluster_iso_its_01_ue[icluster] =
            half(cluster_iso_its_01_ue);
        _branch_cluster_iso_its_02_ue[icluster] =
            half(cluster_iso_its_02_ue);
        _branch_cluster_iso_its_03_ue[icluster] =
            half(cluster_iso_its_03_ue);
        _branch_cluster_iso_its_04_ue[icluster] =
            half(cluster_iso_its_04_ue);

        _branch_cluster_frixione_its_04_02[icluster] =
            half(frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                          0.4, 0.2));
        _branch_cluster_frixione_its_04_05[icluster] =
            half(frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                          0.4, 0.5));
        _branch_cluster_frixione_its_04_10[icluster] =
            half(frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                          0.4, 1.0));
        _branch_cluster_frixione_its_04_02_with_ue[icluster] =
            half(frixione_iso_max_x_e_eps(
                     delta_vs_iso_its_with_ue, 0.4, 0.2));
        _branch_cluster_frixione_its_04_05_with_ue[icluster] =
            half(frixione_iso_max_x_e_eps(
                     delta_vs_iso_its_with_ue, 0.4, 0.5));
        _branch_cluster_frixione_its_04_10_with_ue[icluster] =
            half(frixione_iso_max_x_e_eps(
                     delta_vs_iso_its_with_ue, 0.4, 1.0));
        _branch_cluster_anti_frixione_its_04_02[icluster] =
            half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                               0.4, 0.2));
        _branch_cluster_anti_frixione_its_04_05[icluster] =
            half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                               0.4, 0.5));
        _branch_cluster_anti_frixione_its_04_10[icluster] =
            half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                               0.4, 1.0));

        icluster++;
        if (icluster >= NCLUSTER_MAX) {
            break;
        }
    }
}

void AliAnalysisTaskNTGJ::getIsolationCluster(
    AliClusterContainer *cluster_container,
    std::vector<bool> cell_pass_basic_quality,
    std::map<size_t, size_t> cluster_reco_index,
    std::vector<double> particle_reco_area_cluster,
    std::pair<std::pair<std::vector<double>, std::vector<double>>, double> ue_estimate_cluster)
{
    UInt_t icluster = 0;
    for (auto cluster : cluster_container->accepted()) {
        AliVCluster *c = cluster;
        TLorentzVector p;
        c->GetMomentum(p, _branch_primary_vertex);

        double cluster_iso_cluster_01 = 0;
        double cluster_iso_cluster_02 = 0;
        double cluster_iso_cluster_03 = 0;
        double cluster_iso_cluster_04 = 0;
        double cluster_iso_cluster_01_ue = 0;
        double cluster_iso_cluster_02_ue = 0;
        double cluster_iso_cluster_03_ue = 0;
        double cluster_iso_cluster_04_ue = 0;

        UInt_t jcluster = 0;
        for (auto c : cluster_container->accepted()) {
            if (!cell_pass_basic_quality[jcluster]) {
                continue;
            }

            TLorentzVector p1;
            c->GetMomentum(p1, _branch_primary_vertex);

            const double dpseudorapidity = p1.Eta() - p.Eta();
            const double dazimuth = angular_range_reduce(
                                        angular_range_reduce(p1.Phi()) -
                                        angular_range_reduce(p.Phi()));
            const double dr_2 =
                std::pow(dpseudorapidity, 2) +
                std::pow(dazimuth, 2);
            const double ue =
                dr_2 < 0.4 * 0.4 ?
                cluster_reco_index.find(jcluster) !=
                cluster_reco_index.end() ?
                evaluate_ue(ue_estimate_cluster.first, p1.Eta(),
                            p1.Phi()) *
                particle_reco_area_cluster
                [cluster_reco_index[jcluster]] :
                0 : NAN;
            const double cluster_pt_minus_ue =
                dr_2 < 0.4 * 0.4 ?
                cluster_reco_index.find(jcluster) !=
                cluster_reco_index.end() ?
                p1.Pt() - ue :
                0 : NAN;

            if (dr_2 < 0.1 * 0.1) {
                cluster_iso_cluster_01 += cluster_pt_minus_ue;
                cluster_iso_cluster_01_ue += ue;
            }
            if (dr_2 < 0.2 * 0.2) {
                cluster_iso_cluster_02 += cluster_pt_minus_ue;
                cluster_iso_cluster_02_ue += ue;
            }
            if (dr_2 < 0.3 * 0.3) {
                cluster_iso_cluster_03 += cluster_pt_minus_ue;
                cluster_iso_cluster_03_ue += ue;
            }
            if (dr_2 < 0.4 * 0.4) {
                cluster_iso_cluster_04 += cluster_pt_minus_ue;
                cluster_iso_cluster_04_ue += ue;
            }

            jcluster++;
        }

        cluster_iso_cluster_01 -= p.Pt();
        cluster_iso_cluster_02 -= p.Pt();
        cluster_iso_cluster_03 -= p.Pt();
        cluster_iso_cluster_04 -= p.Pt();

        _branch_cluster_iso_cluster_01[icluster] =
            half(cluster_iso_cluster_01);
        _branch_cluster_iso_cluster_02[icluster] =
            half(cluster_iso_cluster_02);
        _branch_cluster_iso_cluster_03[icluster] =
            half(cluster_iso_cluster_03);
        _branch_cluster_iso_cluster_04[icluster] =
            half(cluster_iso_cluster_04);
        _branch_cluster_iso_cluster_01_ue[icluster] =
            half(cluster_iso_cluster_01_ue);
        _branch_cluster_iso_cluster_02_ue[icluster] =
            half(cluster_iso_cluster_02_ue);
        _branch_cluster_iso_cluster_03_ue[icluster] =
            half(cluster_iso_cluster_03_ue);
        _branch_cluster_iso_cluster_04_ue[icluster] =
            half(cluster_iso_cluster_04_ue);

        icluster++;
        if (icluster >= NCLUSTER_MAX) {
            break;
        }
    }
}

void AliAnalysisTaskNTGJ::getJetsTpc(
    std::vector<fastjet::PseudoJet> particle_reco_tpc,
    std::vector<fastjet::PseudoJet> particle_reco_tagged_ak04tpc,
    std::vector<fastjet::PseudoJet> jet_truth_ak04,
    std::vector<fastjet::PseudoJet> jet_charged_truth_ak04,
    fastjet::ClusterSequenceArea *cluster_sequence_truth,
    std::vector<double> particle_reco_area_tpc,
    std::pair<std::pair<std::vector<double>, std::vector<double>>, double> ue_estimate_tpc,
    AliClusterContainer *cluster_container,
    std::vector<bool> cell_pass_basic_quality,
    AliMCParticleContainer *mc_container,
    std::vector<Int_t> reverse_stored_parton_algorithmic_index)
{
    // lines 1551-1556, 1522-1523, 1529-1530, 1535-1536, 1551-1567, 2322-2327

    static const double jet_antikt_d_04 = 0.4;
    static const double scale_ghost = pow(2.0, -30.0);

    // add ghost particles to reco tagged vector
    UInt_t icluster = 0;
    for (auto c : cluster_container->accepted()) {
        if (cell_pass_basic_quality[icluster]) {
            TLorentzVector p;
            c->GetMomentum(p, _branch_primary_vertex);

            const fastjet::PseudoJet
            pj(p.Px(), p.Py(), p.Pz(), p.P());

            particle_reco_tagged_ak04tpc.push_back(pj * scale_ghost);
            particle_reco_tagged_ak04tpc.back().
            set_user_index(USER_INDEX_EM);
        }
        icluster++;
    }

    TAG_PARTICLE_RECO_JET_TRUTH(particle_reco_tagged_ak04tpc,
                                jet_truth_ak04, 0);
    TAG_PARTICLE_RECO_JET_TRUTH(particle_reco_tagged_ak04tpc,
                                jet_charged_truth_ak04, 1);
    // Insert partons as ghosts into the reco particles
    TAG_PARTICLE_RECO_PARTON(particle_reco_tagged_ak04tpc,
                             algorithmic, ALGORITHMIC);

    // anti-kT clustering
    const fastjet::ClusterSequenceArea
    cluster_sequence_reco_ak04tpc(
        particle_reco_tpc,
        fastjet::JetDefinition(fastjet::JetDefinition(
                                   fastjet::antikt_algorithm, jet_antikt_d_04)),
        fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_reco_ak04tpc =
        cluster_sequence_reco_ak04tpc.inclusive_jets(0);

    const fastjet::ClusterSequenceArea
    cluster_sequence_reco_tagged_ak04tpc(
        particle_reco_tagged_ak04tpc,
        fastjet::JetDefinition(fastjet::JetDefinition(
                                   fastjet::antikt_algorithm, jet_antikt_d_04)),
        fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_reco_tagged_ak04tpc =
        cluster_sequence_reco_tagged_ak04tpc.inclusive_jets(0);

    // fill many branches
    FILL_BRANCH_JET(ak04tpc, jet_reco_ak04tpc,
                    cluster_sequence_reco_ak04tpc,
                    jet_reco_tagged_ak04tpc,
                    cluster_sequence_reco_tagged_ak04tpc,
                    ak04, jet_truth_ak04, jet_charged_truth_ak04,
                    particle_reco_area_tpc, ue_estimate_tpc);
}

void AliAnalysisTaskNTGJ::getJetsIts(
    std::vector<fastjet::PseudoJet> particle_reco_its,
    std::vector<fastjet::PseudoJet> particle_reco_tagged_ak04its,
    std::vector<fastjet::PseudoJet> jet_truth_ak04,
    std::vector<fastjet::PseudoJet> jet_charged_truth_ak04,
    fastjet::ClusterSequenceArea *cluster_sequence_truth,
    std::vector<double> particle_reco_area_its,
    std::pair<std::pair<std::vector<double>, std::vector<double>>, double> ue_estimate_its,
    AliClusterContainer *cluster_container,
    std::vector<bool> cell_pass_basic_quality,
    AliMCParticleContainer *mc_container,
    std::vector<Int_t> reverse_stored_parton_algorithmic_index)
{
    // lines 1573-1580, 1524-1525, 1531-1532, 1537-1538, 1573-1589, 2328-2333

    static const double jet_antikt_d_04 = 0.4;
    static const double scale_ghost = pow(2.0, -30.0);

    // add ghost particles to reco tagged vector
    UInt_t icluster = 0;
    for (auto c : cluster_container->accepted()) {
        if (cell_pass_basic_quality[icluster]) {
            TLorentzVector p;
            c->GetMomentum(p, _branch_primary_vertex);

            const fastjet::PseudoJet
            pj(p.Px(), p.Py(), p.Pz(), p.P());

            particle_reco_tagged_ak04its.push_back(pj * scale_ghost);
            particle_reco_tagged_ak04its.back().
            set_user_index(USER_INDEX_EM);
        }
        icluster++;
    }

    TAG_PARTICLE_RECO_JET_TRUTH(particle_reco_tagged_ak04its,
                                jet_truth_ak04, 0);
    TAG_PARTICLE_RECO_JET_TRUTH(particle_reco_tagged_ak04its,
                                jet_charged_truth_ak04, 1);
    // Insert partons as ghosts into the reco particles
    TAG_PARTICLE_RECO_PARTON(particle_reco_tagged_ak04its,
                             algorithmic, ALGORITHMIC);

    // anti-kT clustering
    const fastjet::ClusterSequenceArea
    cluster_sequence_reco_ak04its(
        particle_reco_its,
        fastjet::JetDefinition(fastjet::JetDefinition(
                                   fastjet::antikt_algorithm, jet_antikt_d_04)),
        fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_reco_ak04its =
        cluster_sequence_reco_ak04its.inclusive_jets(0);

    const fastjet::ClusterSequenceArea
    cluster_sequence_reco_tagged_ak04its(
        particle_reco_tagged_ak04its,
        fastjet::JetDefinition(fastjet::JetDefinition(
                                   fastjet::antikt_algorithm, jet_antikt_d_04)),
        fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_reco_tagged_ak04its =
        cluster_sequence_reco_tagged_ak04its.inclusive_jets(0);

    // fill many branches
    FILL_BRANCH_JET(ak04its, jet_reco_ak04its,
                    cluster_sequence_reco_ak04its,
                    jet_reco_tagged_ak04its,
                    cluster_sequence_reco_tagged_ak04its,
                    ak04, jet_truth_ak04, jet_charged_truth_ak04,
                    particle_reco_area_its, ue_estimate_its);
}

bool AliAnalysisTaskNTGJ::skimJets()
{
    // lines 2335-2387
    if (_skim_jet_min_pt[0] > -INFINITY) {
        // FIXME: JEC?
        std::vector<float>
        pt(_branch_jet_ak04tpc_pt_raw,
           _branch_jet_ak04tpc_pt_raw + _branch_njet_ak04tpc);
        if (pt.size() >= 3) {
            std::partial_sort(pt.begin(), pt.begin() + 3, pt.end(),
                              std::greater<float>());
        }
        else {
            std::sort(pt.begin(), pt.end(), std::greater<float>());
        }
        if (_skim_jet_min_pt.size() >= 3 &&
                (pt.size() < 3 ||
                 !(pt[0] >= _skim_jet_min_pt[0] &&
                   pt[1] >= _skim_jet_min_pt[1] &&
                   pt[2] >= _skim_jet_min_pt[2]))) {
            // Discard this event
            return false;
        }
        else if (_skim_jet_min_pt.size() >= 2 &&
                 (pt.size() < 2 ||
                  !(pt[0] >= _skim_jet_min_pt[0] &&
                    pt[1] >= _skim_jet_min_pt[1]))) {
            // Discard this event
            return false;
        }
        else if (_skim_jet_min_pt.size() >= 1 &&
                 (pt.size() < 1 ||
                  !(pt[0] >= _skim_jet_min_pt[0]))) {
            // Discard this event
            return false;
        }
    }

    if (_skim_jet_average_pt > -INFINITY) {
        // FIXME: JEC?
        std::vector<float>
        pt(_branch_jet_ak04tpc_pt_raw,
           _branch_jet_ak04tpc_pt_raw + _branch_njet_ak04tpc);
        if (pt.size() >= 2) {
            std::partial_sort(pt.begin(), pt.begin() + 2, pt.end(),
                              std::greater<float>());
        }
        else {
            std::sort(pt.begin(), pt.end(), std::greater<float>());
        }
        if (!(pt.size() >= 2 &&
                0.5 * (pt[0] + pt[1]) >= _skim_jet_average_pt)) {
            // Discard this event
            return false;
        }
    }

    return true;
}

void AliAnalysisTaskNTGJ::fillCellBranches(AliVCaloCells *emcal_cell,
        std::vector<size_t> stored_mc_truth_index)
{
    // lines 891-942, 2389-2417
    std::fill(&_branch_cell_position[0][0],
              &_branch_cell_position[0][0] +
              sizeof(_branch_cell_position) /
              sizeof(_branch_cell_position[0][0]), NAN);
    std::fill(_branch_cell_voronoi_area,
              _branch_cell_voronoi_area +
              sizeof(_branch_cell_voronoi_area) /
              sizeof(*_branch_cell_voronoi_area), NAN);

    if (_emcal_geometry != NULL) {
        // this fixes the memory leak
        if (_emcal_cell_position == NULL) {
            _emcal_cell_position = new std::vector<point_2d_t>();
        }
        else {
            reinterpret_cast<std::vector<point_2d_t> *>
            (_emcal_cell_position)->clear();
        }

        for (unsigned int cell_id = 0; cell_id < EMCAL_NCELL;
                cell_id++) {
            TVector3 v;

            _emcal_geometry->GetGlobal(cell_id, v);
            v -= TVector3(_branch_primary_vertex);

            _branch_cell_position[cell_id][0] = v.X();
            _branch_cell_position[cell_id][1] = v.Y();
            _branch_cell_position[cell_id][2] = v.Z();
            reinterpret_cast<std::vector<point_2d_t> *>
            (_emcal_cell_position)->push_back(
                point_2d_t(v.Eta(), v.Phi()));
        }
        // //FIXME: check to see if "computeVoronoiAreas" is just for jet area (is likeley)
        // std::vector<double> emcal_cell_area;
        // std::vector<std::set<size_t> > emcal_cell_incident;

        // voronoi_area_incident(emcal_cell_area, emcal_cell_incident,
        //                       *reinterpret_cast<std::vector<point_2d_t> *>
        //                       (_emcal_cell_position));

        // double sum_area_inside = 0;
        // size_t count_inside = 0;

        // for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
        //     if (inside_edge(cell_id, 1)) {
        //         sum_area_inside += emcal_cell_area[cell_id];
        //         count_inside++;
        //     }
        // }

        // const double mean_area_inside =
        //     sum_area_inside / count_inside;

        // for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
        //     _branch_cell_voronoi_area[cell_id] =
        //         inside_edge(cell_id, 1) ?
        //         emcal_cell_area[cell_id] : mean_area_inside;
        // }
    }

// why is this here of all places? //No idea....
#ifdef WITH_EFP7
#include "fill_efp7.h"
#endif // WITH_EFP7

    std::fill(_branch_cell_e,
              _branch_cell_e + sizeof(_branch_cell_e) /
              sizeof(*_branch_cell_e), NAN);
    std::fill(_branch_cell_tof,
              _branch_cell_tof + sizeof(_branch_cell_tof) /
              sizeof(*_branch_cell_tof), NAN);
    std::fill(_branch_cell_mc_truth_index,
              _branch_cell_mc_truth_index +
              sizeof(_branch_cell_mc_truth_index) /
              sizeof(*_branch_cell_mc_truth_index), USHRT_MAX);

    for (Short_t i = 0; i < emcal_cell->GetNumberOfCells(); i++) {
        Short_t cell_number = -1;
        Double_t cell_energy = NAN;
        Double_t tof = NAN;
        Int_t mc_truth_index = -1;
        Double_t efrac = NAN;

        if (emcal_cell->GetCell(
                    i, cell_number, cell_energy, tof, mc_truth_index, efrac) ==
                kTRUE &&
                cell_number >= 0 && cell_number < EMCAL_NCELL) {
            _branch_cell_e[cell_number]   = half(cell_energy);
            _branch_cell_tof[cell_number] = half(tof * 1e+9);
            _branch_cell_mc_truth_index[cell_number] =
                SAFE_MC_TRUTH_INDEX_TO_USHRT(mc_truth_index);
        }
    }
}

void AliAnalysisTaskNTGJ::fillMuonBranches()
{

}

void AliAnalysisTaskNTGJ::fillEgNtrial()
{
    // lines 2512-2515
    // Now that the event is accepted, copy over the total counted
    // ntrials so far, and reset the ntrial counter
    _branch_eg_ntrial = _skim_sum_eg_ntrial;
    _skim_sum_eg_ntrial = 0;
}

// line 2419; again, not sure where to put this
#undef SAFE_MC_TRUTH_INDEX_TO_USHRT

AliEMCALRecoUtils *AliAnalysisTaskNTGJ::GetEMCALRecoUtils(void)
{
    return _reco_util;
}

void AliAnalysisTaskNTGJ::SetAliROOTVersion(const char *version)
{
    strncpy(_branch_version_aliroot, version, BUFSIZ);
    _branch_version_aliroot[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::SetAliPhysicsVersion(const char *version)
{
    strncpy(_branch_version_aliphysics, version, BUFSIZ);
    _branch_version_aliphysics[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::SetGridDataDir(const char *dir)
{
    strncpy(_branch_grid_data_dir, dir, BUFSIZ);
    _branch_grid_data_dir[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::SetGridDataPattern(const char *pattern)
{
    strncpy(_branch_grid_data_pattern, pattern, BUFSIZ);
    _branch_grid_data_pattern[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::
SetEMCALGeometryFilename(const char *emcal_geometry_filename)
{
    _emcal_geometry_filename = emcal_geometry_filename;
}

void AliAnalysisTaskNTGJ::
SetEMCALLocal2MasterFilename(const char *emcal_local2master_filename)
{
    _emcal_local2master_filename = emcal_local2master_filename;
}

void AliAnalysisTaskNTGJ::
SetForceUESubtraction(bool force_ue_subtraction)
{
    _force_ue_subtraction = force_ue_subtraction;
}

void AliAnalysisTaskNTGJ::SetSkimClusterMinE(double min_e)
{
    _skim_cluster_min_e = min_e;
}

void AliAnalysisTaskNTGJ::SetSkimTrackMinPt(double min_pt)
{
    _skim_track_min_pt = min_pt;
}

void AliAnalysisTaskNTGJ::SetSkimMuonTrackMinPt(double min_pt)
{
    _skim_muon_track_min_pt = min_pt;
}

void AliAnalysisTaskNTGJ::SetSkimJetMinPt(double min_pt_1,
        double min_pt_2,
        double min_pt_3)
{
    _skim_jet_min_pt.clear();

    if (min_pt_1 != -INFINITY) {
        _skim_jet_min_pt.push_back(min_pt_1);
        if (min_pt_2 != -INFINITY) {
            _skim_jet_min_pt.push_back(min_pt_2);
            if (min_pt_3 != -INFINITY) {
                _skim_jet_min_pt.push_back(min_pt_3);
            }
        }
    }
}

void AliAnalysisTaskNTGJ::SetSkimJetAveragePt(double average_pt)
{
    _skim_jet_average_pt = average_pt;
}

void AliAnalysisTaskNTGJ::SetSkimMultiplicityTrackletMinN(int min_n)
{
    _skim_multiplicity_tracklet_min_n = min_n;
}

void AliAnalysisTaskNTGJ::SetStoredTrackMinPt(double min_pt)
{
    _stored_track_min_pt = min_pt;
}

// This is primarily intended for JEC derivation, and therefore
// defined in pT raw

void AliAnalysisTaskNTGJ::SetStoredJetMinPtRaw(double min_pt_raw)
{
    _stored_jet_min_pt_raw = min_pt_raw;
}

void AliAnalysisTaskNTGJ::SetNRandomIsolation(unsigned int nrandom_isolation)
{
    _nrandom_isolation = nrandom_isolation;
}

void AliAnalysisTaskNTGJ::SetIsEmbed(bool is_embed)
{
    _is_embed = is_embed;
}
