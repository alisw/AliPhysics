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
#include <AliStack.h>
#include <AliAODMCParticle.h>
#include <AliESDMuonTrack.h>

#include <AliMultSelection.h>
#include <AliEventplane.h>

#include <AliMagF.h>

#include "AliAnalysisTaskNtuple.h"

#ifndef __CINT__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
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

ClassImp(AliAnalysisTaskNtuple);

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
    _emcal_mask(std::vector<bool>()),                       \
    _emcal_cell_position(NULL),                             \
    _emcal_cell_area(std::vector<double>()),                \
    _emcal_cell_incident(std::vector<std::set<size_t> >()), \
    _load_intel_mklml(false),                               \
    _libiomp5(NULL), _libmklml_gnu(NULL),                   \
    _keras_model_photon_discrimination(NULL),               \
    _alien_plugin(NULL),                                    \
    _metadata_filled(false)

AliAnalysisTaskNtuple::AliAnalysisTaskNtuple(const char *name)
    : AliAnalysisTaskSE(name), CLASS_INITIALIZATION
{
    DefineOutput(1, TTree::Class());
}

AliAnalysisTaskNtuple::AliAnalysisTaskNtuple(
    const AliAnalysisTaskNtuple &x)
    : AliAnalysisTaskSE(), CLASS_INITIALIZATION
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

AliAnalysisTaskNtuple &AliAnalysisTaskNtuple::operator=(
    const AliAnalysisTaskNtuple &x)
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
}

AliAnalysisTaskNtuple::~AliAnalysisTaskNtuple(void)
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

void AliAnalysisTaskNtuple::UserCreateOutputObjects(void)
{
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

void AliAnalysisTaskNtuple::UserExec(Option_t *option)
{
    if (!_emcal_geometry_filename.empty() &&
        !_emcal_local2master_filename.empty() &&
        _emcal_geometry == NULL) {
        if (!gSystem->
            AccessPathName(gSystem->
                           ExpandPathName(_emcal_geometry_filename.
                                          c_str()))) {
            TGeoManager::Import(_emcal_geometry_filename.c_str());
            _emcal_geometry = AliEMCALGeometry::
                GetInstance(_emcal_geometry_name);
        }

        AliOADBContainer emcal_geometry_container("emcal");

        if (!gSystem->
            AccessPathName(gSystem->
                           ExpandPathName(_emcal_local2master_filename.
                                          c_str()))) {
            emcal_geometry_container.
                InitFromFile(_emcal_local2master_filename.c_str(),
                             "AliEMCALgeo");
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
                        geometry_matrix->At(sm)),
                    sm);
             }
             _branch_has_misalignment_matrix = true;
        }
    }

    if (_load_intel_mklml) {
        if (_libiomp5 == NULL) {
            _libiomp5 = dlopen("libiomp5_so", RTLD_NOW |
                               RTLD_GLOBAL | RTLD_NODELETE);
        }
        _branch_debug_libmklml_gnu_error[0] = '\0';
        if (_libmklml_gnu == NULL) {
            _libmklml_gnu = dlopen("libmklml_gnu_so",
                                   RTLD_NOW | RTLD_NODELETE);
            if (_libmklml_gnu == NULL) {
                snprintf(_branch_debug_libmklml_gnu_error, BUFSIZ,
                         "%s:%d: %s\n", __FILE__, __LINE__,
                         dlerror());
                _branch_debug_libmklml_gnu_loaded = false;
            }
            else {
                _branch_debug_libmklml_gnu_loaded = true;
            }
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

    AliVEvent *event = InputEvent();

    if (event == NULL) {
        return;
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

    AliESDEvent *esd_event = dynamic_cast<AliESDEvent *>(event);
    AliAODEvent *aod_event = NULL;

    if (esd_event == NULL) {
        aod_event = dynamic_cast<AliAODEvent *>(event);
        if (aod_event == NULL) {
            return;
        }
    }

    alice_jec_t jec;

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

    if (_f1_ncluster_tpc_linear_pt_dep == NULL) {
        // PWG-JE cuts, see AddTrackCutsLHC10h() in
        // AliPhysics/PWGJE/macros/AddTaskESDFilterPWGJETrain.C

        _f1_ncluster_tpc_linear_pt_dep =
            new TFormula("_f1_ncluster_tpc_linear_pt_dep",
                         "70.+30./20.*x");

        for (size_t i = 0; i < 2; i++) {
            _track_cut.push_back(AliESDtrackCuts("AliESDtrackCuts"));
            _track_cut.back().
                SetMinNClustersTPCPtDep(_f1_ncluster_tpc_linear_pt_dep,
                                        20.0);
            _track_cut.back().SetMinNClustersTPC(70);
            _track_cut.back().SetMaxChi2PerClusterTPC(4);
            _track_cut.back().SetRequireTPCStandAlone(kTRUE);
            _track_cut.back().SetAcceptKinkDaughters(kFALSE);
            _track_cut.back().SetRequireTPCRefit(kTRUE);
            _track_cut.back().SetMaxFractionSharedTPCClusters(0.4);
            _track_cut.back().SetRequireITSRefit(kTRUE);
            _track_cut.back().SetMaxDCAToVertexXY(2.4);
            _track_cut.back().SetMaxDCAToVertexZ(3.2);
            _track_cut.back().SetDCAToVertex2D(kTRUE);
            _track_cut.back().SetMaxChi2PerClusterITS(36);
            _track_cut.back().SetMaxChi2TPCConstrainedGlobal(36);
            _track_cut.back().SetRequireSigmaToVertex(kFALSE);
            _track_cut.back().SetEtaRange(-0.9, 0.9);
            _track_cut.back().SetPtRange(0.15, 1e+15);

            if (i == 0) {
                _track_cut.back().
                    SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                             AliESDtrackCuts::kAny);
            }
            else {
                _track_cut.back().SetRequireITSRefit(kFALSE);
            }
        }

        // "2015 PbPb" cuts, see
        // AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb() in
        // AliRoot/ANALYSIS/ANALYSISalice/AliESDtrackCuts.cxx . Both
        // clusterCut = 0 or 1 cases are kept, but the options
        // selPrimaries = kTRUE, cutAcceptanceEdges = kTRUE,
        // removeDistortedRegions = kFALSE are fixed

        for (size_t i = 0; i < 2; i++) {
            _track_cut.push_back(AliESDtrackCuts("AliESDtrackCuts"));

            if (i == 0) {
                _track_cut.back().SetMinNClustersTPC(50);
            }
            else {
                _track_cut.back().SetMinNCrossedRowsTPC(70);
                _track_cut.back().
                    SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
            }
            // cutAcceptanceEdges == kTRUE
            _track_cut.back().SetCutGeoNcrNcl(2., 130., 1.5, 0.0, 0.0);
            _track_cut.back().SetMaxChi2PerClusterTPC(4);
            _track_cut.back().SetAcceptKinkDaughters(kFALSE);
            _track_cut.back().SetRequireTPCRefit(kTRUE);
            _track_cut.back().SetRequireITSRefit(kTRUE);
            _track_cut.back().
                    SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                             AliESDtrackCuts::kAny);
            // selPrimaries == kTRUE
            _track_cut.back().
                SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
            // selPrimaries == kTRUE
            _track_cut.back().SetMaxChi2TPCConstrainedGlobal(36);
            _track_cut.back().SetMaxDCAToVertexZ(2);
            _track_cut.back().SetDCAToVertex2D(kFALSE);
            _track_cut.back().SetRequireSigmaToVertex(kFALSE);
            _track_cut.back().SetMaxChi2PerClusterITS(36);
        }

        // Relaxed version of the union of
        // AliESDtrackCuts::GetStandardITSPureSATrackCuts2009() and
        // AliESDtrackCuts::GetStandardITSPureSATrackCuts2010() in
        // AliRoot/ANALYSIS/ANALYSISalice/AliESDtrackCuts.cxx for ITS
        // + SDD tracks

        _track_cut.push_back(AliESDtrackCuts("AliESDtrackCuts"));

        _track_cut.back().SetRequireITSPureStandAlone(kTRUE);
        _track_cut.back().SetPtRange(0.15, 1e+15);
        _track_cut.back().SetMinNClustersITS(5);
        _track_cut.back().SetMaxDCAToVertexXY(2.4);
        _track_cut.back().SetMaxDCAToVertexZ(3.2);
        _track_cut.back().SetDCAToVertex2D(kTRUE);
        _track_cut.back().SetMaxChi2PerClusterITS(36);
    }

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

    std::fill(_branch_centrality, _branch_centrality + 9, NAN);
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
              _branch_event_plane_psi_v0 + 3, NAN);

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

    if (_keras_model_photon_discrimination == NULL) {
        _keras_model_photon_discrimination = new KerasModel;
        reinterpret_cast<KerasModel *>(
            _keras_model_photon_discrimination)->
        LoadModel("photon_discr.model");
    }

    AliMCEvent *mc_truth_event = MCEvent();

    if (mc_truth_event != NULL) {
        mc_truth_event->PreReadAll();
    }

    _branch_eg_signal_process_id = INT_MIN;
    _branch_eg_mpi = INT_MIN;
    _branch_eg_pt_hat = NAN;
    _branch_eg_cross_section = NAN;
    _branch_eg_weight = NAN;
    // This should be default to zero to avoid counting mishap
    _branch_eg_ntrial = 0;

    // Not stored by ALICE SW

    _branch_eg_scale_pdf = NAN;
    _branch_eg_alpha_qcd = NAN;
    _branch_eg_alpha_qed = NAN;
    std::fill(_branch_eg_pdf_id, _branch_eg_pdf_id + 2, INT_MIN);
    std::fill(_branch_eg_pdf_x, _branch_eg_pdf_x + 2, NAN);
    std::fill(_branch_eg_pdf_x_pdf, _branch_eg_pdf_x_pdf + 2, NAN);

    // FIXME: Weight is missing, AliGenEventHeader::EventWeight()

    AliGenEventHeader *mc_truth_header = mc_truth_event != NULL ?
        mc_truth_event->GenEventHeader() : NULL;
    AliGenPythiaEventHeader *mc_truth_pythia_header;

    if (mc_truth_header != NULL) {
        _branch_eg_weight = mc_truth_header->EventWeight();

        TArrayF eg_primary_vertex(3);

        mc_truth_header->PrimaryVertex(eg_primary_vertex);

        for (Int_t i = 0; i < 3; i++) {
            _branch_eg_primary_vertex[i] = eg_primary_vertex.At(i);
        }
        mc_truth_pythia_header =
            dynamic_cast<AliGenPythiaEventHeader *>(mc_truth_header);
        if (mc_truth_pythia_header != NULL) {
            _branch_eg_signal_process_id =
                mc_truth_pythia_header->ProcessType();
            _branch_eg_mpi = mc_truth_pythia_header->GetNMPI();
            _branch_eg_pt_hat =
                mc_truth_pythia_header->GetPtHard();
            _branch_eg_cross_section =
                mc_truth_pythia_header->GetXsection();
            // Count ntrial, because the event might get skimmed away
            // by the ntuplizer
            _skim_sum_eg_ntrial +=
                mc_truth_pythia_header->Trials();
        }
    }

    AliStack *stack;

    if (mc_truth_event != NULL) {
        stack = mc_truth_event->Stack();
    }

    if (esd_event != NULL) {
        esd_event->InitMagneticField();
    }
    else if (aod_event != NULL) {
        aod_event->InitMagneticField();
    }

    std::fill(_branch_primary_vertex,
              _branch_primary_vertex + 3, NAN);
    std::fill(_branch_primary_vertex_sigma,
              _branch_primary_vertex_sigma + 3, NAN);
    _branch_primary_vertex_ncontributor = INT_MIN;

    const AliVVertex *primary_vertex = event->GetPrimaryVertex();

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

    std::fill(_branch_primary_vertex_spd,
              _branch_primary_vertex_spd + 3, NAN);
    std::fill(_branch_primary_vertex_spd_sigma,
              _branch_primary_vertex_spd_sigma + 3, NAN);
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
    _branch_event_selected = fInputHandler->IsEventSelected();

    if (_skim_multiplicity_tracklet_min_n > INT_MIN) {
        AliVMultiplicity *multiplicity = event->GetMultiplicity();

        if (multiplicity == NULL ||
            !(multiplicity->GetNumberOfTracklets() >=
              _skim_multiplicity_tracklet_min_n)) {
            return;
        }
    }

    TRefArray calo_cluster;

    event->GetEMCALClusters(&calo_cluster);

    if (_skim_cluster_min_e > -INFINITY) {
        double cluster_e_max = -INFINITY;

        for (Int_t i = 0; i < calo_cluster.GetEntriesFast(); i++) {
            AliVCluster *c =
                static_cast<AliVCluster *>(calo_cluster.At(i));

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
            return;
        }
    }

    if (!_metadata_filled) {
        // Use gitattributes ident mechanism to track the blob object
        // name
        strncpy(_branch_id_git, "$Id: 870eb10aa49957c38be6b311c8f24a8c762b161d $", BUFSIZ);
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
                      _branch_beam_particle + 2, 0);
        }

        std::fill(_branch_cell_eta, _branch_cell_eta + EMCAL_NCELL,
                  NAN);
        std::fill(_branch_cell_phi, _branch_cell_phi + EMCAL_NCELL,
                  NAN);
        if (_emcal_geometry != NULL) {
            _emcal_cell_position = new std::vector<point_2d_t>();
            for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
                TVector3 v;

                _emcal_geometry->GetGlobal(cell_id, v);
                _branch_cell_eta[cell_id] = v.Eta();
                _branch_cell_phi[cell_id] = v.Phi();
                reinterpret_cast<std::vector<point_2d_t> *>
                    (_emcal_cell_position)->push_back(
                        point_2d_t(v.Eta(), v.Phi()));
            }

            voronoi_area_incident(
                _emcal_cell_area, _emcal_cell_incident,
                *reinterpret_cast<std::vector<point_2d_t> *>
                (_emcal_cell_position));

            double sum_area_inside = 0;
            size_t count_inside = 0;

            for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
                if (inside_edge(cell_id, 1)) {
                    sum_area_inside += _emcal_cell_area[cell_id];
                    count_inside++;
                }
            }

            const double mean_area_inside =
                sum_area_inside / count_inside;

            for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
                if (!inside_edge(cell_id, 1)) {
                    _emcal_cell_area[cell_id] = mean_area_inside;
                }
                _branch_cell_voronoi_area[cell_id] =
                    _emcal_cell_area[cell_id];
            }
        }

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
                  _branch_beam_particle + 2, 0);
        _branch_trigger_class.clear();

        for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
            _branch_cell_eta[cell_id] = NAN;
            _branch_cell_phi[cell_id] = NAN;
            _branch_cell_voronoi_area[cell_id] = NAN;
        }
    }

    std::vector<size_t> stored_mc_truth_index;
    std::vector<Int_t> reverse_stored_mc_truth_index;
    std::vector<Int_t> reverse_stored_parton_algorithmic_index;

    if (mc_truth_event != NULL) {
        stored_mc_truth_index.resize(
            mc_truth_event->GetNumberOfTracks(), ULONG_MAX);

        size_t nmc_truth = 0;

        for (Int_t i = 0;
             i < mc_truth_event->GetNumberOfTracks(); i++) {
            // Bookkeeping for primary final state particles
            if (final_state_primary(mc_truth_event, i)) {
                stored_mc_truth_index[i] = nmc_truth;
                reverse_stored_mc_truth_index.push_back(i);
                nmc_truth++;
            }
            // Bookkeeping for partons
            if (parton_cms_algorithmic(mc_truth_event, i)) {
                reverse_stored_parton_algorithmic_index.push_back(i);
            }
        }

        // Assign secondaries to the primaries

        for (Int_t i = 0;
             i < mc_truth_event->GetNumberOfTracks(); i++) {
            // Skip primaries
            if (final_state_primary(mc_truth_event, i)) {
                continue;
            }

            Int_t j = i;
            bool has_physical_primary_ancestor = false;

            // Ensure termination even if there is a loop
            for (size_t k = 0; k < 1000; k++) {
                const AliMCParticle *p =
                    static_cast<AliMCParticle *>(
                        mc_truth_event->GetTrack(j));

                if (p == NULL || p->Particle() == NULL) {
                    break;
                }
                j = p->Particle()->GetFirstMother();
                if (!(j >= 0 &&
                      j < mc_truth_event->GetNumberOfTracks())) {
                    break;
                }
                if (final_state_primary(mc_truth_event, j)) {
                    has_physical_primary_ancestor = true;
                    break;
                }
            }
            if (has_physical_primary_ancestor) {
                stored_mc_truth_index[i] = stored_mc_truth_index[j];
            }
        }
    }

    std::vector<fastjet::PseudoJet> particle_reco_tpc;
    std::vector<fastjet::PseudoJet> particle_reco_its;

    std::fill(_branch_met_tpc, _branch_met_tpc + 2, 0);
    std::fill(_branch_met_its, _branch_met_its + 2, 0);

    std::map<size_t, size_t> track_reco_index_tpc;
    std::vector<size_t> reco_stored_track_index_tpc;
    std::map<size_t, size_t> track_reco_index_its;
    std::vector<size_t> reco_stored_track_index_its;

    double met_tpc_kahan_error[2] = { 0, 0 };
    double met_its_kahan_error[2] = { 0, 0 };

    _branch_ntrack = 0;
    if (esd_event != NULL) {
        for (Int_t i = 0; i < esd_event->GetNumberOfTracks(); i++) {
            AliESDtrack *t = esd_event->GetTrack(i);

            if (t == NULL) {
                continue;
            }

            // Apply PWG-JE cuts (track cuts 0 and 1)

            if (_track_cut[0].AcceptTrack(t) ||
                _track_cut[1].AcceptTrack(t)) {
                track_reco_index_tpc[i] = particle_reco_tpc.size();
                reco_stored_track_index_tpc.push_back(_branch_ntrack);
                particle_reco_tpc.push_back(fastjet::PseudoJet(
                    t->Px(), t->Py(), t->Pz(), t->P()));
                kahan_sum(_branch_met_tpc[0], met_tpc_kahan_error[0],
                          t->Px());
                kahan_sum(_branch_met_tpc[1], met_tpc_kahan_error[1],
                          t->Py());
            }

            // Apply ITS only cut (track cut 4)

            if (_track_cut[4].AcceptTrack(t)) {
                track_reco_index_its[i] = particle_reco_its.size();
                reco_stored_track_index_its.push_back(_branch_ntrack);
                particle_reco_its.push_back(fastjet::PseudoJet(
                    t->Px(), t->Py(), t->Pz(), t->P()));
                kahan_sum(_branch_met_its[0], met_its_kahan_error[0],
                          t->Px());
                kahan_sum(_branch_met_its[1], met_its_kahan_error[1],
                          t->Py());
            }

            // Store tracks passing PWG-JE *or* "2015 PbPb" cuts

            bool store_track = false;

            for (std::vector<AliESDtrackCuts>::iterator iterator =
                     _track_cut.begin();
                 iterator != _track_cut.end(); iterator++) {
                if (iterator->AcceptTrack(t)) {
                    store_track = true;
                    break;
                }
            }

            if (store_track) {
                _branch_track_e[_branch_ntrack] = half(t->E());
                _branch_track_pt[_branch_ntrack] = half(t->Pt());
                _branch_track_eta[_branch_ntrack] = half(t->Eta());
                _branch_track_phi[_branch_ntrack] =
                    half(angular_range_reduce(t->Phi()));
                if (gGeoManager != NULL) {
                    _branch_track_eta_emcal[_branch_ntrack] =
                        half(t->GetTrackEtaOnEMCal());
                    _branch_track_phi_emcal[_branch_ntrack] =
                        half(angular_range_reduce(
                            t->GetTrackPhiOnEMCal()));
                }
                else {
                    _branch_track_eta_emcal[_branch_ntrack] = NAN;
                    _branch_track_phi_emcal[_branch_ntrack] = NAN;
                }
                _branch_track_charge[_branch_ntrack] =
                    std::min(static_cast<Short_t>(CHAR_MAX),
                             std::max(static_cast<Short_t>(CHAR_MIN),
                                      t->Charge()));

                // Shortened track quality bit mask. Here bit 0 and 1
                // are the PWG-JE's bit 4 and 8. Test for
                // ((track_quality[i] & 3) != 0), i being the track
                // index, to get PWG-JE's "272" (= 1 << 4 | 1 << 8)
                // cut. Test for ((track_quality[i] & 4) == 0) and
                // ((track_quality[i] & 8) == 0) for the "2015 PbPb" cut
                // with clusterCut = 0 and 1.

                _branch_track_quality[_branch_ntrack] = 0U;
                for (size_t j = 0; j < _track_cut.size(); j++) {
                    _branch_track_quality[_branch_ntrack] |=
                        _track_cut[j].AcceptTrack(t) ? 1U << j : 0;
                }

                _branch_track_tpc_dedx[_branch_ntrack] =
                    half(t->GetTPCsignal());

                static const Int_t mode_inner_wall = 1;
                static const Double_t dead_zone_width_cm = 2;
                static const Double_t max_z_cm = 220;

                _branch_track_tpc_length_active_zone
                    [_branch_ntrack] = NAN;
                if (t->GetInnerParam() != NULL) {
                    _branch_track_tpc_length_active_zone
                        [_branch_ntrack] =
                        half(t->GetLengthInActiveZone
                             (mode_inner_wall, dead_zone_width_cm,
                              max_z_cm, event->GetMagneticField()));
                }
                _branch_track_tpc_xrow[_branch_ntrack] =
                    std::min(static_cast<Float_t>(UCHAR_MAX),
                             std::max(0.0F, t->GetTPCCrossedRows()));
                _branch_track_tpc_ncluster[_branch_ntrack] =
                    std::min(UCHAR_MAX,
                             std::max(0, t->GetNumberOfTPCClusters()));
                _branch_track_tpc_ncluster_dedx[_branch_ntrack] =
                    std::min(static_cast<UShort_t>(UCHAR_MAX),
                             std::max(static_cast<UShort_t>(0),
                                      t->GetTPCsignalN()));
                _branch_track_tpc_ncluster_findable[_branch_ntrack] =
                    std::min(static_cast<UShort_t>(UCHAR_MAX),
                             std::max(static_cast<UShort_t>(0),
                                      t->GetTPCNclsF()));
                _branch_track_its_ncluster[_branch_ntrack] =
                    t->GetNumberOfITSClusters();
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

#define SAFE_MC_TRUTH_INDEX_TO_USHRT(mc_truth_index)            \
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

                _branch_track_mc_truth_index[_branch_ntrack] =
                    SAFE_MC_TRUTH_INDEX_TO_USHRT(mc_truth_index);
                _branch_track_voronoi_area[_branch_ntrack] = 0;

                _branch_ntrack++;
                if (_branch_ntrack >= NTRACK_MAX) {
                    break;
                }
            }
        }
    }
    // FIXME: AOD not handled

    std::vector<point_2d_t> particle_reco_area_estimation_tpc;

    for (std::vector<fastjet::PseudoJet>::const_iterator iterator =
             particle_reco_tpc.begin();
         iterator != particle_reco_tpc.end(); iterator++) {
        particle_reco_area_estimation_tpc.push_back(
            point_2d_t(iterator->pseudorapidity(),
                       iterator->phi_std()));
    }

    std::vector<point_2d_t> particle_reco_area_estimation_its;

    for (std::vector<fastjet::PseudoJet>::const_iterator iterator =
             particle_reco_its.begin();
         iterator != particle_reco_its.end(); iterator++) {
        particle_reco_area_estimation_its.push_back(
            point_2d_t(iterator->pseudorapidity(),
                       iterator->phi_std()));
    }

    std::vector<double> particle_reco_area_tpc;
    std::vector<std::set<size_t> > particle_reco_incident_tpc;

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

    std::vector<double> particle_reco_area_its;
    std::vector<std::set<size_t> > particle_reco_incident_its;

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

    // Shared by the isolation and jet code

    enum {
        BEAM_PARTICLE_P = 1001
    };

    const bool subtract_ue =
        _force_ue_subtraction ? true :
        esd_event != NULL ?
        !(esd_event->GetBeamParticle(0) == BEAM_PARTICLE_P &&
          esd_event->GetBeamParticle(1) == BEAM_PARTICLE_P) :
        false;

    static const double jet_antikt_d_04 = 0.4;

    std::vector<fastjet::PseudoJet> particle_truth;
    std::vector<fastjet::PseudoJet> jet_truth_ak04;
    fastjet::ClusterSequenceArea *cluster_sequence_truth = NULL;

    _branch_nmc_truth = 0;

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

    std::fill(_branch_met_truth, _branch_met_truth + 2, 0);

    std::vector<fastjet::PseudoJet> particle_reco_tagged_ak04tpc =
        particle_reco_tpc;
    std::vector<fastjet::PseudoJet> particle_reco_tagged_ak04its =
        particle_reco_its;

    // A value of 2^(-30) < 10^(-9) would map a 10 TeV particle to
    // less than 10 MeV, sufficient to remove any significant momentum
    // bias while not being too greedy to limit the exponent range.
    // Ghost scaling factor is chosen as power of two to maintain the
    // multiplication/division being numerically exact (provided px,
    // py, pz >= 2^(-1022+30) which is of the order 10^(-290) eV).

    static const double scale_ghost = pow(2.0, -30.0);

    if (mc_truth_event != NULL) {
        double met_truth_kahan_error[2] = { 0, 0 };

        for (std::vector<Int_t>::const_iterator iterator =
                 reverse_stored_mc_truth_index.begin();
             iterator != reverse_stored_mc_truth_index.end();
             iterator++) {
            const AliMCParticle *p =
                static_cast<AliMCParticle *>(
                    mc_truth_event->GetTrack(*iterator));

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
                    kahan_sum(_branch_met_truth[0],
                              met_truth_kahan_error[0], p->Px());
                    kahan_sum(_branch_met_truth[1],
                              met_truth_kahan_error[1], p->Py());
                }
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
                std::min(static_cast<Int_t>(UCHAR_MAX),
                         std::max(0,
                                  mc_truth_event->
                                  Particle(*iterator)->
                                  GetStatusCode()));
            // _branch_mc_truth_final_state_primary[_branch_nmc_truth] =
            //     final_state_primary(mc_truth_event, *iterator);
            // _branch_mc_truth_first_parent[_branch_nmc_truth] =
            //     p->Particle()->GetFirstMother();
            // _branch_mc_truth_first_child[_branch_nmc_truth] =
            //     p->Particle()->GetFirstDaughter();
            // _branch_mc_truth_first_child[_branch_nmc_truth] =
            //     p->Particle()->GetLastDaughter();
            _branch_mc_truth_generator_index[_branch_nmc_truth] =
                std::min(static_cast<Short_t>(UCHAR_MAX),
                         std::max(static_cast<Short_t>(0),
                                  p->GetGeneratorIndex()));

            if (p->Particle()->GetFirstMother() >= 0 &&
                p->Particle()->GetFirstMother() <
                mc_truth_event->GetNumberOfTracks()) {
                const AliMCParticle *pp =
                    static_cast<AliMCParticle *>(
                        mc_truth_event->GetTrack(
                            p->Particle()->GetFirstMother()));

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
                    pp->Particle()->GetFirstDaughter() ?
                    pp->Particle()->GetLastDaughter() :
                    *iterator ==
                    pp->Particle()->GetLastDaughter() ?
                    pp->Particle()->GetFirstDaughter() :
                    USHRT_MAX;
            }

            _branch_nmc_truth++;
            if (_branch_nmc_truth >= NMC_TRUTH_MAX) {
                break;
            }
        }

        cluster_sequence_truth = new fastjet::ClusterSequenceArea(
            particle_truth,
            fastjet::JetDefinition(fastjet::JetDefinition(
                fastjet::antikt_algorithm, jet_antikt_d_04)),
            fastjet::VoronoiAreaSpec());
        jet_truth_ak04 = cluster_sequence_truth->inclusive_jets(0);
    }

    // The reco jets will contain EMCAL clusters as ghosts. The idea
    // is calculate a CMS L4 jet energy correction-like
    // electromagnetic fraction (EMF). Its value then can be used as a
    // parameter for the secondary correction, similar to an reversed
    // ATLAS global sequential (GS) correction (the tracking in ALICE
    // being the larger acceptance detector).

    AliVCaloCells *emcal_cell = event->GetEMCALCells();

    for (Int_t i = 0; i < calo_cluster.GetEntriesFast(); i++) {
        AliVCluster *c =
            static_cast<AliVCluster *>(calo_cluster.At(i));

        Int_t cell_id_max = -1;
        Double_t cell_energy_max = -INFINITY;
        Double_t cell_cross = NAN;

        cell_max_cross(cell_id_max, cell_energy_max, cell_cross,
                       c, emcal_cell);
        if (c->GetNCells() > 1 &&
            1 - cell_energy_max / cell_cross < 0.95 &&
            !cell_masked(c, _emcal_mask)) {
            TLorentzVector p;

            c->GetMomentum(p, _branch_primary_vertex);

            const fastjet::PseudoJet
                pj(p.Px(), p.Py(), p.Pz(), p.P());

            particle_reco_tagged_ak04tpc.push_back(pj * scale_ghost);
            particle_reco_tagged_ak04tpc.back().
                set_user_index(USER_INDEX_EM);
            particle_reco_tagged_ak04its.push_back(pj * scale_ghost);
            particle_reco_tagged_ak04its.back().
                set_user_index(USER_INDEX_EM);
        }
    }

    FILL_BRANCH_JET_TRUTH(ak04, jet_truth_ak04);
    TAG_PARTICLE_RECO_JET_TRUTH(particle_reco_tagged_ak04tpc,
                                jet_truth_ak04);
    TAG_PARTICLE_RECO_PARTON(particle_reco_tagged_ak04tpc,
                             algorithmic, ALGORITHMIC);
    TAG_PARTICLE_RECO_JET_TRUTH(particle_reco_tagged_ak04its,
                                jet_truth_ak04);
    TAG_PARTICLE_RECO_PARTON(particle_reco_tagged_ak04its,
                             algorithmic, ALGORITHMIC);

    if (cluster_sequence_truth != NULL) {
        delete cluster_sequence_truth;
    }

    for (size_t i = 0; i < particle_reco_tpc.size(); i++) {
        particle_reco_tpc[i].set_user_index(static_cast<int>(i));
    }

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

    for (size_t i = 0; i < particle_reco_its.size(); i++) {
        particle_reco_its[i].set_user_index(static_cast<int>(i));
    }

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

    static const double jet_kt_d_ue_estimation = 0.3;
    const fastjet::ClusterSequenceArea
        cluster_sequence_ue_estimation_tpc(
            particle_reco_tpc,
            fastjet::JetDefinition(fastjet::JetDefinition(
                fastjet::kt_algorithm, jet_kt_d_ue_estimation)),
            fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_ue_estimation_tpc =
        cluster_sequence_ue_estimation_tpc.inclusive_jets(0);
    const fastjet::ClusterSequenceArea
        cluster_sequence_ue_estimation_its(
            particle_reco_its,
            fastjet::JetDefinition(fastjet::JetDefinition(
                fastjet::kt_algorithm, jet_kt_d_ue_estimation)),
            fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_ue_estimation_its =
        cluster_sequence_ue_estimation_its.inclusive_jets(0);

    // FIXME: Maybe also store ITS UE?

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

    std::pair<std::pair<std::vector<double>, std::vector<double> >,
              double> ue_estimate_tpc =
        ue_estimation_median(cluster_sequence_ue_estimation_tpc,
                             particle_reco_area_tpc);

    std::pair<std::pair<std::vector<double>, std::vector<double> >,
              double> ue_estimate_its =
        ue_estimation_median(cluster_sequence_ue_estimation_its,
                             particle_reco_area_its);

    _branch_ue_estimate_tpc_const =
        evaluate_ue_constant(ue_estimate_tpc.first);
    _branch_ue_estimate_tpc_const_se = ue_estimate_tpc.second;
    _branch_ue_estimate_its_const =
        evaluate_ue_constant(ue_estimate_its.first);
    _branch_ue_estimate_its_const_se = ue_estimate_its.second;

    std::fill(_branch_cell_cluster_index,
              _branch_cell_cluster_index + EMCAL_NCELL, USHRT_MAX);
    _branch_ncluster = 0;

    const Int_t ncalo_cluster = _nrandom_isolation > 0 ?
        _nrandom_isolation : calo_cluster.GetEntriesFast();
    AliESDCaloCluster dummy_cluster;

    for (Int_t i = 0; i < ncalo_cluster; i++) {
#if 0
        if (_nrandom_isolation > 0) {
            TLorentzVector p;

            p.SetPtEtaPhiM(1, _random.Uniform(-1.4, 1.4),
                           _random.Uniform(-M_PI, M_PI), 0);
            for (size_t j = 0; j < 3; j++) {
                dummy_cluster.SetPosition(j, p(j) +
                                          _branch_primary_vertex[j]);
            }
        }
#endif

        AliVCluster *c = _nrandom_isolation > 0 ? &dummy_cluster :
            static_cast<AliVCluster *>(calo_cluster.At(i));
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

        cell_max_cross(cell_id_max, cell_energy_max, energy_cross,
                       c, emcal_cell);
        _branch_cluster_cell_id_max[_branch_ncluster] =
            cell_id_max;
        _branch_cluster_e_max[_branch_ncluster] =
            half(cell_energy_max);
        _branch_cluster_e_cross[_branch_ncluster] =
            half(energy_cross);

        _branch_cluster_nmc_truth[_branch_ncluster] =
            c->GetNLabels();

        // Needed for the isolation below

        std::vector<std::pair<float, unsigned short> >
            mc_truth_energy_index;
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

        if (esd_event != NULL) {
            double cluster_iso_tpc_01 = 0;
            double cluster_iso_tpc_02 = 0;
            double cluster_iso_tpc_03 = 0;
            double cluster_iso_tpc_04 = 0;
            double cluster_iso_tpc_01_ue = 0;
            double cluster_iso_tpc_02_ue = 0;
            double cluster_iso_tpc_03_ue = 0;
            double cluster_iso_tpc_04_ue = 0;

            double cluster_iso_its_01 = 0;
            double cluster_iso_its_02 = 0;
            double cluster_iso_its_03 = 0;
            double cluster_iso_its_04 = 0;
            double cluster_iso_its_01_ue = 0;
            double cluster_iso_its_02_ue = 0;
            double cluster_iso_its_03_ue = 0;
            double cluster_iso_its_04_ue = 0;

            std::vector<std::pair<double, double> > delta_vs_iso_tpc;
            std::vector<std::pair<double, double> > delta_vs_iso_its;
            std::vector<std::pair<double, double> >
                delta_vs_iso_tpc_with_ue;
            std::vector<std::pair<double, double> >
                delta_vs_iso_its_with_ue;

            for (Int_t j = 0; j < esd_event->GetNumberOfTracks(); j++) {
                AliESDtrack *t = esd_event->GetTrack(j);

                if (t == NULL) {
                    continue;
                }

                // Apply PWG-JE cuts (track cuts 0 and 1)

                if (_track_cut[0].AcceptTrack(t) ||
                    _track_cut[1].AcceptTrack(t)) {
                    const double dpseudorapidity = t->Eta() - p.Eta();
                    const double dazimuth = angular_range_reduce(
                        angular_range_reduce(t->Phi()) -
                        angular_range_reduce(p.Phi()));
                    const double dr_2 =
                        std::pow(dpseudorapidity, 2) +
                        std::pow(dazimuth, 2);
                    const double ue =
                        dr_2 < 0.4 * 0.4 ?
                        track_reco_index_tpc.find(j) !=
                        track_reco_index_tpc.end() ?
                        evaluate_ue(ue_estimate_tpc.first, t->Eta(),
                                    t->Phi()) *
                        particle_reco_area_tpc
                        [track_reco_index_tpc[j]] :
                        0 : NAN;
                    const double track_pt_minus_ue =
                        dr_2 < 0.4 * 0.4 ?
                        track_reco_index_tpc.find(j) !=
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
                if (_track_cut[4].AcceptTrack(t)) {
                    const double dpseudorapidity = t->Eta() - p.Eta();
                    const double dazimuth = angular_range_reduce(
                        angular_range_reduce(t->Phi()) -
                        angular_range_reduce(p.Phi()));
                    const double dr_2 =
                        std::pow(dpseudorapidity, 2) +
                        std::pow(dazimuth, 2);
                    const double ue =
                        dr_2 < 0.4 * 0.4 ?
                        track_reco_index_its.find(j) !=
                        track_reco_index_its.end() ?
                        evaluate_ue(ue_estimate_its.first, t->Eta(),
                                    t->Phi()) *
                        particle_reco_area_its
                        [track_reco_index_its[j]] :
                        0 : NAN;
                    const double track_pt_minus_ue =
                        dr_2 < 0.4 * 0.4 ?
                        track_reco_index_its.find(j) !=
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
            }
            _branch_cluster_iso_tpc_01[_branch_ncluster] =
                half(cluster_iso_tpc_01);
            _branch_cluster_iso_tpc_02[_branch_ncluster] =
                half(cluster_iso_tpc_02);
            _branch_cluster_iso_tpc_03[_branch_ncluster] =
                half(cluster_iso_tpc_03);
            _branch_cluster_iso_tpc_04[_branch_ncluster] =
                half(cluster_iso_tpc_04);
            _branch_cluster_iso_tpc_01_ue[_branch_ncluster] =
                half(cluster_iso_tpc_01_ue);
            _branch_cluster_iso_tpc_02_ue[_branch_ncluster] =
                half(cluster_iso_tpc_02_ue);
            _branch_cluster_iso_tpc_03_ue[_branch_ncluster] =
                half(cluster_iso_tpc_03_ue);
            _branch_cluster_iso_tpc_04_ue[_branch_ncluster] =
                half(cluster_iso_tpc_04_ue);

            _branch_cluster_iso_its_01[_branch_ncluster] =
                half(cluster_iso_its_01);
            _branch_cluster_iso_its_02[_branch_ncluster] =
                half(cluster_iso_its_02);
            _branch_cluster_iso_its_03[_branch_ncluster] =
                half(cluster_iso_its_03);
            _branch_cluster_iso_its_04[_branch_ncluster] =
                half(cluster_iso_its_04);
            _branch_cluster_iso_its_01_ue[_branch_ncluster] =
                half(cluster_iso_its_01_ue);
            _branch_cluster_iso_its_02_ue[_branch_ncluster] =
                half(cluster_iso_its_02_ue);
            _branch_cluster_iso_its_03_ue[_branch_ncluster] =
                half(cluster_iso_its_03_ue);
            _branch_cluster_iso_its_04_ue[_branch_ncluster] =
                half(cluster_iso_its_04_ue);

            _branch_cluster_frixione_tpc_04_02[_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                              0.4, 0.2));
            _branch_cluster_frixione_tpc_04_05[_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                              0.4, 0.5));
            _branch_cluster_frixione_tpc_04_10[_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                              0.4, 1.0));

            _branch_cluster_frixione_its_04_02[_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                              0.4, 0.2));
            _branch_cluster_frixione_its_04_05[_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                              0.4, 0.5));
            _branch_cluster_frixione_its_04_10[_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                              0.4, 1.0));

            _branch_cluster_frixione_tpc_04_02_with_ue
                [_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(
                    delta_vs_iso_tpc_with_ue, 0.4, 0.2));
            _branch_cluster_frixione_tpc_04_05_with_ue
                [_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(
                    delta_vs_iso_tpc_with_ue, 0.4, 0.5));
            _branch_cluster_frixione_tpc_04_10_with_ue
                [_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(
                    delta_vs_iso_tpc_with_ue, 0.4, 1.0));

            _branch_cluster_frixione_its_04_02_with_ue
                [_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(
                    delta_vs_iso_its_with_ue, 0.4, 0.2));
            _branch_cluster_frixione_its_04_05_with_ue
                [_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(
                    delta_vs_iso_its_with_ue, 0.4, 0.5));
            _branch_cluster_frixione_its_04_10_with_ue
                [_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(
                    delta_vs_iso_its_with_ue, 0.4, 1.0));

            _branch_cluster_anti_frixione_tpc_04_02
                [_branch_ncluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                                   0.4, 0.2));
            _branch_cluster_anti_frixione_tpc_04_05
                [_branch_ncluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                                   0.4, 0.5));
            _branch_cluster_anti_frixione_tpc_04_10
                [_branch_ncluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_tpc,
                                                   0.4, 1.0));

            _branch_cluster_anti_frixione_its_04_02
                [_branch_ncluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                                   0.4, 0.2));
            _branch_cluster_anti_frixione_its_04_05
                [_branch_ncluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                                   0.4, 0.5));
            _branch_cluster_anti_frixione_its_04_10
                [_branch_ncluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso_its,
                                                   0.4, 1.0));
        }
        else {
            _branch_cluster_iso_tpc_01[_branch_ncluster] = NAN;
            _branch_cluster_iso_tpc_02[_branch_ncluster] = NAN;
            _branch_cluster_iso_tpc_03[_branch_ncluster] = NAN;
            _branch_cluster_iso_tpc_04[_branch_ncluster] = NAN;

            _branch_cluster_iso_its_01[_branch_ncluster] = NAN;
            _branch_cluster_iso_its_02[_branch_ncluster] = NAN;
            _branch_cluster_iso_its_03[_branch_ncluster] = NAN;
            _branch_cluster_iso_its_04[_branch_ncluster] = NAN;

            _branch_cluster_frixione_tpc_04_02[_branch_ncluster] =
                NAN;
            _branch_cluster_frixione_tpc_04_05[_branch_ncluster] =
                NAN;
            _branch_cluster_frixione_tpc_04_10[_branch_ncluster] =
                NAN;

            _branch_cluster_frixione_its_04_02[_branch_ncluster] =
                NAN;
            _branch_cluster_frixione_its_04_05[_branch_ncluster] =
                NAN;
            _branch_cluster_frixione_its_04_10[_branch_ncluster] =
                NAN;

            _branch_cluster_anti_frixione_tpc_04_02
                [_branch_ncluster] = NAN;
            _branch_cluster_anti_frixione_tpc_04_05
                [_branch_ncluster] = NAN;
            _branch_cluster_anti_frixione_tpc_04_10
                [_branch_ncluster] = NAN;

            _branch_cluster_anti_frixione_its_04_02
                [_branch_ncluster] = NAN;
            _branch_cluster_anti_frixione_its_04_05
                [_branch_ncluster] = NAN;
            _branch_cluster_anti_frixione_its_04_10
                [_branch_ncluster] = NAN;
        }

        if (mc_truth_event != NULL) {
            double cluster_iso_01_truth = 0;
            double cluster_iso_02_truth = 0;
            double cluster_iso_03_truth = 0;
            double cluster_iso_04_truth = 0;

            std::vector<std::pair<double, double> > delta_vs_iso;

            for (Int_t j = 0;
                 j < mc_truth_event->GetNumberOfTracks(); j++) {
                if (final_state_primary(mc_truth_event, j) &&
                    cluster_mc_truth_index.find(j) ==
                    cluster_mc_truth_index.end()) {
                    const AliMCParticle *t =
                        static_cast<AliMCParticle *>(
                            mc_truth_event->GetTrack(j));

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
            _branch_cluster_iso_01_truth[_branch_ncluster] =
                half(cluster_iso_01_truth);
            _branch_cluster_iso_02_truth[_branch_ncluster] =
                half(cluster_iso_02_truth);
            _branch_cluster_iso_03_truth[_branch_ncluster] =
                half(cluster_iso_03_truth);
            _branch_cluster_iso_04_truth[_branch_ncluster] =
                half(cluster_iso_04_truth);

            _branch_cluster_frixione_04_02_truth[_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso,
                                              0.4, 0.2));
            _branch_cluster_frixione_04_05_truth[_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso,
                                              0.4, 0.5));
            _branch_cluster_frixione_04_10_truth[_branch_ncluster] =
                half(frixione_iso_max_x_e_eps(delta_vs_iso,
                                              0.4, 1.0));

            _branch_cluster_anti_frixione_04_02_truth
                [_branch_ncluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso,
                                                   0.4, 0.2));
            _branch_cluster_anti_frixione_04_05_truth
                [_branch_ncluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso,
                                                   0.4, 0.5));
            _branch_cluster_anti_frixione_04_10_truth
                [_branch_ncluster] =
                half(anti_frixione_iso_max_x_e_eps(delta_vs_iso,
                                                   0.4, 1.0));
        }
        else {
            _branch_cluster_iso_01_truth[_branch_ncluster] = NAN;
            _branch_cluster_iso_02_truth[_branch_ncluster] = NAN;
            _branch_cluster_iso_03_truth[_branch_ncluster] = NAN;
            _branch_cluster_iso_04_truth[_branch_ncluster] = NAN;

            _branch_cluster_frixione_04_02_truth[_branch_ncluster] =
                NAN;
            _branch_cluster_frixione_04_05_truth[_branch_ncluster] =
                NAN;
            _branch_cluster_frixione_04_10_truth[_branch_ncluster] =
                NAN;

            _branch_cluster_anti_frixione_04_02_truth
                [_branch_ncluster] = NAN;
            _branch_cluster_anti_frixione_04_05_truth
                [_branch_ncluster] = NAN;
            _branch_cluster_anti_frixione_04_10_truth
                [_branch_ncluster] = NAN;
        }

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

        _branch_ncluster++;
        if (_branch_ncluster >= NCLUSTER_MAX) {
            break;
        }
    }

    calo_cluster.Delete();

    FILL_BRANCH_JET(ak04tpc, jet_reco_ak04tpc,
                    cluster_sequence_reco_ak04tpc,
                    jet_reco_tagged_ak04tpc,
                    cluster_sequence_reco_tagged_ak04tpc,
                    ak04, jet_truth_ak04,
                    particle_reco_area_tpc, ue_estimate_tpc);
    FILL_BRANCH_JET(ak04its, jet_reco_ak04its,
                    cluster_sequence_reco_ak04its,
                    jet_reco_tagged_ak04its,
                    cluster_sequence_reco_tagged_ak04its,
                    ak04, jet_truth_ak04,
                    particle_reco_area_its, ue_estimate_its);

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
        if (pt.size() >= 3 &&
            !(pt[0] >= _skim_jet_min_pt[0] &&
              pt[1] >= _skim_jet_min_pt[1] &&
              pt[2] >= _skim_jet_min_pt[2])) {
            // Discard this event
            return;
        }
        else if (pt.size() >= 2 &&
                 !(_skim_jet_min_pt[2] > -INFINITY) &&
                 !(pt[0] >= _skim_jet_min_pt[0] &&
                   pt[1] >= _skim_jet_min_pt[1])) {
            // Discard this event
            return;
        }
        else if (pt.size() >= 1 &&
                 !(_skim_jet_min_pt[1] > -INFINITY ||
                   _skim_jet_min_pt[2] > -INFINITY) &&
                 !(pt[0] >= _skim_jet_min_pt[0])) {
            // Discard this event
            return;
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
        if (pt.size() >= 2) {
            fprintf(stdout, "%s:%d: %f %f\n", __FILE__, __LINE__, 0.5 * (pt[0] + pt[1]), _skim_jet_average_pt);
        }
        if (!(pt.size() >= 2 &&
              0.5 * (pt[0] + pt[1]) >= _skim_jet_average_pt)) {
            // Discard this event
            return;
        }
    }

    std::fill(_branch_cell_e, _branch_cell_e + EMCAL_NCELL, NAN);
    std::fill(_branch_cell_tof, _branch_cell_tof + EMCAL_NCELL, NAN);
    std::fill(_branch_cell_mc_truth_index,
              _branch_cell_mc_truth_index + EMCAL_NCELL, USHRT_MAX);
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

#undef SAFE_MC_TRUTH_INDEX_TO_USHRT

        }
    }

    _branch_nmuon_track = 0;
    if (esd_event != NULL) {
        for (Int_t i = 0; i < esd_event->GetNumberOfMuonTracks(); i++) {
            AliESDMuonTrack *t = esd_event->GetMuonTrack(i);

            if (t == NULL) {
                continue;
            }

            _branch_muon_track_e[_branch_nmuon_track] =
                half(t->E());
            _branch_muon_track_pt[_branch_nmuon_track] =
                half(t->Pt());
            _branch_muon_track_eta[_branch_nmuon_track] =
                half(t->Eta());
            _branch_muon_track_phi[_branch_nmuon_track] =
                half(angular_range_reduce(t->Phi()));

            _branch_muon_track_r_abs[_branch_nmuon_track] =
                half(t->GetRAtAbsorberEnd());
            _branch_muon_track_p_dca[_branch_nmuon_track] = NAN;

            static const double c_505_tan_3_pi_180 =
                505 * tan(3 * M_PI / 180);

            // Default values from
            // AliPhysics/PWG/muon/buildMuonTrackCutsOADB.C

            static const double default_sigma_p_dca_23 = 80;
            static const double default_sigma_p_dca_310 = 54;

            _branch_muon_track_sigma_p_dca[_branch_nmuon_track] =
                t->GetRAtAbsorberEnd() < c_505_tan_3_pi_180 ?
                default_sigma_p_dca_23 : default_sigma_p_dca_310;
            if (_muon_track_cut != NULL) {
                _branch_muon_track_p_dca[_branch_nmuon_track] =
                    half(_muon_track_cut->GetAverageMomentum(t) *
                         _muon_track_cut->GetCorrectedDCA(t).Mag());

                // Here, muon_track_sigma_p_dca is the resoultion
                // corrected value, see
                // AliMuonTrackCuts::GetSelectionMask() in
                // AliPhysics/PWG/muon/AliMuonTrackCuts.cxx

                const AliOADBMuonTrackCutsParam param =
                    _muon_track_cut->GetMuonTrackCutsParam();
                const Double_t sigma_p_dca =
                    _muon_track_cut->IsThetaAbs23(t) ?
                    param.GetSigmaPdca23() : param.GetSigmaPdca310();

                _branch_muon_track_sigma_p_dca[_branch_nmuon_track] =
                    half(sigma_p_dca);

                // In AliMuonTrackCuts::GetSelectionMask(), it would
                // be nrp = nsigma * p * dp

                const Double_t delta_sagitta_p =
                    param.GetRelPResolution() * t->P();

                _branch_muon_track_delta_sagitta_p
                    [_branch_nmuon_track] = half(delta_sagitta_p);

                // p_resolution_effect = sigma_p_dca / (1 - nrp / (1 +
                // nrp));

                static const Double_t z_tc12_cm = 535;
                const Double_t distance_sigma_slope_p_meas =
                    z_tc12_cm * param.GetSlopeResolution() * t->P();

                _branch_muon_track_distance_sigma_slope_p
                    [_branch_nmuon_track] =
                    half(distance_sigma_slope_p_meas);
            }
            _branch_nmuon_track++;
            if (_branch_nmuon_track >= NTRACK_MAX) {
                break;
            }
        }
    }
    else if (aod_event != NULL) {
        // FIXME: Not really implemented
        TRefArray muon_track;

        aod_event->GetMuonTracks(&muon_track);
    }

    // Now that the event is accepted, copy over the total counted
    // ntrials so far, and reset the ntrial counter
    _branch_eg_ntrial = _skim_sum_eg_ntrial;
    _skim_sum_eg_ntrial = 0;

    _tree_event->Fill();
}

AliEMCALRecoUtils *AliAnalysisTaskNtuple::GetEMCALRecoUtils(void)
{
    return _reco_util;
}

void AliAnalysisTaskNtuple::SetAliROOTVersion(const char *version)
{
    strncpy(_branch_version_aliroot, version, BUFSIZ);
    _branch_version_aliroot[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNtuple::SetAliPhysicsVersion(const char *version)
{
    strncpy(_branch_version_aliphysics, version, BUFSIZ);
    _branch_version_aliphysics[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNtuple::SetGridDataDir(const char *dir)
{
    strncpy(_branch_grid_data_dir, dir, BUFSIZ);
    _branch_grid_data_dir[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNtuple::SetGridDataPattern(const char *pattern)
{
    strncpy(_branch_grid_data_pattern, pattern, BUFSIZ);
    _branch_grid_data_pattern[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNtuple::
SetEMCALGeometryFilename(const char *emcal_geometry_filename)
{
    _emcal_geometry_filename = emcal_geometry_filename;
}

void AliAnalysisTaskNtuple::
SetEMCALLocal2MasterFilename(const char *emcal_local2master_filename)
{
    _emcal_local2master_filename = emcal_local2master_filename;
}

void AliAnalysisTaskNtuple::
SetForceUESubtraction(bool force_ue_subtraction)
{
    _force_ue_subtraction = force_ue_subtraction;
}

void AliAnalysisTaskNtuple::SetSkimClusterMinE(double min_e)
{
    _skim_cluster_min_e = min_e;
}

void AliAnalysisTaskNtuple::SetSkimTrackMinPt(double min_pt)
{
    _skim_track_min_pt = min_pt;
}

void AliAnalysisTaskNtuple::SetSkimMuonTrackMinPt(double min_pt)
{
    _skim_muon_track_min_pt = min_pt;
}

void AliAnalysisTaskNtuple::SetSkimJetMinPt(double min_pt_1,
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

void AliAnalysisTaskNtuple::SetSkimJetAveragePt(double average_pt)
{
    _skim_jet_average_pt = average_pt;
}

void AliAnalysisTaskNtuple::SetSkimMultiplicityTrackletMinN(int min_n)
{
    _skim_multiplicity_tracklet_min_n = min_n;
}

void AliAnalysisTaskNtuple::SetStoredTrackMinPt(double min_pt)
{
    _stored_track_min_pt = min_pt;
}

// This is primarily intended for JEC derivation, and therefore
// defined in pT raw

void AliAnalysisTaskNtuple::SetStoredJetMinPtRaw(double min_pt_raw)
{
    _stored_jet_min_pt_raw = min_pt_raw;
}

void AliAnalysisTaskNtuple::
SetNRandomIsolation(unsigned int nrandom_isolation)
{
    _nrandom_isolation = nrandom_isolation;
}
