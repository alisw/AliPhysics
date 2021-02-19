
// -*- mode: c++; -*-

#ifndef ALIANALYSISTASKPHOTONDISC_H_
#define ALIANALYSISTASKPHOTONDISC_H_

#include <limits.h>
#include <vector>
#include <set>
#include <TList.h>
#include <TH2D.h>
#include <TFormula.h>
#include <TClonesArray.h>
#include <AliESDtrackCuts.h>
#include <AliEMCALGeometry.h>
#include <AliEMCALRecoUtils.h>
#include <AliMuonTrackCuts.h>
#include <AliAnalysisTaskEmcal.h>
#include <AliAnalysisAlien.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wdeprecated-register"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequenceArea.hh>

#define EMCAL_NCELL         17664
#define NTRIGGER_CLASS_MAX  100 // = AliESDRun::kNTriggerClasses

#define NCLUSTER_MAX        (1U << 17)
#define NTRACK_MAX          (1U << 17)
#define NMC_TRUTH_MAX       (1U << 17)
#define NJET_MAX            (1U << 17)

// #define WITH_EFP7 1

#ifdef WITH_EFP7
#define IF_EFP7(...) __VA_ARGS__
#else // WITH_EFP7
#define IF_EFP7(...)
#endif // WITH_EFP7

#define CLUSTER_NMC_TRUTH_MAX 32

class AliAnalysisTaskNTGJ : public AliAnalysisTaskEmcal {
private:
    TString _emcal_geometry_name; //!
    TTree *_tree_event; //!

#define MEMBER_BRANCH                                               \
    /* */                                                           \
    BRANCH_STR(id_git)                                              \
    BRANCH_STR(version_aliroot)                                     \
    BRANCH_STR(version_aliphysics)                                  \
    BRANCH_STR(version_jec)                                         \
    BRANCH_STR(grid_data_dir)                                       \
    BRANCH_STR(grid_data_pattern)                                   \
    BRANCH(beam_energy, F)                                          \
    BRANCH_ARRAY(beam_particle, 2, I)                               \
    BRANCH(ntrigger_class, b)                                       \
    BRANCH_STR_ARRAY(trigger_class, ntrigger_class)                 \
    BRANCH(run_number, I)                                           \
    BRANCH(period_number, i)                                        \
    BRANCH(orbit_number, i)                                         \
    BRANCH(bunch_crossing_number, s)                                \
    BRANCH(time_stamp, i)                                           \
    /* */                                                           \
    BRANCH_ARRAY(trigger_mask, 2, l)                                \
    BRANCH_ARRAY(multiplicity_v0, 64, F)                            \
    BRANCH(centrality_v0m, F)                                       \
    BRANCH_ARRAY(centrality, 9, F)                                  \
    BRANCH_ARRAY(event_plane_psi_v0, 3, F)                          \
    BRANCH_ARRAY2(event_plane_q_v0, 3, 2, D)                        \
    BRANCH(has_misalignment_matrix, O)                              \
    BRANCH_ARRAY(primary_vertex, 3, D)                              \
    BRANCH_ARRAY(primary_vertex_sigma, 3, D)                        \
    BRANCH(primary_vertex_ncontributor, I)                          \
    BRANCH_ARRAY(primary_vertex_spd, 3, D)                          \
    BRANCH_ARRAY(primary_vertex_spd_sigma, 3, D)                    \
    BRANCH(primary_vertex_spd_ncontributor, I)                      \
    BRANCH(is_pileup_from_spd_3_08, O)                              \
    BRANCH(is_pileup_from_spd_5_08, O)                              \
    BRANCH(npileup_vertex_spd, I)                                   \
    BRANCH(ncluster_tpc, I)                                         \
    BRANCH(event_selected, O)                                       \
    BRANCH(eg_signal_process_id, I)                                 \
    BRANCH(eg_mpi, I)                                               \
    BRANCH(eg_pt_hat, F)                                            \
    BRANCH(eg_cross_section, F)                                     \
    BRANCH(eg_weight, F)                                            \
    BRANCH_ARRAY(eg_primary_vertex, 3, F)                           \
    BRANCH(eg_ntrial, I)                                            \
    BRANCH(eg_scale_pdf, F)                                         \
    BRANCH(eg_alpha_qcd, F)                                         \
    BRANCH(eg_alpha_qed, F)                                         \
    BRANCH_ARRAY(eg_pdf_id, 2, I)                                   \
    BRANCH_ARRAY(eg_pdf_x, 2, F)                                    \
    BRANCH_ARRAY(eg_pdf_x_pdf, 2, F)                                \
    /* */                                                           \
    BRANCH(ncluster, i)                                             \
    BRANCH_ARRAY(cluster_e, ncluster, F)                            \
    BRANCH_ARRAY(cluster_pt, ncluster, F)                           \
    BRANCH_ARRAY(cluster_eta, ncluster, F)                          \
    BRANCH_ARRAY(cluster_phi, ncluster, F)                          \
    BRANCH_ARRAY2(cluster_lambda_square, ncluster, 2, F)            \
    BRANCH_ARRAY2(cluster_lambda_square_angle, ncluster, 2, F)      \
    BRANCH_ARRAY(cluster_tof, ncluster, F)                          \
    BRANCH_ARRAY(cluster_ncell, ncluster, I)                        \
    BRANCH_ARRAY(cluster_nlocal_maxima, ncluster, b)                \
    BRANCH_ARRAY(cluster_distance_to_bad_channel, ncluster, F)      \
    BRANCH_ARRAY(cluster_cell_id_max, ncluster, s)                  \
    BRANCH_ARRAY(cluster_e_max, ncluster, F)                        \
    BRANCH_ARRAY(cluster_e_cross, ncluster, F)                      \
    BRANCH_ARRAY(cluster_nmc_truth, ncluster, i)                    \
    BRANCH_ARRAY2(cluster_mc_truth_index, ncluster, 32, s)          \
    BRANCH_ARRAY(cluster_iso_tpc_01, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_tpc_01_ue, ncluster, F)                \
    BRANCH_ARRAY(cluster_iso_tpc_02, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_tpc_02_ue, ncluster, F)                \
    BRANCH_ARRAY(cluster_iso_tpc_03, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_tpc_03_ue, ncluster, F)                \
    BRANCH_ARRAY(cluster_iso_tpc_04, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_tpc_04_ue, ncluster, F)                \
    BRANCH_ARRAY(cluster_iso_its_01, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_its_01_ue, ncluster, F)                \
    BRANCH_ARRAY(cluster_iso_its_02, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_its_02_ue, ncluster, F)                \
    BRANCH_ARRAY(cluster_iso_its_03, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_its_03_ue, ncluster, F)                \
    BRANCH_ARRAY(cluster_iso_its_04, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_its_04_ue, ncluster, F)                \
    BRANCH_ARRAY(cluster_iso_cluster_01, ncluster, F)               \
    BRANCH_ARRAY(cluster_iso_cluster_01_ue, ncluster, F)            \
    BRANCH_ARRAY(cluster_iso_cluster_02, ncluster, F)               \
    BRANCH_ARRAY(cluster_iso_cluster_02_ue, ncluster, F)            \
    BRANCH_ARRAY(cluster_iso_cluster_03, ncluster, F)               \
    BRANCH_ARRAY(cluster_iso_cluster_03_ue, ncluster, F)            \
    BRANCH_ARRAY(cluster_iso_cluster_04, ncluster, F)               \
    BRANCH_ARRAY(cluster_iso_cluster_04_ue, ncluster, F)            \
    BRANCH_ARRAY(cluster_frixione_tpc_04_02, ncluster, F)           \
    BRANCH_ARRAY(cluster_frixione_tpc_04_05, ncluster, F)           \
    BRANCH_ARRAY(cluster_frixione_tpc_04_10, ncluster, F)           \
    BRANCH_ARRAY(cluster_frixione_its_04_02, ncluster, F)           \
    BRANCH_ARRAY(cluster_frixione_its_04_05, ncluster, F)           \
    BRANCH_ARRAY(cluster_frixione_its_04_10, ncluster, F)           \
    BRANCH_ARRAY(cluster_frixione_tpc_04_02_with_ue, ncluster, F)   \
    BRANCH_ARRAY(cluster_frixione_tpc_04_05_with_ue, ncluster, F)   \
    BRANCH_ARRAY(cluster_frixione_tpc_04_10_with_ue, ncluster, F)   \
    BRANCH_ARRAY(cluster_frixione_its_04_02_with_ue, ncluster, F)   \
    BRANCH_ARRAY(cluster_frixione_its_04_05_with_ue, ncluster, F)   \
    BRANCH_ARRAY(cluster_frixione_its_04_10_with_ue, ncluster, F)   \
    BRANCH_ARRAY(cluster_anti_frixione_tpc_04_02, ncluster, F)      \
    BRANCH_ARRAY(cluster_anti_frixione_tpc_04_05, ncluster, F)      \
    BRANCH_ARRAY(cluster_anti_frixione_tpc_04_10, ncluster, F)      \
    BRANCH_ARRAY(cluster_anti_frixione_its_04_02, ncluster, F)      \
    BRANCH_ARRAY(cluster_anti_frixione_its_04_05, ncluster, F)      \
    BRANCH_ARRAY(cluster_anti_frixione_its_04_10, ncluster, F)      \
    BRANCH_ARRAY(cluster_iso_01_truth, ncluster, F)                 \
    BRANCH_ARRAY(cluster_iso_02_truth, ncluster, F)                 \
    BRANCH_ARRAY(cluster_iso_03_truth, ncluster, F)                 \
    BRANCH_ARRAY(cluster_iso_04_truth, ncluster, F)                 \
    BRANCH_ARRAY(cluster_frixione_04_02_truth, ncluster, F)         \
    BRANCH_ARRAY(cluster_frixione_04_05_truth, ncluster, F)         \
    BRANCH_ARRAY(cluster_frixione_04_10_truth, ncluster, F)         \
    BRANCH_ARRAY(cluster_anti_frixione_04_02_truth, ncluster, F)    \
    BRANCH_ARRAY(cluster_anti_frixione_04_05_truth, ncluster, F)    \
    BRANCH_ARRAY(cluster_anti_frixione_04_10_truth, ncluster, F)    \
    BRANCH_ARRAY2(cluster_s_nphoton, ncluster, 4, F)                \
    BRANCH_ARRAY2(cluster_s_ncharged_hadron, ncluster, 4, F)        \
    BRANCH_ARRAY(cell_e, 17664, F)                                  \
    BRANCH_ARRAY2(cell_position, 17664, 3, F)                       \
    BRANCH_ARRAY(cell_tof, 17664, F)                                \
    BRANCH_ARRAY(cell_cluster_index, 17664, s)                      \
    BRANCH_ARRAY(cell_mc_truth_index, 17664, s)                     \
    BRANCH_ARRAY(cell_voronoi_area, 17664, F)                       \
    /* */                                                           \
    BRANCH(ntrack, i)                                               \
    BRANCH_ARRAY(track_e, ntrack, F)                                \
    BRANCH_ARRAY(track_pt, ntrack, F)                               \
    BRANCH_ARRAY(track_eta, ntrack, F)                              \
    BRANCH_ARRAY(track_phi, ntrack, F)                              \
    BRANCH_ARRAY(track_eta_emcal, ntrack, F)                        \
    BRANCH_ARRAY(track_phi_emcal, ntrack, F)                        \
    BRANCH_ARRAY(track_charge, ntrack, B)                           \
    BRANCH_ARRAY(track_quality, ntrack, b)                          \
    BRANCH_ARRAY(track_tpc_dedx, ntrack, F)                         \
    BRANCH_ARRAY(track_tpc_length_active_zone, ntrack, F)           \
    BRANCH_ARRAY(track_tpc_xrow, ntrack, b)                         \
    BRANCH_ARRAY(track_tpc_ncluster, ntrack, b)                     \
    BRANCH_ARRAY(track_tpc_ncluster_dedx, ntrack, b)                \
    BRANCH_ARRAY(track_tpc_ncluster_findable, ntrack, b)            \
    BRANCH_ARRAY(track_its_ncluster, ntrack, b)                     \
    BRANCH_ARRAY(track_its_chi_square, ntrack, F)                   \
    BRANCH_ARRAY(track_dca_xy, ntrack, F)                           \
    BRANCH_ARRAY(track_dca_z, ntrack, F)                            \
    BRANCH_ARRAY(track_mc_truth_index, ntrack, s)                   \
    BRANCH_ARRAY(track_voronoi_area, ntrack, F)                     \
    /* */                                                           \
    BRANCH(nmuon_track, i)                                          \
    BRANCH_ARRAY(muon_track_e, nmuon_track, F)                      \
    BRANCH_ARRAY(muon_track_pt, nmuon_track, F)                     \
    BRANCH_ARRAY(muon_track_eta, nmuon_track, F)                    \
    BRANCH_ARRAY(muon_track_phi, nmuon_track, F)                    \
    BRANCH_ARRAY(muon_track_r_abs, nmuon_track, F)                  \
    BRANCH_ARRAY(muon_track_p_dca, nmuon_track, F)                  \
    BRANCH_ARRAY(muon_track_sigma_p_dca, nmuon_track, F)            \
    BRANCH_ARRAY(muon_track_delta_sagitta_p, nmuon_track, F)        \
    BRANCH_ARRAY(muon_track_distance_sigma_slope_p, nmuon_track, F) \
    BRANCH_ARRAY(muon_track_mc_truth_index, nmuon_track, s)         \
    /* */                                                           \
    BRANCH(nmc_truth, i)                                            \
    BRANCH_ARRAY(mc_truth_e, nmc_truth, F)                          \
    BRANCH_ARRAY(mc_truth_pt, nmc_truth, F)                         \
    BRANCH_ARRAY(mc_truth_eta, nmc_truth, F)                        \
    BRANCH_ARRAY(mc_truth_phi, nmc_truth, F)                        \
    BRANCH_ARRAY(mc_truth_charge, nmc_truth, B)                     \
    BRANCH_ARRAY(mc_truth_pdg_code, nmc_truth, S)                   \
    BRANCH_ARRAY(mc_truth_status, nmc_truth, b)                     \
    BRANCH_ARRAY(mc_truth_generator_index, nmc_truth, b)            \
    /* BRANCH_ARRAY(mc_truth_physical_primary, nmc_truth, O) */     \
    /* BRANCH_ARRAY(mc_truth_first_parent, nmc_truth, I) */         \
    /* BRANCH_ARRAY(mc_truth_first_child, nmc_truth, I) */          \
    /* BRANCH_ARRAY(mc_truth_second_child, nmc_truth, I) */         \
    BRANCH_ARRAY(mc_truth_first_parent_pdg_code, nmc_truth, S)      \
    BRANCH_ARRAY(mc_truth_first_parent_e, nmc_truth, F)             \
    BRANCH_ARRAY(mc_truth_first_parent_pt, nmc_truth, F)            \
    BRANCH_ARRAY(mc_truth_first_parent_eta, nmc_truth, F)           \
    BRANCH_ARRAY(mc_truth_first_parent_phi, nmc_truth, F)           \
    BRANCH_ARRAY(mc_truth_sibling_index, nmc_truth, s)              \
    BRANCH(debug_njet_ue_estimation, i)                             \
    BRANCH_ARRAY(debug_jet_ue_estimation_pt_raw,                    \
                 debug_njet_ue_estimation, F)                       \
    BRANCH_ARRAY(debug_jet_ue_estimation_eta_raw,                   \
                 debug_njet_ue_estimation, F)                       \
    BRANCH_ARRAY(debug_jet_ue_estimation_phi_raw,                   \
                 debug_njet_ue_estimation, F)                       \
    BRANCH_ARRAY(debug_jet_ue_estimation_area_raw,                  \
                 debug_njet_ue_estimation, F)                       \
    /* */                                                           \
    BRANCH(ue_estimate_tpc_const, F)                                \
    BRANCH(ue_estimate_tpc_const_se, F)                             \
    BRANCH(ue_estimate_its_const, F)                                \
    BRANCH(ue_estimate_its_const_se, F)                             \
    BRANCH(ue_estimate_cluster_const, F)                            \
    BRANCH(ue_estimate_cluster_const_se, F)                         \
    /* */                                                           \
    BRANCH(njet_ak04tpc, i)                                         \
    BRANCH_ARRAY(debug_jet_ak04tpc_tag_dr_square, njet_ak04tpc, F)  \
    BRANCH_ARRAY(jet_ak04tpc_e_raw, njet_ak04tpc, F)                \
    BRANCH_ARRAY(jet_ak04tpc_e, njet_ak04tpc, F)                    \
    BRANCH_ARRAY(jet_ak04tpc_e_charged, njet_ak04tpc, F)            \
    BRANCH_ARRAY(jet_ak04tpc_pt_raw_ue, njet_ak04tpc, F)            \
    BRANCH_ARRAY(jet_ak04tpc_pt_raw, njet_ak04tpc, F)               \
    BRANCH_ARRAY(jet_ak04tpc_pt, njet_ak04tpc, F)                   \
    BRANCH_ARRAY(jet_ak04tpc_pt_charged, njet_ak04tpc, F)           \
    BRANCH_ARRAY(jet_ak04tpc_eta_raw, njet_ak04tpc, F)              \
    BRANCH_ARRAY(jet_ak04tpc_eta, njet_ak04tpc, F)                  \
    BRANCH_ARRAY(jet_ak04tpc_phi, njet_ak04tpc, F)                  \
    BRANCH_ARRAY(jet_ak04tpc_m2_raw, njet_ak04tpc, F)               \
    BRANCH_ARRAY(jet_ak04tpc_m2, njet_ak04tpc, F)                   \
    BRANCH_ARRAY(jet_ak04tpc_area_raw, njet_ak04tpc, F)             \
    BRANCH_ARRAY(jet_ak04tpc_area, njet_ak04tpc, F)                 \
    BRANCH_ARRAY(jet_ak04tpc_emf_raw, njet_ak04tpc, F)              \
    BRANCH_ARRAY(jet_ak04tpc_emf, njet_ak04tpc, F)                  \
    BRANCH_ARRAY(jet_ak04tpc_multiplicity_raw, njet_ak04tpc, s)     \
    BRANCH_ARRAY(jet_ak04tpc_multiplicity, njet_ak04tpc, F)         \
    BRANCH_ARRAY2(jet_ak04tpc_width_sigma_raw, njet_ak04tpc, 2, F)  \
    BRANCH_ARRAY2(jet_ak04tpc_width_sigma, njet_ak04tpc, 2, F)      \
    BRANCH_ARRAY(jet_ak04tpc_ptd_raw, njet_ak04tpc, F)              \
    BRANCH_ARRAY(jet_ak04tpc_ptd, njet_ak04tpc, F)                  \
    IF_EFP7(BRANCH_ARRAY2(jet_ak04tpc_efp_raw,                      \
                          njet_ak04tpc, 489, F))                    \
    BRANCH_STR(debug_blas_version)                                  \
    BRANCH_ARRAY2(jet_ak04tpc_truth_index_z_truth, njet_ak04tpc,    \
                  2, I)                                             \
    BRANCH_ARRAY2(jet_ak04tpc_truth_z_truth, njet_ak04tpc, 2, F)    \
    BRANCH_ARRAY2(jet_ak04tpc_truth_index_z_reco, njet_ak04tpc,     \
                  2, I)                                             \
    BRANCH_ARRAY2(jet_ak04tpc_truth_z_reco, njet_ak04tpc, 2, F)     \
    BRANCH_ARRAY2(jet_ak04tpc_charged_truth_index_z_truth,          \
                  njet_ak04tpc, 2, I)                               \
    BRANCH_ARRAY2(jet_ak04tpc_charged_truth_z_truth, njet_ak04tpc,  \
                  2, F)                                             \
    BRANCH_ARRAY2(jet_ak04tpc_charged_truth_index_z_reco,           \
                  njet_ak04tpc, 2, I)                               \
    BRANCH_ARRAY2(jet_ak04tpc_charged_truth_z_reco, njet_ak04tpc,   \
                  2, F)                                             \
    BRANCH_ARRAY2(jet_ak04tpc_pdg_code_algorithmic, njet_ak04tpc,   \
                  2, I)                                             \
    BRANCH_ARRAY2(jet_ak04tpc_pdg_code_algorithmic_z, njet_ak04tpc, \
                  2, F)                                             \
    BRANCH_ARRAY(jet_ak04tpc_e_truth, njet_ak04tpc, F)              \
    BRANCH_ARRAY(jet_ak04tpc_pt_truth, njet_ak04tpc, F)             \
    BRANCH_ARRAY(jet_ak04tpc_eta_truth, njet_ak04tpc, F)            \
    BRANCH_ARRAY(jet_ak04tpc_phi_truth, njet_ak04tpc, F)            \
    BRANCH_ARRAY(jet_ak04tpc_e_charged_truth, njet_ak04tpc, F)      \
    BRANCH_ARRAY(jet_ak04tpc_pt_charged_truth, njet_ak04tpc, F)     \
    BRANCH_ARRAY(jet_ak04tpc_eta_charged_truth, njet_ak04tpc, F)    \
    BRANCH_ARRAY(jet_ak04tpc_phi_charged_truth, njet_ak04tpc, F)    \
    /* */                                                           \
    BRANCH(njet_ak04its, i)                                         \
    BRANCH_ARRAY(debug_jet_ak04its_tag_dr_square, njet_ak04its, F)  \
    BRANCH_ARRAY(jet_ak04its_e_raw, njet_ak04its, F)                \
    BRANCH_ARRAY(jet_ak04its_e, njet_ak04its, F)                    \
    BRANCH_ARRAY(jet_ak04its_e_charged, njet_ak04its, F)            \
    BRANCH_ARRAY(jet_ak04its_pt_raw_ue, njet_ak04its, F)            \
    BRANCH_ARRAY(jet_ak04its_pt_raw, njet_ak04its, F)               \
    BRANCH_ARRAY(jet_ak04its_pt, njet_ak04its, F)                   \
    BRANCH_ARRAY(jet_ak04its_pt_charged, njet_ak04its, F)           \
    BRANCH_ARRAY(jet_ak04its_eta_raw, njet_ak04its, F)              \
    BRANCH_ARRAY(jet_ak04its_eta, njet_ak04its, F)                  \
    BRANCH_ARRAY(jet_ak04its_phi, njet_ak04its, F)                  \
    BRANCH_ARRAY(jet_ak04its_m2_raw, njet_ak04tpc, F)               \
    BRANCH_ARRAY(jet_ak04its_m2, njet_ak04tpc, F)                   \
    BRANCH_ARRAY(jet_ak04its_area_raw, njet_ak04its, F)             \
    BRANCH_ARRAY(jet_ak04its_area, njet_ak04its, F)                 \
    BRANCH_ARRAY(jet_ak04its_emf_raw, njet_ak04its, F)              \
    BRANCH_ARRAY(jet_ak04its_emf, njet_ak04its, F)                  \
    BRANCH_ARRAY(jet_ak04its_multiplicity_raw, njet_ak04its, s)     \
    BRANCH_ARRAY(jet_ak04its_multiplicity, njet_ak04its, F)         \
    BRANCH_ARRAY2(jet_ak04its_width_sigma_raw, njet_ak04its, 2, F)  \
    BRANCH_ARRAY2(jet_ak04its_width_sigma, njet_ak04its, 2, F)      \
    BRANCH_ARRAY(jet_ak04its_ptd_raw, njet_ak04its, F)              \
    BRANCH_ARRAY(jet_ak04its_ptd, njet_ak04its, F)                  \
    BRANCH_ARRAY2(jet_ak04its_truth_index_z_truth, njet_ak04its,    \
                  2, I)                                             \
    BRANCH_ARRAY2(jet_ak04its_truth_z_truth, njet_ak04its, 2, F)    \
    BRANCH_ARRAY2(jet_ak04its_truth_index_z_reco, njet_ak04its,     \
                  2, I)                                             \
    BRANCH_ARRAY2(jet_ak04its_truth_z_reco, njet_ak04its, 2, F)     \
    BRANCH_ARRAY2(jet_ak04its_charged_truth_index_z_truth,          \
                  njet_ak04its, 2, I)                               \
    BRANCH_ARRAY2(jet_ak04its_charged_truth_z_truth, njet_ak04its,  \
                  2, F)                                             \
    BRANCH_ARRAY2(jet_ak04its_charged_truth_index_z_reco,           \
                  njet_ak04its, 2, I)                               \
    BRANCH_ARRAY2(jet_ak04its_charged_truth_z_reco, njet_ak04its,   \
                  2, F)                                             \
    BRANCH_ARRAY2(jet_ak04its_pdg_code_algorithmic, njet_ak04its,   \
                  2, I)                                             \
    BRANCH_ARRAY2(jet_ak04its_pdg_code_algorithmic_z, njet_ak04its, \
                  2, F)                                             \
    BRANCH_ARRAY(jet_ak04its_e_truth, njet_ak04its, F)              \
    BRANCH_ARRAY(jet_ak04its_pt_truth, njet_ak04its, F)             \
    BRANCH_ARRAY(jet_ak04its_eta_truth, njet_ak04its, F)            \
    BRANCH_ARRAY(jet_ak04its_phi_truth, njet_ak04its, F)            \
    BRANCH_ARRAY(jet_ak04its_e_charged_truth, njet_ak04its, F)      \
    BRANCH_ARRAY(jet_ak04its_pt_charged_truth, njet_ak04its, F)     \
    BRANCH_ARRAY(jet_ak04its_eta_charged_truth, njet_ak04its, F)    \
    BRANCH_ARRAY(jet_ak04its_phi_charged_truth, njet_ak04its, F)    \
    /* */                                                           \
    BRANCH(njet_truth_ak04, i)                                      \
    BRANCH_ARRAY(jet_truth_ak04_e, njet_truth_ak04, F)              \
    BRANCH_ARRAY(jet_truth_ak04_pt, njet_truth_ak04, F)             \
    BRANCH_ARRAY(jet_truth_ak04_eta, njet_truth_ak04, F)            \
    BRANCH_ARRAY(jet_truth_ak04_phi, njet_truth_ak04, F)            \
    BRANCH_ARRAY(jet_truth_ak04_area, njet_truth_ak04, F)           \
    BRANCH_ARRAY(jet_truth_ak04_emf, njet_truth_ak04, F)            \
    BRANCH_ARRAY(jet_truth_ak04_multiplicity, njet_truth_ak04, s)   \
    BRANCH_ARRAY2(jet_truth_ak04_width_sigma, njet_truth_ak04,      \
                  2, F)                                             \
    BRANCH_ARRAY(jet_truth_ak04_ptd, njet_truth_ak04, F)            \
    IF_EFP7(BRANCH_ARRAY2(jet_truth_ak04_efp,                       \
                          njet_truth_ak04, 489, F))                 \
    /* */                                                           \
    BRANCH(njet_charged_truth_ak04, i)                              \
    BRANCH_ARRAY(jet_charged_truth_ak04_e,                          \
                 njet_charged_truth_ak04, F)                        \
    BRANCH_ARRAY(jet_charged_truth_ak04_pt,                         \
                 njet_charged_truth_ak04, F)                        \
    BRANCH_ARRAY(jet_charged_truth_ak04_eta,                        \
                 njet_charged_truth_ak04, F)                        \
    BRANCH_ARRAY(jet_charged_truth_ak04_phi,                        \
                 njet_charged_truth_ak04, F)                        \
    BRANCH_ARRAY(jet_charged_truth_ak04_area,                       \
                 njet_charged_truth_ak04, F)                        \
    BRANCH_ARRAY(jet_charged_truth_ak04_emf,                        \
                 njet_charged_truth_ak04, F)                        \
    BRANCH_ARRAY(jet_charged_truth_ak04_multiplicity,               \
                 njet_charged_truth_ak04, s)                        \
    BRANCH_ARRAY2(jet_charged_truth_ak04_width_sigma,               \
                  njet_charged_truth_ak04, 2, F)                    \
    BRANCH_ARRAY(jet_charged_truth_ak04_ptd,                        \
                 njet_charged_truth_ak04, F)                        \
    IF_EFP7(BRANCH_ARRAY2(jet_charged_truth_ak04_efp,               \
                          njet_charged_truth_ak04, 489, F))         \
    /* */                                                           \
    BRANCH_ARRAY(met_tpc, 2, D)                                     \
    BRANCH_ARRAY(met_its, 2, D)                                     \
    BRANCH_ARRAY(met_truth, 2, D)                                   \


#define B Char_t
#define b UChar_t
#define S Short_t
#define s UShort_t
#define I Int_t
#define i UInt_t
#define F Float_t
#define D Double_t
#define L Long_t
#define l ULong_t
#define O Bool_t

#define ntrigger_class              NTRIGGER_CLASS_MAX
#define ncluster                    NCLUSTER_MAX
#define ntrack                      NTRACK_MAX
#define nmuon_track                 NTRACK_MAX
#define nmc_truth                   NMC_TRUTH_MAX
#define debug_njet_ue_estimation    NJET_MAX
#define njet_ak04tpc                NJET_MAX
#define njet_ak04its                NJET_MAX
#define njet_truth_ak04             NJET_MAX
#define njet_charged_truth_ak04     NJET_MAX
#define full_emcal_ncell            EMCAL_NCELL

#define BRANCH(b, t)                            \
    t _branch_ ## b;
#define BRANCH_ARRAY(b, d, t)                   \
    t _branch_ ## b [(d)];
#define BRANCH_ARRAY2(b, d, e, t)               \
    t _branch_ ## b [(d)][(e)];
#define BRANCH_STR(b)                           \
    char _branch_ ## b[BUFSIZ];
#define BRANCH_STR_ARRAY(b, d)                  \
    std::vector<std::string> _branch_ ## b;

    MEMBER_BRANCH;

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
#undef BRANCH_STR
#undef BRANCH_STR_ARRAY

#undef ntrigger_class
#undef ncluster
#undef ntrack
#undef nmc_truth
#undef debug_njet_ue_estimation
#undef njet_ak04tpc
#undef njet_ak04its
#undef njet_truth_ak04
#undef njet_charged_truth_ak04
#undef full_emcal_ncell

#undef C
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

    Int_t _run_number_current; //!

    TFormula *_f1_ncluster_tpc_linear_pt_dep; //!
    std::vector<AliESDtrackCuts> _track_cut; //!

    AliEMCALRecoUtils *_reco_util; //!
    AliEMCALGeometry *_emcal_geometry; //!

    AliMuonTrackCuts *_muon_track_cut; //!

    size_t _ncell; //!

    std::string _emcal_geometry_filename;
    std::string _emcal_local2master_filename;

    bool _force_ue_subtraction;
    double _skim_cluster_min_e;
    double _skim_track_min_pt;
    double _skim_muon_track_min_pt;
    std::vector<double> _skim_jet_min_pt;
    double _skim_jet_average_pt;
    int _skim_multiplicity_tracklet_min_n;
    int _skim_sum_eg_ntrial;
    double _stored_track_min_pt;
    double _stored_jet_min_pt_raw;
    unsigned int _nrandom_isolation;
    bool _is_embed;

    std::vector<bool> _emcal_mask; //!

    void *_emcal_cell_position; //!

    void *_keras_model_photon_discrimination; //!

    AliAnalysisAlien *_alien_plugin; //!
    bool _metadata_filled; //!

    void loadEmcalGeometry();
    bool getEvent(AliVEvent *&event,
                  AliESDEvent *&esd_event,
                  AliAODEvent *&aod_event);
    void getContainers(AliClusterContainer *&cluster_container,
                       std::vector<AliTrackContainer*> &track_containers,
                       AliMCParticleContainer *&mc_container,
                       AliVCaloCells *&emcal_cell,
                       AliVEvent *event);
    void setTrackCuts();
    void getMultiplicityCentralityEventPlane(AliVEvent *event);
    void loadPhotonNNModel();
    void loadMC(AliAODEvent *aod_event);
    void getBeamProperties(AliVEvent *event,
                           AliESDEvent *esd_event,
                           AliAODEvent *aod_event,
                           const AliVVertex *&primary_vertex);
    void getPrimaryVertex(AliVEvent *event,
                          AliESDEvent *esd_event,
                          const AliVVertex *&primary_vertex);
    void getPrimaryVertexSPD(AliVEvent *event,
                             AliESDEvent *esd_event,
                             AliAODEvent *aod_event);
    void getPileup(AliESDEvent *esd_event,
                   AliAODEvent *aod_event);
    bool skimMultiplicityTracklet(AliVEvent *event);
    bool skimClusterE(AliClusterContainer *calo_cluster);
    void getMetadata(AliESDEvent *esd_event,
                     AliAODEvent *aod_event);
    void getPrimaryMCParticles(AliMCParticleContainer *mc_container,
                               std::vector<size_t> &stored_mc_truth_index,
                               std::vector<Int_t> &reverse_stored_mc_truth_index,
                               std::vector<Int_t> &reverse_stored_parton_algorithmic_index);
    void doTrackLoop(AliVEvent *event,
                     AliAODEvent *aod_event,
                     std::vector<AliTrackContainer*> track_containers,
                     AliMCParticleContainer *mc_container,
                     std::vector<size_t> stored_mc_truth_index,
                     const AliVVertex *primary_vertex);
    void doMCParticleLoop(AliMCParticleContainer *mc_container,
                          AliESDEvent *esd_event,
                          std::vector<Int_t> reverse_stored_mc_truth_index);
    void doClusterLoop(AliVEvent *event,
                       AliVCaloCells *emcal_cell,
                       AliClusterContainer *cluster_container,
                       std::vector<size_t> stored_mc_truth_index,
                       std::vector<bool> &cell_pass_basic_quality);
    void fillClusterBranches(AliVCaloCells *emcal_cell,
                             AliVCluster *c,
                             Int_t i, // i is cluster index in loop
                             std::vector<size_t> stored_mc_truth_index);
    void fillPhotonNNBranches(AliVCluster *c,
                              std::vector<size_t> stored_mc_truth_index,
                              AliVCaloCells *emcal_cell,
                              AliVVZERO *v0);
    void getUEJetsIsolation(AliESDEvent *esd_event,
                            AliMCParticleContainer *mc_container,
                            AliClusterContainer *cluster_container,
                            std::vector<AliTrackContainer*> track_containers,
                            std::vector<size_t> stored_mc_truth_index,
                            std::vector<Int_t> reverse_stored_mc_truth_index,
                            std::vector<Int_t> reverse_stored_parton_algorithmic_index,
                            std::vector<bool> cell_pass_basic_quality);
    void getTruthJetsAndIsolation(AliESDEvent *esd_event,
                                  AliMCParticleContainer *mc_container,
                                  AliClusterContainer *cluster_container,
                                  std::vector<size_t> stored_mc_truth_index,
                                  std::vector<Int_t> reverse_stored_mc_truth_index,
                                  std::vector<fastjet::PseudoJet> &jet_truth_ak04,
                                  std::vector<fastjet::PseudoJet> &jet_charged_truth_ak04,
                                  fastjet::ClusterSequenceArea *&cluster_sequence_truth);
    void getTruthJets(bool subtract_ue,
                      AliMCParticleContainer *mc_container,
                      std::vector<Int_t> reverse_stored_mc_truth_index,
                      std::vector<fastjet::PseudoJet> &jet_truth_ak04,
                      std::vector<fastjet::PseudoJet> &jet_charged_truth_ak04,
                      fastjet::ClusterSequenceArea *&cluster_sequence_truth);
    void getTruthIsolation(AliClusterContainer *cluster_container,
                           AliMCParticleContainer *mc_container,
                           bool subtract_ue,
                           std::vector<size_t> stored_mc_truth_index);
    void getTpcUEJetsIsolation(std::vector<AliTrackContainer*> track_containers,
                               AliClusterContainer *cluster_container,
                               std::vector<bool> cell_pass_basic_quality,
                               AliMCParticleContainer *mc_container,
                               std::vector<Int_t> reverse_stored_parton_algorithmic_index,
                               std::vector<fastjet::PseudoJet> jet_truth_ak04,
                               std::vector<fastjet::PseudoJet> jet_charged_truth_ak04,
                               fastjet::ClusterSequenceArea *cluster_sequence_truth);
    void getItsUEJetsIsolation(std::vector<AliTrackContainer*> track_containers,
                               AliClusterContainer *cluster_container,
                               std::vector<bool> cell_pass_basic_quality,
                               AliMCParticleContainer *mc_container,
                               std::vector<Int_t> reverse_stored_parton_algorithmic_index,
                               std::vector<fastjet::PseudoJet> jet_truth_ak04,
                               std::vector<fastjet::PseudoJet> jet_charged_truth_ak04,
                               fastjet::ClusterSequenceArea *cluster_sequence_truth);
    void getClusterUEIsolation(AliClusterContainer *cluster_container,
                               std::vector<bool> cell_pass_basic_quality);
    void getVoronoiAreaTpc(std::vector<AliTrackContainer*> track_containers,
                           std::vector<fastjet::PseudoJet> &particle_reco_tpc,
                           std::map<size_t, size_t> &track_reco_index_tpc,
                           std::vector<double> &particle_reco_area_tpc);
    void getVoronoiAreaIts(std::vector<AliTrackContainer*> track_containers,
                           std::vector<fastjet::PseudoJet> &particle_reco_its,
                           std::map<size_t, size_t> &track_reco_index_its,
                           std::vector<double> &particle_reco_area_its);
    void getVoronoiAreaCluster(AliClusterContainer *cluster_container,
                               std::vector<bool> cell_pass_basic_quality,
                               std::vector<fastjet::PseudoJet> &particle_reco_cluster,
                               std::map<size_t, size_t> &cluster_reco_index,
                               std::vector<double> &particle_reco_area_cluster);
    void getUETpc(std::vector<fastjet::PseudoJet> particle_reco_tpc,
                  std::vector<double> particle_reco_area_tpc,
                  std::pair<std::pair<std::vector<double>, std::vector<double>>, double> &ue_estimate_tpc);
    void getUEIts(std::vector<fastjet::PseudoJet> particle_reco_its,
                  std::vector<double> particle_reco_area_its,
                  std::pair<std::pair<std::vector<double>, std::vector<double>>, double> &ue_estimate_its);
    void getUECluster(std::vector<fastjet::PseudoJet> particle_reco_cluster,
                      std::vector<double> particle_reco_area_cluster,
                      std::pair<std::pair<std::vector<double>, std::vector<double>>, double> &ue_estimate_cluster);
    void getIsolationTpc(AliClusterContainer *cluster_container,
                         std::vector<AliTrackContainer*> track_containers,
                         std::map<size_t, size_t> track_reco_index_tpc,
                         std::vector<double> particle_reco_area_tpc,
                         std::pair<std::pair<std::vector<double>, std::vector<double>>, double> ue_estimate_tpc);
    void getIsolationIts(AliClusterContainer *cluster_container,
                         std::vector<AliTrackContainer*> track_containers,
                         std::map<size_t, size_t> track_reco_index_its,
                         std::vector<double> particle_reco_area_its,
                         std::pair<std::pair<std::vector<double>, std::vector<double>>, double> ue_estimate_its);
    void getIsolationCluster(AliClusterContainer *cluster_container,
                             std::vector<bool> cell_pass_basic_quality,
                             std::map<size_t, size_t> cluster_reco_index,
                             std::vector<double> particle_reco_area_cluster,
                             std::pair<std::pair<std::vector<double>, std::vector<double>>, double> ue_estimate_cluster);
    void getJetsTpc(std::vector<fastjet::PseudoJet> particle_reco_tpc,
                    std::vector<fastjet::PseudoJet> particle_reco_tagged_ak04_tpc,
                    std::vector<fastjet::PseudoJet> jet_truth_ak04,
                    std::vector<fastjet::PseudoJet> jet_charged_truth_ak04,
                    fastjet::ClusterSequenceArea *cluster_sequence_truth,
                    std::vector<double> particle_reco_area_tpc,
                    std::pair<std::pair<std::vector<double>, std::vector<double>>, double> ue_estimate_tpc,
                    AliClusterContainer *cluster_container,
                    std::vector<bool> cell_pass_basic_quality,
                    AliMCParticleContainer *mc_container,
                    std::vector<Int_t> reverse_stored_parton_algorithmic_index);
    void getJetsIts(std::vector<fastjet::PseudoJet> particle_reco_its,
                    std::vector<fastjet::PseudoJet> particle_reco_tagged_ak04_its,
                    std::vector<fastjet::PseudoJet> jet_truth_ak04,
                    std::vector<fastjet::PseudoJet> jet_charged_truth_ak04,
                    fastjet::ClusterSequenceArea *cluster_sequence_truth,
                    std::vector<double> particle_reco_area_its,
                    std::pair<std::pair<std::vector<double>, std::vector<double>>, double> ue_estimate_its,
                    AliClusterContainer *cluster_container,
                    std::vector<bool> cell_pass_basic_quality,
                    AliMCParticleContainer *mc_container,
                    std::vector<Int_t> reverse_stored_parton_algorithmic_index);
    bool skimJets();
    void fillCellBranches(AliVCaloCells *emcal_cell,
                          std::vector<size_t> stored_mc_truth_index);
    void fillMuonBranches();
    void fillEgNtrial();

public:
    AliAnalysisTaskNTGJ(const char *name = "NTGJ");
    AliAnalysisTaskNTGJ(const AliAnalysisTaskNTGJ &);
    AliAnalysisTaskNTGJ &operator=(const AliAnalysisTaskNTGJ &);
    ~AliAnalysisTaskNTGJ(void);
    virtual void UserCreateOutputObjects(void);
    AliEMCALRecoUtils *GetEMCALRecoUtils(void);
    void SetAliROOTVersion(const char *version);
    void SetAliPhysicsVersion(const char *version);
    void SetGridDataDir(const char *dir);
    void SetGridDataPattern(const char *pattern);
    //
    void SetEMCALGeometryFilename(const char *
                                  emcal_geometry_filename);
    void SetEMCALLocal2MasterFilename(const char *
                                      emcal_local2master_filename);
    //
    void SetForceUESubtraction(bool force_ue_subtraction = true);
    void SetSkimClusterMinE(double min_e = -INFINITY);
    void SetSkimTrackMinPt(double min_pt = -INFINITY);
    void SetSkimMuonTrackMinPt(double min_pt = -INFINITY);
    void SetSkimJetMinPt(double min_pt_1 = -INFINITY,
                         double min_pt_2 = -INFINITY,
                         double min_pt_3 = -INFINITY);
    void SetSkimJetAveragePt(double average_pt = -INFINITY);
    void SetSkimMultiplicityTrackletMinN(int min_n = INT_MIN);
    //
    void SetStoredTrackMinPt(double min_pt = -INFINITY);
    void SetStoredJetMinPtRaw(double min_pt_raw = -INFINITY);
    void SetNRandomIsolation(unsigned int nrandom_isolation = 0);
    void SetIsEmbed(bool is_embed = false);
    ClassDef(AliAnalysisTaskNTGJ, 1);

protected:
    Bool_t Run();

};

#endif // ALIANALYSISTASKPHOTONDISC_H_
