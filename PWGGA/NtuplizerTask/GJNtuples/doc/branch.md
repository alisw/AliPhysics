# Description of Ntuple Branches

## Data Types, Indexing Convention

The following C style notation are used below to describe the less obvious data types in the ntuple.

* C string: `const char *`
* C string array: `const char (*)[]`
* [IEEE 754-2008](https://dx.doi.org/10.1109/IEEESTD.2008.4610935) binary16 (&ldquo;half precision&rdquo;) stored in binary32: `__fp16`

Otherwise the x86-64 standard LP64 convention is used, with:

* 32 bit signed integer: `int`
* 64 bit signed integer: `long`

To efficiently use available space, unsigned data types are used extensively. Invalid or empty entries, e.g. the Monte Carlo truth particle index for secondary tracks, are indexed by the largest possible value of the data type (or what would be the equivalent value for &minus;1 in two&rsquo;s complement signed representation). For example with `unsigned short`, 2<sup>16</sup>&nbsp;&minus;&nbsp;1&nbsp;=&nbsp;65535 means an invalid or empty entry, but 65534 is an actual index.

## Run/Event Metadata

| Name                 | Type               | Description                               |
| -------------------- | ------------------ | ----------------------------------------- |
| `id_git`             | `const char *`     | Git blob SHA1 hash of the analyzer        |
| `version_aliroot`    | `const char *`     | AliRoot version                           |
| `version_aliphysics` | `const char *`     | AliPhysics version                        |
| `version_jec`        | `const char *`     | Jet energy correction version             |
| `grid_data_dir`      | `const char *`     | AliEn data production directory           |
| `grid_data_pattern`  | `const char *`     | AliEn data pattern for the &ldquo;find&rdquo; command |
| `beam_particle`      | `int[2]`           | _Z_&nbsp;+&nbsp;1000&#8239;_A_ for each beam particle |
| `ntrigger_class`     | `unsigned long`    | Number of trigger classes                 |
| `trigger_class`      | `const char (*)[ntrigger_class]` | Trigger class names         |
| `run_number`         | `int`              | ALICE run number                          |
| `trigger_mask`       | `unsigned long[2]` | Trigger mask, each storing 50 bits        |
| `mixed_event`        | `char`             | True if ntuple underwent event mixing     |
| `multiplicity_v0`    | `float[64]`        | V0 multiplicity, per each of the 64 channels |
| `centrality_v0m`     | `float`            | V0M centrality, alias for `centrality[0]` |
| `centrality`         | `float[9]`         | Centrality using: {V0M, CL0, CL1, V0Mplus05, V0Mplus10, V0Mminus05, V0Mminus10, SPDClustersCorr, SPDTracklets} |
| `event_plane_psi_v0` | `float[3]`         | V0 event plane angle &Psi;<sub>_n_</sub>, _n_&nbsp;=&nbsp;1, &hellip; 3, for the directed/elliptic/triangular flows  |
| `event_plane_q_v0`   | `double[3][2]`     | V0 event plane _Q_-vector _Q_<sub>_n_</sub>, _n_&nbsp;=&nbsp;1, &hellip; 3, for the directed/elliptic/triangular flows |
| `has_misalignment_matrix` | `bool`        | True if the Electromagnet Calorimeter (EMCAL) misalignment matrix was loaded |
| `cell_eta`           | `__fp16[ncluster]` | EMCAL cell pseudorapidity _&eta;_                                           |
| `cell_phi`           | `__fp16[ncluster]` | EMCAL cell azimuth _&straightphi;_                                          |
| `cell_voronoi_area`  | `__fp16[ntrack]`   | The Voronoi diagram area occupied by the EMCAL cell in the _&eta;_&ndash;_&straightphi;_-plane |
| `primary_vertex`     | `double[3]`        | The _x_, _y_, and _z_ compoents of the primary vertex position (cm) |
| `primary_vertex_sigma` | `double[3]`      | The _x_, _y_, and _z_ uncertainties of the primary vertex position (cm) |
| `primary_vertex_ncontributor` | `int`     | The number of tracks that contributed to the primary vertex determination |
| `primary_vertex_spd` | `double[3]`        | The _x_, _y_, and _z_ compoents of the Silicon Pixel Detector (SPD) primary vertex position (cm) |
| `primary_vertex_spd_sigma` | `double[3]`  | The _x_, _y_, and _z_ uncertainties of the SPD primary vertex position (cm) |
| `primary_vertex_spd_ncontributor` | `int` | The number of tracklets that contributed to the SPD primary vertex determination |
| `npileup_vertex_spd`   | `int`            | The number of SPD pileup vertices |
| `pileup_vertex_spd_ncontributor` | `int`  | The number tracklets that contributed to the SPD pileup vertices determination |
| `is_pileup_from_spd_3_08` | `bool`        | Event is SPD pileup with multiple 3 contributors vertices and z separation < 0.8&nbsp;cm |
| `is_pileup_from_spd_5_08` | `bool`        | Event is SPD pileup with multiple 5 contributors vertices and z separation < 0.8&nbsp;cm |
| `ncluster_tpc`       | `int`              | The number of TPC clusters |
| `event_selected`     | `bool`             | Event is selected by the input handler (e.g. the ALICE DPG&rsquo;s physics selection, if ntuplized with `AliPhysicsSelectionTask`) |
| `eg_signal_process_id` | `int`            | The code of the process of the current event, HepMC&rsquo;s `GenEvent::signal_process_id()` and e.g. PYTHIA 8&rsquo;s `Info::code()` |
| `eg_mpi`             | `int`              | The number of hard interactions of the current event, HepMC&rsquo;s `GenEvent::mpi()` and e.g. PYTHIA 8&rsquo;s `Info::nMPI()` |
| `eg_pt_hat`          | `float`            | The rest frame transverse momentum _p&#770;_<sub>&perp;</sub> (GeV/_c_) of the current event, e.g. PYTHIA 8&rsquo;s `Info::pTHat()` |
| `eg_cross_section`   | `float`            | Cross section _&sigma;_ (mb) of the current event, e.g. PYTHIA 8&rsquo;s `Info::sigmaGen()` and HERWIG 7&rsquo;s `EventHandler::integratedXSec() / Units::millibarn` |
| `eg_weight`          | `float`            | Weight of the current event, HepMC&rsquo;s `GenEvent::weights().front()` and e.g. PYTHIA 8&rsquo;s `Info::weight()` and HERWIG 7&rsquo;s `EventHandler::weight() / MEBase::reWeight()` |
| `eg_primary_vertex`  | `float[3]`         | The generator truth primary vertex, HepMC&rsquo;s `GenEvent::signal_process_vertex()->point3d()` |
| `eg_ntrial`          | `int`              | Number of trial until the current event is generated (ALICE SW specific) |
| `eg_scale_pdf`       | `float`            | The scale _Q_ for the of the current event, as used for the calculation of parton distribution functions (PDFs), HepMC&rsquo;s `GenEvent::event_scale()` and e.g. PYTHIA 8&rsquo;s `Info::scalup()` (not stored by ALICE SW) |
| `eg_alpha_qcd`       | `float`            | The strong coupling constant, _&alpha;_<sub>S</sub>, of the current event, HepMC&rsquo;s `GenEvent::alphaQCD()` and e.g. PYTHIA 8&rsquo;s `Info::alphaS()`, HERWIG 7&rsquo;s `EventHandler::lastAlphaS()` (not stored by ALICE SW) |
| `eg_alpha_qed`       | `float`            | The electromagnetic coupling constant, _&alpha;_<sub>EM</sub>, of the current event, HepMC&rsquo;s `GenEvent::alphaQED()` and e.g. PYTHIA 8&rsquo;s `Info::alphaEM()`, HERWIG 7&rsquo;s `EventHandler::lastAlphaEM()` (not stored by ALICE SW) |
| `eg_pdf_id`          | `int[2]`           | The flavor code of the two incoming partons, HepMC&rsquo;s `PdfInfo::id1()`/`id2()` and e.g. PYTHIA 8&rsquo;s `Info::id1pdf()`/`id2pdf()` (not stored by ALICE SW) |
| `eg_pdf_x`           | `float[2]`         | The longitudinal momentum fraction _x_ of the two incoming partons, HepMC&rsquo;s `PdfInfo::x1()`/`x2()` and e.g. PYTHIA 8&rsquo;s `Info::x1pdf()`/`x2pdf()`, HERWIG 7&rsquo;s `EventHandler::lastX1()`/`lastX2()` (not stored by ALICE SW) |
| `eg_pdf_x_pdf`       | `float[2]`         | The parton densities _x_&#8239;_f_(_x_, _Q_<sup>2</sup>) of the two incoming partons, HepMC&rsquo;s `PdfInfo::pdf1()`/`pdf2()` and e.g. PYTHIA 8&rsquo;s `Info::pdf1()`/`pdf2()` (not stored by ALICE SW) |

## EMCAL Clusters

| Name                     | Type                           | Description                                                              |
| ------------------------ | ------------------------------ | ------------------------------------------------------------------------ |
| `ncluster`               | `unsigned long`                | Number of EMCAL clusters                                                 |
| `cluster_e`              | `__fp16[ncluster]`             | Cluster energy _E_ (GeV) (with masses from _dE_/_dx_)                    |
| `cluster_pt`             | `__fp16[ncluster]`             | Cluster transverse momentum _p_<sub>T</sub> (GeV/_c_)                    |
| `cluster_eta`            | `__fp16[ncluster]`             | Cluster pseudorapidity _&eta;_                                           |
| `cluster_phi`            | `__fp16[ncluster]`             | Cluster azimuth _&straightphi;_                                          |
| `cluster_lambda`         | `__fp16[ncluster][2]`          | Minor (`[0]`) and major (`[1]`) axes of the cluster dispersion (_&lambda;_<sub>0</sub>)<sup>2</sup> and (_&lambda;_<sub>1</sub>)<sup>2</sup> |
| `cluster_tof`            | `__fp16[ncluster]`             | Cluster time-of-flight _T_<sub>0</sub> (ns)                               |
| `cluster_ncell`          | `int[ncluster]`                | Number of cells (towers) in the cluster                                  |
| `cluster_cell_id_max`    | `unsigned short[ncluster]`     | Index of the cell (tower) in the cluster with the highest energy         |
| `cluster_e_max`          | `__fp16[ncluster]`             | Energy of the cell (tower) in the cluster with the highest energy (GeV)  |
| `cluster_e_cross`        | `__fp16[ncluster]`             | Energy sum of the 4 cross-shaped cells (towers) adjacent to the cell in the cluster with the highest energy (GeV) |
| `cluster_nmc_truth`      | `unsigned short[ncluster]`     | Number of matched Monte Carlo (MC) truth particles, in `mc_truth_`*      |
| `cluster_mc_truth_index` | `unsigned short[ncluster][32]` | Indices of the matched MC truth particles, in `mc_truth_`*               |
| `cluster_iso_tpc_01`     | `__fp16[ncluster]`             | Time Projection Chamber (TPC) isolation transverse momentum _p_<sub>T,iso</sub> (GeV/_c_) with _R_<sub>0</sub>&nbsp;=&nbsp;0.1 |
| `cluster_iso_tpc_02`     | `__fp16[ncluster]`             | TPC isolation transverse momentum _p_<sub>T,iso</sub> (GeV/_c_) with _R_<sub>0</sub>&nbsp;=&nbsp;0.2 |
| `cluster_iso_tpc_03`     | `__fp16[ncluster]`             | TPC isolation transverse momentum _p_<sub>T,iso</sub> (GeV/_c_) with _R_<sub>0</sub>&nbsp;=&nbsp;0.3 |
| `cluster_iso_tpc_04`     | `__fp16[ncluster]`             | TPC isolation transverse momentum _p_<sub>T,iso</sub> (GeV/_c_) with _R_<sub>0</sub>&nbsp;=&nbsp;0.4 |
| `cluster_frixione_tpc_04_02` | `__fp16[ncluster]`         | TPC Frixione max<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;0.2 |
| `cluster_frixione_tpc_04_05` | `__fp16[ncluster]`         | TPC Frixione max<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;0.5 |
| `cluster_frixione_tpc_04_10` | `__fp16[ncluster]`         | TPC Frixione max<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;1.0 |
| `cluster_frixione_its_04_02` | `__fp16[ncluster]`         | ITS Frixione max<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;0.2 |
| `cluster_frixione_its_04_05` | `__fp16[ncluster]`         | ITS Frixione max<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;0.5 |
| `cluster_frixione_its_04_10` | `__fp16[ncluster]`         | ITS Frixione max<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;1.0 |
| `cluster_anti_frixione_tpc_04_02` | `__fp16[ncluster]`    | TPC &ldquo;anti-Frixione&rdquo; min<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;0.2 |
| `cluster_anti_frixione_tpc_04_05` | `__fp16[ncluster]`    | TPC &ldquo;anti-Frixione&rdquo; min<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;0.5 |
| `cluster_anti_frixione_tpc_04_10` | `__fp16[ncluster]`    | TPC &ldquo;anti-Frixione&rdquo; min<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;1.0 |
| `cluster_anti_frixione_its_04_02` | `__fp16[ncluster]`    | ITS &ldquo;anti-Frixione&rdquo; min<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;0.2 |
| `cluster_anti_frixione_its_04_05` | `__fp16[ncluster]`    | ITS &ldquo;anti-Frixione&rdquo; min<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;0.5 |
| `cluster_anti_frixione_its_04_10` | `__fp16[ncluster]`    | ITS &ldquo;anti-Frixione&rdquo; min<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;1.0 |
| `cluster_iso_01_truth`   | `__fp16[ncluster]`             | MC truth isolation transverse momentum _p_<sub>T,iso</sub> (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.1 |
| `cluster_iso_02_truth`   | `__fp16[ncluster]`             | MC truth isolation transverse momentum _p_<sub>T,iso</sub> (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.2 |
| `cluster_iso_03_truth`   | `__fp16[ncluster]`             | MC truth isolation transverse momentum _p_<sub>T,iso</sub> (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.3 |
| `cluster_iso_04_truth`   | `__fp16[ncluster]`             | MC truth isolation transverse momentum _p_<sub>T,iso</sub> (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4 |
| `cluster_frixione_04_02_truth` | `__fp16[ncluster]`       | MC truth Frixione max<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;0.2 |
| `cluster_frixione_04_05_truth` | `__fp16[ncluster]`       | MC truth Frixione max<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;0.5 |
| `cluster_frixione_04_10_truth` | `__fp16[ncluster]`       | MC truth Frixione max<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_) (GeV/_c_) with _&delta;_<sub>0</sub>&nbsp;=&nbsp;0.4, _n_&nbsp;=&nbsp;1.0 |

To apply the Frixione isolation _p_<sub>T,iso</sub>(_&delta;_)&nbsp;&equiv;&nbsp;&#8721;<sub>_i_</sub>&#8239;_p_<sub>T,_i_</sub>&#8239;_&theta;_(_&delta;_&nbsp;&minus;&nbsp;_R_<sub>_i&gamma;_</sub>)&nbsp;&le;&nbsp;&#120039;(_&delta;_) for all _&delta;_&nbsp;&le;&nbsp;_&delta;_<sub>0</sub>, cut on max<sub>_&delta;_</sub>&#8239;_p_<sub>T,iso</sub>(_&delta;_)_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>/&#120039;(_&delta;_)&nbsp;&le;&nbsp;_E_<sub>_&gamma;_</sub>_&straightepsilon;_<sub>_&gamma;_</sub>. See [S. Frixione, Phys. Lett. B 429, 369 (1998)](https://dx.doi.org/10.1016/S0370-2693(98)00454-7) [[arXiv:hep-ph/9801442]](https://arxiv.org/abs/hep-ph/9801442) for details and notation.

## EMCAL Cells

A full 17664 (=&nbsp;10&nbsp;&times;&nbsp;48&nbsp;&times;&nbsp;24 + 6&nbsp;&times;&nbsp;32&nbsp;&times;&nbsp;24 + 4&nbsp;&times;&nbsp;48&nbsp;&times;&nbsp;8, for the 10 full-`iphi` EMCAL, 6 full-`iphi` Di-Jet Calorimeter (DCAL), and 4 one-third-`iphi` EMCAL/DCAL supermodules) entries array is stored per event.

| Name                  | Type                    | Description                                                           |
| --------------------- | ----------------------- | --------------------------------------------------------------------- |
| `cell_energy`         | `__fp16[17664]`         | Cell energy _E_ (GeV)                                                 |
| `cell_tof`            | `__fp16[17664]`         | Cell time-of-flight _T_<sub>0</sub> (ns)                               |
| `cell_mc_truth_index` | `unsigned short[17664]` | Index of the matched Monte Carlo (MC) truth particle, in `mc_truth_`* |

## Central Region Tracks

| Name                           | Type                     | Description                                                          |
| -------------------------------| ------------------------ | -------------------------------------------------------------------- |
| `ntrack`                       | `unsigned long`          | Number of central region tracks                                      |
| `track_e`                      | `__fp16[ntrack]`         | Track energy _E_ (GeV) (with masses from _dE_/_dx_)                  |
| `track_pt`                     | `__fp16[ntrack]`         | Track transverse momentum _p_<sub>T</sub> (GeV/_c_)                  |
| `track_eta`                    | `__fp16[ntrack]`         | Track pseudorapidity _&eta;_                                         |
| `track_phi`                    | `__fp16[ntrack]`         | Track azimuth _&straightphi;_                                        |
| `track_quality`                | `unsigned char[ntrack]`  | Track quality bit mask, |
| `track_tpc_dedx`               | `__fp16[ntrack]`         | Track _dE_/_dx_ (arbitrary units)                                    |
| `track_tpc_length_active_zone` | `__fp16[ntrack]`         | Track length in the active zone of the TPC (cm)                      |
| `track_tpc_xrow`               | `unsigned char[ntrack]`  | Number of TPC rows crossed by the track                              |
| `track_tpc_ncluster`           | `unsigned char[ntrack]`  | Number of (actual) TPC clusters in the track                         |
| `track_tpc_ncluster_dedx`      | `unsigned char[ntrack]`  | Number of TPC clusters used to obtain the _dE_/_dx_                  |
| `track_tpc_ncluster_findable`  | `unsigned char[ntrack]`  | Number of potentially findable TPC clusters                          |
| `track_its_ncluster`           | `unsigned char[ntrack]`  | Number of (actual) Inner Tracking System (ITS) clusters                                      |
| `track_dca_xy`                 | `__fp16[ntrack]`         | Distance of closest approach to the primary vertex in the _x_&ndash;_y_-plane (cm) |
| `track_dca_z`                  | `__fp16[ntrack]`         | Distance of closest approach to the primary vertex in the _z_-direction (cm) |
| `track_mc_truth_index`         | `unsigned short[ntrack]` | Index of the matched Monte Carlo truth particle, in `mc_truth_`*     |
| `track_voronoi_area`           | `__fp16[ntrack]`         | The Voronoi diagram area occupied by the track in the _&eta;_&ndash;_&straightphi;_-plane |

## Muon Tracks

| Name                                | Type                  | Description                                                        |
| ----------------------------------- | --------------------- | ------------------------------------------------------------------ |
| `nmuon_track`                       | `unsigned long`       | Number of muon tracks                                              |
| `muon_track_e`                      | `__fp16[nmuon_track]` | Muon track energy _E_ (GeV)                                        |
| `muon_track_pt`                     | `__fp16[nmuon_track]` | Muon track transverse momentum _p_<sub>T</sub> (GeV/_c_)           |
| `muon_track_eta`                    | `__fp16[nmuon_track]` | Muon track pseudorapidity _&eta;_                                  |
| `muon_track_phi`                    | `__fp16[nmuon_track]` | Muon track azimuth _&straightphi;_                                 |
| `muon_track_p_dca`                  | `__fp16[nmuon_track]` | Muon track _p_&nbsp;&times;&nbsp;DCA (GeV&nbsp;cm/_c_) (see below) |
| `muon_track_sigma_p_dca`            | `__fp16[nmuon_track]` | Uncorrected muon track _&sigma;_(_p_&nbsp;&times;&nbsp;DCA) (see below) |
| `muon_track_sigma_slope_p`          | `__fp16[nmuon_track]` | (see below)                                                        |
| `muon_track_distance_sigma_slope_p` | `__fp16[nmuon_track]` | (see below)                                                        |

To make an arbitrary, resolution corrected &ldquo;n&sigma;&rdquo; cut on _p_&nbsp;&times;&nbsp;DCA, the boolean `pass_cut` is calculated as:

```C++
const float nrp = nsigma_p_dca * muon_track_delta_sagitta_p[i]
const float p_resolution_effect = muon_track_sigma_p_dca[i] / (1 - nrp / (1 + nrp))
const float slope_resolution_effect = muon_track_distance_sigma_slope_p[i]
const float sigma_p_dca_resolution_effect = sqrt(std::pow(p_resolution_effect, 2) + std::pow(slope_resolution_effect, 2))
const bool pass_cut = muon_track_p_dca[i] < sigma_p_dca_resolution_effect
```

## Monte Carlo Truth Particles

| Name                       | Type                | Description                                            |
| -------------------------- | ------------------- | ------------------------------------------------------ |
| `nmc_truth`                | `unsigned long`     | Number of Monte Carlo (MC) truth particles             |
| `mc_truth_e`               | `__fp16[nmc_truth]` | Particle energy _E_ (GeV)                              |
| `mc_truth_pt`              | `__fp16[nmc_truth]` | Particle transverse momentum _p_<sub>T</sub> (GeV/_c_) |
| `mc_truth_eta`             | `__fp16[nmc_truth]` | Particle pseudorapidity _&eta;_                        |
| `mc_truth_phi`             | `__fp16[nmc_truth]` | Particle azimuth _&straightphi;_                       |
| `mc_truth_pdg_code`        | `__fp16[nmc_truth]` | The particle species, according to Particle Data Group&rsquo;s MC Number Scheme |
| `mc_truth_status`          | `__fp16[nmc_truth]` | The HepMC status code of the particle                  |
| `mc_truth_generator_index` | `__fp16[nmc_truth]` | Index of the generator that produced the particle      |
| `mc_truth_voronoi_area`    | `__fp16[nmc_truth]` | The Voronoi diagram area occupied by the particle in the _&eta;_&ndash;_&straightphi;_-plane |
| `charged_mc_truth_voronoi_area` | `__fp16[nmc_truth]` | The Voronoi diagram area occupied by the charged particle (and among only the charged MC truth particles) in the _&eta;_&ndash;_&straightphi;_-plane |

## Jets

The `jet_` prefix is generic, and in reality there are multiple algorithms, e.g. it is replaced with `jet_ak04tpc_` for anti-_k<sub>t</sub>_ _D_&nbsp;=&nbsp;0.4 with TPC tracks.

| Name                      | Type                   | Description                                                         |
| --------------------------| ---------------------- | ------------------------------------------------------------------- |
| `njet`                    | `unsigned long`        | Number of reconstructed jets                                        |
| `jet_e_raw`               | `__fp16[njet]`         | Raw (observed) jet energy _E_ (GeV)                                 |
| `jet_e`                   | `__fp16[njet]`         | Calibrated jet energy _E_ (GeV)                                     |
| `jet_e_charged`           | `__fp16[njet]`         | Calibrated, charged-particles-only jet energy _E_ (GeV)             |
| `jet_pt_raw_ue`           | `__fp16[njet]`         | The estimated underlying event (UE) of the raw (observed) jet transverse momentum _p_<sub>T</sub> (GeV/_c_) |
| `jet_pt_raw`              | `__fp16[njet]`         | Raw (observed) jet transverse momentum _p_<sub>T</sub> (GeV/_c_), after subtracting the UE (`jet_pt_raw_ue`) |
| `jet_pt`                  | `__fp16[njet]`         | Calibrated jet transverse momentum _p_<sub>T</sub> (GeV/_c_)        |
| `jet_pt_charged`          | `__fp16[njet]`         | Calibrated, charged-particles-only jet transverse momentum _p_<sub>T</sub> (GeV/_c_) |
| `jet_eta_raw`             | `__fp16[njet]`         | Raw (observed) jet pseudorapidity _&eta;_                           |
| `jet_eta`                 | `__fp16[njet]`         | Calibrated jet pseudorapidity _&eta;_                               |
| `jet_phi`                 | `__fp16[njet]`         | Raw (observed) jet azimuth _&straightphi;_                          |
| `jet_area_raw`            | `__fp16[njet]`         | Raw (observed) jet area (units of &Delta;_&eta;_&nbsp;&times;&nbsp;&Delta;_&straightphi;_) |
| `jet_area`                | `__fp16[njet]`         | Calibrated jet area (units of &Delta;_&eta;_&nbsp;&times;&nbsp;&Delta;_&straightphi;_) |
| `jet_emf_raw`             | `__fp16[njet]`         | Raw (observed) jet electromagnetic fraction (EMF)                   |
| `jet_emf`                 | `__fp16[njet]`         | Calibrated jet electromagnetic fraction (EMF)                       |
| `jet_multiplicity_raw`    | `unsigned short[njet]` | Raw (observed) jet multiplicity                                     |
| `jet_multiplicity`        | `__fp16[njet]`         | Calibrated jet multiplicity                                         |
| `jet_width_sigma_raw`     | `__fp16[njet][2]`      | Raw (observed) normalized major (`[0]`) and minor (`[1]`) axes of the 2-dimensional &#10216;_j_<sub>T</sub>&#10217; |
| `jet_width_sigma`         | `__fp16[njet][2]`      | Calibrated, normalized major (`[0]`) and minor (`[1]`) axes of the 2-dimensional &#10216;_j_<sub>T</sub>&#10217; |
| `jet_ptd_raw`             | `__fp16[njet]`         | Raw (observed) jet _p_<sub>T</sub>_D_                               |
| `jet_ptd`                 | `__fp16[njet]`         | Calibrated jet _p_<sub>T</sub>_D_                                   |
| `jet_efp`                 | `__fp16[njet][489]`    | Scaled raw (observed) jet energy flow polynomials with order _n_&nbsp;&le;&nbsp;7, see below. |
| `jet_e_truth`             | `__fp16[njet]`         | Alias for `jet_truth_e[jet_truth_index_z_reco[i][0]]`               |
| `jet_pt_truth`            | `__fp16[njet]`         | Alias for `jet_truth_pt[jet_truth_index_z_reco[i][0]]`              |
| `jet_eta_truth`           | `__fp16[njet]`         | Alias for `jet_truth_eta[jet_truth_index_z_reco[i][0]]`             |
| `jet_phi_truth`           | `__fp16[njet]`         | Alias for `jet_truth_phi[jet_truth_index_z_reco[i][0]]`             |
| `jet_area_truth`          | `__fp16[njet]`         | Alias for `jet_truth_area[jet_truth_index_z_reco[i][0]]`            |
| `jet_emf_truth`           | `__fp16[njet]`         | Alias for `jet_truth_emf[jet_truth_index_z_reco[i][0]]`             |
| `jet_multiplicity_truth`  | `unsigned short[njet]` | Alias for `jet_truth_multiplicity[jet_truth_index_z_reco[i][0]]`    |
| `jet_width_sigma_truth`   | `__fp16[njet][2]`      | Alias for `jet_truth_width_sigma[jet_truth_index_z_reco[i][0]]`     |
| `jet_ptd_truth`           | `__fp16[njet]`         | Alias for `jet_truth_ptd[jet_truth_index_z_reco[i][0]]`             |
| `jet_truth_index_z_truth` | `int[njet][2]`         | Index of the truth jet with the two highest longitudinal momentum fraction _z_<sub>truth</sub>, relative to the truth jet 3-momenta, in descending _z_<sub>truth</sub> order |
| `jet_truth_z_truth`       | `__fp16[njet][2]`      | The longitudinal momentum fraction _z_<sub>truth</sub> of the two highest truth jets, relative to the truth jet 3-momenta, in descending _z_<sub>truth</sub> order |
| `jet_truth_index_z_reco`  | `int[njet][2]`         | Index of the truth jet with the two highest longitudinal momentum fraction _z_<sub>reco</sub>, relative to the reconstructed jet 3-momenta, in descending _z_<sub>reco</sub> order |
| `jet_truth_z_reco`        | `__fp16[njet][2]`      | The longitudinal momentum fraction _z_<sub>reco</sub> of the two highest truth jets, relative to the reconstructed jet 3-momenta, in descending _z_<sub>reco</sub> order |
| `njet_truth`              | `unsigned long`        | Number of Monte Carlo (MC) truth jets                               |
| `jet_truth_e`             | `__fp16[njet]`         | MC truth jet energy _E_ (GeV)                                       |
| `jet_truth_pt_ue`         | `__fp16[njet]`         | MC truth underlying event (UE) of the jet momentum _p_<sub>T</sub> (GeV/_c_) |
| `jet_truth_pt`            | `__fp16[njet]`         | MC truth jet transverse momentum _p_<sub>T</sub> (GeV/_c_)          |
| `jet_truth_eta`           | `__fp16[njet]`         | MC truth jet pseudorapidity _&eta;_                                 |
| `jet_truth_phi`           | `__fp16[njet]`         | MC truth jet azimuth _&straightphi;_                                |
| `jet_truth_area`          | `__fp16[njet]`         | MC truth jet area (units of &Delta;_&eta;_&nbsp;&times;&nbsp;&Delta;_&straightphi;_) |
| `jet_truth_emf`           | `__fp16[njet]`         | MC truth jet electromagnetic fraction (EMF)                         |
| `jet_truth_multiplicity`  | `__fp16[njet]`         | MC truth jet multiplicity                                           |
| `jet_truth_width_sigma`   | `__fp16[njet][2]`      | MC truth normalized major (`[0]`) and minor (`[1]`) axes of the 2-dimensional &#10216;_j_<sub>T</sub>&#10217; |
| `jet_truth_ptd`           | `__fp16[njet]`         | MC truth jet _p_<sub>T</sub>_D_                                     |
| `met_tpc`                 | `double[2]`            | The _x_ and _y_ components of the TPC _&#582;_<sub>T</sub> (GeV)    |
| `met_truth`               | `double[2]`            | The _x_ and _y_ components of the MC truth _&#582;_<sub>T</sub> (GeV) |

### Notes on energy flow polynomials

The 489 values per jet are filled according to the order of the `EnergyFlow` package (https://pypi.org/project/EnergyFlow/). To list the form of the Einstein sum and opening angle exponents, install the package and run:

```Python
import energyflow
for e in energyflow.EFPSet('d<=7').efpelems: print(e.einstr, e.weights)
```

To avoid denormalized half precision numbers, the values are scaled by 2<sup>15</sup>&nbsp;=&nbsp;32768, such that an actual EFP value of 1 is stored as 32768. The _&chi;_&nbsp;=&nbsp;4 entries (array indices 105, 106, 235, 239) corresponding those requiring a rank-4 tensor trace are only filled if there are less than 200 particles inside the jet.
