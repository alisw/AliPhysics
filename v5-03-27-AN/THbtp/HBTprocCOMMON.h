#ifndef HBTprocessorCOMMON
#define HBTprocessorCOMMON


enum
  {
    TRK_MAXLEN    = 25000,
    TRK2_MAXLEN   = 25000,
    MAX_TRK_SEC   = 150,
    SEC_MAXLEN    = 28000,
    MAX_TRK_SEC2  = 150,
    SEC_MAXLEN2   = 28000,
    MAX_TRK_SAVE  = 150,
    PART_MAXLEN   = 50,
    MAX_C2_1D     = 100, 
    MAX_C2_3D     = 10,
    MAX_EVENTS    = 100,
    MAX_C2_COUL   = 288,
    MAX_H_1D      = 100,
    MAX_H_3D      = 10
  };


extern "C" {


#ifndef __CINT__
#define f2cFortran
#ifndef __CFORTRAN_LOADED
#include "cfortran.h"
#endif
#endif
 
#ifndef __CINT__

/*       common/parameters/
     1   ref_control,print_full,print_sector_data,n_pid_types,pid(2),
     2   n_events,switch_1d,switch_3d,switch_type,switch_coherence,
     3   switch_coulomb,switch_fermi_bose,n_part_1_trk,n_part_2_trk,
     4   n_part_tot_trk,n_part_1_trk2,n_part_2_trk2,n_part_tot_trk2,
     5   n_part_used_1_trk,n_part_used_2_trk,n_part_used_1_trk2,
     6   n_part_used_2_trk2,n_part_used_1_ref,n_part_used_2_ref,
     7   n_part_used_1_inc,n_part_used_2_inc,num_pairs_like,
     8   num_pairs_unlike,num_pairs_like_ref,num_pairs_unlike_ref,
     9   num_pairs_like_inc,num_pairs_unlike_inc,event_line_counter,
     1   maxit,irand,file10_line_counter,
     2   lambda,R_1d,Rside,Rout,Rlong,Rperp,Rparallel,R0,Q0,deltap,
     3   delchi,pi,rad,hbc,chisq_wt_like_1d,chisq_wt_unlike_1d,
     4   chisq_wt_like_3d_fine,chisq_wt_unlike_3d_fine,
     5   chisq_wt_like_3d_coarse,chisq_wt_unlike_3d_coarse,
     6   chisq_wt_hist1_1,chisq_wt_hist1_2,mass1,mass2,trk_accep

*/ 
 


 typedef struct //PARAMETERS
   {
      Int_t ALICE;
      Int_t errorcode;
      
      Int_t ref_control;
      Int_t print_full;             // Full print out option - each event
      Int_t print_sector_data;      // Print sector overflow diagnostics
      Int_t n_pid_types;            // # particle ID types to correlate
      Int_t pid[2];                    // Geant particle ID #s, max of 2 types
      Int_t n_events ;              // # events in input event text file
      Int_t switch_1d;              // Include 1D correlations
      Int_t switch_3d;              // Include 3D correlations
      Int_t switch_type ;           // For like, unlike or both PID pairs
      Int_t switch_coherence;       // To include incoh/coher mixed source
      Int_t switch_coulomb;         // Coulomb correction selection options
      Int_t switch_fermi_bose;      // For fermions or bosons

//   Numbers of particles and pairs:

      Int_t n_part_1_trk;           // Total # PID #1 in 'trk', all flags
      Int_t n_part_2_trk;           // Total # PID #2 in 'trk', all flags
      Int_t n_part_tot_trk;         // Total # all part. in 'trk', all flgs
      Int_t n_part_used_1_trk;      // # PID#1, used (flag=0) in 'trk'
      Int_t n_part_used_2_trk;      // # PID#2, used (flag=0) in 'trk'

      Int_t n_part_1_trk2;          // Total # PID #1 in 'trk2', all flags
      Int_t n_part_2_trk2;          // Total # PID #2 in 'trk2', all flags
      Int_t n_part_tot_trk2;        // Total # all part. in 'trk2', all flgs
      Int_t n_part_used_1_trk2;     // # PID#1, used (flag=0) in 'trk2'
      Int_t n_part_used_2_trk2;     // # PID#2, used (flag=0) in 'trk2'

      Int_t n_part_used_1_ref;      // # PID#1, used (flag=0) in Reference
      Int_t n_part_used_2_ref;      // # PID#2, used (flag=0) in Reference
      Int_t n_part_used_1_inc;      // # PID#1, used (flag=0) in Inclusive 
      Int_t n_part_used_2_inc;      // # PID#2, used (flag=0) in Inclusive

      Int_t num_pairs_like;         // # like pairs used (flag=0) in fit 
      Int_t num_pairs_unlike;       // # unlike pairs used (flag=0) in fit
      Int_t num_pairs_like_ref;     // # like pairs used (flag=0) in Ref. 
      Int_t num_pairs_unlike_ref;   // # unlike pairs used (flag=0) in Ref. 
      Int_t num_pairs_like_inc;     // # like pairs used (flag=0) in Incl. 
      Int_t num_pairs_unlike_inc;   // # unlike pairs used (flag=0) in Incl. 

//   Counters:

      Int_t event_line_counter;     // Input event text file line counter
      Int_t maxit;                  // Max # iterations in track adjustment
      Int_t irand;                  // Random # starting seed (Def=12345)      
      Int_t file10_line_counter;    // Output, correlated event text file
//                                    //    line counter

//   Correlation Model Parameters:

      Float_t    lambda;               // Chaoticity parameter
      Float_t    R_1d;                   // Spherical source radius (fm)
      Float_t    Rside;                  // 3D Bertsch-Pratt source 'side' R (fm)
      Float_t    Rout;                   // 3D Bertsch-Pratt source 'out'  R (fm)
      Float_t    Rlong;                  // 3D Bertsch-Pratt source 'long' R (fm)
      Float_t    Rperp;                  // 3D YKP source transverse radius  (fm)
      Float_t    Rparallel;              // 3D YKP source longitudinal radius(fm)
      Float_t    R0;                     // 3D YKP source emission time durat(fm)
      Float_t    Q0;                     // NA35 Coulomb parameter (GeV/c) or
//                                    // Coul radius for Pratt finite src (fm)

//   Search Control Parameters:


      Float_t    deltap;                 // Max limit for x,y,z momt shifts(GeV/c)
      Float_t    delchi;                 // Min% change in Chi-Sq to stop iterat.

      Float_t pi;// = 3.141592654; // PI 
      
      Float_t rad;// = 180.0/3.14;                    // radian = 180.0/pi
      Float_t hbc;// = 0.19732891;                       // h-bar-c (GeV*fm)

//   Chi-Square Values:

      Float_t    chisq_wt_like_1d;          // 1D, Like pairs
      Float_t    chisq_wt_unlike_1d;        // 1D, Unlike pairs
      Float_t    chisq_wt_like_3d_fine;     // 3D, Like pairs, Fine Mesh
      Float_t    chisq_wt_unlike_3d_fine;   // 3D, Unlike pairs, Fine Mesh
      Float_t    chisq_wt_like_3d_coarse;   // 3D, Like pairs, Coarse Mesh
      Float_t    chisq_wt_unlike_3d_coarse; // 3D, Unlike pairs, Coarse Mesh
      Float_t    chisq_wt_hist1_1;          // One-body, particle ID type #1
      Float_t    chisq_wt_hist1_2;          // One-body, particle ID type #2

//   Particle Masses:

      Float_t    mass1, mass2;           // Particle ID# 1 and 2 masses (GeV)

//   Constants:


//     parameter (pi = 3.141592654)
//     parameter (hbc = 0.19732891)

//  Random Track Selection Fraction, for low multiplicity particles

      Float_t trk_accep;                 // ranges from 0.0 -> 1.0

  }  HBTprocParamsCommon;
 
#define PARAMETERS COMMON_BLOCK(PARAMETERS,parameters)

COMMON_BLOCK_DEF(HBTprocParamsCommon,PARAMETERS);

//HBTprocParamsCommon PARAMETERS;


/************************************************************************************************/ 
/************************************************************************************************/ 

 typedef struct //MESH COMMON BLOCK
  {
//   One-Body Mesh:

      Int_t n_pt_bins;             // # one-body pt bins
      Int_t n_phi_bins;            // # one-body phi bins
      Int_t n_eta_bins;            // # one-body eta bins
     
      Int_t n_1d_fine;             // # bins for 1D, Fine Mesh
      Int_t n_1d_coarse;           // # bins for 1D, Coarse Mesh
      Int_t n_1d_total;            // Total # bins for 1D
      Int_t n_3d_fine ;            // # bins for 3D, Fine Mesh
      Int_t n_3d_coarse;           // # bins for 3D, Coarse Mesh
      Int_t n_3d_total;            // Total # bins for 3D
      Int_t n_3d_fine_project;     // # 3D fine mesh bins to sum over for

//   Momentum Space Sectors for Track Sorting:

      Int_t n_px_bins;             // # sector bins in px
      Int_t n_py_bins;             // # sector bins in py
      Int_t n_pz_bins;             // # sector bins in pz
      Int_t n_sectors;             // Total # sectors in 3D momentum space

//   Temporary Momentum Space Sector information storage during trk adjust.

      Int_t old_sec_ntrk;          // Old sector # tracks
      Int_t old_sec_flag;         // Old sector flag value
      Int_t old_sec_trkid[MAX_TRK_SAVE];         // Old sector track id array

      Int_t new_sec_ntrk;          // New sector # tracks
      Int_t new_sec_flag;          // New sector flag value
      Int_t new_sec_trkid[MAX_TRK_SAVE];         // New sector track id array
      Int_t new_sec_save;          // New sector ID value
      Int_t old_sec_save;          // Old sector ID value
     
      Float_t    pt_bin_size ;          // One-body pt bin size in (GeV/c)
      Float_t    pt_min;       // One-body pt min/max in (GeV/c)
      Float_t    pt_max;
      
      Float_t    phi_bin_size;          // One-body phi bin size in (degrees)
      Float_t    phi_min;      // One-body phi min/max in (degrees)
      Float_t    phi_max;
      
      Float_t    eta_bin_size ;         // One-body eta bin size
      Float_t    eta_min;       // One-body eta min/max
      Float_t    eta_max;
//   Two-Body Histograms and Correlation Mesh for 1D and 3D distributions:

//                                     // projections onto single axis.

      Float_t    binsize_1d_fine;       // Bin Size - 1D, Fine Mesh in (GeV/c)
      Float_t    binsize_1d_coarse;     // Bin Size - 1D, Coarse Mesh in (GeV/c)
      Float_t    qmid_1d;               // q (GeV/c) at fine-coarse mesh boundary
      Float_t    qmax_1d;               // Max q (GeV/c) for 1D distributions
      Float_t    binsize_3d_fine;       // Bin Size - 3D, Fine Mesh in (GeV/c)
      Float_t    binsize_3d_coarse;     // Bin Size - 3D, Coarse Mesh in (GeV/c)
      Float_t    qmid_3d;               // q (GeV/c) at fine-coarse mesh boundary
      Float_t    qmax_3d;               // Max q (GeV/c) for 3D distributions




      Float_t    px_min;        // Sector range in px in GeV/c
      Float_t    px_max;
      Float_t    delpx;                 // Mom. space sector cell size - px(GeV/c)     
      
      Float_t    py_min;        // Sector range in py in GeV/c
      Float_t    py_max;
      Float_t    delpy;                 // Mom. space sector cell size - py(GeV/c)     

      Float_t    pz_min;        // Sector range in pz in GeV/c
      Float_t    pz_max;
      Float_t    delpz;                 // Mom. space sector cell size - pz(GeV/c)

    
  }HBTprocMeshCommon;

#define MESH COMMON_BLOCK(MESH,mesh)
COMMON_BLOCK_DEF(HBTprocMeshCommon, MESH);
//HBTprocMeshCommon MESH;

/************************************************************************************************/ 
/************************************************************************************************/ 

typedef struct //TRACK1 COMMON
  {
      Int_t   trk_id[TRK_MAXLEN];         // Track ID number
      Int_t   trk_px_sec[TRK_MAXLEN];      // px sector number
      Int_t   trk_py_sec[TRK_MAXLEN];      // py sector number
      Int_t   trk_pz_sec[TRK_MAXLEN];      // pz sector number
      Int_t   trk_sector[TRK_MAXLEN];      // unique sector ID number
      Int_t   trk_flag[TRK_MAXLEN];        // normally=0,if 1 indicates track assigned
//                               // to sector with too many tracks, if = 1
//                               // then track is not used.  See /sec_trk_map/ 
      Int_t   trk_out_flag[TRK_MAXLEN];    // flag indicating track in/out of accept.
//                               // non-zero for track pushed out of accept.
      Int_t   trk_merge_flag[TRK_MAXLEN];  // flag indicating track is merged (not used)
      Int_t   trk_ge_pid[TRK_MAXLEN];      // Geant particle ID code number
      Int_t   trk_start_vertex[TRK_MAXLEN]; // From input event file - track's start vrtx
      Int_t   trk_stop_vertex[TRK_MAXLEN];  // From input event file - track's stop vrtx
      Int_t   trk_event_line[TRK_MAXLEN];  // Line # of track in input event text file
      
      Float_t      trk_px[TRK_MAXLEN];          // x component of track momentum in GeV/c
      Float_t      trk_py[TRK_MAXLEN];          // y component of track momentum in GeV/c
      Float_t      trk_pz[TRK_MAXLEN];          // z component of track momentum in GeV/c
      Float_t      trk_E[TRK_MAXLEN];           // Total energy of track in GeV
      Float_t      trk_pt[TRK_MAXLEN];          // pt of track momentum in GeV/c
      Float_t      trk_phi[TRK_MAXLEN];         // azimuthal angle of track in degrees 
      Float_t      trk_eta[TRK_MAXLEN];         // pseudorapidity of track

    
  }HBTprocTrack1Common;
  
#define TRACK1 COMMON_BLOCK(TRACK1,track1)
COMMON_BLOCK_DEF(HBTprocTrack1Common, TRACK1);

/************************************************************************************************/ 
/************************************************************************************************/ 

  typedef struct //TRACK2  COMMON
   {
      Int_t   trk2_id[TRK2_MAXLEN];          // Track ID number
      Int_t   trk2_px_sec[TRK2_MAXLEN];      // px sector number
      Int_t   trk2_py_sec[TRK2_MAXLEN];      // py sector number
      Int_t   trk2_pz_sec [TRK2_MAXLEN];     // pz sector number
      Int_t   trk2_sector[TRK2_MAXLEN];      // unique sector ID number
      Int_t   trk2_flag[TRK2_MAXLEN];        // normally=0,if 1 indicates track assigned
//                               // to sector with too many tracks, if = 1
//                               // then track is not used.  See /sec_trk2_map/ 
      Int_t   trk2_out_flag[TRK2_MAXLEN];    // flag indicating track in/out of accept.
//                               // non-zero for track pushed out of accept.
      Int_t   trk2_merge_flag[TRK2_MAXLEN];  // flag indicating track is merged (not used)
      Int_t   trk2_ge_pid[TRK2_MAXLEN];      // Geant particle ID code number
      Int_t   trk2_start_vertex[TRK2_MAXLEN]; // From input event file - track's start vrtx
      Int_t   trk2_stop_vertex[TRK2_MAXLEN];  // From input event file - track's stop vrtx
      Int_t   trk2_event_line[TRK2_MAXLEN];  // Line # of track in input event text file
      
      Float_t      trk2_px[TRK2_MAXLEN];          // x component of track momentum in GeV/c
      Float_t      trk2_py[TRK2_MAXLEN];          // y component of track momentum in GeV/c
      Float_t      trk2_pz[TRK2_MAXLEN];          // z component of track momentum in GeV/c
      Float_t      trk2_E[TRK2_MAXLEN];           // Total energy of track in GeV
      Float_t      trk2_pt[TRK2_MAXLEN];          // pt of track momentum in GeV/c
      Float_t      trk2_phi[TRK2_MAXLEN];         // azimuthal angle of track in degrees 
      Float_t      trk2_eta [TRK2_MAXLEN];        // pseudorapidity of track

   }HBTprocTrack2Common;

#define TRACK2 COMMON_BLOCK(TRACK2,track2)
COMMON_BLOCK_DEF(HBTprocTrack2Common, TRACK2);
/************************************************************************************************/ 
/************************************************************************************************/ 
typedef struct //SEC_TRK_MAP COMMON BLOCK
  {
    
      Int_t stm_sec_id [SEC_MAXLEN];       //  unique sector ID number
      Int_t stm_n_trk_sec[SEC_MAXLEN];     //  Number of tracks assigned to sector
      Int_t stm_flag[SEC_MAXLEN];          //  normally=0, if = 1 then more than
//                                         //  max_trk_sec tracks could have been 
//                                         //  assigned to this sector, however the
//                                         //  maximum number that can be assigned is
//                                         //  max_trk_sec.
      Int_t stm_track_id[SEC_MAXLEN][MAX_TRK_SEC];      //  Foreign keys to tracks in /track/ that
//                                                      //  are assigned to this sector.

  }HBTprocSecTrackMapCommon;

#define SEC_TRK_MAP COMMON_BLOCK(SEC_TRK_MAP,sec_trk_map)
COMMON_BLOCK_DEF(HBTprocSecTrackMapCommon, SEC_TRK_MAP);

/************************************************************************************************/ 
/************************************************************************************************/ 
typedef struct //SEC_TRK_MAP2 COMMON BLOCK
  {
    
      Int_t stm_sec_id2 [SEC_MAXLEN2];       //  unique sector ID number
      Int_t stm_n_trk_sec2[SEC_MAXLEN2];     //  Number of tracks assigned to sector
      Int_t stm_flag2[SEC_MAXLEN2];          //  normally=0, if = 1 then more than
//                                         //  max_trk_sec tracks could have been 
//                                         //  assigned to this sector, however the
//                                         //  maximum number that can be assigned is
//                                         //  max_trk_sec.
      Int_t stm_track_id2[SEC_MAXLEN2][MAX_TRK_SEC2];      //  Foreign keys to tracks in /track/ that
//                                                      //  are assigned to this sector.

  }HBTprocSecTrackMap2Common;

#define SEC_TRK_MAP2 COMMON_BLOCK(SEC_TRK_MAP2,sec_trk_map2)
COMMON_BLOCK_DEF(HBTprocSecTrackMap2Common, SEC_TRK_MAP2);

/************************************************************************************************/ 
/************************************************************************************************/ 

typedef struct //PARTICLE COMMON
 {
//   Variable Type Declarations:
 
      Int_t         part_id[PART_MAXLEN];        // Geant particle ID code number; required
//                                  // to be equal to the row number
      Int_t         part_charge[PART_MAXLEN];    // Electric charge in units of |e|
      Float_t       part_mass[PART_MAXLEN];      // Rest mass in GeV/c**2
      Float_t       part_lifetime[PART_MAXLEN];  // Proper lifetime in sec.    
 }HBTprocParticleCommon;

#define PARTICLE COMMON_BLOCK(PARTICLE,particle)
COMMON_BLOCK_DEF(HBTprocParticleCommon, PARTICLE);

/************************************************************************************************/ 
/************************************************************************************************/ 
typedef struct //CORRELATIONS COMMON
 {
//   One-dimensional Functions:

      Float_t c2mod_like_1d[MAX_C2_1D]; 
      Float_t c2mod_unlike_1d[MAX_C2_1D];
      Float_t c2fit_like_1d[MAX_C2_1D];
      Float_t c2fit_unlike_1d[MAX_C2_1D];
      Float_t c2err_like_1d[MAX_C2_1D];
      Float_t c2err_unlike_1d[MAX_C2_1D];

//   Three-dimensional Functions:

      Float_t c2mod_like_3d_fine[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2mod_unlike_3d_fine[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2mod_like_3d_coarse[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2mod_unlike_3d_coarse[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2fit_like_3d_fine[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2fit_unlike_3d_fine[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2fit_like_3d_coarse[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2fit_unlike_3d_coarse[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2err_like_3d_fine[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2err_unlike_3d_fine[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2err_like_3d_coarse[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];
      Float_t c2err_unlike_3d_coarse[MAX_C2_3D][MAX_C2_3D][MAX_C2_3D];

  
 }HBTprocCorrelCommon;

#define CORRELATIONS COMMON_BLOCK(CORRELATIONS,correlations)
COMMON_BLOCK_DEF(HBTprocCorrelCommon, CORRELATIONS);

/************************************************************************************************/ 
/************************************************************************************************/ 
 typedef struct //COULMB
   {
      Float_t c2_coul_like[MAX_C2_COUL];
      Float_t c2_coul_unlike[MAX_C2_COUL];
      Float_t q_coul[MAX_C2_COUL]; 
   }HBTprocCoulumbCommon;
   
#define COULMB COMMON_BLOCK(COULMB,coulmb)

COMMON_BLOCK_DEF(HBTprocCoulumbCommon,COULMB);
   
/************************************************************************************************/ 
/************************************************************************************************/ 

typedef struct 
 {
   Float_t num_iter[MAX_EVENTS]; 
   Float_t niter_mean; 
   Float_t niter_rms;
   Float_t n_part_used_1_store[MAX_EVENTS]; 
   Float_t npart1_mean; 
   Float_t npart1_rms;
   Float_t n_part_used_2_store[MAX_EVENTS]; 
   Float_t npart2_mean; 
   Float_t npart2_rms;
   Float_t n_part_tot_store[MAX_EVENTS]; 
   Float_t npart_tot_mean; 
   Float_t npart_tot_rms;
   Float_t num_sec_flagged_store[MAX_EVENTS]; 
   Float_t nsec_flag_mean; 
   Float_t nsec_flag_rms;
   Float_t frac_trks_out[MAX_EVENTS]; 
   Float_t frac_trks_out_mean;
   Float_t frac_trks_out_rms;
   Float_t frac_trks_flag[MAX_EVENTS];
   Float_t frac_trks_flag_mean;
   Float_t frac_trks_flag_rms;
   Float_t chisq_like_1d_store[MAX_EVENTS];         
   Float_t chi_l1d_mean;
   Float_t chi_l1d_rms;
   Float_t chisq_unlike_1d_store[MAX_EVENTS];       
   Float_t chi_u1d_mean;
   Float_t chi_u1d_rms;
                                                                         
  }HBTprocEventSummaryCommon;

#define EVENT_SUMMARY COMMON_BLOCK(EVENT_SUMMARY,event_summary)

COMMON_BLOCK_DEF(HBTprocEventSummaryCommon,EVENT_SUMMARY);
/************************************************************************************************/ 
/************************************************************************************************/ 
typedef struct//      common/histograms/
  {
      Int_t hist1_pt_1[MAX_H_1D];
      Int_t hist1_phi_1[MAX_H_1D];
      Int_t hist1_eta_1[MAX_H_1D];
      Int_t hist1_pt_2[MAX_H_1D];
      Int_t hist1_phi_2[MAX_H_1D];
      Int_t hist1_eta_2[MAX_H_1D];
      Int_t htmp1_pt_1[MAX_H_1D];
      Int_t htmp1_phi_1[MAX_H_1D];
      Int_t htmp1_eta_1[MAX_H_1D];
      Int_t htmp1_pt_2[MAX_H_1D];
      Int_t htmp1_phi_2[MAX_H_1D];
      Int_t htmp1_eta_2[MAX_H_1D];
      Int_t href1_pt_1[MAX_H_1D];
      Int_t href1_phi_1[MAX_H_1D];
      Int_t href1_eta_1[MAX_H_1D];
      Int_t href1_pt_2[MAX_H_1D];
      Int_t href1_phi_2[MAX_H_1D];
      Int_t href1_eta_2[MAX_H_1D];
      Int_t hinc1_pt_1[MAX_H_1D];
      Int_t hinc1_phi_1[MAX_H_1D];
      Int_t hinc1_eta_1[MAX_H_1D];
      Int_t hinc1_pt_2[MAX_H_1D];
      Int_t hinc1_phi_2[MAX_H_1D];
      Int_t hinc1_eta_2[MAX_H_1D];
      Int_t hist_like_1d[MAX_H_1D];
      Int_t hist_unlike_1d[MAX_H_1D];
      Int_t htmp_like_1d[MAX_H_1D];
      Int_t htmp_unlike_1d[MAX_H_1D];
      Int_t href_like_1d[MAX_H_1D];
      Int_t href_unlike_1d[MAX_H_1D];
      Int_t hinc_like_1d[MAX_H_1D];
      Int_t hinc_unlike_1d[MAX_H_1D];
      
      Int_t hist_like_3d_fine[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t hist_unlike_3d_fine[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t hist_like_3d_coarse[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t hist_unlike_3d_coarse[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t htmp_like_3d_fine[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t htmp_unlike_3d_fine[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t htmp_like_3d_coarse[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t htmp_unlike_3d_coarse[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t href_like_3d_fine[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t href_unlike_3d_fine[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t href_like_3d_coarse[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t href_unlike_3d_coarse[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t hinc_like_3d_fine[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t hinc_unlike_3d_fine[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t hinc_like_3d_coarse[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      Int_t hinc_unlike_3d_coarse[MAX_H_3D][MAX_H_3D][MAX_H_3D];
      
  }HBTprocHistosCommon;
  
#define HISTOGRAMS COMMON_BLOCK(HISTOGRAMS,histograms)

COMMON_BLOCK_DEF(HBTprocHistosCommon,HISTOGRAMS);

/************************************************************************************************/ 
/************************************************************************************************/ 
 
#endif

}

 
#endif                     
