CCC   common/mesh/
C
C     This common block includes the following data:
C
C          o Mesh Bin Definitions for the following:
C            i) One-Body distributions in {pt,phi,eta}
C           ii) Two-Body Histograms and Correlations for 1D and 3D
C          iii) Momentum space sectors for track sorting in {px,py,pz}

      integer*4 max_trk_save
      parameter (max_trk_save = 150)

      common/mesh/
     1   n_pt_bins,n_phi_bins,n_eta_bins,
     2   n_1d_fine,n_1d_coarse,n_1d_total,
     3   n_3d_fine,n_3d_coarse,n_3d_total,n_3d_fine_project,
     4   n_px_bins,n_py_bins,n_pz_bins,n_sectors,
     5   old_sec_ntrk,old_sec_flag,old_sec_trkid(max_trk_save),
     6   new_sec_ntrk,new_sec_flag,new_sec_trkid(max_trk_save),
     7   new_sec_save,old_sec_save,
     8   pt_bin_size,pt_min,pt_max,
     9   phi_bin_size,phi_min,phi_max,
     1   eta_bin_size,eta_min,eta_max,
     2   binsize_1d_fine,binsize_1d_coarse,qmid_1d,qmax_1d,
     3   binsize_3d_fine,binsize_3d_coarse,qmid_3d,qmax_3d,
     4   px_min,px_max,delpx,
     5   py_min,py_max,delpy,
     6   pz_min,pz_max,delpz

CCC   Variable Type Declarations:

CCC   One-Body Mesh:

      integer*4 n_pt_bins             ! # one-body pt bins
      integer*4 n_phi_bins            ! # one-body phi bins
      integer*4 n_eta_bins            ! # one-body eta bins
      real*4    pt_bin_size           ! One-body pt bin size in (GeV/c)
      real*4    phi_bin_size          ! One-body phi bin size in (degrees)
      real*4    eta_bin_size          ! One-body eta bin size
      real*4    pt_min, pt_max        ! One-body pt min/max in (GeV/c)
      real*4    phi_min, phi_max      ! One-body phi min/max in (degrees)
      real*4    eta_min, eta_max      ! One-body eta min/max

CCC   Two-Body Histograms and Correlation Mesh for 1D and 3D distributions:

      integer*4 n_1d_fine             ! # bins for 1D, Fine Mesh
      integer*4 n_1d_coarse           ! # bins for 1D, Coarse Mesh
      integer*4 n_1d_total            ! Total # bins for 1D
      integer*4 n_3d_fine             ! # bins for 3D, Fine Mesh
      integer*4 n_3d_coarse           ! # bins for 3D, Coarse Mesh
      integer*4 n_3d_total            ! Total # bins for 3D
      integer*4 n_3d_fine_project     ! # 3D fine mesh bins to sum over for
C                                     ! projections onto single axis.

      real*4    binsize_1d_fine       ! Bin Size - 1D, Fine Mesh in (GeV/c)
      real*4    binsize_1d_coarse     ! Bin Size - 1D, Coarse Mesh in (GeV/c)
      real*4    qmid_1d               ! q (GeV/c) at fine-coarse mesh boundary
      real*4    qmax_1d               ! Max q (GeV/c) for 1D distributions
      real*4    binsize_3d_fine       ! Bin Size - 3D, Fine Mesh in (GeV/c)
      real*4    binsize_3d_coarse     ! Bin Size - 3D, Coarse Mesh in (GeV/c)
      real*4    qmid_3d               ! q (GeV/c) at fine-coarse mesh boundary
      real*4    qmax_3d               ! Max q (GeV/c) for 3D distributions

CCC   Momentum Space Sectors for Track Sorting:

      integer*4 n_px_bins             ! # sector bins in px
      integer*4 n_py_bins             ! # sector bins in py
      integer*4 n_pz_bins             ! # sector bins in pz
      integer*4 n_sectors             ! Total # sectors in 3D momentum space

      real*4    delpx                 ! Mom. space sector cell size - px(GeV/c)
      real*4    delpy                 ! Mom. space sector cell size - py(GeV/c)
      real*4    delpz                 ! Mom. space sector cell size - pz(GeV/c)
      real*4    px_min, px_max        ! Sector range in px in GeV/c
      real*4    py_min, py_max        ! Sector range in py in GeV/c
      real*4    pz_min, pz_max        ! Sector range in pz in GeV/c

CCC   Temporary Momentum Space Sector information storage during trk adjust.

      integer*4 old_sec_ntrk          ! Old sector # tracks
      integer*4 old_sec_flag          ! Old sector flag value
      integer*4 old_sec_trkid         ! Old sector track id array
      integer*4 old_sec_save          ! Old sector ID value
      integer*4 new_sec_ntrk          ! New sector # tracks
      integer*4 new_sec_flag          ! New sector flag value
      integer*4 new_sec_trkid         ! New sector track id array
      integer*4 new_sec_save          ! New sector ID value

CCC   End common/mesh/
