CCC  Common for bin-by-bin distribution input:
CCC  NOTE:  Include file 'Parameter_values.inc' must accompany and 
CCC         precede this file everywhere it occurs.

      integer n_bins_max
      parameter(n_bins_max = 50) ! maximum # of input pt, eta bins

      Common/dist_bin/ pt_start(npid),eta_start(npid),pt_stop(npid),
     1    eta_stop(npid),delta_pt(npid,n_bins_max),
     2    delta_eta(npid,n_bins_max),pt_bin(npid,n_bins_max),
     3    eta_bin(npid,n_bins_max),
     4    pt_eta_bin(npid,n_bins_max,n_bins_max),
     5    pt_bin_mesh(npid,n_bins_max),eta_bin_mesh(npid,n_bins_max),
     6    n_pt_bins(npid),n_eta_bins(npid) 

      integer n_pt_bins,n_eta_bins
      real*4 pt_start,eta_start,pt_stop,eta_stop
      real*4 delta_pt,delta_eta,pt_bin,eta_bin,pt_eta_bin
      real*4 pt_bin_mesh,eta_bin_mesh
