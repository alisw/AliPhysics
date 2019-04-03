#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TBranch.h"
#include "TIterator.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TError.h"
#include "AliTRDtrendingManager.h"

#endif

void makeTrendingDB(const Char_t *fl)
{
// Make trending of variable list "tl" from trending file list "fl"
// The trending value list should be formated "var1:var2:var3"
// The trending file from the list should be found on a path formated "your path"/runId/TRD.PerformanceTrend.root 
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWGmuon");

  const Int_t nt(322);
  const Char_t *tvn[nt][2] = {
    {"TRDcheckDET_NTracksEvent", "<N_{track}/Event>"},
    {"TRDcheckDET_NTracksEventRMS", "RMS(N_{track}/Event)"},
    {"TRDcheckDET_NTracksSector", "<N_{track}/Sector>"},
    {"TRDcheckDET_NClustersTrack", "<N_{cls}/Track>"},
    {"TRDcheckDET_NClustersTrackRMS", "RMS(N_{cls}/Track)"},
    {"TRDcheckDET_NClustersTracklet", "<N_{cls}/Tracklet>"},
    {"TRDcheckDET_NClustersTrackletRMS", "RMS(N_{cls}/Tracklet)"},
    {"TRDcheckDET_NTrackletsTrack", "<N_{tracklet}/Track>"},
    {"TRDcheckDET_NTrackletsTrackRMS", "RMS(N_{tracklet}/Track)"},
    {"TRDcheckDET_ChargeTracklet", "<dQdl>"},
    {"TRDcheckDET_ChargeTrackletRMS", "RMS(dQdl)"},
    {"TRDcheckDET_PHplateau", "Plateau(<PH>)"},
    {"TRDcheckDET_PHslope", "Slope(<PH>)"},
    {"TRDcheckDET_PHamplificationPeak", "Peak(<PH>)"},
//=======================================================
    {"TRDefficiency_pt", "TPC-TRD efficiency"},
    {"TRDefficiency_pt0", "Eff.| #it{p_{t}[GeV/c]<0.5}"},
    {"TRDefficiency_pt1", "Eff.| #it{0.5<=p_{t}[GeV/c]<0.8}"},
    {"TRDefficiency_pt2", "Eff.| #it{0.8<=p_{t}[GeV/c]<1.5}"},
    {"TRDefficiency_pt3", "Eff.| #it{1.5<=p_{t}[GeV/c]<5.0}"},
    {"TRDefficiency_pt4", "Eff.| #it{5.0<=p_{t}[GeV/c]}"},
//=======================================================
    {"TRDresolution_ClS0", "Cl :: #sigma(#Deltay) [cm] | Ly0"},
    {"TRDresolution_ClS1", "Cl :: #sigma(#Deltay) [cm] | Ly1"},
    {"TRDresolution_ClS2", "Cl :: #sigma(#Deltay) [cm] | Ly2"},
    {"TRDresolution_ClS3", "Cl :: #sigma(#Deltay) [cm] | Ly3"},
    {"TRDresolution_ClS4", "Cl :: #sigma(#Deltay) [cm] | Ly4"},
    {"TRDresolution_ClS5", "Cl :: #sigma(#Deltay) [cm] | Ly5"},
//=======================================================
    {"TRDresolution_TrkltY0", "Trklt :: <#Deltay> [cm] | Ly0"},
    {"TRDresolution_TrkltY1", "Trklt :: <#Deltay> [cm] | Ly1"},
    {"TRDresolution_TrkltY2", "Trklt :: <#Deltay> [cm] | Ly2"},
    {"TRDresolution_TrkltY3", "Trklt :: <#Deltay> [cm] | Ly3"},
    {"TRDresolution_TrkltY4", "Trklt :: <#Deltay> [cm] | Ly4"},
    {"TRDresolution_TrkltY5", "Trklt :: <#Deltay> [cm] | Ly5"},
    {"TRDresolution_TrkltYS0", "Trklt :: #sigma(#Deltay) [cm] | Ly0"},
    {"TRDresolution_TrkltYS1", "Trklt :: #sigma(#Deltay) [cm] | Ly1"},
    {"TRDresolution_TrkltYS2", "Trklt :: #sigma(#Deltay) [cm] | Ly2"},
    {"TRDresolution_TrkltYS3", "Trklt :: #sigma(#Deltay) [cm] | Ly3"},
    {"TRDresolution_TrkltYS4", "Trklt :: #sigma(#Deltay) [cm] | Ly4"},
    {"TRDresolution_TrkltYS5", "Trklt :: #sigma(#Deltay) [cm] | Ly5"},
//=======================================================
    {"TRDresolution_TrkInYnL0", "TrkIn :: <#Deltay>^{e-}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYnl0", "TrkIn :: <#Deltay>^{e-}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYni0", "TrkIn :: <#Deltay>^{e-}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYnh0", "TrkIn :: <#Deltay>^{e-}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYnH0", "TrkIn :: <#Deltay>^{e-}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYn0", "TrkIn :: <#Deltay>^{e-} [cm]"},
    {"TRDresolution_TrkInYSnL0", "TrkIn :: RMS(#Deltay)^{e-}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYSnl0", "TrkIn :: RMS(#Deltay)^{e-}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYSni0", "TrkIn :: RMS(#Deltay)^{e-}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYSnh0", "TrkIn :: RMS(#Deltay)^{e-}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYSnH0", "TrkIn :: RMS(#Deltay)^{e-}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYSn0", "TrkIn :: RMS(#Deltay)^{e-} [cm]"},
    {"TRDresolution_TrkInPhnL0", "TrkIn :: <#Delta#phi>^{e-}| #it{p_{t}[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhnl0", "TrkIn :: <#Delta#phi>^{e-}| #it{0.5<=p_{t}[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhni0", "TrkIn :: <#Delta#phi>^{e-}| #it{0.8<=p_{t}[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhnh0", "TrkIn :: <#Delta#phi>^{e-}| #it{1.5<=p_{t}[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhnH0", "TrkIn :: <#Delta#phi>^{e-}| #it{5.0<=p_{t}[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhn0", "TrkIn :: <#Delta#phi>^{e-} [deg]"},
    {"TRDresolution_TrkInPhSnL0", "TrkIn :: RMS<#Delta#phi>^{e-}| #it{p[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhSnl0", "TrkIn :: RMS<#Delta#phi>^{e-}| #it{0.5<=p[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhSni0", "TrkIn :: RMS<#Delta#phi>^{e-}| #it{0.8<=p[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhSnh0", "TrkIn :: RMS<#Delta#phi>^{e-}| #it{1.5<=p[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhSnH0", "TrkIn :: RMS<#Delta#phi>^{e-}| #it{5.0<=p[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhSn0", "TrkIn :: RMS<#Delta#phi>^{e-} [deg]"},
    {"TRDresolution_TrkInQnL0", "TrkIn :: MPV(dQdl)^{e-}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQnl0", "TrkIn :: MPV(dQdl)^{e-}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQni0", "TrkIn :: MPV(dQdl)^{e-}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQnh0", "TrkIn :: MPV(dQdl)^{e-}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQnH0", "TrkIn :: MPV(dQdl)^{e-}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQn0", "TrkIn :: MPV(dQdl)^{e-} [a.u.]"},
    {"TRDresolution_TrkInQSnL0", "TrkIn :: <dQdl>^{e-}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQSnl0", "TrkIn :: <dQdl>^{e-}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQSni0", "TrkIn :: <dQdl>^{e-}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQSnh0", "TrkIn :: <dQdl>^{e-}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQSnH0", "TrkIn :: <dQdl>^{e-}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQSn0", "TrkIn :: <dQdl>^{e-} [a.u.]"},
    {"TRDresolution_TrkInYnL1", "TrkIn :: <#Deltay>^{#mu#pi-}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYnl1", "TrkIn :: <#Deltay>^{#mu#pi-}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYni1", "TrkIn :: <#Deltay>^{#mu#pi-}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYnh1", "TrkIn :: <#Deltay>^{#mu#pi-}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYnH1", "TrkIn :: <#Deltay>^{#mu#pi-}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYn1", "TrkIn :: <#Deltay>^{#mu#pi-} [cm]"},
    {"TRDresolution_TrkInYSnL1", "TrkIn :: RMS(#Deltay)^{#mu#pi-}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYSnl1", "TrkIn :: RMS(#Deltay)^{#mu#pi-}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYSni1", "TrkIn :: RMS(#Deltay)^{#mu#pi-}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYSnh1", "TrkIn :: RMS(#Deltay)^{#mu#pi-}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYSnH1", "TrkIn :: RMS(#Deltay)^{#mu#pi-}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYSn1", "TrkIn :: RMS(#Deltay)^{#mu#pi-} [cm]"},
    {"TRDresolution_TrkInPhnL1", "TrkIn :: <#Delta#phi>^{#mu#pi-}| #it{p_{t}[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhnl1", "TrkIn :: <#Delta#phi>^{#mu#pi-}| #it{0.5<=p_{t}[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhni1", "TrkIn :: <#Delta#phi>^{#mu#pi-}| #it{0.8<=p_{t}[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhnh1", "TrkIn :: <#Delta#phi>^{#mu#pi-}| #it{1.5<=p_{t}[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhnH1", "TrkIn :: <#Delta#phi>^{#mu#pi-}| #it{5.0<=p_{t}[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhn1", "TrkIn :: <#Delta#phi>^{#mu#pi-} [deg]"},
    {"TRDresolution_TrkInPhSnL1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi-}| #it{p[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhSnl1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi-}| #it{0.5<=p[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhSni1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi-}| #it{0.8<=p[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhSnh1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi-}| #it{1.5<=p[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhSnH1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi-}| #it{5.0<=p[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhSn1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi-} [deg]"},
    {"TRDresolution_TrkInQnL1", "TrkIn :: MPV(dQdl)^{#mu#pi-}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQnl1", "TrkIn :: MPV(dQdl)^{#mu#pi-}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQni1", "TrkIn :: MPV(dQdl)^{#mu#pi-}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQnh1", "TrkIn :: MPV(dQdl)^{#mu#pi-}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQnH1", "TrkIn :: MPV(dQdl)^{#mu#pi-}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQn1", "TrkIn :: MPV(dQdl)^{#mu#pi-} [a.u.]"},
    {"TRDresolution_TrkInQSnL1", "TrkIn :: <dQdl>^{#mu#pi-}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQSnl1", "TrkIn :: <dQdl>^{#mu#pi-}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQSni1", "TrkIn :: <dQdl>^{#mu#pi-}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQSnh1", "TrkIn :: <dQdl>^{#mu#pi-}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQSnH1", "TrkIn :: <dQdl>^{#mu#pi-}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQSn1", "TrkIn :: <dQdl>^{#mu#pi-} [a.u.]"},
    {"TRDresolution_TrkInYnL2", "TrkIn :: <#Deltay>^{Kp-}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYnl2", "TrkIn :: <#Deltay>^{Kp-}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYni2", "TrkIn :: <#Deltay>^{Kp-}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYnh2", "TrkIn :: <#Deltay>^{Kp-}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYnH2", "TrkIn :: <#Deltay>^{Kp-}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYn2", "TrkIn :: <#Deltay>^{Kp-} [cm]"},
    {"TRDresolution_TrkInYSnL2", "TrkIn :: RMS(#Deltay)^{Kp-}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYSnl2", "TrkIn :: RMS(#Deltay)^{Kp-}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYSni2", "TrkIn :: RMS(#Deltay)^{Kp-}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYSnh2", "TrkIn :: RMS(#Deltay)^{Kp-}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYSnH2", "TrkIn :: RMS(#Deltay)^{Kp-}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYSn2", "TrkIn :: RMS(#Deltay)^{Kp-} [cm]"},
    {"TRDresolution_TrkInPhnL2", "TrkIn :: <#Delta#phi>^{Kp-}| #it{p_{t}[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhnl2", "TrkIn :: <#Delta#phi>^{Kp-}| #it{0.5<=p_{t}[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhni2", "TrkIn :: <#Delta#phi>^{Kp-}| #it{0.8<=p_{t}[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhnh2", "TrkIn :: <#Delta#phi>^{Kp-}| #it{1.5<=p_{t}[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhnH2", "TrkIn :: <#Delta#phi>^{Kp-}| #it{5.0<=p_{t}[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhn2", "TrkIn :: <#Delta#phi>^{Kp-} [deg]"},
    {"TRDresolution_TrkInPhSnL2", "TrkIn :: RMS<#Delta#phi>^{Kp-}| #it{p[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhSnl2", "TrkIn :: RMS<#Delta#phi>^{Kp-}| #it{0.5<=p[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhSni2", "TrkIn :: RMS<#Delta#phi>^{Kp-}| #it{0.8<=p[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhSnh2", "TrkIn :: RMS<#Delta#phi>^{Kp-}| #it{1.5<=p[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhSnH2", "TrkIn :: RMS<#Delta#phi>^{Kp-}| #it{5.0<=p[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhSn2", "TrkIn :: RMS<#Delta#phi>^{Kp-} [deg]"},
    {"TRDresolution_TrkInQnL2", "TrkIn :: MPV(dQdl)^{Kp-}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQnl2", "TrkIn :: MPV(dQdl)^{Kp-}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQni2", "TrkIn :: MPV(dQdl)^{Kp-}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQnh2", "TrkIn :: MPV(dQdl)^{Kp-}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQnH2", "TrkIn :: MPV(dQdl)^{Kp-}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQn2", "TrkIn :: MPV(dQdl)^{Kp-} [a.u.]"},
    {"TRDresolution_TrkInQSnL2", "TrkIn :: <dQdl>^{Kp-}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQSnl2", "TrkIn :: <dQdl>^{Kp-}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQSni2", "TrkIn :: <dQdl>^{Kp-}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQSnh2", "TrkIn :: <dQdl>^{Kp-}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQSnH2", "TrkIn :: <dQdl>^{Kp-}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQSn2", "TrkIn :: <dQdl>^{Kp-} [a.u.]"},
    {"TRDresolution_TrkInYnL", "TrkIn :: <#Deltay>^{-}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYSnL", "TrkIn :: RMS(#Deltay)^{-}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInPhnL", "TrkIn :: <#Delta#phi>^{-}| #it{p_{t}[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhSnL", "TrkIn :: RMS<#Delta#phi>^{-}| #it{p_{t}[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInYnl", "TrkIn :: <#Deltay>^{-}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYSnl", "TrkIn :: RMS(#Deltay)^{-}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInPhnl", "TrkIn :: <#Delta#phi>^{-}| #it{0.5<=p_{t}[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhSnl", "TrkIn :: RMS<#Delta#phi>^{-}| #it{0.5<=p_{t}[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInYni", "TrkIn :: <#Deltay>^{-}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYSni", "TrkIn :: RMS(#Deltay)^{-}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInPhni", "TrkIn :: <#Delta#phi>^{-}| #it{0.8<=p_{t}[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhSni", "TrkIn :: RMS<#Delta#phi>^{-}| #it{0.8<=p_{t}[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInYnh", "TrkIn :: <#Deltay>^{-}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYSnh", "TrkIn :: RMS(#Deltay)^{-}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInPhnh", "TrkIn :: <#Delta#phi>^{-}| #it{1.5<=p_{t}[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhSnh", "TrkIn :: RMS<#Delta#phi>^{-}| #it{1.5<=p_{t}[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInYnH", "TrkIn :: <#Deltay>^{-}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYSnH", "TrkIn :: RMS(#Deltay)^{-}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInPhnH", "TrkIn :: <#Delta#phi>^{-}| #it{5.0<=p_{t}[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhSnH", "TrkIn :: RMS<#Delta#phi>^{-}| #it{5.0<=p_{t}[GeV/c]} [deg]"},
    {"TRDresolution_TrkInYn", "TrkIn :: <#Deltay>^{-} [cm]"},
    {"TRDresolution_TrkInYSn", "TrkIn :: RMS(#Deltay)^{-} [cm]"},
    {"TRDresolution_TrkInPhn", "TrkIn :: <#Delta#phi>^{-} [deg]"},
    {"TRDresolution_TrkInPhSn", "TrkIn :: RMS<#Delta#phi>^{-} [deg]"},
    {"TRDresolution_TrkInYpL0", "TrkIn :: <#Deltay>^{e+}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYpl0", "TrkIn :: <#Deltay>^{e+}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYpi0", "TrkIn :: <#Deltay>^{e+}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYph0", "TrkIn :: <#Deltay>^{e+}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYpH0", "TrkIn :: <#Deltay>^{e+}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYp0", "TrkIn :: <#Deltay>^{e+} [cm]"},
    {"TRDresolution_TrkInYSpL0", "TrkIn :: RMS(#Deltay)^{e+}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYSpl0", "TrkIn :: RMS(#Deltay)^{e+}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYSpi0", "TrkIn :: RMS(#Deltay)^{e+}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYSph0", "TrkIn :: RMS(#Deltay)^{e+}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYSpH0", "TrkIn :: RMS(#Deltay)^{e+}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYSp0", "TrkIn :: RMS(#Deltay)^{e+} [cm]"},
    {"TRDresolution_TrkInPhpL0", "TrkIn :: <#Delta#phi>^{e+}| #it{p_{t}[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhpl0", "TrkIn :: <#Delta#phi>^{e+}| #it{0.5<=p_{t}[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhpi0", "TrkIn :: <#Delta#phi>^{e+}| #it{0.8<=p_{t}[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhph0", "TrkIn :: <#Delta#phi>^{e+}| #it{1.5<=p_{t}[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhpH0", "TrkIn :: <#Delta#phi>^{e+}| #it{5.0<=p_{t}[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhp0", "TrkIn :: <#Delta#phi>^{e+} [deg]"},
    {"TRDresolution_TrkInPhSpL0", "TrkIn :: RMS<#Delta#phi>^{e+}| #it{p[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhSpl0", "TrkIn :: RMS<#Delta#phi>^{e+}| #it{0.5<=p[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhSpi0", "TrkIn :: RMS<#Delta#phi>^{e+}| #it{0.8<=p[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhSph0", "TrkIn :: RMS<#Delta#phi>^{e+}| #it{1.5<=p[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhSpH0", "TrkIn :: RMS<#Delta#phi>^{e+}| #it{5.0<=p[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhSp0", "TrkIn :: RMS<#Delta#phi>^{e+} [deg]"},
    {"TRDresolution_TrkInQpL0", "TrkIn :: MPV(dQdl)^{e+}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQpl0", "TrkIn :: MPV(dQdl)^{e+}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQpi0", "TrkIn :: MPV(dQdl)^{e+}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQph0", "TrkIn :: MPV(dQdl)^{e+}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQpH0", "TrkIn :: MPV(dQdl)^{e+}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQp0", "TrkIn :: MPV(dQdl)^{e+} [a.u.]"},
    {"TRDresolution_TrkInQSpL0", "TrkIn :: <dQdl>^{e+}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQSpl0", "TrkIn :: <dQdl>^{e+}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQSpi0", "TrkIn :: <dQdl>^{e+}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQSph0", "TrkIn :: <dQdl>^{e+}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQSpH0", "TrkIn :: <dQdl>^{e+}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQSp0", "TrkIn :: <dQdl>^{e+} [a.u.]"},
    {"TRDresolution_TrkInYpL1", "TrkIn :: <#Deltay>^{#mu#pi+}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYpl1", "TrkIn :: <#Deltay>^{#mu#pi+}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYpi1", "TrkIn :: <#Deltay>^{#mu#pi+}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYph1", "TrkIn :: <#Deltay>^{#mu#pi+}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYpH1", "TrkIn :: <#Deltay>^{#mu#pi+}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYp1", "TrkIn :: <#Deltay>^{#mu#pi+} [cm]"},
    {"TRDresolution_TrkInYSpL1", "TrkIn :: RMS(#Deltay)^{#mu#pi+}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYSpl1", "TrkIn :: RMS(#Deltay)^{#mu#pi+}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYSpi1", "TrkIn :: RMS(#Deltay)^{#mu#pi+}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYSph1", "TrkIn :: RMS(#Deltay)^{#mu#pi+}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYSpH1", "TrkIn :: RMS(#Deltay)^{#mu#pi+}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYSp1", "TrkIn :: RMS(#Deltay)^{#mu#pi+} [cm]"},
    {"TRDresolution_TrkInPhpL1", "TrkIn :: <#Delta#phi>^{#mu#pi+}| #it{p_{t}[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhpl1", "TrkIn :: <#Delta#phi>^{#mu#pi+}| #it{0.5<=p_{t}[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhpi1", "TrkIn :: <#Delta#phi>^{#mu#pi+}| #it{0.8<=p_{t}[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhph1", "TrkIn :: <#Delta#phi>^{#mu#pi+}| #it{1.5<=p_{t}[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhpH1", "TrkIn :: <#Delta#phi>^{#mu#pi+}| #it{5.0<=p_{t}[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhp1", "TrkIn :: <#Delta#phi>^{#mu#pi+} [deg]"},
    {"TRDresolution_TrkInPhSpL1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi+}| #it{p[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhSpl1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi+}| #it{0.5<=p[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhSpi1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi+}| #it{0.8<=p[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhSph1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi+}| #it{1.5<=p[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhSpH1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi+}| #it{5.0<=p[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhSp1", "TrkIn :: RMS<#Delta#phi>^{#mu#pi+} [deg]"},
    {"TRDresolution_TrkInQpL1", "TrkIn :: MPV(dQdl)^{#mu#pi+}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQpl1", "TrkIn :: MPV(dQdl)^{#mu#pi+}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQpi1", "TrkIn :: MPV(dQdl)^{#mu#pi+}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQph1", "TrkIn :: MPV(dQdl)^{#mu#pi+}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQpH1", "TrkIn :: MPV(dQdl)^{#mu#pi+}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQp1", "TrkIn :: MPV(dQdl)^{#mu#pi+} [a.u.]"},
    {"TRDresolution_TrkInQSpL1", "TrkIn :: <dQdl>^{#mu#pi+}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQSpl1", "TrkIn :: <dQdl>^{#mu#pi+}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQSpi1", "TrkIn :: <dQdl>^{#mu#pi+}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQSph1", "TrkIn :: <dQdl>^{#mu#pi+}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQSpH1", "TrkIn :: <dQdl>^{#mu#pi+}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQSp1", "TrkIn :: <dQdl>^{#mu#pi+} [a.u.]"},
    {"TRDresolution_TrkInYpL2", "TrkIn :: <#Deltay>^{Kp+}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYpl2", "TrkIn :: <#Deltay>^{Kp+}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYpi2", "TrkIn :: <#Deltay>^{Kp+}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYph2", "TrkIn :: <#Deltay>^{Kp+}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYpH2", "TrkIn :: <#Deltay>^{Kp+}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYp2", "TrkIn :: <#Deltay>^{Kp+} [cm]"},
    {"TRDresolution_TrkInYSpL2", "TrkIn :: RMS(#Deltay)^{Kp+}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYSpl2", "TrkIn :: RMS(#Deltay)^{Kp+}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYSpi2", "TrkIn :: RMS(#Deltay)^{Kp+}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYSph2", "TrkIn :: RMS(#Deltay)^{Kp+}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYSpH2", "TrkIn :: RMS(#Deltay)^{Kp+}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYSp2", "TrkIn :: RMS(#Deltay)^{Kp+} [cm]"},
    {"TRDresolution_TrkInPhpL2", "TrkIn :: <#Delta#phi>^{Kp+}| #it{p_{t}[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhpl2", "TrkIn :: <#Delta#phi>^{Kp+}| #it{0.5<=p_{t}[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhpi2", "TrkIn :: <#Delta#phi>^{Kp+}| #it{0.8<=p_{t}[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhph2", "TrkIn :: <#Delta#phi>^{Kp+}| #it{1.5<=p_{t}[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhpH2", "TrkIn :: <#Delta#phi>^{Kp+}| #it{5.0<=p_{t}[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhp2", "TrkIn :: <#Delta#phi>^{Kp+} [deg]"},
    {"TRDresolution_TrkInPhSpL2", "TrkIn :: RMS<#Delta#phi>^{Kp+}| #it{p[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhSpl2", "TrkIn :: RMS<#Delta#phi>^{Kp+}| #it{0.5<=p[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhSpi2", "TrkIn :: RMS<#Delta#phi>^{Kp+}| #it{0.8<=p[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhSph2", "TrkIn :: RMS<#Delta#phi>^{Kp+}| #it{1.5<=p[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhSpH2", "TrkIn :: RMS<#Delta#phi>^{Kp+}| #it{5.0<=p[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhSp2", "TrkIn :: RMS<#Delta#phi>^{Kp+} [deg]"},
    {"TRDresolution_TrkInQpL2", "TrkIn :: MPV(dQdl)^{Kp+}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQpl2", "TrkIn :: MPV(dQdl)^{Kp+}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQpi2", "TrkIn :: MPV(dQdl)^{Kp+}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQph2", "TrkIn :: MPV(dQdl)^{Kp+}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQpH2", "TrkIn :: MPV(dQdl)^{Kp+}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQp2", "TrkIn :: MPV(dQdl)^{Kp+} [a.u.]"},
    {"TRDresolution_TrkInQSpL2", "TrkIn :: <dQdl>^{Kp+}| #it{p[GeV/c]<0.5} [a.u.]"},
    {"TRDresolution_TrkInQSpl2", "TrkIn :: <dQdl>^{Kp+}| #it{0.5<=p[GeV/c]<0.8} [a.u.]"},
    {"TRDresolution_TrkInQSpi2", "TrkIn :: <dQdl>^{Kp+}| #it{0.8<=p[GeV/c]<1.5} [a.u.]"},
    {"TRDresolution_TrkInQSph2", "TrkIn :: <dQdl>^{Kp+}| #it{1.5<=p[GeV/c]<5.0} [a.u.]"},
    {"TRDresolution_TrkInQSpH2", "TrkIn :: <dQdl>^{Kp+}| #it{5.0<=p[GeV/c]} [a.u.]"},
    {"TRDresolution_TrkInQSp2", "TrkIn :: <dQdl>^{Kp+} [a.u.]"},
    {"TRDresolution_TrkInYpL", "TrkIn :: <#Deltay>^{+}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInYSpL", "TrkIn :: RMS(#Deltay)^{+}| #it{p_{t}[GeV/c]<0.5} [cm]"},
    {"TRDresolution_TrkInPhpL", "TrkIn :: <#Delta#phi>^{+}| #it{p_{t}[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInPhSpL", "TrkIn :: RMS<#Delta#phi>^{+}| #it{p_{t}[GeV/c]<0.5} [deg]"},
    {"TRDresolution_TrkInYpl", "TrkIn :: <#Deltay>^{+}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInYSpl", "TrkIn :: RMS(#Deltay)^{+}| #it{0.5<=p_{t}[GeV/c]<0.8} [cm]"},
    {"TRDresolution_TrkInPhpl", "TrkIn :: <#Delta#phi>^{+}| #it{0.5<=p_{t}[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInPhSpl", "TrkIn :: RMS<#Delta#phi>^{+}| #it{0.5<=p_{t}[GeV/c]<0.8} [deg]"},
    {"TRDresolution_TrkInYpi", "TrkIn :: <#Deltay>^{+}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInYSpi", "TrkIn :: RMS(#Deltay)^{+}| #it{0.8<=p_{t}[GeV/c]<1.5} [cm]"},
    {"TRDresolution_TrkInPhpi", "TrkIn :: <#Delta#phi>^{+}| #it{0.8<=p_{t}[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInPhSpi", "TrkIn :: RMS<#Delta#phi>^{+}| #it{0.8<=p_{t}[GeV/c]<1.5} [deg]"},
    {"TRDresolution_TrkInYph", "TrkIn :: <#Deltay>^{+}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInYSph", "TrkIn :: RMS(#Deltay)^{+}| #it{1.5<=p_{t}[GeV/c]<5.0} [cm]"},
    {"TRDresolution_TrkInPhph", "TrkIn :: <#Delta#phi>^{+}| #it{1.5<=p_{t}[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInPhSph", "TrkIn :: RMS<#Delta#phi>^{+}| #it{1.5<=p_{t}[GeV/c]<5.0} [deg]"},
    {"TRDresolution_TrkInYpH", "TrkIn :: <#Deltay>^{+}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInYSpH", "TrkIn :: RMS(#Deltay)^{+}| #it{5.0<=p_{t}[GeV/c]} [cm]"},
    {"TRDresolution_TrkInPhpH", "TrkIn :: <#Delta#phi>^{+}| #it{5.0<=p_{t}[GeV/c]} [deg]"},
    {"TRDresolution_TrkInPhSpH", "TrkIn :: RMS<#Delta#phi>^{+}| #it{5.0<=p_{t}[GeV/c]} [deg]"},
    {"TRDresolution_TrkInYp", "TrkIn :: <#Deltay>^{+} [cm]"},
    {"TRDresolution_TrkInYSp", "TrkIn :: RMS(#Deltay)^{+} [cm]"},
    {"TRDresolution_TrkInPhp", "TrkIn :: <#Delta#phi>^{+} [deg]"},
    {"TRDresolution_TrkInPhSp", "TrkIn :: RMS<#Delta#phi>^{+} [deg]"},
    {"TRDresolution_TrkInY0", "TrkIn :: <#Deltay>^{e} [cm]"},
    {"TRDresolution_TrkInPh0", "TrkIn :: <#Delta#phi>^{e} [deg]"},
    {"TRDresolution_TrkInQ0", "TrkIn :: MPV(dQdl)^{e} [a.u.]"},
    {"TRDresolution_TrkInQS0", "TrkIn :: <dQdl>^{e} [a.u.]"},
    {"TRDresolution_TrkInY1", "TrkIn :: <#Deltay>^{#mu#pi} [cm]"},
    {"TRDresolution_TrkInPh1", "TrkIn :: <#Delta#phi>^{#mu#pi} [deg]"},
    {"TRDresolution_TrkInQ1", "TrkIn :: MPV(dQdl)^{#mu#pi} [a.u.]"},
    {"TRDresolution_TrkInQS1", "TrkIn :: <dQdl>^{#mu#pi} [a.u.]"},
    {"TRDresolution_TrkInY2", "TrkIn :: <#Deltay>^{Kp} [cm]"},
    {"TRDresolution_TrkInPh2", "TrkIn :: <#Delta#phi>^{Kp} [deg]"},
    {"TRDresolution_TrkInQ2", "TrkIn :: MPV(dQdl)^{Kp} [a.u.]"},
    {"TRDresolution_TrkInQS2", "TrkIn :: <dQdl>^{Kp} [a.u.]"},
    {"TRDresolution_TrkInRCZ", "TrkIn :: <#Deltaz> [cm]"},
    {"TRDresolution_TrkInRCZS", "TrkIn :: #sigma(#Deltaz) [cm]"},
    {"TRDresolution_TrkInPtn0", "TrkIn :: <p_{t}>^{e-} [GeV/c]"},
    {"TRDresolution_TrkInPtp0", "TrkIn :: <p_{t}>^{e+} [GeV/c]"},
    {"TRDresolution_TrkInPtn1", "TrkIn :: <p_{t}>^{#mu#pi-} [GeV/c]"},
    {"TRDresolution_TrkInPtp1", "TrkIn :: <p_{t}>^{#mu#pi+} [GeV/c]"},
    {"TRDresolution_TrkInPtn2", "TrkIn :: <p_{t}>^{Kp-} [GeV/c]"},
    {"TRDresolution_TrkInPtp2", "TrkIn :: <p_{t}>^{Kp+} [GeV/c]"}
  };
  const char *resName[] = {"Markus Fasel", "Alexandru Bercuci"},
             *resMail[] = {"M.Fasel@gsi.de", "A.Bercuci@gsi.de"};
  const char *notName[] = {"Julian Book", "Hans Beck", "Ionut Arsene", "Raphaelle Bailache", "Christoph Blume"},
             *notMail[] = {"jbook@ikf.uni-frankfurt.de", "hbeck@ikf.uni-frankfurt.de", "I.C.Arsene@gsi.de", "R.Bailhache@gsi.de", "blume@ikf.uni-frankfurt.de"};
  Int_t nDet(0), nEff(0), nRes(0);
  for(Int_t jnt(0); jnt<nt; jnt++){
    //printf("%3d %s %s\n", jnt, tvn[jnt][0], tvn[jnt][1]);
    if(strstr(tvn[jnt][0], "TRDcheckDET")) nDet++;
    else if(strstr(tvn[jnt][0], "TRDefficiency")) nEff++;
    else if(strstr(tvn[jnt][0], "TRDresolution")) nRes++;
    else {
      Error("makeTrendingDB", "Entry \"%s\" not registered as trending task.", tvn[jnt][0]);
      return;
    }
  }
  Info("makeTrendingDB", "Trends :: %3d = %3d[DET] %3d[EFF] %3d[RES]", nt, nDet, nEff, nRes);

  TFile *fDB = TFile::Open("TRD.TrendDB.root", "RECREATE");
  TTree *tDB = new TTree("trend", "Reference Trend Values");
  Double_t val[nt+1000]; Int_t jt(0);
  TBranch *b(NULL);
  for(Int_t it(0); it<nt; it++){  // ALL
    b = tDB->Branch(tvn[it][0], &val[jt++], Form("%s/D", tvn[it][0]));
    b->SetTitle(tvn[it][1]);
  }
  for(Int_t it(nDet); it<nDet+nEff; it++){ // extra EFF (MC)
    TString stn("TRDefficiency_MC"); stn+=&tvn[it][0][14];
    b = tDB->Branch(stn.Data(), &val[jt++], Form("%s/D", stn.Data()));
    b->SetTitle(Form("[MC] %s", tvn[it][1]));
  }
  for(Int_t it(nDet+nEff); it<nDet+nEff+nRes; it++){ // extra RES (MC)
    TString stn("TRDresolution_MC"); stn+=&tvn[it][0][14];
    b= tDB->Branch(stn.Data(), &val[jt++], Form("%s/D", stn.Data()));
    b->SetTitle(Form("[MC] %s", tvn[it][1]));
  }
  for(Int_t it(nDet+nEff); it<nDet+nEff+nRes; it++){ // extra RES (V0)
    TString stn("TRDresolution_V0"); stn+=&tvn[it][0][14];
    b = tDB->Branch(stn.Data(), &val[jt++], Form("%s/D", stn.Data()));
    b->SetTitle(Form("[V0] %s", tvn[it][1]));
  }
  gROOT->cd();
  Info("makeTrendingDB", "Trends :: Combined [%3d]", tDB->GetNbranches());
  TIterator *ib = tDB->GetListOfBranches()->MakeIterator();
  
  AliTRDtrendValue *tv(NULL);
  FILE *fp = fopen(fl, "rt");
  TString sfp;
  while(sfp.Gets(fp)){
    if(!TFile::Open(sfp.Data())) continue;
    Int_t nmiss(0), nbads(0), it(-1); ib->Reset();
    while((b=(TBranch*)ib->Next())){
      val[++it] = -999;
      if(!(tv = (AliTRDtrendValue*)gFile->Get(b->GetName()))) {
        //Warning("makeTrendingDB()", "Missing %s from %s", b->GetName(), sfp.Data());
        nmiss++;
        continue;
      }
      if((strstr(b->GetName(), "TRDcheckDET") ||
          strstr(b->GetName(), "QS") ||
          strstr(b->GetName(), "YS")) && TMath::Abs(tv->GetVal()) < 1.e-5){
        //Info("makeTrendingDB()", "Found bad value for %s[%f] in %s", b->GetName(), tv->GetVal(), sfp.Data());
        nbads++;
        continue;
      }
      val[it] = tv->GetVal();
    }
    Warning("makeTrendingDB()", "%s :: Missing[%d] Bads[%d]", sfp.Data(), nmiss, nbads);
    gFile->Close();
    tDB->Fill();
  }

//   TFile *fDB = TFile::Open("TRD.TrendDB.root");
//   TTree *tDB = (TTree*)gFile->Get("trend");
//   Double_t val[nt];
//   for(Int_t it(0); it<nt; it++) tDB->SetBranchAddress(tvn[it][0], &val[it]);
//   gROOT->cd();

  TString res[] = {Form("%s/%s", resName[0], resMail[0]), Form("%s/%s", resName[1], resMail[1])};
  TString notifiable;
  for(Int_t inot(0); inot<5; inot++){
    notifiable+=notName[inot];
    notifiable+="/";
    notifiable+=notMail[inot];
    if(inot<4) notifiable+=",";
  }
  TF1 f("f", "gaus", -100, 100); TH1 *h(NULL);
  AliTRDtrendingManager *tm = AliTRDtrendingManager::Instance();
  TCanvas *c = new TCanvas("c", "Trend Distrib.", 10, 10, 500, 500);
  Int_t ntr=tDB->GetEntries();
  Int_t it(-1); ib->Reset();
  while((b=(TBranch*)ib->Next())){
    tDB->Draw(b->GetName(), "", "goff"); it++;
    Double_t *v = tDB->GetV1(), xmin(100.), xmax(-100);
    Int_t ntr0(0);
    for(Int_t ir=0; ir<ntr; ir++){
      if(v[ir]<-100) continue;
      ntr0++;
      if(v[ir]<xmin) xmin = v[ir];
      if(v[ir]>xmax) xmax = v[ir];
    }
    Double_t m(0.), s(0.);
    if(ntr0<10){
      Warning("makeTrendingDB", "%s :: Couldn't create reference value out of %d entries.", b->GetName(), ntr0);
    } else {
      if((h =(TH1F*)gROOT->FindObject("hp"))){delete h; h = NULL;}
      h = new TH1F("hp", Form("%s;%s;entries", b->GetTitle(), b->GetName()), 25, 0.5*(3*xmin-xmax), 0.5*(3*xmax - xmin));
      tDB->Draw(Form("%s>>hp", b->GetName()), Form("%s>-100", b->GetName()));
      if(h->Integral() < 1) continue;
      f.SetParameter(0, h->Integral());
      f.SetParameter(1, h->GetMean());
      f.SetParameter(2, h->GetRMS());
      h->Fit(&f, "WQ");
      c->Modified(); c->Update(); c->SaveAs(Form("Fit_%s.gif", b->GetName()));

      // write trending value to manager
      Info("makeTrendingDB", "%s :: %f+-%f [%f - %f]", b->GetName(), f.GetParameter(1), f.GetParameter(2), xmin, xmax);
  /*    if(strstr(tvn[it][0], "TrkInYS")) {
        m=0.4; s=0.06;
      } else if(strstr(tvn[it][0], "TrkInY")) {
        m=0.; s=0.1;*/
  /*    } else if(strstr(tvn[it][0], "TrkInPh")) {
        m=0.; s=0.35;*/
  /*    } else if(strstr(tvn[it][0], "TrkInQ") || strstr(tvn[it][0], "TrkInQS")) {
        m=-2.; s=0.2;*/
  //    } else {
        m=f.GetParameter(1); s=h->GetRMS()/*f.GetParameter(2)*/;
  //    }
    }
    tm->AddValue(b->GetName(), m, s, b->GetTitle(), res[it>13], notifiable);
  }
  tm->Terminate();

  fDB->cd();
  tDB->Write();
  fDB->Close();
}
