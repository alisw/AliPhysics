#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TError.h"
#include "PWGPP/TRD/info/AliTRDtrendingManager.h"

#endif

void makeTrendingDB(const Char_t *fl)
{
// Make trending of variable list "tl" from trending file list "fl"
// The trending value list should be formated "var1:var2:var3"
// The trending file from the list should be found on a path formated "your path"/runId/TRD.PerformanceTrend.root 
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGPP.so");
  gSystem->Load("libPWGmuon.so");

  const Int_t nt(167);
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
    {"TRDresolution_TrkInPhnL0", "TrkIn :: <#Delta#phi>^{e-}(p_{t}[GeV/c]<0.5) [deg]"},
    {"TRDresolution_TrkInPhnl0", "TrkIn :: <#Delta#phi>^{e-}(0.5<=p_{t}[GeV/c]<0.8) [deg]"},
    {"TRDresolution_TrkInPhni0", "TrkIn :: <#Delta#phi>^{e-}(0.8<=p_{t}[GeV/c]<1.5) [deg]"},
    {"TRDresolution_TrkInPhnh0", "TrkIn :: <#Delta#phi>^{e-}(1.5<=p_{t}[GeV/c]<5.0) [deg]"},
    {"TRDresolution_TrkInPhnH0", "TrkIn :: <#Delta#phi>^{e-}(5.0<=p_{t}[GeV/c]) [deg]"},
    {"TRDresolution_TrkInPhn0", "TrkIn :: <#Delta#phi>^{e-} [deg]"},
    {"TRDresolution_TrkInQnL0", "TrkIn :: MPV(dQdl)^{e-}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQnl0", "TrkIn :: MPV(dQdl)^{e-}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQni0", "TrkIn :: MPV(dQdl)^{e-}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQnh0", "TrkIn :: MPV(dQdl)^{e-}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQnH0", "TrkIn :: MPV(dQdl)^{e-}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQn0", "TrkIn :: MPV(dQdl)^{e-} [a.u.]"},
    {"TRDresolution_TrkInQSnL0", "TrkIn :: <dQdl>^{e-}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQSnl0", "TrkIn :: <dQdl>^{e-}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQSni0", "TrkIn :: <dQdl>^{e-}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQSnh0", "TrkIn :: <dQdl>^{e-}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQSnH0", "TrkIn :: <dQdl>^{e-}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQSn0", "TrkIn :: <dQdl>^{e-} [a.u.]"},
    {"TRDresolution_TrkInPhnL1", "TrkIn :: <#Delta#phi>^{#mu#pi-}(p_{t}[GeV/c]<0.5) [deg]"},
    {"TRDresolution_TrkInPhnl1", "TrkIn :: <#Delta#phi>^{#mu#pi-}(0.5<=p_{t}[GeV/c]<0.8) [deg]"},
    {"TRDresolution_TrkInPhni1", "TrkIn :: <#Delta#phi>^{#mu#pi-}(0.8<=p_{t}[GeV/c]<1.5) [deg]"},
    {"TRDresolution_TrkInPhnh1", "TrkIn :: <#Delta#phi>^{#mu#pi-}(1.5<=p_{t}[GeV/c]<5.0) [deg]"},
    {"TRDresolution_TrkInPhnH1", "TrkIn :: <#Delta#phi>^{#mu#pi-}(5.0<=p_{t}[GeV/c]) [deg]"},
    {"TRDresolution_TrkInPhn1", "TrkIn :: <#Delta#phi>^{#mu#pi-} [deg]"},
    {"TRDresolution_TrkInQnL1", "TrkIn :: MPV(dQdl)^{#mu#pi-}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQnl1", "TrkIn :: MPV(dQdl)^{#mu#pi-}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQni1", "TrkIn :: MPV(dQdl)^{#mu#pi-}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQnh1", "TrkIn :: MPV(dQdl)^{#mu#pi-}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQnH1", "TrkIn :: MPV(dQdl)^{#mu#pi-}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQn1", "TrkIn :: MPV(dQdl)^{#mu#pi-} [a.u.]"},
    {"TRDresolution_TrkInQSnL1", "TrkIn :: <dQdl>^{#mu#pi-}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQSnl1", "TrkIn :: <dQdl>^{#mu#pi-}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQSni1", "TrkIn :: <dQdl>^{#mu#pi-}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQSnh1", "TrkIn :: <dQdl>^{#mu#pi-}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQSnH1", "TrkIn :: <dQdl>^{#mu#pi-}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQSn1", "TrkIn :: <dQdl>^{#mu#pi-} [a.u.]"},
    {"TRDresolution_TrkInPhnL2", "TrkIn :: <#Delta#phi>^{Kp-}(p_{t}[GeV/c]<0.5) [deg]"},
    {"TRDresolution_TrkInPhnl2", "TrkIn :: <#Delta#phi>^{Kp-}(0.5<=p_{t}[GeV/c]<0.8) [deg]"},
    {"TRDresolution_TrkInPhni2", "TrkIn :: <#Delta#phi>^{Kp-}(0.8<=p_{t}[GeV/c]<1.5) [deg]"},
    {"TRDresolution_TrkInPhnh2", "TrkIn :: <#Delta#phi>^{Kp-}(1.5<=p_{t}[GeV/c]<5.0) [deg]"},
    {"TRDresolution_TrkInPhnH2", "TrkIn :: <#Delta#phi>^{Kp-}(5.0<=p_{t}[GeV/c]) [deg]"},
    {"TRDresolution_TrkInPhn2", "TrkIn :: <#Delta#phi>^{Kp-} [deg]"},
    {"TRDresolution_TrkInQnL2", "TrkIn :: MPV(dQdl)^{Kp-}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQnl2", "TrkIn :: MPV(dQdl)^{Kp-}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQni2", "TrkIn :: MPV(dQdl)^{Kp-}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQnh2", "TrkIn :: MPV(dQdl)^{Kp-}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQnH2", "TrkIn :: MPV(dQdl)^{Kp-}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQn2", "TrkIn :: MPV(dQdl)^{Kp-} [a.u.]"},
    {"TRDresolution_TrkInQSnL2", "TrkIn :: <dQdl>^{Kp-}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQSnl2", "TrkIn :: <dQdl>^{Kp-}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQSni2", "TrkIn :: <dQdl>^{Kp-}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQSnh2", "TrkIn :: <dQdl>^{Kp-}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQSnH2", "TrkIn :: <dQdl>^{Kp-}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQSn2", "TrkIn :: <dQdl>^{Kp-} [a.u.]"},
    {"TRDresolution_TrkInYnL", "TrkIn :: <#Deltay>^{-}(p_{t}[GeV/c]<0.5) [cm]"},
    {"TRDresolution_TrkInYSnL", "TrkIn :: RMS(#Deltay)^{-}(p_{t}[GeV/c]<0.5) [cm]"},
    {"TRDresolution_TrkInPhnL", "TrkIn :: <#Delta#phi>^{-}(p_{t}[GeV/c]<0.5) [deg]"},
    {"TRDresolution_TrkInYnl", "TrkIn :: <#Deltay>^{-}(0.5<=p_{t}[GeV/c]<0.8) [cm]"},
    {"TRDresolution_TrkInYSnl", "TrkIn :: RMS(#Deltay)^{-}(0.5<=p_{t}[GeV/c]<0.8) [cm]"},
    {"TRDresolution_TrkInPhnl", "TrkIn :: <#Delta#phi>^{-}(0.5<=p_{t}[GeV/c]<0.8) [deg]"},
    {"TRDresolution_TrkInYni", "TrkIn :: <#Deltay>^{-}(0.8<=p_{t}[GeV/c]<1.5) [cm]"},
    {"TRDresolution_TrkInYSni", "TrkIn :: RMS(#Deltay)^{-}(0.8<=p_{t}[GeV/c]<1.5) [cm]"},
    {"TRDresolution_TrkInPhni", "TrkIn :: <#Delta#phi>^{-}(0.8<=p_{t}[GeV/c]<1.5) [deg]"},
    {"TRDresolution_TrkInYnh", "TrkIn :: <#Deltay>^{-}(1.5<=p_{t}[GeV/c]<5.0) [cm]"},
    {"TRDresolution_TrkInYSnh", "TrkIn :: RMS(#Deltay)^{-}(1.5<=p_{t}[GeV/c]<5.0) [cm]"},
    {"TRDresolution_TrkInPhnh", "TrkIn :: <#Delta#phi>^{-}(1.5<=p_{t}[GeV/c]<5.0) [deg]"},
    {"TRDresolution_TrkInYnH", "TrkIn :: <#Deltay>^{-}(5.0<=p_{t}[GeV/c]) [cm]"},
    {"TRDresolution_TrkInYSnH", "TrkIn :: RMS(#Deltay)^{-}(5.0<=p_{t}[GeV/c]) [cm]"},
    {"TRDresolution_TrkInPhnH", "TrkIn :: <#Delta#phi>^{-}(5.0<=p_{t}[GeV/c]) [deg]"},
    {"TRDresolution_TrkInYn", "TrkIn :: <#Deltay>^{-} [cm]"},
    {"TRDresolution_TrkInYSn", "TrkIn :: RMS(#Deltay)^{-} [cm]"},
    {"TRDresolution_TrkInPhn", "TrkIn :: <#Delta#phi>^{-} [deg]"},
    {"TRDresolution_TrkInPhpL0", "TrkIn :: <#Delta#phi>^{e+}(p_{t}[GeV/c]<0.5) [deg]"},
    {"TRDresolution_TrkInPhpl0", "TrkIn :: <#Delta#phi>^{e+}(0.5<=p_{t}[GeV/c]<0.8) [deg]"},
    {"TRDresolution_TrkInPhpi0", "TrkIn :: <#Delta#phi>^{e+}(0.8<=p_{t}[GeV/c]<1.5) [deg]"},
    {"TRDresolution_TrkInPhph0", "TrkIn :: <#Delta#phi>^{e+}(1.5<=p_{t}[GeV/c]<5.0) [deg]"},
    {"TRDresolution_TrkInPhpH0", "TrkIn :: <#Delta#phi>^{e+}(5.0<=p_{t}[GeV/c]) [deg]"},
    {"TRDresolution_TrkInPhp0", "TrkIn :: <#Delta#phi>^{e+} [deg]"},
    {"TRDresolution_TrkInQpL0", "TrkIn :: MPV(dQdl)^{e+}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQpl0", "TrkIn :: MPV(dQdl)^{e+}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQpi0", "TrkIn :: MPV(dQdl)^{e+}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQph0", "TrkIn :: MPV(dQdl)^{e+}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQpH0", "TrkIn :: MPV(dQdl)^{e+}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQp0", "TrkIn :: MPV(dQdl)^{e+} [a.u.]"},
    {"TRDresolution_TrkInQSpL0", "TrkIn :: <dQdl>^{e+}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQSpl0", "TrkIn :: <dQdl>^{e+}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQSpi0", "TrkIn :: <dQdl>^{e+}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQSph0", "TrkIn :: <dQdl>^{e+}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQSpH0", "TrkIn :: <dQdl>^{e+}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQSp0", "TrkIn :: <dQdl>^{e+} [a.u.]"},
    {"TRDresolution_TrkInPhpL1", "TrkIn :: <#Delta#phi>^{#mu#pi+}(p_{t}[GeV/c]<0.5) [deg]"},
    {"TRDresolution_TrkInPhpl1", "TrkIn :: <#Delta#phi>^{#mu#pi+}(0.5<=p_{t}[GeV/c]<0.8) [deg]"},
    {"TRDresolution_TrkInPhpi1", "TrkIn :: <#Delta#phi>^{#mu#pi+}(0.8<=p_{t}[GeV/c]<1.5) [deg]"},
    {"TRDresolution_TrkInPhph1", "TrkIn :: <#Delta#phi>^{#mu#pi+}(1.5<=p_{t}[GeV/c]<5.0) [deg]"},
    {"TRDresolution_TrkInPhpH1", "TrkIn :: <#Delta#phi>^{#mu#pi+}(5.0<=p_{t}[GeV/c]) [deg]"},
    {"TRDresolution_TrkInPhp1", "TrkIn :: <#Delta#phi>^{#mu#pi+} [deg]"},
    {"TRDresolution_TrkInQpL1", "TrkIn :: MPV(dQdl)^{#mu#pi+}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQpl1", "TrkIn :: MPV(dQdl)^{#mu#pi+}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQpi1", "TrkIn :: MPV(dQdl)^{#mu#pi+}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQph1", "TrkIn :: MPV(dQdl)^{#mu#pi+}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQpH1", "TrkIn :: MPV(dQdl)^{#mu#pi+}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQp1", "TrkIn :: MPV(dQdl)^{#mu#pi+} [a.u.]"},
    {"TRDresolution_TrkInQSpL1", "TrkIn :: <dQdl>^{#mu#pi+}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQSpl1", "TrkIn :: <dQdl>^{#mu#pi+}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQSpi1", "TrkIn :: <dQdl>^{#mu#pi+}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQSph1", "TrkIn :: <dQdl>^{#mu#pi+}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQSpH1", "TrkIn :: <dQdl>^{#mu#pi+}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQSp1", "TrkIn :: <dQdl>^{#mu#pi+} [a.u.]"},
    {"TRDresolution_TrkInPhpL2", "TrkIn :: <#Delta#phi>^{Kp+}(p_{t}[GeV/c]<0.5) [deg]"},
    {"TRDresolution_TrkInPhpl2", "TrkIn :: <#Delta#phi>^{Kp+}(0.5<=p_{t}[GeV/c]<0.8) [deg]"},
    {"TRDresolution_TrkInPhpi2", "TrkIn :: <#Delta#phi>^{Kp+}(0.8<=p_{t}[GeV/c]<1.5) [deg]"},
    {"TRDresolution_TrkInPhph2", "TrkIn :: <#Delta#phi>^{Kp+}(1.5<=p_{t}[GeV/c]<5.0) [deg]"},
    {"TRDresolution_TrkInPhpH2", "TrkIn :: <#Delta#phi>^{Kp+}(5.0<=p_{t}[GeV/c]) [deg]"},
    {"TRDresolution_TrkInPhp2", "TrkIn :: <#Delta#phi>^{Kp+} [deg]"},
    {"TRDresolution_TrkInQpL2", "TrkIn :: MPV(dQdl)^{Kp+}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQpl2", "TrkIn :: MPV(dQdl)^{Kp+}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQpi2", "TrkIn :: MPV(dQdl)^{Kp+}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQph2", "TrkIn :: MPV(dQdl)^{Kp+}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQpH2", "TrkIn :: MPV(dQdl)^{Kp+}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQp2", "TrkIn :: MPV(dQdl)^{Kp+} [a.u.]"},
    {"TRDresolution_TrkInQSpL2", "TrkIn :: <dQdl>^{Kp+}(p_{t}[GeV/c]<0.5) [a.u.]"},
    {"TRDresolution_TrkInQSpl2", "TrkIn :: <dQdl>^{Kp+}(0.5<=p_{t}[GeV/c]<0.8) [a.u.]"},
    {"TRDresolution_TrkInQSpi2", "TrkIn :: <dQdl>^{Kp+}(0.8<=p_{t}[GeV/c]<1.5) [a.u.]"},
    {"TRDresolution_TrkInQSph2", "TrkIn :: <dQdl>^{Kp+}(1.5<=p_{t}[GeV/c]<5.0) [a.u.]"},
    {"TRDresolution_TrkInQSpH2", "TrkIn :: <dQdl>^{Kp+}(5.0<=p_{t}[GeV/c]) [a.u.]"},
    {"TRDresolution_TrkInQSp2", "TrkIn :: <dQdl>^{Kp+} [a.u.]"},
    {"TRDresolution_TrkInYpL", "TrkIn :: <#Deltay>^{+}(p_{t}[GeV/c]<0.5) [cm]"},
    {"TRDresolution_TrkInYSpL", "TrkIn :: RMS(#Deltay)^{+}(p_{t}[GeV/c]<0.5) [cm]"},
    {"TRDresolution_TrkInPhpL", "TrkIn :: <#Delta#phi>^{+}(p_{t}[GeV/c]<0.5) [deg]"},
    {"TRDresolution_TrkInYpl", "TrkIn :: <#Deltay>^{+}(0.5<=p_{t}[GeV/c]<0.8) [cm]"},
    {"TRDresolution_TrkInYSpl", "TrkIn :: RMS(#Deltay)^{+}(0.5<=p_{t}[GeV/c]<0.8) [cm]"},
    {"TRDresolution_TrkInPhpl", "TrkIn :: <#Delta#phi>^{+}(0.5<=p_{t}[GeV/c]<0.8) [deg]"},
    {"TRDresolution_TrkInYpi", "TrkIn :: <#Deltay>^{+}(0.8<=p_{t}[GeV/c]<1.5) [cm]"},
    {"TRDresolution_TrkInYSpi", "TrkIn :: RMS(#Deltay)^{+}(0.8<=p_{t}[GeV/c]<1.5) [cm]"},
    {"TRDresolution_TrkInPhpi", "TrkIn :: <#Delta#phi>^{+}(0.8<=p_{t}[GeV/c]<1.5) [deg]"},
    {"TRDresolution_TrkInYph", "TrkIn :: <#Deltay>^{+}(1.5<=p_{t}[GeV/c]<5.0) [cm]"},
    {"TRDresolution_TrkInYSph", "TrkIn :: RMS(#Deltay)^{+}(1.5<=p_{t}[GeV/c]<5.0) [cm]"},
    {"TRDresolution_TrkInPhph", "TrkIn :: <#Delta#phi>^{+}(1.5<=p_{t}[GeV/c]<5.0) [deg]"},
    {"TRDresolution_TrkInYpH", "TrkIn :: <#Deltay>^{+}(5.0<=p_{t}[GeV/c]) [cm]"},
    {"TRDresolution_TrkInYSpH", "TrkIn :: RMS(#Deltay)^{+}(5.0<=p_{t}[GeV/c]) [cm]"},
    {"TRDresolution_TrkInPhpH", "TrkIn :: <#Delta#phi>^{+}(5.0<=p_{t}[GeV/c]) [deg]"},
    {"TRDresolution_TrkInYp", "TrkIn :: <#Deltay>^{+} [cm]"},
    {"TRDresolution_TrkInYSp", "TrkIn :: RMS(#Deltay)^{+} [cm]"},
    {"TRDresolution_TrkInPhp", "TrkIn :: <#Delta#phi>^{+} [deg]"},
    {"TRDresolution_TrkInPh0", "TrkIn :: <#Delta#phi>^{e} [deg]"},
    {"TRDresolution_TrkInQ0", "TrkIn :: MPV(dQdl)^{e} [a.u.]"},
    {"TRDresolution_TrkInQS0", "TrkIn :: <dQdl>^{e} [a.u.]"},
    {"TRDresolution_TrkInPh1", "TrkIn :: <#Delta#phi>^{#mu#pi} [deg]"},
    {"TRDresolution_TrkInQ1", "TrkIn :: MPV(dQdl)^{#mu#pi} [a.u.]"},
    {"TRDresolution_TrkInQS1", "TrkIn :: <dQdl>^{#mu#pi} [a.u.]"},
    {"TRDresolution_TrkInPh2", "TrkIn :: <#Delta#phi>^{Kp} [deg]"},
    {"TRDresolution_TrkInQ2", "TrkIn :: MPV(dQdl)^{Kp} [a.u.]"},
    {"TRDresolution_TrkInQS2", "TrkIn :: <dQdl>^{Kp} [a.u.]"}
  };
  const char *resName[] = {"Markus Fasel", "Alexandru Bercuci"},
             *resMail[] = {"M.Fasel@gsi.de", "A.Bercuci@gsi.de"};
  const char *notName[] = {"Julian Book", "Hans Beck", "Ionut Arsene", "Raphaelle Bailache", "Christoph Blume"},
             *notMail[] = {"jbook@ikf.uni-frankfurt.de", "hbeck@ikf.uni-frankfurt.de", "I.C.Arsene@gsi.de", "R.Bailhache@gsi.de", "blume@ikf.uni-frankfurt.de"};

  TFile *fDB = TFile::Open("TRD.TrendDB.root", "RECREATE");
  TTree *tDB = new TTree("trend", "Reference Trend Values");
  Double_t val[nt];
  for(Int_t it(0); it<nt; it++) tDB->Branch(tvn[it][0], &val[it], Form("%s/D", tvn[it][0]));
  gROOT->cd();

  AliTRDtrendValue *tv(NULL);
  FILE *fp = fopen(fl, "rt");
  TString sfp;
  while(sfp.Gets(fp)){
    if(!TFile::Open(sfp.Data())) continue;
    for(Int_t it(0); it<nt; it++){
      val[it] = -999;
      if(!(tv = (AliTRDtrendValue*)gFile->Get(tvn[it][0]))) {
        Warning("makeTrendingDB()", "Missing %s from %s", tvn[it][0], sfp.Data());
        continue;
      }
      if((strstr(tvn[it][0], "QS") || strstr(tvn[it][0], "YS")) &&tv->GetVal() < 1.e-5){
        Info("makeTrendingDB()", "Found bad value for %s[%f] in %s", tvn[it][0], tv->GetVal(), sfp.Data());
        continue;
      }
      val[it] = tv->GetVal();
    }
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
  for(Int_t it(0); it<nt; it++){
    tDB->Draw(tvn[it][0], "", "goff");
    Double_t *v = tDB->GetV1(), xmin(100.), xmax(-100);
    Int_t ntr0(0);
    for(Int_t ir=0; ir<ntr; ir++){
      if(v[ir]<-100) continue;
      ntr0++;
      if(v[ir]<xmin) xmin = v[ir];
      if(v[ir]>xmax) xmax = v[ir];
    }
    if(ntr0<10){
      Warning("makeTrendingDB", "Couldn't create entry %s. Too few values %d", tvn[it][0], ntr0);
      continue;
    }
    if((h =(TH1F*)gROOT->FindObject("hp"))){delete h; h = NULL;}
    h = new TH1F("hp", Form(";%s;entries", tvn[it][0]), 10, 0.5*(3*xmin-xmax), 0.5*(3*xmax - xmin));
    tDB->Draw(Form("%s>>hp", tvn[it][0]), Form("%s>-100", tvn[it][0]));
    if(h->Integral() < 1) continue;
    f.SetParameter(0, h->Integral());
    f.SetParameter(1, h->GetMean());
    f.SetParameter(2, h->GetRMS());
    h->Fit(&f, "WQ");
    c->Modified(); c->Update(); c->SaveAs(Form("Fit_%s.gif", tvn[it][0]));

    // write trending value to manager
    Info("makeTrendingDB", "%s [%f - %f] %f[%f]", tvn[it][0], xmin, xmax, f.GetParameter(1), f.GetParameter(2));
    Double_t m(0.), s(0.);
    if(strstr(tvn[it][0], "TrkInY")) {
      m=0.; s=0.1;
    } else if(strstr(tvn[it][0], "TrkInPh")) {
      m=0.; s=0.35;
    } else if(strstr(tvn[it][0], "TrkInQ") || strstr(tvn[it][0], "TrkInQS")) {
      m=-2.; s=0.2;
    } else {
      m=f.GetParameter(1); s=f.GetParameter(2);
    }
    tm->AddValue(tvn[it][0], m, s, tvn[it][1], res[it>13], notifiable);
  }
  tm->Terminate();

  fDB->cd();
  tDB->Write();
  fDB->Close();
}
