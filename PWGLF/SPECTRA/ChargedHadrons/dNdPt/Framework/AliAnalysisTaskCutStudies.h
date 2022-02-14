/// \class AliAnalysisTaskCutStudies

#ifndef AliAnalysisTaskCutStudies_H
#define AliAnalysisTaskCutStudies_H

#include "AliAnalysisTaskMKBase.h"
#include "AliAnalysisHelpersHist.h"

class AliESDtrackCuts;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class AliESDtrack;
class AliMCParticle;

class AliAnalysisTaskCutStudies : public AliAnalysisTaskMKBase
{
public:
  static AliESDtrackCuts* GetCutSetting(const std::string& identifier);

  AliAnalysisTaskCutStudies();
  AliAnalysisTaskCutStudies(const char* name);
  virtual ~AliAnalysisTaskCutStudies();

  virtual void AddOutput();                   // called at the beginning
  virtual Bool_t IsEventSelected();           // called for each event
  virtual void AnaEvent();                    // called once for every selected event
  virtual void AnaTrack(Int_t flag = 0);      // called once for every track in DATA+MC event
  virtual void AnaTrackMC(Int_t flag = 0);    // called once for every track in DATA event
  virtual void AnaParticleMC(Int_t flag = 0); // called once for every track in MC event

  static AliAnalysisTaskCutStudies* AddTaskCutStudies(const char* name = "TaskCutStudies");

private:
  AliAnalysisTaskCutStudies(const AliAnalysisTaskCutStudies&);            // not implemented
  AliAnalysisTaskCutStudies& operator=(const AliAnalysisTaskCutStudies&); // not implemented

  // track related properties
  Hist::Hist<TH1D> fHist_x{};         //!<!  x at dca (radial distance to vertex)
  Hist::Hist<TH1D> fHist_y{};         //!<!  local Y-coordinate of track at dca  (cm)
  Hist::Hist<TH1D> fHist_z{};         //!<!  local Z-coordinate of track at dca  (cm)
  Hist::Hist<TH1D> fHist_alpha{};     //!<!  local to global angle
  Hist::Hist<TH1D> fHist_signed1Pt{}; //!<!  signed 1/pt (1/(GeV/c))
  Hist::Hist<TH1D> fHist_snp{};       //!<!  local sine of the track momentum azimuthal angle
  Hist::Hist<TH1D> fHist_tgl{};       //!<!  tangent of the track momentum dip angle
  Hist::Hist<TH1D> fHist_dcaxy{};     //!<!  distance of closest approach in xy plane
  Hist::Hist<TH1D> fHist_dcaz{};      //!<!  distance of closest approach in beam direction z
  Hist::Hist<TH1D> fHist_flag{};      //!<!  flag info assigned to the track
  Hist::Hist<TH1D> fHist_pt{};        //!<!  transverse momentum
  Hist::Hist<TH1D> fHist_eta{};       //!<!  pseudorapidity
  Hist::Hist<TH1D> fHist_phi{};       //!<!  azimuthal angle phi
  Hist::Hist<TH1D> fHist_zInner{};    //!<!  z at inner param

  // its related properties
  Hist::Hist<TH1D> fHist_itsFoundClusters{};  //!<!  found clusters ITS
  Hist::Hist<TH1D> fHist_itsChi2PerCluster{}; //!<!  chi2 per cluster ITS
  Hist::Hist<TH1D> fHist_itsHits{};           //!<!  hitmap ITS

  // tpc related properties
  Hist::Hist<THnF> fHist_tpcFindableClusters{};       //!<!  findable clusters TPC
  Hist::Hist<THnF> fHist_tpcFoundClusters{};          //!<!  found clusters TPC
  Hist::Hist<TH1D> fHist_tpcSharedClusters{};         //!<!  shared clusters TPC
  Hist::Hist<TH1D> fHist_tpcFractionSharedClusters{}; //!<!  fraction of shared clusters TPC
  Hist::Hist<THnF> fHist_tpcCrossedRows{};            //!<!  crossed rows in TPC
  Hist::Hist<THnF> fHist_tpcCrossedRowsOverFindableClusters{}; //!<!  rows / findable clusters TPC
  Hist::Hist<THnF> fHist_tpcChi2PerCluster{};                  //!<!  chi2 per cluster TPC

  Hist::Hist<THnF> fHist_tpcNClustersPID{};      //!<!  number of clusters used for PID in TPC
  Hist::Hist<THnF> fHist_tpcGoldenChi2{};        //!<! chi2 global vs tpc constrained track
  Hist::Hist<THnF> fHist_tpcGeomLength{};        //!<! track length in active volume of the TPC
  Hist::Hist<TH2D> fHist_correlChi2GeomLength{}; //!<! chi2 per cluster vs geom length

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCutStudies, 1);
  /// \endcond
};

#endif
