/// \class AliAnalysisTaskCutStudies

#ifndef AliAnalysisTaskCutStudies_H
#define AliAnalysisTaskCutStudies_H

#include "AliAnalysisTaskMKBase.h"
#include "Hist.h"

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

using namespace AnalysisHelpers; // TODO: remove this at some point to avoid polluting the global namespace!!

class AliAnalysisTaskCutStudies : public AliAnalysisTaskMKBase
{
public:
  AliAnalysisTaskCutStudies();
  AliAnalysisTaskCutStudies(const char *name);
  virtual ~AliAnalysisTaskCutStudies();

  virtual void AddOutput();                     //called at the beginning
  virtual Bool_t IsEventSelected();             //called for each event
  virtual void AnaEvent();                      //called once for every selected event
  virtual void AnaTrack(Int_t flag = 0);        //called once for every track in DATA+MC event
  virtual void AnaTrackMC(Int_t flag = 0);      //called once for every track in DATA event
  virtual void AnaParticleMC(Int_t flag = 0);   //called once for every track in MC event

  static AliAnalysisTaskCutStudies* AddTaskCutStudies(const char* name = "TaskCutStudies");
  
private:
  
  
  AliAnalysisTaskCutStudies(const AliAnalysisTaskCutStudies&); // not implemented
  AliAnalysisTaskCutStudies& operator=(const AliAnalysisTaskCutStudies&); // not implemented

  using RootHist_t = TH1D;
  
  // track related properties
  Hist<RootHist_t> fHist_x;                                     //!<!  x at dca (radial distance to vertex)
  Hist<RootHist_t> fHist_y;                                     //!<!  local Y-coordinate of track at dca  (cm)
  Hist<RootHist_t> fHist_z;                                     //!<!  local Z-coordinate of track at dca  (cm)
  Hist<RootHist_t> fHist_alpha;                                 //!<!  local to global angle
  Hist<RootHist_t> fHist_signed1Pt;                             //!<!  signed 1/pt (1/(GeV/c))
  Hist<RootHist_t> fHist_snp;                                   //!<!  local sine of the track momentum azimuthal angle
  Hist<RootHist_t> fHist_tgl;                                   //!<!  tangent of the track momentum dip angle
  Hist<RootHist_t> fHist_dcaxy;                                 //!<!  distance of closest approach in xy plane
  Hist<RootHist_t> fHist_dcaz;                                  //!<!  distance of closest approach in beam direction z
  Hist<RootHist_t> fHist_flag;                                  //!<!  flag info assigned to the track
  Hist<RootHist_t> fHist_pt;                                    //!<!  transverse momentum
  Hist<RootHist_t> fHist_eta;                                   //!<!  pseudorapidity
  Hist<RootHist_t> fHist_phi;                                   //!<!  azimuthal angle phi

  // its related properties
  Hist<RootHist_t> fHist_itsFoundClusters;                      //!<!  found clusters ITS
  Hist<RootHist_t> fHist_itsChi2PerCluster;                     //!<!  chi2 per cluster ITS
  Hist<RootHist_t> fHist_itsHits;                               //!<!  hitmap ITS

  // tpc related properties
  Hist<RootHist_t> fHist_tpcFindableClusters;                   //!<!  findable clusters TPC
  Hist<RootHist_t> fHist_tpcFoundClusters;                      //!<!  found clusters TPC
  Hist<RootHist_t> fHist_tpcSharedClusters;                     //!<!  shared clusters TPC
  Hist<RootHist_t> fHist_tpcFractionSharedClusters;             //!<!  fraction of shared clusters TPC
  Hist<RootHist_t> fHist_tpcCrossedRows;                        //!<!  crossed rows in TPC
  Hist<RootHist_t> fHist_tpcCrossedRowsOverFindableClusters;    //!<!  crossed rows over findable clusters in TPC
  Hist<RootHist_t> fHist_tpcChi2PerCluster;                     //!<!  chi2 per cluster TPC
  
  Hist<RootHist_t> fHist_tpcNClustersPID;                       //!<!  number of clusters used for PID in TPC
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCutStudies, 1);
  /// \endcond
};

#endif
