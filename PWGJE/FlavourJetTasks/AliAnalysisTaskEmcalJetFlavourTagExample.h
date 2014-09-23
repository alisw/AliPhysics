#ifndef AliAnalysisTaskEmcalJetFlavourTagExample_h
#define AliAnalysisTaskEmcalJetFlavourTagExample_h

class TClonesArray;
class TH1;
class TH2;
class TH3;
class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TLorentzVector;
class TGraph;

class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliESDtrack;
class AliESDtrackCuts;
class AliAnalysisEtCuts;
class AliDetectorPID;
class AliESDCaloCluster;

class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;


// this whole section of includes added 
#include <AliEmcalJet.h>
#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom3.h>
#include <AliLog.h>
#include "AliAnalysisTaskEmcalJet.h"
#include "AliEMCALPIDResponse.h"
#include <AliESDCaloCluster.h>
#include <AliESDtrackCuts.h>
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliTPCdEdxInfo.h"

// PID stuff                                                                                                                                
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliESDpid.h"
#include "AliAnalysisFilter.h"

class AliAnalysisTaskEmcalJetFlavourTagExample : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalJetFlavourTagExample();
  AliAnalysisTaskEmcalJetFlavourTagExample(const char *name);
  virtual ~AliAnalysisTaskEmcalJetFlavourTagExample();
  
  virtual void           UserCreateOutputObjects();
  
    /** Cuts info */
    AliAnalysisEtCuts * GetCuts() const { return fCuts; }
    virtual void SetCuts(const AliAnalysisEtCuts *cuts)
    { fCuts = (AliAnalysisEtCuts *) cuts; }

  // setters
  void                    SetJetPt(Double_t jpt)           { fJetHIpt = jpt; }  // jet threshold pt cut
  void                    SetTrackPtCut(Double_t trpt)     { fTrackPtCut = trpt; } // track pt threshold to do PID on
  virtual void            SetTrackEta(Double_t e)                 { fTrackEta   = e; }  //eta range of the associated tracks
  void                    SetTrackCuts(AliESDtrackCuts *cuts)         { fesdTrackCuts = cuts; }
  
  // eta and phi limits of jets - setters
  virtual void            SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void            SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }
  
  // event no.
  Int_t event;          // event number (processed)
  
 protected:
  Bool_t                      Run();
  virtual void                Terminate(Option_t *);
  virtual Int_t               AcceptMyJet(AliEmcalJet *jet);   // applies basic jet tests/cuts before accepting
  void 			                  ExecOnce();

  
  // data type switch
  
  AliAnalysisEtCuts *fCuts;                            // keeper of basic cuts

  // cuts
  Double_t              fPhimin;
  Double_t              fPhimax;
  Double_t              fEtamin;
  Double_t              fEtamax;
  Double_t              fAreacut;                     //area cut
  Double_t              fJetHIpt;                    // high jet pt 
  Double_t              fTrackPtCut;                 // track pt cut to do PID on
  Double_t              fTrackEta;  
  AliESDtrackCuts       *fesdTrackCuts;           //esd track cuts!!
  
  // PID                                                                                                                                    
  AliPIDResponse        *fPIDResponse;   // PID response object                                                                             
  AliTPCPIDResponse     *fTPCResponse;   // TPC pid response object
  
  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters
  AliParticleContainer       *fTracksJetCont;                 //! JetTracks
  AliClusterContainer        *fCaloClustersJetCont;           //!JetClusters
  

 private:
  AliESDEvent               *fESD;//!  // ESD object
  AliAODEvent               *fAOD;//!  // AOD Object

  
  TH1F			                *fHistJetPhi;//!              // Njets vs Phi
  TH1F			                *fHistCorJetPt;//!            // (Njets) vs Corrected Jet Pt (local rho)
  TH1F			                *fHistJetPt;//!		            // raw jet pt (uncorrected)
  TH1F                      *fHistHighJetPt;//!

  TH2F                      *fHistnSigElecPt;//!  check
  TH2F                      *fHistdEdXvsPt;//! check
  TH2F                      *fHistnJetTrackvnJetClusters;//!  check

  AliAnalysisTaskEmcalJetFlavourTagExample(const AliAnalysisTaskEmcalJetFlavourTagExample & g) ; // cpy ctor
  AliAnalysisTaskEmcalJetFlavourTagExample& operator=(const AliAnalysisTaskEmcalJetFlavourTagExample&); // not implemented
  


  ClassDef(AliAnalysisTaskEmcalJetFlavourTagExample, 4); // Emcal jet Heavy Flavor
};
#endif
