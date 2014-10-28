#ifndef AliAnalysisTaskEmcalJetHF_h
#define AliAnalysisTaskEmcalJetHF_h

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
//#include "AliCaloTrackESDReader.h"
//#include "AliCaloTrackAODReader.h"
//#include "AliCaloTrackReader.h"

// PID stuff                                                                                                                                
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliESDpid.h"
#include "AliAnalysisFilter.h"



class AliAnalysisTaskEmcalJetHF : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalJetHF();
  AliAnalysisTaskEmcalJetHF(const char *name);
  virtual ~AliAnalysisTaskEmcalJetHF();
  
  virtual void           UserCreateOutputObjects();
  
    /** Cuts info */
    AliAnalysisEtCuts * GetCuts() const { return fCuts; }
    virtual void SetCuts(const AliAnalysisEtCuts *cuts)
    { fCuts = (AliAnalysisEtCuts *) cuts; }
    
    
    void SetTPCITSTrackCuts(const AliESDtrackCuts *cuts){ fEsdtrackCutsITSTPC = (AliESDtrackCuts *) cuts;}
    void SetTPCOnlyTrackCuts(const AliESDtrackCuts *cuts){ fEsdtrackCutsTPC = (AliESDtrackCuts *) cuts;}
    void SetITSTrackCuts(const AliESDtrackCuts *cuts){ fEsdtrackCutsITS = (AliESDtrackCuts *) cuts;}
  
  //PID Sparse
  virtual THnSparse*      NewTHnSparseDHF(const char* name, UInt_t entries);
  virtual void            GetDimParamsHF(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  
  //JetQA Sparse
  virtual THnSparse*      NewTHnSparseDJetQA(const char* name, UInt_t entries);
  virtual void            GetDimParamsJetQA(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);


  // setters
  virtual void                  SetGlobalQA(Bool_t PID)                     { fGlobalQA = PID; }  //Global PID fpr all tracks and clusters in event
  void                          SetJetPt(Double_t jpt)                      { fJetHIpt = jpt; }  // jet threshold pt cut
  void                          SetTrackPtCut(Double_t trpt)                { fTrackPtCut = trpt; } // track pt threshold to do PID on
  virtual void                  SetTrackEta(Double_t e)                     { fTrackEta   = e; }  //eta range of the associated tracks
  virtual void                  SetTrackQACut(Double_t trkQAcut)            { fTrkQAcut = trkQAcut; }
  void                          SetTrackCuts(AliESDtrackCuts *cuts)         { fesdTrackCuts = cuts; }
  
  // eta and phi limits of jets - setters
  virtual void            SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void            SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }
  
  // event no.
  Int_t event;          // event number (processed)
  Int_t fillHist;
 protected:
  Bool_t                      Run();
  virtual void                Terminate(Option_t *);
  virtual Int_t               AcceptMyJet(AliEmcalJet *jet);   // applies basic jet tests/cuts before accepting
  virtual Int_t               GetCentBin(Double_t cent) const;  // Get Centrality bin
  void 			                  ExecOnce();
  void                        CheckClusTrackMatchingQA();
  
  // Trigger bit
  Bool_t           IsEventEMCALL1Gamma1()            const { return fEventTrigEMCALL1Gamma1  ; }
	Bool_t           IsEventEMCALL1Gamma2()            const { return fEventTrigEMCALL1Gamma2  ; }
  Bool_t           IsEventEMCALL1Gamma()             const { return (fEventTrigEMCALL1Gamma1 || fEventTrigEMCALL1Gamma2) ; }
  
  Bool_t           fEventTrigEMCALL1Gamma1 ;   // Event is L1-Gamma, threshold 1 on its name,  it should correspond kEMCEGA
	Bool_t           fEventTrigEMCALL1Gamma2 ;   // Event is L1-Gamma, threshold 2 on its name,  it should correspond kEMCEGA
  
  
  // data type switch
  Bool_t           fGlobalQA;
  AliVEvent      * fInputEvent;                //! pointer to esd or aod input
  
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
  Double_t              fTrkQAcut;                    //trkQA cut
  AliESDtrackCuts       *fesdTrackCuts;           //esd track cuts!!
  
  // PID                                                                                                                                    
  AliPIDResponse        *fPIDResponse;   // PID response object                                                                             
  AliTPCPIDResponse     *fTPCResponse;   // TPC pid response object
  
  AliESDtrackCuts* fEsdtrackCutsITSTPC;//esd track cuts for ITS+TPC tracks
  AliESDtrackCuts* fEsdtrackCutsTPC;//esd track cuts for TPC tracks (which may also contain ITS hits)
  AliESDtrackCuts* fEsdtrackCutsITS;//esd track cuts for ITS stand alone tracks
  
  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters
  AliParticleContainer       *fTracksJetCont;                 //! JetTracks
  AliClusterContainer        *fCaloClustersJetCont;           //!JetClusters
  

 private:
  AliESDEvent           *fESD;//!  // ESD object
  AliAODEvent           *fAOD;//!  // AOD Object

  
  TH2F                  *fHistRhovsCent; //!
  TH1F			            *fHistJetPhi;//!              // Njets vs Phi
  TH1F			            *fHistCorJetPt;//!            // (Njets) vs Corrected Jet Pt (local rho)
  TH1F			            *fHistJetPt;//!		            // raw jet pt (uncorrected)
  TH1F                  *fHistHighJetPt;//!

  TH2F                  *fHistJetPtvsTrackPt[6];//!
  TH1F                  *fHistTrackPt[6];//!
  TH1F                  *fHistEP0[6];//!
  TH1F                  *fHistEP0A[6];//!
  TH1F                  *fHistEP0C[6];//!
  TH2F                  *fHistEPAvsC[6];//!

  // PID status histo's                                                                                                                     
  TH2F                  *fHistnSigElecPt;//!  check
  TH2F                  *fHistnJetTrackvnJetClusters;//!  check
  TH3F                  *fHistPtDEtaDPhiTrackClus; //!track pt, delta eta, delta phi to matched cluster
  TH3F                  *fHistPtDEtaDPhiClusTrack; //! cluster pt, delta eta, delta phi to matched track
  
  //Sparse
  THnSparse             *fhnPIDHF;//!                  // PID sparse
  THnSparse             *fhnJetQA;//!                  //Jet QA Sparse
  THnSparse             *fhnClusterTrackQA;//!         // cluster QA sparse
  THnSparse             *fhnTrackClusterQA;//!         // track QA sparse
  THnSparse             *fhnPIDHFTtoC;//!              //Jet PID using track to Cluster matching
  //Declare it private to avoid compilation warning
  
  AliAnalysisTaskEmcalJetHF(const AliAnalysisTaskEmcalJetHF & g) ; // cpy ctor
  AliAnalysisTaskEmcalJetHF& operator=(const AliAnalysisTaskEmcalJetHF&); // not implemented
  


  ClassDef(AliAnalysisTaskEmcalJetHF, 4); // Emcal jet Heavy Flavor
};
#endif
