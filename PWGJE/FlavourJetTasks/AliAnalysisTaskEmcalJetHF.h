#ifndef AliAnalysisTaskEmcalJetHF_h
#define AliAnalysisTaskEmcalJetHF_h

//ROOT
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
class TClonesArray;
class TArrayI;
class TProfile;
//ALIROOT
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
class AliEMCALTriggerPatchInfo;
class AliVCaloTrigger;
//INCLUDES
#include <TRef.h>
#include <TBits.h>
#include <TMath.h>
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
#include "AliEMCALTriggerPatchInfo.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliEmcalTriggerSetupInfo.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliESDpid.h"
#include "AliAnalysisFilter.h"

class AliAnalysisTaskEmcalJetHF : public AliAnalysisTaskEmcalJet {
  public:
    
  enum MainPatchType {
    kManual = 0,    //just select highest energy patch in array
    kEmcalJet = 1   //use functionality of AliAnalysisTaskEmcal
  };
    
  AliAnalysisTaskEmcalJetHF();
  AliAnalysisTaskEmcalJetHF(const char *name);
  virtual                        ~AliAnalysisTaskEmcalJetHF();//!
  virtual void                   UserCreateOutputObjects();//!
  void                           SetMainPatchType(MainPatchType t)    { fMainPatchType      = t;}//!
  void                           SetMainTriggerTypeCat(TriggerCategory cat, Bool_t b) {fMainTrigCat = cat; fMainTrigSimple = b;}//!
  /** Cuts info **/
  AliAnalysisEtCuts * GetCuts()  const { return fCuts; }//!
  virtual void                   SetCuts(const AliAnalysisEtCuts *cuts){ fCuts = (AliAnalysisEtCuts *) cuts; }//!
  void                           SetTPCITSTrackCuts(const AliESDtrackCuts *cuts){ fEsdtrackCutsITSTPC = (AliESDtrackCuts *) cuts;}//!
  void                           SetTPCOnlyTrackCuts(const AliESDtrackCuts *cuts){ fEsdtrackCutsTPC = (AliESDtrackCuts *) cuts;}//!
  void                           SetITSTrackCuts(const AliESDtrackCuts *cuts){ fEsdtrackCutsITS = (AliESDtrackCuts *) cuts;}//!
  //PID Sparse
  virtual THnSparse*              NewTHnSparseDHF(const char* name, UInt_t entries);//!
  virtual void                    GetDimParamsHF(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);//!
  //JetTrigger Sparse
  virtual THnSparse*              NewTHnSparseDJetTrigger(const char* name, UInt_t entries);//!
  virtual void                    GetDimParamsJetTrigger(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);//
  // setters
  virtual void                    SetGlobalQA(Bool_t PID)                { fGlobalQA = PID; }//!               Global PID fpr all tracks and clusters in event
  void                            SetJetPt(Double_t jpt)                 { fJetHIpt = jpt; }//!                jet threshold pt cut
  void                            SetTrackPtCut(Double_t trpt)           { fTrackPtCut = trpt; }//!            track pt threshold to do PID on
  virtual void                    SetTrackEta(Double_t e)                { fTrackEta   = e; }//!               eta range of the associated tracks
  virtual void                    SetTrackQACut(Double_t trkQAcut)       { fTrkQAcut = trkQAcut; }//!
  void                            SetTrackCuts(AliESDtrackCuts *cuts)    { fesdTrackCuts = cuts; }//!
  virtual void                    SetFillHistograms(Bool_t fillhists)    { fFillHists = fillhists; }//!
  virtual void                    SetFillJetPID(Bool_t jetPID)           { fJetPID = jetPID; }//!              Jet PID
  // eta and phi limits of jets - setters
  virtual void                    SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }//!
  virtual void                    SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }//!
    
  protected:
  Bool_t                          Run();
  virtual void                    Terminate(Option_t *);//!                                                                Output message at end of analysis
  virtual Int_t                   AcceptMyJet(AliEmcalJet *jet);//!                                                        applies basic jet tests/cuts before accepting
  virtual Int_t                   AcceptJetforTag(AliVCluster *clusMatch, AliVTrack *AcceptedTrack);//!                    Accept Jet as HFE Jet to pass along to other tasks
  virtual Int_t                   GetCentBin(Double_t cent) const;//!    Get Centrality bin
  void 			                  ExecOnce();//!
  void                            CheckClusTrackMatchingQA();//!
  void                            ExtractMainPatch();//!
  // Trigger bit
  Bool_t                          IsEventEMCALL1Gamma1()                 const { return fEventTrigEMCALL1Gamma1  ; }//!
  Bool_t                          IsEventEMCALL1Gamma2()                 const { return fEventTrigEMCALL1Gamma2  ; }//!
  Bool_t                          IsEventEMCALL1Gamma()                  const { return (fEventTrigEMCALL1Gamma1 || fEventTrigEMCALL1Gamma2) ; }//!
  Bool_t                          fEventTrigEMCALL1Gamma1 ;//!           Event is L1-Gamma, threshold 1 on its name,  it should correspond kEMCEGA
  Bool_t                          fEventTrigEMCALL1Gamma2 ;//!           Event is L1-Gamma, threshold 2 on its name,  it should correspond kEMCEGA
  // data type switch
  AliVEvent                      *fInputEvent;//!                        pointer to esd or aod input
  AliAnalysisEtCuts              *fCuts;//!                              keeper of basic cut
  // cuts
  Double_t                       fPhimin;//!
  Double_t                       fPhimax;//!
  Double_t                       fEtamin;//!
  Double_t                       fEtamax;//!
  Double_t                       fAreacut;//!                            area cut on jet
  Double_t                       fJetHIpt;//!                            high jet pt
  Double_t                       fTrackPtCut;//!                         Track pt cut to do PID on
  Double_t                       fTrackEta;//!
  Double_t                       fTrkQAcut;//!                           Track QA cut
  Double_t                       fM02max;//!                             Cut on max EMCal Shower shape size for Electron ID
  Double_t                       fM02min;//!                             Cut on min EMCal Shower shape size for Electron ID
  AliESDtrackCuts                *fesdTrackCuts;//!                      ESD track cuts!!
  //General Task Switches
  Bool_t                         fGlobalQA;//!                           Switch for all HF Electron Candidates (Non-Jet)
  Bool_t                         fFillHists;//!                          Switch to Fill General Histograms
  Bool_t                         fJetPID;//!                             Switch for basic Jet HF Candidates
  // event no.
  Int_t event;//!                                                        event number (processed)
  // PID                                                                                                                                    
  AliPIDResponse                 *fPIDResponse;//!                       PID response object
  AliTPCPIDResponse              *fTPCResponse;//!                       TPC pid response object
  AliESDtrackCuts*               fEsdtrackCutsITSTPC;//!                 esd track cuts for ITS+TPC tracks
  AliESDtrackCuts*               fEsdtrackCutsTPC;//!                    esd track cuts for TPC tracks (which may also contain ITS hits)
  AliESDtrackCuts*               fEsdtrackCutsITS;//!                    esd track cuts for ITS stand alone tracks
  //Containers
  AliJetContainer                *fJetsCont;//!                          Jets
  AliParticleContainer           *fTracksCont;//!                        Tracks
  AliClusterContainer            *fCaloClustersCont;//!                  Clusters
  AliParticleContainer           *fTracksJetCont;//!                     JetTracks
  AliClusterContainer            *fCaloClustersJetCont;//!               JetClusters

  private:                                                               //Declare it private to avoid compilation warning
  AliESDEvent                   *fESD;//!                                ESD object
  AliAODEvent                   *fAOD;//!                                AOD Object
  AliEMCALTriggerPatchInfo      *fMaxPatch;//!                           main patch
  MainPatchType                 fMainPatchType;//!                       method to select main patch
  TriggerCategory               fMainTrigCat; //!                        trigger category for main trigger from AliAnalysisTaskEmcal::GetMainTriggerPatch
  Bool_t                        fMainTrigSimple;//!                      use offline trigger instead of online from AliAnalysisTaskEmcal::GetMainTriggerPatch
  //Histos
  TH1F                          *fHistMatchedClusJet;//!
  TH1F                          *fHistTriggerBitInfo;//!
  TH1F                          *fHistMaxTriggerBitInfo;//!
  TH1F                          *fHistEventSelection;//!
  TH2F                          *fHistRecalcGASize;//!                   Highest Energy GA Patch size per event (tower)
  TH1F                          *fHistRecalcGAEnergy;//!                 Highest Patch Energy per event
  TH2F                          *fHistCorrJetEvsPatchE;//!               Corrected Jet E v Trigger Patch E
  TH2F                          *fHistClusEvPatchE;//!                   Matched Cluster E vs Patch E
  TH2F                          *fHistdEtaPatchvdPhiPatch;//!
  TH2F                          *fHistRawJetEvPatchE;//!
  //Sparse
  THnSparse                     *fhnTriggerInfo;//!                      correlation between jets, patch energy and event observables
  THnSparse                     *fhnTrackClusterQA;//!                   track QA sparse
  THnSparse                     *fhnPIDHFTtoC;//!                        Jet PID using track to Cluster matching
  THnSparse                     *fhnJetTrigger;//!                       Jet Trigger Sparse
  
  AliAnalysisTaskEmcalJetHF(const AliAnalysisTaskEmcalJetHF & g) ;// cpy ctor
  AliAnalysisTaskEmcalJetHF& operator=(const AliAnalysisTaskEmcalJetHF&);// not implemented
  ClassDef(AliAnalysisTaskEmcalJetHF, 4);// Emcal jet Heavy Flavor
};
#endif
