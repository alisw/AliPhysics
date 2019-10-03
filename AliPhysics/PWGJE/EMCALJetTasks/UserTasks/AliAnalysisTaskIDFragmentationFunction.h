// *************************************************************************
// * Task for Fragmentation Function Analysis in PWG4 Jet Task Force Train *
// *************************************************************************

#ifndef ALIANALYSISTASKIDFRAGMENTATIONFUNCTION_H
#define ALIANALYSISTASKIDFRAGMENTATIONFUNCTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class AliESDEvent;
class AliAODEvent;
class AliEmcalJet;
class AliAODExtension;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class THnSparse; 
class TRandom3;
class TArrayS;
class AliAnalysisUtils;
class AliAODTrack;
class AliAODMCParticle;
class AliMCParticleContainer;
class AliAnalysisTaskPID;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliPID.h"
  
class AliAnalysisTaskIDFragmentationFunction : public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskIDFragmentationFunction(); 
  AliAnalysisTaskIDFragmentationFunction(const char *name);
  virtual ~AliAnalysisTaskIDFragmentationFunction();
  
  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   LocalInit() {Init();}

  Bool_t         FillHistograms();
  virtual void   Terminate(Option_t* );
  virtual Bool_t Notify();

  virtual void   SetNonStdFile(char* c){fNonStdFile = c;} 
  
  virtual void   SetCentralityEstimator(TString c){fCentralityEstimator = c;}
  
  virtual void   SetNameTrackContainer(TString c) {fNameTrackContainer = c;}
  virtual TString GetNameTrackContainer() const { return fNameTrackContainer;}
  
  virtual void   SetNameTrackContainerEfficiency(TString c) {fNameTrackContainerEfficiency = c;}
  virtual TString GetNameTrackContainerEfficiency() const { return fNameTrackContainerEfficiency;}

  virtual void   SetNameMCParticleContainer(TString c) {fNameMCParticleContainer = c;}
  virtual TString GetNameMCParticleContainer() const { return fNameMCParticleContainer;}
  
  virtual void   SetNameJetContainer(TString c) {fNameJetContainer = c;}
  virtual TString GetNameJetContainer() const { return fNameJetContainer;}
  
  virtual void   SetNameMCParticleJetContainer(TString c) {fNameMCParticleJetContainer = c;}
  virtual TString GetNameMCParticleJetContainer() const { return fNameMCParticleJetContainer;}

  virtual void   UseAODInputJets(Bool_t b) {fUseAODInputJets = b;}  
  virtual void   UsePhysicsSelection(Bool_t b) {fUsePhysicsSelection = b;}
  virtual void   SetEventSelectionMask(UInt_t i){fEvtSelectionMask = i;}
  virtual void   SetEventClass(Int_t i){fEventClass = i;}
  virtual void   SetMaxVertexZ(Float_t z){fMaxVertexZ = z;}

  virtual void   SetFFRadius(Float_t r = 0.4) { fFFRadius = r; }
  virtual void   SetMCPtHardCut(Float_t ptHardCut)      { fMCPtHardCut = ptHardCut; }
  
  virtual Bool_t GetStudyTransversalJetStructure() const { return fStudyTransversalJetStructure; }
  virtual void SetStudyTransversalJetStructure(Bool_t flag) { fStudyTransversalJetStructure = flag; }
  
  virtual Bool_t GetFillDCA() const { return fFillDCA; }
  virtual void SetFillDCA(Bool_t flag) { fFillDCA = flag; }  
  
  virtual Bool_t GetUseInclusivePIDtask() const {return fUseInclusivePIDtask; }
  virtual void SetUseInclusivePIDtask(Bool_t flag) {fUseInclusivePIDtask = flag; }
  
  virtual Bool_t GetUseJetPIDtask() const {return fUseJetPIDtask; }
  virtual void SetUseJetPIDtask(Bool_t flag) {fUseJetPIDtask = flag; }
  
  virtual Bool_t GetUseJetUEPIDtask() const {return fUseJetUEPIDtask; }
  virtual void SetUseJetUEPIDtask(Bool_t flag) {fUseJetUEPIDtask = flag; }
  
  //Begin of underlying event calculations
  virtual TList* GetUEJetsWithRandomConeMethod(AliJetContainer* jetContainer, Double_t coneRadius, Double_t maxEtaTrack);
  virtual TList* GetUEJetsWithPerpendicularConeMethod(AliJetContainer* jetContainer);
  virtual AliEmcalJet* GetRandomCone(AliEmcalJet* processedJet, Double_t dEtaConeMax, Double_t dDistance) const;
  virtual AliEmcalJet* GetPerpendicularCone(AliEmcalJet* processedJet, Double_t perpAngle) const;
  virtual TList* GetTracksInCone(const AliEmcalJet* jet, AliParticleContainer* particleContainer = 0x0) const;
  Bool_t IsParticleInCone(const AliVParticle* part1, const AliVParticle* part2, Double_t dRMax) const;
  virtual Bool_t IsRCJCOverlap(const AliVParticle* part, Double_t dDistance) const; 
  
  //Filling the efficiency containers of the specified task with the generated jet yield
  virtual void PerformJetMonteCarloAnalysisGeneratedYield(AliEmcalJet* jet, AliVParticle* trackVP, AliAnalysisTaskPID* task, Double_t centPercent, AliJetContainer* mcJetContainer = 0x0);
  
  //Jet Track Calculations, including filling of the efficiency containers. If no task is specified, the function loops over all Jet tasks in fJetPIDtask, using trackRejectedByTask[] to decide if the track is accepted. If a task is specified, everything is done (without checking further) for the specified task. 
  virtual void AnalyseJetTrack(AliVTrack* track, AliEmcalJet* jet, AliAnalysisTaskPID* task, const Bool_t* trackRejectedByTask, Double_t centPercent, AliMCParticleContainer* mcParticleContainer = 0x0);
  
  //Fill DCA
  virtual void FillDCA(AliVTrack* track, AliMCParticleContainer* mcParticleContainer); 
  
  //End of underlying event calculations

  static  void   SetProperties(TH1* h,const char* x, const char* y);
  static  void   SetProperties(TH1* h,const char* x, const char* y,const char* z);
  static  void   SetProperties(THnSparse* h,Int_t dim, const char** labels);
  
  void SetStoreXi(Bool_t storeXi) { fStoreXi = storeXi; }
  Bool_t GetStoreXi() { return fStoreXi; }
  
  const TString  GetCentralityEstimator() const { return fCentralityEstimator; }
  Float_t  GetFFRadius() const { return fFFRadius; }
  Float_t  GetFFMinLTrackPt() const { return fFFMinLTrackPt; }
  Float_t  GetFFMaxTrackPt() const { return fFFMaxTrackPt; }
  Float_t  GetFFMinNTracks() const { return fFFMinnTracks; }
  Float_t  GetMCPtHardCut() const  { return fMCPtHardCut; }
  
  Double_t GetDistanceJetTrack(const AliEmcalJet* jet, const AliVParticle* track, Bool_t useInternalFlag = kTRUE) const;
  
  Double_t GetMCStrangenessFactor(Double_t pt) const;
  Double_t GetMCStrangenessFactorCMS(AliAODMCParticle* daughter, AliMCParticleContainer* mcParticleContainer) const;
  
  Double_t GetPerpendicularMomentumTrackJet(const AliEmcalJet* jet, const AliVParticle* track, Bool_t useInternalFlag = kTRUE) const;
  
  Bool_t IsSecondaryWithStrangeMotherMC(AliAODMCParticle* part, AliMCParticleContainer* mcParticleContainer);

  const TString* GetNamesOfInclusivePIDtasks() const { return fNameInclusivePIDtask; };
  void SetNamesOfInclusivePIDtasks(Int_t numNames, const TString* names);
  
  const TString* GetNamesOfJetPIDtasks() const { return fNameJetPIDtask; };
  void SetNamesOfJetPIDtasks(Int_t numNames, const TString* names);
  
  const TString* GetNamesOfJetUEPIDtasks() const { return fNameJetUEPIDtask; };
  void SetNamesOfJetUEPIDtasks(Int_t numNames, const TString* names, const TString* methods = 0x0);
	
	
  Bool_t GetIsPP() const { return fIsPP; };
  void SetIsPP(Bool_t flag) { fIsPP = flag; };
  
  UInt_t GetRCTrials() const { return fRCTrials; }
  void SetRCTrials(UInt_t trials) {fRCTrials = trials; }

  const TString* GetUEMethods() const { return fUEMethods; }
  void SetUEMethods(const TString* names);
  
  Bool_t GetOnlyLeadingJets() const { return fOnlyLeadingJets; }
  void SetOnlyLeadingJets(Bool_t onlyLeadingJets) { fOnlyLeadingJets = onlyLeadingJets; }
  
  void RemoveParticleContainer(const char* n) {fParticleCollArray.Remove(GetParticleContainer(n));}
  void RemoveJetContainer(const char* n) {fJetCollArray.Remove(GetJetContainer(n));}

 
 protected:

  AliESDEvent* fESD;      //! ESD event
  AliAODEvent* fAOD;      //! AOD event
  AliAODEvent* fAODJets;  //! AOD event with jet branch (case we have AOD both in input and output)
  AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD
  TString       fNonStdFile; // name of delta aod file to catch the extension
 
  TString fCentralityEstimator;   // Estimator for the Centrality, V0M is default, set to be V0A for pPb-collisions    
  
  TString fNameTrackContainer;
  TString fNameTrackContainerEfficiency;
  TString fNameMCParticleContainer;
  TString fNameJetContainer;
  TString fNameMCParticleJetContainer;

  Bool_t  fUseAODInputJets;     // take jets from in/output - only relevant if AOD event both in input AND output and we want to use output
  Bool_t  fUsePhysicsSelection; // switch for event selection
  UInt_t  fEvtSelectionMask;    // trigger class selection
  Int_t   fEventClass;          // centrality class selection
  Float_t fMaxVertexZ;          // maximum abs(z) position of primiary vertex [cm]

  Float_t fFFRadius;        // if radius > 0 construct FF from tracks within cone around jet axis, otherwise use trackRefs  
  Float_t fFFMinLTrackPt;   // reject jets with leading track with pt smaller than this value
  Float_t fFFMaxTrackPt;    // reject jets containing any track with pt larger than this value
  Int_t   fFFMinnTracks;    // reject jets with less tracks than this value  

  Float_t fAvgTrials;       // average number of trials per event
  
  // histogram bins  

  Bool_t fStoreXi;          // Store Xi = -ln(z) explicitly
  Bool_t fStudyTransversalJetStructure; // Store observables related to transversal jet structure
  
  // Histograms
  TList	*fCommonHistList;         // List of common histos
  
  TH1F  *fh1EvtSelection;         //! event cuts 
  TH1F  *fh1VtxSelection;         //! type of accepted vertices
  TH1F	*fh1VertexNContributors;  //! NContributors to prim vertex
  TH1F	*fh1VertexZ;              //! prim vertex z distribution
  TH1F	*fh1EvtMult;              //! number of reconstructed tracks after cuts 
  TH1F	*fh1EvtCent;              //! centrality percentile 

  TProfile* fh1Xsec;              //! pythia cross section and trials
  TH1F*     fh1Trials;            //! sum of trials
  TH1F*     fh1PtHard;            //! pt hard of the event
  TH1F*     fh1PtHardTrials;      //! pt hard of the event
  
  TH1F*     fh1EvtsPtHardCut;     //! Number events before and after the cut on MC pT hard

  TH1F  *fh1nRecJetsCuts;         //! number of jets from reconstructed tracks per event 
  TH1F  *fh1nGenJets;             //! number of jets from generated tracks per event
  
  TH2F  *fhDCA_XY;                //! DCA XY for all rec. particles
  TH2F  *fhDCA_Z;                 //! DCA Z for all rec. particles
  
  TH2F  *fhJetPtRefMultEta5;      //! Jet pT vs. reference multiplicity (|eta|<0.5)
  TH2F  *fhJetPtRefMultEta8;      //! Jet pT vs. reference multiplicity (|eta|<0.8)
  TH2F  *fhJetPtMultPercent;      //! Jet pT vs. multiplicity percentile (usually V0M)

  TH2F  *fhDCA_XY_prim_MCID[AliPID::kSPECIES];   //! DCA XY for all rec. prim. particles sorted by MC-ID
  TH2F  *fhDCA_Z_prim_MCID[AliPID::kSPECIES];    //! DCA Z for all rec. prim. particles sorted by MC-ID
 
  TH2F  *fhDCA_XY_sec_MCID[AliPID::kSPECIES];    //! DCA XY for all rec. sec. particles sorted by MC-ID
  TH2F  *fhDCA_Z_sec_MCID[AliPID::kSPECIES];     //! DCA Z for all rec. sec. particles sorted by MC-ID

  TRandom3* fRandom;                        //! TRandom3 for background estimation 
  
  Bool_t fOnlyLeadingJets;                  // Flag indicating whether some histos are filled with leading jets only or all jets
  Float_t fMCPtHardCut;                     // Cut on MC pThard (smaller that threshold), if set to non-negative value
  
  AliAnalysisUtils *fAnaUtils;              //! Object to use analysis utils like pile-up rejection
  
  // PID framework
  Int_t fNumInclusivePIDtasks;              // Number of inclusive PID tasks used 
  Int_t fNumJetPIDtasks;                    // Number of jet PID tasks used
  Int_t fNumJetUEPIDtasks;                  // Number of jet UE PID tasks used
  
  TString* fNameInclusivePIDtask;           //[fNumInclusivePIDtasks] Names of the tasks for inclusive PID spectra
  TString* fNameJetPIDtask;                 //[fNumJetPIDtasks] Names of the tasks for jet PID spectra
  TString* fNameJetUEPIDtask;               //[fNumJetUEPIDtasks] Names of the tasks for jet UE PID spectra
  
  AliAnalysisTaskPID** fInclusivePIDtask;   //! Pointer to tasks for inclusive PID spectra
  AliAnalysisTaskPID** fJetPIDtask;         //! Pointer to tasks for jet PID spectra
  AliAnalysisTaskPID** fJetUEPIDtask;       //! Pointer to tasks for jet UE PID spectra
  
  Bool_t fUseInclusivePIDtask;              // Process inclusive PID spectra?
  Bool_t fUseJetPIDtask;                    // Process jet PID spectra?
  Bool_t fUseJetUEPIDtask;                  // Process jet UE PID spectra?
  
  Bool_t fIsPP;                             // Is pp collision system? -> If yes, centrality will be set to -1
  
  Bool_t fFillDCA;                          //Shall the DCA histograms be filled?
  
  //Underlying event members
  UInt_t fRCTrials;
  TString* fUEMethods;                      //[fNumJetUEPIDtasks] Names for the underlying event methods
  
private:
  AliAnalysisTaskIDFragmentationFunction(const  AliAnalysisTaskIDFragmentationFunction&);   //Not implemented in AliAnalysisTaskEmcalJet
  AliAnalysisTaskIDFragmentationFunction& operator=(const  AliAnalysisTaskIDFragmentationFunction);   //Not implemented AliAnalysisTaskEmcalJet
  
  ClassDef(AliAnalysisTaskIDFragmentationFunction, 23);
};


inline void AliAnalysisTaskIDFragmentationFunction::SetNamesOfInclusivePIDtasks(Int_t numNames, const TString* names)
{
  delete [] fNameInclusivePIDtask;
  fNameInclusivePIDtask = 0x0;
  
  if (!names || numNames < 0) {
    fNumInclusivePIDtasks = 0;
    return;
  }
  
  fNumInclusivePIDtasks = numNames;
  
  if (numNames > 0) {
    fNameInclusivePIDtask = new TString[numNames];
    SetUseInclusivePIDtask(kTRUE);
    
    for (Int_t i = 0; i < numNames; i++) {
      fNameInclusivePIDtask[i] = names[i];
    }
  }  
}

inline void AliAnalysisTaskIDFragmentationFunction::SetNamesOfJetPIDtasks(Int_t numNames, const TString* names)
{
  delete [] fNameJetPIDtask;
  fNameJetPIDtask = 0x0;
  
  if (!names || numNames < 0) {
    fNumJetPIDtasks = 0;
    return;
  }
  
  fNumJetPIDtasks = numNames;
  
  if (numNames > 0) {
    fNameJetPIDtask = new TString[numNames];
    SetUseJetPIDtask(kTRUE);
    
    for (Int_t i = 0; i < numNames; i++) {
      fNameJetPIDtask[i] = names[i];
    }
  }  
}

inline void AliAnalysisTaskIDFragmentationFunction::SetNamesOfJetUEPIDtasks(Int_t numNames, const TString* names, const TString* methods)
{
  delete [] fNameJetUEPIDtask;
  fNameJetUEPIDtask = 0x0;
  
  if (!names || numNames < 0) {
    fNumJetUEPIDtasks = 0;
    return;
  }
  
  fNumJetUEPIDtasks = numNames;
  
  if (numNames > 0) {
    fNameJetUEPIDtask = new TString[numNames];
    SetUseJetUEPIDtask(kTRUE);
    
    for (Int_t i = 0; i < numNames; i++) {
      fNameJetUEPIDtask[i] = names[i];
    }

    if(methods) {
      SetUEMethods(methods);
    }
    else {
      TString* defaultMethods = 0x0;
      switch(numNames) {
        case 1:
          defaultMethods = new TString[1];
          defaultMethods[1] = "RC";
          SetUEMethods(defaultMethods);
          break;
        case 2:
          defaultMethods = new TString[2];
          defaultMethods[1] = "RC";
          defaultMethods[2] = "PC";
          SetUEMethods(defaultMethods);
          break;
        default:
          std::cout << "Too many Underlying event PID tasks to specify default methods. Please specify them yourselves for each task by giving TString array as third argument in SetNamesOfJetUEPIDtasks!" << std::endl;
          break;
      }
      delete defaultMethods;
      defaultMethods = 0x0;
    }
  }
}

inline void AliAnalysisTaskIDFragmentationFunction::SetUEMethods(const TString* names) {
  if (fNumJetUEPIDtasks < 1) {
    std::cout << "Number of Underlying event PID tasks has to be greater than 0 to set the methods" << std::endl;
    return;
  }
  delete [] fUEMethods;
  fUEMethods = 0x0;
  if (!names) {
    return;
  }
  
  fUEMethods = new TString[fNumJetUEPIDtasks];
  
  for (Int_t i = 0; i < fNumJetUEPIDtasks; ++i) {
    fUEMethods[i] = names[i];
  }
}

#endif
