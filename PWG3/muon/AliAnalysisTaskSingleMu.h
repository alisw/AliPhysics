/// \ingroup "PWG3muon"
/// \class AliAnalysisTaskSingleMu
/// \brief Analysis task for single muons in the spectrometer
///
//  Author Diego Stocco

class TList;
class AliMCParticle;
class TTree;
class TMap;
class TObjArray;
//class TAxis;
class AliCFManager;

class AliAnalysisTaskSingleMu : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSingleMu(const char *name = "AliAnalysisTaskSingleMu", Int_t fillTreeScaleDown = 0, Bool_t keepAll = kFALSE);
  virtual ~AliAnalysisTaskSingleMu();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);
  virtual void   NotifyRun();
  virtual void   FinishTaskOutput();

  // Please put the CINT1B class as first!
  void SetTriggerClasses(TString triggerClasses = 
			 "CINT1B-ABCE-NOPF-ALL CINT1A-ABCE-NOPF-ALL CINT1C-ABCE-NOPF-ALL CINT1-E-NOPF-ALL CMUS1B-ABCE-NOPF-MUON CMUS1A-ABCE-NOPF-MUON CMUS1C-ABCE-NOPF-MUON CMUS1-E-NOPF-MUON");

  /// Get CORRFW manager
  AliCFManager * GetCFManager() const { return fCFManager; }

  enum {
    kHvarPt,         ///< Pt at vertex
    kHvarY,          ///< Rapidity
    kHvarPhi,        ///< Phi
    kHvarDCA,        ///< DCA
    kHvarVz,         ///< Z vertex position
    kHvarThetaZones, ///< Theta at absorber end (4 zones)
    kHvarCharge,     ///< Particle charge
    kHvarMatchTrig,  ///< Matching trigger
    kHvarTrigClass,  ///< Trigger classes
    kHvarIsGoodVtx,  ///< IP vertex correctly reconstructed
    kHvarMotherType, ///< Mother type (MC only)
    kNvars           ///< THnSparse dimensions
  };  

  enum {
    kStepReconstructed,  ///< Reconstructed tracks
    kStepAcceptance,     ///< Track in acceptance
    kStepGeneratedMC,    ///< Generated tracks (MC)
    kStepAcceptanceMC,   ///< Track in acceptance (MC)
    kNsteps              ///< Number of steps
  };
  
 private:

  AliAnalysisTaskSingleMu(const AliAnalysisTaskSingleMu&);
  AliAnalysisTaskSingleMu& operator=(const AliAnalysisTaskSingleMu&);


  enum {
    kHistoNeventsPerTrig,    ///< Number of events per trigger
    kHistoMuonMultiplicity,  ///< Number of muons per event
    kHistoEventVz,           ///< Vertex z distribution for all events
    kHistoEventVzNorm,       ///< Vertex z distribution normalized to the muon events
    kHistoNeventsPerRun,     ///< Number of triggers per run (for check)
    kHistoNmuonsPerRun,      ///< Number of muons per run (for check)
    kNsummaryHistos          ///< Number of summary histograms
  };

  enum {
    kNoMatchTrig,  ///< No match with trigger
    kAllPtTrig,    ///< Match All Pt
    kLowPtTrig,    ///< Match Low Pt
    kHighPtTrig,   ///< Match High Pt
    kNtrigCuts     ///< Total number of trigger types
  };


  // Histograms for MC
  enum {
    kHistoCheckVzMC,    ///< Check vertex distribution
    kNsummaryHistosMC   ///< Summary histograms for MC
  };

  enum {
    kHistoPtResolutionMC, ///< Pt resolution
    kNhistoTypesMC        ///< Number of MC histograms
  };

  enum {
    kCharmMu,       ///< Mu from charm
    kBeautyMu,      ///< Mu from beauty
    kPrimaryMu,     ///< Primary mu
    kSecondaryMu,   ///< Secondary mu
    kRecoHadron,    ///< Reconstructed hadron
    kUnknownPart,   ///< Particle that fails matching kine
    kNtrackSources  ///< Total number of track sources
  };

  // Trees
  enum {
    kVarPx, ///< Px at vertex
    kVarPy, ///< Py at vertex
    kVarPz, ///< Pz at vertex
    kVarPt, ///< Pt at vertex
    kVarPxAtDCA, ///< Px at DCA
    kVarPyAtDCA, ///< Py at DCA
    kVarPzAtDCA, ///< Pz at DCA
    kVarPtAtDCA, ///< Pt at DCA
    kVarPxUncorrected, ///< Px at first chamber
    kVarPyUncorrected, ///< Py at first chamber
    kVarPzUncorrected, ///< Pz at first chamber
    kVarPtUncorrected, ///< Pt at first chamber
    kVarXUncorrected, ///< X at first chamber
    kVarYUncorrected, ///< Y at at first chamber
    kVarZUncorrected, ///< Z at at first chamber
    kVarXatDCA, ///< X position at DCA
    kVarYatDCA, ///< Y position at DCA
    kVarDCA, ///< DCA
    kVarEta, ///< Eta
    kVarRapidity, ///< Rapidity
    kVarCharge, ///< Charge
    kVarRAtAbsEnd, ///< R at absorber end
    // Global event info
    kVarIPVx, ///< IP x position
    kVarIPVy, ///< IP y position
    kVarIPVz, ///< IP z position
    kNvarFloat
  };

  enum {
    kVarMatchTrig, ///< Match trigger
    kVarIsMuon, ///< Is muon
    kVarIsGhost, ///< Is Ghost (trigger track not matching tracker)
    kVarLocalCircuit, ///< Fired local circuit
    // Global event info
    kVarPassPhysicsSelection, ///< Pass physics selection (Requires ESD)
    kVarNVtxContrib, ///< Vertex contributors
    kVarNspdTracklets, ///< SPD tracklets
    kVarIsPileup, ///< Pileup vertices
    kNvarInt
  };

  enum {
    // Global event info
    kVarTrigMask, ///< Fires triggers mask
    kNvarChar
  };

  // Event ID
  enum {
    // Global event info
    kVarBunchCrossNumber, ///< Bunch crossing number
    kVarOrbitNumber, ///< Orbit number
    kVarPeriodNumber, ///< Period number
    kVarRunNumber, ///< Run number
    kNvarUInt
  };
  
  // Trees MC
  enum {
    kVarPxMC, ///< Px from Kine
    kVarPyMC, ///< Py from Kine
    kVarPzMC, ///< Pz from Kine
    kVarPtMC, ///< Pt from Kine
    kVarEtaMC,  ///< Eta from Kine
    kVarRapidityMC, ///< Rapidity from Kine
    kVarVxMC, ///< Particle production x vertex from Kine
    kVarVyMC, ///< Particle production y vertex from Kine
    kVarVzMC, ///< Particle production z vertex from Kine
    kNvarFloatMC
  };

  enum {
    kVarPdg, ///< PDG
    kVarMotherType, ///< Mother type
    kNvarIntMC
  };

  Int_t GetHistoIndex(Int_t histoTypeIndex, Int_t trigIndex = -1, Int_t srcIndex = -1);
  Float_t GetBinThetaAbsEnd(Float_t RAtAbsEnd, Bool_t isTheta = kFALSE);
  Float_t GetBinTrigClass(const Char_t* trigClass);

  void FillTriggerHistos(Int_t histoIndex, Int_t matchTrig, Int_t motherType,
			 Float_t var1, Float_t var2 = 0. , Float_t var3 = 0.);
  
  Int_t RecoTrackMother(AliMCParticle* mcParticle);

  void Reset(Bool_t keepGlobal = kTRUE);

  void SetAxisLabel(TAxis* axis);

  Int_t fFillTreeScaleDown; ///< Ntuple must be filled each fFillTreeScaleDown events
  Bool_t fKeepAll; ///< Flag indicating to keep all info in tree
  Int_t  fkNvtxContribCut; ///< Number of cuts in vertex contributors
  AliCFManager* fCFManager; //!< Pointer to the CF manager
  TList* fHistoList;   //!< List of histograms for data
  TList* fHistoListMC; //!< List of histograms for MC
  TList* fHistoListQA; //!< List of QA histos
  TTree* fTreeSingleMu; //!< Optional output Tree
  TTree* fTreeSingleMuMC; //!< Optional output Tree for MC
  Float_t* fVarFloat; //!< Reconstructed parameters float
  Int_t* fVarInt; //!< Reconstructed parameters int
  Char_t** fVarChar; //!< Reconstructed parameters string
  UInt_t* fVarUInt; //!< Reconstructed parameters Uint
  Float_t* fVarFloatMC; //!< MC parameters float
  Int_t* fVarIntMC; //!< MC parameters int
  TMap* fAuxObjects; //!< Map of vertex distribution per run
  TObjArray* fTriggerClasses; //!< full trigger class name

  ClassDef(AliAnalysisTaskSingleMu, 2); // Single muon analysis
};

