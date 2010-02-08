/// \ingroup "PWG3muon"
/// \class AliAnalysisTaskSingleMu
/// \brief Analysis task for single muons in the spectrometer
///
//  Author Diego Stocco

class TList;
class AliMCParticle;
class AliMCEvent;
class TTree;

class AliAnalysisTaskSingleMu : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSingleMu(const char *name = "AliAnalysisTaskSingleMu", Bool_t fillTree = kFALSE, Bool_t keepAll = kFALSE);
  virtual ~AliAnalysisTaskSingleMu();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);
  virtual void   NotifyRun();
  
 private:

  AliAnalysisTaskSingleMu(const AliAnalysisTaskSingleMu&);
  AliAnalysisTaskSingleMu& operator=(const AliAnalysisTaskSingleMu&);

  enum {
    kNoMatchTrig,  // No match with trigger
    kAllPtTrig,    // Match All Pt
    kLowPtTrig,    // Match Low Pt
    kHighPtTrig,   // Match High Pt
    kNtrigCuts     // Total number of trigger types
  };
  
  // Histograms for data
  enum {
    // 1D histos
    kHistoPt,         ///< Pt distribution
    kHistoDCA,        ///< DCA distribution
    kHistoVz,         ///< Interaction vertex distribution
    kHistoEta,        ///< Eta distribution
    kHistoRapidity,   ///< Rapidity distribution
    // 2D histos
    kHistoPtDCA,      ///< DCA vs Pt distribution
    kHistoPtVz,       ///< Interaction vertex vs Pt distribution
    kHistoPtRapidity, ///< Rapidity vs Pt distribution
    kNhistoTypes
  };

  enum {
    kHistoMuonMultiplicity,   ///< Number of muons per event
    kNsummaryHistos
  };

  // Histograms for MC
  enum {
    kHistoPtResolutionMC, ///< Pt resolution
    kHistoPtDCAMC,        ///< DCA vs Pt distribution
    kHistoPtVzMC,         ///< Interaction vertex vs Pt distribution
    kNhistoTypesMC
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
    // Global event info
    kVarPassPhysicsSelection, ///< Pass physics selection (Requires ESD)
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
    //kVarIPVxMC, ///< IP x position from Kine
    //kVarIPVyMC, ///< IP y position from Kine
    //kVarIPVzMC, ///< IP z position from Kine
    kNvarFloatMC
  };

  enum {
    kVarPdg, ///< PDG
    kVarMotherType, ///< Mother type
    kNvarIntMC
  };

  Int_t GetHistoIndex(Int_t histoTypeIndex, Int_t trigIndex, Int_t srcIndex = -1);

  void FillTriggerHistos(Int_t histoIndex, Int_t matchTrig, Int_t motherType,
			 Float_t var1, Float_t var2 = 0. , Float_t var3 = 0.);
  
  Int_t RecoTrackMother(AliMCParticle* mcTrack, AliMCEvent* mcEvent);

  void Reset(Bool_t keepGlobal = kTRUE);

  Bool_t fUseMC;  ///< Flag indicating that Monte Carlo information is available
  Bool_t fFillTree; ///< Flag indicating that ntuple must be filled
  Bool_t fKeepAll; ///< Flag indicating to keep all info in tree
  TList* fHistoList;   //!< List of histograms for data
  TList* fHistoListMC; //!< List of histograms for MC
  TTree* fTreeSingleMu; //!< Optional output Tree
  TTree* fTreeSingleMuMC; //!< Optional output Tree for MC
  Float_t* fVarFloat; //!< Reconstructed parameters float
  Int_t* fVarInt; //!< Reconstructed parameters int
  Char_t** fVarChar; //!< Reconstructed parameters string
  UInt_t* fVarUInt; //!< Reconstructed parameters Uint
  Float_t* fVarFloatMC; //!< MC parameters float
  Int_t* fVarIntMC; //!< MC parameters int

  ClassDef(AliAnalysisTaskSingleMu, 2); // Single muon analysis
};

