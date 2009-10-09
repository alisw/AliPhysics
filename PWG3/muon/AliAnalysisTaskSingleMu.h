/// \ingroup "PWG3muon"
/// \class AliAnalysisTaskSingleMu
/// \brief Analysis task for single muons in the spectrometer
///
//  Author Diego Stocco

class TList;
class AliMCParticle;
class AluMCEvent;

class AliAnalysisTaskSingleMu : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSingleMu(const char *name = "AliAnalysisTaskSingleMu");
  virtual ~AliAnalysisTaskSingleMu();

  AliAnalysisTaskSingleMu(const AliAnalysisTaskSingleMu&);
  AliAnalysisTaskSingleMu& operator=(const AliAnalysisTaskSingleMu&);
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);
  
 private:

  void SetFlagsFromHandler ();

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
    kHistoPt,        ///< Pt distribution
    kHistoDCA,       ///< DCA distribution
    kHistoVz,        ///< Interaction vertex distribution
    kHistoEta,       ///< Eta distribution
    kHistoRapidity,  ///< Rapidity distribution
    // 2D histos
    kHistoPtDCA,     ///< DCA vs Pt distribution
    kHistoPtVz,      ///< Interaction vertex vs Pt distribution
    kHistoPtRapidity ///< Rapidity vs Pt distribution
  };

  // Histograms for MC
  enum {
    kHistoPtResolution, ///< Pt resolution
    kHistoPtDCAType,    ///< DCA vs Pt distribution
    kHistoPtVzType      ///< Interaction vertex vs Pt distribution
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

  Int_t GetHistoIndex(Int_t histoTypeIndex, Int_t trigIndex, Int_t srcIndex = -1);

  void FillTriggerHistos(Int_t histoIndex, Int_t matchTrig, Int_t motherType,
			 Float_t var1, Float_t var2 = 0. , Float_t var3 = 0.);
  
  Int_t RecoTrackMother(AliMCParticle* mcTrack, AliMCEvent* mcEvent);

  Bool_t fUseMC;  //!< Flag indicating that Monte Carlo information is available
  TList* fHistoList;   //!< List of histograms for data
  TList* fHistoListMC; //!< List of histograms for MC

  ClassDef(AliAnalysisTaskSingleMu, 1); // Single muon analysis
};

