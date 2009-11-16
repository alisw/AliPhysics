/// \ingroup "PWG3muon"
/// \class AliAnalysisTaskTrigChEff
/// \brief Analysis task for trigger chamber efficiency determination
///
//  Author Diego Stocco

class TList;

class AliAnalysisTaskTrigChEff : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTrigChEff(const char *name = "AliAnalysisTaskTrigChEff");
  virtual ~AliAnalysisTaskTrigChEff() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  /// Use ghost tracks in calculations
  void SetUseGhostTracks(Bool_t useGhosts = kTRUE) { fUseGhosts = useGhosts; }

protected:
  void ResetHistos();

private:
  /// Not implemented
  AliAnalysisTaskTrigChEff(const AliAnalysisTaskTrigChEff& rhs);
  /// Not implemented
  AliAnalysisTaskTrigChEff& operator = (const AliAnalysisTaskTrigChEff& rhs);
    
  Bool_t fUseGhosts; ///< Flag to use also the trigger tracks not matching the tracker in eff. calculation

  TList*  fList; ///<TList output object

  enum {
    kNcathodes = 2,  ///< Number of cathodes
    kNchambers = 4,  ///< Number of chambers
    kNplanes   = 8,  ///< Number of planes
    kNslats    = 18 ///< Number of slats
  };

  enum {kChHit, kChNonHit, kNcounts};

  enum {
    kHtracksInSlat  = 0,  ///< Tracks in slat histogram index
    kHtracksInBoard = 1,  ///< Tracks in board histogram index
    kHchamberEff    = 2,  ///< N44 per cathode histogram index
    kHchamberNonEff = 4,  ///< N33 per cathode histogram index
    kHslatEff       = 6,  ///< N44 per slat histogram index
    kHslatNonEff    = 14, ///< N33 per slat histogram index
    kHboardEff      = 22, ///< N44 per board histogram index
    kHboardNonEff   = 30, ///< N33 per board histogram index
    kHthetaX        = 38, ///< Angular distribution theta_x
    kHthetaY        = 39  ///< Angular distribution theta_y
  };
  
  /// Given cathode and chamber, return plane number
  Int_t GetPlane(Int_t cathode, Int_t chamber) { return kNchambers*cathode + chamber; }

  ClassDef(AliAnalysisTaskTrigChEff, 1); // Trigger chamber efficiency analysis
};

