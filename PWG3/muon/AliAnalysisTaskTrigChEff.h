#include "TH1.h"
#include "TList.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"

class AliAnalysisTaskTrigChEff : public AliAnalysisTask {
 public:
  AliAnalysisTaskTrigChEff(const char *name = "AliAnalysisTaskTrigChEff");
  virtual ~AliAnalysisTaskTrigChEff() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetType(const char* type) {fAnalysisType = type;}

protected:
  void ResetHistos();
  
  /// Getting flag telling which efficiency is performable for current track
  Int_t GetEffFlag(UShort_t pattern) { return (pattern >> 8) & 0x03; }
  
  /// Getting crossed slat
  Int_t GetSlat(UShort_t pattern) { return (pattern >> 10) & 0x1F; }
  
  Int_t IsChInefficient(UShort_t pattern, Int_t cathode);

private:
  /// Not implemented
  AliAnalysisTaskTrigChEff(const AliAnalysisTaskTrigChEff& rhs);
  /// Not implemented
  AliAnalysisTaskTrigChEff& operator = (const AliAnalysisTaskTrigChEff& rhs);
    
  AliESDEvent* fESD; //!< ESDevent object
  AliAODEvent* fAOD; //!< AODevent object
  TString fAnalysisType; //"ESD" or "AOD"

  TList*  fList; //TList output object

  enum {
    kNcathodes = 2,  ///< Number of cathodes
    kNchambers = 4,  ///< Number of chambers
    kNplanes   = 8,  ///< Number of planes
    kNslats    = 18 ///< Number of slats
  };

  enum {kAllChEff, kChNonEff, kNcounts};

  enum {
    kNoEff,
    kChEff,
    kSlatEff,
    kBoardEff
  };

  enum {
    kHtracksInSlat  = 0,  ///< Tracks in slat histogram index
    kHtracksInBoard = 1,  ///< Tracks in board histogram index
    kHchamberAllEff = 2,  ///< N44 per cathode histogram index
    kHchamberNonEff = 4,  ///< N33 per cathode histogram index
    kHslatAllEff    = 6,  ///< N44 per slat histogram index
    kHslatNonEff    = 14, ///< N33 per slat histogram index
    kHboardAllEff   = 22, ///< N44 per board histogram index
    kHboardNonEff   = 30  ///< N33 per board histogram index
  };
  
  /// Given cathode and chamber, return plane number
  Int_t GetPlane(Int_t cathode, Int_t chamber) { return kNchambers*cathode + chamber; }

  ClassDef(AliAnalysisTaskTrigChEff, 0); // Single muon analysis
};

