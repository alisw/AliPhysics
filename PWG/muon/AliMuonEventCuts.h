#ifndef ALIMUONEVENTCUTS_H
#define ALIMUONEVENTCUTS_H

#include "AliAnalysisCuts.h"
#include "TArrayD.h"

class AliVEvent;
class AliVVertex;
class TList;
class TAxis;
class AliTimeStamp;
class TArrayI;

class AliMuonEventCuts : public AliAnalysisCuts
{
 public:
  
  enum {
    kPhysicsSelected = BIT(0),
    kSelectedCentrality = BIT(1),
    kSelectedTrig = BIT(2),
    kGoodVertex = BIT(3)
  };
  
  AliMuonEventCuts();
  AliMuonEventCuts(const char* name, const char* title);
  AliMuonEventCuts(const AliMuonEventCuts& obj);
  AliMuonEventCuts& operator=(const AliMuonEventCuts& obj);
  
  virtual ~AliMuonEventCuts();
  
  virtual UInt_t GetSelectionMask ( const TObject* obj );
  virtual Bool_t IsSelected ( TObject* obj );
  virtual Bool_t IsSelected ( TList* /*list */ );
  
  void SetDefaultFilterMask();
  void SetDefaultParameters();
  
  enum {
    kVertexMinNContributors,
    kVertexVzMin,
    kVertexVzMax,
    kNParameters
  };
  
  // Handle trigger
  void SetTrigClassPatterns ( const TString pattern );
  /// Get default trigger class patterns
  TString GetDefaultTrigClassPatterns() { return fDefaultTrigClassPatterns; };
  void SetTrigClassLevels ( const TString pattern = "MSL:Lpt MSH:Hpt MUL:LptLpt MLL:LptLpt" );
  TArrayI GetTrigClassPtCutLevel ( const TString trigClassName ) const;
  /// Get trigger classes found in run
  TList* GetAllSelectedTrigClasses () const { return fAllSelectedTrigClasses; }
  TObjArray* GetSelectedTrigClassesInEvent ( const AliVEvent* event );

  
  // Handle centrality
  void SetCentralityClasses(Int_t nCentralityBins = -1, Double_t* centralityBins = 0x0);
  /// Get centrality classes
  TAxis* GetCentralityClasses() const { return fCentralityClasses; }
  
  void SetCentralityEstimator ( const TString centralityEstimator = "V0M" );
  TString GetCentralityEstimator () const;
  Double_t GetCentrality ( const AliVEvent* event ) const;
  
  
  /// Set Physics selection mask
  void SetPhysicsSelectionMask ( const UInt_t physicsSelectionMask ) { fPhysicsSelectionMask = physicsSelectionMask; }
  
  
  /// Set minimum number of vertex contributors
  void SetVertexMinNContributors ( const Int_t vertexMinNContributors ) { SetParameter(kVertexMinNContributors, (Double_t)vertexMinNContributors); }
  /// Get minimum number of vertex contributors
  Int_t GetVertexMinNContributors () const { return (Int_t)fParameters[kVertexMinNContributors]; }
  /// Set Vz limits
  void SetVertexVzLimits ( Double_t vzMin = -999., Double_t vzMax = 999. ) { SetParameter(kVertexVzMin, vzMin); SetParameter(kVertexVzMax, vzMax); }
  /// Get Vtx vz min
  Double_t GetVertexVzMin () const { return fParameters[kVertexVzMin]; }
  /// Get Vtx vz max
  Double_t GetVertexVzMax () const { return fParameters[kVertexVzMax]; }
  
  //Bool_t SetRun(Int_t runNumber);
  //void SetUseCustomParam( Bool_t useCustomParam = kTRUE, Int_t runNumber = -1 );
  //void SetIsMC(Bool_t isMC = kTRUE) { fIsMC = isMC; }

  void Print ( Option_t* option = "" ) const;

  //Bool_t StreamParameters ( Int_t runNumber, Int_t maxRun );

 protected:
  
  Bool_t SetParameter ( Int_t iparam, Float_t value );
  void BuildTriggerClasses ( const TString firedTrigClasses );
  Bool_t UpdateEvent( const AliVEvent* event );
  void SetDefaultTrigClassPatterns();
  
  //Bool_t RunMatchesRange ( Int_t runNumber, const Char_t* objName ) const;

  //Bool_t fIsMC;             ///< Monte Carlo analysis
  //Bool_t fUseCustomParam;   ///< Use custom parameters (do not search in OADB)
  
  UInt_t fPhysicsSelectionMask; ///< Physics selection mask
  
  TArrayD fParameters;      ///< List of parameters
  
  TString fDefaultTrigClassPatterns; ///< Default trigger class patterns
  TObjArray* fSelectedTrigPattern; ///< List of triggers to be kept
  TObjArray* fRejectedTrigPattern; ///< List of triggers to be rejected
  TObjArray* fSelectedTrigLevel;   ///< Track-trigger pt cut for selected trigger class
  TList* fAllSelectedTrigClasses;  ///< List of all selected trigger classes found
  TAxis* fCentralityClasses;   ///< Centrality classes
  
  private:
  AliTimeStamp* fTimeStamp; //!< current event time stamp
  TObjArray* fSelectedTrigClassesInEvent; //!< list of selected trigger classes in current event 
  
  ClassDef(AliMuonEventCuts, 1); // Class for muon event filters
};

#endif

