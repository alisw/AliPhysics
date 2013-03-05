#ifndef ALIMUONEVENTCUTS_H
#define ALIMUONEVENTCUTS_H

#include "AliAnalysisCuts.h"

class AliVEvent;
class AliVVertex;
class TList;
class TAxis;
class TArrayI;
class TString;
class TObjString;
class TObjArray;

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
  
  // Handle trigger
  void SetTrigClassPatterns ( const TString trigPattern );
  /// Get default trigger class patterns
  TString GetDefaultTrigClassPatterns() { return fDefaultTrigClassPatterns; };
  void SetTrigClassLevels ( const TString pattern = "MSL:Lpt,MUSL:Lpt,MSH:Hpt,MUSH:Hpt,MUL:LptLpt,MUU:LptLpt,MLL:LptLpt" );
  TArrayI GetTrigClassPtCutLevel ( const TString trigClassName ) const;
  void SetTrigInputsMap ( const TString trigInputsMap );
  /// Get trigger classes found in run
  TList* GetAllSelectedTrigClasses () const { return fAllSelectedTrigClasses; }
  const TObjArray* GetSelectedTrigClassesInEvent ( const AliVEvent* event );

  
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
  void SetVertexMinNContributors ( const Int_t vertexMinNContributors ) { fVertexMinNContributors = vertexMinNContributors; }
  /// Get minimum number of vertex contributors
  Int_t GetVertexMinNContributors () const { return fVertexMinNContributors; }
  /// Set Vz limits
  void SetVertexVzLimits ( Double_t vzMin = -999., Double_t vzMax = 999. ) { fVertexVzMin = vzMin; fVertexVzMax = vzMax; }
  /// Get Vtx vz min
  Double_t GetVertexVzMin () const { return fVertexVzMin; }
  /// Get Vtx vz max
  Double_t GetVertexVzMax () const { return fVertexVzMax; }

  void Print ( Option_t* option = "" ) const;

 protected:
  
  void BuildTriggerClasses ( const TString firedTrigClasses, UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs );
  Bool_t CheckTriggerClassPattern ( const TString& toCheck ) const;
  Bool_t CheckTriggerClassCombination ( const TObjArray* combo, const TString& firedTriggerClasses, UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs ) const;
  void AddToEventSelectedClass ( const TString& toCheck, const TObjString* foundTrig );
  Bool_t UpdateEvent( const AliVEvent* event );
  void SetDefaultTrigClassPatterns();
  void SetDefaultTrigInputsMap();
    
  UInt_t fPhysicsSelectionMask; ///< Physics selection mask
  
  Int_t fVertexMinNContributors;  ///< Minimum number of SPD vertex contributors
  Double_t fVertexVzMin;          ///< SPD vertex Vz min
  Double_t fVertexVzMax;          ///< SPD vertex Vz max
  
  TString fDefaultTrigClassPatterns; ///< Default trigger class patterns
  TObjArray* fSelectedTrigPattern; ///< List of triggers to be kept
  TObjArray* fRejectedTrigPattern; ///< List of triggers to be rejected
  TObjArray* fSelectedTrigLevel;   ///< Track-trigger pt cut for selected trigger class
  TObjArray* fSelectedTrigCombination; ///< Selected trigger combinations
  TList* fTrigInputsMap;       ///< Trigger inputs map
  TList* fAllSelectedTrigClasses;  ///< List of all selected trigger classes found
  TAxis* fCentralityClasses;   ///< Centrality classes
  
  private:
  ULong64_t fEventTriggerMask; //!< Fired trigger mask in the event
  TObjArray* fSelectedTrigClassesInEvent; //!< list of selected trigger classes in current event 
  
  ClassDef(AliMuonEventCuts, 4); // Class for muon event filters
};

#endif

