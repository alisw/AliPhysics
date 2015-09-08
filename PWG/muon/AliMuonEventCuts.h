#ifndef ALIMUONEVENTCUTS_H
#define ALIMUONEVENTCUTS_H

#include "AliAnalysisCuts.h"

class AliVEvent;
class AliVVertex;
class AliInputEventHandler;
class THashList;
class TList;
class TAxis;
class TArrayI;
class TString;
class TObjString;
class TObjArray;
class AliAnalysisUtils;

class AliMuonEventCuts : public AliAnalysisCuts
{
 public:
  
  enum {
    kPhysicsSelected = BIT(0),
    kSelectedCentrality = BIT(1),
    kSelectedTrig = BIT(2),
    kGoodVertex = BIT(3),
    kNoPileup = BIT(4)
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
  
  /// Skip tests which are not active in the filter mask
  void SkipTestsNonInFilterMask ( UInt_t skipMask = 0xFFFF) { fCheckMask = ~skipMask; }
  
  // Handle trigger
  void SetTrigClassPatterns ( TString trigPattern, TString trigInputsMap = "" );
  /// Get default trigger class patterns
  TString GetDefaultTrigClassPatterns() const;
  TString GetDefaultTrigInputsMap() const;
  void SetTrigClassLevels (TString pattern = "MSL:Lpt,MUSL:Lpt,MSH:Hpt,MUSH:Hpt,MUL:LptLpt,MUU:LptLpt,MLL:LptLpt" );
  TArrayI GetTrigClassPtCutLevel (TString trigClassName ) const;
  /// Get trigger classes found in run
  THashList* GetAllSelectedTrigClasses () const { return fAllSelectedTrigClasses; }
  const TObjArray* GetSelectedTrigClassesInEvent ( const TString& firedTriggerClasses,
                                                  UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs, UInt_t physicsSelection );
  const TObjArray* GetSelectedTrigClassesInEvent ( const TString& firedTriggerClasses,
                                                  UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs );
  const TObjArray* GetSelectedTrigClassesInEvent ( const AliVEvent* event );
  const TObjArray* GetSelectedTrigClassesInEvent ( const AliInputEventHandler* eventHandler );

  UInt_t GetTriggerInputBitMaskFromInputName(const char* inputName) const;

  // Handle centrality
  void SetCentralityClasses(Int_t nCentralityBins = -1, Double_t* centralityBins = 0x0);
  /// Get centrality classes
  TAxis* GetCentralityClasses() const { return fCentralityClasses; }
  
  void SetCentralityEstimator (TString centralityEstimator = "V0M" );
  TString GetCentralityEstimator () const;
  Double_t GetCentrality ( const AliVEvent* event ) const;
  
  
  /// Set Physics selection mask
  void SetPhysicsSelectionMask (UInt_t physicsSelectionMask ) { fPhysicsSelectionMask = physicsSelectionMask; }
  
  void SetPhysSelBits();
  
  /// Set minimum number of vertex contributors
  void SetVertexMinNContributors (Int_t vertexMinNContributors ) { fVertexMinNContributors = vertexMinNContributors; }
  /// Get minimum number of vertex contributors
  Int_t GetVertexMinNContributors () const { return fVertexMinNContributors; }
  /// Set Vz limits
  void SetVertexVzLimits ( Double_t vzMin = -999., Double_t vzMax = 999. ) { fVertexVzMin = vzMin; fVertexVzMax = vzMax; }
  /// Get Vtx vz min
  Double_t GetVertexVzMin () const { return fVertexVzMin; }
  /// Get Vtx vz max
  Double_t GetVertexVzMax () const { return fVertexVzMax; }
  
  /// Return pointer to analysis utils (to configure cuts)
  AliAnalysisUtils* GetAnalysisUtils ( ) { return fAnalysisUtils; }

  void Print ( Option_t* option = "" ) const;

 protected:
  
  enum {
    kL0Input,     /// L0input index
    kL1Input,     /// L1input index
    kL2Input,     /// L2input index
    kPhysSelBit,  /// Physcs selection bit index
    kTrigClass,   /// Trigger class index
    kNtypes
  };

  void BuildTriggerClasses ( TString firedTrigClasses, UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs, UInt_t physicsSelection );
  Bool_t CheckTriggerClassPattern ( const TString& toCheck ) const;
  Bool_t CheckTriggerClassCombination ( const TObjArray* combo, const TString& firedTriggerClasses, UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs, UInt_t physicsSelection ) const;
  void AddToEventSelectedClass ( const TString& toCheck, const TObjString* foundTrig, const UInt_t comboType = 0 );
  Bool_t UpdateEvent ( const AliVEvent* event, UInt_t physicsSelection );
  void SetDefaultTrigClassPatterns();
  void SetTrigInputsMap ( TString trigInputsMap );
    
  UInt_t fPhysicsSelectionMask; ///< Physics selection mask
  
  Int_t fVertexMinNContributors;  ///< Minimum number of SPD vertex contributors
  Double_t fVertexVzMin;          ///< SPD vertex Vz min
  Double_t fVertexVzMax;          ///< SPD vertex Vz max
  
  UInt_t fCheckMask;              ///< Mask telling which cuts to check (by default check filter mask)
  
  TObjArray* fSelectedTrigPattern; ///< List of triggers to be kept
  TObjArray* fRejectedTrigPattern; ///< List of triggers to be rejected
  TObjArray* fSelectedTrigLevel;   ///< Track-trigger pt cut for selected trigger class
  TObjArray* fSelectedTrigCombination; ///< Selected trigger combinations
  THashList* fTrigInputsMap;       ///< Trigger inputs map
  THashList* fPhysSelBits;       ///< Physics selection bits
  THashList* fAllSelectedTrigClasses;  ///< List of all selected trigger classes found
  TAxis* fCentralityClasses;   ///< Centrality classes
  AliAnalysisUtils* fAnalysisUtils;    ///< Analysis utility
  
  private:
  ULong64_t fEventTriggerMask; //!< Fired trigger mask in the event
  UInt_t fEventL0Inputs; //!< L0 trigger inputs in the event
  UInt_t fEventL1Inputs; //!< L1 trigger inputs in the event
  UInt_t fEventL2Inputs; //!< L2 trigger inputs in the event
  UInt_t fEventPS; //!< Physics selection for the event
  TObjArray* fSelectedTrigClassesInEvent; //!< list of selected trigger classes in current event
  enum {kComboSimple, kComboFormula, kComboAND, kComboOR}; //!< Trigger combination types
  
  ClassDef(AliMuonEventCuts, 9); // Class for muon event filters
};

#endif

