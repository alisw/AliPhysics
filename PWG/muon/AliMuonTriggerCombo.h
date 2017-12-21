/// \class AliMuonTriggerCombo
/// \brief Combination of triggers classes, trigger inputs or physics selection bits
///
/// The class allows to define a combination of trigger classes,
/// trigger inputs or physics selection bits
/// It is used by AliMuonEventCuts
///
/// \author Diego Stocco <dstocco@cern.ch>, Subatech
/// \date Sep 23, 2015

#ifndef ALIMUONTRIGGERCOMBO_H
#define ALIMUONTRIGGERCOMBO_H

#include "TNamed.h"

class TObjArray;
class THashList;
class TObjString;

class AliMuonTriggerCombo : public TNamed
{
public:
  AliMuonTriggerCombo ();
  AliMuonTriggerCombo ( const char* name, const char* trigInputsString, const char* trigPtMatchLevel );
  virtual ~AliMuonTriggerCombo();

  Bool_t MatchEvent ( const TString& firedTriggerClasses, UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs, UInt_t physicsSelection );

  /// Get trigger matching level
  Int_t GetTrigMatchLevel () const { return fTrigPtMatchLevel; }

  /// Is dimuon trigger
  Bool_t IsDimuTrigger () const { return fIsDimuTrig; }

  /// Get combination type
  Int_t GetType() const { return fComboType; }

  enum {
    kBadPattern,    //!<! Ill formed pattern
    kComboSimple,   //!<! Single element
    kComboFormula,  //!<! Complex formula
    kComboAND,      //!<! Formula with only AND
    kComboOR,       //!<! Formula with onlt OR
    kRegex,         //!<! Match pattern
  };

  /// Combination has trigger inputs
  Bool_t HasTriggerClasses() { return fTriggerClasses ? kTRUE : kFALSE; }


private:
  Bool_t Init ( const char* name, const char* trigInputsString, const char* trigPtMatchLevel );
  TObjArray* GetInputs ( Int_t ilevel );
  THashList* GetTrigInputsMap ( const char* trigInputsString ) const;
  THashList* GetPhysSelBits () const;
  Bool_t CheckElement ( TString& formula, const TObjString* element, Bool_t ok ) const;

  /// not implemented
  AliMuonTriggerCombo& operator=(const AliMuonTriggerCombo &rhs);
  /// not implemented
  AliMuonTriggerCombo(const AliMuonTriggerCombo &iter);

  TObjArray* fL0Inputs; ///< List of L0 inputs required
  TObjArray* fL1Inputs; ///< List of L1 inputs required
  TObjArray* fL2Inputs; ///< List of L2 inputs required
  TObjArray* fPhysSelBits; ///< List of physics selection bits required
  TObjArray* fTriggerClasses; ///< List of trigger classes regexp
  UInt_t fComboType; ///< Type of combination
  Int_t fTrigPtMatchLevel; ///< Trigger pt cut level for this combination
  Bool_t fIsDimuTrig; ///< Is di-muon trigger

  /// \cond CLASSIMP
  ClassDef(AliMuonTriggerCombo, 2); // Class for muon event trigger combination
  /// \endcond
};

#endif
