//
// Class AliRsnCutSet
//
// This is the front-end for cut management and checking.
// It must be prepared by adding all required single cuts,
// and then with a logical expression which combines all cuts
// with the "AND", "OR" and "NOT" operators.
//
// author: M. Vala (martin.vala@cern.ch)
//

#ifndef ALIRSNCUTSET_H
#define ALIRSNCUTSET_H

#include <TNamed.h>
#include <TObjArray.h>

#include "AliRsnCut.h"

class AliRsnDaughter;
class AliRsnExpression;
class AliRsnPairParticle;
class AliRsnEvent;

class AliRsnCutSet : public TNamed
{
  public:

    AliRsnCutSet();
    AliRsnCutSet(TString name);
    AliRsnCutSet(const AliRsnCutSet &copy);
    AliRsnCutSet& operator=(const AliRsnCutSet& /*copy*/) {return (*this);}
    ~AliRsnCutSet();

    void      AddCut(AliRsnCut* cut);

    void      ShowCuts() const;
    Int_t     GetIndexByCutName(TString s);
    Bool_t    Passed();
    Bool_t    IsValidScheme();
    TString   ShowCutScheme();
    Int_t     TestExpression(TString opt="short");
    void      PrintSetInfo();

    Bool_t    IsSelected(AliRsnCut::ETarget type, AliRsnDaughter *daughter);
    Bool_t    IsSelected(AliRsnCut::ETarget type, AliRsnPairParticle *pair);
    Bool_t    IsSelected(AliRsnCut::ETarget type, AliRsnEvent *event);
    Bool_t    IsSelected(AliRsnCut::ETarget type, AliRsnEvent *ev1, AliRsnEvent *ev2);

    void SetBoolValue(Bool_t theValue,Int_t index) { fBoolValues[index] = theValue; }
    Bool_t GetBoolValue(Int_t index) const { return fBoolValues[index]; }

    void SetCutScheme(const TString& theValue);
    TString GetCutScheme() const { return fCutScheme; }

    void SetCutSchemeIndexed(TString theValue);
    TString   GetCutSchemeIndexed();

    TObjArray *GetCuts() { return &fCuts; }

  private:

    TObjArray     fCuts;                  // array of cuts
    Int_t         fNumOfCuts;             // number of cuts
    TString       fCutScheme;             // cut scheme
    TString       fCutSchemeIndexed;      // cut scheme indexed

    Bool_t       *fBoolValues;            //[fNumOfCuts]
    Bool_t        fIsScheme;              // is scheme

    AliRsnExpression  *fExpression;       // pointer to AliRsnExpression

    ClassDef(AliRsnCutSet,1)    // ROOT dictionary
};

#endif
