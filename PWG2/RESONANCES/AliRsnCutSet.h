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

#include "AliRsnTarget.h"

class AliRsnCut;
class AliRsnDaughter;
class AliRsnExpression;
class AliRsnPairParticle;
class AliRsnEvent;

class AliRsnCutSet : public AliRsnTarget {
public:

   AliRsnCutSet();
   AliRsnCutSet(const char *name, RSNTARGET target);
   AliRsnCutSet(const AliRsnCutSet &copy);
   AliRsnCutSet& operator=(const AliRsnCutSet& copy);
   ~AliRsnCutSet();

   void      AddCut(AliRsnCut* cut);

   void      ShowCuts() const;
   Int_t     GetIndexByCutName(TString s);
   Bool_t    Passed();
   Bool_t    IsValidScheme();
   TString   ShowCutScheme() const;
   Int_t     TestExpression(TString opt = "short");
   void      PrintSetInfo();

   Bool_t    IsSelected(TObject *object);

   void SetBoolValue(Bool_t theValue, Int_t index) { fBoolValues[index] = theValue; }
   Bool_t GetBoolValue(Int_t index) const { return fBoolValues[index]; }

   void SetCutScheme(const char *theValue);
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

   ClassDef(AliRsnCutSet, 1)   // ROOT dictionary
};

#endif
