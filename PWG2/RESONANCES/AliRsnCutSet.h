#ifndef ALIRSNCUTSET_H
#define ALIRSNCUTSET_H

#include <TNamed.h>
#include <TObjArray.h>

#include "AliRsnCut.h"
// class AliRsnCut;
class AliRsnCutMgr;

class AliRsnDaughter;
class AliRsnExpression;
class AliRsnPairParticle;

class AliRsnCutSet : public TNamed
{

  public:
    AliRsnCutSet();
    AliRsnCutSet ( TString name );
    AliRsnCutSet ( const AliRsnCutSet &copy );

    ~AliRsnCutSet();

    void      AddCut ( AliRsnCut* cut );
//     void      SetCutScheme ( TString sheme ) {fCutScheme = sheme;}

    void      ShowCuts();
    Int_t     GetIndexByCutName ( TString s );
    Bool_t    Passed();
    Bool_t    IsValidScheme();
    TString   ShowCutScheme();
    Int_t     TestExpression ( TString opt="short" );
    void      PrintSetInfo();

    Bool_t    IsSelected (AliRsnCut::ECutSetType type, AliRsnDaughter *daughter );
    Bool_t    IsSelected (AliRsnCut::ECutSetType type, AliRsnPairParticle *pair );
//     Bool_t    IsSelected ( AliRsnEffectiveParticle *pair );
//     Bool_t    IsSelected(TObject*obj);

    void SetBoolValue ( Bool_t theValue,Int_t index ) { fBoolValues[index] = theValue; }
    Bool_t GetBoolValue ( Int_t index ) const { return fBoolValues[index]; }

    void SetCutScheme ( const TString& theValue );
    TString GetCutScheme() const { return fCutScheme; }

    void SetCutSchemeIndexed ( TString theValue );
    TString   GetCutSchemeIndexed();

    TObjArray *GetCuts() { return &fCuts; }


  private:

    AliRsnCutSet& operator=(const AliRsnCutSet& /*copy*/) {return (*this);}

    TObjArray     fCuts;                  // array of cuts
    Int_t         fNumOfCuts;             // number of cuts
    TString       fCutScheme;             // cut scheme
    TString       fCutSchemeIndexed;      // cut scheme indexed

    Bool_t        *fBoolValues;           //[fNumOfCuts]
    Bool_t        fIsScheme;              // is scheme
    AliRsnExpression  *fExpression;     // pointer to AliRsnExpression

    ClassDef ( AliRsnCutSet,1 );
};

#endif
