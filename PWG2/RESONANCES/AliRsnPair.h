//
// *** Class AliRsnPair ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef ALIRSNPAIR_H
#define ALIRSNPAIR_H

#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TArrayI.h"

#include "AliRsnDaughter.h"
#include "AliRsnPairDef.h"
#include "AliRsnEventBuffer.h"
#include "AliRsnEvent.h"
#include "AliRsnCutMgr.h"
#include "AliRsnHistoDef.h"
#include "AliRsnFunctionDef.h"
#include "AliRsnFunctionNew.h"

class AliRsnFunction;

class AliRsnPair : public TObject
{
  public:

    enum EPairType
    {
      kNoPID = 0,    kNoPIDMix,
      kRealisticPID, kRealisticPIDMix,
      kPerfectPID,   kPerfectPIDMix,
      kPairTypes
    };

    AliRsnPair(EPairType type = kRealisticPID, AliRsnPairDef *def = 0);
    ~AliRsnPair();

    void    Print(Option_t *option = "") const;
    void    ProcessPair(AliRsnEvent *ev1, AliRsnEvent *ev2 = 0);
    void    SetCutMgr(AliRsnCutMgr* theValue) { fCutMgr = theValue; }
    void    AddFunction(AliRsnFunction *fcn);
    //void    AddFunction(AliRsnFunctionDef *fcn);
    TList*  GenerateHistograms(TString prefix = "");
    void    GenerateHistograms(TString prefix, TList *tgt);

    Bool_t  IsMixed() {return fIsMixed;}
    Bool_t  IsPairEqual() {if (fPIDMethod == AliRsnDaughter::kNoPID) return (fPairDef->IsLikeSign());
                           else return (fPairDef->IsLikeSign() && fPairDef->HasEqualTypes());}

    TString GetPairTypeName(EPairType type) const;
    TString GetPairName() const;
    TString GetPairHistName(AliRsnFunction *fcn, TString text = "") const;
    TString GetPairHistTitle(AliRsnFunction *fcn, TString text = "") const;
    //TString GetPairHistName(AliRsnFunctionNew *fcn, TString text = "") const;
    //TString GetPairHistTitle(AliRsnFunctionNew *fcn, TString text = "") const;

  private:

    AliRsnPair (const AliRsnPair &copy) : TObject(copy),
      fIsMixed(kFALSE),fPairType(kPairTypes),fPIDMethod(AliRsnDaughter::kRealistic),
      fPairDef(0x0),fCutMgr(0x0),fFunctions("AliRsnFunction",0) {}
    AliRsnPair& operator=(const AliRsnPair&) {return *this;}

    void     SetUp(EPairType type);
    void     SetAllFlags(AliRsnDaughter::EPIDMethod pid, Bool_t mix) {fPIDMethod = pid; fIsMixed = mix;}
    void     LoopPair(AliRsnEvent *ev1, TArrayI *a1, AliRsnEvent *ev2, TArrayI *a2);

    Bool_t   CutPass(AliRsnDaughter *d);
    Bool_t   CutPass(AliRsnPairParticle *p);
    Bool_t   CutPass(AliRsnEvent *e);

    Bool_t                      fIsMixed;        // doing event-mixing ?
    EPairType                   fPairType;       // pair type (PID + mixing or not)
    AliRsnDaughter::EPIDMethod  fPIDMethod;      // pid type variable for single track

    AliRsnPairDef              *fPairDef;        // pair definition (particles, charges)
    AliRsnCutMgr               *fCutMgr;         // cut manager
    TClonesArray                fFunctions;      // functions

    ClassDef (AliRsnPair, 2)
};

#endif
