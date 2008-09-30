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

#include "AliRsnDaughter.h"
#include "AliRsnPairDef.h"
#include "AliRsnEventBuffer.h"
#include "AliRsnEvent.h"
#include "AliRsnCutMgr.h"
#include "AliRsnHistoDef.h"

class AliRsnFunction;

class AliRsnPair : public TObject
{
  public:

    enum EPairType
    {
      kNoPID = 0,    kNoPIDMix,
      kRealisticPID, kRealisticPIDMix,
      kPerfectPID,   kPerfectPIDMix,
      kTruePairs,
      kPairTypes
    };

    AliRsnPair(EPairType type = kRealisticPID, AliRsnPairDef *def = 0, 
               Int_t mixNum = 1, Double_t mixVzCut = 1.0, Int_t mixMultCut = 10);
    ~AliRsnPair();

    void    Init();
    void    Print();
    void    ProcessPair(AliRsnEventBuffer *buf);
    void    SetCutMgr(AliRsnCutMgr* theValue) { fCutMgr = theValue; }
    void    SetMixingCut(AliRsnCutSet* theValue) { fMixingCut = theValue; }
    void    AddFunction(AliRsnFunction *fcn);
    TList*  GenerateHistograms(TString prefix = "");
    void    GenerateHistograms(TString prefix, TList *tgt);

    TString GetPairTypeName(EPairType type);
    TString GetPairName();
    TString GetPairHistName(AliRsnFunction *fcn, TString text = "");
    TString GetPairHistTitle(AliRsnFunction *fcn, TString text="");

  private:

    AliRsnPair (const AliRsnPair &copy) : TObject(copy),
      fIsMixed(kFALSE),fUseMC(kFALSE),fIsLikeSign(kFALSE),fMixNum(1),fMixingCut(0x0),
      fPairDef(0x0),fPairType(kPairTypes),fTypePID(AliRsnDaughter::kRealistic),
      fCutMgr(0x0),fFunctions("AliRsnFunction",0) {}
    AliRsnPair& operator=(const AliRsnPair&) {return *this;}

    void           SetUp(EPairType type);  // sets up all flags
    void           SetAllFlags(AliRsnDaughter::EPIDMethod pidType,Bool_t isMix, Bool_t useMC);
    AliRsnEvent*   FindEventByEventCut(AliRsnEventBuffer *buf,Int_t & num);
    void           LoopPair(AliRsnEvent *ev1,TArrayI *a1,AliRsnEvent *ev2,TArrayI *a2);

    Bool_t         CutPass(AliRsnDaughter *d);
    Bool_t         CutPass(AliRsnPairParticle *p);
    Bool_t         CutPass(AliRsnEvent *e);
    
    // flags & integer data
    Bool_t         fIsMixed;                 // doing event-mixing ?
    Bool_t         fUseMC;                   // using MC inv. mass ?
    Bool_t         fIsLikeSign;              // is a like-sign pair ?
    Int_t          fMixNum;                  // number of mixed events
    AliRsnCutSet  *fMixingCut;               // cut for event mixing
    
    // work management
    AliRsnPairDef              *fPairDef;                // pair definition (particles, charges)
    EPairType                   fPairType;               // pair type (PID + mixing or not)
    AliRsnDaughter::EPIDMethod  fTypePID;                // pid type variable for single track
    AliRsnCutMgr               *fCutMgr;                 // cut manager
    TClonesArray                fFunctions;              // functions

    ClassDef (AliRsnPair, 1)
};

#endif
