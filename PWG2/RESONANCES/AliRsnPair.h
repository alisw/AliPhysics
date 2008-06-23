#ifndef ALIRSNMVPAIR_H
#define ALIRSNMVPAIR_H

#include <TObject.h>
#include <TH1.h>
#include <TClonesArray.h>
#include <TTree.h>

#include "AliRsnDaughter.h"
#include "AliRsnPairParticle.h"
#include "AliRsnCutSet.h"
#include "AliRsnCutMgr.h"
#include "AliRsnPairDef.h"
#include "AliRsnPID.h"

#include "AliRsnEventBuffer.h"

class AliRsnEvent;
class AliRsnCut;


/**
  @author Martin Vala <Martin.Vala@cern.ch>
*/
class AliRsnPair : public TObject
{
  public:
    enum EPairType
    {
      kESDNoPID = 0,kESDNoPIDMix,
      kESDNormal, kESDMix,
      kMCNoPID,
      kMCNormal, kMCMix,
      kMCSignalOnly, kMCBackgroundOnly,
      kLastIndex
    };

    AliRsnPair();
    AliRsnPair ( AliRsnPair::EPairType type , AliRsnPairDef *pairDef,Int_t numOfMix = 0 );
    ~AliRsnPair();

    TString         GetESDParticleName ( AliRsnPID::EType type );
    TString         GetPairTypeName ( EPairType type );

    TH1F            *GenerateEffMassHist ( Int_t index = 0 );

    void            ProcessPair ( AliRsnEvent *event,TH1F *hist,Int_t index=0 );

    void            DoCleanUpAfterOneEvent();

    void            AddCutMgr ( AliRsnCutMgr* theValue );
    AliRsnCutMgr* GetCutMgr ( Int_t index ) { return ( AliRsnCutMgr* ) fCutMgrs.At ( index );}
    TObjArray*      GetCutMgr () { return &fCutMgrs;}
    Int_t           GetNumOfCutMgr () { return fCutMgrs.GetEntriesFast();}

    void SetMass ( Double_t theValue ) { fMass[0] = theValue; fMass[1] = theValue;}
    void SetMass ( Double_t theValue , Int_t index ) { fMass[index] = theValue; }
    Double_t GetMass ( Int_t index=0 ) const { return fMass[index]; }

    TString GetEffMassHistName ( Int_t index = 0 );
    TString GetEffMassHistTitle ( Int_t index = 0 );

    void SetNumOfMixEvent ( const Int_t& theValue ) { fNumOfMixEvent = theValue; }
    Int_t GetNumOfMixEvent() const { return fNumOfMixEvent;}

    void SetIsFilledOnlyInHistRange ( const Bool_t& theValue ) { fIsFilledOnlyInHistRange = theValue; }

    void                    PrepareMixForPair ( AliRsnEvent * event,TTree *tree );

    void SetRsnMVEventBuffer ( AliRsnEventBuffer* theValue ) { fRsnMVEventBuffer = theValue; }
    AliRsnEventBuffer* GetRsnMVEventBuffer() const { return fRsnMVEventBuffer; }


  private:

    AliRsnPairDef          fPairDef;                // pair definition
    AliRsnPair::EPairType  fPairType;               // pair type
  
    AliRsnCutMgr          *fCurrentCutMgr;          // cut manager
    TObjArray               fCutMgrs;               // array of cuts

    Double_t                fMass[2];               // mass for nopid

    AliRsnPairParticle      fEffMassParticle;

    Int_t                   fNumOfMixEvent;         // number of events to be mix with current one
    AliRsnEventBuffer      *fRsnMVEventBuffer;      // event buffer for event mixing

    Bool_t                  fIsSignSame;            // flag for same sign

    Bool_t                  fIsFilledOnlyInHistRange; // flag filling histogram

    void                    DoLoopPairESD ( AliRsnEvent *event1, TArrayI *array1,AliRsnEvent *event2, TArrayI *array2 ,TH1F *hist,Int_t index=0 );
    void                    DoLoopPairMC ( AliRsnEvent *event1, TArrayI*array1,AliRsnEvent *event2, TArrayI*array2 ,TH1F *hist,Int_t index=0 );

    void                    DoESDNoPID ( AliRsnEvent *event,TH1F *hist ,Int_t index=0 );
    void                    DoESDNoPIDMix ( AliRsnEvent *event,TH1F *hist,Int_t index=0 );
    void                    DoESDNormal ( AliRsnEvent *event,TH1F *hist,Int_t index=0 );
    void                    DoESDMix ( AliRsnEvent *event,TH1F *hist,Int_t index=0 );

    void                    DoMCNoPID ( AliRsnEvent *event,TH1F *hist,Int_t index=0 );
    void                    DoMCNormal ( AliRsnEvent *event,TH1F *hist,Int_t index=0 );

    ClassDef ( AliRsnPair, 1 );
};

#endif
