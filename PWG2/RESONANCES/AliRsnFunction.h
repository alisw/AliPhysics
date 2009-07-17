//
// Class AliRsn Fcn
//
// This class defines a base classe to implement a typical computation
// which uses the internal RSN package event format (AliRsnEvent).
// It contains some default flags which turn out to be useful:
//  - a flag to select only the "true" pairs (tracks from same resonance)
//  - a flag to know if the computation is done over two events (mixing)
//
// Any kind of analysis object should be implemented as inheriting from this
// because the AliRsnAnalyzer which executes the analysis will accept a collection
// of such objects, in order to have a unique format of processing method
//
// The user who implements a kind of computation type should inherit from
// this class and override the virtual functions defined in it, which
// initialize the final output histogram and define how to process data.
//
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNFUNCTION_H
#define ALIRSNFUNCTION_H

#include <TClonesArray.h>
#include <THnSparse.h>
#include <TNamed.h>

#include "AliRsnCut.h"

class AliRsnFunctionAxis;
class AliRsnPairParticle;
class AliRsnPairDef;

class AliRsnFunction : public TNamed
{

  public:

    AliRsnFunction();
    AliRsnFunction(const AliRsnFunction &copy);
    virtual ~AliRsnFunction() { delete fHistogram; }
    const AliRsnFunction& operator=(const AliRsnFunction &copy);

    void                 SetPairDef(AliRsnPairDef * const def) {fPairDef = def;}
    void                 SetTrack(AliRsnDaughter * const track) {fTrack = track;}
    void                 SetPair(AliRsnPairParticle * const pair) {fPair = pair;}
    void                 SetEvent(AliRsnEvent *const event) {fEvent = event;}

    AliRsnPairDef*       GetPairDef() const {return fPairDef;}
    AliRsnDaughter*      GetTrack() const {return fTrack;}
    AliRsnPairParticle*  GetPair() const {return fPair;}
    AliRsnEvent*         GetEvent() const {return fEvent;}
    virtual const char*  GetName() const;

    void                 AddAxis(AliRsnFunctionAxis*const axis);
    Int_t                GetNumberOfAxes() {return fAxisList.GetEntries();}
    THnSparseD*          CreateHistogram(const char *histoName, const char *histoTitle);

    Bool_t               Fill();

  protected:

    AliRsnPairDef      *fPairDef;     // reference to used pair definition
    TClonesArray        fAxisList;    // list of axis

    AliRsnDaughter     *fTrack;       // processed track
    AliRsnPairParticle *fPair;        // processed pair
    AliRsnEvent        *fEvent;       // processed event

    THnSparseD         *fHistogram;   // output histogram

    // ROOT dictionary
    ClassDef(AliRsnFunction, 3)
};

#endif
