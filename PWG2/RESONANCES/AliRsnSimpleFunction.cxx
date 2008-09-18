//
// Class AliRsnSimpleFcn
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

#include <Riostream.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include "AliLog.h"

#include "AliRsnMCInfo.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnCut.h"

#include "AliRsnSimpleFunction.h"

ClassImp(AliRsnSimpleFunction)

//________________________________________________________________________________________
AliRsnSimpleFunction::AliRsnSimpleFunction
(const char *name, AliRsnPairDef *pd, AliRsnHistoDef *hd,
 AliRsnCutMgr *cuts, Option_t *option) :
  TNamed(name, ""),
  fTrueFlag(kFALSE),
  fMixFlag(kFALSE),
  fPair(),
  fPairDef(pd),
  fHistoDef(hd),
  fCuts(cuts),
  fHisto1D(0x0),
  fHisto2D(0x0)
{
//
// Constructor.
// The histogram data member cannot be passed externally,
// its initialization MUST be defined inside the Init() method,
// which must be overridden in any derivate implementation.
// Last argument (option) allows to set in a user-friedly way the flags:
// - TRUE (for true pairs)
// - MIX (for mixing pairs)
//

    TString opt(option);
    opt.ToUpper();
    fTrueFlag = opt.Contains("TRUE");
    fMixFlag = opt.Contains("MIX");
}
//________________________________________________________________________________________
AliRsnSimpleFunction::AliRsnSimpleFunction(const AliRsnSimpleFunction &copy) :
  TNamed(copy),
  fTrueFlag(copy.fTrueFlag),
  fMixFlag(copy.fMixFlag),
  fPair(),
  fPairDef(copy.fPairDef),
  fHistoDef(copy.fHistoDef),
  fCuts(copy.fCuts),
  fHisto1D(0x0),
  fHisto2D(0x0)
{
//
// Copy constructor.
// Pointer data members are not cloned.
// The histogram is not copied, because it is always initialized
// by the Init() method.
//
}
//________________________________________________________________________________________
const AliRsnSimpleFunction& AliRsnSimpleFunction::operator=
(const AliRsnSimpleFunction &copy)
{
//
// Assignment operator.
// Behaves like copy constructor.
// Also in this case, the histogram is not copied, and,
// if it was present, it is destroyed and will need to be recreated.
//

    Clear();
    SetName(copy.GetName());

    fTrueFlag = copy.fTrueFlag;
    fMixFlag = copy.fMixFlag;
    
    fPairDef = copy.fPairDef;
    fHistoDef = copy.fHistoDef;
    fCuts = copy.fCuts;

    return (*this);
}
//________________________________________________________________________________________
void AliRsnSimpleFunction::Clear(Option_t* /*option*/)
{
//
// Clear arrays and histogram.
// For the sake of security, all pointers are also set explicitly to NULL.
//
    if (fHisto1D) delete fHisto1D;
    fHisto1D = 0x0;
    if (fHisto2D) delete fHisto2D;
    fHisto2D = 0x0;
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFunction::Init()
{
//
// Initialization function.
// By default, it initializes the owned histogram using the method
// from AliRsnHistoDef class, giving the same name and title of this.
// A user can override this behaviour, if necessary.
// Before creating, the HistoDef is checked for proper initialization.
//

    Clear();

    Char_t name[200];
    sprintf(name, "h_%s", GetName());

    fHisto1D = (TH1D*)fHistoDef->CreateHistogram(name, "");
    return kTRUE;
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFunction::ProcessOne(AliRsnEvent* /*event*/)
{
//
// Process a single event according to the purposes of the derived class.
// Here the user must put the code which defines how the internal histogram is filled.
// The boolean return value allows to catch a failure in execution.
//

    AliWarning("This method must be defined in the derived class");
    return kFALSE;
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFunction::ProcessTwo(AliRsnEvent* /*event1*/, AliRsnEvent* /*event2*/)
{
//
// Process two events according to relations between them.
// This method must be defined whenever the derived object is used in event mixing.
// The boolean return value allows to catch a failure in execution.
//

    AliWarning("This method must be defined in the derived class");
    return kFALSE;
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFunction::CutPass(AliRsnDaughter *d)
{
//
// Check if the AliRsnDaughter argument pass its cuts.
// If the cut data member is not initialized for it, returns kTRUE.
//
    if (!fCuts) return kTRUE;
    else return fCuts->IsSelected(AliRsnCut::kParticle, d);
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFunction::CutPass(AliRsnPairParticle *p)
{
//
// Check if the AliRsnPairParticle argument pass its cuts.
// If the cut data member is not initialized for it, returns kTRUE.
// In this case, another separate check which could be necessary
// concerns the possibility that the two tracks are a "true pair" of 
// daughters of the same resonance. If the corresponding flag is set,
// this further check is done, and the method returns kTRUE only 
// when also this check is passed.
//

    Bool_t cutPassed, truePair;
    
    cutPassed = fCuts->IsSelected(AliRsnCut::kPair, p);
    truePair = fPair.IsTruePair(fPairDef->GetMotherPDG());
    
    if (fTrueFlag) {
        return (cutPassed && truePair);
    }
    else {
        return cutPassed;
    }
}

//________________________________________________________________________________________
Bool_t AliRsnSimpleFunction::CutPass(AliRsnEvent *e)
{
//
// Check if the AliRsnEvent argument pass its cuts.
// If the cut data member is not initialized for it, returns kTRUE.
//
    if (!fCuts) return kTRUE;
    else return fCuts->IsSelected(AliRsnCut::kEvent, e);
}
