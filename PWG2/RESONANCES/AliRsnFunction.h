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

#include <TArrayD.h>
#include <TString.h>

#include "AliRsnCut.h"
#include "AliRsnHistoDef.h"
#include "AliRsnPairParticle.h"

class TH1D;
class TH2D;
class AliRsnEvent;

class AliRsnFunction : public TObject
{

public:

    enum EFcnType
    {
      kInvMass,
      kInvMassMC,
      kResolution,
      kPtSpectrum,
      kFcnTypes
    };

    AliRsnFunction();
    AliRsnFunction(EFcnType type, AliRsnHistoDef *hd, Bool_t skipOut = kTRUE);
    AliRsnFunction(const AliRsnFunction &copy);
    virtual ~AliRsnFunction() {Clear();}
    virtual void Clear(Option_t *option = "");

    Bool_t           UseBins() {return fUseBins;}
    Bool_t           SkipOut() {return fSkipOutsideInterval;}
    AliRsnHistoDef*  GetHistoDef() {return fHistoDef;}
    TString          GetFcnName();
    TString          GetFcnTitle();
    
    void  SetBinningCut(AliRsnCut::EType type, Double_t min, Double_t max, Double_t step);
    void  SetBinningCut(AliRsnCut::EType type, Int_t nbins, Double_t *bins); 
    void  SetHistoDef(AliRsnHistoDef *def) {fHistoDef = def;}

    // working routines
    TList* Init(const char *histoName, const char *histoTitle);
    void   Init(const char *histoName, const char *histoTitle, TList *tgt);
    Bool_t Fill(AliRsnPairParticle *pair, AliRsnPairDef *ref, Double_t weight = 0.0);
    Double_t FcnValue(AliRsnPairParticle *pair, AliRsnPairDef *ref);

private:
    
    const AliRsnFunction& operator=(const AliRsnFunction &copy);
    
    Double_t    FcnResolution(AliRsnPairParticle *pair, AliRsnPairDef *pd);
    
    EFcnType         fFcnType;       // function type

    Bool_t           fUseBins;       // flag to choose if binning is used
    Bool_t           fSkipOutsideInterval; // skip pairs which fall outside histogram interval
    
    TArrayD          fBins;          // low edge of each bin (upper is the low edge of next bin)
    AliRsnCut        fBinningCut;    // binning cut
    AliRsnCut::EType fBinningCutType;// binning cut type

    AliRsnHistoDef  *fHistoDef;      // definitions for histogram
    TH1D            *fHisto[100];    // binned histograms
    
    // ROOT dictionary
    ClassDef(AliRsnFunction, 1)
};

#endif
