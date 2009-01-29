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

#ifndef ALIRsnEventFunction_H
#define ALIRsnEventFunction_H

#include <TArrayD.h>
#include <TString.h>

#include "AliRsnCut.h"
#include "AliRsnCutSet.h"
#include "AliRsnHistoDef.h"
#include "AliRsnPairParticle.h"

class TH1D;
class TH2D;
class AliRsnEvent;

class AliRsnEventFunction : public TObject
{

  public:

    enum EType
    {
      kMultiplicity,
      kLeadingMomentum,
      kLeadingTheta,
      kAverageMomentum,
      kAngleLeadingMean,
      kAngleLeadingRMS,
      kPtResolution,
      kVtResolution,
      kVzResolution,

      kTypes
    };

    AliRsnEventFunction();
    AliRsnEventFunction(EType type, AliRsnHistoDef *hd,
                        AliRsnDaughter::EPIDMethod pidMethod = AliRsnDaughter::kNoPID,
                        AliRsnPID::EType pidType = AliRsnPID::kUnknown,
                        Char_t sign = '0');
    AliRsnEventFunction(const AliRsnEventFunction &copy);
    virtual ~AliRsnEventFunction() {Clear();}
    virtual void Clear(Option_t *option = "");

    Bool_t           UseBins() {return fUseBins;}
    AliRsnHistoDef*  GetHistoDef() {return fHistoDef;}
    TString          GetFcnName();

    void  SetBinningCut(AliRsnCut::EType type, Double_t min, Double_t max, Double_t step);
    void  SetBinningCut(AliRsnCut::EType type, Int_t nbins, Double_t *bins);
    void  SetHistoDef(AliRsnHistoDef *def) {fHistoDef = def;}
    void  SetEventCuts(AliRsnCutSet *cuts) {fEventCuts = cuts;}
    void  SetTrackCuts(AliRsnCutSet *cuts) {fTrackCuts = cuts;}
    void  SetLeadingPtMin(Double_t value) {fLeadPtMin = value;}

    // working routines
    void     Init(TList *tgt);
    Bool_t   Fill(AliRsnEvent *event);
    Double_t FcnValue(AliRsnEvent *event);

  private:

    const AliRsnEventFunction& operator=(const AliRsnEventFunction& /*copy*/) { return *this; }

    EType                      fType;      // function type

    AliRsnDaughter::EPIDMethod fPIDMethod; // PID method to be used
    AliRsnPID::EType           fPIDType;   // PID species to be used
    Char_t                     fCharge;    // charge sign
    Double_t                   fLeadPtMin; // smallest acceptable momentum of leading particle

    Bool_t                     fAccept;    // internal flag to check a computed value
    Bool_t                     fUseBins;   // flag to choose if binning is used

    TArrayD                    fBins;            // low edge of each bin (upper is the low edge of next bin)
    AliRsnCut                  fBinningCut;      // binning cut
    AliRsnCut::EType           fBinningCutType;  // binning cut type

    AliRsnCutSet              *fEventCuts;       // selection cuts for events
    AliRsnCutSet              *fTrackCuts;       // selection cuts for tracks in each event

    AliRsnHistoDef            *fHistoDef;        // definitions for histogram
    TH1D                      *fHisto[100];      // binned histograms

    // ROOT dictionary
    ClassDef(AliRsnEventFunction, 1)
};

#endif
