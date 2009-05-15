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

#ifndef ALIRSNFunction_H
#define ALIRSNFunction_H

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNamed.h>

#include "AliRsnCut.h"
#include "AliRsnHistoDef.h"
#include "AliRsnPairParticle.h"

class TH1;
class AliRsnPairDef;

class AliRsnFunction : public TNamed
{

  public:

    enum EFcnType
    {
      kTrackPt = 0,
      kTrackEta,
      kInvMass,
      kInvMassMC,
      kResolution,
      kPairPt,
      kPairEta,
      kEventMult,
      kFcnTypes
    };

    enum EBinType
    {
      kNoBins = 0,
      kBinPairPt,
      kBinPairEta,
      kBinEventMult,
      kBinTypes
    };

    AliRsnFunction();
    AliRsnFunction(EFcnType fcnType, AliRsnHistoDef *hd);
    AliRsnFunction(EFcnType fcnType, EBinType binType, AliRsnHistoDef *hdMain, AliRsnHistoDef *hdBin);
    AliRsnFunction(EFcnType fcnType, EBinType binType1, EBinType binType2, AliRsnHistoDef *hdMain, AliRsnHistoDef *hdBin1, AliRsnHistoDef *hdBin2);
    AliRsnFunction(const AliRsnFunction &copy);
    virtual ~AliRsnFunction() { delete fHistogram; }
    const AliRsnFunction& operator=(const AliRsnFunction &copy);

    void                 DefineName();

    void                 SetFcnType(EFcnType value) {fFcnType = value;}
    void                 SetPairDef(AliRsnPairDef *def) {fPairDef = def;}
    void                 SetTrack(AliRsnDaughter *track) {fTrack = track;}
    void                 SetPair(AliRsnPairParticle *pair) {fPair = pair;}
    void                 SetEvent(AliRsnEvent *event) {fEvent = event;}
    void                 SetMainHistoDef(AliRsnHistoDef *hd) {fHistoDef[0] = hd;}
    void                 SetPrimaryBinningHistoDef(AliRsnHistoDef *hd) {fHistoDef[1] = hd;}
    void                 SetSecondaryBinningHistoDef(AliRsnHistoDef *hd) {fHistoDef[2] = hd;}

    EFcnType             GetFcnType() {return fFcnType;}
    AliRsnPairDef*       GetPairDef() {return fPairDef;}
    AliRsnDaughter*      GetTrack() {return fTrack;}
    AliRsnPairParticle*  GetPair() {return fPair;}
    AliRsnEvent*         GetEvent() {return fEvent;}
    AliRsnHistoDef*      GetMainHistoDef(AliRsnHistoDef *hd) {return fHistoDef[0];}
    AliRsnHistoDef*      GetPrimaryBinningHistoDef(AliRsnHistoDef *hd) {return fHistoDef[1];}
    AliRsnHistoDef*      GetSecondaryBinningHistoDef(AliRsnHistoDef *hd) { return fHistoDef[2];}

    TH1*                 CreateHistogram(const char *histoName, const char *histoTitle);

    Double_t             Eval();
    Bool_t               Fill();

  protected:

    Int_t       CheckDim();
    Bool_t      CheckInput(Option_t *option);
    const char* BinName(EBinType binType);
    const char* FcnName();
    Double_t    BinValue(EBinType binType);

    EFcnType            fFcnType;     // function type
    EBinType            fBinType[2];  // binning type

    AliRsnPairDef      *fPairDef;     // reference to used pair definition
    AliRsnHistoDef     *fHistoDef[3]; // histogram definition for each axis

    AliRsnDaughter     *fTrack;       // processed track
    AliRsnPairParticle *fPair;        // processed pair
    AliRsnEvent        *fEvent;       // processed event

    TH1                *fHistogram;   // output histogram

    // ROOT dictionary
    ClassDef(AliRsnFunction, 2)
};

#endif
