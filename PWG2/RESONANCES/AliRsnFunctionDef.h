//
// Class AliRsnFunctionDef
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

#ifndef ALIRSNFUNCTIONDEF_H
#define ALIRSNFUNCTIONDEF_H

#include <TArrayD.h>
#include <TString.h>

class AliRsnFunctionDef : public TObject
{
  public:

    enum EFcnType
    {
      kInvMass,
      kInvMassMC,
      kInvMassRotated,
      kResolution,
      kPtSpectrum,
      kEtaSpectrum,
      kFcnTypes
    };

    enum EFcnBinType
    {
      kPt,            // binning in Pt of the pair
      kEta,           // binning in Eta of the pair
      kBinningTypes
    };

    enum EConstValue
    {
      kBinMax = 5
    };

    AliRsnFunctionDef();
    AliRsnFunctionDef(EFcnType type, Int_t nbins, Double_t min, Double_t max);
    AliRsnFunctionDef(const AliRsnFunctionDef &copy);
    const AliRsnFunctionDef& operator=(const AliRsnFunctionDef &copy);
    virtual ~AliRsnFunctionDef() { }

    // names
    TString  GetFcnName();
    TString  GetFcnTitle();
    TString  GetBinningName(Int_t i);

    // main histogram definition
    EFcnType GetType() {return fFcnType;}
    Int_t    GetNBinsX() {return fNBinsX;}
    Double_t GetXmin() {return fXmin;}
    Double_t GetXmax() {return fXmax;}
    void     SetBinningX(Int_t nbins, Double_t xmin, Double_t xmax) {fNBinsX=nbins;fXmin=xmin;fXmax=xmax;}

    // secondary binning definition
    Int_t    GetMaxYBinUsed() {return fYUsed;}
    Int_t    GetNBinsY(Int_t i) {if (i<=fYUsed) return fNBinsY[i]; else return 0;}
    TArrayD* GetYbins(Int_t i) {if (i<=fYUsed) return &fYbins[i]; else return 0x0;}
    EFcnBinType GetYtype(Int_t i) {if (i<=fYUsed) return fBinType[i]; else return (EFcnBinType)0;}
    void     SetBinningY(Int_t i, EFcnBinType type, Double_t xmin, Double_t xmax, Double_t step);
    void     SetBinningY(Int_t i, EFcnBinType type, Int_t nbins, Double_t *bins);
    Bool_t   AddBinningY(EFcnBinType type, Double_t min, Double_t max, Double_t step);
    Bool_t   AddBinningY(EFcnBinType type, Int_t nbins, Double_t *bins);

    // rotation angle
    void     SetRotationAngle(Double_t rotAngle) {fRotAngle = rotAngle;}

    // computations
    Double_t EvalX(AliRsnPairParticle *pair, AliRsnPairDef *ref);
    Double_t EvalY(Int_t i, AliRsnPairParticle *pair);

    // other
    void     Print(Option_t *opt = "");

  private:

    EFcnType     fFcnType;              // function type

    Int_t        fNBinsX;               // primary histogram bins
    Double_t     fXmin;                 // primary histogram low edge
    Double_t     fXmax;                 // primary histogram up edge

    Int_t        fYUsed;                // index of last secondary binning used
    EFcnBinType  fBinType[kBinMax]; // binning type
    Int_t        fNBinsY[kBinMax];      // secondary bins number
    TArrayD      fYbins[kBinMax];       // secondary histogram bins

    Double_t         fRotAngle;             // rotation angle (for "rotated" invMass)

    ClassDef(AliRsnFunctionDef,1)
};

#endif

