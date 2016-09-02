#ifndef ALIANALYSISMUMUBINNING_H
#define ALIANALYSISMUMUBINNING_H

#include "TNamed.h"
#include "TString.h"
#include "TMath.h"

//
// AliAnalysisMuMuBinning :
// helper class to store the various bins used
//
// author: L. Aphecetche (Subatech)
//

class TMap;
class TObjArray;


class AliAnalysisMuMuBinning : public TNamed
{
public:

  class Range;

  AliAnalysisMuMuBinning(const char* name="", const char* title="");
  AliAnalysisMuMuBinning(const AliAnalysisMuMuBinning& rhs);
  AliAnalysisMuMuBinning& operator=(const AliAnalysisMuMuBinning& rhs);

  virtual ~AliAnalysisMuMuBinning();

  void AddBin(const AliAnalysisMuMuBinning::Range& bin);

  void AddBin(const char* what, const char* quantity,
              Double_t xmin=TMath::Limits<Double_t>::Max(),
              Double_t xmax=TMath::Limits<Double_t>::Max(),
              const char* flavour="")
  { AddBin(what,quantity,xmin,xmax,TMath::Limits<Double_t>::Max(),TMath::Limits<Double_t>::Max(),flavour); }

  void AddBin(const char* what, const char* quantity,
              Double_t xmin,
              Double_t xmax,
              Double_t ymin,
              Double_t ymax,
              const char* flavour="");

  TObjArray* CreateWhatArray() const;

  TObjArray* CreateQuantityArray() const;

  Double_t* CreateBinArray() const;
  Double_t* CreateBinArrayY() const;
  Double_t* CreateBinArrayX() const;

  TObjArray* CreateBinObjArray() const;
  TObjArray* CreateBinObjArray(const char* what) const;
  TObjArray* CreateBinObjArray(const char* what, const char* quantity, const char* flavour) const;

  TObjArray* CreateBinStrArray() const;
  TObjArray* CreateBinStrArray(const char* what) const;
  TObjArray* CreateBinStrArray(const char* what, const char* quantity, const char* flavour) const;

  Int_t GetNBinsX() const;
  Int_t GetNBinsY() const;

  AliAnalysisMuMuBinning* Project(const char* what, const char* quantity, const char* flavour="") const;

  virtual void Print(Option_t* opt="") const;

  void CreateMesh(const char* what, const char* quantity1, const char* quantity2, const char* flavour="", Bool_t remove12=kFALSE);

  Long64_t Merge(TCollection* list);

  Bool_t IsEqual(const TObject* obj) const;

  class Range : public TObject {

  public:

    Range(const char* what="",const char* quantity="",
          Double_t xmin=TMath::Limits<Double_t>::Max(),
          Double_t xmax=TMath::Limits<Double_t>::Max(),
          Double_t ymin=TMath::Limits<Double_t>::Max(),
          Double_t ymax=TMath::Limits<Double_t>::Max(),
          const char* version="");

    virtual Int_t	Compare(const TObject* obj) const;
    Bool_t IsEqual(const TObject* obj) const { return Compare(obj)==0; }
    Bool_t IsSortable() const { return kTRUE; }

    virtual TObject* Clone(const char* /*newname*/ = "") const { return new Range(*this); }

    bool operator==(const Range& other) const { return Compare(&other)==0; }

    bool operator!=(const Range& other) const { return !(*this==other); }

    Bool_t IsIntegrated() const;

    TString Quantity() const { return fQuantity; }
    TString What() const { return fWhat; }
    Double_t Xmin() const { return fXmin; }
    Double_t Xmax() const { return fXmax; }
    Double_t Ymin() const { return fYmin; }
    Double_t Ymax() const { return fYmax; }

    Double_t WidthX() const { return TMath::Abs(fXmin-fXmax); }

    Double_t WidthY() const { return TMath::Abs(fYmin-fYmax); }

    Bool_t Is2D() const { return fYmax > fYmin; }

    const char* GetName() const { return What().Data(); }

    TString AsString() const;

    virtual void Print(Option_t* opt="") const;

    Bool_t IsInRange(Double_t x, Double_t y=TMath::Limits<Double_t>::Max()) const;

    TString Flavour() const { return fFlavour; }

  private:
    TString fWhat; // what this range is about (e.g. J/psi particle, event, etc...)
    TString fQuantity; // binning type (e.g. pt, y, phi)
    Double_t fXmin; // x-min of the range
    Double_t fXmax; // x-max of the range
    Double_t fYmin; // x-min of the range
    Double_t fYmax; // x-max of the range
    TString fFlavour; // flavour (if any) this range, e.g. coarse, fine, etc...

    ClassDef(AliAnalysisMuMuBinning::Range,3)
  };


 private:

  TMap* fBins; // list of bins (what -> list of bins) = (TObjString -> TObjArray)

  ClassDef(AliAnalysisMuMuBinning,2) // custom binning for MuMu analysis
};

#endif
