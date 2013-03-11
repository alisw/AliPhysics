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

  void AddBin(const char* particle, const char* type,
              Double_t xmin=TMath::Limits<Double_t>::Max(),
              Double_t xmax=TMath::Limits<Double_t>::Max(),
              Double_t ymin=TMath::Limits<Double_t>::Max(),
              Double_t ymax=TMath::Limits<Double_t>::Max(),
              const char* flavour="");

  TObjArray* CreateParticleArray() const;

  TObjArray* CreateTypeArray() const;
  
  Double_t* CreateBinArray() const;

  TObjArray* CreateBinObjArray() const;
  TObjArray* CreateBinObjArray(const char* particle) const;
  TObjArray* CreateBinObjArray(const char* particle, const char* type, const char* flavour) const;
  
  AliAnalysisMuMuBinning* Project(const char* particle, const char* type, const char* flavour="") const;
  
  virtual void Print(Option_t* opt="") const;
  
  void CreateMesh(const char* particle, const char* type1, const char* type2, const char* flavour="", Bool_t remove12=kFALSE);

  Long64_t Merge(TCollection* list);

  Bool_t IsEqual(const TObject* obj) const;
  
  class Range : public TObject {
    
  public:
    
    Range(const char* particle="",const char* type="",
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

    Bool_t IsNullObject() const;
    
    TString Type() const { return fType; }
    TString Particle() const { return fParticle; }
    Double_t Xmin() const { return fXmin; }
    Double_t Xmax() const { return fXmax; }
    Double_t Ymin() const { return fYmin; }
    Double_t Ymax() const { return fYmax; }
    
    Double_t WidthX() const { return TMath::Abs(fXmin-fXmax); }
    
    Double_t WidthY() const { return TMath::Abs(fYmin-fYmax); }
    
    Bool_t Is2D() const { return fYmax > fYmin; }
    
    const char* GetName() const { return Type().Data(); }
    
    TString AsString() const;
    
    virtual void Print(Option_t* opt="") const;
    
    Bool_t IsInRange(Double_t x, Double_t y=TMath::Limits<Double_t>::Max()) const;
    
    TString Flavour() const { return fFlavour; }
    
  private:
    TString fParticle; // particle
    TString fType; // binning type (e.g. pt, y, phi)
    Double_t fXmin; // x-min of the range
    Double_t fXmax; // x-max of the range
    Double_t fYmin; // x-min of the range
    Double_t fYmax; // x-max of the range
    TString fFlavour; // flavour (if any) this range, e.g. coarse, fine, etc...
    
    ClassDef(AliAnalysisMuMuBinning::Range,2)
  };
  

 private:

  TMap* fBins; // list of bins (particle -> list of bins) = (TObjString -> TObjArray)
  
  ClassDef(AliAnalysisMuMuBinning,1) // custom binning for MuMu analysis
};

#endif
