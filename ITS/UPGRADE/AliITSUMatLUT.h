#ifndef ALIITSUMATLUT
#define ALIITSUMATLUT

#include <TObject.h>
class TH1;

class AliITSUMatLUT : public TObject
{
 public:
  enum {kParX2X0,kParRhoL,kNParTypes};
  //
  AliITSUMatLUT();
  AliITSUMatLUT(Double_t rmin,Double_t rmax,Int_t nbin);
  AliITSUMatLUT(const AliITSUMatLUT& src);
  AliITSUMatLUT &operator=(const AliITSUMatLUT& src);
  //
  virtual ~AliITSUMatLUT();
  virtual void Print(Option_t* option = "") const;
  //
  void   FillData(Int_t ntest, Double_t zmin,Double_t zmax);
  TH1*   GetHisto(const Option_t* option="", const Char_t *name=0) const;
  //
  Double_t GetMatBudget(const Double_t *pnt0, const Double_t *pnt1, Double_t *ret) const;
  Double_t GetData(Int_t parTyp, Int_t bin)     const {return fData[parTyp][bin];}
  Double_t GetDataDiff(Int_t parTyp, Int_t bin) const {return fData[parTyp][bin] - (bin ? fData[parTyp][bin-1] : 0);}  
  void    GetData(double r, double* dest)        const;
  //
 protected:
  Double_t  fRMin;             // min radius
  Double_t  fRMax;             // max radius
  Double_t  fDRInv;            // inverse bin size
  Double_t  fDR;               // bin size
  Int_t    fNBins;            // number of bins
  Double_t* fData[kNParTypes]; //[fNBins] array per type
  //
  ClassDef(AliITSUMatLUT,1)
};


#endif
