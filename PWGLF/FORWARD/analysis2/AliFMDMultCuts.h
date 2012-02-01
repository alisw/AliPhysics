#ifndef ALIFMDMULTCUTS_H
#define ALIFMDMULTCUTS_H
#include <TObject.h>

class AliFMDMultCuts : public TObject 
{
public:
  AliFMDMultCuts();
  AliFMDMultCuts(const AliFMDMultCuts& o);
  AliFMDMultCuts& operator=(const AliFMDMultCuts& o);
  Double_t GetMultCut(UShort_t d, Char_t r, Double_t eta, Bool_t errors) const;
  Double_t GetMultCut(UShort_t d, Char_t r, Int_t etabin, Bool_t errors) const;
  
  void UnsetMultCuts() { SetMultCuts(-1); }
  void SetMultCuts(Double_t fmd1i, 
		   Double_t fmd2i=-1, 
		   Double_t fmd2o=-1, 
		   Double_t fmd3i=-1, 
		   Double_t fmd3o=-1);
  void SetMPVFraction(Double_t frac=0) { fMPVFraction = frac; }
  void SetNXi(Double_t nXi) { fNXi = nXi; }
  void SetIncludeSigma(Bool_t in) { fIncludeSigma = in; }
  void Print(Option_t* option="") const;
  void Output(TList* l, const char* name=0) const;
  Double_t GetFixedCut(UShort_t d, Char_t r) const;
protected:
  Double_t fMultCuts[5];
  Double_t fMPVFraction;
  Double_t fNXi;
  Bool_t   fIncludeSigma;
  
  ClassDef(AliFMDMultCuts,1); // Cuts on ESD Mult 
};

#endif
// Local Variables: 
//  mode: C++
// End:
