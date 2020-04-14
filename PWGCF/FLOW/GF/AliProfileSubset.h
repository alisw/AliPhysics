/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572)
*/
#ifndef ALIPROFILESUBSET__H
#define ALIPROFILESUBSET__H
//Helper function to select a subrange of a TProfile
#include "TProfile.h"
#include "TProfile2D.h"
#include "TError.h"
class AliProfileSubset : public TProfile2D {
  public:
    AliProfileSubset() {};
    AliProfileSubset(TProfile2D &inpf);
    TProfile *GetSubset(Bool_t onx, const char *name, Int_t fb, Int_t lb, Int_t l_nbins=0, Double_t *l_binarray=0);
    ~AliProfileSubset() {};
    void OverrideBinContent(Double_t x, Double_t y, Double_t x2, Double_t y2, Double_t val);
    void OverrideBinContent(Double_t x, Double_t y, Double_t x2, Double_t y2, TProfile2D *sourceProf);
    Bool_t OverrideBinsWithZero(Int_t xb1, Int_t yb1, Int_t xb2, Int_t yb2);
    ClassDef(AliProfileSubset,1);
};
#endif
