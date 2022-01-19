/*
Author: Vytautas Vislavicius
Container to store the terms required to calculate higher order moments of the pT spectrum.
Extra layer to calculate the moments. Current implementation supports only the first and second moments.
If used, modified, or distributed, please aknowledge the original author of this code.
*/
#ifndef ALICKCONTAINER__H
#define ALICKCONTAINER__H
#include "TCollection.h"
#include "AliProfileBS.h"
#include "TNamed.h"
#include "TH1.h"
class AliCkContainer: public TNamed {
  public:
    AliCkContainer();
    ~AliCkContainer();
    AliCkContainer(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins);
    AliCkContainer(const char* name, const char* title, Int_t nbinsx, Double_t xmin, Double_t xmax);
    Long64_t Merge(TCollection *collist);
    TList *getObsList() { return fObsList; };
    void Initialize(Int_t nbinsx, const Double_t* xbins);
    void Initialize(Int_t nbinsx, Double_t xmin, Double_t xmax);
    void InitializeSubsamples(const Int_t &nrb);
    void FillObs(const Double_t inArr[5], const Double_t &lMult, const Double_t &rn=0);
    TH1 *RecalculateCkHists(TH1 **inh);
    void CalculateCks();
    TList *getCkList() { return fCkList; };
    TH1 *getHist(Int_t ind=-1);
    void RebinMulti(Int_t nbins) {if(!fObsList) return; for(Int_t i=0;i<fObsList->GetEntries();i++) ((AliProfileBS*)fObsList->At(i))->RebinMulti(nbins); };
    void RebinMulti(Int_t nbins, Double_t *binedges) {if(!fObsList) return; for(Int_t i=0;i<fObsList->GetEntries();i++) ((AliProfileBS*)fObsList->At(i))->RebinMulti(nbins, binedges); };
  protected:
    TList *fObsList;
    TList *fCkList; //!
  ClassDef(AliCkContainer,1);
};
#endif
