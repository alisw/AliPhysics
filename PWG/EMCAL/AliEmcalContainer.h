#ifndef AliEmcalContainer_H
#define AliEmcalContainer_H

//
// container with name, TClonesArray
//

class TLorentzVector;
class AliVEvent;

#include <TNamed.h>
#include <TClonesArray.h>

class AliEmcalContainer : public TNamed {
 public:
  AliEmcalContainer();
  AliEmcalContainer(const char *name); 
  virtual ~AliEmcalContainer(){;}

  void SetArrayName(const char *n)                 {fClArrayName = n;}
  TClonesArray               *GetArray()           {return fClArray;}
  Int_t                       GetNEntries() const  {return fClArray->GetEntriesFast();}
  const TString&              GetArrayName() const {return fClArrayName;}
  void                        SortArray() {fClArray->Sort();}
  virtual void                GetMomentum(TLorentzVector &mom, Int_t i) const = 0;

 protected:
  void SetArray(AliVEvent *event, const char *clname=0);

  TClonesArray               *fClArray;                 //!TClonesArray
  Double_t                    fVertex[3];               //!event vertex array
  TString                     fClArrayName;             // name of branch

 private:
  AliEmcalContainer(const AliEmcalContainer& obj); // copy constructor
  AliEmcalContainer& operator=(const AliEmcalContainer& other); // assignment

  ClassDef(AliEmcalContainer,1);

};

#endif

