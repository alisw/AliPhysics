#ifndef AliEmcalContainer_H
#define AliEmcalContainer_H

//
// container with name, TClonesArray
//

class AliVEvent;

#include <iostream>
#include "Rtypes.h"
#include <TArrayS.h>
#include "TString.h"
#include "TNamed.h"
#include "TClonesArray.h"

class AliEmcalContainer : public TNamed {
 public:
  AliEmcalContainer();
  AliEmcalContainer(const char *name); 
  virtual ~AliEmcalContainer();

  void SetArrayName(const char *n)                 {fClArrayName = n;}
  TClonesArray               *GetArray()           {return fClArray;}
  Int_t                       GetNEntries() const  {return fClArray->GetEntriesFast();}
  TString                     GetArrayName()       {return fClArrayName;}

 protected:
  void SetArray(AliVEvent *event, const char *clname=0);

  TClonesArray               *fClArray;                 //!TClonesArray
  TString                     fClArrayName;             // name of branch

 private:
  AliEmcalContainer(const AliEmcalContainer& obj); // copy constructor
  AliEmcalContainer& operator=(const AliEmcalContainer& other); // assignment

  ClassDef(AliEmcalContainer,1);

};

#endif

