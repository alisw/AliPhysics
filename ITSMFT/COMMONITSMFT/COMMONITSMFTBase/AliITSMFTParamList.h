#ifndef AliITSMFTPARAMLIST_H
#define AliITSMFTPARAMLIST_H

#include "AliParamList.h"
#include <TObjArray.h>

class AliITSMFTParamList : public AliParamList
{
 public:
  AliITSMFTParamList(Int_t n=0, const Double_t *parVal=0);
  AliITSMFTParamList(const AliITSMFTParamList& src);
  AliITSMFTParamList& operator=(const AliITSMFTParamList& src);
  virtual ~AliITSMFTParamList();
  //
  TObjArray*  GetParamObjects()                const {return fParamObj;}
  TObject*    GetParamObject(const char* name) const {return fParamObj ? fParamObj->FindObject(name):0;}
  void        AddParamObject(TObject* obj);
  //
  virtual void  Print(Option_t *opt="") const;
  //
 protected:
  TObjArray* fParamObj;            // optional array of parameters objects
  //
  ClassDef(AliITSMFTParamList,1)
};


#endif
