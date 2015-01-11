#ifndef ALIITSUPARAMLIST_H
#define ALIITSUPARAMLIST_H

#include "AliParamList.h"
#include <TObjArray.h>

class AliITSUParamList : public AliParamList
{
 public:
  AliITSUParamList(Int_t n=0, const Double_t *parVal=0);
  AliITSUParamList(const AliITSUParamList& src);
  AliITSUParamList& operator=(const AliITSUParamList& src);
  virtual ~AliITSUParamList();
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
  ClassDef(AliITSUParamList,1)
};


#endif
