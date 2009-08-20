#ifndef AliTPCConfigParser_H
#define AliTPCConfigParser_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
#include <TObject.h>

class TObjArray;
class TList;

class AliTPCConfigParser : public TObject{
  
public:
  AliTPCConfigParser();
  AliTPCConfigParser(const char* cfgfile);
  AliTPCConfigParser(const AliTPCConfigParser &cfg);
  AliTPCConfigParser& operator = (const AliTPCConfigParser &cfg);
 
  virtual ~AliTPCConfigParser();
  
  Int_t ParseConfigFileTxt(const char* cfgfile);
  Float_t GetValue(const char* name, UInt_t position=0); 

  const TList* GetConfigurationMap() const {return fConfigMap;}
  void ResetMap();

  const TList*    operator()() const {return fConfigMap;}
  const TObject* operator()(const char* key) const {return fConfigMap->FindObject(key);}
  const TObject* operator()(TObject* key) const {return fConfigMap->FindObject(key);}
  TObject* NextKey();
  TObject* NextValue(const char *key);
  TObject* NextValue(TObject *key);

private:
  TList *fConfigMap;
  TIterator *fKeyIter;
  TIterator *fValIter;

  TObject *NextValueIter(TObjArray *obj);
  
  ClassDef(AliTPCConfigParser, 1)         // TPC DA configuration file parser
};
#endif
