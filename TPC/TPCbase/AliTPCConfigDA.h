#ifndef ALITPCCONFIGDA_H
#define ALITPCCONFIGDA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
#include <TObject.h>
////////////////////////////////////////////////////////////////////////////
//
// Simple configuration file parser
//
////////////////////////////////////////////////////////////////////////////

class TMap;

class AliTPCConfigDA : public TObject{
  
public:
  AliTPCConfigDA();
  AliTPCConfigDA(const char* cfgfile);
  AliTPCConfigDA(const AliTPCConfigDA &cfg);
  AliTPCConfigDA& operator = (const AliTPCConfigDA &cfg);
 
  virtual ~AliTPCConfigDA();
  
  Int_t ParseConfigFileTxt(const char* cfgfile);
  Float_t GetValue(const char* name) const; 

  const TMap* GetConfigurationMap() const {return fConfigMap;}
  void ResetMap();

private:
  TMap *fConfigMap;                   // Configuration map
  
  ClassDef(AliTPCConfigDA, 1)         // TPC DA configuration file parser
};
#endif
