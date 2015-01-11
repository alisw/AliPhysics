#ifndef ALIITSQASSDREFDATA_H
#define ALIITSQASSDREFDATA_H
/* Copyright(c) 2009-2011, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliITSQASSDRefData
//                     ITS SSD reference values for the QA
//
//         Origin: Panos.Christakoglou@cern.ch, NIKHEF-Utrecht University
//-------------------------------------------------------------------------


#include <TObject.h>
#include <TObjArray.h>
#include <TArrayD.h>


class AliITSQASSDRefData : public TObject {
  
 public:
  AliITSQASSDRefData(); 
  AliITSQASSDRefData(Int_t specie); 
  AliITSQASSDRefData(const char* path);
  AliITSQASSDRefData(const AliITSQASSDRefData& refData);   
  AliITSQASSDRefData& operator = (const AliITSQASSDRefData& refData);
  virtual ~AliITSQASSDRefData(); //destructor
  
  void AddReference(const char* name, Int_t id, Double_t value);
  Int_t GetID(const char*);
  
  Double_t *GetReferenceData() {return fRefList->GetArray();}
  Double_t GetReferenceValue(const char*);
  Double_t GetReferenceValue(Int_t id);

  void SetDefault(Int_t eventSpecie);

  void SetReferenceData(const char* path);
  void SetReferenceValue(const char* name, Double_t value);
  void SetReferenceValue(Int_t id, Double_t value);
  
  void PrintTable();
  
 private:
  TArrayD *fRefList;//* = new TArrayD(11,{0,500,0,50,0,100,0,50,0,100,5});
  TObjArray *fNameList;//* = new TObjArray(11);
  
  ClassDef(AliITSQASSDRefData,1)           // description 
};


#endif
