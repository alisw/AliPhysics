#ifndef ALICOMPARISONOBJECT_H
#define ALICOMPARISONOBJECT_H

//------------------------------------------------------------------------------
// Abstract class to keep information from comparison of 
// reconstructed and MC particle tracks.   
// 
// Author: J.Otwinowski 04/14/2008 
//------------------------------------------------------------------------------

#include "TNamed.h"
#include "TFolder.h"

class AliMCInfo;
class AliESDRecInfo;

class AliComparisonObject : public TNamed {
public :
  AliComparisonObject(); 
  AliComparisonObject(const char* name="AliComparisonObject"); 
  virtual ~AliComparisonObject();

  // Init data members
  // call once before event loop
  virtual void Init() = 0;

  // Execute analysis
  // call in the event loop 
  virtual void Exec(AliMCInfo* infoMC=0, AliESDRecInfo *infoRC=0) = 0;

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list=0) = 0;

  // Analyse output histograms
  virtual void Analyse() = 0;

  // Get output folder for analysed histograms
  virtual TFolder* GetAnalysisFolder() = 0;
  //virtual TFolder* CreateFolder(TString name="",TString title="") = 0;

  ClassDef(AliComparisonObject,1);
};

#endif
