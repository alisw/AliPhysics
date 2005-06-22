#ifndef ALIJETHEADER_H
#define ALIJETHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
 
//---------------------------------------------------------------------
// Jet header base class 
// Stores a comment which describes the jet analysis
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------
 
#include <TNamed.h>
#include <TString.h>
 
class AliJetHeader : public TNamed
{
 public:
 
  AliJetHeader(const char* name);
  AliJetHeader();
  virtual ~AliJetHeader() { }

  // Getters
  virtual TString GetComment() {return fComment;} 
   
  // Setters
  virtual void SetComment(const char* com) {fComment=TString(com);} 

  // others
  
protected:
  TString fComment; // a comment 

  ClassDef(AliJetHeader,1)
};
 
#endif
