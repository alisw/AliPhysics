// @(#) $Id$

#ifndef ALIL3GENERATOR_H
#define ALIL3GENERATOR_H

#include "AliL3RootTypes.h"
#include "AliGenerator.h"

class AliL3Generator : public AliGenerator {

 private:
  
  
  

 public:
  AliL3Generator();
  virtual ~AliL3Generator();
  
  void Generate();
  void Init();
  void ReadParticles(TClonesArray *particles);

  ClassDef(AliL3Generator,1)
};
