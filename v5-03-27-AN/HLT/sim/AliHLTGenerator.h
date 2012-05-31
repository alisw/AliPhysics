// @(#) $Id$

#ifndef ALIL3GENERATOR_H
#define ALIL3GENERATOR_H

#include "AliHLTRootTypes.h"
#include "AliGenerator.h"

class AliHLTGenerator : public AliGenerator {

 private:
  
  
  

 public:
  AliHLTGenerator();
  virtual ~AliHLTGenerator();
  
  void Generate();
  void Init();
  void ReadParticles(TClonesArray *particles);

  ClassDef(AliHLTGenerator,1)
};

typedef AliHLTGenerator AliL3Generator; // for backward compatibility
