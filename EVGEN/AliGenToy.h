#ifndef ALIGENTOY_H
#define ALIGENTOY_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class AliVParticle;
class AliGenToyEventHeader;

#include "AliGenerator.h"

class AliGenToy : public AliGenerator
{
 public:
  AliGenToy(const std::string &name = "toy MC");
  virtual ~AliGenToy();
  virtual void Init();
  virtual void Generate();

  virtual void UserGenerate() = 0;

 protected:
  // AddParticle returns label
  Int_t AddParticle(const TLorentzVector &part);
  Int_t AddParticle(const AliVParticle &part);
  Int_t AddParticle(Double_t px, Double_t py, Double_t pz, Int_t pdg);

  void SetCentrality(Double_t cent);
  void SetVertex(Double_t vx, Double_t vy, Double_t vz);
  void SetValue(const std::string &key, Double_t value);

  AliGenToyEventHeader *fHeader;

  Int_t   fNProduced;
  Float_t fEventWeight;
  TArrayF fEventVertex;

 private:
  AliGenToy(const AliGenToy &para);
  AliGenToy& operator = (const AliGenToy &para);

  ClassDef(AliGenToy, 1) // toy generator
};
#endif
