#ifndef ALIEMCALTRIGGERBINNINGFACTORY_H
#define ALIEMCALTRIGGERBINNINGFACTORY_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerBinningFactory {
public:
  AliEMCalTriggerBinningFactory();
  virtual ~AliEMCalTriggerBinningFactory(){}

  void Create(AliEMCalTriggerBinningComponent * const data);

protected:
  void CreateDefaultPtBinning(TArrayD &binning) const;
  void CreateDefaultEtaBinning(TArrayD& binning) const;
  void CreateDefaultZVertexBinning(TArrayD &binning) const;
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERBINNINGFACTORY_H */
