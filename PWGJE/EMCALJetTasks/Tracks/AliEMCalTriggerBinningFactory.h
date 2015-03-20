#ifndef ALIEMCALTRIGGERBINNINGFACTORY_H
#define ALIEMCALTRIGGERBINNINGFACTORY_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
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
  void CreateLinearBinning(TArrayD &binning, int nbins, double min, double max) const;
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERBINNINGFACTORY_H */
