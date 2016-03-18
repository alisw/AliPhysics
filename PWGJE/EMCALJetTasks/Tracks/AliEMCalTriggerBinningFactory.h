#ifndef ALIEMCALTRIGGERBINNINGFACTORY_H
#define ALIEMCALTRIGGERBINNINGFACTORY_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
namespace EMCalTriggerPtAnalysis {

/**
 * \class AliEMCalTriggerBinningFactory
 * \brief Global binning handler used by several analysis components
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 *
 * This class steers the binning component and set the default binnings for various dimensions,
 * which are used by all analysis components. In case users already defined a binning, this handler
 * does not overwrite this.
 */
class AliEMCalTriggerBinningFactory {
public:
  AliEMCalTriggerBinningFactory();
  /**
   * Destructor, nothing to do
   */
  virtual ~AliEMCalTriggerBinningFactory(){}

  void Create(AliEMCalTriggerBinningComponent * const data);

  static void CreateLinearBinning(TArrayD &binning, int nbins, double min, double max);
protected:
  void CreateMarkusPtBinning(TArrayD &binning) const;
  void CreateRAAPtBinning(TArrayD &binning) const;
  void CreateDefaultEtaBinning(TArrayD& binning) const;
  void CreateDefaultZVertexBinning(TArrayD &binning) const;
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERBINNINGFACTORY_H */
