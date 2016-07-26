#ifndef ALIEMCALTRIGGERBINNINGFACTORY_H
#define ALIEMCALTRIGGERBINNINGFACTORY_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

namespace EMCalTriggerPtAnalysis {

/**
 * @class AliEMCalTriggerBinningFactory
 * @brief Global binning handler used by several analysis components
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @date Dec 12, 2014
 *
 * This class steers the binning component and set the default binnings for various dimensions,
 * which are used by all analysis components. In case users already defined a binning, this handler
 * does not overwrite this.
 */
class AliEMCalTriggerBinningFactory {
public:
  /**
   * Default constructor, nothing to do
   */
  AliEMCalTriggerBinningFactory();
  /**
   * Destructor, nothing to do
   */
  virtual ~AliEMCalTriggerBinningFactory(){}

  /**
   * Initialise binning component with default binning
   * @param[out] data the binning component to be initialised
   */
  void Create(AliEMCalTriggerBinningComponent * const data);

  /**
   * Create any kind of linear binning from given ranges and stores it in the binning array.
   * @param[out] binning output array
   * @param[in] nbins Number of bins
   * @param[in] min lower range
   * @param[in] max upper range
   */
  static void CreateLinearBinning(TArrayD &binning, int nbins, double min, double max);
protected:
  /**
   * Creating the default \f$ p_{t} \f$ binning.
   *
   * Definition used:
   * - from 0 to 2.5 GeV/c: 0.1 GeV/c bins
   * - from 2.5 to 7 GeV/c: 0.25 GeV/c bins
   * - from 7 to 10 GeV/c: 0.5 GeV/c bins
   * - from 10 to 15 GeV/c: 1 GeV/c bins
   * - from 15 to 20 GeV/c: 2.5 GeV/c bins
   * - from 20 to 30 GeV/c: 5 GeV/c bins
   * - from 30 to 100 GeV/c: 10 GeV/c bins
   * - from 100 to 200 GeV/c: 20 GeV/c bins
   *
   * @param[out] binning Array where to store the results.
   */
  void CreateMarkusPtBinning(TArrayD &binning) const;
  /**
   * Create \f$ p_{t} \f$ binning used in the \f$ R_{AA} \f$ analysis:
   *
   * Definitions are:
   * - from 0.15 to 1 GeV/c: 0.05 GeV/c bins
   * - from 1 to 2 GeV/c: 0.1 GeV/c bins
   * - from 2 to 4 GeV/c: 0.2 GeV/c bins
   * - from 4 to 7 GeV/c: 0.5 GeV/c bins
   * - from 7 to 16 GeV/c: 1 GeV/c bins
   * - from 16 to 36 GeV/c: 2 GeV/c bins
   * - from 36 to 40 GeV/c: 4 GeV/c bins
   * - from 40 to 50 GeV/c: 5 GeV/c bins
   * - from 50 to 100 GeV/c: 10 GeV/c bins
   *
   * @param[out] binning Array where to store the results
   */
  void CreateRAAPtBinning(TArrayD &binning) const;
  /**
   * Creating default \f$ \eta \$f  binning. Bin size fixed at 0.1 units of rapidity.
   * @param[out] binning Array where to store the results.
   */
  void CreateDefaultEtaBinning(TArrayD& binning) const;
  /**
   * Creating default z-Vertex binning. Bin size 5 cm.
   * @param[out] binning Array where to store the results.
   */
  void CreateDefaultZVertexBinning(TArrayD &binning) const;
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERBINNINGFACTORY_H */
