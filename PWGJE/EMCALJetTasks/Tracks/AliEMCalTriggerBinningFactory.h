#ifndef ALIEMCALTRIGGERBINNINGFACTORY_H
#define ALIEMCALTRIGGERBINNINGFACTORY_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TLinearBinning.h>
#include <TCustomBinning.h>

namespace PWGJE {
  
namespace EMCALJetTasks {

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
   * @class DefaultZVertexBinning
   * @brief Creating default z-Vertex binning. Bin size 5 cm.
   */
  class DefaultZVertexBinning : public TLinearBinning {
  public:
    DefaultZVertexBinning(): TLinearBinning(4, -10, 10) {}
    virtual ~DefaultZVertexBinning() {}
  };

  /**
   * @class DefaultEtaBinning
   * @brief Creating default \f$ \eta \$f  binning. Bin size fixed at 0.1 units of rapidity.
   */
  class DefaultEtaBinning : public TLinearBinning {
  public:
    DefaultEtaBinning() : TLinearBinning(16, -0.8, 0.8) {}
    virtual ~DefaultEtaBinning() {}
  };

  /**
   * @class DefaultPtBinning
   * @brief Create \f$ p_{t} \f$ binning used in the \f$ R_{AA} \f$ analysis:
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
   */
  class DefaultPtBinning : public TCustomBinning {
  public:
    DefaultPtBinning();
    virtual ~DefaultPtBinning() {}
  };

  /**
   * @class MarkusPtBinning
   * @brief Creating the default \f$ p_{t} \f$ binning.
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
   */
  class MarkusPtBinning : public TCustomBinning {
  public:
    MarkusPtBinning();
    virtual ~MarkusPtBinning() {}
  };

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
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERBINNINGFACTORY_H */
