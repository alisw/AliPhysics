#ifndef TBINNING_H
#define TBINNING_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
 * @class TBinning
 * @brief Interface for binnings used by the histogram handler
 * @author Markus Fasel <markus.fasel@cern.ch>,
 * @ingroup Histmanager
 * @since May 31st, 2016
 *
 * This class is the base class for binning descriptions used by the histogram
 * handler during the creation of histograms. Classes implementing the binning
 * must implement the function
 *
 * ~~~{.cxx}
 * void CreateBinEdges(TArrayD &binedges) const.
 * ~~~
 */
class TBinning {
public:

  /**
   * Constructor
   */
  TBinning() {}

  /**
   * Destructor
   */
  virtual ~TBinning() {}

  /**
   * Function creating bin edges from a descripition
   * To be implemented by classes inheriting from TBinning
   * @param [out] binedges Target array of bin edges
   */
  virtual void CreateBinEdges(TArrayD &binedges) const = 0;
};

#endif // TBINNING_H
