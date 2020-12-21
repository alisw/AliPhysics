#ifndef TBINNING_H
#define TBINNING_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

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
class TBinning : public TObject {
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
   * virtual function creating a copy of the current binning
   * @return Copy of the current binning
   */
  virtual TBinning * MakeCopy() const = 0;

  /**
   * Function creating bin edges from a descripition
   * To be implemented by classes inheriting from TBinning
   * @param [out] binedges Target array of bin edges
   */
  virtual void CreateBinEdges(TArrayD &binedges) const = 0;

  ClassDef(TBinning, 1);
};

#endif // TBINNING_H
