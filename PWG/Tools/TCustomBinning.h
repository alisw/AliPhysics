#ifndef TCUSTOMBINNING_H
#define TCUSTOMBINNING_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TArrayD;

#include <Rtypes.h>
#include <TBinning.h>
#include <exception>
#include <map>

/**
 * @class TCustomBinning
 * @brief Helper class creating user defined custom binning.
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @ingroup Histmanager
 * @since May 31st, 2016
 *
 * A custom binning is defined by a user an can contain bins of different width. The
 * user has to define ranges in which the bin width is equal, bin edges are calculated
 * accordingly. Do do so, the first step is to define a minimum for the range.
 *
 * ~~~{.cxx}
 * TCustomBinning mybinning;
 * mybinning.SetMinimum(-100.);
 * ~~~
 *
 * Now ranges can be added via the AddStep function:
 *
 * ~~~{.cxx}
 * mybinning.AddStep(-10., 10.);
 * mybinning.AddStep(10., 1.);
 * mybinning.AddStep(100., 10);
 * ~~~
 *
 * In this example the binning is done:
 * - from -100 to -10 in steps of 10
 * - from -10 to 10 in steps of 1
 * - from 10 to 100 in steps of 10
 *
 * The binning can be converted to a TArrayD which contains the bin edges in increasing order:
 *
 * ~~~{.cxx}
 * TArrayD binlimits;
 * mybinning.CreateBinEdges(binlimits);
 * ~~~
 *
 * @note In case the binning is used together with the THistManager the last step is done by the
 * THistManager and does not need to be performed by the user.
 */
class TCustomBinning : public TBinning {
public:

  /**
   * @class MinNotSetException
   * @brief Exception thrown in case the minimum is not set
   * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
   * @ingroup Histmanager
   * @since May 31st, 2016
   *
   * Exception thrown in case the minimum is not set. Does not need to carry any information.
   * Only implements the what function.
   */
  class MinNotSetException : public std::exception {
  public:

    /**
     * Constructor
     */
    MinNotSetException() {}

    /**
     * Destructor
     */
    virtual ~MinNotSetException() throw() {}

    /**
     * Provide error message, saying that the minimum is not set
     * @return Error message
     */
    virtual const char *what() const throw() { return "Minimum of the binning not set"; }
  };

  /**
   * Constructor
   */
  TCustomBinning();

  /**
   * Destructor
   */
  virtual ~TCustomBinning() {}

  /**
   * Set the minumum (the lowest bin edge). As the
   * bin edges are calculated based on the minumum this
   * information is cructial to be provided by the user.
   * @param[in] minimum bin edge
   */
  void SetMinimum(Double_t min) { fMinimum = min; fMinimumSet = true; }

  /**
   * Add new step. Will create bins in steps of the binwidth
   * from the previous maximum to the maximum assigned here.
   * @param[in] max New maximum of the binning
   * @param[in] binwidth Width of the bin
   */
  void AddStep(Double_t max, Double_t binwidth);

  /**
   * Create bin edges from the minimum and the ranges. Bin edges
   * will be stored in the array edges in increasing order
   * @param[out] edges Array storing bin edges
   * @throw MinNotSetException in case the minimum is not defined
   */
  virtual void CreateBinEdges(TArrayD &edges) const;

private:
  Double_t                        fMinimum;           ///< Minimum of the binning
  Bool_t                          fMinimumSet;        ///< Define whether minimum is set. Attention: Bin edges will not be created without minimum
  std::map<double, double>        fSteps;             ///< List of ranges with common bin width
};

#endif /* TCUSTOMBINNING_H */
