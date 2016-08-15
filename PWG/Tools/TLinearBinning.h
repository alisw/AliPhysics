#ifndef TLINEARBINNING_H
#define TLINEARBINNING_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TArrayD;

#include <exception>
#include <Rtypes.h>
#include <TBinning.h>

/**
 * @class TLinearBinning
 * @brief Class creating a linear binning, used in the histogram manager
 * @author Markus Fasel
 * @ingroup Histmanager
 * @since May 31st, 2016
 *
 * This class creates a linear binning. For this the user must provide
 * - A minimum
 * - A maximum
 * - The number of bins
 *
 * The information can be set either in the constructor
 *
 * ~~~{.cxx}
 * TLinearBinning mybinning(100, -10., 10.);
 * ~~~
 *
 * or using the set function
 *
 * ~~~{.cxx}
 * TLinearBinning mybinning;
 * mybinning.Set(100, -10, 10);
 * ~~~
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
class TLinearBinning : public TBinning{
public:

  /**
   * @class LimitsNotSetException
   * @brief Exception indicating that the limits are not set
   * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
   * @ingroup Histmanager
   * @since May 31st, 2016
   *
   * Dummy exception thrown in case the limits are not set. Only creating error message
   */
  class LimitsNotSetException : public std::exception{
  public:

    /**
     * Constructor
     */
    LimitsNotSetException() {}

    /**
     * Destructor
     */
    virtual ~LimitsNotSetException() throw() {}

    /**
     * Creating error message that the limits for the binning are not set
     * @return Error message specifying that the limits for the binning are not set
     */
    virtual const char *what() const throw() { return "Limits needed for the linear binning are not set."; }
  };

  /**
   * Constructor
   */
  TLinearBinning();

  /**
   * Constructor, defining the limits and the number of bins
   * @param[in] nbins Number of bins
   * @param[in] min Minimum bin edge of the binning
   * @param[in] max Maximum bin edge of the binning
   */
  TLinearBinning(Int_t nbins, Double_t min, Double_t max);
  virtual ~TLinearBinning() {}

  inline void Set(Int_t nbins, Double_t min, Double_t max);

  /**
   * Converting the linear binning in a set of bin edges
   * @param binedges
   */
  virtual void CreateBinEdges(TArrayD &binedges) const;

private:
  Int_t                                 fNbins;     ///< Number of bins
  Double_t                              fMinimum;   ///< Minimum of the binning
  Double_t                              fMaximum;   ///< Maximum of the binning
  Bool_t                                fLimitsSet; ///< Switch indicating that the binning is initialized
};

void TLinearBinning::Set(Int_t nbins, Double_t min, Double_t max){
  fMinimum = min;
  fMaximum = max;
  fNbins = nbins;
  fLimitsSet = true;
}

#endif /* TLINEARBINNING_H */
