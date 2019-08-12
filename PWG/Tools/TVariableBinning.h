#ifndef TVARIABLEBINNING_H
#define TVARIABLEBINNING_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <exception>
#include <Rtypes.h>
#include <TArrayD.h>
#include <TBinning.h>

/**
 * @class TVariableBinning
 * @brief Class creating a variable binning, used in the histogram manager
 * @author Markus Fasel
 * @ingroup Histmanager
 * @since May 31st, 2016
 *
 * This class creates a varible (non-linear) binning. For this the user must provide
 * - A minimum
 * - A maximum
 * - The number of bins
 *
 * The information can be set either in the constructor
 *
 * ~~~{.cxx}
 * TArrayD binedges;
 * binedges[0] = 0.;
 * binedges[1] = 0.5;
 * binedges[2] = 1;
 * binedges[3] = 2.;
 * binedges[4] = 5.;
 * binedges[5] = 10.;
 * TVariableBinning mybinning(binedges);
 * ~~~
 *
 * or using the set function
 *
 * ~~~{.cxx}
 * TArrayD binedges;
 * binedges[0] = 0.;
 * binedges[1] = 0.5;
 * binedges[2] = 1;
 * binedges[3] = 2.;
 * binedges[4] = 5.;
 * binedges[5] = 10.;
 * TVariableBinning mybinning;
 * mybinning.Set(binedges);
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
class TVariableBinning : public TBinning{
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
  TVariableBinning();

  /**
   * Constructor, defining the bin edges from a c-array
   * @param[in] nbins Number of bins
   * @param[in] binedges Bin edges
   */
  TVariableBinning(Int_t nbins, const Double_t* binedges);

  /**
   * Constructor, defining the bin edges from a ROOT-array
   * @param[in] binedges Bin edges
   */
  TVariableBinning(const TArrayD &binedges);

  /**
   * Destructor
   */
  virtual ~TVariableBinning() {}

  /**
   * Implementation of the copy function for the linear binning class
   * @return Copy of this binning with the exact same bin edges
   */
  virtual TBinning *MakeCopy() const;

  /**
   * Set the binning from a c-array
   * @param[in] nbins Number of bins
   * @param[in] binedges Bin edges
   */
  inline void Set(Int_t nbins, const Double_t *binedges);

  /**
   * Set the binning from a ROOT array
   * @param[in] binedges Bin edges
   */
  void Set(const TArrayD &binedges) { fBinEdges = binedges; }


  /**
   * Converting the variable binning in a set of bin edges
   * @param binedges
   */
  virtual void CreateBinEdges(TArrayD &binedges) const;

private:
  TArrayD                         fBinEdges;

  ClassDef(TVariableBinning, 1);
};

void TVariableBinning::Set(Int_t nbins, const Double_t *binedges){
  fBinEdges.Set(nbins+1);
  for(int i = 0; i < nbins + 1; i++) fBinEdges[i] = binedges[i];
}

#endif /* TVARIABLEBINNING_H */
