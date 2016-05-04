#ifndef ALIEMCALTRIGGERBINNINGCOMPONENT_H
#define ALIEMCALTRIGGERBINNINGCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TArrayD.h>
#include <TNamed.h>
#include <Riosfwd.h>

class TObjArray;

namespace EMCalTriggerPtAnalysis {

/**
 * @class AliEMCalTriggerBinningDimension
 * @brief Binning definition for a certain dimension
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @date Dec 12, 2014
 *
 * This class contains the binning definition for a certain dimension. By construction
 * a variable binning is assumed.
 */
class AliEMCalTriggerBinningDimension : public TNamed{
public:

  /**
   * Dummy Constructor
   */
  AliEMCalTriggerBinningDimension():
    TNamed(),
    fBinning()
  {}

  /**
   * Named constructor
   * @param[in] name Name of the dimension
   */
  AliEMCalTriggerBinningDimension(const char *name):
    TNamed(name, ""),
    fBinning()
  {}

  /**
   * Constructor initializing the dimension from a C-array
   * @param[in] name Name of the dimension
   * @param[in] nbins Number of bins
   * @param[in] binning Array of bin limits
   */
  AliEMCalTriggerBinningDimension(const char *name, int nbins, const double *binning):
    TNamed(name, ""),
    fBinning(nbins+1, binning)
  {}

  /**
   * Constructor initializing the dimension from a ROOT array
   * @param[in] name Name of the dimension
   * @param[in] binning Array of bin limits
   */
  AliEMCalTriggerBinningDimension(const char *name, const TArrayD &binning):
    TNamed(name, ""),
    fBinning(binning.GetSize(), binning.GetArray())
  {}

  /**
   * Destructor
   */
  virtual ~AliEMCalTriggerBinningDimension() {}

  /**
   * Set the bin limits of the dimension from a C-array
   * @param[in] nbins Number of bins
   * @param[in] binning Array of bin limits
   */
  void Set(int nbins, const double *binning) { fBinning.Set(nbins+1, binning); }

  /**
   * Set the bin limits of the dimension from a ROOT array
   * @param[in] binning Array of bin limits
   */
  void Set(const TArrayD &binning) { fBinning = binning; }

  /**
   * Get array of bin limits
   * @return C-array of bin limits
   */
  const double *GetBinLimits() const { return fBinning.GetArray(); }

  /**
   * Get the array of bin limits for this dimension
   * @return Array of bin limits
   */
  const TArrayD &GetBinning() const { return fBinning; }

  /**
   * Initialize output array with binning stored in this dimension
   * @param[out] out Array to initialize with this binning
   */
  void InitializeArray(TArrayD out) const { out = fBinning; }

  /**
   * Get the number of bins of the dimension
   * @return Number of bins
   */
  int GetNumberOfBins() const { return fBinning.GetSize() - 1; }
  virtual void Print(Option_t *option="") const;
  void PrintStream(std::ostream &stream) const;

  friend std::ostream &operator<<(std::ostream &stream, const AliEMCalTriggerBinningDimension &dim);

private:
  TArrayD fBinning;             ///< Bin limits

  /// \cond CLASSIMP
  ClassDef(AliEMCalTriggerBinningDimension, 1);
  /// \endcond
};

/**
 * @class AliEMCalTriggerBinningComponent
 * @brief Global binning definition for the high-\f$ p_{t} \f$ charged particle \f$ p_{t}\f$ analysis
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @date Dec 12, 2014
 *
 * This class contains the binning definition for various dimensions shared globally among
 * analysis components. The dimensions are handled via the class AliEMCalTriggerBinningDimension.
 * Getters and setters are provided.
 */
class AliEMCalTriggerBinningComponent: public TObject {
public:
  AliEMCalTriggerBinningComponent();
  AliEMCalTriggerBinningComponent(const AliEMCalTriggerBinningComponent &ref);
  AliEMCalTriggerBinningComponent &operator=(const AliEMCalTriggerBinningComponent &ref);
  virtual ~AliEMCalTriggerBinningComponent();

  AliEMCalTriggerBinningDimension *GetBinning(const char *name) const;
  void SetBinning(const char *dimname, int nbins, const double *binning);
  void SetBinning(const char *dimname, const TArrayD &binning);
  void SetLinearBinning(const char *dirname, int nbins, double min, double max);

private:
  TObjArray       *fDimensions;           ///< List of binnings (dimensions)

  /// \cond CLASSIMP
  ClassDef(AliEMCalTriggerBinningComponent, 1);
  /// \endcond
};

std::ostream &operator<<(std::ostream& stream, const AliEMCalTriggerBinningDimension &dim);

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERBINNINGCOMPONENT_H */
