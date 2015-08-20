/**
 * \file AliEMCalTriggerBinningComponent.h
 * \brief Declaration of the classes AliEMCalTriggerBinningComponent and AliEMCalTriggerBinningDimension
 *
 * In this header file the classes AliEMCalTriggerBinningComponent and
 * AliEMCalTriggerBinningDimension are declared. Class AliEMCalTriggerBinningComponent is
 * the container for binnings in various dimensions shared among different analysis components
 * while AliEMCalTriggerBinningDimension contains the bin limits for a certain dimension.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIEMCALTRIGGERBINNINGCOMPONENT_H
#define ALIEMCALTRIGGERBINNINGCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TArrayD.h>
#include <TNamed.h>

class TObjArray;

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
namespace EMCalTriggerPtAnalysis {

/**
 * \class AliEMCalTriggerBinningDimension
 * \brief Binning definition for a certain dimension
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
   *
   * \param name Name of the dimension
   */
  AliEMCalTriggerBinningDimension(const char *name):
    TNamed(name, ""),
    fBinning()
  {}
  /**
   * Constructor initializing the dimension from a C-array
   *
   * \param name Name of the dimension
   * \param nbins Number of bins
   * \param binning Array of bin limits
   */
  AliEMCalTriggerBinningDimension(const char *name, int nbins, double *binning):
    TNamed(name, ""),
    fBinning(nbins+1, binning)
  {}
  /**
   * Constructor initializing the dimension from a ROOT array
   *
   * \param name Name of the dimension
   * \param binning Array of bin limits
   */
  AliEMCalTriggerBinningDimension(const char *name, const TArrayD &binning):
    TNamed(name, ""),
    fBinning(binning.GetSize(), binning.GetArray())
  {}
  /**
   * Destructor
   */
  ~AliEMCalTriggerBinningDimension() {}

  /**
   * Set the bin limits of the dimension from a C-array
   *
   * \param nbins Number of bins
   * \param binning Array of bin limits
   */
  void Set(int nbins, double *binning) { fBinning.Set(nbins+1, binning); }
  /**
   * Set the bin limits of the dimension from a ROOT array
   *
   * \param binning Array of bin limits
   */
  void Set(const TArrayD &binning) { fBinning = binning; }
  /**
   * Get array of bin limits
   *
   * \return C-array of bin limits
   */
  const double *GetBinLimits() const { return fBinning.GetArray(); }
  /**
   * Get the number of bins of the dimension
   *
   * \return Number of bins
   */
  int GetNumberOfBins() const { return fBinning.GetSize() - 1; }
  virtual void Print(Option_t *option="") const;

private:
  TArrayD fBinning;             ///< Bin limits

  /// \cond CLASSIMP
  ClassDef(AliEMCalTriggerBinningDimension, 1);
  /// \endcond
};

/**
 * \class AliEMCalTriggerBinningComponent
 * \brief Global binning definition for the high-\f$ p_{t} \f$ charged particle \f$ p_{t}\f$ analysis
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
  void SetBinning(const char *dimname, int nbins, double *binning);
  void SetBinning(const char *dimname, const TArrayD &binning);

private:
  TObjArray       *fDimensions;           ///< List of binnings (dimensions)

  /// \cond CLASSIMP
  ClassDef(AliEMCalTriggerBinningComponent, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERBINNINGCOMPONENT_H */
