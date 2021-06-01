#ifndef ALIEMCALTRIGGERBINNINGCOMPONENT_H
#define ALIEMCALTRIGGERBINNINGCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TArrayD.h>
#include <TBinning.h>
#include <TNamed.h>
#include <iosfwd>

class TObjArray;

namespace PWGJE {
  
namespace EMCALJetTasks {

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
  /**
   * @class AliEMCalTriggerBinningData
   * @brief Wrapper for Binning data, connecting with name
   */
  class AliEMCalTriggerBinningData : public TNamed {
  public:

    /**
     * Default constructor
     */
    AliEMCalTriggerBinningData();

    /**
     * Named constructor, defining name and binning
     */
    AliEMCalTriggerBinningData(const char *name, TBinning *data);

    /**
     * Copy constructor
     * @param[in] data Reference for the copy
     */
    AliEMCalTriggerBinningData(const AliEMCalTriggerBinningData &data);

    /**
     * Assignment operator
     * @param[in] data Reference for assignment
     */
    AliEMCalTriggerBinningData &operator=(const AliEMCalTriggerBinningData &data);

    virtual ~AliEMCalTriggerBinningData();

    /**
     * Set the underlying binning. The data handler is the owner of the
     * binning, so any binning assigned will be deleted
     * @param [in] binning New binning to be set
     */
    void SetBinning(TBinning *binning);

    /**
     * Get the underlying binning
     * @return Underlying binning
     */
    TBinning *GetBinning() const { return fBinning; }

  private:
    TBinning                          *fBinning;        ///< Underlying binning data

    ClassDef(AliEMCalTriggerBinningData, 1);
  };

  /**
   * Main constructor
   */
  AliEMCalTriggerBinningComponent();

  /**
   * Copy constructor, creating a deep copy.
   * @param[in] ref Reference for the copy
   */
  AliEMCalTriggerBinningComponent(const AliEMCalTriggerBinningComponent &ref);

  /**
   * Assignment operator, doing a deep copy.
   * @param[in] ref Reference for the assignment
   */
  AliEMCalTriggerBinningComponent &operator=(const AliEMCalTriggerBinningComponent &ref);

  /**
   * Destructor
   */
  virtual ~AliEMCalTriggerBinningComponent();

  /**
   * Get binning information for a given axis. Return nullpointer if axis is not yet defined
   * @param[in] name axis name
   * @return the axis information
   */
  TBinning *GetBinning(const char *name) const;

  /**
   * Set binning for dimension. If not yet existing, create it
   * @param[in] dimname: axis name
   * @param[in] nbins: Number of bins
   * @param[in] binning: array of bin limits (size nbins+1)
   */
  void SetBinning(const char *dimname, int nbins, const double *binning);

  /**
   * Set binning for dimension. If not yet existing, create it.
   * @param[in] dimname axis name
   * @param[in] binning array of bin limits (size nbins+1)
   */
  void SetBinning(const char *dimname, const TArrayD &binning);

  /**
   * Set pre-defined binning initialized outside of the binning component
   * @param[in] dimname Name of the dimension
   * @param[in] binning Binning for the dimension
   */
  void SetBinning(const char *dimname, TBinning *binning);

  /**
   * Set a linear binning for dimension. If not yet existing, create it.
   * @param[in] dimname axis name
   * @param[in] nbins Number of bins
   * @param[in] min Minimum of the range (= lowest bin limit)
   * @param[in] max Maximum of the range (= highest bin limit)
   */
  void SetLinearBinning(const char *dirname, int nbins, double min, double max);

private:
  /**
   * Find binning for the given dimension in the binning component
   * @param[in] name of the dimension
   * @return Binning data (null if not found)
   */
  AliEMCalTriggerBinningData *FindBinning(const char *dim) const;

  TObjArray       *fDimensions;           ///< List of binnings (dimensions)

  ClassDef(AliEMCalTriggerBinningComponent, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERBINNINGCOMPONENT_H */
