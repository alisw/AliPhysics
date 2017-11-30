#ifndef ALIEMCALAODFILTERBITCUT_H
#define ALIEMCALAODFILTERBITCUT_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliVCuts.h"

class AliAODTrack;

namespace PWG{

namespace EMCAL {

/**
 * @class AliEmcalAODFilterBitCuts
 * @brief Implementation of the AOD filter bit selection as virtual cut class
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since July 1st, 2016
 * @ingroup EMCALCOREFW
 *
 * AliEmcalAODFilterBitCuts implements the track filter bit selection of AOD tracks within
 * the framework of virtual cuts. This means that in filters it can be used as any virtual
 * cut class. Like this no special treatment is for AOD filter bit selection.
 *
 * The cut is applied in the function
 *
 * ~~~{.cxx}
 * AliEmcalAODFilterBitCuts::IsSelected(TObject *o);
 * ~~~
 *
 * which needs to be called in any filter with the AOD track as argument.
 *
 * The filter allows two selection modes:
 * + Any: at least one of the bits needs to be found
 * + All: all bits need to be found.
 */
class AliEmcalAODFilterBitCuts : public AliVCuts {
public:
  /**
   * @enum SelectionMode_t
   * @brief Definition of the mode how multiple filter bits are selected
   */
  enum SelectionMode_t {
    kSelAny = 0,      //!< Any (either of the bits set)
    kSelAll = 1       //!< All (all bits must be found in the AOD track)
  };

  /**
   * Dummy constructor
   */
  AliEmcalAODFilterBitCuts();

  /**
   * Constructor, defining also name and title
   * @param[in] name Name of the cut object
   * @param[in] title Title of the cut object
   */
  AliEmcalAODFilterBitCuts(const char *name, const char *title);

  virtual ~AliEmcalAODFilterBitCuts() {}

  /**
   * Request filter bit at a given bit number (not bit representation).
   * @param[in] bitnumber Number of the bit to be set (integer number)
   */
  void AddFilterBitNumber(ULong_t bitnumber) {if(bitnumber < sizeof(ULong_t)*8) fAODfilterBits |= 1 << bitnumber; }

  void AddStatusBitNumber(ULong_t bitnumber) {if(bitnumber < sizeof(ULong_t)*8) fAODstatusBits |= 1 << bitnumber; }

  /**
   * @brief Set the filter bits to be checked.
   *
   * Function using the bit representation, not the number of then bit.
   * Bits will be added to the existing bits (not if doReset is true).
   * @param[in] filterbits Filter bits requested
   * @param[in] doReset If true existing filter bits will be set to 0
   */
  void SetFilterBits(ULong_t filterbits, Bool_t doReset = false) { if(doReset) fAODfilterBits = 0; fAODfilterBits |= filterbits; }

  /**
   * @brief Set the track status bits to be checked.
   *
   * Function using the bit representation, not the number of then bit.
   * Bits will be added to the existing bits (not if doReset is true).
   * @param[in] filterbits Filter bits requested
   * @param[in] doReset If true existing filter bits will be set to 0
   */
  void SetStatusBits(ULong_t statusbits, Bool_t doReset = false) { if(doReset) fAODstatusBits = 0; fAODstatusBits |= statusbits; }

  /**
   * @brief Set the selection mode
   *
   * Can be any (any of the filter/status) bits set, or all. Note the
   * selection mode is applied to filter and status bits in the same way.
   * Filter and status bit cuts must be fulfilled independently.
   *
   * @param[in] mode Selection mode (any/all) used in the bit selection
   */
  void SetSelectionMode(SelectionMode_t mode) { fSelectionMode = mode; };

  /**
   * Select AOD tracks according which contain any of the bits. The way of the selection
   * is determined by the selection mode (default: any):
   * + Any: At least one bit needs to be found
   * + All: All bits must be found
   *
   * @param[in] o Object to be checked: Must be an AliAODTrack
   * @return True if object is selected. Always false in case object is not an AliAODTrack
   */
  virtual Bool_t IsSelected(TObject *o);

protected:

  /**
   * @brief Select track according to presence of track filter bits
   * @param[in] trk track to check
   * @return True if the track is selected, false otherwise
   */
  Bool_t IsFilterBitsSelected(const AliAODTrack *const trk) const;

  /**
   * @brief Select track according to presence of track status bits
   * @param[in] trk track to check
   * @return True if the track is selected, false otherwise
   */
  Bool_t IsStatusBitsSelected(const AliAODTrack *const trk) const;

  ULong_t                        fAODfilterBits;          ///< Requested filter bits
  ULong_t                        fAODstatusBits;          ///< Requested track status bits
  SelectionMode_t                fSelectionMode;          ///< Mode of the filter bit selection (any or all)

  /// \cond CLASSIMP
  ClassDef(AliEmcalAODFilterBitCuts, 1)
  /// \endcond
};

}

}

#endif /* ALIEMCALAODFILTERBITCUT_H */
