#ifndef ALIQNCORRECTIONS_CUTSSET_H
#define ALIQNCORRECTIONS_CUTSSET_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

#include "AliQnCorrectionsCutsBase.h"

/// \file AliQnCorrectionsCutsSet.h
/// \brief Set of cuts class support for the Q vector correction framework

/// \class AliQnCorrectionsCutsSet
/// \brief Set of cuts to assign to a detector
///
/// Array of cuts that contains the whole set of cuts to assign to
/// a concrete detector or detector configuration to filter the desired
/// events (in the broad sense) it will have to handle.
///
/// Inherits all the methods of TObjArray specially the
/// subscript [] operator and Add method that allows
/// the array to expand.
///
/// Provides IsSelected that goes through the whole set of cuts
/// to check whether the current variables values pass the them.
///
/// The cuts objects are not own by the array so,
/// they are not destroyed when the the set is destroyed. This allows
/// to create several sets with the same cuts.
/// Pay attention to this when you create your cuts,
/// they should live at least the same time you expect the sets to
/// live.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 4, 2016


class AliQnCorrectionsCutsSet : public TObjArray {
public:
  /// Normal constructor
  /// \param n number of cuts in the set
  AliQnCorrectionsCutsSet(Int_t n = TCollection::kInitCapacity) : TObjArray(n) {}
  /// Copy constructor
  /// \param ccs the object instance to be copied
  AliQnCorrectionsCutsSet(const AliQnCorrectionsCutsSet &ccs) : TObjArray(ccs) {}
  /// Default destructor
  virtual ~AliQnCorrectionsCutsSet() {}

  /// Access the event class variable at the passed position
  /// \param i position in the array (starting at zero)
  /// \return the event class variable object a position i
  virtual AliQnCorrectionsCutsBase *At(Int_t i) const { return (AliQnCorrectionsCutsBase *) TObjArray::At(i); }

  Bool_t IsSelected(const Float_t *variableContainer);

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsCutsSet, 1);
/// \endcond
};

/// Checks that the current content of the variableContainer passes
/// the whole set of cuts by going through all the array components
///
/// \param variableContainer the current variables content addressed by var Id
/// \return kTRUE if the actual values pass the set of cuts else kFALSE
inline Bool_t AliQnCorrectionsCutsSet::IsSelected(const Float_t *variableContainer) {
  for (Int_t icut = 0; icut < GetEntriesFast(); icut++) {
    if (!At(icut)->IsSelected(variableContainer)) {
      return kFALSE;
    }
  }
  return kTRUE;
}

#endif // ALIQNCORRECTIONS_CUTSSET_H
