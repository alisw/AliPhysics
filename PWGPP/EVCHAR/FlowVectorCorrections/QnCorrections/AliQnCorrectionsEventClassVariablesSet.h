#ifndef ALIQNCORRECTIONS_EVENTCLASSVARSET_H
#define ALIQNCORRECTIONS_EVENTCLASSVARSET_H
/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsEventClassVariablesSet.h
/// \brief Class that models the set of variables that define an event class for the Q vector correction framework

#include "AliQnCorrectionsEventClassVariable.h"

/// \class AliQnCorrectionsEventClassVariablesSet
/// \brief The set of variables which define an event class
///
/// Array of EventClassVariables that fully define the different
/// event classes considered within the Q vector correction framework.
/// The objects of this class are associated to concrete
/// detectors or detector configurations in such a way that
/// all Q vector corrections are performed according to the
/// event class the involved event is allocated.
///
/// Inherits all the methods of TObjArray specially the
/// subscript [] operator and Add method that allows
/// the array to expand.
///
/// The event class variables objects are not own by the array so,
/// they are not destroyed when the the set is destroyed. This allows
/// to create several sets with the same event class variables.
/// Pay attention to this when you create your event class variables,
/// they should live at least the same time you expect the sets to
/// live.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 4, 2016


class AliQnCorrectionsEventClassVariablesSet : public TObjArray {
public:
  /// Normal constructor
  /// \param n number of variables in the set
  AliQnCorrectionsEventClassVariablesSet(Int_t n = TCollection::kInitCapacity) : TObjArray(n) {}
  /// Copy constructor
  /// \param cecvs the object instance to be copied
  AliQnCorrectionsEventClassVariablesSet(const AliQnCorrectionsEventClassVariablesSet &cecvs) : TObjArray(cecvs) {}
  /// Default destructor
  virtual ~AliQnCorrectionsEventClassVariablesSet() {}

  /// Access the event class variable at the passed position
  /// \param i position in the array (starting at zero)
  /// \return the event class variable object a position i
  virtual AliQnCorrectionsEventClassVariable *At(Int_t i) const { return (AliQnCorrectionsEventClassVariable *) TObjArray::At(i); }

  void GetMultidimensionalConfiguration(Int_t *nbins, Double_t *minvals, Double_t *maxvals);

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsEventClassVariablesSet, 1);
/// \endcond
};

#endif /* ALIQNCORRECTIONS_EVENTCLASSVARSET_H */
