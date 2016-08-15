#ifndef ALIQNCORRECTIONS_CUTS_BASE_H
#define ALIQNCORRECTIONS_CUTS_BASE_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsCutsBase.h
/// \brief Base class for the classes that model the cuts support for the Q vector correction framework

#include <TObject.h>
#include <TObjArray.h>

/// \class AliQnCorrectionsCutsBase
/// \brief Base class for the Q vector correction cuts
///
/// Stores the external variable Id the cut should act on
///
/// Provides the interface for the set of different cuts
/// classes.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 22, 2016
class AliQnCorrectionsCutsBase: public TObject {

 public:
  AliQnCorrectionsCutsBase();
  AliQnCorrectionsCutsBase(const AliQnCorrectionsCutsBase &cut);
  AliQnCorrectionsCutsBase(Int_t varId);
  virtual ~AliQnCorrectionsCutsBase();

  /// Gets the variable Id the cut is applied to
  Int_t           GetVariableId() const { return fVarId; }

  /// Check if the actual variable value passes the cut
  ///
  /// Interface declaration function.
  /// Default behavior. Base class should not be instantiated.
  ///
  /// \param variableContainer the current variables content addressed by var Id
  /// \return kTRUE if the actual value passes the cut else kFALSE
  virtual Bool_t IsSelected(const Float_t *variableContainer) = 0;
 protected:
  Int_t         fVarId;   ///< The external Id for the variable in the data bank

  static const Int_t nHighestBitNumberSupported;       ///< The highest bit number the framework support

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsCutsBase, 1);
/// \endcond
};

#endif // ALIQNCORRECTIONS_CUTS_BASE_H
