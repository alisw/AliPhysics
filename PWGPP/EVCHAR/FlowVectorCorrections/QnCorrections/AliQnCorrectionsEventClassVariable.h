#ifndef ALIQNCORRECTIONS_EVENTCLASSVAR_H
#define ALIQNCORRECTIONS_EVENTCLASSVAR_H
/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsEventClassVariable.h
/// \brief Class that models variables used for defining an event class within the Q vector correction framework

/// \class AliQnCorrectionsEventClassVariable
/// \brief One variable used for defining an event class
///
/// Class defining one variable and its associated binning allowing
/// its use for the definition of event classes within the Q vector
/// correction framework.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 4, 2016
 

#include <TObject.h>
#include <TObjArray.h>


class AliQnCorrectionsEventClassVariable : public TObject {

 public:
  AliQnCorrectionsEventClassVariable();
  AliQnCorrectionsEventClassVariable(const AliQnCorrectionsEventClassVariable &ecv);
  AliQnCorrectionsEventClassVariable(Int_t varId, const char *varname, Int_t nbins, Double_t min, Double_t max);
  AliQnCorrectionsEventClassVariable(Int_t varId, const char *varname, Int_t nbins, Double_t *bins);
  AliQnCorrectionsEventClassVariable(Int_t varId, const char *varname, Double_t binsArray[][2]);
  ~AliQnCorrectionsEventClassVariable();

  /// Gets the variable unique Id
  Int_t           GetVariableId() const { return fVarId; }
  /// Gets the variable name / label
  const char *    GetVariableLabel() const { return (const char *) fLabel; }
  /// Gets the number of bins
  Int_t           GetNBins() const { return fNBins; }
  /// Gets the actual bins edges array
  const Double_t *GetBins() const { return fBins; }
  /// Gets the lower edge for the passed bin number
  /// \param bin bin number starting from one
  Double_t        GetBinLowerEdge(Int_t bin) const { return (((bin < 1) || (bin > fNBins)) ? 0.0 : fBins[bin-1]); }
  /// Gets the upper edge for the passed bin number
  /// \param bin bin number starting from one
  Double_t        GetBinUpperEdge(Int_t bin) const { return (((bin < 1) || (bin > fNBins)) ? 0.0 : fBins[bin]); }

  /// Gets the lowest variable value considered
  Double_t        GetLowerEdge() {return fBins[0]; }
  /// Gets the highest variabel value considered
  Double_t        GetUpperEdge() {return fBins[fNBins]; }

 private:
  Int_t         fVarId;        ///< The external Id for the variable in the data bank
  Int_t         fNBins;        ///< The number of bins for the variable when shown in a histogram
  Int_t         fNBinsPlusOne; ///< the number of bins plus one. Needed for object persistence
  /// Bin edges array for the variable when shown in a histogram
  Double_t     *fBins;         //[fNBinsPlusOne]
  TString       fLabel;        ///< Label to use in an axis that shows the variable

 private:
  /// Assignment operator
  /// Not allowed. Forced private.
  AliQnCorrectionsEventClassVariable& operator= (const AliQnCorrectionsEventClassVariable &);
/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsEventClassVariable, 1);
/// \endcond
};

#endif /* ALIQNCORRECTIONS_EVENTCLASSVAR_H */
