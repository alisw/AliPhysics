#ifndef ALIHFCUTVARFDSUBANALYSISMANAGERD0_H
#define ALIHFCUTVARFDSUBANALYSISMANAGERD0_H
/// \class AliHFCutVarFDsubAnalysisManagerD0
/// \brief Analysis manager with D0 specific functionality for the cut variation feed down method analysis
///
///
///
///
/// \author Felix Reidt <felix.reidt@cern.ch>, CERN
/// \date Aug 17, 2015

#include "AliHFCutVarFDsubAnalysisManager.h"

class AliHFCutVarFDsubAnalysisManagerD0 : public AliHFCutVarFDsubAnalysisManager {
protected:

  Double_t fNevents; // Event count for normalisation

  AliHFCutVarFDsubAnalysisManagerD0(const AliHFCutVarFDsubAnalysisManagerD0& am); /// Copy constructor

public:
  AliHFCutVarFDsubAnalysisManagerD0(); /// Default constructor
  ~AliHFCutVarFDsubAnalysisManagerD0(); /// Destructor

  Int_t GetTHnSparses(const TString strFileData,
                      const TString strFilePrompt,
                      const TString strContNamePrompt,
                      const TString strFileSec,
                      const TString strContNameSec,
                      Bool_t MConly);

  void GetCuts();
  void GetAxes(Int_t version=0); /// get the axes of the THnSparses

  /// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubAnalysisManagerD0, 1);
  /// \endcond
};
#endif // ALIHFCUTVARFDSUBANALYSISMANAGERD0_H
