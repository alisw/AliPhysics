#ifndef ALIHFCUTVARFDSUBANALYSISMANAGERDPLUS_H
#define ALIHFCUTVARFDSUBANALYSISMANAGERDPLUS_H
/// \class AliHFCutVarFDsubAnalysisManagerDplus
/// \brief Analysis manager specifically for the D+ for the cut variation feed down method analysis
///
///
///
///
/// \author Fabrizio Grosa <fabrizio.grosa@to.infn.it>
/// \date Sep 1, 2015

#include "AliHFCutVarFDsubAnalysisManager.h"

class AliHFCutVarFDsubAnalysisManagerDplus : public AliHFCutVarFDsubAnalysisManager {
protected:

  Double_t fNevents; // Event count for normalisation

  AliHFCutVarFDsubAnalysisManagerDplus(const AliHFCutVarFDsubAnalysisManagerDplus& analysisManagerDplus); /// Copy constructor
  AliHFCutVarFDsubAnalysisManagerDplus operator=(const AliHFCutVarFDsubAnalysisManagerDplus& analysisManagerDplus); // Assignment operator

public:
  AliHFCutVarFDsubAnalysisManagerDplus(); /// Default constructor
  ~AliHFCutVarFDsubAnalysisManagerDplus(); /// Destructor

  Int_t GetTHnSparses(const TString strFileData,
                      const TString strFileMC,
		      const TString strDirData,
		      const TString strDirMC,
		      const TString strListData,
		      const TString strListMC,
                      Bool_t MConly);

  void GetCuts(); /// get the list of cuts
  void GetAxes(Int_t version=0); /// get the axes of the THnSparses

  /// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubAnalysisManagerDplus, 1);
  /// \endcond
};
#endif // ALIHFCUTVARFDSUBANALYSISMANAGERDPLUS_H
