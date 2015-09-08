#ifndef ALIHFCUTVARFDSUBANALYSISMANAGERDPLUS_H
#define ALIHFCUTVARFDSUBANALYSISMANAGERDPLUS_H
/// \class AliHFCutVarFDsubAnalysisManagerDplus
/// \brief Pt bin for the cut variation feed down method analysis
///
///
///
///
/// \author Fabrizio Grosa <fabrizio.grosa@to.infn.it>
/// \date Sep 1, 2015
#include "TH1F.h"

#include "AliHFCutVarFDsubAnalysisManager.h"

class AliHFCutVarFDsubAnalysisManagerDplus : public AliHFCutVarFDsubAnalysisManager {
protected:

  Double_t fNevents; // Event count for normalisation

  AliHFCutVarFDsubAnalysisManagerDplus(const AliHFCutVarFDsubAnalysisManagerDplus& analysisManagerDplus); /// Copy constructor

  TH1F* CalculateCrossSection(const TString AccFilePath,
                              const TString GenLimAccHistoName,
                              const TString GenAccHistoName,
                              const TString system,
                              const TString PromptOrFD); ///calculate cross section for pPb or pp system

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

  TH1F* GetCrossSecPrompt(const TString AccFilePath, const TString GenLimAccHistoName,const TString GenAccHistoName,const TString system) {
    return CalculateCrossSection(AccFilePath,GenLimAccHistoName,GenAccHistoName,system,"Prompt"); } /// get the prompt cross section

  TH1F* GetCrossSecFD(const TString AccFilePath,const TString GenLimAccHistoName,const TString GenAccHistoName,const TString system) {
    return CalculateCrossSection(AccFilePath,GenLimAccHistoName,GenAccHistoName,system,"FD"); } /// get the FD cross section

  /// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubAnalysisManagerDplus, 1);
  /// \endcond
};
#endif // ALIHFCUTVARFDSUBANALYSISMANAGERDPLUS_H
