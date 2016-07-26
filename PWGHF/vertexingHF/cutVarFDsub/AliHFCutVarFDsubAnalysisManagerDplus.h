#ifndef ALIHFCUTVARFDSUBANALYSISMANAGERDPLUS_H
#define ALIHFCUTVARFDSUBANALYSISMANAGERDPLUS_H
/// \class AliHFCutVarFDsubAnalysisManagerDplus
/// \brief Pt bin for the cut variation feed down method analysis
///
///
///
///
/// \author Fabrizio Grosa <fabrizio.grosa@to.infn.it>, INFN Torino
/// \date Sep 1, 2015
#include "TH1F.h"
#include "TList.h"

#include "AliHFCutVarFDsubAnalysisManager.h"

class AliHFCutVarFDsubAnalysisManagerDplus : public AliHFCutVarFDsubAnalysisManager {
protected:

  Double_t fNevents; // Event count for normalisation
  Bool_t fPID; ///flag to activate the PID (it adds a cut on THnSparse)
  Int_t fPIDAxis; ///axis on THnSparses which corresponds to PID
  
  AliHFCutVarFDsubAnalysisManagerDplus(const AliHFCutVarFDsubAnalysisManagerDplus& analysisManagerDplus); /// Copy constructor
  AliHFCutVarFDsubAnalysisManagerDplus operator=(const AliHFCutVarFDsubAnalysisManagerDplus& analysisManagerDplus); // Assignment operator

  TH1F* CalculateCrossSection(const TString AccFilePath,
                              const TString GenLimAccHistoName,
                              const TString GenAccHistoName,
                              const TString system,
                              const TString PromptOrFD); ///calculate cross section for pPb or pp system

public:
  AliHFCutVarFDsubAnalysisManagerDplus(); /// Default constructor
  ~AliHFCutVarFDsubAnalysisManagerDplus(); /// Destructor
  
  Int_t GetTHnSparses(const TString strFileMC,
                      const TString strFileData,
                      const TString strDirMC,
                      const TString strDirData,
                      const TString strListMC,
                      const TString strListData,
                      Bool_t MConly);
  
  void GetCuts(Double_t*** cutslowset, ///first dimension: number of set, second: number of pt bin, third: the number of cut variable (pt included)
               Double_t*** cutshighset,
               Double_t** means, ///first dimension: number of set, second: number of pt bin
               Double_t** sigmas, 
               Int_t Rebin,
               Int_t fsig,
               Int_t fbkg,
	       Double_t min,
	       Double_t max,
	       Int_t nSets,
               Int_t nPtBins,
               Int_t nCutVariables); /// get the list of cuts
  
  void GetAxes(UInt_t* dataAxesNo,UInt_t* MCGenAxesNo,UInt_t* MCCutAxesNo,TString* axesName,Int_t nAxes); /// get the axes of the THnSparses

  TH1F* GetCrossSecPrompt(const TString AccFilePath, const TString GenLimAccHistoName,const TString GenAccHistoName,const TString system) {
    return CalculateCrossSection(AccFilePath,GenLimAccHistoName,GenAccHistoName,system,"Prompt"); } /// get the prompt cross section

  TH1F* GetCrossSecFD(const TString AccFilePath,const TString GenLimAccHistoName,const TString GenAccHistoName,const TString system) {
    return CalculateCrossSection(AccFilePath,GenLimAccHistoName,GenAccHistoName,system,"FD"); } /// get the FD cross section
  
  TH1F* GetYieldsPrompt() const {return fCorrYieldPrompt;}
  TH1F* GetYieldsFD() const {return fCorrYieldFD;}
  TList* GetPromptFraction() const {return fFprompt;}
  TList* GetPromptFractionRaw() const {return fFpromptRaw;}
  TList* GetResiduals() const {return fResiduals;}
  TList* GetPulls() const {return fPulls;}

  void SetPID(Bool_t isPIDon=kTRUE, Int_t PIDaxis=7) {fPID = isPIDon; fPIDAxis = PIDaxis;}
  
  /// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubAnalysisManagerDplus, 3);
  /// \endcond
};
#endif // ALIHFCUTVARFDSUBANALYSISMANAGERDPLUS_H
