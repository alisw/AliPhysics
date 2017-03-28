#ifndef ALIANALYSISMUMUPARAMETERS_H
#define ALIANALYSISMUMUPARAMETERS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/**
  @ingroup pwg_muondep_mumu

  @class AliAnalysisMuMuConfig

  @brief helper class to store steering options for other MuMu classes

 Holds some options (e.g. for the AliAnalysisMuMu and AliAnalysisMuMuEvolution classesi)
 like the list of triggers to consider, the fit to be performed, etc...
 both for real data and for simulations (which might differ in e.g.
 the naming of the triggers).
 This class reads an extern file config.mumu. Each line should be written as <key> : <value> <type>

 #DimuonTrigger: CMUL8-S-NOPF-MUON real sim
 MuonTrigger : CMSL7-8-NOPF-MUON sim
 MuonTrigger : CMSL7-7-NOPF-MUON sim

 \author Laurent Aphecetche, Subatech
 @author : Benjamin Audurier (Subatech)

 @todo : make it readeable/writeable from/to a simple ASCII file ?

 */

#include "TObject.h"
#include "TString.h"
#include <map>
#include <string>
#include <set>

class TObjArray;
class TMap;

class AliAnalysisMuMuConfig : public TObject
{

public:

  enum EColor
  {
    kBlue=1500,
    kOrange=1501,
    kGreen=1502
  };

  enum ETriggerType
  {
    kMB=1,
    kMUL=2,
    kMSL=3,
    kMSH=4,
    kMIX=5
  };

  enum EPerRunInfo
  {
    kMBTriggerClassName=0,
    kMULTriggerClassName,
    kMSLTriggerClassName,
    kMSHTriggerClassName,
    kCentralityName,
    kMIXTriggerClassName

  };

  AliAnalysisMuMuConfig();
  virtual ~AliAnalysisMuMuConfig();

  AliAnalysisMuMuConfig(const AliAnalysisMuMuConfig& other);
  AliAnalysisMuMuConfig& operator=(const AliAnalysisMuMuConfig& other);

  void Add(const char* key, const TString& line);

  void ReadFromFile(const char* inputfile);

  /** @name Per Run Information
   *  Methods to retrieve some run-dependent information
   */
  ///@{
  /** Switch board to the GetXXTriggerClassName methods below */
  TString GetTriggerClassName(ETriggerType tt, Int_t runNumber) const;
  /** Get the name of the Minimum Bias (MB) trigger class */
  TString GetMBTriggerClassName(Int_t runNumber) const;
  /** Get the name of the Mix trigger class */
  TString GetMIXTriggerClassName(Int_t runNumber) const;
  /** Get the name of the dimuon unlike sign (MUL) trigger class */
  TString GetMULTriggerClassName(Int_t runNumber) const;
  /** Get the name of the single muon low pt threshold (MSL) trigger class */
  TString GetMSLTriggerClassName(Int_t runNumber) const;
  /** Get the name of the single muon high pt threshold (MSH) trigger class */
  TString GetMSHTriggerClassName(Int_t runNumber) const;
  /** Get the name of the centrality selection used (e.g. V0A, PP) */
  TString GetCentralityName(Int_t runNumber) const;
  /** Do we have run information for all those runs ? */
  Bool_t HasRunInformation(std::set<int>& runs, Bool_t show=kFALSE) const;
  /** Set Run info */
  void SetRunInfo(const TString& runranges, const TString& runinfo);
  ///@}

  ///  Convenience method to retrieve the first element of lists
  TString First(const char* key, Bool_t simulation=kFALSE) const;

  /// Do we have value for the given type and data type ?
  Bool_t Has(const char* key, const char* value, Bool_t simulation=kFALSE) const;

  void Clear(Option_t* = "");

  static TString GetTriggerTypeName(ETriggerType tt);

  void SetColorScheme();

  void SetOCDBPath(const char* ocdbPath) { fOCDBPath = ocdbPath; }

  TString OCDBPath() const { return fOCDBPath; }

  void SetCompactGraphs(Bool_t value=kTRUE) { fIsCompactGraphs = value; }

  Bool_t CompactGraphs() { return fIsCompactGraphs; }

  void Print(Option_t* opt="") const;

  /** @name Configuration keys
   *  Valid keys to be used in the configuration file or in the query methods
   */
  ///@{
  /// Used to specify the list of specific trigger classnames used
  const char* RefMixTriggerKey() const;
  /// Used to specify the list of specific trigger classnames used
  const char* RefMixEventSelectionKey() const;
  /// Used to specify the list of dimuon trigger classnames used
  const char* DimuonTriggerKey() const;
  /// Used to specify the list of single muon trigger classnames used
  const char* MuonTriggerKey() const;
  /// Used to specify the list of minimum bias trigger classnames used
  const char* MinbiasTriggerKey() const;
  /// Used to specify the list of mix trigger classnames used
  const char* MixTriggerKey() const;
  /// Used to specify the list of event selections used
  const char* EventSelectionKey() const;
  /// Used to specify the list of event selections used for mixing
  const char* EventSelectionMixKey() const;
  /// Used to specify the list of pair selections used
  const char* PairSelectionKey() const;
  /// Used to specify the list of centralities used
  const char* CentralitySelectionKey() const;
  /// Used to specify the list of fits used
  const char* FitTypeKey() const;
  /// Used to specify a particular fit used
  const char* FitSingleKey() const;
  /// Used to specify whether or not to use compact graphs
  const char* CompactGraphKey() const;
  /// Used to specify OCDB path
  const char* OCDBPathKey() const;
  ///@}

  TObjArray* GetListElements(const char* type, Bool_t simulation) const;
  TObjArray* GetTriggersList(Bool_t simulation) const;

  void LoadAliceStyles();

private:

  enum EDataType { kSim = 1<<0, kReal = 1<<1 };

  void ShowList(const char* key, Bool_t isSimulation) const;

public:
  TMap* Map() const;

private:

  mutable TMap* fMap; ///< internal map of configuration values
  TString fOCDBPath; ///< OCDB to be used (raw:// by default)
  Bool_t fIsCompactGraphs; ///< whether the graph produced should be compact
  std::map<int,std::map<AliAnalysisMuMuConfig::EPerRunInfo,std::string> > fPerRunInfo; ///< per run information

  ClassDef(AliAnalysisMuMuConfig,4) // class to hold configuration of AliAnalysisMuMu(Evolution) class
};

#endif
