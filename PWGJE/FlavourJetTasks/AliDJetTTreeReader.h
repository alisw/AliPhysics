/**
 * \file AliDJetTTreeReader.h
 * \brief Declaration of class AliDJetTTreeReader
 *
 * In this header file the class AliDJetTTreeReader is declared.
 * Class to extract the invariant mass plots used to extract the raw yield.
 * This implementation takes a TTree as input.
 * This is the output obtained from the task AliAnalysisTaskDmesonJets
 *
 * \author Fabio Colamaria <fabio.colamaria@cern.ch>, INFN Bari
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Mar 7, 2017
 */


#ifndef ALIDJETTTREEREADER_H
#define ALIDJETTTREEREADER_H

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <vector>
#include <string>
#include <TString.h>

#include "AliDJetVReader.h"

class TChain;

/**
 * \class AliDJetTTreeReader
 * \brief Implementation of an abstract class to read the invariant mass histograms used to extract the raw yield.
 *
 * Implementation of an abstract class to read the invariant mass histograms used to extract the raw yield.
 * This implementation takes a TTree histogram as input.
 */
class AliDJetTTreeReader : public AliDJetVReader {

public:
  AliDJetTTreeReader();
  AliDJetTTreeReader(const AliDJetTTreeReader &source);
  virtual ~AliDJetTTreeReader();

  void SetInputTreename(TString treename)           { fTreeName       = treename ; }
  void SetInputDBranchname(TString dname)           { fDBranchName    = dname    ; }
  void SetInputJetBranchname(TString jetname)       { fJetBranchName  = jetname  ; }
  void AddInputFileName(std::string filename)       { fInputFileNames.push_back(filename);}
  void SetMassEdgesAndBinWidthForMassPlot(Double_t mmin, Double_t mmax, Double_t mbinwidth) { fmassmin = mmin; fmassmax = mmax; fmasswidth = mbinwidth; }

  TChain* GenerateChain();

  Bool_t ExtractInputMassPlotEffScale();
  Bool_t ExtractInputMassPlotSideband();

protected:

  std::vector<std::string> fInputFileNames ; ///< Name of input file
  TString                  fTreeName       ; ///< Name of input TTree
  TString                  fDBranchName    ; ///< Name of input branch for D meson
  TString                  fJetBranchName  ; ///< Name of input branch for jet
  Double_t                 fmassmin        ; ///< Mass lower edge of inv.mass plots
  Double_t                 fmassmax        ; ///< Mass upper edge of inv.mass plots
  Double_t                 fmasswidth      ; ///< Mass plots bin width

private:
  ClassDef(AliDJetTTreeReader,1);
};

#endif
