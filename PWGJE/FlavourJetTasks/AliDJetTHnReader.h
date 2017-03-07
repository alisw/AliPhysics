/**
 * \file AliDJetTHnReader.h
 * \brief Declaration of class AliDJetTHnReader
 *
 * In this header file the class AliDJetTHnReader is declared.
 * Class to extract the invariant mass plots used to extract the raw yield.
 * This implementation takes a THn histogram as input.
 * This is the output obtained from the task AliAnalysisTaskFlavourJetCorrelations
 *
 * \author Fabio Colamaria <fabio.colamaria@cern.ch>, INFN Bari
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Mar 7, 2017
 */

#ifndef ALIDJETTHNREADER_H
#define ALIDJETTHNREADER_H

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TString.h>

#include "AliDJetVReader.h"

class TFile;

/**
 * \class AliDJetTHnReader
 * \brief Implementation of an abstract class to read the invariant mass histograms used to extract the raw yield.
 *
 * Implementation of an abstract class to read the invariant mass histograms used to extract the raw yield.
 * This implementation takes a THn histogram as input.
 */
class AliDJetTHnReader : public AliDJetVReader {
public:

  AliDJetTHnReader();
  AliDJetTHnReader(const AliDJetTHnReader &source);
  virtual ~AliDJetTHnReader();

  void SetInputFilename(TString filename)   {fFileNameInput  = filename ; }
  void SetInputDirname(TString dirname)     {fDirName        = dirname  ; }
  void SetInputListname(TString listname)   {fListName       = listname ; }
  void SetInputObjectname(TString objname)  {fObjectName     = objname  ; }

  Bool_t ExtractInputMassPlotEffScale();
  Bool_t ExtractInputMassPlotSideband();

protected:

  TString     fFileNameInput ; ///< Name of input file
  TString     fDirName       ; ///< Name of input directory in the root file
  TString     fListName      ; ///< Name of input list
  TString     fObjectName    ; ///< Name of input container to extract the mass plot
  TFile      *fFileInput     ; //!<!File containing the task output

private:
  ClassDef(AliDJetTHnReader, 1);
};

#endif
