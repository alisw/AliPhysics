/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <iostream>
#include <cstring>
#include <set>

#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TRandom2.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>
#include <TDatabasePDG.h>
#include <Riostream.h>
#include <TLine.h>

#include "AliDJetVReader.h"

/// \cond CLASSIMP
ClassImp(AliDJetVReader);
/// \endcond

/**
 * Default constructor.
 */
AliDJetVReader::AliDJetVReader():
  TObject(),
  fpTmin(0),
  fpTmax(0),
  fzmin(0),
  fzmax(0),
  fnDbins(0),
  fDbinpTedges(nullptr),
  fnJetbins(0),
  fJetbinpTedges(nullptr),
  fDEffValues(nullptr),
  fMassRebin(1),
  fMassPlot(nullptr),
  fMassVsJetPtPlot(nullptr)
{
}

/**
 * Copy constructor.
 * @param[in] source Const reference to an object to copy from
 */
AliDJetVReader::AliDJetVReader(const AliDJetVReader &source):
  TObject(),
  fpTmin(source.fpTmin),
  fpTmax(source.fpTmax),
  fzmin(source.fzmin),
  fzmax(source.fzmax),
  fnDbins(0),
  fDbinpTedges(nullptr),
  fnJetbins(0),
  fJetbinpTedges(nullptr),
  fDEffValues(nullptr),
  fMassRebin(source.fMassRebin),
  fMassPlot(nullptr),
  fMassVsJetPtPlot(nullptr)
{
  if (source.fnDbins > 0) {
    fnDbins = source.fnDbins;
    fDbinpTedges = new Double_t[fnDbins+1];
    memcpy(fDbinpTedges, source.fDbinpTedges, sizeof(Double_t)*(fnDbins+1));
    fDEffValues = new Double_t[fnDbins+1];
    memcpy(fDEffValues, source.fDEffValues, sizeof(Double_t)*(fnDbins+1));
  }
  if (source.fnJetbins > 0) {
    fnJetbins = source.fnJetbins;
    fJetbinpTedges = new Double_t[fnJetbins+1];
    memcpy(fJetbinpTedges, source.fJetbinpTedges, sizeof(Double_t)*(fnJetbins+1));
  }
}

/**
 * Destructor
 */
AliDJetVReader::~AliDJetVReader()
{
  if (fMassPlot) delete fMassPlot;
  if (fMassVsJetPtPlot) delete fMassVsJetPtPlot;
}

/**
 * Set the D meson pt bins
 * @param[in] nbins Number of pt bins
 * @param[in] ptedges Edges of the pt bins
 */
void AliDJetVReader::SetDmesonPtBins(Int_t nbins, Double_t* ptedges)
{
  fnDbins = nbins;
  if (fDbinpTedges) {
    delete[] fDbinpTedges;
    fDbinpTedges = nullptr;
  }
  if (nbins == 0) return;
  fDbinpTedges = new Double_t[fnDbins + 1];
  memcpy(fDbinpTedges, ptedges, sizeof(Double_t) * (fnDbins + 1));
}

/**
 * Set the jet pt bins
 * @param[in] nbins Number of pt bins
 * @param[in] ptedges Edges of the pt bins
 */
void AliDJetVReader::SetJetPtBins(Int_t nbins, Double_t* ptedges)
{
  fnJetbins = nbins;
  if (fJetbinpTedges) {
    delete[] fJetbinpTedges;
    fJetbinpTedges = nullptr;
  }
  if (nbins == 0) return;
  fJetbinpTedges = new Double_t[fnJetbins + 1];
  memcpy(fJetbinpTedges, ptedges, sizeof(Double_t) * (fnJetbins + 1));
}

/**
 * Set the efficiency values in bins of D meson pt
 * @param[in] sigmafix Values of the efficiency
 */
void AliDJetVReader::SetDmesonEfficiency(Double_t* effvalues)
{
  if (fDEffValues) {
    delete[] fDEffValues;
    fDEffValues = nullptr;
  }
  if (fnDbins == 0) return;
  fDEffValues = new Double_t[fnDbins];
  memcpy(fDEffValues, effvalues, sizeof(Double_t) * fnDbins);
}
