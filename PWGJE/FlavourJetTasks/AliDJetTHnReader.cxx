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

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THnSparse.h>
#include <Riostream.h>

#include "AliDJetTHnReader.h"

/// \cond CLASSIMP
ClassImp(AliDJetTHnReader);
/// \endcond

/**
 * Default constructor.
 */
AliDJetTHnReader::AliDJetTHnReader():
  AliDJetVReader(),
  fFileNameInput(),
  fDirName(),
  fListName(),
  fObjectName(),
  fFileInput(nullptr)
{
}

/**
 * Copy constructor.
 * @param[in] source Const reference to an object to copy from
 */
AliDJetTHnReader::AliDJetTHnReader(const AliDJetTHnReader &source):
  AliDJetVReader(source),
  fFileNameInput(source.fFileNameInput),
  fDirName(source.fDirName),
  fListName(source.fListName),
  fObjectName(source.fObjectName),
  fFileInput(nullptr)
{
}

/**
 * Destructor
 */
AliDJetTHnReader::~AliDJetTHnReader()
{
}

/**
 * Extract the input mass plots for the efficiency scaled method.
 */
Bool_t AliDJetTHnReader::ExtractInputMassPlotEffScale()
{
  std::cout << "Extracting input mass plot: " << fpTmin << " to " << fpTmax << std::endl;

  TDirectoryFile* dir = dynamic_cast<TDirectoryFile*>(fFileInput->Get(fDirName.Data()));

  TList *histList = dynamic_cast<TList*>(dir->Get(Form("%s0", fListName.Data())));
  THnSparseF *sparse = dynamic_cast<THnSparseF*>(histList->FindObject(fObjectName.Data()));
  sparse->GetAxis(0)->SetRangeUser(fJetzBinEdges[0], fJetzBinEdges[fnJetzbins]);
  TH3D* hInvMassptD = static_cast<TH3D*>(sparse->Projection(3, 1, 2));
  hInvMassptD->SetName("hInvMassptD");

  if (!hInvMassptD) return kFALSE;

  TH1D* hmassjet[fnDbins];
  TH1D* hmassjet_scale[fnDbins];
  for (int i = 0; i < fnDbins; i++) {
    hmassjet[i] = nullptr;
    hmassjet_scale[i] = nullptr;
  }
  TH1D *hmass = nullptr;

  for (int j = 0; j < fnDbins; j++) {
    hmassjet[j] = hInvMassptD->ProjectionX(Form("hmassjet%d",j), hInvMassptD->GetYaxis()->FindBin(fpTmin), hInvMassptD->GetYaxis()->FindBin(fpTmax)-1, hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[j]), hInvMassptD->GetZaxis()->FindBin(fDbinpTedges[j + 1]) - 1);

    hmassjet_scale[j] = static_cast<TH1D*>(hmassjet[j]->Clone(Form("hmassjet%d_scale", j)));
    hmassjet_scale[j]->Scale(1. / fDEffValues[j]);

    if (!j) hmass = static_cast<TH1D*>(hmassjet_scale[j]->Clone("hmass"));
    else hmass->Add(hmassjet_scale[j]);
  }

  fMassPlot = static_cast<TH1D*>(hmass->Clone("inputSpectrum"));
  if (fMassRebin != 1) fMassPlot->Rebin(fMassRebin);

  return kTRUE;
}

/**
 * Extract the input mass plots for the side band method.
 */
Bool_t AliDJetTHnReader::ExtractInputMassPlotSideband()
{
  std::cout << "Extracting input mass plot: " << fpTmin << " to " << fpTmax << std::endl;

  double jetmin = 5, jetmax = 30;

  TDirectoryFile* dir = dynamic_cast<TDirectoryFile*>(fFileInput->Get(fDirName.Data()));

  TList *histList = dynamic_cast<TList*>(dir->Get(Form("%s0", fListName.Data())));
  THnSparseF *sparse = dynamic_cast<THnSparseF*>(histList->FindObject(fObjectName.Data()));
  sparse->GetAxis(0)->SetRangeUser(fJetzBinEdges[0], fJetzBinEdges[fnJetzbins]);
  sparse->GetAxis(1)->SetRangeUser(jetmin, jetmax);
  sparse->GetAxis(2)->SetRangeUser(fpTmin, fpTmax);

  fMassPlot = sparse->Projection(3);
  fMassPlot->SetName("inputSpectrum");
  if (fMassRebin != 1) fMassPlot->Rebin(fMassRebin);
  fMassVsJetPtPlot = sparse->Projection(3, 1);
  fMassVsJetPtPlot->SetName("hInvMassJetPt");

  return kTRUE;
}
