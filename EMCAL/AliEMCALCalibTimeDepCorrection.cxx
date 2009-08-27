/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: $ */

// Objects of this class contain info on APD bias settings/voltages
//

#include <fstream>
#include <TString.h>
#include <TArrayF.h>

#include "AliEMCALCalibTimeDepCorrection.h"

using namespace std;

ClassImp(AliEMCALCalibTimeDepCorrection)

//____________________________________________________________________________
AliEMCALCalibTimeDepCorrection::AliEMCALCalibTimeDepCorrection() : 
  fNSuperModule(0),
  fSuperModuleData(0),
  fStartTime(0),
  fNTimeBins(0),
  fTimeBinSize(0)
{
  //Default constructor.
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::InitCorrection(Int_t nSM, Int_t nBins, Float_t val)
{
  // This methods assumes that you are using SuperModules 0..nSM-1
  fNSuperModule = nSM;
  if (fSuperModuleData) delete [] fSuperModuleData;
  fSuperModuleData = new AliEMCALSuperModuleCalibTimeDepCorrection[fNSuperModule];

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t nCorr = nBins;
  Float_t Correction = val;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibTimeDepCorrection &t = fSuperModuleData[i];
    t.fSuperModuleNum = iSM; // assume SMs are coming in order

    for (Int_t j=0; j<nAPDPerSM; j++) {

      // set size of TArray
      t.fCorrection[iCol][iRow].Set(nCorr);
      for (Int_t k=0; k<nCorr; k++) {
	// add to TArray
	t.fCorrection[iCol][iRow].AddAt(Correction, k); // AddAt = SetAt..
      }
    }

  } // i, SuperModule

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::SetCorrection(Int_t supModIndex, Int_t iCol, Int_t iRow, Int_t iBin, Float_t val)
{ // if you call for non-existing data, there may be a crash..
  fSuperModuleData[supModIndex].fCorrection[iCol][iRow].AddAt(val, iBin); // AddAt = SetAt..
 return; 
}

//____________________________________________________________________________
Float_t AliEMCALCalibTimeDepCorrection::GetCorrection(Int_t supModIndex, Int_t iCol, Int_t iRow, Int_t iBin)const
{ // if you call for non-existing data, there may be a crash..
  return fSuperModuleData[supModIndex].fCorrection[iCol][iRow].At(iBin);
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::ReadCalibTimeDepCorrectionInfo(Int_t nSM, const TString &txtFileName,
								    Bool_t swapSides)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALCalibTimeDepCorrection::ReadCalibTimeDepCorrectionInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;
  if (fSuperModuleData) delete [] fSuperModuleData;
  fSuperModuleData = new AliEMCALSuperModuleCalibTimeDepCorrection[fNSuperModule];

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t nCorr = 0;
  Float_t Correction = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibTimeDepCorrection &t = fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALCalibTimeDepCorrection::ReadCalibTimeDepCorrectionInfo - Error while reading input file; likely EOF..");
      return;
    }
    inputFile >> iSM;
    t.fSuperModuleNum = iSM;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow >> nCorr;

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      // set size of TArray
      t.fCorrection[iCol][iRow].Set(nCorr);
      for (Int_t k=0; k<nCorr; k++) {
	inputFile >> Correction;
	// add to TArray
	t.fCorrection[iCol][iRow].AddAt(Correction, k);
      }
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::WriteCalibTimeDepCorrectionInfo(const TString &txtFileName,
								     Bool_t swapSides)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALCalibTimeDepCorrection::WriteCalibTimeDepCorrectionInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t nCorr = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibTimeDepCorrection &t = fSuperModuleData[i];
    outputFile << t.fSuperModuleNum << endl;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      nCorr = t.fCorrection[iCol][iRow].GetSize();

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow << " " << nCorr << endl;
      for (Int_t k=0; k<nCorr; k++) {
	outputFile << t.fCorrection[iCol][iRow].At(k) << " ";
      }
      outputFile << endl;

    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
AliEMCALCalibTimeDepCorrection::~AliEMCALCalibTimeDepCorrection()
{
  delete [] fSuperModuleData;
}

//____________________________________________________________________________
AliEMCALCalibTimeDepCorrection::AliEMCALSuperModuleCalibTimeDepCorrection AliEMCALCalibTimeDepCorrection::GetSuperModuleCalibTimeDepCorrectionId(Int_t supModIndex)const
{
  AliEMCALSuperModuleCalibTimeDepCorrection t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fSuperModuleData)
    return t;

  return fSuperModuleData[supModIndex];
}

//____________________________________________________________________________
AliEMCALCalibTimeDepCorrection::AliEMCALSuperModuleCalibTimeDepCorrection AliEMCALCalibTimeDepCorrection::GetSuperModuleCalibTimeDepCorrectionNum(Int_t supModIndex)const
{
  AliEMCALSuperModuleCalibTimeDepCorrection t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fSuperModuleData)
    return t;

  for (int i=0; i<fNSuperModule; i++) {
    if (fSuperModuleData[i].fSuperModuleNum == supModIndex) {
      return fSuperModuleData[i];
    }
  }

  return t;
}

