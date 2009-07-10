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

#include "AliEMCALBiasAPD.h"

using namespace std;

ClassImp(AliEMCALBiasAPD)

//____________________________________________________________________________
AliEMCALBiasAPD::AliEMCALBiasAPD() : 
  fNSuperModule(0),
  fSuperModuleData(0)
{
  //Default constructor.
}

//____________________________________________________________________________
void AliEMCALBiasAPD::ReadBiasAPDInfo(Int_t nSM, const TString &txtFileName,
				      Bool_t swapSides)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALBiasAPD::ReadBiasAPDInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;
  if (fSuperModuleData) delete [] fSuperModuleData;
  fSuperModuleData = new AliEMCALSuperModuleBiasAPD[fNSuperModule];

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t iElecId = 0;
  Int_t iDAC = 0;
  Float_t voltage = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleBiasAPD &t = fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALBiasAPD::ReadBiasAPDInfo - Error while reading input file; likely EOF..");
      return;
    }
    inputFile >> iSM;
    t.fSuperModuleNum = iSM;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow >> iElecId >> iDAC >> voltage;

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      t.fElecId[iCol][iRow] = iElecId;
      t.fDAC[iCol][iRow] = iDAC;
      t.fVoltage[iCol][iRow] = voltage;
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALBiasAPD::WriteBiasAPDInfo(const TString &txtFileName,
				       Bool_t swapSides)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALBiasAPD::WriteBiasAPDInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t iElecId = 0;
  Int_t iDAC = 0;
  Float_t voltage = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleBiasAPD &t = fSuperModuleData[i];
    outputFile << t.fSuperModuleNum << endl;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      iElecId = t.fElecId[iCol][iRow];
      iDAC = t.fDAC[iCol][iRow];
      voltage = t.fVoltage[iCol][iRow];

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow << " " 
		 << iElecId << " " << iDAC << " "
		 << voltage << endl;
    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
AliEMCALBiasAPD::~AliEMCALBiasAPD()
{
  delete [] fSuperModuleData;
}

//____________________________________________________________________________
AliEMCALBiasAPD::AliEMCALSuperModuleBiasAPD AliEMCALBiasAPD::GetSuperModuleBiasAPDId(Int_t supModIndex)const
{
  AliEMCALSuperModuleBiasAPD t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fSuperModuleData)
    return t;

  return fSuperModuleData[supModIndex];
}

//____________________________________________________________________________
AliEMCALBiasAPD::AliEMCALSuperModuleBiasAPD AliEMCALBiasAPD::GetSuperModuleBiasAPDNum(Int_t supModIndex)const
{
  AliEMCALSuperModuleBiasAPD t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fSuperModuleData)
    return t;

  for (int i=0; i<fNSuperModule; i++) {
    if (fSuperModuleData[i].fSuperModuleNum == supModIndex) {
      return fSuperModuleData[i];
    }
  }

  return t;
}

