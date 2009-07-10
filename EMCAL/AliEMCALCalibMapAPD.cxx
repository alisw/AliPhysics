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

// Objects of this class contain info on APD calibration and map info
//

#include <fstream>
#include <TString.h>

#include "AliEMCALCalibMapAPD.h"

using namespace std;

ClassImp(AliEMCALCalibMapAPD)

//____________________________________________________________________________
AliEMCALCalibMapAPD::AliEMCALCalibMapAPD() : 
  fNSuperModule(0),
  fSuperModuleData(0)
{
  //Default constructor.
}

//____________________________________________________________________________
void AliEMCALCalibMapAPD::ReadCalibMapAPDInfo(Int_t nSM, const TString &txtFileName,
				      Bool_t swapSides)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALCalibMapAPD::ReadCalibMapAPDInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;
  if (fSuperModuleData) delete [] fSuperModuleData;
  fSuperModuleData = new AliEMCALSuperModuleCalibMapAPD[fNSuperModule];

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  // list of values to be read
  Int_t iHW = 0;
  Int_t APDNum = 0;
  Float_t V30 = 0;     
  Float_t Par[3] = {0};   
  Float_t ParErr[3] = {0}; 
  Int_t BreakDown = 0;
  Float_t DarkCurrent = 0; 
  // end - all values

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibMapAPD &t = fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALCalibMapAPD::ReadCalibMapAPDInfo - Error while reading input file; likely EOF..");
      return;
    }
    inputFile >> iSM;
    t.fSuperModuleNum = iSM;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow >> iHW 
		>> APDNum >> V30 
		>> Par[0] >> Par[1] >> Par[2]
		>> ParErr[0] >> ParErr[1] >> ParErr[2]
		>> BreakDown >> DarkCurrent;

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      AliEMCALCalibMapAPDVal &v = t.fAPDVal[iCol][iRow];

      v.fHardWareId = iHW;
      v.fAPDNum = APDNum;
      v.fV30 = V30;
      v.fPar[0] = Par[0];
      v.fPar[1] = Par[1];
      v.fPar[2] = Par[2];
      v.fParErr[0] = ParErr[0];
      v.fParErr[1] = ParErr[1];
      v.fParErr[2] = ParErr[2];
      v.fBreakDown = BreakDown;
      v.fDarkCurrent = DarkCurrent;
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibMapAPD::WriteCalibMapAPDInfo(const TString &txtFileName,
				       Bool_t swapSides)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALCalibMapAPD::WriteCalibMapAPDInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibMapAPD &t = fSuperModuleData[i];
    outputFile << t.fSuperModuleNum << endl;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      AliEMCALCalibMapAPDVal &v = t.fAPDVal[iCol][iRow];

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow << " " << v.fHardWareId 
		 << " " << v.fAPDNum << " " << v.fV30 
		 << " " << v.fPar[0] << " " << v.fPar[1] << " " << v.fPar[2]
		 << " " << v.fParErr[0] << " " << v.fParErr[1] << " " << v.fParErr[2]
		 << " " << v.fBreakDown << " " << v.fDarkCurrent << endl;
    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
AliEMCALCalibMapAPD::~AliEMCALCalibMapAPD()
{
  delete [] fSuperModuleData;
}

//____________________________________________________________________________
AliEMCALCalibMapAPD::AliEMCALSuperModuleCalibMapAPD AliEMCALCalibMapAPD::GetSuperModuleCalibMapAPDId(Int_t supModIndex)const
{
  AliEMCALSuperModuleCalibMapAPD t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fSuperModuleData)
    return t;

  return fSuperModuleData[supModIndex];
}

//____________________________________________________________________________
AliEMCALCalibMapAPD::AliEMCALSuperModuleCalibMapAPD AliEMCALCalibMapAPD::GetSuperModuleCalibMapAPDNum(Int_t supModIndex)const
{
  AliEMCALSuperModuleCalibMapAPD t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fSuperModuleData)
    return t;

  for (int i=0; i<fNSuperModule; i++) {
    if (fSuperModuleData[i].fSuperModuleNum == supModIndex) {
      return fSuperModuleData[i];
    }
  }

  return t;
}

