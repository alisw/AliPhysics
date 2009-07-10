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

// Objects of this class contain basis for absolute calibrations
//

#include <fstream>
#include <TString.h>

#include "AliEMCALCalibAbs.h"

using namespace std;

ClassImp(AliEMCALCalibAbs)

//____________________________________________________________________________
AliEMCALCalibAbs::AliEMCALCalibAbs() : 
  fNSuperModule(0),
  fSuperModuleData(0)
{
  //Default constructor.
}

//____________________________________________________________________________
void AliEMCALCalibAbs::ReadCalibAbsInfo(Int_t nSM, const TString &txtFileName,
				      Bool_t swapSides)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALCalibAbs::ReadCalibAbsInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;
  if (fSuperModuleData) delete [] fSuperModuleData;
  fSuperModuleData = new AliEMCALSuperModuleCalibAbs[fNSuperModule];

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t id = 0;

  // list of values to be read
  // first: overall values for the whole SuperModule
  Int_t CalibMethod; 
  Int_t CalibPass; 
  Int_t CalibTime; 
  Float_t AbsoluteGain; 
  // second: additional info for LED Reference and SM temperature
  Float_t LEDRefAmp;
  Float_t LEDRefAmpRMS;
  Float_t LEDRefHighLowRatio;
  Int_t LEDRefHighLow;
  Float_t Temperature;
  Float_t TemperatureRMS;
  // third: info for each tower
  Float_t RelativeGain; // (ADC>GeV relative gain/conversion), value around 1
  Float_t HighLowRatio; // value around 16 or so
  Int_t HighLow; // 
  Float_t LEDAmp; // low gain eq. amplitude
  Float_t LEDAmpRMS; //
  // end - all values

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibAbs &t = fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALCalibAbs::ReadCalibAbsInfo - Error while reading input file; likely EOF..");
      return;
    }
    inputFile >> iSM;
    t.fSuperModuleNum = iSM;

    // first: overall values for the whole SuperModule
    inputFile >> CalibMethod >> CalibPass >> CalibTime >> AbsoluteGain;
    t.fCalibMethod = CalibMethod;
    t.fCalibPass = CalibPass;
    t.fCalibTime = CalibTime;
    t.fAbsoluteGain = AbsoluteGain;

    // second: additional info for LED Reference and SM temperature
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      inputFile >> id >> LEDRefAmp >> LEDRefAmpRMS >> LEDRefHighLowRatio >> LEDRefHighLow;
      t.fLEDRefAmp[id] = LEDRefAmp;
      t.fLEDRefAmpRMS[id] = LEDRefAmpRMS;
      t.fLEDRefHighLowRatio[id] = LEDRefHighLowRatio;
      t.fLEDRefHighLow[id] = LEDRefHighLow;
    }
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      inputFile >> id >> Temperature >> TemperatureRMS;
      t.fTemperature[id] = Temperature;
      t.fTemperatureRMS[id] = TemperatureRMS;
    }

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow 
		>> RelativeGain >> HighLowRatio >> HighLow >> LEDAmp >> LEDAmpRMS;

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      AliEMCALCalibAbsVal &v = t.fAPDVal[iCol][iRow];

      v.fRelativeGain = RelativeGain;
      v.fHighLowRatio = HighLowRatio;
      v.fHighLow = HighLow;
      v.fLEDAmp = LEDAmp;
      v.fLEDAmpRMS = LEDAmpRMS;
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibAbs::WriteCalibAbsInfo(const TString &txtFileName,
				       Bool_t swapSides)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALCalibAbs::WriteCalibAbsInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibAbs &t = fSuperModuleData[i];
    // first: overall values for the whole SuperModule
    outputFile << t.fSuperModuleNum << endl;
    outputFile << t.fCalibMethod << " " 
	       << t.fCalibPass << " " 
	       << t.fCalibTime << " " 
	       << t.fAbsoluteGain << endl;

    // second: additional info for LED Reference and SM temperature
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALLEDRefs; j++) {
      outputFile << j << " " << t.fLEDRefAmp[j] << " " << t.fLEDRefAmpRMS[j] 
		 << " " << t.fLEDRefHighLowRatio[j] << " " << t.fLEDRefHighLow[j] 
		 << endl;
    }
    for (Int_t j=0; j<AliEMCALGeoParams::fgkEMCALTempSensors; j++) {
      outputFile << j << " " << t.fTemperature[j] << " " << t.fTemperatureRMS[j] << endl;
    }

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      AliEMCALCalibAbsVal &v = t.fAPDVal[iCol][iRow];

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow 
		 << " " << v.fRelativeGain
		 << " " << v.fHighLowRatio 
		 << " " << v.fHighLow 
		 << " " << v.fLEDAmp 
		 << " " << v.fLEDAmpRMS << endl;
    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
AliEMCALCalibAbs::~AliEMCALCalibAbs()
{
  delete [] fSuperModuleData;
}

//____________________________________________________________________________
AliEMCALCalibAbs::AliEMCALSuperModuleCalibAbs AliEMCALCalibAbs::GetSuperModuleCalibAbsId(Int_t supModIndex)const
{
  AliEMCALSuperModuleCalibAbs t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fSuperModuleData)
    return t;

  return fSuperModuleData[supModIndex];
}

//____________________________________________________________________________
AliEMCALCalibAbs::AliEMCALSuperModuleCalibAbs AliEMCALCalibAbs::GetSuperModuleCalibAbsNum(Int_t supModIndex)const
{
  AliEMCALSuperModuleCalibAbs t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fSuperModuleData)
    return t;

  for (int i=0; i<fNSuperModule; i++) {
    if (fSuperModuleData[i].fSuperModuleNum == supModIndex) {
      return fSuperModuleData[i];
    }
  }

  return t;
}

