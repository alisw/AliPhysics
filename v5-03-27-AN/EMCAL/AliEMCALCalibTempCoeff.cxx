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

// Objects of this class contain temperature-dependence coefficients
//

#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "AliEMCALCalibTempCoeff.h"

using namespace std;

ClassImp(AliEMCALCalibTempCoeff)

//____________________________________________________________________________
AliEMCALCalibTempCoeff::AliEMCALCalibTempCoeff(const int nSM) : 
  fNSuperModule(nSM),
  fSuperModuleData()
{
  //Default constructor.
  for (int i=0; i<fNSuperModule; i++) {
    fSuperModuleData.Add(new AliEMCALSuperModuleCalibTempCoeff(i));
  }
  fSuperModuleData.Compress(); // compress the TObjArray
  fSuperModuleData.SetOwner(kTRUE); 
}

//____________________________________________________________________________
void AliEMCALCalibTempCoeff::ReadTextCalibTempCoeffInfo(Int_t nSM, const TString &txtFileName,
					    Bool_t swapSides)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALCalibTempCoeff::ReadCalibTempCoeffInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;

  // list of values to be read
  Int_t iSrc = 0; 
  Float_t tempCoeff = 0; 
  // end - all values

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibTempCoeff * t = (AliEMCALSuperModuleCalibTempCoeff*) fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALCalibTempCoeff::ReadCalibTempCoeffInfo - Error while reading input file; likely EOF..\n");
      return;
    }
    inputFile >> iSM;
    t->SetSuperModuleNum(iSM);

    // info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow >> tempCoeff >> iSrc;

      // check that input values are not out bounds
      if (iCol<0 || iCol>(AliEMCALGeoParams::fgkEMCALCols-1) ||
	  iRow<0 || iRow>(AliEMCALGeoParams::fgkEMCALRows-1) ) {
	printf("AliEMCALCalibTempCoeff::ReadCalibTempCoeffInfo - Error while reading input file; j %d iCol %d iRow %d\n", j, iCol, iRow);
      return;
      }

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      t->SetTC(iCol, iRow, tempCoeff);
      t->SetSrc(iCol, iRow, iSrc);
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTempCoeff::WriteTextCalibTempCoeffInfo(const TString &txtFileName,
					     Bool_t swapSides)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALCalibTempCoeff::WriteCalibTempCoeffInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  Float_t tempCoeff = 0;
  Int_t iSrc = 0;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibTempCoeff * t = (AliEMCALSuperModuleCalibTempCoeff*) fSuperModuleData[i];

    outputFile << t->GetSuperModuleNum() << endl;

    // info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      tempCoeff = t->GetTC(iCol, iRow);
      iSrc = t->GetSrc(iCol, iRow);

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow 
		 << " " << tempCoeff 
		 << " " << iSrc << endl;
    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTempCoeff::ReadRootCalibTempCoeffInfo(const TString &rootFileName,
					    Bool_t swapSides)
{
  //Read data from root file. ; coordinates given on SuperModule basis
  TFile inputFile(rootFileName, "read");  

  TTree *tree = (TTree*) inputFile.Get("tree");

  ReadTreeCalibTempCoeffInfo(tree, swapSides);

  inputFile.Close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTempCoeff::ReadTreeCalibTempCoeffInfo(TTree *tree,
					    Bool_t swapSides)
{
  // how many SuperModule's worth of info do we have?
  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  fNSuperModule = tree->GetEntries();

  Int_t iSM = 0; // SuperModule index
  // list of values to be read
  // info for each tower
  Float_t tempCoeff[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  Int_t iSrc[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  // end - all values

  // just to make the initializations of the arrays are done correctly, let's use memset
  memset(tempCoeff, 0, sizeof(tempCoeff)); 
  memset(iSrc, 0, sizeof(iSrc)); 

  // declare the branches
  tree->SetBranchAddress("iSM", &iSM);
  tree->SetBranchAddress("TempCoeff", tempCoeff);
  tree->SetBranchAddress("Src", iSrc);

  // indices for looping over the towers
  Int_t iCol = 0;
  Int_t iRow = 0;

  for (int ient=0; ient<tree->GetEntries(); ient++) {
    tree->GetEntry(ient);

    // assume the index SuperModules come in order: i=iSM
    AliEMCALSuperModuleCalibTempCoeff * t = (AliEMCALSuperModuleCalibTempCoeff*) fSuperModuleData[iSM];

    t->SetSuperModuleNum(iSM);

    // third: info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      // help variables: possibly modified or swapped indices
      int iColMod = iCol;
      int iRowMod = iRow;
      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iColMod = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRowMod = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      t->SetTC(iColMod, iRowMod, tempCoeff[iCol][iRow]);
      t->SetSrc(iColMod, iRowMod, iSrc[iCol][iRow]);
    }

  } // loop over entries

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTempCoeff::WriteRootCalibTempCoeffInfo(const TString &rootFileName,
					     Bool_t swapSides)
{
  // write data to root file. ; coordinates given on SuperModule basis
  TFile destFile(rootFileName, "recreate");  
  if (destFile.IsZombie()) {
    return;
  }  
  destFile.cd();

  TTree *tree = new TTree("tree","");

  // variables for filling the TTree
  Int_t iSM = 0; // SuperModule index
  // list of values to be written

  Float_t tempCoeff[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  Int_t iSrc[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows];   // end - all values

  // just to make the initializations of the arrays are done correctly, let's use memset
  memset(tempCoeff, 0, sizeof(tempCoeff)); 
  memset(iSrc, 0, sizeof(iSrc)); 

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  // for looping over towers
  Int_t iCol = 0;
  Int_t iRow = 0;

  // declare the branches
  // first
  tree->Branch("iSM", &iSM, "iSM/I");
  // info for each tower; see if a 2D array works OK or if we'll have to use 1D arrays instead 
  tree->Branch( "TempCoeff", &tempCoeff, Form("TempCoeff[%d][%d]/F", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );
  tree->Branch( "Src", &iSrc, Form("Src[%d][%d]/I", AliEMCALGeoParams::fgkEMCALCols, AliEMCALGeoParams::fgkEMCALRows) );

  for (iSM = 0; iSM < fNSuperModule; iSM++) {
    AliEMCALSuperModuleCalibTempCoeff * t = (AliEMCALSuperModuleCalibTempCoeff*) fSuperModuleData[iSM];

    iSM = t->GetSuperModuleNum();

    // info for each tower
    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      // help variables: possibly modified or swapped indices
      int iColMod = iCol;
      int iRowMod = iRow;
      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iColMod = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRowMod = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      tempCoeff[iColMod][iRowMod] = t->GetTC(iCol, iRow);
      iSrc[iColMod][iRowMod] = t->GetSrc(iCol, iRow);
    }

    tree->Fill();
  } // i, SuperModule

  tree->Write();
  destFile.Close();

  return;
}

//____________________________________________________________________________
AliEMCALCalibTempCoeff::~AliEMCALCalibTempCoeff()
{
  fSuperModuleData.Delete();
}

//____________________________________________________________________________
AliEMCALSuperModuleCalibTempCoeff * AliEMCALCalibTempCoeff::GetSuperModuleCalibTempCoeffNum(Int_t supModIndex)const
{ // getter via index
  for (int i=0; i<fNSuperModule; i++) {
    AliEMCALSuperModuleCalibTempCoeff * t = (AliEMCALSuperModuleCalibTempCoeff*) fSuperModuleData[i];
    if (t->GetSuperModuleNum() == supModIndex) {
      return t;
    }
  }

  // if we arrived here, then nothing was found.. just return a NULL pointer 
  return NULL;
}

