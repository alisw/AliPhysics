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

// Objects of this class read txt file with APD number data
//

#include <fstream>
#include <TString.h>

#include "AliEMCALMapAPD.h"

const int kFirstAPD = 10213; // dummy number, only used for testing

ClassImp(AliEMCALMapAPD)

//____________________________________________________________________________
AliEMCALMapAPD::AliEMCALMapAPD() : 
  fNSuperModule(0),
  fSuperModuleData(0)
{
  //Default constructor.
}

//____________________________________________________________________________
void AliEMCALMapAPD::ReadMapAPDInfoStripBasis(Int_t nSM, const TString &txtFileName)
{
  //Read data from txt file; coordinates given on StripModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALMapAPD::ReadMapAPDInfoStripBasis - Cannot open the APD info file %s\n",txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;
  if (fSuperModuleData) delete [] fSuperModuleData;
  fSuperModuleData = new AliEMCALSuperModuleMapAPD[fNSuperModule];

  Int_t iSM = 0; // SuperModule index
  Int_t iAPD = 0;
  // info in map is based on Strip Info
  Int_t iStrip = 0;
  Int_t iStripCol = 0;
  Int_t iStripRow = 0;

  // we'll convert this into SuperModule Info
  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerSM = fgkEmCalCols * fgkEmCalRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleMapAPD &t = fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALMapAPD::ReadMapAPDInfoStripBasis - Error while reading input file; likely EOF..\n");
      return;
    }
    inputFile >> iSM;
    t.fSuperModuleNum = iSM;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iStrip >> iStripCol >> iStripRow >> iAPD;
      // iStrip is a number in the range 0..23 (number of StripModules per SuperModule)
      // iStripCol is a number in the range 0..1 (number of tower columns per StripModule)
      // iStripRow is a number in the range 0..23 (number of tower rows per StripModule)
      iCol = iStrip*2 + iStripCol;
      iRow = iStripRow;

      if (iSM%2 == 1) { // C side, oriented differently than A side: swap..
	iCol = fgkEmCalCols-1 - iCol;
	iRow = fgkEmCalRows-1 - iRow;
      }

      t.fAPDNum[iCol][iRow] = iAPD;
    }

  } // i, SuperModule

  inputFile.close();

  return;
}


//____________________________________________________________________________
void AliEMCALMapAPD::ReadMapAPDInfoSingleStripBasis(Int_t iSM, Int_t iStrip, const TString &txtFileName)
  // iSM is the SuperModule number
  // iStrip is a number in the range 0..23 (number of StripModules per SuperModule)
{
  //Read data from txt file; coordinates given on StripModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALMapAPD::ReadMapAPDInfoSingleStripBasis - Cannot open the APD info file %s\n",txtFileName.Data());
    return;
  }

  // see if there is an existing SuperModule with the right index
  Int_t foundSM = -1;
  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleMapAPD &t = fSuperModuleData[i];
    if (t.fSuperModuleNum == iSM) foundSM = i;
  }

  if (foundSM == -1) {
    printf("AliEMCALMapAPD::ReadMapAPDInfoSingleStripBasis - no SuperModule %d found!\n", iSM);
    return;
  }

  AliEMCALSuperModuleMapAPD &t = fSuperModuleData[foundSM];

  Int_t iAPD = 0;
  // info in map is based on Strip Info
  Int_t iStripCol = 0;
  Int_t iStripRow = 0;

  // we'll convert this into SuperModule Info
  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerStrip = 2 * fgkEmCalRows; // 2 columns in a strip

  for (Int_t j=0; j<nAPDPerStrip; j++) {
    if (!inputFile) {
      printf("AliEMCALMapAPD::ReadMapAPDInfoSingleStripBasis - Error while reading input file; likely EOF..\n");
      return;
    }
    inputFile >> iStripCol >> iStripRow >> iAPD;
    // iStripCol is a number in the range 0..1 (number of tower columns per StripModule)
    // iStripRow is a number in the range 0..23 (number of tower rows per StripModule)
    iCol = iStrip*2 + iStripCol;
    iRow = iStripRow;

    /*
// For the SuperModule calibration we will typically use SuperModule 0
// for all - so don't worry about any swaps for now.. 
// I.e. we'll work in a local column and row coord. system, not necessarily ALICE
// May revisit later..
    if (iSM%2 == 1) { // C side, oriented differently than A side: swap..
      iCol = fgkEmCalCols-1 - iCol;
      iRow = fgkEmCalRows-1 - iRow;
    }
    */

    t.fAPDNum[iCol][iRow] = iAPD;
  }

  inputFile.close();

  return;
}


//____________________________________________________________________________
void AliEMCALMapAPD::ReadMapAPDInfo(Int_t nSM, const TString &txtFileName)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALMapAPD::ReadMapAPDInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;
  if (fSuperModuleData) delete [] fSuperModuleData;
  fSuperModuleData = new AliEMCALSuperModuleMapAPD[fNSuperModule];

  Int_t iSM = 0; // SuperModule index
  Int_t iAPD = 0;
  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerSM = fgkEmCalCols * fgkEmCalRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleMapAPD &t = fSuperModuleData[i];
    if (!inputFile) {
      printf("AliEMCALMapAPD::ReadMapAPDInfo - Error while reading input file; likely EOF..");
      return;
    }
    inputFile >> iSM;
    t.fSuperModuleNum = iSM;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow >> iAPD;

      // assume that this info is already swapped and done for this basis?
      /*
      if (iSM%2 == 1) { // C side, oriented differently than A side: swap..
	iCol = fgkEmCalCols-1 - iCol;
	iRow = fgkEmCalRows-1 - iRow;
      }
      */

      t.fAPDNum[iCol][iRow] = iAPD;
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALMapAPD::WriteMapAPDInfo(const TString &txtFileName)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALMapAPD::WriteMapAPDInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerSM = fgkEmCalCols * fgkEmCalRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleMapAPD &t = fSuperModuleData[i];
    outputFile << t.fSuperModuleNum << endl;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / fgkEmCalRows;
      iRow = j % fgkEmCalRows;
      outputFile << iCol << " " << iRow << " " 
		 << t.fAPDNum[iCol][iRow] << endl;
    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
AliEMCALMapAPD::~AliEMCALMapAPD()
{
  delete [] fSuperModuleData;
}

//____________________________________________________________________________
AliEMCALMapAPD::AliEMCALSuperModuleMapAPD AliEMCALMapAPD::GetSuperModuleMapAPDId(Int_t supModIndex)const
{
  AliEMCALSuperModuleMapAPD t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fSuperModuleData)
    return t;

  return fSuperModuleData[supModIndex];
}

//____________________________________________________________________________
AliEMCALMapAPD::AliEMCALSuperModuleMapAPD AliEMCALMapAPD::GetSuperModuleMapAPDNum(Int_t supModIndex)const
{
  AliEMCALSuperModuleMapAPD t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fSuperModuleData)
    return t;

  for (int i=0; i<fNSuperModule; i++) {
    if (fSuperModuleData[i].fSuperModuleNum == supModIndex) {
      return fSuperModuleData[i];
    }
  }

  return t;
}

//____________________________________________________________________________
void AliEMCALMapAPD::GenerateDummyAPDInfo(Int_t nSM, Int_t * iSM)
{
  // just a temporary method to create some info to exercise I/O

  fNSuperModule = nSM;
  if (fSuperModuleData) delete [] fSuperModuleData;
  fSuperModuleData = new AliEMCALSuperModuleMapAPD[fNSuperModule];

  Int_t iAPD = 0;
  Int_t iCol = 0;
  Int_t iRow = 0;

  Int_t nAPDPerSM = fgkEmCalCols * fgkEmCalRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleMapAPD &t = fSuperModuleData[i];
    t.fSuperModuleNum = iSM[i]; // set SuperModules in Order

    for (Int_t j=0; j<nAPDPerSM; j++) {

      iCol = j / fgkEmCalRows;
      iRow = j % fgkEmCalRows;
      iAPD = j + kFirstAPD + i*nAPDPerSM; // just a dummy number; assuming all APDs are assigned in order from Catania.. 

      t.fAPDNum[iCol][iRow] = iAPD;
    }

  } // i, SuperModule

  return;
}

//____________________________________________________________________________
int AliEMCALMapAPD::CheckForDuplicates()
{ 
  // keep it simple: have one big array with place
  // for all APDs from Catania (10000-19999 max) Houston (20000-29999 max)
  // - and see how many times each APD occurs

  const int kMaxAPDNum = 30000;
  int counter[kMaxAPDNum] = {0};
  for (int k=0; k<kMaxAPDNum; k++) { 
    counter[k] = 0; 
  }

  Int_t nAPDPerSM = fgkEmCalCols * fgkEmCalRows;

  // go through all APDs
  int iCol, iRow;
  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleMapAPD &t = fSuperModuleData[i];
    
    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / fgkEmCalRows;
      iRow = j % fgkEmCalRows;
      counter[t.fAPDNum[iCol][iRow]]++;
    }
  } // i, SuperModule

  int nProblems = 0;
  for (int k=0; k<kMaxAPDNum; k++) { 
    if (counter[k] > 1) {
      printf("AliEMCALMapAPD::CheckForDuplicates - APD %d occurs more than once!\n",k);
      nProblems++;
    } 
  }

  return nProblems;
}
