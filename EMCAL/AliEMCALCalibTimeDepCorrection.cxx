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

// Objects of this class contain info on time-dependent corrections
//

#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "AliEMCALCalibTimeDepCorrection.h"

using namespace std;

ClassImp(AliEMCALCalibTimeDepCorrection)

//____________________________________________________________________________
AliEMCALCalibTimeDepCorrection::AliEMCALCalibTimeDepCorrection(const int nSM) : 
  fNSuperModule(nSM),
  fSuperModuleData(),
  fStartTime(0),
  fNTimeBins(0),
  fTimeBinSize(0)
{
  //Default constructor.
  for (int i=0; i<fNSuperModule; i++) {
    fSuperModuleData.Add(new AliEMCALSuperModuleCalibTimeDepCorrection(i));
  }
  fSuperModuleData.Compress(); // compress the TObjArray
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::InitCorrection(Int_t nSM, Int_t nBins, Float_t val)
{
  // This methods assumes that you are using SuperModules 0..nSM-1
  fNSuperModule = nSM;

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t nCorr = nBins;
  Float_t correction = val;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibTimeDepCorrection * t = (AliEMCALSuperModuleCalibTimeDepCorrection*) fSuperModuleData[i];
    iSM = i;
    t->SetSuperModuleNum(iSM); // assume SMs are coming in order

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      // set size of TArray
      t->GetCorrection(iCol,iRow)->Set(nCorr);
      for (Int_t k=0; k<nCorr; k++) {
	// add to TArray
	t->GetCorrection(iCol,iRow)->AddAt(correction, k); // AddAt = SetAt..
      }
    }

  } // i, SuperModule

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::SetCorrection(Int_t supModIndex, Int_t iCol, Int_t iRow, Int_t iBin, Float_t val)
{ // if you call for non-existing data, there may be a crash..
  ((AliEMCALSuperModuleCalibTimeDepCorrection*)fSuperModuleData[supModIndex])->GetCorrection(iCol,iRow)->AddAt(val, iBin); // AddAt = SetAt..
 return; 
}

//____________________________________________________________________________
Float_t AliEMCALCalibTimeDepCorrection::GetCorrection(Int_t supModIndex, Int_t iCol, Int_t iRow, Int_t iBin) const
{ // if you call for non-existing data, there may be a crash..
  return ((AliEMCALSuperModuleCalibTimeDepCorrection*)fSuperModuleData[supModIndex])->GetCorrection(iCol,iRow)->At(iBin);
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::ReadTextInfo(Int_t nSM, const TString &txtFileName,
						  Bool_t swapSides)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALCalibTimeDepCorrection::ReadTextInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNSuperModule = nSM;

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t nCorr = 0;
  Float_t correction = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  // read header info 
  inputFile >> fStartTime >> fNTimeBins >> fTimeBinSize;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibTimeDepCorrection * t = (AliEMCALSuperModuleCalibTimeDepCorrection*) fSuperModuleData[i];

    if (!inputFile) {
      printf("AliEMCALCalibTimeDepCorrection::ReadTextInfo - Error while reading input file; likely EOF..\n");
      return;
    }
    inputFile >> iSM;
    t->SetSuperModuleNum(iSM);

    for (Int_t j=0; j<nAPDPerSM; j++) {
      inputFile >> iCol >> iRow >> nCorr;

      // check that input values are not out bounds
      if (iCol<0 || iCol>(AliEMCALGeoParams::fgkEMCALCols-1) ||
	  iRow<0 || iRow>(AliEMCALGeoParams::fgkEMCALRows-1) || 
	  nCorr<0 ) {
	printf("AliEMCALCalibTimeDepCorrection::ReadTextInfo - Error while reading input file; j %d iCol %d iRow %d nCorr %d\n", j, iCol, iRow, nCorr);
      return;
      }

      // assume that this info is already swapped and done for this basis?
      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      // set size of TArray
      t->GetCorrection(iCol,iRow)->Set(nCorr);
      for (Int_t k=0; k<nCorr; k++) {
	inputFile >> correction;
	// add to TArray
	t->GetCorrection(iCol,iRow)->AddAt(correction, k);
      }
    }

  } // i, SuperModule

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::WriteTextInfo(const TString &txtFileName,
						   Bool_t swapSides)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALCalibTimeDepCorrection::WriteTextInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t nCorr = 0;

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  // write header info 
  outputFile << fStartTime << " " << fNTimeBins << " " << fTimeBinSize << endl;

  for (Int_t i = 0; i < fNSuperModule; i++) {
    AliEMCALSuperModuleCalibTimeDepCorrection * t = (AliEMCALSuperModuleCalibTimeDepCorrection*) fSuperModuleData[i];
    outputFile << t->GetSuperModuleNum() << endl;

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      nCorr = t->GetCorrection(iCol,iRow)->GetSize();

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      outputFile << iCol << " " << iRow << " " << nCorr << endl;
      for (Int_t k=0; k<nCorr; k++) {
	outputFile << t->GetCorrection(iCol,iRow)->At(k) << " ";
      }
      outputFile << endl;

    }

  } // i, SuperModule

  outputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::ReadRootInfo(const TString &rootFileName,
						  Bool_t swapSides)
{
  //Read data from root file. ; coordinates given on SuperModule basis
  TFile inputFile(rootFileName, "read");  

  TTree *treeGlob = (TTree*) inputFile.Get("treeGlob");
  TTree *treeCorr = (TTree*) inputFile.Get("treeCorr");

  ReadTreeInfo(treeGlob, treeCorr, swapSides);

  inputFile.Close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::ReadTreeInfo(TTree *treeGlob, TTree *treeCorr,
						  Bool_t swapSides)
{
 // how many SuperModule's worth of entries / APDs do we have?
  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;
  fNSuperModule = treeCorr->GetEntries() / nAPDPerSM;

  // global variables : only one entry
  treeGlob->SetBranchAddress("fStartTime", &fStartTime);
  treeGlob->SetBranchAddress("fNTimeBins", &fNTimeBins);
  treeGlob->SetBranchAddress("fTimeBinSize", &fTimeBinSize);
  treeGlob->GetEntry(0); 

  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  // list of values to be read
  Int_t nCorr = 0;
  Float_t correction[fgkMaxTimeBins] = {0};
  // make sure it's really initialized correctly
  memset(correction, 0, sizeof(correction)); // better safe than sorry
  // end - all values

  // declare the branches
  treeCorr->SetBranchAddress("iSM", &iSM);
  treeCorr->SetBranchAddress("iCol", &iCol);
  treeCorr->SetBranchAddress("iRow", &iRow);
  treeCorr->SetBranchAddress("nCorr", &nCorr);
  treeCorr->SetBranchAddress("correction", correction);

  for (int ient=0; ient<treeCorr->GetEntries(); ient++) {
    treeCorr->GetEntry(ient);

    // assume the index SuperModules come in order: i=iSM
    AliEMCALSuperModuleCalibTimeDepCorrection * t = (AliEMCALSuperModuleCalibTimeDepCorrection*) fSuperModuleData[iSM];
    t->SetSuperModuleNum(iSM);

    // assume that this info is already swapped and done for this basis?
    if (swapSides) {
      // C side, oriented differently than A side: swap is requested
      iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
      iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
    }


    // set size of TArray
    t->GetCorrection(iCol,iRow)->Set(nCorr);
    for (Int_t k=0; k<nCorr; k++) {
      // add to TArray
      t->GetCorrection(iCol,iRow)->AddAt(correction[k], k);
    }

  } // entry

  return;
}

//____________________________________________________________________________
void AliEMCALCalibTimeDepCorrection::WriteRootInfo(const TString &rootFileName,
						   Bool_t swapSides)
{
  // write data to root file. ; coordinates given on SuperModule basis
  TFile destFile(rootFileName, "recreate");  
  if (destFile.IsZombie()) {
    return;
  }  
  destFile.cd();

  TTree *treeGlob = new TTree("treeGlob","");
  TTree *treeCorr = new TTree("treeCorr","");

  // global part only has one entry
  treeGlob->Branch("fStartTime", &fStartTime, "fStartTime/i"); // unsigned int..
  treeGlob->Branch("fNTimeBins", &fNTimeBins, "fNTimeBins/I");
  treeGlob->Branch("fTimeBinSize", &fTimeBinSize, "fTimeBinSize/I");
  treeGlob->Fill();

  // variables for filling the TTree
  Int_t iSM = 0; // SuperModule index
  Int_t iCol = 0;
  Int_t iRow = 0;
  Int_t nCorr = 0;
  Float_t correction[fgkMaxTimeBins] = {0};
  // make sure it's really initialized correctly
  memset(correction, 0, sizeof(correction)); // better safe than sorry

  // declare the branches
  treeCorr->Branch("iSM", &iSM, "iSM/I");
  treeCorr->Branch("iCol", &iCol, "iCol/I");
  treeCorr->Branch("iRow", &iRow, "iRow/I");
  treeCorr->Branch("nCorr", &nCorr, "nCorr/I");
  treeCorr->Branch("correction", &correction, "correction[nCorr]/F");

  Int_t nAPDPerSM = AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows;

  for (iSM = 0; iSM < fNSuperModule; iSM++) {
    AliEMCALSuperModuleCalibTimeDepCorrection * t = (AliEMCALSuperModuleCalibTimeDepCorrection*) fSuperModuleData[iSM];

    for (Int_t j=0; j<nAPDPerSM; j++) {
      iCol = j / AliEMCALGeoParams::fgkEMCALRows;
      iRow = j % AliEMCALGeoParams::fgkEMCALRows;

      nCorr = t->GetCorrection(iCol,iRow)->GetSize();
      if (nCorr > fgkMaxTimeBins) {
	printf("AliEMCALCalibTimeDepCorrection::WriteRootInfo - too many correction/timebins %d kept\n", nCorr);
	return;
      }

      if (swapSides) {
	// C side, oriented differently than A side: swap is requested
	iCol = AliEMCALGeoParams::fgkEMCALCols-1 - iCol;
	iRow = AliEMCALGeoParams::fgkEMCALRows-1 - iRow;
      }

      for (Int_t k=0; k<nCorr; k++) {
	correction[k] = t->GetCorrection(iCol,iRow)->At(k);
      }

      treeCorr->Fill();
    }

  } // i, SuperModule

  treeGlob->Write();
  treeCorr->Write();
  destFile.Close();

  return;
}

//____________________________________________________________________________
AliEMCALCalibTimeDepCorrection::~AliEMCALCalibTimeDepCorrection()
{
  fSuperModuleData.Delete();
}

//____________________________________________________________________________
AliEMCALSuperModuleCalibTimeDepCorrection * AliEMCALCalibTimeDepCorrection::GetSuperModuleCalibTimeDepCorrectionNum(Int_t supModIndex)const
{ // getter via index
  for (int i=0; i<fNSuperModule; i++) {
    AliEMCALSuperModuleCalibTimeDepCorrection * t = (AliEMCALSuperModuleCalibTimeDepCorrection*) fSuperModuleData[i];
    if (t->GetSuperModuleNum() == supModIndex) {
      return t;
    }
  }

  // if we arrived here, then nothing was found.. just return a NULL pointer 
  return NULL;
}

