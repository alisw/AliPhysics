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

#include "AliEMCALCalibAPD.h"

const int kMaxLen = 1000; // maximum length of single line (# of characters)
// an OK line with complete info should have a certain number of characters
const int kMinLenAPDLine = 135;
const int kMaxLenAPDLine = 1000;

ClassImp(AliEMCALCalibAPD)

//____________________________________________________________________________
AliEMCALCalibAPD::AliEMCALCalibAPD() : 
  fNCalibAPD(0),
  fData(0)
{
  //Default constructor.
}

//____________________________________________________________________________
void AliEMCALCalibAPD::ReadCalibAPDInfo(Int_t nAPD, const TString &txtFileName)
{
  //Read data from txt file. ; coordinates given on SuperModule basis

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    printf("AliEMCALCalibAPD::ReadCalibAPDInfo - Cannot open the APD info file %s\n", txtFileName.Data());
    return;
  }

  fNCalibAPD = nAPD;
  if (fData) delete [] fData;
  fData = new AliEMCALCalibAPDData[fNCalibAPD];

  char line[kMaxLen];

  /* DS: header lines skipped when switched to white-spaced separeted dat files instead of csv:
     conversion from spreadsheet can be done a la
     sed 's/,/ /g' APD-database-Houston.csv | egrep ^2 | awk '{if (NF==19) {print $0}}' > APD-database-Houston.dat
     - meaning "replace , with whitespace, only get the columns that start with the number 2 (Houston APDs),
     and only get rows with the full 19 fields

  // get header lines:
  inputFile.getline(line, kMaxLen);
  //  printf(" 1st header line character count %d\n", inputFile.gcount());
  inputFile.getline(line, kMaxLen);
  //  printf(" 2nd header line character count %d\n", inputFile.gcount());
  */

  // variables for reading
  int i1,i2,i3,i4,i5;
  float f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12;
  char c1[10], c2[20];

  int j = 0; // index for the ones that were read OK
  for (Int_t i = 0; i < fNCalibAPD; i++) {
    AliEMCALCalibAPDData &t = fData[j];
    if (!inputFile) {
      printf("AliEMCALCalibAPD::ReadCalibAPDInfo - Error while reading input file; likely EOF..\n");
      fNCalibAPD = j; // that's how many we actually read succesfully
      printf("AliEMCALCalibAPD::ReadCalibAPDInfo - read %d OK\n", fNCalibAPD);
      return;
    }

    // use some temporary/local values to perhaps help with the order when 
    // trying to read all the many fields in a line..
    inputFile.getline(line, kMaxLen);
    int nchar = inputFile.gcount();
    //printf(" line %d ok %d - character count %d\n", i, j, nchar);

    if (nchar>kMinLenAPDLine && nchar<kMaxLenAPDLine) {
      // looks like the line has about the right number of characters, let's
      // try to decode it now..

      //      printf("input: %s\n",line);
      sscanf(line, "%d %u %s %s %d %d %f %f %f %f %f %f %f %f %f %d %f %f %f",
	     &i1, &i2, c1, c2, &i3, &i4, // header-type info
	     &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, // measurements
	     &i5, &f10, &f11, &f12); // Hamamatsu

      //      printf("after scanf: %s\n",line);
      /*
      printf("%d,%u,%s,%s,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%g,%13.11f,%g\n",
	     i1, i2, c1, c2, i3, i4, // header-type info
	     f1, f2, f3, f4, f5, f6, f7, f8, f9, // measurements
	     i5, f10, f11, f12); // Hamamatsu
      */
      // assign the variables:
      t.fAPDNum = i1; 
      t.fSerialNum = i2; 
      sprintf(t.fStatus,"%s",c1);
      sprintf(t.fLocation,"%s",c2);
      t.fRunNum = i3; 
      t.fTestPos = i4; 

      t.fV30 = f1; 
      t.fV50 = f2; 
      t.fVoltCoeff = f3;
      t.fPar[0] = f4;
      t.fPar[1] = f5;
      t.fPar[2] = f6;
      t.fParErr[0] = f7;
      t.fParErr[1] = f8;
      t.fParErr[2] = f9;

      t.fBreakDown = i5; 
      t.fHamV50 = f10;
      t.fDarkCurrent = f11;
      t.fTestTemp = f12;

      j++; // increment our 'OK' counter..
    } // line length appears OK      
  } // i, APD

  inputFile.close();

  return;
}

//____________________________________________________________________________
void AliEMCALCalibAPD::WriteCalibAPDInfo(const TString &txtFileName)
{
  // write data to txt file. ; coordinates given on SuperModule basis

  std::ofstream outputFile(txtFileName.Data());
  if (!outputFile) {
    printf("AliEMCALCalibAPD::WriteCalibAPDInfo - Cannot open the APD output file %s\n", txtFileName.Data());
    return;
  }

  char *comma = ",";

  for (Int_t i = 0; i < fNCalibAPD; i++) {
    AliEMCALCalibAPDData &t = fData[i];

    outputFile << t.fAPDNum << comma
	       << t.fSerialNum << comma // Serial Number; from Hamamatsu
	       << t.fStatus << comma //
	       << t.fLocation << comma //
	       << t.fRunNum << comma
	       << t.fTestPos << comma
	       << t.fV30 << comma // Catania/Houston Voltage V30 (V) at T = 25 deg C
	       << t.fV50 << comma 
	       << t.fVoltCoeff << comma // 1/M x dM/dV
	       << t.fPar[0] << comma // fit parameters, p0,p1,p2 - for ADC vs bias measurement
	       << t.fPar[1] << comma
	       << t.fPar[2] << comma
	       << t.fParErr[0] << comma // error on fit parameters	
	       << t.fParErr[1] << comma 
	       << t.fParErr[2] << comma 
	       << t.fBreakDown << comma // Hamamatsu Breakdown Voltage (V)	
	       << t.fHamV50 << comma; // Hamamatsu Voltage V50 (V)
    // I wasn't able to quite reproduce the values as they appeared on the 
    // original file: e.g. dark current is not always 11-field - if last digit is zero..
    // the other floats have 6 significant digits, but sometimes switch to 
    // scientific notation - don't know how to handle this for varying length of fields.. -leave it as is for now..
    outputFile << t.fDarkCurrent << comma; // Hamamatsu Dark Current (A)	
    /*
    // some tweaking for the dark-current field to get the output the same
    // as on the original CSV
    outputFile.precision(11);	
    outputFile << fixed << t.fDarkCurrent << comma; // Hamamatsu Dark Current (A)	
    // back to normal..
    outputFile.precision(6); outputFile.unsetf(ios_base::floatfield);	
    */

    outputFile << t.fTestTemp // Hamamatsu Testing Temperature (deg C)	
	       << endl; 

  } // i, APD

  outputFile.close();

  return;
}

//____________________________________________________________________________
AliEMCALCalibAPD::~AliEMCALCalibAPD()
{
  delete [] fData;
}

//____________________________________________________________________________
AliEMCALCalibAPD::AliEMCALCalibAPDData AliEMCALCalibAPD::GetCalibAPDDataId(Int_t apdIndex)const
{
  AliEMCALCalibAPDData t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fData)
    return t;

  return fData[apdIndex];
}

//____________________________________________________________________________
AliEMCALCalibAPD::AliEMCALCalibAPDData AliEMCALCalibAPD::GetCalibAPDDataNum(Int_t apdNum)const
{
  AliEMCALCalibAPDData t;  // just to maybe prevent a crash, but we are returning something not-initialized so maybe not better really..
  if (!fData)
    return t;

  for (int i=0; i<fNCalibAPD; i++) {
    if (fData[i].fAPDNum == apdNum) {
      return fData[i];
    }
  }

  // if we made it to here, then no match was found - just return then
  return t;
}


