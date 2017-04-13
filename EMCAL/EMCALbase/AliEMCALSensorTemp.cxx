/**************************************************************************
 * Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
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

#include <strings.h>

#include "AliEMCALSensorTemp.h"

/// \cond CLASSIMP
ClassImp(AliEMCALSensorTemp) ;
/// \endcond

///
/// Default constructor
//______________________________________________________________________________________________
AliEMCALSensorTemp::AliEMCALSensorTemp(): AliDCSSensor(),
  fSide(0),
  fSector(0),
  fNum(0)
{ }

///
///  Copy constructor
//______________________________________________________________________________________________
AliEMCALSensorTemp::AliEMCALSensorTemp(const AliEMCALSensorTemp& source) :
  AliDCSSensor(source),
   fSide(source.fSide),
   fSector(source.fSector),
   fNum(source.fNum)
{ }

///
/// Assignment operator
//______________________________________________________________________________________________
AliEMCALSensorTemp& AliEMCALSensorTemp::operator=(const AliEMCALSensorTemp& source)
{
  if (&source == this) return *this;
  
  new (this) AliEMCALSensorTemp(source);
  
  return *this;  
}

///
/// Read values from ascii file
//______________________________________________________________________________________________
TClonesArray * AliEMCALSensorTemp::ReadList(const char *fname,
                                            const TString& amandaString) 
{
  TTree * tree = new TTree("asci","asci");
  tree->ReadFile(fname,"");
  TClonesArray *arr = ReadTree(tree, amandaString);
  delete tree;
  return arr;
}

///
/// Read selected info from TTree
//______________________________________________________________________________________________
TClonesArray * AliEMCALSensorTemp::ReadTree(TTree *tree, 
                                            const TString& amandaString) 
{
  Int_t nentries = tree->GetEntries();
  Int_t sensor=0;
  Int_t sector=0;
  char  side[100];
  Int_t num=0;
  Int_t echa=0;
  //Double_t x=0;
  //Double_t y=0;
  //Double_t z=0;
  //String_t namedtp[100];
  
  tree->SetBranchAddress("Sensor",&sensor);
  tree->SetBranchAddress("Side",&side);
  tree->SetBranchAddress("Sec",&sector);
  tree->SetBranchAddress("Num",&num);
  tree->SetBranchAddress("ECha",&echa);
  //tree->SetBranchAddress("X",&x);
  //tree->SetBranchAddress("Y",&y);
  //tree->SetBranchAddress("Z",&z);
  
  // firstSensor = (Int_t)tree->GetMinimum("ECha");
  // lastSensor = (Int_t)tree->GetMaximum("ECha");
  
  TClonesArray * array = new TClonesArray("AliEMCALSensorTemp",nentries);
  
  for (Int_t isensor=0; isensor<nentries; isensor++)
  {
    AliEMCALSensorTemp * temp = new ((*array)[isensor])AliEMCALSensorTemp;
    tree->GetEntry(isensor);
    
    temp->SetId(sensor);
    temp->SetIdDCS(echa);
    
    TString stringID = Form (amandaString.Data(),echa);
    temp->SetStringID(stringID);
    
    if (side[0]=='C') temp->SetSide(1);
    temp->SetSector(sector);
    temp->SetNum(num);
    
    // Don't yet know the local or global coordinates for where the sensors will be placed..
    //temp->SetX(x);
    //temp->SetY(y);
    //temp->SetZ(z);
  }
  
  return array;
}
