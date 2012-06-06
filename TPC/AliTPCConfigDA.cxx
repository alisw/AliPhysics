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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Class for Parsing simple text configuration files                        //
//  It produces a TMap for the Key=>Value pairs found in the                 //
//  Configutation file.                                                      //
//                                                                           //
// The configuration file has a simple structure:                            //
// * Lines starting with a # or empty lines are ignored                      //
// * Key and Value are separated either by a <tab> or <space>es              //
//                                                                           //
//  Currently the class is used in the TPC DAs to allow an adjustment of     //
//  the most relevant parameters without recompiling the DAs                 //
///////////////////////////////////////////////////////////////////////////////


#include <fstream>
//Root includes
#include <TObjString.h>
#include <TObjArray.h>
#include <TString.h>
#include <TIterator.h>
#include <TMap.h>
//AliRoot includes

//header
#include "AliTPCConfigDA.h"

using std::ifstream;

AliTPCConfigDA::AliTPCConfigDA() :
  TObject(),
  fConfigMap(new TMap)
{
 //
 // default constructor
 //
 fConfigMap->SetOwnerKeyValue();
} 
//_____________________________________________________________________
AliTPCConfigDA::AliTPCConfigDA(const AliTPCConfigDA &cfg) :
  TObject(),
  fConfigMap((TMap*)cfg.fConfigMap->Clone())
{
  //
  // copy constructor
  //
  fConfigMap->SetOwnerKeyValue();
}

//_____________________________________________________________________
AliTPCConfigDA::AliTPCConfigDA(const char* cfgfile) :
  TObject(),
  fConfigMap(new TMap)
{
  //
  // default constructor using the config file name as input parameter
  //
  fConfigMap->SetOwnerKeyValue();
  if ( !cfgfile ) return;
  ParseConfigFileTxt(cfgfile);
}
//_____________________________________________________________________
AliTPCConfigDA& AliTPCConfigDA::operator = (const  AliTPCConfigDA &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTPCConfigDA(source);

  return *this;
}
//_____________________________________________________________________
AliTPCConfigDA::~AliTPCConfigDA()
{
 //
 // dtor
 //   
  delete fConfigMap;
}
//_____________________________________________________________________
Int_t AliTPCConfigDA::ParseConfigFileTxt(const char* cfgfile)
{
 //
 // Function to parse a configuration file
 //

 ifstream file(cfgfile);
 if ( !file.is_open() ){
  Error("ParseConfigFileTxt","File %s could not be opened!", cfgfile);
  return 1;
 } 
 TString strFile;
 strFile.ReadFile(file);
 TObjArray *arr=strFile.Tokenize("\n");
 if ( !arr ) {
   file.close(); 
   return 2;
 }
 TIter nextLine(arr);
 while (TObject *l=nextLine()){
  TString line(((TObjString*)l)->GetString());
  //remove whitespcaces 
  line.Remove(TString::kBoth,' ');
  line.Remove(TString::kBoth,'\t');
  if ( line.BeginsWith("#") || line=="" ) continue;
  TObjArray *arrValues=line.Tokenize(" \t");
  //currently only name => Value is supported
  if (arrValues->GetEntries()!=2){
    printf("AliTPCConfigDA::ParseConfigFileTxt: Cannot parse line '%s'\n",line.Data());
    delete arrValues;
    continue;
  }
  fConfigMap->Add(arrValues->At(0)->Clone(),arrValues->At(1)->Clone());
  delete arrValues; 
 } 

 delete arr;
 return 0;
}
//_____________________________________________________________________
Float_t AliTPCConfigDA::GetValue(const char *key) const
{
  //
  //Get value for the speciefied key
  //
  TObject *val=fConfigMap->GetValue(key);
  if ( !val ) return -999.;
  TString sval(((TObjString*)val)->GetString());
  return sval.Atof(); 
}
//_____________________________________________________________________
void AliTPCConfigDA::ResetMap()
{
  //
  // Reset the map with the configuration values
  //
  fConfigMap->DeleteAll();
}
