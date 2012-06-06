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
//  It produces a TList for the TObjArrays with the name of the Key
//    the TObjArray contain the Values, split from kommas, as found in the   //
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
#include <TList.h>
#include <TSystem.h>
//AliRoot includes

//header
#include "AliTPCConfigParser.h"

using std::ifstream;

AliTPCConfigParser::AliTPCConfigParser() :
TObject(),
fConfigMap(new TList),
fKeyIter(0),
fValIter(0)
{
 //
 // default constructor
 //
  fConfigMap->SetOwner();
} 
//_____________________________________________________________________
AliTPCConfigParser::AliTPCConfigParser(const AliTPCConfigParser &cfg) :
TObject(),
fConfigMap((TList*)cfg.fConfigMap->Clone()),
fKeyIter(0),
fValIter(0)
{
  //
  // copy constructor
  //
  fConfigMap->SetOwner();
}

//_____________________________________________________________________
AliTPCConfigParser::AliTPCConfigParser(const char* cfgfile) :
TObject(),
fConfigMap(new TList),
fKeyIter(0),
fValIter(0)
{
  //
  // default constructor using the config file name as input parameter
  //
  fConfigMap->SetOwner();
  if ( !cfgfile ) return;
  ParseConfigFileTxt(cfgfile);
}
//_____________________________________________________________________
AliTPCConfigParser& AliTPCConfigParser::operator = (const  AliTPCConfigParser &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTPCConfigParser(source);
  
  return *this;
}
//_____________________________________________________________________
AliTPCConfigParser::~AliTPCConfigParser()
{
 //
 // dtor
 //
  delete fConfigMap;
  delete fKeyIter;
  delete fValIter;
}
//_____________________________________________________________________
Int_t AliTPCConfigParser::ParseConfigFileTxt(const char* cfgfile)
{
 //
 // Function to parse a configuration file
 //
  ResetMap();
  ifstream file(gSystem->ExpandPathName(cfgfile));
  if ( !file.is_open() ){
    Error("ParseConfigFileTxt","File '%s' could not be opened!", cfgfile);
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
    line.ReplaceAll(", ",",");
    TObjArray *arrValues=line.Tokenize(" \t");
  //currently only name => Value is supported
  //and            name => 'nothing'
  //value may be a comma separated list, in which case a TObjArray
  //of the list will be created and stored as the value
    Int_t nentries=arrValues->GetEntries();
    if (nentries>2){
      Error("AliTPCConfigParser","ParseConfigFileTxt: Cannot parse line '%s'\n",line.Data());
      delete arrValues;
      continue;
    }
    TObjArray  *objArr=0x0;
    if (nentries==2){
      TObject *objVal=arrValues->At(1);
      const TString str=objVal->GetName();
      if (str.Contains(","))
        objArr=str.Tokenize(",");
      else{
        objArr=new TObjArray;
        objArr->Add(objVal->Clone());
      }
      objArr->SetOwner(kTRUE);
    } else {
      objArr=new TObjArray;
    }
    objArr->SetName(arrValues->At(0)->GetName());
    fConfigMap->AddLast(objArr);
    delete arrValues;
  }
  
  delete arr;
  return 0;
}
//_____________________________________________________________________
Float_t AliTPCConfigParser::GetValue(const char *key, UInt_t position)
{
  //
  //Get value for the speciefied key
  //
  TObject *val=((TObjArray*)fConfigMap->FindObject(key))->At(position);
  if ( !val ) return -999.;
  TString sval(((TObjString*)val)->GetString());
  return sval.Atof();
}
//_____________________________________________________________________
const char* AliTPCConfigParser::GetData(const char *key, UInt_t position)
{
  //
  //Get value for the speciefied key
  //
  TObjArray *arr=((TObjArray*)fConfigMap->FindObject(key));
  if (position>=(UInt_t)(arr->GetEntries())) return "";
  TObject *val=arr->At(position);
  if ( !val ) return "";
  return val->GetName();
}
//_____________________________________________________________________
Float_t AliTPCConfigParser::GetValue(const TObject *key, UInt_t position)
{
  //
  //Get value for the speciefied key
  //
  TObject *val=((TObjArray*)fConfigMap->FindObject(key))->At(position);
  if ( !val ) return -999.;
  TString sval(((TObjString*)val)->GetString());
  return sval.Atof();
}
//_____________________________________________________________________
const char* AliTPCConfigParser::GetData(const TObject *key, UInt_t position)
{
  //
  //Get value for the speciefied key
  //
  TObjArray *arr=((TObjArray*)fConfigMap->FindObject(key));
  if (position>=((UInt_t)arr->GetEntries())) return "";
  TObject *val=arr->At(position);
  if ( !val ) return "";
  return val->GetName();
}
//_____________________________________________________________________
Int_t AliTPCConfigParser::GetNumberOfValues(const char* key) const
{
  //
  // return the number of values for key
  //
  return ((TObjArray*)fConfigMap->FindObject(key))->GetEntries();
}
//_____________________________________________________________________
Int_t AliTPCConfigParser::GetNumberOfValues(TObject* key) const
{
  //
  // return the number of values for key
  //
  return ((TObjArray*)fConfigMap->FindObject(key))->GetEntries();
}
//_____________________________________________________________________
TObject* AliTPCConfigParser::NextKey(){
  if (!fKeyIter) fKeyIter=fConfigMap->MakeIterator();
  TObject *obj=fKeyIter->Next();
  if (!obj) {
    delete fKeyIter;
    fKeyIter=0;
  }
  return obj;
}
//_____________________________________________________________________
TObject* AliTPCConfigParser::NextValue(const char *key){
  return NextValueIter((TObjArray*)fConfigMap->FindObject(key));
}
//_____________________________________________________________________
TObject* AliTPCConfigParser::NextValue(TObject *key){
  return NextValueIter((TObjArray*)fConfigMap->FindObject(key));
}
//_____________________________________________________________________
TObject* AliTPCConfigParser::NextValueIter(TObjArray *key){
  if (!key) return 0;
  //check if the collection has changed
  if (fValIter && key!=fValIter->GetCollection()) {
    delete fValIter;
    fValIter=0x0;
  }
  if (!fValIter) fValIter=key->MakeIterator();
  TObject *value=fValIter->Next();
  if (!value) {
    delete fValIter;
    fValIter=0;
  }
  return value;
}
//_____________________________________________________________________
void AliTPCConfigParser::ResetMap()
{
  //
  // Reset the map with the configuration values
  //
  fConfigMap->Delete();
}
