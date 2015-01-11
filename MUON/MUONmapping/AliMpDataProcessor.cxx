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

// $Id$
//
// Category: basic

//-----------------------------------------------------------------------------
// Class AliMpDataProcessor
// -----------------------
// Class for converting ASCII data files in the map of string
// Author: Ivana Hrivnacova, IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpFiles.h"

#include "AliLog.h"

#include <TSystem.h>
#include <TObjString.h>
#include <TFile.h>
#include <Riostream.h>

#include <fstream>
#include <sstream>
#include <map>

//
// static private methods
//

//_____________________________________________________________________________
const TString& AliMpDataProcessor::GetHeaderFileName() 
{
  /// Return the default name for generated header file
  static const TString kHeaderFileName = "add.h";
  return kHeaderFileName;
}    

//_____________________________________________________________________________
const TString& AliMpDataProcessor::GetImplFileName() 
{
  /// Return the default name for generated implementation file
  static const TString kImplFileName = "add.cxx";
  return kImplFileName;
}  

//
// ctors, dtor
//

//_____________________________________________________________________________
AliMpDataProcessor::AliMpDataProcessor()
  : TObject(),
    fCounter(0),
    fHeaderFile(),
    fImplFile()
{
/// Default and standar constructor

  fHeaderFile.open(GetHeaderFileName().Data(), std::ios::out);
  fImplFile.open(GetImplFileName().Data(), std::ios::out);
  
  // Add check

}

//_____________________________________________________________________________
AliMpDataProcessor:: ~AliMpDataProcessor()
{
/// Destructor

}

//
// private methods
//


//_____________________________________________________________________________
void AliMpDataProcessor::ProcessDirectory(const TString& path, 
                                          AliMpDataMap* dataMap) 
{
/// Recursive function to process data directory
  
  AliDebugStream(2) << "ProcessDirectory " << path.Data() << endl;

  // skip some directories
  //if ( path.Contains("manu_serial") ) return;
  //if ( path.Contains("stationTrigger") ) return;
  
  gSystem->cd(path);
  gSystem->Exec("ls > /tmp/mpfiles.txt");
  
  ifstream mpfiles("/tmp/mpfiles.txt");
  TString mpfilesStr;
  TString fileName;
  while ( mpfiles >> fileName ) {
    mpfilesStr += fileName.Data();
    mpfilesStr += " ";
  }  

  istringstream mpFilesStream(mpfilesStr.Data());
  while ( mpFilesStream >> fileName ) {
    TString filePath = path + "/" + fileName;
    if ( gSystem->OpenDirectory(filePath.Data()) )
      ProcessDirectory(filePath, dataMap);
    else
      ProcessFile(filePath, dataMap);
  }  

  AliDebugStream(2) << "ProcessDirectory done " << path.Data() << endl;
} 

//_____________________________________________________________________________
void AliMpDataProcessor::ProcessFile(const TString& path, AliMpDataMap* dataMap) 
{
/// Dump the content of the file specified by its path
/// in a string and fill it in the dataMap

  // cut top path
  string top = AliMpFiles::GetTop().Data();
  string fullDataPath = path.Data();
  string dataPath = fullDataPath;
  if ( dataPath.find(top) != string::npos ) dataPath.erase(0, top.size()+1);
  dataPath.erase(0,dataPath.find('/')+1); 
  
  AliDebugStream(2) << "Processing file " << dataPath << endl;

  // skip macros
  if ( dataPath.find(".C") != std::string::npos || 
       dataPath.find(".h") != std::string::npos ) return;
 
  ifstream input(path.Data(), ios::in);
  if ( ! input.is_open() ) {
    cerr << "Cannot open input file " << path.Data() << endl;
    return;
  }  

  string fileString;
  string line;

  do {
    getline(input, line);
    if ( ! input.eof() ) { 
      fileString += line;
      fileString += '\n';
    }  
  }
  while ( ! input.eof() );

  dataMap->Add(dataPath, fileString);
  // dataMap->Add(new TObjString(dataPath.c_str()), new TObjString(fileString.c_str()));
}  

//_____________________________________________________________________________
void AliMpDataProcessor::GenerateFunction(const TString& path, 
                                          const TString& data) 
{
/// Generate a C++ function which defines a string with the data content 
/// and map it to the given path in the map

  AliDebugStream(2) << "Processing path " << path.Data() << endl;

  // skip motif & padPos
  //if ( dataPath.find("motif") != std::string::npos || 
  //     dataPath.find("padPos") != std::string::npos ) return;
 
  // add function declaration in the header file name
  fHeaderFile << "    void AddData" << ++fCounter << "();" << endl;

  // add function implementation in .cxx
  fImplFile 
    << "//_____________________________________________________________________________" << endl
    << "void AliMpDataMap::AddData" << fCounter << "()" << endl
    << "{" << endl
    << "/// Automatically generated function for data file: " << endl
    << "/// " << path << endl
    << endl
    << "  string path = \"" << path << "\";" << endl << endl;
    
  GenerateFileCode(data);
  
  fImplFile 
    << endl
    << "  fMap[path] = data;" << endl
    << "}" << endl
    << endl;      
} 

//_____________________________________________________________________________
void AliMpDataProcessor::GenerateFileCode(const TString& data) 
{
/// Dump the content of the file specified by its path as a C++ string

  istringstream in(data.Data());

  string line;
  getline(in, line);
  fImplFile << "  string data = \""; 
  //for ( unsigned int i=line.size(); i<58; i++ ) line = line + " ";
  fImplFile << line << " \\n\";";
  if ( in.eof() ) return; 

  do {
    getline(in, line);
    if ( ! in.eof() ) { 
      fImplFile << endl << "  data = data + \"";
      // for ( unsigned int i=line.size(); i<58; i++ ) line = line + " ";
      fImplFile << line << " \\n\";";
    }  
  }
  while ( ! in.eof() );

  fImplFile << endl;
}  

//_____________________________________________________________________________
void AliMpDataProcessor::GenerateFill()
{
/// Generate function which calls all previously generated functions

  // add function declaration in the header file name
  fHeaderFile << endl << "    void Fill();" << endl;

  // add function implementation in .cxx
  fImplFile 
    << "//_____________________________________________________________________________" << endl
    << "void AliMpDataMap::Fill()" << endl
    << "{" << endl
    << "/// Automatically generated function for calling all generated functions " << endl
    << endl;
    
  for ( Int_t i=1; i<=fCounter; ++i ) 
   fImplFile << "  AddData" << i << "();" << endl;   
      
  fImplFile << "}" << endl;
}  
      

//
// public methods
//


//_____________________________________________________________________________
AliMpDataMap* AliMpDataProcessor::CreateDataMap(const TString& dataDir) 
{
/// Process data directory and map a string with the content of each file to
/// the file path.

  TString curDir = gSystem->pwd();
  
  AliMpDataMap* dataMap = new AliMpDataMap();

  TString dataPath = AliMpFiles::GetTop() + "/" + dataDir;
  ProcessDirectory(dataPath, dataMap); 
  
  gSystem->cd(curDir); 
  
  return dataMap;
}
  

//_____________________________________________________________________________
Bool_t AliMpDataProcessor::GenerateData(AliMpDataMap* dataMap,
                                        const TString& outputDataDir) 
{
/// Generate ASCII data files in outputDataDir from dataMap

  TString curDir = gSystem->pwd();  

  TString outputDataPath = AliMpFiles::GetTop() + "/" + outputDataDir;
  if ( gSystem->OpenDirectory(outputDataPath.Data()) ) {
    AliErrorStream() 
      << "Directory " << outputDataPath.Data() << " already exists" << endl;
    return kFALSE;
  }
  else {
    AliDebugStream(2) << "Making directory " <<  outputDataPath.Data() << endl;
    gSystem->mkdir(outputDataPath.Data());
  }  
  
/*
  // std::map implementation

  const map<string, string>& kMap = dataMap->GetMap();
  map<string, string>::const_iterator it;
  
  for ( it = kMap.begin(); it != kMap.end(); it++ ) {
    string path = it->first;
    string data = it->second; 
*/    

  const TMap& kMap = dataMap->GetMap();
  TMapIter it(&kMap);
  TObject* keyObj;
  while ( ( keyObj = it.Next() ) ) {

    TString tpath = ((TObjString*)keyObj)->String();

    TObject* dataObj = kMap.GetValue(keyObj);
    if ( ! dataObj ) {
      AliErrorStream() 
        << "Cannot find value when iterating over map." << endl;
      return kFALSE;
    }      
    TString tdata = ((TObjString*)dataObj)->String();   
  
    string path = tpath.Data();
    string data = tdata.Data(); 
  
    string::size_type slashPos =  path.find_last_of("/");
    if ( slashPos != string::npos ) {
      string dataDir(path, 0, slashPos);
      TString dataDirPath
        = outputDataPath + "/" + dataDir.c_str();
      if ( ! gSystem->OpenDirectory(dataDirPath.Data()) ) {
        AliDebugStream(2) << "Making directory " <<  dataDirPath.Data() << endl;
        gSystem->mkdir(dataDirPath.Data(), kTRUE);
      }
    }     

    TString dataPath
      = outputDataPath + "/" + path.c_str();

    ofstream out(dataPath.Data(), ios::out);
    if ( ! out.good() ) {
      AliErrorStream() 
        << "Cannot open output file  " << outputDataPath.Data() << endl;
   
      return kFALSE;  
    }
    
    out << data;
    out.close();
  }

  delete dataMap;
  return kTRUE;
}
  
//_____________________________________________________________________________
Bool_t AliMpDataProcessor::GenerateCode(AliMpDataMap* dataMap) 
{
/// Generate C++ code from dataMap.
/// <pre>
/// AliMpDataProcessor mp;
/// AliMpDataMap* dataMap = mp.CreateDataMap();
/// mp.GenerateCode(dataMap);
/// </pre>
/// Not really used, but kept for eventual future explorations.

/*
  // std::map implementation

  const map<string, string>& kMap = dataMap->GetMap();
  map<string, string>::const_iterator it;
  
  for ( it = kMap.begin(); it != kMap.end(); it++ ) {
    string path = it->first;
    string data = it->second; 
*/

  const TMap& kMap = dataMap->GetMap();
  TMapIter it(&kMap);
  TObject* keyObj;
  while ( ( keyObj = it.Next() ) ) {

    TString tpath = ((TObjString*)keyObj)->String();

    TObject* dataObj = kMap.GetValue(keyObj);
    if ( ! dataObj ) {
      AliErrorStream() 
        << "Cannot find value when iterating over map." << endl;
      return kFALSE;
    }      
    TString tdata = ((TObjString*)dataObj)->String();   
  
    GenerateFunction(tpath, tdata);
  }
  
  GenerateFill();

  return kTRUE;
}
  
  
