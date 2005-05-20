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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// access class to a DB file inside an organized directory structure         //
// (DBFolder/detector/dbType/detSpecType)                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRegexp.h>

#include "AliLog.h"
#include "AliRunData.h"
#include "AliSelectionMetaData.h"
#include "AliRunDataOrganizedFile.h"


ClassImp(AliRunDataOrganizedFile)

//_____________________________________________________________________________
AliRunDataOrganizedFile::AliRunDataOrganizedFile(const char* DBFolder) :
  AliRunDataStorage(),
  fDBFolder(DBFolder)
{
// constructor
  TString buffer(fDBFolder);
  gSystem->ExpandPathName(buffer);
  if(!gSystem->OpenDirectory(buffer)){
    AliError(Form("Path %s not a directory",fDBFolder.Data()));
  }
}

//_____________________________________________________________________________
AliRunDataOrganizedFile::~AliRunDataOrganizedFile()
{
// destructor

  if (fDBFolder) fDBFolder="";
}

//_____________________________________________________________________________
AliRunDataOrganizedFile::AliRunDataOrganizedFile(const AliRunDataOrganizedFile& /*db*/) :
  AliRunDataStorage(),
  fDBFolder("")
{
// copy constructor

  AliFatal("not implemented");
}

//_____________________________________________________________________________
AliRunDataOrganizedFile& AliRunDataOrganizedFile::operator = (const AliRunDataOrganizedFile& /*db*/)
{
// assignment operator

 AliFatal("not implemented");
 return *this;
}

//_____________________________________________________________________________
AliRunData* AliRunDataOrganizedFile::GetEntry(AliSelectionMetaData& selMetaData, Int_t runNumber)
{
// get an object from the data base   

 TDirectory* saveDir = gDirectory;

// Find the right file in the directory

  TObjArray *objarr = FindDataBaseFile(selMetaData, runNumber);
  if(!objarr || objarr->GetEntries()==0) return NULL;
  if(objarr->GetEntries()>1) 
  AliWarning("Warning: more than 1 file match requirements, I will open the first found!");
   
// Open the file

 TObjString *objstr= (TObjString*) objarr->At(0); // there should be only one item
 TString fileName(objstr->GetName());
 TFile *dbFile = new TFile(fileName.Data(),"READ");
 if (!dbFile || !dbFile->IsOpen()) {
   AliError(Form("could not open file %s", fileName.Data()));
   return NULL;
 }
   
// get the only AliRunData object from the file
// I suppose that the object in the file is a AliRunData entry with
// name="DetSpecType" (set in SelectionMetaData)

 dbFile->cd();

 AliRunData *entry = (AliRunData*) dbFile->Get(selMetaData.GetDetSpecType());
 
 if(!entry || !entry->InheritsFrom(AliRunData::Class())) {
   AliError(Form("No entry named %s found!",selMetaData.GetDetSpecType())); 
   dbFile->Close(); delete dbFile;  
   if (saveDir) saveDir->cd(); else gROOT->cd();
   return NULL;
 }
  
// now set the run range and version got from the filename 
// to the object's metadata!  
  
 TString fileNameShort=fileName;
 fileNameShort.Remove(0,fileNameShort.Last('/')+1);
 int numbers[3]={-1,-1,-1}; // numbers[0]=firstRun, numbers[1]=lastRun, numbers[2]=Version
 GetNumbers(fileNameShort, numbers);
 entry->SetRunRange(numbers[0],numbers[1]);
 entry->SetVersion(numbers[2]);
  
// close file, return retieved entry

 dbFile->Close(); delete dbFile;  
 if (saveDir) saveDir->cd(); else gROOT->cd();
   
 if(selMetaData.GetVersion() > -1 && numbers[2] != selMetaData.GetVersion()) 
     AliWarning(Form("Warning: selected version (%d) not found, got version %d instead",
            selMetaData.GetVersion(),numbers[2]));
 return entry;
  
}


//_____________________________________________________________________________
Bool_t AliRunDataOrganizedFile::PutEntry(AliRunData* entry)
{
// puts an object into the database

// AliRunData entry is composed by the object and its MetaData
// this method takes the metaData, reads the name, runRange and Version
// creates the directory structure and the file name
// looks for runs with same or overlapping runrange, if exist increment version
// (therefore version should not be put in the metadata)
// if the runrange is different (but overlapping) from a preceding version, a warning message
// is issued. 
// sets the runrange and version in the object metadata = -1 (to avoid inconsistencies)
// open the filem, write the entry in the file.
// Note: the key name of the entry is "DetSpecType"
// return result 
   
 if(!entry) return kFALSE;
 TDirectory* saveDir = gDirectory;
  
 Int_t firstRun=entry->GetObjectMetaData().GetFirstRun();
 Int_t lastRun=entry->GetObjectMetaData().GetLastRun();
 if(firstRun<0 || lastRun<0 || lastRun<firstRun) {
   AliError(Form("Run range not set or not valid: %d - %d !", firstRun, lastRun));
   return kFALSE;
 }

 TString name(entry->GetObjectMetaData().GetName()); 
  
 while(name.EndsWith("/")) name.Remove(name.Last('/')); 
 while(name.BeginsWith("/")) name.Remove(name.First('/'),1);

 TString detSpecType(name(name.Last('/')+1, name.Length()-name.Last('/')));
 
 TString buffer(fDBFolder);
 gSystem->ExpandPathName(buffer);
 while(buffer.EndsWith("/")) buffer.Remove(buffer.Last('/'));
 
 void *dir=0;
 Int_t index = -1;
 name+='/'; // name=detector/dbType/detSpecType/

 while ((index = name.Index("/")) >= 0) {
   TString dirName(name(0, index));
   buffer+='/'; buffer+=dirName;
   dir=gSystem->OpenDirectory(buffer);
   if (!dir) {
    AliWarning(Form("Directory %s does not exist! It will be created...",buffer.Data()));
    TString command = "mkdir "+ buffer;
    gSystem->Exec(command.Data());
   }
   name.Remove(0, index+1);
 } 
 
 TString strfName="Run";
 
 TString levelContent="";
 Int_t maxVersion=-1, run1=-1, run2=-1;
 int numbers[3]={-1,-1,-1}; // numbers[0]=firstRun, numbers[1]=lastRun, numbers[2]=Version
 while(levelContent=gSystem->GetDirEntry(dir)){
   if(levelContent=="." || levelContent=="..") continue; if(levelContent=="") break;
   if(levelContent.Contains(strfName)){
     GetNumbers(levelContent, numbers);
     if((firstRun>=numbers[0] && firstRun<=numbers[1]) ||
	(lastRun>=numbers[0] && lastRun<=numbers[1]) ||
	(firstRun<=numbers[0] && lastRun>=numbers[1])) {// overlap!     
        if(numbers[2]>maxVersion) {
	   maxVersion=numbers[2];
	   run1=numbers[0]; run2=numbers[1];
	}
     }
   }
 }
 
 if((run1!=-1 && run2!=-1) && (firstRun!=run1 || lastRun!=run2)) 
    AliWarning(Form("Run range modified w.r.t. preceding version (%d, %d)",run1, run2));
    
 
 if(firstRun==lastRun) {
   strfName+=firstRun;
 }else{
   strfName+=firstRun; strfName+="-"; strfName+=lastRun;
 }
 strfName+="_v"; strfName+=maxVersion+1; strfName+=".root";
 buffer+='/'; buffer+=strfName;
 
 // opening file
 TFile *dbFile = new TFile(buffer.Data(),"NEW");
 if(!dbFile || !dbFile->IsWritable()){
   AliError(Form("The data base file is not writable. "
		  "The object %s was not inserted", entry->GetName()));
   return kFALSE;
 }
 
 dbFile->cd();

 entry->SetRunRange(-1,-1); entry->SetVersion(-1);
 
 // write object
  Bool_t result = (entry->Write(detSpecType) != 0);
  if (saveDir) saveDir->cd(); else gROOT->cd();
  dbFile->Close(); delete dbFile;
  
  if(result) {
    AliInfo(Form("Run object %s",entry->GetName()));
    AliInfo(Form("was successfully written into file %s !",buffer.Data()));
  }
  
  return result;
}

/*****************************************************************************/ 

TObjArray* AliRunDataOrganizedFile::FindDataBaseFile(AliSelectionMetaData& selMetaData, Int_t runNumber){
// Find DataBase file name in a local directory. The filename must be in the form: Run#run1-#run2_v#version.root
// TRegexp allowed: selMetaData's name can be for example: "detector/*" !!

 TObjArray *fileNameColl=new TObjArray();

 TString buffer(fDBFolder);
 if(!(buffer.EndsWith("/"))) buffer+="/";
 gSystem->ExpandPathName(buffer);
  
 TString bufftInit=buffer; // buffInit="$ALICE_ROOT/DB/
 TString levelContent="";
 
 TString detector(selMetaData.GetDetector());
 TString dbType(selMetaData.GetDBType());
 TString detSpecType(selMetaData.GetDetSpecType());
 int selVersion = selMetaData.GetVersion();

 void *dirLevInit = gSystem->OpenDirectory(buffer);
 while(levelContent=gSystem->GetDirEntry(dirLevInit)){ // lev0! In Detector directory (ZDC, TPC...)!!

   if(levelContent=="." || levelContent=="..") continue; if(levelContent=="") break;
   if(!(detector=="*") && !levelContent.Contains(TRegexp(detector)) ) continue;

   buffer=bufftInit+levelContent; buffer+='/'; // buffer="$ALICE_ROOT/DB/detector/
   TString bufft0=buffer; // bufft0="$ALICE_ROOT/DB/detector/

   void *dirLev0 = gSystem->OpenDirectory(buffer);
   while(levelContent=gSystem->GetDirEntry(dirLev0)){ // lev1! dbType directory (Calib, Align)!!

     if(levelContent=="." || levelContent=="..") continue; if(levelContent=="") break;
     if(!(dbType=="*") && !levelContent.Contains(TRegexp(dbType))) continue; 

     buffer=bufft0+levelContent;buffer+='/'; // buffer="$ALICE_ROOT/DB/detector/dbType/
     TString bufft1=buffer; // bufft1="$ALICE_ROOT/DB/detector/dbType/

     void *dirLev1 = gSystem->OpenDirectory(buffer);
     while(levelContent=gSystem->GetDirEntry(dirLev1)){ // lev2! detSpecType directory (Pedestals, gain....)!!

       if(levelContent=="." || levelContent=="..") continue; if(levelContent=="") break;
       if(!(detSpecType=="*") && !levelContent.Contains(TRegexp(detSpecType))) continue; 

       buffer=bufft1+levelContent;buffer+='/'; // buffer="$ALICE_ROOT/DB/detector/dbType/detSpecType/
       TString bufft2=buffer; // bufft2="$ALICE_ROOT/DB/detector/dbType/detSpecType/

       void *dirLev2 = gSystem->OpenDirectory(buffer);
       int oldVers=-1;
       TObjString *str=0;
       while(levelContent=gSystem->GetDirEntry(dirLev2)){ // lev3! Run directory (Run#XXX-#YYY_v#ZZ.root)!!

         if(levelContent=="." || levelContent=="..") continue; if(levelContent=="") break;
         if(!levelContent.BeginsWith("Run")) continue;
         
	 int numbers[3]={-1,-1,-1}; // numbers[0]=firstRun, numbers[1]=lastRun, numbers[2]=Version 
	 GetNumbers(levelContent,numbers);
	 if(numbers[0]<0 || numbers[1]<0 || numbers[2]<0 ) continue; // wrong run filename format!
         if(numbers[0]>runNumber ||  numbers[1]<runNumber) continue; // data not valid for run number

	 if((selVersion == -1 || numbers[2] <= selVersion) && numbers[2] >= oldVers) {
	   buffer=bufft2+levelContent;
	   str=new TObjString(buffer.Data());
	   oldVers=numbers[2];
	   if(numbers[2]==selVersion) break;
	   }
         } // end loop on runs       
       if(str) fileNameColl->Add(str);
       } // end loop in lev1
     } // end loop in lev0
   } // end loop in levInit

 AliInfo(Form("Found %d entries matching requirements", fileNameColl->GetEntriesFast()));
 ToAliInfo(fileNameColl->ls());
 return fileNameColl;
}

//_____________________________________________________________________________

void AliRunDataOrganizedFile::GetNumbers(const TString strName, int *numArray)
{
// Gets the numbers (#Run1, #Run2, #Version) from the filename: Run#Run1-#Run2_v#Version.root 

 int indexMinus=strName.Last('-');
 int indexUScore=strName.Last('_');
 int indexPoint=strName.Last('.');
 
 if(indexUScore<0 || indexPoint<0 )
    {AliError(Form("Check sintax %s",strName.Data())); return;}

 if(indexMinus<0){ // only 1 Run number!
   TString cRun=strName(3,indexUScore-3);
   if(!(cRun.IsDigit()))
      {AliError(Form("%s not a digit! Check sintax %s",cRun.Data(),strName.Data())); return;}
   numArray[0] = (int) strtol(cRun.Data(),0,10);
   numArray[1] = numArray[0];
 }else{
   TString cFirstRun = strName(3,indexMinus-3);	 
   TString cLastRun = strName(indexMinus+1,indexUScore-(indexMinus+1));	 
   if(!(cFirstRun.IsDigit()) || !(cLastRun.IsDigit()))
      {AliError(Form("%s or %s are not digit! Check sintax %s",
         cFirstRun.Data(), cLastRun.Data(), strName.Data())); return;}
   numArray[0] = (int) strtol(cFirstRun.Data(),0,10);
   numArray[1] = (int) strtol(cLastRun.Data(),0,10);
 }
 TString cVersion = strName(indexUScore+2,indexPoint-(indexUScore+2));
 if(!(cVersion.IsDigit())){
   AliError(Form("%s not a digit! Check sintax %s",cVersion.Data(),strName.Data())); return;}
 numArray[2] = (int) strtol(cVersion.Data(),0,10);

 return;
}
