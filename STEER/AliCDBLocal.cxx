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

/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
// access class to a DB file inside an organized directory structure                           //
// file name = "DBFolder/detector/dbType/detSpecType/Run#firstRun-#lastRun _v#version.root"    //                                   
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRegexp.h>

#include "AliLog.h"
#include "AliCDBLocal.h"


ClassImp(AliCDBLocal)

//_____________________________________________________________________________
AliCDBLocal::AliCDBLocal(const char* DBFolder) :
  AliCDBStorage(),
  fDBFolder(DBFolder)
{
// constructor
  gSystem->ExpandPathName(fDBFolder);
  void *dir=0;
  if(!(dir=gSystem->OpenDirectory(fDBFolder))){
    AliFatal(Form("Path %s not a directory",fDBFolder.Data()));
  }
  gSystem->FreeDirectory(dir);
 
  while(fDBFolder.EndsWith("/")) fDBFolder.Remove(fDBFolder.Last('/')); 
  fDBFolder+="/";
}

//_____________________________________________________________________________
AliCDBLocal::~AliCDBLocal()
{
 // destructor

}

//_____________________________________________________________________________
AliCDBLocal::AliCDBLocal(const AliCDBLocal& /*db*/) :
  AliCDBStorage(),
  fDBFolder("")
{
 // copy constructor

 AliFatal("not implemented");
}

//_____________________________________________________________________________
AliCDBLocal& AliCDBLocal::operator = (const AliCDBLocal& /*db*/)
{
// assignment operator

 AliFatal("not implemented");
 return *this;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBLocal::GetEntry(AliCDBMetaDataSelect& selMetaData, Int_t runNumber)
{
// get an object from the data base   

 TDirectory* saveDir = gDirectory;

 // Find the right file in the directory
 TString prefix="_v"; // development mode: fileName=Run#Run1-#Run2_v#Version.root
 if(fStorageMode==kProduction) prefix="_Prod"; // production mode: fileName=Run#Run1-#Run2_Prod#Version.root	 	 

 TString buffer(fDBFolder);
 TString name(selMetaData.GetName());
 buffer+=name; buffer+='/';

 int selVersion = selMetaData.GetVersion();

 void *dir = gSystem->OpenDirectory(buffer);
 if(!dir) {
   AliError(Form("Directory %s not found", name.Data()));
   AliError(Form("in DB folder %s", fDBFolder.Data()));
   return NULL;
 }

 TString fileName="";
 TString levelContent="";
 int oldVers=-1;
 // in this array the "numbers" of the retrieved fileName (Run1, Run2, Version) are stored for later usage
 int fileNumbers[3]={-1,-1,-1}; 
 while(levelContent=gSystem->GetDirEntry(dir)){ 

   if(levelContent=="." || levelContent=="..") continue; if(levelContent=="") break;
   if(!levelContent.Contains("Run") || !levelContent.Contains(prefix)) continue;
 	 
   int numbers[3]={-1,-1,-1}; // numbers[0]=firstRun, numbers[1]=lastRun, numbers[2]=Version 
   // gets the 3 "numbers" in the file name
   if(!DecodeFileName(levelContent,numbers, prefix)) continue; // wrong run filename format
   if(numbers[0]>runNumber ||  numbers[1]<runNumber) continue; // data not valid for run number

   if(selVersion == -1) {
     if(numbers[2] >= oldVers){
        if(numbers[2] == oldVers){ 
	   // more than one file valid for the run -> error!   
           AliError(Form("More than one object valid for run %d, version %d!", runNumber, oldVers));
           AliError(Form("No object will be returned!"));
	   gSystem->FreeDirectory(dir);
	   return NULL;
        }
        fileName=levelContent;
        oldVers=numbers[2];
        fileNumbers[0]=numbers[0]; fileNumbers[1]=numbers[1]; fileNumbers[2]=numbers[2]; 
     }
   } else {
     if(numbers[2] == selVersion){
        if(fileName != ""){ 
	   // filename was already assigned, this means there is more than one file valid for the run -> error!   
           AliError(Form("More than one object valid for run %d, version %d!", runNumber, selVersion));
           AliError(Form("No object will be returned!"));
	   gSystem->FreeDirectory(dir);
	   return NULL;
        }
        fileName=levelContent;
        fileNumbers[0]=numbers[0]; fileNumbers[1]=numbers[1]; fileNumbers[2]=numbers[2]; 
     }
   }   

 } // end loop on runs    
 
 gSystem->FreeDirectory(dir);
 buffer+=fileName;
 if(!buffer.EndsWith(".root")){
    AliError(Form("No DB file matching criteria found!"));
    return NULL;    
 }
 
 TFile *DBFile = new TFile(buffer.Data(),"READ");
 if(!DBFile || !DBFile->IsOpen()) {
    AliError(Form("could not open file %s", buffer.Data()));
    return NULL;
 }
 
 AliInfo(Form("File %s succesfully opened", buffer.Data()));
   
// get the only AliCDBEntry object from the file
// I assume that the object in the file is a AliCDBEntry entry with
// name="detSpecType" (set in CDBMetaDataSelect)

 DBFile->cd();

 AliCDBEntry *entry = (AliCDBEntry*) DBFile->Get(selMetaData.GetDetSpecType());
 
 if(!entry || !entry->InheritsFrom(AliCDBEntry::Class())) {
   AliError(Form("No entry named %s found!",selMetaData.GetDetSpecType())); 
   DBFile->Close(); delete DBFile; DBFile=0; 
   if (saveDir) saveDir->cd(); else gROOT->cd();
   return NULL;
 }
  
// Version 1:
// set the run range and version got from the filename 
// to the object's metadata!  
  
// entry->SetRunRange(fileNumbers[0],fileNumbers[1]);
// entry->SetVersion(fileNumbers[2]);

// Version 2: The object's metadata are not reset during storage
// If object's metadata runRange or version do not match with filename,
// it means that someone renamed file by hand. In this case a warning msg is issued.
 Int_t objFirstRun=(entry->GetCDBMetaData()).GetFirstRun();
 Int_t objLastRun=(entry->GetCDBMetaData()).GetLastRun();
 Int_t objVersion=(entry->GetCDBMetaData()).GetVersion();
 
 if(objFirstRun != fileNumbers[0] || objLastRun != fileNumbers[1] || objVersion != fileNumbers[2]){
    AliWarning(Form("Either RunRange or Version in the object's metadata do noth match with fileName numbers:"));
    AliWarning(Form("someone renamed file by hand!"));
 }

// close file, return retieved entry

 DBFile->Close(); delete DBFile; DBFile=0;
 if (saveDir) saveDir->cd(); else gROOT->cd();
   
//  if(selMetaData.GetVersion() > -1 && fileNumbers[2] != selMetaData.GetVersion()) 
//      AliWarning(Form("Warning: selected version (%d) not found, got version %d instead",
//             selMetaData.GetVersion(),fileNumbers[2]));

 return entry;
  
}


//_____________________________________________________________________________
Bool_t AliCDBLocal::PutEntry(AliCDBEntry* entry)
{
// puts an object into the database

// AliCDBEntry entry is composed by the object and its MetaData
// this method takes the metaData, reads the name, runRange and Version
// creates the directory structure and the file name
// looks for runs with same or overlapping runrange, if exist increment version
// (therefore version should not be put in the metadata)
// if the runrange is different (but overlapping) from a preceding version, a warning message
// is issued. 
// sets the runrange and version in the object metadata = -1 (to avoid inconsistencies)
// open the filem, write the entry in the file.
// Note: the key name of the entry is "detSpecType"
// return result 
   
 if(!entry) return kFALSE;
 TDirectory* saveDir = gDirectory;
  
 Int_t firstRun=entry->GetCDBMetaData().GetFirstRun();
 Int_t lastRun=entry->GetCDBMetaData().GetLastRun();
 if(firstRun<0 || lastRun<0 || lastRun<firstRun) {
    AliError(Form("Run range not set or not valid: %d - %d !", firstRun, lastRun));
    return kFALSE;
 }

 TString name(entry->GetName()); 
 
 
 TString detSpecType(name(name.Last('/')+1, name.Length()-name.Last('/')));
 TString buffer(fDBFolder);
 
 void *dir=0;
 Int_t index = -1;
 name+='/'; // name=detector/dbType/detSpecType/

 while ((index = name.Index("/")) >= 0) {
   TString dirName(name(0, index+1));
   buffer+=dirName;
   dir=gSystem->OpenDirectory(buffer);
   if (!dir) {
     AliWarning(Form("Directory %s does not exist! It will be created...",buffer.Data()));
     gSystem->mkdir(buffer.Data());
   }
   name.Remove(0, index+1);
   gSystem->FreeDirectory(dir);
 } 
 
 dir = gSystem->OpenDirectory(buffer);
 TString levelContent="";
 Int_t maxVersion=-1, run1=-1, run2=-1;
 int numbers[3]={-1,-1,-1}; // numbers[0]=firstRun, numbers[1]=lastRun, numbers[2]=Version
 while(levelContent=gSystem->GetDirEntry(dir)){
   if(levelContent=="." || levelContent=="..") continue; if(levelContent=="") break;
   if(levelContent.Contains("Run")){
     if(levelContent.Contains("_Prod")) continue; //skip "Production" links
     if(!DecodeFileName(levelContent, numbers, "_v")) continue;
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
 gSystem->FreeDirectory(dir);
 
 if((run1!=-1 && run2!=-1) && (firstRun!=run1 || lastRun!=run2)) 
    AliWarning(Form("Run range modified w.r.t. preceding version (%d, %d)",run1, run2));
    
 TString strfName=EncodeFileName(firstRun, lastRun, maxVersion+1);
 buffer+=strfName;
 
 // opening file
 TFile *DBFile = new TFile(buffer.Data(),"NEW");
 if(!DBFile || !DBFile->IsWritable()){
    AliError(Form("The data base file is not writable. "
		"The object %s was not inserted", entry->GetName()));
    if(!DBFile->IsWritable()) DBFile->Close(); DBFile->Delete(); delete DBFile; DBFile=0;
    return kFALSE;
 }
  
 DBFile->cd();
 
 entry->SetVersion(maxVersion+1);
 
 // write object
 Bool_t result = (entry->Write(detSpecType) != 0); 
 if (saveDir) saveDir->cd(); else gROOT->cd();
 DBFile->Close(); DBFile->Delete(); delete DBFile; DBFile=0;
 if(result) {
    AliInfo(Form("Run object %s",entry->GetName()));
    AliInfo(Form("was successfully written into file %s",buffer.Data()));
 }

 return result;

}

/*****************************************************************************/ 
TObjArray* AliCDBLocal::FindDBFiles(const char* name, Int_t runNumber){
// Find DataBase file name in a local directory. The filename must be in the form: Run#run1-#run2_v#version.root
// TRegexp allowed: name can be for example: "detector/*" !!

 TObjArray *FileNameColl=new TObjArray();

 TString prefix="_v"; // development mode: fileName=Run#Run1-#Run2_v#Version.root
 if(fStorageMode==kProduction) prefix="_Prod"; // production mode: fileName=Run#Run1-#Run2_Prod#Version.root	 	 

 TString buffer(fDBFolder);
// gSystem->ExpandPathName(buffer);
  
 TString bufftInit=buffer; // buffInit="$ALICE_ROOT/DB/
 TString levelContent="";
 
 AliCDBMetaDataSelect selMetaData(name);
 
 TString detector(selMetaData.GetDetector());
 TString dbType(selMetaData.GetDBType());
 TString detSpecType(selMetaData.GetDetSpecType());
 int selVersion = selMetaData.GetVersion();

 void *dirLevInit = gSystem->OpenDirectory(buffer);
 while(levelContent=gSystem->GetDirEntry(dirLevInit)){ // lev0! In detector directory (ZDC, TPC...)!!

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
       TObjString *str=0;
       while(levelContent=gSystem->GetDirEntry(dirLev2)){ // lev3! Run directory (Run#XXX-#YYY_v#ZZ.root)!!

         if(levelContent=="." || levelContent=="..") continue; if(levelContent=="") break;
         if(!levelContent.BeginsWith("Run")) continue;
       	 if(!levelContent.Contains(prefix)) continue;
	 	 
	 int numbers[3]={-1,-1,-1}; // numbers[0]=firstRun, numbers[1]=lastRun, numbers[2]=Version 
	 if(!DecodeFileName(levelContent,numbers, prefix)) continue; // wrong run filename format!
         if(numbers[0]>runNumber ||  numbers[1]<runNumber) continue; // data not valid for run number

	 if(numbers[2]==selVersion) {
	    buffer=bufft2+levelContent;
	    str=new TObjString(buffer.Data());
            FileNameColl->Add(str);
	    break;
	 }
	 if(selVersion == -1) { // if version is not specified, collect all versions
	    buffer=bufft2+levelContent;
	    str=new TObjString(buffer.Data());
            FileNameColl->Add(str);
	   }
         } // end loop on runs       
       } // end loop in lev1
     } // end loop in lev0
   } // end loop in levInit

 AliInfo(Form("Found %d entries matching requirements", FileNameColl->GetEntriesFast()));
 ToAliInfo(FileNameColl->ls());
 return FileNameColl;
}

//_____________________________________________________________________________
void AliCDBLocal::TagForProduction(const AliCDBMetaDataSelect& selMetaData, UInt_t prodVers){

TString workingDir=gSystem->pwd();
//Build the file path 
TString buffer(fDBFolder); //gSystem->ExpandPathName() already done in ctor
TString fName="";

buffer+=selMetaData.GetName(); buffer+='/';
//gSystem->ExpandPathName(dirName);

if(!gSystem->cd(buffer))
   {AliError(Form("Directory %s does not exist... check name!", buffer.Data())); gSystem->cd(workingDir.Data()); return;}

// if version is not specified (=-1), then tag the highest version (BE CAREFUL!!)
if(selMetaData.GetVersion() != -1){  
  //Build the filename
   fName = EncodeFileName(selMetaData.GetFirstRun(), selMetaData.GetLastRun(), selMetaData.GetVersion());
} else {
   //look in directory for valid DB files, seek highest version
   void *dir = gSystem->OpenDirectory(buffer);
   TString levelContent="";
   int oldVers=-1;
   while(levelContent=gSystem->GetDirEntry(dir)){ 

      if(levelContent=="." || levelContent=="..") continue; if(levelContent=="") break;
      if(!levelContent.Contains("Run") || !levelContent.Contains("_v")) continue;
 	 
      int numbers[3]={-1,-1,-1}; // numbers[0]=firstRun, numbers[1]=lastRun, numbers[2]=Version 
      // gets the 3 "numbers" in the file name
      if(!DecodeFileName(levelContent,numbers, "_v")) continue; // wrong run filename format!
      if(numbers[0] != selMetaData.GetFirstRun() ||  numbers[1] != selMetaData.GetLastRun()) continue;
      if(numbers[2] >= oldVers) {
         fName=levelContent;
         oldVers=numbers[2];
      }   
   } // end loop on runs    
}

   //check that the flename exists
if(!gSystem->IsFileInIncludePath(fName.Data())){
      AliError(Form("File name %s not found... check!", fName.Data())); 
      gSystem->cd(workingDir.Data()); 
      return;
}

// file exists: make symbolic link! 
TString prodfName=EncodeFileName(selMetaData.GetFirstRun(), selMetaData.GetLastRun(), prodVers, "_Prod");
if(gSystem->Symlink(fName.Data(),prodfName.Data())<0){
    AliError(Form("Link name already existing (%s): linkage failed!",prodfName.Data()));
} else {
    AliError(Form("File %s tagged for production with symlink %s",fName.Data(), prodfName.Data()));
}

gSystem->cd(workingDir);
return; 

}

//_____________________________________________________________________________
Bool_t AliCDBLocal::DecodeFileName(const TString strName, int *numArray, TString prefix)
{
// Gets the numbers (#Run1, #Run2, #Version) 
// from the filename: Run#Run1-#Run2_v#Version.root or  Run#Run1-#Run2_Prod#prodVers.root 

 int indexMinus=strName.Last('-');
 int indexUScore=strName.Last('_');
 int indexPoint=strName.Last('.');
 
 int nSkipChar=prefix.Length(); // prefix can be either "_v" or "_Prod" depending on fStorageMode 
 //if(prefix=="_v") {nSkipChar=2;} // development mode: _v# skip 2 characters
 //else if(prefix=="_Prod") {nSkipChar=5;} // production mode: _Prod# skip 5 characters

 if(indexUScore<0 || indexPoint<0 )
    {AliDebug(2, Form("Check sintax %s",strName.Data())); return kFALSE;}

 if(indexMinus<0){ // only 1 Run number!
   TString cRun=strName(3,indexUScore-3);
   if(!(cRun.IsDigit()))
      {AliDebug(2, Form("%s not a digit! Check sintax %s",cRun.Data(),strName.Data())); return kFALSE;}
   numArray[0] = (int) strtol(cRun.Data(),0,10);
   numArray[1] = numArray[0];
 }else{
   TString cFirstRun = strName(3,indexMinus-3);	 
   TString cLastRun = strName(indexMinus+1,indexUScore-(indexMinus+1));	 
   if(!(cFirstRun.IsDigit()) || !(cLastRun.IsDigit()))
      {AliDebug(2, Form("%s or %s are not digit! Check sintax %s",
         cFirstRun.Data(), cLastRun.Data(), strName.Data())); return kFALSE;}
   numArray[0] = (int) strtol(cFirstRun.Data(),0,10);
   numArray[1] = (int) strtol(cLastRun.Data(),0,10);
 }
// TString cVersion = strName(indexUScore+2,indexPoint-(indexUScore+2));
 TString cVersion = strName(indexUScore+nSkipChar,indexPoint-(indexUScore+nSkipChar));
 if(!(cVersion.IsDigit())){
   AliDebug(2, Form("%s not a digit! Check sintax %s",cVersion.Data(),strName.Data())); return kFALSE;}
 numArray[2] = (int) strtol(cVersion.Data(),0,10);

 return kTRUE;
}


//_____________________________________________________________________________
TString AliCDBLocal::EncodeFileName(int firstRun, int lastRun, int version, TString prefix){
// Builds a file name of the form: Run#firstRun-#lastRun_v#Version.root

TString fName="Run"; 
if(firstRun==lastRun) {
   fName+=firstRun;
}else{
   fName+=firstRun; fName+="-"; fName+=lastRun;
}
fName+=prefix; fName+=version; fName+=".root";

return fName;
}
