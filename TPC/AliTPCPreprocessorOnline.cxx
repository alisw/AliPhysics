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


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Preprocessor class for HLT and DAQ                                        //
//  Possible usage: preprocess TPC calibration data                           //
//  to form needed for online viewing/visualizing TPC calibration data        //
//                                                                            //
//  On HLT or DAQ AliTPCPreprocessorOnline::AddComponent(...) is called,      //
//  whereas ... is either an AliTPCCalPad, or an AliTPCCalROC.                //
//  Their names have to be acording the naming convention for AliTPCCalPads   //
//  Special for CalROCs: Add '_ROC<ROC_Number>' to the name.                  //
//                                                                            //
//  Internal the AliTPCCalPads are stored in a TMap, they are retrieved by    //
//  their name.                                                               //
//                                                                            //
//  Once enough values are accumulated, call ::DumpToFile(fileName).          //
//  A TObjArray is created out of the TMap, which is passed to                //
//  AliTPCCalibViewer::MakeTree(...) and the calibratioTree is written        //
//  is written to "filenName". 
//                               
//  The data flow is as follows: 
/* 
Begin_Html 
   <img src="../TPC/doc/dataFlowCalibViewer.png"> 
End_Html
*/
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

//     
// ROOT includes 
//
#include <TObject.h>
#include <iostream>
#include <TString.h>
#include "TMap.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TIterator.h"


//
// AliRoot includes
//
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalibViewer.h"
#include "AliTPCPreprocessorOnline.h"
#include <fstream>



ClassImp(AliTPCPreprocessorOnline)

AliTPCPreprocessorOnline::AliTPCPreprocessorOnline():TObject(), fMap(0)
{
  //
  // Default constructor
  //
  fMap = new TMap();

}

//_____________________________________________________________________________
AliTPCPreprocessorOnline::AliTPCPreprocessorOnline(const AliTPCPreprocessorOnline &c):TObject(c), fMap(0)
{
  //
  // dummy AliTPCPreprocessorOnline copy constructor
  // not yet working!!!
  //
  fMap = c.fMap;
  printf("AliTPCPreprocessorOnline's copy constructor called, NOT WORKING!!!\n");
}

//_____________________________________________________________________________
AliTPCPreprocessorOnline::AliTPCPreprocessorOnline(TMap *map):TObject(), fMap(0)
{
  //
  // Constructor to "copy" the AliTPCPreprocessorOnline
  //
  fMap = map;
}


//____________________________________________________________________________
AliTPCPreprocessorOnline & AliTPCPreprocessorOnline::operator =(const AliTPCPreprocessorOnline & param){
   //
   // assignment operator - dummy
   // not yet working!!!
   //
  if (this == &param) return (*this);

   fMap = param.fMap;
   std::cout << "AliTPCPreprocessorOnline's assignment operator called, NOT WORKING!!!" << std::endl;
   return (*this);
}


//_____________________________________________________________________________
AliTPCPreprocessorOnline::~AliTPCPreprocessorOnline()
{
   //
   // AliTPCPreprocessorOnline destructor
   //
   printf("AliTPCPreprocessorOnline's destructor called. \n");
   fMap->DeleteValues();
   delete fMap;
}



void AliTPCPreprocessorOnline::AddComponent(TObject *obj) {
   //
   // Add components form HLT or DAQ here 
   // The components are either AliTPCCalPads or AliTPCCalROCs
   // other to be defined
   // 
   // As from HLT they will come ROC wise, they have to be set together to one AliTPCCalPad here
   // Then they are added to the Calibration Tree
   // This calibration tree will be written by AliTPCCalibViewer::MakeTree, once you call ::DumpToFile(fileName)
   // 
   // To distinguish, what kind of component is added, there is a Naming-Convention:
   // The normal, already definded naming-Convention for AliTPCCalPads <PadName>
   // <padName>_ROC<ROC_Number>  for example: "CEQmean_ROC34"
   // 
   // 
   // Get name of obj
   // Check, if it ends with "_ROC<Int_t>" / check, if it contains "ROC"
   // If it contains "ROC", find out, to which calibration CalPad it belongs -> Parse first part of the name
   // Get the corsponding AliTPCCalPad / Make a new AliTPCCalPad 
   // Set "_ROC<Int_t>" to zero /? delete it -> Set the new values, out of obj
   // 

   
   // printf(" ****** AliTPCPreprocessorOnline::AddComponent ****** \n");
  if (!obj) return;
   TString objName = obj->GetName();
   
   // Add a whole AliTPCCalPad to the list
   if (TString(obj->ClassName()) == "AliTPCCalPad") {
      // printf("AliTPCCalPad found\n");
      AliTPCCalPad *calPad = (AliTPCCalPad*)obj;
      TObject *listObj = fMap->GetValue(calPad->GetName());
      if (listObj == 0) {
         // current data not yet written to the list, write it to the list
         fMap->Add(new TObjString(calPad->GetName()), calPad);
         return;
      }
      // current data already present
      if (TString(listObj->ClassName()) != "AliTPCCalPad") 
         Error("AddComponent", "Mismatch in internal list: '%s' is no AliTPCCalPad. \n", listObj->ClassName());
      AliTPCCalPad *listCalPad = (AliTPCCalPad*)listObj;
      listCalPad->Add(listCalPad, -1); // reset the current CalPad
      listCalPad->Add(calPad);
      return;
   }
   
   if (TString(obj->ClassName()) == "AliTPCCalROC") {
      // printf("AliTPCCalROC found\n");
      if (! objName.Contains("_ROC"))
         Warning("AddComponent", "Uncomplete object name: '%s' is missing _ROC<ROC_number>\n", objName.Data());
      TObjArray *stokens = objName.Tokenize("_ROC");
      TString rocNumber(""); 
      if (stokens->GetEntriesFast()>1) rocNumber = ((TObjString*)stokens->At(1))->GetString();
      delete stokens;
      Int_t iroc = -1;
      if (rocNumber.IsAlnum()) iroc = rocNumber.Atoi();
      
      // Extract CalPad Name:
      TString removeString("_ROC");
      removeString += rocNumber;
      objName.ReplaceAll(removeString, "");  // Remove "_ROC<number>" from the end
      
      // objName is cleaned and can now be used to extract the coresponding CalPad from the list
      TObject *listObj = fMap->GetValue(objName.Data());
      AliTPCCalROC *calROC = (AliTPCCalROC*)obj;
      if (iroc == -1) iroc = calROC->GetSector();
      if (iroc != (Int_t)calROC->GetSector()) 
         Warning("AddComponent","Mismatch in  ROC_number: ROC has number %i, in its name %i was specified. Treating now as %i. \n", calROC->GetSector(), iroc, iroc);
      calROC->SetNameTitle(objName.Data(), objName.Data());
      if (listObj == 0) {
         // current data not yet written to the list, write it to the list
         AliTPCCalPad *calPad = new AliTPCCalPad(objName.Data(), objName.Data());
         calPad->SetCalROC(calROC, iroc);
         fMap->Add(new TObjString(objName.Data()), calPad);
         return;
      }
      // current data already present
      if (TString(listObj->ClassName()) != "AliTPCCalPad") 
         Error("AddComponent", "Mismatch in internal list: '%s' is no AliTPCCalPad \n", listObj->ClassName());
      AliTPCCalPad *listCalPad = (AliTPCCalPad*)listObj;
      listCalPad->SetCalROC(calROC, iroc);
      return;
     
      
   }
   
   Warning("AddComponent", "Unknown calibration object '%s', skipped.\n", obj->ClassName());
   
   return;
   
   
}



void AliTPCPreprocessorOnline::DumpToFile(const char* fileName){
   // 
   // Function to dump the tree to file
   // A TObjArray is created out of the TMap, which is passed to   
   // AliTPCCalibViewer::MakeTree(...)  
   // 
   
   printf("AliTPCPreprocessorOnline::DumpToFile\n");
   TObjArray *listOfCalPads = new TObjArray();
   TIterator *iterator = fMap->MakeIterator();
   AliTPCCalPad *calPad = 0x0;
//    while ((calibTracks = (AliTPCcalibTracks*)listIterator->Next()) ){
   TObject * obj=0;
   //   while (( calPad = (AliTPCCalPad*)fMap->GetValue(iterator->Next()) )) {
   while ( ( obj = iterator->Next()) ) {
     calPad = (AliTPCCalPad*)fMap->GetValue(obj);
      if (!calPad) continue;
      calPad = (AliTPCCalPad*)calPad;
      printf("adding the following element to the TObjList: %s \n", calPad->GetName());
      printf("It's type:: %s \n", calPad->ClassName());
      calPad->Print();
      listOfCalPads->Add(calPad);
   }
   printf("writing the tree... \n");
   //   AliTPCCalibViewer::MakeTree(fileName, listOfCalPads, "$ALICE_ROOT/TPC/Calib/MapCalibrationObjects.root");
   AliTPCCalibViewer::MakeTree(fileName, listOfCalPads, 0);
 
}

