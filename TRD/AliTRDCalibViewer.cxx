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

/* $Id: AliTRDCalibViewer.cxx 40390 2010-04-14 09:43:23Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Class which implements AliBaseCalibViewer for the TRD                    //
//  used for the calibration monitor                                         //
//                                                                           //
//  Authors:     Marian Ivanov (Marian.Ivanov@cern.ch)                       //
//               Jens Wiechula (Jens.Wiechula@cern.ch)                       //
//               Ionut Arsene  (iarsene@cern.ch)                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TString.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h> 
#include <TH1F.h>
#include <TMath.h>
#include <TVectorD.h>
#include <THashTable.h>
#include <TObjString.h>
#include <TTimeStamp.h>
#include <TObjString.h>
#include <TTreeStream.h>
#include <TFile.h>
#include <TKey.h>
#include <TGraph.h>
#include <TDirectory.h>
#include <TFriendElement.h>
#include <TGrid.h>
#include <TGeoManager.h>
#include "AliTRDCalDet.h"
#include "AliTRDCalPad.h"
#include "AliTRDCalROC.h"
#include "AliTRDCalChamberStatus.h"
#include "AliTRDCalSingleChamberStatus.h"
#include "AliTRDCalPadStatus.h"
#include "AliTRDCalDCS.h"
#include "AliTRDCalDCSFEE.h"
#include "AliTRDcalibDB.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliTRDalignment.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

#include "AliTRDCalibViewer.h"


using namespace std;

ClassImp(AliTRDCalibViewer)

//_____________________________________________________________________________
AliTRDCalibViewer::AliTRDCalibViewer()
                  :AliBaseCalibViewer()
{
  //
  // Default constructor (just call base class constructor)
  //
}

//_____________________________________________________________________________
AliTRDCalibViewer::AliTRDCalibViewer(const AliTRDCalibViewer &c)
                  :AliBaseCalibViewer(c)
{
  //
  // copy constructor (just call base class copy constructor)
  // 
}

//_____________________________________________________________________________
AliTRDCalibViewer::AliTRDCalibViewer(TTree* tree)
                  :AliBaseCalibViewer(tree)
{
  //
  // Constructor (just call the corresponding base constructor)
  //
}

//_____________________________________________________________________________
AliTRDCalibViewer::AliTRDCalibViewer(const char* fileName, const char* treeName)
                  :AliBaseCalibViewer(fileName, treeName)                   
{
   //
   // Constructor (just call the corresponding base constructor)
   //
}

//_____________________________________________________________________________
AliTRDCalibViewer& AliTRDCalibViewer::operator=(const AliTRDCalibViewer& param) {
  //
  // assignment operator
  //
  fTree = param.fTree;
  fTreeMustBeDeleted = param.fTreeMustBeDeleted;
  fListOfObjectsToBeDeleted = param.fListOfObjectsToBeDeleted;
  fAbbreviation = param.fAbbreviation;
  fAppendString = param.fAppendString;
  return(*this);
}

//_____________________________________________________________________________
AliTRDCalibViewer::~AliTRDCalibViewer()
{
   //
   // AliTRDCalibViewer destructor
   // do nothing, the base class destructor will do the job
}

/*
//_____________________________________________________________________________
void AliTRDCalibViewer::GetTimeInfoOCDB(const Char_t* runList, const Char_t* outFile,
                        	        Int_t firstRun, Int_t lastRun, UInt_t infoFlags,
                                        const Char_t* ocdbStorage) {
//
//  Get time information from OCDB by calling the DumpOCDBtoTree.C macro
//
  DumpOCDBtoTree(runList, outFile, firstRun, lastRun,
                 TESTBIT(infoFlags,1), TESTBIT(infoFlags,2), 
                 TESTBIT(infoFlags,3), TESTBIT(infoFlags,4),
		 ocdbStorage);
}
*/

//_____________________________________________________________________________
const char* AliTRDCalibViewer::AddAbbreviations(char* c, Bool_t printDrawCommand){ 
   // Replace all "<variable>" with "<variable><fAbbreviation>" (Adds forgotten "~")
   // but take care on the statistical information, like "CEQmean_Mean"
   // and also take care on correct given variables, like "CEQmean~"
   // 
   // For each variable out of "listOfVariables":
   // - 'Save' correct items:
   //   - form <replaceString>, take <variable>'s first char, add <removeString>, add rest of <variable>, e.g. "C!#EQmean" (<removeString> = "!#")
   //   - For each statistical information in "listOfNormalizationVariables":
   //     - ReplaceAll <variable><statistical_Information> with <replaceString><statistical_Information>
   //   - ReplaceAll <variable><abbreviation> with <replaceString><abbreviation>, e.g. "CEQmean~" -> "C!#EQmean~"
   //   - ReplaceAll <variable><appendStr> with <replaceString><appendStr>, e.g. "CEQmean.fElements" -> "C!#EQmean.fElements"
   //
   // - Do actual replacing:
   //   - ReplaceAll <variable> with <variable><fAbbreviation>, e.g. "CEQmean" -> "CEQmean~"
   //
   // - Undo saving:
   //   - For each statistical information in "listOfNormalizationVariables":
   //     - ReplaceAll <replaceString><statistical_Information> with <variable><statistical_Information> 
   //   - ReplaceAll <replaceString><abbreviation> with <variable><abbreviation>, e.g. "C!#EQmean~" -> "CEQmean~"
   //   - ReplaceAll <replaceString><appendStr> with <variable><appendStr>, e.g. "C!#EQmean.fElements" -> "CEQmean.fElements"
   // 
   // Now all the missing "~" should be added.
   
   TString str(c);
   TString removeString = "!#";  // very unprobable combination of chars
   TString replaceString = "";
   TString searchString = "";
   TString normString = "";
   TObjArray *listOfVariables = GetListOfVariables();
   // variables used for mapping the pads, mcms, ...
   listOfVariables->Add(new TObjString("SuperModule"));
   listOfVariables->Add(new TObjString("Layer"));
   listOfVariables->Add(new TObjString("Stack"));
   listOfVariables->Add(new TObjString("Channel"));
   listOfVariables->Add(new TObjString("Row"));
   listOfVariables->Add(new TObjString("Column"));
   listOfVariables->Add(new TObjString("Chamber"));
   listOfVariables->Add(new TObjString("PadSuperRow"));
   listOfVariables->Add(new TObjString("PadSuperColumn"));
   listOfVariables->Add(new TObjString("MCMSuperRow"));
   listOfVariables->Add(new TObjString("MCMSuperColumn"));
   listOfVariables->Add(new TObjString("ROB"));
   listOfVariables->Add(new TObjString("MCM"));
   TObjArray *listOfNormalizationVariables = GetListOfNormalizationVariables();
   Int_t nVariables = listOfVariables->GetEntriesFast();
   Int_t nNorm = listOfNormalizationVariables->GetEntriesFast();
   
   Int_t *varLengths = new Int_t[nVariables];
   for (Int_t i = 0; i < nVariables; i++) {
      varLengths[i] = ((TObjString*)listOfVariables->At(i))->String().Length();
   }
   Int_t *normLengths = new Int_t[nNorm];
   for (Int_t i = 0; i < nNorm; i++) {
      normLengths[i] = ((TObjString*)listOfNormalizationVariables->At(i))->String().Length();
   }
   Int_t *varSort = new Int_t[nVariables];
   TMath::Sort(nVariables, varLengths, varSort, kTRUE);
   Int_t *normSort = new Int_t[nNorm];
   TMath::Sort(nNorm, normLengths, normSort, kTRUE);
      
   for (Int_t ivar = 0; ivar < nVariables; ivar++) {
      // ***** save correct tokens *****
      // first get the next variable:
      searchString = ((TObjString*)listOfVariables->At(varSort[ivar]))->String();
      // form replaceString:
      replaceString = "";
      for (Int_t i = 0; i < searchString.Length(); i++) {
         replaceString.Append(searchString[i]);
         if (i == 0) replaceString.Append(removeString);
      }
      // go through normalization:
      for (Int_t inorm = 0; inorm < nNorm; inorm++) {
	normString = ((TObjString*)listOfNormalizationVariables->At(normSort[inorm]))->String();
	str.ReplaceAll(searchString + normString, replaceString + normString);
	// like: str.ReplaceAll("CEQmean_Mean", "C!EQmean_Mean");
      }
      str.ReplaceAll(searchString + fAbbreviation, replaceString + fAbbreviation);
      // like: str.ReplaceAll("CEQmean~", "C!EQmean~");
      str.ReplaceAll(searchString + fAppendString,    replaceString + fAppendString);
      // like: str.ReplaceAll("CEQmean.fElements", "C!EQmean.fElements");
      
      // ***** add missing extensions *****
      str.ReplaceAll(searchString, replaceString + fAbbreviation);
      // like: str.ReplaceAll("CEQmean", "C!EQmean~");
   }
   
   // ***** undo saving *****
   str.ReplaceAll(removeString, "");
  
   if (printDrawCommand) std::cout << "The string looks now like: " << str.Data() << std::endl;
   delete varSort;
   delete normSort;
   return str.Data();
}

//_____________________________________________________________________________
TObjArray* AliTRDCalibViewer::GetListOfVariables(Bool_t printList) {
  //
  // scan the tree  - produces a list of available variables in the tree
  // printList: print the list to the screen, after the scan is done
  //
  TObjArray* arr = new TObjArray();
  TObjString* str = 0;
  if (!fTree) {
    return 0;
  }
  Int_t nentries = fTree->GetListOfBranches()->GetEntries();
  for (Int_t i = 0; i < nentries; i++) {
    str = new TObjString(fTree->GetListOfBranches()->At(i)->GetName());
    str->String().ReplaceAll("_Mean", "");
    str->String().ReplaceAll("_RMS", "");
    str->String().ReplaceAll("_Median", "");
    str->String().ReplaceAll(".", "");
    str->String().ReplaceAll("_Run", "");
    str->String().ReplaceAll("_SuperModule", "");
    str->String().ReplaceAll("_Chamber", "");
    // add all the variables in the tree to a list
    // exception make variables which are used for mapping, specified in AddAbbreviations()
    // These two functions should be kept synchronized with respect to the mapping variables 
    if (!arr->FindObject(str) && 
	!(str->String() == "run" || 
	  str->String() == "SuperModule" || str->String() == "Layer" || str->String() == "Stack" || 
	  str->String() == "Chamber" ||
	  str->String() == "Channel" || str->String() == "Row" || str->String() == "Column" || 
	  str->String() == "PadSuperRow" || str->String() == "PadSuperColumn" || str->String() == "MCMSuperRow"
	  || str->String() == "MCMSuperColumn" || str->String() == "MCM" || str->String() == "ROB")) {
      arr->Add(str);
    }
  }
  
  // loop over all friends (if there are some) and add them to the list
  if (fTree->GetListOfFriends()) {
    for (Int_t ifriend = 0; ifriend < fTree->GetListOfFriends()->GetEntries(); ifriend++){
      // printf("iterating through friendlist, currently at %i\n", ifriend);
      // printf("working with %s\n", fTree->GetListOfFriends()->At(ifriend)->ClassName());
      if (TString(fTree->GetListOfFriends()->At(ifriend)->ClassName()) != "TFriendElement") continue; // no friendElement found
      TFriendElement *friendElement = (TFriendElement*)fTree->GetListOfFriends()->At(ifriend);
      if (friendElement->GetTree() == 0) continue; // no tree found in friendElement
      // printf("friend found \n");
      for (Int_t i = 0; i < friendElement->GetTree()->GetListOfBranches()->GetEntries(); i++) {
	// printf("iterating through friendelement entries, currently at %i\n", i);
	str = new TObjString(friendElement->GetTree()->GetListOfBranches()->At(i)->GetName());
	str->String().ReplaceAll("_Mean", "");
	str->String().ReplaceAll("_RMS", "");
	str->String().ReplaceAll("_Median", "");
	str->String().ReplaceAll(".", "");
	str->String().ReplaceAll("_Run", "");
	str->String().ReplaceAll("_SuperModule", "");
	str->String().ReplaceAll("_Chamber", "");
	if (!(str->String() == "run" || 
	      str->String() == "SuperModule" || str->String() == "Layer" || str->String() == "Stack" || 
	      str->String() == "Chamber" ||
	      str->String() == "Channel" || str->String() == "Row" || str->String() == "Column" || 
	      str->String() == "PadSuperRow" || str->String() == "PadSuperColumn" || str->String() == "MCMSuperRow"
	      || str->String() == "MCMSuperColumn" || str->String() == "MCM" || str->String() == "ROB")){
	  // insert "<friendName>." at the beginning: (<friendName> is per default "R")
	  str->String().Insert(0, ".");
	  str->String().Insert(0, friendElement->GetName());
	  if (!arr->FindObject(str)) {
	    arr->Add(str);
	  }
	  // printf("added string %s \n", str->String().Data());
	}
      }
    }
  } // if (fTree->GetListOfFriends())
  
  arr->Sort();
    
  if (printList) {
    TIterator* iter = arr->MakeIterator();
    iter->Reset();
    TObjString* currentStr = 0;
    while ( (currentStr = (TObjString*)(iter->Next())) ) {
      std::cout << currentStr->GetString().Data() << std::endl;
    }
    delete iter;
  }
  return arr;
}

TObjArray* AliTRDCalibViewer::GetListOfNormalizationVariables(Bool_t printList) const{
  //
  // produces a list of available variables for normalization in the tree
  // printList: print the list to the screen, after the scan is done
  //
   TObjArray* arr = new TObjArray();
   arr->Add(new TObjString("_Mean_Run"));
   arr->Add(new TObjString("_Mean_SuperModule"));
   arr->Add(new TObjString("_Mean_Chamber"));
   arr->Add(new TObjString("_Median_Run"));
   arr->Add(new TObjString("_Median_SuperModule"));
   arr->Add(new TObjString("_Median_Chamber"));
   
   if (printList) {
     TIterator* iter = arr->MakeIterator();
     iter->Reset();
     TObjString* currentStr = 0;
     while ((currentStr = (TObjString*)(iter->Next()))) {
       std::cout << currentStr->GetString().Data() << std::endl;
     }
     delete iter;
   }
   return arr;
}

void AliTRDCalibViewer::GetLayerSectorStack(TString trdString, Int_t& layerNo, Int_t& sectorNo, Int_t& stackNo) const {
  // Get the layer, sector and stack numbers out of a string
  // encoded with the following format:
  // Layer%dSector%dStack%d

  sscanf(trdString.Data(), "Layer%dSector%dStack%d", &layerNo, &sectorNo, &stackNo);

  return;
}

//_____________________________________________________________________________
Int_t AliTRDCalibViewer::EasyDraw(const char* drawCommand, const char* sector, const char* cuts, const char* drawOptions, Bool_t writeDrawCommand) const {
  //
  // easy drawing of data, use '~' for abbreviation of '.fElements'
  // example: EasyDraw("CETmean~-CETmean_mean", "A", "(CETmean~-CETmean_mean)>0")
  // sector: sector-number - only the specified sector will be drwawn
  //         'A'/'C' or 'a'/'c' - side A/C will be drawn
  //         'ALL' - whole TPC will be drawn, projected on one side
  // cuts: specifies cuts
  // drawOptions: draw options like 'same'
  // writeDrawCommand: write the command, that is passed to TTree::Draw
  //

  TString drawStr(drawCommand);

  TString sectorStr(sector);
  Int_t layerNo = -1; 
  Int_t sectorNo = -1; 
  Int_t stackNo = -1;
  GetLayerSectorStack(sectorStr, layerNo, sectorNo, stackNo);
  if(layerNo==-1) {
     Warning("EasyDraw", "The sector string must always contain the Layer number!");
     return -1;
   }
   if(layerNo<0 || layerNo>5) {
     Warning("EasyDraw", "The Layer number must be in the range [0,5] !");
     return -1;
   }
   if(sectorNo!=-1 && (sectorNo<0 || sectorNo>17)) {
     Warning("EasyDraw", "The SuperModule number must be in the range [0,17] !");
     return -1;
   }
   if(stackNo!=-1 && (stackNo<0 || stackNo>4)) {
     Warning("EasyDraw", "The Stack number must be in the range [0,4] !");
     return -1;
   }

   TString cutStr("");

   Bool_t dangerousToDraw = drawStr.Contains(":") || drawStr.Contains(">>");
   if (dangerousToDraw) {
      Warning("EasyDraw", "The draw string must not contain ':' or '>>'. Using only first variable for drawing!");
      drawStr.Resize(drawStr.First(":"));
   }

   TString drawOptionsStr("");
   TRandom rnd(0);
   Int_t rndNumber = rnd.Integer(10000);

   if (drawOptions && strcmp(drawOptions, "") != 0)
      drawOptionsStr += drawOptions;
   else
      drawOptionsStr += "profcolz";

   const Int_t gkNRows[ 5] = {16, 16, 12, 16, 16};  // number of pad rows in the chambers from each of the 5 stacks

   // draw calibration stuff
   if(drawStr.Contains("Status") || drawStr.Contains("Gain") || drawStr.Contains("Noise") ||
      drawStr.Contains("Vdrift") || drawStr.Contains("T0") || 
      drawStr.Contains("gain") || drawStr.Contains("chiSquare")) {
     if(sectorNo==-1 && stackNo==-1) {    // plot the entire layer
       drawStr += Form(":PadSuperColumn%s:PadSuperRow%s>>prof", fAppendString.Data(), fAppendString.Data());
       drawStr += rndNumber;
       drawStr += "(76,-0.5,75.5,2592,-0.5,2591.5)";
       cutStr += Form("Layer==%d", layerNo);
     }
     else if(sectorNo!=-1 && stackNo==-1) {     // plot a sector from a layer
       drawStr += Form(":Column%s:PadSuperRow%s>>prof", fAppendString.Data(), fAppendString.Data());
       drawStr += rndNumber;
       drawStr += "(76,-0.5,75.5,144,-0.5,143.5)";
       cutStr += Form("Layer==%d && SuperModule==%d", layerNo, sectorNo);
     }
     else if(sectorNo==-1 && stackNo!=-1) {       // plot a stack from a layer
       drawStr += Form(":PadSuperColumn%s:Row%s>>prof", fAppendString.Data(), fAppendString.Data());
       drawStr += rndNumber;
       drawStr += Form("(%d,-0.5,%d-0.5,2592,-0.5,2591.5)", gkNRows[stackNo], gkNRows[stackNo]);
       cutStr += Form("Layer==%d && Stack==%d", layerNo, stackNo);
     }
     else {                // the layer, sector and stack are determined -> so plot a chamber
       drawStr += Form(":Column%s:Row%s>>prof", fAppendString.Data(), fAppendString.Data());
       drawStr += rndNumber;
       drawStr += Form("(%d,-0.5,%d-0.5,144,-0.5,143.5)", gkNRows[stackNo], gkNRows[stackNo]);
       cutStr += Form("Layer==%d && SuperModule==%d && Stack==%d", layerNo, sectorNo, stackNo);
     }
   }
   // draw FEE stuff
   else if(drawStr.Contains("SORandEOR") || 
	   drawStr.Contains("gsmSOR") || drawStr.Contains("gsmDelta") ||
	   drawStr.Contains("nimSOR") || drawStr.Contains("nimDelta") ||
	   drawStr.Contains("nevSOR") || drawStr.Contains("nevDelta") ||
	   drawStr.Contains("nptSOR") || drawStr.Contains("nptDelta")) {
     if(sectorNo==-1 && stackNo==-1) {    // plot the entire layer
       drawStr += Form(":MCMSuperColumn%s:MCMSuperRow%s>>prof", fAppendString.Data(), fAppendString.Data());
       drawStr += rndNumber;
       drawStr += "(76,-0.5,75.5,144,-0.5,143.5)";
       cutStr += Form("Layer==%d", layerNo);
     }
     else if(sectorNo!=-1 && stackNo==-1) {     // plot a sector from a layer
       drawStr += Form(":MCMSuperColumn%s:MCMSuperRow%s>>prof", fAppendString.Data(), fAppendString.Data());
       drawStr += rndNumber;
       drawStr += "(76,-0.5,75.5,144,-0.5,143.5)";
       cutStr += Form("Layer==%d && SuperModule==%d", layerNo, sectorNo);
     }
     else if(sectorNo==-1 && stackNo!=-1) {       // plot a stack from a layer
       drawStr += Form(":MCMSuperColumn%s:MCMSuperRow%s>>prof", fAppendString.Data(), fAppendString.Data());
       drawStr += rndNumber;
       //       drawStr += Form("(%d,-0.5,%d-0.5,2592,-0.5,2591.5)", gkNRows[stackNo], gkNRows[stackNo]);
       drawStr += "(76,-0.5,75.5,144,-0.5,143.5)";
       cutStr += Form("Layer==%d && Stack==%d", layerNo, stackNo);
     }
     else {                // the layer, sector and stack are determined -> so plot a chamber
       drawStr += Form(":ROB%s:MCM%s>>prof", fAppendString.Data(), fAppendString.Data());
       drawStr += rndNumber;
       drawStr += Form("(16,-0.5,15.5,%d,-0.5,%d-0.5)", gkNRows[stackNo]/2, gkNRows[stackNo]/2);
       cutStr += Form("Layer==%d && SuperModule==%d && Stack==%d", layerNo, sectorNo, stackNo);
     }
   }
   // draw alignment stuff
   else if(drawStr.Contains("Align")) {
     if(sectorNo==-1 && stackNo==-1) {  // plot the entire layer
       drawStr += ":SuperModule:Stack>>prof";
       drawStr += rndNumber;
       drawStr += "(5,-0.5,4.5,18,-0.5,17.5)";
       cutStr += Form("Layer==%d", layerNo);
     }
     else if(sectorNo!=-1 && stackNo==-1) {     // plot a sector from a layer
       drawStr += ":SuperModule:Stack>>prof";
       drawStr += rndNumber;
       drawStr += Form("(5,-0.5,4.5,1,%f,%f)", sectorNo-0.5, sectorNo+0.5);
       cutStr += Form("Layer==%d && SuperModule==%d", layerNo, sectorNo);
     }
     else if(sectorNo==-1 && stackNo!=-1) {       // plot a stack from a layer
       drawStr += ":SuperModule:Stack>>prof";
       drawStr += rndNumber;
       drawStr += Form("(1,%f,%f,18,-0.5,17.5)", stackNo-0.5, stackNo+0.5);
       cutStr += Form("Layer==%d && Stack==%d", layerNo, stackNo);
     }
     else {                // the layer, sector and stack are determined -> so plot a chamber
       drawStr += ":SuperModule:Stack>>prof";
       drawStr += rndNumber;
       drawStr += Form("(1,%f,%f,1,%f,%f)", stackNo-0.5, stackNo+0.5, sectorNo-0.5, sectorNo+0.5);
       cutStr += Form("Layer==%d && SuperModule==%d && Stack==%d", layerNo, sectorNo, stackNo);
     }
   }


   if (cuts && cuts[0] != 0) {
      if (cutStr.Length() != 0) cutStr += "&& ";
      cutStr += "(";
      cutStr += cuts;
      cutStr += ")";
   }
   drawStr.ReplaceAll(fAbbreviation, fAppendString);
   cutStr.ReplaceAll(fAbbreviation, fAppendString);
   if (writeDrawCommand) std::cout << "fTree->Draw(\"" << drawStr << "\", \"" <<  cutStr << "\", \"" << drawOptionsStr << "\");" << std::endl;
   Int_t returnValue = fTree->Draw(drawStr.Data(), cutStr.Data(), drawOptionsStr.Data());
   TString profName("prof");
   profName += rndNumber;
   TObject *obj = gDirectory->Get(profName.Data());
   // set the names of the axes 
   TH1 *histObj = (TH1*)obj;
   if(drawStr.Contains("Status") || drawStr.Contains("Gain") || drawStr.Contains("Noise") ||
      drawStr.Contains("Vdrift") || drawStr.Contains("T0") ||
      drawStr.Contains("gain") || drawStr.Contains("chiSquare")) {
     histObj->GetXaxis()->SetTitle("Row");
     histObj->GetYaxis()->SetTitle("Column");
   }
   else if(drawStr.Contains("SORandEOR") || 
	   drawStr.Contains("gsmSOR") || drawStr.Contains("gsmDelta") ||
	   drawStr.Contains("nimSOR") || drawStr.Contains("nimDelta") ||
	   drawStr.Contains("nevSOR") || drawStr.Contains("nevDelta") ||
	   drawStr.Contains("nptSOR") || drawStr.Contains("nptDelta")) {
     histObj->GetXaxis()->SetTitle("MCM Row");
     histObj->GetYaxis()->SetTitle("MCM Column");
   }
   else if(drawStr.Contains("Align")) {
     histObj->GetXaxis()->SetTitle("Stack");
     histObj->GetYaxis()->SetTitle("Sector");
   }

   if (obj && obj->InheritsFrom("TH1")) FormatHistoLabels((TH1*)obj);
   return returnValue;
}

//_____________________________________________________________________________
Int_t AliTRDCalibViewer::EasyDraw1D(const char* drawCommand, const char* sector, const char* cuts, const char* drawOptions, Bool_t writeDrawCommand) const {
  //
  // easy drawing of data, use '~' for abbreviation of '.fElements'
  // example: EasyDraw("CETmean~-CETmean_mean", "A", "(CETmean~-CETmean_mean)>0")
  // sector: sector-number - the specified sector will be drwawn
  //         'A'/'C' or 'a'/'c' - side A/C will be drawn
  //         'ALL' - whole TPC will be drawn, projected on one side
  // cuts: specifies cuts
  // drawOptions: draw options like 'same'
  // writeDrawCommand: write the command, that is passed to TTree::Draw
  //

   TString drawStr(drawCommand);

   TString sectorStr(sector);
   Int_t layerNo = -1; 
   Int_t sectorNo = -1; 
   Int_t stackNo = -1;
   GetLayerSectorStack(sectorStr, layerNo, sectorNo, stackNo);
   if(layerNo==-1) {
     Warning("EasyDraw", "The sector string must always contain the Layer number!");
     return -1;
   }
   if(layerNo<0 || layerNo>5) {
     Warning("EasyDraw", "The Layer number must be in the range [0,5] !");
     return -1;
   }
   if(sectorNo!=-1 && (sectorNo<0 || sectorNo>17)) {
     Warning("EasyDraw", "The Sector number must be in the range [0,17] !");
     return -1;
   }
   if(stackNo!=-1 && (stackNo<0 || stackNo>4)) {
     Warning("EasyDraw", "The Stack number must be in the range [0,4] !");
     return -1;
   }

   TString drawOptionsStr(drawOptions);
   TString cutStr("");

   if(sectorNo==-1 && stackNo==-1)     // plot the entire layer
     cutStr += Form("Layer==%d", layerNo);
   else if(sectorNo!=-1 && stackNo==-1)      // plot a sector from a layer
     cutStr += Form("Layer==%d && SuperModule==%d", layerNo, sectorNo);
   else if(sectorNo==-1 && stackNo!=-1)        // plot a stack from a layer
     cutStr += Form("Layer==%d && Stack==%d", layerNo, stackNo);
   else                 // the layer, sector and stack are determined -> so plot a chamber
     cutStr += Form("Layer==%d && SuperModule==%d && Stack==%d", layerNo, sectorNo, stackNo);
   
   if(cuts && cuts[0] != 0) {
      if (cutStr.Length() != 0) cutStr += "&& ";
      cutStr += "(";
      cutStr += cuts;
      cutStr += ")";
   }

   drawStr.ReplaceAll(fAbbreviation, fAppendString);
   cutStr.ReplaceAll(fAbbreviation, fAppendString);
   if (writeDrawCommand) std::cout << "fTree->Draw(\"" << drawStr << "\", \"" <<  cutStr << "\", \"" << drawOptionsStr << "\");" << std::endl;
   Int_t returnValue = fTree->Draw(drawStr.Data(), cutStr.Data(), drawOptionsStr.Data());
   if (returnValue == -1) return -1;
   
   TObject *obj = (gPad) ? gPad->GetPrimitive("htemp") : 0; 
   if (!obj) obj = (TH1F*)gDirectory->Get("htemp");
   if (!obj) obj = gPad->GetPrimitive("tempHist");
   if (!obj) obj = (TH1F*)gDirectory->Get("tempHist");
   if (!obj) obj = gPad->GetPrimitive("Graph");
   if (!obj) obj = (TH1F*)gDirectory->Get("Graph");
   if (obj && obj->InheritsFrom("TH1")) FormatHistoLabels((TH1*)obj);
   return returnValue;
}

//_____________________________________________________________________________
Int_t AliTRDCalibViewer::EasyDraw(const char* drawCommand, Int_t chamber, const char* cuts, const char* drawOptions, Bool_t writeDrawCommand) const {
  //
  // easy drawing of data, use '~' for abbreviation of '.fElements'
  // example: EasyDraw("CETmean~-CETmean_mean", 34, "(CETmean~-CETmean_mean)>0")
  // sector: sector-number - only the specified sector will be drwawn
  // cuts: specifies cuts
  // drawOptions: draw options like 'same'
  // writeDrawCommand: write the command, that is passed to TTree::Draw
  //
  if(chamber >= 0 && chamber < 540) {
    Int_t superModuleNo = chamber/30;
    Int_t stackNo = (chamber%30)/6;
    Int_t layerNo = (chamber%30)%6;
    char sectorChr[22];
    sprintf(sectorChr, "Layer%iSector%iStack%i", layerNo, superModuleNo, stackNo);
    return EasyDraw(drawCommand, sectorChr, cuts, drawOptions, writeDrawCommand);
  }
  Error("EasyDraw","The TRD contains only chamber from 0 to 539");
  return -1;
}

//_____________________________________________________________________________
Int_t AliTRDCalibViewer::EasyDraw1D(const char* drawCommand, Int_t chamber, const char* cuts, const char* drawOptions, Bool_t writeDrawCommand) const {
  //
  // easy drawing of data, use '~' for abbreviation of '.fElements'
  // example: EasyDraw("CETmean~-CETmean_mean", 34, "(CETmean~-CETmean_mean)>0")
  // sector: sector-number - the specified sector will be drwawn
  // cuts: specifies cuts
  // drawOptions: draw options like 'same'
  // writeDrawCommand: write the command, that is passed to TTree::Draw
  //

  if (chamber >= 0 && chamber < 539) {
    Int_t superModuleNo = chamber/30;
    Int_t stackNo = (chamber%30)/6;
    Int_t layerNo = (chamber%30)%6;
    char sectorChr[22];
    sprintf(sectorChr, "Layer%iSector%iStack%i", layerNo, superModuleNo, stackNo);
    return EasyDraw1D(drawCommand, sectorChr, cuts, drawOptions, writeDrawCommand);
  }
  Error("EasyDraw1D","The TRD contains only chambers from 0 to 539");
  return -1;
}

//_____________________________________________________________________________
Bool_t AliTRDCalibViewer::DumpOCDBtoTreeDetails(const Char_t* runListFilename,
						const Char_t* outFilename,
						Int_t firstRun, Int_t lastRun,
						const Char_t* storage,
                                                Int_t version,
                                                Int_t subVersion,
                                                Bool_t getCalibs,
                                                Bool_t getDCS,
                                                Bool_t getAlign) {
  //
  // Retrieve TRD OCDB information for a given run list/range
  //

  if(runListFilename[0]!='\0' && firstRun==-1 && lastRun==-1) {
    cout << "AliTRDCalibViewer::DumpOCDBtoTreeDetails(): You must provide at least a run range or an ascii filename with run numbers" 
	 << endl;
    return kFALSE;
  }
  // initialize the OCDB manager
  TString storageString = storage;
  if(storageString.Contains("alien://")) {
    TGrid::Connect("alien://");
  }
  AliCDBManager *manager = AliCDBManager::Instance();
  if(storage[0]!='\0') {
    manager->SetDefaultStorage(storage);
  }
  else {
    if(!manager->IsDefaultStorageSet()) {
      cout << "AliTRDCalibViewer::DumpOCDBtoTreeDetails(): Default OCDB storage not set!!" << endl;
      return kFALSE;
    }
  }
  manager->SetRun(1);

  // open the ascii file
  ifstream in;
  if(runListFilename[0]!='\0')
    in.open(runListFilename);

  // initialize the tree streamer
  if(outFilename[0]=='\0') outFilename = "trdDetails.root";
  TString calibFilename = outFilename;
  
  TTreeSRedirector *treeStreamer = new TTreeSRedirector(calibFilename.Data());

  Int_t currRun;
  if(runListFilename[0]=='\0' && firstRun!=-1 && lastRun!=-1)
    currRun = firstRun;

  TVectorD runs;

  // loop over runs
  while(1) {
    if(runListFilename[0]!='\0') {
      if(!(in>>currRun)) continue;
      if(currRun < (firstRun==-1 ? 0 : firstRun) ||
	 currRun > (lastRun==-1 ? 999999999 : lastRun))
	continue;
    }
    else {
      if(currRun>lastRun) break;
    }
    cout << "run = " << currRun << endl;
    manager->SetRun(currRun);

    // Get GRP data. If there is no proper GRP object for this run than
    // this run is aborted
    AliCDBEntry *entry = manager->Get("GRP/GRP/Data");
    AliGRPObject* grpObject = 0;
    if(entry) {
      entry->SetOwner(kFALSE);
      grpObject = dynamic_cast<AliGRPObject*>(entry->GetObject());
    }
    else {
      currRun++;
      //      continue;
      //      return kFALSE;
    }
    if(!grpObject)
      cout << "No GRP info available for this run " << endl;

    time_t startTimeGRP = 0;
    TObjString runType("");
    if(grpObject) {
      startTimeGRP = grpObject->GetTimeStart();
      TTimeStamp start(grpObject->GetTimeStart());
      TTimeStamp end(grpObject->GetTimeEnd());
      cout << "Start time: " << start.GetDate()/10000 << "/" 
	   << (start.GetDate()/100)-(start.GetDate()/10000)*100 << "/" 
	   << start.GetDate()%100 << "   "
	   << start.GetTime()/10000 << ":"
	   << (start.GetTime()/100)-(start.GetTime()/10000)*100 << ":" 
	   << start.GetTime()%100 << endl;
      cout << "End time: " << end.GetDate()/10000 << "/" 
	   << (end.GetDate()/100)-(end.GetDate()/10000)*100 << "/" 
	   << end.GetDate()%100 << "   "
	   << end.GetTime()/10000 << ":"
	   << (end.GetTime()/100)-(end.GetTime()/10000)*100 << ":"
	   << end.GetTime()%100 << endl;
      cout << "Run type = " << grpObject->GetRunType().Data() << endl;
      runType = grpObject->GetRunType().Data();
    }

    // gain
    AliTRDCalDet *chamberGainFactor = 0;
    if(getCalibs) {
      entry = manager->Get("TRD/Calib/ChamberGainFactor", currRun, version, subVersion);
      if(entry) {
        entry->SetOwner(kFALSE);
        chamberGainFactor = (AliTRDCalDet*)entry->GetObject();
      }
    }
    AliTRDCalPad *padGainFactor = 0;
    if(getCalibs) {
      entry = manager->Get("TRD/Calib/LocalGainFactor", currRun, version, subVersion);
      if(entry) {
        entry->SetOwner(kFALSE);
        padGainFactor = (AliTRDCalPad*)entry->GetObject();
      }
    }
    Double_t runMeanGain, runRMSGain;
    TVectorD chamberMeanGain(AliTRDcalibDB::kNdet);
    TVectorD chamberRMSGain(AliTRDcalibDB::kNdet);
    TVectorD smMeanGain(AliTRDcalibDB::kNsector);
    TVectorD smRMSGain(AliTRDcalibDB::kNsector);
    for(Int_t iNdet=0; iNdet<AliTRDcalibDB::kNdet; iNdet++) {chamberMeanGain[iNdet] = 0.0; chamberRMSGain[iNdet] = 0.0;}
    for(Int_t iSm=0; iSm<AliTRDcalibDB::kNsector; iSm++) {smMeanGain[iSm] = 0.0; smRMSGain[iSm] = 0.0;}
    TString parName("Gain");
    if(getCalibs)
      ProcessTRDCalibArray(chamberGainFactor, padGainFactor, 
			   parName,
			   runMeanGain, runRMSGain,
			   chamberMeanGain, chamberRMSGain,
			   smMeanGain, smRMSGain);

    // noise/pedestals
    AliTRDCalDet *chamberNoise = 0;
    if(getCalibs) {
      entry = manager->Get("TRD/Calib/DetNoise", currRun, version, subVersion);
      if(entry) {
        entry->SetOwner(kFALSE);
        chamberNoise = (AliTRDCalDet*)entry->GetObject();
      }
    }
    AliTRDCalPad *padNoise = 0;
    if(getCalibs) {
      entry = manager->Get("TRD/Calib/PadNoise", currRun, version, subVersion);
      if(entry) {
        entry->SetOwner(kFALSE);
        padNoise = (AliTRDCalPad*)entry->GetObject();
      }
    }
    Double_t runMeanNoise, runRMSNoise;
    TVectorD chamberMeanNoise(AliTRDcalibDB::kNdet);
    TVectorD chamberRMSNoise(AliTRDcalibDB::kNdet);
    TVectorD smMeanNoise(AliTRDcalibDB::kNsector);
    TVectorD smRMSNoise(AliTRDcalibDB::kNsector);
    for(Int_t iNdet=0; iNdet<AliTRDcalibDB::kNdet; iNdet++) {chamberMeanNoise[iNdet] = 0.0; chamberRMSNoise[iNdet] = 0.0;}
    for(Int_t iSm=0; iSm<AliTRDcalibDB::kNsector; iSm++) {smMeanNoise[iSm] = 0.0; smRMSNoise[iSm] = 0.0;}
    parName = "Noise";
    if(getCalibs)
      ProcessTRDCalibArray(chamberNoise, padNoise, 
			   parName,
			   runMeanNoise, runRMSNoise,
			   chamberMeanNoise, chamberRMSNoise,
			   smMeanNoise, smRMSNoise);

    // vdrift
    AliTRDCalDet *chamberVdrift = 0;
    if(getCalibs) {
      entry = manager->Get("TRD/Calib/ChamberVdrift", currRun, version, subVersion);
      if(entry) {
        entry->SetOwner(kFALSE);
        chamberVdrift = (AliTRDCalDet*)entry->GetObject();
      }
    }
    AliTRDCalPad *padVdrift = 0;
    if(getCalibs) {
      entry = manager->Get("TRD/Calib/LocalVdrift", currRun, version, subVersion);
      if(entry) {
        entry->SetOwner(kFALSE);
        padVdrift = (AliTRDCalPad*)entry->GetObject();
      }
    }
    Double_t runMeanVdrift, runRMSVdrift;
    TVectorD chamberMeanVdrift(AliTRDcalibDB::kNdet);
    TVectorD chamberRMSVdrift(AliTRDcalibDB::kNdet);
    TVectorD smMeanVdrift(AliTRDcalibDB::kNsector);
    TVectorD smRMSVdrift(AliTRDcalibDB::kNsector);
    for(Int_t iNdet=0; iNdet<AliTRDcalibDB::kNdet; iNdet++) {chamberMeanVdrift[iNdet] = 0.0; chamberRMSVdrift[iNdet] = 0.0;}
    for(Int_t iSm=0; iSm<AliTRDcalibDB::kNsector; iSm++) {smMeanVdrift[iSm] = 0.0; smRMSVdrift[iSm] = 0.0;}
    parName = "Vdrift";
    if(getCalibs)
      ProcessTRDCalibArray(chamberVdrift, padVdrift, 
			   parName,
			   runMeanVdrift, runRMSVdrift,
			   chamberMeanVdrift, chamberRMSVdrift,
			   smMeanVdrift, smRMSVdrift);

    // T0
    AliTRDCalDet *chamberT0 = 0;
    if(getCalibs) {
      entry = manager->Get("TRD/Calib/ChamberT0", currRun, version, subVersion);
      if(entry) {
        entry->SetOwner(kFALSE);
        chamberT0 = (AliTRDCalDet*)entry->GetObject();
      }
    }
    AliTRDCalPad *padT0 = 0;
    if(getCalibs) {
      entry = manager->Get("TRD/Calib/LocalT0", currRun, version, subVersion);
      if(entry) {
        entry->SetOwner(kFALSE);
        padT0 = (AliTRDCalPad*)entry->GetObject();
      }
    }
    Double_t runMeanT0, runRMST0;
    TVectorD chamberMeanT0(AliTRDcalibDB::kNdet);
    TVectorD chamberRMST0(AliTRDcalibDB::kNdet);
    TVectorD smMeanT0(AliTRDcalibDB::kNsector);
    TVectorD smRMST0(AliTRDcalibDB::kNsector);
    for(Int_t iNdet=0; iNdet<AliTRDcalibDB::kNdet; iNdet++) {chamberMeanT0[iNdet] = 0.0; chamberRMST0[iNdet] = 0.0;}
    for(Int_t iSm=0; iSm<AliTRDcalibDB::kNsector; iSm++) {smMeanT0[iSm] = 0.0; smRMST0[iSm] = 0.0;}
    parName = "T0";
    if(getCalibs)
      ProcessTRDCalibArray(chamberT0, padT0, 
			   parName,
			   runMeanT0, runRMST0,
			   chamberMeanT0, chamberRMST0,
			   smMeanT0, smRMST0);

    // status
    AliTRDCalChamberStatus* chamberStatus = 0;
    if(getCalibs) {
      entry = manager->Get("TRD/Calib/ChamberStatus", currRun, version, subVersion);
      if(entry) {
        entry->SetOwner(kFALSE);
        chamberStatus = (AliTRDCalChamberStatus*)entry->GetObject();
      }
    }
    AliTRDCalPadStatus *padStatus = 0;
    if(getCalibs) {
      entry = manager->Get("TRD/Calib/PadStatus", currRun, version, subVersion);
      if(entry) {
        entry->SetOwner(kFALSE);
        padStatus = (AliTRDCalPadStatus*)entry->GetObject();
      }
    }

    // DCS FEE information
    TObjArray *dcsArray = 0;
    if(getDCS) {
      entry = manager->Get("TRD/Calib/DCS");
      if(entry) {
        entry->SetOwner(kTRUE);
        dcsArray = (TObjArray*)entry->GetObject();
      }
    }
    AliTRDCalDCS *dcsSOR = 0;
    AliTRDCalDCS *dcsEOR = 0;
    if(getDCS && dcsArray) {
      dcsSOR = (AliTRDCalDCS*)dcsArray->At(0);
      dcsEOR = (AliTRDCalDCS*)dcsArray->At(1);
    }

    // Alignment information
    // get the geometry from OCDB
    TGeoManager *geoMan = 0x0;
    if(getAlign) {
      entry=manager->Get("GRP/Geometry/Data");
      if(entry)
        geoMan=(TGeoManager*)entry->GetObject();
      else
        cout << "Cannot get an entry for the geometry storage" << endl;
    }
    // get the alignment from OCDB
    AliTRDalignment *alignMan=0;
    if(getAlign && geoMan) {
      entry=manager->Get("TRD/Align/Data", currRun, version, subVersion);
      if(entry) {
        alignMan = new AliTRDalignment();
        cout << "storage for alignment = " << manager->GetDefaultStorage()->GetURI().Data() << endl;
        alignMan->ReadDB(manager->GetDefaultStorage()->GetURI().Data(), "TRD/Align/Data", currRun, version, subVersion);
      }
      else {
        cout << "Cannot get an entry for the alignment info" << endl;
      }
    }

    Int_t kSuperModuleStatus[18] = {1, 1, 0, 0, 0, 0,     // super module status (1- installed, 0- not installed)
				    0, 1, 1, 1, 1, 0, 
				    0, 0, 0, 0, 0, 1};
    Int_t kNRows[ 5] = {16, 16, 12, 16, 16};  // number of pad rows in the chambers from each of the 5 stacks
    Int_t kNCols = 144;          // number of pad columns in the chambers from each of the 18 supermodules
    Int_t kROB[5] = {8, 8, 6, 8, 8};   // number of read out boards(ROB) per chamber (6 in stack 2 and 8 in the rest)
    Int_t kMCM = 16;                   // number of MCMs per ROB
    for(Short_t iLayer=0; iLayer<AliTRDgeometry::kNlayer; iLayer++) {   // loop over layers
      for(Short_t iSector=0; iSector<AliTRDgeometry::kNsector; iSector++) {  // loop over supermodules
	if(kSuperModuleStatus[iSector]==0) 
	  continue;
	Double_t alignSMPars[6];
        for(Int_t ipar=0; ipar<6; ipar++) alignSMPars[ipar]=0.0;
        if(getAlign && alignMan)
	  alignMan->GetSm(iSector, alignSMPars);
	for(Short_t iStack=0; iStack<AliTRDgeometry::kNstack; iStack++) {    // loop over stacks
	  Short_t chamberNo = AliTRDgeometry::GetDetector(iLayer, iStack, iSector);
	  AliTRDCalROC *gainROC = 0;
	  if(padGainFactor) gainROC = padGainFactor->GetCalROC(chamberNo);
	  AliTRDCalROC *noiseROC = 0;
	  if(padNoise) noiseROC = padNoise->GetCalROC(chamberNo);
	  AliTRDCalROC *vdriftROC = 0;
	  if(padVdrift) vdriftROC = padVdrift->GetCalROC(chamberNo);
	  AliTRDCalROC *t0ROC = 0;
	  if(t0ROC) t0ROC = padT0->GetCalROC(chamberNo);
	  AliTRDCalSingleChamberStatus *statusROC = 0;
	  if(padStatus) statusROC = padStatus->GetCalROC(chamberNo);
	  TVectorD channelVector(kNRows[iStack]*kNCols);
	  TVectorD rowVector(kNRows[iStack]*kNCols);
	  TVectorD colVector(kNRows[iStack]*kNCols);
	  TVectorD statusVector(kNRows[iStack]*kNCols);
	  TVectorD gainVector(kNRows[iStack]*kNCols);
	  TVectorD noiseVector(kNRows[iStack]*kNCols);
	  TVectorD vdriftVector(kNRows[iStack]*kNCols);
	  TVectorD t0Vector(kNRows[iStack]*kNCols);
	  TVectorD padSuperRowVector(kNRows[iStack]*kNCols);
	  TVectorD padSuperColumnVector(kNRows[iStack]*kNCols);
          for(Int_t ipar=0; ipar<kNRows[iStack]*kNCols; ipar++) {
            channelVector[ipar] = 0; rowVector[ipar] = 0; colVector[ipar] = 0; 
            statusVector[ipar] = 0; gainVector[ipar] = 0; noiseVector[ipar] = 0; 
            vdriftVector[ipar] = 0; t0Vector[ipar] = 0; padSuperRowVector[ipar] = 0; 
            padSuperColumnVector[ipar] = 0;
          }
	  Int_t index = 0;
	  if(getCalibs) {
            for(Short_t iRow=0; iRow<kNRows[iStack]; iRow++) {   // loop over pad rows
	      for(Short_t iCol=0; iCol<kNCols; iCol++) {    // loop over pad columns
	        Short_t padSuperRow = iRow;
	        for(Int_t i=0; i<iStack; i++) padSuperRow = padSuperRow + kNRows[i]; 
	        padSuperRowVector[index] = padSuperRow;
	        Short_t padSuperColumn = iCol;
	        for(Int_t i=0; i<iSector; i++) padSuperColumn = padSuperColumn + kNCols;
	        padSuperColumnVector[index] = padSuperColumn;
	        Short_t channelNo = -1;
	        Float_t gain = -99.;
	        if(gainROC && chamberGainFactor) {
		  channelNo = gainROC->GetChannel(iCol, iRow);
		  gain = chamberGainFactor->GetValue(chamberNo) * gainROC->GetValue(iCol, iRow);
	        }
	        Float_t noise = -99.;
	        if(noiseROC && chamberNoise)
		  noise = chamberNoise->GetValue(chamberNo) * noiseROC->GetValue(iCol, iRow);
	        Float_t vdrift = -99.;
	        if(vdriftROC && chamberVdrift)
		  vdrift = chamberVdrift->GetValue(chamberNo) * vdriftROC->GetValue(iCol, iRow);
	        Float_t t0 = -99.;
	        if(t0ROC && chamberT0)
		  t0 = chamberT0->GetValue(chamberNo) + t0ROC->GetValue(iCol, iRow);
	        Int_t status = -99;
	        if(statusROC)
		  status = statusROC->GetStatus(iCol, iRow);
	        channelVector[index] = channelNo;
	        rowVector[index] = iRow;
	        colVector[index] = iCol;
	        statusVector[index] = status;
	        gainVector[index] = gain;
	        noiseVector[index] = noise;
	        vdriftVector[index] = vdrift;
	        t0Vector[index] = t0;
	        index++;
	      }  // end loop over pad columns
	    }  // end loop over pad rows
          }   // end if(getCalibs)

	  // get the dcs information
	  AliTRDCalDCSFEE *dcsfeeSOR = 0;
	  AliTRDCalDCSFEE *dcsfeeEOR = 0;
          if(getDCS) {
	    if(dcsSOR) dcsfeeSOR = dcsSOR->GetCalDCSFEEObj(chamberNo);
	    if(dcsEOR) dcsfeeEOR = dcsEOR->GetCalDCSFEEObj(chamberNo);
          }
          
	  Bool_t sorAndEor = kFALSE;
	  if(getDCS && dcsfeeSOR && dcsfeeEOR) sorAndEor = kTRUE;
	  if(getDCS && !dcsfeeSOR && dcsfeeEOR) dcsfeeSOR = dcsfeeEOR;
	  TVectorD robVector(kROB[iStack]*kMCM);
	  TVectorD mcmVector(kROB[iStack]*kMCM);
	  TVectorD sorandeorVector(kROB[iStack]*kMCM);
	  TVectorD gsmSorVector(kROB[iStack]*kMCM);
	  TVectorD gsmDeltaVector(kROB[iStack]*kMCM);
	  TVectorD nimSorVector(kROB[iStack]*kMCM);
	  TVectorD nimDeltaVector(kROB[iStack]*kMCM);
	  TVectorD nevSorVector(kROB[iStack]*kMCM);
	  TVectorD nevDeltaVector(kROB[iStack]*kMCM);
	  TVectorD nptSorVector(kROB[iStack]*kMCM);
	  TVectorD nptDeltaVector(kROB[iStack]*kMCM);
	  TVectorD mcmSuperRowVector(kROB[iStack]*kMCM);
	  TVectorD mcmSuperColumnVector(kROB[iStack]*kMCM);
          for(Int_t ipar=0; ipar<kROB[iStack]*kMCM; ipar++) {
            robVector[ipar] = 0; mcmVector[ipar] = 0; sorandeorVector[ipar] = 0; 
            gsmSorVector[ipar] = 0; gsmDeltaVector[ipar] = 0; nimSorVector[ipar] = 0; 
            nimDeltaVector[ipar] = 0; nevSorVector[ipar] = 0; nevDeltaVector[ipar] = 0; 
            nptSorVector[ipar] = 0; nptDeltaVector[ipar] = 0; mcmSuperRowVector[ipar] = 0; 
            mcmSuperColumnVector[ipar] = 0;
          }

	  Int_t robsRowDirection = kNRows[iStack]/4;    // 4 or 3 ROBs per chamber in row direction
	  Int_t index1 = 0; 
          if(getDCS && (dcsfeeSOR || dcsfeeEOR) && dcsfeeSOR->GetStatusBit()==0) {
	    for(Int_t iROB=0; iROB<kROB[iStack]; iROB++) { // loop over ROBs
	      for(Int_t iMCM=0; iMCM<kMCM; iMCM++) {  // loop over MCMs
	        Short_t superRowMCM = iMCM%4;            // 4 MCMs per ROB in row direction
	        superRowMCM += 4*(iROB%robsRowDirection);   // now we have the row of this MCM inside one chamber
	        for(Int_t kk=0; kk<iStack; kk++) superRowMCM += kNRows[kk];   // add number of rows in previous stacks
	        Short_t superColumnMCM = iMCM/4;        // 4 MCMs per ROB in column direction
	        superColumnMCM += 4*(iROB/robsRowDirection);    // should yield 0 or 1 (2 ROBs per chamber in col direction)
	        superColumnMCM += iSector*8;
	        mcmSuperRowVector[index1] = superRowMCM;
	        mcmSuperColumnVector[index1] = superColumnMCM;
	        Int_t gsm = dcsfeeSOR->GetMCMGlobalState(iROB, iMCM);
	        Int_t nim = dcsfeeSOR->GetMCMStateNI(iROB, iMCM);
	        Int_t nev = dcsfeeSOR->GetMCMEventCnt(iROB, iMCM);
	        Int_t npt = dcsfeeSOR->GetMCMPtCnt(iROB, iMCM);
	        Int_t dgsm = -100000;
	        Int_t dnim = -100000;
	        Int_t dnev = -100000;
	        Int_t dnpt = -100000;
	        if(sorAndEor) {
		  dgsm = gsm - dcsfeeEOR->GetMCMGlobalState(iROB, iMCM);
		  dnim = nim - dcsfeeEOR->GetMCMStateNI(iROB, iMCM);
		  dnev = nev - dcsfeeEOR->GetMCMEventCnt(iROB, iMCM);
		  dnpt = npt - dcsfeeEOR->GetMCMPtCnt(iROB, iMCM);
		  if(gsm==-1 && dgsm==0) dgsm = -100000;
		  if(nim==-1 && dnim==0) dnim = -100000;
		  if(nev==-1 && dnev==0) dnev = -100000;
		  if(npt==-1 && dnpt==0) dnpt = -100000;
	        }
	        robVector[index1] = iROB;
	        mcmVector[index1] = iMCM;
	        sorandeorVector[index1] = sorAndEor;
	        gsmSorVector[index1] = gsm;
	        gsmDeltaVector[index1] = dgsm;
	        nimSorVector[index1] = nim;
	        nimDeltaVector[index1] = dnim;
	        nevSorVector[index1] = nev;
	        nevDeltaVector[index1] = dnev;
	        nptSorVector[index1] = npt;
	        nptDeltaVector[index1] = dnpt;
	        index1++;
	      }  // end loop over MCMs
	    }  // end loop over ROBs
          }  // end if(getDCS ...)

	  Double_t alignChamberPars[6];
          for(Int_t ipar=0; ipar<6; ipar++) alignChamberPars[ipar]=0;
          if(getAlign && alignMan)
	    alignMan->GetCh(chamberNo, alignChamberPars);

	  (*treeStreamer)<< "TRDcalibDetails"
			 << "run=" << currRun
			 << "SuperModule=" << iSector
			 << "Stack=" << iStack
                         << "Layer=" << iLayer
			 << "Chamber=" << chamberNo;
          if(getAlign)
            (*treeStreamer)<< "TRDcalibDetails"
			   << "Align_SM_ShiftRphi=" << alignSMPars[0]
			   << "Align_SM_ShiftZ=" << alignSMPars[1]
			   << "Align_SM_ShiftR=" << alignSMPars[2]
			   << "Align_SM_RotRphi=" << alignSMPars[3]
			   << "Align_SM_RotZ=" << alignSMPars[4]
			   << "Align_SM_RotR=" << alignSMPars[5]
			   << "Align_Ch_ShiftRphi=" << alignChamberPars[0]
			   << "Align_Ch_ShiftZ=" << alignChamberPars[1]
			   << "Align_Ch_ShiftR=" << alignChamberPars[2]
			   << "Align_Ch_RotRphi=" << alignChamberPars[3]
			   << "Align_Ch_RotZ=" << alignChamberPars[4]
			   << "Align_Ch_RotR=" << alignChamberPars[5];
          if(getCalibs)
            (*treeStreamer)<< "TRDcalibDetails"
			   << "Gain_Mean_Run=" << runMeanGain
			   << "Gain_RMS_Run=" << runRMSGain
			   << "Gain_Mean_SuperModule=" << smMeanGain[iSector]
			   << "Gain_RMS_SuperModule=" << smRMSGain[iSector]
			   << "Gain_Mean_Chamber=" << chamberMeanGain[chamberNo]
			   << "Gain_RMS_Chamber=" << chamberRMSGain[chamberNo]
			   << "Noise_Mean_Run=" << runMeanNoise
			   << "Noise_RMS_Run=" << runRMSNoise
			   << "Noise_Mean_SuperModule=" << smMeanNoise[iSector]
			   << "Noise_RMS_SuperModule=" << smRMSNoise[iSector]
			   << "Noise_Mean_Chamber=" << chamberMeanNoise[chamberNo]
			   << "Noise_RMS_Chamber=" << chamberRMSNoise[chamberNo]
			   << "Vdrift_Mean_Run=" << runMeanVdrift
			   << "Vdrift_RMS_Run=" << runRMSVdrift
			   << "Vdrift_Mean_SuperModule=" << smMeanVdrift[iSector]
			   << "Vdrift_RMS_SuperModule=" << smRMSVdrift[iSector]
			   << "Vdrift_Mean_Chamber=" << chamberMeanVdrift[chamberNo]
			   << "Vdrift_RMS_Chamber=" << chamberRMSVdrift[chamberNo]
			   << "T0_Mean_Run=" << runMeanT0
			   << "T0_RMS_Run=" << runRMST0
			   << "T0_Mean_SuperModule=" << smMeanT0[iSector]
			   << "T0_RMS_SuperModule=" << smRMST0[iSector]
			   << "T0_Mean_Chamber=" << chamberMeanT0[chamberNo]
			   << "T0_RMS_Chamber=" << chamberRMST0[chamberNo]
                           << "Channel.=" << &channelVector
			   << "Row.=" << &rowVector
			   << "Column.=" << &colVector
			   << "PadSuperRow.=" << &padSuperRowVector
			   << "PadSuperColumn.=" << &padSuperColumnVector
			   << "Status.=" << &statusVector
			   << "Gain.=" << &gainVector
			   << "Noise.=" << &noiseVector
			   << "Vdrift.=" << &vdriftVector
			   << "T0.=" << &t0Vector;
          if(getDCS)
	    (*treeStreamer)<< "TRDcalibDetails"
			   << "ROB.=" << &robVector
			   << "MCM.=" << &mcmVector
			   << "MCMSuperRow.=" << &mcmSuperRowVector
			   << "MCMSuperColumn.=" << &mcmSuperColumnVector
			   << "SORandEOR.=" << &sorandeorVector
			   << "gsmSOR.=" << &gsmSorVector
			   << "gsmDelta.=" << &gsmDeltaVector
			   << "nimSOR.=" << &nimSorVector
			   << "nimDelta.=" << &nimDeltaVector
			   << "nevSOR.=" << &nevSorVector
			   << "nevDelta.=" << &nevDeltaVector
			   << "nptSOR.=" << &nptSorVector
			   << "nptDelta.=" << &nptDeltaVector;
          (*treeStreamer)<< "TRDcalibDetails"
			 << "\n";
	}  // end loop over stacks
      }  // end loop over supermodules
    }  // end loop over layers
    
    // add the run number to the list of runs
    runs.ResizeTo(runs.GetNoElements()+1);
    runs[runs.GetNoElements()-1] = currRun;
    
    // do cleaning
    if(chamberGainFactor) delete chamberGainFactor;
    if(padGainFactor) delete padGainFactor;
    if(chamberNoise) delete chamberNoise;
    if(padNoise) delete padNoise;
    if(chamberVdrift) delete chamberVdrift;
    if(padVdrift) delete padVdrift;
    if(chamberT0) delete chamberT0;
    if(padT0) delete padT0;
    if(chamberStatus) delete chamberStatus;
    if(padStatus) delete padStatus;

    // check if we still have run numbers in the file or provided range
    if(runListFilename[0]=='\0' && firstRun!=-1 && lastRun!=-1) {
      currRun++;
      if(currRun>lastRun) break;
    }
    if(runListFilename[0]!='\0' && in.eof())
      break;
  }   // end loop over runs

  treeStreamer->GetFile()->cd();
  runs.Write("runs");
  delete treeStreamer;
  return kTRUE;
  //  delete treeStreamerDCS;
}

//_________________________________________________________________________
void AliTRDCalibViewer::DumpCalibToTree(const Char_t* inFilename, const Char_t* outFilename)
{
  //
  //  extract info from CalPad objects and dump them into a tree to be viewed
  //
  TTreeSRedirector *treeStreamer = new TTreeSRedirector(outFilename);
  //open file and retrieve list of calPad objects
  TFile f(inFilename);
  TList *l=(TList*)f.GetListOfKeys();

  TObjArray arrCalPads;
  TObjArray arrSMmean;
  TObjArray arrSMrms;
  arrCalPads.SetOwner();
  arrSMmean.SetOwner();
  arrSMrms.SetOwner();

  TIter next(l);
  TKey *k=0x0;
  while ( (k=(TKey*)next()) ){
    AliTRDCalPad *pad=dynamic_cast<AliTRDCalPad*>(k->ReadObj());
    if (!pad) continue;
    arrCalPads.Add(pad);

    TVectorD *smMean=new TVectorD(AliTRDcalibDB::kNsector);
    TVectorD *smRMS=new TVectorD(AliTRDcalibDB::kNsector);

    arrSMmean.Add(smMean);
    arrSMrms.Add(smRMS);

    ProcessTRDCalibArray(pad, *smMean, *smRMS);
  }

  Int_t kSuperModuleStatus[18] = {1, 1, 0, 0, 0, 0,     // super module status (1- installed, 0- not installed)
      0, 1, 1, 1, 1, 0,
      0, 0, 0, 0, 0, 1};

  AliTRDgeometry trdGeom;
  Int_t kNRows[5] = {16, 16, 12, 16, 16};  // number of pad rows in the chambers from each of the 5 stacks
  
  for(Short_t iLayer=0; iLayer<AliTRDgeometry::kNlayer; iLayer++) {   // loop over layers
    for(Short_t iSM=0; iSM<AliTRDgeometry::kNsector; iSM++) {  // loop over supermodules
      if(kSuperModuleStatus[iSM]==0)
        continue;
      
      for(Short_t iStack=0; iStack<AliTRDgeometry::kNstack; iStack++) {    // loop over stacks
        AliTRDpadPlane &plane=*trdGeom.GetPadPlane(iLayer, iStack);

        Short_t chamberNo = AliTRDgeometry::GetDetector(iLayer, iStack, iSM);
        const Int_t nrows=plane.GetNrows();
        const Int_t ncols=plane.GetNcols();
        const Int_t nchannels=nrows*ncols;
//         printf("chamberNo: %d (%03d-%03d-%03d)\n", chamberNo,nrows,ncols,nchannels);
        
        TVectorD channelVector(nchannels);
        TVectorD rowVector(nchannels);
        TVectorD colVector(nchannels);
        
        TVectorD gxVector(nchannels);
        TVectorD gyVector(nchannels);
        TVectorD gzVector(nchannels);
        
        TVectorD padSuperRowVector(nchannels);
        TVectorD padSuperColumnVector(nchannels);
        
        Int_t index = 0;
        for(Short_t iRow=0; iRow<nrows; iRow++) {   // loop over pad rows
          for(Short_t iCol=0; iCol<ncols; iCol++) {    // loop over pad columns
            Short_t padSuperRow = iRow;
            for(Int_t i=0; i<iStack; i++) padSuperRow = padSuperRow + kNRows[i];
            padSuperRowVector.GetMatrixArray()[index] = padSuperRow;
            
            Short_t padSuperColumn = iCol;
            for(Int_t i=0; i<iSM; i++) padSuperColumn = padSuperColumn + ncols;
            padSuperColumnVector.GetMatrixArray()[index] = padSuperColumn;
            
            rowVector.GetMatrixArray()[index] = iRow;
            colVector.GetMatrixArray()[index] = iCol;
            
            index++;
          }  // end loop over pad columns
        }  // end loop over pad rows
        
        (*treeStreamer)<< "TRDcalibDetails"
          << "SuperModule=" << iSM
          << "Stack=" << iStack
          << "Layer=" << iLayer
          << "Chamber=" << chamberNo
            //geographical information
          << "Channel.=" << &channelVector
          << "Row.=" << &rowVector
          << "Column.=" << &colVector
          << "PadSuperRow.=" << &padSuperRowVector
          << "PadSuperColumn.=" << &padSuperColumnVector;
//          << "gx.=" << &gxVector
//          << "gy.=" << &gyVector
//          << "gz.=" << &gzVector;
          
        //
        // pad calibrations
        //
        TObjArray arrTrash;
        arrTrash.SetOwner();
        Int_t ncalib=arrCalPads.GetEntriesFast();
        for (Int_t iCalib=0; iCalib<ncalib; ++iCalib){
          AliTRDCalPad *pad=(AliTRDCalPad*)arrCalPads.UncheckedAt(iCalib);
          AliTRDCalROC *calROC=pad->GetCalROC(chamberNo);
          
          TVectorD &smMean=*((TVectorD*)arrSMmean.UncheckedAt(iCalib));
          TVectorD &smRMS=*((TVectorD*)arrSMrms.UncheckedAt(iCalib));
          
          TString calibName=pad->GetName();
          
          TVectorD *valueVector=new TVectorD(nchannels);
          arrTrash.Add(valueVector);
          
          Double_t rocMean=0;
          Double_t rocRMS=0;
          Double_t rocMedian=0;
          
          if (calROC){
            Int_t index = 0;
            for(Short_t iRow=0; iRow<nrows; iRow++) {
              for(Short_t iCol=0; iCol<ncols; iCol++) {
                valueVector->GetMatrixArray()[index] = calROC->GetValue(iCol,iRow);
                index++;
              }
            }
            rocMean   = calROC->GetMean();
            rocRMS    = calROC->GetRMS();
            rocMedian = calROC->GetMedian();
            //check for NaN
            if ( !(rocMean<1e30) ) rocMean=0;
            if ( !(rocRMS<1e30) ) rocRMS=0;
            
//             printf("mean:   %f\n",rocMean);
//             printf("rms:    %f\n",rocRMS);
//             printf("median: %f\n",rocMedian);
          }
          
          (*treeStreamer)<< "TRDcalibDetails"
            // statistical information
            << (Char_t*)((calibName+"_Mean_SuperModule=").Data()) << smMean[iSM]
            << (Char_t*)((calibName+"_RMS_SuperModule=").Data())  << smRMS[iSM]
            << (Char_t*)((calibName+"_Mean_Chamber=").Data())     << rocMean
            << (Char_t*)((calibName+"_RMS_Chamber=").Data())      << rocRMS
            << (Char_t*)((calibName+"_Median_Chamber=").Data())   << rocMedian
            //pad by pad information
            << (Char_t*)((calibName+".=").Data()) << valueVector;
          
        }   // end loop over calib objects
        
        (*treeStreamer)<< "TRDcalibDetails"
          << "\n";
        arrTrash.Delete();
      }  // end loop over stacks
    }  // end loop over supermodules
  }  // end loop over layers
  delete treeStreamer;
}

//_____________________________________________________________________________
void AliTRDCalibViewer::ProcessTRDCalibArray(AliTRDCalDet* chamberCalib, AliTRDCalPad *padCalib,
					     TString parName,
					     Double_t &runValue, Double_t &runRMS,
					     TVectorD &chamberValues, TVectorD &chamberValuesRMS,
					     TVectorD &superModuleValues, TVectorD &superModuleValuesRMS) {
  // Process the calibrations for a given run.
  // Calculates the run and chamber wise averages
  //
  if(!chamberCalib) return;
  if(!padCalib) return;
  Int_t kSuperModuleStatus[18] = {1, 1, 0, 0, 0, 0,     // super module status (1- installed, 0- not installed)
				  0, 1, 1, 1, 1, 0, 
				  0, 0, 0, 0, 0, 1};

  // initialize the histograms used for extracting the mean and RMS
  TH1F *runWiseHisto = new TH1F("runHisto", "runHisto", 200, -10.0, 10.0);
  TH1F *superModuleWiseHisto = new TH1F("smHisto", "smHisto", 200, -10.0, 10.0);
  TH1F *chamberWiseHisto = new TH1F("chamberHisto", "chamberHisto", 200, -10.0, 10.0);

  // check if the calibration parameter is multiplicative or additive
  Bool_t multiplicative = kTRUE;
  if(!parName.CompareTo("T0")) multiplicative = kFALSE;

  // first iteration (calculate all averages and RMS without discrimination on the SM average)
  runWiseHisto->Reset();
  for(Int_t iSM = 0; iSM<AliTRDcalibDB::kNsector; iSM++) {   // loop over supermodules
    // reset the super module histogram
    superModuleWiseHisto->Reset();
    // check if SM is installed
    if(!kSuperModuleStatus[iSM]) continue;
    for(Int_t iChamber=iSM*AliTRDcalibDB::kNstack*AliTRDcalibDB::kNlayer; 
	iChamber < (iSM+1)*AliTRDcalibDB::kNstack*AliTRDcalibDB::kNlayer; 
	iChamber++) {  // loop over chambers in this supermodule
      // get the chamber value
      Float_t chamberValue = chamberCalib->GetValue(iChamber);
      // get the ROC object
      AliTRDCalROC *chamberROC = padCalib->GetCalROC(iChamber);
      if(!chamberROC) 
	continue;
      chamberWiseHisto->Reset();
      for(Int_t iChannel = 0; iChannel < chamberROC->GetNchannels(); iChannel++){ // loop over channels
	// calculate the calibration parameter for this pad
	Float_t padValue;
	if(multiplicative)
	  padValue = chamberValue * chamberROC->GetValue(iChannel);
	else
	  padValue = chamberValue + chamberROC->GetValue(iChannel);
	// fill the run, SM and chamber wise histograms
	chamberWiseHisto->Fill(padValue);
	// if the parameter is Noise then check if the pad value is not a default one
	// Default value is now 1.2!!!! Check with Raphaelle for more informations
	if(parName.Contains("Noise") &&
	   TMath::Abs(padValue-1.2)<0.00001) continue;
	superModuleWiseHisto->Fill(padValue);
	runWiseHisto->Fill(padValue);
      }  // end loop over channels
      // get the chamber wise mean and RMS
      chamberValues[iChamber] = chamberWiseHisto->GetMean();
      chamberValuesRMS[iChamber] = chamberWiseHisto->GetRMS();
    }  // end loop over chambers
    // SM wise mean and RMS
    superModuleValues[iSM] = superModuleWiseHisto->GetMean();
    superModuleValuesRMS[iSM] = superModuleWiseHisto->GetRMS();
  }  // end loop over supermodules
  // run wise mean and RMS
  runValue = runWiseHisto->GetMean();
  runRMS = runWiseHisto->GetRMS();

  // Noise and Gain are finished processing
  if(parName.Contains("Noise") || parName.Contains("Gain"))
    return;
  // second iteration (calculate SM and run wise averages and RMS for Vdrift and T0)
  // The pads with calib parameter equal to the SM average are discarded (default value)
  runWiseHisto->Reset();
  for(Int_t iSM = 0; iSM<AliTRDcalibDB::kNsector; iSM++) {   // loop over supermodules
    superModuleWiseHisto->Reset();
    // eliminate the uninstalled super modules
    if(!kSuperModuleStatus[iSM]) continue;
    for(Int_t iChamber=iSM*AliTRDcalibDB::kNstack*AliTRDcalibDB::kNlayer; 
	iChamber < (iSM+1)*AliTRDcalibDB::kNstack*AliTRDcalibDB::kNlayer; 
	iChamber++) {  // loop over chambers
      // the chamber value
      Float_t chamberValue = chamberCalib->GetValue(iChamber);
      AliTRDCalROC *chamberROC = padCalib->GetCalROC(iChamber);
      if(!chamberROC) 
	continue;
      
      for(Int_t iChannel = 0; iChannel < chamberROC->GetNchannels(); iChannel++){ // loop over channels in a chamber
	// get the pad value
	Float_t padValue;
	if(multiplicative)
	  padValue = chamberValue * chamberROC->GetValue(iChannel);
	else
	  padValue = chamberValue + chamberROC->GetValue(iChannel);
	// eliminate from the average and RMS calculation all pads which
	// have the calib parameter equal with the SM average
	if((parName.Contains("Vdrift") || parName.Contains("T0")) && 
	   TMath::Abs(padValue-superModuleValues[iSM])<0.00001) continue;
	superModuleWiseHisto->Fill(padValue);
	runWiseHisto->Fill(padValue);
      }   // end loop over channels
    }   // end loop over chambers 
    
    superModuleValues[iSM] = superModuleWiseHisto->GetMean();
    superModuleValuesRMS[iSM] = superModuleWiseHisto->GetRMS();
  }   // end loop over super modules
  runValue = runWiseHisto->GetMean();
  runRMS = runWiseHisto->GetRMS();

  delete chamberWiseHisto;
  delete superModuleWiseHisto;
  delete runWiseHisto;

  return;
}


//_____________________________________________________________________________
void AliTRDCalibViewer::ProcessTRDCalibArray(AliTRDCalPad *padCalib,
                        	             TVectorD &superModuleValues, TVectorD &superModuleValuesRMS)
{
// Process the calibrations for a given run.
// Calculates the run and chamber wise averages
  
  if(!padCalib) return;
  Int_t kSuperModuleStatus[18] = {1, 1, 0, 0, 0, 0,     // super module status (1- installed, 0- not installed)
      0, 1, 1, 1, 1, 0,
      0, 0, 0, 0, 0, 1};

  for(Int_t iSM = 0; iSM<AliTRDcalibDB::kNsector; iSM++) {   // loop over supermodules
    Double_t mean=0;
    Double_t rms=0;
    Double_t count=0;
    if(!kSuperModuleStatus[iSM]) continue;

    // loop over chambers in this supermodule
    const Int_t nchambers=AliTRDcalibDB::kNstack*AliTRDcalibDB::kNlayer;
    for(Int_t iChamber=iSM*nchambers; iChamber < (iSM+1)*nchambers; iChamber++) {
      AliTRDCalROC *chamberROC = padCalib->GetCalROC(iChamber);
      if(!chamberROC)
        continue;

      // loop over channels
      for(Int_t iChannel = 0; iChannel < chamberROC->GetNchannels(); iChannel++){ 
        mean+=chamberROC->GetValue(iChannel);
        rms+=chamberROC->GetValue(iChannel)*chamberROC->GetValue(iChannel);
        ++count;
      }  // end loop over channels
    }  // end loop over chambers

    //calculate mean and rms
    if (count>0){
      mean/=count;
      rms/=count;
      rms=TMath::Sqrt(TMath::Abs(mean*mean-rms));
    }
    // SM wise mean and RMS
    superModuleValues[iSM]    = mean;
    superModuleValuesRMS[iSM] = rms;
  }  // end loop over supermodules
}
