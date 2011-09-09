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
//  Class for viewing/visualizing TPC calibration data                       //
//  base on  TTree functionality for visualization                           //
//                                                                           //
//  Create a list of AliTPCCalPads, arrange them in an TObjArray.            //
//  Pass this TObjArray to MakeTree and create the calibration Tree          //
//  While craating this tree some statistical information are calculated     //
//  Open the viewer with this Tree: AliTPCCalibViewer v("CalibTree.root")    //
//  Have fun!                                                                //
//  EasyDraw("CETmean~-CETmean_mean", "A", "(CETmean~-CETmean_mean)>0")      //
//                                                                           //
//  If you like to click, we recommand you the                               //
//    AliTPCCalibViewerGUI                                                   //
//                                                                           //
//    THE DOCUMENTATION IS STILL NOT COMPLETED !!!!                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//
// ROOT includes 
//

#include <fstream>
#include <iostream>

#include <TFile.h>
#include <TFriendElement.h>
#include <TGraph.h>
#include <TKey.h>
#include <TPad.h>
//#include <TCanvas.h>
#include <TH1.h> 
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjString.h>
//#include <TROOT.h>
#include <TRandom.h>
#include <TString.h>
#include <TStyle.h>
#include <TTreeStream.h>

#include "AliTPCCalibCE.h"
#include "AliMathBase.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalibPedestal.h"
#include "AliTPCCalibPulser.h"

//
// AliRoot includes
//
#include "AliTPCCalibViewer.h"

ClassImp(AliTPCCalibViewer)


AliTPCCalibViewer::AliTPCCalibViewer()
                  :TObject(),
                   fTree(0),
                   fFile(0),
                   fListOfObjectsToBeDeleted(0),
                   fTreeMustBeDeleted(0), 
                   fAbbreviation(0), 
                   fAppendString(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTPCCalibViewer::AliTPCCalibViewer(const AliTPCCalibViewer &c)
                  :TObject(c),
                   fTree(0),
                   fFile(0),
                   fListOfObjectsToBeDeleted(0),
                   fTreeMustBeDeleted(0),
                   fAbbreviation(0), 
                   fAppendString(0)
{
  //
  // dummy AliTPCCalibViewer copy constructor
  // not yet working!!!
  //
  fTree = c.fTree;
  fTreeMustBeDeleted = c.fTreeMustBeDeleted;
  //fFile = new TFile(*(c.fFile));
  fListOfObjectsToBeDeleted = c.fListOfObjectsToBeDeleted;
  fAbbreviation = c.fAbbreviation;
  fAppendString = c.fAppendString;
}

//_____________________________________________________________________________
AliTPCCalibViewer::AliTPCCalibViewer(TTree *const tree)
                  :TObject(),
                   fTree(0),
                   fFile(0),
                   fListOfObjectsToBeDeleted(0),
                   fTreeMustBeDeleted(0),
                   fAbbreviation(0), 
                   fAppendString(0)
{
  //
  // Constructor that initializes the calibration viewer
  //
  fTree = tree;
  fTreeMustBeDeleted = kFALSE;
  fListOfObjectsToBeDeleted = new TObjArray();
  fAbbreviation = "~";
  fAppendString = ".fElements";
}

//_____________________________________________________________________________
AliTPCCalibViewer::AliTPCCalibViewer(const char* fileName, const char* treeName)
                  :TObject(),
                   fTree(0),
                   fFile(0),
                   fListOfObjectsToBeDeleted(0),
                   fTreeMustBeDeleted(0),
                   fAbbreviation(0), 
                   fAppendString(0)
                   
{
   //
   // Constructor to initialize the calibration viewer
   // the file 'fileName' contains the tree 'treeName'
   //
   fFile = new TFile(fileName, "read");
   fTree = (TTree*) fFile->Get(treeName);
   fTreeMustBeDeleted = kTRUE;
   fListOfObjectsToBeDeleted = new TObjArray();
   fAbbreviation = "~";
   fAppendString = ".fElements";
}
                   
//____________________________________________________________________________
AliTPCCalibViewer & AliTPCCalibViewer::operator =(const AliTPCCalibViewer & param)
{
   //
   // assignment operator - dummy
   // not yet working!!!
   //
   fTree = param.fTree;
   fTreeMustBeDeleted = param.fTreeMustBeDeleted;
   //fFile = new TFile(*(param.fFile));
   fListOfObjectsToBeDeleted = param.fListOfObjectsToBeDeleted;
   fAbbreviation = param.fAbbreviation;
   fAppendString = param.fAppendString;
   return (*this);
}

//_____________________________________________________________________________
AliTPCCalibViewer::~AliTPCCalibViewer()
{
   //
   // AliTPCCalibViewer destructor
   // all objects will be deleted, the file will be closed, the pictures will disappear
   //
   if (fTree && fTreeMustBeDeleted) {
      fTree->SetCacheSize(0);
      fTree->Delete();
      //delete fTree;
   }
   if (fFile) {
      fFile->Close();
      fFile = 0;
   }

   for (Int_t i = fListOfObjectsToBeDeleted->GetEntriesFast()-1; i >= 0; i--) {
      //cout << "Index " << i << " trying to delete the following object: " << fListOfObjectsToBeDeleted->At(i)->GetName() << "..."<< endl;
      delete fListOfObjectsToBeDeleted->At(i);
   }
   delete fListOfObjectsToBeDeleted;
}

//_____________________________________________________________________________
void AliTPCCalibViewer::Delete(Option_t* /*option*/) {
   //
   // Should be called from AliTPCCalibViewerGUI class only.
   // If you use Delete() do not call the destructor.
   // All objects (except those contained in fListOfObjectsToBeDeleted) will be deleted, the file will be closed.
   //
   
   if (fTree && fTreeMustBeDeleted) {
      fTree->SetCacheSize(0);
      fTree->Delete();
   }
   delete fFile;
   delete fListOfObjectsToBeDeleted;
}


const char* AliTPCCalibViewer::AddAbbreviations(const Char_t *c, Bool_t printDrawCommand){ 
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
   TString removeString = "!#";  // very unpropable combination of chars
   TString replaceString = "";
   TString searchString = "";
   TString normString = "";
   TObjArray *listOfVariables = GetListOfVariables();
   listOfVariables->Add(new TObjString("channel"));
   listOfVariables->Add(new TObjString("gx"));
   listOfVariables->Add(new TObjString("gy"));
   listOfVariables->Add(new TObjString("lx"));
   listOfVariables->Add(new TObjString("ly"));
   listOfVariables->Add(new TObjString("pad"));
   listOfVariables->Add(new TObjString("row"));
   listOfVariables->Add(new TObjString("rpad"));
   listOfVariables->Add(new TObjString("sector"));
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
      // printf("normLengths[%i] (%s) = %i \n", i,((TObjString*)listOfNormalizationVariables->At(i))->String().Data(), normLengths[i]);
   }
   Int_t *varSort = new Int_t[nVariables];
   TMath::Sort(nVariables, varLengths, varSort, kTRUE);
   Int_t *normSort = new Int_t[nNorm];
   TMath::Sort(nNorm, normLengths, normSort, kTRUE);
   // for (Int_t i = 0; i<nNorm; i++)  printf("normLengths: %i\n", normLengths[normSort[i]]);
   // for (Int_t i = 0; i<nVariables; i++) printf("varLengths: %i\n", varLengths[varSort[i]]);
   
   for (Int_t ivar = 0; ivar < nVariables; ivar++) {
      // ***** save correct tokens *****
      // first get the next variable:
      searchString = ((TObjString*)listOfVariables->At(varSort[ivar]))->String();
      // printf("searchString: %s ++++++++++++++\n", searchString.Data());
      // form replaceString:
      replaceString = "";
      for (Int_t i = 0; i < searchString.Length(); i++) {
         replaceString.Append(searchString[i]);
         if (i == 0) replaceString.Append(removeString);
      }
      // go through normalization:
      // printf("go through normalization\n");
      for (Int_t inorm = 0; inorm < nNorm; inorm++) {
         // printf(" inorm=%i, nNorm=%i, normSort[inorm]=%i \n", inorm, nNorm, normSort[inorm]);
         normString = ((TObjString*)listOfNormalizationVariables->At(normSort[inorm]))->String();
         // printf(" walking in normalization, i=%i, normString=%s \n", inorm, normString.Data());
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
   delete [] varLengths;
   delete [] normLengths;
   delete [] varSort;
   delete [] normSort;
   return str.Data();
}




//_____________________________________________________________________________
Int_t AliTPCCalibViewer::EasyDraw(const char* drawCommand, const char* sector, const char* cuts, const char* drawOptions, Bool_t writeDrawCommand) const {
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
   sectorStr.ToUpper();
   TString cutStr("");
   //TString drawOptionsStr("profcolz ");
   Bool_t dangerousToDraw = drawStr.Contains(":") || drawStr.Contains(">>");
   if (dangerousToDraw) {
      Warning("EasyDraw", "The draw string must not contain ':' or '>>'. Using only first variable for drawing!");
//      return -1;
//      drawStr.Resize(drawStr.First(">"));
      drawStr.Resize(drawStr.First(":"));
   }

   TString drawOptionsStr("");
   TRandom rnd(0);
   Int_t rndNumber = rnd.Integer(10000);

   if (drawOptions && strcmp(drawOptions, "") != 0)
      drawOptionsStr += drawOptions;
   else
      drawOptionsStr += "profcolz";

   if (sectorStr == "A") {
      drawStr += Form(":gy%s:gx%s>>prof", fAppendString.Data(), fAppendString.Data());
      drawStr += rndNumber;
      drawStr += "(330,-250,250,330,-250,250)";
      cutStr += "(sector/18)%2==0 ";
   }
   else if  (sectorStr == "C") {
      drawStr += Form(":gy%s:gx%s>>prof", fAppendString.Data(), fAppendString.Data());
      drawStr += rndNumber;
      drawStr += "(330,-250,250,330,-250,250)";
      cutStr += "(sector/18)%2==1 ";
   }
   else if  (sectorStr == "ALL") {
      drawStr += Form(":gy%s:gx%s>>prof", fAppendString.Data(), fAppendString.Data());
      drawStr += rndNumber;
      drawStr += "(330,-250,250,330,-250,250)";
   }
   else if  (sectorStr.Contains("S")) {
      drawStr += Form(":rpad%s:row%s+(sector>35)*63>>prof", fAppendString.Data(), fAppendString.Data());
      drawStr += rndNumber;
      drawStr += "(159,0,159,140,-70,70)";
      TString sec=sectorStr;
      sec.Remove(0,1);
      cutStr += "sector%36=="+sec+" ";
   }
   else if (sectorStr.IsDigit()) {
      Int_t isec = sectorStr.Atoi();
      drawStr += Form(":rpad%s:row%s>>prof", fAppendString.Data(), fAppendString.Data());
      drawStr += rndNumber;
      if (isec < 36 && isec >= 0)
         drawStr += "(63,0,63,108,-54,54)";
      else if (isec < 72 && isec >= 36)
         drawStr += "(96,0,96,140,-70,70)";
      else {
         Error("EasyDraw","The TPC contains only sectors between 0 and 71.");
         return -1;
      }
      cutStr += "(sector==";
      cutStr += isec;
      cutStr += ") ";
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
   if (obj && obj->InheritsFrom("TH1")) FormatHistoLabels((TH1*)obj);
   return returnValue;
}


Int_t AliTPCCalibViewer::EasyDraw(const char* drawCommand, Int_t sector, const char* cuts, const char* drawOptions, Bool_t writeDrawCommand) const {
  //
  // easy drawing of data, use '~' for abbreviation of '.fElements'
  // example: EasyDraw("CETmean~-CETmean_mean", 34, "(CETmean~-CETmean_mean)>0")
  // sector: sector-number - only the specified sector will be drwawn
  // cuts: specifies cuts
  // drawOptions: draw options like 'same'
  // writeDrawCommand: write the command, that is passed to TTree::Draw
  //
   if (sector >= 0 && sector < 72) {
      return EasyDraw(drawCommand, Form("%i", sector), cuts, drawOptions, writeDrawCommand);
   }
   Error("EasyDraw","The TPC contains only sectors between 0 and 71.");
   return -1;
}


//_____________________________________________________________________________
Int_t AliTPCCalibViewer::EasyDraw1D(const char* drawCommand, const char* sector, const char* cuts, const char* drawOptions, Bool_t writeDrawCommand) const {
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
   TString drawOptionsStr(drawOptions);
   sectorStr.ToUpper();
   TString cutStr("");

   if (sectorStr == "A")
      cutStr += "(sector/18)%2==0 ";
   else if  (sectorStr == "C")
      cutStr += "(sector/18)%2==1 ";
   else if (sectorStr.IsDigit()) {
      Int_t isec = sectorStr.Atoi();
      if (isec < 0 || isec > 71) {
         Error("EasyDraw","The TPC contains only sectors between 0 and 71.");
         return -1;
      }
      cutStr += "(sector==";
      cutStr += isec;
      cutStr += ") ";
   }
   else if  (sectorStr.Contains("S")) {
      TString sec=sectorStr;
      sec.Remove(0,1);
      cutStr += "sector%36=="+sec+" ";
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


Int_t AliTPCCalibViewer::EasyDraw1D(const char* drawCommand, Int_t sector, const char* cuts, const char* drawOptions, Bool_t writeDrawCommand) const {
  //
  // easy drawing of data, use '~' for abbreviation of '.fElements'
  // example: EasyDraw("CETmean~-CETmean_mean", 34, "(CETmean~-CETmean_mean)>0")
  // sector: sector-number - the specified sector will be drwawn
  // cuts: specifies cuts
  // drawOptions: draw options like 'same'
  // writeDrawCommand: write the command, that is passed to TTree::Draw
  //

   if (sector >= 0 && sector < 72) {
      return EasyDraw1D(drawCommand, Form("%i",sector), cuts, drawOptions, writeDrawCommand);
   }
  Error("EasyDraw","The TPC contains only sectors between 0 and 71.");
  return -1;
}


void AliTPCCalibViewer::FormatHistoLabels(TH1 *histo) const {
   // 
   // formats title and axis labels of histo 
   // removes '.fElements'
   // 
   if (!histo) return;
   TString replaceString(fAppendString.Data());
   TString *str = new TString(histo->GetTitle());
   str->ReplaceAll(replaceString, "");
   histo->SetTitle(str->Data());
   delete str;
   if (histo->GetXaxis()) {
      str = new TString(histo->GetXaxis()->GetTitle());
      str->ReplaceAll(replaceString, "");
      histo->GetXaxis()->SetTitle(str->Data());
      delete str;
   }
   if (histo->GetYaxis()) {
      str = new TString(histo->GetYaxis()->GetTitle());
      str->ReplaceAll(replaceString, "");
      histo->GetYaxis()->SetTitle(str->Data());
      delete str;
   }
   if (histo->GetZaxis()) {
      str = new TString(histo->GetZaxis()->GetTitle());
      str->ReplaceAll(replaceString, "");
      histo->GetZaxis()->SetTitle(str->Data());
      delete str;
   }
}


Int_t  AliTPCCalibViewer::DrawHisto1D(const char* drawCommand, Int_t sector, const char* cuts, const char *sigmas, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM) const {
   // 
   // Easy drawing of data, in principle the same as EasyDraw1D
   // Difference: A line for the mean / median / LTM is drawn 
   // in 'sigmas' you can specify in which distance to the mean/median/LTM you want to see a line in sigma-units, separated by ';'
   // example: sigmas = "2; 4; 6;"  at Begin_Latex 2 #sigma End_Latex, Begin_Latex 4 #sigma End_Latex and Begin_Latex 6 #sigma End_Latex  a line is drawn.
   // "plotMean", "plotMedian" and "plotLTM": what kind of lines do you want to see?
   // 
   if (sector >= 0 && sector < 72) {
      return DrawHisto1D(drawCommand, Form("%i", sector), cuts, sigmas, plotMean, plotMedian, plotLTM);
   }
   Error("DrawHisto1D","The TPC contains only sectors between 0 and 71.");
   return -1;
}   


Int_t  AliTPCCalibViewer::DrawHisto1D(const char* drawCommand, const char* sector, const char* cuts, const char *sigmas, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM) const {
   // 
   // Easy drawing of data, in principle the same as EasyDraw1D
   // Difference: A line for the mean / median / LTM is drawn 
   // in 'sigmas' you can specify in which distance to the mean/median/LTM you want to see a line in sigma-units, separated by ';'
   // example: sigmas = "2; 4; 6;"  at Begin_Latex 2 #sigma End_Latex, Begin_Latex 4 #sigma End_Latex and Begin_Latex 6 #sigma End_Latex  a line is drawn.
   // "plotMean", "plotMedian" and "plotLTM": what kind of lines do you want to see?
   // 
   Int_t oldOptStat = gStyle->GetOptStat();
   gStyle->SetOptStat(0000000);
   Double_t ltmFraction = 0.8;
   
   TObjArray *sigmasTokens = TString(sigmas).Tokenize(";");  
   TVectorF nsigma(sigmasTokens->GetEntriesFast());
   for (Int_t i = 0; i < sigmasTokens->GetEntriesFast(); i++) {
      TString str(((TObjString*)sigmasTokens->At(i))->GetString());
      Double_t sig = (str.IsFloat()) ? str.Atof() : 0;
      nsigma[i] = sig;
   }
   
   TString drawStr(drawCommand);
   Bool_t dangerousToDraw = drawStr.Contains(":") || drawStr.Contains(">>");
   if (dangerousToDraw) {
      Warning("DrawHisto1D", "The draw string must not contain ':' or '>>'.");
      return -1;
   }
   drawStr += " >> tempHist";
   Int_t entries = EasyDraw1D(drawStr.Data(), sector, cuts);
   TH1F *htemp = (TH1F*)gDirectory->Get("tempHist");
   // FIXME is this histogram deleted automatically?
   Double_t *values = fTree->GetV1();  // value is the array containing 'entries' numbers
   
   Double_t mean = TMath::Mean(entries, values);
   Double_t median = TMath::Median(entries, values);
   Double_t sigma = TMath::RMS(entries, values);
   Double_t maxY = htemp->GetMaximum();
   
   TLegend * legend = new TLegend(.7,.7, .99, .99, "Statistical information");
   //fListOfObjectsToBeDeleted->Add(legend);

   if (plotMean) {
      // draw Mean
      TLine* line = new TLine(mean, 0, mean, maxY);
      //fListOfObjectsToBeDeleted->Add(line);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      line->SetLineStyle(1);
      line->Draw();
      legend->AddEntry(line, Form("Mean: %f", mean), "l");
      // draw sigma lines
      for (Int_t i = 0; i < nsigma.GetNoElements(); i++) {
         TLine* linePlusSigma = new TLine(mean + nsigma[i] * sigma, 0, mean + nsigma[i] * sigma, maxY);
         //fListOfObjectsToBeDeleted->Add(linePlusSigma);
         linePlusSigma->SetLineColor(kRed);
         linePlusSigma->SetLineStyle(2 + i);
         linePlusSigma->Draw();
         TLine* lineMinusSigma = new TLine(mean - nsigma[i] * sigma, 0, mean - nsigma[i] * sigma, maxY);
         //fListOfObjectsToBeDeleted->Add(lineMinusSigma);
         lineMinusSigma->SetLineColor(kRed);
         lineMinusSigma->SetLineStyle(2 + i);
         lineMinusSigma->Draw();
         legend->AddEntry(lineMinusSigma, Form("%i #sigma = %f",(Int_t)(nsigma[i]), (Float_t)(nsigma[i] * sigma)), "l");
      }
   }
   if (plotMedian) {
      // draw median
      TLine* line = new TLine(median, 0, median, maxY);
      //fListOfObjectsToBeDeleted->Add(line);
      line->SetLineColor(kBlue);
      line->SetLineWidth(2);
      line->SetLineStyle(1);
      line->Draw();
      legend->AddEntry(line, Form("Median: %f", median), "l");
      // draw sigma lines
      for (Int_t i = 0; i < nsigma.GetNoElements(); i++) {
         TLine* linePlusSigma = new TLine(median + nsigma[i] * sigma, 0, median + nsigma[i]*sigma, maxY);
         //fListOfObjectsToBeDeleted->Add(linePlusSigma);
         linePlusSigma->SetLineColor(kBlue);
         linePlusSigma->SetLineStyle(2 + i);
         linePlusSigma->Draw();
         TLine* lineMinusSigma = new TLine(median - nsigma[i] * sigma, 0, median - nsigma[i]*sigma, maxY);
         //fListOfObjectsToBeDeleted->Add(lineMinusSigma);
         lineMinusSigma->SetLineColor(kBlue);
         lineMinusSigma->SetLineStyle(2 + i);
         lineMinusSigma->Draw();
         legend->AddEntry(lineMinusSigma, Form("%i #sigma = %f",(Int_t)(nsigma[i]), (Float_t)(nsigma[i] * sigma)), "l");
      }
   }
   if (plotLTM) {
      // draw LTM
      Double_t ltmRms = 0;
      Double_t ltm = GetLTM(entries, values, &ltmRms, ltmFraction);
      TLine* line = new TLine(ltm, 0, ltm, maxY);
      //fListOfObjectsToBeDeleted->Add(line);
      line->SetLineColor(kGreen+2);
      line->SetLineWidth(2);
      line->SetLineStyle(1);
      line->Draw();
      legend->AddEntry(line, Form("LTM: %f", ltm), "l");
      // draw sigma lines
      for (Int_t i = 0; i < nsigma.GetNoElements(); i++) {
         TLine* linePlusSigma = new TLine(ltm + nsigma[i] * ltmRms, 0, ltm + nsigma[i] * ltmRms, maxY);
         //fListOfObjectsToBeDeleted->Add(linePlusSigma);
         linePlusSigma->SetLineColor(kGreen+2);
         linePlusSigma->SetLineStyle(2+i);
         linePlusSigma->Draw();
   
         TLine* lineMinusSigma = new TLine(ltm - nsigma[i] * ltmRms, 0, ltm - nsigma[i] * ltmRms, maxY);
         //fListOfObjectsToBeDeleted->Add(lineMinusSigma);
         lineMinusSigma->SetLineColor(kGreen+2);
         lineMinusSigma->SetLineStyle(2+i);
         lineMinusSigma->Draw();
         legend->AddEntry(lineMinusSigma, Form("%i #sigma = %f", (Int_t)(nsigma[i]), (Float_t)(nsigma[i] * ltmRms)), "l");
      }
   }
   if (!plotMean && !plotMedian && !plotLTM) return -1;
   legend->Draw();
   gStyle->SetOptStat(oldOptStat);
   return 1;
}


Int_t AliTPCCalibViewer::SigmaCut(const char* drawCommand, Int_t sector, const char* cuts, Float_t sigmaMax, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM, Bool_t pm, const char *sigmas, Float_t sigmaStep) const {
   //
   // Creates a histogram Begin_Latex S(t, #mu, #sigma) End_Latex, where you can see, how much of the data are inside sigma-intervals around the mean value
   // The data of the distribution Begin_Latex f(x, #mu, #sigma) End_Latex are given in 'array', 'n' specifies the length of the array
   // 'mean' and 'sigma' are Begin_Latex #mu End_Latex and  Begin_Latex #sigma End_Latex of the distribution in 'array', to be specified by the user
   // 'nbins': number of bins, 'binLow': first bin, 'binUp': last bin
   // sigmaMax: up to which sigma around the mean/median/LTM the histogram is generated (in units of sigma, Begin_Latex t #sigma End_Latex)
   // sigmaStep: the binsize of the generated histogram
   // Begin_Latex 
   // f(x, #mu, #sigma)     #Rightarrow       S(t, #mu, #sigma) = (#int_{#mu}^{#mu + t #sigma} f(x, #mu, #sigma) dx + #int_{#mu}^{#mu - t #sigma} f(x, #mu, #sigma) dx) / (#int_{-#infty}^{+#infty} f(x, #mu, #sigma) dx)
   // End_Latex
   // 
   //
   // Creates a histogram, where you can see, how much of the data are inside sigma-intervals 
   // around the mean/median/LTM
   // with drawCommand, sector and cuts you specify your input data, see EasyDraw
   // sigmaMax: up to which sigma around the mean/median/LTM the histogram is generated (in units of sigma)
   // sigmaStep: the binsize of the generated histogram
   // plotMean/plotMedian/plotLTM: specifies where to put the center
   //
   if (sector >= 0 && sector < 72) {
      return SigmaCut(drawCommand, Form("%i", sector), cuts, sigmaMax, plotMean, plotMedian, plotLTM, pm, sigmas, sigmaStep);
   }
   Error("SigmaCut","The TPC contains only sectors between 0 and 71.");
   return -1;
}


Int_t AliTPCCalibViewer::SigmaCut(const char* drawCommand, const char* sector, const char* cuts, Float_t sigmaMax, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM, Bool_t pm, const char *sigmas, Float_t sigmaStep) const {
   //
   // Creates a histogram, where you can see, how much of the data are inside sigma-intervals 
   // around the mean/median/LTM
   // with drawCommand, sector and cuts you specify your input data, see EasyDraw
   // sigmaMax: up to which sigma around the mean/median/LTM the histogram is generated (in units of sigma)
   // sigmaStep: the binsize of the generated histogram
   // plotMean/plotMedian/plotLTM: specifies where to put the center
   //
  
   Double_t ltmFraction = 0.8;
   
   TString drawStr(drawCommand);
   Bool_t dangerousToDraw = drawStr.Contains(":") || drawStr.Contains(">>");
   if (dangerousToDraw) {
      Warning("SigmaCut", "The draw string must not contain ':' or '>>'.");
      return -1;
   }
   drawStr += " >> tempHist";
   
   Int_t entries = EasyDraw1D(drawStr.Data(), sector, cuts, "goff");
   TH1F *htemp = (TH1F*)gDirectory->Get("tempHist");
   // FIXME is this histogram deleted automatically?
   Double_t *values = fTree->GetV1();  // value is the array containing 'entries' numbers
   
   Double_t mean = TMath::Mean(entries, values);
   Double_t median = TMath::Median(entries, values);
   Double_t sigma = TMath::RMS(entries, values);
   
   TLegend * legend = new TLegend(.7,.7, .99, .99, "Cumulative");
   //fListOfObjectsToBeDeleted->Add(legend);
   TH1F *cutHistoMean = 0;
   TH1F *cutHistoMedian = 0;
   TH1F *cutHistoLTM = 0;
   
   TObjArray *sigmasTokens = TString(sigmas).Tokenize(";");  
   TVectorF nsigma(sigmasTokens->GetEntriesFast());
   for (Int_t i = 0; i < sigmasTokens->GetEntriesFast(); i++) {
      TString str(((TObjString*)sigmasTokens->At(i))->GetString());
      Double_t sig = (str.IsFloat()) ? str.Atof() : 0;
      nsigma[i] = sig;
   }
  
   if (plotMean) {
      cutHistoMean = AliTPCCalibViewer::SigmaCut(htemp, mean, sigma, sigmaMax, sigmaStep, pm);
      if (cutHistoMean) {
         //fListOfObjectsToBeDeleted->Add(cutHistoMean);
         cutHistoMean->SetLineColor(kRed);
         legend->AddEntry(cutHistoMean, "Mean", "l");
         cutHistoMean->SetTitle(Form("%s, cumulative; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         cutHistoMean->Draw();
         DrawLines(cutHistoMean, nsigma, legend, kRed, pm);
      } // if (cutHistoMean)
       
   }
   if (plotMedian) {
      cutHistoMedian = AliTPCCalibViewer::SigmaCut(htemp, median, sigma, sigmaMax, sigmaStep, pm);
      if (cutHistoMedian) {
         //fListOfObjectsToBeDeleted->Add(cutHistoMedian);
         cutHistoMedian->SetLineColor(kBlue);
         legend->AddEntry(cutHistoMedian, "Median", "l");
         cutHistoMedian->SetTitle(Form("%s, cumulative; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         if (plotMean && cutHistoMean) cutHistoMedian->Draw("same");
            else cutHistoMedian->Draw();
         DrawLines(cutHistoMedian, nsigma, legend, kBlue, pm);
      }  // if (cutHistoMedian)
   }
   if (plotLTM) {
      Double_t ltmRms = 0;
      Double_t ltm = GetLTM(entries, values, &ltmRms, ltmFraction);
      cutHistoLTM = AliTPCCalibViewer::SigmaCut(htemp, ltm, ltmRms, sigmaMax, sigmaStep, pm);
      if (cutHistoLTM) {
         //fListOfObjectsToBeDeleted->Add(cutHistoLTM);
         cutHistoLTM->SetLineColor(kGreen+2);
         legend->AddEntry(cutHistoLTM, "LTM", "l");
         cutHistoLTM->SetTitle(Form("%s, cumulative; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         if ((plotMean && cutHistoMean) || (plotMedian && cutHistoMedian)) cutHistoLTM->Draw("same");
            else cutHistoLTM->Draw();
         DrawLines(cutHistoLTM, nsigma, legend, kGreen+2, pm);
      }
   }
   if (!plotMean && !plotMedian && !plotLTM) return -1;
   legend->Draw();
   return 1;
}


Int_t AliTPCCalibViewer::SigmaCutNew(const char* drawCommand, const char* sector, const char* cuts, Float_t /*sigmaMax*/, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM, Bool_t /*pm*/, const char *sigmas, Float_t /*sigmaStep*/) const {
   //
   // Creates a histogram, where you can see, how much of the data are inside sigma-intervals 
   // around the mean/median/LTM
   // with drawCommand, sector and cuts you specify your input data, see EasyDraw
   // sigmaMax: up to which sigma around the mean/median/LTM the histogram is generated (in units of sigma)
   // sigmaStep: the binsize of the generated histogram
   // plotMean/plotMedian/plotLTM: specifies where to put the center
   //
  
   // Double_t ltmFraction = 0.8;  //unused
   
   TString drawStr(drawCommand);
   drawStr += " >> tempHist";
   
   Int_t entries = EasyDraw1D(drawStr.Data(), sector, cuts, "goff");
   TH1F *htemp = (TH1F*)gDirectory->Get("tempHist");
   TGraph *cutGraphMean   = 0;
   // TGraph *cutGraphMedian = 0;
   // TGraph *cutGraphLTM    = 0;
   Double_t *values = fTree->GetV1();  // value is the array containing 'entries' numbers
   Int_t    *index  = new Int_t[entries];
   Float_t  *xarray = new Float_t[entries];
   Float_t  *yarray = new Float_t[entries];
   TMath::Sort(entries, values, index, kFALSE);
   
   Double_t mean = TMath::Mean(entries, values);
   // Double_t median = TMath::Median(entries, values);
   Double_t sigma = TMath::RMS(entries, values);
   
   TLegend * legend = new TLegend(.7,.7, .99, .99, "Cumulative");
   //fListOfObjectsToBeDeleted->Add(legend);
   
   // parse sigmas string
   TObjArray *sigmasTokens = TString(sigmas).Tokenize(";");  
   TVectorF nsigma(sigmasTokens->GetEntriesFast());
   for (Int_t i = 0; i < sigmasTokens->GetEntriesFast(); i++) {
      TString str(((TObjString*)sigmasTokens->At(i))->GetString());
      Double_t sig = (str.IsFloat()) ? str.Atof() : 0;
      nsigma[i] = sig;
   }
   
   if (plotMean) {
      for (Int_t i = 0; i < entries; i++) {
         xarray[i] = TMath::Abs(values[index[i]] - mean) / sigma; 
         yarray[i] = float(i) / float(entries);
      }
      cutGraphMean = new TGraph(entries, xarray, yarray);
      if (cutGraphMean) {
         //fListOfObjectsToBeDeleted->Add(cutGraphMean);
         cutGraphMean->SetLineColor(kRed);
         legend->AddEntry(cutGraphMean, "Mean", "l");
         cutGraphMean->SetTitle(Form("%s, Cumulative; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         cutGraphMean->Draw("alu");
         DrawLines(cutGraphMean, nsigma, legend, kRed, kTRUE);
      }
   }
  delete [] index;
  delete [] xarray;
  delete [] yarray;
   /*
   if (plotMedian) {
      cutHistoMedian = AliTPCCalibViewer::SigmaCut(htemp, median, sigma, sigmaMax, sigmaStep, pm);
      if (cutHistoMedian) {
         fListOfObjectsToBeDeleted->Add(cutHistoMedian);
         cutHistoMedian->SetLineColor(kBlue);
         legend->AddEntry(cutHistoMedian, "Median", "l");
         cutHistoMedian->SetTitle(Form("%s, cumulative; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         if (plotMean && cutHistoMean) cutHistoMedian->Draw("same");
            else cutHistoMedian->Draw();
         DrawLines(cutHistoMedian, nsigma, legend, kBlue, pm);
      }  // if (cutHistoMedian)
   }
   if (plotLTM) {
      Double_t ltmRms = 0;
      Double_t ltm = GetLTM(entries, values, &ltmRms, ltmFraction);
      cutHistoLTM = AliTPCCalibViewer::SigmaCut(htemp, ltm, ltmRms, sigmaMax, sigmaStep, pm);
      if (cutHistoLTM) {
         fListOfObjectsToBeDeleted->Add(cutHistoLTM);
         cutHistoLTM->SetLineColor(kGreen+2);
         legend->AddEntry(cutHistoLTM, "LTM", "l");
         cutHistoLTM->SetTitle(Form("%s, cumulative; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         if (plotMean && cutHistoMean || plotMedian && cutHistoMedian) cutHistoLTM->Draw("same");
            else cutHistoLTM->Draw();
         DrawLines(cutHistoLTM, nsigma, legend, kGreen+2, pm);
      }
   }*/
   if (!plotMean && !plotMedian && !plotLTM) return -1;
   legend->Draw();
   return 1;
}


Int_t AliTPCCalibViewer::Integrate(const char* drawCommand,       Int_t sector, const char* cuts, Float_t sigmaMax, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM, const char *sigmas, Float_t sigmaStep) const {
   //
   // Creates an integrated histogram Begin_Latex S(t, #mu, #sigma) End_Latex, out of the input distribution distribution Begin_Latex f(x, #mu, #sigma) End_Latex, given in "histogram"   
   // "mean" and "sigma" are Begin_Latex #mu End_Latex and  Begin_Latex #sigma End_Latex of the distribution in "histogram", to be specified by the user
   // sigmaMax: up to which sigma around the mean/median/LTM you want to integrate 
   // if "igma == 0" and "sigmaMax == 0" the whole histogram is integrated
   // "sigmaStep": the binsize of the generated histogram, -1 means, that the maximal reasonable stepsize is used
   // The actual work is done on the array.
   /* Begin_Latex 
         f(x, #mu, #sigma)     #Rightarrow       S(t, #mu, #sigma) = #int_{-#infty}^{#mu + t #sigma} f(x, #mu, #sigma) dx / #int_{-#infty}^{+#infty} f(x, #mu, #sigma) dx
      End_Latex  
   */
   if (sector >= 0 && sector < 72) {
      return Integrate(drawCommand, Form("%i", sector), cuts, sigmaMax, plotMean, plotMedian, plotLTM, sigmas, sigmaStep);
   }
   Error("Integrate","The TPC contains only sectors between 0 and 71.");
   return -1;
   
}


Int_t AliTPCCalibViewer::IntegrateOld(const char* drawCommand, const char* sector, const char* cuts, Float_t sigmaMax, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM, const char *sigmas, Float_t sigmaStep) const {
   //
   // Creates an integrated histogram Begin_Latex S(t, #mu, #sigma) End_Latex, out of the input distribution distribution Begin_Latex f(x, #mu, #sigma) End_Latex, given in "histogram"   
   // "mean" and "sigma" are Begin_Latex #mu End_Latex and  Begin_Latex #sigma End_Latex of the distribution in "histogram", to be specified by the user
   // sigmaMax: up to which sigma around the mean/median/LTM you want to integrate 
   // if "igma == 0" and "sigmaMax == 0" the whole histogram is integrated
   // "sigmaStep": the binsize of the generated histogram, -1 means, that the maximal reasonable stepsize is used
   // The actual work is done on the array.
   /* Begin_Latex 
         f(x, #mu, #sigma)     #Rightarrow       S(t, #mu, #sigma) = #int_{-#infty}^{#mu + t #sigma} f(x, #mu, #sigma) dx / #int_{-#infty}^{+#infty} f(x, #mu, #sigma) dx
      End_Latex  
   */
   
   Double_t ltmFraction = 0.8;
   
   TString drawStr(drawCommand);
   drawStr += " >> tempHist";
   
   Int_t entries = EasyDraw1D(drawStr.Data(), sector, cuts, "goff");
   TH1F *htemp = (TH1F*)gDirectory->Get("tempHist");
   // FIXME is this histogram deleted automatically?
   Double_t *values = fTree->GetV1();  // value is the array containing 'entries' numbers
   
   Double_t mean = TMath::Mean(entries, values);
   Double_t median = TMath::Median(entries, values);
   Double_t sigma = TMath::RMS(entries, values);
    
   TObjArray *sigmasTokens = TString(sigmas).Tokenize(";");  
   TVectorF nsigma(sigmasTokens->GetEntriesFast());
   for (Int_t i = 0; i < sigmasTokens->GetEntriesFast(); i++) {
      TString str(((TObjString*)sigmasTokens->At(i))->GetString());
      Double_t sig = (str.IsFloat()) ? str.Atof() : 0;
      nsigma[i] = sig;
   }
  
   TLegend * legend = new TLegend(.7,.7, .99, .99, "Integrated histogram");
   //fListOfObjectsToBeDeleted->Add(legend);
   TH1F *integralHistoMean = 0;
   TH1F *integralHistoMedian = 0;
   TH1F *integralHistoLTM = 0;
  
   if (plotMean) {
      integralHistoMean = AliTPCCalibViewer::Integrate(htemp, mean, sigma, sigmaMax, sigmaStep);
      if (integralHistoMean) {
         //fListOfObjectsToBeDeleted->Add(integralHistoMean);
         integralHistoMean->SetLineColor(kRed);
         legend->AddEntry(integralHistoMean, "Mean", "l");
         integralHistoMean->SetTitle(Form("%s, integrated; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         integralHistoMean->Draw();
         DrawLines(integralHistoMean, nsigma, legend, kRed, kTRUE);
      }
   }
   if (plotMedian) {
      integralHistoMedian = AliTPCCalibViewer::Integrate(htemp, median, sigma, sigmaMax, sigmaStep);
      if (integralHistoMedian) {
         //fListOfObjectsToBeDeleted->Add(integralHistoMedian);
         integralHistoMedian->SetLineColor(kBlue);
         legend->AddEntry(integralHistoMedian, "Median", "l");
         integralHistoMedian->SetTitle(Form("%s, integrated; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         if (plotMean && integralHistoMean) integralHistoMedian->Draw("same");
            else integralHistoMedian->Draw();
         DrawLines(integralHistoMedian, nsigma, legend, kBlue, kTRUE);
      }
   }
   if (plotLTM) {
      Double_t ltmRms = 0;
      Double_t ltm = GetLTM(entries, values, &ltmRms, ltmFraction);
      integralHistoLTM = AliTPCCalibViewer::Integrate(htemp, ltm, ltmRms, sigmaMax, sigmaStep);
      if (integralHistoLTM) {
         //fListOfObjectsToBeDeleted->Add(integralHistoLTM);
         integralHistoLTM->SetLineColor(kGreen+2);
         legend->AddEntry(integralHistoLTM, "LTM", "l");
         integralHistoLTM->SetTitle(Form("%s, integrated; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         if ((plotMean && integralHistoMean) || (plotMedian && integralHistoMedian)) integralHistoLTM->Draw("same");
            else integralHistoLTM->Draw();
         DrawLines(integralHistoLTM, nsigma, legend, kGreen+2, kTRUE);
      }
   }
   if (!plotMean && !plotMedian && !plotLTM) return -1;
   legend->Draw();
   return 1;
}


Int_t AliTPCCalibViewer::Integrate(const char* drawCommand, const char* sector, const char* cuts, Float_t /*sigmaMax*/, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM, const char *sigmas, Float_t /*sigmaStep*/) const {
   //
   // Creates an integrated histogram Begin_Latex S(t, #mu, #sigma) End_Latex, out of the input distribution distribution Begin_Latex f(x, #mu, #sigma) End_Latex, given in "histogram"   
   // "mean" and "sigma" are Begin_Latex #mu End_Latex and  Begin_Latex #sigma End_Latex of the distribution in "histogram", to be specified by the user
   // sigmaMax: up to which sigma around the mean/median/LTM you want to integrate 
   // if "igma == 0" and "sigmaMax == 0" the whole histogram is integrated
   // "sigmaStep": the binsize of the generated histogram, -1 means, that the maximal reasonable stepsize is used
   // The actual work is done on the array.
   /* Begin_Latex 
         f(x, #mu, #sigma)     #Rightarrow       S(t, #mu, #sigma) = #int_{-#infty}^{#mu + t #sigma} f(x, #mu, #sigma) dx / #int_{-#infty}^{+#infty} f(x, #mu, #sigma) dx
      End_Latex  
   */
   
   Double_t ltmFraction = 0.8;
   
   TString drawStr(drawCommand);
   Bool_t dangerousToDraw = drawStr.Contains(":") || drawStr.Contains(">>");
   if (dangerousToDraw) {
      Warning("Integrate", "The draw string must not contain ':' or '>>'.");
      return -1;
   }
   drawStr += " >> tempHist";
   
   Int_t entries = EasyDraw1D(drawStr.Data(), sector, cuts, "goff");
   TH1F *htemp = (TH1F*)gDirectory->Get("tempHist");
   TGraph *integralGraphMean   = 0;
   TGraph *integralGraphMedian = 0;
   TGraph *integralGraphLTM    = 0;
   Double_t *values = fTree->GetV1();  // value is the array containing 'entries' numbers
   Int_t    *index  = new Int_t[entries];
   Float_t  *xarray = new Float_t[entries];
   Float_t  *yarray = new Float_t[entries];
   TMath::Sort(entries, values, index, kFALSE);
   
   Double_t mean = TMath::Mean(entries, values);
   Double_t median = TMath::Median(entries, values);
   Double_t sigma = TMath::RMS(entries, values);
   
   // parse sigmas string
   TObjArray *sigmasTokens = TString(sigmas).Tokenize(";");  
   TVectorF nsigma(sigmasTokens->GetEntriesFast());
   for (Int_t i = 0; i < sigmasTokens->GetEntriesFast(); i++) {
      TString str(((TObjString*)sigmasTokens->At(i))->GetString());
      Double_t sig = (str.IsFloat()) ? str.Atof() : 0;
      nsigma[i] = sig;
   }
  
   TLegend * legend = new TLegend(.7,.7, .99, .99, "Integrated histogram");
   //fListOfObjectsToBeDeleted->Add(legend);
  
   if (plotMean) {
      for (Int_t i = 0; i < entries; i++) {
         xarray[i] = (values[index[i]] - mean) / sigma; 
         yarray[i] = float(i) / float(entries);
      }
      integralGraphMean = new TGraph(entries, xarray, yarray);
      if (integralGraphMean) {
         //fListOfObjectsToBeDeleted->Add(integralGraphMean);
         integralGraphMean->SetLineColor(kRed);
         legend->AddEntry(integralGraphMean, "Mean", "l");
         integralGraphMean->SetTitle(Form("%s, integrated; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         integralGraphMean->Draw("alu");
         DrawLines(integralGraphMean, nsigma, legend, kRed, kTRUE);
      }
   }
   if (plotMedian) {
      for (Int_t i = 0; i < entries; i++) {
         xarray[i] = (values[index[i]] - median) / sigma; 
         yarray[i] = float(i) / float(entries);
      }
      integralGraphMedian = new TGraph(entries, xarray, yarray);
      if (integralGraphMedian) {
         //fListOfObjectsToBeDeleted->Add(integralGraphMedian);
         integralGraphMedian->SetLineColor(kBlue);
         legend->AddEntry(integralGraphMedian, "Median", "l");
         integralGraphMedian->SetTitle(Form("%s, integrated; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         if (plotMean && integralGraphMean) integralGraphMedian->Draw("samelu");
            else integralGraphMedian->Draw("alu");
         DrawLines(integralGraphMedian, nsigma, legend, kBlue, kTRUE);
      }
   }
   if (plotLTM) {
      Double_t ltmRms = 0;
      Double_t ltm = GetLTM(entries, values, &ltmRms, ltmFraction);
      for (Int_t i = 0; i < entries; i++) {
         xarray[i] = (values[index[i]] - ltm) / ltmRms; 
         yarray[i] = float(i) / float(entries);
      }
      integralGraphLTM = new TGraph(entries, xarray, yarray);
      if (integralGraphLTM) {
         //fListOfObjectsToBeDeleted->Add(integralGraphLTM);
         integralGraphLTM->SetLineColor(kGreen+2);
         legend->AddEntry(integralGraphLTM, "LTM", "l");
         integralGraphLTM->SetTitle(Form("%s, integrated; Multiples of #sigma; Fraction of included data", htemp->GetTitle()));
         if ((plotMean && integralGraphMean) || (plotMedian && integralGraphMedian)) integralGraphLTM->Draw("samelu");
            else integralGraphLTM->Draw("alu");
         DrawLines(integralGraphLTM, nsigma, legend, kGreen+2, kTRUE);
      }
   }
   delete [] index;
   delete [] xarray;
   delete [] yarray;
  
   if (!plotMean && !plotMedian && !plotLTM) return -1;
   legend->Draw();
   return entries;
}


void AliTPCCalibViewer::DrawLines(TH1F *histogram, TVectorF nsigma, TLegend *legend, Int_t color, Bool_t pm) const {
   // 
   // Private function for SigmaCut(...) and Integrate(...)
   // Draws lines into the given histogram, specified by "nsigma", the lines are addeed to the legend
   // 
   
   // start to draw the lines, loop over requested sigmas
   for (Int_t i = 0; i < nsigma.GetNoElements(); i++) {
      if (!pm) { 
         Int_t bin = histogram->GetXaxis()->FindBin(nsigma[i]);
         TLine* lineUp = new TLine(nsigma[i], 0, nsigma[i], histogram->GetBinContent(bin));
         //fListOfObjectsToBeDeleted->Add(lineUp);
         lineUp->SetLineColor(color);
         lineUp->SetLineStyle(2 + i);
         lineUp->Draw();
         TLine* lineLeft = new TLine(nsigma[i], histogram->GetBinContent(bin), 0, histogram->GetBinContent(bin));
         //fListOfObjectsToBeDeleted->Add(lineLeft);
         lineLeft->SetLineColor(color);
         lineLeft->SetLineStyle(2 + i);
         lineLeft->Draw();
         legend->AddEntry(lineLeft, Form("Fraction(%f #sigma) = %f",nsigma[i], histogram->GetBinContent(bin)), "l");
      }
      else { // if (pm)
         Int_t bin = histogram->GetXaxis()->FindBin(nsigma[i]);
         TLine* lineUp1 = new TLine(nsigma[i], 0, nsigma[i], histogram->GetBinContent(bin));
         //fListOfObjectsToBeDeleted->Add(lineUp1);
         lineUp1->SetLineColor(color);
         lineUp1->SetLineStyle(2 + i);
         lineUp1->Draw();
         TLine* lineLeft1 = new TLine(nsigma[i], histogram->GetBinContent(bin), histogram->GetBinLowEdge(0)+histogram->GetBinWidth(0), histogram->GetBinContent(bin));
         //fListOfObjectsToBeDeleted->Add(lineLeft1);
         lineLeft1->SetLineColor(color);
         lineLeft1->SetLineStyle(2 + i);
         lineLeft1->Draw();
         legend->AddEntry(lineLeft1, Form("Fraction(+%f #sigma) = %f",nsigma[i], histogram->GetBinContent(bin)), "l");
         bin = histogram->GetXaxis()->FindBin(-nsigma[i]);
         TLine* lineUp2 = new TLine(-nsigma[i], 0, -nsigma[i], histogram->GetBinContent(bin));
         //fListOfObjectsToBeDeleted->Add(lineUp2);
         lineUp2->SetLineColor(color);
         lineUp2->SetLineStyle(2 + i);
         lineUp2->Draw();
         TLine* lineLeft2 = new TLine(-nsigma[i], histogram->GetBinContent(bin), histogram->GetBinLowEdge(0)+histogram->GetBinWidth(0), histogram->GetBinContent(bin));
         //fListOfObjectsToBeDeleted->Add(lineLeft2);
         lineLeft2->SetLineColor(color);
         lineLeft2->SetLineStyle(2 + i);
         lineLeft2->Draw();
         legend->AddEntry(lineLeft2, Form("Fraction(-%f #sigma) = %f",nsigma[i], histogram->GetBinContent(bin)), "l");
      }
   }  // for (Int_t i = 0; i < nsigma.GetNoElements(); i++)   
}


void AliTPCCalibViewer::DrawLines(TGraph *graph, TVectorF nsigma, TLegend *legend, Int_t color, Bool_t pm) const {
   // 
   // Private function for SigmaCut(...) and Integrate(...)
   // Draws lines into the given histogram, specified by "nsigma", the lines are addeed to the legend
   // 
   
   // start to draw the lines, loop over requested sigmas
   for (Int_t i = 0; i < nsigma.GetNoElements(); i++) {
      if (!pm) { 
         TLine* lineUp = new TLine(nsigma[i], 0, nsigma[i], graph->Eval(nsigma[i]));
         //fListOfObjectsToBeDeleted->Add(lineUp);
         lineUp->SetLineColor(color);
         lineUp->SetLineStyle(2 + i);
         lineUp->Draw();
         TLine* lineLeft = new TLine(nsigma[i], graph->Eval(nsigma[i]), 0, graph->Eval(nsigma[i]));
         //fListOfObjectsToBeDeleted->Add(lineLeft);
         lineLeft->SetLineColor(color);
         lineLeft->SetLineStyle(2 + i);
         lineLeft->Draw();
         legend->AddEntry(lineLeft, Form("Fraction(%f #sigma) = %f",nsigma[i], graph->Eval(nsigma[i])), "l");
      }
      else { // if (pm)
         TLine* lineUp1 = new TLine(nsigma[i], 0, nsigma[i], graph->Eval(nsigma[i]));
         //fListOfObjectsToBeDeleted->Add(lineUp1);
         lineUp1->SetLineColor(color);
         lineUp1->SetLineStyle(2 + i);
         lineUp1->Draw();
         TLine* lineLeft1 = new TLine(nsigma[i], graph->Eval(nsigma[i]), graph->GetHistogram()->GetXaxis()->GetBinLowEdge(0), graph->Eval(nsigma[i]));
         //fListOfObjectsToBeDeleted->Add(lineLeft1);
         lineLeft1->SetLineColor(color);
         lineLeft1->SetLineStyle(2 + i);
         lineLeft1->Draw();
         legend->AddEntry(lineLeft1, Form("Fraction(+%f #sigma) = %f",nsigma[i], graph->Eval(nsigma[i])), "l");
         TLine* lineUp2 = new TLine(-nsigma[i], 0, -nsigma[i], graph->Eval(-nsigma[i]));
         //fListOfObjectsToBeDeleted->Add(lineUp2);
         lineUp2->SetLineColor(color);
         lineUp2->SetLineStyle(2 + i);
         lineUp2->Draw();
         TLine* lineLeft2 = new TLine(-nsigma[i], graph->Eval(-nsigma[i]), graph->GetHistogram()->GetXaxis()->GetBinLowEdge(0), graph->Eval(-nsigma[i]));
         //fListOfObjectsToBeDeleted->Add(lineLeft2);
         lineLeft2->SetLineColor(color);
         lineLeft2->SetLineStyle(2 + i);
         lineLeft2->Draw();
         legend->AddEntry(lineLeft2, Form("Fraction(-%f #sigma) = %f",nsigma[i], graph->Eval(-nsigma[i])), "l");
      }
   }  // for (Int_t i = 0; i < nsigma.GetNoElements(); i++)   
}





/////////////////
// Array tools //
/////////////////


Int_t AliTPCCalibViewer::GetBin(Float_t value, Int_t nbins, Double_t binLow, Double_t binUp){
   // Returns the 'bin' for 'value'
   // The interval between 'binLow' and 'binUp' is divided into 'nbins' equidistant bins
   // avoid index out of bounds error: 'if (bin < binLow) bin = binLow' and vice versa
   /* Begin_Latex
         GetBin(value) = #frac{nbins - 1}{binUp - binLow} #upoint (value - binLow) +1
      End_Latex
   */
   
   Int_t bin =  TMath::Nint( (Float_t)(value - binLow) / (Float_t)(binUp - binLow) * (nbins-1) ) + 1;
   // avoid index out of bounds:   
   if (value < binLow) bin = 0;
   if (value > binUp)  bin = nbins + 1;
   return bin;
   
}   


Double_t AliTPCCalibViewer::GetLTM(Int_t n, const Double_t *const array, Double_t *const sigma, Double_t fraction){
   //
   //  returns the LTM and sigma
   //
   Double_t *ddata = new Double_t[n];
   Double_t mean = 0, lsigma = 0;
   UInt_t nPoints = 0;
   for (UInt_t i = 0; i < (UInt_t)n; i++) {
         ddata[nPoints]= array[nPoints];
         nPoints++;
   }
   Int_t hh = TMath::Min(TMath::Nint(fraction * nPoints), Int_t(n));
   AliMathBase::EvaluateUni(nPoints, ddata, mean, lsigma, hh);
   if (sigma) *sigma = lsigma;
   delete [] ddata;
   return mean;
}


TH1F* AliTPCCalibViewer::SigmaCut(TH1F *const histogram, Float_t mean, Float_t sigma, Float_t sigmaMax, Float_t sigmaStep, Bool_t pm) {
   //
   // Creates a cumulative histogram Begin_Latex S(t, #mu, #sigma) End_Latex, where you can see, how much of the data are inside sigma-intervals around the mean value
   // The data of the distribution Begin_Latex f(x, #mu, #sigma) End_Latex are given in 'histogram'
   // 'mean' and 'sigma' are Begin_Latex #mu End_Latex and  Begin_Latex #sigma End_Latex of the distribution in 'histogram', to be specified by the user
   // sigmaMax: up to which sigma around the mean/median/LTM the histogram is generated (in units of sigma, Begin_Latex t #sigma End_Latex)
   // sigmaStep: the binsize of the generated histogram, -1 means, that the maximal reasonable stepsize is used
   // pm: Decide weather Begin_Latex t > 0 End_Latex (first case) or Begin_Latex t End_Latex arbitrary (secound case)
   // The actual work is done on the array.
   /* Begin_Latex 
         f(x, #mu, #sigma)     #Rightarrow       S(t, #mu, #sigma) = (#int_{#mu}^{#mu + t #sigma} f(x, #mu, #sigma) dx + #int_{#mu}^{#mu - t #sigma} f(x, #mu, #sigma) dx) / (#int_{-#infty}^{+#infty} f(x, #mu, #sigma) dx),    for  t > 0    
         or      
         f(x, #mu, #sigma)     #Rightarrow       S(t, #mu, #sigma) = #int_{#mu}^{#mu + t #sigma} f(x, #mu, #sigma) dx / #int_{-#infty}^{+#infty} f(x, #mu, #sigma) dx
      End_Latex  
      Begin_Macro(source)
      {
         Float_t mean = 0;
         Float_t sigma = 1.5;
         Float_t sigmaMax = 4;
         gROOT->SetStyle("Plain");
         TH1F *distribution = new TH1F("Distrib1", "Distribution f(x, #mu, #sigma)", 1000,-5,5);
         TRandom rand(23);
         for (Int_t i = 0; i <50000;i++) distribution->Fill(rand.Gaus(mean, sigma));
         Float_t *ar = distribution->GetArray();
         
         TCanvas* macro_example_canvas = new TCanvas("cAliTPCCalibViewer1", "", 350, 350);
         macro_example_canvas->Divide(0,3);
         TVirtualPad *pad1 = macro_example_canvas->cd(1);
         pad1->SetGridy();
         pad1->SetGridx();
         distribution->Draw();
         TVirtualPad *pad2 = macro_example_canvas->cd(2);
         pad2->SetGridy();
         pad2->SetGridx();
         
         TH1F *shist = AliTPCCalibViewer::SigmaCut(distribution, mean, sigma, sigmaMax);
         shist->SetNameTitle("Cumulative","Cumulative S(t, #mu, #sigma)");
         shist->Draw();  
         TVirtualPad *pad3 = macro_example_canvas->cd(3);
         pad3->SetGridy();
         pad3->SetGridx();
         TH1F *shistPM = AliTPCCalibViewer::SigmaCut(distribution, mean, sigma, sigmaMax, -1, kTRUE);
         shistPM->Draw();   
         return macro_example_canvas;
      }  
      End_Macro
   */ 
   
   Float_t *array = histogram->GetArray();
   Int_t    nbins = histogram->GetXaxis()->GetNbins();
   Float_t binLow = histogram->GetXaxis()->GetXmin();
   Float_t binUp  = histogram->GetXaxis()->GetXmax();
   return AliTPCCalibViewer::SigmaCut(nbins, array, mean, sigma, nbins, binLow, binUp, sigmaMax, sigmaStep, pm);
}   
   

TH1F* AliTPCCalibViewer::SigmaCut(Int_t n, const Float_t *array, Float_t mean, Float_t sigma, Int_t nbins, Float_t binLow, Float_t binUp, Float_t sigmaMax, Float_t sigmaStep, Bool_t pm){
   //
   // Creates a histogram Begin_Latex S(t, #mu, #sigma) End_Latex, where you can see, how much of the data are inside sigma-intervals around the mean value
   // The data of the distribution Begin_Latex f(x, #mu, #sigma) End_Latex are given in 'array', 'n' specifies the length of the array
   // 'mean' and 'sigma' are Begin_Latex #mu End_Latex and  Begin_Latex #sigma End_Latex of the distribution in 'array', to be specified by the user
   // 'nbins': number of bins, 'binLow': first bin, 'binUp': last bin
   // sigmaMax: up to which sigma around the mean/median/LTM the histogram is generated (in units of sigma, Begin_Latex t #sigma End_Latex)
   // sigmaStep: the binsize of the generated histogram
   // Here the actual work is done.
   
   if (sigma == 0) return 0;
   Float_t binWidth = (binUp-binLow)/(nbins - 1);
   if (sigmaStep <= 0) sigmaStep = binWidth;
   Int_t kbins = (Int_t)(sigmaMax * sigma / sigmaStep) + 1; // + 1  due to overflow bin in histograms
   if (pm) kbins = 2 * (Int_t)(sigmaMax * sigma / sigmaStep) + 1;
   Float_t kbinLow = !pm ? 0 : -sigmaMax;
   Float_t kbinUp  = sigmaMax;
   TH1F *hist = new TH1F("sigmaCutHisto","Cumulative; Multiples of #sigma; Fraction of included data", kbins, kbinLow, kbinUp); 
   hist->SetDirectory(0);
   hist->Reset();
   
   // calculate normalization
   Double_t normalization = 0;
   for (Int_t i = 0; i <= n; i++) {
        normalization += array[i];
   }
   
   // given units: units from given histogram
   // sigma units: in units of sigma
   // iDelta: integrate in interval (mean +- iDelta), given units
   // x:      ofset from mean for integration, given units
   // hist:   needs 
   
//    printf("nbins: %i, binLow: %f, binUp: %f \n", nbins, binLow, binUp);
   // fill histogram
   for (Float_t iDelta = 0; iDelta <= sigmaMax * sigma; iDelta += sigmaStep) {
      // integrate array
      Double_t valueP = array[GetBin(mean, nbins, binLow, binUp)];
      Double_t valueM = array[GetBin(mean-binWidth, nbins, binLow, binUp)];
      // add bin of mean value only once to the histogram
//       printf("++ adding bins: ");
      for (Float_t x = binWidth; x <= iDelta; x += binWidth) {
         valueP += (mean + x <= binUp)  ? array[GetBin(mean + x, nbins, binLow, binUp)] : 0;
         valueM += (mean-binWidth - x >= binLow) ? array[GetBin(mean-binWidth - x, nbins, binLow, binUp)] : 0; 
//          printf("%i, ", GetBin(mean + x, nbins, binLow, binUp));        
      }
//       printf("\n");
      if (valueP / normalization > 100) printf("+++ Error, value to big: %f, normalization with %f will fail  +++ \n", valueP, normalization);
      if (valueP / normalization > 100) return hist;
      if (valueM / normalization > 100) printf("+++ Error, value to big: %f, normalization with %f will fail  +++ \n", valueM, normalization);
      if (valueM / normalization > 100) return hist;
      valueP = (valueP / normalization);
      valueM = (valueM / normalization);
      if (pm) {
         Int_t bin = GetBin(iDelta/sigma, kbins, kbinLow, kbinUp);
         hist->SetBinContent(bin, valueP);
         bin = GetBin(-iDelta/sigma, kbins, kbinLow, kbinUp);
         hist->SetBinContent(bin, valueM);
      }
      else { // if (!pm)
         Int_t bin = GetBin(iDelta/sigma, kbins, kbinLow, kbinUp);
         hist->SetBinContent(bin, valueP + valueM);
//          printf("  first integration bin: %i, last integration bin in + direction: %i \n", GetBin(mean+binWidth, nbins, binLow, binUp), GetBin(iDelta, nbins, binLow, binUp));
//          printf("  first integration bin: %i, last integration bin in - direction: %i \n", GetBin(mean+binWidth, nbins, binLow, binUp), GetBin(-iDelta, nbins, binLow, binUp));
//          printf("  value: %f, normalization: %f, iDelta: %f, Bin: %i \n", valueP+valueM, normalization, iDelta, bin);
      }
   }
   //hist->SetMaximum(0.7);
   if (!pm) hist->SetMaximum(1.2);
   return hist;
}


TH1F* AliTPCCalibViewer::SigmaCut(Int_t /*n*/, const Double_t */*array*/, Double_t /*mean*/, Double_t /*sigma*/, Int_t /*nbins*/, const Double_t */*xbins*/, Double_t /*sigmaMax*/){
   // 
   // SigmaCut for variable binsize
   // NOT YET IMPLEMENTED !!!
   // 
   printf("SigmaCut with variable binsize, Not yet implemented\n");
   
   return 0;
}   


TH1F* AliTPCCalibViewer::Integrate(TH1F *const histogram, Float_t mean, Float_t sigma, Float_t sigmaMax, Float_t sigmaStep){
   //
   // Creates an integrated histogram Begin_Latex S(t, #mu, #sigma) End_Latex, out of the input distribution distribution Begin_Latex f(x, #mu, #sigma) End_Latex, given in "histogram"   
   // "mean" and "sigma" are Begin_Latex #mu End_Latex and  Begin_Latex #sigma End_Latex of the distribution in "histogram", to be specified by the user
   // sigmaMax: up to which sigma around the mean/median/LTM you want to integrate 
   // if "igma == 0" and "sigmaMax == 0" the whole histogram is integrated
   // "sigmaStep": the binsize of the generated histogram, -1 means, that the maximal reasonable stepsize is used
   // The actual work is done on the array.
   /* Begin_Latex 
         f(x, #mu, #sigma)     #Rightarrow       S(t, #mu, #sigma) = #int_{-#infty}^{#mu + t #sigma} f(x, #mu, #sigma) dx / #int_{-#infty}^{+#infty} f(x, #mu, #sigma) dx
      End_Latex  
      Begin_Macro(source)
      {
         Float_t mean = 0;
         Float_t sigma = 1.5;
         Float_t sigmaMax = 4;
         gROOT->SetStyle("Plain");
         TH1F *distribution = new TH1F("Distrib2", "Distribution f(x, #mu, #sigma)", 1000,-5,5);
         TRandom rand(23);
         for (Int_t i = 0; i <50000;i++) distribution->Fill(rand.Gaus(mean, sigma));
         Float_t *ar = distribution->GetArray();
         
         TCanvas* macro_example_canvas = new TCanvas("cAliTPCCalibViewer2", "", 350, 350);
         macro_example_canvas->Divide(0,2);
         TVirtualPad *pad1 = macro_example_canvas->cd(1);
         pad1->SetGridy();
         pad1->SetGridx();
         distribution->Draw();
         TVirtualPad *pad2 = macro_example_canvas->cd(2);
         pad2->SetGridy();
         pad2->SetGridx();
         TH1F *shist = AliTPCCalibViewer::Integrate(distribution, mean, sigma, sigmaMax);
         shist->SetNameTitle("Cumulative","Cumulative S(t, #mu, #sigma)");
         shist->Draw();  
         
         return macro_example_canvas_Integrate;
      }  
      End_Macro
   */ 

   
   Float_t *array = histogram->GetArray();
   Int_t    nbins = histogram->GetXaxis()->GetNbins();
   Float_t binLow = histogram->GetXaxis()->GetXmin();
   Float_t binUp  = histogram->GetXaxis()->GetXmax();
   return AliTPCCalibViewer::Integrate(nbins, array, nbins, binLow, binUp, mean, sigma, sigmaMax, sigmaStep);
}   


TH1F* AliTPCCalibViewer::Integrate(Int_t n, const Float_t *const array, Int_t nbins, Float_t binLow, Float_t binUp, Float_t mean, Float_t sigma, Float_t sigmaMax, Float_t sigmaStep){
   // Creates an integrated histogram Begin_Latex S(t, #mu, #sigma) End_Latex, out of the input distribution distribution Begin_Latex f(x, #mu, #sigma) End_Latex, given in "histogram"   
   // "mean" and "sigma" are Begin_Latex #mu End_Latex and  Begin_Latex #sigma End_Latex of the distribution in "histogram", to be specified by the user
   // sigmaMax: up to which sigma around the mean/median/LTM you want to integrate 
   // if "igma == 0" and "sigmaMax == 0" the whole histogram is integrated
   // "sigmaStep": the binsize of the generated histogram, -1 means, that the maximal reasonable stepsize is used
   // Here the actual work is done.
      
   Bool_t givenUnits = kTRUE;
   if (sigma != 0 && sigmaMax != 0) givenUnits = kFALSE;
   if (givenUnits) {
      sigma = 1;
      sigmaMax = (binUp - binLow) / 2.;
   }
   
   Float_t binWidth = (binUp-binLow)/(nbins - 1);
   if (sigmaStep <= 0) sigmaStep = binWidth;
   Int_t kbins =  (Int_t)(sigmaMax * sigma / sigmaStep) + 1;  // + 1  due to overflow bin in histograms
   Float_t kbinLow = givenUnits ? binLow : -sigmaMax;
   Float_t kbinUp  = givenUnits ? binUp  : sigmaMax;
   TH1F *hist = 0; 
   if (givenUnits)  hist = new TH1F("integratedHisto","Integrated Histogram; Given x; Fraction of included data", kbins, kbinLow, kbinUp); 
   if (!givenUnits) hist = new TH1F("integratedHisto","Integrated Histogram; Multiples of #sigma; Fraction of included data", kbins, kbinLow, kbinUp); 
   hist->SetDirectory(0);
   hist->Reset();
   
   // calculate normalization
 //  printf("calculating normalization, integrating from bin 1 to %i \n", n);
   Double_t normalization = 0;
   for (Int_t i = 1; i <= n; i++) {
        normalization += array[i];
   }
 //  printf("normalization: %f \n", normalization);
   
   // given units: units from given histogram
   // sigma units: in units of sigma
   // iDelta: integrate in interval (mean +- iDelta), given units
   // x:      ofset from mean for integration, given units
   // hist:   needs 
   
   // fill histogram
   for (Float_t iDelta = mean - sigmaMax * sigma; iDelta <= mean + sigmaMax * sigma; iDelta += sigmaStep) {
      // integrate array
      Double_t value = 0;
      for (Float_t x = mean - sigmaMax * sigma; x <= iDelta; x += binWidth) {
         value += (x <= binUp && x >= binLow)  ? array[GetBin(x, nbins, binLow, binUp)] : 0;
      }
      if (value / normalization > 100) printf("+++ Error, value to big: %f, normalization with %f will fail  +++ \n", value, normalization);
      if (value / normalization > 100) return hist;
      Int_t bin = GetBin(iDelta/sigma, kbins, kbinLow, kbinUp);
    //  printf("first integration bin: %i, last integration bin: %i \n", GetBin(mean - sigmaMax * sigma, nbins, binLow, binUp), GetBin(iDelta, nbins, binLow, binUp));
    //  printf("value: %f, normalization: %f, normalized value: %f, iDelta: %f, Bin: %i \n", value, normalization, value/normalization, iDelta, bin);
      value = (value / normalization);
      hist->SetBinContent(bin, value);
   }
   return hist;
}





////////////////////////
// end of Array tools //
////////////////////////



//_____________________________________________________________________________
AliTPCCalPad* AliTPCCalibViewer::GetCalPadOld(const char* desiredData, const char* cuts, const char* calPadName) const {
  //
  // creates a AliTPCCalPad out of the 'desiredData'
  // the functionality of EasyDraw1D is used
  // calPadName specifies the name of the created AliTPCCalPad
  //  - this takes a while -
  //
   TString drawStr(desiredData);
   drawStr.Append(":channel");
   drawStr.Append(fAbbreviation);
   AliTPCCalPad * createdCalPad = new AliTPCCalPad(calPadName, calPadName);
   Int_t entries = 0;
   for (Int_t sec = 0; sec < 72; sec++) {
     AliTPCCalROC * roc = createdCalPad->GetCalROC(sec);
      entries = EasyDraw1D(drawStr.Data(), (Int_t)sec, cuts, "goff");
      if (entries == -1) return 0;
      const Double_t *pchannel = fTree->GetV2();
      const Double_t *pvalue   = fTree->GetV1();
      for (Int_t i = 0; i < entries; i++) 
         roc->SetValue((UInt_t)(pchannel[i]), (Float_t)(pvalue[i]));
   }
   return createdCalPad;   
}


//_____________________________________________________________________________
AliTPCCalPad* AliTPCCalibViewer::GetCalPad(const char* desiredData, const char* cuts, const char* calPadName) const {
  //
  // creates a AliTPCCalPad out of the 'desiredData'
  // the functionality of EasyDraw1D is used
  // calPadName specifies the name of the created AliTPCCalPad
  //  - this takes a while -
  //
   TString drawStr(desiredData);
   drawStr.Append(":channel.fElements:sector");
   AliTPCCalPad * createdCalPad = new AliTPCCalPad(calPadName, calPadName);
   //
   Int_t entries = fTree->Draw(drawStr, cuts,"goff");
   const Double_t *pvalue   = fTree->GetV1();
   const Double_t *pchannel = fTree->GetV2();
   const Double_t *psector  = fTree->GetV3();

   for (Int_t ientry=0; ientry<entries; ientry++){
     Int_t sector= TMath::Nint(psector[ientry]);
     AliTPCCalROC * roc = createdCalPad->GetCalROC(sector);
     if (roc) roc->SetValue((UInt_t)(pchannel[ientry]), (Float_t)(pvalue[ientry]));
   }

  //  for (Int_t sec = 0; sec < 72; sec++) {
//      AliTPCCalROC * roc = createdCalPad->GetCalROC(sec);
//       entries = EasyDraw1D(drawStr.Data(), (Int_t)sec, cuts, "goff");
//       if (entries == -1) return 0;
//       for (Int_t i = 0; i < entries; i++) 
//          roc->SetValue((UInt_t)(pchannel[i]), (Float_t)(pvalue[i]));
//    }
   return createdCalPad;   
}

//_____________________________________________________________________________
AliTPCCalROC* AliTPCCalibViewer::GetCalROC(const char* desiredData, UInt_t sector, const char* cuts) const {
  //
  // creates a AliTPCCalROC out of the desiredData
  // the functionality of EasyDraw1D is used
  // sector specifies the sector of the created AliTPCCalROC
  //
   TString drawStr(desiredData);
   drawStr.Append(":channel");
   drawStr.Append(fAbbreviation);
   Int_t entries = EasyDraw1D(drawStr.Data(), (Int_t)sector, cuts, "goff");
   if (entries == -1) return 0;
   AliTPCCalROC * createdROC = new AliTPCCalROC(sector);
   for (Int_t i = 0; i < entries; i++) 
      createdROC->SetValue((UInt_t)(fTree->GetV2()[i]), fTree->GetV1()[i]);
   return createdROC;
}


TObjArray* AliTPCCalibViewer::GetListOfVariables(Bool_t printList) {
  //
  // scan the tree  - produces a list of available variables in the tree
  // printList: print the list to the screen, after the scan is done
  //
   TObjArray* arr = new TObjArray();
   TObjString* str = 0;
   if (!fTree) return 0;
   Int_t nentries = fTree->GetListOfBranches()->GetEntries();
   for (Int_t i = 0; i < nentries; i++) {
      str = new TObjString(fTree->GetListOfBranches()->At(i)->GetName());
      str->String().ReplaceAll("_Median", "");
      str->String().ReplaceAll("_Mean", "");
      str->String().ReplaceAll("_RMS", "");
      str->String().ReplaceAll("_LTM", "");
      str->String().ReplaceAll("_OutlierCutted", "");
      str->String().ReplaceAll(".", "");
      if (!arr->FindObject(str) && 
          !(str->String() == "channel" || str->String() == "gx" || str->String() == "gy" || 
            str->String() == "lx" || str->String() == "ly" || str->String() == "pad" || 
            str->String() == "row" || str->String() == "rpad" || str->String() == "sector"  ))
         arr->Add(str);
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
            str->String().ReplaceAll("_Median", "");
            str->String().ReplaceAll("_Mean", "");
            str->String().ReplaceAll("_RMS", "");
            str->String().ReplaceAll("_LTM", "");
            str->String().ReplaceAll("_OutlierCutted", "");
            str->String().ReplaceAll(".", "");
            if (!(str->String() == "channel" || str->String() == "gx" || str->String() == "gy" || 
                  str->String() == "lx" || str->String() == "ly" || str->String() == "pad" || 
                  str->String() == "row" || str->String() == "rpad" || str->String() == "sector"  )){
               // insert "<friendName>." at the beginning: (<friendName> is per default "R")
               str->String().Insert(0, ".");
               str->String().Insert(0, friendElement->GetName());
               if (!arr->FindObject(str)) arr->Add(str);
               // printf("added string %s \n", str->String().Data());
            }
         }
      }
   } // if (fTree->GetListOfFriends())
   
   arr->Sort();
//   ((TFriendElement*)gui->GetViewer()->GetTree()->GetListOfFriends()->At(0))->GetTree()->GetListOfBranches()->At(0)->GetName()
// ((TFriendElement*)gui->GetViewer()->GetTree()->GetListOfFriends()->At(0))->GetTree()->GetListOfBranches()


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


TObjArray* AliTPCCalibViewer::GetListOfNormalizationVariables(Bool_t printList) const{
  //
  // produces a list of available variables for normalization in the tree
  // printList: print the list to the screen, after the scan is done
  //
   TObjArray* arr = new TObjArray();
   arr->Add(new TObjString("_Mean"));
   arr->Add(new TObjString("_Mean_OutlierCutted"));
   arr->Add(new TObjString("_Median"));
   arr->Add(new TObjString("_Median_OutlierCutted"));
   arr->Add(new TObjString("_LTM"));
   arr->Add(new TObjString("_LTM_OutlierCutted"));
   arr->Add(new TObjString(Form("LFitIntern_4_8%s", fAppendString.Data())));
   arr->Add(new TObjString(Form("GFitIntern_Lin%s", fAppendString.Data())));
   arr->Add(new TObjString(Form("GFitIntern_Par%s", fAppendString.Data())));
   arr->Add(new TObjString("FitLinLocal"));
   arr->Add(new TObjString("FitLinGlobal"));
   arr->Add(new TObjString("FitParLocal"));
   arr->Add(new TObjString("FitParGlobal"));

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


TFriendElement* AliTPCCalibViewer::AddReferenceTree(const char* filename, const char* treename, const char* refname){
  //
  // add a reference tree to the current tree
  // by default the treename is 'calPads' and the reference treename is 'R'
  //
   TFile *file = new TFile(filename);
   fListOfObjectsToBeDeleted->Add(file);
   TTree * tree = (TTree*)file->Get(treename);
   return AddFriend(tree, refname);
}


TObjArray* AliTPCCalibViewer::GetArrayOfCalPads(){
  //
  // Returns a TObjArray with all AliTPCCalPads that are stored in the tree
  //  - this takes a while - 
  //
   TObjArray *listOfCalPads = GetListOfVariables();
   TObjArray *calPadsArray = new TObjArray();
   Int_t numberOfCalPads = listOfCalPads->GetEntries();
   for (Int_t i = 0; i < numberOfCalPads; i++) {
      std::cout << "Creating calPad " << (i+1) << " of " << numberOfCalPads << "\r" << std::flush;
      char* calPadName = (char*)((TObjString*)(listOfCalPads->At(i)))->GetString().Data();
      TString drawCommand = ((TObjString*)(listOfCalPads->At(i)))->GetString();
      drawCommand.Append(fAbbreviation.Data());
      AliTPCCalPad* calPad = GetCalPad(drawCommand.Data(), "", calPadName); 
      calPadsArray->Add(calPad); 
   }
   std::cout << std::endl;
   listOfCalPads->Delete();
   delete listOfCalPads;
   return calPadsArray;
}


TString* AliTPCCalibViewer::Fit(const char* drawCommand, const char* formula, const char* cuts, Double_t & chi2, TVectorD &fitParam, TMatrixD &covMatrix){
   //
   // fit an arbitrary function, specified by formula into the data, specified by drawCommand and cuts
   // returns chi2, fitParam and covMatrix
   // returns TString with fitted formula
   //
   
   TString formulaStr(formula); 
   TString drawStr(drawCommand);
   TString cutStr(cuts);
   
   // abbreviations:
   drawStr.ReplaceAll(fAbbreviation, fAppendString);
   cutStr.ReplaceAll(fAbbreviation, fAppendString);
   formulaStr.ReplaceAll(fAbbreviation, fAppendString);
   
   formulaStr.ReplaceAll("++", fAbbreviation);
   TObjArray* formulaTokens = formulaStr.Tokenize(fAbbreviation.Data()); 
   Int_t dim = formulaTokens->GetEntriesFast();
   
   fitParam.ResizeTo(dim);
   covMatrix.ResizeTo(dim,dim);
   
   TLinearFitter* fitter = new TLinearFitter(dim+1, Form("hyp%d",dim));
   fitter->StoreData(kTRUE);   
   fitter->ClearPoints();
   
   Int_t entries = Draw(drawStr.Data(), cutStr.Data(), "goff");
   if (entries == -1) return new TString("An ERROR has occured during fitting!");
   Double_t **values = new Double_t*[dim+1] ; 
   
   for (Int_t i = 0; i < dim + 1; i++){
      Int_t centries = 0;
      if (i < dim) centries = fTree->Draw(((TObjString*)formulaTokens->At(i))->GetName(), cutStr.Data(), "goff");
      else  centries = fTree->Draw(drawStr.Data(), cutStr.Data(), "goff");
      
      if (entries != centries) {
        delete [] values;
        return new TString("An ERROR has occured during fitting!");
      }
      values[i] = new Double_t[entries];
      memcpy(values[i],  fTree->GetV1(), entries*sizeof(Double_t)); 
   }
   
   // add points to the fitter
   for (Int_t i = 0; i < entries; i++){
      Double_t x[1000];
      for (Int_t j=0; j<dim;j++) x[j]=values[j][i];
      fitter->AddPoint(x, values[dim][i], 1);
   }

   fitter->Eval();
   fitter->GetParameters(fitParam);
   fitter->GetCovarianceMatrix(covMatrix);
   chi2 = fitter->GetChisquare();
//    chi2 = chi2;
   
   TString *preturnFormula = new TString(Form("( %e+",fitParam[0])), &returnFormula = *preturnFormula;
   
   for (Int_t iparam = 0; iparam < dim; iparam++) {
     returnFormula.Append(Form("%s*(%e)",((TObjString*)formulaTokens->At(iparam))->GetName(),fitParam[iparam+1]));
     if (iparam < dim-1) returnFormula.Append("+");
   }
   returnFormula.Append(" )");
   delete formulaTokens;
   delete fitter;
   for (Int_t i = 0; i < dim + 1; i++){
     delete [] values[i];
   }
   delete [] values;
   return preturnFormula;
}


void AliTPCCalibViewer::MakeTreeWithObjects(const char *fileName, const TObjArray *const array, const char * mapFileName) {
  //
  // Write tree with all available information
  // im mapFileName is speciefied, the Map information are also written to the tree
  // AliTPCCalPad-Objects are written directly to the tree, so that they can be accessd later on
  // (does not work!!!)
  //
   AliTPCROC* tpcROCinstance = AliTPCROC::Instance();

   TObjArray* mapIROCs = 0;
   TObjArray* mapOROCs = 0;
   TVectorF *mapIROCArray = 0;
   TVectorF *mapOROCArray = 0;
   Int_t mapEntries = 0;
   TString* mapNames = 0;
   
   if (mapFileName) {
      TFile mapFile(mapFileName, "read");
      
      TList* listOfROCs = mapFile.GetListOfKeys();
      mapEntries = listOfROCs->GetEntries()/2;
      mapIROCs = new TObjArray(mapEntries*2);
      mapOROCs = new TObjArray(mapEntries*2);
      mapIROCArray = new TVectorF[mapEntries];
      mapOROCArray = new TVectorF[mapEntries];
      
      mapNames = new TString[mapEntries];
      for (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
         TString rocName(((TKey*)(listOfROCs->At(ivalue*2)))->GetName());
         rocName.Remove(rocName.Length()-4, 4);
         mapIROCs->AddAt((AliTPCCalROC*)mapFile.Get((rocName + "IROC").Data()), ivalue);
         mapOROCs->AddAt((AliTPCCalROC*)mapFile.Get((rocName + "OROC").Data()), ivalue);
         mapNames[ivalue].Append(rocName);
      }
      
      for (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
         mapIROCArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(0));
         mapOROCArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(36));
      
         for (UInt_t ichannel = 0; ichannel < tpcROCinstance->GetNChannels(0); ichannel++)
            (mapIROCArray[ivalue])[ichannel] = ((AliTPCCalROC*)(mapIROCs->At(ivalue)))->GetValue(ichannel);
         for (UInt_t ichannel = 0; ichannel < tpcROCinstance->GetNChannels(36); ichannel++)
            (mapOROCArray[ivalue])[ichannel] = ((AliTPCCalROC*)(mapOROCs->At(ivalue)))->GetValue(ichannel);
      }

   } //  if (mapFileName)
  
   TTreeSRedirector cstream(fileName);
   Int_t arrayEntries = array->GetEntries();
   
   // Read names of AliTPCCalPads and save them in names[]
   TString* names = new TString[arrayEntries];
   for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++)
      names[ivalue].Append(((AliTPCCalPad*)array->At(ivalue))->GetName());

   for (UInt_t isector = 0; isector < tpcROCinstance->GetNSectors(); isector++) {
      
      TVectorF *vectorArray = new TVectorF[arrayEntries];
      for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++)
         vectorArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));
            
      
      //
      // fill vectors of variable per pad
      //
      TVectorF *posArray = new TVectorF[8];
      for (Int_t ivalue = 0; ivalue < 8; ivalue++)
         posArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));

      Float_t posG[3] = {0};
      Float_t posL[3] = {0};
      Int_t ichannel = 0;
      for (UInt_t irow = 0; irow < tpcROCinstance->GetNRows(isector); irow++) {
         for (UInt_t ipad = 0; ipad < tpcROCinstance->GetNPads(isector, irow); ipad++) {
            tpcROCinstance->GetPositionLocal(isector, irow, ipad, posL);
            tpcROCinstance->GetPositionGlobal(isector, irow, ipad, posG);
            posArray[0][ichannel] = irow;
            posArray[1][ichannel] = ipad;
            posArray[2][ichannel] = posL[0];
            posArray[3][ichannel] = posL[1];
            posArray[4][ichannel] = posG[0];
            posArray[5][ichannel] = posG[1];
            posArray[6][ichannel] = (Int_t)(ipad - (Double_t)(tpcROCinstance->GetNPads(isector, irow))/2);
            posArray[7][ichannel] = ichannel;
            
            // loop over array containing AliTPCCalPads
            for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
               AliTPCCalPad* calPad = (AliTPCCalPad*) array->At(ivalue);
               AliTPCCalROC* calROC = calPad->GetCalROC(isector);
               if (calROC)
                  (vectorArray[ivalue])[ichannel] = calROC->GetValue(irow, ipad);
               else
                  (vectorArray[ivalue])[ichannel] = 0;
            }
            ichannel++;
         }
      }
      AliTPCCalROC dummyROC(0);
      for  (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
         AliTPCCalROC *roc = ((AliTPCCalPad*)array->At(ivalue))->GetCalROC(isector);
         if (!roc) roc = &dummyROC;
         cstream << "calPads" <<
            (Char_t*)((names[ivalue] + ".=").Data()) << &vectorArray[ivalue];
         cstream << "calPads" << 
            (Char_t*)((names[ivalue] + "Pad.=").Data()) << roc;
      }

      if (mapFileName) {
         for  (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
            if (isector < 36)
               cstream << "calPads" <<
                  (Char_t*)((mapNames[ivalue] + ".=").Data()) << &mapIROCArray[ivalue];
            else
               cstream << "calPads" <<
                  (Char_t*)((mapNames[ivalue] + ".=").Data()) << &mapOROCArray[ivalue];
         }
      }
      
      cstream << "calPads" <<
         "sector=" << isector;

      cstream << "calPads" <<
         "row.=" << &posArray[0] <<
         "pad.=" << &posArray[1] <<
         "lx.=" << &posArray[2] <<
         "ly.=" << &posArray[3] <<
         "gx.=" << &posArray[4] <<
         "gy.=" << &posArray[5] <<
         "rpad.=" << &posArray[6] <<
         "channel.=" << &posArray[7];

      cstream << "calPads" <<
         "\n";

      delete[] posArray;
      delete[] vectorArray;
   } //for (UInt_t isector = 0; isector < tpcROCinstance->GetNSectors(); isector++)

   delete[] names;
   if (mapFileName) {
      delete mapIROCs;
      delete mapOROCs;
      delete[] mapIROCArray;
      delete[] mapOROCArray;
      delete[] mapNames;
   }
}


void AliTPCCalibViewer::MakeTree(const char * fileName, TObjArray * array, const char * mapFileName, AliTPCCalPad *const outlierPad, Float_t ltmFraction) {
  //
  // Write a tree with all available information
  // if mapFileName is speciefied, the Map information are also written to the tree
  // pads specified in outlierPad are not used for calculating statistics
  // The following statistical information on the basis of a ROC are calculated: 
  // "_Median", "_Mean", "_LTM", "_RMS_LTM"
  // "_Median_OutlierCutted", "_Mean_OutlierCutted", "_RMS_OutlierCutted", "_LTM_OutlierCutted", "_RMS_LTM_OutlierCutted"
  // The following position variables are available:
  // "row", "pad", "lx", "ly", "gx", "gy", "rpad", "channel"
  // 
  // The tree out of this function is the basis for the AliTPCCalibViewer and the AliTPCCalibViewerGUI.
   
  AliTPCROC* tpcROCinstance = AliTPCROC::Instance();

  TObjArray* mapIROCs = 0;
  TObjArray* mapOROCs = 0;
  TVectorF *mapIROCArray = 0;
  TVectorF *mapOROCArray = 0;
  Int_t mapEntries = 0;
  TString* mapNames = 0;
  
  if (mapFileName) {
    TFile mapFile(mapFileName, "read");
    
    TList* listOfROCs = mapFile.GetListOfKeys();
    mapEntries = listOfROCs->GetEntries()/2;
    mapIROCs = new TObjArray(mapEntries*2);
    mapOROCs = new TObjArray(mapEntries*2);
    mapIROCArray = new TVectorF[mapEntries];
    mapOROCArray = new TVectorF[mapEntries];
    
    mapNames = new TString[mapEntries];
    for (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
      TString rocName(((TKey*)(listOfROCs->At(ivalue*2)))->GetName());
      rocName.Remove(rocName.Length()-4, 4);
      mapIROCs->AddAt((AliTPCCalROC*)mapFile.Get((rocName + "IROC").Data()), ivalue);
      mapOROCs->AddAt((AliTPCCalROC*)mapFile.Get((rocName + "OROC").Data()), ivalue);
      mapNames[ivalue].Append(rocName);
    }
    
    for (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
      mapIROCArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(0));
      mapOROCArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(36));
      
      for (UInt_t ichannel = 0; ichannel < tpcROCinstance->GetNChannels(0); ichannel++)
        (mapIROCArray[ivalue])[ichannel] = ((AliTPCCalROC*)(mapIROCs->At(ivalue)))->GetValue(ichannel);
      for (UInt_t ichannel = 0; ichannel < tpcROCinstance->GetNChannels(36); ichannel++)
        (mapOROCArray[ivalue])[ichannel] = ((AliTPCCalROC*)(mapOROCs->At(ivalue)))->GetValue(ichannel);
    }
    
  } //  if (mapFileName)
  
  TTreeSRedirector cstream(fileName);
  Int_t arrayEntries = 0;
  if (array) arrayEntries = array->GetEntries();
  
  TString* names = new TString[arrayEntries];
  for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++)
    names[ivalue].Append(((AliTPCCalPad*)array->At(ivalue))->GetName());
  
  for (UInt_t isector = 0; isector < tpcROCinstance->GetNSectors(); isector++) {
      //
      // get statistic for given sector
      //
    TVectorF median(arrayEntries);
    TVectorF mean(arrayEntries);
    TVectorF rms(arrayEntries);
    TVectorF ltm(arrayEntries);
    TVectorF ltmrms(arrayEntries);
    TVectorF medianWithOut(arrayEntries);
    TVectorF meanWithOut(arrayEntries);
    TVectorF rmsWithOut(arrayEntries);
    TVectorF ltmWithOut(arrayEntries);
    TVectorF ltmrmsWithOut(arrayEntries);
    
    TVectorF *vectorArray = new TVectorF[arrayEntries];
    for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++){
      vectorArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));
      vectorArray[ivalue].SetUniqueID(array->UncheckedAt(ivalue)->GetUniqueID());
//       printf("UniqueID: %d\n",vectorArray[ivalue].GetUniqueID());
    }
    
    for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
      AliTPCCalPad* calPad = (AliTPCCalPad*) array->At(ivalue);
      AliTPCCalROC* calROC = calPad->GetCalROC(isector);
      AliTPCCalROC* outlierROC = 0;
      if (outlierPad) outlierROC = outlierPad->GetCalROC(isector);
      if (calROC) {
        median[ivalue] = calROC->GetMedian();
        mean[ivalue] = calROC->GetMean();
        rms[ivalue] = calROC->GetRMS();
        Double_t ltmrmsValue = 0;
        ltm[ivalue] = calROC->GetLTM(&ltmrmsValue, ltmFraction);
        ltmrms[ivalue] = ltmrmsValue;
        if (outlierROC) {
          medianWithOut[ivalue] = calROC->GetMedian(outlierROC);
          meanWithOut[ivalue] = calROC->GetMean(outlierROC);
          rmsWithOut[ivalue] = calROC->GetRMS(outlierROC);
          ltmrmsValue = 0;
          ltmWithOut[ivalue] = calROC->GetLTM(&ltmrmsValue, ltmFraction, outlierROC);
          ltmrmsWithOut[ivalue] = ltmrmsValue;
        }
      }
      else {
        median[ivalue] = 0.;
        mean[ivalue] = 0.;
        rms[ivalue] = 0.;
        ltm[ivalue] = 0.;
        ltmrms[ivalue] = 0.;
        medianWithOut[ivalue] = 0.;
        meanWithOut[ivalue] = 0.;
        rmsWithOut[ivalue] = 0.;
        ltmWithOut[ivalue] = 0.;
        ltmrmsWithOut[ivalue] = 0.;
      }
    }
    
      //
      // fill vectors of variable per pad
      //
    TVectorF *posArray = new TVectorF[8];
    for (Int_t ivalue = 0; ivalue < 8; ivalue++)
      posArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));
    
    Float_t posG[3] = {0};
    Float_t posL[3] = {0};
    Int_t ichannel = 0;
    for (UInt_t irow = 0; irow < tpcROCinstance->GetNRows(isector); irow++) {
      for (UInt_t ipad = 0; ipad < tpcROCinstance->GetNPads(isector, irow); ipad++) {
        tpcROCinstance->GetPositionLocal(isector, irow, ipad, posL);
        tpcROCinstance->GetPositionGlobal(isector, irow, ipad, posG);
        posArray[0][ichannel] = irow;
        posArray[1][ichannel] = ipad;
        posArray[2][ichannel] = posL[0];
        posArray[3][ichannel] = posL[1];
        posArray[4][ichannel] = posG[0];
        posArray[5][ichannel] = posG[1];
        posArray[6][ichannel] = (Int_t)(ipad - (Double_t)(tpcROCinstance->GetNPads(isector, irow))/2);
        posArray[7][ichannel] = ichannel;
        
            // loop over array containing AliTPCCalPads
        for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
          AliTPCCalPad* calPad = (AliTPCCalPad*) array->At(ivalue);
          AliTPCCalROC* calROC = calPad->GetCalROC(isector);
                if (calROC)
                  (vectorArray[ivalue])[ichannel] = calROC->GetValue(irow, ipad);
          else
            (vectorArray[ivalue])[ichannel] = 0;
        }
        ichannel++;
      }
    }
    
    cstream << "calPads" <<
      "sector=" << isector;
    
    for  (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
      cstream << "calPads" <<
        (Char_t*)((names[ivalue] + "_Median=").Data()) << median[ivalue] <<
        (Char_t*)((names[ivalue] + "_Mean=").Data()) << mean[ivalue] <<
        (Char_t*)((names[ivalue] + "_RMS=").Data()) << rms[ivalue] <<
        (Char_t*)((names[ivalue] + "_LTM=").Data()) << ltm[ivalue] <<
        (Char_t*)((names[ivalue] + "_RMS_LTM=").Data()) << ltmrms[ivalue];
      if (outlierPad) {
        cstream << "calPads" <<
          (Char_t*)((names[ivalue] + "_Median_OutlierCutted=").Data()) << medianWithOut[ivalue] <<
          (Char_t*)((names[ivalue] + "_Mean_OutlierCutted=").Data()) << meanWithOut[ivalue] <<
          (Char_t*)((names[ivalue] + "_RMS_OutlierCutted=").Data()) << rmsWithOut[ivalue] <<
          (Char_t*)((names[ivalue] + "_LTM_OutlierCutted=").Data()) << ltmWithOut[ivalue] <<
          (Char_t*)((names[ivalue] + "_RMS_LTM_OutlierCutted=").Data()) << ltmrmsWithOut[ivalue];
      }
        //timestamp and run, if given in title
/*      TString title(((AliTPCCalPad*) array->At(ivalue))->GetTitle());
      TObjArray *arrtitle=title.Tokenize(",");
      Int_t run=-1;
      UInt_t time=0;
      TIter next(arrtitle);
      TObject *o=0;
      while ( (o=next()) ){
        TString &entry=((TObjString*)o)->GetString();
        entry.Remove(TString::kBoth,' ');
        entry.Remove(TString::kBoth,'\t');
        if (entry.BeginsWith("Run:")) {
          run=entry(4,entry.Length()).Atoi();
        } else if (entry.BeginsWith("Time:")) {
          time=entry(6,entry.Length()).Atoi();
        }
      }
      delete arrtitle;*/
      
    }
    
    for  (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
      cstream << "calPads" <<
        (Char_t*)((names[ivalue] + ".=").Data()) << &vectorArray[ivalue];
    }
    
    if (mapFileName) {
          for  (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
            if (isector < 36)
              cstream << "calPads" <<
              (Char_t*)((mapNames[ivalue] + ".=").Data()) << &mapIROCArray[ivalue];
            else
                cstream << "calPads" <<
              (Char_t*)((mapNames[ivalue] + ".=").Data()) << &mapOROCArray[ivalue];
         }
    }

    cstream << "calPads" <<
      "row.=" << &posArray[0] <<
      "pad.=" << &posArray[1] <<
      "lx.=" << &posArray[2] <<
      "ly.=" << &posArray[3] <<
      "gx.=" << &posArray[4] <<
      "gy.=" << &posArray[5] <<
      "rpad.=" << &posArray[6] <<
      "channel.=" << &posArray[7];
    
    cstream << "calPads" <<
      "\n";
    
    delete[] posArray;
    delete[] vectorArray;
  }
  
  
  delete[] names;
  if (mapFileName) {
    delete mapIROCs;
    delete mapOROCs;
    delete[] mapIROCArray;
    delete[] mapOROCArray;
    delete[] mapNames;
  }
}


void AliTPCCalibViewer::MakeTree(const char *outPutFileName, const Char_t *inputFileName, AliTPCCalPad *outlierPad, Float_t ltmFraction, const char *mapFileName ){
   // 
   // Function to create a calibration Tree with all available information.
   // See also documentation to MakeTree   
   // the file "inputFileName" must be in the following format (see also CreateObjectList):
   // (each colum separated by tabs, "dependingOnType" can have 2 or 3 colums)
   // 
   // type	path	dependingOnType
   // 
   // type == "CE":
   // dependingOnType = CETmean	CEQmean	CETrms
   // 
   // type == "Pulser":
   // dependingOnType = PulserTmean	PulsterQmean	PulserTrms
   // 
   // type == "Pedestals":
   // dependingOnType = Pedestals	Noise
   // 
   // type == "CalPad":
   // dependingOnType = NameInFile	NameToWriteToFile
   // 
   // 
   TObjArray objArray;
   CreateObjectList(inputFileName, &objArray);
   MakeTree(outPutFileName, &objArray, mapFileName, outlierPad, ltmFraction);   
}


void AliTPCCalibViewer::CreateObjectList(const Char_t *filename, TObjArray *calibObjects){
   // 
   // Function to create a TObjArray out of a given file
   // the file must be in the following format:
   // (each colum separated by tabs, "dependingOnType" can have 2 or 3 colums)
   // 
   // 
   // type	path	dependingOnType
   // 
   // type == "CE":
   // dependingOnType = CETmean	CEQmean	CETrms
   // 
   // type == "Pulser":
   // dependingOnType = PulserTmean	PulsterQmean	PulserTrms
   // 
   // type == "Pedestals":
   // dependingOnType = Pedestals	Noise
   // 
   // type == "CalPad":
   // dependingOnType = NameInFile	NameToWriteToFile
   // 
   // 
   // 
   if ( calibObjects == 0x0 ) return;
   ifstream in;
   in.open(filename);
   if ( !in.is_open() ){
      fprintf(stderr,"Error: cannot open list file '%s'", filename);
      return;
   }
   
   AliTPCCalPad *calPad=0x0;
   
   TString sFile;
   sFile.ReadFile(in);
   in.close();
   
   TObjArray *arrFileLine = sFile.Tokenize("\n");
   TIter nextLine(arrFileLine);
   
   TObjString *sObjLine = 0x0;
   while ( (sObjLine = (TObjString*)nextLine()) ){
      TString sLine(sObjLine->GetString());
      
      TObjArray *arrCol = sLine.Tokenize("\t");
      Int_t nCols = arrCol->GetEntriesFast();
      
      TObjString *sObjType     = (TObjString*)(arrCol->At(0));
      TObjString *sObjFileName = (TObjString*)(arrCol->At(1));
      TObjString *sObjName = 0x0;
      
      if ( !sObjType || !sObjFileName ) continue;
      TString sType(sObjType->GetString());
      TString sFileName(sObjFileName->GetString());
      printf("Type %s, opening %s \n", sType.Data(), sFileName.Data());
      TFile *fIn = TFile::Open(sFileName);
      if ( !fIn ){
         fprintf(stderr,"File not found: '%s'", sFileName.Data());
         continue;
      }
      
      if ( sType == "CE" ){  // next three colums are the names for CETmean, CEQmean and CETrms
         AliTPCCalibCE *ce = (AliTPCCalibCE*)fIn->Get("AliTPCCalibCE");
         calPad = new AliTPCCalPad((TObjArray*)ce->GetCalPadT0());         
         if (nCols > 2) {  // check, if the name is provided
            sObjName = (TObjString*)(arrCol->At(2));
            calPad->SetNameTitle(sObjName->GetString().Data(), sObjName->GetString().Data());
         }
         else calPad->SetNameTitle("CETmean","CETmean");
         calibObjects->Add(calPad);
         
         calPad = new AliTPCCalPad((TObjArray*)ce->GetCalPadQ());         
         if (nCols > 3) {  // check, if the name is provided
            sObjName = (TObjString*)(arrCol->At(3));
            calPad->SetNameTitle(sObjName->GetString().Data(), sObjName->GetString().Data());
         }
         else calPad->SetNameTitle("CEQmean","CEQmean");
         calibObjects->Add(calPad);        
         
         calPad = new AliTPCCalPad((TObjArray*)ce->GetCalPadRMS());
         if (nCols > 4) {  // check, if the name is provided
            sObjName = (TObjString*)(arrCol->At(4));
            calPad->SetNameTitle(sObjName->GetString().Data(), sObjName->GetString().Data());
         }
         else calPad->SetNameTitle("CETrms","CETrms");
         calibObjects->Add(calPad);         
                  
      } else if ( sType == "Pulser") {
         AliTPCCalibPulser *sig = (AliTPCCalibPulser*)fIn->Get("AliTPCCalibPulser");
         
         calPad = new AliTPCCalPad((TObjArray*)sig->GetCalPadT0());         
         if (nCols > 2) {
            sObjName = (TObjString*)(arrCol->At(2));
            calPad->SetNameTitle(sObjName->GetString().Data(), sObjName->GetString().Data());
         }
         else calPad->SetNameTitle("PulserTmean","PulserTmean");
         calibObjects->Add(calPad);
         
         calPad = new AliTPCCalPad((TObjArray*)sig->GetCalPadQ());         
         if (nCols > 3) {
            sObjName = (TObjString*)(arrCol->At(3));
            calPad->SetNameTitle(sObjName->GetString().Data(), sObjName->GetString().Data());
         }
         else calPad->SetNameTitle("PulserQmean","PulserQmean");
         calibObjects->Add(calPad);        
         
         calPad = new AliTPCCalPad((TObjArray*)sig->GetCalPadRMS());
         if (nCols > 4) {
            sObjName = (TObjString*)(arrCol->At(4));
            calPad->SetNameTitle(sObjName->GetString().Data(), sObjName->GetString().Data());
         }
         else calPad->SetNameTitle("PulserTrms","PulserTrms");
         calibObjects->Add(calPad);         
      
      } else if ( sType == "Pedestals") {
         AliTPCCalibPedestal *ped = (AliTPCCalibPedestal*)fIn->Get("AliTPCCalibPedestal");
         
         calPad = new AliTPCCalPad((TObjArray*)ped->GetCalPadPedestal());         
         if (nCols > 2) {
            sObjName = (TObjString*)(arrCol->At(2));
            calPad->SetNameTitle(sObjName->GetString().Data(), sObjName->GetString().Data());
         }
         else calPad->SetNameTitle("Pedestals","Pedestals");
         calibObjects->Add(calPad);
         
         calPad = new AliTPCCalPad((TObjArray*)ped->GetCalPadRMS());         
         if (nCols > 3) {
            sObjName = (TObjString*)(arrCol->At(3));
            calPad->SetNameTitle(sObjName->GetString().Data(), sObjName->GetString().Data());
         }
         else calPad->SetNameTitle("Noise","Noise");
         calibObjects->Add(calPad);        
     
      } else if ( sType == "CalPad") {
         if (nCols > 2) sObjName = (TObjString*)(arrCol->At(2));
         else continue;
         calPad = new AliTPCCalPad(*(AliTPCCalPad*)fIn->Get(sObjName->GetString().Data()));
         if (nCols > 3) {
            sObjName = (TObjString*)(arrCol->At(3));
            calPad->SetNameTitle(sObjName->GetString().Data(), sObjName->GetString().Data());
         }
         calibObjects->Add(calPad);
      } else {
         fprintf(stderr,"Undefined Type: '%s'",sType.Data());
      }
      delete fIn;
   }
}



