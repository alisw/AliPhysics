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
///////////////////////////////////////////////////////////////////////////////

//
// ROOT includes 
//
#include <iostream>
#include <TString.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1F.h>
#include <THashTable.h>
#include <TObjString.h>
#include "TTreeStream.h"
#include "TFile.h"
#include "TKey.h"


//
// AliRoot includes
//
#include "AliTPCCalibViewer.h"

ClassImp(AliTPCCalibViewer)

AliTPCCalibViewer::AliTPCCalibViewer()
                  :TObject(),
                   fTree(0),
                   fFile(0),
                   fListOfObjectsToBeDeleted(0)
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
                   fListOfObjectsToBeDeleted(0)
{
  //
  // dummy AliTPCCalibViewer copy constructor
  // not yet working!!!
  //
  fTree = c.fTree;
  //fFile = new TFile(*(c.fFile));
  fListOfObjectsToBeDeleted = c.fListOfObjectsToBeDeleted;
}

//_____________________________________________________________________________
AliTPCCalibViewer::AliTPCCalibViewer(TTree* tree)
                  :TObject(),
                   fTree(0),
                   fFile(0),
                   fListOfObjectsToBeDeleted(0)
{
  //
  // Constructor that initializes the calibration viewer
  //
  fTree = tree;
  fListOfObjectsToBeDeleted = new TObjArray();
}

//_____________________________________________________________________________
AliTPCCalibViewer::AliTPCCalibViewer(char* fileName, char* treeName)
                  :TObject(),
                   fTree(0),
                   fFile(0),
                   fListOfObjectsToBeDeleted(0)
{
   //
   // Constructor to initialize the calibration viewer
   // the file 'fileName' contains the tree 'treeName'
   //
   fFile = new TFile(fileName, "read");
   fTree = (TTree*) fFile->Get(treeName);
   fListOfObjectsToBeDeleted = new TObjArray();
}
                   
//____________________________________________________________________________
AliTPCCalibViewer & AliTPCCalibViewer::operator =(const AliTPCCalibViewer & param)
{
   //
   // assignment operator - dummy
   // not yet working!!!
   //
   fTree = param.fTree;
   //fFile = new TFile(*(param.fFile));
   fListOfObjectsToBeDeleted = param.fListOfObjectsToBeDeleted;
   return (*this);
}

//_____________________________________________________________________________
AliTPCCalibViewer::~AliTPCCalibViewer()
{
   //
   // AliTPCCalibViewer destructor
   // all objects will be deleted, the file will be closed, the pictures will disapear
   //
   /*if (fTree) {
      delete fTree;
      fTree = 0;
   }*/
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
   TString drawOptionsStr("profcolz ");
   TRandom rnd(0);
   Int_t rndNumber = rnd.Integer(10000);
   if (drawOptions && drawOptions != "")
      drawOptionsStr += drawOptions;

   if (sectorStr == "A") {
      drawStr += ":gy.fElements:gx.fElements>>prof";
      drawStr += rndNumber;
      drawStr += "(330,-250,250,330,-250,250)";
      cutStr += "(sector/18)%2==0 ";
   }
   else if  (sectorStr == "C") {
      drawStr += ":gy.fElements:gx.fElements>>prof";
      drawStr += rndNumber;
      drawStr += "(330,-250,250,330,-250,250)";
      cutStr += "(sector/18)%2==1 ";
   }
   else if  (sectorStr == "ALL") {
      drawStr += ":gy.fElements:gx.fElements>>prof";
      drawStr += rndNumber;
      drawStr += "(330,-250,250,330,-250,250)";
   }
   else if (sectorStr.IsDigit()) {
      Int_t isec = sectorStr.Atoi();
      drawStr += ":rpad.fElements:row.fElements>>prof";
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
   drawStr.ReplaceAll("~", ".fElements");
   cutStr.ReplaceAll("~", ".fElements");
   if (writeDrawCommand) std::cout << "fTree->Draw(\"" << drawStr << "\", \"" <<  cutStr << "\", \"" << drawOptionsStr << "\");" << std::endl;
   return fTree->Draw(drawStr.Data(), cutStr.Data(), drawOptionsStr.Data());
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
      char sectorChr[3];
      sprintf(sectorChr, "%i", sector);
      return EasyDraw(drawCommand, sectorChr, cuts, drawOptions, writeDrawCommand);
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

   if (cuts && cuts[0] != 0) {
      if (cutStr.Length() != 0) cutStr += "&& ";
      cutStr += "(";
      cutStr += cuts;
      cutStr += ")";
   }

   drawStr.ReplaceAll("~", ".fElements");
   cutStr.ReplaceAll("~", ".fElements");
   if (writeDrawCommand) std::cout << "fTree->Draw(\"" << drawStr << "\", \"" <<  cutStr << "\", \"" << drawOptionsStr << "\");" << std::endl;
   return fTree->Draw(drawStr.Data(), cutStr.Data(), drawOptionsStr.Data());
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
      char sectorChr[3];
      sprintf(sectorChr, "%i", sector);
      return EasyDraw1D(drawCommand, sectorChr, cuts, drawOptions, writeDrawCommand);
   }
  Error("EasyDraw","The TPC contains only sectors between 0 and 71.");
  return -1;
}

//_____________________________________________________________________________
Int_t AliTPCCalibViewer::DrawHisto1D(const char* type, Int_t sector, TVectorF& nsigma, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM) const {
  //
  // draws a 1-dimensional histogram of 'type' for sector 'sector'
  // TVectorF nsigma: Specifies, for which distances from the mean/median/LTM lines should be drawn, in units of sigma
  // example: nsigma={2, 4, 6}: Three lines will be drawn, distance to mean/median/LTM: 2, 3 and 6 sigma
  // plotMean, plotMedian, plotLTM: specifies, if mean, median and LTM should be drawn as lines into the histogram
  //

   TString typeStr(type);
   TString sectorStr("sector==");
   sectorStr += sector;

   TCanvas* canvas = ((TCanvas*)gROOT->GetListOfCanvases()->Last());
   Int_t oldOptStat = gStyle->GetOptStat();
   gStyle->SetOptStat(0000000);

   if (!canvas) {
      canvas = new TCanvas();
      fListOfObjectsToBeDeleted->Add(canvas);
   }
   
   char c[500];
   sprintf(c, "%s, sector: %i", type, sector);
   TLegend * legend = new TLegend(.8,.6, .99, .99, c);
   fListOfObjectsToBeDeleted->Add(legend);

   Int_t nentries = fTree->Draw((typeStr+".fElements").Data(), sectorStr.Data(), "");
   ((TH1F*)canvas->GetPrimitive("htemp"))->SetTitle("");
      
   //****************************************************************
   //!!!!!!!!!!!!!!!!! Needs further investigaton !!!!!!!!!!!!!!!!!!!
   //****************************************************************
   //fListOfObjectsToBeDeleted->Add(canvas->GetPrimitive("htemp"));
/*
  By default the temporary histogram created is called "htemp", but only in
  the one dimensional Draw("e1") it contains the TTree's data points. For
  a two dimensional Draw, the data is filled into a TGraph which is named
  "Graph". They can be retrieved by calling
    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp"); // 1D
    TGraph *graph = (TGraph*)gPad->GetPrimitive("Graph"); // 2D
*/
   
   canvas->Update();
   Double_t sigma = 0;

   if (plotMean) {
      fTree->Draw((typeStr+"_Mean").Data(), sectorStr.Data(), "goff");
      Double_t lineX = fTree->GetV1()[0];
      fTree->Draw((typeStr+"_RMS").Data(), sectorStr.Data(), "goff");
      sigma = fTree->GetV1()[0];
      TLine* line = new TLine(lineX, 0, lineX, canvas->GetUymax());
      fListOfObjectsToBeDeleted->Add(line);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      line->SetLineStyle(1);
      line->Draw();
      sprintf(c, "Mean: %f", lineX);
      legend->AddEntry(line, c, "l");

      for (Int_t i = 0; i < nsigma.GetNoElements(); i++) {
         TLine* linePlusSigma = new TLine(lineX+nsigma[i]*sigma, 0, lineX+nsigma[i]*sigma, canvas->GetUymax());
         fListOfObjectsToBeDeleted->Add(linePlusSigma);
         linePlusSigma->SetLineColor(kRed);
         linePlusSigma->SetLineStyle(2+i);
         linePlusSigma->Draw();
   
         TLine* lineMinusSigma = new TLine(lineX-nsigma[i]*sigma, 0, lineX-nsigma[i]*sigma, canvas->GetUymax());
         fListOfObjectsToBeDeleted->Add(lineMinusSigma);
         lineMinusSigma->SetLineColor(kRed);
         lineMinusSigma->SetLineStyle(2+i);
         lineMinusSigma->Draw();
         sprintf(c, "%i #sigma = %f",(Int_t)(nsigma[i]), (Float_t)(nsigma[i]*sigma));
	 std::cout << "nsigma-char*: " << c << std::endl;
         legend->AddEntry(lineMinusSigma, c, "l");
      }
   }

   if (plotMedian) {
      fTree->Draw((typeStr+"_Median").Data(), sectorStr.Data(), "goff");
      Double_t lineX = fTree->GetV1()[0];
      fTree->Draw((typeStr+"_RMS").Data(), sectorStr.Data(), "goff");
      sigma = fTree->GetV1()[0];
      TLine* line = new TLine(lineX, 0, lineX, canvas->GetUymax());
      fListOfObjectsToBeDeleted->Add(line);
      line->SetLineColor(kBlue);
      line->SetLineWidth(2);
      line->SetLineStyle(1);
      line->Draw();
      sprintf(c, "Median: %f", lineX);
      legend->AddEntry(line, c, "l");
      
      for (Int_t i = 0; i < nsigma.GetNoElements(); i++) {
         TLine* linePlusSigma = new TLine(lineX+nsigma[i]*sigma, 0, lineX+nsigma[i]*sigma, canvas->GetUymax());
         fListOfObjectsToBeDeleted->Add(linePlusSigma);
         linePlusSigma->SetLineColor(kBlue);
         linePlusSigma->SetLineStyle(2+i);
         linePlusSigma->Draw();
   
         TLine* lineMinusSigma = new TLine(lineX-nsigma[i]*sigma, 0, lineX-nsigma[i]*sigma, canvas->GetUymax());
         fListOfObjectsToBeDeleted->Add(lineMinusSigma);
         lineMinusSigma->SetLineColor(kBlue);
         lineMinusSigma->SetLineStyle(2+i);
         lineMinusSigma->Draw();
         sprintf(c, "%i #sigma = %f",(Int_t)(nsigma[i]), (Float_t)(nsigma[i]*sigma));
         legend->AddEntry(lineMinusSigma, c, "l");
      }
   }
   
   if (plotLTM) {
      fTree->Draw((typeStr+"_LTM").Data(), sectorStr.Data(), "goff");
      Double_t lineX = fTree->GetV1()[0];
      fTree->Draw((typeStr+"_RMS_LTM").Data(), sectorStr.Data(), "goff");
      sigma = fTree->GetV1()[0];
      TLine* line = new TLine(lineX, 0, lineX, canvas->GetUymax());
      fListOfObjectsToBeDeleted->Add(line);
      line->SetLineColor(kGreen+2);
      line->SetLineWidth(2);
      line->SetLineStyle(1);
      line->Draw();
      sprintf(c, "LTM: %f", lineX);
      legend->AddEntry(line, c, "l");

      for (Int_t i = 0; i < nsigma.GetNoElements(); i++) {
         TLine* linePlusSigma = new TLine(lineX+nsigma[i]*sigma, 0, lineX+nsigma[i]*sigma, canvas->GetUymax());
         fListOfObjectsToBeDeleted->Add(linePlusSigma);
         linePlusSigma->SetLineColor(kGreen+2);
         linePlusSigma->SetLineStyle(2+i);
         linePlusSigma->Draw();
   
         TLine* lineMinusSigma = new TLine(lineX-nsigma[i]*sigma, 0, lineX-nsigma[i]*sigma, canvas->GetUymax());
         fListOfObjectsToBeDeleted->Add(lineMinusSigma);
         lineMinusSigma->SetLineColor(kGreen+2);
         lineMinusSigma->SetLineStyle(2+i);
         lineMinusSigma->Draw();
         sprintf(c, "%i #sigma = %f", (Int_t)(nsigma[i]), (Float_t)(nsigma[i]*sigma));
         legend->AddEntry(lineMinusSigma, c, "l");
      }
   }
   legend->Draw();
   gStyle->SetOptStat(oldOptStat);
   return nentries;
}

//_____________________________________________________________________________
void AliTPCCalibViewer::SigmaCut(const char* type, Int_t sector, Float_t sigmaMax, Float_t sigmaStep, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM) const {
  //
  // Creates a histogram, where you can see, how much of the data are inside sigma-intervals around the mean/median/LTM
  // type: For which type of data the histogram is generated, e.g. 'CEQmean'
  // sector: For which sector the histogram is generated
  // sigmaMax: up to which sigma around the mean/median/LTM the histogram is generated
  // sigmaStep: the binsize of the generated histogram
  // plotMean/plotMedian/plotLTM: specifies where to put the center
  //
   Int_t oldOptStat = gStyle->GetOptStat();
   gStyle->SetOptStat(0000000);

   TString typeStr(type);
   TString sectorStr("sector==");
   sectorStr += sector;
   Int_t entries = fTree->Draw((typeStr+".fElements").Data(), sectorStr.Data(), "goff");
   char headline[500];
   sprintf(headline, "%s in sector %i; Multiples of #sigma; Fraction of used pads", type, sector);
   TH1F *histMean =   new TH1F("histMean",headline, (Int_t)(sigmaMax/sigmaStep+1), 0, sigmaMax+sigmaStep);
   sprintf(headline, "%s in sector %i; Multiples of #sigma; Fraction of used pads", type, sector);
   TH1F *histMedian = new TH1F("histMedian",headline, (Int_t)(sigmaMax/sigmaStep+1), 0, sigmaMax+sigmaStep);
   sprintf(headline, "%s in sector %i; Multiples of #sigma; Fraction of used pads", type, sector);
   TH1F *histLTM =    new TH1F("histLTM",headline, (Int_t)(sigmaMax/sigmaStep+1), 0, sigmaMax+sigmaStep);                                                                   
   histMean->SetDirectory(0);
   histMedian->SetDirectory(0);
   histLTM->SetDirectory(0);
   fListOfObjectsToBeDeleted->Add(histMean);
   fListOfObjectsToBeDeleted->Add(histMedian);
   fListOfObjectsToBeDeleted->Add(histLTM);

   
   // example-cut: sector==34 && TMath::Abs(CEQmean.fElements - CEQmean_Mean) < nsigma * CEQmean_RMS
   for (Float_t nsigma = 0; nsigma <= sigmaMax; nsigma += sigmaStep) {
     std::cout << "Calculating histograms,  step: " << (Int_t)(nsigma/sigmaStep) << " of: " << (Int_t)(sigmaMax/sigmaStep) << "\r" << std::flush;
      char cuts[5000];
      
      if (plotMean) {
         sprintf(cuts, "sector==%i && ( %s.fElements - %s_Median) < %f * %s_RMS", sector, type, type, nsigma, type );
         sprintf(cuts,         "%s && (-%s.fElements + %s_Median) < %f * %s_RMS",   cuts, type, type, nsigma, type );
         Float_t value = (Float_t)fTree->Draw((typeStr+".fElements").Data(), cuts, "goff")/entries;
         histMean->Fill(nsigma, value);
      }
      if (plotMedian) {
         sprintf(cuts, "sector==%i && ( %s.fElements - %s_Mean) < %f * %s_RMS", sector, type, type, nsigma, type );
         sprintf(cuts,         "%s && (-%s.fElements + %s_Mean) < %f * %s_RMS",   cuts, type, type, nsigma, type );
         Float_t value = (Float_t)fTree->Draw((typeStr+".fElements").Data(), cuts, "goff")/entries;
         histMedian->Fill(nsigma, value);
      }
      if (plotLTM) {
         sprintf(cuts, "sector==%i && ( %s.fElements - %s_LTM) < %f * %s_RMS_LTM", sector, type, type, nsigma, type );
         sprintf(cuts,         "%s && (-%s.fElements + %s_LTM) < %f * %s_RMS_LTM",   cuts, type, type, nsigma, type );
         Float_t value = (Float_t)fTree->Draw((typeStr+".fElements").Data(), cuts, "goff")/entries;
         histLTM->Fill(nsigma, value);
      }
   }
   
   char c[500];
   sprintf(c, "Sigma Cut");
   TLegend * legend = new TLegend(.85,.8, .99, .99, c);
   fListOfObjectsToBeDeleted->Add(legend);
   
   if (plotMean){
      histMean->SetLineColor(kBlack);
      sprintf(c, "Mean");
      legend->AddEntry(histMean, c, "l");
      histMean->Draw();
   }
   if (plotMedian){
      histMedian->SetLineColor(kRed);
      sprintf(c, "Median");
      legend->AddEntry(histMedian, c, "l");
      histMedian->Draw("same");
   }
   if (plotLTM){
      histLTM->SetLineColor(kBlue);
      sprintf(c, "LTM");
      legend->AddEntry(histLTM, c, "l");
      histLTM->Draw("same");
   }   

   legend->Draw();
   gStyle->SetOptStat(oldOptStat);
}


//_____________________________________________________________________________
AliTPCCalPad* AliTPCCalibViewer::GetCalPad(const char* desiredData, char* cuts, char* calPadName) const {
  //
  // creates a AliTPCCalPad out of the 'desiredData'
  // the functionality of EasyDraw1D is used
  // calPadName specifies the name of the created AliTPCCalPad
  //  - this takes a while -
  //
   TString drawStr(desiredData);
   drawStr.Append(":channel~");
   AliTPCCalPad * createdCalPad = new AliTPCCalPad(calPadName, calPadName);
   Int_t entries = 0;
   for (Int_t sec = 0; sec < 72; sec++) {
      entries = EasyDraw1D(drawStr.Data(), (Int_t)sec, cuts, "goff");
      for (Int_t i = 0; i < entries; i++) 
         createdCalPad->GetCalROC(sec)->SetValue((UInt_t)(fTree->GetV2()[i]), (Float_t)(fTree->GetV1()[i]));
   }
   return createdCalPad;   
}

//_____________________________________________________________________________
AliTPCCalROC* AliTPCCalibViewer::GetCalROC(const char* desiredData, UInt_t sector, char* cuts) const {
  //
  // creates a AliTPCCalROC out of the desiredData
  // the functionality of EasyDraw1D is used
  // sector specifies the sector of the created AliTPCCalROC
  //
   TString drawStr(desiredData);
   drawStr.Append(":channel~");
   Int_t entries = EasyDraw1D(drawStr.Data(), (Int_t)sector, cuts, "goff");
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


TObjArray* AliTPCCalibViewer::GetListOfNormalizationVariables(Bool_t printList) {
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
   arr->Add(new TObjString("LFitIntern_4_8.fElements"));
   arr->Add(new TObjString("GFitIntern_Lin.fElements"));
   arr->Add(new TObjString("GFitIntern_Par.fElements"));

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
      drawCommand.Append("~");
      AliTPCCalPad* calPad = GetCalPad(drawCommand.Data(), "", calPadName); 
      calPadsArray->Add(calPad); 
   }
   std::cout << std::endl;
   listOfCalPads->Delete();
   delete listOfCalPads;
   return calPadsArray;
}


void AliTPCCalibViewer::MakeTreeWithObjects(const char * fileName, TObjArray * array, const char * mapFileName) {
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
         TString ROCname(((TKey*)(listOfROCs->At(ivalue*2)))->GetName());
         ROCname.Remove(ROCname.Length()-4, 4);
         mapIROCs->AddAt((AliTPCCalROC*)mapFile.Get((ROCname + "IROC").Data()), ivalue);
         mapOROCs->AddAt((AliTPCCalROC*)mapFile.Get((ROCname + "OROC").Data()), ivalue);
         mapNames[ivalue].Append(ROCname);
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

void AliTPCCalibViewer::MakeTree(const char * fileName, TObjArray * array, const char * mapFileName, AliTPCCalPad* outlierPad, Float_t ltmFraction) {
  //
  // Write a tree with all available information
  // im mapFileName is speciefied, the Map information are also written to the tree
  // pads specified in outlierPad are not used for calculating statistics
  //  - the same function as AliTPCCalPad::MakeTree - 
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
         TString ROCname(((TKey*)(listOfROCs->At(ivalue*2)))->GetName());
         ROCname.Remove(ROCname.Length()-4, 4);
         mapIROCs->AddAt((AliTPCCalROC*)mapFile.Get((ROCname + "IROC").Data()), ivalue);
         mapOROCs->AddAt((AliTPCCalROC*)mapFile.Get((ROCname + "OROC").Data()), ivalue);
         mapNames[ivalue].Append(ROCname);
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
      for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++)
         vectorArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));
      
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

