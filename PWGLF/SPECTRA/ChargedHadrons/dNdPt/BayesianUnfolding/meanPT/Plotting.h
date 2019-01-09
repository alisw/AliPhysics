#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TExec.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TKey.h"
#include "TROOT.h"
#include "TPaletteAxis.h"
#include <sstream>


void PlotArray(TObjArray *arraytoplot, const char *controlstring ,Short_t *colorArray, Short_t *markerArray, Double_t relativeTextSize,Float_t lableSize, Float_t titleSize, Int_t lableFont, Int_t titleFont, Float_t titleOffsetY, Float_t titleOffsetX, Float_t margin,Double_t xMin,Double_t xMax, Float_t titleOffsetZ, Short_t* colorArrayRatios = 0){
  TString control(controlstring);
  Int_t nbrEntries = arraytoplot->GetEntries();

  for(Int_t hh = 0; hh < nbrEntries; hh++){
    if (arraytoplot->At(hh)->InheritsFrom("TH2") && (hh==0)){
      TH2D *hist2d = (TH2D*) arraytoplot->At(hh);
      hist2d->SetStats(0);
      hist2d->GetXaxis()->SetLabelSize(lableSize*relativeTextSize);
      hist2d->GetXaxis()->SetTitleSize(titleSize*relativeTextSize);
      hist2d->GetXaxis()->SetLabelFont(lableFont);
      hist2d->GetXaxis()->SetTitleFont(titleFont);
      hist2d->GetXaxis()->SetTitleOffset(titleOffsetX);
      hist2d->GetYaxis()->SetLabelSize(lableSize*relativeTextSize);
      hist2d->GetYaxis()->SetTitleSize(titleSize*relativeTextSize);
      hist2d->GetYaxis()->SetLabelFont(lableFont);
      hist2d->GetYaxis()->SetTitleFont(titleFont);
      hist2d->GetYaxis()->SetTitleOffset(titleOffsetY);
      hist2d->GetZaxis()->SetLabelSize(lableSize*relativeTextSize);
      hist2d->GetZaxis()->SetTitleSize(titleSize*relativeTextSize);
      hist2d->GetZaxis()->SetLabelFont(lableFont);
      hist2d->GetZaxis()->SetTitleFont(titleFont);
      hist2d->GetZaxis()->SetTitleOffset(titleOffsetZ);
      hist2d->SetTitle("");
      if(control.Contains("Thick")||control.Contains("thick")){
	       hist2d->SetLineWidth(2);
      }
      hist2d->Draw("colz");

    }else if (arraytoplot->At(hh)->InheritsFrom("TH1")){
      TH1D *hist = (TH1D*) arraytoplot->At(hh);
      hist->SetStats(0);
      hist->GetXaxis()->SetLabelSize(lableSize*relativeTextSize);
      hist->GetXaxis()->SetTitleSize(titleSize*relativeTextSize);
      hist->GetXaxis()->SetLabelFont(lableFont);
      hist->GetXaxis()->SetTitleFont(titleFont);
      hist->GetXaxis()->SetTitleOffset(titleOffsetX);
      hist->GetYaxis()->SetLabelSize(lableSize*relativeTextSize);
      hist->GetYaxis()->SetTitleSize(titleSize*relativeTextSize);
      hist->GetYaxis()->SetLabelFont(lableFont);
      hist->GetYaxis()->SetTitleFont(titleFont);
      hist->GetYaxis()->SetTitleOffset(titleOffsetY);
      hist->SetTitle("");
      hist->SetAxisRange(xMin,xMax,"X");
      hist->SetMarkerStyle(markerArray[hh]);
      hist->SetMarkerColor(colorArray[hh]);
      hist->SetLineColor(colorArray[hh]);
      if(control.Contains("Thick")||control.Contains("thick")){
        hist->SetLineWidth(2);
        hist->SetMarkerSize(1.2);

      }
      if(control.Contains("fill")||control.Contains("Fill")){
        hist->SetFillStyle(1001);
        hist->SetFillColor(colorArray[hh]);
      }
      if(control.Contains("fuell")){
      	if(hh < 4) {
      	  hist->SetFillStyle(1001);
      	  hist->SetFillColor(colorArray[hh]);
        }
      }
      if(control.Contains("Lines")||control.Contains("lines")){
      	if(control.Contains("2xthick")){
      	  hist->SetLineWidth(4);
      	}
      	if(control.Contains("4xthick")){
      	  hist->SetLineWidth(8);
      	}
      	if(markerArray){ hist->SetLineStyle(markerArray[hh]);}
      	if(hh == 0) { hist->Draw("hist");}
      	else{hist->Draw("Same hist"); }
      }else{
        if(control.Contains("dpg") && hh == 1){
          hist->SetFillStyle(1001);
          hist->SetFillColor(colorArray[hh]);
          hist->Draw("Same hist");
          continue;
        }
      	TString histName(hist->GetName());
      	if(histName.Contains("SystErr")){
      	  TExec *errorBoxesOn = new TExec("errorBoxesOn","gStyle->SetErrorX(0.48)");
      	  errorBoxesOn->Draw();
      	  hist->SetFillStyle(0);
          hist->SetMarkerStyle(markerArray[hh-1]);
          hist->SetMarkerColor(colorArray[hh-1]);
          hist->SetLineColor(colorArray[hh-1]);
          hist->Draw("SameE2");
      	  TExec *errorBoxesOff = new TExec("errorBoxesOff","gStyle->SetErrorX(0)");
      	  errorBoxesOff->Draw();
      	}else{
      	  if(hh==0) {hist->Draw("");}
      	  else{hist->Draw("Same");}
      	}
      }
    }else if (arraytoplot->At(hh)->InheritsFrom("TLine")){
      TLine *line = (TLine*) arraytoplot->At(hh);
      line->Draw("Same");
    }else if (arraytoplot->At(hh)->InheritsFrom("TPaveText")){
      TPaveText *pt = (TPaveText*)arraytoplot->At(hh);
      if (!pt) cout << "Error with TPaveText" << endl;
      pt->SetFillStyle(0);
      pt->SetFillColor(0);
      pt->SetBorderSize(0);
      pt->SetTextFont(titleFont);
      pt->SetTextSize(relativeTextSize);
      pt->Draw();
    }else if (arraytoplot->At(hh)->InheritsFrom("TLegend")){
      TLegend *leghh = (TLegend*) arraytoplot->At(hh);
      leghh->SetTextFont(titleFont);
      leghh->SetTextSize(relativeTextSize);
      leghh->SetFillColor(kWhite);
      leghh->SetFillStyle(0);
      leghh->SetBorderSize(0);
      leghh->Draw("same");
    }else if (arraytoplot->At(hh)->InheritsFrom("TBox")){
      TBox *box = (TBox*) arraytoplot->At(hh);
      box->Draw();
    }else if (arraytoplot->At(hh)->InheritsFrom("TGraphAsymmErrors")){
      TGraphAsymmErrors *asym = (TGraphAsymmErrors*) arraytoplot->At(hh);
      asym->SetLineColor(colorArray[hh]);
      asym->SetLineStyle(markerArray[hh]);
      asym->SetLineWidth(4);
      asym->Draw("SAME");
    }else if (arraytoplot->At(hh)->InheritsFrom("TGraphErrors")){
      TGraphErrors *sym = (TGraphErrors*) arraytoplot->At(hh);
      sym->SetLineColor(colorArray[hh]);
      sym->SetLineStyle(markerArray[hh]);
      sym->SetLineWidth(4);
      sym->Draw("SAME P");
    }else if (arraytoplot->At(hh)->InheritsFrom("TF1")){
      TF1 *fun = (TF1*) arraytoplot->At(hh);
      fun->SetLineColor(colorArray[hh]);
      fun->SetLineStyle(markerArray[hh]);
      if(control.Contains("Thick")||control.Contains("thick")) fun->SetLineWidth(2);
      if(control.Contains("2xthick"))fun->SetLineWidth(4);
      fun->Draw("SAME");
    }
    gPad->RedrawAxis();
  }
}


TCanvas *makeCanvas(string canvasName, TObjArray *histArray, TObjArray *ratioArray = 0,const char *controlstring="", Short_t *colorArray = 0, Short_t *markerArray = 0, Short_t *colorArrayRatios = 0, Short_t *markerArrayRatios = 0){

  gStyle->SetTextFont(43);
  TString control(controlstring);
  if (control.Contains("CMYK")||control.Contains("cmyk")){gStyle->SetColorModelPS(1);}

  Short_t defaultColorArray[14]={4,2,kGreen+3,kMagenta+2,8,9,11,12,13,14,15,16,17,18};
  Short_t defaultMarkerArray[14]={20,21,34,33,24,25,26,27,28,29,30,2,3,5};

  Int_t nHist = histArray->GetEntries();
  if(nHist == 0){cout << "ERROR: histArrays is empty!" << endl; return NULL;}

  // ----------------------- Settings ---------------------------------------------------------------
  Float_t lableSize       = 1.00;   // (1.0)Factor the label will be larger than the relative textsize
  Float_t titleSize       = 1.20;   // Factor the title will be larger than the relative textsize
  Float_t textSizeFactor  = 17000;  //(15000) 12000
  Int_t lableFont         = 43;
  Int_t titleFont         = 43;

  Float_t titleOffsetY    = 1.5;
  Float_t titleOffsetX    = 1.5;    //3.5
  Float_t titleOffsetZ    = 2.5;
  Float_t margin          = 0.12;

  if(control.Contains("square")||control.Contains("Square")||control.Contains("SQUARE")){
    titleOffsetY  = 1.2; //1.4;
    titleOffsetX  = 1.1; //1.5;
    titleOffsetZ  = 1.4;
    margin        = 0.12;
    if(control.Contains("2D")){titleOffsetY= 1.0; titleOffsetX= 1.1;}
  }else if (control.Contains("horizontal")||control.Contains("Horizontal")||control.Contains("HORIZONTAL")){
    titleOffsetY  = 1.05;
    titleOffsetX  = 1.02;
    margin        = 0.09;
  }else if (control.Contains("A4")||control.Contains("a4")){
    titleOffsetY  = 1.4;
    titleOffsetX  = 1.4;
    margin        = 0.11;
  }
  if(ratioArray){
    titleOffsetY  = 1.3;  //1.5;
    titleOffsetX  = 3.1;  //4.1;
  }
  // ---------------------------------------------------------------------------------------------------


  TPad *upperPad = NULL;
  TRandom *random = new TRandom(histArray->Hash());
  TString title("default");

  Double_t xMin = 0;
  Double_t xMax = 0;
  /// Check if first element of the histogram array is indeed an Histogram.
  if(histArray->At(0)){
    if (!histArray->At(0)->InheritsFrom("TH1")){cout << "ERROR: First entry has to be an Histogram" << endl; return NULL;}
    TH1D *hist0 = (TH1D*) histArray->At(0);
    if(!hist0){cout<<"| ERROR: Histogram empty"<<endl;}
    title += Form("%f",random->Rndm());
    xMin = hist0->GetXaxis()->GetBinLowEdge(hist0->GetXaxis()->GetFirst())+0.0000001;
    xMax = hist0->GetXaxis()->GetBinUpEdge(hist0->GetXaxis()->GetLast())-0.0000001;
    hist0->GetYaxis()->CenterTitle(0);
  }
  /// Check if the first Element of the ratio array is an histogram.
  if(ratioArray){
    if(ratioArray->At(0)){
      if (!ratioArray->At(0)->InheritsFrom("TH1")){cout << "ERROR: First entry has to be an Histogram" << endl; return NULL;}
      TH1D *ratio0 = (TH1D*) ratioArray->At(0);
      ratio0->GetXaxis()->SetTickLength(0.06);
      ratio0->GetYaxis()->SetNdivisions(506);
      ratio0->SetLineWidth(2);
      ratio0->GetYaxis()->CenterTitle(1);
    }
  }


  /// Create Canvas in a given size
  TCanvas *newCanvas = NULL;
  if(control.Contains("square")||control.Contains("Square")||control.Contains("SQUARE")){
    newCanvas = new TCanvas(title,title,10,10,710,710);
  }else if (control.Contains("horizontal")||control.Contains("Horizontal")||control.Contains("HORIZONTAL")){
    newCanvas = new TCanvas(title,title,0,0,600*1.414213562 ,600);
  }else if (control.Contains("A4")||control.Contains("a4")){
    newCanvas = new TCanvas(title,title,0,0,600,600*1.414213562);
  }else{
    newCanvas = new TCanvas(title,title,10,10,710,810);
  }
  newCanvas->SetFillStyle(4000);
  newCanvas->cd();


  /// If ratioArray exist, split the Canvas in two and create a ratio panel
  if(ratioArray){
    /// Create upper pad for the histograms \remark upper pad has 72%
    upperPad = new TPad("upperPad","Distribution",0,0.28,1,1);
    TPad *lowerPad = new TPad("lowerPad","Ratio",0,0,1,0.28);

    upperPad->SetFillStyle(4000); /// \note  \c SetFillStyle(4000) gives transparent canvas
    upperPad->SetTopMargin(0.03);
    upperPad->SetLeftMargin(0.12);
    upperPad->SetRightMargin(0.05);
    upperPad->SetBottomMargin(0.0);
    upperPad->Draw();

    lowerPad->SetFillStyle(4000);
    lowerPad->SetTopMargin(0.0);
    lowerPad->SetLeftMargin(0.12);
    lowerPad->SetRightMargin(0.05);
    lowerPad->SetBottomMargin(0.35); //0.3
    lowerPad->Draw();


    if(control.Contains("logX")||control.Contains("logx")||control.Contains("LogX")||control.Contains("LOGX")){upperPad->SetLogx(1); lowerPad->SetLogx(1);}
    if(control.Contains("logY")||control.Contains("logy")||control.Contains("LogY")||control.Contains("LOGY")){upperPad->SetLogy(1);}
    if(control.Contains("RatioGridY")){lowerPad->SetGridy(1);}

    newCanvas->cd();

    Double_t relativeTextSize;
    Double_t pad_width;
    Double_t pad_height;

    /// Calculate the relative text size, independent of the size of the pad
    pad_width  = gPad->XtoPixel(gPad->GetX2());
    pad_height = gPad->YtoPixel(gPad->GetY1());
    if (pad_width < pad_height){relativeTextSize = 1/pad_width;}
    else{relativeTextSize = 1/pad_height;}
    relativeTextSize = textSizeFactor * relativeTextSize;

    upperPad->cd();

    /// Use the standardtised PlotArray function to loop over the histogram an to draw the elements.
    if(!markerArray && !colorArray){PlotArray(histArray,controlstring,defaultColorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
    else if(!markerArray){PlotArray(histArray,controlstring,colorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
    else if(!colorArray){ PlotArray(histArray,controlstring,defaultColorArray,markerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
    else{PlotArray(histArray,controlstring,colorArray,markerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}

    /// Go to the ratio pad and repeat the procedure
    lowerPad->cd();
    if(!markerArrayRatios && !colorArrayRatios){PlotArray(ratioArray,controlstring,defaultColorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
    else if(!markerArrayRatios){PlotArray(ratioArray,controlstring,colorArrayRatios,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
    else if(!colorArrayRatios){PlotArray(ratioArray,controlstring,defaultColorArray,markerArrayRatios,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
    else{PlotArray(ratioArray,controlstring,colorArrayRatios,markerArrayRatios,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
    upperPad->cd();

  }else{// If no ratios are given:

    if(control.Contains("square")||control.Contains("Square")||control.Contains("SQUARE")){
      newCanvas->SetLeftMargin(margin);
      newCanvas->SetRightMargin(margin-0.05);
      newCanvas->SetTopMargin(margin-0.05);
      newCanvas->SetBottomMargin(margin+0.02);
      if(histArray->At(0)->InheritsFrom("TH2")){
//        newCanvas->SetRightMargin(margin+0.03);
        newCanvas->SetRightMargin(margin+0.06); //TODO
        newCanvas->SetTopMargin(margin-0.05);
        newCanvas->SetBottomMargin(margin+0.02);
      }
    }else if (control.Contains("horizontal")||control.Contains("Horizontal")||control.Contains("HORIZONTAL")){
      newCanvas->SetLeftMargin(margin*1.414213562);
      newCanvas->SetRightMargin(margin*1.414213562);
      newCanvas->SetTopMargin(margin - 0.03);
      newCanvas->SetBottomMargin(margin + 0.03);
    }else if (control.Contains("A4")||control.Contains("a4")){
      newCanvas->SetLeftMargin(margin);
      newCanvas->SetRightMargin(margin);
      newCanvas->SetTopMargin(margin*1.414213562 - 0.05);
      newCanvas->SetBottomMargin(margin*1.414213562 + 0.05);
    }else{
      newCanvas->SetLeftMargin(margin);
      newCanvas->SetRightMargin(margin);
      newCanvas->SetTopMargin(margin-0.1);
      newCanvas->SetBottomMargin(margin+0.1);
    }


    if(control.Contains("logX")||control.Contains("logx")||control.Contains("LogX")||control.Contains("LOGX")){
      newCanvas->SetLogx(1);}
    if(control.Contains("logY")||control.Contains("logy")||control.Contains("LogY")||control.Contains("LOGY")){
      newCanvas->SetLogy(1);}
    if(control.Contains("logZ")||control.Contains("logz")||control.Contains("LogZ")||control.Contains("LOGZ")){
      newCanvas->SetLogz(1);}
    if(control.Contains("gridY")||control.Contains("gridy")||control.Contains("GridY")||control.Contains("GRIDY")){
      newCanvas->SetGridy(1);}

    Double_t relativeTextSize;
    Double_t pad_width  = gPad->XtoPixel(gPad->GetX2());
    Double_t pad_height = gPad->YtoPixel(gPad->GetY1());
    if (pad_width < pad_height){relativeTextSize = 1/pad_width;}
    else{relativeTextSize = 1/pad_height;}
    relativeTextSize = textSizeFactor * relativeTextSize;

    if(!markerArray && !colorArray){PlotArray(histArray,controlstring,defaultColorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
    else if(!markerArray){PlotArray(histArray,controlstring,colorArray,defaultMarkerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
    else if(!colorArray){PlotArray(histArray,controlstring,defaultColorArray,markerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
    else{PlotArray(histArray,controlstring,colorArray,markerArray,relativeTextSize,lableSize, titleSize, lableFont, titleFont, titleOffsetY, titleOffsetX, margin, xMin,xMax, titleOffsetZ);}
  }

  if (histArray->At(0)->InheritsFrom("TH2")){
    TH2D *hist2d = (TH2D*) histArray->At(0);
    hist2d->GetZaxis()->SetTitleOffset(1.6);
    newCanvas->Update();
//    gStyle->SetPalette(1,0);
    TPaletteAxis *palette = (TPaletteAxis*)hist2d->GetListOfFunctions()->FindObject("palette");
//    palette->SetX1NDC(palette->GetX1NDC() * 0.8);
    palette->SetX2NDC(0.88);
//    palette->SetTitleOffset(20);
  //  palette->SetY1NDC(0.1);
  //  palette->SetY2NDC(0.5);
    newCanvas->Update();

  }

  if(ratioArray){upperPad->cd();}
  newCanvas->SetName(canvasName.c_str());
  return newCanvas;
}




///--------------------------- my stuff ----------------------------------------

string appendFileName(TObject* hist, string fileName){
  string newName = "";
  newName = newName + hist->GetName() + "_" + fileName;
  return newName;
}

string appendFileName(string name, string fileName){
  string newName = name + "_" + fileName;
  return newName;
}



TLegend* makeLegend(Double_t x, Double_t y, TObjArray* entries, string* titles, Double_t xWidth = 0.3, Bool_t noErrBars = kFALSE){

  Int_t nEntriesTrue = entries->GetEntries();
  Int_t nEntries = 0;

  for (Int_t i = 0; i < nEntriesTrue; i++){
    TString histName = entries->At(i)->GetName();
    if(titles[i] == "dummy") continue;
    if(histName.Contains("SystErr")) continue;
    if(histName.Contains("dummyHist")) continue;
    nEntries++;
  }

  Double_t yWidth = 0.05*nEntries;
  if(nEntries > 8) yWidth = 0.05*3;

  TLegend* leg = new TLegend(x, y - yWidth, x + xWidth, y);

  if(nEntries > 8) {leg->SetNColumns(5);}

  for (Int_t i = 0; i < nEntriesTrue; i++){
    TString histName = entries->At(i)->GetName();
    if(titles[i] == "dummy") continue;
    if(histName.Contains("SystErr")) continue;
    if(histName.Contains("dummyHist")) continue;
    if((entries->At(i)->InheritsFrom("TF1")) || noErrBars)
    {
      leg->AddEntry((TH1D*)entries->At(i), titles[i].c_str(), "l");
    }
    else leg->AddEntry((TH1D*)entries->At(i), titles[i].c_str(), "ep");
  }

  return leg;

}


TPaveText* makeText(Double_t x, Double_t y, string lines, Int_t nLines = 3, Double_t xWidth = 0.3){
  TPaveText* text = new TPaveText(x, y - 0.04 * nLines, x + xWidth, y, "NDC");

  unsigned long pos = 0;
  string delimiter = "<|>";
  string token;
  while ((pos = lines.find(delimiter)) != string::npos) {
      token = lines.substr(0, pos);
      text->AddText(token.c_str());
      lines.erase(0, pos + delimiter.length());
  }
  text->AddText(lines.c_str());
  text->SetTextAlign(11);
  return text;
}


string getIntegral(string histName, string fileName, TList* histList){

  Double_t nEvents = ((TH1D*)(histList->FindObject(fileName.c_str())->FindObject((appendFileName(histName.c_str(), fileName)).c_str())))->Integral();
  nEvents = floor(nEvents / 1e5)/1e1;

  ostringstream outputString;
  outputString <<  nEvents;
  outputString << "m events";


  return outputString.str();
}

void cutHist(TH1D* hist, Double_t cutoff){

    Int_t cutoffBin = hist->GetXaxis()->FindBin(cutoff);
    for(Int_t i = cutoffBin; i <= hist->GetNbinsX(); i++){
      hist->SetBinContent(i, 0);
      hist->SetBinError(i, 0);
    }
}

void setRange(string axis, Double_t low, Double_t high, TH1D* hist){

  if(axis.compare("X") == 0) hist->GetXaxis()->SetRangeUser(low, high);
  if(axis.compare("Y") == 0) hist->GetYaxis()->SetRangeUser(low, high);
  if(axis.compare("Z") == 0) hist->GetZaxis()->SetRangeUser(low, high);

}
void setTitle(string axis, string title, TH1D* hist){

  if(axis.compare("X") == 0) hist->GetXaxis()->SetTitle(title.c_str());
  if(axis.compare("Y") == 0) hist->GetYaxis()->SetTitle(title.c_str());
  if(axis.compare("Z") == 0) hist->GetZaxis()->SetTitle(title.c_str());

}


TObject* getClone(string histName, string fileName, TList* histList){
  TObject* returnObj = (histList->FindObject(fileName.c_str())->FindObject((appendFileName(histName.c_str(), fileName)).c_str()));
  if(!returnObj) cout << "ERROR: Could not find " << histName << " in " << fileName << endl;
  return returnObj->Clone();
}

TH1D* getRatio(string hist1Name, string hist2Name, string ratioName, string fileName, TList* histList, string fileName2 = ""){

  if(fileName2 == "") fileName2 = fileName;
  TH1D* hist1 = (TH1D*)getClone(hist1Name, fileName, histList);
  TH1D* hist2 = (TH1D*)getClone(hist2Name, fileName2, histList);

  hist1->Divide(hist2);
  hist1->GetYaxis()->SetTitle(ratioName.c_str());

  delete hist2;
  return hist1;
}

TList* getHistosFromFile(string fileName){

  TList* list = new TList();
  list->SetOwner();
  list->SetName(fileName.c_str());

  //TODO ohne stringstream! ueberall
  ostringstream filePath;
  filePath << "Histograms/";
  filePath <<  fileName;
  filePath << ".root";

  TFile* inputFile = TFile::Open((filePath.str()).c_str(),"READ");

  TIter next(inputFile->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {

    TNamed* obj = (TNamed*) key->ReadObj();

    if(obj->IsA() == TList::Class()) {

      TIter nextListItem((TList*)obj);
      TNamed* listItem;
      while ((listItem = (TNamed*)nextListItem())) {
	string listName = listItem->GetName();
	listItem->SetName((appendFileName(obj, listName + "_" + fileName)).c_str());
	list->Add(listItem);
      }

    }else{
      obj->SetName((appendFileName(obj, fileName)).c_str());
      list->Add(obj);
    }


  }
  return list;
}




TList* getHistosFromMacro(){

  gROOT->LoadMacro("Publication/meanpt_pub.C");

  TList* list = new TList();
  list->SetOwner();
  list->SetName("pub");

  TH1D* momentReweighted1SysPP = getHistogramMeanPtPP7TeVSys();
  momentReweighted1SysPP->SetName("momentReweighted1Sys_pp_7TeV_03_pub");
  list->Add(momentReweighted1SysPP);

  TH1D* momentReweighted1SysPPb = getHistogramMeanPtPPb5TeVSys();
  momentReweighted1SysPPb->SetName("momentReweighted1Sys_pPb_5TeV_03_pub");
  list->Add(momentReweighted1SysPPb);

  TH1D* momentReweighted1SysPbPb = getHistogramMeanPtPbPb2TeVSys();
  momentReweighted1SysPbPb->SetName("momentReweighted1Sys_PbPb_2TeV_03_pub");
  list->Add(momentReweighted1SysPbPb);

  return list;
}


/// ---------------------------------------------------------------------------
/// Function to obtain moments from momentum distributions
/// Assumes Pt distribution to be the invariant yield.
/// ---------------------------------------------------------------------------
TH1D* getMomentFromSpectra(TH2D* multPt, Int_t moment){

  TH1D* multMoment = (TH1D*) multPt->ProjectionX("dummyName");
  multMoment->Reset();

  for (Int_t multBin = multPt->GetXaxis()->FindBin(1) ; multBin <= multPt->GetNbinsX() ; multBin++){

    Double_t numerator = 0;
    Double_t denominator = 0;
    Double_t ptBinCenter = 0;
    Double_t ptBinWidth = 0;
    Double_t yield = 0;
    Double_t yieldError = 0;
    Double_t weightedQuantity = 0;

    // calculate moment (startbin at Pt > 150MeV)
    for (Int_t ptBin = multPt->GetYaxis()->FindBin(0.1501); ptBin <= multPt->GetNbinsY() ; ptBin++){

      ptBinCenter = multPt->GetYaxis()->GetBinCenter(ptBin);
      ptBinWidth = multPt->GetYaxis()->GetBinWidth(ptBin);
      yield = multPt->GetBinContent(multBin, ptBin) * ptBinCenter ;

      weightedQuantity = 1.;
      for(Int_t i = 0; i < moment; i++) weightedQuantity *= ptBinCenter;

      numerator += yield * weightedQuantity * ptBinWidth;
      denominator += yield * ptBinWidth;
    }
    Double_t mean =  numerator / denominator;

    // calculate errors
    Double_t errSum = 0;
    for (Int_t ptBin = multPt->GetYaxis()->FindBin(0.1501); ptBin <= multPt->GetNbinsY() ; ptBin++){

      ptBinCenter = multPt->GetYaxis()->GetBinCenter(ptBin);
      ptBinWidth = multPt->GetYaxis()->GetBinWidth(ptBin);
      yieldError = multPt->GetBinError(multBin, ptBin) * ptBinCenter;

      weightedQuantity = 1.;
      for(Int_t i = 0; i < moment; i++) weightedQuantity *= ptBinCenter;

      errSum += TMath::Power((weightedQuantity - mean) * yieldError * ptBinWidth, 2);
    }

    // fill results in histogram
    if(denominator){
      multMoment->SetBinContent(multBin, mean);
      Double_t error = 1/denominator * TMath::Sqrt(errSum);
      multMoment->SetBinError(multBin, error);
    }
  }
  return multMoment;
}


/// ---------------------------------------------------------------------------
/// Function to normalize response Matrix
/// ---------------------------------------------------------------------------
TH2D* normalize(TH2D* matrix){

  TH2D* normalizedMatrix = (TH2D*)matrix->Clone("responseMatrix");
  normalizedMatrix->Reset();

  for (Int_t y = 2 ; y <= matrix->GetNbinsY() ; y++){

    Double_t yIntegral = 0;
    for (Int_t x = 2 ; x <= matrix->GetNbinsX() ; x++) yIntegral += matrix->GetBinContent(x, y);

    for (Int_t x = 2 ; x <= matrix->GetNbinsX() ; x++){

      Double_t value = matrix->GetBinContent(x, y);
      if (value) normalizedMatrix->SetBinContent(x, y, value/yIntegral);

    }
  }
  return normalizedMatrix;
}

TH1D* getRebinRatio(TH1D* fineBinnedHist, TH1D* coarseBinnedHist){

  TH1D* original = (TH1D*)coarseBinnedHist->Clone("original");
  original->Reset();

  for(Int_t i = 1; i <= original->GetNbinsX(); i++){

    Double_t xValue = original->GetXaxis()->GetBinCenter(i);
    Int_t bin = fineBinnedHist->FindBin(xValue);

    Double_t content = fineBinnedHist->GetBinContent(bin);
    Double_t error = fineBinnedHist->GetBinError(bin);
    original->SetBinContent(i, content);
    original->SetBinError(i, error);
  }

  TH1D* ratio = (TH1D*)coarseBinnedHist->Clone("ratio");
  ratio->Divide(original);
  delete original;
  return ratio;
}


TH2D* getRelativeErrorRatio(TH2D* efficiency2D){

  TH2D* relErrors = (TH2D*)efficiency2D->Clone("relErrors");
  relErrors->Clear();

  for (Int_t x = 1 ; x <= efficiency2D->GetNbinsX() ; x++){
    for (Int_t y = 1 ; y <= efficiency2D->GetNbinsY() ; y++){
      Double_t content = efficiency2D->GetBinContent(x,y);
      if(!content) continue;
      Double_t error = efficiency2D->GetBinError(x,y)/content;
      relErrors->SetBinContent(x, y, error);
    }
  }
  return relErrors;
}


/// ---------------------------------------------------------------------------
/// Function to convert Nch on xAxis to dNch/dEta
///
/// ---------------------------------------------------------------------------
TH1D* compressAxis(TH1D* inputSpectrum, Double_t factor){

  ostringstream stringStream;
  stringStream << inputSpectrum->GetName() << " dnchdeta";
  string name = stringStream.str();

  Double_t lowerBound = factor * inputSpectrum->GetBinLowEdge(1);
  Double_t lastBin = inputSpectrum->GetNbinsX();
  Double_t upperBound = factor * inputSpectrum->GetBinLowEdge(lastBin + inputSpectrum->GetBinWidth(lastBin));

  TH1D* outputSpectrum = new TH1D(name.c_str(), name.c_str(), inputSpectrum->GetNbinsX(), lowerBound, upperBound);
  outputSpectrum->GetXaxis()->SetTitle("d#it{N}_{ch} / d#it{#eta}");
  outputSpectrum->GetYaxis()->SetTitle(inputSpectrum->GetYaxis()->GetTitle());

  for(Int_t i = 1; i <= lastBin; i++){
    outputSpectrum->SetBinContent(i, inputSpectrum->GetBinContent(i));
    outputSpectrum->SetBinError(i, inputSpectrum->GetBinError(i));
  }

  return outputSpectrum;
}
