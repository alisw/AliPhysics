/************************************************************
     -- provided by Gamma Conversion Group, PWG4,
     -- Kathrin Koch, kkoch@physi.uni-heidelberg.de
     -- Friederike Bock, friederike.bock@cern.ch   
                                                          
 ************************************************************

 == This macro can be used to display the Photon 
     Characteristics of the conversion method in ALICE, 
     it can be operated  *****

 == on the output of the GammaConversionTask. It can take 
     2 input files, the second one should be MC, if this 
     is not    *****

 == the case all histograms including MC need to be 
     commented out otherwise the running will crash.  *****
 *************************************************************/

/*----------------------------------------------
A small Modificaiton is done by sjena to:
 -- impliment the unique name of each object
 -- output files for comparision purpose
 ----------------------------------------------*/



#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include "TGaxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 

TString textGenerator;
TString collisionSystem;
TString textPeriod;
TString textDate;

void SetStyleTLatex( TLatex* text, Size_t textSize, Width_t lineWidth, Color_t textColor = 1, Bool_t kNDC = kTRUE){
   if (kNDC) {text->SetNDC();}
   text->SetTextColor(textColor);
   text->SetTextSize(textSize);
   text->SetLineWidth(lineWidth);
}


/* DrawAutoGammaHisto is function used for styling a histograma of the gamma conversion group with standart settings
* histo1 - first histogram (Data)
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
* YRangeMax    = kTRUE will scale by Maximum and Minimum Range in Y
*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
*YMinimum - this will be used if YRangeMax is set
*YRange     = kTRUE will Cut y-axis by YMin and YMax 
- will be set to kFAlSE if YRangeMax is set
*YMin - minimum Y
*YMax - maximum Y
*XRange  = kTRUE will Cut x-axis by XMin and XMax
*XMin - minimum Y
*XMax - maximum Y
*/ 
void DrawAutoGammaHisto( TH1* histo1, 
               TString Title, TString XTitle, TString YTitle,
               Bool_t YRangeMax, Float_t YMaxFactor, Float_t YMinimum,
               Bool_t YRange, Float_t YMin ,Float_t YMax,  
               Bool_t XRange, Float_t XMin, Float_t XMax) {
   if (YRangeMax && !XRange){
      YRange = kFALSE;
      Double_t maxRangeR = histo1->GetMaximum();
      Double_t minRangeR = histo1->GetMinimum();      
      if(YMinimum > minRangeR){minRangeR = YMinimum;}
      histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);   
   }
   if (YRangeMax && XRange){
      YRange = kFALSE;
      Double_t maxRangeR = histo1->GetMaximum();
      Double_t minRangeR = histo1->GetMinimum();      
      if(YMinimum > minRangeR){minRangeR = YMinimum;}
      histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);   
      histo1->GetXaxis()->SetRangeUser(XMin, XMax);   
   }
   if (YRange && XRange){
      histo1->GetYaxis()->SetRangeUser(YMin, YMax);   
      histo1->GetXaxis()->SetRangeUser(XMin, XMax);   
   }
   if (!YRangeMax && !YRange && XRange){
      histo1->GetXaxis()->SetRangeUser(XMin, XMax);   
   }
   
   if (YRange && !XRange){
      histo1->GetYaxis()->SetRangeUser(YMin, YMax);
   }
   
   histo1->SetTitle(Title.Data());
   
   if(XTitle.CompareTo("") != 0){
      histo1->SetXTitle(XTitle.Data());
   }
   if(YTitle.CompareTo("") != 0){
      histo1->SetYTitle(YTitle.Data());
   }
   histo1->GetYaxis()->SetLabelSize(0.03);
   histo1->GetYaxis()->SetTitleSize(0.035);  
   histo1->GetYaxis()->SetDecimals();
   histo1->GetYaxis()->SetTitleOffset(1.8);
   histo1->GetXaxis()->SetTitleSize(0.035);
   histo1->GetXaxis()->SetLabelSize(0.03);   
   histo1->SetLineColor(kBlue+2);
   histo1->SetMarkerColor(kBlue+2);
   histo1->SetMarkerStyle(kFullCircle);
   histo1->SetMarkerSize(1.5);
   histo1->DrawCopy("e,p");
}


void DrawLabelsEvents(Float_t startX, Float_t startY, Float_t textHeight, Float_t decrease,  TString collisionSystemDummy, TString textGeneratorDummy, TString textPeriodDummy){
   
   Float_t aliceStartY = startY - textHeight * 1.15;  
   TLatex *pp7 = NULL;
   if( collisionSystemDummy.CompareTo("PbPb @ #sqrt{s_{NN}} = 2.76 TeV") == 0){
      pp7 = new TLatex((startX-2*decrease),(aliceStartY),collisionSystemDummy.Data()); // Bo: this was modified
   } else {
      pp7 = new TLatex((startX+2*decrease),(aliceStartY),collisionSystemDummy.Data()); // Bo: this was modified
   }
   pp7->SetNDC();
   pp7->SetTextColor(1);
   pp7->SetTextFont(62);   
   pp7->SetTextSize(textHeight);
   pp7->SetLineWidth(2);
   pp7->Draw("same");
   if (textGeneratorDummy.CompareTo("")!=0 && textPeriodDummy.CompareTo("")!=0){
      TLatex *generator = new TLatex((startX+decrease),(aliceStartY-1*textHeight*1.15),Form("%s   %s",textGeneratorDummy.Data(),textPeriodDummy.Data())); // Bo: this was modified
      generator->SetNDC();
      generator->SetTextColor(1);
      generator->SetTextFont(62);
      generator->SetTextSize(textHeight);
      generator->SetLineWidth(2);
      generator->Draw("same");   
   } else if (textGeneratorDummy.CompareTo("")!=0) {
      TLatex *generator = new TLatex((startX+decrease),(aliceStartY-1*textHeight*1.15),Form("%s",textGeneratorDummy.Data())); // Bo: this was modified
      generator->SetNDC();
      generator->SetTextColor(1);
      generator->SetTextFont(62);
      generator->SetTextSize(textHeight);
      generator->SetLineWidth(2);
      generator->Draw("same");   
   }
}

/*DrawAutoGammaHisto2D is a function for drawing a 2D-histogram of the gamma conversion group
* histo - histogramm which need to be drawn
* Title - histogram title
* XTitle - X- axis-title
* YTitle - Y-axis-title
* Input - Legend 
* YRange - if kTRUE will scale by YMin and YMay
* YMin  - Y minimum
* YMax - Y maximum
* XRange - if kTRUE will scale by XMin and XMax
* XMin - X minimum
* XMax - X maximum
*/
void DrawAutoGammaHisto2D( TH2 *histo,  
                  TString Title, TString XTitle, TString YTitle, TString Input,
                  Bool_t YRange, Float_t YMin ,Float_t YMax, 
                  Bool_t XRange, Float_t XMin, Float_t XMax,Float_t titleOffsetX=1.4, Float_t titleOffsetY=1.2) {
   
   
   if (YRange && XRange){
      histo->GetYaxis()->SetRangeUser(YMin, YMax); 
      histo->GetXaxis()->SetRangeUser(XMin, XMax); 
   }
   if ( !YRange && XRange){
      histo->GetXaxis()->SetRangeUser(XMin, XMax); 
   }
   
   if (YRange && !XRange){
      histo->GetYaxis()->SetRangeUser(YMin, YMax);
   }
   
//    if(Title.CompareTo("") != 0){
      histo->SetTitle(Title.Data());
//    }
   if(XTitle.CompareTo("") != 0){
      histo->SetXTitle(XTitle.Data());
   }
   if(YTitle.CompareTo("") != 0){
      histo->SetYTitle(YTitle.Data());
   }
   histo->GetYaxis()->SetTitleSize(0.043);   
   histo->GetYaxis()->SetLabelSize(0.035);
   histo->GetXaxis()->SetLabelSize(0.035);
   histo->GetYaxis()->SetDecimals();
   histo->GetYaxis()->SetTitleOffset(titleOffsetY);
   histo->GetXaxis()->SetTitleOffset(titleOffsetX);
   histo->GetXaxis()->SetTitleSize(0.043);   
   histo->DrawCopy("colz");
   if(Input.CompareTo("") != 0){
      TLegend* leg2 = new TLegend(0.6,0.82,0.83,0.9);
      leg2->SetTextSize(0.04);         
      leg2->SetFillColor(0);
      leg2->AddEntry(histo,(Input.Data()));
      leg2->Draw("same");
   }
}

void PlotStandard2D( TH2* histo2D, TString nameOutput, TString title, TString xTitle, TString yTitle, Bool_t kRangeY, Double_t startY, Double_t endY, Bool_t kRangeX, Double_t startX, Double_t endX, Int_t logX, Int_t logY, Int_t logZ, Float_t* floatLogo, Int_t canvasSizeX = 500, Int_t canvasSizeY = 500, TString generator ="" , TString period =""){
   TCanvas * canvasStandard = new TCanvas("canvasStandard","",10,10,canvasSizeX,canvasSizeY);  // gives the page size      
   canvasStandard->SetLogx(logX);
   canvasStandard->SetLogy(logY);
   canvasStandard->SetLogz(logZ);
   canvasStandard->SetRightMargin(0.12);     
   canvasStandard->SetLeftMargin(0.12);      
   canvasStandard->SetBottomMargin(0.1);     
   canvasStandard->SetTopMargin(0.04);       
   canvasStandard->cd();
   histo2D->SetTitle("");
   DrawAutoGammaHisto2D(   histo2D,
                           title.Data(), xTitle.Data(), yTitle.Data(),"",kRangeY, startY, endY, kRangeX, startX, endX);
   histo2D->GetXaxis()->SetTitleOffset(1.05);
//    cout << histo2D->GetYaxis()->GetTitleOffset() << endl;
   histo2D->GetYaxis()->SetTitleOffset(1.35);
   if (logX==1){
//       cout << histo2D->GetXaxis()->GetLabelOffset() << endl;
      histo2D->GetXaxis()->SetLabelOffset(0.);
   }   
      
   histo2D->Draw("colz");
   DrawLabelsEvents(floatLogo[0],floatLogo[1],floatLogo[2], 0.00, collisionSystem, generator, period);
   
   canvasStandard->Update();
   canvasStandard->SaveAs(nameOutput.Data());
   delete canvasStandard;
}

TString GetCentralityString(TString cutNumber){
   TString centralityCutNumberStart = cutNumber(1,1);
   TString centralityCutNumberEnd = cutNumber(2,1);
   TString ppCutNumber = cutNumber(0,1);
   if (ppCutNumber.CompareTo("0") ==0){
           return "pp"; 
   } else if ( ppCutNumber.CompareTo("1") ==0 || ppCutNumber.CompareTo("2") ==0 || ppCutNumber.CompareTo("5") ==0 || ppCutNumber.CompareTo("8") ==0 || ppCutNumber.CompareTo("9") ==0){       
      if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
          return "0-100%"; 
      } else {
          return Form("%i-%i%s", centralityCutNumberStart.Atoi()*10,centralityCutNumberEnd.Atoi()*10,"%");
      }
   } else if (ppCutNumber.CompareTo("3") ==0 || ppCutNumber.CompareTo("6") ==0){
      if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
          return "0-45%"; 
      } else {
          return Form("%i-%i%s", centralityCutNumberStart.Atoi()*5,centralityCutNumberEnd.Atoi()*5,"%");
      }
   } else if (ppCutNumber.CompareTo("4") ==0 || ppCutNumber.CompareTo("7") ==0){
      if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
          return "45-95%"; 
      } else {
          return Form("%i-%i%s",45+centralityCutNumberStart.Atoi()*5,45+centralityCutNumberEnd.Atoi()*5,"%");
      }
   } else return ""; 
}  

TString ReturnFullCollisionsSystem(TString fEnergyFlagOpt){ 
   if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
      return  "pp, #sqrt{s} = 7 TeV";
   } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
      return  "pp, #sqrt{s} = 900 GeV";
   } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
      return  "pp, #sqrt{s} = 2.76 TeV";
   } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
      return "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV";
   } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
      return "p-Pb, #sqrt{s_{NN}} = 5.023 TeV";
   } else {
      cout << "No correct collision system specification, has been given" << endl;
      return "";
   }
}

/* StyleSettingsThesis will make some standard settings for gStyle 
*/
void StyleSettingsThesis(){
   gStyle->SetOptDate(0);   //show day and time
   gStyle->SetOptStat(0);  //show statistic
   gStyle->SetPalette(1,0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameFillColor(0);
   gStyle->SetTitleFillColor(0);
   gStyle->SetTextSize(0.5);
   gStyle->SetLabelSize(0.03,"xyz");
   gStyle->SetLabelOffset(0.002,"xyz");
   gStyle->SetTitleFontSize(0.04);
   gStyle->SetTitleOffset(1,"y");
   gStyle->SetTitleOffset(0.7,"x");    
   gStyle->SetCanvasColor(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   //gStyle->SetLineWidth(0.01);
   
   gStyle->SetPadTopMargin(0.03);
   gStyle->SetPadBottomMargin(0.09);
   gStyle->SetPadRightMargin(0.03);
   gStyle->SetPadLeftMargin(0.13);
   
   
   TGaxis::SetMaxDigits(5);
   gErrorIgnoreLevel=kError;
}


void SetPlotStyle() {
   const Int_t nRGBs = 5;
   const Int_t nCont = 255;
   
   Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   Double_t green[nRGBs] = { 0.31, 0.81, 1.00, 0.20, 0.00 };
   Double_t blue[nRGBs]  = { 0.51, 1., 0.12, 0.00, 0.00};

   TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
   gStyle->SetNumberContours(nCont);
}

               
void processProduceFastQA(TString fileNameInput = "myOutput", 
			  TString cutSelection = "5080001022092970023220000000", 
			  TString suffix = "eps", 
			  TString optEnergy="", 
			  TString optMCGenerator="", 
			  TString optPeriod="", 
			  const char* outfile="ProduceFastQA_output.root"){   
   
   gROOT->Reset();   
   gSystem->Load("libCore");
   gSystem->Load("libTree");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC");
   gSystem->Load("libPhysics");
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libCORRFW");
   gROOT->SetStyle("Plain");
   
   StyleSettingsThesis();   
   SetPlotStyle();
   
   collisionSystem = ReturnFullCollisionsSystem(optEnergy);
   if (collisionSystem.CompareTo("") == 0){
      cout << "No correct collision system specification, has been given" << endl;
      return;
   }
   TString centralityCutNumber = cutSelection(0,3);
   TString centralityString = GetCentralityString(centralityCutNumber);
   cout<< centralityCutNumber.Data() << "\t" << centralityString.Data() << endl;
   if (centralityString.CompareTo("pp") !=0){
      collisionSystem = Form("%s %s", centralityString.Data(), collisionSystem.Data());
   }

   TString textGenerator;

   if(optMCGenerator.CompareTo("") ==0){
      textGenerator = "";
   } else {
      textGenerator = optMCGenerator;
   }
   
   TFile*  fileInput = new TFile(fileNameInput.Data());
fileInput->ls();

   TDirectory* directoryQA = (TDirectory*)fileInput->Get(Form("GammaConvV1_QA_%s",cutSelection.Data())); 
directoryQA->ls();
   TList* listQA = (TList*)directoryQA->Get(Form("GammaConv_V1QA_%s",cutSelection.Data()));
listQA->ls();
   TList* listQAESD = (TList*)listQA->FindObject("ESD QA");


   // Added by sjena
   TFile *fout = TFile::Open(outfile,"UPDATE");
   fout->ls();
   
   TDirectoryFile *cdd = NULL;
   cdd = (TDirectoryFile*)fout->Get("GA");
   if(!cdd) {
     Printf("Warning: GA <dir> doesn't exist, creating a new one");
     cdd = (TDirectoryFile*)fout->mkdir("GA");
   }
   cdd->cd();
   cdd->ls();
   

   TH1F* histoVertexZ =             (TH1F*)listQAESD->FindObject("Vertex_Z");
   histoVertexZ->Sumw2();
   Double_t nEvt = histoVertexZ->GetEntries();
   histoVertexZ->Scale(1./nEvt);
   histoVertexZ->Write(Form("fig_ga_%s",histoVertexZ->GetName()));
   
   TH1I* histoContrVertexZ =         (TH1I*)listQAESD->FindObject("ContrVertex_Z");
   TH1D* histoDContrVertexZ = new TH1D("ContrVertex_Z","ContrVertex_Z",3000,0,3000);
   histoDContrVertexZ->Sumw2();
   for (Int_t i = 1; i < histoContrVertexZ->GetNbinsX(); i++){
      histoDContrVertexZ->SetBinContent(i, histoContrVertexZ->GetBinContent(i)/nEvt);
      histoDContrVertexZ->SetBinError(i, histoContrVertexZ->GetBinError(i)/nEvt);
   }   
   histoDContrVertexZ->Rebin(8);
   histoDContrVertexZ->Write(Form("fig_ga_%s",histoDContrVertexZ->GetName()));
			     
   TH1I* histoGoodESDTracks = (TH1I*)listQAESD->FindObject("GoodESDTracks");
   TH1D* histoDGoodESDTracks = new TH1D("GoodESDTracks","GoodESDTracks",3000,0,3000);
   histoDGoodESDTracks->Sumw2();
   for (Int_t i = 1; i < histoGoodESDTracks->GetNbinsX(); i++){
     histoDGoodESDTracks->SetBinContent(i, histoGoodESDTracks->GetBinContent(i)/nEvt);
     histoDGoodESDTracks->SetBinError(i, histoGoodESDTracks->GetBinError(i)/nEvt);
   }   
   histoDGoodESDTracks->Rebin(8);
   histoDGoodESDTracks->Write(Form("fig_ga_%s",histoDGoodESDTracks->GetName()));
   
   
   TH1I* histoV0Mult = (TH1I*)listQAESD->FindObject("V0 Multiplicity");
   TH1D* histoDV0Mult = new TH1D("V0Multiplicity","V0 Multiplicity",30000,0,30000);
   histoDV0Mult->Sumw2();
   for (Int_t i = 1; i < histoV0Mult->GetNbinsX(); i++){
      histoDV0Mult->SetBinContent(i, histoV0Mult->GetBinContent(i)/nEvt);
      histoDV0Mult->SetBinError(i, histoV0Mult->GetBinError(i)/nEvt);
   }   
   histoDV0Mult->Rebin(10);
   histoDV0Mult->Write(Form("fig_ga_%s",histoDV0Mult->GetName()));
   
      
   TH2F* histoITSClusterPhi = (TH2F*)listQAESD->FindObject("ITSClusterPhi");
   histoITSClusterPhi->Sumw2();
   histoITSClusterPhi->Scale(1./nEvt);
   histoITSClusterPhi->Write(Form("fig_ga_%s",histoITSClusterPhi->GetName()));   


   TH1F* histoGammaPt = (TH1F*)listQAESD->FindObject("Gamma_Pt");
   histoGammaPt->Sumw2();
   histoGammaPt->Scale(1./nEvt);
   histoGammaPt->Write(Form("fig_ga_%s",histoGammaPt->GetName()));
   
   
   TH1F* histoGammaPhi = (TH1F*)listQAESD->FindObject("Gamma_Phi");
   histoGammaPhi->Sumw2();
   histoGammaPhi->Scale(1./nEvt);
   histoGammaPhi->Rebin(2);
   histoGammaPhi->Write(Form("fig_ga_%s",histoGammaPhi->GetName()));   


   TH1F* histoGammaEta = (TH1F*)listQAESD->FindObject("Gamma_Eta");
   histoGammaEta->Sumw2();
   histoGammaEta->Scale(1./nEvt);
   histoGammaEta->Rebin(2);
   histoGammaEta->Write(Form("fig_ga_%s",histoGammaEta->GetName()));

   
   TH1F* histoGammaChi2 = (TH1F*)listQAESD->FindObject("Gamma_Chi2perNDF");
   histoGammaChi2->Sumw2();
   histoGammaChi2->Scale(1./nEvt);
   histoGammaChi2->Write(Form("fig_ga_%s",histoGammaChi2->GetName()));


   
   TH1F* histoGammaPsiPair = (TH1F*)listQAESD->FindObject("Gamma_PsiPair");
   histoGammaPsiPair->Sumw2();
   histoGammaPsiPair->Scale(1./nEvt);
   histoGammaPsiPair->Write(Form("fig_ga_%s",histoGammaPsiPair->GetName()));   


   TH1F* histoGammaCosPoint = (TH1F*)listQAESD->FindObject("Gamma_CosinePointingAngle");
   histoGammaCosPoint->Sumw2();
   histoGammaCosPoint->Scale(1./nEvt);
   histoGammaCosPoint->Write(Form("fig_ga_%s",histoGammaCosPoint->GetName()));


   TH1F* histoGammaInvMass = (TH1F*)listQAESD->FindObject("Gamma_InvMass");
   histoGammaInvMass->Sumw2();
   histoGammaInvMass->Scale(1./nEvt);
   histoGammaInvMass->Write(Form("fig_ga_%s",histoGammaInvMass->GetName()));   
   

   TH2F* histoGammaArmenteros = (TH2F*)listQAESD->FindObject("Gamma_Armenteros");
   histoGammaArmenteros->Sumw2();
   histoGammaArmenteros->Scale(1./nEvt);
   histoGammaArmenteros->Write(Form("fig_ga_%s",histoGammaArmenteros->GetName()));   
   

   
   TH2F* histoEPPt = (TH2F*)listQAESD->FindObject("Electron_Positron_Pt");
   histoEPPt->Sumw2();
   histoEPPt->Scale(1./nEvt);
   histoEPPt->Write(Form("fig_ga_%s",histoEPPt->GetName()));


   TH2F* histoEPEta = (TH2F*)listQAESD->FindObject("Electron_Positron_Eta");
   histoEPEta->Sumw2();
   histoEPEta->Scale(1./nEvt);
   histoEPEta->Write(Form("fig_ga_%s",histoEPEta->GetName()));   


   TH2F* histoEPPhi = (TH2F*)listQAESD->FindObject("Electron_Positron_Phi");
   histoEPPhi->Sumw2();
   histoEPPhi->Scale(1./nEvt);
   histoEPPhi->Write(Form("fig_ga_%s",histoEPPhi->GetName()));   


   TH1F* histoEFindableClusterTPC = (TH1F*)listQAESD->FindObject("Electron_findableClusterTPC");
   histoEFindableClusterTPC->Sumw2();
   histoEFindableClusterTPC->Scale(1./nEvt);
   histoEFindableClusterTPC->Write(Form("fig_ga_%s",histoEFindableClusterTPC->GetName()));   


   TH1F* histoPFindableClusterTPC = (TH1F*)listQAESD->FindObject("Positron_findableClusterTPC");
   histoPFindableClusterTPC->Sumw2();
   histoPFindableClusterTPC->Scale(1./nEvt);
   histoPFindableClusterTPC->Write(Form("fig_ga_%s",histoPFindableClusterTPC->GetName()));   

   
   TH2F* histoEdEdxPTPC = (TH2F*)listQAESD->FindObject("Electron_dEdx_P");
   Double_t nElectronTPC = histoEdEdxPTPC->GetEntries();
   histoEdEdxPTPC->Sumw2();
   histoEdEdxPTPC->Scale(1./nEvt);
   histoEdEdxPTPC->Write(Form("fig_ga_%s",histoEdEdxPTPC->GetName()));
   
   
   TH2F* histoPdEdxPTPC = (TH2F*)listQAESD->FindObject("Positron_dEdx_P");
   Double_t nPositronTPC = histoPdEdxPTPC->GetEntries();
   histoPdEdxPTPC->Sumw2();
   histoPdEdxPTPC->Scale(1./nEvt);
   histoPdEdxPTPC->Write(Form("fig_ga_%s",histoPdEdxPTPC->GetName()));
   
   TH2F* histoENSigmadEdxPTPC = (TH2F*)listQAESD->FindObject("Electron_NSigmadEdx_P");
   histoENSigmadEdxPTPC->Sumw2();
   histoENSigmadEdxPTPC->Scale(1./nEvt);
   histoENSigmadEdxPTPC->Write(Form("fig_ga_%s",histoENSigmadEdxPTPC->GetName()));

   TH2F* histoPNSigmadEdxPTPC = (TH2F*)listQAESD->FindObject("Positron_NSigmadEdx_P");
   histoPNSigmadEdxPTPC->Sumw2();
   histoPNSigmadEdxPTPC->Scale(1./nEvt);
   histoPNSigmadEdxPTPC->Write(Form("fig_ga_%s",histoPNSigmadEdxPTPC->GetName()));

   TH2F* histoENSigmaPiondEdxPTPC = (TH2F*)listQAESD->FindObject("Electron_NSigmaPiondEdx_P");
   histoENSigmaPiondEdxPTPC->Sumw2();
   histoENSigmaPiondEdxPTPC->Scale(1./nEvt);
   histoENSigmaPiondEdxPTPC->Write(Form("fig_ga_%s",histoENSigmaPiondEdxPTPC->GetName()));
   
   
   TH2F* histoPNSigmaPiondEdxPTPC = (TH2F*)listQAESD->FindObject("Positron_NSigmaPiondEdx_P");
   histoPNSigmaPiondEdxPTPC->Sumw2();
   histoPNSigmaPiondEdxPTPC->Scale(1./nEvt);
   histoPNSigmaPiondEdxPTPC->Write(Form("fig_ga_%s",histoPNSigmaPiondEdxPTPC->GetName()));
   
   TH2F* histoETOFP = (TH2F*)listQAESD->FindObject("Electron_TOF_P");
   Double_t nElectronTOF = histoETOFP->GetEntries();
   histoETOFP->Sumw2();
   histoETOFP->Scale(1./nEvt);
   histoETOFP->Write(Form("fig_ga_%s",histoETOFP->GetName()));


   TH2F* histoPTOFP = (TH2F*)listQAESD->FindObject("Positron_TOF_P");
   Double_t nPositronTOF = histoPTOFP->GetEntries();
   histoPTOFP->Sumw2();
   histoPTOFP->Scale(1./nEvt);
   histoPTOFP->Write(Form("fig_ga_%s",histoPTOFP->GetName()));
   
   TH2F* histoENSigmaTOFP = (TH2F*)listQAESD->FindObject("Electron_NSigmaTOF_P");
   histoENSigmaTOFP->Sumw2();
   histoENSigmaTOFP->Scale(1./nEvt);
   histoENSigmaTOFP->Write(Form("fig_ga_%s",histoENSigmaTOFP->GetName()));

   TH2F* histoPNSigmaTOFP = (TH2F*)listQAESD->FindObject("Positron_NSigmaTOF_P");
   histoPNSigmaTOFP->Sumw2();
   histoPNSigmaTOFP->Scale(1./nEvt);
   histoPNSigmaTOFP->Write(Form("fig_ga_%s",histoPNSigmaTOFP->GetName()));   

   TH2F* histoEdEdxPITS = (TH2F*)listQAESD->FindObject("Electron_ITSdEdx_P");
   Double_t nElectronITS = histoEdEdxPITS->GetEntries();
   histoEdEdxPITS->Sumw2();
   histoEdEdxPITS->Scale(1./nEvt);
   histoEdEdxPITS->Write(Form("fig_ga_%s",histoEdEdxPITS->GetName()));

   TH2F* histoPdEdxPITS = (TH2F*)listQAESD->FindObject("Positron_ITSdEdx_P");
   Double_t nPositronITS = histoPdEdxPITS->GetEntries();
   histoPdEdxPITS->Sumw2();
   histoPdEdxPITS->Scale(1./nEvt);
   histoPdEdxPITS->Write(Form("fig_ga_%s",histoPdEdxPITS->GetName()));

   TH2F* histoENSigmadEdxPITS = (TH2F*)listQAESD->FindObject("Electron_NSigmaITS_P");
   histoENSigmadEdxPITS->Sumw2();
   histoENSigmadEdxPITS->Scale(1./nEvt);
   histoENSigmadEdxPITS->Write(Form("fig_ga_%s",histoENSigmadEdxPITS->GetName()));

   TH2F* histoPNSigmadEdxPITS = (TH2F*)listQAESD->FindObject("Positron_NSigmaITS_P");
   histoPNSigmadEdxPITS->Sumw2();
   histoPNSigmadEdxPITS->Scale(1./nEvt);
   histoPNSigmadEdxPITS->Write(Form("fig_ga_%s",histoPNSigmadEdxPITS->GetName()));
   
   //fout->Close();


   TCanvas * canvasEventProp = new TCanvas("canvasEventProp","",0,0,1000,1000);  // gives the page size
   TPad* padEventProp = new TPad("padEventProp","",0.0,0.0,1,1,0);   // gives the size of the histo areas 
   padEventProp->SetFillColor(0);
   padEventProp->GetFrame()->SetFillColor(0);
   padEventProp->SetBorderMode(0);
   
   padEventProp->Divide(2,2);
   padEventProp->Draw();
   padEventProp->cd(1);
   
   DrawAutoGammaHisto( histoVertexZ,
                       "", "Z_{vtx} (cm)","dZ/dN_{evt}",
                       kTRUE, 1.2, 0.,
                       kFALSE, -10,10,
                       kFALSE, -10,10);
   
   
   TLatex *labelDataSet = NULL;
   if (optPeriod.CompareTo("") ){
      labelDataSet = new TLatex(0.18,0.9,Form("%s",optPeriod.Data()));
      SetStyleTLatex( labelDataSet, 0.05,4);
      labelDataSet->Draw();
   }
   padEventProp->cd(2);
   padEventProp->cd(2)->SetLogy(1);
   
   DrawAutoGammaHisto( histoDContrVertexZ,
                       "", "# Contr to prim Vtx","norm counts",
                       kTRUE, 2., 0.5/nEvt,
                       kFALSE, -10,10,
                       kFALSE, -10,10);
    
   padEventProp->cd(3);
   padEventProp->cd(3)->SetLogy(1);
   
   DrawAutoGammaHisto( histoDGoodESDTracks,
                       "", "# Good ESD tracks","norm counts",
                       kTRUE, 2., 0.5/nEvt,
                       kFALSE, -10,10,
                       kFALSE, -10,10);
   
   padEventProp->cd(4);
   padEventProp->cd(4)->SetLogy(1);
   
   DrawAutoGammaHisto( histoDV0Mult,
                       "", "V0 signal","norm counts",
                       kTRUE, 2., 0.5/nEvt,
                       kFALSE, -10,10,
                       kFALSE, -10,10);
   
   canvasEventProp->Update();
   canvasEventProp->SaveAs(Form("fig_ga_EventCharacteristics.%s",suffix.Data()));

   TCanvas * canvasdEdxTPC = new TCanvas("canvasdEdxTPC","",0,0,1000,1500);  // gives the page size
   TPad* paddEdxTPC = new TPad("paddEdxTPC","",0.0,0.0,1,1,0);   // gives the size of the histo areas 
   paddEdxTPC->SetFillColor(0);
   paddEdxTPC->GetFrame()->SetFillColor(0);
   paddEdxTPC->SetBorderMode(0);
   
   paddEdxTPC->Divide(2,3);
   paddEdxTPC->Draw();
   paddEdxTPC->cd(1);
   paddEdxTPC->cd(1)->SetLogx(1);
   paddEdxTPC->cd(1)->SetLogz(1);
   paddEdxTPC->cd(1)->SetTopMargin(0.01);
   paddEdxTPC->cd(1)->SetRightMargin(0.12);
   Double_t maximumEPTPC = 1.2*histoEdEdxPTPC->GetMaximum();
   histoEdEdxPTPC->GetZaxis()->SetRangeUser(1/nEvt,maximumEPTPC);
   histoEdEdxPTPC->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoEdEdxPTPC,  
                  "", "#it{p}_{e^{-}} (GeV/c)", "dE_{e^{-}-cand} /dx",  "",
                  kTRUE, 30., 180., 
                  kTRUE, 0.01, 20.,0.95);
   
   TLatex *labelElectronTPC = new TLatex(0.5,0.9,"Electrons TPC");
   SetStyleTLatex( labelElectronTPC, 0.05,4);
   labelElectronTPC->Draw();

   if (labelDataSet) labelDataSet->Draw();
   
   
   paddEdxTPC->cd(2);
   paddEdxTPC->cd(2)->SetLogx(1);
   paddEdxTPC->cd(2)->SetLogz(1);
   paddEdxTPC->cd(2)->SetTopMargin(0.01);
   paddEdxTPC->cd(2)->SetRightMargin(0.12);
   histoPdEdxPTPC->GetZaxis()->SetRangeUser(1/nEvt,maximumEPTPC);
   histoPdEdxPTPC->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoPdEdxPTPC,  
                  "", "#it{p}_{e^{+}} (GeV/c)", "dE_{e^{+}-cand} /dx",  "",
                  kTRUE, 30., 180., 
                  kTRUE, 0.01, 20.,0.95);
   
   TLatex *labelPositronTPC = new TLatex(0.5,0.9,"Positrons TPC");
   SetStyleTLatex( labelPositronTPC, 0.05,4);
   labelPositronTPC->Draw();

   paddEdxTPC->cd(3);
   paddEdxTPC->cd(3)->SetLogx(1);
   paddEdxTPC->cd(3)->SetLogz(1);
   paddEdxTPC->cd(3)->SetTopMargin(0.01);
   paddEdxTPC->cd(3)->SetRightMargin(0.12);
   histoENSigmadEdxPTPC->GetZaxis()->SetRangeUser(1/nEvt,maximumEPTPC);
   histoENSigmadEdxPTPC->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoENSigmadEdxPTPC,  
                  "", "#it{p}_{e^{-}} (GeV/c)", "n #sigma_{e^{-}} dE/dx",  "",
                  kTRUE, -10., 10., 
                  kTRUE, 0.01, 20.,0.95);
 
   paddEdxTPC->cd(4);
   paddEdxTPC->cd(4)->SetLogx(1);
   paddEdxTPC->cd(4)->SetLogz(1);
   paddEdxTPC->cd(4)->SetTopMargin(0.01);
   paddEdxTPC->cd(4)->SetRightMargin(0.12);
   histoPNSigmadEdxPTPC->GetZaxis()->SetRangeUser(1/nEvt,maximumEPTPC);
   histoPNSigmadEdxPTPC->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoPNSigmadEdxPTPC,  
                  "", "#it{p}_{e^{+}} (GeV/c)", "n #sigma_{e^{+}} dE/dx",  "",
                  kTRUE, -10., 10., 
                  kTRUE, 0.01, 20.,0.95);
 
   paddEdxTPC->cd(5);
   paddEdxTPC->cd(5)->SetLogx(1);
   paddEdxTPC->cd(5)->SetLogz(1);
   paddEdxTPC->cd(5)->SetTopMargin(0.01);
   paddEdxTPC->cd(5)->SetRightMargin(0.12);
   histoENSigmaPiondEdxPTPC->GetZaxis()->SetRangeUser(1/nEvt,maximumEPTPC);
   histoENSigmaPiondEdxPTPC->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoENSigmaPiondEdxPTPC,  
                  "", "#it{p}_{e^{-}} (GeV/c)", "n #sigma_{#pi^{-}} dE/dx",  "",
                  kTRUE, -10., 10., 
                  kTRUE, 0.01, 20.,0.95);
 
   paddEdxTPC->cd(6);
   paddEdxTPC->cd(6)->SetLogx(1);
   paddEdxTPC->cd(6)->SetLogz(1);
   paddEdxTPC->cd(6)->SetTopMargin(0.01);
   paddEdxTPC->cd(6)->SetRightMargin(0.12);
   histoPNSigmaPiondEdxPTPC->GetZaxis()->SetRangeUser(1/nEvt,maximumEPTPC);
   histoPNSigmaPiondEdxPTPC->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoPNSigmaPiondEdxPTPC,  
                  "", "#it{p}_{e^{+}} (GeV/c)", "n #sigma_{#pi^{+}} dE/dx",  "",
                  kTRUE, -10., 10., 
                  kTRUE, 0.01, 20.,0.95);
   
   canvasdEdxTPC->Update();
   canvasdEdxTPC->SaveAs(Form("fig_ga_Electrons_dEdx_TPC.%s",suffix.Data()));

   TCanvas * canvasTOF = new TCanvas("canvasTOF","",0,0,1000,1000);  // gives the page size
   TPad* padTOF = new TPad("padTOF","",0.0,0.0,1,1,0);   // gives the size of the histo areas 
   padTOF->SetFillColor(0);
   padTOF->GetFrame()->SetFillColor(0);
   padTOF->SetBorderMode(0);
   
   padTOF->Divide(2,2);
   padTOF->Draw();
   padTOF->cd(1);
   padTOF->cd(1)->SetLogx(1);
   padTOF->cd(1)->SetLogz(1);
   padTOF->cd(1)->SetTopMargin(0.01);
   padTOF->cd(1)->SetRightMargin(0.12);
   Double_t maximumEPTOF = 1.2*histoETOFP->GetMaximum();
   histoETOFP->GetZaxis()->SetRangeUser(1/nEvt,maximumEPTOF);
   histoETOFP->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoETOFP,  
                  "", "#it{p}_{e^{-}} (GeV/c)", "t_{measured}-t_{expected} e^{-}",  "",
                  kTRUE, -1000, 2000., 
                  kTRUE, 0.01, 20.,0.95);
   
   TLatex *labelElectronTOF = new TLatex(0.5,0.9,"Electrons TOF");
   SetStyleTLatex( labelElectronTOF, 0.05,4);
   labelElectronTOF->Draw();
   Double_t fracElecTOF = nElectronTOF/nElectronTPC*100;
   TLatex *labelFracElectronTOF = new TLatex(0.5,0.845,Form("%4.2f %%",fracElecTOF ));
   SetStyleTLatex( labelFracElectronTOF, 0.05,4);
   labelFracElectronTOF->Draw();
   

   if (labelDataSet) labelDataSet->Draw();
   
   
   padTOF->cd(2);
   padTOF->cd(2)->SetLogx(1);
   padTOF->cd(2)->SetLogz(1);
   padTOF->cd(2)->SetTopMargin(0.01);
   padTOF->cd(2)->SetRightMargin(0.12);
   histoPTOFP->GetZaxis()->SetRangeUser(1/nEvt,maximumEPTOF);
   histoPTOFP->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoPTOFP,  
                  "", "#it{p}_{e^{+}} (GeV/c)", "t_{measured}-t_{expected} e^{+}",  "",
                  kTRUE, -1000, 2000., 
                  kTRUE, 0.01, 20.,0.95);
   
   TLatex *labelPositronTOF = new TLatex(0.5,0.9,"Positrons TOF");
   SetStyleTLatex( labelPositronTOF, 0.05,4);
   labelPositronTOF->Draw();
   Double_t fracPosiTOF = nPositronTOF/nPositronTPC*100;
   TLatex *labelFracPositronTOF = new TLatex(0.5,0.845,Form("%4.2f %%",fracPosiTOF ));
   SetStyleTLatex( labelFracPositronTOF, 0.05,4);
   labelFracPositronTOF->Draw();
   
   padTOF->cd(3);
   padTOF->cd(3)->SetLogx(1);
   padTOF->cd(3)->SetLogz(1);
   padTOF->cd(3)->SetTopMargin(0.01);
   padTOF->cd(3)->SetRightMargin(0.12);
   histoENSigmaTOFP->GetZaxis()->SetRangeUser(1/nEvt,maximumEPTOF);
   histoENSigmaTOFP->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoENSigmaTOFP,  
                  "", "#it{p}_{e^{-}} (GeV/c)", "n #sigma_{e^{-}} TOF",  "",
                  kTRUE, -10., 10., 
                  kTRUE, 0.01, 20.,0.95);
 
   padTOF->cd(4);
   padTOF->cd(4)->SetLogx(1);
   padTOF->cd(4)->SetLogz(1);
   padTOF->cd(4)->SetTopMargin(0.01);
   padTOF->cd(4)->SetRightMargin(0.12);
   histoPNSigmaTOFP->GetZaxis()->SetRangeUser(1/nEvt,maximumEPTOF);
   histoPNSigmaTOFP->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoPNSigmaTOFP,  
                  "", "#it{p}_{e^{+}} (GeV/c)", "n #sigma_{e^{+}} TOF",  "",
                  kTRUE, -10., 10., 
                  kTRUE, 0.01, 20.,0.95);
 
   canvasTOF->Update();
   canvasTOF->SaveAs(Form("fig_ga_Electrons_TOF.%s",suffix.Data()));
   
   
   TCanvas * canvasITS = new TCanvas("canvasITS","",0,0,1000,1000);  // gives the page size
   TPad* padITS = new TPad("padITS","",0.0,0.0,1,1,0);   // gives the size of the histo areas 
   padITS->SetFillColor(0);
   padITS->GetFrame()->SetFillColor(0);
   padITS->SetBorderMode(0);
   
   padITS->Divide(2,2);
   padITS->Draw();
   padITS->cd(1);
   padITS->cd(1)->SetLogx(1);
   padITS->cd(1)->SetLogz(1);
   padITS->cd(1)->SetTopMargin(0.01);
   padITS->cd(1)->SetRightMargin(0.12);
   Double_t maximumEPITS = 1.2*histoEdEdxPITS->GetMaximum();
   histoEdEdxPITS->GetZaxis()->SetRangeUser(1/nEvt,maximumEPITS);
   histoEdEdxPITS->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoEdEdxPITS,  
                  "", "#it{p}_{e^{-}} (GeV/c)", "dE_{e^{-}-cand} /dx ITS",  "",
                  kTRUE, 0, 180., 
                  kTRUE, 0.01, 20.,0.95);
   
   TLatex *labelElectronITS = new TLatex(0.5,0.9,"Electrons ITS");
   SetStyleTLatex( labelElectronITS, 0.05,4);
   labelElectronITS->Draw();
   Double_t fracElecITS = nElectronITS/nElectronTPC*100;
   TLatex *labelFracElectronITS = new TLatex(0.5,0.845,Form("%4.2f %%",fracElecITS ));
   SetStyleTLatex( labelFracElectronITS, 0.05,4);
   labelFracElectronITS->Draw();
   

   if (labelDataSet) labelDataSet->Draw();
   
   
   padITS->cd(2);
   padITS->cd(2)->SetLogx(1);
   padITS->cd(2)->SetLogz(1);
   padITS->cd(2)->SetTopMargin(0.01);
   padITS->cd(2)->SetRightMargin(0.12);
   
   histoPdEdxPITS->GetZaxis()->SetRangeUser(1/nEvt,maximumEPITS);
   histoPdEdxPITS->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoPdEdxPITS,  
                  "", "#it{p}_{e^{+}} (GeV/c)","dE_{e^{+}-cand} /dx ITS",  "",
                  kTRUE, 0, 180., 
                  kTRUE, 0.01, 20.,0.95);
   
   TLatex *labelPositronITS = new TLatex(0.5,0.9,"Positrons ITS");
   SetStyleTLatex( labelPositronITS, 0.05,4);
   labelPositronITS->Draw();
   Double_t fracPosiITS = nPositronITS/nPositronTPC*100;
   TLatex *labelFracPositronITS = new TLatex(0.5,0.845,Form("%4.2f %%",fracPosiITS ));
   SetStyleTLatex( labelFracPositronITS, 0.05,4);
   labelFracPositronITS->Draw();
   
   padITS->cd(3);
   padITS->cd(3)->SetLogx(1);
   padITS->cd(3)->SetLogz(1);
   padITS->cd(3)->SetTopMargin(0.01);
   padITS->cd(3)->SetRightMargin(0.12);
   histoENSigmadEdxPITS->GetZaxis()->SetRangeUser(1/nEvt,maximumEPITS);
   histoENSigmadEdxPITS->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoENSigmadEdxPITS,  
                  "", "#it{p}_{e^{-}} (GeV/c)", "n #sigma_{e^{-}} ITS",  "",
                  kTRUE, -10., 10., 
                  kTRUE, 0.01, 20.,0.95);
 
   padITS->cd(4);
   padITS->cd(4)->SetLogx(1);
   padITS->cd(4)->SetLogz(1);
   padITS->cd(4)->SetTopMargin(0.01);
   padITS->cd(4)->SetRightMargin(0.12);
   histoPNSigmadEdxPITS->GetZaxis()->SetRangeUser(1/nEvt,maximumEPITS);
   histoPNSigmadEdxPITS->GetXaxis()->SetLabelOffset(-0.009);
   DrawAutoGammaHisto2D( histoPNSigmadEdxPITS,  
                  "", "#it{p}_{e^{+}} (GeV/c)", "n #sigma_{e^{+}} ITS",  "",
                  kTRUE, -10., 10., 
                  kTRUE, 0.01, 20.,0.95);
 
   canvasITS->Update();
   canvasITS->SaveAs(Form("fig_ga_Electrons_ITS.%s",suffix.Data()));
   
   
   
   TCanvas * canvasPhoton = new TCanvas("canvasPhoton","",0,0,1000,1000);  // gives the page size
   TPad* padPhoton = new TPad("padPhoton","",0.0,0.0,1,1,0);   // gives the size of the histo areas 
   padPhoton->SetFillColor(0);
   padPhoton->GetFrame()->SetFillColor(0);
   padPhoton->SetBorderMode(0);
   
   padPhoton->Divide(3,3);
   padPhoton->Draw();
   padPhoton->cd(1);
   padPhoton->cd(1)->SetLogy(1);
   padPhoton->cd(1)->SetLogx(1);
   padPhoton->cd(1)->SetTopMargin(0.01);
   DrawAutoGammaHisto( histoGammaPt,
                       "", "#it{p}_{#gamma,T} (GeV/c)","d#it{N}_{#gamma}/d#it{N}_{evt}",
                       kTRUE, 4., 0.5/nEvt,
                       kFALSE, -10,10,
                       kTRUE, 0.1,30);
   histoGammaPt->SetMarkerSize(0.5);
   histoGammaPt->Draw("ep1");

   TLatex *labelPhotons = new TLatex(0.75,0.9,"Photon");
   SetStyleTLatex( labelPhotons, 0.05,4);
   labelPhotons->Draw();
   
   if (labelDataSet) labelDataSet->Draw();
   
   padPhoton->cd(2);
   padPhoton->cd(2)->SetLogy(1);
   padPhoton->cd(2)->SetTopMargin(0.01);
   
   DrawAutoGammaHisto( histoGammaEta,
                       "", "#eta_{#gamma}","d#it{N}_{#gamma}/d#it{N}_{evt}",
                       kTRUE, 2., 0.5/nEvt,
                       kFALSE, -10,10,
                       kFALSE, 0.1,30);
   histoGammaEta->SetMarkerSize(0.5);
   histoGammaEta->Draw("ep1");

   padPhoton->cd(3);
   padPhoton->cd(3)->SetLogy(1);
   padPhoton->cd(3)->SetTopMargin(0.01);

   DrawAutoGammaHisto( histoGammaPhi,
                       "", "#phi_{#gamma}","d#it{N}_{#gamma}/d#it{N}_{evt}",
                       kTRUE, 2., 0.5/nEvt,
                       kFALSE, -10,10,
                       kFALSE, 0.1,30);
   histoGammaPhi->SetMarkerSize(0.5);
   histoGammaPhi->Draw("ep1");

   padPhoton->cd(4);
   padPhoton->cd(4)->SetLogy(1);
   padPhoton->cd(4)->SetTopMargin(0.01);

   DrawAutoGammaHisto( histoGammaInvMass,
                       "", "M_{#gamma#gamma} (GeV/c^{2})","d#it{N}_{#gamma}/d#it{N}_{evt}",
                       kTRUE, 2., 0.5/nEvt,
                       kFALSE, -10,10,
                       kTRUE, 0.,0.1);
   histoGammaInvMass->SetMarkerSize(0.5);
   histoGammaInvMass->Draw("ep1");
   
   padPhoton->cd(5);
   padPhoton->cd(5)->SetLogy(1);
   padPhoton->cd(5)->SetTopMargin(0.01);

   DrawAutoGammaHisto( histoGammaChi2,
                       "", "#chi^{2}/NDF","d#it{N}_{#gamma}/d#it{N}_{evt}",
                       kTRUE, 2., 0.5/nEvt*1e1,
                       kFALSE, -10,10,
                       kFALSE, 0.,0.1);
   histoGammaChi2->SetMarkerSize(0.5);
   histoGammaChi2->Draw("ep1");
   
   padPhoton->cd(6);
   padPhoton->cd(6)->SetLogy(1);
   padPhoton->cd(6)->SetTopMargin(0.01);

   DrawAutoGammaHisto( histoGammaPsiPair,
                       "", "#psi_{Pair}","d#it{N}_{#gamma}/d#it{N}_{evt}",
                       kTRUE, 2., 0.5/nEvt*1e2,
                       kFALSE, -10,10,
                       kTRUE, 0.,0.5);
   histoGammaPsiPair->SetMarkerSize(0.5);
   histoGammaPsiPair->Draw("ep1");
   
   padPhoton->cd(7);
   padPhoton->cd(7)->SetLogy(1);
   padPhoton->cd(7)->SetTopMargin(0.01);

   DrawAutoGammaHisto( histoGammaCosPoint,
                       "", "cos(#theta_{Point})","d#it{N}_{#gamma}/d#it{N}_{evt}",
                       kTRUE, 2., 0.5/nEvt,
                       kFALSE, -10,10,
                       kFALSE, 0.,0.5);
   histoGammaCosPoint->SetMarkerSize(0.5);
   histoGammaCosPoint->Draw("ep1");
 
   padPhoton->cd(8);
   padPhoton->cd(8)->SetLogz(1);
   padPhoton->cd(8)->SetTopMargin(0.01);
   padPhoton->cd(8)->SetRightMargin(0.12);
   Double_t maximumPhotons = 1.2*histoGammaArmenteros->GetMaximum();
   histoGammaArmenteros->GetZaxis()->SetRangeUser(1/nEvt,maximumPhotons);
   DrawAutoGammaHisto2D( histoGammaArmenteros,  
                  "", "#alpha = (p^{+}_{L}-p^{-}_{L})/(p^{+}_{L}+p^{-}_{L})", "q_{T} (GeV/c)",  "",
                  kFALSE, -10., 10., 
                  kTRUE, -1., 1.,0.95);
 
   canvasPhoton->Update();
   canvasPhoton->SaveAs(Form("fig_ga_Photons.%s",suffix.Data()));
   
}

