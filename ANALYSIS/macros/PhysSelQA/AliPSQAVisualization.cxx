/////////////////////////////////////
// Created by: Kevin McDermott     //
// email: kmcderm3@nd.edu          //
// CERN Summer Student 2012        //
// University of Notre Dame du Lac //
//                                 // 
// Revision: 1.0                   //
// Created on: August 6, 2012      //
/////////////////////////////////////

#include "AliPSQAVisualization.h"

ClassImp(AliPSQAVisualization)

AliPSQAVisualization::AliPSQAVisualization(): // Default constructor 
TObject()
{
  // Initialize some private data members from AliPSQAV
  fInDirectory = "";
  fROOTInput = "";
  fSavePDFs = kFALSE;
  fOutDirectory = "";
  fOutPDFName = ""; 
  fOutEPSName = "";
  fDrawOverPlot = kFALSE;
  fOverPlotTitle = "";
  fMaximum = -1000;
  fMinimum =  1000;
  fScaleAuto = kFALSE;
  fUseColorArray = kFALSE;
  fMultMax = 0;
  fDivMin = 0;
  InitializeColorArray("");
  InitializeSelectedPlots("");
}

//________________________________________________________________________________________________
AliPSQAVisualization::~AliPSQAVisualization(){
  delete[] fCanvas;
  delete[] fDrawPlot;
  delete[] fSelectedPlots;
  delete[] fColors;
  delete fOverLegend;
}

//________________________________________________________________________________________________
void AliPSQAVisualization::InitializeColorArray(const Char_t * listOfColors){ // Function to custom set color array used for plots, not essential
  Color_t colors; //color enums from list
  ifstream input_data_col; //text file object
  input_data_col.open(listOfColors, ios::in );  //open the text file
  Int_t Ncol = 0; //number of color names to be used
  
  while ( input_data_col >> colors ) { // count number of color names
    Ncol++;
  }
  
  fNColors = Ncol; // Set private data member to total color names to be used
  input_data_col.close(); // reset the file

  fColors = new Color_t [fNColors]; // initialize private data member with number of trig names
  input_data_col.open(listOfColors, ios::in );
  for (Int_t icol = 0; icol < fNColors; icol++){
    input_data_col >> fColors[icol];
  }

  input_data_col.close();  
}

//________________________________________________________________________________________________
void AliPSQAVisualization::InitializeSelectedPlots(const Char_t * listOfPlots){
  TString plotnames; //plot names from list
  ifstream input_data_pn; //text file object
	if(!listOfPlots) cout << "No list of plots" << endl;
  input_data_pn.open(listOfPlots, ios::in );  //open the text file
  Int_t Npn = 0; //number of plot names to be used
  
  while ( input_data_pn >> plotnames ) { // count number of plot names
    Npn++;
    }

  fNSelectedPlots = Npn; // Set private data member to total plot names to be used
  input_data_pn.close(); // reset the file
	cout << "Number of selected plots = " << fNSelectedPlots << endl;
  fSelectedPlots = new TString [fNSelectedPlots]; // initialize private data member with number of trig names
  input_data_pn.open(listOfPlots, ios::in );
  for (Int_t iplot = 0; iplot < fNSelectedPlots; iplot++){
    input_data_pn >> fSelectedPlots[iplot];
	  
	  cout << "Plots to process :" << fSelectedPlots[iplot] << endl;
	  
  }

  input_data_pn.close();  
}

//________________________________________________________________________________________________
void AliPSQAVisualization::PostProcessQA(){

	if(!fROOTInput) cout << "No .root file found" << endl;
	
	
	//cout << "Sto morendo qua? 1" << endl;
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
   // do not use scientific notation for run number
  TGaxis::SetMaxDigits(6);
  //cout << "Sto morendo qua? 2" << endl;
  ConvertTGraphErrorsToTH1Ds();  // Convert all plots to TH1D's from specified plot list
  //cout << "Sto morendo qua? 3" << endl;
  fCanvas = new TCanvas[fNSelectedPlots]; // Allocate memory for the canvases to draw PostProcessing
  
   
  if ((fDrawOverPlot) || (fSaveOverPlotPDF) || (fSaveOverPlotEPS)){
    fOverLegend = new TLegend(0.55,0.625,.75,.85);
    if (fScaleAuto){
      ScaleMinAndMax();
    }
  }

  if ((fSavePDFs) || (fSaveOverPlotPDF) || (fSaveOverPlotEPS)){
    FileStat_t dummy_filestat;
    if (gSystem->GetPathInfo(fOutDirectory.Data(), dummy_filestat) == 1) { // cache directory does not exist
      MakeDir(fOutDirectory);
    }
  }

  for (Int_t iplot = 0; iplot < fNSelectedPlots; iplot++){ // Loop over the number of plots to be made pretty and drawn
    DrawSelected(iplot); // Draw the Selected plots, decide later whether to have them open/save them
  
    if (fSavePDFs){ // Write the canvas to the PDF file if kTRUE
      SavePDFs(iplot);
    }

    if (!fDrawSelected){ // Close the canvas if you do not want to see the plot
      fCanvas[iplot].Close();
    }

    if ( (fDrawOverPlot) || (fSaveOverPlotPDF) || (fSaveOverPlotEPS) ) { // Set points and draw overplot if any are true
      if (fUseColorArray){
	DrawOverPlotCA(iplot);
      }
      else{
	DrawOverPlotNoCA(iplot);
      }
    }
    else if( (!fDrawOverPlot) && ( (fSaveOverPlotPDF) || (fSaveOverPlotEPS) ) ) { // Set the points to save to and/or pdf/eps but then close it outside the loop
      if (fUseColorArray){
	DrawOverPlotCA(iplot);
      }
      else{
	DrawOverPlotNoCA(iplot);
      }
    }	
  }

  if (fSaveOverPlotPDF){ // Save the overplot to pdf
    SaveOverPlotPDF();
  }
  if (fSaveOverPlotEPS){ // Save the overplot to eps
    SaveOverPlotEPS();
  }
  if (!fDrawOverPlot){ // Close if you do not want to see it
    fOverCanvas.Close();
  }
}

//________________________________________________________________________________________________
void AliPSQAVisualization::ConvertTGraphErrorsToTH1Ds(){
  
	fRootFile = TFile::Open(Form("%s",fROOTInput.Data()));
	
	
  fDrawPlot = new TH1D[fNSelectedPlots]; // Array of TH1D's to be converted
  
  for (Int_t iplot = 0; iplot < fNSelectedPlots; iplot++){
    // Variables from the old TGraphErrors needed to draw TH1Fs
    // Uncomment if seg fault returns TGraph::GetX (this=0x0), as last printout will say which graph is not correct
    cout << fInDirectory.Data() << "/" << fROOTInput.Data() << ": " << fSelectedPlots[iplot].Data() << endl;

	  
	  
	  cout << "Problem is here? " << fSelectedPlots[iplot].Data() << endl;
    TGraphErrors * TGEplot  = (TGraphErrors*) fRootFile->Get(fSelectedPlots[iplot].Data());
	  if(!TGEplot) {
		  cout << "Missing plot : " <<   fSelectedPlots[iplot].Data() << endl;
		  continue;
	  }
	  
    Double_t     * TGEx     = TGEplot->GetX(); // needed to properly index the x-axis (run numbers)
    Double_t     * TGEy     = TGEplot->GetY(); // array of double_t's, qaval
    Int_t        TGEnpoints = TGEplot->GetN(); // number of points used in the TGE
    Double_t      * TGEey   = TGEplot->GetEY(); // array of double_t's corresponding to the errors in the TGE
    TString      TGEname    = TGEplot->GetName();  
    TString      TGEtitle   = TGEplot->GetTitle(); 
    // Parameters for the to be drawn TH1D
    TString     TH1Dname  = TGEname.Data();
    TH1Dname  += "_TH1D";
    TString     TH1Dtitle = TGEtitle.Data();
    TH1Dtitle += " TH1D";
    // See Below for the strange numbering of bins

    Int_t       TH1Dxnbin = TGEnpoints;
    Int_t       TH1Dxlow  = 1;
    Int_t       TH1Dxhigh = TH1Dxnbin + 1;
    // Declare the new plot to be drawn as a TH1D
    
    fDrawPlot[iplot].SetBins(TH1Dxnbin, TH1Dxlow, TH1Dxhigh); // start with bin label 1, up through bin number +1 (TH1D bin numbering has strange conventions, as zero bin is the underflow bin)

    fDrawPlot[iplot].SetNameTitle(TH1Dname.Data(),TH1Dtitle.Data());
    fDrawPlot[iplot].GetXaxis()->SetTitle("Run Number");
    TString TH1Dytitle = SetDrawSelectedYTitle(TGEtitle);
    fDrawPlot[iplot].SetYTitle(TH1Dytitle.Data());

    for (Int_t ibin = 1; ibin < TH1Dxnbin+1; ibin++){ // everything shifted off by one because tgraph starts at zero, th1f starts bins at 1
      fDrawPlot[iplot].SetBinContent(fDrawPlot[iplot].FindBin(ibin),TGEy[ibin-1]); // y[i-1] is the bin in plot
      fDrawPlot[iplot].SetBinError(fDrawPlot[iplot].FindBin(ibin),TGEey[ibin-1]);
      
      Int_t TH1Dx = Int_t(TGEx[ibin-1]); // x[i-1] = actual bin label in plot
      TString TH1DxBinLabel = Form("%i",TH1Dx);
      fDrawPlot[iplot].GetXaxis()->SetBinLabel(fDrawPlot[iplot].FindBin(ibin), TH1DxBinLabel.Data());
      fDrawPlot[iplot].GetXaxis()->SetTitleSize(.05); // default = 4%, made it 5%
    }
    //
    // adding some drawing options
    //
    fDrawPlot[iplot].SetMarkerStyle(20);
    fDrawPlot[iplot].SetMarkerSize(1.4);
    fDrawPlot[iplot].SetMarkerColor(kBlue);
    fDrawPlot[iplot].SetLineColor(kBlue);
    fDrawPlot[iplot].GetYaxis()->SetRangeUser(-0.001,1.05);
    fDrawPlot[iplot].GetXaxis()->LabelsOption("V");
  }
}

//________________________________________________________________________________________________
TString AliPSQAVisualization::SetDrawSelectedYTitle(TString TGEtitle){
  
  Ssiz_t start  = 0;  
  Ssiz_t first_o = TGEtitle.First('o'); // o in over, lowercase o that appears first is the over
  Ssiz_t first_f = TGEtitle.First('f');

  TString label[2];
  label[0] = TString(TGEtitle( start , first_o - start - 1)); // Shift one over from the 'h'
  label[1] = TString(TGEtitle( first_o + 4, first_f - first_o - 4 ));  //Shift five over from the 'O' 

  // Make this substitution to make legend smaller
  
  if (label[1].Contains("Trigger class")){
    label[1].ReplaceAll("Trigger class","All");
  }

  TString yTH1Dtitle = label[0];
  yTH1Dtitle += " / ";
  yTH1Dtitle += label[1];

  return yTH1Dtitle;
}

//________________________________________________________________________________________________
void AliPSQAVisualization::ScaleMinAndMax(){ // Get maximum and minimum value of traghs to set bounds for overplot
  Double_t tmp_max; // max values from tgraphs
  Double_t tmp_min; // min values from tgraphs

  for (Int_t iplot = 0; iplot < fNSelectedPlots; iplot++){ // Loop over the number of plots to be made 
    // Set maximum
    tmp_max = fDrawPlot[iplot].GetMaximum();
    if (tmp_max > fMaximum){
       fMaximum = tmp_max;
    }
    // Set minimum
    tmp_min = fDrawPlot[iplot].GetMinimum();
    if (tmp_max < fMinimum){
      fMinimum = tmp_min;
    }
  }

  // Set the min to be lower and max to be higher so that all points can be seen
  fMaximum *= fMultMax;
  fMinimum /= fDivMin;
}

//________________________________________________________________________________________________
void AliPSQAVisualization::MakeDir(TString dir){ // have to pass a private data member, kinda ugly, but otherwise use sentinel values or copy this function three times... taken from AliPSQA
  TString mkDir = "mkdir -p ";
  mkDir += dir.Data();
  gSystem->Exec(mkDir); // mkdir for the cache/output
}

//________________________________________________________________________________________________
void AliPSQAVisualization::DrawSelected(Int_t iplot){
   fCanvas[iplot].cd(); // write to this canvas first
   fCanvas[iplot].SetBottomMargin(0.15); // default 0.15
   fCanvas[iplot].SetRightMargin(0.02);
   fDrawPlot[iplot].GetXaxis()->SetTitleOffset(1.6);
   fDrawPlot[iplot].Draw("E P X0"); // draw the new plot to canvas
}

//________________________________________________________________________________________________
void AliPSQAVisualization::DrawOverPlotCA(Int_t iplot){
  if (iplot == 0){
    fOverCanvas.cd();
    fDrawPlot[iplot].SetMarkerColor(fColors[iplot]); // iplot == 0 for Color_t == white...
    fDrawPlot[iplot].SetMarkerStyle(Style_t(iplot+20)); // 20 = big circle, afterwards, useful styles, otherwise too small
    fDrawPlot[iplot].SetLineColor(fColors[iplot]);
    fDrawPlot[iplot].SetTitle(fOverPlotTitle);
    fOverLegend->AddEntry(&fDrawPlot[iplot],fDrawPlot[iplot].GetYaxis()->GetTitle(),"lep");
    fDrawPlot[iplot].GetYaxis()->SetTitle("");
    fDrawPlot[iplot].GetYaxis()->SetRangeUser(fMinimum, fMaximum);
    fDrawPlot[iplot].Draw("E P X0");
  }
  else if (iplot == (fNSelectedPlots - 1) ){ 
    fOverCanvas.cd();
    fDrawPlot[iplot].SetMarkerColor(fColors[iplot]);
    fDrawPlot[iplot].SetMarkerStyle(Style_t(iplot+20)); // 20 = big circle, afterwards, useful styles, otherwise too small
    fDrawPlot[iplot].SetLineColor(fColors[iplot]);
    fDrawPlot[iplot].SetTitle(fOverPlotTitle);
    fOverLegend->AddEntry(&fDrawPlot[iplot],fDrawPlot[iplot].GetYaxis()->GetTitle(),"lep");
    fDrawPlot[iplot].GetYaxis()->SetTitle("");
    fDrawPlot[iplot].Draw("E P X0 SAME");
    //    fOverLegend->Draw();
  }
  else{
    fOverCanvas.cd();
    fDrawPlot[iplot].SetMarkerColor(fColors[iplot]);
    fDrawPlot[iplot].SetMarkerStyle(Style_t(iplot+20)); // 20 = big circle, afterwards, useful styles, otherwise too small
    fDrawPlot[iplot].SetLineColor(fColors[iplot]);
    fDrawPlot[iplot].SetTitle(fOverPlotTitle);
    fOverLegend->AddEntry(&fDrawPlot[iplot],fDrawPlot[iplot].GetYaxis()->GetTitle(),"lep");
    fDrawPlot[iplot].GetYaxis()->SetTitle("");
    fDrawPlot[iplot].Draw("E P X0 SAME");
  }
}

//________________________________________________________________________________________________
void AliPSQAVisualization::DrawOverPlotNoCA(Int_t iplot){
  if (iplot == 0){
    fOverCanvas.cd();
    fDrawPlot[iplot].SetMarkerColor(Color_t(iplot+1)); // iplot == 0 for Color_t == white...
    fDrawPlot[iplot].SetMarkerStyle(Style_t(iplot+20)); // 20 = big circle, afterwards, useful styles, otherwise too small
    fDrawPlot[iplot].SetLineColor(Color_t(iplot+1));
    fDrawPlot[iplot].SetTitle(fOverPlotTitle);
    fOverLegend->AddEntry(&fDrawPlot[iplot],fDrawPlot[iplot].GetYaxis()->GetTitle(),"lep");
    fDrawPlot[iplot].GetYaxis()->SetTitle("");
    fDrawPlot[iplot].GetYaxis()->SetRangeUser(fMinimum, fMaximum);
    fDrawPlot[iplot].Draw("E P X0");
  }
  else if (iplot == (fNSelectedPlots - 1) ){ 
    fOverCanvas.cd();
    fDrawPlot[iplot].SetMarkerColor(Color_t(iplot+1));
    fDrawPlot[iplot].SetMarkerStyle(Style_t(iplot+20)); // 20 = big circle, afterwards, useful styles, otherwise too small
    fDrawPlot[iplot].SetLineColor(Color_t(iplot+1));
    fDrawPlot[iplot].SetTitle(fOverPlotTitle);
    fOverLegend->AddEntry(&fDrawPlot[iplot],fDrawPlot[iplot].GetYaxis()->GetTitle(),"lep");
    fDrawPlot[iplot].GetYaxis()->SetTitle("");
    fDrawPlot[iplot].Draw("E P X0 SAME");
    fOverLegend->Draw();
  }
  else{
    fOverCanvas.cd();
    fDrawPlot[iplot].SetMarkerColor(Color_t(iplot+1));
    fDrawPlot[iplot].SetMarkerStyle(Style_t(iplot+20)); // 20 = big circle, afterwards, useful styles, otherwise too small
    fDrawPlot[iplot].SetLineColor(Color_t(iplot+1));
    fDrawPlot[iplot].SetTitle(fOverPlotTitle);
    fOverLegend->AddEntry(&fDrawPlot[iplot],fDrawPlot[iplot].GetYaxis()->GetTitle(),"lep");
    fDrawPlot[iplot].GetYaxis()->SetTitle("");
    fDrawPlot[iplot].Draw("P SAME");
  }
}

//________________________________________________________________________________________________
void AliPSQAVisualization::SavePDFs(Int_t iplot){
  if ( iplot == 0 ){ // Open the PDF file if to be saved and save the first plot
    fCanvas[iplot].Print(Form("%s/%s(",fOutDirectory.Data(),fOutPDFName.Data()),"pdf");
	  
  }
  else if ( (iplot == (fNSelectedPlots - 1)) && (!fSaveOverPlotPDF) ){ // Close the PDF and save the last plot
    fCanvas[iplot].Print(Form("%s/%s)",fOutDirectory.Data(),fOutPDFName.Data()),"pdf");
	    fCanvas[iplot].SaveAs(Form("%s/%s_%d.pdf",fOutDirectory.Data(),fOutPDFName.Data(),iplot));
  }
  else{ // Save the PDF with the inbetween plots
    fCanvas[iplot].Print(Form("%s/%s",fOutDirectory.Data(),fOutPDFName.Data()),"pdf");
  }
}

//________________________________________________________________________________________________
void AliPSQAVisualization::SaveOverPlotPDF(){
  fOverCanvas.cd();
  if (fSavePDFs){
    fOverCanvas.Print(Form("%s/%s)",fOutDirectory.Data(),fOutPDFName.Data()),"pdf");
	
  }
  else{
    fOverCanvas.Print(Form("%s/%s",fOutDirectory.Data(),fOutPDFName.Data()),"pdf");
  }
}

//________________________________________________________________________________________________
void AliPSQAVisualization::SaveOverPlotEPS(){
  fOverCanvas.cd();
  fOverCanvas.Print(Form("%s/%s",fOutDirectory.Data(),fOutEPSName.Data()),"eps");
}

