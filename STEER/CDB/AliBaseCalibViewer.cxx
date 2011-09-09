///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Base class for the AliTPCCalibViewer and AliTRDCalibViewer               //
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
//#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h> 
#include <TH1F.h>
#include <TMath.h>
#include <THashTable.h>
#include <TObjString.h>
#include <TLinearFitter.h>
#include <TFile.h>
#include <TKey.h>
#include <TGraph.h>
#include <TDirectory.h>
#include <TFriendElement.h>

#include "TTreeStream.h"
#include "AliBaseCalibViewer.h"

ClassImp(AliBaseCalibViewer)

AliBaseCalibViewer::AliBaseCalibViewer()
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
AliBaseCalibViewer::AliBaseCalibViewer(const AliBaseCalibViewer &c)
                  :TObject(c),
                   fTree(0),
                   fFile(0),
                   fListOfObjectsToBeDeleted(0),
                   fTreeMustBeDeleted(0),
                   fAbbreviation(0), 
                   fAppendString(0)
{
  //
  // dummy AliBaseCalibViewer copy constructor
  // not yet working!!!
  //
  fTree = c.fTree;
  fTreeMustBeDeleted = c.fTreeMustBeDeleted;
  fListOfObjectsToBeDeleted = c.fListOfObjectsToBeDeleted;
  fAbbreviation = c.fAbbreviation;
  fAppendString = c.fAppendString;
}

//_____________________________________________________________________________
AliBaseCalibViewer::AliBaseCalibViewer(TTree* tree)
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
AliBaseCalibViewer::AliBaseCalibViewer(const Char_t* fileName, const Char_t* treeName)
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
AliBaseCalibViewer & AliBaseCalibViewer::operator =(const AliBaseCalibViewer & param)
{
   //
   // assignment operator - dummy
   // not yet working!!!
   //
   fTree = param.fTree;
   fTreeMustBeDeleted = param.fTreeMustBeDeleted;
   fListOfObjectsToBeDeleted = param.fListOfObjectsToBeDeleted;
   fAbbreviation = param.fAbbreviation;
   fAppendString = param.fAppendString;
   return (*this);
}

//_____________________________________________________________________________
AliBaseCalibViewer::~AliBaseCalibViewer()
{
   //
   // AliBaseCalibViewer destructor
   // all objects will be deleted, the file will be closed, the pictures will disappear
   //
   if (fTree && fTreeMustBeDeleted) {
      fTree->SetCacheSize(0);
      fTree->Delete();
   }
   if (fFile) {
      fFile->Close();
      fFile = 0;
   }

   for (Int_t i = fListOfObjectsToBeDeleted->GetEntriesFast()-1; i >= 0; i--) {
      delete fListOfObjectsToBeDeleted->At(i);
   }
   delete fListOfObjectsToBeDeleted;
}

//_____________________________________________________________________________
void AliBaseCalibViewer::Delete(Option_t* /*option*/) {
   //
   // Should be called from AliBaseCalibViewerGUI class only.
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

//_____________________________________________________________________________
void AliBaseCalibViewer::FormatHistoLabels(TH1 *histo) const {
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

//_____________________________________________________________________________
TFriendElement* AliBaseCalibViewer::AddReferenceTree(const Char_t* filename, const Char_t* treename, const Char_t* refname){
  //
  // add a reference tree to the current tree
  // by default the treename is 'tree' and the reference treename is 'R'
  //
   TFile *file = new TFile(filename);
   fListOfObjectsToBeDeleted->Add(file);
   TTree * tree = (TTree*)file->Get(treename);
   return AddFriend(tree, refname);
}

//_____________________________________________________________________________
TString* AliBaseCalibViewer::Fit(const Char_t* drawCommand, const Char_t* formula, const Char_t* cuts, 
				 Double_t & chi2, TVectorD &fitParam, TMatrixD &covMatrix){
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
   if (entries == -1) {
     delete fitter;
     return new TString("An ERROR has occured during fitting!");
   }
   Double_t **values = new Double_t*[dim+1] ; 
   
   for (Int_t i = 0; i < dim + 1; i++){
      Int_t centries = 0;
      if (i < dim) centries = fTree->Draw(((TObjString*)formulaTokens->At(i))->GetName(), cutStr.Data(), "goff");
      else  centries = fTree->Draw(drawStr.Data(), cutStr.Data(), "goff");
      
      if (entries != centries) {
	delete fitter;
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
   
   TString *preturnFormula = new TString(Form("( %e+",fitParam[0])), &returnFormula = *preturnFormula;
   
   for (Int_t iparam = 0; iparam < dim; iparam++) {
     returnFormula.Append(Form("%s*(%e)",((TObjString*)formulaTokens->At(iparam))->GetName(),fitParam[iparam+1]));
     if (iparam < dim-1) returnFormula.Append("+");
   }
   returnFormula.Append(" )");
   delete formulaTokens;
   delete fitter;
   for (Int_t i = 0; i < dim + 1; i++) delete [] values[i];
   delete[] values;
   return preturnFormula;
}

//_____________________________________________________________________________
Double_t AliBaseCalibViewer::GetLTM(Int_t n, Double_t *array, Double_t *sigma, Double_t fraction){
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

//_____________________________________________________________________________
Int_t AliBaseCalibViewer::GetBin(Float_t value, Int_t nbins, Double_t binLow, Double_t binUp){
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

//_____________________________________________________________________________
TH1F* AliBaseCalibViewer::SigmaCut(TH1F *histogram, Float_t mean, Float_t sigma, Float_t sigmaMax, 
				   Float_t sigmaStep, Bool_t pm) {
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
         TH1F *distribution = new TH1F("Distribution1", "Distribution f(x, #mu, #sigma)", 1000,-5,5);
         TRandom rand(23);
         for (Int_t i = 0; i <50000;i++) distribution->Fill(rand.Gaus(mean, sigma));
         Float_t *ar = distribution->GetArray();
         
         TCanvas* macro_example_canvas = new TCanvas("cAliBaseCalibViewer", "", 350, 350);
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
   return AliBaseCalibViewer::SigmaCut(nbins, array, mean, sigma, nbins, binLow, binUp, sigmaMax, sigmaStep, pm);
}   
   
//_____________________________________________________________________________
TH1F* AliBaseCalibViewer::SigmaCut(Int_t n, Float_t *array, Float_t mean, Float_t sigma, Int_t nbins, Float_t binLow, Float_t binUp, Float_t sigmaMax, Float_t sigmaStep, Bool_t pm){
   //
   // Creates a histogram Begin_Latex S(t, #mu, #sigma) End_Latex, where you can see, how much of the data are inside sigma-intervals around the mean value
   // The data of the distribution Begin_Latex f(x, #mu, #sigma) End_Latex are given in 'array', 'n' specifies the length of the array
   // 'mean' and 'sigma' are Begin_Latex #mu End_Latex and  Begin_Latex #sigma End_Latex of the distribution in 'array', to be specified by the user
   // 'nbins': number of bins, 'binLow': first bin, 'binUp': last bin
   // sigmaMax: up to which sigma around the mean/median/LTM the histogram is generated (in units of sigma, Begin_Latex t #sigma End_Latex)
   // sigmaStep: the binsize of the generated histogram
   // Here the actual work is done.
   
   if (TMath::Abs(sigma) < 1.e-10) return 0;
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
   
   // fill histogram
   for (Float_t iDelta = 0; iDelta <= sigmaMax * sigma; iDelta += sigmaStep) {
      // integrate array
      Double_t valueP = array[GetBin(mean, nbins, binLow, binUp)];
      Double_t valueM = array[GetBin(mean-binWidth, nbins, binLow, binUp)];
      // add bin of mean value only once to the histogram
      for (Float_t x = binWidth; x <= iDelta; x += binWidth) {
         valueP += (mean + x <= binUp)  ? array[GetBin(mean + x, nbins, binLow, binUp)] : 0;
         valueM += (mean-binWidth - x >= binLow) ? array[GetBin(mean-binWidth - x, nbins, binLow, binUp)] : 0; 
      }

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
      }
   }
   if (!pm) hist->SetMaximum(1.2);
   return hist;
}

//_____________________________________________________________________________
TH1F* AliBaseCalibViewer::SigmaCut(Int_t /*n*/, Double_t */*array*/, Double_t /*mean*/, Double_t /*sigma*/, 
				   Int_t /*nbins*/, Double_t */*xbins*/, Double_t /*sigmaMax*/){
   // 
   // SigmaCut for variable binsize
   // NOT YET IMPLEMENTED !!!
   // 
   printf("SigmaCut with variable binsize, Not yet implemented\n");

   return 0;
}

//_____________________________________________________________________________
Int_t  AliBaseCalibViewer::DrawHisto1D(const Char_t* drawCommand, const Char_t* sector, const Char_t* cuts, 
				       const Char_t *sigmas, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM) const
{
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

//_____________________________________________________________________________
Int_t AliBaseCalibViewer::SigmaCut(const Char_t* drawCommand, const Char_t* sector, const Char_t* cuts, 
				   Float_t sigmaMax, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM, Bool_t pm, 
				   const Char_t *sigmas, Float_t sigmaStep) const {
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
      cutHistoMean = SigmaCut(htemp, mean, sigma, sigmaMax, sigmaStep, pm);
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
      cutHistoMedian = SigmaCut(htemp, median, sigma, sigmaMax, sigmaStep, pm);
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
      cutHistoLTM = SigmaCut(htemp, ltm, ltmRms, sigmaMax, sigmaStep, pm);
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

//_____________________________________________________________________________
Int_t AliBaseCalibViewer::Integrate(const Char_t* drawCommand, const Char_t* sector, const Char_t* cuts, 
				    Float_t /*sigmaMax*/, Bool_t plotMean, Bool_t plotMedian, Bool_t plotLTM, 
				    const Char_t *sigmas, Float_t /*sigmaStep*/) const {
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

//_____________________________________________________________________________
TH1F* AliBaseCalibViewer::Integrate(TH1F *histogram, Float_t mean, Float_t sigma, Float_t sigmaMax, Float_t sigmaStep){
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
         TH1F *distribution = new TH1F("Distribution2", "Distribution f(x, #mu, #sigma)", 1000,-5,5);
         TRandom rand(23);
         for (Int_t i = 0; i <50000;i++) distribution->Fill(rand.Gaus(mean, sigma));
         Float_t *ar = distribution->GetArray();
         
         TCanvas* macro_example_canvas = new TCanvas("macro_example_canvas_Integrate", "", 350, 350);
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
   return Integrate(nbins, array, nbins, binLow, binUp, mean, sigma, sigmaMax, sigmaStep);
}

//_____________________________________________________________________________
TH1F* AliBaseCalibViewer::Integrate(Int_t n, Float_t *array, Int_t nbins, Float_t binLow, Float_t binUp, 
				    Float_t mean, Float_t sigma, Float_t sigmaMax, Float_t sigmaStep){
   // Creates an integrated histogram Begin_Latex S(t, #mu, #sigma) End_Latex, out of the input distribution distribution Begin_Latex f(x, #mu, #sigma) End_Latex, given in "histogram"   
   // "mean" and "sigma" are Begin_Latex #mu End_Latex and  Begin_Latex #sigma End_Latex of the distribution in "histogram", to be specified by the user
   // sigmaMax: up to which sigma around the mean/median/LTM you want to integrate 
   // if "igma == 0" and "sigmaMax == 0" the whole histogram is integrated
   // "sigmaStep": the binsize of the generated histogram, -1 means, that the maximal reasonable stepsize is used
   // Here the actual work is done.
      
   Bool_t givenUnits = kTRUE;
   if (TMath::Abs(sigma) < 1.e-10 && TMath::Abs(sigmaMax) < 1.e-10) givenUnits = kFALSE;
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

//_____________________________________________________________________________
void AliBaseCalibViewer::DrawLines(TH1F *histogram, TVectorF nsigma, TLegend *legend, Int_t color, Bool_t pm) const {
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

//_____________________________________________________________________________
void AliBaseCalibViewer::DrawLines(TGraph *graph, TVectorF nsigma, TLegend *legend, Int_t color, Bool_t pm) const {
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
