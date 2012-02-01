#include "AliMassFitControl.h"

// Function for fitting back ground away from peak
Float_t quad(Double_t *x, Double_t *par){
  if (x[0] > par[3] && x[0] < par[4]){
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0]; //a+bx+cx**2
}

TH1F* PtMassAna2(TH2F *PtMass, Int_t mode,Int_t ihist, const Int_t NControl, TObjArray *ControlArray,Char_t* fulllabel){

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1111);
  //TString *tLabel = new TString(fulllabel); //not reqd

  //Arguments TH2F *PtMass - 2D histogram of some variable versus mass
  //          Int_t mode - mode for switching things depending on particle
  //          const Int_t NControl - number of bins to project into which must
  //          be the same as the number of controller objects (see next arg).
  //          TObjArray *ContolArray - holds the controllers objects, 1 per bin

  //const Int_t NDefault=27; // Used in default array (change as required)
  // Trap user errors
  if (ControlArray->GetEntries() != NControl){
    cout << "PtMassAna2: Problem: Specified number of bins and number of bins in supplied TObjArray do not match" << endl;
    return 0;
  }
  //Set up mode mapping to particle and names for histograms
  if(mode == 0){ 
    const Char_t* part = "K0";
    const Char_t* partName = "K^{0}_{s}";
  }
  else if (mode == 1 ){
    const Char_t* part = "LAM";
    const Char_t* partName = "#Lambda";
  }
  else if (mode == 2 ){//since lambda and anti-Lambda have same limits etc.
    const Char_t* part = "LAM"; 
    const Char_t* partName = "#bar{#Lambda}";
  }
  else if (mode ==3) {// This is for combined Lambda + anti-Lambda
    const Char_t* part = "LAM";
    const Char_t* partName = "#Lambda + #bar{#Lambda}";
  }
  else if (mode ==4) {
    const Char_t* part = "XI";
    const Char_t* partName = "#Xi";
  }
  else{
    cout << "Mode not recognized, defaulting to K0" << endl;
    const Char_t* part = "K0";
  }
  cout << "Debug info. Particle: " << part << endl;

  Float_t xxlo, xxhi; //Exclusion limits for fitting
  // FUTURE - these can also be part of AliMassFitControl to allow different regions
  // in different bins
  Float_t NSigmaEx=4.0;  // Number of sigma to exclude from bkgd fit
  Float_t NSigmaInt=3.5; // Number of sigma to integrate histo over

  // Fitting function for initial combined peak + background fit
  TF1 *gausQuad = new TF1("gausQuad","(([1]/([2]*pow(6.2831,0.5)))*[6]*exp(-0.5*pow((x-[0])/[2],2)))+ [3] + [4]*x + [5]*x*x",0.3,1.8);
  //>>>>>>>>> NB Need to set proper range when it is used <<<<<<<<<<<<<<<<<<<<<
  //>>>>>>>>> AND use 'R' in fit option                   <<<<<<<<<<<<<<<<<<<<<
  gausQuad->SetLineWidth(1);
  gausQuad->SetLineColor(kRed);

  // Fitting function for background 
  TF1* quad = new TF1("quad",quad,0.3,1.8,5);
  quad->SetLineWidth(2);
  quad->SetLineColor(kBlue);
  quad->SetLineStyle(2);

  // Function for background under peak - to be integrated and drawn
  TF1* qback = new TF1("qback","[0]+[1]*x+[2]*x*x",0.3,1.8);
  qback->SetLineWidth(2);
  qback->SetLineColor(kRed);
  qback->SetLineStyle(2);
  qback->SetFillStyle(3013);
  qback->SetFillColor(kRed);

  // Set ranges and limits according to particle
  Float_t xmin,xmax;
  Float_t defMean,defArea,defWidth,defa0,defa1,defa2,defBW; //defaults
  Float_t defMnLo, defMnHi;
  if (part=="K0"){
    // FUTURE: xmin, xmax define the range considered in fit so could go into
    // AliMassFitControl
    //xmin=0.418;
    //xmax=0.578;
    //    xxlo=0.46;xxhi=0.53;
    //    gausQuad->SetParameters(0.4977,2000,0.003,1,1,0,0.001*RBFact);
    //gausQuad->SetParameters(0.4977,2000,0.003,1,1,0,0.001);
    defMean=0.4977; defArea=1000; defWidth=0.01;
    defa0=2; defa1=2; defa2=0;
    //gausQuad->SetParLimits(0,0.483,0.523); // Constrain mean
    defMnLo = 0.483; defMnHi=0.523;
  }
  if (part=="LAM"){
    // FUTURE: See above
    //xmin=1.085;
    //xmax=1.17;
    //  xxlo=1.10;xxhi=1.135;
    //    gausQuad->SetParameters(1.115,2000,0.0028,10,100,-100,0.001*RBFact);
    //  gausQuad->SetParameters(1.115,2000,0.0028,10,100,-100,0.001);
    defMean=1.115; defArea=1000; defWidth=0.08;
    defa0=100; defa1=-10; defa2=0;
    //    gausQuad->SetParLimits(0,1.113,1.123); // Constrain mean
    defMnLo=1.113; defMnHi=1.123;
  }
  if (part=="XI"){
    defMean=1.32171; defArea=500; defWidth=0.08;
    defa0=100; defa1=-10; defa2=0;
    defMnLo=1.301; defMnHi=1.341;
  }
  //gausQuad->SetRange(xmin,xmax);
  gausQuad->SetParLimits(1,0.0,1E7); // Constrain area 0 to 10 million
  gausQuad->SetParLimits(2,0.0001,0.04); // Constrain width 0.1 MeV to 40 MeV
  gausQuad->SetParNames("Mean","Area1","Sigma1","p0","p1","p2","BinWid"); 

  // Deciding how to divide up the canvas
  c1->Clear();
  c1->SetWindowSize(1024,768);
  Int_t NDraw = NControl;
  if(NDraw==16){c1->Divide(4,4);}
  elseif(NDraw<=36 && NDraw>30){ c1->Divide(6,6);}
  elseif(NDraw<=30 && NDraw>25) { c1->Divide(6,5);}
  elseif(NDraw<=25 && NDraw>20){ c1->Divide(5,5);}
  elseif(NDraw<=20 && NDraw>16){  c1->Divide(5,4);}
  elseif(NDraw<=15 && NDraw>12){ c1->Divide(5,3);}
  elseif(NDraw<=12 && NDraw>8){ c1->Divide(4,3);}
  elseif(NDraw<=8 && NDraw>6){ c1->Divide(4,2);}
  elseif(NDraw<=6 && NDraw>4){c1->Divide(3,2);}
  elseif(NDraw<=4 && NDraw>2){c1->Divide(2,2);}
  elseif(NDraw==2){c1->Divide(2,1);}
  elseif(NDraw==1){/*do nothing*/;}
  else{c1->Divide(7,6);}

  //
  // Project out the histograms from 2D into 1D 'slices'
  //
  Char_t id[20], title[80];
  //cout << "NControl: " << NControl << endl;
  TH1D* hMassSlice[NControl];
  //Arrays to store various quantities which can later be histogrammed.
  //  const Int_t NBinsArrays = NBins+1-NFirst; // array of 1D projection histograms is counted from 0 but other arrays are used for histograms contents and therefore N+1 are need because element zero is left empty
  const Int_t NBinsArrays = NControl+1;
  Stat_t sigBkgd[NBinsArrays], errSigBkgd[NBinsArrays]; // Signal+background and error on this
  Stat_t bkgd[NBinsArrays], errBkgd[NBinsArrays]; // Background and error on this
  Stat_t signal[NBinsArrays], errSignal[NBinsArrays]; // Net signal and error

  Stat_t chi2PerDOF1[NBinsArrays], errChi2PerDOF1[NBinsArrays]; // Chi-squared per defree of freedom from the initial fit
  Stat_t chi2PerDOF2[NBinsArrays], errChi2PerDOF2[NBinsArrays]; // Chi-squared per defree of freedom from the 2nd (background) fit

  // REDO won't need this?
  // Zero the signal - avoids problems with plotting
  for(Int_t j=0;j<NBinsArrays;j++){
    signal[j]=0.0;
    errSignal[j]=1.0;
  } // <<< Not actually in use a the moment <<<

  Stat_t mean[NBinsArrays], errMean[NBinsArrays]; // Mean and error
  Stat_t sigma[NBinsArrays], errSigma[NBinsArrays]; // Gaussian width and error
  Stat_t area[NBinsArrays], errArea[NBinsArrays]; // Peak area from fit and error
  Stat_t integLow[NBinsArrays], integUp[NBinsArrays]; // Lower/upper limits for integration

  // Float_t DefaultBinPtEdges[NDefault]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0};
  Float_t BinPtEdges[NBinsArrays];
  
  //  Float_t BinPtEdges[NBins+1]={0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0};

  Int_t BinLo, BinHi; //For pt bins
  // **** Main loop over the AliMassFitControllers ****
 
  TIter controlIter(ControlArray);
  AliMassFitControl *controller;
  while (controller=(AliMassFitControl*)controlIter.Next()) { 
    if(ihist == 0)controller->CalcBinLimits(20); //This BinsPerGeV argument should be calculated from the data
    if(ihist == 1)controller->CalcBinLimits(50); //This BinsPerGeV argument should be calculated from the data
    if(ihist == 2)controller->CalcBinLimits(50); //This BinsPerGeV argument should be calculated from the data
    if(ihist == 3)controller->CalcBinLimits(1); //This BinsPerGeV argument should be calculated from the data
    if(ihist == 4)controller->CalcBinLimits(1); //This BinsPerGeV argument should be calculated from the data
    if(ihist == 6)controller->CalcBinLimits(20000); //This BinsPerGeV argument should be calculated from the data
    if(ihist == 5)controller->CalcBinLimits(100); //This BinsPerGeV argument should be calculated from the data
    //Had to introduce this Nint fn otherwise got inconsistencies after type implicit conversion
    //    BinLo=TMath::Nint(1.+BinPtEdges[N]*20.); //cout << "BinLo: " << BinLo << ", " << 1+BinPtEdges[N]*20. << endl;
    //    BinHi=TMath::Nint(BinPtEdges[N+1]*20.); //cout << "BinHi: " << BinHi << ", " << BinPtEdges[N+1]*20. << endl;
    BinLo=controller->BinLower();
    BinHi=controller->BinUpper();
    Int_t N=ControlArray->IndexOf(controller) - 1;
    sprintf(id,"Mass%d",N);
    cout << "Mass histo:" << N << " Firstbin: " << BinLo << " Last bin:" << BinHi << " " << controller->PtLower() << "-" << controller->PtUpper() << " GeV" << endl; 
    //cout << "About to create mass projection " << N << endl;
    if(ihist == 0) hMassSlice[N] = PtMass->ProjectionX(id,BinLo,BinHi);
    if(ihist != 0) hMassSlice[N] = PtMass->ProjectionY(id,BinLo,BinHi);
    //cout << "Mass projection " << N << " created." << endl;
   sprintf(title,"%s Mass, %.2f < p_{t} < %.2f",partName,controller->PtLower(),controller->PtUpper());
    hMassSlice[N]->SetTitle(title);
    hMassSlice[N]->SetXTitle("GeV/c^{2}");
    // cout << "Mass projection " << N << " titles set." << endl;
    hMassSlice[N]->Rebin(controller->RebinFactor());
    // } // End of creating the projections
    Int_t massBinLo,massBinHi; //For mass histogram bins
    Int_t NArray;
  // Do all the fitting
    //REDO this is all one loop (iteration over TObjArray) now
  //for(Int_t N=0;N<NBins-NFirst;N++){
    NArray=N+1; // Arrays numbered from 1 not 0..
    c1->cd(N+1); //Canvas subdivisions numbered from 1 and NOT 0

    // Set range - can be different each time
    xmin=controller->MinMass();
    xmax=controller->MaxMass();
    gausQuad->SetRange(xmin,xmax);
    // Everytime need to set some sensible parameters for both fits
    gausQuad->SetParameters(defMean,defArea,defWidth,defa0,defa1,defa2,0.1);
    quad->SetParameters(defa0,defa1,defa2,0.49,0.51);// Last 2 parameters are
    // the range to exclude but these are fixed later by FixParameter(3,.. etc.
    gausQuad->FixParameter(6,hMassSlice[N]->GetBinWidth(1));
    // Release Linear parameter
    gausQuad->ReleaseParameter(4);
    quad->ReleaseParameter(1);
    // Release Quadratic parameter
    gausQuad->ReleaseParameter(5);
    quad->ReleaseParameter(2);
    gausQuad->SetParLimits(0,defMnLo,defMnHi);
    gausQuad->SetParLimits(1,0.0,1E7); // Constrain area 0 to 10 million
    gausQuad->SetParLimits(2,0.0001,0.04); // Constrain width 0.1 MeV to 40 MeV
//     gausQuad->SetParLimits(4,0,0); // Releases linear parameter
//     quad->SetParLimits(1,0,0); // Releases linear parameter
//     gausQuad->SetParLimits(5,0,0); // Releases quadratic parameter
//     quad->SetParLimits(2,0,0); // Releases quadratic parameter
    if(controller->FixedLin()){
      quad->FixParameter(1,0.0);
      gausQuad->FixParameter(4,0.0);
      cout << "Info: linear part fixed to zero." << endl;
      gausQuad->SetParLimits(3,0.0,1000); // Ensure constant is +ve
      // NB documentation indicates that limits only used when option B
      // is specified. Also how are they can they be freed later?
    }
    if(controller->FixedQuad()){ 
      quad->FixParameter(2,0.0);
      gausQuad->FixParameter(5,0.0);
      cout << "Info: quadratic part fixed to zero." << endl;
    }
    hMassSlice[N]->Fit(gausQuad,"LR");
    hMassSlice[N]->Draw("E,SAME");
    gausQuad->DrawCopy("SAME");
    hMassSlice[N]->Clone("hMySlice");
    //Fit the background:
    xxlo=gausQuad->GetParameter(0)-NSigmaEx*gausQuad->GetParameter(2);
    xxhi=gausQuad->GetParameter(0)+NSigmaEx*gausQuad->GetParameter(2);
    // Limits have to be adjusted to match bins
    massBinLo=hMassSlice[N]->FindBin(xxlo);
    xxlo = hMassSlice[N]->GetBinLowEdge(massBinLo);
    massBinHi=hMassSlice[N]->FindBin(xxhi)+1;
    xxhi = hMassSlice[N]->GetBinLowEdge(massBinHi);
    quad->SetParameters(gausQuad->GetParameter(3),gausQuad->GetParameter(4),gausQuad->GetParameter(5),xxlo,xxhi);
    quad->FixParameter(3,xxlo);
    quad->FixParameter(4,xxhi);
    quad->SetRange(xmin,xmax);
    hMassSlice[N]->Fit(quad,"LRN+");
    quad->SetFillColor(30);
    quad->SetFillStyle(3013);
    quad->SetRange(xmin,xxlo);
    quad->DrawCopy("SAME");
    quad->SetRange(xxhi,xmax);
    quad->DrawCopy("SAME");

    // Draw the background - use a new function since other one not defined everywhere
    xxlo=gausQuad->GetParameter(0)-NSigmaInt*gausQuad->GetParameter(2);
    xxhi=gausQuad->GetParameter(0)+NSigmaInt*gausQuad->GetParameter(2);
    massBinLo=hMassSlice[N]->FindBin(xxlo);
    xxlo = hMassSlice[N]->GetBinLowEdge(massBinLo);
    massBinHi=hMassSlice[N]->FindBin(xxhi)+1;
    xxhi = hMassSlice[N]->GetBinLowEdge(massBinHi);
    qback->SetRange(xxlo,xxhi);
    qback->SetParameter(0,quad->GetParameter(0));
    qback->SetParameter(1,quad->GetParameter(1));
    qback->SetParameter(2,quad->GetParameter(2));
    qback->DrawCopy("SAME");
    if(N==28){
TCanvas *myCan = new TCanvas();
TPad *mypad = new TPad();
myCan->cd();
mypad->Draw();
hMassSlice[N]->Draw();
quad->SetRange(xmin,xxlo);
quad->DrawCopy("SAME");
quad->SetRange(xxhi,xmax);
quad->DrawCopy("SAME");
qback->DrawCopy("SAME");

}
   

    //Integrate the signal+background (i.e. original histo)
    sigBkgd[NArray]=hMassSlice[N]->Integral(massBinLo,massBinHi);
    errSigBkgd[NArray]=TMath::Sqrt(sigBkgd[NArray]);
    //Integrate background - better to do this analytically ?
    bkgd[NArray]=qback->Integral(xxlo,xxhi)/gausQuad->GetParameter(6);//Divide by bin width
    errBkgd[NArray]=TMath::Sqrt(bkgd[NArray]); //Get errors from fit parameters?

    // Fill arrays for diagnostic histograms
    mean[NArray]=gausQuad->GetParameter(0);
    errMean[NArray]=gausQuad->GetParError(0);
    area[NArray]=gausQuad->GetParameter(1);
    errArea[NArray]=gausQuad->GetParError(1);
    sigma[NArray]=gausQuad->GetParameter(2);
    errSigma[NArray]=gausQuad->GetParError(2);
    chi2PerDOF1[NArray]=gausQuad->GetChisquare()/(Double_t)gausQuad->GetNDF();
    errChi2PerDOF1[NArray]= 0.0; // Don't know how to calc. error for this?
    chi2PerDOF2[NArray]=quad->GetChisquare()/(Double_t)quad->GetNDF();
    errChi2PerDOF2[NArray]= 0.0; // Don't know how to calc. error for this?

    BinPtEdges[N]=controller->PtLower();
    //    BinPtEdges(N+1)=controller->PtUpper();
  }
    BinPtEdges[NBinsArrays-1]=((AliMassFitControl*)ControlArray->Last())->PtUpper();
//     for (Int_t jj=0;jj<NBinsArrays;jj++){
//       cout << "BinPtEdges " << jj << " = " << BinPtEdges[jj] << endl;
//       cout << "Mean " << jj << " = " << mean[jj] << endl;
//       cout << "Sigma " << jj << " = " << sigma[jj] << endl;
//       cout << "Signal + bkgd " << jj << " = " << sigBkgd[jj] << endl;
//       cout << "Background " << jj << " = " << bkgd[jj] << endl;
//     }
  // Yields in each bin returned in a histogram
  TH1F* hYields= new TH1F("Yield","Raw Particle Yield",NBinsArrays-1,BinPtEdges);
  //hYields->Reset("ICE");
  hYields->SetXTitle("p_{t} / (GeV/c)");
  //hYields->SetContent(sigBkgd);
  //hYields->SetError(errSigBkgd);
  //c1->cd(4);
  //hYields->Draw();

  // Fill the diagnostic histograms: clone, set title, contents, errors.
  TH1F* hMeans = hYields->Clone("Means");
  hMeans->SetTitle("Mean from Gaussian Fit");
  hMeans->SetContent(mean);
  hMeans->SetError(errMean);
  //hMeans->Sumw2();
  TH1F* hSigmas = hYields->Clone("Sigmas");
  hSigmas->SetTitle("Sigma from Gaussian Fit");
  hSigmas->SetContent(sigma);
  hSigmas->SetError(errSigma);
  TH1F* hAreas = hYields->Clone("Areas");
  hAreas->SetTitle("Area from Gaussian Fit");
  hAreas->SetContent(area);
  hAreas->SetError(errArea);
  TH1F* hSigBkgd = hYields->Clone("SigBkgd");
  hSigBkgd->SetTitle("Signal + Background (Histo. Integral)");
  hSigBkgd->SetContent(sigBkgd);
  hSigBkgd->SetError(errSigBkgd);
  TH1F* hBkgd = hYields->Clone("Bkgd");
  hBkgd->SetTitle("Background (Fn. integral)");
  hBkgd->SetContent(bkgd);
  hBkgd->SetError(errBkgd);

  hYields->Sumw2();
  hYields->Add(hSigBkgd,hBkgd,1.0,-1.0);

  TH1F* hSoverB = hYields->Clone("SoverB");
  hSoverB->SetTitle("Signal to Background ratio");
  hSoverB->Sumw2();
  hSoverB->Divide(hBkgd);

  TH1F* hSoverA = hYields->Clone("SoverA"); // Signal over area
  hSoverA->SetTitle("Ratio of Signal to Area from Fit");
  hSoverA->Sumw2();
  hSoverA->Divide(hAreas);

  TH1F* hChi2PerDOF1 = hYields->Clone("Chi2PerDOF1");
  hChi2PerDOF1->SetTitle("Chi-squared per d.o.f. - initial fit");
  hChi2PerDOF1->SetContent(chi2PerDOF1);
  hChi2PerDOF1->SetError(errChi2PerDOF1);

  TH1F* hChi2PerDOF2 = hYields->Clone("Chi2PerDOF2");
  hChi2PerDOF2->SetTitle("Chi-squared per d.o.f. - background fit");
  hChi2PerDOF2->SetContent(chi2PerDOF2);
  hChi2PerDOF2->SetError(errChi2PerDOF2);

  // Draw the diagnostic histograms on their own canvas
  TCanvas* cDiag = new TCanvas("Diag","Diagnostics",600,600);
  cDiag->Divide(2,3);

  cDiag->cd(1);
  Float_t bookmass;
  if(part=="LAM"){
    bookmass=1.11563;
    hMeans->SetMinimum(bookmass-0.01);
    hMeans->SetMaximum(bookmass+0.01);
  } else if (part=="K0"){
    bookmass=0.4976;
    hMeans->SetMinimum(bookmass-0.01);
    hMeans->SetMaximum(bookmass+0.01);
  } else if (part=="XI") {
    hMeans->SetMaximum(1.34);
    hMeans->SetMinimum(1.30);
    bookmass=1.32171;
  }
  Float_t maxPt = hMeans->GetBinLowEdge(hMeans->GetNbinsX()+1);
  TLine PDGmass;
  PDGmass.SetLineStyle(2);
  hMeans->Draw();
  PDGmass.DrawLine(0,bookmass,maxPt,bookmass);

  cDiag->cd(2);
  if (part=="K0"){
    hSigmas->SetMaximum(0.06);
  }else{
    hSigmas->SetMaximum(0.006);}
  hSigmas->SetMinimum(0.0);
  hSigmas->Draw();

  cDiag->cd(3);
  hSoverB->SetMinimum(0.0);
  hSoverB->Draw();

  cDiag->cd(4);
  hSoverA->SetMinimum(0.5);
  hSoverA->SetMaximum(1.5);
  hSoverA->Draw();

  cDiag->cd(5);
  hChi2PerDOF1->SetMinimum(0);
  hChi2PerDOF1->SetMaximum(6);
  hChi2PerDOF1->Draw();

  cDiag->cd(6);
  hChi2PerDOF2->SetMinimum(0);
  hChi2PerDOF2->SetMaximum(6);
  hChi2PerDOF2->Draw();

  Char_t fileNameBase[80];
  sprintf(fileNameBase,"Diagnostics%s_%s",part,fulllabel);
  Char_t fileNamePng[80];
  sprintf(fileNamePng,"%s.png",fileNameBase);
  Char_t fileNameEps[80];
  sprintf(fileNameEps,"%s.eps",fileNameBase);
  Char_t fileNamePdf[80];
  sprintf(fileNamePdf,"%s.pdf",fileNameBase);
  
  cDiag->SaveAs(fileNamePng);
  cDiag->SaveAs(fileNameEps);
  cDiag->SaveAs(fileNamePdf);
 
  // Scale result by the bin size to get dN/dpT and by bin centre to get 1/pt...
  Float_t y, ey;
  for(Int_t j=0;j<=hYields->GetNbinsX();j++){
    y = hYields->GetBinContent(j)/hYields->GetBinWidth(j);
    ey = hYields->GetBinError(j)/hYields->GetBinWidth(j);
    hYields->SetBinContent(j,y);
    hYields->SetBinError(j,ey);
//    y = hYields->GetBinContent(j)/hYields->GetBinCenter(j);
//    ey = hYields->GetBinError(j)/hYields->GetBinCenter(j);
    hYields->SetBinContent(j,y);
    hYields->SetBinError(j,ey);
  }
  return hYields;

}
