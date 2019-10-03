
#include <TAxis.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLine.h>
#include <TMatrixDSym.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>
#include <cstring>
#include <THnSparse.h>
#include <TLatex.h>
#include <TMinuit.h>
#include <TRandom3.h>

#include <TFractionFitter.h>
#include <AliMCLogLFitter.h>
#include <mystyle.h>


TH2D * gAliFitDiagramData; 
TH2D * gAliFitDiagramDalitz;
TH2D * gAliFitDiagramConversion;
TH2D * gAliFitDiagramD;
TH2D * gAliFitDiagramB;
TH2D * gAliFitDiagramHadrons;
TH2D * gAliFitDiagramLambdaProtons;

TH2D * gAliTempSave;
TH1D * gAliTempBMCguess;

using namespace std;

double avgChisq;

const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

void LowerStatistics2D(TH2D * distribution, double factor)
{
  // Can be used to lower statistics in a sample while keeping the Poisson statistics approximately intact
  if(factor > 1)factor=1.0/factor;
  if(factor<=0)factor=1.0;
  TRandom3 *  rand = new TRandom3(3);
  double * rd = new double[3];
  
  int removed=0, maxremove=int(distribution->Integral()*(1.0-factor));
  int binX=0,binY=0,binNr=0;
  int max = int(distribution->GetMaximum());
  
  while(removed<maxremove)
  {
    rand->RndmArray(3, rd);
    binX=int(rd[0]*distribution->GetNbinsX()+1);
    binY=int(rd[1]*distribution->GetNbinsY()+1);
    binNr=distribution->GetBin(binX,binY);
    
    if(distribution->GetBinContent(binNr)>int(rd[2]*max+1.0))
    {
      distribution->Fill(distribution->GetXaxis()->GetBinCenter(binX),distribution->GetYaxis()->GetBinCenter(binY),-1);
      removed++;
    }
  }  
}



TH2D * MirrorDiagram(TH2D * input)
{
  // Mirrors Histogram along middle of y-axis
  TH2D * mirroredDiagram=(TH2D*)input->Clone();
  mirroredDiagram->SetName(Form("%s-mirrored",input->GetName()));
  mirroredDiagram->Scale(0);double xVal, yVal;
  //int nBins=input->GetNbinsY();
  for(int i=1;i<input->GetNbinsX();i++)
    for(int j=1;j<input->GetNbinsY();j++)
    {
      xVal=input->GetXaxis()->GetBinCenter(i);
      yVal=input->GetYaxis()->GetBinCenter(j);
      mirroredDiagram->Fill(xVal,-yVal,input->GetBinContent(i,j));
    }
    return mirroredDiagram;
}

void ShiftToLeft(TH1D * in, int n)
{
  for(int i=1;i<=in->GetNbinsX()-n;i++)
  {
    in->SetBinContent(i,in->GetBinContent(i+n));
  }
}

void ShiftToRight(TH1D * in, int n)
{
  for(int i=in->GetNbinsX();i>=1+n;i--)
  {
    in->SetBinContent(i,in->GetBinContent(i-n));
  }
}

void CleanUpDiagram(TH1D * input, double min, double max) // Sets everything outside range to 0
{
  int nBins=input->GetNbinsX();
  for(int i=1;i<=nBins;i++)
    if(input->GetBinCenter(i)<min || input->GetBinCenter(i)>max) input->SetBinContent(i,0.);
}

void ReturnChiSq(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) 
{
  // Calculates Chi2 for minimization
  if(gin || npar || iflag==0){};
  //int nBinsY=gAliFitDiagramData->GetYaxis()->GetNbins();
  // 0: Conversion, 1: Dalitz, 2: D, 3: B, 4: pT-Bin
  Double_t Chisq=0, up, low;
  //par[1]=par[0];
  Double_t D,B,Data,Dalitz,Conversion, LP;
  Double_t Fittot,Datatot;
  int lowbin,upbin;
  lowbin=gAliFitDiagramData->GetYaxis()->FindBin(-0.06);
  upbin=gAliFitDiagramData->GetYaxis()->FindBin(0.06); // 0.06 for pp
  int diff= upbin-lowbin;
  Datatot=Fittot=0;
  for(int i=lowbin;i<upbin;i++) 
  {
    LP=gAliFitDiagramLambdaProtons->GetBinContent(int(par[4]),i);
    D=gAliFitDiagramD->GetBinContent(int(par[4]+0.5),i);
    B=gAliFitDiagramB->GetBinContent(int(par[4]+0.5),i);
    Dalitz=gAliFitDiagramDalitz->GetBinContent(int(par[4]+0.5),i);
    Conversion=gAliFitDiagramConversion->GetBinContent(int(par[4]+0.5),i);
    Fittot+=Conversion*par[0]+D*par[2]+B*par[3]+Dalitz*par[1];
    Data=gAliFitDiagramData->GetBinContent(int(par[4]+0.5),i);
    Datatot+=Data;
    up=(Data-Conversion*par[0]-Dalitz*par[1]-D*par[2]-B*par[3]-LP*par[5]*0)*(Data-Conversion*par[0]-Dalitz*par[1]-D*par[2]-B*par[3]-LP*par[5]*0);
    low=(TMath::Max(Data,0.0)+Conversion*par[0]*par[0]+Dalitz*par[1]*par[1]+D*par[2]*par[2]+B*par[3]*par[3]+LP*par[5]*par[5]*0);
    //low=(TMath::Max(Data,0.0)+Conversion*par[0]*par[0]+Dalitz*par[1]*par[1]+D*par[2]*par[2]+B*par[3]*par[3]);
    if(low<7.0) low=7.0;
    //      if(low<1.001) low=1.0;
    //if(gAliFitDiagramData->GetYaxis()->GetBinCenter(i)**2 > 9)
    Chisq+=up/low/double(diff+1-4);
  }
  f=Chisq;
  avgChisq=f;
}

double FitChiSq(int bin, double * sourceElectrons)
{
  int temp;
  double savB, savD, savDa, savCo,savBG, avChisqbest, tempDBL;
  TH1D * Bslice=gAliFitDiagramB->ProjectionY("bSliceTemp",bin,bin);
  TH1D * Dslice=gAliFitDiagramD->ProjectionY("dSliceTemp",bin,bin);
  TH1D * Convslice=gAliFitDiagramConversion->ProjectionY("coSlicTempe",bin,bin);
  TH1D * Daslice=gAliFitDiagramDalitz->ProjectionY("daSliceTemp",bin,bin);
  TH1D * Dataslice=gAliFitDiagramData->ProjectionY("dataSliceTemp",bin,bin);
  
  TMinuit * minimizingObject = new TMinuit(6); // new for Terhi
  minimizingObject->SetPrintLevel(-1);
  minimizingObject->SetFCN(ReturnChiSq); // new for Terhi
  minimizingObject->mnparm(4,"pT Bin", bin, 0.1, 0.1, 26.0, temp);
  minimizingObject->FixParameter(4);
  avChisqbest=1e10;
  for(int run=0; run<20;run++){
    //minimizingObject->mnrset(1);
          minimizingObject->mnparm(0,"Cv Factor", 0.3, 0.3, 0.002, 1.0,temp);  // 1.0 ; 0.6 ; ; 5
	  minimizingObject->mnparm(1,"Da Factor", 0.3, 0.3, 0.002, 1.0,temp);  // 1.0 ; 0.6 ; ; 5
          minimizingObject->mnparm(2,"D Factor", 0.9, (0.01+double(run)/60.0)/2, 0.001, 1.0,temp); // 0.9 ; 0.2 ; 0.7
	  minimizingObject->mnparm(3,"B Factor", 0.2, 0.001, 0.001, 1.0,temp);//*/ 
    
    //Pb-Pb
    /*minimizingObject->mnparm(0,"Cv Factor", 0.6*10.0, 0.3*10.0, 0.0*10.0, 20.0,temp);  // 1.0 ; 0.6 ; ; 5
    minimizingObject->mnparm(1,"Da Factor", 0.6*10.0, 0.3*10.0, 0.0*10.0, 20.0,temp);  // 1.0 ; 0.6 ; ; 5
    minimizingObject->mnparm(2,"D Factor", 0.15/2.0*10.0, (0.01+double(run)/60.0)/2*10.0, 0.05/2.0*10.0, 50,temp); // 0.9 ; 0.2 ; 0.7
    minimizingObject->mnparm(3,"B Factor", 0.04*10.0, 0.001*10.0, 0.001*10.0, 10,temp);//*/
    minimizingObject->mnparm(5,"LP Factor", 0.00, 0.001, 0.0000, 0.04,temp);    
    minimizingObject->FixParameter(5);
    
    for(int seeknr=0; seeknr < run/2.0 ; seeknr++)
      minimizingObject->mnseek();
    minimizingObject->Migrad();
    minimizingObject->mnimpr();
    if(avgChisq<avChisqbest)
    {
      minimizingObject->GetParameter(0, savCo, tempDBL);
      minimizingObject->GetParameter(1, savDa, tempDBL); 
      minimizingObject->GetParameter(2, savD, tempDBL); 
      minimizingObject->GetParameter(3, savB, tempDBL);
      minimizingObject->GetParameter(5, savBG, tempDBL);
      avChisqbest=avgChisq;
    }
    //cout << "try: " << run+1 << " chisq: " << avgChisq << " best: " << avChisqbest << endl;
  }
  minimizingObject->mnparm(0,"Cv Factor", savCo, 0.5, 0.0, 2.5*100.0, temp);  // 1.0 ; 0.6 ; ; 5
  minimizingObject->mnparm(1,"Da Factor", savDa, 0.5, 0.0, 2.55*100.0, temp);  // 1.0 ; 0.6 ; ; 5
  minimizingObject->mnparm(2,"D Factor", savD, 0.5, 0.0, 50, temp); // 0.9 ; 0.2 ; 0.7
  minimizingObject->mnparm(3,"B Factor", savB, 0.1, 0.00, 10, temp);  //  was .2 ; .15*/ 
  minimizingObject->mnparm(5,"BG", savBG, 0.0001, 0.00, 0.25*100.0, temp);  //  was .2 ; .15*/ 
  sourceElectrons[0]=savCo;sourceElectrons[1]=savDa;sourceElectrons[2]=savD;sourceElectrons[3]=savB;sourceElectrons[4]=savBG;
  
  //minimizingObject->mnprin(1, 12345);
  
  delete minimizingObject;
  delete Bslice;
  delete Dslice;
  delete Convslice;
  delete Daslice;
  delete Dataslice;
  return avChisqbest;
}

double FitLikelyhoodOld(int bin, double * sourceElectrons)
{
  double temp;
  TH1D * Bslice=gAliFitDiagramB->ProjectionY("bSliceTemp",bin,bin);
  TH1D * Dslice=gAliFitDiagramD->ProjectionY("dSliceTemp",bin,bin);
  TH1D * Convslice=gAliFitDiagramConversion->ProjectionY("coSliceTemp",bin,bin);
  TH1D * Daslice=gAliFitDiagramDalitz->ProjectionY("daSliceTemp",bin,bin);
  TH1D * Dataslice=gAliFitDiagramData->ProjectionY("dataSliceTemp",bin,bin);
  int lowBin=Dataslice->FindBin(-0.15);
  int upBin=Dataslice->FindBin(0.15);
  double Bnr=Bslice->Integral(lowBin, upBin);
  double Dnr=Dslice->Integral(lowBin, upBin);
  double Convnr=Convslice->Integral(lowBin, upBin);
  double Danr=Daslice->Integral(lowBin, upBin);
  double totNr=Dataslice->Integral(lowBin, upBin);
  TFractionFitter * fitter;
  TObjArray * mc = new TObjArray(1);
  mc->Clear();
  mc->Add(Convslice);
  mc->Add(Daslice);
  mc->Add(Dslice);
  mc->Add(Bslice);
  fitter = new TFractionFitter(Dataslice,mc);

  fitter->SetRangeX(lowBin, upBin);
  fitter->Constrain(0,0.02,1);
  fitter->Constrain(1,0.02,1);
  fitter->Constrain(2,0.02,1);
  fitter->Constrain(3,0.02,1);
  fitter->Fit();  
  fitter->GetResult(0, sourceElectrons[0], temp);
  fitter->GetResult(1, sourceElectrons[1], temp);  
  fitter->GetResult(2, sourceElectrons[2], temp);  
  fitter->GetResult(3, sourceElectrons[3], temp);  
  cout << "from new Fit: " << sourceElectrons[0] << " " <<sourceElectrons[1] << " " <<sourceElectrons[2] << " " <<sourceElectrons[3] << " Sum: "
  << sourceElectrons[0]+sourceElectrons[1]+sourceElectrons[2]+sourceElectrons[3] << endl;
  sourceElectrons[0]=totNr/Convnr*sourceElectrons[0];
  sourceElectrons[1]=totNr/Danr*sourceElectrons[1];
  sourceElectrons[2]=totNr/Dnr*sourceElectrons[2];
  sourceElectrons[3]=totNr/Bnr*sourceElectrons[3];
  cout << "from new Fit: " << sourceElectrons[0] << " " <<sourceElectrons[1] << " " <<sourceElectrons[2] << " " <<sourceElectrons[3] << " Ratio to all: " <<
       (sourceElectrons[0]*Convnr+sourceElectrons[1]*Danr+sourceElectrons[2]*Dnr+sourceElectrons[3]*Bnr)/totNr << endl;
  delete Bslice;
  delete Dslice;
  delete Convslice;
  delete Daslice;
  delete Dataslice;
  return fitter->GetChisquare();
}

double FitLikelyhood(int bin, double * sourceElectrons)
{
  int nFunctions=4;
  
  TH1D ** MCDiagrams = new TH1D*[nFunctions];  // Conversion, Dalitz, Charm, Beauty
  MCDiagrams[3]=gAliFitDiagramB->ProjectionY("bSliceTemp",bin,bin);
  MCDiagrams[2]=gAliFitDiagramD->ProjectionY("dSliceTemp",bin,bin);
  MCDiagrams[0]=gAliFitDiagramConversion->ProjectionY("coSliceTemp",bin,bin);
  //if(bin==11)
    //ShiftToRight(MCDiagrams[0],2);   // Shift! Take out again!
  MCDiagrams[1]=gAliFitDiagramDalitz->ProjectionY("daSliceTemp",bin,bin);
  TH1D * Dataslice=gAliFitDiagramData->ProjectionY("dataSliceTemp",bin,bin);

  AliMCLogLFitter * Fitter = new AliMCLogLFitter(nFunctions ,(TH1**)MCDiagrams,(TH1*)Dataslice);
  Fitter->IfNoMCParticleIsProbablyA(3);
  //Fitter->SetUseChiSq();
  double pT = gAliFitDiagramData->GetXaxis()->GetBinCenter(bin);
  //double range = 0.05+0.15/pT;
  double range = 0.05+3./20./pT;
  
  if(range==0) cout << "strange pT range" << endl;
  //Fitter->SetFitRange(-range, range); // 0.15/x+0.05
  //Fitter->SetFitRange(-0.06, 0.06);
  Fitter->SetFitRange(-0.1, 0.1);
  //Fitter->SetParameter(0, 0.001, 0.25, 0.6);
  //Fitter->SetParameter(1, 0.001, 0.25, 0.6);
  //Fitter->FitVeryCarefully();
  Fitter->Fit();
   
  
  for(int i=0;i<nFunctions;i++) sourceElectrons[i] = Fitter->ReturnParameter(i);
  TH1D * expectedConversion = (TH1D*)Fitter->ReturnExpectedDiagram(0); expectedConversion->Scale(sourceElectrons[0]);
  TH1D * expectedDalitz = (TH1D*)Fitter->ReturnExpectedDiagram(1); expectedDalitz->Scale(sourceElectrons[1]);
  TH1D * expectedCharm = (TH1D*)Fitter->ReturnExpectedDiagram(2); expectedCharm->Scale(sourceElectrons[2]);
  TH1D * expectedBeauty = (TH1D*)Fitter->ReturnExpectedDiagram(3); expectedBeauty->Scale(sourceElectrons[3]);
  //cout << "scales: " << sourceElectrons[0] << " " << sourceElectrons[1] << " " << sourceElectrons[2] << " " << sourceElectrons[3] << endl;
  gAliTempBMCguess = (TH1D*)expectedConversion->Clone(); gAliTempBMCguess->SetName("BMCExpected");
  gAliTempBMCguess->Add(expectedDalitz);
  gAliTempBMCguess->Add(expectedCharm);
  gAliTempBMCguess->Add(expectedBeauty);
  
  delete Fitter;
  
  AliMCLogLFitter * Fitter2 = new AliMCLogLFitter(nFunctions ,(TH1**)MCDiagrams,(TH1*)Dataslice);
  Fitter2->IfNoMCParticleIsProbablyA(3); 
  Fitter2->SetFitRange(-0.06, 0.06);
  //gAliTempSave=Fitter2->ScanVariable(2, 0.21, 0.35);
  //Fitter2->Fit(); // take out again
  //gAliTempSave=Fitter2->ScanVariable2d(2, 0.1, 0.3, 3, 0.01, 0.05);
  //Fitter->ScanAll();
  for(int i=0;i<nFunctions;i++) delete MCDiagrams[i];
  delete MCDiagrams;
  delete Dataslice;
  delete Fitter2;
  return 0;
}


Double_t HistFunction(double x, int nr, int binx=1) // not currently used
{
  int biny = gAliFitDiagramDalitz->GetYaxis()->FindBin(x);
  bool smooth=!false;
  double output=0;
  //if(x*x<1) return 0;
  if(smooth && x*x >-50){
    if(nr==1)
      output = TMath::Sqrt(TMath::Sqrt(gAliFitDiagramDalitz->GetBinContent(binx, biny)*gAliFitDiagramDalitz->GetBinContent(binx, biny-1))*TMath::Sqrt(gAliFitDiagramDalitz->GetBinContent(binx, biny)*gAliFitDiagramDalitz->GetBinContent(binx, biny+1)));
    if(nr==0)
      output = TMath::Sqrt(TMath::Sqrt(gAliFitDiagramConversion->GetBinContent(binx, biny)*gAliFitDiagramConversion->GetBinContent(binx, biny-1))*TMath::Sqrt(gAliFitDiagramConversion->GetBinContent(binx, biny)*gAliFitDiagramConversion->GetBinContent(binx, biny+1)));
    if(nr==2)
      output = TMath::Sqrt(TMath::Sqrt(gAliFitDiagramD->GetBinContent(binx, biny)*gAliFitDiagramD->GetBinContent(binx, biny-1))*TMath::Sqrt(gAliFitDiagramD->GetBinContent(binx, biny)*gAliFitDiagramD->GetBinContent(binx, biny+1)));
    if(nr==3)
      output = TMath::Sqrt(TMath::Sqrt(gAliFitDiagramB->GetBinContent(binx, biny)*gAliFitDiagramB->GetBinContent(binx, biny-1))*TMath::Sqrt(gAliFitDiagramB->GetBinContent(binx, biny)*gAliFitDiagramB->GetBinContent(binx, biny+1)));
  }
  else
  {
    if(nr==1)
      output = gAliFitDiagramDalitz->GetBinContent(binx, biny);
    if(nr==0)
      output = gAliFitDiagramConversion->GetBinContent(binx, biny);
    if(nr==2)
      output = gAliFitDiagramD->GetBinContent(binx, biny);
    if(nr==3)
      output = gAliFitDiagramB->GetBinContent(binx, biny);
  }
  return output;
}


Double_t TotalHistFunction(Double_t * x, Double_t *par)
{
  // combines 0-Conversion, 1-Dalitz, 2-c and 3-b diagrams
  return par[0]*HistFunction(x[0], 0, int(par[4]))+par[1]*HistFunction(x[0], 1, int(par[4]))+par[2]*HistFunction(x[0], 2, int(par[4]))+par[3]*HistFunction(x[0], 3, int(par[4]));
}


TH2D * Remap(TH2D * inDiag, char * name)
{
  int nBinsInX, nBinsInY;
  nBinsInX=inDiag->GetNbinsX(); nBinsInY=inDiag->GetNbinsY(); 
  
  
  Double_t  binLimits[1000];//={-0.2,-0.16,-0.15,-0.14,-0.13,-0.12,-0.11,-0.1,-0.09,-0.085,-0.08,-0.075,-0.07,-0.065,-0.06,-0.055,-0.05,-0.045,-0.04,-0.036,-0.032,-0.03,-0.028,-0.026,-0.024,-0.022,-0.02,-0.018,-0.016,-0.014,-0.013,-0.012,-0.011,-0.01,-0.009,-0.008,-0.007,-0.006,-0.005,-0.004,-0.003,-0.002,-0.001,0.0,
  //0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.016,0.018,0.02,0.024,0.028,0.032,0.036,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.2}; // for absolute  
  for(int i=0;i<=90;i++)
    binLimits[i]=-0.2+double(i)*0.002;
  for(int i=1;i<=80;i++)
    binLimits[i+90]=-0.02+double(i)*0.0005;
  for(int i=1;i<=91;i++)
    binLimits[i+170]=0.02+double(i)*0.002;
  
    Double_t ptbinningX[19] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 12., 16., 20.};
    // Double_t ptbinningX[35] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5,2.75, 3.,3.5, 4., 5., 6., 8., 12., 16., 20.}; 
    
    int tot; int binX, binY; tot=0;
    for(int i=0;binLimits[i]<0.2;i++)tot=i+1;
    //cout << "total number: " << tot << endl;
    
    TH2D * eleData2dRebinned = new TH2D(name,name,18, ptbinningX,tot, binLimits); 
    int nBinsin=nBinsInX*nBinsInY;
    //cout << "Total bin Nr: " << nBinsin << endl;
    for(int j=1; j<nBinsin;j++)
    {
      int temp;
      inDiag->GetBinXYZ(j,binX,binY,temp);
      eleData2dRebinned->Fill(inDiag->GetXaxis()->GetBinCenter(binX),inDiag->GetYaxis()->GetBinCenter(binY), inDiag->GetBinContent(j));
      //eleData2dRebinned->Fill(inDiag->GetXaxis()->GetBinCenter(binX),-inDiag->GetYaxis()->GetBinCenter(binY), inDiag->GetBinContent(j));
    }
    return eleData2dRebinned;
    
}


TH2D * Remap2(TH2D * inDiag, char * name)
{
  int nBinsInX, nBinsInY;
  nBinsInX=inDiag->GetNbinsX(); nBinsInY=inDiag->GetNbinsY(); 
  TH2D * outDiag = (TH2D*)inDiag->Clone();
  outDiag->SetName(name);
  outDiag->Scale(0);
  outDiag->RebinY(20);
  Double_t  binLimits[110];
  for(int j=0;j<21;j++)
    if(j>10)
      binLimits[j]=0.05*TMath::Power(2,(TMath::Abs(j-10.0)-1));
    else if(j<10)
      binLimits[j]=-0.05*TMath::Power(2,(TMath::Abs(j-10.0)-1));
    else binLimits[j]=0.0;
    TH1D * eleData1dRebinned = new TH1D("ds","sd",20, binLimits); 
  
  for(int j=1; j<nBinsInX;j++)
  {
    eleData1dRebinned->Scale(0.0);
    for(int i=1;i<nBinsInY;i++)
    {
      //cout << "Bin Center: " << inDiag->GetYaxis()->GetBinCenter(i) << " Content: " << inDiag->GetBinContent(j,i) << endl;
      eleData1dRebinned->Fill(inDiag->GetYaxis()->GetBinCenter(i),inDiag->GetBinContent(j,i));
    }
    cout << "Test2 " << j << endl;
    for(int i=1; i<21; i++)
      outDiag->SetBinContent(j,i+nBinsInY/2/20-10,eleData1dRebinned->GetBinContent(i));
  }
  delete eleData1dRebinned;
  return outDiag;
}

TH2D * SwitchBinning(TH2D * inDiag, char * name)
{
  Double_t ptbinningX[19] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 12., 16., 20.};
  //Double_t ptbinningX[25] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5, 3., 3.5, 4., 5., 6., 8., 12., 16., 20.}; 
  
  int nxIn= inDiag->GetXaxis()->GetNbins();
  int nyIn= inDiag->GetYaxis()->GetNbins();
  double low =  inDiag->GetYaxis()->GetBinLowEdge(1);
  double up =  inDiag->GetYaxis()->GetBinUpEdge(nyIn);
  double increment = (up - low)/double(nyIn);
  Double_t * IPbinningY= new Double_t[nyIn+1];
  cout << "Switch binning from : " << low << " to " << up <<  " with " << nyIn << "increments of " << increment << endl;
  for(int i=0; i<=nyIn; i++){IPbinningY[i]=low + double(i)*increment;}
  TH2D * outDiag = new TH2D(name, name, 18, ptbinningX, nyIn, IPbinningY);
  for(int i=1;i<2000;i++) // y
    for(int j=1; j<=nxIn; j++) // x
      outDiag->Fill(inDiag->GetXaxis()->GetBinCenter(j),inDiag->GetYaxis()->GetBinCenter(i), inDiag->GetBinContent(j,i));
   delete IPbinningY;
   return outDiag;  
}


void FitBCdc()
{
  mystyle();
  
  
  int rebinY=10; // 5
  int lowbin=7, upbin=14; // 9 (7),14
  bool remap=false;
  bool binWidthCorrection=true; // not anymore done for the file output, only for diagram
  bool errorBarEstimation=true;
  bool useTOF = true;
  int errormax=5;
  TString * OutputFileName = new TString("SpectraFitPbPb-05HalfRAAlastHadron.root");
  
  bool useLikelyhood=true;
  
  double events;
  
  bool useLargerCentralityForDalitz = true;
  
  // Make these different from 1 to test what would happen with more (or less) statistics
  double increaseStatistics=1;
  double increaseStatisticsMC=1;
  //
  
  int DrawParticularBin = 8;
  
  
  TH1D ** distOneBin = new TH1D*[upbin-lowbin+1];
  for(int i=0;i<upbin-lowbin+1;i++) distOneBin[i] = new TH1D(Form("BErrDist%d",i), Form("BErrDist%d",i), 500, 0.0, 1.0);
  TH1D * distOneBinD = new TH1D("oneBinD", "oneBinD", 500, 0.0, 1.0);
  TH1D * distOneBinDB = new TH1D("oneBinD", "oneBinD", 500, 0.0, 5.0);
  TH1D * distOneBinConv = new TH1D("oneBinco", "oneBinco", 200, 0.0, 2.5);
  TH1D * distOneBinDalitz = new TH1D("oneBinda", "oneBinda", 200, 0.0, 2.5);
  
  TH1D * distOneBinBextra = new TH1D("distOneBinBextra", "distOneBinBextra", 500, 0.0, 1.0);
  TH1D * distOneBinDextra = new TH1D("distOneBinDextra", "distOneBinDextra", 500, 0.0, 1.0);
  TH1D * distOneBinDBextra = new TH1D("distOneBinDBextra", "distOneBinDBextra", 500, 0.0, 5.0);
  
  TH1D * BExpectationOneBin;
  bool ExtraErrorTestSuccessful=false;
  bool ExtraErrorTestNeeded=false;
  
  TLine ** parPos = new TLine*[upbin-lowbin+1];
  
  //TList * ppLambdaListNew;  
  //AliHFEcollection * ppLambdaCollectionNew;
  
  TH2D * eleData2d2;
  
  TH2D * eleData2d;
  TH2D * eleD2d;
  TH2D * eleB2d;
  TH2D * eleConversion2d;
  TH2D * eleDalitz2d;
  
  TH2D * eleHadrons2d;
  
  TH1D * whiteLine = new TH1D;
    whiteLine->SetLineColor(kWhite);
  
    
  TFile * inFileTreeData;
  if(useTOF)
    inFileTreeData = new TFile("local/TreesJan20/TOF/outputTreePbPbDataEta-05pos020.root");
  else
    inFileTreeData = new TFile("local/TreesJan20/noTOF/outputTreePbPbDataEta-05pos020.root");
  double EventsTrees = ((TH1D*)inFileTreeData->Get("EventNr"))->Integral();
  double EventsTreesBin0 = EventsTrees;//((TH1D*)inFileTreeData->Get("Bin0Numbers"))->GetBinContent(2);
  cout << "Nr of events in outputTreeData: " << EventsTrees << endl;
  cout << "Nr of events in outputTreeData - with bin0: " << EventsTreesBin0 << endl;
  events=EventsTreesBin0;
  TH3D * DCAwithCut3d = (TH3D*) inFileTreeData->Get("TrackDCApTnoCut");
  //DCAwithCut3d->RebinY(rebinY);
  //DCAwithCut3d->GetZaxis()->SetRange(2,2);
  TH2D * DCAwithCut = (TH2D*) DCAwithCut3d->Project3D("yx"); // This might need to be changed to only use first bin!
  DCAwithCut->RebinY(rebinY);
  
  TH2D * HadronDCA = (TH2D*) inFileTreeData->Get("TrackDCApTnoCutHadron");
  
  TFile * inFileTreeDataPos;
  if(useTOF)
    inFileTreeDataPos = new TFile("local/TreesJan20/TOF/outputTreePbPbDataEta-05neg020.root");
  else
    inFileTreeDataPos = new TFile("local/TreesJan20/noTOF/outputTreePbPbDataEta-05neg020.root");
  //TFile * inFileTreeDataPos = new TFile("local/TreesApr15/outputTreePbPbDataneg020.root");
  EventsTrees += ((TH1D*)inFileTreeDataPos->Get("EventNr"))->Integral();
  EventsTreesBin0 = EventsTrees;
  TH3D * DCAwithCut3dPos = (TH3D*) inFileTreeDataPos->Get("TrackDCApTnoCut");
  //DCAwithCut3dPos->GetZaxis()->SetRange(2,2);
  TH2D * DCAwithCutPos = (TH2D*) DCAwithCut3dPos->Project3D("yx");
  DCAwithCutPos=MirrorDiagram(DCAwithCutPos);
  DCAwithCutPos->RebinY(rebinY);
  DCAwithCut->Add(DCAwithCutPos);
  
  TH2D * HadronDCApos = (TH2D*) inFileTreeDataPos->Get("TrackDCApTnoCutHadron"); HadronDCApos->SetName("HadronDCApos");
  HadronDCApos=MirrorDiagram(HadronDCApos);
  HadronDCA->Add(HadronDCApos);
  HadronDCA->RebinY(rebinY);
  eleHadrons2d=HadronDCA;
  
  TFile * inFileTreeMB;
  //inFileTreeMB = new TFile("local/currentOutput/outputTreeMCMBpos020.root");
  inFileTreeMB = new TFile("local/Mar19Current/outputTreeMCMBPos020.root");
  double EventsTreesMCMB = ((TH1D*)inFileTreeMB->Get("EventNr"))->Integral();
  cout << "Nr of events in outputTreeMCMB: " << EventsTreesMCMB << endl;
  TH3D * DCAwithCut3dTreeMB;
  if(useTOF)
  {DCAwithCut3dTreeMB = (TH3D*) inFileTreeMB->Get("TrackDCApTnoCut"); DCAwithCut3dTreeMB->SetName("DCAwithCut3dTreeMB");}
  else
  {DCAwithCut3dTreeMB = (TH3D*) inFileTreeMB->Get("TrackDCApTnoCutnoTOF"); DCAwithCut3dTreeMB->SetName("DCAwithCut3dTreeMB");}
  //DCAwithCut3dTreeMB->RebinY(rebinY);
  TH2D * DCAwithCutTreeMBAll = (TH2D*) DCAwithCut3dTreeMB->Project3D("yx");  DCAwithCutTreeMBAll->SetName("DCAwithCutTreeMBAll"); DCAwithCutTreeMBAll->RebinY(rebinY);
  DCAwithCut3dTreeMB->GetZaxis()->SetRange(1,1);
  TH2D * DCAwithCutTreeMBD = (TH2D*) DCAwithCut3dTreeMB->Project3D("yx");  DCAwithCutTreeMBD->SetName("DCAwithCutTreeMBD"); DCAwithCutTreeMBD->RebinY(rebinY);
  DCAwithCut3dTreeMB->GetZaxis()->SetRange(2,2);
  TH2D * DCAwithCutTreeMBB = (TH2D*) DCAwithCut3dTreeMB->Project3D("yx");  DCAwithCutTreeMBB->SetName("DCAwithCutTreeMBB"); DCAwithCutTreeMBB->RebinY(rebinY);
  DCAwithCut3dTreeMB->GetZaxis()->SetRange(3,3);
  TH2D * DCAwithCutTreeMBConv = (TH2D*) DCAwithCut3dTreeMB->Project3D("yx");  DCAwithCutTreeMBConv->SetName("DCAwithCutTreeMBConv"); DCAwithCutTreeMBConv->RebinY(rebinY);
  DCAwithCut3dTreeMB->GetZaxis()->SetRange(4,4);
  TH2D * DCAwithCutTreeMBDalitz = (TH2D*) DCAwithCut3dTreeMB->Project3D("yx");  DCAwithCutTreeMBDalitz->SetName("DCAwithCutTreeMBDalitz"); DCAwithCutTreeMBDalitz->RebinY(rebinY);
  
  //TFile * inFileTreeMBPos = new TFile("local/ResolutionCorrection/10correction/outputTreeMCMBneg020.root");
  TFile * inFileTreeMBPos;
  //inFileTreeMBPos = new TFile("local/currentOutput/outputTreeMCMBneg020.root");
  inFileTreeMBPos = new TFile("local/Mar19Current/outputTreeMCMBNeg020.root");
  cout << "Events in positive field MBMC: " << ((TH1D*)inFileTreeMBPos->Get("EventNr"))->Integral() << endl;
  EventsTreesMCMB += ((TH1D*)inFileTreeMBPos->Get("EventNr"))->Integral();
  TH3D * DCAwithCut3dTreeMBPos;
  if(useTOF)
  {DCAwithCut3dTreeMBPos = (TH3D*) inFileTreeMBPos->Get("TrackDCApTnoCut"); DCAwithCut3dTreeMBPos->SetName("DCAwithCut3dTreeMBPos");}
  else
  {DCAwithCut3dTreeMBPos = (TH3D*) inFileTreeMBPos->Get("TrackDCApTnoCutnoTOF"); DCAwithCut3dTreeMBPos->SetName("DCAwithCut3dTreeMBPos");}
  TH2D * DCAwithCutTreeMBAllPos = (TH2D*) DCAwithCut3dTreeMBPos->Project3D("yx");  DCAwithCutTreeMBAllPos->SetName("DCAwithCutTreeMBAllPos");
  DCAwithCut3dTreeMBPos->GetZaxis()->SetRange(1,1);
  TH2D * DCAwithCutTreeMBDPos = (TH2D*) DCAwithCut3dTreeMBPos->Project3D("yx");  DCAwithCutTreeMBDPos->SetName("DCAwithCutTreeMBDPos"); 
  DCAwithCut3dTreeMBPos->GetZaxis()->SetRange(2,2);
  TH2D * DCAwithCutTreeMBBPos = (TH2D*) DCAwithCut3dTreeMBPos->Project3D("yx");  DCAwithCutTreeMBBPos->SetName("DCAwithCutTreeMBBPos"); 
  DCAwithCut3dTreeMBPos->GetZaxis()->SetRange(3,3);
  TH2D * DCAwithCutTreeMBConvPos = (TH2D*) DCAwithCut3dTreeMBPos->Project3D("yx");  DCAwithCutTreeMBConvPos->SetName("DCAwithCutTreeMBConvPos");
  DCAwithCut3dTreeMBPos->GetZaxis()->SetRange(4,4);
  TH2D * DCAwithCutTreeMBDalitzPos = (TH2D*) DCAwithCut3dTreeMBPos->Project3D("yx");  DCAwithCutTreeMBDalitzPos->SetName("DCAwithCutTreeMBDalitzPos"); 
  DCAwithCutTreeMBAllPos=MirrorDiagram(DCAwithCutTreeMBAllPos);
  DCAwithCutTreeMBDPos=MirrorDiagram(DCAwithCutTreeMBDPos);
  DCAwithCutTreeMBBPos=MirrorDiagram(DCAwithCutTreeMBBPos);
  DCAwithCutTreeMBConvPos=MirrorDiagram(DCAwithCutTreeMBConvPos);
  DCAwithCutTreeMBDalitzPos=MirrorDiagram(DCAwithCutTreeMBDalitzPos);
  DCAwithCutTreeMBAllPos->RebinY(rebinY);
  DCAwithCutTreeMBDPos->RebinY(rebinY);
  DCAwithCutTreeMBBPos->RebinY(rebinY);
  DCAwithCutTreeMBConvPos->RebinY(rebinY);
  DCAwithCutTreeMBDalitzPos->RebinY(rebinY);
  DCAwithCutTreeMBAll->Add(DCAwithCutTreeMBAllPos);
  DCAwithCutTreeMBD->Add(DCAwithCutTreeMBDPos);
  DCAwithCutTreeMBB->Add(DCAwithCutTreeMBBPos);
  DCAwithCutTreeMBConv->Add(DCAwithCutTreeMBConvPos);
  DCAwithCutTreeMBDalitz->Add(DCAwithCutTreeMBDalitzPos);
  
  // larger Centrality range for Dalitz here!
  //TFile * inFileTreeMBlargeCent = new TFile("local/ResolutionCorrection/10correction/outputTreeMCMBpos060.root");
  TFile * inFileTreeMBlargeCent;
  //inFileTreeMBlargeCent = new TFile("local/currentOutput/outputTreeMCMBpos060.root");
  inFileTreeMBlargeCent = new TFile("local/Mar19Current/outputTreeMCMBPos060.root");
  TH3D * DCAwithCut3dTreeMBlargeCent;
  if(useTOF)
  {DCAwithCut3dTreeMBlargeCent = (TH3D*) inFileTreeMBlargeCent->Get("TrackDCApTnoCut"); DCAwithCut3dTreeMBlargeCent->SetName("DCAwithCut3dTreeMBlargeCent");}
  else
  {DCAwithCut3dTreeMBlargeCent = (TH3D*) inFileTreeMBlargeCent->Get("TrackDCApTnoCutnoTOF"); DCAwithCut3dTreeMBlargeCent->SetName("DCAwithCut3dTreeMBlargeCent");}
  DCAwithCut3dTreeMBlargeCent->GetZaxis()->SetRange(4,4);
  TH2D * DCAwithCutTreeMBDalitzlargeCent = (TH2D*) DCAwithCut3dTreeMBlargeCent->Project3D("yx");  DCAwithCutTreeMBDalitzlargeCent->SetName("DCAwithCutTreeMBDalitzlargeCent"); DCAwithCutTreeMBDalitzlargeCent->RebinY(rebinY);
  double EventsTreesMCMBlargeCent = ((TH1D*)inFileTreeMBlargeCent->Get("EventNr"))->Integral();
  DCAwithCut3dTreeMBlargeCent->GetZaxis()->SetRange(2,2);
  TH2D * DCAwithCutTreeMBBlargerCentrality = (TH2D*) DCAwithCut3dTreeMBlargeCent->Project3D("yx");  DCAwithCutTreeMBBlargerCentrality->SetName("DCAwithCutTreeMBBlargerCentrality"); DCAwithCutTreeMBBlargerCentrality->RebinY(rebinY);
  
  //TFile * inFileTreeMBPoslargeCent = new TFile("local/ResolutionCorrection/10correction/outputTreeMCMBneg060.root");
  TFile * inFileTreeMBPoslargeCent;
  //inFileTreeMBPoslargeCent = new TFile("local/currentOutput/outputTreeMCMBneg060.root");
  inFileTreeMBPoslargeCent = new TFile("local/Mar19Current/outputTreeMCMBNeg060.root");
  TH3D * DCAwithCut3dTreeMBPoslargeCent;
  if(useTOF)
  {DCAwithCut3dTreeMBPoslargeCent = (TH3D*) inFileTreeMBPoslargeCent->Get("TrackDCApTnoCut"); DCAwithCut3dTreeMBPoslargeCent->SetName("DCAwithCut3dTreeMBPoslargeCent");}
  else
  {DCAwithCut3dTreeMBPoslargeCent = (TH3D*) inFileTreeMBPoslargeCent->Get("TrackDCApTnoCutnoTOF"); DCAwithCut3dTreeMBPoslargeCent->SetName("DCAwithCut3dTreeMBPoslargeCent");}
  DCAwithCut3dTreeMBPoslargeCent->GetZaxis()->SetRange(4,4);
  TH2D * DCAwithCutTreeMBDalitzPoslargeCent = (TH2D*) DCAwithCut3dTreeMBPoslargeCent->Project3D("yx");  DCAwithCutTreeMBDalitzPoslargeCent->SetName("DCAwithCutTreeMBDalitzPoslargeCent");
  DCAwithCut3dTreeMBPoslargeCent->GetZaxis()->SetRange(2,2);
  TH2D * DCAwithCutTreeMBBlargerCentralityPos = (TH2D*) DCAwithCut3dTreeMBPoslargeCent->Project3D("yx");  DCAwithCutTreeMBBlargerCentralityPos->SetName("DCAwithCutTreeMBBlargerCentralityPos");
  
  DCAwithCutTreeMBDalitzPoslargeCent=MirrorDiagram(DCAwithCutTreeMBDalitzPoslargeCent);
  DCAwithCutTreeMBDalitzPoslargeCent->RebinY(rebinY);
  DCAwithCutTreeMBBlargerCentralityPos=MirrorDiagram(DCAwithCutTreeMBBlargerCentralityPos);
  DCAwithCutTreeMBBlargerCentralityPos->RebinY(rebinY);
  EventsTreesMCMBlargeCent += ((TH1D*)inFileTreeMBPoslargeCent->Get("EventNr"))->Integral();
  
  
  DCAwithCutTreeMBDalitzlargeCent->Add(DCAwithCutTreeMBDalitzPoslargeCent);
  DCAwithCutTreeMBBlargerCentrality->Add(DCAwithCutTreeMBBlargerCentralityPos);
  
  if(useLargerCentralityForDalitz){
    DCAwithCutTreeMBDalitz=DCAwithCutTreeMBDalitzlargeCent;
    cout << "Dalitz original statistics: " << DCAwithCutTreeMBDalitz->Integral() << " new: " << DCAwithCutTreeMBDalitzlargeCent->Integral() << endl;}
  // end of larger Centrality range for Dalitz in MB
  
  
  //TFile * inFileTreeenh = new TFile("local/ResolutionCorrection/10correction/outputTreeMCenhpos060.root");
  TFile * inFileTreeenh;
  inFileTreeenh = new TFile("local/Mar19Current/outputTreeMCenhPos060.root");
  cout << "Nr of events in outputTreeMCenh: " << ((TH1D*)inFileTreeenh->Get("EventNr"))->Integral() << endl;
  TH3D * DCAwithCut3dTreeenh;
  if(useTOF)
  {DCAwithCut3dTreeenh = (TH3D*) inFileTreeenh->Get("TrackDCApTnoCut"); DCAwithCut3dTreeenh->SetName("DCAwithCut3dTreeenh");}
  else
  {DCAwithCut3dTreeenh = (TH3D*) inFileTreeenh->Get("TrackDCApTnoCutnoTOF"); DCAwithCut3dTreeenh->SetName("DCAwithCut3dTreeenh");}
  //DCAwithCut3dTreeenh->RebinY(rebinY);
  TH2D * DCAwithCutTreeenhAll = (TH2D*) DCAwithCut3dTreeenh->Project3D("yx");  DCAwithCutTreeenhAll->SetName("DCAwithCutTreeenhAll"); DCAwithCutTreeenhAll->RebinY(rebinY);
  DCAwithCut3dTreeenh->GetZaxis()->SetRange(1,1);
  TH2D * DCAwithCutTreeenhD = (TH2D*) DCAwithCut3dTreeenh->Project3D("yx");  DCAwithCutTreeenhD->SetName("DCAwithCutTreeenhD"); DCAwithCutTreeenhD->RebinY(rebinY);
  DCAwithCut3dTreeenh->GetZaxis()->SetRange(2,2);
  TH2D * DCAwithCutTreeenhB = (TH2D*) DCAwithCut3dTreeenh->Project3D("yx");  DCAwithCutTreeenhB->SetName("DCAwithCutTreeenhB"); DCAwithCutTreeenhB->RebinY(rebinY);
  DCAwithCut3dTreeenh->GetZaxis()->SetRange(3,3);
  TH2D * DCAwithCutTreeenhConv = (TH2D*) DCAwithCut3dTreeenh->Project3D("yx");  DCAwithCutTreeenhConv->SetName("DCAwithCutTreeenhConv"); DCAwithCutTreeenhConv->RebinY(rebinY);
  DCAwithCut3dTreeenh->GetZaxis()->SetRange(4,4);
  TH2D * DCAwithCutTreeenhDalitz = (TH2D*) DCAwithCut3dTreeenh->Project3D("yx");  DCAwithCutTreeenhDalitz->SetName("DCAwithCutTreeenhDalitz"); DCAwithCutTreeenhDalitz->RebinY(rebinY);
  
  
  //TFile * inFileTreeenhPos = new TFile("local/ResolutionCorrection/10correction/outputTreeMCenhneg060.root");
  TFile * inFileTreeenhPos;
  inFileTreeenhPos = new TFile("local/Mar19Current/outputTreeMCenhNeg060.root");
  cout << "Nr of events in outputTreeMCenhPos: " << ((TH1D*)inFileTreeenhPos->Get("EventNr"))->Integral() << endl;
  TH3D * DCAwithCut3dTreeenhPos;
  if(useTOF)
  {DCAwithCut3dTreeenhPos = (TH3D*) inFileTreeenhPos->Get("TrackDCApTnoCut"); DCAwithCut3dTreeenhPos->SetName("DCAwithCut3dTreeenhPos");}
  else
  {DCAwithCut3dTreeenhPos = (TH3D*) inFileTreeenhPos->Get("TrackDCApTnoCutnoTOF"); DCAwithCut3dTreeenhPos->SetName("DCAwithCut3dTreeenhPos");}
  //DCAwithCut3dTreeenhPos->RebinY(rebinY);
  TH2D * DCAwithCutTreeenhAllPos = (TH2D*) DCAwithCut3dTreeenhPos->Project3D("yx");  DCAwithCutTreeenhAllPos->SetName("DCAwithCutTreeenhAllPos"); //DCAwithCutTreeenhAllPos->RebinY(rebinY);
  DCAwithCut3dTreeenhPos->GetZaxis()->SetRange(1,1);
  TH2D * DCAwithCutTreeenhDPos = (TH2D*) DCAwithCut3dTreeenhPos->Project3D("yx");  DCAwithCutTreeenhDPos->SetName("DCAwithCutTreeenhDPos"); //DCAwithCutTreeenhDPos->RebinY(rebinY);
  DCAwithCut3dTreeenhPos->GetZaxis()->SetRange(2,2);
  TH2D * DCAwithCutTreeenhBPos = (TH2D*) DCAwithCut3dTreeenhPos->Project3D("yx");  DCAwithCutTreeenhBPos->SetName("DCAwithCutTreeenhBPos"); //DCAwithCutTreeenhBPos->RebinY(rebinY);
  DCAwithCut3dTreeenhPos->GetZaxis()->SetRange(3,3);
  TH2D * DCAwithCutTreeenhConvPos = (TH2D*) DCAwithCut3dTreeenhPos->Project3D("yx");  DCAwithCutTreeenhConvPos->SetName("DCAwithCutTreeenhConvPos"); //DCAwithCutTreeenhConvPos->RebinY(rebinY);
  DCAwithCut3dTreeenhPos->GetZaxis()->SetRange(4,4);
  TH2D * DCAwithCutTreeenhDalitzPos = (TH2D*) DCAwithCut3dTreeenhPos->Project3D("yx");  DCAwithCutTreeenhDalitzPos->SetName("DCAwithCutTreeenhDalitzPos"); //DCAwithCutTreeenhDalitzPos->RebinY(rebinY);
  DCAwithCutTreeenhAllPos=MirrorDiagram(DCAwithCutTreeenhAllPos);
  DCAwithCutTreeenhDPos=MirrorDiagram(DCAwithCutTreeenhDPos);
  DCAwithCutTreeenhBPos=MirrorDiagram(DCAwithCutTreeenhBPos);
  DCAwithCutTreeenhConvPos=MirrorDiagram(DCAwithCutTreeenhConvPos);
  DCAwithCutTreeenhDalitzPos=MirrorDiagram(DCAwithCutTreeenhDalitzPos);
  DCAwithCutTreeenhAllPos->RebinY(rebinY);
  DCAwithCutTreeenhDPos->RebinY(rebinY);
  DCAwithCutTreeenhBPos->RebinY(rebinY);
  DCAwithCutTreeenhConvPos->RebinY(rebinY);
  DCAwithCutTreeenhDalitzPos->RebinY(rebinY);
  DCAwithCutTreeenhAll->Add(DCAwithCutTreeenhAllPos);
  DCAwithCutTreeenhD->Add(DCAwithCutTreeenhDPos);
  DCAwithCutTreeenhB->Add(DCAwithCutTreeenhBPos);
  DCAwithCutTreeenhConv->Add(DCAwithCutTreeenhConvPos);
  DCAwithCutTreeenhDalitz->Add(DCAwithCutTreeenhDalitzPos);
  TH2D * Bdistenh = (TH2D*)DCAwithCutTreeenhB->Clone();Bdistenh->SetName("Bdistenh");
  
  
  TH2D * dcaBeauty = (TH2D*) ((TH3D*)inFileTreeenh->Get("BeautyMotherCorrelationHalfRAA"))->Project3D("yx"); dcaBeauty->SetName("dcaBeauty");
  TH2D * dcaBeautyNeg = (TH2D*) ((TH3D*)inFileTreeenhPos->Get("BeautyMotherCorrelationHalfRAA"))->Project3D("yx"); dcaBeautyNeg->SetName("dcaBeautyNeg");
  dcaBeauty->Add(MirrorDiagram(dcaBeautyNeg));
  TH2D * dcaCharm = (TH2D*) ((TH3D*)inFileTreeenh->Get("CharmMotherCorrelationHalfRAA"))->Project3D("yx"); dcaCharm->SetName("dcaCharm");
  TH2D * dcaCharmNeg = (TH2D*) ((TH3D*)inFileTreeenhPos->Get("CharmMotherCorrelationHalfRAA"))->Project3D("yx"); dcaCharmNeg->SetName("dcaCharmNeg");
  dcaCharm->Add(MirrorDiagram(dcaCharmNeg));
  
  // Now put here the new data!
  //DCAwithCutTreeMBDalitz=HadronDCA;
  // last bin becomes Hadrons!
  for(int i=1;i<=DCAwithCutTreeMBDalitz->GetNbinsY();i++){
    DCAwithCutTreeMBDalitz->SetBinContent(14, i, HadronDCA->GetBinContent(14,i));
    DCAwithCutTreeMBDalitz->SetBinError(14, i, HadronDCA->GetBinError(14,i));}
  
  DCAwithCutTreeenhB=dcaBeauty;
  DCAwithCutTreeenhD=dcaCharm;
  eleData2d2 = DCAwithCut;
  eleData2d = DCAwithCut;
  eleD2d = DCAwithCutTreeenhD;
  eleB2d = DCAwithCutTreeenhB;
  eleConversion2d = DCAwithCutTreeMBConv;
  eleDalitz2d = DCAwithCutTreeMBDalitz;
  
  //eleD2d->Add(DCAwithCutTreeMBD);
  //eleB2d->Add(DCAwithCutTreeMBB);
  //eleDalitz2d->Add(DCAwithCutTreeenhDalitz);
  //eleConversion2d->Add(DCAwithCutTreeenhConv);
  
  gAliFitDiagramDalitz=eleDalitz2d;
  gAliFitDiagramConversion=eleConversion2d;
  gAliFitDiagramD=eleD2d;
  gAliFitDiagramB=eleB2d;
  gAliFitDiagramData=eleData2d2;
  gAliFitDiagramHadrons=eleHadrons2d;
  
  
  cout << "Fit Diagrams: Dalitz: " << gAliFitDiagramDalitz << " Conv: " << gAliFitDiagramConversion << " B: " << gAliFitDiagramB 
  << " D: " << gAliFitDiagramD << " Data: " << gAliFitDiagramData << endl;
  cout << "contain: Dalitz: " << gAliFitDiagramDalitz->Integral() << " Conv: " << gAliFitDiagramConversion->Integral() << " B: " << gAliFitDiagramB->Integral() 
  << " D: " << gAliFitDiagramD->Integral() << " Data: " << gAliFitDiagramData->Integral() << endl;
  
  TH1D * bOrig =gAliFitDiagramData->ProjectionX(); //////////////////////////////////////////////////
  cout << "B at " << bOrig->GetBinLowEdge(bOrig->FindBin(2.5)) << "GeV-" << bOrig->GetBinLowEdge(bOrig->FindBin(2.5)+1) << " : " <<
  bOrig->GetBinContent(bOrig->FindBin(2.3)) << endl;
  
  
  TH2D * Daown2d = (TH2D*)gAliFitDiagramDalitz->Clone();
  TH2D * Coown2d = (TH2D*)gAliFitDiagramConversion->Clone();
  TH2D * Down2d = (TH2D*)gAliFitDiagramD->Clone();
  TH2D * Bown2d = (TH2D*)gAliFitDiagramB->Clone();
  Daown2d->SetName("Da Own 2d");
  Coown2d->SetName("Co Own 2d");
  Down2d->SetName("D Own 2d");
  Bown2d->SetName("B Own 2d");
  Daown2d->Scale(0);
  Coown2d->Scale(0);
  Down2d->Scale(0);
  Bown2d->Scale(0);
  TH1D * DaownNorm;
  TH1D * CoownNorm;
  TH1D * BownNorm;
  TH1D * DownNorm;
  TH1D * dataown;
  TH1D * dataownDist;
  TH2D * dataown2d=(TH2D*)gAliFitDiagramDalitz->Clone();
  dataown2d->SetName("dataown2d");
  dataown2d->Scale(0);
  

  TH1D * Daown=gAliFitDiagramDalitz->ProjectionY("Da Own",9,18);
  TH1D * Coown=gAliFitDiagramConversion->ProjectionY("Co Own",9,18);
  TH1D * Bown=gAliFitDiagramB->ProjectionY("B Own",9,18);
  TH1D * Down=gAliFitDiagramD->ProjectionY("D Own",9,18);
  
  TH1D * BSpectrum = eleData2d2->ProjectionX("BSpectrum",1, 2);
  TH1D * DSpectrum = eleData2d2->ProjectionX("DSpectrum",1, 2);
  TH1D * DalitzSpectrum = eleData2d2->ProjectionX("DalitzSpectrum",1, 2);
  TH1D * ConversionSpectrum = eleData2d2->ProjectionX("ConversionSpectrum",1, 2);
  TH1D * RatB =  eleData2d2->ProjectionX("ratb",1, 2);
  TH1D * RatD =  eleData2d2->ProjectionX("ratd",1, 2);
  TH1D * RatDalitz =  eleData2d2->ProjectionX("ratdalitz",1, 2);
  TH1D * RatConversion =  eleData2d2->ProjectionX("ratconv",1, 2);
  TH1D * BoverC =  eleData2d2->ProjectionX("b/c",1, 2);
  TH1D * ConvoverDalitz =  eleData2d2->ProjectionX("conv/dalitz",1, 2);
  TH1D * InclusiveoverHFE =  eleData2d2->ProjectionX("inc/hfe",1, 2);
  TH1D * InclusiveoverBG =  eleData2d2->ProjectionX("inc/bg",1, 2);
  
  TH1D * BSpectrumSys = eleData2d2->ProjectionX("BSpectrumSys",1, 2);
  TH1D * DSpectrumSys = eleData2d2->ProjectionX("DSpectrumSys",1, 2);
  TH1D * DBSpectrum = eleData2d2->ProjectionX("DBSpectrum",1, 2);
  TH1D * DBSpectrumSys = eleData2d2->ProjectionX("DBSpectrumSys",1, 2);
  TH1D * DalitzSpectrumSys = eleData2d2->ProjectionX("DalitzSpectrumSys",1, 2);
  TH1D * ConversionSpectrumSys = eleData2d2->ProjectionX("ConversionSpectrumSys",1, 2);
  
  TGraphAsymmErrors * BSpectrumG = new TGraphAsymmErrors(upbin-lowbin+1);
  TGraphAsymmErrors * BSpectrumGsys = new TGraphAsymmErrors(upbin-lowbin+1);
  TGraphAsymmErrors * DSpectrumG = new TGraphAsymmErrors(upbin-lowbin+1);
  TGraphAsymmErrors * DSpectrumGsys = new TGraphAsymmErrors(upbin-lowbin+1);
  TGraphAsymmErrors * DalitzSpectrumG = new TGraphAsymmErrors(upbin-lowbin+1);
  TGraphAsymmErrors * DalitzSpectrumGsys = new TGraphAsymmErrors(upbin-lowbin+1);
  TGraphAsymmErrors * ConversionSpectrumG = new TGraphAsymmErrors(upbin-lowbin+1);
  TGraphAsymmErrors * ConversionSpectrumGsys = new TGraphAsymmErrors(upbin-lowbin+1);
  
  BSpectrumG->SetName("BSpectrumG");
  BSpectrumGsys->SetName("BSpectrumGsys");
  DSpectrumG->SetName("DSpectrumG");
  DSpectrumGsys->SetName("DSpectrumGsys");
  DalitzSpectrumG->SetName("DalitzSpectrumG");
  DalitzSpectrumGsys->SetName("DalitzSpectrumGsys");
  ConversionSpectrumG->SetName("ConversionSpectrumG");
  ConversionSpectrumGsys->SetName("ConversionSpectrumGsys");
  
  
  BSpectrum->Scale(0.0);
  DSpectrum->Scale(0.0);
  DalitzSpectrum->Scale(0.0);
  ConversionSpectrum->Scale(0.0);
  RatB->Scale(0.0);
  RatD->Scale(0.0);
  RatDalitz->Scale(0.0);
  RatConversion->Scale(0.0);
  BoverC->Scale(0.0);
  ConvoverDalitz->Scale(0.0);
  
  cout << "Nbins = " << eleData2d2->GetXaxis()->GetNbins() << endl;
  cout << "Lower Edge: " << eleData2d2->GetXaxis()->GetBinLowEdge(lowbin) << "GeV " << endl;
  cout << "Upper Edge: " << eleData2d2->GetXaxis()->GetBinUpEdge(upbin) << "GeV " << endl;
  
  TH1D **  eleData1d = new TH1D*[upbin-lowbin+1];
  TH1D **  eleD1d = new TH1D*[upbin-lowbin+1];
  TH1D **  eleB1d = new TH1D*[upbin-lowbin+1];
  TH1D **  eleConversion1d = new TH1D*[upbin-lowbin+1];
  TH1D **  eleDalitz1d = new TH1D*[upbin-lowbin+1];
  TH1D **  relErr = new TH1D*[upbin-lowbin+1];
  TH1D **  reconstructed = new TH1D*[upbin-lowbin+1];

  
  TH1D ** eleData1d2 = new TH1D*[upbin-lowbin+1];
  
  TLegend ** FitLegend= new TLegend*[upbin-lowbin+1];
  
  
  //TF1 ** FitFct= new TF1*[upbin-lowbin+1];//necessary?

  
  double * Bparameters = new double[upbin-lowbin+1]; // just to check afterwards
  double * Dparameters = new double[upbin-lowbin+1]; // just to check afterwards
  double * Conversionparameters = new double[upbin-lowbin+1]; // just to check afterwards
  double * Dalitzparameters = new double[upbin-lowbin+1]; // just to check afterwards
  double * ChiSquares = new double[upbin-lowbin+1]; // just to check afterwards
  
  
  
  double Bfactor, Dfactor, Convfactor, Dalitzfactor, TotalNr;
  double Btot, Dtot, Dalitztot, Conversiontot, binWidthpT, binCenter, DBtot;
  
  int nCanv=int(double(upbin-lowbin+1)/4.0+0.8);  
  TCanvas ** CVFits= new TCanvas*[nCanv];
  for(int i=0; i<nCanv; i++)
  {
    CVFits[i] = new TCanvas(Form("Fits%d",i), Form("Fits%d",i), 800, 600);
    CVFits[i]->Divide(2,2);
  }
  TMinuit * minimizingObject = new TMinuit(6); // new for Terhi
  minimizingObject->SetPrintLevel(-1);
  minimizingObject->SetFCN(ReturnChiSq); // new for Terhi

  double expErr;
  double differencetoFct;
  int nBinsIP;
  
  double initTotConv, initTotDalitz, initTotD, initTotB;
  double Berr, Derr, Converr, Dalitzerr, DBerr, DBerrSys;
  double BerrSys, DerrSys, ConverrSys, DalitzerrSys;
  
  // initialize
  Berr = Derr = Converr = Dalitzerr= DBerr = DBerrSys = BerrSys = DerrSys = ConverrSys = DalitzerrSys = 0;
  
  
  // main loop
  for(int i=lowbin;i<=upbin;i++)
  {
    
    
    eleData1d[i-lowbin]=eleData2d->ProjectionY(Form("data%d",i),i, i);
    eleD1d[i-lowbin] = eleD2d->ProjectionY(Form("charm%d",i),i,i);
    eleB1d[i-lowbin] = eleB2d->ProjectionY(Form("beauty%d",i),i,i);
    eleConversion1d[i-lowbin] = eleConversion2d->ProjectionY(Form("conv%d",i),i,i);
    eleDalitz1d[i-lowbin] = eleDalitz2d->ProjectionY(Form("dalitz%d",i),i,i);
    relErr[i-lowbin] = eleDalitz2d->ProjectionY(Form("relErr%d",i),i,i);
    
    // Uncertainties
    for(int j=1;j<=eleD1d[i-lowbin]->GetNbinsX();j++)eleD1d[i-lowbin]->SetBinError(j, TMath::Sqrt(eleD1d[i-lowbin]->GetBinContent(j)));
    for(int j=1;j<=eleD1d[i-lowbin]->GetNbinsX();j++)eleB1d[i-lowbin]->SetBinError(j, TMath::Sqrt(eleB1d[i-lowbin]->GetBinContent(j)));
    for(int j=1;j<=eleD1d[i-lowbin]->GetNbinsX();j++)eleConversion1d[i-lowbin]->SetBinError(j, TMath::Sqrt(eleConversion1d[i-lowbin]->GetBinContent(j)));
    for(int j=1;j<=eleD1d[i-lowbin]->GetNbinsX();j++)eleDalitz1d[i-lowbin]->SetBinError(j, TMath::Sqrt(eleDalitz1d[i-lowbin]->GetBinContent(j)));
    
    //

    
    
    reconstructed[i-lowbin] = eleDalitz2d->ProjectionY(Form("recon%d",i),i,i);
    
    reconstructed[i-lowbin]->Scale(0);
    
    //CleanUpDiagram(eleConversion1d[i-lowbin], -0.03, 0.03);
    //CleanUpDiagram(eleDalitz1d[i-lowbin], -0.03, 0.03);
    
    initTotConv = eleConversion1d[i-lowbin]->Integral();
    initTotDalitz = eleDalitz1d[i-lowbin]->Integral();
    initTotD = eleD1d[i-lowbin]->Integral();
    initTotB = eleB1d[i-lowbin]->Integral();
    
    cout << "\npT from " << eleData2d2->GetXaxis()->GetBinLowEdge(i) << " to " << eleData2d2->GetXaxis()->GetBinUpEdge(i) << endl;
    
    
    
    eleData1d2[i-lowbin] = eleData2d2->ProjectionY(Form("data2%d",i),i, i);
    
    //eleData1d[i-lowbin]->Smooth(50);
    TotalNr=eleData1d2[i-lowbin]->Integral();
    cout << "Total Nr of electrons: " << TotalNr << endl;
    cout << "Total Nr of B in MC: " << eleB1d[i-lowbin]->Integral() << endl;
    cout << "Total Nr of D in MC: " << eleD1d[i-lowbin]->Integral() << endl;
    cout << "Total Nr of conv in MC: " << eleConversion1d[i-lowbin]->Integral() << endl;
    cout << "Total Nr of dalitz in MC: " << eleDalitz1d[i-lowbin]->Integral() << endl;
    
    
    nBinsIP=eleDalitz1d[i-lowbin]->GetNbinsX();
    CVFits[(i-lowbin)/4]->cd((i-lowbin)%4+1);
    
    
    eleData1d2[i-lowbin]->SetAxisRange(-0.1,0.1);
    if(remap)
      eleData1d2[i-lowbin]->SetAxisRange(-17.0,17.0);
    eleData1d2[i-lowbin]->GetXaxis()->SetTitle("DCA (#sigma)");
    eleData1d2[i-lowbin]->GetYaxis()->SetTitle("Counts");
    eleData1d2[i-lowbin]->SetLineWidth(2);
    eleData1d2[i-lowbin]->SetLineColor(1);
    eleData1d2[i-lowbin]->SetMarkerColor(1);
    eleData1d2[i-lowbin]->SetTitle("");
    setOPT_hists((TH1F*) eleData1d2[i-lowbin],"dca #times sgn(charge #times field) (cm)","Entries",505,20,1.0);

    
    
    eleD1d[i-lowbin]->SetAxisRange(-20, 20,"X");
    eleB1d[i-lowbin]->SetLineColor(2);
    
    
    double savB, savD, savDa, savCo;
    minimizingObject->FixParameter(4);
    minimizingObject->mnrset(1);
    double currentChiSq;
    double factors[6];
    
    if(useLikelyhood)currentChiSq = FitLikelyhood(i, factors);   // All the interesting stuff happens here
    else currentChiSq = FitChiSq(i, factors);  // Fit!
    
    
    
    ChiSquares[i-lowbin]=currentChiSq;
    Convfactor=factors[0];
    Dalitzfactor=factors[1];
    Dfactor=factors[2];
    Bfactor=factors[3];
    //LPfactor=factors[4];
    Bparameters[i-lowbin]=Bfactor;
    Dparameters[i-lowbin]=Dfactor;
    Conversionparameters[i-lowbin]=Convfactor;
    Dalitzparameters[i-lowbin]=Dalitzfactor;
    eleConversion1d[i-lowbin]->Scale(Convfactor);
    eleDalitz1d[i-lowbin]->Scale(Dalitzfactor);
    eleD1d[i-lowbin]->Scale(Dfactor);
    eleB1d[i-lowbin]->Scale(Bfactor);
    
    if(DrawParticularBin == i)
    {
      BExpectationOneBin = (TH1D*)gAliTempBMCguess->Clone();
      BExpectationOneBin->SetName("BExpectationOneBin");
      BExpectationOneBin->SetLineColor(kOrange+1);
    }
    
    
    reconstructed[i-lowbin]->Add(eleConversion1d[i-lowbin]);
    reconstructed[i-lowbin]->Add(eleDalitz1d[i-lowbin]);
    reconstructed[i-lowbin]->Add(eleD1d[i-lowbin]);
    reconstructed[i-lowbin]->Add(eleB1d[i-lowbin]);
    
    
    for(int j=1;j<nBinsIP;j++)
    {  // Only one factor for the diagrams because they are already scaled!
      expErr=TMath::Sqrt(eleData1d2[i-lowbin]->GetBinContent(j) + Bfactor * eleB1d[i-lowbin]->GetBinContent(j) + Dfactor * eleD1d[i-lowbin]->GetBinContent(j) + Dalitzfactor * eleDalitz1d[i-lowbin]->GetBinContent(j) + Convfactor * eleConversion1d[i-lowbin]->GetBinContent(j));
      differencetoFct=TMath::Abs(eleData1d2[i-lowbin]->GetBinContent(j)-reconstructed[i-lowbin]->GetBinContent(j));
      //cout << "Expected Error: "  << expErr << "Rel Diff: " << differencetoFct/expErr << endl;
      if(expErr>0)      relErr[i-lowbin]->SetBinContent(j,differencetoFct/TMath::Max(expErr,.01));
				   else relErr[i-lowbin]->SetBinContent(j,0.1);
	relErr[i-lowbin]->SetBinError(j,0.);
    }
    
    parPos[i-lowbin] = new TLine(Bfactor*initTotB/TotalNr,0,Bfactor*initTotB/TotalNr,3);
    cout << "Total ples: " << TotalNr << endl;
    
    
    if(binWidthCorrection)
      binWidthpT=eleData2d2->GetXaxis()->GetBinWidth(i);
    else
      binWidthpT=1.0;
    binCenter=eleData2d2->GetXaxis()->GetBinCenter(i);
    Conversiontot=eleConversion1d[i-lowbin]->Integral();
    Dalitztot=eleDalitz1d[i-lowbin]->Integral();
    Dtot=eleD1d[i-lowbin]->Integral();
    Btot=eleB1d[i-lowbin]->Integral();
    DBtot=Dtot+Btot;
    
    
    //////////////////////////////////////////////////////////////////////////////
    // Start of error calculation code
    //////////////////////////////////////////////////////////////////////////////
    if(errorBarEstimation)
    {
      distOneBin[i-lowbin]->Scale(0);
      distOneBinD->Scale(0);
      distOneBinConv->Scale(0);
      distOneBinDalitz->Scale(0);
      distOneBinDB->Scale(0);
      distOneBinBextra->Scale(0);
      distOneBinDextra->Scale(0);
      distOneBinDBextra->Scale(0);
      
      int ErrEstLowerBin=i;
      DaownNorm=eleDalitz2d->ProjectionY("Da OwnNorm",ErrEstLowerBin,18);
      CoownNorm=eleConversion2d->ProjectionY("Co OwnNorm",ErrEstLowerBin,18);
      BownNorm=eleB2d->ProjectionY("B OwnNorm",ErrEstLowerBin,18);
      DownNorm=eleD2d->ProjectionY("D OwnNorm",ErrEstLowerBin,18);
      dataown=eleDalitz2d->ProjectionY("dataown",ErrEstLowerBin,18);
      dataownDist=eleDalitz2d->ProjectionY("dataownDist",ErrEstLowerBin,18);
      CoownNorm->Scale(1.0/CoownNorm->Integral());
      DownNorm->Scale(1.0/DownNorm->Integral());
      BownNorm->Scale(1.0/BownNorm->Integral());
      DaownNorm->Scale(1.0/DaownNorm->Integral());
      
      for(int errorNr=0;errorNr<errormax;errorNr++)
      {
	cout << "\rError analysis: iteration " << errorNr+1 << "/" << errormax << "   " << flush;
	Coown->Scale(0);
	Bown->Scale(0);
	Daown->Scale(0);
	Down->Scale(0);
	Coown->FillRandom(CoownNorm, int(eleConversion2d->ProjectionY(Form("Copr%d",i),i,i)->Integral())*increaseStatisticsMC); 
	Daown->FillRandom(DaownNorm, int(eleDalitz2d->ProjectionY(Form("Dapr%d",i),i,i)->Integral())*increaseStatisticsMC);
	Down->FillRandom(DownNorm, int(eleD2d->ProjectionY(Form("Dpr%d",i),i,i)->Integral())*increaseStatisticsMC);
	Bown->FillRandom(BownNorm, int(eleB2d->ProjectionY(Form("Bpr%d",i),i,i)->Integral())*increaseStatisticsMC);
	dataownDist->Scale(0);
	dataown->Scale(0);
	dataownDist->Add(CoownNorm,Conversiontot);
	dataownDist->Add(DaownNorm,Dalitztot);
	dataownDist->Add(DownNorm,Dtot);
	dataownDist->Add(BownNorm,Btot);
	dataown->FillRandom(dataownDist,int(TotalNr)*increaseStatistics);
	for(int ownBins=1;ownBins<nBinsIP;ownBins++)
	{
	  dataown2d->SetBinContent(dataown2d->GetBin(i,ownBins),dataown->GetBinContent(ownBins));
	  Coown2d->SetBinContent(dataown2d->GetBin(i,ownBins),Coown->GetBinContent(ownBins));
	  Daown2d->SetBinContent(dataown2d->GetBin(i,ownBins),Daown->GetBinContent(ownBins));
	  Down2d->SetBinContent(dataown2d->GetBin(i,ownBins),Down->GetBinContent(ownBins));
	  Bown2d->SetBinContent(dataown2d->GetBin(i,ownBins),Bown->GetBinContent(ownBins));
	}
	gAliFitDiagramDalitz=Daown2d;
	gAliFitDiagramConversion=Coown2d;
	gAliFitDiagramD=Down2d;
	gAliFitDiagramB=Bown2d;
	gAliFitDiagramData=dataown2d;
	
	if(useLikelyhood)currentChiSq = FitLikelyhood(i, factors);
          else currentChiSq = FitChiSq(i, factors);  // Fit!
	savCo=factors[0];
	savDa=factors[1];
	savD=factors[2];
	savB=factors[3];
	
	distOneBin[i-lowbin]->Fill(savB*initTotB/TotalNr);
	distOneBinD->Fill(savD*initTotD/TotalNr);
	distOneBinDalitz->Fill(savDa*initTotDalitz/TotalNr);
	distOneBinConv->Fill(savCo*initTotConv/TotalNr);
	distOneBinDB->Fill((savB*initTotB+savD*initTotD)/DBtot);
      }
      cout << endl;
      
      
      
      if(distOneBin[i-lowbin]->GetMean()>0){
	Berr = distOneBin[i-lowbin]->GetRMS()/(Bfactor*initTotB/TotalNr);
	BerrSys = TMath::Abs(distOneBin[i-lowbin]->GetMean()-Bfactor*initTotB/TotalNr)/(Bfactor*initTotB/TotalNr);}
	
      if(distOneBinD->GetMean()>0){
	Derr = distOneBinD->GetRMS()/(Dfactor*initTotD/TotalNr);
	DerrSys = TMath::Abs(distOneBinD->GetMean()-Dfactor*initTotD/TotalNr)/(Dfactor*initTotD/TotalNr);}
	
      DBerr=distOneBinDB->GetRMS();   DBerrSys=TMath::Abs(distOneBinDB->GetMean()-1.0);
      
      if(distOneBinDalitz->GetMean()>0){
	Dalitzerr = distOneBinDalitz->GetRMS()/(Dalitzfactor*initTotDalitz/TotalNr);
	DalitzerrSys = TMath::Abs(distOneBinDalitz->GetMean()-Dalitzfactor*initTotDalitz/TotalNr)/(Dalitzfactor*initTotDalitz/TotalNr);}
	
      if(distOneBinConv->GetMean()>0){
	Converr = distOneBinConv->GetRMS()/(Convfactor*initTotConv/TotalNr);
	ConverrSys = TMath::Abs(distOneBinConv->GetMean()-Convfactor*initTotConv/TotalNr)/(Convfactor*initTotConv/TotalNr);}

	double TotalNrExtra=0.;
	cout << "relative total d->e uncertainty: " << TMath::Sqrt(Derr*Derr + DerrSys*DerrSys) << endl;
	if(TMath::Sqrt(Derr*Derr + DerrSys*DerrSys)>0.5)
	{
	  cout << "D->e  error too large! Needs additional error estimation! " << Derr << " , " << DerrSys << endl;
	  cout << "B->e error of " << Berr << " , " << BerrSys << " probably underestimated" << endl;
	  double DoverBExtra = TMath::Max(0.5, 1./BoverC->GetBinContent(i-1));
	  cout << "D/B ratio used: " << DoverBExtra << endl;
	  
	  
	  for(int errorNr=0;errorNr<errormax;errorNr++)
	  {
	    ExtraErrorTestNeeded=true;
	    cout << "\rAdditional error analysis: iteration " << errorNr+1 << "/" << errormax << "   " << flush;
	    Coown->Scale(0);
	    Bown->Scale(0);
	    Daown->Scale(0);
	    Down->Scale(0);
	    Coown->FillRandom(CoownNorm, int(eleConversion2d->ProjectionY(Form("Copr%d",i),i,i)->Integral())); 
	    Daown->FillRandom(DaownNorm, int(eleDalitz2d->ProjectionY(Form("Dapr%d",i),i,i)->Integral()));
	    Down->FillRandom(DownNorm, int(eleD2d->ProjectionY(Form("Dpr%d",i),i,i)->Integral()));
	    Bown->FillRandom(BownNorm, int(eleB2d->ProjectionY(Form("Bpr%d",i),i,i)->Integral()));
	    dataownDist->Scale(0);
	    dataown->Scale(0);
	    dataownDist->Add(CoownNorm,Conversiontot);
	    dataownDist->Add(DaownNorm,Dalitztot);
	    dataownDist->Add(DownNorm,Btot*DoverBExtra);
	    dataownDist->Add(BownNorm,Btot);
	    TotalNrExtra = int(Conversiontot+Dalitztot+Btot*DoverBExtra*0.+Btot);
	    dataown->FillRandom(dataownDist,TotalNrExtra);
	    for(int ownBins=1;ownBins<nBinsIP;ownBins++)
	    {
	      dataown2d->SetBinContent(dataown2d->GetBin(i,ownBins),dataown->GetBinContent(ownBins));
	      Coown2d->SetBinContent(dataown2d->GetBin(i,ownBins),Coown->GetBinContent(ownBins));
	      Daown2d->SetBinContent(dataown2d->GetBin(i,ownBins),Daown->GetBinContent(ownBins));
	      Down2d->SetBinContent(dataown2d->GetBin(i,ownBins),Down->GetBinContent(ownBins));
	      Bown2d->SetBinContent(dataown2d->GetBin(i,ownBins),Bown->GetBinContent(ownBins));
	    }
	    gAliFitDiagramDalitz=Daown2d;
	    gAliFitDiagramConversion=Coown2d;
	    gAliFitDiagramD=Down2d;
	    gAliFitDiagramB=Bown2d;
	    gAliFitDiagramData=dataown2d;
	    
	    if(useLikelyhood)currentChiSq = FitLikelyhood(i, factors);
	    else currentChiSq = FitChiSq(i, factors);  // Fit!
	    savCo=factors[0];
	    savDa=factors[1];
	    savD=factors[2];
	    savB=factors[3];
	    
	    distOneBinBextra->Fill(savB*initTotB/TotalNrExtra);
	    distOneBinDextra->Fill(savD*initTotD/TotalNrExtra);
	    distOneBinDBextra->Fill((savB*initTotB+savD*initTotD)/DBtot);
	  }
	  cout << endl;
	  double BerrExtra=0;
	  double BerrSysExtra=0;
	  if(distOneBinBextra->GetMean()>0){
	    BerrExtra = distOneBinBextra->GetRMS()/(Bfactor*initTotB/TotalNrExtra);
	    BerrSysExtra = TMath::Abs(distOneBinBextra->GetMean()-Bfactor*initTotB/TotalNrExtra)/(Bfactor*initTotB/TotalNrExtra);}
	   
	  double DerrExtra=0;
	  double DerrSysExtra=0;
	  if(distOneBinDextra->GetMean()>0){
	    DerrExtra = distOneBinDextra->GetRMS()/(Bfactor*initTotB*DoverBExtra/TotalNrExtra);
	    DerrSysExtra = TMath::Abs(distOneBinDextra->GetMean()-Bfactor*initTotB*DoverBExtra/TotalNrExtra)/(Bfactor*initTotB*DoverBExtra/TotalNrExtra);}
	
	  double DBerrExtra=0;
	  double DBerrSysExtra=0;
	  if(distOneBinDBextra->GetMean()>0){
	    DBerrExtra=distOneBinDBextra->GetRMS();
	    DBerrSysExtra=TMath::Abs(distOneBinDBextra->GetMean()-1.0);}
	  
	
	
	  cout << "B->e new error of " << BerrExtra << " , " << BerrSysExtra << endl;
	  cout << "New D errors " << DerrExtra << " , " << DerrSysExtra << endl;
	  if(TMath::Sqrt(DerrExtra*DerrExtra + DerrSysExtra*DerrSysExtra)>0.5)
	  {
	    cout << "CAREFUL! Extra Error analysis insufficient!" << endl;
	  }
	  else
	  {
	    ExtraErrorTestSuccessful=true; 
	    cout << "Extra Error analysis has helped!" << endl;
	  }
	    if(BerrExtra>Berr)Berr=BerrExtra;
	    if(BerrSysExtra>BerrSys)BerrSys=BerrSysExtra;
	    if(DBerrExtra>DBerr){DBerr=DBerrExtra;}
	    if(DBerrSysExtra>DBerrSys)DBerrSys=DBerrSysExtra;
	}
      
      
      // cleanup
      delete DaownNorm;
      delete CoownNorm;
      delete BownNorm;
      delete DownNorm;
      delete dataown;
      delete dataownDist;
      
      gAliFitDiagramDalitz=eleDalitz2d;
      gAliFitDiagramConversion=eleConversion2d;
      gAliFitDiagramD=eleD2d;
      gAliFitDiagramB=eleB2d;
      gAliFitDiagramData=eleData2d2;
    }
    
    //////////////////////////////////////////////////////////////////////////////
    // End of error calculation code
    //////////////////////////////////////////////////////////////////////////////
    
    if(!errorBarEstimation){
      Berr = 0.0;
      Derr = 0.0;
    }
    
    BSpectrumG->SetPoint(i-lowbin, binCenter, Btot/binWidthpT/events);
    BSpectrumG->SetPointError(i-lowbin, binWidthpT/2.0, eleData2d2->GetXaxis()->GetBinWidth(i)/2.0, Berr*Btot/binWidthpT/events, Berr*Btot/binWidthpT/events);
    BSpectrumGsys->SetPoint(i-lowbin, binCenter, Btot/binWidthpT/events);
    BSpectrumGsys->SetPointError(i-lowbin, binWidthpT/2.0, eleData2d2->GetXaxis()->GetBinWidth(i)/2.0, BerrSys*Btot/binWidthpT/events, BerrSys*Btot/binWidthpT/events);
    DSpectrumG->SetPoint(i-lowbin, binCenter, Dtot/binWidthpT/events);
    DSpectrumG->SetPointError(i-lowbin, eleData2d2->GetXaxis()->GetBinWidth(i)/2.0, binWidthpT/2.0, Derr*Dtot/binWidthpT/events, Derr*Dtot/binWidthpT/events);
    DSpectrumGsys->SetPoint(i-lowbin, binCenter, Dtot/binWidthpT/events);
    DSpectrumGsys->SetPointError(i-lowbin, eleData2d2->GetXaxis()->GetBinWidth(i)/2.0, binWidthpT/2.0, DerrSys*Dtot/binWidthpT/events, DerrSys*Dtot/binWidthpT/events);
    
    DalitzSpectrumG->SetPoint(i-lowbin, binCenter, Dalitztot/binWidthpT/events);
    DalitzSpectrumG->SetPointError(i-lowbin, binWidthpT/2.0, eleData2d2->GetXaxis()->GetBinWidth(i)/2.0, Dalitzerr*Dalitztot/binWidthpT/events, Dalitzerr*Dalitztot/binWidthpT/events);
    ConversionSpectrumG->SetPoint(i-lowbin, binCenter, Conversiontot/binWidthpT/events);
    ConversionSpectrumG->SetPointError(i-lowbin, binWidthpT/2.0, eleData2d2->GetXaxis()->GetBinWidth(i)/2.0, Converr*Conversiontot/binWidthpT/events, Converr*Conversiontot/binWidthpT/events);
    DalitzSpectrumGsys->SetPoint(i-lowbin, binCenter, Dalitztot/binWidthpT/events);
    DalitzSpectrumGsys->SetPointError(i-lowbin, binWidthpT/2.0, eleData2d2->GetXaxis()->GetBinWidth(i)/2.0, DalitzerrSys*Dalitztot/binWidthpT/events, DalitzerrSys*Dalitztot/binWidthpT/events);
    ConversionSpectrumGsys->SetPoint(i-lowbin, binCenter, Conversiontot/binWidthpT/events);
    ConversionSpectrumGsys->SetPointError(i-lowbin, binWidthpT/2.0, eleData2d2->GetXaxis()->GetBinWidth(i)/2.0, ConverrSys*Conversiontot/binWidthpT/events, ConverrSys*Conversiontot/binWidthpT/events);
    
    InclusiveoverHFE->SetBinContent(i,(TotalNr)/(Btot+Dtot));
    InclusiveoverBG->SetBinContent(i,(TotalNr)/(Dalitztot+Conversiontot));
    
    BSpectrum->Fill(binCenter, Btot/*/binWidthpT*/);
    BSpectrum->SetBinError(i,Berr*Btot/*/binWidthpT*/);
    BSpectrumSys->Fill(binCenter, Btot/*/binWidthpT*/);
    BSpectrumSys->SetBinError(i,BerrSys*Btot/*/binWidthpT*/);
    DSpectrum->Fill(binCenter, Dtot/*/binWidthpT*/);
    DSpectrum->SetBinError(i,Derr*Dtot/*/binWidthpT*/);
    DSpectrumSys->Fill(binCenter, Dtot/*/binWidthpT*/);
    DSpectrumSys->SetBinError(i,DerrSys*Dtot/*/binWidthpT*/);
    DBSpectrum->Fill(binCenter,DBtot);
    DBSpectrum->SetBinError(i,DBerr*DBtot/*/binWidthpT*/);
    DBSpectrumSys->Fill(binCenter,DBtot);
    DBSpectrumSys->SetBinError(i,DBerrSys*DBtot/*/binWidthpT*/);

    ConversionSpectrum->Fill(binCenter, Conversiontot);
    ConversionSpectrum->SetBinError(i,Converr*Conversiontot);
    DalitzSpectrum->Fill(binCenter, Dalitztot);
    DalitzSpectrum->SetBinError(i,Dalitzerr*Dalitztot);
    ConversionSpectrumSys->Fill(binCenter, Conversiontot);
    ConversionSpectrumSys->SetBinError(i,ConverrSys*Conversiontot);
    DalitzSpectrumSys->Fill(binCenter, Dalitztot);
    DalitzSpectrumSys->SetBinError(i,DalitzerrSys*Dalitztot);
    
    
    RatConversion->Fill(binCenter, Conversiontot/TotalNr);
    RatD->Fill(binCenter, Dtot/TotalNr);
    RatB->Fill(binCenter, Btot/TotalNr);
    RatDalitz->Fill(binCenter, Dalitztot/TotalNr);
    if(Dtot==0)Dtot=0.001;if(Dalitztot==0)Dalitztot=0.001;
			 BoverC->Fill(binCenter, Btot/Dtot);
    ConvoverDalitz->Fill(binCenter, Conversiontot/Dalitztot);
    
    eleConversion1d[i-lowbin]->SetLineColor(kGray+2);
    eleDalitz1d[i-lowbin]->SetLineColor(4);
    eleD1d[i-lowbin]->SetLineColor(kGreen+1);
    eleD1d[i-lowbin]->SetFillStyle(0);
    eleD1d[i-lowbin]->SetFillColor(kGreen+1);
    eleB1d[i-lowbin]->SetLineColor(2);
    eleB1d[i-lowbin]->SetFillStyle(0);
    eleB1d[i-lowbin]->SetFillColor(2);
    reconstructed[i-lowbin]->SetLineColor(kYellow+2);
    
    eleConversion1d[i-lowbin]->SetMarkerColor(kGray+2);
    eleDalitz1d[i-lowbin]->SetMarkerColor(4);
    eleD1d[i-lowbin]->SetMarkerColor(kGreen+1);
    eleB1d[i-lowbin]->SetMarkerColor(2);
    reconstructed[i-lowbin]->SetMarkerColor(kYellow+2);
    relErr[i-lowbin]->SetMarkerColor(6);
       
    eleConversion1d[i-lowbin]->SetLineWidth(2);
    eleDalitz1d[i-lowbin]->SetLineWidth(2);
    eleD1d[i-lowbin]->SetLineWidth(2);
    eleB1d[i-lowbin]->SetLineWidth(2);
    reconstructed[i-lowbin]->SetLineWidth(2);
    relErr[i-lowbin]->SetLineWidth(2);
    
    eleConversion1d[i-lowbin]->SetMarkerStyle(markers[7]);
    eleDalitz1d[i-lowbin]->SetMarkerStyle(markers[1]);
    eleD1d[i-lowbin]->SetMarkerStyle(markers[2]);
    eleB1d[i-lowbin]->SetMarkerStyle(markers[3]);
    reconstructed[i-lowbin]->SetMarkerStyle(markers[4]);
    relErr[i-lowbin]->SetMarkerStyle(markers[5]);
    
    eleConversion1d[i-lowbin]->SetMarkerSize(1);
    eleDalitz1d[i-lowbin]->SetMarkerSize(1);
    eleD1d[i-lowbin]->SetMarkerSize(1);
    eleB1d[i-lowbin]->SetMarkerSize(1);
    reconstructed[i-lowbin]->SetMarkerSize(1);
    relErr[i-lowbin]->SetMarkerSize(1);

    eleData1d2[i-lowbin]->Draw("pe");
    gPad->SetLogy();
    
    eleConversion1d[i-lowbin]->Draw("SAMEP");
    eleDalitz1d[i-lowbin]->Draw("SAMEP");
    eleD1d[i-lowbin]->Draw("SAMEP");
    eleB1d[i-lowbin]->Draw("SAMEP");
    
    relErr[i-lowbin]->SetLineColor(6);
    relErr[i-lowbin]->Draw("SAME");
    reconstructed[i-lowbin]->Draw("SAMEP");
    
    
    
    FitLegend[i-lowbin]=plotLegend("right_top"," ",0.8,0.9,+0.08,+0.01,Form("%.1f<#it{p}_{T}<%.1f GeV/#it{c}",eleData2d2->GetXaxis()->GetBinLowEdge(i),eleData2d2->GetXaxis()->GetBinUpEdge(i)),1);
    FitLegend[i-lowbin]->AddEntry(eleData1d2[i-lowbin], "Total electrons", "PE");
    FitLegend[i-lowbin]->AddEntry(reconstructed[i-lowbin], "Fit", "P");
    FitLegend[i-lowbin]->AddEntry(whiteLine, "", "L");
    FitLegend[i-lowbin]->AddEntry(whiteLine, "HIJING+PYTHIA:", "L");
    FitLegend[i-lowbin]->AddEntry(eleConversion1d[i-lowbin], "Conversion electrons", "P");
    FitLegend[i-lowbin]->AddEntry(eleDalitz1d[i-lowbin], "Dalitz electrons", "P");
    FitLegend[i-lowbin]->AddEntry(eleD1d[i-lowbin], "c #rightarrow e", "P");
    FitLegend[i-lowbin]->AddEntry(eleB1d[i-lowbin], "b (#rightarrow c) #rightarrow e", "P");
    FitLegend[i-lowbin]->AddEntry(relErr[i-lowbin], "relative Error", "L");
    
    
    
    FitLegend[i-lowbin]->Draw("SAME");
    
  }
  
  
  TCanvas * CVSpectrum = new TCanvas("spectrumcv", "spectrumcv", 800, 600);
  setOPT_canvas(CVSpectrum);
  TLegend * specLeg = plotLegend("right_top"," ",0.8,0.5,+0.0,+0.00,"",1);//new TLegend(0.6,0.8,0.95,0.95);
  RatB->SetLineColor(2);
  RatD->SetLineColor(3);
  RatConversion->SetLineColor(4);
  RatDalitz->SetLineColor(48);
  
  setOPT_graph(BSpectrumG, "p_{T}", "dN/dp_{T}");
  setOPT_graph(BSpectrumGsys, "p_{T}", "dN/dp_{T}");
  setOPT_graph(DSpectrumG, "p_{T}", "dN/dp_{T}");
  setOPT_graph(DSpectrumGsys, "p_{T}", "dN/dp_{T}");
  setOPT_graph(DalitzSpectrumG, "p_{T}", "dN/dp_{T}");
  setOPT_graph(ConversionSpectrumG, "p_{T}", "dN/dp_{T}");
  setOPT_graph(DalitzSpectrumGsys, "p_{T}", "dN/dp_{T}");
  setOPT_graph(ConversionSpectrumGsys, "p_{T}", "dN/dp_{T}");
  BSpectrumGsys->SetMarkerSize(0);
  BSpectrumGsys->SetFillColor(kRed);
  BSpectrumGsys->SetFillStyle(3002);
  DSpectrumGsys->SetMarkerSize(0);
  DSpectrumGsys->SetFillColor(kGreen);
  DSpectrumGsys->SetFillStyle(3001);
  ConversionSpectrumGsys->SetMarkerSize(0);
  ConversionSpectrumGsys->SetFillColor(kGray+2);
  ConversionSpectrumGsys->SetFillStyle(3102);
  DalitzSpectrumGsys->SetMarkerSize(0);
  DalitzSpectrumGsys->SetFillColor(kOrange);
  DalitzSpectrumGsys->SetFillStyle(3101);
  specLeg->AddEntry(BSpectrumG, "Fit b#rightarrow e Spectrum");
  specLeg->AddEntry(DSpectrumG, "Fit c#rightarrow e Spectrum");
  //specLeg->AddEntry(DalitzSpectrumG, "Dalitz e Spectrum");
  //specLeg->AddEntry(ConversionSpectrumG, "Conversion e Spectrum");
  
  
  
  TH1F *myBlankHisto = new TH1F("myBlankHisto","",20,0,8);
  myBlankHisto->SetMaximum(5E-1);
  myBlankHisto->SetMinimum(5E-5);
  setOPT_hists(myBlankHisto);
  myBlankHisto->SetYTitle("1/N_{evts} dN/dp_{t} (1/GeV/#it{c})");
  myBlankHisto->SetXTitle("p_{t} (GeV/#it{c})");
  myBlankHisto->SetNdivisions(505,"y");
  myBlankHisto->Draw();
  gPad->SetLogy();
  
  
  BSpectrum->SetAxisRange(0.0,10.0);
  BSpectrum->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  BSpectrum->SetLineColor(kBlue);
  DSpectrum->SetLineColor(kBlue);
  
  DSpectrumG->SetMarkerColor(kGreen);
  BSpectrumG->SetMarkerColor(kRed);
  ConversionSpectrumG->SetMarkerColor(kGray+2);
  DalitzSpectrumG->SetMarkerColor(kOrange);
  //DalitzSpectrumG->Draw("SAMEP");
  //ConversionSpectrumG->Draw("SAMEP");
  //DalitzSpectrumGsys->Draw("SAMEP2");
  //ConversionSpectrumGsys->Draw("SAMEP2");
  BSpectrumGsys->Draw("SAMEP2");
  DSpectrumGsys->Draw("SAMEP2");
  BSpectrumG->Draw("SAMEP");
  DSpectrumG->Draw("SAMEP");
  
  
  specLeg->Draw("SAME");
  
  TF1 * compIncHFEfct = new TF1("cihfefct", "0.879432+x*1.32325+x**2*(-0.978663)+0.290138*x**3-0.0385313*x**4+0.00190326*x**5", 1.0, 5.0);
  
  TCanvas * CVboverc = new TCanvas("BoverC", "BoverC", 800, 600);
  TLegend * RatiosLeg =plotLegend("right_top"," ",0.8,0.5,+0.0,+0.00,"",1);
  BoverC->SetAxisRange(0.0,10.0);
  BoverC->SetAxisRange(0.0,5.0, "Y");
  InclusiveoverHFE->SetLineColor(4);
  InclusiveoverHFE->SetLineWidth(2);
  InclusiveoverBG->SetLineColor(kBlack);
  InclusiveoverBG->SetLineWidth(2);
  RatiosLeg->AddEntry(BoverC, "b #rightarrow e / c #rightarrow e");
  RatiosLeg->AddEntry(ConvoverDalitz, "Conversion / Dalitz");
  //RatiosLeg->AddEntry(InclusiveoverHFE, "Inclusive / HFE");
  RatiosLeg->AddEntry(InclusiveoverBG, "Inclusive / Background");
  BoverC->SetLineColor(2);
  ConvoverDalitz->SetLineColor(3);
  setOPT_hists((TH1F*)BoverC, "p_{T} (GeV/c)", "counts", 510, 20, 1.3, kRed+1);
  setOPT_hists((TH1F*)ConvoverDalitz, "p_{T} (GeV/c)", "counts", 510, 20, 1.3, kGreen+1);
  setOPT_hists((TH1F*)InclusiveoverHFE, "p_{T} (GeV/c)", "counts", 510, 20, 1.3, kBlack);
  setOPT_hists((TH1F*)InclusiveoverBG, "p_{T} (GeV/c)", "counts", 510, 20, 1.3, kBlack);
  for(int bn=1;bn<lowbin;bn++)
  {
    BoverC->SetBinContent(bn,0.);
    ConvoverDalitz->SetBinContent(bn,0.);
    InclusiveoverHFE->SetBinContent(bn,0.);
    InclusiveoverBG->SetBinContent(bn,0.);
  }
  BoverC->SetTitle("");
  BoverC->Draw();
  BoverC->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  ConvoverDalitz->Draw("SAME");
  RatiosLeg->Draw("SAME");
  //InclusiveoverHFE->Draw("SAME");
  compIncHFEfct->Draw("SAME");
  InclusiveoverBG->Draw("SAME");
  CVboverc->Flush();
  
  TLegend ** errleg = new TLegend*[upbin-lowbin+1];
  TH1F * redLine = new TH1F();
  setOPT_hists(redLine);
  redLine->SetLineColor(kRed);
  redLine->SetMarkerSize(0.0);
  
  TCanvas ** CVErrors= new TCanvas*[nCanv];
  if(errorBarEstimation)
  {
    for(int i=0; i<nCanv; i++)
    {
      CVErrors[i] = new TCanvas(Form("ParVar%d",i), Form("ParVar%d",i), 800, 600);
    setOPT_canvas(CVErrors[i]);
    CVErrors[i]->Divide(2,2);
    }
    
    for(int i=0;i<upbin-lowbin+1;i++){
      errleg[i] = plotLegend("right_top"," ",1.0,0.8,-0.05,-0.15,Form("%.2f < p_{T} < %.2f GeV/c",BSpectrum->GetXaxis()->GetBinLowEdge(i+lowbin), BSpectrum->GetXaxis()->GetBinUpEdge(i+lowbin)),1);
      errleg[i]->AddEntry(distOneBin[0], "Model Results b #rightarrow e fit parameter");
      errleg[i]->AddEntry(redLine, "Measurement");
      CVErrors[i/4]->cd(i%4+1);
      setOPT_hists((TH1F*)distOneBin[i],"b #rightarrow e / all e","counts");
      distOneBin[i]->SetMarkerSize(0.01);
      distOneBin[i]->SetTitle("");distOneBin[i]->SetName("");
      distOneBin[i]->SetLineColor(kBlue);
      distOneBin[i]->SetLineWidth(2);
      distOneBin[i]->Draw();
      //parPos[i]->SetMarkerSize(0.01);
      parPos[i]->SetLineColor(2);
      parPos[i]->SetLineWidth(2);
      parPos[i]->Draw();
      errleg[i]->Draw("SAME");
    }
  }
  
  
  for(int i=0; i<(upbin-lowbin+1);i++)
    cout << "#" << i << " p_t " << BSpectrum->GetBinCenter(i+lowbin) << " Conv/Dal/D/B " << Conversionparameters[i] << " " << Dalitzparameters[i] << " " << Dparameters[i] << " " << Bparameters[i] << " Chi2: " << ChiSquares[i] << endl;
  cout << "Total: " << endl;
  for(int i=0; i<(upbin-lowbin+1);i++)
    cout << "#" << i << " p_t " << BSpectrum->GetBinCenter(i+lowbin) << " Conv/Dal/D/B " << eleConversion1d[i]->Integral() << " " << eleDalitz1d[i]->Integral() << " " << eleD1d[i]->Integral() << " " << eleB1d[i]->Integral() << " Data: " << eleData1d2[i]->Integral(1,eleData1d2[i]->GetNbinsX()) << endl;

  //cout << "B moment 0: " << momentsB[0] << " B moment 1: " << momentsB[1] << " B moment 2: " << momentsB[2] << " B moment 4: " << momentsB[4] << endl;
  //cout << "D moment 0: " << momentsD[0] << " D moment 1: " << momentsD[1] << " D moment 2: " << momentsD[2] << " D moment 4: " << momentsD[4] << endl;
  //cout << "Dalitz moment 0: " << momentsDalitz[0] << " Daliz moment 1: " << momentsDalitz[1] << " Dalitz moment 2: " << momentsDalitz[2] << " Dalitz moment 4: " << momentsDalitz[4] << endl;
  //cout << "Tot moment 0: " << momentsTot[0] << " Tot moment 1: " << momentsTot[1] << " Tot moment 2: " << momentsTot[2] << " Tot moment 4: " << momentsTot[4] << endl;
  
  
  TH1D * evts = new TH1D("events", "", 1, 0., 1.);
  evts->Fill(0.5,EventsTrees);
  TH1D * evtsb0 = new TH1D("evtsb0", "", 1, 0., 1.);
  evtsb0->Fill(0.5,EventsTreesBin0);
  
  TH1D * totalSpectrum = eleData2d->ProjectionX();totalSpectrum->SetName("totalSpectrum");
  
  
  new TCanvas("DCAFit1520", "ExampleBin", 800, 600);
  eleData1d2[DrawParticularBin-lowbin]->Draw("pe");
  gPad->SetLogy();
  eleConversion1d[DrawParticularBin-lowbin]->Draw("SAMEP");
  eleDalitz1d[DrawParticularBin-lowbin]->Draw("SAMEP");
  eleD1d[DrawParticularBin-lowbin]->Draw("SAMEP");
  eleB1d[DrawParticularBin-lowbin]->Draw("SAMEP");
  relErr[DrawParticularBin-lowbin]->Draw("SAME");
  reconstructed[DrawParticularBin-lowbin]->Draw("SAMEP");
  FitLegend[DrawParticularBin-lowbin]->Draw("SAME");
  TLatex *   tex = new TLatex(0.13, 0.90, 1. ? "ALICE Preliminary" : "ALICE");
  tex->SetNDC();
  tex->SetTextFont(42);
  //tex->Draw("SAME");
  gPad->SetBottomMargin(0.11);
  TLatex *   texSystem1 = new TLatex(0.13, 0.81, "Pb-Pb, 0-20% centrality");
  texSystem1->SetNDC();
  texSystem1->SetTextFont(42);
  texSystem1->SetTextSize(0.045);
  texSystem1->Draw("SAME");
  TLatex *   texSystem2 = new TLatex(0.13, 0.76, "#sqrt{#it{s}}_{NN} = 2.76 TeV");
  texSystem2->SetNDC();
  texSystem2->SetTextFont(42);
  texSystem2->SetTextSize(0.045);
  texSystem2->Draw("SAME");
  //BExpectationOneBin->Draw("same");
  
  
  new TCanvas("Uncert1520", "Uncert1520", 800, 600);
  distOneBin[DrawParticularBin-lowbin]->Draw();
  parPos[DrawParticularBin-lowbin]->Draw();
  TLatex *   tex2 = new TLatex(0.33, 0.90, 1. ? "ALICE Preliminary (requested)" : "ALICE");
  tex2->SetNDC();
  tex2->SetTextFont(42);
  //tex2->Draw("SAME");
  gPad->SetBottomMargin(0.11);
  TLatex *   texSystem3 = new TLatex(0.5, 0.81, "#bf{Pb-Pb}, 0-20% centrality");
  texSystem3->SetNDC();
  texSystem3->SetTextFont(42);
  //texSystem1->SetTextSize(0.05);
  texSystem3->Draw("SAME");
  TLatex *   texSystem4 = new TLatex(0.5, 0.76, "#sqrt{#it{s}}_{NN} = 2.76 TeV");
  texSystem4->SetNDC();
  texSystem4->SetTextFont(42);
  //texSystem2->SetTextSize(0.05);
  texSystem4->Draw("SAME");
  errleg[DrawParticularBin-lowbin]->Draw("SAME");
  
  
  TFile * outFile = new TFile(OutputFileName->Data(), "RECREATE");
  BSpectrum->Write();
  BSpectrumSys->Write();
  DSpectrum->Write();
  DSpectrumSys->Write();
  DBSpectrum->Write();
  DBSpectrumSys->Write();
  ConversionSpectrum->Write();
  ConversionSpectrumSys->Write();
  DalitzSpectrum->Write();
  DalitzSpectrumSys->Write();
  BSpectrumG->Write();
  BSpectrumGsys->Write();
  DSpectrumG->Write();
  DSpectrumGsys->Write();
  DalitzSpectrumG->Write();
  ConversionSpectrumG->Write();
  totalSpectrum->Write();
  evts->Write();
  evtsb0->Write();
  for(int i=lowbin;i<=upbin;i++)
  {
    eleData1d2[i-lowbin]->Write();
    eleConversion1d[i-lowbin]->Write();
    eleDalitz1d[i-lowbin]->Write();
    eleD1d[i-lowbin]->Write();
    eleB1d[i-lowbin]->Write();
    relErr[i-lowbin]->Write();
    relErr[i-lowbin]->Write();
    reconstructed[i-lowbin]->Write();
  }
  outFile->Close();
  
  
  if(ExtraErrorTestNeeded && !ExtraErrorTestSuccessful)
    cout << "Warning! There was a problem with the additional error test!" << endl;

  

}
