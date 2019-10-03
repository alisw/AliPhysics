#include <TAxis.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TMatrixDSym.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <AliLog.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <iostream>
#include <cstring>

#include <AliESDpid.h>

#include "TROOT.h"

#include "AliHFEhadronicbackground.h"


ClassImp(AliHFEhadronicbackground)  

AliHFEhadronicbackground::AliHFEhadronicbackground(char * run):
TObject(),
  fLcutTPC(NULL),
  fUcutTPC(NULL)
{
  GetParametersFromFit=false;
  nSlices=0;
  ranCalculation=false;
  fElectronWidthIsFixed=false;
  fElectronCenterIsFixed=false;
  fLcutTPC=new TF1("lcut","-1.",-10,10);
  fUcutTPC=new TF1("ucut","3.0",-10,10);
  fOutputFileName=new char[300];
  strcpy(fOutputFileName,"output.root");
  
  TFile * functionFile = new TFile("extractedFunctionsAug22.root","READ");
  PiFitHist = (TH1D *)(functionFile->Get("pions")); // must be normalized
  PiFitHist->Scale(1.0/PiFitHist->GetMaximum());
  PiHistCenter=(PiFitHist->GetXaxis())->GetBinCenter(PiFitHist->GetMaximumBin());
  printf("PiHistCenter: %f\n",PiHistCenter);
  PiFitHist->SetDirectory(0);
  functionFile->Close();
  
  useGaussians=0;
  useLandau=0;
  useLandauForElectrons=true;
  electronExponentialParameter=5.5;
  WriteToFile=false;
  FitRangeLow = -12.0;
  FitRangeUp = 5.0;
  
  
  PIDobject=SetupPIDSplinesTPC(run,2,kFALSE,kFALSE,kFALSE);
  
  SetupElectronLandau();
}

AliHFEhadronicbackground::~AliHFEhadronicbackground(){
  printf("destructor.\n");
  if(ranCalculation)
    {
      delete fCenterPi;
      delete fCenterKa;
      delete fCenterEl;
      delete fWidthPi;
      delete fWidthKa;
      delete fWidthEl;
      delete fAmpliPi;
      delete fAmpliEl;
      delete fAmpliKa;
      delete fContamination;
      delete fEfficiency;
      delete fTRDEfficiency;
      delete[] slices;
      delete[] difference;
      delete[] ratio;
      delete[] fElFits;
      delete[] fPiFits;
      delete[] fKaFits;
      delete[] fAllCombined;
      delete[] electronWidth;
      delete[] electronCenter;
      delete[] PionWidth;
      delete[] PionCenter;
      delete[] fAllCombined;
      delete[] momentumArray;
    }
  delete[] fOutputFileName;
}




Double_t AliHFEhadronicbackground::LandauExp(Double_t * x, Double_t * par)
{
  double val=0, dx=0;
  int nrConv= 1;
  for(int i=-nrConv;i<(nrConv+1);i++)
    {
      dx=par[4]*double(i)*3.0/double(nrConv);
      val+=par[0]*TMath::Landau(x[0]+dx, par[1], par[2], true)*TMath::Exp(-par[3]*(x[0]+dx+8.5-par[1]))*TMath::Gaus(dx,0,par[4]);
    } // added par[1] in Exp. Remove if it does not work
  return val;
}


void AliHFEhadronicbackground::SetupElectronLandau(void)
{
  eleLanParas = new Double_t[3];
  eleLanParas[0]=1e20;
  for(int i=1;i<3;i++)eleLanParas[i]=0;
  TF1 * eleLandTest = new TF1("eleLandTest", this, &AliHFEhadronicbackground::ElectronLandauExp, -20.0,20.0,5);
  eleLandTest->SetParameter(0,1);
  eleLandTest->SetParameter(1,0);
  eleLandTest->SetParameter(2,1);
  eleLandTest->SetParameter(3, electronExponentialParameter);
  eleLandTest->SetParameter(4, 0.01);
  
  int n=3;
  eleLanParas[2]=0;
  eleLanParas[1]=10;
  eleLanParas[1]-=eleLandTest->GetMaximumX();
  //cout << "Beginning width: " << eleLandTest->Variance(-20,20) << endl;
  
  while(TMath::Abs(TMath::Sqrt(eleLandTest->Variance(-20,20))-1)>0.02)
  {
    eleLanParas[1]-=eleLandTest->GetMaximumX();
    eleLanParas[0]/=eleLandTest->GetMaximum(-20,20);
     eleLanParas[2]+=(1-TMath::Sqrt(eleLandTest->Variance(-20,20)))/double(n);
     //n++;
     //cout << "eleLanParas[2]: " << eleLanParas[2] << " width: " << TMath::Sqrt(eleLandTest->Variance(-20,20)) << endl;
  }
  eleLanParas[1]-=eleLandTest->GetMaximumX();
  //cout << "Maximum " << eleLandTest->GetMaximum(-5,5) << " at " << eleLandTest->GetMaximumX() << endl;
  eleLanParas[0]/=eleLandTest->GetMaximum(-5,5);
  
}

Double_t AliHFEhadronicbackground::ElectronLandauExp(Double_t * x, Double_t * par) // this is for electrons to have a well defined center
{
  double val=0, dx=0;
  int nrConv= 1;
  for(int i=-nrConv;i<(nrConv+1);i++)
    {
      dx=par[4]*double(i)*3.0/double(nrConv);
      val+=eleLanParas[0]*par[0]*TMath::Landau(x[0]+dx, par[1]+eleLanParas[1], par[2]+eleLanParas[2]+1, true)*TMath::Exp(-par[3]*(x[0]/par[2]-par[1]+dx))*TMath::Gaus(dx,0,par[4]);
    } // added par[1] in Exp. Remove if it does not work
  return val;
}


Double_t  AliHFEhadronicbackground::KaonFitFunction(Double_t * dEdx2, Double_t *par)   // This should not actually be necessary. It is identical to the PionFitFunction
{
  Double_t dEdx=(dEdx2[0]-par[1])/par[2]+PiHistCenter; // For Comparison to Gaussian center on 0
  //  dEdx=dEdx2[0];
  Double_t par0=par[0];
  Int_t binNr = PiFitHist->FindBin(dEdx);
  if(binNr < 2 || binNr>PiFitHist->GetNbinsX()) return 0.001;
  Double_t binCenter = PiFitHist->GetBinCenter(binNr);
  if(binCenter>dEdx) binNr--; // Now x1 gives the lower of the two
  Double_t x2=PiFitHist->GetBinCenter(binNr+1);
  Double_t x1=PiFitHist->GetBinCenter(binNr);
  Double_t y2=PiFitHist->GetBinContent(binNr+1);
  Double_t y1=PiFitHist->GetBinContent(binNr);
  return (y2+(x2-dEdx)/(x2-x1)*(y1-y2))*par0+0.001;
  
}




Double_t  AliHFEhadronicbackground::TotalFitFunction(Double_t * dEdx, Double_t *par)
{
  Double_t electrons, kaons, landauPions;
  
  Double_t ModLandauParam[5];
  ModLandauParam[0]=par[0];
  ModLandauParam[1]=par[1];
  ModLandauParam[2]=par[2];
  ModLandauParam[3]=par[11];
  ModLandauParam[4]=par[12];
  
  Double_t kaonparam[3];
  kaonparam[0]=par[3];
  kaonparam[1]=par[4];
  kaonparam[2]=par[5];
  
  Double_t electronparam[5];
  electronparam[0]=par[6];
  electronparam[1]=par[7];
  electronparam[2]=par[8];
  electronparam[3]=electronExponentialParameter;
  electronparam[4]=0.01;
  
  
  landauPions=LandauExp(dEdx,ModLandauParam);
  if(useLandauForElectrons)
    electrons=ElectronLandauExp(dEdx, electronparam);
    else
      electrons=electronparam[0]*TMath::Gaus(dEdx[0], electronparam[1], electronparam[2]);
  
  kaons = KaonFitFunction(dEdx, kaonparam);
  
  return TMath::Max((electrons + kaons + landauPions)+0.00000001,double(0.0));
  
  // Meaning of Parameters:
  // 0:  Amplitude of Pions
  // 1:  Center of Pions
  // 2:  Width of Pions (deviation from expected width)
  // 3:  Amplitude of Kaons
  // 4:  Center of Kaons (deviation from expected center)
  // 5:  Width of Kaons (deviation from expected width)
  // 6:  Amplitude of Electrons
  // 7:  Center of Electrons
  // 8:  Width of Electrons
  // 9:  Epsilon parameter for Pions
  // 10: Amplitude of Muons
  // 11: Center of Muons (deviation from expected center)
  // 12: Width of Muons
}



void AliHFEhadronicbackground::SetElectronWidth(Double_t Width){fFixedElectronWidth=Width;fElectronWidthIsFixed=true;}

void AliHFEhadronicbackground::SetElectronCenter(Double_t Center){fFixedElectronCenter=Center;fElectronCenterIsFixed=true;}

void AliHFEhadronicbackground::SetOutputFileName(const char * newFileName){strcpy(fOutputFileName,newFileName);}


void AliHFEhadronicbackground::Process(TH2D* hSigmaTPC, Double_t pMin, Double_t pMax){
  if(pMin > pMax){
    // order swapped, swap back
    Double_t tmp = pMin;
    pMin = pMax;
    pMax = tmp;
  }
  
  if(ranCalculation)
    {
      for(int i=0;i<nSlices;i++)
	{
	  delete slices[i];
	}
      nSlices=0;
      delete fLambdaPion;
      delete fCenterPi;
      delete fCenterKa;
      delete fCenterEl;
      delete fWidthPi;
      delete fWidthKa;
      delete fWidthEl;
      delete fAmpliPi;
      delete fAmpliEl;
      delete fAmpliKa;
      delete fContamination;
      delete fEfficiency;
      delete fTRDEfficiency;
      delete[] slices;
      delete[] difference;
      delete[] ratio;
      delete[] fElFits;
      delete[] fPiFits;
      delete[] fKaFits;
      delete[] fAllCombined;
      delete[] electronWidth;
      delete[] electronCenter;
      delete[] PionWidth;
      delete[] PionCenter;
      delete[] fAllCombined;
      delete[] momentumArray;
    }
  
  Int_t firstBin = (hSigmaTPC->GetXaxis())->FindBin(pMin);
  Int_t lastBin = (hSigmaTPC->GetXaxis())->FindBin(pMax);
  Int_t Ybins = hSigmaTPC->GetNbinsY();
  UInt_t nSteps=UInt_t(lastBin-firstBin+1);
  Double_t stepWidth = ((hSigmaTPC->GetXaxis())->GetBinUpEdge(lastBin)-(hSigmaTPC->GetXaxis())->GetBinLowEdge(firstBin))/Double_t(nSteps);
  
  momentumArray=new double[nSteps];
  
  slices=new TH1D*[nSteps];
  difference=new TH1D*[nSteps];
  ratio=new TH1D*[nSteps];
  fElFits=new TF1*[nSteps];
  fPiFits=new TF1*[nSteps];
  fKaFits=new TF1*[nSteps];
  fAllCombined=new TF1*[nSteps];
  
  FitLegend=new TLegend*[nSteps];
  
  electronWidth=new Double_t[nSteps];
  electronCenter=new Double_t[nSteps];
  PionWidth=new Double_t[nSteps];
  PionCenter=new Double_t[nSteps];
  nParameters=nSteps; // this might not work after changing
  parMin=(hSigmaTPC->GetXaxis())->GetBinLowEdge(firstBin)+stepWidth/2.;
  parStep=stepWidth;
  
  Int_t pBin = 0;
  
  double_t tempVal,tempVal2;
  
  for(UInt_t iStep = 0; iStep < nSteps; iStep++){
    pBin = firstBin+iStep;
    momentumArray[iStep]=(hSigmaTPC->GetXaxis())->GetBinCenter(pBin);
    slices[iStep] = hSigmaTPC->ProjectionY(Form("Slice%d", iStep), pBin, pBin);
    nSlices++;
    FitSlice(slices[iStep], momentumArray[iStep]); // Here the fit happens
    fAllCombined[iStep]=tempCombined;
    if(useLandauForElectrons)
      {
    fElFits[iStep]= new TF1("eleLandTest", this, &AliHFEhadronicbackground::ElectronLandauExp, -20.0,20.0,5);
    fElFits[iStep]->SetParameter(0,tempCombined->GetParameter(6));
    fElFits[iStep]->SetParameter(1,tempCombined->GetParameter(7));
    fElFits[iStep]->SetParameter(2,tempCombined->GetParameter(8));  
    fElFits[iStep]->SetParameter(3,electronExponentialParameter);
    fElFits[iStep]->SetParameter(4,0.01);      
      }
      else
      {
    fElFits[iStep]=new TF1(Form("el%d",iStep), "[0] * TMath::Gaus(x, [1], [2])", -20, 20);
    fElFits[iStep]->SetParameter(0,tempCombined->GetParameter(6));
    fElFits[iStep]->SetParameter(1,tempCombined->GetParameter(7));
    fElFits[iStep]->SetParameter(2,tempCombined->GetParameter(8));
      }
    fPiFits[iStep]=new TF1(Form("pi%d",iStep), this, &AliHFEhadronicbackground::LandauExp, -20.0,20.0,5);
    fPiFits[iStep]->SetParameter(0,tempCombined->GetParameter(0));
    fPiFits[iStep]->SetParameter(1,tempCombined->GetParameter(1));
    fPiFits[iStep]->SetParameter(2,tempCombined->GetParameter(2));
    fPiFits[iStep]->SetParameter(3,tempCombined->GetParameter(11));
    fPiFits[iStep]->SetParameter(4,tempCombined->GetParameter(12));
    fKaFits[iStep] = new TF1(Form("ka%d",iStep), this, &AliHFEhadronicbackground::KaonFitFunction, -20.0,20.0,3);
    fKaFits[iStep]->SetParameter(0,tempCombined->GetParameter(3));
    fKaFits[iStep]->SetParameter(1,tempCombined->GetParameter(4));
    fKaFits[iStep]->SetParameter(2,tempCombined->GetParameter(5));

    difference[iStep]=(TH1D*)slices[iStep]->Clone();
    ratio[iStep]=(TH1D*)slices[iStep]->Clone();
    for(int bin=0;bin<Ybins-1;bin++)
        {
	  tempVal=TMath::Abs(difference[iStep]->GetBinContent(bin)-fAllCombined[iStep]->Eval(difference[iStep]->GetBinCenter(bin)))
	    /TMath::Sqrt(0.0001+fAllCombined[iStep]->Eval(difference[iStep]->GetBinCenter(bin)));
	  tempVal2=difference[iStep]->GetBinContent(bin)/(0.00001+fAllCombined[iStep]->Eval(difference[iStep]->GetBinCenter(bin)));
          difference[iStep]->Fill(difference[iStep]->GetBinCenter(bin),-difference[iStep]->GetBinContent(bin)+tempVal);
          ratio[iStep]->Fill(ratio[iStep]->GetBinCenter(bin),-ratio[iStep]->GetBinContent(bin)+tempVal2);
        }    
  }
  
  fLambdaPion= new TH1D("LambdaPar", "lambdaPar",nParameters,parMin-parStep/2.,nParameters*parStep+parMin-parStep/2.);
  
  
  fCenterPi = new TH1D("centerpi", "centerpi",nParameters,parMin-parStep/2.,nParameters*parStep+parMin-parStep/2.);
  fCenterKa = new TH1D("centerka", "centerka",nParameters,parMin-parStep/2.,nParameters*parStep+parMin-parStep/2.);
  fCenterEl = new TH1D("centerel", "centerel",nParameters,parMin-parStep/2.,nParameters*parStep+parMin-parStep/2.);
  
  for(int i=0;i<nParameters;i++)
    {
      fCenterPi->Fill(parMin+parStep*double(i),fAllCombined[i]->GetParameter(1));
      fCenterKa->Fill(parMin+parStep*double(i),0);
      fCenterEl->Fill(parMin+parStep*double(i),fAllCombined[i]->GetParameter(7));    
      fLambdaPion->Fill(parMin+parStep*double(i),fAllCombined[i]->GetParameter(11));
    }
  
  fWidthPi = new TH1D("widthpi", "widthpi",nParameters,parMin-parStep/2.,nParameters*parStep+parMin-parStep/2.);
  fWidthKa = new TH1D("widthka", "widthka",nParameters,parMin-parStep/2.,nParameters*parStep+parMin-parStep/2.);
  fWidthEl = new TH1D("widthel", "widthel",nParameters,parMin-parStep/2.,nParameters*parStep+parMin-parStep/2.);
  for(int i=0;i<nParameters;i++)
    {
      fWidthPi->Fill(parMin+parStep*double(i),fAllCombined[i]->GetParameter(2));
      fWidthKa->Fill(parMin+parStep*double(i),fAllCombined[i]->GetParameter(5));
      fWidthEl->Fill(parMin+parStep*double(i),fAllCombined[i]->GetParameter(8));
    }
  
  fAmpliPi = new TH1D("amplipi", "amplipi",nParameters,parMin-parStep/2.,nParameters*parStep+parMin-parStep/2.);
  fAmpliKa = new TH1D("amplika", "amplika",nParameters,parMin-parStep/2.,nParameters*parStep+parMin-parStep/2.);
  fAmpliEl = new TH1D("ampliel", "ampliel",nParameters,parMin-parStep/2.,nParameters*parStep+parMin-parStep/2.);
  for(int i=0;i<nParameters;i++)
    {
      fAmpliPi->Fill(parMin+parStep*double(i),fAllCombined[i]->GetParameter(0));
      fAmpliKa->Fill(parMin+parStep*double(i),fAllCombined[i]->GetParameter(3));
      fAmpliEl->Fill(parMin+parStep*double(i),fAllCombined[i]->GetParameter(6));
    }
  
  fContamination = new TH1D("Contamination", "contamination",nParameters ,parMin-parStep/2.,parMin+parStep*nParameters-parStep/2.);
  fEfficiency = new TH1D("Complementary Efficiency", "compefficiency",nParameters ,parMin-parStep/2.,parMin+parStep*nParameters-parStep/2.);
  fTRDEfficiency = new TH1D("TRD Efficiency", "efficiencyTRD",nParameters ,parMin-parStep/2.,parMin+parStep*nParameters-parStep/2.);
  double elnr, elAfterCut, piAfterCut;
  elNrDiag = new TH1D("elnr", "elnr", nParameters,parMin-parStep/2.,parMin+parStep*nParameters-parStep/2.);
  
  for(unsigned int i=0;i<nSteps;i++)
    {
      elnr=fElFits[i]->Integral(-10.,5.)/slices[i]->GetBinWidth(5);
      elNrDiag->Fill(parMin+i*parStep,elnr);
      elAfterCut=fElFits[i]->Integral(fLcutTPC->Eval(parMin+i*parStep),fUcutTPC->Eval(parMin+i*parStep))/slices[i]->GetBinWidth(5);
      piAfterCut=fAllCombined[i]->Integral(fLcutTPC->Eval(parMin+i*parStep),fUcutTPC->Eval(parMin+i*parStep))/slices[i]->GetBinWidth(5)-elAfterCut;
      if(elAfterCut<0.1)elAfterCut=0.1;
      if(piAfterCut<0.1)piAfterCut=0.1;
      fContamination->Fill(parMin+i*parStep,piAfterCut/(piAfterCut+elAfterCut));
      fContamination->SetBinError(i+1, 1.0/TMath::Sqrt(piAfterCut)*piAfterCut/(piAfterCut+elAfterCut));
      printf("Nr el at %f GeV: %f\n",parMin+i*parStep,elnr);
      fCenterEl->SetBinError(i+1, 2.002/TMath::Sqrt(elnr));
      fEfficiency->Fill(parMin+i*parStep,1.-elAfterCut/elnr);
    }
  
  
  // Output start
  
  char title[30];
  if(WriteToFile)
    {
      TFile * outFile = new TFile(fOutputFileName, "RECREATE");
      
      for(unsigned int i=0;i<nSteps;i++)
	{
	  sprintf(title, "%iMeVto%iMeVFits", int(1000*(parMin+(i-0.5)*parStep)),int(1000*(parMin+(0.5+i)*parStep)));
	  fAllCombined[i]->SetName(title);
	  fAllCombined[i]->Write();
	  sprintf(title, "%iMeVto%iMeVData", int(1000*(parMin+(i-0.5)*parStep)),int(1000*(parMin+(0.5+i)*parStep)));
	  slices[i]->SetName(title);
	  slices[i]->Write();
	  sprintf(title, "%iMeVto%iMeVElectrons", int(1000*(parMin+(i-0.5)*parStep)),int(1000*(parMin+(0.5+i)*parStep)));
	  fElFits[i]->SetName(title);
	  fElFits[i]->Write();
	  sprintf(title, "%iMeVto%iMeVPions", int(1000*(parMin+(i-0.5)*parStep)),int(1000*(parMin+(0.5+i)*parStep)));
	  fPiFits[i]->SetName(title);
	  fPiFits[i]->Write();
	  sprintf(title, "%iMeVto%iMeVKaons", int(1000*(parMin+(i-0.5)*parStep)),int(1000*(parMin+(0.5+i)*parStep)));
	  fKaFits[i]->SetName(title);
	  fKaFits[i]->Write();
	  
	  sprintf(title, "%iMeVto%iMeVRatio", int(1000*(parMin+(i-0.5)*parStep)),int(1000*(parMin+(0.5+i)*parStep)));
	  ratio[i]->SetName(title);
	  ratio[i]->Write();
	  sprintf(title, "%iMeVto%iMeVrelDiff", int(1000*(parMin+(i-0.5)*parStep)),int(1000*(parMin+(0.5+i)*parStep)));
	  difference[i]->SetName(title);
	  difference[i]->Write();
	}
      
      sprintf(title, "fCenterElectrons%iMeVto%iMeV", int(1000*(parMin-0.0*parStep)),int(1000*(parMin+(-0.0+nParameters)*parStep)));
      fCenterEl->SetName(title);
      fCenterEl->Write();
      sprintf(title, "fCenterPions%iMeVto%iMeV", int(1000*(parMin-0.0*parStep)),int(1000*(parMin+(-0.0+nParameters)*parStep)));
      fCenterPi->SetName(title);
      fCenterPi->Write();
      sprintf(title, "fCenterKaons%iMeVto%iMeV", int(1000*(parMin-0.0*parStep)),int(1000*(parMin+(-0.0+nParameters)*parStep)));
      fCenterKa->SetName(title);
      fCenterKa->Write();
      
      sprintf(title, "fWidthElectrons%iMeVto%iMeV", int(1000*(parMin-0.0*parStep)),int(1000*(parMin+(-0.0+nParameters)*parStep)));
      fWidthEl->SetName(title);
      fWidthEl->Write();
      sprintf(title, "fWidthPions%iMeVto%iMeV", int(1000*(parMin-0.0*parStep)),int(1000*(parMin+(-0.0+nParameters)*parStep)));
      fWidthPi->SetName(title);
      fWidthPi->Write();
      sprintf(title, "fWidthKaons%iMeVto%iMeV", int(1000*(parMin-0.0*parStep)),int(1000*(parMin+(-0.0+nParameters)*parStep)));
      fWidthKa->SetName(title);
      fWidthKa->Write();
      
      sprintf(title, "AmplitudeElectrons%iMeVto%iMeV", int(1000*(parMin-0.0*parStep)),int(1000*(parMin+(-0.0+nParameters)*parStep)));
      fAmpliEl->SetName(title);
      fAmpliEl->Write();
      sprintf(title, "AmplitudePions%iMeVto%iMeV", int(1000*(parMin-0.0*parStep)),int(1000*(parMin+(-0.0+nParameters)*parStep)));
      fAmpliPi->SetName(title);
      fAmpliPi->Write();
      sprintf(title, "AmplitudeKaons%iMeVto%iMeV", int(1000*(parMin-0.0*parStep)),int(1000*(parMin+(-0.0+nParameters)*parStep)));
      fAmpliKa->SetName(title);
      fAmpliKa->Write();
      
      sprintf(title, "TPCCutEfficiency%iMeVto%iMeV", int(1000*(parMin-0.0*parStep)),int(1000*(parMin+(-0.0+nParameters)*parStep)));
      fEfficiency->SetName(title);
      fEfficiency->Write();
      
      
      sprintf(title, "TPCCUTContamination%iMeVto%iMeV", int(1000*(parMin-0.5*parStep)),int(1000*(parMin+(-0.5+nParameters)*parStep)));
      fContamination->SetName(title);
      fContamination->Write();
      
      outFile->Close();
    }
  ranCalculation=true;
  
  
  printf("End of Process()\n");
}




void AliHFEhadronicbackground::FitSlice(TH1 *hslice, Double_t momentum){
  // Make fit for electrons and pions (hSigmaTPC->GetXaxis())->GetBinLowEdge(pBin)separately first to constrain
  // parameters for the combined fit
  //
  if(!fLcutTPC){
    AliError("Cannot fit electrons without lower cut");
    return;
  }
  
  TF1 *combined;
  
  // Prepare use of Splines:
  AliTPCPIDResponse* TPCResponse= &PIDobject->GetTPCResponse();
  
  combined = new TF1("combined", this, &AliHFEhadronicbackground::TotalFitFunction, -20.0,20.0,13);
  combined->SetNpx(10000);
  
  combined->SetParName(0,"pi ampl.");
  combined->SetParName(1,"pi cent.");
  combined->SetParName(2,"pi widt. (r)");
  combined->SetParName(3,"ka ampl.");
  combined->SetParName(4,"ka cent. (r)");
  combined->SetParName(5,"ka widt. (r)");
  combined->SetParName(6,"el ampl.");
  combined->SetParName(7,"el cent.");
  combined->SetParName(8,"el widt.");
  combined->SetParName(9,"momentum");
  combined->SetParName(10,"mu ampl.");
  combined->SetParName(11,"pi lambda");
  combined->SetParName(12,"pi widening");
  
  // Meaning of Parameters:
  // 0:  Amplitude of Pions
  // 1:  Center of Pions
  // 2:  Width of Pions (deviation from expected width)
  // 3:  Amplitude of Kaons
  // 4:  Center of Kaons (deviation from expected center)
  // 5:  Width of Kaons (deviation from expected width)
  // 6:  Amplitude of Electrons
  // 7:  Center of Electrons
  // 8:  Width of Electrons
  // 9:  Momentum
  // 10: Amplitude of Muons
  // 11: Pion Lambda
  // 12: Pion Gaussian widening
  
  // some initial constraints to make the fit converging
  Double_t minEle = 5e-1;
  Double_t maxEle = 1e9;
  Double_t maxPi = 2*hslice->GetMaximum();
  
  printf("\n\nNow fitting slice at momentum %f GeV\n",momentum);
  
  combined->FixParameter(9, momentum);
  combined->SetParLimits(6, minEle, maxEle);
  combined->SetParameter(6, 1.0);  // electron stuff
  combined->SetParLimits(8, 0.5, 1.3);
  combined->SetParameter(8, 1.);
  combined->SetParLimits(7, -0.7, 0.9);
  if(fElectronWidthIsFixed)  combined->FixParameter(8,fFixedElectronWidth);
  if(fElectronCenterIsFixed) combined->FixParameter(7,fFixedElectronCenter);
  
  // sort here all fixed
  
  combined->FixParameter(11, 0);
  combined->FixParameter(10,0.0);
  // end of fixed
  
  
  combined->SetParLimits(3, 1., maxPi);
  combined->SetParLimits(4, -6.5, 4.5);
  combined->SetParLimits(5, 0.5,1.1);
  combined->SetParameter(4, -5.0);
  combined->SetParameter(5, 0.8);  // Kaon Stuff
  
  combined->SetParLimits(0, 0.0, 1e8);
  combined->SetParameter(0, 200.0);
  combined->SetParLimits(1, 3.0, 9.5); 
  combined->SetParameter(1, 6.0); 
  combined->SetParLimits(2, 0.3, 4.0);  // Pion Stuff
  combined->SetParameter(2, 2.4);
  combined->SetParLimits(11, 5.0, 6.3); 
  combined->SetParameter(11, 5.5);
  combined->FixParameter(12,0.002);
  
  TF1 * currentPiFit=new TF1("temppi", this, &AliHFEhadronicbackground::LandauExp, -20.0,20.0,5);
  currentPiFit->SetParameter(0,combined->GetParameter(0));
  currentPiFit->SetParameter(1,combined->GetParameter(1));
  currentPiFit->SetParameter(2,combined->GetParameter(2));
  currentPiFit->SetParameter(3,combined->GetParameter(11));
  currentPiFit->SetParameter(4,combined->GetParameter(12));
  combined->SetParameter(1, combined->GetParameter(1)+TPCResponse->GetNumberOfSigmas(momentum, TPCResponse->GetExpectedSignal(momentum,AliPID::kPion), 125, AliPID::kElectron)
				-currentPiFit->GetMaximumX());
  double kaPos=combined->GetParameter(4)+TPCResponse->GetNumberOfSigmas(momentum, TPCResponse->GetExpectedSignal(momentum,AliPID::kKaon), 125, AliPID::kElectron)
				-currentPiFit->GetMaximumX();
  combined->SetParameter(4,kaPos);
  combined->SetParLimits(4,kaPos-2.3,kaPos+2.3);
  delete currentPiFit;
  
  //combined->FixParameter(0,0.00);
  TFitResultPtr fitresultsCombined = hslice->Fit(combined, "NLM", "", FitRangeLow, FitRangeUp); // for pp 3.0 with "N"
  
  tempCombined=combined;
  
  printf("end of FitSlice\n");
}

void AliHFEhadronicbackground::DrawDiagrams(void)
{
  printf("Start Draw Diagram Function\n");
  char name[30];
  TCanvas ** IndividualSlice=new TCanvas*[nParameters/4+1];
  for(int i=0;i<nParameters/4+1;i++)
    {
      sprintf(name,"cv%i", i);
      IndividualSlice[i]=new TCanvas(name, "Fits",1200, 900);
      IndividualSlice[i]->SetFillColor(10);
      IndividualSlice[i]->Divide(2,2);
      for(int j=0;j<4 && j<nParameters-i*4;j++){
	FitLegend[i]=new TLegend(0.6,0.65,1.,1.);
        IndividualSlice[i]->cd(j+1);
        gPad->SetLogy(1);
        gStyle->SetOptTitle(0);
        gStyle->SetCanvasColor(0);
        gStyle->SetOptStat(0);
        slices[i*4+j]->SetAxisRange(-10, 4, "X");
        slices[i*4+j]->SetLineColor(1);
        slices[i*4+j]->SetLineWidth(1);
        slices[i*4+j]->Draw();
        difference[i*4+j]->SetLineColor(6);
        difference[i*4+j]->SetLineWidth(1);
        difference[i*4+j]->Draw("SAME");
        ratio[i*4+j]->SetLineColor(7);
        ratio[i*4+j]->SetLineWidth(1);
        ratio[i*4+j]->Draw("SAME");
        fPiFits[i*4+j]->SetLineColor(3);
        fPiFits[i*4+j]->SetLineWidth(1);
        fPiFits[i*4+j]->Draw("SAME");
        fKaFits[i*4+j]->SetLineColor(15);
        fKaFits[i*4+j]->SetLineWidth(1);
        fKaFits[i*4+j]->Draw("SAME");  
        fElFits[i*4+j]->SetLineColor(2);
        fElFits[i*4+j]->SetLineWidth(1);
        fElFits[i*4+j]->Draw("SAME");
        fAllCombined[i*4+j]->SetLineColor(50);
        fAllCombined[i*4+j]->SetLineWidth(1);
        fAllCombined[i*4+j]->Draw("SAME");
        FitLegend[i]->AddEntry(slices[i],"Data");
        FitLegend[i]->AddEntry(fAllCombined[i],"Combined Fit");
        FitLegend[i]->AddEntry(fElFits[i],"Electron Fit");
        FitLegend[i]->AddEntry(fPiFits[i],"Pion Fit");
        FitLegend[i]->AddEntry(fKaFits[i],"Kaon Fit");
        FitLegend[i]->AddEntry(ratio[i],"Ratio Data/Fit");
        FitLegend[i]->AddEntry(difference[i],"Deviation relative to expected");
        sprintf(name,"%.1f GeV to %.1f GeV", double(4*i+j-0.5)*parStep+parMin, double(i*4+j+0.5)*parStep+parMin);
        FitLegend[i]->SetHeader(name);
        FitLegend[i]->SetFillColor(0);
        FitLegend[i]->SetLineStyle(1);
        FitLegend[i]->SetLineWidth(1);
        FitLegend[i]->SetBorderSize(1);
        FitLegend[i]->Draw();
	
       }
    }
  
  TCanvas * analysisCanvas=new TCanvas("cvana", "consistency",1200, 900);
  analysisCanvas->Divide(2,2);
  analysisCanvas->SetFillColor(10);
  analysisCanvas->cd(1);
  fCenterPi->SetAxisRange(-10, 4, "Y");
  gStyle->SetCanvasColor(0);
  (fCenterPi->GetXaxis())->SetTitle("p [GeV/c]");
  (fCenterPi->GetYaxis())->SetTitle("(Center-<electron peak>)/#sigma_{e}");
  fCenterPi->SetMarkerColor(3);
  fCenterPi->SetMarkerStyle(3);
  fCenterPi->Draw("P");
  fCenterKa->SetMarkerColor(15);
  fCenterKa->SetMarkerStyle(3);
  fCenterKa->Draw("P SAME");
  fCenterEl->SetMarkerColor(2);
  fCenterEl->SetMarkerStyle(3);
  fCenterEl->Draw("P SAME");
  fLambdaPion->SetMarkerColor(kCyan);
  fLambdaPion->SetMarkerStyle(3);
  fLambdaPion->Draw("P SAME");
  TLegend * consistencyLegend1=new TLegend(0.7,0.7,1.,1.);
  consistencyLegend1->SetHeader("Center of Distributions");
  consistencyLegend1->AddEntry(fCenterEl,"Electrons");
  consistencyLegend1->AddEntry(fCenterPi,"Pions");
  consistencyLegend1->AddEntry(fCenterKa,"Kaons");
  consistencyLegend1->AddEntry(fLambdaPion,"Lambda Parameter");
  consistencyLegend1->SetFillColor(0);
  consistencyLegend1->SetLineStyle(1);
  consistencyLegend1->SetLineWidth(1);
  consistencyLegend1->SetBorderSize(1);
  consistencyLegend1->Draw("SAME");
  
  analysisCanvas->cd(2);
  fWidthPi->SetAxisRange(0, 2, "Y");
  (fWidthPi->GetXaxis())->SetTitle("p [GeV/c]");
  (fWidthPi->GetYaxis())->SetTitle("Width/Electron Width");
  fWidthPi->SetMarkerColor(3);
  fWidthPi->SetMarkerStyle(3);
  fWidthPi->Draw("PL");
  fWidthKa->SetMarkerColor(15);
  fWidthKa->SetMarkerStyle(3);
  fWidthKa->Draw("PL SAME");
  fWidthEl->SetMarkerColor(2);
  fWidthEl->SetMarkerStyle(3);
  fWidthEl->Draw("PL SAME");
  TLegend * consistencyLegend2=new TLegend(0.7,0.7,1.,1.);
  consistencyLegend2->SetHeader("Width of Distributions");
  consistencyLegend2->AddEntry(fWidthEl,"Electrons");
  consistencyLegend2->AddEntry(fWidthPi,"Pions");
  consistencyLegend2->AddEntry(fWidthKa,"Kaons");
  consistencyLegend2->SetFillColor(0);
  consistencyLegend2->SetLineStyle(1);
  consistencyLegend2->SetLineWidth(1);
  consistencyLegend2->SetBorderSize(1);
  consistencyLegend2->Draw("SAME");
  
  analysisCanvas->cd(3);
  
  analysisCanvas->cd(4);
  fCenterEl->Draw("P");
  TF1 * constantCenterFit = new TF1 ("centerFitConst", "[0]", 0,10);
  fCenterEl->Fit(constantCenterFit, "N", "", 0,10);
  constantCenterFit->SetLineColor(2);
  constantCenterFit->SetLineWidth(2);
  constantCenterFit->Draw("SAME");
  
  TLegend * consistencyLegend4=new TLegend(0.5,0.7,0.8,0.8);
  consistencyLegend4->AddEntry(fCenterEl,"Center of Electron Distribution");
  consistencyLegend4->AddEntry(constantCenterFit,"Weighted Average");
  consistencyLegend4->SetFillColor(0);
  consistencyLegend4->SetLineStyle(1);
  consistencyLegend4->SetLineWidth(1);
  consistencyLegend4->SetBorderSize(1);
  consistencyLegend4->Draw("SAME");
  
  TCanvas * CutAnalysis=new TCanvas("cutana", "Analysis of Cuts",1200, 600);
  CutAnalysis->Divide(2,1);
  CutAnalysis->SetFillColor(10);
  CutAnalysis->cd(1);
  gPad->SetLogy(1);
  gStyle->SetCanvasColor(0);
  fContamination->SetMarkerStyle(3);
  fContamination->SetAxisRange(fContamination->GetMinimum()/10.,1.2,"Y");
  (fContamination->GetXaxis())->SetTitle("p [GeV/c]");
  (fContamination->GetYaxis())->SetTitle("Contamination of Electrons after TPC Cut");
  fContamination->Draw("PE");
  CutAnalysis->cd(2);
  gPad->SetLogy(1);
  
  fEfficiency->SetMarkerStyle(3);
  (fEfficiency->GetXaxis())->SetTitle("p [GeV/c]");
  (fEfficiency->GetYaxis())->SetTitle("1 - Electron Efficiency of TPC Cut");
  fEfficiency->SetAxisRange(fEfficiency->GetMinimum()/10.,1.2,"Y");
  fEfficiency->Draw("P");
  CutAnalysis->cd(3);
  gPad->SetLogy(1);
  fTRDEfficiency->SetMarkerStyle(3);
  (fTRDEfficiency->GetXaxis())->SetTitle("p [GeV/c]");
  
  /*TCanvas * testCv=new TCanvas("testCv", "Test Canvas",800, 600);
  TF1 * elemodLandau = new TF1("elemodLandau", this, &AliHFEhadronicbackground::ElectronLandauExp, -20.0,20.0,5);
  elemodLandau->SetParameter(0,1);
  elemodLandau->SetParameter(1,0);
  elemodLandau->SetParameter(2,1);
  elemodLandau->SetParameter(3, electronExponentialParameter);
  elemodLandau->SetParameter(4, 0.01);
  elemodLandau->SetRange(-3,5);
  elemodLandau->Draw();
  gPad->SetLogy();
  cout << "Width of test: " << TMath::Sqrt(elemodLandau->Variance(-5,5)) << "\nMean: " <<  elemodLandau->Mean(-5,5)  << endl;//*/
}




AliESDpid* AliHFEhadronicbackground::SetupPIDSplinesTPC(const TString &periodin, Int_t pass, Bool_t isMC, Bool_t draw, Bool_t mipMom)
{
  //
  // Setup an AliESDpid object with the proper splines attached.
  // Setting the period (data) is needed in order to identify the correct run number and beam type.
  //  This is mandatory and also needs to be set for MC (corresponding data period)
  // In addition the reconstruction pass is needed. In case of MC this will be set to 1
  // In case of MC set isMC to kTRUE
  // If draw is kTRUE, the splines will be drawn into the current canvas (same option)
  // By default the MIP is scaled to 1 and the abscissa is betagamma. In case mipMom is kTRUE (default), it is scaled
  //  to the default mip 50 and the particle momentum.
  //
  // !!! NOTE:  If you intend to use the AliESDpid object to calculate nSigmas you MUST set mipMom to kFALSE !!!
  //
  // The returned AliESDpid object contains the splines for each particle specie.
  // See at the end of the macro (draw part) how to access the splines from the object, if needed.
  //
  TString period(periodin);

  // ESD pid response functions
  AliESDpid *fESDpid=0x0;
 
  //set run number and beamtype according to the period
  TString beamtype="PP";
  if(!period.CompareTo("LHC10H") || !period.CompareTo("LHC11H")) beamtype = "PBPB";
 
  //Get CDB Entry with pid response parametrisations
  //TFile *in = TFile::Open("$ALICE_ROOT/OADB/COMMON/PID/data/TPCPIDResponse.root");
  TFile *in = TFile::Open("$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCPIDResponse.root");
  gROOT->cd();
  TObjArray *fArrPidResponseMaster=0x0;
  if (in){
    fArrPidResponseMaster=dynamic_cast<TObjArray*>(in->Get("TPCPIDResponse"));
  }

  //
  //set the new parametrisation
  //
    //data type
  TString datatype="DATA";
  //in case of mc pass is per default 1
  if (isMC) {
    datatype="MC";
    pass=1;
  }
 
  if (fArrPidResponseMaster){
    fESDpid=new AliESDpid;
    TSpline3 *grAll=0x0;
    //for MC don't use period information
    if (isMC) period="[A-Z0-9]*";
    //pattern for the default entry (valid for all particles)
    TPRegexp reg(Form("TSPLINE3_%s_([A-Z]*)_%s_PASS%d_%s_MEAN",datatype.Data(),period.Data(),pass,beamtype.Data()));
   
    //loop over entries and filter them
    for (Int_t iresp=0; iresp<fArrPidResponseMaster->GetEntriesFast();++iresp){
      TSpline3 *responseFunction=(TSpline3*)fArrPidResponseMaster->At(iresp);
     
      if (responseFunction==0x0) continue;
      TString responseName=responseFunction->GetName();
     
      if (!reg.MatchB(responseName)) continue;
     
      TObjArray *arr=reg.MatchS(responseName);
      TString particleName=arr->At(1)->GetName();
      delete arr;
      if (particleName.IsNull()) continue;
      if (particleName=="ALL") grAll=responseFunction;
      else {
        //find particle id
        for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
          TString particle=AliPID::ParticleName(ispec);
          particle.ToUpper();
          if ( particle == particleName ){
            //scale to MIP and P (instead of 1 and beta gamma)
            if (mipMom){
              Double_t mass=AliPID::ParticleMass(ispec);
              for (Int_t iknot=0; iknot<responseFunction->GetNp(); ++iknot){
                Double_t x=0,y=0;
                responseFunction->GetKnot(iknot,x,y);
                responseFunction->SetPoint(iknot,x*mass,50*y);
              }
            }
           
            //test if there is already a function set. If yes, cleanup
            TObject *old=const_cast<TObject*>(fESDpid->GetTPCResponse().GetResponseFunction((AliPID::EParticleType)ispec));
            if (old) delete old;
            fESDpid->GetTPCResponse().SetResponseFunction((AliPID::EParticleType)ispec,responseFunction);
            fESDpid->GetTPCResponse().SetUseDatabase(kTRUE);
            printf("SetupPIDSplinesTPC Adding graph: %d - %s\n",ispec,(const char*)responseFunction->GetName());
            break;
          }
        }
      }
    }
   
    //set default response function to all particles which don't have a specific one
    if (grAll){
      for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
        TSpline3 *responseFunction=(TSpline3*)grAll;
        if (!fESDpid->GetTPCResponse().GetResponseFunction((AliPID::EParticleType)ispec)){
          responseFunction=(TSpline3*)responseFunction->Clone(AliPID::ParticleName(ispec));
          //scale to MIP and P (instead of 1 and beta gamma)
          if (mipMom){
            Double_t mass=AliPID::ParticleMass(ispec);
            for (Int_t iknot=0; iknot<responseFunction->GetNp(); ++iknot){
              Double_t x=0,y=0;
              responseFunction->GetKnot(iknot,x,y);
              responseFunction->SetPoint(iknot,x*mass,50*y);
            }
          }
         
          fESDpid->GetTPCResponse().SetResponseFunction((AliPID::EParticleType)ispec,responseFunction);
          printf("SetupPIDSplinesTPC Adding graph: %d - %s\n",ispec,(const char *)grAll->GetName());
        }
      }
    }
   
  } else {
    Error("SetupPIDSplinesTPC","Spline Array not set!!!");
  }

  //
  // if requested draw the pid response functions in the current canvas
  //
  if (fESDpid&&draw){
    for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
      TSpline3 *spline=0x0;
      if ( (spline=(TSpline3*)fESDpid->GetTPCResponse().GetResponseFunction((AliPID::EParticleType)ispec)) ){
        spline->SetNpx(500);
        spline->Draw("csame");
      }
    }
  }
 
  return fESDpid;
}



