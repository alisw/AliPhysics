// Code by Martin Völkl:  martin.andreas.volkl@cern.ch

#include <AliMCLogLFitter.h>
#include <TH1.h>
#include <TH2D.h>
#include <TMath.h>
#include <TMinuit.h>
#include <iostream>

ClassImp(AliMCLogLFitter);

AliMCLogLFitter::AliMCLogLFitter(Int_t nMCFunctions, TH1 ** MCDistributions, TH1 * dataHistogram)
{
  // Constructor
  // Sets some constants and makes a guess for the initial values.
  fNPar=nMCFunctions;
  fParameter=new Double_t[fNPar];
  fBestParameters=new Double_t[fNPar];
  fFitDiagsMC=MCDistributions;
  fFitDiagsMCExpected=new TH1*[fNPar];
  fFitDiagData=dataHistogram;
  fStartParameter=new Double_t[fNPar];
  fUpperLimitParameter=new Double_t[fNPar];
  fLowerLimitParameter=new Double_t[fNPar];
  for(int i=0;i<fNPar;i++)
  {
    fStartParameter[i]=1.0/double(nMCFunctions); // was / 10.0
    fUpperLimitParameter[i]=1.0;//1.0;
    fLowerLimitParameter[i]=0.001;  
    fFitDiagsMCExpected[i]= (TH1*)fFitDiagsMC[i]->Clone();
    fFitDiagsMCExpected[i]->SetName(Form("DiagramExpectedNr%d",i));
  }  
  int totalBins=dataHistogram->GetNbinsX();
  for(int i=0;i<fNPar;i++)if(totalBins!=fFitDiagsMC[i]->GetNbinsX())  std::cout << "Careful! Number of Bins does not match between the different distributions." << endl;
  
  flowerFitRange=1;
  fupperFitRange=totalBins;
  
  // Standard constants
  fiterations=10; // was 3
  fDefaultParameterInLowBins=-1;
  fCarefulFitting=false;
  fTotalCalls=0;
  scan2d=false;
  fUseLogL=true;
}

AliMCLogLFitter::~AliMCLogLFitter()
{
  // Destructor
  delete fFitDiagsMCExpected;
  delete fBestParameters;
  delete fParameter;
  delete fStartParameter;
  delete fUpperLimitParameter;
  delete fLowerLimitParameter;
}

void AliMCLogLFitter::SetFitRange(double lowerBound, double upperBound)
{
  // Set Fit range using start and end point of the fit range in corrdinate values
  flowerFitRange=fFitDiagData->FindBin(lowerBound);
  fupperFitRange=fFitDiagData->FindBin(upperBound);
}

void AliMCLogLFitter::SetParameter(int fParameterNr, double startValue, double lowerBound, double upperBound)
{
  // Set Starting value and range of one parameter
  if(startValue>0 && startValue<1.0) fStartParameter[fParameterNr]=startValue;
  if(lowerBound>0 && lowerBound<1.0) fLowerLimitParameter[fParameterNr]=lowerBound; else fLowerLimitParameter[fParameterNr]=0.0;
  if(upperBound>lowerBound && upperBound<1.0) fUpperLimitParameter[fParameterNr]=upperBound; else  fUpperLimitParameter[fParameterNr]=1.0;
}

void AliMCLogLFitter::ChangeFunctions(Int_t nMCFunctions, TH1 ** MCDistributions, TH1 * dataHistogram)
{
  // Switch the fir templates for others, possible with a different number.
  delete fFitDiagsMCExpected;
  delete fParameter;
  
  fNPar=nMCFunctions;
  fParameter=new Double_t[fNPar];
  fFitDiagsMC=MCDistributions;
  fFitDiagsMCExpected=new TH1*[fNPar];
  fFitDiagData=dataHistogram;
  for(int i=0;i<fNPar;i++)
  {
    fFitDiagsMCExpected[i]= (TH1*)fFitDiagsMC[i]->Clone();
  }
  // Here there should be a check to ensure all input distibutions have the same number of bins
}


void AliMCLogLFitter::UpdateExpectationIteratively(Double_t * par, Int_t bin) 
{
  // Iterative procedure to minimize the A_ji within one bins
  // Only used internally
  double fji;
  int i=bin;
  double di=fFitDiagData->GetBinContent(i);
  double * nr = new double[fNPar];
  double p,q;
  double temp=0;
  int largestpar;
  
  
  for(int j=0;j<fNPar;j++)
  {
    nr[j]=fFitDiagsMC[j]->GetBinContent(i);
    fFitDiagsMCExpected[j]->SetBinContent(i,TMath::Max(nr[j],0.1));
    //fFitDiagsMCExpected[j]->SetBinContent(i,2.0);
  }
  
  for(int j=0;j<fNPar;j++) temp+=nr[j];
  if(temp<0.1)
  {
    temp=0;largestpar=0;
   for(int j=0;j<fNPar;j++)  if(par[j]>temp){largestpar=j;temp=par[j];}
   if(fDefaultParameterInLowBins>=0 && fDefaultParameterInLowBins<fNPar)largestpar=fDefaultParameterInLowBins; // this should release problems with low bins
   for(int j=0;j<fNPar;j++) fFitDiagsMCExpected[j]->SetBinContent(i,0.);  // Apr 27
   fFitDiagsMCExpected[largestpar]->SetBinContent(i,di/(par[largestpar]+1.));
  }
  else
  for(int iter=0;iter<fiterations;iter++)
  {
    for(int j=0;j<fNPar;j++)
    {
      fji=0.0;
      for(int k=0;k<fNPar;k++)if(k!=j)fji+=par[k]*fFitDiagsMCExpected[k]->GetBinContent(i);
      p=-(di*par[j]-fji*(par[j]+1)+nr[j]*par[j])/(par[j]*(par[j]+1));
      q=-(fji*nr[j])/(par[j]*(par[j]+1));
      if(nr[j]<0.8 && fDefaultParameterInLowBins!=j)
	fFitDiagsMCExpected[j]->SetBinContent(i,0.);
      else
	fFitDiagsMCExpected[j]->SetBinContent(i,TMath::Max(-p/2.+TMath::Sqrt(p*p/4.0-q),0.0)); // 0.001
    }
  }
  delete nr;
}


void AliMCLogLFitter::ReturnNegLog(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) 
{
  // Calculate the -LogL. Form used to allow as input for TMinuit
  if(fUseLogL){
  if(gin || npar || iflag==0){}
  // par[0]/[1] are amplitudes
  Double_t LogLikelihood=0;
  Double_t fi;
  for(int i=flowerFitRange;i<fupperFitRange;i++)
  {
    //add=0;
    UpdateExpectationIteratively(par, i);
    fi=0;
    for(int p=0;p<npar;p++)
      fi+=par[p]*fFitDiagsMCExpected[p]->GetBinContent(i);
    if(fi<0.00001)fi=0.00001;
    
    LogLikelihood+=fFitDiagData->GetBinContent(i)*TMath::Log(fi)-fi;
    //add+=fFitDiagData->GetBinContent(i)*TMath::Log(fi)-fi;
    for(int p=0;p<npar;p++)
      LogLikelihood+=fFitDiagsMC[p]->GetBinContent(i)*TMath::Log(fFitDiagsMCExpected[p]->GetBinContent(i)+0.0001)-(fFitDiagsMCExpected[p]->GetBinContent(i));
  }
  f=-LogLikelihood;
  
  }
  else // Calculation for Chi Square
  {
    Double_t Chisq=0, up, low;
    int diff= fupperFitRange-flowerFitRange;
    
      for(int i=flowerFitRange;i<fupperFitRange;i++) 
    {
      up=fFitDiagData->GetBinContent(i);
      low=fFitDiagData->GetBinContent(i);
      for(int j=0;j<fNPar;j++)
      {
	up-=fFitDiagsMC[j]->GetBinContent(i)*par[j];
	low+=fFitDiagsMC[j]->GetBinContent(i)*par[j]*par[j];
      }
      up*=up;
      if(low<7.0) low=7.0;
      Chisq+=up/low/double(diff+1-fNPar);
    }
    f=Chisq;
  }
  fTotalCalls++;
  
}

Double_t AliMCLogLFitter::GetLogLikelyhood(void)
{
  // Returns the current LogL.
  double likelyhood;
  ReturnNegLog(fNPar, 0, likelyhood, fParameter, 0); 
  
  return likelyhood;
}

TH1D * AliMCLogLFitter::ScanVariable(Int_t fParameterNr, Double_t min, Double_t max)
{
  // Scans one variable, finding maximum likelyhood for all others
  // This is not a marginal!
  int maxiter=30;
  
  TH1D * scanHist = new TH1D("scanHist","scanHist", maxiter, min, max);
  std::cout << "Scanning Likelyhood for fParameter " << fParameterNr << endl;
  int tempint; double temp;
  
  double bestglobal=1e10;
  double * bestparaGlobal = new double[4];
  
  double start, upper, lower;
  
  //double increment=(max-min)/double(maxiter);
  
  TMinuit * minimizingObject = new TMinuit(fNPar);
  minimizingObject->SetName("scanFit2");
  minimizingObject->SetPrintLevel(-1);
  SomewhereElseInAliMCLogLFitter::fitter=this;  // Workaround
  minimizingObject->SetFCN(SomewhereElseInAliMCLogLFitter::OutsourcedReturnNegLog);
  fBestLikelyhood=1;
  
  

  for(int ScanBin=0;ScanBin<maxiter;ScanBin++)
  {
    if(ScanBin%20==0)std::cout << "scanning bin " << ScanBin << endl;
    for(int fitIter=0;fitIter<2;fitIter++)
    {
      for(int i=0;i<fNPar;i++)
      {
	if(i!=fParameterNr){
	start=  fFitDiagData->Integral()/fFitDiagsMC[i]->Integral()*fStartParameter[i];
	upper=  fFitDiagData->Integral()/fFitDiagsMC[i]->Integral()*fUpperLimitParameter[i];
	lower=  fFitDiagData->Integral()/fFitDiagsMC[i]->Integral()*fLowerLimitParameter[i];
	minimizingObject->mnparm(i,Form("par[%d]",i), start, start/10.0*(fitIter+1), lower, upper, tempint);
	}
      }
      //minimizingObject->mnparm(fParameterNr,Form("par[%d]",fParameterNr), (ScanBin+0.5)*increment+min, start/10.0, lower, upper, tempint);
      minimizingObject->mnparm(fParameterNr,Form("par[%d]",fParameterNr), scanHist->GetXaxis()->GetBinCenter(ScanBin+1), 0.0, 0, 1000, tempint);
	//minimizingObject->FixParameter(fParameterNr);
	//minimizingObject->mnfixp(fParameterNr,temp);
	//Fit happens here!
	for(int nseek=0;nseek<fitIter;nseek++)minimizingObject->mnseek();
	minimizingObject->Migrad();
	// end of Fit - make more involved later on
	for(int i=0;i<fNPar;i++)  minimizingObject->GetParameter(i, fParameter[i], temp);
	if(GetLogLikelyhood()<fBestLikelyhood || fitIter==0)
	{
	  fBestLikelyhood=GetLogLikelyhood();
	  for(int i=0;i<fNPar;i++)  fBestParameters[i]=fParameter[i];
	}
	//scanHist->SetBinContent(ScanBin+1,fBestLikelyhood);
    }
        scanHist->SetBinContent(ScanBin+1,fBestLikelyhood);
        if(bestglobal>fBestLikelyhood){bestglobal=fBestLikelyhood;for(int i=0;i<4;i++)bestparaGlobal[i]=fBestParameters[i];}
  }
  
  std::cout << "Best Likelyhood: " << bestglobal << endl;
  std::cout << "Best Parameters of scan: " << bestparaGlobal[0] << " " << bestparaGlobal[1] << " " << bestparaGlobal[2] << " " << bestparaGlobal[3] << " " << endl;  
  
  std::cout << "Another calculation. Setting fParameters." << endl;
  for(int i=0;i<4;i++)fParameter[i]=bestparaGlobal[i];
  std::cout << "Now: " << GetLogLikelyhood() << endl << "End of insertion." << endl;
  fParameter[0]=0.588;fParameter[1]=0.495;fParameter[0]=0.261;fParameter[0]=0.028;
  std::cout << "Test with values from fit: " << GetLogLikelyhood() << endl;
  
  delete minimizingObject;
  return scanHist;
}


TH2D * AliMCLogLFitter::ScanVariable2d(Int_t fParameterNr, Double_t min, Double_t max, Int_t fParameterNr2, Double_t min2, Double_t max2)
{
  // Scanning two parameters
  int maxiter=20;
  
  TH2D * scanHist = new TH2D("scanHist","scanHist", maxiter, min, max, maxiter, min2, max2);
  std::cout << "Scanning Likelyhood for fParameter " << fParameterNr << " and " << fParameterNr2 << endl;
  int tempint; double temp;
  double start, upper, lower;
  TMinuit * minimizingObject = new TMinuit(fNPar-2);
  minimizingObject->SetName("scanFit");
  minimizingObject->SetPrintLevel(-1);
  SomewhereElseInAliMCLogLFitter::fitter=this;  // Workaround
  minimizingObject->SetFCN(SomewhereElseInAliMCLogLFitter::OutsourcedReturnNegLog);
  fBestLikelyhood=1;
  
  double bestglobal=1e10;
  double * bestparaGlobal = new double[4];
  
  int N=0;
  //int fNParUsed=0;

  for(int ScanBin=1;ScanBin<=maxiter;ScanBin++)
    for(int ScanBin2=1;ScanBin2<=maxiter;ScanBin2++)
  {
    N++;
    if(N%2==0)cout << "\rscanning " << double(N)/double(maxiter*maxiter)*100.0 << "%        " << flush;
    for(int fitIter=0;fitIter<2;fitIter++)
    {
      
      for(int i=0;i<fNPar;i++)
      {
	if(i!=fParameterNr && i!=fParameterNr2){
	start=  fFitDiagData->Integral()/fFitDiagsMC[i]->Integral()*fStartParameter[i];
	upper=  fFitDiagData->Integral()/fFitDiagsMC[i]->Integral()*fUpperLimitParameter[i];
	lower=  fFitDiagData->Integral()/fFitDiagsMC[i]->Integral()*fLowerLimitParameter[i];
	minimizingObject->mnparm(i,Form("par[%d]",i), start*3.0, start/10.0*(fitIter+1), lower, upper, tempint);
	//cout << "Start Par " << i << ": " << start << endl;
	}
      }
      //minimizingObject->mnparm(fParameterNr,Form("par[%d]",fParameterNr), (ScanBin+0.5)*increment+min, start/10.0, lower, upper, tempint);
        minimizingObject->mnparm(fParameterNr,Form("par[%d]",fParameterNr), scanHist->GetXaxis()->GetBinCenter(ScanBin), 0.00001, scanHist->GetXaxis()->GetBinCenter(ScanBin)-0.00001, scanHist->GetXaxis()->GetBinCenter(ScanBin)+0.00001, tempint);
	//minimizingObject->FixParameter(fParameterNr);
        //minimizingObject->mnfixp(fParameterNr,temp);
        minimizingObject->mnparm(fParameterNr2,Form("par[%d]",fParameterNr2), scanHist->GetYaxis()->GetBinCenter(ScanBin2), 0.00001,  scanHist->GetYaxis()->GetBinCenter(ScanBin2)-0.00001,  scanHist->GetYaxis()->GetBinCenter(ScanBin2)+0.00001, tempint);
	//minimizingObject->FixParameter(fParameterNr2);
	//minimizingObject->mnfixp(fParameterNr2,temp);
	//cout << "checking : " <<  scanHist->GetXaxis()->GetBinCenter(ScanBin) << " - " << scanHist->GetYaxis()->GetBinCenter(ScanBin2) << endl;
	//Fit happens here!
	for(int nseek=0;nseek<fitIter;nseek++)minimizingObject->mnseek();
	minimizingObject->Migrad();
	// end of Fit - make more involved later on
	for(int i=0;i<fNPar;i++)  minimizingObject->GetParameter(i, fParameter[i], temp);
	if(GetLogLikelyhood()<fBestLikelyhood || fitIter==0)
	{
	  fBestLikelyhood=GetLogLikelyhood();
	  for(int i=0;i<fNPar;i++)  fBestParameters[i]=fParameter[i];
	}
	scanHist->SetBinContent(ScanBin, ScanBin2,fBestLikelyhood);
	//cout << "B end: " <<  fParameter[3] << endl;
    }
    if(bestglobal>fBestLikelyhood){bestglobal=fBestLikelyhood;for(int i=0;i<4;i++)bestparaGlobal[i]=fBestParameters[i];}
  }
  delete minimizingObject;
  cout << endl;
  cout << "Best Likelyhood: " << bestglobal << endl;
  cout << "Best Parameters of scan: " << bestparaGlobal[0] << " " << bestparaGlobal[1] << " " << bestparaGlobal[2] << " " << bestparaGlobal[3] << " " << endl;
  
  return scanHist;
}

void AliMCLogLFitter::ScanAll(void)
{
  // Start a global maximization
  double * bestparaGlobal = new double[4];
  
  int N=0;
  double min[4]={0.5,0.4,0.2,0.0};
  double max[4]={0.7,0.6,0.35,0.2};
  double steps=20;
  double stepSize[4];
  for(int i=0;i<4;i++) stepSize[i]=(max[i]-min[i])/steps;
  double currentLikelyhood,fBestLikelyhoodLocal=0;
  int ScanBin[4];
  std::cout << endl;
  for(ScanBin[0]=1;ScanBin[0]<=steps;ScanBin[0]++)
    for(ScanBin[1]=1;ScanBin[1]<=steps;ScanBin[1]++)
        for(ScanBin[2]=1;ScanBin[2]<=steps;ScanBin[2]++)
	    for(ScanBin[3]=1;ScanBin[3]<=steps;ScanBin[3]++)
  {
    N++;
    if(N%20==0) std::cout << "\r" << double(N)/double(steps*steps*steps*steps)*100.0 << "%" << flush;
    for(int i=0;i<4;i++)fParameter[i]=min[i]+stepSize[i]*ScanBin[i];
    currentLikelyhood=GetLogLikelyhood();
    
    if(currentLikelyhood<fBestLikelyhoodLocal)
    {
      for(int i=0;i<4;i++)bestparaGlobal[i]=fParameter[i];
      fBestLikelyhoodLocal=currentLikelyhood;
    }
  }
  std::cout << endl;
  std::cout << "Scan of all. Best Likelyhood: " << fBestLikelyhoodLocal << endl
  << "best paramters: " << bestparaGlobal[0] << " " << bestparaGlobal[1] << " " << bestparaGlobal[2] << " " << bestparaGlobal[3] << "\n end of total scan." << endl;

}



void AliMCLogLFitter::Fit(void)
{
  scan2d=false;
  int tempint; double temp;
  
  double start, upper, lower;
  
  TMinuit * minimizingObject = new TMinuit(fNPar);
  minimizingObject->SetName("globalFit");
  minimizingObject->SetPrintLevel(-1);
  SomewhereElseInAliMCLogLFitter::fitter=this;  // Workaround
  minimizingObject->SetFCN(SomewhereElseInAliMCLogLFitter::OutsourcedReturnNegLog);
  fBestLikelyhood=1;
  //std::cout << "Par Min Max" << endl;
  for(int fitIter=0;fitIter<5;fitIter++)
  {
    for(int i=0;i<fNPar;i++)
    {
      start=  fFitDiagData->Integral(1,fFitDiagData->GetNbinsX())/fFitDiagsMC[i]->Integral(1,fFitDiagsMC[i]->GetNbinsX())*fStartParameter[i];
      upper=  fFitDiagData->Integral(1,fFitDiagData->GetNbinsX())/fFitDiagsMC[i]->Integral(1,fFitDiagsMC[i]->GetNbinsX())*fUpperLimitParameter[i];
      lower=  fFitDiagData->Integral(1,fFitDiagData->GetNbinsX())/fFitDiagsMC[i]->Integral(1,fFitDiagsMC[i]->GetNbinsX())*fLowerLimitParameter[i];
      minimizingObject->mnparm(i,Form("par[%d]",i), start, start/10.0*(fitIter+1), lower, upper, tempint);
      //std::cout << i << " " << lower << " " << upper << endl;
    }
    //std::cout << "end of Pars." << endl;
    //std::cout << "Fit range: " << flowerFitRange << " to " << fupperFitRange << endl;
    //Fit happens here!
    for(int nseek=0;nseek<fitIter;nseek++)minimizingObject->mnseek();
    minimizingObject->Migrad();
    //minimizingObject->mnimpr();
    // end of Fit - make more involved later on
    for(int i=0;i<fNPar;i++)  minimizingObject->GetParameter(i, fParameter[i], temp);
    //std::cout << "Iter " << fitIter << " likelyhood: " << GetLogLikelyhood();
    //for(int i=0;i<fNPar;i++) std::cout << " par[" << i << "] = " << fParameter[i];
    //std::cout << endl;
    if(GetLogLikelyhood()<fBestLikelyhood || fitIter==0)
    {
      fBestLikelyhood=GetLogLikelyhood();
      for(int i=0;i<fNPar;i++)  fBestParameters[i]=fParameter[i];
    }
  }
  
  //std::cout << "Likelyhood: " << fBestLikelyhood << endl;
  //std::cout << "fParameters: " ;
  //for(int i=0;i<fNPar;i++) { fParameter[i]=fBestParameters[i];std::cout << " par[" << i << "] = " << fParameter[i];}
  //std::cout << endl;
  
  
  //  fParameter[0]=0.59*2;fParameter[1]=0.49*2;fParameter[0]=0.338*2;fParameter[0]=0.0001*2;
  //std::cout << "Test with values from scan: " << GetLogLikelyhood() << endl;
  
  
  
  /*int atLimits=-1;
  if(fCarefulFitting)
  {
    std::cout << "Fitting with particular care." << endl;
    for(int i=0;i<fNPar;i++)if(fParameter[i]/(fFitDiagData->Integral()/fFitDiagsMC[i]->Integral()*fStartParameter[i])<10.0){atLimits=i;}
    if(atLimits>-1)
    {
      scanHist=ScanVariable(atLimits, 0.0, 0.5);
      scanHist->Smooth(5);
      double bestPar = scanHist->GetBinCenter(scanHist->GetMinimumBin());
      for(int i=0;i<fNPar;i++)
      {
	start=  fFitDiagData->Integral()/fFitDiagsMC[i]->Integral()*fStartParameter[i];
	upper=  fFitDiagData->Integral()/fFitDiagsMC[i]->Integral()*fUpperLimitParameter[i];
	lower=  fFitDiagData->Integral()/fFitDiagsMC[i]->Integral()*fLowerLimitParameter[i];
	minimizingObject->mnparm(i,Form("par[%d]",i), start, start/10.0, lower, upper, tempint);
      }
      minimizingObject->mnparm(atLimits,Form("par[%d]",atLimits), bestPar, start/10.0, lower, upper, tempint);
      minimizingObject->FixParameter(atLimits);
      std::cout << "Improvement found for fParameter " << atLimits << " good value: " << bestPar << endl;
      minimizingObject->Migrad();
      minimizingObject->Release(atLimits);
      minimizingObject->mnparm(atLimits,Form("par[%d]",atLimits), bestPar, start/10.0, bestPar/2.0, bestPar*2.0, tempint);
      minimizingObject->Migrad();
      for(int i=0;i<fNPar;i++)  fBestParameters[i]=fParameter[i];
      for(int i=0;i<fNPar;i++)  fParameter[i]=fBestParameters[i];      
      std::cout << "Final value: " << fParameter[atLimits] << endl;
    }
  }*/
  //std::cout << "fTotalCalls: " << fTotalCalls << endl;
  delete minimizingObject;
}
