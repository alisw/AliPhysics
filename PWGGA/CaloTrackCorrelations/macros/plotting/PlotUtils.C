///
/// \file PlotUtils.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Plotting utilities
///
/// Macro that containes useful small methods
/// that can be used in plotting macros
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraphErrors.h>

#endif

///
/// Crystal ball function. Fitting function.
///
//-----------------------------------------------------------------------------
static Double_t FunctionCrystalBall(Double_t *x, Double_t *par)
{
  Double_t N = par[0];
  Double_t width = par[1];
  
  Double_t mean = par[2];
  Double_t sigma = par[3];
  Double_t alpha = par[4];
  Double_t n = par[5];
  
  Double_t A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  Double_t B = n/fabs(alpha) - fabs(alpha);
  
  if ((x[0]-mean)/sigma>-alpha)
    return N*width*TMath::Gaus(x[0],mean,sigma,1);
  else
    return N/(sqrt(TMath::TwoPi())*sigma)*width*A*pow(B-(x[0]-mean)/sigma,-n);
}

///
/// Gaussian plus polinomial order 3 function. Fitting function.
///
//-----------------------------------------------------------------------------
static Double_t FunctionGaussPol3(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0];
  return gaus+back;
}

///
/// Gaussian plus polinomial order 2 function. Fitting function.
///
//-----------------------------------------------------------------------------
static Double_t FunctionGaussPol2(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0];
  return gaus+back;
}

///
/// Gaussian plus polinomial order 1 function. Fitting function.
///
//-----------------------------------------------------------------------------
static Double_t FunctionGaussPol1(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0];
  return gaus+back;
}

///
/// Gaussian plus constant function. Fitting function.
///
//-----------------------------------------------------------------------------
static Double_t FunctionGaussPol0(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = 0;//par[3];
  return gaus+back;
}

///
/// Truncated function in eta region. Fitting function.
///
//-----------------------------------------------------------------------------
static Double_t FunctionTruncatedPolEta(Double_t *x, Double_t *par)
{
  if ( x[0] > 0.425 && x[0] < 0.65 ) 
  {
    TF1::RejectPoint();
    return 0;
  }
  
  return par[0] + par[1]*x[0];
  // + x[0]*x[0]*par[2]  + x[0]*x[0]*x[0]*par[3];//+ x[0]*x[0]*x[0]*x[0]*par[4];
}

///
/// Truncated function in pi0 region. Fitting function.
///
//-----------------------------------------------------------------------------
static Double_t FunctionTruncatedPolPi0(Double_t *x, Double_t *par)
{
  //if ( x[0] > 0.07 && x[0] < 0.2 )
  if ( x[0] < 0.25 )
  {
    TF1::RejectPoint();
    return 0;
  }
  
  return par[0] + par[1]*x[0];
  // + x[0]*x[0]*par[2] + x[0]*x[0]*x[0]*par[3] ;//+ x[0]*x[0]*x[0]*x[0]*par[4];
}

//-----------------------------------------------------------------------------
/// Calculation of error of a fraction
///
/// \return fraction error
/// \param  num numerator
//-----------------------------------------------------------------------------
static Double_t GetFractionError
(Double_t num   , Double_t den, 
 Double_t numErr, Double_t denErr)
{
  if ( num == 0 || den == 0 ) return 0.;
  
  //printf("\t num %e den %e numErr %e denErr %e\n",num,den,numErr,denErr);
  
  return num/den * TMath::Sqrt( ( numErr * numErr ) / ( num * num ) + 
                                ( denErr * denErr ) / ( den * den )   );
}

//-----------------------------------------------------------------------------
/// Divide 2 TGraphErrors
///
/// \return TGraphError result of division
/// \param  gNum TGraphError numerator
/// \param  gDen TGraphError denominator
//-----------------------------------------------------------------------------
static TGraphErrors * DivideGraphs
(TGraphErrors* gNum, TGraphErrors *gDen)
{
  if ( !gDen || !gNum ) 
  {
    printf("Graph num %p or graph den %p not available\n",gNum,gDen);
    return 0x0;
  }
  const Int_t nBins = gNum->GetN();
  if ( nBins != gDen->GetN() )
  {
    printf("Cannot divide %s with %d bins and %s with %d bins!\n",
           gNum->GetName(),nBins,gDen->GetName(),gDen->GetN());
    return 0x0;
  }
  
  Double_t ratio   [nBins];
  Double_t ratioErr[nBins]; 
  Double_t x       [nBins];
  Double_t xErr    [nBins];
  for (Int_t ibin = 0; ibin < nBins; ibin++) 
  {
    Double_t num    =  gNum->GetY ()[ibin];
    Double_t den    =  gDen->GetY ()[ibin];
    Double_t numErr =  gNum->GetEY()[ibin];
    Double_t denErr =  gDen->GetEY()[ibin];
    
    x   [ibin]      =  gNum->GetX ()[ibin];
    xErr[ibin]      =  gNum->GetEX()[ibin];
    
    if ( num == 0 || den == 0 ) 
    {
      ratio   [ibin] = 0; 
      ratioErr[ibin] = 0; 
      continue;
    }
    
    ratio   [ibin] = num / den ;
    ratioErr[ibin] =  GetFractionError(num,den,numErr,denErr);  
    
    //    printf("bin %d, x %f (%f) num %f (%f), den %f (%f), ratio %f (%f) \n",
    //           ibin,x[ibin],xErr[ibin],num,numErr,den,denErr,ratio[ibin],ratioErr[ibin]);
  } // do the ratio to sum
  
  return new TGraphErrors(nBins,x,ratio,xErr,ratioErr);
}


//-----------------------------------------------------------------------------
/// Get a TH1D integral and its error within a bin range
///
/// \param h TH1D 
/// \param binMin First bin of integration
/// \param binMax Last bin of integration
/// \param integral sum between first and last bin
/// \param integralErr error of integral
//-----------------------------------------------------------------------------
static void GetRangeIntegralAndError
(TH1D* h, Int_t binMin, Int_t binMax, 
 Double_t & integral, Double_t & integralErr )
{
  integral    = 0;
  integralErr = 0;
  for(Int_t ibin = binMin; ibin <= binMax; ibin++)
  {
    //if ( h->GetBinContent(ibin) == 0 ) continue ; 
    
    integral    += h->GetBinContent(ibin);
    integralErr += h->GetBinError  (ibin) * h->GetBinError(ibin);
  }
//  printf("\t bin range [%d,%d], n %d, sum %2.2e err %2.2e\n",
//         binMin,binMax,n,integral,integralErr);
  if ( integralErr > 0 ) integralErr = TMath::Sqrt(integralErr);
}

//-----------------------------------------------------------------------------
/// Get a TH1D integral and its error within a x axis range
///
/// \param h TH1D 
/// \param minX x axis minimum value
/// \param maxX x axis maximum value
/// \param integral sum between first and last bin
/// \param integralErr error of integral
//-----------------------------------------------------------------------------
static void GetRangeIntegralAndError
(TH1D* h, Float_t minX, Float_t maxX, 
 Double_t & integral, Double_t & integralErr )
{
  Int_t binMin = h->FindBin(minX);
  Int_t binMax = h->FindBin(maxX)-1;
  return GetRangeIntegralAndError(h, binMin, binMax, 
                                  integral, integralErr );
}

//-----------------------------------------------------------------------------
/// When plotting multiple graphs, and need to set different ranges
/// Those points with too large error or out of the desired range are set
/// to a default value out of the plot scale
///
/// \param graph TGraphError 
/// \param min x axis minimum value
/// \param max x axis maximum value
/// \param errFraction If error of point is larger than a fraction of the point, it is removed
/// \param value default value for removed point
/// \param valueErr default error value for removed point
//-----------------------------------------------------------------------------
void RemovePointsOutOfRangeOrLargeErrorFromGraph
(TGraphErrors* graph, 
 Double_t min, Double_t max,
 Float_t errFraction = 0.5,
 Float_t value = -1, Float_t valueErr = 0 )
{
  for(Int_t ibin = 0; ibin < graph->GetN(); ibin++)
  {
    Float_t x  = graph->GetX()[ibin];
    Float_t y  = graph->GetY()[ibin];
    Float_t yE = graph->GetEY()[ibin];
    if(x < min || x > max )
    {
      graph->SetPoint     (ibin,x, value   );
      graph->SetPointError(ibin,x, valueErr);
    }
    else 
      if(yE > errFraction*y) // Remove points with more than 50% error
      {
        //printf("ibin %d, y %f, yE %f\n",ibin,y,yE);
        graph->SetPoint     (ibin,x,value   );
        graph->SetPointError(ibin,x,valueErr);
      }
  }
}

//-----------------------------------------------------------------------------
/// Scale histogram bins by its size
//-----------------------------------------------------------------------------
static void ScaleBinBySize(TH1D* h)
{
  for(Int_t ibin = 1; ibin < h->GetNbinsX();ibin++)
  {
    Double_t width   = h->GetBinWidth(ibin);
    Double_t content = h->GetBinContent(ibin);
    Double_t error   = h->GetBinError(ibin);
    
    //printf("bin %d, width %f, content %e\n",ibin,width,content);
    h->SetBinContent(ibin,content/width);
    h->SetBinError  (ibin,  error/width);
  }
}

//-----------------------------------------------------------------------------
/// Scale 2D histogram X bins by its integral
//-----------------------------------------------------------------------------
void ScaleXaxis2D(TH2F* h2D)
{
  Int_t nbinsy = h2D->GetNbinsY();
  
  for(Int_t j = 1; j <= h2D->GetNbinsX(); j++)
  {
    TH1D* temp1 = h2D->ProjectionY(Form("Bin%d",j),j,j);
    
    Float_t scale1 = temp1 -> Integral(-1,-1);
    
    for(Int_t i = 1; i <= nbinsy; i++)
    {
      //printf("i %d, j %d;  content %f / scale %f = %f\n",
      //i,j,h2D->GetBinContent(j,i),scale2,h2D->GetBinContent(j,i)/scale2);
      
      if ( scale1 > 0 )
      {
        h2D->SetBinContent(j,i, h2D->GetBinContent(j,i)/scale1);
        h2D->SetBinError  (j,i, h2D->GetBinError  (j,i)/scale1);
      }
      else
      {
        h2D->SetBinContent(j,i, 0);
        h2D->SetBinError  (j,i, 0);       
      }
      
    } // y bin loop 
  } // x bin loop
  
  h2D->SetZTitle("x bin norm. to integral");
}

//-----------------------------------------------------------------------------
/// Scale 2D histogram X bins by its integral and scale y bin by size
//-----------------------------------------------------------------------------
void ScaleXaxisIntegralYAxisSize2D(TH2F* h2D)
{
  Int_t nbinsy = h2D->GetNbinsY();
  
  for(Int_t j = 1; j <= h2D->GetNbinsX(); j++)
  {
    TH1D* temp1 = h2D->ProjectionY(Form("Bin%d",j),j,j);
    
    Float_t scale1 = temp1 -> Integral(-1,-1);
    
    for(Int_t i = 1; i <= nbinsy; i++)
    {
      //printf("i %d, j %d;  content %f / scale %f = %f\n",
      //i,j,h2D->GetBinContent(j,i),scale2,h2D->GetBinContent(j,i)/scale2);
      
      if ( scale1 > 0 )
      {
        h2D->SetBinContent(j,i, h2D->GetBinContent(j,i)/scale1);
        h2D->SetBinError  (j,i, h2D->GetBinError  (j,i)/scale1);
      }
      else
      {
        h2D->SetBinContent(j,i, 0);
        h2D->SetBinError  (j,i, 0);       
      }
      
    } // y bin loop 
  } // x bin loop
  
  
  // Normalize by y bin size
   
   for(Int_t j = 1; j <= h2D->GetNbinsX(); j++)
   {
     
     for(Int_t i = 1; i <= nbinsy; i++)
     {
       //printf("NLM2: i %d, j %d;  content %f / scale %f = %f\n",i,j,hNLM2->GetBinContent(j,i),scale2,hNLM2->GetBinContent(j,i)/scale2);
         h2D->SetBinContent(j,i, h2D->GetBinContent(j,i)/h2D->GetYaxis()->GetBinWidth(i));
         h2D->SetBinError  (j,i, h2D->GetBinError  (j,i)/h2D->GetYaxis()->GetBinWidth(i));
     } // y bin loop 
   } // x bin loop
  
  h2D->SetZTitle("x bin norm. to integral, y bin by size");

}


//-----------------------------------------------------------------------------
/// When more than 1 frame pad in a canvas, define the number of 
/// columns or rows, depending on the total number of frames
/// \param npad number of pads
/// \param ncol number of assigned columns
/// \param nrow number of assigned rows
/// \param square take the number of columns and rows as the square root of npad
//-----------------------------------------------------------------------------
void GetCanvasColRowNumber(Int_t npad, Int_t & ncol, Int_t & nrow, 
                           Bool_t square = kFALSE)
{
  if(npad <= 0)
  {
    ncol = 0;
    nrow = 0;
    return;
  }
  
  ncol = TMath::Sqrt(npad);
  nrow = ncol; 
  
  if ( square ) 
    return;
  
  if ( npad == 1 )
  {
    ncol = 1;
    nrow = 1;
  }
  else if ( npad == 2 )
  {
    ncol = 2;
    nrow = 1;
  }
  else if ( npad == 3 )
  {
    ncol = 3;
    nrow = 1;
  }
  else if ( npad == 4 )
  {
    ncol = 2;
    nrow = 2;
  }
  else if ( npad < 7 )
  {
    ncol = 3;
    nrow = 2;
  }
  else if ( npad < 10 )
  {
    ncol = 3;
    nrow = 3;
  }  
  else if ( npad < 13 )
  {
    ncol = 4;
    nrow = 3;
  }    
  else if ( npad < 17 )
  {
    ncol = 4;
    nrow = 4;
  }  
  else if ( npad < 21 )
  {
    ncol = 5;
    nrow = 4;
  }  
  else if ( npad < 26 )
  {
    ncol = 5;
    nrow = 5;
  }    
  else if ( npad < 31 )
  {
    ncol = 6;
    nrow = 5;
  }  
  else  if ( npad < 37 )
  {
    ncol = 6;
    nrow = 6; 
  }

}

