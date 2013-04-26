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
//
// Toolkit containing various usefull things
// Usable everywhere in the hfe software package
// For more information see the cxx file
//
// Authors
//   All authors of the HFE group
//
#include <TArrayD.h>
#include <TMath.h>
#include <TParticle.h>
#include <TF1.h>
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TKey.h"
#include "TROOT.h"

#include "AliAODMCParticle.h"
#include "AliAODpidUtil.h"
#include "AliESDpid.h"
#include "AliLog.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"

#include "AliHFEtools.h"

ClassImp(AliHFEtools)

AliPIDResponse *AliHFEtools::fgDefaultPID = NULL;
Int_t AliHFEtools::fgLogLevel = 0;

//__________________________________________
AliHFEtools::AliHFEtools():
  TObject()
{
}

//__________________________________________
Double_t *AliHFEtools::MakeLinearBinning(Int_t nBins, Double_t ymin, Double_t ymax){
  //
  // Helper function for linearly binned array
  //
  Double_t *bins = new Double_t[nBins + 1];
  Double_t stepsize = (ymax - ymin) / static_cast<Double_t>(nBins);
  bins[0] = ymin;
  for(Int_t ibin = 1; ibin <= nBins; ibin++)
    bins[ibin] = bins[ibin-1] + stepsize;
  return bins;
}

//__________________________________________
void AliHFEtools::FillLinearBinning(TArrayD &bins, Int_t nBins, Double_t ymin, Double_t ymax){
  //
  // Helper function for linearly binned array
  //
  Double_t stepsize = (ymax - ymin) / static_cast<Double_t>(nBins);
  bins[0] = ymin;
  for(Int_t ibin = 1; ibin <= nBins; ibin++)
    bins[ibin] = bins[ibin-1] + stepsize;
}

//__________________________________________
Double_t *AliHFEtools::MakeLogarithmicBinning(Int_t nBins, Double_t ymin, Double_t ymax){
  //
  // Helper function for logartimically binned array
  //
  Double_t *bins = new Double_t[nBins+1];
  bins[0] = ymin;
  Double_t factor = TMath::Power(ymax/ymin, 1./nBins);
  for(Int_t ibin = 1; ibin <= nBins; ibin++){
    bins[ibin] = factor * bins[ibin-1];
  }
  return bins;
}

//__________________________________________
void AliHFEtools::FillLogarithmicBinning(TArrayD &bins, Int_t nBins, Double_t ymin, Double_t ymax){
  //
  // Helper function for logartimically binned array
  //
  bins[0] = ymin;
  Double_t factor = TMath::Power(ymax/ymin, 1./nBins);
  for(Int_t ibin = 1; ibin <= nBins; ibin++)
    bins[ibin] = factor * bins[ibin-1];
}

//_________________________________________
Bool_t AliHFEtools::BinLogAxis(TObject *o, Int_t dim){

  // 
  // converts the axis (defined by the dimension) of THx or THnSparse
  // object to Log scale. Number of bins and bin min and bin max are preserved
  //


  if(!o){
    AliError("Input histogram is null pointer");
    return kFALSE;    
  }
  
  TAxis *axis = 0x0;
  if(o->InheritsFrom("TH1")){
    TH1 *h1 = dynamic_cast<TH1F*>(o); 
    if(h1) axis = h1->GetXaxis();
  }
  else if(o->InheritsFrom("TH2")){
    TH2 *h2 = dynamic_cast<TH2F*>(o);
    if(h2 && 0 == dim){
      axis = h2->GetXaxis();
    }
    else if(h2 && 1 == dim){
      axis = h2->GetYaxis();
    }
     else{
       AliError("Only dim = 0 or 1 possible for TH2F");
     }
  }
  else if(o->InheritsFrom("THnSparse")){
    THnSparse *hs = dynamic_cast<THnSparse*>(o);
    if(hs) axis = hs->GetAxis(dim);
  }
  else{
    AliError("Type of input object not recognized, please check your code or update this finction");
    return kFALSE;
  }
  if(!axis){
    AliError(Form("Axis '%d' could not be identified in the object \n", dim));
    return kFALSE;
  }
  
  Int_t bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  if(from <= 0){
    AliError(Form("Log binning not possible for this axis [min = %f]\n", from));
  }
  Double_t to = axis->GetXmax();
  TArrayD newBins(bins+1);
  newBins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i){
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins.GetArray());

  return kTRUE;
}

//__________________________________________
Float_t AliHFEtools::GetRapidity(const TParticle *part){
  //
  // return rapidity
  //
  Float_t rapidity;
  if(!((part->Energy() - part->Pz())*(part->Energy() + part->Pz())>0)) rapidity=-999;
  else rapidity = 0.5*(TMath::Log((part->Energy()+part->Pz()) / (part->Energy()-part->Pz())));
  return rapidity;
}

//__________________________________________
Float_t AliHFEtools::GetRapidity(const AliAODMCParticle *part){
  // return rapidity

  Float_t rapidity;        
  if(!((part->E() - part->Pz())*(part->E() + part->Pz())>0)) rapidity=-999; 
  else rapidity = 0.5*(TMath::Log((part->E()+part->Pz()) / (part->E()-part->Pz()))); 
  return rapidity;
}

//__________________________________________
AliPIDResponse* AliHFEtools::GetDefaultPID(Bool_t isMC, Bool_t isESD){
  //
  // Get the default PID as singleton instance
  //
  if(!fgDefaultPID){

    if(isESD) fgDefaultPID = new AliESDpid;
    else fgDefaultPID = new AliAODpidUtil;
    Double_t tres = isMC ? 80. : 130.;
    fgDefaultPID->GetTOFResponse().SetTimeResolution(tres);

    // TPC Bethe Bloch parameters
    Double_t alephParameters[5];
    if(isMC){
      // simulation
      alephParameters[0] = 2.15898e+00/50.;
      alephParameters[1] = 1.75295e+01;
      alephParameters[2] = 3.40030e-09;
      alephParameters[3] = 1.96178e+00;
      alephParameters[4] = 3.91720e+00;
    } else {
      alephParameters[0] = 0.0283086/0.97;
      //alephParameters[0] = 0.0283086;
      alephParameters[1] = 2.63394e+01;
      alephParameters[2] = 5.04114e-11;
      alephParameters[3] = 2.12543e+00;
      alephParameters[4] = 4.88663e+00;
    }
    fgDefaultPID->GetTPCResponse().SetBetheBlochParameters(alephParameters[0],alephParameters[1],alephParameters[2], alephParameters[3],alephParameters[4]);
    fgDefaultPID->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);

  }
  if(fgLogLevel){
    printf("Error - You are using the default PID: You should use the PID coming from the tender\n");
    printf("Error - Arrrrrrrrr...\n");
    printf("Error - Please rethink your program logic. Using default PID is really dangerous\n");
    printf("Error - TOF PID is adapted to Monte Carlo\n");
  }
  return fgDefaultPID;
}


//__________________________________________
void AliHFEtools::DestroyDefaultPID(){
  //
  // Destroy default PID object if existing
  //
  if(fgDefaultPID) delete fgDefaultPID;
  fgDefaultPID = NULL;
}

//__________________________________________
Int_t AliHFEtools::GetPdg(const AliVParticle *track){
  // 
  // Get MC PDG code (only MC particles supported)
  //
  Int_t pdg = 0;
  if(!TString(track->IsA()->GetName()).CompareTo("AliMCParticle")){
    const AliMCParticle *mctrack = dynamic_cast<const AliMCParticle *>(track);
    pdg = mctrack ? mctrack->Particle()->GetPdgCode() : 0;
  } else if(!TString(track->IsA()->GetName()).CompareTo("AliAODMCParticle")){
    const AliAODMCParticle *aodmctrack = dynamic_cast<const AliAODMCParticle *>(track);
    pdg = aodmctrack ? aodmctrack->GetPdgCode() : 0;
  }
  return pdg;
}

//__________________________________________
Int_t AliHFEtools::PDG2AliPID(Int_t pdg){
  // 
  // Helper function to convert MC PID codes into AliPID codes
  //
  Int_t pid = -1;
  switch(TMath::Abs(pdg)){
    case 11: pid = AliPID::kElectron; break;
    case 13: pid = AliPID::kMuon; break;
    case 211: pid = AliPID::kPion; break;
    case 321: pid = AliPID::kKaon; break;
    case 2212: pid = AliPID::kProton; break;
    default: pid = -1; break;
  };
  return pid;
}

//__________________________________________
TH1D* AliHFEtools::GraphErrorsToHist(TGraphErrors* g, Double_t firstBinWidth, Bool_t exchange, Int_t markerstyle,Int_t markercolor,Float_t markersize)
{
    // Creates a TH1D from TGraph g. The binwidth of the first bin has to
    // specified. The others bins are calculated automatically. Supports also Graphs
    // with non constant x steps. The axis of the Graph can be exchanged if
    // exchange=kTRUE (modifies the Graph).

    TH1D* result = GraphToHist(g,firstBinWidth,exchange, markerstyle,markercolor,markersize);
    if( result == 0) return result;

    //------------------------------------------
    // setup the final hist
    Int_t nBinX = g->GetN();
    Double_t* err = g->GetEY();           // error y is still ok even if exchanged
    for(Int_t i = 0; i < nBinX; i ++){
	result->SetBinError(i + 1, err[i]);
    }
    if(exchange){
        AliHFEtools::ExchangeXYGraph(g);        // undo  what has been done in GraphToHist
	AliHFEtools::ExchangeXYGraphErrors(g);  // now exchange including errors
    }

    return result;
}

//__________________________________________
Bool_t AliHFEtools::ExchangeXYGraph(TGraph* g)
{
    // exchanges x-values and y-values.
    if(g==0) return kFALSE;
    Int_t nbin=g->GetN();
    Double_t x,y;
    for(Int_t i = 0; i < nbin; i ++)
    {
        g->GetPoint(i,x,y);
        g->SetPoint(i,y,x);
    }

    return kTRUE;
}

//__________________________________________
Bool_t AliHFEtools::ExchangeXYGraphErrors(TGraphErrors* g)
{
    // exchanges x-values and y-values and
    // corresponding errors.
    if(g==0) return kFALSE;
    Int_t nbin=g->GetN();
    Double_t x,y;
    Double_t ex,ey;
    for(Int_t i = 0; i < nbin; i ++)
    {
        g->GetPoint(i,x,y);
        ex = g->GetErrorX(i);
        ey = g->GetErrorY(i);
        g->SetPoint(i,y,x);
        g->SetPointError(i,ey,ex);
    }

    return kTRUE;

}

//__________________________________________
TH1D* AliHFEtools::GraphToHist(TGraph* g, Double_t firstBinWidth, Bool_t exchange, Int_t markerstyle,Int_t markercolor,Float_t markersize)
{
    // Creates a TH1D from TGraph g. The binwidth of the first bin has to
    // specified. The others bins are calculated automatically. Supports also Graphs
    // with non constant x steps. The axis of the Graph can be exchanged if
    // exchange=kTRUE (modifies the Graph).


    TH1D* result = 0;
    if(g == 0)              return result;
    if(firstBinWidth == -1) return result;
    TString myname="";
    myname = g->GetName();
    myname += "_graph";
    if(exchange) AliHFEtools::ExchangeXYGraph(g);

    Int_t nBinX = g->GetN();
    Double_t* x = g->GetX();
    Double_t* y = g->GetY();

    if(nBinX < 1) return result;

    //------------------------------------------
    // create the Matrix for the equation system
    // and init the values

    Int_t nDim = nBinX - 1;
    TMatrixD a(nDim,nDim);
    TMatrixD b(nDim,1);

    Double_t* aA = a.GetMatrixArray();
    Double_t* aB = b.GetMatrixArray();
    memset(aA,0,nDim * nDim * sizeof(Double_t));
    memset(aB,0,nDim * sizeof(Double_t));
    //------------------------------------------

    //------------------------------------------
    // setup equation system
    // width for 1st bin is given therefore
    // we shift bin parameter (column) by one to the left
    // to reduce the matrix size

    Double_t* xAxis = new Double_t [nBinX + 1];
    Double_t* binW  = new Double_t [nBinX ];
    binW[0] = firstBinWidth;

    aB[0] = x[1] - x[0]  - 0.5 * binW[0];
    aA[0] = 0.5;

    for(Int_t col = 1; col < nDim ; col ++)
    {
	Int_t row = col;
	aB[col] = x[col + 1] - x[ col ];
	aA[row * nDim + col - 1 ] = 0.5;
	aA[row * nDim + col     ] = 0.5;
    }
    //------------------------------------------

    //------------------------------------------
    // solve the equations
    a.Invert();
    TMatrixD c = a * b;
    //------------------------------------------

    //------------------------------------------
    // calculate the bin boundaries
    xAxis[0] = x[0] - 0.5 * binW[0];
    memcpy(&binW[1],c.GetMatrixArray(),nDim * sizeof(Double_t));
    for(Int_t col = 0; col < nBinX ; col ++) {
	xAxis[col + 1] = x[col] + 0.5 * binW[col];
    }
    //------------------------------------------

    //------------------------------------------
    // setup the final hist
    result = new TH1D(myname.Data(),myname.Data(),nBinX, xAxis);
    for(Int_t i = 0; i < nBinX; i ++){
	result->SetBinContent(i + 1, y[i]);
    }
    result->SetMarkerColor(markercolor);
    result->SetMarkerStyle(markerstyle);
    result->SetMarkerSize(markersize);
    //------------------------------------------

    delete [] xAxis;
    delete [] binW;


    return result;
}

//__________________________________________
void AliHFEtools::BinParameterisation(const TF1 &fun, const TArrayD &xbins, TArrayD &bincontent){
    //
    // Calculate binned version of a function defined as the integral of x*f(x) in
    // the integration range xmin,xmax, where xmin and xmax are the bin limits, divided
    // by the binwidth. The function is important in case of steeply falling functions
    //
    // Parameters
    //   fun:           the function to be binned
    //   xbins:         the bin limits
    //   bincontent:    the binned parameterisation
    //
    TString expression(Form("x*%s", fun.GetName()));
    Double_t xmin(0), xmax(0);
    fun.GetRange(xmin,xmax);
    // check range
    xmin = TMath::Min(xmin, xbins[0]);
    xmax = TMath::Max(xmax, xbins[xbins.GetSize()-1]);
    TF1 helper("helper",expression.Data(),xmin,xmax);   // make function x*f(x)
    if(bincontent.GetSize() != xbins.GetSize()-1)
        bincontent.Set(xbins.GetSize()-1); // Adapt array to number of bins
    //Caclulate Binned
    for(Int_t ib = 0; ib < xbins.GetSize()-1; ib++){
        xmin = xbins[ib];
        xmax = xbins[ib+1];
        bincontent[ib] = (helper.Integral(xmin, xmax))/(xmax - xmin);
    }
}




//_________________________________________________________________________
//Function  AliHFEtools::GetHFEResultList() - opens file from argument and returns TList Object containing String "Results"
//_________________________________________________________________________
TList *AliHFEtools::GetHFEResultList(const TString str){

    TFile *f = TFile::Open(str.Data());
    if(!f || f->IsZombie()){
        printf("Could not read file %s\n",str.Data()); 
        return NULL ;
    }
    gROOT->cd();
    TKey *k;
    TIter next(f->GetListOfKeys());
    while ((k = dynamic_cast<TKey *>(next()))){
        TString s(k->GetName());
        if(s.Contains("Results")) break;
    }
    if(!k){
        printf("Output container not found\n");
        f->Close(); delete f;
        return NULL;
    } 
    TList *returnlist = dynamic_cast<TList *>(k->ReadObj());
    f->Close(); delete f;
    return returnlist;
}

//__________________________________________
void AliHFEtools::NormaliseBinWidth(TH1 *histo){
  //
  // Helper function to correct histograms for the bin width
  //
  Double_t binwidth(0.);
  for(Int_t ipt = 1; ipt <= histo->GetNbinsX(); ipt++){
    binwidth = histo->GetBinWidth(ipt);
    histo->SetBinContent(ipt, histo->GetBinContent(ipt)/binwidth);
    histo->SetBinError(ipt, histo->GetBinError(ipt)/binwidth);
  }
}

//__________________________________________
void AliHFEtools::NormaliseBinWdith(TGraphErrors *graph){
  //
  // Helper function to correct graphs with symmetric errors 
  // for the bin width
  //
  Double_t binwidth(0.);
  Double_t *ypoints = graph->GetY(),
           *yerrors = graph->GetEY();
  for(int ipt = 0; ipt < graph->GetN(); ipt++){
    binwidth = 2*graph->GetEX()[ipt];
    ypoints[ipt] /= binwidth;
    yerrors[ipt] /= binwidth;
  }
}

//__________________________________________
void AliHFEtools::NormaliseBinWdithAsymm(TGraphAsymmErrors *graph){
  //
  // Helper function to correct graphs with asymmetric errors 
  // for the bin width
  //
  Double_t binwidth(0.);
  Double_t *ypoints = graph->GetY(),
           *yerrorslow = graph->GetEYlow(),
           *yerrorshigh = graph->GetEYhigh();
  for(int ipt = 0; ipt < graph->GetN(); ipt++){
    binwidth = graph->GetEXlow()[ipt] + graph->GetEXhigh()[ipt];
    ypoints[ipt] /= binwidth;
    yerrorslow[ipt] /= binwidth;
    yerrorshigh[ipt] /= binwidth;
  }
}
