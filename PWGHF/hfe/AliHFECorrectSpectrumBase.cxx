
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
// Class for spectrum correction
// Subtraction of hadronic background, Unfolding of the data and
// Renormalization done here
// The following containers have to be set:
//  - Correction framework container for real data
//  - Correction framework container for MC (Efficiency Map)
//  - Correction framework container for background coming from data
//  - Correction framework container for background coming from MC
//
//  Author: 
//            Raphaelle Bailhache <R.Bailhache@gsi.de>
//            Markus Fasel <M.Fasel@gsi.de>
//

#include <TArrayD.h>
#include <TH1.h>
#include <TList.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TPad.h>
#include <TH2D.h>
#include <TF1.h>

#include "AliPID.h"
#include "AliCFContainer.h"
#include "AliCFDataGrid.h"
#include "AliCFEffGrid.h"
#include "AliCFGridSparse.h"
#include "AliCFUnfolding.h"
#include "AliLog.h"

#include "AliHFECorrectSpectrumBase.h"
#include "AliHFEcuts.h"
#include "AliHFEcontainer.h"
#include "AliHFEtools.h"

ClassImp(AliHFECorrectSpectrumBase)

//____________________________________________________________
AliHFECorrectSpectrumBase::AliHFECorrectSpectrumBase(const char *name):
  TNamed(name, ""),
  fCFContainers(new TObjArray(kNbCFContainers)),
  fCorrelation(NULL),
  fEfficiencyFunction(NULL),
  fEtaSelected(kFALSE),
  fSetSmoothing(kFALSE),
  fNbDimensions(1),
  fNEvents(0),
  fStepMC(-1),
  fStepTrue(-1),
  fStepData(-1),
  fStepBeforeCutsV0(-1),
  fStepAfterCutsV0(-1),
  fStepGuessedUnfolding(-1),
  fNumberOfIterations(10),
  fChargeChoosen(kAllCharge),
  fTestCentralityLow(-1),
  fTestCentralityHigh(-1)
{
  //
  // Default constructor
  //

  memset(fEtaRange, 0, sizeof(Double_t) * 2);
  memset(fEtaRangeNorm, 0, sizeof(Double_t) * 2);
 
}
//____________________________________________________________
AliHFECorrectSpectrumBase::AliHFECorrectSpectrumBase(const AliHFECorrectSpectrumBase &ref):
  TNamed(ref),
  fCFContainers(NULL),
  fCorrelation(NULL),
  fEfficiencyFunction(NULL),
  fEtaSelected(ref.fEtaSelected),
  fSetSmoothing(ref.fSetSmoothing),
  fNbDimensions(ref.fNbDimensions),
  fNEvents(ref.fNEvents),
  fStepMC(ref.fStepMC),
  fStepTrue(ref.fStepTrue),
  fStepData(ref.fStepData),
  fStepBeforeCutsV0(ref.fStepBeforeCutsV0),
  fStepAfterCutsV0(ref.fStepAfterCutsV0),
  fStepGuessedUnfolding(ref.fStepGuessedUnfolding),
  fNumberOfIterations(ref.fNumberOfIterations),
  fChargeChoosen(ref.fChargeChoosen),
  fTestCentralityLow(ref.fTestCentralityLow),
  fTestCentralityHigh(ref.fTestCentralityHigh)
{
  //
  // Copy constructor
  //
  ref.Copy(*this);

}
//____________________________________________________________
AliHFECorrectSpectrumBase &AliHFECorrectSpectrumBase::operator=(const AliHFECorrectSpectrumBase &ref){
  //
  // Assignment operator
  //
  if(this == &ref) 
    ref.Copy(*this);
  return *this;
}
//____________________________________________________________
void AliHFECorrectSpectrumBase::Copy(TObject &o) const {
  // 
  // Copy into object o
  //
  AliHFECorrectSpectrumBase &target = dynamic_cast<AliHFECorrectSpectrumBase &>(o);
  target.fCFContainers = fCFContainers;
  target.fCorrelation = fCorrelation;
  target.fEfficiencyFunction = fEfficiencyFunction;
  target.fEtaSelected = fEtaSelected;
  target.fSetSmoothing = fSetSmoothing;
  target.fNbDimensions = fNbDimensions;
  target.fNEvents = fNEvents;
  target.fStepMC = fStepMC;
  target.fStepTrue = fStepTrue;
  target.fStepData = fStepData;
  target.fStepBeforeCutsV0 = fStepBeforeCutsV0;
  target.fStepAfterCutsV0 = fStepAfterCutsV0;
  target.fStepGuessedUnfolding = fStepGuessedUnfolding;
  target.fNumberOfIterations = fNumberOfIterations;
  target.fChargeChoosen = fChargeChoosen;
  target.fTestCentralityLow = fTestCentralityLow;
  target.fTestCentralityHigh = fTestCentralityHigh;
  target.fEtaRange[0] = fEtaRange[0];
  target.fEtaRange[1] = fEtaRange[1];
  target.fEtaRangeNorm[0] = fEtaRangeNorm[0];
  target.fEtaRangeNorm[1] = fEtaRangeNorm[1];

}

//____________________________________________________________
AliHFECorrectSpectrumBase::~AliHFECorrectSpectrumBase(){
  //
  // Destructor
  //
  if(fCFContainers) delete fCFContainers;
 
}
//__________________________________________________________________________________
TGraphErrors *AliHFECorrectSpectrumBase::Normalize(THnSparse * const spectrum) const {
  //
  // Normalize the spectrum to 1/(2*Pi*p_{T})*dN/dp_{T} (GeV/c)^{-2}
  // Give the final pt spectrum to be compared
  //

  if(fNEvents > 0) {

    TH1D* projection = spectrum->Projection(0);
    CorrectFromTheWidth(projection);
    TGraphErrors *graphError = NormalizeTH1(projection);
    return graphError;
  
  }
    
  return 0x0;
  

}
//__________________________________________________________________________________
TGraphErrors *AliHFECorrectSpectrumBase::Normalize(AliCFDataGrid * const spectrum) const {
  //
  // Normalize the spectrum to 1/(2*Pi*p_{T})*dN/dp_{T} (GeV/c)^{-2}
  // Give the final pt spectrum to be compared
  //
  if(fNEvents > 0) {

    TH1D* projection = (TH1D *) spectrum->Project(0);
    CorrectFromTheWidth(projection);
    TGraphErrors *graphError = NormalizeTH1(projection);

    return graphError;
    
  }

  return 0x0;
  

}
//__________________________________________________________________________________
TGraphErrors *AliHFECorrectSpectrumBase::NormalizeTH1(TH1 *input) const {
  //
  // Normalize the spectrum to 1/(2*Pi*p_{T})*dN/dp_{T} (GeV/c)^{-2}
  // Give the final pt spectrum to be compared
  //
  Double_t chargecoefficient = 0.5;
  if(fChargeChoosen != 0) chargecoefficient = 1.0;

  Double_t etarange = fEtaSelected ? fEtaRangeNorm[1] - fEtaRangeNorm[0] : 1.6;
  printf("Normalizing Eta Range %f\n", etarange);
  if(fNEvents > 0) {

    TGraphErrors *spectrumNormalized = new TGraphErrors(input->GetNbinsX());
    Double_t p = 0, dp = 0; Int_t point = 1;
    Double_t n = 0, dN = 0;
    Double_t nCorr = 0, dNcorr = 0;
    Double_t errdN = 0, errdp = 0;
    for(Int_t ibin = input->GetXaxis()->GetFirst(); ibin <= input->GetXaxis()->GetLast(); ibin++){
      point = ibin - input->GetXaxis()->GetFirst();
      p = input->GetXaxis()->GetBinCenter(ibin);
      dp = input->GetXaxis()->GetBinWidth(ibin)/2.;
      n = input->GetBinContent(ibin);
      dN = input->GetBinError(ibin);
      // New point
      nCorr = chargecoefficient * 1./etarange * 1./(Double_t)(fNEvents) * 1./(2. * TMath::Pi() * p) * n;
      errdN = 1./(2. * TMath::Pi() * p);
      errdp = 1./(2. * TMath::Pi() * p*p) * n;
      dNcorr = chargecoefficient * 1./etarange * 1./(Double_t)(fNEvents) * TMath::Sqrt(errdN * errdN * dN *dN + errdp *errdp * dp *dp);
      
      spectrumNormalized->SetPoint(point, p, nCorr);
      spectrumNormalized->SetPointError(point, dp, dNcorr);
    }
    spectrumNormalized->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    spectrumNormalized->GetYaxis()->SetTitle("#frac{1}{2 #pi p_{T}} #frac{dN}{dp_{T}} / [GeV/c]^{-2}");
    spectrumNormalized->SetMarkerStyle(22);
    spectrumNormalized->SetMarkerColor(kBlue);
    spectrumNormalized->SetLineColor(kBlue);

    return spectrumNormalized;
    
  }

  return 0x0;
  

}
//____________________________________________________________
void AliHFECorrectSpectrumBase::SetContainer(AliCFContainer *cont, AliHFECorrectSpectrumBase::CFContainer_t type){
  //
  // Set the container for a given step to the 
  //
  if(!fCFContainers) fCFContainers = new TObjArray(kNbCFContainers);
  fCFContainers->AddAt(cont, type);
}

//____________________________________________________________
AliCFContainer *AliHFECorrectSpectrumBase::GetContainer(AliHFECorrectSpectrumBase::CFContainer_t type){
  //
  // Get Correction Framework Container for given type
  //
  if(!fCFContainers) return NULL;
  return dynamic_cast<AliCFContainer *>(fCFContainers->At(type));
}
//____________________________________________________________
AliCFContainer *AliHFECorrectSpectrumBase::GetSlicedContainer(AliCFContainer *container, Int_t nDim, Int_t *dimensions,Int_t source,Chargetype_t charge, Int_t centralitylow, Int_t centralityhigh) {
  //
  // Slice bin for a given source of electron
  // nDim is the number of dimension the corrections are done
  // dimensions are the definition of the dimensions
  // source is if we want to keep only one MC source (-1 means we don't cut on the MC source)
  // positivenegative if we want to keep positive (1) or negative (0) or both (-1)
  // centrality (-1 means we do not cut on centrality)
  //
  
  Double_t *varMin = new Double_t[container->GetNVar()],
           *varMax = new Double_t[container->GetNVar()];

  Double_t *binLimits;
  for(Int_t ivar = 0; ivar < container->GetNVar(); ivar++){
    
    binLimits = new Double_t[container->GetNBins(ivar)+1];
    container->GetBinLimits(ivar,binLimits);
    varMin[ivar] = binLimits[0];
    varMax[ivar] = binLimits[container->GetNBins(ivar)];
    // source
    if(ivar == 4){
      if((source>= 0) && (source<container->GetNBins(ivar))) {
	      varMin[ivar] = binLimits[source];
	      varMax[ivar] = binLimits[source];
      }     
    }
    // charge
    if(ivar == 3) {
      if(charge != kAllCharge) varMin[ivar] = varMax[ivar] = charge;
    }
    // eta
    if(ivar == 1){
      for(Int_t ic = 1; ic <= container->GetAxis(1,0)->GetLast(); ic++) 
        AliDebug(1, Form("eta bin %d, min %f, max %f\n", ic, container->GetAxis(1,0)->GetBinLowEdge(ic), container->GetAxis(1,0)->GetBinUpEdge(ic))); 
      if(fEtaSelected){
        varMin[ivar] = fEtaRange[0];
        varMax[ivar] = fEtaRange[1];
      }
    }
    if(fEtaSelected){
      fEtaRangeNorm[0] = container->GetAxis(1,0)->GetBinLowEdge(container->GetAxis(1,0)->FindBin(fEtaRange[0]));
      fEtaRangeNorm[1] = container->GetAxis(1,0)->GetBinUpEdge(container->GetAxis(1,0)->FindBin(fEtaRange[1]));
      AliInfo(Form("Normalization done in eta range [%f,%f]\n", fEtaRangeNorm[0], fEtaRangeNorm[0]));
    }
    // centrality
    if(ivar == 5){
	if((centralitylow>= 0) && (centralitylow<container->GetNBins(ivar)) && (centralityhigh>= 0) && (centralityhigh<container->GetNBins(ivar))) {
	    varMin[ivar] = binLimits[centralitylow];
	    varMax[ivar] = binLimits[centralityhigh];

	    TAxis *axistest = container->GetAxis(5,0);
	    AliDebug(1, Form("Number of bin in centrality direction %d\n",axistest->GetNbins()));
	    AliDebug(1, Form("Project from %f to %f\n",binLimits[centralitylow],binLimits[centralityhigh]));
	    Double_t lowcentrality = axistest->GetBinLowEdge(axistest->FindBin(binLimits[centralitylow]));
	    Double_t highcentrality = axistest->GetBinUpEdge(axistest->FindBin(binLimits[centralityhigh]));
	    AliDebug(1, Form("Low centrality %f and high centrality %f\n",lowcentrality,highcentrality));
	
	}
    }
    
    delete[] binLimits;
    
  }
  
  AliCFContainer *k = container->MakeSlice(nDim, dimensions, varMin, varMax);
  delete[] varMin; delete[] varMax;

  return k;

}

//_________________________________________________________________________
THnSparseF *AliHFECorrectSpectrumBase::GetSlicedCorrelation(THnSparseF *correlationmatrix, Int_t nDim, Int_t *dimensions,Chargetype_t charge,Int_t centralitylow, Int_t centralityhigh) const {
  //
  // Slice correlation
  //

  Int_t ndimensions = correlationmatrix->GetNdimensions();
  //printf("Number of dimension %d correlation map\n",ndimensions);
  if(ndimensions < (2*nDim)) {
    AliError("Problem in the dimensions");
    return NULL;
  }
  
  // Cut in centrality is centrality > -1
  if((5+((Int_t)(ndimensions/2.))) < ndimensions) {
    if((centralitylow >=0) && (centralityhigh >=0)) {
      
      TAxis *axiscentrality0 = correlationmatrix->GetAxis(5);
      TAxis *axiscentrality1 = correlationmatrix->GetAxis(5+((Int_t)(ndimensions/2.)));
      
      Int_t bins0 = axiscentrality0->GetNbins();
      Int_t bins1 = axiscentrality1->GetNbins();
      
      AliDebug(1, Form("Number of centrality bins: %d and %d\n",bins0,bins1));
      if(bins0 != bins1) {
	AliError("Problem in the dimensions");
	return NULL;
      }
      
      if((centralitylow>= 0) && (centralitylow<bins0) && (centralityhigh>= 0) && (centralityhigh<bins0)) {
	axiscentrality0->SetRangeUser(centralitylow,centralityhigh);
	axiscentrality1->SetRangeUser(centralitylow,centralityhigh);
	
	Double_t lowcentrality0 = axiscentrality0->GetBinLowEdge(axiscentrality0->FindBin(centralitylow));
	Double_t highcentrality0 = axiscentrality0->GetBinUpEdge(axiscentrality0->FindBin(centralityhigh));
	Double_t lowcentrality1 = axiscentrality1->GetBinLowEdge(axiscentrality1->FindBin(centralitylow));
	Double_t highcentrality1 = axiscentrality1->GetBinUpEdge(axiscentrality1->FindBin(centralityhigh));
	AliDebug(1,Form("0 Low centrality %f and high centrality %f\n",lowcentrality0,highcentrality0));
	AliDebug(1,Form("1 Low centrality %f and high centrality %f\n",lowcentrality1,highcentrality1));
	
      }    
    }
  }

  // Cut in eta > -1
  if(fEtaSelected){
    if((1+((Int_t)(ndimensions/2.))) < ndimensions) {
      
      TAxis *axiseta0 = correlationmatrix->GetAxis(1);
      TAxis *axiseta1 = correlationmatrix->GetAxis(1+((Int_t)(ndimensions/2.)));
      
      Int_t bins0 = axiseta0->GetNbins();
      Int_t bins1 = axiseta1->GetNbins();
      
      AliDebug(1, Form("Number of eta bins: %d and %d\n",bins0,bins1));
      if(bins0 != bins1) {
	AliError("Problem in the dimensions");
	return NULL;
      }
      
      axiseta0->SetRangeUser(fEtaRange[0],fEtaRange[1]);
      axiseta1->SetRangeUser(fEtaRange[0],fEtaRange[1]);
	
	Double_t loweta0 = axiseta0->GetBinLowEdge(axiseta0->FindBin(fEtaRange[0]));
	Double_t higheta0 = axiseta0->GetBinUpEdge(axiseta0->FindBin(fEtaRange[1]));
	Double_t loweta1 = axiseta1->GetBinLowEdge(axiseta1->FindBin(fEtaRange[0]));
	Double_t higheta1 = axiseta1->GetBinUpEdge(axiseta1->FindBin(fEtaRange[1]));
	AliInfo(Form("0 Low eta %f and high eta %f\n",loweta0,higheta0));
	AliInfo(Form("1 Low eta %f and high eta %f\n",loweta1,higheta1));
	
    }    
  }

  // Cut in charge
  if(charge != kAllCharge) {
    if((3+((Int_t)(ndimensions/2.))) < ndimensions) {
      
      TAxis *axischarge0 = correlationmatrix->GetAxis(3);
      TAxis *axischarge1 = correlationmatrix->GetAxis(3+((Int_t)(ndimensions/2.)));
      
      Int_t bins0 = axischarge0->GetNbins();
      Int_t bins1 = axischarge1->GetNbins();
      
      AliDebug(1, Form("Number of charge bins: %d and %d\n",bins0,bins1));
      if(bins0 != bins1) {
	AliError("Problem in the dimensions");
	return NULL;
      }
      
      axischarge0->SetRangeUser(charge,charge);
      axischarge1->SetRangeUser(charge,charge);
      
      Double_t lowcharge0 = axischarge0->GetBinLowEdge(axischarge0->FindBin(charge));
      Double_t highcharge0 = axischarge0->GetBinUpEdge(axischarge0->FindBin(charge));
      Double_t lowcharge1 = axischarge1->GetBinLowEdge(axischarge1->FindBin(charge));
      Double_t highcharge1 = axischarge1->GetBinUpEdge(axischarge1->FindBin(charge));
      AliInfo(Form("0 Low charge %f and high charge %f\n",lowcharge0,highcharge0));
      AliInfo(Form("1 Low charge %f and high charge %f\n",lowcharge1,highcharge1));
      
    }
  }
  

  Int_t ndimensionsContainer = (Int_t) ndimensions/2;
  
  Int_t *dim = new Int_t[nDim*2];
  for(Int_t iter=0; iter < nDim; iter++){
    dim[iter] = dimensions[iter];
    dim[iter+nDim] = ndimensionsContainer + dimensions[iter];
  }
    
  THnSparseF *k = (THnSparseF *) correlationmatrix->Projection(nDim*2,dim);

  delete[] dim; 
  return k;
  
}
//___________________________________________________________________________
void AliHFECorrectSpectrumBase::CorrectFromTheWidth(TH1D *h1) const {
  //
  // Correct from the width of the bins --> dN/dp_{T} (GeV/c)^{-1}
  //

  TAxis *axis = h1->GetXaxis();
  Int_t nbinX = h1->GetNbinsX();

  for(Int_t i = 1; i <= nbinX; i++) {

    Double_t width = axis->GetBinWidth(i);
    Double_t content = h1->GetBinContent(i);
    Double_t error = h1->GetBinError(i); 
    h1->SetBinContent(i,content/width); 
    h1->SetBinError(i,error/width);
  }

}

//___________________________________________________________________________
void AliHFECorrectSpectrumBase::CorrectStatErr(AliCFDataGrid *backgroundGrid) const { 
  //
  // Correct statistical error
  //

  TH1D *h1 = (TH1D*)backgroundGrid->Project(0);
  Int_t nbinX = h1->GetNbinsX();
  Int_t bins[1];
  for(Long_t i = 1; i <= nbinX; i++) {
    bins[0] = i;
    Float_t content = h1->GetBinContent(i);
    if(content>0){
      Float_t error = TMath::Sqrt(content);
      backgroundGrid->SetElementError(bins, error);
    }
  }
}
//_________________________________________________________________________
TObject* AliHFECorrectSpectrumBase::GetSpectrum(const AliCFContainer * const c, Int_t step) {
  AliCFDataGrid* data = new AliCFDataGrid("data","",*c, step);
  return data;
}
//_________________________________________________________________________
TObject* AliHFECorrectSpectrumBase::GetEfficiency(const AliCFContainer * const c, Int_t step, Int_t step0){
  // 
  // Create efficiency grid and calculate efficiency
  // of step to step0
  //
  TString name("eff");
  name += step;
  name+= step0;
  AliCFEffGrid* eff = new AliCFEffGrid((const char*)name,"",*c);
  eff->CalculateEfficiency(step,step0);
  return eff;
}
