/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt                        *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to do analysis on conversion pairs
//---------------------------------------------
////////////////////////////////////////////////

#include "AliGammaConversionHistograms.h"
#include "TMath.h"
#include "TObjString.h"
#include "TMap.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"


using namespace std;

ClassImp(AliGammaConversionHistograms)


AliGammaConversionHistograms::AliGammaConversionHistograms() :
  fHistogramMap(new TMap()),
  fNPhiIndex(0),
  fNRIndex(0),
  fMinRadius(0.),
  fMaxRadius(0.),
  fDeltaR(0.),
  fMinPhi(0.),
  fMaxPhi(0.),
  fDeltaPhi(0.),
  fMappingContainer(NULL)
{
  // see header file for documenation
}


AliGammaConversionHistograms::AliGammaConversionHistograms(const AliGammaConversionHistograms & original) :
  fHistogramMap(original.fHistogramMap),
  fNPhiIndex(original.fNPhiIndex),
  fNRIndex(original.fNRIndex),
  fMinRadius(original.fMinRadius),
  fMaxRadius(original.fMaxRadius),
  fDeltaR(original.fDeltaR),
  fMinPhi(original.fMinPhi),
  fMaxPhi(original.fMaxPhi),
  fDeltaPhi(original.fDeltaPhi),
  fMappingContainer(original.fMappingContainer)
{    
  //see header file for documentation
}


AliGammaConversionHistograms & AliGammaConversionHistograms::operator = (const AliGammaConversionHistograms & /*original*/)
{
  // assignment operator
  return *this;
}


AliGammaConversionHistograms::~AliGammaConversionHistograms() {
  //destructor


}

void AliGammaConversionHistograms::AddHistogram(TString histogramName, TString histogramTitle, Int_t nXBins, Double_t firstX,Double_t lastX,TString xAxisTitle, TString yAxisTitle){
  // see header file for documentation
  TH1F *tmp = new TH1F(histogramName, histogramTitle,nXBins,firstX,lastX);
  tmp->GetXaxis()->SetTitle(xAxisTitle);
  tmp->GetYaxis()->SetTitle(yAxisTitle);
  TObjString* tobjstring = new TObjString(histogramName.Data());
  fHistogramMap->Add((TObject*)tobjstring,(TObject*)tmp);
}

void AliGammaConversionHistograms::AddHistogram(TString histogramName, TString histogramTitle, Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  // see header file for documentation
  TH2F *tmp = new TH2F(histogramName, histogramTitle,nXBins,firstX,lastX,nYBins,firstY,lastY);
  tmp->GetXaxis()->SetTitle(xAxisTitle);
  tmp->GetYaxis()->SetTitle(yAxisTitle);
  TObjString *tobjstring = new TObjString(histogramName.Data());
  fHistogramMap->Add((TObject*)tobjstring,(TObject*)tmp);
}

void AliGammaConversionHistograms::FillHistogram(TString histogramName, Double_t xValue) const{
  //see header file for documentation
  TH1 *tmp = (TH1*)fHistogramMap->GetValue(histogramName.Data());
  if(tmp){
      tmp->Fill(xValue);
  }
  else{
    cout<<"Histogram does not exist: "<<histogramName.Data()<<endl;
  }
}

void AliGammaConversionHistograms::FillHistogram(TString histogramName, Double_t xValue, Double_t yValue) const{
  //see header file for documentation
  TH1 *tmp = (TH1*)fHistogramMap->GetValue(histogramName.Data());
  if(tmp){
    tmp->Fill(xValue, yValue);
  }
  else{
    cout<<"Histogram does not exist: "<<histogramName.Data()<<endl;
  }
}

void AliGammaConversionHistograms::GetOutputContainer(TList *fOutputContainer){
  //checking if the container is alrerady created

  if(fOutputContainer == NULL){
    //print warning
    return;
  }

  if(fHistogramMap != NULL){
    TIter iter(fHistogramMap);
    TObjString *histogramName;
    while ((histogramName = (TObjString*) iter.Next())) {
      TString histogramString = histogramName->GetString();
      if(histogramString.Contains("Mapping")){// means it should be put in the mapping folder
	if(fMappingContainer == NULL){
	  fMappingContainer = new TList();
	  fMappingContainer->SetName("Mapping histograms");
	}
	if(fMappingContainer != NULL){
	  fMappingContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
	}
      }
      fOutputContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      histogramName = NULL;
    } // end while
    //   fOutputContainer->Add(fMappingContainer);
  }
}

Int_t AliGammaConversionHistograms::GetRBin(Double_t radius) const{
  // see header file for documentation
  Int_t iResult=0;
  if(fDeltaR>0){
    iResult = (Int_t)((radius - fMinRadius)/fDeltaR);
  }
  return iResult;
}

Int_t AliGammaConversionHistograms::GetPhiBin(Double_t phi) const{
  // see header file for documentation
  Int_t iResult=0;
  if(fDeltaPhi>0){
    if(phi>TMath::Pi()){
      phi-=2*TMath::Pi();
    }
    iResult = (Int_t)((phi - fMinPhi)/fDeltaPhi);
  }
  return iResult;
}



void AliGammaConversionHistograms::InitializeMappingValues(Int_t nPhiIndex, Int_t nRIndex, Int_t nBinsR, Double_t minRadius, Double_t maxRadius,Int_t nBinsPhi, Double_t minPhi, Double_t maxPhi){
  // Initializing the valuse for the mapping

  fNPhiIndex = nPhiIndex;
  fNRIndex   = nRIndex;
  fMinRadius      = minRadius;
  fMaxRadius      = maxRadius;
  if(nBinsR>0 && nRIndex!=0){
    fDeltaR       = (fMaxRadius - fMinRadius)/nRIndex;
  }
  fMinPhi         = minPhi;
  fMaxPhi         = maxPhi;
  if(nBinsPhi>0 && nPhiIndex!=0){
    fDeltaPhi     = (fMaxPhi-fMinPhi)/nPhiIndex;
  }
}


//mapping
void AliGammaConversionHistograms::AddMappingHistograms(Int_t nPhiIndex, Int_t nRIndex,Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  // see header file for documentation
  
  for(Int_t phi =0; phi<=fNPhiIndex;phi++){

    for(Int_t r =0; r<fNRIndex;r++){

      // setting axis to "" changes below
      xAxisTitle="";
      yAxisTitle="";
      //Creating the axis titles
      if(xAxisTitle.Length() == 0){
	xAxisTitle.Form("Phi %02d",phi);
      }
      
      if(yAxisTitle.Length() == 0){
	yAxisTitle.Form("R %02d",phi);
      }

      //MC
      TString nameMC="";
      nameMC.Form("MC_EP_Mapping-Phi%02d-R%02d",phi,r);
      TString titleMC="";
      titleMC.Form("Electron-Positron MC Mapping-Phi%02d-R%02d",phi,r);

      AddHistogram(nameMC, titleMC, nXBins, firstX, lastX, nYBins, firstY, lastY, xAxisTitle, yAxisTitle);

      //ESD
      TString nameESD="";
      nameESD.Form("ESD_EP_Mapping-Phi%02d-R%02d",phi,r);
      TString titleESD="";
      titleESD.Form("Electron-Positron ESD Mapping-Phi%02d-R%02d",phi,r);

      AddHistogram(nameESD, titleESD, nXBins, firstX, lastX, nYBins, firstY, lastY, xAxisTitle, yAxisTitle);
    }
  }


  for(Int_t phi =0; phi<=nPhiIndex;phi++){ 

    // setting axis to "" changes below
    xAxisTitle="";
    yAxisTitle="";
    //Creating the axis titles
    if(xAxisTitle.Length() == 0){
      xAxisTitle.Form("Phi %02d",phi);
    }
    if(yAxisTitle.Length() == 0){
      yAxisTitle = "Counts";
    }
    
    //MC
    TString nameMC="";
    nameMC.Form("MC_EP_Mapping-Phi%02d",phi);
    TString titleMC="";
    titleMC.Form("Electron-Positron MC Mapping-Phi%02d",phi);
    
    AddHistogram(nameMC, titleMC, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);

    //MC
    TString nameESD="";
    nameESD.Form("ESD_EP_Mapping-Phi%02d",phi);
    TString titleESD="";
    titleESD.Form("Electron-Positron ESD Mapping-Phi%02d",phi);
    
    AddHistogram(nameESD, titleESD, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
  }


  for(Int_t r =0; r<=nRIndex;r++){

    // setting axis to "" changes below
    xAxisTitle="";
    yAxisTitle="";
    //Creating the axis titles
    if(xAxisTitle.Length() == 0){
      xAxisTitle.Form("R %02d",r);
    }
    if(yAxisTitle.Length() == 0){
      yAxisTitle = "Counts";
    }
    
    //MC
    TString nameMC="";
    nameMC.Form("MC_EP_Mapping-R%02d",r);
    TString titleMC="";
    titleMC.Form("Electron-Positron MC Mapping-R%02d",r);
    
    AddHistogram(nameMC, titleMC, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);

    //ESD
    TString nameESD="";
    nameESD.Form("ESD_EP_Mapping-R%02d",r);
    TString titleESD="";
    titleESD.Form("Electron-Positron ESD Mapping-R%02d",r);
    
    AddHistogram(nameESD, titleESD, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);

    //Mapping Phi in R
    TString nameMCPhiInR="";
    nameMCPhiInR.Form("MC_EP_Mapping_Phi_R-%02d",r);
    TString titleMCPhiInR="";
    titleMCPhiInR.Form("MC Mapping of Phi in R%02d",r);
    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);

    //Mapping Phi in R
    TString nameESDPhiInR="";
    nameESDPhiInR.Form("ESD_EP_Mapping_Phi_R-%02d",r);
    TString titleESDPhiInR="";
    titleESDPhiInR.Form("ESD Mapping of Phi in R%02d",r);
    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
  }
}
