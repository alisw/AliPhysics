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
  fDeltaPhi(0.)
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
  fDeltaPhi(original.fDeltaPhi)
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
    cout<<"Histogram does not exist"<<histogramName.Data()<<endl;
  }
}

void AliGammaConversionHistograms::FillHistogram(TString histogramName, Double_t xValue, Double_t yValue) const{
  //see header file for documentation
  TH1 *tmp = (TH1*)fHistogramMap->GetValue(histogramName.Data());
  if(tmp){
    tmp->Fill(xValue, yValue);
  }
  else{
    cout<<"Histogram does not exist"<<histogramName.Data()<<endl;
  }
}

void AliGammaConversionHistograms::GetOutputContainer(TList *fOutputContainer) const{
  //checking if the container is alrerady created

  if(fOutputContainer == NULL){
    //print warning
    return;
  }
  cout<<"Creating the histogram output container"<<endl;

  if(fHistogramMap){
    TIter iter(fHistogramMap);
    TObjString *histogramName;
    while ((histogramName = (TObjString*) iter.Next())) {
      cout<<"Histohram name "<<histogramName->GetString().Data()<<endl;
      TString histogramString = histogramName->GetString();
      fOutputContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      histogramName = NULL;
    }  
  }

  //remember mapping stuff!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  /*
  TList*  fMappingContainer = new TList();
  fMappingContainer->SetName("Mapping Histograms");

  if(fMCEPR != NULL){ fOutputContainer->Add(fMCEPR);}
  if(fMCEPZR != NULL){ fOutputContainer->Add(fMCEPZR);}
  if(fMCEPXY != NULL){ fOutputContainer->Add(fMCEPXY);}
  if(fMCEPOpeningAngle != NULL){ fOutputContainer->Add(fMCEPOpeningAngle);}

  if(fMCEEnergy != NULL){ fOutputContainer->Add(fMCEEnergy);}
  if(fMCEPt != NULL){ fOutputContainer->Add(fMCEPt);}
  if(fMCEEta != NULL){ fOutputContainer->Add(fMCEEta);}
  if(fMCEPhi != NULL){ fOutputContainer->Add(fMCEPhi);}

  if(fMCPEnergy != NULL){ fOutputContainer->Add(fMCPEnergy);}
  if(fMCPPt != NULL){ fOutputContainer->Add(fMCPPt);}
  if(fMCPEta != NULL){ fOutputContainer->Add(fMCPEta);}
  if(fMCPPhi != NULL){ fOutputContainer->Add(fMCPPhi);}

  if(fMCGammaEnergy != NULL){ fOutputContainer->Add(fMCGammaEnergy);}
  if(fMCGammaPt != NULL){ fOutputContainer->Add(fMCGammaPt);}
  if(fMCGammaEta != NULL){ fOutputContainer->Add(fMCGammaEta);}
  if(fMCGammaPhi != NULL){ fOutputContainer->Add(fMCGammaPhi);}

  if(fMCDirectGammaEnergy != NULL){ fOutputContainer->Add(fMCDirectGammaEnergy);}
  if(fMCDirectGammaPt != NULL){ fOutputContainer->Add(fMCDirectGammaPt);}
  if(fMCDirectGammaEta != NULL){ fOutputContainer->Add(fMCDirectGammaEta);}
  if(fMCDirectGammaPhi != NULL){ fOutputContainer->Add(fMCDirectGammaPhi);}

  //mapping
  for(UInt_t i=0;i<fMCMapping.size();i++){
    for(UInt_t j=0;j<fMCMapping[i].size();j++){
      if(fMCMapping[i][j] != NULL){fMappingContainer->Add(fMCMapping[i][j]);}
    }
  }
  for(UInt_t i=0;i<fMCMappingPhi.size();i++){
    if(fMCMappingPhi[i] != NULL){fMappingContainer->Add(fMCMappingPhi[i]);}
  }
  for(UInt_t i=0;i<fMCMappingR.size();i++){
    if(fMCMappingR[i] != NULL){fMappingContainer->Add(fMCMappingR[i]);}
  }
  if(fMCMatchGammaEta != NULL){ fOutputContainer->Add(fMCMatchGammaEta);}
  if(fMCMatchGammaPhi != NULL){ fOutputContainer->Add(fMCMatchGammaPhi);}
  if(fMCMatchGammaPt != NULL){ fOutputContainer->Add(fMCMatchGammaPt);}
  if(fMCMatchGammaEnergy != NULL){ fOutputContainer->Add(fMCMatchGammaEnergy);}
  if(fMCMatchGammaMass != NULL){ fOutputContainer->Add(fMCMatchGammaMass);}
  if(fMCMatchGammaOpeningAngle != NULL){ fOutputContainer->Add(fMCMatchGammaOpeningAngle);}
  if(fMCMatchGammaR != NULL){ fOutputContainer->Add(fMCMatchGammaR);}
  if(fMCMatchGammaZR != NULL){ fOutputContainer->Add(fMCMatchGammaZR);}
  if(fMCMatchGammaXY != NULL){ fOutputContainer->Add(fMCMatchGammaXY);}

  if(fMCPi0Eta != NULL){ fOutputContainer->Add(fMCPi0Eta);}
  if(fMCPi0Phi != NULL){ fOutputContainer->Add(fMCPi0Phi);}
  if(fMCPi0Pt != NULL){ fOutputContainer->Add(fMCPi0Pt);}
  if(fMCPi0Energy != NULL){ fOutputContainer->Add(fMCPi0Energy);}
  if(fMCPi0Mass != NULL){ fOutputContainer->Add(fMCPi0Mass);}
  if(fMCPi0OpeningAngleGamma != NULL){ fOutputContainer->Add(fMCPi0OpeningAngleGamma);}
  if(fMCPi0R != NULL){ fOutputContainer->Add(fMCPi0R);}
  if(fMCPi0ZR != NULL){ fOutputContainer->Add(fMCPi0ZR);}
  if(fMCPi0XY != NULL){ fOutputContainer->Add(fMCPi0XY);}
  if(fMCPi0SecondariesXY != NULL){ fOutputContainer->Add(fMCPi0SecondariesXY);}

  if(fMCEtaEta != NULL){ fOutputContainer->Add(fMCEtaEta);}
  if(fMCEtaPhi != NULL){ fOutputContainer->Add(fMCEtaPhi);}
  if(fMCEtaPt != NULL){ fOutputContainer->Add(fMCEtaPt);}
  if(fMCEtaEnergy != NULL){ fOutputContainer->Add(fMCEtaEnergy);}
  if(fMCEtaMass != NULL){ fOutputContainer->Add(fMCEtaMass);}
  if(fMCEtaOpeningAngleGamma != NULL){ fOutputContainer->Add(fMCEtaOpeningAngleGamma);}
  if(fMCEtaR != NULL){ fOutputContainer->Add(fMCEtaR);}
  if(fMCEtaZR != NULL){ fOutputContainer->Add(fMCEtaZR);}
  if(fMCEtaXY != NULL){ fOutputContainer->Add(fMCEtaXY);}
    
  // Histograms from esd tracks
  if(fESDEPR != NULL){ fOutputContainer->Add(fESDEPR);}
  if(fESDEPZR != NULL){ fOutputContainer->Add(fESDEPZR);}
  if(fESDEPXY != NULL){ fOutputContainer->Add(fESDEPXY);}
  if(fESDEPOpeningAngle != NULL){ fOutputContainer->Add(fESDEPOpeningAngle);}

  if(fESDEEnergy != NULL){ fOutputContainer->Add(fESDEEnergy);}
  if(fESDEPt != NULL){ fOutputContainer->Add(fESDEPt);}
  if(fESDEEta != NULL){ fOutputContainer->Add(fESDEEta);}
  if(fESDEPhi != NULL){ fOutputContainer->Add(fESDEPhi);}

  if(fESDPEnergy != NULL){ fOutputContainer->Add(fESDPEnergy);}
  if(fESDPPt != NULL){ fOutputContainer->Add(fESDPPt);}
  if(fESDPEta != NULL){ fOutputContainer->Add(fESDPEta);}
  if(fESDPPhi != NULL){ fOutputContainer->Add(fESDPPhi);}

  if(fESDGammaEnergy != NULL){ fOutputContainer->Add(fESDGammaEnergy);}
  if(fESDGammaPt != NULL){ fOutputContainer->Add(fESDGammaPt);}
  if(fESDGammaEta != NULL){ fOutputContainer->Add(fESDGammaEta);}
  if(fESDGammaPhi != NULL){ fOutputContainer->Add(fESDGammaPhi);}

  //mapping
  for(UInt_t i=0;i<fESDMapping.size();i++){
    for(UInt_t j=0;j<fESDMapping[i].size();j++){
      if(fESDMapping[i][j] != NULL){fMappingContainer->Add(fESDMapping[i][j]);}
    }
  }
  for(UInt_t i=0;i<fESDMappingPhi.size();i++){
    if(fESDMappingPhi[i] != NULL){fMappingContainer->Add(fESDMappingPhi[i]);}
  }
  for(UInt_t i=0;i<fESDMappingR.size();i++){
    if(fESDMappingR[i] != NULL){fMappingContainer->Add(fESDMappingR[i]);}
  }

  fOutputContainer->Add(fMappingContainer);

  if(fESDMatchGammaOpeningAngle != NULL){ fOutputContainer->Add(fESDMatchGammaOpeningAngle);}
  if(fESDMatchGammaEnergy != NULL){ fOutputContainer->Add(fESDMatchGammaEnergy);}
  if(fESDMatchGammaPt != NULL){ fOutputContainer->Add(fESDMatchGammaPt);}
  if(fESDMatchGammaEta != NULL){ fOutputContainer->Add(fESDMatchGammaEta);}
  if(fESDMatchGammaPhi != NULL){ fOutputContainer->Add(fESDMatchGammaPhi);}
  if(fESDMatchGammaMass != NULL){ fOutputContainer->Add(fESDMatchGammaMass);}
  if(fESDMatchGammaWidth != NULL){ fOutputContainer->Add(fESDMatchGammaWidth);}
  if(fESDMatchGammaChi2 != NULL){ fOutputContainer->Add(fESDMatchGammaChi2);}
  if(fESDMatchGammaNDF != NULL){ fOutputContainer->Add(fESDMatchGammaNDF);}
  if(fESDMatchGammaR != NULL){ fOutputContainer->Add(fESDMatchGammaR);}
  if(fESDMatchGammaZR != NULL){ fOutputContainer->Add(fESDMatchGammaZR);}
  if(fESDMatchGammaXY != NULL){ fOutputContainer->Add(fESDMatchGammaXY);}

  if(fESDPi0OpeningAngleGamma != NULL){ fOutputContainer->Add(fESDPi0OpeningAngleGamma);}
  if(fESDPi0Energy != NULL){ fOutputContainer->Add(fESDPi0Energy);}
  if(fESDPi0Pt != NULL){ fOutputContainer->Add(fESDPi0Pt);}
  if(fESDPi0Eta != NULL){ fOutputContainer->Add(fESDPi0Eta);}
  if(fESDPi0Phi != NULL){ fOutputContainer->Add(fESDPi0Phi);}
  if(fESDPi0Mass != NULL){ fOutputContainer->Add(fESDPi0Mass);}
  if(fESDPi0R != NULL){ fOutputContainer->Add(fESDPi0R);}
  if(fESDPi0ZR != NULL){ fOutputContainer->Add(fESDPi0ZR);}
  if(fESDPi0XY != NULL){ fOutputContainer->Add(fESDPi0XY);}

  if(fESDEtaOpeningAngleGamma != NULL){ fOutputContainer->Add(fESDEtaOpeningAngleGamma);}
  if(fESDEtaEnergy != NULL){ fOutputContainer->Add(fESDEtaEnergy);}
  if(fESDEtaPt != NULL){ fOutputContainer->Add(fESDEtaPt);}
  if(fESDEtaEta != NULL){ fOutputContainer->Add(fESDEtaEta);}
  if(fESDEtaPhi != NULL){ fOutputContainer->Add(fESDEtaPhi);}
  if(fESDEtaMass != NULL){ fOutputContainer->Add(fESDEtaMass);}
  if(fESDEtaR != NULL){ fOutputContainer->Add(fESDEtaR);}
  if(fESDEtaZR != NULL){ fOutputContainer->Add(fESDEtaZR);}
  if(fESDEtaXY != NULL){ fOutputContainer->Add(fESDEtaXY);}

  if(fESDBackgroundOpeningAngleGamma != NULL){ fOutputContainer->Add(fESDBackgroundOpeningAngleGamma);}
  if(fESDBackgroundEnergy != NULL){ fOutputContainer->Add(fESDBackgroundEnergy);}
  if(fESDBackgroundPt != NULL){ fOutputContainer->Add(fESDBackgroundPt);}
  if(fESDBackgroundEta != NULL){ fOutputContainer->Add(fESDBackgroundEta);}
  if(fESDBackgroundPhi != NULL){ fOutputContainer->Add(fESDBackgroundPhi);}
  if(fESDBackgroundMass != NULL){ fOutputContainer->Add(fESDBackgroundMass);}
  if(fESDBackgroundR != NULL){ fOutputContainer->Add(fESDBackgroundR);}
  if(fESDBackgroundZR != NULL){ fOutputContainer->Add(fESDBackgroundZR);}
  if(fESDBackgroundXY != NULL){ fOutputContainer->Add(fESDBackgroundXY);}

  if(fResolutiondPt != NULL){ fOutputContainer->Add(fResolutiondPt);}
  if(fResolutiondR != NULL){ fOutputContainer->Add(fResolutiondR);}
  if(fResolutiondZ != NULL){ fOutputContainer->Add(fResolutiondZ);}
  if(fResolutiondRdPt != NULL){ fOutputContainer->Add(fResolutiondRdPt);}
  if(fResolutionMCPt != NULL){ fOutputContainer->Add(fResolutionMCPt);}
  if(fResolutionMCR != NULL){ fOutputContainer->Add(fResolutionMCR);}
  if(fResolutionMCZ != NULL){ fOutputContainer->Add(fResolutionMCZ);}
  if(fResolutionESDPt != NULL){ fOutputContainer->Add(fResolutionESDPt);}
  if(fResolutionESDR != NULL){ fOutputContainer->Add(fResolutionESDR);}
  if(fResolutionESDZ != NULL){ fOutputContainer->Add(fResolutionESDZ);}

  if(fNumberOfV0s != NULL){fOutputContainer->Add(fNumberOfV0s);}
  if(fNumberOfSurvivingV0s != NULL){fOutputContainer->Add(fNumberOfSurvivingV0s);}
  if(fV0MassDebugCut1 != NULL){fOutputContainer->Add(fV0MassDebugCut1);}
  if(fV0MassDebugCut2 != NULL){fOutputContainer->Add(fV0MassDebugCut2);}
  if(fV0MassDebugCut3 != NULL){fOutputContainer->Add(fV0MassDebugCut3);}
  if(fV0MassDebugCut4 != NULL){fOutputContainer->Add(fV0MassDebugCut4);}
  if(fV0MassDebugCut5 != NULL){fOutputContainer->Add(fV0MassDebugCut5);}
  if(fV0MassDebugCut6 != NULL){fOutputContainer->Add(fV0MassDebugCut6);}
  if(fV0MassDebugCut7 != NULL){fOutputContainer->Add(fV0MassDebugCut7);}
  if(fV0MassDebugCut8 != NULL){fOutputContainer->Add(fV0MassDebugCut8);}
  
  return fOutputContainer;
*/
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
    nameMC.Form("ESD_EP_Mapping-Phi%02d",phi);
    TString titleESD="";
    titleMC.Form("Electron-Positron ESD Mapping-Phi%02d",phi);
    
    AddHistogram(nameESD, titleESD, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
  }


  for(Int_t r =0; r<=nRIndex;r++){
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
    nameMCPhiInR.Form("MC_EP_Mapping_Phi_vs_R_R-%02d",r);
    TString titleMCPhiInR="";
    titleMCPhiInR.Form("Electron-Positron MC Mapping of Phi in R%02d",r);
    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, nYBins, firstY, lastY, xAxisTitle, yAxisTitle);

    TString nameESDPhiInR="";
    nameESDPhiInR.Form("ESD_EP_Mapping_Phi_vs_R_R-%02d",r);
    TString titleESDPhiInR="";
    titleESDPhiInR.Form("Electron-Positron ESD Mapping of Phi in R%02d",r);
    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, nYBins, firstY, lastY, xAxisTitle, yAxisTitle);    
  }
}
