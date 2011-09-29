/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt                        *
 * Version 1.1                                                            *
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
#include "TH3F.h"
#include "AliLog.h"

using namespace std;

ClassImp(AliGammaConversionHistograms)


AliGammaConversionHistograms::AliGammaConversionHistograms() :
  fHistogramMap(new TMap()),
  fNPhiIndex(0),
  fNRIndex(0),
  fNZIndex(0),
//  fRBinLimits(0),
//  fZBinLimits(0),
  fMinRadius(0.),
  fMaxRadius(0.),
  fDeltaR(0.),
  fMinPhi(0.),
  fMaxPhi(0.),
  fDeltaPhi(0.),
  fMappingContainer(NULL),
  fBackgroundContainer(NULL),
  fDebugContainer(NULL),
  fResolutionContainer(NULL),
  fMatchContainer(NULL),
  fESDContainer(NULL),
  fMCContainer(NULL),
  fTableContainer(NULL),	
  fOtherContainer(NULL),
  f3DContainer(NULL)
{
  // see header file for documenation
  for(Int_t i=0;i<14;i++){
    fRBinLimits[i]=0.;
  }
  for(Int_t i=0;i<12;i++){
    fZBinLimits[i]=0.;
  }
}


AliGammaConversionHistograms::AliGammaConversionHistograms(const AliGammaConversionHistograms & original) :
  fHistogramMap(original.fHistogramMap),
  fNPhiIndex(original.fNPhiIndex),
  fNRIndex(original.fNRIndex),
  fNZIndex(original.fNZIndex),
  //  fRBinLimits(original.fRBinLimits),
  //  fZBinLimits(original.fZBinLimits),
  fMinRadius(original.fMinRadius),
  fMaxRadius(original.fMaxRadius),
  fDeltaR(original.fDeltaR),
  fMinPhi(original.fMinPhi),
  fMaxPhi(original.fMaxPhi),
  fDeltaPhi(original.fDeltaPhi),
  fMappingContainer(original.fMappingContainer),
  fBackgroundContainer(original.fBackgroundContainer),
  fDebugContainer(original.fDebugContainer),
  fResolutionContainer(original.fResolutionContainer),
  fMatchContainer(original.fMatchContainer),
  fESDContainer(original.fESDContainer),
  fMCContainer(original.fMCContainer),
  fTableContainer(original.fTableContainer), 
  fOtherContainer(original.fOtherContainer),
  f3DContainer(original.f3DContainer)
{    
  //see header file for documentation
  for(Int_t i=0;i<14;i++){
    fRBinLimits[i]= original.fRBinLimits[i];
  }
  for(Int_t i=0;i<12;i++){
    fZBinLimits[i]=original.fZBinLimits[i];
  }
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
  if( fHistogramMap->Contains(histogramName.Data()) ==  kFALSE ){
    TH1F *tmp = new TH1F(histogramName, histogramTitle,nXBins,firstX,lastX);
    tmp->GetXaxis()->SetTitle(xAxisTitle);
    tmp->GetYaxis()->SetTitle(yAxisTitle);
    TObjString* tobjstring = new TObjString(histogramName.Data());
    fHistogramMap->Add((TObject*)tobjstring,(TObject*)tmp);
  }
  else{
    cout << "Warning: Histogram ( "<<histogramName.Data()<<" ) already exists " << endl;
  }
}

void AliGammaConversionHistograms::AddHistogram(TString histogramName, TString histogramTitle, Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle, Int_t logAxis){
  // see header file for documentation
  if( fHistogramMap->Contains(histogramName.Data()) ==  kFALSE ){
    TH2F *tmp = new TH2F(histogramName, histogramTitle,nXBins,firstX,lastX,nYBins,firstY,lastY);
    tmp->GetXaxis()->SetTitle(xAxisTitle);
    tmp->GetYaxis()->SetTitle(yAxisTitle);
    TObjString *tobjstring = new TObjString(histogramName.Data());
    fHistogramMap->Add((TObject*)tobjstring,(TObject*)tmp);

    if(logAxis >= 0){
      BinLogAxis(histogramName.Data(), logAxis);
    }
  }
  else{
    cout << "Warning: Histogram ( "<<histogramName.Data()<<" ) already exists " << endl;
  }
}

void AliGammaConversionHistograms::AddHistogram(TString histogramName, TString histogramTitle, Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, Int_t nZBins, Double_t firstZ, Double_t lastZ, TString xAxisTitle, TString yAxisTitle, TString zAxisTitle, Int_t logAxis ){
  // see header file for documentation
  if( fHistogramMap->Contains(histogramName.Data()) ==  kFALSE ){
    TH3F *tmp = new TH3F(histogramName, histogramTitle,nXBins,firstX,lastX,nYBins,firstY,lastY,nZBins,firstZ,lastZ);
    tmp->GetXaxis()->SetTitle(xAxisTitle);
    tmp->GetYaxis()->SetTitle(yAxisTitle);
    tmp->GetZaxis()->SetTitle(zAxisTitle);
    TObjString *tobjstring = new TObjString(histogramName.Data());
    fHistogramMap->Add((TObject*)tobjstring,(TObject*)tmp);
    if(logAxis >= 0){
      BinLogAxis(histogramName.Data(), logAxis);
    }
  }
  else{
    cout << "Warning: Histogram ( "<<histogramName.Data()<<" ) already exists " << endl;
  }
}


Bool_t AliGammaConversionHistograms::BinLogAxis(const char* name, Int_t dim){

  //
  // converts the axis (defined by the dimension) of THx or THnSparse
  // object to Log scale. Number of bins and bin min and bin max are preserved
  
 
  TObject *o =  fHistogramMap->GetValue(name);
  TAxis *axis = 0x0;
  if(o->InheritsFrom("TH1")){
    axis = (dynamic_cast<TH1F*>(o))->GetXaxis();
  }
  if(o->InheritsFrom("TH2")){
    if(0 == dim){
      axis = (dynamic_cast<TH2F*>(o))->GetXaxis();
    }
    else if(1 == dim){
      axis = (dynamic_cast<TH2F*>(o))->GetYaxis();
    }
     else{
       //  AliError("Only dim = 0 or 1 possible for TH2F");
     }
  }
  //  if(o->InheritsFrom("THnSparse")){
  //  axis = (dynamic_cast<THnSparse*>(o))->GetAxis(dim);
  //}

  if(!axis){
    //AliError(Form("Axis '%d' could not be identified in the object '%s'\n", dim, name));
    return kFALSE;
  }

  Int_t bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  if(from <= 0){
    // AliError(Form(" Log binning not possible for object '%s'because the '%d' axis starts from '%f\n'", name, dim, from));
    return kFALSE;
  }
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins+1];
  newBins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i){
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;

  return kTRUE;


}

void AliGammaConversionHistograms::AddTable(TString tableName,TString tableTitle,Int_t nXBins,const char * axesLabel[]){
  //see header file for documentation

  if( fHistogramMap->Contains(tableName.Data()) ==  kFALSE ){
    TH1F *tmp = new TH1F(tableName,tableTitle,nXBins,0,nXBins);
    for(Int_t xbin=1; xbin<=nXBins; xbin++){
      tmp->GetXaxis()->SetBinLabel(xbin,axesLabel[xbin-1]);
    }
    tmp->SetStats(0);
    
    TObjString *tobjstring = new TObjString(tableName.Data());
    fHistogramMap->Add((TObject*)tobjstring,(TObject*)tmp);
  }
  else{
    cout << "Warning: Table ( "<<tableName.Data()<<" ) already exists " << endl;
  }
}
void AliGammaConversionHistograms::AddTable(TString tableName,TString tableTitle,Int_t nXBins,const char * axesXLabel[],Int_t nYBins,const char * axesYLabel[]){
  //see header file for documentation

  if( fHistogramMap->Contains(tableName.Data()) ==  kFALSE ){
    TH2F *tmp = new TH2F(tableName,tableTitle,nXBins,0,nXBins,nYBins,0,nYBins);
    for(Int_t xbin=1; xbin<=nXBins; xbin++){
      tmp->GetXaxis()->SetBinLabel(xbin,axesXLabel[xbin-1]);
    }
    for(Int_t ybin=1; ybin<=nYBins; ybin++){
      tmp->GetYaxis()->SetBinLabel(ybin,axesYLabel[ybin-1]);
    }
    tmp->SetStats(0);
    
    TObjString *tobjstring = new TObjString(tableName.Data());
    fHistogramMap->Add((TObject*)tobjstring,(TObject*)tmp);
  }
  else{
    cout << "Warning: Table ( "<<tableName.Data()<<" ) already exists " << endl;
  }
}

void AliGammaConversionHistograms::AddTable(TString tableName,TString tableTitle,Int_t nXBins,const char * axesXLabel[],Int_t nYBins,const char * axesYLabel[], Int_t nZBins,const char * axesZLabel[]){
  //see header file for documentation

  if( fHistogramMap->Contains(tableName.Data()) ==  kFALSE ){
    TH3F *tmp = new TH3F(tableName,tableTitle,nXBins,0,nXBins,nYBins,0,nYBins,nZBins,0,nZBins);
    for(Int_t xbin=1; xbin<=nXBins; xbin++){
      tmp->GetXaxis()->SetBinLabel(xbin,axesXLabel[xbin-1]);
    }
    for(Int_t ybin=1; ybin<=nYBins; ybin++){
      tmp->GetYaxis()->SetBinLabel(ybin,axesYLabel[ybin-1]);
    }
    for(Int_t zbin=1; zbin<=nZBins; zbin++){
      tmp->GetZaxis()->SetBinLabel(zbin,axesZLabel[zbin-1]);
    }
    
    tmp->SetStats(0);
    
    TObjString *tobjstring = new TObjString(tableName.Data());
    fHistogramMap->Add((TObject*)tobjstring,(TObject*)tmp);
  }
  else{
    cout << "Warning: Table ( "<<tableName.Data()<<" ) already exists " << endl;
  }
}

void AliGammaConversionHistograms::FillTable(TString tableName,Double_t xValue) const {
  //see header file for documentation
  TH1 *tmp = (TH1*)fHistogramMap->GetValue(tableName.Data());
  if(tmp){
    tmp->Fill(xValue);
  }
}
void AliGammaConversionHistograms::FillTable(TString tableName,Double_t xValue,Double_t yValue) const {
  //see header file for documentation
  TH2 *tmp = (TH2*)fHistogramMap->GetValue(tableName.Data());
  if(tmp){
    tmp->Fill(xValue,yValue);
  }
}
void AliGammaConversionHistograms::FillTable(TString tableName,Double_t xValue,Double_t yValue, Double_t zValue) const {
  //see header file for documentation
  TH3 *tmp = (TH3*)fHistogramMap->GetValue(tableName.Data());
  if(tmp){
    tmp->Fill(xValue,yValue,zValue);
  }
}

void AliGammaConversionHistograms::FillHistogram(TString histogramName, Double_t xValue) const{
  //see header file for documentation
  TH1 *tmp = (TH1*)fHistogramMap->GetValue(histogramName.Data());
  if(tmp){
    tmp->Fill(xValue);
  }
}

void AliGammaConversionHistograms::FillHistogram(TString histogramName, Double_t xValue, Double_t yValue) const{
  //see header file for documentation
  TH1 *tmp = (TH1*)fHistogramMap->GetValue(histogramName.Data());
  if(tmp){
    tmp->Fill(xValue, yValue);
  }
}

void AliGammaConversionHistograms::FillHistogram(TString histogramName, Double_t xValue, Double_t yValue, Double_t zValue) const{
  //see header file for documentation
  TH3 *tmp = (TH3*)fHistogramMap->GetValue(histogramName.Data());
  if(tmp){
    tmp->Fill(xValue, yValue, zValue);
  }
}


TObject* AliGammaConversionHistograms::GetValue(const TString& name)
{ 
  //Get pointer to histogram with name
  return fHistogramMap->GetValue(name.Data());
}

void AliGammaConversionHistograms::GetOutputContainer(TList *fOutputContainer){
  //checking if the container is alrerady created
	
  if(fOutputContainer == NULL){
    cout<<"WARNING: GetOutputContainer: output container object is NULL"<<endl;
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
	  fMappingContainer->SetOwner(kTRUE);
	  fMappingContainer->SetName("Mapping histograms");
	}
	fMappingContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      }
      else if(histogramString.Contains("Background")){// means it should be put in the background folder
	if(fBackgroundContainer == NULL){
	  fBackgroundContainer = new TList();
	  fBackgroundContainer->SetOwner(kTRUE);
	  fBackgroundContainer->SetName("Background histograms");
	}
	fBackgroundContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      }
      else if(histogramString.Contains("Debug")){// means it should be put in the debug folder
	if(fDebugContainer == NULL){
	  fDebugContainer = new TList();
	  fDebugContainer->SetOwner(kTRUE);
	  fDebugContainer->SetName("Debug histograms");
	}
	fDebugContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      }
      else if(histogramString.Contains("Resolution")){// means it should be put in the resolution folder
	if(fResolutionContainer == NULL){
	  fResolutionContainer = new TList();
	  fResolutionContainer->SetOwner(kTRUE);
	  fResolutionContainer->SetName("Resolution histograms");
	}
	fResolutionContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      }
      else if(histogramString.Contains("TrueConv")){// means it should be put in the true conv folder
	if(fMatchContainer == NULL){
	  fMatchContainer = new TList();
	  fMatchContainer->SetOwner(kTRUE);
	  fMatchContainer->SetName("True conversion histograms");
	}
	fMatchContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      }
      else if(histogramString.Contains("ESD")){// means it should be put in the ESD folder
	if(fESDContainer == NULL){
	  fESDContainer = new TList();
	  fESDContainer->SetOwner(kTRUE);
	  fESDContainer->SetName("ESD histograms");
	}
	fESDContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      }
      else if(histogramString.Contains("MC")){// means it should be put in the MC folder
	if(fMCContainer == NULL){
	  fMCContainer = new TList();
	  fMCContainer->SetOwner(kTRUE);
	  fMCContainer->SetName("MC histograms");
	}
	fMCContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      }
      else if(histogramString.Contains("Table")){// means it should be put in the Table Folder
	if(fTableContainer == NULL){
	   fTableContainer = new TList();
	   fTableContainer->SetOwner(kTRUE);
	   fTableContainer->SetName("Tables");
	}
	fTableContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      }
	else if(histogramString.Contains("3DPlots")){// means it should be put in the Table Folder
		if(f3DContainer == NULL){
			f3DContainer = new TList();
			f3DContainer->SetOwner(kTRUE);
			f3DContainer->SetName("3D histograms");
		}
		f3DContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      }
      else{
	if(fOtherContainer == NULL){
	  fOtherContainer = new TList();
	  fOtherContainer->SetOwner(kTRUE);
	  fOtherContainer->SetName("Other histograms");
	}
	fOtherContainer->Add((TH1*)fHistogramMap->GetValue(histogramString.Data()));
      }
      histogramName = NULL;
    } // end while

    if(fMappingContainer != NULL){
      fOutputContainer->Add(fMappingContainer);
    }
    if(fBackgroundContainer != NULL){
      fOutputContainer->Add(fBackgroundContainer);
    }
    if(fDebugContainer != NULL){
      fOutputContainer->Add(fDebugContainer);
    }
    if(fResolutionContainer != NULL){
      fOutputContainer->Add(fResolutionContainer);
    }
    if(fMatchContainer != NULL){
      fOutputContainer->Add(fMatchContainer);
    }
    if(fESDContainer != NULL){
      fOutputContainer->Add(fESDContainer);
    }
    if(fMCContainer != NULL){
      fOutputContainer->Add(fMCContainer);
    }
    if(fTableContainer !=  NULL){
       fOutputContainer->Add(fTableContainer);	
    }
	if(f3DContainer !=  NULL){
       fOutputContainer->Add(f3DContainer);	
    }
    if(fOtherContainer != NULL){
      fOutputContainer->Add(fOtherContainer);
    }
  }
}

Int_t AliGammaConversionHistograms::GetRBin(Double_t radius) const{
  // see header file for documentation
  Int_t iResult=0;
//   if(fDeltaR>0){
//     iResult = (Int_t)((radius - fMinRadius)/fDeltaR);
//   }
  for(Int_t i=0;i<fNRIndex;i++){
    //    cout<<"Test-limits::"<< fRBinLimits[i]<<endl;
    if( radius>=fRBinLimits[i] && radius<fRBinLimits[i+1] ){
      iResult=i;
    }
  }
  return iResult;
}

Int_t AliGammaConversionHistograms::GetZBin(Double_t zPos) const{
  // see header file for documentation
  Int_t iResult=0;

  for(Int_t i=0;i<fNZIndex;i++){
    //    cout<<"Test-limits::"<< fZBinLimits[i]<<endl;
    if( zPos>=fZBinLimits[i] && zPos<fZBinLimits[i+1] ){
      iResult=i;
    }
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
  if(nRIndex<=14){
    fNRIndex   = nRIndex;
  }else{
    fNRIndex=14;
  }

  fNZIndex = 13;

  //  fRBinLimits= new Double_t[8];   Kenneth: moved from pointer to fixed array
  /*
  fRBinLimits[0]=0.;
  fRBinLimits[1]=13.;   //changed from 12 to 13: A. Marin 01.03.10
  fRBinLimits[2]=21.;   //changed from 22 to 21: A. Marin 01.03.10 
  fRBinLimits[3]=35.;
  fRBinLimits[4]=55.;
  fRBinLimits[5]=72.;
  fRBinLimits[6]=90.;
  fRBinLimits[7]=500.;
  */

  fRBinLimits[0]=0.;
  fRBinLimits[1]=3.5;
  fRBinLimits[2]=5.75;
  fRBinLimits[3]=9.5;
  fRBinLimits[4]=13.;
  fRBinLimits[5]=21.;
  fRBinLimits[6]=27.5;
  fRBinLimits[7]=35.;
  fRBinLimits[8]=42.;
  fRBinLimits[9]=55.;
  fRBinLimits[10]=72.;
  fRBinLimits[11]=79.5; // change from 81.5 to 79.5 to have CE in 1 r bin 81.05
  fRBinLimits[12]=90.;
  fRBinLimits[13]=500.;



  //  fZBinLimits= new Double_t[7]; Kenneth: moved from pointer to fixed array
  fZBinLimits[0]=-500.;
  fZBinLimits[1]=-200.;
  fZBinLimits[2]=-100.;
  fZBinLimits[3]=-50.;
  fZBinLimits[4]=-30.;
  fZBinLimits[5]=-15.;
  fZBinLimits[6]= 0.;
  fZBinLimits[7]= 15.;
  fZBinLimits[8]= 30.;
  fZBinLimits[9]= 50.;
  fZBinLimits[10]=100.;
  fZBinLimits[11]=200.;
  fZBinLimits[12]=500.;


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

  Double_t tmptogetridofwarning = firstX + lastX + nYBins + firstY + lastY;
  if(tmptogetridofwarning < 0){
    cout<<"Less than zero"<<endl;
  }
	
  for(Int_t phi =0; phi<fNPhiIndex;phi++){
		
    for(Int_t r =0; r<fNRIndex;r++){
			
      // setting axis to "" changes below
      xAxisTitle="z [cm]";
      yAxisTitle="#eta";
	  
      //Creating the axis titles
      //if(xAxisTitle.Length() == 0){
	//xAxisTitle.Form("Phi %02d",phi);
	 //      }
			
      //if(yAxisTitle.Length() == 0){
	//yAxisTitle.Form("R %02d",phi);
		// }
			
      //MC
      TString nameMC="";
      nameMC.Form("MC_Conversion_Mapping_Phi%02d_R%02d",phi,r);
      TString titleMC="";
      titleMC.Form("Electron-Positron MC Mapping-Phi%02d-R%02d",phi,r);
			
      //AddHistogram(nameMC, titleMC, nXBins, firstX, lastX, nYBins, firstY, lastY, xAxisTitle, yAxisTitle);
			
      //ESD
      TString nameESD="";
      nameESD.Form("ESD_Conversion_Mapping_Phi%02d_R%02d",phi,r);
      TString titleESD="";
      titleESD.Form("Electron-Positron ESD Mapping-Phi%02d-R%02d",phi,r);
			
      //AddHistogram(nameESD, titleESD, nXBins, firstX, lastX, nYBins, firstY, lastY, xAxisTitle, yAxisTitle);
    }
  }
	
	
  for(Int_t phi =0; phi<=nPhiIndex;phi++){ 
		
    // setting axis to "" changes below
     xAxisTitle="z [cm]";
      yAxisTitle="#eta";
   //Creating the axis titles
    //if(xAxisTitle.Length() == 0){
    //  xAxisTitle.Form("Phi %02d",phi);
    //}
    //if(yAxisTitle.Length() == 0){
    //  yAxisTitle = "Counts";
    //}
		
    //MC
    TString nameMC="";
    nameMC.Form("MC_Conversion_Mapping_Phi%02d",phi);
    TString titleMC="";
    titleMC.Form("Electron-Positron MC Mapping-Phi%02d",phi);
		
    //AddHistogram(nameMC, titleMC, nXBins, firstX, lastX, nYBins, firstY, lastY, xAxisTitle, yAxisTitle);
		
    //MC
    TString nameESD="";
    nameESD.Form("ESD_Conversion_Mapping_Phi%02d",phi);
    TString titleESD="";
    titleESD.Form("Electron-Positron ESD Mapping-Phi%02d",phi);
		
    // AddHistogram(nameESD, titleESD, nXBins, firstX, lastX,nYBins, firstY, lastY, xAxisTitle, yAxisTitle);
  }
	
	
  for(Int_t r =0; r<nRIndex;r++){
		
    // setting axis to "" changes below
    xAxisTitle="#phi";
    yAxisTitle="counts";
    //Creating the axis titles
    //if(xAxisTitle.Length() == 0){
    //  xAxisTitle.Form("R %02d",r);
    //}
    //if(yAxisTitle.Length() == 0){
    //  yAxisTitle = "Counts";
    //}
		
    //MC
    TString nameMC="";
    nameMC.Form("MC_Conversion_Mapping_R%02d",r);
    TString titleMC="";
    titleMC.Form("Electron-Positron MC Mapping-R%02d",r);
		
    // AddHistogram(nameMC, titleMC, nXBins, firstX, lastX, nYBins, firstY, lastY, xAxisTitle, yAxisTitle);
		
    //ESD
    TString nameESD="";
    nameESD.Form("ESD_Conversion_Mapping_R%02d",r);
    TString titleESD="";
    titleESD.Form("Electron-Positron ESD Mapping-R%02d",r);
		
    //AddHistogram(nameESD, titleESD, nXBins, firstX, lastX,nYBins, firstY, lastY, xAxisTitle, yAxisTitle);
		
    //Mapping Phi in R
    TString nameMCPhiInR="";
    nameMCPhiInR.Form("MC_Conversion_Mapping_Phi_in_R_%02d",r);
    TString titleMCPhiInR="";
    titleMCPhiInR.Form("MC Mapping of Phi in R%02d",r);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);
		

    //Mapping Z in R
    TString nameMCZInR="";
    nameMCZInR.Form("MC_Conversion_Mapping_Z_in_R_%02d",r);
    TString titleMCZInR="";
    titleMCZInR.Form("MC Mapping of Z in R%02d",r);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCZInR, titleMCZInR, nXBins, -300, 300, xAxisTitle, yAxisTitle);


   //Mapping Phi in R Middle Pt
    TString nameMCMidPtPhiInR="";
    nameMCMidPtPhiInR.Form("MC_Conversion_Mapping_MidPt_Phi_in_R_%02d",r);
    TString titleMCMidPtPhiInR="";
    titleMCMidPtPhiInR.Form("MC Mapping Middle Pt of Phi in R%02d",r);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCMidPtPhiInR, titleMCMidPtPhiInR, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);
		

    //Mapping Z in R Middle Pt
    TString nameMCMidPtZInR="";
    nameMCMidPtZInR.Form("MC_Conversion_Mapping_MidPt_Z_in_R_%02d",r);
    TString titleMCMidPtZInR="";
    titleMCMidPtZInR.Form("MC Mapping Middle Pt of Z in R%02d",r);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCMidPtZInR, titleMCMidPtZInR, nXBins, -300, 300, xAxisTitle, yAxisTitle);




    //Mapping Phi in R
    TString nameESDPhiInR="";
    nameESDPhiInR.Form("ESD_Conversion_Mapping_Phi_in_R_%02d",r);
    TString titleESDPhiInR="";
    titleESDPhiInR.Form("ESD Mapping of Phi in R%02d",r);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);    

   //Mapping Z in R
    TString nameESDZInR="";
    nameESDZInR.Form("ESD_Conversion_Mapping_Z_in_R_%02d",r);
    TString titleESDZInR="";
    titleESDZInR.Form("ESD Mapping of Z in R%02d",r);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDZInR, titleESDZInR, nXBins, -300, 300, xAxisTitle, yAxisTitle);    

    //Mapping Phi in R Middle Pt 
    TString nameESDMidPtPhiInR="";
    nameESDMidPtPhiInR.Form("ESD_Conversion_Mapping_MidPt_Phi_in_R_%02d",r);
    TString titleESDMidPtPhiInR="";
    titleESDMidPtPhiInR.Form("ESD Mapping Middle Pt of Phi in R%02d",r);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDMidPtPhiInR, titleESDMidPtPhiInR, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);    

   //Mapping Z in R Middle Pt
    TString nameESDMidPtZInR="";
    nameESDMidPtZInR.Form("ESD_Conversion_Mapping_MidPt_Z_in_R_%02d",r);
    TString titleESDMidPtZInR="";
    titleESDMidPtZInR.Form("ESD Mapping Middle Pt of Z in R%02d",r);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDMidPtZInR, titleESDMidPtZInR, nXBins, -300, 300, xAxisTitle, yAxisTitle);    


 
  }



  for(Int_t z =0; z<fNZIndex;z++){
    //Mapping Phi in Z
    TString nameMCPhiInZ="";
    nameMCPhiInZ.Form("MC_Conversion_Mapping_Phi_in_Z_%02d",z);
    TString titleMCPhiInZ="";
    titleMCPhiInZ.Form("MC Mapping of Phi in Z%02d",z);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCPhiInZ, titleMCPhiInZ, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);
 
    //Mapping Phi in Z for FMD
    TString nameMCFMDPhiInZ="";
    nameMCFMDPhiInZ.Form("MC_Conversion_Mapping_FMD_Phi_in_Z_%02d",z);
    TString titleMCFMDPhiInZ="";
    titleMCFMDPhiInZ.Form("MC Mapping FMD of Phi in Z%02d",z);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCFMDPhiInZ, titleMCFMDPhiInZ, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);
		
    //Mapping Phi in Z for ITSTPC
    TString nameMCITSTPCPhiInZ="";
    nameMCITSTPCPhiInZ.Form("MC_Conversion_Mapping_ITSTPC_Phi_in_Z_%02d",z);
    TString titleMCITSTPCPhiInZ="";
    titleMCITSTPCPhiInZ.Form("MC Mapping ITSTPC of Phi in Z%02d",z);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCITSTPCPhiInZ, titleMCITSTPCPhiInZ, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);


    //Mapping R in Z
    TString nameMCRInZ="";
    nameMCRInZ.Form("MC_Conversion_Mapping_R_in_Z_%02d",z);
    TString titleMCRInZ="";
    titleMCRInZ.Form("MC Mapping of R in Z%02d",z);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCRInZ, titleMCRInZ, nXBins, fMinRadius, fMaxRadius, xAxisTitle, yAxisTitle);

   //Mapping Phi in Z Middle Pt
    TString nameMCMidPtPhiInZ="";
    nameMCMidPtPhiInZ.Form("MC_Conversion_Mapping_MidPt_Phi_in_Z_%02d",z);
    TString titleMCMidPtPhiInZ="";
    titleMCMidPtPhiInZ.Form("MC Mapping Middle Pt of Phi in Z%02d",z);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCMidPtPhiInZ, titleMCMidPtPhiInZ, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);
		
   //Mapping Phi in Z Middle Pt for FMD
    TString nameMCMidPtFMDPhiInZ="";
    nameMCMidPtFMDPhiInZ.Form("MC_Conversion_Mapping_MidPt_FMD_Phi_in_Z_%02d",z);
    TString titleMCMidPtFMDPhiInZ="";
    titleMCMidPtFMDPhiInZ.Form("MC Mapping Middle Pt of Phi in Z%02d",z);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCMidPtFMDPhiInZ, titleMCMidPtFMDPhiInZ, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);
		


    //Mapping R in Z Middle Pt
    TString nameMCMidPtRInZ="";
    nameMCMidPtRInZ.Form("MC_Conversion_Mapping_MidPt_R_in_Z_%02d",z);
    TString titleMCMidPtRInZ="";
    titleMCMidPtRInZ.Form("MC Mapping Middle Pt of R in Z%02d",z);
    //    AddHistogram(nameMCPhiInR, titleMCPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);
    AddHistogram(nameMCMidPtRInZ, titleMCMidPtRInZ, nXBins, fMinRadius, fMaxRadius, xAxisTitle, yAxisTitle);




    //Mapping Phi in Z
    TString nameESDPhiInZ="";
    nameESDPhiInZ.Form("ESD_Conversion_Mapping_Phi_in_Z_%02d",z);
    TString titleESDPhiInZ="";
    titleESDPhiInZ.Form("ESD Mapping of Phi in Z%02d",z);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDPhiInZ, titleESDPhiInZ, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);    


    //Mapping Phi in Z for FMD
    TString nameESDFMDPhiInZ="";
    nameESDFMDPhiInZ.Form("ESD_Conversion_Mapping_FMD_Phi_in_Z_%02d",z);
    TString titleESDFMDPhiInZ="";
    titleESDFMDPhiInZ.Form("ESD Mapping FMD of Phi in Z%02d",z);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDFMDPhiInZ, titleESDFMDPhiInZ, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);    

    //Mapping Phi in Z for ITSTPC
    TString nameESDITSTPCPhiInZ="";
    nameESDITSTPCPhiInZ.Form("ESD_Conversion_Mapping_ITSTPC_Phi_in_Z_%02d",z);
    TString titleESDITSTPCPhiInZ="";
    titleESDITSTPCPhiInZ.Form("ESD Mapping ITSTPC of Phi in Z%02d",z);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDITSTPCPhiInZ, titleESDITSTPCPhiInZ, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);    


   //Mapping R in Z
    TString nameESDRInZ="";
    nameESDRInZ.Form("ESD_Conversion_Mapping_R_in_Z_%02d",z);
    TString titleESDRInZ="";
    titleESDRInZ.Form("ESD Mapping of R in Z%02d",z);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDRInZ, titleESDRInZ, nXBins, fMinRadius, fMaxRadius, xAxisTitle, yAxisTitle);    


   //Mapping Phi in Z Middle Pt
    TString nameESDMidPtPhiInZ="";
    nameESDMidPtPhiInZ.Form("ESD_Conversion_Mapping_MidPt_Phi_in_Z_%02d",z);
    TString titleESDMidPtPhiInZ="";
    titleESDMidPtPhiInZ.Form("ESD Mapping Middle Ptof Phi in R%02d",z);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDMidPtPhiInZ, titleESDMidPtPhiInZ, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);    

   //Mapping Phi in Z Middle Pt for FMD
    TString nameESDMidPtFMDPhiInZ="";
    nameESDMidPtFMDPhiInZ.Form("ESD_Conversion_Mapping_MidPt_FMD_Phi_in_Z_%02d",z);
    TString titleESDMidPtFMDPhiInZ="";
    titleESDMidPtFMDPhiInZ.Form("ESD Mapping Middle Pt FMD of Phi in Z%02d",z);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDMidPtFMDPhiInZ, titleESDMidPtFMDPhiInZ, nXBins, -TMath::Pi(), TMath::Pi(), xAxisTitle, yAxisTitle);    


   //Mapping R in Z Middle Pt
    TString nameESDMidPtRInZ="";
    nameESDMidPtRInZ.Form("ESD_Conversion_Mapping_MidPt_R_in_Z_%02d",z);
    TString titleESDMidPtRInZ="";
    titleESDMidPtRInZ.Form("ESD Mapping Middle Pt of R in Z%02d",z);
    //    AddHistogram(nameESDPhiInR, titleESDPhiInR, nXBins, firstX, lastX, xAxisTitle, yAxisTitle);    
    AddHistogram(nameESDMidPtRInZ, titleESDMidPtRInZ, nXBins, fMinRadius, fMaxRadius, xAxisTitle, yAxisTitle);    



  }



}
