/**************************************************************************
 * Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
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

// Mean Vertex preprocessor:
// 2) takes data after  pass0 ,
// processes it, and stores either to OCDB .
//
// Davide Caffarri

#include "AliMeanVertexPreprocessorOffline.h"

#include "AliCDBStorage.h"
#include "AliCDBMetaData.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliLog.h"

#include <TTimeStamp.h>
#include <TFile.h>
#include <TObjString.h>
#include <TNamed.h>
#include "TClass.h"
#include <TCanvas.h>
#include "TString.h"

#include "AliESDVertex.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

ClassImp(AliMeanVertexPreprocessorOffline)

//_______________________________________________________

const Char_t *AliMeanVertexPreprocessorOffline::fgkStatusCodeName[AliMeanVertexPreprocessorOffline::kNStatusCodes] = {
  "ok",
  "open file error or missing histos",
  "too low statistics",
  "problems storing OCDB",
  "write MeanVertex computed online",
  "write SPD vtx offline",
  "lumi region or cov matrix computation problems, default values set"
  "problme in the fit of Z coordinate for Cpass1 update"
};


//____________________________________________________
AliMeanVertexPreprocessorOffline::AliMeanVertexPreprocessorOffline():
  TNamed("AliMeanVertexPreprocessorOffline","AliMeanVertexPreprocessorOffline"),
  fStatus(kOk),
  fShowPlots(kFALSE)
{
  //constructor
}
//____________________________________________________

AliMeanVertexPreprocessorOffline::~AliMeanVertexPreprocessorOffline()
{
  //destructor
    
}

//____________________________________________________
void AliMeanVertexPreprocessorOffline::ProcessOutput(const char *filename, AliCDBStorage *db, Int_t runNb) {
    
  TFile *file = TFile::Open(filename);
  if (!file || !file->IsOpen()) {
    AliError(Form("cannot open outputfile %s", filename));
    fStatus=kInputError;
    return;
  }
    
  if (!db) {
    AliError("no OCDB storage found, return");
    fStatus=kInputError;
    return;
  }
    
  AliCDBManager *manCheck = AliCDBManager::Instance();
  manCheck->SetDefaultStorage("raw://");
  manCheck->SetRun(runNb);
    
  TList *list = (TList*)file->Get("MeanVertex");
    
  TH1F *histTRKvtxX = 0x0;
  TH1F *histTRKvtxY = 0x0;
  TH1F *histTRKvtxZ = 0x0;
    
  TH1F *histITSSAvtxX = 0x0;
  TH1F *histITSSAvtxY = 0x0;
  TH1F *histITSSAvtxZ = 0x0;
    
  TH1F *histSPDvtxX = 0x0;
  TH1F *histSPDvtxY = 0x0;
  TH1F *histSPDvtxZ = 0x0;
    
    
  if (!list) {
        
    histTRKvtxX = (TH1F*)file->Get("hTRKVertexX");
    histTRKvtxY = (TH1F*)file->Get("hTRKVertexY");
    histTRKvtxZ = (TH1F*)file->Get("hTRKVertexZ");
        
    histITSSAvtxX = (TH1F*)file->Get("hITSSAVertexX");
    histITSSAvtxY = (TH1F*)file->Get("hITSSAVertexY");
    histITSSAvtxZ = (TH1F*)file->Get("hITSSAVertexZ");
        
    histSPDvtxX = (TH1F*)file->Get("hSPDVertexX");
    histSPDvtxY = (TH1F*)file->Get("hSPDVertexY");
    histSPDvtxZ = (TH1F*)file->Get("hSPDVertexZ");
        
  }else{
        
    histTRKvtxX = (TH1F*)list->FindObject("hTRKVertexX");
    histTRKvtxY = (TH1F*)list->FindObject("hTRKVertexY");
    histTRKvtxZ = (TH1F*)list->FindObject("hTRKVertexZ");
        
    histITSSAvtxX = (TH1F*)list->FindObject("hITSSAVertexX");
    histITSSAvtxY = (TH1F*)list->FindObject("hITSSAVertexY");
    histITSSAvtxZ = (TH1F*)list->FindObject("hITSSAVertexZ");
        
    histSPDvtxX = (TH1F*)list->FindObject("hSPDVertexX");
    histSPDvtxY = (TH1F*)list->FindObject("hSPDVertexY");
    histSPDvtxZ = (TH1F*)list->FindObject("hSPDVertexZ");
        
  }
    
  // Hierachical search for available vertices
  // globbal tracks -> ITS standalone tracks -> SPD tracklets
  Bool_t useTRKvtx = kTRUE;
  Bool_t useITSSAvtx = kFALSE;
  Bool_t useSPDvtx = kFALSE;
  Bool_t spdAvailable = kTRUE;
  Bool_t writeMeanVertexSPD = kFALSE;
  Bool_t vertexerSPD3Doff=kFALSE;
    
  if ( !histTRKvtxX || !histTRKvtxY || !histTRKvtxZ ||
       histTRKvtxX->GetEffectiveEntries()<50 ||
       histTRKvtxY->GetEffectiveEntries()<50 ||
       histTRKvtxZ->GetEffectiveEntries()<50)
    {
      AliError("TRK vertex histograms have too few entries for fitting");
      useTRKvtx=kFALSE;
      useITSSAvtx=kTRUE;
      histTRKvtxX = histITSSAvtxX;
      histTRKvtxY = histITSSAvtxY;
      histTRKvtxZ = histITSSAvtxZ;
    }
  if (useITSSAvtx &&
      (!histTRKvtxX || !histTRKvtxY || !histTRKvtxZ ||
       histTRKvtxX->GetEffectiveEntries()<50 ||
       histTRKvtxY->GetEffectiveEntries()<50 ||
       histTRKvtxZ->GetEffectiveEntries()<50)
      ) {
    AliError("ITSSA vertex histograms have too few entries for fitting");
    useITSSAvtx = kFALSE;
    useSPDvtx=kTRUE;
  }
    
  if(useSPDvtx){
    if (!histSPDvtxX || !histSPDvtxY || !histSPDvtxZ) {
      AliError("cannot find any histograms available from file");
      spdAvailable = kFALSE;
      fStatus=kInputError;
      return;
    }
        
    Float_t nEntriesX = histSPDvtxX->GetEffectiveEntries();
    Float_t nEntriesY = histSPDvtxY->GetEffectiveEntries();
    Float_t nEntriesZ = histSPDvtxZ->GetEffectiveEntries();
        
    if (nEntriesX < 50. || nEntriesY<50. || nEntriesZ<50.) {
      spdAvailable = kFALSE;
      if ((useTRKvtx==kFALSE) && (useITSSAvtx==kFALSE)) {
        AliError("Also SPD vertex histograms have too few entries for fitting, return");
        fStatus=kLowStatistics;
        return;
      }
    }
    if((nEntriesX == 0.)&&(nEntriesY==0.) && (nEntriesZ>0.)) {
      vertexerSPD3Doff = kTRUE;
      AliWarning("Vertexer SPD 3D off");
    }
  }
    
    
  const char* whichpass = gSystem->Getenv("CPASSMODE");
    
  int cPassMode = -1;
  if (whichpass && *whichpass == '0') {
    cPassMode = 0;
  }
  else if (whichpass && *whichpass == '1') {
    cPassMode = 1;
  }
  else {
    AliFatal("Environment variable CPASSMODE must be set to either 0 or 1");
  }
    
  if (cPassMode==1) {
    // check if there was a cpass0 already done, i.e. there is a run-specific version of MeanVertex
    AliCDBEntry* entryMV = manCheck->Get("GRP/Calib/MeanVertex");
    AliCDBId& idOldMV = entryMV->GetId();
    AliESDVertex* oldMVtx = (AliESDVertex*)entryMV->GetObject();
        
    if (idOldMV.GetFirstRun()==0 || idOldMV.GetLastRun()==999999999 ||
        oldMVtx->GetXRes()>2.8 // > pipe radius --> it's a dummy object, don't use it
        ) {
      AliWarningF("In CPass1 old MeanVertex has run range %d:%d and XRes:%.2f "
                  "-> switch to CPass0 mode with full update of meanVertex",
                  idOldMV.GetFirstRun(),idOldMV.GetLastRun(),oldMVtx->GetXRes());
      cPassMode = 0;
    }
  }
    
  if (cPassMode==0) {
        
    Double_t xMeanVtx=0., yMeanVtx=0., zMeanVtx=0.;
    Double_t xSigmaVtx=0., ySigmaVtx=0., zSigmaVtx=0.;
        
    TF1 *fitVtxX= new TF1("fitVtxX", "gaus");
    TF1 *fitVtxY= new TF1("fitVtxY", "gaus");
    TF1 *fitVtxZ= new TF1("fitVtxZ", "gaus");
    Double_t meanh,rmsh,peakh;
    Double_t minfit,maxfit;
    if (useTRKvtx || useITSSAvtx) {
      
      printf("\n------ x-vertex fit ------\n");
      meanh=histTRKvtxX->GetMean();
      rmsh=histTRKvtxX->GetRMS();
      peakh=histTRKvtxX->GetBinCenter(histTRKvtxX->GetMaximumBin());
      printf("mean=%f  peak=%f rms=%f\n",meanh,peakh,rmsh);
      if(TMath::Abs(meanh-peakh)<0.01){
        minfit=meanh-5.*rmsh;
        maxfit=meanh+5.*rmsh;
      }else{
        minfit=-0.3;
        maxfit=0.3;
      }
      printf("fit range: %f %f\n",minfit,maxfit);
      if(fShowPlots){
        TCanvas* cxv=new TCanvas("cxv","xvertex");
        cxv->cd();
        histTRKvtxX->Draw();
      }
      fitVtxX->SetParameters(histTRKvtxX->Integral(), meanh, rmsh);
      histTRKvtxX ->Fit(fitVtxX, "M", "", minfit,maxfit);
      xMeanVtx = fitVtxX -> GetParameter(1);
      if (TMath::Abs(xMeanVtx) > 2.) {
        xMeanVtx = 0.;
        writeMeanVertexSPD=kTRUE;
        fStatus=kWriteMeanVertexSPD;
      }
          
      printf("\n------ y-vertex fit ------\n");
      meanh=histTRKvtxY->GetMean();
      rmsh=histTRKvtxY->GetRMS();
      peakh=histTRKvtxY->GetBinCenter(histTRKvtxY->GetMaximumBin());
      printf("mean=%f  peak=%f rms=%f\n",meanh,peakh,rmsh);
      if(TMath::Abs(meanh-peakh)<0.01){
        minfit=meanh-5.*rmsh;
        maxfit=meanh+5.*rmsh;
      }else{
        minfit=-0.4;
        maxfit=0.7;
      }
      printf("fit range: %f %f\n",minfit,maxfit);
      if(fShowPlots){
        TCanvas* cyv=new TCanvas("cyv","yvertex");
        cyv->cd();
        histTRKvtxY->Draw();
      }
      fitVtxY->SetParameters(histTRKvtxY->Integral(), meanh, rmsh);
      histTRKvtxY ->Fit(fitVtxY, "M", "",minfit,maxfit);
      yMeanVtx = fitVtxY -> GetParameter(1);
      if (TMath::Abs(yMeanVtx) > 2.) {
        yMeanVtx = 0.;
        writeMeanVertexSPD=kTRUE;
        fStatus=kWriteMeanVertexSPD;
      }

      printf("\n------ z-vertex fit ------\n");
      if(fShowPlots){
        TCanvas* czv=new TCanvas("czv","zvertex");
        czv->cd();
        histTRKvtxZ->Draw();
      }
      histTRKvtxZ ->Fit(fitVtxZ, "M", "", -12, 12);
      zMeanVtx = fitVtxZ -> GetParameter(1);
      zSigmaVtx = fitVtxZ -> GetParameter(2);
      if ((TMath::Abs(zMeanVtx) > 20.) || (zSigmaVtx>12.)) {
        zMeanVtx = histTRKvtxZ->GetMean();
        zSigmaVtx = histTRKvtxZ->GetRMS();
        writeMeanVertexSPD=kTRUE;
        fStatus=kWriteMeanVertexSPD;
      }
        
      // Check fits: compare histo mean with fit mean value
      Double_t xHistoMean = histTRKvtxX -> GetMean();
      Double_t xHistoRMS = histTRKvtxX ->GetRMS();
      Double_t yHistoMean = histTRKvtxY ->GetMean();
      Double_t yHistoRMS = histTRKvtxY ->GetRMS();
      Double_t zHistoMean = histTRKvtxZ -> GetMean();
      Double_t zHistoRMS = histTRKvtxZ ->GetRMS();

      if (TMath::Abs(xHistoMean-xMeanVtx) > 0.5 ||
          TMath::Abs(yHistoMean-yMeanVtx) > 0.5 ||
          TMath::Abs(zHistoMean-zMeanVtx) > 1.) {
        AliWarning("Possible problems with the fit mean very different from histo mean... using SPD vertex");
        useTRKvtx = kFALSE;
        useITSSAvtx = kFALSE;
        useSPDvtx = kTRUE;
        fStatus=kUseOfflineSPDvtx;
      }
      AliDebug(2, Form("xHistoRMS = %f, yHistoRMS = %f, zHistoRMS = %f", xHistoRMS, yHistoRMS, zHistoRMS));
    }
       
    if ((useSPDvtx) && (spdAvailable) && (!vertexerSPD3Doff)) {
      if(fShowPlots){
        TCanvas* cxs=new TCanvas("cxs","xSPDvertex");
        cxs->cd();
        histSPDvtxX->Draw();
      }
      histSPDvtxX ->Fit(fitVtxX, "M");
      xMeanVtx = fitVtxX -> GetParameter(1);
      xSigmaVtx = fitVtxX -> GetParameter(2);
      if (TMath::Abs(xMeanVtx) > 2.) {
        xMeanVtx = 0.;
        writeMeanVertexSPD=kTRUE;
      }
      if(fShowPlots){
        TCanvas* cys=new TCanvas("cys","ySPDvertex");
        cys->cd();
        histSPDvtxY->Draw();
      }
      histSPDvtxY ->Fit(fitVtxY, "M");
      yMeanVtx = fitVtxY -> GetParameter(1);
      ySigmaVtx = fitVtxY -> GetParameter(2);
      if (TMath::Abs(yMeanVtx) > 2.) {
        yMeanVtx = 0.;
        writeMeanVertexSPD=kTRUE;
      }
      if(fShowPlots){
        TCanvas* czs=new TCanvas("czs","zSPDvertex");
        czs->cd();
        histSPDvtxZ->Draw();
      }
      histSPDvtxZ ->Fit(fitVtxZ, "M", "", -12, 12);
      zMeanVtx = fitVtxZ -> GetParameter(1);
      zSigmaVtx = fitVtxZ -> GetParameter(2);
      if ((TMath::Abs(zMeanVtx) > 20.) || (zSigmaVtx>12.)) {
        zMeanVtx = histSPDvtxZ ->GetMean();
        zSigmaVtx = histSPDvtxZ->GetRMS();
        writeMeanVertexSPD = kTRUE;
      }
    }
    else if ((useSPDvtx) && (!spdAvailable)) {
      AliError("Difference between trkVtx and online one, SPD histos not enough entry or SPD 3D vertex off. Writing Mean Vertex SPD");
      writeMeanVertexSPD = kTRUE;
    }


    delete fitVtxX;
    delete fitVtxY;
    delete fitVtxZ;
    
    // Check with online position
    Double_t posOnline[3], sigmaOnline[3];
        
    if (useTRKvtx || useITSSAvtx || writeMeanVertexSPD) {
      AliCDBManager *manCheck = AliCDBManager::Instance();
      manCheck->SetDefaultStorage("raw://");
      manCheck->SetRun(runNb);
            
      AliCDBEntry *entr = manCheck->Get("GRP/Calib/MeanVertexSPD");
      if (entr) {
        AliESDVertex *vtxOnline = (AliESDVertex*)entr->GetObject();
                
        posOnline[0] = vtxOnline->GetX();
        posOnline[1] = vtxOnline->GetY();
        posOnline[2] = vtxOnline->GetZ();
                
        sigmaOnline[0] = vtxOnline->GetXRes();
        sigmaOnline[1] = vtxOnline->GetYRes();
        sigmaOnline[2] = vtxOnline->GetZRes();
                
        AliDebug(2, Form("sigmaOnline[0] = %f, sigmaOnline[1] = %f, sigmaOnline[2] = %f", sigmaOnline[0], sigmaOnline[1], sigmaOnline[2]));
        //vtxOnline->GetSigmaXYZ(sigmaOnline);
                
        if ((TMath::Abs(posOnline[0]-xMeanVtx) > 0.1) || (TMath::Abs(posOnline[1]-yMeanVtx) > 0.1) || (TMath::Abs(posOnline[2]-zMeanVtx) > 1.)) {
          AliWarning("vertex offline far from the online one");
        }
      }
    }
    
    if (writeMeanVertexSPD) {
            
      AliWarning("Writing online Mean Vertex SPD, Mean Vertex not available");
      
      Double_t sigma[3] = { 0.0150, 0.0150, zSigmaVtx };
      
      AliESDVertex *vertex = new AliESDVertex(posOnline, sigma, "vertex");
      AliCDBId id("GRP/Calib/MeanVertex", runNb, runNb);
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(0); //check!!!!
      metaData.SetResponsible("Davide Caffarri");
      metaData.SetComment("Mean Vertex object used in reconstruction");
      
      if (!db->Put(vertex, id, &metaData)) {
        AliError(Form("Error while putting object in storage %s", db->GetURI().Data()));
      }
      
      delete vertex;
      return;
    }
    
    Bool_t highMultEnvironment = kFALSE;
    Bool_t highMultppEnvironment = kFALSE;
    Bool_t lowMultppEnvironment = kFALSE;
    
    AliCDBEntry *grpEntry = manCheck->Get("GRP/GRP/Data");
    if (!grpEntry) {
      Printf("GRP entry not found!");
      return;
    }
    
    AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(grpEntry->GetObject());
    if (!grpData) {
      Printf("GRP Data file not found");
      return;
    }
    
    TString beamType(grpData->GetBeamType());
    if (beamType == "A-A") highMultEnvironment = kTRUE;
    else {
      if ((beamType == "p-A") || (beamType == "A-p")) highMultppEnvironment = kTRUE;
      else lowMultppEnvironment = kTRUE;
    }
        
    Float_t meanMult = 38.;
    Float_t p2 = 1.3;
    Float_t resolVtx = 0.04;
        
    Double_t xSigmaMult, ySigmaMult, corrXZ, corrYZ, lumiRegSquaredX, lumiRegSquaredY;
    Double_t covarXZ=0., covarYZ=0.;
        
    TF1 *corrFit;
        
    //TH1F *histTRKdefMultX=0;
    //TH1F *histTRKdefMultY=0;
    // TH1F *histTRKHighMultX=0;
    // TH1F *histTRKHighMultY=0;
    TH2F *histTRKVertexXZ=0;
    TH2F *histTRKVertexYZ=0;
        
    TH2F *histTRKvsMultX=0x0;
    TH2F *histTRKvsMultY=0x0;
        
    if (useTRKvtx) {
      if (list) {
        //histTRKdefMultX = (TH1F*)list->FindObject("hTRKVertexXdefMult");
        //histTRKdefMultY = (TH1F*)list->FindObject("hTRKVertexYdefMult");
        // histTRKHighMultX = (TH1F*)list->FindObject("hTRKVertexXHighMult");
        // histTRKHighMultY = (TH1F*)list->FindObject("hTRKVertexYHighMult");
                
        histTRKvsMultX = (TH2F*)list->FindObject("hTRKVertexXvsMult");
        histTRKvsMultY = (TH2F*)list->FindObject("hTRKVertexYvsMult");
                
        histTRKVertexXZ = (TH2F*)list->FindObject("hTRKVertexXZ");
        histTRKVertexYZ = (TH2F*)list->FindObject("hTRKVertexYZ");
      }
      else {
        //histTRKdefMultX = (TH1F*)file->Get("hTRKVertexXdefMult");
        //histTRKdefMultY = (TH1F*)file->Get("hTRKVertexYdefMult");
        // histTRKHighMultX = (TH1F*)file->Get("hTRKVertexXHighMult");
        // histTRKHighMultY = (TH1F*)file->Get("hTRKVertexYHighMult");
                
        histTRKvsMultX = (TH2F*)file->FindObject("hTRKVertexXvsMult");
        histTRKvsMultY = (TH2F*)file->FindObject("hTRKVertexYvsMult");
                
        histTRKVertexXZ = (TH2F*)file->Get("hTRKVertexXZ");
        histTRKVertexYZ = (TH2F*)file->Get("hTRKVertexYZ");
      }
    }
        
    if (useITSSAvtx) {
      if (list) {
        //histTRKdefMultX = (TH1F*)list->FindObject("hITSSAVertexXdefMult");
        //histTRKdefMultY = (TH1F*)list->FindObject("hITSSAVertexYdefMult");
        // histTRKHighMultX = (TH1F*)list->FindObject("hITSSAVertexXHighMult");
        // histTRKHighMultY = (TH1F*)list->FindObject("hITSSAVertexYHighMult");
                
        histTRKvsMultX = (TH2F*)list->FindObject("hITSSAVertexXvsMult");
        histTRKvsMultY = (TH2F*)list->FindObject("hITSSAVertexYvsMult");
                
        histTRKVertexXZ = (TH2F*)list->FindObject("hITSSAVertexXZ");
        histTRKVertexYZ = (TH2F*)list->FindObject("hITSSAVertexYZ");
      }
            
      else {
        //histTRKdefMultX = (TH1F*)file->Get("hITSSAVertexXdefMult");
        //histTRKdefMultY = (TH1F*)file->Get("hITSSAVertexYdefMult");
        // histTRKHighMultX = (TH1F*)file->Get("hITSSAVertexXHighMult");
        // histTRKHighMultY = (TH1F*)file->Get("hITSSAVertexYHighMult");
                
        histTRKvsMultX = (TH2F*)file->FindObject("hITSSAVertexXvsMult");
        histTRKvsMultY = (TH2F*)file->FindObject("hITSSAVertexYvsMult");
                
        histTRKVertexXZ = (TH2F*)file->Get("hITSSAVertexXZ");
        histTRKVertexYZ = (TH2F*)file->Get("hITSSAVertexYZ");
      }
    }
        
    Int_t minprojbin=1;
    Int_t maxprojbin=300;
    if (highMultEnvironment) {
      minprojbin=15;
      maxprojbin=300;
    } else if (highMultppEnvironment) {
      minprojbin=10;
      maxprojbin=30;
    } else if (lowMultppEnvironment) {
      minprojbin=3;
      maxprojbin=5;
    }
    TH1D *projXvsMult = (TH1D*)histTRKvsMultX->ProjectionX("projXHighMult",minprojbin,maxprojbin);
    TH1D *projYvsMult = (TH1D*)histTRKvsMultY->ProjectionX("projYHighMult",minprojbin,maxprojbin);
    Float_t nEntriesMultX = projXvsMult->GetEffectiveEntries();
    Float_t nEntriesMultY = projYvsMult->GetEffectiveEntries();
            
    if ((nEntriesMultX <100) && (nEntriesMultY<100)) AliWarning("Few entries in High Mulitplicity histos");

    Bool_t doHighMultFits=kTRUE;
    if (lowMultppEnvironment) {            
      TH1D *projMult=(TH1D*)histTRKvsMultX->ProjectionY("projMultHighMultPbPb");
      projMult->GetXaxis()->SetRange(3,5);
      meanMult = projMult->GetMean();
      if ((projXvsMult->GetEntries() < 40.) || (projYvsMult->GetEntries() < 40.)) {
        AliWarning("Histos for lumi reg calculation too few entries, default value set");
        xSigmaVtx=0.0120;
        ySigmaVtx=0.0120;
        fStatus=kLumiRegCovMatrixProblem;
        doHighMultFits=kFALSE;
      }
    }
    
    if(doHighMultFits){
      printf("\n------ Fit high multiplicity histos ------\n");
      if (fShowPlots) {
        TCanvas *cx = new TCanvas("nwCx", "nwCx");
        cx->cd();
        projXvsMult->Draw();
      }
      TF1 *fGausX = new TF1("fGausX", "gaus");
      fGausX->SetParameters(projXvsMult->Integral(), projXvsMult->GetMean(), projXvsMult->GetRMS());
      projXvsMult -> Fit(fGausX, "M", "", projXvsMult->GetMean()-5*projXvsMult->GetRMS(),projXvsMult->GetMean()+5*projXvsMult->GetRMS());

      xSigmaMult = fGausX->GetParameter(2);
      if(lowMultppEnvironment){
        lumiRegSquaredX = (xSigmaMult*xSigmaMult - ((resolVtx*resolVtx)/TMath::Power(meanMult, p2)));
        if (lumiRegSquaredX < 0 || lumiRegSquaredX < 1E-5) {
          AliWarning("Difficult luminous region determination X, keep convoluted sigma");
          xSigmaVtx = xSigmaMult;
          fStatus=kLumiRegCovMatrixProblem;
        } else {
          if (lumiRegSquaredX > 0 && lumiRegSquaredX < 0.0005) {
            xSigmaVtx = TMath::Sqrt(lumiRegSquaredX);
            xSigmaVtx = xSigmaVtx*1.1;
          } else {
            AliWarning("Not possible to define a luminous region X. Default values set");
            xSigmaVtx = 0.0120;
            fStatus=kLumiRegCovMatrixProblem;
          }
        }
      }else{ //highMultEnvironment || highMultppEnvironment
        if ((xSigmaMult <0) || (xSigmaMult>0.03)) {
          AliWarning("Problems with luminosity region determination, update of the postion only");
          xSigmaMult = 0.;
          xSigmaVtx = 0.0120;
          fStatus=kLumiRegCovMatrixProblem;
        } else {
          xSigmaVtx = xSigmaMult;
          xSigmaVtx = xSigmaVtx*1.1;
        }
      }

      if (fShowPlots) {
        TCanvas *cy = new TCanvas("nwCy", "nwCy");
        cy->cd();
        projYvsMult->Draw();
      }
      TF1 *fGausY = new TF1("fGausY", "gaus");
      fGausY->SetParameters(projYvsMult->Integral(), projYvsMult->GetMean(), projYvsMult->GetRMS());
      projYvsMult -> Fit(fGausY, "M", "", projYvsMult->GetMean()-5*projYvsMult->GetRMS(),projYvsMult->GetMean()+5*projYvsMult->GetRMS());
                
      ySigmaMult = fGausY->GetParameter(2);
      if(lowMultppEnvironment){                
        lumiRegSquaredY = (ySigmaMult*ySigmaMult - ((resolVtx*resolVtx)/TMath::Power(meanMult, p2)));
        if (lumiRegSquaredY < 0 || lumiRegSquaredY < 1E-5) {
          AliWarning("Difficult luminous region determination Y, keep convoluted sigma");
          ySigmaVtx = ySigmaMult;
          fStatus=kLumiRegCovMatrixProblem;
        } else {
          if (lumiRegSquaredY > 0 && lumiRegSquaredY < 0.0005) {
            ySigmaVtx = TMath::Sqrt(lumiRegSquaredY);
            ySigmaVtx = ySigmaVtx*1.1;
          } else {
            AliWarning("Not possible to define a luminous region Y. Default values set");
            ySigmaVtx = 0.0120;
            fStatus=kLumiRegCovMatrixProblem;
          }       
        }
      } else { //highMultEnvironment || highMultppEnvironment
        if ((ySigmaMult <0) || (ySigmaMult>0.03)) {
          AliWarning("Problems with luminosiyy region determination, update of the postion only");
          ySigmaMult = 0.;
          ySigmaVtx = 0.0120;
          fStatus=kLumiRegCovMatrixProblem;
        } else {
          ySigmaVtx = ySigmaMult;
          ySigmaVtx = ySigmaVtx*1.1;
        }
      }
      
      delete fGausY;
      delete fGausX;

      TProfile *htrkXZ = histTRKVertexXZ ->ProfileY();
      if (fShowPlots) {
        TCanvas *cxz = new TCanvas("nwCxz", "nwCxz");
        cxz->cd();
        htrkXZ->Draw();
      }
      htrkXZ -> Fit("pol1", "M", "", -10., 10.);
      corrFit = htrkXZ->GetFunction("pol1");
      corrXZ = corrFit->GetParameter(1);
            
      if (TMath::Abs(corrXZ) > 0.01) {
        AliWarning("Problems in the correlation fitting, not update the covariance matrix");
        corrXZ =0.;
        fStatus=kLumiRegCovMatrixProblem;
      }
      else {
        covarXZ = corrXZ * zSigmaVtx*zSigmaVtx;
      }
            
      TProfile *htrkYZ = histTRKVertexYZ ->ProfileY();
      if (fShowPlots) {
        TCanvas *cyz = new TCanvas("nwCyz", "nwCyz");
        cyz->cd();
        htrkYZ->Draw();
      }
      htrkYZ -> Fit("pol1", "M", "", -10., 10.);
      corrFit = htrkYZ->GetFunction("pol1");
      corrYZ = corrFit->GetParameter(1);
      
      if (TMath::Abs(corrYZ) > 0.01) {
        AliWarning("Problems in the correlation fitting, not update the covariance matrix");
        corrYZ =0.;
        fStatus=kLumiRegCovMatrixProblem;
      }
      else {
        covarYZ = corrYZ*zSigmaVtx*zSigmaVtx;
      }
    }
    
    Double_t position[3], covMatrix[6];
    Double_t chi2=1.;
    Int_t nContr=1;
    
    position[0] = xMeanVtx;
    position[1] = yMeanVtx;
    position[2] = zMeanVtx;
    
    covMatrix[0] = xSigmaVtx*xSigmaVtx;
    covMatrix[1] = 0.; //xy
    covMatrix[2] = ySigmaVtx*ySigmaVtx;
    covMatrix[3] = covarXZ;
    covMatrix[4] = covarYZ;
    covMatrix[5] = zSigmaVtx*zSigmaVtx;
    
    //Printf ("sigma x = %f, sigma y = %f, %f, %f", xSigmaVtx, ySigmaVtx, covarXZ, covarYZ);
    
    AliESDVertex  *vertex =  new AliESDVertex(position, covMatrix, chi2, nContr, "vertex");
    
    AliCDBId id("GRP/Calib/MeanVertex", runNb, runNb);
    
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0); //check!!!!
    metaData.SetResponsible("Davide Caffarri");
    metaData.SetComment("Mean Vertex object used in reconstruction");
    
    if (!db->Put(vertex, id, &metaData)) {
      AliError(Form("Error while putting object in storage %s", db->GetURI().Data()));
      fStatus=kStoreError;
    }
        
    delete vertex;
  } // end of pass0 case
  else if (cPassMode == 1) {
    TF1* gs = new TF1("gs", "gaus", -30, 30);
    gs->SetParameters(histSPDvtxZ->GetMaximum(),histSPDvtxZ->GetMean(),histSPDvtxZ->GetRMS());
    TFitResultPtr rSPD = histSPDvtxZ->Fit(gs, "ons");
    gs->SetParameters(histTRKvtxZ->GetMaximum(),histTRKvtxZ->GetMean(),histTRKvtxZ->GetRMS());
    TFitResultPtr rTRK = histTRKvtxZ->Fit(gs, "ons");
    //
    int ndfSPD = rSPD->Ndf(), ndfTRK = rTRK->Ndf();
    //        double chiSPD = rSPD->Chi2(), chiTRK = rTRK->Chi2();
    //
    Bool_t okSPD=kFALSE,okTRK=kFALSE;
    //
    if (ndfSPD>1 && (histSPDvtxZ->GetEntries()>100)) okSPD = kTRUE;
    if (ndfTRK>1 && (histTRKvtxZ->GetEntries()>100)) okTRK = kTRUE;
    //
    if (!okSPD && !okTRK) {
      printf("Neither histos fits have convergecd\n");
      fStatus=kFitUpdateZFailed;
    }
        
    if      (!okSPD) rSPD = rTRK;
    else if (!okTRK) rTRK = rSPD;
        
    if (okTRK || okSPD) {
      ModObject("GRP/Calib/MeanVertex",rTRK->GetParams()[1],rTRK->GetParams()[2], "ZcoordUpdated", db);
      ModObject("GRP/Calib/MeanVertexSPD",rSPD->GetParams()[1],rSPD->GetParams()[2], "ZcoordUpdated", db);
    }
  } // end of pass1 case
    
  Int_t status=GetStatus();
  if (status == 0) {
    AliInfo(Form("MeanVertex calibration successful: %s (status=%d)", fgkStatusCodeName[fStatus], status));
  }
  else if (status > 0) {
    AliInfo(Form("MeanVertex calibration failed: %s (status=%d)", fgkStatusCodeName[fStatus], status));
  }
  else if (status < 0) {
    AliInfo(Form("MeanVertex calibration but not fatal error: %s (status=%d)", fgkStatusCodeName[fStatus], status));
  }
}

//__________________________________________________________________________
Int_t AliMeanVertexPreprocessorOffline::GetStatus() {
  /*
   * get status
   */
    
  switch (fStatus) {
            
    /* OK, return zero */
  case kOk:
    return 0;
    break;
            
    /* non-fatal error, return negative status */
  case kLowStatistics:
  case kWriteMeanVertexSPD:
  case kUseOfflineSPDvtx:
  case kLumiRegCovMatrixProblem:
  case kFitUpdateZFailed:
    return -fStatus;
    break;
            
    /* fatal error, return positive status */
  case kInputError:
  case kStoreError:
    return fStatus;
    break;
            
    /* anything else, return negative large number */
  default:
    return -999;
    break;
  }
    
  /* should never arrive here, anyway return negative large number */
  return -999;
}

//_______________________________________________________________________________________________
void AliMeanVertexPreprocessorOffline::ModObject(const char* url, double zv, double zs, const char* commentAdd, AliCDBStorage *db)
{
  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBEntry* entry = man->Get(url);
  AliESDVertex* vtx = (AliESDVertex*)entry->GetObject();
  double zvOld = vtx->GetZ(), zsOld = vtx->GetZRes(), covm[6];
  vtx->GetCovarianceMatrix(covm);
  covm[5] = zs*zs;
  vtx->SetZv(zv);
  vtx->SetCovarianceMatrix(covm);
  //
  printf("run %d %s %+e/%e -> %+e/%e\n",man->GetRun(),url,zvOld,zsOld,zv,zs);
    
  AliCDBMetaData* mdold = entry->GetMetaData();
  AliCDBMetaData* mdnew = new AliCDBMetaData();
  AliCDBId&idOld = entry->GetId();
  int firstRun = idOld.GetFirstRun(), lastRun = idOld.GetLastRun();
    
  TString commComb = "";
  if (mdold) {
    mdnew->SetResponsible(mdold->GetResponsible());
    mdnew->SetBeamPeriod(mdold->GetBeamPeriod());
    mdnew->SetAliRootVersion(mdold->GetAliRootVersion());
    commComb += mdold->GetComment();
  }
  commComb += commentAdd;
  mdnew->SetComment(commComb.Data());
    
  if (firstRun==0 || lastRun==999999999) {
    printf("Object %s to be modified has run range %d : %d, restricting to %d\n",
           url,firstRun,lastRun,man->GetRun());
    firstRun = lastRun = man->GetRun();
  }
    
  AliCDBId id(url,firstRun,lastRun);
  if (db) {
    db->Put(vtx,id,mdnew); //comment to run locally
  }
  else {
    man->Put(vtx,id,mdnew); //comment to run locally
  }
  //uncomment to run locally
  //AliCDBStorage* stor = man->GetStorage("local://");
  //stor->Put(vtx,id,mdnew);
    
}

