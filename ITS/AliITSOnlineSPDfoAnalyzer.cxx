/**************************************************************************
 * Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

////////////////////////////////////////////////////////////
// Author: A. Mastroserio                                 //
// This class is used in the detector algorithm framework //
// to process the data stored in special container files  //
// (see AliITSOnlineSPDfo). The "good" set of DAC values  //
// is extracted.                                          //
////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TArrayI.h>
#include <TCanvas.h>
#include <THnSparse.h>
#include <TKey.h>
#include <iostream>
#include <fstream>
#include "AliITSOnlineSPDfoChipConfig.h"
#include "AliITSOnlineSPDfoChip.h"
#include "AliITSOnlineSPDfoInfo.h"
#include "AliITSOnlineSPDfo.h"
#include "AliITSOnlineSPDfoAnalyzer.h"
#include "AliLog.h"

AliITSOnlineSPDfoAnalyzer::AliITSOnlineSPDfoAnalyzer(const TString fileName, Bool_t readFromGridFile) :
  fFileName(0),
  fNdims(0),
  fNbins(0x0),
  fXmin(0x0),
  fXmax(0x0),
  fFOHandler(0x0),
  fHighOccupancyCheck(kTRUE)
{
  // constructor
  fFileName = fileName;
  for(Int_t iqual =0; iqual<kNqualityFlags; iqual++) fGeneralThresholds[iqual]=0;
  
  for (UInt_t chipNr=0; chipNr<10; chipNr++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for(Int_t i =0; i< kNqualityFlags; i++) fNh[i][hs][chipNr]=NULL;
    }
  }
  
  Init(readFromGridFile); 
}
//-------------------------------------------------------------------
AliITSOnlineSPDfoAnalyzer::AliITSOnlineSPDfoAnalyzer(const AliITSOnlineSPDfoAnalyzer& foan) :
  fFileName(foan.fFileName),
  fNdims(foan.fNdims),
  fNbins(foan.fNbins),
  fXmin(foan.fXmin),
  fXmax(foan.fXmax),
  fFOHandler(foan.fFOHandler),
  fHighOccupancyCheck(foan.fHighOccupancyCheck)
{
  //
  // copy constructor, only copies the filename and params (not the processed data
  //
  
  for(Int_t iqual =0; iqual<3; iqual++) fGeneralThresholds[iqual] =foan.fGeneralThresholds[iqual];
  for (UInt_t chipNr=0; chipNr<10; chipNr++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for(Int_t i=0; i<kNqualityFlags;i++) fNh[i][hs][chipNr]=NULL;
    }
  }
  
  Init();
}
//-------------------------------------------------------------------
AliITSOnlineSPDfoAnalyzer::~AliITSOnlineSPDfoAnalyzer() 
{
  //
  // destructor
  // 
  
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chipNr=0; chipNr<10; chipNr++) {
      for(Int_t i=0; i<kNqualityFlags ; i++ ) if(fNh[i][hs][chipNr]!=NULL) delete fNh[i][hs][chipNr];          
    }
  }
  delete fFOHandler;
}
//-------------------------------------------------------------------
AliITSOnlineSPDfoAnalyzer& AliITSOnlineSPDfoAnalyzer::operator=(const AliITSOnlineSPDfoAnalyzer& foan) 
{
  // assignment operator, only copies the filename and params (not the processed data)
  if (this!=&foan) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chipNr=0; chipNr<10; chipNr++) {
	for(Int_t i=0; i<kNqualityFlags ; i++ ) if(fNh[i][hs][chipNr]!=NULL) delete fNh[i][hs][chipNr];
      }
    }
    
    fFileName=foan.fFileName;
    fHighOccupancyCheck=foan.fHighOccupancyCheck;
    Init();    
  }
  return *this;
}
//-------------------------------------------------------------------
void AliITSOnlineSPDfoAnalyzer::Init(Bool_t readFromGridFile) 
{
  //
  // first checks type of container and then initializes container obj
  //
  if (!readFromGridFile) {
    TFile* p0 = TFile::Open(fFileName.Data());
    if (p0 == NULL) {
      printf("no file open!"); 
      return;
    }
    else { 
      fFOHandler = new AliITSOnlineSPDfo();
      fFOHandler->SetFile(fFileName);
    }
    p0->Close();   
  }
}
//-------------------------------------------------------------------
void AliITSOnlineSPDfoAnalyzer::SetGeneralThresholds(Float_t thre[3])
{
  // here the settings for the 3 quality flags are defined
  for(Int_t i=0; i< 3; i++) fGeneralThresholds[i]=thre[i];  
}
//-------------------------------------------------------------------
void AliITSOnlineSPDfoAnalyzer::SetNdimensions()
{
  //
  // here the axis of the N-dim histograms are setted 
  // 
  if(!fFOHandler) {
    printf("no fo object. Exiting...\n"); 
    return;
  }
  
  TArrayI array = fFOHandler->GetDACscanParams();
  fNdims = array.GetSize()/3;
  fNbins = new Int_t[fNdims];
  fXmin = new Double_t[fNdims];
  fXmax = new Double_t[fNdims];
  for(Int_t i=0; i< fNdims; i++){
    fXmin[i] = array.At(3*i)-0.25;
    fXmax[i] = array.At(3*i+1)+0.25;
    fNbins[i] = (Int_t)((fXmax[i]-fXmin[i])/0.5)+1;//to avoid Int->Double conversion problems when checking results. 
  } 
}
//-------------------------------------------------------------------
void AliITSOnlineSPDfoAnalyzer::BuildTHnSparse(Int_t ihs, Int_t ichip)
{
  //
  // here the N-dim histogram is booked per chip in one HS 
  //
  
  if(!fNdims) SetNdimensions();
  
  for(Int_t iqual =0; iqual < 3; iqual++) fNh[iqual][ihs][ichip] = new THnSparseI(Form("h%i_HS%i_C%i",iqual,ihs,ichip), Form("h%i_HS%i_C%i",iqual,ihs,ichip), fNdims, fNbins, fXmin, fXmax);
}
//-------------------------------------------------------------------
Int_t AliITSOnlineSPDfoAnalyzer::Select(const AliITSOnlineSPDfoChip *chip) const
{
  //
  // Selects the DAC values if in the chip: the I configuration corresponds to
  // 0% efficiency ( no FO response case). All the others shoud be within the 
  // predefined thresholds
  
  if(!fFOHandler->GetFOscanInfo()) {
    printf("no general information object in the file. Exiting...\n");
    return -1;
  }
  
  Double_t npulses = (fFOHandler->GetFOscanInfo())->GetNumTriggers();
  
  if(npulses == 0. ) {
    Info("AliITSOnlineSPDfoAnalyzer::Select","no trigger pulses set. Exiting...");
    return -999;
  } 
  
  TObjArray *array = chip->GetChipConfigInfo();
  if(!array) {
    printf("No measurement array found in the chip!!\n");
    return 0;
  }
  
  Int_t quality = -1;
  
  Float_t counts = 0;
  
  Int_t processedconfigurations = chip->GetNumberOfChipConfigs();
  
 
  
  
  for(Int_t isteps =0; isteps < processedconfigurations; isteps++){ 
    
    Int_t matrixId = ((AliITSOnlineSPDfoChipConfig*)array->At(isteps))->GetChipConfigMatrixId();
    counts = (Float_t)(((AliITSOnlineSPDfoChipConfig*)array->At(isteps))->GetChipConfigCounter());
    if(matrixId==0 && counts > 0) return -1;
    if(fHighOccupancyCheck &&  matrixId ==6) continue;
    
    Float_t efficiency = counts/npulses;
    
    if(matrixId > 0){
      Int_t response = IsSelected(efficiency);
      if( response >=0) {
	if(quality < response) quality = response;
      }
      else return -1; 
    }
  }
  return quality;
}
//-----------------------------------------------------
Int_t AliITSOnlineSPDfoAnalyzer::IsSelected(Float_t eff) const
{
  //
  // returns the quality of the selection 
  //
  
  for(Int_t i=0; i<3; i++){  
    if(eff <= 1.+ fGeneralThresholds[i] && eff >= 1. - fGeneralThresholds[i]  ) return i;
  }
  return -1;
}
//----------------------------------------------------
void AliITSOnlineSPDfoAnalyzer::Process()
{ 
  //
  // The procedure is the following:
  // - DAC 4-tuples are checked  
  // - if the 4-tuple survives the selection, the chip-related histograms are filled.
  // (- Per each histogram the mean values of each axis are taken)
  //
  
  if(!fFOHandler) { 
    Warning("AliITSOnlineSPDfoAnalyzer::Process","no fo object. Exiting.. \n");
    return;
  } 
  
  TKey *key;
  TIter iter((fFOHandler->GetFile())->GetListOfKeys());  
  while ((key = (TKey*)(iter.Next()))) {
    TString classname = key->GetClassName();
    if(classname.Contains("OnlineSPD")) break;
    
    TObjArray *array = (TObjArray*)(fFOHandler->GetFile())->Get(key->GetName()); // array of chips corresponding to the DACS (size 1-60)
    if(!array){
      printf("no array found! Exiting...\n");
      break; 
    }
    
    Double_t *dacvalues = fFOHandler->GetDACvaluesD(key->GetName(), GetFOHandler()->GetFOscanInfo()->GetNumDACindex());
    
    for(Int_t i=0; i< array->GetSize(); i++){
      AliITSOnlineSPDfoChip *chip = (AliITSOnlineSPDfoChip *)array->At(i); 
          
      if(!chip) continue;
      Int_t hs = chip->GetActiveHS();
      Int_t chipid = chip->GetChipId();     
      Int_t quality = Select(chip);
      if(quality<0) continue; 
      if(!fNh[quality][hs][chipid]) BuildTHnSparse(hs,chipid);
      fNh[quality][hs][chipid]->Fill(dacvalues);       
    } 

    if(dacvalues) delete [] dacvalues;
  } 
}
//---------------------------------------------
void AliITSOnlineSPDfoAnalyzer::WriteToFile(TString outputfile)
{
  //
  // The 4-dim histograms are stored into a file
  //
  
  TFile * f = TFile::Open(outputfile.Data(),"recreate");
  for(Int_t ihs =0; ihs < 6; ihs++) {
    for(Int_t ichip =0; ichip < 10; ichip++){
      for(Int_t i=0; i<kNqualityFlags ; i++ ) {
      if(fNh[i][ihs][ichip]) f->WriteObjectAny(fNh[i][ihs][ichip],"THnSparse",Form("h%i_hs%i_chip%i",i,ihs,ichip));
      }
      
    }
  }
  f->Close();
}
//---------------------------------------------
void AliITSOnlineSPDfoAnalyzer::CheckResults(TString filename, Int_t hs, Int_t ichip, Int_t iqual) const
{
  //  
  //The chip related 4-dim histograms are produced and stored into eps files
  //    
  
  TFile *f = TFile::Open(filename.Data());
  if(!f) {
    Info("AliITSOnlineSPDfoAnalyzer::CehckResults","no file open!. Exiting...\n");
    return; 
  }
  
  THnSparse *hn;
  TObject *obj;
  
  TIter iter(f->GetListOfKeys());
  while((obj=iter.Next())){
    TString name = obj->GetName();
    hn=(THnSparse*)f->Get(name.Data()); 
    if(name.Contains(Form("hs%i",hs)) && name.Contains(Form("chip%i",ichip)) && name.Contains(Form("h%i",iqual))) GetCanvases(hn,hs,ichip,iqual);
  } 
}
//-----------------------------------------------------------------------------------------------
void AliITSOnlineSPDfoAnalyzer::GetCanvases(const THnSparse *hn,Int_t hs, Int_t chip, Int_t iqual) const
{
  //
  // 1-Dim and 2 Dim Projections of 4-dim histograms are visualized in canvases per quality selection
  //
  //
  
  if(!hn) {printf("no thnsparse...exiting!\n"); return;} 
  
  gStyle->SetPalette(1);
  
  TString dacname[4];
  
  if(  (fFOHandler->GetFOscanInfo())->GetDACindex(0) == 20 ) dacname[0] = "FOPOL";
  else dacname[0] = "possible DAC-name/ DAC-register-number mismatch";  
  if(  (fFOHandler->GetFOscanInfo())->GetDACindex(1) == 17 ) dacname[1] = "CONVPOL";
  else dacname[1] = "possible DAC-name/ DAC-register-number mismatch";
  if(  (fFOHandler->GetFOscanInfo())->GetDACindex(2) == 16 ) dacname[2] = "COMPREF";
  else dacname[2] = "possible DAC-name/ DAC-register-number mismatch";
  if(  (fFOHandler->GetFOscanInfo())->GetDACindex(3) == 39 ) dacname[3] = "pre_VTH";
  else dacname[3] = "possible DAC-name/ DAC-register-number mismatch";
  
  TString titles = Form("|eff-1| <= %f   for CHIP %i  in HS  %i",fGeneralThresholds[iqual],hs,chip);
  TCanvas *c[3];
  
  for(Int_t i=0; i<2; i++) {
    c[i] = new TCanvas(Form("c%i",i+1),Form(" %i DIM plots. %s ",i+1,titles.Data()),1200,800); 
    c[i]->Divide(2,2);
  }
  
  
  for(Int_t idim =0; idim<2; idim++){
    for(Int_t k=1; k<5; k++){
      c[idim]->cd(k);
      
      TH1D *h1 =0x0;
      TH2D *h2=0x0;       
      if(idim == 0) {
        h1 = hn->Projection(k-1);
        if(!h1) {
         printf("no histogram!!...\n\n\n");
        } else {
        h1->SetXTitle(Form("DAC %i  ( %s )",k-1,dacname[k-1].Data()));
        h1->SetYTitle("entries (eff within thresholds)");
        h1->Draw();
        }
      } 
      
      if(idim==1) {      
        if(k<4){
          h2  = hn->Projection(k-1,k);
          h2->SetXTitle(Form("DAC %i ( %s )",k,dacname[k].Data())); h2->SetYTitle(Form("DAC %i  ( %s )",k-1,dacname[k-1].Data()));
          h2->SetTitleOffset(2,"Y"); h2->SetTitleOffset(1.5,"X");
	} else {
	  h2 = hn->Projection(0,3);
	  h2->SetXTitle(Form("DAC %i ( %s )",3,dacname[3].Data())); h2->SetYTitle(Form("DAC %i  ( %s )",0, dacname[0].Data()));
	  h2->SetTitleOffset(2,"Y"); h2->SetTitleOffset(1.5,"X"); 
	}    
	
        h2->Draw("lego2");
      }  
    }// k loop
    
    c[idim]->SaveAs(Form("c%iDIM_%i_%i_%i.eps",idim,iqual,hs,chip));   
  }// idim loop  
}
//-----------------------------------------------------
TArrayI AliITSOnlineSPDfoAnalyzer::ChooseDACValues(Int_t hs, Int_t chip) const
{
  //
  // here the "best" 4 dacs are chosen. The most present are taken. 
  // If the n-tuple does not correspond to a meaningful entry, the
  // closest value to the mean point in the n-dimensional histogram
  // is taken.
  
  TH1D *tmp[5];
  for(Int_t jj=0;jj<5;jj++)tmp[jj]=NULL;
  if(fNdims > 5) printf("AliITSOnlineSPDfoAnalyzer::ChooseDACValues -> N. of dimensions are more than expected! Break! \n");
  TArrayI dacs(fNdims+1);
  
  for(Int_t i=0; i<fNdims+1; i++) dacs.AddAt(-1,i);
  
  for(Int_t iqual =0; iqual < kNqualityFlags; iqual++){
    if(!fNh[iqual][hs][chip]) continue;    
    if(fNh[iqual][hs][chip]->GetEntries()==0) continue;
    for(Int_t idim =0; idim<fNdims; idim++){
      if(dacs.At(idim)>=0) continue;
      tmp[idim] = fNh[iqual][hs][chip]->Projection(idim);
      dacs.AddAt((Int_t)tmp[idim]->GetBinLowEdge(tmp[idim]->GetMaximumBin()+1),idim);
      Int_t bin=-1;
      if(fFOHandler->GetFOscanInfo()->GetDACindex(idim)==fFOHandler->kIdPreVTH && CorrectPreVTHChioce(tmp[idim],bin)) {
       dacs.AddAt((Int_t)tmp[idim]->GetBinLowEdge(bin+1),idim);
      }
      dacs.AddAt(iqual,fNdims);
    }//idim
  }//iqual
  
  if(!IsExisting(dacs,hs,chip)  && dacs.At(fNdims)>-1) {   
   TArrayI centraldacs = GetCentralDACS(dacs.At(fNdims),hs,chip,tmp);
    dacs = centraldacs;
  }
  return dacs; 
}
//_____________________________________________________
Bool_t AliITSOnlineSPDfoAnalyzer::IsExisting(TArrayI dacs,Int_t hs, Int_t chip) const
{
  //
  // check the N-tuple corresponds to a real one
  //  
  
  if(dacs.At(fNdims)<0) return kFALSE;
  Bool_t isOk = kFALSE;
  
  Int_t size = dacs.GetSize()-1;
  Double_t *entry = new Double_t[size];
  for(Int_t i=0; i<size; i++) entry[i] = dacs.At(i);
  Int_t checkbin = fNh[dacs.At(dacs.GetSize()-1)][hs][chip]->GetBin(entry,kFALSE); // kFALSE does not allocate another bin
  if(checkbin > -1) isOk = kTRUE;
  delete [] entry;
  return isOk; 
}
//-----------------------------------------------------------
TArrayI AliITSOnlineSPDfoAnalyzer::GetCentralDACS(Int_t qualityflag, Int_t hs, Int_t chip, TH1D **hd) const
{
  //
  // This method gets the DAC values which are closest to the mean point in the N-dim histogram
  // It is called when the most probable DAC set does not correspond to a real entry in the N-dim
  // histogram.
  //
 
  TArrayI dacs(fNdims+1);
   
  Double_t *mean=new Double_t[fNdims];
  Int_t *goodbins=new Int_t[fNdims];
  Double_t distance = 9999999;
    for(Int_t i=0; i<fNdims ;i++){ 
    mean[i]=hd[i]->GetMean();
    goodbins[i]=0;
    dacs.AddAt(-1,i);
  }
  
  Int_t *bins = new Int_t[fNdims];
  Double_t *val=new Double_t[fNdims];
  
  for(Int_t in=0; in< fNh[qualityflag][hs][chip]->GetNbins() ; in++){
    
    fNh[qualityflag][hs][chip]->GetBinContent(in,bins);
    Double_t r2 = 0;  
    for(Int_t j=0; j<fNdims; j++) {
     val[j] = hd[j]->GetBinCenter(bins[j]);
      r2+=TMath::Power(val[j]-mean[j],2);
    }
    
    if(r2<distance) {
      distance = r2;
      fNh[qualityflag][hs][chip]->GetBinContent(in,goodbins);
    }    
  }
  
  
  for(Int_t k=0; k<fNdims; k++) {  
   dacs.AddAt((Int_t)(hd[k]->GetBinCenter(goodbins[k]) + 0.5*hd[k]->GetBinWidth(goodbins[k])),k);
  }
 
  dacs.AddAt(qualityflag,fNdims);
  
  
  delete [] mean;
  delete [] goodbins;
  delete [] bins;
  delete [] val;
  return dacs;
}
//-------------------------------------------------------
void AliITSOnlineSPDfoAnalyzer::ReadParamsFromLocation(const Char_t *dirName) 
{
  //
  // opens file (default name) in dir dirName and reads parameters from it
  // The file should be in database
  //
  
  TString paramsFileName = Form("./%s/focalib_params.txt",dirName);
  
  ifstream paramsFile;
  paramsFile.open(paramsFileName, ifstream::in);
  if (paramsFile.fail()) {
    printf("No config file (%s) present. Using default tuning parameters.\n",paramsFileName.Data());
  }
  else {
    while(1) {
      Char_t paramN[50];
      Char_t paramV[50];
      paramsFile >> paramN;
      if (paramsFile.eof()) break;
      paramsFile >> paramV;
      SetParam(paramN,paramV);
      if (paramsFile.eof()) break;
    }
    paramsFile.close();
  }
}
//---------------------------------------------------------
void AliITSOnlineSPDfoAnalyzer::SetParam(const Char_t *pname, const Char_t *pval) 
{
  //
  // sets a single parameter when reading from a file in the database
  //
  
  TString name = pname;
  TString val = pval;
  
  
  if (name.CompareTo("fGeneralThresholds0")==0) {
    fGeneralThresholds[0] = val.Atof();
  }
  else if (name.CompareTo("fGeneralThresholds1")==0) {
    fGeneralThresholds[1] = val.Atof();
  }
  else if (name.CompareTo("fGeneralThresholds2")==0) {
    fGeneralThresholds[2] = val.Atof();
  }
  else {
    Error("AliITSOnlineSPDscanAnalyzer::SetParam","Parameter %s in configuration file unknown.",name.Data());
  }
}
//--------------------------------------------------------
Bool_t AliITSOnlineSPDfoAnalyzer::CorrectPreVTHChioce(const TH1D *h,Int_t &bin) const
{
  //
  // Checks if more maxima of the same height are present in the pre_VTH case
  //
  
  
  Int_t maxbin = (Int_t)h->GetMaximumBin();
  Double_t maxentries = h->GetBinContent(maxbin);
  
  Int_t binindex = -1;
  Bool_t check=kFALSE;

  for(Int_t i=0; i< h->GetNbinsX(); i++){    
     if(h->GetBinContent(i) == maxentries){
      if(binindex <= i) binindex =i;
    }
  }
  
  
  if(binindex>-1) {
    bin=binindex;
    check = kTRUE; 
  }
  
  
  return check; 
}
