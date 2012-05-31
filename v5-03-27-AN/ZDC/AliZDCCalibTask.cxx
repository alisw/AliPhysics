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

#include <Riostream.h> 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TList.h>
#include <TH1F.h>
#include <TFile.h>

#include "AliZDCCalibTask.h"
#include "AliESDEvent.h"
#include "AliLog.h"

ClassImp(AliZDCCalibTask)

//_______________________________________________________
AliZDCCalibTask::AliZDCCalibTask(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0),	
  fListOfHistos(0x0),
  fZNCHisto(0x0),	
  fZPCHisto(0x0),
  fZNAHisto(0x0),
  fZPAHisto(0x0)	
{
  // Constructor
  DefineOutput(1, TList::Class());

}

//_______________________________________________________
AliZDCCalibTask::AliZDCCalibTask(const AliZDCCalibTask &calibTask):
  AliAnalysisTaskSE("AliZDCCalibTask"),
  fESD(calibTask.fESD),	
  fListOfHistos(calibTask.fListOfHistos),
  fZNCHisto(calibTask.fZNCHisto),	
  fZPCHisto(calibTask.fZPCHisto),
  fZNAHisto(calibTask.fZNAHisto),
  fZPAHisto(calibTask.fZPAHisto)	
{
  // Copy constructor
 
}

//______________________________________________________________________________
AliZDCCalibTask:: ~AliZDCCalibTask() 
{
  // destructor
  delete fESD;        
  delete fListOfHistos;
  delete fZNCHisto;   
  delete fZPCHisto;   
  delete fZNAHisto;   
  delete fZPAHisto;   
}

//______________________________________________________________________________
AliZDCCalibTask& AliZDCCalibTask::operator=(const AliZDCCalibTask &calibtask)  
{ 
  //assignment operator
  
  fESD = calibtask.fESD;	
  fListOfHistos = calibtask.fListOfHistos;
  fZNCHisto = calibtask.fZNCHisto;	
  fZPCHisto = calibtask.fZPCHisto;	
  fZNAHisto = calibtask.fZNAHisto;	
  fZPAHisto = calibtask.fZPAHisto; 
  
  return *this;  
}

//--------------------------------------------------------------------------
void AliZDCCalibTask::BookHistos(){

  //booking histos

  AliInfo(Form("***  Booking Histograms %s", GetName())) ; 
  
  fZNCHisto = new TH1F("fZNCHisto", "Energy in ZNC",100,0.,MAXENVALUE);
  fZNCHisto->SetXTitle("E (GeV)");
  fZPCHisto = new TH1F("fZPCHisto", "Energy in ZPC",100,0.,MAXENVALUE);
  fZPCHisto->SetXTitle("E (GeV)");
  fZNAHisto = new TH1F("fZNAHisto", "Energy in ZNA",100,0.,MAXENVALUE);
  fZNAHisto->SetXTitle("E (GeV)");
  fZPAHisto = new TH1F("fZPAHisto", "Energy in ZPA",100,0.,MAXENVALUE);
  fZPAHisto->SetXTitle("E (GeV)");
  
  if(!fListOfHistos) fListOfHistos = new TList();
  fListOfHistos->SetOwner();
  fListOfHistos->AddAt(fZNCHisto, 0);
  fListOfHistos->AddAt(fZPCHisto, 1);
  fListOfHistos->AddAt(fZNAHisto, 2);
  fListOfHistos->AddAt(fZPAHisto, 3);
  
}

//--------------------------------------------------------------------------
void AliZDCCalibTask::DrawHistos(){

  //drawing output histos
  
  AliInfo(Form("*** Drawing Histograms %s", GetName())) ; 

  if(!gROOT->IsBatch()){
    fListOfHistos = dynamic_cast<TList*>(GetOutputData(1));
    fZNCHisto = (TH1F*) fListOfHistos->At(0);
    fZPCHisto = (TH1F*) fListOfHistos->At(1);
    fZNAHisto = (TH1F*) fListOfHistos->At(2);
    fZPAHisto = (TH1F*) fListOfHistos->At(3);
    
    TCanvas * canvPlot = new TCanvas("canvPlot","ZDC energy",0,0,600,600);
    canvPlot->Divide(2,2);
    canvPlot->cd(1);
    fZNCHisto->SetLineColor(kTeal+1);
    fZNCHisto->DrawCopy("hist");
    canvPlot->cd(2);
    fZPCHisto->SetLineColor(kAzure+1);
    fZPCHisto->DrawCopy("hist");
    canvPlot->cd(3);
    fZNCHisto->SetLineColor(kGreen+1);
    fZNCHisto->DrawCopy("hist");
    canvPlot->cd(4);
    fZPCHisto->SetLineColor(kBlue+1);
    fZPCHisto->DrawCopy("hist");
    
  }
  
}

//_______________________________________________________
void AliZDCCalibTask::UserCreateOutputObjects()
{
  // Create histograms

  if(fDebug>1) printf("ZDCCalibTask::UserCreateOutputObjects() \n");

  //Open histograms
  OpenFile(1);
  BookHistos();
  
}

//_______________________________________________________
void AliZDCCalibTask::UserExec(Option_t */*option*/)
{
    // Execute analysis    
    AliVEvent *event = InputEvent();
    if (!event) {
      printf("ERROR: Could not retrieve event\n");
      return;
    }
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    AliESDZDC *esdZDC = esd->GetESDZDC();

    Float_t *znctow = (Float_t *) (esdZDC->GetZN1TowerEnergyLR());    
    Float_t *zpctow = (Float_t *) (esdZDC->GetZP1TowerEnergyLR());    
    Float_t *znatow = (Float_t *) (esdZDC->GetZN2TowerEnergyLR());    
    Float_t *zpatow = (Float_t *) (esdZDC->GetZP2TowerEnergyLR());    
    
    Float_t energyZNC=0., energyZPC=0.;
    Float_t energyZNA=0., energyZPA=0.;
    for(Int_t j=0; j<5; j++){
      energyZNC += znctow[j];
      energyZPC += zpctow[j];
      energyZNA += znatow[j];
      energyZPA += zpatow[j];
    }
    
    //printf(" Filling histos with values: %f  %f  %f  %f\n",
    //	energyZNC,energyZPC,energyZNA,energyZPA);
    //
    fZNCHisto->Fill(energyZNC);
    fZPCHisto->Fill(energyZPC);
    fZNAHisto->Fill(energyZNA);
    fZPAHisto->Fill(energyZPA);
    
    PostData(1, fListOfHistos);
}

//_______________________________________________________
void AliZDCCalibTask::Terminate(Option_t *)
{
  // Draw results
  TH1::AddDirectory(0);

  fListOfHistos = dynamic_cast<TList*>(GetOutputData(1));
  fZNCHisto = (TH1F*) fListOfHistos->At(0);
  fZPCHisto = (TH1F*) fListOfHistos->At(1);
  fZNAHisto = (TH1F*) fListOfHistos->At(2);
  fZPAHisto = (TH1F*) fListOfHistos->At(3);
  
  Int_t binMax[4], nBinsx[4];
  Float_t yMax[4], meanFitVal[4]={0.,0.,0.,0.};
  TF1 *fitfun[4];

  /*printf("\n # of entries in histos: %d  %d  %d  %d \n",
  	(Int_t) fZNCHisto->GetEntries(),(Int_t) fZPCHisto->GetEntries(),
	(Int_t) fZNAHisto->GetEntries(),(Int_t) fZPAHisto->GetEntries());*/

  if(fZNCHisto->GetEntries() != 0){ 
    binMax[0] = fZNCHisto->GetMaximumBin();
    yMax[0]   = (fZNCHisto->GetXaxis())->GetXmax();
    nBinsx[0] = (fZNCHisto->GetXaxis())->GetNbins();
    fZNCHisto->Fit("gaus","Q","",binMax[0]*yMax[0]/nBinsx[0]*0.7,binMax[0]*yMax[0]/nBinsx[0]*1.25);
    fitfun[0] = fZNCHisto->GetFunction("gaus");
    if(fitfun[0]) meanFitVal[0] = (Float_t) (fitfun[0]->GetParameter(1));
    else printf(" !!! No fit on fZNCHisto !!!\n");
  } 
  else AliWarning(" !!! ZNC empty histo !!!\n");

  if(fZPCHisto->GetEntries() != 0){ 
    binMax[1] = fZPCHisto->GetMaximumBin();
    yMax[1]   = (fZPCHisto->GetXaxis())->GetXmax();
    nBinsx[1] = (fZPCHisto->GetXaxis())->GetNbins();
    fZPCHisto->Fit("gaus","Q","",binMax[1]*yMax[1]/nBinsx[1]*0.7,binMax[1]*yMax[1]/nBinsx[1]*1.25);
    fitfun[1] = fZPCHisto->GetFunction("gaus");
    if(fitfun[1]) meanFitVal[1] = (Float_t) (fitfun[1]->GetParameter(1));
    else printf(" !!! No fit on fZPCHisto !!!\n");
  } 
  else AliWarning(" !!! ZPC empty histo !!!\n");

  if(fZNAHisto->GetEntries() != 0){ 
    binMax[2] = fZNAHisto->GetMaximumBin();
    yMax[2]   = (fZNAHisto->GetXaxis())->GetXmax();
    nBinsx[2] = (fZNAHisto->GetXaxis())->GetNbins();
    fZNAHisto->Fit("gaus","Q","",binMax[2]*yMax[2]/nBinsx[2]*0.7,binMax[2]*yMax[2]/nBinsx[2]*1.25);
    fitfun[2] = fZNAHisto->GetFunction("gaus");
    if(fitfun[2]) meanFitVal[2] = (Float_t) (fitfun[2]->GetParameter(1));
    else printf(" !!! No fit on fZNAHisto !!!\n");
  } 
  else AliWarning(" !!! ZNA empty histo !!!\n");

  if(fZPAHisto->GetEntries() != 0){ 
    binMax[3] = fZPAHisto->GetMaximumBin();
    yMax[3]   = (fZPAHisto->GetXaxis())->GetXmax();
    nBinsx[3] = (fZPAHisto->GetXaxis())->GetNbins();
    fZPAHisto->Fit("gaus","Q","",binMax[3]*yMax[3]/nBinsx[3]*0.7,binMax[3]*yMax[3]/nBinsx[3]*1.25);
    fitfun[3] = fZPAHisto->GetFunction("gaus");
    if(fitfun[3]) meanFitVal[3] = (Float_t) (fitfun[3]->GetParameter(1));
    else printf(" !!! No fit on fZPAHisto !!!\n");
  } 
  else AliWarning(" !!! ZPA empty histo !!!\n");
   
  FILE *fileOut = fopen("EMDCalibData.dat","w");
  for(Int_t j=0; j<6; j++){
     if(j<4) fprintf(fileOut,"\t%f\n",meanFitVal[j]);
     else    fprintf(fileOut,"\t1.0000\n");
  }
  fclose(fileOut);

  DrawHistos();
  
    

}
