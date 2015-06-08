#if !defined( __CINT__) || defined(__MAKECINT__)


#include <Riostream.h>
#include <TSystem.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TStyle.h>

#include "AliTestShuttle.h"
#include "AliShuttleInterface.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"


#include "../TRD/AliTRDarrayF.h"
#include "../TRD/AliTRDCalibPadStatus.h"
#include "../TRD/Cal/AliTRDCalPadStatus.h"
#include "../TRD/Cal/AliTRDCalDet.h"
#include "../TRD/Cal/AliTRDCalPad.h"
#include "../TRD/Cal/AliTRDCalROC.h"


#endif


void AliTRDCheckPreprocessorold()
{
  // load library
  //gSystem->Load("libTestShuttle");


  AliTestShuttle::SetMainCDB("local://TestCDB");
  AliTestShuttle::SetMainRefStorage("local://TestReference");
  printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());

  //Style
  //************************
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

  //Check the reference data 
  //***************************

  //Test reference data gain HLT
  //***************************
  Int_t ErrorRefDataGainHLT = 0;
  AliCDBEntry* entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainRefStorage())->Get("TRD/HLTData/Gain", 8);  
  if(!entry) ErrorRefDataGainHLT = 1;
  else{
    TH2I *histogainhlt = (TH2I *) entry->GetObject();
    if(!histogainhlt) ErrorRefDataGainHLT = 2;
    else{
      Int_t NbinsX = histogainhlt->GetNbinsY();
      if(NbinsX != 540) ErrorRefDataGainHLT = 3;
      TCanvas *cgainhlt = new TCanvas("cgainhlt","",50,50,600,800);
      cgainhlt->cd();
      histogainhlt->Draw("LEGO");
    }
  }


  //Test reference data vdriftt0 HLT
  //***************************
  Int_t ErrorRefDataVdriftT0HLT = 0;
  if(entry) delete entry;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainRefStorage())->Get("TRD/HLTData/VdriftT0", 8);
  if(!entry) ErrorRefDataVdriftT0HLT = 1;
  else{
    TProfile2D *histovdriftt0hlt = (TProfile2D *) entry->GetObject();
    if(!histovdriftt0hlt) ErrorRefDataVdriftT0HLT = 2;
    else{
      Int_t NbinsX = histovdriftt0hlt->GetNbinsY();
      if(NbinsX != 540) ErrorRefDataVdriftT0HLT = 3;
      TCanvas *cvdrifthlt = new TCanvas("cvdrifthlt","",50,50,600,800);
      cvdrifthlt->cd();
      histovdriftt0hlt->Draw("LEGO");
    }
  }
  
  
  //Test reference data PRF HLT
  //***************************
  Int_t ErrorRefDataPRFHLT = 0;
  if(entry) delete entry;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainRefStorage())->Get("TRD/HLTData/PRF", 8);
  if(!entry) ErrorRefDataPRFHLT = 1;
  else{
    TProfile2D *histoprfhlt = (TProfile2D *) entry->GetObject();
    if(!histoprfhlt) ErrorRefDataPRFHLT = 2;
    else{
      Int_t NbinsX = histoprfhlt->GetNbinsY();
      if(NbinsX != 540) ErrorRefDataPRFHLT = 3;
      TCanvas *cprfhlt = new TCanvas("cprfhlt","",50,50,600,800);
      cprfhlt->cd();
      histoprfhlt->Draw("LEGO");
    }
  }
  

  //Test reference data PadStatus DAQ
  //***************************
  Int_t ErrorRefDataPadStatus = 0;

  Int_t nbsm = 0;
  
  for(Int_t k = 0; k < 18; k++){

    TString padstatus("TRD/DAQData/PadStatus");
    padstatus += k;
    
    if(entry) delete entry;
    entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainRefStorage())->Get((const char*)padstatus, 8);
    if(entry){
      AliTRDCalibPadStatus *calpadstatus = (AliTRDCalibPadStatus *) entry->GetObject();
      if(!calpadstatus) ErrorRefDataPadStatus = 2;
      else{
	//Make the AliTRDCalDet correspondant
	AliTRDCalDet calDet = AliTRDCalDet("dummy","dummy for mean");
	for(Int_t l = 0; l < 540; l++){
	  calDet.SetValue(l,10.0);
	}
	
	//Make the AliTRDCalPad correspondant
	AliTRDCalPad calPad1 = AliTRDCalPad("meanpad","dummy for mean");
	AliTRDCalPad calPad2 = AliTRDCalPad("squarespad","dummy for squares");
	AliTRDCalROC *calROC1 = 0x0;
	AliTRDCalROC *calROC2 = 0x0;
	for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)
	  {
	    AliTRDCalROC *calROC11 = calPad1.GetCalROC(det);
	    AliTRDCalROC *calROC22 = calPad2.GetCalROC(det);
	    calROC1                = calpadstatus->GetCalRocMean(det,kTRUE);
	    calROC2                = calpadstatus->GetCalRocRMS(det,kTRUE);
	    for(Int_t k = 0; k < calROC22->GetNchannels(); k++){
	      calROC11->SetValue(k,calROC1->GetValue(k));
	      calROC22->SetValue(k,calROC2->GetValue(k));
	    }
	  }
	TCanvas *cpadstatusm = new TCanvas((const char*)padstatus,(const char*)padstatus,50,50,600,800);
	cpadstatusm->Divide(3,2);
	cpadstatusm->cd(1);
	((TH2F *)calPad1.MakeHisto2DSmPl(k,0,&calDet,0,9.0,11.0,-1))->Draw("colz");
	cpadstatusm->cd(2);
	((TH2F *)calPad1.MakeHisto2DSmPl(k,1,&calDet,0,9.0,11.0,-1))->Draw("colz");
	cpadstatusm->cd(3);
	((TH2F *)calPad1.MakeHisto2DSmPl(k,2,&calDet,0,9.0,11.0,-1))->Draw("colz");
	cpadstatusm->cd(4);
	((TH2F *)calPad1.MakeHisto2DSmPl(k,3,&calDet,0,9.0,11.0,-1))->Draw("colz");
	cpadstatusm->cd(5);
	((TH2F *)calPad1.MakeHisto2DSmPl(k,4,&calDet,0,9.0,11.0,-1))->Draw("colz");
	cpadstatusm->cd(6);
	((TH2F *)calPad1.MakeHisto2DSmPl(k,5,&calDet,0,9.0,11.0,-1))->Draw("colz");
	
	padstatus += 1982;

	TCanvas *cpadstatusrms = new TCanvas((const char*)padstatus,(const char*)padstatus,50,50,600,800);
	cpadstatusrms->Divide(3,2);
	cpadstatusrms->cd(1);
	((TH2F *)calPad2.MakeHisto2DSmPl(k,0,&calDet,0,0.2,2.0,-1))->Draw("colz");
	cpadstatusrms->cd(2);
	((TH2F *)calPad2.MakeHisto2DSmPl(k,1,&calDet,0,0.2,2.0,-1))->Draw("colz");
	cpadstatusrms->cd(3);
	((TH2F *)calPad2.MakeHisto2DSmPl(k,2,&calDet,0,0.2,2.0,-1))->Draw("colz");
	cpadstatusrms->cd(4);
	((TH2F *)calPad2.MakeHisto2DSmPl(k,3,&calDet,0,0.2,2.0,-1))->Draw("colz");
	cpadstatusrms->cd(5);
	((TH2F *)calPad2.MakeHisto2DSmPl(k,4,&calDet,0,0.2,2.0,-1))->Draw("colz");
	cpadstatusrms->cd(6);
	((TH2F *)calPad2.MakeHisto2DSmPl(k,5,&calDet,0,0.2,2.0,-1))->Draw("colz");
	nbsm++;
      }
    }
  }
  printf("there is %d with Padstatus reference entries\n",nbsm);
  if(nbsm==0) ErrorRefDataPadStatus = 1;


  //Test reference data vdriftt0 DAQ
  //***************************
  Int_t ErrorRefDataVdriftT0DAQ = 0;
  if(entry) delete entry;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainRefStorage())->Get("TRD/DAQData/VdriftT0", 8);
  if(!entry) ErrorRefDataVdriftT0DAQ = 1;
  else{
    TProfile2D *histovdriftt0daq = (TProfile2D *) entry->GetObject();
    if(!histovdriftt0daq) ErrorRefDataVdriftT0DAQ = 2;
    else{
      Int_t NbinsX = histovdriftt0daq->GetNbinsY();
      if(NbinsX != 540) ErrorRefDataVdriftT0DAQ = 3;
      TCanvas *cvdriftdaq = new TCanvas("cvdriftdaq","",50,50,600,800);
      cvdriftdaq->cd();
      histovdriftt0daq->Draw("LEGO");
    }
  }


  //Check the detector OCDB values
  //********************************
 
  //Test for pads
  //******************
  //Gain
  //*****
  AliTRDCalPad *calPad = 0x0;
  Int_t ErrorGainPad = 0;
  if(entry) delete entry;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())->Get("TRD/Calib/LocalGainFactor", 8);
  if(!entry) ErrorGainPad++;
  else{
    calPad = (AliTRDCalPad *) entry->GetObject();
    if(!calPad) ErrorGainPad++;
    else{
      for(Int_t det = 0; det < 540; det++){
	AliTRDCalROC *calROC = calPad->GetCalROC(det);
	for(Int_t channel =0; channel < calROC->GetNchannels(); channel++){
	  if(calROC->GetValue(channel) != 1.0) ErrorGainPad++;
	}//channel loop
      }//det loop
    }
  }

  //Vdrift
  //*****
  Int_t ErrorVdriftPad = 0;
  if(entry) delete entry;
  if(calPad) delete calPad;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())->Get("TRD/Calib/LocalVdrift", 8);
  if(!entry) ErrorVdriftPad++;
  else{
    calPad = (AliTRDCalPad *) entry->GetObject();
    if(!calPad) ErrorVdriftPad++;
    else{
      for(Int_t det = 0; det < 540; det++){
	AliTRDCalROC *calROC = calPad->GetCalROC(det);
	for(Int_t channel =0; channel < calROC->GetNchannels(); channel++){
	  if(calROC->GetValue(channel) != 1.0) ErrorVdriftPad++;
	}//channel loop
      }//det loop
    }
  }

  //T0
  //*****
  Int_t ErrorT0Pad = 0;
  if(entry) delete entry;
  if(calPad) delete calPad;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())->Get("TRD/Calib/LocalT0", 8);
  if(!entry) ErrorT0Pad++;
  else{
    calPad = (AliTRDCalPad *) entry->GetObject();
    if(!calPad) ErrorT0Pad++;
    else{
      for(Int_t det = 0; det < 540; det++){
	AliTRDCalROC *calROC = calPad->GetCalROC(det);
	for(Int_t channel =0; channel < calROC->GetNchannels(); channel++){
	  if(calROC->GetValue(channel) != 0.0) ErrorT0Pad++;
	}//channel loop
      }//det loop
    }
  }


  //PRFWidth
  //********
  Int_t ErrorPRFWidthPad = 0;
  if(entry) delete entry;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())->Get("TRD/Calib/PRFWidth", 8);
  if(!entry) ErrorPRFWidthPad++;
  else{
    AliTRDCalPad *calPadPrf = (AliTRDCalPad *) entry->GetObject();
    if(!calPadPrf) ErrorPRFWidthPad++;
    else{
      Float_t value = 0.0;
      
      for(Int_t plane = 0; plane < 6; plane++){
	
	if(plane == 0) value = 0.515;
	if(plane == 1) value = 0.502;
	if(plane == 2) value = 0.491;
	if(plane == 3) value = 0.481;
	if(plane == 4) value = 0.471;
	if(plane == 5) value = 0.463;
	
	for(Int_t chamber = 0; chamber < 5; chamber++){
	  for(Int_t sector = 0; sector < 18; sector++){
	    
	    AliTRDCalROC *calROC = calPadPrf->GetCalROC(plane,chamber,sector);
	    for(Int_t channel =0; channel < calROC->GetNchannels(); channel++){
	      if((calROC->GetValue(channel) > 1.25*value) || (calROC->GetValue(channel) < 0.8*value)) ErrorPRFWidthPad++;
	    }//channel loop
	  }//sector loop
	}//chamber loop
      }//plane loop
      TCanvas *cpadprf = new TCanvas("cpadprf","cpadprf",50,50,600,800);
      cpadprf->cd();
      ((TH1F *)calPadPrf->MakeHisto1D(0,0,0.45,0.59,-1))->Draw();
    }
  }


  
  //Padstatus
  //********
  Int_t ErrorPadStatusPad = 0;
    
  if(entry) delete entry;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())->Get("TRD/Calib/PadStatus", 8);
  if(!entry) ErrorPadStatusPad++;
  else{
    AliTRDCalPadStatus *calPadStatus = (AliTRDCalPadStatus *) entry->GetObject();
    if(!calPadStatus) ErrorPadStatusPad++;
    else{
      for(Int_t k = 0; k < 18; k++){
	
	TString padstatus("PadStatus");
	padstatus += k;

	TCanvas *cpadstatus = new TCanvas((const char*)padstatus,(const char*)padstatus,50,50,600,800);
	cpadstatus->Divide(3,2);
	cpadstatus->cd(1);
	((TH2F *)calPadStatus->MakeHisto2DSmPl(k,0))->Draw("colz");
	cpadstatus->cd(2);
	((TH2F *)calPadStatus->MakeHisto2DSmPl(k,1))->Draw("colz");
	cpadstatus->cd(3);
	((TH2F *)calPadStatus->MakeHisto2DSmPl(k,2))->Draw("colz");
	cpadstatus->cd(4);
	((TH2F *)calPadStatus->MakeHisto2DSmPl(k,3))->Draw("colz");
	cpadstatus->cd(5);
	((TH2F *)calPadStatus->MakeHisto2DSmPl(k,4))->Draw("colz");
	cpadstatus->cd(6);
	((TH2F *)calPadStatus->MakeHisto2DSmPl(k,5))->Draw("colz");
      }
    }
  }

  //Test for detector values
  //*************************

  //Gain
  //******
  Int_t ErrorGainDetector = 0;
  if(entry) delete entry;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())->Get("TRD/Calib/ChamberGainFactor", 8);  
  if(!entry) ErrorGainDetector++;
  else{
    AliTRDCalDet *objectg = (AliTRDCalDet *) entry->GetObject();
    if(!objectg) ErrorGainDetector++;
    else{
      for(Int_t det = 0; det < 540; det++){
	if((objectg->GetValue(det)> 2.0) || (objectg->GetValue(det) < 0.0)) ErrorGainDetector++;
      }
      TCanvas *cdetgain = new TCanvas("cdetgain","",50,50,600,800);
      cdetgain->cd();
      ((TH1F *)objectg->MakeHisto1Distribution(0.0,2.0,-1))->Draw();
    }
  }
 

  //Vdrift
  //******
  Int_t ErrorVdriftDetector = 0;
  if(entry) delete entry;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())->Get("TRD/Calib/ChamberVdrift", 8);  
  if(!entry) ErrorVdriftDetector++;
  else{
    AliTRDCalDet *objectv = (AliTRDCalDet *) entry->GetObject();
    if(!objectv) ErrorVdriftDetector++;
    else{
      for(Int_t det = 0; det < 540; det++){
	if((objectv->GetValue(det)> 2.5) || (objectv->GetValue(det) < 0.6)) ErrorVdriftDetector++;
      }
      TCanvas *cdetv = new TCanvas("cdetv","",50,50,600,800);
      cdetv->cd();
      ((TH1F *)objectv->MakeHisto1Distribution(0.6,2.5,-1))->Draw();
    }
  }


  //T0
  //******
  Int_t ErrorT0Detector = 0;
  if(entry) delete entry;
  entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())->Get("TRD/Calib/ChamberT0", 8);  
  if(!entry) ErrorT0Detector++;
  else{
    AliTRDCalDet *objectt = (AliTRDCalDet *) entry->GetObject();
    if(!objectt) ErrorT0Detector++;
    else{
      for(Int_t det = 0; det < 540; det++){
	if((objectt->GetValue(det)> 1.0) || (objectt->GetValue(det) < -1.0)) ErrorT0Detector++;
      }
      TCanvas *cdett = new TCanvas("cdett","",50,50,600,800);
      cdett->cd();
      ((TH1F *)objectt->MakeHisto1Distribution(-0.5,0.5,-1))->Draw();
    }
  }

  

  //Bilan
  //**************

  //OCDB values
  //************

  printf("For the local gain factor there are %d strange values\n",ErrorGainPad);
  printf("For the local vdrift there are %d strange values\n",ErrorVdriftPad);
  printf("For the local t0 there are %d strange values\n",ErrorT0Pad);
  printf("For the chamber gain factor there are %d strange values\n",ErrorGainDetector);
  printf("For the chamber vdrift there are %d strange values\n",ErrorVdriftDetector);
  printf("For the chamber t0 there are %d strange values\n",ErrorT0Detector);
  printf("For the prf width there are %d strange values\n",ErrorPRFWidthPad);
 
  if(ErrorPadStatusPad > 0) printf("there is no calPadStatus object\n");


 //Reference data
 //****************

  //gain HLT
  
  if(ErrorRefDataGainHLT == 1) printf("There is no reference data entry for the gain HLT!\n");
  if(ErrorRefDataGainHLT == 2) printf("There is no reference data histogram for the gain HLT!\n");
  if(ErrorRefDataGainHLT == 3) printf("The reference data histogram has not the good number of Xbins for the gain HLT!\n");
  // vdrift HLT

  if(ErrorRefDataVdriftT0HLT == 1) printf("There is no reference data entry for the vdriftt0 HLT!\n");
  if(ErrorRefDataVdriftT0HLT == 2) printf("There is no reference data histogram for the vdriftt0 HLT!\n");
  if(ErrorRefDataVdriftT0HLT == 3) printf("The reference data profile has not the good number of Xbins for the gain HLT!\n");
  
  // prf HLT
  if(ErrorRefDataPRFHLT == 1) printf("There is no reference data entry for the prf HLT!\n");
  if(ErrorRefDataPRFHLT == 2) printf("There is no reference data profile for the prf HLT!\n");
  if(ErrorRefDataPRFHLT == 3) printf("The reference data profile has not the good number of Xbins for the prf HLT!\n");
  
  // vdrift DAQ

  if(ErrorRefDataVdriftT0DAQ == 1) printf("There is no reference data entry for the vdriftt0 DAQ!\n");
  if(ErrorRefDataVdriftT0DAQ == 2) printf("There is no reference data histogram for the vdriftt0 DAQ!\n");
  if(ErrorRefDataVdriftT0DAQ == 3) printf("The reference data profile has not the good number of Xbins for the gain DAQ!\n");


  // PadStatus DAQ

  if(ErrorRefDataPadStatus == 1) printf("There is no reference data entry for the pad status DAQ!\n");
  if(ErrorRefDataPadStatus == 2) printf("There is no reference data object for pad status DAQ!\n");
  
 
}
