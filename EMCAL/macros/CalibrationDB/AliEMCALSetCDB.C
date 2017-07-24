///
/// \file AliEMCALSetCDB.C
/// \ingroup EMCAL_CalibDB
/// \brief Set the energy calibration parameters in OCDB
///
/// Script to create energy calibration parameters and store them into CDB
/// Thre sets of calibration parameters can be created:
/// 1) equal parameters
/// 2) randomly distributed parameters for decalibrated detector silumations
/// 3) gausian  distributed parameters for decalibrated detector silumations
/// 
/// Execute like this:
///.x $ALICE_ROOT/EMCAL/macros/CalibrationDB/AliEMCALSetCDB.C
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__)
#include <TControlBar.h>
#include <TString.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>

#include "AliRun.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALGeoParams.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif

///
/// Main method
/// When execution, menu appears
//------------------------------------------------------------------------
void AliEMCALSetCDB()
{
  TControlBar *menu = new TControlBar("vertical","EMCAL CDB");
  menu->AddButton("Help to run EMCAL CDB","Help()",
		  "Explains how to use EMCAL CDS menus");
  menu->AddButton("Equal CC","SetCC(0)",
		  "Set equal calibration coefficients");
  menu->AddButton("Random De-calibration","SetCC(1)",
		  "Set random decalibration calibration coefficients");
  menu->AddButton("Gaussian De-calibration","SetCC(2)",
		  "Set gausian decalibration calibration coefficients");
  menu->AddButton("Read equal CC","GetCC(0)",
		  "Read initial equal calibration coefficients");
  menu->AddButton("Read random CC","GetCC(1)",
		  "Read random decalibration calibration coefficients");
  menu->AddButton("Read gaussian CC","GetCC(2)",
		  "Read gausian decalibration calibration coefficients");
  menu->Show();
}

//------------------------------------------------------------------------
void Help()
{
  printf(
         "\nSet calibration parameters and write them into ALICE CDB. Press button \"Equal CC\" "
         " to create equal pedestals and gain factors. Press button \"Decalibrate\" to create random"
         "pedestals and gain factors to imitate decalibrated detector\n"
         );
}

///
/// Writing calibration coefficients into the Calibration DB
/// Arguments:
///   flag=0: all calibration coefficients are equal
///   flag=1: all calibration coefficients random (decalibration)
///   flag=2: all calibration coefficients have Gaussian random distribution (decalibration)
//------------------------------------------------------------------------
void SetCC(Int_t flag=0)
{
  // Writing calibration coefficients into the Calibration DB
  // Arguments:
  //   flag=0: all calibration coefficients are equal
  //   flag=1: all calibration coefficients random (decalibration)
  //   flag=2: all calibration coefficients have Gaussian random distribution (decalibration)
  
  TString DBFolder;
  Int_t firstRun   =  0; // What is this
  Int_t lastRun    = 10;
  Int_t beamPeriod =  1;
  TString objFormat  = "";
  
  if      (flag == 0)
  {
    DBFolder  ="local://InitCalibDB";
    firstRun  =  0;
    lastRun   =  AliCDBRunRange::Infinity();
    objFormat = "EMCAL initial gain factors and pedestals";
  }
  else if (flag == 1)
  {
    DBFolder  ="local://DeCalibDB";
    firstRun  =  0;
    lastRun   = AliCDBRunRange::Infinity();
    objFormat = "EMCAL random pedestals and ADC gain factors (12x48x24)";
  }
  else if (flag == 2)
  {
    DBFolder  ="local://DeCalibDB"; // create directory DeCalibDB in current directory
    firstRun  =  0;
    lastRun   = AliCDBRunRange::Infinity();
    objFormat = "EMCAL random pedestals and gausian ADC gain factors (12x48x24)";
  }
  
  AliEMCALCalibData *calibda=new AliEMCALCalibData("EMCAL");
  
  Float_t fADCpedestal = 0.000;
  Float_t fADCchannel  = 0.0162;  // 250 GeV / (16*1024)
  Float_t rDecalibration  = 0.02 * fADCchannel; // 2% decalibration - just a guess
  Float_t cc=0, ped;
  
  TRandom rn;
  Int_t nSMod = AliEMCALGeoParams::fgkEMCALModules;
  
  for(Int_t supermodule=0; supermodule < nSMod; supermodule++)
  {
    Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;
    Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;
    
    // Set all the channels even the known to not exist in 1/3 sm and DCAL
    for(Int_t column=0; column< nCol; column++)
    {
      for(Int_t row=0; row< nRow; row++)
      {
        cc  = fADCchannel;
        ped = fADCpedestal;
        if (flag == 1)
        {
          // Decalibration:
          // Spread calibration coefficients uniformly with
          // Cmax/Cmin = 5, (Cmax-Cmin)/2 = 0.0015
          // and pedestals 0.005 +-10%
          //	  fADCchannel  = rn.Uniform(0.00075,0.00375);
          cc  = rn.Uniform(0.00140,0.00160);
          ped = 0;
        }
        else if (flag == 2)
        { // Gaussian
          cc  = rn.Gaus(fADCchannel, rDecalibration);
          ped = rn.Uniform(0.0045,0.0055);
        }
        
        calibda->SetADCchannel      (supermodule,column,row, cc);
        calibda->SetADCchannelOnline(supermodule,column,row, cc);
        calibda->SetADCpedestal     (supermodule,column,row, ped);
        cout<<"Set SM: "<<supermodule<<" col "<<column<<" row "<<row
        <<" cc "<< cc <<" ped "<<ped<<endl;
      }
    }
  }
  
  //Store calibration data into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Gustavo Conesa");
  
  AliCDBId id("EMCAL/Calib/Data",firstRun,lastRun); // create in EMCAL/Calib/Data DBFolder
  
  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBStorage* loc = man->GetStorage(DBFolder.Data());
  loc->Put(calibda, id, &md);
}

///
/// Read calibration coefficients into the Calibration DB
/// Arguments:
///   flag=0: all calibration coefficients are equal
///   flag=1: all calibration coefficients random (decalibration)
//------------------------------------------------------------------------
void GetCC(Int_t flag=0)
{
  TString DBFolder;
  int drawKey=2;
  
  if      (flag == 0)
  {
    DBFolder  ="local://InitCalibDB";
  }
  else if (flag == 1 || flag == 2)
  {
    //    DBFolder  ="local://DeCalibDB"; // Get DB in current folder
    DBFolder  ="local:///data/r22b/ALICE/PROD/CALIBRATION_May_2007/PI0/DeCalibDB"; // absolute name - Aug 31, 2007
    TString HOST(gSystem->Getenv("HOST"));
    if(HOST.Contains("pc")) { // pdsf; May 31, 2007
      DBFolder  ="local:///eliza5/alice/pavlinov/PROD/CALIBRATION_May_2007/PI0/10GEV/DECALIB/DeCalibDB";
    }
  }
  
  const Int_t runNumber = 0;
  AliEMCALCalibData* clb  = (AliEMCALCalibData*)
  (AliCDBManager::Instance()
   ->GetStorage(DBFolder.Data())
   ->Get("EMCAL/Calib/Data",runNumber)->GetObject());
  
  const Int_t nSMod = AliEMCALGeoParams::fgkEMCALModules;
  Int_t nCC   = 0;
  
  TH2F *hPed[nSMod], *hGain[nSMod], *hGainOnline[nSMod];
  TH1F *hCCSum    = new TH1F("hCCSum"," CC summary (in MeV) ", 200, 0.0, 20.);
  TH1F *hCCOnSum  = new TH1F("hCCSumOn"," CC online summary (in MeV) ", 200, 0.0, 20.);
  TH1F *hPedSum   = new TH1F("hPedSum"," pedestal summary (in MeV) ", 100, 4., 6.);
  
  TCanvas *cPed =0, *cGain =0, *cGainOn =0;
  TCanvas *cPed2=0, *cGain2=0, *cGainOn2=0;
  TCanvas *cPed3=0, *cGain3=0, *cGainOn3=0;
  TCanvas *cPed4=0, *cGain4=0, *cGainOn4=0;
  
  if(drawKey>1)
  {
    cPed     = new TCanvas("cPed" ,"Pedestals Mod 0-5"   , 10,10,400,800);
    cGain    = new TCanvas("cGain","Gain factors Mod 0-5", 410,10,400,800);
    cGainOn  = new TCanvas("cGainOn","Gain online factors Mod 0-5", 410,10,400,800);
    cPed    ->Divide(2,3);
    cGain   ->Divide(2,3);
    cGainOn ->Divide(2,3);
    
    cPed2    = new TCanvas("cPed2","Pedestals SMod 6-11", 10,10,400,800);
    cGain2   = new TCanvas("cGain2","Gain factors SMod 6-11", 410,10,400,800);
    cGainOn2 = new TCanvas("cGainOn2","Gain online factors SMod 6-11", 410,10,400,800);
    cPed2   ->Divide(2,3);
    cGain2  ->Divide(2,3);
    cGainOn2->Divide(2,3);
    
    cPed3    = new TCanvas("cPed3","Pedestals SMod 12-17", 10,10,400,800);
    cGain3   = new TCanvas("cGain3","Gain factors SMod 12-17", 410,10,400,800);
    cGainOn3 = new TCanvas("cGainOn3","Gain online factors SMod 7-17", 410,10,400,800);
    cPed3   ->Divide(2,3);
    cGain3  ->Divide(2,3);
    cGainOn3->Divide(2,3);
    
    cPed4    = new TCanvas("cPed4","Pedestals SMod 18-21", 10,10,400,400);
    cGain4   = new TCanvas("cGain4","Gain factors SMod 18-21", 410,10,400,800);
    cGainOn4 = new TCanvas("cGainOn4","Gain online factors SMod 18-21", 410,10,400,400);
    cPed4   ->Divide(2,2);
    cGain4  ->Divide(2,2);
    cGainOn4->Divide(2,2);
  }
  
  TCanvas *cSum   = new TCanvas("cSum" ,"summary"   , 10,10,600,800);
  cSum->Divide(1,2);
  
  cout<<endl;
  
  for (Int_t supermodule=0; supermodule<nSMod; supermodule++)
  {
    Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;
    Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;
    
//    if(supermodule /2 == 5)
//    nRow = nRow/2;
//    if(supermodule > 11)
//    nCol  = nCol*2/3;
    
    TString namePed="hPed";
    namePed+=supermodule;
    TString titlePed="Pedestals in supermodule ";
    titlePed+=supermodule;
    hPed[supermodule] = new TH2F(namePed.Data(),titlePed.Data(),
                                 nCol,1.,1.*nCol,nRow,1.,1.*nRow);
    
    TString nameGain="hGain";
    TString nameGainOnline="hGainOnline";
    nameGain+=supermodule;
    nameGainOnline+=supermodule;
    
    TString titleGain="Gain factors in supermodule ";
    titleGain+=supermodule;
    
    TString titleGainOnline="Gain online factors in supermodule ";
    titleGainOnline+=supermodule;
    
    hGain      [supermodule] = new TH2F(nameGain.Data(),titleGain.Data(),
                                        nCol,1.,1.*nCol,nRow,1.,1.*nRow);
    hGainOnline[supermodule] = new TH2F(nameGainOnline.Data(),titleGainOnline.Data(),
                                        nCol,1.,1.*nCol,nRow,1.,1.*nRow);
    
    for (Int_t column=0; column<nCol; column++)
    {
      for (Int_t row=0; row<nRow; row++)
      {
        Float_t ped    = clb->GetADCpedestal     (supermodule,column,row);
        Float_t gainOn = clb->GetADCchannelOnline(supermodule,column,row);
        Float_t gain   = clb->GetADCchannel      (supermodule,column,row);
        //cout<<"Get SM: "<<supermodule<<" col "<<column<<" row "<<row
        //<<" chan "<<gain<<endl;
        hPed       [supermodule]->SetBinContent(column+1,row+1,ped*1.e+3);  // in mev
        hGain      [supermodule]->SetBinContent(column+1,row+1,gain*1.e+3); // in mev
        hGainOnline[supermodule]->SetBinContent(column+1,row+1,gain*1.e+3); // in mev
        
        hPedSum ->Fill(ped   *1.e+3);
        hCCSum  ->Fill(gain  *1.e+3);
        hCCOnSum->Fill(gainOn*1.e+3);
        
        nCC++;
      }
    }
    
    cout<<" Fill cc for SM "<< supermodule << " nCC "<< nCC << endl;
    
    if(drawKey>1)
    {
      if(supermodule < 6)
      {
        cPed ->cd(supermodule+1);
        hPed[supermodule]->Draw("lego2");
        cGain->cd(supermodule+1);
        hGain[supermodule]->Draw("lego2");
        cGainOn->cd(supermodule+1);
        hGainOnline[supermodule]->Draw("lego2");
      }
      else if(supermodule < 12)
      {
        cPed2 ->cd(supermodule-5);
        hPed[supermodule]->Draw("lego2");
        cGain2->cd(supermodule-5);
        hGain[supermodule]->Draw("lego2");
        cGainOn2->cd(supermodule-5);
        hGainOnline[supermodule]->Draw("lego2");
      }
      else if(supermodule < 18)
      {
        cPed3 ->cd(supermodule-11);
        hPed[supermodule]->Draw("lego2");
        cGain3->cd(supermodule-11);
        hGain[supermodule]->Draw("lego2");
        cGainOn3->cd(supermodule-11);
        hGainOnline[supermodule]->Draw("lego2");
      }
      else
      {
        cPed4 ->cd(supermodule-17);
        hPed[supermodule]->Draw("lego2");
        cGain4->cd(supermodule-17);
        hGain[supermodule]->Draw("lego2");
        cGainOn4->cd(supermodule-17);
        hGainOnline[supermodule]->Draw("lego2");
      }
    }
  }
  cout << " Get "<<nCC<<" calibration coeffs"<<endl;
  
  cSum->cd(1);
  gStyle->SetOptFit(111);
  hCCSum->Fit("gaus");
  hCCSum->SetLineWidth(2);
  hCCSum->GetFunction("gaus")->SetLineColor(2);
  
  cSum->cd(2);
  hPedSum->Draw();
  hPedSum->SetLineWidth(2);
  
  cSum->Update();
  
  if(drawKey>1)
  {
   cPed   ->Print("pedestals_SM_0_6.eps");
   cGain  ->Print("gains_SM_0_5.eps");
   cPed2  ->Print("pedestals_SM_6_11.eps");
   cGain2 ->Print("gains_SM_6_11.eps");
  }
}
