/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for EMCAL cosmic calibration results                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TFile.h>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "TTree.h"
#include <Riostream.h>

#include "AliEMCALSMCalibCosmicResult.h"

ClassImp(AliEMCALSMCalibCosmicResult)


AliEMCALSMCalibCosmicResult::AliEMCALSMCalibCosmicResult()
{
  // Default Constructor  
}


AliEMCALSMCalibCosmicResult::AliEMCALSMCalibCosmicResult(const TString &name)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  cout << " object name: "<< name << endl;
  
}

AliEMCALSMCalibCosmicResult::~AliEMCALSMCalibCosmicResult()
{
  // Destructor
}


void AliEMCALSMCalibCosmicResult::InitArrays()
{
  //initialize arrays
  
  cout << "======== arrays are initialized ========"<< endl;

  for (int j=0;j<48;j++) {
      fLEDRefADC[j]=-1; 
      for (int i=0;i<24;i++) {
	  fMIPPeakADC[j][i] = -1;
	  fLEDPeakADC[j][i] = -1;
	  fAPDVoltage[j][i] = -1;
   	}
  }
  for (int j=0;j<3;j++) {
      fMinTemp[j] =-1; 
      fMaxTemp[j] =-1; 
      fMeanTemp[j]=-1; 
  }

}

void AliEMCALSMCalibCosmicResult::ReadSMPart(int ipart)
{
  // open file and get tree and histograms

  TString fin="./plots/results_part_";
  fin += ipart;
  fin  +=".root";
  
  cout << "ReadSMPart:: input file: "<< fin << endl;
  TFile *f1 = new TFile(fin,"READ");

  TH1D* hev = (TH1D*)f1->Get("ev");
  cout<<"Total number of events: "<< hev->GetEntries() << endl;

  fhref =(TH2D*) f1->Get("hLedRefStrip");
  fTree =(TTree*) f1->Get("calib_tree");
 
  fTree->SetBranchAddress("col",&fCol);
  fTree->SetBranchAddress("row",&fRow);
  fTree->SetBranchAddress("mip",&fMip);
  fTree->SetBranchAddress("led",&fLed);

}



void  AliEMCALSMCalibCosmicResult::PrintAPD() 
{
  // print out APD bias values
  
  for (Int_t col=0; col< fgkEmCalCols; col++){
    for (Int_t row=0; row< fgkEmCalRows; row++){
      printf(" %3.0f",fAPDVoltage[col][row]);
    }
    printf("\n");
  }
}

void  AliEMCALSMCalibCosmicResult::PrintMIP() 
{
  // print out mean MIP values

  for (Int_t col=0; col< fgkEmCalCols; col++){
    for (Int_t row=0; row< fgkEmCalRows; row++){
      printf(" %2.1f",fMIPPeakADC[col][row]);
    }
    printf("\n");
  }
}

void  AliEMCALSMCalibCosmicResult::PrintLED() 
{
  //print out mean  LED values
  
  for (Int_t col=0; col< fgkEmCalCols; col++){
    for (Int_t row=0; row< fgkEmCalRows; row++){
      printf(" %3.0f",fLEDPeakADC[col][row]);
    }
    printf("\n");
  }
}

void  AliEMCALSMCalibCosmicResult::PrintLEDref() 
{
  //print out LED reference

  for (Int_t col=0; col< fgkEmCalCols; col++){
      printf(" %3.0f",fLEDRefADC[col]);
    }
  printf("\n");
    
}

void  AliEMCALSMCalibCosmicResult::PrintTemps() 
{
  //print out temperature
  
  for (Int_t part=0; part< 3; part++){
      printf("min: %2.1f ,",fMinTemp[part]);
      printf("max: %2.1f ,",fMaxTemp[part]);
      printf("ave: %2.1f \n",fMeanTemp[part]);
    }

    
}



void  AliEMCALSMCalibCosmicResult::ReadLEDRefADCValues(int ipart)
{
  //read LED reference from histogram

  cout << "-------- read LED ref values SM part #: "<< ipart << " ---------"<< endl;
  
  int firstStrip = 0;
  int lastStrip  = 0;
  if (ipart==0) {
  	firstStrip = 0;
  	lastStrip  = 8;
  }  
  if (ipart==1) {
  	firstStrip = 8;
  	lastStrip  = 16;
  }  
  if (ipart==2) {
  	firstStrip = 16;
  	lastStrip  = 24;
  }  
  
  char title[30];
  float stripmean[8];
  
  for (int j=firstStrip;j< lastStrip;j++) {
	sprintf(title,"strip%d",j);
  	TH1D* hLedAmp = (TH1D*)fhref->ProjectionY(title,j+1,j+1);
	stripmean[j-firstStrip] = hLedAmp->GetMean();
	cout << title << ", mean-Ref="<< hLedAmp->GetMean() << endl;
        //2 cols per strip
	fLEDRefADC[2*j]   = hLedAmp->GetMean();
	fLEDRefADC[2*j+1] = hLedAmp->GetMean();	
  }

 
  return;
 
}


int  AliEMCALSMCalibCosmicResult::ReadMIPPeakADCValues(int th)
{
  // read MIP values from Tree

  cout << "--------- read MIP peak , threshold "<< th << " ---------"<< endl;
  
  int badMip=0;

  for (int j=0;j<fTree->GetEntries();j++)
    {
	fTree->GetEntry(j);
	if (fMip> th) 	fMIPPeakADC[fCol][fRow] = fMip;
	else  {        	fMIPPeakADC[fCol][fRow] = -1;
 			cout << "col:"<< fCol<< ", row:"<< fRow << ", mean =" << fMip << endl;  
			badMip++; 
	}
    }

  return badMip;

}


int  AliEMCALSMCalibCosmicResult::ReadLEDPeakADCValues(int th)
{
  // read LED values from Tree
  cout << "--------- read LED peaks, threshold "<< th << " ---------"<< endl;
  
  int bad=0;

  for (int j=0;j<fTree->GetEntries();j++)
    {
	fTree->GetEntry(j);
	if (fLed > th) 	fLEDPeakADC[fCol][fRow] = fLed;
	else  {        	fLEDPeakADC[fCol][fRow] = -1;
 			cout << "col:"<< fCol<< ", row:"<< fRow << ", mean =" << fLed << endl; 
			bad++;  
	}
    }

  return bad;


}

int  AliEMCALSMCalibCosmicResult::ReadAPDVoltageValues(int ipart, const TString &fname)
{
  // read APD voltages from text file
  cout << "--------- read APD Bias for sm-part "<< ipart << " ---------"<< endl;
  cout << "--------- filename: "<< fname << " ---------"<< endl;

  FILE *f=fopen(fname.Data(),"r");
  
  int offset=0;
  if (ipart==1) offset=16;
  if (ipart==2) offset=32;

  int dat, row,col;
  float voltage; 
  int bad=0;

  while ( !feof(f) ) 
     { 
       dat = fscanf(f,"%i %i %f",&col,&row,&voltage);

       if (col< offset || col>(offset+15) ) continue;    

       // cout<<col<<" "<<row<<" v="<<voltage<<endl;
       
       fAPDVoltage[col][row] = voltage;
       if (voltage>395) bad++;
       
     }
     
  return bad;
     
}

void AliEMCALSMCalibCosmicResult::ReadTempSensors(int ipart, int run)
{
  // read temperature sensor data 
  cout << "--------- read Temp Sensors for sm-part "<< ipart << " ---------"<< endl;
  cout << "--------- and run number "<< run << " ---------"<< endl;
  
  // take average of left and right sensors.
  
  //correction= measured - reference thermometer
  float tcorr[8] = {-0.4,0.7,0,0.9,-0.2,0,0.1,0.5};

  TString fin="./Temperature/temp_sensors_run_";
  fin += run;
  fin +=".root";
  
  TFile *f = new TFile(fin,"READ");
  
  float max1=0,min1=0,mean1=0;
  float max2=0,min2=0,mean2=0;
  
  if (ipart==0) {
    	TProfile* h1 = (TProfile*)f->Get("sensor1");
   	TProfile* h2 = (TProfile*)f->Get("sensor2");
    	TF1* fit1    = (TF1*)f->Get("fitsensor1");
   	TF1* fit2    = (TF1*)f->Get("fitsensor2");

	max1  = h1->GetMaximum()      - tcorr[0];
	min1  = h1->GetMinimum()      - tcorr[0];
	mean1 = fit1->GetParameter(0) - tcorr[0];
	max2  = h2->GetMaximum()      - tcorr[1];
	min2  = h2->GetMinimum()      - tcorr[1];
	mean2 = fit2->GetParameter(0) - tcorr[1];
  }
  
   if (ipart==1) {
    	TProfile* h1 = (TProfile*)f->Get("sensor3");
   	TProfile* h2 = (TProfile*)f->Get("sensor4");
    	TF1* fit1    = (TF1*)f->Get("fitsensor3");
   	TF1* fit2    = (TF1*)f->Get("fitsensor4");

	max1  = h1->GetMaximum()      - tcorr[2];
	min1  = h1->GetMinimum()      - tcorr[2];
	mean1 = fit1->GetParameter(0) - tcorr[2];
	max2  = h2->GetMaximum()      - tcorr[3];
	min2  = h2->GetMinimum()      - tcorr[3];
	mean2 = fit2->GetParameter(0) - tcorr[3];
  }
   
  if (ipart==2) {
    	TProfile* h1 = (TProfile*)f->Get("sensor7");
   	TProfile* h2 = (TProfile*)f->Get("sensor8");
    	TF1* fit1    = (TF1*)f->Get("fitsensor7");
   	TF1* fit2    = (TF1*)f->Get("fitsensor8");

	max1  = h1->GetMaximum()      - tcorr[6];
	min1  = h1->GetMinimum()      - tcorr[6];
	mean1 = fit1->GetParameter(0) - tcorr[6];
	max2  = h2->GetMaximum()      - tcorr[7];
	min2  = h2->GetMinimum()      - tcorr[7];
	mean2 = fit2->GetParameter(0) - tcorr[7];
  }
      
  fMinTemp[ipart]  = (min1+min2)/2.0;
  fMaxTemp[ipart]  = (max1+max2)/2.0;
  fMeanTemp[ipart] = (mean1+mean2)/2.0;
    
  return;
 
 
}

