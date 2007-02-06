
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

/* $Id$ */

//_________________________________________________________________________
// An analysis task to check the PMD data in simulated data
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h> 
#include <TH1F.h>
#include <TH2F.h>
#include <TLine.h> 
#include <TROOT.h>
#include <TStyle.h> 

#include "AliPMDQATask.h" 
#include "AliPMDUtility.h" 
#include "AliESD.h" 
#include "AliLog.h"

//______________________________________________________________________________
AliPMDQATask::AliPMDQATask(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0), 
  fhPMDP1(0),
  fhPMDC2(0),
  fhPMDP2(0),
  fhPMDC3(0),
  fhPMDP3(0),
  fhPMDP4(0),
  fhPMDC5(0),
  fhPMDP5(0),
  fhPMDCP0(0),
  fhPMDCP1(0),
  fhPMDCP2(0),
  fhPMDCP3(0),
  fhPMDCP4(0),
  fhPMDSM1(0),
  fhPMDSM2(0),
  fhPMDSM3(0),
  fhPMDSM4(0),
  fhPMDSM5(0),
  fhPMDSM6(0),
  fhPMDSM7(0),
  fhPMDSM8(0),
  fhPMDSM9(0),
  fhPMDSM10(0),
  fhPMDSM11(0),
  fhPMDSM12(0),
  fhPMDSM13(0),
  fhPMDSM14(0),
  fhPMDSM15(0),
  fhPMDSM16(0),
  fhPMDSM17(0),
  fhPMDSM18(0),
  fhPMDSM19(0),
  fhPMDSM20(0),
  fhPMDSM21(0),
  fhPMDSM22(0),
  fhPMDSM23(0),
  fhPMDSM24(0),
  fhPMDSM (0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

//______________________________________________________________________________
AliPMDQATask::~AliPMDQATask()
{
  // dtor
  fOutputContainer->Clear() ; 
  delete fOutputContainer ;
  
  delete fhPMDP1  ;
  delete fhPMDC2  ;
  delete fhPMDP2  ;
  delete fhPMDC3  ;
  delete fhPMDP3  ;
  delete fhPMDP4  ;
  delete fhPMDC5  ;
  delete fhPMDP5  ;
  delete fhPMDCP0 ;
  delete fhPMDCP1 ;
  delete fhPMDCP2 ;
  delete fhPMDCP3 ;
  delete fhPMDCP4 ;
  delete fhPMDSM1  ;
  delete fhPMDSM2  ;
  delete fhPMDSM3  ;
  delete fhPMDSM4  ;
  delete fhPMDSM5  ;
  delete fhPMDSM6  ;
  delete fhPMDSM7  ;
  delete fhPMDSM8  ;
  delete fhPMDSM9  ;
  delete fhPMDSM10 ;
  delete fhPMDSM11 ;
  delete fhPMDSM12 ;
  delete fhPMDSM13 ;
  delete fhPMDSM14 ;
  delete fhPMDSM15 ;
  delete fhPMDSM16 ;
  delete fhPMDSM17 ;
  delete fhPMDSM18 ;
  delete fhPMDSM19 ;
  delete fhPMDSM20 ;
  delete fhPMDSM21 ;
  delete fhPMDSM22 ;
  delete fhPMDSM23 ;
  delete fhPMDSM24 ;
  delete fhPMDSM   ;
  
}

//______________________________________________________________________________
void AliPMDQATask::ConnectInputData(const Option_t*)
{
  // Initialisation of branch container and histograms 
    
  AliInfo(Form("*** Initialization of %s", GetName())) ; 
    
  // Get input data
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {
    AliError(Form("Input 0 for %s not found\n", GetName()));
    return ;
  }
  
  // One should first check if the branch address was taken by some other task
  char ** address = (char **)GetBranchAddress(0, "ESD");
  if (address) {
    fESD = (AliESD*)(*address);
  } else {
    fESD = new AliESD();
    SetBranchAddress(0, "ESD", &fESD);
  }
}

//________________________________________________________________________
void AliPMDQATask::CreateOutputObjects()
{  
  // create histograms 
  
  fhPMDP1   = new TH2F("fhPMDP1","XY of Clusters",100,-100.,100.,100,-100.,100.);
  fhPMDC2   = new TH1F("fhPMDC2","CPV  PHI",200,-1,9);
  fhPMDP2   = new TH1F("fhPMDP2","PRE  PHI",200,-1,9);
  fhPMDC3   = new TH1F("fhPMDC3","CPV  Clus",30,0.,500.);
  fhPMDP3   = new TH1F("fhPMDP3","PRE  N-gammalike",20,0.,500.);
  fhPMDP4   = new TH1F("fhPMDP4","PRE  EDEP",30,0.,1000.);
  fhPMDC5   = new TH1F("fhPMDC5","CPV  n-cell",20,0.,100.);
  fhPMDP5   = new TH1F("fhPMDP5","PMD  n-cell",20,0.,100.);
  fhPMDCP0  = new TH2F("fhPMDCP0","PRE CLUS Quad.1 vs 2",150,0.,300.,150,0.,300.);
  fhPMDCP1  = new TH2F("fhPMDCP1","PRE CLUS Quad.3 vs 4",150,0.,300.,150,0.,300.);
  fhPMDCP2  = new TH2F("fhPMDCP2","PRE EDEP Quad.3 vs 4",50,0.,300.,50,0.,300.);
  fhPMDCP3  = new TH2F("fhPMDCP3","PRE EDEP vs Tot Clus ",10,0.,1000.,10,0.,300.);
  fhPMDCP4  = new TH2F("fhPMDCP4","PRE Clus vs CPV Clus ",150,0.,200.,150,0.,200.);

  fhPMDSM1  = new TH2F("fhPMDSM1","PRE Cluster XY",200,-100,100,200,-100,100);
  fhPMDSM2  = new TH2F("fhPMDSM2","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM3  = new TH2F("fhPMDSM3","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM4  = new TH2F("fhPMDSM4","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM5  = new TH2F("fhPMDSM5","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM6  = new TH2F("fhPMDSM6","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM7  = new TH2F("fhPMDSM7","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM8  = new TH2F("fhPMDSM8","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM9  = new TH2F("fhPMDSM9","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM10 = new TH2F("fhPMDSM10","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM11 = new TH2F("fhPMDSM11","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM12 = new TH2F("fhPMDSM12","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM13 = new TH2F("fhPMDSM13","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM14 = new TH2F("fhPMDSM14","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM15 = new TH2F("fhPMDSM15","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM16 = new TH2F("fhPMDSM16","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM17 = new TH2F("fhPMDSM17","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM18 = new TH2F("fhPMDSM18","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM19 = new TH2F("fhPMDSM19","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM20 = new TH2F("fhPMDSM20","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM21 = new TH2F("fhPMDSM21","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM22 = new TH2F("fhPMDSM22","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM23 = new TH2F("fhPMDSM23","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM24 = new TH2F("fhPMDSM24","",999,-100.0,100.0,999,-100.0,100.0);
  fhPMDSM   = new TH1F("fhPMDSM","Plot of all 24 Super Modules",24,0,24);
  
  // create output container
  
  fOutputContainer = new TObjArray(38) ; 
  fOutputContainer->SetName("PMD") ; 

  fOutputContainer->AddAt(fhPMDP1,    0 );
  fOutputContainer->AddAt(fhPMDC2,    1 );
  fOutputContainer->AddAt(fhPMDP2,    2 );
  fOutputContainer->AddAt(fhPMDC3,    3 );
  fOutputContainer->AddAt(fhPMDP3,    4 );
  fOutputContainer->AddAt(fhPMDP4,    5 );
  fOutputContainer->AddAt(fhPMDC5,    6 );
  fOutputContainer->AddAt(fhPMDP5,    7 );
  fOutputContainer->AddAt(fhPMDCP0,   8 );
  fOutputContainer->AddAt(fhPMDCP1,   9);
  fOutputContainer->AddAt(fhPMDCP2,  10 );
  fOutputContainer->AddAt(fhPMDCP3,  11 );
  fOutputContainer->AddAt(fhPMDCP4,  12 );
  
  fOutputContainer->AddAt(fhPMDSM1,  13 );
  fOutputContainer->AddAt(fhPMDSM2,  14 );
  fOutputContainer->AddAt(fhPMDSM3,  15 );
  fOutputContainer->AddAt(fhPMDSM4,  16 );
  fOutputContainer->AddAt(fhPMDSM5,  17 );
  fOutputContainer->AddAt(fhPMDSM6,  18 );
  fOutputContainer->AddAt(fhPMDSM7,  19 );
  fOutputContainer->AddAt(fhPMDSM8,  20 );
  fOutputContainer->AddAt(fhPMDSM9,  21 );
  fOutputContainer->AddAt(fhPMDSM10, 22 );
  fOutputContainer->AddAt(fhPMDSM11, 23 );
  fOutputContainer->AddAt(fhPMDSM12, 24 );
  fOutputContainer->AddAt(fhPMDSM13, 25 );
  fOutputContainer->AddAt(fhPMDSM14, 26 );
  fOutputContainer->AddAt(fhPMDSM15, 27 );
  fOutputContainer->AddAt(fhPMDSM16, 28 );
  fOutputContainer->AddAt(fhPMDSM17, 29 );
  fOutputContainer->AddAt(fhPMDSM18, 30 );
  fOutputContainer->AddAt(fhPMDSM19, 31 );
  fOutputContainer->AddAt(fhPMDSM20, 32 );
  fOutputContainer->AddAt(fhPMDSM21, 33 );
  fOutputContainer->AddAt(fhPMDSM22, 34 );
  fOutputContainer->AddAt(fhPMDSM23, 35 );
  fOutputContainer->AddAt(fhPMDSM24, 36 );
  fOutputContainer->AddAt(fhPMDSM,   37 );
}

//______________________________________________________________________________
void AliPMDQATask::Exec(Option_t *) 
{
  // Processing of one event
  
  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
  
  // ************************ PMD

  AliPMDUtility *cc = new AliPMDUtility(); 
  
  Int_t smn;
  Int_t n=0;
  Float_t totCPVClus     ;
  Float_t totPREClus     ;
  Float_t totPREEdep  ;
  Float_t totCPVCell     ;
  Float_t totPRECell     ;
  Float_t preCluQUAD[4]  ;
  Float_t cpvCluQUAD[4]  ;
  Float_t preADCQUAD[4]  ;
  Float_t cpvADCQUAD[4]  ;
  Float_t preCelQUAD[4]  ;
  Float_t cpvCelQUAD[4]  ;

  Int_t npmdCl = fESD->GetNumberOfPmdTracks();
  
  // ****** The loop over PMD clusters

  for (Int_t kk = 0; kk < 4 ; kk++) {
    cpvCluQUAD[kk] = 0.0 ;
    preCluQUAD[kk] = 0.0 ;
    cpvCelQUAD[kk] = 0.0 ;
    preCelQUAD[kk] = 0.0 ;
    preADCQUAD[kk] = 0.0 ;
  } 
 
  while (npmdCl--) {
    
    AliESDPmdTrack * pmdtr = fESD->GetPmdTrack(npmdCl);
    Int_t   det   = pmdtr->GetDetector();
    Float_t clsX  = pmdtr->GetClusterX();
    Float_t clsY  = pmdtr->GetClusterY();
    Float_t clsZ  = pmdtr->GetClusterZ();
    Float_t ncell = pmdtr->GetClusterCells();
    Float_t adc   = pmdtr->GetClusterADC();
    
    cc->SetXYZ(clsX,clsY,clsZ);
    cc->CalculateEta();
    cc->CalculatePhi();
    Float_t eta = cc->GetEta();
    Float_t phi = cc->GetPhi();
    
   // Calculating S.Module Number from Cluster .
    
    CalculateSMN(clsX, clsY, smn);
    if( det == 1)
      {
	if(smn >= 0 && smn <= 5) {
	  ++cpvCluQUAD[0] ;
	  cpvADCQUAD[0] =+ adc  ;
	  cpvCelQUAD[0] =+ ncell ;
	}
	if(smn >= 6 && smn <=11) {
	  ++cpvCluQUAD[1] ;
	  cpvADCQUAD[1] =+ adc  ;
	  cpvCelQUAD[1] =+ ncell ;
	}
	if(smn >=12 && smn <=17) {
	  ++cpvCluQUAD[2] ;
	  cpvADCQUAD[2] =+ adc  ;
	  cpvCelQUAD[2] =+ ncell ;
	}
	if(smn >=18 && smn <=23) {
	  ++cpvCluQUAD[3] ;
	  cpvADCQUAD[3] =+ adc  ;
	  cpvCelQUAD[3] =+ ncell ;
	}
	
	if(eta >= 2.3 && eta <= 3.5)
	  {
	    fhPMDC2->Fill(phi);
	  }
      }
    if( det == 0)
      {
	if(smn >= 0 && smn <= 5) { 
	  ++preCluQUAD[0] ;
	  preADCQUAD[0] =+ adc  ;    
	  preCelQUAD[0] =+ ncell ;    
	}    
	if(smn >= 6 && smn <=11) { 
	  ++preCluQUAD[1] ;
	  preADCQUAD[1] =+ adc  ;    
	  preCelQUAD[1] =+ ncell ;    
	}    
	if(smn >=12 && smn <=17) { 
	  ++preCluQUAD[2] ;
	  preADCQUAD[2] =+ adc  ;    
	  preCelQUAD[2] =+ ncell ;    
	}    
	if(smn >=18 && smn <=23) { 
	  ++preCluQUAD[3] ;
	  preADCQUAD[3] =+ adc  ;    
	  preCelQUAD[3] =+ ncell ;    
	}    
	if ( n <= 100 ) { 
	  fhPMDSM->Fill(smn);
	  if(smn == 0) fhPMDSM1->Fill(-clsX,clsY);
	  if(smn == 0) fhPMDSM1->Fill(-clsX,clsY);
	  if(smn == 1) fhPMDSM2->Fill(-clsX,clsY);
	  if(smn == 2) fhPMDSM3->Fill(-clsX,clsY);
	  if(smn == 3) fhPMDSM4->Fill(-clsX,clsY);
	  if(smn == 4) fhPMDSM5->Fill(-clsX,clsY);
	  if(smn == 5) fhPMDSM6->Fill(-clsX,clsY);
	  if(smn == 6) fhPMDSM7->Fill(-clsX,clsY);
	  if(smn == 7) fhPMDSM8->Fill(-clsX,clsY);
	  if(smn == 8) fhPMDSM9->Fill(-clsX,clsY);
	  if(smn == 9) fhPMDSM10->Fill(-clsX,clsY);
	  if(smn ==10) fhPMDSM11->Fill(-clsX,clsY);
	  if(smn ==11) fhPMDSM12->Fill(-clsX,clsY);
	  if(smn ==12) fhPMDSM13->Fill(-clsX,clsY);
	  if(smn ==13) fhPMDSM14->Fill(-clsX,clsY);
	  if(smn ==14) fhPMDSM15->Fill(-clsX,clsY);
	  if(smn ==15) fhPMDSM16->Fill(-clsX,clsY);
	  if(smn ==16) fhPMDSM17->Fill(-clsX,clsY);
	  if(smn ==17) fhPMDSM18->Fill(-clsX,clsY);
	  if(smn ==18) fhPMDSM19->Fill(-clsX,clsY);
	  if(smn ==19) fhPMDSM20->Fill(-clsX,clsY);
	  if(smn ==20) fhPMDSM21->Fill(-clsX,clsY);
	  if(smn ==21) fhPMDSM22->Fill(-clsX,clsY);
	  if(smn ==22) fhPMDSM23->Fill(-clsX,clsY);
	  if(smn ==23) fhPMDSM24->Fill(-clsX,clsY);
	}     
	if(eta >= 2.3 && eta <= 3.5)
	  {
	    fhPMDP2->Fill(phi);
	  }
	fhPMDP1->Fill(clsX,clsY);
      }
  } 
  for (Int_t k = 0 ; k < 4 ; k++) {
    totCPVClus =+ cpvCluQUAD [k] ;
    totPREClus =+ preCluQUAD [k] ;
    totCPVCell =+ cpvCelQUAD [k] ;
    totPRECell =+ preCelQUAD [k] ;
    totPREEdep =+ preADCQUAD [k] ;     
    }
  Float_t totCPVpreClus = totPREClus + totCPVClus ;

  //  if(eta >= 2.3 && eta <= 3.5) {
  fhPMDC3->Fill(totCPVClus);
  fhPMDP3->Fill(totPREClus);
  fhPMDP4->Fill(totPREEdep);
  fhPMDP5->Fill(totPRECell);
  fhPMDCP0->Fill(preCluQUAD[0],preCluQUAD[1]);
  fhPMDCP1->Fill(preCluQUAD[2],preCluQUAD[3]);
  fhPMDCP2->Fill(preADCQUAD[2],preADCQUAD[3]);
  fhPMDCP3->Fill(totPREEdep,totCPVpreClus);
  fhPMDCP4->Fill(totPREClus,totCPVClus);
  //    }
  totCPVClus = 0.0; 
  totPREClus = 0.0; 
  totCPVCell = 0.0; 
  totPRECell = 0.0; 
  totPREEdep = 0.0; 

  PostData(0, fOutputContainer);
} 

//______________________________________________________________________________
void AliPMDQATask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  fOutputContainer = (TObjArray*)GetOutputData(0);
  
  fhPMDP1   = (TH2F*)fOutputContainer->At(0);
  fhPMDC2   = (TH1F*)fOutputContainer->At(1);
  fhPMDP2   = (TH1F*)fOutputContainer->At(2);
  fhPMDC3   = (TH1F*)fOutputContainer->At(3);
  fhPMDP3   = (TH1F*)fOutputContainer->At(4);
  fhPMDP4   = (TH1F*)fOutputContainer->At(5);
  fhPMDC5   = (TH1F*)fOutputContainer->At(6);
  fhPMDP5   = (TH1F*)fOutputContainer->At(7);
  fhPMDCP0  = (TH2F*)fOutputContainer->At(8);
  fhPMDCP1  = (TH2F*)fOutputContainer->At(9);
  fhPMDCP2  = (TH2F*)fOutputContainer->At(10);
  fhPMDCP3  = (TH2F*)fOutputContainer->At(11);
  fhPMDCP4  = (TH2F*)fOutputContainer->At(12);

  fhPMDSM1  = (TH2F*)fOutputContainer->At(13);
  fhPMDSM2  = (TH2F*)fOutputContainer->At(14);
  fhPMDSM3  = (TH2F*)fOutputContainer->At(15);
  fhPMDSM4  = (TH2F*)fOutputContainer->At(16);
  fhPMDSM5  = (TH2F*)fOutputContainer->At(17);
  fhPMDSM6  = (TH2F*)fOutputContainer->At(18);
  fhPMDSM7  = (TH2F*)fOutputContainer->At(19);
  fhPMDSM8  = (TH2F*)fOutputContainer->At(20);
  fhPMDSM9  = (TH2F*)fOutputContainer->At(21);
  fhPMDSM10 = (TH2F*)fOutputContainer->At(22);
  fhPMDSM11 = (TH2F*)fOutputContainer->At(23);
  fhPMDSM12 = (TH2F*)fOutputContainer->At(24);
  fhPMDSM13 = (TH2F*)fOutputContainer->At(25);
  fhPMDSM14 = (TH2F*)fOutputContainer->At(26);
  fhPMDSM15 = (TH2F*)fOutputContainer->At(27);
  fhPMDSM16 = (TH2F*)fOutputContainer->At(28);
  fhPMDSM17 = (TH2F*)fOutputContainer->At(29);
  fhPMDSM18 = (TH2F*)fOutputContainer->At(30);
  fhPMDSM19 = (TH2F*)fOutputContainer->At(31);
  fhPMDSM20 = (TH2F*)fOutputContainer->At(32);
  fhPMDSM21 = (TH2F*)fOutputContainer->At(33);
  fhPMDSM22 = (TH2F*)fOutputContainer->At(34);
  fhPMDSM23 = (TH2F*)fOutputContainer->At(35);
  fhPMDSM24 = (TH2F*)fOutputContainer->At(36);
  fhPMDSM   = (TH1F*)fOutputContainer->At(37);

  gStyle->SetOptStat(110000);
  gStyle->SetOptFit(1);

  TCanvas *cPMD0 = new TCanvas("cPMD0","PMD ESD Test #1", 10,10, 600, 600);
  cPMD0->Range(-100, -100,100 ,100 );
  fhPMDSM1->SetMarkerColor(2);
  fhPMDSM1->Draw();
  fhPMDSM1->GetXaxis()->SetTitle("Cluster X");
  fhPMDSM1->GetYaxis()->SetTitle("Cluster Y");
  fhPMDSM2->SetMarkerColor(2);
  fhPMDSM2->Draw("same");
  fhPMDSM3->SetMarkerColor(2);
  fhPMDSM3->Draw("same");
  fhPMDSM4->SetMarkerColor(2);
  fhPMDSM4->Draw("same");
  fhPMDSM5->SetMarkerColor(2);
  fhPMDSM5->Draw("same");
  fhPMDSM6->SetMarkerColor(2);
  fhPMDSM6->Draw("same");
  fhPMDSM7->SetMarkerColor(4);
  fhPMDSM7->Draw("same");
  fhPMDSM8->SetMarkerColor(4);
  fhPMDSM8->Draw("same");
  fhPMDSM9->SetMarkerColor(4);
  fhPMDSM9->Draw("same");
  fhPMDSM10->SetMarkerColor(4);
  fhPMDSM10->Draw("same");
  fhPMDSM11->SetMarkerColor(4);
  fhPMDSM11->Draw("same");
  fhPMDSM12->SetMarkerColor(4);
  fhPMDSM12->Draw("same");
  fhPMDSM13->SetMarkerColor(6);
  fhPMDSM13->Draw("same");
  fhPMDSM14->SetMarkerColor(6);
  fhPMDSM14->Draw("same");
  fhPMDSM15->SetMarkerColor(6);
  fhPMDSM15->Draw("same");
  fhPMDSM16->SetMarkerColor(6);
  fhPMDSM16->Draw("same");
  fhPMDSM17->SetMarkerColor(6);
  fhPMDSM17->Draw("same");
  fhPMDSM18->SetMarkerColor(6);
  fhPMDSM18->Draw("same");
  fhPMDSM19->SetMarkerColor(8);
  fhPMDSM19->Draw("same");
  fhPMDSM20->SetMarkerColor(8);
  fhPMDSM20->Draw("same");
  fhPMDSM21->SetMarkerColor(8);
  fhPMDSM21->Draw("same");
  fhPMDSM22->SetMarkerColor(8);
  fhPMDSM22->Draw("same");
  fhPMDSM23->SetMarkerColor(8);
  fhPMDSM23->Draw("same");
  fhPMDSM24->SetMarkerColor(8);
  fhPMDSM24->Draw("same");

  DrawPMDBoundary();
  DrawPMDBoundarySM1();
  DrawPMDBoundarySM2();
  DrawPMDBoundarySM3();
  DrawPMDBoundarySM4();
  cPMD0->Print("ClusterXY.eps");
  
  TCanvas *cPMD1 = new TCanvas("cPMD1"," PMD ESD Test #2",10, 10, 600,600);
  cPMD1->Divide(1,2);
  cPMD1->cd(1);
  cPMD1->SetFillColor(0);
  fhPMDC2->SetLineColor(4);
  fhPMDC2->Draw();
  cPMD1->cd(2);
  fhPMDP2->SetLineColor(2);
  fhPMDP2->Draw();
  cPMD1->Print("CPVPREphi.eps");

  TCanvas *cPMD2 = new TCanvas("cPMD2","PMD ESD test #3",10, 10, 600, 600);
  cPMD2->cd();
  fhPMDSM->SetFillColor(2);
  fhPMDSM->Draw();
  cPMD2->Print("AllSMN.eps");

  TCanvas *cPMD3 = new TCanvas("cPMD3", "PMD ESD test #4",10, 10, 600, 600);
  cPMD3->Divide(2,2);
  cPMD3->cd(1);
  fhPMDCP0->SetMarkerColor(9);
  fhPMDCP0->Draw();
  cPMD3->cd(2);
  fhPMDCP1->SetMarkerColor(6);
  fhPMDCP1->Draw();
  cPMD3->cd(3);
  fhPMDP3->SetLineColor(2);
  fhPMDP3->Draw();
  cPMD3->cd(4);
  fhPMDCP4->SetMarkerColor(3);
  fhPMDCP4->Draw();
  cPMD3->Print("CPVPREClus.eps");

  TCanvas *cPMD4 = new TCanvas("cPMD4","PMD ESD test #5", 10, 10, 600, 600);
  cPMD4->Divide(1,2);
  cPMD4->cd(1);
  fhPMDC3->SetLineColor(4);
  fhPMDC3->Draw();
  cPMD4->cd(2);
  fhPMDP4->SetLineColor(2);
  fhPMDP4->Draw();
  cPMD4->Print("CPVPREAdc.eps");

  char line[1024] ; 
  sprintf(line, ".!tar -zcvf %s.tar.gz *.eps", GetName()) ; 
  gROOT->ProcessLine(line);
  
  AliInfo(Form("!!! All the eps files are in %s.tar.gz !!! \n", GetName())) ;
  
}

//______________________________________________________________________________
void AliPMDQATask::CalculateSMN( Float_t clsX, Float_t clsY, Int_t & smn)
{
  Double_t xcon[96] = {75.133, 54.204, 53.254, 32.326, 31.376,10.447,
		       75.133, 54.204, 53.254, 32.326, 31.376,10.447,
		       75.133, 54.204, 53.254, 32.326, 31.376,10.447,
		       75.133, 54.204, 53.254, 32.326, 31.376,10.447,
		       -75.133, -54.204, -53.254, -32.326, -31.376,-10.447,
		       -75.133, -54.204, -53.254, -32.326, -31.376,-10.447,
		       -75.133, -54.204, -53.254, -32.326, -31.376,-10.447,
		       -75.133, -54.204, -53.254, -32.326, -31.376,-10.447,
		       9.167, -32.543, -33.493, -75.133,
		       9.167, -32.543, -33.493, -75.133,
		       9.167, -32.543, -33.493, -75.133,
		       9.167, -32.543, -33.493, -75.133,
		       9.167, -32.543, -33.493, -75.133,
		       9.167, -32.543, -33.493, -75.133,
		       -9.167, 32.543, 33.493, 75.133,
		       -9.167, 32.543, 33.493, 75.133,
		       -9.167, 32.543, 33.493, 75.133,
		       -9.167, 32.543, 33.493, 75.133,
		       -9.167, 32.543, 33.493, 75.133,
		       -9.167, 32.543, 33.493, 75.133};
  
  Double_t ycon[96] =  {86.475,  86.475,  86.475, 86.475,  86.475,  86.475,
			38.225,  38.225,  38.225, 38.225,  38.225,  38.225,
			37.325,  37.325,  37.325, 37.325,  37.325,  37.325,
			-10.925, -10.925, -10.925, -10.925, -10.925, -10.925,
			-86.475, -86.475, -86.475, -86.475, -86.475, -86.475,
			-38.225,  -38.225,  -38.225, -38.225, -38.225, -38.225,
			-37.325,  -37.325,  -37.325, -37.325,  -37.325,  -37.325,
			10.925, 10.925, 10.925, 10.925, 10.925, 10.925,
			86.475,  86.475, 86.475,  86.475,
			62.225,  62.225, 62.225,  62.225,
			61.325,  61.325, 61.325,  61.325,
			37.075, 37.075, 37.075, 37.075,
			36.175,  36.175, 36.175,  36.175,
			11.925, 11.925, 11.925 , 11.925,
			-86.475,  -86.475, -86.475,  -86.475,
			-62.225,  -62.225, -62.225,  -62.225,
			-61.325,  -61.325, -61.325,  -61.325,
			-37.075,  -37.075, -37.075,  -37.075,
			-36.175,  -36.175, -36.175,  -36.175,
			-11.925, -11.925,  -11.925 , -11.925 };
  
  if((clsX <= xcon[0]) && (clsX >= xcon[1]) &&
     (clsY <= ycon[0]) && (clsY >= ycon[6])) smn = 0 ;
  
  else if((clsX <=xcon[2]) && (clsX >= xcon[3]) &&
	  (clsY <= ycon[1]) && (clsY >= ycon[7]))smn = 1 ;
  
  else if((clsX <=xcon[4]) && (clsX >= xcon[5]) &&
	  (clsY <= ycon[3]) && (clsY >= ycon[8]))smn = 2 ;
  
  else if((clsX <= xcon[0]) && (clsX >= xcon[1]) &&
	  (clsY <= ycon[12]) && (clsY >= ycon[18])) smn = 3 ;
  
  else if((clsX <=xcon[2]) && (clsX >= xcon[3]) &&
	  (clsY <= ycon[12]) && (clsY >= ycon[18]))smn = 4 ;
  
  else if((clsX <=xcon[4]) && (clsX >= xcon[5]) &&
	  (clsY <= ycon[12]) && (clsY >= ycon[18]))smn = 5 ;
  //------------------------------------------------------------------
  else if((clsX >= xcon[24]) && (clsX <= xcon[25]) &&
	  (clsY >= ycon[24]) && (clsY <= ycon[30])) smn = 6 ;
  
  else if((clsX >=xcon[26]) && (clsX <= xcon[27]) &&
	  (clsY >= ycon[25]) && (clsY <= ycon[31]))smn = 7 ;
  
  else if((clsX >=xcon[28]) && (clsX <= xcon[29]) &&
	  (clsY >= ycon[26]) && (clsY <= ycon[32]))smn = 8 ;
  
  else if((clsX >= xcon[24]) && (clsX <= xcon[25]) &&
	  (clsY >= ycon[36]) && (clsY <= ycon[42])) smn = 9 ;
  
  else if((clsX >=xcon[26]) && (clsX <= xcon[27]) &&
	  (clsY >= ycon[36]) && (clsY <= ycon[42]))smn = 10;
  
  else if((clsX >=xcon[28]) && (clsX <= xcon[29]) &&
	  (clsY >= ycon[36]) && (clsY <= ycon[42]))smn = 11;
  //------------------------------------------------------------------
  else if((clsX <= xcon[48]) && (clsX >= xcon[49]) &&
	  (clsY <= ycon[48]) && (clsY >= ycon[52])) smn = 12 ;
  
  else if((clsX <=xcon[50]) && (clsX >= xcon[51]) &&
	  (clsY <= ycon[48]) && (clsY >= ycon[52]))smn = 13 ;
  
  else if((clsX <=xcon[48]) && (clsX >= xcon[49]) &&
	  (clsY <= ycon[56]) && (clsY >= ycon[60]))smn = 14 ;
  
  else if((clsX <=xcon[50]) && (clsX >= xcon[51]) &&
	  (clsY <= ycon[56]) && (clsY >= ycon[60]))smn = 15 ;
  
  else if((clsX <=xcon[48]) && (clsX >= xcon[49]) &&
	  (clsY <= ycon[64]) && (clsY >= ycon[68]))smn = 16 ;
  
  else if((clsX <=xcon[50]) && (clsX >= xcon[51]) &&
	  (clsY <= ycon[64]) && (clsY >= ycon[68]))smn = 17 ;
  //--------------------------------------------------------------
  else if((clsX >= xcon[72]) && (clsX <= xcon[73]) &&
	  (clsY >= ycon[72]) && (clsY <= ycon[76])) smn = 18 ;
  
  else if((clsX >=xcon[74]) && (clsX <= xcon[75]) &&
	  (clsY >= ycon[72]) && (clsY <= ycon[76]))smn = 19 ;
  
  else if((clsX >=xcon[72]) && (clsX <= xcon[73]) &&
	  (clsY >= ycon[80]) && (clsY <= ycon[84]))smn = 20 ;
  
  else if((clsX >=xcon[74]) && (clsX <= xcon[75]) &&
	  (clsY >= ycon[80]) && (clsY <= ycon[84]))smn = 21;
  
  else if((clsX >= xcon[72]) && (clsX <= xcon[73]) &&
	  (clsY >= ycon[88]) && (clsY <= ycon[92])) smn = 22 ;
  
  else if((clsX >=xcon[74]) && (clsX <= xcon[75]) &&
	  (clsY >= ycon[88]) && (clsY <= ycon[92]))smn = 23 ;
  else smn = 111;
 }

//______________________________________________________________________________
void AliPMDQATask::DrawPMDBoundary()
{
  // Draw PMD boundaries 
  
  gStyle->SetLineWidth(2);
  gStyle->SetLineColor(2);
  TLine * l;
  l = new TLine(75.1333, 86.475, -75.1333, 86.475); l->Draw("same");
  l = new TLine(-75.1333, 86.470,-75.1333, -86.475); l->Draw("same");
  l = new TLine(-75.1333, -86.475,75.1333, -86.475); l->Draw("same");
  l = new TLine(75.1333, -86.475,75.1333, 86.475); l->Draw("same");
}

//______________________________________________________________________________
void AliPMDQATask::DrawPMDBoundarySM1()
{
  // Draw boundaries of Super Module 1 

  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(4);
  TLine * l;
  l = new TLine(-75.1333, 86.475, -10.447,  86.475); l->Draw("same");
  l = new TLine(-10.447, 86.475, -10.446, -10.925); l->Draw("same");
  l = new TLine(-10.446, -10.925, -75.1333,-10.925); l->Draw("same");
  l = new TLine(-75.1333,-10.925, -75.1333, 86.475); l->Draw("same");
}

//______________________________________________________________________________
void AliPMDQATask::DrawPMDBoundarySM2()
{
  // Draw boundaries of Super Module 2 

  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(4);
  TLine * l;
  l = new TLine(75.1333, -86.475, 10.446,  -86.475); l->Draw("same");
  l = new TLine(10.446,  -86.475, 10.446,  10.925); l->Draw("same");
  l = new TLine(10.446,   10.925, 75.1333, 10.925); l->Draw("same");
  l = new TLine(75.1333,  10.925, 75.1333, -86.475); l->Draw("same");
}


//______________________________________________________________________________
void AliPMDQATask::DrawPMDBoundarySM3()
{
  // Draw boundaries of Super Module 3 

  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(1);
  TLine * l;
  l = new TLine(  -9.167, 86.475, 75.1333, 86.475); l->Draw("same");
  l = new TLine(75.1333,86.475, 75.1333, 11.925); l->Draw("same");
  l = new TLine(75.1333,11.925,   -9.167,  11.925); l->Draw("same");
  l = new TLine(  -9.167, 11.925,   -9.167,  86.475); l->Draw("same");
}

//______________________________________________________________________________
void AliPMDQATask::DrawPMDBoundarySM4()
{
  // Draw boundaries of Super Module 4 

  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(1);
  TLine * l;
  l = new TLine(9.167, -86.475, -75.1333,-86.475); l->Draw("same");
  l = new TLine(-75.1333,-86.475, -75.1333,-11.925); l->Draw("same");
  l = new TLine(-75.1333,-11.925, 9.167, -11.925); l->Draw("same");
  l = new TLine(9.167, -11.925, 9.167, -86.475); l->Draw("same");
}
