// Script to create dead channel map and store them into CDB
//  - 4 sets of maps parameters can be created, with now dead channels 
// and 5%, 10%, 20% and 30% of dead channels 
//  - it reads the stored map in a given file. 
// Author: Gustavo Conesa

//.x $ALICE_ROOT/EMCAL/macros/CalibrationDB/AliEMCALSetTowerStatusCDB.C

#if !defined(__CINT__)
#include <TControlBar.h>
#include <TString.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>

#include "AliRun.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALGeoParams.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliEMCALPedestalCDB()
{
  TControlBar *menu = new TControlBar("vertical","EMCAL CDB");
  menu->AddButton("Help to run EMCAL CDB","Help()",
		  "Explains how to use EMCAL CDS menus");
  menu->AddButton("Equal Tower Status Map, all Alive","SetTowerStatusMap(0)",
		  "Set all channels to alive");
  menu->AddButton("Create Random Status Map, 5% dead","SetTowerStatusMap(5)",
		  "Set randomly 5% of the channels dead");
  menu->AddButton("Create Random Status Map, 10% dead","SetTowerStatusMap(10)",
					"Set randomly 10% of the channels dead");
  menu->AddButton("Create Random Status Map, 20% dead","SetTowerStatusMap(20)",
					"Set randomly 20% of the channels dead");
  menu->AddButton("Create Random Status Map, 30% dead","SetTowerStatusMap(30)",
					"Set randomly 30% of the channels dead");
  menu->AddButton("Set Map from txt file","SetTowerStatusMap(\"map.txt\")",
		  "Read bad channels from txt file and set them in root file");
  menu->AddButton("Read Tower Status Map","GetTowerStatusMap()",
		  "Read initial equal calibration coefficients");
  menu->Show();
}

//------------------------------------------------------------------------
void Help()
{
  char *string =
    "\nSet tower status map (dead, hot, alive) and write them into ALICE CDB. Press button \"Equal kAlive\" to set all channels alive. Press button \"Random, 5% dead\" to create random dead channel map at 5%\n";
  printf(string);
}

//------------------------------------------------------------------------
void SetTowerStatusMap(Int_t percent=0)
{
  // Writing status of all the channels in the OCDB with equal value
  // except a percent to be not alive. Right now only "alive" or "dead",
  // we need to implement the other cases like "hot"	

  TString sDBFolder ="local://PedestalsDB";
  Int_t firstRun   =  0; 
  Int_t lastRun    =  999999999;
  Int_t beamPeriod =  1;
  char* objFormat = Form("%d percent of bad channels", percent);
  
  AliCaloCalibPedestal *caloped=new AliCaloCalibPedestal(AliCaloCalibPedestal::kEmCal);
  caloped->Init();

  TObjArray map = caloped->GetDeadMap();
  printf("MAP entries %d\n",map.GetEntries());

  TRandom rn;
  //for(Int_t iSM = 0; iSM < AliEMCALGeoParams::fgkEMCALModules; iSM ++){
	for(Int_t iSM = 0; iSM < map.GetEntries(); iSM ++){
	  Int_t ndead = 0;
	  printf(" >>> SM %d <<< Entries %d, NbinsX %d, NbinsY %d\n",iSM,((TH2D*)map[iSM])->GetEntries(),((TH2D*)map[iSM])->GetNbinsX(),((TH2D*)map[iSM])->GetNbinsY());
		for(Int_t i = 0; i < ((TH2D*)map[iSM])->GetNbinsX() ; i++){
			for(Int_t j = 0; j < ((TH2D*)map[iSM])->GetNbinsY() ; j++){
				//printf("Bin (%d-%d) Content, before: %d ",i,j,((TH2D*)map[iSM])->GetBinContent(i, j));	
			    
				if(rn.Uniform(0,100) > percent)
				  caloped->SetChannelStatus(iSM, i, j, AliCaloCalibPedestal::kAlive);
				else{
					caloped->SetChannelStatus(iSM, i, j, AliCaloCalibPedestal::kDead);
					ndead++;
				}
				//printf("; after: %d \n",((TH2D*)map[iSM])->GetBinContent(i, j));	
			}	
		}
		caloped->SetDeadTowerCount(caloped->GetDeadTowerCount()+ndead);
		printf("--- dead %d\n",ndead);
  }

  printf("--- total dead %d\n",caloped->GetDeadTowerCount());

  //Store map into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Gustavo Conesa");
  
  AliCDBId id("EMCAL/Calib/Pedestals",firstRun,lastRun); // create in EMCAL/Calib/Pedestal sDBFolder 

  AliCDBManager* man = AliCDBManager::Instance();  
  AliCDBStorage* loc = man->GetStorage(sDBFolder.Data());
  loc->Put(caloped, id, &md);

}

//____________________________________________




//------------------------------------------------------------------------
void SetTowerStatusMap(char * file = "map.txt")
{
  // Get the list of dead/hot channels from file and set them in OCDB
  
  TString sDBFolder ="local://PedestalsDB";
  Int_t firstRun   =  0; 
  Int_t lastRun    =  999999999;
  Int_t beamPeriod =  1;
  char* objFormat = Form("bad channels extracted from file %s", file);
  
  AliCaloCalibPedestal *caloped=new AliCaloCalibPedestal(AliCaloCalibPedestal::kEmCal);
  caloped->Init();

  // Read parameter file line-by-line  
  ifstream f;
  f.open(file);
  
  Int_t iSM=-1, icol=-1, irow=-1, istatus=-1, ndead=0 ;
  TString string;
  if (f.good()) {
    while(string.ReadLine(f, kFALSE) && !f.eof()) {
      
      sscanf(string.Data(), "%d %d %d %d",&iSM,&icol,&irow,&istatus);
      cout<<"SM= "<<iSM<<", col= "<<icol<<", row= "<<irow<<", status="<<istatus<<endl;
      if(iSM==-1) continue;
      caloped->SetChannelStatus(iSM, icol, irow, istatus);
      ndead++;
    }
  }
  caloped->SetDeadTowerCount(ndead-2);
  printf("--- dead %d\n",ndead-2);
  
  printf("--- total dead %d\n",caloped->GetDeadTowerCount());
  
  //Store map into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Gustavo Conesa");
  
  AliCDBId id("EMCAL/Calib/Pedestals",firstRun,lastRun); // create in EMCAL/Calib/Pedestal sDBFolder 
  
  AliCDBManager* man = AliCDBManager::Instance();  
  AliCDBStorage* loc = man->GetStorage(sDBFolder.Data());
    loc->Put(caloped, id, &md);
    
}

  
//------------------------------------------------------------------------
void GetTowerStatusMap()
{
  // Read status map
 
  TString sDBFolder ="local://PedestalsDB";
  Int_t runNumber   =  0; 
  
  AliCaloCalibPedestal* caloped  = (AliCaloCalibPedestal*) 
	(AliCDBManager::Instance()
	 ->GetStorage(sDBFolder.Data())
	 ->Get("EMCAL/Calib/Pedestals", 
		   runNumber)->GetObject());

  cout<<endl;
  TObjArray map = caloped->GetDeadMap();
  printf("MAP entries %d\n",map.GetEntries());
  for(Int_t iSM = 0; iSM < map.GetEntries(); iSM ++){ 
	  TCanvas *cMap   = new TCanvas(Form("cMap%d",iSM),Form("SM %d dead map",iSM), 12,12,400,400);
	  cMap->Divide(1,1); 
 	  Int_t ndead = 0;
	  printf(" >>> SM %d <<< Entries %d, NbinsX %d, NbinsY %d\n",iSM,((TH2D*)map[iSM])->GetEntries(),((TH2D*)map[iSM])->GetNbinsX(),((TH2D*)map[iSM])->GetNbinsY());
	  for(Int_t i = 0; i < ((TH2D*)map[iSM])->GetNbinsX() ; i++){
	    for(Int_t j = 0; j < ((TH2D*)map[iSM])->GetNbinsY() ; j++){
	      if(((TH2D*)map[iSM])->GetBinContent(i, j)!=AliCaloCalibPedestal::kAlive)
		printf("Bin (%d-%d) Content: %d \n",i,j,((TH2D*)map[iSM])->GetBinContent(i, j));	
	      
	      if(((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kDead ||
		 ((TH2D*)map[iSM])->GetBinContent(i, j)==AliCaloCalibPedestal::kHot)
		ndead++;
	    }	
	  }
	  printf("--- dead %d\n",ndead);  
	  cMap->cd(iSM);
	  (TH2D*)map[iSM])->Draw("lego2");
}

printf("Total DEAD %d\n", caloped->GetDeadTowerCount());

}
