//
// This macro checks the stored reference data taken with a pedestal run.
// Execute like
//   root -b -q $ALICE_ROOT/TRD/macros/AliTRDCheckPedestal.C+\(234911\)
//
// Works for all pedestal data. In 2009 data only SM17 can be read.
// Author: Hans.Beck@cern.ch
//
#if !defined __CINT__ || defined __MAKECINT__

#include <vector>

#include <TGrid.h>
#include <TCanvas.h>
#include <TH2F.h>

#include <AliCDBManager.h>
#include <AliCDBStorage.h>
#include <AliCDBEntry.h>

#include <AliTRDCalibPadStatus.h>
#include <AliTRDCalDet.h>
#include <AliTRDCalPad.h>
#include <AliTRDCalROC.h>

#endif
//__________________________________________________________

std::vector<Int_t> GetLDCVector(const Int_t year);
std::vector<Int_t> GetSMVector(const Int_t year,const Int_t ildc);

//__________________________________________________________
void AliTRDCheckPedestal(const Int_t runNr){
  // Establish grid connection
  if(!TGrid::Connect("alien://")){printf("F-No grid connection\n");return;}
  // Set the storage to the OCDB of this runNr. It will be like
  // alien://folder=/alice/data/2016/OCDB
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorageFromRun(runNr);
  // We derive the reference storage from it
  const AliCDBStorage *stor =  man->GetDefaultStorage();
  TString folder = stor->GetBaseFolder();
  folder.ReplaceAll("OCDB","Reference");
  man->SetDefaultStorage(Form("alien://folder=%s",folder.Data()));
  // Set the run number
  man->SetRun(runNr);

  // We abuse the folder name to derive the year
  TString yearString(folder);
  yearString.ReplaceAll("/alice/data/","");
  yearString.ReplaceAll("/Reference/","");
  const Int_t year = yearString.Atoi();
  printf("W-Experimental: Derived year %d from storage folder\n",year);
  
  // The reference data is stored per Local Data Concentrator
  std::vector<Int_t> LDCvec = GetLDCVector(year);
  // Loop over LDCs
  for(std::vector<Int_t>::iterator LDCit = LDCvec.begin();LDCit!=LDCvec.end();LDCit++){
    const TString padstatus = Form("TRD/DAQData/PadStatus%d",*LDCit);
    AliCDBEntry *entry = AliCDBManager::Instance()->Get(padstatus.Data());
    AliTRDCalibPadStatus *calpadstatus = dynamic_cast<AliTRDCalibPadStatus *>(entry->GetObject());
    if(!calpadstatus){printf("E-Can not find %s in %s \n",padstatus.Data(),folder.Data());continue;}

    //Create the noise pad with the RMS values of each channel
    AliTRDCalPad *noisePad = calpadstatus->CreateCalPad();
    // LDC -> SM mapping
    std::vector<Int_t> SMvec = GetSMVector(year,*LDCit);
    for(std::vector<Int_t>::iterator SMit=SMvec.begin();SMit!=SMvec.end();SMit++){
      const TString padstatussm = Form("PadNoise-LDC%d-SM%02d",*LDCit,*SMit);
      TCanvas *cpadstatusm = new TCanvas(padstatussm.Data(),padstatussm.Data(),50,50,600,800);
      cpadstatusm->Divide(3,2);
      // Draw each layer (or here plane)
      const Float_t zRange[2]={0.,0.2};
      for(Int_t iLayer = 0;iLayer<6;iLayer++){
  	cpadstatusm->cd(iLayer+1);
	noisePad->MakeHisto2DSmPl(*SMit,iLayer,0,0,zRange[0],zRange[1],-1)->Draw("colz");
      }
      cpadstatusm->SaveAs(Form("%s.pdf",cpadstatusm->GetName()));
    } // Loop over SMs of this LDC
  } // End of loop over LDCs
} // End of void AliTRDcheckPedestal

//__________________________________________________________
std::vector<Int_t> GetLDCVector(const Int_t year){
  //
  // Returns a vector with the used LDCs for a given year
  //
  std::vector<Int_t> LDCvec;
  if(year>=2015){
    // 0..8 c++11 has std::iota which does exactly this
    for(Int_t ildc=0;ildc<9;ildc++){
      LDCvec.push_back(ildc);
    }
  } else if(year==2014){
    // No pedestal taken in 2014 (Long Shutdown)
  } else if (year == 2013 ||
	     year == 2012 ){
    // 0..3 c++11 has std::iota which does exactly this
    for(Int_t ildc=0;ildc<4;ildc++){
      LDCvec.push_back(ildc);
    }
  } else if (year == 2011) {
    // 0..9 c++11 has std::iota which does exactly this
    for(Int_t ildc=0;ildc<10;ildc++){
      LDCvec.push_back(ildc);
    }
  } else if (year == 2010) {
    // 1, 3, 4, 6
    LDCvec.push_back(1);
    LDCvec.push_back(3);
    LDCvec.push_back(4);
    LDCvec.push_back(6);
  } else if (year == 2009) {
    LDCvec.push_back(3);
    // 1, 4, 6 exist but do not work?
    // LDCvec.push_back(1);
    // LDCvec.push_back(4);
    // LDCvec.push_back(6);
  }
  else {
    printf("E-No data earlier than 2009\n");
  }
  return LDCvec;
}
//__________________________________________________________
std::vector<Int_t> GetSMVector(const Int_t year,const Int_t ildc){
  //
  // Returns a vector with the supermodules read out
  // by a given LDC in a given year
  //
  std::vector<Int_t> SMvec;
  
  if(year>=2015){
    // 2015, 2016
    SMvec.push_back(2*ildc);
    SMvec.push_back(2*ildc +1);
  } else if(year==2014){
    // No pedestal taken in 2014 (Long Shutdown)
  } else if (year == 2013 ||
	     year == 2012){
    // 2013, 2012
    if(ildc==0){
      SMvec.push_back(0);
      SMvec.push_back(1);
      SMvec.push_back(2);
      SMvec.push_back(3);
    } else if(ildc==1){
      SMvec.push_back(6);
      SMvec.push_back(7);
      SMvec.push_back(8);
    } else if(ildc==2){
      SMvec.push_back(9);
      SMvec.push_back(10);
      SMvec.push_back(11);
    } else if(ildc==3){
      SMvec.push_back(15);
      SMvec.push_back(16);
      SMvec.push_back(17);
    }
  } else if (year == 2011){
    // 2011
    if(ildc==0){
      SMvec.push_back(0);
    } else if(ildc==1){
      SMvec.push_back(1);
    } else if(ildc==2){
      SMvec.push_back(7);
    } else if(ildc==3){
      SMvec.push_back(8);
    } else if(ildc==4){
      SMvec.push_back(9);
    } else if(ildc==5){
      SMvec.push_back(10);
    } else if(ildc==6){
      SMvec.push_back(11);
    } else if(ildc==7){
      SMvec.push_back(15);
    } else if(ildc==8){
      SMvec.push_back(16);
    } else if(ildc==9){
      SMvec.push_back(17);
    }
  } else if (year == 2010) {
    // 2010
    if(ildc==1){
      SMvec.push_back(0);
      SMvec.push_back(1);
    } else if(ildc==3){
      SMvec.push_back(7);
      SMvec.push_back(8);
    } else if(ildc==4){
      SMvec.push_back(9);
      SMvec.push_back(10);
    } else if(ildc==6){
      SMvec.push_back(17);
    }
  } else if (year == 2009) {
    // 2009
    if(ildc==3){
      SMvec.push_back(17);
    }
    // 1,4,6 exist but do not work?
  } else {
    printf("E-No data earlier than 2009\n");
  }

  return SMvec;
}
