/// \ingroup macros
/// \file MUONReCalcGlobalTrigger.C
/// \brief Re-calculate regional/global trigger response from local response.
///
/// Basic usage is :
///
/// MUONReCalcGlobalTrigger("path_to_reconstruction_galice");
///
/// Starting from local responses the macro will re-calculate regional and
/// global response and print-out the global trigger decision.
/// It is used for comissioning data with cosmics where the global trigger
/// was not written in the raw stream.
/// The purpose is (for the future) to compare the re-calculated answer with 
/// the global trigger decision returned by the CTP.
/// 
/// \author Bogdan Vulpescu

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliCDBManager.h"
#include "AliMpCDB.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONGlobalTriggerBoard.h"
#include "AliMUONDataInterface.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONTriggerCrate.h"
#include "AliMUONTriggerCrateConfig.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONRegionalTriggerBoard.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONTriggerBoard.h"

#include <TArrayS.h>
#include <TObjArray.h>
#include <TMath.h>

#endif

UShort_t locResp[235]; 

AliMUONCalibrationData *calibData;
AliMUONTriggerCrateStore *fCrates;
AliMUONGlobalTriggerBoard *fGlobalTriggerBoard;
AliMUONRegionalTriggerConfig* regionalConfig;
AliMUONGlobalCrateConfig * globalConfig;

TIterator *cratesIterator;

Int_t debug;

//___________________________________________________________________________
void PrintPattBin(Short_t s) {
  /// binary print-out of the strip patterns

  printf("   ");
  Int_t mask = 0;
  for (Int_t i = 15; i >= 0; i--) {
    mask = (Int_t)TMath::Power(2,i);
    printf("%1d",(s & mask) >> i);
  }
  printf(" \n");

}

//___________________________________________________________________________
void PrintGloBin(UShort_t s) {
  /// binary print-out of global trigger decision

  Int_t mask = 0;
  for (Int_t i = 5; i >= 0; i--) {
    mask = (Int_t)TMath::Power(2,i);
    printf("%1d",(s & mask) >> i);
  }
  printf(" \n");

}

//___________________________________________________________________________
Bool_t ReCalcGlobalTrigger(TIter *nextCrates) {
  /// re-calculate regional/global decision from array of local triggers

  Int_t loLpt, loHpt;
  AliMUONTriggerCrate* cr;

  // regional response
  
  nextCrates->Reset();
  
  Int_t irb(0);
  
  while ( ( cr = static_cast<AliMUONTriggerCrate*>(nextCrates->Next()) ) ) {
    
    if (debug) printf("Crate nr = %2d \n",++irb);
    
    TObjArray *boards = cr->Boards();
    
    AliMUONRegionalTriggerBoard *regb = (AliMUONRegionalTriggerBoard*)boards->At(0);
    regb->Reset();
    
    Int_t nrBoard = 0;
    
    UShort_t regLocResp[16]; for (Int_t j=0; j<16; j++) regLocResp[j] = 0;
    
    for (Int_t j = 1; j < boards->GetEntries(); j++) {
      
      TObject *o = boards->At(j);
      
      AliMUONLocalTriggerBoard *board = (AliMUONLocalTriggerBoard*)o;
      
      if (board->GetNumber() == 0) continue;
      
      if (debug) {
	printf("...Board nr = %2d : ",++nrBoard);
	printf("%3d %s in slot %2d of crate %s \n",board->GetNumber(),board->GetName(),j,cr->GetName());
      }
      
      UShort_t response = locResp[board->GetNumber()];
      
      if (debug) printf("......Response = %x \n",response);
      
      if (response != 0) {
	loLpt =  response &  3;
	loHpt = (response & 12) >> 2;
	//printf("Response loLpt = %02b loHpt = %02b \n",loLpt,loHpt);
      }
      
      regLocResp[j-1] = response;
      
    }  // local board loop
    
    AliMUONTriggerCrateConfig* crateConfig = regionalConfig->FindTriggerCrate(cr->GetName());
    UShort_t rmask= crateConfig->GetMask();
    regb->Mask(rmask);
    regb->SetLocalResponse(regLocResp);
    regb->Response();
    //for (Int_t j=0; j<16; j++) printf("%3d ",regLocResp[j]);
    //printf("Reg %2d Response %3d mask %4x\n",irb,regb->GetResponse(),rmask);
    
    irb++;
    
  }  // crate loop
  
  // global response
  
  fGlobalTriggerBoard->Reset();
  
  if (!globalConfig)
    printf("No valid trigger crate configuration in CDB\n");

  UInt_t gmask = 0;

  for (Int_t i = 0; i < 4; i++) {
    gmask = globalConfig->GetGlobalMask(i);
    fGlobalTriggerBoard->Mask(i,gmask);
  }

  nextCrates->Reset();
  
  UShort_t regional[16];
  irb = 0;
  
  if ( !fCrates->NumberOfCrates() >= 16 ) {
    printf("Something is wrong : too many crates %d",fCrates->NumberOfCrates());
    return kFALSE;
  }
  
  for (Int_t iSide = 0; iSide < 2; iSide++) // right & left side
  {
    for (Int_t iReg = 0; iReg < 8; iReg++) // 8 crates/regional boards for each side.
    {
      cr = fCrates->Crate(iSide, iReg);

      AliMUONTriggerBoard* rb =
        static_cast<AliMUONTriggerBoard*>(cr->Boards()->At(0));
      regional[irb] = rb->GetResponse();
      ++irb;
    }
  }
  
  fGlobalTriggerBoard->SetRegionalResponse(regional);
  fGlobalTriggerBoard->Response();
  
  if (fGlobalTriggerBoard->GetResponse() != 0) {
    fGlobalTriggerBoard->Scan("");
    printf("Global trigger response = ");
    PrintGloBin(fGlobalTriggerBoard->GetResponse());
    return kTRUE;
  }
  
  return kFALSE;
  
}

//___________________________________________________________________________
void MUONReCalcGlobalTrigger(const char* input) {
  /// create array of local triggers from the raw data, run the re-calculation
  /// and print-out the results

  debug = 0;

  Int_t runNumber = 0;

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(runNumber);
  AliMpCDB::LoadDDLStore();

  calibData = new AliMUONCalibrationData(runNumber);

  regionalConfig = calibData->RegionalTriggerConfig();
  globalConfig = calibData->GlobalTriggerCrateConfig();

  fCrates = new AliMUONTriggerCrateStore;
  fCrates->ReadFromFile(calibData);
  cratesIterator = fCrates->CreateCrateIterator();
  fGlobalTriggerBoard = new AliMUONGlobalTriggerBoard;
  
  TIter nextCrates(cratesIterator);

  AliMUONDataInterface diRec(input);

  printf("Number of events = %d \n",diRec.NumberOfEvents());
  Int_t nEvents = diRec.NumberOfEvents();

  AliMUONLocalTrigger* localTrig;
  Int_t circ, loLpt, loHpt, lutLpt[2], lutHpt[2];
  TArrayS xPattern[235];
  TArrayS yPattern[235];

  for (Int_t ievent = 0; ievent < nEvents; ++ievent) {

    for (Int_t i = 0; i < 234; i++) {
      locResp[i] = 0;
    }

    AliMUONVTriggerStore* triggerStore = diRec.TriggerStore(ievent,"R");
    TIter nextLocal(triggerStore->CreateLocalIterator());
    while ( (localTrig = static_cast<AliMUONLocalTrigger*>( nextLocal() )) ) {

      if (localTrig->IsNull()) continue;

      circ = localTrig->LoCircuit();
      
      loLpt = localTrig->LoLpt();
      loHpt = localTrig->LoHpt();
      
      lutLpt[0] =  loLpt & 1;
      lutLpt[1] = (loLpt & 2) >> 1;
      lutHpt[0] =  loHpt & 1;
      lutHpt[1] = (loHpt & 2) >> 1;
      
      locResp[circ] = lutLpt[0]              +
	static_cast<int>(lutLpt[1]<<1) +
	static_cast<int>(lutHpt[0]<<2) +
	static_cast<int>(lutHpt[1]<<3);
      
      localTrig->GetXPattern(xPattern[circ]);
      localTrig->GetYPattern(yPattern[circ]);

      if (debug) {
	printf("Event %4d circ %3d loLpt %1d loHpt %1d resp %3d\n",ievent,circ,loLpt,loHpt,locResp[circ]);
      }
      
    } // local trigger loop

    if (ReCalcGlobalTrigger(&nextCrates)) {
      printf("............ for event %5d \n",ievent);
      for (Int_t ic = 1; ic <= 234; ic++) {
	if (locResp[ic] != 0) {
	  UShort_t response = locResp[ic];
	  loLpt =  response &  3;
	  loHpt = (response & 12) >> 2;
	  printf("............ in circuit %3d loLpt %1d loHpt %1d resp %3d\n",ic,loLpt,loHpt,response);
	  
	  printf("   Pattern X:\n");
	  PrintPattBin(xPattern[ic].At(0));
	  PrintPattBin(xPattern[ic].At(1));
	  PrintPattBin(xPattern[ic].At(2));
	  PrintPattBin(xPattern[ic].At(3));
	  printf("   Pattern Y:\n");
	  PrintPattBin(yPattern[ic].At(0));
	  PrintPattBin(yPattern[ic].At(1));
	  PrintPattBin(yPattern[ic].At(2));
	  PrintPattBin(yPattern[ic].At(3));
	  
	}
      }
      printf("\n\n");
    }
    
  } // event loop

  delete fGlobalTriggerBoard;
  delete fCrates;
  delete calibData;

}

