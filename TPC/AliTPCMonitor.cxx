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

/*
$Log$
Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/   

#include "AliTPCMonitor.h"
#ifdef ALI_DATE
#include "AliTPCMonitorDateMonitor.h"
#endif
ClassImp(AliTPCMonitor)

//____________________________________________________________________________
AliTPCMonitor::AliTPCMonitor(char* name, char* title) : AliTPCMonitorConfig(name,title)
{
  // Constructor 
  
  fPad       = new Int_t*[GetMaxHwAddr()];     
  for(Int_t i = 0; i<GetMaxHwAddr(); i++) { fPad[i]       = new Int_t[GetTimeBins()];}
  
  fPadMapRCU = new Int_t*[GetMaxHwAddr()];     
  for(Int_t i = 0; i<GetMaxHwAddr(); i++) { fPadMapRCU[i] = new Int_t[6];}
  
  fPadMapHw  = new Float_t[GetMaxHwAddr()];
  
 

  fReaderROOT                  =   0;  
  fReaderDATE                  =   0; 
  #ifdef ALI_DATE
  fReaderDATEMon               =   0; 
  #endif
  // row and pad settings 
  kNRowsIroc                   =  63;
  kNRowsOroc                   =  96;
  
  kNPadsIroc                   = 110;
  kNPadsOroc                   = 140; 
  kNPadMinIroc                 = -55;
  kNPadMaxIroc                 =  55;
  
  kNPadMinOroc                 = -70;
  kNPadMaxOroc                 =  70;
  
    
    
  // histograms to be used      ////////////////////////////////////////////////////////////////
  fHistIROC                    =0;
  fHistOROC                    =0;
  fHistIROCIndex               =0;
  fHistOROCIndex               =0;
  fHistIROCTime                =0; 
  fHistOROCTime                =0;
  fHistIROCClone               =0;
  fHistOROCClone               =0;
    
  fHistIROCRMS                 =0;
  fHistOROCRMS                 =0;

  fHistIROCBASE                =0;
  fHistOROCBASE                =0;

  fHistIROCSUM                 =0;
  fHistOROCSUM                 =0;
  
  fHistChannelTime             =0;
  fHistAddrMapIndex            =0;
  fHistAddrMaxAdc              =0;
  fHistAddrBaseMean            =0;
  fHistAddrMaxAdcX             =0;
  fHistAddrAdcSum              =0;
  fHistAddrBaseRms             =0;

  fHistDistrSumIROC            =0;
  fHistDistrMaxIROC            =0;
    
  fHistDistrSumOROC            =0;
  fHistDistrMaxOROC            =0;
    
  fHistDistrBase2dIROC         =0;
  fHistDistrBase2dOROC         =0;
  
  fHistDistrBaseRmsIROC        =0;
  fHistDistrBaseMeanIROC       =0;
  
  fHistDistrBaseRmsOROC        =0;
  fHistDistrBaseMeanOROC       =0;

  fHistGlobalMaxA              =0;
  fHistGlobalMaxC              =0;


  // some error flags //////////////////////////////////////////////////////////////////////////////
  fVerb                        =0;
  fMapHand                     =0;


  fLastEv                      =0;
  fEventNumber                 =0;
  fEventNumberOld              =0;
  
  fDisableFit                  =0;
	 
  fExecGlob                    =0;
  fExecPlaneMax                =0;
  fExecPadIrocRms              =0;
  fExecPadOrocRms              =0;
    
  fRunId                       =0;
  fEqId                        =0;
  
  fPadUsedRoc                  =-1;
  fPadUsedHwAddr               =-1;
    
  fMirror                      = 1;
  
  fChannelIter                 = 0;
 

  
  fMapEqidsRcu   = new Int_t[1000];
  fMapEqidsSec   = new Int_t*[36];
  for(Int_t i = 0; i<36; i++) { fMapEqidsSec[i]       = new Int_t[6];}
  
  
  if (gDirectory) 
    {
      if (!gDirectory->GetList()) 
	{
	  Warning("Build","Current directory is not a valid directory");
	  return;
	}
      AliTPCMonitor *hold = (AliTPCMonitor*)gDirectory->GetList()->FindObject(GetName());
      if(hold) 
	{
	  Warning("Build","Replacing existing histogram: %s (Potential memory leak).",GetName());
	  gDirectory->GetList()->Remove(hold);
	}
      gDirectory->Append(this);
     }
  
  fHistList = new TObjArray();

  fGdcId    =  0;
  fLdcId    =  0;
  fLdcIdOld = -1;

  CreateHistos();
  SetEqIds();

}
//____________________________________________________________________________
AliTPCMonitor::~AliTPCMonitor() 
{
  // Destructor
  
  for(Int_t i = 0; i<GetMaxHwAddr(); i++) { delete[] fPad[i] ;}
  for(Int_t i = 0; i<GetMaxHwAddr(); i++) { delete[] fPadMapRCU[i];}
  delete[] fPadMapHw ;
  DeleteHistos();
}

//____________________________________________________________________________
void AliTPCMonitor::CreateHistos() 
{
  // Create histograms to be displayed 

  if(fVerb) cout << " create new ones " << endl;
  fHistIROC            = new TH2F("fHistIROC"            ,"fHistIROC"           ,kNRowsIroc,0,kNRowsIroc,kNPadsIroc, kNPadMinIroc, kNPadMaxIroc);       fHistList->Add(fHistIROC);
  fHistOROC            = new TH2F("fHistOROC"            ,"fHistOROC"           ,kNRowsOroc,0,kNRowsOroc,kNPadsOroc, kNPadMinOroc, kNPadMaxOroc);       fHistList->Add(fHistOROC); 

  fHistIROCIndex       = new TH2S("fHistIROCIndex"       ,"fHistIROCIndex"      ,kNRowsIroc,0,kNRowsIroc,kNPadsIroc, kNPadMinIroc, kNPadMaxIroc);       fHistList->Add(fHistIROCIndex);
  fHistOROCIndex       = new TH2S("fHistOROCIndex"       ,"fHistOROCIndex"      ,kNRowsOroc,0,kNRowsOroc,kNPadsOroc, kNPadMinOroc, kNPadMaxOroc);       fHistList->Add(fHistOROCIndex);

  fHistIROCTime        = new TH2F("fHistIROCTime"        ,"fHistIROCTime"       ,kNRowsIroc,0,kNRowsIroc,kNPadsIroc, kNPadMinIroc, kNPadMaxIroc);       fHistList->Add(fHistIROCTime);
  fHistOROCTime        = new TH2F("fHistOROCTime"        ,"fHistOROCTime"       ,kNRowsOroc,0,kNRowsOroc,kNPadsOroc, kNPadMinOroc, kNPadMaxOroc);       fHistList->Add(fHistOROCTime);
  
  fHistIROCRMS         = new TH2F("fHistIROCRMS"         ,"fHistIROCRMS"        ,kNRowsIroc,0,kNRowsIroc,kNPadsIroc, kNPadMinIroc, kNPadMaxIroc);       fHistList->Add(fHistIROCRMS);
  fHistOROCRMS         = new TH2F("fHistOROCRMS"         ,"fHistOROCRMS"        ,kNRowsOroc,0,kNRowsOroc,kNPadsOroc, kNPadMinOroc, kNPadMaxOroc);       fHistList->Add(fHistOROCRMS);
   
  fHistIROCSUM         = new TH2F("fHistIROCSUM"         ,"fHistIROCSUM"        ,kNRowsIroc,0,kNRowsIroc,kNPadsIroc, kNPadMinIroc, kNPadMaxIroc);       fHistList->Add(fHistIROCSUM);
  fHistOROCSUM         = new TH2F("fHistOROCSUM"         ,"fHistOROCSUM"        ,kNRowsOroc,0,kNRowsOroc,kNPadsOroc, kNPadMinOroc, kNPadMaxOroc);       fHistList->Add(fHistOROCSUM);
 
  fHistIROCBASE        = new TH2F("fHistIROCBASE"        ,"fHistIROCBASE"       ,kNRowsIroc,0,kNRowsIroc,kNPadsIroc, kNPadMinIroc, kNPadMaxIroc);       fHistList->Add(fHistIROCBASE);
  fHistOROCBASE        = new TH2F("fHistOROCBASE"        ,"fHistOROCBASE"       ,kNRowsOroc,0,kNRowsOroc,kNPadsOroc, kNPadMinOroc, kNPadMaxOroc);       fHistList->Add(fHistOROCBASE);


  fHistChannelTime     = new TH2F("fHistChannelTime"     ,"fHistChannelTime"    ,GetNumOfChannels(),0,GetNumOfChannels(),GetTimeBins(),0,GetTimeBins());fHistList->Add(fHistChannelTime);
  fHistAddrMapIndex    = new TH1F("fHistAddrMapIndex"    ,"fHistAddrMapIndex"   ,GetMaxHwAddr()    ,0,GetMaxHwAddr());                                  fHistList->Add(fHistAddrMapIndex);
  fHistAddrMaxAdc      = new TH1F("fHistAddrMaxAdc"      ,"fHistAddrMaxAdc"     ,GetMaxHwAddr(),0,GetMaxHwAddr());                                      fHistList->Add(fHistAddrMaxAdc);
  fHistAddrMaxAdcX     = new TH1F("fHistAddrMaxAdcX"     ,"fHistAddrMaxAdcX"    ,GetMaxHwAddr(),0,GetMaxHwAddr());                                      fHistList->Add(fHistAddrMaxAdcX);
  fHistAddrBaseMean    = new TH1F("fHistAddrBaseMean"    ,"fHistAddrBaseMean"   ,GetMaxHwAddr(),0,GetMaxHwAddr());                                      fHistList->Add(fHistAddrBaseMean);
  fHistAddrAdcSum      = new TH1F("fHistAddrAdcSum"      ,"fHistAddrAdcSum"     ,GetMaxHwAddr(),0,GetMaxHwAddr());                                      fHistList->Add(fHistAddrAdcSum);
  fHistAddrBaseRms     = new TH1F("fHistAddrBaseRms"     ,"fHistAddrBaseRms"    ,GetMaxHwAddr(),0,GetMaxHwAddr());                                      fHistList->Add(fHistAddrBaseRms);
  fHistDistrSumIROC    = new TH1F("fHistDistrSumIROC"    ,"fHistDistrSumIROC"   ,400,0.0,4000.0);                                                       fHistList->Add(fHistDistrSumIROC);
  fHistDistrMaxIROC    = new TH1F("fHistDistrMaxIROC"    ,"fHistDistrMaxIROC"   ,500,0.0,1000.0);                                                       fHistList->Add(fHistDistrMaxIROC);
  fHistDistrSumOROC    = new TH1F("fHistDistrSumOROC"    ,"fHistDistrSumOROC"   ,400,0.0,4000.0);                                                       fHistList->Add(fHistDistrSumOROC);
  fHistDistrMaxOROC    = new TH1F("fHistDistrMaxOROC"    ,"fHistDistrMaxOROC"   ,500,0.0,1000.0);                                                       fHistList->Add(fHistDistrMaxOROC);
   
  fHistDistrBase2dIROC = new TH2F("fHistDistrBase2dIROC" ,"fHistDistrBase2dIROC",100,0.0,100.0,100,0.0,10.0);                                           fHistList->Add(fHistDistrBase2dIROC);
  fHistDistrBase2dOROC = new TH2F("fHistDistrBase2dOROC" ,"fHistDistrBase2dOROC",100,0.0,100.0,100,0.0,10.0);                                           fHistList->Add(fHistDistrBase2dOROC);
  
  fHistGlobalMaxA      = new TH2S("SIDE A"               ,"SIDE A"              ,500,-3000,3000,500,-3000,3000);                                        fHistList->Add(fHistGlobalMaxA);
  fHistGlobalMaxC      = new TH2S("SIDE C"               ,"SIDE C"              ,500,-3000,3000,500,-3000,3000);                                        fHistList->Add(fHistGlobalMaxC);
  
  ResetArrays();
}
//____________________________________________________________________________
Int_t AliTPCMonitor::ProcessEvent()
{
  // Process Event
  // Depending on the value of the sector id all sectors (sectorid == -1) are processed.
  //
  // In this case only the maximum values are calculated per pad and filled to the global histograms
  // In a second loop the last processed(displayed) sector will be processed (sectorid!=-1) 
  // again and the baseline rms and further quantities are calculated
  //
  // If only one sector should be processed SetProcOneSector(1) should be set.
  // In this case only the specified (last/last displayed) sector will be processed.
  //
  // If GetProcNextEvent()==0 the same event will be processed again


  Int_t sectorid = 0;
  Int_t retflag  = 0; // id of last sector + 1000,  or error flag
  if(GetProcNextEvent()==1 && fLastEv) { AliInfo("Last event already processed"); }
  if(GetProcNextEvent()==1) ResetSectorArray();
  
  
  if(GetProcNextEvent()==0 || GetProcOneSector()==1 ) sectorid = GetLastSector();
  else                                                sectorid = -1;
  
  // first iteration 
  retflag = ReadData(sectorid);
  
  SetLastProcFile(GetFile());
  
  if(retflag>=10 && retflag<1000){ AliError("Could not read event properly: Check file name and format or try next event");  return 0  ;}
  
  DrawHists(3);
  
  // second iteration 
  if(sectorid==-1 && retflag >1000) 
    { 
      AliInfo("Second read cycle");
      SetProcNextEvent(0);
      if(GetLastSectorDisplayed()==-1) {sectorid =  GetLastSector()         ;    } 
      else                             {sectorid =  GetLastSectorDisplayed();  SetLastSector(sectorid)           ;  } 
      retflag = ReadData(sectorid);
    }
  
  SetLastSectorDisplayed(sectorid)  ;
  fMapHand->ReadfecHwMap(GetLastSector());
  FillHistsPadPlane();
  DrawHists(1);
  SetEventProcessed(1);
  return retflag;
}

//__________________________________________________________________
Int_t AliTPCMonitor::ReadData(Int_t secid)
{
  // Read Data File/Stream  for specified Format.
  // Payload will be extracted from either ROOT or DATE format 
  // and passed to FillHistsDecode for decoding of the adc information 
  
  Int_t format = GetFormat();
  

  //if(format==2 && !gSystem->Getenv("DATE_ROOT")) {   AliError("DATE not initialized on this system"); return  11;}
  
  if(     format==0) {return  ReadDataROOT(secid);}
  else if(format==1) {return  ReadDataDATEFile(secid);}
  else if(format==2) 
    {
      return  ReadDataDATEStream(secid);
    }
  
  AliWarning("Function should already be left");
  return 11;
} 
//__________________________________________________________________
Int_t AliTPCMonitor::ReadDataDATEFile(Int_t secid) 
{
  // Read Data in Date format.
  // Date file and monitor will be handled in this function.
  
  if(fReaderROOT) { delete fReaderROOT ; fReaderROOT=0;}
  
#ifdef ALI_DATE
  if(fReaderDATEMon) { delete fReaderDATEMon ; fReaderDATEMon=0;}
#endif
  
  Char_t*                  eventPtr = 0;
  AliTPCMonitorDateFormat* DateForm = 0;
  
  // Create objects //
  if( fReaderDATE==0 || ( fReaderDATE!=0 && (strcmp(GetLastProcFile(),GetFile())!=0) )  ) 
    {
      cout << " Create new DATE file "<< endl; 
      if(fReaderDATE)  { delete fReaderDATE ; fReaderDATE=0; }
      fReaderDATE = new AliTPCMonitorDateFile() ;
      fReaderDATE->SetName("fReaderDATE") ;
      fReaderDATE->OpenDateFile(GetFile()) ;
      fEventNumber =0;
    }
  
  // Get Event pointer ///
  if(fReaderDATE->IsDateFileOpen() ==false )  { AliError("Could not open Date File"); return  11 ; }
  // Rewind event if new event number is smaller than old one
  if(fEventNumberOld>fEventNumber)           { fReaderDATE->ResetFilePos(); }
  while(GetProcNextEvent() || fEventNumber==0)
    {
      fReaderDATE->ReadEvent();
      eventPtr  = fReaderDATE->GetMemoryPointer();
      DateForm = new AliTPCMonitorDateFormat(eventPtr);
      Int_t current_ev =  DateForm->GetEventID();
      if(fEventNumber <= current_ev ){  break; }
    }
  
  eventPtr  = fReaderDATE->GetMemoryPointer();
  
  if(DateForm==0) DateForm = new AliTPCMonitorDateFormat(eventPtr); 
  fEventNumber     =  DateForm->GetEventID();
  fEventNumberOld  =  DateForm->GetEventID();
  
  if(fReaderDATE->IsLastEvent()) fLastEv =1;
  if(DateForm->IsEventEndOfRun()) {  AliInfo("Event is end of run: eventType END_OF_RUN "); }
  
  if(fReaderDATE->IsEventValid()       == false )   {AliInfo("Skipping Event because invalid");                     return 10;   } 
  if(DateForm->IsEventPhysicsEvent()   == false )   {AliInfo("Skipping Header/event because not Physics event");    return 12;   } 
  
  ResetHistos() ; 
  
  Int_t lastrcu = ReadDataDATESubEventLoop(DateForm,secid);
  
  delete DateForm;
  if(fVerb) cout << " last rcu " << lastrcu << endl; 
  return lastrcu;
}


#ifdef ALI_DATE
//__________________________________________________________________
Int_t AliTPCMonitor::ReadDataDATEStream(Int_t secid) 
{
  // Read Data from DATE stream.
  // Can also be used for DATE file.
 
  
  if(fReaderROOT) { delete fReaderROOT ; fReaderROOT=0;}
  if(fReaderDATE) { delete fReaderDATE    ; fReaderDATE   =0; }
  
  Char_t*                  eventPtr         = 0;
  
  AliTPCMonitorDateFormat* DateForm =0 ;
  
  // Create objects ///
  if((fReaderDATEMon==0 || ( fReaderDATEMon!=0 && (strcmp(GetLastProcFile(),GetFile())!=0) )))
    {
      if(fReaderDATEMon!=0) 
	{
	  fReaderDATEMon->Logout();
	  delete fReaderDATEMon;
	}
      fReaderDATEMon  = new AliTPCMonitorDateMonitor();
      fReaderDATEMon->SetName("DMon");
      Int_t status = fReaderDATEMon->OpenMonitoring(GetFile());
      if(status) {  AliError(Form("Could not read event online: Error: %s",fReaderDATEMon->DecodeError(status))); return 11; }
    }
  
  if(GetProcNextEvent())
    {
      fReaderDATEMon->Free();
      Int_t status = fReaderDATEMon->GetEvent();
      if(status) {  AliError(Form("Could not read event online: Error: %s",fReaderDATEMon->DecodeError(status))); return 11 ;} 
      //fReaderDATEMon->Logout();
    }
  
  eventPtr  = fReaderDATEMon->GetEventPointerasChar();
  
  if(DateForm==0) DateForm = new AliTPCMonitorDateFormat(eventPtr); 
  fEventNumber     =  DateForm->GetEventID();
  fEventNumberOld  =  DateForm->GetEventID();
  
  if(DateForm->IsEventEndOfRun()) {  AliInfo("Event is end of run: eventType END_OF_RUN "); }
  if(DateForm->IsEventPhysicsEvent()      == false  )  {AliInfo("Skipping Header/event because not Physics event");    return 12;   } 
  
  ResetHistos() ; 
  
  Int_t lastrcu = ReadDataDATESubEventLoop(DateForm,secid);
  
  delete DateForm;
  if(fVerb) cout << " last rcu " << lastrcu << endl; 
  return lastrcu;
}
#else
//__________________________________________________________________
Int_t AliTPCMonitor::ReadDataDATEStream(Int_t /*secid*/) 
{
  // Read Data from DATE stream.
  // Can also be used for DATE file.
  // In case DATE is not install
  // this method is dummy

  AliError("DATE not initialized on this system"); 
  return  11;
}
#endif

Int_t AliTPCMonitor::ReadDataDATESubEventLoop(AliTPCMonitorDateFormat* DateForm, Int_t secid)
{
  // Loop over DATE Subevents  

  Bool_t                   exitSubEventLoop = false;
  Int_t                    last_rcu_id      = 0; 
  Int_t                    last_eq          = 0;
  Int_t                    start            = 1;
  Int_t                    rcu_patch        = 0;
  Char_t*                  eventPtr         = 0;
  UInt_t*                  eventPtrUI       = 0;

  fChannelIter =0;

  while(exitSubEventLoop == false)
    {
      if(start==1) { DateForm->GotoSubEventHeader()    ;  start=0;}
      else         { DateForm->GotoNextSubEventHeader();          }
      if(DateForm->IsLastSubEventHeader()) exitSubEventLoop = true;
      
      if(fVerb) cout << " next subevent LDC " << (Int_t)DateForm->GetSubEventLDC() <<  endl;
      
      last_eq    = 0 ;
      Int_t neq  = 0 ;
      
      while(last_eq==0) 
	{
	  if(neq ==0){ DateForm->GotoFirstEquipment();} 
	  else       { DateForm->GotoNextEquipment();}
	  
	  fGdcId      = DateForm->GetSubEventGDC() ;
	  fLdcId      = DateForm->GetSubEventLDC();
	  fRunId      = DateForm->GetEventRunID(); 
	  fEqId       = DateForm->GetEquipmentID();
	  rcu_patch   = GetRCUPatch(fRunId, fEqId);
	  last_rcu_id = (rcu_patch+1000);
	  neq++;
	  
	  if(fLdcIdOld!=fLdcId &&  fChannelIter!=0) {
	    if(secid==-1)
	      {
		FillGlobal(GetLastSector());
		ResetArrays();	      
		fChannelIter =0;
	      }
	    else 
	      {
		return last_rcu_id;
	      }
	  }
	  fLdcIdOld = fLdcId ;
	  
	  
	  if(DateForm->IsLastEquipment() != false )  last_eq = 1;    
	  
	  eventPtr = DateForm->GetFirstDataPointer();
	  eventPtrUI = (UInt_t *) eventPtr;
	  Int_t payload = DateForm->GetPayloadSize();
	  
	  if(fVerb)DumpHeader(DateForm);	  
	  if(fVerb) cout << "Check sector and fEqId  " << endl;
	  
	  if(fVerb && secid >0 ) 
	    {
	      cout << "secid : "<< secid  << " fEqId "<<  fEqId << " equipment 0 form map  "  << endl;
	      cout <<  fMapEqidsSec[secid][0] <<  endl;
	      cout << " showed map_eqids " << endl;
	    }
	  
	  if(CheckEqId(secid,fEqId)) 
	    {
	      if(fVerb) cout << " init altro " << endl;
	      AliTPCMonitorAltro* altro  = new AliTPCMonitorAltro((UInt_t *)eventPtrUI,(payload/4),1);
	      altro->SetWrite10Bit(GetWrite10Bit());
	      altro->SetActFilename(GetFile());
	      if(fVerb) cout << " allocated 10bit " << endl;
	      altro->Allocate10BitArray();
	      altro->Decodeto10Bit(fEqId);
	      AliInfo(Form("Process eqid %i , patch %i ",fEqId,rcu_patch%6));
	      FillHistsDecode(altro,(rcu_patch%6),secid);
	      if(fChannelIter!=0 && secid==-1 ) SetSectorFilled(rcu_patch/6);
	      delete altro;
	    }// if(CheckID)
	}// while last eq
      SetLastSector(rcu_patch/6);
    }
  if(fChannelIter!=0 && secid==-1){ FillGlobal(GetLastSector()); }  
  return last_rcu_id;
  
}
//__________________________________________________________________
Int_t AliTPCMonitor::ReadDataROOT(Int_t secid) 
{
  // Read in data in ROOT format 
  if(fReaderDATE){ delete fReaderDATE ; fReaderDATE=0;}
  fChannelIter =0;
  Int_t           last_rcu_id   = 0; 
  Int_t           rcu_patch     = 0;
  Int_t           equipmentSize = 0;
  UChar_t*        eventPtr      = 0 ; 
  UInt_t*         eventPtrUI    = 0 ;
  Int_t           evtype        = 0;

  if(fReaderROOT==0 || ( fReaderROOT!=0 && (strcmp(GetLastProcFile(),GetFile())!=0) )  ) 
    {
      if(fVerb) cout << "AliTPCMonitor::ReadDataROOT create new reader " << endl;
      if(fReaderROOT)  { delete fReaderROOT ; fReaderROOT=0; }
      fReaderROOT = new AliRawReaderRoot(GetFile());
      if(!fReaderROOT) { AliError("Could not initiate  AliRawReaderRoo "); return 10;}
      fEventNumber=0;
    }
 
  // reset to beginning of the event
  fReaderROOT->Reset();
  
  // Rewind event if new event number is smaller than old one
  if(fEventNumberOld>fEventNumber) fReaderROOT->RewindEvents(); 
  
  while(GetProcNextEvent() || fEventNumber==0)
    {
      if(fVerb) cout << "AliTPCMonitor::ReadDataROOT get event " << endl;
      if(!fReaderROOT->NextEvent()) { AliError("Could not get next Event"); return 11 ;}
      Int_t current_ev =  *(fReaderROOT->GetEventId());
     
      if(fEventNumber <= current_ev ){  break; }
    }
  
  fEventNumber     =  *(fReaderROOT->GetEventId());
  fEventNumberOld  =  *(fReaderROOT->GetEventId()); 
  
  ResetHistos() ; 
  
  while(fReaderROOT->ReadHeader()) 
    {
      if(fVerb) cout << "AliTPCMonitor::ReadDataROOT read header " << endl;
      fGdcId        = fReaderROOT->GetGDCId() ;
      fLdcId        = fReaderROOT->GetLDCId() ;
      fRunId        = fReaderROOT->GetRunNumber() ;
      equipmentSize = fReaderROOT->GetEquipmentSize();
      fEqId         = fReaderROOT->GetEquipmentId();
      evtype        = fReaderROOT->GetType();
      rcu_patch     = GetRCUPatch(fRunId, fEqId);
      last_rcu_id   = (rcu_patch+1000);
      
      if(evtype==1)        { AliWarning(Form("EventType==1 in event %i ",fEventNumber)); return 10; }
      if(equipmentSize==0) { AliWarning(Form("Equipmentsize ==0  in event %i ",fEventNumber)); return 10; }
      
      if(fVerb) DumpHeader(fReaderROOT);
      
      if(fLdcIdOld!=fLdcId &&  fChannelIter!=0) {
	if(secid==-1)
	  {
	    FillGlobal(GetLastSector());
	    ResetArrays();	      
	    fChannelIter =0;
	  }
	else 
	  {
	    return last_rcu_id;
	  }
      }
      fLdcIdOld = fLdcId ;
      
      if(CheckEqId(secid,fEqId))
	{
	  fReaderROOT->ReadNextData(eventPtr); 
	  eventPtrUI = (UInt_t *) eventPtr;
	  Int_t offset = 16;
	  Int_t fsize = (Int_t)((equipmentSize/4) -offset) +1;
	  AliTPCMonitorAltro* altro  = new AliTPCMonitorAltro((UInt_t *)eventPtrUI,fsize,2);
	  altro->SetWrite10Bit(GetWrite10Bit());
	  if(fVerb) cout << "AliTPCMonitor::ReadDataROOT Set Write10bit to " << GetWrite10Bit() << endl;
	  altro->SetActFilename(GetFile());
	  altro->Allocate10BitArray();
	  altro->Decodeto10Bit(fEqId);
	  AliInfo(Form("Process sector %i eqid %i , patch %i ",rcu_patch/6,fEqId,rcu_patch%6));
	  FillHistsDecode(altro,(rcu_patch%6),secid);
	  if(fChannelIter!=0 && secid==-1 ) SetSectorFilled(rcu_patch/6);
	  delete altro;
	}
      SetLastSector(rcu_patch/6);
    }
  if(fChannelIter!=0 && secid==-1) { FillGlobal(GetLastSector()); }
  return last_rcu_id;
}

 
//____________________________________________________________________________
void AliTPCMonitor::FillHistsDecode(AliTPCMonitorAltro* altro ,Int_t rcu_patch, Int_t secid) 
{
  // Decode Channels, calculate base mean and rms depending on the 
  // value of secid (see ProcessEvent) and fill histograms


  if(fVerb)   cout << "FillHistsDecode : rcu_patch " << rcu_patch << " id " << secid <<  endl;
  Int_t     time_stamp        = 0;
  Int_t     sample_iter       = 0;
  Int_t     sample_length     = 0;
  Int_t     sample_bins       = 0;
  Float_t   max               = 0;
  Float_t   max_x             = 0;
  Float_t   sum               = 0.0;
  Int_t     sumn              = 0;
  Int_t     block_pos         = 0;
  Int_t     hw                = 0;
  Int_t     nwords            = 0; 
  Int_t     nextHwAddress     = 0;
  Int_t     ntime             = 0;
  Int_t     adc               = 0;
  Int_t     nextpos           = 0;
  Int_t     lastpos           = altro->Get10BitArraySize()-1;
  Short_t*  entries           = altro->Get10BitArray();
  Double_t  h_rms             = 0.0;
  Double_t  h_mean            = 0.0;
  Int_t     supnextpos        = 0;
  TH1D*     hbase             =  new TH1D("hbase","hbase",GetTimeBins(),0.5,(GetTimeBins()+0.5));
 
  while(lastpos>0) 
    {
      nextpos    = altro->DecodeTrailer(lastpos);
      supnextpos = altro->GetNextTrailerPos();
      if(nextpos==-1) { break; }
      
      lastpos               = nextpos;
      block_pos             = altro->GetTrailerBlockPos();
      hw                    = altro->GetTrailerHwAddress(); 
      nwords                = altro->GetTrailerNWords();
      nextHwAddress         = ( hw + (rcu_patch<<12) );
      fPad[fChannelIter][0] = nextHwAddress ;
      
      if(fPadMapHw[nextHwAddress]!=-1 ) 
	{
	  //Int_t hw_before1 = fPad[fChannelIter-2][0];
	  //Int_t hw_before2 = fPad[fChannelIter-3][0];
	  
	  if(fVerb){ cout  <<"\n //// Ambiguous hwaddress "   << nextHwAddress << "  write 10bit files and check file for eqid : "  << fEqId << " /////// " << endl;}
	  
	  if( TMath::Abs(fPadMapRCU[nextHwAddress][4] - fChannelIter)==1) 
	    {
	      if(fVerb) cout << "//// Difference to previous channel==1 : reset branch bit of hw from last channel to 1 " << endl;
	      Int_t hwreset =  (nextHwAddress + (1<<11));
	      fPad[fChannelIter-1][0] = hwreset;
	      
	      fPadMapHw[hwreset]    =  fChannelIter-1  ;
	      fPadMapRCU[hwreset][0]=  rcu_patch      ;
	      fPadMapRCU[hwreset][1]=  ((hwreset & AliTPCMonitorAltro::GetHwMaskBranch())    >> 11);
	      fPadMapRCU[hwreset][2]=  ((hwreset & AliTPCMonitorAltro::GetHwMaskFEC())       >>7  );
	      fPadMapRCU[hwreset][3]=  ( hwreset & AliTPCMonitorAltro::GetHwMaskFECChannel()      );
	      fPadMapRCU[hwreset][4]=  fChannelIter-1;
	      fPadMapRCU[hwreset][5]=  altro->GetTrailerPos();
	    }
	}
      
      fPadMapHw[nextHwAddress]    =  fChannelIter  ;
      fPadMapRCU[nextHwAddress][0]=  rcu_patch     ;
      fPadMapRCU[nextHwAddress][1]=  ((nextHwAddress &  AliTPCMonitorAltro::GetHwMaskBranch())>> 11)    ;
      fPadMapRCU[nextHwAddress][2]=  ((nextHwAddress &  AliTPCMonitorAltro::GetHwMaskFEC())   >>7);
      fPadMapRCU[nextHwAddress][3]=  (nextHwAddress  &  AliTPCMonitorAltro::GetHwMaskFECChannel());
      fPadMapRCU[nextHwAddress][4]=  fChannelIter;
      fPadMapRCU[nextHwAddress][5]=  altro->GetTrailerPos();
              
      time_stamp    = 0;
      sample_iter   = 0;
      sample_length = 0;
      sample_bins   = 0;
      
      max           = 0.0;
      max_x         = 0.0;
      sum           = 0.0;
      sumn          = 0;
      
      hbase->Reset();
      
      for(Int_t iter_words = 0 ; iter_words< nwords ; iter_words++) 
	{
	  if(entries[block_pos-iter_words]==682) { continue; }
	  if(entries[block_pos-iter_words]!=682 && sample_iter==0) 
	    {
	      sample_length =  entries[block_pos-iter_words];
	      sample_iter   =  sample_length;
	      sample_bins   =  0;
	      time_stamp    =  entries[block_pos-iter_words-1];
	      iter_words++;
	    }
	  else 
	    {
	      ntime = time_stamp-sample_bins;
	      adc   = entries[block_pos-iter_words];
	      fPad[fChannelIter][ntime]  = adc;
	      if( (adc!=0)  && (ntime>=GetRangeBaseMin()  ) && (ntime<GetRangeBaseMax()    ))  {hbase->Fill(adc)        ;}
	      if( (adc>max) && (ntime>=GetRangeMaxAdcMin()) && (ntime<GetRangeMaxAdcMax()  ))  {max = adc;max_x = ntime ;}
	      if(              (ntime>=GetRangeSumMin())    && (ntime<GetRangeSumMax()     ))  {sum+=adc; sumn++;}
	      sample_bins++;
	    }
	}
      
      h_mean = hbase->GetMean();
      hbase->GetXaxis()->SetRangeUser(h_mean- h_mean/3 , h_mean + h_mean/3);
      h_mean =  hbase->GetMean();
      h_rms  = hbase->GetRMS();

      if(       GetPedestals()==1) fHistAddrMaxAdc->SetBinContent(  nextHwAddress,max- h_mean);
      else                         fHistAddrMaxAdc->SetBinContent(  nextHwAddress,max        );
      
      if(secid!=-1)
	{
	  if(rcu_patch<2)
	    {
	      fHistDistrBase2dIROC->Fill(h_mean,h_rms);
	      fHistDistrSumIROC->Fill(sum);
	      if(     GetPedestals()==1 ) { fHistDistrMaxIROC->Fill(max-h_mean);  fHistDistrSumIROC->Fill(sum -sumn*h_mean);}
	      else                        { fHistDistrMaxIROC->Fill(max);         fHistDistrSumIROC->Fill(sum );}
	    }
	  else
	    {
	      fHistDistrBase2dOROC->Fill(h_mean,h_rms);
	      fHistDistrSumOROC->Fill(sum);
	      if(     GetPedestals()==1 ){ fHistDistrMaxOROC->Fill(max-h_mean); fHistDistrSumOROC->Fill(sum -sumn*h_mean);}
	      else                       { fHistDistrMaxOROC->Fill(max);        fHistDistrSumOROC->Fill(sum)             ;}
	    }
	  
	  fHistAddrAdcSum->SetBinContent(  nextHwAddress,sum);
	  fHistAddrMapIndex->SetBinContent(nextHwAddress,fChannelIter);
	  fHistAddrBaseMean->SetBinContent(nextHwAddress,h_mean);
	  fHistAddrMaxAdcX->SetBinContent( nextHwAddress,max_x);
	  fHistAddrBaseRms->SetBinContent( nextHwAddress,h_rms);
	}
      fChannelIter++;
      if(nextpos<0)  { AliError("Error :  next pos < 0 "); break  ;}
    }
  delete hbase;
  return ;
}


//____________________________________________________________________________
void AliTPCMonitor::FillHistsPadPlane() 
{
  // Fill 2Dim histograms for IROC and OROC (max , rms and sum)
  
  if(fVerb)cout << "AliTPCMonitor::FillHistsPadPlane() Start  " << endl;
  if(fVerb)PrintConfig();

  Int_t pad    = 0;
  Int_t row    = 0;
  Int_t padmax = 0;
  Int_t hwadd  = 0;
  
  for(Int_t ch = 0; ch<fChannelIter; ch++) 
    {
      hwadd=  fPad[ch][0];
      fHistChannelTime->SetCellContent(ch,0,hwadd);
      
      for(Int_t bin = 1; bin <GetTimeBins(); bin++) 
	{ 
	  if( fHistChannelTime->GetCellContent(ch,bin)!=0)  cout << " cellcontent already set " << endl;
	  if(     GetPedestals()==1 ) fHistChannelTime->SetCellContent(ch,bin,(fPad[ch][bin]- fHistAddrBaseMean->GetBinContent(hwadd)));  
	  else                        fHistChannelTime->SetCellContent(ch,bin,(fPad[ch][bin]));
	}
      
      pad    = fMapHand->GetPad(   hwadd);
      row    = fMapHand->GetPadRow(hwadd);
      padmax = fMapHand->GetNumofPads(row);
      
      if(row<63)  
	{
	  fHistIROC->SetCellContent(     row    +1 ,pad +55 -padmax/2 +1,fHistAddrMaxAdc->GetBinContent(  hwadd));
	  fHistIROCIndex->SetCellContent(row    +1 ,pad +55 -padmax/2 +1,ch);
	  fHistIROCRMS->SetCellContent(  row    +1 ,pad +55 -padmax/2 +1,fHistAddrBaseRms->GetBinContent( hwadd));
	  fHistIROCBASE->SetCellContent( row    +1 ,pad +55 -padmax/2 +1,fHistAddrBaseMean->GetBinContent(hwadd));
	  fHistIROCSUM->SetCellContent(  row    +1 ,pad +55 -padmax/2 +1,fHistAddrAdcSum->GetBinContent(  hwadd));
	} 
      else 
	{ 
	  fHistOROC->SetCellContent(     row-63 +1 ,pad +70 -padmax/2 +1,fHistAddrMaxAdc->GetBinContent(  hwadd));
	  fHistOROCIndex->SetCellContent(row-63 +1 ,pad +70 -padmax/2 +1,ch);
	  fHistOROCRMS->SetCellContent(  row-63 +1 ,pad +70 -padmax/2 +1,fHistAddrBaseRms->GetBinContent( hwadd));
	  fHistOROCBASE->SetCellContent( row-63 +1 ,pad +70 -padmax/2 +1,fHistAddrBaseMean->GetBinContent(hwadd));
	  fHistOROCSUM->SetCellContent(  row-63 +1 ,pad +70 -padmax/2 +1,fHistAddrAdcSum->GetBinContent(  hwadd));
	}
    }
  
  fHistChannelTime->GetXaxis()->SetRange(0,fChannelIter);
  fHistChannelTime->GetYaxis()->SetRange(0,GetTimeBins());
}



//____________________________________________________________________________
void AliTPCMonitor::ResetArrays() 
{
  // Reset data arrays 
  for(Int_t row = 0 ; row < kNRowsIroc; row++) 
    {
    for(Int_t pad = 0 ; pad <  kNPadsIroc ; pad++) 
      {
	fHistIROCIndex->SetCellContent(row+1,pad+1,-1);
      }
    }
  for(Int_t row = 0 ; row < kNRowsOroc; row++) 
    {
    for(Int_t pad = 0 ; pad <  kNPadsOroc ; pad++) 
      {
	fHistOROCIndex->SetCellContent(row+1,pad+1,-1);
      }
    }
  
  for(Int_t ch= 0; ch<GetMaxHwAddr(); ch++) 
    { 
      fHistAddrMaxAdcX->SetBinContent(ch,-1);
      fHistAddrMapIndex->SetBinContent(ch,-1);
      fHistAddrMaxAdc->SetBinContent(  ch, 0);
      fHistAddrBaseMean->SetBinContent(  ch, 0);
      fHistAddrAdcSum->SetBinContent(  ch, 0);
      fHistAddrBaseRms->SetBinContent(ch, 0);
    }
  
  for(Int_t ch = 0; ch< GetNumOfChannels(); ch++)
    {
    for(Int_t bin = 0; bin< GetTimeBins(); bin++) 
      {
	fPad[ch][bin] = 0;
      }
    }
  for(Int_t ch = 0; ch< GetMaxHwAddr(); ch++) 
    {
      fPadMapHw[ch]=-1;
      fPadMapRCU[ch][0]=-1;
      fPadMapRCU[ch][1]=-1;
      fPadMapRCU[ch][2]=-1;
      fPadMapRCU[ch][3]=-1;
      fPadMapRCU[ch][4]=-1;
      fPadMapRCU[ch][5]=-1;
    }
  fPadUsedHwAddr=-1;
}


//____________________________________________________________________________
void AliTPCMonitor::ResetHistos() 
{
  // Reset all but 
  for(Int_t i =0; i<fHistList->GetEntries(); i++)
    {
      if(GetProcNextEvent()==0 && strcmp(((TH1*)fHistList->At(i))->GetName(),"SIDE A")==0) continue;
      if(GetProcNextEvent()==0 && strcmp(((TH1*)fHistList->At(i))->GetName(),"SIDE C")==0) continue;
      ((TH1*)fHistList->At(i))->Reset();
    }
  ResetArrays();
}

//____________________________________________________________________________
void AliTPCMonitor::DeleteHistos() 
{
  // Delete histograms 
  for(Int_t i =0; i<fHistList->GetEntries(); i++)
    {
      delete (TH1*)fHistList->At(i);
    }
}


//__________________________________________________________________
Int_t AliTPCMonitor::CheckEqId(Int_t secid,Int_t eqid)
{
  // Check if equipment id corresponds to any  rcu patch in sector 
  // Equipment ids changed during commisioning in 2006 (starting from run 704)
  // However Runids started from 0 again in 2007
  // Now only runids from commissioning in 2006 after runid 704 and all new ones are supported.
  // Comment in equipment check for runids < 704 if old runs should be processed

  if(fVerb) cout << "AliTPCMonitor::CheckEqId  : SectorId  " << secid << " EquipmentId " << eqid  << " runid " << fRunId <<  endl;
  Int_t passed =1;
  if(fRunId<704 && 0) // commented out --> runs with runid < 704 in 2006 are not recognized anymore
    {
      if( (secid>-1) && (secid<36) )   // if ( secid is in range) { take only specific eqids}  else { take all }
	{
	  if(      (secid==13) && ( eqid!=408 && eqid!=409  &&  eqid!=509 && eqid!=512  && eqid!=513 && eqid!=517 )) {passed=0;} 
	  else if( (secid==4)  && ( eqid!=404 && eqid!=504  &&  eqid!=407 && eqid!=503  && eqid!=508 && eqid!=506 )) {passed=0;} 
	}
      else                                                                   {if(fVerb) cout << "passed check "<< endl; }
    }
  else
    {
      if( (secid>-1) && (secid<36) )   // if ( secid is in range) { take only specific eqids}  else { take all }
	{
	  if(eqid!=fMapEqidsSec[secid][0] && eqid!= fMapEqidsSec[secid][1] && eqid!=fMapEqidsSec[secid][2] &&
	     eqid!=fMapEqidsSec[secid][3] && eqid!= fMapEqidsSec[secid][4] && eqid!=fMapEqidsSec[secid][5] )  {passed=0;} 
	}
      else                                                                   {if(fVerb) cout << "passed check "<< endl;}
    }
  return passed;
}

//__________________________________________________________________
void AliTPCMonitor::SetEqIds()
{
  // Set mapping for equipment ids
  for(Int_t i = 0; i<36 ; i++) 
    {     
      for(Int_t j = 0; j<6; j++)
	{
	  if(j<2) fMapEqidsSec[i][j]= 768+i*2+j;
	  else    fMapEqidsSec[i][j]= 840+i*4+j-2;
	}
    }
  
  for(Int_t i = 0; i<36 ; i++)
    {     
      for(Int_t j = 0; j<6; j++)
	{
	  if(j<2) fMapEqidsRcu[768+i*2+j]   = i*6 +j;
	  else    fMapEqidsRcu[840+i*4+j-2] = i*6 +j;
	}
    }
}

//__________________________________________________________________
void AliTPCMonitor::FillGlobal(Int_t sector) 
{

  // Fill global histograms with max adc for each channel 
  
  TH2S* hglob =0;
  if((sector/18) ==0) hglob =  fHistGlobalMaxA;
  else                hglob =  fHistGlobalMaxC;
  
  Float_t  rotsec  = (2*TMath::Pi()/18.0);
  Float_t  rot     = (-rotsec*(sector%18)  +4*rotsec);
  
  Float_t  m11     =    TMath::Cos(rot);
  Float_t  m12     =    TMath::Sin(rot);
  Float_t  m21     = -1*TMath::Sin(rot);
  Float_t  m22     =    TMath::Cos(rot);
  
  Int_t    max     = 0; // use integer for global view
  
  Double_t xval    = 0.0;
  Double_t yval    = 0.0;
  
  Int_t    pad     = 0;
  Int_t    row     = 0;
  Int_t    padmax  = 0;
  
  Float_t  xdr     = 0;
  Float_t  ydr     = 0;
  
  for(Int_t hw = 0; hw<fHistAddrMaxAdc->GetNbinsX(); hw++) 
    {
      max = (Int_t)fHistAddrMaxAdc->GetBinContent(hw);
      if(max!=-1) 
	{
	  pad    = fMapHand->GetPad(   hw);
	  row    = fMapHand->GetPadRow(hw);
	  padmax = fMapHand->GetNumofPads(row);
	  GetXY(xval ,yval , padmax,row ,pad);
	  xdr    =  xval*m11 +yval*m12;
	  ydr    =  xval*m21 +yval*m22;
	  if(hglob->GetBinContent(hglob->GetXaxis()->FindBin(xdr),hglob->GetYaxis()->FindBin(ydr))==0)   hglob->Fill(xdr,ydr,(Int_t)max);
	}
    }
}


//__________________________________________________________________
void AliTPCMonitor::GetXY( Double_t& xval , Double_t& yval , Int_t padmax, Int_t row , Int_t pad) 
{
  // Get x and y position of pad
  
  if(row<63) 
    {
      xval = fMirror*( 2*padmax -4*pad -2);
      yval = 852.25 +7.5*row;
    } 
  else 
    {
      xval = fMirror*( 3*padmax -6*pad -3);
      if((row-63)<63) {	  yval = 10*(row-63) +1351;	} 
      else 	      {	  yval = 15*(row-63-64)+1993.5;	}
    }
  
} 

//__________________________________________________________________
Int_t AliTPCMonitor::GetPadAtX(Float_t xval, Int_t row, Int_t padmax)
{
  // Get pad number at given position in x
  Int_t pad    = 0;
  
  if(row<63) {pad = (Int_t)( ( (xval/fMirror) +2 -2*padmax)/-4);}
  else       {pad = (Int_t)( ( (xval/fMirror) +3 -3*padmax)/-6);}
  
  if(pad>=padmax) return -1;
  else            return pad ;
  
}

//__________________________________________________________________
Int_t AliTPCMonitor::GetPadAtX(Float_t xval, Int_t row)
{

 // Get pad number at given position in x

  Int_t padmax = fMapHand->GetNumofPads(row);
  Int_t pad    = 0;
  
  if(row<63) {pad = (Int_t)( ( (xval/fMirror) +2 -2*padmax)/-4);}
  else       {pad = (Int_t)( ( (xval/fMirror) +3 -3*padmax)/-6);}
  
  if(pad>=padmax) return -1;
  else            return pad ;
  
}

//__________________________________________________________________
void AliTPCMonitor::DrawHists(Int_t histos) 
{

  // Draw sets of histograms
  // histos==1 : 2Dim histos for MAX adc  and add executables
  // histos==2 : distributions max/rms/sum 
  // histos==3 : global max adc for specified SideA/C


  if(fVerb)    cout << " Draw histos " << endl;  
  Char_t cside[10];
  if(GetLastSector()/18==0 ) sprintf(cside,"A");
  else                       sprintf(cside,"C");
  
  Char_t titleSEC[256];   sprintf(titleSEC   ,"Sector %i Side %s Run : %05i EventID %i "       ,GetLastSector()%18,cside,fRunId, fEventNumber);
  Char_t titleEvent[256]; sprintf(titleEvent ,"Time <-> Channles  %s"                          ,titleSEC);
  Char_t titleIROC[256];  sprintf(titleIROC  ,"IROC %s"                                        ,titleSEC);
  Char_t titleOROC[256];  sprintf(titleOROC  ,"OROC %s"                                        ,titleSEC);
  
  Char_t titleMAX[256];   sprintf(titleMAX   ,"Max (timebin: %i,%i) %s"                        ,GetRangeMaxAdcMin(),GetRangeMaxAdcMax(),titleSEC);
  Char_t titleSUM[256];   sprintf(titleSUM   ,"Sum (timebin: %i,%i) %s"                        ,GetRangeSumMin()   ,GetRangeSumMax()   ,titleSEC);
  Char_t titleBASE[256];  sprintf(titleBASE  ,"Baseline RMS<->Mean  (timebin: %i-%i) %s"       ,GetRangeBaseMin()  ,GetRangeBaseMax()  ,titleSEC);
  Char_t titleMEAN[256];  sprintf(titleMEAN  ,"Baseline Mean (timebin: %i-%i) %s"              ,GetRangeBaseMin()  ,GetRangeBaseMax()  ,titleSEC);
  Char_t titleRMS[256] ;  sprintf(titleRMS   ,"Baseline RMS (timebin: %i-%i) %s"               ,GetRangeBaseMin()  ,GetRangeBaseMax()  ,titleSEC);

  if(histos==1) 
    {
      // IROC _______________________________________________________________
      TCanvas* ciroc = 0;
      ciroc = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("ciroc");
      if(!ciroc) 
	{
	  ciroc = CreateCanvas("ciroc");
	  fExecPlaneMax=0;
	}
      ciroc->cd();
      
      fHistIROC->SetXTitle("row");
      fHistIROC->SetYTitle("pad");
      if(GetPedestals()) fHistIROC->SetZTitle("max ADC (baseline sub)");
      else               fHistIROC->SetZTitle("max ADC ");
      fHistIROC->SetTitle(titleIROC);
      fHistIROC->SetMinimum(0.01);
      fHistIROC->Draw("COLZ");
      ciroc->UseCurrentStyle();
      
      
      fHistIROCTime->SetXTitle("row"); fHistIROCTime->SetZTitle("peak time (fit)");       fHistIROCTime->SetYTitle("pad");      fHistIROCTime->SetTitle(titleIROC);
      fHistIROCRMS->SetXTitle("row");  fHistIROCRMS->SetZTitle( "baseline rms (ADC)");    fHistIROCRMS->SetYTitle("pad");       fHistIROCRMS->SetTitle(titleIROC);
      
      // OROC
      TCanvas* coroc = 0;
      coroc =(TCanvas*)gROOT->GetListOfCanvases()->FindObject("coroc");
      if(!coroc) {
	coroc = CreateCanvas("coroc");
	fExecPlaneMax=0;
      }
      coroc->cd();
      
      fHistOROC->SetXTitle("row");
      fHistOROC->SetYTitle("pad");
      if(GetPedestals()) fHistOROC->SetZTitle("max ADC (baseline sub)");
      else               fHistOROC->SetZTitle("max ADC ");
      fHistOROC->SetTitle(titleOROC);
      fHistOROC->SetMinimum(0.01);
      fHistOROC->Draw("COLZ");
      coroc->UseCurrentStyle();
      

      fHistOROCTime->SetXTitle("row"); fHistOROCTime->SetZTitle("peak time (fit) (timebins)"); fHistOROCTime->SetYTitle("pad"); fHistOROCTime->SetTitle(titleOROC);
      fHistOROCRMS->SetXTitle("row");  fHistOROCRMS->SetZTitle("baseline rms (ADC)");          fHistOROCRMS->SetYTitle("pad");  fHistOROCRMS->SetTitle(titleOROC);

      // SUM 
      Char_t name_sum[256] ; sprintf(name_sum,"ADC sum (bins: %i, %i)",GetRangeSumMin() ,GetRangeSumMax() );
      fHistIROCSUM->SetXTitle("row");      fHistIROCSUM->SetZTitle(name_sum);      fHistIROCSUM->SetYTitle("pad");      fHistIROCSUM->SetTitle(titleIROC);
      fHistOROCSUM->SetXTitle("row");      fHistOROCSUM->SetZTitle(name_sum);      fHistOROCSUM->SetYTitle("pad");      fHistOROCSUM->SetTitle(titleOROC);
    
      // BASE
      Char_t name_base[256] ; sprintf(name_base ,"base mean (timbebin: %i, %i )",GetRangeBaseMin(),GetRangeBaseMax());
      fHistIROCBASE->SetXTitle("row"); fHistIROCBASE->SetZTitle(name_base);  fHistIROCBASE->SetYTitle("pad");      fHistIROCBASE->SetTitle(titleIROC);
      fHistOROCBASE->SetXTitle("row"); fHistOROCBASE->SetZTitle(name_base);  fHistOROCBASE->SetYTitle("pad");      fHistOROCBASE->SetTitle(titleOROC);

      if(fHistIROCClone) fHistIROCClone->Delete();
      if(fHistOROCClone) fHistOROCClone->Delete();
      fHistIROCClone = (TH2F*)fHistIROC->Clone("fHistIROCClone");
      fHistOROCClone = (TH2F*)fHistOROC->Clone("fHistOROCClone");
      
      // Executables
      if(fExecPlaneMax==0) 
	{
	  Char_t carry1[100];
	  sprintf(carry1,".x %s/TPC/AliTPCMonitorExec.C(1)",gSystem->Getenv("ALICE_ROOT"));
	  ciroc->AddExec("pad",carry1);
	  coroc->AddExec("pad",carry1);
	  
	  Char_t carry2[100];
	  sprintf(carry2,".x %s/TPC/AliTPCMonitorExec.C(2)",gSystem->Getenv("ALICE_ROOT"));
	  ciroc->AddExec("row",carry2);
	  coroc->AddExec("row",carry2);
	  fExecPlaneMax=1;
	}
      coroc->Update();
      ciroc->Update();
    }
  else if(histos==2) 
    {
      // MAX ADC distribution  ____________________________________________
      TCanvas* cmax = 0;
      cmax = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cmax");
      if(!cmax)  cmax = CreateCanvas("cmax"); 
    
      cmax->cd();
      fHistDistrMaxIROC->GetXaxis()->SetRangeUser(0.0,1000.0);
      fHistDistrMaxIROC->SetXTitle("max ADC (ADC)");
      fHistDistrMaxIROC->SetYTitle("counts");
      fHistDistrMaxIROC->SetTitle(titleMAX);
      fHistDistrMaxIROC->Draw("");
      fHistDistrMaxOROC->SetLineColor(2);
      fHistDistrMaxOROC->Draw("same");
    
      if(fHistDistrMaxOROC->GetMaximum()> fHistDistrMaxIROC->GetMaximum())  fHistDistrMaxIROC->SetMaximum(fHistDistrMaxOROC->GetMaximum()*1.1);
      
      TLegend* legio = new TLegend(0.6,0.6,0.8,0.8);
      legio->SetFillColor(0);
      legio->AddEntry(fHistDistrMaxIROC,"IROC","l");
      legio->AddEntry(fHistDistrMaxOROC,"OROC","l");
      legio->Draw("same");
    
      // ADC sum distribution
      TCanvas* csum = 0;
      csum = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("csum");
      if(!csum)  csum = CreateCanvas("csum") ;
      csum->cd();
    
      fHistDistrSumIROC->SetXTitle("sum ADC (ADC)");
      fHistDistrSumIROC->SetYTitle("counts");
      fHistDistrSumIROC->SetTitle(titleSUM);
      fHistDistrSumIROC->Draw("");
      fHistDistrSumOROC->SetLineColor(2);
      fHistDistrSumOROC->Draw("same");
      if(fHistDistrSumOROC->GetMaximum()> fHistDistrSumIROC->GetMaximum())  fHistDistrSumIROC->SetMaximum(fHistDistrSumOROC->GetMaximum()*1.1);
      legio->Draw("same");
    
      // BASELINE MEAN distribution
      TCanvas* cbasemean = 0;
      cbasemean = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cbasemean");
      if(!cbasemean)  cbasemean = CreateCanvas("cbasemean"); 
      cbasemean->cd();
    
      fHistDistrBaseMeanIROC = fHistDistrBase2dIROC->ProjectionX("fHistDistrBaseMeanIROC");
      fHistDistrBaseMeanIROC->SetXTitle("base mean (ADC)");
      fHistDistrBaseMeanIROC->SetYTitle("counts");
      fHistDistrBaseMeanIROC->SetTitle(titleMEAN);
      fHistDistrBaseMeanIROC->Draw("");
    
      fHistDistrBaseMeanOROC = fHistDistrBase2dOROC->ProjectionX("fHistDistrBaseMeanOROC");
      fHistDistrBaseMeanOROC->SetLineColor(2);
      fHistDistrBaseMeanOROC->Draw("same");
      if(fHistDistrBaseMeanOROC->GetMaximum()>fHistDistrBaseMeanIROC->GetMaximum())   fHistDistrBaseMeanIROC->SetMaximum(fHistDistrBaseMeanOROC->GetMaximum()*1.1);
      legio->Draw("same");

      TCanvas* cbaserms = 0;
      cbaserms = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cbaserms");
      if(!cbaserms) cbaserms = CreateCanvas("cbaserms") ;
      cbaserms->cd();
    
      // BASELINE RMS distribution
      fHistDistrBaseRmsIROC = fHistDistrBase2dIROC->ProjectionY("fHistDistrBaseRmsIROC");
      fHistDistrBaseRmsIROC->SetXTitle("base rms (ADC)");
      fHistDistrBaseRmsIROC->SetYTitle("counts");
      fHistDistrBaseRmsIROC->SetTitle(titleRMS);
      fHistDistrBaseRmsIROC->Draw("");
    
      fHistDistrBaseRmsOROC = fHistDistrBase2dOROC->ProjectionY("fHistDistrBaseRmsOROC");
      fHistDistrBaseRmsOROC->SetLineColor(2);
      fHistDistrBaseRmsOROC->Draw("same");
      if(fHistDistrBaseRmsOROC->GetMaximum()>fHistDistrBaseRmsIROC->GetMaximum())  fHistDistrBaseRmsIROC->SetMaximum(fHistDistrBaseRmsOROC->GetMaximum()*1.1);
      legio->Draw("same");
    
      cmax->Update();
      csum->Update();
      cbasemean->Update();
      cbaserms->Update();
    }
  else  if(histos==3)
    {
      // GLOBAL MAX ADC _________________________________
      if(GetProcNextEvent()==1)
	{
	  TCanvas* cglobA =0;
	  TCanvas* cglobC =0;
	    
	  if(!(cglobC=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("SIDE C all"))) cglobC = CreateCanvas("SIDE C all"); 
	  if(!(cglobA=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("SIDE A all"))) cglobA = CreateCanvas("SIDE A all"); 
	    
	  Char_t gl_title1[256]; sprintf(gl_title1,"SIDE A Run %05i (EventID %i)",fRunId,fEventNumber);
	  Char_t gl_title2[256]; sprintf(gl_title2,"SIDE C Run %05i (EventID %i)",fRunId,fEventNumber);
	    
	  fHistGlobalMaxA->SetTitle(gl_title1);
	  fHistGlobalMaxC->SetTitle(gl_title2);
	  fHistGlobalMaxA->SetXTitle("x/mm");
	  fHistGlobalMaxA->SetYTitle("y/mm");
	  fHistGlobalMaxC->SetXTitle("x/mm");
	  fHistGlobalMaxC->SetYTitle("y/mm");
	    
	  if(GetPedestals()==0)     {	      fHistGlobalMaxA->SetZTitle("max adc (not baseline sub)");   fHistGlobalMaxC->SetZTitle("max adc (not baseline sub)");  }
	  else 	                    {	      fHistGlobalMaxA->SetZTitle("max adc ");	                  fHistGlobalMaxC->SetZTitle("max adc ");	                   }
	  
	  fHistGlobalMaxA->SetMinimum(0.01);
	  fHistGlobalMaxC->SetMinimum(0.01);
	    
	  cglobC->cd() ; fHistGlobalMaxC->Draw("COLZ");
	  cglobA->cd() ; fHistGlobalMaxA->Draw("COLZ");
	    
	  Char_t name_om[256];
	  sprintf(name_om,".x  %s/TPC/AliTPCMonitorExec.C(3)",gSystem->Getenv("ALICE_ROOT"));
	    
	  if(fExecGlob==0) 
	    { 
	      if(fVerb)cout << " set exec " << name_om << endl;
	      cglobC->AddExec("glob",name_om);
	      cglobA->AddExec("glob",name_om);
	      fExecGlob = 1;
	    } 
	  else 
	    {
	      cglobC->DeleteExec("glob");
	      cglobA->DeleteExec("glob");
		
	      if(fVerb)  cout << " set exec " << name_om << endl;
	      cglobC->AddExec("glob",name_om);
	      cglobA->AddExec("glob",name_om);
		
	    }
	  cglobC->Update();
	  cglobA->Update();
	}
      
    }
}



//__________________________________________________________________
void AliTPCMonitor::DrawRMSMap() 
{
  // Draw 2Dim rms histos for IROC and OROC 
  // and set executables for canvases
  
  TCanvas* crmsoroc =0;
  TCanvas* crmsiroc =0;
  if(!(crmsoroc = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("crmsoroc")))    crmsoroc    = CreateCanvas("crmsoroc");
  if(!(crmsiroc = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("crmsiroc")))    crmsiroc    = CreateCanvas("crmsiroc");
   
  crmsiroc->cd();  fHistIROCRMS->Draw("COLZ");
  crmsoroc->cd();  fHistOROCRMS->Draw("COLZ");
 
  Char_t carry1[100];  sprintf(carry1,".x %s/TPC/AliTPCMonitorExec.C(1)",gSystem->Getenv("ALICE_ROOT"));
  Char_t carry2[100];  sprintf(carry2,".x %s/TPC/AliTPCMonitorExec.C(2)",gSystem->Getenv("ALICE_ROOT"));
 
  if(fExecPadIrocRms==0) 
    {
      crmsiroc->AddExec("pad",carry1);
      crmsiroc->AddExec("row",carry2);
      fExecPadIrocRms=1;
    }
  
  if(fExecPadOrocRms==0) 
    {
      crmsoroc->AddExec("pad",carry1);
      crmsoroc->AddExec("row",carry2);
      fExecPadOrocRms=1;
    }
  
  crmsiroc->Update();
  crmsoroc->Update();
  
  DrawHists(2); 

}

//__________________________________________________________________
void AliTPCMonitor::ExecPad() 
{

  // Executable for Pad 
  // Show time profile for channel the mouse is pointing at 
  
  Int_t event = gPad->GetEvent();
  if (event != 51)   return;
  
  TObject *select = gPad->GetSelected();
  if(!select)    return;
  if(!select->InheritsFrom("TH2")) { return;  }
  gPad->GetCanvas()->FeedbackMode(kTRUE);
  
  // get position
  Int_t    px        = gPad->GetEventX();
  Int_t    py        = gPad->GetEventY();
  Float_t  upy       = gPad->AbsPixeltoY(py);
  Float_t  upx       = gPad->AbsPixeltoX(px);
  Float_t  y         = gPad->PadtoY(upy);
  Float_t  x         = gPad->PadtoX(upx);

  Int_t    setrange  = 0;
  
  TCanvas* cpad      = 0;
  //  Char_t   namehist[50];
  Char_t   projhist[60];
  Char_t   namesel[256];
  Char_t   namecanv[256];

  Int_t    xbin_min  = 0;
  Int_t    xbin_max  = 0;
  Float_t  ybin_min  = 0;
  Float_t  ybin_max  = 0;
  Int_t    ROC       = 0;
  
  // Check wich Canvas executed the event 
  TH2S* fHistIndex=0;
  sprintf(namesel,select->GetName());
  if(strcmp(namesel,"fHistOROC")==0 || strcmp(namesel,"fHistOROCRMS")==0 || strcmp(namesel,"fHistOROCTime")==0 ) 
    {
      ROC = 1;
      fPadUsedRoc =1;
      sprintf(projhist,"ProjectionOROC");
      sprintf(namecanv,"coroc_ch");
      fHistIndex = fHistOROCIndex;
    }
  if(strcmp(namesel,"fHistIROC")==0 || strcmp(namesel,"fHistIROCRMS")==0 || strcmp(namesel,"fHistIROCTime")==0 ) 
    {
      ROC = 0;
      fPadUsedRoc=0;
      sprintf(projhist,"ProjectionIROC");
      sprintf(namecanv,"ciroc_ch");
      fHistIndex = fHistIROCIndex;
    }
  
  // Check if Canvas already existed and get Ranges from former Prjection histogram 
  if(cpad=(TCanvas*)gROOT->GetListOfCanvases()->FindObject(namecanv)) 
    {
      cpad->cd();
      if(gROOT->Get(projhist)) 
	{
	  setrange = 1;
	  xbin_min = ((TH1D*)gROOT->Get(projhist))->GetXaxis()->GetFirst();
	  xbin_max = ((TH1D*)gROOT->Get(projhist))->GetXaxis()->GetLast();
	  ybin_min = ((TH1D*)gROOT->Get(projhist))->GetMinimum();
	  ybin_max = ((TH1D*)gROOT->Get(projhist))->GetMaximum();
	  delete gROOT->Get("legfit");
	  delete gROOT->Get("fg");
	}
    } 
  else 
    {
      cpad = CreateCanvas(namecanv); cpad->cd();
    }
  
  // Get Bin 
  Int_t testy = fHistIndex->GetYaxis()->FindBin(y);
  Int_t testx = fHistIndex->GetXaxis()->FindBin(x);
  Int_t binchannel = (Int_t)fHistIndex->GetCellContent(testx,testy);
  if(binchannel>30000 || binchannel<0) return;

  if(gROOT->Get(projhist))  delete gROOT->Get(projhist);
  // Get Projection   
  TH1D *hp = (TH1D*)(((TH1D*)fHistChannelTime->ProjectionY("hp",binchannel,binchannel))->Clone(projhist));
  
  // Make title and Pave for channel Info
  Char_t title[256];
  Int_t npadRow , npad  , nhw , nmax , hwadd;
  
  hwadd   = (Int_t)fHistChannelTime->GetCellContent(binchannel,0);
  fPadUsedHwAddr = hwadd;
  
  if(ROC==0) npadRow = fMapHand->GetPadRow(hwadd);
  else       npadRow = fMapHand->GetPadRow(hwadd)-63;
  npad               = fMapHand->GetPad(hwadd);
  nhw                = hwadd;
  nmax               = (Int_t)hp->GetMaximum();
  

  TPaveText*  legstat = new TPaveText(0.18,0.65,0.3,0.8,"NDC");

  Int_t   connector   =  fMapHand->GetFECconnector(hwadd);
  Int_t   FECnr       =  fMapHand->GetFECfromHw(hwadd);
  Int_t   FECch       =  fMapHand->GetFECchannel(hwadd);
  Int_t   Altro       =  fMapHand->GetAltro(hwadd);
  Int_t   AltroCh     =  (fMapHand->GetAltroChannel(hwadd))%16;
  Int_t   FECloc      =  fMapHand->U2fGetFECinRCU(FECnr) ;
  Int_t   FEClocbran  =  fMapHand->U2fGetFECinBranch(FECnr);
  Int_t   branch      =  fMapHand->U2fGetBranch(FECnr);

  
  Short_t fecget      = (hwadd & AliTPCMonitorAltro::GetHwMaskFEC())   >> 7;
  Short_t branchget   = (hwadd & AliTPCMonitorAltro::GetHwMaskBranch())>> 11;


  Char_t nstat1[100];  Char_t nstat2[100];  Char_t nstat3[100];  Char_t nstat4[100];  
  Char_t nstat5[100];  Char_t nstat6[100];  Char_t nstat7[100];  Char_t nstat8[100];
  
  sprintf(nstat1,"Branch (map) \t %i (%i) \n",branchget,branch);
  sprintf(nstat2,"Fec in patch \t %i \n",FECloc);
  sprintf(nstat8,"Fec in branch (map)\t %i (%i)\n",fecget,FEClocbran);
  sprintf(nstat7,"Connector  \t %i \n",connector);
  sprintf(nstat3,"Fec No.   \t %i \n",FECnr);
  sprintf(nstat4,"Fec chan  \t %i \n",FECch);
  sprintf(nstat5,"Altro chip\t %i \n",Altro);
  sprintf(nstat6,"Altro chan\t %i \n",AltroCh);
  
  legstat->AddText(nstat1); legstat->AddText(nstat2);  legstat->AddText(nstat8);  legstat->AddText(nstat7);  
  legstat->AddText(nstat3); legstat->AddText(nstat4);  legstat->AddText(nstat5);  legstat->AddText(nstat6);

  sprintf(title,"Row=%d Pad=%d Hw =%d maxADC =%d count =%d",npadRow,npad,nhw,nmax,binchannel);
  
  hp->SetName(projhist);
  hp->SetTitleSize(0.04);
  hp->SetTitle(title);
  hp->SetYTitle("ADC");
  hp->SetXTitle("Timebin");
  hp->GetXaxis()->SetTitleColor(1);
  
  if(setrange) 
    {
      hp->GetXaxis()->SetRange(xbin_min,xbin_max);
      hp->SetMinimum(ybin_min);
      hp->SetMaximum(ybin_max);
    }
  else  
    { 
      hp->SetMinimum(0.0);
      hp->SetMaximum(1000.0);
    }
  
  cpad->cd();
  hp->Draw();
  
  // Make Fit to peak
  if(GetPedestals() && fDisableFit==0) 
    {
      Int_t max_x  =    (Int_t)fHistAddrMaxAdcX->GetBinContent(hwadd);
      Float_t max  =  (Float_t)fHistAddrMaxAdc->GetBinContent(hwadd);
      Float_t base =  (Float_t)fHistAddrBaseMean->GetBinContent(hwadd);
      if(base!=0) 
	{
	  if( ((max+base)/base)>1.2) 
	    {
	      TF1* fg = new TF1("fg",AliTPCMonitor::Gamma4,max_x-5,max_x+5,4);
	      fg->SetParName(0,"Normalisation");
	      fg->SetParName(1,"Minimum");
	      fg->SetParName(2,"Width");
	      fg->SetParName(3,"Base");
	      fg->SetParameter(0,max);
	      fg->SetParameter(1,max_x-2);
	      fg->SetParameter(2,1.5);
	      fg->FixParameter(3,0);
	      fg->SetLineColor(4);
	      fg->SetLineWidth(1);
	      hp->Fit("fg","RQ");
	      
	      TLegend* legfit = new TLegend(0.6,0.7,0.7,0.8);
	      legfit->AddEntry("fg","#Gamma 4 fit","l");
	      legfit->SetFillColor(0);
	      legfit->SetName("legfit");
	      legfit->Draw("same");
	    }
	}
    }
  legstat->SetFillColor(0);
  legstat->Draw("same");
  cpad->Update();
  return;
}

//__________________________________________________________________
void AliTPCMonitor::ExecRow() 
{

  // Executable for Pad 
  // Show profile of max adc over given pad row 
  // and 2dim histo  adc(pad-in-row,time bin)  

  Int_t event = gPad->GetEvent();
  if (event != 61)    return;
  gPad->cd();
  TObject *select = gPad->GetSelected();
  if(!select)    return;
  if(!select->InheritsFrom("TH2")) {   return;  }

  Int_t ROC = 0;
  //  Char_t namehist[50];
  Char_t rowhist[60];
  Char_t rowhist_sum[60];
  Char_t rowhist_max[60];
  Char_t rowhist_xmax[60];
  
  sprintf(rowhist,         "hrow_time");
  sprintf(rowhist_xmax    ,"hx_max");
  sprintf(rowhist_max    , "hrow_max");
  
  // get position 
  Int_t   px  = gPad->GetEventX();
  Int_t   py  = gPad->GetEventY();
  Float_t upy = gPad->AbsPixeltoY(py);
  Float_t upx = gPad->AbsPixeltoX(px);
  Float_t y   = gPad->PadtoY(upy);
  Float_t x   = gPad->PadtoX(upx);
  
  TCanvas*crow_time   = 0;
  TCanvas*crow_max    = 0;
  TCanvas*cx_max      = 0;
  
  TH2S*  fHistIndex  = 0;

  // ranges from already existing histos 
  Int_t    row_time_xmin = 0;
  Int_t    row_time_xmax = 0;
  Int_t    row_time_ymin = 0;
  Int_t    row_time_ymax = 0;
  Float_t  row_time_zmin = 0;
  Float_t  row_time_zmax = 0;
  
  Int_t    prof_row_xmin = 0;
  Int_t    prof_row_xmax = 0;
  Double_t prof_row_ymin = 0;
  Double_t prof_row_ymax = 0;

  Int_t    prof_x_xmin   = 0;
  Int_t    prof_x_xmax   = 0;
  Double_t prof_x_ymin   = 0;
  Double_t prof_x_ymax   = 0;


  Int_t    setrange      = 0;


  if(     strcmp(select->GetName(),"fHistIROC")==0 || strcmp(select->GetName(),"fHistIROCRMS")==0 ) { fHistIndex = fHistIROCIndex;     ROC =1;   } 
  else if(strcmp(select->GetName(),"fHistOROC")==0 || strcmp(select->GetName(),"fHistOROCRMS")==0 ) { fHistIndex = fHistOROCIndex;     ROC =2;   }
  else                                                                                  { cout << " not implemented for this histo " << endl; return; }
  
  gPad->GetCanvas()->FeedbackMode(kTRUE);
  
 

  // check if canvases exist //
  crow_time  = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("crow_time");
  crow_max   = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("crow_max");
  cx_max     = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cx_max");
  
  if(!crow_time)   crow_time   = CreateCanvas("crow_time") ;
  if(!crow_max)    crow_max    = CreateCanvas("crow_max")  ;
  if(!cx_max  )    cx_max      = CreateCanvas("cx_max")    ;
  
  // check ranges of already existing histos 
  if(gROOT->Get(rowhist)) 
    {
      row_time_xmin  = ((TH2F*)gROOT->Get(rowhist))->GetXaxis()->GetFirst();
      row_time_xmax  = ((TH2F*)gROOT->Get(rowhist))->GetXaxis()->GetLast();
      row_time_ymin  = ((TH2F*)gROOT->Get(rowhist))->GetYaxis()->GetFirst();
      row_time_ymax  = ((TH2F*)gROOT->Get(rowhist))->GetYaxis()->GetLast();
      row_time_zmin  = ((TH2F*)gROOT->Get(rowhist))->GetMinimum();
      row_time_zmax  = ((TH2F*)gROOT->Get(rowhist))->GetMaximum();
      
      prof_row_xmin  = ((TH1F*)gROOT->Get(rowhist_max))->GetXaxis()->GetFirst();
      prof_row_xmax  = ((TH1F*)gROOT->Get(rowhist_max))->GetXaxis()->GetLast();
      prof_row_ymin  = ((TH1F*)gROOT->Get(rowhist_max))->GetMinimum();
      prof_row_ymax  = ((TH1F*)gROOT->Get(rowhist_max))->GetMaximum();

      prof_x_xmin    = ((TH1F*)gROOT->Get(rowhist_xmax))->GetXaxis()->GetFirst();
      prof_x_xmax    = ((TH1F*)gROOT->Get(rowhist_xmax))->GetXaxis()->GetLast();
      prof_x_ymin    = ((TH1F*)gROOT->Get(rowhist_xmax))->GetMinimum();
      prof_x_ymax    = ((TH1F*)gROOT->Get(rowhist_xmax))->GetMaximum();

      setrange =1;
      
      delete gROOT->Get(rowhist);
      delete gROOT->Get(rowhist_max);
      delete gROOT->Get(rowhist_sum);
      delete gROOT->Get("hx_max");
      delete gROOT->Get("legrow");
    }
  
  // get channel for xy bin -> getRow -> getNrows -> getHw for each Pad in Row -> get channel for each hw -> make proj
  Int_t testy      = fHistIndex->GetYaxis()->FindBin(y);
  Int_t testx      = fHistIndex->GetXaxis()->FindBin(x);
  Int_t binchannel = (Int_t)fHistIndex->GetCellContent(testx,testy);
  
  if(binchannel>30000)    return;
  if(binchannel<=0 ) { crow_time->Update() ;    crow_max->Update() ;    return ;  }
  
  // get hwaddress 
  Int_t hwadd     = (Int_t)fHistChannelTime->GetCellContent(binchannel,0);
  Int_t row       = fMapHand->GetPadRow(hwadd);
  Int_t pad       = fMapHand->GetPad(hwadd)   ;
  Int_t numofpads =  fMapHand->GetNumofPads(row);
  
  // create histos 
  TH2F *hrow_time = new TH2F(rowhist     , ""  ,numofpads,0,numofpads,GetTimeBins(),0.0,GetTimeBins());
  TH1F *hrow_max  = new TH1F(rowhist_max , ""  ,numofpads,0,numofpads);
  TH1F *hx_max    = new TH1F(rowhist_xmax, ""  ,159,0,159      );
  
  // Row profile ///////////
  if(fVerb) cout << " Number of pads " << numofpads << endl;
  for(Int_t padnr = 0; padnr<numofpads;padnr++) 
    {
      Int_t addrinrow = fMapHand->GetPadAddInRow(row,padnr );
      Int_t channel   = (Int_t)fHistAddrMapIndex->GetBinContent(addrinrow);
      if(channel==-1) continue;
      
      hrow_max->SetBinContent(padnr+1,fHistAddrMaxAdc->GetBinContent(addrinrow));
      TH1D *hp = fHistChannelTime->ProjectionY("hp",channel,channel);
      for(Int_t time = 0;time<GetTimeBins();time++) {
	
	Float_t val = hp->GetBinContent(time);
	hrow_time->SetCellContent(padnr+1,time+1,val);
      }
    }
  
  // X profile  /////////////
  Double_t xval  = 0.0;
  Double_t yval  = 0.0;
  GetXY(xval,yval,numofpads,row,pad);
  
  Int_t padnr = 0; 
  Int_t hw    = 0;
  for(Int_t nrow = 0; nrow<159; nrow++)
    {
      padnr = GetPadAtX(xval,nrow);
      if(padnr>=0)
	{
	  hw = fMapHand->GetPadAddInRow(nrow,padnr);
	  if(fPadMapHw[hw]==-1){ continue                                                      ; }
	  else                { hx_max->SetBinContent(nrow+1,fHistAddrMaxAdc->GetBinContent(hw))   ;	}
	}
    }
  
  cx_max->cd();
  Char_t hx_titl[50] ; sprintf(hx_titl,"max adc in pads at x=%5.1f mm",xval);
  hx_max->SetTitle(hx_titl);
  hx_max->SetXTitle("row");
  if(!GetPedestals()) hx_max->SetYTitle("max adc (baseline sub.)");
  else                hx_max->SetYTitle("max adc ");
  hx_max->SetMinimum(0.01);
  hx_max->Draw("l");
  
  if(setrange)
    {
      hx_max->GetXaxis()->SetRange(prof_x_xmin,prof_x_xmax);
      hx_max->SetMinimum(prof_x_ymin);
      hx_max->SetMaximum(prof_x_ymax);
    }
  
  cx_max->Update();
 
  crow_time->cd();
  Char_t title[256];
  Char_t title_max[256];
  if(ROC==1) {sprintf(title,"%s Row=%d",((TH2*)select)->GetTitle(),row)   ;    sprintf(title_max,"IROC  max/sum Row=%d",row   );} 
  else       {sprintf(title,"%s Row=%d",((TH2*)select)->GetTitle(),row-63);    sprintf(title_max,"OROC  max/sum Row=%d",row-63);}
  if(fVerb) cout << " set name " << endl;
  

  // row vs time
  crow_time->cd();
  hrow_time->SetTitleSize(0.04);
  hrow_time->SetTitle(title);
  hrow_time->SetYTitle("timbin");
  hrow_time->SetXTitle("pad in row");
  hrow_time->SetZTitle("signal (ADC)");

  hrow_time->GetXaxis()->SetTitleColor(1);
  hrow_time->SetMaximum(1000.0);
  hrow_time->SetMinimum(0.0);
  
  if(setrange) 
    {
      hrow_time->GetXaxis()->SetRange(row_time_xmin,row_time_xmax);
      hrow_time->GetYaxis()->SetRange(row_time_ymin,row_time_ymax);
      hrow_time->SetMinimum(row_time_zmin);
      hrow_time->SetMaximum(row_time_zmax);
    }

  hrow_time->Draw("COLZ");
  crow_time->UseCurrentStyle();
  crow_time->Update();

  // max and sum /////////////////////////
  crow_max->cd();
  if(setrange) {
    hrow_max->GetXaxis()->SetRange(prof_row_xmin,prof_row_xmax);
    hrow_max->SetMinimum(prof_row_ymin);
    hrow_max->SetMaximum(prof_row_ymax);
  }
  hrow_max->SetTitleSize(0.04);
  hrow_max->SetTitle(title);
  hrow_max->SetYTitle("max adc");
  hrow_max->SetXTitle("pad in row");
  hrow_max->GetXaxis()->SetTitleColor(1);

  hrow_max->SetLineColor(2);
  hrow_max->Draw("l");                  
  crow_max->Update();
  
  return;
}

//__________________________________________________________________
void AliTPCMonitor::Write10bitChannel()
{
  
  // Write 10 bit words form histogram for active(last pointed)  channel 
  
  if(fPadUsedHwAddr==-1){ AliWarning(" No active pad "); return ;}
  
  Int_t  pad     = (Int_t)fMapHand->GetPad(   fPadUsedHwAddr); 
  Int_t  row     = (Int_t)fMapHand->GetPadRow(fPadUsedHwAddr); 
  Int_t  channel = (Int_t)fPadMapHw[fPadUsedHwAddr];
  
  Char_t filename_root[256];
  Char_t filename_dat[256];
  Char_t projhist[256];

  if(fPadUsedRoc==1) { sprintf(projhist,"ProjectionOROC"); }
  if(fPadUsedRoc==0) { sprintf(projhist,"ProjectionIROC"); }
  
  sprintf(filename_dat, "Channel_Run%05i_EventID_%i_Pad_%i_Row_%i.dat"      ,fRunId,fEventNumber,pad,row);
  sprintf(filename_root,"Channel_Run%05i_EventID_%i_Pad_%i_Row_%i.root"     ,fRunId,fEventNumber,pad,row);
  
  TH1D* hpr = 0; 
  if(hpr=(TH1D*)gROOT->Get(projhist)) 
    {
      // root histo 
      TFile f(filename_root,"recreate");
      hpr->Write();
      f.Close();

      // raw singal 
      ofstream datout(filename_dat,ios::out);
      datout <<"Timebin \t ADC value " << endl; 
      for(Int_t i = 1; i <GetTimeBins(); i++)
	{
	  datout << i << " \t \t " << fPad[channel][i] << endl; 
	}
      datout.close();
    }
  else
    {
      AliWarning("No projection histo found ");
    }
}

//__________________________________________________________________
void AliTPCMonitor::ExecTransform() 
{

  // Make Fast Fourier Transformation for active pad
  // fft is only performed for a data sample of size 2^n
  // reduce window according to largest  power of 2 which is smaller than the viewing  range 

  Char_t namecanv[256]; 
  Char_t projhist[256];
  Char_t name_htr_imag[256];
  Char_t name_htr_real[256];
  
  if(fPadUsedRoc==1) {    sprintf(namecanv,"coroc_ch_trans") ;    sprintf(projhist,"ProjectionOROC");  }
  if(fPadUsedRoc==0) {    sprintf(namecanv,"ciroc_ch_trans") ;    sprintf(projhist,"ProjectionIROC");  }

  TH1D*  hproj = 0;
  
  if((TH1D*)gROOT->Get(projhist)==0){AliWarning("Proj histo does not exist \n Move mouse over 2d histo choose channel \n and drag mouse form histo again!");  return ;}
  else      hproj = (TH1D*)gROOT->Get(projhist) ;

  
  if(fPadUsedRoc==1) {  sprintf(name_htr_imag,"htrans_imag_freq_oroc");    sprintf(name_htr_real,"htrans_real_freq_oroc");  } 
  else               {  sprintf(name_htr_imag,"htrans_imag_freq_iroc");    sprintf(name_htr_real,"htrans_real_freq_iroc");  }
  
  if( gROOT->Get(name_htr_imag))  delete  gROOT->Get(name_htr_imag);
  if( gROOT->Get(name_htr_real))  delete  gROOT->Get(name_htr_real);
  
  TCanvas *ctrans = 0;
  if(!(ctrans = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(namecanv))) 
    {
      ctrans = CreateCanvas(namecanv);
      ctrans->Divide(1,2);
    }
  
  Int_t bin_first  =  hproj->GetXaxis()->GetFirst();
  Int_t bin_last   =  hproj->GetXaxis()->GetLast();
  Int_t bins       =  bin_last -bin_first +1;
    
  Int_t power = 0;
  for(Int_t pot = 0; pot<=10 ; pot++) 
    {
      Int_t comp =  (Int_t)TMath::Power(2,pot);
      if(bins>=comp)power = pot;
    }
  
  bins = (Int_t)TMath::Power(2,power);
  
  // sampling frequency ;
  Double_t  delta_t = 1.0/(Float_t)GetSamplingFrequency();
  
  // output histo
  TH1D* htrans_real_freq = new TH1D(name_htr_real,name_htr_real,10000,-1/(2*delta_t),1/(2*delta_t));
  TH1D* htrans_imag_freq = new TH1D(name_htr_imag,name_htr_imag,10000,-1/(2*delta_t),1/(2*delta_t));


  Char_t title_real[256];
  Char_t title_imag[256];
  if(fPadUsedRoc==1) {    sprintf(title_real,"OROC DFT real part");  sprintf(title_imag,"OROC DFT imag part");  } 
  else {                  sprintf(title_real,"IROC DFT real part");  sprintf(title_imag,"IROC DFT imag part");  }
  
  htrans_real_freq->SetTitle(title_real);  htrans_real_freq->SetXTitle("f/hz");  htrans_real_freq->SetYTitle("z_{real}(f)");
  htrans_imag_freq->SetTitle(title_imag);  htrans_imag_freq->SetXTitle("f/hz");  htrans_imag_freq->SetYTitle("z_{imag}(f)");
  
  
  // create complex packed data array  
  const Int_t datasiz = 2*bins;
  Double_t*  data = new Double_t[datasiz];
  for(Int_t i=0;i<2*bins;i++)  { data[i]   =  0.0;}
  for(Int_t i=0;i<bins;i++)    { data[2*i] = (Double_t)hproj->GetBinContent(bin_first+i); }
  
  // make fourier transformation
  AliTPCMonitorFFT* four = new AliTPCMonitorFFT();
  four->ComplexRadix2ForwardWrap(data,1,bins);
  
  // write output  and fill histos forward  
  Double_t freq =  0.0;
  for(Int_t i=0;i<2*bins;i++) 
    {
      if(i<bins) 
	{
	  if(i<(bins/2))  { freq = i/(bins*delta_t)            ; } 
	  else            { freq = -1*((bins-i)/(bins*delta_t)); }
	  htrans_real_freq->Fill( freq,data[2*i]  );
	  htrans_imag_freq->Fill( freq,data[2*i+1]);
	  
	  
	}
    }
  
  ctrans->cd(1);
  htrans_real_freq->Draw();
  ctrans->cd(2);
  htrans_imag_freq->Draw();
  ctrans->Update();
  delete four;
  delete data;
}

//__________________________________________________________________
void AliTPCMonitor::ShowSel(Int_t* comp_val)               
{
  
  // Show only selected components 
  // First restore original histogram from clone 
  // Than remove all not matching pads form histos 
  
  Int_t   connector   =  0;
  Int_t   FECnr       =  0;
  Int_t   Altro       =  0;
  Int_t   FEClocbran  =  0;
  Int_t   branch      =  0;
  Short_t rcuget      =  0;
  Int_t   emptyI      =  1;
  Int_t   index       = -1;
  Int_t   hwadd       =  0;

  Float_t maxiroc     =  fHistIROCClone->GetMaximum();
  Float_t maxoroc     =  fHistOROCClone->GetMaximum();
 
  
  //  restore original histos 
  for(Int_t row = 0; row<kNRowsIroc;  row++) 
    {
    for(Int_t pad = 0; pad<kNPadsIroc;  pad++) 
      {
      index = (Int_t)fHistIROCIndex->GetCellContent(row+1,pad+1);
      if(index==-1)continue;
      else  fHistIROC->SetCellContent(row+1,pad+1,fHistIROCClone->GetCellContent(row+1,pad+1));
      }
    }
  for(Int_t row = 0; row<kNRowsOroc;  row++) 
    {
      for(Int_t pad = 0; pad<kNPadsOroc;  pad++) 
	{
	  index = (Int_t)fHistOROCIndex->GetCellContent(row+1,pad+1);
	  if(index==-1)continue;
	  else    fHistOROC->SetCellContent(row+1,pad+1,fHistOROCClone->GetCellContent(row+1,pad+1));
	}
    }
  
  
  // remove not matching entries from fHistIROC/fHistOROC
  
  TH2F* fHist       =0;
  TH2S* fHistIndex  =0;
  Int_t npads       =0;
  Int_t subrows     =0;
  
  for(Int_t row = 0; row< (kNRowsIroc + kNRowsOroc);  row++) 
    {
      if(row<kNRowsIroc) {  fHist=fHistIROC ; fHistIndex = fHistIROCIndex; npads = kNPadsIroc; subrows =0         ;}
      else               {  fHist=fHistOROC ; fHistIndex = fHistOROCIndex; npads = kNPadsOroc; subrows =kNRowsIroc;}
      
      for(Int_t pad = 0; pad<npads;  pad++) 
	{
	  index    = (Int_t)fHistIndex->GetCellContent(row -subrows +1,pad+1);
	  if(index==-1)  continue ;  
	  hwadd    = (Int_t)fHistChannelTime->GetCellContent(index,0);
	  
	  // global fecnr
	  FECnr     =  fMapHand->GetFECfromHw(hwadd);
	  if(comp_val[0]!=-1 && FECnr!=comp_val[0])      {	  fHist->SetCellContent(row-subrows+1,pad+1,0);	          continue;	}
	  
	  // rcu
	  rcuget      = (hwadd & AliTPCMonitorAltro::GetHwMaskRCU())>> 12;
	  if(comp_val[1]!=-1 && rcuget!=comp_val[1])     {	  fHist->SetCellContent(row-subrows+1,pad+1,0);	          continue;	}
	  
	  // branch
	  branch    =  fMapHand->U2fGetBranch(FECnr) ;
	  if(comp_val[2]!=-1 && branch!=comp_val[2])     {	  fHist->SetCellContent(row-subrows+1,pad+1,0);	          continue;	}
	  
	  // local fec
	  FEClocbran=   fMapHand->U2fGetFECinBranch(FECnr) ;
	  if(comp_val[3]!=-1 && FEClocbran!=comp_val[3]) { 	  fHist->SetCellContent(row-subrows+1,pad+1,0); 	  continue; 	}
	  
	  // connector
	  connector =  fMapHand->GetFECconnector(hwadd);
	  if(comp_val[4]!=-1 && connector!=comp_val[4])  { 	  fHist->SetCellContent(row-subrows+1,pad+1,0); 	  continue; 	}
	  
	  // Altro chip
	  Altro     =  fMapHand->GetAltro(hwadd); 	
	  if(comp_val[5]!=-1 && Altro!=comp_val[5])      { 	  fHist->SetCellContent(row-subrows+1,pad+1,0); 	  continue; 	}
	  emptyI =0;
	}
    }
  
  TCanvas* c1 =0;
  TCanvas* c2 =0;
  if(gROOT->GetListOfCanvases()->FindObject("ciroc")) 
    {
      c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("ciroc");
      c1->cd() ;
      fHistIROC->Draw("COLZ");
      fHistIROC->SetMaximum(maxiroc);
      fHistIROC->SetMinimum(0.0);   
      c1->Update();
    }
  if(gROOT->GetListOfCanvases()->FindObject("coroc")) 
    {
      c2 =  (TCanvas*)gROOT->GetListOfCanvases()->FindObject("coroc");
      c2->cd() ;
      fHistOROC->Draw("COLZ");
      fHistOROC->SetMaximum(maxoroc);
      fHistOROC->SetMinimum(0.0);
      c2->Update();
    }
  return ;
}

//__________________________________________________________________
void AliTPCMonitor::ResizeCanv() 
{
  // Resize canvases and delete some of them
  
  Char_t carry1[100];  
  sprintf(carry1,".x %s/TPC/AliTPCMonitorExec.C(1)",gSystem->Getenv("ALICE_ROOT"));
  Char_t carry3[100];  
  sprintf(carry3,".x %s/TPC/AliTPCMonitorExec.C(2)",gSystem->Getenv("ALICE_ROOT"));
  if(fVerb) cout <<  " canv 1 " << endl;
  
  if(gROOT->GetListOfCanvases()->FindObject(        "coroc_ch")) {  delete gROOT->GetListOfCanvases()->FindObject("coroc_ch") ; }
  if(gROOT->GetListOfCanvases()->FindObject(        "ciroc_ch")) {  delete gROOT->GetListOfCanvases()->FindObject("ciroc_ch") ; }
  
  // for 2dim plots delete create and draw again
  if(gROOT->GetListOfCanvases()->FindObject("ciroc")) 
    {
      delete  gROOT->GetListOfCanvases()->FindObject("ciroc");
      TCanvas* ciroc = CreateCanvas("ciroc");
      ciroc->cd();
      fHistIROC->Draw("COLZ");
      ciroc->AddExec("pad",carry1);
      ciroc->AddExec("row",carry3);
      fExecPlaneMax=1;
      ciroc->Update();
    }
  // for 2dim plots delete create and draw again
  if(gROOT->GetListOfCanvases()->FindObject("coroc")) 
    {
      delete gROOT->GetListOfCanvases()->FindObject("coroc");
      TCanvas* coroc = CreateCanvas("coroc");
      coroc->cd();
      fHistOROC->Draw("COLZ");
      
      coroc->AddExec("pad",carry1);
      coroc->AddExec("row",carry3);
      coroc->Update();
      fExecPlaneMax=1;
    } 
  
  if(gROOT->GetListOfCanvases()->FindObject(       "cbasemean")) {    delete gROOT->GetListOfCanvases()->FindObject("cbasemean"); }
  if(gROOT->GetListOfCanvases()->FindObject(           "cbase")) {    delete gROOT->GetListOfCanvases()->FindObject("cbase");}
  if(gROOT->GetListOfCanvases()->FindObject(        "cbaserms")) {    delete gROOT->GetListOfCanvases()->FindObject("cbaserms");  }
  if(gROOT->GetListOfCanvases()->FindObject(            "cmax")) {    delete gROOT->GetListOfCanvases()->FindObject("cmax");      }
  if(gROOT->GetListOfCanvases()->FindObject(            "csum")) {    delete gROOT->GetListOfCanvases()->FindObject("csum");      }
  if(gROOT->GetListOfCanvases()->FindObject(  "ciroc_ch_trans")) {    delete gROOT->GetListOfCanvases()->FindObject("ciroc_ch_trans");}
  if(gROOT->GetListOfCanvases()->FindObject(  "coroc_ch_trans")) {    delete gROOT->GetListOfCanvases()->FindObject("coroc_ch_trans");}
  if(gROOT->GetListOfCanvases()->FindObject(       "crow_time")) {    delete gROOT->GetListOfCanvases()->FindObject("crow_time"); }
  if(gROOT->GetListOfCanvases()->FindObject(        "crow_max")) {    delete gROOT->GetListOfCanvases()->FindObject("crow_max");  }
  if(gROOT->GetListOfCanvases()->FindObject(          "cx_max")) {    delete gROOT->GetListOfCanvases()->FindObject("cx_max");    }
  if(gROOT->GetListOfCanvases()->FindObject(        "crmsoroc")) {    delete gROOT->GetListOfCanvases()->FindObject("crmsoroc");         fExecPadOrocRms = 0;   }
  if(gROOT->GetListOfCanvases()->FindObject(        "crmsiroc")) {    delete gROOT->GetListOfCanvases()->FindObject("crmsiroc");         fExecPadIrocRms = 0;   }
  if(gROOT->GetListOfCanvases()->FindObject(        "crow_max")) {    delete gROOT->GetListOfCanvases()->FindObject("crow_max");  }
  if(gROOT->GetListOfCanvases()->FindObject(        "crow_max")) {    delete gROOT->GetListOfCanvases()->FindObject("crow_max");  }

}




//__________________________________________________________________
Int_t AliTPCMonitor::ExecProcess() 
{
  // Executable for global Histogram
  // Will be called from /TPC/AliTPCMonitorExec.C(3)
  // Call ProcessEvent for same event and sector pointed at 
  
  Int_t side   = 0;
  Int_t sector = 0;
  
  Int_t event = gPad->GetEvent();
  if(event != 61)  return -1;
  
  TObject *select = gPad->GetSelected();
  if(!select)  return -1;
  if(!select->InheritsFrom("TH2")) {gPad->SetUniqueID(0);    return -1;  }
  if(       strcmp(select->GetName(),"hglobal" )==0 || ( strcmp(select->GetName(),"SIDE A" )==0) ) side = 0;
  else  if( strcmp(select->GetName(),"hglobal2")==0 || ( strcmp(select->GetName(),"SIDE C" )==0) ) side = 1;

  // get position
  Int_t   px    = gPad->GetEventX();
  Int_t   py    = gPad->GetEventY();
  Float_t upy   = gPad->AbsPixeltoY(py);
  Float_t upx   = gPad->AbsPixeltoX(px);
  Float_t y     = gPad->PadtoY(upy);
  Float_t x     = gPad->PadtoX(upx);

  Int_t testy = ((TH2*)select)->GetYaxis()->FindBin(y);
  Int_t testx = ((TH2*)select)->GetXaxis()->FindBin(x);
  if(((TH2*)select)->GetCellContent(testx,testy)==0) return -1 ;

  Float_t alpha = 0.0;
  if(x>0.0 && y > 0.0)    alpha = TMath::Abs(TMath::ATan(TMath::Abs(x/y)));
  if(x>0.0 && y < 0.0)    alpha = TMath::Abs(TMath::ATan(TMath::Abs(y/x)));
  if(x<0.0 && y < 0.0)    alpha = TMath::Abs(TMath::ATan(TMath::Abs(x/y)));
  if(x<0.0 && y > 0.0)    alpha = TMath::Abs(TMath::ATan(TMath::Abs(y/x)));
  
  if(x>0.0 && y < 0.0)    alpha += (    TMath::Pi()/2);
  if(x<0.0 && y < 0.0)    alpha += (    TMath::Pi());
  if(x<0.0 && y > 0.0)    alpha += (1.5*TMath::Pi());
  
  sector =   (Int_t)(alpha/(2*TMath::Pi()/18.0));
  if(alpha> (sector+0.5)*(2*TMath::Pi()/18.0))  sector+=1;
  
  if(sector==18 && side ==0 ) {
    AliWarning("There was a wromg assignment of sector 0 with sector 18. Check sectors");
    sector =0;
  }
  
  sector = (18-sector +4)%18;
  SetLastSector(sector+ side*18);
  SetProcNextEvent(0);
  
  if(fVerb) cout << "AliTPCMonitor::ExecProcess()  next side          " <<   side    << " next sector        " <<    sector  << endl;
  
  return (Int_t)ProcessEvent();

}

//__________________________________________________________________
Int_t AliTPCMonitor::GetRCUPatch(Int_t runid, Int_t eqid)
{
  
  // Return RCU patch index for given equipment id eqid 
  Int_t patch = 0;
  //if(runid>=704)
  if(runid>=0)
    {
      if(eqid>=1000) return 0;
      patch = fMapEqidsRcu[eqid] ;
    }
  else
    {
      if(eqid==408) {patch =  13*6+4 +0;  }
      if(eqid==409) {patch =  13*6+5 +0;  }
      if(eqid==509) {patch =  13*6+0 +0;  }
      if(eqid==512) {patch =  13*6+3 +0;  }
      if(eqid==513) {patch =  13*6+1 +0;  }
      if(eqid==517) {patch =  13*6+2 +0;  }
      
      if(eqid==404) {patch =  4*6+5 +0;   }
      if(eqid==504) {patch =  4*6+4 +0;   }
      if(eqid==407) {patch =  4*6+3 +0;   }  
      if(eqid==503) {patch =  4*6+2 +0;   }
      if(eqid==508) {patch =  4*6+1 +0;   }
      if(eqid==506) {patch =  4*6+0 +0;   }
    }
  return patch;
}

//__________________________________________________________________
void AliTPCMonitor::DumpHeader(AliRawReaderRoot * fReaderROOT)
{
  // Dump Event header for ROOT format
  
  cout << "EventHeader     : fReaderROOT->GetEquipmentSize()            :" << fReaderROOT->GetEquipmentSize()        << endl;
  cout << "EventHeader     : fReaderROOT->GetType()                     :" << fReaderROOT->GetType()                 << endl;
  cout << "EventHeader     : fReaderROOT->GetRunNumber()                :" << fReaderROOT->GetRunNumber()            << endl;
  cout << "EventHeader     : fReaderROOT->GetEventId()                  :" << *(fReaderROOT->GetEventId())           << endl;
  cout << "EventHeader     : fReaderROOT->GetLDCId()                    :" << fReaderROOT->GetLDCId()                << endl;
  cout << "EventHeader     : fReaderROOT->GetGDCId()                    :" << fReaderROOT->GetGDCId()                << endl;
}

//__________________________________________________________________
void AliTPCMonitor::DumpHeader(AliTPCMonitorDateFormat* DateForm)
{
  // Dump Event header for DATE format
  
  cout << "EquipmentHeader : DateForm->GetEquipmentSize()               :" << DateForm->GetEquipmentSize()          << endl;
  cout << "EquipmentHeader : DateForm->GetEquipmentType()               :" << DateForm->GetEquipmentType()          << endl;
  cout << "EquipmentHeader : DateForm->GetEquipmentID()                 :" << DateForm->GetEquipmentID()            << endl;
  cout << "EquipmentHeader : DateForm->GetEquipmentTypeAttribute()      :" << DateForm->GetEquipmentTypeAttribute() << endl;
  cout << "EquipmentHeader : DateForm->GetEquipmentBasicSize()          :" << DateForm->GetEquipmentBasicSize()     << endl;
  cout << "EquipmentHeader : DateForm->GetEquipmentHeaderSize()         :" << DateForm->GetEquipmentHeaderSize()    << endl;
  cout << "EquipmentHeader : DateForm->IsLastSubEventHeader()           :" << DateForm->IsLastSubEventHeader()      << endl;
}

//__________________________________________________________________
Double_t AliTPCMonitor::Gamma4(Double_t* x, Double_t* par) {
  
  // Gamma4 function used to fit signals
  // Defined in sections: diverging branch set to 0 
  
  Double_t val  = 0.0; 
  if(x[0] > par[1])
    val = par[0]*exp(4.0)* pow((x[0]-par[1])/par[2],4)*exp(-4.0*(x[0]-par[1])/par[2])+ par[3];
  else 
    val = 0;
  return val;
}

//__________________________________________________________________
TCanvas* AliTPCMonitor::CreateCanvas(Char_t* name)
{
  // Create Canvases 
  
  TCanvas* canv =0;
  
  Int_t xoffset  = GetCanvasXOffset();
  Int_t xsize    = GetCanvasXSize();
  Int_t ysize    = GetCanvasYSize();
  Int_t xspace   = GetCanvasXSpace();
  Int_t yspace   = GetCanvasYSpace();
  
  // ROC 2dim max distribution
  if(     strcmp(name,"coroc"         )==0) {    canv   = new TCanvas("coroc"         ,"coroc"      ,                   -1+xoffset,(Int_t)(yspace+0.5*ysize) ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(strcmp(name,"ciroc"         )==0) {    canv   = new TCanvas("ciroc"         ,"ciroc"      ,                   -1+xoffset,                     0    ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  // ROC  2dim rms distribution
  else if(strcmp(name,"crmsoroc"      )==0) {    canv   = new TCanvas("crmsoroc"      ,"crmsoroc"   ,                   -1+xoffset,(Int_t)(yspace+0.5*ysize) ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  } 
  else if(strcmp(name,"crmsiroc"      )==0) {    canv   = new TCanvas("crmsiroc"      ,"crmsiroc"   ,                   -1+xoffset,                        0 ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  // Global ADC max Histos
  else if(strcmp(name,"SIDE C all"    )==0) {    canv   = new TCanvas("SIDE C all"    ,"SIDE C all" ,   (Int_t)(3*xspace+ xoffset),(Int_t)(yspace+0.5*ysize) ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(strcmp(name,"SIDE A all"    )==0) {    canv   = new TCanvas("SIDE A all"    ,"SIDE A all" ,   (Int_t)(3*xspace+ xoffset),                        0 ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  // 1 dim max sum basekine distribution
  else if(strcmp(name,"cmax"          )==0) {    canv   = new TCanvas("cmax"          ,"cmax"       ,                   -1+xoffset,                 3*yspace ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else if(strcmp(name,"csum"          )==0) {    canv   = new TCanvas("csum"          ,"csum"       ,               xspace+xoffset,                 3*yspace ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else if(strcmp(name,"cbasemean"     )==0) {    canv   = new TCanvas("cbasemean"     ,"cbasemean"  ,             2*xspace+xoffset,                 3*yspace ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }  
  else if(strcmp(name,"cbaserms"      )==0) {    canv   = new TCanvas("cbaserms"      ,"cbaserms"   ,             3*xspace+xoffset,                 3*yspace ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  // Projections of single channel
  else if(strcmp(name,"coroc_ch"      )==0) {    canv   = new TCanvas("coroc_ch"      ,"coroc_ch"   ,   (Int_t)(1.5*xspace+xoffset),(Int_t)(yspace+0.5*ysize),(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(strcmp(name,"ciroc_ch"      )==0) {    canv   = new TCanvas("ciroc_ch"      ,"ciroc_ch"   ,   (Int_t)(1.5*xspace+xoffset),                       0 ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  // FFT for single channel
  else if(strcmp(name,"coroc_ch_trans")==0) {    canv   = new TCanvas("coroc_ch_trans","coroc_ch_trans",(Int_t)(3.0*xspace+xoffset),(Int_t)(yspace+0.5*ysize),(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  else if(strcmp(name,"ciroc_ch_trans")==0) {    canv   = new TCanvas("ciroc_ch_trans","ciroc_ch_trans",(Int_t)(3.0*xspace+xoffset),                       0 ,(Int_t)(1.5*xsize),(Int_t)(1.5*ysize)); return canv;  }
  // row profile histograms
  else if(strcmp(name,"crow_time"     )==0) {    canv   = new TCanvas("crow_time"     ,"crow_time"  ,              1*xspace+xoffset,         2*yspace +ysize ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else if(strcmp(name,"crow_max"      )==0) {    canv   = new TCanvas("crow_max"      ,"crow_max"   ,              2*xspace+xoffset,         2*yspace +ysize ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else if(strcmp(name,"cx_max"        )==0) {    canv   = new TCanvas("cx_max"        ,"cx_max"     ,              3*xspace+xoffset,         2*yspace +ysize ,(Int_t)(1.0*xsize),(Int_t)(1.0*ysize)); return canv;  }
  else                                      {    cout   << " Warning Canvas name unknown "  << endl;                                                                                                  return 0   ;  }
}

//__________________________________________________________________
void AliTPCMonitor::WriteHistos()  
{
  // Writes all available histograms to a file in current directory
  // File name will be specified by : sector, side, runid and eventnumber 
  
  if(GetEventProcessed())
    {
      AliInfo("Write histos to file");
      Char_t name[256]; 
      sprintf(name,"SIDE_%i_SECTOR_%02i_RUN_%05i_EventID_%06i.root",(GetLastSector()/18),(GetLastSector()%18),fRunId,fEventNumber);
      TFile* f = new TFile(name,"recreate");
      for(Int_t i =0; i<fHistList->GetEntries(); i++)
	{
	  if(((TH1*)fHistList->At(i))!=0)    
	    {
	      ((TH1*)fHistList->At(i))->Write();
	    }
	}
      f->ls();
      f->Close();
    }
  else
    { 
      AliError("No Event Processed : Chose Format , File and push 'Next Event' ");
    }
}


//__________________________________________________________________
TH1* AliTPCMonitor::GetHisto(char* histname) 
{
  
  // Returns histogram specified by histname
  // check available names for histos in CreateHistos()
  
  TH1* hist = 0;
  if((TH1*)fHistList->FindObject(histname))
    {
      hist = (TH1*)fHistList->FindObject(histname);
    }
  else
    {
      cout << " AliTPCMonitor::GetHisto :: Can not find histo with name " << histname << endl;
    }
  return hist ;
}
