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
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
*/

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2I.h> 
#include <TH2D.h> 
#include <TGraph.h> 
#include <TParameter.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliVZEROQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliVZERORawStream.h"
#include "AliVZEROReconstructor.h"
#include "event.h"


ClassImp(AliVZEROQADataMakerRec)
           
//____________________________________________________________________________ 
  AliVZEROQADataMakerRec::AliVZEROQADataMakerRec() : 
	AliQADataMakerRec(AliQA::GetDetName(AliQA::kVZERO), "VZERO Quality Assurance Data Maker"),
	fCalibData(0x0),
    fEvent(0)
    
{
		AliInfo("Construct VZERO QA Object");
  // constructor
  
   for(Int_t i=0; i<64; i++){  
       fEven[i] = 0;   
       fOdd[i]  = 0;  }
  
   for(Int_t i=0; i<128; i++){  
       fADC_Mean[i] = 0.0;   }	
}

//____________________________________________________________________________ 
  AliVZEROQADataMakerRec::AliVZEROQADataMakerRec(const AliVZEROQADataMakerRec& qadm) :
  AliQADataMakerRec(),
	fCalibData(0x0),
    fEvent(0)
  
{
  //copy constructor 
   SetName((const char*)qadm.GetName()) ; 
   SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliVZEROQADataMakerRec& AliVZEROQADataMakerRec::operator = (const AliVZEROQADataMakerRec& qadm )
{
  // Equal operator
  this->~AliVZEROQADataMakerRec();
  new(this) AliVZEROQADataMakerRec(qadm);
  return *this;
}

//____________________________________________________________________________
AliVZEROCalibData* AliVZEROQADataMakerRec::GetCalibData() const

{
  AliCDBManager *man = AliCDBManager::Instance();
  
  //man->SetDefaultStorage("local://$ALICE_ROOT");

  AliCDBEntry *entry=0;
  
  entry = man->Get("VZERO/Calib/Data",fRun);

  // Retrieval of data in directory VZERO/Calib/Data:

  AliVZEROCalibData *calibdata = 0;

  if (entry) calibdata = (AliVZEROCalibData*) entry->GetObject();
  if (!calibdata)  AliFatal("No calibration data from calibration database !");

  return calibdata;
}
 
//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kVZERO, task, list) ;
  
  if(task == AliQA::kRAWS){
  	  int NMaxBin = GetRawsData(kPedestal_Time_Int0)->GetNbinsY();
	  if(fCurrentCycle%NMaxBin==0) {
  		GetRawsData(kPedestal_Time_Int0)->Reset();
  		GetRawsData(kPedestal_Time_Int1)->Reset();
  		GetRawsData(kChargeEoI_Time_Int0)->Reset();
  		GetRawsData(kChargeEoI_Time_Int1)->Reset();
	  }
	  TH1D* hProj;
	  char Name[50];
	  for(Int_t iChannel=0; iChannel<64; iChannel++) {
  		for(Int_t integrator=0;integrator<2;integrator++){
			sprintf(Name,"Ped_%d_%d",iChannel,integrator);
  			hProj = ((TH2I*)GetRawsData((integrator == 0 ? kPedestal_Cycle_Int0 : kPedestal_Cycle_Int1)))->ProjectionY(Name,iChannel+1,iChannel+1);
			((TH2D*)GetRawsData((integrator == 0 ? kPedestal_Time_Int0 : kPedestal_Time_Int1)))->Fill((double)iChannel,(double)(fCurrentCycle%NMaxBin),(double)hProj->GetMean());
			delete hProj;

 			sprintf(Name,"Charge_%d_%d",iChannel,integrator);
 			hProj = ((TH2I*)GetRawsData((integrator == 0 ? kChargeEoI_Cycle_Int0 : kChargeEoI_Cycle_Int1)))->ProjectionY(Name,iChannel+1,iChannel+1);
			((TH2D*)GetRawsData((integrator == 0 ? kChargeEoI_Time_Int0 : kChargeEoI_Time_Int1)))->Fill((double)iChannel,(double)(fCurrentCycle%NMaxBin),hProj->GetMean());
			delete hProj;
  		}
	  }
  } else if (task == AliQA::kESDS) {
  }
  
}

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::InitESDs()
{
  //Create histograms to control ESD
	Bool_t expert   = kTRUE ; 
	
	TH2D * h2d;
	TH1I * h1i;
	TH1D * h1d;
		
  h1i = new TH1I("H1I_Cell_Multiplicity_V0A", "Cell Multiplicity in V0A", 35, 0, 35) ;  
  Add2ESDsList(h1i, kCellMultiV0A, !expert)  ;  
  h1i->GetXaxis()->SetTitle("Multiplicity (Nb of Cell)");
                                                                                                        
  h1i = new TH1I("H1I_Cell_Multiplicity_V0C", "Cell Multiplicity in V0C", 35, 0, 35) ;  
  Add2ESDsList(h1i, kCellMultiV0C, !expert)  ;  
  h1i->GetXaxis()->SetTitle("Multiplicity (Nb of Cell)");
   
  h1d = new TH1D("H1D_MIP_Multiplicity_V0A", "MIP Multiplicity in V0A", 1000, 0, 1000) ;  
  Add2ESDsList(h1d, kMIPMultiV0A, !expert)  ;  
  h1d->GetXaxis()->SetTitle("Multiplicity (Nb of MIP)");
  
  h1d = new TH1D("H1D_MIP_Multiplicity_V0C", "MIP Multiplicity in V0C", 1000, 0, 1000) ;  
  Add2ESDsList(h1d, kMIPMultiV0C, !expert)  ;  
  h1d->GetXaxis()->SetTitle("Multiplicity (Nb of MIP)");

  h2d = new TH2D("H2D_MIP_Multiplicity_Channel", "MIP Multiplicity per Channel",64, 0, 64, 100, 0, 100) ;  
  h2d->GetXaxis()->SetTitle("Channel");
  h2d->GetYaxis()->SetTitle("Multiplicity (Nb of MIP)");
  Add2ESDsList(h2d, kMIPMultiChannel, !expert)  ;  
  
  h1d = new TH1D("H1D_BBFlag_Counters", "BB Flag Counters",64, 0, 64) ;  
  h1d->GetXaxis()->SetTitle("Channel");
  Add2ESDsList(h1d, kBBFlag, !expert)  ;  
  
  h1d = new TH1D("H1D_BGFlag_Counters", "BG Flag Counters",64, 0, 64) ;  
  h1d->GetXaxis()->SetTitle("Channel");
  Add2ESDsList(h1d, kBGFlag, !expert)  ;  
  
  h2d = new TH2D("H2D_Charge_Channel", "ADC Charge per channel",64, 0, 64, 1024, 0, 1024) ;  
  h2d->GetXaxis()->SetTitle("Channel");
  h2d->GetYaxis()->SetTitle("Charge (ADC counts)");
  Add2ESDsList(h2d, kChargeChannel, !expert)  ;  
  
  h2d = new TH2D("H2D_Time_Channel", "Time per channel",64, 0, 64, 820, 0, 410) ;  
  h2d->GetXaxis()->SetTitle("Channel");
  h2d->GetYaxis()->SetTitle("Time (ns)");
  Add2ESDsList(h2d, kTimeChannel, !expert)  ;  
  
  h1d = new TH1D("H1D_V0A_Time", "Mean V0A Time",2048, 0., 409.6);
  h1d->GetXaxis()->SetTitle("Time (ns)");
  Add2ESDsList(h1d,kESDV0ATime, !expert); 
  
  h1d = new TH1D("H1D_V0C_Time", "Mean V0C Time",2048, 0., 409.6);
  h1d->GetXaxis()->SetTitle("Time (ns)");
  Add2ESDsList(h1d,kESDV0CTime, !expert); 
  
  h1d = new TH1D("H1D_Diff_Time", "Diff Time V0A - V0C",2*2048, -409.6, 409.6);
  h1d->GetXaxis()->SetTitle("Diff Time V0A - V0C (ns)");
  Add2ESDsList(h1d,kESDDiffTime, !expert); 
	
}

//____________________________________________________________________________ 
 void AliVZEROQADataMakerRec::InitRaws()
 {
   // create Raws histograms in Raws subdir
	 fCalibData = GetCalibData();

  Bool_t expert   = kTRUE ; 
  Bool_t saveCorr = kTRUE ; 

 char Name[50] , Title[100];
 const Int_t kNintegrator  =    2;
 
 const Int_t kNTdcTimeBins  = 2048;
 const Int_t kTdcTimeMin    =    0;
 const Int_t kTdcTimeMax    = 4096;
 const Int_t kNTdcWidthBins =  128;
 const Int_t kTdcWidthMin   =    0;
 const Int_t kTdcWidthMax   =  128;
 const Int_t kNChargeBins   = 1024;
 const Int_t kChargeMin     =    0;
 const Int_t kChargeMax     = 1024;
 const Int_t kNChannelBins  =   64;
 const Int_t kChannelMin    =    0;
 const Int_t kChannelMax 	=   64;
 const Int_t kNPedestalBins =  200;
 const Int_t kPedestalMin   =    0;
 const Int_t kPedestalMax   =  200;
 const Int_t kTimeMin =   0;
 const Int_t kTimeMax = 100;
 const Int_t kNMIPBins = 200;
 const Int_t kMIPMin   =   0;
 const Int_t kMIPMax   = 200;


 TH2I * h2i;
 TH2D * h2d;
 TH1I * h1i;
 TH1D * h1d;

 int iHisto =0;
 
 
  // Creation Cell Multiplicity Histograms
  h1i = new TH1I("H1I_Multiplicity_V0A", "Cell Multiplicity in V0A", 35, 0, 35) ;  
  Add2RawsList(h1i,kMultiV0A, !expert, saveCorr);   iHisto++;
  h1i = new TH1I("H1I_Multiplicity_V0C", "Cell Multiplicity in V0C", 35, 0, 35) ;  
  Add2RawsList(h1i,kMultiV0C, !expert, saveCorr);   iHisto++;
 
  // Creation Total Charge Histograms
  h1d = new TH1D("H1D_Charge_V0A", "Total Charge in V0A", 2048, 0, 32768) ;  
  Add2RawsList(h1d,kChargeV0A, !expert, saveCorr);   iHisto++;
  h1d = new TH1D("H1D_Charge_V0C", "Total Charge in V0C", 2048, 0, 32768) ;  
  Add2RawsList(h1d,kChargeV0C, !expert, saveCorr);   iHisto++;
  h1d = new TH1D("H1D_Charge_V0", "Total Charge in V0", 2048, 0, 65536) ;  
  Add2RawsList(h1d,kChargeV0, !expert, saveCorr);   iHisto++;
  
  // Creation of MIP Histograms
  h1d = new TH1D("H1D_MIP_V0A", "Total MIP in V0A", 2*kNMIPBins,kMIPMin ,32*kMIPMax) ;  
  Add2RawsList(h1d,kRawMIPV0A, !expert, saveCorr);   iHisto++;
  h1d = new TH1D("H1D_MIP_V0C", "Total MIP in V0C", 2*kNMIPBins,kMIPMin ,32*kMIPMax) ;  
  Add2RawsList(h1d,kRawMIPV0C, !expert, saveCorr);   iHisto++;
  h1d = new TH1D("H1D_MIP_V0", "Total MIP in V0", 2*kNMIPBins,kMIPMin ,32*kMIPMax) ;  
  Add2RawsList(h1d,kRawMIPV0, !expert, saveCorr);   iHisto++;
  h2d = new TH2D("H2D_MIP_Channel", "Nb of MIP per channel", kNChannelBins, kChannelMin, kChannelMax,kNMIPBins,kMIPMin ,kMIPMax) ;  
  Add2RawsList(h2d,kRawMIPChannel, expert, !saveCorr);   iHisto++;
  
 
 for(Int_t iInt=0;iInt<kNintegrator;iInt++){
    // Creation of Pedestal histograms 
    sprintf(Name,"H2I_Pedestal_Int%d",iInt);
    sprintf(Title,"Pedestal (Int%d)",iInt);
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax,kNPedestalBins,kPedestalMin ,kPedestalMax );
	Add2RawsList(h2i,(iInt == 0 ? kPedestal_Int0 : kPedestal_Int1), expert, !saveCorr); iHisto++;
	
	// Creation of temporary Pedestal histo used for the mean versus time histogram. This histogram will be reset at the end of each cycle
    sprintf(Name,"H2I_Pedestal_Cycle_Int%d",iInt);
    sprintf(Title,"One Cycle Pedestal (Int%d)",iInt);
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax,kNPedestalBins,kPedestalMin ,kPedestalMax );
	Add2RawsList(h2i,(iInt == 0 ? kPedestal_Cycle_Int0 : kPedestal_Cycle_Int1), expert, !saveCorr); iHisto++;
		
	// Creation of Pedestal versus time graph.
    sprintf(Name,"H2D_Pedestal_Time_Int%d",iInt);
    sprintf(Title,"Pedestal Versus Time (Int%d)",iInt);
 	h2d = new TH2D(Name, Title,kNChannelBins, kChannelMin, kChannelMax,kTimeMax,kTimeMin ,kTimeMax );
	Add2RawsList(h2d,(iInt == 0 ? kPedestal_Time_Int0 : kPedestal_Time_Int1), expert, !saveCorr); iHisto++;

   // Creation of Charge EoI histograms 
    sprintf(Name,"H2I_ChargeEoI_Int%d",iInt);
    sprintf(Title,"Charge EoI (Int%d)",iInt);
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax, kNChargeBins, kChargeMin, kChargeMax);
	Add2RawsList(h2i,(iInt == 0 ? kChargeEoI_Int0 : kChargeEoI_Int1), !expert, !saveCorr); iHisto++;

   // Creation of temporary Charge EoI histograms used for the mean versus time histogram. This histogram will be reset at the end of each cycle
    sprintf(Name,"H2I_ChargeEoI_Cycle_Int%d",iInt);
    sprintf(Title,"One Cycle Charge EoI (Int%d)",iInt);
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax, kNChargeBins, kChargeMin, kChargeMax);
	Add2RawsList(h2i,(iInt == 0 ? kChargeEoI_Cycle_Int0 : kChargeEoI_Cycle_Int1), expert, !saveCorr); iHisto++;
		
	// Creation of Chercge EoI versus time graph.
    sprintf(Name,"H2D_ChargeEoI_Time_Int%d",iInt);
    sprintf(Title,"Charge EoI Versus Time (Int%d)",iInt);
 	h2d = new TH2D(Name, Title,kNChannelBins, kChannelMin, kChannelMax,kTimeMax,kTimeMin ,kTimeMax );
	Add2RawsList(h2d,(iInt == 0 ? kChargeEoI_Time_Int0 : kChargeEoI_Time_Int1), expert, !saveCorr); iHisto++;
    
    sprintf(Name,"H2I_ChargeEoI_BB_Int%d",iInt);
    sprintf(Title,"Charge EoI w/ BB Flag (Int%d)",iInt);
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax, kNChargeBins, kChargeMin, kChargeMax);
	Add2RawsList(h2i,(iInt == 0 ? kChargeEoI_BB_Int0 : kChargeEoI_BB_Int1), expert, !saveCorr); iHisto++;
    
    sprintf(Name,"H2I_ChargeEoI_BG_Int%d",iInt);
    sprintf(Title,"Charge EoI w/ BG Flag (Int%d)",iInt);
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax, kNChargeBins, kChargeMin, kChargeMax);
	Add2RawsList(h2i,(iInt == 0 ?  kChargeEoI_BG_Int0: kChargeEoI_BG_Int1), expert, !saveCorr); iHisto++;

    // Creation of Charge versus LHC Clock histograms 
    sprintf(Name,"H2D_ChargeVsClock_Int%d",iInt);
    sprintf(Title,"Charge Versus LHC-Clock (Int%d)",iInt);
 	h2d = new TH2D(Name, Title,kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
	Add2RawsList(h2d,(iInt == 0 ? kChargeVsClock_Int0 : kChargeVsClock_Int1 ), expert, !saveCorr); iHisto++;
	
	// Creation of Minimum Bias Charge histograms 
	for(Int_t iBB=0;iBB<2;iBB++){
		for(Int_t iBG=0;iBG<2;iBG++){
			sprintf(Name,"H2I_ChargeMB_BB%d_BG%d_Int%d",iBB,iBG,iInt);
			sprintf(Title,"MB Charge (BB=%d, BG=%d, Int=%d)",iBB,iBG,iInt);
			h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax,kNChargeBins, kChargeMin, kChargeMax);
			int idx;
			if(iInt==0){
				if(iBB==0){
					if(iBG==0) idx = kChargeMB_BB0_BG0_Int0;
					else idx = kChargeMB_BB0_BG1_Int0;
				} else {
					if(iBG==0) idx = kChargeMB_BB1_BG0_Int0;
					else idx = kChargeMB_BB1_BG1_Int0;
				}
			} else {
				if(iBB==0){
					if(iBG==0) idx = kChargeMB_BB0_BG0_Int1;
					else idx = kChargeMB_BB0_BG1_Int1;
				} else {
					if(iBG==0) idx = kChargeMB_BB1_BG0_Int1;
					else idx = kChargeMB_BB1_BG1_Int1;
				}
			}
			Add2RawsList(h2i,idx, expert, !saveCorr); iHisto++;
		}
	}
	
 }
 
// Creation of Time histograms 
	sprintf(Name,"H2I_Width");
	sprintf(Title,"HPTDC Width");
	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax);
 	Add2RawsList(h2i,kWidth, expert, !saveCorr); iHisto++;

 	sprintf(Name,"H2I_Width_BB");
 	sprintf(Title,"HPTDC Width w/ BB Flag condition");
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax);
 	Add2RawsList(h2i,kWidth_BB, expert, !saveCorr); iHisto++;

 	sprintf(Name,"H2I_Width_BG");
 	sprintf(Title,"HPTDC Width w/ BG Flag condition");
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax);
 	Add2RawsList(h2i,kWidth_BG, expert, !saveCorr); iHisto++;

 	sprintf(Name,"H2I_HPTDCTime");
 	sprintf(Title,"HPTDC Time");
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
 	Add2RawsList(h2i,kHPTDCTime, !expert, !saveCorr); iHisto++;

 	sprintf(Name,"H2I_HPTDCTime_BB");
 	sprintf(Title,"HPTDC Time w/ BB Flag condition");
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
 	Add2RawsList(h2i,kHPTDCTime_BB, expert, !saveCorr); iHisto++;

 	sprintf(Name,"H2I_HPTDCTime_BG");
 	sprintf(Title,"HPTDC Time w/ BG Flag condition");
 	h2i = new TH2I(Name, Title,kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
 	Add2RawsList(h2i,kHPTDCTime_BG, expert, !saveCorr); iHisto++;
	
 	sprintf(Name,"H1D_V0A_Time");
 	sprintf(Title,"V0A Time");
 	h1d = new TH1D(Name, Title,kNTdcTimeBins, kTdcTimeMin/10, kTdcTimeMax/10);
 	Add2RawsList(h1d,kV0ATime, !expert, saveCorr); iHisto++;
	
 	sprintf(Name,"H1D_V0C_Time");
 	sprintf(Title,"V0C Time");
 	h1d = new TH1D(Name, Title,kNTdcTimeBins, kTdcTimeMin/10, kTdcTimeMax/10);
 	Add2RawsList(h1d,kV0CTime, !expert, saveCorr); iHisto++;
	
 	sprintf(Name,"H1D_Diff_Time");
 	sprintf(Title,"Diff V0A-V0C Time");
 	h1d = new TH1D(Name, Title,2*kNTdcTimeBins, -kTdcTimeMax/10, kTdcTimeMax/10);
 	Add2RawsList(h1d,kDiffTime, !expert, saveCorr); iHisto++;
	
 
 	// Creation of Flag versus LHC Clock histograms 
 	sprintf(Name,"H2D_BBFlagVsClock");
 	sprintf(Title,"BB-Flags Versus LHC-Clock");
 	h2d = new TH2D(Name, Title,kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
 	Add2RawsList(h2d,kBBFlagVsClock, expert, !saveCorr); iHisto++;
	
 	sprintf(Name,"H2D_BGFlagVsClock");
 	sprintf(Title,"BG-Flags Versus LHC-Clock");
 	h2d = new TH2D(Name, Title,kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
 	Add2RawsList(h2d,kBGFlagVsClock, expert, !saveCorr); iHisto++;
	 
 	AliInfo(Form("%d Histograms has been added to the Raws List",iHisto));
 }

//____________________________________________________________________________
void AliVZEROQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs
  UInt_t EventType = esd->GetEventType();

  switch (EventType){
  	case PHYSICS_EVENT:
  		AliESDVZERO *esdVZERO=esd->GetVZEROData();
   
		if (!esdVZERO) break;
		  
    	GetESDsData(kCellMultiV0A)->Fill(esdVZERO->GetNbPMV0A());
    	GetESDsData(kCellMultiV0C)->Fill(esdVZERO->GetNbPMV0C());  
    	GetESDsData(kMIPMultiV0A)->Fill(esdVZERO->GetMTotV0A());
    	GetESDsData(kMIPMultiV0C)->Fill(esdVZERO->GetMTotV0C());  
    	
		Float_t TimeV0A = 0., TimeV0C =0., DiffTime;
		Int_t iTimeV0A=0, iTimeV0C=0;
		
		for(Int_t i=0;i<64;i++) {
			GetESDsData(kMIPMultiChannel)->Fill((Float_t) i,(Float_t) esdVZERO->GetMultiplicity(i));
		   	GetESDsData(kChargeChannel)->Fill((Float_t) i,(Float_t) esdVZERO->GetAdc(i));
			if(esdVZERO->GetBBFlag(i)) GetESDsData(kBBFlag)->Fill((Float_t) i);
		  	if(esdVZERO->GetBGFlag(i)) GetESDsData(kBGFlag)->Fill((Float_t) i);

			Float_t Time = (Float_t) esdVZERO->GetTime(i)/10.; //Convert in ns:  1 TDC channel = 100ps 
			GetESDsData(kTimeChannel)->Fill((Float_t) i,Time);

			if(Time>0.){
				if (i<32) {
					iTimeV0C++;
					TimeV0C += Time;
				}else{
					iTimeV0A++;
					TimeV0A += Time;
				}
			}
    	}
		if(iTimeV0A>0) TimeV0A /= iTimeV0A; 
		else TimeV0A = -1.;
		if(iTimeV0C>0) TimeV0C /= iTimeV0C;
		else TimeV0C = -1.;
		if(TimeV0A<0. || TimeV0C<0.) DiffTime = -10000.;
		else DiffTime = TimeV0A - TimeV0C;
				
		GetESDsData(kESDV0ATime)->Fill(TimeV0A);
		GetESDsData(kESDV0CTime)->Fill(TimeV0C);
		GetESDsData(kESDDiffTime)->Fill(DiffTime);
		
		break;
	}  
  
}

//____________________________________________________________________________
 void AliVZEROQADataMakerRec::MakeRaws(AliRawReader* rawReader)
 {
  //Fill histograms with Raws, computes average ADC values dynamically (pedestal subtracted)
                  
  rawReader->Reset() ; 
  AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader); 
  rawStream->Next();
  
  eventTypeType EventType = rawReader->GetType();


  	Int_t mulV0A = 0 ; 
  	Int_t mulV0C = 0 ; 
	Double_t TimeV0A =0., TimeV0C = 0.;
	UInt_t iTimeV0A=0, iTimeV0C=0;
	Double_t ChargeV0A=0., ChargeV0C=0.;
	Double_t MIPV0A=0., MIPV0C=0.;

	Double_t DiffTime=-100000.;

  
  switch (EventType){
  	case PHYSICS_EVENT:
       int iFlag=0;
	   int pedestal;
	   int integrator;
	   Bool_t BBFlag;	 
       Bool_t BGFlag;	 
       UInt_t Time, Width;
	   Int_t MBCharge, Charge;
	   Int_t OfflineCh;
	   TH1D * hProj;

	   //printf("----------------------------------------\n");
       for(Int_t iChannel=0; iChannel<64; iChannel++) { // BEGIN : Loop over channels
		   
		   OfflineCh = rawStream->GetOfflineChannel(iChannel);
		   
		   // Fill Pedestal histograms
		   //-------------------------

	       for(Int_t j=15; j<21; j++) {
		       if((rawStream->GetBGFlag(iChannel,j) || rawStream->GetBBFlag(iChannel,j))) iFlag++;
	       }

	       if(iFlag == 0){ //No Flag found
		       for(Int_t j=15; j<21; j++){
	   		       pedestal=rawStream->GetPedestal(iChannel, j);
	   		       integrator = rawStream->GetIntegratorFlag(iChannel, j);

	   		       GetRawsData((integrator == 0 ? kPedestal_Int0 : kPedestal_Int1))->Fill(OfflineCh,pedestal);
	   		       GetRawsData((integrator == 0 ? kPedestal_Cycle_Int0 : kPedestal_Cycle_Int1))->Fill(OfflineCh,pedestal);
		       }
	       }

		   // Fill Charge EoI histograms
		   //---------------------------
		   // Look for the maximum in the LHC clock train
           Charge = 0;
           Int_t iClock  = 0;
           Int_t iCharge = 0;
           for(Int_t iEvent=0; iEvent<21; iEvent++){
               iCharge = rawStream->GetPedestal(iChannel,iEvent);
               if(iCharge>Charge)  {
	 	           Charge = iCharge;
	    	       iClock = iEvent;
	           }
           }   // End of maximum searching procedure

           integrator   = rawStream->GetIntegratorFlag(iChannel,iClock);
           BBFlag	 = rawStream->GetBBFlag(iChannel, iClock);
           BGFlag	 = rawStream->GetBGFlag(iChannel,iClock );

           GetRawsData((integrator == 0 ? kChargeEoI_Int0 : kChargeEoI_Int1))->Fill(OfflineCh,Charge);
		   if(BBFlag) GetRawsData((integrator == 0 ? kChargeEoI_BB_Int0 : kChargeEoI_BB_Int1))->Fill(OfflineCh,Charge);
           if(BGFlag) GetRawsData((integrator == 0 ? kChargeEoI_BG_Int0 : kChargeEoI_BG_Int1))->Fill(OfflineCh,Charge);

		   hProj = ((TH2I*)GetRawsData((integrator == 0 ? kPedestal_Int0 : kPedestal_Int1)))->ProjectionY("",OfflineCh+1,OfflineCh+1);
		   Double_t ped   = hProj->GetMean();
		   Double_t sigma = hProj->GetRMS();
		   delete hProj;

		   Double_t ChargeEoI = Charge - ped;
		   
		   // Calculation of the number of MIP
		   Double_t MIPEoI = ChargeEoI * fCalibData->GetMIPperADC(OfflineCh);

		   
	           if(Charge<1023 && ChargeEoI > 5.*sigma){ 
			   ((TH2I*)GetRawsData((integrator == 0 ? kChargeEoI_Cycle_Int0 : kChargeEoI_Cycle_Int1)))->Fill(OfflineCh,ChargeEoI);
			   ((TH2D*)GetRawsData(kRawMIPChannel))->Fill(OfflineCh,MIPEoI);
        	           if(OfflineCh<32) {
				   mulV0C++;
				   ChargeV0C += ChargeEoI;
				   MIPV0C += MIPEoI;
			   } else {
				   mulV0A++;
				   ChargeV0A += ChargeEoI;
				   MIPV0A += MIPEoI;
			   }
		   }

		   // Fill Charge Minimum Bias Histograms
		   //------------------------------------
		   int idx;
		   for(Int_t iBunch=0; iBunch<10; iBunch++){
			   integrator = rawStream->GetIntMBFlag(iChannel, iBunch);
			   BBFlag	   = rawStream->GetBBMBFlag(iChannel, iBunch);
			   BGFlag	   = rawStream->GetBGMBFlag(iChannel, iBunch);
			   MBCharge   = rawStream->GetChargeMB(iChannel, iBunch);

			   if(integrator==0){
				   if(BBFlag==0){
					   if(BGFlag==0) idx = kChargeMB_BB0_BG0_Int0;
					   else idx = kChargeMB_BB0_BG1_Int0;
				   } else {
					   if(BGFlag==0) idx = kChargeMB_BB1_BG0_Int0;
					   else idx = kChargeMB_BB1_BG1_Int0;
				   }
			   } else {
				   if(BBFlag==0){
					   if(BGFlag==0) idx = kChargeMB_BB0_BG0_Int1;
					   else idx = kChargeMB_BB0_BG1_Int1;
				   } else {
					   if(BGFlag==0) idx = kChargeMB_BB1_BG0_Int1;
					   else idx = kChargeMB_BB1_BG1_Int1;
				   }
			   }
			   GetRawsData(idx)->Fill(OfflineCh,MBCharge);
                   }   

		   // Fill HPTDC Time Histograms
		   //---------------------------

		   BBFlag	 = rawStream->GetBBFlag(iChannel, 10);
           BGFlag	 = rawStream->GetBGFlag(iChannel, 10);
           Time  = rawStream->GetTime(iChannel);
           Width = rawStream->GetWidth(iChannel);

		   if(Time>0.){
			   if (OfflineCh<32) {
				   iTimeV0C++;
				   TimeV0C += Time;
			   }else{
				   iTimeV0A++;
				   TimeV0A += Time;
			   }
		   }
           GetRawsData(kHPTDCTime)->Fill(OfflineCh,Time);
           GetRawsData(kWidth)->Fill(OfflineCh,Width);
           if(BBFlag) {
			   GetRawsData(kHPTDCTime_BB)->Fill(OfflineCh,Time);
			   GetRawsData(kWidth_BB)->Fill(OfflineCh,Width);
           }
		   if(BGFlag) {
			   GetRawsData(kHPTDCTime_BG)->Fill(OfflineCh,Time);
			   GetRawsData(kWidth_BG)->Fill(OfflineCh,Width);
           }

		   // Fill Flag and Charge Versus LHC-Clock histograms
		   //-------------------------------------------------

	   for(Int_t iEvent=0; iEvent<21; iEvent++){
               Charge = rawStream->GetPedestal(iChannel,iEvent);
               integrator   = rawStream->GetIntegratorFlag(iChannel,iEvent);
               BBFlag	 = rawStream->GetBBFlag(iChannel, iEvent);
               BGFlag	 = rawStream->GetBGFlag(iChannel,iEvent );

               ((TH2*) GetRawsData((integrator == 0 ? kChargeVsClock_Int0 : kChargeVsClock_Int1 )))->Fill(OfflineCh,(float)iEvent-10,(float)Charge);
               ((TH2*) GetRawsData(kBBFlagVsClock))->Fill(OfflineCh,(float)iEvent-10,(float)BBFlag);
               ((TH2*) GetRawsData(kBGFlagVsClock))->Fill(OfflineCh,(float)iEvent-10,(float)BGFlag);
           }

       }// END : Loop over channels

		if(iTimeV0A>0) TimeV0A /= (iTimeV0A * 10); // iTimeV0A Channels and divide by 10 to have the result in ns because 1 TDC Channel = 100 ps
		else TimeV0A = -1.;
	    if(iTimeV0C>0) TimeV0C /= (iTimeV0C * 10);
		else TimeV0C = -1.;
		if(TimeV0A<0. || TimeV0C<0.) DiffTime = -10000.;
		else DiffTime = TimeV0A - TimeV0C;
		
		GetRawsData(kV0ATime)->Fill(TimeV0A);
		GetRawsData(kV0CTime)->Fill(TimeV0C);
		GetRawsData(kDiffTime)->Fill(DiffTime);
		
		GetRawsData(kMultiV0A)->Fill(mulV0A);
		GetRawsData(kMultiV0C)->Fill(mulV0C);

		GetRawsData(kChargeV0A)->Fill(ChargeV0A);
		GetRawsData(kChargeV0C)->Fill(ChargeV0C);
		GetRawsData(kChargeV0)->Fill(ChargeV0A + ChargeV0C);
		
		GetRawsData(kRawMIPV0A)->Fill(MIPV0A);
		GetRawsData(kRawMIPV0C)->Fill(MIPV0C);
		GetRawsData(kRawMIPV0)->Fill(MIPV0A + MIPV0C);
		break;
	} //---------> END : SWITCH EVENT TYPE
	
	fEvent++; 
	TParameter<double> * p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kMultiV0A)->GetName()))) ; 
	p->SetVal((double)mulV0A) ; 

	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kMultiV0C)->GetName()))) ; 
	p->SetVal((double)mulV0C) ;                     

	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kChargeV0A)->GetName()))) ; 
	p->SetVal((double)ChargeV0A) ; 

	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kChargeV0C)->GetName()))) ; 
	p->SetVal((double)ChargeV0C) ;                     

	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kChargeV0)->GetName()))) ; 
	p->SetVal((double)(ChargeV0A + ChargeV0C)) ;                     
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kRawMIPV0A)->GetName()))) ; 
	p->SetVal((double)MIPV0A) ; 
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kRawMIPV0C)->GetName()))) ; 
	p->SetVal((double)MIPV0C) ;                     
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kRawMIPV0)->GetName()))) ; 
	p->SetVal((double)(MIPV0A + MIPV0C)) ;                     
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kV0ATime)->GetName()))) ; 
	p->SetVal((double)TimeV0A) ; 
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kV0CTime)->GetName()))) ; 
	p->SetVal((double)TimeV0C) ;                     
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(AliQA::kRAWS).Data(), GetRawsData(kDiffTime)->GetName()))) ; 
	p->SetVal((double)DiffTime) ;                     
	
  	delete rawStream; rawStream = 0x0;      


 }

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::StartOfDetectorCycle()
{
  // Detector specific actions at start of cycle
  
  // Resetting the histogram used to have the trend versus time.
  TH1* h;
  h = GetRawsData(kPedestal_Cycle_Int0);
  if(h) h->Reset();
  h = GetRawsData(kPedestal_Cycle_Int1); 
  if(h) h->Reset();
  h = GetRawsData(kChargeEoI_Cycle_Int0);
  if(h) h->Reset();
  h = GetRawsData(kChargeEoI_Cycle_Int1);
  if(h) h->Reset();

}
