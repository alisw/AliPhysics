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
// L0-Trigger simulator of PHOS
//
//-- Author: Satoshi Yano (Hiroshima University)

#include <iostream>
using std::cout;
using std::endl;
#include "AliCaloTriggerSimulator.h"
#include "TRefArray.h"
#include "TH1F.h"
#include "TF1.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCluster.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEmcCalibData.h"
#include "AliPHOSCalibData.h"
#include "AliCDBStorage.h"
#include "AliAODCaloTrigger.h"

ClassImp(AliCaloTriggerSimulator) 
//===============================================
AliCaloTriggerSimulator::AliCaloTriggerSimulator():
  fAOD(NULL),
  fCluster(NULL),
  fCells(NULL),
  fPHOSGeo(NULL),
  fCalibDataEmc(NULL),
  fPHOSCalibData(NULL),
  fCDBstorage(NULL),
  fThreshold(384)
{
  
  //fAOD  = (AliAODEvent*)fEvent->Clone();

} 

//===============================================
AliCaloTriggerSimulator::AliCaloTriggerSimulator(AliAODEvent* fEvent):
  fAOD(fEvent),
  fCluster(NULL),
  fCells(NULL),
  fPHOSGeo(NULL),
  fCalibDataEmc(NULL),
  fPHOSCalibData(NULL),
  fCDBstorage(NULL),
  fThreshold(384)
{

  fCells = dynamic_cast<AliAODCaloCells*> (fEvent->GetPHOSCells());
  //this->Init();

} 
//===============================================
void AliCaloTriggerSimulator::SetDB(AliPHOSGeometry* fGeo,AliPHOSEmcCalibData* fCalibData1, AliPHOSCalibData* fCalibData2, AliCDBStorage *fCDB)
{
  
  fPHOSGeo             = fGeo;
  fCalibDataEmc        = fCalibData1;
  fPHOSCalibData       = fCalibData2;
  fCDBstorage          = fCDB;

}
//===============================================
void AliCaloTriggerSimulator::SetEvent(AliAODEvent* fEvent)
{
  fAOD = fEvent;
  fCells = dynamic_cast<AliAODCaloCells*> (fEvent->GetPHOSCells());
}
//===============================================
AliAODCaloTrigger* AliCaloTriggerSimulator::CreateTriggerMap()
{

  if(!fPHOSGeo || !fCalibDataEmc || !fPHOSCalibData || !fCDBstorage || !fAOD){
    cout<<"Several DB files are not included. Therefore, the code stop."<<endl;
    return 0;
  }
  
  AliAODCaloCells* cells = dynamic_cast<AliAODCaloCells*> (fAOD->GetPHOSCells());
  
    
  TRefArray* fClusters = new TRefArray();
  fAOD->GetPHOSClusters(fClusters);

  Bool_t RoI[3][8]={};

  for(Int_t index = 0; index < fClusters->GetEntriesFast(); ++index){
    AliAODCluster* cluster = (AliAODCluster*) fClusters->At(index);  
    Float_t  position[3];
    cluster->GetPosition(position);
    TVector3 global(position) ;
    Int_t relId[4] = {0,0,0,0};
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    Int_t    mod = relId[0] ;
    Int_t    col = relId[2] ;
    Int_t    row = relId[3] ;
    Int_t    tru = GetTRUNum(col,row);
    RoI[mod-1][tru-1] = true;
  }

  Int_t size = 0;
  Int_t mod[1000];
  Int_t abs[1000];
  Float_t amp[1000];

  for(Int_t mod1=1; mod1<=3; mod1++){

    if(mod1==2) continue;

    for(Int_t cellX1=1; cellX1<=64; cellX1++){
      for(Int_t cellZ1=1; cellZ1<=56; cellZ1++){

        Int_t trelId[4]={0,0,0,0};
        Int_t tAbsId=0;

        if(cellX1%2==0 || cellZ1%2==0) continue;
        if(cellX1==15)                 continue;
        if(cellX1==31)                 continue;
        if(cellX1==47)                 continue;
        if(cellX1==63)                 continue;
        if(cellZ1==27)                 continue;
        if(cellZ1==55)                 continue;
	
        //Sliding window left bottom corner
        trelId[0]=mod1;
        trelId[1]=0;
        trelId[2]=cellX1;
        trelId[3]=cellZ1;

        Int_t tru1 = GetTRUNum(trelId[2],trelId[3]);

        if(!RoI[mod1-1][tru1-1])
	  continue;
        
	fPHOSGeo->RelToAbsNumbering(trelId,tAbsId);

        Int_t SlidingRelId[4];
        Double_t SumADC4x4=0;
        Double_t SumClust4x4=0;
        Double_t Sum4x4Analogue=0;

        TH1F* h4x4AnalogSum;
        Int_t nTileSum = 0;
	
        Bool_t fisLGFlag = false;

	for(Int_t ix=0; ix<4; ++ix){
          for(Int_t iz=0; iz<4; ++iz){
	
            SlidingRelId[0] = mod1;
            SlidingRelId[1] = 0;
            SlidingRelId[2] = trelId[2]+ix;
            SlidingRelId[3] = trelId[3]+iz;
            
	    Int_t    AbsId=0;
            fPHOSGeo->RelToAbsNumbering(SlidingRelId,AbsId);
            
	    Double_t cellAmp     = cells->GetCellAmplitude(AbsId);
	    if(!(cellAmp>0)) continue;

	    Bool_t   isHG      = 0;
            Double_t CDBcalibEmc = fCalibDataEmc->GetADCchannelEmc(SlidingRelId[0],SlidingRelId[3],SlidingRelId[2]);

	    if(cellAmp/CDBcalibEmc > 1024.){
	      isHG = 0;
	      fisLGFlag = true;
	    }
	    else{
	      isHG = 1;
	    }
	    
	    Double_t adcCh       = 0;
	    
	    //High gain channel
	    adcCh = cellAmp / (CDBcalibEmc);
            
	    if(cellAmp < 1E-4) continue;
	    
            SumADC4x4   += adcCh;
            SumClust4x4 += cellAmp;
	    
            ++nTileSum;
	    
          }//end of loop iz
        }//end of loop ix
	
        if(SumADC4x4<1E-4)
	  continue;
	
	Double_t ratioTRUandFEEADC = 11.2/5.;
        SumADC4x4 /= ratioTRUandFEEADC;

        if(fThreshold<SumADC4x4){
	  amp[size] = SumClust4x4;
	  mod[size] = trelId[0];
	  abs[size] = tAbsId;
	  ++size;
	}
      
      }//end of loop cellZ1
    }//end of loop cellX1
  }//end of loop mod1
  
  AliAODCaloTrigger* trigger = new AliAODCaloTrigger();
  trigger->Allocate(size);
  trigger->Reset();
  
  for(Int_t i=0; i<size; ++i){
    Float_t time = 1.E-9;
    Int_t ntrgtimes = 0;
    Int_t trgtimes[ntrgtimes];
    Int_t trgts     = 0;
    Int_t trgbits   = 0;
    trigger->Add(mod[i],abs[i],amp[i],time,trgtimes,ntrgtimes,trgts,trgbits);
  }
  
  trigger->Reset();

  if(trigger)
    return trigger;
  else 
    return 0;

}

Int_t AliCaloTriggerSimulator::GetTRUNum(Int_t cellX, Int_t cellZ)
{
  //Return TRU region number for given cell.
  //cellX: [1-64], cellZ: [1-56]
  Int_t iTRU=-111;
  //RCU0: TRU 1,2
  if(1<=cellX&&cellX<=16) {
    if(1<=cellZ&&cellZ<=28) iTRU=2;
    else iTRU=1;
  }
  //RCU1: TRU 3,4
  if(17<=cellX&&cellX<=32) {
    if(1<=cellZ&&cellZ<=28) iTRU=4;
    else iTRU=3;
  }
  //RCU2: TRU 5,6
  if(33<=cellX&&cellX<=48) {
    if(1<=cellZ&&cellZ<=28) iTRU=6;
    else iTRU=5;
  }
  //RCU3: TRU 7,8
  if(49<=cellX&&cellX<=64) {
    if(1<=cellZ&&cellZ<=28) iTRU=8;
    else iTRU=7;
  }

  return iTRU;
}
