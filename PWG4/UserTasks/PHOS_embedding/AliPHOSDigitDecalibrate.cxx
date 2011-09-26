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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  PHOS tender, recalibrate PHOS clusters                                   //
//  and do track matching                                                    //
//  Author : Dmitri Peressounko (RRC KI)                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliTender.h>
#include <TH2.h>
#include <TROOT.h>

#include "AliPHOSDigitDecalibrate.h"
#include "AliPHOSGeometry.h"

AliPHOSDigitDecalibrate::AliPHOSDigitDecalibrate() :
  AliTenderSupply()
  ,fPHOSGeo(0x0)
  ,fPHOSCalibData(0x0)
{
	//
	// default ctor
	//
}

//_____________________________________________________
AliPHOSDigitDecalibrate::AliPHOSDigitDecalibrate(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender)
  ,fPHOSGeo(0x0)
  ,fPHOSCalibData(0x0)
{
	//
	// named ctor
	//
}

//_____________________________________________________
AliPHOSDigitDecalibrate::~AliPHOSDigitDecalibrate()
{
  //Destructor
  for(int m=0; m<5; m++){
    if(hDec[m]){
      delete hDec[m];
      hDec[m]=0;
    }
  }
}

//_____________________________________________________
void AliPHOSDigitDecalibrate::Init()
{
  //
  // Initialise PHOS tender
  //
    

  
}

//_____________________________________________________
void AliPHOSDigitDecalibrate::ProcessEvent()
{
  //Choose PHOS clusters and recalibrate them
  //that it recalculate energy, position and distance 
  //to closest track extrapolation	

  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;

  // Init goemetry
  if(!fPHOSGeo){
    fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
  }


  AliESDCaloCells * cells = event->GetPHOSCells() ;

  for (Short_t icell = 0; icell < cells->GetNumberOfCells(); icell++) {
    Short_t id=0;
    Double_t time=0., amp=0. ;
    if (cells->GetCell(icell, id, amp, time) != kTRUE)
      break;

    Int_t relId[4] ;
    fPHOSGeo->AbsToRelNumbering(id,relId);
    Int_t   module = relId[0];
    Int_t   column = relId[3];
    Int_t   row    = relId[2];
   
//    amp=amp*fPHOSCalibData->GetADCchannelEmc(module,column,row);
    amp=amp*hDec[module-1]->GetBinContent(row,column);

    cells->SetCell(icell, id, amp, time);     

  }




}
//_____________________________________________________
void  AliPHOSDigitDecalibrate::SetDecalibration(Int_t mod, TH2F * dec){

  if(mod<0 || mod>4){
    printf("Don't know module %d \n",mod) ;
    return ;
  }
  if(dec==0)
    return ;
  if(hDec[mod])
    delete hDec[mod] ;
  gROOT->cd() ;
  hDec[mod] = new TH2F(*dec) ;
  char key[55] ;
  snprintf(key,55,"DecalibrationModule%d",mod) ;
  hDec[mod]->SetName(key) ;



}

