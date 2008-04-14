/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*********************************************
Purpose:  A macro to generate LookupTable
          in the following form
buspatchId+manuid+channelId  buspatchId Ix  IY  X  Y  B/NB

Created:  7/10/2005
Modified: 22/12/2005
Modified: 09/02/2006
Modified: 09/04/2007
Modified: 24/08/2007 (To adopt to AliRoot v4-06-Release)

Run Info: To run this code copy "rootlogon.C" 
          in the current directory from $ALICE_ROOT/MUON 
           then compile and run using
          .L CreateHitRecoLookupTables.C+

Author:   Indranil Das, HEP, SINP, Kolkata
Email:    indra.das@saha.ac.in
***********************************************/

#include <iostream> 

//STEER 
#include "AliCDBManager.h"
#include "AliGeomManager.h"

//MUON
#include "AliMUONGeometryTransformer.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVCalibParam.h"

//MUON/mapping 
#include "AliMpCDB.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"
#include <AliMpDetElement.h>

using namespace std;

Bool_t CreateHitRecoLookupTables(TString CDBPath = "local://$ALICE_ROOT", Int_t run = 0, Bool_t warn = kTRUE)
{
  Char_t filename1[50];
  Int_t chamberId;
  
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage(CDBPath.Data());
  cdbManager->SetRun(run);

  if (! AliMpCDB::LoadDDLStore(warn)){
    cerr<<__FILE__<<": Failed to Load DDLStore specified for CDBPath "<<CDBPath<<", and Run : "<<run<<endl;
    return kFALSE;
  }

  AliMpSegmentation *mpSegFactory = AliMpSegmentation::Instance(); 
  AliGeomManager::LoadGeometry();
  AliMUONGeometryTransformer* chamberGeometryTransformer = new AliMUONGeometryTransformer();
  if(! chamberGeometryTransformer->LoadGeometryData()){
    cerr<<__FILE__<<": Failed to Load Geomerty Data "<<endl;
    return kFALSE;
  }
  
  int maxDDL = 8;
  FILE *fout[maxDDL];
  for(int iDDL = 0;iDDL<maxDDL; iDDL++){
    sprintf(filename1,"Lut%d.dat",iDDL+13);
    fout[iDDL] = fopen(filename1,"w");
  }
  
  AliMUONCalibrationData calibData(run);

  int totMaxIX = -1;
  int totMaxIY = -1;

  for(Int_t iCh = 6; iCh < 10; iCh++){ // max 4

    chamberId = iCh ;

    AliMpDEIterator it;
    for ( it.First(chamberId); ! it.IsDone(); it.Next() ) {
    
      Int_t detElemId = it.CurrentDEId();
      int iDDL = AliMpDDLStore::Instance()->GetDetElement(detElemId)->GetDdlId() - 12 ;
      for(Int_t iCath = 0 ; iCath <= 1 ; iCath++){
	
	AliMp::CathodType cath;

	if(iCath == 0)
	  cath = AliMp::kCath0 ;
	else
	  cath = AliMp::kCath1 ;

	const AliMpVSegmentation* seg = mpSegFactory->GetMpSegmentation(detElemId, cath);
	AliMp::PlaneType plane = seg->PlaneType(); 
	Int_t maxIX = seg->MaxPadIndexX();  
	Int_t maxIY = seg->MaxPadIndexY(); 
	if(maxIX > totMaxIX)
	  totMaxIX = maxIX;
	if(maxIY > totMaxIY)
	  totMaxIY = maxIY;

	Int_t idManuChannel, manuId, channelId, buspatchId;
	float padSizeX, padSizeY;
	float halfPadSize ;
	Double_t realX, realY, realZ;
	Double_t localX, localY, localZ;
	Float_t calibA0Coeff,calibA1Coeff,pedestal,sigma;
	Int_t thresold,saturation;

// 	cout<<"Running for detElemId :"<<detElemId<<", and plane : "<<plane<<endl;
	//Pad Info of a segment to print in lookuptable
	for(Int_t iX = 0; iX<= maxIX ; iX++){
	  for(Int_t iY = 0; iY<= maxIY ; iY++){
	    if(seg->HasPad(AliMpIntPair(iX,iY))){
	      AliMpPad pad = seg->PadByIndices(AliMpIntPair(iX,iY),kFALSE);

	      // Getting Manu id
	      manuId = pad.GetLocation().GetFirst();
	      manuId &= 0x7FF; // 11 bits 

	      buspatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
	      
	      // Getting channel id
	      channelId =  pad.GetLocation().GetSecond();
	      channelId &= 0x3F; // 6 bits
	      
	      idManuChannel &= 0x0;
	      idManuChannel = (idManuChannel|buspatchId)<<11;  
	      idManuChannel = (idManuChannel|manuId)<<6 ;
	      idManuChannel |= channelId ;
	      
	      localX = pad.Position().X();
	      localY = pad.Position().Y();
	      localZ = 0.0;

 	      chamberGeometryTransformer->Local2Global(detElemId,localX,localY,localZ,
 						       realX,realY,realZ);

	      padSizeX = pad.Dimensions().X();
	      padSizeY = pad.Dimensions().Y();

	      calibA0Coeff = (calibData.Gains(detElemId,manuId))->ValueAsFloat(channelId,0) ;
	      calibA1Coeff = (calibData.Gains(detElemId,manuId))->ValueAsFloat(channelId,1) ;
	      thresold = (calibData.Gains(detElemId,manuId))->ValueAsInt(channelId,2) ;
	      saturation = (calibData.Gains(detElemId,manuId))->ValueAsInt(channelId,4) ;

	      pedestal = (calibData.Pedestals(detElemId,manuId))->ValueAsFloat(channelId,0);
	      sigma = (calibData.Pedestals(detElemId,manuId))->ValueAsFloat(channelId,1);

	      if(plane==0)
		halfPadSize = padSizeX;
	      else
		halfPadSize = padSizeY;

	      fprintf(fout[iDDL],"%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%d\t%d\n",
		      idManuChannel,detElemId,iX,iY,realX,realY,realZ,
		      halfPadSize,plane,pedestal,sigma,calibA0Coeff,calibA1Coeff,thresold,saturation);

	    }// HasPad Condn
	  }// iY loop
	}// iX loop
	
      }// iPlane

    } // detElemId loop

//     fclose(fout1);

  }// ichamber loop

  for(int iDDL = 0;iDDL<maxDDL; iDDL++){
    fclose(fout[iDDL]);
  }

//   cout<<"TotMaxIX : "<<totMaxIX<<", and totMaxIY : "<<totMaxIY<<endl;
  
  return kTRUE;
}
