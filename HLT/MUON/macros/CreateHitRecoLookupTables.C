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
detElemId+manuid+channelId  buspatchId Ix  IY  X  Y  B/NB

Created:  7/10/2005
Modified: 22/12/2005
Modified: 09/02/2006
Modified: 09/04/2007
Modified: 24/08/2007 (To adopt to AliRoot v4-06-Release)

Run Info: To run this code copy "rootlogon.C" 
          in the current directory from $ALICE_ROOT/MUON 
          and specify the 
          transformFileName as "geometry.root" of the 
          simulation directory.Then compile and run using
          .L generateLookupTable.C++

Author:   Indranil Das, HEP, SINP, Kolkata
Email:    indra.das@saha.ac.in
***********************************************/

#include <iostream> 

//STEER 
#include "AliCDBManager.h"

//MUON
#include "AliMUONGeometryTransformer.h"

//MUON/mapping 
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"

class AliMpDDLStore ;

using namespace std;


Bool_t CreateHitRecoLookupTables(TString transformFileName = "geometry.root")
{
  Char_t filename1[20], filename2[20];
  Int_t chamberId;
  
  AliMpSegmentation *mpSegFactory = AliMpSegmentation::ReadData(); 

  AliMpDDLStore* fDDLStore = AliMpDDLStore::ReadData();
  
  AliMUONGeometryTransformer* chamberGeometryTransformer = new AliMUONGeometryTransformer();
  chamberGeometryTransformer->LoadGeometryData(transformFileName);
  
  for(Int_t iCh = 0; iCh < 4; iCh++){ // max 4
    
    sprintf(filename1,"Lut%d.dat",2*(6+iCh));
    FILE *fout1 = fopen(filename1,"w");

    sprintf(filename2,"Lut%d.dat",2*(6+iCh)+1);
    FILE *fout2 = fopen(filename2,"w");

    chamberId = iCh + 6;

    AliMpDEIterator it;
    for ( it.First(chamberId); ! it.IsDone(); it.Next() ) {
    
      Int_t detElemId = it.CurrentDEId();
      
      cout<<"Running for detElemId :"<<detElemId<<endl;
      
      for(Int_t iPlane = 0 ; iPlane <= 1 ; iPlane++){
	
	AliMp::CathodType cath;

	if(iPlane == 0)
	  cath = AliMp::kCath0 ;
	else
	  cath = AliMp::kCath1 ;
	
	AliMpVSegmentation* seg = mpSegFactory->CreateMpSegmentation(detElemId, cath);
	
	Int_t maxIX = seg->MaxPadIndexX();  
	Int_t maxIY = seg->MaxPadIndexY(); 
	Int_t idManuChannel, manuId, channelId,idetElemId;
	Int_t busPatchId;
	Double_t realX, realY, realZ;
	Double_t localX, localY, localZ;
	Double_t padSizeX, padSizeY;
	Int_t pcbType;

	//Pad Info of a segment
	for(Int_t iX = 0; iX<= maxIX ; iX++){
	  for(Int_t iY = 0; iY<= maxIY ; iY++){
	    if(seg->HasPad(AliMpIntPair(iX,iY))){
	      AliMpPad pad = seg->PadByIndices(AliMpIntPair(iX,iY),kFALSE);
	      
	      // Getting Manu id
	      manuId = pad.GetLocation().GetFirst();
	      manuId &= 0x7FF; // 11 bits 
	      
	      busPatchId = fDDLStore->GetBusPatchId(detElemId,manuId);
	      
	      // Getting channel id
	      channelId =  pad.GetLocation().GetSecond();
	      channelId &= 0x3F; // 6 bits
	      
	      idetElemId = detElemId%100;
	      
	      idManuChannel &= 0x0;
	      idManuChannel = (idManuChannel|idetElemId)<<11;  
	      idManuChannel = (idManuChannel|manuId)<<6 ;
	      idManuChannel |= channelId ;
	      
	      localX = pad.Position().X();
	      localY = pad.Position().Y();
	      localZ = 0.0;

 	      chamberGeometryTransformer->Local2Global(detElemId,localX,localY,localZ,
 						       realX,realY,realZ);

	      padSizeX = 2.0*pad.Dimensions().X();
	      padSizeY = 2.0*pad.Dimensions().Y();

	      if(iPlane == 0 ){
		if(padSizeX==2.5)
		  pcbType = 0;
		else if(padSizeX==5.0)
		  pcbType = 1;
		else
		  pcbType = 2;
	      }
	      else{
		if(padSizeY==2.5)
		  pcbType = 0;
		else if(padSizeY==5.0)
		  pcbType = 1;
		else
		  pcbType = 2;
	      }
		
	      if(idetElemId<7 || idetElemId > 19){
  		fprintf(fout2,"%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n",idManuChannel,iX,iY,realX,realY,realZ,pcbType,iPlane);
	      }
	      else{
 		fprintf(fout1,"%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n",idManuChannel,iX,iY,realX,realY,realZ,pcbType,iPlane);
	      }// HasPad Condn
	      
	    }
	    else{
	      cout<<"pad not found for iX :"<<iX<<", iY:"<<iY<<endl;
	    }
	  }// iY loop
	}// iX loop
	
      }// iPlane

    } // detElemId loop

    fclose(fout1);
    fclose(fout2);
  }// ichamber loop
  
  return kTRUE;
}
