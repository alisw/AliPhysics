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

/**********************************************
Purpose:  A modified macro to generate MappingTable
          in the following form

detElemId+iPlane+iX+iY  Ix  IY  X  Y  B/NB pcbtype

          for the trigger stations.

Created:  13/05/2007
Modified: 31/08/2007 (Adopted to AliRoot v4-06-Release)
Modified: 08/09/2007 (To produce mapping on the basis of uniqueID as detElemId+iPlane+localCard+ibitXY )
Author:   Indranil Das, HEP, SINP, Kolkata
Email:    indra.das@saha.ac.in | indra.ehep@gmail.com
***********************************************/

#include <iostream> 

//STEER 
#include "AliCDBManager.h"

//MUON
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVCalibParam.h"

//MUON/mapping 
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"


using namespace std;


Bool_t CreateTriggerRecoLookupTables(TString transformFileName = "geometry.root")
{

  Char_t filename1[20], filename2[20];
  Int_t chamberId;

  Int_t runNumber = 0;
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliMUONCalibrationData cd(runNumber);

  AliMpSegmentation *mpSegFactory = AliMpSegmentation::ReadData(); 

  AliMUONGeometryTransformer* chamberGeometryTransformer = new AliMUONGeometryTransformer();
  chamberGeometryTransformer->LoadGeometryData(transformFileName);
  
    sprintf(filename1,"Lut%d.dat",2*(10+0));
    //ofstream fout1(filename1);
    FILE *fout1 = fopen(filename1,"w");

    sprintf(filename2,"Lut%d.dat",2*(10+0)+1);
    //ofstream fout2(filename2);
    FILE *fout2 = fopen(filename2,"w");

  for(Int_t iCh = 0; iCh < 4; iCh++){ // max 4

    chamberId = iCh + 10;
    
    cout<<"Running for Chamber : "<<chamberId<<endl;

    AliMpDEIterator it;
    for ( it.First(chamberId); ! it.IsDone(); it.Next() ) {
    
      Int_t detElemId = it.CurrentDEId();
      
      Double_t globPos[3] = {0.0, 0.0, 0.0};
      Double_t locPos[3] = {0.0, 0.0, 0.0};

      chamberGeometryTransformer->GetDetElement(detElemId)->Local2Global(locPos[0],locPos[1],locPos[2],globPos[0],globPos[1],globPos[2]);
       
      for(Int_t iPlane = 0 ; iPlane <= 1 ; iPlane++){
	
	AliMp::CathodType cath;

	if(iPlane == 0)
	  cath = AliMp::kCath0 ;
	else
	  cath = AliMp::kCath1 ;

	AliMpVSegmentation* seg = mpSegFactory->CreateMpSegmentation(detElemId, cath);

	Int_t maxIX = seg->MaxPadIndexX();  
	Int_t maxIY = seg->MaxPadIndexY(); 
	Int_t idManuChannel,idetElemId;

	Double_t realX, realY, realZ;
	Double_t localX, localY, localZ;
	Double_t padSizeX, padSizeY;
	Int_t pcbType;
	Int_t locCard, ibitxy;

	//Pad Info of a segment
	for(Int_t iX = 0; iX<= maxIX ; iX++){
	  for(Int_t iY = 0; iY<= maxIY ; iY++){
	    if(seg->HasPad(AliMpIntPair(iX,iY))){
	      AliMpPad pad = seg->PadByIndices(AliMpIntPair(iX,iY),kFALSE);
	      
	      locCard = pad.GetLocation(0).GetFirst();
	      ibitxy = pad.GetLocation(0).GetSecond();

 	      idetElemId = detElemId%1000;
	      idetElemId &= 0x1FF ;
	      iPlane &= 0x1 ;
	      locCard &= 0xFF ;
	      ibitxy &= 0xF ;
	      
	      idManuChannel &= 0x0;
	      idManuChannel = (idManuChannel|idetElemId)<<1;  
 	      idManuChannel = (idManuChannel|iPlane)<<8;  
	      idManuChannel = (idManuChannel|locCard)<<4 ;
	      idManuChannel |= ibitxy;
	      
	      localX = pad.Position().X();
	      localY = pad.Position().Y();
	      localZ = 0.0;

	      chamberGeometryTransformer->Local2Global(detElemId,localX,localY,localZ,
						       realX,realY,realZ);
	      padSizeX = 2.0*pad.Dimensions().X();
	      padSizeY = 2.0*pad.Dimensions().Y();

	      if(iPlane == 0 ){
		if(padSizeX==17.0)
		  pcbType = 0;
		else if(padSizeX==34.0)
		  pcbType = 1;
		else
		  pcbType = 2;
	      }
	      else{
		if(padSizeY==51.0)
		  pcbType = 0;
		else
		  pcbType = 1;
	      }
	      
	      idetElemId %= 100;

	      if(idetElemId<5 || idetElemId > 13){
 		fprintf(fout1,"%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n",idManuChannel,iX,iY,realX,realY,realZ,pcbType,iPlane);
	      }
	      else{
		fprintf(fout2,"%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n",idManuChannel,iX,iY,realX,realY,realZ,pcbType,iPlane);
	      }// HasPad Condn
	      
	    }
	  }// iY loop
	}// iX loop
	
      }// iPlane
    } // detElemId loop

  }// ichamber loop

  fclose(fout1);
  fclose(fout2);
  
  return kTRUE;
}
