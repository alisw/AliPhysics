/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////
//Author : Indranil Das, SINP, INDIA
//         Sukalyan Chattopadhyay, SINP, INDIA
//         
//
//Email :  indra.das@saha.ac.in
//         sukalyan.chattopadhyay@saha.ac.in 
//
// This class implements a hit reconstruction algorithm for the dimuon
// high level trigger.
// The algorithm finds 3 pad clusters by looking for unique pads with a charge
// above a certain threshold. A centre of gravity type calculation is applied
// to the three pads forming the cluster to find the hit's X or Y coordinate
// along the non-bending and bending planes individually.
// The sepperate X and Y coordinates are then merged to give the full coordinate
// of the hit.
/////////////////////////////////////////////////

#include <strings.h>
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTMUONRecHitsBlockStruct.h"


const int AliHLTMUONHitReconstructor::fgkDetectorId = 0xA00;
const int AliHLTMUONHitReconstructor::fgkDDLOffSet = 12 ;
const int AliHLTMUONHitReconstructor::fgkNofDDL = 8 ;

const int AliHLTMUONHitReconstructor::fgkDDLHeaderSize = 8;

const int AliHLTMUONHitReconstructor::fgkEvenLutSize = 3364287 + 1;
const int AliHLTMUONHitReconstructor::fgkOddLutSize = 1645631 + 1;

const int AliHLTMUONHitReconstructor::fgkLutLine[2] = {54208, 59648};

const int AliHLTMUONHitReconstructor::fgkMinIdManuChannel[2] = {64, 917696};
const int AliHLTMUONHitReconstructor::fgkMaxIdManuChannel[2] = {3364351,2563007};

const float AliHLTMUONHitReconstructor::fgkHalfPadSize[3] = {1.25, 2.50, 5.00};



//ClassImp(AliHLTMUONHitReconstructor)

AliHLTMUONHitReconstructor::AliHLTMUONHitReconstructor():
  fkBlockHeaderSize(8),
  fkDspHeaderSize(8),
  fkBuspatchHeaderSize(4),
  fDCCut(0),
  fPadData(NULL),
  fLookUpTableData(NULL),
  fRecPoints(NULL),
  fRecPointsCount(NULL),
  fMaxRecPointsCount(0),
  fCentralCountB(0),
  fCentralCountNB(0),
  fIdOffSet(0),
  fDDLId(0),
  fDigitPerDDL(0),
  fDetManuChannelIdList(NULL),
  fCentralChargeB(NULL),
  fCentralChargeNB(NULL),
  fRecX(NULL),
  fRecY(NULL),
  fAvgChargeX(NULL),
  fAvgChargeY(NULL),
  fNofFiredDetElem(0),
  fDebugLevel(0),
  fBusToDetElem()
{
  // ctor 
  
  if(AliHLTMUONHitReconstructor::fgkEvenLutSize > AliHLTMUONHitReconstructor::fgkOddLutSize){
    fPadData = new DHLTPad[AliHLTMUONHitReconstructor::fgkEvenLutSize];
  }
  else{
    fPadData = new DHLTPad[AliHLTMUONHitReconstructor::fgkOddLutSize];
  }


  fkBlockHeaderSize    = 8;
  fkDspHeaderSize      = 8;
  fkBuspatchHeaderSize = 4;

  bzero(fGetIdTotalData,336*80*2*sizeof(int));
}


AliHLTMUONHitReconstructor::~AliHLTMUONHitReconstructor()
{
  // dtor

  //printf("\nEnd of Run\n");

  delete []fPadData;
  delete []fLookUpTableData;

}

int AliHLTMUONHitReconstructor::GetLutLine(int iDDL) const { return ( iDDL<16 ) ? fgkLutLine[0] : fgkLutLine[1] ;}

bool AliHLTMUONHitReconstructor::LoadLookUpTable(DHLTLut* lookUpTableData, int lookUpTableId)
{
  // function that loads LookUpTable (= position of each pad with electronic channel associated with it)

  if(lookUpTableId<fgkDDLOffSet || lookUpTableId>= fgkDDLOffSet + fgkNofDDL){
    printf("DDL number is out of range (must be %d<=iDDL<%d)\n",fgkDDLOffSet,fgkDDLOffSet+fgkNofDDL);
    return false;
  }
  
  fDDLId = lookUpTableId;

  int lutSize = ((lookUpTableId%2)==0) ? fgkEvenLutSize : fgkOddLutSize ;
  int nofLutLine = GetLutLine(lookUpTableId);
  int idOffSet = fgkMinIdManuChannel[lookUpTableId%2];

  int detManuChannelId;

  fLookUpTableData = new DHLTLut[lutSize];

  fLookUpTableData[0].fIdManuChannel = 0;
  fLookUpTableData[0].fIX = 0 ;
  fLookUpTableData[0].fIY = 0 ;
  fLookUpTableData[0].fRealX = 0.0 ;
  fLookUpTableData[0].fRealY = 0.0 ;
  fLookUpTableData[0].fRealZ = 0.0 ;
  fLookUpTableData[0].fPlane = -1 ;
  fLookUpTableData[0].fPcbZone = -1 ;

  for(int i=0; i<nofLutLine; i++){

    detManuChannelId = lookUpTableData[i].fIdManuChannel - idOffSet + 1;
    fLookUpTableData[detManuChannelId].fIdManuChannel = lookUpTableData[i].fIdManuChannel - idOffSet;
    fLookUpTableData[detManuChannelId].fIX = lookUpTableData[i].fIX ;
    fLookUpTableData[detManuChannelId].fIY = lookUpTableData[i].fIY ;
    fLookUpTableData[detManuChannelId].fRealX = lookUpTableData[i].fRealX ;
    fLookUpTableData[detManuChannelId].fRealY = lookUpTableData[i].fRealY ;
    fLookUpTableData[detManuChannelId].fRealZ = lookUpTableData[i].fRealZ ;
    fLookUpTableData[detManuChannelId].fPcbZone = lookUpTableData[i].fPcbZone ;
    fLookUpTableData[detManuChannelId].fPlane = lookUpTableData[i].fPlane ;
  }
  return true;
}

bool AliHLTMUONHitReconstructor::SetBusToDetMap(BusToDetElem busToDetElem)
{

  // function that loads BusPatch To Detection Element (SlatId) map

  if(busToDetElem.size()==0)
    return false;
  else
    fBusToDetElem = busToDetElem;
  
  return true;
}


bool AliHLTMUONHitReconstructor::Run(int* rawData, int *rawDataSize, AliHLTMUONRecHitStruct recHit[], int *nofHit) 
{  
  // main function called by HLTReconstructor to perform DHLT Hitreconstruction 

  fRecPoints = &recHit[0];
  fMaxRecPointsCount = *nofHit;
  fRecPointsCount = nofHit;
  *fRecPointsCount = 0;

  fPadData[0].fDetElemId = 0;
  fPadData[0].fBuspatchId = 0;
  fPadData[0].fIdManuChannel = 0;
  fPadData[0].fIX = 0 ;
  fPadData[0].fIY = 0 ;
  fPadData[0].fRealX = 0.0 ;
  fPadData[0].fRealY = 0.0 ;
  fPadData[0].fRealZ = 0.0 ;
  fPadData[0].fPlane = -1 ;
  fPadData[0].fPcbZone = -1 ;
  fPadData[0].fCharge = 0 ;
  
  if(!ReadDDL(rawData,rawDataSize)){
    printf("Failed to read the complete DDL file\n");
    return false;
  }

  if(!FindRecHits()){
    printf("Failed to generate RecHits\n");
    return false;
  }
    
  return true;
}



bool AliHLTMUONHitReconstructor::ReadDDL(int* rawData, int *rawDataSize)
{
  //function to read Raw Data files

  int ddlRawDataSize;
  ddlRawDataSize = *rawDataSize;

  int *buffer = new int[ddlRawDataSize]; 
  buffer = (int *)rawData; 

  fIdOffSet= fgkMinIdManuChannel[(fDDLId%2)];
  fDetManuChannelIdList = new int[ddlRawDataSize];

  int index = 0;
  int dataCount = 0;
  fNofFiredDetElem = 0;
  int detElemId = 0 ;
  int buspatchId = 0;
  int prevDetElemId = 0 ;
  int totalBlockSize,blockRawDataSize;
  int totalDspSize,dspRawDataSize;
  int totalBuspatchSize,buspatchRawDataSize;
  int indexDsp,indexBuspatch,indexRawData;
  unsigned int dataWord;
  int charge;
  int idManuChannel;
  
  for(int iBlock = 0; iBlock < 2 ;iBlock++){  // loop over 2 blocks
    totalBlockSize = buffer[index + 1];
    blockRawDataSize = buffer[index + 2];
    indexDsp = index + fkBlockHeaderSize;
    while(blockRawDataSize > 0){
      totalDspSize = buffer[indexDsp + 1];
      dspRawDataSize = buffer[indexDsp + 2];
      dspRawDataSize --;                              // temporary solution to read buspatches 
      indexBuspatch = indexDsp + fkDspHeaderSize + 2; // this extra 2 word comes from the faulty defination of Dsp header size
      while(dspRawDataSize > 0){
	totalBuspatchSize = buffer[indexBuspatch + 1];
	buspatchRawDataSize = buffer[indexBuspatch + 2];
	buspatchId = buffer[indexBuspatch + 3];
	if((detElemId = fBusToDetElem[buspatchId])==0){
	  printf("No Detection element found for buspatch : %d\n",buspatchId);
	  return false;
	}
	indexRawData = indexBuspatch + fkBuspatchHeaderSize;
	while(buspatchRawDataSize > 0){
	  dataWord = buffer[indexRawData];
	  charge = (unsigned short)(dataWord & 0xFFF);
	  
	  idManuChannel = 0x0;
	  idManuChannel = (idManuChannel|(detElemId%100))<<17;
	  idManuChannel |= (dataWord >> 12) & 0x1FFFF;
	  idManuChannel -= fIdOffSet ;
	  
	  if(charge > fDCCut && charge > 5){  // (charge > 4) is due cut out the noise level  			
	    fPadData[idManuChannel].fBuspatchId = buspatchId;
	    fPadData[idManuChannel].fDetElemId = detElemId;
	    fPadData[idManuChannel].fIdManuChannel = idManuChannel;
	    fPadData[idManuChannel].fIX = fLookUpTableData[idManuChannel+1].fIX;
	    fPadData[idManuChannel].fIY = fLookUpTableData[idManuChannel+1].fIY;
	    fPadData[idManuChannel].fRealX = fLookUpTableData[idManuChannel+1].fRealX;
	    fPadData[idManuChannel].fRealY = fLookUpTableData[idManuChannel+1].fRealY;
	    fPadData[idManuChannel].fRealZ = fLookUpTableData[idManuChannel+1].fRealZ;
	    fPadData[idManuChannel].fPcbZone = fLookUpTableData[idManuChannel+1].fPcbZone;
	    fPadData[idManuChannel].fPlane = fLookUpTableData[idManuChannel+1].fPlane;
	    fPadData[idManuChannel].fCharge = charge;
	    
	    fDetManuChannelIdList[dataCount] = idManuChannel;
	    if(detElemId != prevDetElemId){
	      if(fNofFiredDetElem>0){
		fMaxFiredPerDetElem[fNofFiredDetElem-1] = dataCount;
	      }
	      fNofFiredDetElem++;
	      prevDetElemId = detElemId ;
		}
	    dataCount ++;
	  }
	  

	  indexRawData++;
	  buspatchRawDataSize --;
	}
	indexBuspatch += totalBuspatchSize;
	dspRawDataSize -= totalBuspatchSize;
      }// buspatch loop
      indexDsp += totalDspSize;
      blockRawDataSize -= totalDspSize;
    }// DSP loop
    index = totalBlockSize;
  }// Block loop
  
  delete[] buffer;
  
  fDigitPerDDL = dataCount;
  fMaxFiredPerDetElem[fNofFiredDetElem-1] = dataCount;
  
  
  return true;

}

bool AliHLTMUONHitReconstructor::FindRecHits() 
{
  // fuction that calls hit reconstruction detector element-wise   

  for(int iDet=0; iDet<fNofFiredDetElem ; iDet++){
    
    fCentralCountB = 0 ;
    fCentralCountNB = 0 ;
    fCentralChargeB = new int[fMaxFiredPerDetElem[iDet]];
    fCentralChargeNB = new int[fMaxFiredPerDetElem[iDet]];
    
    if(iDet>0)
      FindCentralHits(fMaxFiredPerDetElem[iDet-1],fMaxFiredPerDetElem[iDet]);
    else
      FindCentralHits(0,fMaxFiredPerDetElem[iDet]);
    
    RecXRecY();
    if(!MergeRecHits()){
      printf("Failed to merge hits\n");
      return false;
    }
    
    if(iDet==0)
      for(int i=0;i<fMaxFiredPerDetElem[iDet];i++)
	fGetIdTotalData[fPadData[fDetManuChannelIdList[i]].fIX][fPadData[fDetManuChannelIdList[i]].fIY][fPadData[fDetManuChannelIdList[i]].fPlane] = 0;
    else
      for(int i=fMaxFiredPerDetElem[iDet-1];i<fMaxFiredPerDetElem[iDet];i++)
	fGetIdTotalData[fPadData[fDetManuChannelIdList[i]].fIX][fPadData[fDetManuChannelIdList[i]].fIY][fPadData[fDetManuChannelIdList[i]].fPlane] = 0;

    //fDHLTTree->Fill();

    delete []fCentralChargeB;
    delete []fCentralChargeNB;

  }
    
  //for(int iPad=fDataPerDetElem[i];iPad<fDataPerDetElem[i+1];iPad++){
  for(int iPad=0;iPad<fDigitPerDDL;iPad++){
    fGetIdTotalData[fPadData[fDetManuChannelIdList[iPad]].fIX][fPadData[fDetManuChannelIdList[iPad]].fIY][fPadData[fDetManuChannelIdList[iPad]].fPlane] = 0;
    fPadData[fDetManuChannelIdList[iPad]].fDetElemId = 0;
    fPadData[fDetManuChannelIdList[iPad]].fBuspatchId = 0;
    fPadData[fDetManuChannelIdList[iPad]].fIdManuChannel = 0;
    fPadData[fDetManuChannelIdList[iPad]].fIX = 0 ;
    fPadData[fDetManuChannelIdList[iPad]].fIY = 0 ;
    fPadData[fDetManuChannelIdList[iPad]].fRealX = 0.0 ;
    fPadData[fDetManuChannelIdList[iPad]].fRealY = 0.0 ;
    fPadData[fDetManuChannelIdList[iPad]].fRealZ = 0.0 ;
    fPadData[fDetManuChannelIdList[iPad]].fPlane = -1 ;
    fPadData[fDetManuChannelIdList[iPad]].fPcbZone = -1 ;
    fPadData[fDetManuChannelIdList[iPad]].fCharge = 0 ;
  }  
  
  for(int i=0;i<13;i++)
    fMaxFiredPerDetElem[i] = 0;
  delete []fDetManuChannelIdList;

  return true;
}

void AliHLTMUONHitReconstructor::FindCentralHits(int minPadId, int maxPadId)
{
  // to find central hit associated with each cluster

  int b,nb;
  int idManuChannelCentral;
  bool hasFind;
  int idManuChannel;
  
  for(int iPad=minPadId;iPad<maxPadId;iPad++){
    idManuChannel   = fDetManuChannelIdList[iPad];
    
    
    fGetIdTotalData[fPadData[idManuChannel].fIX]
      [fPadData[idManuChannel].fIY]
      [fPadData[idManuChannel].fPlane] = idManuChannel;
    
    if(fPadData[idManuChannel].fPlane == 0 ){//&& fPadData[idManuChannel].fIY > (0+1) && fPadData[idManuChannel].fIY < (79 - 1)){
      //if(fPadData[idManuChannel].fIY > 0){
      if(fCentralCountB>0){
	hasFind = false;
	for(b = 0;b<fCentralCountB;b++){
	  idManuChannelCentral = fCentralChargeB[b];
	  if(fPadData[idManuChannel].fIX == fPadData[idManuChannelCentral].fIX
	     &&
	     (fPadData[idManuChannel].fIY 
	      == fPadData[idManuChannelCentral].fIY + 1 
	      ||
	      fPadData[idManuChannel].fIY 
	      == fPadData[idManuChannelCentral].fIY + 2 
	      ||
	      fPadData[idManuChannel].fIY 
	      == fPadData[idManuChannelCentral].fIY - 2 
	      ||
	      fPadData[idManuChannel].fIY 
	      == fPadData[idManuChannelCentral].fIY - 1)){
	    
	    hasFind = true;
	    if(fPadData[idManuChannel].fCharge > fPadData[idManuChannelCentral].fCharge){
	      fCentralChargeB[b] = idManuChannel;
	    }// if condn on pad charge
	  }// if condon on pad position
	}// for loop over b
	if(!hasFind){
	  fCentralChargeB[fCentralCountB] = idManuChannel;
	  fCentralCountB++;
	}
      }
      else{
	fCentralChargeB[fCentralCountB] = idManuChannel;
	fCentralCountB++;
      }// check the size of centralHitB
      for(b = 0;b<fCentralCountB;b++){
	idManuChannelCentral = fCentralChargeB[b];
      }
      //}// if cond on iY > 2 (to avoid edge value pb)
    }// B/Nb checking
    else{
      if(fCentralCountNB>0){
	hasFind = false;
	for(nb = 0;nb<fCentralCountNB;nb++){
	  idManuChannelCentral = fCentralChargeNB[nb];
	  if(fPadData[idManuChannel].fIY == fPadData[idManuChannelCentral].fIY
	     &&
	     (fPadData[idManuChannel].fIX 
	      == fPadData[idManuChannelCentral].fIX + 1 
	      ||
	      fPadData[idManuChannel].fIX
	      == fPadData[idManuChannelCentral].fIX + 2
	      ||
	      fPadData[idManuChannel].fIX
	      == fPadData[idManuChannelCentral].fIX - 2
	      ||
	      fPadData[idManuChannel].fIX
	      == fPadData[idManuChannelCentral].fIX - 1)){
	    
	    hasFind = true;	  
	    if(fPadData[idManuChannel].fCharge > fPadData[idManuChannelCentral].fCharge){
	      fCentralChargeNB[nb] = idManuChannel;
	    }// if condn over to find higher charge
	  }// if condn over to find position
	}// for loop over presently all nb values
	if(!hasFind){
	  fCentralChargeNB[fCentralCountNB] = idManuChannel;
	  fCentralCountNB++;
	}
      }// centralHitNB size test
      else{
	fCentralChargeNB[fCentralCountNB] = idManuChannel;
	fCentralCountNB++;
      }// centralHitNB size test
      
    }// fill for bending and nonbending hit
  }// detElemId loop


}

void AliHLTMUONHitReconstructor::RecXRecY()
{
  // find reconstructed X and Y for each plane separately
  int b,nb;
  int idCentral;
  int idLower = 0;
  int idUpper = 0;
  int idRight = 0;
  int idLeft = 0;
  fRecY = new float[fCentralCountB];
  fRecX = new float[fCentralCountNB];
  
  fAvgChargeY = new float[fCentralCountB];
  fAvgChargeX = new float[fCentralCountNB];

  for(b=0;b<fCentralCountB;b++){
    idCentral = fCentralChargeB[b];
    
    if(fPadData[idCentral].fIY==0)
      idLower = 0;
    else
      idLower = fGetIdTotalData[fPadData[idCentral].fIX][fPadData[idCentral].fIY-1][0];
    
    if(fPadData[idCentral].fIY==79)
      idUpper = 0;
    else
      idUpper = fGetIdTotalData[fPadData[idCentral].fIX][fPadData[idCentral].fIY+1][0];

    fRecY[b] = (fPadData[idCentral].fRealY*fPadData[idCentral].fCharge
	       +
	       fPadData[idUpper].fRealY*fPadData[idUpper].fCharge
	       +
	       fPadData[idLower].fRealY*fPadData[idLower].fCharge
	       )/(fPadData[idCentral].fCharge + fPadData[idUpper].fCharge + fPadData[idLower].fCharge) ;

    fAvgChargeY[b] = fPadData[idCentral].fCharge;
  

    fRecY[b] += 0.025*sin(12.56637*(0.25-(fRecY[b] - fPadData[idCentral].fRealY))) ;
  }
      
  for(nb=0;nb<fCentralCountNB;nb++){
    idCentral = fCentralChargeNB[nb];

    if(fPadData[idCentral].fIX==0)
      idLeft = 0;
    else
      idLeft = fGetIdTotalData[fPadData[idCentral].fIX-1][fPadData[idCentral].fIY][1];
    
    if(fPadData[idCentral].fIX==335)
      idRight = 0 ;
    else
      idRight = fGetIdTotalData[fPadData[idCentral].fIX+1][fPadData[idCentral].fIY][1];

    fRecX[nb] = (fPadData[idCentral].fRealX*fPadData[idCentral].fCharge
		 +
		 fPadData[idRight].fRealX*fPadData[idRight].fCharge
		 +
		 fPadData[idLeft].fRealX*fPadData[idLeft].fCharge
		 )/(fPadData[idCentral].fCharge + fPadData[idRight].fCharge + fPadData[idLeft].fCharge);
    

    fAvgChargeX[nb] = fPadData[idCentral].fCharge;
    
  }

}

bool AliHLTMUONHitReconstructor::MergeRecHits()
{
  // Merge reconstructed hits first over same plane then bending plane with non-bending plane

  int idCentralB,idCentralNB ;
  float padCenterXB;
  float padCenterYNB;
  float diffX,diffY;
  float halfPadLengthX,halfPadLengthY;

  // MERGE Bending Plane hits, which are placed side by side
  for(int i=0;i<fCentralCountB-1;i++){
    if(fRecY[i] != 0.0){
      for(int j=i+1;j<fCentralCountB;j++){
		 
	if(fCentralChargeB[i]==fCentralChargeB[j]){
	  fRecY[j] = 0.0;
	  continue;
	}
	else if(
	   (
	    fPadData[fCentralChargeB[i]].fIY == fPadData[fCentralChargeB[j]].fIY
	    )
 	   &&
	   (
	    fPadData[fCentralChargeB[i]].fIX == fPadData[fCentralChargeB[j]].fIX + 1
	    ||
	    fPadData[fCentralChargeB[i]].fIX == fPadData[fCentralChargeB[j]].fIX - 1
	    )
 	   &&
 	   fRecY[j] != 0.0
	   &&
	   fRecY[i] != 0.0
	   ){

	  if(fAvgChargeY[i] > fAvgChargeY[j]){
	    fRecY[i] = (fRecY[i]*fAvgChargeY[i] + fRecY[j]*fAvgChargeY[j]
			)/(fAvgChargeY[i] + fAvgChargeY[j]);
	    fRecY[j] = 0.0;
	  }
	  else{
	    fRecY[j] = (fRecY[i]*fAvgChargeY[i] + fRecY[j]*fAvgChargeY[j]
			)/(fAvgChargeY[i] + fAvgChargeY[j]);
	    fRecY[i] = 0.0;

	  }// search for higher charge
	}//pad position
      }//j for loop
    }//if fRecY[i] != 0.0
  }// i for loop


  
  // MERGE Non Bending Plane hits, which are place side by side
  for(int i=0;i<fCentralCountNB-1;i++){
    if(fRecX[i] != 0.0){
      for(int j=i+1;j<fCentralCountNB;j++){

	if(fCentralChargeNB[i]==fCentralChargeNB[j]){
	  fRecX[j] = 0.0;
	  continue;
	}
	else if(
	   (
	    fPadData[fCentralChargeNB[i]].fIX == fPadData[fCentralChargeNB[j]].fIX
	    )
	   &&
	   (
	    fPadData[fCentralChargeNB[i]].fIY == fPadData[fCentralChargeNB[j]].fIY + 1
	    ||
	    fPadData[fCentralChargeNB[i]].fIY == fPadData[fCentralChargeNB[j]].fIY - 1
	    )
	   &&
	   fRecX[j] != 0.0
	   &&
	   fRecX[i] != 0.0
	   ){

	  if(fAvgChargeX[i] > fAvgChargeX[j]){
	    fRecX[i] = (fRecX[i]*fAvgChargeX[i] + fRecX[j]*fAvgChargeX[j]
		       )/(fAvgChargeX[i] + fAvgChargeX[j]);
	    fRecX[j] = 0.0;
	  }
	  else{
	    fRecX[j] = (fRecX[i]*fAvgChargeX[i] + fRecX[j]*fAvgChargeX[j]
		       )/(fAvgChargeX[i] + fAvgChargeX[j]);
	    fRecX[i] = 0.0;
	  }// search for higher charge
	}//pad position
      }//j for loop
    }//if fRecX[i] != 0.0
  }// i for loop



  // Merge bending Plane hits with Non Bending
  

  for(int b=0;b<fCentralCountB;b++){
    if(fRecY[b]!=0.0){
      idCentralB = fCentralChargeB[b];
      padCenterXB = fPadData[idCentralB].fRealX; 
      
      halfPadLengthX = fgkHalfPadSize[fPadData[idCentralB].fPcbZone] ;

      for(int nb=0;nb<fCentralCountNB;nb++){
	if(fRecX[nb]!=0.0){
	  idCentralNB = fCentralChargeNB[nb];

	  padCenterYNB = fPadData[idCentralNB].fRealY;

	  halfPadLengthY = fgkHalfPadSize[fPadData[idCentralNB].fPcbZone] ;

	  if(fabsf(fRecX[nb]) > fabsf(padCenterXB))
	    diffX = fabsf(fRecX[nb]) -  fabsf(padCenterXB);
	  else
	    diffX = fabsf(padCenterXB) -  fabsf(fRecX[nb]);
	  
	  if(fabsf(padCenterYNB)>fabsf(fRecY[b]))
	    diffY = fabsf(padCenterYNB) - fabsf(fRecY[b]);
	  else
	    diffY =  fabsf(fRecY[b]) - fabsf(padCenterYNB);

	  if(diffX < halfPadLengthX && diffY < halfPadLengthY ){//&& fPadData[idCentralB].fIY != 0){

	    fRecPoints[(*fRecPointsCount)].fX = fRecX[nb];
	    fRecPoints[(*fRecPointsCount)].fY = fRecY[b];
	    fRecPoints[(*fRecPointsCount)].fZ = fPadData[idCentralB].fRealZ;
//	    fRecPoints[(*fRecPointsCount)].fDetElemId = (AliHLTUInt32_t)fPadData[idCentralB].fDetElemId;
	    (*fRecPointsCount)++;
	    if((*fRecPointsCount) == fMaxRecPointsCount){
	      printf("Nof RecHit (i.e. %d) exceeds the max nof RecHit limit %d\n",(*fRecPointsCount),fMaxRecPointsCount);
	      return false;
	    }
	  }//if lies wihtin 5.0 mm
	}// condn over fRecX ! = 0.0
      }// loop over NB side
    }// condn on fRecY[b] !=  0.0
  }// loop over B side;

  delete []fRecX;
  delete []fRecY;
  
  delete []fAvgChargeX;
  delete []fAvgChargeY;

  return true;
}

