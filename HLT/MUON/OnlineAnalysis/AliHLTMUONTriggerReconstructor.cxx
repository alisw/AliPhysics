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

/* $Id$ */

/**********************************************************************
 Created on : 16/05/2007
 Purpose    : This class is supposed to read the tracker DDL files and 
              give the output AliMUONCoreTriggerRecord
 Author     : Indranil Das, HEP Division, SINP
 Email      : indra.das@saha.ac.in | indra.ehep@gmail.com
**********************************************************************/

///*
//
//  The trigger reconstructor class is designed to deal the rawdata inputfiles
//  to findout the the reconstructed hits at the trigger DDL. The output is send
//  to the output block for further processing.
//
//  Author : Indranil Das ( indra.das@saha.ac.in || indra.ehep@gmail.com )
// 
//*/

#include "AliHLTMUONTriggerReconstructor.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONCalculations.h"
#include <vector>
#include <cassert>

const int AliHLTMUONTriggerReconstructor::fgkDetectorId = 0xB00;
const int AliHLTMUONTriggerReconstructor::fgkDDLOffSet = 20 ;
const int AliHLTMUONTriggerReconstructor::fgkNofDDL = 2 ;

const int AliHLTMUONTriggerReconstructor::fgkDDLHeaderSize = 8;
const int AliHLTMUONTriggerReconstructor::fgkEvenLutSize = 2602351 + 1; 
const int AliHLTMUONTriggerReconstructor::fgkOddLutSize = 2528735 + 1;

const int AliHLTMUONTriggerReconstructor::fgkLutLine = 10496;

const int AliHLTMUONTriggerReconstructor::fgkMinIdManuChannel[2] = {819616, 862288};
const int AliHLTMUONTriggerReconstructor::fgkMaxIdManuChannel[2] = {3421966, 3391022};

const float AliHLTMUONTriggerReconstructor::fgkHalfPadSizeXB[3] = {8.5, 17.0, 25.5};
const float AliHLTMUONTriggerReconstructor::fgkHalfPadSizeYNB[2] = {25.5, 34.0};

const int AliHLTMUONTriggerReconstructor::fgkDetElem = 9*4 ; // 9 detele per half chamber


AliHLTMUONTriggerReconstructor::AliHLTMUONTriggerReconstructor()
  :
  fPadData(NULL),
  fLookUpTableData(NULL),
  fRecPoints(NULL),
  fRecPointsCount(NULL),
  fMaxRecPointsCount(0),
  fMaxFiredPerDetElem(),
  fDetElemToDataId(),
  fDDLId(0),
  fIdOffSet(0),
  fTrigRecId(0)
{
  // ctor 
  
  if(AliHLTMUONTriggerReconstructor::fgkEvenLutSize > AliHLTMUONTriggerReconstructor::fgkOddLutSize){
    fPadData = new AliHLTMUONHitReconstructor::DHLTPad[AliHLTMUONTriggerReconstructor::fgkEvenLutSize];
  }
  else{
    fPadData = new AliHLTMUONHitReconstructor::DHLTPad[AliHLTMUONTriggerReconstructor::fgkOddLutSize];
  }

  bzero(fGetIdTotalData,104*64*2*sizeof(int));
}


AliHLTMUONTriggerReconstructor::~AliHLTMUONTriggerReconstructor()
{
  // dtor
  delete []fPadData;
  delete []fLookUpTableData;
}

bool AliHLTMUONTriggerReconstructor::SetRegToLocCardMap(RegToLoc* regToLoc)
{
  if(!memcpy(fRegToLocCard,regToLoc,128*sizeof(RegToLoc)))
    return false;

  for(int i=0;i<128;i++){
    HLTDebug("DDL : %d, reg : %d, loc : %d",fRegToLocCard[i].fTrigDDL,
	    fRegToLocCard[i].fRegId,fRegToLocCard[i].fLocId);
  }

  return true;
}

bool AliHLTMUONTriggerReconstructor::LoadLookUpTable(AliHLTMUONHitReconstructor::DHLTLut* lookUpTableData, int lookUpTableId)
{
  if(lookUpTableId<fgkDDLOffSet || lookUpTableId>= fgkDDLOffSet + fgkNofDDL){
    HLTError("DDL number is out of range (must be %d<=iDDL<%d)",fgkDDLOffSet,fgkDDLOffSet+fgkNofDDL);
    return false;
  }
  
  fDDLId = lookUpTableId;

  int lutSize = ((lookUpTableId%2)==0) ? fgkEvenLutSize : fgkOddLutSize ;
  int nofLutLine = fgkLutLine ;
  fIdOffSet = fgkMinIdManuChannel[lookUpTableId%2];

  int detManuChannelId;

  fLookUpTableData = new AliHLTMUONHitReconstructor::DHLTLut[lutSize];

  memset(fLookUpTableData,-1,lutSize*sizeof(AliHLTMUONHitReconstructor::DHLTLut));

  for(int i=0; i<nofLutLine; i++){

    detManuChannelId = lookUpTableData[i].fIdManuChannel - fIdOffSet + 1;
    fLookUpTableData[detManuChannelId].fIdManuChannel = lookUpTableData[i].fIdManuChannel - fIdOffSet;
    fLookUpTableData[detManuChannelId].fIX = lookUpTableData[i].fIX ;
    fLookUpTableData[detManuChannelId].fIY = lookUpTableData[i].fIY ;
    fLookUpTableData[detManuChannelId].fRealX = lookUpTableData[i].fRealX ;
    fLookUpTableData[detManuChannelId].fRealY = lookUpTableData[i].fRealY ;
    fLookUpTableData[detManuChannelId].fRealZ = lookUpTableData[i].fRealZ ;
    fLookUpTableData[detManuChannelId].fPcbZone = lookUpTableData[i].fPcbZone ;
    fLookUpTableData[detManuChannelId].fPlane = lookUpTableData[i].fPlane ;
  }
  return true;

  return true;
}

bool AliHLTMUONTriggerReconstructor::Run(
		const AliHLTUInt32_t* rawData,
		AliHLTUInt32_t rawDataSize,
		AliHLTMUONTriggerRecordStruct* trigRecord,
		AliHLTUInt32_t& nofTrigRec
	)
{

  fRecPoints = trigRecord;
  fMaxRecPointsCount = nofTrigRec;
  fRecPointsCount = &nofTrigRec;
  *fRecPointsCount = 0;
  fMaxFiredPerDetElem.clear();
  fDetElemToDataId.clear();

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
  
  if(!ReadDDL(rawData, rawDataSize)){
    HLTError("Failed to read the complete DDL file\n");
    return false;
  }

  if(!FindTrigHits()){
    HLTError("Failed to generate RecHits\n");
    return false;
  }
    
  return true;
}


bool AliHLTMUONTriggerReconstructor::ReadDDL(
		const AliHLTUInt32_t* rawData,
		AliHLTUInt32_t /*rawDataSize*/
	)
{

  int idManuChannel ;
  
  int index = 0;
  int dataCount = 0;
  int detElemId = 0 ;
  int reg_output,reg_phys_trig_occur ;
  int iLocIndex,loc,locDec,triggY,sign,loDev,triggX;
  int iRegLoc, locId ;
  short pattern[2][4]; // 2 stands for two cathode planes and 4 stands for 4 chambers

  Int_t offset,ithSwitch,secondLocation,idetElemId;

  int shiftIndex = 10 - 6 - 1; // the one comes due to indexing from zero

  DataIdIndex dataIndex;
#ifdef DEBUG
  int globalcard_data_occurance = (rawData[index]>>10)&0x1; //Set to 1 if global info present in DDL else set to 0 
  int version = (rawData[index]>>12)&0xFF; // software version
  int serial_number =  (rawData[index]>>20)&0xF; // serial number set to 0xF 
#endif
  int phys_trig_occur = (rawData[index]>>30)&0x1; // 1 for physics trigger, 0 for software trigger
  
  HLTDebug("globalcard_data_occurance  %d, version  %d, serial_number  %d, phys_trig_occur  %d",
	 globalcard_data_occurance,version,serial_number,phys_trig_occur);

  if(!phys_trig_occur) // for software trigger
    index += 8 ;// corresponding to scalar words
  
  index += 1 ; // To skip the separator 0xDEADFACE
  
  index += 4 ; // corresponding to global input
  
  index += 1 ; // reaches to global output
  
  if((fDDLId - AliHLTMUONTriggerReconstructor::fgkDDLOffSet) == 0){ //if globalData is present in DDL 0 (presummed may be changed)
#ifdef DEBUG
    int singleLpt = rawData[index] & 0x1;
    int singleHpt = (rawData[index] >> 1) & 0x1;
    
    int pairUnlikeLpt = (rawData[index] >> 4)  & 0x1;
    int pairUnlikeHpt = (rawData[index] >> 5)  & 0x1;
    
    int pairLikeLpt = (rawData[index] >> 2)  & 0x1;
    int pairLikeHpt = (rawData[index] >> 3)  & 0x1;
#endif
    HLTDebug("singleLpt : %x, singleHpt : %x, pairUnlikeLpt : %x, pairUnlikeHpt : %x, pairLikeLpt : %x, pairLikeHpt : %x",
	     singleLpt,singleHpt,pairUnlikeLpt,pairUnlikeHpt,pairLikeLpt,pairLikeHpt);
  }

  if(!phys_trig_occur)
    index += 10 ;// corresponds to scalar words

  index += 1; // separator 0xDEADBEEF 

  for (int iReg = 0; iReg < 8; iReg++) {
    index += 1; // DARC Status Word
    index += 1; // Regeional Word
    reg_output = rawData[index] & 0xFF;
    reg_phys_trig_occur = ( rawData[index] >> 31) & 0x1;
    
    index += 2; // 2 words for regional input
    
    index += 1; // L0 counter
    
    if(!reg_phys_trig_occur)
      index += 10;
    
    index += 1 ; // end of Regeonal header 0xBEEFFACE

    for(int iLoc = 0; iLoc < 16 ; iLoc++){

      iLocIndex = index ;      

      loc = (rawData[index+5] >> 19) &  0xF ;
      
      locDec = (rawData[index+5] >> 15) & 0xF;
      triggY = (rawData[index+5] >> 14) & 0x1 ;
      sign = (rawData[index+5] >> 9) & 0x1;
      loDev = (rawData[index+5] >> 5) & 0xF ;
      triggX = (loDev >> 4 & 0x1 ) && !(loDev & 0xF) ;

      if( locDec != 0x9 ){ // check for Dec
	
	iRegLoc = iReg*16 + iLoc;
	locId = fRegToLocCard[iRegLoc].fLocId ; 
	
	if(locId<=234){ // to avoid the copy locCards
	  
	  index += 1;
	  pattern[0][0] = rawData[index] & 0xFFFF; // x-strip pattern for chamber 0 
	  pattern[0][1] = (rawData[index] >> 16) & 0xFFFF; // x-strip pattern for chamber 1
 	  index += 1; 
	  pattern[0][2] = rawData[index] & 0xFFFF; 
	  pattern[0][3] = (rawData[index] >> 16) & 0xFFFF; 
	  
 	  index += 1;
	  pattern[1][0] = rawData[index] & 0xFFFF; // y-strip pattern for chamber 0
	  pattern[1][1] = (rawData[index] >> 16) & 0xFFFF; // y-strip pattern for chamber 0 
 	  index += 1; 
	  pattern[1][2] = rawData[index] & 0xFFFF; 
	  pattern[1][3] = (rawData[index] >> 16) & 0xFFFF; 
	  
	  if(pattern[0][0] || pattern[0][1] || pattern[0][2] || pattern[0][3]
	     || pattern[1][0] || pattern[1][1] || pattern[1][2] || pattern[1][3]
	     ){

	    HLTDebug("iReg: %d, iLoc :%d, locId : %d,X : %x, %x, %x, %x ...Y : %x, %x, %x, %x",
		    iReg,iLoc,locId,pattern[0][0],pattern[0][1],pattern[0][2],pattern[0][3],
		    pattern[1][0],pattern[1][1],pattern[1][2],pattern[1][3]);

	    for(int iChamber = 0; iChamber < 4 ; iChamber++){ //4 chambers per DDL 
	      for(int iPlane = 0; iPlane < 2 ; iPlane++){// 2 cathode plane
		if(pattern[iPlane][iChamber]){
		  detElemId = fRegToLocCard[iRegLoc].fDetElemId[iChamber];
		  HLTDebug("\tdetElemId : %d\n",detElemId);
		  for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy) {
		    if ((pattern[iPlane][iChamber] >> ibitxy) & 0x1) {
	  
		      // not quite sure about this
		      offset = 0;
		      ithSwitch = (fRegToLocCard[iRegLoc].fSwitch >> shiftIndex) & 0x1;
		      if (iPlane && ithSwitch) offset = -8;
		      
		      secondLocation = ibitxy + offset;
		      
		      idetElemId = detElemId%1000;

		      idetElemId &= 0x1FF ;
		      iPlane &= 0x1 ;
		      locId &= 0xFF ;
		      secondLocation &= 0xF ;
		      
		      idManuChannel &= 0x0;
		      idManuChannel = (idManuChannel|idetElemId)<<1;  
		      idManuChannel = (idManuChannel|iPlane)<<8;  
		      idManuChannel = (idManuChannel|locId)<<4 ;
		      idManuChannel |= secondLocation  ;

		      idManuChannel -= fIdOffSet ;
		      
		      if(fLookUpTableData[idManuChannel+1].fIdManuChannel == -1) //skip uninitialized values
			continue;

	   	      fPadData[idManuChannel].fDetElemId = detElemId;
	   	      fPadData[idManuChannel].fIdManuChannel = idManuChannel;
	   	      fPadData[idManuChannel].fIX = fLookUpTableData[idManuChannel+1].fIX;
	   	      fPadData[idManuChannel].fIY = fLookUpTableData[idManuChannel+1].fIY;
	   	      fPadData[idManuChannel].fRealX = fLookUpTableData[idManuChannel+1].fRealX;
	   	      fPadData[idManuChannel].fRealY = fLookUpTableData[idManuChannel+1].fRealY;
	   	      fPadData[idManuChannel].fRealZ = fLookUpTableData[idManuChannel+1].fRealZ;
	   	      fPadData[idManuChannel].fPcbZone = fLookUpTableData[idManuChannel+1].fPcbZone;
	   	      fPadData[idManuChannel].fPlane = fLookUpTableData[idManuChannel+1].fPlane;
	   	      
	   	      // !!!!!!!!!!!! THIS IS A COMPLETE HACK !!!!!!!!!
	   	      // Have to use fBuspatchId to store the chamber number,
	   	      // and fCharge as the trigger records ID number.
	   	      // This must be fixed in the future by using more appropriate
	   	      // structures.
	   	      fPadData[idManuChannel].fBuspatchId = iChamber;
	   	      fPadData[idManuChannel].fCharge = fTrigRecId;
	   	      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		      HLTDebug("\t Hit Found fo ich : %d, iPlane : %d, detelem %d, id : %d, at (%lf, %lf, %lf) cm"
			      ,iChamber,fLookUpTableData[idManuChannel+1].fPlane,detElemId,fLookUpTableData[idManuChannel+1].fIdManuChannel,
			      fPadData[idManuChannel].fRealX,
			      fPadData[idManuChannel].fRealY,fPadData[idManuChannel].fRealZ);
			      
		      if(fMaxFiredPerDetElem[detElemId] == 0){
			DataIdIndex first;
			first.push_back(idManuChannel);
			fDetElemToDataId[detElemId] = first;
		      }else{
			dataIndex =  fDetElemToDataId[detElemId];
			dataIndex.push_back(idManuChannel);
			fDetElemToDataId[detElemId] = dataIndex;
		      }

		      fMaxFiredPerDetElem[detElemId] = fMaxFiredPerDetElem[detElemId] + 1;

	   	      dataCount ++;
		      
		    }//pattern maching is found 
		  }// loop of ibitxy
		}// if pattern
	      }// iplane
	    }// ichamber
	    
		// Increment trigger record Id and keep it positive.
		//TODO: handle the wrapparound better.
		if (fTrigRecId < 0x7FFFFFFF)
			fTrigRecId++;
		else
			fTrigRecId = 0;
	    
	  }// if any non zero pattern found


	  index += 1 ; // skipping the last word though it is important
	  
	}// if locId <=234
      }// Dec Condn
	
      
      if(!reg_phys_trig_occur)
	 index += 45;
	
      index += 1; // end of local Data 0xCAFEFADE

      HLTDebug("iReg %d, iLoc %d, locId : %d, trigY %x, triggX %x, loDev %x, dec %x, sign %x,rawData : %x",
	       iReg,iLoc,locId,triggY,triggX,loDev,dec,sign, rawData[index]);

      index = iLocIndex + 6 ; //important to reset the index counter for fake locids like 235 
      
     }// iLoc loop
     
  }// iReg Loop

  return true;
}


// !!!!!!!!!!!! THIS IS A COMPLETE HACK !!!!!!!!!
// The following intermediate structures should not be necessary
// in a proper implementation. This will have to be fixed.
struct HitInfo
{
	AliHLTMUONRecHitStruct fHit;
	AliHLTInt32_t fChamber;
	AliHLTInt32_t fTrigRecId;
};

// This array stores the reconstructed hits with extra chamber numbers
// and trigger records numbers.
std::vector<HitInfo> hitInfo;
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


bool AliHLTMUONTriggerReconstructor::FindTrigHits() 
{

  map<int,DataIdIndex>::iterator it;

  hitInfo.clear();
  for(it = fDetElemToDataId.begin(); it != fDetElemToDataId.end(); it++){
    HLTDebug("Nof data found in Detelem : %d = %d",it->first,(it->second).size());
    if(!MergeTrigHits(it->second))
      return false;
  }// loop over detection element
  
// !!!!!!!!!!!! THIS IS A COMPLETE HACK !!!!!!!!!

	// We now have to go through the list of reconstructed hits
	// and collect all the ones that belong to a given trigger record.
	while (hitInfo.size() > 0)
	{
		// Must check if we have not overflowed the buffer.
		if((*fRecPointsCount) == fMaxRecPointsCount){
			HLTError("Nof RecHit (i.e. %d) exceeds the max nof RecHit limit %d\n",(*fRecPointsCount),fMaxRecPointsCount);
			return false;
		}
	
		bool hitset[4] = {false, false, false, false};
		
		// Choose the current trigger record we are collecting hits for.
		AliHLTInt32_t currentTrig = hitInfo[0].fTrigRecId;
		
		// Find the hits that belong to the current trigger record "currentTrig"
		for (int n = (int)hitInfo.size() - 1; n >= 0; n--)
		{
			if (hitInfo[n].fTrigRecId != currentTrig) continue;
			
			assert( 0 <= hitInfo[n].fChamber and hitInfo[n].fChamber < 4 );
			fRecPoints[(*fRecPointsCount)].fHit[hitInfo[n].fChamber] = hitInfo[n].fHit;
			hitset[hitInfo[n].fChamber] = true;
			hitInfo.erase(hitInfo.end());
		}
		
		// Make sure the hits that are not set, get set to a nil value.
		for (int i = 0; i < 4; i++)
		{
			if (not hitset[i])
			{
				fRecPoints[(*fRecPointsCount)].fHit[i]
					= AliHLTMUONConstants::NilRecHitStruct();
			}
		}
		
		fRecPoints[(*fRecPointsCount)].fId = currentTrig;
		
		// Calculate the momentum and fill in the flags and momentum fields.
		AliHLTMUONCalculations::ComputeMomentum(
				fRecPoints[(*fRecPointsCount)].fHit[0].fX,
				fRecPoints[(*fRecPointsCount)].fHit[0].fY,
				fRecPoints[(*fRecPointsCount)].fHit[2].fY,
				fRecPoints[(*fRecPointsCount)].fHit[0].fZ,
				fRecPoints[(*fRecPointsCount)].fHit[2].fZ
			);
		fRecPoints[(*fRecPointsCount)].fPx = AliHLTMUONCalculations::Px();
		fRecPoints[(*fRecPointsCount)].fPy = AliHLTMUONCalculations::Py();
		fRecPoints[(*fRecPointsCount)].fPz = AliHLTMUONCalculations::Pz();
		
		fRecPoints[(*fRecPointsCount)].fFlags =
			AliHLTMUONUtils::PackTriggerRecordFlags(
				AliHLTMUONCalculations::Sign(),
				hitset
			);
		
		// Increment the counter for the number of trigger records found.
		(*fRecPointsCount)++;
	}
  
	hitInfo.clear();
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // Reseting the fPadData to zero for next call to Run.
  DataIdIndex dataIndex;
  for(it = fDetElemToDataId.begin(); it != fDetElemToDataId.end(); it++){
    dataIndex = it->second;
    for(size_t i=0;i<dataIndex.size();i++){
      fPadData[dataIndex.at(i)].fDetElemId = 0;
      fPadData[dataIndex.at(i)].fBuspatchId = 0;
      fPadData[dataIndex.at(i)].fIdManuChannel = 0;
      fPadData[dataIndex.at(i)].fIX = 0 ;
      fPadData[dataIndex.at(i)].fIY = 0 ;
      fPadData[dataIndex.at(i)].fRealX = 0.0 ;
      fPadData[dataIndex.at(i)].fRealY = 0.0 ;
      fPadData[dataIndex.at(i)].fRealZ = 0.0 ;
      fPadData[dataIndex.at(i)].fPlane = -1 ;
      fPadData[dataIndex.at(i)].fPcbZone = -1 ;
      fPadData[dataIndex.at(i)].fCharge = 0 ;
    }// data per detelem loop  
  }//detelem loop
  
  return true;
}


bool AliHLTMUONTriggerReconstructor::MergeTrigHits(DataIdIndex& dataIndex)
{
  int idManuChannelB, idManuChannelNB;
  float halfPadLengthX,halfPadLengthY;
  float diffX,diffY;

  HLTDebug("\tThe bending plane hits are :");
  for(size_t iPad=0;iPad<dataIndex.size();iPad++){
    idManuChannelB   = dataIndex.at(iPad);
    if(fPadData[idManuChannelB].fPlane == 0){
      HLTDebug("\t detelem :%d, pcbzone : %d, (%f, %f, %f) cm",fPadData[idManuChannelB].fDetElemId,fPadData[idManuChannelB].fPcbZone,fPadData[idManuChannelB].fRealX,
	      fPadData[idManuChannelB].fRealY,fPadData[idManuChannelB].fRealZ);
    }
  }

  HLTDebug("\tThe non-bending plane hits are :");
  for(size_t jPad=0;jPad<dataIndex.size();jPad++){
    idManuChannelNB   = dataIndex.at(jPad);
    if(fPadData[idManuChannelNB].fPlane == 1){
      HLTDebug("\t detelem :%d, pcbzone : %d,(%f, %f, %f) cm",fPadData[idManuChannelNB].fDetElemId,fPadData[idManuChannelNB].fPcbZone,fPadData[idManuChannelNB].fRealX,
	      fPadData[idManuChannelNB].fRealY,fPadData[idManuChannelNB].fRealZ);
    }
  }

  for(size_t iPad=0;iPad<dataIndex.size();iPad++){
    idManuChannelB   = dataIndex.at(iPad);
    if(fPadData[idManuChannelB].fPlane == 0){
      
      halfPadLengthX = AliHLTMUONTriggerReconstructor::fgkHalfPadSizeXB[fPadData[idManuChannelB].fPcbZone] ;

      for(size_t jPad=0;jPad<dataIndex.size();jPad++){
	idManuChannelNB   = dataIndex.at(jPad);;
	if(fPadData[idManuChannelNB].fPlane == 1){
	  
	  halfPadLengthY = AliHLTMUONTriggerReconstructor::fgkHalfPadSizeYNB[fPadData[idManuChannelNB].fPcbZone] ;
	  
	  if(fabsf(fPadData[idManuChannelNB].fRealX) > fabsf(fPadData[idManuChannelB].fRealX))
	    diffX = fabsf(fPadData[idManuChannelNB].fRealX) - fabsf(fPadData[idManuChannelB].fRealX);
	  else
	    diffX = fabsf(fPadData[idManuChannelB].fRealX) - fabsf(fPadData[idManuChannelNB].fRealX) ;

	  
	  if(fabsf(fPadData[idManuChannelNB].fRealY) > fabsf(fPadData[idManuChannelB].fRealY))
	    diffY = fabsf(fPadData[idManuChannelNB].fRealY) - fabsf(fPadData[idManuChannelB].fRealY);
	  else
	    diffY = fabsf(fPadData[idManuChannelB].fRealY) - fabsf(fPadData[idManuChannelNB].fRealY) ;
 	  HLTDebug("\tdiffX %f,  halfPadLengthX %f,  diffY  %f, halfPadLengthY  %f\n",diffX,halfPadLengthX,diffY,halfPadLengthY);

 	  if(diffX < halfPadLengthX + 1.0 && diffY < halfPadLengthY + 1.0 ){// added redundancy of 1.0 cm due to the pb of geometrical segmentation 

	    HitInfo hit;
	    hit.fHit.fX = fPadData[idManuChannelNB].fRealX;
	    hit.fHit.fY = fPadData[idManuChannelB].fRealY;
	    hit.fHit.fZ = fPadData[idManuChannelNB].fRealZ;
	    // !!!!!!!!!!!! THIS IS A COMPLETE HACK !!!!!!!!!
	    // Remember we had to use fBuspatchId to store the chamber number,
	    // and fCharge for the trigger record's ID number.
	    hit.fChamber = fPadData[idManuChannelNB].fBuspatchId;
	    hit.fTrigRecId = fPadData[idManuChannelNB].fCharge;
	    hitInfo.push_back(hit);
	    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    //fRecPoints[(*fRecPointsCount)].fHit[0] = hit;
	    //fRecPoints[(*fRecPointsCount)].fId = fPadData[idManuChannelB].fDetElemId ;
	    
	    //(*fRecPointsCount)++;
	    //if((*fRecPointsCount) == fMaxRecPointsCount){
	    //  HLTFatal("Nof RecHit (i.e. %d) exceeds the max nof RecHit limit %d\n",(*fRecPointsCount),fMaxRecPointsCount);
	    //  return false;
	    //}

	    HLTDebug("\t\t\tdetelem : %d, x %f, y %f, z %f\n",fPadData[idManuChannelB].fDetElemId,fPadData[idManuChannelNB].fRealX,
		   fPadData[idManuChannelB].fRealY,fPadData[idManuChannelB].fRealZ);
	  }
	  
	}//condn for non-bending plane
      }//for loop for non-bending plane

    }// condn for bending plane
  }// for loop for bending plane
	
  return true;
}

