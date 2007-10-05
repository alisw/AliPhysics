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

const int AliHLTMUONTriggerReconstructor::fgkMinIdUnique[2] = {819616, 862288};
const int AliHLTMUONTriggerReconstructor::fgkMaxIdUnique[2] = {3421966, 3391022};

const float AliHLTMUONTriggerReconstructor::fgkHalfPadSizeXB[3] = {8.5, 17.0, 25.5};
const float AliHLTMUONTriggerReconstructor::fgkHalfPadSizeYNB[2] = {25.5, 34.0};

const int AliHLTMUONTriggerReconstructor::fgkDetElem = 9*4 ; // 9 detele per half chamber


AliHLTMUONTriggerReconstructor::AliHLTMUONTriggerReconstructor() :
	fLookUpTableData(NULL),
	fMaxRecPointsCount(0),
	fDDLId(0),
	fIdOffSet(0),
	fTrigRecId(0)
{
	// ctor
}


AliHLTMUONTriggerReconstructor::~AliHLTMUONTriggerReconstructor()
{
	// dtor
	if (fLookUpTableData != NULL)
		delete [] fLookUpTableData;
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
  fIdOffSet = fgkMinIdUnique[lookUpTableId%2];

  int detUniqueId;

  fLookUpTableData = new AliHLTMUONHitReconstructor::DHLTLut[lutSize];

  memset(fLookUpTableData,-1,lutSize*sizeof(AliHLTMUONHitReconstructor::DHLTLut));

  for(int i=0; i<nofLutLine; i++){

    detUniqueId = lookUpTableData[i].fIdManuChannel - fIdOffSet + 1;
    fLookUpTableData[detUniqueId].fIdManuChannel = lookUpTableData[i].fIdManuChannel - fIdOffSet;
    fLookUpTableData[detUniqueId].fIX = lookUpTableData[i].fIX ;
    fLookUpTableData[detUniqueId].fIY = lookUpTableData[i].fIY ;
    fLookUpTableData[detUniqueId].fRealX = lookUpTableData[i].fRealX ;
    fLookUpTableData[detUniqueId].fRealY = lookUpTableData[i].fRealY ;
    fLookUpTableData[detUniqueId].fRealZ = lookUpTableData[i].fRealZ ;
    fLookUpTableData[detUniqueId].fPcbZone = lookUpTableData[i].fPcbZone ;
    fLookUpTableData[detUniqueId].fPlane = lookUpTableData[i].fPlane ;
  }
  
  return true;
}

bool AliHLTMUONTriggerReconstructor::Run(
		const AliHLTUInt32_t* rawData,
		// TODO: if we are not checking rawDataSize then it means we are
		// not parsing the raw data safely or checking for corruption carefully.
		// This must be fixed at some point.
		AliHLTUInt32_t /*rawDataSize*/,
		AliHLTMUONTriggerRecordStruct* trigRecord,
		AliHLTUInt32_t& nofTrigRec,
		bool suppressPartialTrigs
	)
{
  fMaxRecPointsCount = nofTrigRec;
  
  // nofTrigRec now becomes the output of how many trigger records were found.
  nofTrigRec = 0;
  
  int lutAddress;
  
  int index = 0;
  int dataCount = 0;
  int detElemId = 0 ;
  int reg_output,reg_phys_trig_occur ;
  int iLocIndex,loc,locDec,triggY,sign,loDev,triggX;
  int iRegLoc = 0, locId = 0;
  short pattern[2][4]; // 2 stands for two cathode planes and 4 stands for 4 chambers

  Int_t offset,ithSwitch,secondLocation,idetElemId;

  int shiftIndex = 10 - 6 - 1; // the one comes due to indexing from zero

#ifdef __DEBUG
  int globalcard_data_occurance = (rawData[index]>>10)&0x1; //Set to 1 if global info present in DDL else set to 0 
  int version = (rawData[index]>>12)&0xFF; // software version
  int serial_number =  (rawData[index]>>20)&0xF; // serial number set to 0xF 
#endif
  int phys_trig_occur = (rawData[index]>>30)&0x1; // 1 for physics trigger, 0 for software trigger
  
  HLTDebug("globalcard_data_occurance  %d, version  %d, serial_number  %d, phys_trig_occur  %d",
	 globalcard_data_occurance,version,serial_number,phys_trig_occur
  );

  if(!phys_trig_occur) // for software trigger
    index += 8 ;// corresponding to scalar words
  
  index += 1 ; // To skip the separator 0xDEADFACE
  
  index += 4 ; // corresponding to global input
  
  index += 1 ; // reaches to global output
  
  if((fDDLId - AliHLTMUONTriggerReconstructor::fgkDDLOffSet) == 0){ //if globalData is present in DDL 0 (presummed may be changed)
#ifdef __DEBUG
    int singleLpt = rawData[index] & 0x1;
    int singleHpt = (rawData[index] >> 1) & 0x1;
    
    int pairUnlikeLpt = (rawData[index] >> 4)  & 0x1;
    int pairUnlikeHpt = (rawData[index] >> 5)  & 0x1;
    
    int pairLikeLpt = (rawData[index] >> 2)  & 0x1;
    int pairLikeHpt = (rawData[index] >> 3)  & 0x1;
#endif
    HLTDebug("singleLpt : %x, singleHpt : %x, pairUnlikeLpt : %x, pairUnlikeHpt : %x, pairLikeLpt : %x, pairLikeHpt : %x",
	     singleLpt,singleHpt,pairUnlikeLpt,pairUnlikeHpt,pairLikeLpt,pairLikeHpt
    );
  }

  if(!phys_trig_occur)
    index += 10; // corresponds to scalar words

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
    
    index += 1; // end of Regeonal header 0xBEEFFACE

    for(int iLoc = 0; iLoc < 16 ; iLoc++){

      iLocIndex = index;      

      loc = (rawData[index+5] >> 19) &  0xF ;
      
      locDec = (rawData[index+5] >> 15) & 0xF;
      triggY = (rawData[index+5] >> 14) & 0x1;
      sign = (rawData[index+5] >> 9) & 0x1;
      loDev = (rawData[index+5] >> 5) & 0xF ;
      triggX = (loDev >> 4 & 0x1 ) && !(loDev & 0xF);

      if( locDec != 0x9 ){ // check for Dec
	
	iRegLoc = iReg*16 + iLoc;
	locId = fRegToLocCard[iRegLoc].fLocId; 
	
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

	  if (pattern[0][0] || pattern[0][1] || pattern[0][2] || pattern[0][3]
	      || pattern[1][0] || pattern[1][1] || pattern[1][2] || pattern[1][3]
	     )
	  {
	    if (nofTrigRec == fMaxRecPointsCount) {
	      HLTError("Output buffer is overflowed maximum assiged arraysize : %d, present array index : %d",
		       fMaxRecPointsCount, nofTrigRec
	      );
	      return false;
	    }

	    HLTDebug("iReg: %d, iLoc :%d, locId : %d,X : %x, %x, %x, %x ...Y : %x, %x, %x, %x",
		    iReg,iLoc,locId,pattern[0][0],pattern[0][1],pattern[0][2],pattern[0][3],
		    pattern[1][0],pattern[1][1],pattern[1][2],pattern[1][3]
	    );

	    // hitset indicates which hits on chambers 7 to 10 have been found and filled.
	    bool hitset[4] = {false, false, false, false};
	    bool Xset[4] = {false, false, false, false};
	    bool Yset[4] = {false, false, false, false};

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

		      idetElemId &= 0x1FF;
		      iPlane &= 0x1;
		      locId &= 0xFF;
		      secondLocation &= 0xF;
		      
		      lutAddress &= 0x0;
		      lutAddress = (lutAddress|idetElemId)<<1;  
		      lutAddress = (lutAddress|iPlane)<<8;  
		      lutAddress = (lutAddress|locId)<<4;
		      lutAddress |= secondLocation;

		      lutAddress -= fIdOffSet;
		      
		      if(fLookUpTableData[lutAddress+1].fIdManuChannel == -1) //skip uninitialized values
			continue;
			
		      if (iPlane == 1)
		      {      
		        trigRecord[nofTrigRec].fHit[iChamber].fX = fLookUpTableData[lutAddress+1].fRealX;
		        Xset[iChamber] = true;
		      }
		      else
		      {
		        trigRecord[nofTrigRec].fHit[iChamber].fY = fLookUpTableData[lutAddress+1].fRealY;
		        trigRecord[nofTrigRec].fHit[iChamber].fZ = fLookUpTableData[lutAddress+1].fRealZ;
		        Yset[iChamber] = true;
		      }	
			
		      HLTDebug("\t Hit Found fo ich : %d, iPlane : %d, detelem %d, id : %d, at (%lf, %lf, %lf) cm",
			       iChamber,fLookUpTableData[lutAddress+1].fPlane,detElemId,fLookUpTableData[lutAddress+1].fIdManuChannel,
			       fLookUpTableData[lutAddress+1].fRealX,
			       fLookUpTableData[lutAddress+1].fRealY,
			       fLookUpTableData[lutAddress+1].fRealZ
		      );
		      
	   	      dataCount++;
		      
		    }//pattern maching is found 
		  }// loop of ibitxy
		}// if pattern
	      }// iplane
	    }// ichamber
	    
	      {
	      
		// Fill the hitset flags and make sure the hit structures that were not
		// filled (set) get set to a nil value.
		for (int i = 0; i < 4; i++)
		{
			if (Xset[i] and Yset[i]) hitset[i] = true;
			
			if (not hitset[i])
			{
				trigRecord[nofTrigRec].fHit[i]
					= AliHLTMUONConstants::NilRecHitStruct();
			}
		}

		trigRecord[nofTrigRec].fId = fTrigRecId;
	    
		// Increment trigger record Id and keep it positive.
		//TODO: handle the wrapparound better.
		if (fTrigRecId < 0x7FFFFFFF)
			fTrigRecId++;
		else
			fTrigRecId = 0;
		
		AliHLTMUONRecHitStruct* hit1 = NULL;
		if (hitset[0])
			hit1 = &trigRecord[nofTrigRec].fHit[0];
		else if (hitset[1])
			hit1 = &trigRecord[nofTrigRec].fHit[1];
		AliHLTMUONRecHitStruct* hit2 = NULL;
		if (hitset[2])
			hit2 = &trigRecord[nofTrigRec].fHit[2];
		else if (hitset[3])
			hit2 = &trigRecord[nofTrigRec].fHit[3];
		
		if (hit1 != NULL and hit2 != NULL)
		{
			// Calculate the momentum and fill in the flags and momentum fields.
			AliHLTMUONCalculations::ComputeMomentum(
					hit1->fX,
					hit1->fY, hit2->fY,
					hit1->fZ, hit2->fZ
				);
			trigRecord[nofTrigRec].fPx = AliHLTMUONCalculations::Px();
			trigRecord[nofTrigRec].fPy = AliHLTMUONCalculations::Py();
			trigRecord[nofTrigRec].fPz = AliHLTMUONCalculations::Pz();

			trigRecord[nofTrigRec].fFlags =
				AliHLTMUONUtils::PackTriggerRecordFlags(
					AliHLTMUONCalculations::Sign(),
					hitset
				);
			
			nofTrigRec++;
		}
		else if ((hit1 != NULL or hit2 != NULL) and not suppressPartialTrigs)
		{
			trigRecord[nofTrigRec].fPx = 0;
			trigRecord[nofTrigRec].fPy = 0;
			trigRecord[nofTrigRec].fPz = 0;

			trigRecord[nofTrigRec].fFlags =
				AliHLTMUONUtils::PackTriggerRecordFlags(
					kSignUnknown,
					hitset
				);
			
			nofTrigRec++;
		}
	    
		//int xPos = rawData[index] & 0x1F;
		//int dev = (rawData[index]>>5) & 0x1F;
		//int yPos = (rawData[index]>>10) & 0xF;
	    }
	  
	  }// if any non zero pattern found

	  index += 1 ; // the last word, important one
	}// if locId <=234
      }// Dec Condn

      if(!reg_phys_trig_occur)
	 index += 45;
	
      index += 1; // end of local Data 0xCAFEFADE

      HLTDebug("iReg %d, iLoc %d, locId : %d, trigY %x, triggX %x, loDev %x, dec %x, sign %x,rawData : %x",
	       iReg,iLoc,locId,triggY,triggX,loDev,dec,sign, rawData[index]);

      index = iLocIndex + 6; //important to reset the index counter for fake locids like 235 
      
     }// iLoc loop
     
  }// iReg Loop

  return true;
}
