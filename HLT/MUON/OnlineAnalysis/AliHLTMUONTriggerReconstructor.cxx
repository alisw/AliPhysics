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
              give the output AliHLTMUONTriggerRecordStruct
 Author     : Indranil Das, HEP Division, SINP
 Email      : indra.das@saha.ac.in | indra.ehep@gmail.com
**********************************************************************/

///*
//
//  The TrigRec class is designed to deal the rawdata inputfiles to findout the 
//  the reconstructed hits at the trigger DDL. The output is send to the output block for further 
//  processing.
//
//  Author : Indranil Das ( indra.das@saha.ac.in || indra.ehep@gmail.com )
// 
//*/

#if __GNUC__ >= 3
using namespace std;
#endif

#include <vector>

#include "TObjArray.h"

#include "AliHLTMUONTriggerReconstructor.h"

#include "AliMUONTriggerCrate.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONTriggerCircuit.h"

#include "AliMpPad.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpDDLStore.h"


const int AliHLTMUONTriggerReconstructor::fgkDetectorId = 0xB00;
const int AliHLTMUONTriggerReconstructor::fgkDDLOffSet = 20 ;
const int AliHLTMUONTriggerReconstructor::fgkNofDDL = 2 ;

const int AliHLTMUONTriggerReconstructor::fgkDDLHeaderSize = 8;

const int AliHLTMUONTriggerReconstructor::fgkEvenLutSize =  5208448+ 1;
const int AliHLTMUONTriggerReconstructor::fgkOddLutSize = 5058432 + 1;

const int AliHLTMUONTriggerReconstructor::fgkLutLine = 10496;

const int AliHLTMUONTriggerReconstructor::fgkMinIdManuChannel[2] = {1638400, 1720320};
const int AliHLTMUONTriggerReconstructor::fgkMaxIdManuChannel[2] = {6846848, 6778752};

const float AliHLTMUONTriggerReconstructor::fgkHalfPadSizeXB[3] = {8.5, 17.0, 25.5};
const float AliHLTMUONTriggerReconstructor::fgkHalfPadSizeYNB[2] = {25.5, 34.0};

const int AliHLTMUONTriggerReconstructor::fgkDetElem = 9*4 ; // 9 detele per half chamber

AliHLTMUONTriggerReconstructor::AliHLTMUONTriggerReconstructor() :
  fPadData(NULL),
  fLookUpTableData(NULL),
  fRecPoints(NULL),
  fRecPointsCount(NULL),
  fMaxRecPointsCount(0),
  fDigitPerDDL(0),
  fNofFiredDetElem(0),
  fMaxFiredPerDetElem(NULL),
  fDetManuChannelIdList(NULL),
  fCentralChargeB(NULL),
  fCentralChargeNB(NULL),
  fDDLId(0),
  fIdOffSet(0),
  fCrateManager(NULL)
{
  // ctor 
  
  if(AliHLTMUONTriggerReconstructor::fgkEvenLutSize > AliHLTMUONTriggerReconstructor::fgkOddLutSize){
    fPadData = new AliHLTMUONHitReconstructor::DHLTPad[AliHLTMUONTriggerReconstructor::fgkEvenLutSize];
  }
  else{
    fPadData = new AliHLTMUONHitReconstructor::DHLTPad[AliHLTMUONTriggerReconstructor::fgkOddLutSize];
  }

  fMaxFiredPerDetElem = new int[fgkDetElem] ;

  fCrateManager = new AliMUONTriggerCrateStore();   
  fCrateManager->ReadFromFile();

  bzero(fGetIdTotalData,104*64*2*sizeof(int));
}


AliHLTMUONTriggerReconstructor::~AliHLTMUONTriggerReconstructor()
{
  // dtor

  //HLTError("\nEnd of Run\n");

  delete []fPadData;
  delete []fLookUpTableData;
  delete fCrateManager ;
  delete []fMaxFiredPerDetElem;
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

  fLookUpTableData[0].fIdManuChannel = 0;
  fLookUpTableData[0].fIX = 0 ;
  fLookUpTableData[0].fIY = 0 ;
  fLookUpTableData[0].fRealX = 0.0 ;
  fLookUpTableData[0].fRealY = 0.0 ;
  fLookUpTableData[0].fRealZ = 0.0 ;
  fLookUpTableData[0].fPlane = -1 ;
  fLookUpTableData[0].fPcbZone = -1 ;

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

bool AliHLTMUONTriggerReconstructor::Run(int *rawData, int *rawDataSize, AliHLTMUONTriggerRecordStruct trigRecord[], int *nofTrigRec)
{

  fRecPoints = &trigRecord[0];
  fMaxRecPointsCount = *nofTrigRec;
  fRecPointsCount = nofTrigRec;
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
    HLTError("Failed to read the complete DDL file\n");
    return false;
  }

  if(!FindTrigHits()){
    HLTError("Failed to generate RecHits\n");
    return false;
  }
    
  return true;
}


bool AliHLTMUONTriggerReconstructor::ReadDDL(int *rawData, int *rawDataSize)
{

  vector<AliHLTMUONHitReconstructor::DHLTPad> padList;
  int idManuChannel ;
  
  int index = 0;
  int dataCount = 0;
  fNofFiredDetElem = 0;
  int detElemId = 0 ;
  int prevDetElemId = 0 ;

  fDetManuChannelIdList = new int[(*rawDataSize)];

#ifdef DEBUG
  int globalcard_data_occurance = (rawData[index]>>10)&0x1; //Set to 1 if global info present in DDL else set to 0 
  int version = (rawData[index]>>12)&0xFF; // software version
  int serial_number =  (rawData[index]>>20)&0xF; // serial number set to 0xF 
#endif
  int phys_trig_occur = (rawData[index]>>30)&0x1; // 1 for physics trigger, 0 for software trigger
  
  // Values not set
//   int regional_structure = (rawData[index])&0xFF ; 
//   int DAQ_interfaced = (rawData[index]>>8)&0x1;
//   int central_or_LTU = (rawData[index]>>9)&0x1;
//   int VME_trigger = (rawData[index]>>11)&0x1;
//   int DARC_type =  (rawData[index]>>24)&0x3;
//   int dimuon_ZDC = (rawData[index]>>27)&0x3;
//   int MBZ = (rawData[index]>>31)&0x1;

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
#endif // DEBUG

    HLTDebug("singleLpt : %x, singleHpt : %x, pairUnlikeLpt : %x, pairUnlikeHpt : %x, pairLikeLpt : %x, pairLikeHpt : %x",
	     singleLpt,singleHpt,pairUnlikeLpt,pairUnlikeHpt,pairLikeLpt,pairLikeHpt);
  }

  if(!phys_trig_occur)
    index += 10 ;// corresponds to scalar words
  
  index += 1; // separator 0xDEADBEEF 

  for (int iReg = 0; iReg < 8; iReg++) {
    index += 1; // DARC Status Word
    index += 1; // Regeional Word
    //int reg_output = rawData[index] & 0xFF;
    int reg_phys_trig_occur = ( rawData[index] >> 31) & 0x1;
//     int reg_version;// = ;
//     int reg_Id ;//= ;
//     int reg_serial_number;// = ;
    
    index += 2; // 2 words for regional input
    
    index += 1; // L0 counter
    
    if(!reg_phys_trig_occur)
      index += 10;
    
    index += 1 ; // end of Regeonal header
    
    AliMUONTriggerCrate* crate = fCrateManager->Crate((fDDLId - AliHLTMUONTriggerReconstructor::fgkDDLOffSet), iReg);
    TObjArray *boards = crate->Boards();


    for(int iLoc = 0; iLoc < 16 ; iLoc++){

      int iLocIndex = index ;      

      int locId = (rawData[index+5] >> 19) &  0xF ;
      
      AliMUONLocalTriggerBoard* localBoard = (AliMUONLocalTriggerBoard*)boards->At(locId + 1);
      int iLocCard = localBoard->GetNumber();

      int dec = (rawData[index+5] >> 15) & 0xF;
      //int triggY = (rawData[index+5] >> 14) & 0x1 ;
      //int sign = (rawData[index+5] >> 9) & 0x1;
      //int loDev = (rawData[index+5] >> 5) & 0xF ;
      //int triggX = (loDev >> 4 & 0x1 ) && !(loDev & 0xF) ;


//       HLTDebug(" \n",
// 	       iLocCard,locId);

//      index += 5;

      if(iLocCard > 0){
	if( dec != 0x9 ){ // check for Dec

	  index += 1;
	  short X1_pattern = rawData[index] & 0xFFFF; 
	  short X2_pattern = (rawData[index] >> 16) & 0xFFFF; 
	  index += 1; 
	  short X3_pattern = rawData[index] & 0xFFFF; 
	  short X4_pattern = (rawData[index] >> 16) & 0xFFFF; 
	  
	  index += 1;
	  short Y1_pattern = rawData[index] & 0xFFFF; 
	  short Y2_pattern = (rawData[index] >> 16) & 0xFFFF; 
	  index += 1; 
	  short Y3_pattern = rawData[index] & 0xFFFF; 
	  short Y4_pattern = (rawData[index] >> 16) & 0xFFFF; 
	  
	  TArrayS xyPattern[2];
	  xyPattern[0].Set(4);
	  xyPattern[1].Set(4);
	  
	  xyPattern[0].AddAt(X1_pattern,0);
	  xyPattern[0].AddAt(X2_pattern,1);
	  xyPattern[0].AddAt(X3_pattern,2);
	  xyPattern[0].AddAt(X4_pattern,3);
	  
	  xyPattern[1].AddAt(Y1_pattern,0);
	  xyPattern[1].AddAt(Y2_pattern,1);
	  xyPattern[1].AddAt(Y3_pattern,2);
	  xyPattern[1].AddAt(Y4_pattern,3);
	  
	  HLTDebug("iLocCard : %d, locId : %d, X : %x, %x, %x, %x .... Y : %x, %x, %x, %x\n",
		   iLocCard,locId,X1_pattern,X2_pattern,X3_pattern,X4_pattern,
		   Y1_pattern,Y2_pattern,Y3_pattern,Y4_pattern);
	
	  index += 1 ; // skipping the last word though it is important
	  
	  padList.clear() ;
	  if( Pattern2Pad(iLocCard, xyPattern, padList) ) {

	    for (UInt_t iEntry = 0; iEntry < padList.size(); iEntry++) {
	      
	      AliHLTMUONHitReconstructor::DHLTPad dPad = padList[iEntry];
	      
	      detElemId = dPad.fDetElemId;
	      int idetElemId = (detElemId)%1000;
	      idetElemId &= 0x1FF ;
	      int iPlane = dPad.fPlane & 0x1 ;
	      int iX = dPad.fIX & 0x7F ;
	      int iY = dPad.fIY & 0x3F ;
	      
	      idManuChannel &= 0x0;
	      idManuChannel = (idManuChannel|idetElemId)<<1;  
	      idManuChannel = (idManuChannel|iPlane)<<7;  
	      idManuChannel = (idManuChannel|iX)<<6 ;
	      idManuChannel |= iY ;
	      idManuChannel -= fIdOffSet ;
	      
	      fPadData[idManuChannel].fDetElemId = dPad.fDetElemId;
	      fPadData[idManuChannel].fIdManuChannel = idManuChannel;
	      fPadData[idManuChannel].fIX = fLookUpTableData[idManuChannel+1].fIX;
	      fPadData[idManuChannel].fIY = fLookUpTableData[idManuChannel+1].fIY;
	      fPadData[idManuChannel].fRealX = fLookUpTableData[idManuChannel+1].fRealX;
	      fPadData[idManuChannel].fRealY = fLookUpTableData[idManuChannel+1].fRealY;
	      fPadData[idManuChannel].fRealZ = fLookUpTableData[idManuChannel+1].fRealZ;
	      fPadData[idManuChannel].fPcbZone = fLookUpTableData[idManuChannel+1].fPcbZone;
	      fPadData[idManuChannel].fPlane = fLookUpTableData[idManuChannel+1].fPlane;

 	      fDetManuChannelIdList[dataCount] = idManuChannel;
	      if(detElemId != prevDetElemId){
		if(fNofFiredDetElem>0){
		  fMaxFiredPerDetElem[fNofFiredDetElem-1] = dataCount;
		}
		fNofFiredDetElem++;
		prevDetElemId = detElemId ;
	      } // if detelem condn
	      dataCount ++;

// 	      printf("detelemId : %d, plane : %d, IX : %d, IY : %d realX : %f, realY : %f , realZ %f\n",
// 		     dPad.fDetElemId,dPad.fPlane,dPad.fIX,dPad.fIY,
// 		     fPadData[idManuChannel].fRealX,fPadData[idManuChannel].fRealY,fPadData[idManuChannel].fRealZ);
	      
	    }// for loop of entry

	  }//pattern2pad 
	  
      }// Dec Condn

      }// iLocCard > 0

      if(!reg_phys_trig_occur)
	index += 45;
	
      index += 1; // end of local Data
      HLTDebug("iReg %d, iLoc %d, iLocCard : %d, locId : %d, trigY %x, triggX %x, loDev %x, dec %x, sign %x,rawData : %x",
	       iReg,iLoc,iLocCard,locId,triggY,triggX,loDev,dec,sign, rawData[index]);

      index = iLocIndex + 6 ;
      
//       delete localBoard;

    }// iLoc loop
  
//     delete crate;
//     delete boards;
    }// iReg Loop

//   fDigitPerDDL = dataCount;
//   fMaxFiredPerDetElem[fNofFiredDetElem-1] = dataCount;
  

  return true;
}

bool AliHLTMUONTriggerReconstructor::Pattern2Pad(int nBoard, TArrayS* xyPattern, vector<AliHLTMUONHitReconstructor::DHLTPad>& padList)
{

  Int_t detElemId;
  Int_t previousDetElemId[4] = {0};
  Int_t previousBoard[4] = {0};


  // loop over x1-4 and y1-4
  for(Int_t iChamber = 0; iChamber < 4; ++iChamber){
    for(Int_t iCath = 0; iCath < 2; ++iCath){
      //int index = 0;  
      Int_t pattern = (Int_t)xyPattern[iCath].At(iChamber); 
      if (!pattern) continue;
      
      // get detElemId
      /*
      AliMUONTriggerCircuit triggerCircuit;
      AliMUONLocalTriggerBoard* localBoard = fCrateManager->LocalBoard(nBoard);
      detElemId = triggerCircuit.DetElemId(iChamber+10, localBoard->GetName());//FIXME +/-10 (should be ok with new mapping)
      */
      AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
      if (ddlStore != NULL) continue;
      AliMpLocalBoard* localBoard = ddlStore->GetLocalBoard(nBoard);
      if (localBoard == NULL) continue;
      detElemId = localBoard->GetDEIdByChamber(iChamber);


      if(iCath == 1){ // FIXME should find a more elegant way
	// Don't save twice the same digit
	// (since strips in non bending plane can cross several boards)
	Int_t prevDetElemId = previousDetElemId[iChamber];
	Int_t prevBoard = previousBoard[iChamber];
	previousDetElemId[iChamber] = detElemId;
	previousBoard[iChamber] = nBoard;

	if(detElemId == prevDetElemId){
	  if(nBoard-prevBoard==1) continue;
	}
      }

      const AliMpVSegmentation* seg 
	= AliMpSegmentation::Instance()
	  ->GetMpSegmentation(detElemId, AliMp::GetCathodType(iCath));  

      // loop over the 16 bits of pattern
      for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy) {
	
	if ((pattern >> ibitxy) & 0x1) {
	  
	  //Int_t temp = (pattern >> ibitxy) & 0x1 ;
	  // not quite sure about this
	  Int_t offset = 0;
	  if (iCath && localBoard->GetSwitch(6)) offset = -8;

	  AliMpPad pad = seg->PadByLocation(AliMpIntPair(nBoard,ibitxy+offset),kTRUE);

	  AliHLTMUONHitReconstructor::DHLTPad dPad;
	  if (!pad.IsValid()) {
	    //AliWarning(Form("No pad for detElemId: %d, nboard %d, ibitxy: %d\n",
	    //	    detElemId, nBoard, ibitxy));
	    continue ;
	  } // 

	  Int_t padX = pad.GetIndices().GetFirst();
	  Int_t padY = pad.GetIndices().GetSecond();
	  // file digit
	  dPad.fIX = padX ;
	  dPad.fIY = padY ;
	  dPad.fPlane = iCath ;
	  dPad.fDetElemId = detElemId ;
          //printf("nBoard : %d, detElemId : %d, chamber %d, iCath %d, ibitxy %d, pattern %d, switch %d, offset %d, temp %x, padX : %d, padY: %d\n",nBoard,detElemId,iChamber,iCath,ibitxy,pattern,localBoard->GetSwitch(6), offset, temp,padX,padY);
	  //dPad.fDetElemId = detElemId ;

	  padList.push_back(dPad);
	
	}// xyPattern
      }// ibitxy
//       delete localBoard;
    }// cath
  } // ichamber

  return true;
}

bool AliHLTMUONTriggerReconstructor::FindTrigHits() 
{
  for(int iDet=0; iDet<fNofFiredDetElem ; iDet++){
    
    if(iDet>0)
      MergeTrigHits(fMaxFiredPerDetElem[iDet-1],fMaxFiredPerDetElem[iDet]);
    else
      MergeTrigHits(0,fMaxFiredPerDetElem[iDet]);

    
//     if(iDet==0)
//       for(int i=0;i<fMaxFiredPerDetElem[iDet];i++)
// 	fGetIdTotalData[fPadData[fDetManuChannelIdList[i]].fIX][fPadData[fDetManuChannelIdList[i]].fIY][fPadData[fDetManuChannelIdList[i]].fPlane] = 0;
//     else
//       for(int i=fMaxFiredPerDetElem[iDet-1];i<fMaxFiredPerDetElem[iDet];i++)
// 	fGetIdTotalData[fPadData[fDetManuChannelIdList[i]].fIX][fPadData[fDetManuChannelIdList[i]].fIY][fPadData[fDetManuChannelIdList[i]].fPlane] = 0;

  }// loop over detection element
    
  //for(int iPad=fDataPerDetElem[i];iPad<fDataPerDetElem[i+1];iPad++){
  for(int iPad=0;iPad<fDigitPerDDL;iPad++){
//     fGetIdTotalData[fPadData[fDetManuChannelIdList[iPad]].fIX][fPadData[fDetManuChannelIdList[iPad]].fIY][fPadData[fDetManuChannelIdList[iPad]].fPlane] = 0;
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
  
  for(int i=0;i<fgkDetElem;i++)
    fMaxFiredPerDetElem[i] = 0;

  delete []fDetManuChannelIdList;

  return true;
}

bool AliHLTMUONTriggerReconstructor::MergeTrigHits(int minPadId, int maxPadId)
{
  int idManuChannelB, idManuChannelNB;
  float halfPadLengthX,halfPadLengthY;
  float diffX,diffY;

  for(int iPad=minPadId;iPad<maxPadId;iPad++){
    idManuChannelB   = fDetManuChannelIdList[iPad];
    //printf("idManuChannelB : %d, fPadData[idManuChannelB].fPlane : %d\n",idManuChannelB,fPadData[idManuChannelB].fPlane);
    if(fPadData[idManuChannelB].fPlane == 0){
      
      halfPadLengthX = AliHLTMUONTriggerReconstructor::fgkHalfPadSizeXB[fPadData[idManuChannelB].fPcbZone] ;

      for(int iPad=minPadId;iPad<maxPadId;iPad++){
	idManuChannelNB   = fDetManuChannelIdList[iPad];
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
	  //printf("diffX %f,  halfPadLengthX %f,  diffY  %f, halfPadLengthY  %f\n",diffX,halfPadLengthX,diffY,halfPadLengthY);

	  if(diffX < halfPadLengthX + 1.0 && diffY < halfPadLengthY + 1.0 ){// added redundancy of 1.0 cm due to the pb of geometrical segmentation 

	    AliHLTMUONRecHitStruct hit;
	    hit.fX = fPadData[idManuChannelNB].fRealX;
	    hit.fY = fPadData[idManuChannelNB].fRealY;
	    hit.fZ = fPadData[idManuChannelNB].fRealZ;

	    int ichamber = int((fPadData[idManuChannelB].fDetElemId - 1000)/100);
	    ichamber--;
	    fRecPoints[(*fRecPointsCount)].fHit[ichamber] = hit;
	    
	    (*fRecPointsCount)++;
	    if((*fRecPointsCount) == fMaxRecPointsCount){
	      printf("Nof RecHit (i.e. %d) exceeds the max nof RecHit limit %d\n",(*fRecPointsCount),fMaxRecPointsCount);
	      return false;
	    }

	    printf("ichamber : %d, detelem : %d, x %f, y %f, z %f\n",ichamber,fPadData[idManuChannelB].fDetElemId,fPadData[idManuChannelNB].fRealX,
		   fPadData[idManuChannelB].fRealY,fPadData[idManuChannelB].fRealZ);
	  }
	  
	  
	  

	}//condn for non-bending plane
      }//for loop for non-bending plane

    }// condn for bending plane
  }// for loop for bending plane

  return true;
}

