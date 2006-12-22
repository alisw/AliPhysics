#include <strings.h>
#include "HLTMUONHitReconstructor.h"


const int HLTMUONHitReconstructor::fgkDetectorId = 0xA00;
const int HLTMUONHitReconstructor::fgkDDLOffSet = 12 ;
const int HLTMUONHitReconstructor::fgkNofDDL = 8 ;

const int HLTMUONHitReconstructor::fgkDDLHeaderSize = 8;

const int HLTMUONHitReconstructor::fgkEvenLutSize = 3364287 + 1;
const int HLTMUONHitReconstructor::fgkOddLutSize = 1645631 + 1;

const int HLTMUONHitReconstructor::fgkLutLine[2] = {54208, 59648};

const int HLTMUONHitReconstructor::fgkMinIdManuChannel[2] = {64, 917696};
const int HLTMUONHitReconstructor::fgkMaxIdManuChannel[2] = {3364351,2563007};

const float HLTMUONHitReconstructor::fgkHalfPadSize[3] = {1.25, 2.50, 5.00};

const int HLTMUONHitReconstructor::GetLutLine(int iDDL){ return ( iDDL<16 ) ? fgkLutLine[0] : fgkLutLine[1] ;}

HLTMUONHitReconstructor::HLTMUONHitReconstructor(): 
  fDCCut(0),
  fDebugLevel(0)
{
  // ctor 
  
  if(HLTMUONHitReconstructor::fgkEvenLutSize > HLTMUONHitReconstructor::fgkOddLutSize){
    fPadData = new DHLTPad[HLTMUONHitReconstructor::fgkEvenLutSize];
  }
  else{
    fPadData = new DHLTPad[HLTMUONHitReconstructor::fgkOddLutSize];
  }


  fkBlockHeaderSize    = 8;
  fkDspHeaderSize      = 8;
  fkBuspatchHeaderSize = 4;

  bzero(fGetIdTotalData,336*80*2*sizeof(int));
}


HLTMUONHitReconstructor::HLTMUONHitReconstructor(const HLTMUONHitReconstructor& /*rhs*/)
{
// Protected copy constructor

  printf("Not implemented.\n");
}


HLTMUONHitReconstructor & 
HLTMUONHitReconstructor::operator=(const HLTMUONHitReconstructor& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  printf("Not implemented.\n");
    
  return *this;  
}

HLTMUONHitReconstructor::~HLTMUONHitReconstructor()
{
  // dtor

  printf("\nEnd of Run\n");

  if(fPadData != NULL)
    delete []fPadData;
  
  if(fLookUpTableData!=NULL)
    delete [] fLookUpTableData;


}

void HLTMUONHitReconstructor::CleanUp()
{

  if(fDetManuChannelIdList != NULL)
    delete [] fDetManuChannelIdList;

  if(fCentralChargeB != NULL)
    delete [] fCentralChargeB;

  if(fCentralChargeNB != NULL)
    delete [] fCentralChargeNB;
  
  if(fRecX != NULL)
    delete [] fRecX;

  if(fRecY != NULL)
    delete [] fRecY;
    
  if(fAvgChargeX != NULL)
    delete [] fAvgChargeX;

  if(fAvgChargeY != NULL)
    delete [] fAvgChargeY;
  
}

bool HLTMUONHitReconstructor::LoadLookUpTable(DHLTLut* lookUpTableData, int lookUpTableId)
{

  if(fLookUpTableData!=NULL)
     delete [] fLookUpTableData;
  

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

bool HLTMUONHitReconstructor::SetBusToDetMap(BusToDetElem busToDetElem)
{

  // function that loads BusPatch To Detection Element (SlatId) map

  if(busToDetElem.size()==0)
    return false;
  else
    fBusToDetElem = busToDetElem;
  
  return true;
}


bool HLTMUONHitReconstructor::Run(const int* rawData, int rawDataSize, DHLTRecPoint recHit[], int *nofHit) 
{  
  // main function called by HLTReconstructor to perform DHLT Hitreconstruction 
#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::Run", "Trace")
		<< "Started Run(0x" << AliHLTLog::kHex << (unsigned long)rawData
		<< ", " << AliHLTLog::kDec << rawDataSize
		<< ", 0x" << AliHLTLog::kHex << (unsigned long)(recHit)
		<< ", " << AliHLTLog::kDec << *nofHit
		<< ") ..." << ENDLOG;
#	endif

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
#	ifdef DEBUG
	LOG(AliHLTLog::kError, "HLTMUONHitReconstructor::Run", "Trace")
		<< "Failed to read the DDL file." << ENDLOG;
	LOG(AliHLTLog::kError, "HLTMUONHitReconstructor::Run", "Trace")
		<< "Please check whether end of event has been reached or probably the input rawdata path is missplled." << ENDLOG;
#	endif
    printf("Failed to read the DDL file\n");
    printf("Please check whether end of event has been reached or probably the input rawdata path is missplled\n");
    CleanUp();
    return false;
  }

  if(!FindRecHits()){
#	ifdef DEBUG
	LOG(AliHLTLog::kError, "HLTMUONHitReconstructor::Run", "Trace")
		<< "Failed to generate RecHits." << ENDLOG;
#	endif
    printf("Failed to generate RecHits\n");
    CleanUp();
    return false;
  }
    
  return true;
}



bool HLTMUONHitReconstructor::ReadDDL(const int* rawData, int rawDataSize)
{
  //function to read Raw Data files
#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::ReadDDL", "Trace")
		<< "Start of ReadDDL(0x" << AliHLTLog::kHex << (unsigned long)rawData << ", "
		<< AliHLTLog::kDec << rawDataSize << ") ..." << ENDLOG;
#	endif

  int ddlRawDataSize;
  ddlRawDataSize = rawDataSize / 4;

//  int *buffer = new int[ddlRawDataSize]; 
  const int* buffer = rawData;

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
#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::ReadDDL", "Trace")
		<< "Reading iBlock = " << AliHLTLog::kDec << iBlock
		<< ", totalBlockSize = " << AliHLTLog::kDec << totalBlockSize
		<< ", blockRawDataSize = " << AliHLTLog::kDec << blockRawDataSize
		<< ", indexDsp = " << AliHLTLog::kDec << indexDsp
		<< ENDLOG;
	cerr << "HLTMUONHitReconstructor::ReadDDL - Trace -"
		<< "Reading iBlock = " << iBlock
		<< ", totalBlockSize = " << totalBlockSize
		<< ", blockRawDataSize = " << blockRawDataSize
		<< ", indexDsp = " << indexDsp
		<< endl;
#	endif
    while(blockRawDataSize > 0){
      totalDspSize = buffer[indexDsp + 1];
      dspRawDataSize = buffer[indexDsp + 2];
      dspRawDataSize --;                              // temporary solution to read buspatches 
      indexBuspatch = indexDsp + fkDspHeaderSize + 2; // this extra 2 word comes from the faulty defination of Dsp header size
#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::ReadDDL", "Trace")
		<< "DSP block: totalDspSize = " << AliHLTLog::kDec << totalDspSize
		<< ", dspRawDataSize = " << AliHLTLog::kDec << dspRawDataSize
		<< ", indexBuspatch = " << AliHLTLog::kDec << indexBuspatch
		<< ENDLOG;
	cerr << "HLTMUONHitReconstructor::ReadDDL - Trace - "
		<< "DSP block: totalDspSize = " << totalDspSize
		<< ", dspRawDataSize = " << dspRawDataSize
		<< ", indexBuspatch = " << indexBuspatch
		<< endl;
#	endif
      while(dspRawDataSize > 0){
	totalBuspatchSize = buffer[indexBuspatch + 1];
	buspatchRawDataSize = buffer[indexBuspatch + 2];
	buspatchId = buffer[indexBuspatch + 3];
	detElemId = fBusToDetElem[buspatchId];
	indexRawData = indexBuspatch + fkBuspatchHeaderSize;
#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::ReadDDL", "Trace")
		<< "Buspatch block: totalBuspatchSize = " << AliHLTLog::kDec << totalBuspatchSize
		<< ", buspatchRawDataSize = " << AliHLTLog::kDec << buspatchRawDataSize
		<< ", buspatchId = " << AliHLTLog::kDec << buspatchId
		<< ", detElemId = " << AliHLTLog::kDec << detElemId
		<< ", indexRawData = " << AliHLTLog::kDec << indexRawData
		<< ENDLOG;
	cerr << "HLTMUONHitReconstructor::ReadDDL - Trace - "
		<< "Buspatch block: totalBuspatchSize = " << totalBuspatchSize
		<< ", buspatchRawDataSize = " << buspatchRawDataSize
		<< ", buspatchId = " << buspatchId
		<< ", detElemId = " << detElemId
		<< ", indexRawData = " << indexRawData
		<< endl;
#	endif
	while(buspatchRawDataSize > 0){
	  dataWord = buffer[indexRawData];
	  charge = (unsigned short)(dataWord & 0xFFF);
	  
	  idManuChannel = 0x0;
	  idManuChannel = (idManuChannel|(detElemId%100))<<17;
	  idManuChannel |= (dataWord >> 12) & 0x1FFFF;
	  idManuChannel -= fIdOffSet ;
	  
#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::ReadDDL", "Trace")
		<< "Found pad [buspatchId = " << AliHLTLog::kDec << buspatchId
		<< ", detElemId = " << AliHLTLog::kDec << detElemId
		<< ", idManuChannel = " << AliHLTLog::kDec << idManuChannel
		<< ", charge = " << AliHLTLog::kDec << charge
		<< "]" << ENDLOG;
	cerr << "HLTMUONHitReconstructor::ReadDDL - Trace - "
		<< "Found pad [buspatchId = " << buspatchId
		<< ", detElemId = " << detElemId
		<< ", idManuChannel = " << idManuChannel
		<< ", charge = " << charge
		<< "]" << endl;
#	endif

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

#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::ReadDDL", "Trace")
		<< "Keeping pad [fBuspatchId = " << AliHLTLog::kDec << fPadData[idManuChannel].fBuspatchId
		<< ", fDetElemId = " << AliHLTLog::kDec << fPadData[idManuChannel].fDetElemId
		<< ", fIdManuChannel = " << AliHLTLog::kDec << fPadData[idManuChannel].fIdManuChannel
		<< ", fIX = " << AliHLTLog::kDec << fPadData[idManuChannel].fIX
		<< ", fIY = " << AliHLTLog::kDec << fPadData[idManuChannel].fIY
		<< ", fRealX = " << AliHLTLog::kDec << fPadData[idManuChannel].fRealX
		<< ", fRealY = " << AliHLTLog::kDec << fPadData[idManuChannel].fRealY
		<< ", fRealZ = " << AliHLTLog::kDec << fPadData[idManuChannel].fRealZ
		<< ", fPcbZone = " << AliHLTLog::kDec << fPadData[idManuChannel].fPcbZone
		<< ", fPlane = " << AliHLTLog::kDec << fPadData[idManuChannel].fPlane
		<< ", fCharge = " << AliHLTLog::kDec << fPadData[idManuChannel].fCharge
		<< "]" << ENDLOG;
	cerr << "HLTMUONHitReconstructor::ReadDDL - Trace - "
		<< "Keeping pad [fBuspatchId = " << fPadData[idManuChannel].fBuspatchId
		<< ", fDetElemId = " << fPadData[idManuChannel].fDetElemId
		<< ", fIdManuChannel = " << fPadData[idManuChannel].fIdManuChannel
		<< ", fIX = " << fPadData[idManuChannel].fIX
		<< ", fIY = " << fPadData[idManuChannel].fIY
		<< ", fRealX = " << fPadData[idManuChannel].fRealX
		<< ", fRealY = " << fPadData[idManuChannel].fRealY
		<< ", fRealZ = " << fPadData[idManuChannel].fRealZ
		<< ", fPcbZone = " << fPadData[idManuChannel].fPcbZone
		<< ", fPlane = " << fPadData[idManuChannel].fPlane
		<< ", fCharge = " << fPadData[idManuChannel].fCharge
		<< "]" << endl;
#	endif
	    
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

#	ifdef DEBUG
	cerr << "HLTMUONHitReconstructor::ReadDDL - Trace - "
		<< "End of Buspatch block: totalBuspatchSize = " << AliHLTLog::kDec << totalBuspatchSize
		<< ", buspatchRawDataSize = " << AliHLTLog::kDec << buspatchRawDataSize
		<< ", buspatchId = " << AliHLTLog::kDec << buspatchId
		<< ", detElemId = " << AliHLTLog::kDec << detElemId
		<< ", indexRawData = " << AliHLTLog::kDec << indexRawData
		<< endl;
#	endif

      }// buspatch loop
      indexDsp += totalDspSize;
      blockRawDataSize -= totalDspSize;

#	ifdef DEBUG
	cerr << "HLTMUONHitReconstructor::ReadDDL - Trace - "
		<< "End of DSP block: totalDspSize = " << AliHLTLog::kDec << totalDspSize
		<< ", dspRawDataSize = " << AliHLTLog::kDec << dspRawDataSize
		<< ", indexBuspatch = " << AliHLTLog::kDec << indexBuspatch
		<< endl;
#	endif

    }// DSP loop
    index = totalBlockSize;

#	ifdef DEBUG
	cerr << "HLTMUONHitReconstructor::ReadDDL - Trace - "
		<< "End of block: iBlock = " << iBlock << endl;
#	endif

  }// Block loop
  
//  delete[] buffer;
  
  fDigitPerDDL = dataCount;
  fMaxFiredPerDetElem[fNofFiredDetElem-1] = dataCount;
  
  
  return true;

}

bool HLTMUONHitReconstructor::FindRecHits() 
{
  // fuction that calls hit reconstruction detector element-wise   
#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::FindRecHits", "Trace")
		<< "Start of FindRecHits() ..." << ENDLOG;
#	endif

  for(int iDet=0; iDet<fNofFiredDetElem ; iDet++){

#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::FindRecHits", "Trace")
		<< "Processing iDet = " << AliHLTLog::kDec << iDet << ENDLOG;
#	endif
    
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

void HLTMUONHitReconstructor::FindCentralHits(int minPadId, int maxPadId)
{
  // to find central hit associated with each cluster
#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::FindCentralHits", "Trace")
		<< "Start of FindCentralHits(" << AliHLTLog::kDec << minPadId << ", "
		<< AliHLTLog::kDec << maxPadId << ") ..." << ENDLOG;
#	endif

  int b,nb;
  int idManuChannelCentral;
  bool hasFind;
  int idManuChannel;
  
  for(int iPad=minPadId;iPad<maxPadId;iPad++){
    idManuChannel   = fDetManuChannelIdList[iPad];
    
    
    fGetIdTotalData[fPadData[idManuChannel].fIX]
      [fPadData[idManuChannel].fIY]
      [fPadData[idManuChannel].fPlane] = idManuChannel;
    
    if(fPadData[idManuChannel].fPlane == 0){
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

void HLTMUONHitReconstructor::RecXRecY()
{
  // find reconstructed X and Y for each plane separately
#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::RecXRecY", "Trace")
		<< "Start of RecXRecY() ..." << ENDLOG;
#	endif

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

bool HLTMUONHitReconstructor::MergeRecHits()
{
  // Merge reconstructed hits first over same plane then bending plane with non-bending plane
#	ifdef DEBUG
	LOG(AliHLTLog::kDebug, "HLTMUONHitReconstructor::MergeRecHits", "Trace")
		<< "Start of MergeRecHits() ..." << ENDLOG;
#	endif

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

	    (fRecPoints + *fRecPointsCount)->X = fRecX[nb];
	    (fRecPoints + *fRecPointsCount)->Y = fRecY[b];
	    (fRecPoints + *fRecPointsCount)->Z = fPadData[idCentralB].fRealZ;
	    (fRecPoints + *fRecPointsCount)->DetElemId = fPadData[idCentralB].fDetElemId;
	    (*fRecPointsCount)++;
	    if((*fRecPointsCount) == fMaxRecPointsCount){
	      printf("Number of RecHit (i.e. %d) exceeds the max RecHit limit of %d.\n", (*fRecPointsCount), fMaxRecPointsCount);
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


