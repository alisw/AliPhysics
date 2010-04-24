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

// $Id$

///////////////////////////////////////////////
//Author : Indranil Das, SINP, INDIA
//         Sukalyan Chattopadhyay, SINP, INDIA
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

#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONClustersBlockStruct.h"
#include "AliHLTMUONChannelsBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include "AliRawDataHeader.h"
#include <cassert>


const AliHLTInt32_t AliHLTMUONHitReconstructor::fgkDetectorId = 0xA00;
const AliHLTInt32_t AliHLTMUONHitReconstructor::fgkDDLOffSet = 0;
const AliHLTInt32_t AliHLTMUONHitReconstructor::fgkNofDDL = 20;
const AliHLTInt32_t AliHLTMUONHitReconstructor::fgkDDLHeaderSize = 8;
const AliHLTInt32_t AliHLTMUONHitReconstructor::fgkLutLine = 59648 + 1;
const AliHLTInt32_t AliHLTMUONHitReconstructor::fgkMaxNofDataPerDetElem = 3000;
const AliHLTInt32_t AliHLTMUONHitReconstructor::fgkNofDetElemInDDL[20] =   {  2,   2,   2,   2, 
									      2,   2,   2,   2,
									      36,  36,  36,  36,
									      26,  26,  26,  26,
									      26,  26,  26,  26};        
const AliHLTInt32_t AliHLTMUONHitReconstructor::fgkMinDetElemIdInDDL[20] = {  100, 101, 200, 201,
									      300, 301, 400, 401,
									      505, 501, 510, 500,
									      707, 700, 807, 800,
									      907, 900,1007,1000};      


AliHLTMUONHitReconstructor::AliHLTMUONHitReconstructor() :
	AliHLTLogging(),
	fHLTMUONDecoder(),
	fkBlockHeaderSize(8),
	fkDspHeaderSize(8),
	fkBuspatchHeaderSize(4),
	fDCCut(-1),
	fPadData(NULL),
	fkLookUpTableData(NULL),
	fRecPoints(NULL),
	fRecPointsCount(NULL),
	fMaxRecPointsCount(0),
	fClusters(NULL),
	fClusterCount(0),
	fMaxClusters(0),
	fGenerateClusterInfo(false),
	fNewClusterId(0),
	fDDL(-1),
	fChannels(NULL),
	fChannelCount(0),
	fMaxChannels(0),
	fGenerateChannelInfo(false),
	fMaxChannelMult(6),
	fCentralCountB(0),
	fCentralCountNB(0),
	fDigitPerDDL(0),
	fDataCountListPerDetElem(NULL),
	fNofDataInDetElem(NULL),
	fCentralChargeB(NULL),
	fCentralChargeNB(NULL),
	fRecX(NULL),
	fRecY(NULL),
	fAvgChargeX(NULL),
	fAvgChargeY(NULL),
	fTotChargeX(NULL),
	fTotChargeY(NULL),
	fNofBChannel(NULL),
	fNofNBChannel(NULL),
	fNofYNeighbour(NULL),
	fkIdToEntry(),
	fkMaxEntryPerBusPatch(0),
	fRecoveryMode(kDontTryRecover)
{
	/// Default constructor
	
	fkBlockHeaderSize    = 8;
	fkDspHeaderSize      = 8;
	fkBuspatchHeaderSize = 4;
	
	try
	{
		fPadData = new AliHLTMUONPad[fgkLutLine];
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Dynamic memory allocation failed for fPadData in constructor.");
		throw;
	}
	
	fPadData[0].fDetElemId = 0;
	fPadData[0].fIX = 0 ;
	fPadData[0].fIY = 0 ;
	fPadData[0].fRealX = 0.0 ;
	fPadData[0].fRealY = 0.0 ;
	fPadData[0].fRealZ = 0.0 ;
	fPadData[0].fHalfPadSize = 0.0 ;
	fPadData[0].fPlane = -1 ;
	fPadData[0].fCharge = 0 ;
	fPadData[0].fBusPatch = -1;
	fPadData[0].fRawData = 0 ;
	
	bzero(fGetIdTotalData, 336*237*2*sizeof(int));
}


AliHLTMUONHitReconstructor::~AliHLTMUONHitReconstructor()
{
	/// Default destructor
	
	if (fPadData)
	{
		delete [] fPadData;
		fPadData = NULL;
	}
	
	if (fClusters != NULL)
	{
		delete [] fClusters;
	}
	if (fChannels != NULL)
	{
		delete [] fChannels;
	}
}


void AliHLTMUONHitReconstructor::SetLookUpTable(
		const AliHLTMUONHitRecoLutRow* lookupTable,
		const IdManuChannelToEntry* idToEntry,
		const MaxEntryPerBusPatch* maxEntryPerBP
	)
{
	/// Sets the Lookup table (LUT) containing the position of each pad with
	/// electronic channel associated with it. Also the appropriate manu
	/// channel ID mapping to LUT row is also set.

	assert( lookupTable != NULL );
	assert( idToEntry != NULL );
	assert( maxEntryPerBP != NULL );
	
	fkLookUpTableData = lookupTable;
	fkIdToEntry = idToEntry;
	fkMaxEntryPerBusPatch = maxEntryPerBP;
}

bool AliHLTMUONHitReconstructor::DeInitDetElemInDDLArray()
{
  /// Deinitialisation

	if (fDataCountListPerDetElem)
	{
		delete [] fDataCountListPerDetElem;
		fDataCountListPerDetElem = NULL;
	}

	if (fNofDataInDetElem)
	{
		delete [] fNofDataInDetElem;
		fNofDataInDetElem = NULL;
	}

  return true;
}

bool AliHLTMUONHitReconstructor::InitDetElemInDDLArray()
{

  ///Initialisation

  if(GetkNofDetElemInDDL(fDDL)==-1){
    HLTError("Check if the DDLNumber(AliHLTInt32_t value) is Set before this method");
    return false;
  }
  
  try{
    fDataCountListPerDetElem = new AliHLTUInt16_t*[GetkNofDetElemInDDL(fDDL)];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fDataCountListPerDetElem");
    throw;
  }

  for( Int_t idet=0;idet<GetkNofDetElemInDDL(fDDL);idet++)
    try{
      fDataCountListPerDetElem[idet] = new AliHLTUInt16_t[fgkMaxNofDataPerDetElem];
    }catch (const std::bad_alloc&){
      HLTError("Dynamic memory allocation failed in constructor : fDataCountListPerDetElem[%d]",idet);
      throw;
    }
  
  try{
    fNofDataInDetElem = new AliHLTUInt16_t[GetkNofDetElemInDDL(fDDL)];
  }catch (const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed in constructor : fDataCountListPerDetElem");
    throw;
  }

  for( Int_t idet=0;idet<GetkNofDetElemInDDL(fDDL);idet++)
    fNofDataInDetElem[idet] = 0;
  
  return true;
}

void AliHLTMUONHitReconstructor::TryRecover(ERecoveryMode mode)
{
	/// Sets if the decoder should enable the error recovery logic.
	
	// Here we setup the various flags to control exactly how the DDL raw data
	// decoder will behave and what output is generated during errors.
	fRecoveryMode = mode;
	switch (mode)
	{
	case kRecoverFull:
		fHLTMUONDecoder.TryRecover(true);
		fHLTMUONDecoder.ExitOnError(false);
		fHLTMUONDecoder.GetHandler().WarnOnly(true);
		fHLTMUONDecoder.GetHandler().PrintParityErrorAsWarning(true);
		break;
	case kRecoverJustSkip:
		fHLTMUONDecoder.TryRecover(false);
		fHLTMUONDecoder.ExitOnError(false);
		fHLTMUONDecoder.GetHandler().WarnOnly(true);
		fHLTMUONDecoder.GetHandler().PrintParityErrorAsWarning(true);
		break;
	case kRecoverFromParityErrorsOnly:
		fHLTMUONDecoder.TryRecover(false);
		fHLTMUONDecoder.ExitOnError(false);
		fHLTMUONDecoder.GetHandler().WarnOnly(false);
		fHLTMUONDecoder.GetHandler().PrintParityErrorAsWarning(true);
		break;
	default:
		fRecoveryMode = kDontTryRecover;
		fHLTMUONDecoder.TryRecover(false);
		fHLTMUONDecoder.ExitOnError(true);
		fHLTMUONDecoder.GetHandler().WarnOnly(false);
		fHLTMUONDecoder.GetHandler().PrintParityErrorAsWarning(false);
		break;
	}
}


bool AliHLTMUONHitReconstructor::Run(
		const AliHLTUInt32_t* rawData,
		AliHLTUInt32_t rawDataSize,
		AliHLTMUONRecHitStruct* const recHit,
		AliHLTUInt32_t& nofHit
	) 
{
  // main function called by HLTReconstructor to perform DHLT Hitreconstruction 

  fRecPoints = recHit;
  fMaxRecPointsCount = nofHit;
  fRecPointsCount = &nofHit;
  *fRecPointsCount = 0;
  fDigitPerDDL = 0;
  fClusterCount = 0;
  fChannelCount = 0;
  HLTDebug("Decoding for DDL : %d",fDDL);
  if(GetkMinDetElemIdInDDL(fDDL) == -1 or GetkNofDetElemInDDL(fDDL)==-1){
    HLTError("DDL value fDDL : %d, out of range",fDDL);
  }
  if (not DecodeDDL(rawData, rawDataSize)) {
    // Dont need to log any message again. Already done so in DecodeDDL.
    return false;
  }

  if (fDigitPerDDL == 1)
  {
    // There are no digits to process so stop here.
    return true;
  }
  
  // Allocate fClusters and fChannels if required to do so and only if the allocated
  // size of the arrays is too small.
  try
  {
    if (fGenerateClusterInfo and fMaxClusters < fMaxRecPointsCount)
    {
      if (fClusters != NULL)
      {
        delete [] fClusters;
        fMaxClusters = 0;
      }
      fClusters = new AliHLTMUONClusterStruct[fMaxRecPointsCount];
      fMaxClusters = fMaxRecPointsCount;
    }
    if (fGenerateChannelInfo and fMaxChannels < fMaxRecPointsCount*fMaxChannelMult)
    {
      if (fChannels != NULL)
      {
        delete [] fChannels;
        fMaxChannels = 0;
      }
      fChannels = new AliHLTMUONChannelStruct[fMaxRecPointsCount*fMaxChannelMult];
      fMaxChannels = fMaxRecPointsCount*fMaxChannelMult;
    }
  }
  catch(const std::bad_alloc&)
  {
    HLTError("Could not allocate memory for the extra cluster and channel information.");
    return false;
  }

  if (not FindRecHits()) {
    HLTError("Failed to generate RecHits");
    return false;
  }

  return true;
}


bool AliHLTMUONHitReconstructor::FillClusterData(
		AliHLTMUONClusterStruct* clusters, AliHLTUInt32_t& nofClusters
	)
{
	/// Fills the output clusters array with extra cluster information.
	
	bool sizeOk = fClusterCount <= nofClusters;
	AliHLTUInt32_t n = sizeOk ? fClusterCount : nofClusters;
	memcpy(clusters, fClusters, sizeof(AliHLTMUONClusterStruct)*n);
	nofClusters = n;
	return sizeOk;
}


bool AliHLTMUONHitReconstructor::FillChannelData(
		AliHLTMUONChannelStruct* channels, AliHLTUInt32_t& nofChannels
	)
{
	/// Fills the output channels array with extra channel information for each cluster.
	
	bool sizeOk = fChannelCount <= nofChannels;
	AliHLTUInt32_t n = sizeOk ? fChannelCount : nofChannels;
	memcpy(channels, fChannels, sizeof(AliHLTMUONChannelStruct)*n);
	nofChannels = n;
	return sizeOk;
}


bool AliHLTMUONHitReconstructor::DecodeDDL(const AliHLTUInt32_t* rawData,AliHLTUInt32_t rawDataSize)
{
  //function to decode Raw Data 

  AliHLTMUONRawDecoder& handler = reinterpret_cast<AliHLTMUONRawDecoder&>(fHLTMUONDecoder.GetHandler());
  UInt_t bufferSize = UInt_t(rawDataSize*sizeof(AliHLTUInt32_t));

  handler.SetDCCut(fDCCut);
  handler.SetPadData(fPadData);
  handler.SetLookUpTable(fkLookUpTableData);
  handler.SetIdManuChannelToEntry(fkIdToEntry);
  handler.DDLNumber(fDDL);
  handler.SetNofFiredDetElemId(fNofDataInDetElem);
  handler.SetMaxFiredPerDetElem(fDataCountListPerDetElem);
  handler.SetMaxEntryPerBusPatch(fkMaxEntryPerBusPatch);
 
  if (not fHLTMUONDecoder.Decode(rawData,bufferSize))
  {
	switch (TryRecover())
	{
	case kRecoverFull:
		// Do not print the following warning for option "-dontprintparityerrors" if there
		// were only parity errors.
		if (fHLTMUONDecoder.GetHandler().NonParityErrorFound() or
		    (fHLTMUONDecoder.GetHandler().ParityErrorFound() and not fHLTMUONDecoder.GetHandler().DontPrintParityErrors())
		   )
		{
			HLTWarning("There was a problem with the raw data."
				" Recovered as much data as possible."
				" Will continue processing the next event."
			);
		}
		break;
	case kRecoverJustSkip:
		if (fHLTMUONDecoder.GetHandler().NonParityErrorFound() or
		    (fHLTMUONDecoder.GetHandler().ParityErrorFound() and not fHLTMUONDecoder.GetHandler().DontPrintParityErrors())
		   )
		{
			HLTWarning("There was a problem with the raw data."
				" Skipped corrupted data structures."
				" Will continue processing the next event."
			);
		}
		break;
	case kRecoverFromParityErrorsOnly:
		if (fHLTMUONDecoder.GetHandler().NonParityErrorFound())
		{
			HLTError("Failed to decode the tracker DDL raw data.");
			return false;
		}
		if (not fHLTMUONDecoder.GetHandler().DontPrintParityErrors())
		{
			assert( fHLTMUONDecoder.GetHandler().ParityErrorFound() );
			HLTWarning("Found parity errors in the raw data,"
				" but will continue processing."
			);
		}
		break;
	default:
		HLTError("Failed to decode the tracker DDL raw data.");
		return false;
	}
  }

  fDigitPerDDL = handler.GetDataCount();

  // fMaxFiredPerDetElem[fNofFiredDetElem-1] = handler.GetDataCount();
  
  // HLTDebug("fNofFiredDetElem : %d, NofDigits %d and max reco point limit is : %d, nofDetElems : %d",
  //         fNofFiredDetElem,fDigitPerDDL,fMaxRecPointsCount,fNofFiredDetElem);

  // for(int iDet=0; iDet<TMath::Max(fNofFiredDetElem,130); iDet++)
  //   HLTDebug("NofCount (fMaxFiredPerDetElem) in iDet %d is : %d", iDet, fMaxFiredPerDetElem[iDet]);
  
  // if(fNofFiredDetElem>129){
  //   HLTError("Number of fired detection elements is %d, which is more than 129.", fNofFiredDetElem);
  //   return false;
  // }
  
  if(fDigitPerDDL == 1){
    HLTDebug("An Empty DDL file was found.");
  }
  
  return true;
}


bool AliHLTMUONHitReconstructor::FindRecHits()
{
  // fuction that calls hit reconstruction detector element-wise.

  assert( fCentralChargeB == NULL );
  assert( fCentralChargeNB == NULL );
  assert( fRecX == NULL );
  assert( fRecY == NULL );
  assert( fAvgChargeX == NULL );
  assert( fAvgChargeY == NULL );
  assert( fTotChargeX == NULL );
  assert( fTotChargeY == NULL );
  assert( fNofBChannel == NULL );
  assert( fNofNBChannel == NULL );
  assert( fNofYNeighbour == NULL );
  
  bool resultOk = false;
  


  for(int iDet=0; iDet< GetkNofDetElemInDDL(fDDL) ; iDet++)
  {
    fCentralCountB = 0;
    fCentralCountNB = 0;

    try
    {
      fCentralChargeB = new int[fNofDataInDetElem[iDet]];
      HLTDebug("Allocated fCentralChargeB with %d elements.", fNofDataInDetElem[iDet]);
      fCentralChargeNB = new int[fNofDataInDetElem[iDet]];
      HLTDebug("Allocated fCentralChargeNB with %d elements.", fNofDataInDetElem[iDet]);
      resultOk = true;
    }
    catch(const std::bad_alloc&)
    {
      HLTError("Dynamic memory allocation failed for fCentralChargeNB and fCentralChargeB");
      resultOk = false;
      //break; Do not break. Might have smaller memory requirements in the next iteration.
    }

    // Continue processing, but check if everything is OK as we do, otherwise
    // do not execute the next steps.
    if (resultOk)
    {
      HLTDebug("Finding central hists for nofDigit  : %d, to in iDet  : %4d and min detelem : %4d",
                 fNofDataInDetElem[iDet],iDet,GetkMinDetElemIdInDDL(fDDL));
      FindCentralHits(iDet);
      
      HLTDebug("For iDet : %d, Found fCentralCountB : %d, fCentralCountNB : %d",iDet,fCentralCountB,fCentralCountNB);
      if(fCentralCountB==0 or fCentralCountNB==0)
      {
        HLTDebug("There is no fired pad in bending/nonbending plane...skipping this detection element");
        if (fCentralChargeB != NULL)
        {
          delete [] fCentralChargeB;
          HLTDebug("Released fCentralChargeB array.");
          fCentralChargeB = NULL;
        }
        if (fCentralChargeNB != NULL)
        {
          delete [] fCentralChargeNB;
          HLTDebug("Released fCentralChargeNB array.");
          fCentralChargeNB = NULL;
        }
        continue;
      }
    }

    if (resultOk)
    {
      try
      {
        fRecY = new float[fCentralCountB];
        HLTDebug("Allocated fRecY with %d elements.", fCentralCountB);
        fRecX = new float[fCentralCountNB];
        HLTDebug("Allocated fRecX with %d elements.", fCentralCountNB);
        fAvgChargeY = new float[fCentralCountB];
        HLTDebug("Allocated fAvgChargeY with %d elements.", fCentralCountB);
        fAvgChargeX = new float[fCentralCountNB];
        HLTDebug("Allocated fAvgChargeX with %d elements.", fCentralCountNB);
        fTotChargeY = new float[fCentralCountB];
        HLTDebug("Allocated fTotChargeY with %d elements.", fCentralCountB);
        fTotChargeX = new float[fCentralCountNB];
        HLTDebug("Allocated fTotChargeX with %d elements.", fCentralCountNB);
        fNofBChannel = new int[fCentralCountB];
        HLTDebug("Allocated fNofBChannel with %d elements.", fCentralCountB);
        fNofNBChannel = new int[fCentralCountNB];
        HLTDebug("Allocated fNofNBChannel with %d elements.", fCentralCountNB);
        fNofYNeighbour = new int[fCentralCountB];
        HLTDebug("Allocated fNofBChannel with %d elements.", fCentralCountB);
        resultOk = true;
      }
      catch(const std::bad_alloc&){
        HLTError("Dynamic memory allocation failed for internal arrays.");
        resultOk = false;
        //break; Must not break, this will prevent calling delete and memory cleanup, i.e. memory leak.
      }
    }

    if (resultOk) RecXRecY();
    
    if (resultOk)
    {
      if(fDDL<8)
	resultOk = MergeQuadRecHits();
      else
	resultOk = MergeSlatRecHits();
    }
    
    // minimum value in loop is 1 because dataCount in ReadDDL starts from 1 instead of 0;
    for(int i=0;i<fNofDataInDetElem[iDet];i++)
      fGetIdTotalData[fPadData[fDataCountListPerDetElem[iDet][i]].fIX]
	[fPadData[fDataCountListPerDetElem[iDet][i]].fIY]
	[fPadData[fDataCountListPerDetElem[iDet][i]].fPlane] = 0;

    // Make sure to release any memory that was allocated.
    if (fCentralChargeB != NULL)
    {
      delete [] fCentralChargeB;
      HLTDebug("Released fCentralChargeB array.");
      fCentralChargeB = NULL;
    }
    if (fCentralChargeNB != NULL)
    {
      delete [] fCentralChargeNB;
      HLTDebug("Released fCentralChargeNB array.");
      fCentralChargeNB = NULL;
    }
    if (fRecX != NULL)
    {
      delete [] fRecX;
      HLTDebug("Released fRecX array.");
      fRecX = NULL;
    }
    if (fRecY != NULL)
    {
      delete [] fRecY;
      HLTDebug("Released fRecY array.");
      fRecY = NULL;
    }
    if (fAvgChargeX != NULL)
    {
      delete [] fAvgChargeX;
      HLTDebug("Released fAvgChargeX array.");
      fAvgChargeX = NULL;
    }
    if (fAvgChargeY != NULL)
    {
      delete [] fAvgChargeY;
      HLTDebug("Released fAvgChargeY array.");
      fAvgChargeY = NULL;
    }
    if (fTotChargeX != NULL)
    {
      delete [] fTotChargeX;
      HLTDebug("Released fTotChargeX array.");
      fTotChargeX = NULL;
    }
    if (fTotChargeY != NULL)
    {
      delete [] fTotChargeY;
      HLTDebug("Released fTotChargeY array.");
      fTotChargeY = NULL;
    }
    if (fNofBChannel != NULL)
    {
      delete [] fNofBChannel;
      HLTDebug("Released fNofBChannel array.");
      fNofBChannel = NULL;
    }
    if (fNofNBChannel != NULL)
    {
      delete [] fNofNBChannel;
      HLTDebug("Released fNofNBChannel array.");
      fNofNBChannel = NULL;
    }
    if (fNofYNeighbour != NULL)
    {
      delete [] fNofYNeighbour;
      HLTDebug("Released fNofYNeighbour array.");
      fNofYNeighbour = NULL;
    }
  }

  Clear();  // clear internal arrays.

  return resultOk;
}


void AliHLTMUONHitReconstructor::FindCentralHits(int iDet)
{
  // to find central hit associated with each cluster

  assert( fCentralChargeB != NULL );
  assert( fCentralChargeNB != NULL );

  int b,nb;
  int idManuChannelCentral;
  bool hasFind;
  int iPad = 0;

  for(int iEntry=0;iEntry<fNofDataInDetElem[iDet];iEntry++){

    iPad = fDataCountListPerDetElem[iDet][iEntry];

    //if(fPadData[iPad].fDetElemId==102)
    HLTDebug("iPad : %d, detElem : %d, fCentralCountB : %d, fCentralCountNB : %d",iPad,fPadData[iPad].fDetElemId,fCentralCountB,fCentralCountNB);

    fGetIdTotalData[fPadData[iPad].fIX]
      [fPadData[iPad].fIY]
      [fPadData[iPad].fPlane] = iPad ;

    if(fPadData[iPad].fCharge <= fDCCut ) continue;
    
    if(fPadData[iPad].fPlane == 0 ){//&& fPadData[iPad].fIY > (0+1) && fPadData[iPad].fIY < (79 - 1)){
      //if(fPadData[iPad].fIY > 0){
      if(fCentralCountB>0){
	hasFind = false;
	for(b = 0;b<fCentralCountB;b++){
	  idManuChannelCentral = fCentralChargeB[b];
	  if(fPadData[iPad].fIX == fPadData[idManuChannelCentral].fIX
	     &&
	     (fPadData[iPad].fIY 
	      == fPadData[idManuChannelCentral].fIY + 1 
	      ||
	      fPadData[iPad].fIY 
	      == fPadData[idManuChannelCentral].fIY + 2 
	      ||
	      fPadData[iPad].fIY 
	      == fPadData[idManuChannelCentral].fIY - 2 
	      ||
	      fPadData[iPad].fIY 
	      == fPadData[idManuChannelCentral].fIY - 1)){
	    
	    hasFind = true;
	    if(fPadData[iPad].fCharge > fPadData[idManuChannelCentral].fCharge){
	      fCentralChargeB[b] = iPad;
	    }// if condn on pad charge
	  }// if condon on pad position
	}// for loop over b
	if(!hasFind){
	  fCentralChargeB[fCentralCountB] = iPad;
	  fCentralCountB++;
	}
      }
      else{
	fCentralChargeB[fCentralCountB] = iPad;
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
	  if(fPadData[iPad].fIY == fPadData[idManuChannelCentral].fIY
	     &&
	     (fPadData[iPad].fIX 
	      == fPadData[idManuChannelCentral].fIX + 1 
	      ||
	      fPadData[iPad].fIX
	      == fPadData[idManuChannelCentral].fIX + 2
	      ||
	      fPadData[iPad].fIX
	      == fPadData[idManuChannelCentral].fIX - 2
	      ||
	      fPadData[iPad].fIX
	      == fPadData[idManuChannelCentral].fIX - 1)){
	    
	    hasFind = true;	  
	    if(fPadData[iPad].fCharge > fPadData[idManuChannelCentral].fCharge){
	      fCentralChargeNB[nb] = iPad;
	    }// if condn over to find higher charge
	  }// if condn over to find position
	}// for loop over presently all nb values
	if(!hasFind){
	  fCentralChargeNB[fCentralCountNB] = iPad;
	  fCentralCountNB++;
	}
      }// centralHitNB size test
      else{
	fCentralChargeNB[fCentralCountNB] = iPad;
	fCentralCountNB++;
      }// centralHitNB size test
      
    }// fill for bending and nonbending hit
  }// detElemId loop

}

void AliHLTMUONHitReconstructor::RecXRecY()
{
  // find reconstructed X and Y for each plane separately

  assert( fRecX != NULL );
  assert( fRecY != NULL );
  assert( fAvgChargeX != NULL );
  assert( fAvgChargeY != NULL );
  assert( fTotChargeX != NULL );
  assert( fTotChargeY != NULL );
  assert( fNofBChannel != NULL );
  assert( fNofNBChannel != NULL );
  assert( fNofYNeighbour  != NULL );

  int b,nb;
  int idCentral;
  int idLower = 0;
  int idUpper = 0;
  int idRight = 0;
  int idLeft = 0;
  
  for(b=0;b<fCentralCountB;b++){
    idCentral = fCentralChargeB[b];

    if(fPadData[idCentral].fIY==0)
      idLower = 0;
    else
      idLower = fGetIdTotalData[fPadData[idCentral].fIX][fPadData[idCentral].fIY-1][0];
    
    if(fPadData[idCentral].fIY==236)
      idUpper = 0;
    else
      idUpper = fGetIdTotalData[fPadData[idCentral].fIX][fPadData[idCentral].fIY+1][0];

    
    fTotChargeY[b] = (fPadData[idCentral].fCharge + fPadData[idUpper].fCharge + fPadData[idLower].fCharge);
    if(fTotChargeY[b]==0.0) continue;
    fAvgChargeY[b] = fTotChargeY[b]/3.0 ;

    fRecY[b] = (fPadData[idCentral].fRealY*fPadData[idCentral].fCharge
	       +
		fPadData[idUpper].fRealY*fPadData[idUpper].fCharge
	       +
		fPadData[idLower].fRealY*fPadData[idLower].fCharge
		)/fTotChargeY[b];//(fPadData[idCentral].fCharge + fPadData[idUpper].fCharge + fPadData[idLower].fCharge) ;
    
    fNofBChannel[b] = 0;
    fNofYNeighbour[b] = 0;
    if(fPadData[idLower].fCharge>0.0){
      fNofBChannel[b]++ ;
      fNofYNeighbour[b]++;
    }
    if(fPadData[idCentral].fCharge>0.0)
      fNofBChannel[b]++ ;
    if(fPadData[idUpper].fCharge>0.0){
      fNofBChannel[b]++ ;
      fNofYNeighbour[b]++;
    }

    HLTDebug("detelem : %d, Y charge : lower : %f, middle : %f, upper : %f",fPadData[idCentral].fDetElemId,fPadData[idLower].fCharge,
	    fPadData[idCentral].fCharge,fPadData[idUpper].fCharge);

    //collect left coloumn
    if((fPadData[idCentral].fIX-1)>=0){
      
      idLeft = fGetIdTotalData[fPadData[idCentral].fIX-1][fPadData[idCentral].fIY][0];
      
      if(fPadData[idLeft].fIY==0)
	idLower = 0;
      else
	idLower = fGetIdTotalData[fPadData[idLeft].fIX][fPadData[idLeft].fIY-1][0];
      
      if(fPadData[idLeft].fIY==236)
	idUpper = 0;
      else
	idUpper = fGetIdTotalData[fPadData[idLeft].fIX][fPadData[idLeft].fIY+1][0];

      fTotChargeY[b] += (fPadData[idLeft].fCharge + fPadData[idUpper].fCharge + fPadData[idLower].fCharge);

      if(fPadData[idLower].fCharge>0.0)
	fNofBChannel[b]++ ;
      if(fPadData[idLeft].fCharge>0.0)
	fNofBChannel[b]++ ;
      if(fPadData[idUpper].fCharge>0.0)
	fNofBChannel[b]++ ;

    }
    ////////////////////////////////////////////////////

    //collect right coloumn
    if((fPadData[idCentral].fIX+1)<=335){

      idRight = fGetIdTotalData[fPadData[idCentral].fIX+1][fPadData[idCentral].fIY][0];
    
      if(fPadData[idRight].fIY==0)
	idLower = 0;
      else
	idLower = fGetIdTotalData[fPadData[idRight].fIX][fPadData[idRight].fIY-1][0];
      
      if(fPadData[idRight].fIY==236)
	idUpper = 0;
      else
	idUpper = fGetIdTotalData[fPadData[idRight].fIX][fPadData[idRight].fIY+1][0];
  
      fTotChargeY[b] += (fPadData[idRight].fCharge + fPadData[idUpper].fCharge + fPadData[idLower].fCharge);

      if(fPadData[idLower].fCharge>0.0)
	fNofBChannel[b]++ ;
      if(fPadData[idRight].fCharge>0.0)
	fNofBChannel[b]++ ;
      if(fPadData[idUpper].fCharge>0.0)
	fNofBChannel[b]++ ;

    }
    //////////////////////////////////////////////////////////////////////////////////
    HLTDebug("RecY[%d] : %f, nofChannel : %d, detelem : %d",b,fRecY[b],fNofBChannel[b],fPadData[idCentral].fDetElemId);
   
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

    fTotChargeX[nb] = (fPadData[idCentral].fCharge + fPadData[idRight].fCharge + fPadData[idLeft].fCharge);
    if(fTotChargeX[nb]==0.0) continue;
    fAvgChargeX[nb] = fTotChargeX[nb]/3.0 ;
    
    fRecX[nb] = (fPadData[idCentral].fRealX*fPadData[idCentral].fCharge
		 +
		 fPadData[idRight].fRealX*fPadData[idRight].fCharge
		 +
		 fPadData[idLeft].fRealX*fPadData[idLeft].fCharge
		 )/fTotChargeX[nb];//(fPadData[idCentral].fCharge + fPadData[idRight].fCharge + fPadData[idLeft].fCharge);
    

    fNofNBChannel[nb] = 0;
    if(fPadData[idLeft].fCharge>0.0)
      fNofNBChannel[nb]++ ;
    if(fPadData[idCentral].fCharge>0.0)
      fNofNBChannel[nb]++ ;
    if(fPadData[idRight].fCharge>0.0)
      fNofNBChannel[nb]++ ;

    HLTDebug("detelem : %d, X charge left : %f, middle : %f, right : %f",fPadData[idCentral].fDetElemId,fPadData[idLeft].fCharge,
	    fPadData[idCentral].fCharge,fPadData[idRight].fCharge);


    // lower row 
    if((fPadData[idCentral].fIY-1)>=0){

      idLower = fGetIdTotalData[fPadData[idCentral].fIX][fPadData[idCentral].fIY-1][1];
      
      if(fPadData[idLower].fIX==0)
	idLeft = 0;
      else
	idLeft = fGetIdTotalData[fPadData[idLower].fIX-1][fPadData[idLower].fIY][1];
    
      if(fPadData[idLower].fIX==335)
	idRight = 0 ;
      else
	idRight = fGetIdTotalData[fPadData[idLower].fIX+1][fPadData[idLower].fIY][1];

      fTotChargeX[nb] += (fPadData[idLower].fCharge + fPadData[idRight].fCharge + fPadData[idLeft].fCharge);

      if(fPadData[idLeft].fCharge>0.0)
	fNofNBChannel[nb]++ ;
      if(fPadData[idLower].fCharge>0.0)
	fNofNBChannel[nb]++ ;
      if(fPadData[idRight].fCharge>0.0)
	fNofNBChannel[nb]++ ;

    }
    ////////////////////////////////////////////////////////////

    // Upper row
    if((fPadData[idCentral].fIY+1)<=236){

      idUpper = fGetIdTotalData[fPadData[idCentral].fIX][fPadData[idCentral].fIY+1][1];

      if(fPadData[idUpper].fIX==0)
	idLeft = 0;
      else
	idLeft = fGetIdTotalData[fPadData[idUpper].fIX-1][fPadData[idUpper].fIY][1];
      
      if(fPadData[idUpper].fIX==335)
	idRight = 0 ;
      else
	idRight = fGetIdTotalData[fPadData[idUpper].fIX+1][fPadData[idUpper].fIY][1];
      
      fTotChargeX[nb] += (fPadData[idUpper].fCharge + fPadData[idRight].fCharge + fPadData[idLeft].fCharge);

      if(fPadData[idLeft].fCharge>0.0)
	fNofNBChannel[nb]++ ;
      if(fPadData[idRight].fCharge>0.0)
	fNofNBChannel[nb]++ ;
      if(fPadData[idRight].fCharge>0.0)
	fNofNBChannel[nb]++ ;

    }
    ////////////////////////////////////////////////////////////

    HLTDebug("RecX[%d] : %f, nofChannel : %d",nb,fRecX[nb],fNofNBChannel[nb]);

   

  }
}


bool AliHLTMUONHitReconstructor::MergeQuadRecHits()
{
  // Merge reconstructed hits first over same plane then bending plane with non-bending plane

  assert( fRecX != NULL );
  assert( fRecY != NULL );
  assert( fAvgChargeX != NULL );
  assert( fAvgChargeY != NULL );
  assert( fTotChargeX != NULL );
  assert( fTotChargeY != NULL );
  assert( fNofBChannel != NULL );
  assert( fNofNBChannel != NULL );
  assert( fNofYNeighbour  != NULL );

  int idCentralB=0,idCentralNB=0 ;
  float padCenterXB;
  float padCenterYNB;
  float diffX,diffY;
  float minPadArea;
  float halfPadLengthX,halfPadLengthY;
  bool *isMergedY = new bool[fCentralCountB];
  bool *isMergedX = new bool[fCentralCountNB];

  // MERGE Bending Plane hits, which are placed side by side
  for(int i=0;i<fCentralCountB-1;i++){
    isMergedY[i] = false;
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
  isMergedY[fCentralCountB-1] = false;

  // MERGE Non Bending Plane hits, which are placed side by side
  for(int i=0;i<fCentralCountNB-1;i++){
    isMergedX[i] = false;
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
  isMergedX[fCentralCountNB-1] = false;
  
  // Merge bending Plane hits with Non Bending
  for(int b=0;b<fCentralCountB;b++){
    if(fRecY[b]!=0.0){
      idCentralB = fCentralChargeB[b];
      padCenterXB = fPadData[idCentralB].fRealX; 
      
      halfPadLengthX = fPadData[idCentralB].fIY ;

      for(int nb=0;nb<fCentralCountNB;nb++){
	if(fRecX[nb]!=0.0){
	  idCentralNB = fCentralChargeNB[nb];

	  padCenterYNB = fPadData[idCentralNB].fRealY;

	  halfPadLengthY = fPadData[idCentralNB].fHalfPadSize ;

	  if(fabsf(fRecX[nb]) > fabsf(padCenterXB))
	    diffX = fabsf(fRecX[nb]) -  fabsf(padCenterXB);
	  else
	    diffX = fabsf(padCenterXB) -  fabsf(fRecX[nb]);
	  
	  if(fabsf(padCenterYNB)>fabsf(fRecY[b]))
	    diffY = fabsf(padCenterYNB) - fabsf(fRecY[b]);
	  else
	    diffY =  fabsf(fRecY[b]) - fabsf(padCenterYNB);

	  if(diffX < halfPadLengthX && diffY < halfPadLengthY ){//&& fPadData[idCentralB].fIY != 0){
	    
	    isMergedY[b] = true;
	    isMergedX[nb] = true;
	    
	    if(fNofYNeighbour[b]==2){
	      if(fPadData[idCentralB].fDetElemId<104)
		fRecY[b] += 0.02*sin(14.5*(fRecY[b] - fPadData[idCentralB].fRealY)) ;
	      else if(fPadData[idCentralB].fDetElemId>=200 && fPadData[idCentralB].fDetElemId<204)
		fRecY[b] += 0.02*sin(14.0*(fRecY[b] - fPadData[idCentralB].fRealY)) ;
	      else
		fRecY[b] += 0.025*sin(12.0*(fRecY[b] - fPadData[idCentralB].fRealY)) ;
	    }
	    
	    if(fPadData[idCentralNB].fDetElemId<204)
	      fRecX[nb] += 0.095*sin(10.5*(fRecX[nb] - fPadData[idCentralNB].fRealX)) ;
	    else //if(fPadData[idCentralNB].fDetElemId>=300 && fPadData[idCentralNB].fDetElemId<404)
	      fRecX[nb] += 0.085*sin(9.0*(fRecX[nb] - fPadData[idCentralNB].fRealX)) ;
	    
	    
	    // First check that we have not overflowed the buffer.
	    if((*fRecPointsCount) == fMaxRecPointsCount){
	      HLTWarning("Number of RecHits (i.e. %d) exceeds the max number of RecHit limit %d."
	                " Output buffer is too small.",
	               (*fRecPointsCount),fMaxRecPointsCount
	      );
	      delete [] isMergedY;
	      delete [] isMergedX;
	      return true;
	    }

	    AliHLTUInt32_t idflags = AliHLTMUONUtils::PackRecHitFlags(
	         (fPadData[idCentralB].fDetElemId / 100) - 1,
	         fPadData[idCentralB].fDetElemId
	      );
	    fRecPoints[(*fRecPointsCount)].fFlags = idflags;
	    fRecPoints[(*fRecPointsCount)].fX = fRecX[nb];
	    fRecPoints[(*fRecPointsCount)].fY = fRecY[b];
	    fRecPoints[(*fRecPointsCount)].fZ = fPadData[idCentralB].fRealZ;
	    
 	    if (fGenerateClusterInfo)
 	    {
 	      if (fClusterCount >= fMaxClusters)
 	      {
	        HLTError("Ran out of space in internal cluster array of size %d.", fMaxClusters);
		delete [] isMergedY;
		delete [] isMergedX;
	        return false;
	      }
	      
	      fClusters[fClusterCount].fId = (fNewClusterId << 5) | fDDL;
	      
	      // Increment the cluster ID and warp it around at 0x03FFFFFF since
	      // the bottom 5 bits are filled with the source DDL number and the
	      // sign bit in fClusters[fClusterCount].fId must be positive.
	      fNewClusterId = (fNewClusterId + 1) & 0x03FFFFFF;
	      
	      fClusters[fClusterCount].fHit = fRecPoints[(*fRecPointsCount)];
	      fClusters[fClusterCount].fDetElemId = fPadData[idCentralB].fDetElemId;
 	      fClusters[fClusterCount].fNchannelsB = fNofBChannel[b];
 	      fClusters[fClusterCount].fNchannelsNB = fNofNBChannel[nb];
 	      fClusters[fClusterCount].fChargeB = fTotChargeY[b];
 	      fClusters[fClusterCount].fChargeNB = fTotChargeX[nb];
 	      fClusterCount++;
 	    }
 	    
 	    if (fGenerateChannelInfo)
 	    {
 	      // 3 by 3 pad structure around the central pad for the 2 planes.
 	      int pad[2][3][3] = {{
 	          {0,      0,      0},
 	          {0, idCentralB,  0},
 	          {0,      0,      0}
 	        },{
 	          {0,      0,      0},
 	          {0, idCentralNB, 0},
 	          {0,      0,      0}
 	      }};
              
              // Find the pad index numbers for the central pads and the pads surrounding them.
              // All these pads would have contributed to the cluster as long as their charge is != 0.
              if (fPadData[idCentralB].fIY > 0)
                pad[0][0][1] = fGetIdTotalData[fPadData[idCentralB].fIX][fPadData[idCentralB].fIY-1][0];
              if (fPadData[idCentralB].fIY < 236)
                pad[0][2][1] = fGetIdTotalData[fPadData[idCentralB].fIX][fPadData[idCentralB].fIY+1][0];
              if (fPadData[idCentralB].fIX > 0)
              {
                int idLeft = fGetIdTotalData[fPadData[idCentralB].fIX-1][fPadData[idCentralB].fIY][0];
                pad[0][1][0] = idLeft;
                if (fPadData[idLeft].fIY > 0)
                  pad[0][0][0] = fGetIdTotalData[fPadData[idLeft].fIX][fPadData[idLeft].fIY-1][0];
                if (fPadData[idLeft].fIY < 236)
                  pad[0][2][0] = fGetIdTotalData[fPadData[idLeft].fIX][fPadData[idLeft].fIY+1][0];
              }
              if (fPadData[idCentralB].fIX < 335)
              {
                int idRight = fGetIdTotalData[fPadData[idCentralB].fIX+1][fPadData[idCentralB].fIY][0];
                pad[0][1][2] = idRight;
                if (fPadData[idRight].fIY > 0)
                  pad[0][0][2] = fGetIdTotalData[fPadData[idRight].fIX][fPadData[idRight].fIY-1][0];
                if (fPadData[idRight].fIY < 236)
                  pad[0][2][2] = fGetIdTotalData[fPadData[idRight].fIX][fPadData[idRight].fIY+1][0];
              }
              
              if (fPadData[idCentralNB].fIX > 0)
                pad[1][1][0] = fGetIdTotalData[fPadData[idCentralNB].fIX-1][fPadData[idCentralNB].fIY][1];
              if (fPadData[idCentralNB].fIX < 335)
                pad[1][1][2] = fGetIdTotalData[fPadData[idCentralNB].fIX+1][fPadData[idCentralNB].fIY][1];
              if (fPadData[idCentralNB].fIY > 0)
              {
                int idLower = fGetIdTotalData[fPadData[idCentralNB].fIX][fPadData[idCentralNB].fIY-1][1];
                pad[1][0][1] = idLower;
                if (fPadData[idLower].fIX > 0)
                  pad[1][0][0] = fGetIdTotalData[fPadData[idLower].fIX-1][fPadData[idLower].fIY][1];
                if (fPadData[idLower].fIX < 335)
                  pad[1][0][2] = fGetIdTotalData[fPadData[idLower].fIX+1][fPadData[idLower].fIY][1];
              }
              if (fPadData[idCentralNB].fIY < 236)
              {
                int idUpper = fGetIdTotalData[fPadData[idCentralNB].fIX][fPadData[idCentralNB].fIY+1][1];
                pad[1][2][1] = idUpper;
                if (fPadData[idUpper].fIX > 0)
                  pad[1][2][0] = fGetIdTotalData[fPadData[idUpper].fIX-1][fPadData[idUpper].fIY][1];
                if (fPadData[idUpper].fIX < 335)
                  pad[1][2][2] = fGetIdTotalData[fPadData[idUpper].fIX+1][fPadData[idUpper].fIY][1];
              }
              
              AliHLTInt32_t clusterId = fGenerateClusterInfo ? fClusters[fClusterCount-1].fId : -1;
              
              // Now generate the pad structures from all the pad indices found above.
              for (int i = 0; i < 2; i++)
              for (int j = 0; j < 3; j++)
              for (int k = 0; k < 3; k++)
              {
                AliHLTMUONPad& p = fPadData[pad[i][j][k]];
                // Skip pads that have zero charge because they would not have
                // contributed to the cluster.
                if (p.fCharge <= 0) continue;
                
                UShort_t manuId; UChar_t channelId; UShort_t adc;
                AliHLTMUONRawDecoder::UnpackADC(p.fRawData, manuId, channelId, adc);
                
                fChannels[fChannelCount].fClusterId = clusterId;
                fChannels[fChannelCount].fBusPatch = p.fBusPatch;
                fChannels[fChannelCount].fManu = manuId;
                fChannels[fChannelCount].fChannelAddress = channelId;
                fChannels[fChannelCount].fSignal = adc;
                fChannels[fChannelCount].fRawDataWord = p.fRawData;
                fChannelCount++;
              }
 	    }

	    HLTDebug("Part 1 : Reconstructed hit (X,Y,Z) : (%f,%f,%f)",
                     fRecPoints[(*fRecPointsCount)].fX,
                     fRecPoints[(*fRecPointsCount)].fY,
                     fRecPoints[(*fRecPointsCount)].fZ
            );
	    (*fRecPointsCount)++;
	  }//if lies wihtin 5.0 mm
	}// condn over fRecX ! = 0.0
      }// loop over NB side
    }// condn on fRecY[b] !=  0.0
  }// loop over B side;

  //Hit only in bending plane and in zone 1 which has not merged and which has number of channels > 2 are considered as valid hits
  for(int b=0;b<fCentralCountB;b++){
    if(fRecY[b]!=0.0 and !isMergedY[b] and fNofBChannel[b]>2){
      idCentralB = fCentralChargeB[b];
      
      minPadArea = (fPadData[idCentralB].fDetElemId < 204) ? 10.0*2.0*0.315*10.0*2.0*0.21 : 10.0*2.0*0.375*10.0*2.0*0.25 ;

      if(TMath::Abs(400.0*fPadData[idCentralB].fHalfPadSize*fPadData[idCentralB].fPadSizeXY - minPadArea)> 1.0e-5) continue;
      
      if(fNofYNeighbour[b]==2){
	if(fPadData[idCentralB].fDetElemId<104)
	  fRecY[b] += 0.02*sin(14.5*(fRecY[b] - fPadData[idCentralB].fRealY)) ;
	else if(fPadData[idCentralB].fDetElemId>=200 && fPadData[idCentralB].fDetElemId<204)
	  fRecY[b] += 0.02*sin(14.0*(fRecY[b] - fPadData[idCentralB].fRealY)) ;
	else
	  fRecY[b] += 0.025*sin(12.0*(fRecY[b] - fPadData[idCentralB].fRealY)) ;
      }
      
      padCenterXB = fPadData[idCentralB].fRealX; 
      if(fPadData[idCentralB].fDetElemId<204)
	padCenterXB += 0.095*sin(10.5*(padCenterXB - fPadData[idCentralB].fRealX)) ;
      else //if(fPadData[idCentralNB].fDetElemId>=300 && fPadData[idCentralNB].fDetElemId<404)
	padCenterXB += 0.085*sin(9.0*(padCenterXB - fPadData[idCentralB].fRealX)) ;
      	    
      // First check that we have not overflowed the buffer.
      if((*fRecPointsCount) == fMaxRecPointsCount){
	HLTWarning("Number of RecHits (i.e. %d) exceeds the max number of RecHit limit %d."
		 " Output buffer is too small.",
		 (*fRecPointsCount),fMaxRecPointsCount
		 );
	delete [] isMergedY;
	delete [] isMergedX;
	return true;
      }
      
      AliHLTUInt32_t idflags = AliHLTMUONUtils::PackRecHitFlags(
	         (fPadData[idCentralB].fDetElemId / 100) - 1,
	         fPadData[idCentralB].fDetElemId
	      );
      fRecPoints[(*fRecPointsCount)].fFlags = idflags;
      fRecPoints[(*fRecPointsCount)].fX = padCenterXB;
      fRecPoints[(*fRecPointsCount)].fY = fRecY[b];
      fRecPoints[(*fRecPointsCount)].fZ = fPadData[idCentralB].fRealZ;
	    
      if (fGenerateClusterInfo)
	{
	  if (fClusterCount >= fMaxClusters)
	    {
	      HLTError("Ran out of space in internal cluster array of size %d.", fMaxClusters);
	      delete [] isMergedY;
	      delete [] isMergedX;
	      return false;
	    }
	      
	  fClusters[fClusterCount].fId = (fNewClusterId << 5) | fDDL;
	      
	  // Increment the cluster ID and warp it around at 0x03FFFFFF since
	  // the bottom 5 bits are filled with the source DDL number and the
	  // sign bit in fClusters[fClusterCount].fId must be positive.
	  fNewClusterId = (fNewClusterId + 1) & 0x03FFFFFF;
	  
	  fClusters[fClusterCount].fHit = fRecPoints[(*fRecPointsCount)];
	  fClusters[fClusterCount].fDetElemId = fPadData[idCentralB].fDetElemId;
	  fClusters[fClusterCount].fNchannelsB = fNofBChannel[b];
	  fClusters[fClusterCount].fNchannelsNB = fNofBChannel[b];
	  fClusters[fClusterCount].fChargeB = fTotChargeY[b];
	  fClusters[fClusterCount].fChargeNB = fTotChargeY[b];
	  fClusterCount++;
	}
      
      if (fGenerateChannelInfo)
	{
	  // 3 by 3 pad structure around the central pad for the 2 planes.
	  int pad[2][3][3] = {{
	      {0,      0,      0},
	      {0, idCentralB,  0},
	      {0,      0,      0}
	    },{
	      {0,      0,      0},
	      {0,      0     , 0},
	      {0,      0,      0}
	    }};
	  
	  // Find the pad index numbers for the central pads and the pads surrounding them.
	  // All these pads would have contributed to the cluster as long as their charge is != 0.
	  if (fPadData[idCentralB].fIY > 0)
	    pad[0][0][1] = fGetIdTotalData[fPadData[idCentralB].fIX][fPadData[idCentralB].fIY-1][0];
	  if (fPadData[idCentralB].fIY < 236)
	    pad[0][2][1] = fGetIdTotalData[fPadData[idCentralB].fIX][fPadData[idCentralB].fIY+1][0];
	  if (fPadData[idCentralB].fIX > 0)
	    {
	      int idLeft = fGetIdTotalData[fPadData[idCentralB].fIX-1][fPadData[idCentralB].fIY][0];
	      pad[0][1][0] = idLeft;
	      if (fPadData[idLeft].fIY > 0)
		pad[0][0][0] = fGetIdTotalData[fPadData[idLeft].fIX][fPadData[idLeft].fIY-1][0];
	      if (fPadData[idLeft].fIY < 236)
		pad[0][2][0] = fGetIdTotalData[fPadData[idLeft].fIX][fPadData[idLeft].fIY+1][0];
	    }
	  if (fPadData[idCentralB].fIX < 335)
	    {
	      int idRight = fGetIdTotalData[fPadData[idCentralB].fIX+1][fPadData[idCentralB].fIY][0];
	      pad[0][1][2] = idRight;
	      if (fPadData[idRight].fIY > 0)
		pad[0][0][2] = fGetIdTotalData[fPadData[idRight].fIX][fPadData[idRight].fIY-1][0];
	      if (fPadData[idRight].fIY < 236)
		pad[0][2][2] = fGetIdTotalData[fPadData[idRight].fIX][fPadData[idRight].fIY+1][0];
	    }

	  // For hits only in bending plane no need to copy the information to the non-bending side
	  // for(int i=0;i<3;i++)
	  //   for(int j=0;j<3;j++)
	  //     pad[1][i][j] = pad[0][i][j];

              
	  AliHLTInt32_t clusterId = fGenerateClusterInfo ? fClusters[fClusterCount-1].fId : -1;
          
	  // Now generate the pad structures from all the pad indices found above.
	  for (int i = 0; i < 1; i++)
	    for (int j = 0; j < 3; j++)
              for (int k = 0; k < 3; k++)
		{
		  AliHLTMUONPad& p = fPadData[pad[i][j][k]];
		  // Skip pads that have zero charge because they would not have
		  // contributed to the cluster.
		  if (p.fCharge <= 0) continue;
		  
		  UShort_t manuId; UChar_t channelId; UShort_t adc;
		  AliHLTMUONRawDecoder::UnpackADC(p.fRawData, manuId, channelId, adc);
		  
		  fChannels[fChannelCount].fClusterId = clusterId;
		  fChannels[fChannelCount].fBusPatch = p.fBusPatch;
		  fChannels[fChannelCount].fManu = manuId;
		  fChannels[fChannelCount].fChannelAddress = channelId;
		  fChannels[fChannelCount].fSignal = adc;
		  fChannels[fChannelCount].fRawDataWord = p.fRawData;
		  fChannelCount++;
		}
	}
      
      HLTDebug("Part 2 : Reconstructed hit (X,Y,Z) : (%f,%f,%f)",
	       fRecPoints[(*fRecPointsCount)].fX,
	       fRecPoints[(*fRecPointsCount)].fY,
	       fRecPoints[(*fRecPointsCount)].fZ
	       );
      (*fRecPointsCount)++;

    }// condn on fRecY[b] !=  0.0
  }// loop over B side;





  //Hit only in non-bending plane and in zone 1 which has not merged and which has number of channels > 2 are considered as valid hits
  for(int nb=0;nb<fCentralCountNB;nb++){

    idCentralNB = fCentralChargeNB[nb];
    if(fRecX[nb]!=0.0 and !isMergedX[nb] and fNofNBChannel[nb]>2){
      
      minPadArea = (fPadData[idCentralNB].fDetElemId < 204) ? 10.0*2.0*0.315*10.0*2.0*0.21 : 10.0*2.0*0.375*10.0*2.0*0.25 ;

      if(TMath::Abs(400.0*fPadData[idCentralNB].fHalfPadSize*fPadData[idCentralNB].fPadSizeXY - minPadArea)> 1.0e-5) continue;
      
      padCenterYNB = fPadData[idCentralNB].fRealY;
      
      
      if(fPadData[idCentralNB].fDetElemId<104)
	padCenterYNB += 0.02*sin(14.5*(padCenterYNB - fPadData[idCentralNB].fRealY)) ;
      else if(fPadData[idCentralNB].fDetElemId>=200 && fPadData[idCentralNB].fDetElemId<204)
	padCenterYNB += 0.02*sin(14.0*(padCenterYNB - fPadData[idCentralNB].fRealY)) ;
      else
	padCenterYNB += 0.025*sin(12.0*(padCenterYNB - fPadData[idCentralNB].fRealY)) ;
      
      
      if(fPadData[idCentralNB].fDetElemId<204)
	fRecX[nb] += 0.095*sin(10.5*(fRecX[nb] - fPadData[idCentralNB].fRealX)) ;
      else //if(fPadData[idCentralNB].fDetElemId>=300 && fPadData[idCentralNB].fDetElemId<404)
	fRecX[nb] += 0.085*sin(9.0*(fRecX[nb] - fPadData[idCentralNB].fRealX)) ;
      
      
      // First check that we have not overflowed the buffer.
      if((*fRecPointsCount) == fMaxRecPointsCount){
	HLTWarning("Number of RecHits (i.e. %d) exceeds the max number of RecHit limit %d."
		 " Output buffer is too small.",
		 (*fRecPointsCount),fMaxRecPointsCount
		 );
	delete [] isMergedY;
	delete [] isMergedX;
	return true;
      }

      AliHLTUInt32_t idflags = AliHLTMUONUtils::PackRecHitFlags(
	         (fPadData[idCentralNB].fDetElemId / 100) - 1,
	         fPadData[idCentralNB].fDetElemId
	      );
      fRecPoints[(*fRecPointsCount)].fFlags = idflags;
      fRecPoints[(*fRecPointsCount)].fX = fRecX[nb];
      fRecPoints[(*fRecPointsCount)].fY = padCenterYNB;
      fRecPoints[(*fRecPointsCount)].fZ = fPadData[idCentralNB].fRealZ;
      
      if (fGenerateClusterInfo)
	{
	  if (fClusterCount >= fMaxClusters)
	    {
	      HLTError("Ran out of space in internal cluster array of size %d.", fMaxClusters);
	      delete [] isMergedY;
	      delete [] isMergedX;
	      return false;
	    }
	  
	  fClusters[fClusterCount].fId = (fNewClusterId << 5) | fDDL;
	  
	  // Increment the cluster ID and warp it around at 0x03FFFFFF since
	  // the bottom 5 bits are filled with the source DDL number and the
	  // sign bit in fClusters[fClusterCount].fId must be positive.
	  fNewClusterId = (fNewClusterId + 1) & 0x03FFFFFF;
	      
	  fClusters[fClusterCount].fHit = fRecPoints[(*fRecPointsCount)];
	  fClusters[fClusterCount].fDetElemId = fPadData[idCentralNB].fDetElemId;
	  fClusters[fClusterCount].fNchannelsB = fNofNBChannel[nb];
	  fClusters[fClusterCount].fNchannelsNB = fNofNBChannel[nb];
	  fClusters[fClusterCount].fChargeB = fTotChargeX[nb];
	  fClusters[fClusterCount].fChargeNB = fTotChargeX[nb];
	  fClusterCount++;
	}
      
      if (fGenerateChannelInfo)
	{
	  // 3 by 3 pad structure around the central pad for the 2 planes.
	  int pad[2][3][3] = {{
	      {0,      0,      0},
	      {0,      0,      0},
	      {0,      0,      0}
 	        },{
	      {0,      0,      0},
	      {0, idCentralNB, 0},
	      {0,      0,      0}
	    }};
              
	  // Find the pad index numbers for the central pads and the pads surrounding them.
	  // All these pads would have contributed to the cluster as long as their charge is != 0.
              
	  if (fPadData[idCentralNB].fIX > 0)
	    pad[1][1][0] = fGetIdTotalData[fPadData[idCentralNB].fIX-1][fPadData[idCentralNB].fIY][1];
	  if (fPadData[idCentralNB].fIX < 335)
	    pad[1][1][2] = fGetIdTotalData[fPadData[idCentralNB].fIX+1][fPadData[idCentralNB].fIY][1];
	  if (fPadData[idCentralNB].fIY > 0)
	    {
	      int idLower = fGetIdTotalData[fPadData[idCentralNB].fIX][fPadData[idCentralNB].fIY-1][1];
	      pad[1][0][1] = idLower;
	      if (fPadData[idLower].fIX > 0)
		pad[1][0][0] = fGetIdTotalData[fPadData[idLower].fIX-1][fPadData[idLower].fIY][1];
	      if (fPadData[idLower].fIX < 335)
		pad[1][0][2] = fGetIdTotalData[fPadData[idLower].fIX+1][fPadData[idLower].fIY][1];
	    }
	  if (fPadData[idCentralNB].fIY < 236)
	    {
	      int idUpper = fGetIdTotalData[fPadData[idCentralNB].fIX][fPadData[idCentralNB].fIY+1][1];
                pad[1][2][1] = idUpper;
                if (fPadData[idUpper].fIX > 0)
                  pad[1][2][0] = fGetIdTotalData[fPadData[idUpper].fIX-1][fPadData[idUpper].fIY][1];
                if (fPadData[idUpper].fIX < 335)
                  pad[1][2][2] = fGetIdTotalData[fPadData[idUpper].fIX+1][fPadData[idUpper].fIY][1];
	    }

	  // For hits only in non-bending plane no need to copy the information to the bending side
	  // for(int i=0;i<3;i++)
	  //   for(int j=0;j<3;j++)
	  //     pad[0][i][j] = pad[1][i][j];
	  
              
	  AliHLTInt32_t clusterId = fGenerateClusterInfo ? fClusters[fClusterCount-1].fId : -1;
	  
	  // Now generate the pad structures from all the pad indices found above.
	  for (int i = 1; i < 2; i++)
	    for (int j = 0; j < 3; j++)
              for (int k = 0; k < 3; k++)
		{
		  AliHLTMUONPad& p = fPadData[pad[i][j][k]];
		  // Skip pads that have zero charge because they would not have
		  // contributed to the cluster.
		  if (p.fCharge <= 0) continue;
		  
		  UShort_t manuId; UChar_t channelId; UShort_t adc;
		  AliHLTMUONRawDecoder::UnpackADC(p.fRawData, manuId, channelId, adc);
		  
		  fChannels[fChannelCount].fClusterId = clusterId;
		  fChannels[fChannelCount].fBusPatch = p.fBusPatch;
		  fChannels[fChannelCount].fManu = manuId;
		  fChannels[fChannelCount].fChannelAddress = channelId;
		  fChannels[fChannelCount].fSignal = adc;
		  fChannels[fChannelCount].fRawDataWord = p.fRawData;
		  fChannelCount++;
		}
	}
      
      HLTDebug("Part 3 : Reconstructed hit (X,Y,Z) : (%f,%f,%f), detelemId : %d",
		 fRecPoints[(*fRecPointsCount)].fX,
		 fRecPoints[(*fRecPointsCount)].fY,
		 fRecPoints[(*fRecPointsCount)].fZ,
		 fPadData[idCentralNB].fDetElemId
		 );
      (*fRecPointsCount)++;
    }//if lies wihtin 5.0 mm
  }// condn over fRecX ! = 0.0

  delete [] isMergedY;
  delete [] isMergedX;
  
  return true;
}

bool AliHLTMUONHitReconstructor::MergeSlatRecHits()
{
  // Merge reconstructed hits first over same plane then bending plane with non-bending plane

  assert( fRecX != NULL );
  assert( fRecY != NULL );
  assert( fAvgChargeX != NULL );
  assert( fAvgChargeY != NULL );
  assert( fTotChargeX != NULL );
  assert( fTotChargeY != NULL );
  assert( fNofBChannel != NULL );
  assert( fNofNBChannel != NULL );
  assert( fNofYNeighbour  != NULL );

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
  
  // MERGE Non Bending Plane hits, which are placed side by side
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
      
      halfPadLengthX = fPadData[idCentralB].fIY ;

      for(int nb=0;nb<fCentralCountNB;nb++){
	if(fRecX[nb]!=0.0){
	  idCentralNB = fCentralChargeNB[nb];

	  padCenterYNB = fPadData[idCentralNB].fRealY;

	  halfPadLengthY = fPadData[idCentralNB].fHalfPadSize ;

	  if(fabsf(fRecX[nb]) > fabsf(padCenterXB))
	    diffX = fabsf(fRecX[nb]) -  fabsf(padCenterXB);
	  else
	    diffX = fabsf(padCenterXB) -  fabsf(fRecX[nb]);
	  
	  if(fabsf(padCenterYNB)>fabsf(fRecY[b]))
	    diffY = fabsf(padCenterYNB) - fabsf(fRecY[b]);
	  else
	    diffY =  fabsf(fRecY[b]) - fabsf(padCenterYNB);

	  if(diffX < halfPadLengthX && diffY < halfPadLengthY ){//&& fPadData[idCentralB].fIY != 0){
	    
	    if(fNofYNeighbour[b]==2){
	      if(fPadData[idCentralB].fDetElemId<104)
		fRecY[b] += 0.02*sin(14.5*(fRecY[b] - fPadData[idCentralB].fRealY)) ;
	      else if(fPadData[idCentralB].fDetElemId>=200 && fPadData[idCentralB].fDetElemId<204)
		fRecY[b] += 0.02*sin(14.0*(fRecY[b] - fPadData[idCentralB].fRealY)) ;
	      else
		fRecY[b] += 0.025*sin(12.0*(fRecY[b] - fPadData[idCentralB].fRealY)) ;
	    }
	    
	    if(fPadData[idCentralNB].fDetElemId<204)
	      fRecX[nb] += 0.095*sin(10.5*(fRecX[nb] - fPadData[idCentralNB].fRealX)) ;
	    else if(fPadData[idCentralNB].fDetElemId>=300 && fPadData[idCentralNB].fDetElemId<404)
	      fRecX[nb] += 0.085*sin(9.0*(fRecX[nb] - fPadData[idCentralNB].fRealX)) ;
	    else
	      fRecX[nb] += 0.075*sin(9.5*(fRecX[nb] - fPadData[idCentralNB].fRealX)) ;
	    
	    
	    // First check that we have not overflowed the buffer.
	    if((*fRecPointsCount) == fMaxRecPointsCount){
	      HLTWarning("Number of RecHits (i.e. %d) exceeds the max number of RecHit limit %d."
	                " Output buffer is too small.",
	               (*fRecPointsCount),fMaxRecPointsCount
	      );
	      return true;
	    }

	    AliHLTUInt32_t idflags = AliHLTMUONUtils::PackRecHitFlags(
	         (fPadData[idCentralB].fDetElemId / 100) - 1,
	         fPadData[idCentralB].fDetElemId
	      );
	    fRecPoints[(*fRecPointsCount)].fFlags = idflags;
	    fRecPoints[(*fRecPointsCount)].fX = fRecX[nb];
	    fRecPoints[(*fRecPointsCount)].fY = fRecY[b];
	    fRecPoints[(*fRecPointsCount)].fZ = fPadData[idCentralB].fRealZ;
	    
 	    if (fGenerateClusterInfo)
 	    {
 	      if (fClusterCount >= fMaxClusters)
 	      {
	        HLTError("Ran out of space in internal cluster array of size %d.", fMaxClusters);
	        return false;
	      }
	      
	      fClusters[fClusterCount].fId = (fNewClusterId << 5) | fDDL;
	      
	      // Increment the cluster ID and warp it around at 0x03FFFFFF since
	      // the bottom 5 bits are filled with the source DDL number and the
	      // sign bit in fClusters[fClusterCount].fId must be positive.
	      fNewClusterId = (fNewClusterId + 1) & 0x03FFFFFF;
	      
	      fClusters[fClusterCount].fHit = fRecPoints[(*fRecPointsCount)];
	      fClusters[fClusterCount].fDetElemId = fPadData[idCentralB].fDetElemId;
 	      fClusters[fClusterCount].fNchannelsB = fNofBChannel[b];
 	      fClusters[fClusterCount].fNchannelsNB = fNofNBChannel[nb];
 	      fClusters[fClusterCount].fChargeB = fTotChargeY[b];
 	      fClusters[fClusterCount].fChargeNB = fTotChargeX[nb];
 	      fClusterCount++;
 	    }
 	    
 	    if (fGenerateChannelInfo)
 	    {
 	      // 3 by 3 pad structure around the central pad for the 2 planes.
 	      int pad[2][3][3] = {{
 	          {0,      0,      0},
 	          {0, idCentralB,  0},
 	          {0,      0,      0}
 	        },{
 	          {0,      0,      0},
 	          {0, idCentralNB, 0},
 	          {0,      0,      0}
 	      }};
              
              // Find the pad index numbers for the central pads and the pads surrounding them.
              // All these pads would have contributed to the cluster as long as their charge is != 0.
              if (fPadData[idCentralB].fIY > 0)
                pad[0][0][1] = fGetIdTotalData[fPadData[idCentralB].fIX][fPadData[idCentralB].fIY-1][0];
              if (fPadData[idCentralB].fIY < 236)
                pad[0][2][1] = fGetIdTotalData[fPadData[idCentralB].fIX][fPadData[idCentralB].fIY+1][0];
              if (fPadData[idCentralB].fIX > 0)
              {
                int idLeft = fGetIdTotalData[fPadData[idCentralB].fIX-1][fPadData[idCentralB].fIY][0];
                pad[0][1][0] = idLeft;
                if (fPadData[idLeft].fIY > 0)
                  pad[0][0][0] = fGetIdTotalData[fPadData[idLeft].fIX][fPadData[idLeft].fIY-1][0];
                if (fPadData[idLeft].fIY < 236)
                  pad[0][2][0] = fGetIdTotalData[fPadData[idLeft].fIX][fPadData[idLeft].fIY+1][0];
              }
              if (fPadData[idCentralB].fIX < 335)
              {
                int idRight = fGetIdTotalData[fPadData[idCentralB].fIX+1][fPadData[idCentralB].fIY][0];
                pad[0][1][2] = idRight;
                if (fPadData[idRight].fIY > 0)
                  pad[0][0][2] = fGetIdTotalData[fPadData[idRight].fIX][fPadData[idRight].fIY-1][0];
                if (fPadData[idRight].fIY < 236)
                  pad[0][2][2] = fGetIdTotalData[fPadData[idRight].fIX][fPadData[idRight].fIY+1][0];
              }
              
              if (fPadData[idCentralNB].fIX > 0)
                pad[1][1][0] = fGetIdTotalData[fPadData[idCentralNB].fIX-1][fPadData[idCentralNB].fIY][1];
              if (fPadData[idCentralNB].fIX < 335)
                pad[1][1][2] = fGetIdTotalData[fPadData[idCentralNB].fIX+1][fPadData[idCentralNB].fIY][1];
              if (fPadData[idCentralNB].fIY > 0)
              {
                int idLower = fGetIdTotalData[fPadData[idCentralNB].fIX][fPadData[idCentralNB].fIY-1][1];
                pad[1][0][1] = idLower;
                if (fPadData[idLower].fIX > 0)
                  pad[1][0][0] = fGetIdTotalData[fPadData[idLower].fIX-1][fPadData[idLower].fIY][1];
                if (fPadData[idLower].fIX < 335)
                  pad[1][0][2] = fGetIdTotalData[fPadData[idLower].fIX+1][fPadData[idLower].fIY][1];
              }
              if (fPadData[idCentralNB].fIY < 236)
              {
                int idUpper = fGetIdTotalData[fPadData[idCentralNB].fIX][fPadData[idCentralNB].fIY+1][1];
                pad[1][2][1] = idUpper;
                if (fPadData[idUpper].fIX > 0)
                  pad[1][2][0] = fGetIdTotalData[fPadData[idUpper].fIX-1][fPadData[idUpper].fIY][1];
                if (fPadData[idUpper].fIX < 335)
                  pad[1][2][2] = fGetIdTotalData[fPadData[idUpper].fIX+1][fPadData[idUpper].fIY][1];
              }
              
              AliHLTInt32_t clusterId = fGenerateClusterInfo ? fClusters[fClusterCount-1].fId : -1;
              
              // Now generate the pad structures from all the pad indices found above.
              for (int i = 0; i < 2; i++)
              for (int j = 0; j < 3; j++)
              for (int k = 0; k < 3; k++)
              {
                AliHLTMUONPad& p = fPadData[pad[i][j][k]];
                // Skip pads that have zero charge because they would not have
                // contributed to the cluster.
                if (p.fCharge <= 0) continue;
                
                UShort_t manuId; UChar_t channelId; UShort_t adc;
                AliHLTMUONRawDecoder::UnpackADC(p.fRawData, manuId, channelId, adc);
                
                fChannels[fChannelCount].fClusterId = clusterId;
                fChannels[fChannelCount].fBusPatch = p.fBusPatch;
                fChannels[fChannelCount].fManu = manuId;
                fChannels[fChannelCount].fChannelAddress = channelId;
                fChannels[fChannelCount].fSignal = adc;
                fChannels[fChannelCount].fRawDataWord = p.fRawData;
                fChannelCount++;
              }
 	    }

	    HLTDebug("Part 4(Slat) : Reconstructed hit (X,Y,Z) : (%f,%f,%f)",
                     fRecPoints[(*fRecPointsCount)].fX,
                     fRecPoints[(*fRecPointsCount)].fY,
                     fRecPoints[(*fRecPointsCount)].fZ
            );
	    (*fRecPointsCount)++;
	  }//if lies wihtin 5.0 mm
	}// condn over fRecX ! = 0.0
      }// loop over NB side
    }// condn on fRecY[b] !=  0.0
  }// loop over B side;
  
  return true;
}


void AliHLTMUONHitReconstructor::Clear()
{
  // function to clear internal arrays.

  HLTDebug("Clearing fPadData and fNofDataInDetElem buffers.");


  for(int iPad=1;iPad<fDigitPerDDL;iPad++){
    fGetIdTotalData[fPadData[iPad].fIX][fPadData[iPad].fIY][fPadData[iPad].fPlane] = 0;
    fPadData[iPad].fDetElemId = 0;
    fPadData[iPad].fIX = 0 ;
    fPadData[iPad].fIY = 0 ;
    fPadData[iPad].fRealX = 0.0 ;
    fPadData[iPad].fRealY = 0.0 ;
    fPadData[iPad].fRealZ = 0.0 ;
    fPadData[iPad].fHalfPadSize = -1 ;
    fPadData[iPad].fPlane = -1 ;
    fPadData[iPad].fCharge = 0 ;
    fPadData[iPad].fBusPatch = -1;
    fPadData[iPad].fRawData = 0;
  }  

  for( Int_t idet=0;idet<GetkNofDetElemInDDL(fDDL);idet++)
    fNofDataInDetElem[idet] = 0;
  
  // for(int i=0;i<130;i++)
  //   fMaxFiredPerDetElem[i] = 0;
}


AliHLTInt32_t AliHLTMUONHitReconstructor::GetkNofDetElemInDDL(Int_t iDDL)
{
	/// Returns the number of detection elements for a DDL.
	
	if(iDDL>=0 && iDDL<=19)
		return fgkNofDetElemInDDL[iDDL];
	else
		return -1;
}


AliHLTInt32_t AliHLTMUONHitReconstructor::GetkMinDetElemIdInDDL(Int_t iDDL)
{
	/// Returns the first detection element ID for a DDL.
	
	if(iDDL>=0 && iDDL<=19)
		return fgkMinDetElemIdInDDL[iDDL];
	else
		return -1;
}


AliHLTMUONHitReconstructor::AliHLTMUONRawDecoder::AliHLTMUONRawDecoder() :
	fkBufferStart(NULL),
	fBusPatchId(0),
	fDCCut(-1),
	fPadData(NULL),
	fkLookUpTableData(NULL),
	fkIdToEntry(),
	fkMaxEntryPerBusPatch(),
	fDDL(-1),
	fDataCount(1),
	fPrevDetElemId(0),
	fPadCharge(0),
	fCharge(0.0),
	fIdManuChannel(0x0),
	fLutEntry(0),
	fDataCountListPerDetElem(NULL),
	fNofDataInDetElem(NULL),
	fWarnOnly(false),
	fSkipParityErrors(false),
	fDontPrintParityErrors(false),
	fPrintParityErrorAsWarning(false),
	fParityErrorFound(false),
	fNonParityErrorFound(false),
	fIsMuchNoisy(false)
{
	// ctor
}


AliHLTMUONHitReconstructor::AliHLTMUONRawDecoder::~AliHLTMUONRawDecoder()
{
	// dtor
}


void AliHLTMUONHitReconstructor::AliHLTMUONRawDecoder::OnNewBuffer(const void* buffer, UInt_t /*bufferSize*/)
{
	/// Called for every new raw DDL data payload being processed.
	/// Just clears internal counters.
	/// \param buffer  The pointer to the raw data buffer.

	assert( buffer != NULL );
	fkBufferStart = buffer;
	// dataCount starts from 1 because the 0-th element of fPadData is used as null value.
	fDataCount = 1;
	for( Int_t idet=0;idet<GetkNofDetElemInDDL(fDDL);idet++)
	  fNofDataInDetElem[idet] = 0;
	fPrevDetElemId = 0;
	fParityErrorFound = false;
	fNonParityErrorFound = false;
};


void AliHLTMUONHitReconstructor::AliHLTMUONRawDecoder::OnError(ErrorCode code, const void* location)
{
	/// Called if there was an error detected in the raw DDL data.
	/// Logs an error message.
	/// \param code  The error code describing the problem.
	/// \param location  A pointer to the location in the raw data buffer
	///      where the problem was found.
	
	if (code == kParityError)
	{
		fParityErrorFound = true;
	}
	else
	{
		fNonParityErrorFound = true;
	}
	if (fDontPrintParityErrors and code == kParityError) return;
	
	long bytepos = long(location) - long(fkBufferStart) + sizeof(AliRawDataHeader);
	if (fWarnOnly or (fPrintParityErrorAsWarning and code == kParityError))
	{
		HLTWarning("There is a problem with decoding the raw data."
			" %s (Error code: %d, at byte %d). Trying to recover from corrupt data.",
			ErrorCodeToMessage(code), code, bytepos
		);
	}
	else
	{
		HLTError("There is a problem with decoding the raw data. %s (Error code: %d, at byte %d)",
			ErrorCodeToMessage(code), code, bytepos
		);
	}
};

void  AliHLTMUONHitReconstructor::AliHLTMUONRawDecoder::OnNewBusPatch(const AliMUONBusPatchHeaderStruct* header, const void* /*data*/) 
{
  // operation to perform on new data
  fBusPatchId = int(header->fBusPatchId);
  MaxEntryPerBusPatch& maxEntryPerBusPatch 
    = * const_cast<MaxEntryPerBusPatch*>(fkMaxEntryPerBusPatch);
  fIsMuchNoisy = false;
  if(AliHLTInt32_t(header->fLength)> maxEntryPerBusPatch[fBusPatchId])
    fIsMuchNoisy = true;
  
};

void AliHLTMUONHitReconstructor::AliHLTMUONRawDecoder::OnData(UInt_t dataWord, bool parityError)
{
  //function to arrange the decoded Raw Data

  if (fSkipParityErrors and parityError) return;
  
  if(fIsMuchNoisy) return;

  fIdManuChannel = 0x0;
  fIdManuChannel = (fIdManuChannel|fBusPatchId)<<17;
  fIdManuChannel |= (dataWord >> 12) & 0x1FFFF;
  
  IdManuChannelToEntry& idToEntry = * const_cast<IdManuChannelToEntry*>(fkIdToEntry);
  fLutEntry = idToEntry[fIdManuChannel];
  if(fLutEntry==0)
  {
    HLTDebug("Failed to find a valid LUT entry.");
    return;
  }
  fPadCharge = int(((unsigned short)(dataWord & 0xFFF)) - fkLookUpTableData[fLutEntry].fPed);
  
  fCharge = 0;	  
  if(fPadCharge > 2.0*fkLookUpTableData[fLutEntry].fSigma){  // (charge > 4) is due cut out the noise level  			
      
    fPadData[fDataCount].fDetElemId = fkLookUpTableData[fLutEntry].fDetElemId;
    fPadData[fDataCount].fIX = fkLookUpTableData[fLutEntry].fIX;
    fPadData[fDataCount].fIY = fkLookUpTableData[fLutEntry].fIY;
    fPadData[fDataCount].fRealX = fkLookUpTableData[fLutEntry].fRealX;
    fPadData[fDataCount].fRealY = fkLookUpTableData[fLutEntry].fRealY;
    fPadData[fDataCount].fRealZ = fkLookUpTableData[fLutEntry].fRealZ;
    fPadData[fDataCount].fHalfPadSize = fkLookUpTableData[fLutEntry].fHalfPadSize;
    fPadData[fDataCount].fPadSizeXY = fkLookUpTableData[fLutEntry].fPadSizeXY;
    fPadData[fDataCount].fPlane = fkLookUpTableData[fLutEntry].fPlane;
    fPadData[fDataCount].fBusPatch = fBusPatchId;
    fPadData[fDataCount].fRawData = dataWord;
    
    if ( fPadCharge < fkLookUpTableData[fLutEntry].fThres ) {
      fCharge = (fkLookUpTableData[fLutEntry].fA0)*fPadCharge;
    }else{
      fCharge = (fkLookUpTableData[fLutEntry].fA0)*(fkLookUpTableData[fLutEntry].fThres) 
	+ (fkLookUpTableData[fLutEntry].fA0)*(fPadCharge-fkLookUpTableData[fLutEntry].fThres) 
	+ (fkLookUpTableData[fLutEntry].fA1)*(fPadCharge-fkLookUpTableData[fLutEntry].fThres)*(fPadCharge-fkLookUpTableData[fLutEntry].fThres);
    }
    
    fPadData[fDataCount].fCharge = fCharge;
    
    if(fkLookUpTableData[fLutEntry].fDetElemId/100 == 6){
      fDataCountListPerDetElem[
			       ((fkLookUpTableData[fLutEntry].fDetElemId-(GetkMinDetElemIdInDDL(fDDL)+100))%GetkNofDetElemInDDL(fDDL) 
				+ GetkNofDetElemInDDL(fDDL)/2)
			       ]
	[fNofDataInDetElem[((fkLookUpTableData[fLutEntry].fDetElemId-(GetkMinDetElemIdInDDL(fDDL)+100))%GetkNofDetElemInDDL(fDDL) 
			    + GetkNofDetElemInDDL(fDDL)/2)
			   ]++] = fDataCount;
    }else{
      fDataCountListPerDetElem[(fkLookUpTableData[fLutEntry].fDetElemId-GetkMinDetElemIdInDDL(fDDL))%GetkNofDetElemInDDL(fDDL)]
	[fNofDataInDetElem[(fkLookUpTableData[fLutEntry].fDetElemId-GetkMinDetElemIdInDDL(fDDL))%GetkNofDetElemInDDL(fDDL)]++] = fDataCount;
    }
    
    // if(fkLookUpTableData[fLutEntry].fDetElemId != fPrevDetElemId){
    //   if((*fNofFiredDetElem)>0){
    // 	fMaxFiredPerDetElem[(*fNofFiredDetElem)-1] = fDataCount;
    //   }
    
    //   HLTDebug("detElem : %d, prevDetElem : %d, datacount : %d, maxFiredPerDetElem[%d] : %d",
    //            fkLookUpTableData[fLutEntry].fDetElemId,fPrevDetElemId,fDataCount,
    //            ((*fNofFiredDetElem)-1),fMaxFiredPerDetElem[(*fNofFiredDetElem)-1]
    // 	       );
    
    //   (*fNofFiredDetElem)++;
    //   fPrevDetElemId =  fkLookUpTableData[fLutEntry].fDetElemId ;
    // }
    
    //if(fPadData[fDataCount].fDetElemId==102)
    
    HLTDebug("%x, fLutEntry : %d, id : %d, (fDDL,buspatch,detele,plane) : (%2d,%4d,%4d,%1d) (manu,channel) : (%4d,%2d) \n(iX,iY) : (%3d,%3d) (X,Y) : (%f, %f, %f), adc : %d, charge : %f, ped : %f, sigma : %f, padsize : %f",
	       fkLookUpTableData,fLutEntry,fIdManuChannel,fDDL,
	       fBusPatchId,fPadData[fDataCount].fDetElemId,fPadData[fDataCount].fPlane,
    	     ((dataWord >> 18) & 0x7FF),((dataWord >> 12) & 0x3F),
    	     fPadData[fDataCount].fIX,fPadData[fDataCount].fIY,
    	     fPadData[fDataCount].fRealX,fPadData[fDataCount].fRealY,fPadData[fDataCount].fRealZ,
    	     (dataWord & 0xFFF),fPadData[fDataCount].fCharge,fkLookUpTableData[fLutEntry].fPed,fkLookUpTableData[fLutEntry].fSigma,
    	     fkLookUpTableData[fLutEntry].fHalfPadSize);
    
    // HLTDebug("%x, fLutEntry : %d, buspatch : %d, detele : %d, id : %d, manu : %d, channel : %d, iX : %d, iY: %d, (X,Y) : (%f, %f, %f), charge : %f, padsize : %f, plane : %d, ped : %f, sigma : %f",
    // 	     fkLookUpTableData,fLutEntry,fBusPatchId,fPadData[fDataCount].fDetElemId,
    // 	     fIdManuChannel,((dataWord >> 18) & 0x7FF),((dataWord >> 12) & 0x3F),
    // 	     fPadData[fDataCount].fIX,fPadData[fDataCount].fIY,
    // 	     fPadData[fDataCount].fRealX,fPadData[fDataCount].fRealY,fPadData[fDataCount].fRealZ,
    // 	     fPadData[fDataCount].fCharge,fkLookUpTableData[fLutEntry].fHalfPadSize,fPadData[fDataCount].fPlane,fkLookUpTableData[fLutEntry].fPed,fkLookUpTableData[fLutEntry].fSigma);
    
    fDataCount++;
  }// if charge is more than DC Cut limit condition
  
}
