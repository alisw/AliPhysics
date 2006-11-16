// $Id$

/** \class AliHLTTPCDisplayMain
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayMain
//
// Display class for the HLT events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#define DEBUG 0


#if defined(HAVE_HOMERREADER) 
#include "HOMERReader.h"
#endif // defined(HAVE_HOMERREADER) 
//-----------
#include "AliHLTTPCDisplayMain.h"
#include "AliHLTTPCDisplayCharge.h"
#include "AliHLTTPCDisplayPadRow.h"
#include "AliHLTTPCDisplayPad.h"
#include "AliHLTTPCDisplay3D.h"
#include "AliHLTTPCDisplayResiduals.h"
#include "AliHLTTPCDisplayFront.h"
//-----------
#include "AliHLTStdIncludes.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCDigitReaderRaw.h"
#include "AliHLTTPCPad.h"
#include "AliHLT_C_Component_WrapperInterface.h"
#include "AliHLTTPCLogging.h"

#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDigitReaderPacked.h"

#ifdef use_aliroot
#include <TClonesArray.h>
#include <AliRun.h>
#include <AliSimDigits.h>
#include <AliTPCParam.h>
#endif

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCDisplayMain)

//____________________________________________________________________________________________________
AliHLTTPCDisplayMain::AliHLTTPCDisplayMain(void* pt2GUI, void (*pt2Function)(void*, Int_t)) {
  //constructor

  // set N of TimeBins
  fgNTimeBins = 1024;   //  446 or 1024

  //callback handler
  fPt2Gui = pt2GUI;
  fPadCallback = pt2Function;
	
  fReader = NULL;
  
  fConnect = kFALSE;
  fEventID = 0;


  fDisplayCharge = NULL; 
  fDisplayPadRow = NULL;
  fDisplayPad = NULL;
  fDisplay3D = NULL;
  fDisplayResiduals = NULL;
  fDisplayFront = NULL;

  fCanvasCharge = NULL;
  fCanvasPadRow = NULL;
  fCanvasPad = NULL;
  fCanvas3D = NULL;
  fCanvasResiduals = NULL;
  fCanvasFront = NULL;
  fCanvasHits_S = NULL;
  fCanvasQ_Track = NULL;
  fCanvasQ_S = NULL;
  fCanvasPadRow_Pad = NULL;


  fExistsRawData = kFALSE;
  fExistsClusterData = kFALSE;
  fExistsTrackData = kFALSE;

  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
  memset(fNcl, 0, 36*6*sizeof(UInt_t)); 

  fTracks = NULL;

  fZeroSuppression = kTRUE;
  fNPads = AliHLTTPCTransform::GetNPads(0);
  fPad = -1;
  fPadRow = 0;
  fSlicePadRow = 0; 
  fSplitPadRow = kFALSE;
  
  fFrontDataSwitch = 2;
  fTimeBinMin = 0;
  fTimeBinMax = GetNTimeBins() -1;
  fSplitFront = kFALSE;
    
  fSelectTrack = -1;
  fSelectTrackSlice = 0;
  fSelectTrackSwitch = kFALSE;
  
  fSelectCluster = 0;

  fMinSlice = 0;
  fMaxSlice = 35;
  fSlicePair = kFALSE;

  SetSliceArray();
  
  fTheta = 90.;
  fPhi = 0.;

  fBackColor = 1;
  fLineColor = 0;
  fKeepView = kTRUE;

  fSwitch3DCluster = kFALSE;
  fSwitch3DTracks = kFALSE;
  fSwitch3DPadRow = kFALSE;
  fSwitch3DGeometry = kTRUE;
  fSwitch3DRaw = 0;

  fCutHits = 0;
  fCutPt = 0.;
  fCutPsi = 0.;
  fCutLambda = 0.;
  fCut_S = 0.;
  fCutPadrow = 159;
  

  fTrackParam.kappa = 0.;
  fTrackParam.nHits = 0;
  fTrackParam.charge = 0;
  fTrackParam.radius = 0.;
  fTrackParam.slice = 0;
  fTrackParam.phi0 = 0.;
  fTrackParam.psi = 0.;
  fTrackParam.lambda = 0.;
  fTrackParam.pt = 0.;
  fTrackParam.id = 0;
  fTrackParam.bfield = 0.;
  fTrackParam.s = 0.;
}

//____________________________________________________________________________________________________
AliHLTTPCDisplayMain::~AliHLTTPCDisplayMain() {
    //destructor
    if ( fTracks )  delete fTracks;
    fTracks = NULL;

#if defined(HAVE_HOMERREADER) 
    HOMERReader* reader= (HOMERReader*) fReader;

    if (reader)
      delete reader;
    fReader = NULL;
#endif // defined(HAVE_HOMERREADER) 

    if (fDisplayResiduals)
      delete fDisplayResiduals;
    fDisplayResiduals = NULL;

    if (fDisplay3D)
      delete fDisplay3D;
    fDisplay3D = NULL;

    if (fDisplayPad)
      delete fDisplayPad;
    fDisplayPad = NULL;

    if (fDisplayPadRow)
      delete fDisplayPadRow;
    fDisplayPadRow = NULL;

    if (fDisplayFront)
      delete fDisplayFront;
    fDisplayFront = NULL;

    if (fDisplayCharge)
      delete fDisplayCharge;
    fDisplayCharge = NULL;
}

//____________________________________________________________________________________________________
Int_t AliHLTTPCDisplayMain::Connect( unsigned int cnt, const char** hostnames, unsigned short* ports, Char_t *gfile){

#if defined(HAVE_HOMERREADER) 
  // -- input datatypes , reverse
  Char_t* spptID="SRETSULC";       // CLUSTERS
  Char_t* trkID = "SGESKART";      // TRAKSEGS
  Char_t* padrowID = "KPWR_LDD";   // DDL_RWPK
  
  Int_t ret;
  
  while(1) {

    // -- CONNECT to TCPDumpSubscriber via HOMER
    // ---------------------------------------------
    HOMERReader* reader = new HOMERReader( cnt, (const char**) hostnames, ports );
    ret=reader->GetConnectionStatus();
    
    if ( ret ) {	
      Int_t ndx = reader->GetErrorConnectionNdx();
      if ( ndx < (Int_t) cnt) HLTError("Error establishing connection to TCP source %s:%hu: %s (%d)\n", hostnames[ndx], ports[ndx], strerror(ret), ret );
      else HLTError("Error establishing connection to unknown source with index %d: %s (%d)\n",ndx, strerror(ret), ret );
      
      delete reader;
      reader = NULL;
      return -1;
    }

    // -- ERROR HANDLING for HOMER (error codes and empty blocks)
    // --------------------------------------------------------------
    Int_t  ret1 = reader->ReadNextEvent(3000000); // timeout in microseconds (3s)
    if ( ret1 ) {	  
      Int_t ndx = reader->GetErrorConnectionNdx();
      
      if ( ret1 == 111 || ret1 == 32 || ret1 == 6 ) {
	HLTError( "Error, No Connection to source %d: %s (%d)\n", ndx, strerror(ret1), ret1 );
	
	delete reader;
	reader = NULL;
	
	return -1; 
      }
      else if (ret1 == 110){
	HLTError( "Timout occured, reading event from source %d: %s (%d)\n", ndx, strerror(ret1), ret1 );
	
	delete reader;
	reader = NULL;
	
	return -2;
      }
      else if ( ret1 == 56) {
	HLTError( "Error reading event from source %d: %s (%d) -- RETRY\n", ndx, strerror(ret1), ret1 );
	
	delete reader;
	reader = NULL;
	
	continue; 
      }
      else {
	HLTError( "General Error reading event from source %d: %s (%d)\n", ndx, strerror(ret1), ret1 );
	
	delete reader;
	reader = NULL;
	
	return -3;
      }
    }

    unsigned long blockCnt = reader->GetBlockCnt();  

    HLTDebug( "Event 0x%016LX (%Lu) with %lu blocks\n", (ULong64_t)reader->GetEventID(), (ULong64_t)reader->GetEventID(), blockCnt );

    if (blockCnt == 0){
      delete reader;
      reader = NULL;
      continue;
    }
    
    // SWITCH on the Display parts according what DATATAYPES are present
    // ---------------------------------------------------------------------
    for ( unsigned long i = 0; i < blockCnt; i++ ) {
      char tmp1[9];
      memset( tmp1, 0, 9 );
      void *tmp11 = tmp1;
      ULong64_t* tmp12 = (ULong64_t*)tmp11;
      *tmp12 = reader->GetBlockDataType( i );
      cout << "XXX=" << tmp1 << endl;
      if (!strcmp(tmp1,padrowID)) fExistsRawData = kTRUE;
      else if (!strcmp(tmp1,spptID)) fExistsClusterData = kTRUE;
      else if (!strcmp(tmp1,trkID)) fExistsTrackData = kTRUE;   	  
    }
  
    // -- Set reader and eventID
    // -------------------------------
    fEventID = (ULong64_t)reader->GetEventID();
    fReader = (void*) reader;
    fConnect = kTRUE; 

    break; // leave while(1), if connected to source
  }

  
  // -- Initialize TPC Display Classes  -- IMPORTANT... don't change the order of them
  //    Check if necessary data types are present
  // ------------------------------------------------------------------------------------
  fDisplay3D = new AliHLTTPCDisplay3D(this, gfile);
  
  if (ExistsRawData()){
    fDisplayPadRow = new AliHLTTPCDisplayPadRow(this);
    fDisplayPad = new AliHLTTPCDisplayPad(this);
    fDisplayFront = new AliHLTTPCDisplayFront(this);
  }
  
  if (ExistsClusterData()){
    fDisplayCharge = new AliHLTTPCDisplayCharge(this);
  }
  
  if (ExistsTrackData() && ExistsClusterData() ){
    fDisplayResiduals = new AliHLTTPCDisplayResiduals(this);
  }

  return 0;

#else
    HLTFatal("HOMER reader not available");
#endif // defined(HAVE_HOMERREADER) 
}

//____________________________________________________________________________________________________
Int_t AliHLTTPCDisplayMain::Disconnect(){
#if defined(HAVE_HOMERREADER) 
    // READ CLUSTER and TRACK data
    HOMERReader* reader = (HOMERReader*)fReader;

    if (reader)
      delete reader;
    fReader = NULL;

#else
    HLTFatal("HOMER raeder not available");
#endif // defined(HAVE_HOMERREADER) 

    return 0;
}


//____________________________________________________________________________________________________
Int_t AliHLTTPCDisplayMain::ReadData(Bool_t nextSwitch){
// READ CLUSTER and TRACK data

#if defined(HAVE_HOMERREADER) 
  HOMERReader* reader = (HOMERReader*)fReader;
  Char_t* spptID="SRETSULC";       // CLUSTERS
  Char_t* trkID = "SGESKART";      // TRAKSEGS

  // -- reset TracksPerSlice Array
  for(Int_t ii=0;ii<36;ii++) fTracksPerSlice[ii] = 0;

  Int_t ret;
  ULong_t blk;

  // -- READ next event data and ERROR HANDLING for HOMER (error codes and empty blocks)
  // ---------------------------------------------------------------------------------------
  if (nextSwitch) {
    
    while(1){
      ret = reader->ReadNextEvent();
      
      if ( ret == 111 || ret == 32 || ret == 6 ) {
	Int_t ndx = reader->GetErrorConnectionNdx();
	HLTError( "Error, No Connection to source %d: %s (%d)\n", ndx, strerror(ret), ret );
	return -1; 
      }
      else if ( ret == 110 ) {
	Int_t ndx = reader->GetErrorConnectionNdx();
	HLTError( "Timout occured, reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
	return -2; 
      }
      else if ( ret == 56) {
	Int_t ndx = reader->GetErrorConnectionNdx();
	HLTError( "Error reading event from source %d: %s (%d) -- RETRY\n", ndx, strerror(ret), ret );
	continue; 
      }
      else if ( ret ) {
	Int_t ndx = reader->GetErrorConnectionNdx();
	HLTError( "General Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
	return -3;
      }
      else {
	break;
      }
    }
  } // end if (nextSwitch)

  // -- Set blockCnt and eventID
  // -------------------------------
  ULong_t blockCnt = reader->GetBlockCnt();
  fEventID = (ULong64_t)reader->GetEventID();

  HLTDebug( "Event 0x%016LX (%Lu) with %lu blocks\n", (ULong64_t)reader->GetEventID(), (ULong64_t)reader->GetEventID(), blockCnt );

  // Loop for Debug only
  for ( ULong_t i = 0; i < blockCnt; i++ ) {
    Char_t tmp1[9], tmp2[5];
    memset( tmp1, 0, 9 );
    memset( tmp2, 0, 5 );
    void *tmp11 = tmp1;
    ULong64_t* tmp12 = (ULong64_t*)tmp11;
    *tmp12 = reader->GetBlockDataType( i );
    void *tmp21 = tmp2;
    ULong_t* tmp22 = (ULong_t*)tmp21;
    *tmp22 = reader->GetBlockDataOrigin( i );

    HLTDebug( "Block %lu length: %lu - type: %s - origin: %s\n",i, reader->GetBlockDataLength( i ), tmp1, tmp2 );
  } // end for ( ULong_t i = 0; i < blockCnt; i++ ) {


  // -- READ RAW DATA
  //---------------------
  ReadRawData();
  
  // -- READ CLUSTER DATA
  //-------------------------
  ReadClusterData();

  // -- READ TRACKS
  //-------------------
  ReadTrackData();

#else
  HLTFatal("HOMER raeder not available");
#endif // defined(HAVE_HOMERREADER) 

  return 0;
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::DisplayEvent(Bool_t newRawSlice){
#if defined(HAVE_HOMERREADER) 
    HOMERReader* reader = (HOMERReader*)fReader;

    //--------------------------------------------------------------------------------------------
    // READ RAW DATA for PADROW/PAD HISTOGRAM
    //--------------------------------------------------------------------------------------------
    if ( ExistsRawData() ){	
	fDisplayPadRow->Reset();
	fDisplayPad->Reset();
	fDisplayFront->Reset();

	char* rawID = "KPWR_LDD";
	Int_t padrowpatch =  AliHLTTPCTransform::GetPatch(fPadRow);
	Int_t maxpadrowpatch =  padrowpatch;


	if (newRawSlice) {
	  printf("=!=!=!=!=!= READ NEW ARRAY!\n");
	  ReadRawData();
	}
	
	if (fPadRow == 30 || fPadRow == 90 || fPadRow == 139) maxpadrowpatch++;

	for (Int_t patch=0; patch < 6; patch++){
	    unsigned long blk;
	    AliHLTUInt32_t rawSpec = AliHLTTPCDefinitions::EncodeDataSpecification( (AliHLTUInt8_t) fSlicePadRow, (AliHLTUInt8_t) fSlicePadRow,
										    (AliHLTUInt8_t) patch,(AliHLTUInt8_t) patch );
	    // READ RAW DATA BLOCK
	    blk = reader->FindBlockNdx( rawID, " CPT", rawSpec );	

	    if ( ~(unsigned long)0 != blk ){

#if HOMER_VERSION >= 2
	      // Check for corrupt data
	      AliHLTUInt64_t corruptFlag = reader->GetBlockStatusFlags( blk );
	      if (corruptFlag & 0x00000001) {
		LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplayMain::ReadData","Block status flags") << "Data block is corrupt"<<ENDLOG; 
		continue;
	      }
#endif

#if DEBUG
		printf( "Raw Data found for slice %u/patch %u\n", fSlicePadRow, patch );
#endif
		unsigned long rawDataBlock = (unsigned long) reader->GetBlockData( blk );
		unsigned long rawDataLen = reader->GetBlockDataLength( blk );
	    
		fDisplayFront->Fill(patch, rawDataBlock, rawDataLen);

		if (patch == padrowpatch || patch == maxpadrowpatch ) {
		    fDisplayPadRow->Fill(patch, rawDataBlock, rawDataLen);
		    fDisplayPad->Fill(patch, rawDataBlock, rawDataLen);
		}
	    }
	}

#if 0
	// -- READ raw data of one sector into 3D arrays fRawData and fRawDataZeroSuppressed
	// -------------------------------------------------------------------------------------
	
	// -- reset arrays
	memset(fRawData, 0, 159*140*1024*sizeof(UInt_t));
	memset(fRawDataZeroSuppressed, 0, 159*140*1024*sizeof(UInt_t));

	// -- read raw data blocks
	// -----------------------
	ULong_t blk;
	blk = reader->FindBlockNdx( rawID, " CPT",0xFFFFFFFF );

	while ( blk != ~(ULong_t)0 ) {
	  
	  HLTDebug( "Found raw data block %lu\n", blk );

	  // -- Check for corrupt data
	  AliHLTUInt64_t corruptFlag = reader->GetBlockStatusFlags( blk );
	  if (corruptFlag & 0x00000001) {
	    HLTError( "Data block %lu is corrupt\n",blk );
	    blk = reader->FindBlockNdx( rawID, " CPT", 0xFFFFFFFF, blk+1 );
	    continue;
	  }

	  void* rawDataBlock = (void*) reader->GetBlockData( blk );
	  unsigned long rawDataLen = reader->GetBlockDataLength( blk );

	  ULong_t spec = reader->GetBlockDataSpec( blk );
	  Int_t patch = AliHLTTPCDefinitions::GetMinPatchNr( spec );
	  Int_t slice = AliHLTTPCDefinitions::GetMinSliceNr( spec );

	  HLTDebug( "Raw data found for slice %u - patch %u\n", slice, patch );

	  // -- Wrong slice for raw data 
	  if ( GetSlicePadRow() != slice) {
	    blk = reader->FindBlockNdx( rawID, " CPT", 0xFFFFFFFF, blk+1 );
	    continue;
	  } 

#if defined(HAVE_TPC_MAPPING)

	  // -- Read data out of block
	  AliHLTTPCDigitReaderRaw digitReader(0);

    	  // Initialize RAW DATA
	  Int_t firstRow = AliHLTTPCTransform::GetFirstRow(patch);
	  Int_t lastRow = AliHLTTPCTransform::GetLastRow(patch);
	  
	  // Outer sector, patches 2, 3, 4, 5 -  start counting in patch 2 with row 0
	  Int_t rowOffset = 0;
	  if ( patch >= 2 ) rowOffset = AliHLTTPCTransform::GetFirstRow( 2 );

	  // Initialize block for reading packed data
	  digitReader.InitBlock(rawDataBlock,rawDataLen,firstRow,lastRow,patch,slice);

	  Bool_t readValue = digitReader.Next();

	  if (!readValue){	
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplayPadRow::Fill","Read first value") << "No value in data block" << ENDLOG;
	    blk = reader->FindBlockNdx( rawID, " CPT", 0xFFFFFFFF, blk+1 );
	    continue;
	  }

	  // LOOP over raw data and fill arrays
	  while ( readValue ){ 
	    
	    Int_t row = digitReader.GetRow() + rowOffset;
	    Int_t pad = (Int_t) digitReader.GetPad();
	    Int_t time = (Int_t) digitReader.GetTime();
	    UInt_t signal = (UInt_t) digitReader.GetSignal();

	    fRawData[row][pad][time] = signal;
	    fRawDataZeroSuppressed[row][pad][time] = signal;
	    
	    // -- read next value
	    readValue = digitReader.Next();
      
	    // -- No more value
	    if(!readValue) break; 
	  } // end  while ( readValue ){ 


#else //! defined(HAVE_TPC_MAPPING)
	  HLTFatal("DigitReaderRaw not available - check your build");
#endif //defined(HAVE_TPC_MAPPING)
	  
	  blk = reader->FindBlockNdx( rawID, " CPT", 0xFFFFFFFF, blk+1 );

	} // end while ( blk != ~(ULong_t)0 ) {
#endif

	fDisplayPad->Draw();
	fDisplayPadRow->Draw();
	fDisplayFront->Draw();
    }

    //--------------------------------------------------------------------------------------------
    // DRAW RESIDUALS
    //--------------------------------------------------------------------------------------------
    if ( ExistsClusterData() && ExistsTrackData() ){
	fDisplayResiduals->Reset();
	fDisplayResiduals->Fill();
	fDisplayResiduals->Draw();
    }

    //--------------------------------------------------------------------------------------------
    // DRAW CHARGE 
    //--------------------------------------------------------------------------------------------
    if ( ExistsClusterData() ){
	fDisplayCharge->Reset();
	fDisplayCharge->Fill();
	fDisplayCharge->Draw();
    }

    //--------------------------------------------------------------------------------------------
    // DRAW 3D
    //--------------------------------------------------------------------------------------------

    // TAKE CARE !!! EXISTSxxxData() HANDLING of 3D will be done IN this class !!!
    fDisplay3D->Draw(); 
#else
    HLTFatal("HOMER raeder not available");
#endif // defined(HAVE_HOMERREADER) 

}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::SaveHistograms(){    
  // Save histograms as eps
  fDisplayCharge->Save();
  fDisplayPadRow->Save();
  fDisplayPad->Save();
  fDisplay3D->Save();
  fDisplayResiduals->Save();
  fDisplayFront->Save();
  fDisplayCharge->Save();
}

//----------------------------------------------------------------------------------------------------
//                 SETUP
//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::SetupCluster(Int_t slice, Int_t patch, UInt_t nofClusters, AliHLTTPCSpacePointData* data)  {  

    if (data && slice>=0 && slice<36 && patch>=0 && patch<AliHLTTPCTransform::GetNPatches()) {
	if (fClusters[slice][patch]!=NULL) {
	    delete(fClusters[slice][patch]);
	    fClusters[slice][patch]=NULL;
	}
	Int_t arraysize=nofClusters*sizeof(AliHLTTPCSpacePointData);
	fClusters[slice][patch] = (AliHLTTPCSpacePointData*)new Byte_t[arraysize];
	if (fClusters[slice][patch]) {
	    memcpy(fClusters[slice][patch], data, arraysize);
	    fNcl[slice][patch]=nofClusters;
	} else {
	    fNcl[slice][patch]=nofClusters;
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupCluster","memory allocation") << "memory allocation failed "<<ENDLOG; 
	}
    } else LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupCluster","argument check") << "invalid argument "<<ENDLOG; 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::SetupTracks() {

    // Set USED cluster
    Int_t ntracks = fTracks->GetNTracks();

    for(Int_t j=0; j<ntracks; j++) { 	
	AliHLTTPCTrack *gtrack = fTracks->GetCheckedTrack(j); 
	if(!gtrack) continue;

	Int_t nHits = gtrack->GetNHits();
	UInt_t *hitnum = gtrack->GetHitNumbers();

	for(Int_t h=0; h<nHits; h++){
	  
	    UInt_t id=hitnum[h];
	    Int_t slice = (id>>25) & 0x7f;
	    Int_t patch = (id>>22) & 0x7;
	    UInt_t pos = id&0x3fffff;	    
		
	    AliHLTTPCSpacePointData *points = fClusters[slice][patch];
	    
	    if(!points) {
		LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplayMain::SetupTracks","Clusterarray") 
		  <<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
		continue;
	    }
 
	    if(pos>=fNcl[slice][patch]) {
		LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplayMain::SetupTracks","Clusterarray") 
		  <<"Pos is too large: pos "<<pos <<" ncl "<<fNcl[slice][patch]<<ENDLOG;
		continue;
	    }
	    points[pos].fUsed = kTRUE;
	    points[pos].fTrackN = j;
	}
    }
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::SetSliceArray() {
    Int_t slice=0;
    Int_t minSlice = fMinSlice; 
    Int_t maxSlice = fMaxSlice; 
    Int_t realslice = 0;

    for (slice=0;slice< 36;slice++){
	fSliceArray[slice] = kFALSE;
    }

    // Single Slice, or Range
    if (minSlice > maxSlice) maxSlice += 17;
	
    for (slice=minSlice;slice<=maxSlice;slice++){
	realslice = slice % 18;
	fSliceArray[realslice] = kTRUE;
	fSliceArray[realslice+18] = kTRUE;
    }

    // Pair of Slices
    if (fSlicePair) {
	minSlice = fMinSlice + 9;
	maxSlice = fMaxSlice + 9;
	
	if (minSlice > maxSlice) maxSlice += 17;
	
	for (slice=minSlice;slice<=maxSlice;slice++){
	    realslice = slice % 18;
	    fSliceArray[realslice] = kTRUE;	
	    fSliceArray[realslice+18] = kTRUE;
	}
    }
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::ReadRawData(){
  // -- READ raw data of one sector into 3D arrays fRawData and fRawDataZeroSuppressed

  // -- reset arrays
  memset(fRawData, 0, 159*140*1024*sizeof(UInt_t));
  memset(fRawDataZeroSuppressed, 0, 159*140*1024*sizeof(UInt_t));

#if defined(HAVE_HOMERREADER) 
  HOMERReader* reader = (HOMERReader*)fReader;
  
  // -- read raw data blocks
  // ---------------------------
  ULong_t blk;
  Char_t* rawID = "KPWR_LDD";    
  blk = reader->FindBlockNdx( rawID, " CPT",0xFFFFFFFF );
  
  while ( blk != ~(ULong_t)0 ) {
    HLTDebug( "Found raw data block %lu\n", blk );

#if HOMER_VERSION >= 2
    // -- Check for corrupt data
    AliHLTUInt64_t corruptFlag = reader->GetBlockStatusFlags( blk );
    if (corruptFlag & 0x00000001) {
      HLTError( "Data block %lu is corrupt\n",blk );
      blk = reader->FindBlockNdx( rawID, " CPT", 0xFFFFFFFF, blk+1 );
      continue;
    }
#endif

    void* rawDataBlock = (void*) reader->GetBlockData( blk );
    unsigned long rawDataLen = reader->GetBlockDataLength( blk );
    
    ULong_t spec = reader->GetBlockDataSpec( blk );
    Int_t patch = AliHLTTPCDefinitions::GetMinPatchNr( spec );
    Int_t slice = AliHLTTPCDefinitions::GetMinSliceNr( spec );
    
    HLTDebug( "Raw data found for slice %u - patch %u\n", slice, patch );

    // -- Wrong slice for raw data 
    if ( GetSlicePadRow() != slice) {
      blk = reader->FindBlockNdx( rawID, " CPT", 0xFFFFFFFF, blk+1 );
      continue;
    } 
    
#if defined(HAVE_TPC_MAPPING)

    // -- Read data out of block
    AliHLTTPCDigitReaderRaw digitReader(0);
    
    // Initialize RAW DATA
    Int_t firstRow = AliHLTTPCTransform::GetFirstRow(patch);
    Int_t lastRow = AliHLTTPCTransform::GetLastRow(patch);
    
    // Outer sector, patches 2, 3, 4, 5 -  start counting in patch 2 with row 0
    Int_t rowOffset = 0;
    if ( patch >= 2 ) rowOffset = AliHLTTPCTransform::GetFirstRow( 2 );
    
    // Initialize block for reading packed data
    digitReader.InitBlock(rawDataBlock,rawDataLen,firstRow,lastRow,patch,slice);
    
    Bool_t readValue = digitReader.Next();

    if (!readValue){	
      HLTError ( "No value in data block %lu \n", blk);
      blk = reader->FindBlockNdx( rawID, " CPT", 0xFFFFFFFF, blk+1 );
      continue;
    }

    // -- Fill Zero Suppressed Array
    // ---------------------------------

    Int_t newRow = 0;    
    UShort_t time=0,newTime=0;
    UInt_t pad=0,newPad=0;
    Int_t signal=0;

    Int_t gatingGridOffset=50;
    Int_t signalThreshold = 10;
    AliHLTTPCPad baseline(gatingGridOffset, GetNTimeBins() );

    // just to make later conversion to a list of objects easier
    AliHLTTPCPad* pCurrentPad=NULL;
    if ( signalThreshold >= 0) {
      pCurrentPad=&baseline;
      baseline.SetThreshold(signalThreshold);
    }

    while ( readValue ) { // Reads through all digits in block
     
      while(1){ //Loop over time bins of current pad
	// read all the values for one pad at once to calculate the base line
	if (pCurrentPad) {
	  
	  if (!pCurrentPad->IsStarted()) {
	    
	    HLTDebug("reading data for pad %d, padrow %d", digitReader.GetPad(), digitReader.GetRow()+rowOffset);
	    
	    pCurrentPad->SetID(digitReader.GetRow()+rowOffset,digitReader.GetPad());
	    
	    if ( (pCurrentPad->StartEvent()) >= 0) {
	      // loop over data of one pad and read it into AliHLTTPCPad class
	      do {
		if ( (digitReader.GetRow()+rowOffset) != pCurrentPad->GetRowNumber() ) break;
		if ( digitReader.GetPad() != pCurrentPad->GetPadNumber() ) break;
		pCurrentPad->SetRawData( digitReader.GetTime(), digitReader.GetSignal() );
		
		HLTDebug("set raw data to pad: bin %d charge %d", digitReader.GetTime(), digitReader.GetSignal());
		
	      } while ( (readValue = digitReader.Next()) != 0 );
	    }

	    // calculate baseline of pad
	    pCurrentPad->CalculateBaseLine( GetNTimeBins() / 2);
	    
	    if ( pCurrentPad->Next(kTRUE/*do zero suppression*/) == 0 ) {
	      HLTDebug("no data available after zero suppression");
	      
	      pCurrentPad->StopEvent();
	      pCurrentPad->ResetHistory();
	      break;
	    }

	    Int_t time = pCurrentPad->GetCurrentPosition();
	    if ( time > pCurrentPad->GetSize() ) {
	      HLTError("invalid time bin for pad");
	      break;
	    }	  
	  } // end - if (!pCurrentPad->IsStarted()) {
	} // end - if (pCurrentPad) {
	
	if (pCurrentPad) {
	  Float_t occupancy=pCurrentPad->GetOccupancy();
	  if ( occupancy < 0.01 ) {
	    signal = pCurrentPad->GetCorrectedData();
	  } 
	  else {
	    signal = 0;
	    HLTInfo("ignoring pad %d with occupancy level %f %%", pCurrentPad->GetPadNumber(), occupancy);
	  }
	  //	  signal = pCurrentPad->GetCorrectedData();
	} else {
	  signal = digitReader.GetSignal();
	}

	fRawDataZeroSuppressed[ pCurrentPad->GetRowNumber() ][ pCurrentPad->GetPadNumber() ][ pCurrentPad->GetCurrentPosition() ] = (UInt_t)signal;

	
	HLTDebug("get next charge value: position %d charge %d", time, signal);
	
	if( time >= AliHLTTPCTransform::GetNTimeBins() ){
	  HLTWarning("Timebin (%d) out of range (%d)", time, AliHLTTPCTransform::GetNTimeBins());
	  break;
	}
 	  
	if (pCurrentPad) {
	
	  if( (pCurrentPad->Next(kTRUE/*do zero suppression*/)) == 0 ) {
	    pCurrentPad->StopEvent();
	    pCurrentPad->ResetHistory();

	    if(readValue) {
	      newPad = digitReader.GetPad();
	      newTime = digitReader.GetTime();
	      newRow = digitReader.GetRow() + rowOffset;
	    }
	    break;
	  }
	  
	  newPad=pCurrentPad->GetPadNumber();
	  newTime=pCurrentPad->GetCurrentPosition();
	  newRow=pCurrentPad->GetRowNumber();
	} 
	else {
	  
	  readValue = digitReader.Next();
	  
	  //Check where to stop:
	  if(!readValue) break; //No more value
	  
	  newPad = digitReader.GetPad();
	  newTime = digitReader.GetTime();
	  newRow = digitReader.GetRow() + rowOffset;
	}
	
	if(newPad != pad) break; //new pad
	if(newTime != time+1) break; //end of sequence
	
	time = newTime;
	
      } // end - while(1){ // Loop over time bins of current pad

      if ( !readValue ) break;

    } // end  while ( readValue ) {

    // -- END ZERO SUPPRESSION

    // Rewind block for reading packed data
    digitReader.InitBlock(rawDataBlock,rawDataLen,firstRow,lastRow,patch,slice);

    // -- Fill Raw Array
    // ---------------------
    readValue = digitReader.Next();
    
    // LOOP over raw data and fill arrays
    while ( readValue ){ 
	    
      Int_t row = digitReader.GetRow() + rowOffset;
      Int_t pad = (Int_t) digitReader.GetPad();
      Int_t time = (Int_t) digitReader.GetTime();
      UInt_t signal = (UInt_t) digitReader.GetSignal();

      fRawData[row][pad][time] = signal;
	    
      // -- read next value
      readValue = digitReader.Next();
      
      // -- No more value
      if(!readValue) break; 

    } // end  while ( readValue ){ 

#else //! defined(HAVE_TPC_MAPPING)
    HLTFatal("DigitReaderRaw not available - check your build");
#endif //defined(HAVE_TPC_MAPPING)
	  
    blk = reader->FindBlockNdx( rawID, " CPT", 0xFFFFFFFF, blk+1 );
    
  } // end while ( blk != ~(ULong_t)0 ) {
  
#else
  HLTFatal("HOMER raeder not available");
#endif // defined(HAVE_HOMERREADER) 
  
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::ReadClusterData(){
 // -- READ cluster data 
#if defined(HAVE_HOMERREADER) 
  HOMERReader* reader = (HOMERReader*)fReader;

  ULong_t blk;
  Char_t* spptID="SRETSULC";       // CLUSTERS
  blk = reader->FindBlockNdx( spptID, " CPT",0xFFFFFFFF );
   
  while ( blk != ~(ULong_t)0 ) {	    
    HLTDebug( "Found clusters block %lu\n", blk );
    
    const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*)reader->GetBlockData( blk );
    if ( !clusterData ) {
      HLTError( "No track data for block %lu\n", blk );
      blk = reader->FindBlockNdx( spptID, " CPT", 0xFFFFFFFF, blk+1 );
      continue;
    }
	
    ULong_t spec = reader->GetBlockDataSpec( blk );
    Int_t patch = AliHLTTPCDefinitions::GetMinPatchNr( spec );
    Int_t slice = AliHLTTPCDefinitions::GetMinSliceNr( spec );

    HLTDebug( "%lu Clusters found for slice %u - patch %u\n", clusterData->fSpacePointCnt, slice, patch );

    void* tmp30 = (void*)clusterData;
    Byte_t* tmp31 = (Byte_t*)tmp30;
    ULong_t offset;
    offset = sizeof(clusterData->fSpacePointCnt);
    if ( offset <= reader->GetBlockTypeAlignment( blk, 1 ) ) offset = reader->GetBlockTypeAlignment( blk, 1 );
    tmp31 += offset;
    tmp30 = tmp31;
    AliHLTTPCSpacePointData* tmp32 = (AliHLTTPCSpacePointData*)tmp30;
	
    SetupCluster( slice, patch, clusterData->fSpacePointCnt, tmp32 );
    
    blk = reader->FindBlockNdx( spptID, " CPT", 0xFFFFFFFF, blk+1 );
  }
  
#else
  HLTFatal("HOMER raeder not available");
#endif // defined(HAVE_HOMERREADER) 
}
//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::ReadTrackData(){
  // -- READ track data
  
  if ( fTracks ) delete fTracks; 
  fTracks = new AliHLTTPCTrackArray;
  
#if defined(HAVE_HOMERREADER) 
  HOMERReader* reader = (HOMERReader*)fReader;
  
  ULong_t blk;
  Char_t* trkID = "SGESKART";      // TRAKSEGS
  blk = reader->FindBlockNdx( trkID, " CPT", 0xFFFFFFFF );
  
  while ( blk != ~(ULong_t)0 ) {

    HLTDebug( "Found tracks in block %lu\n", blk );

    const AliHLTTPCTrackletData* trackData = (const AliHLTTPCTrackletData*)reader->GetBlockData( blk );

    if ( !trackData ) {
      HLTError( "No track data for block %lu\n", blk );
      blk = reader->FindBlockNdx( trkID, " CPT", 0xFFFFFFFF, blk+1 );
      continue;
    }
    
    ULong_t spec = reader->GetBlockDataSpec( blk );	
    Int_t slice = AliHLTTPCDefinitions::GetMinSliceNr( spec );  
    Int_t patchmin = AliHLTTPCDefinitions::GetMinPatchNr( spec );
    Int_t patchmax = AliHLTTPCDefinitions::GetMaxPatchNr( spec );
    
    HLTDebug( "%lu tracks found for slice %u - patch %u-%u\n", trackData->fTrackletCnt, slice, patchmin, patchmax );

    fTracksPerSlice[slice] = trackData->fTrackletCnt; // Fill number if tracks per slice
	
    AliHLTTPCTrackSegmentData* tmp42 = (AliHLTTPCTrackSegmentData*) &(trackData->fTracklets[0]);
    fTracks->FillTracks( trackData->fTrackletCnt, tmp42, slice );
    
    blk = reader->FindBlockNdx( trkID, " CPT", 0xFFFFFFFF, blk+1 );	
  }
    
  SetupTracks(); 

#else
  HLTFatal("HOMER raeder not available");
#endif // defined(HAVE_HOMERREADER) 
}

//----------------------------------------------------------------------------------------------------
//                 GETTER
//____________________________________________________________________________________________________
Int_t AliHLTTPCDisplayMain::GetGlobalTrack(Int_t slice){
    // get global track out of the "selectTrack" parameters
    Int_t currenttrack= -1;
    Int_t trackcounter = 0;
    Int_t ntracks = fTracks->GetNTracks();
    
    if ( slice == fSelectTrackSlice) {
	for(Int_t j=0; j<ntracks; j++) { 	
	    
	    AliHLTTPCTrack *gtrack = fTracks->GetCheckedTrack(j); 
	    if(!gtrack) continue;
	    
	    // --- CHECK if track is should be drawn
	    // select Single Track
	    if(gtrack->GetSector() != fSelectTrackSlice) continue;
	    
	    if (trackcounter != fSelectTrack){
		trackcounter++;  
		continue;
	    }
	    trackcounter++;
	    // +++
	    //if((fPtThreshold > 0) && (gtrack->GetPt()< fPtThreshold)) continue;
	    //if(gtrack->GetNHits() < fMinHits) continue;
	    
	    currenttrack = j;
	    break;
	}
    }

    return currenttrack;
}

//----------------------------------------------------------------------------------------------------
//                 EVENTS
//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::ExecPadEvent(Int_t event, Int_t px, Int_t py, TObject *selected){
   TCanvas *c = (TCanvas *) gTQSender;

   if (event == 11 &&selected->InheritsFrom("TH2F")) {
	TH2F *hist = (TH2F*) selected;

	Int_t binx = hist->GetXaxis()->FindBin(c->AbsPixeltoX(px)) -1 ;
//	Int_t biny = hist->GetYaxis()->FindBin(c->AbsPixeltoY(py)) -1;

        fPadCallback(fPt2Gui, binx);
  }
}


