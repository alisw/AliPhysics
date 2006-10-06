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

    // Set B-Field
    AliHLTTPCTransform::SetBField( 0.8 );

    //callback handler
    fPt2Gui = pt2GUI;
    fPadCallback = pt2Function;
	
    fReader = NULL;

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

    fConnect = kFALSE;

    fExistsRawData = kFALSE;
    fExistsClusterData = kFALSE;
    fExistsTrackData = kFALSE;

    AliHLTTPCTransform::SetBField(0.4);

    memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
    memset(fNcl, 0, 36*6*sizeof(UInt_t)); 

    fTracks = NULL;

    fTheta = 90.;
    fPhi = 0.;

    fSwitch3DCluster = kFALSE;
    fSwitch3DTracks = kFALSE;
    fSwitch3DPadRow = kFALSE;
    fSwitch3DGeometry = kTRUE;


    // fMinHits = 0;
    // fPtThreshold = 0.;

    fCutHits = 0;
    fCutPt = 0.;
    fCutPsi = 0.;
    fCutLambda = 0.;
    fCut_S = 0.;
    fCutPadrow = 159;

    fNPads = AliHLTTPCTransform::GetNPads(0);
    fPad = -1;
    fPadRow = 0;
    fTimebin = 0;

    fAllTimebins = kTRUE;
    fSplitPadRow = kFALSE;
    
    fSlicePadRow = 0; 
    fSelectTrack = -1;
    fSelectTrackSlice = 0;
    fSelectTrackSwitch = kFALSE;
    fSelectCluster = 0;

    fMinSlice = 0;
    fMaxSlice = 35;
    fSlicePair = kFALSE;

    SetSliceArray();

    fBackColor = 1;
    fLineColor = 0;
    fKeepView = kTRUE;

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

    Char_t* defaultGeometry=NULL;
#if defined(DEFAULT_GEOMETRY)
    defaultGeometry=DEFAULT_GEOMETRY;
#endif
    if (gfile!=NULL) {
      HLTDebug("probing geometry file %s", gfile);
      ifstream test(gfile);
      if (test.fail()) {
	HLTWarning("unable to find geometry file %s, using default file", gfile);
	gfile=defaultGeometry;
      }
      test.close();
    } else {
      HLTDebug("using default geometry file %s", gfile, defaultGeometry);
      gfile=defaultGeometry;
    }
    if (gfile==NULL) {
      HLTError("geometry file missing");
      return -EINVAL;
    }
#if defined(HAVE_HOMERREADER) 
    // -- input datatypes , reverse
    Char_t* spptID="SRETSULC";
    Char_t* trkID = "SGESKART";
    Char_t* padrowID = "KPWR_LDD";
     
    Int_t ret;
    Int_t tries = 0;
    
    while(1) {
	HOMERReader* reader = new HOMERReader( cnt, (const char**) hostnames, ports );

	ret=reader->GetConnectionStatus();

	if ( ret ) {	
	    Int_t ndx = reader->GetErrorConnectionNdx();
	    if ( ndx < (Int_t) cnt) printf( "Error establishing connection to TCP source %s:%hu: %s (%d)\n", hostnames[ndx], ports[ndx], strerror(ret), ret );
	    else printf( "Error establishing connection to unknown source with index %d: %s (%d)\n",ndx, strerror(ret), ret );

	    delete reader;
	    reader = NULL;
	    return -1;
	}

	//cout << "TRY = " << tries << " || " ;
	tries++;

	Int_t  ret1 = reader->ReadNextEvent();
	if ( ret1 ) {
	    Int_t ndx = reader->GetErrorConnectionNdx();
	    printf( "Error reading event from source %d: %s (%d)\n", ndx, strerror(ret1), ret1 );
	    delete reader;
	    reader = NULL;
	    continue;
	}

	unsigned long blockCnt = reader->GetBlockCnt();  // ####
#if DEBUG
	printf( "Event 0x%016LX (%Lu) with %lu blocks\n", (ULong64_t)reader->GetEventID(), (ULong64_t)reader->GetEventID(), blockCnt );
#endif
	if (blockCnt == 0){
	    delete reader;
	    reader = NULL;
	    continue;
	}
	
	fReader = (void*) reader;

	// Switch on the Display parts according what data types are present
	for ( unsigned long i = 0; i < blockCnt; i++ ) {
	    char tmp1[9];
	    memset( tmp1, 0, 9 );
	    void *tmp11 = tmp1;
	    ULong64_t* tmp12 = (ULong64_t*)tmp11;
	    *tmp12 = reader->GetBlockDataType( i );
	    if (!strcmp(tmp1,padrowID)) fExistsRawData = kTRUE;
	    if (!strcmp(tmp1,spptID)) fExistsClusterData = kTRUE;
	    if (!strcmp(tmp1,trkID)) fExistsTrackData = kTRUE;   	  
	}

	break;
    }
    fConnect = kTRUE;  

    // INITIALIZE TPC DISPLAY Classes  -- IMPORTANT... don't change the order of them
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
    HLTFatal("HOMER raeder not available during package configuration, libAliHLTTPCDisplay compiled without HOMER support");
    return -EFAULT;
#endif // defined(HAVE_HOMERREADER) 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::ReadData(Bool_t nextSwitch){
#if defined(HAVE_HOMERREADER) 
    // READ CLUSTER and TRACK data
    HOMERReader* reader = (HOMERReader*)fReader;

    for(Int_t ii=0;ii<36;ii++) fTracksPerSlice[ii] = 0;

    // -- input datatypes , reverse
    Char_t* spptID="SRETSULC";
    Char_t* trkID = "SGESKART";
    Int_t ret;
    ULong_t blk;

    if (nextSwitch) {
	while(1){
	    ret = reader->ReadNextEvent();
	    
	    if ( ret ) {
		Int_t ndx = reader->GetErrorConnectionNdx();
		printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
		continue; 
	    }
	    break;
	}
    }

    ULong_t blockCnt = reader->GetBlockCnt();
#if DEBUG
    printf( "Event 0x%016LX (%Lu) with %lu blocks\n", (ULong64_t)reader->GetEventID(), (ULong64_t)reader->GetEventID(), blockCnt );
#endif    
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
#if DEBUG
	printf( "Block %lu length: %lu - type: %s - origin: %s\n",i, reader->GetBlockDataLength( i ), tmp1, tmp2 );
#endif
    }
    
    //--------------------------------------------------------------------------------------------
    // READ CLUSTER DATA
    //-------------------------------------------------------------------------------------------- 
    blk = reader->FindBlockNdx( spptID, " CPT",0xFFFFFFFF );
    
    while ( blk != ~(ULong_t)0 ) {	    
#if DEBUG
	printf( "Found clusters block %lu\n", blk );
#endif
	const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*)reader->GetBlockData( blk );
	if ( !clusterData ) {
#if DEBUG
	    printf( "No track data for block %lu\n", blk );
#endif
	    continue;
	}
	
	ULong_t spec = reader->GetBlockDataSpec( blk );
	Int_t patch = AliHLTTPCDefinitions::GetMinPatchNr( spec );
	Int_t slice = AliHLTTPCDefinitions::GetMinSliceNr( spec );
#if DEBUG
	printf( "%lu Clusters found for slice %u - patch %u\n", clusterData->fSpacePointCnt, slice, patch );
#endif	
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

    //--------------------------------------------------------------------------------------------
    // READ TRACKS
    //-------------------------------------------------------------------------------------------- 
    if ( fTracks ) delete fTracks; 
    fTracks = new AliHLTTPCTrackArray;
    
    blk = reader->FindBlockNdx( trkID, " CPT", 0xFFFFFFFF );
	
    while ( blk != ~(ULong_t)0 ) {
#if DEBUG
	printf( "Found tracks in block %lu\n", blk );
#endif
	const AliHLTTPCTrackletData* trackData = (const AliHLTTPCTrackletData*)reader->GetBlockData( blk );
	if ( !trackData ) {
#if DEBUG
	    printf( "No track data for block %lu\n", blk );
#endif
	    continue;
	}
	
	ULong_t spec = reader->GetBlockDataSpec( blk );	
	Int_t slice = AliHLTTPCDefinitions::GetMinSliceNr( spec );  
#if DEBUG
	Int_t patchmin = AliHLTTPCDefinitions::GetMinPatchNr( spec );
	Int_t patchmax = AliHLTTPCDefinitions::GetMaxPatchNr( spec );

	printf( "%lu tracks found for slice %u - patch %u-%u\n", trackData->fTrackletCnt, slice, patchmin, patchmax );
#endif	
	fTracksPerSlice[slice] = trackData->fTrackletCnt; // Fill number if tracks per slice
	
	AliHLTTPCTrackSegmentData* tmp42 = (AliHLTTPCTrackSegmentData*) &(trackData->fTracklets[0]);
	fTracks->FillTracks( trackData->fTrackletCnt, tmp42, slice );
	
	blk = reader->FindBlockNdx( trkID, " CPT", 0xFFFFFFFF, blk+1 );	
    }
    
    SetupTracks();  
#else
    HLTFatal("HOMER raeder not available during package configuration, libAliHLTTPCDisplay compiled without HOMER support");
#endif // defined(HAVE_HOMERREADER) 
}

//____________________________________________________________________________________________________
void AliHLTTPCDisplayMain::DisplayEvent(){
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
	
	if (fPadRow == 30 || fPadRow == 90 || fPadRow == 139) maxpadrowpatch++;

	for (Int_t patch=0; patch < 6; patch++){
	    unsigned long blk;
	    AliHLTUInt32_t rawSpec = AliHLTTPCDefinitions::EncodeDataSpecification( (AliHLTUInt8_t) fSlicePadRow, (AliHLTUInt8_t) fSlicePadRow,
										    (AliHLTUInt8_t) patch,(AliHLTUInt8_t) patch );
	    // READ RAW DATA BLOCK
	    blk = reader->FindBlockNdx( rawID, " CPT", rawSpec );	

	    if ( ~(unsigned long)0 != blk ){
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
    HLTFatal("HOMER raeder not available during package configuration, libAliHLTTPCDisplay compiled without HOMER support");
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
