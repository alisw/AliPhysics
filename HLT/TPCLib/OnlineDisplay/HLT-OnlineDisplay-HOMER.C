/*
 * Requirements: Running HLT Analysis Framework
 *               with TCP-Dump Subscribe listining on 
 *                  - Cluster: Hostname, Port
 *                  - Tracker: Hostname, Port
 *                  - Controlled FilePublisher:  Hostname, Port
 * Usage :
 * 
 * (1) USE the StartDisplayMacro.C , by adding 
 *             gROOT->Macro("StartDisplayMacro.C"); 
 *       to the rootlogon.C file      
 *
 * (2) IN  StartDisplayMacro.C subject to change are: 
 * ODH_INIT("<PATH TO GEO FILE>", "<PATH TO LIBRARIES>"); 
 *      - <PATH TO LIBRARIES> : The path where libHOMERReader_ROOT.so and libHOMERReader.so are located.
 *      - <PATH TO GEO FILE>  : The path where alice.geom is located.
 * ODH_CONNECT("<NODE>",<PORT>,"<NODE>",<PORT>,"<NODE>",<PORT>);
 *      - <NODE> : The port on which the TCPDumpSubscriber is running.
 *      - <PORT> : The port specified in the TCPDumpSubscriber in the XML configuration file. 
 *    use first pair for Cluster Data, second for Tracks Data and third for Raw Data...
 *    but it doesn't matter if you change it.
 *
 *
 * When starting root, Display and Control Bar will pop up.
 * Nothing will be displayed at the beginning.
 * - In order to display the geometry. Click <Show Geometry> first.
 *     Different sets of sectors can be selected also.
 * - In order to dislpay event's click <Next clusters>, <Next tracks> or <Next clusters and tracks>
 * Displaying PadRow's:
 * - as Histogram: 
 *         once : <Setup PadRow 20 with Histogram> 
 *         next event : <Display PadRow>
 * - in Geometry:
 *         once : <Setup PadRow 20 with Geometry> 
 *         next event : <Display PadRow>, <Display PadRow with Clusters>, <Display PadRow with Tracks> or <Display PadRow with Clusters and Tracks>
 */


// gROOT->GetInterpreter()->AddIncludePath( "$ALIHLT_TOPDIR/BASE" );
// gROOT->GetInterpreter()->AddIncludePath( "$ALIHLT_TOPDIR/TPCLib" );

#include "AliHLTDataTypes.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCDigitReader.h"

// include because of connecttofile()
#include <stdio.h>
#include <stdlib.h>
//#include <sys/time.h>
//#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <vector>
// #include "MLUCString.hpp"
// #include "MLUCCmdlineParser.hpp"
// #include "MLUCDynamicLibrary.hpp"
#include "AliHLTDataTypes.h"
#include "AliHLT_C_Component_WrapperInterface.h"



//#include "AliHLTTPCTrackArray.h"

    void* gODH_Reader; /* really HOMERReader* */
    void* gODH_Display; /* really AliHLTTPCDisplay* */


/* Dummy function to work around the problem that the interpreter does not load the first function properly... */
int ODH_Dummy() {
    return 0;
}

// ########################################################################################################################################
int HLT_OnlineDisplay_HOMER() {
    gODH_Display = NULL;
    return 0;
}

// ########################################################################################################################################
int ODH_Init( char* path_to_geom_file, char* path_to_homer_lib = NULL ) {
    // --- Load LIBS
    TString lib;
    cout << "Loading ROOT libraries (ROOTSYS has to be set)" << endl;

    gSystem->Load("libPhysics");
    gSystem->Load("libEG");
//     if(gSystem->Load("libMC")==-1) {
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
// 	}

    cout << "Loading ALICE TPC libraries (ALICE_ROOT & ALICE_TARGET have to be set)" << endl;
    gSystem->Load("libESD");
    gSystem->Load("libSTEER");
    gSystem->Load("libRAWDatabase");
    gSystem->Load("libRAWDatarec");
//     gSystem->Load("libCONTAINERS");
//     if(gSystem->Load("libTPC")!=0) {
    gSystem->Load("libTPCbase");
    gSystem->Load("libTPCrec");
    gSystem->Load("libTPCsim");
    gSystem->Load("libTPCfast");
//  	}

    cout << "Loading HLT libraries (ALIHLT_LIBDIR has to be set)" << endl;
    gSystem->Load("$(ALIHLT_LIBDIR)/libHLTbase.so");
    gSystem->Load("$(ALIHLT_LIBDIR)/libAliHLTTPC.so");
    cout << "Loading HOMER library" << endl;

    if ( path_to_homer_lib ) {
	lib = path_to_homer_lib;
	lib += "/";
    }
    else {
	lib = "";
    }
    
    lib += "libHOMERReader_ROOT.so";
    gSystem->Load( lib );


    // --- Create DISPLAY
    cout << "Creating display" << endl;

    Char_t *gfile="alice.geom"; // geometrie file
    TString geoPath;
    if ( path_to_geom_file ){	
	geoPath = path_to_geom_file;
	geoPath += "/";
    }
    else {
	geoPath = "";
    }
    geoPath += gfile;

    Int_t slices[] = { 0, 35 };
    AliHLTTPCDisplay* display = new AliHLTTPCDisplay( slices, geoPath.Data() );

    display->SetSlices();
    display->DisplayAll(0,kFALSE,kFALSE,kFALSE,kFALSE);

    gODH_Display = (void*)display;
    return 0;
}

int ODH_SetSliceRange(int minslice,int maxslice){
  // sets a range of Slices
  if (!gODH_Display ) return -1;
  AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)gODH_Display;

  display->SetSlices(minslice,maxslice);
  display->DisplayAll(0,kFALSE,kFALSE,kFALSE,kFALSE);
}

int ODH_SetSliceRange(int slice){
  // sets one slice
  if (!gODH_Display ) return -1;
  AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)gODH_Display;

  display->SetSlices(slice); 
  display->DisplayAll(0,kFALSE,kFALSE,kFALSE,kFALSE);
}

int ODH_SetSliceRange(){
  // sets all Slices
  if (!gODH_Display ) return -1;
  AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)gODH_Display;

  display->SetSlices();
  display->DisplayAll(0,kFALSE,kFALSE,kFALSE,kFALSE);
}

int ODH_SetSlicePair(int slice){
  // sets a pair of slices
  if (!gODH_Display ) return -1;
  AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)gODH_Display;

  display->SetSlicesPair(slice);
  display->DisplayAll(0,kFALSE,kFALSE,kFALSE,kFALSE);
}

int ODH_SetSlicePair(int minslice, int maxslice){
  // sets a pair of slices
  if (!gODH_Display ) return -1;
  AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)gODH_Display;

  display->SetSlicesPair(minslice,maxslice);
  display->DisplayAll(0,kFALSE,kFALSE,kFALSE,kFALSE);
}

int ODH_SetInvert(){
  // sets a pair of slices
  if (!gODH_Display ) return -1;
  AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)gODH_Display;

  display->SetInvert();
  display->DisplayAll(0,kFALSE,kFALSE,kFALSE,kFALSE);
}

int ODH_SetDrawGeo(){
  // sets a pair of slices
  if (!gODH_Display ) return -1;
  AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)gODH_Display;

  display->SetDrawGeo();
  display->DisplayAll(0,kFALSE,kFALSE,kFALSE,kFALSE);
}

// ########################################################################################################################################
int ODH_Connect( char* hostname_clusters, short port_clusters, char* hostname_tracks, short port_tracks, char* hostname_raw, short port_raw ) {
    HOMERReader* reader;
    char* hostnames[] = { "NULL", "NULL","NULL" };
    unsigned short ports[] = { 0, 0, 0 };
    int cnt=0;
    if ( hostname_clusters ){
	hostnames[cnt] = hostname_clusters;
	ports[cnt] = port_clusters;
	cnt++;
    }
    if ( hostname_tracks ){
	hostnames[cnt] = hostname_tracks;
	ports[cnt] = port_tracks;
	cnt++;
    }
    if ( hostname_raw ){
	hostnames[cnt] = hostname_raw;
	ports[cnt] = port_raw;
	cnt++;
    }

    reader = new HOMERReader( cnt, hostnames, ports );
    int ret;
    ret=reader->GetConnectionStatus();
    if ( ret ) {
	int ndx = reader->GetErrorConnectionNdx();
	if ( ndx < cnt ) {
	    printf( "Error establishing connection to TCP source %s:%hu: %s (%d)\n", 
		    hostnames[ndx], ports[ndx], 
		    strerror(ret), ret );
	}
	else {	    
	    printf( "Error establishing connection to unknown source with index %d: %s (%d)\n", 
		    ndx, 
		    strerror(ret), ret );
	}
	delete reader;
	reader = NULL;
    }
    gODH_Reader = (void*)reader;
    return ret;
}

// ########################################################################################################################################
int ODH_Disconnect() {
    if ( !gODH_Reader )
	return 0;
    HOMERReader* reader = (HOMERReader*)gODH_Reader;
    delete reader;
    gODH_Reader = NULL;
    return 0;
}

// ########################################################################################################################################
Int_t ODH_SetupPadRow(Int_t histswitch, Int_t slice, Int_t padrow){  
    // PADROW 0 - 158
    if (!gODH_Display ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)gODH_Display;

    // Setup Geometry:  histswitch = 0;
    // Setup Histogram: histswitch = 1|2|3;
    display->SetupPadRow(histswitch,slice,padrow);

    return 0;
}

// ########################################################################################################################################
int ODH_DisplayNextEvent( Bool_t clusterswitch, Bool_t trackswitch, Bool_t padrowswitch, Float_t* clustersEtaRange = NULL ) {
    
    if ( !gODH_Reader || !gODH_Display ) return -1;
    
    // -- input datatypes , reverse
    char* spptID="SRETSULC";
    //char* spptID="STPECAPS";
    char* trkID = "SGESKART";
    //char* trkID = "SGESCART";
    char* padrowID = "KPWR_LDD";

    Int_t minHits = 0;
    Bool_t x3don = kFALSE;

    HOMERReader* reader = (HOMERReader*)gODH_Reader;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)gODH_Display;
    int ret = reader->ReadNextEvent();

    if ( ret ) {
	int ndx = reader->GetErrorConnectionNdx();
	printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
	ODH_DisplayNextEvent( clusterswitch, trackswitch, padrowswitch, clustersEtaRange);
	return ret;
    }

    unsigned long blockCnt = reader->GetBlockCnt();
    printf( "Event 0x%016LX (%Lu) with %lu blocks\n", (ULong64_t)reader->GetEventID(), (ULong64_t)reader->GetEventID(), blockCnt );

    for ( unsigned long i = 0; i < blockCnt; i++ ) {
	char tmp1[9], tmp2[5];
	memset( tmp1, 0, 9 );
	memset( tmp2, 0, 5 );
	void *tmp11 = tmp1;
	ULong64_t* tmp12 = (ULong64_t*)tmp11;
	*tmp12 = reader->GetBlockDataType( i );
	void *tmp21 = tmp2;
	ULong_t* tmp22 = (ULong_t*)tmp21;
	*tmp22 = reader->GetBlockDataOrigin( i );
 	printf( "Block %lu length: %lu - type: %s - origin: %s\n",
 		i, reader->GetBlockDataLength( i ), tmp1, tmp2 );
    }
    
    display->ResetDisplay();

    // -------------- CLUSTER
    if ( clusterswitch ) {
	unsigned long blk = reader->FindBlockNdx( spptID, " CPT",0xFFFFFFFF );
	printf( "blk: %lu\n", blk );
	while ( blk != ~(unsigned long)0 ) {	    
	    printf( "Found clusters block %lu\n", blk );
	    const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*)reader->GetBlockData( blk );
	    unsigned long dataLen = reader->GetBlockDataLength( blk );
	    
	    ULong_t spec = reader->GetBlockDataSpec( blk );
	    Int_t patch = AliHLTTPCDefinitions::GetMinPatchNr( spec );
	    Int_t slice = AliHLTTPCDefinitions::GetMinSliceNr( spec );
	    printf( "%lu Clusters found for slice %u - patch %u\n", clusterData->fSpacePointCnt, slice, patch );

	    void* tmp30 = clusterData;
	    Byte_t* tmp31 = (Byte_t*)tmp30;
	    unsigned long offset;
	    offset = sizeof(clusterData->fSpacePointCnt);
	    if ( offset <= reader->GetBlockTypeAlignment( blk, 1 ) )
		offset = reader->GetBlockTypeAlignment( blk, 1 );
	    tmp31 += offset;
	    tmp30 = tmp31;
	    AliHLTTPCSpacePointData* tmp32 = (AliHLTTPCSpacePointData*)tmp30;

	    display->SetupClusterDataForPatch( slice, patch, clusterData->fSpacePointCnt, tmp32 );

	    blk = reader->FindBlockNdx( spptID, " CPT", 0xFFFFFFFF, blk+1 );
	    printf( "blk: %lu\n", blk );
	}
    }

    // -------------- PADROW
    if ( padrowswitch ) {
	Int_t padrow = display->GetPadrow();
	Int_t patch =  AliHLTTPCTransform::GetPatch(padrow);
	Int_t maxpatch = patch;
	AliHLTUInt8_t padslice = (AliHLTUInt8_t) display->GetSlice();

	if (padrow == 30 || padrow == 90 || padrow == 139) maxpatch++;
      
	for (Int_t tpatch=patch;tpatch <= maxpatch;tpatch++){

	  
	    AliHLTUInt32_t padrowSpec = AliHLTTPCDefinitions::EncodeDataSpecification( padslice, padslice,(AliHLTUInt8_t) tpatch,(AliHLTUInt8_t) tpatch );
	    
	    unsigned long blk;
	    
	    // READ RAW DATA BLOCK
	    blk = reader->FindBlockNdx( padrowID, " CPT", padrowSpec );
	    printf( "Raw Data found for slice %u/patch %u\n", padslice, tpatch );
	    
	    unsigned long rawDataBlock = reader->GetBlockData( blk );
	    unsigned long rawDataLen = reader->GetBlockDataLength( blk );
	    
	    // READ CLUSTERS BLOCK
	    blk = reader->FindBlockNdx( spptID, " CPT", padrowSpec );
	    
	    const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*)reader->GetBlockData( blk );
	    
	    printf( "%lu Clusters found for slice %u - patch %u\n", clusterData->fSpacePointCnt, padslice, tpatch );
	    
	    void* tmp30 = clusterData;
	    Byte_t* tmp31 = (Byte_t*)tmp30;
	    unsigned long offset;
	    offset = sizeof(clusterData->fSpacePointCnt);
	    if ( offset <= reader->GetBlockTypeAlignment( blk, 1 ) )
		offset = reader->GetBlockTypeAlignment( blk, 1 );
	    tmp31 += offset;
	    tmp30 = tmp31;
	    AliHLTTPCSpacePointData* tmp32 = (AliHLTTPCSpacePointData*)tmp30;
	    
	    // DISPLAY RAW DATA AND CLUSTERS
	    display->FillPadRow( tpatch, rawDataBlock, rawDataLen, clusterData->fSpacePointCnt, tmp32);
	}


    }

    AliHLTTPCTrackArray* trackArray = NULL;
    trackArray = new AliHLTTPCTrackArray;
    if ( !trackArray ) {
	printf( "No track array\n" );
	return -1;
    }

    // -------------- TRACKS
    if ( trackswitch ) {
	unsigned long blk = reader->FindBlockNdx( trkID, " CPT", 0xFFFFFFFF );

	while ( blk != ~(unsigned long)0 ) {
	    printf( "Found tracks in block %lu\n", blk );
	    const AliHLTTPCTrackletData* trackData = (const AliHLTTPCTrackletData*)reader->GetBlockData( blk );
	    if ( !trackData ) {
		printf( "No track data for block %lu\n", blk );
		continue;
	    }
	    
	    ULong_t spec = reader->GetBlockDataSpec( blk );
	    Int_t patchmin = AliHLTTPCDefinitions::GetMinPatchNr( spec );
	    Int_t patchmax = AliHLTTPCDefinitions::GetMaxPatchNr( spec );
 	    Int_t slice = AliHLTTPCDefinitions::GetMinSliceNr( spec );  

	    printf( "%lu tracks found for slice %u - patch %u-%u\n", trackData->fTrackletCnt, slice, patchmin, patchmax );

	    void* tmp40 = trackData;
	    Byte_t* tmp41 = (Byte_t*)tmp40;
	    unsigned long offset;

	    offset = sizeof(trackData->fTrackletCnt);

	    if ( offset <= reader->GetBlockTypeAlignment( blk, 1 ) ) offset = reader->GetBlockTypeAlignment( blk, 1 );

	    tmp41 += offset;
	    tmp40 = tmp41;
	    AliHLTTPCTrackSegmentData* tmp42 = (AliHLTTPCTrackSegmentData*)tmp40;

	    trackArray->FillTracks( trackData->fTrackletCnt, tmp42, slice );

 	    blk = reader->FindBlockNdx( trkID, " CPT", 0xFFFFFFFF, blk+1 );
	}


    }

    display->SetTracks( trackArray );
    
    // DISPLAY
    display->DisplayAll( minHits, clusterswitch, trackswitch, x3don, 0. /*, clustersEtaRange */ );
    if ( padrowswitch) display->DrawPadRow(x3don);
    
    if ( trackArray )  delete trackArray;   

    return 0;
}

// ########################################################################################################################################
// ########################################################################################################################################
// ########################################################################################################################################
int ODH_SimpleDisplay(char * infiles[7],unsigned int slice){
    if ( !gODH_Display ) return -1;

    FILE * inFH = 1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)gODH_Display;	
    int i = 1;
    char* error = NULL;
    int errorArg = -1;
    int errorParam = -1;

    AliHLTComponent_DataType dataType;
    AliHLTUInt32_t dataSpec = ~(AliHLTUInt32_t)0;
    dataType.fStructSize = sizeof(dataType);
    memset( dataType.fID, '*', 8 );
    memset( dataType.fOrigin, '*', 4 );
    AliHLTComponent_BlockData blocks[7];
    unsigned long totalRead = 0;
  
    // --- get INFILE
    for (int ii = 0;ii<7;ii++){

	inFH = fopen( infiles[ii], "r");
	if ( inFH == -1 ) {
	    error = "Unable to open input file for reading";
	    errorArg = i;
	    errorParam = i+1;
	    break;
	}
    
	AliHLTComponent_BlockData newBlock;
	newBlock.fStructSize = sizeof(AliHLTComponent_BlockData);

	newBlock.fShmKey.fStructSize = sizeof(AliHLTComponent_ShmData);
	newBlock.fShmKey.fShmType = gkAliHLTComponent_InvalidShmType;
	newBlock.fShmKey.fShmID = gkAliHLTComponent_InvalidShmID;
	
	fseek( inFH, 0, SEEK_END );
	long sz = ftell(inFH);
	
	fseek( inFH, 0, SEEK_SET);
	newBlock.fPtr = new uint8[ sz ];
	
	if ( !newBlock.fPtr ) {
	    fprintf( stderr, "Out of memory trying to allocate memory for input file '%s' of %lu bytes.\n",
		     infile, sz );
	    return -1;
	}

	unsigned long curSize = 0;
	int ret=1;
	void * tmp50 = newBlock.fPtr;
	uint8* tmp51 = (uint8*)tmp50;
	while ( ret>0 ) {
	    ret = fread(tmp51+curSize,1, sz-curSize,inFH);

 	    if ( ret >= 0 ) {
		curSize += (unsigned long)ret;
		if ( curSize >= (unsigned long)sz ) {
		    newBlock.fSize = sz;
		    blocks[ii] =  newBlock;
		    break;
		}
	    }
	    else {
		fprintf( stderr, "%s error reading data from input file after %lu bytes.\n", infile, curSize );
		return -1;
	    }
	}
    }
  
    //  --- Setup CLUSTERS
    for (int ii = 1;ii<7;ii++){
	unsigned int patch = ii - 1;
	const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*) blocks[ii].fPtr;
	
	if ( !clusterData )	{
	    printf( "No ClusterData\n" );
	    return -1;
	}
	
		void* tmp30 = clusterData;
	Byte_t* tmp31 = (Byte_t*)tmp30;
	unsigned long offset;
	
	offset = sizeof(clusterData->fSpacePointCnt);
	
	tmp31 += offset;
	tmp30 = tmp31;
	AliHLTTPCSpacePointData* tmp32 = (AliHLTTPCSpacePointData*)tmp30;
	
	display->SetupClusterDataForPatch( slice, patch, clusterData->fSpacePointCnt, tmp32 );
    }

    // --- Setup TRACKS
    AliHLTTPCTrackArray* trackArray = NULL;
    trackArray = new AliHLTTPCTrackArray;
    
    if ( !trackArray )	{
	printf( "No track array\n" );
	return -1;
    }

    const AliHLTTPCTrackletData* trackData = (const AliHLTTPCTrackletData*) blocks[0].fPtr;
    if ( !trackData )
    {
	printf( "No track data for block %lu\n", blk );
	continue;
    }
   
    void* tmp40 = trackData;
    Byte_t* tmp41 = (Byte_t*)tmp40;
    unsigned long offset;
    
    offset = sizeof(trackData->fTrackletCnt);
    
    tmp41 += offset;
    tmp40 = tmp41;
    AliHLTTPCTrackSegmentData* tmp42 = (AliHLTTPCTrackSegmentData*)tmp40;
    
    trackArray->FillTracks( trackData->fTrackletCnt, tmp42, slice );
    
    display->SetTracks( trackArray );

    // --- Display 
    Int_t minHits = 0; 

    display->SetSlices(0,4);
    display->DisplayAll( minHits, clusterswitch, trackswitch, kFALSE, 0. /*, clustersEtaRange */ );
    if ( trackArray ) delete trackArray;
    
    return 0;

} // END ODH_SimpleDisplay
