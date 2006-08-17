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

#include "AliHLTDataTypes.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCDigitReader.h"
//#include "AliHLTTPCTrackArray.h"

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


void* fODHReader;  /* really HOMERReader* */
void* fODHDisplay; /* really AliHLTTPCDisplay* */

Int_t fTracksPerSlice[36];   // TrackCount per slice

Float_t fTheta;
Float_t fPhi;

//vector<QString> fHostnames;
//vector<QString> fPorts;

class AliHLTTPCTrackArray;
AliHLTTPCTrackArray* fTrackArray;

TCanvas *fCanvasPad;
TCanvas *fCanvasPadRow;
TCanvas *fCanvas3D;
TCanvas *fCanvasCharge;
TCanvas *fCanvasResiduals;

/* Dummy function to work around the problem that the interpreter does not load the first function properly... */
int ODH_Dummy() {
    return 0;
}

// ########################################################################################################################################
int HLT_OnlineDisplay_HOMER() {
    fODHDisplay = NULL;
    fODHReader = NULL; 
    fTrackArray = NULL;  
    fTheta = 90.;
    fPhi = 0.;

    fCanvasPad = NULL;
    fCanvasPadRow = NULL;
    fCanvas3D = NULL;
    fCanvasCharge = NULL;
    fCanvasResiduals = NULL;
    return 0;
}

// ########################################################################################################################################
int ODH_Init() {
    // --- Load LIBS
    cout << "Loading ROOT libraries (ROOTSYS has to be set)" << endl;
    gSystem->Load("libPhysics");
    gSystem->Load("libEG");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libGpad");

    cout << "Loading ALICE TPC libraries (ALICE_ROOT & ALICE_TARGET have to be set)" << endl;
    gSystem->Load("libESD");
    gSystem->Load("libSTEER");
    // AliRoot versions v4-04-Release or higher have a changed RAW library layout
    // TODO: check this at configure time
#ifndef ALIRAW_4_04
    gSystem->Load("libRAWData");
#else //ALIRAW_4_04
    gSystem->Load("libRAWDatabase");
    gSystem->Load("libRAWDatarec");
#endif //ALIRAW_4_04
//     gSystem->Load("libCONTAINERS");
//     if(gSystem->Load("libTPC")!=0) {
    gSystem->Load("libTPCbase");
    gSystem->Load("libTPCrec");
    gSystem->Load("libTPCsim");
    gSystem->Load("libTPCfast");

    cout << "Loading HLT libraries (ALIHLT_LIBDIR has to be set)" << endl;
    gSystem->Load("$(ALIHLT_LIBDIR)/libHLTbase.so");
    gSystem->Load("$(ALIHLT_LIBDIR)/libAliHLTTPC.so");

    cout << "Loading HOMER library" << endl;
    gSystem->Load("${ALIHLT_DC_DIR}/lib/Linux-i686/libHOMERReader_ROOT.so" );

    // --- Create DISPLAY
    cout << "Creating display" << endl;
    fCanvasPad = new TCanvas("fCanvasPad","HLT Online Display - Pad Display",900,900);
    fCanvasPad->Divide(1,3);
    fCanvasPadRow = new TCanvas("fCanvasPadRow","HLT Online Display - PadRow Display",900,900);
    fCanvas3D = new TCanvas("fCanvas3D","HLT Online Display - 3D Display",900,900);
    fCanvasCharge = new TCanvas("fCanvasCharge","HLT Online Display - Charge Display",900,900);
    fCanvasResiduals = new TCanvas("fCanvasResiduals","HLT Online Display - Residuals Display",900,900);

    fTheta = 90.;
    fPhi = 0.;

    // --- Load Geometry
    TString geoPath;

    Char_t * geometryPath = getenv("ALIHLT_DC_DIR");
    Char_t * geometryFile = "alice.geom";

    if (geometryPath) {
	geoPath = geometryPath;
	geoPath += "/share/TPC-OnlineDisplayHOMER";
    }
    else geoPath = ".";
    
    geoPath += "/";
    geoPath += geometryFile;

    AliHLTTPCDisplay* display = new AliHLTTPCDisplay( (Char_t*)geoPath.Data() );

    fODHDisplay = (void*)display;

    display->SetupHist();

    return 0;
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
    fODHReader = (void*)reader;

    printf( "CONNECTED\n");
    return ret;
}

// ########################################################################################################################################
Int_t ODH_Disconnect() {
    if ( !fODHReader )
	return 0;
    HOMERReader* reader = (HOMERReader*)fODHReader;
    delete reader;
    fODHReader = NULL;
  
    if ( fTrackArray )  delete fTrackArray;   
    fTrackArray = NULL;  

    return 0; 
} 


// ########################################################################################################################################
Int_t ODH_DisplayEvent(Bool_t nextSwitch=kTRUE, Bool_t threeDSwitch=kTRUE, Bool_t PadRowSwitch= kTRUE){

    if ( !fODHReader || !fODHDisplay ) {
    	printf( "ERROR no READER or DISPLAY" );
	return -1;
    }
    
    HOMERReader* reader = (HOMERReader*)fODHReader;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;

    // -- input datatypes , reverse
    char* spptID="SRETSULC";
    char* trkID = "SGESKART";
    char* padrowID = "KPWR_LDD";
    Int_t ret;

    if (nextSwitch) {
	ret = reader->ReadNextEvent();
    
	if ( ret ) {
	    int ndx = reader->GetErrorConnectionNdx();
	    printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
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
	    printf( "Block %lu length: %lu - type: %s - origin: %s\n",i, reader->GetBlockDataLength( i ), tmp1, tmp2 );
	}
	
	for(Int_t ii=0;ii<36;ii++) fTracksPerSlice[ii] = 0;
    }
    
    //--------------------------------------------------------------------------------------------
    // READ CLUSTER DATA
    //-------------------------------------------------------------------------------------------- 
    if (nextSwitch) {
	unsigned long blk = reader->FindBlockNdx( spptID, " CPT",0xFFFFFFFF );
    
	while ( blk != ~(unsigned long)0 ) {	    
	    printf( "Found clusters block %lu\n", blk );
	    const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*)reader->GetBlockData( blk );
	    if ( !clusterData ) {
		printf( "No track data for block %lu\n", blk );
		continue;
	    }
	    
	    ULong_t spec = reader->GetBlockDataSpec( blk );
	    Int_t patch = AliHLTTPCDefinitions::GetMinPatchNr( spec );
	    Int_t slice = AliHLTTPCDefinitions::GetMinSliceNr( spec );
	    printf( "%lu Clusters found for slice %u - patch %u\n", clusterData->fSpacePointCnt, slice, patch );
	    
	    void* tmp30 = (void*)clusterData;
	    Byte_t* tmp31 = (Byte_t*)tmp30;
	    unsigned long offset;
	    offset = sizeof(clusterData->fSpacePointCnt);
	    if ( offset <= reader->GetBlockTypeAlignment( blk, 1 ) ) offset = reader->GetBlockTypeAlignment( blk, 1 );
	    tmp31 += offset;
	    tmp30 = tmp31;
	    AliHLTTPCSpacePointData* tmp32 = (AliHLTTPCSpacePointData*)tmp30;
	    
	    display->SetupCluster( slice, patch, clusterData->fSpacePointCnt, tmp32 );
	    
	    blk = reader->FindBlockNdx( spptID, " CPT", 0xFFFFFFFF, blk+1 );
	}
    }
    //--------------------------------------------------------------------------------------------
    // READ TRACKS
    //-------------------------------------------------------------------------------------------- 
    if (nextSwitch) {
	if ( fTrackArray ) delete fTrackArray; 
	fTrackArray = new AliHLTTPCTrackArray;
	
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
	    
	    fTracksPerSlice[slice] = trackData->fTrackletCnt;
	    
	    void* tmp40 = (void*)trackData;
	    Byte_t* tmp41 = (Byte_t*)tmp40;
	    unsigned long offset;
	    offset = sizeof(trackData->fTrackletCnt);
	    if ( offset <= reader->GetBlockTypeAlignment( blk, 1 ) ) offset = reader->GetBlockTypeAlignment( blk, 1 );
	    tmp41 += offset;
	    tmp40 = tmp41;
	    AliHLTTPCTrackSegmentData* tmp42 = (AliHLTTPCTrackSegmentData*)tmp40;
	    
	    fTrackArray->FillTracks( trackData->fTrackletCnt, tmp42, slice );
	    
	    blk = reader->FindBlockNdx( trkID, " CPT", 0xFFFFFFFF, blk+1 );	
	}
	
	display->SetupTracks( fTrackArray );  
    }

    //--------------------------------------------------------------------------------------------
    // READ RAW DATA for PADROW HISTOGRAM
    //--------------------------------------------------------------------------------------------
    if (PadRowSwitch) {
	AliHLTUInt8_t padslice = (AliHLTUInt8_t) display->GetSlicePadRow();
	Int_t padrow = display->GetPadRow();
	Int_t patch =  AliHLTTPCTransform::GetPatch(padrow);
	Int_t maxpatch = patch;
	
	display->ResetHistPadRow();
	
	if (padrow == 30 || padrow == 90 || padrow == 139) maxpatch++;
	
	for (Int_t tpatch=patch;tpatch <= maxpatch;tpatch++) {
	    unsigned long blk;
	    AliHLTUInt32_t padrowSpec = AliHLTTPCDefinitions::EncodeDataSpecification( padslice, padslice,(AliHLTUInt8_t) tpatch,(AliHLTUInt8_t) tpatch );
	    
	    // READ RAW DATA BLOCK - READ CLUSTERS BLOCK
	    blk = reader->FindBlockNdx( padrowID, " CPT", padrowSpec );	

	    if ( ~(unsigned long)0 != blk ){
		printf( "Raw Data found for slice %u/patch %u\n", padslice, tpatch );
		unsigned long rawDataBlock = (unsigned long) reader->GetBlockData( blk );
		unsigned long rawDataLen = reader->GetBlockDataLength( blk );

		display->FillPadRow( tpatch, rawDataBlock, rawDataLen);
	    }
	}
    
	if (fCanvasPadRow){
	    fCanvasPadRow->cd();
	    display->DrawHistPadRow();
	    fCanvasPadRow->Update();
	}
	
	if (fCanvasPad){
	    fCanvasPad->cd(1);
	    display->DrawHistPad1();
	    fCanvasPad->GetCanvas()->Update();
	    
	    fCanvasPad->cd(2);
	    display->DrawHistPad2();
	    fCanvasPad->GetCanvas()->Update();
	    
	    fCanvasPad->cd(3);
	    display->DrawHistPad3();
	    fCanvasPad->GetCanvas()->Update();
	}
    }

    //--------------------------------------------------------------------------------------------
    // RESET HISTOGRAMS
    //--------------------------------------------------------------------------------------------
    if ( display->Get3DSwitchCluster() ) display->ResetHistCharge();
    if ( display->Get3DSwitchTracks() ) display->ResetHistResiduals();

    //--------------------------------------------------------------------------------------------
    // DRAW 3D
    //--------------------------------------------------------------------------------------------
    if ( (threeDSwitch || (PadRowSwitch && display->Get3DSwitchPadRow())) && fCanvas3D ){

	fCanvas3D->cd(); 
	fCanvas3D->Clear();
	fCanvas3D->SetFillColor(display->GetBackColor());
	
	if ( !display->GetKeepView() ){
	    fCanvas3D->SetTheta(fTheta);
	    fCanvas3D->SetPhi(fPhi);
	}
	
	display->Draw3D();
	
	fCanvas3D->Modified();
	fCanvas3D->Update();		    
    }
    
    //--------------------------------------------------------------------------------------------
    // DRAW RESIDUALS
    //--------------------------------------------------------------------------------------------
    if ( display->Get3DSwitchTracks() && fCanvasResiduals) {
	fCanvasResiduals->cd();  
	fCanvasResiduals->Clear();
	fCanvasResiduals->Divide(1,2);
	fCanvasResiduals->cd(1);
	display->DrawHistResiduals(kTRUE);
	fCanvasResiduals->cd(2);
	display->DrawHistResiduals(kFALSE);
	fCanvasResiduals->Update();
    }
    //--------------------------------------------------------------------------------------------
    // DRAW CHARGE
    //--------------------------------------------------------------------------------------------
    if ( display->Get3DSwitchCluster() && fCanvasCharge != NULL) {
	fCanvasCharge->cd();       
	fCanvasCharge->Clear();

	display->DrawHistCharge();

	fCanvasCharge->Update();
    }

    return 0;
}

// ####################################################################################################
//       SETTER - minor functions
// ####################################################################################################
Int_t ODH_SetSliceRange(){
    // -- sets all Slices
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;
    
    display->SetSlices();
    fTheta = 90.;

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_SetSliceRange(Int_t slice){
    // -- sets one slice
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;
    
    display->SetSlices(slice);  
    fTheta = 0.;

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_SetSliceRange(Int_t minslice,Int_t maxslice){
    // -- sets a range of Slices
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;
    
    display->SetSlices(minslice,maxslice);  
    fTheta = 90.;

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_SetSlicePair(Int_t slice){
    // -- sets a pair of slices
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;
    
    display->SetSlicesPair(slice); 
    fTheta = 90.;
    
    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_SetSlicePair(Int_t minslice, Int_t maxslice){
    // -- sets a pair of slices
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;
    
    display->SetSlicesPair(minslice,maxslice);
    fTheta = 90.;

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_SetCluster(Bool_t used=kTRUE, Bool_t unused=kTRUE){    
    // -- set used/unuse/all cluster
    if ( !fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;
    
    //  ALL Cluster = 0  | USED Cluster = 1 | UNUSED Cluster = 2
    if (used && unused) display->SetSelectCluster(0);
    else if (used) display->SetSelectCluster(1);
    else if (unused) display->SetSelectCluster(2);

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_SetInvert(){    
    // -- set invert 3D display
    if ( !fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;

    display->SetInvert();

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_SetKeepView(Bool_t keepview){ 
    // -- keep angle view 
    if ( !fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;

    if (keepview) display->SetKeepView(kTRUE);
    else display->SetKeepView(kFALSE);

    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_SetPadRow(Int_t slice, Int_t padrow, Int_t pad = 0) {
    // -- set padrow and pad histograms
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;
    
    if( padrow < 0 || padrow > 158){
	printf("Padrow %d out of range, has to be in [0,158].",padrow);
	return -2;
    }

    display->SetSlicePadRow(slice);
    display->SetPadRow(padrow); 
    display->SetHistPadRowAxis();
    display->SetPad(pad);

    ODH_DisplayEvent(kFALSE,kTRUE,kTRUE);
    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_SetSelectTrack(Bool_t trackswitch=kFALSE, Int_t slice=0, Int_t track=0){
    // -- trackswitch : turn on /off of single track
    // -- slice : select slice for single track
    // -- track :select single track in slice
    
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;

    display->SetSelectTrackSwitch(trackswitch);
    if(!trackswitch) display->SetSelectTrack(0); 
    else {
	if ( track < 0 || track >  fTracksPerSlice[slice]) {
	    printf ("No Track %d, Maximum of Tracks in Slice %d is %d ! ",track,slice,fTracksPerSlice[slice] );
	    return -2;
	}
	
	display->SetSelectTrackSlice(slice); 
	display->SetSelectTrack(track);
    }

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_SetTrack(Int_t minhits = 0, Float_t ptthreshold = 0.){
    // -- set track properties
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;

    display->SetMinHits(minhits);
    display->SetPtThreshold(ptthreshold );
    return 0;
}
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같
Int_t ODH_Set3D(Bool_t tracks=kFALSE, Bool_t cluster=kFALSE, Bool_t padrow=kFALSE, Bool_t geometry=kFALSE){    
    // -- set 3D Display
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;
    
    if  (padrowCheckBox->isChecked()) setPadRow();
    
    display->SetSwitches(tracks, cluster, padrow, geometry);

    if (display->Get3DSwitchPadRow()) ODH_DisplayEvent(kFALSE,kTRUE,kTRUE);
    else ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}  
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같  
Int_t ODH_Set3DTracks(Bool_t on=kFALSE){
    // --  Set 3D tracks
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;

    display->Set3DSwitchTracks(on);

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}  
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같  
Int_t ODH_Set3DCluster(Bool_t on=kFALSE){
    // --  Set 3D cluster
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;

    display->Set3DSwitchCluster(on);

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}  
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같  
Int_t ODH_Set3DPadRow(Bool_t on=kFALSE){
    // --  Set 3D padrow
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;

    display->Set3DSwitchPadRow(on);

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}  
// 같같같같같같같같같같같같같같같같같같같같같같같같같같같같  
Int_t ODH_Set3DGeometry(Bool_t on=kFALSE){
    // --  Set 3D geometry
    if (!fODHDisplay ) return -1;
    AliHLTTPCDisplay* display = (AliHLTTPCDisplay*)fODHDisplay;

    display->Set3DSwitchGeometry(on);

    ODH_DisplayEvent(kFALSE,kTRUE,kFALSE);
    return 0;
}  


