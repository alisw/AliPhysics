//tree setting
Bool_t isMC=0;
Bool_t is13TeV=0;
Bool_t storeTracks        =1; // if kTRUE - Barrel tracks will be stored as AliCFParticles
Bool_t storeTracklets     =0; // if kTRUE - SPD tracklets will be stored as AliCFParticles
Bool_t storeMuons         =0; // if kTRUE - muon tracks will be stored as AliCFParticles
Bool_t storeMcTracks      =0; // if kTRUE - mc particles will be stored as AliCFParticles
Bool_t storeMcTracklets   =0; // if kTRUE - Store Monte-Carlo info for tracklets
Bool_t storeMcMuons       =0; // if kTRUE - Store Monte-Carlo info for muons
Bool_t storeTrackInfo     =0; // if kTRUE - Store additional track info on tracks
Bool_t storePidInfo       =0; // if kTRUE - Store PID info for tracks

//input
TString dataDir = "/alice/data/2010/LHC10d/";
TString dataPattern = "/ESDs/pass2/AOD147/*/AliAOD.root";
TString workingDir = "lhc10d";

Int_t nRuns = 53;
Int_t runList[] = {
126432, 126425, 126424, 126422, 126409, 126408, 126407, 126406, 126405, 126404, 126403, 126359, 126352, 126351, 126285, 126284, 126283, 126168, 126160, 126158, 126097, 126090, 126088, 126082, 126081, 126078, 126073, 126008, 126007, 126004, 125855, 125851, 125850, 125849, 125848, 125847, 125844, 125843, 125842, 125633, 125632, 125630, 125628, 125296, 125134, 125101, 125100, 125097, 125085, 125023, 124751, 122375, 122374
};

//Int_t nRuns = 1;
//Int_t runList[] = {126088};

//event selection
UInt_t classbit=1<<0|1<<2|1<<3|1<<4|1<<5|1<<6;
UInt_t eventselectionbit=AliVEvent::kMB;

Double_t zvertex=10.;
Bool_t applyPhysicsSelectionCut =1; // skip events not passing selectionBit mask
Bool_t storeOnlyEventsWithMuons =0; // if kTRUE store only events with at least one muon
Bool_t storeCutBitsInTrackMask  =0; // if kTRUE modify additional bits in track mask

//track selection
Double_t ptmin=0.7;
UInt_t trackfilterbit=1<<5|1<<6|1<<8|1<<9|1<<4;
Double_t tracketacut=1.4;

//tracklet selection
Double_t dphicut=0.005;
Double_t trackletetacut=1.2;

//muon selection


//masterjob setting
Int_t ntestfiles=5;
Int_t maxinputfilenumber=15;
Int_t maxmergestage=1;
Bool_t checkcopy=kFALSE;
Bool_t keeplogs=kFALSE;
Int_t NrunsPerMaster=1;

//AliPhysics version
char *AliPhysicsVersion="vAN-20160224-1";
char *AliROOTVersion="v5-08-00-3";
char *APIVersion="V1.1x";
