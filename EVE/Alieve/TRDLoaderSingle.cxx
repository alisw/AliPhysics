#include "TRDLoaderSingle.h"

//#include "AliTRDv1.h"

#include <Reve/RGTopFrame.h>

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TFile.h"

using namespace Reve;
using namespace Alieve;
using namespace std;

ClassImp(Alieve::TRDLoaderSingle)

///////////////////////////////////////////////////////////
/////////////     TRDLoaderSingle     /////////////////////
///////////////////////////////////////////////////////////


//________________________________________________________
TRDLoaderSingle::TRDLoaderSingle(const Text_t* n, const Text_t* t) : TRDLoader(n, t)
{
}


//________________________________________________________
Bool_t	TRDLoaderSingle::GoToEvent(const int ev)
{
	fEvent = ev;

	Unload();
	
	TTree *t = 0x0;
	TFile *f = new TFile(Form("%s/%s", fDir.Data(), fFilename.Data()));
	if(! f->cd(Form("Event%d", ev))){
		Error("GoToEvent()", Form("Could not find event %d in file %s.", ev, fFilename.Data()));
		return kFALSE;
	}
	
	if(kLoadHits){
		t = (TTree*)gDirectory->Get("TreeH");
		if(!t) return kFALSE;
		if(!LoadHits(t)) return kFALSE;
	}
	if(kLoadDigits){
		t = (TTree*)gDirectory->Get("TreeD");
		if(!t) return kFALSE;
		if(!LoadDigits(t)) return kFALSE;
	}
	if(kLoadClusters){
		t = (TTree*)gDirectory->Get("TreeR");
		if(!t) return kFALSE;
		if(!LoadClusters(t)) return kFALSE;
	}
	if(kLoadTracks){
		t = (TTree*)gDirectory->Get("TreeT");
		if(!t) return kFALSE;
		if(!LoadTracklets(t)) return kFALSE;
	}
	f->Close();
	
	gReve->Redraw3D();
	
	return kTRUE;
}

	
//________________________________________________________
Bool_t	TRDLoaderSingle::Open(const char *filename, const char *dir)
{
	fFilename = filename;
	fDir = dir;
	
	TObjArray *so = fFilename.Tokenize(".");

	if(((TObjString*)(*so)[0])->GetString().CompareTo("TRD") != 0){
		Error("Open()", "Filename didn't fulfill naming conventions. No TRD data.");
		return kFALSE;
	}
	if(((TObjString*)(*so)[1])->GetString().CompareTo("Hits") == 0) kLoadHits = kTRUE;
	else if(((TObjString*)(*so)[1])->GetString().CompareTo("Digits") == 0) kLoadDigits = kTRUE;
	else if(((TObjString*)(*so)[1])->GetString().CompareTo("RecPoints") == 0) kLoadClusters = kTRUE;
	else if(((TObjString*)(*so)[1])->GetString().CompareTo("Tracks") == 0) kLoadTracks = kTRUE;
	else{
		Error("Open()", "Filename didn't fulfill naming conventions. No data type specified.");
		return kFALSE;
	}
	
	return kTRUE;
}

