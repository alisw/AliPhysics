#ifndef ALI_MUON_DATA_INTERFACE_H
#define ALI_MUON_DATA_INTERFACE_H

// Author: Artur Szostak
//  email: artur@alice.phy.uct.ac.za

#include <TObject.h>
#include <TString.h>
#include <TParticle.h>

#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliMUONData.h"

#include "AliMUONHit.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONLocalTrigger.h"


// An easy to use interface to the MUON module data stored in TreeK, TreeH, TreeS, TreeD and TreeR
// One can fetch any of the data objects with all the calls to runloader, muon loader and AliMUONData
// done behind the scenes and automatically.
// This interface in not necessarily the fastest way to fetch the data but it is the easiest.
// Note: If independant calls to the run loader, muon loader or AliMUONData objects are interspersed
// with calls to the AliMUONDataInterface to fetch data, one might need to call the Reset method
// between these method calls at some point to prevent AliMUONDataInterface from getting confused.
// This is necessary since this object assumes the state of runloader, muon loader nor AliMUONData
// has not changed between calls. If the state has changes then one must call Reset so that 
// AliMUONDataInterface refreshes what it knows about the state of the loader and AliMUONData objects.
//
class AliMUONDataInterface : public TObject
{
public:

	AliMUONDataInterface();
	~AliMUONDataInterface();
	
	// Sets all internal pointers to NULL without releasing the current runloader.
	void Reset();
	
	Int_t NumberOfEvents(TString filename, TString foldername);

	Int_t NumberOfParticles(TString filename, TString foldername, Int_t event);
	TParticle* Particle(TString filename, TString foldername, Int_t event, Int_t particle);

	Int_t NumberOfTracks(TString filename, TString foldername, Int_t event);
	Int_t NumberOfHits(TString filename, TString foldername, Int_t event, Int_t track);
	AliMUONHit* Hit(TString filename, TString foldername, Int_t event, Int_t track, Int_t hit);

	Int_t NumberOfSDigits(TString filename, TString foldername, Int_t event, Int_t chamber, Int_t cathode);
	AliMUONDigit* SDigit(TString filename, TString foldername, Int_t event, Int_t chamber, Int_t cathode, Int_t sdigit);

	Int_t NumberOfDigits(TString filename, TString foldername, Int_t event, Int_t chamber, Int_t cathode);
	AliMUONDigit* Digit(TString filename, TString foldername, Int_t event, Int_t chamber, Int_t cathode, Int_t digit);

	Int_t NumberOfRawClusters(TString filename, TString foldername, Int_t event, Int_t chamber);
	AliMUONRawCluster* RawCluster(TString filename, TString foldername, Int_t event, Int_t chamber, Int_t cluster);

	Int_t NumberOfLocalTriggers(TString filename, TString foldername, Int_t event);
	AliMUONLocalTrigger* LocalTrigger(TString filename, TString foldername, Int_t event, Int_t trigger);
	
	Bool_t SetFile(TString filename = "galice.root", TString foldername = "MUONFolder");
	Bool_t GetEvent(Int_t event = 0);

	Int_t NumberOfEvents();
	
	Int_t NumberOfParticles();
	TParticle* Particle(Int_t particle);

	Int_t NumberOfTracks();
	Int_t NumberOfHits(Int_t track);
	AliMUONHit* Hit(Int_t track, Int_t hit);

	Int_t NumberOfSDigits(Int_t chamber, Int_t cathode);
	AliMUONDigit* SDigit(Int_t chamber, Int_t cathode, Int_t sdigit);

	Int_t NumberOfDigits(Int_t chamber, Int_t cathode);
	AliMUONDigit* Digit(Int_t chamber, Int_t cathode, Int_t digit);

	Int_t NumberOfRawClusters(Int_t chamber);
	AliMUONRawCluster* RawCluster(Int_t chamber, Int_t cluster);

	Int_t NumberOfLocalTriggers();
	AliMUONLocalTrigger* LocalTrigger(Int_t trigger);
	
	
	// Returns the name of the currently selected file.
	TString CurrentFile()     { return fFilename;    };
	
	// Returns the name of the currently selected folder.
	TString CurrentFolder()   { return fFoldername;  };
	
	// Returns the number of the currently selected event.
	Int_t   CurrentEvent()    { return fEventnumber; };
	
	// Returns the currently selected track.
	Int_t   CurrentTrack()    { return fTrack;       };
	
	// Returns the currently selected cathode in TreeS.
	Int_t   CurrentSCathode() { return fSCathode;    };
	
	// Returns the currently selected cathode in TreeD.
	Int_t   CurrentDCathode() { return fCathode;     };

private:

	Bool_t LoadLoaders(TString filename, TString foldername);
	Bool_t FetchLoaders(TString filename, TString foldername);
	Bool_t FetchEvent(Int_t event);
	Bool_t FetchTreeK();
	Bool_t FetchTreeH();
	Bool_t FetchTreeS();
	Bool_t FetchTreeD();
	Bool_t FetchTreeR();
	
	
	Bool_t fHitAddressSet;     //! Flag specifying if the TTree address for the hit tree was set.
	Bool_t fSDigitAddressSet;  //! Flag specifying if the TTree address for the s-digit tree was set.
	Bool_t fDigitAddressSet;   //! Flag specifying if the TTree address for the digit tree was set.
	Bool_t fClusterAddressSet; //! Flag specifying if the TTree address for the cluster tree was set.
	Bool_t fTriggerAddressSet; //! Flag specifying if the TTree address for the trigger tree was set.

	AliRunLoader* fRunloader;  //! Pointer to the runloader object used.
	AliLoader* fMuonloader;    //! Pointer to the muon loader object used.
	AliMUONData fData;         //! Pointer to the muon raw data interface.
	TString fFilename;         //! The file name from which we are fetching data.
	TString fFoldername;       //! The folder name from which we are fetching data.
	Int_t fEventnumber;        //! The currently selected event.
	Int_t fTrack;              //! The currently selected track.
	Int_t fSCathode;           //! The currently selected cathode in TreeS.
	Int_t fCathode;            //! The currently selected cathode in TreeD.

	ClassDef(AliMUONDataInterface, 0)  // A easy to use interface to data in the MUON module.
};


#endif // ALI_MUON_DATA_INTERFACE_H
