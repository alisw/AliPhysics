/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
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

#include <cassert>
#include <TError.h>
#include <TParticle.h>

#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliMUONDataInterface.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONHit.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliLog.h"

#include <iostream>
using std::endl;
using std::cout;

///
/// \class AliMUONDataInterface
///
/// An easy to use interface to the MUON module data stored in
/// TreeK, TreeH, TreeS, TreeD and TreeR
/// One can fetch any of the data objects with all the calls to runloader, 
/// muon loader and AliMUONData done behind the scenes and automatically.
///
/// This interface in not necessarily the fastest way to fetch the data but
/// it is the easiest.
/// Note: If independant calls to the run loader, muon loader or 
/// AliMUONData objects are interspersed with calls to the 
/// AliMUONDataInterface to fetch data, one might need to call the Reset 
/// method between these method calls at some point to prevent 
/// AliMUONDataInterface from getting confused.
/// This is necessary since this object assumes the state of runloader,
/// muon loader nor AliMUONData has not changed between calls.
/// If the state has changes then one must call Reset so that 
/// AliMUONDataInterface refreshes what it knows about the state
/// of the loader and AliMUONData objects.
///
/// \deprecated We have to revisit all this AliMUONData stuff anyway,
/// and probably make a real AliMUONLoader instead...
///
/// \author Artur Szostak
/// email: artur@alice.phy.uct.ac.za


/// \cond CLASSIMP
ClassImp(AliMUONDataInterface)
/// \endcond

AliMUONDataInterface::AliMUONDataInterface()
	: TObject(), 
	  fCreatedRunLoader(kFALSE),
	  fCreatedRunLoaderSim(kFALSE),
	  fHitAddressSet(kFALSE),
	  fSDigitAddressSet(kFALSE),
	  fDigitAddressSet(kFALSE),
	  fClusterAddressSet(kFALSE),
	  fTriggerAddressSet(kFALSE),
	  fRecTracksAddressSet(kFALSE),
	  fRunloader(0x0),
	  fRunloaderSim(0x0),
	  fRecLoader(0x0),
	  fSimLoader(0x0),
	  fRecData(0x0, "MUON", "MUON"),
	  fSimData(0x0, "MUON", "MUON"),
	  fFilename(),
	  fFoldername("MUONLoader"),
	  fFoldernameSim("MUONLoaderSim"),
	  fEventnumber(-1),
	  fTrack(-1),
	  fSCathode(-1),
	  fCathode(-1)
{
/// Set all internal pointers to 0x0 and indices to -1.

	Reset();
}

AliMUONDataInterface::~AliMUONDataInterface()
{
/// Delete the runloader if we created it.
/// If the runloader is not to be deleted then call Reset just before 
/// the destructor is called.

	if (fRunloader != NULL && fCreatedRunLoader)
		delete fRunloader;

	if (fRunloaderSim != NULL && fCreatedRunLoaderSim)
		delete fRunloaderSim;
}

void AliMUONDataInterface::Reset()
{
/// Sets all internal pointers to NULL and indices to -1.
/// Note: No resources are released!
/// Specificaly AliRunLoader is not deleted.

	fCreatedRunLoader = kFALSE;
	fCreatedRunLoaderSim = kFALSE;
	fRunloader = NULL;
	fRunloaderSim = NULL;
	fRecLoader = NULL;
	fSimLoader = NULL;
	fEventnumber = -1;
	fTrack = -1;
	fSCathode = -1;
	fCathode = -1;
	fHitAddressSet = kFALSE;
	fSDigitAddressSet = kFALSE;
	fDigitAddressSet = kFALSE;
	fClusterAddressSet = kFALSE;
	fTriggerAddressSet = kFALSE;
	fRecTracksAddressSet = kFALSE;
}


Bool_t AliMUONDataInterface::UseCurrentRunLoader()
{
/// Tries to fetch the current runloader with AliRunLoader::GetRunLoader. If nothing is
/// currently loaded then kFALSE is returned and AliMUONDataInterface is reset.

	Reset();
	fRunloader = AliRunLoader::GetRunLoader();
	if (fRunloader == NULL) return kFALSE;
	// Fetch the current file name, folder name and event number.
	fFilename = fRunloader->GetFileName();
        // fFoldername = fRunloader->GetEventFolder()->GetName();
	fEventnumber = fRunloader->GetEventNumber();

	if ( ! FetchMuonLoader(fFilename.Data()) )
	{
		Reset();
		return kFALSE;
	}		

	return kTRUE;
}


Bool_t AliMUONDataInterface::FetchMuonLoader(TString filename)
{
/// Fetches the muon loader for the given filename/foldername


	fRecLoader = fRunloader->GetLoader("MUONLoader");
	if (fRecLoader == NULL)
	{
		AliError(Form("Could not find the MUON loader in file: %s and folder: %s", 
			(const char*)filename, fFoldername.Data() ));
		return kFALSE;
	}
	
	// Need to connect the muon loader to the AliMUONData object,
	// else class to fRecData will return NULL.
	fRecData.SetLoader(fRecLoader);

	fSimLoader = fRunloaderSim->GetLoader("MUONLoader");
	if (fSimLoader == NULL)
	{
		AliError(Form("Could not find the MUON loader in file: %s and folder: %s", 
			(const char*)filename, fFoldernameSim.Data()));
		return kFALSE;
	}
	
	// Need to connect the muon loader to the AliMUONData object,
	// else class to fSimData will return NULL.
	fSimData.SetLoader(fSimLoader);
	return kTRUE;
}


Bool_t AliMUONDataInterface::LoadLoaders(TString filename)
{
/// Load the run and muon loaders from the specified file and folder.
/// kTRUE is returned on success and kFALSE on failure.

	fRunloader = AliRunLoader::Open(filename, "MUONFolder", "READ");
	if (fRunloader == NULL)
	{
		AliError(Form("Could not find or load the run loader for the file: %s and folder: MUONFolder", 
			(const char*)filename));
		return kFALSE;
	}
	fCreatedRunLoader = kTRUE;

	fRunloaderSim = AliRunLoader::Open(filename, "MUONFolderSim", "READ");
	if (fRunloaderSim == NULL)
	{
		AliError(Form("Could not find or load the run loader for the file: %s and folder: MUONFolderSim", 
			(const char*)filename));
		return kFALSE;
	}
	fCreatedRunLoaderSim = kTRUE;

	if ( ! FetchMuonLoader(filename) )
	{
		fRunloader = NULL;
		fRunloaderSim = NULL;
		return kFALSE;
	}
	
	fFilename = filename;
	fEventnumber = -1;  // Reset the event number to force the event to be loaded.
	return kTRUE;
}


Bool_t AliMUONDataInterface::FetchLoaders(TString filename)
{
/// Fetch the run loader and muon loader objects from memory if they already exist,
/// or from memory if they do not. 
/// If the currently loaded run loader (if any) is not refering to the file and folder
/// we are interested in then it is deleted and reopened with the required file and
/// folder names.

	if (fRunloader == NULL)
	{
		fRunloader = AliRunLoader::GetRunLoader();
		if (fRunloader == NULL)
			return LoadLoaders(filename);
		else
		{
			if (fRecLoader == NULL)
			{
				if ( ! FetchMuonLoader(filename) )
				{
					fRunloader = NULL;
					return kFALSE;
				}
			}
		}
		
		// Fetch the current file and folder names.
		fFilename = fRunloader->GetFileName();
		// fFoldername = fRunloader->GetEventFolder()->GetName();
	}

	// If filename or foldername are not the same as the ones currently selected then
	// reopen the file.
	if ( filename.CompareTo(fFilename) != 0 )
	{
		delete fRunloader;
		return LoadLoaders(filename);
	}
	return kTRUE;
}


Bool_t AliMUONDataInterface::FetchEvent(Int_t event)
{
/// Fetch the specified event from the runloader and reset all the track, cathode
/// and address flags to force them to be reloaded.
/// If a negative event number is specified then the current runloader event
/// number is used.

	if (fEventnumber < 0)
	{
		fEventnumber = fRunloader->GetEventNumber();
		fTrack = -1;
		fSCathode = -1;
		fCathode = -1;
		fHitAddressSet = kFALSE;
		fSDigitAddressSet = kFALSE;
		fDigitAddressSet = kFALSE;
		fClusterAddressSet = kFALSE;
		fTriggerAddressSet = kFALSE;
		fRecTracksAddressSet = kFALSE;
	}
	if ( event != fEventnumber )
	{
		if ( fRunloader->GetEvent(event) < 0 ) return kFALSE;
		fEventnumber = event;
		fTrack = -1;
		fSCathode = -1;
		fCathode = -1;
		fHitAddressSet = kFALSE;
		fSDigitAddressSet = kFALSE;
		fDigitAddressSet = kFALSE;
		fClusterAddressSet = kFALSE;
		fTriggerAddressSet = kFALSE;
		fRecTracksAddressSet = kFALSE;
	}
	return kTRUE;
}


Bool_t AliMUONDataInterface::FetchTreeK()
{
/// Fetch the Kine tree from the current run loader.

	if (fRunloaderSim->TreeK() == NULL)
	{
		fRunloaderSim->LoadKinematics("READ");
		if (fRunloaderSim->TreeK() == NULL)
		{
			AliError("Could not load TreeK.");
			return kFALSE;
		}
	}
	return kTRUE;
}


Bool_t AliMUONDataInterface::FetchTreeH()
{
/// Fetch the Hits tree from the current muon loader.
/// Set all the required addresses etc...

	if (fSimLoader->TreeH() == NULL)
	{
		fSimLoader->LoadHits("READ");
		if (fSimLoader->TreeH() == NULL)
		{
			AliError("Could not load TreeH.");
			return kFALSE;
		}
		fSimData.SetTreeAddress("H");
		fHitAddressSet = kTRUE;
	}
	else if ( ! fHitAddressSet )
	{
		fSimData.SetTreeAddress("H");
		fHitAddressSet = kTRUE;
	}
	return kTRUE;
}


Bool_t AliMUONDataInterface::FetchTreeS()
{
/// Fetch the S-Digits tree from the current muon loader.
/// Set all the required addresses etc...

	if (fSimLoader->TreeS() == NULL)
	{
		fSimLoader->LoadSDigits("READ");
		if (fSimLoader->TreeS() == NULL)
		{
			AliError("Could not load TreeS.");
			return kFALSE;
		}
		fSimData.SetTreeAddress("S");
		fSDigitAddressSet = kTRUE;
	}
	else if ( ! fSDigitAddressSet )
	{
		fSimData.SetTreeAddress("S");
		fSDigitAddressSet = kTRUE;
	}
	return kTRUE;
}


Bool_t AliMUONDataInterface::FetchTreeD()
{
/// Fetch the digits tree from the current muon loader.
/// Set all the required addresses etc...

	if (fSimLoader->TreeD() == NULL)
	{
		fSimLoader->LoadDigits("READ");
		if (fSimLoader->TreeD() == NULL)
		{
			AliError("Could not load TreeD.");
			return kFALSE;
		}
		fSimData.SetTreeAddress("D");
		fDigitAddressSet = kTRUE;
	}
	else if ( ! fDigitAddressSet )
	{
		fSimData.SetTreeAddress("D");
		fDigitAddressSet = kTRUE;
	}
	return kTRUE;
}


Bool_t AliMUONDataInterface::FetchTreeR()
{
/// Fetch the reconstructed objects tree from the current muon loader.
/// Note: The addresses must still be set. 
  
  if (fRecLoader->TreeR() == NULL)
    {
      fRecLoader->LoadRecPoints("READ");
      if (fRecLoader->TreeR() == NULL)
	{
	  AliError("Could not load TreeR.");
	  return kFALSE;
	}
      
      // Need to reset these flags so that the cluster and trigger address
      // gets reset after this method. 
      fClusterAddressSet = kFALSE;
      fTriggerAddressSet = kFALSE;
    }
  return kTRUE;
}

Bool_t AliMUONDataInterface::FetchTreeT()
{
/// fetch the reconstructed tracks tree from the current muon loader
/// note : the addresses must still be set.
  if (fRecLoader->TreeT() == NULL)
    {
      fRecLoader->LoadTracks("READ");
      if (fRecLoader->TreeT() == NULL)
	{
	  AliError("Could not load TreeT.");
	  return kFALSE;
	}
      
      // Need to reset these flags so that the rec tracks address
      // gets reset after this method. 
      fRecTracksAddressSet = kFALSE;
    }
  return kTRUE;
}
  
Int_t AliMUONDataInterface::NumberOfEvents(TString filename, TString foldername)
{
/// Returns the number of events in the specified file/folder, and -1 on error.

	if ( ! FetchLoaders(filename) ) return -1;
	return fRunloader->GetNumberOfEvents();
}


Int_t AliMUONDataInterface::NumberOfParticles(TString filename, TString foldername, Int_t event)
{
/// Returns the number of events in the specified file/folder, and -1 on error.

	if ( ! FetchLoaders(filename) ) return -1;
	if ( ! FetchEvent(event) ) return -1;
	if ( ! FetchTreeK() ) return -1;
	return (Int_t) fRunloader->TreeK()->GetEntriesFast();
}


TParticle* AliMUONDataInterface::Particle(
		TString filename, TString foldername, Int_t event, Int_t particle
	)
{
/// Returns the specified particle in the given file, folder and event.
/// NULL is returned on error.

	if ( ! FetchLoaders(filename) ) return NULL;
	if ( ! FetchEvent(event) ) return NULL;
	if ( ! FetchTreeK() ) return NULL;
	
	TTree* treeK = fRunloader->TreeK();
	TParticle* p = NULL;
	treeK->GetBranch("Particles")->SetAddress(&p);
	treeK->GetEvent(particle);
	return p;
}


Int_t AliMUONDataInterface::NumberOfTracks(TString filename, TString foldername, Int_t event)
{
/// Returns the number of tracks in the specified file/folder and event.
/// -1 is returned on error.

	if ( ! FetchLoaders(filename) ) return -1;
	if ( ! FetchEvent(event) ) return -1;
	if ( ! FetchTreeH() ) return -1;
	return fSimData.GetNtracks();
}


Int_t AliMUONDataInterface::NumberOfHits(
		TString filename, TString foldername, Int_t event, Int_t track
	)
{
/// Returns the number of hits in the specified file/folder, event and track.
/// -1 is returned on error.

	if ( ! FetchLoaders(filename) ) return -1;
	if ( ! FetchEvent(event) ) return -1;
	if ( ! FetchTreeH() ) return -1;

	if (fTrack < 0 || fTrack != track)
	{
		fSimData.ResetHits();
		fSimData.GetTrack(track);
		fTrack = track;
	}
	return fSimData.Hits()->GetEntriesFast();
}


AliMUONHit* AliMUONDataInterface::Hit(
		TString filename, TString foldername, Int_t event,
		Int_t track, Int_t hit
	)
{
/// Returns the specified hit in the given file, folder, event and track.
/// NULL is returned on error.

	if ( ! FetchLoaders(filename) ) return NULL;
	if ( ! FetchEvent(event) ) return NULL;
	if ( ! FetchTreeH() ) return NULL;

	if (fTrack < 0 || fTrack != track)
	{
		fSimData.ResetHits();
		fSimData.GetTrack(track);
		fTrack = track;
	}
	return static_cast<AliMUONHit*>( fSimData.Hits()->At(hit) );
}


Int_t AliMUONDataInterface::NumberOfSDigits(
		TString filename, TString foldername, Int_t event,
		Int_t chamber, Int_t cathode
	)
{
/// Returns the number of s-digits in the given file, folder, event,
/// chamber and cathode. -1 is returned on error.

	assert( 0 <= chamber && chamber <= 13 );
	assert( 0 <= cathode && cathode <= 1 );
	
	if ( ! FetchLoaders(filename) ) return -1;
	if ( ! FetchEvent(event) ) return -1;
	if ( ! FetchTreeS() ) return -1;

	if ( fSCathode != cathode )
	{
		fSimData.ResetSDigits();
		fSimData.GetSDigits();
		fSCathode = cathode;
	}
	return fSimData.SDigits(chamber)->GetEntriesFast();
}


AliMUONDigit* AliMUONDataInterface::SDigit(
		TString filename, TString foldername, Int_t event,
		Int_t chamber, Int_t cathode, Int_t sdigit
	)
{
/// Returns the specified s-digit in the given file, folder, event,
/// chamber and cathode. NULL is returned on error.

	assert( 0 <= chamber && chamber <= 13 );
	assert( 0 <= cathode && cathode <= 1 );
	
	if ( ! FetchLoaders(filename) ) return NULL;
	if ( ! FetchEvent(event) ) return NULL;
	if ( ! FetchTreeS() ) return NULL;

	if ( fSCathode != cathode )
	{
		fSimData.ResetSDigits();
		fSimData.GetSDigits();
		fSCathode = cathode;
	}
	return static_cast<AliMUONDigit*>( fSimData.SDigits(chamber)->At(sdigit) );
}


Int_t AliMUONDataInterface::NumberOfDigits(
		TString filename, TString foldername, Int_t event,
		Int_t chamber, Int_t cathode
	)
{
/// Returns the number of digits in the given file, folder, event,
/// chamber and cathode. -1 is returned on error.
	assert( 0 <= chamber && chamber <= 13 );
	assert( 0 <= cathode && cathode <= 1 );
	
	if ( ! FetchLoaders(filename) ) return -1;
	if ( ! FetchEvent(event) ) return -1;
	if ( ! FetchTreeD() ) return -1;

	if ( fCathode != cathode )
	{
		fSimData.ResetDigits();
		fSimData.GetDigits();
		fCathode = cathode;
	}
	return fSimData.Digits(chamber)->GetEntriesFast();
}


AliMUONDigit* AliMUONDataInterface::Digit(
		TString filename, TString foldername, Int_t event,
		Int_t chamber, Int_t cathode, Int_t digit
	)
{
/// Returns the specified digit in the given file, folder, event,
/// chamber and cathode. NULL is returned on error.

	assert( 0 <= chamber && chamber <= 13 );
	assert( 0 <= cathode && cathode <= 1 );
	
	if ( ! FetchLoaders(filename) ) return NULL;
	if ( ! FetchEvent(event) ) return NULL;
	if ( ! FetchTreeD() ) return NULL;

	if ( fCathode != cathode )
	{
		fSimData.ResetDigits();
		fSimData.GetDigits();
		fCathode = cathode;
	}
	return static_cast<AliMUONDigit*>( fSimData.Digits(chamber)->At(digit) );
}


Int_t AliMUONDataInterface::NumberOfRawClusters(
		TString filename, TString foldername, Int_t event, Int_t chamber
	)
{
/// Returns the number of raw clusters in the specified file, folder, event and chamber.
/// -1 is returned or error.

	assert( 0 <= chamber && chamber <= 13 );
	if ( ! FetchLoaders(filename) ) return -1;
	if ( ! FetchEvent(event) ) return -1;
	if ( ! FetchTreeR() ) return -1;
	if ( ! fClusterAddressSet )
	{
		// If the raw cluster address in TreeR is not set yet then set it now.
		fRecData.SetTreeAddress("RC");
		fRecData.ResetRawClusters();
		fRecData.GetRawClusters();
		fClusterAddressSet = kTRUE;
	}
	return fRecData.RawClusters(chamber)->GetEntriesFast();
}


AliMUONRawCluster* AliMUONDataInterface::RawCluster(
		TString filename, TString foldername, Int_t event,
		Int_t chamber, Int_t cluster
	)
{
/// Fetch the specified raw cluster from the given file, folder, event and chamber number.
/// NULL is returned on error.

	assert( 0 <= chamber && chamber <= 13 );
	if ( ! FetchLoaders(filename) ) return NULL;
	if ( ! FetchEvent(event) ) return NULL;
	if ( ! FetchTreeR() ) return NULL;
	if ( ! fClusterAddressSet )
	{
		// If the raw cluster address in TreeR is not set yet then set it now.
		fRecData.SetTreeAddress("RC");
		fRecData.ResetRawClusters();
		fRecData.GetRawClusters();
		fClusterAddressSet = kTRUE;
	}
	return static_cast<AliMUONRawCluster*>( fRecData.RawClusters(chamber)->At(cluster) );
}


Int_t AliMUONDataInterface::NumberOfLocalTriggers(TString filename, TString foldername, Int_t event)
{
/// Return the number of local trigger objects in the specified file, folder and
/// event number. -1 is returned on error.

	if ( ! FetchLoaders(filename) ) return -1;
	if ( ! FetchEvent(event) ) return -1;
	if ( ! FetchTreeD() ) return -1;
	if ( ! fTriggerAddressSet )
	{
		// If the local trigger address in TreeR is not set yet then set it now.
		fRecData.SetTreeAddress("GLT");
		fRecData.ResetTrigger();
		fRecData.GetTriggerD();
		fTriggerAddressSet = kTRUE;
	}
	return fRecData.LocalTrigger()->GetEntriesFast();
}


AliMUONLocalTrigger* AliMUONDataInterface::LocalTrigger(
		TString filename, TString foldername, Int_t event, Int_t trigger
	)
{
/// Fetch the specified local trigger object from the given file, folder and event number.
/// NULL is returned on error.

	if ( ! FetchLoaders(filename) ) return NULL;
	if ( ! FetchEvent(event) ) return NULL;
	if ( ! FetchTreeD() ) return NULL;
	if ( ! fTriggerAddressSet )
	{
		// If the local trigger address in TreeR is not set yet then set it now.
		fRecData.SetTreeAddress("GLT");
		fRecData.ResetTrigger();
		fRecData.GetTriggerD();
		fTriggerAddressSet = kTRUE;
	}
	return static_cast<AliMUONLocalTrigger*>( fRecData.LocalTrigger()->At(trigger) );
}

Bool_t AliMUONDataInterface::SetFile(TString filename, TString foldername)
{
/// Set the current file and folder from which to fetch data.
/// kTRUE is returned if the run and muon loaders were found, else kFALSE. 

	return FetchLoaders(filename);
}


Bool_t AliMUONDataInterface::GetEvent(Int_t event)
{
/// Select the current event from which to fetch data.
/// kTRUE is returned if the event was found, else kFALSE is returned.

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return kFALSE;
	}
	else
		return FetchEvent(event);
}


Int_t AliMUONDataInterface::NumberOfEvents()
{
/// Get the number of events in the currently selected file.
/// -1 is returned on error.

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return -1;
	}
	return fRunloader->GetNumberOfEvents();
}


Int_t AliMUONDataInterface::NumberOfParticles()
{
/// Get the number of particles in the current event.
/// -1 is returned on error.

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return -1;
	}
	if ( ! FetchTreeK() ) return -1;
	return (Int_t) fRunloader->TreeK()->GetEntriesFast();
}


TParticle* AliMUONDataInterface::Particle(Int_t particle)
{
/// Fetch the specified particle from the current event.
/// NULL is returned on error.

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return NULL;
	}
	if (fEventnumber < 0)
	{
		AliError("Event not chosen.");
		return NULL;
	}
	if ( ! FetchTreeK() ) return NULL;
	TTree* treeK = fRunloader->TreeK();
	TParticle* p = NULL;
	treeK->GetBranch("Particles")->SetAddress(&p);
	treeK->GetEvent(particle);
	return p;
}


Int_t AliMUONDataInterface::NumberOfTracks()
{
/// Get the number of tracks in the current event.
/// -1 is returned on error.

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return -1;
	}
	if (fEventnumber < 0)
	{
		AliError( "Event not chosen.");
		return -1;
	}
	if ( ! FetchTreeH() ) return -1;
	return fSimData.GetNtracks();
}


Int_t AliMUONDataInterface::NumberOfHits(Int_t track)
{
/// Get the number of hits for the given track in the current event.
/// -1 is returned on error.

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return -1;
	}
	if (fEventnumber < 0)
	{
		AliError("Event not chosen.");
		return -1;
	}
	if ( ! FetchTreeH() ) return -1;
	if (fTrack < 0 || fTrack != track)
	{
		fSimData.ResetHits();
		fSimData.GetTrack(track);
		fTrack = track;
	}
	return fSimData.Hits()->GetEntriesFast();
}


AliMUONHit* AliMUONDataInterface::Hit(Int_t track, Int_t hit)
{
/// Fetch the specified hit from the current event.
/// NULL is returned on error.

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return NULL;
	}
	if (fEventnumber < 0)
	{
		AliError("Event not chosen.");
		return NULL;
	}
	if ( ! FetchTreeH() ) return NULL;
	if (fTrack < 0 || fTrack != track)
	{
		fSimData.ResetHits();
		fSimData.GetTrack(track);
		fTrack = track;
	}
	return static_cast<AliMUONHit*>( fSimData.Hits()->At(hit) );
}


Int_t AliMUONDataInterface::NumberOfSDigits(Int_t chamber, Int_t cathode)
{
/// Get the number of s-digits on the chamber, cathode in the current event.
/// -1 is returned on error.

	assert( 0 <= chamber && chamber <= 13 );
	assert( 0 <= cathode && cathode <= 1 );
	
	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return -1;
	}
	if (fEventnumber < 0)
	{
		AliError("Event not chosen.");
		return -1;
	}

	if ( ! FetchTreeS() ) return -1;
	if ( fSCathode != cathode )
	{
		fSimData.ResetSDigits();
		fSimData.GetSDigits();
		fSCathode = cathode;
	}
	return fSimData.SDigits(chamber)->GetEntriesFast();
}


AliMUONDigit* AliMUONDataInterface::SDigit(Int_t chamber, Int_t cathode, Int_t sdigit)
{
/// Fetch the specified s-digits on the chamber, cathode from the current event.
/// NULL is returned on error.

	assert( 0 <= chamber && chamber <= 13 );
	assert( 0 <= cathode && cathode <= 1 );
	
	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return NULL;
	}
	if (fEventnumber < 0)
	{
		AliError("Event not chosen.");
		return NULL;
	}

	if ( ! FetchTreeS() ) return NULL;
	if ( fSCathode != cathode )
	{
		fSimData.ResetSDigits();
		fSimData.GetSDigits();
		fSCathode = cathode;
	}
	return static_cast<AliMUONDigit*>( fSimData.SDigits(chamber)->At(sdigit) );
}


Int_t AliMUONDataInterface::NumberOfDigits(Int_t chamber, Int_t cathode)
{
/// Get the number of digits on the chamber, cathode in the current event.
/// -1 is returned on error.

	assert( 0 <= chamber && chamber <= 13 );
	assert( 0 <= cathode && cathode <= 1 );
	
	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return -1;
	}
	if (fEventnumber < 0)
	{
		AliError("Event not chosen.");
		return -1;
	}
	
	if ( ! FetchTreeD() ) return -1;
	if ( fCathode != cathode )
	{
		fSimData.ResetDigits();
		fSimData.GetDigits();
		fCathode = cathode;
	}
	return fSimData.Digits(chamber)->GetEntriesFast();
}


AliMUONDigit* AliMUONDataInterface::Digit(Int_t chamber, Int_t cathode, Int_t digit)
{
/// Fetch the specified digits on the chamber, cathode from the current event.
/// NULL is returned on error.

	assert( 0 <= chamber && chamber <= 13 );
	assert( 0 <= cathode && cathode <= 1 );
	
	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return NULL;
	}
	if (fEventnumber < 0)
	{
		AliError("Event not chosen.");
		return NULL;
	}

	if ( ! FetchTreeD() ) return NULL;
	if ( fCathode != cathode )
	{
		fSimData.ResetDigits();
		fSimData.GetDigits();
		fCathode = cathode;
	}
	return static_cast<AliMUONDigit*>( fSimData.Digits(chamber)->At(digit) );
}


Int_t AliMUONDataInterface::NumberOfRawClusters(Int_t chamber)
{
/// Get the number of raw clusters on the given chamber in the current event.
/// -1 is returned on error.

	assert( 0 <= chamber && chamber <= 13 );

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return -1;
	}
	if (fEventnumber < 0)
	{
		AliError("Event not chosen.");
		return -1;
	}

	if ( ! FetchTreeR() ) return -1;
	if ( ! fClusterAddressSet )
	{
		fRecData.SetTreeAddress("RC");
		fRecData.ResetRawClusters();
		fRecData.GetRawClusters();
		fClusterAddressSet = kTRUE;
	}
	return fRecData.RawClusters(chamber)->GetEntriesFast();
}


AliMUONRawCluster* AliMUONDataInterface::RawCluster(Int_t chamber, Int_t cluster)
{
/// Fetch the specified raw cluster on the given chamber from the current event.
/// NULL is returned on error.

	assert( 0 <= chamber && chamber <= 13 );

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return NULL;
	}
	if (fEventnumber < 0)
	{
		AliError("Event not chosen.");
		return NULL;
	}

	if ( ! FetchTreeR() ) return NULL;
	if ( ! fClusterAddressSet )
	{
		fRecData.SetTreeAddress("RC");
		fRecData.ResetRawClusters();
		fRecData.GetRawClusters();
		fClusterAddressSet = kTRUE;
	}
	return static_cast<AliMUONRawCluster*>( fRecData.RawClusters(chamber)->At(cluster) );
}


Int_t AliMUONDataInterface::NumberOfLocalTriggers()
{
/// Get the number of local trigger objects in the current event.
/// -1 is returned on error.

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return -1;
	}
	if (fEventnumber < 0)
	{
		AliError("Event not chosen.");
		return -1;
	}

	if ( ! FetchTreeD() ) return -1;
	if ( ! fTriggerAddressSet )
	{
		fRecData.SetTreeAddress("GLT");
		fRecData.ResetTrigger();
		fRecData.GetTriggerD();
		fTriggerAddressSet = kTRUE;
	}
	return fRecData.LocalTrigger()->GetEntriesFast();
}


AliMUONLocalTrigger* AliMUONDataInterface::LocalTrigger(Int_t trigger)
{
/// Fetch the specified local trigger object from the current event.
/// NULL is returned on error.

	if (fRunloader == NULL)
	{
		AliError("File not set.");
		return NULL;
	}
	if (fEventnumber < 0)
	{
		AliError( "Event not chosen.");
		return NULL;
	}

	if ( ! FetchTreeD() ) return NULL;
	if ( ! fTriggerAddressSet )
	{
		fRecData.SetTreeAddress("GLT");
		fRecData.ResetTrigger();
		fRecData.GetTriggerD();
		fTriggerAddressSet = kTRUE;
	}
	return static_cast<AliMUONLocalTrigger*>( fRecData.LocalTrigger()->At(trigger) );
}

Int_t AliMUONDataInterface::NumberOfGlobalTriggers()
{
/// Get the number of local trigger objects in the current event.
/// -1 is returned on error.
  
  if (fRunloader == NULL)
    {
      AliError("File not set.");
      return -1;
    }
  if (fEventnumber < 0)
    {
      AliError("Event not chosen.");
      return -1;
    }
  
  if ( ! FetchTreeD() ) return -1;
  if ( ! fTriggerAddressSet )
    {
      fRecData.SetTreeAddress("GLT");
      fRecData.ResetTrigger();
      fRecData.GetTriggerD();
      fTriggerAddressSet = kTRUE;
    }
  return fRecData.GlobalTrigger()->GetEntriesFast();
}

AliMUONGlobalTrigger* AliMUONDataInterface::GlobalTrigger(Int_t trigger)
{
/// Fetch the specified local trigger object from the current event.
/// NULL is returned on error.
  
  if (fRunloader == NULL)
    {
      AliError("File not set.");
      return NULL;
    }
  if (fEventnumber < 0)
    {
      AliError( "Event not chosen.");
      return NULL;
    }
  
  if ( ! FetchTreeD() ) return NULL;
  if ( ! fTriggerAddressSet )
    {
      fRecData.SetTreeAddress("GLT");
      fRecData.ResetTrigger();
      fRecData.GetTriggerD();
      fTriggerAddressSet = kTRUE;
    }
  return static_cast<AliMUONGlobalTrigger*>( fRecData.GlobalTrigger()->At(trigger) );
}

Int_t AliMUONDataInterface::NumberOfRecTracks()
{
/// Fetch the number of reconstructed tracks from the current event.
/// NULL is returned on error.
  
  if (fRunloader == NULL)
    {
      AliError("File not set.");
      return -1;
    }
  if (fEventnumber < 0)
    {
      AliError( "Event not chosen.");
      return -1;
    }
  
  if ( ! FetchTreeT() ) return -1;
  if ( ! fRecTracksAddressSet )
    {
      fRecData.SetTreeAddress("RT");
      fRecData.ResetRecTracks();
      fRecData.GetRecTracks();
      fRecTracksAddressSet = kTRUE;
    }
  return fRecData.RecTracks()->GetEntriesFast();
}

AliMUONTrack* AliMUONDataInterface::RecTrack(Int_t rectrack)
{
/// Fetch the specified reconstructed track object from the current event.
/// NULL is returned on error.
  
  if (fRunloader == NULL)
    {
      AliError("File not set.");
      return NULL;
    }
  if (fEventnumber < 0)
    {
      AliError( "Event not chosen.");
      return NULL;
    }
  
  if ( ! FetchTreeT() ) return NULL;
  if ( ! fRecTracksAddressSet )
    {
      fRecData.SetTreeAddress("RT");
      fRecData.ResetRecTracks();
      fRecData.GetRecTracks();
      fRecTracksAddressSet = kTRUE;
    }
  return static_cast<AliMUONTrack*>( fRecData.RecTracks()->At(rectrack) );
  // return (AliMUONTrack*)(fRecData.RecTracks()->At(rectrack));
}
