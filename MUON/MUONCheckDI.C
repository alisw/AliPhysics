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

//
// Macro for checking AliMUONDataInterface and AliMUONMCDataInterface.
// By Bruce Becker, DAPNIA/SPhN/CEA Saclay
//
// Modified to updated versions of data interfaces.
//  Artur Szostak <artursz@iafrica.com> (University of Cape Town)
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <Riostream.h>
#include <TObjArray.h>
#include <TIterator.h>
#include "AliMUONHit.h"
#include "AliMUONVDigit.h"
#include "AliMUONVCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONDataInterface.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMpConstants.h"
#include "AliMpDEManager.h"
#include <cstdlib>
#endif


/**
 * This routine implements a the comparison functionality which is missing for
 * classes like AliMUONTrack and the various AliMUONxxxTrigger classes.
 * A result of -1 is returned if a < b, 0 if a == b and +1 if a > b.
 */
Int_t Compare(const TObject* a, const TObject* b)
{
	int result = -999;
	if (a->IsA() == AliMUONTrack::Class() && b->IsA() == AliMUONTrack::Class())
	{
		const AliMUONTrack* ta = static_cast<const AliMUONTrack*>(a);
		const AliMUONTrack* tb = static_cast<const AliMUONTrack*>(b);
		if (ta->GetNClusters() < tb->GetNClusters()) return -1;
		if (ta->GetNClusters() > tb->GetNClusters()) return 1;
		if (ta->GetMatchTrigger() < tb->GetMatchTrigger()) return -1;
		if (ta->GetMatchTrigger() > tb->GetMatchTrigger()) return 1;
		if (ta->GetLoTrgNum() < tb->GetLoTrgNum()) return -1;
		if (ta->GetLoTrgNum() > tb->GetLoTrgNum()) return 1;
		if (ta->GetChi2MatchTrigger() < tb->GetChi2MatchTrigger()) return -1;
		if (ta->GetChi2MatchTrigger() > tb->GetChi2MatchTrigger()) return 1;
		const AliMUONTrackParam* tpa = static_cast<const AliMUONTrackParam*>(ta->GetTrackParamAtCluster()->First());
		const AliMUONTrackParam* tpb = static_cast<const AliMUONTrackParam*>(tb->GetTrackParamAtCluster()->First());
		if (tpa->GetNonBendingCoor() < tpb->GetNonBendingCoor()) return -1;
		if (tpa->GetNonBendingCoor() > tpb->GetNonBendingCoor()) return 1;
		if (tpa->GetNonBendingSlope() < tpb->GetNonBendingSlope()) return -1;
		if (tpa->GetNonBendingSlope() > tpb->GetNonBendingSlope()) return 1;
		if (tpa->GetBendingCoor() < tpb->GetBendingCoor()) return -1;
		if (tpa->GetBendingCoor() > tpb->GetBendingCoor()) return 1;
		if (tpa->GetBendingSlope() < tpb->GetBendingSlope()) return -1;
		if (tpa->GetBendingSlope() > tpb->GetBendingSlope()) return 1;
		if (tpa->GetInverseBendingMomentum() < tpb->GetInverseBendingMomentum()) return -1;
		if (tpa->GetInverseBendingMomentum() > tpb->GetInverseBendingMomentum()) return 1;
		if (tpa->GetCharge() < tpb->GetCharge()) return -1;
		if (tpa->GetCharge() > tpb->GetCharge()) return 1;
		return 0;
	}
	else if (a->IsA() == AliMUONLocalTrigger::Class() && b->IsA() == AliMUONLocalTrigger::Class())
	{
		result = memcmp(a, b, sizeof(AliMUONLocalTrigger));
	}
	else if (a->IsA() == AliMUONRegionalTrigger::Class() && b->IsA() == AliMUONRegionalTrigger::Class())
	{
		result = memcmp(a, b, sizeof(AliMUONRegionalTrigger));
	}
	else if (a->IsA() == AliMUONGlobalTrigger::Class() && b->IsA() == AliMUONGlobalTrigger::Class())
	{
		result = memcmp(a, b, sizeof(AliMUONGlobalTrigger));
	}
	else if (a->IsA() == AliMUONTriggerTrack::Class() && b->IsA() == AliMUONTriggerTrack::Class())
	{
		const AliMUONTriggerTrack* ta = static_cast<const AliMUONTriggerTrack*>(a);
		const AliMUONTriggerTrack* tb = static_cast<const AliMUONTriggerTrack*>(b);
		if (ta->GetX11() < tb->GetX11()) return -1;
		if (ta->GetX11() > tb->GetX11()) return 1;
		if (ta->GetY11() < tb->GetY11()) return -1;
		if (ta->GetY11() > tb->GetY11()) return 1;
		if (ta->GetThetax() < tb->GetThetax()) return -1;
		if (ta->GetThetax() > tb->GetThetax()) return 1;
		if (ta->GetThetay() < tb->GetThetay()) return -1;
		if (ta->GetThetay() > tb->GetThetay()) return 1;
		if (ta->GetLoTrgNum() < tb->GetLoTrgNum()) return -1;
		if (ta->GetLoTrgNum() > tb->GetLoTrgNum()) return 1;
		if (ta->GetGTPattern() < tb->GetGTPattern()) return -1;
		if (ta->GetGTPattern() > tb->GetGTPattern()) return 1;
		return 0;
	}
	else
	{
		result = memcmp(a, b, sizeof(TObject));
	}
	
	if (result < 0) return -1;
	if (result > 0) return 1;
	return 0;
}

/**
 * This method fills internal arrays with local and regional triggers returned
 * by AliMUONMCDataInterface. For each set of interface methods available in
 * AliMUONMCDataInterface a TObjArray is created for local and another for regional
 * triggers. These arrays are filled with copies of the trigger objects.
 * The global trigger object is also copied out using the 2 different methods.
 * The arrays and objects are then compared to each other. The arrays and objects
 * should contain the same information if everything is working correctly with
 * AliMUONMCDataInterface. If not then the difference is printed together with an
 * error message and kFALSE is returned.
 */
bool SimTriggersOk()
{
	AliMUONMCDataInterface data;
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray localsFromStore, regionalsFromStore;
		localsFromStore.SetOwner(kTRUE);
		regionalsFromStore.SetOwner(kTRUE);
		AliMUONVTriggerStore* store = data.TriggerStore(event);
		AliMUONGlobalTrigger* globalFromStore = static_cast<AliMUONGlobalTrigger*>(store->Global()->Clone());
		TIter nextLocal(store->CreateLocalIterator());
		AliMUONLocalTrigger* localTrig;
		while ( (localTrig = static_cast<AliMUONLocalTrigger*>( nextLocal() )) != NULL )
		{
			localsFromStore.Add(localTrig->Clone());
		}
		TIter nextRegional(store->CreateRegionalIterator());
		AliMUONRegionalTrigger* regionalTrig;
		while ( (regionalTrig = static_cast<AliMUONRegionalTrigger*>( nextRegional() )) != NULL )
		{
			regionalsFromStore.Add(regionalTrig->Clone());
		}
		
		TObjArray localsByIndex, regionalsByIndex;
		localsByIndex.SetOwner(kTRUE);
		regionalsByIndex.SetOwner(kTRUE);
		data.GetEvent(event);
		AliMUONGlobalTrigger* globalByMethod = static_cast<AliMUONGlobalTrigger*>(data.GlobalTrigger()->Clone());
		Int_t nlocals = data.NumberOfLocalTriggers();
		for (Int_t i = 0; i < nlocals; i++)
		{
			localTrig = data.LocalTrigger(i);
			localsByIndex.Add(localTrig->Clone());
		}
		Int_t nregionals = data.NumberOfRegionalTriggers();
		for (Int_t i = 0; i < nregionals; i++)
		{
			regionalTrig = data.RegionalTrigger(i);
			regionalsByIndex.Add(regionalTrig->Clone());
		}
		
		// Now check that all the lists of local, regional and global triggers
		// contain the same results.
		// They must. If they do not then something is wrong with the implementation
		// of AliMUONMCDataInterface.
		if (Compare(globalFromStore, globalByMethod) != 0)
		{
			Error(	"SimTriggersOk",
				"The AliMUONMCDataInterface does not return identical global"
				  " triggers through all its user interface methods."
			);
			globalFromStore->Print();
			globalByMethod->Print();
			return false;
		}
		delete globalFromStore;
		delete globalByMethod;
		if (localsFromStore.GetEntriesFast() != localsByIndex.GetEntriesFast())
		{
			Error(	"SimTriggersOk",
				"The AliMUONMCDataInterface does not return all the local triggers"
				  " correctly through all its user interface methods. We got the"
				  " following numbers of local triggers: %d and %d",
				localsFromStore.GetEntriesFast(),
				localsByIndex.GetEntriesFast()
			);
			return false;
		}
		if (regionalsFromStore.GetEntriesFast() != regionalsByIndex.GetEntriesFast())
		{
			Error(	"SimTriggersOk",
				"The AliMUONMCDataInterface does not return all the regional triggers"
				  " correctly through all its user interface methods. We got the"
				  " following numbers of regional triggers: %d and %d",
				regionalsFromStore.GetEntriesFast(),
				regionalsByIndex.GetEntriesFast()
			);
			return false;
		}
		for (Int_t i = 0; i < localsFromStore.GetEntriesFast(); i++)
		{
			if (Compare(localsFromStore[i], localsByIndex[i]) != 0)
			{
				Error(	"SimTriggersOk",
					"The AliMUONMCDataInterface does not return identical local"
					  " triggers through all its user interface methods. The"
					  " incorrect local trigger has index %d.",
					i
				);
				localsFromStore[i]->Print();
				localsByIndex[i]->Print();
				return false;
			}
		}
		for (Int_t i = 0; i < regionalsFromStore.GetEntriesFast(); i++)
		{
			if (Compare(regionalsFromStore[i], regionalsByIndex[i]) != 0)
			{
				Error(	"SimTriggersOk",
					"The AliMUONMCDataInterface does not return identical regional"
					  " triggers through all its user interface methods. The"
					  " incorrect regional trigger has index %d.",
					i
				);
				regionalsFromStore[i]->Print();
				regionalsByIndex[i]->Print();
				return false;
			}
		}
	}
	return true;
}

/**
 * This method fills internal arrays with digits returned by the AliMUONDataInterface.
 * For each set of interface methods available a TObjArray is filled with copies of 
 * the digits. These arrays are sorted and then compared to each other. The arrays
 * should contain the same digit information if everything is working correctly with
 * AliMUONDataInterface. If not then the difference is printed together with an
 * error message and kFALSE is returned.
 */
bool RecDigitsOk()
{
	AliMUONDataInterface data;
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray digitsFromStore;
		digitsFromStore.SetOwner(kTRUE);
		AliMUONVDigitStore* store = data.DigitStore(event);
		TIter next(store->CreateIterator());
		AliMUONVDigit* digit;
		while ( (digit = static_cast<AliMUONVDigit*>( next() )) != NULL )
		{
			digitsFromStore.Add(digit->Clone());
		}
		digitsFromStore.Sort();
		
		TObjArray digitsByDetElem;
		digitsByDetElem.SetOwner(kTRUE);
		data.GetEvent(event);
		for (Int_t detElem = 0; detElem < 1500; detElem++)
		{
			if (! AliMpDEManager::IsValidDetElemId(detElem)) continue;
			Int_t ndigits = data.NumberOfDigits(detElem);
			for (Int_t i = 0; i < ndigits; i++)
			{
				AliMUONVDigit* digit = data.Digit(detElem, i);
				digitsByDetElem.Add(digit->Clone());
			}
		}
		digitsByDetElem.Sort();
		
		TObjArray digitsByChamber;
		digitsByChamber.SetOwner(kTRUE);
		data.GetEvent(event);
		for (Int_t chamber = 0; chamber < AliMpConstants::NofChambers(); chamber++)
		for (Int_t cathode = 0; cathode < 2; cathode++)
		{
			Int_t ndigits = data.NumberOfDigits(chamber, cathode);
			for (Int_t i = 0; i < ndigits; i++)
			{
				AliMUONVDigit* digit = data.Digit(chamber, cathode, i);
				digitsByChamber.Add(digit->Clone());
			}
		}
		digitsByChamber.Sort();
		
		// Now check that all the lists of digits contain the same results.
		// They must. If they do not then something is wrong with the implementation
		// of AliMUONDataInterface.
		if (digitsFromStore.GetEntriesFast() != digitsByDetElem.GetEntriesFast()
		    || digitsFromStore.GetEntriesFast() != digitsByChamber.GetEntriesFast())
		{
			Error(	"RecDigitsOk",
				"The AliMUONDataInterface does not return all the digits correctly"
				  " through all its user interface methods. We got the following"
				  " numbers of digits: %d, %d and %d",
				digitsFromStore.GetEntriesFast(),
				digitsByDetElem.GetEntriesFast(),
				digitsByChamber.GetEntriesFast()
			);
			return false;
		}
		for (Int_t i = 0; i < digitsFromStore.GetEntriesFast(); i++)
		{
			if (digitsFromStore[i]->Compare(digitsByDetElem[i]) != 0
			    || digitsFromStore[i]->Compare(digitsByChamber[i]) != 0)
			{
				Error(	"RecDigitsOk",
					"The AliMUONDataInterface does not return identical digits"
					  " through all its user interface methods. The incorrect"
					  " digit has index %d after sorting.",
					i
				);
				digitsFromStore[i]->Print();
				digitsByChamber[i]->Print();
				digitsByDetElem[i]->Print();
				return false;
			}
		}
	}
	return true;
}

/**
 * This method fills internal arrays with raw clusters returned by AliMUONDataInterface.
 * For each set of interface methods available in AliMUONDataInterface a TObjArray is
 * filled with copies of the raw clusters. These arrays are sorted and then compared
 * to each other. The arrays should contain the same information if everything is
 * working correctly with AliMUONDataInterface. If not then the difference is printed
 * together with an error message and kFALSE is returned.
 */
bool RawClustersOk()
{
	AliMUONDataInterface data;
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray clustersFromStore;
		clustersFromStore.SetOwner(kTRUE);
		AliMUONVClusterStore* store = data.ClusterStore(event);
		TIter next(store->CreateIterator());
		AliMUONVCluster* cluster;
		while ( (cluster = static_cast<AliMUONVCluster*>( next() )) != NULL )
		{
			clustersFromStore.Add(cluster->Clone());
		}
		clustersFromStore.Sort();
		
		TObjArray clustersByChamber;
		clustersByChamber.SetOwner(kTRUE);
		data.GetEvent(event);
		for (Int_t chamber = 0; chamber < AliMpConstants::NofChambers(); chamber++)
		{
			Int_t nclusters = data.NumberOfRawClusters(chamber);
			for (Int_t i = 0; i < nclusters; i++)
			{
				AliMUONVCluster* cluster = data.RawCluster(chamber, i);
				clustersByChamber.Add(cluster->Clone());
			}
		}
		clustersByChamber.Sort();
		
		// Now check that all the lists of clusters contain the same results.
		// They must. If they do not then something is wrong with the implementation
		// of AliMUONDataInterface.
		if (clustersFromStore.GetEntriesFast() != clustersByChamber.GetEntriesFast())
		{
			Error(	"RawClustersOk",
				"The AliMUONDataInterface does not return all the clusters correctly"
				  " through all its user interface methods. We got the following"
				  " numbers of clusters: %d and %d",
				clustersFromStore.GetEntriesFast(),
				clustersByChamber.GetEntriesFast()
			);
			return false;
		}
		for (Int_t i = 0; i < clustersFromStore.GetEntriesFast(); i++)
		{
			if (clustersFromStore[i]->Compare(clustersByChamber[i]) != 0)
			{
				Error(	"RawClustersOk",
					"The AliMUONDataInterface does not return identical clusters"
					  " through all its user interface methods. The incorrect"
					  " cluster has index %d after sorting.",
					i
				);
				clustersFromStore[i]->Print();
				clustersByChamber[i]->Print();
				return false;
			}
		}
	}
	return true;
}

/**
 * This method fills internal arrays with tracks returned by AliMUONDataInterface.
 * For each set of interface methods available in AliMUONDataInterface a TObjArray is
 * filled with copies of the raw clusters. These arrays are then compared to each
 * other. The arrays should contain the same information if everything is working
 * correctly with AliMUONDataInterface. If not then the difference is printed
 * together with an error message and kFALSE is returned.
 */
bool TracksOk()
{
	AliMUONDataInterface data;
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray tracksFromStore;
		tracksFromStore.SetOwner(kTRUE);
		AliMUONVTrackStore* store = data.TrackStore(event);
		TIter next(store->CreateIterator());
		AliMUONTrack* track;
		while ( (track = static_cast<AliMUONTrack*>( next() )) != NULL )
		{
			tracksFromStore.Add(track->Clone());
		}
		
		TObjArray tracksByIndex;
		tracksByIndex.SetOwner(kTRUE);
		data.GetEvent(event);
		Int_t ntracks = data.NumberOfTracks();
		for (Int_t i = 0; i < ntracks; i++)
		{
			AliMUONTrack* track = data.Track(i);
			tracksByIndex.Add(track->Clone());
		}
		
		// Now check that all the lists of tracks contain the same results.
		// They must. If they do not then something is wrong with the implementation
		// of AliMUONDataInterface.
		if (tracksFromStore.GetEntriesFast() != tracksByIndex.GetEntriesFast())
		{
			Error(	"TracksOk",
				"The AliMUONDataInterface does not return all the tracks correctly"
				  " through all its user interface methods. We got the following"
				  " numbers of tracks: %d and %d",
				tracksFromStore.GetEntriesFast(),
				tracksByIndex.GetEntriesFast()
			);
			return false;
		}
		for (Int_t i = 0; i < tracksFromStore.GetEntriesFast(); i++)
		{
			if (Compare(tracksFromStore[i], tracksByIndex[i]) != 0)
			{
				Error(	"TracksOk",
					"The AliMUONDataInterface does not return identical tracks"
					  " through all its user interface methods. The incorrect"
					  " track has index %d.",
					i
				);
				tracksFromStore[i]->Print();
				tracksByIndex[i]->Print();
				return false;
			}
		}
	}
	return true;
}

/**
 * This method fills internal arrays with local and regional triggers returned
 * by AliMUONDataInterface. For each set of interface methods available in
 * AliMUONDataInterface a TObjArray is created for local and another for regional
 * triggers. These arrays are filled with copies of the trigger objects.
 * The global trigger object is also copied out using the 2 different methods.
 * The arrays and objects are then compared to each other. The arrays and objects
 * should contain the same information if everything is working correctly with
 * AliMUONDataInterface. If not then the difference is printed together with an
 * error message and kFALSE is returned.
 */
bool TriggersOk()
{
	AliMUONDataInterface data;
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray localsFromStore, regionalsFromStore;
		localsFromStore.SetOwner(kTRUE);
		regionalsFromStore.SetOwner(kTRUE);
		AliMUONVTriggerStore* store = data.TriggerStore(event);
		AliMUONGlobalTrigger* globalFromStore = static_cast<AliMUONGlobalTrigger*>(store->Global()->Clone());
		TIter nextLocal(store->CreateLocalIterator());
		AliMUONLocalTrigger* localTrig;
		while ( (localTrig = static_cast<AliMUONLocalTrigger*>( nextLocal() )) != NULL )
		{
			localsFromStore.Add(localTrig->Clone());
		}
		TIter nextRegional(store->CreateRegionalIterator());
		AliMUONRegionalTrigger* regionalTrig;
		while ( (regionalTrig = static_cast<AliMUONRegionalTrigger*>( nextRegional() )) != NULL )
		{
			regionalsFromStore.Add(regionalTrig->Clone());
		}
		
		TObjArray localsByIndex, regionalsByIndex;
		localsByIndex.SetOwner(kTRUE);
		regionalsByIndex.SetOwner(kTRUE);
		data.GetEvent(event);
		AliMUONGlobalTrigger* globalByMethod = static_cast<AliMUONGlobalTrigger*>(data.GlobalTrigger()->Clone());
		Int_t nlocals = data.NumberOfLocalTriggers();
		for (Int_t i = 0; i < nlocals; i++)
		{
			localTrig = data.LocalTrigger(i);
			localsByIndex.Add(localTrig->Clone());
		}
		Int_t nregionals = data.NumberOfRegionalTriggers();
		for (Int_t i = 0; i < nregionals; i++)
		{
			regionalTrig = data.RegionalTrigger(i);
			regionalsByIndex.Add(regionalTrig->Clone());
		}
		
		// Now check that all the lists of local, regional and global triggers
		// contain the same results.
		// They must. If they do not then something is wrong with the implementation
		// of AliMUONDataInterface.
		if (Compare(globalFromStore, globalByMethod) != 0)
		{
			Error(	"TriggersOk",
				"The AliMUONDataInterface does not return identical global"
				  " triggers through all its user interface methods."
			);
			globalFromStore->Print();
			globalByMethod->Print();
			return false;
		}
		delete globalFromStore;
		delete globalByMethod;
		if (localsFromStore.GetEntriesFast() != localsByIndex.GetEntriesFast())
		{
			Error(	"TriggersOk",
				"The AliMUONDataInterface does not return all the local triggers"
				  " correctly through all its user interface methods. We got the"
				  " following numbers of local triggers: %d and %d",
				localsFromStore.GetEntriesFast(),
				localsByIndex.GetEntriesFast()
			);
			return false;
		}
		if (regionalsFromStore.GetEntriesFast() != regionalsByIndex.GetEntriesFast())
		{
			Error(	"TriggersOk",
				"The AliMUONDataInterface does not return all the regional triggers"
				  " correctly through all its user interface methods. We got the"
				  " following numbers of regional triggers: %d and %d",
				regionalsFromStore.GetEntriesFast(),
				regionalsByIndex.GetEntriesFast()
			);
			return false;
		}
		for (Int_t i = 0; i < localsFromStore.GetEntriesFast(); i++)
		{
			if (Compare(localsFromStore[i], localsByIndex[i]) != 0)
			{
				Error(	"TriggersOk",
					"The AliMUONDataInterface does not return identical local"
					  " triggers through all its user interface methods. The"
					  " incorrect local trigger has index %d.",
					i
				);
				localsFromStore[i]->Print();
				localsByIndex[i]->Print();
				return false;
			}
		}
		for (Int_t i = 0; i < regionalsFromStore.GetEntriesFast(); i++)
		{
			if (Compare(regionalsFromStore[i], regionalsByIndex[i]) != 0)
			{
				Error(	"TriggersOk",
					"The AliMUONDataInterface does not return identical regional"
					  " triggers through all its user interface methods. The"
					  " incorrect regional trigger has index %d.",
					i
				);
				regionalsFromStore[i]->Print();
				regionalsByIndex[i]->Print();
				return false;
			}
		}
	}
	return true;
}

/**
 * This method fills internal arrays with trigger tracks returned by AliMUONDataInterface.
 * For each set of interface methods available in AliMUONDataInterface a TObjArray is
 * filled with copies of the trigger tracks. These arrays are then compared to each
 * other. The arrays should contain the same information if everything is working
 * correctly with AliMUONDataInterface. If not then the difference is printed
 * together with an error message and kFALSE is returned.
 */
bool TriggerTracksOk()
{
	AliMUONDataInterface data;
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray tracksFromStore;
		tracksFromStore.SetOwner(kTRUE);
		AliMUONVTriggerTrackStore* store = data.TriggerTrackStore(event);
		TIter next(store->CreateIterator());
		AliMUONTriggerTrack* track;
		while ( (track = static_cast<AliMUONTriggerTrack*>( next() )) != NULL )
		{
			tracksFromStore.Add(track->Clone());
		}
		
		TObjArray tracksByIndex;
		tracksByIndex.SetOwner(kTRUE);
		data.GetEvent(event);
		Int_t ntracks = data.NumberOfTriggerTracks();
		for (Int_t i = 0; i < ntracks; i++)
		{
			AliMUONTriggerTrack* track = data.TriggerTrack(i);
			tracksByIndex.Add(track->Clone());
		}
		
		// Now check that all the lists of trigger tracks contain the same results.
		// They must. If they do not then something is wrong with the implementation
		// of AliMUONDataInterface.
		if (tracksFromStore.GetEntriesFast() != tracksByIndex.GetEntriesFast())
		{
			Error(	"TriggerTracksOk",
				"The AliMUONDataInterface does not return all the trigger tracks"
				  " correctly through all its user interface methods. We got the"
				  " following numbers of tracks: %d and %d",
				tracksFromStore.GetEntriesFast(),
				tracksByIndex.GetEntriesFast()
			);
			return false;
		}
		for (Int_t i = 0; i < tracksFromStore.GetEntriesFast(); i++)
		{
			if (Compare(tracksFromStore[i], tracksByIndex[i]) != 0)
			{
				Error(	"TriggerTracksOk",
					"The AliMUONDataInterface does not return identical trigger"
					  " tracks through all its user interface methods. The"
					  " incorrect track has index %d.",
					i
				);
				tracksFromStore[i]->Print();
				tracksByIndex[i]->Print();
				return false;
			}
		}
	}
	return true;
}

/**
 * This method performs a check of the AliMUONDataInterface and AliMUONMCDataInterface
 * classes. Basically there are at least 2 ways to fetch data using these interfaces:
 * The expert way using the store objects returned by these interface classes or
 * the much slower but easier way of using the NumberOfxxx and Digit(...),
 * RawCluster(...), Track(...) etc. methods to fetch individual data objects.
 * The MUONCheckDI will check that all these various ways of fetching data results
 * in the same information being returned. If yes then kTRUE is returned and a
 * confirmation message is printed, if not then kFALSE is returned with the failure
 * reason printed to screen.
 */
bool MUONCheckDI()
{
	// TODO: complete checking AliMUONMCDataInterface.
	//cout << "Checking simulated triggers..." << endl;
	//if (! SimTriggersOk()) return false;
	//cout << "Simulated triggers look OK." << endl;
	
	cout << "Checking reconstructed digits..." << endl;
	if (! RecDigitsOk()) return false;
	cout << "Reconstructed digits look OK." << endl;
	
	cout << "Checking raw clusters..." << endl;
	if (! RawClustersOk()) return false;
	cout << "Raw clusters look OK." << endl;

	cout << "Checking reconstructed tracks..." << endl;
	if (! TracksOk()) return false;
	cout << "Reconstructed tracks look OK." << endl;

	cout << "Checking reconstructed triggers..." << endl;
	if (! TriggersOk()) return false;
	cout << "Reconstructed triggers look OK." << endl;

	cout << "Checking reconstructed trigger tracks..." << endl;
	if (! TriggerTracksOk()) return false;
	cout << "Reconstructed trigger tracks look OK." << endl;
	
	return true;
}
