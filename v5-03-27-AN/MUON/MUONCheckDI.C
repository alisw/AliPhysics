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

/// \ingroup macros
/// \file MUONCheckDI.C
/// \brief Macro for checking AliMUONDataInterface and AliMUONMCDataInterface.
///
/// \author Bruce Becker, DAPNIA/SPhN/CEA Saclay
///
/// Modified to updated versions of data interfaces
/// by Artur Szostak <artursz@iafrica.com> (University of Cape Town)

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliMUONHit.h"
#include "AliMUONVDigit.h"
#include "AliMUONVCluster.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONDataInterface.h"
#include "AliMUONVHitStore.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMpConstants.h"
#include "AliMpDEManager.h"

#include "AliCDBManager.h"

#include <Rtypes.h>
#include <Riostream.h>
#include <TObjArray.h>
#include <TIterator.h>
#include <TMatrixD.h>

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
	if (a->IsA() == AliMUONLocalTrigger::Class() && b->IsA() == AliMUONLocalTrigger::Class())
	{
		result = memcmp(a, b, sizeof(AliMUONLocalTrigger));
	}
	else if (a->IsA() == AliMUONRegionalTrigger::Class() && b->IsA() == AliMUONRegionalTrigger::Class())
	{
		const AliMUONRegionalTrigger* ra = static_cast<const AliMUONRegionalTrigger*>(a);
		const AliMUONRegionalTrigger* rb = static_cast<const AliMUONRegionalTrigger*>(b);
		if (ra->GetId() < rb->GetId()) return -1;
		if (ra->GetId() > rb->GetId()) return 1;
		if (ra->GetLocalOutput(0) < rb->GetLocalOutput(0)) return -1;
		if (ra->GetLocalOutput(0) > rb->GetLocalOutput(0)) return 1;
		if (ra->GetLocalOutput(1) < rb->GetLocalOutput(1)) return -1;
		if (ra->GetLocalOutput(1) > rb->GetLocalOutput(1)) return 1;
		if (ra->GetLocalMask() < rb->GetLocalMask()) return -1;
		if (ra->GetLocalMask() > rb->GetLocalMask()) return 1;
		if (ra->GetOutput() < rb->GetOutput()) return -1;
		if (ra->GetOutput() > rb->GetOutput()) return 1;
		return 0;
	}
	else if (a->IsA() == AliMUONGlobalTrigger::Class() && b->IsA() == AliMUONGlobalTrigger::Class())
	{
		result = memcmp(a, b, sizeof(AliMUONGlobalTrigger));
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
 * error message and false is returned.
 */
bool SimTriggersOk()
{
	AliMUONMCDataInterface data("generated/galice.root");
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray localsFromStore, regionalsFromStore;
		localsFromStore.SetOwner(kTRUE);
		regionalsFromStore.SetOwner(kTRUE);
		AliMUONVTriggerStore* store = data.TriggerStore(event);
		if (store == NULL) return false;
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
				regionalsFromStore[i]->Print();
				regionalsByIndex[i]->Print();
		}
	}
	return true;
}

/**
 * This method fills internal arrays with s-digits returned by the AliMUONMCDataInterface.
 * For each set of interface methods available a TObjArray is filled with copies of 
 * the s-digits. These arrays are sorted and then compared to each other. The arrays
 * should contain the same s-digit information if everything is working correctly with
 * AliMUONMCDataInterface. If not then the difference is printed together with an
 * error message and false is returned.
 */
bool SimSDigitsOk()
{
	AliMUONMCDataInterface data("generated/galice.root");
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray digitsFromStore;
		digitsFromStore.SetOwner(kTRUE);
		AliMUONVDigitStore* store = data.SDigitStore(event);
		if (store == NULL) return false;
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
			Int_t ndigits = data.NumberOfSDigits(detElem);
			for (Int_t i = 0; i < ndigits; i++)
			{
				AliMUONVDigit* digit = data.SDigit(detElem, i);
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
			Int_t ndigits = data.NumberOfSDigits(chamber, cathode);
			for (Int_t i = 0; i < ndigits; i++)
			{
				AliMUONVDigit* digit = data.SDigit(chamber, cathode, i);
				digitsByChamber.Add(digit->Clone());
			}
		}
		digitsByChamber.Sort();
		
		// Now check that all the lists of s-digits contain the same results.
		// They must. If they do not then something is wrong with the implementation
		// of AliMUONMCDataInterface.
		if (digitsFromStore.GetEntriesFast() != digitsByDetElem.GetEntriesFast()
		    || digitsFromStore.GetEntriesFast() != digitsByChamber.GetEntriesFast())
		{
			Error(	"SimSDigitsOk",
				"The AliMUONMCDataInterface does not return all the s-digits correctly"
				  " through all its user interface methods. We got the following"
				  " numbers of s-digits: %d, %d and %d",
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
				Error(	"SimSDigitsOk",
					"The AliMUONMCDataInterface does not return identical s-digits"
					  " through all its user interface methods. The incorrect"
					  " s-digit has index %d after sorting.",
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
 * This method fills internal arrays with digits returned by the AliMUONMCDataInterface.
 * For each set of interface methods available a TObjArray is filled with copies of 
 * the digits. These arrays are sorted and then compared to each other. The arrays
 * should contain the same digit information if everything is working correctly with
 * AliMUONMCDataInterface. If not then the difference is printed together with an
 * error message and false is returned.
 */
bool SimDigitsOk()
{
	AliMUONMCDataInterface data("generated/galice.root");
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray digitsFromStore;
		digitsFromStore.SetOwner(kTRUE);
		AliMUONVDigitStore* store = data.DigitStore(event);
		if (store == NULL) return false;
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
		// of AliMUONMCDataInterface.
		if (digitsFromStore.GetEntriesFast() != digitsByDetElem.GetEntriesFast()
		    || digitsFromStore.GetEntriesFast() != digitsByChamber.GetEntriesFast())
		{
			Error(	"SimDigitsOk",
				"The AliMUONMCDataInterface does not return all the digits correctly"
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
				Error(	"SimDigitsOk",
					"The AliMUONMCDataInterface does not return identical digits"
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
 * This method fills internal arrays with hits returned by the AliMUONMCDataInterface.
 * For each set of interface methods available a TObjArray is filled with copies of 
 * the hits. These arrays are then compared to each other. The arrays should contain
 * the same hit information if everything is working correctly with AliMUONMCDataInterface.
 * If not then the difference is printed together with an error message and false is returned.
 */
bool SimHitsOk()
{
	AliMUONMCDataInterface data("generated/galice.root");
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		if (data.NumberOfTracks(event) != data.NumberOfTracks())
		{
			Error(	"SimHitsOk",
				"The AliMUONMCDataInterface does not return the same number of tracks"
				  " through all its user interface methods. We got the following"
				  " numbers of tracks: %d and %d",
				data.NumberOfTracks(event),
				data.NumberOfTracks()
			);
			return false;
		}
		
		for (Int_t track = 0; track < data.NumberOfTracks(); track++)
		{
			TObjArray hitsFromStore;
			hitsFromStore.SetOwner(kTRUE);
			AliMUONVHitStore* store = data.HitStore(event, track);
			if (store == NULL) return false;
			TIter next(store->CreateIterator());
			AliMUONHit* hit;
			while ( (hit = static_cast<AliMUONHit*>( next() )) != NULL )
			{
				hitsFromStore.Add(hit->Clone());
			}
			//hitsFromStore.Sort();  // Unfortunately hits do not implement the Compare method.
			
			TObjArray hitsByMethod;
			hitsByMethod.SetOwner(kTRUE);
			data.GetEvent(event);
			for (Int_t i = 0; i < data.NumberOfHits(track); i++)
			{
				AliMUONHit* hit = data.Hit(track, i);
				hitsByMethod.Add(hit->Clone());
			}
			//hitsByMethod.Sort();  // Unfortunately hits do not implement the Compare method.
			
			// Now check that both lists of hits contain the same results.
			// They must. If they do not then something is wrong with the implementation
			// of AliMUONMCDataInterface.
			if (hitsFromStore.GetEntriesFast() != hitsByMethod.GetEntriesFast())
			{
				Error(	"SimHitsOk",
					"The AliMUONMCDataInterface does not return all the hits correctly"
					" through all its user interface methods. We got the following"
					" numbers of hits: %d and %d",
					hitsFromStore.GetEntriesFast(),
					hitsByMethod.GetEntriesFast()
				);
				return false;
			}
			for (Int_t i = 0; i < hitsFromStore.GetEntriesFast(); i++)
			{
				if (Compare(hitsFromStore[i], hitsByMethod[i]) != 0)
				{
					Error(	"SimHitsOk",
						"The AliMUONMCDataInterface does not return identical hits"
						" through all its user interface methods. The incorrect"
						" hit has index %d after sorting, for track %d.",
						i, track
					);
					hitsFromStore[i]->Print();
					hitsByMethod[i]->Print();
					return false;
				}
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
 * error message and false is returned.
 */
bool RecDigitsOk()
{
	AliMUONDataInterface data;
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray digitsFromStore;
		digitsFromStore.SetOwner(kTRUE);
		AliMUONVDigitStore* store = data.DigitStore(event);
		if (store == NULL) return false;
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
 * together with an error message and false is returned.
 */
bool RawClustersOk()
{
	AliMUONDataInterface data;
	for (Int_t event = 0; event < data.NumberOfEvents(); event++)
	{
		TObjArray clustersFromStore;
		clustersFromStore.SetOwner(kTRUE);
		AliMUONVClusterStore* store = data.ClusterStore(event);
		if (store == NULL) return false;
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
		for (Int_t chamber = 0; chamber < AliMpConstants::NofTrackingChambers(); chamber++)
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
 * This method fills internal arrays with local and regional triggers returned
 * by AliMUONDataInterface. For each set of interface methods available in
 * AliMUONDataInterface a TObjArray is created for local and another for regional
 * triggers. These arrays are filled with copies of the trigger objects.
 * The global trigger object is also copied out using the 2 different methods.
 * The arrays and objects are then compared to each other. The arrays and objects
 * should contain the same information if everything is working correctly with
 * AliMUONDataInterface. If not then the difference is printed together with an
 * error message and false is returned.
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
		if (store == NULL) return false;
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
bool MUONCheckDI(bool checkSim = true, bool checkRec = true)
{
	AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

	// Note: we do not bother checking the AliMUONMCDataInterface::Particle,
	// AliMUONMCDataInterface::Stack and AliMUONMCDataInterface::TrackRefs methods
	// because they are trivial enough to validate from a quick inspecition of
	// the source code.
	
	if (checkSim)
	{
		cout << "Checking simulated hits..." << endl;
		if (! SimHitsOk()) return false;
		cout << "Simulated hits look OK." << endl;
		
		cout << "Checking simulated s-digits..." << endl;
		if (! SimSDigitsOk()) return false;
		cout << "Simulated s-digits look OK." << endl;
		
		cout << "Checking simulated digits..." << endl;
		if (! SimDigitsOk()) return false;
		cout << "Simulated digits look OK." << endl;
		
		cout << "Checking simulated triggers..." << endl;
		if (! SimTriggersOk()) return false;
		cout << "Simulated triggers look OK." << endl;
	}
	
	if (checkRec)
	{
		cout << "Checking reconstructed digits..." << endl;
		if (! RecDigitsOk()) return false;
		cout << "Reconstructed digits look OK." << endl;
		
		cout << "Checking raw clusters..." << endl;
		if (! RawClustersOk()) return false;
		cout << "Raw clusters look OK." << endl;
	
		cout << "Checking reconstructed triggers..." << endl;
		if (! TriggersOk()) return false;
		cout << "Reconstructed triggers look OK." << endl;
	}

	return true;
}
