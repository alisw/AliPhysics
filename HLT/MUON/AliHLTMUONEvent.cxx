/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id: $

///
/// @file   AliHLTMUONEvent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Implementation of the AliHLTMUONEvent class.
///
/// The class is used to store all ROOTified data objects from the dHLT chain
/// for a single event together.

#include "AliHLTMUONEvent.h"
#include "AliHLTMUONRecHit.h"
#include "AliHLTMUONTriggerRecord.h"
#include "AliHLTMUONMansoTrack.h"
#include "AliHLTMUONTrack.h"
#include "AliHLTMUONDecision.h"
#include "AliLog.h"
#include <iostream>
#include <map>

ClassImp(AliHLTMUONEvent);


AliHLTMUONEvent::AliHLTMUONEvent(const AliHLTMUONEvent& event) :
	TObject(event),
	fEventId(event.fEventId),
	fArray()
{
	/// Copy constructor performs a deep copy of the object.
	
	fArray.SetOwner(kTRUE);
	DeepCopy(event);
}


AliHLTMUONEvent& AliHLTMUONEvent::operator = (const AliHLTMUONEvent& event)
{
	/// The assignment operator performs a deep copy of the object.
	
	if (this == &event) return *this;
	fArray.Clear();
	TObject::operator = (event);
	DeepCopy(event);
	return *this;
}


const AliHLTMUONDecision* AliHLTMUONEvent::FindDecision() const
{
	/// Finds the decision object in the event from the list of dHLT objects.
	/// There should only be one such object in the event. If not, then only
	/// the first object found is returned. You will need to manually search
	/// for the other objects.
	/// \returns  The AliHLTMUONDecision object in the event or NULL if none exists.
	
	for (Int_t i = 0; i < fArray.GetEntriesFast(); i++)
	{
		if (fArray[i]->IsA() == AliHLTMUONDecision::Class())
		{
			return static_cast<const AliHLTMUONDecision*>(fArray[i]);
		}
	}
	
	return NULL;
}


void AliHLTMUONEvent::Print(Option_t* option) const
{
	///
	/// Inherited from TObject. Prints the contents of the event objects in fArray.
	/// \param option  This is an option string that is just passed on to individual
	///     objects in the event's fArray list of objects.
	///
	
	std::cout << "################## EVENT: " << fEventId << " ##################" << std::endl;
	for (Int_t i = 0; i < fArray.GetEntriesFast(); i++)
		if (fArray[i] != NULL) fArray[i]->Print(option);
}


void AliHLTMUONEvent::Clear(Option_t* option)
{
	/// Clears the internal array of event objects.
	
	fEventId = AliHLTEventID_t(-1);
	fArray.Clear(option);
}


void AliHLTMUONEvent::Copy(TObject& object) const
{
	/// Deep copy this object to the target object.
	/// \param object  The target object to copy to.
	
	if (object.IsA() != this->IsA())
	{
		AliError(Form("Cannot copy an object of type %s to a type of %s.",
			this->ClassName(), object.ClassName()
		));
		return;
	}
	
	TObject::Copy(object);
	AliHLTMUONEvent& event = static_cast<AliHLTMUONEvent&>(object);
	event.DeepCopy(*this);
}


void AliHLTMUONEvent::DeepCopy(const AliHLTMUONEvent& event)
{
	/// Performs a deep copy of the event.
	
	fEventId = event.fEventId;
	TObjArray tocopy = event.fArray;
	std::map<const TObject*, TObject*> objmap;
	
	// We need to copy all the objects from the old event while maintaining the
	// pointer cross references contained in the track and decision objects:
	// Start by looping over all objects and first copy the trigger records and hits.
	TIter next(&tocopy);
	TObject* obj = NULL;
	while ( (obj = next()) != NULL )
	{
		if (obj->IsA() == AliHLTMUONTriggerRecord::Class() or
		    obj->IsA() == AliHLTMUONRecHit::Class()
		   )
		{
			TObject* newobj = obj->Clone();
			objmap[obj] = newobj;
			fArray.Add(newobj);
			tocopy.Remove(obj);
		}
	}
	
	// Now loop over all objects that still need to be copied and copy the tracks.
	next.Reset();
	while ( (obj = next()) != NULL )
	{
		if (obj->IsA() == AliHLTMUONMansoTrack::Class())
		{
			AliHLTMUONMansoTrack* track = static_cast<AliHLTMUONMansoTrack*>(obj);
			AliHLTMUONMansoTrack* newtrack = new AliHLTMUONMansoTrack(
				track->Id(), track->Sign(),
				track->Px(), track->Py(), track->Pz(),
				track->Chi2(),
				static_cast<const AliHLTMUONTriggerRecord*>( objmap[track->TriggerRecord()] ),
				static_cast<const AliHLTMUONRecHit*>( objmap[track->Hit(7)] ),
				static_cast<const AliHLTMUONRecHit*>( objmap[track->Hit(8)] ),
				static_cast<const AliHLTMUONRecHit*>( objmap[track->Hit(9)] ),
				static_cast<const AliHLTMUONRecHit*>( objmap[track->Hit(10)] ),
				track->Zmiddle(), track->QBL()
			);
			for (int i = 7; i <= 10; ++i)
			{
				const TVector3& b = track->RoICentre(i);
				newtrack->SetRoI(i, b.X(), b.Y(), b.Z(), track->RoIRadius(i));
			}
			objmap[obj] = newtrack;
			fArray.Add(newtrack);
			tocopy.Remove(obj);
		}
		else if (obj->IsA() == AliHLTMUONTrack::Class())
		{
			AliHLTMUONTrack* track = static_cast<AliHLTMUONTrack*>(obj);
			
			const AliHLTMUONRecHit* newhits[16];
			for (int i = 0; i < 16; ++i)
			{
				newhits[i] = static_cast<const AliHLTMUONRecHit*>( objmap[track->Hit(i)] );
			}
			AliHLTMUONTrack* newtrack = new AliHLTMUONTrack(
				track->Id(), track->Sign(),
				track->Px(), track->Py(), track->Pz(),
				track->InverseBendingMomentum(),
				track->ThetaX(), track->ThetaY(),
				track->X(), track->Y(), track->Z(),
				track->Chi2(),
				static_cast<const AliHLTMUONTriggerRecord*>( objmap[track->TriggerRecord()] ),
				newhits
			);
			objmap[obj] = newtrack;
			fArray.Add(newtrack);
			tocopy.Remove(obj);
		}
	}
	
	// Finally copy over the decision objects.
	next.Reset();
	while ( (obj = next()) != NULL )
	{
		if (obj->IsA() == AliHLTMUONDecision::Class())
		{
			AliHLTMUONDecision* dec = static_cast<AliHLTMUONDecision*>(obj);
			AliHLTMUONDecision* newdec = new AliHLTMUONDecision(
				dec->NumberOfLowPtTriggers(),
				dec->NumberOfHighPtTriggers(),
				dec->NumberOfUnlikePairs(),
				dec->NumberOfUnlikeLowPtPairs(),
				dec->NumberOfUnlikeHighPtPairs(),
				dec->NumberOfLikePairs(),
				dec->NumberOfLikeLowPtPairs(),
				dec->NumberOfLikeHighPtPairs(),
				dec->NumberOfMassTriggers(),
				dec->NumberOfLowMassTriggers(),
				dec->NumberOfHighMassTriggers()
			);
			for (Int_t i = 0; i < dec->NumberOfTracks(); ++i)
			{
				const AliHLTMUONDecision::AliTrackDecision* d = dec->SingleTrackDecision(i);
				newdec->AddDecision(
					d->Pt(), d->PassedLowPtCut(), d->PassedHighPtCut(),
					objmap[d->Track()]
				);
			}
			for (Int_t i = 0; i < dec->NumberOfPairs(); ++i)
			{
				const AliHLTMUONDecision::AliPairDecision* d = dec->TrackPairDecision(i);
				newdec->AddDecision(
					d->Mass(), d->PassedLowMassCut(),
					d->PassedHighMassCut(), d->UnlikeSign(),
					d->NumberPassedLowPtCut(), d->NumberPassedHighPtCut(),
					objmap[d->TrackA()], objmap[d->TrackB()]
				);
			}
			objmap[obj] = newdec;
			fArray.Add(newdec);
			tocopy.Remove(obj);
		}
	}
	
	// Copy all the remaining objects that we do not handle in a special way.
	next.Reset();
	while ( (obj = next()) != NULL )
	{
		fArray.Add(obj->Clone());
	}
}
