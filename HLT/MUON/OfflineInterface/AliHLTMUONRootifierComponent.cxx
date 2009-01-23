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

/* $Id$ */

///
/// @file   AliHLTMUONRootifierComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Implementation of the AliHLTMUONRootifierComponent component.
///

#include "AliHLTMUONRootifierComponent.h"
#include "AliHLTMUONEvent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONRecHit.h"
#include "AliHLTMUONTriggerRecord.h"
#include "AliHLTMUONMansoTrack.h"
#include "AliHLTMUONDecision.h"
#include "TClonesArray.h"
#include <cassert>

ClassImp(AliHLTMUONRootifierComponent);


AliHLTMUONRootifierComponent::AliHLTMUONRootifierComponent() :
	AliHLTMUONProcessor(),
	fWarnForUnexpecedBlock(false)
{
	///
	/// Default constructor.
	///
}


AliHLTMUONRootifierComponent::~AliHLTMUONRootifierComponent()
{
	///
	/// Default destructor.
	///
}


bool AliHLTMUONRootifierComponent::IgnoreArgument(const char* arg) const
{
	/// Return true if the argument is one of -cdbpath -run or -delaysetup
	/// to prevent the parent class from parsing these arguments in DoInit.
	
	if (strcmp(arg, "-cdbpath") == 0 or strcmp(arg, "-run") == 0 or
	    strcmp(arg, "-delaysetup") == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}


int AliHLTMUONRootifierComponent::DoInit(int argc, const char** argv)
{
	///
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	///
	
	HLTInfo("Initialising dHLT rootifier component.");

	// Inherit the parents functionality.
	int result = AliHLTMUONProcessor::DoInit(argc, argv);
	if (result != 0) return result;
	
	fWarnForUnexpecedBlock = false;
	
	for (int i = 0; i < argc; i++)
	{
		if (ArgumentAlreadyHandled(i, argv[i])) continue;

		if (strcmp(argv[i], "-warn_on_unexpected_block") == 0)
		{
			fWarnForUnexpecedBlock = true;
			continue;
		}

		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}
	
	return 0;
}


int AliHLTMUONRootifierComponent::DoDeinit()
{
	///
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	///
	
	HLTInfo("Deinitialising dHLT rootifier component.");
	return 0;
}


const char* AliHLTMUONRootifierComponent::GetComponentID()
{
	///
	/// Inherited from AliHLTComponent. Returns the component ID.
	///
	
	return AliHLTMUONConstants::RootifierComponentId();
}


AliHLTComponentDataType AliHLTMUONRootifierComponent::GetOutputDataType()
{
	/// Inherited from AliHLTComponent. Returns kAliHLTMultipleDataType
	/// refer to GetOutputDataTypes for all returned data types.
	
	return kAliHLTMultipleDataType;
}


int AliHLTMUONRootifierComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
	/// Inherited from AliHLTComponent. Returns the output data types.
	
	tgtList.push_back(kAliHLTAnyDataType);
	return tgtList.size();
}


void AliHLTMUONRootifierComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
	///
	/// Inherited from AliHLTProcessor. Returns the list of expected input data types.
	///
	
	list.push_back(kAliHLTAnyDataType);
}


void AliHLTMUONRootifierComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	///
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	///
	
	constBase = 1024*1024;
	inputMultiplier = 100;
}


AliHLTComponent* AliHLTMUONRootifierComponent::Spawn()
{
	///
	/// Inherited from AliHLTComponent. Creates a new object instance.
	///
	
	return new AliHLTMUONRootifierComponent();
}


int AliHLTMUONRootifierComponent::DoEvent(
		const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData
	)
{
	///
	/// Inherited from AliHLTProcessor. Processes the new event data.
	///
	
	AliHLTMUONEvent event(evtData.fEventID);
	const AliHLTComponentBlockData* block = NULL;
	AliHLTUInt32_t specification = 0;  // Contains the output data block spec bits.

	// First process the blocks of reconstructed hits and trigger records.
	for (int i = 0; i < GetNumberOfInputBlocks(); i++)
	{
		block = GetInputBlock(i);
		assert( block != NULL );
		
		HLTDebug("Handling block: %u, with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
			i, DataType2Text(block->fDataType).c_str(), block->fPtr, block->fSize
		);
		
		if (block->fDataType == AliHLTMUONConstants::RecHitsBlockDataType())
		{
			specification |= block->fSpecification;
			AliHLTMUONRecHitsBlockReader inblock(block->fPtr, block->fSize);
			if (not BlockStructureOk(inblock))
			{
				if (DumpDataOnError()) DumpEvent(evtData, trigData);
				continue;
			}
			
			// Decode the source DDL from the specification bits.
			Int_t sourceDDL = -1;
			bool ddl[22];
			AliHLTMUONUtils::UnpackSpecBits(block->fSpecification, ddl);
			for (int k = 0; k < 22; k++)
			{
				if (ddl[k])
				{
					if (sourceDDL == -1)
					{
						sourceDDL = k+1;
					}
					else
					{
						HLTWarning("The input data block %d contains"
							" data from multiple DDL sources.", i
						);
					}
				}
			}
			if (sourceDDL > 20)
			{
				HLTWarning("The source DDL for input data block %d is %d."
					" The expected range for the DDL is [1..20].",
					i, sourceDDL
				);
			}
			
			for (AliHLTUInt32_t n = 0; n < inblock.Nentries(); n++)
			{
				const AliHLTMUONRecHitStruct& h = inblock[n];
				AliHLTUInt8_t chamber;
				AliHLTUInt16_t detElemId;
				AliHLTMUONUtils::UnpackRecHitFlags(h.fFlags, chamber, detElemId);
				event.Add(new AliHLTMUONRecHit(h.fX, h.fY, h.fZ, sourceDDL, detElemId, chamber));
			}
		}
		else if (block->fDataType == AliHLTMUONConstants::TriggerRecordsBlockDataType())
		{
			specification |= block->fSpecification;
			AliHLTMUONTriggerRecordsBlockReader inblock(block->fPtr, block->fSize);
			if (not BlockStructureOk(inblock))
			{
				if (DumpDataOnError()) DumpEvent(evtData, trigData);
				continue;
			}
			
			// Decode the source DDL from the specification bits.
			Int_t sourceDDL = -1;
			bool ddl[22];
			AliHLTMUONUtils::UnpackSpecBits(block->fSpecification, ddl);
			for (int k = 0; k < 22; k++)
			{
				if (ddl[k])
				{
					if (sourceDDL == -1)
					{
						sourceDDL = k+1;
					}
					else
					{
						HLTWarning("The input data block %d contains"
							" data from multiple DDL sources.", i
						);
					}
				}
			}
			if (sourceDDL != -1 and (sourceDDL < 21 or sourceDDL > 22))
			{
				HLTWarning("The source DDL for input data block %d is %d."
					" The expected range for the DDL is [21..22].",
					i, sourceDDL
				);
			}
			
			for (AliHLTUInt32_t n = 0; n < inblock.Nentries(); n++)
			{
				const AliHLTMUONTriggerRecordStruct& t = inblock[n];
				
				AliHLTMUONParticleSign sign;
				bool hitset[4];
				AliHLTMUONUtils::UnpackTriggerRecordFlags(
						t.fFlags, sign, hitset
					);
			
				AliHLTMUONTriggerRecord* tr = new AliHLTMUONTriggerRecord(
						t.fId, sign, t.fPx, t.fPy, t.fPz, sourceDDL
					);
				for (int k = 0; k < 4; k++)
					tr->SetHit(k+11, t.fHit[k].fX, t.fHit[k].fY, t.fHit[k].fZ);
				event.Add(tr);
			}
		}
		else
		{
			if (block->fDataType != AliHLTMUONConstants::MansoTracksBlockDataType() and
			    block->fDataType != AliHLTMUONConstants::SinglesDecisionBlockDataType() and
			    block->fDataType != AliHLTMUONConstants::PairsDecisionBlockDataType()
			   )
			{
				// Log a message indicating that we got a data block that we
				// do not know how to handle.
				if (fWarnForUnexpecedBlock)
					HLTWarning("Received a data block of a type we cannot handle: '%s', spec: 0x%X",
						DataType2Text(block->fDataType).c_str(), block->fSpecification
					);
#ifdef __DEBUG
				else
					HLTDebug("Received a data block of a type we cannot handle: '%s', spec: 0x%X",
						DataType2Text(block->fDataType).c_str(), block->fSpecification
					);
#endif
			}
		}
	}
	
	// Now we can look for tracks to add. We needed the ROOT trigger records
	// and reco hits created before we can create track objects.
	for (block = GetFirstInputBlock(AliHLTMUONConstants::MansoTracksBlockDataType());
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		specification |= block->fSpecification;
		AliHLTMUONMansoTracksBlockReader inblock(block->fPtr, block->fSize);
		if (not BlockStructureOk(inblock))
		{
			if (DumpDataOnError()) DumpEvent(evtData, trigData);
			continue;
		}
		
		for (AliHLTUInt32_t n = 0; n < inblock.Nentries(); n++)
		{
			const AliHLTMUONMansoTrackStruct& t = inblock[n];
			
			AliHLTMUONParticleSign sign;
			bool hitset[4];
			AliHLTMUONUtils::UnpackMansoTrackFlags(
					t.fFlags, sign, hitset
				);
			
			// Try find the trigger record in 'event'.
			const AliHLTMUONTriggerRecord* trigrec = NULL;
			for (Int_t k = 0; k < event.Array().GetEntriesFast(); k++)
			{
				if (event.Array()[k]->IsA() != AliHLTMUONTriggerRecord::Class())
					continue;
				const AliHLTMUONTriggerRecord* tk =
					static_cast<const AliHLTMUONTriggerRecord*>(event.Array()[k]);
				if (tk->Id() == t.fTrigRec)
				{
					trigrec = tk;
					break;
				}
			}
			
			// Now try find the hits in 'event'.
			// If they cannot be found then create new ones.
			const AliHLTMUONRecHit* hit7 = NULL;
			const AliHLTMUONRecHit* hit8 = NULL;
			const AliHLTMUONRecHit* hit9 = NULL;
			const AliHLTMUONRecHit* hit10 = NULL;
			for (Int_t k = 0; k < event.Array().GetEntriesFast(); k++)
			{
				if (event.Array()[k]->IsA() != AliHLTMUONRecHit::Class())
					continue;
				const AliHLTMUONRecHit* h =
					static_cast<const AliHLTMUONRecHit*>(event.Array()[k]);
				
				if (hitset[0] and h->X() == t.fHit[0].fX and h->Y() == t.fHit[0].fY
					and h->Z() == t.fHit[0].fZ)
				{
					hit7 = h;
				}
				if (hitset[1] and h->X() == t.fHit[1].fX and h->Y() == t.fHit[1].fY
					and h->Z() == t.fHit[1].fZ)
				{
					hit8 = h;
				}
				if (hitset[2] and h->X() == t.fHit[2].fX and h->Y() == t.fHit[2].fY
					and h->Z() == t.fHit[2].fZ)
				{
					hit9 = h;
				}
				if (hitset[3] and h->X() == t.fHit[3].fX and h->Y() == t.fHit[3].fY
					and h->Z() == t.fHit[3].fZ)
				{
					hit10 = h;
				}
			}
			AliHLTMUONRecHit* newhit;
			if (hitset[0] and hit7 == NULL)
			{
				newhit = new AliHLTMUONRecHit(t.fHit[0].fX, t.fHit[0].fY, t.fHit[0].fZ);
				event.Add(newhit);
				hit7 = newhit;
			}
			if (hitset[1] and hit8 == NULL)
			{
				newhit = new AliHLTMUONRecHit(t.fHit[1].fX, t.fHit[1].fY, t.fHit[1].fZ);
				event.Add(newhit);
				hit8 = newhit;
			}
			if (hitset[2] and hit9 == NULL)
			{
				newhit = new AliHLTMUONRecHit(t.fHit[2].fX, t.fHit[2].fY, t.fHit[2].fZ);
				event.Add(newhit);
				hit9 = newhit;
			}
			if (hitset[3] and hit10 == NULL)
			{
				newhit = new AliHLTMUONRecHit(t.fHit[3].fX, t.fHit[3].fY, t.fHit[3].fZ);
				event.Add(newhit);
				hit10 = newhit;
			}
		
			AliHLTMUONMansoTrack* tr = new AliHLTMUONMansoTrack(
					t.fId, sign, t.fPx, t.fPy, t.fPz, t.fChi2,
					trigrec, hit7, hit8, hit9, hit10
				);
			event.Add(tr);
		}
	}
	
	bool decisionBlockFound = false;
	UInt_t numLowPt = 0;
	UInt_t numHighPt = 0;
	TClonesArray singlesDecisions("AliHLTMUONDecision::AliTrackDecision");
	
	// Find the single tracks decision blocks and add their information.
	// We just sum the trigger scalars and single decisions.
	for (block = GetFirstInputBlock(AliHLTMUONConstants::SinglesDecisionBlockDataType());
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		decisionBlockFound = true;
		specification |= block->fSpecification;
		AliHLTMUONSinglesDecisionBlockReader inblock(block->fPtr, block->fSize);
		if (not BlockStructureOk(inblock))
		{
			if (DumpDataOnError()) DumpEvent(evtData, trigData);
			continue;
		}
		
		numLowPt += inblock.BlockHeader().fNlowPt;
		numHighPt += inblock.BlockHeader().fNhighPt;
		
		for (AliHLTUInt32_t n = 0; n < inblock.Nentries(); n++)
		{
			const AliHLTMUONTrackDecisionStruct& t = inblock[n];
			
			bool highPt, lowPt;
			AliHLTMUONUtils::UnpackTrackDecisionBits(t.fTriggerBits, highPt, lowPt);
			
			// Try find the corresponding track in the 'event'.
			const AliHLTMUONMansoTrack* track = NULL;
			for (Int_t k = 0; k < event.Array().GetEntriesFast(); k++)
			{
				if (event.Array()[k]->IsA() != AliHLTMUONMansoTrack::Class())
					continue;
				const AliHLTMUONMansoTrack* tk =
					static_cast<const AliHLTMUONMansoTrack*>(event.Array()[k]);
				if (tk->Id() == t.fTrackId)
				{
					track = tk;
					break;
				}
			}
			
			// If the track was not found then create a dummy one.
			if (track == NULL)
			{
				AliHLTMUONMansoTrack* tr = new AliHLTMUONMansoTrack(t.fTrackId);
				event.Add(tr);
				track = tr;
			}
			
			new (singlesDecisions[singlesDecisions.GetEntriesFast()])
				AliHLTMUONDecision::AliTrackDecision(t.fPt, lowPt, highPt, track);
		}
	}
	
	UInt_t numUnlikeAnyPt = 0;
	UInt_t numUnlikeLowPt = 0;
	UInt_t numUnlikeHighPt = 0;
	UInt_t numLikeAnyPt = 0;
	UInt_t numLikeLowPt = 0;
	UInt_t numLikeHighPt = 0;
	UInt_t numAnyMass = 0;
	UInt_t numLowMass = 0;
	UInt_t numHighMass = 0;
	TClonesArray pairsDecisions("AliHLTMUONDecision::AliPairDecision");
	
	// Find the track pairs decision blocks and add their information.
	// We just sum the trigger scalars and track pair decisions.
	for (block = GetFirstInputBlock(AliHLTMUONConstants::PairsDecisionBlockDataType());
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		decisionBlockFound = true;
		specification |= block->fSpecification;
		AliHLTMUONPairsDecisionBlockReader inblock(block->fPtr, block->fSize);
		if (not BlockStructureOk(inblock))
		{
			if (DumpDataOnError()) DumpEvent(evtData, trigData);
			continue;
		}
		
		numUnlikeAnyPt += inblock.BlockHeader().fNunlikeAnyPt;
		numUnlikeLowPt += inblock.BlockHeader().fNunlikeLowPt;
		numUnlikeHighPt += inblock.BlockHeader().fNunlikeHighPt;
		numLikeAnyPt += inblock.BlockHeader().fNlikeAnyPt;
		numLikeLowPt += inblock.BlockHeader().fNlikeLowPt;
		numLikeHighPt += inblock.BlockHeader().fNlikeHighPt;
		numAnyMass += inblock.BlockHeader().fNmassAny;
		numLowMass += inblock.BlockHeader().fNmassLow;
		numHighMass += inblock.BlockHeader().fNmassHigh;
		
		for (AliHLTUInt32_t n = 0; n < inblock.Nentries(); n++)
		{
			const AliHLTMUONPairDecisionStruct& t = inblock[n];
			
			bool highMass, lowMass, unlike;
			AliHLTUInt8_t highPtCount, lowPtCount;
			AliHLTMUONUtils::UnpackPairDecisionBits(
					t.fTriggerBits, highMass, lowMass, unlike,
					highPtCount, lowPtCount
				);
			
			// Try find the corresponding tracks in the 'event'.
			const AliHLTMUONMansoTrack* trackA = NULL;
			const AliHLTMUONMansoTrack* trackB = NULL;
			for (Int_t k = 0; k < event.Array().GetEntriesFast(); k++)
			{
				if (event.Array()[k]->IsA() != AliHLTMUONMansoTrack::Class())
					continue;
				const AliHLTMUONMansoTrack* tk =
					static_cast<const AliHLTMUONMansoTrack*>(event.Array()[k]);
				if (tk->Id() == t.fTrackAId) trackA = tk;
				if (tk->Id() == t.fTrackBId) trackB = tk;
				if (trackA != NULL and trackB != NULL) break;
			}
			
			// If either of the tracks was not found then create a dummy one.
			if (trackA == NULL)
			{
				AliHLTMUONMansoTrack* tr = new AliHLTMUONMansoTrack(t.fTrackAId);
				event.Add(tr);
				trackA = tr;
			}
			if (trackB == NULL)
			{
				AliHLTMUONMansoTrack* tr = new AliHLTMUONMansoTrack(t.fTrackBId);
				event.Add(tr);
				trackB = tr;
			}
			
			new (pairsDecisions[pairsDecisions.GetEntriesFast()])
				AliHLTMUONDecision::AliPairDecision(
					t.fInvMass, lowMass, highMass, unlike,
					lowPtCount, highPtCount, trackA, trackB
				);
		}
	}
	
	
	// Do not add the decision if no decision blocks were found.
	if (decisionBlockFound)
	{
		AliHLTMUONDecision* triggerDecision = new AliHLTMUONDecision(
				numLowPt, numHighPt, numUnlikeAnyPt, numUnlikeLowPt,
				numUnlikeHighPt, numLikeAnyPt, numLikeLowPt,
				numLikeHighPt, numAnyMass, numLowMass, numHighMass
			);
		for (Int_t i = 0; i < singlesDecisions.GetEntriesFast(); i++)
		{
			AliHLTMUONDecision::AliTrackDecision* decision =
				static_cast<AliHLTMUONDecision::AliTrackDecision*>( singlesDecisions[i] );
			triggerDecision->AddDecision(decision);
		}
		for (Int_t j = 0; j < pairsDecisions.GetEntriesFast(); j++)
		{
			AliHLTMUONDecision::AliPairDecision* decision =
				static_cast<AliHLTMUONDecision::AliPairDecision*>( pairsDecisions[j] );
			triggerDecision->AddDecision(decision);
		}
		
		event.Add(triggerDecision);
	}
	
	PushBack(&event, AliHLTMUONConstants::RootifiedEventDataType(), specification);
	
	return 0;
}

