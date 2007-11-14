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

/**
 * @file   AliHLTMUONRootifierComponent.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of the AliHLTMUONRootifierComponent component.
 */

#include "AliHLTMUONRootifierComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONRecHit.h"
#include "AliHLTMUONTriggerRecord.h"
#include "AliHLTMUONMansoTrack.h"
#include <cassert>

ClassImp(AliHLTMUONEvent);
ClassImp(AliHLTMUONRootifierComponent);


AliHLTMUONRootifierComponent::AliHLTMUONRootifierComponent() :
	AliHLTProcessor()
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


int AliHLTMUONRootifierComponent::DoInit(int /*argc*/, const char** /*argv*/)
{
	///
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	///
	
	return 0;
}


int AliHLTMUONRootifierComponent::DoDeinit()
{
	///
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	///
	
	return 0;
}


const char* AliHLTMUONRootifierComponent::GetComponentID()
{
	///
	/// Inherited from AliHLTComponent. Returns the component ID.
	///
	
	return "MUONRootifier";
}


AliHLTComponentDataType AliHLTMUONRootifierComponent::GetOutputDataType()
{
	///
	/// Inherited from AliHLTComponent. Returns the output data type.
	///
	
	return kAliHLTAnyDataType;
}


void AliHLTMUONRootifierComponent::GetInputDataTypes(
		vector<AliHLTComponentDataType>& list
	)
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
		AliHLTComponentTriggerData& /*trigData*/
	)
{
	///
	/// Inherited from AliHLTProcessor. Processes the new event data.
	///
	
	AliHLTMUONEvent event(evtData.fEventID);

	// First process the blocks of reconstructed hits and trigger records.
	for (int i = 0; i < GetNumberOfInputBlocks(); i++)
	{
		const AliHLTComponentBlockData* block = GetInputBlock(i);
		assert( block != NULL );
		if (block->fDataType == AliHLTMUONConstants::RecHitsBlockDataType())
		{
			AliHLTUInt8_t* ptr = reinterpret_cast<AliHLTUInt8_t*>(block->fPtr);
			ptr += block->fOffset;
			AliHLTMUONRecHitsBlockReader inblock(ptr, block->fSize);
			if (not inblock.BufferSizeOk())
			{
				size_t headerSize = sizeof(AliHLTMUONRecHitsBlockReader::HeaderType);
				if (block->fSize < headerSize)
				{
					HLTError("Received a reconstructed hits data block with a size of %d bytes,"
						" which is smaller than the minimum valid header size of %d bytes."
						" The block must be corrupt.",
						block->fSize, headerSize
					);
					continue;
				}
				
				size_t expectedWidth = sizeof(AliHLTMUONRecHitsBlockReader::ElementType);
				if (inblock.CommonBlockHeader().fRecordWidth != expectedWidth)
				{
					HLTError("Received a reconstructed hits data block with a record"
						" width of %d bytes, but the expected value is %d bytes."
						" The block might be corrupt.",
						block->fSize, headerSize
					);
					continue;
				}
				
				HLTError("Received a reconstructed hits data block with a size of %d bytes,"
					" but the block header claims the block should be %d bytes."
					" The block might be corrupt.",
					block->fSize, inblock.BytesUsed()
				);
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
				//AliHLTMUONRecHit rh(h.fX, h.fY, h.fZ, sourceDDL);
				//PushBack(&rh, "ROOTHITS", "MUON");
				event.Add(new AliHLTMUONRecHit(h.fX, h.fY, h.fZ, sourceDDL));
			}
		}
		else if (block->fDataType == AliHLTMUONConstants::TriggerRecordsBlockDataType())
		{
			AliHLTUInt8_t* ptr = reinterpret_cast<AliHLTUInt8_t*>(block->fPtr);
			ptr += block->fOffset;
			AliHLTMUONTriggerRecordsBlockReader inblock(ptr, block->fSize);
			if (not inblock.BufferSizeOk())
			{
				size_t headerSize = sizeof(AliHLTMUONTriggerRecordsBlockReader::HeaderType);
				if (block->fSize < headerSize)
				{
					HLTError("Received a trigger records data block with a size of %d bytes,"
						" which is smaller than the minimum valid header size of %d bytes."
						" The block must be corrupt.",
						block->fSize, headerSize
					);
					continue;
				}
				
				size_t expectedWidth = sizeof(AliHLTMUONTriggerRecordsBlockReader::ElementType);
				if (inblock.CommonBlockHeader().fRecordWidth != expectedWidth)
				{
					HLTError("Received a trigger records data block with a record"
						" width of %d bytes, but the expected value is %d bytes."
						" The block might be corrupt.",
						block->fSize, headerSize
					);
					continue;
				}
				
				HLTError("Received a trigger records data block with a size of %d bytes,"
					" but the block header claims the block should be %d bytes."
					" The block might be corrupt.",
					block->fSize, inblock.BytesUsed()
				);
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
			if (sourceDDL < 21 or sourceDDL > 22)
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
			// TODO: ignore for now, but should log an optional message.
		}
	}
	
	// Now we can look for tracks to add. We needed the ROOT trigger records
	// and reco hits created before we can create track objects.
	for (int i = 0; i < GetNumberOfInputBlocks(); i++)
	{
		const AliHLTComponentBlockData* block = GetInputBlock(i);
		assert( block != NULL );
		if (block->fDataType == AliHLTMUONConstants::MansoTracksBlockDataType())
		{
			AliHLTUInt8_t* ptr = reinterpret_cast<AliHLTUInt8_t*>(block->fPtr);
			ptr += block->fOffset;
			AliHLTMUONMansoTracksBlockReader inblock(ptr, block->fSize);
			if (not inblock.BufferSizeOk())
			{
				size_t headerSize = sizeof(AliHLTMUONMansoTracksBlockReader::HeaderType);
				if (block->fSize < headerSize)
				{
					HLTError("Received a Manso tracks data block with a size of %d bytes,"
						" which is smaller than the minimum valid header size of %d bytes."
						" The block must be corrupt.",
						block->fSize, headerSize
					);
					continue;
				}
				
				size_t expectedWidth = sizeof(AliHLTMUONMansoTracksBlockReader::ElementType);
				if (inblock.CommonBlockHeader().fRecordWidth != expectedWidth)
				{
					HLTError("Received a Manso tracks data block with a record"
						" width of %d bytes, but the expected value is %d bytes."
						" The block might be corrupt.",
						block->fSize, headerSize
					);
					continue;
				}
				
				HLTError("Received a Manso tracks data block with a size of %d bytes,"
					" but the block header claims the block should be %d bytes."
					" The block might be corrupt.",
					block->fSize, inblock.BytesUsed()
				);
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
		else
		{
			// TODO: ignore for now, but should log an optional message.
		}
	}
	
	PushBack(&event, "ROOTEVNT", "MUON");
	
	return 0;
}


void AliHLTMUONEvent::Print(Option_t* option) const
{
	///
	/// Inherited from TObject. Prints the contents of the event objects in fArray.
	///
	
	cout << "################## EVENT: " << fEventId << " ##################" << endl;
	for (Int_t i = 0; i < fArray.GetEntriesFast(); i++)
		if (fArray[i] != NULL) fArray[i]->Print(option);
}

