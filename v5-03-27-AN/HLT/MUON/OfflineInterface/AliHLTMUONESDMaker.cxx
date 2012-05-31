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
/// @file   AliHLTMUONESDMaker.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   30 June 2008
/// @brief  Implementation of the AliHLTMUONESDMaker component.
///
/// The ESD maker component converts dHLT raw internal reconstructed information
/// into AliESDEvent objects.
///

#include "AliHLTMUONESDMaker.h"
#include "AliHLTMUONEvent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONRecHit.h"
#include "AliHLTMUONTriggerRecord.h"
#include "AliHLTMUONMansoTrack.h"
#include "AliHLTMUONDecision.h"
#include "AliMUONConstants.h"
#include "AliMUONVCluster.h"
#include "AliESDEvent.h"
#include "AliESDRun.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "TClonesArray.h"
#include "TList.h"
#include <cmath>
#include <cassert>


ClassImp(AliHLTMUONESDMaker);


AliHLTMUONESDMaker::AliHLTMUONESDMaker() :
	AliHLTMUONProcessor(),
	fWarnForUnexpecedBlock(false),
	fMakeMinimalESD(false),
	fAddCustomData(false),
	fMakeClonesArray(false),
	fMakeESDDataBlock(true),
	fClusterIndex(0)
{
	/// Default constructor.
}


AliHLTMUONESDMaker::~AliHLTMUONESDMaker()
{
	/// Default destructor.
}


bool AliHLTMUONESDMaker::IgnoreArgument(const char* arg) const
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


int AliHLTMUONESDMaker::DoInit(int argc, const char** argv)
{
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	
	HLTInfo("Initialising dHLT ESD maker component.");

	// Inherit the parents functionality.
	int result = AliHLTMUONProcessor::DoInit(argc, argv);
	if (result != 0) return result;
	
	fWarnForUnexpecedBlock = false;
	fMakeMinimalESD = false;
	fMakeClonesArray = false;
	fMakeESDDataBlock = true;
	
	for (int i = 0; i < argc; i++)
	{
		if (ArgumentAlreadyHandled(i, argv[i])) continue;

		if (strcmp(argv[i], "-make_minimal_esd") == 0)
		{
			fMakeMinimalESD = true;
			continue;
		}
		
		if (strcmp(argv[i], "-warn_on_unexpected_block") == 0)
		{
			fWarnForUnexpecedBlock = true;
			continue;
		}
		
		if (strcmp(argv[i], "-add_rootified_objects") == 0)
		{
			fAddCustomData = true;
			continue;
		}

		if (strcmp(argv[i], "-makeclonesarray") == 0)
		{
			fMakeClonesArray = true;
			continue;
		}

		if (strcmp(argv[i], "-makeonlyclonesarray") == 0)
		{
			fMakeMinimalESD = true;
			fMakeClonesArray = true;
			fMakeESDDataBlock = false;
			continue;
		}

		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}
	
	return 0;
}


int AliHLTMUONESDMaker::DoDeinit()
{
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	
	HLTInfo("Deinitialising dHLT ESD maker component.");
	return 0;
}


const char* AliHLTMUONESDMaker::GetComponentID()
{
	/// Inherited from AliHLTComponent. Returns the component ID.
	
	return AliHLTMUONConstants::ESDMakerId();
}


AliHLTComponentDataType AliHLTMUONESDMaker::GetOutputDataType()
{
	/// Inherited from AliHLTComponent. Returns kAliHLTMultipleDataType.
	/// Refer to GetOutputDataTypes for all returned data types.
	
	return kAliHLTMultipleDataType;
}


int AliHLTMUONESDMaker::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
	/// Inherited from AliHLTComponent.
	/// Returns the ESD object and kAliHLTDataTypeTObject data type with MUON origin.
	
	list.push_back(AliHLTMUONConstants::ESDDataType());
	list.push_back(kAliHLTDataTypeTObject | kAliHLTDataOriginMUON);
	return list.size();
}


void AliHLTMUONESDMaker::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
	/// Inherited from AliHLTProcessor.
	/// Returns the list of expected input data types.
	
	list.push_back(AliHLTMUONConstants::TriggerRecordsBlockDataType());
	list.push_back(AliHLTMUONConstants::MansoTracksBlockDataType());
	list.push_back(AliHLTMUONConstants::TracksBlockDataType());
}


void AliHLTMUONESDMaker::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	/// Inherited from AliHLTComponent.
	/// Returns an estimate of the expected output data size.
	
	constBase = sizeof(AliESDEvent) + 1024*1024;  // The extra 1 MByte is for auxilary objects created in AliESDEvent.
	inputMultiplier = 10;
}


AliHLTComponent* AliHLTMUONESDMaker::Spawn()
{
	/// Inherited from AliHLTComponent. Creates a new object instance.
	
	return new AliHLTMUONESDMaker();
}


int AliHLTMUONESDMaker::DoEvent(
		const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData
	)
{
	/// Inherited from AliHLTProcessor. Processes the new event data.
	
	AliESDEvent event;
	fClusterIndex = 0;  // for the cluster unique ID.
	
	// Create and fill in the standard ESD objects or just create the muon
	// tracks array if so requested.
	if (fMakeMinimalESD)
	{
		TClonesArray* muonTracks = new TClonesArray("AliESDMuonTrack",2);
		muonTracks->SetName("MuonTracks");
		event.AddObject(muonTracks);
		TClonesArray* muonClusters = new TClonesArray("AliESDMuonCluster",0);
		muonClusters->SetName("MuonClusters");
		event.AddObject(muonClusters);
		event.GetStdContent();
	}
	else
	{
		event.CreateStdContent();
		event.SetRunNumber(GetRunNo());
	}
	
	const AliHLTComponentBlockData* block = NULL;
	AliHLTUInt32_t specification = 0;  // Contains the output data block spec bits.
	std::vector<const AliHLTMUONTriggerRecordStruct*> triggerRecords;

	// First process the blocks of trigger records. We simply mark each trigger
	// record in the triggerRecords array.
	for (int i = 0; i < GetNumberOfInputBlocks(); i++)
	{
		block = GetInputBlock(i);
		assert( block != NULL );
		
		HLTDebug("Handling block: %u, with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
			i, DataType2Text(block->fDataType).c_str(), block->fPtr, block->fSize
		);
		
		if (block->fDataType == AliHLTMUONConstants::TriggerRecordsBlockDataType())
		{
			specification |= block->fSpecification;
			AliHLTMUONTriggerRecordsBlockReader inblock(block->fPtr, block->fSize);
			if (not BlockStructureOk(inblock))
			{
				if (DumpDataOnError()) DumpEvent(evtData, trigData);
				continue;
			}
			
			for (AliHLTUInt32_t n = 0; n < inblock.Nentries(); n++)
			{
				triggerRecords.push_back(&inblock[n]);
			}
		}
		else if (block->fDataType == AliHLTMUONConstants::RootifiedEventDataType() and fAddCustomData)
		{
			// Do nothing for now, will handle this later.
		}
		else
		{
			if (block->fDataType != AliHLTMUONConstants::MansoTracksBlockDataType() and
			    block->fDataType != AliHLTMUONConstants::TracksBlockDataType() )
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
	
	// If we were requested to add all dHLT rootified data objects then do so.
	if (fAddCustomData)
	{
		const AliHLTComponentDataType& type = AliHLTMUONConstants::RootifiedEventDataType();
		const char* classname = AliHLTMUONEvent::Class_Name();
		const TObject* obj = NULL;
		for (obj = GetFirstInputObject(type, classname); obj != NULL; obj = GetNextInputObject())
		{
			// Clone the object since the ESD takes ownership of it.
			event.AddObject(obj->Clone());
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
			AliESDMuonTrack *muTrack = event.NewMuonTrack();
			
			AliHLTMUONParticleSign sign;
			bool hitset[4];
			AliHLTMUONUtils::UnpackMansoTrackFlags(t.fFlags, sign, hitset);
			
			double signVal = 0;
			switch (sign)
			{
			case kSignMinus:   signVal = +1.; break;
			case kSignUnknown: signVal =  0.; break;
			case kSignPlus:    signVal = -1.; break;
			default:
				HLTWarning("Got a track with an invalid sign value: %d", int(sign));
			}
			
			TVector3 mom(t.fPx, t.fPy, t.fPz);
			if (mom.Mag() != 0)
				muTrack->SetInverseBendingMomentum(signVal/mom.Mag());
			else
				muTrack->SetInverseBendingMomentum(0.);
			muTrack->SetThetaX(atan2(t.fPx, t.fPz));
			muTrack->SetThetaY(atan2(t.fPy, t.fPz));
			muTrack->SetZ(0.);
			muTrack->SetBendingCoor(0.);
			muTrack->SetNonBendingCoor(0.);
			
			// The Manso algorithm assumes the information at the
			// Distance of Closest Approach and chamber 1 is the same
			// as the vertex.
			if (mom.Mag() != 0)
				muTrack->SetInverseBendingMomentumAtDCA(1./mom.Mag());
			else
				muTrack->SetInverseBendingMomentumAtDCA(0.);
			muTrack->SetThetaXAtDCA(atan2(t.fPx, t.fPz));
			muTrack->SetThetaYAtDCA(atan2(t.fPy, t.fPz));
			muTrack->SetBendingCoorAtDCA(0.);
			muTrack->SetNonBendingCoorAtDCA(0.);
			
			if (mom.Mag() != 0)
				muTrack->SetInverseBendingMomentumUncorrected(1./mom.Mag());
			else
				muTrack->SetInverseBendingMomentumUncorrected(0.);
			muTrack->SetThetaXUncorrected(atan2(t.fPx, t.fPz));
			muTrack->SetThetaYUncorrected(atan2(t.fPy, t.fPz));
			muTrack->SetZUncorrected(0.);
			muTrack->SetBendingCoorUncorrected(0.);
			muTrack->SetNonBendingCoorUncorrected(0.);
			
			muTrack->SetChi2(t.fChi2);
			
			// Fill in the track hit points.
			for (int i = 0; i < 4; i++)
			{
				if (not hitset[i]) continue;
				
				AliHLTUInt8_t chamber;
				AliHLTUInt16_t detElemId;
				AliHLTMUONUtils::UnpackRecHitFlags(t.fHit[i].fFlags, chamber, detElemId);
				
				AliESDMuonCluster *cluster = event.NewMuonCluster();
				cluster->SetUniqueID(AliMUONVCluster::BuildUniqueID(chamber, detElemId, fClusterIndex++));
				cluster->SetXYZ(t.fHit[i].fX, t.fHit[i].fY, t.fHit[i].fZ);
				cluster->SetErrXY(    // Use nominal values.
						  AliHLTMUONConstants::DefaultNonBendingReso(),
						  AliHLTMUONConstants::DefaultBendingReso()
						  );
				cluster->SetCharge(-1.);   // Indicate no total charge calculated.
				cluster->SetChi2(-1.);   // Indicate no fit made.
				muTrack->AddClusterId(cluster->GetUniqueID());
				muTrack->AddInMuonClusterMap(i+6);
			}
			
			FillTriggerInfo(triggerRecords, t.fTrigRec, t.fId, *muTrack, event);
		}
	}

	for (block = GetFirstInputBlock(AliHLTMUONConstants::TracksBlockDataType());
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		specification |= block->fSpecification;
		AliHLTMUONTracksBlockReader inblock(block->fPtr, block->fSize);
		if (not BlockStructureOk(inblock))
		{
			if (DumpDataOnError()) DumpEvent(evtData, trigData);
			continue;
		}
		
		for (AliHLTUInt32_t n = 0; n < inblock.Nentries(); n++)
		{
			const AliHLTMUONTrackStruct& t = inblock[n];
			AliESDMuonTrack *muTrack = event.NewMuonTrack();
			
			AliHLTMUONParticleSign sign;
			bool hitset[16];
			AliHLTMUONUtils::UnpackTrackFlags(t.fFlags, sign, hitset);
			
			muTrack->SetInverseBendingMomentum(t.fInverseBendingMomentum);
			muTrack->SetThetaX(t.fThetaX);
			muTrack->SetThetaY(t.fThetaY);
			muTrack->SetZ(t.fZ);
			muTrack->SetBendingCoor(t.fY);
			muTrack->SetNonBendingCoor(t.fX);
			
			// The full tracker assumes the information at the
			// Distance of Closest Approach and chamber 1 is the same
			// as the vertex.
			muTrack->SetInverseBendingMomentumAtDCA(t.fInverseBendingMomentum);
			muTrack->SetThetaXAtDCA(t.fThetaX);
			muTrack->SetThetaYAtDCA(t.fThetaY);
			muTrack->SetBendingCoorAtDCA(t.fY);
			muTrack->SetNonBendingCoorAtDCA(t.fX);
				
			muTrack->SetInverseBendingMomentumUncorrected(t.fInverseBendingMomentum);
			muTrack->SetThetaXUncorrected(t.fThetaX);
			muTrack->SetThetaYUncorrected(t.fThetaY);
			muTrack->SetZUncorrected(t.fZ);
			muTrack->SetBendingCoorUncorrected(t.fY);
			muTrack->SetNonBendingCoorUncorrected(t.fX);
				
			muTrack->SetChi2(t.fChi2);
			
			// Fill in the track hit points.
			for (int i = 0; i < 16; i++)
			{
				if (not hitset[i]) continue;
				
				AliHLTUInt8_t chamber;
				AliHLTUInt16_t detElemId;
				AliHLTMUONUtils::UnpackRecHitFlags(t.fHit[i].fFlags, chamber, detElemId);
				
				AliESDMuonCluster *cluster = event.NewMuonCluster();
				cluster->SetUniqueID(AliMUONVCluster::BuildUniqueID(chamber, detElemId, fClusterIndex++));
				cluster->SetXYZ(t.fHit[i].fX, t.fHit[i].fY, t.fHit[i].fZ);
				cluster->SetErrXY(    // Use nominal values.
						  AliHLTMUONConstants::DefaultNonBendingReso(),
						  AliHLTMUONConstants::DefaultBendingReso()
						  );
				cluster->SetCharge(-1.);   // Indicate no total charge calculated.
				cluster->SetChi2(-1.);   // Indicate no fit made.
				muTrack->AddClusterId(cluster->GetUniqueID());
				muTrack->AddInMuonClusterMap(chamber);
			}
			
			FillTriggerInfo(triggerRecords, t.fTrigRec, t.fId, *muTrack, event);
		}
	}

	if (fMakeClonesArray and event.GetList() != NULL)
	{
		PushBack(
			event.GetList()->FindObject("MuonTracks"),
			kAliHLTDataTypeTObject | kAliHLTDataOriginMUON,
			specification
		);
	        PushBack(
			 event.GetList()->FindObject("MuonClusters"),
			 kAliHLTDataTypeTObject | kAliHLTDataOriginMUON,
			 specification
			 );
	}
	if (fMakeESDDataBlock)
	{
		PushBack(&event, AliHLTMUONConstants::ESDDataType(), specification);
	}
	return 0;
}


void AliHLTMUONESDMaker::FillTriggerInfo(
		const AliTriggerRecordList& triggerRecords,
		AliHLTInt32_t trigRecId, AliHLTInt32_t trackId,
		AliESDMuonTrack& muTrack, AliESDEvent& event
	)
{
	/// Finds the trigger record with ID = 'id' and fills the output ESD track structure.
	
	// Find the corresponding trigger record.
	const AliHLTMUONTriggerRecordStruct* trigrec = NULL;
	for (size_t k = 0; k < triggerRecords.size(); k++)
	{
		if (triggerRecords[k]->fId == trigRecId)
		{
			trigrec = triggerRecords[k];
			break;
		}
	}
	// If the trigger record was found then fill its hit information also.
	if (trigrec != NULL)
	{
		AliHLTMUONParticleSign trsign;
		bool trhitset[4];
		AliHLTMUONUtils::UnpackTriggerRecordFlags(
				trigrec->fFlags, trsign, trhitset
			);
		
		for (int i = 0; i < 4; i++)
		{
			if (not trhitset[i]) continue;
			
			AliHLTUInt8_t chamber;
			AliHLTUInt16_t detElemId;
			AliHLTMUONUtils::UnpackRecHitFlags(trigrec->fHit[i].fFlags, chamber, detElemId);
		
			AliESDMuonCluster *cluster = event.NewMuonCluster();
			cluster->SetUniqueID(AliMUONVCluster::BuildUniqueID(chamber, detElemId, fClusterIndex++));
			cluster->SetXYZ(
					trigrec->fHit[i].fX,
					trigrec->fHit[i].fY,
					trigrec->fHit[i].fZ
					);
			cluster->SetErrXY(    // Use nominal values.
					  AliMUONConstants::TriggerNonBendingReso(),
					  AliMUONConstants::TriggerBendingReso()
					  );
			cluster->SetCharge(-1.);   // Indicate no total charge calculated.
			cluster->SetChi2(-1.);   // Indicate no fit made.
			muTrack.AddClusterId(cluster->GetUniqueID());
			muTrack.AddInMuonClusterMap(i+10);
		}
	}
	else
	{
		HLTWarning("Trigger record (ID = %d) not found for track ID = %d.",
			trigRecId, trackId
		);
	}
}
