// $Id: $
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

/// @file   AliHLTMuonSpectroTriggerComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   9 Nov 2009
/// @brief  Implementation of the muon spectrometer trigger component.
///
/// The AliHLTMuonSpectroTriggerComponent component is used to generate the HLT
/// trigger for the muon spectrometer.

#include "AliHLTMuonSpectroTriggerComponent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTMuonSpectroScalars.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliMUONTriggerDDLDecoderEventHandler.h"

ClassImp(AliHLTMuonSpectroTriggerComponent);


AliHLTMuonSpectroTriggerComponent::AliHLTMuonSpectroTriggerComponent() :
	AliHLTTrigger(),
	fBufferSizeConst(1024*16),
	fBufferSizeMultiplier(1),
	fMakeStats(false),
	fTriggerDDLs(false),
	fTriggerHits(false),
	fTriggerTrigRecs(false),
	fTriggerTracks(true),
	fTriggerDimuons(true)
{
	// Default constructor.
}


AliHLTMuonSpectroTriggerComponent::~AliHLTMuonSpectroTriggerComponent()
{
	// Default destructor.
}


void AliHLTMuonSpectroTriggerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list) const
{
	// Returns the list of input types expected.
	list.push_back(AliHLTMUONConstants::DDLRawDataType());	
	list.push_back(AliHLTMUONConstants::TriggerRecordsBlockDataType());
	list.push_back(AliHLTMUONConstants::RecHitsBlockDataType());
	list.push_back(AliHLTMUONConstants::MansoTracksBlockDataType());
	list.push_back(AliHLTMUONConstants::TracksBlockDataType());
	list.push_back(AliHLTMUONConstants::SinglesDecisionBlockDataType());
	list.push_back(AliHLTMUONConstants::PairsDecisionBlockDataType());
}


void AliHLTMuonSpectroTriggerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list) const
{
	// Return the output data types generated.
	
	list.push_back(kAliHLTDataTypeTriggerDecision);
	list.push_back(kAliHLTDataTypeEventStatistics|kAliHLTDataOriginHLT);
}


void AliHLTMuonSpectroTriggerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
	// Returns the output data size estimate.
	
	constBase = fBufferSizeConst;
	inputMultiplier = fBufferSizeMultiplier;
}


AliHLTComponent* AliHLTMuonSpectroTriggerComponent::Spawn()
{
	// Creates and returns a new instance.
	
	return new AliHLTMuonSpectroTriggerComponent;
}


Int_t AliHLTMuonSpectroTriggerComponent::DoInit(int argc, const char** argv)
{
	// Initialise the component.
	
	fMakeStats = false;
	fTriggerDDLs = false;
	fTriggerHits = false;
	fTriggerTrigRecs = false;
	fTriggerTracks = false;
	fTriggerDimuons = false;
 
	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-makestats") == 0)
		{
			fMakeStats = true;
			continue;
		}
		if (strcmp(argv[i], "-triggerddls") == 0)
		{
			fTriggerDDLs = true;
			continue;
		}
		if (strcmp(argv[i], "-triggerhits") == 0)
		{
			fTriggerHits = true;
			continue;
		}
		if (strcmp(argv[i], "-triggertrigrecs") == 0)
		{
			fTriggerTrigRecs = true;
			continue;
		}
		if (strcmp(argv[i], "-triggertracks") == 0)
		{
			fTriggerTracks = true;
			continue;
		}
		if (strcmp(argv[i], "-triggerdimuons") == 0)
		{
			fTriggerDimuons = true;
			continue;
		}
		if (strcmp(argv[i], "-triggerany") == 0)
		{
			fTriggerHits = true;
			fTriggerTrigRecs = true;
			fTriggerTracks = true;
			fTriggerDimuons = true;
			continue;
		}
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	} // for loop
	
	// If nothing was turned on to trigger then turn everything on by default.
	if (fTriggerHits == false and fTriggerTrigRecs == false
	    and fTriggerTracks == false and fTriggerDimuons == false
	   )
	{
		fTriggerTracks = true;
		fTriggerDimuons = true;
	}
	
	return 0;
}


Int_t AliHLTMuonSpectroTriggerComponent::DoDeinit()
{
	// Cleans up the component.
	
	return 0;
}


int AliHLTMuonSpectroTriggerComponent::DoTrigger()
{
	// Applies the trigger for the HLT.
	
	int result = 0;

	bool gotddls = false;
	bool gothits = false;
	bool gottrigrecs = false;
	bool gottracks = false;
	bool gotsingles = false;
	bool gotpairs = false;
	UInt_t nL0 = 0;
	UInt_t nhits = 0;
	UInt_t nhitsMTR = 0;
	UInt_t nhitsMCH = 0;
	UInt_t nhitsCh[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	UInt_t ntrigrecs = 0;
	UInt_t nL0plus = 0;
	UInt_t nL0minus = 0;
	UInt_t ntracksT = 0;  // from track blocks.
	UInt_t ntracksSD = 0; // from singles decision blocks
	UInt_t nplus = 0;
	UInt_t nminus = 0;
	UInt_t nlowpt = 0;
	UInt_t nhighpt = 0;
	Double_t minpt = -1;
	Double_t maxpt = -1;
	UInt_t nlikeany = 0;
	UInt_t nlikelow = 0;
	UInt_t nlikehigh = 0;
	UInt_t nunlikeany = 0;
	UInt_t nunlikelow = 0;
	UInt_t nunlikehigh = 0;
	UInt_t nlowmass = 0;
	UInt_t nhighmass = 0;
	Double_t minmass = -1;
	Double_t maxmass = -1;
	
	AliHLTComponentDataType blockType = AliHLTMUONConstants::DDLRawDataType();
	for (const AliHLTComponentBlockData* block = GetFirstInputBlock(blockType);
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		gotddls = true;
		nL0 = 1;
	}
	
	blockType = AliHLTMUONConstants::TriggerRecordsBlockDataType();
	for (const AliHLTComponentBlockData* block = GetFirstInputBlock(blockType);
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		AliHLTMUONTriggerRecordsBlockReader trigRecsBlock(block->fPtr, block->fSize);
		if (not IsBlockOk(trigRecsBlock, blockType))
		{
			HLTWarning("Skipping problematic block '%s'", DataType2Text(blockType).c_str());
			continue;
		}
		gothits = true;
		
		ntrigrecs += trigRecsBlock.Nentries();
		for (AliHLTUInt32_t i = 0; i < trigRecsBlock.Nentries(); ++i)
		{
			AliHLTMUONParticleSign sign;
			bool hitset[4];
			AliHLTMUONUtils::UnpackTriggerRecordFlags(trigRecsBlock[i].fFlags, sign, hitset);
			switch (sign)
			{
			case kSignPlus: ++nL0plus; break;
			case kSignMinus: ++nL0minus; break;
			default: break;
			}
			for (int j = 0; j < 4; ++j)
			{
				if (hitset[j])
				{
					++nhitsCh[j+10];
					++nhits;
					++nhitsMTR;
				}
			}
		}
	}
	
	blockType = AliHLTMUONConstants::RecHitsBlockDataType();
	for (const AliHLTComponentBlockData* block = GetFirstInputBlock(blockType);
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		AliHLTMUONRecHitsBlockReader hitsBlock(block->fPtr, block->fSize);
		if (not IsBlockOk(hitsBlock, blockType))
		{
			HLTWarning("Skipping problematic block '%s'", DataType2Text(blockType).c_str());
			continue;
		}
		gottrigrecs = true;
		
		nhits += hitsBlock.Nentries();
		nhitsMCH += hitsBlock.Nentries();
		for (AliHLTUInt32_t i = 0; i < hitsBlock.Nentries(); ++i)
		{
			AliHLTUInt8_t chamber;
			AliHLTUInt16_t detElemId;
			AliHLTMUONUtils::UnpackRecHitFlags(hitsBlock[i].fFlags, chamber, detElemId);
			if (chamber < 10)
			{
				++nhitsCh[chamber];
			}
			else
			{
				HLTWarning("Received a reconstructed hit which indicates"
					" an invalid chamber number of %d. The expected"
					" range is [0..9]. The data block is probably corrupt.",
					int(chamber)
				);
			}
		}
	}
	
	blockType = AliHLTMUONConstants::MansoTracksBlockDataType();
	for (const AliHLTComponentBlockData* block = GetFirstInputBlock(blockType);
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		AliHLTMUONMansoTracksBlockReader tracksBlock(block->fPtr, block->fSize);
		if (not IsBlockOk(tracksBlock, blockType))
		{
			HLTWarning("Skipping problematic block '%s'", DataType2Text(blockType).c_str());
			continue;
		}
		gottracks = true;
		
		ntracksT += tracksBlock.Nentries();
		for (AliHLTUInt32_t i = 0; i < tracksBlock.Nentries(); ++i)
		{
			AliHLTMUONParticleSign sign;
			bool hitset[4];
			AliHLTMUONUtils::UnpackMansoTrackFlags(tracksBlock[i].fFlags, sign, hitset);
			switch (sign)
			{
			case kSignPlus: ++nplus; break;
			case kSignMinus: ++nminus; break;
			default: break;
			}
		}
	}
	
	blockType = AliHLTMUONConstants::TracksBlockDataType();
	for (const AliHLTComponentBlockData* block = GetFirstInputBlock(blockType);
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		AliHLTMUONTracksBlockReader tracksBlock(block->fPtr, block->fSize);
		if (not IsBlockOk(tracksBlock, blockType))
		{
			HLTWarning("Skipping problematic block '%s'", DataType2Text(blockType).c_str());
			continue;
		}
		gottracks = true;
		
		ntracksT += tracksBlock.Nentries();
		for (AliHLTUInt32_t i = 0; i < tracksBlock.Nentries(); ++i)
		{
			AliHLTMUONParticleSign sign;
			bool hitset[16];
			AliHLTMUONUtils::UnpackTrackFlags(tracksBlock[i].fFlags, sign, hitset);
			switch (sign)
			{
			case kSignPlus: ++nplus; break;
			case kSignMinus: ++nminus; break;
			default: break;
			}
		}
	}
	
	blockType = AliHLTMUONConstants::SinglesDecisionBlockDataType();
	for (const AliHLTComponentBlockData* block = GetFirstInputBlock(blockType);
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		AliHLTMUONSinglesDecisionBlockReader singlesBlock(block->fPtr, block->fSize);
		if (not IsBlockOk(singlesBlock, blockType))
		{
			HLTWarning("Skipping problematic block '%s'", DataType2Text(blockType).c_str());
			continue;
		}
		gotsingles = true;
		
		ntracksSD += singlesBlock.Nentries();
		nlowpt += singlesBlock.BlockHeader().fNlowPt;
		nhighpt += singlesBlock.BlockHeader().fNhighPt;
		for (AliHLTUInt32_t i = 0; i < singlesBlock.Nentries(); ++i)
		{
			if (singlesBlock[i].fPt < minpt or minpt == -1) minpt = singlesBlock[i].fPt;
			if (singlesBlock[i].fPt > maxpt or maxpt == -1) maxpt = singlesBlock[i].fPt;
		}
	}
	
	blockType = AliHLTMUONConstants::PairsDecisionBlockDataType();
	for (const AliHLTComponentBlockData* block = GetFirstInputBlock(blockType);
	     block != NULL;
	     block = GetNextInputBlock()
	    )
	{
		AliHLTMUONPairsDecisionBlockReader pairsBlock(block->fPtr, block->fSize);
		if (not IsBlockOk(pairsBlock, blockType))
		{
			HLTWarning("Skipping problematic block '%s'", DataType2Text(blockType).c_str());
			continue;
		}
		gotpairs = true;
		
		nunlikeany += pairsBlock.BlockHeader().fNunlikeAnyPt;
		nlikeany += pairsBlock.BlockHeader().fNlikeAnyPt;
		// Dont use the other scalars from the pair decisions block header because
		// they are much more restrictive. They count pairs where both tracks pass
		// the low or high pT cut. But we want to relax this to just one track needs
		// to pass the pT cut.
		for (AliHLTUInt32_t i = 0; i < pairsBlock.Nentries(); ++i)
		{
			if (pairsBlock[i].fInvMass < minmass or minmass == -1) minmass = pairsBlock[i].fInvMass;
			if (pairsBlock[i].fInvMass > maxmass or maxmass == -1) maxmass = pairsBlock[i].fInvMass;
			bool highMass, lowMass, unlike;
			AliHLTUInt8_t highPtCount, lowPtCount;
			AliHLTMUONUtils::UnpackPairDecisionBits(
				pairsBlock[i].fTriggerBits, highMass, lowMass, unlike,
				highPtCount, lowPtCount
			);
			if (unlike)
			{
				if (lowPtCount >= 1) ++nunlikelow;
				if (highPtCount >= 1) ++nunlikehigh;
				if (lowMass) ++nlowmass;
				if (highMass) ++nhighmass;
			}
			else
			{
				if (lowPtCount >= 1) ++nlikelow;
				if (highPtCount >= 1) ++nlikehigh;
			}
		}
	}
	
	// Select the largest value for nTracks since we might only get this information
	// from singles decision blocks.
	UInt_t ntracks = ntracksSD > ntracksT ? ntracksSD : ntracksT;
	
	bool triggeredOnDDLs = fTriggerDDLs and nL0 > 0;
	bool triggeredOnHits = fTriggerHits and nhitsMCH > 0;
	bool triggeredOnTrigRecs = fTriggerTrigRecs and ntrigrecs > 0;
	bool triggeredOnTracks = fTriggerTracks and ntracks > 0;
	bool triggeredOnDimuons = fTriggerDimuons and (nunlikeany > 0
		or nunlikelow > 0 or nunlikehigh > 0 or nlowmass > 0 or nhighmass > 0);
	
	if (triggeredOnDimuons)
	{
		SetDescription("Dimuon in muon spectrometer");
		SetTriggerDomain(AliHLTTriggerDomain("*******:MUON"));
	}
	else if (triggeredOnTracks)
	{
		SetDescription("Tracks in muon spectrometer");
		SetTriggerDomain(AliHLTTriggerDomain("*******:MUON"));
	}
	else if (triggeredOnTrigRecs)
	{
		SetDescription("Muon trigger chambers triggered");
		if (triggeredOnHits)
		{
			SetReadoutList(AliHLTReadoutList(AliHLTReadoutList::kMUONTRG | AliHLTReadoutList::kMUONTRK));
			SetTriggerDomain(AliHLTTriggerDomain("TRIGRECS:MUON,RECHITS :MUON"));
		}
		else
		{
			SetReadoutList(AliHLTReadoutList(AliHLTReadoutList::kMUONTRG));
			SetTriggerDomain(AliHLTTriggerDomain("TRIGRECS:MUON"));
		}
	}
	else if (triggeredOnHits)
	{
		SetDescription("Hits in muon tracking chambers");
		SetReadoutList(AliHLTReadoutList(AliHLTReadoutList::kMUONTRK));
		SetTriggerDomain(AliHLTTriggerDomain("RECHITS :MUON"));
	}
	else if (triggeredOnDDLs)
	{
		SetDescription("DDL in muon tracking chambers");
		SetReadoutList(AliHLTReadoutList(AliHLTReadoutList::kMUONTRG | AliHLTReadoutList::kMUONTRK));
		SetTriggerDomain(AliHLTTriggerDomain("DDL_RAW :MUON"));
	}
	else
	{
		SetDescription("Not triggered");
		SetTriggerDomain(AliHLTTriggerDomain());
	}
	
	if (triggeredOnDimuons or triggeredOnTracks or triggeredOnTrigRecs or triggeredOnHits or triggeredOnDDLs)
	{
		result = TriggerEvent();
		if (result == -ENOSPC) goto increaseBuffer;
		if (result != 0) return result;
	}
	
	if (fMakeStats)
	{
		AliHLTMuonSpectroScalars scalars;
		if (gotddls) scalars.Add("NL0", "Number of L0 triggered event", nL0);
		if (gothits and gottrigrecs) scalars.Add("NHits", "Total number of hits", nhits);
		if (gottrigrecs) scalars.Add("NHitsMTR", "Number of hits in trigger chambers", nhitsMTR);
		if (gothits)
		{
			scalars.Add("NHitsMCH", "Number of hits in tracking chambers", nhitsMCH);
			for (int i = 0; i < 10; i++)
			{
				scalars.Add(Form("NHitsCh%d", i+1), Form("Number of hits in chamber %d", i+1), nhitsCh[i]);
			}
		}
		if (gottrigrecs)
		{
			for (int i = 10; i < 14; i++)
			{
				scalars.Add(Form("NHitsCh%d", i+1), Form("Number of hits in chamber %d", i+1), nhitsCh[i]);
			}
			scalars.Add("NTrigRecs", "Total number of trigger records", ntrigrecs);
			scalars.Add("NL0+", "Number of positive sign tracks in L0", nL0plus);
			scalars.Add("NL0-", "Number of negative sign tracks in L0", nL0minus);
		}
		if (gottracks or gotsingles) scalars.Add("NTracks", "Total number of tracks", ntracks);
		if (gottracks)
		{
			scalars.Add("N+", "Number of positive sign tracks", nplus);
			scalars.Add("N-", "Number of negative sign tracks", nminus);
		}
		if (gotsingles)
		{
			scalars.Add("NLowPt", "Number of low pT tracks", nlowpt);
			scalars.Add("NHighPt", "Number of high pT tracks", nhighpt);
			scalars.Add("MinPt", "Minimum pT found (GeV/c)", minpt);
			scalars.Add("MaxPt", "Maximum pT found (GeV/c)", maxpt);
		}
		if (gotpairs)
		{
			scalars.Add("NLikeAny", "Number of like sign track pairs", nlikeany);
			scalars.Add("NLikeLow", "Number of like sign pairs with at least 1 low pT track.", nlikelow);
			scalars.Add("NLikeHigh", "Number of like sign pairs with at least 1 high pT track.", nlikehigh);
			scalars.Add("NUnlikeAny", "Number of unlike sign track pairs", nunlikeany);
			scalars.Add("NUnlikeLow", "Number of unlike sign pairs with at least 1 low pT track.", nunlikelow);
			scalars.Add("NUnlikeHigh", "Number of unlike sign pairs with at least 1 high pT track.", nunlikehigh);
			scalars.Add("NLowMass", "Number of low mass dimuons", nlowmass);
			scalars.Add("NHighMass", "Number of high mass dimuons", nhighmass);
			scalars.Add("MinMass", "Minimum invariant mass found for dimuon (GeV/c^2)", minmass);
			scalars.Add("MaxMass", "Maximum invariant mass found for dimuon (GeV/c^2)", maxmass);
		}
		
		result = PushBack(&scalars, kAliHLTDataTypeEventStatistics|kAliHLTDataOriginHLT);
		if (result == -ENOSPC) goto increaseBuffer;
	}
	return result;
	
increaseBuffer:
	// Increase the estimated buffer space required since the PushBack
	// or TriggerEvent methods indicate that they ran out of buffer space.
	fBufferSizeConst += 1024*1024;
	fBufferSizeMultiplier *= 2.;
	return -ENOSPC;
}


template <typename BlockReader>
bool AliHLTMuonSpectroTriggerComponent::IsBlockOk(
		const BlockReader& reader, const AliHLTComponentDataType& type
	) const
{
	// Method for checking the block structure.
	
	if (not reader.BufferSizeOk())
	{
		string name = DataType2Text(type).c_str();
		size_t headerSize = sizeof(typename BlockReader::HeaderType);
		if (reader.BufferSize() < headerSize)
		{
			HLTError("Received a '%s' data block with a size of %d bytes,"
				" which is smaller than the minimum valid size of %d bytes."
				" The block must be corrupt.",
				name.c_str(), reader.BufferSize(), headerSize
			);
		}
		
		size_t expectedWidth = sizeof(typename BlockReader::ElementType);
		if (reader.CommonBlockHeader().fRecordWidth != expectedWidth)
		{
			HLTError("Received a '%s' data block with a record"
				" width of %d bytes, but the expected value is %d bytes."
				" The block might be corrupt.",
				name.c_str(),
				reader.CommonBlockHeader().fRecordWidth,
				expectedWidth
			);
		}
		
		HLTError("Received a '%s' data block with a size of %d bytes,"
			" but the block header claims the block should be %d bytes."
			" The block might be corrupt.",
			name.c_str(), reader.BufferSize(), reader.BytesUsed()
		);
		return false;
	}
	return true;
}
