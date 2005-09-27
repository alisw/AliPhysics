////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/TriggerSource.hpp"
#include "AliRoot/Base.hpp"
#include "Tracking/Calculations.hpp"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliModule.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONHit.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONDataInterface.h"
#include "TDatabasePDG.h"
#ifndef __alpha
#include <math.h>
#else
#include <float.h>
#endif

ClassImp(AliMUONHLT::TriggerSource)
ClassImp(AliMUONHLT::TriggerSource::EventData)

namespace AliMUONHLT
{

TriggerSource::TriggerSource()
	: TObject(), fEventList(TriggerSource::EventData::Class())
{
	fAreaToUse = FromWholePlane;
	fDataToUse = FromLocalTriggers;
	fMaxBlockSize = 0xFFFFFFFF;
	fUseLookupTable = kTRUE;
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fHadToLoadgAlice = kFALSE;
}


TriggerSource::TriggerSource(AliMUONDataInterface* data)
	: TObject(), fEventList(TriggerSource::EventData::Class())
{
	fAreaToUse = FromWholePlane;
	fDataToUse = FromLocalTriggers;
	fMaxBlockSize = 0xFFFFFFFF;
	fUseLookupTable = kTRUE;
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fHadToLoadgAlice = kFALSE;
	FillFrom(data);
}


TriggerSource::~TriggerSource()
{
	fEventList.Delete();
}


void TriggerSource::FillFrom(AliMUONDataInterface* data)
{
	DebugMsg(1, "FillFrom(AliMUONDataInterface*)");
	
	if (FileAndFolderOk(data))
	{
		AliMUON* module = NULL;
		if ( not FetchAliMUON(module) ) return;
		
		for (Int_t i = 0; i < data->NumberOfEvents(); i++)
		{
			AddEventFrom(data, module, i);
		}
		
		FinishedWithAliMUON();
	}
}


void TriggerSource::FillFrom(AliMUONDataInterface* data, Int_t event)
{
	DebugMsg(1, "FillFrom(AliMUONDataInterface*, Int_t)");
	
	if (FileAndFolderOk(data))
	{
		AliMUON* module = NULL;
		if ( not FetchAliMUON(module) ) return;
		AddEventFrom(data, module, event);
		FinishedWithAliMUON();
	}
}


void TriggerSource::FillFrom(
		AliMUONDataInterface* data,
		Int_t event, Int_t trigger, Bool_t newblock
	)
{
	DebugMsg(1, "FillFrom(AliMUONDataInterface*, Int_t, Int_t, Bool_t)");
	
	if (FileAndFolderOk(data))
	{
		data->GetEvent(event);
		AliMUON* module = NULL;
		if ( not FetchAliMUON(module) ) return;

		// Check if the current event corresponds to the event number we are
		// attempting to add to. If they do not or no event is selected then
		// try find the event or create a new one.
		if ( fCurrentEvent == NULL )
		{
			Bool_t found = GetEvent(event);
			if (not found) AddEvent(event);
		}
		else
		{
			if (fCurrentEvent->fEventNumber != event)
			{
				Bool_t found = GetEvent(event);
				if (not found) AddEvent(event);
			};
		};
		
		if ( fCurrentBlock != NULL )
		{
			Assert( fCurrentEvent != NULL );
			// If the newblock flag is set then force a new block.
			if (newblock) AddBlock();
		}
		else
			AddBlock();  // No block selected so we need to create a new block.

		AddTriggerFrom(data, module, trigger);
		FinishedWithAliMUON();
	}
}


void TriggerSource::Clear(Option_t* /*option*/)
{
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fEventList.Clear("C");
}


Bool_t TriggerSource::GetEvent(Int_t eventnumber) const
{
	DebugMsg(1, "TriggerSource::GetEvent(" << eventnumber << ")" );
	
	// Try find the corresponding event in the list of events.
	for (Int_t i = 0; i < fEventList.GetEntriesFast(); i++)
	{
		EventData* current = (EventData*) fEventList[i];
		if (current->fEventNumber == eventnumber)
		{
			fEventIndex = i;
			fCurrentEvent = current;
			GetFirstBlock();
			DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
				<< " , fTriggerIndex = " << fTriggerIndex
			);
			return kTRUE;
		};
	};
	return kFALSE;
}


Bool_t TriggerSource::GetFirstEvent() const
{
	DebugMsg(1, "TriggerSource::GetFirstEvent()");
	if (fEventList.GetEntriesFast() > 0)
	{
		fEventIndex = 0;
		fCurrentEvent = (EventData*) fEventList[0];
		GetFirstBlock();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kTRUE;
	}
	else
	{
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kFALSE;
	}
}


Bool_t TriggerSource::MoreEvents() const
{
	return 0 <= fEventIndex and fEventIndex < fEventList.GetEntriesFast();
}


Bool_t TriggerSource::GetNextEvent() const
{
	DebugMsg(1, "TriggerSource::GetNextEvent()");
	if (fEventIndex < fEventList.GetEntriesFast() - 1)
	{
		fCurrentEvent = (EventData*) fEventList[ ++fEventIndex ];
		GetFirstBlock();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kTRUE;
	}
	else
	{
		ResetAllPointers();
		return kFALSE;
	}
}


Int_t TriggerSource::CurrentEvent() const
{
	if (fCurrentEvent != NULL)
		return fCurrentEvent->fEventNumber;
	else
		return -1;
}


Int_t TriggerSource::NumberOfBlocks() const
{
	DebugMsg(1, "TriggerSource::NumberOfBlocks()");
	if (fCurrentEvent == NULL)
	{
		Error("NumberOfBlocks", "No event selected.");
		return -1;
	}
	else
		return fCurrentEvent->fBlocks.GetEntriesFast();
}


Bool_t TriggerSource::GetBlock(Int_t index) const
{
	DebugMsg(1, "TriggerSource::GetBlock(" << index << ")");
	
	// Note NumberOfBlocks() also checks if the event was selected.
	Int_t numberofblocks = NumberOfBlocks();
	if (numberofblocks < 0) return kFALSE;

	if ( 0 <= index and index < numberofblocks )
	{
		fBlockIndex = index;
		fCurrentBlock = (TClonesArray*) fCurrentEvent->fBlocks[index];
		GetFirstTrigger();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kTRUE;
	}
	else
	{
		// The index is out of bounds so inform the user.
		if (numberofblocks > 0)
			Error(	"GetBlock",
				"The block index (%d) is out of bounds. Valid range is [0, %d]",
				index, numberofblocks - 1
			);
		else
			Error(	"GetBlock",
				"The block index (%d) is out of bounds. No blocks found.",
				index
			);
		return kFALSE;
	}
}


Bool_t TriggerSource::GetFirstBlock() const
{
	DebugMsg(1, "TriggerSource::GetFirstBlock()");
	// Note: NumberOfBlocks() also checks if fCurrentEvent != NULL.
	if (NumberOfBlocks() > 0)
	{
		fBlockIndex = 0;
		fCurrentBlock = (TClonesArray*) fCurrentEvent->fBlocks[fBlockIndex];
		GetFirstTrigger();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kTRUE;
	}
	else
		return kFALSE;
}


Bool_t TriggerSource::MoreBlocks() const
{
	return 0 <= fBlockIndex and fBlockIndex < NumberOfBlocks();
}


Bool_t TriggerSource::GetNextBlock() const
{
	DebugMsg(1, "TriggerSource::GetNextBlock()");

	// Note: NumberOfBlocks() checks if fCurrentEvent != NULL. If it is then it returns -1
	// and since fBlockIndex is always >= -1 the if statement must go to the else part.
	if (fBlockIndex < NumberOfBlocks() - 1)
	{
		fCurrentBlock = (TClonesArray*) fCurrentEvent->fBlocks[ ++fBlockIndex ];
		GetFirstTrigger();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kTRUE;
	}
	else
	{
		ResetBlockPointers();
		return kFALSE;
	}
}


Int_t TriggerSource::NumberOfTriggers() const
{
	DebugMsg(1, "TriggerSource::NumberOfTriggers()");
	if (fCurrentBlock == NULL)
	{
		Error("NumberOfTriggers", "No block selected.");
		return -1;
	}
	else
		return fCurrentBlock->GetEntriesFast();
}


const TriggerRecord* TriggerSource::GetTrigger(Int_t triggernumber) const
{
	DebugMsg(1, "TriggerSource::GetTrigger(" << triggernumber << ")");

	if (fCurrentBlock == NULL)
	{
		Error("GetTrigger", "No block selected.");
		return NULL;
	};
	
	// Try find the corresponding trigger record in the list of events.
	for (Int_t i = 0; i < fCurrentBlock->GetEntriesFast(); i++)
	{
		TriggerRecord* current = (TriggerRecord*) fCurrentBlock->At(i);
		if (current->TriggerNumber() == triggernumber)
		{
			fTriggerIndex = i;
			fCurrentTrigger = current;
			DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
				<< " , fTriggerIndex = " << fTriggerIndex
			);
			return current;
		}
	}
	return NULL;
}


const TriggerRecord* TriggerSource::GetFirstTrigger() const
{
	DebugMsg(1, "TriggerSource::GetFirstTrigger()");
	// Note: NumberOfTriggers() also checks if fCurrentBlock != NULL.
	if (NumberOfTriggers() > 0)
	{
		fTriggerIndex = 0;
		fCurrentTrigger = (TriggerRecord*) fCurrentBlock->At(0);
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return fCurrentTrigger;
	}
	else
		return NULL;
};


Bool_t TriggerSource::MoreTriggers() const
{
	return 0 <= fTriggerIndex and fTriggerIndex < NumberOfTriggers();
}


const TriggerRecord* TriggerSource::GetNextTrigger() const
{
	DebugMsg(1, "TriggerSource::GetNextTrigger()");
	
	// Note: NumberOfTriggers() checks if fCurrentBlock != NULL. If it is then it returns -1
	// and since fTriggerIndex is always >= -1 the if statement must go to the else part.
	if (fTriggerIndex < NumberOfTriggers() - 1)
	{
		fCurrentTrigger = (TriggerRecord*) fCurrentBlock->At( ++fTriggerIndex );
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return fCurrentTrigger;
	}
	else
	{
		ResetTriggerPointers();
		return NULL;
	};
};


Int_t TriggerSource::CurrentTrigger() const
{
	if (fCurrentTrigger != NULL)
	{
		return fCurrentTrigger->TriggerNumber();
	}
	else
	{
		Error("CurrentTrigger", "No trigger record selected.");
		return -1;
	}
}


void TriggerSource::AddEvent(Int_t eventnumber)
{
	DebugMsg(1, "TriggerSource::AddEvent(" << eventnumber << ")");
	Assert( eventnumber >= 0 );

	// Assume the eventnumber does not already exist in the event list.
	fEventIndex = fEventList.GetEntriesFast();
	new ( fEventList[fEventIndex] ) EventData(eventnumber);
	fCurrentEvent = (EventData*) fEventList[fEventIndex];
	
	// Remember to reset the other pointers because the new event is empty.
	ResetBlockPointers();
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


void TriggerSource::AddBlock()
{
	DebugMsg(1, "TriggerSource::AddBlock()");
	
	if (fCurrentEvent == NULL)
	{
		Error("AddBlock", "No event selected.");
		return;
	};
	
	fBlockIndex = fCurrentEvent->fBlocks.GetEntriesFast();
	new ( fCurrentEvent->fBlocks[fBlockIndex] ) TClonesArray(TriggerRecord::Class());
	fCurrentBlock = (TClonesArray*) fCurrentEvent->fBlocks[fBlockIndex];
	
	// Remember to reset the trigger pointer because the new block is empty.
	ResetTriggerPointers();

	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


void TriggerSource::AddTrigger(const TriggerRecord& data)
{
	DebugMsg(1, "TriggerSource::AddTrigger(" << (void*)&data << ")");

	if (fCurrentBlock == NULL)
	{
		Error("AddTrigger", "No block selected.");
		return;
	}
	
	fTriggerIndex = fCurrentBlock->GetEntriesFast();
	new ( (*fCurrentBlock)[fTriggerIndex] ) TriggerRecord(data);
	fCurrentTrigger = (TriggerRecord*) (*fCurrentBlock)[fTriggerIndex];
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


Bool_t TriggerSource::FileAndFolderOk(AliMUONDataInterface* data)
{
	if (fFilename == "")
	{
		// Nothing filled yet so set the file and folder names.
		fFilename = data->CurrentFile();
		fFoldername = data->CurrentFolder();
		return kTRUE;
	};

	if ( fFilename != data->CurrentFile() )
	{
		Error(	"FileAndFolderOk",
			"The Trigger source already contains data from file '%s', cannot add data from file '%s'",
			fFilename.Data(), data->CurrentFile().Data()
		);
		return kFALSE;
	};
	
	if ( fFoldername != data->CurrentFolder() )
	{
		Error(	"FileAndFolderOk",
			"The Trigger source already contains data from folder '%s', cannot add data from folder '%s'",
			fFoldername.Data(), data->CurrentFolder().Data()
		);
		return kFALSE;
	};
	
	return kTRUE;
}


void TriggerSource::AddEventFrom(AliMUONDataInterface* data, AliMUON* module, Int_t event)
{
	if ( data->GetEvent(event) )
	{
		AddEvent(event);

		AddBlock();
		UInt_t currentblocksize = 0;
		TriggerRecord trigdata;

		switch (fDataToUse)
		{
		case FromHits:
			for (Int_t track = 0; track < data->NumberOfTracks(); track++)
			{
				if ( not FillTriggerFromHits(data, track, trigdata) )
					continue;  // Continue if unable to find hits.

				if (InFillRegion(trigdata))
				{
					AddTrigger(trigdata);

					// Create a new block if we reached the maximum block size.
					if ( ++currentblocksize == fMaxBlockSize )
					{
						AddBlock();
						currentblocksize = 0;
					};
				};
			};
			break;

		case FromLocalTriggers:
			Assert( module != NULL );
			DebugMsg(4, "Taking FromLocalTriggers branch...");
			for (Int_t i = 0; i < data->NumberOfLocalTriggers(); i++)
			{
				DebugMsg(4, "for loop: i = " << i);
				AliMUONLocalTrigger* lt = data->LocalTrigger(i);
				FillTriggerFromLocalTrigger(lt, module, trigdata);
				trigdata.TriggerNumber(i);

				if (InFillRegion(trigdata))
				{
					AddTrigger(trigdata);

					// Create a new block if we reached the maximum block size.
					if ( ++currentblocksize == fMaxBlockSize )
					{
						AddBlock();
						currentblocksize = 0;
					};
				};
			};
			break;

		default:
			Error("AddChamberFrom", "fDataToUse is not set to a valid value.");
		}
	}  // Loop on events.
}


void TriggerSource::AddTriggerFrom(AliMUONDataInterface* data, AliMUON* module, Int_t trigger)
{
	DebugMsg(1, "Entering AddTriggerFrom");

	TriggerRecord trigdata;

	switch (fDataToUse)
	{
	case FromHits:
		{
		// Note: in this case we treat the trigger parameter as a track number.
		if ( not FillTriggerFromHits(data, trigger, trigdata) )
			return;  // Continue if unable to find hits.
		}
		break;

	case FromLocalTriggers:
		{
		Assert( module != NULL );
		AliMUONLocalTrigger* lt = data->LocalTrigger(trigger);
		FillTriggerFromLocalTrigger(lt, module, trigdata);
		trigdata.TriggerNumber(trigger);
		}
		break;

	default:
		Error("AddTriggerFrom", "fDataToUse is not set to a valid value.");
		return;
	};
	
	AddTrigger(trigdata);
	
	DebugMsg(1, "Leaving AddTriggerFrom");
}


Bool_t TriggerSource::InFillRegion(const TriggerRecord& data)
{
	switch (fAreaToUse)
	{
	case FromWholePlane:     return kTRUE;
	case FromLeftHalfPlane:  return data.Station1Point().fX <= 0;
	case FromRightHalfPlane: return data.Station1Point().fX > 0;

	default:
		Error("InFillRegion", "fAreaToUse is not set to a valid value.");
		return kFALSE;
	}
}


void TriggerSource::FillTriggerFromLocalTrigger(
		AliMUONLocalTrigger* trigger, AliMUON* module, TriggerRecord& record
	)
{
	DebugMsg(2, "Creating TriggerRecord from AliMUONLocalTrigger object: " << (void*)trigger );
	AliMUONTriggerCircuit& circuit = module->TriggerCircuit(trigger->LoCircuit());

	// Get the sign of the particle the sign of the muon.
	if (trigger->LoLpt() == 1 || trigger->LoHpt() == 1 || trigger->LoApt() == 1)
	{
		record.ParticleSign(-1);
	}
	else
	if (trigger->LoLpt() == 2 || trigger->LoHpt() == 2 || trigger->LoApt() == 2)
	{
		record.ParticleSign(+1);
	}
	else
	{
		record.ParticleSign(0);
	};
	DebugMsg(2, "Particle sign = " << record.ParticleSign() );

	// Compute the transverse momentum.
	if (fUseLookupTable)
	{
		// TODO: implement use of the L0 lookup table.
		Error("FillTriggerFromLocalTrigger", "Use of L0 lookup table is not yet implemented!");
	}
	else
	{
		Float_t pt = circuit.PtCal( trigger->LoStripX(), trigger->LoDev(), trigger->LoStripY() );
		record.Pt(pt);
	};
	DebugMsg(2, "Pt = " << record.Pt() );

	// Build the impact points.
	record.Station1Point().fX = circuit.GetX11Pos(trigger->LoStripY());
	record.Station1Point().fY = circuit.GetY11Pos(trigger->LoStripX());
	record.Station2Point().fY = circuit.GetY21Pos(trigger->LoStripX() + trigger->LoDev() + 1);  // Why + 1?
	record.Station2Point().fX = module->Chamber(12).Z() * record.Station1Point().fX / module->Chamber(10).Z();
	DebugMsg(2, "fStation1x = " << record.Station1Point().fX);
	DebugMsg(2, "fStation1y = " << record.Station1Point().fY);
	DebugMsg(2, "fStation2x = " << record.Station2Point().fX);
	DebugMsg(2, "fStation2y = " << record.Station2Point().fY);
}


Bool_t TriggerSource::FillTriggerFromHits(AliMUONDataInterface* data, Int_t track, TriggerRecord& record)
{
	DebugMsg(2, "Creating TriggerRecord from hits on track: " << track );
	
	Float_t x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
#ifndef __alpha
#ifndef __sun
	x1 = y1 = z1 = x2 = y2 = z2 = x3 = y3 = z3 = x4 = y4 = z4 = NAN;
#else
	x1 = y1 = z1 = x2 = y2 = z2 = x3 = y3 = z3 = x4 = y4 = z4 = 0;
#endif
#else
	x1 = y1 = z1 = x2 = y2 = z2 = x3 = y3 = z3 = x4 = y4 = z4 = FLT_QNAN;
#endif
	// Find the hit that corresponds to chambers. 11 to 14. We can ignore any
	// hits above the first 14. If there are that many it means the particle
	// is cycling in the detector.
	for (Int_t i = 0; i < data->NumberOfHits(track) and i < 14; i++)
	{
		AliMUONHit* h = data->Hit(track, i);

		// Note AliMUONHit::Chamber() returns a value in the range 1..14
		// It is also important to have positive Z coordinates and under the
		// new version of Aliroot we use GEANT coordinates which return 
		// negatives. So we use the fabs routine.
		switch ( h->Chamber() )
		{
		case 11: x1 = h->X(); y1 = h->Y(); z1 = (Float_t)fabs(h->Z()); break;
		case 12: x2 = h->X(); y2 = h->Y(); z2 = (Float_t)fabs(h->Z()); break;
		case 13: x3 = h->X(); y3 = h->Y(); z3 = (Float_t)fabs(h->Z()); break;
		case 14: x4 = h->X(); y4 = h->Y(); z4 = (Float_t)fabs(h->Z()); break;
		};
	};
	DebugMsg(4, "Found: x1 = " << x1 << ", y1 = " << y1 << ", z1 = " << z1);
	DebugMsg(4, "Found: x2 = " << x2 << ", y2 = " << y2 << ", z2 = " << z2);
	DebugMsg(4, "Found: x3 = " << x3 << ", y3 = " << y3 << ", z3 = " << z3);
	DebugMsg(4, "Found: x4 = " << x4 << ", y4 = " << y4 << ", z4 = " << z4);

	// Get a coordinate for station 1, perferably from chamber 11 otherwise
	// use hits from chamber 12.
	if (not TMath::IsNaN(x1))
	{
		record.Station1Point().fX = x1;
		record.Station1Point().fY = y1;
		DebugMsg(3, "Using value from chamber 11: x1 = " << x1 << ", y1 = " << y1 << ", z1 = " << z1 );
	}
	else if (not TMath::IsNaN(x2))
	{
		record.Station1Point().fX = x2;
		record.Station1Point().fY = y2;
		z1 = z2;
		DebugMsg(3, "Using value from chamber 12: x2 = " << x2 << ", y2 = " << y2 << ", z2 = " << z2 );
	}
	else
	{
		// Return false if we could not find any hits on chambers 11 or 12.
		Warning("FillTriggerFromHits", "Could not find any hits on chambers 11 and 12.");
		return kFALSE;
	};

	// Get a coordinate for station 2, perferably from chamber 13 otherwise
	// use hits from chamber 14.
	if (not TMath::IsNaN(x3))
	{
		record.Station2Point().fX = x3;
		record.Station2Point().fY = y3;
		z2 = z3;
		DebugMsg(3, "Using value from chamber 13: x3 = " << x3 << ", y3 = " << y3 << ", z3 = " << z3 );
	}
	else if (not TMath::IsNaN(x4))
	{
		record.Station2Point().fX = x4;
		record.Station2Point().fY = y4;
		z2 = z4;
		DebugMsg(3, "Using value from chamber 14: x4 = " << x4 << ", y4 = " << y4 << ", z4 = " << z4 );
	}
	else
	{
		// Return false if we could not find any hits on chambers 13 or 14.
		Warning("FillTriggerFromHits", "Could not find any hits on chambers 13 and 14.");
		return kFALSE;
	};
	
	record.TriggerNumber(track);
	
	// Get the sign of the particle.
	Int_t particlecode = (Int_t) data->Hit(track, 0)->Particle();
	DebugMsg(3, "particle code = " << particlecode);
	TDatabasePDG* pdb = TDatabasePDG::Instance();
	TParticlePDG* pdata = pdb->GetParticle(particlecode);
	if (pdata->Charge() < 0)
		record.ParticleSign(-1);
	else if (pdata->Charge() > 0)
		record.ParticleSign(+1);
	else
		record.ParticleSign(0);
	DebugMsg(3, "Particle sign = " << record.ParticleSign());
	
	DebugMsg(3, "Calculating Pt: x1 = " << record.Station1Point().fX
			<< ", y1 = " << record.Station1Point().fY
			<< ", y2 = " << record.Station2Point().fY
			<< ", z1 = " << z1
			<< ", z2 = " << z2
		);
	// Calculate and assign the transverse momentum.
	Float_t pt = dHLT::Tracking::CalculatePt(
				record.Station1Point().fX,
				record.Station1Point().fY, record.Station2Point().fY,
				z1, z2
			);
	record.Pt(pt);
	
	DebugMsg(3, "Pt = " << record.Pt());
	
	return kTRUE;
}


Bool_t TriggerSource::FetchAliMUON(AliMUON*& module)
{
	// Check if we even need the MUON module. Not having to load it will
	// save a lot of loading time for AliRoot.
	if (fDataToUse == FromHits)
	{
		// Make sure we do not attempt to unload gAlice in FinishedWithAliMUON,
		// by setting the fHadToLoadgAlice to false.
		fHadToLoadgAlice = kFALSE;
		return kTRUE;
	};
	
	AliRunLoader* runloader = AliRunLoader::GetRunLoader();
	if ( runloader == NULL )
	{
		Error("FetchAliMUON", "AliRunLoader not initialised!");
		return kFALSE;
	};

	// Try fetch the AliRun object. If it is not found then try load it using
	// the runloader.
	AliRun* alirun = runloader->GetAliRun();
	if (alirun == NULL)
	{
		if (runloader->LoadgAlice() != 0)
		{
			// Error.
			DebugMsg(1, "Leaving FillFrom(AliMUONDataInterface*)");
			return kFALSE;
		};
		fHadToLoadgAlice = kTRUE;
		alirun = runloader->GetAliRun();
	}
	else
		fHadToLoadgAlice = kFALSE;

	// Get the MUON module pointer and return it.
	module = dynamic_cast<AliMUON*>( alirun->GetModule("MUON") );
	return kTRUE;
}


void TriggerSource::FinishedWithAliMUON()
{
	// Only unload the gAlice object if we had to load it ourselves.
	if (fHadToLoadgAlice)
		AliRunLoader::GetRunLoader()->UnloadgAlice();
}


void TriggerSource::ResetAllPointers() const
{
	fEventIndex = -1;
	fCurrentEvent = NULL;
	fBlockIndex = -1;
	fCurrentBlock = NULL;
	fTriggerIndex = -1;
	fCurrentTrigger = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


void TriggerSource::ResetBlockPointers() const
{
	fBlockIndex = -1;
	fCurrentBlock = NULL;
	fTriggerIndex = -1;
	fCurrentTrigger = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


void TriggerSource::ResetTriggerPointers() const
{
	fTriggerIndex = -1;
	fCurrentTrigger = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


TriggerSource::EventData::EventData() : fBlocks(TClonesArray::Class())
{
	fEventNumber = -1;
}


TriggerSource::EventData::EventData(Int_t eventnumber)
	: fBlocks(TClonesArray::Class())
{
	fEventNumber = eventnumber;

	// If the following is not set then we do not write the fBlocks properly.
	fBlocks.BypassStreamer(kFALSE);
}


TriggerSource::EventData::~EventData()
{
	//fBlocks.Clear("C");  // Done in fBlocks destructor 
}


} // AliMUONHLT
