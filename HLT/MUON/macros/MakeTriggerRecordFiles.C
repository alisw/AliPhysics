////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "../src/TriggerRecord.hpp"


// TODO: fix the endian encoding of the data format.

/* Generates .trigrec binary files with the following structure:
      
     Byte #   Field

   // Header
         00   size              - Size of the structure in 32-bit words (excluding this word).

         04   Trigger ID offset - The offset to add to the trigger record index number
                                  to get full trigger ID of any trigger record.
      
   // First TriggerRecord structure.
         08   sign              - Signed integer, possible values: -1, 0, +1
         12   pt                - Transverse momentum as a 32bit floating point number.
         16   station1impact.x  - 
         20   station1impact.y  - 
         24   station2impact.x  -
         28   station2impact.y  -

   // Next TriggerRecord structure.
   // ...
 */
void MakeTriggerRecordFiles(const char* outputpath = ".")
{
	//gSystem->Load("../lib/Linux-debug/libMUONHLT.so");
	gSystem->Load("../lib/Linux/libMUONHLT.so");

	using namespace AliMUONHLT;
	
	AliMUONDataInterface data;
	data.SetFile();
	
	// Load the trigger data.
	AliMUONHLT::TriggerSource ts;
	ts.DataToUse(AliMUONHLT::TriggerSource::FromHits);
	ts.FillFrom(&data);
	
	dHLT::UInt triggerIDoffset;
	
	for (Int_t event = 0; event < ts.NumberOfEvents(); event++)
	{
		ts.GetEvent(event);
		
		// For every event create a new file.
		char buffer[1024];
		char* filename = &buffer[0];
		sprintf(filename, "%s/event%d.trigrec", outputpath, event);
		FILE* file = fopen(filename, "w+b");
		if (file == NULL)
		{
			Error("MakeTriggerRecordFiles", "Could not create/open file: %s", filename);
			return;
		};
		
		// Compute and write the size of the whole data structure in 32bit words:
		dHLT::UInt size = 0;
		for (Int_t block = 0; block < ts.NumberOfBlocks(); block++)
		{
			ts.GetBlock(block);
			size += (ts.NumberOfTriggers() * sizeof(dHLT::TriggerRecord)) / 4;
		};
		size += sizeof(triggerIDoffset) / 4; // Add word count of triggerIDoffset.
		
		fwrite(&size, sizeof(size), 1, file);
		if (ferror(file))
		{
			Error("MakeTriggerRecordFiles", "Could not write to file: %s", filename);
			return;
		};
		
		// Compute and write the trigger ID offset.
		triggerIDoffset = 0;  // very simple case.
		fwrite(&triggerIDoffset, sizeof(triggerIDoffset), 1, file);
		if (ferror(file))
		{
			Error("MakeTriggerRecordFiles", "Could not write to file: %s", filename);
			return;
		};
		
		// Write all the TriggerRecord structures into the file.
		for (Int_t block = 0; block < ts.NumberOfBlocks(); block++)
		{
			ts.GetBlock(block);
			for (	const TriggerRecord* trigger = ts.GetFirstTrigger();
				ts.MoreTriggers();
				trigger = ts.GetNextTrigger()
			    )
			{
				dHLT::TriggerRecord record;
				record.sign = trigger->ParticleSign();
				record.pt = trigger->Pt();
				record.station1impact.x = trigger->Station1Point().fX;
				record.station1impact.y = trigger->Station1Point().fY;
				record.station2impact.x = trigger->Station2Point().fX;
				record.station2impact.y = trigger->Station2Point().fY;
				
				fwrite(&record, sizeof(record), 1, file);
				if (ferror(file))
				{
					Error("MakeTriggerRecordFiles", "Could not write to file: %s", filename);
					return;
				};
			};
		};
		
		fclose(file);
	};
};

