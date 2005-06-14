////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "../src/TriggerRecord.hpp"


// Formated printing of floating point numbers.
void Print(float value)
{
	char buffer[1024];
	char* str = &buffer[0];
	sprintf(str, "%f", value);
	int length = strlen(str);
	if (str[0] == '-')
	{
		cout << str;
		for (int i = length; i < 15; i++)
			cout << " ";
	}
	else
	{
		cout << " " << str;
		for (int i = length+1; i < 15; i++)
			cout << " ";
	};
};


// TODO: fix the endian encoding of the data format.

/* Reads the contents of a file generated with MakeTriggerRecordFiles.C
   and prints the contents to screen.
 */
void DumpTriggerRecordFile(const char* filename)
{
	FILE* file = fopen(filename, "r");
	if (file == NULL)
	{
		Error("DumpTriggerRecordFile", "Could not open file: %s", filename);
		return;
	};
	
	// Read the first 32 bits which is the size of the data structure in
	// the file in 32bit words.
	dHLT::UInt size;
	fread(&size, sizeof(size), 1, file);
	if (ferror(file))
	{
		Error("DumpTriggerRecordFile", "Could not read from file: %s", filename);
		return;
	};
	cout << "Size of structure: " << size << " (4 byte) words." << endl;
	
	// Read the Trigger ID offset.
	dHLT::UInt triggerIDoffset;
	fread(&triggerIDoffset, sizeof(triggerIDoffset), 1, file);
	if (ferror(file))
	{
		Error("DumpTriggerRecordFile", "Could not read from file: %s", filename);
		return;
	};
	cout << "Trigger ID offset: " << triggerIDoffset << endl;
	
	cout << "Sign\t Pt             MT1 X          MT1 Y          MT2 X          MT2 Y" << endl;
	
	// Go through all the records and dump them to screen.
	UInt_t recordcount = (size - sizeof(dHLT::UInt)) / (sizeof(dHLT::TriggerRecord) / 4);
	for (UInt_t i = 0; i < recordcount; i++)
	{
		dHLT::TriggerRecord record;
		fread(&record, sizeof(record), 1, file);
		if (ferror(file))
		{
			Error("DumpTriggerRecordFile", "Could not read from file: %s", filename);
			return;
		};
		
		cout << record.sign << "\t";
		Print(record.pt);
		Print(record.station1impact.x);
		Print(record.station1impact.y);
		Print(record.station2impact.x);
		Print(record.station2impact.y);
		cout << endl;
	};
	
	fclose(file);
};

