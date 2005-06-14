////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Utils.hpp"
#include "Error.hpp"
#include "TriggerRecord.hpp"
#include "Cluster.hpp"
#include "Track.hpp"
#include "System/File.hpp"
#include "RegionOfInterest.hpp"

#include <endian.h>
// TODO: fix the endian encoding of the data format.

#include <stdlib.h>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

using namespace dHLT;


// Formated printing of floating point numbers.
void Print(Float value)
{
	char buffer[1024];
	char* str = &buffer[0];
	sprintf(str, "%f", value);
	Int length = strlen(str);
	if (str[0] == '-')
	{
		cout << str;
		for (Int i = length; i < 15; i++)
			cout << " ";
	}
	else
	{
		cout << " " << str;
		for (Int i = length+1; i < 15; i++)
			cout << " ";
	};
};


using System::File;


int PrintTriggerRecords(const char* filename)
{
	File file(filename, File::READ);
	
	// Read the first 32 bits which is the size of the data structure in
	// the file in 32bit words.
	UInt size;
	file.Read(&size, sizeof(size));
	
	// Check that the size value corresponds to the size of the file.
	if (size * 4 + sizeof(size) != file.Size())
	{
		cerr << "Error: The file is corrupt. Word count in file does not correspond to file size." << endl;
		return 3;
	};
	
	cout << "Size of structure: " << size << " (4 byte) words." << endl;
	
	UInt triggerIDoffset;
	file.Read(&triggerIDoffset, sizeof(triggerIDoffset));
	
	cout << "Trigger ID offset: " << triggerIDoffset << endl;
	cout << "Sign\t Pt             MT1 X          MT1 Y          MT2 X          MT2 Y" << endl;
	
	UInt recordcount = (size * 4 - sizeof(UInt)) / sizeof(TriggerRecord);
	
	// Check the record count value is correct.
	if ( (recordcount * sizeof(TriggerRecord) + sizeof(UInt) + 4) != file.Size() )
	{
		cerr << "Error: The file is corrupt or an unknown data format." << endl;
		return 3;
	};
	
	// Go through all the records and dump them to screen.
	for (UInt i = 0; i < recordcount; i++)
	{
		TriggerRecord record;
		file.Read(&record, sizeof(record));
		
		cout << record.sign << "\t";
		Print(record.pt);
		Print(record.station1impact.x);
		Print(record.station1impact.y);
		Print(record.station2impact.x);
		Print(record.station2impact.y);
		cout << endl;
	};
	
	return 0;
};


int PrintClusterPoints(const char* filename)
{
	File file(filename, File::READ);
	
	// Read the first 32 bits which is the size of the data structure in
	// the file in 32bit words.
	UInt size;
	file.Read(&size, sizeof(size));
	
	// Check that the size value corresponds to the size of the file.
	if (size * 4 + sizeof(size) != file.Size())
	{
		cerr << "Error: The file is corrupt. Word count in file does not correspond to file size." << endl;
		return 3;
	};
	
	cout << "Size of structure: " << size << " (4 byte) words." << endl;
	cout << " X              Y" << endl;
	
	UInt recordcount = size / (sizeof(ClusterPoint) / 4);
	
	// Check the record count value is correct.
	if ( (recordcount * sizeof(ClusterPoint) + 4) != file.Size() )
	{
		cerr << "Error: The file is corrupt or an unknown data format." << endl;
		return 3;
	};
	
	// Go through all the records and dump them to screen.
	for (UInt i = 0; i < recordcount; i++)
	{
		ClusterPoint point;
		file.Read(&point, sizeof(point));
		
		Print(point.x);
		Print(point.y);
		cout << endl;
	};
	
	return 0;
};


int PrintTracks(const char* filename)
{
	File file(filename, File::READ);
	
	// Read the first 32 bits which is the size of the data structure in
	// the file in 32bit words.
	UInt size;
	file.Read(&size, sizeof(size));

	// Check that the size value corresponds to the size of the file.
	if (size * 4 + sizeof(size) != file.Size())
	{
		cerr << "Error: The file is corrupt. Word count in file does not correspond to file size." << endl;
		return 3;
	};

	cout << "Size of structure: " << size << " (4 byte) words." << endl;
	
	UInt recordcount = size / (sizeof(Track) / 4);
	
	// Check the record count value is correct.
	if ( (recordcount * sizeof(Track) + 4) != file.Size() )
	{
		cerr << "Error: The file is corrupt or an unknown data format." << endl;
		return 3;
	};
	
	// Go through all the records and dump them to screen.
	for (UInt i = 0; i < recordcount; i++)
	{
		Track track;
		file.Read(&track, sizeof(track));
		
		cout << "Trigger ID: " << track.triggerid;
		switch (track.sign)
		{
		case Plus:  cout << "  Sign: +1  "; break;
		case Minus: cout << "  Sign: -1  "; break;
		default:    cout << "  Sign: X   "; break;
		};
		cout << "Momentum: ";
		Print(track.p);
		cout << "  Pt: ";
		Print(track.pt);
		cout << endl;
		
		cout << "\tChamber\t X              Y             ROI" << endl;
		for (UInt j = 0; j < 10; j++)
		{
			UInt chamber = j+1;
			
			// Check the ROI number if it is supposed to be valid.
			if ( track.region[j] != INVALID_ROI )
			{
				RegionOfInterest roi( track.region[j] );
				if ( roi.Chamber() != (ChamberID) j )
				{
					cerr << "Error: The file is corrupt. ROI number is invalid for chamber "
					<< chamber << ", got: " << roi.Chamber() << endl;
					return 3;
				};
			};
		
			cout << "\t" << chamber << "\t";
			Print(track.point[j].x);
			Print(track.point[j].y);
			printf("0x%8X", track.region[j]);
			cout << endl;
		};
	};
	
	return 0;
};


int main(const int argc, const char** argv)
{
	if (argc != 2)
	{
		cerr << "Usage: > printrawdata <filename>" << endl;
		cerr << "Where <filename> is the name of a file containing:" << endl;
		cerr << "    TriggerRecord data *.trigrec" << endl;
		cerr << "    ClusterPoint data *.clusters" << endl;
		cerr << "    Track data *.tracks" << endl;
		return 1;
	};
	
	int returncode = 0;
	
	try
	{
		// Find the first dot.
		int length = strlen(argv[1]);
		int i;
		for (i = length; i >= 0; i--)
			if (argv[1][i] == '.') break;
		if (argv[1][i] == '.')
		{
			if ( strcmp(&argv[1][i], ".trigrec") == 0 )
			{
				returncode = PrintTriggerRecords(argv[1]);
			}
			else if ( strcmp(&argv[1][i], ".clusters") == 0 )
			{
				returncode = PrintClusterPoints(argv[1]);
			}
			else if ( strcmp(&argv[1][i], ".tracks") == 0 )
			{
				returncode = PrintTracks(argv[1]);
			}
			else
			{
				cerr << "Error: File extension not recognised." << endl;
				returncode = 2;
			};
		}
		else
		{
			cerr << "Error: File extension not recognised." << endl;
			returncode = 2;
		};
	}
	catch (const Error& e)
	{
		cerr << "Error: " << e << endl;
		returncode = 2;
	}
	catch (...)
	{
		cerr << "Error: Unknown exception!" << endl;
		returncode = 2;
	};

	return returncode;
};

