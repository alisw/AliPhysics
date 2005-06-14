////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "../src/Cluster.hpp"


// TODO: fix the endian encoding of the data format.

/* Generates .clusters binary files with the following structure:
      
     Byte #   Field

   // Header
         00   size              - Size of the structure in 32-bit words (excluding this word).
      
   // First ClusterPoint structure
         04   x                 - X coordinate of point.
         08   y                 - Y coordinate of point.

   // Next ClusterPoint structure
         08   x                 - X coordinate of point.
         12   y                 - Y coordinate of point.

   // etc...
 */
void MakeClusterPointFiles(const char* outputpath = ".")
{
	//gSystem->Load("../lib/Linux-debug/libMUONHLT.so");
	gSystem->Load("../lib/Linux/libMUONHLT.so");

	using namespace AliMUONHLT;
	
	AliMUONDataInterface data;
	data.SetFile();
	
	// Load the cluster data.
	AliMUONHLT::ClusterSource cs;
	cs.DataToUse(AliMUONHLT::ClusterSource::FromHits);
	cs.FillFrom(&data);
	
	for (Int_t event = 0; event < cs.NumberOfEvents(); event++)
	{
		cs.GetEvent(event);
		
		for (UInt_t chamber = 0; chamber < 10; chamber++)
		{
			// For every event create a new file.
			char buffer[1024];
			char* filename = &buffer[0];
			sprintf(filename, "%s/event%d_chamber%d.clusters", outputpath, event, chamber+1);
			FILE* file = fopen(filename, "w+b");
			if (file == NULL)
			{
				Error("MakeClusterPointFiles", "Could not create/open file: %s", filename);
				return;
			};

			// Compute and write the size of the whole data structure in 32bit words:
			dHLT::UInt size = 0;
			for (Int_t block = 0; block < cs.NumberOfBlocks(); block++)
			{
				cs.GetBlock(block);
				if (cs.Chamber() != chamber) continue;
				size += (cs.NumberOfClusters() * sizeof(dHLT::ClusterPoint)) / 4;
			};

			fwrite(&size, sizeof(dHLT::UInt), 1, file);
			if (ferror(file))
			{
				Error("MakeClusterPointFiles", "Could not write to file: %s", filename);
				return;
			};

			// Write all the TriggerRecord structures into the file.
			for (Int_t block = 0; block < cs.NumberOfBlocks(); block++)
			{
				cs.GetBlock(block);
				if (cs.Chamber() != chamber) continue;

				for (	const Point* cluster = cs.GetFirstCluster();
					cs.MoreClusters();
					cluster = cs.GetNextCluster()
				)
				{
					dHLT::ClusterPoint point;
					point.x = cluster->fX;
					point.y = cluster->fY;

					fwrite(&point, sizeof(point), 1, file);
					if (ferror(file))
					{
						Error("MakeClusterPointFiles", "Could not write to file: %s", filename);
						return;
					};
				};
			};

			fclose(file);
		};
	};
};

