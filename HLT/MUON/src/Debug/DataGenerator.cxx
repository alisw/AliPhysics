////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Debug/DataGenerator.hpp"
#include "RegionOfInterest.hpp"
#include "Error.hpp"

#include <vector>

namespace dHLT
{
namespace Debug
{


#define PI 3.14159265358979311599796346854418516159057

#define BINNING_GRID_SIZE 5


DataGenerator::DataGenerator()
{
	BL = 3.0;
	chamber_z[0] = 533.5;
	chamber_z[1] = 546.5;
	chamber_z[2] = 678.5;
	chamber_z[3] = 693.5;
	chamber_z[4] = 964.0;
	chamber_z[5] = 986.0;
	chamber_z[6] = 1251.5;
	chamber_z[7] = 1278.5;
	chamber_z[8] = 1416.5;
	chamber_z[9] = 1443.5;
	chamber_z[10] = 1603.5;
	chamber_z[11] = 1620.5;
	chamber_z[12] = 1703.5;
	chamber_z[13] = 1720.5;
	zf = 975.0;
};


void DataGenerator::GenerateData(
		const EventID eventid, const UInt number_of_tracks,
		TriggerSource& triggers, ClusterSource& clusters
	)
{
	typedef std::vector<ClusterPoint> ClusterList;
	ClusterList clusterlists[10][BINNING_GRID_SIZE][BINNING_GRID_SIZE];

	typedef std::vector<TriggerRecord> TriggerList;
	TriggerList triggerlists[BINNING_GRID_SIZE][BINNING_GRID_SIZE];

	for (UInt i = 0; i < number_of_tracks; i++)
	{
		ClusterPoint track[14];
		//Float theta = (Float)( Random() * 9.0 * PI / 180.0 );
		Float theta = (Float)( Random() * 4.0 * PI / 180.0 );
		Float phi = (Float)( Random() * 2.0 * PI );
		Float momentum = (Float)( Random() * 5.0 + 5.0 );
		Float charge = Random() > 0.5 ? 1.0f : -1.0f;
		ComputeTrack(theta, phi, momentum, charge, track);

		Int j;
		for (j = 0; j < 10; j++)
		{
			Double scale = chamber_z[j] * tan( 9.0 * PI / 180.0 );
			Int xi = (Int)( (track[j].x / (scale * 2) + 0.5) * BINNING_GRID_SIZE );
			if (xi < 0) xi = 0;
			if (xi > 9) xi = 9;
			Int yi = (Int)( (track[j].y / (scale * 2) + 0.5) * BINNING_GRID_SIZE );
			if (yi < 0) yi = 0;
			if (yi > 9) yi = 9;

			// Make sure the track coordinate is in the detection region before adding it.
			if ( track[j].x * track[j].x + track[j].y * track[j].y < scale * scale )
				clusterlists[j][xi][yi].push_back(track[j]);
		};
		
		j = 10;
		{
			Double scale = chamber_z[j] * tan( 9.0 * PI / 180.0 );
			Double scale12 = chamber_z[12] * tan( 9.0 * PI / 180.0 );
			Int xi = (Int)( (track[j].x / (scale * 2) + 0.5) * BINNING_GRID_SIZE );
			if (xi < 0) xi = 0;
			if (xi > 9) xi = 9;
			Int yi = (Int)( (track[j].y / (scale * 2) + 0.5) * BINNING_GRID_SIZE );
			if (yi < 0) yi = 0;
			if (yi > 9) yi = 9;

			TriggerRecord rec(charge < 0 ? Minus : Plus, momentum*sin(theta) + GausRandom(), track[10], track[12]);

			// Make sure the track coordinates are in the detection region before
			// adding the trigger record.
			if ( track[10].x * track[10].x + track[10].y * track[10].y < scale * scale
				and track[12].x * track[12].x + track[12].y * track[12].y < scale12 * scale12)
				triggerlists[xi][yi].push_back(rec);
		}
	};

	clusters.NewEvent(eventid);
	for (Int k = 0; k < 10; k++)
	{
		for (Int x = 0; x < BINNING_GRID_SIZE; x++)
		{
			for (Int y = 0; y < BINNING_GRID_SIZE; y++)
			{
				UInt count = clusterlists[k][x][y].size();
				if (count == 0) continue;

				// Create the region of interest of these cluster points.
				Assert( 0 < count and count <= clusterlists[k][x][y].size() );
				RegionOfInterest roi;
				roi.CreateToContain( clusterlists[k][x][y][0], (ChamberID)k );
				for (UInt n = 1; n < count; n++)
					roi.ExpandToContain( clusterlists[k][x][y][n] );

				clusters.NewClusterBlock( roi );
				for (UInt n = 0; n < count; n++)
					clusters.AddCluster( clusterlists[k][x][y][n] );
			};
		}
	};

	triggers.NewEvent(eventid);
	for (Int x = 0; x < BINNING_GRID_SIZE; x++)
	{
		for (Int y = 0; y < BINNING_GRID_SIZE; y++)
		{
			if (triggerlists[x][y].size() == 0) continue;
			
			triggers.NewTriggerBlock();
			for (UInt l = 0; l < triggerlists[x][y].size(); l++)
				triggers.AddTrigger( triggerlists[x][y][l] );
		};
	};
};


void DataGenerator::ComputeTrack(Float theta, Float phi, Float momentum, Float charge, ClusterPoint track[14])
{
	Double px = ((Double)momentum)*sin(theta)*cos(phi);
	Double py = ((Double)momentum)*sin(theta)*sin(phi);
	Double pz = ((Double)momentum)*cos(theta);
	Double pb = sqrt(py*py + pz*pz);

	// px / pz = tan(k)
	// sin(a0) = py / momentum;
	// 0.3 * charge * BL / pb = sin(a0) - sin(a);
	// sin(angle) = py / momentum - 0.3 * charge * BL / pb;
	// y = z * tan(angle)
	Double angle0 = asin(py / momentum);
	Double yf = zf * tan(angle0);
	Double angle = asin(py / momentum - 0.3 * charge * BL / pb);
	for (Int i = 0; i < 14; i++)
	{
		track[i].x = (Float)( px / pz * chamber_z[i] );
		if (chamber_z[i] < zf)
			track[i].y = (Float)( tan(angle0) * chamber_z[i] );
		else
			track[i].y = (Float)( yf + tan(angle) * chamber_z[i] );
	};
};


Float DataGenerator::Random()
{
	return (Float)rand() / (Float)RAND_MAX;
};


Float DataGenerator::GausRandom()
{
	// Using polar form of Box-Muller transformation for Gaussian distributed numbers.
	Float x1, x2, w, y1, y2;

	do {
		x1 = 2.0f * Random() - 1.0f;
		x2 = 2.0f * Random() - 1.0f;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0f );

	w = sqrt( (-2.0f * log( w ) ) / w );
	y1 = x1 * w;
	y2 = x2 * w;

	return y1;
};


} // Debug
} // dHLT
