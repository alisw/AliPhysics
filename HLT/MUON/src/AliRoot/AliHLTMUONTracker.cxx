// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com

#include "AliHLTMUONTracker.h"
#include "Tracking/MansoTracker.hpp"
#include "AliMUONDataInterface.h"
#include "AliRunLoader.h"
#include "AliLog.h"
#include "AliESD.h"
#include "AliESDMuonTrack.h"

ClassImp(AliHLTMUONTracker)


AliHLTMUONTracker::AliHLTMUONTracker(AliRunLoader* runloader) :
	AliTracker(), fdHLT(NULL), fTriggers(NULL), fClusters(NULL), fTracks(NULL)
{
// Creates the the AliHLTMUONMicrodHLT object and its associated data source and sink
// objects. The AliHLTMUONMicrodHLT object is then initialised by hooking to these objects.

	AliDebug(2, Form("Called for object 0x%X and with runloader = 0x%X", (ULong_t)this, (ULong_t)runloader));

	// Create the dHLT objects.
	fdHLT = new AliHLTMUONMicrodHLT();
	fTriggers = new AliHLTMUONTriggerSource();
	fClusters = new AliHLTMUONClusterSource();
	fTracks = new AliHLTMUONTrackSink();

	// Hook up all the objects.
	fdHLT->SetTriggerSource(fTriggers);
	fdHLT->SetClusterSource(fClusters);
	fdHLT->SetTrackSink(fTracks);
}


AliHLTMUONTracker::~AliHLTMUONTracker()
{
// Deletes all dHLT objects created in the constructor.

	AliDebug(2, Form("Called for object 0x%X", (ULong_t)this));
	delete fTracks;
	delete fClusters;
	delete fTriggers;
	delete fdHLT;
}


Int_t AliHLTMUONTracker::LoadClusters(TTree* data)
{
// Fills the trigger and cluster data source from the current runloader.
// The runloader must be open and initialised properly.

	AliDebug(1, Form("Called for object 0x%X and with data tree = 0x%X", (ULong_t)this, (ULong_t)data));

	// Initialise the MUON data interface to use current runloader.
	AliMUONDataInterface di;
	di.UseCurrentRunLoader();
	AliDebug(2, Form("Loading for event %d", di.CurrentEvent()));

	// Load the trigger records.
	fTriggers->DataToUse(AliHLTMUONTriggerSource::kFromLocalTriggers);
	fTriggers->FillFrom(&di, di.CurrentEvent());

#ifndef LOG_NO_DEBUG
	TString str = "====================== Loaded Triggers ========================\nTrig#\tSign\tPt\t\tX1\t\tY1\t\tX2\t\tY2\n";
	for (fTriggers->GetFirstEvent(); fTriggers->MoreEvents(); fTriggers->GetNextEvent())
	for (fTriggers->GetFirstBlock(); fTriggers->MoreBlocks(); fTriggers->GetNextBlock())
	for (fTriggers->GetFirstTrigger(); fTriggers->MoreTriggers(); fTriggers->GetNextTrigger())
	{
		const AliHLTMUONTriggerRecord* trig = fTriggers->GetTrigger();
		str += trig->TriggerNumber();
		str += "\t";
		str += trig->ParticleSign();
		str += "\t";
		str += trig->Pt();
		str += "\t";
		str += trig->Station1Point().X();
		str += "\t";
		str += trig->Station1Point().Y();
		str += "\t";
		str += trig->Station2Point().X();
		str += "\t";
		str += trig->Station2Point().Y();
		str += "\n";
	}
	AliDebug(5, str);
#endif // LOG_NO_DEBUG

	// Load cluster points (reconstructed hits)
	fClusters->DataToUse(AliHLTMUONClusterSource::kFromRawClusters);
	fClusters->FillFrom(&di, di.CurrentEvent());

#ifndef LOG_NO_DEBUG
	str = "====================== Loaded Clusters ========================\nChamber\tX\t\tY\n";
	for (fClusters->GetFirstEvent(); fClusters->MoreEvents(); fClusters->GetNextEvent())
	for (fClusters->GetFirstBlock(); fClusters->MoreBlocks(); fClusters->GetNextBlock())
	for (fClusters->GetFirstCluster(); fClusters->MoreClusters(); fClusters->GetNextCluster())
	{
		Float_t x, y;
		fClusters->FetchCluster(x, y);
		str += fClusters->Chamber();
		str += "\t";
		str += x;
		str += "\t";
		str += y;
		str += "\n";
	}
	AliDebug(5, str);
#endif // LOG_NO_DEBUG

	return 0;
}


void AliHLTMUONTracker::UnloadClusters()
{
// Frees the triggers and clusters loaded in LoadClusters().

	AliDebug(1, Form("Called for object 0x%X", (ULong_t)this));

	// Release internal arrays.
	fTriggers->Clear();
	fClusters->Clear();

	return;
}


const AliHLTMUONTriggerRecord*
AliHLTMUONTracker::FindTriggerRecord(const AliHLTMUONTrack* track) const
{
// Finds the corresponding trigger record object for the given track.
// This is doen by matching up the trigger record number with the trigger
// record ID number in the track.

	fTriggers->GetEvent( fTracks->CurrentEvent() );
	for (fTriggers->GetFirstBlock(); fTriggers->MoreBlocks(); fTriggers->GetNextBlock())
	for (fTriggers->GetFirstTrigger(); fTriggers->MoreTriggers(); fTriggers->GetNextTrigger())
	{
		const AliHLTMUONTriggerRecord* trigrec = fTriggers->GetTrigger();
		if ( trigrec->TriggerNumber() == track->TriggerID() )
			return trigrec;
	}
	return NULL;
}


void AliHLTMUONTracker::LeastSquaresFit(const Double_t x[4], const Double_t y[4], Double_t& m, Double_t& c) const
{
// Least squares fit for a 4 point line: y = m*x + c

	Double_t n = 4.;  // Number of data points.

	// Compute the sums.
	Double_t sumx, sumy, sumxx, sumxy;
	sumx = sumy = sumxx = sumxy = 0.;
	for (Int_t i = 0; i < 4; i++)
	{
		sumx += x[i];
		sumy += y[i];
		sumxx += x[i] * x[i];
		sumxy += x[i] * y[i];
	}
	// Now compute the denominator then m and c parameters.
	Double_t denom = n*sumxx - sumx*sumx;
	m = (n*sumxy - sumx*sumy) / denom;
	c = (sumy*sumxx - sumx*sumxy) / denom;
}


Double_t AliHLTMUONTracker::ComputeChi2(const AliHLTMUONTrack* track) const
{
// Computes the Chi^2 for the line fit of the found track fragment.

	const AliHLTMUONTriggerRecord* trigrec = FindTriggerRecord(track);
	if (trigrec == NULL) return -1.;
	AliDebug(10, Form("Found trigger #%d, particle sign = %d, pt = %f",
			trigrec->TriggerNumber(),
			trigrec->ParticleSign(),
			trigrec->Pt()
		)
	);

	// Initialise X, Y and Z coordinate arrays.	
	Double_t st4x, st4y, st4z, st5x, st5y, st5z;
	if (track->Hit(6).X() == 0. && track->Hit(6).Y() == 0.)
	{
		st4x = track->Hit(7).X();
		st4y = track->Hit(7).Y();
		st4z = AliHLTMUONCoreMansoTracker::GetZ8();
	}
	else
	{
		st4x = track->Hit(6).X();
		st4y = track->Hit(6).Y();
		st4z = AliHLTMUONCoreMansoTracker::GetZ7();
	}
	if (track->Hit(8).X() == 0. && track->Hit(8).Y() == 0.)
	{
		st5x = track->Hit(9).X();
		st5y = track->Hit(9).Y();
		st5z = AliHLTMUONCoreMansoTracker::GetZ10();
	}
	else
	{
		st5x = track->Hit(8).X();
		st5y = track->Hit(8).Y();
		st5z = AliHLTMUONCoreMansoTracker::GetZ9();
	}

	Double_t x[4] = {st4x, st5x,
			trigrec->Station1Point().X(), trigrec->Station2Point().X()
		};
	Double_t y[4] = {st4y, st5y,
			trigrec->Station1Point().Y(), trigrec->Station2Point().Y()
		};
	Double_t z[4] = {st4z, st5z,
			AliHLTMUONCoreMansoTracker::GetZ11(),
			AliHLTMUONCoreMansoTracker::GetZ13()
		};

	// Fit a line to the x, y, z data points.
	Double_t mx, cx;  // for x = mx*z + cx
	Double_t my, cy;  // for y = my*z + cy 
	LeastSquaresFit(z, x, mx, cx);
	AliDebug(11, Form("Fitted line to xz plane: m = %f c = %f", mx, cx));
	LeastSquaresFit(z, y, my, cy);
	AliDebug(11, Form("Fitted line to yz plane: m = %f c = %f", my, cy));

	// Now compute the residual, square it and add it to the total chi2.
	Double_t chi2 = 0.;
	for (Int_t i = 0; i < 4; i++)
	{
		Double_t rx = x[i] - (mx*z[i] + cx);
		Double_t ry = y[i] - (my*z[i] + cy);
		AliDebug(15, Form("Real x = %f, fitted x = %f for z = %f", x[i], (mx*z[i] + cx), z[i]));
		AliDebug(15, Form("Real y = %f, fitted y = %f for z = %f", y[i], (my*z[i] + cy), z[i]));
		chi2 += rx*rx + ry*ry;
		AliDebug(15, Form("Running total for Chi2 = %f", chi2));
	}

	return chi2;
}


Int_t AliHLTMUONTracker::Clusters2Tracks(AliESD* event)
{
// Runs the dHLT tracking algorithm via the AliHLTMUONMicrodHLT fdHLT object. The found tracks are
// filled by AliHLTMUONMicrodHLT in fTracks, which then need to be unpacked into the AliRoot ESD.

	AliDebug(1, Form("Called for object 0x%X and ESD object = 0x%X", (ULong_t)this, (ULong_t)event));

	// Run the dHLT tracking algorithm.
	fdHLT->Run();
	AliDebug(1, "Finished running dHLT.");

	// Unpack the output and fill the ESD.
#ifndef LOG_NO_DEBUG
	TString str = "====================== Found Tracks ========================\nID\tSign\tP\t\tPt\t\tChamber\tX\t\tY\n";
#endif // LOG_NO_DEBUG
	for (fTracks->GetFirstEvent(); fTracks->MoreEvents(); fTracks->GetNextEvent())
	for (fTracks->GetFirstBlock(); fTracks->MoreBlocks(); fTracks->GetNextBlock())
	for (fTracks->GetFirstTrack(); fTracks->MoreTracks(); fTracks->GetNextTrack())
	{
		const AliHLTMUONTrack* track = fTracks->GetTrack();

#ifndef LOG_NO_DEBUG
		// Build some debug logging information.
		str += track->TriggerID();
		str += "\t";
		str += track->ParticleSign();
		str += "\t";
		str += track->P();
		str += "\t";
		str += track->Pt();
		str += "\t7\t";
		str += track->Hit(6).X();
		str += "\t";
		str += track->Hit(6).Y();
		str += "\n\t\t\t\t\t\t8\t";
		str += track->Hit(7).X();
		str += "\t";
		str += track->Hit(7).Y();
		str += "\n\t\t\t\t\t\t9\t";
		str += track->Hit(8).X();
		str += "\t";
		str += track->Hit(8).Y();
		str += "\n\t\t\t\t\t\t10\t";
		str += track->Hit(9).X();
		str += "\t";
		str += track->Hit(9).Y();
		str += "\n";
#endif // LOG_NO_DEBUG

		Double_t z = AliHLTMUONCoreMansoTracker::GetZ7();
		Double_t x = track->Hit(6).X();
		if (track->Hit(6).X() == 0. && track->Hit(6).Y() == 0.)
		{
			z = AliHLTMUONCoreMansoTracker::GetZ8();
			x = track->Hit(7).X();
		};

		Double_t p = track->P();
		Double_t pt = track->Pt();

		// Using the approximation: px / pz = x / z , and pz^2 = p^2 - pt^2
		Double_t pz2 = TMath::Abs(p*p - pt*pt);
		Double_t px2 = -1.;
		if (z != 0.)
			px2 = pz2 * x*x / (z*z);
		Double_t pyz = -1.;
		if (p*p - px2 >= 0.)
			pyz = TMath::Sqrt(p*p - px2);
		Double_t px = -1.;
		if (px2 >= 0.)
			px = TMath::Sqrt(px2);
		Double_t py = -1.;
		if (pt*pt - px2 >= 0.)
			py = TMath::Sqrt(pt*pt - px2);   // pt^2 = px^2 + py^2
		Double_t pz = -1.;
		if (pz2 >= 0.)
			pz = TMath::Sqrt(pz2);

		Double_t invmom = -1.0;
		if (pyz != 0.)
			invmom = 1. / pyz * track->ParticleSign();
		Double_t thetax = TMath::ATan2(px, pz);
		Double_t thetay = TMath::ATan2(py, pz);

		Double_t chi2 = ComputeChi2(track);

		// Create MUON track, fill it and add it to the ESD.
		AliESDMuonTrack mt;
		mt.SetInverseBendingMomentum(invmom);
		mt.SetThetaX(thetax);
		mt.SetThetaY(thetay);

		// Our algorithm assumes the particle originates from the origin.
		mt.SetZ(0.0);
		mt.SetBendingCoor(0.0);
		mt.SetNonBendingCoor(0.0);

		mt.SetChi2(chi2);
		mt.SetNHit(2);  // Always 2, thats the way Manso's algorithm works.
		mt.SetMatchTrigger( chi2 != -1. ? 1 : 0 );
		mt.SetChi2MatchTrigger(chi2);

		AliDebug(2, Form(
			"Adding AliESDMuonTrack with inverse bending momentum"
			" = %f, theta X = %f, theta Y = %f, chi2 = %f",
			invmom, thetax, thetay, chi2)
		);

		event->AddMuonTrack(&mt);
	}

#ifndef LOG_NO_DEBUG
	AliDebug(4, str);
#endif // LOG_NO_DEBUG

	return 0;
}

