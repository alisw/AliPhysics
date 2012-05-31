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

// Simple macro to generate N-tuple for performance analysis.
// This macro must be compiled and run like so from within the aliroot command
// prompt:
//   root [0] gSystem->Load("libAliHLTMUON.so");
//   root [0] .x MakeHitsTable.C+

#include "TVector3.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TParticle.h"
#include "TArrayF.h"
#include "TStopwatch.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TError.h"
#include "TMath.h"

#include "AliCDBManager.h"
#include "AliLoader.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONHit.h"
#include "AliHLTMUONEvent.h"
#include "AliHLTMUONRootifierComponent.h"
#include "AliHLTMUONMansoTrack.h"

#include <cstdlib>
#include <vector>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;


void MakeTrackTable(
		Int_t firstEvent = 0,
		Int_t lastEvent = -1,
		const char* dHLToutputfile = "output.root",
		Float_t maxSigma = 4., // 4 standard deviations
		Float_t sigmaX = 0.1,  // 1 mm resolution
		Float_t sigmaY = 0.01, // 100 micron resolution
		Float_t sigmaZ = 0.02,  // 200 microns resolution
		Float_t sigmaXtrg = 0.5,  // 5 mm resolution
		Float_t sigmaYtrg = 0.5,  // 5 mm resolution
		Float_t sigmaZtrg = 0.02  // 2 microns resolution
	)
{
	// Setup the CDB default storage and run number if nothing was set.
	AliCDBManager* cdbManager = AliCDBManager::Instance();
	if (cdbManager == NULL)
	{
		cerr << "ERROR: Global CDB manager object does not exist." << endl;
		return;
	}
	if (cdbManager->GetDefaultStorage() == NULL)
	{
		cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	}
	if (cdbManager->GetRun() == -1)
	{
		cdbManager->SetRun(0);
	}
	
	gSystem->Load("libAliHLTMUON.so");
	// Must pree load libAliHLTMUON.so before loading this macro and running it in compiled mode.

	TString fieldnames = "event:isprimary:hastrack:cantrigger:pdgcode:sign:px:py:pz";
	for (Int_t i = 0; i <= 14; i++)
	{
		fieldnames += ":x";
		fieldnames += i;
		fieldnames += ":y";
		fieldnames += i;
		fieldnames += ":z";
		fieldnames += i;
	}
	fieldnames += ":dHLTsign:dHLTpx:dHLTpy:dHLTpz:dHLTfit";
	for (Int_t i = 1; i <= 14; i++)
	{
		fieldnames += ":dHLTx";
		fieldnames += i;
		fieldnames += ":dHLTy";
		fieldnames += i;
		fieldnames += ":dHLTz";
		fieldnames += i;
	}
	TNtuple* table = new TNtuple("tracktable", "Table of tracks.", fieldnames);
	
	TFile dHLTfile(dHLToutputfile, "READ");

	AliMUONMCDataInterface data;
	AliStack* stack;
	
	Int_t nEvents = data.NumberOfEvents();
	if (lastEvent < 0) lastEvent = nEvents-1;
	for (Int_t event = firstEvent; event <= lastEvent; event++)
	{
		cout << "Processing event: " << event;
		
		char buf[1024];
		char* str = &buf[0];
		sprintf(str, "AliHLTMUONEvent;%d", event+1);
		AliHLTMUONEvent* dHLTevent = dynamic_cast<AliHLTMUONEvent*>( dHLTfile.Get(str) );
		if (dHLTevent == NULL)
		{
			cout << endl;
			cerr << "ERROR: Could not find " << str << " in " << dHLToutputfile << endl;
			continue;
		}
		
		data.GetEvent(event);
		
		stack = data.Stack(event);
		if (stack == NULL)
		{
			cout << endl;
			cerr << "ERROR: Stack is NULL" << endl;
			continue;
		}
		
		Int_t trackcount = data.NumberOfParticles();
		cout << "\ttrackcount = " << trackcount << endl;
		
		TArrayF isPrimary(trackcount);
		TArrayF hasTrack(trackcount);
		TArrayI hitsInTrigger(trackcount);
		TArrayF pdgcode(trackcount);
		TArrayF sign(trackcount);
		TArrayF px(trackcount);
		TArrayF py(trackcount);
		TArrayF pz(trackcount);
		TArrayF vx(trackcount);
		TArrayF vy(trackcount);
		TArrayF vz(trackcount);
		TArrayF dHLTsign(trackcount);
		TArrayF dHLTpx(trackcount);
		TArrayF dHLTpy(trackcount);
		TArrayF dHLTpz(trackcount);
		TArrayF dHLTfit(trackcount);
		
		//cout << "Fill particle data" << endl;
		for (Int_t i = 0; i < trackcount; i++)
		{
			TParticle* p = data.Particle(i);
			//isPrimary[i] = (i < stack->GetNprimary()) ? 1 : 0;
			isPrimary[i] = stack->IsPhysicalPrimary(i) ? 1 : 0;
			hasTrack[i] = 0;
			hitsInTrigger[i] = 0;
			pdgcode[i] = p->GetPdgCode();
			if (p->GetPDG() != NULL)
			{
				Double_t charge = p->GetPDG()->Charge();
				if (charge > 0)
					sign[i] = 1;
				else if (charge < 0)
					sign[i] = -1;
				else
					sign[i] = 0;
			}
			else
			{
				sign[i] = 0;
				cerr << "ERROR: Could not get the PDG information." << endl;
			}
			px[i] = p->Px();
			py[i] = p->Py();
			pz[i] = p->Pz();
			vx[i] = p->Vx();
			vy[i] = p->Vy();
			vz[i] = p->Vz();
			dHLTsign[i] = 0;
			dHLTpx[i] = 0;
			dHLTpy[i] = 0;
			dHLTpz[i] = 0;
			dHLTfit[i] = -1;
		}
		
		TArrayF hitX[14];
		TArrayF hitY[14];
		TArrayF hitZ[14];
		TArrayF dHLThitX[14];
		TArrayF dHLThitY[14];
		TArrayF dHLThitZ[14];
		for (Int_t i = 0; i < 14; i++)
		{
			hitX[i].Set(trackcount);
			hitY[i].Set(trackcount);
			hitZ[i].Set(trackcount);
			dHLThitX[i].Set(trackcount);
			dHLThitY[i].Set(trackcount);
			dHLThitZ[i].Set(trackcount);
			for (Int_t j = 0; j < trackcount; j++)
			{
				hitX[i][j] = 0;
				hitY[i][j] = 0;
				hitZ[i][j] = 0;
				dHLThitX[i][j] = 0;
				dHLThitY[i][j] = 0;
				dHLThitZ[i][j] = 0;
			}
		}

		//cout << "Fill hits" << endl;
		for (Int_t i = 0; i < data.NumberOfTracks(); i++)
		for (Int_t j = 0; j < data.NumberOfHits(i); j++)
		{
			AliMUONHit* hit = data.Hit(i, j);
			Int_t chamber = hit->Chamber() - 1;
			if (chamber < 0 || chamber >= 14)
			{
				cerr << "ERROR: Chamber number " << chamber << " is out of range." << endl;
				continue;
			}
			
			//cout << "hit->Track() = " << hit->Track() << endl;
			if (hit->Track() < 0 || hit->Track() >= trackcount)
			{
				cerr << "ERROR: Track number " << hit->Track() << " is out of range. [0.."
					<< trackcount << ")" << endl;
				continue;
			}
			Int_t particleindex = hit->Track();

			hitX[chamber][particleindex] = hit->X();
			hitY[chamber][particleindex] = hit->Y();
			hitZ[chamber][particleindex] = hit->Z();
					
			if (chamber == 0)
			{
				hasTrack[particleindex] = 1;
			}
			if (chamber >= 10)
			{
				hitsInTrigger[particleindex]++;
			}
		} // hits loop
		
		// Ok now go through the dHLT information for this event and correlate
		// the dHLT tracks to the kine data.
		// This is done as follows: take a dHLT track. Then for each kine track calculate
		// the distance between each hit of the kine track and corresponding dHLT track.
		// Then sum up the distances and divide by the number of hits matched.
		// This gives us a fit quality /Chi2 type parameter. The kine track that
		// has the smallest sum of differences is then chosen as the kine track
		// correlating to the current dHLT track being processed.
		for (Int_t i = 0; i < dHLTevent->Array().GetEntriesFast(); i++)
		{
			if (dHLTevent->Array().At(i)->IsA() != AliHLTMUONMansoTrack::Class())
				continue;
			AliHLTMUONMansoTrack* mtrack = static_cast<AliHLTMUONMansoTrack*>(dHLTevent->Array().At(i));
			
			// Find the best fitting kine track.
			Double_t bestFitQuality = 1e30;
			Int_t bestFitTrack = -1;
			for (Int_t kinetrack = 0; kinetrack < trackcount; kinetrack++)
			{
				Double_t sumOfDiffs = 0;
				Double_t hitsMatched = 0;
				for (Int_t j = 7; j <= 10; j++)
				{
					const AliHLTMUONRecHit* hit = mtrack->Hit(j);
					if (hit == NULL) continue;
					if (hitX[j-1][kinetrack] == 0 && hitY[j-1][kinetrack] == 0 && hitZ[j-1][kinetrack] == 0)
						continue;
					
					TVector3 hV(hitX[j-1][kinetrack], hitY[j-1][kinetrack], hitZ[j-1][kinetrack]);
					TVector3 diff = hV - hit->Coordinate();
					Double_t diffX = diff.X() / sigmaX;
					Double_t diffY = diff.Y() / sigmaY;
					Double_t diffZ = diff.Z() / sigmaZ;
					if (diffX > maxSigma) continue;
					if (diffY > maxSigma) continue;
					if (diffZ > maxSigma) continue;
					sumOfDiffs += diffX*diffX + diffY*diffY + diffZ*diffZ;
					hitsMatched++;
				}
				if (mtrack->TriggerRecord() != NULL)
				{
					for (Int_t j = 11; j <= 14; j++)
					{
						const TVector3& hit = mtrack->TriggerRecord()->Hit(j);
						if (hit == TVector3(0,0,0)) continue;
						if (hitX[j-1][kinetrack] == 0 && hitY[j-1][kinetrack] == 0 && hitZ[j-1][kinetrack] == 0)
							continue;
						
						TVector3 hV(hitX[j-1][kinetrack], hitY[j-1][kinetrack], hitZ[j-1][kinetrack]);
						TVector3 diff = hV - hit;
						Double_t diffX = diff.X() / sigmaXtrg;
						Double_t diffY = diff.Y() / sigmaYtrg;
						Double_t diffZ = diff.Z() / sigmaZtrg;
						if (diffX > maxSigma) continue;
						if (diffY > maxSigma) continue;
						if (diffZ > maxSigma) continue;
						sumOfDiffs += diffX*diffX + diffY*diffY + diffZ*diffZ;
						hitsMatched++;
					}
				}
				
				// Now check the fit quality.
				if (hitsMatched <= 0) continue;
				Double_t fitQuality = TMath::Sqrt(sumOfDiffs) / hitsMatched;
				if (fitQuality < bestFitQuality)
				{
					bestFitQuality = fitQuality;
					bestFitTrack = kinetrack;
				}
			}
			
			if (bestFitTrack < 0)
			{
				// No track was matched so we need to add a fake track record.
				Float_t args[12+14*3+5+14*3];
				args[0] = event;
				args[1] = -1;
				args[2] = 0;
				args[3] = 0;
				args[4] = 0;
				args[5] = 0;
				args[6] = 0;
				args[7] = 0;
				args[8] = 0;
				args[9] = 0;
				args[10] = 0;
				args[11] = 0;
				for (Int_t j = 0; j < 14; j++)
				{			
					args[12+j*3+0] = 0;
					args[12+j*3+1] = 0;
					args[12+j*3+2] = 0;
				}
				args[12+14*3+0] = mtrack->Sign();
				args[12+14*3+1] = mtrack->Px();
				args[12+14*3+2] = mtrack->Py();
				args[12+14*3+3] = mtrack->Pz();
				args[12+14*3+4] = -1;
				for (Int_t j = 0; j < 14; j++)
				{			
					args[12+14*3+5+j*3+0] = 0;
					args[12+14*3+5+j*3+1] = 0;
					args[12+14*3+5+j*3+2] = 0;
				}
				for (Int_t j = 6; j < 10; j++)
				{
					const AliHLTMUONRecHit* hit = mtrack->Hit(j+1);
					if (hit == NULL) continue;
					args[12+14*3+5+j*3+0] = hit->X();
					args[12+14*3+5+j*3+1] = hit->Y();
					args[12+14*3+5+j*3+2] = hit->Z();
				}
				if (mtrack->TriggerRecord() != NULL)
				{
					for (Int_t j = 10; j < 14; j++)
					{
						const TVector3& hit = mtrack->TriggerRecord()->Hit(j+1);
						args[12+14*3+5+j*3+0] = hit.X();
						args[12+14*3+5+j*3+1] = hit.Y();
						args[12+14*3+5+j*3+2] = hit.Z();
					}
				}
				table->Fill(args);
			}
			else
			{
				// Fill the details about the dHLT track to the best fitting track info.
				dHLTsign[bestFitTrack] = mtrack->Sign();
				dHLTpx[bestFitTrack] = mtrack->Px();
				dHLTpy[bestFitTrack] = mtrack->Py();
				dHLTpz[bestFitTrack] = mtrack->Pz();
				dHLTfit[bestFitTrack] = bestFitQuality;
				for (Int_t j = 7; j <= 10; j++)
				{
					const AliHLTMUONRecHit* hit = mtrack->Hit(j);
					if (hit == NULL) continue;
					dHLThitX[j-1][bestFitTrack] = hit->X();
					dHLThitY[j-1][bestFitTrack] = hit->Y();
					dHLThitZ[j-1][bestFitTrack] = hit->Z();
				}
				if (mtrack->TriggerRecord() != NULL)
				{
					for (Int_t j = 11; j <= 14; j++)
					{
						const TVector3& hit = mtrack->TriggerRecord()->Hit(j);
						if (hit == TVector3(0,0,0)) continue;
						dHLThitX[j-1][bestFitTrack] = hit.X();
						dHLThitY[j-1][bestFitTrack] = hit.Y();
						dHLThitZ[j-1][bestFitTrack] = hit.Z();
					}
				}
			}
		}
		
		//cout << "Fill table" << endl;
		for (Int_t i = 0; i < trackcount; i++)
		{
			Float_t args[12+14*3+5+14*3];
			args[0] = event;
			args[1] = isPrimary[i];
			args[2] = hasTrack[i];
			args[3] = (hitsInTrigger[i] > 2) ? 1 : 0;
			args[4] = pdgcode[i];
			args[5] = sign[i];
			args[6] = px[i];
			args[7] = py[i];
			args[8] = pz[i];
			args[9] = vx[i];
			args[10] = vy[i];
			args[11] = vz[i];
			for (Int_t j = 0; j < 14; j++)
			{			
				args[12+j*3+0] = hitX[j][i];
				args[12+j*3+1] = hitY[j][i];
				args[12+j*3+2] = hitZ[j][i];
			}
			args[12+14*3+0] = dHLTsign[i];
			args[12+14*3+1] = dHLTpx[i];
			args[12+14*3+2] = dHLTpy[i];
			args[12+14*3+3] = dHLTpz[i];
			args[12+14*3+4] = dHLTfit[i];
			for (Int_t j = 0; j < 14; j++)
			{			
				args[12+14*3+5+j*3+0] = dHLThitX[j][i];
				args[12+14*3+5+j*3+1] = dHLThitY[j][i];
				args[12+14*3+5+j*3+2] = dHLThitZ[j][i];
			}
			table->Fill(args);
		}
		
	} // event loop

	TFile file("tracktable.root", "RECREATE");
	file.cd();
	table->Write(table->GetName(), TObject::kOverwrite);

	cout << "Done." << endl;
}
