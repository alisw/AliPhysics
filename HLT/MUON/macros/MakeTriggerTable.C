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


void MakeTriggerTable(
		Int_t firstEvent = 0,
		Int_t lastEvent = -1,
		const char* L0outputfile = "output.root",
		Float_t maxSigma = 4., // 4 standard deviations
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

	TString fieldnames = "event:isprimary:pdgcode:sign:px:py:pz";
	for (Int_t i = 11; i <= 14; i++)
	{
		fieldnames += ":x";
		fieldnames += i;
		fieldnames += ":y";
		fieldnames += i;
		fieldnames += ":z";
		fieldnames += i;
	}
	fieldnames += ":L0sign:L0px:L0py:L0pz:L0fit";
	for (Int_t i = 11; i <= 14; i++)
	{
		fieldnames += ":L0x";
		fieldnames += i;
		fieldnames += ":L0y";
		fieldnames += i;
		fieldnames += ":L0z";
		fieldnames += i;
	}
	TNtuple* table = new TNtuple("triggertable", "Table of triggers.", fieldnames);
	
	TFile L0file(L0outputfile, "READ");

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
		AliHLTMUONEvent* L0event = dynamic_cast<AliHLTMUONEvent*>( L0file.Get(str) );
		if (L0event == NULL)
		{
			cout << endl;
			cerr << "ERROR: Could not find " << str << " in " << L0outputfile << endl;
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
		TArrayF pdgcode(trackcount);
		TArrayF sign(trackcount);
		TArrayF px(trackcount);
		TArrayF py(trackcount);
		TArrayF pz(trackcount);
		TArrayF L0sign(trackcount);
		TArrayF L0px(trackcount);
		TArrayF L0py(trackcount);
		TArrayF L0pz(trackcount);
		TArrayF L0fit(trackcount);
		
		//cout << "Fill particle data" << endl;
		for (Int_t i = 0; i < trackcount; i++)
		{
			TParticle* p = data.Particle(i);
			isPrimary[i] = stack->IsPhysicalPrimary(i) ? 1 : 0;
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
			L0sign[i] = 0;
			L0px[i] = 0;
			L0py[i] = 0;
			L0pz[i] = 0;
			L0fit[i] = -1;
		}
		
		TArrayF hitX[4];
		TArrayF hitY[4];
		TArrayF hitZ[4];
		TArrayF L0hitX[4];
		TArrayF L0hitY[4];
		TArrayF L0hitZ[4];
		for (Int_t i = 0; i < 4; i++)
		{
			hitX[i].Set(trackcount);
			hitY[i].Set(trackcount);
			hitZ[i].Set(trackcount);
			L0hitX[i].Set(trackcount);
			L0hitY[i].Set(trackcount);
			L0hitZ[i].Set(trackcount);
			for (Int_t j = 0; j < trackcount; j++)
			{
				hitX[i][j] = 0;
				hitY[i][j] = 0;
				hitZ[i][j] = 0;
				L0hitX[i][j] = 0;
				L0hitY[i][j] = 0;
				L0hitZ[i][j] = 0;
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

			if (chamber >= 10)
			{
				hitX[chamber-10][particleindex] = hit->X();
				hitY[chamber-10][particleindex] = hit->Y();
				hitZ[chamber-10][particleindex] = hit->Z();
			}
		} // hits loop
		
		// Ok now go through the L0 information for this event and correlate
		// the L0 trigger records to the kine data.
		// This is done as follows: take a trigger record. Then for each kine track
		// calculate the distance between each hit of the kine track and corresponding
		// L0 trigger track fragment hit.
		// Then sum up the distances and divide by the number of hits matched.
		// This gives us a fit quality /Chi2 type parameter. The kine track that
		// has the smallest sum of differences is then chosen as the kine track
		// correlating to the current L0 trigger record being processed.
		for (Int_t i = 0; i < L0event->Array().GetEntriesFast(); i++)
		{
			if (L0event->Array().At(i)->IsA() != AliHLTMUONTriggerRecord::Class())
				continue;
			AliHLTMUONTriggerRecord* trigrec = static_cast<AliHLTMUONTriggerRecord*>(L0event->Array().At(i));
			
			// Find the best fitting kine track.
			Double_t bestFitQuality = 1e30;
			Int_t bestFitTrack = -1;
			for (Int_t kinetrack = 0; kinetrack < trackcount; kinetrack++)
			{
				Double_t sumOfDiffs = 0;
				Double_t hitsMatched = 0;
				for (Int_t j = 11; j <= 14; j++)
				{
					const TVector3& hit = trigrec->Hit(j);
					if (hit == TVector3(0,0,0)) continue;
					TVector3 hV(hitX[j-11][kinetrack], hitY[j-11][kinetrack], hitZ[j-11][kinetrack]);
					if (hV == TVector3(0,0,0)) continue;
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
				
				// Now check the fit quality.
				if (hitsMatched <= 0) continue;
				Double_t fitQuality = sumOfDiffs / hitsMatched;
				if (fitQuality < bestFitQuality)
				{
					bestFitQuality = fitQuality;
					bestFitTrack = kinetrack;
				}
			}
			
			if (bestFitTrack < 0)
			{
				// No track was matched so we need to add a fake track record.
				Float_t args[7+4*3+5+4*3];
				args[0] = event;
				args[1] = -1;
				args[2] = 0;
				args[3] = 0;
				args[4] = 0;
				args[5] = 0;
				args[6] = 0;
				for (Int_t j = 0; j < 4; j++)
				{			
					args[7+j*3+0] = 0;
					args[7+j*3+1] = 0;
					args[7+j*3+2] = 0;
				}
				args[7+4*3+0] = trigrec->Sign();
				args[7+4*3+1] = trigrec->Px();
				args[7+4*3+2] = trigrec->Py();
				args[7+4*3+3] = trigrec->Pz();
				args[7+4*3+4] = -1;
				for (Int_t j = 0; j < 4; j++)
				{
					const TVector3& hit = trigrec->Hit(j+11);
					args[7+4*3+5+j*3+0] = hit.X();
					args[7+4*3+5+j*3+1] = hit.Y();
					args[7+4*3+5+j*3+2] = hit.Z();
				}
				table->Fill(args);
			}
			else
			{
				// Fill the details about the L0 track to the best fitting track info.
				L0sign[bestFitTrack] = trigrec->Sign();
				L0px[bestFitTrack] = trigrec->Px();
				L0py[bestFitTrack] = trigrec->Py();
				L0pz[bestFitTrack] = trigrec->Pz();
				L0fit[bestFitTrack] = bestFitQuality;
				for (Int_t j = 11; j <= 14; j++)
				{
					const TVector3& hit = trigrec->Hit(j);
					if (hit == TVector3(0,0,0)) continue;
					L0hitX[j-11][bestFitTrack] = hit.X();
					L0hitY[j-11][bestFitTrack] = hit.Y();
					L0hitZ[j-11][bestFitTrack] = hit.Z();
				}
			}
		}
		
		//cout << "Fill table" << endl;
		for (Int_t i = 0; i < trackcount; i++)
		{
			Float_t args[7+4*3+5+4*3];
			args[0] = event;
			args[1] = isPrimary[i];
			args[2] = pdgcode[i];
			args[3] = sign[i];
			args[4] = px[i];
			args[5] = py[i];
			args[6] = pz[i];
			for (Int_t j = 0; j < 4; j++)
			{			
				args[7+j*3+0] = hitX[j][i];
				args[7+j*3+1] = hitY[j][i];
				args[7+j*3+2] = hitZ[j][i];
			}
			args[7+4*3+0] = L0sign[i];
			args[7+4*3+1] = L0px[i];
			args[7+4*3+2] = L0py[i];
			args[7+4*3+3] = L0pz[i];
			args[7+4*3+4] = L0fit[i];
			for (Int_t j = 0; j < 4; j++)
			{			
				args[7+4*3+5+j*3+0] = L0hitX[j][i];
				args[7+4*3+5+j*3+1] = L0hitY[j][i];
				args[7+4*3+5+j*3+2] = L0hitZ[j][i];
			}
			table->Fill(args);
		}
		
	} // event loop

	TFile file("triggertable.root", "RECREATE");
	file.cd();
	table->Write(table->GetName(), TObject::kOverwrite);

	cout << "Done." << endl;
}
