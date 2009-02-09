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
#include "AliHLTMUONRecHit.h"
#include "AliHLTMUONTriggerRecord.h"

#include <cstdlib>
#include <vector>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;


void MakeHitsTable(
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

	TString fieldnames = "event:isprimary:pdgcode:detelem:chamber:x:y:z:dHLTx:dHLTy:dHLTz";
	TNtuple* table = new TNtuple("hittable", "Table of hits.", fieldnames);
	
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
		TArrayF pdgcode(trackcount);
		
		//cout << "Fill particle data" << endl;
		for (Int_t i = 0; i < trackcount; i++)
		{
			TParticle* p = data.Particle(i);
			isPrimary[i] = stack->IsPhysicalPrimary(i) ? 1 : 0;
			pdgcode[i] = p->GetPdgCode();
		}
		
		TArrayF detElem[14];
		TArrayF chamberNum[14];
		TArrayF hitX[14];
		TArrayF hitY[14];
		TArrayF hitZ[14];
		TArrayF dHLThitX[14];
		TArrayF dHLThitY[14];
		TArrayF dHLThitZ[14];
		for (Int_t i = 0; i < 14; i++)
		{
			detElem[i].Set(trackcount);
			chamberNum[i].Set(trackcount);
			hitX[i].Set(trackcount);
			hitY[i].Set(trackcount);
			hitZ[i].Set(trackcount);
			dHLThitX[i].Set(trackcount);
			dHLThitY[i].Set(trackcount);
			dHLThitZ[i].Set(trackcount);
			for (Int_t j = 0; j < trackcount; j++)
			{
				detElem[i][j] = 0;
				chamberNum[i][j] = 0;
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

			detElem[chamber][particleindex] = hit->DetElemId();
			chamberNum[chamber][particleindex] = hit->Chamber();
			hitX[chamber][particleindex] = hit->X();
			hitY[chamber][particleindex] = hit->Y();
			hitZ[chamber][particleindex] = hit->Z();
		} // hits loop
		
		for (Int_t i = 0; i < dHLTevent->Array().GetEntriesFast(); i++)
		{
			if (dHLTevent->Array().At(i)->IsA() != AliHLTMUONRecHit::Class())
				continue;
			AliHLTMUONRecHit* hit = static_cast<AliHLTMUONRecHit*>(dHLTevent->Array().At(i));
			
			// Find the best fitting hit.
			Double_t bestFitQuality = 1e30;
			Int_t bestFitTrack = -1;
			Int_t bestFitCh = -1;
			for (Int_t ti = 0; ti < trackcount; ti++)
			for (Int_t ch = 0; ch < 14; ch++)
			{
				if (hit->Coordinate() == TVector3(0,0,0)) continue;
				if (hitX[ch][ti] == 0 && hitY[ch][ti] == 0 && hitZ[ch][ti] == 0)
					continue;
				TVector3 hV(hitX[ch][ti], hitY[ch][ti], hitZ[ch][ti]);
				TVector3 diff = hV - hit->Coordinate();
				Double_t diffX = diff.X() / sigmaX;
				Double_t diffY = diff.Y() / sigmaY;
				Double_t diffZ = diff.Z() / sigmaZ;
				if (diffX > maxSigma) continue;
				if (diffY > maxSigma) continue;
				if (diffZ > maxSigma) continue;
				Double_t fitQuality = diffX*diffX + diffY*diffY + diffZ*diffZ;

				// Now check the fit quality.
				if (fitQuality < bestFitQuality)
				{
					bestFitQuality = fitQuality;
					bestFitTrack = ti;
					bestFitCh = ch;
				}
			}
			
			if (bestFitTrack < 0 || bestFitCh < 0)
			{
				// No hit was matched so we need to add a fake hit.
				table->Fill(event, -1, 0, -1, -1, 0, 0, 0, hit->X(), hit->Y(), hit->Z());
			}
			else
			{
				// Fill the details about the dHLT track to the best fitting track info.
				dHLThitX[bestFitCh][bestFitTrack] = hit->X();
				dHLThitY[bestFitCh][bestFitTrack] = hit->Y();
				dHLThitZ[bestFitCh][bestFitTrack] = hit->Z();
			}
		}
		
		for (Int_t i = 0; i < dHLTevent->Array().GetEntriesFast(); i++)
		{
			if (dHLTevent->Array().At(i)->IsA() != AliHLTMUONTriggerRecord::Class())
				continue;
			AliHLTMUONTriggerRecord* trigrec = static_cast<AliHLTMUONTriggerRecord*>(dHLTevent->Array().At(i));
			
			for (Int_t n = 11; n <= 14; n++)
			{
				const TVector3& hit = trigrec->Hit(n);
				if (hit == TVector3(0,0,0)) continue;
				Int_t ch = n-1;
				
				// Find the best fitting hit.
				Double_t bestFitQuality = 1e30;
				Int_t bestFitTrack = -1;
				Int_t bestFitCh = -1;
				for (Int_t ti = 0; ti < trackcount; ti++)
				{
					if (hit == TVector3(0,0,0)) continue;
					if (hitX[ch][ti] == 0 && hitY[ch][ti] == 0 && hitZ[ch][ti] == 0)
						continue;
					TVector3 hV(hitX[ch][ti], hitY[ch][ti], hitZ[ch][ti]);
					TVector3 diff = hV - hit;
					Double_t diffX = diff.X() / sigmaXtrg;
					Double_t diffY = diff.Y() / sigmaYtrg;
					Double_t diffZ = diff.Z() / sigmaZtrg;
					if (diffX > maxSigma) continue;
					if (diffY > maxSigma) continue;
					if (diffZ > maxSigma) continue;
					Double_t fitQuality = diffX*diffX + diffY*diffY + diffZ*diffZ;

					// Now check the fit quality.
					if (fitQuality < bestFitQuality)
					{
						bestFitQuality = fitQuality;
						bestFitTrack = ti;
						bestFitCh = ch;
					}
				}
				
				if (bestFitTrack < 0 || bestFitCh < 0)
				{
					// No hit was matched so we need to add a fake hit.
					table->Fill(event, -1, 0, -1, -1, 0, 0, 0, hit.X(), hit.Y(), hit.Z());
				}
				else
				{
					// Fill the details about the dHLT track to the best fitting track info.
					dHLThitX[bestFitCh][bestFitTrack] = hit.X();
					dHLThitY[bestFitCh][bestFitTrack] = hit.Y();
					dHLThitZ[bestFitCh][bestFitTrack] = hit.Z();
				}
			}
		}
		
		//cout << "Fill table" << endl;
		for (Int_t i = 0; i < trackcount; i++)
		for (Int_t j = 0; j < 14; j++)
		{
			//event:isprimary::pdgcode:detelem:chamber:x:y:z:dHLTx:dHLTy:dHLTz
			table->Fill(
					event,
					isPrimary[i],
					pdgcode[i],
					detElem[j][i],
					chamberNum[j][i],
					hitX[j][i],
					hitY[j][i],
					hitZ[j][i],
					dHLThitX[j][i],
					dHLThitY[j][i],
					dHLThitZ[j][i]
				);
		}
		
	} // event loop

	TFile file("hittable.root", "RECREATE");
	file.cd();
	table->Write(table->GetName(), TObject::kOverwrite);

	cout << "Done." << endl;
}
