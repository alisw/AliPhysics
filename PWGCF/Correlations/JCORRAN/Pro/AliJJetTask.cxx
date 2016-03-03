/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
//
// Jet fragmentation transverse momentum (j_T) analysis task
//
// Author: Beomkyu Kim, Beomsu Chang, Dongjo Kim

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TSystem.h>
#include <TFile.h>

#include "AliCentrality.h"



#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"
#include "AliJBaseTrack.h"
#include "AliJJet.h"
#include "AliJJetTask.h"


ClassImp(AliJJetTask);

//________________________________________________________________________
AliJJetTask::AliJJetTask() : 
	AliAnalysisTaskEmcalJet("AliJJetTask", kTRUE),
	fJetsCont(),
	fTracksCont(),
	fCaloClustersCont(),
	fJTracks("AliJBaseTrack",1000),
	fJJets(),
	fTaskEntry(-1),
	fJetFinderString(),
	fTrackArrayName("nonejk"),
	fNJetFinder(0),
	debug(0)

{
	// Default constructor.


	SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJJetTask::AliJJetTask(const char *name, const int nJetFinder) : 
	AliAnalysisTaskEmcalJet(name, kTRUE),
	fJetsCont(nJetFinder),
	fTracksCont(nJetFinder),
	fCaloClustersCont(nJetFinder),
	fJTracks("AliJBaseTrack",1000),
	fJJets(),
	fTaskEntry(-1),
	fJetFinderString(nJetFinder),
	fTrackArrayName("nonejk"),
	fNJetFinder(nJetFinder),
	debug(0)
{
	SetMakeGeneralHistograms(kTRUE);
}

AliJJetTask::AliJJetTask(const AliJJetTask& ap) :
	AliAnalysisTaskEmcalJet(ap.fName, kTRUE),
	fJetsCont(ap.fJetsCont),
	fTracksCont(ap.fTracksCont),
	fCaloClustersCont(ap.fCaloClustersCont),
	fJTracks(ap.fJTracks),
	fJJets(ap.fJJets),
	fTaskEntry(ap.fTaskEntry),
	fJetFinderString(ap.fJetFinderString),
	fTrackArrayName(ap.fTrackArrayName),
	fNJetFinder(ap.fNJetFinder),
	debug(ap.debug)
{

}


AliJJetTask& AliJJetTask::operator = (const AliJJetTask& ap)
{

	this->~AliJJetTask();
	new(this) AliJJetTask(ap);
	return *this;
}

//________________________________________________________________________
AliJJetTask::~AliJJetTask()
{
	// Destructor.

	//delete[] fJJets;
}




//________________________________________________________________________
void AliJJetTask::UserCreateOutputObjects()
{
	// Create user output.

	fJetsCont.resize(fNJetFinder); 
	fJJets.clear();
	fJJets.resize(fNJetFinder, TClonesArray("AliJJet",1000));
	fTracksCont.resize(fNJetFinder); 
	fCaloClustersCont.resize(fNJetFinder); 



	AliAnalysisTaskEmcalJet::UserCreateOutputObjects();


	//fJJets = new TClonesArray[fNJetFinder];


	fJetFinderString.clear();
	for (int i=0; i<fNJetFinder; i++){
		fJetsCont[i]           = GetJetContainer(i);
		fJetFinderString.push_back(fJetsCont[i]->GetArrayName());
		cout << i <<"\t" << fJetFinderString[i] << endl;
		fTracksCont[i]       = GetParticleContainer(0);
		fCaloClustersCont[i] = GetClusterContainer(0);
		fTracksCont[i]->SetClassName("AliVTrack");
		fCaloClustersCont[i]->SetClassName("AliAODCaloCluster");
		//fJJets.push_back(TClonesArray("AliJJet",1000));
	}

}

//________________________________________________________________________
Bool_t AliJJetTask::FillHistograms()
{



    for (int itrack = 0; itrack<fTracksCont[0]->GetNParticles(); itrack++){
        AliVTrack *track = static_cast<AliVTrack*>(fTracksCont[0]->GetParticle(itrack));
        new (fJTracks[itrack]) AliJBaseTrack(track->Px(),track->Py(), track->Pz(), track->E(), itrack,0,0);

    }     




	for (int i=0; i<fNJetFinder; i++){
		fJetsCont[i]->ResetCurrentID();
		AliEmcalJet *jet = fJetsCont[i]->GetNextAcceptJet();
		int iJet =0; 

		//fills fJJets[icontainer][ijet] and histograms        
		while(jet){
			TClonesArray & jets = fJJets[i];
			new (jets[iJet]) AliJJet(jet->Px(),jet->Py(), jet->Pz(), jet->E(), jet->GetLabel(),0,0);
			AliJJet * j = (AliJJet*) fJJets[i][iJet];   
			j->SetArea(jet->Area() );


			Int_t nConstituents = jet->GetNumberOfTracks();

			for (int it=0; it<nConstituents; it++){
				int iTrack = jet->TrackAt(it);
				j->AddConstituent(fJTracks[iTrack]);

			}

			if (debug>0) { cout<<"    iContainer    : "<<i<< 
				"    iJet          : "<<iJet<<
					"    nConstituents : "<<nConstituents<<endl;
			}

			//Goes to the next jet
			jet = fJetsCont[i]->GetNextAcceptJet();
			if (debug>0) {
				cout<<"  fJJets N lists : "<<fJJets[i].GetEntries()<<
					"  fJJets constituents : "<<((AliJJet*)fJJets[i][iJet])->GetConstituents()->GetEntries()<<endl;
			}

			iJet++;
		}


	}


	return kTRUE;

}




//________________________________________________________________________
void AliJJetTask::ExecOnce() {

	if(debug > 0){
		cout << "AliJJetTask::ExecOnce(): " << endl;
	}
	AliAnalysisTaskEmcalJet::ExecOnce();

	for (int i=0; i<fNJetFinder; i++){
		if (fJetsCont[i] && fJetsCont[i]->GetArray() == 0) fJetsCont[i] = 0;
		if (fTracksCont[i] && fTracksCont[i]->GetArray() == 0) fTracksCont[i] = 0;
		if (fCaloClustersCont[i] && fCaloClustersCont[i]->GetArray() == 0) fCaloClustersCont[i] = 0;
	}
}

//________________________________________________________________________
Bool_t AliJJetTask::Run()
{
	// Run analysis code here, if needed. It will be executed before FillHistograms().
	fTaskEntry = fEntry;
	for (int i=0; i<fNJetFinder; i++){
		fJJets[i]. Clear();
	}
	fJTracks.Clear(); 

	if (debug >0 && Entry()%1000 ==0 ) cout<<Entry()<<endl;
	return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliJJetTask::Terminate(Option_t *) 
{
	// Called once at the end of the analysis.
}

