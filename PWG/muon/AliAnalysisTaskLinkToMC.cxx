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

/* $Id$ */

///
/// @file   AliAnalysisTaskLinkToMC.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   27 Oct 2008
/// @brief  Implementation of the AliAnalysisTaskLinkToMC task.
///
/// The AliAnalysisTaskLinkToMC connects reconstructed tracks found in the ESD
/// to their corresponding Monte Carlo (MC) tracks via the proximity of their
/// hit coordinates. Each ESD track is copied over to an AOD track and put into
/// the AOD event container, but the AOD track label is filled such that the
/// AOD track label is the same as the MC track label.
/// If a MC track could not be found for a given ESD track then -1 is filled
/// into the AOD track label.
/// The corresponding MC track is found by correlating the reconstructed track
/// clusters with the MC track references. For each cluster the closest track
/// reference coordinate is found and a chi-squared value is computed for each
/// X and Y dimension independently. The chi-square's are summed up and the MC
/// track with the smallest total chi-squared is chosen as the best match to
/// the ESD track.
///
/// The following cuts control the conditions under which we consider a MC track
/// to match a given ESD track:
/// 1) SigmaCut(Double_t value) - This is maximum distance expressed in standard
///      deviations that is allowed between a reconstructed cluster coordinate
///      and the MC track reference coordinate. The condition is expressed as:
///        (x_{cluster} - x_{trackref})^2 / (sigma_{x})^2 < #sigma_{cut}
///      This conditions is applied to each X and Y dimension respectively and
///      cluster / track reference pairs that do not pass this cut are considered
///      not to match.
/// 2) MinClusters(Int_t value) - Indicates the minimum number of clusters that
///      must be matched to a track reference from a MC track for the ESD and
///      MC track to be considered to match.
/// The following are absolute cuts which are used mainly as a speed optimisation
/// so as to skip over cluster / track reference pairs that are far away from
/// each other early.
/// 3) HardCutLimitX(Double_t value) - The absolute maximum distance in centimetres
///      allowed between a cluster's X coordinate and a track reference's X coordinate.
///      If a cluster / track reference pair does not pass this condition then the
///      pair is considered not to match.
/// 4) HardCutLimitY(Double_t value) - Similar to HardCutLimitX but applied to the
///      Y dimension of the cluster and track reference coordinate.
/// 5) HardCutLimitZ(Double_t value) - Similar to HardCutLimitX but applied to the
///      Z dimension.
///
/// The analysis task also optionally generates histograms which show can be used
/// to check the quality of the correlation between ESD and MC tracks.
/// The histograms are generated in 2D for the pT and rapidity variable,
/// and are called:
/// 1) findableTracksHist - Filled with the MC tracks which should in principle
///      be findable by the reconstruction software. A findable track is defined
///      by the IsFindable method.
/// 2) foundTracksHistMC - Filled with the MC information for each muon track that
///      was found in the ESD and for which we could match its corresponding MC track.
/// 3) foundTracksHist - Filled just like the foundTracksHistMC histogram, but with
///      the momentum information taken from the ESD track and not the corresponding
///      MC track.
/// 4) fakeTracksHist - Filled with all the ESD tracks for which we could not find
///      its matching MC track and assume it is a fake track.
/// These histograms are filled by default, but can be disabled by using the
/// GenerateHistograms(false) method call.
///

#include "AliAnalysisTaskLinkToMC.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODVertex.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliTrackReference.h"
#include "AliStack.h"
#include "AliLog.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TMap.h"
#include "TList.h"
#include "TH2D.h"
#include "TStyle.h"

ClassImp(AliAnalysisTaskLinkToMC)


void AliAnalysisTaskLinkToMC::HardCutLimitX(Double_t value)
{
	/// Sets the absolute maximum distance in centimetres between the cluster
	/// and track reference that is allowed in the X direction.
	
	if (value >= 0)
	{
		fHardCutLimitX = value;
	}
	else
	{
		AliWarning("Setting value for fHardCutLimitX which is negative."
			" Will set it to a positive value."
			);
		fHardCutLimitX = TMath::Abs(value);
	}
}


void AliAnalysisTaskLinkToMC::HardCutLimitY(Double_t value)
{
	/// Sets the absolute maximum distance in centimetres between the cluster
	/// and track reference that is allowed in the Y direction.
	
	if (value >= 0)
	{
		fHardCutLimitY = value;
	}
	else
	{
		AliWarning("Setting value for fHardCutLimitY which is negative."
			" Will set it to a positive value."
			);
		fHardCutLimitY = TMath::Abs(value);
	}
}


void AliAnalysisTaskLinkToMC::HardCutLimitZ(Double_t value)
{
	/// Sets the absolute maximum distance in centimetres between the cluster
	/// and track reference that is allowed in the Z direction.
	
	if (value >= 0)
	{
		fHardCutLimitZ = value;
	}
	else
	{
		AliWarning("Setting value for fHardCutLimitZ which is negative."
			" Will set it to a positive value."
			);
		fHardCutLimitZ = TMath::Abs(value);
	}
}


Double_t AliAnalysisTaskLinkToMC::SigmaCut() const
{
	/// Returns the maximum sigma value to allow for each X and Y dimension
	/// between clusters and track references.
	
	return TMath::Sqrt(fSigmaCut2);
}


Bool_t AliAnalysisTaskLinkToMC::ChamberMustMatch(Int_t chamber) const
{
	/// Returns the flag indicating if the ESD track cluster on the given
	/// chamber must match a corresponding track reference to be valid.
	/// \param chamber The chamber number must be in the range 1 to 14.
	
	if (chamber >= 1 and chamber <= 14)
	{
		return fChamberMustMatch[chamber-1];
	}
	else
	{
		AliError(Form("The chamber number %d is invalid. Expected a value"
			" in the range [1..14]. Will return a value of 'false'.",
			chamber
			));
		return false;
	}
}


void AliAnalysisTaskLinkToMC::ChamberMustMatch(Int_t chamber, Bool_t value)
{
	/// Set the flag indicating if the ESD track cluster on the given chamber
	/// must match a corresponding track reference to be valid.
	/// \param chamber The chamber number must be in the range 1 to 14.
	/// \param value  The new value to set to.
	
	if (chamber >= 1 and chamber <= 14)
	{
		fChamberMustMatch[chamber-1] = value;
	}
	else
	{
		AliError(Form("The chamber number %d is invalid. Expected a value"
			" in the range [1..14]. Will not set any fChamberMustMatch flag.",
			chamber
			));
	}
}


Bool_t AliAnalysisTaskLinkToMC::StationMustMatch(Int_t station) const
{
	/// Returns the flag indicating if at least one ESD track cluster on the
	/// given station must match a corresponding track reference to be valid.
	/// \param station The station number must be in the range 1..7. Station
	///     1 to 5 are the tracking stations while 6 is MT1 and 7 is MT2.
	
	if (station >= 1 and station <= 7)
	{
		return fStationMustMatch[station-1];
	}
	else
	{
		AliError(Form("The station number %d is invalid. Expected a value"
			" in the range [1..7]. Will return a value of 'false'.",
			station
			));
		return false;
	}
}


void AliAnalysisTaskLinkToMC::StationMustMatch(Int_t station, Bool_t value)
{
	/// Sets the flag indicating if at least one ESD track cluster on the
	/// given station must match a corresponding track reference to be valid.
	/// \param station The station number must be in the range 1..7. Station
	///     1 to 5 are the tracking stations while 6 is MT1 and 7 is MT2.
	/// \param value  The new value to set to.
	
	if (station >= 1 and station <= 7)
	{
		fStationMustMatch[station-1] = value;
	}
	else
	{
		AliError(Form("The station number %d is invalid. Expected a value"
			" in the range [1..7]. Will not set any fStationMustMatch flag.",
			station
			));
	}
}


AliAnalysisTaskLinkToMC::AliAnalysisTaskLinkToMC() :
	AliAnalysisTaskSE(),
	fHistos(NULL),
	fFindableHist(NULL),
	fFoundHistMC(NULL),
	fFoundHist(NULL),
	fFakeHist(NULL),
	fHardCutLimitX(4.), // cm
	fHardCutLimitY(4.), // cm
	fHardCutLimitZ(4.), // cm
	fVarX(0.144*0.144), // cm^2
	fVarY(0.01*0.01), // cm^2
	fSigmaCut2(5*5),  // sigma^2
	fMinClusters(6),
	fMinClustersInSt12(0),  // no minimum requirement for stations 1 and 2.
	fMinClustersInSt45(3),  // 3 hits required on station 4 and 5 by default.
	fMakeHists(true),  // Generate histograms for monitoring by default
	fUseZCoordinate(false)  // Just use detElemId to match in the Z direction.
{
	/// Default constructor.
	/// Sets input slot 0 to TChain input for ESDs, output slot 0 to TTree
	/// for AODs and output slot 1 to TList for the list of output histograms.
	/// The output histograms generated are to cross check the quality of the
	/// correlation and their generation can be turned off with
	/// AliAnalysisTaskLinkToMC::GenerateHistograms(false).
	
	// Do not force any matching for particular chambers by default.
	for (Int_t i = 0; i < 14; i++)
	{
		fChamberMustMatch[i] = false;
	}
	// Only force a minimum of 1 hit to be matched on stations 1, 2 and 3 by default.
	for (Int_t i = 0; i < 3; i++)
	{
		fStationMustMatch[i] = true;
	}
	for (Int_t i = 3; i < 7; i++)
	{
		fStationMustMatch[i] = false;
	}
	
	DefineOutput(1, TList::Class());
}


AliAnalysisTaskLinkToMC::AliAnalysisTaskLinkToMC(const char* name) :
	AliAnalysisTaskSE(name),
	fHistos(NULL),
	fFindableHist(NULL),
	fFoundHistMC(NULL),
	fFoundHist(NULL),
	fFakeHist(NULL),
	fHardCutLimitX(4.), // cm
	fHardCutLimitY(4.), // cm
	fHardCutLimitZ(4.), // cm
	fVarX(0.144*0.144), // cm^2
	fVarY(0.01*0.01), // cm^2
	fSigmaCut2(5*5),  // sigma^2
	fMinClusters(6),
	fMinClustersInSt12(0),  // no minimum requirement for stations 1 and 2.
	fMinClustersInSt45(3),  // 3 hits required on station 4 and 5 by default.
	fMakeHists(true),  // Generate histograms for monitoring by default
	fUseZCoordinate(false)  // Just use detElemId to match in the Z direction.
{
	/// Constructor with a task name.
	/// Sets input slot 0 to TChain input for ESDs, output slot 0 to TTree
	/// for AODs and output slot 1 to TList for the list of output histograms
	/// storing histograms generated to cross check the quility of the
	/// correlation.
	/// \param name  The name of the task.
	
	// Do not force any matching for particular chambers by default.
	for (Int_t i = 0; i < 14; i++)
	{
		fChamberMustMatch[i] = false;
	}
	// Only force a minimum of 1 hit to be matched on stations 1, 2 and 3 by default.
	for (Int_t i = 0; i < 3; i++)
	{
		fStationMustMatch[i] = true;
	}
	for (Int_t i = 3; i < 7; i++)
	{
		fStationMustMatch[i] = false;
	}
	
	DefineOutput(1, TList::Class());
}


void AliAnalysisTaskLinkToMC::UserCreateOutputObjects()
{
	/// Creates the output histograms containing findable tracks, found tracks
	/// and fake tracks. The binning range for all histograms is 30 bins along
	/// the pT dimension (X direction) from 0 to 20 GeV/c and 30 bins in the
	/// rapidity dimension (Y direction) form -4 to -2.4 units.
	
	fHistos = new TList;
	
	char name[1024];
	char title[1024];
	
	Int_t nBinsX = 30;
	Double_t minX = 0.;
	Double_t maxX = 20.;
	Int_t nBinsY = 30;
	Double_t minY = -4.;
	Double_t maxY = -2.4;
	
	snprintf(name, 1024, "findableTracksHist");
	snprintf(title, 1024, "Findable tracks");
	fFindableHist = new TH2D(name, title, nBinsX, minX, maxX, nBinsY, minY, maxY);
	fFindableHist->SetXTitle("p_{T} [GeV/c]");
	fFindableHist->SetYTitle("Rapidity (Y)");
	fHistos->Add(fFindableHist);
	
	snprintf(name, 1024, "foundTracksHistMC");
	snprintf(title, 1024, "Found tracks (filled with Monte Carlo kinematic information)");
	fFoundHistMC = new TH2D(name, title, nBinsX, minX, maxX, nBinsY, minY, maxY);
	fFoundHistMC->SetXTitle("p_{T} [GeV/c]");
	fFoundHistMC->SetYTitle("Rapidity (Y)");
	fHistos->Add(fFoundHistMC);
	
	snprintf(name, 1024, "foundTracksHist");
	snprintf(title, 1024, "Found tracks (filled with reconstructed kinematic information)");
	fFoundHist = new TH2D(name, title, nBinsX, minX, maxX, nBinsY, minY, maxY);
	fFoundHist->SetXTitle("p_{T} [GeV/c]");
	fFoundHist->SetYTitle("Rapidity (Y)");
	fHistos->Add(fFoundHist);
	
	snprintf(name, 1024, "fakeTracksHist");
	snprintf(title, 1024, "Fake tracks");
	fFakeHist = new TH2D(name, title, nBinsX, minX, maxX, nBinsY, minY, maxY);
	fFakeHist->SetXTitle("p_{T} [GeV/c]");
	fFakeHist->SetYTitle("Rapidity (Y)");
	fHistos->Add(fFakeHist);
}


void AliAnalysisTaskLinkToMC::UserExec(Option_t* /*option*/)
{
	/// This is the top level method that performs the work of the task.
	/// It will fill the histograms for quality checking with the Monte
	/// Carlo (MC) information if they are to be generated.
	/// The input event is then processed which must be of type AliESDEvent.
	/// For all the ESD muon tracks we try to find the corresponding MC
	/// track and build a TMap for the mapping.
	/// If the AOD output event is available then it is filled with AOD
	/// tracks which copy the information form the ESD track and the AOD
	/// track label is filled with the corresponding MC track label.
	/// If a corresponding MC track could not be found then the AOD track
	/// label is filled with -1.
	/// Finally the checking histograms are further filled with the ESD
	/// information based on the ESD to MC track map.
	
	if (fMakeHists and
	    (fFindableHist == NULL or fFoundHistMC == NULL or fFoundHist == NULL or fFakeHist == NULL)
	   )
	{
		AliError("Histograms have not been created.");
		return;
	}
	
	if (fMakeHists and MCEvent() != NULL) FillHistsFromMC();
	
	if (InputEvent() != NULL)
	{
		if (InputEvent()->IsA() == AliESDEvent::Class())
		{
			if (MCEvent() != NULL)
			{
				TMap links;
				FindLinksToMC(links);
				if (AODEvent() != NULL) CreateAODTracks(links);
				if (fMakeHists) FillHistsFromLinks(links);
			}
			else
			{
				AliError("To process the ESD event we must have the Monte Carlo data also.");
			}
		}
		else
		{
			AliError(Form("The input event is of type \"%s\". Do not know how to handle it so it will be skipped.",
					InputEvent()->ClassName()
				));
		}
	}
	
	if (fMakeHists) PostData(1, fHistos);
}


void AliAnalysisTaskLinkToMC::Terminate(Option_t* /*option*/)
{
	/// The terminate method will draw the histograms filled by this task if
	/// they were filled by setting appropriate flag with the
	/// GenerateHistograms(true) method, which is the default.
	
	if (fMakeHists)
	{
		gStyle->SetPalette(1);
		new TCanvas;
		fFindableHist->DrawCopy("colz");
		new TCanvas;
		fFoundHistMC->DrawCopy("colz");
		new TCanvas;
		fFoundHist->DrawCopy("colz");
		new TCanvas;
		fFakeHist->DrawCopy("colz");
		
		AliInfo(Form("Findable tracks = %d", Int_t(fFindableHist->GetEntries()) ));
		AliInfo(Form("Found tracks    = %d", Int_t(fFoundHistMC->GetEntries()) ));
		AliInfo(Form("Fake tracks     = %d", Int_t(fFakeHist->GetEntries()) ));
	}
}


bool AliAnalysisTaskLinkToMC::IsFindable(AliMCParticle* track) const
{
	/// This method is used to identify if a Monte Carlo (MC) track is in principle
	/// findable in the muon spectrometer.
	/// This means the MC track must be a muon particle and leave at least one track
	/// reference point on either chamber in each tracking station.
	/// \param track  The MC track to check.
	/// \returns  true if the track is findable by the muon spectrometer
	///    and false otherwise.
	/// \note  The definition of track find-ability does not affect the correlation
	///    between ESD and MC tracks nor the generation of the AOD output tracks,
	///    only the histogram of findable tracks is effected if it is generated at all.

	// Select only muons.
	if (TMath::Abs(track->Particle()->GetPdgCode()) != 13) return false;
	
	// Build hit mask for chambers.
	bool hitOnCh[14];
	for (Int_t i = 0; i < 14; i++)
	{
		hitOnCh[i] = false;
	}
	for (Int_t i = 0; i < track->GetNumberOfTrackReferences(); i++)
	{
		AliTrackReference* ref = track->GetTrackReference(i);
		if (ref->DetectorId() != AliTrackReference::kMUON) continue;
		Int_t chamberId = ref->UserId() / 100 - 1;
		if (0 <= chamberId and chamberId < 14)
		{
			hitOnCh[chamberId] = true;
		}
	}
	
	// Check that enough hits are available on all tracking chambers.
	bool hitOnSt1 = hitOnCh[0] or hitOnCh[1];
	bool hitOnSt2 = hitOnCh[2] or hitOnCh[3];
	bool hitOnSt3 = hitOnCh[4] or hitOnCh[5];
	bool hitOnSt4 = hitOnCh[6] or hitOnCh[7];
	bool hitOnSt5 = hitOnCh[8] or hitOnCh[9];
	
	return (hitOnSt1 and hitOnSt2 and hitOnSt3 and hitOnSt4 and hitOnSt5);
}


void AliAnalysisTaskLinkToMC::FillHistsFromMC()
{
	/// Fills the histogram containing findable tracks from the MC event information.
	/// Only MC tracks for which the IsFindable method returns true are added
	/// to the histogram.
	
	for (Int_t i = 0; i < MCEvent()->GetNumberOfTracks(); i++)
	{
		AliMCParticle* mcTrack = (AliMCParticle*) MCEvent()->GetTrack(i);
		
		// Select only reconstructible tracks to fill into findable hist.
		if (IsFindable(mcTrack))
		{
			fFindableHist->Fill(mcTrack->Pt(), mcTrack->Y());
		}
	}
}


void AliAnalysisTaskLinkToMC::FillHistsFromLinks(TMap& links)
{
	/// Fills the histograms containing found tracks and fake tracks.
	/// Found tracks are all those that have a ESD muon track in the ESD event
	/// and for which a corresponding MC track could be found.
	/// Fake tracks are ESD muon tracks for which no corresponding MC track
	/// could be matched so it is considered a fake track.
	/// \param links  This contains the mapping between ESD muon tracks and
	///    MC tracks. ESD tracks for which no MC track could be matched should
	///    have the "value" of the key/value pair in the mapping set to NULL.
	
	TIter itmap(links.GetTable());
	TPair* pair = NULL;
	while ( (pair = static_cast<TPair*>(itmap())) != NULL )
	{
		if (pair->Value() != NULL)
		{
			AliESDMuonTrack* esdTrack = static_cast<AliESDMuonTrack*>(pair->Key());
			fFoundHist->Fill(esdTrack->Pt(), esdTrack->Y());
			
			AliMCParticle* mcTrack = static_cast<AliMCParticle*>(pair->Value());
			fFoundHistMC->Fill(mcTrack->Pt(), mcTrack->Y());
		}
		else
		{
			AliESDMuonTrack* esdTrack = static_cast<AliESDMuonTrack*>(pair->Key());
			fFakeHist->Fill(esdTrack->Pt(), esdTrack->Y());
		}
	}
}


void AliAnalysisTaskLinkToMC::CreateAODTracks(TMap& links)
{
	/// This code is copied from AliAnalysisTaskESDMuonFilter::ConvertESDtoAOD with
	/// modifications to fill the MC track label using the ESD to MC track mapping.
	/// \param links  This is the mapping between ESD tracks and MC tracks.
	
	// ESD Muon Filter analysis task executed for each event
	AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
        // CHECK
        if ( ! esd ) {
          AliError("Cannot get input event");
          return;
        }  
	
	// Define arrays for muons
	Double_t pos[3];
	Double_t p[3];
	Double_t pid[10];
	
	// has to be changed once the muon pid is provided by the ESD
	for (Int_t i = 0; i < 10; pid[i++] = 0.) {}
	pid[AliAODTrack::kMuon]=1.;
	
	AliAODHeader* header = AODEvent()->GetHeader();
	AliAODTrack *aodTrack = 0x0;
	AliESDMuonTrack *esdMuTrack = 0x0;
	
	// Access to the AOD container of tracks
	TClonesArray &tracks = *(AODEvent()->GetTracks());
	Int_t jTracks = tracks.GetEntriesFast();
	
	// Read primary vertex from AOD event 
	AliAODVertex *primary = AODEvent()->GetPrimaryVertex();
	if (primary != NULL) primary->Print();
	
	// Loop on muon tracks to fill the AOD track branch
	Int_t nMuTracks = esd->GetNumberOfMuonTracks();
	printf("Number of Muon Tracks=%d\n",nMuTracks);
	
	// Update number of positive and negative tracks from AOD event (M.G.)
	Int_t nPosTracks = header->GetRefMultiplicityPos();
	Int_t nNegTracks = header->GetRefMultiplicityNeg();
	
	for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack)
	{
		esdMuTrack = esd->GetMuonTrack(nMuTrack);
		
		//if (!esdMuTrack->ContainTrackerData()) continue;
		
		UInt_t selectInfo = 0;
		
		p[0] = esdMuTrack->Px(); 
		p[1] = esdMuTrack->Py(); 
		p[2] = esdMuTrack->Pz();
		
		pos[0] = esdMuTrack->GetNonBendingCoor(); 
		pos[1] = esdMuTrack->GetBendingCoor(); 
		pos[2] = esdMuTrack->GetZ();

		Int_t label = -1;
		TPair* pair = static_cast<TPair*>( links.FindObject(esdMuTrack) );
		if (pair != NULL and pair->Value() != NULL)
		{
			AliMCParticle* mcTrack = static_cast<AliMCParticle*>(pair->Value());
			label = mcTrack->GetLabel();
			if (label == -1)
			{
				AliWarning("The MC track was set to -1.");
			}
		}
		
		aodTrack = new(tracks[jTracks++]) AliAODTrack(
				esdMuTrack->GetUniqueID(), // ID
				label, // label
				p, // momentum
				kTRUE, // cartesian coordinate system
				pos, // position
				kFALSE, // isDCA
				0x0, // covariance matrix
				esdMuTrack->Charge(), // charge
				0, // ITSClusterMap
				pid, // pid
				primary, // primary vertex
				kFALSE, // used for vertex fit?
				kFALSE, // used for primary vertex fit?
				AliAODTrack::kPrimary,// track type
				selectInfo
			);
		
		aodTrack->SetXYAtDCA(esdMuTrack->GetNonBendingCoorAtDCA(), esdMuTrack->GetBendingCoorAtDCA());
		aodTrack->SetPxPyPzAtDCA(esdMuTrack->PxAtDCA(), esdMuTrack->PyAtDCA(), esdMuTrack->PzAtDCA());
		aodTrack->ConvertAliPIDtoAODPID();
		aodTrack->SetChi2perNDF(esdMuTrack->GetChi2() / (2.*esdMuTrack->GetNHit() - 5.));
		aodTrack->SetChi2MatchTrigger(esdMuTrack->GetChi2MatchTrigger());
		aodTrack->SetHitsPatternInTrigCh(esdMuTrack->GetHitsPatternInTrigCh());
		aodTrack->SetMuonClusterMap(esdMuTrack->GetMuonClusterMap());
		aodTrack->SetMatchTrigger(esdMuTrack->GetMatchTrigger());
		
		if (primary != NULL) primary->AddDaughter(aodTrack);
		
		if (esdMuTrack->Charge() > 0) nPosTracks++;
		else nNegTracks++;
	}
	
	header->SetRefMultiplicity(jTracks); 
	header->SetRefMultiplicityPos(nPosTracks);
	header->SetRefMultiplicityNeg(nNegTracks);
	
	// Note: we do not have to call PostData for the AOD because the AliAnalysisTaskSE
	// base class already does so just after the UserExec.
}


void AliAnalysisTaskLinkToMC::FindLinksToMC(TMap& links)
{
	/// Finds all MC tracks that correspond to their ESD tracks and creates a
	/// mapping between the ESD tracks and MC tracks which is filled into the
	/// 'links' TMap object.
	/// \param links  This will be filled with the ESD and MC track pairs which
	///    have been matched. Thus, for each ESD track which is set as the TMap key,
	///    the corresponding MC track is identified by the TMap value. If no
	///    MC track could be matched then the MC track value is set to NULL.

	// Algorithm description for connecting reconstructed ESD muon tracks to MC:
	//
	// Build ESD track list
	// Build MC track list
	// while ESD track list not empty and MC track list not empty:
	//   Find ESD/MC track pair with best compatibility.
	//   if no pair is compatible then break loop.
	//   Store found pair in return list.
	//   Remove pair from ESD and MC track list.
	// Mark remaining tracks in ESD track list as fakes.
	//
	// To find ESD/MC track pair with best compatibility:
	//   bestQuality := 1e100
	//   bestESDTrack := nil
	//   bestMCTrack := nil
	//   for each ESD track:
	//     for each MC track:
	//       chi^2 / numOfClustersMatched := calculate compatibility
	//       if tracks are compatible and chi^2 / numOfClustersMatched < bestQuality then:
	//         bestQuality := chi^2 / noOfClustersMatched
	//         bestESDTrack := current ESD track
	//         bestMCTrack := current MC track
	//
	// To calculate compatibility (ESD track, MC track):
	//   chi^2 := 1e100
	//   numberOfClustersAssociated := 0
	//   for each cluster in ESD track:
	//     find nearest 2 track refs from MC track
	//     if 2 track refs found then find average of them into MChit
	//     else MChit := nearest ref track.
	//     k := calculate chi^2 distance between ESD cluster and MChit
	//     if k <= maxValue then:
	//       chi^2 := chi^2 + k
	//       numberOfClustersAssociated := numberOfClustersAssociated + 1
	//     else:
	//       chi^2 := chi^2 + maxAllowedChi2
	//   if numberOfClustersAssociated >= minGoodClusters then:
	//     return (compatible, chi^2)
	//   else:
	//     return (not compatible, 1e100)
	
	AliESDEvent* esd = static_cast<AliESDEvent*>(InputEvent());

	TObjArray esdTracks;
	TObjArray mcTracks;
	for (Int_t i = 0; i < esd->GetNumberOfMuonTracks(); i++)
	{
		esdTracks.Add(esd->GetMuonTrack(i));
	}
	for (Int_t i = 0; i < MCEvent()->GetNumberOfTracks(); i++)
	{
		mcTracks.Add(MCEvent()->GetTrack(i));
	}
	
	while (not esdTracks.IsEmpty() and not mcTracks.IsEmpty())
	{
		Double_t bestQuality = 1e100;
		Int_t bestNoOfClusters = 0;
		AliESDMuonTrack* bestEsdTrack = NULL;
		AliMCParticle* bestMcTrack = NULL;
		TIter itesd(&esdTracks);
		TIter itmc(&mcTracks);
		AliESDMuonTrack* esdTrack = NULL;
		AliMCParticle* mcTrack = NULL;
		
		// Find the ESD and MC track pair that matches the best.
		while ( (esdTrack = static_cast<AliESDMuonTrack*>(itesd())) != NULL )
		{
			while ( (mcTrack = static_cast<AliMCParticle*>(itmc())) != NULL )
			{
				Double_t chi2 = 1e100;
				bool tracksCompatible = false;
				Int_t noOfClusters = 0;
				CalculateTrackCompatibility(esdTrack, mcTrack, tracksCompatible, chi2, noOfClusters);
				Double_t quality = noOfClusters > 0 ? chi2 / Double_t(noOfClusters) : 1e100;
				if (tracksCompatible and quality < bestQuality and noOfClusters >= bestNoOfClusters)
				{
					bestQuality = quality;
					bestNoOfClusters = noOfClusters;
					bestEsdTrack = esdTrack;
					bestMcTrack = mcTrack;
				}
			}
			itmc.Reset();
		}
		
		// If no track pair could be matched then we are done because we
		// will not be able to match anything more next time.
		if (bestMcTrack == NULL) break;
		
		links.Add(bestEsdTrack, bestMcTrack);
		
		esdTracks.Remove(bestEsdTrack);
		mcTracks.Remove(bestMcTrack);
	}
	
	// Add all remaining ESD tracks which could not be matched to any MC track
	// to the mapping but with NULL for the MC track.
	TIter itesd(&esdTracks);
	AliESDMuonTrack* esdTrack = NULL;
	while ( (esdTrack = static_cast<AliESDMuonTrack*>(itesd())) != NULL )
	{
		links.Add(esdTrack, NULL);
	}
}


void AliAnalysisTaskLinkToMC::CalculateTrackCompatibility(
		AliESDMuonTrack* esdTrack, AliMCParticle* mcTrack,
		bool& tracksCompatible, Double_t& chi2, Int_t& noOfClusters
	)
{
	/// Calculates the compatibility between a ESD and MC track pair.
	/// For each cluster of the ESD track, two track references are found
	/// which are spatially closest to it. If the distance between the cluster
	/// and track reference is larger than the HardCutLimit*() parameters then
	/// the track reference is not used. The average coordinate is calculated
	/// from the two track references found and the chi-square between the cluster
	/// and average coordinate is calculated. If only one track reference is
	/// found then it is used for the calculation without taking an average.
	/// If no track references are found then the cluster is not matched and
	/// we continue to the next cluster.
	/// If the chi-squared value for the cluster / track reference pair is
	/// larger than that allowed by the SigmaCut() parameter then the cluster
	/// is also considered not to be matched.
	/// The tracks are compatible if the ESD track has a minimum number of
	/// clusters matched with MC track references from the MC track.
	/// [in] \param esdTrack  The ESD track we are trying to match.
	/// [in] \param mcTrack  The MC track we are trying to match the ESD track to.
	/// [out] \param tracksCompatible  Filled with true if the given esdTrack
	///     and mcTrack are compatible. 
	/// [out] \param chi2  Filled with the total chi-squared calculated for all
	///     matched clusters.
	/// [out] \param noOfClusters  Filled with the number of clusters that were
	///     matched to at least one track reference.
	
	chi2 = 0;
	noOfClusters = 0;
	tracksCompatible = false;
	Int_t noOfClustersSt12 = 0;
	Int_t noOfClustersSt45 = 0;
	bool chamberMatched[14] = {
			false, false, false, false, false, false, false,
			false, false, false, false, false, false, false
		};
	bool stationMatched[7] = {
			false, false, false, false, false, false, false
		};
	
	for (Int_t i = 0; i < esdTrack->GetNClusters(); i++)
	{
		AliESDMuonCluster* cluster = static_cast<AliESDMuonCluster*>( esdTrack->GetClusters()[i] );
		Double_t varX = cluster->GetErrX2();
		Double_t varY = cluster->GetErrY2();
		// If the variance is zero then use the default one or just 1
		// if the default is also not set.
		if (varX == 0) varX = fVarX != 0 ? fVarX : 1.;
		if (varY == 0) varY = fVarY != 0 ? fVarY : 1.;
		
		Double_t chi2A = 1e100;
		AliTrackReference* refA = NULL;
		Double_t chi2B = 1e100;
		AliTrackReference* refB = NULL;
		
		for (Int_t j = 0; j < mcTrack->GetNumberOfTrackReferences(); j++)
		{
			AliTrackReference* ref = mcTrack->GetTrackReference(j);
			
			if (ref->DetectorId() != AliTrackReference::kMUON) continue;
			if (ref->UserId() != cluster->GetDetElemId()) continue;
			
			Double_t dx = ref->X() - cluster->GetX();
			if (TMath::Abs(dx) > fHardCutLimitX) continue;
			Double_t dy = ref->Y() - cluster->GetY();
			if (TMath::Abs(dy) > fHardCutLimitY) continue;
			if (fUseZCoordinate)
			{
				Double_t dz = ref->Z() - cluster->GetZ();
				if (TMath::Abs(dz) > fHardCutLimitZ) continue;
			}
			
			Double_t chi2sum = dx*dx / varX + dy*dy / varY;
			if (chi2sum < chi2A)
			{
				// Copy track reference A to B and set new value for A
				chi2B = chi2A;
				refB = refA;
				chi2A = chi2sum;
				refA = ref;
			}
			else if (chi2sum < chi2B)
			{
				// Track reference A still has smaller chi2 so just replace B.
				chi2B = chi2sum;
				refB = ref;
			}
		}
		
		Double_t x = 0, y = 0;
		if (refA != NULL and refB != NULL)
		{
			// Average the track references if 2 were found.
			x = (refA->X() + refB->X()) * 0.5;
			y = (refA->Y() + refB->Y()) * 0.5;
		}
		else if (refA != NULL)
		{
			x = refA->X();
			y = refA->Y();
		}
		else
		{
			continue;  // no track reference hit found.
		}
		
		Double_t dx = x - cluster->GetX();
		Double_t dy = y - cluster->GetY();
		
		Double_t chi2x = dx*dx / varX;
		if (chi2x > fSigmaCut2) continue;
		Double_t chi2y = dy*dy / varY;
		if (chi2y > fSigmaCut2) continue;
		
		chi2 += chi2x + chi2y;
		noOfClusters++;
		if (cluster->GetChamberId() >= 0 and cluster->GetChamberId() < 4)
		{
			noOfClustersSt12++;
		}
		if (cluster->GetChamberId() >= 6 and cluster->GetChamberId() < 10)
		{
			noOfClustersSt45++;
		}
		
		if (0 <= cluster->GetChamberId() and cluster->GetChamberId() < 14)
		{
			chamberMatched[cluster->GetChamberId()] = true;
			stationMatched[cluster->GetChamberId()/2] = true;
		}
	}
	
	// Need to check that all the chambers and stations that had to match we
	// actually found matching track references.
	
	bool forcedChambersMatched = true;
	for (Int_t i = 0; i < 14; i++)
	{
		if (fChamberMustMatch[i] and not chamberMatched[i])
		{
			forcedChambersMatched = false;
		}
	}
	bool forcedStationsMatched = true;
	for (Int_t i = 0; i < 7; i++)
	{
		if (fStationMustMatch[i] and not stationMatched[i])
		{
			forcedStationsMatched = false;
		}
	}
	
	if (noOfClusters >= fMinClusters and
	    noOfClustersSt12 >= fMinClustersInSt12 and
	    noOfClustersSt45 >= fMinClustersInSt45 and
	    forcedChambersMatched and forcedStationsMatched
	   )
	{
		tracksCompatible = true;
	}
}

