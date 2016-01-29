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

/* $Id$ */

//-----------------------------------------------------------------------------
// Class AliMFTTrack
//-------------------
// Description of an ALICE Standalone MFT track
//-------------------
// Contact author: raphael.tieulent@cern.ch
//-----------------------------------------------------------------------------

#include "TObjArray.h"
#include "TMath.h"

#include "AliLog.h"

#include "AliMFTTrack.h"
#include "AliMFTTrackParam.h"
#include "AliMFTCATrack.h"
#include "AliMFTCACell.h"


/// \cond CLASSIMP
ClassImp(AliMFTTrack); // Class implementation in ROOT context
											 /// \endcond


//=============================================================================================

AliMFTTrack::AliMFTTrack():TObject(),
fChi2(0.),
fTrackParamAtCluster(NULL),
fTrackParamAtVertex(NULL),
fCATrack(NULL),
fTrackID(-1),
fP(0),
fTheta(0),
fPhi(0),
fPt(0)
{
	/// Default constructor
	
}

//=============================================================================================

AliMFTTrack::AliMFTTrack(AliMFTCATrack *catrack):TObject(),
fChi2(0.),
fTrackParamAtCluster(new TObjArray(10)),
fTrackParamAtVertex(NULL),
fCATrack(catrack),
fTrackID(-1),
fP(0),
fTheta(0),
fPhi(0),
fPt(0)
{
	catrack->Print();
	/// Constructor from a AliMFTCATrack
	fTrackParamAtCluster->SetOwner(kTRUE);
	
	Int_t nCells = catrack->GetNcells();
	// Create Empty Track Parameter objects
	Double_t *caHit;
	AliMFTCACell * caCell;
	AliInfo(Form("Nb of Cells in the track = %d ",nCells));
	for (Int_t iCell = 0 ; iCell < nCells  ; iCell++) {
		caCell = catrack->GetCell(iCell);
		caHit = caCell->GetHit2();
		
		AliMFTTrackParam trackParam;
		trackParam.SetInverseTransverseMomentum(1.e-6); // infinite momentum
		trackParam.SetClusterPos(caHit[0], caHit[1], caHit[2]);
		
		AddTrackParamAtCluster(trackParam);
		
		if(iCell==nCells-1) {
			AliMFTTrackParam trackParam2;
			caHit = caCell->GetHit1();
			trackParam2.SetInverseTransverseMomentum(1.e-6); // infinite momentum
			trackParam2.SetClusterPos(caHit[0], caHit[1], caHit[2]);
			
			AddTrackParamAtCluster(trackParam2);
		}
		
	}
	
	AliInfo(Form("Nb of Track param objects = %d ",fTrackParamAtCluster->GetEntries()));
}
//__________________________________________________________________________
AliMFTTrack::AliMFTTrack(const AliMFTTrack& track)
: TObject(track),
fChi2(track.fChi2),
fTrackParamAtCluster(NULL),
fTrackParamAtVertex(NULL),
fCATrack(track.fCATrack),
fTrackID(track.fTrackID),
fP(track.fP),
fTheta(track.fTheta),
fPhi(track.fPhi),
fPt(track.fPt)

{
	///copy constructor
	
	// necessary to make a copy of the objects and not only the pointers in TObjArray.
	if (track.fTrackParamAtCluster) {
		fTrackParamAtCluster = new TObjArray(track.fTrackParamAtCluster->GetSize());
		fTrackParamAtCluster->SetOwner(kTRUE);
		for (Int_t i = 0; i < track.GetNClusters(); i++)
			fTrackParamAtCluster->AddLast(new AliMFTTrackParam(*static_cast<AliMFTTrackParam*>(track.fTrackParamAtCluster->UncheckedAt(i))));
	}
	
	
	// copy track parameters at vertex if any
	if (track.fTrackParamAtVertex) fTrackParamAtVertex = new AliMFTTrackParam(*(track.fTrackParamAtVertex));
	
}
//__________________________________________________________________________
AliMFTTrack & AliMFTTrack::operator=(const AliMFTTrack& track)
{
	/// Asignment operator
	// check assignement to self
	if (this == &track)
		return *this;
	
	// base class assignement
	TObject::operator=(track);
	
	// clear memory
	Clear();
	
	// necessary to make a copy of the objects and not only the pointers in TObjArray
	if (track.fTrackParamAtCluster) {
		fTrackParamAtCluster = new TObjArray(track.fTrackParamAtCluster->GetSize());
		fTrackParamAtCluster->SetOwner(kTRUE);
		for (Int_t i = 0; i < track.GetNClusters(); i++)
			fTrackParamAtCluster->AddLast(new AliMFTTrackParam(*static_cast<AliMFTTrackParam*>(track.fTrackParamAtCluster->UncheckedAt(i))));
	}
	
	
	// copy track parameters at vertex if any
	if (track.fTrackParamAtVertex) {
		if (fTrackParamAtVertex) *fTrackParamAtVertex = *(track.fTrackParamAtVertex);
		else fTrackParamAtVertex = new AliMFTTrackParam(*(track.fTrackParamAtVertex));
	}
	
	fChi2         =  track.fChi2;
	fTrackID =  track.fTrackID;
	fCATrack = track.fCATrack;
	fP=track.fP;
	fTheta=track.fTheta;
	fPhi=track.fPhi;
	fPt=track.fPt;
	
	return *this;
}

//=============================================================================================


AliMFTTrack::~AliMFTTrack() {
	
	delete fTrackParamAtCluster;
	delete fTrackParamAtVertex;
	delete fCATrack;
	
}

//__________________________________________________________________________
void AliMFTTrack::AddTrackParamAtCluster(const AliMFTTrackParam &trackParam)
{
	
	//	// check whether track parameters are given at the correct cluster z position
	//	if (TMath::Abs(cluster.GetZ() - trackParam.GetZ())>1.e-5) {   // AU
	//		AliError("track parameters are given at a different z position than the one of the associated cluster");
	//		return;
	//	}
	
	// add parameters to the array of track parameters
	if (!fTrackParamAtCluster) {
		fTrackParamAtCluster = new TObjArray(10);
		fTrackParamAtCluster->SetOwner(kTRUE);
	}
	AliMFTTrackParam* trackParamAtCluster = new AliMFTTrackParam(trackParam);
	fTrackParamAtCluster->AddLast(trackParamAtCluster);
	
	// sort the array of track parameters
	//	fTrackParamAtCluster->Sort();
}

//__________________________________________________________________________
TObjArray* AliMFTTrack::GetTrackParamAtCluster() const
{
	/// return array of track parameters at cluster (create it if needed)
	if (!fTrackParamAtCluster) {
		fTrackParamAtCluster = new TObjArray(10);
		fTrackParamAtCluster->SetOwner(kTRUE);
	}
	return fTrackParamAtCluster;
}

//__________________________________________________________________________
Int_t AliMFTTrack::GetNDF() const
{
	/// return the number of degrees of freedom
	
	Int_t ndf = 2 * GetNClusters() - 5;
	return (ndf > 0) ? ndf : 0;
}

//__________________________________________________________________________
Double_t AliMFTTrack::GetNormalizedChi2() const
{
	/// return the chi2 value divided by the number of degrees of freedom (or FLT_MAX if ndf <= 0)
	
	Double_t ndf = (Double_t) GetNDF();
	return (ndf > 0.) ? fChi2 / ndf : 1e6;
}

//_____________________________________________-
void AliMFTTrack::Print(Option_t*) const
{
	/// Printing Track information
	
	AliInfo(Form("Chi2 = %.1f", fChi2));
	AliInfo(Form("Chi2/ndf = %.2f", GetNormalizedChi2()));
	AliInfo(Form("No.Clusters = %d", GetNClusters()));
	AliInfo(Form("MC label = %d", fTrackID));
	if (fTrackParamAtCluster) fTrackParamAtCluster->Last()->Print("FULL");
}

//__________________________________________________________________________
void AliMFTTrack::Clear(Option_t* /*opt*/)
{
	/// Clear arrays
	//	delete fTrackParamAtCluster; fTrackParamAtCluster = 0x0;
	//	delete fTrackParamAtVertex; fTrackParamAtVertex = 0x0;
}



