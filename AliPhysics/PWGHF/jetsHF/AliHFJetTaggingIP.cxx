/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

// ************************************************************************
// Impact parameter based b-jet tagger class for HF jet analysis
// To be used with the Emcal Jet Framework
// Author: linus.feldkamp@cern.ch
// ***********************************************************************
#include "TF1.h"
#include "AliHFJetTaggingIP.h"
#include "AliEmcalJet.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "TRefArray.h"
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODv0.h"
#include "AliParticleContainer.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "TVector3.h"
#include "AliPicoTrack.h"
#include "AliExternalTrackParam.h"
// c std lib
#include <vector>
#include <utility>
#include <algorithm>
ClassImp(AliHFJetTaggingIP)
//_____________________________________________________________________________
AliHFJetTaggingIP::AliHFJetTaggingIP(void)
: fUseThresholdFuction(kFALSE)
, fAnaTypeAOD(kFALSE)
, fUseSignAtlas(kFALSE)
, fUseSignificance(kFALSE)
, fUse3DsIP(kFALSE)
, fCurrentDCA(0.)
, fThreshold(-99.)
, fDiscriminators()
, fTrackDCA(0)
, fJetIndices()
, fThresholdFuction(NULL)
, fEvent(NULL)
, fVertex(NULL)
, fVertexRecalculated(NULL)
, fJet(NULL)
, fParticles(NULL)
{
	//========================================================================
	// default constructor
	//========================================================================
	memset(fSelectionCuts, 0, sizeof fSelectionCuts);
	memset(fCurrentTrack, 0, sizeof fCurrentTrack);
	memset(fCurrentTrackDCAz, 0, sizeof fCurrentTrackDCAz);
	fAnaTypeAOD = kFALSE;
	this->InitTrackSelectionParams(0x0);
}
//_________________________________________________________________________
AliHFJetTaggingIP::~AliHFJetTaggingIP()
{
	//========================================================================'
	// default destructor
	//========================================================================
	if(this->fEvent)
		this->fEvent = 0x0; // delete pointers only since owned by Manager
	if(this->fVertex)
		this->fVertex = 0x0; //
	if(this->fJet)
		this->fJet = 0x0; //
}

Bool_t AliHFJetTaggingIP::DoTagging(AliEmcalJet* jetrec, Double_t& n0, Double_t& n1, Double_t& n2)
{
	if(!jetrec)
		return kFALSE;
	AliEmcalJet* jet = jetrec;
	Bool_t bTagged[3] = { kFALSE, kFALSE, kFALSE };
	Double_t bParam[3] = { -99., -99., -99 };
	if(!(this->GetJetDiscriminator(jet, bParam, bTagged)))
		return kFALSE;
	if(bTagged[2]) {
		if(this->fUseThresholdFuction) {
			if(bParam[2] >= this->fThresholdFuction->Eval(jet->Pt())) {
				return kTRUE;
			}
		} else {
			if(bParam[2] >= this->fThreshold) {
				n0 = bParam[0];
				n1 = bParam[1];
				n2 = bParam[2];
				return kTRUE;
			}
		}
	}
	return kFALSE;
}

Double_t CalculateTrackProbability(AliVTrack* bTrack)
{
	// Introduce definition of quality classes
	// To be checked in general track QA
	/*
       Possible criteria:
       Number of ITS and TPC clusters

       Quality Good:

       >=6 ITS clusters
       2   Hits in  SPD
       >   100 TPC clusters
       <   10 cm   d0 2d
       V0 Veto (>=2.6cm conversion point)
	 */
	//

	return 0.;
}

Bool_t AliHFJetTaggingIP::GetJetDiscriminator(AliEmcalJet* jet, Double_t* discriminator, Bool_t* check_discr)
{
	//========================================================================
	// Calculates overall discriminator
	//========================================================================
	std::vector<std::pair<Int_t, Double_t> > bSignedImpactParameter;
	if(!jet)
		return kFALSE;
	if(!fParticles)
		return kFALSE;
	fJet = jet;

	AliVTrack* bTrack = 0x0;
	fCurrentTrack[0] = -1;
	fCurrentTrack[1] = -1;
	fCurrentTrack[2] = -1;

	for(Int_t j = 0; j < fJet->GetNumberOfTracks(); ++j) {
		Double_t bSign = 0.;
		Double_t bIp2d = -999.;
		Double_t bDCAZ = -999.;

		bTrack = (AliVTrack*)((AliPicoTrack*)fParticles->GetParticle(jet->TrackAt((int)j)))->GetTrack();
		if(!bTrack)
			continue;
		if(!GetImpactParameter(bTrack, &bSign, &bIp2d, &bDCAZ))
			continue;
		if(!PassedCuts(bTrack, bIp2d, bDCAZ))
			continue;
		if(fVertexRecalculated) {
			delete fVertexRecalculated;
			fVertexRecalculated = NULL;
		}
		bSignedImpactParameter.push_back(std::make_pair(j, bSign * bIp2d));
		this->fTrackDCA.push_back(std::make_pair(j, bDCAZ));
	}
	discriminator[0] = -99.;
	discriminator[1] = -99.;
	discriminator[2] = -99.;
	Int_t numoftracks = (Int_t)bSignedImpactParameter.size();
	if(numoftracks < this->fSelectionCuts[AliHFJetTaggingIP::S_MINNTRACKS])
		return kFALSE;
	if(numoftracks == 0) {
		bSignedImpactParameter.clear();
		return kFALSE;
	}
	std::sort(bSignedImpactParameter.begin(), bSignedImpactParameter.end(), AliHFJetTaggingIP::mysort);
	if(numoftracks > 2) {
		discriminator[2] = bSignedImpactParameter.at(2).second;
		fCurrentTrack[2] = jet->TrackAt(bSignedImpactParameter.at(2).first);

		check_discr[2] = kTRUE;
	}
	if(numoftracks > 1) {
		discriminator[1] = bSignedImpactParameter.at(1).second;
		fCurrentTrack[1] = jet->TrackAt(bSignedImpactParameter.at(1).first);
		check_discr[1] = kTRUE;
	}
	if(numoftracks > 0) {
		discriminator[0] = bSignedImpactParameter.at(0).second;
		fCurrentTrack[0] = jet->TrackAt(bSignedImpactParameter.at(0).first);
		check_discr[0] = kTRUE;
	}
	SetCurrentDCAz();

	return kTRUE;
}

Bool_t AliHFJetTaggingIP::GetJetDiscriminatorQualityClass(Int_t qtyclass,
		AliEmcalJet* jet,
		Double_t* discriminator,
		Bool_t* check_discr)
{
	//========================================================================
	// Calculates overall discriminator
	//========================================================================
	std::vector<std::pair<Int_t, Double_t> > bSignedImpactParameter;
	if(!jet)
		return kFALSE;
	if(!fParticles)
		return kFALSE;

	fJet = jet;
	fCurrentTrack[0] = -1;
	fCurrentTrack[1] = -1;
	fCurrentTrack[2] = -1;

	AliVTrack* bTrack = 0x0;
	for(Int_t j = 0; j < fJet->GetNumberOfTracks(); ++j) {
		Double_t bSign = 0.;
		Double_t bIp2d = -999.;
		Double_t bDCAZ = -999.;
		bTrack = (AliVTrack*)((AliPicoTrack*)fParticles->GetParticle(jet->TrackAt((int)j)))->GetTrack();
		if(!bTrack)	continue;
		if(!GetImpactParameter(bTrack, &bSign, &bIp2d, &bDCAZ)) continue;
		if(this->fVertexRecalculated) {
			delete fVertexRecalculated;
			fVertexRecalculated = NULL;
		}
		if(!IsInQualityClass((AliAODTrack*)bTrack, qtyclass,bIp2d,bDCAZ)) continue;
		bSignedImpactParameter.push_back(std::make_pair(j, bSign * bIp2d));
		this->fTrackDCA.push_back(std::make_pair(j, bDCAZ));
		//	Printf("%i , pt %f" , j,bTrack->Pt());
	}
	//	Printf("Done with loop\n");
	discriminator[0] = -99.;
	discriminator[1] = -99.;
	discriminator[2] = -99.;
	Int_t numoftracks = (Int_t)bSignedImpactParameter.size();
	if(numoftracks < this->fSelectionCuts[AliHFJetTaggingIP::S_MINNTRACKS])
		return kFALSE;
	if(numoftracks == 0) {
		bSignedImpactParameter.clear();
		return kFALSE;
	}
	std::sort(bSignedImpactParameter.begin(), bSignedImpactParameter.end(), AliHFJetTaggingIP::mysort);
	if(numoftracks > 2) {
		discriminator[2] = bSignedImpactParameter.at(2).second;
		fCurrentTrack[2] = jet->TrackAt((int)bSignedImpactParameter.at(2).first);
		check_discr[2] = kTRUE;
	//	Printf("--> %i , pt %f" , 2,this->GetCurrentTrack(2)->Pt());

	}
	if(numoftracks > 1) {
		discriminator[1] = bSignedImpactParameter.at(1).second;
		fCurrentTrack[1] = jet->TrackAt((int)bSignedImpactParameter.at(1).first);
		check_discr[1] = kTRUE;
		//	Printf("--> %i , pt %f" , 1,this->GetCurrentTrack(1)->Pt());

	}
	if(numoftracks > 0) {
		discriminator[0] = bSignedImpactParameter.at(0).second;
		fCurrentTrack[0] = jet->TrackAt((int)bSignedImpactParameter.at(0).first);
		check_discr[0] = kTRUE;
		//	Printf("--> %i , pt %f" , 0,this->GetCurrentTrack(0)->Pt());

	}
	SetCurrentDCAz();



	return kTRUE;
}

Bool_t AliHFJetTaggingIP::GetImpactParameter(AliVTrack* bTrack, Double_t* bSign, Double_t* bIp2d, Double_t* zDCA)
{
	//========================================================================
	// Calculates the 2d impact parameter significance using the re-calculated event vertex
	//========================================================================
	if(!bTrack || !this->fEvent || !this->fJet)
		return kFALSE;
	Int_t bSkipped[2];
	Float_t bDiamondcovxy[3];
	Double_t bPosAtDCA[2] = { -999, -999 };
	Double_t bCovar[3] = { -999, -999, -999 };
	Double_t bpV[3] = { 0., 0., 0. };
	Double_t bpTrack[3] = { 0., 0., 0. };
	Double_t bpTrackP[3] = { 0., 0., 0. };

	const Double_t kBeampiperadius = 2.6;

	this->fVertex = (AliVVertex*)this->fEvent->GetPrimaryVertex();
	AliVertexerTracks* bVertexer = new AliVertexerTracks(fEvent->GetMagneticField());

	bVertexer->SetITSMode();
	bVertexer->SetMinClusters(4);
	bSkipped[0] = bTrack->GetID();
	bVertexer->SetSkipTracks(1, bSkipped);
	bVertexer->SetConstraintOn();
	fEvent->GetDiamondCovXY(bDiamondcovxy);

	Double_t bpos[3] = { this->fEvent->GetDiamondX(), this->fEvent->GetDiamondY(), 0. };
	Double_t bcov[6] = { bDiamondcovxy[0], bDiamondcovxy[1], bDiamondcovxy[2], 0., 0., 10. };
	AliESDVertex* bDiamond = new AliESDVertex(bpos, bcov, 1., 1);

	bVertexer->SetVtxStart(bDiamond);

	this->fVertexRecalculated = bVertexer->FindPrimaryVertex(fEvent);
	delete bDiamond;
	bDiamond = 0x0;
	delete bVertexer;
	bVertexer = NULL;
	if(this->fVertexRecalculated)
		this->fVertex = this->fVertexRecalculated;

	AliExternalTrackParam betp;
	betp.CopyFromVTrack(bTrack);
	if(!(betp.PropagateToDCA(fVertex, fEvent->GetMagneticField(), kBeampiperadius, bPosAtDCA, bCovar))) {
		if(this->fVertexRecalculated) {
			delete this->fVertexRecalculated;
			this->fVertexRecalculated = NULL;
		}
		return kFALSE;
	}
	fVertex->GetXYZ(bpV);
	betp.GetXYZ(bpTrack);
	betp.GetPxPyPz(bpTrackP);

	Double_t bIPVector[3] = { bpTrack[0] - bpV[0], bpTrack[1] - bpV[1], bpTrack[2] - bpV[2] };
	Double_t absIP =
			sqrt((bpTrack[0] - bpV[0]) * (bpTrack[0] - bpV[0]) + (bpTrack[1] - bpV[1]) * (bpTrack[1] - bpV[1]) +
					(bpTrack[2] - bpV[2]) * (bpTrack[2] - bpV[2]));
	Double_t bVar =
			(bIPVector[0] * this->fJet->Px() + bIPVector[1] * this->fJet->Py() + bIPVector[2] * this->fJet->Pz()) /
			(absIP * this->fJet->P());
	*zDCA = bPosAtDCA[1];
	bVar >= 0 ? *bSign = 1. : *bSign = -1.;

	*bSign = bVar;
	Double_t ptrIP = fabs(bPosAtDCA[0]);
	if(fUse3DsIP) {
		ptrIP = TMath::Sqrt(bPosAtDCA[0] * bPosAtDCA[0] + bPosAtDCA[1] * bPosAtDCA[1]);
	}
	Double_t pJetArray[3];
	fJet->PxPyPz(pJetArray);
	if(fUseSignAtlas) {
		*bSign = GetSignAtlasDefinition(bpTrack, bpTrackP, bpV, pJetArray);
	}
	*bIp2d = ptrIP;

	if(fUseSignificance) {
		if(!fUse3DsIP)
			*bIp2d = ptrIP / bCovar[0];
		else {
			double v3derror = 0.;
			v3derror = bPosAtDCA[0] * (bCovar[0] * bCovar[0]) + bPosAtDCA[1] * (bCovar[2] * bCovar[2]) +
					2 * bPosAtDCA[0] * bPosAtDCA[1] * bCovar[1];
			*bIp2d = ptrIP / v3derror;
		}
	}

	return kTRUE;
}
Double_t AliHFJetTaggingIP::GetDecayLength(AliVTrack* bTrack)
{
	//========================================================================
	// Calculates decay length i.e. distance from primary vertex to the position
	// of closest aproach to the jet and the distance in that point
	//========================================================================
	if(!bTrack)
		return -1.;
	Double_t bcv[21] = { 0 };
	Double_t bpxpypz[3] = { fJet->Px(), fJet->Py(), fJet->Pz() };
	Double_t bpos[3] = { 0x0 };
	Double_t xa = 0., xb = 0.;
	Double_t xyz[3] = { 0., 0., 0. };
	Double_t xyzb[3] = { 0., 0., 0. };
	this->fVertex->GetXYZ(bpos);
	AliExternalTrackParam bjetparam(bpos, bpxpypz, bcv, (Short_t)0);
	AliExternalTrackParam betp;
	betp.CopyFromVTrack(bTrack);
	if(!bTrack || !(this->fEvent->GetMagneticField()))
		return -1.;
	this->fCurrentDCA = bjetparam.GetDCA(&betp, this->fEvent->GetMagneticField(), xa, xb);
	bjetparam.GetXYZAt(xa, this->fEvent->GetMagneticField(), xyz);
	betp.GetXYZAt(xb, this->fEvent->GetMagneticField(), xyzb);

	Double_t bdecaylength =
			TMath::Sqrt((bpos[0] - xyz[0]) * (bpos[0] - xyz[0]) + (bpos[1] - xyz[1]) * (bpos[1] - xyz[1]) +
					(bpos[2] - xyz[2]) * (bpos[2] - xyz[2]));
	if(bdecaylength > 0)
		return bdecaylength;
	return -1.;
}
Bool_t AliHFJetTaggingIP::PassedCuts(AliVTrack* bTrack, Double_t ip, Double_t ipz)
{
	if(!bTrack)
		return kFALSE;
	if(bTrack->Pt() < this->fSelectionCuts[AliHFJetTaggingIP::S_PTTRACK])
		return kFALSE;

	UInt_t status = ((AliVTrack*)bTrack)->GetStatus();
	if((status & AliAODTrack::kITSrefit) == 0)
		return kFALSE;
	if((status & AliAODTrack::kTPCrefit) == 0)
		return kFALSE;
	if(fAnaTypeAOD) {
		if(!((((AliAODTrack*)bTrack)->HasPointOnITSLayer(0)) || (((AliAODTrack*)bTrack)->HasPointOnITSLayer(1))))
			return kFALSE;
	} else {
		if(!((((AliESDtrack*)bTrack)->HasPointOnITSLayer(0)) || (((AliESDtrack*)bTrack)->HasPointOnITSLayer(1))))
			return kFALSE;
	}
	if(fAnaTypeAOD) {
		if(((AliAODTrack*)bTrack)->GetITSNcls() < (int)this->fSelectionCuts[AliHFJetTaggingIP::S_ITSNCLS])
			return kFALSE;
	} else {
		if(((AliESDtrack*)bTrack)->GetITSNcls() < (int)this->fSelectionCuts[AliHFJetTaggingIP::S_ITSNCLS])
			return kFALSE;
	}

	Double_t dl = this->GetDecayLength(bTrack);
	if(fabs(ip) > this->fSelectionCuts[AliHFJetTaggingIP::S_TRANSVERSEIP])
		return kFALSE;
	if(sqrt(dl * dl + this->fCurrentDCA * this->fCurrentDCA) > fSelectionCuts[AliHFJetTaggingIP::S_DECAYLENGTH])
		return kFALSE;
	if(dl > fSelectionCuts[AliHFJetTaggingIP::S_DECAYLENGTH])
		return kFALSE;
	if(this->fCurrentDCA > this->fSelectionCuts[AliHFJetTaggingIP::S_MAXDCAJETTRACK])
		return kFALSE;
	if(fabs(ipz) > this->fSelectionCuts[AliHFJetTaggingIP::S_MAXDCAZ])
		return kFALSE;


	return kTRUE;
}
bool AliHFJetTaggingIP::mysort(const std::pair<Int_t, Double_t>& i, const std::pair<Int_t, Double_t>& j)
{
	if(i.second <= j.second)
		return false;
	else
		return true;
}
void AliHFJetTaggingIP::InitTrackSelectionParams(Double_t* params)
{
	if(!params) {
		// Set to default values
		this->fSelectionCuts[AliHFJetTaggingIP::S_PTTRACK] = 1.0;   // GeV/c
		this->fSelectionCuts[AliHFJetTaggingIP::S_ITSNCLS] = 1.0;   // (int)
		this->fSelectionCuts[AliHFJetTaggingIP::S_TRACKCHI2] = 5.0; //
		this->fSelectionCuts[AliHFJetTaggingIP::S_MINNTRACKS] = 1;
		this->fSelectionCuts[AliHFJetTaggingIP::S_DECAYLENGTH] = 10.;     // cm (max)
		this->fSelectionCuts[AliHFJetTaggingIP::S_MAXDCAJETTRACK] = 0.07; // cm
		this->fSelectionCuts[AliHFJetTaggingIP::S_TRANSVERSEIP] = 0.4;    // cm
		this->fSelectionCuts[AliHFJetTaggingIP::S_MAXDCAZ] = 10.;         // cm
	} else {
		this->fSelectionCuts[AliHFJetTaggingIP::S_PTTRACK] = params[AliHFJetTaggingIP::S_PTTRACK];     // GeV/c
		this->fSelectionCuts[AliHFJetTaggingIP::S_ITSNCLS] = params[AliHFJetTaggingIP::S_ITSNCLS];     // (int)
		this->fSelectionCuts[AliHFJetTaggingIP::S_TRACKCHI2] = params[AliHFJetTaggingIP::S_TRACKCHI2]; //
		this->fSelectionCuts[AliHFJetTaggingIP::S_MINNTRACKS] = params[AliHFJetTaggingIP::S_MINNTRACKS];
		this->fSelectionCuts[AliHFJetTaggingIP::S_DECAYLENGTH] = params[AliHFJetTaggingIP::S_DECAYLENGTH]; // cm (max)
		this->fSelectionCuts[AliHFJetTaggingIP::S_MAXDCAJETTRACK] = params[AliHFJetTaggingIP::S_MAXDCAJETTRACK]; // cm
		this->fSelectionCuts[AliHFJetTaggingIP::S_TRANSVERSEIP] = params[AliHFJetTaggingIP::S_TRANSVERSEIP];     // cm
		this->fSelectionCuts[AliHFJetTaggingIP::S_MAXDCAZ] = params[AliHFJetTaggingIP::S_MAXDCAZ];
	}
}
void AliHFJetTaggingIP::SetImpactParameterThreshold(Double_t threshold)
{
	this->fThreshold = threshold;
}
void AliHFJetTaggingIP::SetImpactParameterThreshold(TF1* threshld_fct)
{
	this->fUseThresholdFuction = kTRUE;
	this->fThresholdFuction = threshld_fct;
}

Bool_t AliHFJetTaggingIP::IsV0DaughterRadius(AliVTrack* track, Double_t& Radius)
{
	AliESDv0* v0esd = 0x0;
	AliAODv0* v0aod = 0x0;
	for(int i = 0; i < fEvent->GetNumberOfV0s(); ++i) {
		if(!fAnaTypeAOD) {
			v0esd = ((AliESDEvent*)fEvent)->GetV0(i);
			int posid = v0esd->GetPindex();
			int negid = v0esd->GetNindex();
			int trackid = track->GetID();
			if(posid == trackid || negid == trackid) {
				Double_t P[3];
				v0esd->XvYvZv(P);
				Radius = sqrt(P[0] * P[0] + P[1] * P[1]);
				return kTRUE;
			}
		} else if(fAnaTypeAOD) {
			v0aod = ((AliAODEvent*)fEvent)->GetV0(i);
			int posid = v0aod->GetPosID();
			int negid = v0aod->GetNegID();
			int trackid = track->GetID();
			if(posid == trackid || negid == trackid) {
				Double_t P[3];
				P[0] = v0aod->DecayVertexV0X();
				P[1] = v0aod->DecayVertexV0Y();
				P[2] = v0aod->DecayVertexV0Z();
				Radius = sqrt(P[0] * P[0] + P[1] * P[1]);
				return kTRUE;
			}
		}
	}
	return kFALSE;
}

Bool_t AliHFJetTaggingIP::IsInQualityClass(AliAODTrack* track, int qclass, double ip, double ipz)
{

	if(!track)
		return kFALSE;
	if(track->Pt() < 1.)
		return kFALSE;
	ULong_t status = track->GetStatus();

	if(!(status & AliAODTrack::kTPCrefit))
		return kFALSE;
	if(!(status & AliAODTrack::kITSrefit))
		return kFALSE;

	int nSPDHits = 0;
	if(track->HasPointOnITSLayer(0))
		nSPDHits++;
	if(track->HasPointOnITSLayer(1))
		nSPDHits++;
	int nITSHits = nSPDHits;
	for(int j = 2; j < 6; ++j)
		if(track->HasPointOnITSLayer(j))
			nITSHits++;
	int nTPCcls = 0;
	nTPCcls = ((AliAODTrack*)track)->GetTPCNcls();
	Float_t cRatioTPC = track->GetTPCNclsF() > 0. ?
			static_cast<Float_t>(track->GetTPCNcls()) / static_cast<Float_t>(track->GetTPCNclsF()) :
			1.;
	Bool_t isV0Daughter = kFALSE;
	Double_t v0Radius = 0.;
	isV0Daughter = IsV0DaughterRadius(track, v0Radius);

	if(fabs(ipz) > this->fSelectionCuts[AliHFJetTaggingIP::S_MAXDCAZ])
		return kFALSE;
	if(fabs(ip) > this->fSelectionCuts[AliHFJetTaggingIP::S_TRANSVERSEIP])
		return kFALSE;

	switch(qclass) {
	case 1:
		if(nSPDHits < 2)
			return kFALSE;
		if(nITSHits < 4)
			return kFALSE;
		if(nTPCcls < 90)
			return kFALSE;
		if(cRatioTPC < 0.6)
			return kFALSE;
		if(isV0Daughter && v0Radius > 2.0)
			return kFALSE;
		break;
	case 2:
		if(nSPDHits < 1)
			return kFALSE;
		if(nITSHits < 4)
			return kFALSE;
		if(nTPCcls < 90)
			return kFALSE;
		if(cRatioTPC < 0.6)
			return kFALSE;
		if(isV0Daughter && v0Radius > 2.0)
			return kFALSE;
		break;
	case 3:
		if(nSPDHits < 1)
			return kFALSE;
		if(nITSHits < 3)
			return kFALSE;
		if(nTPCcls < 90)
			return kFALSE;
		if(cRatioTPC < 0.6)
			return kFALSE;
		break;
	case 4:
		if(nITSHits < 3)
			return kFALSE;
		if(nTPCcls < 80)
			return kFALSE;
		if(cRatioTPC < 0.6)
			return kFALSE;
		break;
	default:
		return kFALSE;
		break;
	}

	return kTRUE;
}
Double_t AliHFJetTaggingIP::GetSignAtlasDefinition(Double_t* xDCA, Double_t* pDCA, Double_t* xVtx, Double_t* pJet)
{
	// Get vectors to make calculations easier
	TVector3 vXAtdca(xDCA[0], xDCA[1], xDCA[2]);
	TVector3 vPAtdca(pDCA[0], pDCA[1], pDCA[2]);
	TVector3 pVJet(pJet[0], pJet[1], pJet[2]);
	TVector3 pPV(xVtx[0], xVtx[1], xVtx[2]);

	TVector3 tmpVec = pVJet.Cross(vPAtdca);
	TVector3 tmpVec2 = pPV - vXAtdca;
	TVector3 tmpVec3 = vPAtdca.Cross(tmpVec2);
	Double_t tmpVec4 = tmpVec.Dot(tmpVec3);
	Double_t N = TMath::Abs(tmpVec4);

	return tmpVec4 / N;
}

AliVTrack* AliHFJetTaggingIP::GetCurrentTrack(int i)
{
	AliVTrack* res = 0x0;
	res = (AliVTrack*)(((AliPicoTrack*)(fParticles->GetParticle(this->fCurrentTrack[i])))->GetTrack());
	return res;
}

Double_t AliHFJetTaggingIP::GetCurrentTrackDCAz(int i)
{
	if(i < 0 || i > 2)
		return -9999.;
	else
		return this->fCurrentTrackDCAz[i];
}

void AliHFJetTaggingIP::SetCurrentDCAz()
{
	for(int i = 0; i < 3; ++i) {
		int curIdx = this->fCurrentTrack[i];
		if(curIdx < 0)
			continue;
		for(int j = 0; j < (int)this->fTrackDCA.size(); ++j) {
			if(fTrackDCA.at(j).first == curIdx)
				this->fCurrentTrackDCAz[i] = fTrackDCA.at(j).second;
		}
	}
}
