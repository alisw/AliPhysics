/*************************************************************************
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
// Apply special trackcuts to find out N-track events
// in central-productive events.
//
// Author:
//  Martin Poghosyan <Martin.Poghosyan@cern.ch>
//  continued by
//  Taesoo Kim <taesoo.kim@cern.ch>

#include "TMath.h"
#include "TObject.h"
#include "TIterator.h"
#include "TList.h"

#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliSPDUtils.h"
#include "AliITSsegmentationSPD.h"
#include "AliESDv0.h"

#include "AliMultiplicitySelectionCPPWA.h"


ClassImp(AliMultiplicitySelectionCPPWA)

//Class Constructor
AliMultiplicitySelectionCPPWA::AliMultiplicitySelectionCPPWA():TObject(),
	fkCheckReferenceMultiplicity(0)
{
	//TrackCuts
	fTrackCutListPrim = new TList();
	fTrackCutListPrim->SetOwner();
	fTrackCutListPrim->SetName("PrimaryTrackCut");

	//Set basic parameters
	SetTPCnclsS();
	SetTrackDCAz();
	SetTrackEtaRange();

	IgnoreV0s();
	for(Int_t i = 0; i< fkNtrackMax; i++)
		fkIsTrackSec[i]= kFALSE;
}
//Class Destructor
AliMultiplicitySelectionCPPWA::~AliMultiplicitySelectionCPPWA()
{
	if(fTrackCutListPrim)
	{
		fTrackCutListPrim->Delete();
		delete fTrackCutListPrim;
	}
	fTrackCutListPrim = 0;
}
//Member function
void AliMultiplicitySelectionCPPWA::InitDefaultTrackCuts(Int_t clusterCut, Bool_t ITSSACut, Bool_t IsRun2, Int_t nCluster)
{
	//Important message for 7TeV analysis (LHC10b,c,d,e)
	/*
	   Alexander Kalweit 
	   Email to PWG conveners on 22 Apr 2014
	   LHC10b&c (pass2):
	   ==================

	   Default cut which is currently recommended: AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE/kFALSE, 1)
	   Important is the second argument (=1) which replaces the cut on 70 clusters with a crossed rows cuts.  
	   !!! Please note, that a cut on 70 clusters is strongly discouraged in LHC10b&c pass2 data analysis!!! 
	   Changing to number of clusters (=0) and variations of the cut to 60 or 80 should be included in the systematic studies

	   LHC10deh (pass2):
	   ==================
	   Default cut which is currently recommended: AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE/kFALSE, 0)
	   In this period, a cut on 70 clusters should be okay, however, changing to a crossed rows cut and lowering the cut to 60 clusters should be included in the systematic error study.
	   */

	//Will be used as standard trackcuts
	if (IsRun2) {//Run2
		if (ITSSACut == kFALSE) {//ITS+TPC
			AliESDtrackCuts *fcutITSTPC_P = new AliESDtrackCuts;// = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, 0);
			//fcutITSTPC_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
			fcutITSTPC_P -> SetMaxDCAToVertexXYPtDep("(0.0182+0.0350/pt^1.01)");
			fcutITSTPC_P -> SetMinNCrossedRowsTPC(nCluster);
			fcutITSTPC_P -> SetMaxDCAToVertexZ(2);
			fcutITSTPC_P -> SetEtaRange(-0.9,0.9);
			fcutITSTPC_P -> SetMaxChi2PerClusterTPC(4);
			fcutITSTPC_P -> SetRequireTPCRefit(kTRUE);
			fcutITSTPC_P -> SetRequireITSRefit(kTRUE);
			fcutITSTPC_P -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
			fcutITSTPC_P -> SetAcceptKinkDaughters(kFALSE);
			fcutITSTPC_P -> SetMaxChi2PerClusterITS(36);
			fcutITSTPC_P -> SetMaxChi2TPCConstrainedGlobal(36);
			fcutITSTPC_P -> SetPtRange(0.15);
			fcutITSTPC_P -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fcutITSTPC_P -> SetMaxFractionSharedTPCClusters(0.4);
			fcutITSTPC_P->SetName("ITSTPC");
			AddPrimaryTrackCut(fcutITSTPC_P);
			AliESDtrackCuts *fcutITSSA_P = AliESDtrackCuts::GetStandardITSSATrackCuts2010(kTRUE, 0);
			fcutITSSA_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
			fcutITSSA_P->SetName("ITSSA");
			AddPrimaryTrackCut(fcutITSSA_P);
		}
		else {//ITSSA
			AliESDtrackCuts *fcutITSSA_P = AliESDtrackCuts::GetStandardITSSATrackCuts2010(kTRUE, 0);
			fcutITSSA_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
			fcutITSSA_P->SetName("ITSSA");
			AddPrimaryTrackCut(fcutITSSA_P);
		}
		return;
	}
	else {//Run1
		if (ITSSACut == kFALSE) {//ITS+TPC
			AliESDtrackCuts *fcutITSTPC_P = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, clusterCut);
			fcutITSTPC_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
			if (clusterCut == 1) {//For 10b,c
				fcutITSTPC_P->SetMinNCrossedRowsTPC(nCluster);
			}
			else {
				fcutITSTPC_P->SetMinNClustersTPC(nCluster);
			}

			fcutITSTPC_P->SetName("ITSTPC");
			AddPrimaryTrackCut(fcutITSTPC_P);
			AliESDtrackCuts *fcutITSSA_P = AliESDtrackCuts::GetStandardITSSATrackCuts2010(kTRUE, 0);
			fcutITSSA_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
			fcutITSSA_P->SetName("ITSSA");
			AddPrimaryTrackCut(fcutITSSA_P);
		}
		else {//ITSSA
			AliESDtrackCuts *fcutITSSA_P = AliESDtrackCuts::GetStandardITSSATrackCuts2010(kTRUE, 0);
			fcutITSSA_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
			fcutITSSA_P->SetName("ITSSA");
			AddPrimaryTrackCut(fcutITSSA_P);
		}
		return;
	}
}
//Member function
void AliMultiplicitySelectionCPPWA::AddPrimaryTrackCut(AliESDtrackCuts *cut)
{
	fTrackCutListPrim->Add(cut);
}
//Member function
Int_t AliMultiplicitySelectionCPPWA::GetNumberOfITSTPCtracks(AliESDEvent *esd)
{
	TArrayI indices;
	return GetNumberOfITSTPCtracks(esd, indices);
}
//Member function: Include V0 particles??
Bool_t AliMultiplicitySelectionCPPWA::InitV0Daughters(AliESDEvent *esd)
{
	if(fkNtrackMax < esd->GetNumberOfTracks() )
	{
		AliFatal(" fkNtrackMax < esd->GetNumberOfTracks() !!!\n");
	}

	for(Int_t i=0; i< esd->GetNumberOfTracks(); i++)
	{
		fkIsTrackSec[i] = kFALSE;
	}

	//  if(!fkIgnoreV0s) return kTRUE;

	Int_t Nv0  = esd->GetNumberOfV0s();

	for(Int_t iv0 = 0; iv0<Nv0; iv0++)
	{
		AliESDv0 *v0 = esd->GetV0(iv0);
		if(!v0) continue;

		fkIsTrackSec[v0->GetPindex()] = kTRUE;
		fkIsTrackSec[v0->GetNindex()] = kTRUE;
	}

	return kTRUE;
}
//Return Number of ITSTPC tracks with Standard + Special trackcuts
Int_t AliMultiplicitySelectionCPPWA::GetNumberOfITSTPCtracks(AliESDEvent *esd, TArrayI &indices)
{
	//Initialize TArrayI
	indices.Set(esd->GetNumberOfTracks());
	indices.Reset(-1);

	fIndicesN.Set(esd->GetNumberOfTracks());
	fIndicesN.Reset(-1);
	fIndicesP.Set(esd->GetNumberOfTracks());
	fIndicesP.Reset(-1);

	const AliESDVertex *vtxESD = esd->GetPrimaryVertex();

	Int_t NtracksSel = 0;
	Int_t NpureITStracks = 0;

	Int_t NtracksSelN = 0;
	Int_t NtracksSelP = 0;

	Double_t bfield = esd->GetMagneticField();
	Double_t dca[2], cov[3];

	esd->ConnectTracks();

	if(fkIgnoreV0s) InitV0Daughters(esd); 

	//Start trackcuts
	for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++)
	{
		AliESDtrack* track = esd->GetTrack(iTrack);
		track->SetESDEvent(esd);

		if(track->GetTPCnclsS()>fTPCnclsS) return -1;

		if(!track->PropagateToDCA(vtxESD, bfield, 500., dca, cov))
			continue;

		if(fkIgnoreV0s && fkIsTrackSec[iTrack])
			continue;

		Bool_t isITSpureSA = ((track->GetStatus() & AliESDtrack::kITSpureSA) != 0);
		if(isITSpureSA) 
		{
			NpureITStracks++;
			continue;
		}

		if(TMath::Abs(track->Zv() - vtxESD->GetZ())>fTrackDCAz) 
			continue;

		indices.AddAt(iTrack, NtracksSel);

		NtracksSel++;

		if(track->GetSign()<0)
		{
			fIndicesN.AddAt(iTrack, NtracksSelN);
			NtracksSelN++;
		}
		else if(track->GetSign()>0)
		{
			fIndicesP.AddAt(iTrack, NtracksSelP);
			NtracksSelP++;
		}


	}

	indices.Set(NtracksSel);

//	printf("NtracksSelN = %d   NtracksSelP = %d   ***************\n",NtracksSelN,NtracksSelP);

	fIndicesN.Set(NtracksSelN);
	fIndicesP.Set(NtracksSelP);

	for(Int_t i = 0; i < NtracksSel; i++)
	{
		AliESDtrack* tr = esd->GetTrack(indices.At(i));
		if(tr->Eta() < fTrackEtaMin || tr->Eta() > fTrackEtaMax)
			return -2;
		if(!AcceptTrack(tr, kTRUE))
			return -3;
	}

	const AliMultiplicity *mult = esd->GetMultiplicity();

	if(NpureITStracks>NtracksSel || mult->GetNumberOfTracklets() > NtracksSel)
		return -4;

	if(!TestFiredChips(esd, indices))
		return -5;

	if(fkCheckReferenceMultiplicity)
	{
		Int_t NRefMult = AliESDtrackCuts::GetReferenceMultiplicity(esd, AliESDtrackCuts::kTrackletsITSTPC, 3);

		if(NRefMult > NtracksSel)
			return -6;
	}

	return NtracksSel; 

}
//TODO :: Still under discussion
//Return Number of ITSSA tracks with Standard + Special trackcuts
Int_t AliMultiplicitySelectionCPPWA::GetNumberOfITSSAtracks(AliESDEvent *esd, TArrayI &indices)
{
	indices.Set(esd->GetNumberOfTracks());
	indices.Reset(-1);

	fIndicesN.Set(esd->GetNumberOfTracks());
	fIndicesN.Reset(-1);
	fIndicesP.Set(esd->GetNumberOfTracks());
	fIndicesP.Reset(-1);

	const AliESDVertex *vtxESD = esd->GetPrimaryVertex();

	Int_t NtracksSel = 0;
	Int_t NpureITStracks = 0;

	Int_t NtracksSelN = 0;
	Int_t NtracksSelP = 0;

	Double_t bfield = esd->GetMagneticField();
	Double_t dca[2], cov[3];


	esd->ConnectTracks();

	if(fkIgnoreV0s) InitV0Daughters(esd); 

	for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++)
	{
		AliESDtrack* track = esd->GetTrack(iTrack);
		track->SetESDEvent(esd);

		//     if(track->GetTPCnclsS()>fTPCnclsS) return -1;

		if(!track->PropagateToDCA(vtxESD, bfield, 500., dca, cov))
			continue;

		if(fkIgnoreV0s && fkIsTrackSec[iTrack])
			continue;

		Bool_t isITSpureSA = ((track->GetStatus() & AliESDtrack::kITSpureSA) != 0);
		if(isITSpureSA) 
		{
			NpureITStracks++;
			continue;
		}

		if(TMath::Abs(track->Zv() - vtxESD->GetZ())>fTrackDCAz) 
			continue;

		indices.AddAt(iTrack, NtracksSel);

		NtracksSel++;

		if(track->GetSign()<0)
		{
			fIndicesN.AddAt(iTrack, NtracksSelN);
			NtracksSelN++;
		}
		else if(track->GetSign()>0)
		{
			fIndicesP.AddAt(iTrack, NtracksSelP);
			NtracksSelP++;
		}
	}

	indices.Set(NtracksSel);

//	printf("NtracksSelN = %d   NtracksSelP = %d   ***************\n",NtracksSelN,NtracksSelP);

	fIndicesN.Set(NtracksSelN);
	fIndicesP.Set(NtracksSelP);

	for(Int_t i = 0; i< NtracksSel; i++)
	{
		AliESDtrack* tr = esd->GetTrack(indices.At(i));
		if(tr->Eta() < fTrackEtaMin || tr->Eta() > fTrackEtaMax)
			return -2;
		if(!AcceptTrack(tr, kTRUE))
			return -3;
	}

	const AliMultiplicity *mult = esd->GetMultiplicity();

	if(NpureITStracks>NtracksSel || mult->GetNumberOfTracklets() > NtracksSel)
		return -4;

	if(!TestFiredChips(esd, indices))
		return -5;

	if(fkCheckReferenceMultiplicity)
	{
		Int_t NRefMult = AliESDtrackCuts::GetReferenceMultiplicity(esd, AliESDtrackCuts::kTrackletsITSTPC, 3);

		if(NRefMult > NtracksSel)
			return -6;
	}

	return NtracksSel; 

}

Bool_t AliMultiplicitySelectionCPPWA::AcceptTrack(AliESDtrack *track, Bool_t asPrimary)
{
	if(asPrimary)
	{
		TIter next(fTrackCutListPrim);
		AliESDtrackCuts *cut;
		while ((cut=(AliESDtrackCuts*)next()))
		{
			if(cut->AcceptTrack(track))
				return kTRUE;
		}
		return kFALSE;
	}

	else{
		Bool_t isITSrefit = ((track->GetStatus() & AliESDtrack::kITSrefit) != 0);
		Bool_t isTPCrefit = ((track->GetStatus() & AliESDtrack::kTPCrefit) != 0);

		if(isITSrefit || isTPCrefit) return kTRUE;
		else return kFALSE;
	}

	return kFALSE;
}


Bool_t AliMultiplicitySelectionCPPWA::IsTrackSelected(Int_t index)
{

	for(Int_t i = 0; i< fIndicesN.GetSize(); i++)
	{
		if(fIndicesN.At(i)==index) 
			return kTRUE;
	}

	for(Int_t i = 0; i< fIndicesP.GetSize(); i++)
	{
		if(fIndicesP.At(i)==index) 
			return kTRUE;
	}

	return  kFALSE;
}


Bool_t AliMultiplicitySelectionCPPWA::TestFiredChips(AliESDEvent *esd, TArrayI indices)
{

	const AliMultiplicity *mult = esd->GetMultiplicity();
	Int_t Ntracks = indices.GetSize();
	UInt_t *Modules = new UInt_t[2*Ntracks];

	for(Int_t iT = 0; iT< Ntracks; iT++)
	{
//		printf("AliMultiplicitySelectionCPPWA::TestFiredChips:  indices.At(%d) = %d \n", iT, indices.At(iT));

		Int_t statusLay;
		Int_t idet = -1;
		Float_t xloc,zloc;
		AliESDtrack* track = esd->GetTrack(indices.At(iT));
		Bool_t retc=track->GetITSModuleIndexInfo(0,idet,statusLay,xloc,zloc);
		if(retc && statusLay!=5) Modules[2*iT] = idet;
		retc=track->GetITSModuleIndexInfo(1,idet,statusLay,xloc,zloc);
		if(retc && statusLay!=5) Modules[2*iT+1] = idet;
	}

	UInt_t eq, hs, chip;
	for (Int_t i=0; i<1200; i++)
	{
		if (!mult->TestFiredChipMap(i)) continue;
		AliSPDUtils::GetOnlineFromOfflineChipKey(i, eq, hs,  chip);
		UInt_t module = AliSPDUtils::GetOfflineModuleFromOnline(eq, hs, chip);

		Bool_t ktmp = kFALSE;
		for(Int_t iM = 0; iM<2*Ntracks; iM++)
		{
			if(Modules[iM]==module)
				ktmp=kTRUE;
		}
		if(!ktmp) 
		{
			delete[] Modules;
			return kFALSE;
		}
	}

	delete[] Modules;
	return kTRUE;
}
