/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include "AliEmcalTriggerMakerPart.h"

ClassImp(PWG::EMCAL::TriggerPart::AliEmcalTriggerMakerPart)

using namespace PWG::EMCAL::TriggerPart;

/**
 * Constructor, initializing channel maps for EMCAL and DCAL-PHOS
 */
AliEmcalTriggerMakerPart::AliEmcalTriggerMakerPart() :
	TObject(),
	fJetTrigger(),
	fGammaTrigger(),
	fTriggerChannelsEMCAL(48, 64),
	fTriggerChannelsDCALPHOS(48, 40),
	fTriggerMapping(),
	fBadChannelsEMCAL(),
	fBadChannelsDCALPHOS(),
	fHasRun(false),
	fAcceptPHOSPatches(true),
	fGammaEMCAL(),
	fGammaDCALPHOS(),
	fJetEMCAL(),
	fJetDCALPHOS(),
	fJetEMCAL8x8(),
	fJetDCALPHOS8x8()
{
	fJetTrigger.SetTriggerSetup(&fTriggerSetup);
	fGammaTrigger.SetTriggerSetup(&fTriggerSetup);
}

/**
 * Destructor
 */
AliEmcalTriggerMakerPart::~AliEmcalTriggerMakerPart() {
}

/**
 * Reset channel maps
 */
void AliEmcalTriggerMakerPart::Reset() {
	fTriggerChannelsEMCAL.Reset();
	fTriggerChannelsDCALPHOS.Reset();
	fGammaEMCAL.clear();
	fGammaDCALPHOS.clear();
	fJetEMCAL.clear();
	fJetDCALPHOS.clear();
	fJetEMCAL8x8.clear();
	fJetDCALPHOS8x8.clear();
	fHasRun = false;
}

/**
 * Main function to reconstruct trigger patches in the EMCAL and in the DCAL-PHOS.
 * Patches are found using the trigger channels map, which has to be filled from outside.
 * The order in the output list is:
 * #- EMCAL gamma
 * #- DCAL-PHOS gamma
 * #- EMCAL jet
 * #- DCAL-PHOS jet
 * @return vector with all trigger patches.
 */

void AliEmcalTriggerMakerPart::FindPatches() {

	fGammaEMCAL     = fGammaTrigger.FindPatches 	(	&fTriggerChannelsEMCAL		);
	fGammaDCALPHOS  = fGammaTrigger.FindPatches 	(	&fTriggerChannelsDCALPHOS	);
	fJetEMCAL       = fJetTrigger.  FindPatches 	(	&fTriggerChannelsEMCAL		);
	fJetDCALPHOS    = fJetTrigger.  FindPatches 	(	&fTriggerChannelsDCALPHOS	);
	fJetEMCAL8x8    = fJetTrigger.  FindPatches8x8	(	&fTriggerChannelsEMCAL		);
	fJetDCALPHOS8x8 = fJetTrigger.  FindPatches8x8	(	&fTriggerChannelsDCALPHOS	);

	fHasRun = true;
}

std::vector<AliEmcalTriggerPartRawPatch> AliEmcalTriggerMakerPart::GetPatches(const int what) {
	if (fHasRun == false)
		FindPatches();

	std::vector<AliEmcalTriggerPartRawPatch> result;

	if (what == AliEmcalTriggerPartRawPatch::kAny || what == AliEmcalTriggerPartRawPatch::kEMCALpatchGA ) {
		for (std::vector<AliEmcalTriggerPartRawPatch>::iterator patchiter = fGammaEMCAL.begin(); patchiter != fGammaEMCAL.end(); patchiter++) {
			patchiter->SetPatchType(AliEmcalTriggerPartRawPatch::kEMCALpatch);
			result.push_back(*patchiter);
		}
	}

	if (what == AliEmcalTriggerPartRawPatch::kAny || what == AliEmcalTriggerPartRawPatch::kDCALpatchGA ) {
		for (std::vector<AliEmcalTriggerPartRawPatch>::iterator patchiter = fGammaDCALPHOS.begin(); patchiter != fGammaDCALPHOS.end(); patchiter++) {
			if(!fAcceptPHOSPatches && IsPHOSPatch(patchiter->GetColStart(), patchiter->GetRowStart(), 2)) continue;
			patchiter->SetPatchType(AliEmcalTriggerPartRawPatch::kDCALPHOSpatch);
			result.push_back(*patchiter);
		}
	}

	if (what == AliEmcalTriggerPartRawPatch::kAny || what == AliEmcalTriggerPartRawPatch::kEMCALpatchJE ) {
		for (std::vector<AliEmcalTriggerPartRawPatch>::iterator patchiter = fJetEMCAL.begin(); patchiter != fJetEMCAL.end(); patchiter++) {
			patchiter->SetPatchType(AliEmcalTriggerPartRawPatch::kEMCALpatch);
			result.push_back(*patchiter);
		}
	}

	if (what == AliEmcalTriggerPartRawPatch::kAny || what == AliEmcalTriggerPartRawPatch::kDCALpatchJE ) {
		for (std::vector<AliEmcalTriggerPartRawPatch>::iterator patchiter = fJetDCALPHOS.begin(); patchiter != fJetDCALPHOS.end(); patchiter++) {
			if(!fAcceptPHOSPatches && IsPHOSPatch(patchiter->GetColStart(), patchiter->GetRowStart(), 16)) continue;
			patchiter->SetPatchType(AliEmcalTriggerPartRawPatch::kDCALPHOSpatch);
			result.push_back(*patchiter);
		}
	}

	if (what == AliEmcalTriggerPartRawPatch::kAny || what == AliEmcalTriggerPartRawPatch::kEMCALpatchJE8x8 ) {
		for (std::vector<AliEmcalTriggerPartRawPatch>::iterator patchiter = fJetEMCAL8x8.begin(); patchiter != fJetEMCAL8x8.end(); patchiter++) {
			patchiter->SetPatchType(AliEmcalTriggerPartRawPatch::kEMCALpatch);
			result.push_back(*patchiter);
		}
	}

	if (what == AliEmcalTriggerPartRawPatch::kAny || what == AliEmcalTriggerPartRawPatch::kDCALpatchJE8x8 ) {
		for (std::vector<AliEmcalTriggerPartRawPatch>::iterator patchiter = fJetDCALPHOS8x8.begin(); patchiter != fJetDCALPHOS8x8.end(); patchiter++) {
			if(!fAcceptPHOSPatches && IsPHOSPatch(patchiter->GetColStart(), patchiter->GetRowStart(), 8)) continue;
			patchiter->SetPatchType(AliEmcalTriggerPartRawPatch::kDCALPHOSpatch);
			result.push_back(*patchiter);
		}
	}

	return result;
}

bool AliEmcalTriggerMakerPart::IsPHOSPatch(int col, int row, int size){
	return (col >= kMinEtaPHOS) && (col+size <kMaxEtaPHOS) && (row >= kMinRowPHOS) && (row+size < kMaxRowPHOS);
}

AliEmcalTriggerPartRawPatch AliEmcalTriggerMakerPart::GetMaxGammaEMCAL()
{
	if (fHasRun == false)
		FindPatches();
	if (fGammaEMCAL.empty())
		return AliEmcalTriggerPartRawPatch();
	else
		return fGammaEMCAL.back();
}

AliEmcalTriggerPartRawPatch AliEmcalTriggerMakerPart::GetMaxGammaDCALPHOS()
{
	if (fHasRun == false)
		FindPatches();
	if (fGammaDCALPHOS.empty())
		return AliEmcalTriggerPartRawPatch();
	else
		return fGammaDCALPHOS.back();
}

AliEmcalTriggerPartRawPatch AliEmcalTriggerMakerPart::GetMaxJetEMCAL()
{
	if (fHasRun == false)
		FindPatches();
	if (fJetEMCAL.empty())
		return AliEmcalTriggerPartRawPatch();
	else
		return fJetEMCAL.back();
}

AliEmcalTriggerPartRawPatch AliEmcalTriggerMakerPart::GetMaxJetDCALPHOS()
{
	if (fHasRun == false)
		FindPatches();
	if (fJetDCALPHOS.empty())
		return AliEmcalTriggerPartRawPatch();
	else
		return fJetDCALPHOS.back();
}

AliEmcalTriggerPartRawPatch AliEmcalTriggerMakerPart::GetMaxJetEMCAL8x8()
{
	if (fHasRun == false)
		FindPatches();
	if (fJetEMCAL8x8.empty())
		return AliEmcalTriggerPartRawPatch();
	else
		return fJetEMCAL8x8.back();
}

AliEmcalTriggerPartRawPatch AliEmcalTriggerMakerPart::GetMaxJetDCALPHOS8x8()
{
	if (fHasRun == false)
		FindPatches();
	if (fJetDCALPHOS8x8.empty())
		return AliEmcalTriggerPartRawPatch();
	else
		return fJetDCALPHOS8x8.back();
}

double AliEmcalTriggerMakerPart::GetMedian(std::vector<AliEmcalTriggerPartRawPatch> v)
{
	double median = 0;
	size_t size = v.size();

	if (size > 0)
	{
		size_t halfsize = v.size() / 2;
		if (size % 2 == 0)
		{
			median = (v[halfsize - 1].GetADC() + v[halfsize].GetADC()) / 2;
		}
		else
		{
			median = v[halfsize].GetADC();
		}
	}
	return median;
}

double AliEmcalTriggerMakerPart::GetMedianGammaEMCAL()
{
	if (fHasRun == false)
		FindPatches();
	return GetMedian(fGammaEMCAL);
}

double AliEmcalTriggerMakerPart::GetMedianGammaDCALPHOS()
{
	if (fHasRun == false)
		FindPatches();
	return GetMedian(fGammaDCALPHOS);
}

double AliEmcalTriggerMakerPart::GetMedianJetEMCAL()
{
	if (fHasRun == false)
		FindPatches();
	return GetMedian(fJetEMCAL);
}

double AliEmcalTriggerMakerPart::GetMedianJetDCALPHOS()
{
	if (fHasRun == false)
		FindPatches();
	return GetMedian(fJetDCALPHOS);
}

double AliEmcalTriggerMakerPart::GetMedianJetEMCAL8x8()
{
	if (fHasRun == false)
		FindPatches();
	return GetMedian(fJetEMCAL8x8);
}

double AliEmcalTriggerMakerPart::GetMedianJetDCALPHOS8x8()
{
	if (fHasRun == false)
		FindPatches();
	return GetMedian(fJetDCALPHOS8x8);
}

/**
 * Fill trigger channel map depending on where the particle hits the detector in the EMCAL DCAL-PHOS surface. Adds the charge
 * to the already existing charge. Doesn't do anything if the particle is outside of the detector acceptance of either of the
 * tow subsystems.
 * @param eta Track/Particle eta
 * @param phi Track/Particle phi
 * @param energy Track/Particle energy
 */
void AliEmcalTriggerMakerPart::FillChannelMap(double eta, double phi, double energy) {
	AliEmcalTriggerPartChannel position = fTriggerMapping.GetPositionFromEtaPhi(eta, phi);
	if (position.IsEMCAL()) {
		if (!fBadChannelsEMCAL.HasChannel(position.GetCol(), position.GetRow()))
			fTriggerChannelsEMCAL.AddADC(position.GetCol(), position.GetRow(), energy);
	} else if (position.IsDCALPHOS()) {
		if (!fBadChannelsDCALPHOS.HasChannel(position.GetCol(), position.GetRow()))
			fTriggerChannelsDCALPHOS.AddADC(position.GetCol(), position.GetRow(), energy);
	}
}
