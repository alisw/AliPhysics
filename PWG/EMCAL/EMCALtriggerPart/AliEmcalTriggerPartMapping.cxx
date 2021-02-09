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
#include <cstdlib>

#include "AliEmcalTriggerPartMapping.h"

ClassImp(PWG::EMCAL::TriggerPart::AliEmcalTriggerPartMapping);
ClassImp(PWG::EMCAL::TriggerPart::AliEmcalTriggerPartChannel);

using namespace PWG::EMCAL::TriggerPart;

/**
 * Constructor, setting EMCAL and DCAL+PHOS dimensions
 *
 * Phi Limits:
 * EMCAL
 * Sector 0: min(1.40413), max(1.73746), dist(0.333328)
 * Sector 1: min(1.7532),  max(2.08653), dist(0.333328)
 * Sector 2: min(2.10226), max(2.43559), dist(0.333328)
 * Sector 3: min(2.45133), max(2.78466), dist(0.333328)
 * Sector 4: min(2.8004),  max(3.13372), dist(0.333328)
 * Sector 5: min(3.14946), max(3.26149), dist(0.112032)
 *
 * DCAL + PHOS
 * Sector 0: min(4.54573), max(4.87905), dist(0.333328)
 * Sector 1: min(4.89479), max(5.22812), dist(0.333328)
 * Sector 2: min(5.24386), max(5.57718), dist(0.333328)
 * Sector 3: min(5.59292), max(5.70495), dist(0.112032)
 */
AliEmcalTriggerPartMapping::AliEmcalTriggerPartMapping():
	TObject(),
	fPhiLimitsEMCAL(),
	fPhiLimitsDCALPHOS(),
	fEtaMin(-0.668305),
	fEtaMax(0.668305),
	fEtaSizeFOR()
{
	fEtaSizeFOR = (fEtaMax - fEtaMin) / 48;
	fPhiLimitsEMCAL.push_back(SectorPhi(0, 1.40413, 1.73746, 12));
	fPhiLimitsEMCAL.push_back(SectorPhi(1, 1.7532, 2.08653, 12));
	fPhiLimitsEMCAL.push_back(SectorPhi(2, 2.10226, 2.43559, 12));
	fPhiLimitsEMCAL.push_back(SectorPhi(3, 2.45133, 2.78466, 12));
	fPhiLimitsEMCAL.push_back(SectorPhi(4, 2.8004, 3.13372, 12));
	fPhiLimitsEMCAL.push_back(SectorPhi(5, 3.14946, 3.26149, 4));

	fPhiLimitsDCALPHOS.push_back(SectorPhi(0, 4.54573, 4.87905, 12));
	fPhiLimitsDCALPHOS.push_back(SectorPhi(1, 4.89479, 5.22812, 12));
	fPhiLimitsDCALPHOS.push_back(SectorPhi(2, 5.24386, 5.57718, 12));
	fPhiLimitsDCALPHOS.push_back(SectorPhi(3, 5.59292, 5.70495, 4));
}

/**
 * Destructor, nothing to do
 */
AliEmcalTriggerPartMapping::~AliEmcalTriggerPartMapping() {
}

/**
 * Map the position of the trigger channel. Always returns a trigger channel postition (also if outside the EMCAL
 * or DCAL+PHOS acceptance), which can however be undefined. In case of accessing row or col of an undefined trigger
 * position an exception is thrown.
 * @param eta Eta of the particle
 * @param phi Phi of the particle
 * @return The trigger channel corresponding to the Eta-Phi position of the particle.
 */
AliEmcalTriggerPartChannel AliEmcalTriggerPartMapping::GetPositionFromEtaPhi(double eta, double phi) const{
	if(IsEMCAL(eta, phi)){
		return GetPositionFromEtaPhiEMCAL(eta, phi);
	} else if(IsDCALPHOS(eta, phi)){
		return GetPositionFromEtaPhiDCALPHOS(eta, phi);
	}
	return AliEmcalTriggerPartChannel();
}

/**
 * Check whether particle is in the EMCAL trigger active area
 * @param eta Particle Eta
 * @param phi Particle Phi
 * @return True if the particle is in the active area of the EMCAL trigger, false otherwise
 */
bool AliEmcalTriggerPartMapping::IsEMCAL(double eta, double phi) const{
	if(eta < fEtaMin || eta > fEtaMax) return false;
	bool hasfound(false);
	for(std::vector<SectorPhi>::const_iterator phiit = fPhiLimitsEMCAL.begin(); phiit != fPhiLimitsEMCAL.end(); ++phiit){
		if(phiit->IsInSector(phi)){
			hasfound = true;
			break;
		}
	}
	return hasfound;
}

/**
 * Check whether particle is in the DCAL+PHOS trigger active area
 * @param eta Particle Eta
 * @param phi Particle Phi
 * @return True if the particle is in the active area of the EMCAL trigger, false otherwise
 */
bool AliEmcalTriggerPartMapping::IsDCALPHOS(double eta, double phi) const{
	if(eta < fEtaMin || eta > fEtaMax) return false;
	bool hasfound(false);
	for(std::vector<SectorPhi>::const_iterator phiit = fPhiLimitsDCALPHOS.begin(); phiit != fPhiLimitsDCALPHOS.end(); ++phiit){
		if(phiit->IsInSector(phi)){
			hasfound = true;
			break;
		}
	}
	return hasfound;
}

/**
 * Map eta and phi to trigger position in EMCAL. Mapping is done in the following way:
 * - linear model between min and max phi using a constant size in phi of the FastOR
 * @param eta Particle Eta
 * @param phi Particle Phi
 * @return Trigger position in the EMCAL
 */
AliEmcalTriggerPartChannel AliEmcalTriggerPartMapping::GetPositionFromEtaPhiEMCAL(double eta, double phi) const{
	AliEmcalTriggerPartChannel result;
	int row(-1), col(-1);
	// first get the row
	const SectorPhi *emcsec = FindSectorEMCAL(phi);
	if(!emcsec) 		// dead area
		return result;
	if(emcsec->GetSectorID() == 0){
		row = emcsec->GetRowNumberInSector(phi);
	} else {
		int rowsec = emcsec->GetRowNumberInSector(phi);
		row = 0;
		for(int isec = 0; isec < 6; isec++){
			if(isec < emcsec->GetSectorID()) row += fPhiLimitsEMCAL[isec].GetNumberOfRows();
			else break;
		}
		row += rowsec;
	}
	// then get the column
	// assume linear model, mapping from positive to negative eta as obtained from fastor + cell mapping from AliEMCALGeometry
	for(int coliter = 0; coliter < 48; coliter++){
		if(eta > fEtaMax - (coliter+1) * fEtaSizeFOR && eta < fEtaMax - coliter * fEtaSizeFOR){
			col = coliter;
			break;
		}
	}

	if(col >= 0)
		result.Set(row, col, AliEmcalTriggerPartChannel::kEMCAL);
	return result;
}

/**
 * Map eta and phi to trigger position in EMCAL. Mapping is done in the following way:
 * - linear model between min and max phi using a constant size in phi of the FastOR
 * @param eta Particle Eta
 * @param phi Particle Phi
 * @return Trigger position in the EMCAL
 */
AliEmcalTriggerPartChannel AliEmcalTriggerPartMapping::GetPositionFromEtaPhiDCALPHOS(double eta, double phi) const{
	AliEmcalTriggerPartChannel result;
	int row(-1), col(-1);
	// first get the row
	// first get the row
	const SectorPhi *emcsec = FindSectorDCALPHOS(phi);
	if(!emcsec) 		// dead area
		return result;
	if(emcsec->GetSectorID() == 0){
		row = emcsec->GetRowNumberInSector(phi);
	} else {
		int rowsec = emcsec->GetRowNumberInSector(phi);
		row = 0;
		for(int isec = 0; isec < 4; isec++){
			if(isec < emcsec->GetSectorID()) row += fPhiLimitsDCALPHOS[isec].GetNumberOfRows();
			else break;
		}
		row += rowsec;
	}
	// then get the column
	// assume linear model, mapping from positive to negative eta as obtained from fastor + cell mapping from AliEMCALGeometry
	for(int coliter = 0; coliter < 48; coliter++){
		if(eta > fEtaMax - (coliter+1) * fEtaSizeFOR && eta < fEtaMax - coliter * fEtaSizeFOR){
			col = coliter;
			break;
		}
	}

	if(col >= 0)
		result.Set(row, col, AliEmcalTriggerPartChannel::kDCALPHOS);
	return result;
}

/**
 * Find the sector according to the given phi angle, either in the EMCAL or in the DCAL/PHOS
 * @param phi Phi of the track/particle
 * @param isEMCAL Switch whether to test EMCAL or DCAL/PHOS
 * @return Full sector information (NULL if not found)
 */
const AliEmcalTriggerPartMapping::SectorPhi *AliEmcalTriggerPartMapping::FindSector(double phi, bool isEMCAL) const {
	if(isEMCAL) return FindSectorEMCAL(phi);
	else return FindSectorDCALPHOS(phi);
}

/**
 * Find the sector in the EMCAL for a given phi of the particle/track
 * @param phi Particle/Track phi
 * @return Sector information according to the phi (NULL if not found)
 */
const AliEmcalTriggerPartMapping::SectorPhi *AliEmcalTriggerPartMapping::FindSectorEMCAL(double phi) const {
	const SectorPhi *result = NULL;
	for(std::vector<SectorPhi>::const_iterator secit = fPhiLimitsEMCAL.begin(); secit != fPhiLimitsEMCAL.end(); ++secit){
		if(secit->IsInSector(phi)){
			result = &(*secit);
			break;
		}
	}
	return result;
}

/**
 * Find the sector in the DCAL/PHOS for a given phi of the particle/track
 * @param phi Particle/Track phi
 * @return Sector information according to the phi (NULL if not found)
 */
const AliEmcalTriggerPartMapping::SectorPhi *AliEmcalTriggerPartMapping::FindSectorDCALPHOS(double phi) const {
	const SectorPhi *result = NULL;
	for(std::vector<SectorPhi>::const_iterator secit = fPhiLimitsDCALPHOS.begin(); secit != fPhiLimitsDCALPHOS.end(); ++secit){
		if(secit->IsInSector(phi)){
			result = &(*secit);
			break;
		}
	}
	return result;
}

/**
 * Dummy constructor for ROOT I/O
 */
AliEmcalTriggerPartMapping::SectorPhi::SectorPhi():
	fSectorID(-1),
	fMinimum(-100),
	fMaximum(-100),
	fNRows(-1)
{}

/**
 * Initialize the Sector with mandatory information to map the position into a row number
 * @param sectorID ID of the sector, starting from 0
 * @param phiMin Min. phi of the sector
 * @param phiMax Max. phi of the sector
 * @param nrow Number of rows in the sector
 */
AliEmcalTriggerPartMapping::SectorPhi::SectorPhi(int sectorID, double phiMin, double phiMax, int nrow):
	fSectorID(sectorID),
	fMinimum(phiMin),
	fMaximum(phiMax),
	fNRows(nrow)
{}

/**
 * Calculate row number in sector from the phi of the particle / track
 * @param phi
 * @return Row number of the FastOR within the chamber (-1 if not found)
 */
int AliEmcalTriggerPartMapping::SectorPhi::GetRowNumberInSector(double phi) const {
	int rownumber = -1;
	double phiwidth = (fMaximum - fMinimum) / fNRows;
	int rowcounter = 0;
	for(double phiiter = fMinimum; phiiter < fMaximum; phiiter += phiwidth){
		if(phi > phiiter && phi < phiiter + phiwidth){
			rownumber = rowcounter;
			break;
		}
		rowcounter++;
	}
	return rownumber;
}
