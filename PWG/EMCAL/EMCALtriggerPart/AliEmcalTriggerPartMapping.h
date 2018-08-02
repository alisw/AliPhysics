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
#ifndef ALIEMCALTRIGGERPARTMAPPING_H
#define ALIEMCALTRIGGERPARTMAPPING_H

#include <exception>
#include <vector>
#include <TObject.h>

namespace PWG {

namespace EMCAL {

namespace TriggerPart {

class AliEmcalTriggerPartChannel : public TObject {
public:
	enum Detector{
		kEMCAL,
		kDCALPHOS,
		kUndefined,
	};
	class TriggerChannelException : public std::exception{
	public:
		TriggerChannelException() {}
		virtual ~TriggerChannelException() throw() {}
		const char *what() const throw() {
			return "Trigger channel not existing";
		}
	};

	AliEmcalTriggerPartChannel():
		TObject(),
		fDetector(kUndefined),
		fRow(0),
		fCol(0)
	{}
	virtual ~AliEmcalTriggerPartChannel() {}

	void Set(int row, int col, Detector det){
		fRow = row;
		fCol = col;
		fDetector = det;
	}

	int GetRow() const { if(fDetector == kUndefined) throw TriggerChannelException(); return fRow; }
	int GetCol() const { if(fDetector == kUndefined) throw TriggerChannelException(); return fCol; }

	bool IsEMCAL() const { if(fDetector == kEMCAL) return true; return false; }
	bool IsDCALPHOS() const { if(fDetector == kDCALPHOS) return true; return false; }

private:
	Detector		fDetector;
	int 			fRow;
	int 			fCol;

	ClassDef(AliEmcalTriggerPartChannel, 1);
};

class AliEmcalTriggerPartMapping : public TObject {
public:
	AliEmcalTriggerPartMapping();
	virtual ~AliEmcalTriggerPartMapping();

	AliEmcalTriggerPartChannel GetPositionFromEtaPhi(double eta, double phi) const;
	bool IsEMCAL(double eta, double phi) const;
	bool IsDCALPHOS(double eta, double phi) const;

	class SectorPhi{
	public:
		SectorPhi();
		SectorPhi(int sectorID, double phiMin, double phiMax, int nrow);
		~SectorPhi() {}

		int	 		GetSectorID() const { return fSectorID; }
		int	 		GetNumberOfRows() const { return fNRows; }
		double	 	GetPhiMin() const { return fMinimum; }
		double 		GetPhiMax() const { return fMaximum; }

		bool IsInSector(double phi)  const { return phi > fMinimum && phi < fMaximum; }
		int GetRowNumberInSector(double phi) const;

	private:
		int 			fSectorID;			///< ID of the sector, starting from 0 (Indices separate for EMCAL and DCAL)
		double 			fMinimum;			///< Min. phi of the sector (0 - 2 pi)
		double			fMaximum;			///< Max. phi of the sector (0 - 2 pi)
		int 			fNRows;  			///< Number of rows in a sector in phi (12 for big supermodules, 4 for small)
	};

protected:
	const SectorPhi *FindSector(double phi, bool isEMCAL) const;
	const SectorPhi *FindSectorEMCAL(double phi) const;
	const SectorPhi *FindSectorDCALPHOS(double phi) const;
	AliEmcalTriggerPartChannel GetPositionFromEtaPhiEMCAL(double eta, double phi) const;
	AliEmcalTriggerPartChannel GetPositionFromEtaPhiDCALPHOS(double eta, double phi) const;

	std::vector<SectorPhi>		fPhiLimitsEMCAL;
	std::vector<SectorPhi> 		fPhiLimitsDCALPHOS;
	double						fEtaMin;
	double						fEtaMax;
	double						fEtaSizeFOR;

	ClassDef(AliEmcalTriggerPartMapping, 1);
};

}
}
}

#endif /* ALIEMCALTRIGGERPARTMAPPING_H_*/
