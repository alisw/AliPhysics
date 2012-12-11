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
// AliCDMesonUtilsStripped
//
//  Author:
//  Felix Reidt <Felix.Reidt@cern.ch>

#ifndef ALICDMESONUTILSSTRIPPED_H
#define ALICDMESONUTILSSTRIPPED_H

class TH1;
class TH2;
class TLorentzVector;

class AliESDEvent;
class AliAODEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliVTrack;
class AliPIDResponse;

class AliCDMesonTracks;

class AliCDMesonUtilsStripped
{
public:
	enum{
		kInnerPixel = 0,
		kOuterPixel,
		kIPA,
		kIPC,
		kOPA,
		kOPC
	};

	// ESD only
	//---------

	// cuts for ESD analysis
	static Bool_t CutEvent(const AliESDEvent *ESDEvent);

	// gap determination
	static Int_t GetGapConfig(const AliESDEvent *ESDEvent);

	static void SPDLoadGeom(Int_t run); // only needed for ESDs, not in AODs

	// AOD only
	//---------
	static Bool_t CutEvent(const AliAODEvent *AODEvent);

private:
	// ESD only
	//---------

	// Gap determination functions
	static Int_t GetV0(const AliESDEvent *ESDEvent);
	static Int_t GetFMD(const AliESDEvent *ESDEvent);
	static Int_t GetSPD(const AliESDEvent *ESDEvent);
	static Int_t GetTPC(const AliESDEvent *ESDEvent);
	static Int_t GetZDC(const AliESDEvent *ESDEvent); // not used so far

	// helpers for the SPD gap determination
	static Bool_t SPDLoc2Glo(Int_t id, const Double_t *loc, Double_t *glo);
	static Int_t CheckChipEta(Int_t chipKey, TString scut,
	                          const Double_t vtxPos[]);
	static void GetNFO(const AliESDEvent *ESDEvent, TString etacut,
	                   Int_t ctr[]);
	// AOD only
	//---------


	// independent
	//----------
};

#endif
