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
// AliDDMesonBase
// for
// AliAnalysisTaskDDMeson
//
//  Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>

#include "TH2.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TAxis.h"
#include "THnSparse.h"
#include "TString.h"
#include "TList.h"

#include "AliPBBase.h"


//------------------------------------------------------------------------------
//== INVARIANT MASS DISTRIBUTIONS (ThnMother) ================================
//-- event characteristics
// number of charged primary particles - combined => can be w/ soft tracks
const Int_t AliPBBase::fgkNNcombined      = 2; // number of bins
const Double_t AliPBBase::fgkMinNcombined = 2; // lower border
const Double_t AliPBBase::fgkMaxNcombined = 4; // upper border (#bins + lower border)

// unlike sign or like sign charged?
const Int_t AliPBBase::fgkNCombCh      = 2; // kBinPPMM = number of bins
const Double_t AliPBBase::fgkMinCombCh = 1.; // kBinPM = lower border
const Double_t AliPBBase::fgkMaxCombCh = 3.; // kBinPPMM + kBinPM = upper border

// two track PID, take care the enum changed from time to time
const Int_t AliPBBase::fgkNCombPID      = 13; // kBinPIDUnknown = number of bins
const Double_t AliPBBase::fgkMinCombPID = 1.; // kBinPionE = lower border
const Double_t AliPBBase::fgkMaxCombPID = 14.;// ...  = upper border

// gap configuration, used for different detectors
const Int_t AliPBBase::fgkNGapConfig      = 4; // kBinNG = number of bins
const Double_t AliPBBase::fgkMinGapConfig = 1.; // kBinDG = lower border
const Double_t AliPBBase::fgkMaxGapConfig = 5.; // kBinNG + kBinDG = upper border

//-- mother kinematics
// invariant mass of the two-track system
const Int_t AliPBBase::fgkNMass      = 1024; // number of bins
const Double_t AliPBBase::fgkMinMass = 0.; // lower border
const Double_t AliPBBase::fgkMaxMass = 5.12; // upper border

// transverse momentum of the two-track system
const Int_t AliPBBase::fgkNMotherPt      = 128;//10; // number of bins
const Double_t AliPBBase::fgkMinMotherPt = 0.; // lower border
const Double_t AliPBBase::fgkMaxMotherPt = .64; //6.4 //2; // upper border

// cosine theta* (opening angle of the two daugther tracks in the
// centre-of-mass system of the two-track/mother system)
// **no meaning full variable**
const Int_t AliPBBase::fgkNCTS      = 2; // number of bins
const Double_t AliPBBase::fgkMinCTS = -1.; // lower border
const Double_t AliPBBase::fgkMaxCTS = -0.9; // upper border

// opening angle in the lab frame
const Int_t AliPBBase::fgkNOA      = 20; // number of bins
const Double_t AliPBBase::fgkMinOA = -1.; // lower border
const Double_t AliPBBase::fgkMaxOA = 1.; // upper border

//-- daughter kinematics
// transverse momentum of one of the two daughter particles
// (randomly selected)
const Int_t AliPBBase::fgkNDaughterPt      = 128; // number of bins
const Double_t AliPBBase::fgkMinDaughterPt = 0.; // lower border
const Double_t AliPBBase::fgkMaxDaughterPt = 6.4; // upper border

// pseudo rapidity of one of the two daughter particles
// (randomly selected)
//const Int_t AliPBBase::fgkNDaughterEta      = 64; // number of bins
//const Double_t AliPBBase::fgkMinDaughterEta = -1.28; // lower border
//const Double_t AliPBBase::fgkMaxDaughterEta =  1.28; // upper border

//-- Event quality information
// boolean values to reduce output size

// are there tracks in addition to the ones selected using AliPBTracks
// (ITSTPC, ITSsa, ITSpureSA) (0 = no, 1 = yes)
const Int_t AliPBBase::fgkNTrackResiduals      = 2; // number of bins
const Double_t AliPBBase::fgkMinTrackResiduals = 0.; // lower border
const Double_t AliPBBase::fgkMaxTrackResiduals = 2.; // upper border

// vertex with in +/-4cm (0 = no, 1 = yes)
const Int_t AliPBBase::fgkNVertexZinRng      = 2; // number of bins
const Double_t AliPBBase::fgkMinVertexZinRng = 0.; // lower border
const Double_t AliPBBase::fgkMaxVertexZinRng = 2.; // upper border

// are the vertices from SPD and tracks within 0.5cm? (0 = no, 1 = yes)
const Int_t AliPBBase::fgkNVertexCoincidence      = 2; // number of bins
const Double_t AliPBBase::fgkMinVertexCoincidence = 0.; // lower border
const Double_t AliPBBase::fgkMaxVertexCoincidence = 2.; // upper border

// are there SPD tracklets which are not assigned to tracks? (0 = no, 1 = yes)
const Int_t AliPBBase::fgkNTrackletResiduals      = 2; // number of bins
const Double_t AliPBBase::fgkMinTrackletResiduals = 0.; // lower border
const Double_t AliPBBase::fgkMaxTrackletResiduals = 2.; // upper border

//-- MC event information
const Int_t AliPBBase::fgkNProcessType      = 4; // kBinDD = number of bins
const Double_t AliPBBase::fgkMinProcessType = 0.; // kBinND = lower border
const Double_t AliPBBase::fgkMaxProcessType = 4.; // kBinDD = upper border


//== EMPTY EVENT STUDY (ThnEmptyEvents) ======================================
// event type
const Int_t AliPBBase::fgkNEventType      = 5; // kBinEventE = number of bins (5)
const Double_t AliPBBase::fgkMinEventType = 1.; // kBinEventI = lower border (1)
const Double_t AliPBBase::fgkMaxEventType = 6.; // kBinEventE+kBinEventI = u. b.

// multiplicities (reused for different detectors and ways of counting)
const Int_t AliPBBase::fgkNMult      = 32; // number of bins
const Double_t AliPBBase::fgkMinMult = 0.; // lower border
const Double_t AliPBBase::fgkMaxMult = 31.; // upper border

// multplicities - extended range
// (reused for different detectors and ways of counting)
const Int_t AliPBBase::fgkNMultW      = 64; // number of bins
const Double_t AliPBBase::fgkMinMultW = 0; // lower border
const Double_t AliPBBase::fgkMaxMultW = 63; // upper border

//== MULTIPLICITY STUDY (TnnMultiplicity) ====================================
// number of ITSTPC tracks in event
const Int_t AliPBBase::fgkNNch      = 51; // number of bins
const Double_t AliPBBase::fgkMinNch = 0.; // lower border
const Double_t AliPBBase::fgkMaxNch = 51.; // upper border

// number of ITS standalone tracks in event
const Int_t AliPBBase::fgkNNsoft      = 11; // number of bins
const Double_t AliPBBase::fgkMinNsoft = 0.; // lower border
const Double_t AliPBBase::fgkMaxNsoft = 11.; // upper border

// combined multiplicity
const Int_t AliPBBase::fgkNNcomb      = 61; // number of bins
const Double_t AliPBBase::fgkMinNcomb = 0.; // lower border
const Double_t AliPBBase::fgkMaxNcomb = 61.; // upper border

// gap configuration is reused from THnMother

// number of residual tracks
const Int_t AliPBBase::fgkNNresidualTracks      = 11; // number of bins
const Double_t AliPBBase::fgkMinNresidualTracks = 0.; // lower border
const Double_t AliPBBase::fgkMaxNresidualTracks = 11.; // upper border

// number of residual tracklets
const Int_t AliPBBase::fgkNNresidualTracklets      = 21; // number of bins
const Double_t AliPBBase::fgkMinNresidualTracklets = 0.; // lower border
const Double_t AliPBBase::fgkMaxNresidualTracklets = 21.; // upper border

// vertex z-position
const Int_t AliPBBase::fgkNVertexZ      = 20; // number of bins
const Double_t AliPBBase::fgkMinVertexZ = -10.; // lower border
const Double_t AliPBBase::fgkMaxVertexZ = 10.; // upper border

// SPD and track vertex distance
const Int_t AliPBBase::fgkNVerticesDistance      = 10; // number of bins
const Double_t AliPBBase::fgkMinVerticesDistance = 0.; // lower border
const Double_t AliPBBase::fgkMaxVerticesDistance = 5.; // upper border


//------------------------------------------------------------------------------
Int_t AliPBBase::GetGapBin(TString tag, Int_t gapcg,
                                Bool_t checkCentralActivity /* = kTRUE */)
{
	//
	// retrieve gap topology for a given string
	//

	tag.ToUpper();

	Bool_t ka = kFALSE, kc = kFALSE;

	if(tag.Contains("V0")){
		ka = ka || (gapcg & kBitV0A);
		kc = kc || (gapcg & kBitV0C);
	}
	if(tag.Contains("FMD")){
		ka = ka || (gapcg & kBitFMDA);
		kc = kc || (gapcg & kBitFMDC);
	}
	if(tag.Contains("AD")){
		ka = ka || (gapcg & kBitADA);
		kc = kc || (gapcg & kBitADC);
	}
	if(tag.Contains("SPD")){
		ka = ka || (gapcg & kBitSPDA);
		kc = kc || (gapcg & kBitSPDC);
	}
	if(tag.Contains("TPC")){
		ka = ka || (gapcg & kBitTPCA);
		kc = kc || (gapcg & kBitTPCC);
	}
	if(tag.Contains("ZDC")){
		ka = ka || (gapcg & kBitZDCA);
		kc = kc || (gapcg & kBitZDCC);
	}

	if(ka && kc)
		return kBinNG;
	else{
		if(!ka && !kc)
			if (((gapcg & kBitCentAct) && checkCentralActivity) ||
			    !checkCentralActivity) {
				return kBinDG; // central activity seen (or not required)
			}
			else {
				return kBinNG; // no central activity
			}
		else if(!kc)
			return kBinGC;
		else
			return kBinGA;
	}
}


//------------------------------------------------------------------------------
void AliPBBase::CheckRange(Double_t &var, Double_t min,
                                Double_t max)
{
	//
	// check whether the value is with in the range specified with min and max
	//

	const Double_t eps = 1e-3;
	if( var >= max ) var = max - eps;
	if( var <= min ) var = min + eps;
}


//------------------------------------------------------------------------------
Int_t AliPBBase::GetAxis(TString thntit, TString name)
{
	//
	// get the number of an axis, derived from the ThnSparse title (list of axes)
	//

	thntit.ToUpper();
	thntit.ReplaceAll(" ","");
	const Int_t nmax = 20;
	TString tits[nmax];
	Int_t counter = 0;
	while(thntit.Contains(",")){
		const Int_t pos = thntit.First(",");
		tits[counter] = thntit(0, pos);
		thntit = thntit(pos+1, thntit.Length()-pos);
		counter++;
		if(counter>=nmax-1){
			printf("AliPBBase AliPBBase::GetAxis too small nmax! %d %d\n",
			       counter, nmax);
			return -1; //exit(1); // TODO
		}
	}
	tits[counter++] = thntit;

	//----------------

	name.ToUpper();

	for(Int_t ii=0; ii<counter; ii++){
	  if( tits[ii] == name )
		  return ii;
	}
	printf("AliPBBase AliPBBase::GetAxis !%s! %s not found!\n",
	       name.Data(), thntit.Data());
	for(Int_t ii=0; ii<counter; ii++){
		printf("*************** AliPBBase::GetAxis *****************\n");
		printf("AliPBBase AliPBBase::GetAxis %d !%s!\n", ii,
		       tits[ii].Data());
		printf("\n");
	}
	return -1; //exit(1); // TODO
}


//------------------------------------------------------------------------------
TString AliPBBase::GetTitleMother()
{
	//
	// get the title of ThnMother (= list of axes)
	//

	TString title = "Ncombined, CombCh, CombPID, V0, FMD, SPD, TPC,";
	title += " Mass, Pt, CTS, OA, DaughterPt, TrackResiduals, VertexZinRng,";
	title += " ProcessType, VertexCoincidence, TrackletResiduals";
	return title;
}


//------------------------------------------------------------------------------
THnSparseD * AliPBBase::GetThnMother(TString name /* = "PB_Mother" */)
{
	// creates the THnSparse called Mother (used for the generation of invariant
	// mass distributions)
	//

	const Int_t nbin[] = {
		fgkNNcombined, fgkNCombCh, fgkNCombPID, fgkNGapConfig, fgkNGapConfig,
		fgkNGapConfig, fgkNGapConfig, fgkNMass, fgkNMotherPt, fgkNCTS, fgkNOA,
		fgkNDaughterPt, fgkNTrackResiduals, fgkNVertexZinRng, fgkNProcessType,
		fgkNVertexCoincidence, fgkNTrackletResiduals
	};
	const Double_t binmin[] = {
		fgkMinNcombined, fgkMinCombCh, fgkMinCombPID, fgkMinGapConfig,
		fgkMinGapConfig, fgkMinGapConfig, fgkMinGapConfig, fgkMinMass,
		fgkMinMotherPt, fgkMinCTS, fgkMinOA, fgkMinDaughterPt,
		fgkMinTrackResiduals, fgkMinVertexZinRng, fgkMinProcessType,
		fgkMinVertexCoincidence, fgkMinTrackletResiduals
	};
	const Double_t binmax[] = {
		fgkMaxNcombined, fgkMaxCombCh, fgkMaxCombPID, fgkMaxGapConfig,
		fgkMaxGapConfig, fgkMaxGapConfig, fgkMaxGapConfig, fgkMaxMass,
		fgkMaxMotherPt, fgkMaxCTS, fgkMaxOA, fgkMaxDaughterPt,
		fgkMaxTrackResiduals, fgkMaxVertexZinRng, fgkMaxProcessType,
		fgkMaxVertexCoincidence, fgkMaxTrackletResiduals
	};
	
	const Int_t npar = sizeof(nbin)/sizeof(Int_t);

	return new THnSparseD(name.Data(), GetTitleMother().Data(), npar, nbin, binmin,
	                      binmax);
}


//------------------------------------------------------------------------------
void AliPBBase::FillThnMother(THnSparseD *thn, Double_t vNch,
                                   Double_t vCombCh,
                                   Double_t vCombPID, Double_t vV0,
                                   Double_t vFMD, Double_t vSPD,
                                   Double_t vTPC, Double_t vMass,
                                   Double_t vPt,  Double_t vOA,
                                   Double_t vCTS,
                                   Double_t vDaughterPt,
                                   Double_t vTrackResiduals,
                                   Double_t vVertexZ,
                                   Double_t vProcessType,
                                   Double_t vVertexCoincidence,
                                   Double_t vTrkltResiduals)
{
	//
	// fill ThnMother
	//

	Double_t var[]={
		vNch, vCombCh, vCombPID, vV0, vFMD, vSPD, vTPC, vMass, vPt, vOA, vCTS,
		vDaughterPt, vTrackResiduals, vVertexZ, vProcessType, vVertexCoincidence,
		vTrkltResiduals
	};
	const Int_t nv = sizeof(var)/sizeof(Double_t);
	if(nv!=thn->GetNdimensions()){
		printf("AliPBBase::FillThnMother nv error!! %d %d\n", nv,
		       thn->GetNdimensions());
		return; //exit(1);
	}

	CheckRange(var[7], fgkMinMass, fgkMaxMass);
	CheckRange(var[8], fgkMinMotherPt, fgkMaxMotherPt);
	CheckRange(var[11], fgkMinDaughterPt, fgkMaxDaughterPt);

	thn->Fill(var);
}


//------------------------------------------------------------------------------
Int_t AliPBBase::GetAxisMother(TString name)
{
	//
	// return axis number corresponding to the name
	//

	return GetAxis(GetTitleMother(), name);
}


//------------------------------------------------------------------------------
TString AliPBBase::GetTitleEmptyEvents()
{
	//
	// title / list of axes of the empty event study THnSparse
	//

	TString title = "EventType, FMD-A, FMD-C, SPD-I-A, SPD-I-C, SPD-O-A, SPD-O-C";
	title += ", SPDtrkltA, SPDtrkltC";
	title += ", fmdSum1I, fmdSum2I, fmdSum2O, fmdSum3I, fmdSum3O";
	//title += ", TPC_all, TPC_diffVertex";
	return title;
}


//------------------------------------------------------------------------------
THnSparseD* AliPBBase::GetThnEmptyEvents()
{
	// creates the THnSparse called Empty Events
	// EventType, FMD-A, FMD-C, SPD-I-A, SPD-I-C, SPD-O-A, SPD-O-C, FMD1I..FMD3O
	// TPC_all, TPC_diffVertex (not used so far)

	const Int_t nbin[] = {
		fgkNEventType, fgkNMultW, fgkNMultW, fgkNMult, fgkNMult, fgkNMult, fgkNMult,
		fgkNMult, fgkNMult, fgkNMult, fgkNMult, fgkNMult, fgkNMult, fgkNMult
		//, fgkNMult, fgkNMult
	};
	const Double_t binmin[] = {
		fgkMinEventType, fgkMinMultW, fgkMinMultW, fgkMinMult, fgkMinMult,
		fgkMinMult, fgkMinMult, fgkMinMult, fgkMinMult, fgkMinMult, fgkMinMult,
		fgkMinMult, fgkMinMult, fgkMinMult //, fgkMinMult, fgkMinMult
	};
	const Double_t binmax[] = {
		fgkMaxEventType, fgkMaxMultW, fgkMaxMultW, fgkMaxMult, fgkMaxMult,
		fgkMaxMult, fgkMaxMult, fgkMaxMult, fgkMaxMult, fgkMaxMult, fgkMaxMult,
		fgkMaxMult, fgkMaxMult, fgkMaxMult //, fgkMaxMult, fgkMaxMult
	};
	
	const Int_t npar = sizeof(nbin)/sizeof(Int_t);

	return new THnSparseD("PB_EmptyEvents", GetTitleEmptyEvents(), npar,
	                      nbin, binmin, binmax);
}


//------------------------------------------------------------------------------
void AliPBBase::FillThnEmptyEvents(THnSparseD * thn, Int_t eventType,
                                        Int_t multFMDA,
                                        Int_t multFMDC,
                                        Int_t multSPDIA,
                                        Int_t multSPDIC,
                                        Int_t multSPDOA,
                                        Int_t multSPDOC,
                                        Int_t multSPDtrkltA,
                                        Int_t multSPDtrkltC,
                                        Int_t fmdSum1I,
                                        Int_t fmdSum2I,
                                        Int_t fmdSum2O,
                                        Int_t fmdSum3I,
                                        Int_t fmdSum3O/*,
                                        Int_t multTPC,
                                        Int_t multTPCdiffVertex */)
{
	//
	// Fill ThnEmptyEvents
	//

	Double_t var[]={
	  static_cast<Double_t>(eventType), static_cast<Double_t>(multFMDA), static_cast<Double_t>(multFMDC), static_cast<Double_t>(multSPDIA), static_cast<Double_t>(multSPDIC), static_cast<Double_t>(multSPDOA),
	  static_cast<Double_t>(multSPDOC), static_cast<Double_t>(multSPDtrkltA), static_cast<Double_t>(multSPDtrkltC), static_cast<Double_t>(fmdSum1I), static_cast<Double_t>(fmdSum2I), static_cast<Double_t>(fmdSum2O),
	  static_cast<Double_t>(fmdSum3I), static_cast<Double_t>(fmdSum3O)//, multTPC, multTPCdiffVertex
	};
	const Int_t nv = sizeof(var)/sizeof(Double_t);
	if(nv!=thn->GetNdimensions()){
		printf("AliPBBase::FillThnEmptyEvents nv error!! %d %d\n", nv,
		       thn->GetNdimensions());
		return; //exit(1); // TODO
	}

	thn->Fill(var);
}


//------------------------------------------------------------------------------
Int_t AliPBBase::GetAxisEmptyEvents(TString name)
{
	//
	// return axis number corresponding to the name
	//

	return GetAxis(GetTitleEmptyEvents(), name);
}


//------------------------------------------------------------------------------
TString AliPBBase::GetTitleMultiplicity()
{
	//
	// get title of the multiplicity study ThnSparse
	//

	TString title = "Nch, Nsoft, Ncombined, V0, FMD, SPD, TPC, NresidualTracks";
	title += ", NresidualTracklets, VertexZ, VerticesDistance, ProcessType";
	return title;
}


//------------------------------------------------------------------------------
THnSparseD* AliPBBase::GetThnMultiplicity()
{
	//
	// creates the THnSparse called Multiplicity
	//

	const Int_t nbin[] = {
		fgkNNch, fgkNNsoft, fgkNNcomb,fgkNGapConfig, fgkNGapConfig,
		fgkNGapConfig, fgkNGapConfig, fgkNNresidualTracks, fgkNNresidualTracklets,
		fgkNVertexZ, fgkNVerticesDistance, fgkNProcessType
	};
	const Double_t binmin[] = {
		fgkMinNch, fgkMinNsoft, fgkMinNcomb, fgkMinGapConfig, fgkMinGapConfig,
		fgkMinGapConfig, fgkMinGapConfig, fgkMinNresidualTracks,
		fgkMinNresidualTracklets, fgkMinVertexZ, fgkMinVerticesDistance,
		fgkMinProcessType
	};
	const Double_t binmax[] = {
		fgkMaxNch, fgkMaxNsoft, fgkMaxNcomb, fgkMaxGapConfig, fgkMaxGapConfig,
		fgkMaxGapConfig, fgkMaxGapConfig, fgkMaxNresidualTracks,
		fgkMaxNresidualTracklets, fgkMaxVertexZ, fgkMaxVerticesDistance,
		fgkMaxProcessType
	};

	const Int_t npar = sizeof(nbin)/sizeof(Int_t);

	return new THnSparseD("PB_Multiplicity", GetTitleMultiplicity().Data(), npar,
	                      nbin, binmin, binmax);
}


//------------------------------------------------------------------------------
void AliPBBase::FillThnMultiplicity(THnSparseD *thn, Double_t vNch,
                                         Double_t vNsoft,
                                         Double_t vNcombined,
                                         Double_t vV0,
                                         Double_t vFMD,
                                         Double_t vSPD,
                                         Double_t vTPC,
                                         Double_t vNresidualTracks,
                                         Double_t vNresidualTracklets,
                                         Double_t vVertexZ,
                                         Double_t vVerticesDistance,
                                         Double_t vProcessType)
{
	// fill ThnMultiplicity
	// input list copied from GetTitle
	// var[] copied from GetTitle

	Double_t var[]={
		vNch, vNsoft, vNcombined, vV0, vFMD, vSPD, vTPC, vNresidualTracks,
		vNresidualTracklets, vVertexZ, vVerticesDistance, vProcessType
	};
	const Int_t nv = sizeof(var)/sizeof(Double_t);
	if(nv!=thn->GetNdimensions()){
		printf("AliPBBase::FillThnMultiplicity nv error!! %d %d\n", nv,
		       thn->GetNdimensions());
		return; //exit(1);
	}

	thn->Fill(var);
}


//------------------------------------------------------------------------------
Int_t AliPBBase::GetAxisMultiplicity(TString name)
{
	//
	// return axis number corresponding to the name
	//

	return GetAxis(GetTitleMultiplicity(), name);
}


//------------------------------------------------------------------------------
TH1F* AliPBBase::GetHistStatsFlow()
{
	//
	// setup the stats flow histogram
	//

	TH1F *hist = new TH1F("c00_statsFlow", "", AliPBBase::kBinLastValue,
	                0, AliPBBase::kBinLastValue);
	TAxis* axis = hist->GetXaxis();
	axis->SetBinLabel(AliPBBase::kBinTotalInput+1, "total Input");
	axis->SetBinLabel(AliPBBase::kBinGoodInput+1, "good ESDs");
	axis->SetBinLabel(AliPBBase::kBinV0OR+1, "V0-OR");
	axis->SetBinLabel(AliPBBase::kBinV0AND+1, "V0-AND");
	axis->SetBinLabel(AliPBBase::kBinEventsAfterCuts+1, "after cuts");
	axis->SetBinLabel(AliPBBase::kBinEventsWithOutPileUp+1, "w/o pile up");
	axis->SetBinLabel(AliPBBase::kBinv0Gap+1, "with V0 DG gap");
	axis->SetBinLabel(AliPBBase::kBinv0fmdGap+1, "with V0-FMD DG gap");
	axis->SetBinLabel(AliPBBase::kBinv0fmdspdGap+1,
	                  "with V0-FMD-SPD DG gap");
	axis->SetBinLabel(AliPBBase::kBinv0fmdspdtpcGap+1,
	                  "with V0-FMD-SPD-TPC DG gap");
	axis->SetBinLabel(AliPBBase::kBinv0fmdspdtpczdcGap+1,
	                  "with V0-FMD-SPD-TPC-ZDC DG gap");
	axis->SetBinLabel(AliPBBase::kBinfmdGap+1, "with FMD DG gap");
	axis->SetBinLabel(AliPBBase::kBinspdGap+1, "with SPD DG gap");
	axis->SetBinLabel(AliPBBase::kBintpcGap+1, "with TPC DG gap");
	axis->SetBinLabel(AliPBBase::kBintpcspdGap+1, "width TPC-SPD DG gap");
	axis->SetBinLabel(AliPBBase::kBintpcspdfmdGap+1,
	                  "width TPC-SPD-FMD DG gap");
	axis->SetBinLabel(AliPBBase::kBintpcspdfmdv0Gap+1,
	                  "width TPC-SPD-FMD-V0 DG gap");
	axis->SetBinLabel(AliPBBase::kBinspdfmdGap+1, "with SPD FMD gap");
	axis->SetBinLabel(AliPBBase::kBinspdfmdv0Gap+1, "with SPD FMD V0 gap");
	axis->SetBinLabel(AliPBBase::kBinTwoTrackEvents+1, "with two tracks");
	axis->SetBinLabel(AliPBBase::kBinThreeTrackEvents+1,
	                  "with three tracks");
	axis->SetBinLabel(AliPBBase::kBinPionEvents+1, "with two pions");
	axis->SetBinLabel(AliPBBase::kBinKaonEvents+1, "with two kaons");
	axis->SetBinLabel(AliPBBase::kBinProtonEvents+1, "with two proton");
	axis->SetBinLabel(AliPBBase::kBinElectronEvents+1, "with two electron");
	axis->SetBinLabel(AliPBBase::kBinUnknownPIDEvents+1, "with unknown PID");
	axis->SetBinLabel(AliPBBase::kBinResidualTracks+1,
	                  "without residual tracks");
	axis->SetBinLabel(AliPBBase::kBinResidualTracklets+1,
	                  "without residual tracklets");
	axis->SetBinLabel(AliPBBase::kBinCDonlyEvents+1, "CD only events");
	return hist;
}


//------------------------------------------------------------------------------
TH2F* AliPBBase::GetHistPIDStudies(TString name)
{
	//
	// setup the PID studies histogram
	//

	TH2F *hist = new TH2F(name.Data(), ";particle 1;particle 2",
	                      AliPBBase::kBinPIDUnknown,
	                      AliPBBase::kBinPionE,
	                      AliPBBase::kBinPIDUnknown+1,
	                      AliPBBase::kBinPIDUnknown,
	                      AliPBBase::kBinPionE,
	                      AliPBBase::kBinPIDUnknown+1);
	TAxis* x = hist->GetXaxis();
	TAxis* y = hist->GetYaxis();
	x->SetBinLabel(AliPBBase::kBinPionE, "#pi (ex)");
	x->SetBinLabel(AliPBBase::kBinPion, "#pi");
	x->SetBinLabel(AliPBBase::kBinSinglePion, "-");
	x->SetBinLabel(AliPBBase::kBinKaonE, "K (ex)");
	x->SetBinLabel(AliPBBase::kBinKaon, "K");
	x->SetBinLabel(AliPBBase::kBinSingleKaon, ",");
	x->SetBinLabel(AliPBBase::kBinProtonE, "p (ex)");
	x->SetBinLabel(AliPBBase::kBinProton, "p");
	x->SetBinLabel(AliPBBase::kBinSingleProton, "_");
	x->SetBinLabel(AliPBBase::kBinElectronE, "e (ex)");
	x->SetBinLabel(AliPBBase::kBinElectron, "e");
	x->SetBinLabel(AliPBBase::kBinSingleElectron, ".");
	x->SetBinLabel(AliPBBase::kBinPIDUnknown, "X");
	y->SetBinLabel(AliPBBase::kBinPionE, "#pi (ex)");
	y->SetBinLabel(AliPBBase::kBinPion, "#pi");
	y->SetBinLabel(AliPBBase::kBinSinglePion, "-");
	y->SetBinLabel(AliPBBase::kBinKaonE, "K (ex)");
	y->SetBinLabel(AliPBBase::kBinKaon, "K");
	y->SetBinLabel(AliPBBase::kBinSingleKaon, ",");
	y->SetBinLabel(AliPBBase::kBinProtonE, "p (ex)");
	y->SetBinLabel(AliPBBase::kBinProton, "p");
	y->SetBinLabel(AliPBBase::kBinSingleProton, "_");
	y->SetBinLabel(AliPBBase::kBinElectronE, "e (ex)");
	y->SetBinLabel(AliPBBase::kBinElectron, "e");
	y->SetBinLabel(AliPBBase::kBinSingleElectron, ".");
	y->SetBinLabel(AliPBBase::kBinPIDUnknown, "X");
	return hist;
}


//------------------------------------------------------------------------------
TObjArray* AliPBBase::GetHistVZEROStudies(TList* l)
{
	//
	// Create Histograms for the VZERO trigger studies
	//

	TObjArray* arr = new TObjArray(130);

	TObject* o = 0x0;
	for (Int_t iPMT = 0; iPMT < 64; ++iPMT) {
		o = (TObject*)new TH2F(Form("h00_%02d_ADC_TriggerThr", iPMT),
		                       ";ADC Counts;Trigger Threshold (ADC Counts)",
		                       400, 0., 4000., 48, 2., 50.);
		arr->Add(o);
		l->Add(o);
	}
	for (Int_t iPMT = 0; iPMT < 64; ++iPMT) {
		o = (TObject*)new TH2F(Form("h01_%02d_ADC_Multiplicity", iPMT),
		                       ";ADC Counts;Multiplicity",
		                       400, 0., 4000., 250, 0., 250.);
		arr->Add(o);
		l->Add(o);
	}
	/*
	  // not of use for pp data in 2010
	o = (TObject*)new TH2F("h02_TriggerChargeA_Trigger",
	                       ";Trigger Charge A;Trigger Decision A",
	                       250, 0., 5000., 2, 0., 2.);
	arr->Add(o);
	l->Add(o);
	o = (TObject*)new TH2F("h03__TriggerChargeC_Trigger",
	                       ";Trigger Charge C;Trigger Decision C",
	                       250, 0., 5000., 2, 0., 2.);
	arr->Add(o);
	l->Add(o);
	*/

	return arr;
}


//------------------------------------------------------------------------------
void AliPBBase::GetGapTriggers(THnSparseI* gaprun, Int_t gapCondition,
                                    Int_t run, Double_t& triggers,
                                    Double_t& total)
{
	// determine the number of certain triggers with a gap in the detectors
	// specified by gapCondition in run and the total number of events
	// surviving all cuts (including pile-up rejection)

	triggers = 0;
	total = 0;
	Int_t nTuple[] = { gaprun->GetAxis(0)->FindBin(run), 0 };
	for (Int_t i = 0; i< kBitGapMax; ++i) {
		nTuple[1] = i;
		Double_t temp = gaprun->GetBinContent(nTuple);
		if (!(i & gapCondition)) {
			triggers += temp;
		}
		total += temp;
	}
}


//------------------------------------------------------------------------------
void AliPBBase::GetNoGapTriggers(THnSparseI* gaprun, Int_t gapCondition,
                                      Int_t run, Double_t& triggers,
                                      Double_t& total)
{
	// determine the number of certain triggers with a NO gap in the detectors
	// specified by the gapCondition in run and the total number of events
	// surviving all cuts (including pile-up rejection)

	triggers = 0;
	total = 0;
	Int_t nTuple[] = { gaprun->GetAxis(0)->FindBin(run), 0 };
	for (Int_t i = 0; i< kBitGapMax; ++i) {
		nTuple[1] = i;
		Double_t temp = gaprun->GetBinContent(nTuple);
		if (i & gapCondition) {
			triggers += temp;
		}
		total += temp;
	}
}
