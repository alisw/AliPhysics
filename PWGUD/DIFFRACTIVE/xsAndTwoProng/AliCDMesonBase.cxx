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

#include "AliCDMesonBase.h"


//------------------------------------------------------------------------------
Int_t AliCDMesonBase::GetGapBin(TString tag, const Int_t gapcg,
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
void AliCDMesonBase::CheckRange(Double_t &var, const Double_t min,
                                const Double_t max)
{
	//
	// check whether the value is with in the range specified with min and max
	//

	const Double_t eps = 1e-3;
	if( var >= max ) var = max - eps;
	if( var <= min ) var = min + eps;
}


//------------------------------------------------------------------------------
Int_t AliCDMesonBase::GetAxis(TString thntit, TString name)
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
			printf("AliCDMesonBase AliCDMesonBase::GetAxis too small nmax! %d %d\n",
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
	printf("AliCDMesonBase AliCDMesonBase::GetAxis !%s! %s not found!\n",
	       name.Data(), thntit.Data());
	for(Int_t ii=0; ii<counter; ii++){
		printf("*************** AliCDMesonBase::GetAxis *****************\n");
		printf("AliCDMesonBase AliCDMesonBase::GetAxis %d !%s!\n", ii,
		       tits[ii].Data());
		printf("\n");
	}
	return -1; //exit(1); // TODO
}


//------------------------------------------------------------------------------
TString AliCDMesonBase::GetTitleMother()
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
THnSparseD * AliCDMesonBase::GetThnMother(TString name /* = "CDMeson_Mother" */)
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

	return new THnSparseD(name.Data(), GetTitleMother(), npar, nbin, binmin,
	                      binmax);
}


//------------------------------------------------------------------------------
void AliCDMesonBase::FillThnMother(THnSparseD *thn, const Double_t vNch,
                                   const Double_t vCombCh,
                                   const Double_t vCombPID, const Double_t vV0,
                                   const Double_t vFMD, const Double_t vSPD,
                                   const Double_t vTPC, const Double_t vMass,
                                   const Double_t vPt,  const Double_t vOA,
                                   const Double_t vCTS,
                                   const Double_t vDaughterPt,
                                   const Double_t vTrackResiduals,
                                   const Double_t vVertexZ,
                                   const Double_t vProcessType,
                                   const Double_t vVertexCoincidence,
                                   const Double_t vTrkltResiduals)
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
		printf("AliCDMesonBase::FillThnMother nv error!! %d %d\n", nv,
		       thn->GetNdimensions());
		return; //exit(1);
	}

	CheckRange(var[7], fgkMinMass, fgkMaxMass);
	CheckRange(var[8], fgkMinMotherPt, fgkMaxMotherPt);
	CheckRange(var[11], fgkMinDaughterPt, fgkMaxDaughterPt);

	thn->Fill(var);
}


//------------------------------------------------------------------------------
Int_t AliCDMesonBase::GetAxisMother(const TString name)
{
	//
	// return axis number corresponding to the name
	//

	return GetAxis(GetTitleMother(), name);
}


//------------------------------------------------------------------------------
TString AliCDMesonBase::GetTitleEmptyEvents()
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
THnSparseD* AliCDMesonBase::GetThnEmptyEvents()
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

	return new THnSparseD("CDMeson_EmptyEvents", GetTitleEmptyEvents(), npar,
	                      nbin, binmin, binmax);
}


//------------------------------------------------------------------------------
void AliCDMesonBase::FillThnEmptyEvents(THnSparseD * thn, const Int_t eventType,
                                        const Int_t multFMDA,
                                        const Int_t multFMDC,
                                        const Int_t multSPDIA,
                                        const Int_t multSPDIC,
                                        const Int_t multSPDOA,
                                        const Int_t multSPDOC,
                                        const Int_t multSPDtrkltA,
                                        const Int_t multSPDtrkltC,
                                        const Int_t fmdSum1I,
                                        const Int_t fmdSum2I,
                                        const Int_t fmdSum2O,
                                        const Int_t fmdSum3I,
                                        const Int_t fmdSum3O/*,
                                        const Int_t multTPC,
                                        const Int_t multTPCdiffVertex */)
{
	//
	// Fill ThnEmptyEvents
	//

	Double_t var[]={
		eventType, multFMDA, multFMDC, multSPDIA, multSPDIC, multSPDOA,
		multSPDOC, multSPDtrkltA, multSPDtrkltC, fmdSum1I, fmdSum2I, fmdSum2O,
		fmdSum3I, fmdSum3O//, multTPC, multTPCdiffVertex
	};
	const Int_t nv = sizeof(var)/sizeof(Double_t);
	if(nv!=thn->GetNdimensions()){
		printf("AliCDMesonBase::FillThnEmptyEvents nv error!! %d %d\n", nv,
		       thn->GetNdimensions());
		return; //exit(1); // TODO
	}

	thn->Fill(var);
}


//------------------------------------------------------------------------------
Int_t AliCDMesonBase::GetAxisEmptyEvents(const TString name)
{
	//
	// return axis number corresponding to the name
	//

	return GetAxis(GetTitleEmptyEvents(), name);
}


//------------------------------------------------------------------------------
TString AliCDMesonBase::GetTitleMultiplicity()
{
	//
	// get title of the multiplicity study ThnSparse
	//

	TString title = "Nch, Nsoft, Ncombined, V0, FMD, SPD, TPC, NresidualTracks";
	title += ", NresidualTracklets, VertexZ, VerticesDistance, ProcessType";
	return title;
}


//------------------------------------------------------------------------------
THnSparseD* AliCDMesonBase::GetThnMultiplicity()
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

	return new THnSparseD("CDMeson_Multiplicity", GetTitleMultiplicity(), npar,
	                      nbin, binmin, binmax);
}


//------------------------------------------------------------------------------
void AliCDMesonBase::FillThnMultiplicity(THnSparseD *thn, const Double_t vNch,
                                         const Double_t vNsoft,
                                         const Double_t vNcombined,
                                         const Double_t vV0,
                                         const Double_t vFMD,
                                         const Double_t vSPD,
                                         const Double_t vTPC,
                                         const Double_t vNresidualTracks,
                                         const Double_t vNresidualTracklets,
                                         const Double_t vVertexZ,
                                         const Double_t vVerticesDistance,
                                         const Double_t vProcessType)
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
		printf("AliCDMesonBase::FillThnMultiplicity nv error!! %d %d\n", nv,
		       thn->GetNdimensions());
		return; //exit(1);
	}

	thn->Fill(var);
}


//------------------------------------------------------------------------------
Int_t AliCDMesonBase::GetAxisMultiplicity(const TString name)
{
	//
	// return axis number corresponding to the name
	//

	return GetAxis(GetTitleMultiplicity(), name);
}


//------------------------------------------------------------------------------
TH1F* AliCDMesonBase::GetHistStatsFlow()
{
	//
	// setup the stats flow histogram
	//

	TH1F *hist = new TH1F("c00_statsFlow", "", AliCDMesonBase::kBinLastValue,
	                0, AliCDMesonBase::kBinLastValue);
	TAxis* axis = hist->GetXaxis();
	axis->SetBinLabel(AliCDMesonBase::kBinTotalInput+1, "total Input");
	axis->SetBinLabel(AliCDMesonBase::kBinGoodInput+1, "good ESDs");
	axis->SetBinLabel(AliCDMesonBase::kBinV0OR+1, "V0-OR");
	axis->SetBinLabel(AliCDMesonBase::kBinV0AND+1, "V0-AND");
	axis->SetBinLabel(AliCDMesonBase::kBinEventsAfterCuts+1, "after cuts");
	axis->SetBinLabel(AliCDMesonBase::kBinEventsWithOutPileUp+1, "w/o pile up");
	axis->SetBinLabel(AliCDMesonBase::kBinv0Gap+1, "with V0 DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBinv0fmdGap+1, "with V0-FMD DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBinv0fmdspdGap+1,
	                  "with V0-FMD-SPD DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBinv0fmdspdtpcGap+1,
	                  "with V0-FMD-SPD-TPC DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBinv0fmdspdtpczdcGap+1,
	                  "with V0-FMD-SPD-TPC-ZDC DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBinfmdGap+1, "with FMD DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBinspdGap+1, "with SPD DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBintpcGap+1, "with TPC DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBintpcspdGap+1, "width TPC-SPD DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBintpcspdfmdGap+1,
	                  "width TPC-SPD-FMD DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBintpcspdfmdv0Gap+1,
	                  "width TPC-SPD-FMD-V0 DG gap");
	axis->SetBinLabel(AliCDMesonBase::kBinspdfmdGap+1, "with SPD FMD gap");
	axis->SetBinLabel(AliCDMesonBase::kBinspdfmdv0Gap+1, "with SPD FMD V0 gap");
	axis->SetBinLabel(AliCDMesonBase::kBinTwoTrackEvents+1, "with two tracks");
	axis->SetBinLabel(AliCDMesonBase::kBinThreeTrackEvents+1,
	                  "with three tracks");
	axis->SetBinLabel(AliCDMesonBase::kBinPionEvents+1, "with two pions");
	axis->SetBinLabel(AliCDMesonBase::kBinKaonEvents+1, "with two kaons");
	axis->SetBinLabel(AliCDMesonBase::kBinProtonEvents+1, "with two proton");
	axis->SetBinLabel(AliCDMesonBase::kBinElectronEvents+1, "with two electron");
	axis->SetBinLabel(AliCDMesonBase::kBinUnknownPIDEvents+1, "with unknown PID");
	axis->SetBinLabel(AliCDMesonBase::kBinResidualTracks+1,
	                  "without residual tracks");
	axis->SetBinLabel(AliCDMesonBase::kBinResidualTracklets+1,
	                  "without residual tracklets");
	axis->SetBinLabel(AliCDMesonBase::kBinCDonlyEvents+1, "CD only events");
	return hist;
}


//------------------------------------------------------------------------------
TH2F* AliCDMesonBase::GetHistPIDStudies(TString name)
{
	//
	// setup the PID studies histogram
	//

	TH2F *hist = new TH2F(name.Data(), ";particle 1;particle 2",
	                      AliCDMesonBase::kBinPIDUnknown,
	                      AliCDMesonBase::kBinPionE,
	                      AliCDMesonBase::kBinPIDUnknown+1,
	                      AliCDMesonBase::kBinPIDUnknown,
	                      AliCDMesonBase::kBinPionE,
	                      AliCDMesonBase::kBinPIDUnknown+1);
	TAxis* x = hist->GetXaxis();
	TAxis* y = hist->GetYaxis();
	x->SetBinLabel(AliCDMesonBase::kBinPionE, "#pi (ex)");
	x->SetBinLabel(AliCDMesonBase::kBinPion, "#pi");
	x->SetBinLabel(AliCDMesonBase::kBinSinglePion, "-");
	x->SetBinLabel(AliCDMesonBase::kBinKaonE, "K (ex)");
	x->SetBinLabel(AliCDMesonBase::kBinKaon, "K");
	x->SetBinLabel(AliCDMesonBase::kBinSingleKaon, ",");
	x->SetBinLabel(AliCDMesonBase::kBinProtonE, "p (ex)");
	x->SetBinLabel(AliCDMesonBase::kBinProton, "p");
	x->SetBinLabel(AliCDMesonBase::kBinSingleProton, "_");
	x->SetBinLabel(AliCDMesonBase::kBinElectronE, "e (ex)");
	x->SetBinLabel(AliCDMesonBase::kBinElectron, "e");
	x->SetBinLabel(AliCDMesonBase::kBinSingleElectron, ".");
	x->SetBinLabel(AliCDMesonBase::kBinPIDUnknown, "X");
	y->SetBinLabel(AliCDMesonBase::kBinPionE, "#pi (ex)");
	y->SetBinLabel(AliCDMesonBase::kBinPion, "#pi");
	y->SetBinLabel(AliCDMesonBase::kBinSinglePion, "-");
	y->SetBinLabel(AliCDMesonBase::kBinKaonE, "K (ex)");
	y->SetBinLabel(AliCDMesonBase::kBinKaon, "K");
	y->SetBinLabel(AliCDMesonBase::kBinSingleKaon, ",");
	y->SetBinLabel(AliCDMesonBase::kBinProtonE, "p (ex)");
	y->SetBinLabel(AliCDMesonBase::kBinProton, "p");
	y->SetBinLabel(AliCDMesonBase::kBinSingleProton, "_");
	y->SetBinLabel(AliCDMesonBase::kBinElectronE, "e (ex)");
	y->SetBinLabel(AliCDMesonBase::kBinElectron, "e");
	y->SetBinLabel(AliCDMesonBase::kBinSingleElectron, ".");
	y->SetBinLabel(AliCDMesonBase::kBinPIDUnknown, "X");
	return hist;
}


//------------------------------------------------------------------------------
TObjArray* AliCDMesonBase::GetHistVZEROStudies(TList* l)
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
void AliCDMesonBase::GetGapTriggers(THnSparseI* gaprun, Int_t gapCondition,
                                    const Int_t run, Double_t& triggers,
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
void AliCDMesonBase::GetNoGapTriggers(THnSparseI* gaprun, Int_t gapCondition,
                                      const Int_t run, Double_t& triggers,
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
