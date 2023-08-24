/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
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
/*  
 *  
 *  AliAnalysisTaskGenUeSpherocity.cxx
 *  
 *
 *  Author: Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)
 *  Original analysistask provided by Gyula BENCEDI  <Gyula.Bencedi@cern.ch>, WIGNER RCP
 * 
 */

//_____ ROOT headers
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include <TF1.h>
#include <TRandom.h>
#include <TTreeStream.h>
#include "TChain.h"
#include <THnSparse.h>
#include <TDatabasePDG.h>
#include "TObjArray.h"
#include <TClonesArray.h>


//_____ ALIROOT headers
#include "AliAnalysisTask.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h" 
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

//_____ Additional includes
#include "AliVEvent.h"
#include "AliGenEventHeader.h"
// #include "AliAnalysisUtils.h"

//_____ AnalysisTask headers
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskGenUeSpherocity.h"

//_____ STL includes
#include <iostream>
using namespace std;

const Char_t * estimators[3]={"Mid05","Mid08","V0M"};
const Int_t NchPercBin=7;
const Int_t NchBin08=122;// multiplicity |eta|<0.8
const Int_t nTSBins=100;
Double_t nchBin_gen0[NchPercBin+1]={0x0};
Double_t nchBin_gen1[NchPercBin+1]={0x0};
Double_t nchBin_gen2[NchPercBin+1]={0x0};
Double_t nchBin_rec0[NchPercBin+1]={0x0};
Double_t nchBin_rec1[NchPercBin+1]={0x0};
Double_t nchBin_rec2[NchPercBin+1]={0x0};

TF1* ch_Eff; // efficiency for pions
TF1* k_neg_Eff;
TF1* k_pos_Eff;
TF1* p_neg_Eff;
TF1* p_pos_Eff;
TF1* xi_Eff;
TF1* phi_Eff;
TF1* k0s_Eff;
TF1* la_Eff;
TF1* labar_Eff;
TF1* k0star_Eff;

const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
const Char_t * pidNames[11] = { "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar", "KStarPM", "SigmaZero" };
Bool_t isPrimary[11] = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE };
const Int_t nPtBins = 66;
Double_t PtBins[nPtBins+1] = {
	0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
	0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
	1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
	2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
	4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
	11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
	26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0 };
const Int_t nPtBins08 = 46;
Double_t PtBins08[nPtBins08+1] = {
	0.15,  0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
	0.5 ,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
	1.0 ,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
	2.0 ,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
	4.0 ,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0};

const Int_t nSoBins = 200;
Double_t SoBins[nSoBins+1]={
	0.000, 0.005, 0.01,  0.015, 0.02,  0.025, 0.03,  0.035, 0.04,  0.045, 
	0.05,  0.055, 0.06,  0.065, 0.07,  0.075, 0.08,  0.085, 0.09,  0.095,
	0.1,   0.105, 0.11,  0.115, 0.12,  0.125, 0.13,  0.135, 0.14,  0.145,
	0.15,  0.155, 0.16,  0.165, 0.17,  0.175, 0.18,  0.185, 0.19,  0.195,
	0.2,   0.205, 0.21,  0.215, 0.22,  0.225, 0.23,  0.235, 0.24,  0.245,
	0.25,  0.255, 0.26,  0.265, 0.27,  0.275, 0.28,  0.285, 0.29,  0.295,
	0.3,   0.305, 0.31,  0.315, 0.32,  0.325, 0.33,  0.335, 0.34,  0.345,
	0.35,  0.355, 0.36,  0.365, 0.37,  0.375, 0.38,  0.385, 0.39,  0.395,
	0.4,   0.405, 0.41,  0.415, 0.42,  0.425, 0.43,  0.435, 0.44,  0.445,
	0.45,  0.455, 0.46,  0.465, 0.47,  0.475, 0.48,  0.485, 0.49,  0.495,
	0.5,   0.505, 0.51,  0.515, 0.52,  0.525, 0.53,  0.535, 0.54,  0.545,
	0.55,  0.555, 0.56,  0.565, 0.57,  0.575, 0.58,  0.585, 0.59,  0.595,
	0.6,   0.605, 0.61,  0.615, 0.62,  0.625, 0.63,  0.635, 0.64,  0.645,
	0.65,  0.655, 0.66,  0.665, 0.67,  0.675, 0.68,  0.685, 0.69,  0.695,
	0.7,   0.705, 0.71,  0.715, 0.72,  0.725, 0.73,  0.735, 0.74,  0.745,
	0.75,  0.755, 0.76,  0.765, 0.77,  0.775, 0.78,  0.785, 0.79,  0.795,
	0.8,   0.805, 0.81,  0.815, 0.82,  0.825, 0.83,  0.835, 0.84,  0.845,
	0.85,  0.855, 0.86,  0.865, 0.87,  0.875, 0.88,  0.885, 0.89,  0.895,
	0.9,   0.905, 0.91,  0.915, 0.92,  0.925, 0.93,  0.935, 0.94,  0.945,
	0.95,  0.955, 0.96,  0.965, 0.97,  0.975, 0.98,  0.985, 0.99,  0.995,
	1.0
};


const Int_t nDphiBins=64;
Double_t DphiBins[nDphiBins+1]={ 
	-1.5708,   -1.47262,  -1.37445,  -1.27627,  -1.1781,   -1.07992,  -0.981748, -0.883573,
	-0.785398, -0.687223, -0.589049, -0.490874, -0.392699, -0.294524, -0.19635,  -0.0981748, 
	0.0000000, 0.0981748, 0.19635,    0.294524,  0.392699,  0.490874,  0.589049,  0.687223,
	0.785398,  0.883573,  0.981748,   1.07992,   1.1781,    1.27627,   1.37445,   1.47262,
	1.5708,    1.66897,   1.76715,    1.86532,   1.9635,    2.06167,   2.15984,   2.25802,
	2.35619,   2.45437,   2.55254,    2.65072,   2.74889,   2.84707,   2.94524,   3.04342,
	3.14159,   3.23977,   3.33794,    3.43612,   3.53429,   3.63247,   3.73064,   3.82882,
	3.92699,   4.02517,   4.12334,    4.22152,   4.31969,   4.41786,   4.51604,   4.61421, 4.7123};

ClassImp( AliAnalysisTaskGenUeSpherocity )

	//_____________________________________________________________________________

	AliAnalysisTaskGenUeSpherocity::AliAnalysisTaskGenUeSpherocity():
		AliAnalysisTaskSE(),
		fMcEvent(0x0),
		fMcHandler(0x0),
		fStack(0),
		fGenerator("Pythia8"),
		fIndexLeadingGen(-1),
		fIndexLeadingRec(-1),
		fMinPtLeading(5.0),
		fSizeStep(0.1),
		//fNso_gen(3),
		//fNso_rec(3),
		fspherocity_gen_ptWeighted(-1),
		fspherocity_gen(-1),
		fspherocity_rec_ptWeighted(-1),
		fspherocity_rec(-1),
		fbinPerc0_gen(0),
		fbinPerc1_gen(0),
		fbinPerc2_gen(0),
		fbinPerc0_rec(0),
		fbinPerc1_rec(0),
		fbinPerc2_rec(0),
		fY(0.5),
		fHistEvt(0x0),
		fHistEta(0x0),
		fHistPart(0x0),
		fDphiNS(0x0),
		fDphiAS(0x0),
		fDphiTS(0x0),
		fMultTS(0x0),
		fDphiNSRec(0x0),
		fDphiASRec(0x0),
		fDphiTSRec(0x0),
		fMultTSRec(0x0),
		fSoWeighedVsNchPtL(0x0),
		fListOfObjects(0)
{

	for(Int_t i_pid=0; i_pid<11; ++i_pid){

		fHistPtVsNchNS[i_pid]=0;
		fHistPtVsNchAS[i_pid]=0;
		fHistPtVsNchTS[i_pid]=0;

		fHistPtVsNchNSRec[i_pid]=0;
		fHistPtVsNchASRec[i_pid]=0;
		fHistPtVsNchTSRec[i_pid]=0;

	}



	for(Int_t i_pid=0; i_pid<11; ++i_pid){

		fHistPt[i_pid]=0;
		fHistPtRec[i_pid]=0;

		for(Int_t i_nch=0; i_nch<NchPercBin; ++i_nch){// loop over mult classes
			// spherocity pt-weighted
			fHistPtVsSoPtW0[i_pid][i_nch]=0;// eta05
			fHistPtVsSoPtWRec0[i_pid][i_nch]=0;
			fHistPtVsSoPtW1[i_pid][i_nch]=0;// eta08
			fHistPtVsSoPtWRec1[i_pid][i_nch]=0;
			fHistPtVsSoPtW2[i_pid][i_nch]=0;// v0m
			fHistPtVsSoPtWRec2[i_pid][i_nch]=0;
			// normalization considering pT=1
			fHistPtVsSo0[i_pid][i_nch]=0;// eta05
			fHistPtVsSoRec0[i_pid][i_nch]=0;
			fHistPtVsSo1[i_pid][i_nch]=0;// eta08
			fHistPtVsSoRec1[i_pid][i_nch]=0;
			fHistPtVsSo2[i_pid][i_nch]=0;// v0m
			fHistPtVsSoRec2[i_pid][i_nch]=0;


		}

	}
	for(Int_t i=0; i<3; ++i){// loop over mult estimators

		fMult[i] = 0;
		fMultRec[i] = 0;

		fNch[i]=0;
		fNchSoSel[i]=0;
		fSoVsNch[i] = 0;
		fSoWeighedVsNch[i] = 0;

		fNchRec[i]=0;
		fNchSoSelRec[i]=0;
		fSoVsNchRec[i] = 0;
		fSoWeighedVsNchRec[i] = 0;

	}

	for(Int_t i_nch08=0; i_nch08<NchBin08; ++i_nch08){

		hDphiSoIS[i_nch08]=0;
		hPtVsSoIS[i_nch08]=0;
		hPtVsSoNS[i_nch08]=0;
		hPtVsSoAS[i_nch08]=0;
		hPtVsSoTS[i_nch08]=0;

	}

	// Default constructor (should not be used)
}

//______________________________________________________________________________

AliAnalysisTaskGenUeSpherocity::AliAnalysisTaskGenUeSpherocity(const char *name):
	AliAnalysisTaskSE(name),
	fMcEvent(0x0),
	fMcHandler(0x0),
	fStack(0),
	fGenerator("Pythia8"),
	fIndexLeadingGen(-1),
	fIndexLeadingRec(-1),
	fMinPtLeading(5.0),
	fSizeStep(0.1),
	//fNso_gen(3),
	//fNso_rec(3),
	fspherocity_gen_ptWeighted(-1),
	fspherocity_gen(-1),
	fspherocity_rec_ptWeighted(-1),
	fspherocity_rec(-1),
	fbinPerc0_gen(0),
	fbinPerc1_gen(0),
	fbinPerc2_gen(0),
	fbinPerc0_rec(0),
	fbinPerc1_rec(0),
	fbinPerc2_rec(0),
	fY(0.5),
	fHistEvt(0x0),
	fHistEta(0x0),
	fHistPart(0x0),
	fDphiNS(0x0),
	fDphiAS(0x0),
	fDphiTS(0x0),
	fMultTS(0x0),
	fDphiNSRec(0x0),
	fDphiASRec(0x0),
	fDphiTSRec(0x0),
	fMultTSRec(0x0),
	fSoWeighedVsNchPtL(0x0),
	fListOfObjects(0)
{

	for(Int_t i_pid=0; i_pid<11; ++i_pid){

		fHistPtVsNchNS[i_pid]=0;
		fHistPtVsNchAS[i_pid]=0;
		fHistPtVsNchTS[i_pid]=0;

		fHistPtVsNchNSRec[i_pid]=0;
		fHistPtVsNchASRec[i_pid]=0;
		fHistPtVsNchTSRec[i_pid]=0;

	}

	for(Int_t i_pid=0; i_pid<11; ++i_pid){

		fHistPt[i_pid]=0;
		fHistPtRec[i_pid]=0;
		for(Int_t i_nch=0; i_nch<NchPercBin; ++i_nch){// loop over mult classes

			fHistPtVsSoPtW0[i_pid][i_nch]=0;// eta05
			fHistPtVsSoPtWRec0[i_pid][i_nch]=0;
			fHistPtVsSoPtW1[i_pid][i_nch]=0;// eta08
			fHistPtVsSoPtWRec1[i_pid][i_nch]=0;
			fHistPtVsSoPtW2[i_pid][i_nch]=0;// v0m
			fHistPtVsSoPtWRec2[i_pid][i_nch]=0;

			fHistPtVsSo0[i_pid][i_nch]=0;// eta05
			fHistPtVsSoRec0[i_pid][i_nch]=0;
			fHistPtVsSo1[i_pid][i_nch]=0;// eta08
			fHistPtVsSoRec1[i_pid][i_nch]=0;
			fHistPtVsSo2[i_pid][i_nch]=0;// v0m
			fHistPtVsSoRec2[i_pid][i_nch]=0;

		}

	}
	for(Int_t i=0; i<3; ++i){
		fMult[i] = 0;
		fMultRec[i] = 0;

		fNch[i]=0;
		fNchSoSel[i]=0;
		fSoVsNch[i] = 0;
		fSoWeighedVsNch[i] = 0;

		fNchRec[i]=0;
		fNchSoSelRec[i]=0;
		fSoVsNchRec[i] = 0;
		fSoWeighedVsNchRec[i] = 0;


	}
	for(Int_t i_nch08=0; i_nch08<NchBin08; ++i_nch08){

		hDphiSoIS[i_nch08]=0;
		hPtVsSoIS[i_nch08]=0;
		hPtVsSoNS[i_nch08]=0;
		hPtVsSoAS[i_nch08]=0;
		hPtVsSoTS[i_nch08]=0;

	}
	DefineInput( 0, TChain::Class());
	DefineOutput(1, TList::Class() ); // Basic output slot 
}


//_____________________________________________________________________________

AliAnalysisTaskGenUeSpherocity::~AliAnalysisTaskGenUeSpherocity(){
	// Destructor
	// histograms are in the output list and deleted when the output
	// list is deleted by the TSelector dtor
	if (fListOfObjects) { delete fListOfObjects; fListOfObjects=0x0; }
}

//______________________________________________________________________________

void AliAnalysisTaskGenUeSpherocity::UserCreateOutputObjects(){

	ch_Eff = 0;
	ch_Eff = new TF1("ch_Eff",
			"(x>=0.15&&x<[0])*([1]+x*[2])+(x>=[0]&&x<[3])*([4]+[5]*x+[6]*x*x)+(x>=[3])*([7])", 0.0, 1e2);
	ch_Eff->SetParameters(0.4,0.559616,0.634754,1.5,0.710963,0.312747,-0.163094,0.813976);

	k_neg_Eff = 0;
	k_neg_Eff = new TF1("k_neg_Eff_Eff",
			"(x>=0.25&&x<[0])*([1]+[2]*x+[3]*x*x+[4]*x*x*x+[5]*x*x*x*x)+(x>=[0]&&x<[6])*([7]+[8]*x)+(x>=[6])*([9])", 0.0, 1e3);
	k_neg_Eff->SetParameters(1.4, -0.437094, 4.03935, -6.10465, 4.41681, -1.21593, 3.6, 0.693548, 0.0206902,
			0.784834);

	k_pos_Eff = 0;
	k_pos_Eff = new TF1("k_pos_Eff_Eff",
			"(x>=0.25&&x<[0])*([1]+[2]*x+[3]*x*x+[4]*x*x*x+[5]*x*x*x*x)+(x>=[0]&&x<[6])*([7]+[8]*x)+(x>=[6])*([9])", 0.0, 1e3);
	k_pos_Eff->SetParameters(1.0, -0.437094, 4.03935, -6.10465, 4.41681, -1.21593, 3.0, 0.694729, 0.0238144,
			0.784834);

	p_neg_Eff = 0;
	p_neg_Eff = new TF1("p_neg_Eff_Eff",
			"(x>=0.3&&x<[0])*([1]+[2]*x)+(x>=[0]&&x<[3])*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)+(x>=[3])*([9])", 0.0, 1e3);
	p_neg_Eff->SetParameters(0.4, -1.49693, 5.55626, 2.4, 0.21432, 2.0683, -2.28951, 1.04733, -0.171858,
			0.802333);

	p_pos_Eff = 0; 
	p_pos_Eff = new TF1("p_pos_Eff_Eff",
			"(x>=0.3&&x<[0])*([1]+[2]*x)+(x>=[0]&&x<[3])*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)+(x>=[3])*([9])", 0.0, 1e3);
	p_pos_Eff->SetParameters(0.35, -1.49693, 5.55626, 2.2, 0.477826, 1.10708, -1.08169, 0.419493, -0.0565612,
			0.802333);

	xi_Eff = 0;
	xi_Eff = new TF1("xi_Eff", "0*(x<[0]) + ([1]*(x-[0])+[2]*(x-[0])*(x-[0]))*([0]<=x&&x<[3])+[4]*(1-[5]/x)*([3]<=x)", 0.0, 1e3);
	xi_Eff->SetParameters(0.643, 0.114, -0.00594, 3.06, 0.283, 0.489);

	phi_Eff = 0;
	phi_Eff = new TF1("phi_Eff","([0]*x^[1]*exp(-x) - [2]*x +[3])*(x>[4])",0,1e3);
	phi_Eff->SetParameters(-5.25411e-01,5.42076e-02,-2.99137e-03,3.03124e-01,2.77814e-02);

	k0s_Eff = 0;
	k0s_Eff = new TF1("k0s_Eff","([0]*x^[1]*exp(-x) - [2]*x +[3])*(x>[4])",0,20);
	k0s_Eff->SetParameters(-7.17967e-01,1.27198e-01,1.51269e-02,5.19441e-01,0.1);

	la_Eff = 0;
	la_Eff = new TF1("la_Eff","([0]*x^[1]*exp(-x) - [2]*x +[3])*(x>[4])",0,20);
	la_Eff->SetParameters(-5.46924e-01,8.50191e-03,1.49558e-02,3.96383e-01,0.35);

	labar_Eff = 0;
	labar_Eff = new TF1("labar_Eff","([0]*x^[1]*exp(-x) - [2]*x +[3])*(x>[4])",0,20);
	labar_Eff->SetParameters(-5.42213e-01,2.16541e-02,1.32033e-02,3.79549e-01,0.35);

	k0star_Eff = 0;
	k0star_Eff = new TF1("k0star_Eff","(([0]+[1]*x+[2]*x*x+[3]*x*x*x)*(x>0.4&&x<3.0))+(([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)*(x>3.0&&x<15.0))",0,20);
	k0star_Eff->SetParameters(-0.06334,0.3354,-0.03492,-0.003904,0.2856,0.1316,-0.02146,0.001467,-3.563e-05);

	// ### Analysis output
	fListOfObjects = new TList();
	fListOfObjects->SetOwner(kTRUE);

	// ### Create histograms
	fHistEvt = 0;	
	fHistEvt = new TH1I("fHistEvt","fHistEvt",2,0,2) ;
	fHistEvt->GetXaxis()->SetBinLabel(1,"All events");
	fHistEvt->GetXaxis()->SetBinLabel(2,"All particles");
	fHistEvt->Sumw2();
	fListOfObjects->Add(fHistEvt);

	fHistEta = 0;
	fHistEta = new TH1F("fHistEta","Eta Distr.; #eta; N_{part}", 200, -1., 1.);
	fListOfObjects->Add(fHistEta);

	InitHisto<TH1F>("fHistY", "Y Distr.", 200, -1., 1., "#it{y}", "N_{part}");
	InitHisto<TH1F>("fHistYRec", "Y Distr. rec", 200, -1., 1., "#it{y}", "N_{part}");

	fDphiNS = 0;
	fDphiNS = new TH1D("hDphiNS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiNS);

	fDphiAS = 0;
	fDphiAS = new TH1D("hDphiAS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiAS);

	fDphiTS = 0;
	fDphiTS = new TH1D("hDphiTS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiTS);

	fMultTS = 0;
	fMultTS = new TH1D("fMultTS","",100,-0.5,99.5);
	fListOfObjects->Add(fMultTS);

	fDphiNSRec = 0;
	fDphiNSRec = new TH1D("hDphiNSRec","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiNSRec);

	fDphiASRec = 0;
	fDphiASRec = new TH1D("hDphiASRec","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiASRec);

	fDphiTSRec = 0;
	fDphiTSRec = new TH1D("hDphiTSRec","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiTSRec);

	fMultTSRec = 0;
	fMultTSRec = new TH1D("fMultTSRec","",100,-0.5,99.5);
	fListOfObjects->Add(fMultTSRec);

	fSoWeighedVsNchPtL = 0;
	fSoWeighedVsNchPtL = new TH2D("fSoWeighedVsNchPtL","So vs Nch |eta|<0.8 pTL>5",300,-0.5,299.5,200,0.0,1.0);
	fListOfObjects->Add(fSoWeighedVsNchPtL);


	Double_t TSBins[nTSBins+1]={0x0};
	for(Int_t i=0;i<nTSBins;++i){
		TSBins[i]=i*1.0-0.5;
	}
	TSBins[nTSBins]=99.5;

	for(Int_t i_pid=0; i_pid<11; ++i_pid){

		fHistPtVsNchNS[i_pid]=0;
		fHistPtVsNchNS[i_pid]=new TH2D(Form("fHistPtVsNchNS_%s",pidNames[i_pid]), "Generated #it{p}_{T} distribution NS",nTSBins, TSBins, nPtBins,PtBins);
		fListOfObjects->Add(fHistPtVsNchNS[i_pid]);

		fHistPtVsNchAS[i_pid]=0;
		fHistPtVsNchAS[i_pid]=new TH2D(Form("fHistPtVsNchAS_%s",pidNames[i_pid]), "Generated #it{p}_{T} distribution AS",nTSBins, TSBins, nPtBins,PtBins);
		fListOfObjects->Add(fHistPtVsNchAS[i_pid]);

		fHistPtVsNchTS[i_pid]=0;
		fHistPtVsNchTS[i_pid]=new TH2D(Form("fHistPtVsNchTS_%s",pidNames[i_pid]), "Generated #it{p}_{T} distribution TS",nTSBins, TSBins, nPtBins,PtBins);
		fListOfObjects->Add(fHistPtVsNchTS[i_pid]);

		fHistPtVsNchNSRec[i_pid]=0;
		fHistPtVsNchNSRec[i_pid]=new TH2D(Form("fHistPtVsNchNSRec_%s",pidNames[i_pid]), "Rec #it{p}_{T} distribution NS",nTSBins, TSBins, nPtBins,PtBins);
		fListOfObjects->Add(fHistPtVsNchNSRec[i_pid]);

		fHistPtVsNchASRec[i_pid]=0;
		fHistPtVsNchASRec[i_pid]=new TH2D(Form("fHistPtVsNchASRec_%s",pidNames[i_pid]), "Rec #it{p}_{T} distribution AS",nTSBins, TSBins, nPtBins,PtBins);
		fListOfObjects->Add(fHistPtVsNchASRec[i_pid]);

		fHistPtVsNchTSRec[i_pid]=0;
		fHistPtVsNchTSRec[i_pid]=new TH2D(Form("fHistPtVsNchTSRec_%s",pidNames[i_pid]), "Rec #it{p}_{T} distribution TS",nTSBins, TSBins, nPtBins,PtBins);
		fListOfObjects->Add(fHistPtVsNchTSRec[i_pid]);

	}



	for(Int_t i_pid=0; i_pid<11; ++i_pid){

		fHistPt[i_pid]=0;
		fHistPt[i_pid]=new TH1F(Form("fHistPt_%s",pidNames[i_pid]), "Generated #it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c}); Entries",nPtBins,PtBins);
		fListOfObjects->Add(fHistPt[i_pid]);

		fHistPtRec[i_pid]=0;
		fHistPtRec[i_pid]=new TH1F(Form("fHistPtRec_%s",pidNames[i_pid]), "Rec #it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c}); Entries",nPtBins,PtBins);
		fListOfObjects->Add(fHistPtRec[i_pid]);


		for(Int_t i_nch=0; i_nch<NchPercBin; ++i_nch){// loop over mult classes

			fHistPtVsSoPtW0[i_pid][i_nch]=0;
			fHistPtVsSoPtW0[i_pid][i_nch]=new TH2D(Form("fHistPtVsSoPtW0_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSoPtW0[i_pid][i_nch]);

			fHistPtVsSoPtWRec0[i_pid][i_nch]=0;
			fHistPtVsSoPtWRec0[i_pid][i_nch]=new TH2D(Form("fHistPtVsSoPtWRec0_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSoPtWRec0[i_pid][i_nch]);

			fHistPtVsSoPtW1[i_pid][i_nch]=0;
			fHistPtVsSoPtW1[i_pid][i_nch]=new TH2D(Form("fHistPtVsSoPtW1_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSoPtW1[i_pid][i_nch]);

			fHistPtVsSoPtWRec1[i_pid][i_nch]=0;
			fHistPtVsSoPtWRec1[i_pid][i_nch]=new TH2D(Form("fHistPtVsSoPtWRec1_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSoPtWRec1[i_pid][i_nch]);

			fHistPtVsSoPtW2[i_pid][i_nch]=0;
			fHistPtVsSoPtW2[i_pid][i_nch]=new TH2D(Form("fHistPtVsSoPtW2_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSoPtW2[i_pid][i_nch]);

			fHistPtVsSoPtWRec2[i_pid][i_nch]=0;
			fHistPtVsSoPtWRec2[i_pid][i_nch]=new TH2D(Form("fHistPtVsSoPtWRec2_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSoPtWRec2[i_pid][i_nch]);

			// normalization considering pT=1
			fHistPtVsSo0[i_pid][i_nch]=0;
			fHistPtVsSo0[i_pid][i_nch]=new TH2D(Form("fHistPtVsSo0_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSo0[i_pid][i_nch]);

			fHistPtVsSoRec0[i_pid][i_nch]=0;
			fHistPtVsSoRec0[i_pid][i_nch]=new TH2D(Form("fHistPtVsSoRec0_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSoRec0[i_pid][i_nch]);

			fHistPtVsSo1[i_pid][i_nch]=0;
			fHistPtVsSo1[i_pid][i_nch]=new TH2D(Form("fHistPtVsSo1_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSo1[i_pid][i_nch]);

			fHistPtVsSoRec1[i_pid][i_nch]=0;
			fHistPtVsSoRec1[i_pid][i_nch]=new TH2D(Form("fHistPtVsSoRec1_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSoRec1[i_pid][i_nch]);

			fHistPtVsSo2[i_pid][i_nch]=0;
			fHistPtVsSo2[i_pid][i_nch]=new TH2D(Form("fHistPtVsSo2_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSo2[i_pid][i_nch]);

			fHistPtVsSoRec2[i_pid][i_nch]=0;
			fHistPtVsSoRec2[i_pid][i_nch]=new TH2D(Form("fHistPtVsSoRec2_%s_%d",pidNames[i_pid],i_nch),"",nSoBins,SoBins,nPtBins,PtBins);
			fListOfObjects->Add(fHistPtVsSoRec2[i_pid][i_nch]);

		}



	}

	for(Int_t i=0; i<3; ++i){

		fMult[i] = 0;
		fMult[i] = new TH2D(Form("fMult_%s",estimators[i]),Form("Selection bias %s; Nch; #eta",estimators[i]),300,-0.5,299.5,100,-5,5); 
		fListOfObjects->Add(fMult[i]);

		fMultRec[i] = 0;
		fMultRec[i] = new TH2D(Form("fMultRec_%s",estimators[i]),Form("Selection bias %s; Nch; #eta",estimators[i]),300,-0.5,299.5,100,-5,5); 
		fListOfObjects->Add(fMultRec[i]);

		fSoVsNch[i] = 0;
		fSoVsNch[i] = new TH2D(Form("fSoVsNch_%s",estimators[i]),Form("So vs Nch %s",estimators[i]),300,-0.5,299.5,200,0.0,1.0);
		fListOfObjects->Add(fSoVsNch[i]);

		fSoWeighedVsNch[i] = 0;
		fSoWeighedVsNch[i] = new TH2D(Form("fSoWeighedVsNch_%s",estimators[i]),Form("So vs Nch %s",estimators[i]),300,-0.5,299.5,200,0.0,1.0);
		fListOfObjects->Add(fSoWeighedVsNch[i]);

		fNch[i]=0;
		fNch[i]=new TH1F(Form("fNch_%s",estimators[i]),Form("Nch %s",estimators[i]),300,-0.5,299.5);
		fListOfObjects->Add(fNch[i]);

		fNchSoSel[i]=0;
		fNchSoSel[i]=new TH1F(Form("fNchSoSel_%s",estimators[i]),Form("Nch %s (after So cuts)",estimators[i]),300,-0.5,299.5);
		fListOfObjects->Add(fNchSoSel[i]);



		fSoVsNchRec[i] = 0;
		fSoVsNchRec[i] = new TH2D(Form("fSoVsNchRec_%s",estimators[i]),Form("So vs Nch %s",estimators[i]),300,-0.5,299.5,200,0.0,1.0);
		fListOfObjects->Add(fSoVsNchRec[i]);

		fSoWeighedVsNchRec[i] = 0;
		fSoWeighedVsNchRec[i] = new TH2D(Form("fSoWeighedVsNchRec_%s",estimators[i]),Form("So vs Nch %s",estimators[i]),300,-0.5,299.5,200,0.0,1.0);
		fListOfObjects->Add(fSoWeighedVsNchRec[i]);

		fNchRec[i]=0;
		fNchRec[i]=new TH1F(Form("fNchRec_%s",estimators[i]),Form("Nch %s",estimators[i]),300,-0.5,299.5);
		fListOfObjects->Add(fNchRec[i]);

		fNchSoSelRec[i]=0;
		fNchSoSelRec[i]=new TH1F(Form("fNchSoSelRec_%s",estimators[i]),Form("Nch %s (after So cuts)",estimators[i]),300,-0.5,299.5);
		fListOfObjects->Add(fNchSoSelRec[i]);


	}

	for(Int_t i_nch08=0; i_nch08<NchBin08; ++i_nch08){

		hDphiSoIS[i_nch08]=0;
		hDphiSoIS[i_nch08]=new TH2D(Form("hDphiVsSoIS_Nch08_%d",i_nch08+3),Form("hDphiVsSoIS_Nch08_%d",i_nch08+3),nSoBins,SoBins,nDphiBins,DphiBins);
		fListOfObjects->Add(hDphiSoIS[i_nch08]);

		hPtVsSoIS[i_nch08]=0;
		hPtVsSoIS[i_nch08]=new TH2D(Form("hPtVsSoIS_Nch08_%d",i_nch08+3),Form("hPtVsSoIS_Nch08_%d",i_nch08+3),nSoBins,SoBins,nPtBins08,PtBins08);
		fListOfObjects->Add(hPtVsSoIS[i_nch08]);

		hPtVsSoNS[i_nch08]=0;
		hPtVsSoNS[i_nch08]=new TH2D(Form("hPtVsSoNS_Nch08_%d",i_nch08+3),Form("hPtVsSoNS_Nch08_%d",i_nch08+3),nSoBins,SoBins,nPtBins08,PtBins08);
		fListOfObjects->Add(hPtVsSoNS[i_nch08]);

		hPtVsSoAS[i_nch08]=0;
		hPtVsSoAS[i_nch08]=new TH2D(Form("hPtVsSoAS_Nch08_%d",i_nch08+3),Form("hPtVsSoAS_Nch08_%d",i_nch08+3),nSoBins,SoBins,nPtBins08,PtBins08);
		fListOfObjects->Add(hPtVsSoAS[i_nch08]);

		hPtVsSoTS[i_nch08]=0;
		hPtVsSoTS[i_nch08]=new TH2D(Form("hPtVsSoTS_Nch08_%d",i_nch08+3),Form("hPtVsSoTS_Nch08_%d",i_nch08+3),nSoBins,SoBins,nPtBins08,PtBins08);
		fListOfObjects->Add(hPtVsSoTS[i_nch08]);

	}


	// ### List of outputs
	PostData(1, fListOfObjects);

}

//______________________________________________________________________________

inline void AliAnalysisTaskGenUeSpherocity::FillHisto(const char* objkey, Double_t x)
{
	TH1* hTmp = 0;
	hTmp = static_cast<TH1*>(fListOfObjects->FindObject(objkey));
	if(!hTmp){
		AliError(Form("Cannot find histogram: %s",objkey)) ;
		return;
	}
	hTmp->Fill(x);
}

inline void AliAnalysisTaskGenUeSpherocity::FillHisto(const char* objkey, Double_t x, Double_t y)
{
	TH2* hTmp = 0;
	hTmp = static_cast<TH2*>(fListOfObjects->FindObject(objkey));
	if(!hTmp){
		AliError(Form("Cannot find histogram: %s",objkey)) ;
		return;
	}
	hTmp->Fill(x,y);
}

//______________________________________________________________________________

template <class T> T* AliAnalysisTaskGenUeSpherocity::InitHisto(const char* hname, const char* htitle, Int_t nxbins, Double_t xmin, Double_t xmax, const char* xtitle, const char* ytitle)
{
	T* hTmp = 0;
	hTmp = new T(hname, htitle, nxbins, xmin, xmax);
	hTmp->GetXaxis()->SetTitle(xtitle);
	hTmp->GetYaxis()->SetTitle(ytitle);
	//	hTmp->SetMarkerStyle(kFullCircle);
	//	hTmp->Sumw2();
	fListOfObjects->Add(hTmp);

	return hTmp;
}
template <class T> T* AliAnalysisTaskGenUeSpherocity::InitHisto(const char* hname, const char* htitle, Int_t nxbins, Double_t xmin, Double_t xmax, Int_t nybins, Double_t ymin, Double_t ymax, const char* xtitle, const char* ytitle)
{
	T* hTmp = 0;
	hTmp = new T(hname, htitle, nxbins, xmin, xmax, nybins, ymin, ymax);
	hTmp->GetXaxis()->SetTitle(xtitle);
	hTmp->GetYaxis()->SetTitle(ytitle);
	//	hTmp->Sumw2();
	fListOfObjects->Add(hTmp);

	return hTmp;
}

//______________________________________________________________________________

void AliAnalysisTaskGenUeSpherocity::Init(){
	//
	fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

	if(fGenerator == "Pythia6"){
		Double_t nchBin_gen0tmp[NchPercBin+1]={-0.5, 1.5, 3.5, 5.5, 7.5, 10.5, 15.5, 299.5};
		Double_t nchBin_gen1tmp[NchPercBin+1]={-0.5, 3.5, 6.5, 9.5, 13.5, 17.5, 25.5, 299.5};
		Double_t nchBin_gen2tmp[NchPercBin+1]={-0.5, 10.5, 18.5, 24.5, 31.5, 42.5, 57.5, 299.5};
		Double_t nchBin_rec0tmp[NchPercBin+1]={-0.5, 1.5, 2.5, 4.5, 5.5, 8.5, 11.5, 299.5};
		Double_t nchBin_rec1tmp[NchPercBin+1]={-0.5, 2.5, 4.5, 6.5, 9.5, 13.5, 18.5, 299.5};
		Double_t nchBin_rec2tmp[NchPercBin+1]={-0.5, 11.5, 19.5, 25.5, 32.5, 43.5, 57.5, 299.5};
		for(Int_t i=0;i<NchPercBin+1;++i){
			nchBin_gen0[i]=nchBin_gen0tmp[i];
			nchBin_gen1[i]=nchBin_gen1tmp[i];
			nchBin_gen2[i]=nchBin_gen2tmp[i];

			nchBin_rec0[i]=nchBin_rec0tmp[i];
			nchBin_rec1[i]=nchBin_rec1tmp[i];
			nchBin_rec2[i]=nchBin_rec2tmp[i];
		}
	}
	else if(fGenerator == "Epos"){
		Double_t nchBin_gen0tmp[NchPercBin+1]={-0.5, 1.5, 4.5, 5.5, 7.5, 10.5, 16.5, 299.5};
		Double_t nchBin_gen1tmp[NchPercBin+1]={-0.5, 3.5, 6.5, 9.5, 12.5, 17.5, 25.5, 299.5};
		Double_t nchBin_gen2tmp[NchPercBin+1]={-0.5, 12.5, 21.5, 27.5, 36.5, 49.5, 68.5, 299.5};
		Double_t nchBin_rec0tmp[NchPercBin+1]={-0.5, 1.5, 2.5, 3.5, 5.5, 7.5, 11.5, 299.5};
		Double_t nchBin_rec1tmp[NchPercBin+1]={-0.5, 2.5, 4.5, 6.5, 9.5, 12.5, 18.5, 299.5};
		Double_t nchBin_rec2tmp[NchPercBin+1]={-0.5, 13.5, 21.5, 28.5, 37.5, 50.5, 69.5, 299.5};
		for(Int_t i=0;i<NchPercBin+1;++i){
			nchBin_gen0[i]=nchBin_gen0tmp[i];
			nchBin_gen1[i]=nchBin_gen1tmp[i];
			nchBin_gen2[i]=nchBin_gen2tmp[i];
			nchBin_rec0[i]=nchBin_rec0tmp[i];
			nchBin_rec1[i]=nchBin_rec1tmp[i];
			nchBin_rec2[i]=nchBin_rec2tmp[i];
		}
	}
	else{
		Double_t nchBin_gen0tmp[NchPercBin+1]={-0.5, 1.5, 4.5, 5.5, 7.5, 10.5, 15.5, 299.5};
		Double_t nchBin_gen1tmp[NchPercBin+1]={-0.5, 3.5, 6.5, 9.5, 12.5, 17.5, 25.5, 299.5};
		Double_t nchBin_gen2tmp[NchPercBin+1]={-0.5, 12.5, 19.5, 26.5, 34.5, 46.5, 64.5, 299.5};
		Double_t nchBin_rec0tmp[NchPercBin+1]={-0.5, 1.5, 2.5, 3.5, 5.5, 7.5, 11.5, 299.5};
		Double_t nchBin_rec1tmp[NchPercBin+1]={-0.5, 2.5, 4.5, 6.5, 9.5, 12.5, 18.5, 299.5};
		Double_t nchBin_rec2tmp[NchPercBin+1]={-0.5, 12.5, 20.5, 26.5, 35.5, 47.5, 65.5, 299.5};
		for(Int_t i=0;i<NchPercBin+1;++i){
			nchBin_gen0[i]=nchBin_gen0tmp[i];
			nchBin_gen1[i]=nchBin_gen1tmp[i];
			nchBin_gen2[i]=nchBin_gen2tmp[i];
			nchBin_rec0[i]=nchBin_rec0tmp[i];
			nchBin_rec1[i]=nchBin_rec1tmp[i];
			nchBin_rec2[i]=nchBin_rec2tmp[i];
		}
	}
}

//______________________________________________________________________________

void AliAnalysisTaskGenUeSpherocity::UserExec(Option_t *){

	// ### Initialize
	Init();

	// ### MC handler
	if(fMcHandler)
		fMcEvent = fMcHandler->MCEvent();
	else { if(fDebug > 1) printf("AliAnalysisTaskGenUeSpherocity::Handler() fMcHandler = NULL\n"); return; }

	// ### MC event
	if( !fMcEvent ) { if(fDebug > 1) printf("AliAnalysisTaskGenUeSpherocity::UserExec() fMcEvent = NULL \n"); return; }

	fStack = ((AliMCEvent*)fMcEvent)->Stack();
	if (!fStack) {
		Printf("ERROR: Could not retrieve MC stack \n");
		cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
		return;
	}

	// ### MC event selection
	Bool_t isEventMCSelected = IsMCEventSelected(fMcEvent);
	if( !isEventMCSelected ) return;


	// GENERATOR LEVEL
	vector<Int_t> mult_estimators_gen;
	vector<Float_t> pt_so_gen;
	vector<Float_t> eta_so_gen;
	vector<Float_t> phi_so_gen;


	Int_t fNso_gen = -1;
	fNso_gen = GetMultipliciy( kFALSE, mult_estimators_gen, pt_so_gen, eta_so_gen, phi_so_gen );

	// PSEUDO REC LEVEL
	vector<Int_t> mult_estimators_rec;
	vector<Float_t> pt_so_rec;
	vector<Float_t> eta_so_rec;
	vector<Float_t> phi_so_rec;

	// RT analysis
	fIndexLeadingGen = -1;
	fIndexLeadingGen = GetIndexLeading(kFALSE);
	TParticle* mcPartLeadingGen         = 0x0;
	if(fIndexLeadingGen>=0){
		mcPartLeadingGen                    = (TParticle *)fMcEvent->Particle(fIndexLeadingGen);
		if(mcPartLeadingGen->Pt()>=fMinPtLeading){
			MakeRTAnalysis(kFALSE);
		}
	}

	fIndexLeadingRec = -1;
	fIndexLeadingRec = GetIndexLeading(kTRUE);
	TParticle* mcPartLeadingRec         = 0x0;
	if(fIndexLeadingRec>=0){
		mcPartLeadingRec                    = (TParticle *)fMcEvent->Particle(fIndexLeadingRec);
		if(mcPartLeadingRec->Pt()>=fMinPtLeading){
			MakeRTAnalysis(kTRUE);
		}
	}


	Int_t fNso_rec = -1;
	fNso_rec = GetMultipliciy( kTRUE, mult_estimators_rec, pt_so_rec, eta_so_rec, phi_so_rec );

	// mult_estimators[3]: mult in |eta|<1
	Bool_t fIsInel0_gen = kFALSE;
	if(mult_estimators_gen[3]>0)
		fIsInel0_gen = kTRUE;// is INEL>0


	Bool_t fIsInel0_rec = kFALSE;
	if(mult_estimators_rec[3]>0)
		fIsInel0_rec = kTRUE;// is INEL>0
	//////////////////////////////////////////////////////////////////////
	///////    Temporal solutution, Ntrk>10, pT>0.15 (to match the data analysis)
	if(fNso_gen<10)
		return;
	///////////////////////////////////////////////////////////////////////

	if(fIsInel0_gen)
		MakeAnaGen(fNso_gen, mult_estimators_gen, pt_so_gen, eta_so_gen, phi_so_gen);

	if(fIsInel0_rec)
		MakeAnaRec(fNso_rec, mult_estimators_rec, pt_so_rec, eta_so_rec, phi_so_rec);

	fbinPerc0_gen=GetMultBin(kFALSE,mult_estimators_gen[0],0);
	fbinPerc1_gen=GetMultBin(kFALSE,mult_estimators_gen[1],1);
	fbinPerc2_gen=GetMultBin(kFALSE,mult_estimators_gen[2],2);


	fbinPerc0_rec=GetMultBin(kTRUE,mult_estimators_rec[0],0);
	fbinPerc1_rec=GetMultBin(kTRUE,mult_estimators_rec[1],1);
	fbinPerc2_rec=GetMultBin(kTRUE,mult_estimators_rec[2],2);


	// Only for Sphero (pT=1) Analysis
	if(fIsInel0_gen&&fIsInel0_rec){
		ParticleSel( kTRUE, mult_estimators_rec );
		ParticleSel( kFALSE, mult_estimators_gen );
	}

	// To check charged particle spherocity analysis
	if(fIndexLeadingGen>=0){
		if(mcPartLeadingGen->Pt()>=0.15){
			MakeUeSoNch08Analysis(mult_estimators_gen);
		}
	}
	//MemInfo_t memInfo;
	//Int_t memUsage = 0;
	//gSystem->GetMemInfo(&memInfo);
	//memUsage = memInfo.fMemUsed;
	//cout<<"mem usage="<<memUsage<<endl;
	// ### Post data for all output slots
	mult_estimators_rec.clear();
	pt_so_rec.clear();
	eta_so_rec.clear();
	phi_so_rec.clear();

	mult_estimators_gen.clear();
	pt_so_gen.clear();
	eta_so_gen.clear();
	phi_so_gen.clear();


	PostData(1, fListOfObjects);

	return;
}

//______________________________________________________________________________

Bool_t AliAnalysisTaskGenUeSpherocity::IsMCEventSelected(TObject* obj){

	Bool_t isSelected = kTRUE;

	AliMCEvent *event = 0x0;
	event = dynamic_cast<AliMCEvent*>(obj);
	if( !event ) 
		isSelected = kFALSE;

	if( isSelected ) 
		FillHisto("fHistEvt",0.5);

	return isSelected;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGenUeSpherocity::GetIndexLeading(Bool_t fIsPseudoRec){

	Double_t ptleading = 0;
	Int_t index_leading = -1;
	Bool_t isPhysPrim = kFALSE;
	Double_t qPart = 0;
	Double_t etaPart = -10;
	Int_t pidCodeMC = 0;
	Int_t pPDG = -10;

	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;
		//selection of primary charged particles
		if(!(mcPart->GetPDG())) continue;
		qPart = mcPart->GetPDG()->Charge()/3.;
		if(TMath::Abs(qPart)<0.001) continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		if(!isPhysPrim)
			continue;

		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);

		if(TMath::Abs(pidCodeMC)==5)
			continue;

		etaPart = mcPart -> Eta();
		if(fIsPseudoRec)
			if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),mcPart->Pt())) continue;
		if(TMath::Abs(etaPart) > 0.8)continue;
		if(mcPart -> Pt()<1.0)continue;
		if(mcPart -> Pt()>ptleading){
			ptleading = mcPart->Pt();
			index_leading = ipart;
		}


	} // particle loop
	cout<<"\n";
	return index_leading;
}

Int_t AliAnalysisTaskGenUeSpherocity::GetMultipliciy(Bool_t fIsPseudoRec, vector<Int_t> &multArray, vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray){


	Int_t mult_so=0;
	multArray.clear();
	ptArray.clear();
	etaArray.clear();
	phiArray.clear();

	Bool_t isPhysPrim = kFALSE;

	Int_t mult_Eta5   = 0;
	Int_t mult_Eta8   = 0;
	Int_t mult_Eta1   = 0;
	Int_t mult_VZEROM = 0;
	Double_t qPart = 0;
	Double_t etaPart = -10;


	Int_t pidCodeMC = 0;
	Int_t pPDG = -10;

	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;
		//selection of primary charged particles
		if(!(mcPart->GetPDG())) continue;
		qPart = mcPart->GetPDG()->Charge()/3.;
		if(TMath::Abs(qPart)<0.001) continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		if(!isPhysPrim)
			continue;

		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);

		etaPart = mcPart -> Eta();
		if( (2.8 < etaPart && etaPart < 5.1) || (-3.7 < etaPart && etaPart <-1.7) ) mult_VZEROM++;

		if(fIsPseudoRec)
			if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),mcPart->Pt())) continue;

		if( TMath::Abs(etaPart) < 0.5 ) mult_Eta5++;
		if( TMath::Abs(etaPart) < 0.8 ){ 
			mult_Eta8++;

			if(mcPart -> Pt()>0.15){

				ptArray.push_back(mcPart->Pt());
				etaArray.push_back(mcPart->Eta());
				phiArray.push_back(mcPart->Phi());
				mult_so++;
			}

		}
		if( TMath::Abs(etaPart) < 1 )
			if(mcPart -> Pt()>0)
				mult_Eta1++;// for INEL>0n

	} // particle loop

	multArray.push_back(mult_Eta5);
	multArray.push_back(mult_Eta8);
	multArray.push_back(mult_VZEROM);
	multArray.push_back(mult_Eta1);

	return mult_so;
}

//______________________________________________________________________________
void AliAnalysisTaskGenUeSpherocity::MakeRTAnalysis(Bool_t fIsPseudoRec){

	// Properties leading particle
	TParticle* mcPartTmp         = 0x0;
	if(fIsPseudoRec)
		mcPartTmp                    = (TParticle *)fMcEvent->Particle(fIndexLeadingRec);
	else
		mcPartTmp                    = (TParticle *)fMcEvent->Particle(fIndexLeadingGen);

	Double_t phiL = mcPartTmp->Phi();
	// Multiplicity transverse side
	Int_t multTS = 0;

	// Get multiplicity in transverse side
	Int_t pidCodeMC = 0;
	Double_t ipt = 0.;
	Double_t etaPart = -10.0;
	Double_t phiPart = -10.0;
	Bool_t isPhysPrim = kFALSE;
	Double_t qPart = 0;
	Int_t pPDG = -10;
	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		if(fIsPseudoRec){
			if(ipart == fIndexLeadingRec)continue;
		}
		else{ 
			if(ipart == fIndexLeadingGen)continue;
		}

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;
		//selection of primary charged particles
		if(!(mcPart->GetPDG())) continue;
		qPart = mcPart->GetPDG()->Charge()/3.;
		if(TMath::Abs(qPart)<0.001) continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		if(!isPhysPrim)
			continue;
		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);
		etaPart = mcPart -> Eta();
		if(TMath::Abs(etaPart)>0.8)continue;
		ipt = mcPart->Pt();
		if(fIsPseudoRec)
			if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),ipt)) continue;
		if(ipt<0.15)continue;
		phiPart = mcPart -> Phi();
		Double_t DPhi = DeltaPhi(phiL,phiPart);

		if(fIsPseudoRec){
			if(TMath::Abs(DPhi)<pi/3.0){
				fDphiNSRec->Fill(DPhi);
			}
			// away side
			else if(TMath::Abs(DPhi-pi)<pi/3.0){
				fDphiASRec->Fill(DPhi);
			}
			// transverse side
			else{
				multTS++;
				fDphiTSRec->Fill(DPhi);
			}
		}
		else{
			if(TMath::Abs(DPhi)<pi/3.0){
				fDphiNS->Fill(DPhi);
			}
			// away side
			else if(TMath::Abs(DPhi-pi)<pi/3.0){
				fDphiAS->Fill(DPhi);
			}
			// transverse side
			else{
				multTS++;
				fDphiTS->Fill(DPhi);
			}
		}


	}

	if(fIsPseudoRec)
		fMultTSRec->Fill(multTS);
	else
		fMultTS->Fill(multTS);


	// selecting topological regions
	pidCodeMC = 0;
	ipt = 0.;
	etaPart = -10.0;
	phiPart = -10.0;
	isPhysPrim = kFALSE;
	qPart = 0;
	pPDG = -10;

	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		if(fIsPseudoRec){
			if(ipart == fIndexLeadingRec)continue;
		}
		else{
			if(ipart == fIndexLeadingGen)continue;
		}


		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;

		if(!mcPart->GetPDG())continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		qPart = mcPart->GetPDG()->Charge()/3.;
		// only primary charged particles
		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);
		Bool_t isSelectedPart = kTRUE;
		for(Int_t i=0; i<11; ++i) 
			if( pidCodeMC == i ) 
				isSelectedPart = kFALSE;
		if ( isSelectedPart ) continue;
		ipt = mcPart->Pt();
		if(ipt<0.15)continue;
		phiPart = mcPart -> Phi();
		Double_t DPhi = DeltaPhi(phiL,phiPart);
		if(TMath::Abs(mcPart->Eta()) >= 0.8)continue;

		for(Int_t i=0; i<11; ++i)
		{
			if( pidCodeMC == i )
			{
				if( isPrimary[i] == kTRUE && isPhysPrim == kFALSE ) 
					continue;

				if(fIsPseudoRec){
					if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),ipt)) continue;
					if(TMath::Abs(DPhi)<pi/3.0){// near side
						fHistPtVsNchNSRec[i]->Fill(1.0*multTS,ipt);
					}
					else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
						fHistPtVsNchASRec[i]->Fill(1.0*multTS,ipt);
					}
					else{// transverse side
						fHistPtVsNchTSRec[i]->Fill(1.0*multTS,ipt);
					}


				}else{

					if(TMath::Abs(DPhi)<pi/3.0){// near side
						fHistPtVsNchNS[i]->Fill(1.0*multTS,ipt);
					}       
					else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
						fHistPtVsNchAS[i]->Fill(1.0*multTS,ipt);
					}       
					else{// transverse side
						fHistPtVsNchTS[i]->Fill(1.0*multTS,ipt);
					}       
				}

			}
		}

	} // particle loop

}

void AliAnalysisTaskGenUeSpherocity::MakeUeSoNch08Analysis(vector<Int_t> &mult){



	Int_t BinNchForSpherocity=mult[1];
	if(fspherocity_gen_ptWeighted<0)
		return;

	fSoWeighedVsNchPtL->Fill(1.0*BinNchForSpherocity,fspherocity_gen_ptWeighted);

	if(BinNchForSpherocity<3||BinNchForSpherocity>122)
		return;

	TParticle* mcPartTmp         = 0x0;
	mcPartTmp                    = (TParticle *)fMcEvent->Particle(fIndexLeadingGen);
	Double_t phiL = mcPartTmp->Phi();

	//Int_t pidCodeMC = 0;
	Double_t ipt = 0.;
	Double_t etaPart = -10.0;
	Double_t phiPart = -10.0;
	Bool_t isPhysPrim = kFALSE;
	Double_t qPart = 0;
	//Int_t pPDG = -10;

	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		if(ipart == fIndexLeadingGen)continue;

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;
		//selection of primary charged particles
		if(!(mcPart->GetPDG())) continue;
		qPart = mcPart->GetPDG()->Charge()/3.;
		if(TMath::Abs(qPart)<0.001) continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		if(!isPhysPrim)
			continue;
		//pPDG = TMath::Abs(mcPart->GetPdgCode());
		//pidCodeMC = GetPidCode(pPDG);
		etaPart = mcPart -> Eta();
		if(TMath::Abs(etaPart)>0.8)continue;
		ipt = mcPart->Pt();
		if(ipt<0.15)continue;
		phiPart = mcPart -> Phi();
		Double_t DPhi = DeltaPhi(phiL,phiPart);

		// inclusive
		hDphiSoIS[BinNchForSpherocity-3]->Fill(fspherocity_gen_ptWeighted,DPhi);
		hPtVsSoIS[BinNchForSpherocity-3]->Fill(fspherocity_gen_ptWeighted,ipt);
		// near side
		if(TMath::Abs(DPhi)<pi/3.0){
			hPtVsSoNS[BinNchForSpherocity-3]->Fill(fspherocity_gen_ptWeighted,ipt);
		}
		// away side
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			hPtVsSoAS[BinNchForSpherocity-3]->Fill(fspherocity_gen_ptWeighted,ipt);
		}
		// transverse side
		else{
			hPtVsSoTS[BinNchForSpherocity-3]->Fill(fspherocity_gen_ptWeighted,ipt);
		}

	} // particle loop

}


void AliAnalysisTaskGenUeSpherocity::ParticleSel(Bool_t fIsPseudoRec, const vector<Int_t> &mult){


	Int_t pidCodeMC = 0;
	Double_t ipt = 0.;

	Bool_t isPhysPrim = kFALSE;
	Double_t qPart = 0;
	Int_t pPDG = -10;
	Double_t y = -10;
	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;

		FillHisto("fHistEvt",1.5);
		if(!mcPart->GetPDG())continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		qPart = mcPart->GetPDG()->Charge()/3.;
		// only primary charged particles
		if( isPhysPrim ){
			if(TMath::Abs(qPart)>0.001)
				for(Int_t i=0;i<3;++i){
					if(fIsPseudoRec)
						fMultRec[i]->Fill(1.0*mult[i],mcPart->Eta());
					//FillHisto(Form("fMultRec_%s",estimators[i]),1.0*mult[i],mcPart->Eta());
					else
						fMult[i]->Fill(1.0*mult[i],mcPart->Eta());
					//	FillHisto(Form("fMult_%s",estimators[i]),1.0*mult[i],mcPart->Eta());
				}
		}

		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);
		Bool_t isSelectedPart = kTRUE;
		for(Int_t i=0; i<11; ++i) 
			if( pidCodeMC == i ) 
				isSelectedPart = kFALSE;
		if ( isSelectedPart ) continue;

		fHistEta->Fill(mcPart->Eta());

		if (!(TMath::Abs(mcPart->Energy()-mcPart->Pz())>0.)) continue;
		Double_t myY = (mcPart->Energy()+mcPart->Pz())/(mcPart->Energy()-mcPart->Pz());
		if( myY <= 0 ) continue;
		y = 0.5*TMath::Log(myY);
		//y = mcPart->Y(); 
		ipt = mcPart->Pt();


		for(Int_t i=0; i<11; ++i)
		{
			if( pidCodeMC == i && TMath::Abs(y) < fY)
			{
				if( isPrimary[i] == kTRUE && isPhysPrim == kFALSE ) 
					continue;

				if(fIsPseudoRec){
					if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),mcPart->Pt())) continue;
					if(i==0) FillHisto("fHistYRec",y);
					fHistPtRec[i]->Fill(ipt);
					if(fbinPerc0_rec>=0){
						fHistPtVsSoPtWRec0[i][fbinPerc0_rec]->Fill(fspherocity_rec_ptWeighted,ipt);
						fHistPtVsSoRec0[i][fbinPerc0_rec]->Fill(fspherocity_rec,ipt);
					}
					if(fbinPerc1_rec>=0){
						fHistPtVsSoPtWRec1[i][fbinPerc1_rec]->Fill(fspherocity_rec_ptWeighted,ipt);
						fHistPtVsSoRec1[i][fbinPerc1_rec]->Fill(fspherocity_rec,ipt);
					}
					if(fbinPerc2_rec>=0){
						fHistPtVsSoPtWRec2[i][fbinPerc2_rec]->Fill(fspherocity_rec_ptWeighted,ipt);
						fHistPtVsSoRec2[i][fbinPerc2_rec]->Fill(fspherocity_rec,ipt);
					}


				}else{
					if(i==0) FillHisto("fHistY",y);
					fHistPt[i]->Fill(ipt);
					if(fbinPerc0_gen>=0){
						fHistPtVsSoPtW0[i][fbinPerc0_gen]->Fill(fspherocity_gen_ptWeighted,ipt);
						fHistPtVsSo0[i][fbinPerc0_gen]->Fill(fspherocity_gen,ipt);
					}
					if(fbinPerc1_gen>=0){
						fHistPtVsSoPtW1[i][fbinPerc1_gen]->Fill(fspherocity_gen_ptWeighted,ipt);
						fHistPtVsSo1[i][fbinPerc1_gen]->Fill(fspherocity_gen,ipt);
					}
					if(fbinPerc2_gen>=0){
						fHistPtVsSoPtW2[i][fbinPerc2_gen]->Fill(fspherocity_gen_ptWeighted,ipt);
						fHistPtVsSo2[i][fbinPerc2_gen]->Fill(fspherocity_gen,ipt);

					}

				}

			}
		}

	} // particle loop

}

Float_t AliAnalysisTaskGenUeSpherocity::GetSpherocity(Int_t nch_so, const vector<Float_t> &pt, const vector<Float_t> &eta, const vector<Float_t> &phi, const Bool_t isPtWeighted ){


	Float_t spherocity = -10.0;
	Float_t pFull = 0;
	Float_t Spherocity = 2;

	//computing total pt
	Float_t sumapt = 0;
	if(isPtWeighted)
		for(Int_t i1 = 0; i1 < nch_so; ++i1){
			sumapt += pt[i1];
		}
	else
		sumapt = 1.0*nch_so;
	//Getting thrust
	for(Int_t i = 0; i < 360/(fSizeStep); ++i){
		Float_t numerador = 0;
		Float_t phiparam  = 0;
		Float_t nx = 0;
		Float_t ny = 0;
		phiparam=( (TMath::Pi()) * i * fSizeStep ) / 180; // parametrization of the angle
		nx = TMath::Cos(phiparam);            // x component of an unitary vector n
		ny = TMath::Sin(phiparam);            // y component of an unitary vector n
		for(Int_t i1 = 0; i1 < nch_so; ++i1){

			Float_t pxA = 0;
			Float_t pyA = 0;
			if(isPtWeighted){
				pxA = pt[i1] * TMath::Cos( phi[i1] );
				pyA = pt[i1] * TMath::Sin( phi[i1] );
			}
			else{
				pxA = 1.0 * TMath::Cos( phi[i1] );
				pyA = 1.0 * TMath::Sin( phi[i1] );
			}



			numerador += TMath::Abs( ny * pxA - nx * pyA );//product between p  proyection in XY plane and the unitary vector
		}
		pFull=TMath::Power( (numerador / sumapt),2 );
		if(pFull < Spherocity)//maximization of pFull
		{
			Spherocity = pFull;
		}
	}

	spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;


	return spherocity;

}


//______________________________________________________________________________

void AliAnalysisTaskGenUeSpherocity::Terminate(Option_t*){

	fListOfObjects = dynamic_cast<TList*> (GetOutputData(1));
	if (!fListOfObjects) { Printf("ERROR: Output list not available"); return; }

	return;
}

//_____________________________________________________________________________

Int_t AliAnalysisTaskGenUeSpherocity::GetPidCode(Int_t pdgCode) const  {

	Int_t pidCode = 999;

	switch (TMath::Abs(pdgCode)) {
		case 211:
			pidCode = 0; // pion
			break;
		case 321:
			pidCode = 1; // kaon
			break;
		case 2212:
			pidCode = 2; // proton
			break;
		case 310:
			pidCode = 3; // K0s
			break;
		case 3122:
			pidCode = 4; // Lambda
			break;
		case 3312:
			pidCode = 5; // Xi-
			break;
		case 3334:
			pidCode = 6; // Omega-
			break;
		case 333:
			pidCode = 7; // phi(1020)
			break;
		case 313:
			pidCode = 8; // K*(892)0
			break;
		case 323:
			pidCode = 9; // K*(892) +-
			break;
		case 3212:
			pidCode = 10; // Sigma 0
			break;    
		default:
			break;
	};

	return pidCode;
}

Bool_t AliAnalysisTaskGenUeSpherocity::IsGoodTrack(Int_t pid, Int_t pdgcode, Double_t pt){
	Double_t efficiency;
	if(pid==0)
		efficiency = ch_Eff->Eval(pt);
	else if(pid==1){
		if(pdgcode>0)
			efficiency = k_pos_Eff->Eval(pt);
		else
			efficiency = k_neg_Eff->Eval(pt);
	}
	else if(pid==2){
		if(pdgcode>0)
			efficiency = p_pos_Eff->Eval(pt);
		else
			efficiency = p_neg_Eff->Eval(pt);
	}
	else if(pid==3)
		efficiency = k0s_Eff->Eval(pt);
	else if(pid==4){
		if(pdgcode>0)
			efficiency = la_Eff->Eval(pt);
		else
			efficiency = labar_Eff->Eval(pt);
	}
	else if(pid==5)
		efficiency = xi_Eff->Eval(pt);
	else if(pid==7)
		efficiency = phi_Eff->Eval(pt);
	else if(pid==8)
		efficiency = k0star_Eff->Eval(pt);
	else
		efficiency = ch_Eff->Eval(pt);


	if(gRandom->Uniform(1.0)<efficiency)
		return kTRUE;
	else
		return kFALSE;
}
void AliAnalysisTaskGenUeSpherocity::MakeAnaGen(Int_t fNso_gen, vector<Int_t> &mult_estimators, vector<Float_t> &pt_so,  vector<Float_t> &eta_so, vector<Float_t> &phi_so){



	fspherocity_gen_ptWeighted=-0.5;
	fspherocity_gen=-0.5;

	if(fNso_gen>2){
		fspherocity_gen_ptWeighted = GetSpherocity(fNso_gen, pt_so, eta_so, phi_so, kTRUE);
		fspherocity_gen = GetSpherocity(fNso_gen, pt_so, eta_so, phi_so, kFALSE);
	}
	// ### Event and particle selection
	for(Int_t i=0;i<3;++i){
		fNch[i]->Fill(1.0*mult_estimators[i]);
	}
	if(fspherocity_gen>=0&&fspherocity_gen<=1){

		for(Int_t i=0;i<3;++i){
			fNchSoSel[i]->Fill(1.0*mult_estimators[i]);
			fSoVsNch[i]->Fill(1.0*mult_estimators[i],fspherocity_gen);
		}

	}
	if(fspherocity_gen_ptWeighted>=0&&fspherocity_gen_ptWeighted<=1){

		for(Int_t i=0;i<3;++i){
			fSoWeighedVsNch[i]->Fill(1.0*mult_estimators[i],fspherocity_gen_ptWeighted);
		}
	}

}
void AliAnalysisTaskGenUeSpherocity::MakeAnaRec(Int_t fNso_gen, vector<Int_t> &mult_estimators, vector<Float_t> &pt_so,  vector<Float_t> &eta_so, vector<Float_t> &phi_so){



	fspherocity_rec_ptWeighted=-0.5;
	fspherocity_rec=-0.5;

	if(fNso_gen>2){
		fspherocity_rec_ptWeighted = GetSpherocity(fNso_gen, pt_so, eta_so, phi_so, kTRUE);
		fspherocity_rec = GetSpherocity(fNso_gen, pt_so, eta_so, phi_so, kFALSE);
	}
	// ### Event and particle selection
	for(Int_t i=0;i<3;++i){
		fNchRec[i]->Fill(1.0*mult_estimators[i]);
	}
	if(fspherocity_rec>=0&&fspherocity_rec<=1){

		for(Int_t i=0;i<3;++i){
			fNchSoSelRec[i]->Fill(1.0*mult_estimators[i]);
			fSoVsNchRec[i]->Fill(1.0*mult_estimators[i],fspherocity_rec);
		}

	}
	if(fspherocity_rec_ptWeighted>=0&&fspherocity_rec_ptWeighted<=1){

		for(Int_t i=0;i<3;++i){
			fSoWeighedVsNchRec[i]->Fill(1.0*mult_estimators[i],fspherocity_rec_ptWeighted);
		}
	}

}
Int_t AliAnalysisTaskGenUeSpherocity::GetMultBin(Bool_t fIsPseudoRec,Int_t mult_int, Int_t mult_select){

	if(fIsPseudoRec){
		for(Int_t j=0;j<NchPercBin;++j){
			if(mult_select==0){
				if(mult_int>=nchBin_rec0[j] && mult_int<nchBin_rec0[j+1])
					return j;
			}
			if(mult_select==1){
				if(mult_int>=nchBin_rec1[j] && mult_int<nchBin_rec1[j+1])
					return j;
			}
			if(mult_select==2){
				if(mult_int>=nchBin_rec2[j] && mult_int<nchBin_rec2[j+1])
					return j;
			}
		}
	}
	else{
		for(Int_t j=0;j<NchPercBin;++j){
			if(mult_select==0){
				if(mult_int>=nchBin_gen0[j] && mult_int<nchBin_gen0[j+1])
					return j;
			}
			if(mult_select==1){
				if(mult_int>=nchBin_gen1[j] && mult_int<nchBin_gen1[j+1])
					return j;
			}
			if(mult_select==2){
				if(mult_int>=nchBin_gen2[j] && mult_int<nchBin_gen2[j+1])
					return j;
			}
		}
	}

	return -1;
}
Double_t AliAnalysisTaskGenUeSpherocity::DeltaPhi(Double_t phia, Double_t phib,
		Double_t rangeMin, Double_t rangeMax)
{
	Double_t dphi = -999;
	Double_t pi = TMath::Pi();

	if (phia < 0)         phia += 2*pi;
	else if (phia > 2*pi) phia -= 2*pi;
	if (phib < 0)         phib += 2*pi;
	else if (phib > 2*pi) phib -= 2*pi;
	dphi = phib - phia;
	if (dphi < rangeMin)      dphi += 2*pi;
	else if (dphi > rangeMax) dphi -= 2*pi;

	return dphi;
}
