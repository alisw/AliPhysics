/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

//==================================================================
// class for two-particle angular correlations analyses.
// by Beomkyu KIM, Junlee Kim.
//==================================================================

#include "TFile.h"
#include "TSystem.h"
#include "TGrid.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliCentrality.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliVMultiplicity.h"
#include "AliVVZERO.h"
#include "AliJetContainer.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <TVector3.h>
#include <TVectorT.h>
#include "AliJJet.h"
#include "AliAnalysisTaskRidge.h"
#include <AliDirList.h>
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalContainer.h"
#include "AliStack.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalJet.h"
#include <TClonesArray.h>
#include <TList.h>
#include <TProfile.h>
#include "AliAnalysisTaskEmcalJet.h"
#include "AliEmcalJet.h"
#include "AliAnalysisTaskRhoSparse.h"
#include "AliRhoParameter.h"

using namespace std;

const Double_t pi = TMath::Pi();

AliAnalysisTaskRidgeRunTable::AliAnalysisTaskRidgeRunTable() :
    fCollisionType(kUnknownCollType)
{;}

AliAnalysisTaskRidgeRunTable::AliAnalysisTaskRidgeRunTable(Int_t runnumber) 
{
    if (runnumber>=114737 && runnumber<=130850) fCollisionType = kPP; //LHC10bcde
    else if (runnumber>=144871 && runnumber<=146860) fCollisionType=kPP;//LHC11a
    else if (runnumber>=136851 && runnumber<=139517) fCollisionType=kAA;//LHC10h
    else if (runnumber>=167813 && runnumber<=170595) fCollisionType=kAA;//LHC11h
    else if (runnumber>=188356 && runnumber<=188503) fCollisionType=kPA;//LHC12g
    else if (runnumber>=189122 && runnumber<=192732) fCollisionType=kPA;//LHC12h
    else if (runnumber>=195344 && runnumber<=195483) fCollisionType=kPA;//LHC13b
    else if (runnumber>=195529 && runnumber<=195677) fCollisionType=kPA;//LHC13c
    else if (runnumber>=195724 && runnumber<=195872) fCollisionType=kPA;//LHC13d
    else if (runnumber>=195955 && runnumber<=195872) fCollisionType=kPA;//LHC13e
    else if (runnumber>=197669 && runnumber<=200000) fCollisionType=kPA;//LHC13g
    else if (runnumber>=256504 && runnumber<=260014) fCollisionType=kPP;//LHC16kl
    else fCollisionType=kPP;
}
AliAnalysisTaskRidgeRunTable::~AliAnalysisTaskRidgeRunTable()
{;}

//___________________________________________________________________
AliAnalysisTaskRidge::AliAnalysisTaskRidge()
:AliAnalysisTaskEmcalJet("AliAnalysisTaskRidge")
    , fOption() 
    , goodtrackindices()
	, fEMpool ()
	, fEMpooltracklet() 
	, fEMpoolMCALICE ()
	, fEMpoolMCCMS ()

{
}
//___________________________________________________________________
AliAnalysisTaskRidge::AliAnalysisTaskRidge
( 
      const char *name
    , const char *option
)
:AliAnalysisTaskEmcalJet(name)
    , fOption(option)
    , goodtrackindices()
	, fEMpool () 
	, fEMpooltracklet() 
	, fEMpoolMCALICE ()
	, fEMpoolMCCMS ()
{
    DefineOutput (1, AliDirList::Class());
}
//___________________________________________________________________
AliAnalysisTaskRidge::AliAnalysisTaskRidge
(
      const AliAnalysisTaskRidge& ap
)
    : fOption(ap.fOption)
    , goodtrackindices(ap.goodtrackindices)
	, fEMpool (ap.fEMpool)
	, fEMpooltracklet(ap.fEMpooltracklet) 
	, fEMpoolMCALICE(ap.fEMpoolMCALICE)
	, fEMpoolMCCMS(ap.fEMpoolMCCMS)
{
    DefineOutput (1, AliDirList::Class());
}
//___________________________________________________________________
AliAnalysisTaskRidge& AliAnalysisTaskRidge::operator = 
(
      const AliAnalysisTaskRidge& ap
)
{
    // assignment operator

    DefineOutput (1, AliDirList::Class());
    this->~AliAnalysisTaskRidge();
    new(this) AliAnalysisTaskRidge(ap);
    return *this;
}
//___________________________________________________________________
AliAnalysisTaskRidge::~AliAnalysisTaskRidge()
{
    delete fTrigger;
    delete fTrackCuts;
    delete fRunTable; 
    delete fOutput;
}

//___________________________________________________________________
void AliAnalysisTaskRidge::UserCreateOutputObjects()
{
	// Histograms container
	fOutput = new AliDirList();
	fOutput->SetOwner();

	// Offline triggers -----------------------------------------------------
	fTrigger = new AliTriggerAnalysis; // offline trigger
//	fTrigger -> SetFMDThreshold(0.3,0.5); // FMD threshold
	//-----------------------------------------------------------------------

	// TrackCuts for strangeness measure-------------------------------------
	fHistos = new THistManager("Ridgehists");

	Double1D varcentbinHigh = {
		   0, 0.001, 0.01, 0.02,
		0.05,   0.1,  0.5,    1};

        Double1D varcentbin = {0,1,2,5,10,20,50,100};
	Double1D varcentbinHeavy = {0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};

	binCent = AxisVar("Cent",varcentbinHigh);

	if( IsAA ) binCent = AxisVar("Cent",varcentbinHeavy);
	if( fOption.Contains("HighMult") ){ binCent = AxisVar("Cent",varcentbinHigh); }
	else if( !fOption.Contains("HighMult") ){ binCent = AxisVar("Cent",varcentbin); }

	Double1D ptbin = {0.1,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,14.0,100};

	binTPt = AxisVar("TPt",ptbin); //trig
	binAPt = AxisVar("APt",ptbin); //associate
	binNtrig = AxisFix("Ntrig",1,0.5,1.5);
	binUnipT = AxisFix("pt",200,0,20);

	binTrig = AxisFix("Trig",2,-0.5,1.5);
	binV0Amp = AxisFix("binV0Amp",3000,0,3e3);

	binPhi = AxisFix("phi",32,-0.5*pi-pi/32.0, 1.5*pi-pi/32.0);
	binEta = AxisFix("eta",40,-2.0,2.0);
        binMCEta = AxisFix("eta",80,-4.0,4.0);
	binTrkEff = AxisFix("Trkbin",5,0.5,5.5);

	Double1D pttrackbin = {
	0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
	0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,   1,  1.1, 1.2,
	 1.3, 1.4,  1.5, 1.6,  1.7, 1.8,  1.9,   2,  2.2, 2.4,
	 2.6, 2.8,    3, 3.2,  3.4, 3.6,  3.8,   4,  4.5,   5,
	20 };
//	 5.5,   6,  6.5,   7,   8 ,  10,   13,  20};
	binPt = AxisVar("Pt",pttrackbin);

        Double1D pttrackbin1 = {
        0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
        0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,   1,  1.1, 1.2,
         1.3, 1.4,  1.5, 1.6,  1.7, 1.8,  1.9,   2,  2.2, 2.4,
         2.6, 2.8,    3, 3.2,  3.4, 3.6,  3.8,   4,  4.5,   5,
	5.5,   6,  6.5,   7,   8 ,  10,   13,  20};	
	binPt1 = AxisVar("Pt",pttrackbin1);

	Double1D ltpttrackbin = {
	0.2, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 13.0, 20.0};
	Double1D jetptbin = {
	0, 10, 20, 30, 40, 50, 60, 80, 100, 1e5 };

	binRho = AxisFix("Rho",300,0,30);

	binLtpt = AxisVar("LPPt",ltpttrackbin);
	binJetpT = AxisVar("JetPt",jetptbin);

	Double1D verzbin = {-10,-8,-6,-4,-2,0,2,4,6,8,10};

	Double1D verzbinFine = {
	-10, -9, 
	 -8, -7, -6, -5, -4, -3, -2,
	  2,  3,  4,  5,  6,  7,  8,
	  9, 10 };

	binZ = AxisVar("Z",verzbin);
	if( fOption.Contains("FineBkg") ){
		binZ = AxisVar("Z",verzbinFine);
	}

        binPhiTrack = AxisFix("PHI",180,0,2*pi);
        binEtaTrack = AxisFix("ETA",80,-4,4);

	CreateTHnSparse("hRidgeLT","RidgeLT",6,{binCent,binPhi,binEta,binTPt,binAPt,binLtpt},"s");
	CreateTHnSparse("hRidgeMixingSLT","RidgeMixingSLT",6,{binCent,binPhi,binEta,binTPt,binAPt,binLtpt},"s");
	CreateTHnSparse("hNtrig","hNtrig",4,{binCent,binTPt,binNtrig,binLtpt},"s");

        CreateTHnSparse("hRidgeJet","RidgeJet",6,{binCent,binPhi,binEta,binTPt,binAPt,binJetpT},"s");
        CreateTHnSparse("hRidgeMixingSJet","RidgeMixingSJet",6,{binCent,binPhi,binEta,binTPt,binAPt,binJetpT},"s");
	CreateTHnSparse("hNtrigJet","hNtrigJet",4,{binCent,binTPt,binNtrig,binJetpT},"s");

	CreateTHnSparse("hTrackData","hTrackData",6,{binPt,binPhiTrack,binEtaTrack,binZ,binTrkEff,binCent},"s");
	CreateTHnSparse("hTrackDataCor","hTrackDataCor",5,{binPt,binPhiTrack,binEtaTrack,binCent,binZ},"s");

	CreateTHnSparse("hTrackDataLTRaw","hTrackDataLTRaw",5,{binCent,binPt,binPhiTrack,binEtaTrack,binZ},"s");

        CreateTHnSparse("nevtForMult","nevtForMult",2,{binCent,binV0Amp},"s");

	CreateTHnSparse("hRho","hRho",3,{binCent,binJetpT,binRho},"s");

	if( fOption.Contains("MC") ){
	        CreateTHnSparse("hRidgeMCALICELT","hRidgeMCALICELT",6,{binCent,binPhi,binMCEta,binTPt,binAPt,binLtpt},"s");
	        CreateTHnSparse("hRidgeMixingSMCALICELT","hRidgeMixingSMCALICELT",6,{binCent,binPhi,binMCEta,binTPt,binAPt,binLtpt},"s");
	        CreateTHnSparse("hNtrigMCALICE","hNtrigMCALICE",4,{binCent,binTPt,binNtrig,binLtpt},"s");

	        CreateTHnSparse("hRidgeMCCMSLT","hRidgeMCCMSLT",6,{binCent,binPhi,binMCEta,binTPt,binAPt,binLtpt},"s");
	        CreateTHnSparse("hRidgeMixingSMCCMSLT","hRidgeMixingSMCCMSLT",6,{binCent,binPhi,binMCEta,binTPt,binAPt,binLtpt},"s");
	        CreateTHnSparse("hNtrigMCCMS","hNtrigMCCMS",4,{binCent,binTPt,binNtrig,binLtpt},"s");

	        CreateTHnSparse("hTrackDataTrue","hTrackDataTrue",5,{binPt,binPhiTrack,binEtaTrack,binCent,binZ},"s");
	        CreateTHnSparse("hTrackMCallcut","hTrackMCallcut",5,{binPt,binPhiTrack,binEtaTrack,binZ,binTrkEff},"s");

	        CreateTHnSparse("hTrackMCCMS","hTrackMCCMS",5,{binPt,binPhiTrack,binEtaTrack,binCent,binZ},"s");
	        CreateTHnSparse("hTrackMCALICE","hTrackMCALICE",5,{binPt,binPhiTrack,binEtaTrack,binCent,binZ},"s");

	        CreateTHnSparse("hTrackMCCMSLT","hTrackMCCMSLT",5,{binCent,binPt,binPhiTrack,binEtaTrack,binZ},"s");
	        CreateTHnSparse("hTrackMCALICELT","hTrackMCALICELT",5,{binCent,binPt,binPhiTrack,binEtaTrack,binZ},"s");

	        CreateTHnSparse("TrigEffMult","TrigEffMult",2,{binCent,binTrig},"s");

	        CreateTHnSparse("MultiplicityStudy","MultiplicityStudy",4,
	                {binEta,binUnipT,binCent,binV0Amp},"s");
	}


	vector<TString> ent ={
	        "All","IsTriggered","IsNotPileup",
	        "IsValidVtx","IsGoodVtx","IsSelectedFromAliMultSelection",
	        "IsMultiplicityInsideBin" };
	
	auto h = fHistos->CreateTH1("hEventNumbers","",ent.size(), 0, ent.size());
	for(auto i=0u;i<ent.size();i++) h->GetXaxis()->SetBinLabel(i+1,ent.at(i).Data());

	fHistos->CreateTH1("hJetPt","",240,0,120);
        fHistos->CreateTH1("hJetEta","",100,-1.0,1.0);
        fHistos->CreateTH1("hJetPhi","",100,-4,4);

        fHistos->CreateTH1("hJetPtCor","",240,0,120);


        fHistos->CreateTH1("hLJetPt","",240,0,120);
        fHistos->CreateTH1("hLJetEta","",100,-1.0,1.0);
        fHistos->CreateTH1("hLJetPhi","",100,-4,4);

	fHistos->CreateTH1("hLHPt","",240,0,120);

	fHistos->CreateTH1("hPtCons","",200,0,5);

	fHistos->CreateTH2("hLHPt_JetpT","",240,0,120,240,0,120);

	fHistos->CreateTH1("hHMT","",1000,0,1,"s");
	fHistos->CreateTH1("hMB","",100,0,100,"s");

	fHistos->CreateTH2("hMB_V0M","",100,0,100,1000,0,3000,"s");
	fHistos->CreateTH2("hHMT_V0M","",1000,0,1,1000,0,3000,"s");

        fHistos->CreateTH1("hZvtx","",620,-15.5,15.5,"s");

	fHistos->CreateTH2("hPhiEta","",180,0,2*pi,40,-2,2);
        fHistos->CreateTH2("hPhiEtaCor","",180,0,2*pi,40,-2,2);

	fEMpool.resize(binCent.GetNbins(),vector<eventpool> (binZ.GetNbins()));
	fEMpooltracklet.resize(binCent.GetNbins(),vector<eventpooltracklet> (binZ.GetNbins()));
	fEMpoolMCALICE.resize(binCent.GetNbins(),vector<eventpoolMC> (binZ.GetNbins()));
	fEMpoolMCCMS.resize(binCent.GetNbins(),vector<eventpoolMC> (binZ.GetNbins()));	

	fOutput -> Add( fHistos->GetListOfHistograms() );
	PostData(1, fOutput);

	if( !fOption.Contains("ITS") ){
		for(int i=0;i<fEff_npT_step;i++){
			std::vector<double> elem;
			elem.resize(fEff_neta_step);
			Eff.push_back(elem);
		}
	}
	else if( fOption.Contains("ITS") ){
                for(int i=0;i<fEff_npT_step;i++){
                        std::vector<double> elem;
                        elem.resize(ITS_fEff_neta_step);
                        Eff.push_back(elem);
                }
	}

	if( fOption.Contains("Add3DEff") ){
		for(int i=0;i<fEff_npT_step;i++){
			std::vector< std::vector<double> > elem2D;
			for(int j=0;j<fEff_neta_step;j++){
				std::vector<double> elem;
				elem.resize(fEff_nphi_step);
				elem2D.push_back(elem);
			}
			Eff3D.push_back(elem2D);
		}
	}

	if( fOption.Contains("GRID") ){
		TGrid::Connect("alien://");
		fefficiencyFile = TFile::Open("alien:///alice/cern.ch/user/j/junlee/Efficiency_RIDGE/EffOut.root","read");
		if( fOption.Contains("Add3DEff") )fefficiency3DFile = TFile::Open("alien:///alice/cern.ch/user/j/junlee/Efficiency_RIDGE/Eff3DOut.root","read");
	}


	V0M_mean = 120;

}
//___________________________________________________________________
void AliAnalysisTaskRidge::Exec(Option_t* )
{

	// Pointer to a event----------------------------------------------------
	AliVEvent *event = InputEvent();
	if( !event ){ Printf("ERROR: Could not retrieve event");  return; }
	// ----------------------------------------------------------------------

	// connect to ESD tree --------------------------------------------------
	Int_t runnumber;
//	TString foption = option;
	TString Period;

	event->IsA()==AliESDEvent::Class() 
		? fEvt = dynamic_cast<AliESDEvent*>(event)
		: fEvt = dynamic_cast<AliAODEvent*>(event);
	if( !fEvt ) return;

	if( fOption.Contains("MC") ){ IsMC = kTRUE; }



	if( IsFirstEvent ){
		runnumber = fEvt->GetRunNumber();
        	fRunTable = new AliAnalysisTaskRidgeRunTable(runnumber);

//		if( !fefficiencyFile ) return;

		     if( runnumber >= 252235 && runnumber <= 252330 ) Period = "LHC16d";
		else if( runnumber >= 253437 && runnumber <= 253591 ) Period = "LHC16e";
		else if( runnumber >= 253659 && runnumber <= 253978 ) Period = "LHC16f";
		else if( runnumber >= 254128 && runnumber <= 254332 ) Period = "LHC16g";
		else if( runnumber >= 254604 && runnumber <= 255467 ) Period = "LHC16h";
		else if( runnumber >= 255539 && runnumber <= 255618 ) Period = "LHC16i";
		else if( runnumber >= 256219 && runnumber <= 256418 ) Period = "LHC16j";
		else if( runnumber >= 256941 && runnumber <= 256219 ) Period = "LHC16k";
		else if( runnumber >= 258962 && runnumber <= 259888 ) { Period = "LHC16l"; V0M_mean=89.9003; }
		else if( runnumber >= 262424 && runnumber <= 264035 ) { Period = "LHC16o"; V0M_mean=86.3912; }
		else if( runnumber >= 264076 && runnumber <= 264347 ) { Period = "LHC16p"; V0M_mean=138.814; }

		else if( runnumber >= 270581 && runnumber <= 270667 ) Period = "LHC17c";
		else if( runnumber >= 270822 && runnumber <= 270830 ) Period = "LHC17e";
		else if( runnumber >= 270854 && runnumber <= 270865 ) Period = "LHC17f";
		else if( runnumber >= 270882 && runnumber <= 271777 ) Period = "LHC17g";
		else if( runnumber >= 271870 && runnumber <= 273103 ) { Period = "LHC17h"; V0M_mean=127.895; }
		else if( runnumber >= 273591 && runnumber <= 274442 ) { Period = "LHC17i"; V0M_mean=124.276; }
		else if( runnumber >= 274593 && runnumber <= 274671 ) Period = "LHC17j";
		else if( runnumber >= 274690 && runnumber <= 276508 ) { Period = "LHC17k"; V0M_mean=121.31; }
		else if( runnumber >= 276551 && runnumber <= 278216 ) { Period = "LHC17l"; V0M_mean=119.144; }
		else if( runnumber >= 278914 && runnumber <= 280140 ) { Period = "LHC17m"; V0M_mean=117.165; }
		else if( runnumber >= 280282 && runnumber <= 281961 ) { Period = "LHC17o"; V0M_mean=113.45; }
		else if( runnumber >= 282528 && runnumber <= 282704 ) { Period = "LHC17r"; V0M_mean=111.462; }

		else if( runnumber >= 285009 && runnumber <= 285396 ) Period = "LHC18b";
		else if( runnumber >= 285978 && runnumber <= 286350 ) { Period = "LHC18d"; V0M_mean=131.868; }
		else if( runnumber >= 286380 && runnumber <= 286937 ) { Period = "LHC18e"; V0M_mean=131.397; }
		else if( runnumber >= 287000 && runnumber <= 287658 ) { Period = "LHC18f"; V0M_mean=130.591; }
		else if( runnumber >= 288750 && runnumber <= 288619 ) Period = "LHC18g";
		else if( runnumber >= 288806 && runnumber <= 288804 ) { Period = "LHC18h"; V0M_mean=130.86; }
		else if( runnumber >= 288909 && runnumber <= 288861 ) Period = "LHC18i";
		else if( runnumber >= 288943 && runnumber <= 288943 ) { Period = "LHC18j"; V0M_mean=131.17; }
		else if( runnumber >= 289240 && runnumber <= 289971 ) { Period = "LHC18l"; V0M_mean=131.59; }
		else if( runnumber >= 290323 && runnumber <= 292839 ) { Period = "LHC18m"; V0M_mean=130.467; }
		else if( runnumber >= 293359 && runnumber <= 293357 ) Period = "LHC18n";
		else if( runnumber >= 289201 && runnumber <= 289165 ) { Period = "LHC18k"; V0M_mean=127.642; }
		else if( runnumber >= 293475 && runnumber <= 293898 ) { Period = "LHC18o"; V0M_mean=124.973; }
		else if( runnumber >= 294009 && runnumber <= 294925 ) Period = "LHC18p";


		TH2D* hEfficiencyHist;
		TH3D* hEfficiency3DHist;

		if( !fefficiencyFile ){
                	for(int i=0;i<fEff_npT_step;i++){
                	        for(int j=0;j<fEff_neta_step;j++){
                	                Eff[i][j] = 0.5;
                	        }
                	}
		}

		if( fOption.Contains("Add3DEff") && !fefficiency3DFile ){
			for(int i=0;i<fEff_npT_step;i++){
				for(int j=0;j<fEff_neta_step;j++){
					for(int k=0;k<fEff_nphi_step;k++){
						Eff3D[i][j][k] = 1.0;
					}
				}
			}
		}
        	if( fefficiencyFile ){
//			cout << (bool)fefficiencyFile->FindObject(Form("%s_Hyb8cm",Period.Data())) << endl;
			hEfficiencyHist = (TH2D*)fefficiencyFile->Get(Form("%s_Hyb8cm",Period.Data()));
        	        if( fOption.Contains("Glb") ){ hEfficiencyHist = (TH2D*)fefficiencyFile->Get(Form("%s_Glb8cm",Period.Data())); }
        	        if( fOption.Contains("SDD") ){ hEfficiencyHist = (TH2D*)fefficiencyFile->Get(Form("%s_GlbSDD8cm",Period.Data())); }
        	        if( fOption.Contains("TightVtx") ){ hEfficiencyHist = (TH2D*)fefficiencyFile->Get(Form("%s_Hyb6cm",Period.Data())); }
	
	                if( !hEfficiencyHist ){ hEfficiencyHist = (TH2D*)fefficiencyFile->Get("LHC16l_Hyb8cm"); }

			if( fOption.Contains("MC") ){ hEfficiencyHist = (TH2D*)fefficiencyFile->Get("LHC16l_Hyb8cm"); }

			if( hEfficiencyHist ){
	                	for(int i=0;i<fEff_npT_step;i++){
	                	        for(int j=0;j<fEff_neta_step;j++){
	                	                if( i<hEfficiencyHist->GetNbinsY() ) Eff[i][j] = hEfficiencyHist->GetBinContent(j+1,i+1);
	                	                else{ Eff[i][j] = hEfficiencyHist->GetBinContent(j+1, hEfficiencyHist->GetNbinsY() ); }	
	                	                if( Eff[i][j] < 0.01 ){ Eff[i][j] = 1.0; }
	                	        }
	                	}
			}
			else if( !hEfficiencyHist ){
				for(int i=0;i<fEff_npT_step;i++){
					for(int j=0;j<fEff_neta_step;j++){
						Eff[i][j] = 0.25;
					}
				}
			}
	        }
		if( fOption.Contains("Add3DEff") && fefficiency3DFile ){
			hEfficiency3DHist = (TH3D*)fefficiency3DFile->Get(Form("%s_Hyb8cm",Period.Data()));
			if( fOption.Contains("Glb") ) hEfficiency3DHist = (TH3D*)fefficiency3DFile->Get(Form("%s_Glb8cm",Period.Data()));
                        if( fOption.Contains("SDD") ) hEfficiency3DHist = (TH3D*)fefficiency3DFile->Get(Form("%s_GlbSDD8cm",Period.Data()));
                        if( fOption.Contains("TightVtx") ) hEfficiency3DHist = (TH3D*)fefficiency3DFile->Get(Form("%s_Hyb6cm",Period.Data()));

                        if( !hEfficiency3DHist ){ hEfficiency3DHist = (TH3D*)fefficiency3DFile->Get("LHC16l_Hyb8cm"); }

                        if( fOption.Contains("MC") ){ hEfficiency3DHist = (TH3D*)fefficiency3DFile->Get("LHC16l_Hyb8cm"); }

                        if( hEfficiency3DHist ){
                                for(int i=0;i<fEff_npT_step;i++){
                                        for(int j=0;j<fEff_neta_step;j++){
						for(int k=0;k<fEff_nphi_step;k++){
                                                	if( i<hEfficiency3DHist->GetNbinsZ() ) Eff3D[i][j][k] = hEfficiency3DHist->GetBinContent(k+1,j+1,i+1);
                                                	else{ Eff3D[i][j][k] = hEfficiency3DHist->GetBinContent(k+1,j+1, hEfficiency3DHist->GetNbinsZ() ); }
                                                	if( Eff3D[i][j][k] < 0.01 ){ Eff3D[i][j][k] = 1.0; }
						}
                                        }
                                }
                        }
                        else if( !hEfficiency3DHist ){
                                for(int i=0;i<fEff_npT_step;i++){
                                        for(int j=0;j<fEff_neta_step;j++){
						for(int k=0;k<fEff_nphi_step;k++){
                                                	Eff3D[i][j][k] = 1.0;
						}
                                        }
                                }
                        }	
		}
		IsFirstEvent = kFALSE;
        }

	fCent = 200;
	fZ = 0.0;

	AliInputEventHandler* inputHandler;
	inputHandler = (AliInputEventHandler*)
		AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();

	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

	fHistos -> FillTH1("hEventNumbers","All",1);

	double v0amplitude=0;

	double JetConeSize = 0.4;
	if( fOption.Contains("SmallCone") ) JetConeSize = 0.2;
	double JetEtaAccpetance = 0.8 - JetConeSize;


	if( !fOption.Contains("HighMult") && !fOption.Contains("SmallCone") )		fJetTask = (AliJJetTask*)(AliAnalysisManager::GetAnalysisManager()->GetTask( "AliJJetTask" ));
	else if( fOption.Contains("HighMult") && !fOption.Contains("SmallCone") )	fJetTask = (AliJJetTask*)(AliAnalysisManager::GetAnalysisManager()->GetTask( "AliJJetTaskHighMult" ));
	else if( !fOption.Contains("HighMult") && fOption.Contains("SmallCone") )	fJetTask = (AliJJetTask*)(AliAnalysisManager::GetAnalysisManager()->GetTask( "AliJJetTaskSmallCone" ));
	else if( fOption.Contains("HighMult") && fOption.Contains("SmallCone") )	fJetTask = (AliJJetTask*)(AliAnalysisManager::GetAnalysisManager()->GetTask( "AliJJetTaskSmallConeHighMult" ));

	sel = (AliMultSelection*) fEvt -> FindListObject("MultSelection");
	if( sel ){ 
		fCent = sel->GetMultiplicityPercentile("V0M");
		AliVVZERO* lVV0 = fEvt->GetVZEROData();
        	for(int i=0;i<64;i++){ v0amplitude += lVV0->GetMultiplicity(i); }
	}
	TObjArray* fjets = (TObjArray*)fJetTask->GetAliJJetList(1);

	auto jettest0 = fJetTask->GetJetContainer(0);
	auto jettest1 = fJetTask->GetJetContainer(1);
	auto ktjets = fJetTask->GetJetContainer(2);

	RHO = 0.0; RHOM = 0.0;

	this->RhoSparse(ktjets, jettest1, 2);

	AliJJet* ktJet;
	AliJBaseTrack* kttrack;
	AliAODTrack* ktt_aod;
	
	TLorentzVector* LorentzTrack;

	double rho = RHO;
	double area = 0.0;
	AliJJet *Ljet = dynamic_cast<AliJJet*>( fjets->At(0) );

	fJetPt = 0.0;
	double JetEta = -10.0;
	double JetPhi = -10.0;


	if( Ljet ){
		fJetPt = Ljet->Pt();
		JetEta = Ljet->Eta();
		JetPhi = Ljet->Phi();
		area = Ljet->GetArea();
	}

	AliJJet* Cjet;
	for(int i=0;i<fjets->GetEntries();i++){
		Cjet = dynamic_cast<AliJJet*>( fjets->At(i) );
		if( fabs( Cjet->Eta() ) > JetEtaAccpetance ){ continue; }
		if( ( fJetPt < Cjet->Pt() ) ){
			fJetPt = Cjet->Pt();
			JetEta = Cjet->Eta();
			JetPhi = Cjet->Phi();
			area = Cjet->GetArea();
		}
	}
        const AliVVertex* trackVtx = fEvt->GetPrimaryVertexTPC() ;
        const AliVVertex* spdVtx   = fEvt->GetPrimaryVertexSPD() ;

        if( IsMC ){
                AliAODMCHeader *cHeaderAOD  = dynamic_cast<AliAODMCHeader*>
                        (fEvt->FindListObject(AliAODMCHeader::StdBranchName()));
                if( cHeaderAOD ) fZ_gen = cHeaderAOD -> GetVtxZ();
        }

	AbsZmax = 8.0;
	if( fOption.Contains("TightVtx") ) AbsZmax = 6.0;
	
	Bool_t IsTriggered = kFALSE;
	Bool_t IsNotPileup = kFALSE;
	Bool_t IsValidVtx = kFALSE;
	Bool_t IsGoodVtx = kFALSE;
	Bool_t IsSelectedFromAliMultSelection = kFALSE;
	Bool_t IsMultiplicityInsideBin = kFALSE;

	Bool_t IsSelectedFromAliMultSelectionForSysZ = kFALSE;
//IsTriggered*************************************
	if(fRunTable->IsAA() || fRunTable->IsPA()){
                IsTriggered = (inputHandler -> IsEventSelected()) & (AliVEvent::kINT7);
        }
        else if(fRunTable->IsPP() ){
                IsTriggered = (inputHandler -> IsEventSelected()) & AliVEvent::kINT7;
                if( fOption.Contains("HighMult") ){
                        IsTriggered = (inputHandler -> IsEventSelected()) & AliVEvent::kHighMultV0;
                }
	}
	if( IsMC ){
                IsTriggered = (inputHandler -> IsEventSelected()) & AliVEvent::kINT7;
//                if( fOption.Contains("HighMult") ){
//                        IsTriggered = inputHandler -> IsEventSelected() & AliVEvent::kHighMultV0;
//                }
        }        
//***********************************************



//IsNotPileup***********************************
	if( IsMC ) IsNotPileup = kTRUE;
	else if( !fRunTable->IsPP() ) IsNotPileup = kTRUE;
	else if( !IsMC && fRunTable->IsPP() && !event->IsPileupFromSPDInMultBins() ) IsNotPileup = kTRUE;
	if( fOption.Contains("PileupTest") ){ IsNotPileup = kTRUE; }
	if( fOption.Contains("PileupMV") ){ IsNotPileup = sel->GetThisEventIsNotPileupMV(); }

//*********************************************


//IsGoodVtx************************************
	zbin = -1;
	if( spdVtx ){
	        if( spdVtx->GetNContributors() > 0.5 ) IsValidVtx = kTRUE;
	        fZ = spdVtx->GetZ();
	        zbin = binZ.FindBin(fZ) -1;
	        if( fabs(fZ) < AbsZmax && !(zbin < 0 )) IsGoodVtx = kTRUE;
	}
//********************************************



//IsSelectedFromAliMultSelection**************
	if( !IsMC && fRunTable->IsPP() && !fOption.Contains("HighMult") ){
		if( sel->IsEventSelected() ) IsSelectedFromAliMultSelection = kTRUE;

		if( sel->GetThisEventIsNotPileupInMultBins() &&
	        sel->GetThisEventINELgtZERO() &&
	        sel->GetThisEventPassesTrackletVsCluster() &&
	        sel->GetThisEventHasNoInconsistentVertices() &&
	        fabs(fZ) < 8.0 ){
	                IsSelectedFromAliMultSelectionForSysZ = kTRUE;
	        }
	}
	else if( fOption.Contains("HighMult") ){
	        if( sel->GetThisEventIsNotPileup() &&
	        sel->GetThisEventIsNotPileupInMultBins() &&
	        sel->GetThisEventHasNoInconsistentVertices() &&
	        sel->GetThisEventPassesTrackletVsCluster() ){
	                IsSelectedFromAliMultSelection = kTRUE;
	                if( fabs(fZ) < 8.0 ){
	                        IsSelectedFromAliMultSelectionForSysZ = kTRUE;
	                }
	        }
		if( fOption.Contains("PileupTest") || fOption.Contains("PileupMV") ){
			if( sel->GetThisEventHasNoInconsistentVertices() &&
			sel->GetThisEventPassesTrackletVsCluster() ){
				IsSelectedFromAliMultSelection = kTRUE;
			}
		}
	}
	else if( IsMC ) IsSelectedFromAliMultSelection = kTRUE;
//*******************************************




//IsMultiplicityInsideBin Flag Configuration********
	centbin = binCent.FindBin(fCent) -1;
	if( centbin >= 0 && centbin < binCent.GetNbins() ) IsMultiplicityInsideBin = kTRUE;
	if( fOption.Contains("HighMult") && fOption.Contains("CUTwithV0M") ){
		if( v0amplitude < 5.0 * V0M_mean || v0amplitude > 9.0 * V0M_mean ){
			IsMultiplicityInsideBin = kFALSE;
		}
	}
//******************************************
	
	if( IsTriggered ) fHistos->FillTH1("hEventNumbers","IsTriggered",1);
	if( IsTriggered && IsNotPileup ) fHistos->FillTH1("hEventNumbers","IsNotPileup",1);
	if( IsTriggered && IsNotPileup && IsValidVtx ) fHistos->FillTH1("hEventNumbers","IsValidVtx",1);
	if( IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx ) fHistos->FillTH1("hEventNumbers","IsGoodVtx",1);
	if( IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx && IsSelectedFromAliMultSelection ) fHistos->FillTH1("hEventNumbers","IsSelectedFromAliMultSelection",1);
	if( IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin ){
	        fHistos->FillTH1("hEventNumbers","IsMultiplicityInsideBin",1);

		if( fJetPt>0.1 ){
			fHistos->FillTH1("hLJetPt",fJetPt,1.0);
			fHistos->FillTH1("hLJetEta",JetEta,1.0);
			fHistos->FillTH1("hLJetPhi",JetPhi,1.0);
		}

        	for(int i=0;i<fjets->GetEntries();i++){
        	        Cjet = dynamic_cast<AliJJet*>( fjets->At(i) );
        	        if( fabs( Cjet->Eta() ) > JetEtaAccpetance ){ continue; }
/*
        	        if( ( fJetPt <  Cjet->Pt() ) ){
        	                fJetPt = Cjet->Pt();
        	                JetEta = Cjet->Eta();
        	                JetPhi = Cjet->Phi();
        	        }
*/
                        fHistos->FillTH1("hJetPt",Cjet->Pt(),1.0);
                        fHistos->FillTH1("hJetEta",Cjet->Eta(),1.0);
                        fHistos->FillTH1("hJetPhi",Cjet->Phi(),1.0);

			FillTHnSparse("hRho",{fCent,Cjet->Pt(),rho},1.0);

			fHistos->FillTH1("hJetPtCor",Cjet->Pt() - rho*Cjet->GetArea(),1.0);
        	}

		fJetPt -= rho*area;

	        if( !fOption.Contains("HighMult") ){
	                fHistos->FillTH1("hMB",fCent,1);
	                fHistos->FillTH2("hMB_V0M",fCent,v0amplitude,1);
	        }
	        else if( fOption.Contains("HighMult") ){
	                fHistos->FillTH1("hHMT",fCent,1);
	                fHistos->FillTH2("hHMT_V0M",fCent,v0amplitude,1);
	        }
	        fHistos->FillTH1("hZvtx",fZ,1);
	}

	if( fabs( JetEta ) > JetEtaAccpetance ){
		fJetPt = 0.0;
	}

	if( IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin && fOption.Contains("EvtSelStudy") && fOption.Contains("jtptstudy") ){
		fsetmixing = kFALSE; this -> GoodTracksSelection( 2 );
	}

//	if( !fOption.Contains("EvtSelStudy") && IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin ){
	if( !fOption.Contains("EvtSelStudy") && IsTriggered && IsNotPileup && IsValidVtx && fabs(fZ) < 10.0 && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin ){
		if( IsMC && fEvt->IsA()==AliAODEvent::Class() ){
			fMCArray = (TClonesArray*) fEvt->FindListObject("mcparticles");
			if (fabs(fZ_gen)<10.0){
	        	        Int_t nTracksMC = fMCArray->GetEntries();
	        	        for(Int_t iTracks = 0; iTracks < nTracksMC; iTracks++){
	        	                AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
	        	                if( !trackMC ) continue;
	        	                Int_t pdgCode = trackMC->PdgCode();
	        	                if( !(trackMC->IsPhysicalPrimary()) ) continue;
					if( trackMC->Charge() == 0 ) continue;

	        	                FillTHnSparse("hTrackMCallcut",{trackMC->Pt(),trackMC->Phi(),trackMC->Eta(),fZ_gen,1.0},1.0);
					if( fabs(fZ) < 8 ) FillTHnSparse("hTrackMCallcut",{trackMC->Pt(),trackMC->Phi(),trackMC->Eta(),fZ_gen,2.0},1.0);
					if( fabs(fZ) < 6 ) FillTHnSparse("hTrackMCallcut",{trackMC->Pt(),trackMC->Phi(),trackMC->Eta(),fZ_gen,3.0},1.0);
	        	        }
	        	}
		}
	}


	if( !fOption.Contains("EvtSelStudy") && IsTriggered && IsNotPileup && IsValidVtx && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin ){
		if( !fOption.Contains("EffCorrection") && IsGoodVtx ){
			if( fOption.Contains("Glb") && !fOption.Contains("SDD") ){ 
								if( this -> GoodTracksSelection( 4 ) ){ this -> FillTracks(); } }
			else if( fOption.Contains("GlbSDD") ){ 	if( this -> GoodTracksSelection( 5 ) ){ this -> FillTracks(); } }
			else{ 					if( this -> GoodTracksSelection( 2 ) ){ this -> FillTracks(); } }
//			else{  this -> GoodTracksSelection( 2 ); }
			if( fOption.Contains("MC") ){ if( this -> GoodTracksSelectionMC() ){ this -> FillTracksMC(); } }

		}
		else if( fOption.Contains("EffCorrection") ){
			fsetmixing = kFALSE;
			if( fabs(fZ) < 10 ){ this -> GoodTracksSelection( 1 ); }
			if( fabs(fZ) < 8  ){ this -> GoodTracksSelection( 2 ); }
			if( fabs(fZ) < 6  ){ this -> GoodTracksSelection( 3 ); }
			if( IsGoodVtx ){ this -> GoodTracksSelection( 4 ); }
			if( IsGoodVtx ){ this -> GoodTracksSelection( 5 ); }
		}
	}	

//Trig Efficiency ************** kIINT7 to INEL>0
	bool IsINEL=false;
	if( IsMC ){
        	if( fEvt->IsA()==AliAODEvent::Class() ){
        	        fMCArray = (TClonesArray*) fEvt->FindListObject("mcparticles");
        	        AliAODMCHeader *cHeaderAOD  = dynamic_cast<AliAODMCHeader*>
        	                (fEvt->FindListObject(AliAODMCHeader::StdBranchName()));
        	        const Int_t nTracksMC = fMCArray->GetEntriesFast();
        	        for(Int_t iTracks = 0; iTracks < nTracksMC; iTracks++){
        	                AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
        	                if( !trackMC ) continue;
        	                if( !(trackMC->IsPhysicalPrimary()) ) continue;
        	                if( trackMC->Charge() == 0 ) continue;
        	                if( fabs( trackMC->Eta() ) > 1.0 ) continue;
        	                IsINEL = true;
        	        }
        	}
	}	

	if( IsINEL ){
	        if( (inputHandler -> IsEventSelected()) & (AliVEvent::kINT7) ){
	                FillTHnSparse("TrigEffMult",{fCent,1.0},1.0 );
	        }
	        else{
	                FillTHnSparse("TrigEffMult",{fCent,0.0},1.0 );
	        }
	}
//*******************************************


//******
//        if( !fOption.Contains("EvtSelStudy") && IsNotPileup && IsValidVtx && fabs(fZ) < 10 && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin && IsINEL ){
	if( !fOption.Contains("EvtSelStudy") && IsNotPileup && IsValidVtx && fabs(fZ) < 10 && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin && IsINEL && IsTriggered ){
		if( IsMC && fEvt->IsA()==AliAODEvent::Class() ){
                        fMCArray = (TClonesArray*) fEvt->FindListObject("mcparticles");
                        if (fabs(fZ_gen)<10.0){
                                Int_t nTracksMC = fMCArray->GetEntries();
                                FillTHnSparse("nevtForMult",{fCent,v0amplitude},1.0);
                                for(Int_t iTracks = 0; iTracks < nTracksMC; iTracks++){
                                        AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
                                        if( !trackMC ) continue;
                                        Int_t pdgCode = trackMC->PdgCode();
                                        if( !(trackMC->IsPhysicalPrimary()) ) continue;
                                        if( trackMC->Charge() == 0 ) continue;

                                        FillTHnSparse("MultiplicityStudy",{trackMC->Eta(),trackMC->Pt(),fCent,v0amplitude},1.0);

				}
			}
		}
	}
//****

	PostData(1, fOutput);
}

Bool_t AliAnalysisTaskRidge::GoodTracksSelection(int trk){

	fFilterBit = 0x300;
        if( trk == 4 ) fFilterBit = 0x20;
        else if( trk == 5 ) fFilterBit = 0x60;

	const UInt_t ntracks = fEvt ->GetNumberOfTracks();
	goodtrackindices.clear();
	
	AliVTrack * track;

	tracklist *etl;
	eventpool *ep;

	//Event mixing pool
	if (fsetmixing && centbin>=0 && zbin>=0 &&
	centbin < binCent.GetNbins() && zbin < binZ.GetNbins() ){
		ep = &fEMpool[centbin][zbin];
		ep -> push_back( tracklist() ); //
//		ep -> push_back( etl );
		etl = &(ep->back());
	}

	fNTracks = 0;
	double LHPt = 0;
	for (UInt_t it = 0; it<ntracks; it++){
		if (fEvt->IsA()==AliESDEvent::Class()){
			track = (AliESDtrack*) fEvt ->GetTrack(it);
			if (!track) continue;
			if (!fTrackCuts->AcceptTrack((AliESDtrack*) track)) continue;
			fHistos->FillTH2("hPhiEta",track->Phi(),track->Eta());
		}
		else {
			track = (AliAODTrack*) fEvt ->GetTrack(it);
			if (!track) continue;

			if( trk < 3.5 ){ if( ! ((AliAODTrack*) track)->TestFilterBit(fFilterBit)) continue; } //for hybrid
			else if( trk > 3.5 ){ if( ! ((AliAODTrack*) track)->TestFilterMask(fFilterBit)) continue; } //for global


                        if( IsMC && fabs( dynamic_cast<AliAODMCHeader*>(fEvt->FindListObject(AliAODMCHeader::StdBranchName()))->GetVtxZ() ) < 10 ){
                                if( track->GetLabel()>-1 && dynamic_cast<AliAODMCParticle*>(fMCArray->At( track->GetLabel() ))->IsPhysicalPrimary() ){
                                        FillTHnSparse("hTrackDataTrue",{track->Pt(),track->Phi(),track->Eta(),fCent,fZ},1.0);
                                }
			}

			if( track->Pt()<fptcut ) continue;
			if( !fOption.Contains("ITS") && fabs(track->Eta())>fetacut ) continue;
			else if( fOption.Contains("ITS") && fabs(track->Eta())>1.3 ) continue;			
			FillTHnSparse("hTrackData",{track->Pt(),track->Phi(),track->Eta(),fZ,(double)trk,fCent},1.0);

			if( LHPt < track->Pt() ){ LHPt = track->Pt(); }

			fHistos->FillTH2("hPhiEta",track->Phi(),track->Eta(),1.0);
			if( !fOption.Contains("ITS") ){
				if( fOption.Contains("AddEff") ){
					if( track->Pt() > fEff_pT_max ){
						fHistos->FillTH2("hPhiEtaCor",track->Phi(),track->Eta(),1.0);
						FillTHnSparse("hTrackDataCor",{track->Pt(),track->Phi(),track->Eta(),fCent,fZ},1.0);
					}
					else{
						fHistos->FillTH2("hPhiEtaCor",track->Phi(),track->Eta(),1.0/
						Eff[ binPt.FindBin(track->Pt())-1 ][ (int)((track->Eta()-fEff_eta_min)/fEff_eta_l) ] );
						FillTHnSparse("hTrackDataCor",{track->Pt(),track->Phi(),track->Eta(),fCent,fZ},1.0/
						Eff[ binPt.FindBin(track->Pt())-1 ][ (int)((track->Eta()-fEff_eta_min)/fEff_eta_l) ] );
					}
				}
				else if( fOption.Contains("Add3DEff") ){
					if( track->Pt() > fEff_pT_max ){
						fHistos->FillTH2("hPhiEtaCor",track->Phi(),track->Eta(),1.0);
						FillTHnSparse("hTrackDataCor",{track->Pt(),track->Phi(),track->Eta(),fCent,fZ},1.0);
					}
					else{
						fHistos->FillTH2("hPhiEtaCor",track->Phi(),track->Eta(),1.0/
						Eff3D[ binPt.FindBin(track->Pt())-1 ][ (int)((track->Eta()-fEff_eta_min)/fEff_eta_l) ][ (int)(track->Phi()/(2.0*pi/fEff_nphi_step)) ] );
						FillTHnSparse("hTrackDataCor",{track->Pt(),track->Phi(),track->Eta(),fCent,fZ},1.0/
						Eff3D[ binPt.FindBin(track->Pt())-1 ][ (int)((track->Eta()-fEff_eta_min)/fEff_eta_l) ][ (int)(track->Phi()/(2.0*pi/fEff_nphi_step)) ] );
					}
				}
				else if( fOption.Contains("AddpTEff") ){
					if( track->Pt() > fEff_pT_max ){
						fHistos->FillTH2("hPhiEtaCor",track->Phi(),track->Eta(),1.0);
	                                        FillTHnSparse("hTrackDataCor",{track->Pt(),track->Phi(),track->Eta(),fCent,fZ},1.0);
					}
					else{
						fHistos->FillTH2("hPhiEtaCor",track->Phi(),track->Eta(),1.0/
	                                        EffpT[ binPt.FindBin(track->Pt())-2 ] );
	                                        FillTHnSparse("hTrackDataCor",{track->Pt(),track->Phi(),track->Eta(),fCent,fZ},1.0/
	                                        EffpT[ binPt.FindBin(track->Pt())-2 ] );
					}
				}
			}
			else if( fOption.Contains("ITS") ){
				if( fOption.Contains("AddEff") ){
					if( track->Pt() > fEff_pT_max ){
						FillTHnSparse("hTrackDataCor",{track->Pt(),track->Phi(),track->Eta(),fCent,fZ},1.0);
					}
					else{
						FillTHnSparse("hTrackDataCor",{track->Pt(),track->Phi(),track->Eta(),fCent,fZ},1.0/
						Eff[ binPt.FindBin(track->Pt())-1 ][ (int)((track->Eta()-ITS_fEff_eta_min)/ITS_fEff_eta_l) ] );
					}
				}				
			}

		}
		fNTracks++;
		goodtrackindices.push_back(it);

		//Event mixing pool		
		if (fsetmixing){
			etl->push_back( (AliVTrack*) track -> Clone() );
		}
	}//track

	fHistos->FillTH1("hLHPt",LHPt,1.0);
	fHistos->FillTH2("hLHPt_JetpT",LHPt,fJetPt,1.0);
	if (fsetmixing){
		if (!goodtrackindices.size()) ep->pop_back();
		if ( ep->size() > bookingsize ){
			for (auto it: ep->front()) delete it;
			ep->pop_front();
		}
	}
	return goodtrackindices.size();
}

Bool_t AliAnalysisTaskRidge::GoodTracksSelectionMC(){

	goodtrackindicesMCALICE.clear();
	goodtrackindicesMCCMS.clear();

        trackMClist *etlMCALICE;
        eventpoolMC *epMCALICE;

        trackMClist *etlMCCMS;
        eventpoolMC *epMCCMS;

	if (fsetmixing && centbin>=0 && zbin>=0 && IsMC){
		epMCALICE = &fEMpoolMCALICE[centbin][zbin];
		epMCALICE->push_back( trackMClist() );
		etlMCALICE = &(epMCALICE->back());

		epMCCMS = &fEMpoolMCCMS[centbin][zbin];
		epMCCMS->push_back( trackMClist() );
		etlMCCMS = &(epMCCMS->back());
	}

        if( IsMC && fEvt->IsA()==AliAODEvent::Class() ){
		if( fabs(fZ_gen)<AbsZmax ){
                        const Int_t nTracksMC = fMCArray->GetEntriesFast();
                        for(Int_t iTracks = 0; iTracks < nTracksMC; iTracks++){
                                AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
                                if( !trackMC ) continue;
                                if( !(trackMC->IsPhysicalPrimary()) ) continue;
                                if( trackMC->Charge() == 0 ) continue;
				if( !fOption.Contains("ITS") ){
					if( ALICEAccp( trackMC->Pt(),trackMC->Eta() ) ){
						goodtrackindicesMCALICE.push_back(iTracks);
						if( fsetmixing ){
							Particlelet par;
							par.eta = trackMC->Eta();
							par.phi = trackMC->Phi();
							par.pt  = trackMC->Pt();
							par.IsTrackRecon = (bool)fEvt->GetTrack( trackMC->GetLabel() );
							etlMCALICE->push_back( par );
						}
						if( !fOption.Contains("CMSOnly") ) FillTHnSparse("hTrackMCALICE",{trackMC->Pt(),trackMC->Phi(),trackMC->Eta(),fCent,fZ_gen},1.0);
					}
				}
				if( fOption.Contains("ITS") ){
					if( trackMC->Pt() > fptcut && fabs( trackMC->Eta() ) < 1.3 ){
                                                goodtrackindicesMCALICE.push_back(iTracks);
                                                if( fsetmixing ){
                                                        Particlelet par;
                                                        par.eta = trackMC->Eta();
                                                        par.phi = trackMC->Phi();
                                                        par.pt  = trackMC->Pt();
                                                        par.IsTrackRecon = (bool)fEvt->GetTrack( trackMC->GetLabel() );
                                                        etlMCALICE->push_back( par );
                                                }
                                                if( !fOption.Contains("CMSOnly") ) FillTHnSparse("hTrackMCALICE",{trackMC->Pt(),trackMC->Phi(),trackMC->Eta(),fCent,fZ_gen},1.0);
                                        }
				}
				
				if( CMSAccp( trackMC->Pt(),trackMC->Eta() ) ){
					goodtrackindicesMCCMS.push_back(iTracks);
					if( fsetmixing ){
						Particlelet par;
						par.eta = trackMC->Eta();
						par.phi = trackMC->Phi();
						par.pt  = trackMC->Pt();
						par.IsTrackRecon = (bool)fEvt->GetTrack( trackMC->GetLabel() );
						etlMCCMS->push_back( par );
					}
					if( fOption.Contains("AddCMS") ) FillTHnSparse("hTrackMCCMS",{trackMC->Pt(),trackMC->Phi(),trackMC->Eta(),fCent,fZ_gen},1.0);
				}
                        }
			if(fsetmixing){
				if( !goodtrackindicesMCALICE.size() ) epMCALICE->pop_back();
				if( epMCALICE->size() > bookingsizeMC ){
					epMCALICE->pop_front();
				}

				if( !goodtrackindicesMCCMS.size() ) epMCCMS->pop_back();
				if( epMCCMS->size() > bookingsizeMC ){
					epMCCMS->pop_front();
				}

			}
                }
        }
	if( IsMC ) return goodtrackindicesMCCMS.size();
	else{ return 0; }
}
Bool_t AliAnalysisTaskRidge::GoodTrackletSelection(){

	trackletlist *etl;
	eventpooltracklet *ep;
	if (fsetmixing && centbin>=0 && zbin>=0){
		ep = &fEMpooltracklet[centbin][zbin];
		ep -> push_back( trackletlist() ); //
		etl = &(ep->back());
	}


	goodtrackindices.clear();
 
	fMultiplicity = fEvt -> GetMultiplicity();
	Int_t ntraklets = fMultiplicity->GetNumberOfTracklets();
	Double_t eta = 0, phi=0; 
	for (UInt_t it = 0; it<ntraklets; it++){
		eta = fMultiplicity->GetEta(it);
		phi = fMultiplicity->GetPhi(it);
		if( fabs(eta)>1.99 ) continue;
		if( TVector2::Phi_0_2pi(phi) >2*pi-0.01 ) continue;
		if( TVector2::Phi_0_2pi(phi) < 0.01 ) continue;
		goodtrackindices.push_back(it);
		fHistos->FillTH2("hPhiEta",phi,eta);
		if( fOption.Contains("SPDEff") )
                fHistos->FillTH2("hPhiEtaCor",phi,eta,
		1.0/Eff[ (int)( TVector2::Phi_0_2pi(phi)*180.0/(2*pi) ) ][ (int)( ( eta+2.0 )*40.0/4.0 ) ]
                );
		if (fsetmixing) {
			Tracklet tracklet;
			tracklet.eta = eta;
			tracklet.phi = phi;
			etl->push_back( tracklet );
		}
	}
	if (fsetmixing){
		if (!goodtrackindices.size()) ep->pop_back();
		if ( ep->size() > bookingsize ){
			ep->pop_front();
		}
	}
	return goodtrackindices.size();

}


void AliAnalysisTaskRidge::FillTracks(){


	NTracksPerPtBin.clear();
	NTracksPerPtBin.resize( binTPt.GetNbins() );

        for(int i=0;i<binTPt.GetNbins();i++){
                NTracksPerPtBin[i] =0 ;
        }

	AliVTrack *track1, *track2;
	TLorentzVector temp1,temp2;
	TLorentzVector vecsum;
	const UInt_t ntracks = goodtrackindices.size();

	tracklist trackpool;
/*
        int epsize=1;	
	if (fsetmixing && centbin>=0 && zbin>=0){
		eventpool &ep = fEMpool[centbin][zbin];
		epsize = ep.size();
		if (ep.size()<bookingsize  ) return;
		int n = 0;
		for (auto pool: ep){
			if (n == (ep.size() -1 )) continue;			
			for (auto track: pool) trackpool.push_back((AliVTrack*)track);
			n++;
		}
	}
*/
/*
	std::random_device rd;
	std::default_random_engine engine{rd()};
	std::shuffle ( goodtrackindices.begin(), goodtrackindices.end(), engine );
*/

	double MaxPhi=0;
	double MaxPt=0;
	double MaxEta=0;

	double effi;

	if( ntracks > 0 )
	for (UInt_t  it = 0; it < ntracks; it++) {
		track1 = (AliVTrack*) fEvt->GetTrack(goodtrackindices[it]);
		if( !track1 ) continue;

		effi = 1.0;
		if( !fOption.Contains("ITS") ){
			if( fOption.Contains("AddEff") ){
				if( binTPt.FindBin( track1->Pt() )-1 >= 0 ){
					if( track1->Pt() < fEff_pT_max ){
						effi = 1.0/Eff[ binPt.FindBin(track1->Pt())-1 ][ (int)((track1->Eta()-fEff_eta_min)/fEff_eta_l) ];
					}
					else{ effi = 1.0; }
				}
			}
			else if( fOption.Contains("Add3DEff") ){
				if( binTPt.FindBin( track1->Pt() )-1 >= 0 ){
					if( track1->Pt() < fEff_pT_max ){
						effi = 1.0/Eff3D[ binPt.FindBin(track1->Pt())-1 ][ (int)((track1->Eta()-fEff_eta_min)/fEff_eta_l) ][ (int)(track1->Phi()/(2.0*pi/fEff_nphi_step)) ];
					}
					else{ effi = 1.0; }
				}
			}
			else{ effi = 1.0; }
		}
		else if( fOption.Contains("ITS") ){
			if( fOption.Contains("AddEff") ){
				if( binTPt.FindBin( track1->Pt() )-1 >= 0 ){
					if( track1->Pt() < fEff_pT_max ){
						effi = 1.0/Eff[ binPt.FindBin(track1->Pt())-1 ][ (int)((track1->Eta()-ITS_fEff_eta_min)/ITS_fEff_eta_l) ];
					}
					else{ effi = 1.0; }
				}
			}
		}

		if( fOption.Contains("Asso") ){ effi = 1.0; }
//			effi = 1.0;
		if( binTPt.FindBin( track1->Pt() )-1 >= 0 ){
			NTracksPerPtBin[ binTPt.FindBin( track1->Pt() )-1 ] += effi;
		}
//		if( binTPt.FindBin( track1->Pt() )-1 >= 0 ){ NTracksPerPtBin[ binTPt.FindBin( track1->Pt() )-1 ]++; }
		if( MaxPt < track1->Pt() ){
			MaxPt = track1->Pt();
			MaxPhi = track1->Phi();
			MaxEta = track1->Eta();
		}
	}

	fLT_pT = MaxPt;
	if( ntracks > 1 ){
		FillTHnSparse("hTrackDataLTRaw",{fCent,MaxPt,MaxPhi,MaxEta,fZ},1.0);
	}
        for(int i=0;i<binTPt.GetNbins();i++){
//		if( NTracksPerPtBin[i] > 0.5 )
		FillTHnSparse("hNtrig",{fCent,binTPt.GetBinCenter(i+1),1.0,MaxPt},NTracksPerPtBin[i]);
		FillTHnSparse("hNtrigJet",{fCent,binTPt.GetBinCenter(i+1),1.0,fJetPt},NTracksPerPtBin[i]);
        }

	double PhiThres = MaxPhi;

	PhiThres = TVector2::Phi_0_2pi(PhiThres);
	// count
	if( PhiThres > pi*0.5 && PhiThres < pi*1.0 ){ PhiThres = pi - PhiThres; }
	else if( PhiThres > pi*1.0 && PhiThres < pi*1.5 ){ PhiThres = PhiThres - pi; }
	else if( PhiThres > 1.5*pi ){ PhiThres = 2*pi - PhiThres; }

        double eff1;
        double eff2;

//	double addPhi = 15.0 * pi / 180.0;	
	double addPhi = 0;
	if( ntracks > 0 )
	for (UInt_t  it = 0; it < ntracks-1; it++) {
		track1 = (AliVTrack*) fEvt->GetTrack(goodtrackindices[it]);
		if (!track1) continue;
		if( fOption.Contains("MyTrack") ){
			if( ( track1->Phi() > PhiThres && track1->Phi() < pi - PhiThres ) ||
			( track1->Phi() > pi + PhiThres && track1->Phi() < 2.0*pi - PhiThres ) ){
				continue;
			}
		}
		if( fOption.Contains("UETrack") ){
			if( PhiThres > pi/4.0-addPhi ){
				if( ( track1->Phi() > PhiThres-pi/4.0+addPhi && track1->Phi() < PhiThres+pi/4.0-addPhi ) ||
				( track1->Phi() > pi+PhiThres-pi/4.0+addPhi && track1->Phi() < pi+PhiThres+pi/4.0-addPhi ) ){
					continue;
				}
			}
			else if( PhiThres < pi/4.0-addPhi ){
				if( ( track1->Phi() > 0 && track1->Phi() < PhiThres+pi/4.0-addPhi ) ||
				( track1->Phi() > pi+PhiThres-pi/4.0+addPhi && track1->Phi() < pi+PhiThres+pi/4.0-addPhi ) || 
				( track1->Phi() > 2.0*pi+PhiThres-pi/4.0+addPhi ) ){
					continue;
				}
			}
		}
		for (UInt_t jt = it+1; jt < ntracks; jt++) {
			track2 = (AliVTrack*) fEvt->GetTrack(goodtrackindices[jt]);
			if (!track2) continue;
			if( fOption.Contains("MyTrack") ){
				if( ( track2->Phi() > PhiThres && track2->Phi() < pi - PhiThres ) ||
				( track2->Phi() > pi + PhiThres && track2->Phi() < 2.0*pi - PhiThres ) ){
					continue;
				}
			}
                	if( fOption.Contains("UETrack") ){
                	        if( PhiThres > pi/4.0-addPhi ){
                	                if( ( track2->Phi() > PhiThres-pi/4.0+addPhi && track2->Phi() < PhiThres+pi/4.0-addPhi ) ||
                	                ( track2->Phi() > pi+PhiThres-pi/4.0+addPhi && track2->Phi() < pi+PhiThres+pi/4.0-addPhi ) ){
                	                        continue;
                	                }
                	        }
                	        else if( PhiThres < pi/4.0-addPhi ){
                	                if( ( track2->Phi() > 0 && track2->Phi() < PhiThres+pi/4.0-addPhi ) ||
                	                ( track2->Phi() > pi+PhiThres-pi/4.0+addPhi && track2->Phi() < pi+PhiThres+pi/4.0-addPhi ) ||
                	                ( track2->Phi() > 2.0*pi+PhiThres-pi/4.0+addPhi ) ){
                	                        continue;
                	                }
                	        }
                	}
			double deltaeta = track1->Eta() - track2->Eta();
			double deltaphi = track1->Phi() - track2->Phi();
			if( track1-> Pt() < track2-> Pt()){
				deltaeta *= -1.0;
				deltaphi *= -1.0;
			}
			deltaphi = TVector2::Phi_0_2pi(deltaphi);

			if( deltaphi > 1.5*pi ) deltaphi -= 2.0*pi; 

			if( !fOption.Contains("ITS") ){
                        	if( fOption.Contains("AddEff") ){
                        	        if( track1->Pt() < fEff_pT_max ){
                        	                eff1 = Eff[ binPt.FindBin(track1->Pt())-1 ][ (int)((track1->Eta()-fEff_eta_min)/fEff_eta_l) ];
                        	        }
                        	        else{ eff1 = 1.0; }
                        	        if( track2->Pt() < fEff_pT_max ){
						eff2 = Eff[ binPt.FindBin(track2->Pt())-1 ][ (int)((track2->Eta()-fEff_eta_min)/fEff_eta_l) ];
                        	        }
                        	        else{ eff2 = 1.0; }
					if( fOption.Contains("Asso") ){
						if( track1-> Pt() > track2-> Pt() ){
							eff1 = 1.0;
						}
						else if( track1-> Pt() < track2-> Pt() ){
							eff2 = 1.0;
						}
					}
				}
				else if( fOption.Contains("Add3DEff") ){
					if( track1->Pt() < fEff_pT_max ){
						eff1 = Eff3D[ binPt.FindBin(track1->Pt())-1 ][ (int)((track1->Eta()-fEff_eta_min)/fEff_eta_l) ][ (int)(track1->Phi()/(2.0*pi/fEff_nphi_step)) ];
					}
					else{ eff1 = 1.0; }
					if( track2->Pt() < fEff_pT_max ){
						eff2 = Eff3D[ binPt.FindBin(track2->Pt())-1 ][ (int)((track2->Eta()-fEff_eta_min)/fEff_eta_l) ][ (int)(track2->Phi()/(2.0*pi/fEff_nphi_step)) ];
					}
				}
				else if( fOption.Contains("AddpTEff") ){
					if( track1->Pt() < fEff_pT_max ){
						eff1 = EffpT[ binPt.FindBin(track1->Pt())-2 ];
        	                        }
        	                        else{ eff1 = 1.0; }
        	                        if( track2->Pt() < fEff_pT_max ){
						eff2 = EffpT[ binPt.FindBin(track2->Pt())-2 ];
        	                        }
        	                        else{ eff2 = 1.0; }
				}
				else{ eff1 = 1.0; eff2 = 1.0; }
			}
                        else if( fOption.Contains("ITS") ){
                                if( fOption.Contains("AddEff") ){
                                        if( track1->Pt() < fEff_pT_max ){
                                                eff1 = Eff[ binPt.FindBin(track1->Pt())-1 ][ (int)((track1->Eta()-ITS_fEff_eta_min)/ITS_fEff_eta_l) ];
                                        }
                                        else{ eff1 = 1.0; }
                                        if( track2->Pt() < fEff_pT_max ){
                                                eff2 = Eff[ binPt.FindBin(track2->Pt())-1 ][ (int)((track2->Eta()-ITS_fEff_eta_min)/ITS_fEff_eta_l) ];
                                        }
                                        else{ eff2 = 1.0; }
                                        if( fOption.Contains("Asso") ){
                                                if( track1-> Pt() > track2-> Pt() ){
                                                        eff1 = 1.0;
                                                }
                                                else if( track1-> Pt() < track2-> Pt() ){
                                                        eff2 = 1.0;
                                                }
                                        }
                                }
                        }


			if( NTracksPerPtBin[binTPt.FindBin( max(track1->Pt(),track2-> Pt()) )-1] > 0 ){
//				FillTHnSparse("hRidge",{fCent, deltaphi, deltaeta,
//					max(track1-> Pt(),track2-> Pt()),
//					min(track1-> Pt(),track2-> Pt()) },
//					1.0/ ( NTracksPerPtBin[binTPt.FindBin( max(track1->Pt(),track2-> Pt()) )-1]*eff1*eff2 ) );
				FillTHnSparse("hRidgeLT",{fCent, deltaphi, deltaeta,
					max(track1-> Pt(),track2-> Pt()),
					min(track1-> Pt(),track2-> Pt()),MaxPt},
					1.0/ ( eff1*eff2 ) );

                                FillTHnSparse("hRidgeJet",{fCent, deltaphi, deltaeta,
                                        max(track1-> Pt(),track2-> Pt()),
                                        min(track1-> Pt(),track2-> Pt()),fJetPt},
                                        1.0/ ( eff1*eff2 ) );
			}
//			FillTHnSparse("hRidgeNTrig",{fCent, deltaphi, deltaeta,
  //                              max(track1-> Pt(),track2-> Pt()),
    //                            min(track1-> Pt(),track2-> Pt()) },
      //                          1.0/ ( eff1*eff2 ) );

		}        
	}


        int epsize=1;
        if (fsetmixing && centbin>=0 && zbin>=0 &&
	centbin < binCent.GetNbins() && zbin < binZ.GetNbins() ){
                eventpool &ep = fEMpool[centbin][zbin];
                epsize = ep.size();
                if (ep.size()<bookingsize  ) return;
                int n = 0;
                for (auto pool: ep){
                        if (n == (ep.size() -1 )) continue;
                        for (auto track: pool) trackpool.push_back((AliVTrack*)track);
                        n++;
                }
        }

	if (fsetmixing){
		for (UInt_t  it = 0; it < ntracks; it++) {
			track1 =  (AliVTrack*) fEvt->GetTrack(goodtrackindices.at(it)) ;
			for (UInt_t jt = 0; jt < trackpool.size(); jt++) {
				track2 = trackpool.at(jt);
				if (track1-> Pt() < track2->Pt()) continue; //trigger pT > associated pT
                		if( fOption.Contains("MyTrack") ){
                        		if( ( track1->Phi() > PhiThres && track1->Phi() < pi - PhiThres ) ||
                        		( track1->Phi() > pi + PhiThres && track1->Phi() < 2.0*pi - PhiThres ) ){ continue; } }
				if( fOption.Contains("MyTrack") ){
                                        if( ( track2->Phi() > PhiThres && track2->Phi() < pi - PhiThres ) ||
                                        ( track2->Phi() > pi + PhiThres && track2->Phi() < 2.0*pi - PhiThres ) ){ continue; } }
		                if( fOption.Contains("UETrack") ){
		                        if( PhiThres > pi/4.0-addPhi ){
		                                if( ( track1->Phi() > PhiThres-pi/4.0+addPhi && track1->Phi() < PhiThres+pi/4.0-addPhi ) ||
		                                ( track1->Phi() > pi+PhiThres-pi/4.0+addPhi && track1->Phi() < pi+PhiThres+pi/4.0-addPhi ) ){
		                                        continue;
		                                }
		                        }
		                        else if( PhiThres < pi/4.0-addPhi ){
	 	       	                        if( ( track1->Phi() > 0 && track1->Phi() < PhiThres+pi/4.0-addPhi ) ||
	        	                        ( track1->Phi() > pi+PhiThres-pi/4.0+addPhi && track1->Phi() < pi+PhiThres+pi/4.0-addPhi ) ||
	        	                        ( track1->Phi() > 2.0*pi+PhiThres-pi/4.0+addPhi ) ){
	        	                                continue;
	        	                        }
	        	                }
				}
	                        if( fOption.Contains("UETrack") ){
	                                if( PhiThres > pi/4.0-addPhi ){
	                                        if( ( track2->Phi() > PhiThres-pi/4.0+addPhi && track2->Phi() < PhiThres+pi/4.0-addPhi ) ||
	                                        ( track2->Phi() > pi+PhiThres-pi/4.0+addPhi && track2->Phi() < pi+PhiThres+pi/4.0-addPhi ) ){
	                                                continue;
	                                        }
	                                }
	                                else if( PhiThres < pi/4.0-addPhi ){
	                                        if( ( track2->Phi() > 0 && track2->Phi() < PhiThres+pi/4.0-addPhi ) ||
	                                        ( track2->Phi() > pi+PhiThres-pi/4.0+addPhi && track2->Phi() < pi+PhiThres+pi/4.0-addPhi ) ||
	                                        ( track2->Phi() > 2.0*pi+PhiThres-pi/4.0+addPhi ) ){
	                                                continue;
	                                        }
	                                }
	                        }

				double deltaeta = track1->Eta() - track2->Eta();
				double deltaphi = TVector2::Phi_0_2pi(track1->Phi() - track2->Phi());
				if( deltaphi > 1.5*pi ) deltaphi -= 2.0*pi; 

				if( !fOption.Contains("ITS") ){
                                	if( fOption.Contains("AddEff") ){
                                	        if( track1->Pt() < fEff_pT_max ){
							eff1 = Eff[ binPt.FindBin(track1->Pt())-1 ][ (int)((track1->Eta()-fEff_eta_min)/fEff_eta_l) ];
                                	        }
                                	        else{ eff1 = 1.0; }
	
	                                        if( track2->Pt() < fEff_pT_max ){
							eff2 = Eff[ binPt.FindBin(track2->Pt())-1 ][ (int)((track2->Eta()-fEff_eta_min)/fEff_eta_l) ];
	                                        }
	                                        else{ eff2 = 1.0; }
	                                	if( fOption.Contains("Asso") ){
	                                	        if( track1-> Pt() > track2-> Pt() ){
	                                	                eff1 = 1.0;
	                                	                eff2 = Eff[ binPt.FindBin(track2->Pt())-1 ][ (int)((track2->Eta()-fEff_eta_min)/fEff_eta_l) ];
	                                	        }
	                                	        else if( track1-> Pt() < track2-> Pt() ){
	                                	                eff1 = Eff[ binPt.FindBin(track1->Pt())-1 ][ (int)((track1->Eta()-fEff_eta_min)/fEff_eta_l) ];
	                                	                eff2 = 1.0;
	                                	        }
	                                	}
	                                }
					else if( fOption.Contains("Add3DEff") ){
						if( track1->Pt() < fEff_pT_max ){
							eff1 = Eff3D[ binPt.FindBin(track1->Pt())-1 ][ (int)((track1->Eta()-fEff_eta_min)/fEff_eta_l) ][ (int)(track1->Phi()/(2.0*pi/fEff_nphi_step)) ];
						}
						else{ eff1 = 1.0; }
						if( track2->Pt() < fEff_pT_max ){
							eff2 = Eff3D[ binPt.FindBin(track2->Pt())-1 ][ (int)((track2->Eta()-fEff_eta_min)/fEff_eta_l) ][ (int)(track2->Phi()/(2.0*pi/fEff_nphi_step)) ];
						}
					}
					else if( fOption.Contains("AddpTEff") ){
	                                        if( track1->Pt() < fEff_pT_max ){
							eff1 = EffpT[ binPt.FindBin(track1->Pt())-2 ];
	                                        }
	                                        else{ eff1 = 1.0; }
	
	                                        if( track2->Pt() < fEff_pT_max ){
	   						eff2 = EffpT[ binPt.FindBin(track2->Pt())-2 ];
	                                        }
	                                        else{ eff2 = 1.0; }
					}
					else{ eff1 = 1.0; eff2 = 1.0; }
				}


                                else if( fOption.Contains("ITS") ){
	                                if( fOption.Contains("AddEff") ){
	                                        if( track1->Pt() < fEff_pT_max ){
	                                                eff1 = Eff[ binPt.FindBin(track1->Pt())-1 ][ (int)((track1->Eta()-ITS_fEff_eta_min)/ITS_fEff_eta_l) ];
	                                        }
	                                        else{ eff1 = 1.0; }
	                                        if( track2->Pt() < fEff_pT_max ){
	                                                eff2 = Eff[ binPt.FindBin(track2->Pt())-1 ][ (int)((track2->Eta()-ITS_fEff_eta_min)/ITS_fEff_eta_l) ];
	                                        }
	                                        else{ eff2 = 1.0; }
	                                        if( fOption.Contains("Asso") ){
	                                                if( track1-> Pt() > track2-> Pt() ){
	                                                        eff1 = 1.0;
	                                                }
	                                                else if( track1-> Pt() < track2-> Pt() ){
	                                                        eff2 = 1.0;
	                                                }
	                                        }
	                                }
                                        else{ eff1 = 1.0; eff2 = 1.0; }
                                }
//				FillTHnSparse("hRidgeMixing",{fCent, deltaphi, deltaeta,
//					max(track1-> Pt(),track2-> Pt()),
  //                              	min(track1-> Pt(),track2-> Pt()) }, 1.0/(epsize-1)/ntracks/eff1/eff2 );
				if( NTracksPerPtBin[binTPt.FindBin( max(track1->Pt(),track2-> Pt()) )-1] > 0 ){
//					FillTHnSparse("hRidgeMixingS",{fCent, deltaphi, deltaeta,
//						max(track1-> Pt(),track2-> Pt()),
//						min(track1-> Pt(),track2-> Pt()) },
//						1.0/(NTracksPerPtBin[binTPt.FindBin( max(track1->Pt(),track2-> Pt()) )-1]*eff1*eff2) );
					FillTHnSparse("hRidgeMixingSLT",{fCent, deltaphi, deltaeta,
						max(track1-> Pt(),track2-> Pt()),
						min(track1-> Pt(),track2-> Pt()),MaxPt},
						1.0/(eff1*eff2) );


                                        FillTHnSparse("hRidgeMixingSJet",{fCent, deltaphi, deltaeta,
                                                max(track1-> Pt(),track2-> Pt()),
                                                min(track1-> Pt(),track2-> Pt()),fJetPt},
                                                1.0/(eff1*eff2) );
				}
//				FillTHnSparse("hRidgeMixingSNTrig",{fCent, deltaphi, deltaeta,
//					max(track1-> Pt(),track2-> Pt()),
//					min(track1-> Pt(),track2-> Pt()) }, 1.0/(eff1*eff2) );
			}
		}
	}
}

void AliAnalysisTaskRidge::FillTracklets(){


	trackletlist trackpool;
	int epsize=1;	
	if (fsetmixing && centbin>=0 && zbin>=0){
		eventpooltracklet &ep = fEMpooltracklet[centbin][zbin];
		epsize = ep.size();
		if (ep.size()< bookingsize/20 ) return;
		int n = 0;
		for (auto pool: ep){
			if (n == (ep.size() -1 )) continue;			
			for (auto track: pool) trackpool.push_back(track);
			n++;
		}
	}

	//eta is in ordering, so detaleta is always negative. solution -> Shuffle
/*
	std::random_device rd;
	std::default_random_engine engine{rd()};
	std::shuffle ( goodtrackindices.begin(), goodtrackindices.end(), engine );
*/

	fMultiplicity = fEvt -> GetMultiplicity();
	const UInt_t ntracks = goodtrackindices.size();
	Double_t eta1 = 0, phi1=0; 
	Double_t eta2 = 0, phi2=0; 
	for (UInt_t it = 0; it<ntracks-1; it++){
		eta1 = fMultiplicity -> GetEta(goodtrackindices[it]);
		phi1 = fMultiplicity -> GetPhi(goodtrackindices[it]);
		for (UInt_t jt = it+1; jt<ntracks; jt++){
			eta2 = fMultiplicity -> GetEta(goodtrackindices[jt]);
			phi2 = fMultiplicity -> GetPhi(goodtrackindices[jt]);
			double deltaeta = eta1 - eta2;
			double deltaphi = TVector2::Phi_0_2pi(phi1 - phi2);
			if( deltaphi > 1.5*pi ) deltaphi -= 2.0*pi; 
			if( !fOption.Contains("SPDEff") )
			FillTHnSparse("hRidge",{fCent, deltaphi, deltaeta,2,2 },1.0/ntracks );
			else if( fOption.Contains("SPDEff") )
			FillTHnSparse("hRidge",{fCent, deltaphi, deltaeta,2,2 },
			1.0/( ntracks* 
			Eff[ (int)( TVector2::Phi_0_2pi(phi1)*180.0/(2*pi) ) ][ (int)( ( eta1+2.0 )*40.0/4.0 ) ]*
			Eff[ (int)( TVector2::Phi_0_2pi(phi2)*180.0/(2*pi) ) ][ (int)( ( eta2+2.0 )*40.0/4.0 ) ] ) );
		}
	}

	if (fsetmixing){
		for (UInt_t it = 0; it<ntracks; it++){
			eta1 = fMultiplicity -> GetEta(goodtrackindices[it]);
			phi1 = fMultiplicity -> GetPhi(goodtrackindices[it]);
			for (UInt_t jt = 0; jt<trackpool.size(); jt++){
				eta2 = trackpool.at(jt).eta;
				phi2 = trackpool.at(jt).phi;
			
				double deltaeta = eta1 - eta2;
				double deltaphi = TVector2::Phi_0_2pi(phi1 - phi2);
				if( deltaphi > 1.5*pi ) deltaphi -= 2.0*pi; 
				if( !fOption.Contains("SPDEff") )
				FillTHnSparse("hRidgeMixing",{fCent, deltaphi, deltaeta,1,1 },1.0/(ntracks*(epsize-1)) );
				else if( fOption.Contains("SPDEff") )
				FillTHnSparse("hRidgeMixing",{fCent, deltaphi, deltaeta,1,1 },
                        	1.0/( ntracks*(epsize-1)*
                        	Eff[ (int)( TVector2::Phi_0_2pi(phi1)*180.0/(2*pi) ) ][ (int)( ( eta1+2.0 )*40.0/4.0 ) ]*
                        	Eff[ (int)( TVector2::Phi_0_2pi(phi2)*180.0/(2*pi) ) ][ (int)( ( eta2+2.0 )*40.0/4.0 ) ] ) );

			}

		}
	}
}


void AliAnalysisTaskRidge::FillTracksMC(){


        NTracksPerPtBinMCALICE.clear();
        NTracksPerPtBinMCALICE.resize( binTPt.GetNbins() );
        NTracksPerPtBinMCCMS.clear();
        NTracksPerPtBinMCCMS.resize( binTPt.GetNbins() );

	for(int i=0;i<binTPt.GetNbins();i++){
		NTracksPerPtBinMCALICE[i] =0 ;
		NTracksPerPtBinMCCMS[i] = 0;
	}

        AliAODMCParticle *par1;
	AliAODMCParticle *par2;

	Particlelet parMixing;

        const UInt_t nTracksMCALICE = goodtrackindicesMCALICE.size();
        const UInt_t nTracksMCCMS = goodtrackindicesMCCMS.size();

	trackMClist trkpoolMCALICE;
	trackMClist trkpoolMCCMS;

        int epsizeALICE=0, epsizeCMS=0;
	int n=0;
        if( fsetmixing && centbin>=0 && zbin>=0 && IsMC ){

		eventpoolMC &epMCALICE = fEMpoolMCALICE[centbin][zbin];
		eventpoolMC &epMCCMS = fEMpoolMCCMS[centbin][zbin];

                epsizeALICE = epMCALICE.size();
		if( epMCALICE.size() < bookingsizeMC ) return;
		n=0;
                for( auto pool: epMCALICE ){
                        if( n == (epMCALICE.size() -1 ) ) continue;
                        for( auto track: pool ) trkpoolMCALICE.push_back(track);
                        n++;
                }

		epsizeCMS = epMCCMS.size();
		if( epMCCMS.size() < bookingsizeMC ) return;
		n=0;
		for( auto pool: epMCCMS ){
			if( n == (epMCCMS.size() -1 ) ) continue;
			for( auto track: pool ) trkpoolMCCMS.push_back(track);
			n++;
		}
        }

	double deltaeta;
	double deltaphi;


	double MaxPtALICE=0;
	double MaxPhiALICE;
	double MaxEtaALICE;
	double MaxPtCMS=0;
	double MaxPhiCMS;
	double MaxEtaCMS;

        if( IsMC && fEvt->IsA()==AliAODEvent::Class() ){
                if( fabs(fZ_gen)<AbsZmax ){
			for(Int_t iTracks = 0; iTracks < nTracksMCALICE; iTracks++){
				par1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(goodtrackindicesMCALICE[iTracks]));
				if( !par1 ) continue;
				NTracksPerPtBinMCALICE[ binTPt.FindBin( par1->Pt() )-1 ]++;
				if( MaxPtALICE < par1->Pt() ){
					MaxPtALICE = par1->Pt();
					MaxPhiALICE = par1->Phi();
					MaxEtaALICE = par1->Eta();
				}
			}

			for(Int_t iTracks = 0; iTracks < nTracksMCCMS; iTracks++){
				par1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(goodtrackindicesMCCMS[iTracks]));
				if( !par1 ) continue;
				NTracksPerPtBinMCCMS[ binTPt.FindBin( par1->Pt() )-1 ]++;
				if( MaxPtCMS < par1->Pt() ){
					MaxPtCMS = par1->Pt();
					MaxPhiCMS = par1->Phi();
					MaxEtaCMS = par1->Eta();
				}
			}
		}
	}


	if( nTracksMCALICE > 1 && !fOption.Contains("CMSOnly") ) FillTHnSparse("hTrackMCALICELT",{fCent,MaxPtALICE,MaxPhiALICE,MaxEtaALICE,fZ_gen},1.0);
	if( nTracksMCCMS > 1 && fOption.Contains("AddCMS") ) FillTHnSparse("hTrackMCCMSLT",{fCent,MaxPtCMS,MaxPhiCMS,MaxEtaCMS,fZ_gen},1.0);

        for(int i=0;i<binTPt.GetNbins();i++){
//                if( NTracksPerPtBinMCALICE[i] > 0.5 )
		if( !fOption.Contains("CMSOnly") ) FillTHnSparse("hNtrigMCALICE",{fCent,binTPt.GetBinCenter(i+1),1.0,MaxPtALICE},NTracksPerPtBinMCALICE[i]);
//                if( NTracksPerPtBinMCCMS[i] > 0.5 )
		if( fOption.Contains("AddCMS") ) FillTHnSparse("hNtrigMCCMS",{fCent,binTPt.GetBinCenter(i+1),1.0,MaxPtCMS},NTracksPerPtBinMCCMS[i]);
        }


	if( IsMC && fEvt->IsA()==AliAODEvent::Class() ){
		if( fabs(fZ_gen)<AbsZmax ){
			if( nTracksMCALICE > 0 )
                        for(Int_t iTracks = 0; iTracks < nTracksMCALICE-1; iTracks++){
                                par1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(goodtrackindicesMCALICE[iTracks]));
                                if( !par1 ) continue;
				for(Int_t jTracks = iTracks+1; jTracks < nTracksMCALICE; jTracks++){
					par2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(goodtrackindicesMCALICE[jTracks]));
					if( !par2 ) continue;

                                	deltaeta = par1->Eta() - par2->Eta();
                                	deltaphi = par1->Phi() - par2->Phi();
                                	if( par1->Pt() < par2->Pt() ){
                                	        deltaeta *= -1.0;
                                	        deltaphi *= -1.0;
                                	}
                                	deltaphi = TVector2::Phi_0_2pi(deltaphi);
                                	if( deltaphi > 1.5*pi ) deltaphi -= 2.0*pi;

					if( NTracksPerPtBinMCALICE[binTPt.FindBin( max(par1->Pt(),par2->Pt()) )-1] > 0 ){
//						FillTHnSparse("hRidgeMCALICE",{fCent,deltaphi,deltaeta,
  //                                      	        max(par1-> Pt(),par2-> Pt()), min(par1-> Pt(),par2-> Pt()) },
    //                                    	        1.0 / ( NTracksPerPtBinMCALICE[binTPt.FindBin( max(par1->Pt(),par2-> Pt()) )-1] ) );
						if( !fOption.Contains("CMSOnly") )
						FillTHnSparse("hRidgeMCALICELT",{fCent,deltaphi,deltaeta,
							max(par1-> Pt(),par2-> Pt()), min(par1-> Pt(),par2-> Pt()),MaxPtALICE},
							1.0 );
					}
//                                        FillTHnSparse("hRidgeMCALICENTrig",{fCent,deltaphi,deltaeta,
  //                                              max(par1-> Pt(),par2-> Pt()), min(par1-> Pt(),par2-> Pt()) },
    //                                            1.0 / ( 1.0 ) );

					if( (bool)(fEvt->GetTrack( par1->GetLabel() ) ) && (bool)(fEvt->GetTrack( par2->GetLabel() ) ) ){
//						FillTHnSparse("hRidgeMCALICET",{fCent,deltaphi,deltaeta,
//							max(par1-> Pt(),par2-> Pt()), min(par1-> Pt(),par2-> Pt()) },
//							1.0 / ( NTracksPerPtBinMCALICE[binTPt.FindBin( max(par1->Pt(),par2-> Pt()) )-1] ) );
					}
				}
			}



			for(Int_t iTracks = 0; iTracks < nTracksMCALICE; iTracks++){
				par1 = (AliAODMCParticle*)(fMCArray->At(goodtrackindicesMCALICE[iTracks]));
				if( !par1 ) continue;
				for(Int_t jTracks = 0; jTracks < trkpoolMCALICE.size(); jTracks++){
					parMixing = trkpoolMCALICE.at(jTracks);
					if( par1->Pt() < parMixing.pt ) continue;

					deltaeta = par1->Eta() - parMixing.eta;
					deltaphi = par1->Phi() - parMixing.phi;
					deltaphi = TVector2::Phi_0_2pi(deltaphi);
					if( deltaphi > 1.5*pi ) deltaphi -= 2.0*pi;

					if( NTracksPerPtBinMCALICE[binTPt.FindBin( max(par1->Pt(),parMixing.pt) )-1] > 0 ){
//						FillTHnSparse("hRidgeMixingSMCALICE",{fCent,deltaphi,deltaeta,
//							max(par1-> Pt(),parMixing.pt), min(par1-> Pt(),parMixing.pt) },
//							1.0 / ( NTracksPerPtBinMCALICE[binTPt.FindBin( max(par1->Pt(),parMixing.pt) )-1] ) );
						if( !fOption.Contains("CMSOnly") )
						FillTHnSparse("hRidgeMixingSMCALICELT",{fCent,deltaphi,deltaeta,
							max(par1-> Pt(),parMixing.pt), min(par1-> Pt(),parMixing.pt),MaxPtALICE},
							1.0 );
					}
  //                                      FillTHnSparse("hRidgeMixingSMCALICENTrig",{fCent,deltaphi,deltaeta,
    //                                            max(par1-> Pt(),parMixing.pt), min(par1-> Pt(),parMixing.pt) },
      //                                          1.0 / (1.0) );
					if( (bool)(fEvt->GetTrack( par1->GetLabel() ) ) && parMixing.IsTrackRecon ){
	//					FillTHnSparse("hRidgeMixingSMCALICET",{fCent,deltaphi,deltaeta,
	//						max(par1-> Pt(),parMixing.pt), min(par1-> Pt(),parMixing.pt) },
	//						1.0 / nTracksMCALICE );
					}
				}
                        }

			if( nTracksMCCMS > 0 && fOption.Contains("AddCMS") )
			for(Int_t iTracks = 0; iTracks < nTracksMCCMS-1; iTracks++){
				par1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(goodtrackindicesMCCMS[iTracks]));
				if( !par1 ) continue;

				for(Int_t jTracks = iTracks+1; jTracks < nTracksMCCMS; jTracks++){
					par2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(goodtrackindicesMCCMS[jTracks]));
					if( !par2 ) continue;

                                        deltaeta = par1->Eta() - par2->Eta();
                                        deltaphi = par1->Phi() - par2->Phi();
                                        if( par1->Pt() < par2->Pt() ){
                                                deltaeta *= -1.0;
                                                deltaphi *= -1.0;
                                        }
                                        deltaphi = TVector2::Phi_0_2pi(deltaphi);
                                        if( deltaphi > 1.5*pi ) deltaphi -= 2.0*pi;
					if( NTracksPerPtBinMCCMS[binTPt.FindBin( max(par1->Pt(),par2-> Pt()) )-1] > 0 ){
	//					FillTHnSparse("hRidgeMCCMS",{fCent,deltaphi,deltaeta,
	  //                                              max(par1-> Pt(),par2-> Pt()), min(par1-> Pt(),par2-> Pt()) },
	    //                                            1.0 / ( NTracksPerPtBinMCCMS[binTPt.FindBin( max(par1->Pt(),par2-> Pt()) )-1] ) );

						FillTHnSparse("hRidgeMCCMSLT",{fCent,deltaphi,deltaeta,
							max(par1-> Pt(),par2-> Pt()), min(par1-> Pt(),par2-> Pt()),MaxPtCMS},
							1.0  );
					}
//                                        FillTHnSparse("hRidgeMCCMSNTrig",{fCent,deltaphi,deltaeta,
  //                                              max(par1-> Pt(),par2-> Pt()), min(par1-> Pt(),par2-> Pt()) },
    //                                            1.0 / (1.0) );
					if( (bool)(fEvt->GetTrack( par1->GetLabel() ) ) && (bool)(fEvt->GetTrack( par2->GetLabel() ) ) ){
//                                                FillTHnSparse("hRidgeMCCMST",{fCent,deltaphi,deltaeta,
  //                                                      max(par1-> Pt(),par2-> Pt()), min(par1-> Pt(),par2-> Pt()) },
    //                                                    1.0 / ( NTracksPerPtBinMCCMS[binTPt.FindBin( max(par1->Pt(),par2-> Pt()) )-1] ) );
                                        }
				}
			}


			if( fOption.Contains("AddCMS") )
			for(Int_t iTracks = 0; iTracks < nTracksMCCMS; iTracks++){
				par1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(goodtrackindicesMCCMS[iTracks]));	
				if( !par1 ) continue;
                                for(Int_t jTracks = 0; jTracks < trkpoolMCCMS.size(); jTracks++){
					parMixing = trkpoolMCCMS.at(jTracks);
                                        if( par1->Pt() < parMixing.pt ) continue;

                                        deltaeta = par1->Eta() - parMixing.eta;
                                        deltaphi = par1->Phi() - parMixing.phi;
                                        deltaphi = TVector2::Phi_0_2pi(deltaphi);
                                        if( deltaphi > 1.5*pi ) deltaphi -= 2.0*pi;
					if( NTracksPerPtBinMCCMS[binTPt.FindBin( max(par1->Pt(),parMixing.pt) )-1] > 0 ){
//	                                        FillTHnSparse("hRidgeMixingSMCCMS",{fCent,deltaphi,deltaeta,
//	                                                max(par1-> Pt(),parMixing.pt), min(par1-> Pt(),parMixing.pt) },
//	                                                1.0 / ( NTracksPerPtBinMCCMS[binTPt.FindBin( max(par1->Pt(),parMixing.pt) )-1] ) );
						FillTHnSparse("hRidgeMixingSMCCMSLT",{fCent,deltaphi,deltaeta,
							max(par1-> Pt(),parMixing.pt), min(par1-> Pt(),parMixing.pt),MaxPtCMS},
							1.0  );
					}
//					FillTHnSparse("hRidgeMixingSMCCMSNTrig",{fCent,deltaphi,deltaeta,
  //                                              max(par1-> Pt(),parMixing.pt), min(par1-> Pt(),parMixing.pt) },
    //                                            1.0 / (1.0) );
                                        if( (bool)(fEvt->GetTrack( par1->GetLabel() ) ) && parMixing.IsTrackRecon ){
//                                                FillTHnSparse("hRidgeMixingSMCCMST",{fCent,deltaphi,deltaeta,
  //                                                      max(par1-> Pt(),parMixing.pt), min(par1-> Pt(),parMixing.pt) },
    //                                                    1.0 / nTracksMCCMS );
                                        } 
                               }
			}
		}
	}
}

void AliAnalysisTaskRidge::RhoSparse(AliJetContainer *ktContainer, AliJetContainer * aktContainer , Int_t numberofexcludingjets) {
	// Lets exclude a dijet
	AliEmcalJet *leading = nullptr;
	AliEmcalJet *subleading = nullptr;
	Int_t n = 0;
	//for (auto ij : aktContainer->accepted_momentum())
	//	{
	//	auto j = ij.second;
	//	cout << "sg jet pt : " << j->Pt() << endl;
	//	n++;
	//}

	n = 0;
	Int_t njetacc = 0;
	static Double_t rhovec[999];
 	Double_t TotaljetAreaPhys=0;
  	Double_t TotalTPCArea=2*TMath::Pi()*0.9;
//	Double_t TotalTPCArea = 0.0;
	for (auto iBg : ktContainer->accepted_momentum())
	{
		if (n < numberofexcludingjets) {
			n++;
			continue;
		}
		auto bgjet = iBg.second;

		Bool_t matched = false;
		for (auto iSg : aktContainer->accepted_momentum())
		{
			auto sgjet = iSg.second;
			matched = (isOverlapping(bgjet, sgjet)) ? true : false;
		}
		//cout << "n = "<<n<< " kt jet pt : " << bgjet->Pt() << " matched : "<< matched << " jet eta = "<<bgjet->Eta()<< endl;

//		TotalTPCArea += bgjet->Area();
		if (bgjet -> GetNumberOfTracks()>0 && bgjet->Pt()>0.1 ){
			rhovec[njetacc] = bgjet->Pt() / bgjet->Area();
			TotaljetAreaPhys += bgjet->Area();
			njetacc++;
		}
			
		n++;
	}

  	Double_t OccCorr=1;
    OccCorr = TotaljetAreaPhys/TotalTPCArea;
	if (njetacc > 0 ) RHO = TMath::Median(njetacc, rhovec);
	else
		RHO = 0;

	RHO *= OccCorr;
	//cout << "jet pt end\n\n"
	//	 << endl;
	cout << 2*TMath::Pi()*0.9 << "< " << TotalTPCArea << ", " << TotaljetAreaPhys << endl;
}

Bool_t AliAnalysisTaskRidge::isOverlapping(AliEmcalJet *jet1, AliEmcalJet *jet2)
{
	for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i)
	{
		Int_t jet1Track = jet1->TrackAt(i);
		for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j)
		{
			Int_t jet2Track = jet2->TrackAt(j);
			if (jet1Track == jet2Track)
				return kTRUE;
		}
	}
	return kFALSE;
}

void AliAnalysisTaskRidge::MeasureBgDensity(AliJetContainer* ktContainer){
	using TMath::Abs;
        RHO=0;
        RHOM=0;
        int n = 0;
        TLorentzVector1D rhoarray;
        Double1D Sumpt;
        Double1D Summ;
        TLorentzVector leadingkt;
        Bool_t isfirstdijet = true;

//        for( auto j : ktContainer->accepted() ){
	for(int i=0;i<ktContainer->GetNJets();i++){
		AliEmcalJet* j = (AliEmcalJet*)ktContainer->GetJet(i);
//		cout << "j->Pt() : " << j->Pt() << endl;
                if (fabs(j->Eta())>0.7)  continue;
		if( j->Pt() <1e-10 ) continue;
//		if( !fJetTask->AcceptJet(j) ) continue;
                double lpt = 0;
                double sumpt = 0;
                double summ = 0;
                TLorentzVector sumkt (0,0,0,0);
                for( int it=0; it<j->GetNumberOfTracks(); it++ ) {
                        auto trk =  j->Track(it);
                        if( ! ((AliAODTrack*) trk)->TestFilterBit(768)) continue;
                        TLorentzVector temp;
                        temp.SetXYZM(trk->Px(),trk->Py(),trk->Pz(),AliPID::ParticleMass(AliPID::kPion));
                        if( lpt < temp.Pt() )  lpt = temp.Pt();
			fHistos->FillTH1("hPtCons",temp.Pt(),1.0);
                        sumkt += temp;
//                        sumpt += temp.Pt();
			sumpt += trk->Pt();
                        summ += (sqrt(AliPID::ParticleMass(AliPID::kPion)*AliPID::ParticleMass(AliPID::kPion)+temp.Pt()*temp.Pt())-temp.Pt());
                }
                sumpt /= j->Area();
                summ /= j->Area();

                if (n==0) { //remove leading kt jet
                        leadingkt = sumkt;
                        n++;
                        continue;
                }
                //remove back-subleading kt jet
                else if (fabs(sumkt.DeltaPhi(leadingkt))>pi/2. && isfirstdijet){
                        n++;
                        isfirstdijet = false;
                        continue;
                }
                Sumpt.push_back(sumpt);
                Summ.push_back(summ);
                n++;
        }

        Double_t rhopt[Sumpt.size()];
        Double_t rhom[Summ.size()];
        int i=0;
        for (auto k : Sumpt){
                rhopt[i] = k;
                i++;
        }
        i=0;
        for (auto k : Summ){
                rhom[i] = k;
                i++;
        }
        RHO = TMath::Median(Sumpt.size(),rhopt);
        RHOM = TMath::Median(Summ.size(),rhom);
//        FillTHnSparse("hRho",{fCent,RHO, pthardbin},sf);
}

//___________________________________________________________________
void AliAnalysisTaskRidge::FinishTaskOutput()
{
    //fOutput->Write();
}
//___________________________________________________________________
void AliAnalysisTaskRidge::Terminate(Option_t*)
{

}
//___________________________________________________________________

bool AliAnalysisTaskRidge::ALICEAccp(double pt, double eta){
	if( pt > fptcut && fabs( eta ) < fetacut && !TMath::IsNaN(pt) ){ return true; }
	return 0; 
}

bool AliAnalysisTaskRidge::CMSAccp(double pt, double eta){
	if( pt > 0.1 && fabs( eta ) < 2.4 && !TMath::IsNaN(pt)){ return true; }
        return 0;
}



//Int_t AliAnalysisTaskRidge::GetPID(AliPIDResponse *pid, const AliVTrack *trk){
//    if (!pid) return -1; // no pid available

    /*Double_t sigmas[] ={-999,-999,-999,-999};

    Int_t ipid = kUnknown;
    Double_t lsigma = 3;
    sigmas[kPion] = pid -> NumberOfSigmasTPC(trk,AliPID::kPion);
    sigmas[kKaon] = pid -> NumberOfSigmasTPC(trk,AliPID::kKaon);
    sigmas[kProton] = pid -> NumberOfSigmasTPC(trk,AliPID::kProton);
    sigmas[kElectron] = pid -> NumberOfSigmasTPC(trk,AliPID::kElectron);
    for (int i=0; i<kUnknown; i++){
        if (fabs(sigmas[i]) < lsigma) {
            lsigma = fabs(sigmas[i]);
            ipid = i;
        }
    }

    // derive information, whether tof pid is available
    if (0){
        const Bool_t ka = !(trk->GetStatus() & AliESDtrack::kTOFmismatch);
        const Bool_t kb =  (trk->GetStatus() & AliESDtrack::kTOFpid);
        const Bool_t ktof = ka && kb;
    }

   if (lsigma>3 ) return kUnknown;
   else  return ipid; 

	*/
/*
	Double_t prob[AliPID::kSPECIES];
	fPIDCombined->ComputeProbabilities(trk,pid,prob);
	Int_t ipid = AliPID::kUnknown;
	Double_t iprob = 0;
	for (int i=0; i<AliPID::kSPECIES; i++){
		if (prob[i]>0.6 && prob[i]>iprob) {
			iprob = prob[i];
			ipid = i;
		}
	}
  if (ipid == AliPID::kUnknown) ipid = AliPID::kPion;
	return ipid;

}
*/

THnSparse * AliAnalysisTaskRidge::CreateTHnSparse(TString name
	, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt){
  const TAxis * axises[bins.size()];
  for( UInt_t i=0;i<bins.size();i++ ) axises[i]= &bins[i];
  THnSparse * h= fHistos->CreateTHnSparse(name, title, ndim, axises,opt );
  return h;
}

THnSparse * AliAnalysisTaskRidge::CreateTHnSparse(TString name
	, TString title, TString templ, Option_t * opt){
  auto o = fHistos->FindObject(templ);
  if( !o ) {
    cout<<"ERROR: no "<<templ<<endl;
    gSystem->Exit(1);
  }
  auto ht = dynamic_cast<THnSparse*>( o );
  const TAxis * axises[ht->GetNdimensions()];
  for( int i=0;i<ht->GetNdimensions();i++ ) axises[i]= ht->GetAxis(i);
  auto h= fHistos->CreateTHnSparse(name, title, ht->GetNdimensions(), axises,opt );
  return h;
}

Long64_t AliAnalysisTaskRidge::FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w ){
  auto hsparse = dynamic_cast<THnSparse*>( fHistos->FindObject(name) );
  if(! hsparse ){
    cout<<"ERROR : no "<<name<<endl;
    exit(1);
  }
  return FillTHnSparse( hsparse, x, w );
}

Long64_t AliAnalysisTaskRidge::FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w ){
  if( int(x.size()) != h->GetNdimensions() ){
    cout<<"ERROR : wrong sized of array while Fill "<<h->GetName()<<endl;
    exit(1);
  }
  return h->Fill( &x.front(), w );
}

TAxis AliAnalysisTaskRidge::AxisFix
	( TString name, int nbin, Double_t xmin, Double_t xmax ){
		TAxis axis(nbin, xmin, xmax);axis.SetName(name);
		return axis;
}

TAxis AliAnalysisTaskRidge::AxisStr( TString name, std::vector<TString> bin ){
  TAxis ax = AxisFix( name, bin.size(), 0.5, bin.size()+0.5);
  UInt_t i=1;
  for( auto blabel : bin )
    ax.SetBinLabel( i++, blabel );
  return ax;
}

TAxis AliAnalysisTaskRidge::AxisVar( TString name, std::vector<Double_t> bin ){
  TAxis axis( bin.size()-1, &bin.front() ) ;axis.SetName(name);
  return axis;
}

TAxis AliAnalysisTaskRidge::AxisLog
	( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0){
  int binoffset = ( xmin0<0 || (xmin-xmin0)<1e-9) ? 0 : 1;
  std::vector<Double_t> bin(nbin+1+binoffset,0);
  double logBW3 = (log(xmax)-log(xmin))/nbin;
  for(int ij=0;ij<=nbin;ij++) bin[ij+binoffset]=xmin*exp(ij*logBW3);
  TAxis axis( nbin, &bin.front() ) ;
  axis.SetName(name);
  return axis;
}


