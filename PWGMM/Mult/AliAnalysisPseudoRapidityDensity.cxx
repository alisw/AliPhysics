
 /*************************************************************************

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
// Simple class for dn/deta analyses.
// by Beomkyu KIM
//==================================================================

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TGeoGlobalMagField.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliCentrality.h"
#include "AliVMultiplicity.h"
#include "AliPWG0Helper.h"
#include "AliMultSelection.h"
#include "AliESDtrack.h"
///#include "AliGenPythiaPlus.h"
#include "AliGenHijingEventHeader.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliITSMultReconstructor.h"
#include "AliITSsegmentationSPD.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGeomManager.h"
#include "AliESDInputHandlerRP.h"
#include "AliMagF.h"
#include "AliGRPObject.h"
#include "AliMultiplicity.h"
#include "AliVVZERO.h"
#include "AliAnalysisPseudoRapidityDensity.h"
const Double_t        pi = TMath::Pi();
//_________________________________________________________________
AliAnalysisPseudoRapidityDensity::AliAnalysisPseudoRapidityDensity()
	:AliAnalysisTaskSE("AliAnalysisPseudoRapidityDensity")
	 , fOption()
{
	DefineInput(0, TChain::Class());
	//DefineOutput (0, TTree::Class());
	DefineOutput (1, TList::Class());
}
//___________________________________________________________________
AliAnalysisPseudoRapidityDensity::AliAnalysisPseudoRapidityDensity
(
 const char *name
 , const char *option
 )
	:AliAnalysisTaskSE(name)
	 , fOption(option)
{
	DefineInput(0, TChain::Class());
	//DefineOutput(1, TTree::Class());
  DefineOutput (1, TList::Class());
}
//___________________________________________________________________
AliAnalysisPseudoRapidityDensity::AliAnalysisPseudoRapidityDensity
(
 const AliAnalysisPseudoRapidityDensity& ap
)
: fOption(ap.fOption)
{
}
//___________________________________________________________________
AliAnalysisPseudoRapidityDensity& AliAnalysisPseudoRapidityDensity::operator =
(
 const AliAnalysisPseudoRapidityDensity& ap
 )
{
	// assignment operator

	this->~AliAnalysisPseudoRapidityDensity();
	new(this) AliAnalysisPseudoRapidityDensity(ap);
	return *this;
}
//___________________________________________________________________
AliAnalysisPseudoRapidityDensity::~AliAnalysisPseudoRapidityDensity()
{
	//delete fHistos;
	delete fTrigger;
	delete fPIDResponse;
	delete fRandom;
}

//___________________________________________________________________
void AliAnalysisPseudoRapidityDensity::UserCreateOutputObjects()

{
	fRandom = new TRandom3;
	fRandom->SetSeed();
	// Histograms container

	// Offline triggers -----------------------------------------------------
	fTrigger = new AliTriggerAnalysis; // offline trigger
	//-----------------------------------------------------------------------


	// TrackCuts for strangeness measure------------------------------------i

	AliESDtrackCuts* tempcutset = new AliESDtrackCuts();
	fTrackCuts.resize(kTrackCutend);
	fTrackCuts.at(0) = *tempcutset;
	tempcutset ->GetStandardITSTPCTrackCuts2010(1,0);

	{
		TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
		tempcutset->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
		tempcutset->SetMinNClustersTPC(70);
		tempcutset->SetMaxChi2PerClusterTPC(4);
		tempcutset->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
		tempcutset->SetAcceptKinkDaughters(kFALSE);
		tempcutset->SetRequireTPCRefit(kTRUE);
		tempcutset->SetMaxFractionSharedTPCClusters(0.4);
		// ITS
		tempcutset->SetRequireITSRefit(kTRUE);
		//accept secondaries
		tempcutset->SetMaxDCAToVertexXY(2.4);
		tempcutset->SetMaxDCAToVertexZ(3.2);
		tempcutset->SetDCAToVertex2D(kTRUE);
		//reject fakes
		tempcutset->SetMaxChi2PerClusterITS(36);
		tempcutset->SetMaxChi2TPCConstrainedGlobal(36);

		tempcutset->SetRequireSigmaToVertex(kFALSE);

		tempcutset->SetEtaRange(-0.9,0.9);
		tempcutset->SetPtRange(0.15, 1E+15);
		tempcutset->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	}
	fTrackCuts.at(kHybrid) = *tempcutset;

	tempcutset ->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
	tempcutset ->SetRequireITSRefit(kTRUE);

	fTrackCutGC = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	fTrackCuts.at(kTPConly) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSPureSATrackCuts2010();
	fTrackCuts.at(kITSSA) = *tempcutset;

	tempcutset = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
	tempcutset -> SetMinNCrossedRowsTPC(120);
	tempcutset ->	SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	fTrackCuts.at(kITSTPC2010) = *tempcutset;

	/*
	AliESDtrackCuts dcaz(fTrackCuts.at(kITSTPC2010));
	dcaz.SetMaxDCAToVertexZ(1);
	fTrackCuts.at(kITSTPC2010dcazdw) = dcaz;
	dcaz.SetMaxDCAToVertexZ(5);
	fTrackCuts.at(kITSTPC2010dcazup) = dcaz;

	AliESDtrackCuts dcar(fTrackCuts.at(kITSTPC2010));
	dcar.SetMaxDCAToVertexXYPtDep("10*(0.0026+0.0050/pt^1.01)");
	fTrackCuts.at(kITSTPC2010dcarup) = dcar;
	dcar.SetMaxDCAToVertexXYPtDep("3*(0.0026+0.0050/pt^1.01)");
	fTrackCuts.at(kITSTPC2010dcardw) = dcar;

	AliESDtrackCuts nclutpc(fTrackCuts.at(kITSTPC2010));
	nclutpc.SetMinNCrossedRowsTPC(100);
	fTrackCuts.at(kITSTPC2010nclutpcdw) = nclutpc;
	nclutpc.SetMinNCrossedRowsTPC(130);
	fTrackCuts.at(kITSTPC2010nclutpcup) = nclutpc;

	AliESDtrackCuts chitpc(fTrackCuts.at(kITSTPC2010));
	chitpc.SetMaxChi2PerClusterTPC(3);
	fTrackCuts.at(kITSTPC2010chitpcdw) = chitpc;
	chitpc.SetMaxChi2PerClusterTPC(5);
	fTrackCuts.at(kITSTPC2010chitpcup) = chitpc;

	AliESDtrackCuts globalcons(fTrackCuts.at(kITSTPC2010));
	globalcons.SetMaxChi2TPCConstrainedGlobal(25);
	fTrackCuts.at(kITSTPC2010globalconsdw) = globalcons;
	globalcons.SetMaxChi2TPCConstrainedGlobal(49);
	fTrackCuts.at(kITSTPC2010globalconsup) = globalcons;
	*/


	fHistos = new THistManager("dndeta");

	//auto binType = AxisStr("Type",{"PN","PP","NN","Mixing"});
	Double1D varcentbin = {0,0.001,0.01,0.05,0.1,0.5,1,5,10,15,20,30,40,50,70,100};
	Double1D varcentbinHeavy = {0,2.5,5,7.5,10,20,30,40,50,60,70,80,90,100};
	//for (auto i=1; i<=100; i++) varcentbin.push_back(i);
	auto binCent = AxisVar("Cent", IsAA ? varcentbinHeavy : varcentbin);

	auto binEta = AxisFix("Eta",40,-2,2);
	auto binPhi = AxisVar("Phi",{0,pi*2/3,pi*4/3,2*pi});
	auto binCentClass = AxisStr("CentClass"
		,{"V0M","V0A","V0C","SPDMult"});
	auto binEventClass   = AxisStr( "EventClass", { "DATA", "INEL","NSD","INELg0" } );
	auto binTriggClass  = AxisStr("TriggClass",{"MBOR","MBAND","MBORg0","MBANDg0","HighMult","HighMultg0"});
	auto binParType = AxisStr("ParticleType"
			,{"Tracks","MotherStrange","Bkg","Pion","Kaon","Proton","Opar"});
	auto binV0Type = AxisStr("V0Type",{"K0s","Lambda","AntiLambda"});
	auto binPtVar = AxisStr("PtVar",{"NoPtVar","PtUp","PtDw","Hybrid"
		,"ITSTPC2010","TPConly","ITSSA"
		//,"ITSTPC2010dcazdw", "ITSTPC2010dcazup"
		//,"ITSTPC2010dcardw", "ITSTPC2010dcarup"
		//,"ITSTPC2010nclutpcdw", "ITSTPC2010nclutpcup"
		//,"ITSTPC2010chitpcdw", "ITSTPC2010chitpcup"
		//,"ITSTPC2010globalconsdw", "ITSTPC2010globalconsup"
		});

	Double1D binmul = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 25.5, 30.5, 40.5, 50.5, 100.5, 300.5};
	auto binMul = AxisVar("multiplicity",binmul);

	//auto binMass = AxisFix("Mass",200,0,2);

	std::vector<TString> ent ={"All","PS","PSpileup","Goodz","Goodzcut","INEL015","INEL050","INEL100", "INEL200", "Goodzcut7", "INELg010", "MBANDg010","INELg0","PSINELg0", "MBANDg0"};
	auto h = fHistos->CreateTH1("hEventNumbers","",ent.size(), 0, ent.size());
	fHistos->CreateTH1("quickcheck","",10, 0, 10);
	for(auto i=0u;i<ent.size();i++) h->GetXaxis()->SetBinLabel(i+1,ent.at(i).Data());
	fHistos -> CreateTH2("hPhiEta","",180,0,2*pi,40,-2,2);
	fHistos -> CreateTH2("hPhiEtaCut","",180,0,2*pi,40,-2,2);
	//Double1D zbins = {-20,-15,-10,-7,-3,0,3,7,10,15,20};
	//Double1D zbins = {-20,-15,-10,-7,-3,0,3,7,10,15,20};
	Double1D zbins = {-30,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,30};
	TAxis binZ = AxisVar("Z",zbins);

	fHistos -> CreateTH1("hdndeta","",60,-6,6);
	fHistos -> CreateTH1("hdndetacut","",60,-6,6);
	fHistos -> CreateTH1("hdndetamc","",60,-6,6);
	fHistos -> CreateTH1("hdndetaINEL015","",1200,-6,6);
	fHistos -> CreateTH1("hdndetaINEL050","",1200,-6,6);
	fHistos -> CreateTH1("hdndetaINEL100","",1200,-6,6);
	fHistos -> CreateTH1("hdndetaINEL200","",1200,-6,6);
	fHistos -> CreateTH1("hGenMultV0AGeo","",100000,0,1000);
	fHistos -> CreateTH1("hGenMultV0CGeo","",100000,0,1000);
	fHistos -> CreateTH1("hGenMultV0MGeo","",100000,0,1000);
	fHistos -> CreateTH1("hGenMultSPDGeo","",100000,0,1000);
	fForward2d = fHistos ->CreateTH2 ("hfmddetadphi","",200, -4, 6, 20, 0, 2*pi);
	fHistos -> CreateTH1("hcentTracklets","cent",100,0,100);
	fHistos -> CreateTH1("hcentselection","cent",100,0,100);
	fHistos -> CreateTH1("ztruth","ztuth",60,-30,30);
	fHistos -> CreateTH1("zdata","zdata",60,-30,30);
	fHistos -> CreateTH1("hchi2","SPD chi2",600,0,6);
	fHistos -> CreateTH1("hcentV0M","cent",100,0,100);
	fHistos -> CreateTH1("hcentV0MHighMult","cent",1000,0,1);
	fHistos -> CreateTH1("hcentSPD","cent",100,0,100);
	fHistos -> CreateTH1("hcentV0Mzcut7","cent",100,0,100);
	fHistos -> CreateTH1("hcentSPDzcut7","cent",100,0,100);
	fHistos -> CreateTH1("hcentV0Mzcut10","cent",100,0,100);
	fHistos -> CreateTH1("hcentSPDzcut10","cent",100,0,100);
	fHistos -> CreateTH1("hcentV0Mzcut20","cent",100,0,100);
	fHistos -> CreateTH1("hcentSPDzcut20","cent",100,0,100);
	fHistos -> CreateTH1("hcentV0Mzcut15","cent",100,0,100);
	fHistos -> CreateTH1("hcentSPDzcut15","cent",100,0,100);
	fHistos -> CreateTH1("hINELg0Nch","hINELg0Nch",300,0,300);
	fHistos -> CreateTH1("hINT7g0Nch","hINT7g0Nch",300,0,300);
	fHistos -> CreateTH1("hEtTruth","hEtTruth",100,0,100,"s");
	fHistos -> CreateTH1("hEtMCrec","hEtMCrec",100,0,100,"s");
	fHistos -> CreateTH2("hMultResponseV0M","hMultResponseV0M", 200,0,200,800,0,800,"s");
	fHistos -> CreateTH2("hMultResponseSPD","hMultResponseSPD", 200,0,200,200,0,200,"s");


	const int nbins=100;
	Double_t logbins[nbins+1];
	Double_t low= 0.01;
	Double_t high=500;
	Double_t logbw= (log(high)-log(low))/nbins;
	for(int ij=0;ij<=nbins;ij++) logbins[ij]=low*exp(ij*logbw);
	fHistos -> CreateTH1 ("hsdmass","SD mass",nbins, logbins);
	low= 0.001;
	high=50;
	logbw= (log(high)-log(low))/nbins;
	for(int ij=0;ij<=nbins;ij++) logbins[ij]=low*exp(ij*logbw);
	fHistos -> CreateTH1 ("hkinept","Kine only pt",nbins, logbins);
	fHistos -> CreateTH1 ("hmcrecpt","MC rec pt",nbins, logbins);

	CreateTHnSparse( "hrecdndeta", "rec dndeta", 9, {
			binEventClass, binTriggClass, binCent, binZ, binParType, binPtVar, binEta
			,binCentClass, binPhi },"s");
	CreateTHnSparse( "hreczvtx", "rec zvtx", 5, {
			binEventClass, binTriggClass, binCent,binZ, binCentClass},"s");
	CreateTHnSparse( "heventcount", "event count", 3, {
			binTriggClass, binCent, binCentClass},"s");

	CreateTHnSparse( "hkinedndeta", "kine only dndeta", 8, {
			binEventClass, binCent, binZ, binParType, binPtVar, binEta
		, binCentClass, binPhi },"s");
	CreateTHnSparse( "hkinezvtx", "kine zvtx", 4, { binEventClass, binCent,binZ, binCentClass},"s");

	CreateTHnSparse( "hv0mass", "V0 inv mass", 5,
			{ binTriggClass, binCent,binV0Type, AxisFix("v0mass",1200,0.3,1.5), binCentClass},"s");
	CreateTHnSparse( "hv0eta", "V0 daughter eta", 5,
			{ binTriggClass, binCent,binV0Type, binEta, binCentClass},"s");

	CreateTHnSparse( "hkineGeoMult","",2,{binCentClass,AxisFix("GeoMult",30000,0,300)},"s");
	CreateTHnSparse( "hkinedNdEtaGeoCent","",3,{binCent,binEta,binCentClass},"s");
	CreateTHnSparse( "hcentGeo","",2,{binCentClass,binCent},"s");

	CreateTHnSparse( "hMult","",2,{binCentClass,AxisFix("Mult",20000,0,2000)},"s");
	CreateTHnSparse( "hMultHigh","",2,{binCentClass,AxisFix("Mult",20000,0,2000)},"s");
	
	CreateTHnSparse( "hMultcent","",3,{binCentClass,binCent,AxisFix("Multcent",20000,0,2000)},"s");
	CreateTHnSparse( "hMultHighcent","",3,{binCentClass,binCent,AxisFix("highMultcent",20000,0,2000)},"s");

	CreateTHnSparse( "hdNdEtaCent","",3,{binCent,binEta,binCentClass},"s");
	CreateTHnSparse( "hcent","",2,{binCentClass,binCent},"s");

	CreateTHnSparse("hPhiEtaTracks","",3,{binPtVar
		,AxisFix("phi",180,0,2*pi),AxisFix("eta",40,-2,2)},"s");


	PostData(1, fHistos->GetListOfHistograms());
	//PID Combined
	fPIDCombined = new AliPIDCombined;
	fPIDCombined->SetDefaultTPCPriors();//Need more update..
	fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);
	fPIDCombined->SetDetectorMask(
			AliPIDResponse::kDetTPC |
			AliPIDResponse::kDetTOF |
			AliPIDResponse::kDetITS |
			AliPIDResponse::kDetTRD);//Do we need??

	Double1D fcent(kCentClassBinEnd,-1);
	fCent = fcent;


}




//___________________________________________________________________
void AliAnalysisPseudoRapidityDensity::UserExec(Option_t* )
{


	// Pointer to a event----------------------------------------------------
	AliVEvent *event = InputEvent();
	if (!event) {
		Printf("ERROR: Could not retrieve event");
		return;
	}

	// connect to ESD tree --------------------------------------------------
	Int_t runnumber;
	event->IsA()==AliESDEvent::Class()
		? fEvt = dynamic_cast<AliESDEvent*>(event)
		: fEvt = dynamic_cast<AliAODEvent*>(event);
	if (!fEvt) return;


	if( IsFirstEvent ) {
		runnumber = fEvt->GetRunNumber();
		IsFirstEvent = false;
	}



	// ----------------------------------------------------------------------
	// centrality
	Double1D fcent(kCentClassBinEnd,-1);
	fCent = fcent;
	sel = (AliMultSelection*) fEvt -> FindListObject("MultSelection");
	//if ( sel->GetEvSelCode() <= 0 ) return;
	//if (! sel ->	GetThisEventPassesTrackletVsCluster() && !fOption.Contains("MC")) return;
	//sel->SetThisEventIsNotAsymmetricInVZERO(true);

	if (sel) {
		fCent[kV0M] = sel->GetMultiplicityPercentile("V0M");
		fCent[kV0A] = sel->GetMultiplicityPercentile("V0A");
		fCent[kV0C] = sel->GetMultiplicityPercentile("V0C");
		fCent[kSPDMult] = sel->GetMultiplicityPercentile("SPDTracklets");

	}
	//if(! sel->IsEventSelected() && !fOption.Contains("MC") && !fOption.Contains("LHC10d")) {
	//	fCent[kV0M] = sel->GetEvSelCode();
	//	fCent[kV0A] = sel->GetEvSelCode();
	//	fCent[kV0C] = sel->GetEvSelCode();
	//	fCent[kSPDMult] = sel->GetEvSelCode();
	//} 

	if(!fOption.Contains("MC") ) {
		fHistos -> FillTH1("quickcheck",1);
		if ( sel->GetThisEventIsNotPileup()) fHistos -> FillTH1("quickcheck",2); 
		if ( sel->GetThisEventIsNotPileupInMultBins()) fHistos -> FillTH1("quickcheck",3); 
		if ( sel->GetThisEventHasNoInconsistentVertices()) fHistos -> FillTH1("quickcheck",4); 
		if ( sel->GetThisEventPassesTrackletVsCluster()) fHistos -> FillTH1("quickcheck",5); 
		if ( !sel->GetThisEventIsNotPileup() 
				|| !sel->GetThisEventIsNotPileupInMultBins()
				|| !sel->GetThisEventHasNoInconsistentVertices()
				|| !sel->GetThisEventPassesTrackletVsCluster() 

			 ){
			fCent[kV0M] = -1;
			fCent[kV0A] = -1;
			fCent[kV0C] = -1;
			fCent[kSPDMult] = -1;;
		}
	}

	// Pointer to a MC event-------------------------------------------------
	AliMCEvent *mcEvent = MCEvent();
	// ----------------------------------------------------------------------

	Bool_t IsMC = kFALSE;
	TArrayF vtxMC3d(3);

	Double_t eta = -10, phi = -10, pt = -1.;
	AliGenPythiaEventHeader* pythiaGenHeader = NULL;
	AliGenDPMjetEventHeader* dpmHeader = NULL;
	AliGenHijingEventHeader* hijing = NULL;

	Bool_1d bevtc(kECend,false);
	sdweightingfactor = 1.;
	Bool_t issd = 0;
	Bool_t isdd = 0;
	Int_t Nch_INELg0_eta1 = 0;
	Double_t nv0mgeo = 0;
	Double_t nspdgeo = 0;
	if (mcEvent) {
		ismc = true;
		pythiaGenHeader =
			dynamic_cast<AliGenPythiaEventHeader*>(mcEvent->GenEventHeader());
		dpmHeader =
			dynamic_cast<AliGenDPMjetEventHeader*>(mcEvent->GenEventHeader());
		hijing =
			dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());

 
		if(pythiaGenHeader) { // Pythia6
			//92 SD1, 93 SD2, 94 DD, -2 CD , -1 ND, 91 EL
			Int_t proc = pythiaGenHeader->ProcessType();
			if (proc == 92 || proc == 93) issd = true;
			else if (proc == 94) isdd  = true;
			//else if (proc == -1 || proc == -2) bevtc[kND] = true;

			if (fOption.Contains("PYTHIA8")){
				// 103 SD1, 104 SD2, 105 DD
				if (proc == 103 || proc == 104) issd = true;
				else if (proc == 105) isdd = true;
			}
		}
		else if (dpmHeader) {
			Int_t proc = dpmHeader->ProcessType();
			if (proc == 5 || proc == 6) issd = true;
			else if (proc == 7) isdd = true;
		}

		if (!issd )  bevtc[kNSD]  = kTRUE;
		bevtc[kINEL] = true;


		if (fOption.Contains("EPOS")) bevtc[kINEL] = true; //EPOS doesn't have diff info



		stack = mcEvent -> Stack();
		Int_t nPrim  = stack->GetNprimary();
		mcEvent->GenEventHeader()->PrimaryVertex(vtxMC3d);
		if (bevtc[kINEL]) fHistos -> FillTH1("ztruth",vtxMC3d[2]);
		for (Int_t i = 0; i < nPrim; i++){
			TParticle* part = stack->Particle(i);
			if (AliPWG0Helper::IsPrimaryCharged(part, nPrim) == kFALSE) continue;
			eta = part->Eta();
			pt = part->Pt();
		}

		Double1D  mulGeo (kCentClassBinEnd,0);
		for (Int_t i = 0; i < nPrim; i++){
			TParticle* part = stack->Particle(i);
			if (AliPWG0Helper::IsPrimaryCharged(part, nPrim) == kFALSE) continue;
			eta = part->Eta();
			phi = part ->Phi();
			pt  = part->Pt();
			if ( bevtc[kINEL] && fabs(eta) < 1.0 ) bevtc[kINELg0] = true;
			if (abs(eta)<1) Nch_INELg0_eta1++;
			fHistos->FillTH1("hdndetamc",eta);
			fHistos->FillTH1("hkinept",pt,1./pt);
			if (abs(eta)<2) nspdgeo++;
			if ( (eta >-3.7 && eta< -1.7) || ( eta >2.8 && eta<5.1 )) nv0mgeo++;
		}
		Double1D geocent(kCentClassBinEnd,-1);
		if ( abs(vtxMC3d[2])<10 && bevtc [kINELg0]) fHistos ->FillTH1("hEtTruth",fCent[kV0M]);
	} else bevtc[kDATA] = true;



	if (fOption.Contains("MC") && bevtc[kINELg0] && abs(vtxMC3d[2])<10){
			fHistos -> FillTH1("hEventNumbers","INELg010",1);
			fHistos -> FillTH1("hINELg0Nch", Nch_INELg0_eta1, 1 );
	}
	if (fOption.Contains("MC") && bevtc[kINELg0] )
		fHistos -> FillTH1("hEventNumbers","INELg0",1);
 
	AliInputEventHandler* inputHandler = (AliInputEventHandler*) 
		AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
	Bool_t IsMinimumBias = kFALSE;
	fHistos -> FillTH1("hEventNumbers","All",1);
	Bool_1d btrigc(kTrigend,false);

	IsMinimumBias = inputHandler -> IsEventSelected() & AliVEvent::kINT7;
	btrigc[kMBAND] = IsMinimumBias;
	btrigc[kMBOR] = false; //let's regard CINT7 as a minimum bias trigger.
	if (!fOption.Contains("MC"))
		btrigc[kHighMult] = inputHandler -> IsEventSelected() & AliVEvent::kHighMultV0;
	else btrigc[kHighMult] = IsMinimumBias; 


	if (fOption.Contains("LHC10")){ 
		IsMinimumBias = inputHandler -> IsEventSelected() & AliVEvent::kMB;
		btrigc[kMBOR] = IsMinimumBias;
		btrigc[kMBAND] = IsMinimumBias
			&& fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kV0A) 
			&& fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kV0C);
	}



	if (IsMinimumBias)  fHistos -> FillTH1("hEventNumbers","PS",1);
	// Reject pile-up events--------------------------------------------------
	if (!mcEvent && event->IsPileupFromSPD(3.,0.8,3.,2.,5.)) return;



	if (IsMinimumBias) fHistos -> FillTH1("hEventNumbers","PSpileup",1);
	// -----------------------------------------------------------------------

	const AliVVertex* trackVtx  = fEvt->GetPrimaryVertexTracks() ;
	const AliVVertex* spdVtx	  = fEvt->GetPrimaryVertexSPD() ;

	Bool_t IsGoodVertex = kFALSE;
	Bool_t IsGoodVertexCut = kFALSE;

	double trackz=-50;
	fZ = -50;

	float vtxf[3] =  {-50,-50,-50};
	if (spdVtx) {
		fZ = spdVtx->GetZ();
		vtxf[0] = spdVtx->GetX();
		vtxf[1] = spdVtx->GetY();
		vtxf[2] = spdVtx->GetZ();
		IsGoodVertex = kTRUE;
		if (spdVtx->GetNContributors()<1) IsGoodVertex = kFALSE;
		else IsGoodVertex = true;
	} else IsGoodVertex = kFALSE;


	if (sel && !sel -> GetThisEventHasNoInconsistentVertices()) IsGoodVertex = false;
	if (sel && !sel -> GetThisEventHasGoodVertex2016 () && !fOption.Contains("LHC10d")) IsGoodVertex = false;
	
	if (IsMinimumBias && IsGoodVertex) fHistos -> FillTH1("hEventNumbers","Goodz",1);
	if (IsGoodVertex) fHistos -> FillTH1("zdata",vtxf[2]);

	if ( IsGoodVertex && fabs(fZ)<20.) {
		IsGoodVertexCut = kTRUE;
		if (IsMinimumBias) {
			fHistos -> FillTH1("hEventNumbers","Goodzcut",1);
			if (fabs(fZ)< 7)fHistos -> FillTH1("hEventNumbers","Goodzcut7",1);
		}
	}



	if (IsGoodVertexCut){
		fMultiplicity = fEvt -> GetMultiplicity();
		Double_t eta,pt;
		Double1D  mul (kCentClassBinEnd,0);
		Int_t ntrks = fMultiplicity->GetNumberOfTracklets();
		Int_t ntrks_int7g0 = 0;
		for (auto it = 0; it<ntrks; it++) {
			Double_t eta = fMultiplicity->GetEta(it);
			Double_t phi = fMultiplicity->GetPhi(it);
			if ( abs(eta) < 1.0 && btrigc[ kMBOR ]) btrigc[ kMBORg0 ] = true;
			if ( abs(eta) < 1.0 && btrigc[ kMBAND ]) btrigc[ kMBANDg0 ] = true;
			if ( abs(eta) < 1.0 && btrigc[ kHighMult ]) btrigc[ kHighMultg0 ] = true;
			if ( abs(eta) < 1.0 ) ntrks_int7g0++;
			mul.at(kSPDMult)++;

		}
		//if ( abs(fZ)<10 && btrigc [kMBANDg0] &&  abs(vtxMC3d[2])<10 && bevtc [kINELg0]) fHistos ->FillTH1("hEtMCrec",fCent[kV0M]);
		
		//if (btrigc[kMBANDg0] && bevtc[kINELg0] && fOption.Contains("MC") )			fHistos -> FillTH1("hEventNumbers","MBANDg0",1);
		if (btrigc[kMBANDg0] && abs(fZ) < 10 && bevtc[kINELg0] && fOption.Contains("MC") && bevtc[kINELg0] && abs(vtxMC3d[2])<10) {
			fHistos -> FillTH1("hEventNumbers","MBANDg010",1);
			fHistos -> FillTH1("hINT7g0Nch", ntrks_int7g0, 1);
		}
		//mul.at(kSPDMult) = fTrackCuts.at(0).GetReferenceMultiplicity((AliESDEvent*)fEvt,  AliESDtrackCuts::kTrackletsITSTPC,2);

		Double_t intensity=0, intensityA=0, intensityC=0;
		AliVVZERO* lVV0 = fEvt->GetVZEROData();
		for (int i=0; i<64; i++) {
			intensity += lVV0->GetMultiplicity(i);
			if (i<32) intensityC += lVV0->GetMultiplicity(i); 
			else intensityA += lVV0->GetMultiplicity(i); 
		}
		mul.at(kV0M) = intensity;
		mul.at(kV0A) = intensityA;
		mul.at(kV0C)  = intensityC;
		Bool_t isfilling = false;
		for (auto i= 0u; i<kCentClassBinEnd; i++){
			Double_t rand = fRandom->Uniform(0,10);
			if (fabs(fZ)<10 && btrigc[kMBAND] && IsGoodVertex 
				&& sel->GetThisEventIsNotPileupInMultBins()
				&& sel->GetThisEventIsNotPileup() 
			  && sel->GetThisEventHasNoInconsistentVertices()
				&& sel->GetThisEventPassesTrackletVsCluster() ){
				if (i<=kV0C) FillTHnSparse("hMult",{double(i),mul.at(i)});
				else FillTHnSparse("hMult",{double(i),mul.at(i)+rand*0.1});
				if (i == kV0M) fHistos-> FillTH2( "hMultResponseV0M", nv0mgeo, mul.at(i)  );
				if (i == kSPDMult) fHistos-> FillTH2( "hMultResponseSPD", nspdgeo, mul.at(i)  );
				FillTHnSparse("hMultcent",{double(i),fCent.at(i),mul.at(i)});
				//else FillTHnSparse("hMult",{double(i),mul.at(i)+rand*0.01});
			}

			if (fabs(fZ)<10 && btrigc[kHighMultg0] && IsGoodVertex && sel->GetThisEventIsNotPileupInMultBins()){
				if (i<=kV0C) FillTHnSparse("hMultHigh",{double(i),mul.at(i)});
				FillTHnSparse("hMultHighcent",{double(i),fCent.at(i),mul.at(i)});
				//else FillTHnSparse("hMult",{double(i),mul.at(i)+rand*0.01});
			}
			

			if (fOption.Contains("LHC10d") || fOption.Contains("LHC12h")){
				if (i<=kV0C) fCent.at(i) = fMult[i].Integral(fMult[i].GetXaxis()->FindBin(mul.at(i)),fMult[i].GetNbinsX())*100.;
				else fCent.at(i) = fMult[i].Integral(fMult[i].GetXaxis()->FindBin(mul.at(i)+rand*0.1),fMult[i].GetNbinsX())*100.;
				     //fCent.at(i) = fMult[i].Integral(fMult[i].GetXaxis()->FindBin(mul.at(i)),fMult[i].GetNbinsX())*100.;
						//else fCent.at(i) = fMult[i].Integral(fMult[i].GetXaxis()->FindBin(mul.at(i)+rand*0.01),fMult[i].GetNbinsX())*100.;
				if (!sel->GetThisEventIsNotPileupInMultBins() && !fOption.Contains("MC")) fCent.at(i) = -1;
				if (!sel->GetThisEventIsNotPileup() && !fOption.Contains("MC")) fCent.at(i) = -1;
				if (!sel->GetThisEventHasNoInconsistentVertices() && !fOption.Contains("MC")) fCent.at(i) = -1;
				if (!sel->GetThisEventPassesTrackletVsCluster() && !fOption.Contains("MC")) fCent.at(i) = -1; 
				FillTHnSparse("hcent",{ double(i),fCent.at(i) });
			}
		}
		if (fOption.Contains("Flat") && fCent.at(kSPDMult)<0.002) fCent.at(kSPDMult) = -1; 

		if (fabs(fZ)< 7 && btrigc[ kMBANDg0 ]) {
			fHistos->FillTH1("hcentV0Mzcut7",fCent[kV0M]);
			fHistos->FillTH1("hcentSPDzcut7",fCent[kSPDMult]);
		}
		if (fabs(fZ)< 10 && btrigc[ kMBANDg0 ]) {
			fHistos->FillTH1("hcentV0Mzcut10",fCent[kV0M]);
			fHistos->FillTH1("hcentSPDzcut10",fCent[kSPDMult]);
		}
		if (fabs(fZ)< 20 && btrigc[ kMBANDg0 ]) {
			fHistos->FillTH1("hcentV0Mzcut20",fCent[kV0M]);
			fHistos->FillTH1("hcentSPDzcut20",fCent[kSPDMult]);
		}
		if (fabs(fZ)< 15 && btrigc[ kMBANDg0 ]) {
			fHistos->FillTH1("hcentV0Mzcut15",fCent[kV0M]);
			fHistos->FillTH1("hcentSPDzcut15",fCent[kSPDMult]);
		}

		const int nESDTracks = fEvt->GetNumberOfTracks();
		for(auto it = 0; it < fEvt->GetNumberOfTracks(); it++){
			AliESDtrack* track = (AliESDtrack*)fEvt->GetTrack(it);
			if (!track) continue;
			if (!fTrackCuts[kHybrid].AcceptTrack(track) && !fTrackCutGC.AcceptTrack(track)) continue;
			//if (!fTrackCuts[kITSTPC2010].AcceptTrack(track) ) continue;
			eta = track->Eta();
			pt  = track->Pt();
		}
		this -> FillTracklets(bevtc,btrigc);
		this -> StrangenessMeasure(btrigc);
	}

	//if (IsMinimumBias ){
		for (auto itrigc=1u ; itrigc<kTrigend; itrigc++){
			if (btrigc[itrigc]){
				for (auto icentc = UInt_t(kCentClassBinBegin); icentc <kCentClassBinEnd; icentc++){
					FillTHnSparse("heventcount",{Double_t(itrigc),fCent[icentc],Double_t(icentc)});
				}
			}
		}
	//}

	if (mcEvent) {
		for (auto ievtc=1u ; ievtc<kECend; ievtc++) {
			if (bevtc[ievtc])
				for (auto icentc = UInt_t(kCentClassBinBegin); icentc <kCentClassBinEnd; icentc++){
					FillTHnSparse("hkinezvtx"
							,{Double_t(ievtc),fCent[icentc],vtxMC3d[2],Double_t(icentc)},sdweightingfactor);
				}
		}
		Int_t nPrim  = stack->GetNprimary();

		Int_t pid ;
		for (Int_t i = 0; i < nPrim; i++){
			TParticle* part = stack->Particle(i);
			if (AliPWG0Helper::IsPrimaryCharged(part, nPrim) == kFALSE) continue;
			eta = part->Eta();
			phi = part ->Phi();
			pt  = part->Pt();
			switch (TMath::Abs(part->GetPdgCode())){
				case 211:   pid = kPion;   break;
				case 321:   pid = kKaon;   break;
				case 2212:  pid = kProton; break;
				default:    pid = kOPar; break;
			}
			for (auto ievtc=1u ; ievtc<=kECend; ievtc++) {
				if (bevtc[ievtc]){
					for (auto icentc = UInt_t(kCentClassBinBegin); icentc <kCentClassBinEnd; icentc++){
						FillTHnSparse( "hkinedndeta",
								{ Double_t(ievtc),fCent[icentc],vtxMC3d[2], Double_t(pid),double(kNoPtVar), eta,Double_t(icentc), phi},sdweightingfactor);
						if (pt<0.05 ){
							FillTHnSparse( "hkinedndeta",
									{ Double_t(ievtc),fCent[icentc],vtxMC3d[2], Double_t(pid),double(kPtUp),eta,Double_t(icentc),phi}, 2.*sdweightingfactor);
							FillTHnSparse( "hkinedndeta",
									{ Double_t(ievtc),fCent[icentc],vtxMC3d[2], Double_t(pid),double(kPtDw),eta,Double_t(icentc),phi}, 0.5*sdweightingfactor);
						} else {
							FillTHnSparse( "hkinedndeta",
									{ Double_t(ievtc),fCent[icentc],vtxMC3d[2], Double_t(pid),double(kPtUp),eta,Double_t(icentc),phi}, sdweightingfactor);
							FillTHnSparse( "hkinedndeta",
									{ Double_t(ievtc),fCent[icentc],vtxMC3d[2], Double_t(pid),double(kPtDw),eta,Double_t(icentc),phi}, sdweightingfactor);
						}
					}
				}
			}
		}




		
	}
	if (IsMinimumBias){
		fHistos->FillTH1("hcentV0M",fCent[kV0M]);
		fHistos->FillTH1("hcentV0MHighMult",fCent[kV0M]);
		fHistos->FillTH1("hcentSPD",fCent[kSPDMult]);
		if (fOption.Contains("MC") && bevtc[kINELg0] )
			fHistos -> FillTH1("hEventNumbers","PSINELg0",1);
	}


	PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisPseudoRapidityDensity::FillTracklets(Bool_1d bevtc, Bool_1d btrigc){
	UInt_t nmul = fMultiplicity->GetNumberOfTracklets();
	Double_t eta, phi, pt = -1., chi2 = 0;

	for (auto it=0u; it<nmul; it++) {
		eta = fMultiplicity->GetEta(it);
		phi = fMultiplicity->GetPhi(it);
		fHistos->FillTH2("hPhiEta",phi,eta);
		fHistos->FillTH1("hdndeta",eta);
		if (
			(phi>0 && phi<1.25) ||
			(phi>2.4 && phi<2.6) ||
			(phi>2.8 && phi<3.6) ||
			(phi>4.5 && phi<5.3) ||
			(phi>5.7 && phi<2*pi) 
			) 
		{
			fHistos->FillTH2("hPhiEtaCut",phi,eta);
			fHistos->FillTH1("hdndetacut",eta);
		}	

		chi2 = ((AliMultiplicity*) fMultiplicity ) -> CalcDist(it);
		fHistos -> FillTH1 ("hchi2",chi2);
		if (chi2>1.6) continue ;

		Int_t pid = 0;

		if (ismc && isbginjection) {
			pid = kBkg;
			pt = 1.;

		}
		else if (ismc){
			TParticle *particle = NULL;
			TParticle *mother = NULL;
			if (fMultiplicity->GetLabel(it,0)<0) continue;
			if (stack->GetNtrack() < fMultiplicity->GetLabel(it,0)
					|| stack->GetNtrack() <fMultiplicity->GetLabel(it,1))
				continue;

			particle = stack->Particle(fMultiplicity->GetLabel(it,0));
			if (!particle) continue;
			mother=AliPWG0Helper::FindPrimaryMother(stack,
					fMultiplicity->GetLabel(it,0));

			if (fMultiplicity->GetLabel(it,0) != fMultiplicity->GetLabel(it,1)
					&& fMultiplicity->GetLabel(it,1) >=0) {
				if (stack->IsPhysicalPrimary(fMultiplicity->GetLabel(it,0))==kFALSE
						&& stack->IsPhysicalPrimary(fMultiplicity->GetLabel(it,1))){
					particle = stack->Particle(fMultiplicity->GetLabel(it,1));
					if(!particle) continue;
					mother=AliPWG0Helper::FindPrimaryMother(stack,
							fMultiplicity->GetLabel(it,1));
				}
			}


			if (fMultiplicity->GetLabel(it,0) == fMultiplicity->GetLabel(it,1)){
				eta = particle -> Eta();
				phi = particle -> Phi();
			}
			else {
				eta = fMultiplicity->GetEta(it);
				phi  = fMultiplicity->GetPhi(it);
			}
			pt = particle -> Pt();
			fHistos -> FillTH1("hmcrecpt",pt,1./pt);

			switch (TMath::Abs(particle->GetPdgCode())){
				case 211:   pid = kPion;   break;
				case 321:   pid = kKaon;   break;
				case 2212:  pid = kProton; break;
				default:    pid = kOPar; break;
			}
			if (IsMotherStrangeParticle(mother)) pid = kMotherStrange;
			else if (!stack->IsPhysicalPrimary(fMultiplicity->GetLabel(it,0))
					&& !stack->IsPhysicalPrimary(fMultiplicity->GetLabel(it,1)))
				pid = kBkg;
		}

		if (!ismc) {
			pid = kParDATA;
		}

		for (auto ievtc=1u ; ievtc<kECend; ievtc++) {
			for (auto itrigc=1u ; itrigc<kTrigend; itrigc++) {
				if (bevtc[ievtc] && btrigc[itrigc] ){
					for (auto icentc = UInt_t(kCentClassBinBegin); icentc <kCentClassBinEnd; icentc++){
						FillTHnSparse( "hrecdndeta",
								{ Double_t(ievtc),Double_t(itrigc),fCent[icentc],fZ,Double_t(pid),double(kNoPtVar),eta,Double_t(icentc), phi}, sdweightingfactor);
						if (ismc){
							if (pt<0.05){
								FillTHnSparse( "hrecdndeta",
										{ Double_t(ievtc),Double_t(itrigc),fCent[icentc],fZ,Double_t(pid),double(kPtUp),eta,Double_t(icentc),phi }, 2.*sdweightingfactor);
								FillTHnSparse( "hrecdndeta",
										{ Double_t(ievtc),Double_t(itrigc),fCent[icentc],fZ,Double_t(pid),double(kPtDw),eta,Double_t(icentc),phi}, 0.5*sdweightingfactor);
							} else {
								FillTHnSparse( "hrecdndeta",
										{ Double_t(ievtc),Double_t(itrigc),fCent[icentc],fZ,Double_t(pid),double(kPtUp),eta,Double_t(icentc),phi}, sdweightingfactor);
								FillTHnSparse( "hrecdndeta",
										{ Double_t(ievtc),Double_t(itrigc),fCent[icentc],fZ,Double_t(pid),double(kPtDw),eta,Double_t(icentc),phi}, sdweightingfactor);
							}
						}
					}
				}
			}
		}
	}
	for (auto ievtc=1u ; ievtc<kECend; ievtc++) {
		for (auto itrigc=1u ; itrigc<kTrigend; itrigc++){
			if (bevtc[ievtc] && btrigc[itrigc] ){
				for (auto icentc = UInt_t(kCentClassBinBegin); icentc <kCentClassBinEnd; icentc++){
					FillTHnSparse("hreczvtx",{Double_t(ievtc),Double_t(itrigc)
							,fCent[icentc],fZ,Double_t(icentc)}, sdweightingfactor);
				}
			}
		}
	}
}

//refered from analysis2AliAODForwardMult.cxx
//refered from analysis2/AliForwardMultiplicityTask.cxx


//___________________________________________________________________
//void AliAnalysisPseudoRapidityDensity::FinishTaskOutput()
//{
	//OpenFile(1);
	//TH1D *fForward1d = (TH1D*)fForward2d->ProjectionX("_x",1,-1,"e");
	//TH1D* norm   = fForward2d->ProjectionX("norm", 0, 1, "");
	//fForward1d -> Divide(norm);
	//fForward1d -> Write();
//}
//___________________________________________________________________
void AliAnalysisPseudoRapidityDensity::Terminate(Option_t*)
{

}
//___________________________________________________________________

Int_t AliAnalysisPseudoRapidityDensity::GetPID(AliPIDResponse *pid, const AliVTrack *trk){
	if (!pid) return -1; // no pid available

	Double_t prob[AliPID::kSPECIES];
	fPIDCombined->ComputeProbabilities(trk,pid,prob);
	Int_t ipid = -1;
	Double_t iprob = 0;
	for (int i=0; i<AliPID::kSPECIES; i++){
		if (prob[i]>iprob) {
			iprob = prob[i];
			ipid = i;
		}
	}

	return ipid;

}


void AliAnalysisPseudoRapidityDensity::StrangenessMeasure( Bool_1d btrigc){
	Double_t tPrimaryVtxPosition[3];
	Int_t nv0s = 0;
	AliESDEvent *evt = dynamic_cast<AliESDEvent*> (fEvt);
	nv0s = evt->GetNumberOfV0s();
	Int_t    lOnFlyStatus = 0, nv0sOn = 0, nv0sOff = 0;
	Double_t lChi2V0 = 0;
	Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
	Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
	Double_t lV0CosineOfPointingAngle = 0;
	Double_t lV0Radius = 0;
	Double_t lV0DecayLength = 0;
	Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
	Double_t lPt       = 0, lRapK0s = 0, lRapLambda = 0;
	Double_t lAlphaV0  = 0, lPtArmV0 = 0;
	Double_t tV0Position[3];
	Double_t lMagneticField      = 999;
	const AliESDVertex *primaryVtx = evt->GetPrimaryVertex();
	tPrimaryVtxPosition[0] = primaryVtx->GetX();
	tPrimaryVtxPosition[1] = primaryVtx->GetY();
	tPrimaryVtxPosition[2] = primaryVtx->GetZ();
	if (TMath::Abs(tPrimaryVtxPosition[2])>10) return;

	for(Int_t iV0 = 0; iV0 < nv0s; iV0++){
		AliESDv0 *v0 = evt->GetV0(iV0);
		UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
		UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());
		lMagneticField = ((AliESDEvent*)evt)->GetMagneticField();


		AliESDtrack *pTrack = evt->GetTrack(lKeyPos);
		AliESDtrack *nTrack = evt->GetTrack(lKeyNeg);
		if (!pTrack || !nTrack) {
			Printf("ERROR: Could not retreive one of the daughter track");
			continue;
		}

		// Remove like-sign
		if ( pTrack->GetSign() == nTrack->GetSign()){
			//cout<< "like sign, continue"<< endl;
			continue;
		}
		// Tracks quality cuts
		if (((pTrack->GetTPCNcls())<80)||((nTrack->GetTPCNcls())<80))
			continue;

		// TPC refit condition (done during reconstruction for
		// Offline but not for On-the-fly)
		if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
		if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;

		if (pTrack) lDcaPosToPrimVertex =
			TMath::Abs(pTrack->GetD(
						tPrimaryVtxPosition[0]
						, tPrimaryVtxPosition[1]
						, lMagneticField) );

		if (nTrack) lDcaNegToPrimVertex =
			TMath::Abs(nTrack->GetD(
						tPrimaryVtxPosition[0]
						, tPrimaryVtxPosition[1]
						, lMagneticField) );

		lOnFlyStatus             = v0->GetOnFlyStatus();
		lChi2V0                  = v0->GetChi2V0();
		lDcaV0Daughters          = v0->GetDcaV0Daughters();
		lDcaV0ToPrimVertex       = v0->GetD(
				tPrimaryVtxPosition[0]
				, tPrimaryVtxPosition[1]
				, tPrimaryVtxPosition[2]);

		lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(
				tPrimaryVtxPosition[0]
				, tPrimaryVtxPosition[1]
				, tPrimaryVtxPosition[2]);
		v0->GetXYZ(tV0Position[0], tV0Position[1], tV0Position[2]);
		lV0Radius      = TMath::Sqrt(
				tV0Position[0]*tV0Position[0]
				+ tV0Position[1]*tV0Position[1]);
		lV0DecayLength = TMath::Sqrt(
				TMath::Power(tV0Position[0] - tPrimaryVtxPosition[0],2) +
				TMath::Power(tV0Position[1] - tPrimaryVtxPosition[1],2) +
				TMath::Power(tV0Position[2] - tPrimaryVtxPosition[2],2));

		// Pt:
		lPt = v0->Pt();

		// Armenteros variables: !!
		lAlphaV0      = v0->AlphaV0();
		lPtArmV0      = v0->PtArmV0();

		// Selections:
		if (1) {
			if ( (lDcaPosToPrimVertex      < 0.05 )||
					(lDcaNegToPrimVertex      < 0.05 )||
					(lDcaV0Daughters          > 0.5 )  ||
					(lV0CosineOfPointingAngle < 0.99)
				 ) continue;
		}



		for (auto icentc = UInt_t(kCentClassBinBegin); icentc <kCentClassBinEnd; icentc++){
			v0->ChangeMassHypothesis(310);
			for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
				if (btrigc[itrigc])
					FillTHnSparse("hv0mass",
							{Double_t(itrigc),fCent[icentc],Double_t(kK0s),v0->GetEffMass(),Double_t(icentc)});
			if (v0->GetEffMass() >0.482 && v0->GetEffMass()<0.509) {
				for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
					if (btrigc[itrigc]) {

						FillTHnSparse("hv0eta",
								{Double_t(itrigc),fCent[icentc],Double_t(kK0s),pTrack->Eta(),Double_t(icentc)});
						FillTHnSparse("hv0eta",
								{Double_t(itrigc),fCent[icentc],Double_t(kK0s),nTrack->Eta(),Double_t(icentc)});
					}
			}
			v0->ChangeMassHypothesis(3122);
			for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
				if (btrigc[itrigc])
					FillTHnSparse("hv0mass",
							{Double_t(itrigc),fCent[icentc],Double_t(kLambda),v0->GetEffMass(),Double_t(icentc)});
			if (v0->GetEffMass() >1.11 && v0->GetEffMass()<1.12) {
				for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
					if (btrigc[itrigc]) {
						FillTHnSparse("hv0eta",
								{Double_t(itrigc),fCent[icentc],Double_t(kLambda),pTrack->Eta(),Double_t(icentc)});
						FillTHnSparse("hv0eta",
								{Double_t(itrigc),fCent[icentc],Double_t(kLambda),nTrack->Eta(),Double_t(icentc)});
					}
			}

			v0->ChangeMassHypothesis(-3122);
			for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
				if (btrigc[itrigc])
					FillTHnSparse("hv0mass",
							{Double_t(itrigc),fCent[icentc],Double_t(kAntiLambda),v0->GetEffMass(),Double_t(icentc)});
			if (v0->GetEffMass() >1.11 && v0->GetEffMass()<1.12) {
				for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
					if (btrigc[itrigc]) {

						FillTHnSparse("hv0eta",
								{Double_t(itrigc),fCent[icentc],Double_t(kAntiLambda),pTrack->Eta(),Double_t(icentc)});
						FillTHnSparse("hv0eta",
								{Double_t(itrigc),fCent[icentc],Double_t(kAntiLambda),nTrack->Eta(),Double_t(icentc)});
					}
			}
		}
	}
}

THnSparse * AliAnalysisPseudoRapidityDensity::CreateTHnSparse(TString name
		, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt){
	const TAxis * axises[bins.size()];
	for( UInt_t i=0;i<bins.size();i++ ) axises[i]= &bins[i];
	THnSparse * h= fHistos->CreateTHnSparse(name, title, ndim, axises,opt );
	return h;
}

THnSparse * AliAnalysisPseudoRapidityDensity::CreateTHnSparse(TString name
		, TString title, TString templ, Option_t * opt){
	auto o = fHistos->FindObject(templ);
	if( !o ) {
		std::cout<<"ERROR: no "<<templ<<std::endl;
		gSystem->Exit(1);
	}
	auto ht = dynamic_cast<THnSparse*>( o );
	const TAxis * axises[ht->GetNdimensions()];
	for( int i=0;i<ht->GetNdimensions();i++ ) axises[i]= ht->GetAxis(i);
	auto h= fHistos->CreateTHnSparse(name, title, ht->GetNdimensions(), axises,opt );
	return h;
}

Long64_t AliAnalysisPseudoRapidityDensity::FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w ){
	auto hsparse = dynamic_cast<THnSparse*>( fHistos->FindObject(name) );
	if(! hsparse ){
		std::cout<<"ERROR : no "<<name<<std::endl;
		exit(1);
	}
	return FillTHnSparse( hsparse, x, w );
}

Long64_t AliAnalysisPseudoRapidityDensity::FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w ){
	if( int(x.size()) != h->GetNdimensions() ){
		std::cout<<"ERROR : wrong sized of array while Fill "<<h->GetName()<<std::endl;
		exit(1);
	}
	return h->Fill( &x.front(), w );
}

TAxis AliAnalysisPseudoRapidityDensity::AxisFix
( TString name, int nbin, Double_t xmin, Double_t xmax ){
	TAxis axis(nbin, xmin, xmax);axis.SetName(name);
	return axis;
}

TAxis AliAnalysisPseudoRapidityDensity::AxisStr( TString name, std::vector<TString> bin ){
	TAxis ax = AxisFix( name, bin.size(), 0.5, bin.size()+0.5);
	UInt_t i=1;
	for( auto blabel : bin )
		ax.SetBinLabel( i++, blabel );
	return ax;
}

TAxis AliAnalysisPseudoRapidityDensity::AxisVar( TString name, std::vector<Double_t> bin ){
	TAxis axis( bin.size()-1, &bin.front() ) ;axis.SetName(name);
	return axis;
}

TAxis AliAnalysisPseudoRapidityDensity::AxisLog
( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0){
	int binoffset = ( xmin0<0 || (xmin-xmin0)<1e-9) ? 0 : 1;
	std::vector<Double_t> bin(nbin+1+binoffset,0);
	double logBW3 = (log(xmax)-log(xmin))/nbin;
	for(int ij=0;ij<=nbin;ij++) bin[ij+binoffset]=xmin*exp(ij*logBW3);
	TAxis axis( nbin, &bin.front() ) ;
	axis.SetName(name);
	return axis;
}
