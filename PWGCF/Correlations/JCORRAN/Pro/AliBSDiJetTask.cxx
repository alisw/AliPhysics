/*************************************************************************
 *
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
// Simple class for di-charged jet analyses.
// by Beomkyu KIM
//==================================================================
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TList.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <THnSparse.h>
#include <TRandom.h>
#include "THistManager.h"
#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"
#include "AliCentrality.h"
#include "AliBSDiJetTask.h"
#include <TRandom3.h>
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliEmcalPythiaInfo.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliAnalysisUtils.h"
//#include "AliAnalysisTaskCounter.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliCalorimeterUtils.h"
#include "AliEMCALGeometry.h" 
#include "AliMultSelection.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
const Double_t pionmass = AliPID::ParticleMass(AliPID::kPion);
const Double_t pi = TMath::Pi();

AliBSDiJetTask::AliBSDiJetTask() 
  : AliAnalysisTaskEmcalJet("AliBSDiJetTask",kTRUE)
    , fOutput(0)
{
  DefineOutput (1, TList::Class());
}

  AliBSDiJetTask::AliBSDiJetTask(const char *name) 
  : AliAnalysisTaskEmcalJet(name,kTRUE)
  , fOutput(0)
{
  DefineOutput (1, TList::Class());
}
  AliBSDiJetTask::AliBSDiJetTask(const char *name, const char *option) 
  : AliAnalysisTaskEmcalJet(name,kTRUE)
  , fOutput(0)
  , fOption(option)
{
  DefineOutput (1, TList::Class());
}

  AliBSDiJetTask::AliBSDiJetTask(const AliBSDiJetTask& ap) 
  : AliAnalysisTaskEmcalJet(ap.fName,kTRUE)
  , fOutput(ap.fOutput)
{
}
AliBSDiJetTask& AliBSDiJetTask::operator = (const AliBSDiJetTask& ap)
{

  this->~AliBSDiJetTask();
  new(this) AliBSDiJetTask(ap);
  return *this;
}

AliBSDiJetTask::~AliBSDiJetTask()
{
  delete fOutput;
  delete fBSRandom;
}

void AliBSDiJetTask::UserCreateOutputObjects(){
  fBSRandom = new TRandom3;
  fBSRandom->SetSeed();
  //===============================
  // BINS
  //===============================
  auto binJetMass    = AxisFix( "JetMass", 50, 0, 100 );
  auto binInvM       = AxisVar( "InvM", {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,160,170,180,190,200,220,240,270,300,500,1000,5000});
  auto binCent       = AxisFix( "Cent", 1, 0, 100 );
  auto binDiJetSel   = AxisStr( "DiJetSel", { "NoCut"
		,"B1","B2","B3","C1","C2","C3","D1","D2","D3"
		,"M1","M2","M3","M4","M5","M6","M7" } );
  fNDiJetSelection   = binDiJetSel.GetNbins();
  if( fIsAA  ){
    binCent = AxisFix( "Cent", 20, 0, 100 );
  }

	if (fOption.Contains("LHC16l")){
		binCent = AxisVar("Cent",{0,0.001,0.01,0.1,0.5,1,5,10,15,20,30,40,50,70,100});

	}
  auto binLog1k     = AxisLog("Log1k",500,0.1,1000,0);
  //auto binLog3c     = AxisVar("Log3c",{0,5,7,9,12,16,21,28,36,45,57,70,85,99,115,132,150,169,190,212,235,500});
  auto binLog3c = AxisVar("Log3c",{0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,160,170,180,190,200,220,240,270     ,300,500,1000,5000});
  auto bintpt = AxisVar("JetTPt",{0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,160,170,180,190,200,220,240,270     ,300,500,1000,5000});
  auto bin1c = AxisFix("Fix1c",100,0,100);
  auto binAsim = AxisFix("Asim", 100 ,0,1);
  auto binM = AxisFix("BinM",600,-300,300);
	//auto binjetpt = AxisVar("binjetpt",{0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,35,40,45,50,55,60,65,70,80,90,100,120,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000});

	//auto binjetpt = AxisVar("binjetpt",{0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,500,1000});
	auto binjetpt = AxisLog("binjetpt",100,1,1000,1);;

  //===============================
  // HISTOGRAMS
  //===============================
  fHistos = new THistManager("jethists");

  const int nbins=100;
  Double_t logbins[nbins+1];
  Double_t low= 0.1;
  Double_t high=300;
  Double_t logbw= (log(high)-log(low))/nbins;
  for(int ij=0;ij<=nbins;ij++) logbins[ij]=low*exp(ij*logbw);

  fHistos -> CreateTH1("hJetMulti","",101,-1,100,"s");
  fHistos -> CreateTH1("jetpt","jetpt",nbins,logbins,"s");
  fHistos -> CreateTH1("jetptresum","jetptresum",nbins,logbins,"s");
  fHistos -> CreateTProfile("jetptresumRatio","jetptresumRatio",nbins,logbins,"s");
  fHistos -> CreateTH1("jetinv","jetinv",600,-300,300,"s");
  fHistos -> CreateTH1("jetinvcorr","jetinvcorr",600,-300,300,"s");
  fHistos -> CreateTH1("zvtx","zvtx",60,-30,30,"s");
  fHistos -> CreateTH1("mczvtx","mczvtx",60,-30,30,"s");
  fHistos -> CreateTH1("nCentrality","nCentrality", 10,0,100,"s");
  fHistos -> CreateTH2("trketaphi","trketaphi",20,-1,1,90,0,TMath::TwoPi(),"s");
  fHistos -> CreateTH2("jetetaphi","jetetaphi",20,-1,1,90,0,TMath::TwoPi(),"s");
  fHistos -> CreateTH2("trketaphifulljet","trketaphifulljet",20,-1,1,90,0,TMath::TwoPi(),"s");
  fHistos -> CreateTH2("cluetaphifulljet","cluetaphifulljet",20,-1,1,90,0,TMath::TwoPi(),"s");
	fHistos -> CreateTH1("hGammaEnergy","EMCAL gamma energy",600,0,60,"s");
	fHistos -> CreateTH1("hGammaEnergyBeIso","EMCAL gamma energy",600,0,60,"s");
	fHistos -> CreateTH1("hGammaEnergyAfIso","EMCAL gamma energy",600,0,60,"s");
	fHistos -> CreateTH2("hLambda0vsEnergy","EMCAL lambda0 vs energy",350 ,5,40,300,0,3,"s");
	fHistos -> CreateTH2("hLambda1vsEnergy","EMCAL lambda1 vs energy",350 ,5,40,300,0,3,"s");
	fHistos -> CreateTH2("hLambda0vsPt","EMCAL lambda0 vs Pt",350 ,5,40,300,0,3,"s");
	fHistos -> CreateTH2("hLambda1vsPt","EMCAL lambda1 vs Pt",350 ,5,40,300,0,3,"s");
	fHistos -> CreateTH2("hDispersionvsEnergy","EMCAL Disp vs energy",350 ,5,40,300,0,3,"s");
	fHistos -> CreateTH2("hDispersionvsPt","EMCAL Disp vs Pt",350 ,5,40,300,0,3,"s");
	fHistos -> CreateTH2("hNCellsvsEnergy","N cells vs energy", 400,0,40,30,0,30,"s");
	fHistos -> CreateTH1("hNLM","Local maxima",30,0,30,"s");
	fHistos -> CreateTH1("hXE","XE",20,0,1,"s");
	fHistos -> CreateTH1("hXEunderlying","XE",20,0,1,"s");
	fHistos -> CreateTH1("hdPhiGamma","dPhi_Gamma_hpm",100,0,2*pi,"s");
	fHistos -> CreateTH2("hIsoEtaPhi","Eta phi of iso gamma",100,0,2*pi,100,-1,1,"s");
	fHistos -> CreateTH2("hGammaJetdetadphi","",100,-1,1,100,-0.5*pi,2.5*pi,"s");
	fHistos -> CreateTH2("hJetPtLeadingSubLeading","LeadingSubleading",100,0,100,100,0,100,"s");
 	fHistos -> CreateTH1("hDijetMassTest","DijetMassTest",300,0,300,"s");
 	fHistos -> CreateTH1("hDijetDPhiTest","DijetDPhiTest",100,0,2*pi,"s");
	fHistos -> CreateTH1("hInclJetPtBefore","InclDeltaJetPt",300,-100,200,"s");
	fHistos -> CreateTH1("hInclJetPtAfter","InclDeltaJetPt",300,-100,200,"s");
  vector<TString> ent ={"All","PassPileUp","PassPileUpGoodz","GoodzNTrials","GoodzNX","NXoNTrials"};
  auto h = fHistos->CreateTH1("hEventNumbers","",ent.size(), 0, ent.size());
  for(auto i=0u;i<ent.size();i++) h->GetXaxis()->SetBinLabel(i+1,ent.at(i).Data());

  CreateTHnSparse( "hJetPtLeading","",3,{
      binDiJetSel,binCent,bintpt}, "s" );
  CreateTHnSparse( "hRho","Rho dist",2,{binCent,bin1c},"s");
  CreateTHnSparse( "hDiJetInvM", "DiJetMass", 4, {
      binDiJetSel, binCent,bintpt,binInvM }, "s" );
  CreateTHnSparse( "hDiJetDPhi_0_2pi", "DiJet #Delta#Phi", 5, {
      binDiJetSel, binCent, bintpt,binInvM, AxisFix( "",100, 0, 2*pi ) },"s");
  CreateTHnSparse( "hDiJetPtPair", "DiJet PtPair", 5, {
      binDiJetSel, binCent, bintpt, binInvM, binLog3c  },"s");
  CreateTHnSparse( "hDiJetPtPairH", "DiJet PtPair", 5, {
      binDiJetSel, binCent, bintpt,binInvM, binLog3c  },"s");
  CreateTHnSparse( "hDiJetPtPairL", "DiJet PtPair", 5, {
      binDiJetSel, binCent, bintpt,binInvM, binLog3c  },"s");
  CreateTHnSparse( "hDiJetInvMPtPairRes", "DiJet InvM PtPair Res Matrix", 7, {
      binDiJetSel, binCent, bintpt,binInvM,binInvM,binLog3c,binLog3c},"s");
  CreateTHnSparse( "hDiJetInvMPtPairMiss", "DiJet InvM PtPair missing tracks", 5, {
      binDiJetSel, binCent, bintpt,binInvM, binLog3c},"s");
  CreateTHnSparse( "hDiJetInvMPtPairFake", "DiJet InvM PtPair fake tracks", 5, {
      binDiJetSel, binCent, bintpt,binInvM, binLog3c},"s");
 
	
  CreateTHnSparse( "hInclJetPt", "Inclusive jet pt", 2, {
      binCent, binjetpt  },"s");
  CreateTHnSparse( "hInclJetPtCorrectedScaledToMB", "Inclusive jet pt", 2, {
      binCent, binjetpt  },"s");
  CreateTHnSparse( "hInclJetPtCorrectedScaledToMBH", "Inclusive jet pt", 2, {
      binCent, binjetpt  },"s");
  CreateTHnSparse( "hInclJetPtCorrectedScaledToMBL", "Inclusive jet pt", 2, {
      binCent, binjetpt  },"s");
  CreateTHnSparse( "hInclJetPtScaledToMB", "Inclusive jet pt", 2, {
      binCent, binjetpt  },"s");
  CreateTHnSparse( "hInclJetPtCorrected", "Inclusive jet pt", 2, {
      binCent, binjetpt  },"s");
  CreateTHnSparse( "hInclJetPtCorrectedAndMatched", "Inclusive jet pt", 2, {
      binCent, binjetpt  },"s");
  CreateTHnSparse( "hInclJetMass", "Inclusive jet Mass", 2, {
      binCent, binjetpt  },"s");
  CreateTHnSparse( "hInclJetMassScaledToMB", "Inclusive jet Mass", 2, {
      binCent, binjetpt  },"s");
  CreateTHnSparse( "hInclJetMassCorrected", "Inclusive jet mass", 2, {
      binCent, binjetpt  },"s");
  CreateTHnSparse( "hInclJetMassCorrectedScaledToMB", "Inclusive jet mass", 2, {
      binCent, binjetpt  },"s");

  CreateTHnSparse( "hInclJetPtRes", "Inclusive jet pt Res Matrix", 3, {
      binCent, binjetpt,binjetpt},"s");
  CreateTHnSparse( "hInclJetPtFake", "Inclusive jet pt fake tracks", 2, {
      binCent, binjetpt},"s");
  CreateTHnSparse( "hInclJetPtMiss", "Inclusive jet pt  missing tracks", 2, {
      binCent, binjetpt},"s");


  PostData(1, fHistos->GetListOfHistograms());

  //===============================
  // For Sure
  //===============================
  std::cout<< "DEBUG4 IsAA?"<< (fIsAA?"AA":"pp")<<std::endl;
  std::cout<<"NBins of Cent : "<<binCent.GetNbins()<<"\t"<<binCent.GetXmin()<<"\t"<<binCent.GetXmax()<<endl;
  if( fNDiJetSelection != kBDiJetSelEnd-1 ){
    cout<<"fNDiJetSelection("<<fNDiJetSelection
      <<") is not match with kBDiJetSelEnd("<<kBDiJetSelEnd<<")"<<endl;
    gSystem->Exit(1);
  }
  for (auto i=0u; i<kBDiJetSelEnd; i++) {
		fDijetSelectionCut.push_back(false);
		fDijetPtPair.push_back(0.);
		fDijetInvM.push_back(0.);
  }
  fJets.clear();
  fUtils = new AliAnalysisUtils();

  fUtils -> 	SetMaxVtxZ(10);
	fCaloUtils = new AliCalorimeterUtils();
	fCaloUtils->SetNumberOfCellsFromEMCALBorder(1);
	fCaloUtils->SetNumberOfCellsFromPHOSBorder (2);
	fCaloUtils->SetNumberOfSuperModulesUsed(10);
	const Int_t year = 2011;
	if     (year == 2010) fCaloUtils->SetNumberOfSuperModulesUsed(4);
	else if(year <= 2013) fCaloUtils->SetNumberOfSuperModulesUsed(10);
	else if(year >  2013) fCaloUtils->SetNumberOfSuperModulesUsed(20);
	else                  fCaloUtils->SetNumberOfSuperModulesUsed(10);
	if(!fIsAA)
	{
		fCaloUtils->SetLocalMaximaCutE(0.1);
		fCaloUtils->SetLocalMaximaCutEDiff(0.03);
	}
	else
	{
		fCaloUtils->SetLocalMaximaCutE(0.2);
		fCaloUtils->SetLocalMaximaCutEDiff(0.03);
	}
	fCaloUtils->SwitchOffRecalculateClusterTrackMatching();
	fCaloUtils->SwitchOffBadChannelsRemoval() ;
	if(!fIsMC)
		fCaloUtils->SwitchOnLoadOwnEMCALGeometryMatrices();
	
	fCaloUtils->InitEMCALGeometry();

	tsf = new TF1("tsf","[0]/pow([1]+x*x*x*x,[2])+[3]",5,1000) ;
	tsfh = new TF1("tsfh","[0]/pow([1]+x*x*x*x,[2])+[3]",5,1000) ;
	tsfl = new TF1("tsfl","[0]/pow([1]+x*x*x*x,[2])+[3]",5,1000) ;


	if (fOption.Contains("LHC12d")) { //ch jet
		if (fOption.Contains("ScaleTo5TeV")){
			tsf->SetParameters(20.8877,197.93,0.813918,0.000309556);
			tsfh->SetParameters(21.9551,207.686,0.818395,0.000332579);
			tsfl->SetParameters(19.864,188.184,0.809403,0.000286503);
		}
		else {
			tsf->SetParameters(27.9874,419.65,0.800532,0.000613293);
			tsfh->SetParameters(29.7364,436.329,0.806079,0.000655603);
			tsfl->SetParameters(26.3315,403.03,0.794951,0.000570913);
		}
	} 
	else if (fOption.Contains("LHC12f")){
		if (fOption.Contains("ScaleTo5TeV")){
			tsf->SetParameters(21.768,224.487,0.811922,0.000256762);
			tsfh->SetParameters(22.937,233.167,0.816485,0.00027829);
			tsfl->SetParameters(20.6467,215.804,0.807307,0.000235205);
		} else {
			tsf->SetParameters(29.2384,450.817,0.798505,0.000520484);
			tsfh->SetParameters(30.9866,465.512,0.803686,0.000557798);
			tsfl->SetParameters(27.5756,436.145,0.793281,0.000483107);
		}
	}
	else if (fOption.Contains("LHC12h")){
		if (fOption.Contains("ScaleTo5TeV")){
			if (fOption.Contains("Trigger1")){
				tsf->SetParameters(23.2189,216.902,0.823087,0.000208479);
				tsfh->SetParameters(24.2083,222.941,0.826666,0.000224221);
				tsfl->SetParameters(22.2599,210.825,0.819468,0.000192719);
			}
			if (fOption.Contains("Trigger2")){
				tsf->SetParameters(17.5794,91.6111,0.794505,0.00028183);
				tsfh->SetParameters(18.4342,99.0281,0.798727,0.000303145);
				tsfl->SetParameters(16.7568,84.1957,0.790243,0.000260487);
			}

		} else { 
			if (fOption.Contains("Trigger1")){
				tsf->SetParameters(30.7108,432.104,0.807976,0.000424779);
				tsfh->SetParameters(32.3394,444.05,0.812517,0.000455315);
				tsfl->SetParameters(29.1507,420.144,0.803393,0.000394201);
			}
			if (fOption.Contains("Trigger2")){
				tsf->SetParameters(22.888,277.488,0.778264,0.000561091);
				tsfh->SetParameters(24.2841,291.415,0.78364,0.000602462);
				tsfl->SetParameters(21.563,263.625,0.772848,0.000519645);
			}
			if (fOption.Contains("All")){
				tsf->SetParameters(26.2369,880.724,0.736723,0.000845877);
				tsfh->SetParameters(28.0387,913.16,0.742835,0.000932078);
				tsfl->SetParameters(24.5484,848.55,0.7306,0.000759487);
			}
		}
	}
	else if (fOption.Contains("LHC13d")) {
		tsf->SetParameters(250.954,1918.07,0.865423,0.000799134);
		tsfh->SetParameters(266.475,1932.12,0.870167,0.000851788);
		tsfl->SetParameters(236.116,1903.93,0.860602,0.000746376);
	}
 	else if (fOption.Contains("LHC13e")) {
		tsf->SetParameters(199.456,1931.83,0.838126,0.00114861);
		tsfh->SetParameters(212.014,1955.28,0.843172,0.00122398);
		tsfl->SetParameters(187.52,1908.36,0.833025,0.00107308);
	}
	else if (fOption.Contains("LHC13f")) {
		tsf->SetParameters(143.084,1774.16,0.800964,0.00128229);
		tsfh->SetParameters(152.963,1803.91,0.806621,0.00139113);
		tsfl->SetParameters(133.764,1744.48,0.795255,0.0011732);
	} 

	else {
		tsf->SetParameters(37.8667,490.515,0.795349,0.000840197);
		tsfh->SetParameters(41.4233,520.858,0.803635,0.00094383);
		tsfl->SetParameters(34.594,460.346,0.787003,0.000736331);
	}

	
	/*if (fName.Contains("FullJet")){
		tsf->SetParameter(0,263.817);
		tsf->SetParameter(1,4146.28);
		tsf->SetParameter(2,1.01805);
		tsf->SetParameter(3,0.000403723);

	} 
  else if (fOption.Contains("ScaleTo5")){
		tsf->SetParameter(0,22.3377);
		tsf->SetParameter(1,401.885);
		tsf->SetParameter(2,0.777928);
		tsf->SetParameter(3,0.000394575);
	}
	else if (fOption.Contains("ScaleTo13bc") && fOption.Contains("13d")){
		tsf->SetParameter(0,248.968);
		tsf->SetParameter(1,1901.53);
		tsf->SetParameter(2,0.864777);
		tsf->SetParameter(3,0.000799705);
	}
	else if (fOption.Contains("ScaleTo13bc") && fOption.Contains("13e")){
		tsf->SetParameter(0,199.376);
		tsf->SetParameter(1,1927.48);
		tsf->SetParameter(2,0.838168);
		tsf->SetParameter(3,0.00115039);
	}

	else if (fOption.Contains("ScaleTo13bc") && fOption.Contains("13f")){
		tsf->SetParameter(0,143.121);
		tsf->SetParameter(1,1771.21);
		tsf->SetParameter(2,0.801047);
		tsf->SetParameter(3,0.00128567);
	}*/


}

//________________________________________________________________________

Bool_t AliBSDiJetTask::Run(){
	using TMath::Abs;
  fCent = -1;
	//Cent = InputEvent()->GetCentrality();
	sel = (AliMultSelection*) InputEvent() -> FindListObject("MultSelection");
	if (sel) {
		if (fOption.Contains("13f")) fCent = sel -> GetMultiplicityPercentile("V0C");
		else if (
				fOption.Contains("13b") ||
				fOption.Contains("13c") ||
				fOption.Contains("13d") ||
				fOption.Contains("13e") 
				) fCent = sel -> GetMultiplicityPercentile("V0A");
		else if (fOption.Contains("15o")) fCent = sel -> GetMultiplicityPercentile("V0M");
		else if (fOption.Contains("16l")) fCent = sel -> GetMultiplicityPercentile("V0M");
		else if (fOption.Contains("17n")) fCent = sel -> GetMultiplicityPercentile("V0M");
		else	fCent = 50.;
	}

  //pt hard bin scaling-----------------------------------------------------------
	//https://twiki.cern.ch/twiki/bin/viewauth/ALICE/JetMCProductionsCrossSections
  Double_t sf = 1.;
	Int_t NTrials = -1;
	Double_t XSection =-1;
  TLorentzVector p6,p7,p67;
  Double_t genzvtx = -30;
	if(fIsMC){
		if (fOption.Contains("AOD")){
			AliVEvent *event = InputEvent();
			if (!event) {
				Printf("ERROR: Could not retrieve event");
				return kFALSE;
			}
			AliAODMCHeader *cHeaderAOD  = dynamic_cast<AliAODMCHeader*>
				(event->FindListObject(AliAODMCHeader::StdBranchName()));
      genzvtx = cHeaderAOD -> GetVtxZ();
      fHistos->FillTH1("mczvtx",genzvtx);
			TList *genHeaders         = cHeaderAOD->GetCocktailHeaders();
			NTrials = -1;
			XSection =-1;
			AliGenEventHeader* gh       = 0;
			for(Int_t i = 0; i<genHeaders->GetEntries();i++){
				gh   = (AliGenEventHeader*)genHeaders->At(i);
				TString GeneratorName   = gh->GetName();
				if (GeneratorName.CompareTo("Pythia") == 0){
					AliGenPythiaEventHeader* gPythia = dynamic_cast<AliGenPythiaEventHeader*>(gh);
					NTrials = gPythia->Trials();
					XSection = gPythia->GetXsection();
				}
			}
			sf = XSection/NTrials;

			fMCArray = (TClonesArray*) event->FindListObject("mcparticles");
			const Int_t nTracksMC = fMCArray->GetEntriesFast();
			for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
				AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
				Int_t pdgCode = trackMC->PdgCode();
				//cout<<"iTracks  : "<<iTracks<<"   pdgcode : "<<pdgCode
				//<< "  Pt " <<trackMC->Pt()
				//<< "   Phi  " <<trackMC->Phi()
				//<<endl;
				if (iTracks == 6) p6.SetXYZT(trackMC->Px(), trackMC->Py(), trackMC->Pz(), trackMC->E());
				if (iTracks == 7) p7.SetXYZT(trackMC->Px(), trackMC->Py(), trackMC->Pz(), trackMC->E());
				// Select only primaries
				//if(!(trackMC->IsPhysicalPrimary())) continue;
				//Float_t ptMC = trackMC->Pt();
			}	
		}
    else {
			AliMCEvent *mcEvent = MCEvent();
			AliStack *stack = mcEvent->Stack();
      AliGenPythiaEventHeader*  gPythia =
				dynamic_cast<AliGenPythiaEventHeader*>(mcEvent->GenEventHeader());
			NTrials = gPythia->Trials();
			XSection = gPythia->GetXsection();
			sf = XSection/NTrials;
      Int_t nPrim  = stack->GetNprimary();
			for (Int_t i = 0; i < nPrim; i++){
				TParticle* p = stack->Particle(i);
        if (!p) continue;
			}
		}
	}
  //The method above doesn't work for MC 13b4_fix and plus AOD files
  //For these MC productions, the method below is used..
	if ( fOption.Contains("13plus") 
    || fOption.Contains("12a15e")
    || fOption.Contains("12a15f")
    || fOption.Contains("13e4")
		|| fOption.Contains("13fix") ){
		TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
		if(!tree) return kFALSE;
		TFile *curfile = tree->GetCurrentFile();
		if(!curfile) return kFALSE;
		TString fCurrFileName = TString(curfile->GetName());
		if (!filename.EqualTo(fCurrFileName)) {
      filename = fCurrFileName.Data();
			fCurrFileName.ReplaceAll(gSystem->BaseName(fCurrFileName.Data()),"");
			TFile *fxsec = TFile::Open(Form("%s%s",fCurrFileName.Data()
						,"pyxsec_hists.root"));
			TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
			TList *list = dynamic_cast<TList*>(key->ReadObj());
			XSection = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
			NTrials = ((TH1F * ) list->FindObject("h1Trials"))->GetBinContent(1);
			Int_t entries = ((TH1F * ) list->FindObject("h1Trials"))->GetEntries();
			NTrials /= entries;
			pyxsechistsf = XSection/NTrials;
      fxsec->Close();
		}
		sf = pyxsechistsf;
	}

  //pT hard bin scaling end---------------------------------------------------

  
  fHistos->FillTH1("hEventNumbers","All",1);
	AliVEvent *event = InputEvent();
	event->IsA()==AliESDEvent::Class()
		? fEvt = dynamic_cast<AliESDEvent*>(event)
		: fEvt = dynamic_cast<AliAODEvent*>(event);
	if (!fEvt) return false;


  if (!fOption.Contains("MC") && fUtils->IsPileUpSPD(InputEvent())) {
    PostData(1, fHistos->GetListOfHistograms());
    return false;
  }
  
	fHistos->FillTH1("hEventNumbers","PassPileUp",1);

  // Fill z_vertex and cut z_vertex range
  Bool_t IsGoodVertex = false;
  IsGenGoodVtx = false;
	Double_t vertex[3] = {0,0,0};             
  if (fName.Contains("kinefull") || fName.Contains("kinech")){
		fHistos -> FillTH1 ("zvtx",genzvtx);
  	if (abs(genzvtx)<=10) {
			IsGoodVertex = true;  
      IsGenGoodVtx = true;
    }
  } else {
  	const AliVVertex* vtx = InputEvent()->GetPrimaryVertex();
    if (fUtils->IsVertexSelected2013pA(InputEvent())){
    //if (vtx->	GetNContributors()>1){
			double zvtx = vtx->GetZ();
			vtx ->GetXYZ(vertex);
			fHistos -> FillTH1 ("zvtx",zvtx);
      if (abs(zvtx)<=10) IsGoodVertex = true; 
    } 
	}
  
  
	if (IsGoodVertex) {
		fHistos->FillTH1("hEventNumbers","PassPileUpGoodz",1);
		fHistos->FillTH1("hEventNumbers","GoodzNTrials",NTrials);
		fHistos->FillTH1("hEventNumbers","GoodzNX",XSection);
		fHistos->FillTH1("hEventNumbers","NXoNTrials",XSection/NTrials);
		fHistos -> FillTH1 ("nCentrality",fCent);
  }
	
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	AliBSDiJetTask *task = nullptr;
  if (fName.Contains("FullJet"))  task = (AliBSDiJetTask*)man -> GetTask("kinefullFullJet");
	else task = (AliBSDiJetTask*)man -> GetTask("kinech");
	//-----------------------------------



  //===============================
  // RESUM JETS
  //===============================

  auto jetContainer = GetJetContainer(0);
  auto trkContainer = GetParticleContainer(0);
	AliClusterContainer* cluContainer = nullptr;
	if (fName.Contains("FullJet")) cluContainer = GetClusterContainer(0);
  auto ktContainer = GetJetContainer(1);

  // Underlying background calculation
  // this will be done by four-momentum sum of constituents
  
  Double_t RHO=0;
	Double_t RHOM=0;
	if (!fName.Contains("kinefull") && !fName.Contains("kinech")){
    int n = 0;
		TLorentzVector1D rhoarray;
		Double1D Sumpt;
		Double1D Summ;
		for( auto j : ktContainer->accepted() ){
      if (Abs(j->Eta())>0.7)  continue;
			double lpt = 0;
      double sumpt = 0;
      double summ = 0;

			for( int it=0; it<j->GetNumberOfTracks(); it++ ) {
				auto trk =  j->Track(it);
				TLorentzVector temp;
				temp.SetXYZM(trk->Px(),trk->Py(),trk->Pz(),pionmass);
        if( lpt < temp.Pt() )  lpt = temp.Pt();
        sumpt += temp.Pt();
        summ += (sqrt(pionmass*pionmass+temp.Pt()*temp.Pt())-temp.Pt());
			}
      sumpt /= j->Area();
      summ /= j->Area();
      //if( lpt < fLeadingParticlePtMin ) continue;
			if (n>=2) {//remove leading and subleading jets
				Sumpt.push_back(sumpt);
				Summ.push_back(summ);
			}
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
	}


  // =============================
  // QA plots
  // ============================
  for (auto trk : trkContainer->accepted()){
    fHistos->FillTH2("trketaphi",trk->Eta(),trk->Phi());
  }
  //cout<<((AliAODEvent*)AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalEvent())->GetCentrality()->GetCentralityPercentile("V0M")<<endl;


  //=== NOTE : We will not use sj[0] because bin in THnSparse is always begin with 1.
  TLorentzVector2D sj( fNDiJetSelection+1, TLorentzVector1D(2));
  TLorentzVector1D sumJets, sumJetsCorr;
  std::vector<Double_t> lpts;
  fJets.clear();
  int ij=0;
  FillTHnSparse("hRho",{fCent,RHO},sf); 
  double rho = 0;
  rho = jetContainer->GetRhoVal();
  for( auto j : jetContainer->accepted() ){
		fHistos->FillTH2("jetetaphi",j->Eta(),j->Phi(),sf);
    TLorentzVector sum (0,0,0,0);
    double sumpt=0;
    double lpt = 0;
    //cout<<"jet pt : "<<j->Pt()<<"  total tracks : "<<j->GetNumberOfTracks()<<endl;

    for( int it=0; it<j->GetNumberOfTracks(); it++ ) {
      auto trk =  j->Track(it);
		 	//cout<<"trk eta : "<<trk->Eta()<<"    trk phi : "<<trk->Phi()<< "   Label "<<trk->GetLabel()<<endl;
      if (fName.Contains("MBTR")){
				if (fBSRandom->Uniform(0,100)<5.) continue;
      }
      TLorentzVector temp;
      temp.SetXYZM(trk->Px(),trk->Py(),trk->Pz(),pionmass);
			fHistos->FillTH2("trketaphifulljet",trk->Eta(),trk->Phi(),sf);
			//if (fName.Contains("FullJet")) cout<<"track Pt : "<<trk->Pt()<<endl;
      sum+=temp;
      sumpt+=trk->Pt();
      if( lpt < temp.Pt() )  lpt = temp.Pt();
      
    }
		if (fName.Contains("FullJet")){
			for( int it=0; it<j->GetNumberOfClusters(); it++ ) {
				auto clu =  j->Cluster(it);
				TLorentzVector temp;
				clu->GetMomentum(temp,vertex);
				fHistos->FillTH2("cluetaphifulljet",temp.Eta(),temp.Phi(),sf);
				//cout<<"Clu Pt : "<<temp.Pt()<<endl;
				sum+=temp;
				sumpt+=temp.Pt();
				if( lpt < temp.Pt() )  lpt = temp.Pt();

			}
		}

    TLorentzVector avec(0,0,0,0);
    avec.SetPtEtaPhiE(j->AreaPt(),j->AreaEta(),j->AreaPhi(),j->AreaE());
    TLorentzVector sumcorr(
			  sum.Px()-RHO*avec.Px()
			, sum.Py()-RHO*avec.Py()
			, sum.Pz()-(RHO+RHOM)*avec.Pz()
			, sum.E() -(RHO+RHOM)*avec.E()
		);
		//MC Pythia jet / pT-hard > 4 cut
		if( fOption.Contains("MC")) {
			if (TMath::Abs(sum.DeltaR(p6))<0.4 && (sumpt/p6.Pt())>4) return false;
			if (TMath::Abs(sum.DeltaR(p7))<0.4 && (sumpt/p7.Pt())>4) return false;
		}
		
    //=== Jet geometric and kinematic CUT 
		//=== MCTRUTH=0, MCREC & DATA = 5GeV
    if( lpt < fLeadingParticlePtMin ) continue;
		if (IsGoodVertex && abs(sum.Eta())<0.5){
			fHistos->FillTH1("hInclJetPtBefore",sum.Pt());
		}
		if (IsGoodVertex && abs(sumcorr.Eta())<0.5){
			fHistos->FillTH1("hInclJetPtAfter",sumcorr.Pt());
		}

		if (IsGoodVertex){
			fHistos->FillTH1("jetptresum",sumpt,sf);
			fHistos->FillProfile("jetptresumRatio",sumpt, j->Pt()/sumpt);
			Double_t tratio = 1.;
			Double_t tratioh = 1.;
			Double_t tratiol = 1.;
			if (fName.Contains("EMCEJE")) {
				tratio = tsf->Eval(sum.Pt());
				tratioh = tsfh->Eval(sum.Pt());
				tratiol = tsfl->Eval(sum.Pt());
			}
			if (abs(sum.Eta())<0.5){
				if (fOption.Contains("LHC16j5") && (sum.DeltaR(p6)<0.1 || sum.DeltaR(p7)<0.1))  
					FillTHnSparse("hInclJetPt",{fCent,sum.Pt()},sf);
				else FillTHnSparse("hInclJetPt",{fCent,sum.Pt()},sf);
				FillTHnSparse("hInclJetPtScaledToMB",{fCent,sum.Pt()},sf*tratio);
				FillTHnSparse("hInclJetMass",{fCent,sum.M()},sf);
				FillTHnSparse("hInclJetMassScaledToMB",{fCent,sum.M()},sf*tratio);
			}	
			if (abs(sumcorr.Eta())<0.5){
				//if (fOption.Contains("LHC16j5") && (sumcorr.DeltaR(p6)<0.1 || sumcorr.DeltaR(p7)<0.1)) 
				FillTHnSparse("hInclJetPtCorrected",{fCent,sumcorr.Pt()},sf);
				FillTHnSparse("hInclJetPtCorrectedScaledToMB",{fCent,sumcorr.Pt()},sf*tratio);
				FillTHnSparse("hInclJetPtCorrectedScaledToMBH",{fCent,sumcorr.Pt()},sf*tratioh);
				FillTHnSparse("hInclJetPtCorrectedScaledToMBL",{fCent,sumcorr.Pt()},sf*tratiol);
				FillTHnSparse("hInclJetMassCorrected",{fCent,sumcorr.M()},sf);
				FillTHnSparse("hInclJetMassCorrectedScaledToMB",{fCent,sumcorr.M()},sf*tratio);
			}
		}

    if( Abs(sum.Eta()) < 0.5 ) sumJets.push_back(sum);
		if( Abs(sumcorr.Eta()) < 0.5) sumJetsCorr.push_back(sumcorr);
		// Using corrected jets for any kind of reconstructed jets
		if (!fName.Contains("kinech")) sum = sumcorr; 
    lpts.push_back(lpt);

    //=== DiJetSelection 1 : LS = Leading-SubLeading
    if (sj[0][0].Pt() < sum.Pt()) {
 			sj[0][0]=sum;
    }
    else if( sj[0][1].Pt() < sum.Pt()){
			sj[0][1]=sum;
    }
    ij++;
  }
	if (IsGoodVertex) fHistos->FillTH2("hJetPtLeadingSubLeading",sj[0][0].Pt(),sj[0][1].Pt(),sf);

	if (IsGoodVertex){
		if (sj[0][0].Pt()>20 && sj[0][1].Pt()>20) fHistos->FillTH1("hDijetMassTest",(sj[0][0]+sj[0][1]).M(),sf);
		if (sj[0][0].Pt()>20 && sj[0][1].Pt()>20) fHistos->FillTH1("hDijetDPhiTest",TVector2::Phi_0_2pi(sj[0][0].DeltaPhi(sj[0][1])),sf);
	}
	// chjet isolated jet correlation



	// Using corrected jets for any kind of reconstructed jets
	if (!fName.Contains("kinech")) sumJets = sumJetsCorr;   



	TLorentzVector1D matchedjets;
	if (task && task->GetIsGenGoodVtx()){ // Inclusive jet pt response matrix
		TLorentzVector1D pjets = task->GetJets();
		TLorentzVector1D rjets = sumJets;
		if (IsGoodVertex){
			for (auto pj : pjets){
				double maxjpt = 0;
				TLorentzVector maxjet(0,0,0,0);
				for (auto rj : rjets){
					//if (rj.Pt()>pj.Pt()*0.5 && maxjet.Pt()<rj.Pt() && rj.DeltaR(pj)<0.1)  {
					if (rj.Pt()>pj.Pt()*0.5 && rj.Pt()<pj.Pt()*1.5 && maxjet.Pt()<rj.Pt())  {
						maxjet = rj;
					}
				}
				if (maxjet.Pt()>0) {
					matchedjets.push_back(maxjet);
					FillTHnSparse("hInclJetPtCorrectedAndMatched",{fCent,maxjet.Pt()},sf);	
					FillTHnSparse("hInclJetPtRes",{fCent,maxjet.Pt(),pj.Pt()},sf);
				}
				if (maxjet.Pt()==0) FillTHnSparse("hInclJetPtMiss",{fCent,pj.Pt()},sf);
				}
				for (auto rj : rjets){
					bool matched = false;
					for (auto mj : matchedjets){
						if (rj.Pt() == mj.Pt()) matched = true;
					}
					if (!matched ) FillTHnSparse("hInclJetPtFake",{fCent,rj.Pt()},sf); 
				}
			}
			else {
				for (auto pj : pjets){
					FillTHnSparse("hInclJetPtMiss",{fCent,pj.Pt()},sf);
				}
			}
		} 

		// do test
		//cout<<"sumJets"<<endl;
		//for (auto j : sumJets)	cout<<j.Pt()<<endl;
		//cout<<"sumJetCorr"<<endl;
		//for (auto j : sumJetsCorr)	cout<<j.Pt()<<endl;
		//cout<<"matched jets"<<endl;
		//for (auto j : matchedjets)	cout<<j.Pt()<<endl;

		//cout<<"sj[0] before "<<endl;
		//cout<<sj[0][0].Pt()<<endl;
		//cout<<sj[0][1].Pt()<<endl;
	
		// Re-searching leading and subleading for matched jets. 
		if (fOption.Contains("16j5") && !fName.Contains("kinech") ){
			sj[0][0] = TLorentzVector(0,0,0,0);
			sj[0][1] = TLorentzVector(0,0,0,0);
			for (auto j : matchedjets){
				if (sj[0][0].Pt() < j.Pt()) {
					sj[0][0]=j;
				}
				else if( sj[0][1].Pt() < j.Pt()){
					sj[0][1]=j;
				}
			}
			//sumjets are redirected to matched jets for embedded data samples
			sumJets = matchedjets;
		}

		//cout<<"sj[0] after "<<endl;
		//cout<<sj[0][0].Pt()<<endl;
		//cout<<sj[0][1].Pt()<<endl;





	/*
  if (task && IsGoodVertex){ // Find a probe jet
    TLorentzVector1D pjets = task->GetJets();
    TLorentzVector1D rjets = sumJets;
    TLorentzVector1D rjetscorr = sumJetsCorr;
		for (auto pj : pjets){
 			double maxjpt = 0;
      TLorentzVector maxjet(0,0,0,0);
			for (auto rj : rjets){
				if (rj.Pt()>pj.Pt()/2. && maxjpt<rj.Pt() && rj.DeltaR(pj)<0.1)	{
					maxjet = rj;
				}	
			}
      rjets.erase(std::remove_if(rjets.begin(), rjets.end(),
						[&](const TLorentzVector& x) { return (x.Pt() > pj.Pt()/2)  && (x.DeltaR(pj)<0.1) ; }), rjets.end());
			
		  if (maxjet.Pt()>0) {
      }
   
			maxjpt = 0;
			maxjet.SetXYZT(0,0,0,0);
			for (auto rj : rjetscorr){
				if (rj.Pt()>pj.Pt()/2. && maxjpt<rj.Pt() && rj.DeltaR(pj)<0.1) {
					maxjet = rj;
				}
			}
			rjetscorr.erase(std::remove_if(rjetscorr.begin(), rjetscorr.end(),
						[&](const TLorentzVector& x)  { return (x.Pt() > pj.Pt()/2)  && (x.DeltaR(pj)<0.1); }), rjetscorr.end());

			if (maxjet.Pt()>0) {
			}
		}
	}
  */

  //cout<<"\n\n"<<endl; 
  for (auto i=0u; i<kBDiJetSelEnd; i++) {
    fDijetPtPair[i] = 0;
    fDijetInvM[i] = 0;
		fDijetSelectionCut[i] = false;
  }
  fJets = sumJets;

  //0 lorentz vector
  TLorentzVector1D zsj(2,TLorentzVector(0,0,0,0));

  fDijetSelectionCut[kNoCut] = true;
  sj[kNoCut] = sj[0];
	Double_t kinecut = 20;
	
  //DiJetSelection 1 : TpTJET>20, ApTJET>20, pT_pair
	if (sj[0][0].Pt()>kinecut && sj[0][1].Pt()>kinecut) {
		sj[kB1] =sj[0]; 
		fDijetSelectionCut[kB1] = true;
	}
  //DiJetSelection 2 : TpTJET>20, ApTJET>20, pT_pair 
  //Leading-Subreading & dPhi>pi/2
	if(sj[0][0].Pt()>kinecut && sj[0][1].Pt()>kinecut &&
  	Abs(sj[0][0].DeltaPhi(sj[0][1])) > pi/2. ){ 
    sj[kB2] = sj[0]; fDijetSelectionCut[kB2] = true;
	}

	//DiJetSelection 3 : TpTJET>20, ApTJET>20, pT_pair 
	//Leading - Subleading in opposite hemisphere
	if(sj[0][0].Pt()>kinecut ){
		auto tj = sj[kB3][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kB3][1] = TLorentzVector(0,0,0,0);
		for( auto jet: sumJets ){
			if( Abs(tj.DeltaPhi( jet )) < pi/2. ) continue; 
			if( sj[kB3][1].Pt() < jet.Pt() ) sj[kB3][1] = jet;
		}
    if (sj[kB3][1].Pt() >kinecut) fDijetSelectionCut[kB3] = true;
	}

	if (sj[0][0].Pt()>30 && sj[0][1].Pt()>30) {
		sj[kC1] =sj[0];
		fDijetSelectionCut[kC1] = true;
	}
	//DiJetSelection 2 : TpTJET>20, ApTJET>20, pT_pair
	//Leading-Subreading & dPhi>pi/2
	if(sj[0][0].Pt()>30 && sj[0][1].Pt()>30 &&
			Abs(sj[0][0].DeltaPhi(sj[0][1])) > pi/2. ){
		sj[kC2] = sj[0]; fDijetSelectionCut[kC2] = true;
	}

	//DiJetSelection 3 : TpTJET>20, ApTJET>20, pT_pair
	//Leading - Subleading in opposite hemisphere
	if(sj[0][0].Pt()>30 ){
		auto tj = sj[kC3][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kC3][1] = TLorentzVector(0,0,0,0);
		for( auto jet: sumJets ){
			if( Abs(tj.DeltaPhi( jet )) < pi/2. ) continue;
			if( sj[kC3][1].Pt() < jet.Pt() ) sj[kC3][1] = jet;
		}
		if (sj[kC3][1].Pt() >30) fDijetSelectionCut[kC3] = true;
	} 

	if (sj[0][0].Pt()>40 && sj[0][1].Pt()>40) {
		sj[kD1] =sj[0];
		fDijetSelectionCut[kD1] = true;
	}
	//DiJetSelection 2 : TpTJET>20, ApTJET>20, pT_pair
	//Leading-Subreading & dPhi>pi/2
	if(sj[0][0].Pt()>40 && sj[0][1].Pt()>40 &&
			Abs(sj[0][0].DeltaPhi(sj[0][1])) > pi/2. ){
		sj[kD2] = sj[0]; fDijetSelectionCut[kD2] = true;
	}

	//DiJetSelection 3 : TpTJET>20, ApTJET>20, pT_pair
	//Leading - Subleading in opposite hemisphere
	if(sj[0][0].Pt()>30 ){
		auto tj = sj[kD3][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kD3][1] = TLorentzVector(0,0,0,0);
		for( auto jet: sumJets ){
			if( Abs(tj.DeltaPhi( jet )) < pi/2. ) continue;
			if( sj[kD3][1].Pt() < jet.Pt() ) sj[kD3][1] = jet;
		}
		if (sj[kD3][1].Pt() >30) fDijetSelectionCut[kD3] = true;
	} 


  //DiJetSelection 4 : 20<TpTJET<40, 20<ApTJET, kTy
	if (sj[0][0].Pt()>20 && sj[0][0].Pt()<40 && sj[0][1].Pt()>20 &&
			Abs(TVector2::Phi_0_2pi(sj[0][0].DeltaPhi(sj[0][1]))-pi) < pi/3.){ 	
		sj[kM1] =sj[0]; fDijetSelectionCut[kM1] = true;}
  else sj[kM1] = zsj;

  //DiJetSelection 5 : 40<TpTJET<60, 20<ApTJET, kTy
	if (sj[0][0].Pt()>40 && sj[0][0].Pt()<60 && sj[0][1].Pt()>20 &&
			Abs(TVector2::Phi_0_2pi(sj[0][0].DeltaPhi(sj[0][1]))-pi) < pi/3.){	
		sj[kM2] =sj[0];fDijetSelectionCut[kM2] = true;}
  else sj[kM2] = zsj;

  //DiJetSelection 6 : 60<TpTJET<80, 20<ApTJET, kTy
	if (sj[0][0].Pt()>60 && sj[0][0].Pt()<80 && sj[0][1].Pt()>20 &&	
			Abs(TVector2::Phi_0_2pi(sj[0][0].DeltaPhi(sj[0][1]))-pi) < pi/3.){	
		sj[kM3] =sj[0]; fDijetSelectionCut[kM3] = true;}
  else sj[kM3] = zsj;

  //DiJetSelection 7 : 80<TpTJET<110, 20<ApTJET, kTy
	if (sj[0][0].Pt()>80 && sj[0][0].Pt()<110 && sj[0][1].Pt()>20 && 	
			Abs(TVector2::Phi_0_2pi(sj[0][0].DeltaPhi(sj[0][1]))-pi) < pi/3.){	
		sj[kM4] =sj[0];fDijetSelectionCut[kM4] = true;}
  else sj[kM4] = zsj;

  //DiJetSelection 8 : 70<TpTJET<110, 20<ApTJET<30, kTy
	if (sj[0][0].Pt()>70 && sj[0][0].Pt()<110 &&
			sj[0][1].Pt()>20 && sj[0][1].Pt()<30 && 	
			Abs(TVector2::Phi_0_2pi(sj[0][0].DeltaPhi(sj[0][1]))-pi) < pi/3.){	
		sj[kM5] =sj[0];fDijetSelectionCut[kM5] = true;}
  else sj[kM5] = zsj;

  //DiJetSelection 9 : 70<TpTJET<110, 30<ApTJET<40, kTy
	if (sj[0][0].Pt()>70 && sj[0][0].Pt()<110 &&
			sj[0][1].Pt()>30 && sj[0][1].Pt()<40 && 	
			Abs(TVector2::Phi_0_2pi(sj[0][0].DeltaPhi(sj[0][1]))-pi) < pi/3.){	
		sj[kM6] =sj[0]; fDijetSelectionCut[kM6] = true;}
  else sj[kM6] = zsj;

  //DiJetSelection 10 : 70<TpTJET<110, 40<ApTJET, kTy
	if (sj[0][0].Pt()>70 && sj[0][0].Pt()<110 && sj[0][1].Pt()>40  && 	
			Abs(TVector2::Phi_0_2pi(sj[0][0].DeltaPhi(sj[0][1]))-pi) < pi/3.){	
		sj[kM7] =sj[0]; fDijetSelectionCut[kM7] = true;}
  else sj[kM7] = zsj;


  //==== FILL
  fHistos->FillTH1("hJetMulti", sumJets.size() );

  //===============================
  // ptPair : leading - subleading
  //===============================


	for( int ids=kBDiJetSelBegin;ids<kBDiJetSelEnd;ids++ ){ // NOTE : begin with 1
		auto j = sj[ids];
		//=== SKIP Empty DiJet
		if( j[0].E() < 1e-4 || j[1].E() < 1e-4 ) {
			//fDijetSelectionCut[ids] = false;
			if (task){
				Double_t truthptpair =  (task->GetDijetPtPair()).at(ids);
				Double_t truthinvM =  (task->GetDijetInvM()).at(ids);
        Bool_t  truthdijetcut = (task->GetDijetSelectionCut()).at(ids);
				if (task->GetIsGenGoodVtx() && truthdijetcut ) {
					FillTHnSparse( "hDiJetInvMPtPairMiss", {(double)ids,fCent,j[0].Pt(),truthinvM,truthptpair},sf);
				}
			}
		}
    if (!fDijetSelectionCut[ids]) continue;
		Double_t nphi=0;
		auto diJetSel = Double_t( ids );
		auto dijet    = j[0] + j[1];
		auto invM     = dijet.M();
		auto dPhi     = j[0].DeltaPhi(j[1]);
		auto dPhiA    = Abs(dPhi);
		auto dPhi_0_2pi  = TVector2::Phi_0_2pi(dPhi);
		auto tpt      = j[0].Pt();
		auto apt      = j[1].Pt();
		auto ptpair   = dijet.Pt(); 
		auto testdphi = nphi-j[0].Phi();
		auto ptAsim   = (tpt-apt)/(tpt+apt);
		auto eAsim    = (j[0].E()-j[1].E())/(j[0].E()+j[1].E());
    auto kty      = j[0].Pt() * TMath::Sin(dPhi_0_2pi);
    if (ids>=kM1) ptpair = kty;

		double tratio = 1.;
		double tratioh = 1.;
		double tratiol = 1.;
		if (fName.Contains("EMCEJE")){
			tratio = tsf->Eval(j[0].Pt());
			tratioh = tsfh->Eval(j[0].Pt());
			tratiol = tsfl->Eval(j[0].Pt());
		}
		
		fDijetPtPair[ids] = ptpair;
		fDijetInvM[ids]=invM;
		if (task ){	
			Double_t truthptpair =  (task->GetDijetPtPair()).at(ids);
			Double_t truthinvM =  (task->GetDijetInvM()).at(ids);
			Bool_t truthdijetcut = (task->GetDijetSelectionCut().at(ids));
			if ( task->GetIsGenGoodVtx()) {
				if (IsGoodVertex){
					if (truthdijetcut){
						FillTHnSparse( "hDiJetInvMPtPairRes", { diJetSel,fCent,tpt,invM,truthinvM,ptpair,truthptpair},sf);
					} else {
						FillTHnSparse( "hDiJetInvMPtPairFake", {diJetSel,fCent,tpt,invM,ptpair},sf);
						FillTHnSparse( "hDiJetInvMPtPairMiss", {diJetSel,fCent,tpt,truthinvM,truthptpair},sf);
					}
				} else FillTHnSparse( "hDiJetInvMPtPairMiss", {diJetSel,fCent,tpt,truthinvM,truthptpair},sf);

			}
		}

		if (IsGoodVertex){
			FillTHnSparse( "hJetPtLeading",        { diJetSel, fCent, tpt},sf*tratio );
			FillTHnSparse( "hDiJetInvM",           { diJetSel, fCent, tpt, invM },sf*tratio);
			FillTHnSparse( "hDiJetDPhi_0_2pi",     { diJetSel, fCent, tpt, invM, dPhi_0_2pi },sf*tratio);
			FillTHnSparse( "hDiJetPtPair",         { diJetSel, fCent, tpt, invM, ptpair },sf*tratio); 
			FillTHnSparse( "hDiJetPtPairH",   { diJetSel, fCent, tpt, invM, ptpair },sf*tratioh); 
			FillTHnSparse( "hDiJetPtPairL",   { diJetSel, fCent, tpt, invM, ptpair },sf*tratiol); 
		}
	}

  PostData(1, fHistos->GetListOfHistograms());;
  //PostData(1, fOutput);
  return kTRUE;
}
void AliBSDiJetTask::FinishTaskOutput()
{
}
//___________________________________________________________________

THnSparse * AliBSDiJetTask::CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt){
  const TAxis * axises[bins.size()];
  for( UInt_t i=0;i<bins.size();i++ ) axises[i]= &bins[i];
  THnSparse * h= fHistos->CreateTHnSparse(name, title, ndim, axises,opt );
  return h;
}

THnSparse * AliBSDiJetTask::CreateTHnSparse(TString name, TString title, TString templ, Option_t * opt){
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

Long64_t AliBSDiJetTask::FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w ){
  auto hsparse = dynamic_cast<THnSparse*>( fHistos->FindObject(name) );
  if(! hsparse ){
    cout<<"ERROR : no "<<name<<endl;
    exit(1);
  }
  return FillTHnSparse( hsparse, x, w );
}

/*
   Long64_t AliBSDiJetTask::Fill( TString name, std::vector<Double_t> x, Double_t w ){
   return FillTHnSparse( name, x, w );
   }
   */

Long64_t AliBSDiJetTask::FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w ){
  if( int(x.size()) != h->GetNdimensions() ){
    cout<<"ERROR : wrong sized of array while Fill "<<h->GetName()<<endl;
    exit(1);
  }
  return h->Fill( &x.front(), w );
}

TAxis AliBSDiJetTask::AxisFix( TString name, int nbin, Double_t xmin, Double_t xmax ){ 
  TAxis axis(nbin, xmin, xmax);axis.SetName(name);
  return axis; 
}

TAxis AliBSDiJetTask::AxisStr( TString name, std::vector<TString> bin ){ 
  TAxis ax = AxisFix( name, bin.size(), 0.5, bin.size()+0.5);
  UInt_t i=1;
  for( auto blabel : bin )
    ax.SetBinLabel( i++, blabel );
  return ax; 
}

TAxis AliBSDiJetTask::AxisVar( TString name, std::vector<Double_t> bin ){
  TAxis axis( bin.size()-1, &bin.front() ) ;axis.SetName(name);
  return axis;
}

TAxis AliBSDiJetTask::AxisLog( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0){
  int binoffset = ( xmin0<0 || (xmin-xmin0)<1e-9) ? 0 : 1;
  std::vector<Double_t> bin(nbin+1+binoffset,0);
  double logBW3 = (log(xmax)-log(xmin))/nbin;
  for(int ij=0;ij<=nbin;ij++) bin[ij+binoffset]=xmin*exp(ij*logBW3);
  TAxis axis( nbin, &bin.front() ) ;
  axis.SetName(name);
  return axis;
}
