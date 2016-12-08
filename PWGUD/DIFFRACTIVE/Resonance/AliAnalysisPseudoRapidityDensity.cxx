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
// Simple class for dn/deta analyses. 
// by Beomkyu KIM
//==================================================================

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include <TRandom3.h>
#include "AliAnalysisPseudoRapidityDensity.h"
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
const Double_t        pi = TMath::Pi();
const Int_t kNstable=20;
const Int_t pdgStable[20] = {
	22,             // Photon
	11,             // Electron
	12,             // Electron Neutrino
	13,             // Muon
	14,             // Muon Neutrino
	15,             // Tau
	16,             // Tau Neutrino
	211,            // Pion
	321,            // Kaon
	311,            // K0
	130,            // K0s
	310,            // K0l
	2212,           // Proton
	2112,           // Neutron
	3122,           // Lambda_0
	3112,           // Sigma Minus
	3222,           // Sigma Plus
	3312,           // Xsi Minus
	3322,           // Xsi0
	3334            // Omega
};

AliAnalysisPseudoRapidityDensityRunTable::AliAnalysisPseudoRapidityDensityRunTable() :
    fCollisionType(kUnknownCollType)
{;}

AliAnalysisPseudoRapidityDensityRunTable::AliAnalysisPseudoRapidityDensityRunTable(Int_t runnumber) 
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
    else if (runnumber>=244340 && runnumber<=244628) fCollisionType=kPP;//LHC15n
    else fCollisionType=kUnknownCollType;
}
AliAnalysisPseudoRapidityDensityRunTable::~AliAnalysisPseudoRapidityDensityRunTable()
{;}

//___________________________________________________________________
AliAnalysisPseudoRapidityDensity::AliAnalysisPseudoRapidityDensity()
:AliAnalysisTaskSE("AliAnalysisPseudoRapidityDensity")
    , fOption() 
    , goodtrackindices()
{
    DefineOutput (1, TDirectory::Class());
}
//___________________________________________________________________
AliAnalysisPseudoRapidityDensity::AliAnalysisPseudoRapidityDensity
( 
      const char *name
    , const char *option
)
:AliAnalysisTaskSE(name)
    , fOption(option)
    , goodtrackindices()
{
    DefineOutput (1, TDirectory::Class());
}
//___________________________________________________________________
AliAnalysisPseudoRapidityDensity::AliAnalysisPseudoRapidityDensity
(
      const AliAnalysisPseudoRapidityDensity& ap
)
    : fOption(ap.fOption)
    , goodtrackindices(ap.goodtrackindices)
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
    delete fHistos; 
    delete fOutput;
    delete fTrigger;
    delete fTrackCuts;
    delete fPIDResponse;
    delete fRunTable;
    delete fRandom;
}

//___________________________________________________________________
void AliAnalysisPseudoRapidityDensity::UserCreateOutputObjects()

{
	fRandom = new TRandom3;
	fRandom->SetSeed();
	// Histograms container
	fOutput = new TList();
	fOutput->SetOwner(kTRUE);

	// Offline triggers -----------------------------------------------------
	fTrigger = new AliTriggerAnalysis; // offline trigger
	fTrigger -> SetFMDThreshold(0.3,0.5); // FMD threshold
	//-----------------------------------------------------------------------
  

	// TrackCuts for strangeness measure-------------------------------------
	fTrackCuts = new AliESDtrackCuts();
	{
		fTrackCuts -> GetStandardITSTPCTrackCuts2010(1,0);
	}
	fHistos = new THistManager("dndeta");

	//auto binType = AxisStr("Type",{"PN","PP","NN","Mixing"});
	binCent = AxisFix("Cent",10,0,100);
  binEta = AxisFix("Eta",60,-6,6);
  binEventClass   = AxisStr( "EventClass", { "DATA", "INEL","NSD"
		,"INELg0","SD","DD" ,"ND" } );
  binTriggClass  = AxisStr("TriggClass",{"MBOR","MBAND","MBORg0"});
  binParType = AxisStr("ParticleType",{"Pion","Kaon","Proton","Opar"});
  binV0Type = AxisStr("V0Type",{"K0s","Lambda","AntiLambda"});
 
	//auto binPt   = AxisFix("Pt",200,0,20);
	//auto binMass = AxisFix("Mass",200,0,2);
	//CreateTHnSparse("hInvMass","InvMass",4,{binType,binCent,binMass,binPt},"s");

	std::vector<TString> ent ={"All","PS","PSpileup","Goodz","Goodzcut"};
	auto h = fHistos->CreateTH1("hEventNumbers","",ent.size(), 0, ent.size());
	for(auto i=0u;i<ent.size();i++) h->GetXaxis()->SetBinLabel(i+1,ent.at(i).Data());
	fHistos -> CreateTH2("hPhiEta","",180,0,2*pi,20,-1,1);
	binZ = AxisFix("Z",100,-50,50);

  fHistos -> CreateTH1("hdndeta","",60,-6,6);
	fHistos -> CreateTH1("hdndetamc","",60,-6,6);
  fForward2d = fHistos ->CreateTH2 ("hfmddetadphi","",200, -4, 6, 20, 0, 2*pi);
  fHistos -> CreateTH1("hcent","cent",100,0,100);
  
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
 
	CreateTHnSparse( "hrecdndeta", "rec dndeta", 6, {
			binEventClass, binTriggClass, binCent, binZ, binParType, binEta  },"s");
	CreateTHnSparse( "hreczvtx", "rec zvtx", 4, { 
			binEventClass, binTriggClass,binCent,binZ},"s");
	CreateTHnSparse( "hkinedndeta", "kine only dndeta", 5, {
			binEventClass, binCent, binZ, binParType ,binEta  },"s");
	CreateTHnSparse( "hkinezvtx", "kine zvtx", 3, { binEventClass, binCent,binZ},"s");
  CreateTHnSparse( "hv0mass", "V0 inv mass", 4, 
		{ binTriggClass, binCent,binV0Type, AxisFix("v0mass",150,0,3)},"s");
  CreateTHnSparse( "hv0eta", "V0 daughter eta", 4, 
		{ binTriggClass, binCent,binV0Type, binEta},"s");
  
	CreateTHnSparse( "hrecdndetaptup", "rec dndeta pt up", 6, {
			binEventClass, binTriggClass, binCent, binZ, binParType, binEta  },"s");
	CreateTHnSparse( "hrecdndetaptdw", "rec dndeta pt dw", 6, {
			binEventClass, binTriggClass, binCent, binZ, binParType, binEta  },"s");
	CreateTHnSparse( "hkinedndetaptup", "kine only dndeta pt up", 5, {
			binEventClass, binCent, binZ, binParType ,binEta  },"s");
	CreateTHnSparse( "hkinedndetaptdw", "kine only dndeta pt dw", 5, {
			binEventClass, binCent, binZ, binParType ,binEta  },"s");


	Double_t xAxis1[26] = {0.01, 0.0154155, 0.0237639, 0.0366333, 0.0564723, 0.0870551, 0.1342, 0.206877, 0.318912, 0.49162, 0.757858, 1.16828, 1.80097, 2.77629, 4.2798, 6.59754, 10.1705, 15.6783, 24.169, 37.2578, 57.4349, 88.539, 136.488, 210.403, 324.348, 500};

	PostData(1, fHistos->GetListOfHistograms());

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
		fRunTable = new AliAnalysisPseudoRapidityDensityRunTable(runnumber);
	}



	// ----------------------------------------------------------------------
	// centrality
	fCent = -999;
	if (fRunTable->IsAA() || fRunTable->IsPA()){
		AliCentrality *cent = event->GetCentrality();
		if( ! cent ) return;
		fCent = cent->GetCentralityPercentile("V0M");
	} else {
		//Make sure naming convention is followed!
		AliMultSelection* sel = (AliMultSelection*) fEvt -> FindListObject("MultSelection");
		if (sel) {
			fCent = sel->GetMultiplicityPercentile("V0M");
		}
  }
  fHistos->FillTH1("hcent",fCent);


	
	// Pointer to a MC event-------------------------------------------------
	AliMCEvent *mcEvent = MCEvent();
	// ----------------------------------------------------------------------

	Bool_t IsMC = kFALSE;
	TArrayF vtxMC3d(3);

	Double_t eta, phi, pt;
	AliGenPythiaEventHeader* pythiaGenHeader = NULL;
	AliGenDPMjetEventHeader* dpmHeader = NULL;
  Bool_1d bevtc(kECend,false);
	sdweightingfactor = 1;
	if (mcEvent) {
    ismc = true;
		pythiaGenHeader =
			dynamic_cast<AliGenPythiaEventHeader*>(mcEvent->GenEventHeader());
		dpmHeader =
			dynamic_cast<AliGenDPMjetEventHeader*>(mcEvent->GenEventHeader());

		if(pythiaGenHeader) { // Pythia6
      //92 SD1, 93 SD2, 94 DD, -2 CD , -1 ND, 91 EL
			Int_t proc = pythiaGenHeader->ProcessType();
      if (proc == 92 || proc == 93) bevtc[kSD] = true;
      else if (proc == 94) bevtc[kDD] = true;
      //else if (proc == -1 || proc == -2) bevtc[kND] = true;

			if (fOption.Contains("PYTHIA8")){
        // 103 SD1, 104 SD2, 105 DD
				if (proc == 103 || proc == 104) bevtc[kSD] = true;
				else if (proc == 105) bevtc[kDD] = true;
				else bevtc[kND] = true;
			}
      if (!bevtc[kSD] && !bevtc[kDD] ) bevtc[kND] = kTRUE;
			if (!bevtc[kSD] )  bevtc[kNSD]  = kTRUE;
      if (bevtc[kSD] || bevtc[kDD] || bevtc[kND]) bevtc[kINEL] = true;
		}
 
		stack = mcEvent -> Stack();
		if (bevtc[kSD] && stack) this->MeasureDiffMass();
		Int_t nPrim  = stack->GetNprimary();
		mcEvent->GenEventHeader()->PrimaryVertex(vtxMC3d);
		zmcbin = binZ.FindBin(vtxMC3d[2])-1 ;
		for (Int_t i = 0; i < nPrim; i++){
			TParticle* part = stack->Particle(i);
			if (AliPWG0Helper::IsPrimaryCharged(part, nPrim) == kFALSE) continue;
			eta = part->Eta();
      if ( bevtc[kINEL] && fabs(eta) < 1.0 ) bevtc[kINELg0] = true;
     }
 
    Int_t pid ;
		for (Int_t i = 0; i < nPrim; i++){
			TParticle* part = stack->Particle(i);
			if (AliPWG0Helper::IsPrimaryCharged(part, nPrim) == kFALSE) continue;
			eta = part->Eta();
			phi = part ->Phi();
      pt  = part->Pt();
			fHistos->FillTH1("hdndetamc",eta);
      fHistos->FillTH1("hkinept",pt,1./pt);
			switch (TMath::Abs(part->GetPdgCode())){
				case 211:   pid = kPion;   break;
				case 321:   pid = kKaon;   break;
				case 2212:  pid = kProton; break;
				default:    pid = kOPar; break;
			} 
			for (auto ievtc=1u ; ievtc<kECend; ievtc++) {
				if (bevtc[ievtc]){ 
        	FillTHnSparse( "hkinedndeta", 
						{ Double_t(ievtc),fCent,vtxMC3d[2], Double_t(pid),eta}); 
          if (pt<0.05){
						FillTHnSparse( "hkinedndetaptup", 
								{ Double_t(ievtc),fCent,vtxMC3d[2], Double_t(pid),eta}, 2.);
						FillTHnSparse( "hkinedndetaptdw", 
								{ Double_t(ievtc),fCent,vtxMC3d[2], Double_t(pid),eta}, 0.5);
          } 
        }
      }
		}
		for (auto ievtc=1u ; ievtc<kECend; ievtc++) {
			if (bevtc[ievtc]) 
				FillTHnSparse("hkinezvtx",{Double_t(ievtc),fCent,vtxMC3d[2]});
		}
	} else bevtc[kDATA] = true;

	// Pointer to a MC event-------------------------------------------------
	//AliMCEvent *mcEvent = MCEvent();
	// ----------------------------------------------------------------------


	// Load InputHandler for each event---------------------------------------
	AliInputEventHandler* inputHandler = (AliInputEventHandler*)
		AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (IsFirstEvent){
		//fTree = inputHandler -> GetTree(); 
		IsFirstEvent = kFALSE;
	}
	// -----------------------------------------------------------------------

	fPIDResponse = (AliPIDResponse*) inputHandler->GetPIDResponse();
	if(!fPIDResponse){
		printf("AliAnalysisPseudoRapidityDensity No PIDd\n");
	}

	//PID Combined
	fPIDCombined = new AliPIDCombined;
	fPIDCombined->SetDefaultTPCPriors();//Need more update..
	fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);
	fPIDCombined->SetDetectorMask(
			AliPIDResponse::kDetTPC |
			AliPIDResponse::kDetTOF |
			AliPIDResponse::kDetITS |
			AliPIDResponse::kDetTRD);//Do we need??


	Bool_t IsMinimumBias = kFALSE;
	fHistos -> FillTH1("hEventNumbers","All",1);
  Bool_1d btrigc(kTrigend,false);

	IsMinimumBias = inputHandler -> IsEventSelected() & AliVEvent::kMB ;
  
  btrigc[kMBOR]  = IsMinimumBias
		&& (
				fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kV0A) ||
				fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kV0C) ||
				fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kADA) ||
				fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kADC) ||
				fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kSPDGFO))
		;
	btrigc[kMBAND] = IsMinimumBias
		&& (    fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kADA) ||
				fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kV0A) )
		&& (    fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kADC) ||
				fTrigger->IsOfflineTriggerFired(fEvt,AliTriggerAnalysis::kV0C))
		;


	if (IsMinimumBias)  fHistos -> FillTH1("hEventNumbers","PS",1);
	// Reject pile-up events--------------------------------------------------
	//if (!IsMC && event->IsPileupFromSPD(3.,0.8,3.,2.,5.)) return;
	if ((!IsAA || !mcEvent) && event->IsPileupFromSPD(3.,0.8,3.,2.,5.)) return;
	if (IsMinimumBias) fHistos -> FillTH1("hEventNumbers","PSpileup",1);
	// -----------------------------------------------------------------------

	const AliVVertex* trackVtx  = fEvt->GetPrimaryVertexTPC() ;
	const AliVVertex* spdVtx	  = fEvt->GetPrimaryVertexSPD() ;

	Bool_t IsGoodVertex = kFALSE;
	Bool_t IsGoodVertexCut = kFALSE;

	if (spdVtx) {
		fZ = spdVtx->GetZ();
		if (spdVtx->GetNContributors()<1) IsGoodVertex = kFALSE;
		else {	
			fHistos -> FillTH1("hEventNumbers","Goodz",1);
			IsGoodVertex = kTRUE;
		}
	} else IsGoodVertex = kFALSE;


	if ( IsGoodVertex && fabs(fZ)<30.) {
		IsGoodVertexCut = kTRUE;
		if (IsMinimumBias) {
			fHistos -> FillTH1("hEventNumbers","Goodzcut",1);
		}
	}

	if (IsMinimumBias && IsGoodVertexCut){
		fMultiplicity = fEvt -> GetMultiplicity();
    Double_t eta;
		for (auto it=0u; it<fMultiplicity->GetNumberOfTracklets(); it++) {
			Double_t eta = fMultiplicity->GetEta(it);
      if ( btrigc[kMBOR] && fabs(eta) < 1.0 ) btrigc[ kMBORg0 ] = true;
		} 
		this -> FillTracklets(bevtc,btrigc) ;
    this -> StrangenessMeasure(btrigc);
	}
	PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisPseudoRapidityDensity::FillTracklets(Bool_1d bevtc, Bool_1d btrigc){
	UInt_t nmul = fMultiplicity->GetNumberOfTracklets();
	Double_t eta, phi, pt;
	for (auto it=0u; it<nmul; it++) {
		eta = fMultiplicity->GetEta(it);
		phi = fMultiplicity->GetPhi(it);
		fHistos->FillTH2("hPhiEta",phi,eta);
		fHistos->FillTH1("hdndeta",eta);
    Int_t pid;
    if (ismc){
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

			//if (IsMotherStrangeParticle(mother)) IsMotherStrangeness = kTRUE;
			//else if (!stack->IsPhysicalPrimary(fMultiplicity->GetLabel(i,0))
			//		&& !stack->IsPhysicalPrimary(fMultiplicity->GetLabel(i,1)))
			//	IsBackground = kTRUE;

			if (fMultiplicity->GetLabel(it,0) == fMultiplicity->GetLabel(it,1)){
				eta = particle -> Eta();
				phi  = particle -> Phi();
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
    } else pid = kParDATA;
		for (auto ievtc=1u ; ievtc<kECend; ievtc++) {
			for (auto itrigc=1u ; itrigc<kTrigend; itrigc++) {
				if (bevtc[ievtc] && btrigc[itrigc] ){
					FillTHnSparse( "hrecdndeta", 
						{ Double_t(ievtc),Double_t(itrigc),fCent,fZ,Double_t(pid),eta}); 
          if (ismc && pt<0.05){
						FillTHnSparse( "hrecdndetaptup", 
								{ Double_t(ievtc),Double_t(itrigc),fCent,fZ,Double_t(pid),eta},2.); 
						FillTHnSparse( "hrecdndetaptdw", 
								{ Double_t(ievtc),Double_t(itrigc),fCent,fZ,Double_t(pid),eta},0.5); 
          }
        }
      }
		}
	}
	for (auto ievtc=1u ; ievtc<kECend; ievtc++) 
		for (auto itrigc=1u ; itrigc<kTrigend; itrigc++) 
			if (bevtc[ievtc] && btrigc[itrigc] )
				FillTHnSparse("hreczvtx",{Double_t(ievtc),Double_t(itrigc), fCent,fZ}); 
}



//___________________________________________________________________
void AliAnalysisPseudoRapidityDensity::FinishTaskOutput()
{
		OpenFile(1);
		TH1D *fForward1d = (TH1D*)fForward2d->ProjectionX("_x",1,-1,"e");
		TH1D* norm   = fForward2d->ProjectionX("norm", 0, 1, "");
		fForward1d -> Divide(norm);
		fForward1d -> Write();
}
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

void AliAnalysisPseudoRapidityDensity::MeasureDiffMass(){
	
	Int_t np = stack->GetNprimary();

	Int_t iPart1=-1;
	Int_t iPart2=-1;
  Double_t cms = 5020;

	Double_t y1 = 1e10;
	Double_t y2 = -1e10;
	for (Int_t i = 0; i < np; ++i){
		TParticle* part = stack->Particle(i);

		Int_t statusCode = part->GetStatusCode();

		// Initial state particle
		if (statusCode != 1)
			continue;

		Int_t pdg = TMath::Abs(part->GetPdgCode());
		Bool_t isStable = kFALSE;
		for (Int_t i1 = 0; i1 < kNstable; i1++) {
			if (pdg == pdgStable[i1]) {
				isStable = kTRUE;
				break;
			}
		}
		if(!isStable)
			continue;

		Double_t y = part->Y();

		if (y < y1)
		{
			y1 = y;
			iPart1 = i;
		}
		if (y > y2)
		{
			y2 = y;
			iPart2 = i;
		}
		if(iPart1>=0 && iPart2 >=0) {
			y1=TMath::Abs(y1);
			y2=TMath::Abs(y2);

			TParticle *  part1 = (TParticle *) stack->Particle(iPart1);
			TParticle *  part2 = (TParticle *) stack->Particle(iPart2);

			Int_t pdg1 = part1->GetPdgCode();
			Int_t pdg2 = part2->GetPdgCode();


			Int_t iPart = -1;
			if (pdg1 == 2212 && pdg2 == 2212)
			{
				if(y1 > y2)
					iPart = iPart1;
				else if(y1 < y2)
					iPart = iPart2;
				else {
					iPart = iPart1;
					if(fRandom->Uniform(0.,1.)>0.5) iPart = iPart2;
				}
			}
			else if (pdg1 == 2212)
				iPart = iPart1;
			else if (pdg2 == 2212)
				iPart = iPart2;

			Double_t M=-1.;
			if(iPart>0) {
				TParticle *  part = (TParticle *) stack->Particle(iPart);
				Double_t E= part->Energy();
				Double_t P= part->P();
				//energy 13 TeV = 13000
				Double_t M2 = (cms-E-P)*(cms-E+P);
				if(M2>0) {
					M= TMath::Sqrt(M2);
					fHistos->FillTH1("hsdmass",M);
				}
			}
		} // end of iPart1>=0 && iPart2 >=0
	}
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


		v0->ChangeMassHypothesis(310);
    for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
    	if (btrigc[itrigc])
    		FillTHnSparse("hv0mass",
					{Double_t(itrigc),fCent,Double_t(kK0s),v0->GetEffMass()});
		if (v0->GetEffMass() >0.482 && v0->GetEffMass()<0.509) {
			for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
				if (btrigc[itrigc]) {
					FillTHnSparse("hv0eta",
							{Double_t(itrigc),fCent,Double_t(kK0s),pTrack->Eta()});
					FillTHnSparse("hv0eta",
							{Double_t(itrigc),fCent,Double_t(kK0s),nTrack->Eta()});
				}
		}
		v0->ChangeMassHypothesis(3122);
    for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
    	if (btrigc[itrigc])
    		FillTHnSparse("hv0mass",
					{Double_t(itrigc),fCent,Double_t(kLambda),v0->GetEffMass()});
		if (v0->GetEffMass() >1.11 && v0->GetEffMass()<1.12) {
			for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
				if (btrigc[itrigc]) {
					FillTHnSparse("hv0eta",
							{Double_t(itrigc),fCent,Double_t(kLambda),pTrack->Eta()});
					FillTHnSparse("hv0eta",
							{Double_t(itrigc),fCent,Double_t(kLambda),nTrack->Eta()});
				}
		}

		v0->ChangeMassHypothesis(-3122);
    for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
    	if (btrigc[itrigc])
    		FillTHnSparse("hv0mass",
					{Double_t(itrigc),fCent,Double_t(kAntiLambda),v0->GetEffMass()});
		if (v0->GetEffMass() >1.11 && v0->GetEffMass()<1.12) { 
			for (auto itrigc=1u ; itrigc<kTrigend; itrigc++)
				if (btrigc[itrigc]) {
					FillTHnSparse("hv0eta",
							{Double_t(itrigc),fCent,Double_t(kAntiLambda),pTrack->Eta()});
					FillTHnSparse("hv0eta",
							{Double_t(itrigc),fCent,Double_t(kAntiLambda),nTrack->Eta()});
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


