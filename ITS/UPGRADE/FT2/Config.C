//
// Configuration for ITS Upgrade TDR simulations
//
// 1 PbPb HIJING event 5.5 TeV with b<5 fm (0-10%)
// +
// 15 per event per type of
//   B+->D0pi
//
// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++")
//
//
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "PYTHIA6/AliGenPythia.h"
#include "TDPMjet/AliGenDPMjet.h"
#include "STEER/AliMagFCheb.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv3.h"
#include "STRUCT/AliDIPOv3.h"
#include "STRUCT/AliHALLv3.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv3.h"
#include "STRUCT/AliPIPEv3.h"
#include "ITS/AliITSv11.h"
#include "ITS/UPGRADE/AliITSUv0.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv6T0.h"
#include "HMPID/AliHMPIDv3.h"
#include "ZDC/AliZDCv3.h"
#include "TRD/AliTRDv1.h"
#include "TRD/AliTRDgeometry.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PHOS/AliPHOSSimParam.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "ACORDE/AliACORDEv1.h"
#include "VZERO/AliVZEROv7.h"
//
//#include "ITS/UPGRADE/AliITSUv11.h"

#endif


enum PDC06Proc_t
{
	kPythia6, kPythia6D6T, kPythia6ATLAS, kPythia6ATLAS_Flat, kPythiaPerugia0, kPhojet, kPythiaPerugia0chadr, kPythiaPerugia0bchadr, kPythiaPerugia0cele, kPythiaPerugia0bele, kPythiaPerugia0Jpsi2e, kPythiaPerugia0BtoJpsi2e, kHijing, kHijing2500, kHijing2500HF, kHydjet, kDpmjet, kAmptHF, kAmpt, kRunMax
};

const char * pprRunName[] = {
	"kPythia6", "kPythia6D6T", "kPythia6ATLAS", "kPythia6ATLAS_Flat", "kPythiaPerugia0", "kPhojet",  "kPythiaPerugia0chadr", "kPythiaPerugia0bchadr", "kPythiaPerugia0cele", "kPythiaPerugia0bele", "kPythiaPerugia0Jpsi2e", "kPythiaPerugia0BtoJpsi2e", "kHijing", "kHijing2500", "kHijing2500HF", "kHydjet", "kDpmjet", "kAmptHF", "kAmpt"
};

enum Mag_t
{
	kNoField, k5kG, kFieldMax
};

const char * pprField[] = {
	"kNoField", "k5kG"
};

enum PprTrigConf_t
{
	kDefaultPPTrig, kDefaultPbPbTrig
};

const char * pprTrigConfName[] = {
	"p-p","Pb-Pb"
};


//--- Functions ---
class AliGenPythia;
//
AliGenerator *Hijing();
AliGenerator *Hijing2500();
AliGenerator *Hijing2500HF(Int_t typeHF);
AliDecayer *decayer = 0;

// Geterator, field, beam energy
static PDC06Proc_t   proc     = kHijing2500HF;
static Mag_t         mag      = k5kG;
static Float_t       energy   = 5500.;				// energy in CMS
static Float_t       bMin     = 0.;
static Float_t       bMax     = 0.;					// 0-5 fm corresponds to around 0-10% (see https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies#Tables_with_centrality_bins_for)
static PprTrigConf_t strig = kDefaultPbPbTrig;		// default pp trigger configuration
static Double_t      JpsiPol  = 0;					// Jpsi polarisation
static Bool_t        JpsiHarderPt = kFALSE;			// Jpsi harder pt spectrum (8.8 TeV)

void ProcessEnvironmentVars();

//========================//
// Set Random Number seed //
//========================//
TDatime dt;
static UInt_t seed    = dt.Get();

// Comment line
static TString comment;

void Config()
{
	
	printf("Config.C: Initializing Macro\n");
	// Get settings from environment variables
	ProcessEnvironmentVars();
	
	gRandom->SetSeed(seed);
	cerr<<"Seed for random number generation= "<<seed<<endl;
	gSystem->Load("libITSUpgradeBase.so");
	gSystem->Load("libITSUpgradeSim.so");
	gSystem->Load("libEVGEN");
	// Libraries required by geant321
#if defined(__CINT__)
	gSystem->Load("liblhapdf");      // Parton density functions
	gSystem->Load("libEGPythia6");   // TGenerator interface
	if (proc == kPythia6 || proc == kPhojet || proc == kDpmjet) {
		gSystem->Load("libpythia6");        // Pythia 6.2
		gSystem->Load("libAliPythia6");     // ALICE specific implementations
	} else if (proc != kHydjet) {
		gSystem->Load("libpythia6.4.21");   // Pythia 6.4
		gSystem->Load("libAliPythia6");     // ALICE specific implementations
	}
	
	gSystem->Load("libhijing");
	gSystem->Load("libTHijing");
	
	gSystem->Load("libgeant321");

#endif
	
	TGeant3TGeo *g3 = new TGeant3TGeo("C++ Interface to Geant3");
	g3->SetSWIT(4,10000);

	
	//=======================================================================
	//  Create the output file
	
	
	AliRunLoader* rl=0x0;
	
	cout<<"Config.C: Creating Run Loader ..."<<endl;
	rl = AliRunLoader::Open("galice.root",
													AliConfig::GetDefaultEventFolderName(),
													"recreate");
	if (rl == 0x0)
	{
		gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
		return;
	}
	rl->SetCompressionLevel(2);
	rl->SetNumberOfEventsPerFile(1000);
	gAlice->SetRunLoader(rl);
	
	// Set the trigger configuration
	AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
	cout<<"Trigger configuration is set to  "<<pprTrigConfName[strig]<<endl;
	
	//
	//=======================================================================
	// ************* STEERING parameters FOR ALICE SIMULATION **************
	// --- Specify event type to be tracked through the ALICE setup
	// --- All positions are in cm, angles in degrees, and P and E in GeV
	
	
	gMC->SetProcess("DCAY",1);
	gMC->SetProcess("PAIR",0); // 1 pair production
	gMC->SetProcess("COMP",0); // 1 compton scattering
	gMC->SetProcess("PHOT",0); // 1 photoelectric effect
	gMC->SetProcess("PFIS",0); // 0 nuclear fission by photon
	gMC->SetProcess("DRAY",0); // 0 delta ray production
	gMC->SetProcess("ANNI",0); // 1 positron annihilation
	gMC->SetProcess("BREM",0); // 1 bremsstrahlung
	gMC->SetProcess("MUNU",0); // 1 muon nucleus interaction
	gMC->SetProcess("CKOV",0); // 1 cherenkov
	gMC->SetProcess("HADR",0); // 1 hadronic interactions
	gMC->SetProcess("LOSS",0); // 2 energy loss
	gMC->SetProcess("MULS",0); // 1 mutliple scattering
	gMC->SetProcess("RAYL",0); // 1 Rayleigh effect
	
	Float_t cut = 1.e-3;        // 1MeV cut by default
	Float_t tofmax = 1.e10;
	
	gMC->SetCut("CUTGAM", cut);
	gMC->SetCut("CUTELE", cut);
	gMC->SetCut("CUTNEU", cut);
	gMC->SetCut("CUTHAD", cut);
	gMC->SetCut("CUTMUO", cut);
	gMC->SetCut("BCUTE",  cut);
	gMC->SetCut("BCUTM",  cut);
	gMC->SetCut("DCUTE",  cut);
	gMC->SetCut("DCUTM",  cut);
	gMC->SetCut("PPCUTM", cut);
	gMC->SetCut("TOFMAX", tofmax);
	
	// RANDOM SELECTION OF ONE OF THE SEVEN GENERATION TYPES
	//
	Int_t typeHF  = 0.;
	
	//======================//
	// Set External decayer //
	//======================//
	TVirtualMCDecayer* decayer = new AliDecayerPythia();
	decayer->SetForceDecay(kHadronicDWithout4Bodies);
	decayer->Init();
	gMC->SetExternalDecayer(decayer);
	
	//=========================//
	// Generator Configuration //
	//=========================//
	AliGenerator* gener = 0x0;
	gener = Hijing2500HF(typeHF);
	
	// abort GEANT propagation outside of some range
	// x,z in cm
	gAlice->GetMCApp()->TrackingLimits(50,250);
	

	
	//
	//
	// Size of the interaction diamond
	// Longitudinal
	Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
	
	//
	// Transverse
	Float_t betast  = 3.5;                      // beta* [m]
	Float_t eps     = 3.75e-6;                   // emittance [m]
	Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
	Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
	
	printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
	
	gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
	gener->SetVertexSmear(kPerEvent);
	gener->Init();
	
	printf("\n \n Comment: %s \n \n", comment.Data());
	
	//
	// FIELD
	//
	TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG,
																											 AliMagF::kBeamTypeAA, 1380.));
	
	
	rl->CdGAFile();
	
	//=================== Alice BODY parameters =============================
	AliBODY *BODY = new AliBODY("BODY", "Alice envelop");
	AliPIPE *PIPE = new AliPIPEupgrade("PIPE", "Beam Pipe",1.8,0.08);
	gROOT->ProcessLine(".x $ALICE_ROOT/ITS/UPGRADE/testITSU/CreateITSU.C");
	
}
void ProcessEnvironmentVars()
{

	cout << "#######################################" << endl;
	cout << "# Simulation Setup" << endl;
	cout << "#######################################" << endl;
	cout << "# Energy set to "<< energy <<" GeV"<<endl;
	cout << "# Seed set to "<< seed <<endl;
	cout << "# Impact parameter in ["<< bMin << "," << bMax << "]" <<endl;
	cout << "#######################################" << endl;

}

AliGenerator* Hijing()
{
	AliGenHijing *gener = new AliGenHijing(-1);
	// centre of mass energy
	gener->SetEnergyCMS(energy);
	gener->SetImpactParameterRange(bMin, bMax);
	// reference frame
	gener->SetReferenceFrame("CMS");
	// projectile
	gener->SetProjectile("A", 208, 82);
	gener->SetTarget    ("A", 208, 82);
	// tell hijing to keep the full parent child chain
	gener->KeepFullEvent();
	// enable jet quenching
	gener->SetJetQuenching(1);
	// enable shadowing
	gener->SetShadowing(1);
	// Don't track spectators
	gener->SetSpectators(0);
	// kinematic selection
	gener->SetSelectAll(0);
	return gener;
}

AliGenerator* Hijing2500()
{
	AliGenHijing *gener = (AliGenHijing*) Hijing();
	gener->SetJetQuenching(0);
	gener->SetPtHardMin (3.7);
	return gener;
}

AliGenerator* Hijing2500HF(Int_t typeHF)
{
	comment = comment.Append("Modified by Johannes: PbPb: Hjing2500 at 5.5 + ITS Upgrade signals");
	AliPDG::AddParticlesToPdgDataBase();
	
	AliGenCocktail *cocktail = new AliGenCocktail();
	
	cocktail->SetProjectile("A", 208, 82);
	cocktail->SetTarget    ("A", 208, 82);
	cocktail->SetEnergyCMS(energy);
	//
	// 1 Hijing event
	TFormula* one    = new TFormula("one",    "1.");
	// provides underlying event and collision geometry
	AliGenHijing *hijing = Hijing2500();
	cocktail->AddGenerator(hijing,"hijing",1);
	Float_t thminH = (180./TMath::Pi())*2.*atan(exp(-2.5));
	Float_t thmaxH = (180./TMath::Pi())*2.*atan(exp( 2.5));
	hijing->SetChildThetaRange(thminH,thmaxH);
	//	hijing->SetDecaysOff(0);  // HIJIJNG decays all particles, including weak strange decayss
	
	Float_t thmin          = (180./TMath::Pi())*2.*atan(exp(-1.));
	Float_t thmax          = (180./TMath::Pi())*2.*atan(exp( 1.));
	
	
	AliGenParam *gen[2];
	UInt_t partId[1] = {AliGenITSULib::kBplus};
	for(Int_t iPart=0; iPart<2 ; iPart++){
		if(iPart%2==0) gen[iPart] = new AliGenParam(15,new AliGenITSULib(),partId[iPart/2],"DIST");
		if(iPart%2==1) gen[iPart]= new AliGenParam(15,new AliGenITSULib(),-partId[iPart/2],"DIST");
	 gen[iPart]->SetDecayer(decayer);
	 gen[iPart]->SetPtRange(0.,999.);
		gen[iPart]->SetPhiRange(0., 360.);
		gen[iPart]->SetYRange(-1.,1.);
		gen[iPart]->SetCutOnChild(1);
		gen[iPart]->SetChildThetaRange(thmin,thmax);
		gen[iPart]->SetSelectAll(kTRUE);
		gen[iPart]->SetForceDecay(kBeautyUpgrade);
		//	gen[iPart]->SetPreserveFullDecayChain(kTRUE);
		cocktail->AddGenerator(gen[iPart], Form("Generator_%i_%i",partId[iPart/2],iPart%2), 1);
	}

	return cocktail;
}
//---------------------------------------
AliGenerator* MbPythiaTunePerugia0chadr()
{
	comment = comment.Append(" pp: Pythia (Perugia0) chadr (1 ccbar per event, 1 c-hadron in |y|<1.5, chadrons decay to hadrons");
	//
	//    Pythia
	AliGenPythia* pythia = new AliGenPythia(-1);
	pythia->SetMomentumRange(0, 999999.);
	pythia->SetThetaRange(0., 180.);
	pythia->SetYRange(-1.,1.);
	pythia->SetPtRange(0,1000.);
	pythia->SetProcess(kPyCharmppMNRwmi);
	pythia->SetEnergyCMS(energy);
	//    Tune
	//    320     Perugia 0
	pythia->SetTune(320);
	pythia->UseNewMultipleInteractionsScenario();
	//
	//    decays
	pythia->SetForceDecay(kHadronicDWithout4Bodies);
	
	//    write only HF sub event
	pythia->SetStackFillOpt(AliGenPythia::kHeavyFlavor);
	
	return pythia;
}
