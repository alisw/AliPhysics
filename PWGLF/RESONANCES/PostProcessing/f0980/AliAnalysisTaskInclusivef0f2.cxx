/*

auther : JunLee Kim

*/
#include "Riostream.h"
#include <iomanip>

#include "TFile.h"
#include "TSystem.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskInclusivef0f2.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliVMultiplicity.h"
#include "AliMCEventHandler.h"
#include "AliVEventHandler.h"
#include "AliEventCuts.h"


using namespace std;


AliAnalysisTaskInclusivef0f2RunTable::AliAnalysisTaskInclusivef0f2RunTable() :
    fCollisionType(kUnknownCollType)
{;}

AliAnalysisTaskInclusivef0f2RunTable::AliAnalysisTaskInclusivef0f2RunTable(Int_t runnumber)
{   
    if (runnumber>=114737 && runnumber<=130850) fCollisionType = kPP; //LHC10bcde
    else if (runnumber>=144871 && runnumber<=146860) fCollisionType=kPP;//LHC11a
    else if (runnumber>=136851 && runnumber<=139517) fCollisionType=kAA;//LHC10h
    else if (runnumber>=167813 && runnumber<=170595) fCollisionType=kAA;//LHC11h
    else if (runnumber>=244917 && runnumber<=246392) fCollisionType=kAA;//LHC15o
    else if (runnumber>=188356 && runnumber<=188503) fCollisionType=kPA;//LHC12g
    else if (runnumber>=189122 && runnumber<=192732) fCollisionType=kPA;//LHC12h
    else if (runnumber>=195344 && runnumber<=195483) fCollisionType=kPA;//LHC13b
    else if (runnumber>=195529 && runnumber<=195677) fCollisionType=kPA;//LHC13c
    else if (runnumber>=195724 && runnumber<=195872) fCollisionType=kPA;//LHC13d
    else if (runnumber>=195955 && runnumber<=195872) fCollisionType=kPA;//LHC13e
    else if (runnumber>=197669 && runnumber<=200000) fCollisionType=kPA;//LHC13g
    else if (runnumber>=244340 && runnumber<=244628) fCollisionType=kPP;//LHC15n
    else if (runnumber>=235245 && runnumber<=236557) fCollisionType=kPP;//LHC15i
    else if (runnumber>=256504 && runnumber<=260014) fCollisionType=kPP;//LHC16kl
    else if (runnumber>=271870 && runnumber<=280140) fCollisionType=kPP;//LHC17hijklm
    else if (runnumber>=280282 && runnumber<=281961) fCollisionType=kPP;//LHC17o
    else if (runnumber>=282008 && runnumber<=282343) fCollisionType=kPP;//LHC17p
    else if (runnumber>=285471 && runnumber<=289971) fCollisionType=kPP;//LHC18
    else if (runnumber>=265309 && runnumber<=265525) fCollisionType=kPA;//LHC16q
    else if (runnumber>=295585 && runnumber<=296623) fCollisionType=kAA;//LHC18q
    else if (runnumber>=296690 && runnumber<=297595) fCollisionType=kAA;//LHC18r
//    else fCollisionType=kUnknownCollType;
    else fCollisionType=kPP;
}   
AliAnalysisTaskInclusivef0f2RunTable::~AliAnalysisTaskInclusivef0f2RunTable()
{;}


AliAnalysisTaskInclusivef0f2::AliAnalysisTaskInclusivef0f2() : AliAnalysisTaskSE("AliAnalysisTaskInclusivef0f2"),
	fEvt(), fOption(),
	fEMpool(),
	fEMpooltrk()
{
}


AliAnalysisTaskInclusivef0f2::AliAnalysisTaskInclusivef0f2(const char* name, const char *option) : AliAnalysisTaskSE(name),
	fEvt(), fOption(option),
	fEMpool(),
	fEMpooltrk()
{
 DefineOutput(1, TList::Class());
}


AliAnalysisTaskInclusivef0f2::AliAnalysisTaskInclusivef0f2(const AliAnalysisTaskInclusivef0f2& ap) :
	fOption(ap.fOption), goodtrackindices(ap.goodtrackindices),
	fEMpool(ap.fEMpool),
	fEMpooltrk(ap.fEMpooltrk)
{
 DefineOutput(1, TList::Class());
}


AliAnalysisTaskInclusivef0f2& AliAnalysisTaskInclusivef0f2::operator =(const AliAnalysisTaskInclusivef0f2& ap)
{
 DefineOutput(1, TList::Class());
 this->~AliAnalysisTaskInclusivef0f2();
 new(this) AliAnalysisTaskInclusivef0f2(ap);
 return *this;
}


AliAnalysisTaskInclusivef0f2::~AliAnalysisTaskInclusivef0f2()
{
// if(fOutputList) {
//	delete fOutputList; 
// }
 delete fTrigger;
 delete fPIDResponse;
// delete fPIDCombined;
 delete fTrackCuts;
 delete fRunTable;
}


void AliAnalysisTaskInclusivef0f2::UserCreateOutputObjects()
{
 fOutput = new TList();    
 fOutput->SetOwner();

/*
 fTrackCuts = new AliESDtrackCuts();
 fTrackCuts -> GetStandardITSTPCTrackCuts2010(1,0);
 fTrackCuts -> SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.01");
 fTrackCuts -> SetMaxDCAToVertexZ(0.5);
 fTrackCuts -> SetEtaRange(-1*fetacut,fetacut);
 fTrackCuts -> SetPtRange(0.15, 1e10);
*/


 Double1D varcentbinHeavy = {0,0.001,0.0033,0.01,0.02,0.033,0.05,0.1,0.2,0.5,1,5,10,15,20,30,40,50,60,70,80,100};
 binCent = AxisVar("Cent",varcentbinHeavy);

 Double1D verzbin = {-15,-10,-7,7,10,15};
 binZ = AxisVar("Z",verzbin);
 
 Double1D verptbin = {0.0, 0.3, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 13.0};
 binPt = AxisVar("Pt",verptbin);
 binPtGen = AxisFix("Pt",200,0.0,20.0);

 binType = AxisStr("Type",{"PN","PP","NN"});

 binMass = AxisFix("Mass",1000,0,5);
 if( fOption.Contains("Fine") ) binMass = AxisFix("Mass",100000,0,5);

 binCharge = AxisFix("Charge",3,-1.5,1.5);

 binTrackCutBit = AxisFix("TrackCut",7,0.5,7.5);

 binPID = AxisFix("PIDBit",4,0.5,4.5);

 binExKaonNum = AxisFix("KNum",5,0.5,5.5);

 binTrackPt = AxisFix("TrackPt",200,0,10);

 binSigma = AxisFix("Sigma",100,-10,10);

 binSwitch = AxisFix("Switch",2,-0.5,1.5);

 binEta = AxisFix("eta",32,-0.8,0.8);

 auto binV0Amp = AxisFix("binV0Amp",3000,0,3e3);
 auto binTrig = AxisFix("Trig",2,-0.5,1.5);
 auto binParType = AxisFix("ParType",2,-0.5,1.5);

 fHistos = new THistManager("Inclusivef0f2hists");

//Event Selection ****************
 vector<TString> ent ={
	"All","IsTriggered","IsNotPileup",
	"IsValidVtx","IsGoodVtx","IsSelectedFromAliMultSelection",
	"IsMultiplicityInsideBin" };

 auto h = fHistos->CreateTH1("hEventNumbers","",ent.size(), 0, ent.size());
 for(auto i=0u;i<ent.size();i++) h->GetXaxis()->SetBinLabel(i+1,ent.at(i).Data());

 fHistos->CreateTH1("hZvtx","",620,-15.5,15.5,"s");

 fHistos->CreateTH1("hHMT","",1000,0,1,"s");
 fHistos->CreateTH1("hMB","",100,0,100,"s");

 fHistos->CreateTH2("hMB_V0M","",100,0,100,1000,0,3000,"s");
 fHistos->CreateTH2("hHMT_V0M","",1000,0,1,1000,0,3000,"s");

 binCentForMC = AxisFix("CentMC",100,0,100);
//*****************************

//Distributions for correction in the event selection ****************

 if( fOption.Contains("MC") ){
 CreateTHnSparse("TrigEffMult","TrigEffMult",2,
	{binCentForMC,binTrig},"s");
 CreateTHnSparse("TrigEffMult0","TrigEffMult0",2,
        {binCentForMC,binTrig},"s");

 CreateTHnSparse("hRhoGenParticle","hRhoGenParticle",4,
        {binZ,binCentForMC,binPtGen,binMass},"s");
 CreateTHnSparse("hF0GenParticle","hF0GenParticle",4,
        {binZ,binCentForMC,binPtGen,binMass},"s");
 CreateTHnSparse("hF2GenParticle","hF2GenParticle",4,
        {binZ,binCentForMC,binPtGen,binMass},"s");

 CreateTHnSparse("VtxSelection","VtxSelection",2,
	{binCentForMC,binSwitch},"s");
 CreateTHnSparse("VtxSelection0","VtxSelection0",2,
        {binCentForMC,binSwitch},"s");

 CreateTHnSparse("SignalLoss","SignalLoss",3,
	{binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0","SignalLoss0",3,
        {binCentForMC,binPt,binSwitch},"s"); 

 CreateTHnSparse("SignalLossPion","SignalLossPion",4,
        {binCentForMC,binPt,binSwitch,binCharge},"s");
 CreateTHnSparse("SignalLoss0Pion","SignalLoss0Pion",4,
        {binCentForMC,binPt,binSwitch,binCharge},"s");

 CreateTHnSparse("SignalLossKaon","SignalLossKaon",3,
        {binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0Kaon","SignalLoss0Kaon",3,
        {binCentForMC,binPt,binSwitch},"s");

 CreateTHnSparse("SignalLossKaonpipi","SignalLossKaonpipi",3,
        {binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0Kaonpipi","SignalLoss0Kaonpipi",3,
        {binCentForMC,binPt,binSwitch},"s");

 CreateTHnSparse("SignalLossPhi","SignalLossPhi",3,
        {binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0Phi","SignalLoss0Phi",3,
        {binCentForMC,binPt,binSwitch},"s");

 CreateTHnSparse("SignalLossPionMt","SignalLossPionMt",4,
        {binCentForMC,binPt,binSwitch,binCharge},"s");
 CreateTHnSparse("SignalLoss0PionMt","SignalLoss0PionMt",4,
        {binCentForMC,binPt,binSwitch,binCharge},"s");

 CreateTHnSparse("SignalLossKaonMt","SignalLossKaonMt",3,
        {binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0KaonMt","SignalLoss0KaonMt",3,
        {binCentForMC,binPt,binSwitch},"s");

 CreateTHnSparse("SignalLossKaonpipiMt","SignalLossKaonpipiMt",3,
        {binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0KaonpipiMt","SignalLoss0KaonpipiMt",3,
        {binCentForMC,binPt,binSwitch},"s");

 CreateTHnSparse("SignalLossPhiMt","SignalLossPhiMt",3,
        {binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0PhiMt","SignalLoss0PhiMt",3,
        {binCentForMC,binPt,binSwitch},"s");



 CreateTHnSparse("SignalLossRho","SignalLossRho",3,
        {binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0Rho","SignalLoss0Rho",3,
        {binCentForMC,binPt,binSwitch},"s");

 CreateTHnSparse("SignalLossRhopipi","SignalLossRhopipi",3,
        {binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0Rhopipi","SignalLoss0Rhopipi",3,
        {binCentForMC,binPt,binSwitch},"s");

 CreateTHnSparse("SignalLossRhoMt","SignalLossRhoMt",3,
        {binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0RhoMt","SignalLoss0RhoMt",3,
        {binCentForMC,binPt,binSwitch},"s");
        
 CreateTHnSparse("SignalLossRhopipiMt","SignalLossRhopipiMt",3,
        {binCentForMC,binPt,binSwitch},"s");
 CreateTHnSparse("SignalLoss0RhopipiMt","SignalLoss0RhopipiMt",3,
        {binCentForMC,binPt,binSwitch},"s");


 CreateTHnSparse("hF0GenParticleFromPion","hF0GenParticleFromPion",4,
        {binZ,binCentForMC,binPt,binMass},"s");
 }
//********************************************

//Used Event Number**************
 fHistos->CreateTH1("hEvtNumberUsed","",1,0.5,1.5,"s");
 CreateTHnSparse("EvtSelector","EvtSelector",2,
	{binZ,binCent},"s");
//******************************



//Track Selection ****************
 fHistos->CreateTH2("PID_TPC_NSIG","",200,0,10,100,-10,10,"s");
 fHistos->CreateTH2("PID_TOF_NSIG","",200,0,10,100,-10,10,"s");
 fHistos->CreateTH2("PID_TPCalone_NSIG","",100,-10,10,200,0,10,"s");
//******************************


//Distribution for tracking efficiency correction
 if( fOption.Contains("MC") ){
 CreateTHnSparse("hRhoTrueParticle","hRhoTrueParticle",5,
        {binZ,binCentForMC,binPt,binMass,binTrackCutBit},"s");
 CreateTHnSparse("hF0TrueParticle","hF0TrueParticle",5,
        {binZ,binCentForMC,binPt,binMass,binTrackCutBit},"s");
 CreateTHnSparse("hF2TrueParticle","hF2TrueParticle",5,
        {binZ,binCentForMC,binPt,binMass,binTrackCutBit},"s");

 CreateTHnSparse("hRhoTrueParticleADDPID","hRhoTrueParticleADDPID",5,
        {binZ,binCentForMC,binPt,binMass,binTrackCutBit},"s");
 CreateTHnSparse("hF0TrueParticleADDPID","hF0TrueParticleADDPID",5,
        {binZ,binCentForMC,binPt,binMass,binTrackCutBit},"s");
 CreateTHnSparse("hF2TrueParticleADDPID","hF2TrueParticleADDPID",5,
        {binZ,binCentForMC,binPt,binMass,binTrackCutBit},"s");

 CreateTHnSparse("hRhoTrueParticleADDPIDTUNE","hRhoTrueParticleADDPIDTUNE",5,
        {binZ,binCentForMC,binPt,binMass,binTrackCutBit},"s");
 CreateTHnSparse("hF0TrueParticleADDPIDTUNE","hF0TrueParticleADDPIDTUNE",5,
        {binZ,binCentForMC,binPt,binMass,binTrackCutBit},"s");
 CreateTHnSparse("hF2TrueParticleADDPIDTUNE","hF2TrueParticleADDPIDTUNE",5,
        {binZ,binCentForMC,binPt,binMass,binTrackCutBit},"s");
 }
//************************************




//Distribution of MC Response**********
 fHistos->CreateTH2("PID_TPC_NSIG_MC","",200,0,10,100,-10,10,"s");
 fHistos->CreateTH2("PID_TOF_NSIG_MC","",200,0,10,100,-10,10,"s");

 fHistos->CreateTH2("PID_TPC_NSIG_MC_TUNE","",200,0,10,100,-10,10,"s");
 fHistos->CreateTH2("PID_TOF_NSIG_MC_TUNE","",200,0,10,100,-10,10,"s");

 CreateTHnSparse("hKSTrueParticleADDPID","hKSTrueParticleADDPID",4,
        {binZ,binCent,binPt,binMass},"s");
 CreateTHnSparse("hOmgTrueParticleADDPID","hOmgTrueParticleADDPID",4,
        {binZ,binCent,binPt,binMass},"s");

 CreateTHnSparse("hKSTrueParticleADDPIDTUNE","hKSTrueParticleADDPIDTUNE",4,
        {binZ,binCent,binPt,binMass},"s");
 CreateTHnSparse("hOmgTrueParticleADDPIDTUNE","hOmgTrueParticleADDPIDTUNE",4,
        {binZ,binCent,binPt,binMass},"s");
//***************************************



//PID Study**************************
 CreateTHnSparse("TPC_PID_After_TOF","TPC_PID_After_TOF",4,
	{binSigma,binTrackPt,binTrig,binParType},"s");
 CreateTHnSparse("TOF_PID_After_TPC","TOF_PID_After_TPC",4,
	{binSigma,binTrackPt,binTrig,binParType},"s");
//*************************************




//Fill Tracks*********************************
 CreateTHnSparse("hInvMass","InvMass",6,
	{binType,binZ,binCent,binPt,binMass,binTrackCutBit},"s");
 CreateTHnSparse("hInvMassMixing","InvMassMixing",6,
	{binType,binZ,binCent,binPt,binMass,binTrackCutBit},"s");
 CreateTHnSparse("hInvMassUnpair","InvMassUnpair",5,
        {binType,binZ,binCent,binPt,binMass},"s");


 CreateTHnSparse("KSTARRecParticle","KSTARRecParticle",5,
	{binType,binZ,binCent,binPt,binMass},"s");
 CreateTHnSparse("OmegaRecParticle","OmegaRecParticle",5,
	{binType,binZ,binCent,binPt,binMass},"s");
 CreateTHnSparse("PhiRecParticle","PhiRecParticle",5,
	{binType,binZ,binCent,binPt,binMass},"s");
 CreateTHnSparse("SigmaRecParticle","SigmaRecParticle",5,
        {binType,binZ,binCent,binPt,binMass},"s");
 CreateTHnSparse("LambdaRecParticle","LambdaRecParticle",5,
        {binType,binZ,binCent,binPt,binMass},"s");
 CreateTHnSparse("ExKaonRecParticle","ExKaonRecParticle",6,
	{binType,binZ,binCent,binPt,binMass,binExKaonNum},"s");
//**********************************************

//QA plots**************************************
 CreateTHnSparse("hSinglePion","hSinglePion",4,
	{binCent,binTrackPt,binCharge,binEta},"s");
//**********************************************

// fEMpool.resize(binCent.GetNbins(),
//	vector<eventpool> (binZ.GetNbins()));

 fEMpooltrk.resize(binTrackCutBit.GetNbins(),
	vector<vector<eventpool>>( binCent.GetNbins(),
	vector<eventpool>( binZ.GetNbins() ) ) );

 fOutput->Add(fHistos); 

 if( fOption.Contains("2018") ){
	fEventCuts.SetupPbPb2018();
 }
 else if( fOption.Contains("PbPb") ){
	fEventCuts.SetupRun2PbPb();
 }

 if( fOption.Contains("PbPb") ){
	fEventCuts.AddQAplotsToList(fHistos->GetListOfHistograms());
 }


 PostData(1, fHistos->GetListOfHistograms());

 fPIDCombined = new AliPIDCombined;
 fPIDCombined->SetDefaultTPCPriors();//Need more update..
 fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);
 fPIDCombined->SetDetectorMask(
	AliPIDResponse::kDetTPC |
        AliPIDResponse::kDetTOF |
        AliPIDResponse::kDetITS |
        AliPIDResponse::kDetTRD);//Do we need??
}


void AliAnalysisTaskInclusivef0f2::UserExec(Option_t *option)
{

 AliVEvent *event = InputEvent();
 if (!event) { Printf("ERROR: Could not retrieve event"); return; }

 event->IsA()==AliESDEvent::Class()
	? fEvt = dynamic_cast<AliESDEvent*>(event)
	: fEvt = dynamic_cast<AliAODEvent*>(event);
 if(!fEvt) return;

 bool IsEventSelectedPbPb = kFALSE;

 IsMC = kFALSE;
// if( IsFirstEvent ){
	Int_t runnumber = fEvt->GetRunNumber();
//	cout << runnumber << endl;
	if( fOption.Contains("MC") ) IsMC = kTRUE;
	fRunTable = new AliAnalysisTaskInclusivef0f2RunTable(runnumber);
	fRunTable->SetColl(0);
	if( fOption.Contains("pPb") ) fRunTable->SetColl(1);
	if( fOption.Contains("PbPb") ){
		fRunTable->SetColl(2);
	}

//	IsFirstEvent = kFALSE;
// }

 AliInputEventHandler* inputHandler = (AliInputEventHandler*)
        AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();

 fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
 if(!fPIDResponse){ printf("AliAnalysisTaskInclusivef0f2 No PIDd\n"); }

 fHistos -> FillTH1("hEventNumbers","All",1);

 fCent = 200;
 sel = (AliMultSelection*) fEvt -> FindListObject("MultSelection");
 AliVVZERO* lVV0 = fEvt->GetVZEROData();

 if( sel ){ fCent = sel->GetMultiplicityPercentile("V0M"); }
 if( fRunTable->IsPA() ) { fCent = sel->GetMultiplicityPercentile("V0A"); }
 if( fOption.Contains("UseZNA") ){ fCent = sel->GetMultiplicityPercentile("ZNA"); } 
 
 fEventCuts.OverrideAutomaticTriggerSelection( (AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral) );
 if( fOption.Contains("PbPb") && (
	( fCent > 10  && fCent < 30 ) || ( fCent > 50 ) ) ){
	fEventCuts.OverrideAutomaticTriggerSelection( (AliVEvent::kINT7) );
 }
 if( fOption.Contains("PbPb") ) IsEventSelectedPbPb = fEventCuts.AcceptEvent( event );

 double v0amplitude=0;
 for(int i=0;i<64;i++){ v0amplitude += lVV0->GetMultiplicity(i); }
 fMultiplicity = fEvt -> GetMultiplicity();


// const AliVVertex* trackVtx = fEvt->GetPrimaryVertexTPC(); //for ESD
 const AliVVertex* trackVtx = fEvt->GetPrimaryVertex();
 const AliVVertex* spdVtx = fEvt->GetPrimaryVertexSPD();
 fZ = -15.5;

// Generated True Particle distributions for efficiecny and acceptance correction
 Double_t genzvtx;
 if( IsMC && fEvt->IsA()==AliAODEvent::Class() && (inputHandler -> IsEventSelected()) & (AliVEvent::kINT7) ){
	fMCArray = (TClonesArray*) fEvt->FindListObject("mcparticles");
	AliAODMCHeader *cHeaderAOD  = dynamic_cast<AliAODMCHeader*>(fEvt->FindListObject(AliAODMCHeader::StdBranchName()));
	genzvtx = cHeaderAOD -> GetVtxZ();
	if( ( fabs(genzvtx)<10 && !fOption.Contains("SysZ") ) ||
	    ( fabs(genzvtx)<15 &&  fOption.Contains("SysZ") ) ){
		const Int_t nTracksMC = fMCArray->GetEntriesFast();
	        for(Int_t iTracks = 0; iTracks < nTracksMC; iTracks++){
	        	AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
	        	if( !trackMC ) continue;
			if( !trackMC->IsPrimary() ) continue;
	        	Int_t pdgCode = trackMC->PdgCode();
			if( pdgCode == 113 ){
				if( fabs( trackMC->Y() ) > 0.5 ) continue;
//				if( fRunTable->IsPA() && trackMC->Y() > 0 ) continue;
				if( fRunTable->IsPA() && ( trackMC->Y() > 0.035 || trackMC->Y() < -0.465 ) ) continue;
				if( trackMC->GetNDaughters() != 2 ) continue;
			        AliAODMCParticle* trackd1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ));
                                AliAODMCParticle* trackd2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ));
				if( !trackd2 ) continue;
				if( !trackd1 ) continue;
				if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->PdgCode() ) != 211 ) continue;
				if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->PdgCode() ) != 211 ) continue;
				FillTHnSparse("hRhoGenParticle",
					{genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass()},1.0 );
			}
			else if( pdgCode == 9010221 ){
				if( fabs( trackMC->Y() ) > 0.5 ) continue;
//				if( fRunTable->IsPA() && trackMC->Y() > 0 ) continue;
				if( fRunTable->IsPA() && ( trackMC->Y() > 0.035 || trackMC->Y() < -0.465 ) ) continue;
	                        if( trackMC->GetNDaughters() != 2 ) continue;
                                AliAODMCParticle* trackd1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ));
                                AliAODMCParticle* trackd2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ));
                                if( !trackd2 ) continue;
                                if( !trackd1 ) continue;
                                if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->PdgCode() ) != 211 ) continue;
                                if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->PdgCode() ) != 211 ) continue;
				FillTHnSparse("hF0GenParticle",
					{genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass()},1.0 );

				FillTHnSparse("hF0GenParticleFromPion",
					{genzvtx,fCent,
					sqrt( pow( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->Px() +
						   dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->Px(), 2) +
					      pow( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->Py() +
						   dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->Py(), 2) ),

					sqrt( pow( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->E() +
						   dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->E(), 2) - 
					      pow( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->Px() +
						   dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->Px(), 2) -
					      pow( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->Py() +
						   dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->Py(), 2) -
					      pow( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->Pz() +
						   dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->Pz(), 2) )}, 1.0 );

			}
			else if( pdgCode == 225 ){
				if( fabs( trackMC->Y() ) > 0.5 ) continue;
//				if( fRunTable->IsPA() && trackMC->Y() > 0 ) continue;
				if( fRunTable->IsPA() && ( trackMC->Y() > 0.035 || trackMC->Y() < -0.465 ) ) continue;  
                               if( trackMC->GetNDaughters() != 2 ) continue;
                                AliAODMCParticle* trackd1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ));
                                AliAODMCParticle* trackd2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ));
                                if( !trackd2 ) continue;
                                if( !trackd1 ) continue;
                                if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->PdgCode() ) != 211 ) continue;
                                if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->PdgCode() ) != 211 ) continue;
				FillTHnSparse("hF2GenParticle",
					{genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass()},1.0 );
			}
	        }
	}
 }
//***********************************************
//***********************************************

 Bool_t IsTriggered = kFALSE;
 Bool_t IsNotPileup = kFALSE;
 Bool_t IsValidVtx = kFALSE;
 Bool_t IsGoodVtx = kFALSE;
 Bool_t IsINT7LTZERO = kFALSE;
 Bool_t IsSelectedFromAliMultSelection = kFALSE;
 Bool_t IsMultiplicityInsideBin = kFALSE;
//Event Selection Criteria


 Bool_t IsSelectedFromAliMultSelectionForSysZ = kFALSE;
//

//IsTriggered Flag Configuration

 IsTriggered = (inputHandler -> IsEventSelected()) & (AliVEvent::kINT7);
 if( fOption.Contains("kMB") ) IsTriggered = (inputHandler -> IsEventSelected()) & (AliVEvent::kMB);
 if( fOption.Contains("HighMult") ) IsTriggered = inputHandler -> IsEventSelected() & AliVEvent::kHighMultV0;
// if( fOption.Contains("2018") ) IsTriggered = fEvt -> GetTriggerMask() & (AliVEvent::kINT7); 
 if( fOption.Contains("2018") ) IsTriggered = fEventCuts.PassedCut(AliEventCuts::kTrigger);
 if( IsMC ) IsTriggered = inputHandler -> IsEventSelected() & AliVEvent::kINT7;
//*****************************


//IsNotPileup Flag Configuration
 if( IsMC ) IsNotPileup = kTRUE;
 else if( fRunTable->IsAA() ) IsNotPileup = kTRUE;
 else if( !IsMC && !event->IsPileupFromSPDInMultBins() &&
	(  fRunTable->IsPP() || fRunTable->IsPA() ) ) IsNotPileup = kTRUE;

 if( fOption.Contains("NoPileupCut") ) IsNotPileup = kTRUE;
//*****************************


//IsGoodVtx Flag Configuration
 if( spdVtx ){
	if( spdVtx->GetNContributors() > 0.5 ) IsValidVtx = kTRUE;
	fZ = spdVtx->GetZ();
	zbin = binZ.FindBin(fZ) -1;
	if( fabs(fZ) < 10.0 && !(zbin < 0 )) IsGoodVtx = kTRUE;
 }
//***************************************



//IsINT7LTZERO *******************
 for(int i=0;i<fMultiplicity->GetNumberOfTracklets();i++){
	if( fMultiplicity->GetEta(i) < 1.0 ) IsINT7LTZERO = kTRUE;
 }
//******************************



//IsSelectedFromAliMultSelection Flag Configuration
 if( !IsMC && !fOption.Contains("HighMult") && !fOption.Contains("kMB") &&
//	( fRunTable->IsPP() || fRunTable->IsPA() || fRunTable->IsAA() ) ){
	( fRunTable->IsPP() || fRunTable->IsPA() ) ){
	if( sel->IsEventSelected() ) IsSelectedFromAliMultSelection = kTRUE;

	if( sel->GetThisEventIsNotPileupInMultBins() &&
	sel->GetThisEventINELgtZERO() &&
	sel->GetThisEventPassesTrackletVsCluster() &&
	sel->GetThisEventHasNoInconsistentVertices() &&
	fabs(fZ) < 15 ){
		IsSelectedFromAliMultSelectionForSysZ = kTRUE;
	}
 }
 if( fOption.Contains("HighMult") ){
        if( sel->GetThisEventIsNotPileup() &&
        sel->GetThisEventIsNotPileupInMultBins() &&
        sel->GetThisEventHasNoInconsistentVertices() &&
        sel->GetThisEventPassesTrackletVsCluster() ){
		IsSelectedFromAliMultSelection = kTRUE;
		if( fabs(fZ) < 15 ){
			IsSelectedFromAliMultSelectionForSysZ = kTRUE;
		}
	}
 }
 if( fRunTable->IsAA() ) IsSelectedFromAliMultSelection = kTRUE;
 if( fRunTable->IsAA() ) IsSelectedFromAliMultSelectionForSysZ = kTRUE;

 if( IsMC ) IsSelectedFromAliMultSelection = kTRUE;
 if( IsMC && fabs(genzvtx)<15 ) IsSelectedFromAliMultSelectionForSysZ = kTRUE;
 if( fOption.Contains("kMB") ) IsSelectedFromAliMultSelection = kTRUE;
 if( fOption.Contains("INEL")) IsSelectedFromAliMultSelection = kTRUE;
//***************************************************


//IsMultiplicityInsideBin Flag Configuration********
 centbin = binCent.FindBin(fCent) -1;
 if( centbin >= 0 || fOption.Contains("INEL")) IsMultiplicityInsideBin = kTRUE; 
//************************************



 if( IsTriggered ) fHistos->FillTH1("hEventNumbers","IsTriggered",1);
 if( IsTriggered && IsNotPileup ) fHistos->FillTH1("hEventNumbers","IsNotPileup",1);
 if( IsTriggered && IsNotPileup && IsValidVtx ) fHistos->FillTH1("hEventNumbers","IsValidVtx",1);
 if( IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx ) fHistos->FillTH1("hEventNumbers","IsGoodVtx",1);
 if( IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx && IsSelectedFromAliMultSelection ) fHistos->FillTH1("hEventNumbers","IsSelectedFromAliMultSelection",1);
 if( ( IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin && !fOption.Contains("2018") ) ||
	( IsEventSelectedPbPb && fOption.Contains("2018") ) ){
	fHistos->FillTH1("hEventNumbers","IsMultiplicityInsideBin",1);
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

 if( fOption.Contains("QAMode") && IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin ){
	if( this -> GoodTracksSelection(0x20, 5, 3, 2,0.01) ) FillTHnSparse("EvtSelector",{fZ,fCent},1.0);
 }

 if( !fOption.Contains("EvtSelStudy") ){
	if( !fOption.Contains("Sys") ){
		if( ( IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin && !fOption.Contains("2018") ) ||
			( IsEventSelectedPbPb && fOption.Contains("2018") ) ){
			if( fOption.Contains("MismatchCheck") ){
				if(this -> GoodTracksSelection(0x20, 5, 3, 2, 100)) this -> FillTracks();
			} else{ if(this -> GoodTracksSelection(0x20, 5, 3, 2, 0.01)) this -> FillTracks(); }
			fHistos->FillTH1("hEvtNumberUsed",1,1);
			FillTHnSparse("EvtSelector",{fZ,fCent},1.0);
		}
	}
	else if( fOption.Contains("SysZ") ){
		if( !fOption.Contains("SysTrk") ){
			if( IsTriggered && IsNotPileup && IsValidVtx && ( fabs(fZ) < 15 ) &&
			IsSelectedFromAliMultSelectionForSysZ && IsMultiplicityInsideBin ){
				if(this -> GoodTracksSelection(0x20, 5, 3, 2,0.01)) this -> FillTracks();
				fHistos->FillTH1("hEvtNumberUsed",1,1);
				FillTHnSparse("EvtSelector",{fZ,fCent},1.0);
			}
		}
		else if( fOption.Contains("SysTrk") ){
			if( !fOption.Contains("SysPID") ){
				if( IsTriggered && IsNotPileup && IsValidVtx && ( fabs(fZ) < 15 ) &&
				IsSelectedFromAliMultSelectionForSysZ && IsMultiplicityInsideBin ){
					if(this -> GoodTracksSelection(0x20, 5, 3, 2,0.01)) this -> FillTracks();
					if(this -> GoodTracksSelection(0x60, 5, 3, 2,0.01)) this -> FillTracks();
					if(this -> GoodTracksSelection(0x300, 5, 3, 2,0.01)) this -> FillTracks();
					fHistos->FillTH1("hEvtNumberUsed",1,1);
					FillTHnSparse("EvtSelector",{fZ,fCent},1.0);
				}
			}
			else if( fOption.Contains("SysPID") ){
				if( IsTriggered && IsNotPileup && IsValidVtx && ( fabs(fZ) < 15 ) &&
				IsSelectedFromAliMultSelectionForSysZ && IsMultiplicityInsideBin ){
					if(this -> GoodTracksSelection(0x20, 5, 3, 2,0.01)) this -> FillTracks();
					if(this -> GoodTracksSelection(0x60, 5, 3, 2,0.01)) this -> FillTracks();
					if(this -> GoodTracksSelection(0x300, 5, 3, 2,0.01)) this -> FillTracks();
					if(this -> GoodTracksSelection(0x20, 5, 3.5, 2,0.01)) this -> FillTracks();
					if(this -> GoodTracksSelection(0x20, 5, 3, 2.5,0.01)) this -> FillTracks();
					if(this -> GoodTracksSelection(0x20, 5, 2.5, 2,0.01)) this -> FillTracks();
					if(this -> GoodTracksSelection(0x20, 5, 3, 1.5,0.01)) this -> FillTracks();

					fHistos->FillTH1("hEvtNumberUsed",1,1);
					FillTHnSparse("EvtSelector",{fZ,fCent},1.0);
				}
			}
		}
	}
 }


//Trigger Efficiency**********************
 bool IsINEL = false;
 bool IsINEL0 = true;
 if( IsMC ){
        if( fEvt->IsA()==AliAODEvent::Class() ){
                fMCArray = (TClonesArray*) fEvt->FindListObject("mcparticles");
                AliAODMCHeader *cHeaderAOD  = dynamic_cast<AliAODMCHeader*>
                        (fEvt->FindListObject(AliAODMCHeader::StdBranchName()));
		const Int_t nTracksMC = fMCArray->GetEntries();  
                for(Int_t iTracks = 0; iTracks < nTracksMC; iTracks++){
                        AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
                        if( !trackMC ) continue;
                        if( trackMC->Charge() == 0 ) continue;
                       	if( fabs( trackMC->Eta() ) > 1.0 ) continue;
                        IsINEL = true;
                }
        }
 }      

 if( IsINEL && IsMC ){
//        if( (inputHandler -> IsEventSelected()) & (AliVEvent::kINT7) ){
	if( IsTriggered ){
                FillTHnSparse("TrigEffMult",{fCent,1.0},1.0 );
        }
        else{
                FillTHnSparse("TrigEffMult",{fCent,0.0},1.0 );
        }
 }

 if( IsINEL0 && IsMC ){
//        if( (inputHandler -> IsEventSelected()) & (AliVEvent::kINT7) ){
	if( IsTriggered ){
                FillTHnSparse("TrigEffMult0",{fCent,1.0},1.0 );
        }
        else{
                FillTHnSparse("TrigEffMult0",{fCent,0.0},1.0 );
        }
 }
//************************************8

//Corrections for event selection******************
 if( IsMC && IsTriggered && IsINEL ){
//	if( IsValidVtx && sel->GetThisEventHasGoodVertex2016() && sel->GetThisEventHasNoInconsistentVertices() && IsGoodVtx ){
	if( IsValidVtx && IsGoodVtx ){
		FillTHnSparse("VtxSelection",{fCent,1.0}, 1.0);
	}
	else{ FillTHnSparse("VtxSelection",{fCent,0.0}, 1.0); }
 }

 if( IsMC && IsTriggered && IsINEL0 ){
//        if( IsValidVtx && sel->GetThisEventHasGoodVertex2016() && sel->GetThisEventHasNoInconsistentVertices() && IsGoodVtx ){
	if( IsValidVtx && IsGoodVtx ){
                FillTHnSparse("VtxSelection0",{fCent,1.0}, 1.0);
        }
        else{ FillTHnSparse("VtxSelection0",{fCent,0.0}, 1.0); }
 }

 if( IsMC ){
	const Int_t nTracksMC = fMCArray->GetEntriesFast();
	for(Int_t iTracks = 0; iTracks < nTracksMC; iTracks++){
		AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
		if( !trackMC ) continue;
		Int_t pdgCode = trackMC->PdgCode();
		if( pdgCode == 9010221 ){
			if( fabs( trackMC->Y() ) > 0.5 ) continue;
//			if( fRunTable->IsPA() && trackMC->Y() > 0 ) continue;
			if( fRunTable->IsPA() && ( trackMC->Y() > 0.035 || trackMC->Y() < -0.465 ) ) continue;
			if( trackMC->GetNDaughters() != 2 ) continue;
			AliAODMCParticle* trackd1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ));
			AliAODMCParticle* trackd2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ));
			if( !trackd2 ) continue; if( !trackd1 ) continue;
			if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->PdgCode() ) != 211 ) continue;
			if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->PdgCode() ) != 211 ) continue;
			FillTHnSparse("SignalLoss0",{fCent,trackMC->Pt(),0.0},1.0);
 			if( IsINEL>0 )FillTHnSparse("SignalLoss",{fCent,trackMC->Pt(),0.0},1.0);
			if( IsTriggered && IsValidVtx && sel->GetThisEventHasGoodVertex2016() && sel->GetThisEventHasNoInconsistentVertices() ){
//			if( IsTriggered && IsNotPileup && IsValidVtx && IsGoodVtx && IsSelectedFromAliMultSelection && IsMultiplicityInsideBin ){
				if( IsINEL>0 ) FillTHnSparse("SignalLoss",{fCent,trackMC->Pt(),1.0},1.0);
				FillTHnSparse("SignalLoss0",{fCent,trackMC->Pt(),1.0},1.0);
			}
		}
		else if( abs(pdgCode) == 211 ){
			if( fabs( trackMC->Y() ) > 0.5 ) continue;
//			if( fRunTable->IsPA() && trackMC->Y() > 0 ) continue;
			if( fRunTable->IsPA() && ( trackMC->Y() > 0.035 || trackMC->Y() < -0.465 ) ) continue;
			FillTHnSparse("SignalLoss0Pion",{fCent,trackMC->Pt(),0.0,(double)trackMC->Charge()/3.0},1.0 );
			FillTHnSparse("SignalLoss0PionMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0,(double)trackMC->Charge()/3.0},1.0 );
			if( IsINEL>0 ){
				FillTHnSparse("SignalLossPion",{fCent,trackMC->Pt(),0.0,(double)trackMC->Charge()/3.0},1.0 );	
				FillTHnSparse("SignalLossPionMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0,(double)trackMC->Charge()/3.0},1.0 );
			}
			if( IsTriggered && IsValidVtx && sel->GetThisEventHasGoodVertex2016() && sel->GetThisEventHasNoInconsistentVertices() ){
				FillTHnSparse("SignalLoss0Pion",{fCent,trackMC->Pt(),1.0,(double)trackMC->Charge()/3.0},1.0 );
				FillTHnSparse("SignalLoss0PionMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0,(double)trackMC->Charge()/3.0},1.0 );
				if( IsINEL>0 ){
					FillTHnSparse("SignalLossPion",{fCent,trackMC->Pt(),1.0,(double)trackMC->Charge()/3.0},1.0 );
					FillTHnSparse("SignalLossPionMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0,(double)trackMC->Charge()/3.0},1.0 );
				}
			}
		}
/*
		else if( pdgCode == 310 || pdgCode == 113 ){
                        if( fabs( trackMC->Y() ) > 0.5 ) continue;
                        if( fRunTable->IsPA() && trackMC->Y() > 0 ) continue;
			if( pdgCode == 310 ){
                        	FillTHnSparse("SignalLoss0Kaon",{fCent,trackMC->Pt(),0.0},1.0 );
				FillTHnSparse("SignalLoss0KaonMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0},1.0 );
                        	if( IsINEL>0 ){
					FillTHnSparse("SignalLossKaon",{fCent,trackMC->Pt(),0.0},1.0 );
					FillTHnSparse("SignalLossKaonMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0},1.0 );
				}
                        	if( IsTriggered && IsValidVtx && sel->GetThisEventHasGoodVertex2016() && sel->GetThisEventHasNoInconsistentVertices() ){
                        	        FillTHnSparse("SignalLoss0Kaon",{fCent,trackMC->Pt(),1.0},1.0 );
					FillTHnSparse("SignalLoss0KaonMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0},1.0 );
                        	        if( IsINEL>0 ){
						FillTHnSparse("SignalLossKaon",{fCent,trackMC->Pt(),1.0},1.0 );
						FillTHnSparse("SignalLossKaonMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0},1.0 );
					}
                        	}
			}
			else if( pdgCode == 113 ){
                                FillTHnSparse("SignalLoss0Rho",{fCent,trackMC->Pt(),0.0},1.0 );
                                FillTHnSparse("SignalLoss0RhoMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0},1.0 );
                                if( IsINEL>0 ){
                                        FillTHnSparse("SignalLossRho",{fCent,trackMC->Pt(),0.0},1.0 );
                                        FillTHnSparse("SignalLossRhoMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0},1.0 );
                                }
                                if( IsTriggered && IsValidVtx && sel->GetThisEventHasGoodVertex2016() && sel->GetThisEventHasNoInconsistentVertices() ){
                                        FillTHnSparse("SignalLoss0Rho",{fCent,trackMC->Pt(),1.0},1.0 );
                                        FillTHnSparse("SignalLoss0RhoMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0},1.0 );
                                        if( IsINEL>0 ){
                                                FillTHnSparse("SignalLossRho",{fCent,trackMC->Pt(),1.0},1.0 );
                                                FillTHnSparse("SignalLossRhoMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0},1.0 );
                                        }
                                }
                        }
                        if( trackMC->GetNDaughters() != 2 ) continue;
                        AliAODMCParticle* trackd1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ));
                        AliAODMCParticle* trackd2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ));
                        if( !trackd2 ) continue; if( !trackd1 ) continue;
                        if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(0) ))->PdgCode() ) != 211 ) continue;
                        if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( trackMC->GetDaughterLabel(1) ))->PdgCode() ) != 211 ) continue;
			if( pdgCode == 310 ){
				FillTHnSparse("SignalLoss0Kaonpipi",{fCent,trackMC->Pt(),0.0}, 1.0 );
				FillTHnSparse("SignalLoss0KaonpipiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0},1.0 );
				if( IsINEL>0 ){
                        	        FillTHnSparse("SignalLossKaonpipi",{fCent,trackMC->Pt(),0.0},1.0 );
                        	        FillTHnSparse("SignalLossKaonpipiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0},1.0 );
				}
                        	if( IsTriggered && IsValidVtx && sel->GetThisEventHasGoodVertex2016() && sel->GetThisEventHasNoInconsistentVertices() ){
                        	        FillTHnSparse("SignalLoss0Kaonpipi",{fCent,trackMC->Pt(),1.0},1.0 );
                        	        FillTHnSparse("SignalLoss0KaonpipiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0},1.0 );
                        	        if( IsINEL>0 ){
                        	                FillTHnSparse("SignalLossKaonpipi",{fCent,trackMC->Pt(),1.0},1.0 );
                        	                FillTHnSparse("SignalLossKaonpipiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0},1.0 );
                        	        }
                        	}
			}
			else if( pdgCode == 113 ){
				FillTHnSparse("SignalLoss0Rhopipi",{fCent,trackMC->Pt(),0.0}, 1.0 );
                                FillTHnSparse("SignalLoss0RhopipiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0},1.0 );
                                if( IsINEL>0 ){
                                        FillTHnSparse("SignalLossRhopipi",{fCent,trackMC->Pt(),0.0},1.0 );
                                        FillTHnSparse("SignalLossRhopipiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0},1.0 );
                                }
                                if( IsTriggered && IsValidVtx && sel->GetThisEventHasGoodVertex2016() && sel->GetThisEventHasNoInconsistentVertices() ){
                                        FillTHnSparse("SignalLoss0Rhopipi",{fCent,trackMC->Pt(),1.0},1.0 );
                                        FillTHnSparse("SignalLoss0RhopipiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0},1.0 );
                                        if( IsINEL>0 ){
                                                FillTHnSparse("SignalLossRhopipi",{fCent,trackMC->Pt(),1.0},1.0 );
                                                FillTHnSparse("SignalLossRhopipiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0},1.0 );
                                        }
                                }
			}
		}
		else if( pdgCode == 333 ){
                        if( fabs( trackMC->Y() ) > 0.5 ) continue;
                        if( fRunTable->IsPA() && trackMC->Y() > 0 ) continue;
                        FillTHnSparse("SignalLoss0Phi",{fCent,trackMC->Pt(),0.0},1.0 );
			FillTHnSparse("SignalLoss0PhiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0},1.0 );
                        if( IsINEL>0 ){
				FillTHnSparse("SignalLossPhi",{fCent,trackMC->Pt(),0.0},1.0 );
				FillTHnSparse("SignalLossPhiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),0.0},1.0 );
			}
                        if( IsTriggered && IsValidVtx && sel->GetThisEventHasGoodVertex2016() && sel->GetThisEventHasNoInconsistentVertices() ){
                                FillTHnSparse("SignalLoss0Phi",{fCent,trackMC->Pt(),1.0},1.0 );
				FillTHnSparse("SignalLoss0PhiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0},1.0 );
                                if( IsINEL>0 ){
					FillTHnSparse("SignalLossPhi",{fCent,trackMC->Pt(),1.0},1.0 );
					FillTHnSparse("SignalLossPhiMt",{fCent,sqrt( pow(trackMC->Pt(),2) + pow(trackMC->M(),2 ) ),1.0},1.0 );
				}
                        }
		}
*/
	}
 }
//************************************************


 PostData(1, fHistos->GetListOfHistograms());
}


bool AliAnalysisTaskInclusivef0f2::GoodTracksSelection(int trkcut, double TPCsig, double TOFsig, double TPCalonesig, double TOFMismatchRatio=0.01){

 const UInt_t ntracks = fEvt ->GetNumberOfTracks();
 goodtrackindices.clear();

 AliVTrack* track;

 tracklist* etl;
 eventpool* ep;

 fFilterBit = trkcut;
 trkbin = 0;

 if( fOption.Contains("SysTrk") ){
	if( trkcut == 0x20 ){
		trkbin = 0;
	}
	else if( trkcut == 0x60 ){
		trkbin = 1;
	}
	else if( trkcut == 0x300 ){
		trkbin = 2;
	}
	fFilterBit = trkcut;
 }
 if( fOption.Contains("SysPID") ){  
	if( TOFsig == 3.5 ){
		trkbin = 3;
	}
	if( TPCalonesig == 2.5 ){
		trkbin = 4;
	}
	if( TOFsig == 2.5 ){
		trkbin = 5;
	}
	if( TPCalonesig == 1.5 ){
		trkbin = 6;
	}
 }


 if( centbin>=0 && zbin>=0 && fRunTable->IsPP() && fOption.Contains("AddMixing") ){
	ep = &fEMpooltrk[trkbin][centbin][zbin];
	ep -> push_back( tracklist() ); 
	etl = &(ep->back());
 }
 for(UInt_t it=0;it<ntracks;it++){
	if( fEvt->IsA()==AliESDEvent::Class() ){
//
//
	}
	else{
		track = (AliAODTrack*)fEvt->GetTrack(it);
		if( !track ) continue;
		if( trkbin != 1 && trkbin != 2 ) if( !((AliAODTrack*) track)->TestFilterMask(fFilterBit) ) continue; //for global
		if( trkbin == 1 || trkbin == 2 ) if( !((AliAODTrack*) track)->TestFilterBit(fFilterBit) ) continue; //for hybrid, SSD global

		if( trkbin == 0 ){
			fHistos -> FillTH2("PID_TPC_NSIG",track->Pt(),fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion),1.0);
			fHistos -> FillTH2("PID_TOF_NSIG",track->Pt(),fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion),1.0);
		}
/*
		if( fPIDResponse->GetTOFMismatchProbability( track ) < 0.01 ){
			if( fabs( fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion) ) > TOFsig ) continue;
		}
		if( fabs( fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion) ) > TOFsig ||
		fPIDResponse->GetTOFMismatchProbability( track ) > 0.01 ){
			if( fabs( fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion) ) > TPCalonesig ) continue;
		}
*/
		if( !(
			( fPIDResponse->GetTOFMismatchProbability( track ) < TOFMismatchRatio &&
			fabs( fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion) ) < TOFsig ) ||

//			( fPIDResponse->GetTOFMismatchProbability( track ) > 0.01 &&
			( fabs( fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion) ) < TPCalonesig ) )
		) continue;

		if( track->Pt() < fptcut ) continue;
		if( fabs( track->Eta() ) > fetacut ) continue;

		FillTHnSparse("hSinglePion", {fCent,track->Pt(),static_cast<double>(track->Charge()),track->Eta()},1.0 );
	}

	goodtrackindices.push_back(it);
		
	if( fRunTable->IsPP() && fOption.Contains("AddMixing") ) etl->push_back( (AliVTrack*) track -> Clone() );
 }

 double genzvtx;
 int Label1, Label2, Label3;
 int trkl1, trkl2, trkl3;
 int trkid1, trkid2, trkid3;
 int PIDcut1, PIDcut2;
 int trk1count=0;
 int trk2count=0;
 int trk3count=0;

 if( IsMC && fEvt->IsA()==AliAODEvent::Class() ){
	fMCArray = (TClonesArray*) fEvt->FindListObject("mcparticles");
	AliAODMCHeader *cHeaderAOD  = dynamic_cast<AliAODMCHeader*>
		(fEvt->FindListObject(AliAODMCHeader::StdBranchName()));
	genzvtx = cHeaderAOD -> GetVtxZ();
        if( ( fabs(genzvtx)<10 && !fOption.Contains("SysZ") ) ||
            ( fabs(genzvtx)<15 &&  fOption.Contains("SysZ") ) ){
		const Int_t nTracksMC = fMCArray->GetEntriesFast();
		for(Int_t iTracks = 0; iTracks < nTracksMC; iTracks++){
			AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
			if( !trackMC ) continue;
			if( !trackMC->IsPrimary() ) continue;
			Int_t pdgCode = trackMC->GetPdgCode();
			if( fabs( trackMC->Y() ) > 0.5 ) continue;
//			if( fRunTable->IsPA() && trackMC->Y() > 0 ) continue;
			if( fRunTable->IsPA() && ( trackMC->Y() > 0.035 || trackMC->Y() < -0.465 ) ) continue;
			if( pdgCode == 113 || pdgCode == 9010221 || pdgCode == 225 ){
				if( trackMC->GetNDaughters() != 2 ) continue;
				Label1 = trackMC->GetDaughterLabel(0);
				Label2 = trackMC->GetDaughterLabel(1);
//				Label1 = trackMC->GetDaughter(0);
  //                              Label2 = trackMC->GetDaughter(1);

//				if( Label1<0 || Label2<0 ) continue;

				AliAODMCParticle* trackd1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( Label1 ));
				AliAODMCParticle* trackd2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( Label2 ));

				if( !trackd1 ) continue;
				if( !trackd2 ) continue;

				if( abs( trackd1->PdgCode() ) != 211 ) continue;
				if( abs( trackd2->PdgCode() ) != 211 ) continue;

                                if( trackd1->Pt() < fptcut ) continue;
                                if( fabs( trackd1->Eta() ) > fetacut ) continue;
                                if( trackd1->Pt() < fptcut ) continue;
                                if( fabs( trackd1->Eta() ) > fetacut ) continue;

                                if( pdgCode == 113 )
                                FillTHnSparse("hRhoTrueParticle",
                                        {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass(),(double)(trkbin+1)},1.0 );
                                else if( pdgCode == 9010221 )
                                FillTHnSparse("hF0TrueParticle",
                                        {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass(),(double)(trkbin+1)},1.0 );
                                else if( pdgCode == 225 )
                                FillTHnSparse("hF2TrueParticle",
                                        {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass(),(double)(trkbin+1)},1.0 );


				trkl1 = trackd1->GetLabel();
				trkl2 = trackd2->GetLabel();
				
				trk1count=0;
				trk2count=0;

				AliAODTrack* trackd1Recon;
				AliAODTrack* trackd2Recon;

				trkid1 = -1;
				trkid2 = -1;


				for(UInt_t it=0;it<ntracks;it++){
					track = (AliAODTrack*)fEvt->GetTrack(it);
					if( !track ) continue;
					if( pdgCode == 9010221 ) cout << track->GetLabel() << ", " << track->GetMother() << endl;
					if( (int)track->GetLabel() == Label1 || (int)track->GetLabel() == Label2 ){
						trackd1Recon = (AliAODTrack*)fEvt->GetTrack( it );
						if( trackd1Recon ){
							trkid1=it;
							break;				
						}		
					}
				}
				if( trkid1 == -1 ){
					continue;
				}

				for(UInt_t it=0;it<ntracks;it++){
					track = (AliAODTrack*)fEvt->GetTrack(it);
					if( !track ) continue;
					if( (int)track->GetLabel() == Label1 || (int)track->GetLabel() == Label2 ){
					if( it == trkid1) continue;
						trackd2Recon = (AliAODTrack*)fEvt->GetTrack( it );
						if( trackd2Recon ){
                                                        trkid2=it;
                                                        break;
                                                }
					}
				}
				if( trkid2 == -1 ){
                                        continue;
				}


				if( trackd1Recon->Pt() < fptcut ) continue;
				if( fabs( trackd1Recon->Eta() ) > fetacut ) continue;
				if( trkbin != 1 && trkbin != 2 ) if( !((AliAODTrack*) trackd1Recon)->TestFilterMask(fFilterBit) ) continue;
				if( trkbin == 1 || trkbin == 2 ) if( !((AliAODTrack*) trackd1Recon)->TestFilterBit(fFilterBit) ) continue;

                                if( trackd2Recon->Pt() < fptcut ) continue;
                                if( fabs( trackd2Recon->Eta() ) > fetacut ) continue;
				if( trkbin != 1 && trkbin != 2 ) if( !((AliAODTrack*) trackd2Recon)->TestFilterMask(fFilterBit) ) continue;
				if( trkbin == 1 || trkbin == 2 ) if( !((AliAODTrack*) trackd2Recon)->TestFilterBit(fFilterBit) ) continue;
	

				if( fPIDResponse->GetTOFResponse().IsA() ){
					fPIDResponse->SetTunedOnData(false);
				       	fHistos -> FillTH2("PID_TOF_NSIG_MC",trackd1Recon->Pt(),fPIDResponse->NumberOfSigmasTOF(trackd1Recon, AliPID::kPion),1.0);
	                                fHistos -> FillTH2("PID_TOF_NSIG_MC",trackd2Recon->Pt(),fPIDResponse->NumberOfSigmasTOF(trackd2Recon, AliPID::kPion),1.0);
					fPIDResponse->SetTunedOnData(true);
	                                fHistos -> FillTH2("PID_TOF_NSIG_MC_TUNE",trackd1Recon->Pt(),fPIDResponse->NumberOfSigmasTOF(trackd1Recon, AliPID::kPion),1.0);
	                                fHistos -> FillTH2("PID_TOF_NSIG_MC_TUNE",trackd2Recon->Pt(),fPIDResponse->NumberOfSigmasTOF(trackd2Recon, AliPID::kPion),1.0);
				}
				fPIDResponse->SetTunedOnData(false);
				fHistos -> FillTH2("PID_TPC_NSIG_MC",trackd1Recon->Pt(),fPIDResponse->NumberOfSigmasTPC(trackd1Recon, AliPID::kPion),1.0);
				fHistos -> FillTH2("PID_TPC_NSIG_MC",trackd2Recon->Pt(),fPIDResponse->NumberOfSigmasTPC(trackd2Recon, AliPID::kPion),1.0);
				fPIDResponse->SetTunedOnData(true);
				fHistos -> FillTH2("PID_TPC_NSIG_MC_TUNE",trackd1Recon->Pt(),fPIDResponse->NumberOfSigmasTPC(trackd1Recon, AliPID::kPion),1.0);
				fHistos -> FillTH2("PID_TPC_NSIG_MC_TUNE",trackd2Recon->Pt(),fPIDResponse->NumberOfSigmasTPC(trackd2Recon, AliPID::kPion),1.0);

				fPIDResponse->SetTunedOnData(false);

				PIDcut1 = 0;
				PIDcut2 = 0;

				if( 
				( fPIDResponse->GetTOFMismatchProbability( trackd1Recon ) < TOFMismatchRatio
				&& fabs( fPIDResponse->NumberOfSigmasTOF(trackd1Recon, AliPID::kPion) ) < TOFsig ) ||
//				( fPIDResponse->GetTOFMismatchProbability( trackd1Recon ) > 0.01
//				&&
				( fabs( fPIDResponse->NumberOfSigmasTPC(trackd1Recon, AliPID::kPion) ) < TPCalonesig ) ){
					PIDcut1=1;
				}

                                if(
				( fPIDResponse->GetTOFMismatchProbability( trackd2Recon ) < TOFMismatchRatio
                                && fabs( fPIDResponse->NumberOfSigmasTOF(trackd2Recon, AliPID::kPion) ) < TOFsig ) ||
//				( fPIDResponse->GetTOFMismatchProbability( trackd2Recon ) > 0.01
//                                &&
				( fabs( fPIDResponse->NumberOfSigmasTPC(trackd2Recon, AliPID::kPion) ) < TPCalonesig ) ){
                                        PIDcut2=1;
                                }

				if( PIDcut1 && PIDcut2 ){

                                	if( pdgCode == 113 ){
                                		FillTHnSparse("hRhoTrueParticleADDPID",
                                		        {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass(),(double)(trkbin+1)},1.0 );
					}
                                	else if( pdgCode == 9010221 ){
                                		FillTHnSparse("hF0TrueParticleADDPID",
                                		        {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass(),(double)(trkbin+1)},1.0 );
					}
                                	else if( pdgCode == 225 ){
                                		FillTHnSparse("hF2TrueParticleADDPID",
                                		        {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass(),(double)(trkbin+1)},1.0 );
					}
				}
                                fPIDResponse->SetTunedOnData(true);

                                PIDcut1 = 0;
                                PIDcut2 = 0;

                                if(
				( fPIDResponse->GetTOFMismatchProbability( trackd1Recon ) < TOFMismatchRatio
                                && fabs( fPIDResponse->NumberOfSigmasTOF(trackd1Recon, AliPID::kPion) ) < TOFsig ) ||
//				( fPIDResponse->GetTOFMismatchProbability( trackd1Recon ) > 0.01
//                                &&
				( fabs( fPIDResponse->NumberOfSigmasTPC(trackd1Recon, AliPID::kPion) ) < TPCalonesig ) ){
                                        PIDcut1=1;
                                }

                                if(
				( fPIDResponse->GetTOFMismatchProbability( trackd2Recon ) < TOFMismatchRatio
                                && fabs( fPIDResponse->NumberOfSigmasTOF(trackd2Recon, AliPID::kPion) ) < TOFsig ) ||
//				( fPIDResponse->GetTOFMismatchProbability( trackd2Recon ) > 0.01
  //                              &&
				( fabs( fPIDResponse->NumberOfSigmasTPC(trackd2Recon, AliPID::kPion) ) < TPCalonesig ) ){
                                        PIDcut2=1;
                                }

                                if( PIDcut1 && PIDcut2 ){
                                        if( pdgCode == 113 )
                                        FillTHnSparse("hRhoTrueParticleADDPIDTUNE",
                                                {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass(),(double)(trkbin+1)},1.0 );
                                        else if( pdgCode == 9010221 )
                                        FillTHnSparse("hF0TrueParticleADDPIDTUNE",
                                                {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass(),(double)(trkbin+1)},1.0 );
                                        else if( pdgCode == 225 )
                                        FillTHnSparse("hF2TrueParticleADDPIDTUNE",
                                                {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass(),(double)(trkbin+1)},1.0 );
                                }
			}
/*
			else if( pdgCode == 223 ){
				if( trackMC->GetNDaughters() != 3 ) continue; //omega
				Label1 = trackMC->GetDaughterLabel(0);
				Label2 = trackMC->GetDaughterLabel(0)+1;
				Label3 = trackMC->GetDaughterLabel(1);

//				if( Label1<0 || Label2<0 ) continue;

				AliAODMCParticle* trackd1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( Label1 ));
                                AliAODMCParticle* trackd2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( Label2 ));
				AliAODMCParticle* trackd3 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( Label3 ));

                                if( !trackd1 ) continue;
                                if( !trackd2 ) continue;
				if( !trackd3 ) continue;

                                if( abs( trackd1->PdgCode() ) != 211 ) continue;
                                if( abs( trackd2->PdgCode() ) != 211 ) continue;
				if( abs( trackd3->PdgCode() ) != 111 ) continue;

                                if( trackd1->Pt() < fptcut ) continue;
                                if( fabs( trackd1->Eta() ) > fetacut ) continue;
                                if( trackd1->Pt() < fptcut ) continue;
                                if( fabs( trackd1->Eta() ) > fetacut ) continue;
				if( trackd3->Pt() < fptcut ) continue; 
				if( fabs( trackd3->Eta() ) > fetacut ) continue;

                                trkl1 = trackd1->GetLabel();
                                trkl2 = trackd2->GetLabel();
				trkl3 = trackd3->GetLabel();

                                trk1count=0;
                                trk2count=0;
                                for(UInt_t it=0;it<ntracks;it++){
                                        track = (AliAODTrack*)fEvt->GetTrack(it);
                                        if( !track ) continue;
                                        if( (int)track->GetLabel() == trkl1 && track->GetID() > -1){
                                                trk1count++;
                                                trkid1=it;
                                        }
                                        if( (int)track->GetLabel() == trkl2 && track->GetID() > -1){
                                                trk2count++;
                                                trkid2=it;
                                        }
				}
				if( trk2count == 0 || trk1count == 0 ) continue;

				AliAODTrack* trackd1Recon = (AliAODTrack*)fEvt->GetTrack( trkid1 );
				AliAODTrack* trackd2Recon = (AliAODTrack*)fEvt->GetTrack( trkid2 );

                                if( !trackd1Recon ) continue;
                                if( trackd1Recon->Pt() < fptcut ) continue;
                                if( fabs( trackd1Recon->Eta() ) > fetacut ) continue;

                                if( !trackd2Recon ) continue;
                                if( trackd2Recon->Pt() < fptcut ) continue;
                                if( fabs( trackd2Recon->Eta() ) > fetacut ) continue;

				fPIDResponse->SetTunedOnData(false);
				if( ( fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 2.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 2.0 ) ||
                                    ( fabs( fPIDResponse->NumberOfSigmasTOF((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 3.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTOF((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 3.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 5.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 5.0 ) ){
					FillTHnSparse("hOmgTrueParticleADDPID",
					{genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass()},1.0 );
				}

                                fPIDResponse->SetTunedOnData(true);
                                if( ( fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 2.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 2.0 ) ||
                                    ( fabs( fPIDResponse->NumberOfSigmasTOF((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 3.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTOF((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 3.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 5.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 5.0 ) ){
                                        FillTHnSparse("hOmgTrueParticleADDPIDTUNE",
					{genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass()},1.0 );
                                }
			}
			else if( pdgCode ==313 ){
				if( trackMC->GetNDaughters() != 2 ) continue; //kstart

                                Label1 = trackMC->GetDaughterLabel(0);
                                Label2 = trackMC->GetDaughterLabel(1);

//				if( Label1<0 || Label2<0 ) continue;

                                AliAODMCParticle* trackd1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( Label1 ));
                                AliAODMCParticle* trackd2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At( Label2 ));

                                if( !trackd1 ) continue;
                                if( !trackd2 ) continue;

				if( abs( trackd1->PdgCode() ) != 211 && abs( trackd2->PdgCode() ) != 211 ) continue;
				if( abs( trackd1->PdgCode() ) != 321 && abs( trackd2->PdgCode() ) != 321 ) continue;

                                if( trackd1->Pt() < fptcut ) continue;
                                if( fabs( trackd1->Eta() ) > fetacut ) continue;
                                if( trackd1->Pt() < fptcut ) continue;
                                if( fabs( trackd1->Eta() ) > fetacut ) continue;

                                trkl1 = trackd1->GetLabel();
                                trkl2 = trackd2->GetLabel();

                                trk1count=0;
                                trk2count=0;
                                for(UInt_t it=0;it<ntracks;it++){
                                        track = (AliAODTrack*)fEvt->GetTrack(it);
                                        if( !track ) continue;
                                        if( (int)track->GetLabel() == trkl1 && track->GetID() > -1){
                                                trk1count++;
                                                trkid1=it;
                                        }
                                        if( (int)track->GetLabel() == trkl2 && track->GetID() > -1){
                                                trk2count++;
                                                trkid2=it;
                                        }
                                }
                                if( trk2count == 0 || trk1count == 0 ) continue;

                                AliAODTrack* trackd1Recon = (AliAODTrack*)fEvt->GetTrack( trkid1 );
                                AliAODTrack* trackd2Recon = (AliAODTrack*)fEvt->GetTrack( trkid2 );

                                if( !trackd1Recon ) continue;
                                if( trackd1Recon->Pt() < fptcut ) continue;
                                if( fabs( trackd1Recon->Eta() ) > fetacut ) continue;

                                if( !trackd2Recon ) continue;
                                if( trackd2Recon->Pt() < fptcut ) continue;
                                if( fabs( trackd2Recon->Eta() ) > fetacut ) continue;

                                fPIDResponse->SetTunedOnData(false);
                                if( ( fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 2.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 2.0 ) ||
                                    ( fabs( fPIDResponse->NumberOfSigmasTOF((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 3.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTOF((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 3.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 5.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 5.0 ) ){
                                        FillTHnSparse("hKSTrueParticleADDPID",
                                        {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass()},1.0 );
                                }

                                fPIDResponse->SetTunedOnData(true);
                                if( ( fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 2.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 2.0 ) ||
                                    ( fabs( fPIDResponse->NumberOfSigmasTOF((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 3.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTOF((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 3.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd1Recon, AliPID::kPion) ) < 5.0 &&
                                      fabs( fPIDResponse->NumberOfSigmasTPC((AliVParticle*)trackd2Recon, AliPID::kPion) ) < 5.0 ) ){
                                        FillTHnSparse("hKSTrueParticleADDPIDTUNE",
                                        {genzvtx,fCent,trackMC->Pt(),trackMC->GetCalcMass()},1.0 );
                                }
			}
*/
		}
	}
 }
 if( fRunTable->IsPP() && fOption.Contains("AddMixing") ){
	if (!goodtrackindices.size()) ep->pop_back();
	if ( ep->size() > bookingsize ){
	        for (auto it: ep->front()) delete it;
		ep->pop_front();
	}
 }
 return goodtrackindices.size();

}



void AliAnalysisTaskInclusivef0f2::FillTracks(){
 AliVTrack *track1, *track2;
 int ntracks = goodtrackindices.size();

 tracklist trackpool;


 int MotherID;

 double e1, e2;
 double PiPiMass;
 double PiPipT;
 double Rap_pair;

 TLorentzVector temp1,temp2;
 TLorentzVector vecsum;

 for(UInt_t it=0;it<ntracks-1;it++){
	track1 = (AliVTrack*)fEvt->GetTrack(goodtrackindices[it]);
	if(!track1) continue;
	for(UInt_t jt=it+1;jt<ntracks;jt++){
		track2 = (AliVTrack*)fEvt->GetTrack(goodtrackindices[jt]);
		if(!track2) continue;

		e1 = sqrt( pow(track1->P(),2)+pow(AliPID::ParticleMass(AliPID::kPion),2) );
		e2 = sqrt( pow(track2->P(),2)+pow(AliPID::ParticleMass(AliPID::kPion),2) );

		if( e1+e2-track1->Pz()-track2->Pz() > 0 )
		Rap_pair = 0.5*TMath::Log( (e1+e2+track1->Pz()+track2->Pz())/(e1+e2-track1->Pz()-track2->Pz()) );
		else{ Rap_pair = -999; }
		if( fabs( Rap_pair ) > 0.5 ) continue;
//                if( fRunTable->IsPA() && Rap_pair > 0 ) continue;

		if( fRunTable->IsPA() && ( Rap_pair > 0.035 || Rap_pair < -0.465 ) ) continue;
		PiPiMass = sqrt( pow(e1+e2,2) -
			pow(track1->Px()+track2->Px(),2) -
			pow(track1->Py()+track2->Py(),2) -
			pow(track1->Pz()+track2->Pz(),2) );

		PiPipT = sqrt( pow( track1->Px()+track2->Px(),2 ) +
			       pow( track1->Py()+track2->Py(),2 ) );
/*
		if( IsMC ){
			if( track1->GetLabel() < 0 || track2->GetLabel() < 0 ) continue;
			if( dynamic_cast<AliAODMCParticle*>(fMCArray->At( track1->GetLabel() ))->GetMother() ==
				dynamic_cast<AliAODMCParticle*>(fMCArray->At( track2->GetLabel() ))->GetMother() ){
				MotherID = dynamic_cast<AliAODMCParticle*>(fMCArray->At( track1->GetLabel() ))->GetMother();

				if( dynamic_cast<AliAODMCParticle*>(fMCArray->At( MotherID ))->PdgCode() == 223 ){

					if( track1->Charge()*track2->Charge() == -1 ){
						FillTHnSparse("OmegaRecParticle",{1,fZ,fCent,PiPipT, PiPiMass},1.0 );
					}
					else if( track1->Charge() + track2->Charge() == 2 ){
						FillTHnSparse("OmegaRecParticle",{2,fZ,fCent,PiPipT, PiPiMass},1.0 );
					}
					else if( track1->Charge() + track2->Charge() == -2 ){
						FillTHnSparse("OmegaRecParticle",{3,fZ,fCent,PiPipT, PiPiMass},1.0 );
					}
				}
				else if( dynamic_cast<AliAODMCParticle*>(fMCArray->At( MotherID ))->PdgCode() == 313 ){
                                        	if( track1->Charge()*track2->Charge() == -1 ){
                                        	        FillTHnSparse("KSTARRecParticle",{1,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        	}
                                        	else if( track1->Charge() + track2->Charge() == 2 ){
                                        	        FillTHnSparse("KSTARRecParticle",{2,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        	}
                                        	else if( track1->Charge() + track2->Charge() == -2 ){
                                        	        FillTHnSparse("KSTARRecParticle",{3,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        	}
					
				}
				else if( dynamic_cast<AliAODMCParticle*>(fMCArray->At( MotherID ))->PdgCode() == 333 ){
					if( track1->Charge()*track2->Charge() == -1 ){
                                        	FillTHnSparse("PhiRecParticle",{1,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == 2 ){
                                        	FillTHnSparse("PhiRecParticle",{2,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == -2 ){
                                        	FillTHnSparse("PhiRecParticle",{3,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        }
				}
                                else if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( MotherID ))->PdgCode() ) == 3212 ){
                                        if( track1->Charge()*track2->Charge() == -1 ){
                                        	FillTHnSparse("SigmaRecParticle",{1,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == 2 ){
                                                FillTHnSparse("SigmaRecParticle",{2,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == -2 ){
                                                FillTHnSparse("SigmaRecParticle",{3,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        }
                                }
                                else if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( MotherID ))->PdgCode() ) == 3122 ){
                                        if( track1->Charge()*track2->Charge() == -1 ){
                                        	FillTHnSparse("LambdaRecParticle",{1,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == 2 ){
                                                FillTHnSparse("LambdaRecParticle",{2,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == -2 ){
                                                FillTHnSparse("LambdaRecParticle",{3,fZ,fCent,PiPipT, PiPiMass},1.0 );
                                        }
                                }
				else if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( MotherID ))->PdgCode() ) == 20313 ){
					if( track1->Charge()*track2->Charge() == -1 ){
						FillTHnSparse("ExKaonRecParticle",{1,fZ,fCent,PiPipT, PiPiMass,1},1.0 );
					}
					else if( track1->Charge() + track2->Charge() == 2 ){
						FillTHnSparse("ExKaonRecParticle",{2,fZ,fCent,PiPipT, PiPiMass,1},1.0 );
					}
					else if( track1->Charge() + track2->Charge() == -2 ){
						FillTHnSparse("ExKaonRecParticle",{3,fZ,fCent,PiPipT, PiPiMass,1},1.0 );
					}	
				}
				else if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( MotherID ))->PdgCode() ) == 100313 ){
                                        if( track1->Charge()*track2->Charge() == -1 ){
                                                FillTHnSparse("ExKaonRecParticle",{1,fZ,fCent,PiPipT, PiPiMass,2},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == 2 ){
                                                FillTHnSparse("ExKaonRecParticle",{2,fZ,fCent,PiPipT, PiPiMass,2},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == -2 ){
                                                FillTHnSparse("ExKaonRecParticle",{3,fZ,fCent,PiPipT, PiPiMass,2},1.0 );
                                        }
				}
				else if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( MotherID ))->PdgCode() ) == 10311 ){
                                        if( track1->Charge()*track2->Charge() == -1 ){
                                                FillTHnSparse("ExKaonRecParticle",{1,fZ,fCent,PiPipT, PiPiMass,3},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == 2 ){
                                                FillTHnSparse("ExKaonRecParticle",{2,fZ,fCent,PiPipT, PiPiMass,3},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == -2 ){
                                                FillTHnSparse("ExKaonRecParticle",{3,fZ,fCent,PiPipT, PiPiMass,3},1.0 );
                                        }
				}
				else if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( MotherID ))->PdgCode() ) == 315 ){
                                        if( track1->Charge()*track2->Charge() == -1 ){
                                                FillTHnSparse("ExKaonRecParticle",{1,fZ,fCent,PiPipT, PiPiMass,4},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == 2 ){
                                                FillTHnSparse("ExKaonRecParticle",{2,fZ,fCent,PiPipT, PiPiMass,4},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == -2 ){
                                                FillTHnSparse("ExKaonRecParticle",{3,fZ,fCent,PiPipT, PiPiMass,4},1.0 );
                                        }
				}
				else if( abs( dynamic_cast<AliAODMCParticle*>(fMCArray->At( MotherID ))->PdgCode() ) == 30313 ){
                                        if( track1->Charge()*track2->Charge() == -1 ){
                                                FillTHnSparse("ExKaonRecParticle",{1,fZ,fCent,PiPipT, PiPiMass,5},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == 2 ){
                                                FillTHnSparse("ExKaonRecParticle",{2,fZ,fCent,PiPipT, PiPiMass,5},1.0 );
                                        }
                                        else if( track1->Charge() + track2->Charge() == -2 ){
                                                FillTHnSparse("ExKaonRecParticle",{3,fZ,fCent,PiPipT, PiPiMass,5},1.0 );
                                        }
				}
			}
			else if( dynamic_cast<AliAODMCParticle*>(fMCArray->At( track1->GetLabel() ))->GetMother() != 
				dynamic_cast<AliAODMCParticle*>(fMCArray->At( track2->GetLabel() ))->GetMother() ){
				if( track1->Charge()*track2->Charge() == -1 ){
					FillTHnSparse("hInvMassUnpair",{1,fZ,fCent,PiPipT, PiPiMass},1.0 );
				}
				else if( track1->Charge() + track2->Charge() == 2 ){
                                        FillTHnSparse("hInvMassUnpair",{2,fZ,fCent,PiPipT, PiPiMass},1.0 );
				}
				else if( track1->Charge() + track2->Charge() == -2 ){
                                        FillTHnSparse("hInvMassUnpair",{3,fZ,fCent,PiPipT, PiPiMass},1.0 );
				}
			}
		}
*/

		if( track1->Charge()*track2->Charge() == -1 ){
			FillTHnSparse("hInvMass",{1,fZ,fCent,
				PiPipT, PiPiMass,(double)(trkbin+1)},1.0 );
		}
		else if( track1->Charge() + track2->Charge() == 2 ){
			FillTHnSparse("hInvMass",{2,fZ,fCent,
				PiPipT, PiPiMass,(double)(trkbin+1)},1.0 );
		}
		else if( track1->Charge() + track2->Charge() == -2 ){
			FillTHnSparse("hInvMass",{3,fZ,fCent,
				PiPipT, PiPiMass,(double)(trkbin+1)},1.0 );
		}
	}
 }

 int epsize=1;
 if( centbin>=0 && zbin>=0 && fRunTable->IsPP() && fOption.Contains("AddMixing") ){
        eventpool &ep = fEMpooltrk[trkbin][centbin][zbin];
        epsize = ep.size();
        if (ep.size()< 10 ) return;
        int n = 0;
        for (auto pool: ep){
                if (n == (ep.size() -1 )) continue;
                for (auto track: pool) trackpool.push_back((AliVTrack*)track);
                n++;
        }
 }


 if( fRunTable->IsPP() && fOption.Contains("AddMixing") ){
 	for (UInt_t  it = 0; it < ntracks; it++) {
		track1 =  (AliVTrack*) fEvt->GetTrack(goodtrackindices.at(it)) ;
		if(!track1) continue;
	        for (UInt_t jt = 0; jt < trackpool.size(); jt++) {
	        	track2 = trackpool.at(jt);
			if(!track2) continue;
 
	                e1 = sqrt( pow(track1->P(),2)+pow(AliPID::ParticleMass(AliPID::kPion),2) );
	                e2 = sqrt( pow(track2->P(),2)+pow(AliPID::ParticleMass(AliPID::kPion),2) );

	                PiPiMass = sqrt( pow(e1+e2,2) -
	                        pow(track1->Px()+track2->Px(),2) -
	                        pow(track1->Py()+track2->Py(),2) -
	                        pow(track1->Pz()+track2->Pz(),2) );
	
			PiPipT = sqrt( pow( track1->Px()+track2->Px(),2 ) +
	                               pow( track1->Py()+track2->Py(),2 ) );

	                if( track1->Charge()*track2->Charge() == -1 ){
	                        FillTHnSparse("hInvMassMixing",{1,fZ,fCent,
	                                PiPipT, PiPiMass,(double)(trkbin+1)},1.0);
	                }
	                else if( track1->Charge() + track2->Charge() == 2 ){
	                        FillTHnSparse("hInvMassMixing",{2,fZ,fCent,
	                                PiPipT, PiPiMass,(double)(trkbin+1)},1.0);
	                }
	                else if( track1->Charge() + track2->Charge() == -2 ){
	                        FillTHnSparse("hInvMassMixing",{3,fZ,fCent,
	                                PiPipT, PiPiMass,(double)(trkbin+1)},1.0);
	                }
	        }
	}
 }

}
void AliAnalysisTaskInclusivef0f2::FinishTaskOutput()
{
 cout << "FinishTaskOutput" << endl;
}


void AliAnalysisTaskInclusivef0f2::Terminate(Option_t *)
{
 cout << "Terminate" << endl;
}


int AliAnalysisTaskInclusivef0f2::GetPID(AliPIDResponse *pid, const AliVTrack *trk){
 if (!pid) return -1; 

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


THnSparse * AliAnalysisTaskInclusivef0f2::CreateTHnSparse(TString name
        , TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt){
  const TAxis * axises[bins.size()];
  for( UInt_t i=0;i<bins.size();i++ ) axises[i]= &bins[i];
  THnSparse * h= fHistos->CreateTHnSparse(name, title, ndim, axises,opt );
  return h;
}

THnSparse * AliAnalysisTaskInclusivef0f2::CreateTHnSparse(TString name
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

Long64_t AliAnalysisTaskInclusivef0f2::FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w ){
  auto hsparse = dynamic_cast<THnSparse*>( fHistos->FindObject(name) );
  if(! hsparse ){
    cout<<"ERROR : no "<<name<<endl;
    exit(1);
  }
  return FillTHnSparse( hsparse, x, w );
}

Long64_t AliAnalysisTaskInclusivef0f2::FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w ){
  if( int(x.size()) != h->GetNdimensions() ){
    cout<<"ERROR : wrong sized of array while Fill "<<h->GetName()<<endl;
    exit(1);
  }
  return h->Fill( &x.front(), w );
}

TAxis AliAnalysisTaskInclusivef0f2::AxisFix
        ( TString name, int nbin, Double_t xmin, Double_t xmax ){
                TAxis axis(nbin, xmin, xmax);axis.SetName(name);
                return axis;
}

TAxis AliAnalysisTaskInclusivef0f2::AxisStr( TString name, std::vector<TString> bin ){
  TAxis ax = AxisFix( name, bin.size(), 0.5, bin.size()+0.5);
  UInt_t i=1;
  for( auto blabel : bin )
    ax.SetBinLabel( i++, blabel );
  return ax;
}

TAxis AliAnalysisTaskInclusivef0f2::AxisVar( TString name, std::vector<Double_t> bin ){
  TAxis axis( bin.size()-1, &bin.front() ) ;axis.SetName(name);
  return axis;
}

TAxis AliAnalysisTaskInclusivef0f2::AxisLog
        ( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0){
  int binoffset = ( xmin0<0 || (xmin-xmin0)<1e-9) ? 0 : 1;
  std::vector<Double_t> bin(nbin+1+binoffset,0);
  double logBW3 = (log(xmax)-log(xmin))/nbin;
  for(int ij=0;ij<=nbin;ij++) bin[ij+binoffset]=xmin*exp(ij*logBW3);
  TAxis axis( nbin, &bin.front() ) ;
  axis.SetName(name);
  return axis;
}

