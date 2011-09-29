#include "AliV0ReaderV1.h"
#include "AliKFParticle.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliKFParticle.h"
#include "AliKFConversionPhoton.h"
#include "AliAODConversionPhoton.h"
#include "AliConversionPhotonBase.h"
#include "TVector.h"
#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "AliStack.h"

class iostream;


using namespace std;

ClassImp(AliV0ReaderV1)

//________________________________________________________________________
AliV0ReaderV1::AliV0ReaderV1(const char *name) : AliAnalysisTaskSE(name),
fConversionGammas(NULL),
fESDEvent(NULL),
fAODEvent(NULL),
fMCStack(NULL),
fOutputList(NULL),
fCurrentMotherKFCandidate(NULL),
fCurrentPositiveKFParticle(NULL),
fCurrentNegativeKFParticle(NULL),
fCurrentTrackLabels(NULL),
fCurrentV0Index(-1),
fNCentralityBins(5),
    fCentralityBin(0),
    fBGHandler(NULL),
    fVertexZ(-999),
    fCentrality(-1),
    fEPAngle(-1),
fMaxVertexZ(10),
fMaxR(180),// 100 meter(outside of ALICE)
fMinR(5),// 100 meter(outside of ALICE)
fEtaCut(0.9),
fEtaCutMin(-0.1),
fPtCut(0.),
fSinglePtCut(0.),
fMaxZ(240),
fMinClsTPC(0.),
fMinClsTPCToF(0.),
fLineCutZRSlope(0.),
fLineCutZValue(7),
fLineCutZRSlopeMin(0.),
fLineCutZValueMin(-2),
fChi2CutConversion(30),
fPIDProbabilityCutNegativeParticle(0),
fPIDProbabilityCutPositiveParticle(0),
fDodEdxSigmaCut(kTRUE),
fDoTOFsigmaCut(kFALSE), // RRnewTOF
fPIDTRDEfficiency(0.95),
fDoTRDPID(kFALSE),
fPIDnSigmaAboveElectronLine(5),
fPIDnSigmaBelowElectronLine(-3),
fTofPIDnSigmaAboveElectronLine(100), // RRnewTOF
fTofPIDnSigmaBelowElectronLine(-100), // RRnewTOF
fPIDnSigmaAbovePionLine(0),
fPIDnSigmaAbovePionLineHighPt(-100),
fPIDMinPnSigmaAbovePionLine(1),
fPIDMaxPnSigmaAbovePionLine(3),
fDoKaonRejectionLowP(kTRUE),
    fDoProtonRejectionLowP(kTRUE),
    fDoPionRejectionLowP(kTRUE),
    fPIDnSigmaAtLowPAroundKaonLine(0),
    fPIDnSigmaAtLowPAroundProtonLine(0),
    fPIDnSigmaAtLowPAroundPionLine(0),
    fPIDMinPKaonRejectionLowP(1.5),
    fPIDMinPProtonRejectionLowP(2),
    fPIDMinPPionRejectionLowP(0.5),
    fDoQtGammaSelection(kTRUE),
    fDoHighPtQtGammaSelection(kTRUE), // RRnew
    fQtMax(0.05),
    fHighPtQtMax(0.06), // RRnew
    fPtBorderForQt(2.5), // RRnew
	fNSigmaMass(0.),
	fUseImprovedVertex(kTRUE),
	fUseOwnXYZCalculation(kTRUE),
	fUseConstructGamma(kFALSE),
	fUseEtaMinCut(kFALSE),
	fUseOnFlyV0Finder(kTRUE),
	fDoPhotonAsymmetryCut(kTRUE),
	fMinPPhotonAsymmetryCut(100.),
fMinPhotonAsymmetry(0.),
    kUseAODConversionPhoton(kFALSE),
    fIsHeavyIon(kTRUE),
    fCreateAOD(kFALSE),
    fDeltaAODFilename("AliAODGammaConversion.root")
{

    fLineCutZRSlope=tan(2*atan(exp(-fEtaCut)));
        // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());

  fCurrentTrackLabels=new Int_t[2];
}


//________________________________________________________________________
AliV0ReaderV1::~AliV0ReaderV1()
{
    if(fConversionGammas){
	fConversionGammas->Delete();// Clear Objects
        delete fConversionGammas;
	fConversionGammas=0x0;
    }

    if(fCurrentTrackLabels){
	delete[] fCurrentTrackLabels;
	fCurrentTrackLabels=NULL;}

}

//________________________________________________________________________
void AliV0ReaderV1::UserCreateOutputObjects()
{

if(fOutputList != NULL){
    delete fOutputList;
    fOutputList = NULL;
  }
  if(fOutputList == NULL){
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
  }

  TList *fCutList=new TList();
  fCutList->SetName("GammaReconstruction");
  fCutList->SetOwner(kTRUE);
  fOutputList->Add(fCutList);

//GammaMass-plots
Int_t kGCnXBinsGammaMass 		= 4000;
Double_t kGCfirstXBinGammaMass= 0.;
Double_t kGClastXBinGammaMass = 1.;

  Int_t kGCnYBinsSpectra = 250;
  Double_t kGCfirstYBinSpectra = 0.;
  Double_t kGClastYBinSpectra = 25.;

// Process Gammas Histograms

hV0CurrentFinder=new TH1F("ESD_V0sCurrentFinder_InvMass","V0sCurrentFinder",kGCnXBinsGammaMass,kGCfirstXBinGammaMass,kGClastXBinGammaMass);
fCutList->Add(hV0CurrentFinder);
hV0AllArmenteros=new TH2F("ESD_V0sCurrentFinder_Armenteros","Armenteros Alpha Qt",200,-1,1,250,0,0.25);
fCutList->Add(hV0AllArmenteros);
hV0Good=new TH1F("ESD_GoodV0s_InvMass","GoodV0s",kGCnXBinsGammaMass,kGCfirstXBinGammaMass,kGClastXBinGammaMass);
fCutList->Add(hV0Good);
hV0GoodArmenteros=new TH2F("ESD_GoodV0s_Armenteros","Armenteros Alpha Qt",200,-1,1,250,0,0.25);
fCutList->Add(hV0GoodArmenteros);


// Track Cuts
hV0CutLikeSign=new TH1F("ESD_CutLikeSign_InvMass","LikeSign",kGCnXBinsGammaMass,kGCfirstXBinGammaMass,kGClastXBinGammaMass);
fCutList->Add(hV0CutLikeSign);
hV0CutRefit=new TH1F("ESD_CutRefit_InvMass","No TPC refit",kGCnXBinsGammaMass,kGCfirstXBinGammaMass,kGClastXBinGammaMass);
fCutList->Add(hV0CutRefit);
hV0CutKinks=new TH1F("ESD_CutKink_InvMass","Kinks",kGCnXBinsGammaMass,kGCfirstXBinGammaMass,kGClastXBinGammaMass);
fCutList->Add(hV0CutKinks);
hV0CutMinNclsTPCToF=new TH1F("ESD_CutMinNClsTPCToF_InvMass","Min Ncls TPC ToF",kGCnXBinsGammaMass,kGCfirstXBinGammaMass,kGClastXBinGammaMass);
fCutList->Add(hV0CutMinNclsTPCToF);

// Event Cuts
/*hV0CutNContributors=new TH1F("ESD_CutNContributors_InvMass","NContributors<=0",kGCnXBinsGammaMass,kGCfirstXBinGammaMass,kGClastXBinGammaMass);
fCutList->Add(hV0CutNContributors);
hV0CutVertexZ=new TH1F("ESD_CutVertexZ_InvMass","VertexZ",kGCnXBinsGammaMass,kGCfirstXBinGammaMass,kGClastXBinGammaMass);
fCutList->Add(hV0CutVertexZ);
*/


// dEdx Cuts

hV0CutdEdxElectron=new TH1F("ESD_CutdEdxSigmaElectronLine_InvMass" ,"dedx ElectronLine" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutdEdxElectron);
hV0CutdEdxPion=new TH1F("ESD_CutdEdxSigmaPionLine_InvMass" ,"dedx PionLine" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutdEdxPion);
hV0CutdEdxKaonLowP=new TH1F("ESD_CutKaonRejectionLowP_InvMass" ,"dedx KaonRejection LowP" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutdEdxKaonLowP);
hV0CutdEdxProtonLowP=new TH1F("ESD_CutProtonRejectionLowP_InvMass" ,"dedx ProtonRejection LowP" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutdEdxProtonLowP);
hV0CutdEdxPionLowP=new TH1F("ESD_CutPionRejectionLowP_InvMass" ,"dedx PionRejection LowP" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutdEdxPionLowP);
hV0CutdEdxTOFElectron=new TH1F("ESD_CutTOFsigmaElec_InvMass", "ESD_CutTOFsigmaElec_InvMass",kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutdEdxTOFElectron);
hV0CutdEdxTRD=new TH1F("ESD_CutTRD_InvMass", "ESD_CutTRD_InvMass",kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutdEdxTRD);
hGammadEdxbefore=new TH2F("Gamma_dEdx_before","dEdx Gamma before" ,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,400, 0,200);
fCutList->Add(hGammadEdxbefore);
hGammadEdxafter=new TH2F("Gamma_dEdx_after","dEdx Gamma after" ,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,400, 0,200);
fCutList->Add(hGammadEdxafter);


// Armenteros

hV0CutQt=new TH1F("ESD_CutQt_InvMass","ESD_CutQt_InvMass",kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutQt);

// Kinematic Cuts

hV0CutR=new TH1F("ESD_CutR_InvMass" ,"Above RMax" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutR);
hV0CutMinR=new TH1F("ESD_CutMinR_InvMass" ,"Above RMax" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutMinR);
hV0CutLine=new TH1F("ESD_CutLine_InvMass" ,"Out of reconstruction area" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutLine);
hV0CutZ=new TH1F("ESD_CutZ_InvMass" ,"Out of reconstruction area" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutZ);
hV0CutEta=new TH1F("ESD_CutEta_InvMass" ,"Above #eta max" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutEta);


hV0CutSinglePt=new TH1F("ESD_CutSinglePt_InvMass" ,"Below p_{t} min" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutSinglePt);
hV0CutNDF=new TH1F("ESD_CutNDF_InvMass" ,"#chi^{2} > Max" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutNDF);
hV0CutChi2=new TH1F("ESD_CutChi2_InvMass" ,"#chi^{2} > Max" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutChi2);
hV0CutPt= new TH1F("ESD_CutPt_InvMass" ,"Below p_{t} min" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutPt);

// Asymmetry Cut

hV0CutAsymmetry=new TH1F("ESD_CutPhotonAsymmetry_InvMass" ,"Out of reconstruction area" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutAsymmetry);

// PID Prob

hV0CutPIDProb=new TH1F("ESD_CutPIDProb_InvMass" ,"wrong TPC PID" , kGCnXBinsGammaMass, kGCfirstXBinGammaMass, kGClastXBinGammaMass);
fCutList->Add(hV0CutPIDProb);

// Event Info

   // Other
  TList *fOtherList=new TList();
  fOtherList->SetName("EventInfo");
  fOtherList->SetOwner(kTRUE);
  fOutputList->Add(fOtherList);


  hV0EventCuts=new TH1F("ESD_EventCuts","Event Cuts",10,-0.5,9.5);
  fOtherList->Add(hV0EventCuts);
  hNEvents=new TH1F("NEvents_vs_Centrality","NEvents vs Centrality",fNCentralityBins,-0.5,fNCentralityBins-0.5);
  fOtherList->Add(hNEvents);
  hCentrality=new TH1F("Centrality","Centrality",100,0,100);
  fOtherList->Add(hCentrality);
  hVertexZ=new TH1F("VertexZ","VertexZ",1000,-50,50);
  fOtherList->Add(hVertexZ);

// QA

  TList *fGammaList=new TList();
  fGammaList->SetName("GammaInfo");
  fGammaList->SetOwner(kTRUE);
  fOutputList->Add(fGammaList);

  hGammaPt=new TH1F*[fNCentralityBins];

  for(int i=0;i<fNCentralityBins;i++){
      hGammaPt[i]=new TH1F(Form("GammaSpectrum_RECO_Pt_%d",i),"Reco Gamma Pt",kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
      fGammaList->Add(hGammaPt[i]);
  }

  hGammaPhi=new TH1F("Phi_Gamma","Phi Gamma" ,36, 0, 2*TMath::Pi());
  fGammaList->Add(hGammaPhi);
  hGammaConversionMapXY=new TH2F("Gamma_ConversionMap_XY","Conversion Point xy",400,-200,200,400,-200,200);
  fGammaList->Add(hGammaConversionMapXY);
  hGammaConversionMapZR=new TH2F("Gamma_ConversionMap_ZR","Conversion Point zr",500,-250,250,180,0,180);
  fGammaList->Add(hGammaConversionMapZR);


  // FILL MC PART only if MC is available
  if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){

    hMCPtTRUE=new TH1F*[fNCentralityBins];
    hMCPtRECOTRUE=new TH1F*[fNCentralityBins];

    for(int i=0;i<fNCentralityBins;i++){
        hMCPtTRUE[i]=new TH1F(Form("GammaSpectrum_MC_Pt_%d",i),"TRUE Gamma Pt",kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
	fGammaList->Add(hMCPtTRUE[i]);
	hMCPtRECOTRUE[i]=new TH1F(Form("GammaSpectrum_RECOTRUE_Pt_%d",i),"True Reco Gamma Pt",kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
	fGammaList->Add(hMCPtRECOTRUE[i]);
    }

    // Call Sumw2 Option
    for (Int_t i=0; i<fGammaList->GetEntries(); i++) {
	TH1 *h1 = dynamic_cast<TH1*>(fGammaList->At(i));
	if (h1){h1->Sumw2();}
    }

    hMCPtResolution=new TH2F("Resolution_Gamma_dPt_Pt","dPt vs Pt", kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra,200,-10,10);
    fGammaList->Add(hMCPtResolution);
    hMCPtResolutionPhi=new TH2F("Resolution_Gamma_dPt_Phi","dPt vs Phi",180,0,2*TMath::Pi(),200,-10,10);
    fGammaList->Add(hMCPtResolutionPhi);
    hMCRResolutionvsR=new TH2F("Resolution_dRAbs_VS_R","dR vs R", 720,0,360,100,-5,5);
    fGammaList->Add(hMCRResolutionvsR);
    hMCZResolutionvsZ=new TH2F("Resolution_dZAbs_VS_Z","dZ vs Z", 200,-50,50,100,-5,5);
    fGammaList->Add(hMCZResolutionvsZ);

}
// Gamma Output

if(fCreateAOD){kUseAODConversionPhoton=kTRUE;}

if(fConversionGammas == NULL){
    if(kUseAODConversionPhoton){
	fConversionGammas = new TClonesArray("AliAODConversionPhoton",100);}
    else{
	fConversionGammas = new TClonesArray("AliKFConversionPhoton",100);}
}
fConversionGammas->Delete();//Reset the TClonesArray

// Create AODs

if(fCreateAOD){
    fConversionGammas->SetName(Form("GammaConv_gamma"));

    AddAODBranch("TClonesArray", &fConversionGammas, fDeltaAODFilename.Data());
    AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile(fDeltaAODFilename.Data());
}

  PostData(1, fOutputList);

}
//________________________________________________________________________
void AliV0ReaderV1::UserExec(Option_t *){

    fAODEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
    fESDEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
    if(!fAODEvent&&!fESDEvent) {
	AliError("No Input event");
	return;
    }

    fMCStack=NULL;
    if(fMCEvent){
	fMCStack = fMCEvent->Stack();}


    if(fESDEvent)AliKFParticle::SetField(fESDEvent->GetMagneticField());
    if(fAODEvent)AliKFParticle::SetField(fAODEvent->GetMagneticField());

    fConversionGammas->Delete();//Reset the TClonesArray

    // Event Cuts

    EventCuts();

    if(EventIsSelected()){

	// Process V0s
	for(fCurrentV0Index=0;fCurrentV0Index<fInputEvent->GetNumberOfV0s();fCurrentV0Index++){
	    if(CheckV0Status()){

		ProcessV0();
	    }
	}
	ProcessMCGammasForEfficiency();

	// Set AOD Output

	///Make sure delta aod is filled if standard aod is filled (for synchronization when reading aod with standard aod)
	if(fCreateAOD) {
	    AliAODHandler * aodhandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
	    if (aodhandler && aodhandler->GetFillAOD()) {
		AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillExtension(kTRUE);
	    }}
    }

  PostData(1, fOutputList);
}
 ///________________________________________________________________________
Bool_t AliV0ReaderV1::IsGammaCandidate(AliConversionPhotonBase *fPhotonCandidate)
{

    // Fill Histos before Cuts
    hV0CurrentFinder->Fill(fPhotonCandidate->GetPhotonMass());
    hV0AllArmenteros->Fill(fPhotonCandidate->GetArmenterosAlpha(),fPhotonCandidate->GetArmenterosQt());

    Bool_t passcuts=kTRUE;

    // Gamma selection based on QT from Armenteros
    if(fDoQtGammaSelection == kTRUE){
        if(!ArmenterosQtCut(fPhotonCandidate))return kFALSE;//passcuts=kFALSE;
    }

    // Chi Cut

    if(fPhotonCandidate->GetChi2perNDF() > fChi2CutConversion || fPhotonCandidate->GetChi2perNDF() <=0){
	hV0CutChi2->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
    }

    // Reconstruction Acceptance Cuts
    if(!AcceptanceCuts(fPhotonCandidate))return kFALSE;//passcuts=kFALSE;


    // Track Cuts
    if(!TrackCuts(fPhotonCandidate))return kFALSE;//passcuts=kFALSE;

   
    // PID Cuts
    if(!dEdxCuts(fPhotonCandidate))return kFALSE;//passcuts=kFALSE;

    // Asymmetry Cut
    if(fDoPhotonAsymmetryCut == kTRUE){
        if(!AsymmetryCut(fPhotonCandidate))return kFALSE;//passcuts=kFALSE;
    }

    //Check the pid probability

   if(!PIDProbabilityCut(fPhotonCandidate))return kFALSE;//passcuts=kFALSE;

return passcuts;

}

///________________________________________________________________________
const AliExternalTrackParam *AliV0ReaderV1::GetExternalTrackParam(Int_t charge){


    if(!(charge==1||charge==-1)){AliError("Charge not defined");return 0x0;}

    Int_t label;
    if(charge>0)label=0;
    else label=1;
    // Check for sign flip

    if(fESDEvent){
	AliESDv0 *fCurrentV0=dynamic_cast<AliESDv0*>(fESDEvent->GetV0(fCurrentV0Index));
	if(fCurrentV0){
	    if(!fCurrentV0->GetParamN()||!fCurrentV0->GetParamP())return 0x0;
	    if(!fESDEvent->GetTrack(fCurrentV0->GetNindex())||!fESDEvent->GetTrack(fCurrentV0->GetPindex()))return 0x0;
	    if((fESDEvent->GetTrack(fCurrentV0->GetPindex()))->Charge()==charge){
		fCurrentTrackLabels[label]=fCurrentV0->GetPindex();
		return fCurrentV0->GetParamP();}
	    if((fESDEvent->GetTrack(fCurrentV0->GetNindex()))->Charge()==charge){
                fCurrentTrackLabels[label]=fCurrentV0->GetNindex();
		return fCurrentV0->GetParamN();}
	}

    }

    if(fAODEvent){
/*	AliAODv0 *fCurrentV0=dynamic_cast<AliAODv0*>(fAODEvent->GetV0(fCurrentV0Index));
	if(fCurrentV0){
	    if(!fCurrentV0->GetParamN()||!fCurrentV0->GetParamP())return 0x0;
	    if(!fAODEvent->GetTrack(fCurrentV0->GetNegID())||!fAODEvent->GetTrack(fCurrentV0->GetPosID()))return 0x0;
	    if((fAODEvent->GetTrack(fCurrentV0->GetPosID()))->Charge()==1){
                fCurrentTrackLabels[label]=fCurrentV0->GetPosID();
		return fCurrentV0->GetParamP();}
	    if((fAODEvent->GetTrack(fCurrentV0->GetNegID()))->Charge()==1){
		fCurrentTrackLabels[label]=fCurrentV0->GetNegID();
		return fCurrentV0->GetParamN();}
        }
  */  }
    return 0x0;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::CheckV0Status()
{
    if(fESDEvent){
	AliESDv0 *fCurrentV0=(AliESDv0*)(fESDEvent->GetV0(fCurrentV0Index));
	if(!fCurrentV0){
	    printf("Requested V0 does not exist");
	    return kFALSE;}
	//checks if on the fly mode is set
	if(fCurrentV0->GetOnFlyStatus()==fUseOnFlyV0Finder)return kTRUE;
    }

    if(fAODEvent){
	     AliAODv0 *fCurrentV0=dynamic_cast<AliAODv0*>(fAODEvent->GetV0(fCurrentV0Index));
	     if(!fCurrentV0){
		 AliWarning("Requested V0 does not exist");
		 return kFALSE;}

	     //checks if on the fly mode is set
	     if(fCurrentV0->GetOnFlyStatus()==fUseOnFlyV0Finder)return kTRUE;
    }
    return kFALSE;
}




///________________________________________________________________________
void AliV0ReaderV1::ProcessV0(){

    // Reset TrackLabels
    fCurrentTrackLabels[0]=-1;
    fCurrentTrackLabels[1]=-1;


   //  cout<<"V0ReaderV1 ProcessV0 "<<fCurrentV0Index<<endl;

    // Get Daughter KF Particles

    const AliExternalTrackParam *fCurrentExternalTrackParamPositive=GetExternalTrackParamP();
    const AliExternalTrackParam *fCurrentExternalTrackParamNegative=GetExternalTrackParamN();

    if(fCurrentExternalTrackParamPositive&&fCurrentExternalTrackParamNegative){

	fCurrentNegativeKFParticle=new AliKFParticle(*(fCurrentExternalTrackParamNegative),11);
	fCurrentPositiveKFParticle=new AliKFParticle(*(fCurrentExternalTrackParamPositive),-11);

    }
    //else{hV0CutLikeSign->Fill(fCurrentMotherKFCandidate->M());}// Like Sign error is already catched here


    // Reconstruct Gamma

    if(fCurrentNegativeKFParticle&&fCurrentPositiveKFParticle){

	if(fUseConstructGamma==kTRUE){
         
	    fCurrentMotherKFCandidate = new AliKFConversionPhoton();
		fCurrentMotherKFCandidate->ConstructGamma(*fCurrentNegativeKFParticle,*fCurrentPositiveKFParticle);
	}else{
		fCurrentMotherKFCandidate = new AliKFConversionPhoton(*fCurrentNegativeKFParticle,*fCurrentPositiveKFParticle);
		fCurrentMotherKFCandidate->SetMassConstraint(0,fNSigmaMass);
	}

	if(fCurrentNegativeKFParticle){delete fCurrentNegativeKFParticle;
	    fCurrentNegativeKFParticle=0x0;}
	if(fCurrentPositiveKFParticle){ delete fCurrentPositiveKFParticle;
	    fCurrentPositiveKFParticle=0x0;}


        // Update Vertex
        if(fUseImprovedVertex == kTRUE){
		AliKFVertex primaryVertexImproved(*GetPrimaryVertex());
		primaryVertexImproved+=*fCurrentMotherKFCandidate;
		fCurrentMotherKFCandidate->SetProductionVertex(primaryVertexImproved);
	}

	// Set Track Labels

	fCurrentMotherKFCandidate->SetV0Index(fCurrentV0Index);
	fCurrentMotherKFCandidate->SetTrackLabels(fCurrentTrackLabels[0],fCurrentTrackLabels[1]);

	//Set MC Label

	if(fMCStack){
            Int_t labeln=TMath::Abs(GetTrack(fCurrentMotherKFCandidate->GetTrackLabelPositive())->GetLabel());
	    Int_t labelp=TMath::Abs(GetTrack(fCurrentMotherKFCandidate->GetTrackLabelNegative())->GetLabel());

	       TParticle *fNegativeMCParticle = fMCStack->Particle(labeln);
	       TParticle *fPositiveMCParticle = fMCStack->Particle(labelp);

	       if(fPositiveMCParticle&&fNegativeMCParticle){
		   fCurrentMotherKFCandidate->SetMCLabelPositive(labelp);
                   fCurrentMotherKFCandidate->SetMCLabelNegative(labeln);
	       }
	}



	//Add PID information with ESD tender (AOD implementation is not complete)


	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
	AliPIDResponse *fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();

	if(fESDEvent){
	    Int_t labelp=((AliESDv0*)fESDEvent->GetV0(fCurrentMotherKFCandidate->GetV0Index()))->GetNindex();
	    Int_t labeln=((AliESDv0*)fESDEvent->GetV0(fCurrentMotherKFCandidate->GetV0Index()))->GetNindex();

	AliESDtrack *trackpos=fESDEvent->GetTrack(labelp);
	AliESDtrack *trackneg=fESDEvent->GetTrack(labeln);

	if(trackpos&&trackneg){

	    Float_t fNSigmadEdxPositive[5];
	    Float_t fNSigmadEdxNegative[5];

	    fNSigmadEdxPositive[0]=fPIDResponse->NumberOfSigmasTPC(trackpos,AliPID::kElectron);
	    fNSigmadEdxPositive[1]=fPIDResponse->NumberOfSigmasTPC(trackpos,AliPID::kMuon);
	    fNSigmadEdxPositive[2]=fPIDResponse->NumberOfSigmasTPC(trackpos,AliPID::kPion);
	    fNSigmadEdxPositive[3]=fPIDResponse->NumberOfSigmasTPC(trackpos,AliPID::kKaon);
	    fNSigmadEdxPositive[4]=fPIDResponse->NumberOfSigmasTPC(trackpos,AliPID::kProton);

	    fNSigmadEdxNegative[0]=fPIDResponse->NumberOfSigmasTPC(trackneg,AliPID::kElectron);
	    fNSigmadEdxNegative[1]=fPIDResponse->NumberOfSigmasTPC(trackneg,AliPID::kMuon);
	    fNSigmadEdxNegative[2]=fPIDResponse->NumberOfSigmasTPC(trackneg,AliPID::kPion);
	    fNSigmadEdxNegative[3]=fPIDResponse->NumberOfSigmasTPC(trackneg,AliPID::kKaon);
	    fNSigmadEdxNegative[4]=fPIDResponse->NumberOfSigmasTPC(trackneg,AliPID::kProton);

	fCurrentMotherKFCandidate->SetNSigmadEdx(fNSigmadEdxPositive,fNSigmadEdxNegative);
	}
    }

	// Calculate ConversionPoint

	if(fUseOwnXYZCalculation){
        
	    Double_t convpos[3]={0,0,0};
	    GetConversionPoint(fCurrentExternalTrackParamPositive,fCurrentExternalTrackParamNegative,convpos);
	    fCurrentMotherKFCandidate->SetConversionPoint(convpos);

	}
	
       // if(kTRUE){
       if(IsGammaCandidate(fCurrentMotherKFCandidate)){
	   // Fill Histos after Cuts

	   // Process MC

	   ProcessMC(fCurrentMotherKFCandidate);


	    hV0Good->Fill(fCurrentMotherKFCandidate->M());
	    hV0GoodArmenteros->Fill(fCurrentMotherKFCandidate->GetArmenterosAlpha(),fCurrentMotherKFCandidate->GetArmenterosQt());


	    // Set Mass Zero for Gammas
	    SetGammaMassZero();

	    // Add Gamma to the TClonesArray

	    if(kUseAODConversionPhoton){
		new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliAODConversionPhoton(fCurrentMotherKFCandidate);
	    }
	    else{
		new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliKFConversionPhoton(*fCurrentMotherKFCandidate);
	    }
	    // Fill QA Histos

	    hGammaPhi->Fill(fCurrentMotherKFCandidate->Phi());
	    hGammaPt[fCentralityBin]->Fill(fCurrentMotherKFCandidate->Pt());
	    hGammaConversionMapXY->Fill(fCurrentMotherKFCandidate->GetConversionX(),fCurrentMotherKFCandidate->GetConversionY());
	    hGammaConversionMapZR->Fill(fCurrentMotherKFCandidate->GetConversionZ(),fCurrentMotherKFCandidate->GetConversionRadius());

       }
	delete fCurrentMotherKFCandidate;
        fCurrentMotherKFCandidate=NULL;
    }

}

///________________________________________________________________________
Bool_t AliV0ReaderV1::ArmenterosQtCut(AliConversionPhotonBase *fPhotonCandidate)
{
    if(fDoHighPtQtGammaSelection){
	if(fPhotonCandidate->GetPhotonPt() < fPtBorderForQt){
	    if(fPhotonCandidate->GetArmenterosQt()>fQtMax){
		hV0CutQt->Fill(fPhotonCandidate->GetPhotonMass());
                return kFALSE;
	    }
	} else {
	    if(fPhotonCandidate->GetArmenterosQt()>fHighPtQtMax){
		hV0CutQt->Fill(fPhotonCandidate->GetPhotonMass());
                return kFALSE;
	    }
	}
    } else {

	if(fPhotonCandidate->GetArmenterosQt()>fQtMax){
	   hV0CutQt->Fill(fPhotonCandidate->GetPhotonMass());
	   return kFALSE;
	}
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::AcceptanceCuts(AliConversionPhotonBase *fPhotonCandidate)
{
    AliVTrack *fCurrentNegativeTrack=GetTrack(fPhotonCandidate->GetTrackLabelNegative());
    AliVTrack *fCurrentPositiveTrack=GetTrack(fPhotonCandidate->GetTrackLabelPositive());

    if(fPhotonCandidate->GetConversionRadius()>fMaxR){ // cuts on distance from collision point
	hV0CutR->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
	     }

    if(fPhotonCandidate->GetConversionRadius()<fMinR){ // cuts on distance from collision point
	hV0CutMinR->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
    }

   if(fPhotonCandidate->GetConversionRadius() <= ((TMath::Abs(fPhotonCandidate->GetConversionZ())*fLineCutZRSlope)-fLineCutZValue)){
	hV0CutLine->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
   }
   /*else if (fUseEtaMinCut &&  fPhotonCandidate->GetConversionRadius() >= ((TMath::Abs(fPhotonCandidate->GetConversionZ())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
	hV0CutLine->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
    }*/
   

    if(TMath::Abs(fPhotonCandidate->GetConversionZ()) > fMaxZ ){ // cuts out regions where we do not reconstruct
	hV0CutZ->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
    }

    if(TMath::Abs(fPhotonCandidate->GetPhotonEta())> fEtaCut || TMath::Abs(fPhotonCandidate->GetPhotonEta())< fEtaCutMin){
	hV0CutEta->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
    }

    if(TMath::Abs(fCurrentNegativeTrack->Eta())> fEtaCut || TMath::Abs(fCurrentNegativeTrack->Eta())< fEtaCutMin){
	hV0CutEta->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
    }

    if(TMath::Abs(fCurrentPositiveTrack->Eta())> fEtaCut || TMath::Abs(fCurrentPositiveTrack->Eta())< fEtaCutMin){
	hV0CutEta->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
    }

    if( fCurrentNegativeTrack->Pt()< fSinglePtCut ||	fCurrentNegativeTrack->Pt()< fSinglePtCut){
	hV0CutSinglePt->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
    }
		    

    if(fPhotonCandidate->GetPhotonPt()<fPtCut){
	hV0CutPt->Fill(fPhotonCandidate->GetPhotonMass());
	return kFALSE;
    }
return kTRUE;
}



///________________________________________________________________________
Bool_t AliV0ReaderV1::TrackCuts(AliConversionPhotonBase *fPhotonCandidate){

    Bool_t passtrackcuts=kTRUE;

  

    if(fESDEvent){

	AliESDtrack *fCurrentNegativeESDTrack=(AliESDtrack*)fESDEvent->GetTrack(fPhotonCandidate->GetTrackLabelNegative());
	AliESDtrack *fCurrentPositiveESDTrack=(AliESDtrack*)fESDEvent->GetTrack(fPhotonCandidate->GetTrackLabelPositive());

          if(!fCurrentNegativeESDTrack||!fCurrentPositiveESDTrack)return kFALSE;

       	// avoid like sign
        if(fCurrentNegativeESDTrack->Charge() == fCurrentPositiveESDTrack->Charge()){

	    hV0CutLikeSign->Fill(fPhotonCandidate->GetPhotonMass());
	    passtrackcuts=kFALSE;
	}

	if( (!(fCurrentNegativeESDTrack->IsOn(AliESDtrack::kTPCrefit))||(!(fCurrentPositiveESDTrack->IsOn(AliESDtrack::kTPCrefit))))){
	    hV0CutRefit->Fill(fPhotonCandidate->GetPhotonMass());
	    passtrackcuts=kFALSE;
         
	}

	if( fCurrentNegativeESDTrack->GetKinkIndex(0) > 0 ||
	   fCurrentPositiveESDTrack->GetKinkIndex(0) > 0) {
	    passtrackcuts=kFALSE;
	}

	if(fCurrentNegativeESDTrack->GetNcls(1) < fMinClsTPC ||	fCurrentPositiveESDTrack->GetNcls(1) < fMinClsTPC ){
	    passtrackcuts=kFALSE;
	    hV0CutKinks->Fill(fPhotonCandidate->GetPhotonMass());
	}

             /*	Double_t negclsToF = 0.;
		if (!fUseCorrectedTPCClsInfo ){
			if(fCurrentNegativeESDTrack->GetTPCNclsF()!=0	){
				negclsToF = (Double_t)fCurrentNegativeESDTrack->GetNcls(1)/(Double_t)fCurrentNegativeESDTrack->GetTPCNclsF();
			}
		} else {
			negclsToF = fCurrentNegativeESDTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(GetXYRadius()));
		}

		Double_t posclsToF = 0.;
		if (!fUseCorrectedTPCClsInfo ){
			if(fCurrentTrack->GetTPCNclsF()!=0	){
				posclsToF = (Double_t)fCurrentTrack->GetNcls(1)/(Double_t)fCurrentTrack->GetTPCNclsF();
			}
		}else{
			posclsToF = fCurrentTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(GetXYRadius()));
		}

		if( negclsToF < fMinClsTPCToF ||	posclsToF < fMinClsTPCToF ){
		    hV0CutMinNclsTPCToF->Fill(fPhotonCandidate->GetPhotonMass());
		    passtrackcuts=kFALSE; }
               */

		  
			    

    }

    if(fAODEvent){

	AliAODTrack *fCurrentNegativeESDTrack=(AliAODTrack*)fAODEvent->GetTrack(fPhotonCandidate->GetTrackLabelNegative());
	AliAODTrack *fCurrentPositiveESDTrack=(AliAODTrack*)fAODEvent->GetTrack(fPhotonCandidate->GetTrackLabelPositive());

        if(!fCurrentNegativeESDTrack||!fCurrentPositiveESDTrack)return kFALSE;

	// avoid like sign
	if(fCurrentNegativeESDTrack->Charge() == fCurrentPositiveESDTrack->Charge()){

	    hV0CutLikeSign->Fill(fPhotonCandidate->GetPhotonMass());
	    passtrackcuts=kFALSE;
	}

	if( !(fCurrentNegativeESDTrack->IsOn(AliESDtrack::kTPCrefit))){
	    hV0CutRefit->Fill(fPhotonCandidate->GetPhotonMass());
	    passtrackcuts=kFALSE;
	}

	if( !(fCurrentPositiveESDTrack->IsOn(AliESDtrack::kTPCrefit))){
	    hV0CutRefit->Fill(fPhotonCandidate->GetPhotonMass());
	    passtrackcuts=kFALSE;
	}

	// to be implemented
	/*
	 if( fCurrentNegativeESDTrack->GetKinkIndex(0) > 0 ||
	 fCurrentPositiveESDTrack->GetKinkIndex(0) > 0) {
	 }*/

	if(fCurrentNegativeESDTrack->GetNcls(1) < fMinClsTPC ||	fCurrentPositiveESDTrack->GetNcls(1) < fMinClsTPC ){
	    passtrackcuts=kFALSE;}


    }

    return passtrackcuts;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::dEdxCuts(AliConversionPhotonBase *fPhotonCandidate){

    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    AliPIDResponse *fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();

    AliVTrack *fCurrentTrack=0x0;

    for(Int_t ilabel=0;ilabel<2;ilabel++){

	fCurrentTrack=GetTrack(fPhotonCandidate->GetTrackLabel(ilabel));

	if(!fCurrentTrack)return kFALSE;

	// Fill dEdx before cuts

	hGammadEdxbefore->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());

     	 if(fDodEdxSigmaCut == kTRUE){
	     if( fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaBelowElectronLine ||
		fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaAboveElectronLine){

		hV0CutdEdxElectron->Fill(fPhotonCandidate->GetPhotonMass());
		 return kFALSE;
	     }
	     
	     if( fCurrentTrack->P()>fPIDMinPnSigmaAbovePionLine && fCurrentTrack->P()<fPIDMaxPnSigmaAbovePionLine ){
		 if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
		    fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
		    fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){

                    hV0CutdEdxPion->Fill(fPhotonCandidate->GetPhotonMass());
		     return kFALSE;
		 }
	     }
			
     	     // High Pt
	     if( fCurrentTrack->P()>fPIDMaxPnSigmaAbovePionLine ){
		 if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
		    fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
		    fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLineHighPt){

		    hV0CutdEdxPion->Fill(fPhotonCandidate->GetPhotonMass());
		     return kFALSE;
		 }
	     }
	 }

	 if(fDoKaonRejectionLowP == kTRUE){
	     if(fCurrentTrack->P()<fPIDMinPKaonRejectionLowP ){
		 if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kKaon))<fPIDnSigmaAtLowPAroundKaonLine){
                    hV0CutdEdxKaonLowP->Fill(fPhotonCandidate->GetPhotonMass());
                    return kFALSE;
	         }
	     }
	 }

	 if(fDoProtonRejectionLowP == kTRUE){
	     if( fCurrentTrack->P()<fPIDMinPProtonRejectionLowP ){
		 if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kProton))<fPIDnSigmaAtLowPAroundProtonLine){
		    hV0CutdEdxProtonLowP->Fill(fPhotonCandidate->GetPhotonMass());
		     return kFALSE;
		 }
	     }
	 }

	 if(fDoPionRejectionLowP == kTRUE){
	     if( fCurrentTrack->P()<fPIDMinPPionRejectionLowP ){
		 if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion))<fPIDnSigmaAtLowPAroundPionLine){
                    hV0CutdEdxPionLowP->Fill(fPhotonCandidate->GetPhotonMass());
                 return kFALSE;
		 }
	     }
	 }


	 if( fDoTOFsigmaCut == kTRUE ){ // RRnewTOF start /////////////////////////////////////////////////////////////////////////////

	    if((fPIDResponse->NumberOfSigmasTOF(fCurrentTrack,AliPID::kElectron)>fTofPIDnSigmaAboveElectronLine) || (fPIDResponse->NumberOfSigmasTOF(fCurrentTrack,AliPID::kElectron)<fTofPIDnSigmaBelowElectronLine)){
              hV0CutdEdxTOFElectron->Fill(fPhotonCandidate->GetPhotonMass());
              return kFALSE;
	    }
	 } /////////////////////////////// RRnewTOF end ///////////////////////////////////////////////////////////////////////////////

         // Apply TRD PID
	 if(fDoTRDPID){
	     if(!fPIDResponse->IdentifiedAsElectronTRD(fCurrentTrack,fPIDTRDEfficiency)){
		 hV0CutdEdxTRD->Fill(fPhotonCandidate->GetPhotonMass());
		 return kFALSE;
	     }
	 }

	 // Fill dEdx Histogram after Cuts

         hGammadEdxafter->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());

    }

    return kTRUE;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::AsymmetryCut(AliConversionPhotonBase *fPhotonCandidate)
{
   
    for(Int_t ii=0;ii<2;ii++){
	AliVTrack *fCurrentTrack=GetTrack(fPhotonCandidate->GetTrackLabel(ii));

	if( fCurrentTrack->P()>fMinPPhotonAsymmetryCut ){
	    Double_t trackNegAsy=0;
	    if (fPhotonCandidate->GetPhotonP()!=0.){
		trackNegAsy= fCurrentTrack->P()/fPhotonCandidate->GetPhotonP();
	    }
	    if( trackNegAsy<fMinPhotonAsymmetry ||trackNegAsy>(1.- fMinPhotonAsymmetry)){
		hV0CutAsymmetry->Fill(fPhotonCandidate->GetPhotonMass());
		return kFALSE;
	    }
	}
    }

    return kTRUE;
}


///________________________________________________________________________

Int_t AliV0ReaderV1::GetNumberOfContributorsVtx(){
	// returns number of contributors to the vertex
    if(fESDEvent){
	if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors()>0) {
	    return fESDEvent->GetPrimaryVertexTracks()->GetNContributors();
	}
     
	if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors()<1) {
	    //		return 0;
	    //-AM test pi0s without SPD only vertex
	    if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
		return fESDEvent->GetPrimaryVertexSPD()->GetNContributors();

	    }
	    if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()<1) {
	     //   cout<<"number of contributors from bad vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
		return 0;
	    }
	}
    }
    if(fAODEvent){
          if(fAODEvent->GetPrimaryVertex()->GetNContributors()>0) {
	    return fESDEvent->GetPrimaryVertex()->GetNContributors();
	}
        if(fAODEvent->GetPrimaryVertex()->GetNContributors()<1) {
	    if(fAODEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
		return fAODEvent->GetPrimaryVertexSPD()->GetNContributors();

	    }
	    if(fAODEvent->GetPrimaryVertexSPD()->GetNContributors()<1) {
		AliWarning(Form("Number of contributors from bad vertex type:: %s",fAODEvent->GetPrimaryVertex()->GetName()));
		return 0;
	    }
	}
    }

    return 0;

}

///________________________________________________________________________
Double_t AliV0ReaderV1::GetCentrality(){
    if(fAODEvent){
	if(fAODEvent->GetHeader()){return fAODEvent->GetHeader()->GetCentrality();}
    }

    if(fESDEvent){
	return ((AliCentrality *)fESDEvent->GetCentrality())->GetCentralityPercentile("V0M");
    }

    return -1;
}

///________________________________________________________________________

AliEventplane *AliV0ReaderV1::GetEventPlane(){
  if(fESDEvent){
	   return fESDEvent->GetEventplane();
  }
  if(fAODEvent){
      return fAODEvent->GetEventplane();
  }
return 0x0;
}

//________________________________________________________________________
Bool_t AliV0ReaderV1::SetEventPlane()
{
    if(fMCStack||!fIsHeavyIon){
	// NO EP in MC mode and pp mode
	fEPAngle=0;
       return kTRUE;
    }

    if(fIsHeavyIon){
	AliEventplane *fEP=GetEventPlane();

	if(!fEP)return kFALSE;

	fEPAngle=fEP->GetEventplane("Q");

        return kTRUE;
    }
    return kFALSE;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::EventCuts(){

    fEventIsSelected=kTRUE;

    // Z Vertex Position Cut

    if(!VertexZCut()){
	hV0EventCuts->Fill(0);
	fEventIsSelected=kFALSE;;
    }
    // Number of Contributors Cut

    if(GetNumberOfContributorsVtx()<=0) {
	hV0EventCuts->Fill(1);
	fEventIsSelected=kFALSE;;
    }
    // Centrality Selection

    if(!CentralitySelection()){
	hV0EventCuts->Fill(2);
	fEventIsSelected=kFALSE;;
    }

    // Event Plane
    if(!SetEventPlane()){hV0EventCuts->Fill(4);}

    // Fill Event Histograms

    if(fEventIsSelected){
	  // Fill Event Histograms
	  hV0EventCuts->Fill(9);
	  hVertexZ->Fill(fVertexZ);
	  hCentrality->Fill(fCentrality);
	  hNEvents->Fill(fCentralityBin);
    }
    else{hV0EventCuts->Fill(8);}

    return fEventIsSelected;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::VertexZCut(){

    fVertexZ=GetPrimaryVertex()->GetZ();

    if(fBGHandler){
	if(fBGHandler->GetZBinIndex(fVertexZ)<0)return kFALSE;
    }
    else{
	if(fVertexZ>fMaxVertexZ)return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
AliVTrack *AliV0ReaderV1::GetTrack(Int_t label){
    if(fESDEvent){
	return (AliESDtrack*)fESDEvent->GetTrack(label);
    }
    if(fAODEvent)return (AliAODTrack*)fAODEvent->GetTrack(label);
    return 0x0;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::PIDProbabilityCut(AliConversionPhotonBase *fPhotonCandidate){

    if(fESDEvent){

	Bool_t iResult=kFALSE;

	Double_t *posProbArray = new Double_t[AliPID::kSPECIES];
	Double_t *negProbArray = new Double_t[AliPID::kSPECIES];

	AliESDtrack* negTrack	= (AliESDtrack*)fESDEvent->GetTrack(fPhotonCandidate->GetTrackLabelNegative());
	AliESDtrack* posTrack	= (AliESDtrack*)fESDEvent->GetTrack(fPhotonCandidate->GetTrackLabelPositive());

	if(negProbArray && posProbArray){

	    negTrack->GetTPCpid(negProbArray);
	    posTrack->GetTPCpid(posProbArray);

	    if(negProbArray[AliPID::kElectron]>=fPIDProbabilityCutNegativeParticle && posProbArray[AliPID::kElectron]>=fPIDProbabilityCutPositiveParticle){
		iResult=kTRUE;
	    }
	    else{hV0CutPIDProb->Fill(fPhotonCandidate->GetPhotonMass());}
	}

	delete [] posProbArray;
	delete [] negProbArray;
	return iResult;

    }
    if(fAODEvent){
	// not possible to implement
	return kTRUE;}
    return kFALSE;
}

///________________________________________________________________________
const AliVertex *AliV0ReaderV1::GetPrimaryVertex() {
    if(fESDEvent)return fESDEvent->GetPrimaryVertex();
    if(fAODEvent)return const_cast<const AliVertex*>(dynamic_cast<AliVertex*>(fAODEvent->GetPrimaryVertex()));
return 0x0;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::GetHelixCenter(const AliExternalTrackParam *track, Double_t b,Int_t charge, Double_t center[2]){
	// see header file for documentation
	
	Double_t	helix[6];
	track->GetHelixParameters(helix,b);
	
	Double_t xpos =	helix[5];
	Double_t ypos =	helix[0];
	Double_t radius = TMath::Abs(1./helix[4]);
	Double_t phi = helix[2];

	if(phi < 0){
		phi = phi + 2*TMath::Pi();
	}

	phi -= TMath::Pi()/2.;
	Double_t xpoint =	radius * TMath::Cos(phi);
	Double_t ypoint =	radius * TMath::Sin(phi);

	if(b<0){
		if(charge > 0){
			xpoint = - xpoint;
			ypoint = - ypoint;
		}

		if(charge < 0){
			xpoint =	xpoint;
			ypoint =	ypoint;
		}
	}
	if(b>0){
		if(charge > 0){
			xpoint =	xpoint;
			ypoint =	ypoint;
		}

		if(charge < 0){
			xpoint = - xpoint;
			ypoint = - ypoint;
		}
	}
	center[0] =	xpos + xpoint;
	center[1] =	ypos + ypoint;

	return 1;
}
///________________________________________________________________________
Bool_t AliV0ReaderV1::GetConversionPoint(const AliExternalTrackParam *pparam,const AliExternalTrackParam *nparam,Double_t convpos[3]){

    if(!pparam||!nparam)return kFALSE;

	Double_t helixcenterpos[2];
	GetHelixCenter(pparam,GetMagneticField(),pparam->Charge(),helixcenterpos);

	Double_t helixcenterneg[2];
	GetHelixCenter(nparam,GetMagneticField(),nparam->Charge(),helixcenterneg);

	Double_t helixpos[6];
	pparam->GetHelixParameters(helixpos,GetMagneticField());
	Double_t posradius = TMath::Abs(1./helixpos[4]);

	Double_t helixneg[6];
	nparam->GetHelixParameters(helixneg,GetMagneticField());
	Double_t negradius = TMath::Abs(1./helixneg[4]);

        // Calculate xy-position

	Double_t xpos = helixcenterpos[0];
	Double_t ypos = helixcenterpos[1];
	Double_t xneg = helixcenterneg[0];
	Double_t yneg = helixcenterneg[1];

	convpos[0] = (xpos*negradius + xneg*posradius)/(negradius+posradius);
	convpos[1] = (ypos*negradius+	yneg*posradius)/(negradius+posradius);


	// Calculate z-position

	Double_t deltaXPos = convpos[0] -	xpos;
	 Double_t deltaYPos = convpos[1] -	ypos;

	 Double_t deltaXNeg = convpos[0] -	xneg;
	 Double_t deltaYNeg = convpos[1] -	yneg;

	 Double_t alphaPos =	TMath::Pi() + TMath::ATan2(-deltaYPos,-deltaXPos);
	 Double_t alphaNeg =	TMath::Pi() + TMath::ATan2(-deltaYNeg,-deltaXNeg);

	 Double_t vertexXNeg =	xneg +	TMath::Abs(negradius)*
	 TMath::Cos(alphaNeg);
	 Double_t vertexYNeg =	yneg +	TMath::Abs(negradius)*
	 TMath::Sin(alphaNeg);

	 Double_t vertexXPos =	xpos +	TMath::Abs(posradius)*
	 TMath::Cos(alphaPos);
	 Double_t vertexYPos =	ypos +	TMath::Abs(posradius)*
	 TMath::Sin(alphaPos);

	 Double_t x0neg =	 helixneg[5];
	 Double_t y0neg =	 helixneg[0];

	 Double_t x0pos =	 helixpos[5];
	 Double_t y0pos =	 helixpos[0];

	 Double_t dNeg = TMath::Sqrt((vertexXNeg -	x0neg)*(vertexXNeg - x0neg)
															 +(vertexYNeg -	y0neg)*(vertexYNeg - y0neg));

	 Double_t dPos = TMath::Sqrt((vertexXPos -	x0pos)*(vertexXPos - x0pos)
															 +(vertexYPos -	y0pos)*(vertexYPos - y0pos));

	 Double_t rNeg =	TMath::Sqrt(negradius*negradius -
	 dNeg*dNeg/4.);

	 Double_t rPos = TMath::Sqrt(posradius*posradius -
	 dPos*dPos/4.);

	 Double_t deltabetaNeg =	2*(TMath::Pi() +	 TMath::ATan2(-dNeg/2.,-rNeg));
	 Double_t deltabetaPos = 2*(TMath::Pi() + TMath::ATan2(-dPos/2.,-rPos));

	 Double_t deltaUNeg = negradius*deltabetaNeg;
	 Double_t deltaUPos = posradius*deltabetaPos;

	 Double_t zphaseNeg = nparam->GetZ() +	deltaUNeg * nparam->GetTgl();
	 Double_t zphasePos = pparam->GetZ() +	deltaUPos * pparam->GetTgl();

	 convpos[2] = (zphasePos*negradius+zphaseNeg*posradius)/(negradius+posradius);

return kTRUE;
}
///________________________________________________________________________
void AliV0ReaderV1::ProcessMC(AliKFConversionPhoton *fCurrentReconstructedGamma){

    if(!fMCStack)return;

    TParticle *fMotherMCParticle=NULL;
    TParticle *fNegativeMCParticle=NULL;
    TParticle *fPositiveMCParticle=NULL;
    
    // Get MC Particles
    fMotherMCParticle   = fCurrentReconstructedGamma->GetMCParticle(fMCStack);
    fNegativeMCParticle = fCurrentReconstructedGamma->GetNegativeMCDaughter(fMCStack);
    fPositiveMCParticle = fCurrentReconstructedGamma->GetPositiveMCDaughter(fMCStack);

    if(fPositiveMCParticle&&fNegativeMCParticle&&fMotherMCParticle){

	// Check if it is a true photon

	if(fMotherMCParticle->GetPdgCode()==22){

     	    hMCPtRECOTRUE[fCentralityBin]->Fill(fCurrentReconstructedGamma->GetPt());

	    // Pt Resolution
          
	    Double_t mcpt	 = fMotherMCParticle->Pt();
	    Double_t esdpt	= fCurrentReconstructedGamma->GetPt();
	    Double_t resdPt = 0.;
	    if(mcpt > 0){
		resdPt = ((esdpt - mcpt)/mcpt)*100.;
	    } else if(mcpt < 0){
		AliWarning("Pt of MC particle is negative, this will cause wrong calculation of resPt");
	    }


	    hMCPtResolution->Fill(mcpt,resdPt);
	    hMCPtResolutionPhi->Fill(fMotherMCParticle->Phi(),resdPt);

	    // Conversion Point Resolution

	    Double_t resdR = 0.;
	    if(fNegativeMCParticle->R() != 0){
		resdR = ((fCurrentReconstructedGamma->GetConversionRadius() - fNegativeMCParticle->R())/fNegativeMCParticle->R())*100.;
	    }
	    Double_t resdRAbs = 0.;
	    resdRAbs = (fCurrentReconstructedGamma->GetConversionRadius() - fNegativeMCParticle->R());

	    hMCRResolutionvsR->Fill(fNegativeMCParticle->R(),resdRAbs);

	    // 		fHistograms->FillHistogram("Resolution_dR", fV0Reader->GetNegativeMCParticle()->R(), resdR);
	    // 		fHistograms->FillHistogram("Resolution_MC_R", fV0Reader->GetNegativeMCParticle()->R());
	    // 		fHistograms->FillHistogram("Resolution_ESD_R", fV0Reader->GetXYRadius());
	    // 		fHistograms->FillHistogram("Resolution_R_dPt", fV0Reader->GetNegativeMCParticle()->R(), resdPt);

	    Double_t resdZ = 0.;
	    if(fNegativeMCParticle->Vz() != 0){
		resdZ = ((fCurrentReconstructedGamma->GetZ() -fNegativeMCParticle->Vz())/fNegativeMCParticle->Vz())*100.;
	    }
	    Double_t resdZAbs = 0.;
	    resdZAbs = fCurrentReconstructedGamma->GetZ() -fNegativeMCParticle->Vz();

	    hMCZResolutionvsZ->Fill( fNegativeMCParticle->Vz(), resdZAbs);
	}
    }

}

///________________________________________________________________________
void AliV0ReaderV1::ProcessMCGammasForEfficiency(){

    if(!fMCStack)return;

    for (Int_t iTracks = 0; iTracks < fMCStack->GetNprimary(); iTracks++) {
	TParticle* particle = (TParticle *)fMCStack->Particle(iTracks);

	//process the gammas

	if(IsMCConversionGammaInAcceptance(particle)){

	    hMCPtTRUE[fCentralityBin]->Fill(particle->Pt());

	}
    }
}


///________________________________________________________________________
Bool_t AliV0ReaderV1::IsMCConversionGammaInAcceptance(TParticle *particle){
    if(!fMCStack)return kFALSE;

    if (particle->GetPdgCode() == 22){
	if(TMath::Abs(particle->Eta())> fEtaCut || TMath::Abs(particle->Eta())< fEtaCutMin)	return kFALSE;

	if(particle->GetMother(0) >-1 && fMCStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
	    return kFALSE; // no photon as mothers!
	}

	if(particle->GetMother(0) >= fMCStack->GetNprimary()){
	    return kFALSE; // the gamma has a mother, and it is not a primary particle
	}

	// looking for conversion (electron + positron from pairbuilding (= 5) )
	TParticle* ePos = NULL;
	TParticle* eNeg = NULL;

	if(particle->GetNDaughters() >= 2){
	    for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
		TParticle *tmpDaughter = fMCStack->Particle(daughterIndex);
		if(tmpDaughter->GetUniqueID() == 5){
		    if(tmpDaughter->GetPdgCode() == 11){
			eNeg = tmpDaughter;
		    } else if(tmpDaughter->GetPdgCode() == -11){
			ePos = tmpDaughter;
		    }
		}
	    }
	}

	if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
	    return kFALSE;
	}

	if(AcceptanceCut(particle,ePos,eNeg))return kTRUE;
    }
    return kFALSE;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::AcceptanceCut(TParticle *particle, TParticle * ePos,TParticle* eNeg){

    // cuts on distance from collision point

    if(particle->R()>fMaxR){
	return kFALSE;}

    if(ePos->R()>fMaxR){
	return kFALSE;
    }

    if(ePos->R()<fMinR){
	return kFALSE;
    }

    if( ePos->R() <= ((TMath::Abs(ePos->Vz())*fLineCutZRSlope)-fLineCutZValue)){
	return kFALSE;
    }
    /*else if (fUseEtaMinCut &&  ePos->R() >= ((TMath::Abs(ePos->Vz())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
	return kFALSE;
    } */

    if(TMath::Abs(eNeg->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
	return kFALSE;
    }

    if(eNeg->Vz()!=ePos->Vz()||eNeg->R()!=ePos->R()){
	return kFALSE;
    }

    if(TMath::Abs(ePos->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
	return kFALSE;
    }

    if(TMath::Abs(particle->Eta())> fEtaCut || TMath::Abs(particle->Eta())< fEtaCutMin){
	return kFALSE;
    }

    if(TMath::Abs(ePos->Eta())> fEtaCut || TMath::Abs(ePos->Eta())< fEtaCutMin){
	return kFALSE;
    }

    if(TMath::Abs(eNeg->Eta())> fEtaCut || TMath::Abs(eNeg->Eta())< fEtaCutMin){
	return kFALSE;
    }

    if( ePos->Pt()< fSinglePtCut ||  eNeg->Pt()< fSinglePtCut){
	return kFALSE;
    }

    if(particle->Pt()<fPtCut){
	return kFALSE;
    }

    return kTRUE;
}
///________________________________________________________________________
void AliV0ReaderV1::PrintCuts(){

    cout<<"V0 Reader initialized with following settings"<<endl;


    cout<<"Acceptance Eta:"<<endl;
    cout<<fEtaCutMin<<" < eta < "<<fEtaCut<<endl;
    cout<<"Conversion Point"<<endl;
    cout<<"Z <"<<fMaxZ<<endl;
    cout<<fMinR<<" < R < "<<fMaxR<<endl;
    cout<<"Line Cut Slope"<<fLineCutZRSlope<<" ZValue "<<fLineCutZValue<<endl;

    cout<<"Pt Gamma > "<<fPtCut<<endl;
    cout<<"Pt Daughters > "<<fSinglePtCut<<endl;

    cout<<"Armenteros Qt Cut"<<endl;

    if(fDoHighPtQtGammaSelection){
	cout<<" qt < "<<fQtMax<<" for pt < "<<fPtBorderForQt<<endl;
	cout<<" qt < "<<fHighPtQtMax<<" for pt > "<<fPtBorderForQt<<endl;
    }
    else{
         cout<<" qt < "<<fQtMax<<endl;
    }
    cout<<"Chi2perNDF > "<<fChi2CutConversion<<endl;



}

//_______________________________________________________________________

Bool_t AliV0ReaderV1::CentralitySelection(){

    fCentralityBin=-1;

    if(!fIsHeavyIon){fCentralityBin=0;fCentrality=0;return kTRUE;}

    fCentrality=GetCentrality();

    if(fIsHeavyIon){

	if(fBGHandler){
	    fCentralityBin=fBGHandler->GetCentralityBinIndex(Int_t(fCentrality));
	}
	else{
	    Double_t fCentralityBins[fNCentralityBins+1];
	    for(int i=0;i<fNCentralityBins;i++){
		fCentralityBins[i]=i*100/Double_t(fNCentralityBins);
	    }

	    for(int i=0;i<fNCentralityBins;i++){
		if(fCentrality>fCentralityBins[i]&&fCentrality<fCentralityBins[i+1]){
		    fCentralityBin=i;
		    return kTRUE;}
	    }
	}

    }
    if(fCentralityBin>=0&&fCentrality>=0){

        return kTRUE;
    }

    AliWarning("Centrality not defined");
    return kFALSE;

}

//________________________________________________________________________
void AliV0ReaderV1::Terminate(Option_t *)
{
 
    printf(Form("V0ReaderV1: V0s processed: %4.0f reconstructed photons: %4.0f \n",hV0CurrentFinder->GetEntries(),hV0Good->GetEntries()));

}
