#include "TList.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TFile.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "AliEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliAnalysisUtils.h"
#include "AliExternalTrackParam.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODPidHF.h"
#include "AliAODMCParticle.h"
#include "AliAODTracklets.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliESDUtils.h"
#include "AliMCEventHandler.h"
#include "AliPIDResponse.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEventHandler.h"
#include "AliVVertex.h"
#include "AliVParticle.h"
#include "AliAODMCHeader.h"
#include "AliJetContainer.h"
#include "AliGenEventHeader.h"
#include "AliVertexerTracks.h"
#include "AliEmcalList.h"
#include "THnSparse.h"
#include "TObjectTable.h"
#include "AliAnalysisTaskEmcalJet.h"


//***********************************//
#include "AliAnalysisTaskHFJetIPQA.h"
//***********************************//
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include "AliAnalysisHelperJetTasks.h"
#include "AliGenPythiaEventHeader.h"
#include "TChain.h"
#include <map>
using std::min;
using std::cout;
using std::endl;
using std::vector;
using std::pair;
using std::map;
ClassImp(AliAnalysisTaskHFJetIPQA)

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<bool>+;
#endif


AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA():
AliAnalysisTaskEmcalJet(),
fEventCuts(0),
fHistManager(),
fEventVertex(nullptr),
fPidResponse(nullptr),
tJetTree(nullptr),     //! tree containing jet properties
fJetRecPt(-99),
fJetArea(-99),
fMatchedJetPt(-99),
fJetProb(-99),
fMeanLNKt(-99),
fMeanTheta(-99),
fMeanLNKtSD(-99),
fMeanThetaSD(-99),
fJetMass(-99),
fJetFlavour(-99),
nTracks(0),
fNEvent(0),
bMatched(kFALSE),
bIsTrueGenV0Jet(0),
fTrackIPs{0},
fTrackIPSigs{0},
fTrackProb{0},
fTrackChi2OverNDF{0},
fTrackPt{0},
fDeltaRij{0},
fV0MotherPt{0},
fV0MotherPtMC{0},
fV0MotherEta{0},
fV0MotherEtaMC{0},
iTrackITSHits{0},
iV0MCID{0},
iV0RecID{0},
bTrackIsV0{0},
bPassedSD{0},
//bFull{0},
bSingle1st{0},
//bSingle2nd{0},
//bSingle3rd{0},
bDouble{0},
//bTriple{0},
jetconrec(nullptr),
jetcongen(nullptr),
fRunSmearing(kFALSE),
fUsePIDJetProb(kFALSE),
fDoMCCorrection(kFALSE),
fDoUnderlyingEventSub(kFALSE),
fDoFlavourMatching(kFALSE),
fParam_Smear_Sigma(1.),
fParam_Smear_Mean(0.),
fGlobalVertex(kFALSE),
fDoNotCheckIsPhysicalPrimary(kFALSE),
fDoJetProb(kFALSE),
fFillCorrelations(kFALSE),
fDoLundPlane(kFALSE),
fDoTCTagging(0),
fDoProbTagging(0),
fDoMCEffs(0),
fUseSignificance(kTRUE),
fResponseMode(kFALSE),
kTagLevel(3),
fFracs(0),
fXsectionWeightingFactor(1),
fProductionNumberPtHard(-1),
fNThresholds(1),
fNTrackTypes(3),
sTemplateFlavour(0),
fUnfoldPseudeDataFrac(50),
sTaskName(""),
fJetRadius(0.4),
fDaughtersRadius(1),
fNoJetConstituents(0),
fTCThresholdPtFixed(0.008),
h1DThresholdsFirst(0),
h1DThresholdsSecond(0),
h1DThresholdsThird(0),
h2DProbLookup(0),
h1DProbThresholds(0),
cCuts(0),
fh1DCutInclusive(0),
fh1dCutudg(0),
fh1dCutc(0),
fh1dCutb(0),
fh1dCuts(0),
fh1dTracksAccepeted(0),
fh1dCutsPrinted(0),
fHLundIterative(nullptr),
fhnV0InJetK0s(nullptr),
fhnV0InJetLambda(nullptr),
fhnV0InJetALambda(nullptr),
fh1V0CounterCentK0s(nullptr),
fh1V0CounterCentLambda(nullptr),
fh1V0CounterCentALambda(nullptr),
fh2dKshortMassVsPt(nullptr),
fh2dLamdaMassVsPt(nullptr),
fh2dAnLamdaMassVsPt(nullptr),
h1DV0FalseRec(nullptr),
h1DV0TrueRec(nullptr),
h1DV0TrueDataDef(nullptr),
h1DV0TrueMCDef(nullptr),
fh1dKshortPtMC(nullptr),
fh1dLamdaPtMC(nullptr),
fh1dAnLamdaPtMC(nullptr),
fh1dKshortEtaMC(nullptr),
fh1dLamdaEtaMC(nullptr),
fh1dAnLamdaEtaMC(nullptr),
fh2dKshortPtVsJetPtMC(nullptr),
fh2dLamdaPtVsJetPtMC(nullptr),
fh2dAnLamdaPtVsJetPtMC(nullptr),
fMCArray(nullptr),
fMCEvent(nullptr),
fESDTrackCut(nullptr),
fVertexer(nullptr),
fV0CandidateArray(nullptr),
fMcEvtSampled(kFALSE),
fBackgroundFactorLinus{0},
fPUdsgJet(100),fPSJet(100),fPCJet(100),fPBJet(100),
fJetCont(10),
fAnalysisCuts{0},
fCombined(nullptr),
fMCglobalDCAxyShift(0.0008),
fMCglobalDCASmear(1),
fVertexRecalcMinPt(1.0),
fHardCutOff(0),
fn1_mix(-999.),
fn2_mix(-999.),
fn3_mix(-999.),
fIsMixSignalReady_n1(kFALSE),
fIsMixSignalReady_n2(kFALSE),
fIsMixSignalReady_n3(kFALSE),
fIsSameEvent_n1(kFALSE),
fIsSameEvent_n2(kFALSE),
fIsSameEvent_n3(kFALSE),
fUseTreeForCorrelations(kFALSE),
fCorrelationCrossCheck(nullptr),
fTREE_n1(-99.),
fTREE_n2(-99.),
fTREE_n3(-99.),
fTREE_pt(-1.)

{
    SetMakeGeneralHistograms(kTRUE);
    SetDefaultAnalysisCuts();
    SetDefaultV0Cuts();
    SetNeedEmcalGeom(kFALSE);
    SetOffTrigger(AliVEvent::kINT7);
    SetVzRange(-10,10);
    DefineOutput(1,  AliEmcalList::Class());
    DefineOutput(2, TTree::Class());
}
AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA(const char *name):
AliAnalysisTaskEmcalJet(name, kTRUE),
fEventCuts(0),
fHistManager(name),
fEventVertex(nullptr),
fPidResponse(nullptr),
tJetTree(nullptr),   //! tree containing jet properties
fJetRecPt(-99),
fJetArea(-99),
fMatchedJetPt(-99),
fJetProb(-99),
fMeanLNKt(-99),
fMeanTheta(-99),
fMeanLNKtSD(-99),
fMeanThetaSD(-99),
fJetMass(-99),
fJetFlavour(-99),
nTracks(0),
fNEvent(0),
bMatched(kFALSE),
bIsTrueGenV0Jet(0),
fTrackIPs{0},
fTrackIPSigs{0},
fTrackProb{0},
fTrackChi2OverNDF{0},
fTrackPt{0},
fDeltaRij{0},
fV0MotherPt{0},
fV0MotherPtMC{0},
fV0MotherEta{0},
fV0MotherEtaMC{0},
iTrackITSHits{0},
iV0MCID{0},
iV0RecID{0},
bTrackIsV0{0},
bPassedSD{0},
//bFull{0},
bSingle1st{0},
//bSingle2nd{0},
//bSingle3rd{0},
bDouble{0},
//bTriple{0},
jetconrec(nullptr),
jetcongen(nullptr),
fRunSmearing(kFALSE),
fUsePIDJetProb(kFALSE),
fDoMCCorrection(kFALSE),
fDoUnderlyingEventSub(kFALSE),
fDoFlavourMatching(kFALSE),
fParam_Smear_Sigma(1.),
fParam_Smear_Mean(0.),
fGlobalVertex(kFALSE),
fDoNotCheckIsPhysicalPrimary(kFALSE),
fDoJetProb(kFALSE),
fFillCorrelations(kFALSE),
fDoLundPlane(kFALSE),
fDoTCTagging(0),
fDoProbTagging(0),
fDoMCEffs(0),
fUseSignificance(kTRUE),
fResponseMode(kFALSE),
kTagLevel(3),
fFracs(0),
fXsectionWeightingFactor(1.),
fProductionNumberPtHard(-1),
fNThresholds(1),
fNTrackTypes(3),
sTemplateFlavour(0),
fUnfoldPseudeDataFrac(50),
sTaskName(""),
fJetRadius(0.4),
fDaughtersRadius(1),
fNoJetConstituents(0),
fTCThresholdPtFixed(0.008),
h1DThresholdsFirst(0),
h1DThresholdsSecond(0),
h1DThresholdsThird(0),
h2DProbLookup(0),
h1DProbThresholds(0),
cCuts(0),
fh1DCutInclusive(0),
fh1dCutudg(0),
fh1dCutc(0),
fh1dCutb(0),
fh1dCuts(0),
fh1dTracksAccepeted(0),
fh1dCutsPrinted(0),
fHLundIterative(nullptr),
fhnV0InJetK0s(nullptr),
fhnV0InJetLambda(nullptr),
fhnV0InJetALambda(nullptr),
fh1V0CounterCentK0s(nullptr),
fh1V0CounterCentLambda(nullptr),
fh1V0CounterCentALambda(nullptr),
fh2dKshortMassVsPt(nullptr),
fh2dLamdaMassVsPt(nullptr),
fh2dAnLamdaMassVsPt(nullptr),
h1DV0FalseRec(nullptr),
h1DV0TrueRec(nullptr),
h1DV0TrueDataDef(nullptr),
h1DV0TrueMCDef(nullptr),
fh1dKshortPtMC(nullptr),
fh1dLamdaPtMC(nullptr),
fh1dAnLamdaPtMC(nullptr),
fh1dKshortEtaMC(nullptr),
fh1dLamdaEtaMC(nullptr),
fh1dAnLamdaEtaMC(nullptr),
fh2dKshortPtVsJetPtMC(nullptr),
fh2dLamdaPtVsJetPtMC(nullptr),
fh2dAnLamdaPtVsJetPtMC(nullptr),
fMCArray(nullptr),
fMCEvent(nullptr),
fESDTrackCut(nullptr),
fVertexer(nullptr),
fV0CandidateArray(nullptr),
fMcEvtSampled(kFALSE),
fBackgroundFactorLinus{0},
fPUdsgJet(100),fPSJet(100),fPCJet(100),fPBJet(100),
fJetCont(10),
fAnalysisCuts{0},
fCombined(nullptr),
fMCglobalDCAxyShift(0.000668),
fMCglobalDCASmear(1),
fVertexRecalcMinPt(1.0),
fHardCutOff(0),
fn1_mix(-999.),
fn2_mix(-999.),
fn3_mix(-999.),
fIsMixSignalReady_n1(kFALSE),
fIsMixSignalReady_n2(kFALSE),
fIsMixSignalReady_n3(kFALSE),
fIsSameEvent_n1(kFALSE),
fIsSameEvent_n2(kFALSE),
fIsSameEvent_n3(kFALSE),
fUseTreeForCorrelations(kFALSE),
fCorrelationCrossCheck(nullptr),
fTREE_n1(-99.),
fTREE_n2(-99.),
fTREE_n3(-99.),
fTREE_pt(-1.)
{
    SetNeedEmcalGeom(kFALSE);
    SetOffTrigger(AliVEvent::kINT7);
    SetVzRange(-10,10);
    SetMakeGeneralHistograms(kTRUE);
    SetDefaultAnalysisCuts();
    SetDefaultV0Cuts();
    DefineOutput(1,  AliEmcalList::Class()) ;
    DefineOutput(2, TTree::Class());
}

/*! \brief ChangeDefaultCutTo
 *
 *
 * Modify default analysis cuts
 */
void AliAnalysisTaskHFJetIPQA::ChangeDefaultCutTo(AliAnalysisTaskHFJetIPQA::bCuts cutname, Double_t newcutvalue){
    fAnalysisCuts[cutname] =newcutvalue;
}
/*! \brief SetDefaultAnalysisCuts
 *
 *
 * Set default analysis cuts
 */
void AliAnalysisTaskHFJetIPQA::SetDefaultAnalysisCuts(){
    //DCA
    fAnalysisCuts[bAnalysisCut_DCAJetTrack]     = 0.07;
    fAnalysisCuts[bAnalysisCut_MaxDecayLength]  = 5.;
    fAnalysisCuts[bAnalysisCut_MaxDCA_XY]       = 1.;
    fAnalysisCuts[bAnalysisCut_MaxDCA_Z]        = 2.;

    //Vertex
    fAnalysisCuts[bAnalysisCut_NContibutors]    = 3 ;
    //fAnalysisCuts[bAnalysisCut_RelError_Y]      = 0.2;
    //fAnalysisCuts[bAnalysisCut_RelError_Z]      = 0.2;
    //fAnalysisCuts[bAnalysisCut_Sigma_Y]         = 0.3;
    //fAnalysisCuts[bAnalysisCut_Sigma_Z]         = 0.3;
    //fAnalysisCuts[bAnalysisCut_SigmaDiamond]    = 2.;
    fAnalysisCuts[bAnalysisCut_MaxVtxZ]         = 10.;
    fAnalysisCuts[bAnalysisCut_Z_Chi2perNDF]    =3.5*3.5;
    fAnalysisCuts[bAnalysisCut_MinNewVertexContrib] = 1;

    //Tracks
    fAnalysisCuts[bAnalysisCut_MinTrackPt]      =.5;
    //fAnalysisCuts[bAnalysisCut_MinTrackPtMC]    =.5;
    fAnalysisCuts[bAnalysisCut_MinTPCClus]      =100;
    fAnalysisCuts[bAnalysisCut_MinITSLayersHit] =4;
    fAnalysisCuts[bAnalysisCut_MinTrackChi2] =5;
    fAnalysisCuts[bAnalysisCut_HasSDD]=1;
    fAnalysisCuts[bAnalysisCut_HasSSD]=1;
    fAnalysisCuts[bAnalysisCut_HasSPD]=1;
    fAnalysisCuts[bAnalysisCut_HasTPCrefit]=1;
    fAnalysisCuts[bAnalysisCut_HasITSrefit]=1;
    fAnalysisCuts[bAnalysisCut_KinkCand]=1;

    //Jet Cuts
    fAnalysisCuts[bAnalysisCut_MinJetPt]        =0;  //only for settings output. Not really used as cuts are done in .C file
    fAnalysisCuts[bAnalysisCut_MaxJetPt]        =1000;
    fAnalysisCuts[bAnalysisCut_MinJetEta]       =-0.9;
    fAnalysisCuts[bAnalysisCut_MaxJetEta]       =0.9;
    fAnalysisCuts[bAnalysisCut_MinJetArea] =0.6*TMath::Pi()*0.4*0.4;
    fAnalysisCuts[bAnalysisCut_SDz]=  0.1;
    fAnalysisCuts[bAnalysisCut_SDbeta]=0;
    fAnalysisCuts[bAnalysisCut_MaxIPLNJP]=25;

    //Events
    fAnalysisCuts[bAnalysisCut_PtHardAndJetPtFactor] =3;
}

void AliAnalysisTaskHFJetIPQA::SetDefaultV0Cuts(){
    fV0Cuts[DaughMaxEta]=0.8;
    fV0Cuts[DaughMinPt]=0.15;
    fV0Cuts[MinDCADaughWrtPV]=0.06;
    fV0Cuts[MaxDCADaughvsDaugh]=1;
    fV0Cuts[IsTPCRefitOn]=1;
    fV0Cuts[DoPosNoTPCClusters]=1;
    fV0Cuts[MinNoCrossedTPCRows]=70;
    fV0Cuts[NoCrossedOverNoTPCClustersMin]=0.8;
    fV0Cuts[NoCrossedOverNoTPCClustersMax]=1000;
    fV0Cuts[IsKinkCand]=1;
    fV0Cuts[MaxSigmadEdxTPC]=3;

    fV0Cuts[MaxV0Eta]=0;
    fV0Cuts[MaxV0Rap]=0.5;
    fV0Cuts[MinDecayRadius]=0.5;
    fV0Cuts[MaxDecayRadius]=1000;
    fV0Cuts[MaxCosPALambda]=0.995;
    fV0Cuts[MinCosPAK0]=0.97;
    fV0Cuts[MaxLifeTimeK0]=5;
    fV0Cuts[MaxLifeTimeLambda]=5;
    fV0Cuts[DoArmenteros]=1;
    fV0Cuts[DoMassWindow]=1;
    fV0Cuts[InvarMassWindowK0]=0.01;
    fV0Cuts[InvarMassWindowLambda]=0.005;

    fV0Cuts[fAV0Cut]=0;
    fV0Cuts[fBV0Cut]=0;
    fV0Cuts[fCV0Cut]=9999;
}

void AliAnalysisTaskHFJetIPQA::SmearTrack(AliAODTrack *track) {
    if(!fIsPythia) return;
    printf("Run Track Smearing.\n");
    // Get reconstructed track parameters
    AliExternalTrackParam et; et.CopyFromVTrack(track);
    Double_t *param=const_cast<Double_t*>(et.GetParameter());
   // Double_t *covar=const_cast<Double_t*>(et.GetCovariance());
    // Get MC info
    Int_t imc=track->GetLabel();
    if (imc<=0) return;
    const AliAODMCParticle *mc=static_cast<AliAODMCParticle*>(fMCArray->At(imc));
    Double_t mcx[3];
    Double_t mcp[3];
    Double_t mccv[36]={0.};
    Short_t  mcc;
    mc->XvYvZv(mcx);
    mc->PxPyPz(mcp);
    mcc=mc->Charge();
    AliExternalTrackParam mct(mcx,mcp,mccv,mcc);
    const Double_t *parammc=mct.GetParameter();
    AliVertex vtx(mcx,1.,1);
    // Correct reference points and frames according to MC
    // TODO: B-Field correct?
    // TODO: failing propagation....
    et.PropagateToDCA(&vtx,track->GetBz(),10.);
    et.Rotate(mct.GetAlpha());
    // Select appropriate smearing functions
    Double_t sd0rpn=fParam_Smear_Sigma;
    Double_t sd0mrpn=fParam_Smear_Mean;//mu m
    Double_t sd0zn =1.;
    Double_t spt1n =1.;
    Double_t sd0rpo=1;
    Double_t sd0mrpo=0;
    Double_t sd0zo =1.;
    Double_t spt1o =1.;
    // Use the same units (i.e. cm and GeV/c)! TODO: pt!
    sd0rpo*=1.e-4;
    sd0zo *=1.e-4;
    sd0rpn*=1.e-4;
    sd0zn *=1.e-4;
    sd0mrpo*=1.e-4;
    sd0mrpn*=1.e-4;
    // Apply the smearing
    Double_t d0zo  =param  [1];
    Double_t d0zmc =parammc[1];
    Double_t d0rpo =param  [0];
    Double_t d0rpmc=parammc[0];
    Double_t pt1o  =param  [4];
    Double_t pt1mc =parammc[4];
    Double_t dd0zo =d0zo-d0zmc;
    Double_t dd0zn =dd0zo *(sd0zo >0. ? (sd0zn /sd0zo ) : 1.);
    Double_t d0zn  =d0zmc+dd0zn;
    Double_t dd0rpo=d0rpo-d0rpmc;
    Double_t dd0rpn=dd0rpo*(sd0rpo>0. ? (sd0rpn/sd0rpo) : 1.);
    Double_t dd0mrpn=sd0mrpn-sd0mrpo;
    Double_t d0rpn =d0rpmc+dd0rpn-dd0mrpn;
    Double_t dpt1o =pt1o-pt1mc;
    Double_t dpt1n =dpt1o *(spt1o >0. ? (spt1n /spt1o ) : 1.);
    Double_t pt1n  =pt1mc+dpt1n;
    param[0]=d0rpn;
    param[1]=d0zn ;
    param[4]=pt1n ;
    // Copy the smeared parameters to the AOD track
    Double_t x[3];
    Double_t p[3];
    et.GetXYZ(x);
    et.GetPxPyPz(p);
    Double_t cv[21];
    et.GetCovarianceXYZPxPyPz(cv);
    track->SetPosition(x,kFALSE);
    track->SetP(p,kTRUE);
    track->SetCovMatrix(cv);
    // Mark the track as "improved" with a trick (this is done with a trick using layer 7 (ie the 8th))
    UChar_t itsClusterMap = track->GetITSClusterMap();
    SETBIT(itsClusterMap,7);
    track->SetITSClusterMap(itsClusterMap);
}


/*int AliAnalysisTaskHFJetIPQA::GetMCTruth(AliAODTrack * track, int &motherpdg){
    if(!fIsPythia) return 0;
    int pdg = 0;
    AliAODMCParticle *pMCAOD = nullptr;
    if(track->GetLabel()< 0) return pdg;
    pMCAOD = static_cast<AliAODMCParticle*>(fMCArray->At(track->GetLabel()));
    if(!(pMCAOD))  return pdg;
    pdg = pMCAOD->PdgCode();
    motherpdg=0;
    AliAODMCParticle *pMCAODmother = nullptr;
    pMCAODmother = static_cast<AliAODMCParticle*>(fMCArray->At(pMCAOD->GetMother()));
    if(!(pMCAOD))  return pdg;
    motherpdg =pMCAODmother->PdgCode();
    return pdg;
}*/



Bool_t AliAnalysisTaskHFJetIPQA::FillTrackHistograms(AliVTrack *track, double *dca, double *cov, double weight)
{
    FillHist("fh2dTracksImpParXY",GetValImpactParameter(kXY,dca,cov),track->Pt(),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh2dTracksImpParXYZ",GetValImpactParameter(kXYZ,dca,cov),track->Pt(),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh2dTracksImpParZ",dca[1],track->Pt(),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh2dTracksImpParXYSignificance",GetValImpactParameter(kXYSig,dca,cov),track->Pt(),1);     //*this->fXsectionWeightingFactor );
    //FillHist("fh2dTracksImpParXYZSignificance",GetValImpactParameter(kXYZSig,dca,cov),track->Pt(),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh2dTracksImpParZSignificance",GetValImpactParameter(kZSig,dca,cov),track->Pt(),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh1dTracksImpParXY",GetValImpactParameter(kXY,dca,cov),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh1dTracksImpParXYZ",GetValImpactParameter(kXYZ,dca,cov),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh1dTracksImpParXYSignificance",GetValImpactParameter(kXYSig,dca,cov),1);     //*this->fXsectionWeightingFactor );
    //FillHist("fh1dTracksImpParXYZSignificance",GetValImpactParameter(kXYZSig,dca,cov),1);     //*this->fXsectionWeightingFactor );
    //if(fIsPythia){
        //FillHist("fh1dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        //FillHist("fh1dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        //FillHist("fh1dTracksImpParXYSignificance_McCorr",GetValImpactParameter(kXYSig,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        //FillHist("fh1dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        //FillHist("fh2dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),track->Pt(),weight);     //*this->fXsectionWeightingFactor );
        //FillHist("fh2dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),track->Pt(),weight);     //*this->fXsectionWeightingFactor );
        //FillHist("fh2dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),track->Pt(),weight);     //*this->fXsectionWeightingFactor );
    //}
    return kTRUE;
}

/*! \brief Transforms local to global coordinates
 *
 *
 * Transforms local to global coordinates
 */
void AliAnalysisTaskHFJetIPQA::localtoglobal(Double_t alpha ,Double_t* local,Double_t* global)
{
    global[0] = local[0]*  TMath::Sin(alpha) + local[1] * TMath::Cos(alpha);
    global[1] = local[0]*  TMath::Cos(alpha) + local[1] * TMath::Sin(-alpha);
    global[2] = local[2];
    return;
}
/*! \brief Cleanup
 *
 *
 * Cleanup of the event-wise globals
 */
/*void AliAnalysisTaskHFJetIPQA::EventwiseCleanup(){
    fEtaBEvt.clear();
    fPhiBEvt.clear();
    fEtaCEvt.clear();
    fPhiCEvt.clear();
    fEtaUdsgEvt.clear();
    fPhiUdsgEvt.clear();
    fEtaSEvt.clear();
    fPhiSEvt.clear();
    fMcEvtSampled = kFALSE;
}
*/

void AliAnalysisTaskHFJetIPQA::FillRecHistograms(Int_t jetflavour, Double_t recjetpt, Double_t fJetGenPt,Double_t fJetRecEta, Double_t fJetGenEta, Double_t fJetRecPhi, Int_t fUnfoldFracCalc){
  FillHist("fh1dJetRecPt",recjetpt, 1);  //this->fXsectionWeightingFactor );
  FillHist("fh1dJetRecEtaPhiAccepted",fJetRecEta,fJetRecPhi, 1);   //this->fXsectionWeightingFactor );
  //FillHist("fh1dJetRecPtAccepted",recjetpt, 1);  //this->fXsectionWeightingFactor );

  if(fJetGenEta<-99) return;
  if(fIsPythia){
    if(jetflavour==0)     FillHist("fh1dJetRecPtUnidentified",recjetpt, 1);    //this->fXsectionWeightingFactor );
      else if(jetflavour==1)FillHist("fh1dJetRecPtudsg",        recjetpt, 1);    //this->fXsectionWeightingFactor );
      else if(jetflavour==2)FillHist("fh1dJetRecPtc",           recjetpt, 1);    //this->fXsectionWeightingFactor );
      else if(jetflavour==3){
        //
        //Filling matched b jet pt histograms:
        //particle level jets without acceptance cuts, detector level jets with acceptance cuts. This is taken care of by using DefineCutsTaskpp cuts only for detlevel jets!
        if(fUnfoldFracCalc<fUnfoldPseudeDataFrac) {
          if((recjetpt>5)&&(recjetpt<120)){
              //printf("Filling fh1dJetTrueMatchedPtb_PseudoData b hist fNEvent=%i, fUnfoldFracCalc=%i\n", fNEvent, fUnfoldFracCalc);
              FillHist("fh1dJetTrueMatchedPtb_PseudoData",fJetGenPt, 1);
          }
          if(PerformGenLevAcceptanceCuts(fJetGenEta)){
              //printf("Filling fh2dJetGenPtVsJetRecPt_PseudoData b hist fNEvent=%i, fUnfoldFracCalc=%i\n", fNEvent, fUnfoldFracCalc);
              FillHist("fh2dJetGenPtVsJetRecPt_PseudoData",recjetpt, fJetGenPt,1);
              FillHist("fh2dJetGenPtVsJetWideRecPt_PseudoData",recjetpt, fJetGenPt,1);
          }
        }
        else{
           if((recjetpt>5)&&(recjetpt<120)){
              //printf("Filling fh1dJetTrueMatchedPtb_Response b hist fNEvent=%i, fUnfoldFracCalc=%i\n", fNEvent, fUnfoldFracCalc);
              FillHist("fh1dJetTrueMatchedPtb_Response", fJetGenPt, 1);
           }
          if(PerformGenLevAcceptanceCuts(fJetGenEta)){
              //printf("Filling fh2dJetGenPtVsJetRecPt_Response b hist fNEvent=%i, fUnfoldFracCalc=%i\n", fNEvent, fUnfoldFracCalc);
              FillHist("fh2dJetGenPtVsJetRecPt_Response",recjetpt, fJetGenPt,1);
              FillHist("fh2dJetGenPtVsJetWideRecPt_Response",recjetpt, fJetGenPt,1);
          }
        }
      }
      else if(jetflavour==4)FillHist("fh1dJetRecPts", recjetpt, 1);
  }
}

Bool_t AliAnalysisTaskHFJetIPQA::PerformGenLevAcceptanceCuts(Double_t fJetGenEta){
  if(TMath::Abs(fJetGenEta)>(0.9-fJetRadius)){
    return kFALSE;
  }
  return kTRUE;
}

void AliAnalysisTaskHFJetIPQA::FillGenHistograms(Int_t jetflavour,Double_t jetgenpt, Int_t fUnfoldFracCalc){
    FillHist("fh1dJetGenPt",jetgenpt, 1);
    if(jetflavour ==0)      FillHist("fh1dJetGenPtUnidentified",jetgenpt, 1);
    else if(jetflavour ==1) FillHist("fh1dJetGenPtudsg",jetgenpt, 1);
    else if(jetflavour ==2) FillHist("fh1dJetGenPtc",jetgenpt, 1);
    else if(jetflavour ==3){
        if(fUnfoldFracCalc<fUnfoldPseudeDataFrac) {
            //printf("Filling fh1dJetGenPtb_PseudoData b hist fNEvent=%i, fUnfoldFracCalc=%i\n", fNEvent, fUnfoldFracCalc);
            FillHist("fh1dJetGenPtb_PseudoData",jetgenpt, 1);
        }
        else{
            //printf("Filling fh1dJetGenPtb_Response b hist fNEvent=%i, fUnfoldFracCalc=%i\n", fNEvent, fUnfoldFracCalc);
            FillHist("fh1dJetGenPtb_Response",jetgenpt, 1);
        }
    }
    else if(jetflavour ==4) FillHist("fh1dJetGenPts",jetgenpt, 1);
}

double AliAnalysisTaskHFJetIPQA::DoUESubtraction(AliJetContainer* &jetcongen, AliJetContainer* &jetconrec, AliEmcalJet* &jetrec, double jetpt){
    //_________________________
    //Underlying Event Subtraction
    if((!(jetconrec->GetRhoParameter() == nullptr)))
    {
        printf("Correct for Underlying Event.\n");
        jetpt = jetpt - jetconrec->GetRhoVal() * jetrec->Area();
    }
    if(fIsPythia){
        if (jetrec->MatchedJet()) {
            Double_t genpt = jetrec->MatchedJet()->Pt();
            if((!(jetcongen->GetRhoParameter() == nullptr)))
            {
                printf("Correct for Underlying Event.\n");
                genpt = genpt - jetcongen->GetRhoVal() * jetrec->MatchedJet()->Area();
            }
            FillHist("fh2dJetGenPtVsJetRecPt",genpt,jetpt,1);    // this->fXsectionWeightingFactor );
        }
    }
    return jetpt;
}

Bool_t AliAnalysisTaskHFJetIPQA::IsEventAccepted(AliAODEvent *ev){
    if(!fEventCuts.AcceptEvent(ev)){
        return kFALSE;
    }

    //if(!fMCRejectFilter) return true;
    if(!(fIsPythia)) return true; // Only relevant for pt-hard production
    AliDebugStream(1) << "Using custom MC outlier rejection" << std::endl;
    //auto partjets =GetJetContainer("mcparticles");
    AliJetContainer * partjets=static_cast<AliJetContainer*>(fJetCollArray.At(1));
    if(!partjets){
      printf("No particle container found\n");
      return true;
    }

      // Check whether there is at least one particle level jet with pt above n * event pt-hard
      auto jetiter = partjets->accepted();
      auto max = std::max_element(jetiter.begin(), jetiter.end(), [](const AliEmcalJet *lhs, const AliEmcalJet *rhs ) { return lhs->Pt() < rhs->Pt(); });
      if(max != jetiter.end()) {
        // At least one jet found with pt > n * pt-hard
        AliDebugStream(1) << "Found max jet with pt " << (*max)->Pt() << " GeV/c" << std::endl;
        if((*max)->Pt() > fAnalysisCuts[bAnalysisCut_PtHardAndJetPtFactor] * fPtHard){
            //printf("Refuse jet with jetpt=%f, fPtHard=%f, fPtHardAndJetPtFactor=%f\n",(*max)->Pt(), fPtHard,fAnalysisCuts[bAnalysisCut_PtHardAndJetPtFactor]);
            return false;
        }
      }
     return true;
}



void AliAnalysisTaskHFJetIPQA::GetV0Properties(SV0Cand* & sV0, AliAODv0* &v0){
  // particle masses from PDG
  Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
  Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
  sV0->fPt = TMath::Sqrt(v0->Pt2V0()); // transverse momentum of V0

  Double_t dSecVtxPos[3]; // V0 vertex position {x,y,z}
  v0->GetSecondaryVtx(dSecVtxPos);
  Double_t dPrimVtxPos[3]; // primary vertex position {x,y,z}
  fEventVertex->GetXYZ(dPrimVtxPos);
  Double_t dDecayPath[3];
  for(Int_t iPos = 0; iPos < 3; iPos++)
    dDecayPath[iPos] = dSecVtxPos[iPos] - dPrimVtxPos[iPos]; // vector of the V0 path
  //Double_t dDecLen = TMath::Sqrt(dDecayPath[0] * dDecayPath[0] + dDecayPath[1] * dDecayPath[1] + dDecayPath[2] * dDecayPath[2]); // path length L
  Double_t dDecLen2D = TMath::Sqrt(dDecayPath[0] * dDecayPath[0] + dDecayPath[1] * dDecayPath[1]); // transverse path length R
  Double_t dROverPt = dDecLen2D / sV0->fPt; // R/pT

  const AliAODTrack* trackPos = (AliAODTrack*)v0->GetDaughter(0); // positive daughter track
  const AliAODTrack* trackNeg = (AliAODTrack*)v0->GetDaughter(1); // negative daughter track
  if(!trackPos || !trackPos){
      sV0->bDaughsMissing=kTRUE;
  }
  //printf("GetV0Properties: Pointers DaughPos=%p, DaughNeg=%p\n", trackPos, trackNeg);
  //set v0 parameters
  sV0->bOnFly=v0->GetOnFlyStatus();
  sV0->fDCAV0DaughvsDaugh=v0->DcaV0Daughters();
  sV0->fPA=v0->CosPointingAngle(fEventVertex);
  sV0->fDecayRadius=TMath::Sqrt(dSecVtxPos[0] * dSecVtxPos[0] + dSecVtxPos[1] * dSecVtxPos[1]);
  sV0->fLifetimeK0 = dMassPDGK0s * dROverPt; // m*R/pT
  sV0->fLifetimeLambda = dMassPDGLambda * dROverPt; // m*R/pT
  sV0->fEta=v0->Eta();
  sV0->fRapK0=v0->RapK0Short();
  sV0->fRapLambda=v0->RapLambda();
  sV0->fDecayLength3D=TMath::Sqrt(dDecayPath[0] * dDecayPath[0] + dDecayPath[1] * dDecayPath[1] + dDecayPath[2] * dDecayPath[2]);
  sV0->fDecayLength2D=TMath::Sqrt(dDecayPath[0] * dDecayPath[0] + dDecayPath[1] * dDecayPath[1]);
  sV0->fArmenterosAlpha=v0->AlphaV0();
  sV0->fArmenterosPt=v0->PtArmV0();
  sV0->fMassK0=v0->MassK0Short();
  sV0->fMassLambda=v0->MassLambda();
  sV0->fMassAntilambda=v0->MassAntiLambda();
  sV0->fSigmaPosPion   = (fPidResponse ? TMath::Abs(fPidResponse->NumberOfSigmasTPC(trackPos, AliPID::kPion)) : 0.); // difference between measured and expected signal of the dE/dx in the TPC
  sV0->fSigmaPosProton = (fPidResponse ? TMath::Abs(fPidResponse->NumberOfSigmasTPC(trackPos, AliPID::kProton)) : 0.);
  sV0->fSigmaNegPion   = (fPidResponse ? TMath::Abs(fPidResponse->NumberOfSigmasTPC(trackNeg, AliPID::kPion)) : 0.);
  sV0->fSigmaNegProton = (fPidResponse ? TMath::Abs(fPidResponse->NumberOfSigmasTPC(trackNeg, AliPID::kProton)) : 0.);
}

void AliAnalysisTaskHFJetIPQA::GetV0DaughProperties(SV0Daugh* & sTrack,AliAODv0* &v0, bool isPos){
    const AliAODTrack* vTrack =0x0;
    const AliAODVertex* prodVtxDaughter =0x0;
    if(isPos){
        vTrack =(AliAODTrack*)v0->GetDaughter(0); // positive daughter track
    }
    else{
        vTrack =(AliAODTrack*)v0->GetDaughter(1); // positive daughter track
    }
    if(!vTrack) AliError(" There should be daughter tracks but GetV0DaughProperties does not get them!\n");
    prodVtxDaughter = (AliAODVertex*)(vTrack->GetProdVertex()); // production vertex of the positive daughter track
    Int_t cTypeVtxProd = prodVtxDaughter->GetType(); // type of the production vertex

    sTrack->fPt= vTrack->Pt();
    sTrack->fEta=vTrack->Eta();
    sTrack->iCharge=vTrack->Charge();
    sTrack->iCrossedTPC=vTrack->GetTPCClusterInfo(2, 1);
    sTrack->iNoTPCCluster=Double_t(vTrack->GetTPCNclsF());
    if(isPos){
        sTrack->fDCAtoPV=TMath::Abs(v0->DcaPosToPrimVertex());
        //printf("v0: dcapos=%f\n",v0->DcaPosToPrimVertex());
    }
    else{
        sTrack->fDCAtoPV=TMath::Abs(v0->DcaNegToPrimVertex());
        //printf("v0: dcaneg=%f\n",v0->DcaNegToPrimVertex());
    }
    sTrack->bTPCRefitOn=vTrack->IsOn(AliAODTrack::kTPCrefit);
    sTrack->bIsKink=(cTypeVtxProd == AliAODVertex::kKink);
}
/*!
 * IsParticleInCone:
 * decides whether a particle is inside a jet cone
 */
Bool_t AliAnalysisTaskHFJetIPQA::IsParticleInCone(const AliVParticle* part, const AliEmcalJet* jet, Double_t dRMax) {
  if(!part || !jet) AliError(Form("Particle or Jet missing: part=%p, jet=%p\n", part,jet));

  TVector3 vecJet(jet->Px(), jet->Py(), jet->Pz());
  TVector3 vecPart(part->Px(), part->Py(), part->Pz());
  Double_t dR = vecJet.DeltaR(vecPart);
  if(dR <= dRMax){
    //printf("IsParticleInCone:: Accepting as dR=%f, maxDR=%f\n",dR, dRMax);
    return kTRUE;
  }
  //printf("IsParticleInCone:: Rejecting as dR=%f, maxDR=%f\n",dR, dRMax);
  return kFALSE;
}//end trackloop


/*!
 * returns number of V0 daughter within generated jet cone plus the maximal ip of the two
 *
 * - ask whether V0 daughters within cone
 * - obtain the ip
 * - return -999 if daughters not in cone
 */
Int_t AliAnalysisTaskHFJetIPQA::NDaughterInCone(vector<Int_t>& vecDaughLabels, const AliEmcalJet* jet, const AliAODEvent* event, Double_t dRMax, Double_t& ipsig) {
// decides whether a particle is inside a jet cone
  if(!jet) AliError(Form("NDaughterInCone:: Jet missing p=%p\n", jet));

  Int_t iDaughInCone=0;
  Bool_t bDebug=kFALSE;
  Double_t fIPSigMax=-999;
  AliAODMCParticle* pDaugh=0x0;

  for(long unsigned iDaugh=0;iDaugh<vecDaughLabels.size();iDaugh++){
      pDaugh=dynamic_cast<AliAODMCParticle*>(GetMCTrack(vecDaughLabels[iDaugh]));
      if(!pDaugh) printf("NDaughterInCone:: Daugh %lu, label %i not accessible!\n",iDaugh,vecDaughLabels[iDaugh]);

      if(bDebug)printf("NDaughtersInCone:: Checking daugh %lu, label %i\n",iDaugh,vecDaughLabels[iDaugh]);

      Double_t fip=-999;
      //if(IsParticleInCone(pDaugh, jet,dRMax)){
      iDaughInCone++;
      GetMCIP(pDaugh,event, jet, fip);
      if(bDebug)printf("NDaughtersInCone:: daugh %lu in cone: ipsig=%f\n",iDaugh,fip);
      //}
      if(fip>fIPSigMax) fIPSigMax=fip;
      pDaugh=0x0;
  }

  ipsig=fIPSigMax;
  if(bDebug)printf("NDaughtersInCone:: iDaughInCone=%i, ipsig=%f\n",iDaughInCone, ipsig);

  return iDaughInCone;
}//end trackloop

Double_t AliAnalysisTaskHFJetIPQA::GetIPSign(Double_t *XYZatDCA, Double_t* jetp,Double_t* VxVyVz){
  TVector3 jetP3(jetp);
  TVector3 JetDir =jetP3.Unit();
  //printf("GetIPSign:: jetx=%f, jety=%f, jetz=%f\n", JetDir.x(), JetDir.y(),JetDir.z());
  TVector3 D0(XYZatDCA);
  TVector3 vertex(VxVyVz);
  TVector3 DD0(D0.x()-vertex.x(),D0.y()-vertex.y(),0.);   //track impact parameter
  //printf("GetIPSign:: DD0x=%f, DD0y=%f, DD0z=%f\n", DD0.x(), DD0.y(), DD0.z());
  double ps =DD0.Dot(JetDir);
  double value = DD0.Mag()*(ps/fabs(ps));                 //absolut track impact parameter
  double jetsign  = TMath::Sign(1.,value);                       //sign
  //printf("GetIPSign:: ps=%f, value=%f, jetsign=%f\n", ps, value, jetsign);

  return jetsign;
}
/*!
 * Calculates IP for generated MC tracks
 *
 * - generate new AliVertex object from MC header information
 * - generate new AliExternalTrackParam object from AliAODMCParticle Pt, xyz, charge
 * - use PropagateToDCA to calculate the xy dca
 * - returns signed xy dca
 */
Bool_t AliAnalysisTaskHFJetIPQA::GetMCIP(const AliAODMCParticle* track,const AliAODEvent *event, const AliEmcalJet* jetgen, Double_t& ipsig){
    Bool_t bDebug=kFALSE;
    Bool_t success=kFALSE;

    AliAODMCHeader* headerMC = (AliAODMCHeader*)event->FindListObject(AliAODMCHeader::StdBranchName());
    if(!headerMC){
      AliError("No MC header found!");
    }

    // create vertex
    Double_t VxVyVz[3];
    headerMC->GetVertex(VxVyVz);   //checked, passing in following order: x,y,z
    double pos[3] ={VxVyVz[0],VxVyVz[1],VxVyVz[2]};
    AliVertex* vtx = new AliVertex(pos,0,0);  // dispersion,nContributors set to 0
    if(bDebug)vtx->Print();
    if(bDebug)printf("GetMCIP:: Check Vx=%f, Vy=%f, Vz=%f\n", VxVyVz[0],VxVyVz[1],VxVyVz[2]);

    // create track, fetch momentum and position of the AliAODMCParticle...
    double mcmom[3] = {track->Px(),track->Py(),track->Pz()};
    double mcpos[3] = {track->Xv(),track->Yv(),track->Zv()};
    double mccov[21] = {0};              // fake cov.matrix:
    mccov[0]=mccov[2]=mccov[5]=mccov[9]=mccov[14]=mccov[20]=1e-4;
    int charge = track->Charge() ; // get the charge of the particle from pdg code
    AliExternalTrackParam* tr = new AliExternalTrackParam(mcpos,mcmom,mccov, charge);
    if(bDebug)printf("GetMCIP:: Before: Px=%f, Py=%f, Pz=%f\n Xv=%f, Yv=%f, Zv=%f\n Charge=%i\n", mcmom[0],mcmom[1],mcmom[2],mcpos[0],mcpos[1],mcpos[2],charge);
    if(bDebug)printf("GetMCIP:: From track: Px=%f, Py=%f, Pz=%f\n Xv=%f, Yv=%f, Zv=%f\n Charge=%i\n", tr->Px(), tr->Py(), tr->Pz(), tr->Xv(), tr->Yv(), tr->Zv(), tr->Charge());

    //Classic dca calculation
    Double_t XYZatDCA[3];
    Double_t dca[2];
    if(tr->PropagateToDCA(vtx, event->GetMagneticField(), 200,dca)){
        success = kTRUE;
        tr->GetXYZ(XYZatDCA);
        //tr.GetPxPyPz(p_at_dca);
            //         if(fIsPythia)   dca[0] *= fMCglobalDCASmear;
            //         if(fIsPythia)   dca[0] += fMCglobalDCAxyShift; // generic mean offset in LHC10e default is 0.007 == 7 µm

    } else return kFALSE;

    if(bDebug)printf("GetMCIP:: XatDCA=%f, YatDCA=%f, ZatDCA=%f, xydca=%f, zdca=%f\n", XYZatDCA[0],XYZatDCA[1],XYZatDCA[2], dca[0],dca[1]);

    //get IP sign
    double jetp[3];
    jetgen->PxPyPz(jetp);
    Double_t jetsign=GetIPSign(XYZatDCA,jetp,VxVyVz);

    ipsig =TMath::Abs(GetValImpactParameter(kXY,dca,mccov))*jetsign;
    if(bDebug)printf("GetMCIP:: Sign=%f, Ipsig=%f\n",jetsign,ipsig);  //ipsig has to be like xydca!

    delete vtx;
    delete tr;
    return success;
}


//==========================================
void AliAnalysisTaskHFJetIPQA::FillV0Candidates(Bool_t isK, Bool_t isL, Bool_t isAL, Int_t iCut/*cut index*/){
  if(isK){
    fh1V0CounterCentK0s->Fill(iCut);
  }
  if(isL){
    fh1V0CounterCentLambda->Fill(iCut);
  }
  if(isAL){
    fh1V0CounterCentALambda->Fill(iCut);
  }
}


/*!
 * \brief AliAnalysisTaskHFJetIPQA::GetGeneratedV0
 *
 * - loops through fMCArray and stores the properties of all V0 particles within GENERATED jets
 * - cuts on:
 *      *physical primary
 *      *V0 ID
 *      *acceptance
 *      *within jet radius
 *
 * TODO:
 *
 * include also Lambda from cascades
 * Double_t DcaPosToPrimVertex() const;
 * Double_t DcaNegToPrimVertex() const;
 */
/*
void AliAnalysisTaskHFJetIPQA::GetGeneratedV0(){
  AliAODMCParticle* pAOD=NULL;
  AliEmcalJet* jetMC=NULL;
  Double_t fJetPt=-99;
  Bool_t bDebug=kFALSE;
  Int_t iDaughInCone=0;

  //for (Int_t i=0; i<fMCArray->GetEntriesFast(); i++) {
  for(Int_t i=0;i<fMCEvent->GetNumberOfTracks();i++){
    iDaughInCone=0;

    pAOD = dynamic_cast<AliAODMCParticle*>(fMCEvent->GetTrack(i));
    if (!pAOD) continue;

    // Get the distance between the production point of the MC V0 particle and the primary vertex !
    //Double_t dx = dPrimVtxMCX - pAOD->Xv();
    //Double_t dy = dPrimVtxMCY - pAOD->Yv();
    //Double_t dz = dPrimVtxMCZ - pAOD->Zv();
    //Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
    //Bool_t bV0MCIsPrimaryDist = (dDistPrimary < 0.01); // Is close enough to be considered primary-like?
    Bool_t bV0MCIsPrimaryDist = pAOD->IsPhysicalPrimary();
    if(!bV0MCIsPrimaryDist){
      if(bDebug)printf("GetGeneratedV0:: Rejecting as not primari: bV0MCIsPrimaryDist=%i\n",bV0MCIsPrimaryDist);
      continue;
    }

    //Ask for PDG !
    Int_t id = pAOD->GetPdgCode();
    Bool_t bV0 = ((id==3122) || (id==-3122) || (id==310));
    if (!bV0) {
      pAOD=0;
      if(bDebug)printf("GetGeneratedV0:: Rejecting as wrong ID=%i\n",id);
      continue;
    }

    //Daughter Properties
    int itrackPos = pAOD->GetDaughterLabel(0); // positive daughter track
    int itrackNeg = pAOD->GetDaughterLabel(1); // negative daughter track

    AliAODMCParticle* trackPos=(AliAODMCParticle*) GetMCTrack(itrackPos);
    AliAODMCParticle* trackNeg=(AliAODMCParticle*) GetMCTrack(itrackNeg);

    if((!trackPos)||(!trackNeg)) continue;
    Double_t fPosPt=trackPos->Pt();
    Double_t fNegPt=trackNeg->Pt();
    Double_t fPosEta=trackPos->Eta();
    Double_t fNegEta=trackNeg->Eta();

    //acceptance cuts !
    if (TMath::Abs(pAOD->Y()) > fV0Cuts[MaxV0Rap]) { pAOD=0; continue; }
    if (pAOD->Pt() <0.15) continue;
    if ((TMath::Abs(fPosEta) > fV0Cuts[DaughMaxEta])||(TMath::Abs(fNegEta) > fV0Cuts[DaughMaxEta])) continue;
    if ((fPosPt < fV0Cuts[DaughMinPt])||(fNegPt < fV0Cuts[DaughMinPt])) continue;

    if(bDebug)printf("GetGeneratedV0:: Passed acceptance cuts: V0y=%f, V0pt=%f, poseta=%f, negeta=%f, pospt%f, negpt=%f\n", pAOD->Y(), pAOD->Pt(), fPosEta,fNegEta,fPosPt, fNegPt);

    if(!jetcongen) AliError("No MC jet container available!\n");
    jetcongen->ResetCurrentID();

    if(bDebug)printf("GetGeneratedV0:: Having %i jets in container!\n", jetcongen->GetNJets());

    while ((jetMC = jetcongen->GetNextAcceptJet())){
      fJetPt= jetMC->Pt();
      //asking whether v0 within jet cone !
      iDaughInCone=kFALSE;//NDaughterInCone(trackPos, trackNeg, jetMC, fJetRadius);

      double thnentries[4]={static_cast<double>(id), pAOD->Pt(), pAOD->Eta(), fJetPt};
      if(iDaughInCone==0){
          continue;
      }
      if(id==310) {
        for(int iDaugh=0;iDaugh<iDaughInCone;iDaugh++){
          fhnV0InJetK0s->Fill(thnentries);
          if(bDebug)printf("Fount MCTrue K0s: id=%i, eta=%f, pt=%f, jetpt=%f\n",id,pAOD->Eta(),pAOD->Pt(), fJetPt);
        }
      }
      if(id==3122) {
        for(int iDaugh=0;iDaugh<iDaughInCone;iDaugh++){
          fhnV0InJetLambda->Fill(thnentries);
          if(bDebug)printf("Fount MCTrue Lambda: id=%i, eta=%f, pt=%f, jetpt=%f\n",id,pAOD->Eta(),pAOD->Pt(), fJetPt);
        }
      }
      if(id==-3122) {
        for(int iDaugh=0;iDaugh<iDaughInCone;iDaugh++){
          fhnV0InJetALambda->Fill(thnentries);
          if(bDebug)printf("Fount MCTrue ALambda: id=%i, eta=%f, pt=%f, jetpt=%f\n",id,pAOD->Eta(),pAOD->Pt(), fJetPt);
        }
      }
    }//end jet const. loop
  }//end fMCArray

    jetMC=NULL;
    pAOD=NULL;
}
*/

Int_t AliAnalysisTaskHFJetIPQA::FindAllV0Daughters(AliAODMCParticle* pAOD, const AliAODEvent* event, const AliEmcalJet* jetgen, const vector<Int_t>& iTrackLabels, const vector<Double_t>& fTrackRecIPs,Int_t iCount, Int_t iLevel){
    Int_t nDaughters=pAOD->GetNDaughters();
    AliAODMCParticle* pDaugh=0x0;
    Int_t iLabel=-99;
    Bool_t bDebug=kFALSE;
    Double_t fMaxIP=-999;
    Int_t iInVectorInx=-999;
    Int_t iInVectorInxMaxIP=-999;

    Int_t iFirstLabel=pAOD->GetDaughterLabel(0);
    if(bDebug)printf("FindAllV0Daughters:: n=%i, pdgmother=%i, iFirstLabel=%i, \n", nDaughters,pAOD->GetPdgCode(),iFirstLabel);

    for(int iDaugh=0;iDaugh<nDaughters;iDaugh++){
      iLabel=iFirstLabel+iDaugh;
      iInVectorInx=-999;
      pDaugh=dynamic_cast<AliAODMCParticle*>(GetMCTrack(iLabel));
      if(!pDaugh){
        AliError(Form("FindAllV0Daughters:: iDaugh=%i, iCount=%i not available!\n", iDaugh, iCount));
        continue;
       }

      ULong64_t nstatus=pDaugh->GetStatus();
      Int_t nDaughDaughters=pDaugh->GetNDaughters();
      //UInt_t flag=pDaugh->GetFlag();

      std::string spaces;
      spaces.resize(2*iLevel,'  ');
      if(bDebug)printf("%s pdg=%i, iDaugh=%i,  nstatus=%llu, ndaughdaughters=%i, label=%i\n",spaces.c_str(),pDaugh->GetPdgCode(), iDaugh,nstatus,nDaughDaughters    ,iLabel);

      iInVectorInx=IsInVector(iTrackLabels, iLabel,__FUNCTION__);
      if((nDaughDaughters==0)||!(iInVectorInx<0)){
        if(!(iInVectorInx<0)){
          //vecDaughLabels.push_back(iLabel);
          if(fTrackRecIPs[iInVectorInx]>fMaxIP){
            fMaxIP=fTrackRecIPs[iInVectorInx];
            iInVectorInxMaxIP=iInVectorInx;
            if(bDebug)printf("FindAllV0Daughters:: Found final state daughter label =%i!, iInVectorInxMaxIP=%i, IP=%f, IPMax=%f\n",  iLabel,iInVectorInxMaxIP, fTrackRecIPs[iInVectorInx], fMaxIP);
          }//
        }
        else{
          //if(bDebug)printf("FindAllV0Daughters:: Rejecting Daughter as not reconstructedin jet lael =%i!\n", iLabel);
        }
      }
      else{
        if(iCount==100) return -999;
        iCount++;
        //if(bDebug)printf("FindAllV0Daughters:: iDaugh %i not final state, going recursive for the %i'time!\n", iDaugh, iCount);
        FindAllV0Daughters(pDaugh, event, jetgen,iTrackLabels,fTrackRecIPs, iCount, iLevel+1);
      }

      pDaugh=NULL;
    }
    return iInVectorInxMaxIP;
}

/*!
 * Identifies V0 daughter particles and returns maximum ip of V0 daughters
 *
 * - cut on primary particles
 * - cut on K0s and Lamda id
 * - acceptance cuts
 * - return -999 if not V0 candidate
 */
void  AliAnalysisTaskHFJetIPQA::GetGenV0DaughterIP(AliAODMCParticle *pAOD, const AliEmcalJet* jetgen, const AliAODEvent* event, const vector<Int_t>& iTrackLabels, const vector<Double_t>& fTrackRecIPs, Int_t& iInVectorInxMaxIP){
    //Int_t iDaughInCone=0;
    Bool_t bDebug=kFALSE;
    Bool_t bV0MCIsK0s=kFALSE;
    Bool_t bV0MCIsLambda=kFALSE;
    Bool_t bV0MCIsALambda=kFALSE;

    Bool_t bV0MCIsPrimaryDist = pAOD->IsPhysicalPrimary();
    if(!bV0MCIsPrimaryDist){
      return;
    }

    //Ask for PDG !
    Int_t id = pAOD->GetPdgCode();
    Bool_t bV0 = ((id==3122) || (id==-3122) || (id==310));
    if (!bV0) {
      pAOD=0;
      return;
    }

    //Daughter Properties
    int itrackPos = pAOD->GetDaughterLabel(0); // positive daughter track
    int itrackNeg = pAOD->GetDaughterLabel(1); // negative daughter track


    AliAODMCParticle* trackPos=(AliAODMCParticle*) GetMCTrack(itrackPos);
    AliAODMCParticle* trackNeg=(AliAODMCParticle*) GetMCTrack(itrackNeg);


    //MC acceptance cuts
    if(!PerformV0MCAcceptanceCuts(pAOD, trackPos,trackNeg,bV0MCIsK0s,bV0MCIsLambda,bV0MCIsALambda)) return;
    if(bDebug)printf("%s:: bV0MCIsK0s=%i, bV0MCIsLambda=%i, bV0MCIsALambda=%i\n",__FUNCTION__, bV0MCIsK0s, bV0MCIsLambda,bV0MCIsALambda);
    if((!bV0MCIsK0s)&&(!bV0MCIsLambda)&&(!bV0MCIsALambda)) return;
    if(bDebug)printf("%s:: Passed PerformV0MCAcceptanceCuts!\n",__FUNCTION__);

    //acceptance cuts !
    if(!PerformV0AcceptanceCuts(pAOD->Pt(), pAOD->Y(), trackPos->Pt(), trackPos->Eta(),trackNeg->Pt(), trackNeg->Eta())) return;
    if(bDebug)printf("%s:: Passed acceptance cuts",__FUNCTION__);

    //asking whether v0 within jet cone !
    iInVectorInxMaxIP=FindAllV0Daughters(pAOD,event, jetgen,iTrackLabels,fTrackRecIPs,0,0);
    if(bDebug)printf("%s: iInVectorInxMaxIP=%i\n",__FUNCTION__, iInVectorInxMaxIP);
}

Int_t AliAnalysisTaskHFJetIPQA::IsInVector(const vector<Int_t>& vec, Int_t iLabel, TString sFunc){
  Int_t iVecFind=-999;
  for(long unsigned iVec=0;iVec<vec.size();iVec++){
    if(iLabel==vec[iVec]){
      iVecFind=iVec;
      //printf("%s: Found Label for %s: %i at iTrackLabel position %i\n",__FUNCTION__,sFunc.Data(),iLabel, iVecFind);
      return iVecFind;
    }
  }
  return iVecFind;
}

/*!
 * Fills thnsparse with MC truth information about V0 pt, eta, jetpt of generated V0 jet
 *
 * - loops over all MC particles in the event
 * - finds final state particles (MC status=1) within jet cone, calculates their IP
 * - finds V0 particles and calculates the maximum IP of their daughters
 * - stores id, pt and eta of V0 mother if one of their daughter has the maximum IP within the jet
 * - fills thnsparse object if largest IP track within jet is V0 daughter and if jetflavour is not b
 */
Int_t AliAnalysisTaskHFJetIPQA::GetGenV0Jets(const AliEmcalJet* jetgen, const AliAODEvent* event, const std::vector<Int_t>& iTrackLabels, const std::vector<Double_t>& fTrackRecIPs, const std::vector<Double_t>& fTrackRecPts,Int_t fGenJetFlavour, Bool_t **kTagDec, Double_t fLNJP){
    Bool_t bDebug=kFALSE;
    if(fGenJetFlavour==B){
      if(bDebug)printf("GetGenV0Jets:: Returning as B jet jetflavour=%i\n",fGenJetFlavour);
      return kFALSE;
    }
    if(iTrackLabels.size()!=fTrackRecIPs.size()) AliError(Form("%s: size of iTrackLabels (%lu) and fTrackRecIPs (%lu) not matching!\n",__FUNCTION__, iTrackLabels.size(), fTrackRecIPs.size()));
    if(iTrackLabels.size()!=fTrackRecPts.size()) AliError(Form("%s: size of iTrackLabels (%lu) and fTrackRecPts (%lu) not matching!\n",__FUNCTION__, iTrackLabels.size(), fTrackRecIPs.size()));

    AliAODMCParticle* pAOD=NULL;

    Double_t fIPSigPart=-999;
    Int_t iPartIsInJet=-999;
    Bool_t bMaxIPIsV0=kFALSE;
    Int_t iMaxIPTrack=-999;
    Int_t fMaxV0ID=-999;
    Int_t iIsV0Jet=V0Untagged;
    Double_t fMaxIP=-999;
    Double_t fMaxV0Pt=-999;
    Double_t fMaxV0Eta=-999;
    Double_t fJetPt=jetgen->Pt();

    if(bDebug)printf("GetGenV0Jets:: Starting up....................................................\n");

    /*for(long unsigned iTrack=0;iTrack<iTrackLabels.size();iTrack++){

      pAOD = dynamic_cast<AliAODMCParticle*>(fMCEvent->GetMCTrack(iTrackLabels[iTrack]));
      if (!pAOD) continue;

      printf("iTrack=%i, iLabel=%i, pt=%f, rap=%f\n", iTrack, iTrackLabels[iTrack], pAOD->Pt(), pAOD->Eta());
    }*/

    for(Int_t i=0;i<fMCEvent->GetNumberOfTracks();i++){
      pAOD = dynamic_cast<AliAODMCParticle*>(GetMCTrack(i));
      if (!pAOD) continue;
      fIPSigPart=-999;
      iPartIsInJet=-999;

      Int_t id = pAOD->GetPdgCode();
      Bool_t bV0 = ((id==3122) || (id==-3122) || (id==310));

      iPartIsInJet=IsInVector(iTrackLabels,i,__FUNCTION__);
      if((!bV0)&&(iPartIsInJet<0)){
          continue;
      }

      if(!(iPartIsInJet<0)){
        fIPSigPart=fTrackRecIPs[iPartIsInJet];
        //GetMCIP(pAOD,event, jetgen, fIPSigPart);
       // printf("%s:: Got fIPSigPart=%f for track with label %i\n",__FUNCTION__,fIPSigPart, iTrackLabels[iPartIsInJet]);
      }

      if(fIPSigPart>fMaxIP){
        Double_t fOldMaxIP=fMaxIP;
        bMaxIPIsV0=kFALSE;
        iMaxIPTrack=i;
        fMaxIP=fIPSigPart;
        if(bDebug)printf("GetGeneratedV0:: Accepting Part as OldMaxIP=%f, NewMaxIP=%f, bMaxIPIsV0=%i, iMaxIPTrack=%i, bPartIsInJet=%i, pt=%f, rap=%f, \n",fOldMaxIP,fMaxIP, bMaxIPIsV0, iMaxIPTrack,iPartIsInJet, pAOD->Pt(), pAOD->Eta());
        //if(nDaughters!=0)printf("GetGeneratedV0:: Accepting Part with nDaughters=%i as OldMaxIP=%f, NewMaxIP=%f, bMaxIPIsV0=%i, iMaxIPTrack=%i, bPartIsInJet=%i, pt=%f, rap=%f, \n",nDaughters,fOldMaxIP,fMaxIP, bMaxIPIsV0, iMaxIPTrack,iPartIsInJet, pAOD->Pt(), pAOD->Eta());
      }


      //Get V0 IP
      Int_t iInVectorInxMaxIP=-999;
      Double_t fIPV0Daugh=-999;
      Double_t fPtV0Daugh=-999;
      GetGenV0DaughterIP(pAOD, jetgen,  event, iTrackLabels,fTrackRecIPs, iInVectorInxMaxIP);
      if(iInVectorInxMaxIP>=0){
        fIPV0Daugh=fTrackRecIPs[iInVectorInxMaxIP];
        fPtV0Daugh=fTrackRecPts[iInVectorInxMaxIP];
        if(bDebug)printf("%s: Overtaking fIPV0Daugh=%f, fPtV0Daugh=%f\n",__FUNCTION__, fIPV0Daugh, fPtV0Daugh);
      }

      if(fIPV0Daugh>fMaxIP){
        Double_t fOldMaxIP=fMaxIP;
        fMaxIP=fIPV0Daugh;
        bMaxIPIsV0=kTRUE;
        iMaxIPTrack=i;
        fMaxV0Pt=fPtV0Daugh;
        fMaxV0Eta=pAOD->Eta();
        fMaxV0ID=pAOD->GetPdgCode();
        if(bDebug)printf("GetGeneratedV0:: Accepting V0 as OldMaxIP=%f, NewMaxIP=%f, fIPV0=%f, bMaxIPIsV0=%i, iMaxIPTrack=%i, fMaxV0Pt=%f, fMaxV0Eta=%f, fMaxV0ID=%i\n",fOldMaxIP,fMaxIP, fIPV0Daugh,bMaxIPIsV0, iMaxIPTrack,fMaxV0Pt, fMaxV0Eta, fMaxV0ID);
      }
      else{
        //if(bDebug)printf("GetGeneratedV0:: Returning V0 daug as fIPSigV0 %f<fMaxIP %f\n",fIPV0, fMaxIP);
        continue;
      }
    }//end fMCArray

    double thnentries[6]={static_cast<double>(fMaxV0ID), fMaxV0Pt, fLNJP, fJetPt, static_cast<double>(kTagDec[1][Double]),fMaxIP};
    if(bMaxIPIsV0==0){
      if(bDebug)printf("GetGeneratedV0:: Returning as no V0 jet\n");
      return kFALSE;
    }
    if(fMaxV0ID==310) {
        fhnV0InJetK0s->Fill(thnentries); //
        if(bDebug)printf("Found MCTrue K0s: id=%i, eta=%f, pt=%f, jetpt=%f, jetflavour=%i, tagging=%i - %f, fIPV0Max=%f, fLNJP=%f\n",fMaxV0ID, fMaxV0Eta, fMaxV0Pt, fJetPt, fGenJetFlavour,kTagDec[1][Double],static_cast<double>(kTagDec[1][Double]),fMaxIP, fLNJP);
        iIsV0Jet=V0K0s;
    }
    if(fMaxV0ID==3122) {
        fhnV0InJetLambda->Fill(thnentries);
        if(bDebug)printf("Found MCTrue Lambda: id=%i, eta=%f, pt=%f, jetpt=%f,jetflavour=%i, tagging=%i - %f, fIPV0Max=%f, fLNJP=%f\n",fMaxV0ID, fMaxV0Eta, fMaxV0Pt, fJetPt,fGenJetFlavour,kTagDec[1][Double],static_cast<double>(kTagDec[1][Double]),fMaxIP, fLNJP);
        iIsV0Jet=V0Lambda;
    }
    if(fMaxV0ID==-3122) {
        fhnV0InJetALambda->Fill(thnentries);
        if(bDebug)printf("Found MCTrue ALambda: id=%i, eta=%f, pt=%f, jetpt=%f,jetflavour=%i, tagging=%i - %f, fIPV0Max=%f, fLNJP=%f\n",fMaxV0ID, fMaxV0Eta, fMaxV0Pt, fJetPt,fGenJetFlavour,kTagDec[1][Double],static_cast<double>(kTagDec[1][Double]),fMaxIP, fLNJP);
        iIsV0Jet=V0AntiLambda;
    }

      pAOD=NULL;
    if((!bMaxIPIsV0)||((fMaxV0ID!=310)&&(fMaxV0ID!=3122)&&(fMaxV0ID!=-3122))) AliError(Form("%s: Inconsistent decision of V0 jets: MaxIPIsV0=%i, fMaxV0ID=%i\n",__FUNCTION__, bMaxIPIsV0,fMaxV0ID));
    return iIsV0Jet;
}


/*!
 *
 * Reject as MC V0 candidate if:
 * - is not primary (ask for IsPhysicalPrimary)
 * - IDs not matching with K0s or (Anti-)lambda
 * - does not fullfill v0 and v0 daughter rapidity, eta and minpt cuts
 * - cuts all based on generated information?!? Check...
 */

  //AliAODMCHeader* headerMC = (AliAODMCHeader*)fAODIn->FindListObject(AliAODMCHeader::StdBranchName());
  //Double_t dPrimVtxMCX = 0., dPrimVtxMCY = 0., dPrimVtxMCZ = 0.; // position of the MC primary vertex

  //if(!headerMC){
  //  AliError("No MC header found!");
  //}
  // get position of the MC primary vertex
  //dPrimVtxMCX = headerMC->GetVtxX();
  //dPrimVtxMCY = headerMC->GetVtxY();
  //dPrimVtxMCZ = headerMC->GetVtxZ();



  // Get the distance between production point of the MC mother particle and the primary vertex
  //Double_t dx = dPrimVtxMCX - particleMCMother->Xv();
  //Double_t dy = dPrimVtxMCY - particleMCMother->Yv();
  //Double_t dz = dPrimVtxMCZ - particleMCMother->Zv();
  //Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
  //Bool_t bV0MCIsPrimaryDist = (dDistPrimary < fV0Cuts[MinDecayRadius]); // Is close enough to be considered primary-like?



//////////////////////////////////////////////////
/// - In data:
///         - go through V0 candidates and compare track labels to V0 daughter labels.
///         - If they match: tagged as V0 daughter tracks, store reconstructed pt
/// - In MC:
///         - go through V0 candidates and compare track labels, check decay channels, store reconstructed and generated pt
///         - (simulation of survived V0 daughter tracks, mainly need generated pt)
///
Int_t AliAnalysisTaskHFJetIPQA::IsV0Daughter(const AliAODEvent* fAODIn,const AliAODTrack* track, Int_t iTrack){
        fflush(stdout);
        if(!track)return kFALSE;
        AliAODv0* v0aod = 0x0;
        int posid = -1;
        int negid = -1;
        int trackid = -1;
        bool bDebug=kFALSE;

        if(!fV0CandidateArray) {AliError("No fV0CandidateArray available!\n");}

        const AliAODTrack* trackPos =0x0;
        const AliAODTrack* trackNeg =0x0; // negative daughter track
        int tracklabel=track->GetLabel();
        int posdaughlabel=-99;
        int negdaughlabel=-99;
        int iV0Tag=V0No;

        Int_t iV0RecIDTemp=V0Untagged;
        Int_t iV0MCIDTemp=V0Untagged;
        Int_t iV0TagTemp=V0No;
        Bool_t bInvalidMCInfo=kFALSE;

        Bool_t isK0=kTRUE;
        Bool_t isLambda=kTRUE;
        Bool_t isAntiLambda=kTRUE;

        if(bDebug)printf("--------------------------IsV0Daughter: V0 data candidates=%i ------------------------\n",fV0CandidateArray->GetEntriesFast());
        for(int i = 0; i < fV0CandidateArray->GetEntriesFast(); ++i) {
          v0aod = dynamic_cast<AliAODv0*>(fV0CandidateArray->At(i));
          if(!v0aod){
            if(bDebug)printf("v0aod not valid? %p\n", v0aod);
            continue;
          }
          posid=-99;
          negid=-99;
          isK0=kTRUE;
          isLambda=kTRUE;
          isAntiLambda=kTRUE;

          //****************************
          //Data: matching of track labels
          trackPos=(AliAODTrack*) (AliAODTrack*)v0aod->GetDaughter(0); // positive daughter track
          trackNeg=(AliAODTrack*)(AliAODTrack*)v0aod->GetDaughter(1); // positive daughter track
          posid = trackPos->GetID();
          negid = trackNeg->GetID();
          trackid = track->GetID();
          if(bDebug)printf("Investigating Data candidate %i, tracklabel=%i\n",i,tracklabel);

          if(posid == trackid || negid == trackid) {
            //set iV0TagTemp to reconstructed, if no true reconstructed V0 has been found before
            if(iV0TagTemp!=V0TrueRec) iV0TagTemp=V0Rec;//iV0Tag=V0Rec;

            IdentifyRecV0PDG(v0aod->MassK0Short(),v0aod->MassLambda(), v0aod->MassAntiLambda(),isK0, isLambda, isAntiLambda, __FUNCTION__);

            //prohibit second reconstructed V0 with different ID
            if(isK0&&(iV0RecIDTemp!=isK0)&&(iV0RecIDTemp!=V0Untagged)){
              if(bDebug)printf("%s: Found differeing second K0s (isK0=%i, iV0RecIDTemp=%i)",__FUNCTION__, isK0,iV0RecIDTemp);
              return V0No;
            }
            if(isLambda&&(iV0RecIDTemp!=isLambda)&&(iV0RecIDTemp!=V0Untagged)){
              if(bDebug)printf("%s: Found differeing second Lambda (isK0=%i, iV0RecIDTemp=%i)",__FUNCTION__, isLambda,iV0RecIDTemp);
              return V0No;
            }
            if(isAntiLambda&&(iV0RecIDTemp!=isAntiLambda)&&(iV0RecIDTemp!=V0Untagged)){
              if(bDebug)printf("%s: Found differeing second ALambda (isK0=%i, iV0RecIDTemp=%i)",__FUNCTION__, isAntiLambda,iV0RecIDTemp);
              return V0No;
            }

            if(isK0) iV0RecIDTemp=V0K0s; //iV0RecID[iTrack]=V0K0s;
            if(isLambda) iV0RecIDTemp=V0Lambda; //iV0RecID[iTrack]=V0Lambda;
            if(isAntiLambda) iV0RecIDTemp=V0AntiLambda; //iV0RecID[iTrack]=V0AntiLambda;

            if(bDebug)printf("Found V0 candidate %i: posid=%i, negid=%i, trackid=%i, tracketa=%f, V0pt=%f, isK0=%i, isLambda=%i, isAntiLambda=%i, iV0RecID=%i\n", i, posid, negid, trackid,track->Eta(),v0aod->Pt(),isK0, isLambda, isAntiLambda, iV0RecID[iTrack]);
          }
          if(trackid<0)printf("Track ID<0: trackid=%i, tracklabel=%i\n",trackid, tracklabel);

          //Check Matching of Labels not IDs
          posdaughlabel=trackPos->GetLabel();
          negdaughlabel=trackNeg->GetLabel();
          if((i == posdaughlabel)&&(trackid!=posid)&&(fIsPythia)) {
            AliError(Form("Mismatch of trackid=%i and posdaughter id=%i! (tracklabel=%i, daughlabel=%i) Are you using hybrid tracks?",trackid, posid, tracklabel, posdaughlabel));
          }
          if((tracklabel == negdaughlabel)&&(trackid!=negid)&&(fIsPythia)){
            AliError(Form("Mismatch of trackid=%i and negdaughter id=%i! (tracklabel=%i, daughlabel=%i) Are you using hybrid tracks?",trackid, negid, tracklabel, negdaughlabel));
          }

        //****************************
        //test the matching of IDs and compare to MC truth
          if(fIsPythia){
            if(bDebug)printf("Investigating MC candidate %i, fV0ptData=%f, tracklabel=%i\n", i, fV0MotherPt[iTrack], tracklabel);
            if(bInvalidMCInfo){
              continue;
            }
            Int_t iV0MCIDTempTemp=GetV0MCVeto(fAODIn, v0aod,tracklabel);
            if((iV0MCIDTempTemp!=iV0MCIDTemp)&&(iV0MCIDTemp!=V0Untagged)&&(iV0MCIDTempTemp!=V0Untagged)){
              AliError(Form("%s: Setting bInvlaidMCInfo as iV0MCIDTempTemp=%i, iV0MCIDTemp=%i\n",__FUNCTION__, iV0MCIDTempTemp,iV0MCIDTemp));
              bInvalidMCInfo=kTRUE;
              continue;
            }

            if((iV0MCIDTempTemp!=V0Untagged)&&(iV0TagTemp!=V0Rec)){
              iV0TagTemp=V0MC;
              iV0MCIDTemp=iV0MCIDTempTemp;
              if(bDebug)printf("%s: Setting iV0TagTemp to V0MC, iV0MCIDTemp=%i, iV0TagTemp=%i\n", __FUNCTION__, iV0MCIDTemp, iV0TagTemp);
            }
            if((iV0MCIDTempTemp!=V0Untagged)&&(iV0TagTemp==V0Rec)){
              iV0TagTemp=V0TrueRec;
              iV0MCIDTemp=iV0MCIDTempTemp;
              if(bDebug)printf("%s: Setting iV0TagTemp to TrueRec, iV0MCIDTemp=%i, iV0TagTemp=%i\n", __FUNCTION__, iV0MCIDTemp, iV0TagTemp);
            }
          }
        }

        if(iV0MCIDTemp>0){
          iV0MCID[iTrack]=iV0MCIDTemp;
          //printf("%s: Setting iV0MCID to %i\n", __FUNCTION__, iV0MCID[iTrack]);
        }
        if(iV0RecIDTemp>0){
          iV0RecID[iTrack]=iV0RecIDTemp;
          //printf("%s: Setting iV0RecID to %i\n", __FUNCTION__, iV0RecID[iTrack]);
        }
        if(iV0TagTemp>0){
          iV0Tag=iV0TagTemp;
          //printf("%s: Setting iV0Tag to %i\n", __FUNCTION__, iV0Tag);
        }

        if(((iV0MCID[iTrack]>0)||(iV0RecID[iTrack]>0))&&iV0Tag==V0No) throw invalid_argument(Form("%s: MC(%i) and Rec(%i) Ids inconsistent with iV0Tag (%i)!\n",__FUNCTION__, iV0MCID[iTrack], iV0RecID[iTrack], iV0Tag));

        if(bDebug)printf("--------------------------End candidates------------------------\n");
        fflush(stdout);

        return iV0Tag;
}

void AliAnalysisTaskHFJetIPQA::SV0Cand::Print() const{
  printf("----------------- V0  -----------------\n");
  printf("bOnFly=%i, bDaughsMissing=%i\n", bOnFly, bDaughsMissing);
  printf("fDCAV0DaughvsDaugh=%f, fPA=%f, fDecayRadius=%f, fLifetimeK0=%f, fLifetimeLambda=%f\n", fDCAV0DaughvsDaugh,fPA,fDecayRadius,fLifetimeK0,fLifetimeLambda);
  printf("fEta=%f, fPt=%f, fRapK0=%f, fRapLambda=%f, fDecayLength3D=%f, fDecayLength2D=%f\n",fEta, fPt, fRapK0, fRapLambda, fDecayLength3D, fDecayLength3D);
  printf("fArmenterosAlpha=%f, fArmenterosPt=%f, fMassK0=%f, fMassLambda=%f, fMassAntilambda=%f\n",fArmenterosAlpha,fArmenterosPt,fMassK0,fMassLambda,fMassAntilambda);
  printf("fSigmaPosPion=%f fSigmaPosProton=%f fSigmaNegPion=%f fSigmaNegProton=%f\n", fSigmaPosPion,fSigmaPosProton,fSigmaNegPion,fSigmaNegProton);
  printf("bIsCandidateK0s=%i ,bIsCandidateLambda=%i, bIsCandidateALambda=%i, bIsInPeakK0s=%i bIsInPeakLambda=%i bIsInPeakALambda=%i\n", bIsCandidateK0s,bIsCandidateLambda,bIsCandidateALambda,bIsInPeakK0s,bIsInPeakLambda,bIsInPeakALambda);
  printf("bIsInConeJet=%i, bIsInConePerp=%i, bIsInConeRnd=%i, bIsInConeMed=%i, bIsOutsideCones=%i\n",bIsInConeJet,bIsInConePerp,bIsInConeRnd,bIsInConeMed,bIsOutsideCones);
  printf("-----------------------------------------\n");
}

void AliAnalysisTaskHFJetIPQA::SV0Daugh::Print() const{
  printf("----------------- V0 Daugh -----------------\n");
  printf("fPT=%f, fEta=%f, iCharge=%i\n", fPt, fEta, iCharge);
  printf("iCrossedTPC=%i iNoTPCCluster=%i  fDCAtoPV=%f  bTPCRefitOn=%i bIsKink=%i\n", iCrossedTPC, iNoTPCCluster, fDCAtoPV, bTPCRefitOn, bIsKink);
  printf("-----------------------------------------\n");
}

/*AliAODMCParticle* AliAnalysisTaskHFJetIPQA::GetMCTrack( const AliAODTrack* track){
  //
  // return MC track
  //
  if(!fIsPythia) return NULL;
  if(!fMCArray) { AliError("No fMCArray"); return NULL;}
  Int_t nStack = fMCArray->GetEntriesFast();
  Int_t iLabel  = track->GetLabel(); // negative label indicate poor matching quality
  if(iLabel<0||(iLabel > nStack) ){
      AliError(Form("Stupid Label given iLabel=%i\n",iLabel));
      return NULL;
  }
  printf("MCTrack: iLabel=%i\n",iLabel);
  AliAODMCParticle *mctrack =  dynamic_cast<AliAODMCParticle *>(fMCArray->At(iLabel));
  return mctrack;*/
//}

AliAODMCParticle* AliAnalysisTaskHFJetIPQA::GetMCTrack(int iLabel){
  //
  // return MC track
  //
  if(!fIsPythia) return NULL;
  if(!fMCEvent){ AliError("No fMCEvent"); return NULL;}
  Int_t nStack = fMCEvent->GetNumberOfTracks();

  if((iLabel < 0) || (iLabel >= nStack)){
      //printf("Daugh not in array range: iLabel=%i\n", iLabel);
      return NULL;
  }

  AliAODMCParticle *mctrack =  dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(iLabel));

  return mctrack;
}

Bool_t AliAnalysisTaskHFJetIPQA::PerformV0AcceptanceCuts(Double_t V0pt, Double_t V0y, Double_t V0PosDaughpt, Double_t V0PosDaughEta,Double_t V0NegDaughpt, Double_t V0NegDaughEta){
  //printf("%s: Before acceptance cuts! V0pt=%f, V0y=%f, pospt=%f, poseta=%f, negpt=%f, negeta=%f\n",__FUNCTION__, V0pt, V0y, V0PosDaughpt, V0PosDaughEta, V0NegDaughpt, V0NegDaughEta);
  if (TMath::Abs(V0y) > fV0Cuts[MaxV0Rap]) return kFALSE;
  if (V0pt <0.15) return kFALSE;
  if ((TMath::Abs(V0PosDaughEta) > fV0Cuts[DaughMaxEta])||(TMath::Abs(V0NegDaughEta) > fV0Cuts[DaughMaxEta])) return kFALSE;
  if ((V0PosDaughpt < fV0Cuts[DaughMinPt])||(V0NegDaughpt < fV0Cuts[DaughMinPt])) return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskHFJetIPQA::PerformV0MCAcceptanceCuts(const AliAODMCParticle* pAODMother, AliAODMCParticle* pAODPosDaugh, AliAODMCParticle* pAODNegDaugh,Bool_t& bV0MCIsK0s,Bool_t& bV0MCIsLambda,Bool_t& bV0MCIsALambda){
    Bool_t bDebug=kFALSE;

    //Existance of V0 daughters
    if(!pAODPosDaugh || !pAODNegDaugh || !pAODMother){
      if(bDebug)printf("%s: posdaugh=%p, negdaugh=%p, mother=%p not existing",__FUNCTION__, pAODPosDaugh, pAODNegDaugh, pAODMother);
      return kFALSE;
    }

    //Consistency of V0 mothers
    Int_t iIndexMotherPos = pAODPosDaugh->GetMother();
    Int_t iIndexMotherNeg = pAODNegDaugh->GetMother();
    if(iIndexMotherNeg != iIndexMotherPos){
      if(bDebug)printf("%s:: Mothers are different: iIndexMotherNeg=%i, iIndexMotherPos=%i\n", __FUNCTION__,iIndexMotherNeg,iIndexMotherPos);
     return kFALSE;
    }

    //IsPhysicalPrimary
    if(!pAODMother->IsPhysicalPrimary()){
      if(bDebug)printf("%s:: Mother is not physical primary!\n",__FUNCTION__);
      return kFALSE;
    }

    //PDGCode Matching
    Int_t iPdgCodeMother = pAODMother->GetPdgCode();
    Int_t iPdgCodeDaughterPos = pAODPosDaugh->GetPdgCode();
    Int_t iPdgCodeDaughterNeg = pAODNegDaugh->GetPdgCode();


    if((iPdgCodeDaughterPos<0)&&(iPdgCodeDaughterNeg>0)){
      /*printf("%s: Positive dump before:\n",__FUNCTION__);
      pAODPosDaugh->Dump();
      printf("%s: Negative dump before:\n",__FUNCTION__);
      pAODNegDaugh->Dump();*/

      AliAODMCParticle* pDummyPos;
      pDummyPos=pAODPosDaugh;
      pAODPosDaugh=pAODNegDaugh;
      pAODNegDaugh=pDummyPos;
      pDummyPos=NULL;

      iPdgCodeDaughterPos = pAODPosDaugh->GetPdgCode();
      iPdgCodeDaughterNeg = pAODNegDaugh->GetPdgCode();
      /*printf("%s: Positive dump after:\n",__FUNCTION__);
      pAODPosDaugh->Dump();
      printf("%s: Negative dump after:\n",__FUNCTION__);
      pAODNegDaugh->Dump();*/

      if((iPdgCodeDaughterPos<0)&&(iPdgCodeDaughterNeg>0)){
        AliError(Form("%s: Now I completely screwed the signs posID=%i, negID=%i\n",__FUNCTION__,iPdgCodeDaughterPos, iPdgCodeDaughterNeg));
      }
    }

    Int_t iPdgCodePion = 211;
    Int_t iPdgCodeProton = 2212;
    Int_t iPdgCodeK0s = 310;
    Int_t iPdgCodeLambda = 3122;
    bV0MCIsK0s = ((iPdgCodeMother == iPdgCodeK0s) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodePion));
    bV0MCIsLambda = ((iPdgCodeMother == +iPdgCodeLambda) && (iPdgCodeDaughterPos == +iPdgCodeProton) && (iPdgCodeDaughterNeg == -iPdgCodePion));
    bV0MCIsALambda = ((iPdgCodeMother == -iPdgCodeLambda) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodeProton));
    if(bDebug)printf("%s:: iPdgCodeDaughterPos=%i, iPdgCodeDaughterNeg=%i, iPdgCodeMother=%i, bV0MCIsK0s=%i, bV0MCIsLambda=%i, bV0MCIsALambda=%i\n",__FUNCTION__, iPdgCodeDaughterPos, iPdgCodeDaughterNeg,iPdgCodeMother,bV0MCIsK0s, bV0MCIsLambda,bV0MCIsALambda);
    if(bDebug)printf("%s:: iIndexMotherPos=%i, iIndexMotherNeg=%i, isPhysicalPrimary=%i\n",__FUNCTION__,iIndexMotherPos,iIndexMotherNeg,pAODMother->IsPhysicalPrimary());

    return kTRUE;
}

int AliAnalysisTaskHFJetIPQA::GetV0MCVeto(const AliAODEvent* fAODIn, AliAODv0* v0, Int_t tracklabel){
  Bool_t bDebug=kFALSE;
  // PDG codes of used particles

  Bool_t bV0MCIsK0s=kFALSE;
  Bool_t bV0MCIsLambda=kFALSE;
  Bool_t bV0MCIsALambda=kFALSE;

  //***********************
  //Get V0 daughter particle
  const AliAODTrack* postrack = (AliAODTrack*)v0->GetDaughter(0); // positive daughter track
  const AliAODTrack* negtrack = (AliAODTrack*)v0->GetDaughter(1); // positive daughter track
  Int_t iposLabel = postrack->GetLabel();
  Int_t inegLabel = negtrack->GetLabel();
  if((iposLabel!=tracklabel)&&(inegLabel!=tracklabel)){
    //if(bDebug)printf("GetV0MCVeto:: labels of pos. %i and neg. %i not matching track label %i!\n", iposLabel, inegLabel, tracklabel);
    return V0Untagged;
  }
  AliAODMCParticle* particleMCDaughter=0x0;
  Int_t iSearchedLabel=-99;
  if(iposLabel==tracklabel)iSearchedLabel=iposLabel;
  if(inegLabel==tracklabel)iSearchedLabel=inegLabel;
  //if(bDebug) printf("%s: iposLabel=%i, inegLabel=%i, iSearchedLabel=%i, tracklabel=%i\n",__FUNCTION__, iposLabel, inegLabel, iSearchedLabel,tracklabel);
  particleMCDaughter=(AliAODMCParticle*) GetMCTrack(iSearchedLabel);
  if(!particleMCDaughter){
    //if(bDebug)printf("GetV0MCVeto:: daughter particle not existing %p!\n", particleMCDaughter);
    return V0Untagged;
  }

  //************************
  //Get V0 Mother
  Int_t iIndexMother = particleMCDaughter->GetMother();
  AliAODMCParticle* particleMCMother = (AliAODMCParticle*)fMCEvent->GetTrack(iIndexMother);
  if(!particleMCMother){
    //if(bDebug)printf("GetV0MCVeto:: particleMCMother not existing!\n");
    return V0Untagged;
  }

  //************************
  //Get V0 daughter particles
  int itrackPos = particleMCMother->GetDaughterLabel(0); // positive daughter track
  int itrackNeg = particleMCMother->GetDaughterLabel(1); // negative daughter track
  if((itrackPos!=iSearchedLabel)&&(itrackNeg!=iSearchedLabel)){
    if(bDebug)printf("GetV0MCVeto:: neither itrackPos=%i nor itrackNeg=%i equal to iSearchedLabel=%i\n", itrackPos, itrackNeg, iSearchedLabel);
    return V0Untagged;
  }
  AliAODMCParticle* particleMCDaughterPos=(AliAODMCParticle*) GetMCTrack(itrackPos);
  AliAODMCParticle* particleMCDaughterNeg=(AliAODMCParticle*) GetMCTrack(itrackNeg);
  if(!particleMCDaughterNeg || !particleMCDaughterPos){
    if(bDebug)printf("GetV0MCVeto:: negative (%p) or positive (%p) daughter not existing!\n", particleMCDaughterPos, particleMCDaughterNeg);
    return V0Untagged;
  }
  if(bDebug)printf("GetV0MCVeto:: poslabelV0=%i, neglabelV0=%i, motherlabel=%i, poslabelMother=%i, neglabelMother=%i!\n", iposLabel, inegLabel, iIndexMother,itrackPos, itrackNeg);
  //if((itrackPos!=iposLabel)||(itrackNeg!=inegLabel)) printf("%s: Labels not matching!\n");

  //************************
  //Perform Cuts
  if(!PerformV0MCAcceptanceCuts(particleMCMother, particleMCDaughterPos,particleMCDaughterNeg, bV0MCIsK0s, bV0MCIsLambda, bV0MCIsALambda)) return V0Untagged;
  if(bDebug)printf("%s: passed V0MCAcceptanceCuts!\n", __FUNCTION__);
  if(!PerformV0AcceptanceCuts(particleMCMother->Pt(),particleMCMother->Y(),particleMCDaughterPos->Pt(), particleMCDaughterPos->Eta(),particleMCDaughterNeg->Pt(), particleMCDaughterNeg->Eta()))  return V0Untagged;
  if(bDebug)printf("%s: passed acceptance Cuts!\n",__FUNCTION__);

  // K0s
  if(bV0MCIsK0s){ // well reconstructed candidates
      if(bDebug)printf("Found K0s iposlabel=%i, ineglabel=%i, tracklabel=%i,bV0MCIsK0=%i, bV0MCIsLambda=%i, bV0MCIsALambda=%i",
        iposLabel, inegLabel, tracklabel, bV0MCIsK0s, bV0MCIsLambda, bV0MCIsALambda);
      return V0K0s;
  }
  // Lambda
  if(bV0MCIsLambda){ // well reconstructed candidates
      if(bDebug)printf("Found Lambda iposlabel=%i, ineglabel=%i, tracklabel=%i, bV0MCIsK0=%i, bV0MCIsLambda=%i, bV0MCIsALambda=%i",
        iposLabel, inegLabel, tracklabel, bV0MCIsK0s, bV0MCIsLambda, bV0MCIsALambda);
      return V0Lambda;
  }

  // anti-Lambda
  if(bV0MCIsALambda){ // well reconstructed candidates
      if(bDebug)printf("Found ALambda iposlabel=%i, ineglabel=%i, tracklabel=%i, bV0MCIsK0=%i, bV0MCIsLambda=%i, bV0MCIsALambda=%i",
        iposLabel, inegLabel, tracklabel, bV0MCIsK0s, bV0MCIsLambda, bV0MCIsALambda);
    return V0AntiLambda;
  }

  return V0Untagged;
}

void AliAnalysisTaskHFJetIPQA::IdentifyRecV0PDG(Double_t fMassK0, Double_t fMassLambda, Double_t fMassAntiLambda, Bool_t& isK0, Bool_t& IsLambda, Bool_t& IsAntiLambda, TString sIsCalledBy){
    Double_t dMassPeakWindowK0s = fV0Cuts[InvarMassWindowK0]; //0.010; // LF p-p
    Double_t dMassPeakWindowLambda = fV0Cuts[InvarMassWindowLambda]; //0.005; // LF p-p
    Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
    Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
    Bool_t doPrintCuts=kFALSE;
    Bool_t IsInPeakK0=kFALSE;
    Bool_t IsInPeakLambda=kFALSE;
    Bool_t IsInPeakAntiLambda=kFALSE;

    if(doPrintCuts)printf("%s: Is called by %s\n",__FUNCTION__, sIsCalledBy.Data());

    if(TMath::Abs(fMassK0 - dMassPDGK0s) < dMassPeakWindowK0s){
      if(doPrintCuts)printf("Found K0s! mass =%f, window=%f\n",fMassK0,dMassPeakWindowK0s);
      IsInPeakK0 = kTRUE;
    }
    if(TMath::Abs(fMassLambda - dMassPDGLambda) < dMassPeakWindowLambda){
      if(doPrintCuts)printf("Found Lambda! mass =%f, window=%f\n",fMassLambda,dMassPeakWindowLambda);
      IsInPeakLambda = kTRUE;
    }
    if(TMath::Abs(fMassAntiLambda - dMassPDGLambda) < dMassPeakWindowLambda){
      if(doPrintCuts)printf("Found AntiLambda! mass =%f, window=%f\n",fMassAntiLambda,dMassPeakWindowLambda);
      IsInPeakAntiLambda = kTRUE;
    }

    if(!IsInPeakK0&&!IsInPeakLambda&&!IsInPeakAntiLambda){
      isK0 = kFALSE;
      IsLambda = kFALSE;
      IsAntiLambda = kFALSE;
      if(doPrintCuts)printf("Setting all candidates to 0!\n");
    }

    if(IsInPeakK0){
      IsLambda = kFALSE;
      IsAntiLambda = kFALSE;
      if(doPrintCuts)printf("Setting lambdas to 0!\n");
    }

    if(IsInPeakLambda){
      isK0 = kFALSE;
      IsAntiLambda = kFALSE;
      if(doPrintCuts)printf("Setting K0 and AntiLambda to 0!\n");
    }

    if(IsInPeakAntiLambda){
      isK0 = kFALSE;
      IsLambda = kFALSE;
      if(doPrintCuts)printf("Setting K0 and Lambda to 0!\n");
    }

    //printf("%s: isK0=%i, isLambda=%i, isALambda=%i, peakk0=%i, peakLambda=%i, peakALambda=%i\n",__FUNCTION__, isK0, IsLambda, IsAntiLambda, IsInPeakK0,IsInPeakLambda,IsInPeakAntiLambda);
}

void AliAnalysisTaskHFJetIPQA::SelectV0Candidates(const AliAODEvent *fAODIn){
    fflush(stdout);

    AliAODv0* v0 = 0; // pointer to V0 candidates
    SV0Cand* sV0=new SV0Cand();
    SV0Daugh* sPosDaugh=new SV0Daugh();
    SV0Daugh* sNegDaugh=new SV0Daugh();
    Bool_t doPrintCuts=kFALSE;

    fV0CandidateArray->Delete();//Reset the TClonesArray
    if(doPrintCuts)printf("-------------- Entries tcarray=%i\n",fV0CandidateArray->GetEntriesFast());

    // Mean lifetime
    Int_t iNV0s = fAODIn->GetNumberOfV0s(); // get the number of V0 candidates

    //printf("################## Select candidates tcarray=%i ########################\n", fV0CandidateArray->GetEntriesFast());

    for(Int_t iV0 = 0; iV0 < iNV0s; iV0++){
      v0 = fAODIn->GetV0(iV0); // get next candidate from the list in AOD
      if(!v0) continue;

      sV0->Reset();
      sPosDaugh->Reset();
      sNegDaugh->Reset();

      //sV0->Print();
      //sPosDaugh->Print();

      GetV0Properties(sV0,  v0);
      bool hasNoDaughters=sV0->bDaughsMissing;
      if(hasNoDaughters) continue;

      GetV0DaughProperties(sPosDaugh,v0, kTRUE);
      GetV0DaughProperties(sNegDaugh, v0, kFALSE);
      //v0->Print();

      Int_t iCutIndex = 0; // indicator of current selection step
      // 1
      // All V0 candidates
      //printf("Found V0 cand\n");
      FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      iCutIndex++;

      // Start of global cuts
      // 2
      // Reconstruction method
      if(sV0->bOnFly) continue;
      FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      iCutIndex++;

      // 3
      // Tracks TPC OK
      if(sPosDaugh->iCharge == sNegDaugh->iCharge)  continue;// daughters have different charge?
      if(sNegDaugh->iCharge != -1)  continue;// daughters have expected charge?
      if(sPosDaugh->iCharge != 1)  continue;// daughters have expected charge?

      if(fV0Cuts[IsTPCRefitOn]){
        if(doPrintCuts)printf("SelectV0Candidates:: Applying Cut TPC refit on.\n");
        if(!sPosDaugh->bTPCRefitOn) continue;// TPC refit is ON?
        if(!sNegDaugh->bTPCRefitOn)  continue;
      }

      if(fV0Cuts[IsKinkCand]){
        if(doPrintCuts)printf("SelectV0Candidates:: Applying kink cut.\n");
        if(sPosDaugh->bIsKink) continue;// kink daughter rejection
        if(sNegDaugh->bIsKink) continue;
      }

      if(fV0Cuts[DoPosNoTPCClusters]){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring positive TPC clusters!\n");
        if(sPosDaugh->iNoTPCCluster <= 0.)  continue;
        if(sNegDaugh->iNoTPCCluster <= 0.)  continue;
      }

      if(fV0Cuts[MinNoCrossedTPCRows] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring minimum number of TPC Rows %f!\n", fV0Cuts[MinNoCrossedTPCRows]);
        if(sPosDaugh->iCrossedTPC < fV0Cuts[MinNoCrossedTPCRows]) continue;// Crossed TPC padrows
        if(sNegDaugh->iCrossedTPC < fV0Cuts[MinNoCrossedTPCRows]) continue;
      }

      if(fV0Cuts[NoCrossedOverNoTPCClustersMin] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring min crossed rows/tpc clusters=%f\n",fV0Cuts[NoCrossedOverNoTPCClustersMin]);
        if(sPosDaugh->iCrossedTPC / sPosDaugh->iNoTPCCluster < fV0Cuts[NoCrossedOverNoTPCClustersMin]) continue;
        if(sNegDaugh->iCrossedTPC / sNegDaugh->iNoTPCCluster < fV0Cuts[NoCrossedOverNoTPCClustersMin]) continue;
      }

      if(fV0Cuts[NoCrossedOverNoTPCClustersMax] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring max crossed rows/tpc clusters=%f\n",fV0Cuts[NoCrossedOverNoTPCClustersMax]);
        if(sPosDaugh->iCrossedTPC / sPosDaugh->iNoTPCCluster > fV0Cuts[NoCrossedOverNoTPCClustersMax])   continue;
        if(sNegDaugh->iCrossedTPC / sNegDaugh->iNoTPCCluster > fV0Cuts[NoCrossedOverNoTPCClustersMax])   continue;
      }

      FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      iCutIndex++;

      // 4
      // Daughters: transverse momentum cut
      if(fV0Cuts[DaughMinPt] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring daughpt >%f\n",fV0Cuts[DaughMinPt]);
        if((sPosDaugh->fPt < fV0Cuts[DaughMinPt]) || (sNegDaugh->fPt < fV0Cuts[DaughMinPt]))   continue;
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 5
      // Daughters: Impact parameter of daughters to prim vtx
      if(fV0Cuts[MinDCADaughWrtPV] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring daugh dca to pv >%f\n",fV0Cuts[MinDCADaughWrtPV]);
        if((sPosDaugh->fDCAtoPV < fV0Cuts[MinDCADaughWrtPV] ) || (sNegDaugh->fDCAtoPV < fV0Cuts[MinDCADaughWrtPV] ))   continue;
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 6
      // Daughters: DCA
      if(fV0Cuts[MaxDCADaughvsDaugh] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring daugh dca daugh >%f\n",fV0Cuts[MaxDCADaughvsDaugh]);
        if(sV0->fDCAV0DaughvsDaugh > fV0Cuts[MaxDCADaughvsDaugh])   continue;
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 7
      // V0: Cosine of the pointing angle
      if(fV0Cuts[MinCosPAK0] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring K0 pointing angle >%f\n",fV0Cuts[MinCosPAK0]);
        if(sV0->fPA < fV0Cuts[MinCosPAK0]){
          sV0->bIsCandidateK0s = kFALSE;
        }
      }
      if(fV0Cuts[MaxCosPALambda] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring Lambda pointing angle >%f\n",fV0Cuts[MaxCosPALambda]);
        if(sV0->fPA < fV0Cuts[MaxCosPALambda]){
          sV0->bIsCandidateLambda = kFALSE;
          sV0->bIsCandidateALambda = kFALSE;
        }
      }
      if(fV0Cuts[MinCosPAK0] > 0. || fV0Cuts[MaxCosPALambda] > 0.){
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;


      // 8
      // V0: Fiducial volume
      /*if(fdCutRadiusDecayMin > 0. && fdCutRadiusDecayMax > 0.)
      {
        if((dRadiusDecay < fdCutRadiusDecayMin) || (dRadiusDecay > fdCutRadiusDecayMax))
          continue;
            FillV0Candidates(bIsCandidateK0s, bIsCandidateLambda, bIsCandidateALambda, iCutIndex);
      }*/

      if(fV0Cuts[MinDecayRadius]>0){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring min decay radius =%f\n",fV0Cuts[MinDecayRadius]);
        if(sV0->fDecayRadius<fV0Cuts[MinDecayRadius])   {
            continue;
        }
      }
      FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      iCutIndex++;

      // 9
      // Daughters: pseudorapidity cut
      if(fV0Cuts[DaughMaxEta] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring daughter rap <%f\n",fV0Cuts[DaughMaxEta]);
        if((TMath::Abs(sPosDaugh->fEta) > fV0Cuts[DaughMaxEta]) || (TMath::Abs(sNegDaugh->fEta) > fV0Cuts[DaughMaxEta]))
          continue;
            FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 10
      // V0: rapidity cut
      if(fV0Cuts[MaxV0Rap] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring V0 rap <%f\n",fV0Cuts[MaxV0Rap]);
        if(TMath::Abs(sV0->fRapK0) > fV0Cuts[MaxV0Rap])
          sV0->bIsCandidateK0s = kFALSE;
        if(TMath::Abs(sV0->fRapLambda) > fV0Cuts[MaxV0Rap]){
          sV0->bIsCandidateLambda = kFALSE;
          sV0->bIsCandidateALambda = kFALSE;
        }
      }
      FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      iCutIndex++;

      // 11
      // Lifetime cut

      // Mean lifetime
      Double_t dCTauK0s = 2.6844; // [cm] c*tau of K0S
      Double_t dCTauLambda = 7.89; // [cm] c*tau of Lambda
      if(fV0Cuts[MaxLifeTimeK0] > 0.){
        if(doPrintCuts)printf("SelectV0Candidates:: Requiring V0 max lifetime K0<%f, Lambda<%f\n",fV0Cuts[MaxLifeTimeK0],fV0Cuts[MaxLifeTimeLambda]);
        if(sV0->fLifetimeK0 > fV0Cuts[MaxLifeTimeK0] * dCTauK0s) sV0->bIsCandidateK0s = kFALSE;
        if(sV0->fLifetimeLambda > fV0Cuts[MaxLifeTimeLambda] * dCTauLambda)
        {
          sV0->bIsCandidateLambda = kFALSE;
          sV0->bIsCandidateALambda = kFALSE;
        }
      }
      if(fV0Cuts[MaxLifeTimeK0] > 0.){
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 12
      // Daughter PID
      if(fV0Cuts[MaxSigmadEdxTPC] > 0.){
          if(doPrintCuts)printf("SelectV0Candidates:: Requiring max sigma tpc %f \n",fV0Cuts[MaxSigmadEdxTPC]);

          if(sV0->fSigmaPosPion > fV0Cuts[MaxSigmadEdxTPC] || sV0->fSigmaNegPion > fV0Cuts[MaxSigmadEdxTPC]){ // pi+, pi-
            sV0->bIsCandidateK0s = kFALSE;
          }
          if(sV0->fSigmaPosProton > fV0Cuts[MaxSigmadEdxTPC] || sV0->fSigmaNegPion > fV0Cuts[MaxSigmadEdxTPC]){ // p+, pi-
            sV0->bIsCandidateLambda = kFALSE;
          }
          if(sV0->fSigmaNegProton > fV0Cuts[MaxSigmadEdxTPC] || sV0->fSigmaPosPion > fV0Cuts[MaxSigmadEdxTPC]){ // p-, pi+
            sV0->bIsCandidateALambda = kFALSE;
          }
          FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 13
      // Armenteros-Podolanski cut
      if(fV0Cuts[DoArmenteros]){
        if(doPrintCuts)printf("SelectV0Candidates:: Doing AP cut \n");
        if(sV0->fArmenterosPt < TMath::Abs(0.2 * sV0->fArmenterosAlpha)){
          sV0->bIsCandidateK0s = kFALSE;
        }
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 14
      // Invariant mass peak selection
      if(fV0Cuts[DoMassWindow]){
        if(doPrintCuts)printf("SelectCandidates:: Before IsK0=%i, IsLambda=%i, IsALamda=%i, massK0=%f, massLambda=%f, massALambda=%f\n",sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda,sV0->fMassK0, sV0->fMassLambda,sV0->fMassAntilambda);
        IdentifyRecV0PDG(sV0->fMassK0, sV0->fMassLambda,sV0->fMassAntilambda, sV0->bIsCandidateK0s,sV0->bIsCandidateLambda, sV0->bIsCandidateALambda);
        if(doPrintCuts)printf("SelectCandidates:: After IsK0=%i, IsLambda=%i, IsALamda=%i\n",sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda);
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      //if(sV0->bIsCandidateK0s) {fh2dKshortMassVsPt->Fill(sV0->fPt, sV0->fMassK0,1); }
      //if(sV0->bIsCandidateLambda) {fh2dLamdaMassVsPt->Fill(sV0->fPt, sV0->fMassLambda,1); }
      //if(sV0->bIsCandidateALambda) {fh2dAnLamdaMassVsPt->Fill(sV0->fPt, sV0->fMassAntilambda,1); }

      /*if(fIsPythia){
        if(!SelectV0CandidatesMC(fAODIn, v0)) continue;
        if(doPrintCuts)printf("Passed all MC cuts!\n");
      }*/

      //Store V0 candidate for later rejection of daughter tracks
      if(sV0->bIsCandidateK0s || sV0->bIsCandidateLambda || sV0->bIsCandidateALambda){
        int nBins=fV0CandidateArray->GetEntriesFast();
        //v0->Print();
        new((*fV0CandidateArray)[nBins]) AliAODv0(*v0);
      }
    }
    delete sV0;
    delete sPosDaugh;
    delete sNegDaugh;
    //printf("################## Select candidates End Dataarray=%i########################\n", fV0CandidateArray->GetEntriesFast());
    fflush(stdout);

}

void AliAnalysisTaskHFJetIPQA::DefaultInitTreeVars(){
  fJetRecPt=-99;
  fJetFlavour=-99;
  nTracks=0;
  fJetArea=-99;
  fMatchedJetPt=-99;
  fJetProb=-99;
  fJetMass=-99;
  fMeanLNKt=-99;
  fMeanTheta=-99;
  fMeanLNKtSD=-99;
  fMeanThetaSD=-99;
  bMatched=kFALSE;
  bIsTrueGenV0Jet=0;
  std::fill( std::begin( fTrackIPs ), std::end( fTrackIPs ), -99 );
  std::fill( std::begin( fTrackIPSigs ), std::end( fTrackIPSigs ), -99 );
  std::fill( std::begin( fTrackProb ), std::end( fTrackProb ), -99 );
  std::fill( std::begin( fTrackPt ), std::end( fTrackPt ), -99 );
  std::fill( std::begin( fTrackChi2OverNDF ), std::end( fTrackChi2OverNDF ), -99 );
  std::fill( std::begin( fDeltaRij), std::end( fDeltaRij), -99);
  std::fill( std::begin( iTrackITSHits ), std::end( iTrackITSHits ), -99 );
  std::fill( std::begin( iV0MCID ), std::end( iV0MCID ), V0Untagged );
  std::fill( std::begin( iV0RecID ), std::end( iV0RecID ), V0Untagged );
  std::fill( std::begin( bTrackIsV0 ), std::end( bTrackIsV0 ), -99 );
  std::fill( std::begin( fV0MotherPt ), std::end( fV0MotherPt ), -99);
  std::fill( std::begin( fV0MotherPtMC ), std::end( fV0MotherPtMC ), -99);
  std::fill( std::begin( fV0MotherEta ), std::end( fV0MotherEta ), -99);
  std::fill( std::begin( fV0MotherEtaMC ), std::end( fV0MotherEtaMC ), -99);
  std::fill( std::begin( bPassedSD ), std::end( bPassedSD ), kFALSE );
  //std::fill( std::begin( bFull ), std::end( bFull ), kFALSE );
  std::fill( std::begin( bSingle1st ), std::end( bSingle1st ), kFALSE );
  //std::fill( std::begin( bSingle2nd ), std::end( bSingle2nd ), kFALSE );
  //std::fill( std::begin( bSingle3rd ), std::end( bSingle3rd ), kFALSE );
  std::fill( std::begin( bDouble ), std::end( bDouble ), kFALSE );
  //std::fill( std::begin( bTriple ), std::end( bTriple ), kFALSE );
}

/*
 * - Sort after storing track parameters as otherwise assignment wrt. IP and IPSig to different tracks!
 * */

void AliAnalysisTaskHFJetIPQA::DetermineIPVars(std::vector<AliAnalysisTaskHFJetIPQA::SJetIpPati>& sImpParXY, std::vector<AliAnalysisTaskHFJetIPQA::SJetIpPati> sImpParXYSig, vector<Float_t> &ipvalsig, vector<Float_t> &ipval, Int_t &nGoodIPTracks)
{
    nGoodIPTracks=sImpParXY.size();
    if(nGoodIPTracks==0) return;

    std::sort(sImpParXY.begin(),sImpParXY.end(),        AliAnalysisTaskHFJetIPQA::mysort);
    std::sort(sImpParXYSig.begin(),sImpParXYSig.end(),  AliAnalysisTaskHFJetIPQA::mysort);
    //std::sort(sImpParXYZ.begin(),sImpParXYZ.end(),      AliAnalysisTaskHFJetIPQA::mysort);
    //std::sort(sImpParXYZSig.begin(),sImpParXYZSig .end(),AliAnalysisTaskHFJetIPQA::mysort);

    for(int iTrack=0;iTrack<nGoodIPTracks;iTrack++){
      ipvalsig.push_back(sImpParXYSig.at(iTrack).first);
      ipval.push_back(sImpParXY.at(iTrack).first);
      //printf("HasIP0, ipval[%i]=%f, ipvalsig[%i]=%f\n", iTrack, ipval[iTrack], iTrack, ipvalsig[iTrack]);
    }
    if(((int)ipvalsig.size()!=nGoodIPTracks)||((int)ipval.size()!=nGoodIPTracks)) AliError("Size of IP vector not valid!\n");
    //if(hasIPs[0])printf("N=1: cursImParXY=%f, TrackWeight=%f,corridx=%i, pt=%f\n",sImpParXYSig.at(0).first, sImpParXYSig.at(0).second, sImpParXYSig.at(0).trackLabel, sImpParXYSig.at(0).trackpt);
    //if(hasIPs[1])printf("N=2: cursImParXY=%f, TrackWeight=%f, corridx=%i, pt=%f\n",sImpParXYSig.at(1).first, sImpParXYSig.at(1).second, sImpParXYSig.at(1).trackLabel, sImpParXYSig.at(1).trackpt);
    //if(hasIPs[2])printf("N=3: cursImParXY=%f, TrackWeight=%f, corridx=%i, pt=%f\n",sImpParXYSig.at(2).first, sImpParXYSig.at(2).second, sImpParXYSig.at(2).trackLabel, sImpParXYSig.at(2).trackpt);
    //printf("*********************************************************\n");
}

void AliAnalysisTaskHFJetIPQA::PrintAllTreeVars(){
  fflush(stdout);
  printf("-------------------------------------\n");
  printf("Printing all tree vars\n");
  printf("fJetRecPt %f\n",fJetRecPt);
  printf("fJetMass=%f\n",fJetMass);
  printf("fJetFlavour %i\n",fJetFlavour);
  printf("nTracks=%i\n", nTracks);
  printf("fJetArea %f\n",fJetArea);
  printf("fMeanLNKt=%f, fMeanLNKtSD=%f\n",fMeanLNKt,fMeanLNKtSD);
  printf("fMeanTheta=%f, fMeanThetaSD=%f\n",fMeanTheta, fMeanThetaSD);
  printf("fMatchedJetPt %f\n", fMatchedJetPt);
  printf("fJetProb=%f\n",fJetProb);
  printf("bMatched %i\n",bMatched);
  printf("bIsTrueGenV0Jet=%i\n", bIsTrueGenV0Jet);
  printf("Tagging Results:\n");
  /*for(int iThresh=0;iThresh<fNThresholds;iThresh++){
     printf("    Thresh %i,\n     bFull=%i,\n     bSingle1st=%i,\n     bSingle2nd=%i,\n     bSingle3rd=%i,\n     bDouble=%i\n    bTriple=%i\n",
     iThresh, bFull[iThresh],bSingle1st[iThresh],bSingle2nd[iThresh],bSingle3rd[iThresh],bDouble[iThresh],bTriple[iThresh]);
  }*/
  for(int iThresh=0;iThresh<fNThresholds;iThresh++){
     printf("    Thresh %i,\n     bSingle1st=%i,\n     bDouble=%i\n    ",iThresh, bSingle1st[iThresh],bDouble[iThresh]);
  }
  printf("TrackProperties\n");
  for(int iTrack=0;iTrack<nTracks;iTrack++){
    printf("     iTrack=%i,\n     fTrackIP=%f,\n     fTrackIPSigs=%f,\n     fTrackProb=%f,\n     fTrackPt=%f,\n     fTrackChi2OverNdf=%f,\n     fDeltaRij=%f\n,     iTrackITSHits=%i,\n     iV0MCID=%i,\n     iV0RecID=%i,\n     bTrackIsV0=%i,\n     bPassedSD=%i,\n     fV0MotherPt=%f,\n     fV0MotherPtMC=%f,     fV0motherEta=%f\n",
    iTrack, fTrackIPs[iTrack], fTrackIPSigs[iTrack], fTrackProb[iTrack], fTrackPt[iTrack], fTrackChi2OverNDF[iTrack], fDeltaRij[iTrack],iTrackITSHits[iTrack], iV0MCID[iTrack], iV0RecID[iTrack], bTrackIsV0[iTrack], bPassedSD[iTrack], fV0MotherPt[iTrack], fV0MotherPtMC[iTrack], fV0MotherEta[iTrack]);
  }
  printf("-------------------------------------\n");
  fflush(stdout);
}

Bool_t AliAnalysisTaskHFJetIPQA::Run(){
    fNEvent++;
    Int_t fUnfoldFracCalc=fNEvent%100;

    /*if(fUnfoldFracCalc<fUnfoldPseudeDataFrac){ printf("Generating Pseudo Data, fNEvent=%i, fUnfoldFracCalc=%i\n", fNEvent, fUnfoldFracCalc);}
    else{printf("Generating Response matrix, fNEvent=%i, fUnfoldFracCalc=%i\n", fNEvent, fUnfoldFracCalc);}*/

    if(fResponseMode&&fDoTCTagging)AliFatal("Response mode and TC tagging turned on. Both modes cannot be run at the same time\n");
    if(fDoTCTagging!=TCNo&&fDoProbTagging!=ProbNo) AliFatal("Don't do track counting and probability tagging simultaneously!");

    //*******************************
    //Selection
    FillGeneralHistograms();

    AliAODEvent *ev = dynamic_cast<AliAODEvent*>(InputEvent());
    fEventVertex = dynamic_cast<AliAODVertex*>(ev->GetPrimaryVertex());

    if(!fResponseMode)IsEventAccepted(ev);

    fIsEsd =  (InputEvent()->IsA()==AliESDEvent::Class())? kTRUE : kFALSE;
   // EventwiseCleanup();
    if(fIsPythia){
        //if(fIsEsd){
            fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent()) ;
            if (!fMCEvent){
                AliError("Could not retrieve  MC particles! Returning");
                return kFALSE;
            }
        //}
        //else{
        //    fMCArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
        //    if (!fMCArray){
        //        AliError("Could not retrieve AOD MC particles! Returning");
        //        return kFALSE;
        //    }
        //}
      jetcongen = static_cast<AliJetContainer*>(fJetCollArray.At(1));
      jetcongen->ResetCurrentID();
    }
    jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    jetconrec->ResetCurrentID();

    SelectV0Candidates(ev);
    //if(fIsPythia)GetGeneratedV0();

    IncHist("fh1dEventsAcceptedInRun",1);
    FillHist("fh1dNoParticlesPerEvent",InputEvent()->GetNumberOfTracks(),1);

    //If fPythia=kTRUE: only events with suitable pthat make it until here

    //************************
    /*Investigate Tracks:
        * (Smearing & particle weighting could take place here)
        * IP dists of all tracks in event which have a valid IP are filled
    */
    Bool_t HasImpactParameter = kFALSE;
    Double_t dca[2] = {-99999,-99999};
    Double_t cov[3] = {-99999,-99999,-99999};
    Double_t TrackWeight       = 1;
    AliVTrack* trackV = NULL;
    for(long itrack= 0; itrack<InputEvent()->GetNumberOfTracks();++itrack){
        trackV = static_cast<AliVTrack*>(InputEvent()->GetTrack(itrack));
               if(!trackV) {
                   AliInfo("Could not retrieve Track");
                   continue;
        }
        IncHist("fh1dTracksAccepeted",1);
        if(!IsTrackAccepted(trackV,-1)) {
            IncHist("fh1dTracksAccepeted",3);
            continue;
        }
        //if(fRunSmearing)SmearTrack((AliAODTrack*)trackV);
        IncHist("fh1dTracksAccepeted",2);
        FillHist("fh2dAcceptedTracksEtaPhi",trackV->Eta(),trackV->Phi(),1); //this->fXsectionWeightingFactor );
        TrackWeight =1;dca[0]=-9999;dca[1]=-9999;cov[0]=-9999;cov[1]=-9999;cov[2]=-9999;
        HasImpactParameter =kFALSE;
        Double_t xyzatdca[3];
        if (GetImpactParameter(((AliAODTrack*)trackV),(AliAODEvent *)InputEvent(), dca, cov,xyzatdca))HasImpactParameter =kTRUE;
        if(fEventVertex) {
            delete fEventVertex;
            fEventVertex =nullptr;
        }
        if(HasImpactParameter)FillTrackHistograms(trackV,dca,cov,TrackWeight);
        /*if(fIsPythia){
            Int_t corrpartidx =-1;
            double ppt;
            //if(fDoMCCorrection) TrackWeight *= GetMonteCarloCorrectionFactor(trackV,corrpartidx,ppt);
        }*/
    }

    //**********************************
    /* JetMatching between generated and reconstructed Jets
     * - track all reconstructed jets per event
     * - match reconstructed and generated jets
     * - fill histograms for generated hists
     */
    //printf("In Program %f < jetpt <%f, %f < jeteta < %f\n",fAnalysisCuts[bAnalysisCut_MinJetPt],fAnalysisCuts[bAnalysisCut_MaxJetPt],fAnalysisCuts[bAnalysisCut_MinJetEta],fAnalysisCuts[bAnalysisCut_MaxJetEta] );
    if (!jetconrec) return kFALSE;
    FillHist("fh1dNoJetsPerEvent",jetconrec->GetNJets(),1);
    AliEmcalJet * jetgen  = nullptr;
    if(fIsPythia){
      if(!MatchJetsGeometricDefault()) AliInfo("Jet matching did not succeed!");
      jetcongen->ResetCurrentID();
      Int_t jetflavour =0;
      Bool_t is_udgjet = kFALSE;
      while ((jetgen = jetcongen->GetNextAcceptJet())){
        if (!jetgen) continue;
          jetflavour=0; is_udgjet=kFALSE;
          jetflavour =IsMCJetPartonFast(jetgen,fJetRadius,is_udgjet);
          if(PerformGenLevAcceptanceCuts(jetgen->Eta())) FillGenHistograms(jetflavour, jetgen->Pt(),fUnfoldFracCalc);
      }
        jetcongen->ResetCurrentID();
        jetconrec->ResetCurrentID();
    }

    //*************************************************
    // Loop over reconstructed/matched jets for template creation and analysis
    AliEmcalJet * jetrec  = nullptr;
    AliEmcalJet * jetmatched  = nullptr;
    Bool_t is_udgjet = kFALSE;
    Double_t fMatchedJetEta=-999;
    Int_t NJetParticles=0;  //Used for counting particles per jet
    Int_t nGoodIPTracks=-1;
    AliVParticle* vp=0x0;

    std::vector<SJetIpPati> sImpParXY,sImpParXYZ,sImpParXYSig,sImpParXYZSig;
    std::vector<Float_t> ipval;
    std::vector<Float_t> ipvalsig;
    std::vector <Int_t> fJetConstTrackID;
    std::vector <Int_t> iTrackLabels;
    std::vector <Double_t> fTrackRecPt;
    std::vector <Double_t> fTrackRecIPs;

    Int_t isV0=V0No;
    Float_t fIPValue=999;
    Int_t corridx=-1;double ppt;

    jetconrec->ResetCurrentID();
    while ((jetrec = jetconrec->GetNextAcceptJet())){//start jetloop
      if(!jetrec) continue;
      if(jetrec->GetNumberOfTracks()==0) continue;
      is_udgjet=kFALSE;
      fMatchedJetEta=-999;
      DefaultInitTreeVars();
      sImpParXY.clear();
      sImpParXYSig.clear();
      sImpParXYZ.clear();
      sImpParXYZSig.clear();
      ipval.clear();
      ipvalsig.clear();
      fJetConstTrackID.clear();
      fTrackRecIPs.clear();
      iTrackLabels.clear();
      fTrackRecPt.clear();
      vp=0x0;
      NJetParticles=0;
      isV0=V0No;
      nGoodIPTracks=-1;

      fJetRecPt=jetrec->Pt();
      fJetMass=jetrec->M();
      //printf("Generated: fJetRecPt=%f, fJetMass=%f\n",fJetRecPt, fJetMass);
      if(fDoUnderlyingEventSub)fJetRecPt=DoUESubtraction(jetcongen, jetconrec,jetrec, fJetRecPt);
      fJetArea=jetrec->Area();
      //printf("Generated: fJetArea=%f\n", fJetArea);
      //________________________
      //Determination of Jet Flavour
        if(fIsPythia){
          jetmatched = nullptr;
          jetmatched =jetrec->MatchedJet();
          if(jetmatched){
            fJetFlavour = IsMCJetPartonFast(jetmatched,fJetRadius,is_udgjet); //Event based association to save memory
            //printf("Generated: fJetFlavour %i\n",fJetFlavour);
            bMatched=kTRUE;
            //printf("Generated: bMatched %i\n",bMatched);
            fMatchedJetPt=jetmatched->Pt();
            //printf("Generated: bMatchedPt %f\n",fMatchedJetPt);
            fMatchedJetEta=jetmatched->Eta();
          }
          else{
            fJetFlavour=Unid;
          }
        }

        FillRecHistograms(fJetFlavour, fJetRecPt,fMatchedJetPt, jetrec->Eta(),fMatchedJetEta,jetrec->Phi(), fUnfoldFracCalc);
        if(fDoLundPlane)RecursiveParents(jetrec, jetconrec,fJetConstTrackID);


        //_____________________________
        //Determination of impact parameters
        Double_t dca[2] = {-99999,-99999};
        Double_t cov[3] = {-99999,-99999,-99999};
        Double_t sign=0;

        for(UInt_t i = 0; i < jetrec->GetNumberOfTracks(); i++) {//start trackloop
          isV0=V0No;
          Double_t xyzatcda[3];
          dca[0]=-99999; dca[1]=-99999;
          cov[0]=-99999; cov[1]=-99999; cov[2]=-99999;
          sign=0;

          vp = static_cast<AliVParticle*>(jetrec->TrackAt(i, jetconrec->GetParticleContainer()->GetArray()));
          if (!vp){
            AliError("AliVParticle associated to constituent not found");
            continue;
          }
          AliVTrack *vtrack = dynamic_cast<AliVTrack*>(vp);
          if (!vtrack) {
            AliError(Form("Could not receive track%d\n", i));
            continue;
          }
          AliAODTrack *trackV = dynamic_cast<AliAODTrack*>(vtrack);

          if (!trackV || !jetrec)            continue;
          if(!IsTrackAccepted((AliAODTrack*)trackV,fJetFlavour)) continue;
          if(!GetImpactParameterWrtToJet((AliAODTrack*)trackV,(AliAODEvent*)InputEvent(),jetrec,dca,cov,xyzatcda,sign, fJetFlavour)) continue;
          if(fEventVertex) {
           delete fEventVertex;
           fEventVertex =nullptr;
          }

          iTrackLabels.push_back(trackV->GetLabel());

          //printf("dca[0] =%f, dca[1]=%f\n",dca[0], dca[1]);
          isV0=IsV0Daughter(ev,trackV, NJetParticles);

          TrackWeight=1;
          corridx=-1;ppt=-1;
          fIPValue=999;

          fTrackPt[NJetParticles]=trackV->Pt();
          fDeltaRij[NJetParticles]=jetrec->DeltaR(vp);

          for(long unsigned iConst=0;iConst<fJetConstTrackID.size();iConst++){
            //printf("const=%i, userindex=%i, TrackID=%i\n", iConst, fJetConstTrackID[iConst],jetrec->TrackAt(i));
            if(fJetConstTrackID[iConst]==jetrec->TrackAt(i)){
              //printf("Setting SD to true!\n");
              bPassedSD[NJetParticles]=kTRUE;
            }
          }

          //printf("Generated: fTrackPt=%f\n",fTrackPt[NJetParticles]);
          //(fIsPythia&&fDoMCCorrection) ? TrackWeight = GetMonteCarloCorrectionFactor(trackV,corridx,ppt) : TrackWeight =1;
          dca[0]=fabs(dca[0]);

          Double_t cursImParXY     =TMath::Abs(GetValImpactParameter(   kXY,dca,cov))*sign;
          Double_t cursImParXYSig  =TMath::Abs(GetValImpactParameter(kXYSig,dca,cov))*sign;
          Double_t cursImParXYZ    =TMath::Abs(GetValImpactParameter(   kXYZ,dca,cov))*sign;
          Double_t cursImParXYZSig =TMath::Abs(GetValImpactParameter(kXYZSig,dca,cov))*sign;
          //printf("cursImParXY=%f, cursImParXYSig=%f, cursImParXYZ=%f, cursImParXYZSig=%f\n", cursImParXY,cursImParXYSig, cursImParXYZ, cursImParXYZSig);
          fTrackRecIPs.push_back(cursImParXY);
          fTrackRecPt.push_back(fTrackPt[NJetParticles]);
          //printf("%s: Track with label %i, IP=%f, pt=%f\n",__FUNCTION__, trackV->GetLabel(), cursImParXY, fTrackPt[NJetParticles]);

          fTrackIPs[NJetParticles]=cursImParXY;
          fTrackIPSigs[NJetParticles]=cursImParXYSig;
          fTrackChi2OverNDF[NJetParticles]=((AliAODTrack*)vtrack)->Chi2perNDF();
          bTrackIsV0[NJetParticles]=isV0;
          iTrackITSHits[NJetParticles]=(int) vtrack->HasPointOnITSLayer(0) + (int) vtrack->HasPointOnITSLayer(1)+(int) vtrack->HasPointOnITSLayer(2) + (int) vtrack->HasPointOnITSLayer(3) + (int) vtrack->HasPointOnITSLayer(4) + (int) vtrack->HasPointOnITSLayer(5);

          //printf("Generated:fTrackChi2OverNDF=%f,  bTrackIsV0=%i, iTrackITSHits=%i\n",  fTrackChi2OverNDF[NJetParticles],bTrackIsV0[NJetParticles], iTrackITSHits[NJetParticles]);
          SJetIpPati a(cursImParXY, TrackWeight,isV0,kFALSE,corridx,fTrackPt[NJetParticles], iV0MCID[NJetParticles]); sImpParXY.push_back(a);
          SJetIpPati b(cursImParXYZ, TrackWeight,isV0,kFALSE,corridx,fTrackPt[NJetParticles], iV0MCID[NJetParticles]); sImpParXYZ.push_back(b);
          SJetIpPati c(cursImParXYSig, TrackWeight,isV0,kFALSE,corridx,fTrackPt[NJetParticles], iV0MCID[NJetParticles]);sImpParXYSig.push_back(c);
          SJetIpPati d(cursImParXYZSig, TrackWeight,isV0,kFALSE,corridx,fTrackPt[NJetParticles], iV0MCID[NJetParticles]);sImpParXYZSig.push_back(d);
          //printf("curImParXY=%f, isV0=%i, pt=%f\n",sImpParXYSig.back().first,sImpParXYSig.back().is_V0, sImpParXYSig.back().trackpt);


          //for(int iTrack=0;iTrack<NJetParticles;iTrack++){
          //  printf("Normaltrack=%i, fNJetParticles=%i : %f\n",iTrack, NJetParticles, fTrackIPs[iTrack]);
          //  if(TMath::Abs(fTrackIPs[NJetParticles-1]+99.)<0.0001) printf("iTrack=%i, fNJetParticles=%i\n",iTrack, NJetParticles);
          //}
          //if(isV0!=V0Untagged){ printf("Found true reconstructed pt=%f, ip=%f, MCID=%i, DataID=%i\n", fTrackPt[NJetParticles],fTrackIPs[NJetParticles], iV0MCID[NJetParticles], iV0RecID[NJetParticles]);}
          ++NJetParticles;
         }//end trackloop
        nTracks=NJetParticles;
        if(nTracks!=(int)sImpParXY.size()) printf("!!!!!!!!!!!!! nTracks=%i, nJetParticles=%i\n",nTracks, NJetParticles);
        if(nTracks>40){
          AliError(Form("Have nTracks>40: %i -> Exceeding capacity of Tree!\n",nTracks));
          FillHist("fh2dnTracksvsJetPt",nTracks,fJetRecPt,1);
          continue;
        }

        //FillHist("fh1dParticlesPerJet",NJetParticles,1);

        DetermineIPVars(sImpParXY, sImpParXYSig, ipvalsig, ipval, nGoodIPTracks);

        //_____________________________
        //V0 tag decisions
        //printf("isV0=%i, jetflaovur=%i, jetpt=%f", );

        //if(fIsPythia){
        //  if(nGoodIPTracks>0)FillV0EfficiencyHists(sImpParXYSig[0].is_V0, fJetFlavour, fJetRecPt, isV0Jet);
        //}

        //_________________________________________
        //TAGGING
        Bool_t ** kTagDec=new Bool_t*[fNThresholds];
        for(int iThresh=0;iThresh<fNThresholds;iThresh++){
          kTagDec[iThresh]=new Bool_t[6];
          for(int iType=0;iType<6;iType++){
            kTagDec[iThresh][iType]=0;
          }
        }

       if(fDoTCTagging!=TCNo){
           if(fUseSignificance){DoTCTagging(fJetRecPt, nGoodIPTracks,ipvalsig, kTagDec);}
           else{DoTCTagging(fJetRecPt, nGoodIPTracks,ipval, kTagDec);}
       }

       //**************
       //Probability Dists
       fJetProb=GetTrackProbability(fJetRecPt,nGoodIPTracks, ipvalsig);

       Double_t fLNJP=-999;
       if(fJetProb>0)fLNJP=-TMath::Log(fJetProb);
       if((fIsPythia)&&(fJetProb>0))bIsTrueGenV0Jet=GetGenV0Jets(jetrec, ev, iTrackLabels,fTrackRecIPs,fTrackRecPt,fJetFlavour, kTagDec, fLNJP);
       //if(bIsTrueGenV0Jet) printf("%s: Found true V0 jet!\n",__FUNCTION__);

       if(sImpParXY.size()>0){
         if((bIsTrueGenV0Jet!=sImpParXY[0].iv0MCID)&&((sImpParXY[0].is_V0==V0MC)||(sImpParXY[0].is_V0==V0TrueRec))&&(fLNJP>0)&&(fJetFlavour!=B)){
           AliError(Form("Found Inconsistency! bIsTrueGenV0Jet=%i, sImpParIsV0=%i\n", bIsTrueGenV0Jet, sImpParXY[0].is_V0));
           PrintAllTreeVars();
         }
       }


       //PrintAllTreeVars();
       tJetTree->Fill();
       for(int iThresh=0;iThresh<fNThresholds;iThresh++){
         delete kTagDec[iThresh];
       }
       delete [] kTagDec;
     }//end jetloop*/
   return kTRUE;
  }



                       /*Double_t AliAnalysisTaskHFJetIPQA::GetLocalAlphaAOD(AliAODTrack * track)
                        {
                            AliExternalTrackParam etp; etp.CopyFromVTrack(track);
                            return etp.GetAlpha();
                        }
                        Double_t AliAnalysisTaskHFJetIPQA::GetLocalThetaAOD(AliAODTrack * track)
                        {
    // convert to AliExternalTrackParam
                            AliExternalTrackParam etp; etp.CopyFromVTrack(track);
    // propagate

                            Double_t dv[2],dcov[3];
                            const Double_t kBeampiperadius=3;
                            AliVEvent * eev = (AliVEvent*)InputEvent();
                            const  AliVVertex *vtxESDSkip =(const  AliVVertex *) (InputEvent()->GetPrimaryVertex())  ;
                            if(!vtxESDSkip) return -9999;
                            if(!(etp.PropagateToDCA(vtxESDSkip, eev->GetMagneticField(), kBeampiperadius, dv, dcov)))return -9999.;
    // update track position and momentum
                            return etp.Theta();
                        }*/

                        /*Double_t AliAnalysisTaskHFJetIPQA::GetTrackCurvature(AliAODTrack * track)
                        {
    // convert to AliExternalTrackParam
                            AliVEvent * eev = (AliVEvent*)InputEvent();

                            AliExternalTrackParam etp; etp.CopyFromVTrack(track);

                            return etp.GetC(eev->GetMagneticField());
                        }*/




void AliAnalysisTaskHFJetIPQA::UserCreateOutputObjects(){
  Printf("Analysing Jets with Radius: R=%f\n",fJetRadius);

  TString BJetCuts[29] = {
    "#sigma_{Dia}",  //0
    "#sigma_{z}",       //1
    "#sigma_{y}",       //2
    "RelError_{z}",     //3
    "RelError_{y}",     //4
    "N_{Cont}",       //5
    "MaxVtxZ",        //6
    "d_{JetTrack}",   //7
    "d_{z}",          //8
    "d_{xy}",         //9
    "DecayLength",    //10
    "d_{chi2/ndf}", //11
    "p_{T,Track}^{min}",      //12
    "p_{T,TrackMC}^{min}",    //13
    "MinTPCClus",     //14
    "n_{ITS Hits}",   //15
    "chi2_{track}", //16
    "p_{T,Jet}^{min}",//17
    "p_{T,Jet}^{max}",//18
    "#eta_{Jet}^{min}",       //19
    "#eta_{Jet}^{max}",        //20
    "SPD Hits", //21
    "SDD Hits",//22
    "SSD Hits",//23
    "Kink",//24
    "TPC Refit",//25
    "ITS Refit", //26
    "PtHard", //27
    "MinNewVtx" //28
  };

  TString sV0Cuts[18] = {
    "all"/*0*/,
    "reconstr."/*1*/,
    "tracks TPC"/*2*/,
    "Daugh pt"/*3*/,
    "DCA(Daugh<->PV)"/*4*/,
    "DCA(Daugh1<->Daugh2)"/*5*/,
    "cos(PA)"/*6*/,
    "R Decay"/*7*/,
    "Daugh Eta"/*8*/,
    "V0 rap, "/*9*/,
    "lifetime"/*10*/,
    "sigma dEdx"/*11*/,
    "pt Arm"/*12*/,
    "mass wind"/*13*/,
    "in jet"/*14*/,
    "a"/*15*/,
    "b"/*16*/,
    "c"/*17*/,
  };

  TString sTemp[7]={"Unidentified","udsg","c","b","udsgV0","cV0",""};
  for(int iS=0;iS<7;iS++){
      sTemplateFlavour.push_back(sTemp[iS]);
  }

  fIsMixSignalReady_n1=kFALSE;
  fIsMixSignalReady_n2=kFALSE;
  fIsMixSignalReady_n3=kFALSE;
  fn1_mix =-1;
  fn2_mix =-1;
  fn3_mix =-1;
  fIsSameEvent_n1=kFALSE;
  fIsSameEvent_n2=kFALSE;
  fIsSameEvent_n3=kFALSE;

  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) AliError("Analysis manager not found!");
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) AliError("Event handler not found!");
  if (evhand->InheritsFrom("AliESDInputHandler"))  fIsEsd = kTRUE;
  else fIsEsd = kFALSE;
  OpenFile(1);

  Double_t lowIPxy =-1.;  //ranges of xy axis of fh2dTracksImpParXY and fh2dTracksImpParZ
  Double_t highIPxy =1.;

  OpenFile(1);
  if(fIsPythia){

    if(!fPidResponse){
        AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
        AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
        fPidResponse = inputHandler->GetPIDResponse();
    }
    if (!fPidResponse) {
      AliFatal("NULL PID response");
    }
    if(!fCombined) fCombined = new AliPIDCombined();
  }

  //AddHistogramm("fh1dTracksAccepeted","# tracks before/after cuts;;",3,0,3);

  //****************************************
  //QA Plots
  fh1dCutsPrinted =(TH1D*)AddHistogramm("fh1dCutsPrinted","",3,0,3);

  fh1dTracksAccepeted =(TH1D*)AddHistogramm("fh1dTracksAccepeted","# tracks before/after cuts;;",3,0,3);
  fh1dTracksAccepeted->GetXaxis()->SetBinLabel(1,"total");
  fh1dTracksAccepeted->GetXaxis()->SetBinLabel(2,"accepted");
  fh1dTracksAccepeted->GetXaxis()->SetBinLabel(3,"rejected");
  //Event Properties
  fHistManager.CreateTH1("fh1dNoParticlesPerEvent","fh1dNoParticlesvsEvent;#;No Particles/Event",5000, 0, 5000,"s");
  fHistManager.CreateTH1("fh1dNoJetsPerEvent","fh1dNoJetsPerEvent;#;No Jets/Event",400, 0, 100,"s");
  fHistManager.CreateTH1("fh1dEventsAcceptedInRun","fh1dEventsAcceptedInRun;Events Accepted;count",1,0,1,"s");
  fHistManager.CreateTH1("fh1dPtHardMonitor","fh1dPtHardMonitor;ptHard;",500,0,250,"s");
  //Jet Properties
  fHistManager.CreateTH2("fh1dJetRecEtaPhiAccepted","detector level jet;#eta;phi",1,-0.5,0.5,1,0.,TMath::TwoPi(),"s");
  fHistManager.CreateTH2("fh2dAcceptedTracksEtaPhi","accepted tracks;#eta;phi",200,-0.9,0.9,200,0.,TMath::TwoPi(),"s");
  fHistManager.CreateTH1("fh1dJetRecPt","detector level jets;pt (GeV/c); count",500,0,250,"s");
  //fHistManager.CreateTH1("fh1dJetRecPtAccepted","accepted detector level jets;pt (GeV/c); count",500,0,250,"s");
  //fHistManager.CreateTH1("fh1dJetArea","fh1dJetArea;# Jet Area",100,0,1,"s");
  //fHistManager.CreateTH1("fh1dParticlesPerJet","fh1dParticlesPerJet;#, Particles/Jet",100,0,100,"s");
  //fHistManager.CreateTH2("fh2dNoAcceptedTracksvsJetArea","fh2dNoAcceptedTracksvsJetArea;No Accepted Tracks;JetArea",20,0,20,100,0,1);
  //MC properties
  if(fIsPythia){
    fHistManager.CreateTH1("fh1dJetGenPt","generator level jets;pt (GeV/c); count",250,0,250,"s");
    fHistManager.CreateTH1("fh1dJetGenPtUnidentified","generator level jets (no flavour assigned);pt (GeV/c); count",250,0,250,"s");
    fHistManager.CreateTH1("fh1dJetGenPtudsg","generator level udsg jets;pt (GeV/c); count",250,0,250,"s");
    fHistManager.CreateTH1("fh1dJetGenPtc","generator level c jets;pt (GeV/c); count",250,0,250,"s");
    fHistManager.CreateTH1("fh1dJetGenPts","generator level s jets;pt (GeV/c); count",250,0,250,"s");

    fHistManager.CreateTH2("fh2dJetGenPtVsJetRecPt_PseudoData","detector momentum response;red pt;gen pt",115,5,120,200,0,200,"s");
    fHistManager.CreateTH2("fh2dJetGenPtVsJetRecPt_Response","detector momentum response;rec pt;gen pt",115,5,120,200,0,200,"s");
    fHistManager.CreateTH2("fh2dJetGenPtVsJetWideRecPt_PseudoData","detector momentum response;rec pt;gen pt",300,0,300,200,0,200,"s");;
    fHistManager.CreateTH2("fh2dJetGenPtVsJetWideRecPt_Response","detector momentum response;rec pt;gen pt",300,0,300,200,0,200,"s");

    fHistManager.CreateTH1("fh1dJetGenPtb_PseudoData","generator level b jets;pt (GeV/c); count",200,0,200,"s");
    fHistManager.CreateTH1("fh1dJetGenPtb_Response","generator level b jets;pt (GeV/c); count",200,0,200,"s");
    fHistManager.CreateTH1("fh1dJetGenMatchedPtb_PseudoData","generator level b jets;pt (GeV/c); count",200,0,200,"s");
    fHistManager.CreateTH1("fh1dJetGenMatchedPtb_Response","generator level b jets;pt (GeV/c); count",200,0,200,"s");

    fHistManager.CreateTH1("fh1dJetRecPtudsg","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPtUnidentified","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPtc","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPts","detector level jets;pt (GeV/c); count",500,0,250,"s");

    fHistManager.CreateTH1("fh1dJetTrueMatchedPtb_PseudoData","detector level jets;pt (GeV/c); count",200,0,200,"s");
    fHistManager.CreateTH1("fh1dJetTrueMatchedPtb_Response","detector level jets;pt (GeV/c); count",200,0,200,"s");
  }

  fh1DCutInclusive=(TH1D*)AddHistogramm("fh1DCutInclusive","fh1DCutInclusive",30,0,30);
  fh1dCutudg=(TH1D*)AddHistogramm("fh1dCutudg","fh1dCutudg",30,0,30);
  fh1dCutc=(TH1D*)AddHistogramm("fh1dCutc","fh1dCutc",30,0,30);
  fh1dCutb=(TH1D*)AddHistogramm("fh1dCutb","fh1dCutb",30,0,30);
  fh1dCuts=(TH1D*)AddHistogramm("fh1dCuts","fh1dCuts",30,0,30);
  fh1dCuts->GetXaxis()->LabelsOption("v");

  for(Int_t iBin = 0; iBin < 29; iBin++){
          fh1DCutInclusive->GetXaxis()->SetBinLabel(iBin + 1, BJetCuts[iBin].Data());
          if(fIsPythia){
                  fh1dCutudg->GetXaxis()->SetBinLabel(iBin + 1, BJetCuts[iBin].Data());
                  fh1dCutb->GetXaxis()->SetBinLabel(iBin + 1, BJetCuts[iBin].Data());
                  fh1dCutc->GetXaxis()->SetBinLabel(iBin + 1, BJetCuts[iBin].Data());
                  fh1dCuts->GetXaxis()->SetBinLabel(iBin + 1, BJetCuts[iBin].Data());
          }
  }

  //****************************************
  //Lund Plane
  const Int_t dimSpec   = 6;
  const Int_t nBinsSpec[6]     = {50,100,100,20,100,2};
  const Double_t lowBinSpec[6] = {0.,-10,0,0,0,0};
  const Double_t hiBinSpec[6]  = {5.,10.,100,20,100,2};
  fHLundIterative = new THnSparseF("fHLundIterative",
                  "LundIterativePlot [log(1/theta),log(z*theta),pTjet,algo]",
                  dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
  if(fDoLundPlane)fOutput->Add(fHLundIterative);

  //****************************************
  //Histograms for Probability Tagging
  //const char * tagtype[6] = {"Full","Single1st","Single2nd","Single3rd","Double","Triple"};

    //****************************************
    //Track Impact Parameter Distributions
    fHistManager.CreateTH2("fh2dnTracksvsJetPt",";;",200,0,200,240,0,120,"s");
    fHistManager.CreateTH2("fh2dTracksImpParXY","radial imp. parameter ;impact parameter xy (cm);a.u.",2000,lowIPxy,highIPxy,500,0,100.,"s");
    fHistManager.CreateTH2("fh2dTracksImpParXYZ","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-1,1,500,0,100.,"s");
    //fHistManager.CreateTH2("fh2dTracksImpParXYZSignificance","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.,"s");
    fHistManager.CreateTH2("fh2dTracksImpParZ","z imp. parameter ;impact parameter xy (cm);a.u.",2000,lowIPxy,highIPxy,500,0,10.,"s");
    fHistManager.CreateTH2("fh2dTracksImpParXYSignificance","radial imp. parameter sig;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.,"s");
    fHistManager.CreateTH2("fh2dTracksImpParZSignificance","z imp. parameter ;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.,"s");
    fHistManager.CreateTH1("fh1dTracksImpParXY","2d imp. parameter ;impact parameter 2d (cm);a.u.",400,-0.2,0.2,"s");
    fHistManager.CreateTH1("fh1dTracksImpParXYZ","3d imp. parameter ;impact parameter 3d (cm);a.u.",2000,0,1.,"s");
    fHistManager.CreateTH1("fh1dTracksImpParXYSignificance","radial imp. parameter ;impact parameter xy significance;a.u.",200,-30,30,"s");
    //fHistManager.CreateTH1 ("fh1dTracksImpParXYZSignificance","3d imp. parameter ;impact parameter 3d significance;a.u.",2000,0.,100.,"s");

    //****************************************
    //Pt Distributions for N1,N2,N3 Tracks

    //V0Cuts
      //V0s from reconstruction
      const Int_t iNDimInJC = 6;
      Int_t binsKInJC[iNDimInJC] = {200, 200, 200, 200,2,400};
      Double_t xminKInJC[iNDimInJC] = {0, 0., 0., 0.,0,-0.5};
      Double_t xmaxKInJC[iNDimInJC] = {3000, 100., 40., 200.,2,1.5};
      Int_t binsLInJC[iNDimInJC] = {200, 200, 200, 200,2,400};
      Double_t xminLInJC[iNDimInJC] = {0, 0., 0, 0.,0,-0.5};
      Double_t xmaxLInJC[iNDimInJC] = {3000, 100., 40., 200.,2,1.5};

      fh2dKshortMassVsPt=(TH2D*)AddHistogramm("fh2dKshortMassVsPt","KShort Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,0.35, 0.65);
      fh2dLamdaMassVsPt =(TH2D*)AddHistogramm("fh2dLamdaMassVsPt","Lamda Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,1.05,1.25);
      fh2dAnLamdaMassVsPt =(TH2D*)AddHistogramm("fh2dAnLamdaMassVsPt","Anti Lamda Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,1.05,1.25);

      fhnV0InJetK0s = new THnSparseD("fhnV0InJetK0s", ";K0s[ID;#it{p}_{T}^{V0} (GeV/#it{c});-Ln(JP);#it{p}_{T}^{jet} (GeV/#it{c})]; DoubleTag; d_{0} (cm)", iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);
      fhnV0InJetLambda = new THnSparseD("fhnV0InJetLambda", ";Lambda[ID;#it{p}_{T}^{V0} (GeV/#it{c});-Ln(JP);#it{p}_{T}^{jet} (GeV/#it{c})]; DoubleTag; d_{0} (cm)", iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
      fhnV0InJetALambda = new THnSparseD("fhnV0InJetALambda", ";ALambda[ID;#it{p}_{T}^{V0} (GeV/#it{c});-Ln(JP);#it{p}_{T}^{jet} (GeV/#it{c})];DoubleTag; d_{0} (cm)", iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
      fhnV0InJetK0s->Sumw2();
      fhnV0InJetLambda->Sumw2();
      fhnV0InJetALambda->Sumw2();

      fh1V0CounterCentK0s = (TH1D*)AddHistogramm(Form("fh1V0CounterCentK0s_%d", 0), Form("Number of K0s candidates after cuts, cent %s;cut;counts","0-100%"), 18, 0, 18);
      fh1V0CounterCentLambda = (TH1D*)AddHistogramm(Form("fh1V0CounterCentLambda_%d", 0), Form("Number of Lambda candidates after cuts, cent %s;cut;counts", "0-100%"), 18, 0, 18);
      fh1V0CounterCentALambda = (TH1D*)AddHistogramm(Form("fh1V0CounterCentALambda_%d", 0), Form("Number of ALambda candidates after cuts, cent %s;cut;counts", "0-100%"), 18, 0, 18);

      fOutput->Add(fhnV0InJetK0s);
      fOutput->Add(fhnV0InJetLambda);
      fOutput->Add(fhnV0InJetALambda);

      for(Int_t j = 0; j < 18; j++){
            fh1V0CounterCentK0s->GetXaxis()->SetBinLabel(j + 1, sV0Cuts[j].Data());
            fh1V0CounterCentLambda->GetXaxis()->SetBinLabel(j + 1, sV0Cuts[j].Data());
            fh1V0CounterCentALambda->GetXaxis()->SetBinLabel(j + 1, sV0Cuts[j].Data());
      }
      if(fV0CandidateArray == NULL){
          fV0CandidateArray = new TClonesArray("AliAODv0",100);
      }
      fV0CandidateArray->Delete();//Reset the TClonesArray


      //V0s from MC
      fh1dKshortPtMC = new TH1D("fh1dKshortPtMC","KShort Pt MC;p_{T} (GeV/c)",200,0,50);
      fh1dLamdaPtMC = new TH1D("fh1dLamdaPtMC","Lamda Pt MC;p_{T} (GeV/c)",200,0,50);
      fh1dAnLamdaPtMC = new TH1D("fh1dAnLamdaPtMC","Anti Lamda Pt MC;p_{T} (GeV/c)",200,0,50);
      fh1dKshortEtaMC = new TH1D("fh1dKshortEtaMC","KShort Eta MC;p_{T} (GeV/c)",200,-10,10);
      fh1dLamdaEtaMC = new TH1D("fh1dLamdaEtaMC","Lamda Eta MC;p_{T} (GeV/c)",200,-10,10);
      fh1dAnLamdaEtaMC = new TH1D("fh1dAnLamdaEtaMC","Anti Lamda Eta MC;p_{T} (GeV/c)",200,-10,10);

      fh2dKshortPtVsJetPtMC = new TH2D("fh2dKshortPtVsJetPtMC","KShort Pt Vs Jet Pt MC;p_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/#it{c})",200,0,50,200,0, 200);
      fh2dLamdaPtVsJetPtMC = new TH2D("fh2dLamdaPtVsJetPtMC","Lamda Pt Vs Jet Pt MC;p_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/#it{c})",200,0,50,200,0,200);
      fh2dAnLamdaPtVsJetPtMC = new TH2D("fh2dAnLamdaPtVsJetPtMC","Anti Lamda Pt Vs Jet Pt MC;p_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/#it{c})",200,0,50,200,0,200);

      //V0 Tagging efficiency
      h1DV0FalseRec=(TH1D*)AddHistogramm("h1DV0FalseRec",";jetpt (GeV/c); N_{V0}",500,0,250);
      h1DV0TrueDataDef =(TH1D*)AddHistogramm("h1DV0TrueDataDef",";jetpt;N_{V0}",500,0,250);
      h1DV0TrueMCDef =(TH1D*)AddHistogramm("h1DV0TrueMCDef",";jetpt;N_{V0}",500,0,250);
      h1DV0TrueRec =(TH1D*)AddHistogramm("h1DV0TrueRec",";jetpt;N_{V0}",500,0,250);

      fOutput->Add(fh1dKshortPtMC);
      fOutput->Add(fh1dLamdaPtMC);
      fOutput->Add(fh1dAnLamdaPtMC);
      fOutput->Add(fh2dKshortPtVsJetPtMC);
      fOutput->Add(fh2dLamdaPtVsJetPtMC);
      fOutput->Add(fh2dAnLamdaPtVsJetPtMC);

    TIter next(fHistManager.GetListOfHistograms());
    TObject* obj = 0;
    while ((obj = next())) {
      printf("Adding Object %s\n",obj->GetName());
      fOutput->Add(obj);
    }

    //***********************
    //Initialise TTree
    tJetTree = new TTree(Form("tJetTree_R%0.2f_%s",fJetRadius,sTaskName.Data()), Form("tJetTree_R%0.2f_%s",fJetRadius,sTaskName.Data()));
    tJetTree->Branch("fJetRecPt", &fJetRecPt, "fJetRecPt/F");
    tJetTree->Branch("fJetMass",&fJetMass, "fJetMass/F");
    tJetTree->Branch("fJetFlavour", &fJetFlavour, "fJetFlavour/I");
    tJetTree->Branch("nTracks", &nTracks, "nTracks/I");
    tJetTree->Branch("fNEvent",&fNEvent,"fNEvent/I");
    tJetTree->Branch("fNThresholds", &fNThresholds, "fNThresholds/I");
    //tJetTree->Branch("fJetArea", &fJetArea, "fJetArea/F");
    //tJetTree->Branch("fMatchedJetPt", &fMatchedJetPt, "fMatchedJetPt/F");
    tJetTree->Branch("fJetProb",&fJetProb, "fJetProb/F");
    tJetTree->Branch("fMeanLNKt",&fMeanLNKt,"fMeanLNKt/F");
    tJetTree->Branch("fMeanTheta",&fMeanTheta,"fMeanTheta/F");
    tJetTree->Branch("fMeanLNKtSD",&fMeanLNKtSD,"fMeanLNKtSD/F");
    tJetTree->Branch("fMeanThetaSD",&fMeanThetaSD,"fMeanThetaSD/F");
    //tJetTree->Branch("bMatched",&bMatched, "bMatched/O");
    tJetTree->Branch("bIsTrueGenV0Jet",&bIsTrueGenV0Jet, "bIsTrueGenV0Jet/I");
    tJetTree->Branch("fTrackIPs",&fTrackIPs,"fTracksIPs[nTracks]/F");
    tJetTree->Branch("fTrackIPSigs",&fTrackIPSigs,"fTrackIPSigs[nTracks]/F");
    tJetTree->Branch("fTrackProb",&fTrackProb,"fTrackProb[nTracks]/F");
    tJetTree->Branch("fTrackChi2OverNDF",&fTrackChi2OverNDF,"fTrackChi2OverNDF[nTracks]/F");
    tJetTree->Branch("fTrackPt",&fTrackPt,"fTrackP[nTracks]/F");
    tJetTree->Branch("fDeltaRij",&fDeltaRij, "fDeltaRij[nTracks]/F");
    tJetTree->Branch("iTrackITSHits",&iTrackITSHits,"iTrackITSHits[nTracks]/I");
    tJetTree->Branch("iV0MCID",&iV0MCID, "iV0MCID[nTracks]/I");
    tJetTree->Branch("iV0RecID",&iV0RecID, "iV0RecID[nTracks]/I");
    tJetTree->Branch("bTrackIsV0",&bTrackIsV0,"bTrackIsV0[nTracks]/I");
    tJetTree->Branch("fV0MotherPt",&fV0MotherPt, "fV0MotherPt[nTracks]/F");
    tJetTree->Branch("fV0MotherPtMC",&fV0MotherPtMC, "fV0MotherPtMC[nTracks]/F");
    tJetTree->Branch("fV0MotherEta",&fV0MotherEta,"fV0MotherEta[nTracks]/F");
    tJetTree->Branch("fV0MotherEtaMC",&fV0MotherEtaMC,"fV0MotherEtaMC[nTracks]/F");
    tJetTree->Branch("bPassedSD",&bPassedSD,"bPassedSD[nTracks]/O");
    //tJetTree->Branch("bFull",&bFull,"bFull[fNThresholds]/O");
    //tJetTree->Branch("bSingle1st",&bSingle1st,"bSingle1st[fNThresholds]/O");
    //tJetTree->Branch("bSingle2nd",&bSingle2nd,"bSingle2nd[fNThresholds]/O");
    //tJetTree->Branch("bSingle3rd",&bSingle3rd,"bSingle3rd[fNThresholds]/O");
    tJetTree->Branch("bDouble",&bDouble,"bDouble[fNThresholds]/O");
    //tJetTree->Branch("bTriple",&bTriple,"bTriple[fNThresholds]/O");

    PostData(1, fOutput);
    PostData(2, tJetTree);
}

void AliAnalysisTaskHFJetIPQA::UserExecOnce(){
    AliJetContainer *  jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    fAnalysisCuts[bAnalysisCut_MinJetPt]=jetconrec->GetJetPtCut();
    fAnalysisCuts[bAnalysisCut_MaxJetPt]=jetconrec->GetJetPtCutMax();
    fAnalysisCuts[bAnalysisCut_MinJetEta]=jetconrec->GetMinEta();
    fAnalysisCuts[bAnalysisCut_MaxJetEta]=jetconrec->GetMaxEta();

    printf("--------------------------------------------------------------------------------\n");
    printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    printf("XXXXXXXXXX Code version 09.02.21 XXXXXXXXXX\n");
    printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    PrintSettings();
    PrintV0Settings();
}

void AliAnalysisTaskHFJetIPQA::PrintV0Settings(){
    TString v0cuts="";
    Int_t version=1;
    int sForm[24]={1,2,2,0,0,0,0,1,1,0,
        0,1,0,1,0,3,3,0,0,0,3,3,2,2};

    v0cuts=version;
    for(int iV0Cut=0;iV0Cut<24;iV0Cut++){
      v0cuts+="+";
      switch (sForm[iV0Cut]){
          case 0:
            v0cuts+=Form("%0.f",fV0Cuts[iV0Cut]);
            break;

          case 1:
            v0cuts+=Form("%0.1f",fV0Cuts[iV0Cut]);
            break;

          case 2:
            v0cuts+=Form("%0.2f",fV0Cuts[iV0Cut]);
            break;

          case 3:
            v0cuts+=Form("%0.3f",fV0Cuts[iV0Cut]);
            break;
      }
    }

    fh1V0CounterCentK0s->SetTitle(v0cuts.Data());
}

void AliAnalysisTaskHFJetIPQA::PrintSettings(){
    TString jetcuts="";
    TString trackcuts="";
    TString vertexcuts="";
    Int_t version=7;

    jetcuts+=version;
    jetcuts+="+";
    jetcuts+=fAnalysisCuts[bAnalysisCut_MinJetPt];
    jetcuts+="+";
    jetcuts+=fAnalysisCuts[bAnalysisCut_MaxJetPt];
    jetcuts+="+";
    jetcuts+=fAnalysisCuts[bAnalysisCut_MinJetEta];
    jetcuts+="+";
    jetcuts+=fAnalysisCuts[bAnalysisCut_MaxJetEta];
    jetcuts+="+";
    jetcuts+=fAnalysisCuts[bAnalysisCut_MinJetArea];
    jetcuts+="+";
    jetcuts+=fNoJetConstituents;
    jetcuts+="+";
    jetcuts+=fDaughtersRadius;
    jetcuts+="+";
    jetcuts+=Form("%0.1f",fJetRadius);
    jetcuts+="+";
    jetcuts+=Form("%0.1f", fAnalysisCuts[bAnalysisCut_SDz]);
    jetcuts+="+";
    jetcuts+=Form("%0.f", fAnalysisCuts[bAnalysisCut_SDbeta]);
    jetcuts+="+";
    jetcuts+=Form("%0.f", fAnalysisCuts[bAnalysisCut_MaxIPLNJP]);

    printf("Cut Jet Settings: %s\n",jetcuts.Data());

    trackcuts+=version;
    trackcuts+="+";
    trackcuts+=Form("%0.2f",fAnalysisCuts[bAnalysisCut_DCAJetTrack]);
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_MaxDecayLength];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_MinTrackPt];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_MinTrackPtMC];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_MaxDCA_Z];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_MaxDCA_XY];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_MinTPCClus];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_MinITSLayersHit];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_MinTrackChi2];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_HasSPD];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_HasSDD];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_HasSSD];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_KinkCand];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_HasTPCrefit];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_HasITSrefit];

    printf("Cut Track Settings: %s\n",trackcuts.Data());

    vertexcuts+=version;
    vertexcuts+="+";
    vertexcuts+=Form("%0.f",fAnalysisCuts[bAnalysisCut_NContibutors]);
    vertexcuts+="+";
    vertexcuts+=Form("%0.f",fAnalysisCuts[bAnalysisCut_Sigma_Z]);
    vertexcuts+="+";
    vertexcuts+=Form("%0.f",fAnalysisCuts[bAnalysisCut_Sigma_Y]);
    vertexcuts+="+";
    vertexcuts+=Form("%0.f",fAnalysisCuts[bAnalysisCut_RelError_Z]);
    vertexcuts+="+";
    vertexcuts+=Form("%0.f",fAnalysisCuts[bAnalysisCut_RelError_Y]);
    vertexcuts+="+";
    vertexcuts+=fAnalysisCuts[bAnalysisCut_SigmaDiamond];
    vertexcuts+="+";
    vertexcuts+=fAnalysisCuts[bAnalysisCut_PtHardAndJetPtFactor];
    vertexcuts+="+";
    vertexcuts+=fAnalysisCuts[bAnalysisCut_MinNewVertexContrib];
    vertexcuts+="+";
    vertexcuts+=fAnalysisCuts[bAnalysisCut_MaxVtxZ];
    vertexcuts+="+";
    vertexcuts+=fAnalysisCuts[bAnalysisCut_Z_Chi2perNDF];
    vertexcuts+="+";
    vertexcuts+=fVertexRecalcMinPt;
    vertexcuts+="+";
    vertexcuts+=fDoMCCorrection;
    vertexcuts+="+";
    vertexcuts+=fRunSmearing;
    vertexcuts+="+";
    vertexcuts+=fDoUnderlyingEventSub;
    vertexcuts+="+";
    vertexcuts+=fDoFlavourMatching;
    vertexcuts+="+";
    vertexcuts+=fResponseMode;
    vertexcuts+="+";
    vertexcuts+=fDoTCTagging;
    vertexcuts+="+";
    vertexcuts+=fDoProbTagging;
    vertexcuts+="+";
    vertexcuts+=fUseSignificance;
    vertexcuts+="+";
    vertexcuts+=fDoJetProb;
    vertexcuts+="+";
    vertexcuts+=fNThresholds;
    vertexcuts+="+";
    vertexcuts+=fNTrackTypes;
    vertexcuts+="+";
    vertexcuts+=Form("%0.3f",fTCThresholdPtFixed);

    fh1dCutsPrinted->SetTitle(jetcuts.Data());
    fh1dCutsPrinted->GetXaxis()->SetTitle(vertexcuts.Data());
    fh1dCutsPrinted->GetYaxis()->SetTitle(trackcuts.Data());

    printf("Vertex Cuts: %s\n",vertexcuts.Data());
}

//NotInUse
/*void AliAnalysisTaskHFJetIPQA::GetMaxImpactParameterCutR(const AliVTrack * const track, Double_t &maximpactRcut){
    //
    // Get max impact parameter cut r (pt dependent)
    //
    Double_t pt = track->Pt();
    if(pt > 0.15) {
            maximpactRcut = 0.0182 + 0.035/TMath::Power(pt,1.01);  // abs R cut
        }
        else maximpactRcut = 9999999999.0;
    }*/


    bool AliAnalysisTaskHFJetIPQA::GetPIDCombined(AliAODTrack * track, double  * prob, int &nDetectors,UInt_t &usedDet ,AliPID::EParticleType &MostProbablePID, bool setTrackPID ){
        AliPIDResponse::EDetPidStatus status[AliPIDResponse::kNdetectors] = {AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal};
        unsigned int nGoodDet = 0;
        for (int j =0; j<AliPIDResponse::kNdetectors;j++)
        {
            for (int i =AliPID::kElectron; i<AliPID::kSPECIES;i++)
            {
                double val = 0;
                status[j] =  fPidResponse->NumberOfSigmas(static_cast <AliPIDResponse::EDetector>(j), track, static_cast <AliPID::EParticleType>(i), val);
                if (status[j] == AliPIDResponse::kDetPidOk ){
                    nGoodDet++;}
                }
            }
            if( nGoodDet/7 <2 ) return false;
            nDetectors = nGoodDet/7;
            Double_t probTPCTOF[AliPID::kSPECIES]={-1.,-1.,-1.,-1.,-1.};
            fCombined->SetDefaultTPCPriors();
            fCombined->SetDetectorMask(AliPIDResponse::kDetITS|AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF|AliPIDResponse::kDetTRD|AliPIDResponse::kDetEMCAL);
            usedDet  = fCombined->ComputeProbabilities((AliVTrack*)track, fPidResponse   , probTPCTOF);
            int maxpid=14;
            double maxpidv=0;
            for (int j =0 ;j< AliPID::kSPECIES;++j ){
                prob[j] =probTPCTOF[j];
                if(prob[j] >maxpidv ){
                    maxpid=j;
                    maxpidv=prob[j] ;
                }
            }
            MostProbablePID=static_cast <AliPID::EParticleType>(maxpid);
            return true;
        }


        bool AliAnalysisTaskHFJetIPQA::IsFromElectron(AliAODTrack*track){

            Double_t p[5] ={0};
            Int_t nDet=0;
            UInt_t nDetUCom=0;
            AliPID::EParticleType mpPID=AliPID::kUnknown;
            if(GetPIDCombined( track, p, nDet,nDetUCom ,mpPID, kFALSE )){
                if(mpPID==AliPID::kElectron && p[0] >0.90 && nDet >1) return kTRUE;
            }
            return false;
        }

        bool AliAnalysisTaskHFJetIPQA::IsFromPion(AliAODTrack*track){
            Double_t p[5] ={0};
            Int_t nDet=0;
            UInt_t nDetUCom=0;
            AliPID::EParticleType mpPID=AliPID::kUnknown;

            if(GetPIDCombined( track, p, nDet,nDetUCom ,mpPID, kFALSE )){
                if(mpPID==AliPID::kPion && p[2] >0.90 && nDet >1) return kTRUE;
            }
            return false;
        }
        bool AliAnalysisTaskHFJetIPQA::IsFromKaon(AliAODTrack*track){
            Double_t p[5] ={0};
            Int_t nDet=0;
            UInt_t nDetUCom=0;
            AliPID::EParticleType mpPID=AliPID::kUnknown;
            if(GetPIDCombined( track, p, nDet,nDetUCom ,mpPID, kFALSE )){
                if(mpPID==AliPID::kKaon && p[3] >0.90 && nDet >1) return kTRUE;
            }
            return false;
        }
        bool AliAnalysisTaskHFJetIPQA::IsFromProton(AliAODTrack*track){
            Double_t p[5] ={0};
            Int_t nDet=0;
            UInt_t nDetUCom=0;
            AliPID::EParticleType mpPID=AliPID::kUnknown;
            if(GetPIDCombined( track, p, nDet,nDetUCom ,mpPID, kFALSE )){
                if(mpPID==AliPID::kProton && p[3] >0.90 && nDet >1) return kTRUE;
            }
            return false;
        }

        int AliAnalysisTaskHFJetIPQA::DetermineUnsuitableVtxTracks(int *skipped, AliAODEvent * const aod, AliVTrack * const track){
            Int_t nTracs=aod->GetNumberOfTracks();
            AliAODTrack * t = nullptr;
            AliExternalTrackParam etp_at_r39_old; etp_at_r39_old.CopyFromVTrack(track);
            etp_at_r39_old.PropagateTo(3.9,InputEvent()->GetMagneticField());
            double angle0 = TMath::ATan2(etp_at_r39_old.Yv(),etp_at_r39_old.Xv());
            double zz0    = etp_at_r39_old.GetZ();
            int nTrksToSkip=1;

            for(Int_t i=0; i<nTracs; i++){
                t = (AliAODTrack *)(aod->GetTrack(i));
                if(!((((AliAODTrack*)t)->TestFilterBit(4))))continue;
                int id = (Int_t)t->GetID();
                AliExternalTrackParam etp_at_r39; etp_at_r39.CopyFromVTrack(t);
                etp_at_r39.PropagateTo(3.9,InputEvent()->GetMagneticField());
                double angle = TMath::ATan2(etp_at_r39.Yv(),etp_at_r39.Xv());
                double zz    = etp_at_r39.GetZ();
                bool doskip=false;
                if(t->Pt()<fVertexRecalcMinPt)                      doskip=true;
                if(fabs(TVector2::Phi_mpi_pi(angle-angle0))>TMath::Pi()/6.) {
                    doskip=true;
                }
                if(fabs(zz-zz0)>0.5) {
                    doskip=true;
                }
                if(doskip && !fGlobalVertex){
                    skipped[nTrksToSkip++] = id;
                }
            }
            return nTrksToSkip;
        }
        /*!
 * \brief AliAnalysisTaskHFJetIPQA::RemoveDaughtersFromPrimaryVtx
 * - discard vertex with title vertexer tracks
 * - perform cuts on number of vertex contributors + chi2 of vertex
 *
 */


AliAODVertex *AliAnalysisTaskHFJetIPQA::RemoveDaughtersFromPrimaryVtx( const AliVTrack * const track) {
   //Initialisation of vertexer
   const AliAODEvent * aod =  ((AliAODEvent*)InputEvent());
   AliAODVertex *vtxAOD =aod ->GetPrimaryVertex();
   if(!vtxAOD) return 0;
   //printf("Before remove:\n");
   //UvtxAOD->Print();

   TString title=vtxAOD->GetTitle();
   if(!title.Contains("VertexerTracks")) return 0;

   AliVertexerTracks vertexer(aod->GetMagneticField());
   vertexer.SetITSMode();
   vertexer.SetMinClusters(3);
   if(title.Contains("WithConstraint")) {
     Float_t diamondcovxy[3];
     aod->GetDiamondCovXY(diamondcovxy);
     Double_t pos[3]={aod->GetDiamondX(),aod->GetDiamondY(),0.};
     Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
     AliESDVertex diamond(pos,cov,1.,1);
     vertexer.SetVtxStart(&diamond);
   }

   //_____________________________
   //Determination of unsuited tracks
   Int_t skipped[5000]; for(Int_t i=0;i<5000;i++) skipped[i]=-1;
   Int_t id = (Int_t)track->GetID();
   if(!(id<0)) skipped[0] = id;  //remove track under investigation from vertex
   int nTrksToSkip=1;
   //nTrksToSkip=DetermineUnsuitableVtxTracks(skipped, aod, track);
   vertexer.SetSkipTracks(nTrksToSkip,skipped);

   //________________________________
   //Determination of new ESD vertex
   AliESDVertex *vtxESDNew = vertexer.FindPrimaryVertex(aod);
   if(!vtxESDNew) return 0;
   Int_t nContrib =vtxESDNew->GetNContributors();
   if(vtxESDNew->GetNContributors()<fAnalysisCuts[bAnalysisCut_MinNewVertexContrib]) {
     //printf("Thrown away with %i contributors\n",vtxESDNew->GetNContributors());
     delete vtxESDNew; vtxESDNew=nullptr;
     return 0;
   }
   if(vtxESDNew->GetChi2toNDF()>fAnalysisCuts[bAnalysisCut_Z_Chi2perNDF]) {
     //printf("Thrown away with chi2 = %f\n",fAnalysisCuts[bAnalysisCut_Z_Chi2perNDF]);
     delete vtxESDNew; vtxESDNew=nullptr;
     return 0;
   }

   //________________________________
   //Conversion to AOD vertex
   Double_t pos[3];
   Double_t cov[6];
   Double_t chi2perNDF;
   vtxESDNew->GetXYZ(pos); // position
   vtxESDNew->GetCovMatrix(cov); //covariance matrix
   chi2perNDF = vtxESDNew->GetChi2toNDF(); //chisquare
   if(vtxESDNew) delete vtxESDNew;
   vtxESDNew=NULL;
   AliAODVertex *vtxAODNew = new AliAODVertex(pos,cov,chi2perNDF);
   vtxAODNew->SetNContributors(nContrib);  //contributors
   return vtxAODNew;
}
/*void AliAnalysisTaskHFJetIPQA::FillParticleCompositionSpectra(AliEmcalJet * jet,const char * histname ){
    if(!jet) return;
    AliVTrack* tr=0x0;
    for(Int_t j = 0; j < jet->GetNumberOfTracks(); ++j) {
      tr = (AliVTrack*)GetParticleContainer(0)->GetParticle((jet->TrackAt(j)));
      if(!tr) continue;
      double pT=0x0;
      Int_t pCorr_indx=-1;
      GetMonteCarloCorrectionFactor(tr,pCorr_indx,pT);
      if(pCorr_indx<0 ) continue;
      FillHist(histname,pCorr_indx+0.5,pT, 1);     //*this->fXsectionWeightingFactor );
  }
  return;
}*/

AliEmcalJet *  AliAnalysisTaskHFJetIPQA::GetPerpendicularPseudoJet (AliEmcalJet *jet_in  , bool rev ){
    TVector3 j(jet_in->Px(), jet_in->Py(), jet_in->Pz());
    TVector3 p1(j);
    std::vector <int > track_inc;
    p1.RotateZ(rev ? -1*TMath::Pi()/2. :TMath::Pi()/2. );
    Double_t sumAllPt1 = 0;
    int nconst_1 =0;
    std::vector <int> const_idx1;
    for(long itrack= 0; itrack<GetParticleContainer(0)->GetNParticles();++itrack){
       AliVTrack *  tr = static_cast<AliVTrack*>(GetParticleContainer(0)->GetParticle(itrack));
       if(!tr) continue;
       TVector3 v(tr->Px(), tr->Py(), tr->Pz());
       Double_t dR1 = v.DrEtaPhi(p1);
       if(v.Pt()>0.150){
        if(dR1 < fJetRadius) {
            sumAllPt1+=v.Pt();
            nconst_1++;
            const_idx1.push_back(itrack);
        }
    }
}
AliEmcalJet* jet1 =0;


if (sumAllPt1>0) {
    jet1 = new AliEmcalJet(sumAllPt1, p1.Eta(), TVector2::Phi_0_2pi (p1.Phi()), 0);
    jet1->SetArea(fJetRadius*fJetRadius*TMath::Pi());
    jet1->SetNumberOfTracks(nconst_1);
    jet1->SetNumberOfClusters(0);
    for (int i = 0 ; i < (int) const_idx1.size();++i) {
        jet1->AddTrackAt(const_idx1.at(i), i);
    }}


    return jet1;
}










/*! \brief CalculateTrackImpactParameter
 *
 *
 * Calculate track impact parameter
 */

Double_t AliAnalysisTaskHFJetIPQA::GetValImpactParameter(TTypeImpPar type,Double_t *impar, Double_t * cov)
{
    Float_t result = -999999;
    Float_t dFdx   = 0;
    Float_t dFdy   = 0;
    switch(type){
        case kXY:
        result = impar[0];
        break;
        case kXYSig:
        result =  impar[0]/sqrt(cov[0]);
        break;
        case kXYZ:
        result = TMath::Sqrt(impar[0]*impar[0]+impar[1]*impar[1]);
        break;
        case kXYZSig:
        result =  TMath::Sqrt(impar[0]*impar[0]+impar[1]*impar[1]);
        dFdx = impar[0]/result;
        dFdy = impar[1]/result;
        result /=TMath::Sqrt(cov[0]*dFdx*dFdx + cov[2]*dFdy*dFdy + 2* cov[1] *dFdx*dFdy);
        break;
        case kZSig:
        result =  impar[1];
        result /=TMath::Sqrt(cov[2]);
        break;
        case kXYZSigmaOnly:
        result =  TMath::Sqrt(impar[0]*impar[0]+impar[1]*impar[1]);
        dFdx = impar[0]/result;
        dFdy = impar[1]/result;
        result =TMath::Sqrt(cov[0]*dFdx*dFdx + cov[2]*dFdy*dFdy + 2* cov[1] *dFdx*dFdy);
        break;
        default:
        break;
    }
    return result;
}

//____________________________________________________
void AliAnalysisTaskHFJetIPQA::FillCandidateJet(Int_t CutIndex, Int_t JetFlavor){
        if(JetFlavor>0) fh1DCutInclusive->Fill(CutIndex);
        if(fIsPythia){
                if(JetFlavor==1)fh1dCutudg->Fill(CutIndex);
                if(JetFlavor==2)fh1dCutc->Fill(CutIndex);
                if(JetFlavor==3)fh1dCutb->Fill(CutIndex);
                if(JetFlavor==4)fh1dCuts->Fill(CutIndex);
        }

}

//Bool_t AliAnalysisTaskHFJetIPQA::IsTrackAcceptedMC(Double_t pt, Double_t eta){
//  if(TMath::Abs(eta)>0.9) return kFALSE;
//  if(pt<fAnalysisCuts[bAnalysisCut_MinTrackPt]) return kFALSE;
//  return kTRUE;
//}

Bool_t AliAnalysisTaskHFJetIPQA::IsTrackAccepted(const AliVTrack* track , int jetflavour){
    if(!track) return kFALSE;
    if(fIsEsd){
        fESDTrackCut->SetMinNClustersITS(fAnalysisCuts[bAnalysisCut_MinITSLayersHit]);
        fESDTrackCut->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
        if(!(fESDTrackCut->AcceptTrack((AliESDtrack*)track))) return kFALSE;
        return kTRUE;
    }
    else {
            //HasMatchedGoodTracklet((AliAODTrack*)track);
        if(!(((AliAODTrack*)track)->TestFilterBit(9) || ((AliAODTrack*)track)->TestFilterBit(4)))return kFALSE;

        if(TMath::Abs(track->Eta())>0.9) return kFALSE;

        //SPD hits
        if(fAnalysisCuts[bAnalysisCut_HasSPD]){
          if(!(((AliAODTrack*)track->HasPointOnITSLayer(0))&&(AliAODTrack*)track->HasPointOnITSLayer(1))){
            FillCandidateJet(bAnalysisCut_HasSPD,jetflavour);
            //printf("Throw away as SPD hits missing: 0=%i, 1=%i\n",track->HasPointOnITSLayer(0),track->HasPointOnITSLayer(1));
            return kFALSE;
          }
        }

        //SDD hits
        if(fAnalysisCuts[bAnalysisCut_HasSDD]){
          if(!(((AliAODTrack*)track->HasPointOnITSLayer(2))&&(AliAODTrack*)track->HasPointOnITSLayer(3))){
            FillCandidateJet(bAnalysisCut_HasSDD,jetflavour);
            //printf("Throw away as SDD hits missing: 2=%i, 3=%i\n",track->HasPointOnITSLayer(2),track->HasPointOnITSLayer(3));
            return kFALSE;
          }
        }

        //SSD hits
        if(fAnalysisCuts[bAnalysisCut_HasSSD]){
          if(!(((AliAODTrack*)track->HasPointOnITSLayer(4))&&(AliAODTrack*)track->HasPointOnITSLayer(5))){
            FillCandidateJet(bAnalysisCut_HasSSD,jetflavour);
            //printf("Throw away as SSD hits missing: 4=%i, 5=%i\n",track->HasPointOnITSLayer(4),track->HasPointOnITSLayer(5));
            return kFALSE;
          }
        }

        //n=hits in ITS layers
        Int_t SPDSSDHits = (int) track->HasPointOnITSLayer(0) + (int) track->HasPointOnITSLayer(1) + (int) track->HasPointOnITSLayer(2) + (int) track->HasPointOnITSLayer(3)+(int) track->HasPointOnITSLayer(4) + (int) track->HasPointOnITSLayer(5);
        if(SPDSSDHits<abs(fAnalysisCuts[bAnalysisCut_MinITSLayersHit])){
            //printf("Throw away due flav=%i ssd points 0 = %i, 1=%i, 2=%i, 3=%i, 4=%i, 5=%i\n",jetflavour,track->HasPointOnITSLayer(0),track->HasPointOnITSLayer(1),track->HasPointOnITSLayer(2),track->HasPointOnITSLayer(3),track->HasPointOnITSLayer(4),track->HasPointOnITSLayer(5));
            FillCandidateJet(bAnalysisCut_MinITSLayersHit,jetflavour);
            return kFALSE;
        }

        //TPC clusters
        if(((AliAODTrack*)track)->GetNcls(1)<fAnalysisCuts[bAnalysisCut_MinTPCClus]){
            //printf("Throw away due flav=%i TPCClus %i, cutvalue=%f\n",jetflavour,((AliAODTrack*)track)->GetNcls(1),fAnalysisCuts[bAnalysisCut_MinTPCClus]);
            FillCandidateJet(bAnalysisCut_MinTPCClus,jetflavour);
            return kFALSE;
        }

        if(track->Pt()<fAnalysisCuts[bAnalysisCut_MinTrackPt]){
            //printf("Throw away due flav=%i pt %f, cutvalue=%f\n",jetflavour,track->Pt(),fAnalysisCuts[bAnalysisCut_MinTrackPt]);
            FillCandidateJet(bAnalysisCut_MinTrackPt,jetflavour);
            return kFALSE;
        }

        AliAODVertex *aodvertex = (( AliAODTrack *)track)->GetProdVertex();
        if(!aodvertex) return kFALSE;
        if(fAnalysisCuts[bAnalysisCut_KinkCand]){
          if(aodvertex->GetType()==AliAODVertex::kKink){
            FillCandidateJet(bAnalysisCut_KinkCand,jetflavour);
            return kFALSE;
          }
        }

        ULong_t status=track->GetStatus();
        if(fAnalysisCuts[bAnalysisCut_HasTPCrefit]){
          if(!(status & AliAODTrack::kTPCrefit)){
            FillCandidateJet(bAnalysisCut_HasTPCrefit,jetflavour);
            return kFALSE;
          }
        }

        if(fAnalysisCuts[bAnalysisCut_HasITSrefit]){
          if(!(status & AliAODTrack::kITSrefit)){
            FillCandidateJet(bAnalysisCut_HasITSrefit,jetflavour);
            return kFALSE;
          }
        }
        if(((AliAODTrack*)track)->Chi2perNDF()>=fAnalysisCuts[bAnalysisCut_MinTrackChi2]){
            FillCandidateJet(bAnalysisCut_MinTrackChi2,jetflavour);
            //printf("Throw away due flav=%i chi2 %f, cutvalue=%f\n",jetflavour,((AliAODTrack*)track)->Chi2perNDF(),fAnalysisCuts[bAnalysisCut_MinTrackChi2]);
            return kFALSE;
        }

        return kTRUE;
    }
    return kTRUE;
}


bool AliAnalysisTaskHFJetIPQA::IsDCAAccepted(double decaylength, double ipwrtjet, Double_t * dca, int jetflavour){
    if(dca[0]>fAnalysisCuts[bAnalysisCut_MaxDCA_XY]){
        FillCandidateJet(bAnalysisCut_MaxDCA_XY,jetflavour);
        //printf("Throw away due to xy cut = %f, cutvalue=%f",dca[0],fAnalysisCuts[bAnalysisCut_MaxDCA_XY]);
        return kFALSE;
    }
    if(dca[1]>fAnalysisCuts[bAnalysisCut_MaxDCA_Z]){
        FillCandidateJet(bAnalysisCut_MaxDCA_Z,jetflavour);
        //printf("Throw away due to z cut = %f, cutvalue=%f",dca[1],fAnalysisCuts[bAnalysisCut_MaxDCA_Z]);
        return kFALSE;
    }
    if(decaylength>fAnalysisCuts[bAnalysisCut_MaxDecayLength]){
        FillCandidateJet(bAnalysisCut_MaxDecayLength,jetflavour);
        //printf("Throw away due to decaylength cut = %f, cutvalue=%f",decaylength,fAnalysisCuts[bAnalysisCut_MaxDecayLength]);
        return kFALSE;
    }
    if(ipwrtjet>fAnalysisCuts[bAnalysisCut_DCAJetTrack]){
        FillCandidateJet(bAnalysisCut_DCAJetTrack,jetflavour);
        //printf("Throw away due to dcajettrack = %f, cutvalue=%f",ipwrtjet,fAnalysisCuts[bAnalysisCut_DCAJetTrack]);
        return kFALSE;
    }

    return kTRUE;
}

/*! \brief FillTrackingEfficiencyDCA
 *
 *
 * Calculates tracking efficiency as a function of the dca xy z
 */

/*! \brief MatchJetsGeometricDefault
 *
 *
 * geometric jet matching
 */
Bool_t AliAnalysisTaskHFJetIPQA::MatchJetsGeometricDefault()
{
    AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));
    Double_t matchingpar1 =0.3;
    Double_t matchingpar2 =0.3;
    if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return kFALSE;
    DoJetLoop();
    AliEmcalJet* jet1 = 0;
    jets1->ResetCurrentID();
    while ((jet1 = jets1->GetNextJet())) {
        AliEmcalJet *jet2 = jet1->ClosestJet();
        if (!jet2) continue;
        if (jet2->ClosestJet() != jet1) continue;
        if (jet1->ClosestJetDistance() > matchingpar1 || jet2->ClosestJetDistance() > matchingpar2) continue;
            // Matched jet found
        jet1->SetMatchedToClosest(1);
        jet2->SetMatchedToClosest(1);
    }
    return kTRUE;
}

/*! \brief DoJetLoop
 *
 *
 * jet matching loop:
 * - reset previous matching and jet IDs
 * - redo jet matching
 */
void AliAnalysisTaskHFJetIPQA::DoJetLoop()
{
    // Do the jet loop.
    Double_t minjetpt =1.;
    AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));
    if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return;
    AliEmcalJet* jet1 = 0;
    AliEmcalJet* jet2 = 0;
    jets2->ResetCurrentID();
    while ((jet2 = jets2->GetNextJet())) jet2->ResetMatching();
    jets1->ResetCurrentID();
    while ((jet1 = jets1->GetNextJet())) {
        jet1->ResetMatching();
        if (jet1->MCPt() < minjetpt) continue;
        jets2->ResetCurrentID();
        while ((jet2 = jets2->GetNextJet())) {
            SetMatchingLevel(jet1, jet2, 1);
                } // jet2 loop
        } // jet1 loop
    }
/*! \brief IsTruePrimary
 *
 *
 *
 */
    Bool_t AliAnalysisTaskHFJetIPQA::IsTruePrimary(AliVParticle * mcpart){
        if(!mcpart) return kFALSE;
        AliVParticle * mcmother = GetVParticleMother(mcpart);
        if(!mcmother) return kTRUE;
        Int_t istatus =-1;
        if(!fIsEsd) {
            istatus =   ( (AliMCParticle*)mcmother)->Particle()->GetStatusCode();
        }
        else {
            istatus =    ((AliAODMCParticle*)mcmother)->GetStatus();
        }
        if(istatus >11)return kTRUE;
        return kFALSE;
    }


/*! \brief GetBMesonWeight
 *
 *
 */
Bool_t AliAnalysisTaskHFJetIPQA::GetBMesonWeight( AliVParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  )
{
    pT = mcpart->Pt();
    switch(pdg){
        case bBPlus:
        idx = bIdxBPlus;
        return kTRUE;
        case bB0:
        idx = bIdxB0;
        return kTRUE;
        case bLambdaB:
        idx = bIdxLambdaB;
        return kTRUE;
        break;
        case bBStarPlus:
        idx = bIdxBStarPlus;
        return kTRUE;
        break;
        default:
        break;
    }
    return kFALSE;
}
/*! \brief IsSelectionParticle
 *
 *
 */
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticle( AliVParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
    pT 	= mcpart->Pt();
    Int_t pdg2 = abs(mcpart->PdgCode());
    idx = -1;

    switch(pdg2){
        case bPi0:
        idx = bIdxPi0;
        if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
        break;
        case bEta:
        idx = bIdxEta;
        if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
        break;
        case bEtaPrime:
        idx = bIdxEtaPrime;
        if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
        break;
        case bOmega:
        idx = bIdxOmega;
        if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
        break;
        case bPhi:
        idx = bIdxPhi;
        if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
        break;
        case bRho:
        idx = bIdxRho;
        if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
        break;
        case bRhoPlus:
        idx = bIdxRho;
        if(!IsSecondaryFromWeakDecay(mcpart))return kTRUE;
        break;
        default:
        break;
    }
    return kFALSE;
}
/*! \brief IsSelectionParticleALICE
 *
 *
 */
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleALICE( AliVParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
    pT 	= mcpart->Pt();
    AliVParticle * mother = nullptr;
    idx = -1;
    pdg = abs(mcpart->PdgCode());
    Bool_t pIsSecStrANGE= IsSecondaryFromWeakDecay(mcpart);
    if(pIsSecStrANGE ) return kFALSE;
    if(!IsPhysicalPrimary(mcpart))return kFALSE;

    switch(pdg){
        case bProton:
        mother = GetVParticleMother(mcpart);
        if(mother){
            if((abs(mother->PdgCode()) ==  3222) )
                return kFALSE;
        }
        idx=bIdxProton;
        return kTRUE;
        break;
        case bPi:
        idx=bIdxPi;
        return kTRUE;
        break;
        case bKaon:
        idx=bIdxKaon;
        return kTRUE;
        break;
        default:
        return kFALSE;
        break;
    }
    return kFALSE;
}
/*! \brief IsSelectionParticleMeson
 *
 *
 */
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleMeson( AliVParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
    pT 	= mcpart->Pt();
    idx = -1;
    pdg = abs(mcpart->PdgCode());
    if(TMath::Abs(mcpart->Y()) >=.5) return kFALSE;
    if(IsSecondaryFromWeakDecay(mcpart))return kFALSE;
    switch(pdg){
        case bD0:
        idx = bIdxD0;
        if(IsPromptDMeson(mcpart))return kTRUE;
        break;
        case bDPlus:
        idx = bIdxDPlus;
        if(IsPromptDMeson(mcpart))return kTRUE;
        break;
        case bDSPlus:
        idx = bIdxDSPlus;
        if(IsPromptDMeson(mcpart))return kTRUE;
        break;
        case bDStarPlus:
        idx = bIdxDStarPlus;
        if(IsPromptDMeson(mcpart))return kTRUE;
        break;
        case bLambdaC:
        idx = bIdxLambdaC;
        if(IsPromptDMeson(mcpart))return kTRUE;
        break;
        case bBPlus:
        idx = bIdxBPlus;
        if(IsPromptBMeson(mcpart))return kTRUE;
        break;
        case bB0:
        idx = bIdxB0;
        if(IsPromptBMeson(mcpart))return kTRUE;
        break;
        case bLambdaB:
        idx = bIdxLambdaB;
        if(IsPromptBMeson(mcpart))return kTRUE;
        break;
        case bBStarPlus:
        idx = bIdxBStarPlus;
        if(IsPromptBMeson(mcpart))return kTRUE;
        break;
        default:
        return kFALSE;
        break;
    }
    return kTRUE;
}
/*! \brief IsSelectionParticleOmegaXiSigmaP
 *
 *
 */
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleOmegaXiSigmaP( AliVParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
    pT 	= mcpart->Pt();
    idx = -1;
    pdg = abs(mcpart->PdgCode());

    if (!IsPhysicalPrimary(mcpart)) return kFALSE;
    switch(pdg){
        case bSigmaMinus:
            idx = bIdxSigmaMinus; // use lambda as proxy
            return kTRUE;
            break;
            case bSigmaPlus:
            idx = bIdxSigmaPlus; //use lambda as proxy
            return kTRUE;
            break;
            case bXiBaryon:
            idx = bIdxBStarPlus;//dummy! Externally sotlved
            return kTRUE;
            break;
            case bOmegaBaryon:
            idx = bIdxBStarPlus;//dummy! Externally solved
            return kTRUE;
            break;
            default:
            return kFALSE;
            break;
        }
        return kFALSE;
    }
/*! \brief GetVParticleMother
 *
 *
 */
    AliVParticle * AliAnalysisTaskHFJetIPQA::GetVParticleMother(AliVParticle * part){
        AliVParticle * mother = nullptr;
        if (part->GetMother()<0) return nullptr;
        if(fIsEsd){
            mother = static_cast <AliMCParticle*>(MCEvent()->GetTrack(part->GetMother()));
        }
        else{
            mother =static_cast<AliAODMCParticle*>(fMCArray->At(part->GetMother()));
        }
        if(!mother) return nullptr;
        return mother;
    }
/*! \brief IsPhysicalPrimary
 *
 *
 */
    Bool_t  AliAnalysisTaskHFJetIPQA::IsPhysicalPrimary( AliVParticle * part){
        return  fIsEsd ? MCEvent()->IsPhysicalPrimary(part->GetLabel()) : static_cast<AliAODMCParticle*>(part)->IsPhysicalPrimary() ;
    }
/*! \brief IsSecondaryFromWeakDecay
 *
 *
 */
    Bool_t  AliAnalysisTaskHFJetIPQA::IsSecondaryFromWeakDecay( AliVParticle * part){
        return fIsEsd ? MCEvent()->IsSecondaryFromWeakDecay(part->GetLabel()) : static_cast<AliAODMCParticle*>(part)->IsSecondaryFromWeakDecay() ;

    }
/*! \brief IsSelectionParticleStrange
 *
 *
 */
    Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleStrange( AliVParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
        pT 	= mcpart->Pt();
        AliVParticle * mother = nullptr;
        idx = -1;
        pdg = abs(mcpart->PdgCode());
        if(!IsPhysicalPrimary(mcpart))return kFALSE;
        switch(pdg){
            case bPhi:
            idx = bIdxPhi;//dummy will be overwritten in posrpos
            return kTRUE;
            break;
            case bK0S892:
            idx = bIdxK0s;//dummy will be overwritten in posrpos
            return kTRUE;
            break;
            case bK0S892plus:
            idx = bIdxK0s; //dummy will be overwritten in posrpos
            return kTRUE;
            break;
            case bK0s:
            idx = bIdxK0s;
            mother = GetVParticleMother(mcpart);
            if(mother){
                if((abs(mother->PdgCode()) == bPhi))
                    return kFALSE;
            }
            return kTRUE;
            break;
            case bK0l:
            idx = bIdxK0s;
            mother = GetVParticleMother(mcpart);
            if(mother){
                if((abs(mother->PdgCode()) == bPhi))
                    return kFALSE;
            }
            return kTRUE;
            break;

            case bLambda:
            idx = bIdxLambda;
            mother = GetVParticleMother(mcpart);
            if(mother){
                if((abs(mother->PdgCode()) ==  3312) || (abs(mother->PdgCode()) ==  3322) || (abs(mother->PdgCode()) ==  3334))
                    return kFALSE;
            }
            return kTRUE;
            break;
            default:
            return kFALSE;
            break;
        }
        return kFALSE;
    }
/*! \brief IsPromptBMeson
 *
 *
 */
    Bool_t AliAnalysisTaskHFJetIPQA::IsPromptBMeson(AliVParticle * part )
    {
        if(!part) return kFALSE;
        Int_t pdg = TMath::Abs(part->PdgCode());
        if ((pdg >= 500 && pdg < 600 )||(pdg >= 5000 && pdg < 6000 ))
        {
            AliVParticle* pm = GetVParticleMother(part);
            Int_t mpdg = TMath::Abs(pm->PdgCode());
            if (!(mpdg >5000 && mpdg <6000) && !(mpdg >500 && mpdg <600))
                return kTRUE;
        }
        return kFALSE;
    }
/*! \brief IsPromptDMeson
 *
 *
 */
    Bool_t AliAnalysisTaskHFJetIPQA::IsPromptDMeson(AliVParticle * part )
    {
        if(!part) return kFALSE;
        Int_t pdg = TMath::Abs(part->PdgCode());
        if ((pdg >= 400 && pdg < 500 )||(pdg >= 4000 && pdg < 5000 ))
        {
            AliVParticle* pm = GetVParticleMother(part);
            if(!pm) return kTRUE;
            Int_t mpdg = TMath::Abs(pm->PdgCode());
            if (!(mpdg >4000 && mpdg <6000) && !(mpdg >400 && mpdg <600))
                return kTRUE;
        }

        return kFALSE;
    }
/*! \brief ParticleIsPossibleSource
 *
 *
 */
    Bool_t AliAnalysisTaskHFJetIPQA::ParticleIsPossibleSource(Int_t pdg){
        Int_t pos[22] = {bPi0,bEta,bEtaPrime,bPhi,bRho,bOmega,bK0s,bLambda,bOmegaBaryon,bXiBaryon,bD0,bPi,bKaon,bProton,bDPlus,bDStarPlus,bDSPlus,bLambdaB,bLambdaC,bBPlus,bB0,bBStarPlus};
        for (Int_t i =0 ;i<22 ;++i){
            if (abs(pdg)==pos[i] ) return kTRUE;
        }
        return kFALSE;
    }
/*! \brief SetMatchingLevel
 *
 * jet matching helper:
 * - define closest and second closest jet
 */
    void AliAnalysisTaskHFJetIPQA::SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Int_t matching)
    {
        Double_t d1 = - 1;
        Double_t d2 = -1;

        switch (matching) {
            case 1:
            GetGeometricalMatchingLevel(jet1,jet2,d1);
            d2 = d1;
            break;
            default:
            break;
        }
        if (d1 >= 0) {
            if (d1 < jet1->ClosestJetDistance()) {
                jet1->SetSecondClosestJet(jet1->ClosestJet(), jet1->ClosestJetDistance());
                jet1->SetClosestJet(jet2, d1);
            }
            else if (d1 < jet1->SecondClosestJetDistance()) {
                jet1->SetSecondClosestJet(jet2, d1);
            }
        }
        if (d2 >= 0) {
            if (d2 < jet2->ClosestJetDistance()) {
                jet2->SetSecondClosestJet(jet2->ClosestJet(), jet2->ClosestJetDistance());
                jet2->SetClosestJet(jet1, d2);
            }
            else if (d2 < jet2->SecondClosestJetDistance()) {
                jet2->SetSecondClosestJet(jet1, d2);
            }
        }
    }
/*! \brief GetGeometricalMatchingLevel
 *
 * jet matching helper
 */
    void AliAnalysisTaskHFJetIPQA::GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const
    {
        Double_t deta = jet2->Eta() - jet1->Eta();
        Double_t dphi = jet2->Phi() - jet1->Phi();
        dphi = TVector2::Phi_mpi_pi(dphi);
        d = sqrt(deta * deta + dphi * dphi);
    }

/*! \brief mysort
 *
 * custom strcut sorter function
 */
    Bool_t AliAnalysisTaskHFJetIPQA::mysort(const SJetIpPati& i, const SJetIpPati& j)
    {
        if(i.first <= j.first)
            return kFALSE;
        else
            return kTRUE;
    }
/*! \brief GetPtCorrected
 *
 *
 */
    Double_t AliAnalysisTaskHFJetIPQA::GetPtCorrected(const AliEmcalJet *jet)
    {
        AliJetContainer * jetconrec = nullptr;
        jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
        if(jet && jetconrec&&fDoUnderlyingEventSub){
            printf("Correct for Underlying Event.\n");
            return jet->Pt() - jetconrec->GetRhoVal() * jet->Area();
        }
        return -1.;
    }
/*! \brief GetPtCorrectedMC
 *
 *
 */
   /* Double_t AliAnalysisTaskHFJetIPQA::GetPtCorrectedMC(const AliEmcalJet *jet)
    {
        /*AliJetContainer * jetconrec = nullptr;
        jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(1));
        if(jet && jetconrec&&fDoUnderlyingEventSub){
            printf("Correct for Underlying Event.\n");
            return jet->Pt() - jetconrec->GetRhoVal() * jet->Area();
        }
        return -1.;
        return jet->Pt();
    }*/


/*! \brief IsJetTaggedTC
 * unused
 *
 */
    Bool_t AliAnalysisTaskHFJetIPQA::IsJetTaggedTC(Int_t n, Double_t thres)
    {
        return kTRUE;
    }


/*! \brief IsParton
 *
 *
 */
    Bool_t AliAnalysisTaskHFJetIPQA::IsParton(int pdg){
      return ((pdg==1)||(pdg==2)||(pdg==3)||(pdg==4)||(pdg==5)||(pdg==21));
    }

/*! \brief IsJetTaggedJetProb
 *
 * unused
 */
    Bool_t AliAnalysisTaskHFJetIPQA::IsJetTaggedJetProb(Double_t thresProb)
    {
        return kTRUE;
    }
    /*! \brief IsMCJetPartonFast
     *
     * Fast jet parton MC matcher
     */
        Int_t  AliAnalysisTaskHFJetIPQA::IsMCJetPartonFast(const AliEmcalJet *jet, Double_t radius,Bool_t &is_udg)
        {
            fJetCont.clear();
            fPUdsgJet.clear();
            fPSJet.clear();
            fPCJet.clear();
            fPBJet.clear();
            daughtermother.clear();

            double p_udsg_max=-999;
            double p_s_max=-999;
            AliVParticle *vp=0x0;
            Int_t pdg=0;
            Double_t p=0;
            int kJetOrigin=-999;

            //printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
            //printf("Starting loop\n");
            if(!jet) return 0;
            if(!(jet->GetNumberOfTracks()>fNoJetConstituents)){
              //printf("Throwing away jets with only too few constituent!\n");
              return 0;
            }

            if(fDoFlavourMatching){
              for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {//start trackloop jet
                vp = static_cast<AliVParticle*>(jet->Track(i));
                if (!vp){
                  AliError("AliVParticle associated to constituent not found\n");
                  continue;
                }

                AliAODMCParticle * part = static_cast<AliAODMCParticle*>(vp);

                if(!part){
                    AliError("Finding no Part!\n");
                    return 0;
                }       // if(!part->IsPrimary()) continue;
                pdg = (abs(part->PdgCode()));

                fJetCont.push_back(part->Label());
                //printf("Daugther pdg=%i, Label=%i, Mother =%i, p=%f, MCStatusCode=%i\n",pdg, part->GetLabel(), part->GetMother(), p, part->MCStatusCode());
              }//end trackloop jet
            }

            //for(Int_t iPrim = 0 ; iPrim<fMCArray->GetEntriesFast();iPrim++){//start trackloop MC
            for(Int_t iPrim=0; iPrim<fMCEvent->GetNumberOfTracks();iPrim++){

                        //AliAODMCParticle * part = static_cast<AliAODMCParticle*>(fMCArray->At(iPrim));
                        AliAODMCParticle * part = static_cast<AliAODMCParticle*>(fMCEvent->GetTrack(iPrim));
                        if(!part) return 0;
                        if(!part->IsPrimary()) continue;
                        Double_t eta = part->Eta();
                        Double_t phi= part->Phi();
                        Double_t etajet = jet->Eta();
                        Double_t phijet = jet->Phi();
                         p=part->P();

                        Int_t pdg = (abs(part->PdgCode()));
                        Double_t deta = etajet - eta;
                        Double_t dphi = phijet-phi;
                        dphi = TVector2::Phi_mpi_pi(dphi);
                        Double_t  d = sqrt(deta * deta + dphi * dphi);
                      //   printf("LINE %i, deta%f, dphi%f, d=%f, fDoFlavourMatching=%i\n",__LINE__, deta,dphi,d,fDoFlavourMatching);
                        if(!fDoFlavourMatching) {
                          //if(!((part->GetStatus()==11) ||(part->GetStatus()==12))) continue;
                          if(!IsParton(pdg)) continue;
                          if(d > radius) continue;
                          kJetOrigin=pdg;
                        }
                        else{

                          if(IsParton(pdg)){
                            if(d > fDaughtersRadius) continue;
                          }

                          //printf("i=%i, Mother Label %i, Mother PdG %i,Daughter:%i, Last Daughter: %i, MCStatusCode=%i, d=%f, p=%f\n", iPrim, part->Label(), part->PdgCode(), part->GetDaughterLabel(0),part->GetDaughterLabel(1),part->MCStatusCode(), d,p);

                          int kFirstDaugh=part->GetDaughterLabel(0);
                          int NDaugh=part->GetNDaughters();
                          for(int idaugh=0;idaugh<NDaugh;idaugh++){
                            int kDaughLabel=kFirstDaugh+idaugh;
                            //printf("Dauglabel=%i, kFirstDaugh=%i, kLastDaugh=%i\n",kDaughLabel, kFirstDaugh,kLastDaugh);

                            bool IsDaughter=std::find(fJetCont.begin(), fJetCont.end(),kDaughLabel) != fJetCont.end();
                            if(IsDaughter){
                              if(IsParton(pdg)){
                                //printf("Directly matched %i with daughter =%i\n",part->GetLabel(), kDaughLabel);
                                kJetOrigin=part->PdgCode();
                              }
                              else{
                                  bool Is2ndDaughter=daughtermother.find(part->Label()) != daughtermother.end();
                                  if(Is2ndDaughter){
                                      kJetOrigin=daughtermother.find(part->Label())->second;
                                      //printf("Matched with %i with 2nd daughter =%i\n",part->GetLabel(), kDaughLabel);
                                  }
                              }
                            }//end is jet daughter
                            else{
                              if(IsParton(pdg)){
                                //printf("Writing Quarks in map: partlabel=%i, daughlabel=%i, status=%i\n", part->GetLabel(),kDaughLabel, part->MCStatusCode());
                                daughtermother.emplace(kDaughLabel, part->PdgCode());
                              }
                              else{
                                bool Is2ndDaughter=daughtermother.find(part->Label()) != daughtermother.end();
                                if(Is2ndDaughter){
                                  //printf("Writing Daughters in map: partlabel=%i, daughlabel=%i, status=%i\n", part->GetLabel(),kDaughLabel, part->MCStatusCode());
                                  int kOriginalQuark=daughtermother.find(part->Label())->second;
                                  daughtermother.emplace(kDaughLabel, kOriginalQuark);
                                }
                              //printf("Asking for daughlabel%i\n",part->Label());
                              }
                            }//end related to jet?
                          }//end daughterloop
                        }//end else DoMatchFlavours

                        //printf("i=%i, Mother Label %i, Mother PdG %i,Daughter:%i, Last Daughter: %i, MCStatusCode=%i, d=%f, p=%f\n", iPrim, part->Label(), part->PdgCode(), part->GetDaughterLabel(0),part->GetLastDaughter(),part->MCStatusCode(), d,p);

                        if(abs(kJetOrigin) == 5) {
                            fPBJet.push_back(p);
                        }
                        else if(abs(kJetOrigin)== 4) {
                            fPCJet.push_back(p);
                        }
                        else if(abs(kJetOrigin) == 3 ) {
                            fPSJet.push_back(p);
                            //printf("Strange pushed with p=%f",p);
                        }
                        else if(abs(kJetOrigin)== 1 ||abs(kJetOrigin)== 2 ||  abs(kJetOrigin) == 21) {
                            fPUdsgJet.push_back(p);
                            //printf("Light pushed with p=%f",p);
                        }

        }//end trackloop MC

        /*printf("Inside JetCont:\n");
          for(int i=0;i<fJetCont.size();i++){
          printf("%f\n", fJetCont[i]);
        }
        printf("Inside map:\n");
          for (auto& x: daughtermother) {
          std::cout << x.first << ": " << x.second << '\n';
        }*/
        if(fPCJet.size() ==0&& fPBJet.size()==0&& fPSJet.size()==0&&fPUdsgJet.size()==0) return 0; //udsg
        //check for b jet
        for (Int_t icj = 0 ; icj <(Int_t)fPBJet.size();++icj ){
            //printf("Bottom Flavour Jet!\n");
            return B;
        }
        //check for c jet
        for (Int_t icj = 0 ; icj <(Int_t)fPCJet.size();++icj ){
            //printf("Charm Flavour Jet!\n");
            return C;
        }
        //check for s and light jet
        if(fPUdsgJet.size()!=0){
          std::sort(fPUdsgJet.begin(), fPUdsgJet.end());
          p_udsg_max=fPUdsgJet[fPUdsgJet.size()-1];
        }
        if(fPSJet.size()!=0){
          std::sort(fPSJet.begin(), fPSJet.end());
          p_s_max=fPSJet[fPSJet.size()-1];
        }

        if(p_s_max>p_udsg_max){
          //printf("S prefered with psmax=%f, pudsgmax=%f\n", p_s_max,p_udsg_max);
          return UDSG;
        }
        else{
            if(fPUdsgJet.size()!=0){
                //printf("Light prefered with psmax=%f, pudsgmax=%f\n", p_s_max,p_udsg_max);
                is_udg =kTRUE;
                return UDSG;
            }
        }
        return 0;
    }


/*! \brief RecursiveParents
 *
 * function which is declustering jets via Camebridge Aachen algorithm and from subjets filling the Lund plane
  */
//_________________________________________________________________________
void AliAnalysisTaskHFJetIPQA::RecursiveParents(const AliEmcalJet *fJet,const AliJetContainer *fJetCont, vector<Int_t> &fJetConstTrackID){
   //printf("Entering recursive parents!\n");
   Int_t nall=0;
   double delta_R=-99;
   double z=-99;
   double zcut=-99;
   double y=-99;
   double lnpt_rel=-99;
   double yh=-99;
   std::vector<fastjet::PseudoJet>  fInputVectors;
   std::vector<fastjet::PseudoJet>   fOutputJets;
   fastjet::PseudoJet  PseudoTracks;
   fastjet::PseudoJet jj;
   fastjet::PseudoJet j1;
   fastjet::PseudoJet j2;

   AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();
   fInputVectors.clear();

   //Fill InputVector, set user index to track ID + 100
   if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
     AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
     AliVTrack *vtrack = dynamic_cast<AliVTrack*>(fTrk);
     if (!vtrack) {
       AliError(Form("Could not receive track%d\n", i));
       continue;
     }
     AliAODTrack *trackV = dynamic_cast<AliAODTrack*>(vtrack);

     if(!trackV)            continue;
     if(!IsTrackAccepted((AliAODTrack*)trackV,-1)) continue;
     PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
     PseudoTracks.set_user_index(fJet->TrackAt(i));   //user index ist erstmal irrelevant für mich...
     //printf("Track %i, px=%f, py=%f, pz=%f, e=%f, userindex=%i\n", i,fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E(),fJet->TrackAt(i));
     fInputVectors.push_back(PseudoTracks);
   }
   if(fInputVectors.size()==0){return;}
   fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
   fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);  //jet algorithm, jet radius, recomb.scheme=0 is E-scheme, clustering strategy

   //try declustering
   try {
      fastjet::ClusterSequence fClustSeqSA(fInputVectors, fJetDef);
      fOutputJets.clear();
      jj.reset(0,0,0,0);
      j1.reset(0,0,0,0);
      j2.reset(0,0,0,0);
      fOutputJets=fClustSeqSA.inclusive_jets(0);
      fOutputJets=sorted_by_pt(fOutputJets);
      if(fOutputJets.size()==0){return;}
      jj=fOutputJets[0];

      fMeanLNKt=0;
      fMeanTheta=0;
      fMeanLNKtSD=0;
      fMeanThetaSD=0;

      while((jj.has_parents(j1,j2))){  //&&(z<fAnalysisCuts[bAnalysisCut_SDz])
          delta_R=-99;
          z=-99;
          zcut=-99;
          y=-99;
          lnpt_rel=-99;
          yh=-99;

          nall=nall+1;
          if(j1.perp() < j2.perp()) swap(j1,j2);

          delta_R=j1.delta_R(j2);
          z=j2.perp()/(j1.perp()+j2.perp());
          y =log(1.0/delta_R);
          lnpt_rel=log(j2.perp()*delta_R);
          yh=j1.e()+j2.e();
          zcut=fAnalysisCuts[bAnalysisCut_SDz]*pow((delta_R/fJetRadius),fAnalysisCuts[bAnalysisCut_SDbeta]);
          //printf("Recursive Parents:: Cut decision zcut=%f, z=%f, DeltaRij=%f\n",zcut, z, delta_R);

          fMeanLNKt=fMeanLNKt+lnpt_rel;
          fMeanTheta=fMeanTheta+delta_R;
          //printf("Not SDped: j1pt=%f, j2pt=%f, delta_R=%f, z=%f, y=%f, lnpt_rel=%f, yh=%f, fMeanLNKt=%f, fMeanTheta=%f\n",
          //          j1.perp(),j2.perp(), delta_R, z,y,lnpt_rel, yh,  fMeanLNKt, fMeanTheta);
          if(z>zcut){
            //printf("Stop reclustering!\n");
            break;
          }

          fMeanLNKtSD=fMeanLNKtSD+lnpt_rel;
          fMeanThetaSD=fMeanThetaSD+delta_R;
          //Double_t LundEntries[6] = {y,lnpt_rel,fOutputJets[0].perp(),nall,yh,flagSubjet};
          //fHLundIterative->Fill(LundEntries);
          //printf("SDped: j1pt=%f, j2pt=%f, delta_R=%f, z=%f, y=%f, lnpt_rel=%f, yh=%f,  fMeanLNKt=%f, fMeanTheta=%f\n",
          //j1.perp(),j2.perp(), delta_R, z,y,lnpt_rel, yh, fMeanLNKt, fMeanTheta);
          jj=j1;
        }
        vector<fastjet::PseudoJet> fGroomedJetConstit=sorted_by_pt(jj.constituents());

        for(long unsigned iConst=0;iConst<fGroomedJetConstit.size();iConst++){
          //printf("Pushing const=%i, userindex=%i\n", iConst, fGroomedJetConstit[iConst].user_index());
          fJetConstTrackID.push_back(fGroomedJetConstit[iConst].user_index());
        }
      } catch (fastjet::Error) {
        AliError(" [w] FJ Exception caught.");
        //return -1;
      }

      return;
}
/*! \brief FillHist
 *
 * 1d
 */
void AliAnalysisTaskHFJetIPQA::FillHist(const char *name, Double_t x, Double_t w){
    TH1D * h1 =GetHist1D(name);
    if(h1){  h1->Fill(x,w);}
    else{
      AliError(Form("FillHist: Histogram with name %s not  available!\n", name));
    }
}
/*! \brief FillHist
 *
 * 2d
 */
void AliAnalysisTaskHFJetIPQA::FillHist(const char *name, Double_t x, Double_t y, Double_t w){
    TH2D * h2 =GetHist2D(name);
    if(h2){ h2->Fill(x,y,w);}
    else{
      AliError(Form("FillHist: Histogram with name %s not  available!\n", name));
    }
}
/*! \brief IncHist
 *
 * increase 1d hist bin
 */
void AliAnalysisTaskHFJetIPQA::IncHist(const char *name, Int_t bin){
    TH1D * h1 =GetHist1D(name);
    h1->SetBinContent(bin,h1->GetBinContent(bin)+1);
}
/*! \brief AddHistogramm
 *
 *add a historgram
 */
TH1 *AliAnalysisTaskHFJetIPQA::AddHistogramm(const char *name, const char *title, Int_t x, Double_t xlow, Double_t xhigh, Int_t y, Double_t ylow, Double_t yhigh){
    TObject * res = nullptr;
    res = fOutput->FindObject(name);
    if((res)) return nullptr;

    TH1 * phist=nullptr;
    if(y==0){ //TH1D*
        phist = new TH1D (name,title,x,xlow,xhigh);
    }
    else  {
        phist = new TH2D(name,title,x,xlow,xhigh,y,ylow,yhigh);
    }
    phist->Sumw2();

    TString sName(name);
    if(sName.EqualTo("fh1dCutsPrinted")){
        fOutput->AddFirst(phist);
    }
    else{
        fOutput->Add(phist);
    }
    return (TH1*)phist;
}

void AliAnalysisTaskHFJetIPQA::setFFillCorrelations(const Bool_t &value)
{
    fFillCorrelations = value;
}



void AliAnalysisTaskHFJetIPQA::setFMCglobalDCASmear(const Double_t value)
{
    fMCglobalDCASmear = value;
}

Double_t AliAnalysisTaskHFJetIPQA::getFVertexRecalcMinPt() const
{
    return fVertexRecalcMinPt;
}

void AliAnalysisTaskHFJetIPQA::setFVertexRecalcMinPt(const Double_t &value)
{
    fVertexRecalcMinPt = value;
}

Double_t AliAnalysisTaskHFJetIPQA::getFMCglobalDCAxyShift() const
{
    return fMCglobalDCAxyShift;
}

void AliAnalysisTaskHFJetIPQA::setFMCglobalDCAxyShift(const Double_t &value)
{
    fMCglobalDCAxyShift = value;
}





/*! \brief SubtractMean
 *
 *
 */


//Currently not in use
/*void AliAnalysisTaskHFJetIPQA::SubtractMean(Double_t val[], AliVTrack *track){
    Double_t  deltamean=fGraphMean->Eval(track->Pt() < 3. ? track->Pt() : 3 );//(fCurrentMeanFactors[0]-fCurrentMeanFactors[1]*TMath::Exp(-1*fCurrentMeanFactors[2]*(track->Pt()-fCurrentMeanFactors[3]))) *1e-4;

    val[0] -=deltamean*1e-4;;
}
//Helpers*/


Bool_t AliAnalysisTaskHFJetIPQA::GetImpactParameter(const AliAODTrack *track, const AliAODEvent *event, Double_t *dca, Double_t *cov, Double_t *XYZatDCA)
{
    if(!track || !event) return kFALSE;
    if(dca==0 || cov ==0 ||XYZatDCA ==0 ) return kFALSE;
    AliExternalTrackParam etp;
    etp.CopyFromVTrack((AliVTrack*)track);

    const Double_t kBeampiperadius=3;  //maximal dca used for track propagation
    fEventVertex = RemoveDaughtersFromPrimaryVtx(track);
    if(!fEventVertex)return kFALSE;
    //printf("After remove:\n");
    //fEventVertex->Print();
    const AliVVertex *vtxESDSkip =fEventVertex;//RemoveDaughtersFromPrimaryVtx(track);
    if(!vtxESDSkip) return kFALSE;
    if(!etp.PropagateTo(1,event->GetMagneticField())) return kFALSE;
    Bool_t success = kFALSE;
    Double_t x_at_dca[3];
    Double_t p_at_dca[3];

    //Classic dca calculation
    if(etp.PropagateToDCA(vtxESDSkip, event->GetMagneticField(), kBeampiperadius, dca, cov)){
        success = kTRUE;
        etp.GetXYZ(XYZatDCA);
        etp.GetXYZ(x_at_dca);
        etp.GetPxPyPz(p_at_dca);
            //         if(fIsPythia)   dca[0] *= fMCglobalDCASmear;
            //         if(fIsPythia)   dca[0] += fMCglobalDCAxyShift; // generic mean offset in LHC10e default is 0.007 == 7 µm

    } else return kFALSE;
    return success;
}

/*AliExternalTrackParam AliAnalysisTaskHFJetIPQA::GetExternalParamFromJet(const AliEmcalJet *jet, const AliAODEvent *event)
{
    double vtx[3]= {0.};
    double cov [21] = {0.};
    double pxpypz[3] = {0.};
    jet->PxPyPz(pxpypz);
    (event->GetPrimaryVertex())->GetXYZ(vtx);
    AliExternalTrackParam etp     (vtx, pxpypz, cov, (Short_t)0);
    return etp;
}*/



Bool_t AliAnalysisTaskHFJetIPQA::getJetVtxMass(AliEmcalJet *jet,double &value ){

    return true;



    const AliVVertex * vertex = InputEvent()->GetPrimaryVertex();

    Printf("Primary vertex %e %e %e",vertex->GetX(),vertex->GetY(),vertex->GetY() );
    /*
    if(!amvf)  amvf = new AliHFAdaptiveMVF(InputEvent()->GetMagneticField());
    amvf->_fitter->_seeder->_vertex_penalty((AliAODVertex*)InputEvent()->GetPrimaryVertex());
    amvf->_fitter->_r_fitter_jet(jet,(AliAODEvent*)InputEvent());
    return kTRUE;*/
}


Bool_t AliAnalysisTaskHFJetIPQA:: GetImpactParameterWrtToJet(const AliAODTrack *track, const AliAODEvent *event, const AliEmcalJet *jet, Double_t *dca, Double_t *cov, Double_t *XYZatDCA, Double_t &jetsign, int jetflavour)
{
    if(!track || !event || !jet)return kFALSE;
    if(dca==0 || cov ==0 ||XYZatDCA ==0 ) return kFALSE;

    if(!GetImpactParameter(track,event,dca,cov,XYZatDCA)) return kFALSE;
    //vertex properties
    Double_t VxVyVz[3]= {0.,0.,0.};
    const  AliVVertex *vtxESDSkip =fEventVertex;    //GetImpactParameter does set fEventVertex to the recalculated one
    if(!fEventVertex) return kFALSE;
    vtxESDSkip->GetXYZ(VxVyVz);
    //printf("Vertex in wrt jet:\n");
    //vtxESDSkip->Print();

    //jet properties
    double jetp[3];
    jet->PxPyPz(jetp);
    TVector3 jetP3(jetp);
    double covjet [21] = {0.};
    double pxpypz[3] = {0.};
    jet->PxPyPz(pxpypz);
    AliExternalTrackParam etp_jet(VxVyVz, pxpypz, covjet, (Short_t)0);

    //Calculation of sign
    jetsign=GetIPSign(XYZatDCA, jetp,VxVyVz);

    //track properties
    AliExternalTrackParam etp_track;    etp_track.CopyFromVTrack(track);
    Double_t xa,xb,xyz_jet_global[3],xyz_track_global[3];

    etp_jet.GetDCA(&etp_track, event->GetMagneticField(), xa, xb);
    etp_jet.GetXYZAt(xa, event->GetMagneticField(),xyz_jet_global);
    etp_track.GetXYZAt(xb, event->GetMagneticField(),xyz_track_global);
    if(!etp_track.PropagateTo(xb,event->GetMagneticField())) return kFALSE;

    if(fEventVertex) {
        delete fEventVertex;
        fEventVertex =nullptr;

    }

    double val = ((VxVyVz[0] - xyz_track_global[0]) * (VxVyVz[0] - xyz_track_global[0]) +
        (VxVyVz[1] - xyz_track_global[1]) * (VxVyVz[1] - xyz_track_global[1])+
        (VxVyVz[2] - xyz_track_global[2]) * (VxVyVz[2] - xyz_track_global[2]));


    double  bdecaylength = val >0 ? sqrt(val) : 1000;
    //printf("decaylength:\n");
    //for(int i=0;i<3;i++){
    //    printf("VxVyVZ=%f, xyztrackglobal=%f\n",VxVyVz[i], xyz_track_global[i]);
    //}

    double dcajetrack = sqrt((xyz_jet_global[0] - xyz_track_global[0]) * (xyz_jet_global[0] - xyz_track_global[0]) +
        (xyz_jet_global[1] - xyz_track_global[1]) * (xyz_jet_global[1] - xyz_track_global[1])+
        (xyz_jet_global[2] - xyz_track_global[2]) * (xyz_jet_global[2]- xyz_track_global[2]));

    //printf("decaylength:\n");
    //for(int i=0;i<3;i++){
    //    printf("xyzjetglobal=%f, xyztrackglobal=%f\n",xyz_jet_global[i], xyz_track_global[i]);
    //}

    if(!(IsDCAAccepted(bdecaylength, dcajetrack, dca, jetflavour))) return kFALSE;
    //printf("decaylength=%f, pi=%f\n", bdecaylength, dcajetrack);

    return kTRUE;
}


void AliAnalysisTaskHFJetIPQA::SetTCThresholds(TObjArray** &threshs){
  for(int iProbSet=0;iProbSet<fNThresholds;iProbSet++){
    TObjArray* oa=(TObjArray*)threshs[iProbSet];

    printf("Pointer oa=%p\n",oa);

    h1DThresholdsFirst.push_back((TH1D*)oa->At(0));
    h1DThresholdsSecond.push_back((TH1D*)oa->At(1));
    h1DThresholdsThird.push_back((TH1D*)oa->At(2));

    TString sFrac=h1DThresholdsFirst[iProbSet]->GetTitle();
    fFracs.push_back(sFrac.Atof());
  }
  //checking fFracs
  //for(int iFrac=0;iFrac<fFracs.size();iFrac++){
  //    printf("iFrac=%i, %f\n",iFrac,fFracs[iFrac]);
  //}
  /*  int nPoints=h1DThresholdsFirst[0]->GetNbinsX();

    printf("LargestIP: 2nd Histogram bins:\n");
    for(int iPoint=0;iPoint<nPoints;iPoint++){
      printf("iPoint=%i, xval=%f, yval=%f\n",iPoint, h1DThresholdsFirst[1]->GetBinCenter(iPoint),h1DThresholdsFirst[1]->GetBinContent(iPoint));
    }
    printf("2ndLargestIP: 3rd Histogram bins:\n");
    for(int iPoint=0;iPoint<nPoints;iPoint++)
      printf("iPoint=%i, xval=%f, yval=%f\n",iPoint, h1DThresholdsSecond[2]->GetBinCenter(iPoint),h1DThresholdsSecond[2]->GetBinContent(iPoint));
    }
    printf("3rdLargestIP: 1st Histogram bins:\n");
    for(int iPoint=0;iPoint<nPoints;iPoint++){
      printf("iPoint=%i, xval=%f, yval=%f\n",iPoint, h1DThresholdsThird[0]->GetBinCenter(iPoint),h1DThresholdsThird[0]->GetBinContent(iPoint));
    }*/
}

void AliAnalysisTaskHFJetIPQA::SetProbThresholds(TObjArray** &threshs){
  for(int iProbSet=0;iProbSet<fNThresholds;iProbSet++){
    TObjArray* oa=(TObjArray*)threshs[iProbSet];
    if(!oa) AliError(Form(" No %i'th Probability Threshold object array!\n",iProbSet));
    printf("Pointer oa=%p\n",oa);
    h1DProbThresholds.push_back((TH1D*)oa->At(0));
    if(!h1DProbThresholds.back()) AliError(Form(" No %i'th Probability Threshold hist!\n",iProbSet));
  }

  /*int nPoints=h1DProbThresholds[0]->GetNbinsX();
  for(int iPoint=0;iPoint<nPoints;iPoint++){
    printf("iPoint=%i, xval=%f, yval=%f\n",iPoint, h1DProbThresholds[0]->GetXaxis()->GetBinLowEdge(iPoint),h1DProbThresholds[0]->GetBinContent(iPoint));
  }*/
}

void AliAnalysisTaskHFJetIPQA::SetTagSettings(int iTagSetting, int ntracktypes){
  this->setDoTCTagging(iTagSetting);

  switch(fDoTCTagging){
    case TCIPSigPtDep:
      fNTrackTypes=3;
      fUseSignificance=kTRUE;
      break;

    case TCIPFixedPt:
      fNTrackTypes=1;
      fUseSignificance=kFALSE;
      break;

    case TCIPSigPtDepVarNTemps:
      fNTrackTypes=ntracktypes;
      fUseSignificance=kTRUE;
      break;
  }
}

// Read Threshold Histograms
//==============================================================================
void AliAnalysisTaskHFJetIPQA::ReadThresholdHists(TString PathToThresholds, TString taskname, int nTCThresh, int iTagSetting, int ntracktypes){
    TFile* fileThresholds=TFile::Open(PathToThresholds.Data());
    if(!fileThresholds ||(fileThresholds&& !fileThresholds->IsOpen())){AliError(Form("%s :: File with threshold values not found",taskname.Data()));}

    Printf("%s :: File %s successfully loaded, setting up threshold functions for %i q-values.",taskname.Data(),PathToThresholds.Data(),nTCThresh);

    SetTagSettings(iTagSetting, ntracktypes);

    if(fileThresholds){
        printf("Reading threshold histograms for track counting...\n");

        //TC Thresholds
        TObjArray** oaTCThresh=new TObjArray*[nTCThresh];
        for(int iThresh=0;iThresh<nTCThresh;iThresh++){
          fileThresholds->GetObject(Form("TCThres_%i",iThresh),oaTCThresh[iThresh]);
        }

        //ProbLookup hists
        TObjArray* oLookup;
        fileThresholds->GetObject("ProbLookup",oLookup);

       /* TObjArray** oaProbThresh=new TObjArray*[nTCThresh];
        for(int iThresh=0;iThresh<nTCThresh;iThresh++){
          fileThresholds->GetObject(Form("ProbThres_%i",iThresh),oaProbThresh[iThresh]);
        }*/

        this->setfNThresholds(nTCThresh);
        this->SetTCThresholds(oaTCThresh);
        //this->SetProbThresholds(oaProbThresh);
        this->ReadProbvsIPLookup(oLookup);
    }
}

void AliAnalysisTaskHFJetIPQA::ReadProbvsIPLookup(TObjArray*& oLookup){
  for(int iN=0;iN<fNTrackTypes;iN++){
    h2DProbLookup.push_back((TH2D*)oLookup->At(iN));
  }

  /*for(int iN=0;iN<fNTrackTypes;iN++){
    printf("Reading hist for %i, %p\n", iN, h2DProbLookup[iN]);
    int nBinsX=h2DProbLookup[iN]->GetNbinsX();
    int nBinsY=h2DProbLookup[iN]->GetNbinsY();
    for(int iBinsX=1;iBinsX<=nBinsX;iBinsX++){
      for(int iBinsY=1;iBinsY<=nBinsY;iBinsY++){
         printf("iBinx=%i, iBiny=%i, %f\n",iBinsX, iBinsY, h2DProbLookup[iN]->GetBinContent(iBinsX,iBinsY));
      }
    }
  }*/
}

void AliAnalysisTaskHFJetIPQA::DoTCTagging(Float_t jetpt, Int_t nGoodIPTracks, const vector<Float_t>& ipval, Bool_t **kTagDec){
  //printf("Start TCTagging!\n");
  //threshold values for tracks with largest, second and third largest IP
  //always takes larger threshold if jetpt = binboundary
  int iJetPtBin=h1DThresholdsFirst[0]->FindBin(jetpt);
  double IPthresN1[fNThresholds];  //IP threshold values for individual separation power
  double IPthresN2[fNThresholds];
  double IPthresN3[fNThresholds];

  if(fDoTCTagging==TCIPSigPtDep){
    for(int iFrac=0;iFrac<fNThresholds;iFrac++){
      IPthresN1[iFrac]=h1DThresholdsFirst[iFrac]->GetBinContent(iJetPtBin);
      IPthresN2[iFrac]=h1DThresholdsSecond[iFrac]->GetBinContent(iJetPtBin);
      IPthresN3[iFrac]=h1DThresholdsThird[iFrac]->GetBinContent(iJetPtBin);
    }
  }
  if(fDoTCTagging==TCIPFixedPt){
    for(int iFrac=0;iFrac<fNThresholds;iFrac++){
      IPthresN1[iFrac]=fTCThresholdPtFixed;
      IPthresN2[iFrac]=fTCThresholdPtFixed;
      IPthresN3[iFrac]=fTCThresholdPtFixed;
    }
  }
  //printf("Starting DoTCTagging with nGoodIPTracks=%i\n",nGoodIPTracks);

  for(int iThresh=0;iThresh<fNThresholds;iThresh++){
    if(nGoodIPTracks==0) continue;
    //printf("DoTCTagging:\n");
    //printf("      iJetPtBin=%i, IPthresN1=%f, IPthresN2=%f, IPthresN3=%f\n", iJetPtBin, IPthresN1[iThresh],IPthresN2[iThresh], IPthresN3[iThresh]);


    if(nGoodIPTracks>2){
      //tripple tag
      //printf("ipval[0]=%f, ipval[1]=%f, ipval[2]=%f\n", ipval[0],ipval[1],ipval[2]);
      if(ipval[0]>IPthresN1[iThresh]&&ipval[1]>IPthresN2[iThresh]&&ipval[2]>IPthresN3[iThresh]){
          //printf("Triple %f!\n",fFracs[iThresh]);
          kTagDec[iThresh][Full]=kTRUE; kTagDec[iThresh][Triple]=kTRUE;
      }

      //single tag
      if(kTagLevel<2){
        //printf("Single catch\n");
        if(ipval[2]>IPthresN3[iThresh]) {
        //  printf("Single3rd %f!\n",fFracs[iThresh]);
          kTagDec[iThresh][Full]=kTRUE; kTagDec[iThresh][Single3rd]=kTRUE;
        }
      }
    }

    if(nGoodIPTracks>1){
      //printf("ipval[0]=%f, ipval[1]=%f\n", ipval[0],ipval[1]);
      //double tag
      if(kTagLevel<3){
        //printf("Double catch\n");
        if(ipval[0]>IPthresN1[iThresh]&&ipval[1]>IPthresN2[iThresh]) {
          //printf("Double %f!\n",fFracs[iThresh]);
          kTagDec[iThresh][Full]=kTRUE; kTagDec[iThresh][Double]=kTRUE;
        }
      }
      //single tag
      if(kTagLevel<2){
        //printf("Single catch\n");
        if(ipval[1]>IPthresN2[iThresh]) {
          //printf("Single2nd %f!\n",fFracs[iThresh]);
          kTagDec[iThresh][Full]=kTRUE; kTagDec[iThresh][Single2nd]=kTRUE;
        }
      }
    }

    //single tag
    if(nGoodIPTracks>0){
      //printf("ipval[0]=%f", ipval[0]);
      if(kTagLevel<2){
        //printf("Single catch\n");
        if(ipval[0]>IPthresN1[iThresh]) {
          //printf("Single1st %f, %f!\n",fFracs[iThresh], IPthresN1[iThresh]);
          kTagDec[iThresh][Full]=kTRUE; kTagDec[iThresh][Single1st]=kTRUE;
        }
      }
    }
    /*printf("Testing kTagLevel:: Line %i\n ",__LINE__);
    for(int iThresh=0;iThresh<fNThresholds;iThresh++){
      for(int iType=0;iType<6;iType++){
        printf("iThresh=%f, %i, kTagDec=%i\n",fFracs[iThresh],iType,kTagDec[iThresh][iType]);
      }
    }*/
    //bFull[iThresh]=kTagDec[iThresh][Full];
    bSingle1st[iThresh]=kTagDec[iThresh][Single1st];
    //bSingle2nd[iThresh]=kTagDec[iThresh][Single2nd];
    //bSingle3rd[iThresh]=kTagDec[iThresh][Single3rd];
    bDouble[iThresh]=kTagDec[iThresh][Double];
    //bTriple[iThresh]=kTagDec[iThresh][Triple];
  }
}

void AliAnalysisTaskHFJetIPQA::DoProbTagging(double probval, double jetpt,bool **kTagDec){
  int iJetPtBin=h1DProbThresholds[0]->FindBin(jetpt);
  for(int iThresh=0;iThresh<fNThresholds;iThresh++){
    double threshval=h1DProbThresholds[iThresh]->GetBinContent(iJetPtBin);
    //printf(" iThres=%i, iJetPtBin=%i, jetpt=%f, iThreshold=%f, probval=%f", iThresh, iJetPtBin, jetpt, threshval, probval);

    if(probval>threshval){
        kTagDec[iThresh][Full]=kTRUE;
        //printf("Tagging condition fullfilled %i!\n",jetflavour);
    }
  }
}

void AliAnalysisTaskHFJetIPQA::FillTaggedJetPtDistribution(bool** kTagDec, double jetpt){
    const char * tagtype[6] = {"Full","Single1st","Single2nd","Single3rd","Double","Triple"};

    for(int iThresh=0;iThresh<fNThresholds;iThresh++){
      for(int iType=0;iType<6;iType++){
        if(kTagDec[iThresh][iType]){
          //printf("Filling fracs=%f, iType=%s\n",fFracs[iThresh],tagtype[iType]);
          FillHist(Form("h1DTaggedJetPt_%s_%0.2f",tagtype[iType],fFracs[iThresh]),jetpt,1);
        }
      }
    }
}

void AliAnalysisTaskHFJetIPQA::FillV0EfficiencyHists(int isV0, int &jetflavour, double jetpt, bool &isV0Jet){
  //printf("isV0=%i, jetflavour=%i, jetpt=%f\n",isV0, jetflavour, jetpt);
    switch (isV0){
      case V0Rec:
          FillHist(Form("h1DV0FalseRec"), jetpt, 1);
          //printf("Found false Tag\n");
          break;

      case V0MC:
        FillHist(Form("h1DV0TrueDataDef"), jetpt, 1);
        if((jetflavour!=B)){
          isV0Jet=kTRUE;
          FillHist(Form("h1DV0TrueMCDef"), jetpt, 1);
          //printf("Found MC def true V0\n");
        }
        //printf("Found MC true V0 Jet: jetflavour=%i\n",jetflavour);
        if(jetflavour==UDSG) jetflavour=UDSGV0;
        if(jetflavour==C) jetflavour=CV0;
          break;

      case V0TrueRec:
        FillHist(Form("h1DV0TrueDataDef"), jetpt, 1);
        FillHist(Form("h1DV0TrueRec"), jetpt, 1);
        if(jetflavour==UDSG) jetflavour=UDSGV0;
        if(jetflavour==C) jetflavour=CV0;
        if(jetflavour!=B){
          FillHist(Form("h1DV0TrueMCDef"), jetpt, 1);
          isV0Jet=kTRUE;
          //printf("Found MC def true V0\n");
        }
        //printf("Found Rec true V0 Jet: jetflavour=%i\n",jetflavour);
          break;
    }
}

void AliAnalysisTaskHFJetIPQA::GetLowUpperBinNo(int &iLowerBin, int &iUpperBin, double min, double max, TString type, Int_t iN){
    double fLowLowerBound=999;
    double fLowUpperBound=999;
    double fUpUpperBound=999;
    if(type.Contains("x")){
      iLowerBin=h2DProbLookup[iN]->GetXaxis()->FindBin(min);
      iUpperBin=h2DProbLookup[iN]->GetXaxis()->FindBin(max);
      fLowLowerBound=h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iLowerBin);
      fLowUpperBound=h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iUpperBin);
      fUpUpperBound=h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iUpperBin+1);
    }
    if(type.Contains("y")){
      iLowerBin=h2DProbLookup[iN]->GetYaxis()->FindBin(min);
      iUpperBin=h2DProbLookup[iN]->GetYaxis()->FindBin(max);
      fLowLowerBound=h2DProbLookup[iN]->GetYaxis()->GetBinLowEdge(iLowerBin);
      fLowUpperBound=h2DProbLookup[iN]->GetYaxis()->GetBinLowEdge(iUpperBin);
      fUpUpperBound=h2DProbLookup[iN]->GetYaxis()->GetBinLowEdge(iUpperBin+1);
    }

    //printf("iLowerBin=%i, iUpperBin=%i, fLowUpperBound=%f, fUpUpperBound=%f, fLowLowerbound=%f\n",iLowerBin, iUpperBin, fLowUpperBound, fUpUpperBound,fLowLowerBound);
    if(max==fLowUpperBound){
        iUpperBin=iUpperBin-1;
        //printf("Decreasing upper bound: iUpperBin=%i\n",iUpperBin);
    }
    else{
        //printf("WARNING: Extending upper limit from %f to bin boundary %f\n",max, fUpUpperBound);
    }
    if(min!=fLowLowerBound) printf("WARNING: Extending lower limit from %f to bin boundary %f\n",min, fLowLowerBound);
}

Float_t AliAnalysisTaskHFJetIPQA::IntegrateIP(Float_t jetpt, Float_t IP, Int_t iN){
  Int_t iZeroIPBin=-99;//=h2DProbLookup[iN]->GetXaxis()->FindBin(0.);
  Int_t iStartIPBin=-99;//=h2DProbLookup[iN]->GetXaxis()->FindBin(-25);
  Int_t iJetPtBin=-99;
  Int_t iIPBin=-99;

  //expand upper and lower values to bin boundaries. if value=bin boundary: expand to higher for upper limit and to lower for lower limit
  GetLowUpperBinNo(iStartIPBin, iZeroIPBin, -25,0,"x",iN);
  GetLowUpperBinNo(iStartIPBin, iIPBin, -25, -IP,"x",iN);
  //if jetpt=binboundary->take higher bin
  iJetPtBin=h2DProbLookup[iN]->GetYaxis()->FindBin(jetpt);

  Float_t probnomi=h2DProbLookup[iN]->Integral(iStartIPBin,iIPBin,iJetPtBin,iJetPtBin);
  Float_t probdenomi=h2DProbLookup[iN]->Integral(iStartIPBin,iZeroIPBin,iJetPtBin,iJetPtBin);
  //printf("probnomi=%f, probdenomi=%f\n",probnomi, probdenomi);
  Float_t prob=probnomi/probdenomi;
  //printf("Integrate: Zero=%f, StartIP(-25)=%f, IPValue=%f, lowy=%f, upy=%f, prob=%f\n", h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iZeroIPBin), h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iStartIPBin),h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iIPBin),h2DProbLookup[iN]->GetYaxis()->GetBinLowEdge(iJetPtBin),h2DProbLookup[iN]->GetYaxis()->GetBinLowEdge(iJetPtBin+1),prob);

  return prob;
}

/* Calculates track probability
 *
 * - takes lookup TH2D where jetpt vs. IP is stored; integrates it according to IntegrateIP()
 * - discards all tracks with negative IP and IP > then reach of fit functions for IP templates
 * - if number of tracks within jet larger than available amount of templates -> take template with largest tracknumber
 * */

Float_t AliAnalysisTaskHFJetIPQA::GetTrackProbability(Float_t jetpt, Int_t nGoodIPTracks, const vector<Float_t>& ipval){
  Float_t prob=1;
  Float_t probval=-999;
  Int_t nIPTracksAboveZero=0;

  for(long unsigned iN=0;iN<ipval.size();iN++){
    if(ipval[iN]<0) continue;
    if(ipval[iN]>fAnalysisCuts[bAnalysisCut_MaxIPLNJP]){
      continue;
    }
    probval=-999;
    if((int)iN>=fNTrackTypes){ probval=IntegrateIP(jetpt,ipval[iN], fNTrackTypes-1);}    //probval[iN]=h2DProbLookup[iN]->GetBinContent(iIPBin[iN],iJetPtBin);}
    else{probval=IntegrateIP(jetpt,ipval[iN], iN);}    //probval[iN]=h2DProbLookup[iN]->GetBinContent(iIPBin[iN],iJetPtBin);}
    //printf("iN=%i: jetpt=%f, ipval[%i]=%f, prob=%f\n",iN, jetpt, iN, ipval[iN],probval);
    fTrackProb[iN]=probval;
    prob=prob*probval;
    nIPTracksAboveZero++;
  }

  Float_t LNExpo=1;
  Float_t Faculty=1;
  Float_t LNFunc=TMath::Log(prob);
  Float_t JP=0;

  //printf("nIPTracksAboveZero=%i, prob=%f, LNExpo=%f, Faculty=%f, LNFunc=%f, JP=%f\n", nIPTracksAboveZero, prob, LNExpo, Faculty, LNFunc, JP);
  if(nIPTracksAboveZero>0)JP=prob;

  for(int iN=1;iN<nIPTracksAboveZero;iN++){
    LNExpo*=-LNFunc;
    Faculty*=iN;
    JP+=prob*LNExpo/Faculty;
    //printf("LNExpo=%f, Faculty=%f, JP=%f\n", LNExpo, Faculty, JP);
  }

  return JP;
}


void AliAnalysisTaskHFJetIPQA::Terminate(Option_t *){

    printf("\n*********************************\n");
    printf("Corrections:\n");
    printf("    MC Corrections (Data/MC+Fluka):%i\n",fDoMCCorrection);
    printf("    Track Smearing:%i\n",fRunSmearing );
    printf("    Underlying Event Subtraction:%i\n", fDoUnderlyingEventSub);
    printf("*********************************\n");
}
