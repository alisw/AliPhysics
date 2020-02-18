#include "TList.h"
#include "TMatrixD.h"
#include "TParticle.h"
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
#include "AliStack.h"
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

AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA():
AliAnalysisTaskEmcalJet(),
fEventCuts(0),
fHistManager(),
fEventVertex(nullptr),
fPidResponse(nullptr),
jetconrec(nullptr),
jetcongen(nullptr),
fRunSmearing(kFALSE),
fUsePIDJetProb(kFALSE),
fDoMCCorrection(kFALSE),
fDoUnderlyingEventSub(kFALSE),
fApplyV0Rej(kFALSE),
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
kTagLevel(3),
fFracs(0),
fXsectionWeightingFactor(1),
fProductionNumberPtHard(-1),
fNThresholds(1),
sTemplateFlavour(0),
fJetRadius(0.4),
fDaughtersRadius(1),
fNoJetConstituents(0),
fTCThresholdPtFixed(0.008),
fGraphMean(nullptr),
fGraphSigmaData(nullptr),
fGraphSigmaMC(nullptr),
fGraphXi(nullptr),
fGraphOmega(nullptr),
fK0Star(nullptr),
fPhi(nullptr),
fGeant3FlukaProton(nullptr),
fGeant3FlukaAntiProton(nullptr),
fGeant3FlukaLambda(nullptr),
fGeant3FlukaAntiLambda(nullptr),
fGeant3FlukaKMinus(nullptr),
h1DThresholdsFirst(0),
h1DThresholdsSecond(0),
h1DThresholdsThird(0),
h2DProbLookup(0),
h2DProbDistsUnid(0),
h2DProbDistsudsg(0),
h2DProbDistsc(0),
h2DProbDistsb(0),
h2DProbDistsudsgV0(0),
h2DProbDistscV0(0),
h2DLNProbDistsUnid(0),
h2DLNProbDistsudsg(0),
h2DLNProbDistsc(0),
h2DLNProbDistsb(0),
h2DLNProbDistsudsgV0(0),
h2DLNProbDistscV0(0),
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
    DefineOutput(1,  AliEmcalList::Class()) ;
}
AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA(const char *name):
AliAnalysisTaskEmcalJet(name, kTRUE),
fEventCuts(0),
fHistManager(name),
fEventVertex(nullptr),
fPidResponse(nullptr),
jetconrec(nullptr),
jetcongen(nullptr),
fRunSmearing(kFALSE),
fUsePIDJetProb(kFALSE),
fDoMCCorrection(kFALSE),
fDoUnderlyingEventSub(kFALSE),
fApplyV0Rej(kFALSE),
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
kTagLevel(3),
fFracs(0),
fXsectionWeightingFactor(1.),
fProductionNumberPtHard(-1),
fNThresholds(1),
sTemplateFlavour(0),
fJetRadius(0.4),
fDaughtersRadius(1),
fNoJetConstituents(0),
fTCThresholdPtFixed(0.008),
fGraphMean(nullptr),
fGraphSigmaData(nullptr),
fGraphSigmaMC(nullptr),
fGraphXi(nullptr),
fGraphOmega(nullptr),
fK0Star(nullptr),
fPhi(nullptr),
fGeant3FlukaProton(nullptr),
fGeant3FlukaAntiProton(nullptr),
fGeant3FlukaLambda(nullptr),
fGeant3FlukaAntiLambda(nullptr),
fGeant3FlukaKMinus(nullptr),
h1DThresholdsFirst(0),
h1DThresholdsSecond(0),
h1DThresholdsThird(0),
h2DProbLookup(0),
h2DProbDistsUnid(0),
h2DProbDistsudsg(0),
h2DProbDistsc(0),
h2DProbDistsb(0),
h2DProbDistsudsgV0(0),
h2DProbDistscV0(0),
h2DLNProbDistsUnid(0),
h2DLNProbDistsudsg(0),
h2DLNProbDistsc(0),
h2DLNProbDistsb(0),
h2DLNProbDistsudsgV0(0),
h2DLNProbDistscV0(0),
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
    fAnalysisCuts[bAnalysisCut_HasTPCrefit]=1;
    fAnalysisCuts[bAnalysisCut_HasITSrefit]=1;
    fAnalysisCuts[bAnalysisCut_KinkCand]=1;

    //Jet Cuts
    fAnalysisCuts[bAnalysisCut_MinJetPt]        =0;  //only for settings output. Not really used as cuts are done in .C file
    fAnalysisCuts[bAnalysisCut_MaxJetPt]        =1000;
    fAnalysisCuts[bAnalysisCut_MinJetEta]       =-0.9;
    fAnalysisCuts[bAnalysisCut_MaxJetEta]       =0.9;

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
    fV0Cuts[MaxLifeTime]=5;
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


int AliAnalysisTaskHFJetIPQA::GetMCTruth(AliAODTrack * track, int &motherpdg){
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
}



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
    if(fIsPythia){
        FillHist("fh1dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        FillHist("fh1dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        FillHist("fh1dTracksImpParXYSignificance_McCorr",GetValImpactParameter(kXYSig,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        //FillHist("fh1dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        FillHist("fh2dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),track->Pt(),weight);     //*this->fXsectionWeightingFactor );
        FillHist("fh2dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),track->Pt(),weight);     //*this->fXsectionWeightingFactor );
        //FillHist("fh2dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),track->Pt(),weight);     //*this->fXsectionWeightingFactor );
    }
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

void AliAnalysisTaskHFJetIPQA::FillRecHistograms(int jetflavour, double jetpt, double eta, double phi){
  FillHist("fh1dJetRecPt",jetpt, 1);  //this->fXsectionWeightingFactor );
  FillHist("fh1dJetRecEtaPhiAccepted",eta,phi, 1);   //this->fXsectionWeightingFactor );
  FillHist("fh1dJetRecPtAccepted",jetpt, 1);  //this->fXsectionWeightingFactor );

  if(fIsPythia){
    if(jetflavour==0)     FillHist("fh1dJetRecPtUnidentified",jetpt, 1);    //this->fXsectionWeightingFactor );
      else if(jetflavour==1)FillHist("fh1dJetRecPtudsg",        jetpt, 1);    //this->fXsectionWeightingFactor );
      else if(jetflavour==2)FillHist("fh1dJetRecPtc",           jetpt, 1);    //this->fXsectionWeightingFactor );
      else if(jetflavour==3)FillHist("fh1dJetRecPtb",           jetpt, 1);    //this->fXsectionWeightingFactor );
      else if(jetflavour==4)FillHist("fh1dJetRecPts",           jetpt, 1);    //this->fXsectionWeightingFactor );
  }
}

void AliAnalysisTaskHFJetIPQA::FillGenHistograms(int jetflavour, AliEmcalJet* jetgen){
    FillHist("fh1dJetGenPt",GetPtCorrectedMC(jetgen), 1); //this->fXsectionWeightingFactor);
    if(jetflavour ==0)      FillHist("fh1dJetGenPtUnidentified",GetPtCorrectedMC(jetgen), 1); // this->fXsectionWeightingFactor );
    else if(jetflavour ==1) FillHist("fh1dJetGenPtudsg",GetPtCorrectedMC(jetgen), 1);   //this->fXsectionWeightingFactor );
    else if(jetflavour ==2) FillHist("fh1dJetGenPtc",GetPtCorrectedMC(jetgen), 1);  //this->fXsectionWeightingFactor );
    else if(jetflavour ==3) FillHist("fh1dJetGenPtb",GetPtCorrectedMC(jetgen), 1);  //this->fXsectionWeightingFactor );
    else if(jetflavour ==4) FillHist("fh1dJetGenPts",GetPtCorrectedMC(jetgen), 1);  //this->fXsectionWeightingFactor );*/
}


void AliAnalysisTaskHFJetIPQA::FillIPTypePtHists(int jetflavour, double jetpt, bool* nTracks){
    //Fill histograms for jets which have largest, second largest and third largest impact parameter
    //tracks passing the selection criterion

    for (Int_t iN = 1 ; iN <=3 ;++iN){
      if(!nTracks[iN]) continue;
      FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",iN,"all"),jetpt,1);     //*this->fXsectionWeightingFactor );

      if(jetflavour==0) continue;
      FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",iN,sTemplateFlavour[jetflavour].Data()),jetpt,1);     //*this->fXsectionWeightingFactor );
    }
}

void AliAnalysisTaskHFJetIPQA::FillIPTemplateHists(double jetpt, int iN,int jetflavour, double* params){
    const char * stype  [4] = {"fh2dJetSignedImpParXY","fh2dJetSignedImpParXYSignificance","fh2dJetSignedImpParXYZ","fh2dJetSignedImpParXYZSignificance"};
    const char * subord [3] = {"First","Second","Third"};

    for (Int_t iType = 0 ;iType <2 ;++iType){
        TString hname = Form("%s%s",stype[iType],subord[iN]);
        FillHist(hname.Data(),jetpt,params[iType],1);
        if(fIsPythia){
          TString hnameflav = Form("%s%s%s",stype[iType],sTemplateFlavour[jetflavour].Data(),subord[iN]);
          FillHist(hnameflav.Data(),jetpt,params[iType],1);
        }
    }
}

void AliAnalysisTaskHFJetIPQA::FillTrackIPvsPt(int isV0, double pt, double IP, int jetflavour){
  if(jetflavour==CV0||jetflavour==UDSGV0){
    FillHist("fh1dTracksIPvsPt_V0JetTracks", pt, IP,1);
    //printf("Filling fh1dTracksIPvsPt_V0JetTracks: isV0=%i, pt=%f, IP=%f, jetflavour=%i",isV0, pt, IP, jetflavour);
    if((isV0==V0MC)||(isV0==V0TrueRec)){
      FillHist("fh1dTracksIPvsPt_V0inV0Jet", pt, IP,1);
      //printf("Filling fh1dTracksIPvsPt_V0inV0Jet: isV0=%i, pt=%f, IP=%f, jetflavour=%i",isV0, pt, IP, jetflavour);
    }
  }
  if(jetflavour==B){
    FillHist("fh1dTracksIPvsPt_B", pt, IP,1);
    //printf("Filling fh1dTracksIPvsPt_B: isV0=%i, pt=%f, IP=%f, jetflavour=%i",isV0, pt, IP, jetflavour);
    if((isV0==V0MC)||(isV0==V0TrueRec)){
      FillHist("fh1dTracksIPvsPt_V0inBJet", pt, IP,1);
      //printf("Filling fh1dTracksIPvsPt_V0inBJet: isV0=%i, pt=%f, IP=%f, jetflavour=%i",isV0, pt, IP, jetflavour);
    }
  }
}


void AliAnalysisTaskHFJetIPQA::FillTrackTypeResHists(){
   printf("Filling track type resolution hists");

   /* if(GetImpactParameterWrtToJet((AliAODTrack*)trackV,(AliAODEvent*)InputEvent(),jetrec,dca,cov,xyzatcda,sign)){
      if(fEventVertex) {
        delete fEventVertex;
        fEventVertex =nullptr;
      }
      dca[0]=fabs(dca[0]);
      Double_t cursImParXYSig  =TMath::Abs(GetValImpactParameter(kXYSig,dca,cov))*sign;
      Double_t cursImParXYZSig =TMath::Abs(GetValImpactParameter(kXYZSig,dca,cov))*sign;

      Int_t corridx=-1;double ppt;
      (fIsPythia&&fDoMCCorrection) ? TrackWeight = GetMonteCarloCorrectionFactor(trackV,corridx,ppt) : TrackWeight =1;
                  Double_t cursImParXY     =TMath::Abs(GetValImpactParameter(   kXY,dca,cov))*sign;
                  Double_t cursImParXYZ    =TMath::Abs(GetValImpactParameter(   kXYZ,dca,cov))*sign;

                if(fIsPythia){
                  if(is_udgjet){
                      if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                          FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShits",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                      }
                      else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                          FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShits",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                      }
                      else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                          FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShits",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                      }
                      else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                          FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShits",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                      }
                  }
                  if(IsFromElectron((AliAODTrack*)trackV)){
                      if(is_udgjet){
                          if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                      }
                  }
                  else if(IsFromPion((AliAODTrack*)trackV)){
                      if(is_udgjet){
                          if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                      }
                  }
                  else if(IsFromKaon((AliAODTrack*)trackV)){
                      if(is_udgjet){
                          if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                      }
                  }
                  else if(IsFromProton((AliAODTrack*)trackV)){
                      if(is_udgjet){
                          if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                          else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                              FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight);     //this->fXsectionWeightingFactor );
                              //FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //this->fXsectionWeightingFactor );
                          }
                      }
                  }
                      //Fill jet probability ipsig histograms for template fitting
                  const char * subtype_jp [4] = {"","udsg","c","b"};*/
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

/*void AliAnalysisTaskHFJetIPQA::SetIPVals(vector <SJetIpPati > sImpPar, bool* hasIPs, double* ipval){
  if((int)sImpPar.size()>0) hasIPs[0]=kTRUE;
  if((int)sImpPar.size()>1) hasIPs[1]=kTRUE;
  if((int)sImpPar.size()>2) hasIPs[2]=kTRUE;

  if(hasIPs[0]){
    ipval[0] =sImpPar.at(0).first;
     //printf("HasIP0, ipval[0]=%f\n", ipval[0]);
  }
  if(hasIPs[1]){
    ipval[1] =sImpPar.at(1).first;
      //printf("HasIP1, ipval[1]=%f\n",ipval[1]);
  }
  if(hasIPs[2]){
    ipval[2] =sImpPar.at(2).first;
      //printf("HasIP2, ipval[2]=%f\n", ipval[2]);
  }
}*/

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

Bool_t AliAnalysisTaskHFJetIPQA::IsParticleInCone(const AliVParticle* part, const AliEmcalJet* jet, Double_t dRMax) {
// decides whether a particle is inside a jet cone
  if(!part || !jet) AliError(Form("Particle or Jet missing: part=%p, jet=%p\n", part,jet));

  TVector3 vecMom2(jet->Px(), jet->Py(), jet->Pz());
  TVector3 vecMom1(part->Px(), part->Py(), part->Pz());
  Double_t dR = vecMom2.DeltaR(vecMom1); // = sqrt(dEta*dEta+dPhi*dPhi)
  if(dR <= dRMax) return kTRUE;// momentum vectors of part1 and part2 are closer than dRMax
  return kFALSE;
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


void AliAnalysisTaskHFJetIPQA::GetV0MCTrueCandidates(AliAODEvent *fAODIn){

  AliAODMCHeader* headerMC = (AliAODMCHeader*)fAODIn->FindListObject(AliAODMCHeader::StdBranchName());

  if(!headerMC){
    AliError("No MC header found!");
  }
  // get position of the MC primary vertex
  Double_t dPrimVtxMCX = headerMC->GetVtxX();
  Double_t dPrimVtxMCY = headerMC->GetVtxY();
  Double_t dPrimVtxMCZ = headerMC->GetVtxZ();

  AliAODMCParticle *pAOD = 0;
  AliEmcalJet * jetMC  = 0x0;
  double fJetPt=0;

  for (Int_t i=0; i<fMCArray->GetEntriesFast(); i++) {
    pAOD = dynamic_cast<AliAODMCParticle*>(fMCArray->At(i));
    if (!pAOD) continue;

    // Get the distance between the production point of the MC V0 particle and the primary vertex
    Double_t dx = dPrimVtxMCX - pAOD->Xv();
    Double_t dy = dPrimVtxMCY - pAOD->Yv();
    Double_t dz = dPrimVtxMCZ - pAOD->Zv();
    Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
    Bool_t bV0MCIsPrimaryDist = (dDistPrimary < 0.01); // Is close enough to be considered primary-like?
    if(!bV0MCIsPrimaryDist) continue;

    //Ask for PDG
    Int_t id = pAOD->GetPdgCode();
    Bool_t bV0 = ((id==3122) || (id==-3122) || (id==310));
    if (!bV0) { pAOD=0; continue; }

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

    //acceptance cuts
    if (TMath::Abs(pAOD->Y()) > fV0Cuts[MaxV0Rap]) { pAOD=0; continue; }
    if (pAOD->Pt() <0.15) continue;
    if ((TMath::Abs(fPosEta) > fV0Cuts[DaughMaxEta])||(TMath::Abs(fNegEta) > fV0Cuts[DaughMaxEta])) continue;
    if ((fPosPt < fV0Cuts[DaughMinPt])||(fNegPt < fV0Cuts[DaughMinPt])) continue;


    if(id==310) {fh1dKshortPtMC->Fill(pAOD->Pt());}
    if(id==3122) { fh1dLamdaPtMC->Fill(pAOD->Pt());}
    if(id==-3122) { fh1dAnLamdaPtMC->Fill(pAOD->Pt());}

    if(!jetcongen) AliError("No MC jet container available!\n");
    jetcongen->ResetCurrentID();
    while ((jetMC = jetcongen->GetNextAcceptJet())){
      fJetPt= jetMC->Pt();
      //if(!(jetcongen->GetRhoParameter() == 0x0)){
      //  fJetPt = fJetPt - jetcongen->GetRhoVal() * jetMC->Area();
      //}
      if(fJetPt < 5.) continue;
      bool bIsInCone=IsParticleInCone(pAOD, jetMC, fJetRadius);
      if(!bIsInCone) continue;

      if(id==310) {
          fh2dKshortPtVsJetPtMC->Fill(pAOD->Pt(), fJetPt);
          //printf("Fount MCTrue K0s: id=%i, y=%f, pt=%f\n",id,pAOD->Y(),pAOD->Pt());
      }
      if(id==3122) {
          fh2dLamdaPtVsJetPtMC->Fill(pAOD->Pt(), fJetPt);
          //printf("Fount MCTrue Lambda: id=%i, y=%f, pt=%f\n",id,pAOD->Y(),pAOD->Pt());
      }
      if(id==-3122) {
          fh2dAnLamdaPtVsJetPtMC->Fill(pAOD->Pt(), fJetPt);
          //printf("Fount MCTrue ALambda: id=%i, y=%f, pt=%f\n",id,pAOD->Y(),pAOD->Pt());
      }
    }
  }
  //delete jetMC;
  jetMC=NULL;
  pAOD=NULL;
  //delete pAOD;
  headerMC=NULL;
  //delete headerMC;
}

//////////////////////////////////////////////////
Int_t AliAnalysisTaskHFJetIPQA::IsV0Daughter(const AliAODTrack* track){
        if(!track)return kFALSE;
        AliAODv0* v0aod = 0x0;
        int posid = -1;
        int negid = -1;
        int trackid = -1;

        if(!fV0CandidateArray) {return kFALSE;}

        const AliAODTrack* trackPos =0x0;
        const AliAODTrack* trackNeg =0x0; // negative daughter track
        AliAODMCParticle *v0MC=0x0;
        AliAODMCParticle *v0Daugh=0x0;
        int tracklabel=track->GetLabel();
        int posdaughlabel=-99;
        int negdaughlabel=-99;
        int iMCLabelMother=-99;
        int iMCPdgMother=-99;
        int iMCPdgDaughter=-99;
        int iV0Tag=V0No;

        //printf("--------------------------Start with %i candidated------------------------\n",fV0CandidateArray->GetEntriesFast());
        for(int i = 0; i < fV0CandidateArray->GetEntriesFast(); ++i) {
                v0aod = dynamic_cast<AliAODv0*>(fV0CandidateArray->At(i));
                if(!v0aod){
                  continue;
                }
                trackPos=(AliAODTrack*) (AliAODTrack*)v0aod->GetDaughter(0); // positive daughter track
                trackNeg=(AliAODTrack*)(AliAODTrack*)v0aod->GetDaughter(1); // positive daughter track
                posid = trackPos->GetID();
                negid = trackNeg->GetID();
                trackid = track->GetID();

                //test the matching of IDs and compare to MC truth
                if(fIsPythia){
                  posdaughlabel=trackPos->GetLabel();
                  negdaughlabel=trackNeg->GetLabel();

                  if((tracklabel==posdaughlabel)||(tracklabel==negdaughlabel)){
                    //printf("Matched track to daughter: tracklabel=%i, posdaughlabel=%i, negdaughlabel=%i, trackID=%i, posID=%i, negIP=%i\n", tracklabel, posdaughlabel, negdaughlabel, trackid, posid, negid);
                  }
                  if((tracklabel == posdaughlabel)&&(trackid!=posid) ) {
                    AliError(Form("Mismatch of trackid=%i and posdaughter id=%i! (tracklabel=%i, daughlabel=%i) Are you using hybrid tracks?",trackid, posid, tracklabel, posdaughlabel));
                  }
                  if((tracklabel == negdaughlabel)&&(trackid!=negid)){
                    AliError(Form("Mismatch of trackid=%i and negdaughter id=%i! (tracklabel=%i, daughlabel=%i) Are you using hybrid tracks?",trackid, negid, tracklabel, negdaughlabel));
                  }
                }

                if(posid == trackid || negid == trackid) {
                    //printf("Reject V0 candidate: posid=%i, negid=%i, trackid=%i\n", posid, negid, trackid);
                    iV0Tag=V0Rec;
                }
        }
        if(fIsPythia){
          if(tracklabel>(fMCArray->GetEntriesFast())||(tracklabel<0)) return 0;
          v0Daugh=dynamic_cast<AliAODMCParticle*>(fMCArray->At(tracklabel));
          if(!v0Daugh) return 0;
          iMCLabelMother=v0Daugh->GetMother();
          v0MC=dynamic_cast<AliAODMCParticle*>(fMCArray->At(iMCLabelMother));
          if(!v0MC){
            return 0;
          }
          iMCPdgMother=v0MC->GetPdgCode();
          iMCPdgDaughter=v0Daugh->GetPdgCode();
          //printf("TrackLabel=%i, iMCLabelMother=%i,  pdgduaghter=%i,pdgmother=%i\n",tracklabel,iMCLabelMother, iMCPdgDaughter,iMCPdgMother);
          //printf("Daughter print!\n");
          //v0Daugh->Print();
          //printf("Mother print!\n");
          //v0MC->Print();
          if(((iMCPdgMother==310)||(iMCPdgMother==3122)||(iMCPdgMother==-3122))&&(iV0Tag!=V0Rec)){
              iV0Tag=V0MC;
              //printf("V0MC\n");
          }
          if(((iMCPdgMother==310)||(iMCPdgMother==3122)||(iMCPdgMother==-3122))&&(iV0Tag==V0Rec)){
              iV0Tag=V0TrueRec;
              //printf("TrueRec!\n");
          }
        }
        //printf("--------------------------End candidates------------------------\n");
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
  if(!fMCArray) { AliError("No fMCArray"); return NULL;}
  Int_t nStack = fMCArray->GetEntriesFast();

  if((iLabel < 0) || (iLabel >= nStack)){
      //printf("Daugh not in array range: iLabel=%i\n", iLabel);
      return NULL;
  }

  AliAODMCParticle *mctrack =  dynamic_cast<AliAODMCParticle *>(fMCArray->At(iLabel));
  return mctrack;
}

int AliAnalysisTaskHFJetIPQA::GetV0MCVeto(AliAODEvent* fAODIn, AliAODv0* v0, bool bIsCandidateK0s,bool bIsCandidateLambda, bool bIsCandidateALambda){
  // PDG codes of used particles
  Int_t iPdgCodePion = 211;
  Int_t iPdgCodeProton = 2212;
  Int_t iPdgCodeK0s = 310;
  Int_t iPdgCodeLambda = 3122;

  AliAODMCHeader* headerMC = (AliAODMCHeader*)fAODIn->FindListObject(AliAODMCHeader::StdBranchName());
  Double_t dPrimVtxMCX = 0., dPrimVtxMCY = 0., dPrimVtxMCZ = 0.; // position of the MC primary vertex

  if(!headerMC){
    AliError("No MC header found!");
  }
  // get position of the MC primary vertex
  dPrimVtxMCX = headerMC->GetVtxX();
  dPrimVtxMCY = headerMC->GetVtxY();
  dPrimVtxMCZ = headerMC->GetVtxZ();

  Int_t iNTracksMC = fMCArray->GetEntriesFast();
  if(!(bIsCandidateK0s) && !(bIsCandidateLambda)  && !(bIsCandidateALambda)) return 0; // chosen candidates with any mass

  const AliAODTrack* postrack = (AliAODTrack*)v0->GetDaughter(0); // positive daughter track
  const AliAODTrack* negtrack = (AliAODTrack*)v0->GetDaughter(1); // positive daughter track
  Int_t iposLabel = postrack->GetLabel();
  Int_t inegLabel = negtrack->GetLabel();
  AliAODMCParticle* particleMCDaughterPos=(AliAODMCParticle*) GetMCTrack(iposLabel);
  AliAODMCParticle* particleMCDaughterNeg=(AliAODMCParticle*) GetMCTrack(inegLabel);
  if(!particleMCDaughterNeg || !particleMCDaughterPos) return 0;

  Int_t iPdgCodeDaughterPos = particleMCDaughterPos->GetPdgCode();
  Int_t iPdgCodeDaughterNeg = particleMCDaughterNeg->GetPdgCode();
  Int_t iIndexMotherPos = particleMCDaughterPos->GetMother();
  Int_t iIndexMotherNeg = particleMCDaughterNeg->GetMother();
  Double_t fPosEta=particleMCDaughterPos->Eta();
  Double_t fPosPt=particleMCDaughterPos->Pt();
  Double_t fNegEta=particleMCDaughterNeg->Eta();
  Double_t fNegPt=particleMCDaughterNeg->Pt();

  if((iIndexMotherNeg < 0) || (iIndexMotherNeg >= iNTracksMC) || (iIndexMotherPos < 0) || (iIndexMotherPos >= iNTracksMC)){
    //printf("Mother not in array range: iIndexMotherNeg=%i, iIndexMotherPos=%i,iNTracksMC=%i\n", iIndexMotherNeg, iIndexMotherPos, iNTracksMC);
    return 0;
  }
  AliAODMCParticle* particleMCMotherPos = (AliAODMCParticle*)fMCArray->At(iIndexMotherNeg);
  AliAODMCParticle* particleMCMotherNeg = (AliAODMCParticle*)fMCArray->At(iIndexMotherPos);
  int posidmother=particleMCMotherPos->GetPdgCode();
  int negidmother=particleMCMotherNeg->GetPdgCode();
  /*int posndaughs=particleMCMotherPos->GetNDaughters();
  int negndaughs=particleMCMotherNeg->GetNDaughters();
  int pos0daughs=((AliAODMCParticle*)fMCArray->At(particleMCMotherPos->GetDaughterLabel(0)))->GetPdgCode();
  int pos1daughs=((AliAODMCParticle*)fMCArray->At(particleMCMotherPos->GetDaughterLabel(1)))->GetPdgCode();
  int neg0daughs=((AliAODMCParticle*)fMCArray->At(particleMCMotherNeg->GetDaughterLabel(0)))->GetPdgCode();
  int neg1daughs=((AliAODMCParticle*)fMCArray->At(particleMCMotherNeg->GetDaughterLabel(1)))->GetPdgCode();

  int pos0daughslabel=((AliAODMCParticle*)fMCArray->At(particleMCMotherPos->GetDaughterLabel(0)))->GetLabel();
  int pos1daughslabel=((AliAODMCParticle*)fMCArray->At(particleMCMotherPos->GetDaughterLabel(1)))->GetLabel();
  int neg0daughslabel=((AliAODMCParticle*)fMCArray->At(particleMCMotherNeg->GetDaughterLabel(0)))->GetLabel();
  int neg1daughslabel=((AliAODMCParticle*)fMCArray->At(particleMCMotherNeg->GetDaughterLabel(1)))->GetLabel();*/

  if((posidmother != iPdgCodeK0s) && (TMath::Abs(posidmother) != iPdgCodeLambda)&&(negidmother != iPdgCodeK0s) && (TMath::Abs(negidmother) != iPdgCodeLambda)) return 0;

  if(iIndexMotherNeg != iIndexMotherPos){
    //printf("Mothers are different: iIndexMotherNeg=%i, iIndexMotherPos=%i, neglab=%i, poslab=%i, negid=%i, posid=%i, v0idneg=%i, v0idpos=%i\n", iIndexMotherNeg,iIndexMotherPos,inegLabel, iposLabel ,iPdgCodeDaughterNeg,iPdgCodeDaughterPos,negidmother,posidmother);
    //printf("posndaughs=%i, pos0daughs=%i, pos1daughs=%i, negndaughs=%i, neg0daughs=%i, neg1daughs=%i\n",posndaughs, pos0daughs, pos1daughs, negndaughs, neg0daughs, neg1daughs);
    //printf("pos0daughslabel=%i, pos1daughslabel=%i, neg0daughslabel=%i, neg1daughslabel=%i\n",pos0daughslabel, pos1daughslabel, neg0daughslabel, neg1daughslabel);
    return 0;
  }
  //_____________________
  // Mother Properties
  AliAODMCParticle* particleMCMother = (AliAODMCParticle*)fMCArray->At(iIndexMotherPos);
  if(!particleMCMother)return 0;

  Int_t iPdgCodeMother = particleMCMother->GetPdgCode();
  Double_t dPtV0Gen = particleMCMother->Pt();
  Double_t dRapV0Gen = particleMCMother->Y();

  // Acceptance Cuts
  if((TMath::Abs(dRapV0Gen) > fV0Cuts[MaxV0Rap])) return 0;
  if(dPtV0Gen <0.15) return 0;
  if((TMath::Abs(fPosEta) > fV0Cuts[DaughMaxEta])||(TMath::Abs(fNegEta) > fV0Cuts[DaughMaxEta])) return 0;
  if((fPosPt<fV0Cuts[DaughMinPt])||(fNegPt<fV0Cuts[DaughMinPt])) return 0;

  // Skip not interesting particles
  if((iPdgCodeMother != iPdgCodeK0s) && (TMath::Abs(iPdgCodeMother) != iPdgCodeLambda)) return 0;

  // Check identity of the MC mother particle and the decay channel
  Bool_t bV0MCIsK0s = ((iPdgCodeMother == iPdgCodeK0s) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodePion));
  Bool_t bV0MCIsLambda = ((iPdgCodeMother == +iPdgCodeLambda) && (iPdgCodeDaughterPos == +iPdgCodeProton) && (iPdgCodeDaughterNeg == -iPdgCodePion));
  Bool_t bV0MCIsALambda = ((iPdgCodeMother == -iPdgCodeLambda) && (iPdgCodeDaughterPos == +iPdgCodePion) && (iPdgCodeDaughterNeg == -iPdgCodeProton));

  // Get the distance between production point of the MC mother particle and the primary vertex
  Double_t dx = dPrimVtxMCX - particleMCMother->Xv();
  Double_t dy = dPrimVtxMCY - particleMCMother->Yv();
  Double_t dz = dPrimVtxMCZ - particleMCMother->Zv();
  Double_t dDistPrimary = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
  Bool_t bV0MCIsPrimaryDist = (dDistPrimary < 0.01); // Is close enough to be considered primary-like?

  // K0s
  if(bIsCandidateK0s){ // selected candidates with any mass
    if(bV0MCIsK0s && bV0MCIsPrimaryDist){ // well reconstructed candidates
      //printf("Found K0s iPdgCodeMother=%i, iPdgCodeDaughterPos=%i, iPdgCodeDaughterNeg=%i, y=%f, ptv0=%f, etadaughpos=%f, ptdaughpos=%f,etadaughneg=%f, ptdaughneg=%f\n",
      //iPdgCodeMother,iPdgCodeDaughterPos,iPdgCodeDaughterNeg, dRapV0Gen,dPtV0Gen,fPosEta, fPosPt, fNegEta,fNegPt);
      return 1;
    }
  }
  // Lambda
  if(bIsCandidateLambda){ // selected candidates with any mass
    if(bV0MCIsLambda && bV0MCIsPrimaryDist){ // well reconstructed candidates
      //printf("Found Lambda iPdgCodeMother=%i, iPdgCodeDaughterPos=%i, iPdgCodeDaughterNeg=%i, y=%f, ptv0=%f, etadaughpos=%f, ptdaughpos=%f,etadaughneg=%f, ptdaughneg=%f\n",
      //iPdgCodeMother,iPdgCodeDaughterPos,iPdgCodeDaughterNeg, dRapV0Gen,dPtV0Gen,fPosEta, fPosPt, fNegEta,fNegPt);
      return 2;
    }
  }

  // anti-Lambda
  if(bIsCandidateALambda){ // selected candidates with any mass
    if(bV0MCIsALambda && bV0MCIsPrimaryDist){ // well reconstructed candidates
      //printf("Found ALambda iPdgCodeMother=%i, iPdgCodeDaughterPos=%i, iPdgCodeDaughterNeg=%i, y=%f, ptv0=%f, etadaughpos=%f, ptdaughpos=%f,etadaughneg=%f, ptdaughneg=%f\n",
      //iPdgCodeMother,iPdgCodeDaughterPos,iPdgCodeDaughterNeg, dRapV0Gen,dPtV0Gen,fPosEta, fPosPt, fNegEta,fNegPt);
      return 3;
    }
  }

  return 0;
}

void AliAnalysisTaskHFJetIPQA::SelectV0Candidates(AliAODEvent *fAODIn){
    AliAODv0* v0 = 0; // pointer to V0 candidates
    SV0Cand* sV0=new SV0Cand();
    SV0Daugh* sPosDaugh=new SV0Daugh();
    SV0Daugh* sNegDaugh=new SV0Daugh();


    fV0CandidateArray->Delete();//Reset the TClonesArray

    // Mean lifetime
    Int_t iNV0s = fAODIn->GetNumberOfV0s(); // get the number of V0 candidates

    //printf("################## Select candidates ########################\n");

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
        if(!sPosDaugh->bTPCRefitOn) continue;// TPC refit is ON?
        if(!sNegDaugh->bTPCRefitOn)  continue;
      }

      if(fV0Cuts[IsKinkCand]){
        if(sPosDaugh->bIsKink) continue;// kink daughter rejection
        if(sNegDaugh->bIsKink) continue;
      }

      if(fV0Cuts[DoPosNoTPCClusters]){
        if(sPosDaugh->iNoTPCCluster <= 0.)  continue;
        if(sNegDaugh->iNoTPCCluster <= 0.)  continue;
      }

      if(fV0Cuts[MinNoCrossedTPCRows] > 0.){
        if(sPosDaugh->iCrossedTPC < fV0Cuts[MinNoCrossedTPCRows]) continue;// Crossed TPC padrows
        if(sNegDaugh->iCrossedTPC < fV0Cuts[MinNoCrossedTPCRows]) continue;
      }

      if(fV0Cuts[NoCrossedOverNoTPCClustersMin] > 0.){
        if(sPosDaugh->iCrossedTPC / sPosDaugh->iNoTPCCluster < fV0Cuts[NoCrossedOverNoTPCClustersMin]) continue;
        if(sNegDaugh->iCrossedTPC / sNegDaugh->iNoTPCCluster < fV0Cuts[NoCrossedOverNoTPCClustersMin]) continue;
      }

      if(fV0Cuts[NoCrossedOverNoTPCClustersMax] > 0.){
        if(sPosDaugh->iCrossedTPC / sPosDaugh->iNoTPCCluster > fV0Cuts[NoCrossedOverNoTPCClustersMax])   continue;
        if(sNegDaugh->iCrossedTPC / sNegDaugh->iNoTPCCluster > fV0Cuts[NoCrossedOverNoTPCClustersMax])   continue;
      }

      FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      iCutIndex++;

      // 4
      // Daughters: transverse momentum cut
      if(fV0Cuts[DaughMinPt] > 0.){
        if((sPosDaugh->fPt < fV0Cuts[DaughMinPt]) || (sNegDaugh->fPt < fV0Cuts[DaughMinPt]))   continue;
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 5
      // Daughters: Impact parameter of daughters to prim vtx
      if(fV0Cuts[MinDCADaughWrtPV] > 0.){
        if((sPosDaugh->fDCAtoPV < fV0Cuts[MinDCADaughWrtPV] ) || (sNegDaugh->fDCAtoPV < fV0Cuts[MinDCADaughWrtPV] ))   continue;
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 6
      // Daughters: DCA
      if(fV0Cuts[MaxDCADaughvsDaugh] > 0.){
        if(sV0->fDCAV0DaughvsDaugh > fV0Cuts[MaxDCADaughvsDaugh])   continue;
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 7
      // V0: Cosine of the pointing angle
      if(fV0Cuts[MinCosPAK0] > 0.){
        if(sV0->fPA < fV0Cuts[MinCosPAK0]){
          sV0->bIsCandidateK0s = kFALSE;
        }
      }
      if(fV0Cuts[MaxCosPALambda] > 0.){
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
        if(sV0->fDecayRadius<fV0Cuts[MinDecayRadius])   {
            continue;
        }
      }
      FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      iCutIndex++;

      // 9
      // Daughters: pseudorapidity cut
      if(fV0Cuts[DaughMaxEta] > 0.){
        if((TMath::Abs(sPosDaugh->fEta) > fV0Cuts[DaughMaxEta]) || (TMath::Abs(sNegDaugh->fEta) > fV0Cuts[DaughMaxEta]))
          continue;
            FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 10
      // V0: rapidity cut
      if(fV0Cuts[MaxV0Rap] > 0.){
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
      if(fV0Cuts[MaxLifeTime] > 0.){
        if(sV0->fLifetimeK0 > fV0Cuts[MaxLifeTime] * dCTauK0s)
          sV0->bIsCandidateK0s = kFALSE;
      }
      if(fV0Cuts[MaxLifeTime] > 0.){
        if(sV0->fLifetimeLambda > fV0Cuts[MaxLifeTime] * dCTauLambda)
        {
          sV0->bIsCandidateLambda = kFALSE;
          sV0->bIsCandidateALambda = kFALSE;
        }
      }
      if(fV0Cuts[MaxLifeTime] > 0.){
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 12
      // Daughter PID
      if(fV0Cuts[MaxSigmadEdxTPC] > 0.){
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
        if(sV0->fArmenterosPt < TMath::Abs(0.2 * sV0->fArmenterosAlpha)){
          sV0->bIsCandidateK0s = kFALSE;
        }
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      // 14
      // Invariant mass peak selection
      Double_t dMassPeakWindowK0s = fV0Cuts[InvarMassWindowK0]; //0.010; // LF p-p
      Double_t dMassPeakWindowLambda = fV0Cuts[InvarMassWindowLambda]; //0.005; // LF p-p
      Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
      Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

      if(TMath::Abs(sV0->fMassK0 - dMassPDGK0s) < dMassPeakWindowK0s)
        sV0->bIsInPeakK0s = kTRUE;
      if(TMath::Abs(sV0->fMassLambda - dMassPDGLambda) < dMassPeakWindowLambda)
        sV0->bIsInPeakLambda = kTRUE;
      if(TMath::Abs(sV0->fMassAntilambda - dMassPDGLambda) < dMassPeakWindowLambda)
        sV0->bIsInPeakALambda = kTRUE;

      if(fV0Cuts[DoMassWindow]){
        if(sV0->bIsInPeakK0s){
          sV0->bIsCandidateLambda = kFALSE;
          sV0->bIsCandidateALambda = kFALSE;
        }
        if(sV0->bIsInPeakLambda){
          sV0->bIsCandidateK0s = kFALSE;
        }
        if(sV0->bIsInPeakALambda){
          sV0->bIsCandidateK0s = kFALSE;
       }
        FillV0Candidates(sV0->bIsCandidateK0s, sV0->bIsCandidateLambda, sV0->bIsCandidateALambda, iCutIndex);
      }
      iCutIndex++;

      if(sV0->bIsCandidateK0s) {fh2dKshortMassVsPt->Fill(sV0->fPt, sV0->fMassK0,1); }
      if(sV0->bIsCandidateLambda) {fh2dLamdaMassVsPt->Fill(sV0->fPt, sV0->fMassLambda,1); }
      if(sV0->bIsCandidateALambda) {fh2dAnLamdaMassVsPt->Fill(sV0->fPt, sV0->fMassAntilambda,1); }

      //Is true MC V0?
      if(fIsPythia){
        int iVetoDec=0;
        iVetoDec=(int) GetV0MCVeto(fAODIn,v0,sV0->bIsCandidateK0s,sV0->bIsCandidateLambda,sV0->bIsCandidateALambda);
        Double_t fIsMCTrueK0=0.;
        Double_t fIsMCTrueLambda=0.;
        Double_t fIsMCTrueALambda=0.;
        if(iVetoDec==1) fIsMCTrueK0=1.;
        if(iVetoDec==2) fIsMCTrueLambda=1.;
        if(iVetoDec==3) fIsMCTrueALambda=1.;

        //printf("VetoDecision=%i, fIsMCTrueK0=%f, fIsMCTrueLambda=%f, fIsMCTrueALamda=%f\n",iVetoDec,fIsMCTrueK0,fIsMCTrueLambda,fIsMCTrueALambda);

        AliEmcalJet * jetrec  = 0x0;
        double fJetPt=0;
        if(!jetconrec)printf("No jet container with reconstructed jets!\n");

        while ((jetrec = jetconrec->GetNextAcceptJet())){
            fJetPt= jetrec->Pt();
            //if(!(fJetContainerData->GetRhoParameter() == 0x0)){
            //        fJetPt = fJetPt - fJetContainerData->GetRhoVal() * jetrec->Area();
            //}
            if(fJetPt < 5.) continue;
            bool bIsInCone=IsParticleInCone(v0, jetrec, fJetRadius);
            //printf("VetoDecision=%i, fIsMCTrueK0=%f, fIsMCTrueLambda=%f, fIsMCTrueALamda=%f, bIsINCone=%i\n",iVetoDec,fIsMCTrueK0,fIsMCTrueLambda,fIsMCTrueALambda,bIsInCone);

            // make inclusive signed imp. parameter constituent histograms
            if(sV0->bIsCandidateK0s &&bIsInCone) {
              Double_t valueKInJC[5] = {sV0->fMassK0, sV0->fPt, sV0->fEta, fJetPt, fIsMCTrueK0};
              fhnV0InJetK0s->Fill(valueKInJC);
              //printf("massk0=%f, pt=%f, eta=%f,  jetpt=%f, istrue=%f\n",sV0->fMassK0, sV0->fPt, sV0->fEta,fJetPt,fIsMCTrueK0);
            }
            if(sV0->bIsCandidateLambda && bIsInCone) {
              Double_t valueLInJC[5] = {sV0->fMassLambda, sV0->fPt, sV0->fEta, fJetPt,fIsMCTrueLambda};
              fhnV0InJetLambda->Fill(valueLInJC);
              //printf("masslambda=%f, pt=%f, eta=%f,  jetpt=%f, istrue=%f\n",sV0->fMassK0, sV0->fPt, sV0->fEta,fJetPt,fIsMCTrueLambda);
            }
            if(sV0->bIsCandidateALambda && bIsInCone) {
              Double_t valueLInJC[5] = {sV0->fMassAntilambda, sV0->fPt, sV0->fEta, fJetPt,fIsMCTrueALambda};
              fhnV0InJetALambda->Fill(valueLInJC);
              //printf("massalambda=%f, pt=%f, eta=%f, jetpt=%f, istrue=%f\n",sV0->fMassK0, sV0->fPt, sV0->fEta,fJetPt,fIsMCTrueALambda);
            }

        }
        jetconrec->ResetCurrentID();
        jetrec=NULL;
        //delete jetrec;

      }
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
    //printf("################## Select candidates End########################\n");

}

Bool_t AliAnalysisTaskHFJetIPQA::Run(){

    //*******************************
    //Selection
    FillGeneralHistograms();
    /*Vertex Pos Selection*/

    AliAODEvent *ev = dynamic_cast<AliAODEvent*>(InputEvent());
    fEventVertex = dynamic_cast<AliAODVertex*>(ev->GetPrimaryVertex());
    fIsSameEvent_n1 = kFALSE;
    fIsSameEvent_n2 = kFALSE;
    fIsSameEvent_n3 = kFALSE;

    IsEventAccepted(ev);


    Bool_t HasImpactParameter = kFALSE;
    Double_t dca[2] = {-99999,-99999};
    Double_t cov[3] = {-99999,-99999,-99999};
    Double_t TrackWeight       = 1;
    AliVTrack* trackV = NULL;
    fIsEsd =  (InputEvent()->IsA()==AliESDEvent::Class())? kTRUE : kFALSE;
   // EventwiseCleanup();
    if(fIsPythia){
        if(fIsEsd){
            fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent()) ;
            if (!fMCEvent){
                AliError("Could not retrieve  MC particles! Returning");
                return kFALSE;
            }
        }
        else{
            fMCArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
            if (!fMCArray){
                AliError("Could not retrieve AOD MC particles! Returning");
                return kFALSE;
            }
        }
      jetcongen = static_cast<AliJetContainer*>(fJetCollArray.At(1));
      jetcongen->ResetCurrentID();
    }

    jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    jetconrec->ResetCurrentID();


    if(fApplyV0Rej){
      SelectV0Candidates(ev);
      if(fIsPythia) GetV0MCTrueCandidates(ev);
    }
    IncHist("fh1dEventsAcceptedInRun",1);

    //if(fIsPythia)FillParticleCompositionEvent(); //Added  for in cone vs inclusive comparison

    FillHist("fh1dNoParticlesPerEvent",InputEvent()->GetNumberOfTracks(),1);

    //************************
    //Investigate Track Corrections: Smearing & particle weighting
    for(long itrack= 0; itrack<InputEvent()->GetNumberOfTracks();++itrack)
    {
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
        if(fRunSmearing)SmearTrack((AliAODTrack*)trackV);
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
        if(!HasImpactParameter)  continue;
        /*if(fIsPythia){
            Int_t corrpartidx =-1;
            double ppt;
            //if(fDoMCCorrection) TrackWeight *= GetMonteCarloCorrectionFactor(trackV,corrpartidx,ppt);
        }*/
       FillTrackHistograms(trackV,dca,cov,TrackWeight);
    }
    //**********************************
    //JetMatching between generated and reconstructed Jets

    //printf("In Program %f < jetpt <%f, %f < jeteta < %f\n",fAnalysisCuts[bAnalysisCut_MinJetPt],fAnalysisCuts[bAnalysisCut_MaxJetPt],fAnalysisCuts[bAnalysisCut_MinJetEta],fAnalysisCuts[bAnalysisCut_MaxJetEta] );
    FillHist("fh1dNoJetsPerEvent",jetconrec->GetNJets(),1);
    if (!jetconrec) return kFALSE;
    AliEmcalJet * jetgen  = nullptr;

    if(fIsPythia){
        if(!MatchJetsGeometricDefault()) AliInfo("Jet matching did not succeed!");
        jetcongen->ResetCurrentID();
        while ((jetgen = jetcongen->GetNextAcceptJet()))
        {
            if (!jetgen) continue;
            //Int_t jetflavour =0;
            //Bool_t is_udgjet = kFALSE;
            //jetflavour =IsMCJetPartonFast(jetgen,fJetRadius,is_udgjet);
            //FillGenHistograms(jetflavour, jetgen);
        }
        jetcongen->ResetCurrentID();
        jetconrec->ResetCurrentID();
    }

    //*************************************************
    // Loop over reconstructed/matched jets for template creation and analysis
    AliEmcalJet * jetrec  = nullptr;
    AliEmcalJet * jetmatched  = nullptr;
    jetconrec->ResetCurrentID();
    Double_t jetpt=0;
    while ((jetrec = jetconrec->GetNextAcceptJet()))
    {//start jetloop
        if(!jetrec) continue;
        jetpt = jetrec->Pt();
        if(fDoUnderlyingEventSub)jetpt=DoUESubtraction(jetcongen, jetconrec,jetrec, jetpt);

        FillHist("fh1dJetArea",jetrec->Area(),1);

        //________________________
        //Determination of Jet Flavour
        Int_t jetflavour=0;
        Bool_t is_udgjet = kFALSE;
        if(fIsPythia){
          jetmatched = nullptr;
          jetmatched =jetrec->MatchedJet();
          if(jetmatched){
            jetflavour = IsMCJetPartonFast(jetmatched,fJetRadius,is_udgjet); //Event based association to save memory
          }
          else{
            jetflavour=Unid;
          }
        }
        FillRecHistograms( jetflavour, jetpt, jetrec->Eta(),jetrec->Phi());
        if(fDoLundPlane)RecursiveParents(jetrec, jetconrec);

        //_____________________________
        //Determination of impact parameters
        std::vector<SJetIpPati> sImpParXY,sImpParXYZ,sImpParXYSig,sImpParXYZSig;
        AliVParticle *vp=0x0;
        int NJetParticles=0;  //Used for counting particles per jet
        int isV0=kFALSE;
        Double_t dca[2] = {-99999,-99999};
        Double_t cov[3] = {-99999,-99999,-99999};
        Double_t sign=0;

        for(UInt_t i = 0; i < jetrec->GetNumberOfTracks(); i++) {//start trackloop
          TrackWeight=1;
          isV0=kFALSE;
          Double_t xyzatcda[3];

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
          if (fIsPythia&&!IsTrackAccepted((AliAODTrack*)trackV,jetflavour))   continue;

          isV0=IsV0Daughter(trackV);
          ++NJetParticles;

          //FillTrackTypeResHists();

          if(GetImpactParameterWrtToJet((AliAODTrack*)trackV,(AliAODEvent*)InputEvent(),jetrec,dca,cov,xyzatcda,sign, jetflavour)){
            if(fEventVertex) {
              delete fEventVertex;
              fEventVertex =nullptr;
            }
            Int_t corridx=-1;double ppt;
            //(fIsPythia&&fDoMCCorrection) ? TrackWeight = GetMonteCarloCorrectionFactor(trackV,corridx,ppt) : TrackWeight =1;
            dca[0]=fabs(dca[0]);
            Double_t cursImParXY     =TMath::Abs(GetValImpactParameter(   kXY,dca,cov))*sign;
            Double_t cursImParXYSig  =TMath::Abs(GetValImpactParameter(kXYSig,dca,cov))*sign;
            Double_t cursImParXYZ    =TMath::Abs(GetValImpactParameter(   kXYZ,dca,cov))*sign;
            Double_t cursImParXYZSig =TMath::Abs(GetValImpactParameter(kXYZSig,dca,cov))*sign;
            FillHist("fh2dJetSignedImpParXY"            ,jetpt,cursImParXY,TrackWeight);     //*this->fXsectionWeightingFactor );
            FillHist("fh2dJetSignedImpParXYSignificance",jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );

            const char * subtype [5] = {"Unidentified","udsg","c","b","s"};
            if(fIsPythia){
              FillHist(Form("fh2dJetSignedImpParXY%s",subtype[jetflavour]),jetpt,cursImParXY,TrackWeight);     //*this->fXsectionWeightingFactor );
              FillHist(Form("fh2dJetSignedImpParXYSignificance%s",subtype[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
            }
            double fTrackPt=trackV->Pt();
            double fIPValue=fV0Cuts[fAV0Cut]*TMath::Exp(fV0Cuts[fBV0Cut]*fTrackPt)+fV0Cuts[fCV0Cut];
            //printf("trackpt=%f, IPValue=%f, TrueIP=%f, a=%f, b=%f, c=%f\n", fTrackPt, fIPValue,cursImParXYSig, fV0Cuts[fAV0Cut], fV0Cuts[fBV0Cut], fV0Cuts[fCV0Cut]);
            if(cursImParXYSig>fIPValue){
              //printf("Going into switch!\n");
              switch (isV0){
                case V0No:
                  isV0=V0Rec;
                  //printf("V0No set to V0Rec!\n");
                  break;
                case V0MC:
                  isV0=V0TrueRec;
                  //printf("V0MC set to V0TrueRec!\n");
                  break;
              }
            }

            SJetIpPati a(cursImParXY, TrackWeight,isV0,kFALSE,corridx,trackV->Pt()); sImpParXY.push_back(a);
            SJetIpPati b(cursImParXYZ, TrackWeight,isV0,kFALSE,corridx,trackV->Pt()); sImpParXYZ.push_back(b);
            SJetIpPati c(cursImParXYSig, TrackWeight,isV0,kFALSE,corridx,trackV->Pt());sImpParXYSig.push_back(c);
            SJetIpPati d(cursImParXYZSig, TrackWeight,isV0,kFALSE,corridx,trackV->Pt());sImpParXYZSig.push_back(d);
            //printf("curImParXY=%f, isV0=%i, pt=%f\n",sImpParXYSig.back().first,sImpParXYSig.back().is_V0, sImpParXYSig.back().trackpt);

           }
         }//end trackloop

                FillHist("fh1dParticlesPerJet",NJetParticles,1);
                //_________________________
                //Sorting of Impact Parameters
                bool hasIPs[4] ={kFALSE,kFALSE,kFALSE, kFALSE};
                Double_t ipval[4] = {-9999.,-9999.,-9999., -9999.};
                Double_t ipvalsig[4] = {-9999.,-9999.,-9999., -9999.};

                std::sort(sImpParXY.begin(),sImpParXY.end(),        AliAnalysisTaskHFJetIPQA::mysort);
                std::sort(sImpParXYSig.begin(),sImpParXYSig.end(),  AliAnalysisTaskHFJetIPQA::mysort);
                std::sort(sImpParXYZ.begin(),sImpParXYZ.end(),      AliAnalysisTaskHFJetIPQA::mysort);
                std::sort(sImpParXYZSig.begin(),sImpParXYZSig.end(),AliAnalysisTaskHFJetIPQA::mysort);

                if((int)sImpParXYSig.size()>0) hasIPs[0]=kTRUE;
                if((int)sImpParXYSig.size()>1) hasIPs[1]=kTRUE;
                if((int)sImpParXYSig.size()>2) hasIPs[2]=kTRUE;
                if((int)sImpParXYSig.size()>3) hasIPs[3]=kTRUE;

                if(hasIPs[0]){
                  ipvalsig[0] =sImpParXYSig.at(0).first;
                  ipval[0] =sImpParXY.at(0).first;
                  //printf("HasIP0, ipval[0]=%f, ipvalsig[0]=%f\n", ipval[0], ipvalsig[0]);
                }
                if(hasIPs[1]){
                  ipvalsig[1] =sImpParXYSig.at(1).first;
                  ipval[1] =sImpParXY.at(1).first;
                  //printf("HasIP1, ipval[1]=%f, ipvalsig[1]=%f\n",ipval[1], ipvalsig[1]);
                }
                if(hasIPs[2]){
                  ipvalsig[2] =sImpParXYSig.at(2).first;
                  ipval[2] =sImpParXY.at(2).first;
                  //printf("HasIP2, ipval[2]=%f, ipvalsig[2]=%f\n", ipval[2], ipvalsig[2]);
                }
                if(hasIPs[3]){
                  ipvalsig[3] =sImpParXYSig.at(3).first;
                  ipval[3] =sImpParXY.at(3).first;
                  //printf("HasIP3, ipval[3]=%f, ipvalsig[3]=%f\n", ipval[3], ipvalsig[3]);
                }

                //if(hasIPs[0])printf("N=1: cursImParXY=%f, TrackWeight=%f,corridx=%i, pt=%f\n",sImpParXYSig.at(0).first, sImpParXYSig.at(0).second, sImpParXYSig.at(0).trackLabel, sImpParXYSig.at(0).trackpt);
                //if(hasIPs[1])printf("N=2: cursImParXY=%f, TrackWeight=%f, corridx=%i, pt=%f\n",sImpParXYSig.at(1).first, sImpParXYSig.at(1).second, sImpParXYSig.at(1).trackLabel, sImpParXYSig.at(1).trackpt);
                //if(hasIPs[2])printf("N=3: cursImParXY=%f, TrackWeight=%f, corridx=%i, pt=%f\n",sImpParXYSig.at(2).first, sImpParXYSig.at(2).second, sImpParXYSig.at(2).trackLabel, sImpParXYSig.at(2).trackpt);
                //printf("*********************************************************\n");

                //_____________________________
                //V0 tag decisions
                //printf("isV0=%i, jetflaovur=%i, jetpt=%f", );
                bool isV0Jet=kFALSE;
                if((hasIPs[0])&&(!fIsPythia)&&(sImpParXYSig[0].is_V0==V0Rec))isV0Jet=kTRUE;
                //printf("New jetflavour=%i, isV0Jet=%i\n",jetflavour, isV0Jet);

                if(fIsPythia){
                  if(hasIPs[0])FillV0EfficiencyHists(sImpParXYSig[0].is_V0, jetflavour, jetpt, isV0Jet);
                  for(long unsigned iTrack=0;iTrack<sImpParXYSig.size();iTrack++){
                    FillTrackIPvsPt(sImpParXYSig[iTrack].is_V0,sImpParXYSig[iTrack].trackpt,sImpParXYSig[iTrack].first,jetflavour);
                  }

                  //_______________________________
                  //IP Template Generation
                  FillIPTypePtHists(jetflavour, jetpt, hasIPs);
                }
                for(int iN=0;iN<3;iN++){
                  if(!hasIPs[iN]) continue;
                  //printf("iN=%i, jetflavour=%i xy=%f, xysig=%f\n",iN,jetflavour,sImpParXY.at(iN).first,sImpParXYSig.at(iN).first);
                  Double_t params [4] ={sImpParXY.at(iN).first,sImpParXYSig.at(iN).first,sImpParXYZ.at(iN).first,sImpParXYZSig.at(iN).first};
                  FillIPTemplateHists(jetpt,iN,jetflavour, params);
                }

                //____________________________________________
                //TAGGING
                if(fDoTCTagging!=TCNo&&fDoProbTagging!=ProbNo) AliFatal("Don't do track counting and probability tagging simultaneously!");

                bool ** kTagDec=new bool*[fNThresholds];
                for(int iThresh=0;iThresh<fNThresholds;iThresh++){
                  kTagDec[iThresh]=new bool[6];
                  for(int iType=0;iType<6;iType++){
                    kTagDec[iThresh][iType]=0;
                  }
                }
                //**************
                //MC Track Counting
                if(fDoTCTagging!=TCNo){
                  if(fIsPythia||((!fIsPythia)&&(!isV0Jet))){
                      //printf("isV0Jet=%i\n", isV0Jet);
                      if(fUseSignificance){DoTCTagging(jetpt, hasIPs,ipvalsig, kTagDec);}
                      else{DoTCTagging(jetpt, hasIPs,ipval, kTagDec);}
                  }
                  if(fDoMCEffs){
                    //printf("Filling Efficiency hists\n");
                    FillEfficiencyHists(kTagDec, jetflavour, jetpt,hasIPs[0]);
                  }
                  FillTaggedJetPtDistribution(kTagDec,jetpt);
                }
                //**************
                //Probability Dists
                double probval=0;
                probval=GetTrackProbability(jetpt,hasIPs, ipvalsig);
                //Generation of Track Probability Hists
                if(fDoJetProb&&(probval>0)){
                  //printf("Doing Jet Probability!\n");
                  FillProbabilityHists(jetpt,  probval, jetflavour, kTagDec);
                  FillProbThreshHists(probval, ipvalsig, jetpt, jetflavour, hasIPs, kTagDec);
                }
                //Probability Tagging
                /*if(fDoProbTagging!=ProbNo){
                  DoProbTagging(probval, jetpt,kTagDec);
                  FillEfficiencyHists(kTagDec, jetflavour, jetpt,hasIPs[0]);
                }*/

                if(sImpParXY.size()!=0){
                  FillHist("fh2dNoAcceptedTracksvsJetArea",(int)sImpParXY.size(),jetrec->Area(),1);
                }
                sImpParXY.clear();
                sImpParXYSig.clear();
                sImpParXYZ.clear();
                sImpParXYZSig.clear();

                for(int iThresh=0;iThresh<fNThresholds;iThresh++){
                  delete kTagDec[iThresh];
                }
                delete [] kTagDec;
              }//end jetloop*/
            return kTRUE;
          }



                        Double_t AliAnalysisTaskHFJetIPQA::GetLocalAlphaAOD(AliAODTrack * track)
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
                        }

                        Double_t AliAnalysisTaskHFJetIPQA::GetTrackCurvature(AliAODTrack * track)
                        {
    // convert to AliExternalTrackParam
                            AliVEvent * eev = (AliVEvent*)InputEvent();

                            AliExternalTrackParam etp; etp.CopyFromVTrack(track);

                            return etp.GetC(eev->GetMagneticField());
                        }


/*! \brief IsSelected
 *
 *
 * Enable event selection
 */




/*! \brief SetUseMonteCarloWeighingLinus
 *
 *
 * Setter for MC composition correction factors
 */
void AliAnalysisTaskHFJetIPQA::SetUseMonteCarloWeighingLinus(TH1F *Pi0, TH1F *Eta, TH1F *EtaP, TH1F *Rho, TH1F *Phi, TH1F *Omega, TH1F *K0s, TH1F *Lambda, TH1F *ChargedPi, TH1F *ChargedKaon, TH1F *Proton, TH1F *D0, TH1F *DPlus, TH1F *DStarPlus, TH1F *DSPlus, TH1F *LambdaC, TH1F *BPlus, TH1F *B0, TH1F *LambdaB, TH1F *BStarPlus)
  {
  for(Int_t i =1 ; i< Pi0->GetNbinsX()+1;++i){
  fBackgroundFactorLinus[bIdxPi0][i-1] =Pi0->GetBinContent(i);
  fBackgroundFactorLinus[bIdxEta][i-1] =Eta->GetBinContent(i);
  fBackgroundFactorLinus[bIdxEtaPrime][i-1] =EtaP->GetBinContent(i);
  fBackgroundFactorLinus[bIdxRho][i-1] =Rho->GetBinContent(i);
  fBackgroundFactorLinus[bIdxPhi][i-1] =Phi->GetBinContent(i);
  fBackgroundFactorLinus[bIdxOmega][i-1] =Omega->GetBinContent(i);
  fBackgroundFactorLinus[bIdxK0s][i-1] =K0s->GetBinContent(i);
  fBackgroundFactorLinus[bIdxLambda][i-1] =Lambda->GetBinContent(i);
  fBackgroundFactorLinus[bIdxPi][i-1] =ChargedPi->GetBinContent(i);
  fBackgroundFactorLinus[bIdxKaon][i-1] =ChargedKaon->GetBinContent(i);
  fBackgroundFactorLinus[bIdxProton][i-1] =Proton->GetBinContent(i);
  fBackgroundFactorLinus[bIdxD0][i-1] =D0->GetBinContent(i);
  fBackgroundFactorLinus[bIdxDPlus][i-1] =DPlus->GetBinContent(i);
  fBackgroundFactorLinus[bIdxDStarPlus][i-1] =DStarPlus->GetBinContent(i);
  fBackgroundFactorLinus[bIdxDSPlus][i-1] =DSPlus->GetBinContent(i);
  fBackgroundFactorLinus[bIdxLambdaC][i-1] =LambdaC->GetBinContent(i);
  fBackgroundFactorLinus[bIdxBPlus][i-1] =BPlus->GetBinContent(i);
  fBackgroundFactorLinus[bIdxB0][i-1] =B0->GetBinContent(i);
  fBackgroundFactorLinus[bIdxLambdaB][i-1] =LambdaB->GetBinContent(i);
  fBackgroundFactorLinus[bIdxBStarPlus][i-1] =BStarPlus->GetBinContent(i);
  }
  return;
}

void AliAnalysisTaskHFJetIPQA::SetFlukaFactor(TGraph* GraphOmega, TGraph* GraphXi, TGraph* K0Star, TGraph* Phi)
  {
   fGraphOmega=(TGraph*)GraphOmega;
   fGraphXi=(TGraph*)GraphXi;
   fK0Star=(TGraph*)K0Star;
   fPhi=(TGraph*)Phi;

   return;
  }



/*! \brief SetResFunction
 *
 *
 * Setter for resolution function (currently unused)
 */

                        Bool_t AliAnalysisTaskHFJetIPQA::SetResFunctionPID(const char * filename){
                            TFile  * jetProbfile=TFile::Open(filename,"READ");

                            int bins_low[5]  = {0,1,2,4,6};
                            int bins_high[5]  = {1,2,4,6,255};
                            int its_hits[4]  = {6,5,4,3};

                            const char * type[5] ={"Electron","Pion","Kaon","Proton",""};
                            for (int k=0;k<5;++k){
                                for (int i=0;i<4;++i){
                                    for (int j=0;j<5;++j){
                                        TGraph *fResulFkt = 0x0;
                                        if(k==4) fResulFkt = (TGraph*)jetProbfile->Get(Form("fResulFkt_ITS_%i_PT_%i_to_%i",its_hits[i],bins_low[j],bins_high[j]));
                                        else fResulFkt = (TGraph*)jetProbfile->Get(Form("fResulFkt_ITS_%i_PT_%i_to_%i_%s",its_hits[i],bins_low[j],bins_high[j],type[k]));
                                        if(!fResulFkt){
                                            return kFALSE;
                                        }
                                        else {
                                            fResolutionFunction [20*k + 4*j +i] = *fResulFkt;
                                            Printf("Added %i %i %i -> [%i][%i][%i]->[%i]",j,k,i,j,i,k,20*k + 4*j +i  );
                                            delete fResulFkt;
                                        }
                                    }
                                }
                            }
                            if(jetProbfile)jetProbfile->Close();

                            fUsePIDJetProb =kTRUE;
                            return kTRUE;

                        }

                        void AliAnalysisTaskHFJetIPQA::FillCorrelations(bool bn[3],double v[3], double jetpt ){
    //Fill all possible same jet distributions
                             fTREE_n1 = -999;
                             fTREE_n2 = -999;
                             fTREE_n3 = -999;
                             fTREE_pt = -999;
                                     
                             if(fUseTreeForCorrelations){
                                if(bn[0] && bn[1]){
                                   fTREE_n1 = v[0];
                                   fTREE_n2 = v[1];
                                   fTREE_n3 = v[2];
                                   fTREE_pt = jetpt;
                                   fCorrelationCrossCheck->Fill();                                   
                                }
                             }   
                            else if (fFillCorrelations && !fUseTreeForCorrelations )
                            {
                            if (bn[0] && bn[1]) {
                                FillHist("fh2dInclusiveCorrelationN1N2",v[0],v[1]);
                                if(jetpt>10 && jetpt <20)    FillHist("fh2dGreater10_20GeVCorrelationN1N2",v[0],v[1]);
                                if(jetpt>20 && jetpt <30)    FillHist("fh2dGreater20_30GeVCorrelationN1N2",v[0],v[1]);
                                if(jetpt>30 && jetpt <100)    FillHist("fh2dGreater30_100GeVCorrelationN1N2",v[0],v[1]);

                            }
                            if (bn[0] && bn[2]) {
                                FillHist("fh2dInclusiveCorrelationN1N3",v[1],v[2]);
                                if(jetpt>10 && jetpt <20)    FillHist("fh2dGreater10_20GeVCorrelationN1N3",v[0],v[2]);
                                if(jetpt>20 && jetpt <30)    FillHist("fh2dGreater20_30GeVCorrelationN1N3",v[0],v[2]);
                                if(jetpt>30 && jetpt <100)    FillHist("fh2dGreater30_100GeVCorrelationN1N3",v[0],v[2]);
                                FillHist("fh2dInclusiveCorrelationN2N3",v[0],v[2]);
                                if(jetpt>10 && jetpt <20)    FillHist("fh2dGreater10_20GeVCorrelationN2N3",v[1],v[2]);
                                if(jetpt>20 && jetpt <30)    FillHist("fh2dGreater20_30GeVCorrelationN2N3",v[1],v[2]);
                                if(jetpt>30 && jetpt <100)    FillHist("fh2dGreater30_100GeVCorrelationN2N3",v[1],v[2]);
                            }
    //Fill if possible mix distributions
                            bool n3wasReady = false;
                            double storedn3=-999;
                            if(bn[0]){

                                if(fIsMixSignalReady_n2) {
                                    double n2 =0;
                                    if(GetMixDCA(2,n2)) {
                                        FillHist("fh2dInclusiveCorrelationN1N2mix",v[0],n2);
                                        if(jetpt>10 && jetpt <20)    FillHist("fh2dGreater10_20GeVCorrelationN1N2mix",v[0],n2);
                                        if(jetpt>20 && jetpt <30)    FillHist("fh2dGreater20_30GeVCorrelationN1N2mix",v[0],n2);
                                        if(jetpt>30 && jetpt <100)    FillHist("fh2dGreater30_100GeVCorrelationN1N2mix",v[0],n2);
                                    }
                                }
                                if(fIsMixSignalReady_n3) {
                                    double n3 =0;
                                    if(GetMixDCA(3,n3)) {
                                        n3wasReady=true;
                                        storedn3=n3;
                                        FillHist("fh2dInclusiveCorrelationN1N3mix",v[0],n3);
                                        if(jetpt>10 && jetpt <20)    FillHist("fh2dGreater10_20GeVCorrelationN1N3mix",v[0],n3);
                                        if(jetpt>20 && jetpt <30)    FillHist("fh2dGreater20_30GeVCorrelationN1N3mix",v[0],n3);
                                        if(jetpt>30 && jetpt <100)    FillHist("fh2dGreater30_100GeVCorrelationN1N3mix",v[0],n3);
                                    }
                                }
                            }

                            if(bn[1]) {
                                if(n3wasReady) {
                                    double n3 =0;
                                    n3=storedn3;
                                    FillHist("fh2dInclusiveCorrelationN2N3mix",v[1],n3);
                                    if(jetpt>10 && jetpt <20)    FillHist("fh2dGreater10_20GeVCorrelationN2N3mix",v[1],n3);
                                    if(jetpt>20 && jetpt <30)    FillHist("fh2dGreater20_30GeVCorrelationN2N3mix",v[1],n3);
                                    if(jetpt>30 && jetpt <100)    FillHist("fh2dGreater30_100GeVCorrelationN2N3mix",v[1],n3);
                                }
                            }
                            }
                            return;
                        }



void AliAnalysisTaskHFJetIPQA::UserCreateOutputObjects(){
  Printf("Analysing Jets with Radius: R=%f\n",fJetRadius);

  TString BJetCuts[25] = {
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
    "Kink",//22
    "TPC Refit",//23
    "ITS Refit" //24
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
  }/*
    //Graphs currently not in use
                                               const Int_t gfProtonN = 9;
                                                const Int_t gfAntiProtonN = 18;
                                                const Int_t gfAntiLambdaN = 34;
                                                const Int_t gfLambdaN = 2;
                                                const Int_t gfKMinusN =13 ;
                                                Double_t gfProtonX [gfProtonN] ={0,0.534483,1.29741,2.21552,3.0819,3.92241,4.5819,5.39655,1000};
                                                Double_t gfProtonY [gfProtonN] ={0.990964,0.990964,0.990964,0.990964,0.990964,0.990964,0.990964,0.990964,0.990964};
                                                Double_t gfAntiProtonX [gfAntiProtonN]  = {0,0.806034,0.922414,1.09052,1.28448,1.5431,1.73707,1.89224,2.17672,2.43534,2.74569,3.06897,
                                                   3.52155,3.88362,4.38793,5.03448,5.38362, 1000};
                                                   Double_t gfAntiProtonY [gfAntiProtonN]  = {0.922892,0.922892,	0.930723,	0.939157,0.94397,0.95241,0.956627,0.959639,0.964458,
                                                       0.966867,0.971084,0.974096,0.978313,0.98012,0.983735,0.986747,0.989157,0.989157};
                                                       Double_t gfAntiLambdaX [gfAntiLambdaN]  = {0.,0.55555,0.64646,0.75757,	0.84848,0.94949,1.06061,1.15152,1.24242,1.35354,1.44444,
                                                           1.54545,1.66667,1.75758,1.84848,1.9596,2.09091,2.30303,2.50505,2.68687,2.90909,3.11111,
                                                           3.31313,3.51515,3.69697,3.89899,4.20202,4.66667,5.21212,5.74747,6.50505,7.51515,9.0101,1000};
                                                           Double_t gfAntiLambdaY [gfAntiLambdaN]  = {0.864925,0.864925,0.895896,0.908209,0.915672,0.921269,0.926866,0.931343,0.935821,0.938806,0.942164,
                                                               0.945149,0.947761,0.95,0.952612,0.954478,0.957836,0.960821,0.96306,0.965672,0.968657,0.970149,
                                                               0.972015,0.973507,0.975,0.976493,0.978358,0.981343,0.983955,0.986194,0.988433,0.991045,0.991045,0.991045};
                                                               Double_t gfLambdaX [gfLambdaN]          =	{0.,1000};
                                                               Double_t gfLambdaY [gfLambdaN]          = {0.991045,0.991045};
                                                               Double_t gfKMinusX [gfKMinusN]          =	{0,0.54741,0.74137,1.03879,1.36207,1.96983,2.52586,3.0819,3.67672,4.19397,5.03448,5.44828,1000};
                                                               Double_t gfKMinusY [gfKMinusN]          = {0,0.979518,0.983133,0.987349,0.989759,0.992169,0.993976,0.996386,0.995783,0.998193,0.99759,1,1000};
                                                               fGeant3FlukaProton 	   = new TGraph(gfProtonN,gfProtonX,gfProtonY);
                                                               fGeant3FlukaAntiProton = new TGraph(gfAntiProtonN,gfAntiProtonX,gfAntiProtonY);
                                                               fGeant3FlukaLambda     = new TGraph(gfLambdaN,gfLambdaX,gfLambdaY);
                                                               fGeant3FlukaAntiLambda = new TGraph(gfAntiLambdaN,gfAntiLambdaX,gfAntiLambdaY);
                                                               fGeant3FlukaKMinus 	   = new TGraph(gfKMinusN,gfKMinusX,gfKMinusY);



  //General Information
  /*fh1dEventRejectionRDHFCuts = (TH1D*)AddHistogramm("fh1dEventRejectionRDHFCuts","fh1dEventRejectionRDHFCuts;reason;count",12,0,12);
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(1,"Event accepted");
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(2,"Event rejected");
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(3,"Wrong physics selection");
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(4,"No vertex");
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(5,"No contributors");
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(6,"Less than 10 contributors");
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(7,">10cm vertex Z distance");
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(8,"Bad diamond X distance");
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(9,"Bad diamond Y distance");
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(10,"Bad diamond Z distance");
  fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(11,"Chi2 vtx >1.5 ");*/
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
  fHistManager.CreateTH1("fh1dJetRecPtAccepted","accepted detector level jets;pt (GeV/c); count",500,0,250,"s");
  fHistManager.CreateTH1("fh1dJetArea","fh1dJetArea;# Jet Area",100,0,1,"s");
  fHistManager.CreateTH1("fh1dParticlesPerJet","fh1dParticlesPerJet;#, Particles/Jet",100,0,100,"s");
  fHistManager.CreateTH2("fh2dNoAcceptedTracksvsJetArea","fh2dNoAcceptedTracksvsJetArea;No Accepted Tracks;JetArea",20,0,20,100,0,1);
  //MC properties
  if(fIsPythia){
    /*fHistManager.CreateTH1("fh1dJetGenPt","generator level jets;pt (GeV/c); count",250,0,250,"s");
    fHistManager.CreateTH1("fh1dJetGenPtUnidentified","generator level jets (no flavour assigned);pt (GeV/c); count",250,0,250,"s");
    fHistManager.CreateTH1("fh1dJetGenPtudsg","generator level udsg jets;pt (GeV/c); count",250,0,250,"s");
    fHistManager.CreateTH1("fh1dJetGenPtc","generator level c jets;pt (GeV/c); count",250,0,250,"s");
    fHistManager.CreateTH1("fh1dJetGenPtb","generator level b jets;pt (GeV/c); count",250,0,250,"s");
    fHistManager.CreateTH1("fh1dJetGenPts","generator level s jets;pt (GeV/c); count",250,0,250,"s");
    fHistManager.CreateTH2("fh2dJetGenPtVsJetRecPt","detector momentum response;gen pt;rec pt",500,0,250,500,0,250,"s");*/
    fHistManager.CreateTH1("fh1dJetRecPtudsg","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPtUnidentified","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPtc","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPtb","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPts","detector level jets;pt (GeV/c); count",500,0,250,"s");
  }

  fh1DCutInclusive=(TH1D*)AddHistogramm("fh1DCutInclusive","fh1DCutInclusive",30,0,30);
  fh1dCutudg=(TH1D*)AddHistogramm("fh1dCutudg","fh1dCutudg",30,0,30);
  fh1dCutc=(TH1D*)AddHistogramm("fh1dCutc","fh1dCutc",30,0,30);
  fh1dCutb=(TH1D*)AddHistogramm("fh1dCutb","fh1dCutb",30,0,30);
  fh1dCuts=(TH1D*)AddHistogramm("fh1dCuts","fh1dCuts",30,0,30);
  fh1dCuts->GetXaxis()->LabelsOption("v");

  for(Int_t iBin = 0; iBin < 25; iBin++){
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
  const char * tagtype[6] = {"Full","Single1st","Single2nd","Single3rd","Double","Triple"};

  for(int iThresh=0;iThresh<fNThresholds;iThresh++){
    for(int iType=0;iType<6;iType++){
      if(fDoJetProb){
          fHistManager.CreateTH2(Form("h2DProbDistsTag_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",200, 0, 1,500, 0, 250);
          fHistManager.CreateTH2(Form("h2DLNProbDistsTag_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",400, 0, 20,500, 0, 250);
      }
      if(fDoTCTagging!=TCNo)fHistManager.CreateTH1(Form("h1DTagged_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",500, 0, 250);
    }
  }
  if(fDoJetProb){
    fHistManager.CreateTH2("h2DProbDists","h2DProbDistsAll",200, 0, 1,500, 0, 250);
    fHistManager.CreateTH2("h2DLNProbDists","h2DProbDistsAll",400, 0, 20,500, 0, 250);

    fHistManager.CreateTH2("h2DProb1Above0","h2DProb1Above0", 400,0,20,500,0,250);
    fHistManager.CreateTH2("h2DProb2Above0","h2DProb2Above0", 400,0,20,500,0,250);
    fHistManager.CreateTH2("h2DProb3Above0","h2DProb3Above0", 400,0,20,500,0,250);

    if(fIsPythia){
      fHistManager.CreateTH2("h2DProb1Above0_UDSG","h2DProb1Above0_UDSG", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb2Above0_UDSG","h2DProb2Above0_UDSG", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb3Above0_UDSG","h2DProb3Above0_UDSG", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb1Above0_C","h2DProb1Above0_C", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb2Above0_C","h2DProb2Above0_C", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb3Above0_C","h2DProb3Above0_C", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb1Above0_B","h2DProb1Above0_B", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb2Above0_B","h2DProb2Above0_B", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb3Above0_B","h2DProb3Above0_B", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb1Above0_V0","h2DProb1Above0_V0", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb2Above0_V0","h2DProb2Above0_V0", 400,0,20,500,0,250);
      fHistManager.CreateTH2("h2DProb3Above0_V0","h2DProb3Above0_V0", 400,0,20,500,0,250);
    }

    fHistManager.CreateTH2("h2DProb1AboveThresh","h2DProb1AboveThresh", 400,0,20,500,0,250);
    fHistManager.CreateTH2("h2DProb1AbThresh1Ab0","h2DProb1AbThresh1Ab0", 400,0,20,500,0,250);
    fHistManager.CreateTH2("h2DProb1AbThresh2Ab0","h2DProb1AbThresh2Ab0", 400,0,20,500,0,250);

    fHistManager.CreateTH2("h2DProb2AboveThresh","h2DProb2AboveThresh", 400,0,20,500,0,250);
    fHistManager.CreateTH2("h2DProb2AbThresh1Ab0","h2DProb2AbThresh1Ab0", 400,0,20,500,0,250);
    fHistManager.CreateTH2("h2DProb2AbThresh2Ab0","h2DProb2AbThresh2Ab0", 400,0,20,500,0,250);
  }
  if (fIsPythia){
    if(fDoMCEffs){
      for(int iThresh=0;iThresh<fNThresholds;iThresh++){
        for(int iType=0;iType<6;iType++){
          fHistManager.CreateTH1(Form("h1DTrueBTagged_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";jet pt; #",500,0,250);
          fHistManager.CreateTH1(Form("h1DFalseCTagged_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";jet pt; #",500,0,250);
          fHistManager.CreateTH1(Form("h1DFalseUDSGTagged_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";jet pt; #",500,0,250);
          fHistManager.CreateTH1(Form("h1DFalseV0Tagged_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";jet pt; #",500,0,250);
        }
      }
    }
    if(fDoJetProb){
      fHistManager.CreateTH2("h2DLNProbDists_B",";; #",400, 0, 20,500, 0, 250);
      fHistManager.CreateTH2("h2DLNProbDists_C",";; #",400, 0, 20,500, 0, 250);
      fHistManager.CreateTH2("h2DLNProbDists_UDSG",";; #",400, 0, 20,500, 0, 250);
      fHistManager.CreateTH2("h2DLNProbDists_V0",";; #",400, 0, 20,500, 0, 250);

      fHistManager.CreateTH2("h2DProbDists_B",";; #",200, 0, 1,500, 0, 250);
      fHistManager.CreateTH2("h2DProbDists_C",";; #",200, 0, 1,500, 0, 250);
      fHistManager.CreateTH2("h2DProbDists_UDSG",";; #",200, 0, 1,500, 0, 250);
      fHistManager.CreateTH2("h2DProbDists_V0",";; #",200, 0, 1,500, 0, 250);

      for(int iThresh=0;iThresh<fNThresholds;iThresh++){
        for(int iType=0;iType<6;iType++){
            fHistManager.CreateTH2(Form("h2DLNProbDistsTag_B_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",400, 0, 20,500, 0, 250);
            fHistManager.CreateTH2(Form("h2DLNProbDistsTag_C_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",400, 0, 20,500, 0, 250);
            fHistManager.CreateTH2(Form("h2DLNProbDistsTag_UDSG_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",400, 0, 20,500, 0, 250);
            fHistManager.CreateTH2(Form("h2DLNProbDistsTag_V0_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",400, 0, 20,500, 0, 250);

            fHistManager.CreateTH2(Form("h2DProbDistsTag_B_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",200, 0, 1,500, 0, 250);
            fHistManager.CreateTH2(Form("h2DProbDistsTag_C_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",200, 0, 1,500, 0, 250);
            fHistManager.CreateTH2(Form("h2DProbDistsTag_UDSG_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",200, 0, 1,500, 0, 250);
            fHistManager.CreateTH2(Form("h2DProbDistsTag_V0_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";; #",200, 0, 1,500, 0, 250);
        }
      }
    }
  }

    //****************************************
    //Track Impact Parameter Distributions
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

    fHistManager.CreateTH2("fh1dTracksIPvsPt_B","Track IP vs Track Pt; p_{T,Track} (GeV/c); d_{0}",500,0,250,300,0,30);
    fHistManager.CreateTH2("fh1dTracksIPvsPt_V0JetTracks","Track IP vs Track Pt; p_{T,Track} (GeV/c); d_{0}",500,0,250,300,0,30);
    fHistManager.CreateTH2("fh1dTracksIPvsPt_V0inV0Jet","Track IP vs Track Pt; p_{T,Track} (GeV/c); d_{0}",500,0,250,300,0,30);
    fHistManager.CreateTH2("fh1dTracksIPvsPt_V0inBJet","Track IP vs Track Pt; p_{T,Track} (GeV/c); d_{0}",500,0,250,300,0,30);



    //****************************************
    //Pt Distributions for N1,N2,N3 Tracks
    if(fIsPythia){
      for(int iN=1;iN<=3;iN++){
        fHistManager.CreateTH1(Form("fh1dJetRecPt_n_%i_all_Accepted",iN),"detector level jets;pt (GeV/c); count",500,0,250,"s");
        fHistManager.CreateTH1(Form("fh1dTrackPt_n_%i_all_Accepted",iN),"detector level jets;pt (GeV/c); count",500,0,200,"s");
        for(long unsigned iFlav=0;iFlav<sTemplateFlavour.size();iFlav++){
          fHistManager.CreateTH1(Form("fh1dJetRecPt_n_%i_%s_Accepted",iN,sTemplateFlavour[iFlav].Data()),"detector level jets;pt (GeV/c); count",500,0,250,"s");
        }
      }
    }

    //V0Cuts
    if(fApplyV0Rej){
      //V0s from reconstruction
      const Int_t iNDimInJC = 5;
      Int_t binsKInJC[iNDimInJC] = {200, 200, 200, 200,5};
      Double_t xminKInJC[iNDimInJC] = {0.35, 0., -1., 0.,0};
      Double_t xmaxKInJC[iNDimInJC] = {0.65, 50., 1., 200.,5};
      Int_t binsLInJC[iNDimInJC] = {200, 200, 200, 200,5};
      Double_t xminLInJC[iNDimInJC] = {1.05, 0., -1., 0.,0};
      Double_t xmaxLInJC[iNDimInJC] = {1.25, 50., 1., 200.,5};

      fh2dKshortMassVsPt=(TH2D*)AddHistogramm("fh2dKshortMassVsPt","KShort Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,0.35, 0.65);
      fh2dLamdaMassVsPt =(TH2D*)AddHistogramm("fh2dLamdaMassVsPt","Lamda Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,1.05,1.25);
      fh2dAnLamdaMassVsPt =(TH2D*)AddHistogramm("fh2dAnLamdaMassVsPt","Anti Lamda Mass Vs Pt;p_{T} (GeV/c);Mass (GeV)",200,0,50,200,1.05,1.25);

      fhnV0InJetK0s = new THnSparseD("fhnV0InJetK0s", "K0s: Mass vs Pt in jets;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c}); IsMCTrueK0", iNDimInJC, binsKInJC, xminKInJC, xmaxKInJC);
      fhnV0InJetLambda = new THnSparseD("fhnV0InJetLambda", "Lambda: Mass vs Pt in jets;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c}); IsMCTrueLamda", iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);
      fhnV0InJetALambda = new THnSparseD("fhnV0InJetALambda", "ALambda: Mass vs Pt in jets;#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{V0} (GeV/#it{c});#it{#eta}_{V0};#it{p}_{T}^{jet} (GeV/#it{c}); IsMCTrueALamda", iNDimInJC, binsLInJC, xminLInJC, xmaxLInJC);

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
    }

    //Final Tagged JetPt Spectra
    for(int iThresh=0;iThresh<fNThresholds;iThresh++){
      for(int iType=0;iType<6;iType++){
        fHistManager.CreateTH1(Form("h1DTaggedJetPt_%s_%0.2f",tagtype[iType],fFracs[iThresh]),";jet pt; #",500, 0, 250);
      }
    }

    //Impact Parameter Template Generation
    const char * base = "fh2dJetSignedImpPar";
    const char * dim[2]  = {"XY","XYZ"};
    const char * typ[2]  = {"","Significance"};
    const char * ordpar [4] = {"","First","Second","Third"};
    const char * special [1] = {"",/*"McCorr"*/};

    Int_t ptbins = 250;
    Double_t ptlow = 0;
    Double_t pthigh = 250;
    Int_t ipbins = 1000;
    Double_t iplow = -.5;
    Double_t iphigh = .5;
    for (Int_t id = 0;id<1;++id)  // XY or XY/
      for (unsigned int ifl = 0;ifl<sTemplateFlavour.size();++ifl)  //flavour
        for (Int_t io = 0;io<4;++io)        //order parameter
          for (Int_t is = 0;is<1;++is)          //special comment
            for (Int_t it = 1;it<2;++it){           //significance or not
              if(it==1) {
                iplow=-30;
                iphigh=30; //from 30
                if(io==0 && ifl==4) ipbins = 1000;//2000;
                  else  ipbins =1000;//2000;\
              }else {
                iplow=-0.5;
                iphigh=0.5;
                ipbins =1000;//;2000;
              }
              if(id==0)  ipbins =1000;//2000;
                if((fIsPythia||(!fIsPythia && ifl==6))){
                  fHistManager.CreateTH2(Form("%s%s%s%s%s%s",base,dim[id],typ[it],sTemplateFlavour[ifl].Data(),ordpar[io],special[is]),
                                Form("%s%s%s%s%s%s;;",base,dim[id],typ[it],sTemplateFlavour[ifl].Data(),ordpar[io],special[is]),
                                ptbins,ptlow,pthigh,ipbins,iplow,iphigh,"s");
                  //printf("Generating%s%s%s%s%s%s",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]);

                  }
              }

    TIter next(fHistManager.GetListOfHistograms());
    TObject* obj = 0;
    while ((obj = next())) {
      printf("Adding Object %s\n",obj->GetName());
      fOutput->Add(obj);
    }

    PostData(1, fOutput);
}

void AliAnalysisTaskHFJetIPQA::UserExecOnce(){
    AliJetContainer *  jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    fAnalysisCuts[bAnalysisCut_MinJetPt]=jetconrec->GetJetPtCut();
    fAnalysisCuts[bAnalysisCut_MaxJetPt]=jetconrec->GetJetPtCutMax();
    fAnalysisCuts[bAnalysisCut_MinJetEta]=jetconrec->GetMinEta();
    fAnalysisCuts[bAnalysisCut_MaxJetEta]=jetconrec->GetMaxEta();

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
    Int_t version=2;


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
    jetcuts+=fNoJetConstituents;
    jetcuts+="+";
    jetcuts+=fDaughtersRadius;
    jetcuts+="+";
    jetcuts+=Form("%0.1f",fJetRadius);

    printf("Cut Track Settings: %s\n",jetcuts.Data());

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
    trackcuts+=fAnalysisCuts[bAnalysisCut_HasSDD];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_KinkCand];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_HasTPCrefit];
    trackcuts+="+";
    trackcuts+=fAnalysisCuts[bAnalysisCut_HasITSrefit];


    printf("Cut Vertex Settings %s\n", trackcuts.Data());

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
    vertexcuts+=fDoTCTagging;
    vertexcuts+="+";
    vertexcuts+=fDoProbTagging;
    vertexcuts+="+";
    vertexcuts+=fUseSignificance;
    vertexcuts+="+";
    vertexcuts+=Form("%0.3f",fTCThresholdPtFixed);

    fh1dCutsPrinted->SetTitle(jetcuts.Data());
    fh1dCutsPrinted->GetXaxis()->SetTitle(trackcuts.Data());
    fh1dCutsPrinted->GetYaxis()->SetTitle(vertexcuts.Data());

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
            Int_t nTracks=aod->GetNumberOfTracks();
            AliAODTrack * t = nullptr;
            AliExternalTrackParam etp_at_r39_old; etp_at_r39_old.CopyFromVTrack(track);
            etp_at_r39_old.PropagateTo(3.9,InputEvent()->GetMagneticField());
            double angle0 = TMath::ATan2(etp_at_r39_old.Yv(),etp_at_r39_old.Xv());
            double zz0    = etp_at_r39_old.GetZ();
            int nTrksToSkip=1;

            for(Int_t i=0; i<nTracks; i++){
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
void AliAnalysisTaskHFJetIPQA::FillParticleCompositionSpectra(AliEmcalJet * jet,const char * histname ){
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
}

/*! \brief FillParticleCompositionEvent
 *
 *
 * Fill Histogram with Correction Factors vs. pT
 */
void AliAnalysisTaskHFJetIPQA::FillParticleCompositionEvent( ){
    AliVTrack* tr=0x0; 
    for(Int_t j = 0; j < InputEvent()->GetNumberOfTracks(); ++j) {
      tr = (AliVTrack*)(InputEvent()->GetTrack(j));
      if(!tr) continue;
      double pT=0x0;
      Int_t pCorr_indx=-1;
      GetMonteCarloCorrectionFactor(tr,pCorr_indx,pT);
      if(pCorr_indx<0) continue;
      FillHist("fh2dParticleSpectra_Event",pCorr_indx+0.5,pT, 1);     //*this->fXsectionWeightingFactor );
  }
  return;
}


void AliAnalysisTaskHFJetIPQA::GetOutOfJetParticleComposition(AliEmcalJet * jet, int flavour){
    if(!jet) return;
    AliEmcalJet * perp_jet1 =   GetPerpendicularPseudoJet(jet,false);
    AliEmcalJet * perp_jet2 =   GetPerpendicularPseudoJet(jet,true);
    FillParticleCompositionSpectra(jet,"fh2dParticleSpectra_InCone");
    if(flavour==3)  FillParticleCompositionSpectra(jet,"fh2dParticleSpectra_InCone_bjet");
    else if(flavour==2)  FillParticleCompositionSpectra(jet,"fh2dParticleSpectra_InCone_cjet");
    else if(flavour==1)  FillParticleCompositionSpectra(jet,"fh2dParticleSpectra_InCone_lfjet");
    FillParticleCompositionSpectra(perp_jet1,"fh2dParticleSpectra_OutOfCone");
    FillParticleCompositionSpectra(perp_jet2,"fh2dParticleSpectra_OutOfCone");
    if(perp_jet1) delete perp_jet1;
    if(perp_jet2) delete perp_jet2;
}

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

Bool_t AliAnalysisTaskHFJetIPQA::IsTrackAccepted(AliVTrack* track , int jetflavour){
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
        if(fAnalysisCuts[bAnalysisCut_HasSDD]){
          if(!(((AliAODTrack*)track->HasPointOnITSLayer(0))&&(AliAODTrack*)track->HasPointOnITSLayer(1))){
            FillCandidateJet(bAnalysisCut_HasSDD,jetflavour);
            return kFALSE;
          }
        }

        //n=hits in ITS layers
        Int_t SPDSSDHits = (int) track->HasPointOnITSLayer(0) + (int) track->HasPointOnITSLayer(1) + (int) track->HasPointOnITSLayer(4) + (int) track->HasPointOnITSLayer(5);
        if(SPDSSDHits<abs(fAnalysisCuts[bAnalysisCut_MinITSLayersHit])){
            //printf("Throw away due flav=%i ssd points 0 = %i, 1=%i, 2=%i, 3=%i, SPDSSDHits=%i cutvalue=%f\n",jetflavour,track->HasPointOnITSLayer(0),track->HasPointOnITSLayer(1),track->HasPointOnITSLayer(4),track->HasPointOnITSLayer(5),SPDSSDHits, abs(fAnalysisCuts[bAnalysisCut_MinITSLayersHit]));
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
/*! \brief Composition correction factor  getter
 *
 * finds the corresponding re-weighing factor for a certain track:
 *      - find mother particle
 *      - define fluca factor according to mother particle
 */
    Double_t AliAnalysisTaskHFJetIPQA::GetWeightFactor( AliVTrack * track,Int_t &pCorr_indx, double &ppt){
        AliMCParticle *pMCESD = nullptr;
        AliAODMCParticle *pMCAOD = nullptr;
        if(track->GetLabel()< 0) return 1;
        if(!fIsPythia) return 1;
        if(fIsEsd){
            pMCESD = ((AliMCParticle*)MCEvent()->GetTrack(abs(track->GetLabel())));
            if(!(pMCESD)) return 1;
        }
        else {
            pMCAOD = static_cast<AliAODMCParticle*>(fMCArray->At(abs(track->GetLabel())));
            if(!(pMCAOD)) return 1;
        }

        AliVParticle * mcpart = fIsEsd ? (   AliVParticle * )pMCESD:(   AliVParticle * )pMCAOD;
        Bool_t _particlesourcefound(kFALSE);
        Int_t  _particlesourcepdg(mcpart->PdgCode());
        Int_t  _particlesourceidx(25);
        Double_t _particlesourcept(0);

        AliVParticle * mcpartclone = mcpart;
    while(mcpart){//Strangenss
        if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;
        _particlesourcept = mcpart->Pt();
        _particlesourcepdg = abs(mcpart->PdgCode());
        if (IsSelectionParticleStrange(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
            _particlesourcefound = kTRUE;
            break;
        }
        mcpart = GetVParticleMother(mcpart);
    }
    if (!_particlesourcefound) { //heavy mesons to improve templates
        mcpart = mcpartclone;

        while(mcpart){
            if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;
            _particlesourcept = mcpart->Pt();
            _particlesourcepdg = abs(mcpart->PdgCode());
            if (IsSelectionParticleMeson(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
                _particlesourcefound = kTRUE;
                break;
            }
            mcpart = GetVParticleMother(mcpart);
        }
    }
    if (!_particlesourcefound) { //charged hadrons
        mcpart = mcpartclone;
        while(mcpart){
            if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;
            _particlesourcept = mcpart->Pt();
            _particlesourcepdg = abs(mcpart->PdgCode());
            if (IsSelectionParticleOmegaXiSigmaP(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
                _particlesourcefound = kTRUE;
                break;
            }
            mcpart = GetVParticleMother(mcpart);
        }
    }
    if (!_particlesourcefound) {
        mcpart = mcpartclone;
        while(mcpart){
            if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;
            _particlesourcept = mcpart->Pt();
            _particlesourcepdg = abs(mcpart->PdgCode());
            if (IsSelectionParticle(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
                _particlesourcefound = kTRUE;
                break;
            }
            mcpart = GetVParticleMother(mcpart);
        }
    }
    if (!_particlesourcefound) {
        mcpart = mcpartclone;
        while(mcpart){
            if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;
            _particlesourcept = mcpart->Pt();
            if (IsSelectionParticleALICE(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
                _particlesourcefound = kTRUE;
                break;
            }
            mcpart = GetVParticleMother(mcpart);
        }
    }

    if (!_particlesourcefound) return 1.;
    // Do the weighting
    Double_t factor = 1;
    if (_particlesourceidx <0) return 1;
    if (_particlesourceidx >19) return 1;
    Double_t wpos = ((_particlesourcept - 0.15)/ 0.05);
    Double_t  fractpart, intpart;
    fractpart = modf (wpos , &intpart);
    if (fractpart > 0) intpart = intpart + 1;
    Int_t  bin = floor(intpart);
    if (bin > 497) bin = 497;// above weight definition
    if (_particlesourcept < 0.1+ 1E-5) bin = 0; //below weight definition
    if((_particlesourceidx == bIdxSigmaMinus) || (_particlesourceidx == bIdxSigmaPlus))    factor = fBackgroundFactorLinus[bIdxLambda][bin];
    else factor = fBackgroundFactorLinus[_particlesourceidx][bin];
    pCorr_indx = _particlesourceidx;
    Double_t flucafactor = 1;

    switch(mcpart->PdgCode())
    {
        case -bPhi:
        factor=1;
        flucafactor =1./fPhi->Eval(_particlesourcept);
        pCorr_indx=bIdxPhi;
        break;
        case -bK0S892:
        factor=1;
        flucafactor =1./fK0Star->Eval(_particlesourcept);
        pCorr_indx=bIdxK0S892;
        break;
        case -bK0S892plus:
        factor=1;
        flucafactor =1./fK0Star->Eval(_particlesourcept);
        pCorr_indx=bK0S892plus;
        break;
        case -bOmegaBaryon:
        pCorr_indx=bIdxOmegaBaryon;
        factor=1;
        flucafactor =1./fGraphOmega->Eval(_particlesourcept);
        break;
        case -bXiBaryon:
        pCorr_indx=bIdxXiBaryon;
        factor=1;
        flucafactor =1./fGraphXi->Eval(_particlesourcept);
        break;
        case bPhi:
        factor=1;
        flucafactor =1./fPhi->Eval(_particlesourcept);
        pCorr_indx=bIdxPhi;
        break;
        case bK0S892:
        factor=1;
        pCorr_indx=bIdxK0S892;
        flucafactor =1./fK0Star->Eval(_particlesourcept);
        break;
        case bK0S892plus:
        pCorr_indx=bK0S892plus;
        factor=1;
        flucafactor =1./fK0Star->Eval(_particlesourcept);
        break;
        case bOmegaBaryon:
        pCorr_indx=bIdxOmegaBaryon;
        factor=1;
        flucafactor =fGraphOmega->Eval(_particlesourcept);
        break;
        case bXiBaryon:
        pCorr_indx=bIdxXiBaryon;
        factor=1;
        flucafactor =fGraphXi->Eval(_particlesourcept);
        break;

        default:
        break;
    }
    factor*=flucafactor;
    if (factor <= 0 || factor > 100.)  {
        return 1;
    }
    ppt = _particlesourcept;
    return factor ;
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
/*! \brief GetMonteCarloCorrectionFactor
 *
 * Composition  correction base function caller
 */
    Double_t AliAnalysisTaskHFJetIPQA::GetMonteCarloCorrectionFactor(AliVTrack* track,Int_t &pCorr_indx, double &ppt){
        printf("Doing MC Correction.\n");
        double val=  GetWeightFactor(track,pCorr_indx,ppt);
        if(val > 0 ) return val;
        return 1.;
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
    Double_t AliAnalysisTaskHFJetIPQA::GetPtCorrectedMC(const AliEmcalJet *jet)
    {
        AliJetContainer * jetconrec = nullptr;
        jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(1));
        if(jet && jetconrec&&fDoUnderlyingEventSub){
            printf("Correct for Underlying Event.\n");
            return jet->Pt() - jetconrec->GetRhoVal() * jet->Area();
        }
        return -1.;
    }


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

            for(Int_t iPrim = 0 ; iPrim<fMCArray->GetEntriesFast();iPrim++){//start trackloop MC

                        AliAODMCParticle * part = static_cast<AliAODMCParticle*>(fMCArray->At(iPrim));
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
void AliAnalysisTaskHFJetIPQA::RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont){

      std::vector<fastjet::PseudoJet>  fInputVectors;
      fInputVectors.clear();
      fastjet::PseudoJet  PseudoTracks;

      AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();

        if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
          AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
          if (!fTrk) continue;
          //if(fDoTwoTrack==kTRUE && CheckClosePartner(i,fJet,fTrk,fTrackCont)) continue;
          PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
          PseudoTracks.set_user_index(fJet->TrackAt(i)+100);
          fInputVectors.push_back(PseudoTracks);

        }
        fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);



      fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 );

      try {
        fastjet::ClusterSequence fClustSeqSA(fInputVectors, fJetDef);
        std::vector<fastjet::PseudoJet>   fOutputJets;
        fOutputJets.clear();
        fOutputJets=fClustSeqSA.inclusive_jets(0);

       fastjet::PseudoJet jj;
       fastjet::PseudoJet j1;
       fastjet::PseudoJet j2;
       jj=fOutputJets[0];
       double ktaverage=0;
       double thetaverage=0;
       double nall=0;
       double flagSubjet=0;
        while(jj.has_parents(j1,j2)){
          nall=nall+1;
        if(j1.perp() < j2.perp()) swap(j1,j2);
        flagSubjet=0;
        double delta_R=j1.delta_R(j2);
        double z=j2.perp()/(j1.perp()+j2.perp());
        double y =log(1.0/delta_R);
        double lnpt_rel=log(j2.perp()*delta_R);
        double yh=j1.e()+j2.e();
         vector < fastjet::PseudoJet > constitj1 = sorted_by_pt(j1.constituents());
         if(constitj1[0].perp()>fAnalysisCuts[bAnalysisCut_MinTrackPt]) flagSubjet=1;
        if(z>fHardCutOff){
          ktaverage=ktaverage+lnpt_rel;
          thetaverage=thetaverage+delta_R;
        Double_t LundEntries[6] = {y,lnpt_rel,fOutputJets[0].perp(),nall,yh,flagSubjet};
        fHLundIterative->Fill(LundEntries);}
        jj=j1;}
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
    if(h1)  h1->Fill(x,w);
}
/*! \brief FillHist
 *
 * 2d
 */
void AliAnalysisTaskHFJetIPQA::FillHist(const char *name, Double_t x, Double_t y, Double_t w){
    TH2D * h2 =GetHist2D(name);
    if(h2) h2->Fill(x,y,w);
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


Bool_t AliAnalysisTaskHFJetIPQA::GetImpactParameterWrtToJet(const AliAODTrack *track, const AliAODEvent *event, const AliEmcalJet *jet, Double_t *dca, Double_t *cov, Double_t *XYZatDCA, Double_t &jetsign, int jetflavour)
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
    TVector3 JetDir =jetP3.Unit();
    TVector3 D0(XYZatDCA);
    TVector3 vertex(VxVyVz);
    TVector3 DD0(D0.x()-vertex.x(),D0.y()-vertex.y(),0.);   //track impact parameter
    double ps =DD0.Dot(JetDir);
    double value = DD0.Mag()*(ps/fabs(ps));                 //absolut track impact parameter
    jetsign  = TMath::Sign(1.,value);                       //sign
    TVector3 dd0 =DD0.Unit();

    //track properties
    AliExternalTrackParam etp_track;    etp_track.CopyFromVTrack(track);
    Double_t xa,xb,xyz_jet_global[3],xyz_track_global[3];

    etp_jet.GetDCA(&etp_track, event->GetMagneticField(), xa, xb);
    etp_jet.GetXYZAt(xa, event->GetMagneticField(),xyz_jet_global);
    etp_track.GetXYZAt(xb, event->GetMagneticField(),xyz_track_global);
    etp_track.PropagateTo(xb,event->GetMagneticField());

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
  if(fDoTCTagging==TCIPFixedPt){
      fNThresholds=1;
  }
  for(int iProbSet=0;iProbSet<fNThresholds;iProbSet++){
    TObjArray* oa=(TObjArray*)threshs[iProbSet];

    //printf("Pointer oa=%p\n",oa);

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

// Read Threshold Histograms
//==============================================================================
void AliAnalysisTaskHFJetIPQA::ReadThresholdHists(TString PathToThresholds, TString taskname, int nTCThresh){
    TFile* fileThresholds=TFile::Open(PathToThresholds.Data());
    if(!fileThresholds ||(fileThresholds&& !fileThresholds->IsOpen())){AliError(Form("%s :: File with threshold values not found",taskname.Data()));}

    Printf("%s :: File %s successfully loaded, setting up threshold functions.",taskname.Data(),PathToThresholds.Data());

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

  for(int iN=0;iN<3;iN++){
    h2DProbLookup.push_back((TH2D*)oLookup->At(iN));
  }
}

void AliAnalysisTaskHFJetIPQA::DoTCTagging(double jetpt, bool* hasIPs, double* ipval, bool **kTagDec){    
  //printf("Start TCTagging!\n");
  //threshold values for tracks with largest, second and third largest IP
  int iJetPtBin=h1DThresholdsFirst[0]->FindBin(jetpt);
  double IPthresN1[fNThresholds];  //IP threshold values for individual separation power
  double IPthresN2[fNThresholds];
  double IPthresN3[fNThresholds];

  if(fDoTCTagging==TCIPSig){
    for(int iN=0;iN<fNThresholds;iN++){
      IPthresN1[iN]=h1DThresholdsFirst[iN]->GetBinContent(iJetPtBin);
      IPthresN2[iN]=h1DThresholdsSecond[iN]->GetBinContent(iJetPtBin);
      IPthresN3[iN]=h1DThresholdsThird[iN]->GetBinContent(iJetPtBin);
    }
  }
  if(fDoTCTagging==TCIPFixedPt){
    for(int iN=0;iN<fNThresholds;iN++){
      IPthresN1[iN]=fTCThresholdPtFixed;
      IPthresN2[iN]=fTCThresholdPtFixed;
      IPthresN3[iN]=fTCThresholdPtFixed;
    }
  }

  for(int iThresh=0;iThresh<fNThresholds;iThresh++){
    if(!hasIPs[0]) continue;
    //printf("DoTCTagging:\n");
    //printf("      iJetPtBin=%i, IPthresN1=%f, IPthresN2=%f, IPthresN3=%f\n", iJetPtBin, IPthresN1[iThresh],IPthresN2[iThresh], IPthresN3[iThresh]);


    if(hasIPs[2]){
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
        //printf("Single3rd %f!\n",fFracs[iThresh]);
          kTagDec[iThresh][Full]=kTRUE; kTagDec[iThresh][Single3rd]=kTRUE;
        }
      }
    }

    if(hasIPs[1]){
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
    if(hasIPs[0]){
      //printf("ipval[0]=%f", ipval[0]);
      if(kTagLevel<2){
        //printf("Single catch\n");
        if(ipval[0]>IPthresN1[iThresh]) {
          //printf("Single1st %f!\n",fFracs[iThresh]);
          kTagDec[iThresh][Full]=kTRUE; kTagDec[iThresh][Single1st]=kTRUE;
        }
      }
    }
    /*printf("Testing kTagLevel::\n ");
    for(int iThresh=0;iThresh<fNThresholds;iThresh++){
      for(int iType=0;iType<6;iType++){
        printf("iThresh=%f, %i, kTagDec=%i\n",fFracs[iThresh],iType,kTagDec[iThresh][iType]);
      }
    }*/
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

void AliAnalysisTaskHFJetIPQA::FillEfficiencyHists(bool** kTagDec, int jetflavour, double jetpt, bool hasIPs){
  //printf("Receiving BTagged decision: %i\n", kTagDec[Full]);
  for(int iThresh=0;iThresh<fNThresholds;iThresh++){
    //printf("kTagDec=%i, jetflavour=%i, hasIPs=%i\n",kTagDec[iThresh][Full],jetflavour, hasIPs);
    if(!fIsPythia&&kTagDec[iThresh][Full]){
      FillHist(Form("h1DTagged_Full_%0.2f",fFracs[iThresh]),jetpt,1);
      if(kTagDec[iThresh][Single1st]) FillHist(Form("h1DTagged_Single1st_%0.2f",fFracs[iThresh]), jetpt, 1);
      if(kTagDec[iThresh][Single2nd]) FillHist(Form("h1DTagged_Single2nd_%0.2f",fFracs[iThresh]), jetpt, 1);
      if(kTagDec[iThresh][Single3rd])FillHist(Form("h1DTagged_Single3rd_%0.2f",fFracs[iThresh]), jetpt, 1);
      if(kTagDec[iThresh][Double])FillHist(Form("h1DTagged_Double_%0.2f",fFracs[iThresh]), jetpt, 1);
      if(kTagDec[iThresh][Triple])FillHist(Form("h1DTagged_Triple_%0.2f",fFracs[iThresh]), jetpt, 1);
    }

    if(fIsPythia){
      if(kTagDec[iThresh][Full]&&(jetflavour!=Unid)){
          FillHist(Form("h1DTagged_Full_%0.2f",fFracs[iThresh]),jetpt,1);
          if(kTagDec[iThresh][Single1st]) FillHist(Form("h1DTagged_Single1st_%0.2f",fFracs[iThresh]), jetpt, 1);
          if(kTagDec[iThresh][Single2nd]) FillHist(Form("h1DTagged_Single2nd_%0.2f",fFracs[iThresh]), jetpt, 1);
          if(kTagDec[iThresh][Single3rd])FillHist(Form("h1DTagged_Single3rd_%0.2f",fFracs[iThresh]), jetpt, 1);
          if(kTagDec[iThresh][Double])FillHist(Form("h1DTagged_Double_%0.2f",fFracs[iThresh]), jetpt, 1);
          if(kTagDec[iThresh][Triple])FillHist(Form("h1DTagged_Triple_%0.2f",fFracs[iThresh]), jetpt, 1);
      }
      if(kTagDec[iThresh][Full]&&(jetflavour==B)&&hasIPs){
        //printf("################################ Before: FoundJet with tagindex=%i!\n",kTagDec[iThresh][Full]);

        FillHist(Form("h1DTrueBTagged_Full_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single1st]) FillHist(Form("h1DTrueBTagged_Single1st_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single2nd]) FillHist(Form("h1DTrueBTagged_Single2nd_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single3rd])FillHist(Form("h1DTrueBTagged_Single3rd_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Double])FillHist(Form("h1DTrueBTagged_Double_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Triple])FillHist(Form("h1DTrueBTagged_Triple_%0.2f",fFracs[iThresh]), jetpt, 1);

       //printf("################################ FoundJet with tagindex=%i!\n",kTagDec[iThresh][Full]);
      }
      if(kTagDec[iThresh][Full]&&(jetflavour==C)&&hasIPs){
        FillHist(Form("h1DFalseCTagged_Full_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single1st]) FillHist(Form("h1DFalseCTagged_Single1st_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single2nd]) FillHist(Form("h1DFalseCTagged_Single2nd_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single3rd]) FillHist(Form("h1DFalseCTagged_Single3rd_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Double])FillHist(Form("h1DFalseCTagged_Double_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Triple])FillHist(Form("h1DFalseCTagged_Triple_%0.2f",fFracs[iThresh]), jetpt, 1);
        //printf("################################ CMistagged: flavour is=%i with tagindex=%i!\n",jetflavour,kTagDec[iThresh][Full]);
      }
      if(kTagDec[iThresh][Full]&&(jetflavour==UDSG)&&hasIPs){
        FillHist(Form("h1DFalseUDSGTagged_Full_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single1st]) FillHist(Form("h1DFalseUDSGTagged_Single1st_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single2nd]) FillHist(Form("h1DFalseUDSGTagged_Single2nd_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single3rd]) FillHist(Form("h1DFalseUDSGTagged_Single3rd_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Double])FillHist(Form("h1DFalseUDSGTagged_Double_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Triple])FillHist(Form("h1DFalseUDSGTagged_Triple_%0.2f",fFracs[iThresh]), jetpt, 1);
        //printf("################################ LightMistagged: flavour is=%i with tagindex=%i!\n",jetflavour,kTagDec[iThresh][Full]);
      }
      if((kTagDec[iThresh][Full])&&((jetflavour==UDSGV0)||(jetflavour==CV0))&&hasIPs){
        FillHist(Form("h1DFalseV0Tagged_Full_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single1st]) FillHist(Form("h1DFalseV0Tagged_Single1st_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single2nd]) FillHist(Form("h1DFalseV0Tagged_Single2nd_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Single3rd]) FillHist(Form("h1DFalseV0Tagged_Single3rd_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Double])FillHist(Form("h1DFalseV0Tagged_Double_%0.2f",fFracs[iThresh]), jetpt, 1);
        if(kTagDec[iThresh][Triple])FillHist(Form("h1DFalseV0Tagged_Triple_%0.2f",fFracs[iThresh]), jetpt, 1);
        //printf("################################ v0Mistagged: flavour is=%i with tagindex=%i!\n",jetflavour,kTagDec[iThresh][Full]);
      }
      if(!kTagDec[iThresh][Full]&&jetflavour==B&&hasIPs){
        //printf("################################ Missed one: flavour is=%i\n", jetflavour);
      }
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

double AliAnalysisTaskHFJetIPQA::IntegrateIP(int iJetPtBin, int iIPBin, int iN){
  int iZeroIPBin=h2DProbLookup[iN]->GetXaxis()->FindBin(0.);
  int iStartIPBin=h2DProbLookup[iN]->GetXaxis()->FindBin(-25);

  double prob=h2DProbLookup[iN]->Integral(iStartIPBin,iIPBin,iJetPtBin,iJetPtBin);
  prob=prob/(h2DProbLookup[iN]->Integral(iStartIPBin,iZeroIPBin,iJetPtBin,iJetPtBin));

  //printf("Integrate: 0x=%f, lowx=%f, upx=%f, lowy=%f, upy=%f, prob=%f\n", h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iZeroIPBin), h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iStartIPBin),h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iIPBin),h2DProbLookup[iN]->GetYaxis()->GetBinLowEdge(iJetPtBin),h2DProbLookup[iN]->GetYaxis()->GetBinLowEdge(iJetPtBin+1),prob);

  return prob;
}

double AliAnalysisTaskHFJetIPQA::GetTrackProbability(double jetpt, bool* hasIPs, double* ipval){
  //printf("Printing h2DProbLook\n");
  double prob=1;
  double probval[3]={0};
  int iJetPtBin=h2DProbLookup[0]->GetYaxis()->FindBin(jetpt);;
  int iIPBin[3]={0};
  //printf("ipval1=%f, ipval2=%f, ipval3=%f\n",ipval[0],ipval[1],ipval[2]);

  for(int iN=0;iN<3;iN++){
    if(!hasIPs[iN])continue;
    if(ipval[iN]<0) continue;
    iIPBin[iN]=h2DProbLookup[iN]->GetXaxis()->FindBin(-ipval[iN]);
    probval[iN]=IntegrateIP(iJetPtBin,iIPBin[iN], iN);
    //probval[iN]=h2DProbLookup[iN]->GetBinContent(iIPBin[iN],iJetPtBin);
    //printf("iN=%i, iIPBin=%i, ipval=%f, lowerIP=%f, higherIP=%f, || iJetPtBin=%i, jetpt=%f, lowerjetpt=%f, higherjetpt=%f, prob=%f\n", iN, iIPBin[iN],ipval[iN],h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iIPBin[iN]),
    //        h2DProbLookup[iN]->GetXaxis()->GetBinLowEdge(iIPBin[iN]+1), iJetPtBin,jetpt, h2DProbLookup[iN]->GetYaxis()->GetBinLowEdge(iJetPtBin),  h2DProbLookup[iN]->GetYaxis()->GetBinLowEdge(iJetPtBin+1),probval[iN]);
    prob=prob*probval[iN];
  }

  double prob3=prob-prob*TMath::Log(prob)+prob*(TMath::Log(prob)*TMath::Log(prob))*0.5;
  double prob2=prob-prob*TMath::Log(prob);
  double prob1=prob;


  if(hasIPs[2]&&ipval[2]>0){
      //printf("3 tracks with prob1=%f, prob2=%f, prob3=%f, prob=%f, prob3=%f\n",probval[0],probval[1],probval[2],prob,prob3);
      return prob3;
  }
  if(hasIPs[1]&&ipval[1]>0){
      //printf("2 tracks with prob1=%f, prob2=%f prob=%f, prob2=%f\n",probval[0],probval[1],prob,prob2);
      return prob2;
  }
  if(hasIPs[0]&&ipval[0]>0){
      //printf("1 track with prob1=%f, prob=%f, prob1=%f\n",probval[0],prob,prob1);
      return prob1;
  }
  return 0;
}

void AliAnalysisTaskHFJetIPQA::FillProbabilityHists(double jetpt,double probval,int jetflavour,bool **kTagDec){
  //  printf("Filling iflavou=%i, jetpt=%f, probval=%f into histogram\n", jetflavour,jetpt,probval);
  double lnprobval=-TMath::Log(probval);

  //printf("logval=%f\n",lnprobval);
  if(fIsPythia){
    FillHist(Form("h2DProbDists%s",sTemplateFlavour[jetflavour].Data()),probval,jetpt,1);     //*this->fXsectionWeightingFactor );
    FillHist(Form("h2DLNProbDists%s",sTemplateFlavour[jetflavour].Data()),lnprobval,jetpt,1);     //*this->fXsectionWeightingFactor );
  }
  //untagged Probability hist
  FillHist(Form("h2DProbDists"),probval,jetpt,1);     //*this->fXsectionWeightingFactor );
  FillHist(Form("h2DLNProbDists"),lnprobval,jetpt,1);     //*this->fXsectionWeightingFactor );
  const char * tagtype[6] = {"Full","Single1st","Single2nd","Single3rd","Double","Triple"};

  //tagged Probability hists
  if(fDoTCTagging!=TCNo){
    for(int iThresh=0;iThresh<fNThresholds;iThresh++){
      for(int iType=0;iType<6;iType++){
          if(kTagDec[iThresh][iType]){
            FillHist(Form("h2DProbDistsTag_%s_%0.2f",tagtype[iType],fFracs[iThresh]),probval,jetpt,1);
            FillHist(Form("h2DLNProbDistsTag_%s_%0.2f",tagtype[iType],fFracs[iThresh]),lnprobval,jetpt,1);
          }
          switch(jetflavour){
            case UDSG:
              if(iThresh==0&&iType==0){
                FillHist(Form("h2DProbDists_UDSG"),probval,jetpt,1);
                FillHist(Form("h2DLNProbDists_UDSG"),lnprobval,jetpt,1);
                //printf("Filling total UDSG %i\n",jetflavour);
              }
              if(kTagDec[iThresh][iType]){
                FillHist(Form("h2DProbDistsTag_UDSG_%s_%0.2f",tagtype[iType],fFracs[iThresh]),probval,jetpt,1);
                FillHist(Form("h2DLNProbDistsTag_UDSG_%s_%0.2f",tagtype[iType],fFracs[iThresh]),lnprobval,jetpt,1);
                //printf("Filling UDSG %i %s %0.2f\n", jetflavour, tagtype[iType],fFracs[iThresh]);
              }
              break;

            case B:
              if(iThresh==0&&iType==0){
                FillHist(Form("h2DProbDists_B"),probval,jetpt,1);
                FillHist(Form("h2DLNProbDists_B"),lnprobval,jetpt,1);
                //printf("Filling total B %i\n",jetflavour);
              }
              if(kTagDec[iThresh][iType]){
                FillHist(Form("h2DProbDistsTag_B_%s_%0.2f",tagtype[iType],fFracs[iThresh]),probval,jetpt,1);
                FillHist(Form("h2DLNProbDistsTag_B_%s_%0.2f",tagtype[iType],fFracs[iThresh]),lnprobval,jetpt,1);
                //printf("Filling B %i %s %0.2f\n,", jetflavour,tagtype[iType],fFracs[iThresh]);
              }
              break;

            case C:
              if(iThresh==0&&iType==0){
                FillHist(Form("h2DProbDists_C"),probval,jetpt,1);
                FillHist(Form("h2DLNProbDists_C"),lnprobval,jetpt,1);
                //printf("Filling total C %i\n",jetflavour);
              }
              if(kTagDec[iThresh][iType]){
                FillHist(Form("h2DProbDistsTag_C_%s_%0.2f",tagtype[iType],fFracs[iThresh]),probval,jetpt,1);
                FillHist(Form("h2DLNProbDistsTag_C_%s_%0.2f",tagtype[iType],fFracs[iThresh]),lnprobval,jetpt,1);
                //printf("Filling C %i %s %0.2f\n",jetflavour, tagtype[iType],fFracs[iThresh]);
              }
              break;

            case UDSGV0:
              if(iThresh==0&&iType==0){
                FillHist(Form("h2DProbDists_V0"),probval,jetpt,1);
                FillHist(Form("h2DLNProbDists_V0"),lnprobval,jetpt,1);
                //printf("Filling total V0 %i\n",jetflavour);
              }
              if(kTagDec[iThresh][iType]){
                FillHist(Form("h2DProbDistsTag_V0_%s_%0.2f",tagtype[iType],fFracs[iThresh]),probval,jetpt,1);
                FillHist(Form("h2DLNProbDistsTag_V0_%s_%0.2f",tagtype[iType],fFracs[iThresh]),lnprobval,jetpt,1);
                //printf("Filling V0 %i %s %0.2f\n",jetflavour, tagtype[iType],fFracs[iThresh]);
              }
              break;

            case CV0:
              if(iThresh==0&&iType==0){
                FillHist(Form("h2DProbDists_V0"),probval,jetpt,1);
                FillHist(Form("h2DLNProbDists_V0"),lnprobval,jetpt,1);
                //printf("Filling total V0 %i\n",jetflavour);
              }
              if(kTagDec[iThresh][iType]){
                FillHist(Form("h2DProbDistsTag_V0_%s_%0.2f",tagtype[iType],fFracs[iThresh]),probval,jetpt,1);
                FillHist(Form("h2DLNProbDistsTag_V0_%s_%0.2f",tagtype[iType],fFracs[iThresh]),lnprobval,jetpt,1);
                //printf("Filling V0 %i %s %0.2f\n", jetflavour,tagtype[iType],fFracs[iThresh]);
              }
              break;
          }

      }
    }
  }
}

void AliAnalysisTaskHFJetIPQA::FillProbThreshHists(double probval, double* ipval, double jetpt, int jetflavour, bool* hasIPs, bool** kTagDec){
  double lnprobval=-TMath::Log(probval);
  TString sFlavour[6]={"Unid","UDSG","C","B","V0","V0"};

  //printf("ipval0=%f, ipval1=%f, ipval2=%f, single1st=%i, double=%i\n", ipval[0],ipval[1],ipval[2],kTagDec[0][Single1st],kTagDec[0][Double]);

  //Untagged
  if(ipval[0]>0){FillHist("h2DProb1Above0",lnprobval,jetpt,1);}
  if((ipval[0]>0)&&(ipval[1]>0)){FillHist("h2DProb2Above0",lnprobval,jetpt,1);}
  if((ipval[0]>0)&&(ipval[1]>0)&&(ipval[2]>0)){FillHist("h2DProb3Above0",lnprobval,jetpt,1);}
  if(fIsPythia){
    if(ipval[0]>0) {FillHist(Form("h2DProb1Above0_%s",sFlavour[jetflavour].Data()),lnprobval,jetpt,1);}
    if((ipval[0]>0)&&(ipval[1]>0)) {FillHist(Form("h2DProb2Above0_%s",sFlavour[jetflavour].Data()),lnprobval,jetpt,1);}
    if((ipval[0]>0)&&(ipval[1]>0)&&(ipval[2]>0))  {FillHist(Form("h2DProb3Above0_%s",sFlavour[jetflavour].Data()),lnprobval,jetpt,1);}
  }
  //Single1st
  if(kTagDec[0][Single1st]) FillHist("h2DProb1AboveThresh",lnprobval,jetpt,1);
  if(kTagDec[0][Single1st]&&(ipval[1]>0)) FillHist("h2DProb1AbThresh1Ab0",lnprobval,jetpt,1);
  if(kTagDec[0][Single1st]&&(ipval[1]>0)&&(ipval[2]>0)) FillHist("h2DProb1AbThresh2Ab0",lnprobval,jetpt,1);

  //Single2nd
  if(kTagDec[0][Double]) FillHist("h2DProb2AboveThresh",lnprobval,jetpt,1);
  if(kTagDec[0][Double]&&(ipval[2]>0)) FillHist("h2DProb2AbThresh1Ab0",lnprobval,jetpt,1);
  if(kTagDec[0][Double]&&(ipval[2]>0)&&(ipval[3]>0)) FillHist("h2DProb2AbThresh2Ab0",lnprobval,jetpt,1);
}

void AliAnalysisTaskHFJetIPQA::Terminate(Option_t *){

    printf("\n*********************************\n");
    printf("Corrections:\n");
    printf("    MC Corrections (Data/MC+Fluka):%i\n",fDoMCCorrection);
    printf("    Track Smearing:%i\n",fRunSmearing );
    printf("    Underlying Event Subtraction:%i\n", fDoUnderlyingEventSub);
    printf("*********************************\n");
}
