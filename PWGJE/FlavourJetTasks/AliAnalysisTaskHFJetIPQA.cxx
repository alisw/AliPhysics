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
fh1dTracksAccepeted(nullptr),
fh1dCuts(nullptr),
fh2dManifoldParton(nullptr),
fh2dLightNotContrib(nullptr),
fh2dCharmNotContrib(nullptr),
fh2dBottomNotContrib(nullptr),
fh2dLightNMatch(nullptr),
fh2dCharmNMatch(nullptr),
fh2dBottomNMatch(nullptr),
fh2dLightDeta(nullptr),
fh2dCharmDeta(nullptr),
fh2dBottomDeta(nullptr),
fHLundIterative(nullptr),
fHistManager(),
fEventVertex(nullptr),
fPidResponse(nullptr),
fRunSmearing(kFALSE),
fUsePIDJetProb(kFALSE),
fDoMCCorrection(kFALSE),
fDoUnderlyingEventSub(kFALSE),
fDoFlavourMatching(kFALSE),
fPythia6(kFALSE),
fPythia8(kFALSE),
fFillCorrelations(kFALSE),
fParam_Smear_Sigma(1.),
fParam_Smear_Mean(0.),
fGlobalVertex(kFALSE),
fDoNotCheckIsPhysicalPrimary(kFALSE),
fDoJetProb(kFALSE),
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
cCuts(0),
fMCArray(nullptr),
fMCEvent(nullptr),
fESDTrackCut(nullptr),
fVertexer(nullptr),
fMcEvtSampled(kFALSE),
fBackgroundFactorLinus{0},
fPUdsgJet(100),fPSJet(100),fPCJet(100),fPBJet(100),
fJetCont(10),
fHardProcess(0,SQuarks(0,0,0,0)),
fQuarkVec(0,SQuarks(0,0,0,0)),
fAnalysisCuts{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
fCombined(nullptr),
fXsectionWeightingFactor(1),
fProductionNumberPtHard(-1),
fJetRadius(0.4),
fDaughtersRadius(1),
fNoJetConstituents(2),
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
    for(Int_t i =0 ; i<498;++i)for(Int_t j =0 ; j<19;++j)  fBackgroundFactorLinus[j][i]=1.;
        SetNeedEmcalGeom(kFALSE);
    SetOffTrigger(AliVEvent::kINT7);
    SetUseAliAnaUtils(kTRUE,kTRUE);
    SetVzRange(-10,10);
    SetUseSPDTrackletVsClusterBG(kTRUE);
    for(Int_t i =0 ; i<200;++i)this->fResolutionFunction[i].Set(1000);
    DefineOutput(1,  AliEmcalList::Class()) ;
}
AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA(const char *name):
AliAnalysisTaskEmcalJet(name, kTRUE),
fh1dTracksAccepeted(nullptr),
fh1dCuts(nullptr),
fh2dManifoldParton(nullptr),
fh2dLightNotContrib(nullptr),
fh2dCharmNotContrib(nullptr),
fh2dBottomNotContrib(nullptr),
fh2dLightNMatch(nullptr),
fh2dCharmNMatch(nullptr),
fh2dBottomNMatch(nullptr),
fh2dLightDeta(nullptr),
fh2dCharmDeta(nullptr),
fh2dBottomDeta(nullptr),
fHLundIterative(nullptr),
fHistManager(name),
fEventVertex(nullptr),
fPidResponse(nullptr),
fRunSmearing(kFALSE),
fUsePIDJetProb(kFALSE),
fDoMCCorrection(kFALSE),
fDoUnderlyingEventSub(kFALSE),
fDoFlavourMatching(kFALSE),
fPythia6(kFALSE),
fPythia8(kFALSE),
fFillCorrelations(kFALSE),
fParam_Smear_Sigma(1.),
fParam_Smear_Mean(0.),
fGlobalVertex(kFALSE),
fDoNotCheckIsPhysicalPrimary(kFALSE),
fDoJetProb(kFALSE),
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
cCuts(0),
fMCArray(nullptr),
fMCEvent(nullptr),
fESDTrackCut(nullptr),
fVertexer(nullptr),
fMcEvtSampled(kFALSE),
fBackgroundFactorLinus{0},
fPUdsgJet(100),fPSJet(100),fPCJet(100),fPBJet(100),
fJetCont(10),
fHardProcess(0,SQuarks(0,0,0,0)),
fQuarkVec(0,SQuarks(0,0,0,0)),
fAnalysisCuts{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
fCombined(nullptr),
fXsectionWeightingFactor(1.),
fProductionNumberPtHard(-1),
fJetRadius(0.4),
fDaughtersRadius(1),
fNoJetConstituents(2),
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
    SetUseAliAnaUtils(kTRUE,kTRUE);
    SetVzRange(-10,10);
    SetUseSPDTrackletVsClusterBG(kTRUE);
    SetMakeGeneralHistograms(kTRUE);
    SetDefaultAnalysisCuts();
    DefineOutput(1,  AliEmcalList::Class()) ;

    for(Int_t i =0 ; i<498;++i)for(Int_t j =0 ; j<19;++j)  fBackgroundFactorLinus[j][i]=1;
        for(Int_t i =0 ; i<200;++i)this->fResolutionFunction[i].Set(1000);

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
    fAnalysisCuts[bAnalysisCut_DCAJetTrack]     = 0.07;
    fAnalysisCuts[bAnalysisCut_MaxDecayLength]  = 10.;
    fAnalysisCuts[bAnalysisCut_MaxDCA_XY]       = 0.5 ;
    fAnalysisCuts[bAnalysisCut_MaxDCA_Z]        = 0.5;
    fAnalysisCuts[bAnalysisCut_NContibutors]    = 3 ;
    fAnalysisCuts[bAnalysisCut_RelError_Y]      = 0.2;
    fAnalysisCuts[bAnalysisCut_RelError_Z]      = 0.2;
    fAnalysisCuts[bAnalysisCut_Sigma_Y]         = 0.3;
    fAnalysisCuts[bAnalysisCut_Sigma_Z]         = 0.3;
    fAnalysisCuts[bAnalysisCut_SigmaDiamond]    = 2.;
    fAnalysisCuts[bAnalysisCut_MaxVtxZ]         = 10.;
    fAnalysisCuts[bAnalysisCut_Z_Chi2perNDF]    =3.5*3.5;
    fAnalysisCuts[bAnalysisCut_MinTrackPt]      =1.;
    fAnalysisCuts[bAnalysisCut_MinTrackPtMC]    =.5;
    fAnalysisCuts[bAnalysisCut_MinTPCClus]      =100;
    fAnalysisCuts[bAnalysisCut_MinJetPt]        =0;
    fAnalysisCuts[bAnalysisCut_MaxJetPt]        =1000;
    fAnalysisCuts[bAnalysisCut_MinJetEta]       =-0.9;
    fAnalysisCuts[bAnalysisCut_MaxJetEta]       =0.9;
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
    FillHist("fh2dTracksImpParXYZSignificance",GetValImpactParameter(kXYZSig,dca,cov),track->Pt(),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh2dTracksImpParZSignificance",GetValImpactParameter(kZSig,dca,cov),track->Pt(),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh1dTracksImpParXY",GetValImpactParameter(kXY,dca,cov),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh1dTracksImpParXYZ",GetValImpactParameter(kXYZ,dca,cov),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh1dTracksImpParXYSignificance",GetValImpactParameter(kXYSig,dca,cov),1);     //*this->fXsectionWeightingFactor );
    FillHist("fh1dTracksImpParXYZSignificance",GetValImpactParameter(kXYZSig,dca,cov),1);     //*this->fXsectionWeightingFactor );
    if(fIsPythia){
        FillHist("fh1dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        FillHist("fh1dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        FillHist("fh1dTracksImpParXYSignificance_McCorr",GetValImpactParameter(kXYSig,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        FillHist("fh1dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),weight);     //*this->fXsectionWeightingFactor );
        FillHist("fh2dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),track->Pt(),weight);     //*this->fXsectionWeightingFactor );
        FillHist("fh2dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),track->Pt(),weight);     //*this->fXsectionWeightingFactor );
        FillHist("fh2dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),track->Pt(),weight);     //*this->fXsectionWeightingFactor );
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


Bool_t AliAnalysisTaskHFJetIPQA::Run(){
  

    FillGeneralHistograms();
    /*Vertex Pos Selection*/
    fEventVertex = dynamic_cast<const AliAODVertex*>(InputEvent()->GetPrimaryVertex());
    fIsSameEvent_n1 = kFALSE;
    fIsSameEvent_n2 = kFALSE;
    fIsSameEvent_n3 = kFALSE;

    if(!fEventVertex) {
        return kFALSE;
    }
    if(fEventVertex->GetNContributors()<1) {
        return kFALSE;
    }
    if(fEventVertex->GetChi2perNDF()>fAnalysisCuts[bAnalysisCut_Z_Chi2perNDF]) {
        return kFALSE;
    }
    if(fEventVertex->GetNContributors()<(int)(fAnalysisCuts[bAnalysisCut_NContibutors])) {
        return kFALSE;
    }
    if(TMath::Abs(fEventVertex->GetZ())>=fAnalysisCuts[bAnalysisCut_MaxVtxZ]) {
        return kFALSE;
    }
    IncHist("fh1dEventsAcceptedInRun",1);

    Bool_t HasImpactParameter = kFALSE;
    Double_t dca[2] = {-99999,-99999};
    Double_t cov[3] = {-99999,-99999,-99999};
    Double_t TrackWeight       = 1;
    AliVTrack* trackV = NULL;
    fIsEsd =  (InputEvent()->IsA()==AliESDEvent::Class())? kTRUE : kFALSE;
   // EventwiseCleanup();
    if(fIsPythia)
    {
        if(fIsEsd)
        {
            fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent()) ;
            if (!fMCEvent)
            {
                AliError("Could not retrieve  MC particles! Returning");
                return kFALSE;
            }
        }
        else{
            fMCArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
            if (!fMCArray)
            {
                AliError("Could not retrieve AOD MC particles! Returning");
                return kFALSE;
            }
        }
    }
    //if(fIsPythia)FillParticleCompositionEvent(); //Added  for in cone vs inclusive comparison

    FillHist("fh1dNoParticlesPerEvent",InputEvent()->GetNumberOfTracks(),1);

    for(long itrack= 0; itrack<InputEvent()->GetNumberOfTracks();++itrack)
    {
        trackV = static_cast<AliVTrack*>(InputEvent()->GetTrack(itrack));
        if(!trackV) {
            AliInfo("Could not retrieve Track");
            continue;
        }
        IncHist("fh1dTracksAccepeted",1);
        if(!IsTrackAccepted(trackV,6)) {
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
        if(fIsPythia){
            Int_t corrpartidx =-1;
            double ppt;
            //if(fDoMCCorrection) TrackWeight *= GetMonteCarloCorrectionFactor(trackV,corrpartidx,ppt);
        }
        FillTrackHistograms(trackV,dca,cov,TrackWeight);
    }
    AliJetContainer *  jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    //fAnalysisCuts[bAnalysisCut_MinJetPt]=jetconrec->GetJetPtCut();
    //fAnalysisCuts[bAnalysisCut_MaxJetPt]=jetconrec->GetJetPtCutMax();
    //fAnalysisCuts[bAnalysisCut_MinJetEta]=jetconrec->GetMinEta();
    //fAnalysisCuts[bAnalysisCut_MaxJetEta]=jetconrec->GetMaxEta();

    //printf("In Program %f < jetpt <%f, %f < jeteta < %f\n",fAnalysisCuts[bAnalysisCut_MinJetPt],fAnalysisCuts[bAnalysisCut_MaxJetPt],fAnalysisCuts[bAnalysisCut_MinJetEta],fAnalysisCuts[bAnalysisCut_MaxJetEta] );
    FillHist("fh1dNoJetsPerEvent",jetconrec->GetNJets(),1);

    if (!jetconrec) return kFALSE;
    AliJetContainer * jetcongen = nullptr;
    AliEmcalJet * jetgen  = nullptr;
    if(fIsPythia){
        jetcongen = static_cast<AliJetContainer*>(fJetCollArray.At(1));
        if(!MatchJetsGeometricDefault()) AliInfo("Jet matching did not succeed!");
        jetcongen->ResetCurrentID();
        while ((jetgen = jetcongen->GetNextAcceptJet()))
        {
            if (!jetgen) continue;
            //Int_t jetflavour =0;
            //Bool_t is_udgjet = kFALSE;
            /*jetflavour =IsMCJetPartonFast(jetgen,fJetRadius,is_udgjet);
            FillHist("fh1dJetGenPt",GetPtCorrectedMC(jetgen), 1); //this->fXsectionWeightingFactor);
            if(jetflavour ==0)      FillHist("fh1dJetGenPtUnidentified",GetPtCorrectedMC(jetgen), 1); // this->fXsectionWeightingFactor );
            else if(jetflavour ==1) FillHist("fh1dJetGenPtudsg",GetPtCorrectedMC(jetgen), 1);   //this->fXsectionWeightingFactor );
            else if(jetflavour ==2) FillHist("fh1dJetGenPtc",GetPtCorrectedMC(jetgen), 1);  //this->fXsectionWeightingFactor );
            else if(jetflavour ==3) FillHist("fh1dJetGenPtb",GetPtCorrectedMC(jetgen), 1);  //this->fXsectionWeightingFactor );
             else if(jetflavour ==4) FillHist("fh1dJetGenPts",GetPtCorrectedMC(jetgen), 1);  //this->fXsectionWeightingFactor );*/
        }
        jetcongen->ResetCurrentID();
        jetconrec->ResetCurrentID();
    }

    // Loop over reconstructed/matched jets for template creation and analysis
    AliEmcalJet * jetrec  = nullptr;
    AliEmcalJet * jetmatched  = nullptr;
    jetconrec->ResetCurrentID();
    Double_t jetpt=0;
    while ((jetrec = jetconrec->GetNextAcceptJet()))
    {//start jetloop
        if(!jetrec) continue;
        jetpt = jetrec->Pt();
        if((!(jetconrec->GetRhoParameter() == nullptr))&&fDoUnderlyingEventSub)
        {
            printf("Correct for Underlying Event.\n");
            jetpt = jetpt - jetconrec->GetRhoVal() * jetrec->Area();
        }
        if(fIsPythia){
            if (jetrec->MatchedJet()) {
                Double_t genpt = jetrec->MatchedJet()->Pt();
                if((!(jetcongen->GetRhoParameter() == nullptr))&&fDoUnderlyingEventSub)
                {
                    printf("Correct for Underlying Event.\n");
                    genpt = genpt - jetcongen->GetRhoVal() * jetrec->MatchedJet()->Area();
                }
                FillHist("fh2dJetGenPtVsJetRecPt",genpt,jetpt,1);    // this->fXsectionWeightingFactor );
            }
        }
        FillHist("fh1dJetArea",jetrec->Area(),1);

        Double_t dca[2] = {-99999,-99999};
        Double_t cov[3] = {-99999,-99999,-99999};
        Double_t sign=0;
        Int_t jetflavour=0;
        Bool_t is_udgjet = kFALSE;
        	if(fIsPythia){
                  jetmatched = nullptr;
                  jetmatched =jetrec->MatchedJet();
                  if(jetmatched){
                    jetflavour = IsMCJetPartonFast(jetmatched,fJetRadius,is_udgjet); //Event based association to save memory
                  }
                  else{
                    jetflavour=0;
                  }

                }
                FillHist("fh1dJetRecPt",jetpt, 1);  //this->fXsectionWeightingFactor );
                if(fIsPythia){
                    if(jetflavour==0)     FillHist("fh1dJetRecPtUnidentified",jetpt, 1);    //this->fXsectionWeightingFactor );
                    else if(jetflavour==1)FillHist("fh1dJetRecPtudsg",        jetpt, 1);    //this->fXsectionWeightingFactor );
                    else if(jetflavour==2)FillHist("fh1dJetRecPtc",           jetpt, 1);    //this->fXsectionWeightingFactor );
                    else if(jetflavour==3)FillHist("fh1dJetRecPtb",           jetpt, 1);    //this->fXsectionWeightingFactor );
                    else if(jetflavour==4)FillHist("fh1dJetRecPts",           jetpt, 1);    //this->fXsectionWeightingFactor );
                }

                RecursiveParents(jetrec, jetconrec);
               

                FillHist("fh1dJetRecEtaPhiAccepted",jetrec->Eta(),jetrec->Phi(), 1);   //this->fXsectionWeightingFactor );
                FillHist("fh1dJetRecPtAccepted",jetpt, 1);  //this->fXsectionWeightingFactor );
                if(fIsPythia){

                    if(jetflavour==0)     FillHist("fh1dJetRecPtUnidentifiedAccepted",jetpt,1);  //this->fXsectionWeightingFactor );
                    else if(jetflavour==1)FillHist("fh1dJetRecPtudsgAccepted",        jetpt,1);  //this->fXsectionWeightingFactor );
                    else if(jetflavour==2)FillHist("fh1dJetRecPtcAccepted",           jetpt,1);  //this->fXsectionWeightingFactor );
                    else if(jetflavour==3)FillHist("fh1dJetRecPtbAccepted",           jetpt,1);  //this->fXsectionWeightingFactor );
                    else if(jetflavour==4)FillHist("fh1dJetRecPtsAccepted",           jetpt,1);  //this->fXsectionWeightingFactor );
                    //GetOutOfJetParticleComposition(jetrec,jetflavour);
                }
                std::vector<SJetIpPati> sImpParXY,sImpParXYZ,sImpParXYSig,sImpParXYZSig;

                double jetprob =  -1;
                if(fDoJetProb){
                    jetprob = CalculateJetProb(jetrec);
                    const char * subtype_jp [4] = {"","udsg","c","b"};
                    if(jetflavour>0 && fIsPythia){
                        if(jetflavour==1) FillHist("fh2d_jetprob_light",jetpt,jetprob,1);  //this->fXsectionWeightingFactor );
                        else if(jetflavour==2) FillHist("fh2d_jetprob_charm",jetpt,jetprob,1);  //this->fXsectionWeightingFactor );
                        else if(jetflavour==3) FillHist("fh2d_jetprob_beauty",jetpt,jetprob,1);    //*this->fXsectionWeightingFactor ); );
                        if(jetprob >0.5)FillHist(Form("fh1dJetRecPt_0_5JP_%sAccepted",subtype_jp[jetflavour]),jetpt,1);  //this->fXsectionWeightingFactor );
                        if (jetprob >0.6)FillHist(Form("fh1dJetRecPt_0_6JP_%sAccepted",subtype_jp[jetflavour]),jetpt,1);  //this->fXsectionWeightingFactor );
                        if (jetprob >0.7)FillHist(Form("fh1dJetRecPt_0_7JP_%sAccepted",subtype_jp[jetflavour]),jetpt,1);  //this->fXsectionWeightingFactor );
                        if (jetprob >0.8)FillHist(Form("fh1dJetRecPt_0_8JP_%sAccepted",subtype_jp[jetflavour]),jetpt,1);  //this->fXsectionWeightingFactor );
                        if (jetprob >0.9)FillHist(Form("fh1dJetRecPt_0_9JP_%sAccepted",subtype_jp[jetflavour]),jetpt,1);  //this->fXsectionWeightingFactor );
                        if (jetprob >0.95)FillHist(Form("fh1dJetRecPt_0_95JP_%sAccepted",subtype_jp[jetflavour]),jetpt,1);  //this->fXsectionWeightingFactor );
                    }
                    if(jetprob >0.5)FillHist(Form("fh1dJetRecPt_0_5JP_%sAccepted","all"),jetpt,1);  //this->fXsectionWeightingFactor );
                    if (jetprob >0.6)FillHist(Form("fh1dJetRecPt_0_6JP_%sAccepted","all"),jetpt,1);  //this->fXsectionWeightingFactor );
                    if (jetprob >0.7)FillHist(Form("fh1dJetRecPt_0_7JP_%sAccepted","all"),jetpt,1);  //this->fXsectionWeightingFactor );
                    if (jetprob >0.8)FillHist(Form("fh1dJetRecPt_0_8JP_%sAccepted","all"),jetpt,1);  //this->fXsectionWeightingFactor );
                    if (jetprob >0.9)FillHist(Form("fh1dJetRecPt_0_9JP_%sAccepted","all"),jetpt,1);  //this->fXsectionWeightingFactor );
                    if (jetprob >0.95)FillHist(Form("fh1dJetRecPt_0_95JP_%sAccepted","all"),jetpt,1);  //this->fXsectionWeightingFactor );
                }
                 AliVParticle *vp=0x0;

                int NJetParticles=0;  //Used for counting particles per jet
                for(UInt_t i = 0; i < jetrec->GetNumberOfTracks(); i++) {//start trackloop
                    TrackWeight=1;
                    Double_t xyzatcda[3];

                    vp = static_cast<AliVParticle*>(jetrec->TrackAt(i, jetconrec->GetParticleContainer()->GetArray()));
                    if (!vp){
                      Printf("ERROR: AliVParticle associated to constituent not found");
                      continue;
                    }

                    AliVTrack *vtrack = dynamic_cast<AliVTrack*>(vp);
                    if (!vtrack) {
                      printf("ERROR: Could not receive track%d\n", i);
                      continue;
                    }

                    AliAODTrack *trackV = dynamic_cast<AliAODTrack*>(vtrack);

                    if (!trackV || !jetrec)            continue;

                    if (!IsTrackAccepted((AliAODTrack*)trackV,3))   continue;
                    ++NJetParticles;


                    if(GetImpactParameterWrtToJet((AliAODTrack*)trackV,(AliAODEvent*)InputEvent(),jetrec,dca,cov,xyzatcda,sign)){
                        if(fEventVertex) {
                            delete fEventVertex;
                            fEventVertex =nullptr;
                        }


                        Int_t corridx=-1;double ppt;
                        //(fIsPythia&&fDoMCCorrection) ? TrackWeight = GetMonteCarloCorrectionFactor(trackV,corridx,ppt) : TrackWeight =1;
                        dca[0]=fabs(dca[0]);
                      //  Double_t cursImParXY     =TMath::Abs(GetValImpactParameter(   kXY,dca,cov))*sign;
                        Double_t cursImParXYSig  =TMath::Abs(GetValImpactParameter(kXYSig,dca,cov))*sign;
                       // Double_t cursImParXYZ    =TMath::Abs(GetValImpactParameter(   kXYZ,dca,cov))*sign;
                        Double_t cursImParXYZSig =TMath::Abs(GetValImpactParameter(kXYZSig,dca,cov))*sign;

                      if(fIsPythia){
                        if(is_udgjet){
                            if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                                FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShits",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                            }
                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                                FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShits",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                            }
                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                                FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShits",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                            }
                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                                FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShits",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                            }
                        }
                        if(IsFromElectron((AliAODTrack*)trackV)){
                            if(is_udgjet){
                                if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                            }
                        }
                        else if(IsFromPion((AliAODTrack*)trackV)){
                            if(is_udgjet){
                                if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                            }
                        }
                        else if(IsFromKaon((AliAODTrack*)trackV)){
                            if(is_udgjet){
                                if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                            }
                        }
                        else if(IsFromProton((AliAODTrack*)trackV)){
                            if(is_udgjet){
                                if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                            }
                        }
                      }


                        if (!IsTrackAccepted((AliAODTrack*)trackV,6))   continue;

                            //Fill jet probability ipsig histograms for template fitting
                        const char * subtype_jp [4] = {"","udsg","c","b"};

                        
                        if(fDoJetProb){
                            if(jetflavour>0 && fIsPythia ){
                                if(jetprob >0.5){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_5JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_5JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                if (jetprob >0.6){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_6JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_6JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                if (jetprob >0.7){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_7JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_7JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                if (jetprob >0.8){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_8JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_8JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                if (jetprob >0.9){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_9JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_9JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                                if (jetprob >0.95){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_95JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_95JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                }
                            }

                            if(jetprob >0.5){
                                FillHist(Form("fh2d_ImpSigXY_%s_0_5JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                FillHist(Form("fh2d_ImpSigXYZ_%s_0_5JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                            }
                            if (jetprob >0.6){
                                FillHist(Form("fh2d_ImpSigXY_%s_0_6JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                FillHist(Form("fh2d_ImpSigXYZ_%s_0_6JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                            }
                            if (jetprob >0.7){
                                FillHist(Form("fh2d_ImpSigXY_%s_0_7JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                FillHist(Form("fh2d_ImpSigXYZ_%s_0_7JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                            }
                            if (jetprob >0.8){
                                FillHist(Form("fh2d_ImpSigXY_%s_0_8JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                FillHist(Form("fh2d_ImpSigXYZ_%s_0_8JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                            }
                            if (jetprob >0.9){
                                FillHist(Form("fh2d_ImpSigXY_%s_0_9JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                FillHist(Form("fh2d_ImpSigXYZ_%s_0_9JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                            }
                            if (jetprob >0.95){
                                FillHist(Form("fh2d_ImpSigXY_%s_0_95JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                                FillHist(Form("fh2d_ImpSigXYZ_%s_0_95JP","all"),jetpt,cursImParXYSig,TrackWeight);     //*this->fXsectionWeightingFactor );
                            }
                        }
                    }

                    if(GetImpactParameterWrtToJet((AliAODTrack*)trackV,(AliAODEvent*)InputEvent(),jetrec,dca,cov,xyzatcda,sign)){
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
                        SJetIpPati a(cursImParXY, TrackWeight,kFALSE,kFALSE,corridx,trackV->Pt()); sImpParXY.push_back(a);
                        SJetIpPati b(cursImParXYZ, TrackWeight,kFALSE,kFALSE,corridx,trackV->Pt()); sImpParXYZ.push_back(b);
                        SJetIpPati c(cursImParXYSig, TrackWeight,kFALSE,kFALSE,corridx,trackV->Pt());sImpParXYSig.push_back(c);
                        SJetIpPati d(cursImParXYZSig, TrackWeight,kFALSE,kFALSE,corridx,trackV->Pt());sImpParXYZSig.push_back(d);
                        //printf("curImParXY=%f, TrackWeight=%f,  corridx=%i, pt=%f\n",cursImParXYSig,TrackWeight,corridx, trackV->Pt());

                    }
                }//end trackloop

                FillHist("fh1dParticlesPerJet",NJetParticles,1);

                std::sort(sImpParXY.begin(),sImpParXY.end(),        AliAnalysisTaskHFJetIPQA::mysort);
                std::sort(sImpParXYSig.begin(),sImpParXYSig.end(),  AliAnalysisTaskHFJetIPQA::mysort);
                std::sort(sImpParXYZ.begin(),sImpParXYZ.end(),      AliAnalysisTaskHFJetIPQA::mysort);
                std::sort(sImpParXYZSig.begin(),sImpParXYZSig.end(),AliAnalysisTaskHFJetIPQA::mysort);
                const char * subtype[5] = {"Unidentified","udsg","c","b","s"};
                const char * subord [3] = {"First","Second","Third"};
                const char * stype  [4] = {"fh2dJetSignedImpParXY","fh2dJetSignedImpParXYSignificance","fh2dJetSignedImpParXYZ","fh2dJetSignedImpParXYZSignificance"};

                bool hasIPs [3] ={kFALSE,kFALSE,kFALSE};
                if((int)sImpParXYSig.size()>0) hasIPs[0]=kTRUE;
                if((int)sImpParXYSig.size()>1) hasIPs[1]=kTRUE;
                if((int)sImpParXYSig.size()>2) hasIPs[2]=kTRUE;

                Double_t ipval [3] = {-9999};
                if(hasIPs[0])ipval[0] =sImpParXYSig.at(0).first;
                if(hasIPs[1])ipval[1] =sImpParXYSig.at(1).first;
                if(hasIPs[2])ipval[2] =sImpParXYSig.at(2).first;
                //if(hasIPs[0])printf("N=1: cursImParXY=%f, TrackWeight=%f,corridx=%i, pt=%f\n",sImpParXYSig.at(0).first, sImpParXYSig.at(0).second, sImpParXYSig.at(0).trackLabel, sImpParXYSig.at(0).trackpt);
                //if(hasIPs[1])printf("N=2: cursImParXY=%f, TrackWeight=%f, corridx=%i, pt=%f\n",sImpParXYSig.at(1).first, sImpParXYSig.at(1).second, sImpParXYSig.at(1).trackLabel, sImpParXYSig.at(1).trackpt);
                //if(hasIPs[2])printf("N=3: cursImParXY=%f, TrackWeight=%f, corridx=%i, pt=%f\n",sImpParXYSig.at(2).first, sImpParXYSig.at(2).second, sImpParXYSig.at(2).trackLabel, sImpParXYSig.at(2).trackpt);
                //printf("*********************************************************\n");

                if(hasIPs[0])FillHist("fh1dTrackPt_n_1_all_Accepted",sImpParXYSig.at(0).trackpt,1);
                if(hasIPs[1])FillHist("fh1dTrackPt_n_2_all_Accepted",sImpParXYSig.at(1).trackpt,1);
                if(hasIPs[2])FillHist("fh1dTrackPt_n_3_all_Accepted",sImpParXYSig.at(2).trackpt,1);

                if(fFillCorrelations || fUseTreeForCorrelations){
                    FillCorrelations(hasIPs,ipval,jetpt);
                    if(fFillCorrelations && ! fUseTreeForCorrelations){
                        if(!fIsMixSignalReady_n1 && hasIPs[0]) SetMixDCA(1,ipval[0] );
                        if(!fIsMixSignalReady_n2 && hasIPs[1]) SetMixDCA(2,ipval[1] );
                        if(!fIsMixSignalReady_n3 && hasIPs[2]) SetMixDCA(3,ipval[2] );
                    }
                }


                if(sImpParXY.size()!=0){
                  FillHist("fh2dNoAcceptedTracksvsJetArea",(int)sImpParXY.size(),jetrec->Area(),1);
                }

                for (Int_t ot = 0 ; ot <3 ;++ot){
                    if ((int)sImpParXY.size()>ot){

                        if(ot==0) {
                            if(jetflavour >0){
                                FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,subtype[jetflavour]),jetpt,1);     //*this->fXsectionWeightingFactor );
                            }
                            FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,"all"),jetpt,1);     //*this->fXsectionWeightingFactor );
                        }
                        else if (ot==1)
                        {
                            if(jetflavour >0){
                                FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,subtype[jetflavour]),jetpt,1);     //*this->fXsectionWeightingFactor );
                            }
                            FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,"all"),jetpt,1);     //*this->fXsectionWeightingFactor );
                        }
                        else if (ot==2){
                            if(jetflavour >0){
                                FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,subtype[jetflavour]),jetpt,1);     //*this->fXsectionWeightingFactor );
                            }
                            FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,"all"),jetpt,1);     //*this->fXsectionWeightingFactor );
                        }


                        Double_t params [4] ={sImpParXY.at(ot).first,sImpParXYSig.at(ot).first,sImpParXYZ.at(ot).first,sImpParXYZSig.at(ot).first};
                        Double_t weights[4] ={sImpParXY.at(ot).second,sImpParXYSig.at(ot).second,sImpParXYZ.at(ot).second,sImpParXYZSig.at(ot).second};
                        Int_t    correctionwindex[4] ={sImpParXY.at(ot).trackLabel,sImpParXYSig.at(ot).trackLabel,sImpParXYZ.at(ot).trackLabel,sImpParXYZSig.at(ot).trackLabel};

                        for (Int_t ost = 0 ; ost <4 ;++ost){
                            TString hname = Form("%s%s",stype[ost],subord[ot]);
                            if(fIsPythia)   FillHist(hname.Data(),jetpt,params[ost],weights[ost] *  1);     //*this->fXsectionWeightingFactor );
                            else  FillHist(hname.Data(),jetpt,params[ost], this->fXsectionWeightingFactor );


                        }
                        if(fIsPythia){

                                        if(ot ==0){//N=1
                                          FillHist("fh2dNMCWeightSpeciesPerJetPtN1_IP_all",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                          FillHist("fh2dNMCWeightSpeciesPerJetPtN1_SIP_all",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                            if(jetflavour==3) {//b
                                                FillHist("fh2dNMCWeightSpeciesPerJetPtN1_IP_b",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                FillHist("fh2dNMCWeightSpeciesPerJetPtN1_SIP_b",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                            }
                                        else   if(jetflavour==2) {//c
                                            FillHist("fh2dNMCWeightSpeciesPerJetPtN1_IP_c",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                            FillHist("fh2dNMCWeightSpeciesPerJetPtN1_SIP_c",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                        }
                                        else   if(jetflavour==1) {//lf
                                            FillHist("fh2dNMCWeightSpeciesPerJetPtN1_IP_lf",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                            FillHist("fh2dNMCWeightSpeciesPerJetPtN1_SIP_lf",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                        }

                                    }
                                        else  if(ot ==1){//N=2
                                            FillHist("fh2dNMCWeightSpeciesPerJetPtN2_IP_all",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                            FillHist("fh2dNMCWeightSpeciesPerJetPtN2_SIP_all",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                if(jetflavour==3) {//b
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN2_IP_b",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN2_SIP_b",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                }
                                                else   if(jetflavour==2) {//c
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN2_IP_c",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN2_SIP_c",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                }
                                                else   if(jetflavour==1) {//lf
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN2_IP_lf",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN2_SIP_lf",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                }

                                            }
                                        else  if(ot==2){//N=3
                                            FillHist("fh2dNMCWeightSpeciesPerJetPtN3_IP_all",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                            FillHist("fh2dNMCWeightSpeciesPerJetPtN3_SIP_all",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                if(jetflavour==3) {//b
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN3_IP_b",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN3_SIP_b",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                }
                                                else   if(jetflavour==2) {//c
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN3_IP_c",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN3_SIP_c",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                }
                                                else   if(jetflavour==1) {//lf
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN3_IP_lf",correctionwindex[0]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dNMCWeightSpeciesPerJetPtN3_SIP_lf",correctionwindex[1]+0.5,jetpt, 1);     //*this->fXsectionWeightingFactor );
                                                }
                                            }
                                            for (Int_t ost = 0 ; ost <4 ;++ost){
                                                TString hname = Form("%s%s%s",stype[ost],subtype[jetflavour],subord[ot]);
                                                FillHist(hname.Data(),jetpt,params[ost],weights[ost]* 1);     //*this->fXsectionWeightingFactor );
                                            }
                                        }
                                    }
                                }
                                sImpParXY.clear();
                                sImpParXYSig.clear();
                                sImpParXYZ.clear();
                                sImpParXYZSig.clear();
                            }//end jetloop

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
                        Bool_t AliAnalysisTaskHFJetIPQA::IsVertexSelected(const AliVVertex *vertex){
                            if(!vertex) return kFALSE;
    //Printf("GetNContributors %i",((AliAODVertex*)vertex)->GetNContributors());

                            if(((AliAODVertex*)vertex)->GetNContributors()<1) {
                                return kFALSE;
                            }

    ///Printf("GetNContributors %i",((AliAODVertex*)vertex)->GetNContributors());
                            if(((AliAODVertex*)vertex)->GetNContributors()<(int)(fAnalysisCuts[bAnalysisCut_NContibutors])) {
                                return kFALSE;
                            }
                            if(TMath::Abs(((AliAODVertex*)vertex)->GetZ())>=fAnalysisCuts[bAnalysisCut_MaxVtxZ]) {
                                return kFALSE;
                            }
                            if((TMath::Abs(((AliAODVertex*)vertex)->GetX() - ((AliAODEvent*)InputEvent())->GetDiamondX()) > fAnalysisCuts[bAnalysisCut_SigmaDiamond]*sqrt(((AliAODEvent*)InputEvent())->GetSigma2DiamondX()))) {
                                return kFALSE;
                            }
                            if((TMath::Abs(((AliAODVertex*)vertex)->GetY() - ((AliAODEvent*)InputEvent())->GetDiamondY()) > fAnalysisCuts[bAnalysisCut_SigmaDiamond]*sqrt(((AliAODEvent*)InputEvent())->GetSigma2DiamondY()))) {
                                return kFALSE;
                            }
                            if((TMath::Abs(((AliAODVertex*)vertex)->GetZ() - ((AliAODEvent*)InputEvent())->GetDiamondZ()) > fAnalysisCuts[bAnalysisCut_SigmaDiamond]*sqrt(((AliAODEvent*)InputEvent())->GetSigma2DiamondZ()))) {
                                return kFALSE;
                            }
                            return kTRUE;
                        }


                        Bool_t AliAnalysisTaskHFJetIPQA::IsSelected(AliVEvent *event, Int_t &WhyRejected,ULong_t &RejectionBits){
                            WhyRejected =0;
                            Bool_t accept=kTRUE;
                            RejectionBits=000;
                            Bool_t isSelected = kFALSE;
                            isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny);

                            if(!isSelected) {
                                WhyRejected=kPhysicsSelection;
                                return kFALSE;
                            }
                            if(!fEventVertex) {
                                WhyRejected=kNoVertex;
                                return kFALSE;
                            }
                            if(fEventVertex->GetNContributors()<1) {
                                WhyRejected=kNoContributors;
                                return kFALSE;
                            }

                            if(fEventVertex->GetChi2perNDF()>fAnalysisCuts[bAnalysisCut_Z_Chi2perNDF]) {
                                WhyRejected=kVertexChi2NDF;
                                return kFALSE;
                            }
                            if(fEventVertex->GetNContributors()<(int)(fAnalysisCuts[bAnalysisCut_NContibutors])) {
                                WhyRejected=kTooFewVtxContrib;
                                return kFALSE;
                            }

                            if(TMath::Abs(fEventVertex->GetZ())>=fAnalysisCuts[bAnalysisCut_MaxVtxZ]) {
                                WhyRejected=kZVtxOutFid;
                                return kFALSE;
                            }
                            return accept;
                        }
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
    if(!fPidResponse)   fPidResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPidResponse) {
      AliFatal("NULL PID response");
    }
    if(!fCombined) fCombined = new AliPIDCombined();
  }
    //Graphs currently not in use
    /*                                            const Int_t gfProtonN = 9;
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

*/

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

  fh1dCuts =(TH1D*)AddHistogramm("fh1dCuts","",3,0,3);

  fh1dTracksAccepeted =(TH1D*)AddHistogramm("fh1dTracksAccepeted","# tracks before/after cuts;;",3,0,3);
  fh1dTracksAccepeted->GetXaxis()->SetBinLabel(1,"total");
  fh1dTracksAccepeted->GetXaxis()->SetBinLabel(2,"accepted");
  fh1dTracksAccepeted->GetXaxis()->SetBinLabel(3,"rejected");

  const Int_t dimSpec   = 6;
  const Int_t nBinsSpec[6]     = {50,100,100,20,100,2};
  const Double_t lowBinSpec[6] = {0.,-10,0,0,0,0};
  const Double_t hiBinSpec[6]  = {5.,10.,100,20,100,2};
  fHLundIterative = new THnSparseF("fHLundIterative",
                  "LundIterativePlot [log(1/theta),log(z*theta),pTjet,algo]",
                  dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
  fOutput->Add(fHLundIterative);

  //Histograms for vertexing factor quicktest
  //fHistManager.CreateTH1("fh1dVERTEXFACTOR_VERTEXZ_FULL","fh1dVERTEXFACTOR_VERTEXZ_FULL;;",400,-100,100,"s");
  //fHistManager.CreateTH1("fh1dVERTEXFACTOR_VERTEXZ_10CM","fh1dVERTEXFACTOR_VERTEXZ_10CM;",400,-100,100,"s");

  //TH1D * hvt = (TH1D*)AddHistogramm("fh1dVERTEXFACTOR_EVENTS","fh1dVERTEXFACTOR_EVENTS;;",3,0,3);
  //hvt->GetXaxis()->SetBinLabel(1,"Event in");
  //hvt->GetXaxis()->SetBinLabel(2,"Event kMB");
  //hvt->GetXaxis()->SetBinLabel(3,"Event kMB with Vertex");

  fHistManager.CreateTH1("fh1dNoParticlesPerEvent","fh1dNoParticlesvsEvent;#;No Particles/Event",5000, 0, 5000,"s");
  fHistManager.CreateTH1("fh1dNoJetsPerEvent","fh1dNoJetsPerEvent;#;No Jets/Event",400, 0, 100,"s");
  fHistManager.CreateTH1("fh1dEventsAcceptedInRun","fh1dEventsAcceptedInRun;Events Accepted;count",1,0,1,"s");

  fh2dLightNotContrib=(TH2D*)AddHistogramm("fh2dLightNotContrib",";particle pt (GeV/c);jet pt (GeV/c)",500,0,100,250,0,250);
  fh2dCharmNotContrib=(TH2D*)AddHistogramm("fh2dCharmNotContrib",";particle pt (GeV/c);jet pt (GeV/c)",500,0,100,250,0,250);
  fh2dBottomNotContrib=(TH2D*)AddHistogramm("fh2dBottomNotContrib",";particle pt (GeV/c);jet pt (GeV/c)",500,0,100,250,0,250);
  fh2dManifoldParton=(TH2D*)AddHistogramm("fh2dManifoldParton",";Number of Flavours per Jet;jet pt GeV/c",50,0,50,250,0,250);
  fh2dLightNMatch=(TH2D*)AddHistogramm("fh2dLightNMatch",";Number of particles within jet;particle pt (GeV/c)",50,0,50,200,0,50);
  fh2dCharmNMatch=(TH2D*)AddHistogramm("fh2dCharmNMatch",";Number of particles within jet;particle pt (GeV/c)",50,0,50,200,0,50);
  fh2dBottomNMatch=(TH2D*)AddHistogramm("fh2dBottomNMatch",";Number of particles within jet;particle pt (GeV/c)",50,0,50,200,0,50);
  fh2dLightDeta=(TH2D*)AddHistogramm("fh2dLightDeta",";particle pt (GeV/c);#Delta R",200,0,50,200,0,10);
  fh2dCharmDeta=(TH2D*)AddHistogramm("fh2dCharmDeta",";particle pt (GeV/c);#Delta R",200,0,50,200,0,10);
  fh2dBottomDeta=(TH2D*)AddHistogramm("fh2dBottomDeta",";particle pt (GeV/c);#Delta R",200,0,50,200,0,10);

  fOutput->Add(fh2dLightNotContrib);
  fOutput->Add(fh2dCharmNotContrib);
  fOutput->Add(fh2dBottomNotContrib);
  fOutput->Add(fh2dLightNMatch);
  fOutput->Add(fh2dCharmNMatch);
  fOutput->Add(fh2dBottomNMatch);
  fOutput->Add(fh2dManifoldParton);
  fOutput->Add(fh2dLightDeta);
  fOutput->Add(fh2dCharmDeta);
  fOutput->Add(fh2dBottomDeta);




  //In-Out of jet particle composition in pythia events
  /*if(fIsPythia){
    fHistManager.CreateTH2("fh2dParticleSpectra_InCone","fh2dParticleSpectra_InCone ;species i;p_t particle (GeV/c)",22,0,22,200,0,50,"s");
    fHistManager.CreateTH2("fh2dParticleSpectra_InCone_bjet","fh2dParticleSpectra_InCone_bjet ;species i;p_t particle (GeV/c)",22,0,22,200,0,50,"s");
    fHistManager.CreateTH2("fh2dParticleSpectra_InCone_cjet","fh2dParticleSpectra_InCone_cjet ;species i;p_t particle (GeV/c)",22,0,22,200,0,50,"s");
    fHistManager.CreateTH2("fh2dParticleSpectra_InCone_lfjet","fh2dParticleSpectra_InCone_lfjet ;species i;p_t particle (GeV/c)",22,0,22,200,0,50,"s");
    fHistManager.CreateTH2("fh2dParticleSpectra_OutOfCone","fh2dParticleSpectra_OutOfCone ;species i;p_t particle (GeV/c)",22,0,22,200,0,50,"s");
    fHistManager.CreateTH2("fh2dParticleSpectra_Event","fh2dParticleSpectra_Event ;species i;p_t particle (GeV/c)",22,0,22,200,0,50,"s");
  }*/
  //MC weights
  /*fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN1_IP_all","fh2dNMCWeightSpeciesPerJetPtN1_all (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN1_IP_b","fh2dNMCWeightSpeciesPerJetPtN1_b (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN1_IP_c","fh2dNMCWeightSpeciesPerJetPtN1_c (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN1_IP_lf","fh2dNMCWeightSpeciesPerJetPtN1_lf (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN2_IP_all","fh2dNMCWeightSpeciesPerJetPtN2_all (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN2_IP_b","fh2dNMCWeightSpeciesPerJetPtN2_b (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN2_IP_c","fh2dNMCWeightSpeciesPerJetPtN2_c (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN2_IP_lf","fh2dNMCWeightSpeciesPerJetPtN2_lf (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN3_IP_all","fh2dNMCWeightSpeciesPerJetPtN3_all (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN3_IP_b","fh2dNMCWeightSpeciesPerJetPtN3_b (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN3_IP_c","fh2dNMCWeightSpeciesPerJetPtN3_c (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN3_IP_lf","fh2dNMCWeightSpeciesPerJetPtN3_lf (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN1_SIP_all","fh2dNMCWeightSpeciesPerJetPtN1_all (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN1_SIP_b","fh2dNMCWeightSpeciesPerJetPtN1_b (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN1_SIP_c","fh2dNMCWeightSpeciesPerJetPtN1_c (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN1_SIP_lf","fh2dNMCWeightSpeciesPerJetPtN1_lf (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN2_SIP_all","fh2dNMCWeightSpeciesPerJetPtN2_all (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN2_SIP_b","fh2dNMCWeightSpeciesPerJetPtN2_b (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN2_SIP_c","fh2dNMCWeightSpeciesPerJetPtN2_c (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN2_SIP_lf","fh2dNMCWeightSpeciesPerJetPtN2_lf (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN3_SIP_all","fh2dNMCWeightSpeciesPerJetPtN3_all (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN3_SIP_b","fh2dNMCWeightSpeciesPerJetPtN3_b (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN3_SIP_c","fh2dNMCWeightSpeciesPerJetPtN3_c (N used weights) ;species i;pt",22,0,22,150,0,150,"s");
  fHistManager.CreateTH2("fh2dNMCWeightSpeciesPerJetPtN3_SIP_lf","fh2dNMCWeightSpeciesPerJetPtN3_lf (N used weights) ;species i;pt",22,0,22,150,0,150,"s");*/
                                                            
  //Check N1N2N3 correlations
  /*if(fUseTreeForCorrelations){
    Printf("Adding Tree to output file");
    OpenFile(1);
    fCorrelationCrossCheck = new TTree("fCorrelationCrossCheck","fCorrelationCrossCheck");
    fCorrelationCrossCheck->Branch("n1",&fTREE_n1,"px/F");
    fCorrelationCrossCheck->Branch("n2",&fTREE_n2,"py/F");
    fCorrelationCrossCheck->Branch("n3",&fTREE_n3,"pz/F");
    fCorrelationCrossCheck->Branch("pt",&fTREE_pt,"pz/F");
    fOutput->Add(fCorrelationCrossCheck);
  }*/

  /*if (fFillCorrelations && !fUseTreeForCorrelations){
    fHistManager.CreateTH2("fh2dInclusiveCorrelationN1N2","fh2dInclusiveCorrelationN1N2 ;N1 impact parameter xy (cm);N2impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dInclusiveCorrelationN1N3","fh2dInclusiveCorrelationN1N3 ;N1 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dInclusiveCorrelationN2N3","fh2dInclusiveCorrelationN2N3 ;N2 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater10_20GeVCorrelationN1N2","fh2dGreater10_20GeVCorrelationN1N2 ;N1 impact parameter xy (cm);N2impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater10_20GeVCorrelationN1N3","fh2dGreater10_20GeVCorrelationN1N3 ;N1 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater10_20GeVCorrelationN2N3","fh2dGreater10_20GeVCorrelationN2N3 ;N2 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater20_30GeVCorrelationN1N2","fh2dGreater20_30GeVCorrelationN1N2 ;N1 impact parameter xy (cm);N2impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater20_30GeVCorrelationN1N3","fh2dGreater20_30GeVCorrelationN1N3 ;N1 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater20_30GeVCorrelationN2N3","fh2dGreater20_30GeVCorrelationN2N3 ;N2 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater30_100GeVCorrelationN1N2","fh2dGreater30_100GeVCorrelationN1N2 ;N1 impact parameter xy (cm);N2impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater30_100GeVCorrelationN1N3","fh2dGreater30_100GeVCorrelationN1N3 ;N1 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater30_100GeVCorrelationN2N3","fh2dGreater30_100GeVCorrelationN2N3 ;N2 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dInclusiveCorrelationN1N2mix","fh2dInclusiveCorrelationN1N2mix ;N1 impact parameter xy (cm);N2impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dInclusiveCorrelationN1N3mix","fh2dInclusiveCorrelationN1N3mix ;N1 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dInclusiveCorrelationN2N3mix","fh2dInclusiveCorrelationN2N3mix ;N2 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater10_20GeVCorrelationN1N2mix","fh2dGreater10_20GeVCorrelationN1N2mix ;N1 impact parameter xy (cm);N2impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater10_20GeVCorrelationN1N3mix","fh2dGreater10_20GeVCorrelationN1N3mix ;N1 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater10_20GeVCorrelationN2N3mix","fh2dGreater10_20GeVCorrelationN2N3mix ;N2 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater20_30GeVCorrelationN1N2mix","fh2dGreater20_30GeVCorrelationN1N2mix ;N1 impact parameter xy (cm);N2impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater20_30GeVCorrelationN1N3mix","fh2dGreater20_30GeVCorrelationN1N3mix ;N1 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater20_30GeVCorrelationN2N3mix","fh2dGreater20_30GeVCorrelationN2N3mix ;N2 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater30_100GeVCorrelationN1N2mix","fh2dGreater30_100GeVCorrelationN1N2mix ;N1 impact parameter xy (cm);N2impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater30_100GeVCorrelationN1N3mix","fh2dGreater30_100GeVCorrelationN1N3mix ;N1 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
    fHistManager.CreateTH2("fh2dGreater30_100GeVCorrelationN2N3mix","fh2dGreater30_100GeVCorrelationN2N3mix ;N2 impact parameter xy (cm);N3impact parameter xy (cm)",1000,-50,50,1000,-50,50,"s");
  }*/

    //Track Impact Parameter Distributions
    fHistManager.CreateTH1("fh1dPtHardMonitor","fh1dPtHardMonitor;ptHard;",500,0,250,"s");
    fHistManager.CreateTH2("fh2dTracksImpParXY","radial imp. parameter ;impact parameter xy (cm);a.u.",2000,lowIPxy,highIPxy,500,0,100.,"s");
    fHistManager.CreateTH2("fh2dTracksImpParXYZ","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-1,1,500,0,100.,"s");
    fHistManager.CreateTH2("fh2dTracksImpParXYZSignificance","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.,"s");
    fHistManager.CreateTH2("fh2dTracksImpParZ","z imp. parameter ;impact parameter xy (cm);a.u.",2000,lowIPxy,highIPxy,500,0,10.,"s");
    fHistManager.CreateTH2("fh2dTracksImpParXYSignificance","radial imp. parameter sig;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.,"s");
    fHistManager.CreateTH2("fh2dTracksImpParZSignificance","z imp. parameter ;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.,"s");
    fHistManager.CreateTH1("fh1dTracksImpParXY","2d imp. parameter ;impact parameter 2d (cm);a.u.",400,-0.2,0.2,"s");
    fHistManager.CreateTH1("fh1dTracksImpParXYZ","3d imp. parameter ;impact parameter 3d (cm);a.u.",2000,0,1.,"s");
    fHistManager.CreateTH1("fh1dTracksImpParXYSignificance","radial imp. parameter ;impact parameter xy significance;a.u.",200,-30,30,"s");
    fHistManager.CreateTH1 ("fh1dTracksImpParXYZSignificance","3d imp. parameter ;impact parameter 3d significance;a.u.",2000,0.,100.,"s");

    //General Jet Properties
    fHistManager.CreateTH2("fh1dJetRecEtaPhiAccepted","detector level jet;#eta;phi",1,-0.5,0.5,1,0.,TMath::TwoPi(),"s");
    fHistManager.CreateTH2("fh2dAcceptedTracksEtaPhi","accepted tracks;#eta;phi",200,-0.9,0.9,200,0.,TMath::TwoPi(),"s");
    fHistManager.CreateTH1("fh1dJetRecPt","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPtAccepted","accepted detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetArea","fh1dJetArea;# Jet Area",100,0,1,"s");
    fHistManager.CreateTH1("fh1dParticlesPerJet","fh1dParticlesPerJet;#, Particles/Jet",100,0,100,"s");
    fHistManager.CreateTH2("fh2dNoAcceptedTracksvsJetArea","fh2dNoAcceptedTracksvsJetArea;No Accepted Tracks;JetArea",20,0,20,100,0,1);




    //IP Distributions for different species for different numbers of ITShits
    if (fIsPythia){
    /*  fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShits", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShits;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShits","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShits;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShits", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShits;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShits","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShits;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShits", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShits;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShits","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShits;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShits", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShits;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShits","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShits;pt (GeV/c); count",200,0,100,500,-30,0,"s");
            //Electrons
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsElectrons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsElectrons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsElectrons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsElectrons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsElectrons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsElectrons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsElectrons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsElectrons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsElectrons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsElectrons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsElectrons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsElectrons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsElectrons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsElectrons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsElectrons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsElectrons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
            //Pions
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsPions", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsPions;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsPions","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsPions;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsPions", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsPions;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsPions","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsPions;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsPions", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsPions;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsPions","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsPions;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsPions", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsPions;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsPions","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsPions;pt (GeV/c); count",200,0,100,500,-30,0,"s");
            //Kaons
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsKaons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsKaons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsKaons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsKaons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsKaons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsKaons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsKaons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsKaons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsKaons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsKaons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsKaons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsKaons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsKaons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsKaons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsKaons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsKaons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
            //Protons
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsProtons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsProtons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsProtons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsProtons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsProtons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsProtons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsProtons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsProtons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsProtons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsProtons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsProtons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsProtons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsProtons", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsProtons;pt (GeV/c); count",200,0,100,500,-30,0,"s");
      fHistManager.CreateTH2("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsProtons","fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsProtons;pt (GeV/c); count",200,0,100,500,-30,0,"s");

      //Jet Probability QA  plots for different particle species
      if(fDoJetProb){
         fHistManager.CreateTH2("fh2d_jetprob_beauty", "fh2d_jetprob_beauty;pt (GeV/c); count",500,0,250,500,0,2.5,"s");
         fHistManager.CreateTH2("fh2d_jetprob_charm", "fh2d_jetprob_charm;pt (GeV/c); count",500,0,250,500,0,2.5,"s");
         fHistManager.CreateTH2("fh2d_jetprob_light", "fh2d_jetprob_light;pt (GeV/c); count",500,0,250,500,0,2.5,"s");
            //Templates for different probabilities
            //50%
         fHistManager.CreateTH2("fh2d_ImpSigXY_b_0_5JP", "fh2d_ImpSigXY_b_0_5JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_b_0_5JP", "fh2d_ImpSigXYZ_b_0_5JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_c_0_5JP", "fh2d_ImpSigXY_c_0_5JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_c_0_5JP", "fh2d_ImpSigXYZ_c_0_5JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_udsg_0_5JP", "fh2d_ImpSigXY_udsg_0_5JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_udsg_0_5JP", "fh2d_ImpSigXYZ_udsg_0_5JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_5JP_bAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_5JP_cAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_5JP_udsgAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
            //60%
         fHistManager.CreateTH2("fh2d_ImpSigXY_b_0_6JP", "fh2d_ImpSigXY_b_0_6JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_b_0_6JP", "fh2d_ImpSigXYZ_b_0_6JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_c_0_6JP", "fh2d_ImpSigXY_c_0_5JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_c_0_6JP", "fh2d_ImpSigXYZ_c_0_6JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_udsg_0_6JP", "fh2d_ImpSigXY_udsg_0_6JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_udsg_0_6JP", "fh2d_ImpSigXYZ_udsg_0_6JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_6JP_bAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_6JP_cAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_6JP_udsgAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
            //70%
         fHistManager.CreateTH2("fh2d_ImpSigXY_b_0_7JP", "fh2d_ImpSigXY_b_0_7JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_b_0_7JP", "fh2d_ImpSigXYZ_b_0_7JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_c_0_7JP", "fh2d_ImpSigXY_c_0_7JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_c_0_7JP", "fh2d_ImpSigXYZ_c_0_7JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_udsg_0_7JP", "fh2d_ImpSigXY_udsg_0_7JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_udsg_0_7JP", "fh2d_ImpSigXYZ_udsg_0_7JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_7JP_bAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_7JP_cAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_7JP_udsgAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
            //80%
         fHistManager.CreateTH2("fh2d_ImpSigXY_b_0_8JP", "fh2d_ImpSigXY_b_0_8JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_b_0_8JP", "fh2d_ImpSigXYZ_b_0_8JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_c_0_8JP", "fh2d_ImpSigXY_c_0_8JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_c_0_8JP", "fh2d_ImpSigXYZ_c_0_8JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_udsg_0_8JP", "fh2d_ImpSigXY_udsg_0_8JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_udsg_0_8JP", "fh2d_ImpSigXYZ_udsg_0_8JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_8JP_bAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_8JP_cAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_8JP_udsgAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
            //90%
         fHistManager.CreateTH2("fh2d_ImpSigXY_b_0_9JP", "fh2d_ImpSigXY_b_0_9JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_b_0_9JP", "fh2d_ImpSigXYZ_b_0_9JP;pt (GeV/c); sig ",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_c_0_9JP", "fh2d_ImpSigXY_c_0_9JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_c_0_9JP", "fh2d_ImpSigXYZ_c_0_9JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_udsg_0_9JP", "fh2d_ImpSigXY_udsg_0_9JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_udsg_0_9JP", "fh2d_ImpSigXYZ_udsg_0_9JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_9JP_bAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_9JP_cAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_9JP_udsgAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
            //95%
         fHistManager.CreateTH2("fh2d_ImpSigXY_b_0_95JP", "fh2d_ImpSigXY_b_0_95JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_b_0_95JP", "fh2d_ImpSigXYZ_b_0_95JP;pt (GeV/c); sig ",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_c_0_95JP", "fh2d_ImpSigXY_c_0_95JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_c_0_95JP", "fh2d_ImpSigXYZ_c_0_95JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_udsg_0_95JP", "fh2d_ImpSigXY_udsg_0_95JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXYZ_udsg_0_95JP", "fh2d_ImpSigXYZ_udsg_0_95JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_95JP_bAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_95JP_cAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_95JP_udsgAccepted","detector level jets;pt (GeV/c); count",500,0,250,"s");

         //IP distributions: all jet types combined
         fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_5JP", "fh2d_ImpSigXY_all_0_5JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_6JP", "fh2d_ImpSigXY_all_0_6JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_7JP", "fh2d_ImpSigXY_all_0_7JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_8JP", "fh2d_ImpSigXY_all_0_8JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_9JP", "fh2d_ImpSigXY_all_0_9JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
         fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_95JP", "fh2d_ImpSigXY_all_0_95JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");

         //Pt distributions: all jet types combined
         fHistManager.CreateTH1("fh1dJetRecPt_0_5JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_6JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_7JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_8JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_9JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
         fHistManager.CreateTH1("fh1dJetRecPt_0_95JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
       }//EndDoJetProb
*/
     //MC Impact Parameter Correction
    /* fHistManager.CreateTH2("fh2dTracksImpParXY_McCorr","radial imp. parameter (after correction);impact parameter xy (cm);a.u.",2000,-1,1,500,0,100,"s");
     fHistManager.CreateTH1("fh1dTracksImpParXY_McCorr","radial imp. parameter (after correction);impact parameter xy (cm);a.u.",400,-0.2,0.2,"s");
     fHistManager.CreateTH1("fh1dTracksImpParXYZ_McCorr","3d imp. parameter (after correction);impact parameter 3d (cm);a.u.",2000,0,100.,"s");
     fHistManager.CreateTH2("fh2dTracksImpParXYZ_McCorr","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-1,1,500,0,100.,"s");
     fHistManager.CreateTH2("fh2dTracksImpParXYZSignificance_McCorr","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.,"s");
     fHistManager.CreateTH1("fh1dTracksImpParXYSignificance_McCorr","radial imp. parameter (after correction);impact parameter xy significance;a.u.",2000,-30,30.,"s");
     fHistManager.CreateTH1("fh1dTracksImpParXYZSignificance_McCorr","3d imp. parameter (after correction);impact parameter 3d significance;a.u.",2000,0.,100.,"s");*/

     //MC General Information
     fHistManager.CreateTH1("fh1dJetGenPt","generator level jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetGenPtUnidentified","generator level jets (no flavour assigned);pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetGenPtudsg","generator level udsg jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetGenPtc","generator level c jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetGenPtb","generator level b jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetGenPts","generator level s jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH2("fh2dJetGenPtVsJetRecPt","detector momentum response;gen pt;rec pt",500,0,250,500,0,250,"s");
     fHistManager.CreateTH1("fh1dJetRecPtudsg","detector level jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetRecPtUnidentified","detector level jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetRecPtc","detector level jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetRecPtb","detector level jets;pt (GeV/c); count",250,0,250,"s"); 
     fHistManager.CreateTH1("fh1dJetRecPts","detector level jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetRecPtUnidentifiedAccepted","detector level jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetRecPtudsgAccepted","detector level jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetRecPtcAccepted","detector level jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetRecPtbAccepted","detector level jets;pt (GeV/c); count",250,0,250,"s");
     fHistManager.CreateTH1("fh1dJetRecPtsAccepted","detector level jets;pt (GeV/c); count",250,0,250,"s");
   }//EndPythiaLoop

    //Pt Distributions for N1,N2,N3 Tracks
    if(fIsPythia){
      fHistManager.CreateTH1("fh1dJetRecPt_n_1_b_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
      fHistManager.CreateTH1("fh1dJetRecPt_n_2_b_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
      fHistManager.CreateTH1("fh1dJetRecPt_n_3_b_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");

      fHistManager.CreateTH1("fh1dJetRecPt_n_1_c_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
      fHistManager.CreateTH1("fh1dJetRecPt_n_2_c_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
      fHistManager.CreateTH1("fh1dJetRecPt_n_3_c_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");

      fHistManager.CreateTH1("fh1dJetRecPt_n_1_udsg_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
      fHistManager.CreateTH1("fh1dJetRecPt_n_2_udsg_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
      fHistManager.CreateTH1("fh1dJetRecPt_n_3_udsg_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");

      fHistManager.CreateTH1("fh1dJetRecPt_n_1_s_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
      fHistManager.CreateTH1("fh1dJetRecPt_n_2_s_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
      fHistManager.CreateTH1("fh1dJetRecPt_n_3_s_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    }
    fHistManager.CreateTH1("fh1dJetRecPt_n_1_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPt_n_2_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPt_n_3_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dTrackPt_n_1_all_Accepted","detector level jets;pt (GeV/c); count",500,0,200,"s");
    fHistManager.CreateTH1("fh1dTrackPt_n_2_all_Accepted","detector level jets;pt (GeV/c); count",500,0,200,"s");
    fHistManager.CreateTH1("fh1dTrackPt_n_3_all_Accepted","detector level jets;pt (GeV/c); count",500,0,200,"s");

    //Template Generation
    const char * flavour[6]  = {"Unidentified","udsg","c","b","s",""};
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
    for (Int_t id = 0;id<2;++id)  // XY or XY/
      for (Int_t ifl = 0;ifl<6;++ifl)  //flavour
        for (Int_t io = 0;io<4;++io)        //order parameter
          for (Int_t is = 0;is<1;++is)          //special comment
            for (Int_t it = 0;it<2;++it){           //significance or not
              if(it==1) {
                iplow=-30;
                iphigh=30; //from 30
                if(io==0 && ifl==4) ipbins = 1000;//2000;
                  else  ipbins =1000;//2000;
              }else {
                iplow=-0.5;
                iphigh=0.5;
                ipbins =1000;//;2000;
              }
              if(id==0)  ipbins =1000;//2000;
                if((fIsPythia||(!fIsPythia && ifl==5))){
                  fHistManager.CreateTH2(Form("%s%s%s%s%s%s",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]),
                                Form("%s%s%s%s%s%s;;",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]),
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
}

void AliAnalysisTaskHFJetIPQA::PrintSettings(){
    TString jetcuts="";
    TString trackcuts="";
    TString vertexcuts="";
    Int_t version=1;

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
    printf("Cut Settings: %s\n",jetcuts.Data());

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
    trackcuts+=fAnalysisCuts[bAnalysisCut_MinTPCClus];

    printf("Cut Track Settings %s\n", trackcuts.Data());

    vertexcuts+=version;
    vertexcuts+="+";
    vertexcuts+=Form("%0.f",fAnalysisCuts[bAnalysisCut_NContibutors]);
    vertexcuts+="+";
    vertexcuts+=fAnalysisCuts[bAnalysisCut_SigmaDiamond];
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

    fh1dCuts->SetTitle(jetcuts.Data());
    fh1dCuts->GetXaxis()->SetTitle(trackcuts.Data());
    fh1dCuts->GetYaxis()->SetTitle(vertexcuts.Data());

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
        AliAODVertex *AliAnalysisTaskHFJetIPQA::RemoveDaughtersFromPrimaryVtx( const AliVTrack * const track) {
            const AliAODEvent * aod =  ((AliAODEvent*)InputEvent());
            AliAODVertex *vtxAOD =aod ->GetPrimaryVertex();
            if(!vtxAOD) return 0;
            TString title=vtxAOD->GetTitle();
            if(!title.Contains("VertexerTracks")) return 0;
            AliVertexerTracks vertexer(aod->GetMagneticField());
            vertexer.SetITSMode();
            vertexer.SetMinClusters(4);
            if(title.Contains("WithConstraint")) {
                Float_t diamondcovxy[3];
                aod->GetDiamondCovXY(diamondcovxy);
                Double_t pos[3]={aod->GetDiamondX(),aod->GetDiamondY(),0.};
                Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
                AliESDVertex diamond(pos,cov,1.,1);
                vertexer.SetVtxStart(&diamond);
            }

            Int_t skipped[5000]; for(Int_t i=0;i<5000;i++) skipped[i]=-1;
            Int_t id = (Int_t)track->GetID();
            if(!(id<0)) skipped[0] = id;
            int nTrksToSkip=1;

            Int_t nTracks=aod->GetNumberOfTracks();
            AliAODTrack * t = nullptr;
            AliExternalTrackParam etp_at_r39_old; etp_at_r39_old.CopyFromVTrack(track);
            etp_at_r39_old.PropagateTo(3.9,InputEvent()->GetMagneticField());
            double angle0 = TMath::ATan2(etp_at_r39_old.Yv(),etp_at_r39_old.Xv());
            double zz0    = etp_at_r39_old.GetZ();

            for(Int_t i=0; i<nTracks; i++){
                t = (AliAODTrack *)(aod->GetTrack(i));
                if(!((((AliAODTrack*)t)->TestFilterBit(4))))continue;
                id = (Int_t)t->GetID();
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
            vertexer.SetSkipTracks(nTrksToSkip,skipped);
            AliESDVertex *vtxESDNew = vertexer.FindPrimaryVertex(aod);
            if(!vtxESDNew) return 0;
            Int_t nContrib =vtxESDNew->GetNContributors();
            if(vtxESDNew->GetNContributors()<=2) {

                delete vtxESDNew; vtxESDNew=nullptr;
                return 0;
            }

            if(vtxESDNew->GetChi2toNDF()>fAnalysisCuts[bAnalysisCut_Z_Chi2perNDF]) {
                delete vtxESDNew; vtxESDNew=nullptr;
                return 0;
            }

            Double_t pos[3];
            Double_t cov[6];

            Double_t chi2perNDF;
    vtxESDNew->GetXYZ(pos); // position
    vtxESDNew->GetCovMatrix(cov); //covariance matrix
    chi2perNDF = vtxESDNew->GetChi2toNDF();
    if(vtxESDNew) delete vtxESDNew;
    vtxESDNew=NULL;
    AliAODVertex *vtxAODNew = new AliAODVertex(pos,cov,chi2perNDF);
    vtxAODNew->	SetNContributors(nContrib);
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

Bool_t AliAnalysisTaskHFJetIPQA::IsTrackAccepted(AliVTrack* track ,Int_t n){
    if(!track) return kFALSE;
    if(fIsEsd){
        fESDTrackCut->SetMinNClustersITS(abs(n));
        fESDTrackCut->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
        if(!(fESDTrackCut->AcceptTrack((AliESDtrack*)track))) return kFALSE;
        return kTRUE;
    }
    else {
            //HasMatchedGoodTracklet((AliAODTrack*)track);
        if(!(((AliAODTrack*)track)->TestFilterBit(9) || ((AliAODTrack*)track)->TestFilterBit(4)))return kFALSE;
        if(!(((AliAODTrack*)track->HasPointOnITSLayer(0))&&(AliAODTrack*)track->HasPointOnITSLayer(1)))  return kFALSE;
        if(((AliAODTrack*)track)->GetNcls(0)<abs(n)) return kFALSE;
        if(((AliAODTrack*)track)->GetNcls(1)<fAnalysisCuts[bAnalysisCut_MinTPCClus]) return kFALSE;
        if(track->Pt()<fAnalysisCuts[bAnalysisCut_MinTrackPt])return kFALSE;
        AliAODVertex *aodvertex = (( AliAODTrack *)track)->GetProdVertex();
        if(!aodvertex) return kFALSE;
        if(aodvertex->GetType()==AliAODVertex::kKink) return kFALSE;
        return kTRUE;
    }
    return kTRUE;
}

Bool_t AliAnalysisTaskHFJetIPQA::IsTrackAcceptedJP(AliVTrack* track ,Int_t n){
    // min pt = 0.5 instead of 1.
    if(!track) return kFALSE;
    if(fIsEsd){
        fESDTrackCut->SetMinNClustersITS(abs(n));
        fESDTrackCut->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
        if(!(fESDTrackCut->AcceptTrack((AliESDtrack*)track))) return kFALSE;
        return kTRUE;
    }
    else {
        if(!(((AliAODTrack*)track)->TestFilterBit(9) || ((AliAODTrack*)track)->TestFilterBit(4))){
            return kFALSE;
        }
        if(!(((AliAODTrack*)track->HasPointOnITSLayer(0))&&(AliAODTrack*)track->HasPointOnITSLayer(1)))  {
            return kFALSE;
        }
        if(((AliAODTrack*)track)->GetNcls(0)<abs(n)) return kFALSE;
        if(((AliAODTrack*)track)->GetNcls(1)<fAnalysisCuts[bAnalysisCut_MinTPCClus]) return kFALSE;
        if(track->Pt()<fAnalysisCuts[bAnalysisCut_MinTrackPtMC])return kFALSE;
        AliAODVertex *aodvertex = (( AliAODTrack *)track)->GetProdVertex();
        if(!aodvertex) return kFALSE;
        if(aodvertex->GetType()==AliAODVertex::kKink) return kFALSE;
        return kTRUE;
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
    Double_t matchingpar1 =0.25;
    Double_t matchingpar2 =0.25;
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
/*! \brief mysort
 *
 * custom strcut sorter function
*/
Bool_t AliAnalysisTaskHFJetIPQA::myquarksort(const SQuarks& i, const SQuarks& j)
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
/*! \brief IsJetTaggedJetProb
 *
 * unused
 */
    Bool_t AliAnalysisTaskHFJetIPQA::IsJetTaggedJetProb(Double_t thresProb)
    {
        return kTRUE;
    }

/*! \brief IsParton
 *
 * Return kTRUE if particle is parton
 */
Bool_t AliAnalysisTaskHFJetIPQA::IsParton(int pdg){
  return ((pdg==1)||(pdg==2)||(pdg==3)||(pdg==4)||(pdg==5)||(pdg==21));
}

/*! \brief StoreDaughters
 *
 * Store all daughters of MC particles
 */
void  AliAnalysisTaskHFJetIPQA::StoreDaughters(AliAODMCParticle* part, int kFirstMotherLabel){
    int kDaughLabel1=part->GetDaughterFirst();
    int kDaughLabel2=part->GetDaughterLast();

    for(int idaughter=kDaughLabel1;idaughter<=kDaughLabel2;idaughter++){
      if(idaughter==-1)continue;
      daughtermother.emplace(idaughter,kFirstMotherLabel);
      //printf("Pushing daughters to vector %i, Mother Label=%i, kFirstMother=%i\n",idaughter, part->GetMother(), kFirstMotherLabel);
    }
}

/*! \brief DoJetPartonMatching
 *
 * Performs matching of parton daughters to jet constituents. Returns the number of parton daughters that are matched with constituents.
 */
int AliAnalysisTaskHFJetIPQA::DoJetPartonMatching(const AliEmcalJet *jet, int partonlabel){
  int nFit=0;
  int kPartLabel=0;
  int kFirstMotherLabel=0;
  AliVParticle *vp=0x0;
  std::map<int ,int >::iterator it = daughtermother.begin();

  //loop over daughtermother
  while(it != daughtermother.end()){
    kPartLabel=it->first;
    kFirstMotherLabel=it->second;

    if(kFirstMotherLabel!=partonlabel){
        it++;
        continue;
    }
    //printf("Withing daughtermother Label %i, firstmotherlabel %i\n", kPartLabel, kFirstMotherLabel);

    //loop over jet constituents
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {//start trackloop jet
      vp = static_cast<AliVParticle*>(jet->Track(i));
      if (!vp){
        Printf("ERROR: AliVParticle associated to constituent not found\n");
        continue;
      }
      AliAODMCParticle * part = static_cast<AliAODMCParticle*>(vp);
      if(!part){
        printf("ERROR: Finding no Part!\n");
        return 0;
      }
      //printf("Jet Constitutent: %i\n", part->Label());
      if(kPartLabel==part->Label()){
          //printf("Found matching for label %i\n",part->Label());
          nFit++;
      }
    }//end trackloop jet
    it++;
  }
  //printf("Number of fits: %i\n",nFit);
  return nFit;
}

/*! \brief DoFlavourVectorFill
 *
 * Fills parton vectors. If fDoFlavourMatching=kTRUE, only partons that have been matched to jet constituents are filled. If kFALSE, alls partons are accepted.
 */
void AliAnalysisTaskHFJetIPQA::DoFlavourVectorFill(int kJetOrigin, double partonpt, double jetpt, int kPartonsInJet, double deta){
    if(abs(kJetOrigin) == 5) {
      if(kPartonsInJet==0){
          fh2dBottomNotContrib->Fill(partonpt,jetpt);
          //printf("Entry in fh2dBottomNotContrib\n");
      }
      if(kPartonsInJet>0||!fDoFlavourMatching){
          fh2dBottomNMatch->Fill(kPartonsInJet,partonpt);
          fh2dBottomDeta->Fill(partonpt,deta);
          fPBJet.push_back(partonpt);
          //printf("Pushing back bottom\n");
      }
    }
    else if(abs(kJetOrigin)== 4) {
      if(kPartonsInJet==0){
          fh2dCharmNotContrib->Fill(partonpt,jetpt);
          //printf("Entry in fh2dCharmNotContrib\n");
      }
      if(kPartonsInJet>0||!fDoFlavourMatching){
          fh2dCharmNMatch->Fill(kPartonsInJet,partonpt);
          fh2dCharmDeta->Fill(partonpt,deta);
          fPCJet.push_back(partonpt);
          //printf("Pushing back charm\n");
      }
    }
    else if(abs(kJetOrigin) == 3 ) {
      if(kPartonsInJet==0){
          fh2dLightNotContrib->Fill(partonpt,jetpt);
          //printf("Entry in fh2dLightNotContrib\n");
      }
      if(kPartonsInJet>0||!fDoFlavourMatching){
          fh2dLightNMatch->Fill(kPartonsInJet,partonpt);
          fh2dLightDeta->Fill(partonpt,deta);
          fPSJet.push_back(partonpt);
          //printf("Strange pushed with p=%f",partonpt);
      }
    }
    else if(abs(kJetOrigin)== 1 ||abs(kJetOrigin)== 2 ||  abs(kJetOrigin) == 21) {
      if(kPartonsInJet==0){
          fh2dLightNotContrib->Fill(partonpt,jetpt);
          //printf("Entry in fh2dLightNotContrib\n");
      }
      if(kPartonsInJet>0||!fDoFlavourMatching){
          fh2dLightNMatch->Fill(kPartonsInJet,partonpt);
          fh2dLightDeta->Fill(partonpt,deta);
          fPUdsgJet.push_back(partonpt);
          //printf("Light pushed with p=%f",partonpt);
      }
    }
}
/*! \brief DoFlavourDecision
 *
 * Perform final decision on flavor assigned to jet depending on particle flavour ordering (bottom accepted
 *  above charm accepted above light flavours) and pt ordering (of light flavours, the hardest is accepted)
 */
int AliAnalysisTaskHFJetIPQA::DoFlavourDecision(){
      double p_udsg_max=-999;  // maximal pt of light quarks
      double p_s_max=-999;      // maximal pt of strange quarks
      bool is_udg=kFALSE;

      if(fPCJet.size() ==0&& fPBJet.size()==0&& fPSJet.size()==0&&fPUdsgJet.size()==0) return 0; //udsg
      //check for c jet
      for (Int_t icj = 0 ; icj <(Int_t)fPCJet.size();++icj ){
          //printf("Charm Flavour Jet!\n");
          return 2;
      }
      //check for b jet
      for (Int_t icj = 0 ; icj <(Int_t)fPBJet.size();++icj ){
          //printf("Bottom Flavour Jet!\n");
          return 3;
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

      //take hardest
      if(p_s_max>p_udsg_max){
        //printf("S prefered with psmax=%f, pudsgmax=%f\n", p_s_max,p_udsg_max);
        return 4;
      }
      else{
          if(fPUdsgJet.size()!=0){
              //printf("Light prefered with psmax=%f, pudsgmax=%f\n", p_s_max,p_udsg_max);
              is_udg =kTRUE;
              return 1;
          }
      }
      return 0;
}

/*! \brief IsMCJetPartonFast
 *
 * Fast jet parton MC matcher
 */
Int_t  AliAnalysisTaskHFJetIPQA::IsMCJetPartonFast(const AliEmcalJet *jet, Double_t radius,Bool_t &is_udg){
  fJetCont.clear();
  fPUdsgJet.clear();
  fPSJet.clear();
  fPCJet.clear();
  fPBJet.clear();
  daughtermother.clear();
  fHardProcess.clear();
  fQuarkVec.clear();


  int kPartonsInJet=-1;      //number of partons contributed to the jet by one flavour
  int kPartonsPerJet=0;     //how many flavors matched with jet
  int kJetOrigin=-999;      //label of flavor that is matched with jet
  int kStoreRun=0;          //if find hardest process go back 2 steps in the particle listing and store daughters of hardest process particles

  if(!jet) return 0;
  if(!(jet->GetNumberOfTracks()>fNoJetConstituents)){
    return 0;
  }
  //************************************
  //loop over all partons within event
  for(Int_t iPrim = 0 ; iPrim<fMCArray->GetEntriesFast();iPrim++){//start trackloop filling fQuarkVec
    AliAODMCParticle * part = static_cast<AliAODMCParticle*>(fMCArray->At(iPrim));
    if(!part) return 0;

    //initialise variables
    Double_t eta = part->Eta();
    Double_t phi= part->Phi();
    Double_t etajet = jet->Eta();
    Double_t phijet = jet->Phi();
    Double_t deta = etajet - eta;
    Double_t dphi = phijet-phi;
    Double_t pt=part->Pt();
    Int_t pdg = (abs(part->GetPdgCode()));
    Int_t kPartLabel=iPrim;
    dphi = TVector2::Phi_mpi_pi(dphi);
    Double_t  d = sqrt(deta * deta + dphi * dphi);
    Int_t kMotherLabel=part->GetMother();

    SQuarks q(pt,  kPartLabel,pdg,d);

    //if(kStoreRun>0)printf("Doing kStorerun %i, iPrim=%i\n",kStoreRun,iPrim);

    if(fPythia6){
      if(!IsParton(pdg)) continue;
      //printf("Going into Pythia 6 loop\n");
      if(!(part->MCStatusCode()==11)&&(part->MCStatusCode()==12)){
        continue;
      }
      else{
        fQuarkVec.push_back(q);
      }
    }
    if(fPythia8){
      //printf("Going into Pythia 8 loop\n");
      for(int iquark=0;iquark<fHardProcess.size();iquark++){
        if(!IsParton(pdg)) continue;
        if((pt==fHardProcess[iquark].first)&&(TMath::Abs(pt)>0.0000000001)&&(kMotherLabel!=fHardProcess[iquark].second)&&(kStoreRun==0)){
          //printf("Found pairing of quarks! q1_label=%i, q2_label=%i, q1_pdg=%i, q2_pdg=%i, q1pt=%.10f, q2pt=%.10f\n", iPrim, fHardProcess[iquark].second, pdg, fHardProcess[iquark].HasPDG, pt, fHardProcess[iquark].first);

          fQuarkVec.push_back(q);
          fQuarkVec.push_back(fHardProcess[iquark]);

          daughtermother.emplace(q.second, q.second);
          daughtermother.emplace(fHardProcess[iquark].second,fHardProcess[iquark].second);
          //printf("Pushing: q1_label=%i, q2_label=%i\n", q.second, fHardProcess[iquark].second);
          kStoreRun=3;
        }
      }

      if(kStoreRun==0)fHardProcess.push_back(q);
      if(kStoreRun!=3){  //Only store daughters if not only just now hardest process particles have been found. Otherwise double counting.
          int kMotherLabel=-999;
          std::map<int,int>::iterator it;
          it = daughtermother.find(kPartLabel);
          if(it != daughtermother.end()){
            kMotherLabel=it->second;
            //printf("Found daughter with %i, pdg =%i, kMotherLabel=%i\n",iPrim,pdg,kMotherLabel);
          }
          if(kMotherLabel!=-999)StoreDaughters(part,kMotherLabel);
      }
      if(kStoreRun>0){  //if found hardest process particles go back 2 steps in eventlisting to store daughter particles
        kStoreRun=kStoreRun-1;
        if(kStoreRun==2)iPrim=iPrim-2;
      }
    }
  }//end trackloop filling fQuarkVec

  //printf("Inside QuarkVec:\n");
  //for(int i=0;i<fQuarkVec.size();i++){
  //  printf("pt =%f, label =%i, pdg=%i, deta=%f\n", fQuarkVec[i].first, fQuarkVec[i].second, fQuarkVec[i].HasPDG, fQuarkVec[i].HasDETA);
  //}

 //**************************************
 //loop over jet constituents
 for(int iquark=0;iquark<fQuarkVec.size();iquark++){
    kJetOrigin=-999;
    kPartonsInJet=-1;
    if(fQuarkVec[iquark].HasDETA<fDaughtersRadius){
        kJetOrigin=fQuarkVec[iquark].HasPDG;
        //printf("Found quark with right radius: Label=%i, deta=%f, pt=%f, pdg=%i\n",fQuarkVec[iquark].second, fQuarkVec[iquark].HasDETA,fQuarkVec[iquark].first, fQuarkVec[iquark].HasPDG);

        if(fDoFlavourMatching){
          kPartonsInJet=DoJetPartonMatching(jet, fQuarkVec[iquark].second);
          if(kPartonsInJet==0){
              //printf("No Jet matched!\n");
          }
          else{
              kPartonsPerJet=kPartonsPerJet+1;
          }
        }//fJetLoop for DoFlavourMatching
    }
    DoFlavourVectorFill(kJetOrigin, fQuarkVec[iquark].first,jet->Pt(),kPartonsInJet, fQuarkVec[iquark].HasDETA);
  }//end trackloop MC

  if(kPartonsPerJet>1){
      //printf("More than one, exactly %i partons matched to jet\n",kPartonsPerJet);
      fh2dManifoldParton->Fill(kPartonsPerJet,jet->Pt());
  }

  return DoFlavourDecision();
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
    if(sName.EqualTo("fh1dCuts")){
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
    const Double_t kBeampiperadius=5;
    fEventVertex = RemoveDaughtersFromPrimaryVtx(track);
    if(!fEventVertex)return kFALSE;
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

AliExternalTrackParam AliAnalysisTaskHFJetIPQA::GetExternalParamFromJet(const AliEmcalJet *jet, const AliAODEvent *event)
{
    double vtx[3]= {0.};
    double cov [21] = {0.};
    double pxpypz[3] = {0.};
    jet->PxPyPz(pxpypz);
    (event->GetPrimaryVertex())->GetXYZ(vtx);
    AliExternalTrackParam etp     (vtx, pxpypz, cov, (Short_t)0);
    return etp;
}



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


Bool_t AliAnalysisTaskHFJetIPQA::GetImpactParameterWrtToJet(const AliAODTrack *track, const AliAODEvent *event, const AliEmcalJet *jet, Double_t *dca, Double_t *cov, Double_t *XYZatDCA, Double_t &jetsign)
{
    if(!track || !event || !jet)return kFALSE;
    if(dca==0 || cov ==0 ||XYZatDCA ==0 ) return kFALSE;
    if(!GetImpactParameter(track,event,dca,cov,XYZatDCA)) return kFALSE;
    Double_t VxVyVz[3]= {0.,0.,0.};
    const  AliVVertex *vtxESDSkip =fEventVertex;//GetKFPrimaryVertex();//;RemoveDaughtersFromPrimaryVtx(track);//THIS IS THE WORKING RECALC OPTION
    if(!fEventVertex) return kFALSE;
    vtxESDSkip->GetXYZ(VxVyVz);
    double jetp[3];
    jet->PxPyPz(jetp);
    TVector3 jetP3(jetp);

    TVector3 JetDir =jetP3.Unit();
    TVector3 D0(XYZatDCA);
    TVector3 vertex(VxVyVz);
    TVector3 DD0(D0.x()-vertex.x(),D0.y()-vertex.y(),0.);
    double ps =DD0.Dot(JetDir);
    double value = DD0.Mag()*(ps/fabs(ps));
    jetsign  = TMath::Sign(1.,value);
    TVector3 dd0 =DD0.Unit();
    AliExternalTrackParam etp;
    etp.CopyFromVTrack((AliVTrack*)track);
    double c[3],d[2];
    if(!etp.PropagateToDCA(vtxESDSkip, InputEvent()->GetMagneticField(), 10., d, c))return kFALSE;
    double cm[21];
    etp.GetCovarianceXYZPxPyPz(cm);
    double cve[6];
    vtxESDSkip->GetCovarianceMatrix(cve);
    TMatrix vertexError(3,3);
    vertexError(0,0) =cve[0];
    vertexError(0,1) =cve[1];
    vertexError(1,0) =cve[1];
    vertexError(1,1) =cve[2];
    vertexError(0,2) =cve[3];
    vertexError(2,0) =cve[3];

    vertexError(1,2) =cve[4];
    vertexError(1,2) =cve[4];
    vertexError(2,2) =cve[5];

    TMatrix trackError(6,6);
    trackError(0,0) =cm[0];
    trackError(1,0) =cm[1];
    trackError(2,0) =cm[3];
    trackError(3,0) =cm[6];
    trackError(4,0) =cm[10];
    trackError(5,0) =cm[15];

    trackError(0,1) =cm[1];
    trackError(1,1) =cm[2];
    trackError(2,1) =cm[4];
    trackError(3,1) =cm[7];
    trackError(4,1) =cm[11];
    trackError(5,1) =cm[16];

    trackError(0,2) =cm[3];
    trackError(1,2) =cm[4];
    trackError(2,2) =cm[5];
    trackError(3,2) =cm[8];
    trackError(4,2) =cm[12];
    trackError(5,2) =cm[17];

    trackError(0,3) =cm[6];
    trackError(1,3) =cm[7];
    trackError(2,3) =cm[8];
    trackError(3,3) =cm[9];
    trackError(4,3) =cm[13];
    trackError(5,3) =cm[18];

    trackError(0,4) =cm[10];
    trackError(1,4) =cm[11];
    trackError(2,4) =cm[12];
    trackError(3,4) =cm[13];
    trackError(4,4) =cm[14];
    trackError(5,4) =cm[19];

    trackError(0,5) =cm[15];
    trackError(1,5) =cm[16];
    trackError(2,5) =cm[17];
    trackError(3,5) =cm[18];
    trackError(4,5) =cm[19];
    trackError(5,5) =cm[20];
    ROOT::Math::SMatrix<double,3> E33Vertex;
    ROOT::Math::SMatrix<double,6> E66Track;

    for (int i = 0;i<3;++i) {
        for (int j = 0;j<3;++j) {
            E33Vertex (i,j) = vertexError(i,j);
        }
    }
    for (int i = 0;i<6;++i) {
        for (int j = 0;j<6;++j) {
            E66Track (i,j) = trackError(i,j);
        }
    }

    ROOT::Math::SVector<double,6>  deriv;
    ROOT::Math::SVector<double,3>  deriv_v;

    deriv_v[0] = -dd0.x();
    deriv_v[1] = -dd0.y();
    deriv_v[2] = -dd0.z();

    deriv[0] = dd0.x();
    deriv[1] = dd0.y();
    deriv[2] = dd0.z();
    deriv[3] = 0.;
    deriv[4] = 0.;
    deriv[5] = 0.;

  
    AliExternalTrackParam etp_jet = GetExternalParamFromJet(jet,event);
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

    double dcajetrack = sqrt((xyz_jet_global[0] - xyz_track_global[0]) * (xyz_jet_global[0] - xyz_track_global[0]) +
        (xyz_jet_global[1] - xyz_track_global[1]) * (xyz_jet_global[1] - xyz_track_global[1])+
        (xyz_jet_global[2] - xyz_track_global[2]) * (xyz_jet_global[2]- xyz_track_global[2]));

    if(bdecaylength>bAnalysisCut_MaxDecayLength) return kFALSE;
    if(dcajetrack  >bAnalysisCut_DCAJetTrack) return kFALSE;
    return kTRUE;
}

Double_t AliAnalysisTaskHFJetIPQA::CalculateJetProb(AliEmcalJet *jet)
{
    if(!jet) return -1.;
    //Loop over all tracks calculate P(s) for all accepted later add looser cuts also
    Int_t ntracks = (Int_t)jet->GetNumberOfTracks();
    AliJetContainer * jetconrec = 0x0;
    jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    // Int_t jcounter =0;
    double SumJet =0;
    double ProbJet = 0;
    double Loginvlog=0;
    int ngoodtracks=0;

    double TrackProb=0;


    std::vector <double> probabilities;

    for(Int_t itrack = 0; itrack < ntracks; ++itrack)
    {
        TrackProb=0;
        AliVTrack * trackV = (((AliVTrack*)((jetconrec->GetParticleContainer())->GetParticle(jet->TrackAt(itrack)))));
        if(!trackV) continue;
        Bool_t isAccepted=kFALSE;
        Int_t track_class=0;
        if(IsTrackAcceptedJP(trackV,6))     {isAccepted=kTRUE;track_class=0;}
        else if(IsTrackAcceptedJP(trackV,5)) {isAccepted=kTRUE;track_class=1;}
        else if(IsTrackAcceptedJP(trackV,4)) {isAccepted=kTRUE;track_class=2;}
        else if(IsTrackAcceptedJP(trackV,3)) {isAccepted=kTRUE;track_class=3;}

        int PIDClass = -1;
        if(IsFromElectron((AliAODTrack*)trackV)) PIDClass=0;
        else if(IsFromPion((AliAODTrack*)trackV)) PIDClass=1;
        else if(IsFromKaon((AliAODTrack*)trackV)) PIDClass=2;
        else if(IsFromProton((AliAODTrack*)trackV)) PIDClass=3;
        if(isAccepted){
            Double_t dca[2] = {0.,0.};
            Double_t cov[3] = {0.,0.,0.};
            Double_t xyz[3] = {0.,0.,0.};
            Double_t sign = 1.;
            if(!GetImpactParameterWrtToJet((AliAODTrack*)trackV,(const AliAODEvent*)InputEvent(),jet,dca,cov,xyz,sign)) continue;
            double tmpProb;
            if(!fUsePIDJetProb || PIDClass <0)  tmpProb =CalculatePSTrack(sign,GetValImpactParameter(kXYSig,dca,cov),trackV->Pt() ,track_class);
            else tmpProb =CalculatePSTrackPID(sign,GetValImpactParameter(kXYSig,dca,cov),trackV->Pt() ,track_class,PIDClass);

            if       (tmpProb<0. && sign <0.)TrackProb=-tmpProb;
            else if  (tmpProb>0. && sign >0.)TrackProb=tmpProb;
            else if  (tmpProb>0. && sign ==0.)TrackProb=0.5*tmpProb;
            else if  (tmpProb<0. && sign ==0.)TrackProb=1.+0.5*tmpProb;
            else continue;
            probabilities.push_back(TrackProb);
            ngoodtracks++;
        }
    }
    std::sort (probabilities.begin(), probabilities.end());
    std::reverse(probabilities.begin(),probabilities.end());
    double m_minTrackProb = 5e-2;
    for(std::vector<double>::const_iterator q = probabilities.begin(); q != probabilities.end(); q++){
        SumJet+=(*q>m_minTrackProb)?log(*q):log(m_minTrackProb);
    }
    if(SumJet<0.){
        if(ngoodtracks>=2){
            Loginvlog=log(-SumJet);
        }
        double Prob=1.;
        double lfact=1.;
        for(int l=1; l!=ngoodtracks; l++){
            lfact*=l;
            Prob+=exp(l*Loginvlog-log(1.*lfact));
        }
        double LogProb=log(Prob);
        ProbJet=
        std::min(exp(std::max(LogProb+SumJet,-30.)),1.);
    }else{
        ProbJet=1.;
    }
    return 1.-ProbJet;
}

Double_t AliAnalysisTaskHFJetIPQA::CalculatePSTrack(Double_t sign, Double_t significance ,Double_t trackPt,Int_t trclass)
{
    Double_t retval = 0;
    //switch resolution function based on track pt;
    int ptbin = 0;
    if(trackPt >0.5 && trackPt<1.)    ptbin=0;
    else if(trackPt >1. && trackPt<2.)ptbin=1;
    else if(trackPt >2. && trackPt<4.)ptbin=2;
    else if(trackPt >4. && trackPt<6.)ptbin=3;
    else if(trackPt >6.)              ptbin=4;
    if(TMath::Abs(significance) >99) significance =99; //Limit to function definition range
    retval = sign * ((fResolutionFunction[20*4 + 4*trclass +ptbin])).Eval(TMath::Abs(significance));
    return retval;
}

Double_t AliAnalysisTaskHFJetIPQA::CalculatePSTrackPID(Double_t sign, Double_t significance ,Double_t trackPt,Int_t trclass,Int_t species)
{

    Double_t retval = 0;
    //switch resolution function based on track pt;
    int ptbin = 0;
    if(trackPt >0.5 && trackPt<1.)ptbin=0;
    else if(trackPt >1. && trackPt<2.)ptbin=1;
    else if(trackPt >2. && trackPt<4.)ptbin=2;
    else if(trackPt >4. && trackPt<6.)ptbin=3;
    else if(trackPt >6.)ptbin=4;
    if(TMath::Abs(significance) >99) significance =99; //Limit to function definition range
    retval = sign * ((fResolutionFunction[20*species + 4*trclass +ptbin])).Eval(TMath::Abs(significance));
    return retval;
}


void AliAnalysisTaskHFJetIPQA::Terminate(Option_t *){

    printf("\n*********************************\n");
    printf("Corrections:\n");
    printf("    MC Corrections (Data/MC+Fluka):%i\n",fDoMCCorrection);
    printf("    Track Smearing:%i\n",fRunSmearing );
    printf("    Underlying Event Subtraction:%i\n", fDoUnderlyingEventSub);
    printf("*********************************\n");
}
