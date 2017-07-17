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
#include "AliRDHFJetsCuts.h"
#include "AliGenEventHeader.h"
#include "AliVertexerTracks.h"
#include "AliEmcalList.h"
#include "AliHFJetsTagging.h"
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
using std::min;
using std::cout;
using std::endl;
using std::vector;
using std::pair;
ClassImp(AliAnalysisTaskHFJetIPQA)

AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA():
    AliAnalysisTaskEmcalJet(),fHistManager(),
    fEventVertex(nullptr),
    fPidResponse(nullptr),
    fUsePIDJetProb(kFALSE),
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
    fMCArray(nullptr),
    fJetCutsHF(new AliRDHFJetsCuts()),
    fMCEvent(nullptr),
    fESDTrackCut(nullptr),
    fVertexer(nullptr),
    fMcEvtSampled(kFALSE),
    fBackgroundFactorLinus{0},
    fEtaSEvt(100),fPhiSEvt(100),fEtaBEvt(100),fPhiBEvt(100),fEtaCEvt(100),fPhiCEvt(100),fEtaUdsgEvt(100),fPhiUdsgEvt(100),
    fAnalysisCuts{0,0,0,0,0,0,0,0,0,0,0},
    fCombined(nullptr),
    fXsectionWeightingFactor(1),
    fProductionNumberPtHard(-1),
    fMCglobalDCAxyShift(0.007),
    fVertexRecalcMinPt(1.0)
{
    SetMakeGeneralHistograms(kTRUE);
    SetDefaultAnalysisCuts();
    for(Int_t i =0 ; i<498;++i)for(Int_t j =0 ; j<19;++j)  fBackgroundFactorLinus[j][i]=1.;
    SetNeedEmcalGeom(kFALSE);
    SetOffTrigger(AliVEvent::kMB);
    SetUseAliAnaUtils(kTRUE,kTRUE);
    SetVzRange(-10,10);
    SetUseSPDTrackletVsClusterBG(kTRUE);
    for(Int_t i =0 ; i<200;++i)this->fResolutionFunction[i].Set(1000);


    DefineOutput(1,  TList::Class()) ;
}
AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA(const char *name):
    AliAnalysisTaskEmcalJet(name, kTRUE),fHistManager(name),
    fEventVertex(nullptr),
    fPidResponse(nullptr),
    fUsePIDJetProb(kFALSE),
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
    fMCArray(nullptr),
    fJetCutsHF(new AliRDHFJetsCuts()),
    fMCEvent(nullptr),
    fESDTrackCut(nullptr),
    fVertexer(nullptr),
    fMcEvtSampled(kFALSE),
    fBackgroundFactorLinus{0},
    fEtaSEvt(100),fPhiSEvt(100),fEtaBEvt(100),fPhiBEvt(100),fEtaCEvt(100),fPhiCEvt(100),fEtaUdsgEvt(100),fPhiUdsgEvt(100),
    fAnalysisCuts{0,0,0,0,0,0,0,0,0,0,0},
    fCombined(nullptr),
    fXsectionWeightingFactor(1.),
    fProductionNumberPtHard(-1),
    fMCglobalDCAxyShift(0.007),
    fVertexRecalcMinPt(1.0)

{
    SetNeedEmcalGeom(kFALSE);
    SetOffTrigger(AliVEvent::kMB);
    SetUseAliAnaUtils(kTRUE,kTRUE);
    SetVzRange(-10,10);
    SetUseSPDTrackletVsClusterBG(kTRUE);
    SetMakeGeneralHistograms(kTRUE);
    SetDefaultAnalysisCuts();
    DefineOutput(1,  TList::Class()) ;
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
}
/*! \brief setN_ITSClusters_Input_global
 *
 *
 * Set default number of its clusters cuts
 */

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
    //2d
    FillHist("fh2dTracksImpParXY",GetValImpactParameter(kXY,dca,cov),track->Pt(),1.*this->fXsectionWeightingFactor );
    FillHist("fh2dTracksImpParXYZ",GetValImpactParameter(kXYZ,dca,cov),track->Pt(),1.*this->fXsectionWeightingFactor );
    FillHist("fh2dTracksImpParZ",dca[1],track->Pt(),1.*this->fXsectionWeightingFactor );
    FillHist("fh2dTracksImpParXYSignificance",GetValImpactParameter(kXYSig,dca,cov),track->Pt(),1.*this->fXsectionWeightingFactor );
    FillHist("fh2dTracksImpParXYZSignificance",GetValImpactParameter(kXYZSig,dca,cov),track->Pt(),1.*this->fXsectionWeightingFactor );
    FillHist("fh2dTracksImpParZSignificance",GetValImpactParameter(kZSig,dca,cov),track->Pt(),1.*this->fXsectionWeightingFactor );
    //1d
    FillHist("fh1dTracksImpParXY",GetValImpactParameter(kXY,dca,cov),1.*this->fXsectionWeightingFactor );
    FillHist("fh1dTracksImpParXYZ",GetValImpactParameter(kXYZ,dca,cov),1.*this->fXsectionWeightingFactor );
    FillHist("fh1dTracksImpParXYSignificance",GetValImpactParameter(kXYSig,dca,cov),1.*this->fXsectionWeightingFactor );
    FillHist("fh1dTracksImpParXYZSignificance",GetValImpactParameter(kXYZSig,dca,cov),1*this->fXsectionWeightingFactor );
    //mc
    if(fIsPythia){
            FillHist("fh1dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),weight*this->fXsectionWeightingFactor );
            FillHist("fh1dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),weight*this->fXsectionWeightingFactor );
            FillHist("fh1dTracksImpParXYSignificance_McCorr",GetValImpactParameter(kXYSig,dca,cov),weight*this->fXsectionWeightingFactor );
            FillHist("fh1dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),weight*this->fXsectionWeightingFactor );
            FillHist("fh2dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),track->Pt(),weight*this->fXsectionWeightingFactor );
            FillHist("fh2dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),track->Pt(),weight*this->fXsectionWeightingFactor );
            FillHist("fh2dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),track->Pt(),weight*this->fXsectionWeightingFactor );
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
void AliAnalysisTaskHFJetIPQA::EventwiseCleanup(){
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



Bool_t AliAnalysisTaskHFJetIPQA::Run(){
    FillGeneralHistograms();
    /*Vertex Pos Selection*/
    fEventVertex = dynamic_cast<const AliAODVertex*>(InputEvent()->GetPrimaryVertex());

    if(!fEventVertex) {
            return kFALSE;
        }
    if(fEventVertex->GetNContributors()<1) {
            return kFALSE;
        }
    if(fEventVertex->GetChi2perNDF()>3.5*3.5) {
            return kFALSE;
        }
    if(fEventVertex->GetNContributors()<(int)(fAnalysisCuts[bAnalysisCut_NContibutors])) {
            return kFALSE;
        }
    if(TMath::Abs(fEventVertex->GetZ())>=fAnalysisCuts[bAnalysisCut_MaxVtxZ]) {
            return kFALSE;
        }


    this->fXsectionWeightingFactor=1;
    if(fIsPythia){
            if(fProductionNumberPtHard>0){
                    Int_t pTHard=fPtHard;

                    if(fProductionNumberPtHard==1){
                            //LHC15g6D
                            if((pTHard >= 5) && (pTHard < 7) )
                                fXsectionWeightingFactor = 1.389571e+00;
                            else if((pTHard >= 7) && (pTHard < 9) )
                                fXsectionWeightingFactor = 1.246220e+00;
                            else if((pTHard >= 9) && (pTHard < 12) )
                                fXsectionWeightingFactor = 1.134698e+00;
                            else if((pTHard >= 12) && (pTHard < 16) )
                                fXsectionWeightingFactor = 6.498058e-01;
                            else if((pTHard >= 16) && (pTHard < 21) )
                                fXsectionWeightingFactor = 2.806916e-01;
                            else if((pTHard >= 21) && (pTHard < 28) )
                                fXsectionWeightingFactor = 1.177605e-01;
                            else if((pTHard >= 28) && (pTHard < 36) )
                                fXsectionWeightingFactor = 3.852459e-02;
                            else if((pTHard >= 36) && (pTHard < 45) )
                                fXsectionWeightingFactor = 1.381213e-02;
                            else if((pTHard >= 45) && (pTHard < 57) )
                                fXsectionWeightingFactor = 5.940404e-03;
                            else if((pTHard >= 57) && (pTHard < 70) )
                                fXsectionWeightingFactor = 2.087004e-03;
                            else if((pTHard >= 70) && (pTHard < 85) )
                                fXsectionWeightingFactor = 8.527722e-04;
                            else if((pTHard >= 85) && (pTHard < 99) )
                                fXsectionWeightingFactor = 3.151402e-04;
                            else if((pTHard >= 99) && (pTHard < 115) )
                                fXsectionWeightingFactor = 1.595649e-04;
                            else if((pTHard >= 115) && (pTHard < 132) )
                                fXsectionWeightingFactor = 7.689924e-05;
                            else if((pTHard >= 132) && (pTHard < 150) )
                                fXsectionWeightingFactor = 3.878501e-05;
                            else if((pTHard >= 150) && (pTHard < 169) )
                                fXsectionWeightingFactor =2.032300e-05 ;
                            else if((pTHard >= 169) && (pTHard < 190) )
                                fXsectionWeightingFactor = 1.137054e-05;
                            else if((pTHard >= 190) && (pTHard < 212) )
                                fXsectionWeightingFactor = 6.136061e-06;
                            else if((pTHard >= 212) && (pTHard < 235) )
                                fXsectionWeightingFactor = 3.411389e-06;
                            else if((pTHard >= 235) && (pTHard < 1000000) )
                                fXsectionWeightingFactor = 4.847459e-06;
                            else
                                fXsectionWeightingFactor = 0;
                        }else  if(fProductionNumberPtHard==2){
                            //LHC15g6C
                            if((pTHard >= 5) && (pTHard < 7) )
                                fXsectionWeightingFactor = 1.052070e+00;
                            else if((pTHard >= 7) && (pTHard < 9) )
                                fXsectionWeightingFactor = 1.046607e+00;
                            else if((pTHard >= 9) && (pTHard < 12) )
                                fXsectionWeightingFactor = 1.023076e+00;
                            else if((pTHard >= 12) && (pTHard < 16) )
                                fXsectionWeightingFactor = 5.988252e-01;
                            else if((pTHard >= 16) && (pTHard < 21) )
                                fXsectionWeightingFactor = 2.616957e-01;
                            else if((pTHard >= 21) && (pTHard < 28) )
                                fXsectionWeightingFactor = 1.095912e-01;
                            else if((pTHard >= 28) && (pTHard < 36) )
                                fXsectionWeightingFactor =3.551833e-02;
                            else if((pTHard >= 36) && (pTHard < 45) )
                                fXsectionWeightingFactor = 1.280720e-02;
                            else if((pTHard >= 45) && (pTHard < 57) )
                                fXsectionWeightingFactor = 5.475954e-03;
                            else if((pTHard >= 57) && (pTHard < 70) )
                                fXsectionWeightingFactor = 1.926764e-03;
                            else if((pTHard >= 70) && (pTHard < 85) )
                                fXsectionWeightingFactor = 7.886817e-04;
                            else if((pTHard >= 85) && (pTHard < 99) )
                                fXsectionWeightingFactor = 2.908735e-04;
                            else if((pTHard >= 99) && (pTHard < 115) )
                                fXsectionWeightingFactor = 1.457728e-04;
                            else if((pTHard >= 115) && (pTHard < 132) )
                                fXsectionWeightingFactor = 7.072972e-05;
                            else if((pTHard >= 132) && (pTHard < 150) )
                                fXsectionWeightingFactor = 3.543739e-05;
                            else if((pTHard >= 150) && (pTHard < 169) )
                                fXsectionWeightingFactor = 1.859674e-05;
                            else if((pTHard >= 169) && (pTHard < 190) )
                                fXsectionWeightingFactor = 1.042401e-05;
                            else if((pTHard >= 190) && (pTHard < 212) )
                                fXsectionWeightingFactor = 5.609724e-06;
                            else if((pTHard >= 212) && (pTHard < 235) )
                                fXsectionWeightingFactor = 3.103228e-06;
                            else if((pTHard >= 235) && (pTHard < 1000000) )
                                fXsectionWeightingFactor = 4.408599e-06;
                            else
                                fXsectionWeightingFactor = 0;
                        }else  if(fProductionNumberPtHard==3){
                            //LHC15G6E
                            if((pTHard >= 5) && (pTHard < 7) )
                                fXsectionWeightingFactor = 1.389406e+00;
                            else if((pTHard >= 7) && (pTHard < 9) )
                                fXsectionWeightingFactor = 1.246014e+00;
                            else if((pTHard >= 9) && (pTHard < 12) )
                                fXsectionWeightingFactor = 1.135135e+00;
                            else if((pTHard >= 12) && (pTHard < 16) )
                                fXsectionWeightingFactor = 6.494695e-01;
                            else if((pTHard >= 16) && (pTHard < 21) )
                                fXsectionWeightingFactor = 2.801671e-01;
                            else if((pTHard >= 21) && (pTHard < 28) )
                                fXsectionWeightingFactor = 1.177731e-01 ;
                            else if((pTHard >= 28) && (pTHard < 36) )
                                fXsectionWeightingFactor = 3.851400e-02;
                            else if((pTHard >= 36) && (pTHard < 45) )
                                fXsectionWeightingFactor = 1.379875e-02;
                            else if((pTHard >= 45) && (pTHard < 57) )
                                fXsectionWeightingFactor = 5.941354e-03;
                            else if((pTHard >= 57) && (pTHard < 70) )
                                fXsectionWeightingFactor =2.086477e-03;
                            else if((pTHard >= 70) && (pTHard < 85) )
                                fXsectionWeightingFactor = 8.531533e-04;
                            else if((pTHard >= 85) && (pTHard < 99) )
                                fXsectionWeightingFactor = 3.152257e-04;
                            else if((pTHard >= 99) && (pTHard < 115) )
                                fXsectionWeightingFactor = 1.595340e-04;
                            else if((pTHard >= 115) && (pTHard < 132) )
                                fXsectionWeightingFactor = 7.689516e-05;
                            else if((pTHard >= 132) && (pTHard < 150) )
                                fXsectionWeightingFactor = 3.876765e-05;
                            else if((pTHard >= 150) && (pTHard < 169) )
                                fXsectionWeightingFactor = 2.078694e-05;
                            else if((pTHard >= 169) && (pTHard < 190) )
                                fXsectionWeightingFactor = 1.196729e-05;
                            else if((pTHard >= 190) && (pTHard < 212) )
                                fXsectionWeightingFactor =6.165574e-06;
                            else if((pTHard >= 212) && (pTHard < 235) )
                                fXsectionWeightingFactor = 3.415262e-06;
                            else if((pTHard >= 235) && (pTHard < 1000000) )
                                fXsectionWeightingFactor = 4.853493e-06;
                            else
                                fXsectionWeightingFactor = 0;


                        }

                }
        }

    Printf("fXsectionWeightingFactor %e",fXsectionWeightingFactor);
    FillHist("fh1dPtHardMonitor",fPtHard,fXsectionWeightingFactor);

    Bool_t HasImpactParameter = kFALSE;
    Double_t dca[2] = {-99999,-99999};
    Double_t cov[3] = {-99999,-99999,-99999};
    Double_t TrackWeight       = 1;
    AliVTrack* trackV = NULL;
    fIsEsd =  (InputEvent()->IsA()==AliESDEvent::Class())? kTRUE : kFALSE;
    EventwiseCleanup();
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

            IncHist("fh1dTracksAccepeted",2);
            FillHist("fh2dAcceptedTracksEtaPhi",trackV->Eta(),trackV->Phi(), this->fXsectionWeightingFactor );
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
                    TrackWeight *= GetMonteCarloCorrectionFactor(trackV,corrpartidx);
                }
            FillTrackHistograms(trackV,dca,cov,TrackWeight);
        }
    AliJetContainer *  jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    if (!jetconrec) return kFALSE;
    AliJetContainer * jetcongen = nullptr;
    AliEmcalJet * jetgen  = nullptr;
    if(fIsPythia){
            jetcongen = static_cast<AliJetContainer*>(fJetCollArray.At(1));
            if(!MatchJetsGeometricDefault()) AliInfo("Jet matching did not succeed!");
            jetcongen->ResetCurrentID();
            while ((jetgen = jetcongen->GetNextJet()))
                {
                    if (!jetgen) continue;
                    Int_t jetflavour =0;
                    Bool_t is_udgjet = kFALSE;
                    //Event based association to save memory
                    jetflavour =IsMCJetPartonFast(jetgen,0.4,is_udgjet);
                    FillHist("fh1dJetGenPt",GetPtCorrectedMC(jetgen), this->fXsectionWeightingFactor);
                    if(jetflavour ==0)      FillHist("fh1dJetGenPtUnidentified",GetPtCorrectedMC(jetgen), this->fXsectionWeightingFactor );
                    else if(jetflavour ==1) FillHist("fh1dJetGenPtudsg",GetPtCorrectedMC(jetgen), this->fXsectionWeightingFactor );
                    else if(jetflavour ==2) FillHist("fh1dJetGenPtc",GetPtCorrectedMC(jetgen), this->fXsectionWeightingFactor );
                    else if(jetflavour ==3) FillHist("fh1dJetGenPtb",GetPtCorrectedMC(jetgen), this->fXsectionWeightingFactor );
                }
            jetcongen->ResetCurrentID();
            jetconrec->ResetCurrentID();
        }
    // Loop over reconstructed/matched jets for template creation and analysis
    AliEmcalJet * jetrec  = nullptr;
    AliEmcalJet * jetmatched  = nullptr;
    jetconrec->ResetCurrentID();
    Double_t jetpt=0;
    while ((jetrec = jetconrec->GetNextJet()))
        {
            if(!jetrec) break;
            double val;
            jetpt = jetrec->Pt();
            if(!(jetconrec->GetRhoParameter() == nullptr))
                {
                    jetpt = jetpt - jetconrec->GetRhoVal() * jetrec->Area();
                }
            if(fIsPythia){
                    if (jetrec->MatchedJet()) {
                            Double_t genpt = jetrec->MatchedJet()->Pt();
                            if(!(jetcongen->GetRhoParameter() == nullptr))
                                {
                                    genpt = genpt - jetcongen->GetRhoVal() * jetrec->MatchedJet()->Area();
                                }
                            FillHist("fh2dJetGenPtVsJetRecPt",genpt,jetpt, this->fXsectionWeightingFactor );
                        }
                }
            // make inclusive signed imp. parameter constituent histograms
            Double_t dca[2] = {-99999,-99999};
            Double_t cov[3] = {-99999,-99999,-99999};
            Double_t sign=0;
            Int_t jetflavour=0;
            Bool_t is_udgjet = kFALSE;
            if(fIsPythia){
                    jetmatched = nullptr;
                    jetmatched =jetrec->MatchedJet();
                    if(jetmatched)jetflavour = IsMCJetPartonFast(jetmatched,0.4,is_udgjet); //Event based association to save memory
                }
            FillHist("fh1dJetRecPt",jetrec->Pt(), this->fXsectionWeightingFactor );
            if(fIsPythia){
                    if(jetflavour==0)     FillHist("fh1dJetRecPtUnidentified",jetpt, this->fXsectionWeightingFactor );
                    else if(jetflavour==1)FillHist("fh1dJetRecPtudgs",          jetpt, this->fXsectionWeightingFactor );
                    else if(jetflavour==2)FillHist("fh1dJetRecPtc",           jetpt, this->fXsectionWeightingFactor );
                    else if(jetflavour==3)FillHist("fh1dJetRecPtb",           jetpt, this->fXsectionWeightingFactor );
                }
            if(!(fJetCutsHF->IsJetSelected(jetrec))) break;
            FillHist("fh1dJetRecEtaPhiAccepted",jetrec->Eta(),jetrec->Phi(), this->fXsectionWeightingFactor );
            FillHist("fh1dJetRecPtAccepted",jetpt, this->fXsectionWeightingFactor );
            if(fIsPythia){
                    if(jetflavour==0)     FillHist("fh1dJetRecPtUnidentifiedAccepted",jetpt, this->fXsectionWeightingFactor );
                    else if(jetflavour==1)FillHist("fh1dJetRecPtudsgAccepted",        jetpt, this->fXsectionWeightingFactor );
                    else if(jetflavour==2)FillHist("fh1dJetRecPtcAccepted",           jetpt, this->fXsectionWeightingFactor );
                    else if(jetflavour==3)FillHist("fh1dJetRecPtbAccepted",           jetpt, this->fXsectionWeightingFactor );
                }

            std::vector<SJetIpPati> sImpParXY,sImpParXYZ,sImpParXYSig,sImpParXYZSig;

            double jetprob = CalculateJetProb(jetrec);



            const char * subtype_jp [4] = {"","udsg","c","b"};
            //Printf("JetProbability %e (flavour %s)", jetprob,subtype_jp[jetflavour]);

            if(jetflavour>0 && fIsPythia){
                    if(jetflavour==1) FillHist("fh2d_jetprob_light",jetpt,jetprob,this->fXsectionWeightingFactor);
                    else if(jetflavour==2) FillHist("fh2d_jetprob_charm",jetpt,jetprob,this->fXsectionWeightingFactor);
                    else if(jetflavour==3) FillHist("fh2d_jetprob_beauty",jetpt,jetprob,this->fXsectionWeightingFactor);


                    if(jetprob >0.5){
                            FillHist(Form("fh1dJetRecPt_0_5JP_%sAccepted",subtype_jp[jetflavour]),jetpt,this->fXsectionWeightingFactor);
                        }
                    if (jetprob >0.6){
                            FillHist(Form("fh1dJetRecPt_0_6JP_%sAccepted",subtype_jp[jetflavour]),jetpt,this->fXsectionWeightingFactor);
                        }
                    if (jetprob >0.7){
                            FillHist(Form("fh1dJetRecPt_0_7JP_%sAccepted",subtype_jp[jetflavour]),jetpt,this->fXsectionWeightingFactor);
                        }
                    if (jetprob >0.8){
                            FillHist(Form("fh1dJetRecPt_0_8JP_%sAccepted",subtype_jp[jetflavour]),jetpt,this->fXsectionWeightingFactor);
                        }
                    if (jetprob >0.9){
                            FillHist(Form("fh1dJetRecPt_0_9JP_%sAccepted",subtype_jp[jetflavour]),jetpt,this->fXsectionWeightingFactor);
                        }
                    if (jetprob >0.95){
                            FillHist(Form("fh1dJetRecPt_0_95JP_%sAccepted",subtype_jp[jetflavour]),jetpt,this->fXsectionWeightingFactor);
                        }
                }



            if(jetprob >0.5){
                    FillHist(Form("fh1dJetRecPt_0_5JP_%sAccepted","all"),jetpt,this->fXsectionWeightingFactor);
                }
            if (jetprob >0.6){
                    FillHist(Form("fh1dJetRecPt_0_6JP_%sAccepted","all"),jetpt,this->fXsectionWeightingFactor);
                }
            if (jetprob >0.7){
                    FillHist(Form("fh1dJetRecPt_0_7JP_%sAccepted","all"),jetpt,this->fXsectionWeightingFactor);
                }
            if (jetprob >0.8){
                    FillHist(Form("fh1dJetRecPt_0_8JP_%sAccepted","all"),jetpt,this->fXsectionWeightingFactor);
                }
            if (jetprob >0.9){
                    FillHist(Form("fh1dJetRecPt_0_9JP_%sAccepted","all"),jetpt,this->fXsectionWeightingFactor);
                }
            if (jetprob >0.95){
                    FillHist(Form("fh1dJetRecPt_0_95JP_%sAccepted","all"),jetpt,this->fXsectionWeightingFactor);
                }


            for(Int_t itrack = 0; itrack < InputEvent()->GetNumberOfTracks(); ++itrack)
                {
                    TrackWeight=1;
                    Double_t xyzatcda[3];
                    AliAODTrack * trackV = (AliAODTrack *) InputEvent()->GetTrack(itrack);
                    if (!trackV || !jetrec)            continue;
                    if (jetrec->DeltaR(trackV) > 0.4) continue;


                    //////////////////////BEGIN jet probability part
                    ///
                    /// Stores the signed ip distributions for UDG jets as input for the JP resolution function
                    if (!IsTrackAccepted((AliAODTrack*)trackV,3))   continue;
                    if(GetImpactParameterWrtToJet((AliAODTrack*)trackV,(AliAODEvent*)InputEvent(),jetrec,dca,cov,xyzatcda,sign)){
                            if(fEventVertex) {
                                    delete fEventVertex;
                                    fEventVertex =nullptr;
                                }


                            Int_t corridx=-1;
                            fIsPythia ? TrackWeight = GetMonteCarloCorrectionFactor(trackV,corridx) : TrackWeight =1;
                            dca[0]=fabs(dca[0]);
                            Double_t cursImParXY     =TMath::Abs(GetValImpactParameter(   kXY,dca,cov))*sign;
                            Double_t cursImParXYSig  =TMath::Abs(GetValImpactParameter(kXYSig,dca,cov))*sign;
                            Double_t cursImParXYZ    =TMath::Abs(GetValImpactParameter(   kXYZ,dca,cov))*sign;
                            Double_t cursImParXYZSig =TMath::Abs(GetValImpactParameter(kXYZSig,dca,cov))*sign;

                            if(is_udgjet){
                                    if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                                            FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShits",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                            FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                        }
                                    else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                                            FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShits",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                            FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                        }
                                    else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                                            FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShits",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                            FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                        }
                                    else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                                            FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShits",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                            FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShits",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                        }
                                }
                            if(IsFromElectron((AliAODTrack*)trackV)){
                                    if(is_udgjet){
                                            if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsElectrons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsElectrons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                        }
                                }
                            else if(IsFromPion((AliAODTrack*)trackV)){
                                    if(is_udgjet){
                                            if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsPions",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsPions",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                        }
                                }
                            else if(IsFromKaon((AliAODTrack*)trackV)){
                                    if(is_udgjet){
                                            if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsKaons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsKaons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                        }
                                }
                            else if(IsFromProton((AliAODTrack*)trackV)){
                                    if(is_udgjet){
                                            if (IsTrackAcceptedJP((AliAODTrack*)trackV,6)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_6ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,5)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_5ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_5ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,4)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_4ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_4ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                            else if (IsTrackAcceptedJP((AliAODTrack*)trackV,3)){
                                                    FillHist("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_3ITShitsProtons",trackV->Pt(),cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                                    FillHist("fh2dJetSignedImpParXYZSignificanceudg_light_resfunction_3ITShitsProtons",trackV->Pt(),cursImParXYZSig,TrackWeight*this->fXsectionWeightingFactor );
                                                }
                                        }
                                }

                            if (!IsTrackAccepted((AliAODTrack*)trackV,6))   continue;

                            //Fill jet probability ipsig histograms for template fitting
                            const char * subtype_jp [4] = {"","udsg","c","b"};
                            if(jetflavour>0 && fIsPythia){
                                    if(jetprob >0.5){
                                            FillHist(Form("fh2d_ImpSigXY_%s_0_5JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                            FillHist(Form("fh2d_ImpSigXYZ_%s_0_5JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                        }
                                    if (jetprob >0.6){
                                            FillHist(Form("fh2d_ImpSigXY_%s_0_6JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                            FillHist(Form("fh2d_ImpSigXYZ_%s_0_6JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                        }
                                    if (jetprob >0.7){
                                            FillHist(Form("fh2d_ImpSigXY_%s_0_7JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                            FillHist(Form("fh2d_ImpSigXYZ_%s_0_7JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                        }
                                    if (jetprob >0.8){
                                            FillHist(Form("fh2d_ImpSigXY_%s_0_8JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                            FillHist(Form("fh2d_ImpSigXYZ_%s_0_8JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                        }
                                    if (jetprob >0.9){
                                            FillHist(Form("fh2d_ImpSigXY_%s_0_9JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                            FillHist(Form("fh2d_ImpSigXYZ_%s_0_9JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                        }
                                    if (jetprob >0.95){
                                            FillHist(Form("fh2d_ImpSigXY_%s_0_95JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                            FillHist(Form("fh2d_ImpSigXYZ_%s_0_95JP",subtype_jp[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                        }
                                }

                            if(jetprob >0.5){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_5JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_5JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                }
                            if (jetprob >0.6){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_6JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_6JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                }
                            if (jetprob >0.7){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_7JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_7JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                }
                            if (jetprob >0.8){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_8JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_8JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                }
                            if (jetprob >0.9){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_9JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_9JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                }
                            if (jetprob >0.95){
                                    FillHist(Form("fh2d_ImpSigXY_%s_0_95JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                    FillHist(Form("fh2d_ImpSigXYZ_%s_0_95JP","all"),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor);
                                }
                        }
                    //////////////////////END jet probability part

                    if(GetImpactParameterWrtToJet((AliAODTrack*)trackV,(AliAODEvent*)InputEvent(),jetrec,dca,cov,xyzatcda,sign)){
                            if(fEventVertex) {
                                    delete fEventVertex;
                                    fEventVertex =nullptr;
                                }

                            Int_t corridx=-1;
                            fIsPythia ? TrackWeight = GetMonteCarloCorrectionFactor(trackV,corridx) : TrackWeight =1;
                            dca[0]=fabs(dca[0]);
                            Double_t cursImParXY     =TMath::Abs(GetValImpactParameter(   kXY,dca,cov))*sign;
                            Double_t cursImParXYSig  =TMath::Abs(GetValImpactParameter(kXYSig,dca,cov))*sign;
                            Double_t cursImParXYZ    =TMath::Abs(GetValImpactParameter(   kXYZ,dca,cov))*sign;
                            Double_t cursImParXYZSig =TMath::Abs(GetValImpactParameter(kXYZSig,dca,cov))*sign;


                            FillHist("fh2dJetSignedImpParXY"            ,jetpt,cursImParXY,TrackWeight*this->fXsectionWeightingFactor );
                            FillHist("fh2dJetSignedImpParXYSignificance",jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );

                            const char * subtype [4] = {"Unidentified","udsg","c","b"};
                            if(fIsPythia){
                                    FillHist(Form("fh2dJetSignedImpParXY%s",subtype[jetflavour]),jetpt,cursImParXY,TrackWeight*this->fXsectionWeightingFactor );
                                    FillHist(Form("fh2dJetSignedImpParXYSignificance%s",subtype[jetflavour]),jetpt,cursImParXYSig,TrackWeight*this->fXsectionWeightingFactor );
                                }
                            SJetIpPati a(cursImParXY, TrackWeight,kFALSE,kFALSE,trackV->GetLabel()); sImpParXY.push_back(a);
                            SJetIpPati b(cursImParXYZ, TrackWeight,kFALSE,kFALSE,trackV->GetLabel()); sImpParXYZ.push_back(b);
                            SJetIpPati c(cursImParXYSig, TrackWeight,kFALSE,kFALSE,trackV->GetLabel());sImpParXYSig.push_back(c);
                            SJetIpPati d(cursImParXYZSig, TrackWeight,kFALSE,kFALSE,trackV->GetLabel());sImpParXYZSig.push_back(d);

                        }
                }
            std::sort(sImpParXY.begin(),sImpParXY.end(),        AliAnalysisTaskHFJetIPQA::mysort);
            std::sort(sImpParXYSig.begin(),sImpParXYSig.end(),  AliAnalysisTaskHFJetIPQA::mysort);
            std::sort(sImpParXYZ.begin(),sImpParXYZ.end(),      AliAnalysisTaskHFJetIPQA::mysort);
            std::sort(sImpParXYZSig.begin(),sImpParXYZSig.end(),AliAnalysisTaskHFJetIPQA::mysort);
            const char * subtype[4] = {"Unidentified","udsg","c","b"};
            const char * subord [3] = {"First","Second","Third"};
            const char * stype  [4] = {"fh2dJetSignedImpParXY","fh2dJetSignedImpParXYSignificance","fh2dJetSignedImpParXYZ","fh2dJetSignedImpParXYZSignificance"};
            for (Int_t ot = 0 ; ot <3 ;++ot){
                    if ((int)sImpParXY.size()>ot){

                            if(ot==0) {
                                    if(jetflavour >0){
                                            FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,subtype[jetflavour]),jetpt,this->fXsectionWeightingFactor);
                                        }
                                    FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,"all"),jetpt,this->fXsectionWeightingFactor);
                                }
                            else if (ot==1)
                                {
                                    if(jetflavour >0){
                                            FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,subtype[jetflavour]),jetpt,this->fXsectionWeightingFactor);
                                        }
                                    FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,"all"),jetpt,this->fXsectionWeightingFactor);
                                }
                            else if (ot==2){
                                    if(jetflavour >0){
                                            FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,subtype[jetflavour]),jetpt,this->fXsectionWeightingFactor);
                                        }
                                    FillHist(Form("fh1dJetRecPt_n_%i_%s_Accepted",ot+1,"all"),jetpt,this->fXsectionWeightingFactor);
                                }
                            Double_t params [4] ={sImpParXY.at(ot).first,sImpParXYSig.at(ot).first,sImpParXYZ.at(ot).first,sImpParXYZSig.at(ot).first};
                            Double_t weights[4] ={sImpParXY.at(ot).second,sImpParXYSig.at(ot).second,sImpParXYZ.at(ot).second,sImpParXYZSig.at(ot).second};
                            for (Int_t ost = 0 ; ost <4 ;++ost){
                                    TString hname = Form("%s%s",stype[ost],subord[ot]);
                                    if(fIsPythia)   FillHist(hname.Data(),jetpt,params[ost],weights[ost] *  this->fXsectionWeightingFactor);
                                    else  FillHist(hname.Data(),jetpt,params[ost], this->fXsectionWeightingFactor );


                                }
                            if(fIsPythia){
                                    for (Int_t ost = 0 ; ost <4 ;++ost){
                                            TString hname = Form("%s%s%s",stype[ost],subtype[jetflavour],subord[ot]);
                                            FillHist(hname.Data(),jetpt,params[ost],weights[ost]* this->fXsectionWeightingFactor  );
                                        }
                                }
                        }
                }
            sImpParXY.clear();
            sImpParXYSig.clear();
            sImpParXYZ.clear();
            sImpParXYZSig.clear();
        }

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
    isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);

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

    if(fEventVertex->GetChi2perNDF()>3.5*3.5) {
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
/*
Bool_t AliAnalysisTaskHFJetIPQA::IsEventSelected()	{

    Int_t WhyRejected =0;
    ULong_t RejectionBits=0;
    if(fIsPythia && fIsEsd){
            AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
            if(!mcH ){
                    AliError("No MC Event Handler available");
                    return kFALSE;
                }
            if(!mcH->InitOk()) return kFALSE;
            if(!mcH->TreeK()) return kFALSE ;
            if(!mcH->TreeTR()) return kFALSE;
        }
    fEventVertex = nullptr;
    fEventVertex = dynamic_cast<const AliAODVertex*>(InputEvent()->GetPrimaryVertex());//GetKFPrimaryVertex());
    if(!fEventVertex) return kFALSE;
    AliVEvent* eev = InputEvent();
    if(eev && fEventVertex && fEventVertex->GetNContributors()>0){
            FillHist("fh1dVertexZ",eev->GetPrimaryVertex()->GetZ(),1);
            Double_t vtxx =fEventVertex->GetX();
            Double_t vtxy = fEventVertex->GetY();
            FillHist("fh1dVertexR",vtxx,vtxy,1);
        }else return kFALSE;
    if(!(IsSelected(eev,WhyRejected,RejectionBits)))
        {
            IncHist("fh1dEventRejectionRDHFCuts",2);
            if(WhyRejected==kPhysicsSelection)         IncHist("fh1dEventRejectionRDHFCuts",3);
            else if(WhyRejected==kNoVertex)            IncHist("fh1dEventRejectionRDHFCuts",4);
            else if(WhyRejected==kNoContributors)      IncHist("fh1dEventRejectionRDHFCuts",5);
            else if(WhyRejected==kTooFewVtxContrib)    IncHist("fh1dEventRejectionRDHFCuts",6);
            else if(WhyRejected==kZVtxOutFid)          IncHist("fh1dEventRejectionRDHFCuts",7);
            else if(WhyRejected==kBadDiamondXDistance) IncHist("fh1dEventRejectionRDHFCuts",8);
            else if(WhyRejected==kBadDiamondYDistance) IncHist("fh1dEventRejectionRDHFCuts",9);
            else if(WhyRejected==kBadDiamondZDistance) IncHist("fh1dEventRejectionRDHFCuts",10);
            else if(WhyRejected==kVertexChi2NDF)       IncHist("fh1dEventRejectionRDHFCuts",11);
            return kFALSE;
        }else {



            Double_t vtxx =fEventVertex->GetX();
            Double_t vtxy = fEventVertex->GetY();
            Double_t vtxz = fEventVertex->GetZ();
            FillHist("fh1dVertexXvsMultiplicity",vtxx,eev->GetNumberOfTracks(),1);
            FillHist("fh1dVertexYvsMultiplicity",vtxy,eev->GetNumberOfTracks(),1);
            FillHist("fh1dVertexZvsMultiplicity",vtxz,eev->GetNumberOfTracks(),1);
            IncHist("fh1dEventRejectionRDHFCuts",1);
            FillHist("fh1dVertexZAccepted",fEventVertex->GetZ(),1);
            FillHist("fh1dVertexRAccepted",vtxx,vtxy,1);
            FillHist("fh2dVertexChi2NDFNESDTracks",fEventVertex->GetChi2perNDF(),eev->GetNumberOfTracks(),1);
            return kTRUE;
        }


    return kTRUE;
}
*/

void AliAnalysisTaskHFJetIPQA::UserCreateOutputObjects(){
    AliAnalysisTaskEmcal::UserCreateOutputObjects();
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if(!mgr) AliError("Analysis manager not found!");
    AliVEventHandler *evhand = mgr->GetInputEventHandler();
    if (!evhand) AliError("Event handler not found!");
    if (evhand->InheritsFrom("AliESDInputHandler"))  fIsEsd = kTRUE;
    else fIsEsd = kFALSE;
    OpenFile(1);
    Double_t xomega[10] ={0,1.00993,1.40318,1.7428,2.04667,2.35055,2.77061,3.39623,4.40616,1000};
    Double_t yomega[10] ={2.74011,2.74011,3.16949,3.37288,4.02825,4.38983,4.23164,4.75141,3.62147,3.62147};
    Double_t xxi[18] ={0,0.723932,0.849057,0.947368,1.04568,1.14399,1.25124,1.35849,1.43893,1.60874,1.80536,2.04667,2.41311,
                       2.87786,3.50348,4.39722,5.45184,1000};
    Double_t yxi[18] ={1.20339,1.20339,1.45198,1.54237,1.76836,1.81356,1.85876,1.97175,2.10734,2.10734,2.15254,2.19774,2.22034,
                       2.19774,2.12994,1.83616,1.36158,1.36158};
    Double_t xK0s[24]= {0,0.0628272,0.162304,0.26178,0.356021,0.465969,0.570681,0.675393,0.743455,0.863874,0.947644,1.10471,1.29843,
                        1.50785,1.70681,1.91099,2.19895,2.60733,3.01571,3.40838,3.82199,4.51832,5.49215,1000};
    Double_t yK0s[24]= {1.31496,1.31496,1.19685,1.11024,1.16535,1.11811,1.07874,1.05512,
                        1.01575,0.976378,0.92126,0.889764,0.858268,0.811024,0.84252,0.858268,
                        0.811024,0.811024,0.779528,0.787402,0.795276,0.811024,0.88189,0.88189};
    Double_t xPhi[24] ={0,0.456294,0.54021,0.645105,0.76049,0.833916,0.938811,1.03322,1.15385,
                        1.23776,1.34266,1.46853,1.5472,1.6521,1.75175,1.86189,1.96154,2.09266,
                        2.27622,2.51748,2.71678,2.91084,3.2465,1000};
    Double_t yPhi[24] ={0.699725,0.699725,0.699725,0.721763,0.661157,0.683196,0.650138,
                        0.639118,0.639118,0.62259,0.61157,0.62259,0.61708,0.61708,0.606061,
                        0.644628,0.694215,0.699725,0.721763,0.787879,0.77686,0.809917,0.870523,0.870523};
    //This is MC /Data not Data over MC
    fGraphOmega = new TGraph(10,xomega,yomega);
    fGraphXi    = new TGraph(18,xxi,yxi);
    fK0Star     = new TGraph(24,xK0s,yK0s);
    fPhi        = new TGraph(24,xPhi,yPhi);

    const Int_t nBins3dSignificance = 500;
    const Int_t nBins3d = 250;

    Double_t lowIPxy =-1.;
    Double_t highIPxy =1.;
    OpenFile(1);

    if(!fPidResponse)   fPidResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPidResponse) {
            AliFatal("NULL PID response");
        }
    if(!fCombined) fCombined = new AliPIDCombined();
    //Make Graphs
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
    //ADD HISTOGRAMS
    TH1D * h = (TH1D*)AddHistogramm("fh1dEventRejectionRDHFCuts","fh1dEventRejectionRDHFCuts;reason;count",12,0,12);
    h->GetXaxis()->SetBinLabel(1,"Event accepted");
    h->GetXaxis()->SetBinLabel(2,"Event rejected");
    h->GetXaxis()->SetBinLabel(3,"Wrong physics selection");
    h->GetXaxis()->SetBinLabel(4,"No vertex");
    h->GetXaxis()->SetBinLabel(5,"No contributors");
    h->GetXaxis()->SetBinLabel(6,"Less than 10 contributors");
    h->GetXaxis()->SetBinLabel(7,">10cm vertex Z distance");
    h->GetXaxis()->SetBinLabel(8,"Bad diamond X distance");
    h->GetXaxis()->SetBinLabel(9,"Bad diamond Y distance");
    h->GetXaxis()->SetBinLabel(10,"Bad diamond Z distance");
    h->GetXaxis()->SetBinLabel(11,"Chi2 vtx >1.5 ");
    AddHistogramm("fh1dTracksAccepeted","# tracks before/after cuts;;",3,0,3);

    TH1D * h1 = GetHist1D("fh1dTracksAccepeted");
    h1->GetXaxis()->SetBinLabel(1,"total");
    h1->GetXaxis()->SetBinLabel(2,"accepted");
    h1->GetXaxis()->SetBinLabel(3,"rejected");
    // Tracks impact parameter histograms
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
    fHistManager.CreateTH2("fh1dJetRecEtaPhiAccepted","detector level jet;#eta;phi",1,-0.5,0.5,1,0.,TMath::TwoPi(),"s");
    fHistManager.CreateTH2("fh2dAcceptedTracksEtaPhi","accepted tracks;#eta;phi",200,-0.9,0.9,200,0.,TMath::TwoPi(),"s");
    fHistManager.CreateTH1("fh1dJetRecPt","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPtAccepted","accepted detector level jets;pt (GeV/c); count",500,0,250,"s");
    if (fIsPythia){
            //Histograms for jet-probability tagger
            fHistManager.CreateTH2("fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShits", "fh2dJetSignedImpParXYSignificanceudg_light_resfunction_6ITShits;pt (GeV/c); count",200,0,100,500,-30,0,"s");
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

            fHistManager.CreateTH2("fh2dTracksImpParXY_McCorr","radial imp. parameter (after correction);impact parameter xy (cm);a.u.",2000,-1,1,500,0,100,"s");
            fHistManager.CreateTH1("fh1dTracksImpParXY_McCorr","radial imp. parameter (after correction);impact parameter xy (cm);a.u.",400,-0.2,0.2,"s");
            fHistManager.CreateTH1("fh1dTracksImpParXYZ_McCorr","3d imp. parameter (after correction);impact parameter 3d (cm);a.u.",2000,0,100.,"s");
            fHistManager.CreateTH2("fh2dTracksImpParXYZ_McCorr","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-1,1,500,0,100.,"s");
            fHistManager.CreateTH2("fh2dTracksImpParXYZSignificance_McCorr","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.,"s");
            fHistManager.CreateTH1("fh1dTracksImpParXYSignificance_McCorr","radial imp. parameter (after correction);impact parameter xy significance;a.u.",2000,-30,30.,"s");
            fHistManager.CreateTH1("fh1dTracksImpParXYZSignificance_McCorr","3d imp. parameter (after correction);impact parameter 3d significance;a.u.",2000,0.,100.,"s");
            fHistManager.CreateTH1("fh1dJetGenPt","generator level jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetGenPtUnidentified","generator level jets (no flavour assigned);pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetGenPtudsg","generator level udsg jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetGenPtc","generator level c jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetGenPtb","generator level b jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH2("fh2dJetGenPtVsJetRecPt","detector momentum response;gen pt;rec pt",500,0,250,500,0,250,"s");
            fHistManager.CreateTH1("fh1dJetRecPtudsg","detector level jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetRecPtUnidentified","detector level jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetRecPtc","detector level jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetRecPtb","detector level jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetRecPtUnidentifiedAccepted","detector level jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetRecPtudsgAccepted","detector level jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetRecPtcAccepted","detector level jets;pt (GeV/c); count",250,0,250,"s");
            fHistManager.CreateTH1("fh1dJetRecPtbAccepted","detector level jets;pt (GeV/c); count",250,0,250,"s");
        }
    fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_5JP", "fh2d_ImpSigXY_all_0_5JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
    fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_6JP", "fh2d_ImpSigXY_all_0_6JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
    fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_7JP", "fh2d_ImpSigXY_all_0_7JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
    fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_8JP", "fh2d_ImpSigXY_all_0_8JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
    fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_9JP", "fh2d_ImpSigXY_all_0_9JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");
    fHistManager.CreateTH2("fh2d_ImpSigXY_all_0_95JP", "fh2d_ImpSigXY_all_0_95JP;pt (GeV/c); sig",500,0,250,1000,-30,30,"s");

    fHistManager.CreateTH1("fh1dJetRecPt_0_5JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPt_0_6JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPt_0_7JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPt_0_8JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPt_0_9JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPt_0_95JP_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");

    //This is for the default baseline analysis

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
        }
    fHistManager.CreateTH1("fh1dJetRecPt_n_1_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPt_n_2_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");
    fHistManager.CreateTH1("fh1dJetRecPt_n_3_all_Accepted","detector level jets;pt (GeV/c); count",500,0,250,"s");




    const char * flavour[5]  = {"Unidentified","udsg","c","b",""};
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
    for (Int_t id = 0;id<2;++id)
        for (Int_t ifl = 0;ifl<5;++ifl)
            for (Int_t io = 0;io<4;++io)
                for (Int_t is = 0;is<1;++is)
                    for (Int_t it = 0;it<2;++it){
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
                            if(id==0)  ipbins =1000;//2000;
                            if((fIsPythia||(!fIsPythia && ifl==4)))  fHistManager.CreateTH2(Form("%s%s%s%s%s%s",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]),
                                                                                            Form("%s%s%s%s%s%s;;",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]),
                                                                                            ptbins,ptlow,pthigh,ipbins,iplow,iphigh,"s");
                        }
    TIter next(fHistManager.GetListOfHistograms());
    TObject* obj = 0;
    while ((obj = next())) {
            fOutput->Add(obj);
        }
    PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}
void AliAnalysisTaskHFJetIPQA::GetMaxImpactParameterCutR(const AliVTrack * const track, Double_t &maximpactRcut){
    //
    // Get max impact parameter cut r (pt dependent)
    //
    Double_t pt = track->Pt();
    if(pt > 0.15) {
            maximpactRcut = 0.0182 + 0.035/TMath::Power(pt,1.01);  // abs R cut
        }
    else maximpactRcut = 9999999999.0;
}


bool AliAnalysisTaskHFJetIPQA::GetPIDCombined(AliAODTrack * track, double  * prob, int &nDetectors,UInt_t &usedDet ,AliPID::EParticleType &MostProbablePID, bool setTrackPID ){
    //Initialize if neede

    //fPidResponse->SetITSPIDmethod(AliPIDResponse::kITSLikelihood);
    double nSigma[AliPIDResponse::kNdetectors][AliPID::kSPECIES];
    AliPIDResponse::EDetPidStatus status[AliPIDResponse::kNdetectors] = {AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,AliPIDResponse::kDetNoSignal,};

    unsigned int nGoodDet = 0;
    for (int j =0; j<AliPIDResponse::kNdetectors;j++)
        {
            for (int i =AliPID::kElectron; i<AliPID::kSPECIES;i++)
                {
                    double val = 0;

                    status[j] =  fPidResponse->NumberOfSigmas(static_cast <AliPIDResponse::EDetector>(j), track, static_cast <AliPID::EParticleType>(i), val);
                    if (status[j] == AliPIDResponse::kDetPidOk ){
                            nSigma[j][i] =val;
                            nGoodDet++;}
                    else nSigma[j][i] =-9999.;
                }
        }
    if( nGoodDet/7 <2 ) return false;

    nDetectors = nGoodDet/7;
    //Setup Combined PID

    Double_t probTPCTOF[AliPID::kSPECIES]={-1.,-1.,-1.,-1.,-1.};

    fCombined->SetDefaultTPCPriors();

    //Combine all available detectors ITS TPC TRD TOF EMCAL
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
    if(fIsPythia){
            int pdgm =0;
            int pdg =                    GetMCTruth(track,pdgm);
        }
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
    if(fIsPythia){
            int pdgm =0;
            int pdg =                    GetMCTruth(track,pdgm);
        }
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
    if(fIsPythia){
            int pdgm =0;
            int pdg =                    GetMCTruth(track,pdgm);
        }
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
    if(fIsPythia){
            int pdgm =0;
            int pdg =                    GetMCTruth(track,pdgm);
        }
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
    int tt=0;

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
            if(doskip){
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

    if(vtxESDNew->GetChi2toNDF()>3.5*3.5) {
            delete vtxESDNew; vtxESDNew=nullptr;
            return 0;
        }

    Double_t pos[3];
    Double_t cov[6];

    Double_t chi2perNDF;
    vtxESDNew->GetXYZ(pos); // position
    vtxESDNew->GetCovMatrix(cov); //covariance matrix
    chi2perNDF = vtxESDNew->GetChi2toNDF();
    if(vtxESDNew) delete vtxESDNew; vtxESDNew=NULL;
    AliAODVertex *vtxAODNew = new AliAODVertex(pos,cov,chi2perNDF);
    vtxAODNew->	SetNContributors(nContrib);
    return vtxAODNew;
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
            if(((AliAODTrack*)track)->GetNcls(1)<100) return kFALSE;
            ULong_t status = track->GetStatus();
            if(track->Pt()<1.)return kFALSE;
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
            //HasMatchedGoodTracklet((AliAODTrack*)track);
            if(!(((AliAODTrack*)track)->TestFilterBit(9) || ((AliAODTrack*)track)->TestFilterBit(4))){
                    return kFALSE;
                }
            if(!(((AliAODTrack*)track->HasPointOnITSLayer(0))&&(AliAODTrack*)track->HasPointOnITSLayer(1)))  {

                    //Printf ("SPD %i  i%",(int)((AliAODTrack*)track->HasPointOnITSLayer(0)),(int)((AliAODTrack*)track->HasPointOnITSLayer(1)));
                    return kFALSE;
                }
            if(((AliAODTrack*)track)->GetNcls(0)<abs(n)) return kFALSE;
            if(((AliAODTrack*)track)->GetNcls(1)<100) return kFALSE;
            ULong_t status = track->GetStatus();
            if(track->Pt()<.5)return kFALSE;
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
 * jet matching loop
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
    Int_t ipdg =abs(mcmother->PdgCode())	;
    return kFALSE;
}
/*! \brief Composition correction factor  getter
 *
 * finds the corresponding re-weighing factor for a certain track
 *
 */
Double_t AliAnalysisTaskHFJetIPQA::GetWeightFactor( AliVTrack * track,Int_t &pCorr_indx){
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
    Double_t val = 1;

    AliVParticle * mcpart = fIsEsd ? (   AliVParticle * )pMCESD:(   AliVParticle * )pMCAOD;
    Bool_t _particlesourcefound(kFALSE);
    Int_t  _particlesourcepdg(mcpart->PdgCode());
    Int_t  _particlesourceidx(25);
    Double_t _particlesourcept(0);

    AliVParticle * mcpartclone = mcpart;
    while(mcpart){//omega and xi test
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

    if (!_particlesourcefound) { //heavy mesons to improve templates
            mcpart = mcpartclone;
            while(mcpart){
                    if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;
                    _particlesourcept = mcpart->Pt();
                    _particlesourcepdg = abs(mcpart->PdgCode());
                    if (IsSelectionParticleALICE(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
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
                    if (IsSelectionParticleOmegaXiSigmaP(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
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


    /*if(abs(mcpart->PdgCode()) == bProton) factor *=0.80;
    if(abs(mcpart->PdgCode()) == bPi) factor *=0.80;
*/


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



    if (factor <= 0 || factor > 6.)  return 1;


    return factor ;




    /*if (factor <= 0 || factor > 6.)  factor = 1.;
    if(IsSecondaryFromWeakDecay(mcpartclone)){
            factor *= ff;
            return factor;
        }
    mother = GetVParticleMother(mcpart);
    if(mother){
            if(IsSecondaryFromWeakDecay(mother)){
                    factor *= ff;
                    return factor;
                }
            mother = GetVParticleMother(mcpart);
            if(mother){
                    if(IsSecondaryFromWeakDecay(mother)){
                            factor *= ff;
                            return factor;
                        }
                }
        }*/
    return factor;
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
            if(IsSecondaryFromWeakDecay(mcpart))return kTRUE;
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
        case bRhoPlus: //Experimental assume same shape correction for neutral and charged rho
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
    AliMCParticle * mother = nullptr;
    idx = -1;
    pdg = abs(mcpart->PdgCode());
    if(!IsPhysicalPrimary(mcpart))return kFALSE;
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
            Printf("bSigmaMinus");
            idx = bIdxSigmaMinus; // use lambda as proxy
            return kTRUE;
            break;
        case bSigmaPlus:
            idx = bIdxSigmaPlus; //use lambda as proxy
            return kTRUE;
            break;
        case bXiBaryon:
            idx = bIdxBStarPlus;//dummy! Externally solved
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
 * jet matching helper
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
Double_t AliAnalysisTaskHFJetIPQA::GetMonteCarloCorrectionFactor(AliVTrack* track,Int_t &pCorr_indx){
    double val=  GetWeightFactor(track,pCorr_indx);
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
    if(jet && jetconrec)
        return jet->Pt() - jetconrec->GetRhoVal() * jet->Area();
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
    if(jet && jetconrec)
        return jet->Pt() - jetconrec->GetRhoVal() * jet->Area();
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
/*! \brief IsMCJetPartonFast
 *
 * Fast jet parton MC matcher
 */
Int_t  AliAnalysisTaskHFJetIPQA::IsMCJetPartonFast(const AliEmcalJet *jet, Double_t radius,Bool_t &is_udg)
{
    if(!jet) return 0;
    if(! fMcEvtSampled){
            if(fIsEsd){
                    //Sample MC stack upward once per event and oly if there are jets
                    AliStack * Mstack = MCEvent()->Stack();
                    for(Int_t iPrim = 0 ; iPrim<Mstack->GetNprimary();iPrim++){
                            TParticle * part = (TParticle*)Mstack->Particle(iPrim);
                            if(!part) return 0;
                            if(!((part->GetStatusCode()==11) ||(part->GetStatusCode()==12))) continue;
                            Double_t etap = part->Eta();
                            Double_t phip = part->Phi();
                            Int_t pdg = (abs(part->GetPdgCode()));
                            if(pdg == 5) {
                                    fEtaBEvt.push_back(etap);
                                    fPhiBEvt.push_back(phip);
                                }
                            else if(pdg== 4) {
                                    fEtaCEvt.push_back(etap);
                                    fPhiCEvt.push_back(phip);
                                }
                            else if(pdg == 3 ) {
                                    fEtaSEvt.push_back(etap);
                                    fPhiSEvt.push_back(phip);
                                }
                            else if(pdg== 1 ||pdg== 2 || pdg== 3 || pdg == 21) {
                                    fEtaUdsgEvt.push_back(etap);
                                    fPhiUdsgEvt.push_back(phip);
                                }
                        }
                    fMcEvtSampled= kTRUE;
                }
            else
                {
                    for(Int_t iPrim = 0 ; iPrim<fMCArray->GetEntriesFast();iPrim++){
                            AliAODMCParticle * part = static_cast<AliAODMCParticle*>(fMCArray->At(iPrim));
                            if(!part) return 0;
                            if(!part->IsPrimary()) continue;
                            if(!((part->GetStatus()==11) ||(part->GetStatus()==12))) continue;
                            Double_t etap = part->Eta();
                            Double_t phip = part->Phi();
                            Int_t pdg = (abs(part->PdgCode()));
                            if(pdg == 5) {
                                    fEtaBEvt.push_back(etap);
                                    fPhiBEvt.push_back(phip);
                                }
                            else if(pdg== 4) {
                                    fEtaCEvt.push_back(etap);
                                    fPhiCEvt.push_back(phip);
                                }
                            else if(pdg == 3 ) {
                                    fEtaSEvt.push_back(etap);
                                    fPhiSEvt.push_back(phip);
                                }
                            else if(pdg== 1 ||pdg== 2 || pdg== 3 || pdg == 21) {
                                    fEtaUdsgEvt.push_back(etap);
                                    fPhiUdsgEvt.push_back(phip);
                                }
                        }
                    fMcEvtSampled= kTRUE;


                }
        }
    if(fEtaBEvt.size() ==0&& fEtaCEvt.size()==0&& fEtaSEvt.size()==0&&fEtaUdsgEvt.size()==0) return 0; //udsg
    Double_t etajet = jet->Eta();
    Double_t phijet = jet->Phi();
    //check for c jet
    for (Int_t icj = 0 ; icj <(Int_t)fEtaCEvt.size();++icj ){
            Double_t eta =fEtaCEvt.at(icj);
            Double_t phi =fPhiCEvt.at(icj);
            Double_t deta = etajet - eta;
            Double_t dphi = phijet-phi;
            dphi = TVector2::Phi_mpi_pi(dphi);
            Double_t  d = sqrt(deta * deta + dphi * dphi);
            if(d < radius) return 2;
        }
    //check for b jet
    for (Int_t icj = 0 ; icj <(Int_t)fEtaBEvt.size();++icj ){
            Double_t eta =fEtaBEvt.at(icj);
            Double_t phi =fPhiBEvt.at(icj);
            Double_t deta = etajet - eta;
            Double_t dphi = phijet - phi;
            dphi = TVector2::Phi_mpi_pi(dphi);
            Double_t  d = sqrt(deta * deta + dphi * dphi);
            if(d < radius) return 3;
        }
    //check for s jet
    for (Int_t icj = 0 ; icj <(Int_t)fEtaSEvt.size();++icj ){
            Double_t eta =fEtaSEvt.at(icj);
            Double_t phi =fPhiSEvt.at(icj);

            Double_t deta = etajet - eta;
            Double_t dphi = phijet - phi;
            dphi = TVector2::Phi_mpi_pi(dphi);
            Double_t  d = sqrt(deta * deta + dphi * dphi);
            if(d < radius) return 1;
        }
    for (Int_t icj = 0 ; icj <(Int_t)fEtaUdsgEvt.size();++icj ){
            Double_t eta =fEtaUdsgEvt.at(icj);
            Double_t phi =fPhiUdsgEvt.at(icj);

            Double_t deta = etajet - eta;
            Double_t dphi = phijet - phi;
            dphi = TVector2::Phi_mpi_pi(dphi);
            Double_t  d = sqrt(deta * deta + dphi * dphi);
            if(d < radius){
                    is_udg =kTRUE;
                    return 1;
                }
        }

    return 0;
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

    fOutput->Add(phist);
    Printf("Adding %s to output list",phist->GetName());
    return (TH1*)phist;
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

void AliAnalysisTaskHFJetIPQA::setFProductionNumberPtHard(const Int_t &value)
{
    fProductionNumberPtHard = value;
}



/*! \brief SubtractMean
 *
 *
 */



void AliAnalysisTaskHFJetIPQA::SubtractMean(Double_t val[], AliVTrack *track){
    Double_t  deltamean=fGraphMean->Eval(track->Pt() < 3. ? track->Pt() : 3 );//(fCurrentMeanFactors[0]-fCurrentMeanFactors[1]*TMath::Exp(-1*fCurrentMeanFactors[2]*(track->Pt()-fCurrentMeanFactors[3]))) *1e-4;

    val[0] -=deltamean*1e-4;;
}
//Helpers


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
    Double_t pv_pos[3];
    Double_t p_at_dca[3];
    //Classic dca calculation
    if(etp.PropagateToDCA(vtxESDSkip, event->GetMagneticField(), kBeampiperadius, dca, cov)){
            success = kTRUE;
            etp.GetXYZ(XYZatDCA);
            etp.GetXYZ(x_at_dca);
            etp.GetPxPyPz(p_at_dca);
         if(fIsPythia)   dca[0] += fMCglobalDCAxyShift; // generic mean offset in LHC10e default is 0.007 == 7 µm
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
    amvf->_fitter->_run_fitter_jet(jet,(AliAODEvent*)InputEvent());
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

    double E1  = ROOT::Math::Similarity(deriv,E66Track);
    double E2  = ROOT::Math::Similarity(deriv_v,E33Vertex);
    double theError =sqrt(E1+E2);

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
    double  bdecaylength =  TMath::Sqrt((VxVyVz[0] - xyz_track_global[0]) * (VxVyVz[0] - xyz_track_global[0]) +
            (VxVyVz[1] - xyz_track_global[1]) * (VxVyVz[1] - xyz_track_global[1])+
            (VxVyVz[2] - xyz_track_global[2]) * (VxVyVz[2] - xyz_track_global[2]));

    double dcajetrack = TMath::Sqrt((xyz_jet_global[0] - xyz_track_global[0]) * (xyz_jet_global[0] - xyz_track_global[0]) +
            (xyz_jet_global[1] - xyz_track_global[1]) * (xyz_jet_global[1] - xyz_track_global[1])+
            (xyz_jet_global[2] - xyz_track_global[2]) * (xyz_jet_global[2]- xyz_track_global[2]));

    if(bdecaylength>bAnalysisCut_MaxDecayLength) return kFALSE;
    if(dcajetrack  >bAnalysisCut_DCAJetTrack) return kFALSE;
    return kTRUE;
}

Double_t AliAnalysisTaskHFJetIPQA::CalculateJetProb(AliEmcalJet *jet)
{
    if(!jet) return -1.;
    Double_t retval = -1;
    //Loop over all tracks calculate P(s) for all accepted later add looser cuts also
    Int_t ntracks = (Int_t)jet->GetNumberOfTracks();
    Double_t prodPS = 1;
    Double_t curps=-1;
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
