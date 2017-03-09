#include "TList.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliKFVertex.h"
#include "AliAnalysisUtils.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "TParticle.h"
#include "AliAODMCParticle.h"
#include "AliStack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliVEventHandler.h"
#include "AliStack.h"
#include "TGraph.h"
#include "AliTracker.h"
#include "AliAODMCHeader.h"
#include "AliJetContainer.h"
#include "AliPicoTrack.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TTree.h"
#include "AliRDHFJetsCuts.h"
#include "AliAnalysisTaskHFJetIPQA.h"
#include "AliExternalTrackParam.h"
#include "AliGenEventHeader.h"
#include "AliVertexerTracks.h"
#include "AliEmcalList.h"
#include "AliHFJetsTagging.h"
#include "AliESDUtils.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include "AliTriggerAnalysis.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliGenPythiaEventHeader.h"
#include "AliOADBContainer.h"
#include "TFile.h"
#include "TFormula.h"
#include "TGrid.h"
using std::min;
using std::cout;
using std::endl;
using std::vector;
using std::pair;
ClassImp(AliAnalysisTaskHFJetIPQA)
// ######################################################################################## CONSTRUCTORS
AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA(): AliAnalysisTaskEmcalJet("AliAnalysisTaskHFJetIPQA", kFALSE),
    fGraphMean(0x0),
    fGraphSigmaData(0x0),
    fGraphSigmaMC(0x0),
    fGraphXi(0x0),
    fGraphOmega(0x0),
    fK0Star(0x0),
    fPhi(0x0),
    fGeant3FlukaProton(0x0),
    fGeant3FlukaAntiProton(0x0),
    fGeant3FlukaLambda(0x0),
    fGeant3FlukaAntiLambda(0x0),
    fGeant3FlukaKMinus(0x0),
    fOutput2(0x0),
    fMCArray(0x0),
    fJetCutsHF(new AliRDHFJetsCuts()),
    fMCEvent(0x0),
    fESDTrackCut(0x0),
    fUtils(new AliAnalysisUtils()),
    fVertexer(0x0),
    fESD(kFALSE),
    fMcEvtSampled(kFALSE),
    fCorrrectionSamplingMode(kFALSE),
    fItsClustersInputGlobal(6),
    fBackgroundFactorLinus{0},
    fEtaSEvt(100),
    fPhiSEvt(100),
    fEtaBEvt(100),
    fPhiBEvt(100),
    fEtaCEvt(100),
    fPhiCEvt(100),
    fEtaUdsgEvt(100),
    fPhiUdsgEvt(100),
    fRandom(0x0),
    fh2dAcceptedTracksEtaPhiPerLayer{0,0,0,0,0,0}
{

    fOutput2 =0x0;
    fDisableWeightingMC= kFALSE;
    fSkipMeanSigmaCorrection= kFALSE;
    DefineOutput(1,  TList::Class()) ;

    for(Int_t i =0 ; i<498;++i)for(Int_t j =0 ; j<19;++j)  fBackgroundFactorLinus[j][i]=1.;
}
//###############################################################################################################
AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA(const char *name): AliAnalysisTaskEmcalJet(name,kFALSE),
    fGraphMean(0x0),
    fGraphSigmaData(0x0),
    fGraphSigmaMC(0x0),
    fGraphXi(0x0),
    fGraphOmega(0x0),
    fK0Star(0x0),
    fPhi(0x0),
    fGeant3FlukaProton(0x0),
    fGeant3FlukaAntiProton(0x0),
    fGeant3FlukaLambda(0x0),
    fGeant3FlukaAntiLambda(0x0),
    fGeant3FlukaKMinus(0x0),
    fOutput2(0x0),
    fMCArray(0x0),
    fJetCutsHF(new AliRDHFJetsCuts()),
    fMCEvent(0x0),
    fESDTrackCut(0x0),
    fUtils(new AliAnalysisUtils()),
    fVertexer(0x0),
    fESD(kFALSE),
    fMcEvtSampled(kFALSE),
    fCorrrectionSamplingMode(kFALSE),
    fItsClustersInputGlobal(6),
    fBackgroundFactorLinus{0},
    fEtaSEvt(100),
    fPhiSEvt(100),
    fEtaBEvt(100),
    fPhiBEvt(100),
    fEtaCEvt(100),
    fPhiCEvt(100),
    fEtaUdsgEvt(100),
    fPhiUdsgEvt(100),
    fRandom(0x0),
    fh2dAcceptedTracksEtaPhiPerLayer{0,0,0,0,0,0}
{

    fDisableWeightingMC= kFALSE;
    fSkipMeanSigmaCorrection= kTRUE;
    fOutput2 =0x0;
    DefineOutput(1,  TList::Class()) ;

    for(Int_t i =0 ; i<498;++i)for(Int_t j =0 ; j<19;++j)  fBackgroundFactorLinus[j][i]=1; //set default to 1
}


void AliAnalysisTaskHFJetIPQA::setN_ITSClusters_Input_global( Int_t value)
{
    fItsClustersInputGlobal = value;
}


Bool_t AliAnalysisTaskHFJetIPQA::Notify()
{
    return AliAnalysisTaskEmcalJet::Notify();
}


void AliAnalysisTaskHFJetIPQA::localtoglobal(Double_t alpha ,Double_t* local,Double_t* global)
{
    global[0] = local[0]*  TMath::Sin(alpha) + local[1] * TMath::Cos(alpha);
    global[1] = local[0]*  TMath::Cos(alpha) + local[1] * TMath::Sin(-alpha);
    global[2] = local[2];
    return;

}

void AliAnalysisTaskHFJetIPQA::getTrueDCAfromPropagatedTrack(AliESDtrack* track ,Double_t* newdca,Double_t*newposatdca)
{
    Double_t alpha= track->GetAlpha();
    Double_t y = track->GetX() *  TMath::Sin(alpha) +  track->GetY()  * TMath::Cos(alpha);
    Double_t x =track->GetX() *  TMath::Cos(alpha) +   track->GetY()  * TMath::Sin(-alpha);

    Double_t v[3] = {0x0,0x0,0x0};
    const AliESDVertex *vtx =    (const AliESDVertex *)InputEvent()->GetPrimaryVertex();
    vtx->GetXYZ(v);
    Double_t z =   track->GetZ();
    newposatdca[0] = x;
    newposatdca[1] = y;
    newposatdca[2] = z;
    Double_t dcaxy = sqrt((x-v[0])*(x-v[0]) + (y-v[1])*(y-v[1]));
    Double_t dcaxyz = sqrt((x-v[0])*(x-v[0]) + (y-v[1])*(y-v[1])+(z-v[2])*(z-v[2]) );
    newdca[0] = dcaxy;
    newdca[1] = dcaxyz;

}

Bool_t AliAnalysisTaskHFJetIPQA::Run(){
    if(fESD) fIsEsd = kTRUE;
    this->fDisableWeightingMC=kFALSE;
    fEtaBEvt.clear(); fPhiBEvt.clear(); fEtaCEvt.clear();fPhiCEvt.clear();fEtaUdsgEvt.clear(); fPhiUdsgEvt.clear();
    fEtaSEvt.clear();fPhiSEvt.clear(); fMcEvtSampled = kFALSE; fMCArray= NULL;fMCEvent = NULL;
    AliESDtrack* trackV = NULL;

    Float_t dca[2] = {-99999,-99999};
    Float_t cov[3] = {-99999,-99999,-99999};

    Double_t weight       = 1;
    Int_t nTracksInEvent  = 0 ;
    Bool_t hasIPSuccess = kFALSE;

    if(fIsPythia)     fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    nTracksInEvent                  = (Int_t)InputEvent()->GetNumberOfTracks() ;
    if(nTracksInEvent<1) return kTRUE;

    for(long itrack= 0; itrack<nTracksInEvent;++itrack)
        {


            trackV=(AliESDtrack*)(InputEvent()->GetTrack(itrack));
            const AliVVertex * vtx= InputEvent()->GetPrimaryVertex();
            if(!trackV) continue;

            IncHist("fh1dTracksAccepeted",1);

            if(!IsTrackAccepted((AliESDtrack*)trackV,6)) {
                    IncHist("fh1dTracksAccepeted",3);
                    continue;
                }

            IncHist("fh1dTracksAccepeted",2);
            FillHist("fh2dAcceptedTracksEtaPhi",trackV->Eta(),trackV->Phi(),1);
            weight =1;dca[0]=-9999;dca[1]=-9999;cov[0]=-9999;cov[1]=-9999;cov[2]=-9999;
            hasIPSuccess =kFALSE;
            if (CalculateTrackImpactParameter((AliESDtrack*)trackV,dca,cov))hasIPSuccess =kTRUE;
            if(!hasIPSuccess)  continue;
            else {
                    weight =1;
                    Double_t posdcatrack[3]= {0.,0.,0.};
                    Double_t truetrackvtx[3]= {0.,0.,0.};

                    FillHist("fh2dTracksImpParXY",GetValImpactParameter(kXY,dca,cov)*-1,trackV->Pt(),weight);
                    FillHist("fh2dTracksImpParXY",GetValImpactParameter(kXY,dca,cov)*1,trackV->Pt(),weight);
                    FillHist("fh2dTracksImpParZ",dca[1],trackV->Pt(),weight);
                    FillHist("fh2dTracksImpParXYSignificance",GetValImpactParameter(kXYSig,dca,cov),trackV->Pt(),weight);
                    FillHist("fh2dTracksImpParZSignificance",GetValImpactParameter(kZSig,dca,cov),trackV->Pt(),weight);
                    FillHist("fh2dTracksImpParXYZ",GetValImpactParameter(kXYZ,dca,cov),trackV->Pt(),weight);
                    FillHist("fh2dTracksImpParXYZSignificance",GetValImpactParameter(kXYZSig,dca,cov),trackV->Pt(),weight);
                    FillHist("fh1dTracksImpParXY",GetValImpactParameter(kXY,dca,cov),weight);
                    FillHist("fh1dTracksImpParXYZ",GetValImpactParameter(kXYZ,dca,cov),weight);
                    FillHist("fh1dTracksImpParXYSignificance",GetValImpactParameter(kXYSig,dca,cov),weight);
                    FillHist("fh1dTracksImpParXYZSignificance",GetValImpactParameter(kXYZSig,dca,cov),weight);


                    if(fIsPythia){
                            Int_t corrpartidx =-1;
                            weight = GetMonteCarloCorrectionFactor(trackV,corrpartidx);

                            //Run geometric track reconstruction efficiency correction
                            FillHist("fh2dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov)*-1.,trackV->Pt(),weight);
                            FillHist("fh2dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov)*1.,trackV->Pt(),weight);

                            FillHist("fh1dTracksImpParXY_McCorr",GetValImpactParameter(kXY,dca,cov),weight);
                            FillHist("fh1dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),weight);
                            FillHist("fh1dTracksImpParXYSignificance_McCorr",GetValImpactParameter(kXYSig,dca,cov),weight);
                            FillHist("fh1dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),weight);
                            FillHist("fh2dTracksImpParXYZ_McCorr",GetValImpactParameter(kXYZ,dca,cov),trackV->Pt(),weight);
                            FillHist("fh2dTracksImpParXYZSignificance_McCorr",GetValImpactParameter(kXYZSig,dca,cov),trackV->Pt(),weight);

                        }
                }
        }

    if(fCorrrectionSamplingMode) return kTRUE;


    AliJetContainer *  jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    if (!jetconrec) return kFALSE;
    AliJetContainer * jetcongen = 0x0;


    AliEmcalJet * jetgen  = 0x0;
    if(fIsPythia){
            jetcongen = static_cast<AliJetContainer*>(fJetCollArray.At(1));
            if(!MatchJetsGeometricDefault()) cout << "Error running jet matching!" << endl;
            jetcongen->ResetCurrentID();
            // Fill gen. level jet histograms
            while ((jetgen = jetcongen->GetNextJet()))
                {
                    if (!jetgen) continue;
                    Int_t jetflavour =0;
                    Bool_t is_udgjet = kFALSE;
                    jetflavour =IsMCJetPartonFast(jetgen,0.4,is_udgjet); //Event based association to save memory
                    FillHist("fh1dJetGenPt",GetPtCorrectedMC(jetgen),1);
                    if(jetflavour ==0)
                        FillHist("fh1dJetGenPtUnidentified",GetPtCorrectedMC(jetgen),1);
                    else if(jetflavour ==1)
                        FillHist("fh1dJetGenPtudsg",GetPtCorrectedMC(jetgen),1);
                    else if(jetflavour ==2)
                        FillHist("fh1dJetGenPtc",GetPtCorrectedMC(jetgen),1);
                    else if(jetflavour ==3)
                        FillHist("fh1dJetGenPtb",GetPtCorrectedMC(jetgen),1);
                }


            jetcongen->ResetCurrentID();
            jetconrec->ResetCurrentID();
        }
    // loop rec level jets
    AliEmcalJet * jetrec  = 0x0;
    AliEmcalJet * jetmatched  = 0x0;
    jetconrec->ResetCurrentID();
    Double_t jetpt=0;

    while ((jetrec = jetconrec->GetNextJet()))
        {
            if(!jetrec) break;
            jetpt = jetrec->Pt();
            if(!(jetconrec->GetRhoParameter() == 0x0)){
                    jetpt = jetpt - jetconrec->GetRhoVal() * jetrec->Area();
                }
            if(fIsPythia){
                    if (jetrec->MatchedJet()) {
                            Double_t genpt = jetrec->MatchedJet()->Pt();
                            if(!(jetcongen->GetRhoParameter() == 0x0)){
                                    genpt = genpt - jetcongen->GetRhoVal() * jetrec->MatchedJet()->Area();
                                }
                            FillHist("fh2dJetGenPtVsJetRecPt",genpt,jetpt,1);
                        }
                }


            // make inclusive signed imp. parameter constituent histograms
            Int_t ntracks = (Int_t)jetrec->GetNumberOfTracks();
            Float_t dca[2] = {-99999,-99999};
            Float_t cov[3] = {-99999,-99999,-99999};
            Double_t sign=0;
            Int_t jetflavour=0;
            Bool_t is_udgjet = kFALSE;
            if(fIsPythia){
                    jetmatched = 0x0;
                    jetmatched =jetrec->MatchedJet();
                    if(jetmatched && fESD)jetflavour =
                            IsMCJetPartonFast(jetmatched,0.4,is_udgjet); //Event based association to save memory



                }
            FillHist("fh1dJetRecPt",jetrec->Pt(),1);
            if(fIsPythia){
                    if(jetflavour==0)     FillHist("fh1dJetRecPtUnidentified",jetpt,1);
                    else if(jetflavour==1)FillHist("fh1dJetRecPtudsg",        jetpt,1);
                    else if(jetflavour==2)FillHist("fh1dJetRecPtc",           jetpt,1);
                    else if(jetflavour==3)FillHist("fh1dJetRecPtb",           jetpt,1);
                }

            if(!(fJetCutsHF->IsJetSelected(jetrec))) break;
            FillHist("fh1dJetRecEtaPhiAccepted",jetrec->Eta(),jetrec->Phi(),1);
            FillHist("fh1dJetRecPtAccepted",jetpt,1);
            if(fIsPythia){
                    if(jetflavour==0)     FillHist("fh1dJetRecPtUnidentifiedAccepted",jetpt,1);
                    else if(jetflavour==1)FillHist("fh1dJetRecPtudsgAccepted",        jetpt,1);
                    else if(jetflavour==2)FillHist("fh1dJetRecPtcAccepted",           jetpt,1);
                    else if(jetflavour==3)FillHist("fh1dJetRecPtbAccepted",           jetpt,1);
                }


            std::vector<SJetIpPati> sImpParXY,sImpParXYZ,sImpParXYSig,sImpParXYZSig;
            for(Int_t itrack = 0; itrack < InputEvent()->GetNumberOfTracks(); ++itrack)
                {
                    Bool_t hasSIP=kFALSE;
                    Double_t dcatrackjet =999,lineardecaylenth = 999;
                    AliESDtrack * trackV = (AliESDtrack *) InputEvent()->GetTrack(itrack);
                    if (!trackV || !jetrec)            continue;
                    if (jetrec->DeltaR(trackV) > 0.4) continue;
                    if (!IsTrackAccepted((AliESDtrack*)trackV,6))   continue;



                    if(CalculateJetSignedTrackImpactParameter((AliESDtrack*)trackV,jetrec,dca,cov,sign,dcatrackjet,lineardecaylenth))hasSIP =kTRUE;


                    if(hasSIP)
                        {






                            if(lineardecaylenth>0.5) continue;
                            if(dcatrackjet > 0.07) continue;

                            Int_t corridx=-1;
                            Double_t cursImParXY    =TMath::Abs(GetValImpactParameter(   kXY,dca,cov))*sign;
                            Double_t cursImParXYSig =TMath::Abs(GetValImpactParameter(kXYSig,dca,cov))*sign;
                            Double_t cursImParXYZ    =TMath::Abs(GetValImpactParameter(   kXYZ,dca,cov))*sign;
                            Double_t cursImParXYZSig =TMath::Abs(GetValImpactParameter(kXYZSig,dca,cov))*sign;

                            fIsPythia ? weight = GetMonteCarloCorrectionFactor(trackV,corridx) : weight =1;
                            //FILL Resolution function for non heavy/strange jets
                              if(is_udgjet){
                                 FillHist("fh2dJetSignedImpParXYSignificanceudg_6",jetpt,cursImParXYSig,weight);
                                 FillHist("fh2dJetSignedImpParXYZSignificanceudg_6",jetpt,cursImParXYZSig,weight);
                              }



                            FillHist("fh2dJetSignedImpParXY"            ,jetpt,cursImParXY,weight);
                            FillHist("fh2dJetSignedImpParXYSignificance",jetpt,cursImParXYSig,weight);
                            const char * subtype [4] = {"Unidentified","udsg","c","b"};
                            if(fIsPythia){
                                    FillHist(Form("fh2dJetSignedImpParXY%s",subtype[jetflavour]),jetpt,cursImParXY,weight);
                                    FillHist(Form("fh2dJetSignedImpParXYSignificance%s",subtype[jetflavour]),jetpt,cursImParXYSig,weight);
                                }
                            SJetIpPati a(cursImParXY, weight,kFALSE,kFALSE); sImpParXY.push_back(a);
                            SJetIpPati b(cursImParXYZ, weight,kFALSE,kFALSE); sImpParXYZ.push_back(b);
                            SJetIpPati c(cursImParXYSig, weight,kFALSE,kFALSE);sImpParXYSig.push_back(c);
                            SJetIpPati d(cursImParXYZSig, weight,kFALSE,kFALSE);sImpParXYZSig.push_back(d);

                        }
                }
            std::sort(sImpParXY.begin(),sImpParXY.end(),       AliAnalysisTaskHFJetIPQA::mysort);
            std::sort(sImpParXYSig.begin(),sImpParXYSig.end(), AliAnalysisTaskHFJetIPQA::mysort);
            std::sort(sImpParXYZ.begin(),sImpParXYZ.end(), AliAnalysisTaskHFJetIPQA::mysort);
            std::sort(sImpParXYZSig.begin(),sImpParXYZSig.end(), AliAnalysisTaskHFJetIPQA::mysort);
            const char * subtype[4] = {"Unidentified","udsg","c","b"};
            const char * subord [3] = {"First","Second","Third"};
            const char * stype  [4] = {"fh2dJetSignedImpParXY","fh2dJetSignedImpParXYSignificance","fh2dJetSignedImpParXYZ","fh2dJetSignedImpParXYZSignificance"};
            for (Int_t ot = 0 ; ot <3 ;++ot){
                    if ((int)sImpParXY.size()>ot){
                            Double_t params [4] ={sImpParXY.at(ot).first,sImpParXYSig.at(ot).first,sImpParXYZ.at(ot).first,sImpParXYZSig.at(ot).first};
                            Double_t weights[4] ={sImpParXY.at(ot).second,sImpParXYSig.at(ot).second,sImpParXYZ.at(ot).second,sImpParXYZSig.at(ot).second};
                            //non flavour dependent histograms
                            for (Int_t ost = 0 ; ost <4 ;++ost){
                                    TString hname = Form("%s%s",stype[ost],subord[ot]);
                                    if(fIsPythia&& !fDisableWeightingMC)   FillHist(hname.Data(),jetpt,params[ost],weights[ost] );
                                    else  FillHist(hname.Data(),jetpt,params[ost],1);
                                }
                            if(fIsPythia){
                                    for (Int_t ost = 0 ; ost <4 ;++ost){
                                            TString hname = Form("%s%s%s",stype[ost],subtype[jetflavour],subord[ot]);
                                            FillHist(hname.Data(),jetpt,params[ost],weights[ost] );
                                            TString hname2 = Form("%s%s%suw",stype[ost],subtype[jetflavour],subord[ot]);
                                            FillHist(hname2.Data(),jetpt,params[ost],1);
                                        }
                                }
                        }
                }
            sImpParXY.clear();
            sImpParXYSig.clear();
        }
    return kTRUE;
}
//###############################################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSelected(AliVEvent *event, Int_t &WhyRejected,ULong_t &RejectionBits){
    WhyRejected =0;
    Bool_t accept=kTRUE;
    Double_t fMaxVtxZ = 10;
    RejectionBits=000;
    Bool_t isSelected = kFALSE;
    isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    if(!isSelected) {
            WhyRejected=kPhysicsSelection;
            return kFALSE;
        }
    const AliVVertex* spdVtx = dynamic_cast<const AliVVertex*>(((AliESDEvent*)event)->GetPrimaryVertex());
    if(!spdVtx) {
            WhyRejected=kNoVertex;
            return kFALSE;
        }
    if(spdVtx->GetNContributors()<1) {
            WhyRejected=kNoContributors;
            return kFALSE;
        }
    if(spdVtx->GetNContributors()<4) {
            WhyRejected=kTooFewVtxContrib;
            return kFALSE;
        }
    if(TMath::Abs(spdVtx->GetZ())>=fMaxVtxZ) {
            WhyRejected=kZVtxOutFid;
            return kFALSE;
        }
    if((TMath::Abs(spdVtx->GetX() - ((AliESDEvent*)event)->GetDiamondX()) > 2.*sqrt(((AliESDEvent*)event)->GetSigma2DiamondX()))) {
            WhyRejected=kBadDiamondXDistance;
            return kFALSE;
        }
    if((TMath::Abs(spdVtx->GetY() - ((AliESDEvent*)event)->GetDiamondY()) > 2.*sqrt(((AliESDEvent*)event)->GetSigma2DiamondY()))) {
            WhyRejected=kBadDiamondYDistance;
            return kFALSE;
        }
    if((TMath::Abs(spdVtx->GetZ() - ((AliESDEvent*)event)->GetDiamondZ()) > 2.*sqrt(((AliESDEvent*)event)->GetSigma2DiamondZ()))) {
            WhyRejected=kBadDiamondZDistance;
            return kFALSE;
        }

    Double_t x2ndf = (((AliESDEvent*)event)->GetPrimaryVertex())->GetChi2perNDF();

    return accept;
}
//###############################################################################################################
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
//###############################################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::SetResFunction(TGraph *f, Int_t j){
    fResolutionFunction[j] = *f;
    return kTRUE;
}
// ######################################################################################## Event Selection
Bool_t AliAnalysisTaskHFJetIPQA::IsEventSelected()	{
    if(fIsPythia){
            AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
            if(!mcH ){
                    AliError("No MC Event Handler available");
                    return kFALSE;
                }
            if(!mcH->InitOk()) return kFALSE;
            if(!mcH->TreeK()) return kFALSE ;
            if(!mcH->TreeTR()) return kFALSE;
        }

    AliAODEvent* aev = NULL;
    Int_t WhyRejected =0;
    ULong_t RejectionBits=0;
    if(!fESD)
        {
            aev = dynamic_cast<AliAODEvent*>(InputEvent());
            if(aev && aev->GetPrimaryVertex() && aev->GetPrimaryVertex()->GetNContributors()>0){
                    FillHist("fh1dVertexZ",aev->GetPrimaryVertex()->GetZ(),1);
                    Double_t vtxx = aev->GetPrimaryVertex()->GetX();
                    Double_t vtxy = aev->GetPrimaryVertex()->GetY();
                    FillHist("fh1dVertexR",vtxx,vtxy,1);
                }else return kFALSE;


            if(!(IsSelected(aev,WhyRejected,RejectionBits)))
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
                    IncHist("fh1dEventRejectionRDHFCuts",1);
                    FillHist("fh1dVertexZAccepted",aev->GetPrimaryVertex()->GetZ(),1);
                    Double_t vtxx = aev->GetPrimaryVertex()->GetX();
                    Double_t vtxy = aev->GetPrimaryVertex()->GetY();
                    FillHist("fh1dVertexRAccepted",vtxx,vtxy,1);
                    FillHist("fh2dVertexChi2NDFNESDTracks",aev->GetPrimaryVertex()->GetChi2perNDF(),aev->GetNumberOfESDTracks(),1);
                    return kTRUE;
                }
        }
    AliESDEvent* eev = NULL;
    if(fESD)
        {
            eev = dynamic_cast<AliESDEvent*>(InputEvent());
            if(eev && eev->GetPrimaryVertex() && eev->GetPrimaryVertex()->GetNContributors()>0){
                    FillHist("fh1dVertexZ",eev->GetPrimaryVertex()->GetZ(),1);
                    Double_t vtxx = eev->GetPrimaryVertex()->GetX();
                    Double_t vtxy = eev->GetPrimaryVertex()->GetY();
                    FillHist("fh1dVertexR",vtxx,vtxy,1);
                }else {
                    return kFALSE;
                }

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
                    Double_t vtxx = eev->GetPrimaryVertex()->GetX();
                    Double_t vtxy = eev->GetPrimaryVertex()->GetY();
                    Double_t vtxz = eev->GetPrimaryVertex()->GetZ();
                    FillHist("fh1dVertexXvsMultiplicity",vtxx,eev->GetNumberOfTracks(),1);
                    FillHist("fh1dVertexYvsMultiplicity",vtxy,eev->GetNumberOfTracks(),1);
                    FillHist("fh1dVertexZvsMultiplicity",vtxz,eev->GetNumberOfTracks(),1);
                    IncHist("fh1dEventRejectionRDHFCuts",1);
                    FillHist("fh1dVertexZAccepted",eev->GetPrimaryVertex()->GetZ(),1);
                    FillHist("fh1dVertexRAccepted",vtxx,vtxy,1);
                    FillHist("fh2dVertexChi2NDFNESDTracks",eev->GetPrimaryVertex()->GetChi2perNDF(),eev->GetNumberOfTracks(),1);
                    return kTRUE;
                }
        }
    return kFALSE;
}
// ######################################################################################## Init histograms
void AliAnalysisTaskHFJetIPQA::UserCreateOutputObjects(){

    if(!fRandom) fRandom = new TRandom3(0);

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (mgr) {
            AliVEventHandler *evhand = mgr->GetInputEventHandler();
            if (evhand) {
                    if (evhand->InheritsFrom("AliESDInputHandler")) {
                            fIsEsd = kTRUE;
                        }
                    else {
                            fIsEsd = kFALSE;
                        }
                }
            else {
                    AliError("Event handler not found!");
                }
        }
    else {
            AliError("Analysis manager not found!");
        }
    OpenFile(1);



    // AliAnalysisTaskEmcal::UserCreateOutputObjects();
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

    const Int_t nBins3dSignificance =500;
    const Int_t nBins3d =250;

    Double_t lowIPxy =-1.;
    Double_t highIPxy =1.;
    OpenFile(1);
    if (!fOutput2) fOutput2 = new TList ();
    fOutput2->SetOwner(kTRUE);

    //Make Graphs
    const Int_t gfProtonN = 9;
    Double_t gfProtonX [gfProtonN] ={0,0.534483,1.29741,2.21552,3.0819,3.92241,4.5819,5.39655,1000};
    Double_t gfProtonY [gfProtonN] ={0.990964,0.990964,0.990964,0.990964,0.990964,0.990964,0.990964,0.990964,0.990964};

    const Int_t gfAntiProtonN = 18;
    Double_t gfAntiProtonX [gfAntiProtonN] ={0,0.806034,0.922414,1.09052,1.28448,1.5431,1.73707,1.89224,2.17672,2.43534,2.74569,3.06897,
                                             3.52155,3.88362,4.38793,5.03448,5.38362, 1000};
    Double_t gfAntiProtonY [gfAntiProtonN] ={0.922892,0.922892,	0.930723,	0.939157,0.94397,0.95241,0.956627,0.959639,0.964458,
                                             0.966867,0.971084,0.974096,0.978313,0.98012,0.983735,0.986747,0.989157,0.989157};


    const Int_t gfAntiLambdaN = 34;
    Double_t gfAntiLambdaX [gfAntiLambdaN] ={0.,0.55555,0.64646,0.75757,	0.84848,0.94949,1.06061,1.15152,1.24242,1.35354,1.44444,
                                             1.54545,1.66667,1.75758,1.84848,1.9596,2.09091,2.30303,2.50505,2.68687,2.90909,3.11111,
                                             3.31313,3.51515,3.69697,3.89899,4.20202,4.66667,5.21212,5.74747,6.50505,7.51515,9.0101,1000};
    Double_t gfAntiLambdaY [gfAntiLambdaN] = {0.864925,0.864925,0.895896,0.908209,0.915672,0.921269,0.926866,0.931343,0.935821,0.938806,0.942164,
                                              0.945149,0.947761,0.95,0.952612,0.954478,0.957836,0.960821,0.96306,0.965672,0.968657,0.970149,
                                              0.972015,0.973507,0.975,0.976493,0.978358,0.981343,0.983955,0.986194,0.988433,0.991045,0.991045,0.991045};

    const Int_t gfLambdaN = 2;
    Double_t gfLambdaX [gfLambdaN] =	{0.,1000};
    Double_t gfLambdaY [gfLambdaN] = {0.991045,0.991045};
    const Int_t gfKMinusN =13 ;
    Double_t gfKMinusX [gfKMinusN] =	{0,0.54741,0.74137,1.03879,1.36207,1.96983,2.52586,3.0819,3.67672,4.19397,5.03448,5.44828,1000};
    Double_t gfKMinusY [gfKMinusN] = {0,0.979518,0.983133,0.987349,0.989759,0.992169,0.993976,0.996386,0.995783,0.998193,0.99759,1,1000};

    fGeant3FlukaProton 	 = new TGraph(gfProtonN,gfProtonX,gfProtonY);
    fGeant3FlukaAntiProton = new TGraph(gfAntiProtonN,gfAntiProtonX,gfAntiProtonY);
    fGeant3FlukaLambda   	 = new TGraph(gfLambdaN,gfLambdaX,gfLambdaY);
    fGeant3FlukaAntiLambda = new TGraph(gfAntiLambdaN,gfAntiLambdaX,gfAntiLambdaY);
    fGeant3FlukaKMinus 	 = new TGraph(gfKMinusN,gfKMinusX,gfKMinusY);


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




    AddHistogramm("fh1dVertexXvsMultiplicity",";cm;# ESD Tracks",300,-1,1,200,0,200);
    AddHistogramm("fh1dVertexYvsMultiplicity",";cm;# ESD Tracks",300,-1,1,200,0,200);
    AddHistogramm("fh1dVertexZvsMultiplicity",";cm;# ESD Tracks",300,-10,10 ,200,0,200);

    //Vertex Z before and after
    AddHistogramm("fh1dVertexZ","Vertex Z before Event selection;primary vertex z (cm);count",500,-30,30);
    AddHistogramm("fh1dVertexZAccepted","Vertex Z after Event selection;primary vertex z (cm);count",500,-30,30);

    AddHistogramm("fh1dVertexR","Vertex R before Event selection;primary vertex xy (cm);x;y",500,-0.5,0.5,500,-0.5,0.5);
    AddHistogramm("fh1dVertexRAccepted","Vertex R after Event selection;primary vertex xy (cm);x;y",500,-0.5,0.5,500,-0.5,0.5);

    AddHistogramm("fh2dVertexChi2NDFNESDTracks","Vertex Chi2/NDF vs # tracks ESD;vertex #chi^{2}/NDF;# tracks esd",200,0,10,500,0,500);
    // AOD tracks accepted
    AddHistogramm("fh1dTracksAccepeted","# tracks before/after cuts;;",3,0,3);
    TH1D * h1 = GetHist1D("fh1dTracksAccepeted");
    h1->GetXaxis()->SetBinLabel(1,"total");
    h1->GetXaxis()->SetBinLabel(2,"accepted");
    h1->GetXaxis()->SetBinLabel(3,"rejected");
    // Tracks impact parameter histograms
    AddHistogramm("fh2dTracksImpParXY","radial imp. parameter ;impact parameter xy (cm);a.u.",2000,lowIPxy,highIPxy,500,0,100.);
    AddHistogramm("fh2dTracksImpParXYZ","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-1,1,500,0,100.);
    AddHistogramm("fh2dTracksImpParXYZSignificance","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.);
    AddHistogramm("fh2dTracksImpParZ","z imp. parameter ;impact parameter xy (cm);a.u.",2000,lowIPxy,highIPxy,500,0,10.);
    AddHistogramm("fh2dTracksImpParXYSignificance","radial imp. parameter sig;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.);
    AddHistogramm("fh2dTracksImpParZSignificance","z imp. parameter ;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.);
    AddHistogramm("fh1dTracksImpParXY","2d imp. parameter ;impact parameter 2d (cm);a.u.",10001,-1,1.);
    AddHistogramm("fh1dTracksImpParXYZ","3d imp. parameter ;impact parameter 3d (cm);a.u.",2000,0,1.);
    AddHistogramm("fh1dTracksImpParXYSignificance","radial imp. parameter ;impact parameter xy significance;a.u.",4000,-30,30);
    AddHistogramm ("fh1dTracksImpParXYZSignificance","3d imp. parameter ;impact parameter 3d significance;a.u.",2000,0.,100.);
    AddHistogramm("fh1dJetRecEtaPhiAccepted","detector level jet;#eta;phi",200,-0.5,0.5,200,0.,TMath::TwoPi());
    AddHistogramm("fh2dAcceptedTracksEtaPhi","accepted tracks;#eta;phi",200,-0.9,0.9,200,0.,TMath::TwoPi());

    for(Int_t i = 0; i<6;++i){
            fh2dAcceptedTracksEtaPhiPerLayer[i] = new TH2D(Form("fh2dAcceptedTracksEtaPhiPerLayer_%i",i),"accepted tracks;#eta;phi",200,-0.9,0.9,200,0.,TMath::TwoPi());
        }

    AddHistogramm("fh1dJetRecPt","detector level jets;pt (GeV/c); count",500,0,250);
    AddHistogramm("fh1dJetRecPtAccepted","accepted detector level jets;pt (GeV/c); count",500,0,250);




    if (fIsPythia){


            AddHistogramm("fh2dJetSignedImpParXYSignificanceudg_6","fh2dJetSignedImpParXYSignificanceudg;pt (GeV/c); count",200,0,100,2000,-100,100);
            AddHistogramm("fh2dJetSignedImpParXYZSignificanceudg_6","fh2dJetSignedImpParXYZSignificanceudg;pt (GeV/c); count",200,0,100,2000,-100,100);



            AddHistogramm("fh2dTracksImpParXY_McCorr","radial imp. parameter (after correction);impact parameter xy (cm);a.u.",2000,-1,1,500,0,100);
            AddHistogramm("fh1dTracksImpParXY_McCorr","radial imp. parameter (after correction);impact parameter xy (cm);a.u.",10000,-1,1.);
            AddHistogramm("fh1dTracksImpParXYZ_McCorr","3d imp. parameter (after correction);impact parameter 3d (cm);a.u.",2000,0,100.);
            AddHistogramm("fh2dTracksImpParXYZ_McCorr","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-1,1,500,0,100.);
            AddHistogramm("fh2dTracksImpParXYZSignificance_McCorr","XYZ imp. parameter ;impact parameter xy (cm);a.u.",2000,-30,30,500,0,100.);

            AddHistogramm("fh1dTracksImpParXYSignificance_McCorr","radial imp. parameter (after correction);impact parameter xy significance;a.u.",10000,-30,30.);
            AddHistogramm("fh1dTracksImpParXYZSignificance_McCorr","3d imp. parameter (after correction);impact parameter 3d significance;a.u.",2000,0.,100.);
            AddHistogramm("fh1dJetGenPt","generator level jets;pt (GeV/c); count",500,0,250);
            AddHistogramm("fh1dJetGenPtUnidentified","generator level jets (no flavour assigned);pt (GeV/c); count",500,0,250);
            AddHistogramm("fh1dJetGenPtudsg","generator level udsg jets;pt (GeV/c); count",500,0,250);
            AddHistogramm("fh1dJetGenPtc","generator level c jets;pt (GeV/c); count",500,0,250);
            AddHistogramm("fh1dJetGenPtb","generator level b jets;pt (GeV/c); count",500,0,250);
            AddHistogramm("fh2dJetGenPtVsJetRecPt","detector momentum response;gen pt;rec pt",500,0,250,500,0,250);

            AddHistogramm("fh1dJetRecPtudsg","detector level jets;pt (GeV/c); count",500,0,250);
            AddHistogramm("fh1dJetRecPtUnidentified","detector level jets;pt (GeV/c); count",500,0,250);

            AddHistogramm("fh1dJetRecPtc","detector level jets;pt (GeV/c); count",500,0,250);
            AddHistogramm("fh1dJetRecPtb","detector level jets;pt (GeV/c); count",500,0,250);
            AddHistogramm("fh1dJetRecPtUnidentifiedAccepted","detector level jets;pt (GeV/c); count",500,0,250);
            AddHistogramm("fh1dJetRecPtudsgAccepted","detector level jets;pt (GeV/c); count",500,0,250);
            AddHistogramm("fh1dJetRecPtcAccepted","detector level jets;pt (GeV/c); count",500,0,250);
            AddHistogramm("fh1dJetRecPtbAccepted","detector level jets;pt (GeV/c); count",500,0,250);
        }

    const char * flavour[5]  = {"Unidentified","udsg","c","b",""};

    for (Int_t i = 0 ;i < 5;++i){
            if (!fIsPythia && (i<4)) continue;
            AddHistogramm(Form("JetProbability_%s",flavour[i]),Form("JetProbability_%s;jept",flavour[i]),500,0,250,500,0,1);
        }
    const char * base = "fh2dJetSignedImpPar";
    const char * dim[2]  = {"XY","XYZ"};
    const char * typ[2]  = {"","Significance"};
    const char * ordpar [4] = {"","First","Second","Third"};
    const char * special [1] = {"",/*"McCorr"*/};

    Int_t ptbins = 250;
    Double_t ptlow = 0;
    Double_t pthigh = 250;

    Int_t ipbins = 1000;
    Double_t iplow = -.3;
    Double_t iphigh = .3;

    for (Int_t id = 0;id<2;++id)
        for (Int_t ifl = 0;ifl<5;++ifl)
            for (Int_t io = 0;io<4;++io)
                for (Int_t is = 0;is<1;++is)
                    for (Int_t it = 0;it<2;++it){
                            if(it==1) {
                                    iplow=-100;
                                    iphigh=100; //from 30
                                    if(io==0 && ifl==4) ipbins =2000;
                                    else  ipbins =2000;
                                }else {
                                    iplow=-1.;
                                    iphigh=1.;
                                    ipbins =2000;

                                }
                            //High-res 2d templates/  low-res 3d
                            if(id==0)  ipbins =2000;
                            if(id==0)  ipbins =2000;
                            if((fIsPythia||(!fIsPythia && ifl==4)))  AddHistogramm(Form("%s%s%s%s%s%s",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]),
                                                                                                                    Form("%s%s%s%s%s%s;;",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]),
                                                                                                                    ptbins,ptlow,pthigh,ipbins,iplow,iphigh);
                            if((fIsPythia||(!fIsPythia && ifl==4)))  AddHistogramm(Form("%s%s%s%s%s%suw",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]),
                                                                                                                    Form("%s%s%s%s%s%suw;;",base,dim[id],typ[it],flavour[ifl],ordpar[io],special[is]),
                                                                                                                    ptbins,ptlow,pthigh,1,iplow,iphigh);
                        }
    PostData(1, fOutput2); // Post data for ALL output slots > 0 here.

}

Bool_t AliAnalysisTaskHFJetIPQA::CalculateTrackImpactParameter(AliESDtrack * track,Float_t *impar, Float_t * cov,Bool_t useTRUEvtx)
{
    Double_t dv[2],dcov[3];
    const Double_t kBeampiperadius=3;
    AliESDEvent * eev = (AliESDEvent*)InputEvent();
    const  AliVVertex *vtxESDSkip =(const  AliVVertex *) (InputEvent()->GetPrimaryVertex())  ;
    if(!vtxESDSkip) return kFALSE;
    if(track){
            if(track->PropagateToDCA(vtxESDSkip, eev->GetMagneticField(), kBeampiperadius, dv, dcov)){
                    Double_t xyz_at_dca_track[3]= {0.,0.,0.};
                    track->GetXYZ(xyz_at_dca_track );
                    Double_t xyz_vertex[3]= {0.,0.,0.};
                    vtxESDSkip->GetXYZ(xyz_vertex);
                    impar[0] = (Float_t) (TMath::Abs(dv[0]) );
                    // if(fIsPythia)impar[0] *=1.036;
                    impar[1] = (Float_t) (TMath::Abs(sqrt((impar[0] *impar[0] + (xyz_at_dca_track[2]-xyz_vertex[2])*(xyz_at_dca_track[2]-xyz_vertex[2])))));
                    cov[0] = (Float_t) dcov[0];
                    cov[1] = (Float_t) dcov[1];
                    cov[2] = (Float_t) dcov[2];
                    return kTRUE;
                }
        }
    return kFALSE;
}

Bool_t AliAnalysisTaskHFJetIPQA::CalculateJetSignedTrackImpactParameter(AliESDtrack * track,AliEmcalJet * jet ,Float_t *impar, Float_t * cov, Double_t &sign, Double_t &dcajetrack, Double_t &lineardecaylength){
    AliESDEvent* esdevent = dynamic_cast<AliESDEvent*>(InputEvent());
    const Double_t kBeampiperadius=3;
    const AliESDVertex *vtxESDSkip =    (const AliESDVertex *)esdevent->GetPrimaryVertex();
    if(track){

            Double_t dv[2]={0},dcov[3]={0},bcv[21]={0x0};
            if(!(track->PropagateToDCA(vtxESDSkip, esdevent->GetMagneticField(), kBeampiperadius, dv, dcov))) return kFALSE;


            Double_t xyz_at_dca_track[3]= {0.,0.,0.};
            impar[0] = (Float_t) (TMath::Abs(dv[0]));
            // if(fIsPythia)  impar[0] *=1.036;
            track->GetXYZ(xyz_at_dca_track );
            cov[0] = (Float_t) dcov[0];
            cov[1] = (Float_t) dcov[1];
            cov[2] = (Float_t) dcov[2];

            Double_t xyz_vertex[3]= {0.,0.,0.};
            vtxESDSkip->GetXYZ(xyz_vertex);
            Double_t ipvector3a[3] = { xyz_at_dca_track[0] - xyz_vertex[0], xyz_at_dca_track[1] - xyz_vertex[1], xyz_at_dca_track[2] - xyz_vertex[2]};
            impar[1] = (Float_t) (TMath::Abs(sqrt((impar[0] *impar[0] + (xyz_at_dca_track[2]-xyz_vertex[2])*(xyz_at_dca_track[2]-xyz_vertex[2])))));

            Double_t bpxpypz[3] = { jet->Px(), jet->Py(), jet->Pz() };


            AliExternalTrackParam bjetparam(xyz_vertex, bpxpypz, bcv, (Short_t)0);
            Double_t xa,xb,xyz_jet_global[3],xyz_track_global[3];
            bjetparam.GetDCA(track, esdevent->GetMagneticField(), xa, xb);
            bjetparam.GetXYZAt(xa, esdevent->GetMagneticField(),xyz_jet_global);
            track->   GetXYZAt(xb, esdevent->GetMagneticField(),xyz_track_global);
            track->PropagateTo(xb,esdevent->GetMagneticField());
            track->PropagateToDCA(vtxESDSkip, 0, kBeampiperadius, dv, dcov);
            track->GetXYZ(xyz_at_dca_track);


            Double_t ipvector3[3] = { xyz_at_dca_track[0] - xyz_vertex[0], xyz_at_dca_track[1] - xyz_vertex[1], xyz_at_dca_track[2] - xyz_vertex[2]};
            sign = TMath::Sign(1.,ipvector3[0]*jet->Px() +ipvector3[1]*jet->Py()+ipvector3[2]*jet->Pz());


            Double_t  bdecaylength =  TMath::Sqrt((xyz_at_dca_track[0] - xyz_track_global[0]) * (xyz_at_dca_track[0] - xyz_track_global[0]) +
                    (xyz_at_dca_track[1] - xyz_track_global[1]) * (xyz_at_dca_track[1] - xyz_track_global[1])+
                    (xyz_at_dca_track[2] - xyz_track_global[2]) * (xyz_at_dca_track[2] - xyz_track_global[2]));
            dcajetrack = TMath::Sqrt((xyz_jet_global[0] - xyz_track_global[0]) * (xyz_jet_global[0] - xyz_track_global[0]) +
                    (xyz_jet_global[1] - xyz_track_global[1]) * (xyz_jet_global[1] - xyz_track_global[1])+
                    (xyz_jet_global[2] - xyz_track_global[2]) * (xyz_jet_global[2]- xyz_track_global[2]));
            if(bdecaylength>0) lineardecaylength=bdecaylength;
            if(bdecaylength >3. ) return kFALSE;
            if(dcajetrack >0.07 ) return kFALSE;
            if(abs(impar[1])>2.) return kFALSE;
            if(abs(impar[0])>1.) return kFALSE;
            return kTRUE;
        }
    return kFALSE;
}

Double_t AliAnalysisTaskHFJetIPQA::GetValImpactParameter(TTypeImpPar type,Float_t *impar, Float_t * cov)
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
}// ########################################################################################Track Selection
Bool_t AliAnalysisTaskHFJetIPQA::IsTrackAccepted(const AliESDtrack* track ,Int_t n){
    if(track->GetNcls(0)<abs(n)) return kFALSE;
    fESDTrackCut->SetMinNClustersITS(abs(n));
    fESDTrackCut->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);

    Double_t relY = TMath::Abs(sqrt(((AliESDtrack*)track)->GetSigmaY2())/((AliESDtrack*)track)->GetY());
    Double_t relZ = TMath::Abs(sqrt(((AliESDtrack*)track)->GetSigmaZ2())/((AliESDtrack*)track)->GetZ());
    Double_t pxyz[3]={0} ;
    ((AliESDtrack*)track)->GetPxPyPz(pxyz);
    if(sqrt(((AliESDtrack*)track)->GetSigmaY2())>5E-2)return kFALSE;
    if(sqrt(((AliESDtrack*)track)->GetSigmaZ2())>5E-2)return kFALSE;
    if(relY>0.20) return kFALSE;
    if(relZ>0.20) return kFALSE;
    if(!(fESDTrackCut->AcceptTrack((AliESDtrack*)track))) return kFALSE;
    else  return kTRUE;
    return kFALSE;
}
// ######################################################################################## Jet matching 1/4
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
// ######################################################################################## Jet matching 2/4
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

Bool_t AliAnalysisTaskHFJetIPQA::IsTruePrimary(AliMCParticle * mcpart){
    if(!mcpart) return kFALSE;
    AliMCParticle * mcmother = (AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother());
    if(!mcmother) return kTRUE;
    Int_t istatus = mcmother->Particle()->GetStatusCode();
    if(istatus >11)return kTRUE;
    Int_t ipdg =abs(mcmother->PdgCode())	;
    return kFALSE;
}
// ######################################################################################## Jet matching 3/4
Double_t AliAnalysisTaskHFJetIPQA::GetWeightFactor( AliMCParticle * mcpart,Int_t &pCorr_indx){
    if(!mcpart)  return 1;
    Bool_t _particlesourcefound(kFALSE);
    Int_t  _particlesourcepdg(abs(mcpart->PdgCode()));
    Int_t  _particlesourceidx(-1);
    Double_t _particlesourcept(0);

    AliMCParticle * mcpartclone = mcpart;
    while(mcpart){//omega and xi test
            if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;
            _particlesourcept = mcpart->Pt();
            _particlesourcepdg = abs(mcpart->PdgCode());
            if (IsSelectionParticleOmegaXiSigmaP(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
                    _particlesourcefound = kTRUE;
                    break;
                }
            mcpart->GetMother() >0 ? mcpart =(AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother()) :mcpart =  0x0;
        }
    if (!_particlesourcefound) { //heavy mesons to improve templates
            mcpart = mcpartclone;

            while(mcpart){
                    if((abs(mcpart->PdgCode()) >0 && abs(mcpart->PdgCode()) <7)|| (abs(mcpart->PdgCode())  == 21)) break;

                    _particlesourcept = mcpart->Pt();
                    _particlesourcepdg = abs(mcpart->PdgCode());
                    if (IsSelectionParticleStrange(mcpart,_particlesourcepdg,_particlesourcept,_particlesourceidx)){
                            _particlesourcefound = kTRUE;
                            break;
                        }
                    mcpart->GetMother() >0 ? mcpart =(AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother()) :mcpart =  0x0;
                }
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
                    mcpart->GetMother() >0 ? mcpart =(AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother()) :mcpart =  0x0;
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
                    mcpart->GetMother() >0 ? mcpart =(AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother()) :mcpart =  0x0;
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
                    mcpart->GetMother() >0 ? mcpart =(AliMCParticle*)MCEvent()->GetTrack(mcpart->GetMother()) :mcpart =  0x0;
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
    factor = fBackgroundFactorLinus[_particlesourceidx][bin];
    pCorr_indx = mcpart->GetLabel();
    Double_t flucafactor = 1;
    switch(mcpart->PdgCode())
        {


        case -bPhi:
            factor=1;
            flucafactor =1./fPhi->Eval(_particlesourcept);
            break;
        case -bK0S892:
            factor=1;
            flucafactor =1./fK0Star->Eval(_particlesourcept);
            break;
        case -bK0S892plus:
            factor=1;
            flucafactor =1./fK0Star->Eval(_particlesourcept);
            break;
        case -bOmegaBaryon:
            factor=1;
            flucafactor =1./fGraphOmega->Eval(_particlesourcept);
            break;
        case -bXiBaryon:
            factor=1;
            flucafactor =1./fGraphXi->Eval(_particlesourcept);
            break;
        case bPhi:
            factor=1;
            flucafactor =1./fPhi->Eval(_particlesourcept);
            break;

        case bK0S892:
            factor=1;
            flucafactor =1./fK0Star->Eval(_particlesourcept);
            break;
        case bK0S892plus:
            factor=1;
            flucafactor =1./fK0Star->Eval(_particlesourcept);
            break;
        case bOmegaBaryon:
            factor=1;
            flucafactor =fGraphOmega->Eval(_particlesourcept);
            break;
        case bXiBaryon:
            factor=1;
            flucafactor =fGraphXi->Eval(_particlesourcept);
            break;
            /*
        case bLambda:
       factor=1;
            flucafactor =fGeant3FlukaLambda->Eval(_particlesourcept);
            break;
        case -bLambda:
            flucafactor =fGeant3FlukaAntiLambda->Eval(_particlesourcept);
            break;
        case bProton:
            flucafactor =fGeant3FlukaProton->Eval(_particlesourcept);
            break;
        case -bProton:
            flucafactor =fGeant3FlukaAntiProton->Eval(_particlesourcept);
            break;
        case -bKaon:
            flucafactor =fGeant3FlukaKMinus->Eval(_particlesourcept);
            break;
*/
        default:
            break;
        }

    factor*=flucafactor;

    if (factor <= 0 || factor > 6.)  factor = 1.;
    return factor;
}

Bool_t AliAnalysisTaskHFJetIPQA::IsSecondaryFromWeakDecay( AliMCParticle * particle ) {
    // If a particle is not a physical primary, check if it comes from weak decay
    Int_t mfl = 0;
    Int_t indexMoth = particle->GetMother();
    if(indexMoth < 0) return kFALSE; // if index mother < 0 and not a physical primary, is a non-stable product or one of the beams
    AliMCParticle* moth =(AliMCParticle*) MCEvent()->GetTrack(indexMoth);
    if(!moth) return kFALSE;
    Int_t codemoth = TMath::Abs(moth->PdgCode());
    Int_t checklist [9] ={310,130,3122,3212,3222,3112,3322,3312,3334};
    for (Int_t i = 0; i<9; ++i) {
            if(codemoth == checklist[i])return kTRUE;

        }
    return kFALSE;
}

Bool_t AliAnalysisTaskHFJetIPQA::GetBMesonWeight( AliMCParticle * mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx  )
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

Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticle( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
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

Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleALICE( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
    pT 	= mcpart->Pt();
    AliMCParticle * mother = 0x0;
    idx = -1;
    pdg = abs(mcpart->PdgCode());
    Bool_t pIsSecStrANGE= MCEvent()->Stack()->IsSecondaryFromWeakDecay(mcpart->GetLabel());

    if(IsSecondaryFromWeakDecay(mcpart) && pIsSecStrANGE ) return kFALSE;

    switch(pdg){
        case bProton:
            mother = (AliMCParticle*)(MCEvent()->GetTrack(mcpart->GetMother()));
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

Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleMeson( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
    pT 	= mcpart->Pt();
    AliMCParticle * mother = 0x0;
    idx = -1;
    pdg = abs(mcpart->PdgCode());
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

Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleOmegaXiSigmaP( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
    pT 	= mcpart->Pt();
    idx = -1;
    pdg = abs(mcpart->PdgCode());
    Bool_t pIsPhysicalPrimary= MCEvent()->IsPhysicalPrimary(mcpart->GetLabel());
    if (!pIsPhysicalPrimary) return kFALSE;
    switch(pdg){
        case bXiBaryon:
            idx = bIdxBStarPlus;//dummy!
            return kTRUE;
            break;
        case bOmegaBaryon:
            idx = bIdxBStarPlus;//dummy!
            return kTRUE;
            break;
        default:
            return kFALSE;
            break;
        }
    return kFALSE;
}

Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticleStrange( AliMCParticle *  mcpart ,Int_t &pdg,Double_t &pT,Int_t &idx ){
    pT 	= mcpart->Pt();
    AliMCParticle * mother = 0x0;
    idx = -1;
    pdg = abs(mcpart->PdgCode());
    Bool_t pIsPhysicalPrimary= MCEvent()->IsPhysicalPrimary(mcpart->GetLabel());
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
            mother = (AliMCParticle*)(MCEvent()->GetTrack(mcpart->GetMother()));
            if(mother){
                    if((abs(mother->PdgCode()) == bPhi))
                        return kFALSE;
                }
            return kTRUE;
            break;
        case bK0l:
            idx = bIdxK0s;
            mother = (AliMCParticle*)(MCEvent()->GetTrack(mcpart->GetMother()));
            if(mother){
                    if((abs(mother->PdgCode()) == bPhi))
                        return kFALSE;
                }
            return kTRUE;
            break;
        case bLambda:
            idx = bIdxLambda;

            mother = (AliMCParticle*)(MCEvent()->GetTrack(mcpart->GetMother()));
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

Bool_t AliAnalysisTaskHFJetIPQA::IsPromptBMeson(AliMCParticle * part )
{
    if(!part) return kFALSE;
    Int_t pdg = TMath::Abs(part->PdgCode());
    if ((pdg >= 500 && pdg < 600 )||(pdg >= 5000 && pdg < 6000 ))
        {
            Int_t imo =  part->GetMother();
            AliMCParticle* pm = dynamic_cast<AliMCParticle *>(MCEvent()->GetTrack(imo));
            Int_t mpdg = TMath::Abs(pm->PdgCode());
            if (!(mpdg >5000 && mpdg <6000) && !(mpdg >500 && mpdg <600))
                return kTRUE;
        }
    return kFALSE;
}

Bool_t AliAnalysisTaskHFJetIPQA::IsPromptDMeson(AliMCParticle * part )
{
    if(!part) return kFALSE;
    Int_t pdg = TMath::Abs(part->PdgCode());
    if ((pdg >= 400 && pdg < 500 )||(pdg >= 4000 && pdg < 5000 ))
        {
            Int_t imo =  part->GetMother();
            if(imo<0) return kTRUE;
            AliMCParticle* pm = ((AliMCParticle*)fMCEvent->GetTrack(abs(imo)));
            if(!pm) return kTRUE;
            Int_t mpdg = TMath::Abs(pm->PdgCode());
            if (!(mpdg >4000 && mpdg <6000) && !(mpdg >400 && mpdg <600))
                return kTRUE;
        }

    return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::ParticleIsPossibleSource(Int_t pdg){
    Int_t pos[22] = {bPi0,bEta,bEtaPrime,bPhi,bRho,bOmega,bK0s,bLambda,bOmegaBaryon,bXiBaryon,bD0,bPi,bKaon,bProton,bDPlus,bDStarPlus,bDSPlus,bLambdaB,bLambdaC,bBPlus,bB0,bBStarPlus};
    for (Int_t i =0 ;i<22 ;++i){
            if (abs(pdg)==pos[i] ) return kTRUE;
        }
    return kFALSE;
}

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

void AliAnalysisTaskHFJetIPQA::GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const
{
    Double_t deta = jet2->Eta() - jet1->Eta();
    Double_t dphi = jet2->Phi() - jet1->Phi();
    dphi = TVector2::Phi_mpi_pi(dphi);
    d = sqrt(deta * deta + dphi * dphi);
}

Double_t AliAnalysisTaskHFJetIPQA::GetMonteCarloCorrectionFactor(AliVTrack* track,Int_t &pCorr_indx){
    if(fDisableWeightingMC) return 1;

    AliMCParticle *pMCESD = 0x0;
    if(track->GetLabel()< 0) return 1;
    pMCESD = ((AliMCParticle*)MCEvent()->GetTrack(abs(track->GetLabel())));
    if(!(pMCESD)) return 1;
    Double_t val = 1;
    val=  GetWeightFactor(pMCESD,pCorr_indx);
    if(val > 0 ) return val;

    return 1.;
}

Bool_t AliAnalysisTaskHFJetIPQA::mysort(const SJetIpPati& i, const SJetIpPati& j)
{
    if(i.first <= j.first)
        return kFALSE;
    else
        return kTRUE;
}

Double_t AliAnalysisTaskHFJetIPQA::GetPtCorrected(const AliEmcalJet *jet)
{
    AliJetContainer * jetconrec = 0x0;
    jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    if(jet && jetconrec)
        return jet->Pt() - jetconrec->GetRhoVal() * jet->Area();
    return -1.;
}

Double_t AliAnalysisTaskHFJetIPQA::GetPtCorrectedMC(const AliEmcalJet *jet)
{
    AliJetContainer * jetconrec = 0x0;
    jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(1));
    if(jet && jetconrec)
        return jet->Pt() - jetconrec->GetRhoVal() * jet->Area();
    return -1.;
}


//###############################################################################################################
//Functions for System 8
Bool_t AliAnalysisTaskHFJetIPQA::IsJetTaggedTC(Int_t n, Double_t thres)
{
    return kTRUE;
}

Bool_t AliAnalysisTaskHFJetIPQA::IsJetTaggedJetProb(Double_t thresProb)
{
    return kTRUE;
}

Int_t  AliAnalysisTaskHFJetIPQA::IsMCJetPartonFast(const AliEmcalJet *jet, Double_t radius,Bool_t &is_udg)
{
    if(!jet) return 0;
    if(! fMcEvtSampled){
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

void AliAnalysisTaskHFJetIPQA::FillHist(const char *name, Double_t x, Double_t w){
    TH1D * h1 =GetHist1D(name);
    h1->Fill(x,w);
}
void AliAnalysisTaskHFJetIPQA::FillHist(const char *name, Double_t x, Double_t y, Double_t w){
    TH2D * h2 =GetHist2D(name);
    h2->Fill(x,y,w);
}
void AliAnalysisTaskHFJetIPQA::IncHist(const char *name, Int_t bin){
    TH1D * h1 =GetHist1D(name);
    h1->SetBinContent(bin,h1->GetBinContent(bin)+1);
}
TH1 *AliAnalysisTaskHFJetIPQA::AddHistogramm(const char *name, const char *title, Int_t x, Double_t xlow, Double_t xhigh, Int_t y, Double_t ylow, Double_t yhigh){
    TObject * res = 0x0;
    res = fOutput2->FindObject(name);
    if((res)) return 0x0;

    TH1 * phist=0x0;
    if(y==0){ //TH1D*
            phist = new TH1D (name,title,x,xlow,xhigh);
        }
    else  {
            phist = new TH2D(name,title,x,xlow,xhigh,y,ylow,yhigh);
        }
    phist->Sumw2();

    fOutput2->Add(phist);
    Printf("Adding %s to output list",phist->GetName());
    return (TH1*)phist;
}

void AliAnalysisTaskHFJetIPQA::SubtractMean(Double_t val[], AliVTrack *track){
    Double_t  deltamean=fGraphMean->Eval(track->Pt() < 3. ? track->Pt() : 3 );//(fCurrentMeanFactors[0]-fCurrentMeanFactors[1]*TMath::Exp(-1*fCurrentMeanFactors[2]*(track->Pt()-fCurrentMeanFactors[3]))) *1e-4;

    val[0] -=deltamean*1e-4;;
}
