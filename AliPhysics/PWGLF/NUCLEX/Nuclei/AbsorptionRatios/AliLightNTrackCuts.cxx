/*
 * AliLightNTrackCuts.cxx
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#include "AliLightNTrackCuts.h"
ClassImp(AliLightNTrackCuts)
AliLightNTrackCuts::AliLightNTrackCuts()
:fMCHists(0)
,fHists(0)
,fMCData(false)
,fDCAPlots(false)
,fCombSigma(false)
,f3DPlotPmass2dca(false)
,fContribSplitting(false)
,fFillQALater(false)
,fCheckFilterBit(false)
,fCheckPileUpITS(false)
,fCheckPileUpTOF(false)
,fCheckPileUp(false)
,fFilterBit(0)
,fpTmin(0.)
,fpTmax(0.)
,fcutPt(false)
,fetamin(0.)
,fetamax(0.)
,fcutEta(false)
,fCutRapidity(false)
,fcutCharge(false)
,fCharge(0)
,fnTPCCls(0)
,fcutnTPCCls(false)
,fDCAProp(false)
,fDCAToVertexXY(0)
,fCutDCAToVtxXY(false)
,fdoITSnSigmaCut(false)
,fDCAToVertexZ(0)
,fCutDCAToVtxZ(false)
,fMinMass(0.)
,fMaxMass(0.)
,fCutSharedCls(false)
,fCheckTPCRefit(false)
,fCutTPCCrossedRows(false)
,fCutPID(false)
,fCutHighPtSig(false)
,fParticleID(AliPID::kUnknown)
,fNSigValueTPC()
,fNSigValueTOF()
,fNSigValueITSmin()
,fNSigValueITSmax()
,fPIDPTPCThreshold(0)
,fRejectPions(false)
{}

AliLightNTrackCuts::~AliLightNTrackCuts() {
    if (fMCHists) {
        delete fMCHists;
    }
    if (fHists) {
        delete fHists;
    }
}

bool AliLightNTrackCuts::isSelected(AliLightNTrack *Track) {
    if (!Track) {
        AliFatal("No Input Track recieved");
    }
    bool pass=true;
    bool ForMassFitPass=true;
    bool PIDEffPass=true;
    if (!Track->IsSet()) {
        pass=false;
        ForMassFitPass=false;
        PIDEffPass=false;
    } else {
        fHists->FillTrackCounter(0);
    }
    if (pass) {
        if (!TrackingCuts(Track)) {
            pass=false;
            ForMassFitPass=false;
            PIDEffPass=false;
        }
    }
    
    if (pass) {
        if (!TPCPIDAODCuts(Track)) {
            pass=false;
            ForMassFitPass=false;
        }
    }
    
    if (pass) {
        if (!PIDAODCuts(Track)) {
            pass=false;
            ForMassFitPass=false;
        } else {
            fHists->FillTrackCounter(20);
        }
    }
    
    if (pass && fdoITSnSigmaCut) {
        if (!ITSPIDAODCuts(Track)) {
            pass=false;
            ForMassFitPass=false;
            PIDEffPass=false;
        }
    }
    
    if (pass) {
        if (!DCACuts(Track)) {
            pass=false;
            ForMassFitPass=false;
            PIDEffPass=false;
        }
    }
    
    //The mass distribution without a TOF PID cut, (to make corrections with a fit later)
    //And a 3D histogram p,mass,dca
    double p =Track->GetP();
    for (int i=0;i<2;++i) {
        if (i==0||(i==1&&ForMassFitPass)) {
            fHists->FillMass2sq(i,p,Track->GetMassSquare());
            
            //Careful: uses high diskspace, better use combination of FillMass2sq and
            //FillDCAXYPBinningTot
            if(f3DPlotPmass2dca){
                if (fDCAProp) {
                    fHists->FillP_mass2_DCAxy(p,Track->GetMassSquare(),Track->GetDCAXYProp());
                }else {
                    fHists->FillP_mass2_DCAxy(p,Track->GetMassSquare(),Track->GetDCAXY());
                }
            }
        }
    }
    
    Track->SetUse(pass);
    if (!fFillQALater) {
        BookQA(Track);
        if (fMCData) {
            BookMC(Track);
            //The momentum distribution of correct identified particles without a PID
            //cut but with the requirement that there is a detector signal (particle reached the detector)
            double p =Track->GetP();
            if (PIDEffPass){
                bool TPCisthere=false;
                bool TOFisthere=false;
                if (Track->GetstatusTPC()==AliPIDResponse::kDetPidOk) {
                    TPCisthere=true;
                }
                if (Track->GetstatusTOF()==AliPIDResponse::kDetPidOk) {
                    TOFisthere=true;
                }
                Int_t PDGcode[6] = {11,13,211,321,2212,1000010020};
                if(TMath::Abs(Track->GetMCPDGCode())==PDGcode[(int)(fParticleID)]){
                    if (p < fPIDPTPCThreshold) {
                        if(TPCisthere&&(Track->GetParticleOrigin()==0)) fMCHists->FillMCCorrPrimNoPID(p);
                    }else{
                        if(TPCisthere&&TOFisthere&&(Track->GetParticleOrigin()==0)) fMCHists->FillMCCorrPrimNoPID(p);
                    }
                }
            }
        }
    }
    return pass;
}

bool AliLightNTrackCuts::TrackingCuts(AliLightNTrack *Track) {
    bool pass=true;
    std::vector<double> eta=Track->GetEta();
    std::vector<int> charge=Track->GetCharge();
    double rapidity = Track->GetRapidity(fParticleID);
    double p =Track->GetP();
    
    if (fCheckFilterBit) {
        if (!Track->TestFilterBit(fFilterBit)) {
            //if there is Filterbit -1 .. see Prong cuts! Daughters don't check for FB
            pass=false;
        } else {
            fHists->FillTrackCounter(1);
        }
    }
    if (pass && fcutPt) {
        if (p < fpTmin||p > fpTmax) {
            pass=false;
        } else {
            fHists->FillTrackCounter(2);
        }
    }
    if (pass && fcutEta) {
        if (eta[0]< fetamin||eta[0]>fetamax) {
            pass=false;
        } else {
            fHists->FillTrackCounter(3);
        }
    }
    if (pass && fcutCharge) {
        if (!(charge[0]==fCharge)) {
            pass=false;
        } else {
            fHists->FillTrackCounter(4);
        }
    }
    if (pass&&fCheckPileUpITS) {
        if (!Track->GetHasITSHit()) {
            pass=false;
        } else {
            fHists->FillTrackCounter(5);
        }
    }
    if (pass&&fCheckPileUpTOF) {
        if (!Track->GetTOFTimingReuqirement()) {
            pass=false;
        } else {
            fHists->FillTrackCounter(6);
        }
    }
    if (pass&&fCheckPileUp) {
        if (!(Track->GetTOFTimingReuqirement()||Track->GetHasITSHit())) {
            pass=false;
        } else {
            fHists->FillTrackCounter(7);
        }
    }
    if (pass && fcutnTPCCls) {
        if (Track->GetNClsTPC()<fnTPCCls) {
            pass=false;
        } else {
            fHists->FillTrackCounter(8);
        }
    }
    if (pass && fcutnITSCls) {
        if (Track->GetNClsITS()<fnITSCls) {
            pass=false;
        } else {
            //fHists->FillTrackCounter(8);
        }
    }
    if (pass && fCutSharedCls) {
        if (!Track->isnoSharedClst()) {
            pass=false;
        } else {
            fHists->FillTrackCounter(9);
        }
    }
    if (pass && fCheckTPCRefit) {
        if (!Track->GetHasTPCRefit()) {
            pass=false;
        } else {
            fHists->FillTrackCounter(10);
        }
    }
    if (pass && fCutTPCCrossedRows) {
        if (Track->GetTPCCrossedRows()<fTPCCrossedRowsCut) {
            pass=false;
        } else {
            fHists->FillTrackCounter(11);
        }
        if (pass) {
            if (Track->GetRatioCr()<fTPCRatio) {
                pass=false;
            } else {
                fHists->FillTrackCounter(12);
            }
        }
    }
    if (pass && fCutRapidity){
        if(rapidity<fRapMin || rapidity>fRapMax) {
            pass=false;
        } else {
            fHists->FillTrackCounter(13);
        }
    }
    if (pass && fCutChi2){
        if(Track->GetTrackChi2perNDF()>fChi2perNDF) {
            pass=false;
        }
    }
    
    return pass;
}


bool AliLightNTrackCuts::ITSPIDAODCuts(AliLightNTrack *Track) {
    //ITS PID cut for (anti-)deuterons in the momentum region 0 < p < 1.4 GeV/c
    bool pass=true;
    bool ITSisthere=false;
    
    if (Track->GetstatusITS()==AliPIDResponse::kDetPidOk) {
        ITSisthere=true;
    }
    
    double p =Track->GetP();
    if (p<1.4) {
        double nSigITS=(Track->GetnSigmaITS((int)(fParticleID)));
        if (nSigITS < fNSigValueITSmin || nSigITS > fNSigValueITSmax) {
            pass=false;
        }
    }
    return pass;
}

bool AliLightNTrackCuts::TPCPIDAODCuts(AliLightNTrack *Track) {
    bool pass=true;
    //PID Method with an nSigma cut, just use the TPC below threshold,
    //and above TPC and TOF Combined
    //there must be TPC & TOF signal (TOF for P>0.7 GeV/c)
    bool TPCisthere=false;
    
    if (Track->GetstatusTPC()==AliPIDResponse::kDetPidOk) {
        TPCisthere=true;
        fHists->FillTrackCounter(14);
    }
    //Below a threshold where the bands are well seperated in the TPC use only
    //TPC for PID, since the TOF has only limited matching efficiency. Above
    //threshold use both detectors and perform a purity check, if another
    //particle species doesn't have a smaller sigma value
    
    if (!TPCisthere) {
        pass=false;
    }else{
        double nSigTPC=(Track->GetnSigmaTPC((int)(fParticleID)));
        if (!(TMath::Abs(nSigTPC)<fNSigValueTPC)) {
            pass=false;
        }
    }
    return pass;
}


bool AliLightNTrackCuts::PIDAODCuts(AliLightNTrack *Track) {
    bool pass=true;
    //PID Method with an nSigma cut, just use the TPC below threshold,
    //and above TPC and TOF Combined
    //there must be TPC & TOF signal (TOF for P>0.75 GeV/c)
    bool TPCisthere=false;
    bool TOFisthere=false;
    
    if (Track->GetstatusTPC()==AliPIDResponse::kDetPidOk) {
        TPCisthere=true;
    }
    if (Track->GetstatusTOF()==AliPIDResponse::kDetPidOk) {
        TOFisthere=true;
    }
    //Below a threshold where the bands are well seperated in the TPC use only
    //TPC for PID, since the TOF has only limited matching efficiency. Above
    //threshold use both detectors and perform a purity check, if another
    //particle species doesn't have a smaller sigma value
    double p =Track->GetP();
    if (p < fPIDPTPCThreshold) {
        if (!TPCisthere) {
            pass=false;
        } else {
            if (pass) {
                double nSigTPC=(Track->GetnSigmaTPC((int)(fParticleID)));
                if (!(TMath::Abs(nSigTPC)<fNSigValueTPC)) {
                    pass=false;
                } else {
                    fHists->FillTrackCounter(15);
                }
            }
            if (fRejectPions&&TOFisthere) {
                double nSigTOFPion=(Track->GetnSigmaTOF((int)(AliPID::kPion)));
                if (TMath::Abs(nSigTOFPion)<fNSigValueTPC) {
                    if (fParticleID==AliPID::kPion) {
                        AliWarning("Sure you want to use this method? Propably want to set"
                                   " SetRejLowPtPionsTOF(kFALSE), since you are selecting Pions");
                    }
                    //if the particle is a Pion according to the TOF, reject it!
                    pass=false;
                } else {
                    fHists->FillTrackCounter(16);
                }
            }
        }
        
    } else {
        if (!(TPCisthere&&TOFisthere)) {
            pass=false;
        } else {
            fHists->FillTrackCounter(17);
            double nSigTPC=(Track->GetnSigmaTPC((int)(fParticleID)));
            double nSigTOF=(Track->GetnSigmaTOF((int)(fParticleID)));
            double nSigComb=TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
            if (!(nSigComb < TMath::Sqrt(fNSigValueTPC*fNSigValueTPC + fNSigValueTOF*fNSigValueTOF))) {
                pass=false;
            } else {
                fHists->FillTrackCounter(18);
                if (fCutHighPtSig) {
                    if (!SmallestNSig(Track)) {
                        pass=false;
                    } else {
                        fHists->FillTrackCounter(19);
                    }
                }
            }
        }
    }
    return pass;
}

bool AliLightNTrackCuts::SmallestNSig(AliLightNTrack *Track) {
    bool pass=true;
    //check before if TPC and TOF PID are available
    //This should just be for PID of high pT particles
    AliPID::EParticleType type[6]={AliPID::kElectron,AliPID::kMuon,AliPID::kPion,
        AliPID::kKaon,AliPID::kProton,AliPID::kDeuteron};
    double nSigmaComb[6];
    //Form the combination:
    for (int i=0; i<6; ++i) {
        nSigmaComb[i]=TMath::Sqrt(pow((Track->GetnSigmaTPC(i)),2.)+
                                  pow((Track->GetnSigmaTOF(i)),2.));
    }
    int index=0;
    for (int i=0; i<6; ++i) {
        if (nSigmaComb[index]>nSigmaComb[i]) {
            index=i;
        }
    }
    if (!(type[index]==fParticleID)) {
        pass=false;
    }
    return pass;
}


bool AliLightNTrackCuts::MassCut_ForDCA(AliLightNTrack *Track) {
    bool pass = true;
    double p =Track->GetP();
    if (!(Track->GetMassSquare()>fMinMass && Track->GetMassSquare()<fMaxMass)){
        if(p > fPIDPTPCThreshold){
            pass = false;
        }
    }
    return pass;
}


bool AliLightNTrackCuts::DCACuts(AliLightNTrack *Track) {
    bool pass=true;
    double p =Track->GetP();
    if (fCutDCAToVtxZ) {
        if (fDCAProp) {
            if (!(TMath::Abs(Track->GetDCAZProp())<fDCAToVertexZ)) {
                pass=false;
            } else {
                fHists->FillTrackCounter(21);
            }
        } else {
            if (!(TMath::Abs(Track->GetDCAZ())<fDCAToVertexZ)) {
                pass=false;
            } else {
                fHists->FillTrackCounter(21);
            }
        }
    }
    
    //Making lose massCut as a first selection to avoid contamination from
    //other particles: MassCut_ForDCA(Track)
    if (pass&&fDCAPlots&&MassCut_ForDCA(Track)) {
        if (fDCAProp) {
            fHists->FillDCAXYPBins(p,Track->GetDCAXYProp());
        } else {
            fHists->FillDCAXYPBins(p,Track->GetDCAXY());
        }
        if (fMCData) {
            if (fDCAPlots) {
                if (fDCAProp) {
                    fMCHists->FillMCDCAXYPtBins(Track->GetParticleOrigin(),
                                                Track->GetMotherWeak(),
                                                p,
                                                Track->GetDCAXYProp());
                } else {
                    fMCHists->FillMCDCAXYPtBins(Track->GetParticleOrigin(),
                                                Track->GetMotherWeak(),
                                                p,
                                                Track->GetDCAXY());
                }
            }
        }
    }
    if (pass&&fCutDCAToVtxXY) {
        if (fDCAProp) {
            if (!(TMath::Abs(Track->GetDCAXYProp())<fDCAToVertexXY)) {
                pass=false;
            } else {
                fHists->FillTrackCounter(22);
            }
        } else {
            if (!(TMath::Abs(Track->GetDCAXY())<fDCAToVertexXY)) {
                pass=false;
            } else {
                fHists->FillTrackCounter(22);
            }
        }
    }
    return pass;
}
void AliLightNTrackCuts::Init() {
    fHists=new AliLightNTrackHist(fDCAPlots,fCombSigma,f3DPlotPmass2dca);
    if (fMCData) {
        fMCHists=new AliLightNTrackMCHist(fContribSplitting,fDCAPlots);
    }
    BookTrackCuts();
}

void AliLightNTrackCuts::BookQA(AliLightNTrack *Track) {
    std::vector<double> eta=Track->GetEta();
    std::vector<double> phi=Track->GetPhi();
    // double pT = Track->GetPt();
    double pTPC = Track->GetMomTPC();
    double p =Track->GetP();
    
    for (int i=0;i<2;++i) {
        if (i==0||(i==1&&Track->UseParticle())) {
            fHists->FilletaCut(i,eta.at(0));
            fHists->FillphiCut(i,phi.at(0));
            fHists->FillpCut(i,p);
            fHists->FillpTPCCut(i,pTPC);
            fHists->FillpDiff_p_pTPC(i,p,p-pTPC);
            fHists->FillTPCclsCut(i,Track->GetNClsTPC());
            if (fDCAProp) {
                fHists->FillDCAxyCut(i,p,Track->GetDCAXYProp());
                fHists->FillDCAzCut(i,p,Track->GetDCAZProp());
            } else {
                fHists->FillDCAxyCut(i,p,Track->GetDCAXY());
                fHists->FillDCAzCut(i,p,Track->GetDCAZ());
            }
            fHists->FillTPCCrossedRowCut(i,Track->GetTPCCrossedRows());
            fHists->FillTPCRatioCut(i,Track->GetRatioCr());
            fHists->FillTPCClsS(i,Track->GetTPCClsC());
            for (int j=0;j<6;++j) {
                if (Track->GetITSHit(j)) {
                    fHists->FillTPCClsCPileUp(i,j,Track->GetTPCClsC());
                } else if (Track->GetHasITSHit()||Track->GetTOFTimingReuqirement()){
                    fHists->FillTPCClsCPileUp(i,j+7,Track->GetTPCClsC());
                }
            }
            if (Track->GetTOFTimingReuqirement()) {
                fHists->FillTPCClsCPileUp(i,6,Track->GetTPCClsC());
            } else if (Track->GetHasITSHit()) {
                fHists->FillTPCClsCPileUp(i,13,Track->GetTPCClsC());
            } else {
                fHists->FillTPCClsCPileUp(i,14,Track->GetTPCClsC());
            }
            
            for (int j=0;j<6;++j) {
                if (Track->GetSharedClusterITS(j)) {
                    fHists->FillHasSharedClsITS(i,j+1,0);
                } else {
                    fHists->FillHasSharedClsITS(i,j+1,1);
                }
                for (int k=0;k<6;++k) {
                    if (Track->GetITSHit(k)) {
                        if (Track->GetSharedClusterITS(j)) {
                            fHists->FillITSSharedPileUp(i,k,j);
                        } else {
                            fHists->FillITSSharedPileUp(i,k,j+6);
                        }
                    } else if (Track->GetHasITSHit()||Track->GetTOFTimingReuqirement()) {
                        if (Track->GetSharedClusterITS(j)) {
                            fHists->FillITSSharedPileUp(i,k+7,j);
                        } else {
                            fHists->FillITSSharedPileUp(i,k+7,j+6);
                        }
                    }
                }
                if (Track->GetTOFTimingReuqirement()) {
                    if (Track->GetSharedClusterITS(j)) {
                        fHists->FillITSSharedPileUp(i,6,j);
                    } else {
                        fHists->FillITSSharedPileUp(i,6,j+6);
                    }
                } else if (Track->GetHasITSHit()) {
                    if (Track->GetSharedClusterITS(j)) {
                        fHists->FillITSSharedPileUp(i,13,j);
                    } else {
                        fHists->FillITSSharedPileUp(i,13,j+6);
                    }
                } else {
                    if (Track->GetSharedClusterITS(j)) {
                        fHists->FillITSSharedPileUp(i,14,j);
                    } else {
                        fHists->FillITSSharedPileUp(i,14,j+6);
                    }
                }
            }
            if(p < fPIDPTPCThreshold){
                fHists->FillEtaPhiTPConlyPID(i,eta.at(0),phi.at(0));
            }else{
                fHists->FillEtaPhiTPCTOFPID(i,eta.at(0),phi.at(0));
            }
            fHists->FillTPCdedx(i,p,Track->GetdEdxTPC());
            fHists->FillTOFbeta(i,p,Track->GetbetaTOF());
            fHists->FillITSdedx(i,p,Track->GetdEdxITS());
            fHists->FillNSigTPC(i,p,(Track->GetnSigmaTPC(fParticleID)));
            fHists->FillNSigTOF(i,p,(Track->GetnSigmaTOF(fParticleID)));
            fHists->FillNSigITS(i,p,(Track->GetnSigmaITS(fParticleID)));
            fHists->FillTPCStatus(i,Track->GetstatusTPC());
            fHists->FillTOFStatus(i,Track->GetstatusTOF());
            fHists->FillITSStatus(i,Track->GetstatusITS());
            fHists->FillrapidityCut(i,Track->GetRapidity(fParticleID));
            fHists->FillTrackChi2perNDF(i,Track->GetTrackChi2perNDF());
            //Fill These Before
            if (i==0&&fCombSigma) {
                fHists->FillNSigComb(p,Track->GetnSigmaTPC(fParticleID),Track->GetnSigmaTOF(fParticleID));
            }
            //Fill These only after
            if(i==1&&Track->UseParticle()){
                if(p>1 && p<2 && TMath::Abs(Track->GetDCAXYProp())<0.05 && TMath::Abs(Track->GetDCAZProp())<0.08 ){
                    fHists->FillTPCclsHighPur(i,Track->GetNClsTPC());
                }
            }
        }
    }
    return;
}

void AliLightNTrackCuts::BookMC(AliLightNTrack *Track) {
    if(!Track->TestFilterBit(1))return;
    double p =Track->GetP();
    double RAPIDITY = Track->GetRapidity(fParticleID);
    Int_t PDGcode[6] = {11,13,211,321,2212,1000010020};
    if (fpTmin<p && p<fpTmax) { 															//to be in same p range
        if (fetamin<Track->GetEta().at(0)&&Track->GetEta().at(0)<fetamax) { 				//to be in same eta range
            if(RAPIDITY>fRapMin && RAPIDITY<fRapMax){ 										//to be in same rapidity range
                if(!fcutCharge){
                    if(TMath::Abs(Track->GetMCPDGCode())==PDGcode[(int)(fParticleID)]){
                        fMCHists->FillMCGen(p);
                        if(Track->GetParticleOrigin() == 0){ 								//Primary
                            fMCHists->FillMCGenPrim(p);
                        }
                        
                    }
                }else{
                    Int_t sign = fCharge/TMath::Abs(fCharge);
                    if(Track->GetMCPDGCode() == sign*PDGcode[(int)(fParticleID)]){
                        fMCHists->FillMCGen(p);
                        if(Track->GetParticleOrigin() == 0){ 								//Primary
                            fMCHists->FillMCGenPrim(p);
                        }
                    }
                }
            }
        }
    }
    if (Track->UseParticle()) {
        // double pT = Track->GetPt();
        int PDGcode[6] = {11,13,211,321,2212,1000010020};
        //Fill Identified
        fMCHists->FillMCIdent(p);
        if (!fcutCharge) {
            if (TMath::Abs(Track->GetMCPDGCode())==TMath::Abs(PDGcode[(int)(fParticleID)])) {
                fMCHists->FillMCCorr(p);
                if(Track->GetParticleOrigin() == 0){                                        //Primary
                    fMCHists->FillMCCorrPrim(p);
                }
            } else {
                Track->SetParticleOrigin(AliLightNBasePart::kContamination);
            }
        } else {
            Int_t sign = fCharge/TMath::Abs(fCharge);
            if (Track->GetMCPDGCode() == sign*PDGcode[(int)(fParticleID)]) {
                fMCHists->FillMCCorr(p);
                if(Track->GetParticleOrigin() == 0){                                        //Primary
                    fMCHists->FillMCCorrPrim(p);
                }
            } else {
                Track->SetParticleOrigin(AliLightNBasePart::kContamination);
            }
        }
        if (fContribSplitting) {
            FillMCContributions(Track);
        }
    }
}

void AliLightNTrackCuts::FillMCContributions(
                                             AliLightNTrack *Track)
{
    
    double p =Track->GetP();
    // double pT=Track->GetPt();
    AliLightNBasePart::PartOrigin org=Track->GetParticleOrigin();
    Int_t iFill = -1;
    switch(org) {
        case AliLightNBasePart::kPhysPrimary:
            fMCHists->FillMCPrimary(p);
            iFill = 0;
            break;
        case AliLightNBasePart::kWeak:
            fMCHists->FillMCFeeddown(p,TMath::Abs(Track->GetMotherWeak()));
            iFill = 1;
            break;
        case AliLightNBasePart::kMaterial:
            fMCHists->FillMCMaterial(p);
            iFill = 2;
            break;
        case AliLightNBasePart::kContamination:
            fMCHists->FillMCCont(p);
            iFill = 3;
            break;
        default:
            AliFatal("Type Not implemented");
            break;
    }
    if (iFill >= 0 && iFill < 4) {
        std::vector<double> eta=Track->GetEta();
        std::vector<double> phi=Track->GetPhi();
        fMCHists->FillMCpTPCCut(iFill,p);
        fMCHists->FillMCetaCut(iFill,eta[0]);
        fMCHists->FillMCphiCut(iFill,phi[0]);
        fMCHists->FillMCTPCclsCut(iFill,p,Track->GetNClsTPC());
        if (fDCAProp) {
            fMCHists->FillMCDCAxyCut(iFill,p,Track->GetDCAXYProp());
            fMCHists->FillMCDCAzCut(iFill,p,Track->GetDCAZProp());
        } else {
            fMCHists->FillMCDCAxyCut(iFill,p,Track->GetDCAXY());
            fMCHists->FillMCDCAzCut(iFill,p,Track->GetDCAZ());
        }
        fMCHists->FillMCTPCCrossedRowCut(iFill,p,Track->GetTPCCrossedRows());
        fMCHists->FillMCTPCRatioCut(iFill,p,Track->GetRatioCr());
        fMCHists->FillMCTPCdedx(iFill,p,Track->GetdEdxTPC());
        fMCHists->FillMCTOFbeta(iFill,p,Track->GetbetaTOF());
        fMCHists->FillMCNSigTPC(iFill,p,Track->GetnSigmaTPC(fParticleID));
        fMCHists->FillMCNSigTOF(iFill,p,Track->GetnSigmaTOF(fParticleID));
    } else {
        TString errMSG =  Form("iFill = %d", iFill);
        AliFatal(errMSG.Data());
    }
    return;
}

void AliLightNTrackCuts::BookTrackCuts() {
    if (!fHists) {
        AliFatal("AliFemtoPPbpbLamSpTrackCuts::BookTrackCuts No Histograms to work with");
    }
    if (fcutPt) {
        fHists->FillConfig(0, fpTmin);
        fHists->FillConfig(1, fpTmax);
    }
    if (fcutEta) {
        fHists->FillConfig(2, fetamin);
        fHists->FillConfig(3, fetamax);
    }
    if (fcutCharge) {
        fHists->FillConfig(4, fCharge);
    }
    
    if (fcutnTPCCls) {
        fHists->FillConfig(5, fnTPCCls);
    }
    
    if (fCheckFilterBit) {
        fHists->FillConfig(6, fFilterBit);
    } else {
        fHists->FillConfig(6, -1);
    }
    
    if (fMCData) {
        fHists->FillConfig(7, 1);
    }
    if (fCutDCAToVtxXY) {
        fHists->FillConfig(8, fDCAToVertexXY);
    }
    if (fCutDCAToVtxZ) {
        fHists->FillConfig(9, fDCAToVertexZ);
    }
    if (fCutSharedCls) {
        fHists->FillConfig(10, 1);
    }
    if (fCutTPCCrossedRows) {
        fHists->FillConfig(11, 1);
    } else {
        fHists->FillConfig(11, 0);
    }
    if (fCutPID) {
        fHists->FillConfig(12, fPIDPTPCThreshold);
        fHists->FillConfig(13, fNSigValueTPC);
        if (fRejectPions) {
            fHists->FillConfig(14, 1);
        } else {
            fHists->FillConfig(14, 0);
        }
        if (fCutHighPtSig) {
            fHists->FillConfig(15, 1);
        } else {
            fHists->FillConfig(15,0);
        }
    } else {
        fHists->FillConfig(12, 0);
        fHists->FillConfig(13, 0);
        fHists->FillConfig(14, 0);
        fHists->FillConfig(15, 0);
    }
    if (fCheckPileUpITS) {
        fHists->FillConfig(16,1);
    }
    if (fCheckPileUpTOF) {
        fHists->FillConfig(17,1);
    }
    if (fCheckPileUp) {
        fHists->FillConfig(18,1);
    }
    if (fCheckTPCRefit) {
        fHists->FillConfig(19,1);
    }
    if (fCutRapidity) {
        fHists->FillConfig(20, fRapMin);
        fHists->FillConfig(21, fRapMax);
    }
}


AliLightNTrackCuts* AliLightNTrackCuts::PrimProtonCuts(bool isMC,bool DCAPlots,bool CombSigma,bool ContribSplitting)
{
    AliLightNTrackCuts *trackCuts = new AliLightNTrackCuts();
    //you can leave DCA cut active, this will still be filled
    //over the whole DCA_xy range
    trackCuts->SetPlotDCADist(DCAPlots);
    trackCuts->SetPlotCombSigma(CombSigma);
    trackCuts->SetPlotContrib(ContribSplitting);
    trackCuts->SetIsMonteCarlo(isMC);
    
    trackCuts->SetFilterBit(256);
    trackCuts->SetPtRange(0.1,1e30);
    trackCuts->SetEtaRange(-0.8, 0.8);
    trackCuts->SetNClsTPC(70);
    trackCuts->SetNClsITS(70);
    trackCuts->SetNClsITS(2); //Min. value
    trackCuts->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
    trackCuts->SetDCAVtxZ(0.2);
    trackCuts->SetDCAVtxXY(0.1);
    trackCuts->SetCutSharedCls(true);
    trackCuts->SetTPCRatioCut(0.8);
    trackCuts->SetChi2perNDFCut(4);
    trackCuts->SetCutTPCCrossedRows(true);
    trackCuts->SetTPCCrossedRowsCut(70);
    trackCuts->SetPID(AliPID::kProton, 0.7,3.,3.);
    trackCuts->SetRapidityRange(-1, 1);
    trackCuts->SetRejLowPtPionsTOF(false);
    trackCuts->SetCutSmallestSig(false);
    trackCuts->SetMassCut_ForDCA(0.3,1.8);
    
    return trackCuts;
}

AliLightNTrackCuts* AliLightNTrackCuts::DecayProtonCuts(bool isMC,bool ContribSplitting) {
    AliLightNTrackCuts *trackCuts = new AliLightNTrackCuts();
    trackCuts->SetPlotDCADist(false);
    trackCuts->SetPlotCombSigma(false);
    trackCuts->SetPlotContrib(ContribSplitting);
    trackCuts->SetIsMonteCarlo(isMC);
    trackCuts->SetFillQALater(true);
    
    trackCuts->SetCheckPileUp(true);
    trackCuts->SetCheckFilterBit(false);
    trackCuts->SetEtaRange(-0.8, 0.8);
    trackCuts->SetNClsTPC(70);
    trackCuts->SetDCAReCalculation(true);
    trackCuts->SetCutCharge(1);
    trackCuts->SetPID(AliPID::kProton, 999.,5,5);
    return trackCuts;
}



AliLightNTrackCuts* AliLightNTrackCuts::PrimDeuteronCuts(bool isMC,bool DCAPlots,bool CombSigma,bool ContribSplitting)
{
    AliLightNTrackCuts *trackCuts = new AliLightNTrackCuts();
    //you can leave DCA cut active, this will still be filled
    //over the whole DCA_xy range
    trackCuts->SetPlotDCADist(DCAPlots);
    trackCuts->SetPlotCombSigma(CombSigma);
    trackCuts->SetPlotContrib(ContribSplitting);
    trackCuts->SetIsMonteCarlo(isMC);
    
    trackCuts->SetFilterBit(256);
    trackCuts->SetPtRange(0.1,1e30);
    trackCuts->SetEtaRange(-0.8, 0.8);
    trackCuts->SetNClsTPC(70);
    trackCuts->SetNClsITS(2); //Min. value
    trackCuts->SetDCAReCalculation(true);//Get the dca from the PropagateToVetex
    trackCuts->SetDCAVtxZ(0.2);
    trackCuts->SetDCAVtxXY(0.1);
    trackCuts->SetCutSharedCls(true);
    trackCuts->SetTPCRatioCut(0.8);
    trackCuts->SetChi2perNDFCut(4);
    trackCuts->SetCutTPCCrossedRows(true);
    trackCuts->SetTPCCrossedRowsCut(70);
    trackCuts->SetPID(AliPID::kDeuteron, 1.4,3.,1e30);
    trackCuts->SetCutITSPID(-2.,1e30,true);
    trackCuts->SetRapidityRange(-1, 1);
    trackCuts->SetRejLowPtPionsTOF(false);
    trackCuts->SetCutSmallestSig(false);
    trackCuts->SetMassCut_ForDCA(3.0,5.0);
    
    return trackCuts;
}




AliLightNTrackCuts* AliLightNTrackCuts::DecayPionCuts(bool isMC,bool ContribSplitting) {
    AliLightNTrackCuts *trackCuts = new AliLightNTrackCuts();
    trackCuts->SetPlotDCADist(false);
    trackCuts->SetPlotCombSigma(false);
    trackCuts->SetPlotContrib(ContribSplitting);
    trackCuts->SetIsMonteCarlo(isMC);
    trackCuts->SetFillQALater(true);
    
    trackCuts->SetCheckPileUp(true);
    trackCuts->SetCheckFilterBit(kFALSE);
    trackCuts->SetEtaRange(-0.8, 0.8);
    trackCuts->SetNClsTPC(70);
    trackCuts->SetDCAReCalculation(kTRUE);
    trackCuts->SetCutCharge(-1);
    trackCuts->SetPID(AliPID::kPion, 999.,5,5);
    return trackCuts;
}

AliLightNTrackCuts* AliLightNTrackCuts::Xiv0PionCuts(bool isMC,bool ContribSplitting) {
    AliLightNTrackCuts *trackCuts = new AliLightNTrackCuts();
    trackCuts->SetPlotDCADist(false);
    trackCuts->SetPlotCombSigma(false);
    trackCuts->SetCheckPileUp(false);
    trackCuts->SetPlotContrib(ContribSplitting);
    trackCuts->SetIsMonteCarlo(isMC);
    trackCuts->SetFillQALater(true);
    
    trackCuts->SetCheckFilterBit(kFALSE);
    trackCuts->SetEtaRange(-0.8, 0.8);
    trackCuts->SetPtRange(0.3,999);
    trackCuts->SetCutTPCCrossedRows(true);
    trackCuts->SetDCAReCalculation(kTRUE);
    trackCuts->SetCutCharge(-1);
    trackCuts->SetCheckTPCRefit(true);
    trackCuts->SetPID(AliPID::kPion, 999.,4,4);
    return trackCuts;
}

AliLightNTrackCuts* AliLightNTrackCuts::Xiv0ProtonCuts(bool isMC,bool ContribSplitting) {
    AliLightNTrackCuts *trackCuts = new AliLightNTrackCuts();
    trackCuts->SetPlotDCADist(false);
    trackCuts->SetPlotCombSigma(false);
    trackCuts->SetCheckPileUp(false);
    trackCuts->SetPlotContrib(ContribSplitting);
    trackCuts->SetIsMonteCarlo(isMC);
    trackCuts->SetFillQALater(true);
    
    trackCuts->SetCheckFilterBit(kFALSE);
    trackCuts->SetEtaRange(-0.8, 0.8);
    trackCuts->SetPtRange(0.3,999);
    trackCuts->SetCutTPCCrossedRows(true);
    trackCuts->SetDCAReCalculation(kTRUE);
    trackCuts->SetCutCharge(1);
    trackCuts->SetCheckTPCRefit(true);
    trackCuts->SetPID(AliPID::kProton, 999.,4,4);
    return trackCuts;
}

AliLightNTrackCuts* AliLightNTrackCuts::XiBachPionCuts(bool isMC,bool ContribSplitting) {
    AliLightNTrackCuts *trackCuts = new AliLightNTrackCuts();
    trackCuts->SetPlotDCADist(false);
    trackCuts->SetPlotCombSigma(false);
    trackCuts->SetCheckPileUp(false);
    trackCuts->SetPlotContrib(ContribSplitting);
    trackCuts->SetIsMonteCarlo(isMC);
    trackCuts->SetFillQALater(true);
    
    trackCuts->SetCheckFilterBit(kFALSE);
    trackCuts->SetEtaRange(-0.8, 0.8);
    trackCuts->SetPtRange(0.3,999);
    trackCuts->SetCutTPCCrossedRows(true);
    trackCuts->SetDCAReCalculation(kTRUE);
    trackCuts->SetCutCharge(-1);
    trackCuts->SetCheckTPCRefit(true);
    trackCuts->SetPID(AliPID::kPion, 999.,4,4);
    return trackCuts;
}
