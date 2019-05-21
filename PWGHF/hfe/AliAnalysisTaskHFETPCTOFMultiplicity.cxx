
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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


//////////////////////////////////////////////
//    Task for Measurement of Heavy Flavour //
//    electron as a function of charged     //
//    particle multiplicity                 //
//    Author: Preeti Dhankher               //
//////////////////////////////////////////////


#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "TVector3.h"
#include <TRandom3.h>
#include "TProfile.h"
#include "TGeoManager.h"

#include "stdio.h"
#include "iostream"
#include "fstream"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskHFETPCTOFMultiplicity.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliEventPoolManager.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDpid.h"
#include "AliAODPid.h"

#include "AliSelectNonHFE.h" // non hfe
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliAODTracklets.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "AliCentrality.h"
#include "AliMagF.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "AliESDUtils.h"
#include "AliVertexingHFUtils.h"
#include "AliNormalizationCounter.h"
#include "AliAnalysisVertexingHF.h"

#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"

class AliAnalysisTaskHFETPCTOFMultiplicity;
using namespace std;
ClassImp(AliAnalysisTaskHFETPCTOFMultiplicity)

AliAnalysisTaskHFETPCTOFMultiplicity::AliAnalysisTaskHFETPCTOFMultiplicity() : AliAnalysisTaskSE(),

//Event Cut
fCutNcontV(2),
//Track Cut
fCutTPCMaxCls(100.),
fCutTPCchi2perNDF(4.),
fCutTPCNCls(80.),
fCutITSNCls(3.),
fCutDCAxy(2.4),
fCutDCAz(3.2),
fCutTrackEta(0.7),
fCutNsigmaTOF(3.),
//PID Cut

fCutNsigmaEMin(-1.),
fCutNsigmaEMax(3.),


//Loose cuts for photonic electron pair
fAssoTPCCluster(80.),
fAssoITSCluster(3.),
fCutAssoEPt(0.1),
fCutAssoEEta(0.9),
fCutAssoENsigma(3),
fAssoITSRefit(kTRUE),
//Mass Cut for photonic electron pair
fCutInvmass(0.14),


// events
fAOD(0),
fNevents(0),
fOutputList(0),
fListProfiles(0),
fpidResponse(0),
// emcal correction
// MC
fMCArray(0),
fMCHeader(0),
fMCparticle(0),
// zvtx
fVtxZ(0),
fVtxX(0),
fVtxY(0),
// multi. estimation
fRejectPUFromSPD(kTRUE),

fSPDBoth(kTRUE),
fSPDAny(kFALSE),
fSPDFirst(kFALSE),

fRefMult(61.26),
gRandom(new TRandom3(0)),
fSparseMulti(0),
fvalueMulti(0),
// TPC info and PID
fTPCdEdx(0x0),
fTPCnsigma(0x0),
fTOFnsigma(0x0),
fTPCdEdxVsPTOFCut(0x0),
fTPCnSigmaVsPTOFCut(0x0),
fTPCnSigmaProtonsVsPTOFCut(0x0),
fTPCnSigmaKaonsVsPTOFCut(0x0),
fLandau(0),
fErr(0),
fNewHadFunc(0),
fHadCont_Landau(0),
fHadCont_Err(0),
fPt_incl_e_Landau(0),
fPt_incl_e_Err(0),
fTPCnSigma(-999.0),
fTrkPt(0),
fTrketa(0),
fTrkphi(0),
// electron info
// Non-Hfe
fNonHFE(new AliSelectNonHFE()),
// photonic electron info
fInvmassLS(0),
fInvmassULS(0),
fInvmassLSPt(0),
fInvmassULSPt(0),
fULSElecPt(0),
fLSElecPt(0),
fSparseLSElectron(0),
fSparseULSElectron(0),
fvaluePHElectron(0),
// MC info
fReadMC(kFALSE),
fNTotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),
// enhance weight cal.
ftype(-1),
fWeight(1),
fCalculateWeight(kFALSE),
fCalculateElectronEffi(kTRUE),
fSprsPi0EtaWeightCal(0),
fPi0Weight(0),
fEtaWeight(0),
fPi0Weight2(0),
fEtaWeight2(0),

fWeightEta(1),
fWeightPi0(1),

// inclusive e-
fInclsElecPt(0),
fInclsElecPtAll(0),
fInclsElecPtReco(0),
// hfe e-
fHFElecPtAll(0),
fHFElecPtReco_wtrkcuts(0),
fHFElecPtReco_wTPCPID(0),
fHFElecPtReco_wtrkCalocuts(0),
fHFElecPtReco_wTPCTOFPID(0),
//  nonhfe MC
fIsFrmPi0(kFALSE),
fIsFrmEta(kFALSE),

fNonHFeTrkPt(0),
fNonHFeWeightTrkPt(0),
fPi0eWeightTrkPt(0),
fEtaeWeightTrkPt(0),
fPi0EtaEleTrkPt(0),
//nonhfe invariant mass
fNonHFePairInvmassLS(0),
fNonHFePairInvmassULS(0),
// non hfe reconstructed
fRecoNonHFeTrkPt(0),
fRecoNonHFeWeightTrkPt(0),
fRecoPi0eWeightTrkPt(0),
fRecoEtaeWeightTrkPt(0),
fRecoPi0EtaEleTrkPt(0),

fRecoLSeTrkPt(0),
fRecoLSeWeightTrkPt(0),
fRecoPi0LSeWeightTrkPt(0),
fRecoEtaLSeWeightTrkPt(0),


fRecoULSeTrkPt(0),
fRecoULSeWeightTrkPt(0),
fRecoPi0ULSeWeightTrkPt(0),
fRecoEtaULSeWeightTrkPt(0)



{
    fvaluePHElectron = new Double_t[3];
    fvalueMulti = new Double_t[6];
    for(Int_t i=0; i<2; i++) fMultEstimatorAvg[i]=0;
    
}

//_____________________________________________________________________________
AliAnalysisTaskHFETPCTOFMultiplicity::AliAnalysisTaskHFETPCTOFMultiplicity(const char* name) : AliAnalysisTaskSE(name),


//Event Cut
fCutNcontV(2),
//Track Cut
fCutTPCMaxCls(100.),
fCutTPCchi2perNDF(4.),
fCutTPCNCls(80.),
fCutITSNCls(3.),
fCutDCAxy(2.4),
fCutDCAz(3.2),
fCutTrackEta(0.7),
fCutNsigmaTOF(3.),
//PID Cut

fCutNsigmaEMin(-1.),
fCutNsigmaEMax(3.),

//Loose cuts for photonic electron pair
fAssoTPCCluster(80.),
fAssoITSCluster(3.),
fCutAssoEPt(0.1),
fCutAssoEEta(0.9),
fCutAssoENsigma(3),
fAssoITSRefit(kTRUE),

fSPDBoth(kTRUE),
fSPDAny(kFALSE),
fSPDFirst(kFALSE),

//Mass Cut for photonic electron pair
fCutInvmass(0.14),


// events
fAOD(0),
fNevents(0),
fOutputList(0),
fListProfiles(0),
fpidResponse(0),
// emcal correction
// MC
fMCArray(0),
fMCHeader(0),
fMCparticle(0),
// zvtx
fVtxZ(0),
fVtxX(0),
fVtxY(0),
// multi. estimation
fRejectPUFromSPD(kTRUE),
fRefMult(61.26),
gRandom(new TRandom3(0)),
fSparseMulti(0),
fvalueMulti(0),
// TPC info and PID
fTPCdEdx(0x0),
fTPCnsigma(0x0),
fTOFnsigma(0x0),
fTPCdEdxVsPTOFCut(0x0),
fTPCnSigmaVsPTOFCut(0x0),
fTPCnSigmaProtonsVsPTOFCut(0x0),
fTPCnSigmaKaonsVsPTOFCut(0x0),
fLandau(0),
fErr(0),
fNewHadFunc(0),
fHadCont_Landau(0),
fHadCont_Err(0),
fPt_incl_e_Landau(0),
fPt_incl_e_Err(0),

fTPCnSigma(-999.0),
fTrkPt(0),
fTrketa(0),
fTrkphi(0),
// electron info

// Non-Hfe
fNonHFE(new AliSelectNonHFE()),

// photonic electron info
fInvmassLS(0),
fInvmassULS(0),
fInvmassLSPt(0),
fInvmassULSPt(0),
fULSElecPt(0),
fLSElecPt(0),
fSparseLSElectron(0),
fSparseULSElectron(0),
fvaluePHElectron(0),
// MC info
fReadMC(kFALSE),
fNTotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),
// enhance weight cal.
ftype(-1),
fWeight(1),
fCalculateWeight(kFALSE),
fCalculateElectronEffi(kTRUE),
fSprsPi0EtaWeightCal(0),
fPi0Weight(0),
fEtaWeight(0),
fPi0Weight2(0),
fEtaWeight2(0),

fWeightEta(1),
fWeightPi0(1),

// inclusive e-
fInclsElecPt(0),
fInclsElecPtAll(0),
fInclsElecPtReco(0),
// hfe e-
fHFElecPtAll(0),
fHFElecPtReco_wtrkcuts(0),
fHFElecPtReco_wTPCPID(0),
fHFElecPtReco_wtrkCalocuts(0),
fHFElecPtReco_wTPCTOFPID(0),
//  nonhfe MC
fIsFrmPi0(kFALSE),
fIsFrmEta(kFALSE),

fNonHFeTrkPt(0),
fNonHFeWeightTrkPt(0),
fPi0eWeightTrkPt(0),
fEtaeWeightTrkPt(0),
fPi0EtaEleTrkPt(0),
//nonhfe invariant mass
fNonHFePairInvmassLS(0),
fNonHFePairInvmassULS(0),
// non hfe reconstructed
fRecoNonHFeTrkPt(0),
fRecoNonHFeWeightTrkPt(0),
fRecoPi0eWeightTrkPt(0),
fRecoEtaeWeightTrkPt(0),
fRecoPi0EtaEleTrkPt(0),

fRecoLSeTrkPt(0),
fRecoLSeWeightTrkPt(0),
fRecoPi0LSeWeightTrkPt(0),
fRecoEtaLSeWeightTrkPt(0),


fRecoULSeTrkPt(0),
fRecoULSeWeightTrkPt(0),
fRecoPi0ULSeWeightTrkPt(0),
fRecoEtaULSeWeightTrkPt(0)


{
    // constructor
    
    fvaluePHElectron = new Double_t[3];
    fvalueMulti = new Double_t[6];
    for(Int_t i=0; i<2; i++) fMultEstimatorAvg[i]=0;
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}


//_____________________________________________________________________________
AliAnalysisTaskHFETPCTOFMultiplicity::~AliAnalysisTaskHFETPCTOFMultiplicity()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
        
        delete fSparseLSElectron;
        delete fSparseULSElectron;
        
        delete []fvaluePHElectron;
        delete fSparseMulti;
        delete []fvalueMulti;
        delete fSprsPi0EtaWeightCal;
        
    }
    
    for(Int_t i=0; i<2; i++) {
        if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
    }
    delete fListProfiles;
    delete gRandom;
}
//_____________________________________________________________________________
void AliAnalysisTaskHFETPCTOFMultiplicity::Init()
{
    // Initialization
    
    
    if(fDebug > 1) printf("AliAnalysisTaskHFETPCTOFMultiplicity::Init() \n");
    
    
    fListProfiles = new TList();
    fListProfiles->SetOwner();
    TString period[2];
    Int_t nProfiles=2;
    period[0]="LHC16s";
    period[1]="LHC16r";
    
    
    for(Int_t i=0; i<nProfiles; i++){
        if(fMultEstimatorAvg[i]){
            TProfile* hprof=new TProfile(*fMultEstimatorAvg[i]);
            hprof->SetName(Form("ProfileTrkVsZvtx%s\n",period[i].Data()));
            fListProfiles->Add(hprof);
        }
    }
    
    fListProfiles->Add(fPi0Weight);
    fListProfiles->Add(fPi0Weight2);
    fListProfiles->Add(fEtaWeight);
    fListProfiles->Add(fEtaWeight2);
    PostData(2,fListProfiles);
    
    
    return;
}

//_____________________________________________________________________________

void AliAnalysisTaskHFETPCTOFMultiplicity::UserCreateOutputObjects()
{
    
    
    fErr=new TF1("fErr","[0]+[1]*TMath::Erf([2]*x-[3])",0,4);
    fErr->SetParameters(3.89999e-01,3.90001e-01,4.81539e-01,3.53973e+00);
    
    fLandau=new TF1("fLandau","[0]*TMath::Landau(x,[1],[2])",0,4);
    fLandau->SetParameters(1.58523e-02,4.34068e+00,1.02726e+00);
    
    fNewHadFunc=new TF1("fNewHadFunc","gaus(0)+gaus(3)+gaus(6)",0.4,5);
    
    //Double_t param[9] = { 7.21109e+23,-3.35277e+00,3.58273e-01,6.86201e-02,9.60855e-01,4.24254e-02,7.45881e+00,1.10354e+01,2.04162e+00};
    Double_t param[9] ={1.07583e-01,4.74117e-01,5.78342e-02,5.17024e-02,9.94191e-01,7.16580e-02,8.48815e+00,1.06248e+01,1.92413e+00};
    fNewHadFunc->SetParameters(param);
    // fNewHadFunc->Draw();
    
    fPi0Weight = new TF1("fPi0Weight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    
    fPi0Weight2 = new TF1("fPi0Weight2","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    
    fEtaWeight = new TF1("fEtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    
    fEtaWeight2 = new TF1("fEtaWeight2","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    
    fPi0Weight->SetParameters(1.26444e+01,-2.50844e+01,9.14740+00 ,2.23096e-06,1.49298e-01);
    fPi0Weight2->SetParameters(2.07376e+00,-2.42689e+00, 9.23739e-01,4.15003e-01,3.60326e-01);
    
    fEtaWeight->SetParameters( 8.51015e+00,-3.11127e+01, 1.15518e+01, 1.47895e-07, 1.08564e-01);
    fEtaWeight2->SetParameters(6.94009e-01,7.94502e-01, 1.59036e-01,2.80932e+01,1.79031e-01);
    
    printf("\nseed = %u\n\n", gRandom->GetSeed());
    
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    
    fNevents = new TH1F ("fNevents","Number of events",6,-0.5,5.5);
    fNevents->GetYaxis()->SetTitle("counts");
    fNevents->GetXaxis()->SetBinLabel(1,"nEvents Total");
    fNevents->GetXaxis()->SetBinLabel(2,"nEvents With nContributor>2 Trks");
    fNevents->GetXaxis()->SetBinLabel(3,"nEvents with |Zvtx|<10cm cut");
    fNevents->GetXaxis()->SetBinLabel(4,"nEvents with Trigger");
    fNevents->GetXaxis()->SetBinLabel(5,"nEvents with pileup cut");
    fNevents->GetXaxis()->SetBinLabel(6,"nEvents with vertex QA cut");
    
    
    fVtxZ         = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
    fVtxY         = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",300,-15,15);
    fVtxX         = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",300,-15,15);
    fTPCdEdx         = new TH2F("fTPCdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",300,0,15,750,10,160);
    fTPCnsigma         = new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",300,0,15,300,-15,15);
    fTOFnsigma         = new TH2F("fTOFnsig","All Track TOF Nsigma distribution;p (GeV/c);#sigma_{TOF}",500,0,50,200,-10,10);
    fTrkPt         = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",500,0,50);
    fTrketa         = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-0.5,0.5);
    fTrkphi         = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,2*3.141);
    
    //-----------------Incl e Spectrum, inv mass and Hadron Contamination--------------------------------------------
    
    
    fInvmassLS         = new TH1F("fInvmassLS", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 50,0,0.5);
    fInvmassULS         = new TH1F("fInvmassULS", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 50,0,0.5);
    fInvmassLSPt         = new TH2F("fInvmassLSPt", "Inv mass of LS (e,e) vs pT; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,50,0,0.5);
    fInvmassULSPt        = new TH2F("fInvmassULSPt", "Inv mass of ULS (e,e) vs pT; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,50,0,0.5);
    
    fULSElecPt          = new TH1F("fULSElecPt","p_{T} distribution of ULS electrons;p_{T} (GeV/c);counts",500,0,50);
    fLSElecPt         = new TH1F("fLSElecPt","p_{T} distribution of LS electrons;p_{T} (GeV/c);counts",500,0,50);
    fTPCdEdxVsPTOFCut     = new TH2F("fTPCdEdxVsPTOFCut", "fTPCdEdxVsPTOFCut distribution",300,0,15,750,10,160);
    fTPCnSigmaVsPTOFCut    = new TH2F("fTPCnSigmaVsPTOFCut","fTPCnSigmaVsPTOFCut distribution",300,0,15,300,-15,15);
    fTPCnSigmaProtonsVsPTOFCut    = new TH2F("fTPCnSigmaProtonsVsPTOFCut","fTPCnSigmaProtonsVsPTOFCut distribution",300,0,15,300,-15,15);
    fTPCnSigmaKaonsVsPTOFCut    = new TH2F("fTPCnSigmaKaonsVsPTOFCut","fTPCnSigmaKaonsVsPTOFCut distribution",300,0,15,300,-15,15);
    fHadCont_Landau    = new TH2F("fHadCont_Landau","fHadCont_Landau ",300,0,15,300,0,300);
    fHadCont_Err        = new TH2F("fHadCont_Err","fHadCont_Err ",300,0,15,300,0,300);
    fPt_incl_e_Landau    = new TH2F("fPt_incl_e_Landau","fPt_incl_e_Landau ",300,0,15,300,0,300);
    fPt_incl_e_Err    = new TH2F("fPt_incl_e_Err","fPt_incl_e_Err ",300,0,15,300,0,300);
    
    //---------------THnSparse------------
    
    
    Int_t binsls[2]    =          {250, 300};
    Double_t xminls[2]    =    { 0,  0};
    Double_t xmaxls[2]    =    { 25, 300};
    
    fSparseLSElectron     = new THnSparseD ("LSElectron","LSElectron;pT;SPDTracklets;",2 ,binsls,xminls,xmaxls);
    fSparseULSElectron     = new THnSparseD ("ULSElectron","ULSElectron;pT;SPDTracklets;",2 ,binsls,xminls,xmaxls);
    
    Int_t binsm[4]    =          {300,300,300,400};
    Double_t xminm[4]    =    {-15,0,0,0};
    Double_t xmaxm[4]    =    { 15,300,300,400};
    fSparseMulti         = new THnSparseD ("Multiplicity","Multiplicity;zvtx;SPDTracklets_data;Corrected_SPDTracklets;ncharge;",4,binsm,xminm,xmaxm);
    
    
    
    if(fReadMC){
        
        fInclsElecPt         = new TH1F("fInclsElecPt","p_{T} distribution of inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
        
        fInclsElecPtAll     = new TH1F("fInclsElecPtAll","p_{T} distribution of all inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
        
        fInclsElecPtReco     = new TH1F("fInclsElecPtReco","p_{T} distribution of reconstructed inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
        
        
        fHFElecPtAll     = new TH2F("fHFElecPtAll","p_{T} distribution of all HFe ;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fHFElecPtReco_wtrkcuts = new TH2F("fHFElecPtReco_wtrkcuts","p_{T} distribution of HF electrons with track cuts;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fHFElecPtReco_wTPCPID = new TH2F("fHFElecPtReco_wTPCPID","p_{T} distribution of all HF electrons with TPC PID;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fHFElecPtReco_wtrkCalocuts = new TH2F("fHFElecPtReco_wtrkCalocuts","p_{T} distribution of all HF electrons with Track matching;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fHFElecPtReco_wTPCTOFPID = new TH2F("fHFElecPtReco_wTPCTOFPID","p_{T} distribution of all HF electrons with EMCal;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        
        fNonHFeTrkPt = new TH2F("fNonHFeTrkPt","Non-HF electrons from all generators;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fNonHFeWeightTrkPt = new TH2F("fNonHFeWeightTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom with weight + No mom;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fPi0eWeightTrkPt = new TH2F("fPi0eWeightTrkPt","Non-HF electrons from embedded #pi^{0} + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fEtaeWeightTrkPt = new TH2F("fEtaeWeightTrkPt","Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fPi0EtaEleTrkPt  = new TH2F("fPi0EtaEleTrkPt","Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        
        fNonHFePairInvmassLS= new TH1F("fNonHFePairInvmassLS", "Inv mass of LS (e,e) if both e- are Non-HFE; mass(GeV/c^2); counts",  50,0,0.5);
        fNonHFePairInvmassULS = new TH1F("fNonHFePairInvmassULS", "Inv mass of ULS (e,e) if both e- are Non-HFE; mass(GeV/c^2); counts",  50,0,0.5);
        
        
        fRecoNonHFeTrkPt = new TH2F("fRecoNonHFeTrkPt"," Reco Non-HF electrons from all generators;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fRecoNonHFeWeightTrkPt = new TH2F("fRecoNonHFeWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fRecoPi0eWeightTrkPt = new TH2F("fRecoPi0eWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fRecoEtaeWeightTrkPt = new TH2F("fRecoEtaeWeightTrkPt","Reco Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fRecoPi0EtaEleTrkPt  = new TH2F("fRecoPi0EtaEleTrkPt","Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        
        
        fRecoLSeTrkPt = new TH2F("fRecoLSeTrkPt"," Reco LS electrons from all generators;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fRecoLSeWeightTrkPt = new TH2F("fRecoLSeWeightTrkPt","Reco LS electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fRecoPi0LSeWeightTrkPt = new TH2F("fRecoPi0LSeWeightTrkPt","Reco LS electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fRecoEtaLSeWeightTrkPt = new TH2F("fRecoEtaLSeWeightTrkPt","Reco LS electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        
        fRecoULSeTrkPt = new TH2F("fRecoULSeTrkPt"," Reco ULS electrons from all generators;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fRecoULSeWeightTrkPt = new TH2F("fRecoULSeWeightTrkPt","Reco ULS electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fRecoPi0ULSeWeightTrkPt = new TH2F("fRecoPi0ULSeWeightTrkPt","Reco ULS electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        fRecoEtaULSeWeightTrkPt = new TH2F("fRecoEtaULSeWeightTrkPt","Reco ULS electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts;SPDTracklets",250,0,25,300,0,300);
        
        
        
        
        
        
        Int_t binw[4] = {500,3,2,7}; //pT, PDG, HijingOrNot, pi0etaType
        Double_t xminw[4] = {0,0,0,-1};
        Double_t xmaxw[4] = {50,3,2,6};
        
        fSprsPi0EtaWeightCal = new THnSparseD("fSprsPi0EtaWeightCal","Sparse to calculate #pi^{0} and #eta weight;p_{T};PDGID;HijingOrNot;pi0etaType;",4,binw,xminw,xmaxw);
        
        
    }
    
    
    fSparseLSElectron->Sumw2();
    fSparseULSElectron->Sumw2();
    fSparseMulti->Sumw2();
    fTPCdEdxVsPTOFCut->Sumw2();
    fTPCnSigmaVsPTOFCut->Sumw2();
    fHadCont_Landau->Sumw2();
    fHadCont_Err->Sumw2();
    fPt_incl_e_Landau->Sumw2();
    fPt_incl_e_Err->Sumw2();
    fTPCnSigmaProtonsVsPTOFCut->Sumw2();
    fTPCnSigmaKaonsVsPTOFCut->Sumw2();
    
    if(fReadMC){
        
        fHFElecPtAll->Sumw2();
        fHFElecPtReco_wtrkcuts->Sumw2();
        fHFElecPtReco_wTPCPID->Sumw2();
        fHFElecPtReco_wtrkCalocuts->Sumw2();
        fHFElecPtReco_wTPCTOFPID->Sumw2();
        fSprsPi0EtaWeightCal->Sumw2();
        
        fNonHFeTrkPt->Sumw2();
        fNonHFeWeightTrkPt->Sumw2();
        fPi0eWeightTrkPt->Sumw2();
        fEtaeWeightTrkPt->Sumw2();
        fPi0EtaEleTrkPt->Sumw2();
        
        fRecoNonHFeTrkPt->Sumw2();
        fRecoNonHFeWeightTrkPt->Sumw2();
        fRecoPi0eWeightTrkPt->Sumw2();
        fRecoEtaeWeightTrkPt->Sumw2();
        fRecoPi0EtaEleTrkPt->Sumw2();
        
        
        fRecoLSeTrkPt->Sumw2();
        fRecoLSeWeightTrkPt->Sumw2();
        fRecoPi0LSeWeightTrkPt->Sumw2();
        fRecoEtaLSeWeightTrkPt->Sumw2();
        
        fRecoULSeTrkPt->Sumw2();
        fRecoULSeWeightTrkPt->Sumw2();
        fRecoPi0ULSeWeightTrkPt->Sumw2();
        fRecoEtaULSeWeightTrkPt->Sumw2();
        
        
    }
    
    // Output list
    fOutputList->Add(fNevents);
    fOutputList->Add(fVtxZ);
    //fOutputList->Add(fVtxY);
    // fOutputList->Add(fVtxX);
    fOutputList->Add(fTPCdEdx);
    fOutputList->Add(fTPCnsigma);
    fOutputList->Add(fTOFnsigma);
    fOutputList->Add(fHadCont_Landau);
    //fOutputList->Add(fHadCont_Err);
    fOutputList->Add(fPt_incl_e_Landau);
    //fOutputList->Add(fPt_incl_e_Err);
    fOutputList->Add(fTrkPt);
    fOutputList->Add(fTrketa);
    fOutputList->Add(fTrkphi);
    fOutputList->Add(fTPCdEdxVsPTOFCut);
    fOutputList->Add(fTPCnSigmaVsPTOFCut);
    fOutputList->Add(fTPCnSigmaProtonsVsPTOFCut);
    fOutputList->Add(fTPCnSigmaKaonsVsPTOFCut);
    
    fOutputList->Add(fInvmassLS);
    fOutputList->Add(fInvmassULS);
    fOutputList->Add(fInvmassLSPt);
    fOutputList->Add(fInvmassULSPt);
    fOutputList->Add(fULSElecPt);
    fOutputList->Add(fLSElecPt);
    
    
    fOutputList->Add(fSparseLSElectron);
    fOutputList->Add(fSparseULSElectron);
    fOutputList->Add(fSparseMulti);
    
    if(fReadMC){
        fOutputList->Add(fSprsPi0EtaWeightCal);
        fOutputList->Add(fInclsElecPt);
        fOutputList->Add(fInclsElecPtAll);
        fOutputList->Add(fInclsElecPtReco);
        fOutputList->Add(fHFElecPtAll);
        fOutputList->Add(fHFElecPtReco_wtrkcuts);
        fOutputList->Add(fHFElecPtReco_wTPCPID);
        fOutputList->Add(fHFElecPtReco_wtrkCalocuts);
        fOutputList->Add(fHFElecPtReco_wTPCTOFPID);
        
        fOutputList->Add(fNonHFeTrkPt);
        fOutputList->Add(fNonHFeWeightTrkPt);
        fOutputList->Add(fPi0eWeightTrkPt);
        fOutputList->Add(fEtaeWeightTrkPt);
        fOutputList->Add(fPi0EtaEleTrkPt);
        
        fOutputList->Add(fRecoNonHFeTrkPt);
        fOutputList->Add(fRecoNonHFeWeightTrkPt);
        fOutputList->Add(fRecoPi0eWeightTrkPt);
        fOutputList->Add(fRecoEtaeWeightTrkPt);
        fOutputList->Add(fRecoPi0EtaEleTrkPt);
        
        fOutputList->Add(fNonHFePairInvmassLS);
        fOutputList->Add(fNonHFePairInvmassULS);
        
        fOutputList->Add(fRecoLSeTrkPt);
        fOutputList->Add(fRecoLSeWeightTrkPt);
        fOutputList->Add(fRecoPi0LSeWeightTrkPt);
        fOutputList->Add(fRecoEtaLSeWeightTrkPt);
        
        fOutputList->Add(fRecoULSeTrkPt);
        fOutputList->Add(fRecoULSeWeightTrkPt);
        fOutputList->Add(fRecoPi0ULSeWeightTrkPt);
        fOutputList->Add(fRecoEtaULSeWeightTrkPt);
        
        
    }
    
    
    PostData(1,fOutputList);
    PostData(2,fListProfiles);
    
}
//_____________________________________________________________________________
void AliAnalysisTaskHFETPCTOFMultiplicity::UserExec(Option_t *)
{
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;
    
    if(!PassEventSelect(fAOD)) return;
    
    fNevents->Fill(3);
    
    
    //////////////////
    // Multiplicity //
    /////////////////
    
    
    //---------------------multiplicity ------------------------------
    Double_t Zvertex1 = -100;
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Zvertex1 =pVtx->GetZ();
    
    if (fRejectPUFromSPD && fAOD->IsPileupFromSPDInMultBins()) return; // pile-up cut
    fNevents->Fill(4);
    //--------------------vertex selection cuts-----------------------
    AliAODVertex* vtxSPD = fAOD->GetPrimaryVertexSPD();
    
    if (!vtxSPD || vtxSPD->GetNContributors() < 1) return;
    Double_t cov[6]={0};
    vtxSPD->GetCovarianceMatrix(cov);
    if (TMath::Sqrt(cov[5]) > 0.25) return;
    if (TMath::Abs(vtxSPD->GetZ() - pVtx->GetZ())>0.5) return;
    
    fNevents->Fill(5);
    
    //----------V0M Multiplicity------------------
    AliAODVZERO *vzeroAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(fAOD)->GetVZEROData());
    Int_t V0AMult = static_cast<Int_t>(vzeroAOD->GetMTotV0A());
    Int_t V0CMult = static_cast<Int_t>(vzeroAOD->GetMTotV0C());
    Int_t V0Mult=V0AMult+V0CMult;
    
    //------------SPDTracklets--------------------
    Int_t nTracklets = 0;
    Int_t nAcc = 0;
    Double_t etaRange = 1.0;
    
    AliAODTracklets *tracklets = static_cast<const AliAODEvent*>(fAOD)->GetTracklets();
    nTracklets = tracklets->GetNumberOfTracklets();
    for (Int_t nn = 0; nn < nTracklets; nn++) {
        Double_t theta = tracklets->GetTheta(nn);
        Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
        if (TMath::Abs(eta) < etaRange) nAcc++;
    }
    // Data driven multiplicity z-vertex correction
    
    
    Int_t countMult = nAcc;
    
    //SPDTracklets correction
    Double_t correctednAcc   = nAcc;
    Double_t correctednAcc1   = nAcc;
    Double_t countCorr       = countMult;
    TProfile* estimatorAvg = GetEstimatorHistogram(fAOD);
    if(estimatorAvg){
        correctednAcc=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,nAcc,Zvertex1,fRefMult));
        correctednAcc1 = AliAnalysisTaskHFETPCTOFMultiplicity::GetTrackletsMeanCorrection(estimatorAvg,nAcc,Zvertex1,fRefMult);
        countCorr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countMult,Zvertex1,fRefMult));
    }
    
    
    //V0M Correction
    Int_t vzeroMultACorr=V0AMult, vzeroMultCCorr=V0CMult, vzeroMultCorr=V0Mult;
    vzeroMultACorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(V0AMult,Zvertex1));
    vzeroMultCCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0C(V0CMult,Zvertex1));
    vzeroMultCorr = vzeroMultACorr + vzeroMultCCorr; // corrected V0M
    
    ////////////////////////////////
    //Monte Carlo information    //
    ////////////////////////////////
    
    if(fReadMC){
        
        fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        if(!fMCArray){
            AliError("Array of MC particles not found");
            return;
        }
        fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        if (!fMCHeader) {
            AliError("Could not find MC Header in AOD");
            return;
        }
        
        ////////////////////////////////
        //Get number of Gen particles //
        ////////////////////////////////
        GetNMCPartProduced();
        
        
        /////////////////////////
        //Electrons in MC stack//
        /////////////////////////
        
        
        //electrons
        for(Int_t imcArrayL=0; imcArrayL< fMCArray->GetEntries(); imcArrayL++){
            AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCArray->At(imcArrayL);
            Int_t PDGcode = TMath::Abs(AODMCtrack->GetPdgCode());
            
            if(TMath::Abs(AODMCtrack->Eta()) > fCutTrackEta) continue;
            if(!AODMCtrack->IsPhysicalPrimary()) continue;
            
            
            if(PDGcode == 11) { // selecting only electrons
                fInclsElecPtAll->Fill(AODMCtrack->Pt());
                
                Int_t IsHFe=GetHFE(AODMCtrack,fMCArray);
                if((IsHFe==kBeauty) || (IsHFe==kCharm)){
                    fHFElecPtAll->Fill(AODMCtrack->Pt(),correctednAcc1);
                }
            }
        }
        /////////////////////////////////
        //Calculate Pi0 and Eta weight //
        /////////////////////////////////
        
        if(fCalculateWeight) GetPi0EtaWeight(fSprsPi0EtaWeightCal);
    }
    
    fvalueMulti[0] = Zvertex1;
    fvalueMulti[1] = nAcc;
    fvalueMulti[2] = correctednAcc1;
    fvalueMulti[3] = GetNcharged();
    
    
    fSparseMulti->Fill(fvalueMulti);    // multiplicity from tracklets
    
    
    fpidResponse = fInputHandler->GetPIDResponse();
    if(!fpidResponse) return;
    
    
    ///////////////////////
    // Track information //
    //////////////////////
    
    for (Int_t iTracks=0; iTracks<fAOD->GetNumberOfTracks(); iTracks++) {
        
        
        AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTracks);
        if(!track) continue;
        if(!Passtrackcuts(track)) continue;
        
        
        Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999,dEdx = -999, nsigmaTPC = -999.0,nsigmaTOF = -999.0, TOFBeta=-999., nsigmaTPCProtons= -999., nsigmaTPCKaons = -999., nsigmaTOFProtons = -999., nsigmaTOFKaons= -999.;
        Int_t pidM = -1;
        
        TOFBeta=Beta(track);
        Double_t weight_lan=0.,weight_lan_inv=0.;
        Double_t weight_err=0.,weight_err_inv=0.;
        
        nsigmaTPC=fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        nsigmaTPCProtons=fpidResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        nsigmaTPCKaons=fpidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        
        nsigmaTOF = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
        nsigmaTOFProtons = fpidResponse->NumberOfSigmasTOF(track, AliPID::kProton);
        nsigmaTOFKaons = fpidResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
        
        dEdx = track->GetTPCsignal();
        
        TrkPhi = track->Phi();
        fTrkphi->Fill(TrkPhi);
        TrkPt = track->Pt();
        fTrkPt->Fill(TrkPt);
        
        TrkEta = track->Eta();
        fTrketa->Fill(TrkEta);
        TrkP = track->P();
        
        
        
        fTPCdEdx->Fill(TrkP,dEdx);
        fTPCnsigma->Fill(TrkP,nsigmaTPC);
        fTOFnsigma->Fill(TrkP,nsigmaTOF);
        
        if(fReadMC){
            Int_t iTrklabel = TMath::Abs(track->GetLabel());
            if(iTrklabel == 0) continue;
            AliAODMCParticle *MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
            if(TMath::Abs(MCPart->Eta()) > fCutTrackEta) continue;
            if(!MCPart->IsPhysicalPrimary()) continue;
            if(TMath::Abs(MCPart->GetPdgCode())==11) { //only electrons
                Int_t IsElecHf=GetHFE(MCPart,fMCArray);
                if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){ //HF electrons
                    fHFElecPtReco_wtrkcuts -> Fill(TrkPt,correctednAcc1);
                    
                }
            }
        }
        
        if(TMath::Abs(nsigmaTOFProtons) < 3.)
        {  fTPCnSigmaProtonsVsPTOFCut->Fill(TrkP,nsigmaTPC); }
        
        if(TMath::Abs(nsigmaTOFKaons) < 3.)
        {fTPCnSigmaKaonsVsPTOFCut->Fill(TrkP,nsigmaTPC); }
        
        if(TMath::Abs(nsigmaTOF) < fCutNsigmaTOF)
        {
            fTPCdEdxVsPTOFCut->Fill(TrkP,dEdx);
            fTPCnSigmaVsPTOFCut->Fill(TrkP,nsigmaTPC);
            
        }
        weight_lan=fNewHadFunc->Eval(TrkP);
        //weight_err=fErr->Eval(TrkP);
        weight_lan_inv=1-weight_lan;
        //weight_err_inv=1-weight_err;
        
        
        
        if(nsigmaTPC>fCutNsigmaEMin && nsigmaTPC<fCutNsigmaEMax)
        {
            
            if(fReadMC){
                Int_t iTrklabel = TMath::Abs(track->GetLabel());
                if(iTrklabel == 0) continue;
                AliAODMCParticle *MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
                if(TMath::Abs(MCPart->Eta()) > fCutTrackEta) continue;
                if(!MCPart->IsPhysicalPrimary()) continue;
                if(TMath::Abs(MCPart->GetPdgCode())==11) { //only electrons
                    Int_t IsElecHf=GetHFE(MCPart,fMCArray);
                    if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){ //HF electrons
                        
                        
                        fHFElecPtReco_wTPCPID -> Fill(TrkPt,correctednAcc1); }
                    
                }
            }
            
            if(TMath::Abs(nsigmaTOF) < fCutNsigmaTOF)
            {
                
                fHadCont_Landau->Fill(TrkPt,correctednAcc1,weight_lan);
                // fHadCont_Err->Fill(TrkPt,correctednAcc1,weight_err);
                fPt_incl_e_Landau->Fill(TrkPt,correctednAcc1,weight_lan_inv);
                // fPt_incl_e_Err->Fill(TrkPt,correctednAcc1,weight_err_inv);
                
                
                // Heavy-flavour electron reconstruction
                if(fReadMC)
                    
                {Int_t iTrklabel = TMath::Abs(track->GetLabel());
                    if(iTrklabel == 0) continue;
                    AliAODMCParticle *MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
                    if(TMath::Abs(MCPart->Eta()) > fCutTrackEta) continue;
                    if(!MCPart->IsPhysicalPrimary()) continue;
                    if(TMath::Abs(MCPart->GetPdgCode())!=11) continue;
                    Int_t IsElecHf=GetHFE(MCPart,fMCArray);
                    if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){
                        
                        fHFElecPtReco_wTPCTOFPID -> Fill(TrkPt,correctednAcc1); }
                }
                
                
                //////////////////////////////////////
                //Reconst NonHFE with invmass method//
                //////////////////////////////////////
                
                //////////////////////////////////
                //Non-HFE efficiency calculation//
                //////////////////////////////////
                Bool_t EffiDenom = kFALSE;
                Bool_t EffiNumTag = kFALSE;
                if(fReadMC){
                    
                    EffiDenom = GetNonHFEEffiDenomGenPurMC(track,correctednAcc1);
                    
                }
                
                ////////////////////
                //NonHFE selection//
                ////////////////////
                Bool_t fFlagNonHFE=kFALSE;
                fvaluePHElectron[0] = TrkPt;
                fvaluePHElectron[1] = correctednAcc1; //SPD Tracklets
                
                SelectPhotonicElectron(iTracks,track,fFlagNonHFE,pidM, correctednAcc1);
                
                //////////////////////////////////
                //Non-HFE efficiency calculation//
                //////////////////////////////////
                if(fReadMC){
                    if(fFlagNonHFE){
                        EffiNumTag = GetNonHFEEffiRecoTagGenPurMC(track,correctednAcc1);
                    }
                }
                
                
            } // TOF nsigma
            
        } //TPC nsignma
        
        
        
        
    }//track loop
    
    PostData(1,fOutputList);
    PostData(2,fListProfiles);
}

//______________________________________________________________________

Bool_t  AliAnalysisTaskHFETPCTOFMultiplicity::IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromHijing, Int_t &type, Int_t &iMCmom, Int_t &MomPDG, Double_t &MomPt)
{
    //Is electron from pi0, eta and gamma
    
    iMCmom = MCPart->GetMother();
    if(iMCmom<0) return kFALSE;
    else{
        AliAODMCParticle *MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
        MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
        MomPt = MCPartMom->Pt();
        
        if((MomPDG == 111) || (MomPDG == 221) || (MomPDG == 22)){
            if(iMCmom >= fNpureMC)fFromHijing = kFALSE;
            type = GetPi0EtaType(MCPartMom);
            return kTRUE;
        }
        else return kFALSE;
    }
}

//______________________________________________________________________
Bool_t AliAnalysisTaskHFETPCTOFMultiplicity::Passtrackcuts(AliAODTrack *atrack)
{
    Double_t d0z0[2]={-999,-999}, cov[3];
    Double_t DCAxyCut = fCutDCAxy, DCAzCut = fCutDCAz;
    Double_t dEdx =-999;
    Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    
    //kink daughters
    Int_t numberofvertices = 100;
    numberofvertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[numberofvertices];
    Int_t numberofmotherkink = 0;
    for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
        AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
        if(!aodvertex) continue;
        if(aodvertex->GetType()==AliAODVertex::kKink) {
            AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
            if(!mother) continue;
            Int_t idmother = mother->GetID();
            listofmotherkink[numberofmotherkink] = idmother;
            numberofmotherkink++;
        }
    }
    
    //reject kink
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
        if(atrack->GetID() == listofmotherkink[kinkmother]) {
            kinkmotherpass = kFALSE;
            continue;
        }
    }
    if(!kinkmotherpass) return kFALSE;
    
    //other track cuts
    TrkPt = atrack->Pt();
    TrkEta = atrack->Eta();
    if(TrkPt < 0.1 ) return kFALSE;
    Double_t nclusF = atrack->GetTPCNclsF();
    Double_t nclusN = atrack->GetTPCsignalN();  // TPC cluster information findable
    if(nclusF > 0.){
        Double_t RatioTPCclusters = (Double_t)nclusN/nclusF;
        if(RatioTPCclusters<0.6) return kFALSE;
        // cout<<"ratio of tpc cluster"<<RatioTPCclusters<<endl;
    }
    
    
    if (TMath::Abs(TrkEta)>0.7) return kFALSE;
    if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE; //minimum cuts- filter bit 4
    
    if(atrack->GetTPCCrossedRows() < fCutTPCMaxCls) return kFALSE;
    if(atrack->Chi2perNDF() >= fCutTPCchi2perNDF) return kFALSE;
    if(atrack->GetTPCNcls() < fCutTPCNCls) return kFALSE;
    if(atrack->GetITSNcls() < fCutITSNCls) return kFALSE;
    if(nclusN< 80.) return kFALSE;
    
    
    if((!(atrack->GetStatus()&AliAODTrack::kITSrefit)|| (!(atrack->GetStatus()&AliAODTrack::kTPCrefit)))) return kFALSE; //TPC refit
    
    if(fSPDBoth){ if(!(atrack->HasPointOnITSLayer(0) && atrack->HasPointOnITSLayer(1))) return kFALSE;} //Hit on first and second SPD layer
    else if(fSPDAny){ if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) return kFALSE;} //Hit on any layer
    else if(fSPDFirst){ if(!(atrack->HasPointOnITSLayer(0))) return kFALSE;} //Hit on first and second SPD layer

    
    if(atrack->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov))
        if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) return kFALSE;
    
    return kTRUE;
    
}

//______________________________________________________________________
Bool_t AliAnalysisTaskHFETPCTOFMultiplicity::PassEventSelect(AliAODEvent *fAOD)
{
    Int_t ntracks = -999;
    fNevents->Fill(0);
    Double_t Zvertex=-100, Xvertex=-100, Yvertex=-100;
    
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    
    Double_t NContV = pVtx->GetNContributors();
    
    if(NContV<fCutNcontV) return kFALSE;
    fNevents->Fill(1);
    
    Zvertex =pVtx->GetZ();
    Yvertex =pVtx->GetY();
    Xvertex =pVtx->GetX();
    
    if(TMath::Abs(Zvertex)>10.0) return kFALSE;
    fNevents->Fill(2);
    
    
    fVtxZ->Fill(Zvertex);
    fVtxX->Fill(Xvertex);
    fVtxY->Fill(Yvertex);
    
    return kTRUE;
    
}
//______________________________________________________________________
Int_t AliAnalysisTaskHFETPCTOFMultiplicity::GetHFE(AliAODMCParticle *electron, TClonesArray *mcArray)
{
    Int_t motherindex=electron->GetMother(); //Getting Electron Mother
    if(motherindex<0) kNoMother;
    
    else {
        AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherindex);
        Int_t motherpdg = mother->GetPdgCode();
        if ( (motherindex >= fNpureMC) &&(int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
        if ( (motherindex >= fNpureMC) && (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
    }
}

//______________________________________________________________________
void AliAnalysisTaskHFETPCTOFMultiplicity::SelectPhotonicElectron(Int_t itrack, AliAODTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC, Int_t correctednAcc1)
{
    //Photonic electron selection
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t d0z0[2]={-999,-999}, cov[3];
    Double_t DCAxyCut = fCutDCAxy, DCAzCut = fCutDCAz;
    
    Bool_t flagPhotonicElec = kFALSE, flagLSElec = kFALSE;
    Double_t ptAsso=-999., nsigmaAsso=-999.0;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    
    for(Int_t jTracks = 0; jTracks < fAOD->GetNumberOfTracks(); jTracks++){
        if(jTracks==itrack) continue;
        
        
        AliAODTrack* atrackAsso = 0x0;
        atrackAsso = (AliAODTrack*)fAOD->GetTrack(jTracks);
        
        if(!atrackAsso) continue;
        
        if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
        if(atrackAsso->GetTPCCrossedRows() < fCutTPCMaxCls) continue;
        if(atrackAsso->Chi2perNDF() >= fCutTPCchi2perNDF) continue;
        if(atrackAsso->GetTPCNcls() < fCutTPCNCls) continue;
        if(atrackAsso->GetITSNcls() < fCutITSNCls) continue;
        if((!(atrackAsso->GetStatus()&AliAODTrack::kITSrefit)|| (!(atrackAsso->GetStatus()&AliAODTrack::kTPCrefit)))) continue; //refit required
        //  if(!(atrackAsso->HasPointOnITSLayer(0) || atrackAsso->HasPointOnITSLayer(1))) continue; //ITS //kAny
        //if(!(atrackAsso->HasPointOnITSLayer(0) && atrackAsso->HasPointOnITSLayer(1))) continue; //ITS //kBoth
        
        if(fSPDBoth){ if(!(atrackAsso->HasPointOnITSLayer(0) && atrackAsso->HasPointOnITSLayer(1)))  continue;} //Hit on first and second SPD layer
        else if(fSPDAny){ if(!(atrackAsso->HasPointOnITSLayer(0) || atrackAsso->HasPointOnITSLayer(1))) continue;} //Hit on any layer
        else if(fSPDFirst){ if(!(atrackAsso->HasPointOnITSLayer(0))) continue;} //Hit on first and second SPD layer

        
        nsigmaAsso = fpidResponse->NumberOfSigmasTPC(atrackAsso, AliPID::kElectron);
        ptAsso = atrackAsso->Pt();
        Int_t chargeAsso = atrackAsso->Charge();
        Int_t charge = track->Charge();
        
        if(ptAsso < fCutAssoEPt) continue; //partner electron pt 0.3
        if(atrackAsso->Eta()<-fCutAssoEEta || atrackAsso->Eta()>fCutAssoEEta) continue;
        if(nsigmaAsso < -fCutAssoENsigma || nsigmaAsso > fCutAssoENsigma) continue;
        
        
        if(atrackAsso->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov))
            if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
        
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        
        fFlagLS=kFALSE; fFlagULS=kFALSE;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        AliKFParticle::SetField(fAOD->GetMagneticField());
        
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*atrackAsso, fPDGe2);
        AliKFParticle recg(ge1,ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        Double_t mass=-999., width = -999.;
        Int_t MassCorrect;
        MassCorrect = recg.GetMass(mass,width);
        
        if(fFlagLS && track->Pt()>1) fInvmassLS->Fill(mass);
        if(fFlagULS && track->Pt()>1) fInvmassULS->Fill(mass);
        if(fFlagLS) fInvmassLSPt->Fill(track->Pt(),mass);
        if(fFlagULS) fInvmassULSPt->Fill(track->Pt(),mass);
        
        //////////////////////////////////
        //Non-HFE efficiency calculation//
        //////////////////////////////////
        Bool_t EffiNumULSLS = kFALSE;
        
        if(fReadMC){
            EffiNumULSLS = GetNonHFEEffiULSLSGenPurMC(track, atrackAsso, fFlagLS, fFlagULS, mass,correctednAcc1);
        }
        Double_t TrkPt = track->Pt();
        if(mass < fCutInvmass){
            if(fFlagLS){
                fLSElecPt->Fill(track->Pt()); //Reco LS e TrkPt
                fSparseLSElectron->Fill(fvaluePHElectron); }
            
            
            if(fFlagULS){
                fULSElecPt->Fill(track->Pt());  //Reco ULS e TrkPt
                fSparseULSElectron->Fill(fvaluePHElectron);
                
            }
        }
        
        
        if(mass<fCutInvmass && fFlagULS && !flagPhotonicElec){
            flagPhotonicElec = kTRUE;
        }
        
    }
    fFlagPhotonicElec = flagPhotonicElec;
    
    
}
//____________________________________________________________________________
TProfile* AliAnalysisTaskHFETPCTOFMultiplicity::GetEstimatorHistogram(const AliAODEvent* fAOD)
{
    
    Int_t runNo  = fAOD->GetRunNumber();
    //cout<<"run number"<<runNo<<endl;
    Int_t period = -1;
    
    
    if (runNo>266436 && runNo<267111) period = 0;
    if (runNo>265593 && runNo<266319) period = 1;
    if (period < 0 || period > 1) return 0;
    
    //cout<<"period ="<<period<<endl;
    
    return fMultEstimatorAvg[period];
}
//____________________________________________________________________________
Double_t AliAnalysisTaskHFETPCTOFMultiplicity::Beta(AliAODTrack *track)
{
    Double_t stoptime=track->GetTOFsignal();
    
    Double_t c=TMath::C()*1.E-9;// m/ns
    Float_t startTime= fpidResponse->GetTOFResponse().GetStartTime(((AliAODTrack*)track)->P());//in ps
    Double_t length= fpidResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron)*1E-3*c;
    stoptime -= startTime;
    Double_t scaleStopTime= stoptime*1E-3;
    scaleStopTime=scaleStopTime*c;
    return length/scaleStopTime;
}
//____________________________________________________________________________


//______________________________________________________________________________
Int_t AliAnalysisTaskHFETPCTOFMultiplicity::GetNcharged(){
    //counts all tracks in eta<1 with charge!=0
    
    Int_t Nch = 0;
    
    if(!fReadMC) return Nch; // if no MC info return 0
    
    // loop over all tracks
    for (Int_t igen = 0; igen < fMCArray->GetEntriesFast(); igen++){
        AliAODMCParticle *mctrack=(AliAODMCParticle*)fMCArray->UncheckedAt(igen);
        Int_t charge = mctrack->Charge();
        Double_t eta = mctrack->Eta();
        Bool_t isPhysPrim = mctrack->IsPhysicalPrimary();
        if(charge!=0){
            if(eta > -1.0 && eta < 1.0){
                if(isPhysPrim){
                    Nch++;
                }
            }
        }
    }
    return Nch;
}
//______________________________________________________________________________
void AliAnalysisTaskHFETPCTOFMultiplicity::GetPi0EtaWeight(THnSparse *SparseWeight)
{
    //Get pi0 and eta information for weight calculation
    
    Double_t fvalue[4] = {-999,-999,-999,-999};
    
    for(int imc=0; imc< fNTotMCpart; imc++)
    {
        AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCArray->At(imc);
        if(TMath::Abs(AODMCtrack->Eta()) > 1.2) continue;
        
        //-------Get PDG
        Int_t TrackPDG = TMath::Abs(AODMCtrack->GetPdgCode());
        if((TrackPDG != 111) && (TrackPDG != 221) && (TrackPDG != 22)) continue;
        
        Double_t fPartPDGid = -999;
        if (TrackPDG == 111) fPartPDGid = 0.2;
        if (TrackPDG == 221) fPartPDGid = 1.2;
        if (TrackPDG == 22) fPartPDGid = 2.2;
        
        Double_t fTrkPt = AODMCtrack->Pt();
        
        //-------Check if the particle is from hijing or not
        Bool_t fFromHijing = kHijing;
        if(imc >= fNpureMC)fFromHijing = kElse;
        
        //------Get type of the particle
        Int_t fType = GetPi0EtaType(AODMCtrack);
        
        fvalue[0] = fTrkPt;
        fvalue[1] = fPartPDGid;
        fvalue[2] = fFromHijing;
        fvalue[3] = fType;
        
        SparseWeight->Fill(fvalue);
    }
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFETPCTOFMultiplicity::GetPi0EtaType(AliAODMCParticle *part)
{
    // Return the type of particle
    
    // IsPrimary
    Bool_t primMC = part->IsPrimary();
    if(!primMC) return kNotIsPrimary;
    
    // Mother
    Int_t motherlabel = part->GetMother();
    if(motherlabel<0) return kNoMother;
    
    else {
        AliAODMCParticle *mother = (AliAODMCParticle*)fMCArray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        
        if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213 || motherpdg == 313 || motherpdg == 323) return kLightMesons;
        
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
        return kNoFeedDown;
    }
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskHFETPCTOFMultiplicity::GetNMCPartProduced()
{
    //Get number of MC particles produced by generators.
    
    TList *lh = fMCHeader->GetCocktailHeaders();
    fNTotMCpart = 0;
    fNembMCpi0 = 0;
    fNembMCeta = 0;
    fNpureMC = 0;
    TString MCgen;
    TString embpi0("pi");
    TString embeta("eta");
    
    if(!lh){
        AliError("no MC header");
        return (0);
    }
    
    for(int igene=0; igene<lh->GetEntries(); igene++)
    {
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
        if(!gh) continue;
        
        MCgen =  gh->GetName();
        
        if(igene==0) fNpureMC = gh->NProduced();  // generated by HIJING
        
        
        if(MCgen.Contains(embpi0))fNembMCpi0 = fNTotMCpart;
        if(MCgen.Contains(embeta))fNembMCeta = fNTotMCpart;
        fNTotMCpart += gh->NProduced();
        //cout << "fNpureMC, fNembMCpi0, fNembMCeta, fNTotMCpart : " <<fNpureMC << ", " << fNembMCpi0 << ", " << fNembMCeta << ", " << fNTotMCpart << endl;
    }
    
    
    return kTRUE;
}
//______________________________________________________________________________

Bool_t AliAnalysisTaskHFETPCTOFMultiplicity::GetNonHFEEffiDenomGenPurMC(AliVTrack *track,Int_t correctednAcc1)
{
    //Calculate Non-HFE efficiency demoninator for General purpose MC
    fIsFrmPi0 = kFALSE, fIsFrmEta = kFALSE;
    ftype = -1, fWeightPi0 = 1.0, fWeightEta = 1.0, fWeight=1.0;
    Bool_t fFromMB = kTRUE;
    
    Int_t MomPDG = -999, GMomPDG=-999, GGMomPDG=-999, GGGMomPDG=-999;
    Int_t iMCmom = -999, iMCgmom = -999, iMCggmom = -999, iMCgggmom = -999;
    Double_t MomPt =-999.0, ptmotherw = -999.0, ptgmotherw= -999.0, ptggmotherw = -999.0;
    
    AliAODMCParticle *MCPart = 0;
    AliAODMCParticle *MCPartMom = 0;
    AliAODMCParticle *MCPartGMom = 0;
    AliAODMCParticle *MCPartGGMom = 0;
    AliAODMCParticle *MCPartGGGMom = 0;
    
    Double_t TrkPt = track->Pt();
    Int_t iTrklabel = TMath::Abs(track->GetLabel());
    if(iTrklabel == 0) return kFALSE;
    
    MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
    if(TMath::Abs(MCPart->GetPdgCode())!=11) return kFALSE;
    fInclsElecPt->Fill(TrkPt);
    
    
    Bool_t fNonHFE = IsNonHFE(MCPart, fFromMB, ftype, iMCmom, MomPDG, MomPt);
    if(!fNonHFE) return kFALSE;
    fNonHFeTrkPt->Fill(TrkPt, correctednAcc1);
    
    MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
    ptmotherw = MCPartMom->Pt();
    iMCgmom = MCPartMom->GetMother();
    
    if(iMCgmom > 0){
        MCPartGMom = (AliAODMCParticle*)fMCArray->At(iMCgmom);
        GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());
        ptgmotherw = MCPartGMom->Pt();
        
        iMCggmom = MCPartGMom->GetMother();
        if(iMCggmom > 0){
            MCPartGGMom = (AliAODMCParticle*)fMCArray->At(iMCggmom);
            GGMomPDG = TMath::Abs(MCPartGGMom->GetPdgCode());
            ptggmotherw = MCPartGGMom->Pt();
            
            iMCgggmom = MCPartGGMom->GetMother();
            if(iMCgggmom > 0){
                MCPartGGGMom = (AliAODMCParticle*)fMCArray->At(iMCgggmom);
                GGGMomPDG = TMath::Abs(MCPartGGGMom->GetPdgCode());
            }
        }
    }
    
    //cases to consider: eta->e, eta->pi0->e, eta->gamma->e, eta->pi0->gamma->e, pi0->e, pi0->gamma->e
    if(MomPDG == 221){
        fIsFrmEta = kTRUE; //eta->e
        fWeightEta = fEtaWeight->Eval(MCPartMom->Pt());
        if(ptmotherw >= 0.4 && ptmotherw <= 1.45) fWeightEta = 1./(fEtaWeight->Eval(ptmotherw));
        if(ptmotherw > 1.45 && ptmotherw <=15) fWeightEta = 1./(fEtaWeight2->Eval(ptmotherw));
        
        
    }
    
    if(MomPDG == 111) {
        fIsFrmPi0 = kTRUE; //pi0 -> e
        if(ptmotherw >= 0.3 && ptmotherw <= 1.3) fWeightPi0 =1./( fPi0Weight->Eval(ptmotherw));
        if(ptmotherw > 1.3 && ptmotherw <=15) fWeightPi0 =1./( fPi0Weight2->Eval(ptmotherw));
        
        if(GMomPDG == 221){
            fIsFrmEta = kTRUE; //eta->pi0-> e
            if(ptgmotherw >= 0.4 && ptgmotherw <= 1.45) fWeightEta =1./( fEtaWeight->Eval(ptgmotherw));
            if(ptgmotherw > 1.45 && ptgmotherw <=15) fWeightEta = 1./(fEtaWeight2->Eval(ptgmotherw));
        }
    }
    
    if(MomPDG == 22){
        if(GMomPDG == 221){
            fIsFrmEta = kTRUE; //eta->gamma-> e
            if(ptgmotherw >= 0.4 && ptgmotherw < 1.45) fWeightEta = 1./(fEtaWeight->Eval(ptgmotherw));
            if(ptgmotherw > 1.45 && ptgmotherw <= 15) fWeightEta = 1./(fEtaWeight2->Eval(ptgmotherw));
            
        }
        
        if(GMomPDG == 111){
            fIsFrmPi0 = kTRUE; //pi0-> gamma-> e
            if(ptgmotherw >= 0.3 && ptgmotherw <= 1.3) fWeightPi0 = 1./(fPi0Weight->Eval(ptgmotherw));
            if(ptgmotherw > 1.3 && ptgmotherw <= 15) fWeightPi0 = 1./(fPi0Weight2->Eval(ptgmotherw));
            
            if(GGMomPDG == 221){
                fIsFrmEta = kTRUE; //eta->pi0->gamma-> e
                if(ptggmotherw >= 0.4 && ptggmotherw < 1.45) fWeightEta = 1./(fEtaWeight->Eval(ptggmotherw));
                if(ptggmotherw > 1.45 && ptggmotherw <= 15) fWeightEta = 1./(fEtaWeight2->Eval(ptggmotherw));
            }
        }
    }
    
    //   cout << "PDG of M, GM, GGM, GGGM of ele: "<< MomPDG << ", " << GMomPDG << ", " << GGMomPDG << ", " << GGGMomPDG << endl;
    
    if(fIsFrmPi0 || fIsFrmEta){
        fPi0EtaEleTrkPt->Fill(TrkPt, correctednAcc1);
        
        if(fIsFrmPi0) {
            fWeight = fWeightPi0;
            
            fPi0eWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightPi0);
            fNonHFeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightPi0);
        }
        if(fIsFrmEta){
            fWeight = fWeightEta;
            fEtaeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightEta);
            fNonHFeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightEta);
        }
    }
    
    return kTRUE;
}

//______________________________________________________________________________

Bool_t AliAnalysisTaskHFETPCTOFMultiplicity::GetNonHFEEffiRecoTagGenPurMC(AliVTrack *track, Int_t correctednAcc1)
{
    //Tagging method for General purpose MC
    
    Double_t TrkPt = track->Pt();
    
    fRecoNonHFeTrkPt->Fill(TrkPt, correctednAcc1);
    if(fIsFrmPi0 || fIsFrmEta){
        fRecoPi0EtaEleTrkPt->Fill(TrkPt,correctednAcc1);
        
        if(fIsFrmPi0) {
            fRecoPi0eWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightPi0);
            fRecoNonHFeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightPi0);
        }
        if(fIsFrmEta){
            fRecoEtaeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightEta);
            fRecoNonHFeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightEta);
        }
    }
    
    return kTRUE;
}
//______________________________________________________________________________

Bool_t AliAnalysisTaskHFETPCTOFMultiplicity::GetNonHFEEffiULSLSGenPurMC(AliAODTrack *track, AliAODTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass, Int_t correctednAcc1)
{
    //ULS-LS method for General purpose MC
    
    Double_t TrkPt = track->Pt();
    
    //Track information
    Int_t iTrklabel = TMath::Abs(track->GetLabel());
    if(iTrklabel == 0) return kFALSE;
    AliAODMCParticle *MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
    
    if(TMath::Abs(MCPart->GetPdgCode())!=11) return kFALSE;
    Bool_t fFromMB = kTRUE;
    Int_t iMCmom=-999, MomPDG = -999, type=-1;
    Double_t MomPt =-999;
    Bool_t fNonHFE = IsNonHFE(MCPart, fFromMB, type, iMCmom, MomPDG, MomPt);
    
    //Associated partner information
    Int_t iTrkAssolabel = TMath::Abs(Assotrack->GetLabel());
    if(iTrkAssolabel == 0) return kFALSE;
    AliAODMCParticle *MCPartAsso = (AliAODMCParticle*)fMCArray->At(iTrkAssolabel);
    
    if(TMath::Abs(MCPartAsso->GetPdgCode())!=11) return kFALSE; // check origin of asso elec
    
    Bool_t fAssoFromMB = kTRUE;
    Int_t iMCAssomom=-999, AssoMomPDG = -999, fAssotype=-1;
    Double_t AssoMomPt =-999;
    Bool_t fAssoNonHFE = IsNonHFE(MCPartAsso, fAssoFromMB, fAssotype, iMCAssomom, AssoMomPDG, AssoMomPt);
    
    //cout << "Asso ele mom : " << iMCAssomom << ", " << AssoMomPDG << ", " << iMCmom << ", " << MomPDG << ", " << fIsFrmPi0 << ", " << fIsFrmEta << ", " << type << endl;
    
    if(!fAssoNonHFE) return kFALSE;
    if(iMCmom != iMCAssomom) return kFALSE; //ensure electron and partner comes from same mother
    
    if(fFlagLS) fNonHFePairInvmassLS->Fill(mass);
    if(fFlagULS) fNonHFePairInvmassULS->Fill(mass);
    
    if(mass < fCutInvmass){
        if(fFlagLS){
            //new method
            if(fIsFrmPi0 || fIsFrmEta) {
                fRecoLSeTrkPt->Fill(TrkPt, correctednAcc1);
                
                if(fIsFrmPi0) {
                    fRecoLSeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightPi0);
                    fRecoPi0LSeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightPi0);
                }
                if(fIsFrmEta){
                    fRecoLSeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightEta);
                    fRecoEtaLSeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightEta);
                }
            }
            
        }
        
        if(fFlagULS){
            //new method
            if(fIsFrmPi0 || fIsFrmEta) {
                fRecoULSeTrkPt->Fill(TrkPt, correctednAcc1);
                
                if(fIsFrmPi0) {
                    fRecoULSeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightPi0);
                    fRecoPi0ULSeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightPi0);
                }
                
                if(fIsFrmEta){
                    fRecoULSeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightEta);
                    fRecoEtaULSeWeightTrkPt->Fill(TrkPt,correctednAcc1,fWeightEta);
                }
            }
            
        }
    }
    
    return kTRUE;
}
//______________________________________________________________________________
Double_t AliAnalysisTaskHFETPCTOFMultiplicity::GetTrackletsMeanCorrection(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult)
{
    if(TMath::Abs(vtxZ)>10.0){
        //    printf("ERROR: Z vertex out of range for correction of multiplicity\n");
        return uncorrectedNacc;
    }
    
    if(!estimatorAvg){
        printf("ERROR: Missing TProfile for correction of multiplicity\n");
        return uncorrectedNacc;
    }
    
    Double_t localAvg = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxZ));
    Double_t deltaM = 0;
    deltaM = uncorrectedNacc*(refMult/localAvg - 1);
    
    Double_t correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->PoissonD(TMath::Abs(deltaM));
    
    if(correctedNacc<0) correctedNacc=0;
    
    return correctedNacc;
    
}


//______________________________________________________________________________

void AliAnalysisTaskHFETPCTOFMultiplicity::Terminate(Option_t *)
{
    
}
//_________________________________________________________________________


