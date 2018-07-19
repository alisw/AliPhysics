//
// HFE tagged jet shape analysis task.
//
// Authors: D. Caffarri, L. Cunqueiro (jet); D. Godoy (HFE)
//just testing
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "TVector2.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliEmcalPythiaInfo.h"
#include "TRandom3.h"
#include "AliAODInputHandler.h"
#include "AliPIDResponse.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliMultiEventInputHandler.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"

#include "AliAODEvent.h"
#include "AliAnalysisTaskEmcalHfeTagging.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalHfeTagging)

//________________________________________________________________________
AliAnalysisTaskEmcalHfeTagging::AliAnalysisTaskEmcalHfeTagging() :
AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalHfeTagging", kTRUE),
fAOD(0),
fVevent(0),
fpidResponse(0),
fTracksTender(0),
pVtx(0),
spdVtx(0),
fMC(0),
fStack(0),
fMCparticle(0),
fMCarray(0),
fMCheader(0),
fContainer(0),
fMinFractionShared(0),
fJetShapeType(kData),
fJetShapeSub(kNoSub),
fJetSelection(kInclusive),
fPtThreshold(-9999.),
fRMatching(0.2),
fSelectedShapes(0),
fminpTTrig(20.),
fmaxpTTrig(50.),
fangWindowRecoil(0.6),
fSemigoodCorrect(0),
fHolePos(0),
fHoleWidth(0),
fCentSelectOn(kTRUE),
fCentMin(0),
fCentMax(10),
fOneConstSelectOn(kFALSE),
fDerivSubtrOrder(0),
fMCweight(0),
fAssPtCut(0.1),
fITSncut(3),
fAssTPCnCut(60),
fTPCnCut(100),
fAssITSrefitCut(kTRUE),
fUseTender(kFALSE),
fSigmaTOFcut(3.),
fSigmaTPCcut(-1.),
fDcaXYcut(1.),
fDcaZcut(2.),
fIMcut(0.1),
fh2ResponseUW(0x0),
fh2ResponseW(0x0),
fPhiJetCorr6(0x0),
fPhiJetCorr7(0x0),
fEtaJetCorr6(0x0),
fEtaJetCorr7(0x0),
fPtJetCorr(0x0),
fPtJet(0x0),
fPtGenJet(0),
fPhiJet(0x0),
fEtaJet(0x0),
fJetEfficiency(0x0),
fhpTjetpT(0x0),
fhPt(0x0),
fhPhi(0x0),
fNbOfConstvspT(0x0),
fnTPCnTOFnocut(0),
fnTPCnocut(0),
fnTOFnocut(0),
fnTPCcut(0),
fnULSmLSpairsPerElectron(0),
fnPartPerJet(0),
fnElecOverPartPerJet(0),
fnInclElecPerJet(0),
fnPhotElecPerJet(0),
fnIncSubPhotElecPerJet(0),
fnTrueElecPerJet(0),
fnTrueHFElecPerJet(0),
fnTruePElecPerJet(0),
fPi0Pt(0),
fEtaPt(0),
fGenHfePt(0),
fGenPePt(0),
fPtP(0),
fptJetIE(0),
fptJetPE(0),
fptJetHFE(0),
fptRecPE(0),
fptTruePE(0),
fptWrongPE(0),
fPhiTrack(0x0),
fEtaTrack(0x0),
fPhiElec(0x0),
fEtaElec(0x0),
fTreeObservableTagging(0)

{
    for(Int_t i=0;i<21;i++){
        fShapesVar[i]=0;
    }
    
    for(Int_t i=0;i<4;i++){
        fptTrueHFE[i] = NULL;
    }
    
    for(Int_t i=0;i<5;i++){
        fInvmassLS[i] = NULL;
        fInvmassULS[i] = NULL;
        fUlsLsElecPt[i] = NULL;
        fTotElecPt[i] = NULL;
        for(Int_t j=0;j<18;j++){
            fnTPCTrueParticles[i][j] = NULL;
        }
    }
    SetMakeGeneralHistograms(kTRUE);
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    
}

//________________________________________________________________________
AliAnalysisTaskEmcalHfeTagging::AliAnalysisTaskEmcalHfeTagging(const char *name) :
AliAnalysisTaskEmcalJet(name, kTRUE),
fAOD(0),
fVevent(0),
fpidResponse(0),
fTracksTender(0),
pVtx(0),
spdVtx(0),
fMC(0),
fStack(0),
fMCparticle(0),
fMCarray(0),
fMCheader(0),
fContainer(0),
fMinFractionShared(0),
fJetShapeType(kData),
fJetShapeSub(kNoSub),
fJetSelection(kInclusive),
fPtThreshold(-9999.),
fRMatching(0.2),
fSelectedShapes(0),
fminpTTrig(20.),
fmaxpTTrig(50.),
fangWindowRecoil(0.6),
fSemigoodCorrect(0),
fHolePos(0),
fHoleWidth(0),
fCentSelectOn(kTRUE),
fCentMin(0),
fCentMax(10),
fOneConstSelectOn(kFALSE),
fDerivSubtrOrder(0),
fMCweight(0),
fAssPtCut(0.1),
fITSncut(3),
fAssTPCnCut(60),
fTPCnCut(100),
fAssITSrefitCut(kTRUE),
fUseTender(kFALSE),
fSigmaTOFcut(3.),
fSigmaTPCcut(-1.),
fDcaXYcut(1.),
fDcaZcut(2.),
fIMcut(0.1),
fh2ResponseUW(0x0),
fh2ResponseW(0x0),
fPhiJetCorr6(0x0),
fPhiJetCorr7(0x0),
fEtaJetCorr6(0x0),
fEtaJetCorr7(0x0),
fPtJetCorr(0x0),
fPtJet(0x0),
fPtGenJet(0),
fPhiJet(0x0),
fEtaJet(0x0),
fJetEfficiency(0x0),
fhpTjetpT(0x0),
fhPt(0x0),
fhPhi(0x0),
fNbOfConstvspT(0x0),
fnTPCnTOFnocut(0),
fnTPCnocut(0),
fnTOFnocut(0),
fnTPCcut(0),
fnULSmLSpairsPerElectron(0),
fnPartPerJet(0),
fnElecOverPartPerJet(0),
fnInclElecPerJet(0),
fnPhotElecPerJet(0),
fnIncSubPhotElecPerJet(0),
fnTrueElecPerJet(0),
fnTrueHFElecPerJet(0),
fnTruePElecPerJet(0),
fPi0Pt(0),
fEtaPt(0),
fGenHfePt(0),
fGenPePt(0),
fPtP(0),
fptJetIE(0),
fptJetPE(0),
fptJetHFE(0),
fptRecPE(0),
fptTruePE(0),
fptWrongPE(0),
fPhiTrack(0x0),
fEtaTrack(0x0),
fPhiElec(0x0),
fEtaElec(0x0),
fTreeObservableTagging(0)

{
    // Standard constructor.
    for(Int_t i=0;i<21;i++){
        fShapesVar[i]=0;
    }
    
    for(Int_t i=0;i<4;i++){
        fptTrueHFE[i] = NULL;
    }
    
    for(Int_t i=0;i<5;i++){
        fInvmassLS[i] = NULL;
        fInvmassULS[i] = NULL;
        fUlsLsElecPt[i] = NULL;
        fTotElecPt[i] = NULL;
        for(Int_t j=0;j<18;j++){
            fnTPCTrueParticles[i][j] = NULL;
        }
    }
    SetMakeGeneralHistograms(kTRUE);
    
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    
}

//________________________________________________________________________
AliAnalysisTaskEmcalHfeTagging::~AliAnalysisTaskEmcalHfeTagging()
{
    // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHfeTagging::UserCreateOutputObjects()
{
    // Create user output.
    
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    
    fh2ResponseUW= new TH2F("fh2ResponseUW", "fh2ResponseUW", 100, 0, 200,  100, 0, 200);
    fOutput->Add(fh2ResponseUW);
    fh2ResponseW= new TH2F("fh2ResponseW", "fh2ResponseW", 100, 0, 200,  100, 0, 200);
    fOutput->Add(fh2ResponseW);
    fPhiJetCorr6= new TH2F("fPhiJetCorr6", "fPhiJetCorr6", 50, 0, 2*TMath::Pi(), 50, 0, 2*TMath::Pi());
    fOutput->Add(fPhiJetCorr6);
    fEtaJetCorr6= new TH2F("fEtaJetCorr6", "fEtaJetCorr6", 50, -1.5, 1.5, 50, -1.5, 1.5);
    fOutput->Add(fEtaJetCorr6);
    
    fPhiJetCorr7= new TH2F("fPhiJetCorr7", "fPhiJetCorr7", 50, 0, 2*TMath::Pi(), 50, 0, 2*TMath::Pi());
    fOutput->Add(fPhiJetCorr7);
    fEtaJetCorr7= new TH2F("fEtaJetCorr7", "fEtaJetCorr7", 50, -1.5, 1.5, 50, -1.5, 1.5);
    fOutput->Add(fEtaJetCorr7);
    
    fPtJetCorr= new TH2F("fPtJetCorr", "fPtJetCorr", 100, 0, 200,  100, 0, 200);
    fOutput->Add(fPtJetCorr);
    
    fPtJet= new TH1F("fPtJet", "fPtJet", 100, 0, 200);
    fOutput->Add(fPtJet);
    
    fPtGenJet= new TH1F("fPtGenJet", "fPtGenJet", 100, 0, 200);
    fOutput->Add(fPtGenJet);
    
    fPhiJet= new TH2F("fPhiJet", "fPhiJet", 100, 0, 200, 100, 0, TMath::TwoPi());
    fOutput->Add(fPhiJet);
    
    fEtaJet= new TH2F("fEtaJet", "fEtaJet", 100, 0, 200, 100, -1,1);
    fOutput->Add(fEtaJet);
    
    fhpTjetpT= new TH2F("fhpTjetpT", "fhpTjetpT", 200, 0, 200,  200, 0, 200);
    fOutput->Add(fhpTjetpT);
    
    fhPt= new TH1F("fhPt", "fhPt", 200, 0, 200);
    fOutput->Add(fhPt);
    
    fhPhi= new TH1F("fhPhi", "fhPhi", 100, -TMath::Pi(), TMath::Pi());
    fOutput->Add(fhPhi);
    
    fNbOfConstvspT=new TH2F("fNbOfConstvspT", "fNbOfConstvspT", 100, 0, 100, 200, 0, 200);
    fOutput->Add(fNbOfConstvspT);
    
    fnTPCnTOFnocut=new TH2F("fnTPCnTOFnocut", "fnTPCnTOFnocut", 200, -10, 10, 200, -10, 10);
    fOutput->Add(fnTPCnTOFnocut);
    
    fnTPCnocut=new TH2F("fnTPCnocut", "fnTPCnocut", 50, 0, 5, 200, -10, 10);
    fOutput->Add(fnTPCnocut);
    
    fnTOFnocut=new TH2F("fnTOFnocut", "fnTOFnocut", 50, 0, 5, 200, -10, 10);
    fOutput->Add(fnTOFnocut);
    
    fnTPCcut=new TH2F("fnTPCcut", "fnTPCcut", 50, 0, 5, 200, -10, 10);
    fOutput->Add(fnTPCcut);
    
    fnULSmLSpairsPerElectron =new TH2F("fnULSmLSpairsPerElectron", "fnULSmLSpairsPerElectron", 50, 0, 5, 100, 0, 100);
    fOutput->Add(fnULSmLSpairsPerElectron);
    
    fnPartPerJet=new TH1F("fnPartPerJet", "fnPartPerJet", 500,0,500);
    fOutput->Add(fnPartPerJet);
    
    fnElecOverPartPerJet=new TH1F("fnElecOverPartPerJet", "fnElecOverPartPerJet", 100,0,0.1);
    fOutput->Add(fnElecOverPartPerJet);
    
    fnInclElecPerJet=new TH1F("fnInclElecPerJet", "fnInclElecPerJet", 100,0,100);
    fOutput->Add(fnInclElecPerJet);
    
    fnPhotElecPerJet=new TH1F("fnPhotElecPerJet", "fnPhotElecPerJet", 100,0,100);
    fOutput->Add(fnPhotElecPerJet);
    
    fnIncSubPhotElecPerJet=new TH1F("fnIncSubPhotElecPerJet", "fnIncSubPhotElecPerJet", 101,-1,100);
    fOutput->Add(fnIncSubPhotElecPerJet);
    
    fnTrueElecPerJet=new TH1F("fnTrueElecPerJet", "fnTrueElecPerJet", 100,0,100);
    fOutput->Add(fnTrueElecPerJet);
    
    fnTrueHFElecPerJet=new TH1F("fnTrueHFElecPerJet", "fnTrueHFElecPerJet", 100,0,100);
    fOutput->Add(fnTrueHFElecPerJet);
    
    fnTruePElecPerJet=new TH1F("fnTruePElecPerJet", "fnTruePElecPerJet", 100,0,100);
    fOutput->Add(fnTruePElecPerJet);
    
    int nbin = 59;
    double xbins[60] =  {0.01,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,
        0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.5,5,
        5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20};
    
    fPi0Pt = new TH2F("fPi0Pt","fPi0Pt",4,0,4,nbin,xbins);
    fOutput->Add(fPi0Pt);
    
    fEtaPt = new TH2F("fEtaPt", "fEtaPt",4,0,4,nbin,xbins);
    fOutput->Add(fEtaPt);
    
    fGenHfePt = new TH1F("fGenHfePt","fGenHfePt",nbin,xbins);
    fOutput->Add(fGenHfePt);
    
    fGenPePt = new TH1F("fGenPePt", "fGenPePt",nbin,xbins);
    fOutput->Add(fGenPePt);
    
    Double_t ptRange[19] = {0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,5};
    
    fPtP=new TH2F("fPtP", "fPtP", 18,ptRange,18,ptRange);
    fOutput->Add(fPtP);
    
    for(Int_t i=0;i<5;i++){
        
        fInvmassLS[i] = new TH2F(Form("fInvmassLS%d",i), Form("fInvmassLS%d",i), 50, 0, 5, 100, 0, 0.5);
        fOutput->Add(fInvmassLS[i]);
        
        fInvmassULS[i] = new TH2F(Form("fInvmassULS%d",i), Form("fInvmassULS%d",i), 50, 0, 5, 100, 0, 0.5);
        fOutput->Add(fInvmassULS[i]);
        
        fUlsLsElecPt[i] = new TH1F(Form("fUlsLsElecPt%d",i), Form("fUlsLsElecPt%d",i), 18,ptRange);
        fOutput->Add(fUlsLsElecPt[i]);
        
        fTotElecPt[i] = new TH1F(Form("fTotElecPt%d",i), Form("fTotElecPt%d",i), 18,ptRange);
        fOutput->Add(fTotElecPt[i]);
        
        
        for(Int_t j=0;j<18;j++){
            fnTPCTrueParticles[i][j] = new TH2F(Form("fnTPCTrueParticles%d%d",i,j),Form("fnTPCTrueParticles%d%d",i,j), 7,0,7,100,-15,15);
            fOutput->Add(fnTPCTrueParticles[i][j]);
        }
    }
    
    for(Int_t i=0;i<4;i++){
        fptTrueHFE[i] = new TH1F(Form("fptTrueHFE%d",i), Form("fptTrueHFE%d",i), 18,ptRange);
        fOutput->Add(fptTrueHFE[i]);
    }
    
    fptJetIE= new TH1F("fptJetIE", "fptJetIE", 100, 0, 200);
    fOutput->Add(fptJetIE);
    
    fptJetPE= new TH1F("fptJetPE", "fptJetPE", 100, 0, 200);
    fOutput->Add(fptJetPE);
    
    fptJetHFE= new TH1F("fptJetHFE", "fptJetHFE", 100, 0, 200);
    fOutput->Add(fptJetHFE);
    
    fptRecPE= new TH1F("fptRecPE", "fptRecPE", 100, 0, 50);
    fOutput->Add(fptRecPE);
    
    fptTruePE= new TH1F("fptTruePE", "fptTruePE", 100, 0, 50);
    fOutput->Add(fptTruePE);
    
    fptWrongPE= new TH1F("fptWrongPE", "fptWrongPE", 100, 0, 50);
    fOutput->Add(fptWrongPE);
    
    fPhiTrack= new TH2F("fPhiTrack", "fPhiTrack", 100, 0, 50, 100, 0, TMath::TwoPi());
    fOutput->Add(fPhiTrack);
    
    fEtaTrack= new TH2F("fEtaTrack", "fEtaTrack", 100, 0, 50, 100, -1,1);
    fOutput->Add(fEtaTrack);
    
    fPhiElec= new TH2F("fPhiElec", "fPhiElec", 100, 0, 50, 100, 0, TMath::TwoPi());
    fOutput->Add(fPhiElec);
    
    fEtaElec= new TH2F("fEtaElec", "fEtaElec", 100, 0, 50, 100, -1,1);
    fOutput->Add(fEtaElec);
    
    // =========== Switch on Sumw2 for all histos ===========
    for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
        TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
        if (h1){
            h1->Sumw2();
            continue;
        }
        THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
        if(hn)hn->Sumw2();
    }
    
    
    TH1::AddDirectory(oldStatus);
    const Int_t nVar = 24;
    
    fTreeObservableTagging = new TTree(GetOutputSlot(2)->GetContainer()->GetName(), GetOutputSlot(2)->GetContainer()->GetName());
    
    TString *fShapesVarNames = new TString [nVar];
    
    fShapesVarNames[0] = "partonCode";
    fShapesVarNames[1] = "ptJet";
    fShapesVarNames[2] = "ptDJet";
    fShapesVarNames[3] = "mJet";
    // fShapesVarNames[4] = "nbOfConst";
    fShapesVarNames[4] = "angularity";
    fShapesVarNames[5] = "circularity";
    fShapesVarNames[6] = "lesub";
    fShapesVarNames[7] = "coronna";
    
    fShapesVarNames[8] = "ptJetMatch";
    fShapesVarNames[9] = "ptDJetMatch";
    fShapesVarNames[10] = "mJetMatch";
    // fShapesVarNames[12] = "nbOfConstMatch";
    fShapesVarNames[11] = "angularityMatch";
    fShapesVarNames[12] = "circularityMatch";
    fShapesVarNames[13] = "lesubMatch";
    fShapesVarNames[14] = "coronnaMatch";
    fShapesVarNames[15]="weightPythia";
    //fShapesVarNames[14]="ntrksEvt";
    fShapesVarNames[16]="rhoVal";
    fShapesVarNames[17]="rhoMassVal";
    fShapesVarNames[18]="ptUnsub";
    fShapesVarNames[19]="pElec";
    fShapesVarNames[20]="ptElec";
    fShapesVarNames[21]="nInclElec";
    fShapesVarNames[22]="nPhotElec";
    fShapesVarNames[23]="hasElec";
    
    
    for(Int_t ivar=0; ivar < nVar; ivar++){
        //cout<<"looping over variables"<<endl;
        fTreeObservableTagging->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));
    }
    
    PostData(1,fOutput);
    PostData(2,fTreeObservableTagging);
    
    delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::Run()
{
    //  // Run analysis code here, if needed. It will be executed before FillHistograms().
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD){
        printf("ERROR: fAOD not available\n");
        return kFALSE;
    }
    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if(!fVevent){
        printf("ERROR: fVEvent not available\n");
        return kFALSE;
    }
    
    //event selection
    pVtx = fVevent->GetPrimaryVertex();
    spdVtx = fVevent->GetPrimaryVertexSPD();
    
    fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    
    Int_t NpureMC = -1;
    
    if (fMCheader){
        TList *lh=fMCheader->GetCocktailHeaders();
        if(lh){
            AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(0);
            NpureMC = gh->NProduced();
        }
    }
    // pileup rejection with SPD multiple vertices and track multivertexer (using default configuration)
    
    Bool_t isPileupfromSPDmulbins = kFALSE, isPileupFromMV = kFALSE;
    isPileupfromSPDmulbins = fAOD->IsPileupFromSPDInMultBins();
    if (isPileupfromSPDmulbins) return kFALSE;
    
    AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(5);
    utils.SetMaxPlpChi2MV(5.);
    utils.SetMinWDistMV(15);
    utils.SetCheckPlpFromDifferentBCMV(kFALSE);
    isPileupFromMV = utils.IsPileUpMV(fAOD);
    if (isPileupFromMV) return kFALSE;
    
    // Selection of pi0 and eta in MC to compute the weight
    
    Bool_t isPrimary = kFALSE, isFromLMdecay = kTRUE, isFromHFdecay=kTRUE;
    
    Double_t MCweight = 1.;
    Int_t iDecay = 0;
    
    if(fMCarray){
        Int_t nParticles = fMCarray->GetEntries();
        for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
            AliAODMCParticle* particle = (AliAODMCParticle*) fMCarray->At(iParticle);
            int fPDG = particle->GetPdgCode();
            double pTMC = particle->Pt();
            
            Bool_t iEnhance = kFALSE;
            if(iParticle>=NpureMC)iEnhance = kTRUE;
            
            Double_t etaMC = particle->Eta();
            //if (TMath::Abs(etaMC)>1.2)continue;
            
            isPrimary = IsPrimary(particle);
            isFromLMdecay = IsFromLMdecay(particle);
            isFromHFdecay = IsFromHFdecay(particle);
            
            MCweight = 1.;
            iDecay = 0;
            
            GetWeightAndDecay(particle,iDecay,MCweight);
            
            
            if (TMath::Abs(etaMC)<0.8 && TMath::Abs(fPDG)==11){
                
                if (isFromHFdecay) fGenHfePt->Fill(pTMC);
                if (iDecay>0 && iDecay<7) fGenPePt->Fill(pTMC);
            }
            
            
            Double_t yMC = particle->Y();
            if (TMath::Abs(yMC)>1.0)continue;
            
            if (isPrimary){
                if(fPDG==111) fPi0Pt->Fill(iEnhance,pTMC); //pi0
                if(fPDG==221) fEtaPt->Fill(iEnhance,pTMC); //eta
            }
            if (!isFromHFdecay && !isFromLMdecay){
                if(fPDG==111) fPi0Pt->Fill(iEnhance+2,pTMC); //pi0
                if(fPDG==221) fEtaPt->Fill(iEnhance+2,pTMC); //eta
            }
        }
    }//MC
    
    
    return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::FillHistograms()
{
    // Fill histograms.
    //cout<<"base container"<<endl;
    AliEmcalJet* jet1 = NULL;
    AliJetContainer *jetCont = GetJetContainer(0);
    
    Float_t kWeight=1;
    if (fCentSelectOn)
        if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
    
    AliAODTrack *triggerHadron = 0x0;
    
    if (fJetSelection == kRecoil) {
        //Printf("Recoil jets!!!, fminpTTrig = %f, fmaxpTTrig = %f", fminpTTrig, fmaxpTTrig);
        Int_t triggerHadronLabel = SelectTrigger(fminpTTrig, fmaxpTTrig);
        
        
        if (triggerHadronLabel==-99999) {
            //Printf ("Trigger Hadron not found, return");
            return 0;}
        
        AliTrackContainer *PartCont =NULL;
        AliParticleContainer *PartContMC=NULL;
        
        if (fJetShapeSub==kConstSub){
            if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) PartContMC = GetParticleContainer(1);
            else PartCont = GetTrackContainer(1);
        }
        else{
            if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) PartContMC = GetParticleContainer(0);
            else PartCont = GetTrackContainer(0);
        }
        TClonesArray *TrackArray = NULL;
        TClonesArray *TrackArrayMC = NULL;
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) TrackArrayMC = PartContMC->GetArray();
        else TrackArray = PartCont->GetArray();
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) triggerHadron = static_cast<AliAODTrack*>(TrackArrayMC->At(triggerHadronLabel));
        else triggerHadron = static_cast<AliAODTrack*>(TrackArray->At(triggerHadronLabel));
        
        
        
        if (!triggerHadron) {
            //Printf("No Trigger hadron with the found label!!");
            return 0;
        }
        
        if(fSemigoodCorrect){
            Double_t disthole=RelativePhi(triggerHadron->Phi(),fHolePos);
            if(TMath::Abs(disthole)+fHoleWidth>TMath::Pi()-fangWindowRecoil){
                return 0;}
        }
        
        fhPt->Fill(triggerHadron->Pt());
        
    }
    
    // list of kink mothers
    Int_t nVertices = 1;
    nVertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[nVertices];
    Int_t nMotherKink = 0;
    for(Int_t ivertex=0; ivertex < nVertices; ivertex++) {
        AliAODVertex *vertex = fAOD->GetVertex(ivertex);
        if(!vertex) continue;
        if(vertex->GetType()==AliAODVertex::kKink) {
            AliAODTrack *mother = (AliAODTrack *) vertex->GetParent();
            if(!mother) continue;
            Int_t idmother = mother->GetID();
            listofmotherkink[nMotherKink] = idmother;
            nMotherKink++;
        }
    }
    
    //AliParticleContainer *partContAn = GetParticleContainer(0);
    //TClonesArray *trackArrayAn = partContAn->GetArray();
    //Int_t ntracksEvt = trackArrayAn->GetEntriesFast();
    
    Float_t rhoVal=0, rhoMassVal = 0.;
    if(jetCont) {
        
        jetCont->ResetCurrentID();
        if ((fJetShapeSub==kConstSub) || (fJetShapeSub==kDerivSub)){
            //rho
            AliRhoParameter* rhoParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoSparseR020"));
            if (!rhoParam) {
                Printf("%s: Could not retrieve rho %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoName().Data());
            } else rhoVal = rhoParam->GetVal();
            //rhom
            AliRhoParameter* rhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoMassSparseR020"));
            if (!rhomParam) {
                Printf("%s: Could not retrieve rho_m %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoMassName().Data());
            } else rhoMassVal = rhomParam->GetVal();
        }
        
        while((jet1 = jetCont->GetNextAcceptJet())) {
            if (!jet1) continue;
            AliEmcalJet* jet2 = 0x0;
            AliEmcalJet* jet3 = 0x0;
            fPtJet->Fill(jet1->Pt());
            fPhiJet->Fill(jet1->Pt(),jet1->Phi());
            fEtaJet->Fill(jet1->Pt(),jet1->Eta());
            AliEmcalJet *jetUS = NULL;
            Int_t ifound=0, jfound=0;
            Int_t ilab=-1, jlab=-1;
            
            if(fSemigoodCorrect && (fJetSelection != kRecoil)){
                Double_t disthole=RelativePhi(jet1->Phi(),fHolePos);
                if(TMath::Abs(disthole)<fHoleWidth){
                    continue;}
            }
            
            Float_t dphiRecoil = 0.;
            if (fJetSelection == kRecoil){
                dphiRecoil = RelativePhi(triggerHadron->Phi(), jet1->Phi());
                if (TMath::Abs(dphiRecoil) < (TMath::Pi() - fangWindowRecoil)) {
                    // Printf("Recoil jets back to back not found! continuing");
                    continue;
                }
                
                fhpTjetpT->Fill(triggerHadron->Pt(), jet1->Pt());
                //Printf(" ************ FILLING HISTOS****** shapeSub = %d, triggerHadron = %f, jet1 = %f", fJetShapeSub, triggerHadron->Pt(), jet1->Pt());
                fhPhi->Fill(RelativePhi(triggerHadron->Phi(), jet1->Phi()));
                
            }
            
            
            fShapesVar[0] = 0.;
            if(fJetShapeType == kDetEmbPartPythia){
                //AliJetContainer *jetContTrue = GetJetContainer(1);
                AliJetContainer *jetContUS = GetJetContainer(2);
                
                if(fJetShapeSub==kConstSub){
                    for(Int_t i = 0; i<jetContUS->GetNJets(); i++) {
                        jetUS = jetContUS->GetJet(i);
                        if(jetUS->GetLabel()==jet1->GetLabel()) {
                            ifound++;
                            if(ifound==1) ilab = i;
                        }
                    }
                    if(ilab==-1) continue;
                    jetUS=jetContUS->GetJet(ilab);
                    jet2=jetUS->ClosestJet();
                }
                
                if(!(fJetShapeSub==kConstSub)) jet2 = jet1->ClosestJet();
                if (!jet2) {
                    Printf("jet2 does not exist, returning");
                    continue;
                }
                
                //AliJetContainer *jetContPart=GetJetContainer(3);
                jet3=jet2->ClosestJet();
                
                if(!jet3){
                    Printf("jet3 does not exist, returning");
                    continue;
                }
                //cout<<"jet 3 exists"<<jet3->Pt()<<endl;
                
                
                fh2ResponseUW->Fill(jet1->Pt(),jet2->Pt());
                
                Double_t fraction=0;
                if(!(fJetShapeSub==kConstSub))  fraction = jetCont->GetFractionSharedPt(jet1);
                if(fJetShapeSub==kConstSub) fraction = jetContUS->GetFractionSharedPt(jetUS);
                //if (fraction > 0.1) cout<<"***** hey a jet matched with fraction"<<fraction<<"  "<<jet1->Pt()<<" "<<jet2->Pt()<<" "<<fCent<<endl;
                
                if(fraction<fMinFractionShared) continue;
                //InputEvent()->Print();
                
            }
            
            
            
            if (fJetShapeType == kPythiaDef){
                
                //AliJetContainer *jetContTrue = GetJetContainer(1);
                AliJetContainer *jetContUS = GetJetContainer(2);
                AliJetContainer *jetContPart = GetJetContainer(3);
                
                if(fJetShapeSub==kConstSub){
                    
                    for(Int_t i = 0; i<jetContUS->GetNJets(); i++) {
                        jetUS = jetContUS->GetJet(i);
                        if(jetUS->GetLabel()==jet1->GetLabel()) {
                            ifound++;
                            if(ifound==1) ilab = i;
                        }
                    }
                    if(ilab==-1) continue;
                    jetUS=jetContUS->GetJet(ilab);
                    jet2=jetUS->ClosestJet();
                    
                    if (!jet2) {
                        Printf("jet2 does not exist, returning");
                        continue;
                    }
                    
                    for(Int_t j=0; j<jetContPart->GetNJets(); j++) {
                        
                        jet3 = jetContPart->GetJet(j);
                        if(!jet3) continue;
                        if(jet3->GetLabel()==jet2->GetLabel()) {
                            jfound++;
                            if(jfound==1) jlab = j;
                        }
                    }
                    if(jlab==-1) continue;
                    jet3=jetContPart->GetJet(jlab);
                    if(!jet3){
                        Printf("jet3 does not exist, returning");
                        continue;
                    }
                }
                if(!(fJetShapeSub==kConstSub)) jet3 = jet1->ClosestJet();
                if (!jet3) {
                    Printf("jet3 does not exist, returning");
                    continue;
                }
                
                
                fh2ResponseUW->Fill(jet1->Pt(),jet3->Pt());
                
                Double_t eff = (jet3->Pt()-jet1->Pt())/jet1->Pt();
                fJetEfficiency->Fill(eff,jet1->Pt());
                
                
            }
            
            
            if (fJetShapeType == kGenOnTheFly){
                const AliEmcalPythiaInfo *partonsInfo = 0x0;
                partonsInfo = GetPythiaInfo();
                Double_t jp1=RelativePhi(jet1->Phi(),partonsInfo->GetPartonPhi6());
                Double_t detap1=(jet1->Eta())-(partonsInfo->GetPartonEta6());
                kWeight=partonsInfo->GetPythiaEventWeight();
                fh2ResponseW->Fill(jet1->Pt(),jet1->Pt(),kWeight);
                
                Float_t dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
                fEtaJetCorr6->Fill(jet1->Eta(), partonsInfo->GetPartonEta6());
                fPhiJetCorr6->Fill(jet1->Phi(), partonsInfo->GetPartonPhi6());
                if(dRp1 < fRMatching) {
                    fShapesVar[0] = partonsInfo->GetPartonFlag6();
                    fPtJetCorr ->Fill(partonsInfo->GetPartonPt6(), jet1->Pt());
                }
                else {
                    jp1=RelativePhi(jet1->Phi(),partonsInfo->GetPartonPhi7());
                    detap1=(jet1->Eta())-(partonsInfo->GetPartonEta7());
                    dRp1 = TMath::Sqrt(jp1 * jp1 + detap1 * detap1);
                    fEtaJetCorr7->Fill(jet1->Eta(), partonsInfo->GetPartonEta7());
                    fPhiJetCorr7->Fill(jet1->Phi(), partonsInfo->GetPartonPhi7());
                    if(dRp1 < fRMatching) {
                        fShapesVar[0] = partonsInfo->GetPartonFlag7();
                        fPtJetCorr->Fill(partonsInfo->GetPartonPt7(), jet1->Pt());
                    }
                    else fShapesVar[0]=0;
                }
            }
            
            
            
            
            Double_t ptSubtracted = 0;
            if (fJetShapeSub==kConstSub) ptSubtracted= jet1->Pt();
            
            else if (fJetShapeSub==kDerivSub)  {
                ptSubtracted=jet1->Pt()-GetRhoVal(0)*jet1->Area();
            }
            
            else if (fJetShapeSub==kNoSub) {
                if ((fJetShapeType==kData) || (fJetShapeType==kDetEmbPartPythia)) ptSubtracted=jet1->Pt()-GetRhoVal(0)*jet1->Area();
                else if ((fJetShapeType==kPythiaDef) || (fJetShapeType==kMCTrue) || (fJetShapeType==kGenOnTheFly)) ptSubtracted= jet1->Pt();
            }
            
            //Printf("ptSubtracted=%f,fPtThreshold =%f ", ptSubtracted, fPtThreshold);
            if (ptSubtracted < fPtThreshold) continue;
            
            if (fOneConstSelectOn == kTRUE) fNbOfConstvspT->Fill(GetJetNumberOfConstituents(jet1,0), ptSubtracted);
            if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1)) continue;
            
            
            fShapesVar[1] = ptSubtracted;
            fShapesVar[2] = GetJetpTD(jet1,0);
            fShapesVar[3] = GetJetMass(jet1,0);
            fShapesVar[4] = GetJetAngularity(jet1,0);
            fShapesVar[5] = GetJetCircularity(jet1,0);
            fShapesVar[6] = GetJetLeSub(jet1,0);
            fShapesVar[6] = GetJetCoronna(jet1,0);
            
            
            
            Float_t ptMatch=0., ptDMatch=0., massMatch=0., angulMatch=0.,circMatch=0., lesubMatch=0., coronnaMatch=0;
            //Float constMatch=0., sigma2Match=0.;
            Int_t kMatched = 0;
            
            if (fJetShapeType==kPythiaDef) {
                kMatched =1;
                if(fJetShapeSub==kConstSub) kMatched = 3;
                
                ptMatch=jet3->Pt();
                ptDMatch=GetJetpTD(jet3, kMatched);
                massMatch=GetJetMass(jet3,kMatched);
                //constMatch=1.*GetJetNumberOfConstituents(jet2,kMatched);
                angulMatch=GetJetAngularity(jet3, kMatched);
                circMatch=GetJetCircularity(jet3, kMatched);
                lesubMatch=GetJetLeSub(jet3, kMatched);
                coronnaMatch=GetJetCoronna(jet3,kMatched);
                //sigma2Match = GetSigma2(jet2, kMatched);
            }
            
            if (fJetShapeType==kDetEmbPartPythia) {
                if(fJetShapeSub==kConstSub) kMatched = 3;
                if(fJetShapeSub==kDerivSub) kMatched = 2;
                ptMatch=jet3->Pt();
                ptDMatch=GetJetpTD(jet3, kMatched);
                massMatch=GetJetMass(jet3,kMatched);
                // constMatch=1.*GetJetNumberOfConstituents(jet3,kMatched);
                angulMatch=GetJetAngularity(jet3, kMatched);
                circMatch=GetJetCircularity(jet3, kMatched);
                lesubMatch=GetJetLeSub(jet3, kMatched);
                coronnaMatch = GetJetCoronna(jet3, kMatched);
            }
            
            
            
            if (fJetShapeType == kMCTrue || fJetShapeType == kData || fJetShapeType == kGenOnTheFly) {
                kMatched = 0;
                ptMatch=0.;
                ptDMatch=0.;
                massMatch=0.;
                //constMatch=0.;
                angulMatch=0.;
                circMatch=0.;
                lesubMatch=0.;
                coronnaMatch =0.;
                
            }
            
            
            
            fShapesVar[8] = ptMatch;
            fShapesVar[9] = ptDMatch;
            fShapesVar[10] = massMatch;
            fShapesVar[11] = angulMatch;
            fShapesVar[12] = circMatch;
            fShapesVar[13] = lesubMatch;
            fShapesVar[14] = coronnaMatch;
            fShapesVar[15] = kWeight;
            //fShapesVar[16] = ntracksEvt;
            fShapesVar[16] = rhoVal;
            fShapesVar[17] = rhoMassVal;
            fShapesVar[18] = jet1->Pt();
            
            Int_t nInclusiveElectrons = 0, nPhotonicElectrons = 0, nTrueElectronsMC, nTrueHFElecMC;
            Double_t pElec = 0., ptElec = 0.;
            Bool_t hasElectrons = kFALSE;
            
            
            if(fJetShapeType != kMCTrue){
                GetNumberOfElectrons(jet1, 0,nMotherKink,listofmotherkink,nInclusiveElectrons,nPhotonicElectrons,pElec,ptElec,hasElectrons);
                GetNumberOfTrueElectrons(jet1, 0,nMotherKink,listofmotherkink,nTrueElectronsMC,nTrueHFElecMC);
            }
            
            
            // generated HFE jets
            
            AliVParticle *vp1 = 0x0;
            Int_t pdgCode = 0;
            Int_t elecCounter = 0;
            
            if(fMCarray){
                for(UInt_t i = 0; i < jet1->GetNumberOfTracks(); i++) {
                    vp1 = static_cast<AliVParticle*>(jet1->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
                    
                    if (!vp1){
                        Printf("AliVParticle associated to constituent not found");
                        continue;
                    }
                    
                    pdgCode = vp1->PdgCode();
                    if (TMath::Abs(pdgCode)==11) elecCounter++;
                }
                if (elecCounter==1) fPtGenJet->Fill(jet1->Pt());
            }
            
            fShapesVar[19] = pElec;
            fShapesVar[20] = ptElec;
            
            fShapesVar[21] = nInclusiveElectrons;
            fShapesVar[22] = nPhotonicElectrons;
            fShapesVar[23] = hasElectrons;
            
            fTreeObservableTagging->Fill();
            
            
        }//jet loop
    }
    return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetMass(AliEmcalJet *jet,Int_t jetContNb=0){
    //calc subtracted jet mass
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtracted();
        else return jet->GetShapeProperties()->GetSecondOrderSubtracted();
        else
            return jet->M();
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHfeTagging::GetNumberOfElectrons(AliEmcalJet *jet, Int_t jetContNb , Int_t nMother, Double_t listMother[] ,  Int_t &nIncElec,  Int_t &nPhotElec, Double_t &pElec, Double_t &ptElec, Bool_t &hasElec){
    // count the number of inclusive electrons per jet and per event
    
    AliVParticle *vp1 = 0x0;
    Int_t nIE=0, nPairs=0, nPE=0, sub = -1;
    Float_t ratioElec = -1.;
    Double_t pte=0., pe=0.;
    Bool_t hasElecCand = kFALSE;
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (jet->GetNumberOfTracks()){
        
        fnPartPerJet->Fill(jet->GetNumberOfTracks());
        fpidResponse = fInputHandler->GetPIDResponse();
        if(!fpidResponse){
            printf("ERROR: pid response not available\n");
        }
        
        for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
            vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
            
            if (!vp1){
                Printf("AliVParticle associated to constituent not found");
                continue;
            }
            
            AliVTrack *vtrack = dynamic_cast<AliVTrack*>(vp1);
            if (!vtrack) {
                printf("ERROR: Could not receive track%d\n", i);
                continue;
            }
            
            AliAODTrack *track = dynamic_cast<AliAODTrack*>(vtrack);
            
            if (!track) continue;
            
            // track cuts
            Bool_t passTrackCut = kFALSE;
            passTrackCut = InclElecTrackCuts(pVtx,track,nMother,listMother);
            if (!passTrackCut) continue;
            
            nPairs=0;
            
            Double_t p=-9., pt=-9., fTPCnSigma=-99., fTOFnSigma=-99., phi = -9., eta = -99.;
            p = track->P();
            pt = track->Pt();
            phi = track->Phi();
            eta = track->Eta();
            
            fPhiTrack->Fill(pt,phi);
            fEtaTrack->Fill(pt,eta);
            
            fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            fTOFnSigma = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
            
            fnTPCnTOFnocut->Fill(fTPCnSigma,fTOFnSigma);
            fnTPCnocut->Fill(p,fTPCnSigma);
            fnTOFnocut->Fill(p,fTOFnSigma);
            if (TMath::Abs(fTOFnSigma)<fSigmaTOFcut) fnTPCcut->Fill(p,fTPCnSigma);
            
            if(TMath::Abs(fTPCnSigma)<3.5){
                hasElecCand = kTRUE;
            }
            
            if ((fTPCnSigma>fSigmaTPCcut)  && (fTPCnSigma<3) && TMath::Abs(fTOFnSigma)<fSigmaTOFcut){
                fPtP->Fill(pt,p);
                fPhiElec->Fill(pt,phi);
                fEtaElec->Fill(pt,eta);
                
                nIE++;
                pte=pt;
                pe=p;
                nPairs = GetNumberOfPairs(jet,track,pVtx,nMother,listMother);
                if (nPairs>0) nPE++;
            }
            
        }//tagged track
        sub = nIE - nPE;
        if (jet->GetNumberOfTracks()>0) ratioElec = nIE/jet->GetNumberOfTracks();
        fnInclElecPerJet->Fill(nIE);
        fnPhotElecPerJet->Fill(nPE);
        fnIncSubPhotElecPerJet->Fill(sub);
        fnElecOverPartPerJet->Fill(ratioElec);
    }
    
    nIncElec = nIE;
    nPhotElec = nPE;
    pElec = pe;
    ptElec = pte;
    hasElec = hasElecCand;
}
//________________________________________________________________________
void AliAnalysisTaskEmcalHfeTagging::GetNumberOfTrueElectrons(AliEmcalJet *jet, Int_t jetContNb ,  Int_t nMother, Double_t listMother[] ,  Int_t &nTrueElec,  Int_t &nTrueHFElec){
    // count the number of inclusive and HF electrons per jet and per event (true MC)
    
    AliVParticle *vp1 = 0x0;
    Int_t nIE=0, nHFE=0, nPE=0, PartId=0, nPairs=0, iDecay = 0;
    Double_t p=-9., pt=-9., fTPCnSigma=-99., fTOFnSigma=-99., MCweight = 1.;
    Double_t ptRange[19] = {0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,5};
    Double_t ptJetRange[6] = {5,15,25,35,45,60};
    
    Bool_t isFromHFdecay=kFALSE;
    Bool_t isFromLMdecay=kFALSE;
    Bool_t passTrackCut = kFALSE;
    
    Double_t jetPt = jet->Pt();
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (jet->GetNumberOfTracks()){
        for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
            vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
            
            if (!vp1){
                Printf("AliVParticle associated to constituent not found");
                continue;
            }
            
            AliVTrack *vtrack = dynamic_cast<AliVTrack*>(vp1);
            if (!vtrack) {
                printf("ERROR: Could not receive track%d\n", i);
                continue;
            }
            
            AliAODTrack *track = dynamic_cast<AliAODTrack*>(vtrack);
            
            passTrackCut = kFALSE;
            isFromHFdecay=kFALSE;
            isFromLMdecay=kFALSE;
            
            p = track->P();
            pt = track->Pt();
            
            fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
            fTOFnSigma = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
            
            nPairs=0;
            MCweight = 1.;
            iDecay = 0;
            
            if(fMCarray){
                Int_t label = track->GetLabel();
                if(label!=0){
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(label);
                    if(fMCparticle){
                        Int_t partPDG = fMCparticle->GetPdgCode();
                        
                        GetWeightAndDecay(fMCparticle,iDecay,MCweight);
                        isFromHFdecay = IsFromHFdecay(fMCparticle);
                        isFromLMdecay = IsFromLMdecay(fMCparticle);
                        
                        if (TMath::Abs(partPDG)==11 && isFromHFdecay) fptTrueHFE[0]->Fill(pt);
                        
                        // track cuts
                        passTrackCut = InclElecTrackCuts(pVtx,track,nMother,listMother);
                        if (!passTrackCut) continue;
                        
                        if (TMath::Abs(partPDG)==11 && isFromHFdecay) fptTrueHFE[1]->Fill(pt);
                        
                        //check sigma_TPC in MC for different particles
                        if (TMath::Abs(fTOFnSigma)<fSigmaTOFcut){
                            if (TMath::Abs(partPDG)==11) PartId = 1; // electrons
                            if (TMath::Abs(partPDG)==13) PartId = 2; // muons
                            if (TMath::Abs(partPDG)==321) PartId = 3; // kaons
                            if (TMath::Abs(partPDG)==2212) PartId = 4; // protons
                            if (TMath::Abs(partPDG)==211) PartId = 5; // pions
                            
                            for (Int_t l=0;l<5;l++){// pt jet range
                                for(Int_t k=0;k<18;k++){// pt electron range
                                    if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1] && p>=ptRange[k] && p<ptRange[k+1]) fnTPCTrueParticles[l][k]->Fill(PartId,fTPCnSigma);
                                }
                            }
                        }
                        
                        if (TMath::Abs(partPDG)!=11) continue;
                        
                        if (isFromHFdecay && (fTPCnSigma>fSigmaTPCcut) && (fTPCnSigma<3)) fptTrueHFE[2]->Fill(pt);
                        if (isFromHFdecay && (fTPCnSigma>fSigmaTPCcut) && (fTPCnSigma<3) && (TMath::Abs(fTOFnSigma)<fSigmaTOFcut)) fptTrueHFE[3]->Fill(pt);
                        
                        
                        nIE++;
                        if (isFromHFdecay) nHFE++;
                        if (isFromLMdecay) nPE++;
                        
                        if ((fTPCnSigma>fSigmaTPCcut) && (fTPCnSigma<3) && (TMath::Abs(fTOFnSigma)<fSigmaTOFcut)){
                            
                            nPairs = GetNumberOfPairs(jet,track,pVtx,nMother,listMother);
                            if (nPairs>0) fptRecPE->Fill(pt,MCweight);
                            if (nPairs>0 && (iDecay<1 || iDecay>6)) fptWrongPE->Fill(pt,MCweight);
                            
                            if (iDecay>0 && iDecay<7){
                                fptTruePE->Fill(pt,MCweight);
                                for (Int_t l=0;l<5;l++){// pt jet range
                                    
                                    if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1]){
                                        fTotElecPt[l]->Fill(pt,MCweight);
                                        if (nPairs>0) fUlsLsElecPt[l]->Fill(pt,MCweight);
                                    }
                                }//jet pt
                            } // decay channels
                        }// PID cuts
                    }
                }
            }
        }//track loop
        
        if (nIE==1) fptJetIE->Fill(jet->Pt());
        if (nIE==1 && nPE==1) fptJetPE->Fill(jet->Pt());
        if (nIE==1 && nHFE==1) fptJetHFE->Fill(jet->Pt());
        
        fnTrueElecPerJet->Fill(nIE);
        fnTrueHFElecPerJet->Fill(nHFE);
        fnTruePElecPerJet->Fill(nPE);
    }
    
    nTrueElec = nIE;
    nTrueHFElec = nHFE;
}

//_________________________________________
Int_t AliAnalysisTaskEmcalHfeTagging::GetNumberOfPairs(AliEmcalJet *jet, AliAODTrack *track,const AliVVertex *pVtx, Int_t nMother, Double_t listMother[])
{
    //Get number of ULS and LS pairs per event
    
    Int_t ntracks = 0;
    ntracks = fVevent->GetNumberOfTracks();
    
    Int_t nULSpairs = 0;
    Int_t nLSpairs = 0;
    Int_t sub = -1;
    
    Double_t ptJetRange[6] = {5,15,25,35,45,60};
    Double_t jetPt = jet->Pt();
    
    for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
        
        AliVParticle* Vassotrack = fVevent->GetTrack(jTracks);
        
        if (!Vassotrack) {
            printf("ERROR: Could not receive associated track %d\n", jTracks);
            continue;
        }
        AliAODTrack *trackAsso = dynamic_cast<AliAODTrack*>(Vassotrack);
        if(!trackAsso) continue;
        
        if((track->GetID())==(trackAsso->GetID())) continue;
        
        // track cuts
        Bool_t passAssoTrackCut = kFALSE;
        passAssoTrackCut = PhotElecTrackCuts(pVtx,trackAsso,nMother,listMother);
        if (!passAssoTrackCut) continue;
        
        // tagged particle
        Double_t p=-9.,pt=-9.,eta =-9.,phi=-9.;
        Int_t charge = 0;
        
        pt = track->Pt();
        p = track->P();
        phi = track->Phi();
        eta = track->Eta();
        charge = track->Charge();
        
        
        // associated particle variables
        Double_t pAsso=-9.,ptAsso=-9.,etaAsso =-9.,fITSnSigmaAsso=-9.,fTPCnSigmaAsso=-9.,phiAsso=-9.;
        Int_t chargeAsso = 0;
        
        ptAsso = trackAsso->Pt();
        pAsso = trackAsso->P();
        phiAsso = trackAsso->Phi();
        etaAsso = trackAsso->Eta();
        chargeAsso = trackAsso->Charge();
        
        
        // looser PID cuts
        fITSnSigmaAsso = fpidResponse->NumberOfSigmasITS(trackAsso, AliPID::kElectron);
        fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(trackAsso, AliPID::kElectron);
        
        if(TMath::Abs(fTPCnSigmaAsso)>3) continue;
        
        // invariant mass
        Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
        Double_t openingAngle = -999., mass=999., width = -999;
        
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        AliKFParticle::SetField(fVevent->GetMagneticField());
        AliKFParticle ge1(*track, fPDGe1);
        AliKFParticle ge2(*trackAsso, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        openingAngle = ge1.GetAngle(ge2);
        //if(openingAngle > fOpeningAngleCut) continue;
        
        Int_t MassCorrect=-9;
        MassCorrect = recg.GetMass(mass,width);
        
        for (Int_t l=0;l<5;l++){// pt jet range
            if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1] && fFlagULS) fInvmassULS[l]->Fill(pt,mass);
            if (jetPt>=ptJetRange[l] && jetPt<ptJetRange[l+1] && fFlagLS)  fInvmassLS[l]->Fill(pt,mass);
        }
        
        if(mass<fIMcut && fFlagULS)
            nULSpairs++;
        
        if(mass<fIMcut && fFlagLS)
            nLSpairs++;
        
    }//track loop
    
    sub = nULSpairs-nLSpairs;
    fnULSmLSpairsPerElectron->Fill(track->Pt(),sub);
    
    return sub;
}


//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::IsFromHFdecay(AliAODMCParticle *particle)
{
    // Check if the mother comes from heavy-flavour decays
    Bool_t isHFdecay = kFALSE;
    //Int_t partPDG = particle->GetPdgCode();
    
    Int_t idMother = particle->GetMother();
    if (idMother>0){
        AliAODMCParticle* mother = (AliAODMCParticle*) fMCarray->At(idMother);
        Int_t motherPDG = TMath::Abs(mother->GetPdgCode());
        
        // c decays
        if((motherPDG % 1000) / 100 == 4) isHFdecay = kTRUE;
        if(motherPDG / 1000 == 4) isHFdecay = kTRUE;
        
        // b decays
        if((motherPDG % 1000) / 100 == 5) isHFdecay = kTRUE;
        if(motherPDG / 1000 == 5) isHFdecay = kTRUE;
        
        // all particles related to  HF
        if(motherPDG==4 || motherPDG==5 || motherPDG == 211 || motherPDG ==13 || motherPDG ==2112 || motherPDG ==130 || motherPDG ==3122 ||
           motherPDG ==310 || motherPDG ==3222 || motherPDG ==2212 || motherPDG ==3112 || motherPDG ==321 ||
           motherPDG ==11 || motherPDG ==3212 || motherPDG ==311 || motherPDG ==20213 || motherPDG ==3312 ||
           motherPDG ==3334 || motherPDG ==3324 || motherPDG ==3322 || motherPDG ==1000010020 || motherPDG ==15
           || motherPDG ==10323 || motherPDG ==2114 || motherPDG ==1000010030 || motherPDG ==2214 || motherPDG ==2224
           || motherPDG ==1114 || motherPDG == 2214 || motherPDG == 3114 || motherPDG ==3224 || motherPDG ==3124
           || motherPDG ==3314 || motherPDG ==10323 || motherPDG == 3214) isHFdecay = kTRUE;
        
    }
    
    return isHFdecay;
}
//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::IsFromLMdecay(AliAODMCParticle *particle)
{
    // Check if  mother comes from light-meson decays
    Bool_t isLMdecay = kFALSE;
    //Int_t partPDG = particle->GetPdgCode();
    
    Int_t idMother = particle->GetMother();
    if (idMother>0){
        AliAODMCParticle* mother = (AliAODMCParticle*) fMCarray->At(idMother);
        Int_t motherPDG = TMath::Abs(mother->GetPdgCode());
        
        if(motherPDG == 111 || motherPDG == 221 || motherPDG==223 || motherPDG==333 || motherPDG==331 ||
           motherPDG==113 || motherPDG==213 || motherPDG==313 || motherPDG==323) isLMdecay = kTRUE;
    }
    
    return isLMdecay;
}
//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::IsPrimary(AliAODMCParticle *particle)
{
    // Check if  particle is primary
    Bool_t isprimary = kFALSE;
    
    //Int_t idMother = particle->GetMother();
    //if (idMother==-1) isprimary = kTRUE;
    if (particle->IsPrimary()) isprimary = kTRUE;
    
    return isprimary;
}
//_________________________________________
Double_t AliAnalysisTaskEmcalHfeTagging::GetPi0weight(Double_t mcPi0pT) const
{
    //Get Pi0 weight
    double weight = 1.0;
    
    double parPi0_enh[5] = {0.530499,0.732775,0.000997414,3.46894,4.84342};
    if (fMCweight==1) weight = parPi0_enh[0] / TMath::Power((exp(parPi0_enh[1]*mcPi0pT - parPi0_enh[2]*mcPi0pT*mcPi0pT) + (mcPi0pT/parPi0_enh[3])) , parPi0_enh[4]);
    
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskEmcalHfeTagging::GetEtaweight(Double_t mcEtapT) const
{
    //Get Pi0 weight
    double weight = 1.0;
    
    double parEta_enh[5] = {0.78512,-0.606822,0.0326254,3.13959,4.83715};
    if (fMCweight==1) weight = parEta_enh[0] / TMath::Power((exp(parEta_enh[1]*mcEtapT - parEta_enh[2]*mcEtapT*mcEtapT) + (mcEtapT/parEta_enh[3])), parEta_enh[4]);
    
    return weight;
}
//________________________________________________________________________
void AliAnalysisTaskEmcalHfeTagging::GetWeightAndDecay(AliAODMCParticle *particle, Int_t &decay, Double_t &weight)
{
    //Get decay channel and weight for pi0 and eta (MC with enhanced signal)
    Double_t w = 1.;
    Int_t d = 0;
    Int_t partPDG = particle->GetPdgCode();
    
    if (TMath::Abs(partPDG)==11){
        Int_t idMother = particle->GetMother();
        
        if (idMother>0){
            AliAODMCParticle* mother = (AliAODMCParticle*) fMCarray->At(idMother);
            Int_t motherPDG = mother->GetPdgCode();
            Double_t motherPt = mother->Pt();
            
            Bool_t isMotherPrimary = IsPrimary(mother);
            Bool_t isMotherFromHF = IsFromHFdecay(mother);
            Bool_t isMotherFromLM = IsFromLMdecay(mother);
            
            if (motherPDG==111 && (isMotherPrimary || (!isMotherFromHF && !isMotherFromLM))){ // pi0 -> e
                d = 1;
                w = GetPi0weight(motherPt);
            }
            
            if (motherPDG==221  && (isMotherPrimary || (!isMotherFromHF && !isMotherFromLM))){ // eta -> e
                d = 2;
                w = GetEtaweight(motherPt);
            }
            
            //Int_t idSecondMother = particle->GetSecondMother();
            Int_t idSecondMother = mother->GetMother();
            
            if (idSecondMother>0){
                AliAODMCParticle* secondMother = (AliAODMCParticle*) fMCarray->At(idSecondMother);
                Int_t secondMotherPDG = secondMother->GetPdgCode();
                Double_t secondMotherPt = secondMother->Pt();
                
                Bool_t isSecondMotherPrimary = IsPrimary(secondMother);
                Bool_t isSecondMotherFromHF = IsFromHFdecay(secondMother);
                Bool_t isSecondMotherFromLM = IsFromLMdecay(secondMother);
                
                if (motherPDG==22 && secondMotherPDG==111 && (isSecondMotherPrimary || (!isSecondMotherFromHF && !isSecondMotherFromLM))){ //pi0 -> g -> e
                    d = 3;
                    w = GetPi0weight(secondMotherPt);
                }
                
                if (motherPDG==22 && secondMotherPDG==221  && (isSecondMotherPrimary || (!isSecondMotherFromHF && !isSecondMotherFromLM))){ //eta -> g -> e
                    d = 4;
                    w = GetEtaweight(secondMotherPt);
                }
                
                if (motherPDG==111 && secondMotherPDG==221  && (isSecondMotherPrimary || (!isSecondMotherFromHF && !isSecondMotherFromLM))){ //eta -> pi0 -> e
                    d = 5;
                    w = GetEtaweight(secondMotherPt);
                }
                
                Int_t idThirdMother = secondMother->GetMother();
                if (idThirdMother>0){
                    AliAODMCParticle* thirdMother = (AliAODMCParticle*) fMCarray->At(idThirdMother);
                    Int_t thirdMotherPDG = thirdMother->GetPdgCode();
                    Double_t thirdMotherPt = thirdMother->Pt();
                    
                    Bool_t isThirdMotherPrimary = IsPrimary(thirdMother);
                    Bool_t isThirdMotherFromHF = IsFromHFdecay(thirdMother);
                    Bool_t isThirdMotherFromLM = IsFromLMdecay(thirdMother);
                    
                    if (motherPDG==22 && secondMotherPDG==111 && thirdMotherPDG==221 && (isThirdMotherPrimary || (!isThirdMotherFromHF && !isThirdMotherFromLM))){//eta->pi0->g-> e
                        d = 6;
                        w = GetEtaweight(thirdMotherPt);
                    }
                }//third mother
            }//second mother
        }//mother
    }// if electron
    decay = d;
    weight = w;
}
//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::Angularity(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        Double_t dphi = RelativePhi(vp1->Phi(),jet->Phi());
        Double_t dr2 = (vp1->Eta()-jet->Eta())*(vp1->Eta()-jet->Eta()) + dphi*dphi;
        Double_t dr = TMath::Sqrt(dr2);
        num=num+vp1->Pt()*dr;
        den=den+vp1->Pt();
    }
    return num/den;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    if((fJetShapeSub==kDerivSub) && (jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity();
        else
            return Angularity(jet, jetContNb);
    
}

//____________________________________________________________________________

Float_t AliAnalysisTaskEmcalHfeTagging::Coronna(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    AliTrackContainer *PartCont = NULL;
    AliParticleContainer *PartContMC = NULL;
    
    
    if (fJetShapeSub==kConstSub){
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) PartContMC = GetParticleContainer(1);
        else PartCont = GetTrackContainer(1);
    }
    else{
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) PartContMC = GetParticleContainer(0);
        else PartCont = GetTrackContainer(0);
    }
    
    TClonesArray *TracksArray = NULL;
    TClonesArray *TracksArrayMC = NULL;
    
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) TracksArrayMC = PartContMC->GetArray();
    else TracksArray = PartCont->GetArray();
    
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly){
        if(!PartContMC || !TracksArrayMC) return -2;
    }
    else {
        if(!PartCont || !TracksArray) return -2;
    }
    
    
    AliAODTrack *Track = 0x0;
    Float_t sumpt=0;
    Int_t NTracks=0;
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) NTracks = TracksArrayMC->GetEntriesFast();
    else NTracks = TracksArray->GetEntriesFast();
    
    for(Int_t i=0; i < NTracks; i++){
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly){
            if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
                if (!Track) continue;
                if(TMath::Abs(Track->Eta())>0.9) continue;
                Double_t dphi = RelativePhi(Track->Phi(),jet->Phi());
                Double_t dr2 = (Track->Eta()-jet->Eta())*(Track->Eta()-jet->Eta()) + dphi*dphi;
                Double_t dr = TMath::Sqrt(dr2);
                if((dr>=0.8) && (dr<1)) sumpt=sumpt+Track->Pt();
            }
        }
        else{
            if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
                if (!Track) continue;
                if(TMath::Abs(Track->Eta())>0.9) continue;
                if (Track->Pt()<0.15) continue;
                Double_t dphi = RelativePhi(Track->Phi(),jet->Phi());
                Double_t dr2 = (Track->Eta()-jet->Eta())*(Track->Eta()-jet->Eta()) + dphi*dphi;
                Double_t dr = TMath::Sqrt(dr2);
                if((dr>=0.8) && (dr<1)) sumpt=sumpt+Track->Pt();
                
            }
        }
    }
    return sumpt;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetCoronna(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    if((fJetShapeSub==kDerivSub) && (jetContNb==0)) return -2;
    else
        return Coronna(jet, jetContNb);
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::PTD(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        num=num+vp1->Pt()*vp1->Pt();
        den=den+vp1->Pt();
    }
    return TMath::Sqrt(num)/den;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetpTD(AliEmcalJet *jet, Int_t jetContNb = 0){
    //calc subtracted jet mass
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedpTD();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedpTD();
        else
            return PTD(jet, jetContNb);
    
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::Circularity(AliEmcalJet *jet, Int_t jetContNb = 0){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t mxx    = 0.;
    Double_t myy    = 0.;
    Double_t mxy    = 0.;
    int  nc     = 0;
    Double_t sump2  = 0.;
    Double_t pxjet=jet->Px();
    Double_t pyjet=jet->Py();
    Double_t pzjet=jet->Pz();
    
    
    //2 general normalized vectors perpendicular to the jet
    TVector3  ppJ1(pxjet, pyjet, pzjet);
    TVector3  ppJ3(- pxjet* pzjet, - pyjet * pzjet, pxjet * pxjet + pyjet * pyjet);
    ppJ3.SetMag(1.);
    TVector3  ppJ2(-pyjet, pxjet, 0);
    ppJ2.SetMag(1.);
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        TVector3 pp(vp1->Px(), vp1->Py(), vp1->Pz());
        
        //local frame
        TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
        TVector3 pPerp = pp - pLong;
        //projection onto the two perpendicular vectors defined above
        
        Float_t ppjX = pPerp.Dot(ppJ2);
        Float_t ppjY = pPerp.Dot(ppJ3);
        Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);
        if(ppjT<=0) return 0;
        
        mxx += (ppjX * ppjX / ppjT);
        myy += (ppjY * ppjY / ppjT);
        mxy += (ppjX * ppjY / ppjT);
        nc++;
        sump2 += ppjT;}
    
    if(nc<2) return 0;
    if(sump2==0) return 0;
    // Sphericity Matrix
    Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};
    TMatrixDSym m0(2,ele);
    
    // Find eigenvectors
    TMatrixDSymEigen m(m0);
    TVectorD eval(2);
    TMatrixD evecm = m.GetEigenVectors();
    eval  = m.GetEigenValues();
    // Largest eigenvector
    int jev = 0;
    //  cout<<eval[0]<<" "<<eval[1]<<endl;
    if (eval[0] < eval[1]) jev = 1;
    TVectorD evec0(2);
    // Principle axis
    evec0 = TMatrixDColumn(evecm, jev);
    Double_t compx=evec0[0];
    Double_t compy=evec0[1];
    TVector2 evec(compx, compy);
    Double_t circ=0;
    if(jev==1) circ=2*eval[0];
    if(jev==0) circ=2*eval[1];
    
    return circ;
    
    
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb =0 ){
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedCircularity();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedCircularity();
        else
            return Circularity(jet, jetContNb);
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::LeSub(AliEmcalJet *jet, Int_t jetContNb =0 ){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t den=0.;
    Double_t num = 0.;
    AliVParticle *vp1 = 0x0;
    AliVParticle *vp2 = 0x0;
    std::vector<int> ordindex;
    ordindex=jet->GetPtSortedTrackConstituentIndexes(jetCont->GetParticleContainer()->GetArray());
    //Printf("Nbof const = %d", jet->GetNumberOfTracks());
    //Printf("ordindex[0] = %d, ordindex[1] = %d", ordindex[0], ordindex[1]);
    
    if(ordindex.size()<2) return -1;
    
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[0], jetCont->GetParticleContainer()->GetArray()));
    if (!vp1){
        Printf("AliVParticle associated to Leading constituent not found");
        return -1;
    }
    
    vp2 = static_cast<AliVParticle*>(jet->TrackAt(ordindex[1], jetCont->GetParticleContainer()->GetArray()));
    if (!vp2){
        Printf("AliVParticle associated to Subleading constituent not found");
        return -1;
    }
    
    
    num=vp1->Pt();
    den=vp2->Pt();
    //Printf("vp1->Pt() =%f, vp2->Pt() =%f", vp1->Pt(), vp2->Pt());
    
    return num-den;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb =0) {
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedLeSub();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub();
        else
            return LeSub(jet, jetContNb);
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb=0){
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedConstituent();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedConstituent();
        else
            return jet->GetNumberOfTracks();
    
}


//______________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::Sigma2(AliEmcalJet *jet, Int_t jetContNb=0){
    
    AliJetContainer *jetCont = GetJetContainer(jetContNb);
    if (!jet->GetNumberOfTracks())
        return 0;
    Double_t mxx    = 0.;
    Double_t myy    = 0.;
    Double_t mxy    = 0.;
    int  nc     = 0;
    Double_t sump2  = 0.;
    
    AliVParticle *vp1 = 0x0;
    for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
        vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
        
        if (!vp1){
            Printf("AliVParticle associated to constituent not found");
            continue;
        }
        
        Double_t ppt=vp1->Pt();
        Double_t dphi = RelativePhi(vp1->Phi(),jet->Phi());
        
        Double_t deta = vp1->Eta()-jet->Eta();
        mxx += ppt*ppt*deta*deta;
        myy += ppt*ppt*dphi*dphi;
        mxy -= ppt*ppt*deta*TMath::Abs(dphi);
        nc++;
        sump2 += ppt*ppt;
        
    }
    if(nc<2) return 0;
    if(sump2==0) return 0;
    // Sphericity Matrix
    Double_t ele[4] = {mxx , mxy , mxy , myy };
    TMatrixDSym m0(2,ele);
    
    // Find eigenvectors
    TMatrixDSymEigen m(m0);
    TVectorD eval(2);
    TMatrixD evecm = m.GetEigenVectors();
    eval  = m.GetEigenValues();
    // Largest eigenvector
    int jev = 0;
    //  cout<<eval[0]<<" "<<eval[1]<<endl;
    if (eval[0] < eval[1]) jev = 1;
    TVectorD evec0(2);
    // Principle axis
    evec0 = TMatrixDColumn(evecm, jev);
    Double_t compx=evec0[0];
    Double_t compy=evec0[1];
    TVector2 evec(compx, compy);
    Double_t sig=0;
    if(jev==1) sig=TMath::Sqrt(TMath::Abs(eval[0])/sump2);
    if(jev==0) sig=TMath::Sqrt(TMath::Abs(eval[1])/sump2);
    
    return sig;
    
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalHfeTagging::GetSigma2(AliEmcalJet *jet, Int_t jetContNb=0){
    //calc subtracted jet mass
    
    if((fJetShapeSub==kDerivSub)&&(jetContNb==0))
        if (fDerivSubtrOrder == 1) return jet->GetShapeProperties()->GetFirstOrderSubtractedSigma2();
        else return jet->GetShapeProperties()->GetSecondOrderSubtractedSigma2();
        else
            return Sigma2(jet, jetContNb);
    
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalHfeTagging::SelectTrigger(Float_t minpT, Float_t maxpT){
    
    AliTrackContainer *PartCont = NULL;
    AliParticleContainer *PartContMC = NULL;
    
    
    if (fJetShapeSub==kConstSub){
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) PartContMC = GetParticleContainer(1);
        else PartCont = GetTrackContainer(1);
    }
    else{
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) PartContMC = GetParticleContainer(0);
        else PartCont = GetTrackContainer(0);
    }
    
    TClonesArray *TracksArray = NULL;
    TClonesArray *TracksArrayMC = NULL;
    
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) TracksArrayMC = PartContMC->GetArray();
    else TracksArray = PartCont->GetArray();
    
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly){
        if(!PartContMC || !TracksArrayMC) return -99999;
    }
    else {
        if(!PartCont || !TracksArray) return -99999;
    }
    
    
    AliAODTrack *Track = 0x0;
    
    
    
    //TList *trackList = new TList();
    Int_t triggers[100];
    for (Int_t iTrigger=0; iTrigger<100; iTrigger++) triggers[iTrigger] = 0;
    Int_t iTT = 0;
    Int_t NTracks=0;
    if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly) NTracks = TracksArrayMC->GetEntriesFast();
    else NTracks = TracksArray->GetEntriesFast();
    
    for(Int_t i=0; i < NTracks; i++){
        if (fJetShapeType == AliAnalysisTaskEmcalHfeTagging::kGenOnTheFly){
            if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
                if (!Track) continue;
                if(TMath::Abs(Track->Eta())>0.9) continue;
                if (Track->Pt()<0.15) continue;
                if ((Track->Pt() >= minpT) && (Track->Pt()< maxpT)) {
                    triggers[iTT] = i;
                    iTT++;
                }
            }
        }
        else{
            if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
                if (!Track) continue;
                if(TMath::Abs(Track->Eta())>0.9) continue;
                if (Track->Pt()<0.15) continue;
                if ((Track->Pt() >= minpT) && (Track->Pt()< maxpT)) {
                    triggers[iTT] = i;
                    iTT++;
                }
            }
        }
    }
    
    
    if (iTT == 0) return -99999;
    Int_t nbRn = 0, index = 0 ;
    TRandom3* random = new TRandom3(0);
    nbRn = random->Integer(iTT);
    index = triggers[nbRn];
    //Printf("iTT Total= %d, nbRn = %d, Index = %d",iTT, nbRn, index );
    return index;
    
}

//__________________________________________________________________________________
Double_t AliAnalysisTaskEmcalHfeTagging::RelativePhi(Double_t mphi,Double_t vphi){
    
    if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
    else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
    if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
    else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
    double dphi = mphi-vphi;
    if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
    else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
    return dphi;//dphi in [-Pi, Pi]
}


//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::RetrieveEventObjects() {
    //
    // retrieve event objects
    //
    if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
        return kFALSE;
    
    return kTRUE;
}

//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::InclElecTrackCuts(const AliVVertex *pVietx,AliAODTrack *ietrack,Int_t nMother, Double_t listMother[])
{
    // track cuts for inclusive electrons
    
    if(!ietrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
    if(TMath::Abs(ietrack->Eta())>0.8) return kFALSE;
    if(ietrack->GetTPCNcls() < fTPCnCut) return kFALSE;
    if (ietrack->GetITSNcls() < fITSncut) return kFALSE;
    if(!ietrack->IsOn(AliAODTrack::kITSrefit)) return kFALSE;
    if(!ietrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
    if(!(ietrack->HasPointOnITSLayer(0) && ietrack->HasPointOnITSLayer(1))) return kFALSE;
    if ((ietrack->Pt()<0.5) || (ietrack->Pt()>4)) return kFALSE;
    
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < nMother; kinkmother++) {
        if(ietrack->GetID() == listMother[kinkmother]) {
            kinkmotherpass = kFALSE;
            continue;
        }
    }
    if(!kinkmotherpass) return kFALSE;
    
    Double_t d0z0[2]={-999,-999}, cov[3];
    
    if(ietrack->PropagateToDCA(pVietx, fVevent->GetMagneticField(), 20., d0z0, cov))
        if(TMath::Abs(d0z0[0]) > fDcaXYcut || TMath::Abs(d0z0[1]) > fDcaZcut) return kFALSE;
    
    return kTRUE;
    
}
//_________________________________________
Bool_t AliAnalysisTaskEmcalHfeTagging::PhotElecTrackCuts(const AliVVertex *pVaetx,AliAODTrack *aetrack,Int_t nMother, Double_t listMother[])
{
    // track cuts for associate tracks of photonic electrons
    
    if(!aetrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) return kFALSE;
    if(aetrack->Pt() < fAssPtCut) return kFALSE;
    if(TMath::Abs(aetrack->Eta())>0.9) return kFALSE;
    if(aetrack->GetTPCNcls() < fAssTPCnCut) return kFALSE;
    if (fAssITSrefitCut && !(aetrack->GetStatus()&AliAODTrack::kITSrefit)) return kFALSE;
    if(!(aetrack->GetStatus()&AliAODTrack::kTPCrefit)) return kFALSE;
    
    return kTRUE;
    
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalHfeTagging::Terminate(Option_t *)
{
    // Called once at the end of the analysis.
    
    // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
    // if (!fTreeObservableTagging){
    //   Printf("ERROR: fTreeObservableTagging not available");
    //   return;
    // }
    
}

