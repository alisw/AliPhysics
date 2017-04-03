//includes
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "TVector3.h"
#include "Math/Vector3Dfwd.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TDirectory.h"
#include <vector>
#include <iostream>
#include "TParticle.h"
#include "TPDGCode.h"
#include "TDatabasePDG.h"
#include "AliPID.h"
#include "AliLog.h"
#include "AliPIDResponse.h"
#include "AliConversionPhotonCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliExternalTrackParam.h"
#include "AliV0ReaderV1.h"
#include "AliESDInputHandler.h"
#include "AliVEventHandler.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDtrackCuts.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliAnalysisTaskDirectPhotons.h"

using namespace std;

ClassImp( AliAnalysisTaskDirectPhotons );

//--------------------------------------------------------------------------------------------------
AliAnalysisTaskDirectPhotons::AliAnalysisTaskDirectPhotons():
    AliAnalysisTaskSE(),
    fOutputList(NULL),
    fDataIsElectron(NULL),
    fDataIsPositron(NULL),
    fPhotonCandidateFirst(NULL),
    fGoodPhotonCandidates(NULL),
    fGammaCutArray(NULL),
    fEventCutArray(NULL),
    fNeutralPionMesonCutArray(NULL),
    pSigmaTPC(NULL),
    fnCuts(0),
    magField(0),
    xPrimVertex(0),
    yPrimVertex(0),
    zPrimVertex(0),
    sigmaVertex(0),
    dPointToSecVertexLine(0),
    nTracks(0),
    fisMC(kTRUE),
    nPhotonCandidates(0),
    fEvent(NULL),
    fPrimVertex(NULL),
    fESDEvent(NULL),
    fAODEvent(NULL),
    fTrackCuts(new AliESDtrackCuts),
    fV0Reader(NULL),
    fMCEvent(0),
    fMCStack(NULL),
    fPIDResponse(NULL),
    trueTrack(NULL),
    myParticle(NULL),
    fConversionCuts(NULL),
    fHistZVertex(NULL),
    fHistXY(NULL),
    fHistPtNegativeParticles(NULL),
    fHistPtPositiveParticles(NULL),
    fHistInvMassPosElec(NULL),
    fHistInvMassPosElecBack(NULL),
    fHistInvMassPosElecAfterDCut(NULL),
    fHistDistanceToPrimVertex(NULL),
    fHistTrueDCA(NULL),
    fHistTrueInvMassPosElecNoCut(NULL),
    fHistTrueInvMassPosElecDCACut(NULL),
    fHistTrueInvMassPosElecDCAvertexCut(NULL),
    fHistNSigmasTPCElectron1(NULL),
    fHistNSigmasTPCElectron2(NULL),
    fHistNSigmasTPCElectron3(NULL),
    fHistNSigmasTPCElectron4(NULL),
    fHistNSigmasTPCElectron5(NULL),
    fHistNSigmasTPCElectron6(NULL),
    fHistInvMassNeutralMeson(NULL),
    fHistTrueMesonPt(NULL),
    fHistInvMassBackNeutralMeson(NULL),
    fHistNSigmasMomentumTPC(NULL)

{

}

//--------------------------------------------------------------------------------------------------
AliAnalysisTaskDirectPhotons::AliAnalysisTaskDirectPhotons(const char* name):
    AliAnalysisTaskSE(name),
    fOutputList(NULL),
    fDataIsElectron(NULL),
    fDataIsPositron(NULL),
    fPhotonCandidateFirst(NULL),
    fGoodPhotonCandidates(NULL),
    fGammaCutArray(NULL),
    fEventCutArray(NULL),
    fNeutralPionMesonCutArray(NULL),
    pSigmaTPC(NULL),
    fnCuts(0),
    magField(0),
    xPrimVertex(0),
    yPrimVertex(0),
    zPrimVertex(0),
    sigmaVertex(0),
    dPointToSecVertexLine(0),
    nTracks(0),
    fisMC(kTRUE),
    nPhotonCandidates(0),
    fEvent(NULL),
    fPrimVertex(NULL),
    fESDEvent(NULL),
    fAODEvent(NULL),
    fTrackCuts(new AliESDtrackCuts),
    fV0Reader(NULL),
    fMCEvent(0),
    fMCStack(NULL),
    fPIDResponse(NULL),
    trueTrack(NULL),
    myParticle(NULL),
    fConversionCuts(NULL),
    fHistZVertex(NULL),
    fHistXY(NULL),
    fHistPtNegativeParticles(NULL),
    fHistPtPositiveParticles(NULL),
    fHistInvMassPosElec(NULL),
    fHistInvMassPosElecBack(NULL),
    fHistInvMassPosElecAfterDCut(NULL),
    fHistDistanceToPrimVertex(NULL),
    fHistTrueDCA(NULL),
    fHistTrueInvMassPosElecNoCut(NULL),
    fHistTrueInvMassPosElecDCACut(NULL),
    fHistTrueInvMassPosElecDCAvertexCut(NULL),
    fHistNSigmasTPCElectron1(NULL),
    fHistNSigmasTPCElectron2(NULL),
    fHistNSigmasTPCElectron3(NULL),
    fHistNSigmasTPCElectron4(NULL),
    fHistNSigmasTPCElectron5(NULL),
    fHistNSigmasTPCElectron6(NULL),
    fHistInvMassNeutralMeson(NULL),
    fHistTrueMesonPt(NULL),
    fHistInvMassBackNeutralMeson(NULL),
    fHistNSigmasMomentumTPC(NULL)

{
  //  DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());

/*    fTrackCuts->SetMinNClustersTPC(70);

    fTrackCuts->SetMaxChi2PerClusterTPC(4);

    fTrackCuts->SetAcceptKinkDaughters(kFALSE);

    //fTrackCuts->SetRequireTPCRefit(kTRUE);

    //fTrackCuts->SetRequireITSRefit(kTRUE);

    fTrackCuts->SetMaxDCAToVertexXY(3.0);

    fTrackCuts->SetMaxDCAToVertexZ(3.0);

    fTrackCuts->SetDCAToVertex2D(kFALSE);

    fTrackCuts->SetRequireSigmaToVertex(kFALSE);

    //fTrackCuts->SetPtRange(0.2, 1e10);

    fTrackCuts->SetEtaRange(-0.8,0.8);*/

}

//--------------------------------------------------------------------------------------------------
AliAnalysisTaskDirectPhotons::~AliAnalysisTaskDirectPhotons()
{

    //--------------------------> Do que faço o delete ou não no destructor?

    if(fOutputList){
        delete fOutputList;
        fOutputList = 0x0;
    }

    if(fPIDResponse){
        delete fPIDResponse;
        fPIDResponse = 0x0;
    }

    /*if(pSigmaTPC){//----------------------------------->faço esse delete ou não?
        delete[] pSigmaTPC;
        pSigmaTPC = 0x0;
    }*/
}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);

    fDataIsElectron = new TClonesArray("AliVTrack");

    fDataIsPositron = new TClonesArray("AliVTrack");

    fPhotonCandidateFirst = new TClonesArray("AliAODConversionPhoton");

    fGoodPhotonCandidates = new TClonesArray("AliAODConversionPhoton");

    fHistZVertex = new TH1D("fHistZVertex","fHistZVertex", 400, -20, 20);
    fOutputList->Add(fHistZVertex);

    fHistXY = new TH2D("fHistXY","fHistXY", 250, -5, 5, 250, -5, 5);
    fOutputList->Add(fHistXY);

    fHistPtNegativeParticles = new TH1F("fHistPtNeg","fHistPtNeg", 500, 0, 50);
    fOutputList->Add(fHistPtNegativeParticles);

    fHistPtPositiveParticles = new TH1F("fHistPtPos","fHistPtPos", 500, 0, 50);
    fOutputList->Add(fHistPtPositiveParticles);   

    fHistInvMassPosElec = new TH2D("fHistInvMassPosElec","fHistInvMassPosElec", 125, 0, 0.25, 100, 0, 10);
    fOutputList->Add(fHistInvMassPosElec);

    fHistInvMassPosElecBack = new TH1F("fHistInvMassPosElecBack","fHistInvMassPosElecBack", 125, 0, 0.25);
    fOutputList->Add(fHistInvMassPosElecBack);

    fHistInvMassPosElecAfterDCut = new TH1F("fHistInvMassPosElecAfterDCut","fHistInvMassPosElecAfterDCut", 125, 0, 0.25);
    fOutputList->Add(fHistInvMassPosElecAfterDCut);

    fHistTrueDCA = new TH1D("fHistTrueDCA","fHistTrueDCA", 100, 0, 5);
    fOutputList->Add(fHistTrueDCA);

    fHistTrueInvMassPosElecNoCut = new TH1F("fHistTrueInvMassPosElecNoCut","fHistTrueInvMassPosElecNoCut", 125, 0, 0.25);
    fOutputList->Add(fHistTrueInvMassPosElecNoCut);

    fHistTrueInvMassPosElecDCACut = new TH1F("fHistTrueInvMassPosElecDCACut","fHistTrueInvMassPosElecDCACut", 125, 0, 0.25);
    fOutputList->Add(fHistTrueInvMassPosElecDCACut);

    fHistTrueInvMassPosElecDCAvertexCut = new TH1F("fHistTrueInvMassPosElecDCAvertexCut","fHistTrueInvMassPosElecDCAvertexCut", 125, 0, 0.25);
    fOutputList->Add(fHistTrueInvMassPosElecDCAvertexCut);

    fHistDistanceToPrimVertex= new TH1D("fHistDistanceToPrimVertex","fHistDistanceToPrimVertex", 1500, 0, 15);
    fOutputList->Add(fHistDistanceToPrimVertex);

    fHistNSigmasTPCElectron1 = new TH1F("fHistnSigmaElec1","fHistnSigmaElec1", 200, -8, 8);
    fOutputList->Add(fHistNSigmasTPCElectron1);

    fHistNSigmasTPCElectron2 = new TH1F("fHistnSigmaElec2","fHistnSigmaElec2", 200, -8, 8);
    fOutputList->Add(fHistNSigmasTPCElectron2);

    fHistNSigmasTPCElectron3 = new TH1F("fHistnSigmaElec3","fHistnSigmaElec3", 200, -8, 8);
    fOutputList->Add(fHistNSigmasTPCElectron3);

    fHistNSigmasTPCElectron4 = new TH1F("fHistnSigmaElec4","fHistnSigmaElec4", 200, -8, 8);
    fOutputList->Add(fHistNSigmasTPCElectron4);

    fHistNSigmasTPCElectron5 = new TH1F("fHistnSigmaElec5","fHistnSigmaElec5", 200, -8, 8);
    fOutputList->Add(fHistNSigmasTPCElectron5);

    fHistNSigmasTPCElectron6 = new TH1F("fHistnSigmaElec6","fHistnSigmaElec6", 200, -8, 8);
    fOutputList->Add(fHistNSigmasTPCElectron6);

    fHistInvMassNeutralMeson = new TH2D("fHistInvMassNeutralMeson","fHistInvMassNeutralMeson", 200, 0, 1.0, 100, 0, 10);
    fOutputList->Add(fHistInvMassNeutralMeson);

    fHistTrueMesonPt = new TH1D("fHistTrueMesonPt","fHistTrueMesonPt", 500, 0, 50);
    fOutputList->Add(fHistTrueMesonPt);

    fHistInvMassBackNeutralMeson = new TH1D("fHistInvMassBackNeutralMeson","fHistInvMassBackNeutralMeson", 200, 0, 1.0);
    fOutputList->Add(fHistInvMassBackNeutralMeson);

    fHistNSigmasMomentumTPC = new TH2F("fHistNSigmaP","fHistNSigmaP", 250, 0, 5, 400, -8, 8);
    fOutputList->Add(fHistNSigmasMomentumTPC);

    PostData(1, fOutputList);

}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::UserExec(Option_t *)
{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (mgr == NULL) return;

    AliESDInputHandler *esdHandler = (AliESDInputHandler*)mgr->GetInputEventHandler();
    if (esdHandler == NULL) return;

    AliMCEventHandler *mcHandler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
    if (mcHandler == NULL) return;

    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
    if(!fV0Reader) {
        cout<<"v0 Reader is missing"<<endl;
        return;
    }

    fConversionCuts = fV0Reader->GetConversionCuts();

    fPIDResponse = (AliPIDResponse*)esdHandler->GetPIDResponse();
    if (fPIDResponse == NULL) return;

    fESDEvent = dynamic_cast<AliESDEvent*>(InputEvent());
    if (fESDEvent == NULL) return;

    if(!fESDEvent){
        fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
        if (fAODEvent == NULL) return;
    }

    if (fESDEvent) fEvent = dynamic_cast<AliVEvent*>(fESDEvent);
    if (fAODEvent) fEvent = dynamic_cast<AliVEvent*>(fAODEvent);

    fMCEvent = MCEvent();

    if(fMCEvent){
        fMCStack = (AliStack*)fMCEvent->Stack();
    }

    fPrimVertex = (AliVVertex*)fEvent->GetPrimaryVertex();

    //Double_t* covmatrix = new Double_t[6];
    //fPrimVertex->GetCovarianceMatrix(covmatrix);

    xPrimVertex = fPrimVertex->GetX();
    yPrimVertex = fPrimVertex->GetY();
    zPrimVertex = fPrimVertex->GetZ();

    fHistXY->Fill(xPrimVertex,yPrimVertex);
    fHistZVertex->Fill(zPrimVertex);

    /*sigmaVertex = TMath::Sqrt((xPrimVertex*xPrimVertex*covmatrix[0]*covmatrix[0]) + (yPrimVertex*yPrimVertex*covmatrix[2]*covmatrix[2]) +
            (zPrimVertex*zPrimVertex*covmatrix[5]*covmatrix[5]) + (2*xPrimVertex*yPrimVertex*covmatrix[1]*covmatrix[1])+
            (2*xPrimVertex*zPrimVertex*covmatrix[3]*covmatrix[3]) + (2*yPrimVertex*zPrimVertex*covmatrix[4]*covmatrix[4]))/
            TMath::Sqrt(xPrimVertex*xPrimVertex + yPrimVertex*yPrimVertex + zPrimVertex*zPrimVertex);

    cout<<sigmaVertex<<endl;*/

    magField = fEvent->GetMagneticField();

    nTracks = fEvent->GetNumberOfTracks();

    pSigmaTPC = new Float_t[nTracks];

    trueTrack = new AliVTrack*[nTracks]; //--------------------------------->Este método está correto?

    /*if (fMCEvent && fMCStack) {
        CalculateEfficiency();
    }*/

    ProcessAllTracks();
    FindGammaCandidates();
    FindNeutralMesonCandidates();
    //FindTrueNeutralMesonCandidates();

    fDataIsElectron->Clear();
    fDataIsPositron->Clear();
    fPhotonCandidateFirst->Clear();
    fGoodPhotonCandidates->Clear();

    PostData(1, fOutputList);

    delete[] pSigmaTPC;
    delete[] trueTrack;
    //delete fV0Reader;

}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::CalculateEfficiency(){

    ProcessStackParticles();
    ProcessDataEfficiency();

}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::ProcessStackParticles(){

    Int_t stackNTracks = fMCStack->GetNtrack();
    //cout<<"n tracks stack"<<stackNTracks<<endl;

    for (Int_t iStackTrack = 0; iStackTrack < stackNTracks; iStackTrack++){

        myParticle = new AliMyTParticleMotherAndDaughter();

        TParticle* particle = fMCStack->Particle(iStackTrack);
        if (!particle) continue;

        Int_t pdg = particle->GetPdgCode();
        if(pdg != 22) continue;

        Int_t labelDaughterOne = particle->GetFirstDaughter();

        Int_t labelDaughterTwo = particle->GetLastDaughter();

        if(labelDaughterTwo != labelDaughterOne+1) continue;

        TParticle* daughterOne = fMCStack->Particle(labelDaughterOne);
        if(!daughterOne) continue;

        Int_t pdgDaughterOne = daughterOne->GetPdgCode();
        if(TMath::Abs(pdgDaughterOne)!= 11) continue;

        TParticle* daughterTwo = fMCStack->Particle(labelDaughterTwo);
        if(!daughterTwo) continue;

        Int_t pdgDaughterTwo = daughterTwo->GetPdgCode();
        if(TMath::Abs(pdgDaughterTwo)!= 11) continue;

        if (pdgDaughterOne == 11 && pdgDaughterTwo == -11) {

            myParticle->labelDaughterElectron= labelDaughterOne;
            myParticle->labelDaughterPositron= labelDaughterTwo;
            myParticle->labelMotherPhoton= iStackTrack; 
        }

        if (pdgDaughterTwo == 11 && pdgDaughterOne == -11) {

            myParticle->labelDaughterElectron= labelDaughterTwo;
            myParticle->labelDaughterPositron= labelDaughterOne;
            myParticle->labelMotherPhoton= iStackTrack;

        }

        ProcessStackInfos();

        delete myParticle;

    }

}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::ProcessStackInfos(){

    Int_t idxelec=0, idxpos=0;

    for (Int_t iTracks = 0 ; iTracks < nTracks ; iTracks++){

        trueTrack[iTracks] = (AliVTrack*)fEvent->GetTrack(iTracks);
        if (!trueTrack[iTracks]) continue;

        Int_t labelTrack = trueTrack[iTracks]->GetLabel();

        if ( labelTrack == myParticle->labelDaughterElectron ){

            //fDataIsElectron->Add(trueTrack[iTracks]);

            AliVTrack *elec = static_cast<AliVTrack*>(fDataIsElectron->New(idxelec++));
            elec = trueTrack[iTracks];

        }

        if ( labelTrack == myParticle->labelDaughterPositron ){

            //fDataIsPositron->Add(trueTrack[iTracks]);

            AliVTrack *pos = static_cast<AliVTrack*>(fDataIsPositron->New(idxpos++));
            pos = trueTrack[iTracks];

        }

    }

}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::ProcessDataEfficiency(){

    Double_t maxd=2, E1=-999., E2=-999., *dz=0, *cov=0;
    Double_t x1=-999., x2=-999., elecMass=0.511e-3, photonMass=-999;

    TLorentzVector firstParticle(0,0,0,0), secondParticle(0,0,0,0), firstPhotonCandidate(0,0,0,0);

    Int_t nElectrons = fDataIsElectron->GetEntries();
    //cout<<"n electrons"<<nElectrons<<endl;

    Int_t nPositrons = fDataIsPositron->GetEntries();
    //cout<<"n positrons"<<nPositrons<<endl;

    for (Int_t iElectrons=0; iElectrons < nElectrons; iElectrons++){

        AliVTrack* electronTrack = (AliVTrack*)fDataIsElectron->At(iElectrons);
        if (!electronTrack) continue;

        AliExternalTrackParam* aliExternal1 = new AliExternalTrackParam();
        Double_t* p1 = new Double_t[3];

        aliExternal1 -> CopyFromVTrack(electronTrack);
        if(!aliExternal1) continue;

        for(Int_t iPositrons = 0; iPositrons < nPositrons; iPositrons++){

            AliVTrack* positronTrack = (AliVTrack*)fDataIsPositron->At(iPositrons);
            if (!positronTrack) continue;

            AliExternalTrackParam* aliExternal2 = new AliExternalTrackParam();
            Double_t* p2 = new Double_t[3];

            if(electronTrack == positronTrack) continue;

            aliExternal2 -> CopyFromVTrack(positronTrack);
            if(!aliExternal2) continue;

            //if(positronTrack->GetMother(0) != electronTrack->GetMother(0)) continue; ----------->NAO POSSO PQ NAO EH TPARTICLE

            //if(positronTrack->GetLabel() != electronTrack->GetLabel()+1 || electronTrack->GetLabel() != positronTrack->GetLabel()+1) continue;

            Double_t dca = aliExternal1->GetDCA(aliExternal2, magField, x1, x2);
            fHistTrueDCA->Fill(dca);

            aliExternal1->GetPxPyPzAt(x1, magField, p1);
            aliExternal2->GetPxPyPzAt(x2, magField, p2);

            E1= TMath::Sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]+elecMass*elecMass);
            E2= TMath::Sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]+elecMass*elecMass);

            firstParticle.SetPxPyPzE(p1[0],p1[1],p1[2],E1);
            secondParticle.SetPxPyPzE(p2[0],p2[1],p2[2],E2);

            firstPhotonCandidate = firstParticle + secondParticle;
            photonMass = firstPhotonCandidate.M();

            fHistTrueInvMassPosElecNoCut->Fill(photonMass);

            if (dca > 2) continue;
            fHistTrueInvMassPosElecDCACut->Fill(photonMass);

            if (aliExternal1->PropagateToDCA(fPrimVertex, magField, maxd, dz, cov)) continue;
            fHistTrueInvMassPosElecDCAvertexCut->Fill(photonMass);

            delete aliExternal2;
            delete[] p2;

        }

        delete aliExternal1;
        delete[] p1;

    }


}


//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::ProcessAllTracks(){

    for (Int_t iTracks = 0; iTracks < nTracks ; iTracks++){

        AliVTrack* track = (AliVTrack*)fEvent->GetTrack(iTracks);
        if (!track) continue;
        //if (!fTrackCuts->AcceptTrack(track)) continue;

        if (track->Charge() < 0){
            fHistPtNegativeParticles->Fill(track->Pt());
        }

        if (track->Charge() > 0){
            fHistPtPositiveParticles->Fill(track->Pt());
        }

        pSigmaTPC[iTracks] =fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        if(!&pSigmaTPC[iTracks]) continue;
        if(pSigmaTPC[iTracks] == -999) continue;

        Double_t trackMom=track->P();

        fHistNSigmasMomentumTPC->Fill(track->P(),pSigmaTPC[iTracks]);

        if (trackMom>0.25 && trackMom<0.5) fHistNSigmasTPCElectron1->Fill(pSigmaTPC[iTracks]);

        if (trackMom>0.5 && trackMom<0.75) fHistNSigmasTPCElectron2->Fill(pSigmaTPC[iTracks]);

        if (trackMom>0.75 && trackMom<1.0) fHistNSigmasTPCElectron3->Fill(pSigmaTPC[iTracks]);

        if (trackMom>1.0 && trackMom<1.5) fHistNSigmasTPCElectron4->Fill(pSigmaTPC[iTracks]);

        if (trackMom>1.5 && trackMom<2.0) fHistNSigmasTPCElectron5->Fill(pSigmaTPC[iTracks]);

        if (trackMom>2.0 && trackMom<5.0) fHistNSigmasTPCElectron6->Fill(pSigmaTPC[iTracks]);

    }

}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::FindGammaCandidates(){

    Int_t idx = 0;

    Double_t photonMass = -999., photonPt = -999.;

    fPhotonCandidateFirst = (TClonesArray*)fV0Reader->GetReconstructedGammas();

    nPhotonCandidates = fPhotonCandidateFirst->GetEntries();

    for (Int_t iPhotons=0; iPhotons < nPhotonCandidates; iPhotons++){

        AliAODConversionPhoton* firstPhotonCandidate = static_cast<AliAODConversionPhoton*>(fPhotonCandidateFirst->At(iPhotons));
        if (!firstPhotonCandidate) continue;

        if(!(AliConversionPhotonCuts*)fConversionCuts->PhotonIsSelected(firstPhotonCandidate,fEvent)) continue;

        photonMass = firstPhotonCandidate->M();
        photonPt = firstPhotonCandidate->Pt();

        if (photonMass > 0.06) continue; //------------->CORTE NO LUGAR CORRETO?

        AliAODConversionPhoton* photonCandidate = static_cast<AliAODConversionPhoton*>(fGoodPhotonCandidates->New(idx++));

        fHistInvMassPosElec->Fill(photonMass, photonPt);

        //fHistInvMassPosElecBack->Fill(photonMass);

    }

}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::FindNeutralMesonCandidates(){

    Double_t neutralMesonMass=-999., neutralMesonPt=-999.;

    for (Int_t iPhotonEntries=0; iPhotonEntries < nPhotonCandidates-1; iPhotonEntries++){

        if (fPhotonCandidateFirst->At(iPhotonEntries)==NULL) continue;

        AliAODConversionPhoton* firstPhotonCandidate = dynamic_cast<AliAODConversionPhoton*>(fPhotonCandidateFirst->At(iPhotonEntries));
        if(!firstPhotonCandidate) continue;

        for(Int_t jPhotonEntries = iPhotonEntries; jPhotonEntries < nPhotonCandidates; jPhotonEntries++){

            if (fPhotonCandidateFirst->At(jPhotonEntries)==NULL) continue;

            AliAODConversionPhoton* secondPhotonCandidate = dynamic_cast<AliAODConversionPhoton*>(fPhotonCandidateFirst->At(jPhotonEntries));
            if(!secondPhotonCandidate) continue;

            if(firstPhotonCandidate->GetTrackLabelPositive() == secondPhotonCandidate->GetTrackLabelPositive() ||
            firstPhotonCandidate->GetTrackLabelNegative() == secondPhotonCandidate->GetTrackLabelNegative() ||
            firstPhotonCandidate->GetTrackLabelNegative() == secondPhotonCandidate->GetTrackLabelPositive() ||
            firstPhotonCandidate->GetTrackLabelPositive() == secondPhotonCandidate->GetTrackLabelNegative() ) continue;

            AliAODConversionMother* neutralMesonCandidate = new AliAODConversionMother(firstPhotonCandidate,secondPhotonCandidate);

            neutralMesonMass = neutralMesonCandidate->M();
            neutralMesonPt = neutralMesonCandidate->Pt();

            if (neutralMesonMass > 1.0) continue;

            fHistInvMassNeutralMeson->Fill(neutralMesonMass, neutralMesonPt);

            delete neutralMesonCandidate;

        }

    }

}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::FindTrueNeutralMesonCandidates(){

    Int_t stackNTracks = fMCStack->GetNtrack();

    Double_t particlePt = -999.;

    for (Int_t iStackTrack = 0; iStackTrack < stackNTracks; iStackTrack++){

        TParticle* particle = fMCStack->Particle(iStackTrack);
        if (!particle) continue;

        Int_t pdg = particle->GetPdgCode();
        if(pdg != 111) continue; //-------->pdg do meson (pi0 ou eta)

        /*|| pdg != 221*/

        particlePt = particle->Pt();

        fHistTrueMesonPt->Fill(particlePt);

        /*Int_t labelDaughterOne = particle->GetFirstDaughter();

        Int_t labelDaughterTwo = particle->GetLastDaughter();

        if(labelDaughterTwo != labelDaughterOne+1) continue;

        TParticle* daughterOne = fMCStack->Particle(labelDaughterOne);
        if(!daughterOne) continue;

        Int_t pdgDaughterOne = daughterOne->GetPdgCode();
        if(TMath::Abs(pdgDaughterOne)!= 22) continue;  //--------------> filha tem que ser fóton

        TParticle* daughterTwo = fMCStack->Particle(labelDaughterTwo);
        if(!daughterTwo) continue;

        Int_t pdgDaughterTwo = daughterTwo->GetPdgCode();
        if(TMath::Abs(pdgDaughterTwo)!= 22) continue;

        AliMyTParticleMotherAndDaughter* myParticle = new AliMyTParticleMotherAndDaughter();

        if (pdg == 111) {

            myParticle->daughterElectron = daughterOne;
            myParticle->labelDaughterElectron= labelDaughterOne;

            myParticle->daughterPositron = daughterTwo;
            myParticle->labelDaughterPositron= labelDaughterTwo;

            myParticle->motherPhoton = particle;
            myParticle->labelMotherPhoton= iStackTrack; //------------------------------------->é de fato istacktrack?

        }

        if (pdg == 221) {

            myParticle->daughterElectron = daughterTwo;
            myParticle->labelDaughterElectron= labelDaughterTwo;

            myParticle->daughterPositron = daughterOne;
            myParticle->labelDaughterPositron= labelDaughterOne;

            myParticle->motherPhoton = particle;
            myParticle->labelMotherPhoton= iStackTrack;

        }*/

    }

}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::FindSecVertexPosition(TLorentzVector *firstPhotonCandidate){

    Double_t t;

    //Vetor Diretor da Reta (momento dos TLorentz)
    Double_t px = firstPhotonCandidate->Px();
    Double_t py = firstPhotonCandidate->Py();
    Double_t pz = firstPhotonCandidate->Pz();
    TVector3 v(px,py,pz);

    //Coordenadas vértice secundário
    Double_t xSecVertex = firstPhotonCandidate->X();
    Double_t ySecVertex = firstPhotonCandidate->Y();
    Double_t zSecVertex = firstPhotonCandidate->Z();
    TVector3 xyzSec(xSecVertex, ySecVertex, zSecVertex);

    //Coordenadas vértice primário
    TVector3 xyzPrim(xPrimVertex, yPrimVertex, zPrimVertex);

    //Vetor do vértice primário ao vértice secundário
    TVector3 PSecPPrim(xyzSec-xyzPrim);

    dPointToSecVertexLine = ((PSecPPrim.Cross(v)).Dot((PSecPPrim.Cross(v))))/(v.Dot(v));
    fHistDistanceToPrimVertex->Fill(dPointToSecVertexLine);

}

//--------------------------------------------------------------------------------------------------
void AliAnalysisTaskDirectPhotons::Terminate(Option_t *){
        AliAnalysisTaskSE::Terminate();
}

