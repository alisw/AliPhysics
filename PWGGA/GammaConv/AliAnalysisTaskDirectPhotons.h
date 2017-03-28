#ifndef AliAnalysisTaskDirectPhotons_h
#define AliAnalysisTaskDirectPhotons_h

/*#include <iostream>
#include <vector>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "Math/Vector3D.h"*/

class TH1;
class TH1F;
class TH1D;
class TH2;
class TH2F;
class TH2D;
class TLorentzVector;
class TObject;
class TVector3;
class TChain;
class TString;
class TClonesArray;
class TList;
class TObjArray;
class TParticle;
class AliPID;
class AliStack;
class AliPIDResponse;
class AliConversionPhotonCuts;
class AliVEvent;
class AliVVertex;
class AliVTrack;
class AliVParticle;
class AliExternalTrackParam;
class AliESDInputHandler;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDv0;
class AliVEventHandler;
class AliV0ReaderV1;
class AliAODInputHandler;
class AliAODEvent;
class AliAODTrack;
class AliAODv0;
class AliMCEventHandler;
class AliMCEvent;
class AliKFParticle;
class AliKFVertex;
class AliMyTParticleMotherAndDaughter;

//Header File
#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskDirectPhotons: public AliAnalysisTaskSE
{
public:

    AliAnalysisTaskDirectPhotons();
    AliAnalysisTaskDirectPhotons(const char *name);
    virtual             ~AliAnalysisTaskDirectPhotons();

    virtual void        UserCreateOutputObjects();
    virtual void        UserExec(Option_t* option);
    virtual void        Terminate(Option_t* );

    void                CalculateEfficiency();
    void                ProcessStackParticles();
    void                ProcessStackInfos();
    void                ProcessDataEfficiency();
    void                ProcessAllTracks();
    void                FindGammaCandidates();
    void                FindSecVertexPosition(TLorentzVector* firstPhotonCandidate);
    void                FindNeutralMesonCandidates();
    void                FindTrueNeutralMesonCandidates();


    void                SetMCEvent(Bool_t isMC){fisMC=isMC;}
    void                SetV0Reader(AliV0ReaderV1 * reader) { fV0Reader = reader; }
    void                SetEventCutList(Int_t nCuts, TList *CutArray){
                        fnCuts= nCuts;
                        fEventCutArray = CutArray;
                        }
    void                SetConversionCutList(TList *CutArray){ fGammaCutArray = CutArray;}
    void                SetNeutralPionCutList(TList *CutArray){ fNeutralPionMesonCutArray = CutArray; }
    //void                SetV0ReaderName(TString name){fV0ReaderName=name; return;}


private:

    //all lists
    TList*                         fOutputList;
    /*TObjArray*                  fNegativeParticles;
    TObjArray*                  fPositiveParticles;
    TObjArray*                  fNeutralParticles;
    TObjArray*                  fTrueElectrons;
    TObjArray*                  fTruePositrons;*/
    TClonesArray*                  fDataIsElectron;
    TClonesArray*                  fDataIsPositron;
    TClonesArray*                  fPhotonCandidateFirst;
    TClonesArray*                  fGoodPhotonCandidates;
    TList*                         fEventCutArray;
    TList* 	       				fGammaCutArray;
    TList* 	       				fNeutralPionMesonCutArray;


//    TClonesArray*               fParticlesArray;
//    TList*                  fNumSigmasTPCTrack;
//    TObjArray*                  fMyParticles;
//    TList*                  fLabelDaughterElectron;
//    TList*                  fLabelDaughterPositron;


    Float_t*               pSigmaTPC;
    Int_t                  fnCuts;
    Double_t               magField;
    Double_t               zPrimVertex;
    Double_t               yPrimVertex;
    Double_t               xPrimVertex;
    Double_t               sigmaVertex;
    Double_t               dPointToSecVertexLine;
    Int_t                  nTracks;
    Bool_t                 fisMC;
//    Int_t                  stackNTracksNoCut;
//    Int_t                  stackNTracks;
    Int_t                  nPhotonCandidates;

    AliVEvent*              fEvent;
    const AliVVertex*       fPrimVertex;
    AliESDEvent*            fESDEvent;
    AliAODEvent*            fAODEvent;
    AliESDtrackCuts*        fTrackCuts;
    //TString                 fV0ReaderName;
    AliV0ReaderV1*          fV0Reader;
    AliMCEvent*             fMCEvent;
    AliStack*               fMCStack;
    AliPIDResponse*         fPIDResponse; //!PID Response
    AliVTrack**             trueTrack;
    AliMyTParticleMotherAndDaughter* myParticle;
    AliConversionPhotonCuts* fConversionCuts;


    TH1D*        fHistZVertex;
    TH2D*        fHistXY;
    TH1F*        fHistPtNegativeParticles;
    TH1F*        fHistPtPositiveParticles;
    TH1D*        fHistTrueDCA;
//    TH1D*        fHistTrueDCAToPrimVtx;
    TH2D*        fHistInvMassPosElec;
 //   TH1F*        fHistInvMassPosElecAfterDCACut;
    TH1F*        fHistInvMassPosElecAfterDCut;
    TH1F*        fHistInvMassPosElecBack;
    TH1D*        fHistDistanceToPrimVertex;
    TH1F*        fHistTrueInvMassPosElecNoCut;
    TH1F*        fHistTrueInvMassPosElecDCACut;
    TH1F*        fHistTrueInvMassPosElecDCAvertexCut;
    TH1F*        fHistNSigmasTPCElectron1;
    TH1F*        fHistNSigmasTPCElectron2;
    TH1F*        fHistNSigmasTPCElectron3;
    TH1F*        fHistNSigmasTPCElectron4;
    TH1F*        fHistNSigmasTPCElectron5;
    TH1F*        fHistNSigmasTPCElectron6;
    TH2D*        fHistInvMassNeutralMeson;
    TH1D*        fHistTrueMesonPt;
    TH1D*        fHistInvMassBackNeutralMeson;
    TH2F*        fHistNSigmasMomentumTPC;

    AliAnalysisTaskDirectPhotons(const AliAnalysisTaskDirectPhotons&); // not implemented
    AliAnalysisTaskDirectPhotons& operator=(const AliAnalysisTaskDirectPhotons&); // not implemented


ClassDef(AliAnalysisTaskDirectPhotons, 1);
};

#endif

#ifndef AliMyTParticle_h
#define AliMyTParticle_h


class AliMyTParticleMotherAndDaughter : public TObject
{
    public: AliMyTParticleMotherAndDaughter() : 

    TObject(), labelDaughterElectron(0), labelDaughterPositron(0), labelMotherPhoton(0) {;}

    public:

    //TParticle* daughterElectron;
    Int_t labelDaughterElectron;

    //TParticle* daughterPositron;
    Int_t labelDaughterPositron;
	
    //TParticle* motherPhoton;
    Int_t labelMotherPhoton;

    ClassDef(AliMyTParticleMotherAndDaughter,1); //primeira versão, mudar a cada atualização

};

#endif
