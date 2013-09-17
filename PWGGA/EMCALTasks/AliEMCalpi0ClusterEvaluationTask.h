#ifndef AliEMCalpi0ClusterEvaluationTask_h
#define AliEMCalpi0ClusterEvaluationTask_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \class AliEMCalpi0ClusterEvaluationTask
/// \Task to analyse ESDs for cluster studies
///
//  \Astrid Morreale -

#include "AliAnalysisTaskSE.h"
class TTree;
class TString;
class TList;
class TH1F;
class TH2F;
class AliESDEvent;



class AliEMCalpi0ClusterEvaluationTask : public AliAnalysisTaskSE
{
    public:
    AliEMCalpi0ClusterEvaluationTask(const char *name = "AliEMCalpi0ClusterEvaluationTask");

    virtual ~AliEMCalpi0ClusterEvaluationTask();

    // virtual void GetMom(TLorentzVector& p, const AliVCluster *c1, const AliVCluster *c2, Double_t *vertex);
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(const Option_t*);
    virtual void FillMixed(const TLorentzVector& p1, const TLorentzVector& p2);
    virtual Double_t GetMaxCellEnergy(const AliVCluster *cluster, Int_t &id) const;

    private:
    void InitHistPointers();
    /// copy constructor (not implemented )
    AliEMCalpi0ClusterEvaluationTask( const AliEMCalpi0ClusterEvaluationTask& );

    /// assignment operator (not implemented )
    AliEMCalpi0ClusterEvaluationTask& operator = ( const AliEMCalpi0ClusterEvaluationTask& );

    enum {kNtype = 10};
    TH2F * fMasspi0EGA[kNtype];
    TH2F * fMasspi0MB[kNtype];
    TH2F * fMasspi0AllMB[kNtype];
    TH2F * fMasspi0Cent[kNtype];
    TH2F * fMasspi0SemiCent[kNtype];

    TH2F * fMassMixedEGA[kNtype];
    TH2F * fMassMixedMB[kNtype];
    TH2F * fMassMixedAllMB[kNtype];
    TH2F * fMassMixedCent[kNtype];
    TH2F * fMassMixedSemiCent[kNtype];


    TH1F * fEventsEGA[kNtype];
    TH1F * fEventsMB[kNtype];
    TH1F * fEventsAllMB[kNtype];
    TH1F * fEventsCent[kNtype];
    TH1F * fEventsSemiCent[kNtype];


    TH1F * fCentrality;
    TH1F * fCentralityMB;
    TH1F * fCentralityEGA;
    TH1F * fCentralityCent;
    TH1F * fCentralitySemiCent;

    TH1F * fTriggers;

    TH1F * fpTMB[kNtype];
    TH1F * fpTAllMB[kNtype];
    TH1F * fpTEGA[kNtype];
    TH1F * fpTkCent[kNtype];
    TH1F * fpTkSemiCent[kNtype];
    TH1F * fDispersion;
    TH1F * fexo;
    TH1F * fshower;

    /// local event counter
    Int_t  fEvent;
    Int_t   ega0,   ega1,  ega2,  ega3,  ega4,  ega5,  ega6,   ega7,  ega8, ega9;
    Int_t    mb0,    mb1,   mb2,   mb3,   mb4,   mb5,   mb6,    mb7,   mb8, mb9;
    Int_t allmb0, allmb1,allmb2,allmb3,allmb4,allmb5,allmb6, allmb7,allmb8, allmb9;
    Int_t cent0,   cent1, cent2, cent3, cent4, cent5, cent6,  cent7, cent8, cent9;
    Int_t semicent0,   semicent1, semicent2, semicent3, semicent4, semicent5, semicent6,  semicent7, semicent8, semicent9;
    Int_t all, allmb, mb, central, semicentral, ega;



    Bool_t  kAllMB;
    Bool_t  isPileup;
    Bool_t  isMB;
    Bool_t  isAnyINT;
    Bool_t  isCentral;
    Bool_t  isSemiCentral;
    Bool_t  isEga;

    Bool_t  isMBmx;
    Bool_t  isAnyINTmx;
    Bool_t  isCentralmx;
    Bool_t  isSemiCentralmx;
    Bool_t  isEgamx;
    Bool_t  kAllMBmx;

    Int_t   trigger;
    Float_t CentralityVZERO;
    Float_t CentralitySPD;
    Int_t   runNumber;
    Int_t   selectionMask;
    Float_t vX; Float_t vY; Float_t vZ;




    //characteristiques cluster
    Float_t   Ecluster;
    Int_t     NCellscluster;
    Float_t   M20cluster;
    Float_t   M02cluster;
    Int_t     NCluscluster;
    Bool_t    isEMCALcluster;
    Float_t   dispersioncluster;
    Float_t   chi2cluster;
    Double_t     distBadChannelcluster;
    Float_t   phicluster;
    Float_t   etacluster;
    Float_t   ptcluster;
    Double_t crossEnergy;
    //characteristics pion
    Float_t   piE;
    Float_t   piphi;
    Float_t   pieta;
    Float_t   ptpi;
    Float_t   pipx;
    Float_t   pipy;
    Float_t   pipz;
    Float_t   asympi;
    Float_t   masspi;

    TList *fHistList;
    TObjArray  *fPool[kNtype];

    ClassDef(AliEMCalpi0ClusterEvaluationTask,1)

};

#endif

