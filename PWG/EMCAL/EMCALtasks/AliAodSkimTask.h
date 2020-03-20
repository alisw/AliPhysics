#ifndef AliAodSkimTask_H
#define AliAodSkimTask_H

/// \class AliAodSkimTask
/// \brief Use to skim AOD files
///
/// Class to skim AOD files with the idea to keep the skimmed file as close as possible to the original AOD.
///
/// \author C.Loizides

#include <AliAnalysisTaskSE.h>
#include <TString.h>
class AliAODMCHeader;
class TH1F;

class AliAodSkimTask: public AliAnalysisTaskSE
{
  public:
    AliAodSkimTask(const char *name=0);
    virtual              ~AliAodSkimTask();
    void                  SetCleanTracklets(Bool_t b)         {fDoCleanTracklets=b;}
    void                  SetCleanTracks(Bool_t b)            {fDoCleanTracks=b;}
    void                  SetClusMinE(Double_t v)             {fClusMinE=v;}
    void                  SetTrackMinPt(Double_t v)           {fTrackMinPt=v;}
    void                  SetTrackMaxPt(Double_t v)           {fTrackMaxPt=v;}
    void                  SetDoBothMinTrackAndClus(Bool_t b)  {fDoBothMinTrackAndClus=b;}
    void                  SetCopyCascades(Bool_t b)           {fDoCopyCascades=b;}
    void                  SetCopyCells(Bool_t b)              {fDoCopyCells=b;}
    void                  SetCopyClusters(Bool_t b)           {fDoCopyClusters=b;}
    void                  SetCopyConv(Bool_t b)               {fDoCopyConv=b;}
    void                  SetCopyDiMuons(Bool_t b)            {fDoCopyDiMuons=b;}
    void                  SetCopyHeader(Bool_t b)             {fDoCopyHeader=b;}
    void                  SetCopyMC(Bool_t b)                 {fDoCopyMC=b;}
    void                  SetCopyMCHeader(Bool_t b)           {fDoCopyMCHeader=b;}
    void                  SetCopyKinks(Bool_t b)              {fDoCopyKinks=b;}
    void                  SetCopyPCells(Bool_t b)             {fDoCopyPCells=b;}
    void                  SetCopyPTrigger(Bool_t b)           {fDoCopyPTrigger=b;}
    void                  SetCopyTOF(Bool_t b)                {fDoCopyTOF=b;}
    void                  SetCopyTZERO(Bool_t b)              {fDoCopyTZERO=b;}
    void                  SetCopyTracklets(Bool_t b)          {fDoCopyTracklets=b;}
    void                  SetCopyTracks(Bool_t b)             {fDoCopyTracks=b;}
    void                  SetCopyTrdTracks(Bool_t b)          {fDoCopyTrdTracks=b;}
    void                  SetCopyTrigger(Bool_t b)            {fDoCopyTrigger=b;}
    void                  SetCopyV0s(Bool_t b)                {fDoCopyV0s=b;}
    void                  SetCopyUserTree(Bool_t b)           {fDoCopyUserTree=b;}
    void                  SetCopyVZERO(Bool_t b)              {fDoCopyVZERO=b;}
    void                  SetCopyVertices(Bool_t b)           {fDoCopyVertices=b;}
    void                  SetCopyZDC(Bool_t b)                {fDoCopyZDC=b;}
    void                  SetCutFilterBit(UInt_t b)           {fCutFilterBit=b;}
    void                  SetCutMcIsPrimary(Bool_t b)         {fCutMcIsPrimary=b;}
    void                  SetCutMcIsPhysicalPrimary(Bool_t b) {fCutMcIsPhysPrimary=b;}
    void                  SetCutMcPt(Double_t pt)             {fCutMcPt=pt;}
    void                  SetCutMcY(Double_t y)               {fCutMcY=y;}
    void                  SetCutMcPhos(Bool_t b)              {fCutMcPhos=b;}
    void                  SetCutMcEmcal(Bool_t b)             {fCutMcEmcal=b;}
    void                  SetDoVertMain(Bool_t b)             {fDoVertMain=b;}
    void                  SetDoVertWoRefs(Bool_t b)           {fDoVertWoRefs=b;}
    void                  SetDoPhosFilter(Bool_t b)           {fDoPhosFilt=b;}
    void                  SetGammaBrName(TString s)           {fGammaBr=s;}
    void                  SetMinCutPt(Double_t pt)            {fCutMinPt=pt;}
    void                  SetRemCovMat(Bool_t b)              {fDoRemCovMat=b;}
    void                  SetRemPid(Bool_t b)                 {fDoRemPid=b;}
    void                  SetRemoveMcParts(Bool_t b)          {fDoRemoveMcParts=b;}
    void                  SetRemoveTracks(Bool_t b)           {fDoRemoveTracks=b;}
    const char           *Str() const;
  protected:
    virtual void          CleanTrack(AliAODTrack *t);
    const char           *GetVersion() const { return "1.6"; }
    virtual Bool_t        KeepMcPart(AliAODMCParticle *p);
    Bool_t                IsDcalAcc(Double_t phi, Double_t eta);
    Bool_t                IsPhosAcc(Double_t phi, Double_t eta);
    Bool_t                IsEmcalAcc(Double_t phi, Double_t eta);
    virtual Bool_t        KeepTrack(AliAODTrack *t);
    Bool_t                PythiaInfoFromFile(const char *currFile, Float_t &xsec, Float_t &trials, Int_t &pthard);
    virtual Bool_t        SelectEvent();
    void                  Terminate(Option_t* option);
    void                  UserCreateOutputObjects();
    void                  UserExec(Option_t* option);
    Bool_t                UserNotify();

    virtual void          CopyCascades();
    virtual void          CopyCells();
    virtual void          CopyCellsP();
    virtual void          CopyConv();
    virtual void          CopyClusters();
    virtual void          CopyDimuons();
    virtual void          CopyHeader();
    virtual void          CopyKinks();
    virtual void          CopyMc();
    virtual void          CopyMcHeader();
    virtual void          CopyTof();
    virtual void          CopyTracklets();
    virtual void          CopyTracks();
    virtual void          CopyTrdTracks();
    virtual void          CopyTrigger();
    virtual void          CopyTriggerP();
    virtual void          CopyTZero();
    virtual void          CopyUserTree();
    virtual void          CopyVertices();
    virtual void          CopyV0s();
    virtual void          CopyVZero();
    virtual void          CopyZdc();
    virtual void          CopyMore() {;}
  
    Double_t              fClusMinE;              //  minimum cluster energy to accept event
    Double_t              fTrackMinPt;            //  minimum track pt to accept event
    Double_t              fTrackMaxPt;            //  maximum track pt to accept event
    Bool_t                fDoBothMinTrackAndClus; // switch to enable simultaneous filtering for minimum track and cluster cuts
    Double_t              fCutMinPt;              //  minimum pT to keep track
    UInt_t                fCutFilterBit;          //  filter bit(s) to select
    TString               fGammaBr;               //  gamma branch name
    Bool_t                fDoCopyHeader;          //  if true copy header
    Bool_t                fDoCopyVZERO;           //  if true copy VZERO
    Bool_t                fDoCopyTZERO;           //  if true copy TZERO
    Bool_t                fDoCopyVertices;        //  if true copy vertices
    Bool_t                fDoCopyTOF;             //  if true copy TOF
    Bool_t                fDoCopyTracklets;       //  if true copy tracklets
    Bool_t                fDoCopyTracks;          //  if true copy tracks
    Bool_t                fDoRemoveTracks;        //  if true remove tracks
    Bool_t                fDoCleanTracks;         //  if true clean tracks
    Bool_t                fDoRemCovMat;           //  if true remove cov matrix from tracks
    Bool_t                fDoRemPid;              //  if true remove PID object from tracks
    Bool_t                fDoCopyTrigger;         //  if true copy trigger (EMC)
    Bool_t                fDoCopyPTrigger;        //  if true copy trigger (PHS)
    Bool_t                fDoCopyCells;           //  if true copy cells (EMC)
    Bool_t                fDoCopyPCells;          //  if true copy cells (PHS)
    Bool_t                fDoCopyClusters;        //  if true copy clusters
    Bool_t                fDoCopyDiMuons;         //  if true copy dimuons
    Bool_t                fDoCopyTrdTracks;       //  if true copy trd tracks
    Bool_t                fDoCopyV0s;             //  if true copy v0s
    Bool_t                fDoCopyCascades;        //  if true copy cascades
    Bool_t                fDoCopyZDC;             //  if true copy zdc
    Bool_t                fDoCopyConv;            //  if true copy conversions
    Bool_t                fDoCopyKinks;           //  if true copy kinks
    Bool_t                fDoCopyMC;              //  if true copy MC particles
    Bool_t                fDoCopyMCHeader;        //  if true copy MC header
    Bool_t                fDoVertWoRefs;          //  if true then do not copy TRefs in vertices
    Bool_t                fDoVertMain;            //  if true then only copy main vertices
    Bool_t                fDoCleanTracklets;      //  if true then clean tracklets
    Bool_t                fDoCopyUserTree;        //  if true copy input user tree
    Bool_t                fDoPhosFilt;            //  if true filter on phos (clusters and min energy cut)
    Bool_t                fDoRemoveMcParts;       //  if true remove mc particles
    Bool_t                fCutMcIsPrimary;        //  if true cut on primary mc particles
    Bool_t                fCutMcIsPhysPrimary;    //  if true cut on physical primary mc particles
    Bool_t                fCutMcPt;               //  if true cut on mc particles with pT<fCutMcPt
    Bool_t                fCutMcY;                //  if true cut on mc particles with |y|<fCutMcY
    Bool_t                fCutMcPhos;             //  if true cut particles not in PHOS
    Bool_t                fCutMcEmcal;            //  if true cut particles not in Emcal
    UInt_t                fTrials;                //! events seen since last acceptance
    Float_t               fPyxsec;                //! pythia xsection
    Float_t               fPytrials;              //! pythia trials
    Int_t                 fPypthardbin;           //! pythia pthard bin
    AliAODEvent          *fAOD;                   //! input event
    AliAODMCHeader       *fAODMcHeader;           //! MC header
    TList                *fOutputList;            //! output list
    TH1F                 *fHevs;                  //! events processed/accepted
    TH1F                 *fHclus;                 //! cluster distribution
    TH1F                 *fHtrack;                //! track distribution

    AliAodSkimTask(const AliAodSkimTask&);             // not implemented
    AliAodSkimTask& operator=(const AliAodSkimTask&);  // not implemented
    ClassDef(AliAodSkimTask, 10); // AliAodSkimTask
};
#endif
