#ifndef AliAodSkimTask_H
#define AliAodSkimTask_H

/// \class AliAodSkimTask
/// \brief Use to skim AOD files
///
/// Class to skim AOD files with the idea to keep the skimmed file as close as possible to the original AOD.
///
/// \author C.Loizides

#include <AliAnalysisTaskSE.h>
class AliAODMCHeader;
class TH1F;

class AliAodSkimTask: public AliAnalysisTaskSE  
{
  public:
    AliAodSkimTask();
    AliAodSkimTask(const char *name);
    virtual              ~AliAodSkimTask();
    void                  SetClusMinE(Double_t v)   {fClusMinE=v;}
    void                  SetCutMC(Bool_t b)        {fCutMC=b;}
    void                  SetYCutMC(Double_t v)     {fYCutMC=v;}
    void                  SetCopyHeader(Bool_t b)   {fDoCopyHeader=b;}
    void                  SetCopyVZERO(Bool_t b)    {fDoCopyVZERO=b;}
    void                  SetCopyTZERO(Bool_t b)    {fDoCopyTZERO=b;}
    void                  SetCopyVertices(Bool_t b) {fDoCopyVertices=b;}
    void                  SetCopyTOF(Bool_t b)      {fDoCopyTOF=b;}
    void                  SetCopyTracklets(Bool_t b){fDoCopyTracklets=b;}
    void                  SetCopyTracks(Bool_t b)   {fDoCopyTracks=b;}
    void                  SetCopyTrigger(Bool_t b)  {fDoCopyTrigger=b;}
    void                  SetCopyPTrigger(Bool_t b) {fDoCopyPTrigger=b;}
    void                  SetCopyCells(Bool_t b)    {fDoCopyCells=b;}
    void                  SetCopyPCells(Bool_t b)   {fDoCopyPCells=b;}
    void                  SetCopyClusters(Bool_t b) {fDoCopyClusters=b;}
    void                  SetCopyDiMuons(Bool_t b)  {fDoCopyDiMuons=b;}
    void                  SetCopyZDC(Bool_t b)      {fDoCopyZDC=b;}
    void                  SetCopyMC(Bool_t b)       {fDoCopyMC=b;}
    void                  SetCopyMCHeader(Bool_t b) {fDoCopyMCHeader=b;}
    const char           *Str() const;
  protected:
    void                  UserCreateOutputObjects();
    void                  UserExec(Option_t* option);
    Bool_t                UserNotify();
    void                  Terminate(Option_t* option);
    Bool_t                PythiaInfoFromFile(const char *currFile, Float_t &xsec, Float_t &trials, Int_t &pthard);
    Double_t              fClusMinE;        //  minimum cluster energy to accept event
    Bool_t                fCutMC;           //  if true cut MC particles with |Y|>fYCutMC
    Double_t              fYCutMC;          //  cut for MC particles (default = 0.7)
    Bool_t                fDoCopyHeader;    //  if true copy header
    Bool_t                fDoCopyVZERO;     //  if true copy VZERO
    Bool_t                fDoCopyTZERO;     //  if true copy TZERO
    Bool_t                fDoCopyVertices;  //  if true copy vertices
    Bool_t                fDoCopyTOF;       //  if true copy TOF
    Bool_t                fDoCopyTracklets; //  if true copy tracklets
    Bool_t                fDoCopyTracks;    //  if true copy tracks
    Bool_t                fDoCopyTrigger;   //  if true copy trigger (EMC)
    Bool_t                fDoCopyPTrigger;  //  if true copy trigger (PHS)
    Bool_t                fDoCopyCells;     //  if true copy cells (EMC)
    Bool_t                fDoCopyPCells;    //  if true copy cells (PHS)
    Bool_t                fDoCopyClusters;  //  if true copy clusters
    Bool_t                fDoCopyDiMuons;   //  if true copy dimuons
    Bool_t                fDoCopyZDC;       //  if true copy zdc
    Bool_t                fDoCopyMC;        //  if true copy MC particles
    Bool_t                fDoCopyMCHeader;  //  if true copy MC header
    UInt_t                fTrials;          //! events seen since last acceptance 
    Float_t               fPyxsec;          //! pythia xsection
    Float_t               fPytrials;        //! pythia trials
    Int_t                 fPypthardbin;     //! pythia pthard bin
    AliAODEvent          *fAOD;             //! input event
    AliAODMCHeader       *fAODMcHeader;     //! MC header
    TList                *fOutputList;      //! output list
    TH1F                 *fHevs;            //! events processed/accepted
    TH1F                 *fHclus;           //! cluster distribution
    const char           *GetVersion() const { return "1.1"; }

    AliAodSkimTask(const AliAodSkimTask&);             // not implemented
    AliAodSkimTask& operator=(const AliAodSkimTask&);  // not implemented
    ClassDef(AliAodSkimTask, 3); // AliAodSkimTask
};
#endif
