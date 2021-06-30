#ifndef AliAnalysisTaskV0multspec_H
#define AliAnalysisTaskV0multspec_H

//class AliESDpid;
//class AliESDEvent;

#include "TString.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"

enum PartType {kk0s, klam, kalam, kxim, kxip, komm, komp, kallparts};

struct V0filler{

    unsigned char Mult;
    Double32_t Pt;
    Double32_t InvMass;
    Double32_t DcaV0Daught;     //[0,2.54,8]
    Double32_t DcaPosToPV;      //[0,2.54,8]
    Double32_t DcaNegToPV;      //[0,2.54,8]
    Double32_t V0Rad;           //[0,25.4,8]
    Double32_t V0CosPA;         //[0.95,1,16]
    Double32_t LeastCRawsOvF;   //[0,2.54,8]
    Double32_t DistOverTotP;    //[0,254,8]
    Double32_t NSigPos;         //[-10,10,8]
    Double32_t NSigNeg;         //[-10,10,8]
    int    LeastCRaws;          //[0,254,8]
    bool   TOFmatch;
    bool   ITSmatch;
    bool   IsNewEvt;
    int    nGen;                //[0,254,8]
    
};

struct Cascfiller{
    
    uint8_t Mult;
    Double32_t Pt;
    Double32_t InvMassXi;
    Double32_t InvMassOm;
    Double32_t CascCosPA;          //[0.95,1,16]  
    Double32_t V0CosPA;            //[0.95,1,16]    
    Double32_t CascRad;            //[0,25.4,8]  
    Double32_t V0Rad;              //[0,25.4,8]
    Double32_t NSigPosProton;      //[-5,5,8]
    Double32_t NSigPosPion;        //[-5,5,8]
    Double32_t NSigNegProton;      //[-5,5,8]
    Double32_t NSigNegPion;        //[-5,5,8]
    Double32_t NSigBacPion;        //[-5,5,8]
    Double32_t NSigBacKaon;        //[-5,5,8]
    Double32_t LeastCRawsOvF;
    Double32_t InvMassLam;   
    Double32_t DcaCascDaught;      //[0,2.54,8] 
    Double32_t DcaV0Daught;        //[0,2.54,8] 
    Double32_t DcaV0ToPV;          //[0,2.54,8]  
    Double32_t DcaBachToPV;        //[0,2.54,8] 
    Double32_t DcaPosToPV;         //[0,2.54,8]  
    Double32_t DcaNegToPV;         //[0,2.54,8]
    Double32_t DistOverTotP;
    Double32_t BacBarCosPA;
    int    LeastCRaws;
    bool   IsNewEvt;
    bool   charge;
    bool   TOFmatch;
    bool   ITSmatch;
    bool   yXi;
    bool   yOm;
    
};


class AliAnalysisTaskV0multspec : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskV0multspec();
    AliAnalysisTaskV0multspec(const char *name, int particle, bool AddCasc, TString lExtraOptions = "", bool ismc = kTRUE);
    virtual ~AliAnalysisTaskV0multspec();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
    
    void SetAddCasc(bool AddCasc) {fAddCasc=AddCasc;}
    void SetIsMC(bool mc) {fisMC=mc;}
    void SetParticle(int part) {if(part<kallparts) fpart = part; else{printf("Invalid particle. Aborting..."); return;}}

  private:
    
    //outputs
    THistManager *fHistos_misc;                               //!<! Output histos
    TTree        *fTreeV0;                                    //!<! Output Tree, V0s
    TTree        *fTreeCascade;                               //!<! Output Tree, Cascades

    //fillers
    V0filler *ffillV0 = nullptr;                               //!<! Transient V0filler
    Cascfiller *ffillCasc = nullptr;                           //!<! Transient Cascfiller
    
    //objects retreived from input handler
    AliPIDResponse *fPIDResponse;                              //!

    //decide if to store cascade tree
    bool fAddCasc;                                             //!
    int fpart;                                                 //! which particle are we analyzing?

    //flag for MC handeling
    bool fisMC;                                                 //

    //functions to allow flushing part of code out of UserExec
    void DataPosting();
    
    AliAnalysisTaskV0multspec(const AliAnalysisTaskV0multspec&);            // not implemented
    AliAnalysisTaskV0multspec& operator=(const AliAnalysisTaskV0multspec&); // not implemented

    ClassDef(AliAnalysisTaskV0multspec, 1); 
    //version 1: first beta version
};

#endif
