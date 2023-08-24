#ifndef AliAnalysisTaskV0multspec_H
#define AliAnalysisTaskV0multspec_H

//class AliESDpid;
//class AliESDEvent;

#include "TString.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"
#include "AliMCEvent.h"

enum PartType {kk0s, klam, kalam, kxi, kom, kallparts};
const int PDGcodes[5] = {310, 3122, -3122, 3312, 3334};

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
    //MC-related
    int    nGen;                //[0,254,8]
    bool   IsPDGmatched;
    
};

struct Cascfiller{
    
    unsigned char Mult;
    Double32_t Pt;
    Double32_t InvMass;
    Double32_t CascCosPA;          //[0.95,1,16]  
    Double32_t V0CosPA;            //[0.95,1,16]    
    Double32_t CascRad;            //[0,25.4,8]  
    Double32_t V0Rad;              //[0,25.4,8]
    Double32_t NSigPos;            //[-10,10,8]
    Double32_t NSigNeg;            //[-10,10,8]
    Double32_t NSigBac;            //[-10,10,8]
    Double32_t LeastCRawsOvF;      //[0,2.54,8]
    Double32_t DcaCascDaught;      //[0,2.54,8] 
    Double32_t DcaV0Daught;        //[0,2.54,8] 
    Double32_t DcaV0ToPV;          //[0,2.54,8]  
    Double32_t DcaBachToPV;        //[0,2.54,8] 
    Double32_t DcaPosToPV;         //[0,2.54,8]  
    Double32_t DcaNegToPV;         //[0,2.54,8]
    Double32_t DistOverTotP;       //[0,254,8]
    Double32_t BacBarCosPA;        //[0.999,1,8]  
    int    LeastCRaws;             //[0,254,8]
    bool   IsNewEvt;
    bool   charge; //kTRUE==negative, kFALSE==positive
    bool   TOFmatch;
    bool   ITSmatch;
    bool   GoodInvMassLam; //kTRUE if ( PDGmass-0.008 GeV/c2 <ImassLam< PDGmass+0.008 GeV/c2)
    //MC-related
    int    nGenP;                    //[0,254,8]
    int    nGenM;                    //[0,254,8]
    bool   IsPDGmatched;
};


class AliAnalysisTaskV0multspec : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskV0multspec();
    AliAnalysisTaskV0multspec(const char *name, int particle, TString lExtraOptions = "", bool ismc = kTRUE);
    virtual ~AliAnalysisTaskV0multspec();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
    
    void SetIsMC(bool mc) {fisMC=mc;}
    void SetParticle(int part) {if(part<kallparts) fpart = part; else{printf("Invalid particle. Aborting..."); return;}}

  private:
    
    //outputs
    THistManager *fHistos_misc;                               //!<! Output histos
    TTree        *fTree;                                      //!<! Output Tree

    //fillers
    V0filler *ffillV0 = nullptr;                               //!<! Transient V0filler
    Cascfiller *ffillCasc = nullptr;                           //!<! Transient Cascfiller
    
    //objects retreived from input handler
    AliPIDResponse *fPIDResponse;                              //!

    //decide if to store cascade tree
    int fpart;                                                 // which particle are we analyzing?

    //flag for MC handeling
    bool fisMC;                                                 //

    //functions to allow flushing part of code out of UserExec
    void DataPosting();
    
    AliAnalysisTaskV0multspec(const AliAnalysisTaskV0multspec&);            // not implemented
    AliAnalysisTaskV0multspec& operator=(const AliAnalysisTaskV0multspec&); // not implemented

    ClassDef(AliAnalysisTaskV0multspec, 2); 
    //version 1: first beta version
    //version 2: updated casc part and correct handeling of "zeros"
};

#endif
