#ifndef AliAnalysisTaskMultspec_H
#define AliAnalysisTaskMultspec_H

//class AliESDpid;
//class AliESDEvent;

#include "TString.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"
#include "TH1D.h"
#include "TH2D.h"
#include "AliMCEvent.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"

enum PartType_Mult {kk0s_Mult, klam_Mult, kalam_Mult, kxi_Mult, kom_Mult, kallparts_Mult};
const int PDGcodes_Mult[5] = {310, 3122, -3122, 3312, 3334};

struct V0filler_Mult{

    unsigned char Mult;
    Double32_t Pt;
    Double32_t InvMass;
    Double32_t CompMass;
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
    bool   ITSTOFtwo;
    bool   IsNewEvt;
    //MC-related
    int    nGen;                   //[0,254,8]
    bool   IsPDGmatched;
    bool   PhysicalPrimary;
    
};

struct Cascfiller_Mult{
    
    unsigned char Mult;
    Double32_t Pt;
    Double32_t InvMass;
    Double32_t CompMass;
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
    bool   ITSTOFtwo;
    bool   GoodInvMassLam; //kTRUE if ( PDGmass-0.008 GeV/c2 <ImassLam< PDGmass+0.008 GeV/c2)
    //MC-related
    int    nGenP;                    //[0,254,8]
    int    nGenM;                    //[0,254,8]
    bool   IsPDGmatched;
    bool   PhysicalPrimary;
};


class AliAnalysisTaskMultspec : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskMultspec();
    AliAnalysisTaskMultspec(const char *name, int particle, TString lExtraOptions = "", bool ismc = kTRUE, bool removePythiaGen = kTRUE);
    virtual ~AliAnalysisTaskMultspec();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
    
    void SetIsMC(bool mc) {fisMC=mc;}
    void SetParticle(int part) {if(part<kallparts_Mult) fpart = part; else{printf("Invalid particle. Aborting..."); return;}}
    void SetRemovePythiaGen(bool pythiagen) {fRemovePythiaGen=pythiagen;}

  private:
    
    //outputs
    THistManager *fHistos_misc;                               //!<! Output histos
    TTree        *fTree;                                      //!<! Output Tree
    
    //fillers
    V0filler_Mult *ffillV0 = nullptr;                               //!<! Transient V0filler
    Cascfiller_Mult *ffillCasc = nullptr;                           //!<! Transient Cascfiller
    
    //objects retreived from input handler
    AliPIDResponse *fPIDResponse;                              //!
    
    //AliEventCuts object
    AliEventCuts fEventCuts;                                   //

    //pile-up rejection flag
    bool fRejectPileupEvts;                                    //

    //decide if to store cascade tree
    int fpart;                                                 // which particle are we analyzing?

    //flag for MC handeling
    bool fisMC;                                                //

    //flag for remove Pythia events
    bool fRemovePythiaGen;  
    
    //functions to allow flushing part of code out of UserExec
    void DataPosting();
    
    AliAnalysisTaskMultspec(const AliAnalysisTaskMultspec&);            // not implemented
    AliAnalysisTaskMultspec& operator=(const AliAnalysisTaskMultspec&); // not implemented

    ClassDef(AliAnalysisTaskMultspec, 1); 
    //version 1: pt-dependence implemented (in order to extend the analysis as a function of the multiplicity) and add a data-member in order to remove Pythia events 
};

#endif
