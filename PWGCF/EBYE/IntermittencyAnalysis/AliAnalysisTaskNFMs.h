#ifndef AliAnalysisTaskNFMs_H
#define AliAnalysisTaskNFMs_H

class TH1F;
class TH1D;
class TH2D;
class TH2F;
class TList;
class TNtuple;
class TString;
class TList;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAnalysisTaskNFMs : public AliAnalysisTaskSE 
  {
  public:
    AliAnalysisTaskNFMs();
    AliAnalysisTaskNFMs(const char *name);
    virtual ~AliAnalysisTaskNFMs();
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *option);

    //Analysis setters
    void SetFilterBit(UInt_t Bit) {filterBit = Bit;}
    void SetCentLim(int f_minCent, int f_maxCent){minCent = f_minCent; maxCent = f_maxCent;}
    void SetpTbins(float bin1, float bin2, float bin3, float bin4, float bin5, float bin6, float bin7, float bin8)
    {ptbin[0] = bin1; ptbin[1] = bin2; ptbin[2] = bin3; ptbin[3] = bin4; ptbin[4] = bin5; ptbin[5] = bin6; ptbin[6] = bin7; ptbin[7] = bin8;}
    void SetMmax(bool mflag) {fIsMmaxET = mflag;}
    void SetPileUpRead(bool flag) {fIsPileUpCuts = flag;}
    void SetSysflag(bool clflag, bool crflag, float clcut, float crcut) {fClflag=clflag; fCrflag=crflag; if(clflag){fNcls = clcut;} if(crflag){fNcrows = crcut;}}
    void SetDataset(bool dflag){if(dflag) {fiskMB=dflag; fisINT7=kFALSE;} else {fiskMB=kFALSE; fisINT7=kTRUE;} fRun1data=dflag;}

  protected:
    AliEventCuts fEventCuts;
  private:
    AliAODEvent *fAOD;                //AOD object
    TList *fOutHList;                 //Output List
    TList *fQAList; 
    TList *fNtupleListBin1;           // Ntuple to store Fqes
    TList *fNtupleListBin2;       
    TList *fNtupleListBin3;       
    TList *fNtupleListBin4;       
    UInt_t filterBit;             
    double fVxMax, fVyMax, fVzMax, minCent, maxCent, minEta, maxEta, NoOfBins;
    static const int M = 40, Q = 6, D = 2; 
    float ptbin[8], fNcls, fNcrows;
    bool fIsPileUpCuts; //for pileup cut 
    bool fIsMmaxET; //flag for Mmax def
    bool fClflag, fCrflag;
    bool fRun1data;
    bool fiskMB, fisINT7;

    TNtuple *fntpMBin1[M];  //!
    TNtuple *fntpMBin2[M];  //!
    TNtuple *fntpMBin3[M];  //!
    TNtuple *fntpMBin4[M];  //!

    TH1F *fHistPtBin[4];     //! Pt spectrum
    TH1F *fEtaBin[4];        //! Eta spectrum
    TH1F *fPhiBin[4];        //! Phi spectrum
    TH1F *fHistMulBin[4];    //! Histogram to register event multiplicity

    //QA hists
    TH1F *fHistQAVx;//!
    TH1F *fHistQAVy;//!
    TH1F *fHistQAVz;//!
    TH1D *fEventCounter;  //! Histogram to track events etc
    TH1F *fHistQACent;   //!   Histogram for centrality Distribution
    TH1F *fHistQAEta;//!
    TH1F *fHistQAPhi;//!
    TH1F *fHistQAMult;//!
    TH1F *fHistQAMultwpc;//!

    TH2D *fHEtaPhiBin1[M];    //! Eta-Phi for 10 different bins  distribution
    TH2D *fHEtaPhiBin2[M];    //! Eta-Phi for 10 different bins  distribution
    TH2D *fHEtaPhiBin3[M];    //! Eta-Phi for 10 different bins  distribution
    TH2D *fHEtaPhiBin4[M];    //! Eta-Phi for 10 different bins  distribution

    //Hists for Efficiency
    TH2D *fhMapEtaPhiBin1M[M];//!
    TH2D *fhMapEtaPhiBin2M[M];//!
    TH2D *fhMapEtaPhiBin3M[M];//!
    TH2D *fhMapEtaPhiBin4M[M];//!

    AliAnalysisTaskNFMs(const AliAnalysisTaskNFMs &);            // not implemented
    AliAnalysisTaskNFMs &operator=(const AliAnalysisTaskNFMs &); // not implemented

    ClassDef(AliAnalysisTaskNFMs, 1); // example of analysis
  };

#endif
