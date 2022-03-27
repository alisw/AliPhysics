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
    virtual void SetFilterBit(UInt_t Bit)
    {
      filterBit = Bit;
    }
    virtual void SetCentLim(int f_minCent, int f_maxCent)
    {
      minCent = f_minCent;
      maxCent = f_maxCent;
    }

private:
AliAODEvent *fAOD;                //AOD object
TList *fOutHList;                 //Output List
TList *fNtupleListBin1;           // Ntuple to store Fqes
TList *fNtupleListBin2;       
TList *fNtupleListBin3;       
TList *fNtupleListBin4;       
UInt_t filterBit;             
double fVxMax, fVyMax, fVzMax, minCent, maxCent, minEta, maxEta;
static const int M = 40, Q = 6, D = 2; 

TNtuple *fntpMBin1[M];  //!
TNtuple *fntpMBin2[M];  //!
TNtuple *fntpMBin3[M];  //!
TNtuple *fntpMBin4[M];  //!

TH1F *fHistPtBin[4];        //! Pt spectrum
TH1F *fEtaBin[4];        //! Eta spectrum
TH1F *fPhiBin[4];        //! Phi spectrum
TH1F *fHistMulBin[4];       //! Histogram to register event multiplicity

//QA hists
TH1F *fHistQAVx;//!
TH1F *fHistQAVy;//!
TH1F *fHistQAVz;//!
TH1D *fEventCounter;  //! Histogram to track events etc
TH1F *fHistQACent;   //!   Histogram for centrality Distribution
TH1F *fHistQAEta;//!
TH1F *fHistQAPhi;//!
                 
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
