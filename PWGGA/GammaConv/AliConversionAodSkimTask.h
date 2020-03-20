#ifndef AliConversionAodSkimTask_H
#define AliConversionAodSkimTask_H

/// \class AliConversionAodSkimTask
/// \brief Use to skim AOD files
///
/// Class to skim AOD files with the idea to keep the skimmed file as close as possible to the original AOD.
///
/// \author C.Loizides

#include <AliAodSkimTask.h>
#include <TString.h>
class AliAODMCHeader;
class TH1F;

class AliConversionAodSkimTask: public AliAodSkimTask
{
  public:
    AliConversionAodSkimTask(const char *name=0);
    virtual              ~AliConversionAodSkimTask();
    void                  SetConvMinPt(Double_t v)            {fConvMinPt=v;}
    void                  SetConvMinEta(Double_t eta)         {fConvMinEta=eta;}
    void                  SetConvMaxEta(Double_t eta)         {fConvMaxEta=eta;}
    void                  SetConvMinPhi(Double_t phi)         {fConvMinPhi=phi;}
    void                  SetConvMaxPhi(Double_t phi)         {fConvMaxPhi=phi;}
    void                  SetDoBothConvPtAndAcc(Bool_t b)     {fDoBothConvPtAndAcc=b;}
    void                  SetDoQA(Bool_t b)                   {fDoQA=b;}
  protected:
    void                  UserCreateOutputObjects();
    Bool_t                SelectEvent();
    Double_t              fConvMinPt;             //  minimum conversion photon pT to accept event
    Double_t              fConvMinEta;            //  minimum eta of photon to accept event
    Double_t              fConvMaxEta;            //  maximum eta of photon to accept event
    Double_t              fConvMinPhi;            //  minimum phi of photon to accept event
    Double_t              fConvMaxPhi;            //  maximum phi of photon to accept event
    Bool_t                fDoBothConvPtAndAcc;    //  switch to enable simultaneous filtering for minimum conv pt and conv acceptance cut
    Bool_t                fDoQA;                  //  do QA
    TH1F                 *fHconvPtBeforeCuts;     //! conv photon distribution
    TH1F                 *fHconvPtAfterCuts;      //! conv photon distribution
    TH2F                 *fHconvAccBeforeCuts;    //! conv photon distribution
    TH2F                 *fHconvAccAfterCuts;     //! conv photon distribution

    AliConversionAodSkimTask(const AliConversionAodSkimTask&);             // not implemented
    AliConversionAodSkimTask& operator=(const AliConversionAodSkimTask&);  // not implemented
    ClassDef(AliConversionAodSkimTask, 2); // AliConversionAodSkimTask
};
#endif
