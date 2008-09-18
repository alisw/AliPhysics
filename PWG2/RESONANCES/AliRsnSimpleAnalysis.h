/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnSimpleAnalysis
//             Reconstruction and analysis of K* Rsn
// ........................................
// ........................................
// ........................................
// ........................................
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef AliRsnSimpleAnalysis_H
#define AliRsnSimpleAnalysis_H

class TTree;
class AliRsnPID;
class AliRsnSimpleAnalyzer;

class AliRsnSimpleAnalysis : public TObject
{
  public:

    AliRsnSimpleAnalysis(AliRsnSimpleAnalyzer *ana = 0x0, AliRsnPID *pid = 0x0);
    virtual ~AliRsnSimpleAnalysis() {Clear();}
    virtual void Clear(Option_t *option = "");

    void    SetPID(AliRsnPID *pid) {fPID = pid;}
    void    SetAnalyzer(AliRsnSimpleAnalyzer *analyzer) {fAnalyzer = analyzer;}
    void    SetEventsTree(TTree *tree);
    void    SetFileName(char *fname) {strcpy(fFileName, fname);}
    void    SetStep(Int_t step) {fStep = step;}

    Bool_t  Initialize();
    Stat_t  Process();
    void    SaveOutput() const;

  private:

    AliRsnSimpleAnalysis(const AliRsnSimpleAnalysis &copy) :
        TObject(copy),fInitialized(kFALSE),fStep(1000),fTree(0x0),fPID(0x0),fAnalyzer(0x0) { }
    AliRsnSimpleAnalysis& operator=(const AliRsnSimpleAnalysis & /*copy*/)
    { return (*this); }

    Bool_t                fInitialized;     // flag to check initialization
    Int_t                 fStep;            // progress step
    Char_t                fFileName[250];   // output file name
    TTree                *fTree;            //! TTree of events
    AliRsnPID            *fPID;             //! PID manager
    AliRsnSimpleAnalyzer *fAnalyzer;        //! analyzer

    ClassDef(AliRsnSimpleAnalysis,1)        // dictionary
};

#endif
