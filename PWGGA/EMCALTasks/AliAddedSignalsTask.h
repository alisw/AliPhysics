#ifndef AliAddedSignalsTask_h
#define AliAddedSignalsTask_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */
// $Id$
/// \class AliAddedSignalsTask
/// \Task to analyse AODs for pion/eta studies
///
//  \Astrid Morreale -subatech

#include "AliAnalysisTaskSE.h"
#include "Alipi0EventStatStruct.h"
class TTree;
class TString;
class TList;
class AliESDEvent;



class AliAddedSignalsTask : public AliAnalysisTaskSE
{
    public:
    AliAddedSignalsTask(const char *name = "AliAddedSignalsTask");

    virtual ~AliAddedSignalsTask();

    // virtual void GetMom(TLorentzVector& p, const AliVCluster *c1, const AliVCluster *c2, Double_t *vertex);
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(const Option_t*);


    virtual Double_t GetMaxCellEnergy(const AliVCluster *cluster, Int_t &id) const;
  //  virtual void TsallisWeight(Double_t  &ptMC);

    private:

    /// copy constructor (not implemented )
    AliAddedSignalsTask( const AliAddedSignalsTask& );

    /// assignment operator (not implemented )
    AliAddedSignalsTask& operator = ( const AliAddedSignalsTask& );



    /// local event counter
    Int_t fEvent;

    //Branch with event stats
    EventStatStruct* fEventStatStruct;

    /// list of track records
    TClonesArray* fClusterStatArray;
    TClonesArray* fClusterMCStatArray;
    TClonesArray* fClusterpiMCStatArray;
    TClonesArray* fClusterHijingMCStatArray;

     /// number of entries in TClonesArray
    Int_t fClusterStatCount;
    Int_t fClusterMCStatCount;
    Int_t fClusterpiMCStatCount;
    Int_t fClusterHijingMCStatCount;



    ClassDef(AliAddedSignalsTask,1)

};

#endif

