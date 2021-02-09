#ifndef AliEMCalpi0AddedSignalsTask_h
#define AliEMCalpi0AddedSignalsTask_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */
// $Id$
/// \class AliEMCalpi0AddedSignalsTask
/// \Task to analyse AODs for pion/eta studies
///
//  \Astrid Morreale -subatech

#include "AliAnalysisTaskSE.h"
#include "Alipi0EventStatStruct.h"
class TTree;
class TString;
class TList;
class AliESDEvent;



class AliEMCalpi0AddedSignalsTask : public AliAnalysisTaskSE
{
    public:
    AliEMCalpi0AddedSignalsTask(const char *name = "AliEMCalpi0AddedSignalsTask");

    virtual ~AliEMCalpi0AddedSignalsTask();

    // virtual void GetMom(TLorentzVector& p, const AliVCluster *c1, const AliVCluster *c2, Double_t *vertex);
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(const Option_t*);


    virtual Double_t GetMaxCellEnergy(const AliVCluster *cluster, Int_t &id) const;
  //  virtual void TsallisWeight(Double_t  &ptMC);

    private:

    /// copy constructor (not implemented )
    AliEMCalpi0AddedSignalsTask( const AliEMCalpi0AddedSignalsTask& );

    /// assignment operator (not implemented )
    AliEMCalpi0AddedSignalsTask& operator = ( const AliEMCalpi0AddedSignalsTask& );



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

    ClassDef(AliEMCalpi0AddedSignalsTask,1)

};

#endif

