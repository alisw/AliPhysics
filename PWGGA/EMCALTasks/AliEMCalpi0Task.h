#ifndef AliEMCalpi0Task_h
#define AliEMCalpi0Task_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */
// $Id$
/// \class AliEMCalpi0Task
/// \brief Task to analyse AODs for pion/eta in PbPb
///
/// \author Astrid Morreale astridmorreale@cern.ch subatech
/// \date April 10 2015
#include "AliAnalysisTaskSE.h"
#include "Alipi0EventStatStruct.h"
class TTree;
class TString;
class TList;
class AliESDEvent;



class AliEMCalpi0Task : public AliAnalysisTaskSE
{
    public:
    AliEMCalpi0Task(const char *name = "AliEMCalpi0Task");

    virtual ~AliEMCalpi0Task();

    // virtual void GetMom(TLorentzVector& p, const AliVCluster *c1, const AliVCluster *c2, Double_t *vertex);
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(const Option_t*);

    virtual void FillMixed(const TLorentzVector& p1, const TLorentzVector& p2);
    virtual Double_t GetMaxCellEnergy(const AliVCluster *cluster, Int_t &id) const;

    private:

    /// copy constructor (not implemented )
    AliEMCalpi0Task( const AliEMCalpi0Task& );

    /// assignment operator (not implemented )
    AliEMCalpi0Task& operator = ( const AliEMCalpi0Task& );

    enum {kNtype = 10};

    /// local event counter
    Int_t  fEvent;



    //Branch with event stats
    EventStatStruct* fEventStatStruct;

    /// list of track records
    TClonesArray* fClusterStatArray;
    TClonesArray* fdiClusterStatArray;
    TClonesArray* fmixedDiClusterStatArray;

     /// number of entries in TClonesArray
    Int_t fClusterStatCount;
    Int_t fdiClusterStatCount;
    Int_t fmixedDiClusterStatCount;


    TObjArray  *fPool[kNtype];

    ClassDef(AliEMCalpi0Task,1)

};

#endif

