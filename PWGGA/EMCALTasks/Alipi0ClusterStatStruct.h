#ifndef ALIPI0CLUSTERSTATSTRUCT_H
#define ALIPI0CLUSTERSTATSTRUCT_H
#include <TObject.h>

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */
/// $Id$
/// \class Alipi0ClusterStatStruct
/// \brief cluster information structure pion/eta in PbPb
///
/// \author Astrid Morreale astridmorreale@cern.ch subatech
/// \date April 10 2015

class ClusterStatStruct: public TObject
{

    public:

    /// object name (re-implemented)
    virtual const char*	GetName() const
    { return "clusterStatStruct"; }

    /// default contructor
    ClusterStatStruct( void ):
        clusterpdg(0),
        clusterMCenergy(0),
        clusterMotherpid(0),
        EMax(0),
        Tcell(0),
        ECross(0),
        ExoC(0),
        RecoMomPDG(0),
        RecoPDG(0),
        fraction(0),
        NCells(0),
        TrackFlag(0),
        TrackFlag1(0),
        TrackFlag2(0),
        Npi(0),
        Neta(0),
        Nhijing(0),
        Nlabels(0),
        distanceBadChannel(0),
        Ecluster(0),
        NCellscluster(0),
        M20cluster(0),
        M02cluster(0),
        NCluscluster(0),
        isEMCALcluster(0),
        dispersioncluster(0),
        chi2cluster(0),
        distBadChannelcluster(0),
        phicluster(0),
        etacluster(0),
        ptcluster(0),
        crossEnergy(0),
        validFlagHijing(kFALSE),
        validFlagpi(kFALSE),
        validFlageta(kFALSE)
    {}

    //characteristics pion
    Int_t     clusterpdg;
    Float_t   clusterMCenergy;
    Int_t     clusterMotherpid;
    Float_t   EMax;
    Float_t   Tcell;
    Float_t   ECross;
    Float_t   ExoC;
    Int_t     RecoMomPDG;
    Int_t     RecoPDG;
    Double_t  fraction;
    Int_t     NCells;
    Bool_t    TrackFlag;
    Bool_t    TrackFlag1;
    Bool_t    TrackFlag2;
    Int_t     Npi;
    Int_t     Neta;
    Int_t     Nhijing;
    Int_t     Nlabels;
    Int_t     distanceBadChannel;
    Float_t   Ecluster;
    Int_t     NCellscluster;
    Float_t   M20cluster;
    Float_t   M02cluster;
    Int_t     NCluscluster;
    Bool_t    isEMCALcluster;
    Float_t   dispersioncluster;
    Float_t   chi2cluster;
    Int_t     distBadChannelcluster;
    Float_t   phicluster;
    Float_t   etacluster;
    Float_t   ptcluster;
    Double_t  crossEnergy;
    Bool_t    validFlagHijing;
    Bool_t    validFlagpi;
    Bool_t    validFlageta;
    //
    ClassDef(ClusterStatStruct,1)

};

#endif
