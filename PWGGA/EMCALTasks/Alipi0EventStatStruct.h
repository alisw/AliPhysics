#ifndef ALIPI0EVENTSTATSTRUCT_H
#define ALIPI0EVENTSTATSTRUCT_H

#include <TObject.h>
/// \class Alipi0EventStatStruct
/// \brief Event information structure pion/eta in PbPb
///
/// \author Astrid Morreale astridmorreale@cern.ch subatech
/// \date April 10 2015

class EventStatStruct: public TObject
{
    public:

    /// object name (re-implemented)
    virtual const char*	GetName() const
    { return "eventStatStruct"; }

    /// constructor
    EventStatStruct( void):

        kAllMB(kFALSE),
        isPileup(kFALSE),
         isMB(kFALSE),
        isAnyINT(kFALSE),
        isCentral(kFALSE),
        isSemiCentral(kFALSE),
        isEga(kFALSE),
        CentralityVZERO(0),
        CentralitySPD(0),
        multiplicity(0),
        runNumber(0),
        nV(0),
        vX(0),
        vY(0),
        vZ(0),
        EMCalClusters(0),
        EMCalHijingClusters(0),
        NTotalClusters(0)
    {}

    Bool_t  kAllMB;
    Bool_t  isPileup;
    Bool_t  isMB;
    Bool_t  isAnyINT;
    Bool_t  isCentral;
    Bool_t  isSemiCentral;
    Bool_t  isEga;

    Float_t CentralityVZERO;
    Float_t CentralitySPD;
    Int_t   multiplicity;
    Int_t   runNumber;
    Int_t   nV;
    Float_t vX;
    Float_t vY;
    Float_t vZ;

    Int_t   EMCalClusters;
    Int_t   EMCalHijingClusters;
    Int_t   NTotalClusters;

    ClassDef(EventStatStruct,1)

};

#endif
