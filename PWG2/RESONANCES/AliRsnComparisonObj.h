//
// *** Class AliRsnComparisonObj ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef ALIRSNCOMPARISONOBJ_H
#define ALIRSNCOMPARISONOBJ_H

#include <TNamed.h>
#include <TList.h>
#include <TH1.h>

#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"

#include "AliRsnPID.h"

class AliRsnComparisonObj : public TNamed
{
  public:

    enum EFormat
    {
      kRSN=0,
      kESD,
      kMC,
      kLastFormat
    };

    enum EComparisonType
    {
      kParticleInfo=0,
      kPid,
      kLastComparisonType
    };

    enum EParameterType
    {
      kP = 0,
      kPt,
      kEta,
      kY,
      kLastParameterType
    };
    
    enum EHistoType
    {
      kIndent=0,
      kGood,
      kFake,
      kTrue,
      kLastHistoType
    };

    enum EPIDType
    {
      kEsd=0,
      kITS,
      kTPC,
      kTOF,
      kITS_TPC,
      kITS_TOF,
      kTPC_TOF,
      kITS_TPC_TOF,
      kITS_TPC_TOF_SP,
      kLastPIDType
    };

    AliRsnComparisonObj(const char*name="RSN");
    ~AliRsnComparisonObj();

    TList     *GenerateParticleInfoHistogramList(TString prefix="");
    TList     *GeneratePIDHistogramList(TString prefix="");
    void FillPIDHistograms(AliRsnDaughter *daughter);
    void FillPIDHistograms(AliESDtrack *track,AliMCEvent *mc=0);
    void FillPIDHistograms(AliMCParticle *mctrack);
    
    void FillHistograms(AliMCParticle *mctrack);

    void SetCurrentESDPID(const EPIDType& type,const Double_t&divValue = 0.0);
    void SetPriorProbs(Double_t* pid) { for (Int_t i=0; i<5; i++) fPriorProbs[i]=pid[i]; }

    void SetESDstatus(const ULong_t status);
    void SetESDTrackQualityCuts(const Int_t& its=-1,const Int_t& tpc=-1,const Int_t& trd=-1);
  private:
  
    AliRsnComparisonObj(const AliRsnComparisonObj& copy) 
    : TNamed(copy),fCurrentComparisonType(kParticleInfo),fCurrentESDPID(kEsd),
    fESDstatus(0),fITSClusters(0),fTPCClusters(0),fTRDClusters(0),fPIDDivValue(0.) {}
    const AliRsnComparisonObj& operator=(const AliRsnComparisonObj&) {return *this;}

    EComparisonType   fCurrentComparisonType;
    EPIDType          fCurrentESDPID;
    Double_t          fPriorProbs[5];
    ULong_t           fESDstatus;
    Int_t             fITSClusters;
    Int_t             fTPCClusters;
    Int_t             fTRDClusters;
    
    Double_t          fPIDDivValue;

    TH1D              *fHistosPartInfo[kLastParameterType][2][AliRsnPID::kSpeciesAll];
    TH1D              *fHistosPID[kLastFormat][kLastHistoType][AliRsnPID::kSpecies+1];
    
    TString     GetFormatName(EFormat type);
    TString     GetHistoTypeName(EHistoType type);

    void        GetESDPID(AliESDtrack *track,Double_t *pid,Double_t p=-1.0);
    TString     GetParameterName(EParameterType type);
    Double_t    GetParameterNameValue(EParameterType type,AliMCParticle * mctrack);

    ClassDef(AliRsnComparisonObj, 1)
};

#endif
