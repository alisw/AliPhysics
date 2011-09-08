// -*- mode: C++ -*- 
#ifndef ALIESDEVENT_H
#define ALIESDEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliESDEvent
//   This is the class to deal with during the physics analysis of data.
//   It also ensures the backward compatibility with the old ESD format.
//      
// Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------

#include <TClonesArray.h>
#include <TObject.h>
#include <TTree.h>
#include <TArrayF.h>


#include "AliVEvent.h"
// some includes for delegated methods
#include "AliESDCaloTrigger.h"
#include "AliESDRun.h"
#include "AliESDHeader.h"
#include "AliESDTZERO.h"
#include "AliESDZDC.h"
#include "AliESDACORDE.h"

// AliESDtrack has to be included so that the compiler 
// knows its inheritance tree (= that it is a AliVParticle).
#include "AliESDtrack.h"
// same for AliESDVertex (which is a AliVVertex)
#include "AliESDVertex.h"
// same for CaloCells and CaloClusters (which is a AliVCaloCells, AliVCluster)
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"

#include "AliTOFHeader.h"

#include "AliESDVZERO.h"

class AliESDfriend;
class AliESDHLTtrack;
class AliESDVertex;
class AliESDPmdTrack;
class AliESDFMD;
class AliESDkink;
class AliESDv0;
class AliMultiplicity;
class AliRawDataErrorLog;
class AliESDRun;
class AliESDTrdTrack;
class AliESDTrdTracklet;
class AliESDMuonTrack;
class AliESD;
class AliESDcascade;
class AliCentrality;
class AliEventplane;
class TRefArray;
class AliESDACORDE;
class AliESDHLTDecision;

class TList;
class TString;


class AliESDEvent : public AliVEvent {
public:


  enum ESDListIndex   {kESDRun,
		       kHeader,
		       kESDZDC,
		       kESDFMD,
		       kESDVZERO,
		       kESDTZERO,
		       kTPCVertex,
		       kSPDVertex,
		       kPrimaryVertex,
		       kSPDMult,
		       kPHOSTrigger,
		       kEMCALTrigger,
		       kSPDPileupVertices,
		       kTrkPileupVertices,
		       kTracks,
		       kMuonTracks,
		       kPmdTracks,
		       kTrdTracks,
		       kTrdTracklets,
		       kV0s,
		       kCascades,
		       kKinks,
		       kCaloClusters,
		       kEMCALCells,
		       kPHOSCells,
		       kErrorLogs,
                       kESDACORDE,
		       kTOFHeader,
		       kESDListN
  };

  AliESDEvent();
  virtual ~AliESDEvent();
  AliESDEvent &operator=(const AliESDEvent& source); // or make private and use only copy? 
  virtual void Copy(TObject& obj) const;

  // RUN
  // move this to the UserData!!!
  const AliESDRun*    GetESDRun() const {return fESDRun;}

  // Delegated methods for fESDRun
  void     SetRunNumber(Int_t n) {if(fESDRun) fESDRun->SetRunNumber(n);}
  Int_t    GetRunNumber() const {return fESDRun?fESDRun->GetRunNumber():-1;}
  void     SetPeriodNumber(UInt_t n){
    if(fESDRun) fESDRun->SetPeriodNumber(n);
    if(fHeader) fHeader->SetPeriodNumber(n);
  }
  UInt_t   GetPeriodNumber() const {return fESDRun?fESDRun->GetPeriodNumber():0;}
  void     SetMagneticField(Double_t mf){if(fESDRun) fESDRun->SetMagneticField(mf);}
  Double_t GetMagneticField() const {return fESDRun?fESDRun->GetMagneticField():0;}
  void     SetDiamond(const AliESDVertex *vertex) { if(fESDRun) fESDRun->SetDiamond(vertex);}
  Double_t  GetDiamondX() const {return fESDRun?fESDRun->GetDiamondX():0;}
  Double_t  GetDiamondY() const {return fESDRun?fESDRun->GetDiamondY():0;}
  Double_t  GetDiamondZ() const {return fESDRun?fESDRun->GetDiamondZ():0;}
  Double_t  GetSigma2DiamondX() const {return  fESDRun?fESDRun->GetSigma2DiamondX():0;}
  Double_t  GetSigma2DiamondY() const {return  fESDRun?fESDRun->GetSigma2DiamondY():0;}
  Double_t  GetSigma2DiamondZ() const {return  fESDRun?fESDRun->GetSigma2DiamondZ():0;}
  void      GetDiamondCovXY(Float_t cov[3]) const {if(fESDRun) fESDRun->GetDiamondCovXY(cov);}   
  void     SetTriggerClass(const char*name, Int_t index) {if(fESDRun) fESDRun->SetTriggerClass(name,index);}
  void     SetPHOSMatrix(TGeoHMatrix*matrix, Int_t i) {if(fESDRun) fESDRun->SetPHOSMatrix(matrix,i);}
  const TGeoHMatrix* GetPHOSMatrix(Int_t i) const {return fESDRun?fESDRun->GetPHOSMatrix(i):0x0;}
  void     SetEMCALMatrix(TGeoHMatrix*matrix, Int_t i) {if(fESDRun) fESDRun->SetEMCALMatrix(matrix,i);}
  const TGeoHMatrix* GetEMCALMatrix(Int_t i) const {return fESDRun?fESDRun->GetEMCALMatrix(i):0x0;}
  void     SetCaloTriggerType(const Int_t* type) {if (fESDRun) fESDRun->SetCaloTriggerType(type);}
  Int_t*   GetCaloTriggerType() const {return fESDRun?fESDRun->GetCaloTriggerType():0x0;};
	
  //
  void        SetCurrentL3(Float_t cur)           const  {if(fESDRun) fESDRun->SetCurrentL3(cur);}
  void        SetCurrentDip(Float_t cur)          const  {if(fESDRun) fESDRun->SetCurrentDip(cur);}
  void        SetBeamEnergy(Float_t be)           const  {if(fESDRun) fESDRun->SetBeamEnergy(be);}
  void        SetBeamType(const char* bt)         const  {if(fESDRun) fESDRun->SetBeamType(bt);}
  void        SetUniformBMap(Bool_t val=kTRUE)    const  {if(fESDRun) fESDRun->SetBit(AliESDRun::kUniformBMap,val);}
  void        SetBInfoStored(Bool_t val=kTRUE)    const  {if(fESDRun) fESDRun->SetBit(AliESDRun::kBInfoStored,val);}
  //
  Float_t     GetCurrentL3()                      const  {return fESDRun?fESDRun->GetCurrentL3():0;}
  Float_t     GetCurrentDip()                     const  {return fESDRun?fESDRun->GetCurrentDip():0;}
  Float_t     GetBeamEnergy()                     const  {return fESDRun?fESDRun->GetBeamEnergy():0;}
  const char* GetBeamType()                       const  {return fESDRun?fESDRun->GetBeamType():0;}
  Bool_t      IsUniformBMap()                     const  {return fESDRun?fESDRun->TestBit(AliESDRun::kUniformBMap):kFALSE;}
  //
  Bool_t      InitMagneticField()                 const  {return fESDRun?fESDRun->InitMagneticField():kFALSE;} 
  void      SetT0spread(Float_t *t)               const  {if(fESDRun) fESDRun->SetT0spread(t);} 
  // HEADER
  AliESDHeader* GetHeader() const {return fHeader;}

  // Delegated methods for fHeader
  void      SetTriggerMask(ULong64_t n) {if(fHeader) fHeader->SetTriggerMask(n);}
  void      SetOrbitNumber(UInt_t n) {if(fHeader) fHeader->SetOrbitNumber(n);}
  void      SetTimeStamp(UInt_t timeStamp){if(fHeader) fHeader->SetTimeStamp(timeStamp);}
  void      SetEventType(UInt_t eventType){if(fHeader) fHeader->SetEventType(eventType);}
  void      SetEventSpecie(UInt_t eventSpecie){if(fHeader) fHeader->SetEventSpecie(eventSpecie);}
  void      SetEventNumberInFile(Int_t n) {if(fHeader) fHeader->SetEventNumberInFile(n);}
  //  void     SetRunNumber(Int_t n) {if(fHeader) fHeader->SetRunNumber(n);}
  void      SetBunchCrossNumber(UShort_t n) {if(fHeader) fHeader->SetBunchCrossNumber(n);}
  void      SetTriggerCluster(UChar_t n) {if(fHeader) fHeader->SetTriggerCluster(n);}
  
  ULong64_t GetTriggerMask() const {return fHeader?fHeader->GetTriggerMask():0;}
  TString   GetFiredTriggerClasses() const {return (fESDRun&&fHeader)?fESDRun->GetFiredTriggerClasses(fHeader->GetTriggerMask()):"";}
  Bool_t    IsTriggerClassFired(const char *name) const {return (fESDRun&&fHeader)?fESDRun->IsTriggerClassFired(fHeader->GetTriggerMask(),name):kFALSE;}
  Bool_t    IsEventSelected(const char *trigExpr) const;
  TObject*  GetHLTTriggerDecision() const;
  TString   GetHLTTriggerDescription() const;
  Bool_t    IsHLTTriggerFired(const char* name=NULL) const;
  UInt_t    GetOrbitNumber() const {return fHeader?fHeader->GetOrbitNumber():0;}
  UInt_t    GetTimeStamp()  const { return fHeader?fHeader->GetTimeStamp():0;}
  UInt_t    GetEventType()  const { return fHeader?fHeader->GetEventType():0;}
  UInt_t    GetEventSpecie()  const { return fHeader?fHeader->GetEventSpecie():0;}
  Int_t     GetEventNumberInFile() const {return fHeader?fHeader->GetEventNumberInFile():-1;}
  UShort_t  GetBunchCrossNumber() const {return fHeader?fHeader->GetBunchCrossNumber():0;}
  UChar_t   GetTriggerCluster() const {return fHeader?fHeader->GetTriggerCluster():0;}

  // ZDC CKB: put this in the header?
  AliESDZDC*    GetESDZDC()  const {return fESDZDC;}
  AliESDZDC*    GetZDCData() const {return fESDZDC;}

  void SetZDCData(AliESDZDC * obj);

  // Delegated methods for fESDZDC
  Double_t GetZDCN1Energy() const {return fESDZDC?fESDZDC->GetZDCN1Energy():0;}
  Double_t GetZDCP1Energy() const {return fESDZDC?fESDZDC->GetZDCP1Energy():0;}
  Double_t GetZDCN2Energy() const {return fESDZDC?fESDZDC->GetZDCN2Energy():0;}
  Double_t GetZDCP2Energy() const {return fESDZDC?fESDZDC->GetZDCP2Energy():0;}
  Double_t GetZDCEMEnergy(Int_t i=0) const {return fESDZDC?fESDZDC->GetZDCEMEnergy(i):0;}
  Int_t    GetZDCParticipants() const {return fESDZDC?fESDZDC->GetZDCParticipants():0;}
  AliCentrality* GetCentrality();
  AliEventplane* GetEventplane();
    

  void     SetZDC(Float_t n1Energy, Float_t p1Energy, Float_t em1Energy, Float_t em2Energy,
                  Float_t n2Energy, Float_t p2Energy, Int_t participants, Int_t nPartA,
	 	  Int_t nPartC, Double_t b, Double_t bA, Double_t bC, UInt_t recoflag)
  {if(fESDZDC) fESDZDC->SetZDC(n1Energy, p1Energy, em1Energy, em2Energy, n2Energy, p2Energy, 
            participants, nPartA, nPartC, b, bA, bC,  recoflag);}
    // FMD
  void SetFMDData(AliESDFMD * obj);
  AliESDFMD *GetFMDData() const { return fESDFMD; }


  // TZERO CKB: put this in the header?
  const AliESDTZERO*    GetESDTZERO() const {return fESDTZERO;}
  void SetTZEROData(AliESDTZERO * obj);
 // delegetated methods for fESDTZERO

  Double32_t GetT0zVertex() const {return fESDTZERO?fESDTZERO->GetT0zVertex():0;}
  void SetT0zVertex(Double32_t z) {if(fESDTZERO) fESDTZERO->SetT0zVertex(z);}
  Double32_t GetT0() const {return fESDTZERO?fESDTZERO->GetT0():0;}
  void SetT0(Double32_t timeStart) {if(fESDTZERO) fESDTZERO->SetT0(timeStart);}
  Double32_t GetT0clock() const {return fESDTZERO?fESDTZERO->GetT0clock():0;}
  void SetT0clock(Double32_t timeStart) {if(fESDTZERO) fESDTZERO->SetT0clock(timeStart);}
  Double32_t GetT0TOF(Int_t icase) const {return fESDTZERO?fESDTZERO->GetT0TOF(icase):0;}
  const Double32_t * GetT0TOF() const {return fESDTZERO?fESDTZERO->GetT0TOF():0x0;}
  void SetT0TOF(Int_t icase,Double32_t timeStart) {if(fESDTZERO) fESDTZERO->SetT0TOF(icase,timeStart);}
  const Double32_t * GetT0time() const {return fESDTZERO?fESDTZERO->GetT0time():0x0;}
  void SetT0time(Double32_t time[24]) {if(fESDTZERO) fESDTZERO->SetT0time(time);}
  const Double32_t * GetT0amplitude() const {return fESDTZERO?fESDTZERO->GetT0amplitude():0x0;}
  void SetT0amplitude(Double32_t amp[24]){if(fESDTZERO) fESDTZERO->SetT0amplitude(amp);}
  Int_t GetT0Trig() const { return fESDTZERO?fESDTZERO->GetT0Trig():0;}
  void SetT0Trig(Int_t tvdc) {if(fESDTZERO) fESDTZERO->SetT0Trig(tvdc);}

  // VZERO 
  AliESDVZERO *GetVZEROData() const { return fESDVZERO; }
  void SetVZEROData(AliESDVZERO * obj);

 // ACORDE
  AliESDACORDE *GetACORDEData() const { return fESDACORDE;}
  void SetACORDEData(AliESDACORDE * obj);

  void SetESDfriend(const AliESDfriend *f) const;
  void GetESDfriend(AliESDfriend *f) const;



  void SetPrimaryVertexTPC(const AliESDVertex *vertex); 
  const AliESDVertex *GetPrimaryVertexTPC() const {return fTPCVertex;}

  void SetPrimaryVertexSPD(const AliESDVertex *vertex); 
  const AliESDVertex *GetPrimaryVertexSPD() const {return fSPDVertex;}
  const AliESDVertex *GetVertex() const {
    //For the backward compatibily only
     return GetPrimaryVertexSPD();
  }

  void SetPrimaryVertexTracks(const AliESDVertex *vertex);
  const AliESDVertex *GetPrimaryVertexTracks() const {return fPrimaryVertex;}
  AliESDVertex *PrimaryVertexTracksUnconstrained() const;

  const AliESDVertex *GetPrimaryVertex() const;



  void SetTOFHeader(const AliTOFHeader * tofEventTime);
  const AliTOFHeader *GetTOFHeader() const {return fTOFHeader;}
  Float_t GetEventTimeSpread() const {if (fTOFHeader) return fTOFHeader->GetT0spread(); else return 0.;}
  Float_t GetTOFTimeResolution() const {if (fTOFHeader) return fTOFHeader->GetTOFResolution(); else return 0.;}


  void SetMultiplicity(const AliMultiplicity *mul);

  const AliMultiplicity *GetMultiplicity() const {return fSPDMult;}
  void   EstimateMultiplicity(Int_t &tracklets,Int_t &trITSTPC,Int_t &trITSSApure,
			      Double_t eta=1.,Bool_t useDCAFlag=kTRUE,Bool_t useV0Flag=kTRUE) const;

  Bool_t Clean(Float_t *cleanPars);
  Bool_t RemoveKink(Int_t i)   const;
  Bool_t RemoveV0(Int_t i)     const;
  Bool_t RemoveTrack(Int_t i)  const;

  const AliESDVertex *GetPileupVertexSPD(Int_t i) const {
    return (const AliESDVertex *)(fSPDPileupVertices?fSPDPileupVertices->UncheckedAt(i):0x0);
  }
  Char_t  AddPileupVertexSPD(const AliESDVertex *vtx);
  const AliESDVertex *GetPileupVertexTracks(Int_t i) const {
    return (const AliESDVertex *)(fTrkPileupVertices?fTrkPileupVertices->UncheckedAt(i):0x0);
  }
  Char_t  AddPileupVertexTracks(const AliESDVertex *vtx);
  TClonesArray* GetPileupVerticesTracks() const {return (TClonesArray*)fTrkPileupVertices;}
  TClonesArray* GetPileupVerticesSPD()    const {return (TClonesArray*)fSPDPileupVertices;}

  virtual Bool_t  IsPileupFromSPD(Int_t minContributors=3, 
				  Double_t minZdist=0.8, 
				  Double_t nSigmaZdist=3., 
				  Double_t nSigmaDiamXY=2., 
				  Double_t nSigmaDiamZ=5.) const;
  
  virtual Bool_t IsPileupFromSPDInMultBins() const;

  AliESDtrack *GetTrack(Int_t i) const {
    return (AliESDtrack *)(fTracks?fTracks->UncheckedAt(i):0x0);
  }
  Int_t  AddTrack(const AliESDtrack *t);

  /// add new track at the end of tracks array and return instance
  AliESDtrack* NewTrack();
  
  AliESDHLTtrack *GetHLTConfMapTrack(Int_t /*i*/) const {
    //    return (AliESDHLTtrack *)fHLTConfMapTracks->UncheckedAt(i);
    return 0;
  }
  void AddHLTConfMapTrack(const AliESDHLTtrack */*t*/) {
    printf("ESD:: AddHLTConfMapTrack do nothing \n");
    //    TClonesArray &fhlt = *fHLTConfMapTracks;
    //  new(fhlt[fHLTConfMapTracks->GetEntriesFast()]) AliESDHLTtrack(*t);
  }
  

  AliESDHLTtrack *GetHLTHoughTrack(Int_t /*i*/) const {
    //    return (AliESDHLTtrack *)fHLTHoughTracks->UncheckedAt(i);
    return 0;
  }
  void AddHLTHoughTrack(const AliESDHLTtrack */*t*/) {
    printf("ESD:: AddHLTHoughTrack do nothing \n");
    //    TClonesArray &fhlt = *fHLTHoughTracks;
    //     new(fhlt[fHLTHoughTracks->GetEntriesFast()]) AliESDHLTtrack(*t);
  }
  
  AliESDMuonTrack *GetMuonTrack(Int_t i) const {
    return (AliESDMuonTrack *)(fMuonTracks?fMuonTracks->UncheckedAt(i):0x0);
  }

  void AddMuonTrack(const AliESDMuonTrack *t);

  AliESDPmdTrack *GetPmdTrack(Int_t i) const {
    return (AliESDPmdTrack *)(fPmdTracks?fPmdTracks->UncheckedAt(i):0x0);
  }

  void AddPmdTrack(const AliESDPmdTrack *t);


  AliESDTrdTrack *GetTrdTrack(Int_t i) const {
    return (AliESDTrdTrack *)(fTrdTracks?fTrdTracks->UncheckedAt(i):0x0);
  }

  
  void AddTrdTrack(const AliESDTrdTrack *t);

  AliESDTrdTracklet* GetTrdTracklet(Int_t idx) const {
    return (AliESDTrdTracklet*)(fTrdTracklets?fTrdTracklets->UncheckedAt(idx):0x0);
  }

  void AddTrdTracklet(const AliESDTrdTracklet *trkl);

  AliESDv0 *GetV0(Int_t i) const {
    return (AliESDv0*)(fV0s?fV0s->UncheckedAt(i):0x0);
  }
  Int_t AddV0(const AliESDv0 *v);

  AliESDcascade *GetCascade(Int_t i) const {
    return (AliESDcascade *)(fCascades?fCascades->UncheckedAt(i):0x0);
  }

  void AddCascade(const AliESDcascade *c);

  AliESDkink *GetKink(Int_t i) const {
    return (AliESDkink *)(fKinks?fKinks->UncheckedAt(i):0x0);
  }
  Int_t AddKink(const AliESDkink *c);

  AliESDCaloCluster *GetCaloCluster(Int_t i) const {
    return (AliESDCaloCluster *)(fCaloClusters?fCaloClusters->UncheckedAt(i):0x0);
  }

  Int_t AddCaloCluster(const AliESDCaloCluster *c);

  AliESDCaloCells *GetEMCALCells() const {return fEMCALCells; }  
  AliESDCaloCells *GetPHOSCells() const {return fPHOSCells; }  

  AliESDCaloTrigger* GetCaloTrigger(TString calo) const 
  {
	  if (calo.Contains("EMCAL")) return fEMCALTrigger;
	  else
		  return fPHOSTrigger;
  }
	
  AliRawDataErrorLog *GetErrorLog(Int_t i) const {
    return (AliRawDataErrorLog *)(fErrorLogs?fErrorLogs->UncheckedAt(i):0x0);
  }
  void  AddRawDataErrorLog(const AliRawDataErrorLog *log) const;

  Int_t GetNumberOfErrorLogs()   const {return fErrorLogs?fErrorLogs->GetEntriesFast():0;}

  Int_t GetNumberOfPileupVerticesSPD() const {
    return (fSPDPileupVertices?fSPDPileupVertices->GetEntriesFast():0);
  }
  Int_t GetNumberOfPileupVerticesTracks() const {
    return (fTrkPileupVertices?fTrkPileupVertices->GetEntriesFast():0);
  }
  Int_t GetNumberOfTracks()     const {return fTracks?fTracks->GetEntriesFast():0;}
  Int_t GetNumberOfHLTConfMapTracks()     const {return 0;} 
  // fHLTConfMapTracks->GetEntriesFast();}
  Int_t GetNumberOfHLTHoughTracks()     const {return  0;  }
  //  fHLTHoughTracks->GetEntriesFast();  }

  Int_t GetNumberOfMuonTracks() const {return fMuonTracks?fMuonTracks->GetEntriesFast():0;}
  Int_t GetNumberOfPmdTracks() const {return fPmdTracks?fPmdTracks->GetEntriesFast():0;}
  Int_t GetNumberOfTrdTracks() const {return fTrdTracks?fTrdTracks->GetEntriesFast():0;}
  Int_t GetNumberOfTrdTracklets() const {return fTrdTracklets?fTrdTracklets->GetEntriesFast():0;}
  Int_t GetNumberOfV0s()      const {return fV0s?fV0s->GetEntriesFast():0;}
  Int_t GetNumberOfCascades() const {return fCascades?fCascades->GetEntriesFast():0;}
  Int_t GetNumberOfKinks() const {return fKinks?fKinks->GetEntriesFast():0;}
  
  Int_t GetEMCALClusters(TRefArray *clusters) const;
  Int_t GetPHOSClusters(TRefArray *clusters) const;
  Int_t GetNumberOfCaloClusters() const {return fCaloClusters?fCaloClusters->GetEntriesFast():0;}

  void SetUseOwnList(Bool_t b){fUseOwnList = b;}
  Bool_t GetUseOwnList() const {return fUseOwnList;}

  void ResetV0s() { if(fV0s) fV0s->Clear(); }
  void ResetCascades() { if(fCascades) fCascades->Clear(); }
  void Reset();

  void  Print(Option_t *option="") const;

  void AddObject(TObject* obj);
  void ReadFromTree(TTree *tree, Option_t* opt = "");
  TObject* FindListObject(const char *name) const;
  AliESD *GetAliESDOld(){return fESDOld;}
  void WriteToTree(TTree* tree) const;
  void GetStdContent();
  void ResetStdContent();
  void CreateStdContent();
  void CreateStdContent(Bool_t bUseThisList);
  void SetStdNames();
  void CopyFromOldESD();
  TList* GetList() const {return fESDObjects;}
  
    //Following needed only for mixed event
  virtual Int_t        EventIndex(Int_t)       const {return 0;}
  virtual Int_t        EventIndexForCaloCluster(Int_t) const {return 0;}
  virtual Int_t        EventIndexForPHOSCell(Int_t)    const {return 0;}
  virtual Int_t        EventIndexForEMCALCell(Int_t)   const {return 0;} 
  

protected:
  AliESDEvent(const AliESDEvent&);
  static Bool_t ResetWithPlacementNew(TObject *pObject);

  TList *fESDObjects;             // List of esd Objects

  AliESDRun       *fESDRun;           //! Run information tmp put in the Userdata
  AliESDHeader    *fHeader;           //! ESD Event Header
  AliESDZDC       *fESDZDC;           //! ZDC information
  AliESDFMD       *fESDFMD;           //! FMD object containing rough multiplicity
  AliESDVZERO     *fESDVZERO;         //! VZERO object containing rough multiplicity
  AliESDTZERO     *fESDTZERO;         //! TZEROObject
  AliESDVertex    *fTPCVertex;        //! Primary vertex estimated by the TPC
  AliESDVertex    *fSPDVertex;        //! Primary vertex estimated by the SPD
  AliESDVertex    *fPrimaryVertex;    //! Primary vertex estimated using ESD tracks
  AliMultiplicity *fSPDMult;          //! SPD tracklet multiplicity
  AliESDCaloTrigger* fPHOSTrigger;     //! PHOS Trigger information
  AliESDCaloTrigger* fEMCALTrigger;    //! PHOS Trigger information
  AliESDACORDE    *fESDACORDE;        //! ACORDE ESD object caontaining bit pattern

  TClonesArray *fSPDPileupVertices;//! Pileup primary vertices reconstructed by SPD 
  TClonesArray *fTrkPileupVertices;//! Pileup primary vertices reconstructed using the tracks 
  TClonesArray *fTracks;           //! ESD tracks 
  TClonesArray *fMuonTracks;       //! MUON ESD tracks
  TClonesArray *fPmdTracks;        //! PMD ESD tracks
  TClonesArray *fTrdTracks;        //! TRD ESD tracks (triggered)
  TClonesArray *fTrdTracklets;     //! TRD tracklets (for trigger)
  TClonesArray *fV0s;              //! V0 vertices
  TClonesArray *fCascades;         //! Cascade vertices
  TClonesArray *fKinks;            //! Kinks
  TClonesArray *fCaloClusters;     //! Calorimeter clusters for PHOS/EMCAL
  AliESDCaloCells *fEMCALCells;     //! EMCAL cell info
  AliESDCaloCells *fPHOSCells;     //! PHOS cell info
  TClonesArray *fErrorLogs;        //! Raw-data reading error messages
 


  AliESD       *fESDOld;           //! Old esd Structure
  AliESDfriend *fESDFriendOld;     //! Old friend esd Structure
  Bool_t    fConnected;            //! flag if leaves are alreday connected
  Bool_t    fUseOwnList;           //! Do not use the list from the esdTree but use the one created by this class 

  static const char* fgkESDListName[kESDListN]; //!

  AliTOFHeader *fTOFHeader;  //! event times (and sigmas) as estimated by TOF
			     //  combinatorial algorithm.
                             //  It contains also TOF time resolution
                             //  and T0spread as written in OCDB
  AliCentrality *fCentrality; //! Centrality for AA collision
  AliEventplane *fEventplane; //! Event plane for AA collision
  ClassDef(AliESDEvent,14)  //ESDEvent class 
};
#endif 

