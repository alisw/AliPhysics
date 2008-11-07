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

class TList;

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

class AliESDfriend;
class AliESDVZERO;
class AliESDHLTtrack;
class AliESDVertex;
class AliESDPmdTrack;
class AliESDFMD;
class AliESDkink;
class AliESDCaloCluster;
class AliESDCaloCells;
class AliESDv0;
class AliMultiplicity;
class AliRawDataErrorLog;
class AliESDRun;
class AliESDTrdTrack;
class AliESDMuonTrack;
class AliESD;
class AliESDcascade;
class TRefArray;
class AliESDACORDE;

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
		       kTracks,
		       kMuonTracks,
		       kPmdTracks,
		       kTrdTracks,
		       kV0s,
		       kCascades,
		       kKinks,
		       kCaloClusters,
		       kEMCALCells,
		       kPHOSCells,
		       kErrorLogs,
                       kESDACORDE,
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
  void     SetRunNumber(Int_t n) {fESDRun->SetRunNumber(n);}
  Int_t    GetRunNumber() const {return fESDRun->GetRunNumber();}
  void     SetPeriodNumber(UInt_t n){fESDRun->SetPeriodNumber(n);}
  UInt_t   GetPeriodNumber() const {return fESDRun->GetPeriodNumber();}
  void     SetMagneticField(Double_t mf){fESDRun->SetMagneticField(mf);}
  Double_t GetMagneticField() const {return fESDRun->GetMagneticField();}
  void     SetDiamond(const AliESDVertex *vertex) { fESDRun->SetDiamond(vertex);}
  Double_t  GetDiamondX() const {return fESDRun->GetDiamondX();}
  Double_t  GetDiamondY() const {return fESDRun->GetDiamondY();}
  Double_t  GetSigma2DiamondX() const {return  fESDRun->GetSigma2DiamondX();}
  Double_t  GetSigma2DiamondY() const {return  fESDRun->GetSigma2DiamondY();}
  void      GetDiamondCovXY(Float_t cov[3]) const {fESDRun->GetDiamondCovXY(cov);}   
  void     SetTriggerClass(const char*name, Int_t index) {fESDRun->SetTriggerClass(name,index);}
  

  // HEADER
  AliESDHeader* GetHeader() const {return fHeader;}

  // Delegated methods for fHeader
  void      SetTriggerMask(ULong64_t n) {fHeader->SetTriggerMask(n);}
  void      SetOrbitNumber(UInt_t n) {fHeader->SetOrbitNumber(n);}
  void      SetTimeStamp(UInt_t timeStamp){fHeader->SetTimeStamp(timeStamp);}
  void      SetEventType(UInt_t eventType){fHeader->SetEventType(eventType);}
  void      SetEventNumberInFile(Int_t n) {fHeader->SetEventNumberInFile(n);}
  //  void     SetRunNumber(Int_t n) {fHeader->SetRunNumber(n);}
  void      SetBunchCrossNumber(UShort_t n) {fHeader->SetBunchCrossNumber(n);}
  void      SetTriggerCluster(UChar_t n) {fHeader->SetTriggerCluster(n);}
  
  ULong64_t GetTriggerMask() const {return fHeader->GetTriggerMask();}
  TString   GetFiredTriggerClasses() const {return fESDRun->GetFiredTriggerClasses(fHeader->GetTriggerMask());}
  Bool_t    IsTriggerClassFired(const char *name) const {return fESDRun->IsTriggerClassFired(fHeader->GetTriggerMask(),name);}
  UInt_t    GetOrbitNumber() const {return fHeader->GetOrbitNumber();}
  UInt_t    GetTimeStamp()  const { return fHeader->GetTimeStamp();}
  UInt_t    GetEventType()  const { return fHeader->GetEventType();}
  Int_t     GetEventNumberInFile() const {return fHeader->GetEventNumberInFile();}
  UShort_t  GetBunchCrossNumber() const {return fHeader->GetBunchCrossNumber();}
  UChar_t   GetTriggerCluster() const {return fHeader->GetTriggerCluster();}

  // ZDC CKB: put this in the header?
  AliESDZDC*    GetESDZDC() const {return fESDZDC;}

  // Delegated methods for fESDZDC
  Double_t GetZDCN1Energy() const {return fESDZDC->GetZDCN1Energy();}
  Double_t GetZDCP1Energy() const {return fESDZDC->GetZDCP1Energy();}
  Double_t GetZDCN2Energy() const {return fESDZDC->GetZDCN2Energy();}
  Double_t GetZDCP2Energy() const {return fESDZDC->GetZDCP2Energy();}
  Double_t GetZDCEMEnergy(Int_t i=0) const {return fESDZDC->GetZDCEMEnergy(i);}
  Int_t   GetZDCParticipants() const {return fESDZDC->GetZDCParticipants();}
  Int_t   GetZDCParticipants2() const {return fESDZDC->GetZDCParticipants2();}
  void    SetZDC(Float_t n1Energy, Float_t p1Energy, Float_t em1Energy, Float_t em2Energy,
                 Float_t n2Energy, Float_t p2Energy, 
		 Int_t participants, Int_t participants2)
  {fESDZDC->SetZDC(n1Energy, p1Energy, em1Energy, em2Energy, n2Energy, p2Energy, 
   participants, participants2);}


  // FMD
  void SetFMDData(AliESDFMD * obj);
  AliESDFMD *GetFMDData(){ return fESDFMD; }


  // TZERO CKB: put this in the header?
  const AliESDTZERO*    GetESDTZERO() const {return fESDTZERO;}
  // delegetated methods for fESDTZERO

  Double_t GetT0zVertex() const {return fESDTZERO->GetT0zVertex();}
  void SetT0zVertex(Float_t z) {fESDTZERO->SetT0zVertex(z);}
  Double_t GetT0() const {return fESDTZERO->GetT0();}
  void SetT0(Float_t timeStart) {fESDTZERO->SetT0(timeStart);}
  const Double_t * GetT0time() const {return fESDTZERO->GetT0time();}
  void SetT0time(Float_t time[24]) {fESDTZERO->SetT0time(time);}
  const Double_t * GetT0amplitude() const {return fESDTZERO->GetT0amplitude();}
  void SetT0amplitude(Float_t amp[24]){fESDTZERO->SetT0amplitude(amp);}

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

  const AliESDVertex *GetPrimaryVertex() const;

  void SetMultiplicity(const AliMultiplicity *mul);

  const AliMultiplicity *GetMultiplicity() const {return fSPDMult;}


  Bool_t Clean(Float_t *cleanPars);
  Bool_t RemoveKink(Int_t i)   const;
  Bool_t RemoveV0(Int_t i)     const;
  Bool_t RemoveTrack(Int_t i)  const;

  AliESDtrack *GetTrack(Int_t i) const {
    return (AliESDtrack *)fTracks->UncheckedAt(i);
  }
  Int_t  AddTrack(const AliESDtrack *t);

  
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
    return (AliESDMuonTrack *)fMuonTracks->UncheckedAt(i);
  }

  void AddMuonTrack(const AliESDMuonTrack *t);

  AliESDPmdTrack *GetPmdTrack(Int_t i) const {
    return (AliESDPmdTrack *)fPmdTracks->UncheckedAt(i);
  }

  void AddPmdTrack(const AliESDPmdTrack *t);


  AliESDTrdTrack *GetTrdTrack(Int_t i) const {
    return (AliESDTrdTrack *)fTrdTracks->UncheckedAt(i);
  }

  
  void AddTrdTrack(const AliESDTrdTrack *t);

  AliESDv0 *GetV0(Int_t i) const {
    return (AliESDv0*)fV0s->UncheckedAt(i);
  }
  Int_t AddV0(const AliESDv0 *v);

  AliESDcascade *GetCascade(Int_t i) const {
    return (AliESDcascade *)fCascades->UncheckedAt(i);
  }

  void AddCascade(const AliESDcascade *c);

  AliESDkink *GetKink(Int_t i) const {
    return (AliESDkink *)fKinks->UncheckedAt(i);
  }
  Int_t AddKink(const AliESDkink *c);

  AliESDCaloCluster *GetCaloCluster(Int_t i) const {
    return (AliESDCaloCluster *)fCaloClusters->UncheckedAt(i);
  }

  Int_t AddCaloCluster(const AliESDCaloCluster *c);

  AliESDCaloCells *GetEMCALCells() const {return fEMCALCells; }  
  AliESDCaloCells *GetPHOSCells() const {return fPHOSCells; }  

  AliRawDataErrorLog *GetErrorLog(Int_t i) const {
    return (AliRawDataErrorLog *)fErrorLogs->UncheckedAt(i);
  }
  void  AddRawDataErrorLog(const AliRawDataErrorLog *log) const;

  Int_t GetNumberOfErrorLogs()   const {return fErrorLogs->GetEntriesFast();}

    
  void AddPHOSTriggerPosition(TArrayF array)   { fPHOSTrigger->AddTriggerPosition(array); }
  void AddPHOSTriggerAmplitudes(TArrayF array) { fPHOSTrigger->AddTriggerAmplitudes(array);}
  void AddEMCALTriggerPosition(TArrayF array)  { fEMCALTrigger->AddTriggerPosition(array); }
  void AddEMCALTriggerAmplitudes(TArrayF array){ fEMCALTrigger->AddTriggerAmplitudes(array); }

  Int_t GetNumberOfTracks()     const {return fTracks->GetEntriesFast();}
  Int_t GetNumberOfHLTConfMapTracks()     const {return 0;} 
  // fHLTConfMapTracks->GetEntriesFast();}
  Int_t GetNumberOfHLTHoughTracks()     const {return  0;  }
  //  fHLTHoughTracks->GetEntriesFast();  }

  Int_t GetNumberOfMuonTracks() const {return fMuonTracks->GetEntriesFast();}
  Int_t GetNumberOfPmdTracks() const {return fPmdTracks->GetEntriesFast();}
  Int_t GetNumberOfTrdTracks() const {return fTrdTracks->GetEntriesFast();}
  Int_t GetNumberOfV0s()      const {return fV0s->GetEntriesFast();}
  Int_t GetNumberOfCascades() const {return fCascades->GetEntriesFast();}
  Int_t GetNumberOfKinks() const {return fKinks->GetEntriesFast();}
  
  Int_t GetEMCALClusters(TRefArray *clusters) const;
  Int_t GetPHOSClusters(TRefArray *clusters) const;
  Int_t GetNumberOfCaloClusters() const {return fCaloClusters->GetEntriesFast();}

  void SetUseOwnList(Bool_t b){fUseOwnList = b;}
  Bool_t GetUseOwnList(){return fUseOwnList;}
  
  // Remove this stuff CKB?
  //---------------------------------------------------
  Int_t GetNumberOfEMCALClusters() const {return fEMCALClusters;}
  void  SetNumberOfEMCALClusters(Int_t clus) {fEMCALClusters = clus;}
  Int_t GetFirstEMCALCluster() const {return fFirstEMCALCluster;}
  void  SetFirstEMCALCluster(Int_t index) {fFirstEMCALCluster = index;}
 
  Int_t GetNumberOfPHOSClusters() const {return fPHOSClusters;}
  void  SetNumberOfPHOSClusters(Int_t part) { fPHOSClusters = part ; }
  void  SetFirstPHOSCluster(Int_t index) { fFirstPHOSCluster = index ; } 
  Int_t GetFirstPHOSCluster() const  { return fFirstPHOSCluster ; }
  //-------------------------------------------------------

  TArrayF *GetEMCALTriggerPosition() const {return  fEMCALTrigger->GetTriggerPosition();}
  TArrayF *GetEMCALTriggerAmplitudes() const {return  fEMCALTrigger->GetTriggerAmplitudes();}
  TArrayF *GetPHOSTriggerPosition() const {return  fPHOSTrigger->GetTriggerPosition();}
  TArrayF *GetPHOSTriggerAmplitudes() const {return  fPHOSTrigger->GetTriggerAmplitudes();}

  void ResetV0s() { fV0s->Clear(); }
  void ResetCascades() { fCascades->Clear(); }
  void Reset();

  void  Print(Option_t *option="") const;

  void AddObject(TObject* obj);
  void ReadFromTree(TTree *tree, Option_t* opt = "");
  TObject* FindListObject(const char *name);
  AliESD *GetAliESDOld(){return fESDOld;}
  void WriteToTree(TTree* tree) const;
  void GetStdContent();
  void ResetStdContent();
  void CreateStdContent();
  void CreateStdContent(Bool_t bUseThisList);
  void SetStdNames();
  void CopyFromOldESD();
  TList* GetList() const {return fESDObjects;}

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

  TClonesArray *fTracks;           //! ESD tracks 
  TClonesArray *fMuonTracks;       //! MUON ESD tracks
  TClonesArray *fPmdTracks;        //! PMD ESD tracks
  TClonesArray *fTrdTracks;        //! TRD ESD tracks (triggered)
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

  // Remove this stuff CKB
  Int_t        fEMCALClusters;   // Number of EMCAL clusters (subset of caloclusters)
  Int_t        fFirstEMCALCluster; // First EMCAL cluster in the fCaloClusters list 

  Int_t        fPHOSClusters;     // Number of PHOS clusters (subset of caloclusters)
  Int_t        fFirstPHOSCluster; // First PHOS cluster in the fCaloClusters list 

  ClassDef(AliESDEvent,9)  //ESDEvent class 
};
#endif 

