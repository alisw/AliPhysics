// -*- mode: C++ -*- 
#ifndef ALIESDEVENT_H
#define ALIESDEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliESDEvent
//   This is the class to deal with during the physical analysis of data.
//   It also ensures the backward compatibility with the old ESD format.
//      
// Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------

#include <TClonesArray.h>
#include <TObject.h>
#include <TTree.h>
#include <TArrayF.h>

class TList;


#include "AliESDMuonTrack.h"
#include "AliESDPmdTrack.h"
#include "AliESDTrdTrack.h"
#include "AliESDVertex.h"
#include "AliESDcascade.h"
#include "AliESDkink.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "AliESDv0.h"
#include "AliESDFMD.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"
#include "AliRawDataErrorLog.h"
#include "AliESDRun.h"
#include "AliESDHeader.h"
#include "AliESDZDC.h"
#include "AliESDTZERO.h"
#include "AliESDCaloTrigger.h"

class AliESDfriend;
class AliESDVZERO;
class AliESDHLTtrack;
class AliESDFMD;
class AliESD;

class AliESDEvent : public TObject {
public:


  enum ESDListIndex_t   {kESDRun,
		       kHeader,
		       kESDZDC,
		       kESDFMD,
		       kESDVZERO,
		       kESDTZERO,
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
		       kErrorLogs,
		       kESDListN
  };

  AliESDEvent();
  virtual ~AliESDEvent(); 


  // RUN
  // move this to the UserData!!!
  const AliESDRun*    GetESDRun() const {return fESDRun;}

  // Delegated methods for fESDRun
  void    SetRunNumber(Int_t n) {fESDRun->SetRunNumber(n);}
  Int_t   GetRunNumber() const {return fESDRun->GetRunNumber();}
  void    SetPeriodNumber(Int_t n){fESDRun->SetPeriodNumber(n);}
  Int_t   GetPeriodNumber() const {return fESDRun->GetPeriodNumber();}
  void    SetMagneticField(Float_t mf){fESDRun->SetMagneticField(mf);}
  Float_t GetMagneticField() const {return fESDRun->GetMagneticField();}
  void SetDiamond(const AliESDVertex *vertex) { fESDRun->SetDiamond(vertex);}
  Float_t GetDiamondX() const {return fESDRun->GetDiamondX();}
  Float_t GetDiamondY() const {return fESDRun->GetDiamondY();}
  Float_t GetSigma2DiamondX() const {return  fESDRun->GetSigma2DiamondX();}
  Float_t GetSigma2DiamondY() const {return  fESDRun->GetSigma2DiamondY();}
  void GetDiamondCovXY(Float_t cov[3]) const {fESDRun->GetDiamondCovXY(cov);}   
  

  // HEADER
  const AliESDHeader* GetHeader() const {return fHeader;}

  // Delegated methods for fHeader
  void      SetTriggerMask(ULong64_t n) {fHeader->SetTriggerMask(n);}
  void      SetOrbitNumber(UInt_t n) {fHeader->SetOrbitNumber(n);}
  void      SetTimeStamp(UInt_t timeStamp){fHeader->SetTimeStamp(timeStamp);}
  void      SetEventType(UInt_t eventType){fHeader->SetEventType(eventType);}
  void      SetEventNumberInFile(Int_t n) {fHeader->SetEventNumberInFile(n);}
  //  void      SetRunNumber(Int_t n) {fHeader->SetRunNumber(n);}
  void      SetBunchCrossNumber(UShort_t n) {fHeader->SetBunchCrossNumber(n);}
  void      SetTriggerCluster(UChar_t n) {fHeader->SetTriggerCluster(n);}
  ULong64_t GetTriggerMask() const {return fHeader->GetTriggerMask();}
  UInt_t    GetOrbitNumber() const {return fHeader->GetOrbitNumber();}
  UInt_t    GetTimeStamp()  const { return fHeader->GetTimeStamp();}
  UInt_t    GetEventType()  const { return fHeader->GetEventType();}
  Int_t     GetEventNumberInFile() const {return fHeader->GetEventNumberInFile();}
  UShort_t  GetBunchCrossNumber() const {return fHeader->GetBunchCrossNumber();}
  UChar_t   GetTriggerCluster() const {return fHeader->GetTriggerCluster();}

  // ZDC CKB: put this in the header?
  const AliESDZDC*    GetESDZDC() const {return fESDZDC;}

  // Delegated methods for fESDZDC
  Float_t GetZDCN1Energy() const {return fESDZDC->GetZDCN1Energy();}
  Float_t GetZDCP1Energy() const {return fESDZDC->GetZDCP1Energy();}
  Float_t GetZDCN2Energy() const {return fESDZDC->GetZDCN2Energy();}
  Float_t GetZDCP2Energy() const {return fESDZDC->GetZDCP2Energy();}
  Float_t GetZDCEMEnergy() const {return fESDZDC->GetZDCEMEnergy();}
  Int_t   GetZDCParticipants() const {return fESDZDC->GetZDCParticipants();}
  void    SetZDC(Float_t n1Energy, Float_t p1Energy, Float_t emEnergy,
                 Float_t n2Energy, Float_t p2Energy, Int_t participants)
  {fESDZDC->SetZDC(n1Energy, p1Energy, emEnergy, n2Energy, p2Energy, participants);}


  // FMD
  void SetFMDData(AliESDFMD * obj);
  AliESDFMD *GetFMDData(){ return fESDFMD; }


  // TZERO CKB: put this in the header?
  const AliESDTZERO*    GetESDTZERO() const {return fESDTZERO;}
  // delegetated methods for fESDTZERO

  Float_t GetT0zVertex() const {return fESDTZERO->GetT0zVertex();}
  void SetT0zVertex(Float_t z) {fESDTZERO->SetT0zVertex(z);}
  Float_t GetT0() const {return fESDTZERO->GetT0();}
  void SetT0(Float_t timeStart) {fESDTZERO->SetT0(timeStart);}
  const Float_t * GetT0time() const {return fESDTZERO->GetT0time();}
  void SetT0time(Float_t time[24]) {fESDTZERO->SetT0time(time);}
  const Float_t * GetT0amplitude() const {return fESDTZERO->GetT0amplitude();}
  void SetT0amplitude(Float_t amp[24]){fESDTZERO->SetT0amplitude(amp);}

  // VZERO 
  AliESDVZERO *GetVZEROData() const { return fESDVZERO; }
  void SetVZEROData(AliESDVZERO * obj);


  void SetESDfriend(const AliESDfriend *f);
  void GetESDfriend(AliESDfriend *f) const;



  void SetVertex(const AliESDVertex *vertex); 
  const AliESDVertex *GetVertex() const {return fSPDVertex;}

  void SetPrimaryVertex(const AliESDVertex *vertex);
  const AliESDVertex *GetPrimaryVertex() const {return fPrimaryVertex;}

  void SetMultiplicity(const AliMultiplicity *mul) {
    *fSPDMult = *mul;
    // CKB 
    //     new (&fSPDMult) AliMultiplicity(*mul);
  }
  const AliMultiplicity *GetMultiplicity() const {return fSPDMult;}
  
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
  void AddMuonTrack(const AliESDMuonTrack *t) {
    TClonesArray &fmu = *fMuonTracks;
    new(fmu[fMuonTracks->GetEntriesFast()]) AliESDMuonTrack(*t);
  }

  AliESDPmdTrack *GetPmdTrack(Int_t i) const {
    return (AliESDPmdTrack *)fPmdTracks->UncheckedAt(i);
  }
  void AddPmdTrack(const AliESDPmdTrack *t) {
    TClonesArray &fpmd = *fPmdTracks;
    new(fpmd[fPmdTracks->GetEntriesFast()]) AliESDPmdTrack(*t);
  }

  AliESDTrdTrack *GetTrdTrack(Int_t i) const {
    return (AliESDTrdTrack *)fTrdTracks->UncheckedAt(i);
  }
  void AddTrdTrack(const AliESDTrdTrack *t) {
    TClonesArray &ftrd = *fTrdTracks;
    new(ftrd[fTrdTracks->GetEntriesFast()]) AliESDTrdTrack(*t);
  }

  AliESDv0 *GetV0(Int_t i) const {
    return (AliESDv0*)fV0s->UncheckedAt(i);
  }
  Int_t AddV0(const AliESDv0 *v);

  AliESDcascade *GetCascade(Int_t i) const {
    return (AliESDcascade *)fCascades->UncheckedAt(i);
  }
  void AddCascade(const AliESDcascade *c) {
    TClonesArray &fc = *fCascades;
    new(fc[fCascades->GetEntriesFast()]) AliESDcascade(*c);
  }

  AliESDkink *GetKink(Int_t i) const {
    return (AliESDkink *)fKinks->UncheckedAt(i);
  }
  Int_t AddKink(const AliESDkink *c);

  AliESDCaloCluster *GetCaloCluster(Int_t i) const {
    return (AliESDCaloCluster *)fCaloClusters->UncheckedAt(i);
  }
  Int_t AddCaloCluster(const AliESDCaloCluster *c);

  AliRawDataErrorLog *GetErrorLog(Int_t i) const {
    return (AliRawDataErrorLog *)fErrorLogs->UncheckedAt(i);
  }
  void  AddRawDataErrorLog(const AliRawDataErrorLog *log) {
    // CKB inline this??
    TClonesArray &errlogs = *fErrorLogs;
    new(errlogs[errlogs.GetEntriesFast()])  AliRawDataErrorLog(*log);
  }
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
  Int_t GetNumberOfCaloClusters() const {return fCaloClusters->GetEntriesFast();}

  Int_t GetNumberOfEMCALClusters() const {return fEMCALClusters;}
  void  SetNumberOfEMCALClusters(Int_t clus) {fEMCALClusters = clus;}
  Int_t GetFirstEMCALCluster() const {return fFirstEMCALCluster;}
  void  SetFirstEMCALCluster(Int_t index) {fFirstEMCALCluster = index;}
  TArrayF *GetEMCALTriggerPosition() const {return  fEMCALTrigger->GetTriggerPosition();}
  TArrayF *GetEMCALTriggerAmplitudes() const {return  fEMCALTrigger->GetTriggerAmplitudes();}

  Int_t GetNumberOfPHOSClusters() const {return fPHOSClusters;}
  void  SetNumberOfPHOSClusters(Int_t part) { fPHOSClusters = part ; }
  void  SetFirstPHOSCluster(Int_t index) { fFirstPHOSCluster = index ; } 
  Int_t GetFirstPHOSCluster() const  { return fFirstPHOSCluster ; }
  TArrayF *GetPHOSTriggerPosition() const {return  fPHOSTrigger->GetTriggerPosition();}
  TArrayF *GetPHOSTriggerAmplitudes() const {return  fPHOSTrigger->GetTriggerAmplitudes();}

  void ResetV0s() { fV0s->Clear(); }
  void ResetCascades() { fCascades->Clear(); }
  void Reset();

  void  Print(Option_t *option="") const;

  void AddObject(TObject* obj);
  void ReadFromTree(TTree *tree);
  TObject* FindListObject(const char *name);
  AliESD *GetAliESDOld(){return fESDOld;}
  const void WriteToTree(TTree* tree) const {tree->Branch(fESDObjects);}
  void GetStdContent();
  void ResetStdContent();
  void CreateStdContent();
  void SetStdNames();
  void CopyFromOldESD();
  TList* GetList(){return fESDObjects;}

protected:
  AliESDEvent(const AliESDEvent&);
  AliESDEvent &operator=(const AliESDEvent& source);


  TList *fESDObjects;             // List of esd Objects

  AliESDRun       *fESDRun;           //! Run information tmp put in the Userdata
  AliESDHeader    *fHeader;           //! ESD Event Header
  AliESDZDC       *fESDZDC;           //! ZDC information
  AliESDFMD       *fESDFMD;           //! FMD object containing rough multiplicity
  AliESDVZERO     *fESDVZERO;         //! VZERO object containing rough multiplicity
  AliESDTZERO     *fESDTZERO;         //! TZEROObject
  AliESDVertex    *fSPDVertex;        //! Primary vertex estimated by the SPD
  AliESDVertex    *fPrimaryVertex;    //! Primary vertex estimated using ESD tracks
  AliMultiplicity *fSPDMult;          //! SPD tracklet multiplicity
  AliESDCaloTrigger* fPHOSTrigger;     //! PHOS Trigger information
  AliESDCaloTrigger* fEMCALTrigger;    //! PHOS Trigger information

  TClonesArray *fTracks;           //! ESD tracks 
  TClonesArray *fMuonTracks;       //! MUON ESD tracks
  TClonesArray *fPmdTracks;        //! PMD ESD tracks
  TClonesArray *fTrdTracks;        //! TRD ESD tracks (triggered)
  TClonesArray *fV0s;              //! V0 vertices
  TClonesArray *fCascades;         //! Cascade vertices
  TClonesArray *fKinks;            //! Kinks
  TClonesArray *fCaloClusters;     //! Calorimeter clusters for PHOS/EMCAL
  TClonesArray *fErrorLogs;        //! Raw-data reading error messages
 


  AliESD    *fESDOld;              //! Old esd Structure
  Bool_t    fConnected;            //! flag if leaves are alreday connected

  static const char* fESDListName[kESDListN];

  // Remove this stuff CKB
  Int_t        fEMCALClusters;   // Number of EMCAL clusters (subset of caloclusters)
  Int_t        fFirstEMCALCluster; // First EMCAL cluster in the fCaloClusters list 

  Int_t        fPHOSClusters;     // Number of PHOS clusters (subset of caloclusters)
  Int_t        fFirstPHOSCluster; // First PHOS cluster in the fCaloClusters list 

  ClassDef(AliESDEvent,2)  //ESDEvent class 
};
#endif 

