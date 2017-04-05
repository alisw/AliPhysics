// Classes used for creating a reduced information tree
// Author: Paul Christoph Baetzing (pbatzing@cern.ch)
// 
//  Basic structure:
//  1. Event wise information
//  2. List of tracks in the event

#ifndef ALIFILTEREDEVENT_H
#define ALIFILTEREDEVENT_H

#include <TClonesArray.h>
#include <AliVEvent.h>
#include <TBits.h>
#include <TMath.h>
#include <AliLog.h>
class AliFilteredTrack;






//_____________________________________________________________________
class AliFilteredEvent : public AliVEvent {
  public:
  AliFilteredEvent();
  AliFilteredEvent(const AliFilteredEvent& other);
  AliFilteredEvent(const char* name);
  AliFilteredEvent& operator=(const AliFilteredEvent& other);
  
  ~AliFilteredEvent();
  
  //Getting and Setting the event Characteristics:
  void SetVertex(Double_t vert[3]){
    fVertex[0]=vert[0];
    fVertex[1]=vert[1];
    fVertex[2]=vert[2];    
  }
  Double_t GetfVertexX(){return fVertex[0];}
  Double_t GetfVertexY(){return fVertex[1];}
  Double_t GetfVertexZ(){return fVertex[2];}
  void SetCentrality(Double_t cent){fCentrality=cent;}
  Double_t GetCentralityP(){return fCentrality;}
  void SetRunNr(Int_t run){fRunNr=run;}
  void     SetRunNumber(Int_t n){SetRunNr(n);}
  Int_t GetRunNr(){return fRunNr;}
  Int_t    GetRunNumber() const{return fRunNr;}
  void SetEventPlane(Double_t plane){fEventPlaneAngle = plane;}
  Double_t GetEventPlaneAngle(){return fEventPlaneAngle;}
  void   SetConnected(Bool_t conn=kTRUE) {fConnected=conn;}
  Bool_t GetConnected() const {return fConnected;}
  void ClearEvent();
  void Reset(){ClearEvent();}
  
  //Setting and getting the tracks:
  TClonesArray* GetTracks()   const {return fTracks;}
  AliVParticle* GetTrack(Int_t i)const {return (i<fNTracks ? (AliVParticle*)fTracks->At(i) : 0x0);}
  void SetNtrks(Int_t nt){fNTracks = nt;}
  Int_t GetNtrks(){return fNTracks;}
  Int_t  GetNumberOfTracks() const{return fNTracks;}
  void ReadFromTree(TTree *tree, Option_t* opt = "");
  
  //empty implementations of AliVEvent virtual functions to stop compiler from complaining:
  virtual AliVHeader* GetHeader() const {return NULL;}
  void AddObject(TObject* obj){;}
  TObject* FindListObject(const char *name) const {return NULL;}
  TList* GetList() const {return NULL;}
  void CreateStdContent(){;}
  void GetStdContent(){;}
  void WriteToTree(TTree* tree) const{;}
  void SetStdNames(){;}
  void Print(Option_t *option="") const{;}
  void     SetPeriodNumber(UInt_t n){;}
  void     SetMagneticField(Double_t mf){;}
  UInt_t   GetPeriodNumber() const {return 0;}
  Double_t GetMagneticField() const{return 0;}
  void      SetOrbitNumber(UInt_t n){;}
  void      SetBunchCrossNumber(UShort_t n){;}
  void      SetEventType(UInt_t eventType){;}
  void      SetTriggerMask(ULong64_t n){;}
  void      SetTriggerCluster(UChar_t n){;}
  UInt_t    GetOrbitNumber() const{return 0;}
  UShort_t  GetBunchCrossNumber() const{return 0;}
  UInt_t    GetEventType()  const {return 0;}
  ULong64_t GetTriggerMask() const{return 0;}
  UChar_t   GetTriggerCluster() const{return UChar_t(0);}
  TString   GetFiredTriggerClasses() const {return TString("");}
  Double_t  GetZDCN1Energy() const{return 0;}
  Double_t  GetZDCP1Energy() const{return 0;}
  Double_t  GetZDCN2Energy() const{return 0;}
  Double_t  GetZDCP2Energy() const{return 0;}
  Double_t  GetZDCEMEnergy(Int_t i) const{return 0;}
  AliVVZERO *GetVZEROData() const {return NULL;}
  AliVZDC   *GetZDCData() const{return NULL;}
  Int_t        GetNumberOfV0s() const{return 0;}
  Int_t        GetNumberOfCascades() const{return 0;}
  AliCentrality* GetCentrality(){return NULL;}
  AliEventplane* GetEventplane(){return NULL;}
  Int_t        EventIndex(Int_t itrack)             const{return 0;}
  Int_t        EventIndexForCaloCluster(Int_t iclu) const{return 0;}
  Int_t        EventIndexForPHOSCell(Int_t icell)   const{return 0;}
  Int_t        EventIndexForEMCALCell(Int_t icell)  const{return 0;}
  EDataLayoutType GetDataLayoutType() const{return EDataLayoutType();}
  
private:
  //Event Characteristics:
  Double_t fCentrality;//Centrality of the event
  Double_t fVertex[3];//vertex
  Int_t    fRunNr;
  Double_t fEventPlaneAngle;
  Bool_t   fConnected;  //! flag if leaves are alreday connected 
  //Track Characteristics:
  Int_t fNTracks;// number of tracks in event
  TClonesArray* fTracks;            //->   array containing global tracks
  static TClonesArray* fgTracks;    //       global tracks  
  ClassDef(AliFilteredEvent, 1)
};
#endif