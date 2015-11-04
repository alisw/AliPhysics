// Classes used for creating a reduced information tree
// Author: Paul Christoph Baetzing (pbatzing@cern.ch)
// 
//  Basic structure:
//  1. Event wise information
//  2. List of tracks in the event

#ifndef ALIFILTEREDEVENT_H
#define ALIFILTEREDEVENT_H

#include <TClonesArray.h>
#include <TBits.h>
#include <TMath.h>
class AliFilteredTrack;






//_____________________________________________________________________
class AliFilteredEvent : public TObject {
  public:
  AliFilteredEvent();
  AliFilteredEvent(const char* name);
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
  Double_t GetCentrality(){return fCentrality;}
  void SetRunNr(Int_t run){fRunNr=run;}
  Int_t GetRunNr(){return fRunNr;}
  void SetEventPlane(Double_t plane){fEventPlaneAngle = plane;}
  Double_t GetEventPlaneAngle(){return fEventPlaneAngle;}
  
  void ClearEvent();
  
  //Setting and getting the tracks:
  TClonesArray* GetTracks()   const {return fTracks;}
  AliFilteredTrack* GetTrack(Int_t i)         const 
    {return (i<fNTracks ? (AliFilteredTrack*)fTracks->At(i) : 0x0);}
  void SetNtrks(Int_t nt){fNTracks = nt;}
  Int_t GetNtrks(){return fNTracks;}
private:
  //Event Characteristics:
  Double_t fCentrality;//Centrality of the event
  Double_t fVertex[3];//vertex
  Int_t    fRunNr;
  Double_t fEventPlaneAngle;
  //Track Characteristics:
  Int_t fNTracks;// number of tracks in event
  TClonesArray* fTracks;            //->   array containing global tracks
  static TClonesArray* fgTracks;    //       global tracks  
  ClassDef(AliFilteredEvent, 1)
};
#endif