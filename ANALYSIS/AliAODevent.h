#ifndef ALIAODEVENT_H
#define ALIAODEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     Implementation of the Analysis Oriented Data (AOD) event summary
//     Purpose : container of event important information for soft analysis
//     Author : Renaud Vernet, IPHC, Strasbourg
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TClonesArray.h>

class AliESD;
class AliAODv0;
class AliAODxi;

class AliAODevent : public TObject {

 public :
  AliAODevent();
  ~AliAODevent();
  AliAODevent(AliESD* esd);
  AliAODevent(const AliAODevent& aodevent); 
  AliAODevent& operator=(const AliAODevent& aodevent);
  
  void AddV0(AliAODv0* aodv0);
  void AddCascade(AliAODxi* aodxi);
  TClonesArray*  GetV0s() const              {return fV0s;}
  TClonesArray*  GetCascades() const         {return fCascades;}
  AliAODv0*      GetV0    (UInt_t idx) const {return ((AliAODv0*)fV0s->UncheckedAt(idx));}
  AliAODxi*      GetCascade(UInt_t idx) const {return ((AliAODxi*)fCascades->UncheckedAt(idx));}

  UInt_t         GetNumberOfTracks() const   {return fNumberOfTracks;}
  UInt_t         GetNumberOfV0s() const      {return fV0s->GetEntries();}
  UInt_t         GetNumberOfCascades() const {return fCascades->GetEntries();}
  UInt_t         GetRunNumber() const        {return fRunNumber;}
  UInt_t         GetEventNumber() const      {return fEventNumber;}

  Double_t       GetPrimVertexX() const      {return fPrimVertexX;}
  Double_t       GetPrimVertexY() const      {return fPrimVertexY;}
  Double_t       GetPrimVertexZ() const      {return fPrimVertexZ;}

 private :
  TClonesArray* fV0s;         // List of V0's ?
  TClonesArray* fCascades;    // List of cascades ?
  Double_t      fPrimVertexX; // Vertex X coordinate ?
  Double_t      fPrimVertexY; // Vertex Y coordinate ?
  Double_t      fPrimVertexZ; // Vertex Z coordinate ?
  
  UInt_t        fRunNumber;   // Run number
  UInt_t        fEventNumber; // Event number
  UInt_t        fNumberOfTracks; // Number of tracks

  ClassDef(AliAODevent,1);
};

#endif
