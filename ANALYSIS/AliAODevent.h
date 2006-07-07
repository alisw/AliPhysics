#ifndef ALIAODEVENT_H
#define ALIAODEVENT_H

#include <TObject.h>
#include <TClonesArray.h>

class AliESD;
class AliAODv0;
class AliAODxi;

class AliAODevent : public TObject {

 private :
  TClonesArray* fV0s;
  TClonesArray* fCascades;
  Double_t      fPrimVertexX; // here put the whole AliESDVertex ?
  Double_t      fPrimVertexY;
  Double_t      fPrimVertexZ;
  
  UInt_t        fRunNumber;
  UInt_t        fEventNumber;
  UInt_t        fNumberOfTracks;

 public :
  AliAODevent();
  ~AliAODevent();
  AliAODevent(AliESD*);
  AliAODevent(const AliAODevent&); 
  AliAODevent& operator=(const AliAODevent&);
  
  void AddV0(AliAODv0*);
  void AddCascade(AliAODxi*);
  inline TClonesArray*  GetV0s() {return fV0s;}
  inline TClonesArray*  GetCascades() {return fCascades;}
  inline AliAODv0*      GetV0     (UInt_t idx) {return ((AliAODv0*)fV0s->UncheckedAt(idx));}
  inline AliAODxi*      GetCascade(UInt_t idx) {return ((AliAODxi*)fCascades->UncheckedAt(idx));}

  inline UInt_t         GetNumberOfTracks()   {return fNumberOfTracks;}
  inline UInt_t         GetNumberOfV0s()      {return fV0s->GetEntries();}
  inline UInt_t         GetNumberOfCascades() {return fCascades->GetEntries();}
  inline UInt_t         GetRunNumber()        {return fRunNumber;}
  inline UInt_t         GetEventNumber()      {return fEventNumber;}

  inline Double_t       GetPrimVertexX()      {return fPrimVertexX;}
  inline Double_t       GetPrimVertexY()      {return fPrimVertexY;}
  inline Double_t       GetPrimVertexZ()      {return fPrimVertexZ;}

  ClassDef(AliAODevent,1);
};

#endif
