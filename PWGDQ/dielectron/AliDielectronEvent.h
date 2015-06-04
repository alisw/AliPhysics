#ifndef ALIDIELECTRONEVENT_H
#define ALIDIELECTRONEVENT_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronEvent                     #
//#                                                           #
//#  Authors:                                                 #
//#   Jens      Wiechula, Uni Tübingen / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

#include <TNamed.h>
#include <TClonesArray.h>

#include "AliDielectronVarManager.h"


class TObjArray;
class TProcessID;

class AliDielectronEvent : public TNamed {
public:
  AliDielectronEvent();
  AliDielectronEvent(const char*name, const char* title);

  virtual ~AliDielectronEvent();

  void SetESD(Int_t sizeP=1000, Int_t sizeN=1000);
  void SetAOD(Int_t sizeP=1000, Int_t sizeN=1000);
  Bool_t IsAOD() const { return fIsAOD; }

  void SetTracks(const TObjArray &arrP, const TObjArray &arrN, const TObjArray &arrPairs);
  void SetEventData(const Double_t data[AliDielectronVarManager::kNMaxValues]);
  const Double_t* GetEventData() const {return fEventData;}
  
  const TClonesArray* GetTrackArrayP() const { return &fArrTrackP; }
  const TClonesArray* GetTrackArrayN() const { return &fArrTrackN; }
  const TClonesArray* GetVertexArray() const { return &fArrVertex; }

  Int_t GetNTracksP() const { return fNTracksP; }
  Int_t GetNTracksN() const { return fNTracksN; }

  void SetProcessID(TProcessID *pid) { fPID=pid;    }
  const TProcessID* GetProcessID()   { return fPID; }
  
virtual void Clear(Option_t *opt="C");


private:
  TClonesArray fArrTrackP;      //positive tracks
  TClonesArray fArrTrackN;      //negative tracks
  TClonesArray fArrVertex;      //track vertices

  TClonesArray fArrPairs;       //Pair array

  Int_t fNTracksP;              //number of positive tracks
  Int_t fNTracksN;              //number of negative tracks

  Bool_t fIsAOD;                // if we deal with AODs

  Double_t fEventData[AliDielectronVarManager::kNMaxValues]; // event informaion from the var manager

  TProcessID *fPID;             //! internal PID for references to buffered objects
  UInt_t      fPIDIndex;        //! index of PID

  AliDielectronEvent(const AliDielectronEvent &c);
  AliDielectronEvent &operator=(const AliDielectronEvent &c);

  void AssignID(TObject *obj);
  
  ClassDef(AliDielectronEvent,1)         // Dielectron Event
};



#endif
