#ifndef ALIFOURPIONEVENTCOLLECTION
#define ALIFOURPIONEVENTCOLLECTION

#include <iostream>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TBits.h"
#include "TObject.h"
#include "TVector2.h"
#include "AliESDtrack.h"

using namespace std;


class AliFourPionTrackStruct{// detector tracks
  
 public:
  AliFourPionTrackStruct();
  virtual ~AliFourPionTrackStruct();
  AliFourPionTrackStruct(const AliFourPionTrackStruct &obj); 
  AliFourPionTrackStruct &operator=(const AliFourPionTrackStruct &obj);


  UInt_t fStatus;
  UInt_t fFiltermap;
  Int_t fId;
  Double_t fPhi;
  Float_t fPt;
  Float_t fMom;
  Double_t fP[3];
  Int_t fCharge;
  Float_t fEta;
  Float_t fMass;
  Float_t fDCAXY;
  Float_t fDCAZ;
  Float_t fDCA;
  Float_t fEaccepted;
  Short_t fKey;
  TBits fClusterMap;
  TBits fSharedMap;
  Double_t fX[3];
  Bool_t fTOFhit;
  Bool_t fElectron;
  Bool_t fMuon;
  Bool_t fPion;
  Bool_t fKaon;
  Bool_t fProton;
  Int_t fLabel;// MC

  ClassDef(AliFourPionTrackStruct, 1);
};

class AliFourPionMCStruct{// MC info

 public:
  AliFourPionMCStruct();
  virtual ~AliFourPionMCStruct();
  AliFourPionMCStruct(const AliFourPionMCStruct &obj); 
  AliFourPionMCStruct &operator=(const AliFourPionMCStruct &obj);

  Float_t fPx;
  Float_t fPy;
  Float_t fPz;
  Float_t fPtot;
  Int_t fPdgCode;
  Int_t fMotherLabel;

  ClassDef(AliFourPionMCStruct, 1);
};

class AliFourPionEventStruct{// like particle_event
  
 public:
  AliFourPionEventStruct();
  virtual ~AliFourPionEventStruct();
  AliFourPionEventStruct(const AliFourPionEventStruct &obj); 
  AliFourPionEventStruct &operator=(const AliFourPionEventStruct &obj);


  Int_t fFillStatus;
  Int_t fNtracks;
  Int_t fMCarraySize;
  AliFourPionTrackStruct *fTracks;
  AliFourPionMCStruct *fMCtracks;

  ClassDef(AliFourPionEventStruct, 1);
};



class AliFourPionEventCollection {
  
  public:
    AliFourPionEventCollection();
    AliFourPionEventCollection(Short_t,Int_t,Int_t,Bool_t);
    virtual ~AliFourPionEventCollection();
    AliFourPionEventCollection(const AliFourPionEventCollection &obj); 
    AliFourPionEventCollection &operator=(const AliFourPionEventCollection &obj);
   
    void FIFOShift();
    void SetBuffSize(Short_t a){fFIFO = a;}
 
    Short_t fFIFO; //Size of the Event Storage buffer.
    Int_t fLimit; //Max number of tracks
    Int_t fMCLimit; //Max number of MC tracks
    AliFourPionEventStruct *fEvtStr;

   ClassDef(AliFourPionEventCollection, 1);
};
#endif
