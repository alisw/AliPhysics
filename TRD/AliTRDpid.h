#ifndef ALITRDPID_H
#define ALITRDPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>

class TObjArray;
class TFile;

class AliTRDgeometry;
class AliTRDtrack;

class AliTRDpid : public TNamed {

 public:

  AliTRDpid();
  AliTRDpid(const char* name, const char* title);
  AliTRDpid(const AliTRDpid &p);
  virtual ~AliTRDpid();
  AliTRDpid &operator=(const AliTRDpid &p);

  virtual void          Copy(TObject &p);
  virtual Bool_t        Init();
  virtual Bool_t        AssignLQ(TObjArray *tarray);
  virtual Bool_t        AssignLQ(AliTRDtrack *t);
  virtual Bool_t        FillQspectra();
  virtual Bool_t        FillQspectra(const AliTRDtrack *t);
  virtual Bool_t        CreateHistograms(const Int_t nmom
                                       , const Float_t minmom
                                       , const Float_t maxmom);
  virtual Float_t       LQPion(const Float_t *charge);
  virtual Float_t       LQElectron(const Float_t *charge);
  virtual Bool_t        Open(const Char_t *name, Int_t event = 0);
  virtual Bool_t        Open(const Char_t *namekine
                           , const Char_t *namecluster
                           , const Char_t *nametracks, Int_t event = 0); 
  virtual Int_t         Pid(const AliTRDtrack *t);
  virtual Bool_t        ReadCluster(const Char_t *name);
  virtual Bool_t        ReadTracks(const Char_t *name);
  virtual Bool_t        ReadKine(const Char_t *name, Int_t event);
  virtual Bool_t        SumCharge(const AliTRDtrack *t, Float_t *charge);

  inline  Int_t         GetIndexLQ(const Int_t imom, const Int_t ipid);
  inline  Int_t         GetIndexLQ(const Float_t mom, const Int_t ipid);
  inline  Int_t         GetIndexQ(const Int_t imom, const Int_t ipla, const Int_t ipid);   
  inline  Int_t         GetIndexQ(const Float_t mom, const Int_t ipla, const Int_t ipid);

          TObjArray*    GetQHist() const                    { return fQHist;  };
          TObjArray*    GetLQHist() const                   { return fLQHist; };

          void          SetGeometry(AliTRDgeometry *geo)    { fGeometry     = geo;    };
          void          SetTrackArray(TObjArray *tarray)    { fTrackArray   = tarray; };
          void          SetClusterArray(TObjArray *carray)  { fClusterArray = carray; };

 protected:

  enum { 
    kNpid     = 2,                   //  Number of pid types (pion + electron)
    kElectron = 0,                   //  Electron pid
    kPion     = 1                    //  Pion pid
  };

  Int_t           fNMom;             //  Number of momentum bins
  Float_t         fMinMom;           //  Lower momentum
  Float_t         fMaxMom;           //  Upper momentum
  Float_t         fWidMom;           //  Width of the momentum bins
  TObjArray      *fLQHist;           //  Array of L-Q histograms
  TObjArray      *fQHist;            //  Array of Q histograms
  TObjArray      *fTrackArray;       //! Array containing the tracks
  TObjArray      *fClusterArray;     //! Array containing the cluster
  AliTRDgeometry *fGeometry;         //! The TRD geometry

  ClassDef(AliTRDpid,1)              //  Assigns e/pi propability to the tracks 

};
#endif
