#ifndef ALITRDPID_H
#define ALITRDPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */                   
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   The TRD particle identification base class                              //
//                                                                           //
//   Its main purposes are:                                                  //
//      - Provide I/O framework for all neccessary files                     //
//      - Assignment of a e/pi propability to a given track                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

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

  virtual void          Copy(TObject &p) const;
  virtual Bool_t        Init();
  virtual Bool_t        AssignLikelihood();
  virtual Bool_t        AssignLikelihood(TObjArray *tarray);
  virtual Bool_t        AssignLikelihood(AliTRDtrack *t) = 0;
  virtual Bool_t        FillSpectra();
  virtual Bool_t        FillSpectra(TObjArray *tarray);
  virtual Bool_t        FillSpectra(const AliTRDtrack *t) = 0;
  virtual Bool_t        Open(const Char_t *name, Int_t event = 0);
  virtual Bool_t        Open(const Char_t *namekine
                           , const Char_t *namecluster
                           , const Char_t *nametracks, Int_t event = 0); 
  virtual Int_t         MCpid(const AliTRDtrack *t);
  virtual Int_t         MCpid(const AliTRDtrack *t, Int_t *pdg, Int_t *nFound, Int_t *indices);
  virtual Bool_t        ReadCluster(const Char_t *name);
  virtual Bool_t        ReadTracks(const Char_t *name);
  virtual Bool_t        ReadKine(const Char_t *name, Int_t event);
  virtual Bool_t        SumCharge(const AliTRDtrack *t, Float_t *charge, Int_t *nCluster);

  virtual Int_t         GetIndex(const AliTRDtrack *t) = 0;

          void          SetGeometry(AliTRDgeometry *geo)    { fGeometry      = geo;    };
          void          SetTrackArray(TObjArray *tarray)    { fTrackArray    = tarray; };
          void          SetClusterArray(TObjArray *carray)  { fClusterArray  = carray; };

          void          SetPIDratioMin(Float_t min)         { fPIDratioMin   = min;    };
          void          SetPIDpurePoints(Bool_t pure)       { fPIDpurePoints = pure;   };
          void          SetPIDindexMin(Int_t min)           { fPIDindexMin   = min;    };
          void          SetPIDindexMax(Int_t max)           { fPIDindexMax   = max;    };

          void          SetThreePadOnly(Bool_t only)        { fThreePadOnly  = only;   };
          void          SetEvent(Int_t event)               { fEvent         = event;  };

          TObjArray    *GetTrackArray()               const { return fTrackArray;      }; 
          TObjArray    *GetClusterArray()             const { return fClusterArray;    };

          Float_t       GetPIDratioMin()              const { return fPIDratioMin;     };
          Bool_t        GetPIDpurePoints()            const { return fPIDpurePoints;   };
          Float_t       GetPIDindexMin()              const { return fPIDindexMin;     };
          Float_t       GetPIDindexMax()              const { return fPIDindexMax;     };

          Bool_t        GetThreePadOnly()             const { return fThreePadOnly;    };

 protected:

  enum { 
    kNpid     = 2,                   //  Number of pid types (pion + electron)
    kElectron = 0,                   //  Electron pid
    kPion     = 1                    //  Pion pid
  };

  Float_t         fPIDratioMin;      //  Minimum fraction of cluster from one particle
  Bool_t          fPIDpurePoints;    //  Require pure (nono overlapping) cluster
  Int_t           fPIDindexMin;      //  Lower index MC particles to be considered
  Int_t           fPIDindexMax;      //  Upper index MC particles to be considered

  Bool_t          fThreePadOnly;     //  Use only three pad cluster in the charge sum

  Int_t           fEvent;            //  Event number

  TObjArray      *fTrackArray;       //! Array containing the tracks
  TObjArray      *fClusterArray;     //! Array containing the cluster
  AliTRDgeometry *fGeometry;         //! The TRD geometry
  TFile          *fFileKine;         //! The kine input file

  ClassDef(AliTRDpid,1)              //  Assigns the e/pi propability to the tracks 

};
#endif
