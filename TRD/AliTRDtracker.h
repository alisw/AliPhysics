#ifndef ALITRDTRACKER_H
#define ALITRDTRACKER_H   

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

#include <TNamed.h>

class TFile;
class TObjArray;

class AliTRDgeometry;
// class AliTRDtrackingSector;
class AliTRDtrack;
class AliTRDseed;


class AliTRDtracker : public TNamed { 

 public:

  AliTRDtracker();
  AliTRDtracker(const Text_t* name, const Text_t* title);
  ~AliTRDtracker(); 

  virtual void  GetEvent(const Char_t *name, Int_t nEvent = 0);
  virtual void  SetUpSectors(AliTRDtrackingSector *sec);
  virtual void  MakeSeeds(Int_t inner, Int_t outer);
  virtual void  FindTracks();
  virtual void  UseClusters(AliTRDseed t);
  virtual Int_t GetTrackLabel(AliTRDseed t);
  virtual Int_t WriteTracks(); 
  virtual void  ReadClusters(TObjArray *array, const Char_t *filename, Int_t nEvent = 0, Int_t option = 1);

 protected:

  TFile            *fInputFile;       // AliROOT input file
  AliTRDgeometry   *fGeom;            // Pointer to TRD geometry
  Int_t            fEvent;            // Event number

  Int_t            fNclusters;        // Number of clusters in TRD 
  TObjArray        *fClusters;        // List of clusters for all sectors

  Int_t            fNseeds;           // Number of track seeds  
  TObjArray        *fSeeds;           // List of track seeds
   
  Int_t            fNtracks;          // Number of reconstructed tracks 
  TObjArray        *fTracks;          // List of reconstructed tracks   

  ClassDef(AliTRDtracker,1)           // manager base class  

};

#endif 
