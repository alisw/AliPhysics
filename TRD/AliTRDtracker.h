#ifndef ALITRDTRACKER_H
#define ALITRDTRACKER_H   

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

#include <TNamed.h>
#include <TH1.h>   

class TFile;
class TParticle;
class TParticlePDG;
class TObjArray;

class AliTRDgeometry;
class AliTRDtrack;
class AliTRDmcTrack;
class AliTRDtrackingSector;

class AliTRDtracker : public TNamed { 

 public:

  AliTRDtracker();
  AliTRDtracker(const Text_t* name, const Text_t* title);
  ~AliTRDtracker(); 

  virtual void  Clusters2Tracks(TH1F *hs, TH1F *hd); 
  Double_t      ExpectedSigmaY2(Double_t r, Double_t tgl, Double_t pt);
  Double_t      ExpectedSigmaZ2(Double_t r, Double_t tgl);
  Int_t         FindProlongation(AliTRDtrack& t, AliTRDtrackingSector *sec,
                              Int_t s, Int_t rf=0, Int_t matched_index = -1,
				 TH1F *hs=0, TH1F *hd=0);
  void          GetEvent(const Char_t *hitfile, const Char_t *clusterfile);
  void          SetUpSectors(AliTRDtrackingSector *sec);
  virtual void  MakeSeeds(Int_t inner, Int_t outer, AliTRDtrackingSector *sec,
			  Int_t turn, TH1F *hs, TH1F *hd);
  virtual void  FindTracks(AliTRDtrackingSector *sec, TH1F *hs, TH1F *hd);
  virtual void  UseClusters(AliTRDtrack t);
  virtual Int_t GetTrackLabel(AliTRDtrack t);
  Int_t         WriteTracks(const Char_t *filename); 
  void          ReadClusters(TObjArray *array, const Char_t *filename);

  Float_t  GetSeedGap()       const {return fSeedGap;}   
  Float_t  GetSeedStep()      const {return fSeedStep;}
  Float_t  GetSeedDepth()     const {return fSeedDepth;}
  Float_t  GetSkipDepth()     const {return fSkipDepth;}
  Double_t GetMaxChi2()       const {return fMaxChi2;}
  Float_t  GetMaxSeedC()      const {return fMaxSeedC;}
  Float_t  GetMaxSeedTan()    const {return fMaxSeedTan;}
  Double_t GetSeedErrorSY()   const {return fSeedErrorSY;}
  Double_t GetSeedErrorSY3()  const {return fSeedErrorSY3;}
  Double_t GetSeedErrorSZ()   const {return fSeedErrorSZ;}
  Float_t  GetLabelFraction() const {return fLabelFraction;}
  Float_t  GetWideRoad()      const {return fWideRoad;}

  Float_t  GetMinClustersInTrack() const {return fMinClustersInTrack;}
  Float_t  GetMinClustersInSeed()  const {return fMinClustersInSeed;} 
  Float_t  GetMaxSeedDeltaZ()      const {return fMaxSeedDeltaZ;}
  Float_t  GetMaxSeedVertexZ()     const {return fMaxSeedVertexZ;}

  void     SetSY2corr(Float_t w)    {fSY2corr = w;}

 protected:

  Int_t            fEvent;            // Event number

  AliTRDgeometry   *fGeom;            // Pointer to TRD geometry

  Int_t            fNclusters;        // Number of clusters in TRD 
  TObjArray        *fClusters;        // List of clusters for all sectors

  Int_t            fNseeds;           // Number of track seeds  
  TObjArray        *fSeeds;           // List of track seeds
   
  Int_t            fNtracks;          // Number of reconstructed tracks 
  TObjArray        *fTracks;          // List of reconstructed tracks   

  Float_t          fSY2corr;          // Correction coefficient for
                                      // cluster SigmaY2 

  static const Float_t  fSeedGap;   // Distance between inner and outer
                                    // time bin in seeding 
				    // (fraction of all time bins) 
  
  static const Float_t  fSeedStep;    // Step in iterations
  static const Float_t 	fSeedDepth;   // Fraction of TRD allocated for seeding
  static const Float_t  fSkipDepth;   // Fraction of TRD which can be skipped
                                      // in track prolongation		   
  static const Double_t fMaxChi2;     // max increment in track chi2 
 	
  static const Float_t  fMinClustersInTrack; // min fraction of clusters in track
  static const Float_t  fMinClustersInSeed;  // min fraction of clusters in seed
  static const Float_t  fMaxSeedDeltaZ;   // max dZ in MakeSeeds
  static const Float_t  fMaxSeedDeltaZ12; // max abs(z1-z2) in MakeSeeds
  static const Float_t  fMaxSeedC;       // max initial curvature in MakeSeeds
  static const Float_t  fMaxSeedTan;     // max initial Tangens(lambda) in MakeSeeds
  static const Float_t  fMaxSeedVertexZ; // max vertex Z in MakeSeeds
  static const Double_t fSeedErrorSY;    // sy parameter in MakeSeeds
  static const Double_t fSeedErrorSY3;   // sy3 parameter in MakeSeeds
  static const Double_t fSeedErrorSZ;    // sz parameter in MakeSeeds
  static const Float_t  fLabelFraction;  // min fraction of clusters in GetTrackLabel
  static const Float_t  fWideRoad;       // max road width in FindProlongation
 
  ClassDef(AliTRDtracker,1)           // manager base class  

};

#endif 
