
#ifndef DEBUGCLASSESMULTESA2013_H
#define DEBUGCLASSESMULTESA2013_H

class DeDxTrack : public TObject
{
 public:
  Float_t   p;
  Float_t   pt;
  //  Float_t   ptcon;
  Float_t   pttrue;
  //  Float_t   tpcchi;
  Float_t   eta;
  Float_t   phi;
  Float_t   dedx;
  Float_t   protNSigma;
  Float_t   pionNSigma;

  Float_t   dcaxy;
  Float_t   dcaz;
  Int_t     mother; // pdg of mother (can be same particle)
  Short_t   q;
  Short_t   filter;
  Short_t   ncl;
  Short_t   neff;
  Short_t   pid;
  Short_t   primary;  
  Short_t   order;
  Int_t     tpcnclS; //number of shared TPC clusters

  DeDxTrack();
  void Copy(TObject& object) const;

  ClassDef(DeDxTrack, 2);    // Help class
};
//_________________________________________________________
class VZEROCell : public TObject
{
 public:

  Int_t   cellmult;
  Float_t cellindex;
  VZEROCell();
  void Copy(TObject& object) const;
  
  ClassDef(VZEROCell, 2);    // Help class
};


//_____________________________________________________________________________
class DeDxV0 : public TObject
{
 public:
  Float_t   p;
  Float_t   pt;
  Float_t   eta;
  Float_t   phi;
  Float_t   pdca;     // Distance of Closest Approach for positive track
  Float_t   ndca;     // Distance of Closest Approach for positive track
  Float_t   dmassG;
  Float_t   dmassK0;
  Float_t   dmassL;
  Float_t   dmassAL;
  Float_t   alpha;
  Float_t   ptarm;
  Float_t   decayr;
  Float_t   decayl;
  Float_t   chi2;
  Float_t   cospt;
  Float_t   dcav0;
  Float_t   dcadaughters;
  Int_t     pdg;
  Int_t     pdgmother;
  Short_t   oobPileupFlag;
  Short_t   primary;  
  Short_t   status;
  DeDxTrack ptrack;
  DeDxTrack ntrack;
   Float_t   y;
  DeDxV0();
  void Copy(TObject& object) const;

  
  ClassDef(DeDxV0, 3);    // Help class
};


//_____________________________________________________________________________
class DeDxTrackMC : public TObject
{
 public:
  Float_t pMC;
  Float_t ptMC;
  Float_t etaMC;
  Float_t phiMC;
  Float_t yMC;
  Short_t qMC;
  Short_t pidMC;
  Short_t orderMC;
  Int_t   pdgMC;

  DeDxTrackMC();
  void Copy(TObject& object) const;

  ClassDef(DeDxTrackMC, 2);    // Help class for MC track debug info
};

//_____________________________________________________________________________
class DeDxEvent : public TObject
{
 public:
  ULong64_t eventid;     // unique event id
  Int_t     run;         // run number
  UInt_t    time;        // time of event
  Float_t   cent;        // centrality V0A+V0C, default
  Float_t   mag;         // magnetic field
  Float_t   zvtx;        // rec vertex
  Float_t   zvtxMC;      // MC true vertes
  Float_t   ptmax;       // Max pt of tracks for this event
  Float_t   ptmaxMC;     // Max pt of MC tracks
  Short_t   vtxstatus;   // Vtx status (-1=no vtx, 0 = outside, 1 = inside cuts)
  Short_t   trackmult;   // Track mult (no cuts)
  Short_t   n;           // Number of added tracks 
  Short_t   nTracks;
  UInt_t    refMult;     // reference multiplicity
  Short_t   trackmultMC; // MC track mult (primary tracks)
  Short_t   nMC;         // MC number of added tracks 
  Short_t   process;     // MC process: -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
  Short_t   trig;        // 0=untriggered, &1 = MB, &2=V0 AND
  Int_t     triggerInt;  // 0 = kMB, 1 = kCent, 2 = kSemiCent
  Int_t     v0Finder;    // 0 = oldFinder, 1 = newFinder
  Int_t     centFramework; // 0 = AliCentrality, 1 = AliMultSelection
 
  DeDxEvent();
  void Copy(TObject& object) const;

  ClassDef(DeDxEvent, 4);    // Help class
};

#endif
