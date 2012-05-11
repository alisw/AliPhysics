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
  Float_t   beta;
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
  Bool_t filterset1;//TPC  
  Bool_t filterset2;//2010 old
  Bool_t filterset3;//2010 golden


  DeDxTrack();
  
  ClassDef(DeDxTrack, 1);    // Help class
};
//_________________________________________________________
class VZEROCell : public TObject
{
 public:

  Int_t   cellmult;
  Float_t cellindex;
  VZEROCell();
  
  ClassDef(VZEROCell, 1);    // Help class
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
  // new
  Float_t   chi2;
  Float_t   cospt;
  Float_t   dcav0;
  Float_t   dcadaughters;
  Int_t     pdg;
  Short_t   primary;  
  Short_t   status;  
  // old
  DeDxTrack ptrack;
  DeDxTrack ntrack;
  
  DeDxV0();
  
  ClassDef(DeDxV0, 1);    // Help class
};


//_____________________________________________________________________________
class DeDxTrackMC : public TObject
{
 public:
  Float_t pMC;
  Float_t ptMC;
  Float_t etaMC;
  Float_t phiMC;
  Short_t qMC;
  Short_t pidMC;
  Short_t orderMC;
  Int_t   pdgMC;

  DeDxTrackMC();
  
  ClassDef(DeDxTrackMC, 1);    // Help class for MC track debug info
};

//_____________________________________________________________________________
class DeDxEvent : public TObject
{
 public:
  ULong64_t eventid;     // unique event id
  Int_t     run;         // run number
  UInt_t    time;        // time of event
  Float_t   cent;        // centrality
  Float_t   mag;         // magnetic field
  Float_t   zvtx;        // rec vertex
  Float_t   zvtxMC;      // MC true vertes
  Float_t   ptmax;       // Max pt of tracks for this event
  Float_t   ptmaxMC;     // Max pt of MC tracks
  Short_t   vtxstatus;   // Vtx status (-1=no vtx, 0 = outside, 1 = inside cuts)
  Short_t   trackmult;   // Track mult (no cuts)
  Short_t   n;           // Number of added tracks 
  Short_t   trackmultMC; // MC track mult (primary tracks)
  Short_t   nMC;         // MC number of added tracks 
  Short_t   process;     // MC process: -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
  Short_t   trig;        // 0=untriggered, &1 = MB, &2=V0 AND
  Short_t   pileup;      // Is the event marked as pileup?
  
  DeDxEvent();
  
  ClassDef(DeDxEvent, 1);    // Help class
};

//_____________________________________________________________________________
ClassImp(DeDxTrack)

DeDxTrack::DeDxTrack():
TObject(),
  p(-1),
  pt(-1),
//  ptcon(-1),
  pttrue(-1),
//  tpcchi(-1),
  eta(-999),
  phi(-999),
  dedx(-999),
  beta(-999),
  dcaxy(-999),
  dcaz(-999),
  mother(0),
  q(-999),
  filter(-999),
  ncl(-999),
  neff(-999),
  pid(-999),
  primary(-999),
  order(-1),
  filterset1(0),
  filterset2(0),
  filterset3(0)

{
  // default constructor
}
//_____________________________________________________________________________
ClassImp(VZEROCell)

VZEROCell::VZEROCell():
TObject(),
  cellmult(-1),
  cellindex(-999)

{
  // default constructor
}

//_____________________________________________________________________________
ClassImp(DeDxV0)

DeDxV0::DeDxV0():
TObject(),
  p(-1),
  pt(-1),
  eta(-999),
  phi(-999),
  pdca(-1),
  ndca(-1),
  dmassG(-1),
  dmassK0(-1),
  dmassL(-1),
  dmassAL(-1),
  alpha(-999),
  ptarm(-999),
  decayr(-999),
  decayl(-999),
  chi2(-1),
  cospt(-999),
  dcav0(999),
  dcadaughters(999),
  pdg(0),
  primary(-1),  
  status(),  
  ptrack(),
  ntrack()
{
  // default constructor
}

//_____________________________________________________________________________
ClassImp(DeDxTrackMC)

DeDxTrackMC::DeDxTrackMC():
TObject(),
  pMC(-1),
  ptMC(-1),
  etaMC(-999),
  phiMC(-999),
  qMC(-999),
  pidMC(-999),
  orderMC(-1),
  pdgMC(0)
{
  // default constructor
}

//_____________________________________________________________________________
ClassImp(DeDxEvent)

DeDxEvent::DeDxEvent():
TObject(),
  eventid(0),      // unique event id
  run(-1),         // run number
  time(-1),        // time of event
  cent(1000),      // centrality
  mag(+999),       // magnetic field
  zvtx(+999),      // rec vertex
  zvtxMC(+999),    // MC true vertes
  ptmax(-1),       // Max pt of tracks for this event
  ptmaxMC(-1),     // Max pt of MC tracks
  vtxstatus(-2),   // Vtx status (-1=no vtx, 0 = outside, 1 = inside cuts)
  trackmult(-1),   // Track mult (no cuts)
  n(-1),           // Number of added tracks 
  trackmultMC(-1), // MC track mult (primary tracks)
  nMC(-1),         // MC number of added tracks 
  process(-2),     // MC process: -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
  trig(-1),        // Was the event triggered
  pileup(-1)       // Is the event marked as pileup?
{
  // default constructor
}
