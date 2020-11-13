// ROOT includes
#include <TObject.h>

#include "DebugClassesMultESA2013.h"

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
  protNSigma(-999),
  pionNSigma(-999),
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
  tpcnclS(0)

{
  // default constructor
}

void DeDxTrack::Copy(TObject& object) const
{
  TObject::Copy(object);

  DeDxTrack* track = (DeDxTrack*)(&object);
  if(!track)
    return;
  
  track->p	    = p;	       
  track->pt	    = pt;	       
  track->pttrue     = pttrue;    
  track->eta        = eta;       
  track->phi        = phi;       
  track->dedx       = dedx;
  track->protNSigma = protNSigma;
  track->pionNSigma = pionNSigma;    
  track->dcaxy      = dcaxy;     
  track->dcaz       = dcaz;      
  track->mother     = mother;    
  track->q	    = q;	       
  track->filter     = filter;    
  track->ncl        = ncl;       
  track->neff       = neff;      
  track->pid        = pid;       
  track->primary    = primary;   
  track->order      = order;     
  track->tpcnclS    = tpcnclS;
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

void VZEROCell::Copy(TObject& object) const
{
  TObject::Copy(object);

  VZEROCell* cellv0 = (VZEROCell*)(&object);
  if(!cellv0)
    return;
  
  cellv0->cellmult  = cellmult;	       
  cellv0->cellindex = cellindex;	       
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
  pdgmother(0),
  oobPileupFlag(-1),
  primary(-1),  
  status(),  
  ptrack(),
  ntrack(),
  y(-999)
{
  // default constructor
}

void DeDxV0::Copy(TObject& object) const
{
  TObject::Copy(object);

  DeDxV0* v0 = (DeDxV0*)(&object);
  if(!v0)
    return;

  v0->p	            = p;		 
  v0->pt	    = pt;		 
  v0->eta	    = eta;		 
  v0->phi	    = phi;		 
  v0->pdca	    = pdca;		 
  v0->ndca	    = ndca;		 
  v0->dmassG	    = dmassG;	 
  v0->dmassK0	    = dmassK0;	 
  v0->dmassL	    = dmassL;	 
  v0->dmassAL	    = dmassAL;	 
  v0->alpha	    = alpha;	 
  v0->ptarm	    = ptarm;	 
  v0->decayr	    = decayr;	 
  v0->decayl	    = decayl;	 
  v0->chi2	    = chi2;		 
  v0->cospt	    = cospt;	 
  v0->dcav0	    = dcav0;	 
  v0->dcadaughters = dcadaughters;	 
  v0->pdg	    = pdg;		 
  v0->pdgmother	    = pdgmother;
  v0->oobPileupFlag = oobPileupFlag;
  v0->primary       = primary;  	 
  v0->status  	    = status;  	 

  ptrack.Copy(v0->ptrack);	 
  ntrack.Copy(v0->ntrack);        
  v0->y = y;
}

//_____________________________________________________________________________
ClassImp(DeDxTrackMC)

DeDxTrackMC::DeDxTrackMC():
TObject(),
  pMC(-1),
  ptMC(-1),
  etaMC(-999),
  phiMC(-999),
  yMC(-999),
  qMC(-999),
  pidMC(-999),
  orderMC(-1),
  pdgMC(0)
{
  // default constructor
}

void DeDxTrackMC::Copy(TObject& object) const
{
  TObject::Copy(object);

  DeDxTrackMC* trackmc = (DeDxTrackMC*)(&object);
  if(!trackmc)
    return;
  
  trackmc->pMC	    = pMC;	       
  trackmc->ptMC	    = ptMC;	       
  trackmc->etaMC    = etaMC;
  trackmc->phiMC    = phiMC;
  trackmc->yMC      = yMC;
  trackmc->qMC      = qMC;
  trackmc->pidMC    = pidMC;
  trackmc->orderMC  = orderMC;
  trackmc->pdgMC    = pdgMC;
 
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
  nTracks(-1),
  refMult(-1),     // reference multiplicity
  trackmultMC(-1), // MC track mult (primary tracks)
  nMC(-1),         // MC number of added tracks 
  process(-2),     // MC process: -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
  trig(-1),        // Was the event triggered
  triggerInt(-1),     // 0 = kMB, 1 = kCent, 2 = kSemiCent
  v0Finder(-1),    // 0 = oldFinder, 1 = newFinder
  centFramework(-1)// 0 = AliCentrality, 1 = AliMultSelection 

{
  // default constructor
}

void DeDxEvent::Copy(TObject& object) const
{
  TObject::Copy(object);

  DeDxEvent* eventIn = (DeDxEvent*)(&object);
  if(!eventIn)
    return;

  eventIn->eventid       = eventid     ; 
  eventIn->run           = run         ; 
  eventIn->time          = time        ; 
  eventIn->cent          = cent        ;
  eventIn->mag           = mag         ; 
  eventIn->zvtx          = zvtx        ; 
  eventIn->zvtxMC        = zvtxMC      ; 
  eventIn->ptmax         = ptmax       ; 
  eventIn->ptmaxMC       = ptmaxMC     ; 
  eventIn->vtxstatus     = vtxstatus   ; 
  eventIn->trackmult     = trackmult   ; 
  eventIn->n             = n           ;
  eventIn->nTracks       = nTracks     ;
  eventIn->refMult       = refMult     ;
  eventIn->trackmultMC   = trackmultMC ; 
  eventIn->nMC           = nMC         ; 
  eventIn->process       = process     ; 
  eventIn->trig          = trig        ; 
  eventIn->triggerInt    = triggerInt  ; 
  eventIn->v0Finder      = v0Finder    ;
  eventIn->centFramework = centFramework;
}
