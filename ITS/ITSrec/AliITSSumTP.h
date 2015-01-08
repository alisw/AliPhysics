#ifndef ALIITSSUMTP_H
#define ALIITSSUMTP_H

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for ITS trackpoints summary + some aux. info  )         //
// Author: Ruben Shahoian                                        //
//                                                               //
///////////////////////////////////////////////////////////////////

/* $Id$ */

class AliTrackPointArray;
#include <TObjArray.h>
#include "AliESDVertex.h"

class AliITSSumTP : public TObject
{
 public:
  enum {kCrvTPC,kCrvTPCErr,kCrvGlo,kCrvGloErr,kNVarPerTrack};
  enum {kTPCInX,kTPCInY,kTPCInZ,kTPCSpare};
  //
  AliITSSumTP() : fTracks(0),fVertex(),fNVars(0),fCrvVars(0),fTPCVars(0) {fTracks.SetOwner(kTRUE);}
  AliITSSumTP(const AliITSSumTP& src);
  virtual ~AliITSSumTP() {Reset();}
  AliITSSumTP &operator=(const  AliITSSumTP& src);
  virtual void Print(Option_t* opt="") const;
  //
  AliESDVertex& GetVertex()             const {return (AliESDVertex&)fVertex;}
  AliTrackPointArray* GetTrack(Int_t i) const {return (AliTrackPointArray*)fTracks[i];}
  TObjArray&   GetTracks()              const {return (TObjArray&) fTracks;}
  Int_t        GetNTracks()             const {return fNVars/kNVarPerTrack;}
  Double_t     GetCrvTPC(Int_t i)       const {return fCrvVars[i*kNVarPerTrack+kCrvTPC];}
  Double_t     GetCrvTPCErr(Int_t i)    const {return fCrvVars[i*kNVarPerTrack+kCrvTPCErr];}
  Double_t     GetCrvGlo(Int_t i)       const {return fCrvVars[i*kNVarPerTrack+kCrvGlo];}
  Double_t     GetCrvGloErr(Int_t i)    const {return fCrvVars[i*kNVarPerTrack+kCrvGloErr];}
  void         GetTPCInnerXYZ(int i, double* xyz) const;
  void         SetTPCInnerXYZ(int i, const double* xyz);
  //
  void         Reset();
  void         BookNTracks(Int_t n);
  void         SetCrvTPC(Int_t i, Double_t v)       {fCrvVars[i*kNVarPerTrack+kCrvTPC] = v;}
  void         SetCrvTPCErr(Int_t i, Double_t v)    {fCrvVars[i*kNVarPerTrack+kCrvTPCErr] = v;}
  void         SetCrvGlo(Int_t i, Double_t v)       {fCrvVars[i*kNVarPerTrack+kCrvGlo] = v;}
  void         SetCrvGloErr(Int_t i, Double_t v)    {fCrvVars[i*kNVarPerTrack+kCrvGloErr] = v;}
  void         AddTrack(AliTrackPointArray* trc)    {fTracks.AddLast((TObject*)trc);}
  void         SetVertex(const AliESDVertex* vtx)   {fVertex = *vtx;}
  //
 protected:
  //
  TObjArray    fTracks;               // TrackPoints
  AliESDVertex fVertex;               // ESD Vertex
  Int_t        fNVars;                // Ntracks*kNVarPerTrack
  Double32_t*  fCrvVars;               //[fNVars];
  Double32_t*  fTPCVars;               //[fNVars];


  ClassDef(AliITSSumTP,2)
};


//---------------------------------------------------------------
inline void AliITSSumTP::GetTPCInnerXYZ(int i, double* xyz) const
{
  // get tpc inner coordinates
  if (!fTPCVars) {xyz[0]=xyz[1]=xyz[2]=0; return;}  // for BWD compatibility
  xyz[0] = fTPCVars[i*kNVarPerTrack+kTPCInX];
  xyz[1] = fTPCVars[i*kNVarPerTrack+kTPCInY];
  xyz[2] = fTPCVars[i*kNVarPerTrack+kTPCInZ];
}

//---------------------------------------------------------------
inline void AliITSSumTP::SetTPCInnerXYZ(int i, const double* xyz)
{
  // Set tpc inner coordinates
  if (!fTPCVars) return; // for BWD compatibility
  fTPCVars[i*kNVarPerTrack+kTPCInX] = xyz[0];
  fTPCVars[i*kNVarPerTrack+kTPCInY] = xyz[1];
  fTPCVars[i*kNVarPerTrack+kTPCInZ] = xyz[2];
}


#endif

