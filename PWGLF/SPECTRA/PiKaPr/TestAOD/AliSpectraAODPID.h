
#ifndef ALISPECTRAAODPID_H
#define ALISPECTRAAODPID_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliSpectraAODPID
//
//
//
//
// Authors: Michele Floris, CERN, Leonardo Milano, Torino
//-------------------------------------------------------------------------

class AliAODEvent;
class TH1F;
class TH2F;
class TList;
class AliAODTrack;
class AliAODMCParticle;
class AliPIDResponse;  
class AliSpectraAODTrackCuts; 

#include "TNamed.h"
#include "AliSpectraAODHistoManager.h" 

namespace AliSpectraNameSpace {

  enum AODPIDType_t
   {
       kNSigmaTPC,
       kNSigmaTOF,
       kNSigmaTPCTOF, // squared sum
   };



}

using namespace AliSpectraNameSpace;

class AliSpectraAODPID : public TNamed
{
public:
  AliSpectraAODPID() ;
  AliSpectraAODPID(AODPIDType_t pidType);
  virtual  ~AliSpectraAODPID() {}

  void FillQAHistos(AliSpectraAODHistoManager * hman, AliAODTrack * track, AliSpectraAODTrackCuts * trackCuts) ;
  void SetNSigmaCut(Float_t nsigma) { fNSigmaPID = nsigma; }

  Int_t GetParticleSpecie(AliSpectraAODHistoManager * hman,AliAODTrack      * trk, AliSpectraAODTrackCuts * trackCuts);
  Int_t GetParticleSpecie(AliAODMCParticle * trk);
  
  Long64_t Merge(TCollection* list);


private:

  AODPIDType_t fPIDType; // PID type
  Float_t fNSigmaPID; // number of sigma for PID cut
  AliPIDResponse   *fPIDResponse;     // ! PID response object


  AliSpectraAODPID(const AliSpectraAODPID&);
  AliSpectraAODPID& operator=(const AliSpectraAODPID&);

  ClassDef(AliSpectraAODPID, 1);

};
#endif

