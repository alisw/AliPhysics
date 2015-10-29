
#ifndef ALISPECTRABOTHPID_H
#define ALISPECTRABOTHPID_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliSpectraBothPID
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
class AliVTrack;
class AliAODMCParticle;
class TParticle;
class AliPIDResponse;  
class AliSpectraBothTrackCuts; 

#include "TNamed.h"
#include "TParticle.h"
#include "AliSpectraBothHistoManager.h" 

/*namespace AliSpectraNameSpaceBoth {

  enum BothPIDType_t
   {
       kNSigmaTPC,
       kNSigmaTOF,
       kNSigmaTPCTOF, // squared sum
   };



}*/

using namespace AliSpectraNameSpaceBoth;

class AliSpectraBothPID : public TNamed
{
public:
   enum BothPIDType_t {kNSigmaTPC,kNSigmaTOF,kNSigmacircleTPCTOF,kNSigmasquareTPCTOF,kNSigmaTPCorTOF};

  AliSpectraBothPID() ;
  AliSpectraBothPID(BothPIDType_t pidType);
  virtual  ~AliSpectraBothPID() {}

  void FillQAHistos(AliSpectraBothHistoManager * hman, AliVTrack * track, AliSpectraBothTrackCuts * trackCuts) ;
  void SetNSigmaCut(Float_t nsigma) { fNSigmaPID = nsigma; }
  void SetPIDtype(BothPIDType_t pidType){fPIDType=pidType;}
  void SetShiftTPC(Float_t shift){fshiftTPC=shift;}
  void SetShiftTOF(Float_t shift){fshiftTOF=shift;}
  Float_t GetNSigmaCut() {return fNSigmaPID; }

  Int_t  GetParticleSpecie(AliSpectraBothHistoManager * hman,AliVTrack * trk, AliSpectraBothTrackCuts * trackCuts, Bool_t* rec);
  Int_t GetParticleSpecie(AliAODMCParticle * trk);
  Int_t GetParticleSpecie(TParticle * trk);
  
  
  
  Long64_t Merge(TCollection* list);


private:

  BothPIDType_t fPIDType; // PID type
  Float_t fNSigmaPID; // number of sigma for PID cut
  AliPIDResponse   *fPIDResponse;     // ! PID response object
  Float_t fshiftTPC; // shift of the nsigma TPC
 Float_t fshiftTOF; // shift of the nsigma TPC

  AliSpectraBothPID(const AliSpectraBothPID&);
  AliSpectraBothPID& operator=(const AliSpectraBothPID&);

  ClassDef(AliSpectraBothPID, 3);

};
#endif

