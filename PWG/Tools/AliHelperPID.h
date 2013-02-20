#ifndef ALIHELPERPID_H
#define ALIHELPERPID_H

class AliAODEvent;
class AliVEvent;
class TH1F;
class TH2F;
class TList;
class AliVTrack;
class AliVParticle;
class AliStack;
class TParticle;
class AliPIDResponse;  

#include "TNamed.h"

namespace AliHelperPIDNameSpace {
  
  enum PIDType_t
  {
    kNSigmaTPC,
    kNSigmaTOF,
    kNSigmaTPCTOF,  // squared sum
    kNSigmaPIDType=kNSigmaTPCTOF
  };
  
  
  enum AliHelperDetectorType_t
  {
    kTPC,
    kTOF,
    kNDetectors
  };
  
  
  enum AliHelperParticleSpecies_t
  {
    kSpPion,
    kSpKaon,
    kSpProton,
    kNSpecies,
    kSpUndefined=999
  }; // Particle species used in plotting
  
  
  enum AliHelperCharge_t
  {
    kChPos,
    kChNeg,
    kNCharge
  };
}

using namespace AliHelperPIDNameSpace;

class AliHelperPID : public TNamed
{
 public:
  
  AliHelperPID();
  virtual  ~AliHelperPID() {}
  
  //MC or data
  Bool_t GetisMC(){return   fisMC;}
  void SetisMC(Bool_t mc){fisMC=mc;}
  //PID Type
  void SetPIDType(PIDType_t PIDType) { fPIDType = PIDType; }
  PIDType_t GetPIDType() {return fPIDType; }
  //NSigma cut
  void SetNSigmaCut(Double_t nsigma) { fNSigmaPID = nsigma; }
  Double_t GetNSigmaCut() {return fNSigmaPID; }
  //TOF PID
  void SetfRequestTOFPID(Bool_t tof){fRequestTOFPID=tof;}//fRequestTOFPID
  Bool_t GetfRequestTOFPID(){return   fRequestTOFPID;}//fRequestTOFPID
  //lower pt fot TOF PID
  Double_t SetPtTOFPID(){return   fPtTOFPID;}
  void SetfPtTOFPID(Double_t pttof){fPtTOFPID=pttof;}
  //getters of the other data members
  TList * GetOutputList() {return fOutputList;}//get the TList with histos
  Double_t* GetNSigmas() {return *fnsigmas;}//get nsigma[ipart][idet], calculated in CalculateNSigmas(trk)
  Bool_t* GetfHasDoubleCounting() {return fHasDoubleCounting;}//get fHasDoubleCounting[ipart], calculated in GetDoubleCounting(trk)
  //getter of histo "name" from fOutput
  TH2F* GetHistogram2D(const char * name);//return histogram "name" from fOutputList
  
  //PID functions
  Int_t GetParticleSpecies(AliVTrack * trk, Bool_t FIllQAHistos);//calculate the PID according to the minimum sigma
  void CalculateNSigmas(AliVTrack * trk, Bool_t FIllQAHistos);//Calcuate nsigma[ipart][idet], fill NSigma histos
  Int_t FindMinNSigma(AliVTrack * trk, Bool_t FIllQAHistos);//retun the minimum Nsigma
  Bool_t* GetDoubleCounting(AliVTrack * trk, Bool_t FIllQAHistos);//if a particle has double counting set fHasDoubleCounting[ipart]=kTRUE
  Int_t GetMCParticleSpecie(AliVEvent* event, AliVTrack * trk, Bool_t FIllQAHistos);//calculate the PID according to MC truth
  void CheckTOF(AliVTrack * trk);//check the TOF matching and set fHasTOFPID
  Long64_t Merge(TCollection* list);
  
 private:
  
  Bool_t fisMC;
  PIDType_t fPIDType; // PID type
  Double_t fNSigmaPID; // number of sigma for PID cut
  AliPIDResponse   *fPIDResponse;     // ! PID response object
  TList     *fOutputList;  // List Histo's
  Double_t fnsigmas[kNSpecies][kNSigmaPIDType+1]; //nsigma values
  Bool_t fHasDoubleCounting[kNSpecies];//array with compatible identities
  Bool_t fRequestTOFPID;//if true returns kSpUndefined if the TOF signal is missing
  Double_t fPtTOFPID; //lower pt bound for the TOF pid
  Bool_t fHasTOFPID;
  
  
  AliHelperPID(const AliHelperPID&);
  AliHelperPID& operator=(const AliHelperPID&);
  
  ClassDef(AliHelperPID, 1);
  
};
#endif

 
