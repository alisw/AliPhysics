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
class AliPIDCombined;  

#include "TNamed.h"

namespace AliHelperPIDNameSpace {
  
  enum PIDType_t
  {
    kNSigmaTPC = 0,
    kNSigmaTOF,
    kNSigmaTPCTOF,  // squared sum
    kBayes
  };
  
  const Int_t kNSigmaPIDType=kNSigmaTPCTOF;//number of Nsigma PID types
  
  enum AliHelperDetectorType_t
  {
    kITS = 0,
    kTPC,
    kTOF,
    kNDetectors
  };
  
  
  enum AliHelperParticleSpecies_t
  {
    kSpPion = 0,
    kSpKaon,
    kSpProton,
    kNSpecies,
    kSpUndefined=999
  }; // Particle species used in plotting
  
  
  enum AliHelperCharge_t
  {
    kChPos = 0,
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
  void SetfRemoveTracksT0Fill(Bool_t tof){fRemoveTracksT0Fill=tof;}//fRemoveTracksT0Fill
  Bool_t GetfRemoveTracksT0Fill(){return   fRemoveTracksT0Fill;}//fRemoveTracksT0Fill
  //Exclusive NSIgma
  void SetfUseExclusiveNSigma(Bool_t nsigEx){fUseExclusiveNSigma=nsigEx;}//fUseExclusiveNSigma
  Bool_t GetfUseExclusiveNSigma(){return   fUseExclusiveNSigma;}//fUseExclusiveNSigma
  //lower pt fot TOF PID
  Double_t GetPtTOFPID(){return   fPtTOFPID;}
  void SetfPtTOFPID(Double_t pttof){fPtTOFPID=pttof;}
  //set PID Combined
  void SetPIDCombined(AliPIDCombined *obj){fPIDCombined=obj;}
  //void SetPIDCombined(AliPIDCombined *obj){Printf("void AliHelperPID::SetPIDCombined(AliPIDCombined *obj) not implemented");}  //FIXME Left for backward compatibility, not the PIDCombined onject is created in the constructor as done in /ANALYSIS/AliAnalysisTaskPIDCombined.cxx (Jul 15th 2014)
  AliPIDCombined *GetPIDCombined(){return fPIDCombined;}
  //set cut on beyesian probability
  void SetBayesCut(Double_t cut){fBayesCut=cut;}
  Double_t GetBayesCut(){return fBayesCut;}
  
  //getters of the other data members
  TList * GetOutputList() {return fOutputList;}//get the TList with histos
  Double_t* GetNSigmas(AliHelperParticleSpecies_t species) {return fnsigmas[species];}//get nsigma[ipart][idet], calculated in CalculateNSigmas(trk)
  Bool_t* GetfHasDoubleCounting() {return fHasDoubleCounting;}//get fHasDoubleCounting[ipart], calculated in GetDoubleCounting(trk)
  //getter of histo "name" from fOutput
  TH2F* GetHistogram2D(const char * name);//return histogram "name" from fOutputList
  
  //PID functions
  // User should call ONLY the function GetParticleSpecies and set the PID strategy in the steering macro!
  Int_t GetParticleSpecies(AliVTrack * trk, Bool_t FIllQAHistos);//calculate the PID according to the slected method.
  Int_t GetParticleSpecies(AliVParticle * part);

  Int_t GetIDBayes(AliVTrack * trk, Bool_t FIllQAHistos);//calculate the PID according to bayesian PID
  UInt_t CalcPIDCombined(const AliVTrack *track,const AliPIDResponse *PIDResponse, Int_t detMask, Double_t* prob) const;
  void CalculateNSigmas(AliVTrack * trk, Bool_t FIllQAHistos);//Calcuate nsigma[ipart][idet], fill NSigma histos
  Int_t FindMinNSigma(AliVTrack * trk, Bool_t FIllQAHistos);//retun the minimum Nsigma
  Bool_t* GetDoubleCounting(AliVTrack * trk, Bool_t FIllQAHistos);//if a particle has double counting set fHasDoubleCounting[ipart]=kTRUE (only for the second or third identity)
  Bool_t* GetAllCompatibleIdentitiesNSigma(AliVTrack * trk, Bool_t FIllQAHistos);//All the identities are true
  Int_t GetMCParticleSpecie(AliVEvent* event, AliVTrack * trk, Bool_t FIllQAHistos);//calculate the PID according to MC truth
  void CheckTOF(AliVTrack * trk);//check the TOF matching and set fHasTOFPID
  Double_t TOFBetaCalc(AliVTrack *track) const;
  Double_t GetMass(AliHelperParticleSpecies_t id) const;
  Long64_t Merge(TCollection* list);
  
 private:
  
  Bool_t fisMC;
  PIDType_t fPIDType; // PID type
  Double_t fNSigmaPID; // number of sigma for PID cut
  Double_t fBayesCut; // Cut on Bayesian probability
  AliPIDResponse   *fPIDResponse;     //! PID response object
  AliPIDCombined   *fPIDCombined;     // PIDCombined
  TList     *fOutputList;  // List Histo's
  Double_t fnsigmas[kNSpecies][kNSigmaPIDType+1]; //nsigma values
  Bool_t fHasDoubleCounting[kNSpecies];//array with compatible identities
  Bool_t fRequestTOFPID;//if true returns kSpUndefined if the TOF signal is missing
  Bool_t fRemoveTracksT0Fill;//if true remove tracks for which only StartTime from To-Fill is available (worst resolution)
  Bool_t fUseExclusiveNSigma;//if true returns the identity only if no double counting
  Double_t fPtTOFPID; //lower pt bound for the TOF pid
  Bool_t fHasTOFPID;
  
  AliHelperPID(const AliHelperPID&);
  AliHelperPID& operator=(const AliHelperPID&);
  
  ClassDef(AliHelperPID, 8);
  
};
#endif

 
