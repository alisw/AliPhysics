#ifndef ALIANALYSISTASKTAGGEDPHOTONS_H
#define ALIANALYSISTASKTAGGEDPHOTONS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis for PHOS Tagged Photons 
// marks photons making pi0 with any other photon
// and calculates necessary corrections for fake pairs and
// decay partners escaped acceptance. If MC info is present 
// fills set of controll histograms.
//
//*-- Dmitry Peresunko
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"  

class AliStack ; 
class AliAODEvent ; 
class THashList ; 
class TH2I ;
class AliPHOSGeometry;
class AliAODPWG4Particle;
class AliVCluster ;
class AliTriggerAnalysis ;
class TParticle ;

class AliAnalysisTaskTaggedPhotons : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskTaggedPhotons() ;
  AliAnalysisTaskTaggedPhotons(const char *name) ;
  AliAnalysisTaskTaggedPhotons(const AliAnalysisTaskTaggedPhotons& ap) ;   
  AliAnalysisTaskTaggedPhotons& operator = (const AliAnalysisTaskTaggedPhotons& ap) ;
  virtual ~AliAnalysisTaskTaggedPhotons() ;
   
  virtual void UserCreateOutputObjects(); 
  virtual void Init() ; 
  virtual void LocalInit() { Init() ; }
  virtual void UserExec(Option_t * opt = "") ;
  virtual void Terminate(Option_t * opt = "") ;

  void SetPHOSBadMap(Int_t mod,TH2I * h)
  {
    if(fPHOSBadMap[mod]) delete fPHOSBadMap[mod] ;
    fPHOSBadMap[mod]=new TH2I(*h) ;
    printf("Set %s \n",fPHOSBadMap[mod]->GetName());
  }

protected:
  void    FillMCHistos() ;
  void    FillTaggingHistos() ;
  Int_t   GetFiducialArea(const Float_t * pos)const ; //what kind of fiducial area hit the photon
  Bool_t  IsSamePi0(const AliAODPWG4Particle *p1, const AliAODPWG4Particle *p2) const; //Check MC genealogy
  Bool_t  IsInPi0Band(Double_t m, Double_t pt,Int_t type)const; //Check if invariant mass is within pi0 peak
  Bool_t  TestDisp(Double_t l0, Double_t l1, Double_t e)const  ;
  Bool_t  TestTOF(Double_t /*t*/,Double_t /*en*/)const{return kTRUE;} 
  Bool_t  TestCharged(Double_t dr,Double_t en)const ;
  void    InitGeometry() ;  //read reotation matrixes from AOD/AOD
  Int_t   EvalIsolation(TLorentzVector * ph) ;
  Bool_t  TestLambda(Double_t pt,Double_t l1,Double_t l2) ;
  Bool_t  TestPID(Int_t iPID, AliAODPWG4Particle* part) ;
  Double_t PrimaryParticleWeight(TParticle * particle) ;
  Int_t   FindPrimary(AliVCluster*, Bool_t&);
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void FillPIDHistograms(const char * name, const AliAODPWG4Particle * p) const ;
  void FillPIDHistograms(const char * name, const AliAODPWG4Particle * p ,Double_t y) const ;

private:

  AliPHOSGeometry  *fPHOSgeom;   //!PHOS geometry
  THashList *   fOutputContainer ;   //! List of output histograms
  AliStack        *fStack ;      //!Pointer to MC stack
  TClonesArray * fTrackEvent ;   //!List of tracks in the event
  TClonesArray * fPHOSEvent ;    //!List of tracks in the event
  TList   * fPHOSEvents[1][5] ; //!Previous events for mixing
  TList   * fCurrentMixedList;   //! list of previous evetns for given centrality
  AliTriggerAnalysis * fTriggerAnalysis ; //!
 
  //Fiducial area parameters
  Float_t fZmax ;               //Rectangular
  Float_t fZmin ;               //area
  Float_t fPhimax ;             //covered by
  Float_t fPhimin ;             //full calorimeter

  //
  Double_t fCentrality;
  Int_t fCentBin ;
  
  // Histograms
  TH2I * fPHOSBadMap[6] ; 
  
  ClassDef(AliAnalysisTaskTaggedPhotons, 2);   // a PHOS photon analysis task 
};
#endif // ALIANALYSISTASKTAGGEDPHOTONS_H
