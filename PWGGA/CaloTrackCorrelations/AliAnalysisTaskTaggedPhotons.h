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
//*-- Dmitry Blau 
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"  

class AliStack ; 
class AliESDEvent ; 
//class AliAODEvent ; 
class TH1D ; 
class TH1I ; 
class TH2D ;
class TH2F ;
class TH3D ;
class AliPHOSGeoUtils;
class AliEMCALGeometry;
class AliAODPWG4Particle;

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

  void SetDebugLevel(Int_t level) { fDebug = level ; }
  void SetPHOS(Bool_t isPHOS=kTRUE){fPHOS=isPHOS ;} //Which calirimeter to analyze
  
  void SetPhotonId(Float_t threshold) { fPhotonId = threshold ; }
  void SetMinEnergyCut(Float_t threshold) { fMinEnergyCut = threshold; }

  //Parameterization of pi0 peak position
  void SetPi0MeanParameters(Float_t p0, Float_t p1, Float_t p2, Float_t p3)
      { fPi0MeanP0 = p0; fPi0MeanP1 = p1; fPi0MeanP2 = p2; fPi0MeanP3 = p3;}
  //Parameterization of pi0 peak width    
  void SetPi0SigmaParameters(Float_t p0, Float_t p1, Float_t p2){ fPi0SigmaP0 = p0; fPi0SigmaP1 = p1;fPi0SigmaP2 = p2; }

protected:
  Float_t GetPhotonId() const { return fPhotonId ; }
  Int_t   GetFiducialArea(Float_t * pos)const ; //what kind of fiducial area hit the photon
  Bool_t  IsSamePi0(const AliAODPWG4Particle *p1, const AliAODPWG4Particle *p2) const; //Check MC genealogy
  Bool_t  IsInPi0Band(Double_t m, Double_t pt)const; //Check if invariant mass is within pi0 peak
  Bool_t  TestDisp(Double_t l0, Double_t l1, Double_t e)const  ;
  Bool_t  TestTOF(Double_t /*t*/,Double_t /*en*/)const{return kTRUE;} 
  Bool_t  TestCharged(Double_t dr,Double_t en)const ;
  void    InitGeometry() ;  //read reotation matrixes from ESD/AOD

private:

  AliPHOSGeoUtils  *fPHOSgeom;   //!PHOS geometry
  AliEMCALGeometry *fEMCALgeom;  //!EMCAL geometry

  AliStack        *fStack ;      //Pointer to MC stack
  Bool_t           fPHOS ;       //Choose Calorimeter: PHOS/EMCAL

  // task parameters
  Float_t   fPhotonId ;          // threshold for photon identification (Bayesian)
  Float_t   fMinEnergyCut;       // min energy of partner photon
  Float_t   fPi0MeanP0;          // Parameterization of pi0 mass:
  Float_t   fPi0MeanP1;          // m_mean_pi0 = p[0] + p[1]*pt + p[2]*pt^2 + p[3]*pt^3
  Float_t   fPi0MeanP2;          // m_mean_pi0 = p[0] + p[1]*pt + p[2]*pt^2 + p[3]*pt^3
  Float_t   fPi0MeanP3;          // m_mean_pi0 = p[0] + p[1]*pt + p[2]*pt^2 + p[3]*pt^3

  Float_t   fPi0SigmaP0;        // sigma_m_pi0 = sqrt ( p0*p0/x + p1*p1 + p2*p2/x/x)
  Float_t   fPi0SigmaP1;      // sigma_m_pi0 = sqrt ( p0*p0/x + p1*p1 + p2*p2/x/x)
  Float_t   fPi0SigmaP2;      // sigma_m_pi0 = sqrt ( p0*p0/x + p1*p1 + p2*p2/x/x)

  //Fiducial area parameters
  Float_t fZmax ;               //Rectangular
  Float_t fZmin ;               //area
  Float_t fPhimax ;             //covered by
  Float_t fPhimin ;             //full calorimeter

  TList   * fOutputList ;        // output data list
  TList   * fEventList ;         //  event list for mixed InvMass

  // Histograms
  //Reconstructed spectra
  TH1D    * fhRecAll[4];               // Spectrum of all reconstructed particles
  TH1D    * fhRecAllArea1;             // Spectrum of rec particles in Fid. Area 1
  TH1D    * fhRecAllArea2;             // Spectrum of rec particles in Fid. Area 2
  TH1D    * fhRecAllArea3;             // Spectrum of rec particles in Fid. Area 3

  //Sort registered particles spectra according MC information
  TH1D    * fhRecPhoton;               // Spectrum of rec. with primary==22 and no PID criteria
  TH1D    * fhRecOther;                // Spectrum of rec. with primary!=22 and no PID criteria
  TH1D    * fhRecPhotonPID[4];         // Spectrum of rec. with primary==22 and different PID criteria
  TH1D    * fhRecOtherPID[4];          // Spectrum of rec. with primary!=22 and different PID criteria
  TH1D    * fhRecPhotPi0 ;             // Spectrum of rec. photons from pi0 decays
  TH1D    * fhRecPhotEta ;             // Spectrum of rec. photons from eta decays
  TH1D    * fhRecPhotOmega ;           // Spectrum of rec. photons from omega decays
  TH1D    * fhRecPhotEtapr ;           // Spectrum of rec. photons from eta prime decays
  TH1D    * fhRecPhotConv ;            // Spectrum of rec. photons from conversion
  TH1D    * fhRecPhotHadron ;          // Spectrum of rec. photons from hadron-matter interactions
  TH1D    * fhRecPhotDirect ;          // Spectrum of rec. photons direct or no primary
  TH1D    * fhRecPhotOther ;           // Spectrum of rec. photons from other hadron decays

  //MC tagging: reasons of partner loss etc.
  TH1D  * fhDecWMCPartner ;         //pi0 decay photon which partner should be registered according to MC
  TH1D  * fhDecWMissedPartnerNotPhoton ; //Decay photon with missed non-photon partner
  TH1D  * fhDecWMissedPartnerAll ;  //Decay photons with partner missed due to some reason (sum of below)
  TH1D  * fhDecWMissedPartnerEmin;  //Decay photons with partner missed due to low energy
  TH1D  * fhDecWMissedPartnerConv;  //Decay photons with partner missed due to conversion
  TH1D  * fhDecWMissedPartnerStack;  //Decay photons with partner not present in Stack
  TH1D  * fhDecWMissedPartnerGeom0; //Decay photons with partner missed due geometry
  TH1D  * fhDecWMissedPartnerGeom1; //Decay photons with partner missed due geometry Fid. area. 1
  TH1D  * fhDecWMissedPartnerGeom2; //Decay photons with partner missed due geometry Fid. area. 2
  TH1D  * fhDecWMissedPartnerGeom3; //Decay photons with partner missed due geometry Fid. area. 3

  //MC tagging: Decay partners spectra
  TH1D  * fhPartnerMCReg ;       //Spectrum of decay partners which should be registered (MC)
  TH1D  * fhPartnerMissedEmin ;  //Spectrum of decay partners lost due to Emin cut  
  TH1D  * fhPartnerMissedConv ;  //Spectrum of decay partners lost due to conversion
  TH1D  * fhPartnerMissedGeo  ;  //Spectrum of decay partners lost due to acceptance

  //Tagging
  TH1D  * fhTaggedAll ;      //Spectrum of all tagged photons
  TH1D  * fhTaggedArea1 ;    //Spectrum of all tagged photons Fid. area1
  TH1D  * fhTaggedArea2 ;    //Spectrum of all tagged photons Fid. area2
  TH1D  * fhTaggedArea3 ;    //Spectrum of all tagged photons Fid. area3
  TH1D  * fhTaggedPID[4] ;   //Spectrum of tagged photons for different PID criteria
  TH1D  * fhTaggedMult ;     //Spectrum of multiply tagged photons

  //Tagging: use MC information if available
  TH1D *  fhTaggedMCTrue ;     //Spectrum of correctly tagged pi0 decay photons
  TH1D *  fhMCMissedTagging ;  //Spectrum of pi0 decay photons missed tagging due to wrong pair mass
  TH1D *  fhMCFakeTagged ;     //Spectrum of photons wrongly tagged according to MC

  //Invariant mass distributions for fake corrections
  TH2D * fhInvMassReal[4] ;    //Two-photon inv. mass vs first photon pt
  TH2D * fhInvMassMixed[4] ;   //Two-photon inv. mass vs first photon pt
  TH2D * fhMCMissedTaggingMass ; //Inv mass of pairs missed tagging
  
  //Conversion and annihilation radius distributions
  TH1D * fhConversionRadius ;      // Radis of photon production (conversion)
  TH1D * fhInteractionRadius ;     // Radis of photon production (hadron interaction)

  TH1D    * fhEvents;                  //  number of processed Events
 
  ClassDef(AliAnalysisTaskTaggedPhotons, 1);   // a PHOS photon analysis task 
};
#endif // ALIANALYSISTASKTAGGEDPHOTONS_H
