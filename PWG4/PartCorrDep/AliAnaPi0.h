#ifndef ALIANAPI0_H
#define ALIANAPI0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Class to fill two-photon invariant mass hisograms 
// to be used to extract pi0 raw yield.
//
//-- Author: Dmitri Peressounko (RRC "KI")
//-- Adapted to PartCorr frame by Lamia Benhabib (SUBATECH)
//-- and Gustavo Conesa (INFN-Frascati)

//Root
class TList;
class TH3D ;
class TH2D ;

//Analysis
class AliAODEvent ;
class AliESDEvent ;
#include "AliAnaPartCorrBaseClass.h"

class AliAnaPi0 : public AliAnaPartCorrBaseClass {
  
public:   
  AliAnaPi0() ; // default ctor
  AliAnaPi0(const char *name) ; // default ctor
  AliAnaPi0(const AliAnaPi0 & g) ; // cpy ctor
  virtual ~AliAnaPi0() ;//virtual dtor

private:
  AliAnaPi0 & operator = (const AliAnaPi0 & api0) ;//cpy assignment
  
public:
  TList * GetCreateOutputObjects(); 
  
  void Print(const Option_t * opt) const;
  
  //void Init();
  void InitParameters();
  
  //void MakeAnalysisFillAOD() {;} //Not needed
  void MakeAnalysisFillHistograms();
  
  //    void SetBadRunsList(){;} ;     //Set list of runs which can be used for this analysis
  //To be defined in future.
  
  //void SetEtalonHisto(TH3D * h);//Provide etalon of binning for histograms
  
  //Setters for parameters of event buffers
  void SetNCentrBin(Int_t n=5) {fNCentrBin=n ;} //number of bins in centrality 
  void SetNZvertBin(Int_t n=5) {fNZvertBin=n ;} //number of bins for vertex position
  void SetNRPBin(Int_t n=6)    {fNrpBin=n ;}    //number of bins in reaction plain
  void SetNMaxEvMix(Int_t n=20){fNmaxMixEv=n ;} //Maximal number of events for mixing
  
  //Setters for event selection
  void SetZvertexCut(Float_t zcut=40.){fZvtxCut=zcut ;} //cut on vertex position
  
  TString GetCalorimeter()   const {return fCalorimeter ; }
  void SetCalorimeter(TString det)    {fCalorimeter = det ; }
  	
  void Terminate(TList* outputList);
  void ReadHistograms(TList * outputList); //Fill histograms with histograms in ouput list, needed in Terminate.
	
  void SetNumberOfModules(Int_t nmod) {fNModules = nmod;}
	
  Int_t GetNPID()   const   {return fNPID ; }
  void  SetNPID(Int_t n)    {fNPID = n ; }
	
  void SwitchOnAngleSelection()    {fUseAngleCut = kTRUE ; }
  void SwitchOffAngleSelection()   {fUseAngleCut = kFALSE ; }

  private:
  Bool_t IsBadRun(Int_t /*iRun*/) const {return kFALSE;} //Tests if this run bad according to private list
  
  private:
  Int_t    fNCentrBin ;	  // Number of bins in event container for centrality
  Int_t    fNZvertBin ;	  // Number of bins in event container for vertex position
  Int_t    fNrpBin ;	  // Number of bins in event container for reaction plain
  Int_t    fNPID ;		  // Number of possible PID combinations
  Int_t    fNmaxMixEv ;	  // Maximal number of events stored in buffer for mixing
  Float_t  fZvtxCut ;	  // Cut on vertex position
  TString  fCalorimeter ; // Select Calorimeter for IM
  Int_t    fNModules ;    // Number of EMCAL/PHOS modules, set as many histogras as modules 
  Bool_t   fUseAngleCut ; // Select pairs depending on their opening angle
  TList ** fEventsList ;  //! Containers for photons in stored events
  
  //Histograms
  
  //TH3D * fhEtalon ; //Etalon histo, all distributions will have same binning as this one
  
  TH3D ** fhReMod ;  //!REAL two-photon invariant mass distribution for different calorimeter modules.
	
  TH3D ** fhRe1 ;  //!REAL two-photon invariant mass distribution for different centralities and PID 
  TH3D ** fhMi1 ;  //!MIXED two-photon invariant mass distribution for different centralities and PID
  TH3D ** fhRe2 ;  //!REAL two-photon invariant mass distribution for different centralities and PID 
  TH3D ** fhMi2 ;  //!MIXED two-photon invariant mass distribution for different centralities and PID
  TH3D ** fhRe3 ;  //!REAL two-photon invariant mass distribution for different centralities and PID 
  TH3D ** fhMi3 ;  //!MIXED two-photon invariant mass distribution for different centralities and PID
  TH3D * fhEvents; //!Number of events per centrality, RP, zbin

  TH2D * fhRealOpeningAngle ;    //! Opening angle of pair versus pair energy
  TH2D * fhRealCosOpeningAngle ; //! Cosinus of opening angle of pair version pair energy

  //Acceptance
  TH1D * fhPrimPt ;    //! Spectrum of Primary 
  TH1D * fhPrimAccPt ; //! Spectrum of primary with accepted daughters 
  TH1D * fhPrimY ;     //! Rapidity distribution of primary particles
  TH1D * fhPrimAccY ;  //! Rapidity distribution of primary with accepted daughters
  TH1D * fhPrimPhi ;   //! Azimutal distribution of primary particles
  TH1D * fhPrimAccPhi; //! Azimutal distribution of primary with accepted daughters	
  TH2D * fhPrimOpeningAngle ;    //! Opening angle of pair versus pair energy, primaries
  TH2D * fhPrimCosOpeningAngle ; //! Cosinus of opening angle of pair version pair energy, primaries
	
  ClassDef(AliAnaPi0,8)
} ;


#endif //ALIANAPI0_H



