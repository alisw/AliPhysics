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

//Analysis
class AliAODEvent ;
class AliESDEvent ;
#include "AliAnaPartCorrBaseClass.h"

#ifdef __PHOSGEO__
	class AliPHOSGeoUtils;
#endif

class AliAnaPi0 : public AliAnaPartCorrBaseClass {
	
public: 
	
    AliAnaPi0() ; // default ctor
    AliAnaPi0(const char *name) ; // default ctor
    AliAnaPi0(const AliAnaPi0 & g) ; // cpy ctor
    AliAnaPi0 & operator = (const AliAnaPi0 & api0) ;//cpy assignment
    virtual ~AliAnaPi0() ;//virtual dtor
	
    TList * GetCreateOutputObjects(); 
	
    void Print(const Option_t * opt) const;

	void Init();
	void InitParameters();

    //void MakeAnalysisFillAOD() {;} //Not needed
    void MakeAnalysisFillHistograms();
	
//    void SetBadRunsList(){;} ;     //Set list of runs which can be used for this analysis
									 //To be defined in future.
    void SetEtalonHisto(TH3D * h);//Provide etalon of binning for histograms
	
    //Setters for parameters of event buffers
    void SetNCentrBin(Int_t n=5){fNCentrBin=n ;} //number of bins in centrality 
    void SetNZvertBin(Int_t n=5){fNZvertBin=n ;} //number of bins for vertex position
    void SetNRPBin(Int_t n=6)   {fNrpBin=n ;}    //number of bins in reaction plain
    void SetNMaxEvMix(Int_t n=20){fNmaxMixEv=n ;}//Maximal number of events for mixing
	
    //Setters for event selection
    void SetZvertexCut(Float_t zcut=40.){fZvtxCut=zcut ;} //cut on vertex position
	
	TString GetCalorimeter()   const {return fCalorimeter ; }
	void SetCalorimeter(TString det)    {fCalorimeter = det ; }
	
private:
	Bool_t IsBadRun(Int_t /*iRun*/) const {return kFALSE;} //Tests if this run bad according to private list
	
private:
	
	Int_t fNCentrBin ;		//number of bins in event container for centrality
	Int_t fNZvertBin ;		//number of bins in event container for vertex position
	Int_t fNrpBin ;			//number of bins in event container for reaction plain
	Int_t fNPID ;			//Number of possible PID combinations
	Int_t fNmaxMixEv ;		//Maximal number of events stored in buffer for mixing
	Float_t fZvtxCut ;		//Cut on vertex position
	TString fCalorimeter ;     //Select Calorimeter for IM
	
	TList ** fEventsList ;	//! containers for photons in stored events
			
	//Histograms
	
	TH3D * fhEtalon ; //Etalon histo, all distributions will have same binning as this one
	
	TH3D ** fhRe1 ;  //!REAL two-photon invariant mass distribution for different centralities and PID 
	TH3D ** fhMi1 ;  //!MIXED two-photon invariant mass distribution for different centralities and PID
	TH3D ** fhRe2 ;  //!REAL two-photon invariant mass distribution for different centralities and PID 
	TH3D ** fhMi2 ;  //!MIXED two-photon invariant mass distribution for different centralities and PID
	TH3D ** fhRe3 ;  //!REAL two-photon invariant mass distribution for different centralities and PID 
	TH3D ** fhMi3 ;  //!MIXED two-photon invariant mass distribution for different centralities and PID
	TH3D * fhEvents;  //!Number of events per centrality, RP, zbin
	
	//Acceptance
	TH1D * fhPrimPt ;    //! Spectrum of Primary 
	TH1D * fhPrimAccPt ; //! Spectrum of primary with accepted daughters 
	TH1D * fhPrimY ;     //! Rapidity distribution of primary particles
	TH1D * fhPrimAccY ;  //! Rapidity distribution of primary with accepted daughters
	TH1D * fhPrimPhi ;   //! Azimutal distribution of primary particles
	TH1D * fhPrimAccPhi; //! Azimutal distribution of primary with accepted daughters	
	
#ifdef __PHOSGEO__
       AliPHOSGeoUtils * fPHOSGeo ; //! PHOS geometry pointer
#endif	

	ClassDef(AliAnaPi0,3)
} ;


#endif //ALIANAPI0_H



