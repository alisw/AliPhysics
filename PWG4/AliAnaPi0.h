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

class TH3D ;
class AliAODEvent ;
class AliESDEvent ;

#include "AliAnalysisTaskSE.h"
 
class AliAnaPi0 : public AliAnalysisTaskSE {
       
  public: 
       
    AliAnaPi0() ; // default ctor
    AliAnaPi0(const char *name) ; // default ctor
    AliAnaPi0(const AliAnaPi0 & g) ; // cpy ctor
    AliAnaPi0 & operator = (const AliAnaPi0 & api0) ;//cpy assignment
    virtual ~AliAnaPi0() ;//virtual dtor
              
    //obligatory methods
    virtual void UserCreateOutputObjects();
    virtual void Init() ;
    virtual void LocalInit() { Init() ; }
    virtual void UserExec(Option_t * opt = "") ;
    virtual void Terminate(Option_t * opt = "") ;
 
    void InitParameters();
    void Print(const Option_t * opt) const;

    void SetBadRunsList(){} ;     //Set list of runs which can be used for this analysis
    void SetEtalonHisto(TH3D * h);//Provide etalon of binning for histograms

    //Setters for parameters of event buffers
    void SetNCentrBin(Int_t n=5){fNCentrBin=n ;} //number of bins in centrality 
    void SetNZvertBin(Int_t n=5){fNZvertBin=n ;} //number of bins for vertex position
    void SetNRPBin(Int_t n=6)   {fNrpBin=n ;}    //number of bins in reaction plain
    void SetNMaxEvMix(Int_t n=20){fNmaxMixEv=n ;}//Maximal number of events for mixing
 
    //Setters for event selection
    void SetZvertexCut(Float_t zcut=40.){fZvtxCut=zcut ;} //cut on vertex position
    void SetPtMin(Float_t pt = 0.2) {fPtMin=pt ;} //Cut on minimal pt of photon

 private:
    void FillHistograms() ; //fill invariant mass distributions
    Bool_t FillFromAOD(AliAODEvent * aod) ;
    Bool_t FillFromESD(AliESDEvent * esd) ;

    Bool_t IsBadRun(Int_t /*iRun*/) {return kFALSE;} //Tests if this run bad according to private list

 private:
 
       Int_t fNCentrBin ;   //number of bins in event container for centrality
       Int_t fNZvertBin ;   //number of bins in event container for vertex position
       Int_t fNrpBin ;      //number of bins in event container for reaction plain
       Int_t fNPID ;        //Number of possible PID combinations
       Int_t fNmaxMixEv ;   //Maximal number of events stored in buffer for mixing

       Int_t fCurCentrBin ; //! bin for current event (centrality)
       Int_t fCurZvertBin ; //! bin for current event (z position)
       Int_t fCurRPBin ;    //! bin for current event (reaction plain)

       Float_t fPtMin ;     //Cut on minimum photon Pt
       Float_t fZvtxCut ;   //Cut on vertex position
       Double_t fVert[3] ;   //Currecnt vertex position
       TList ** fEventsList ;        //! containers for photons in stored events
       TClonesArray * fCurrentEvent ; //!Container  with latest event

       TList * fOutputList ; //list of output objects

       //Histograms
       TH3D * fhEtalon ; //!Etalon histo, all distributions will have same binning as this one

       TH3D ** fhRe1 ;  //!REAL two-photon invariant mass distribution for different centralities and PID 
       TH3D ** fhMi1 ;  //!MIXED two-photon invariant mass distribution for different centralities and PID
       TH3D ** fhRe2 ;  //!REAL two-photon invariant mass distribution for different centralities and PID 
       TH3D ** fhMi2 ;  //!MIXED two-photon invariant mass distribution for different centralities and PID
       TH3D ** fhRe3 ;  //!REAL two-photon invariant mass distribution for different centralities and PID 
       TH3D ** fhMi3 ;  //!MIXED two-photon invariant mass distribution for different centralities and PID
       TH3D * fhEvents;  //!Number of events per centrality, RP, zbin

       //Calo Cells
       ClassDef(AliAnaPi0,1)
 } ;


#endif //ALIANAPI0_H



