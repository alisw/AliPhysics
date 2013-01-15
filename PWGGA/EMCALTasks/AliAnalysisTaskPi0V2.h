/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskPi0V2.h 45956 2010-12-10 12:55:37Z agheata $ */
/* AliAnalysisTaskPi0V2.h
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 */
#ifndef ALIANALYSISTASKPI0V2_H
#define ALIANALYSISTASKPI0V2_H

class TH1F;
class TH2F;
class TH3F;
class TList;
class AliESDCaloCluster;
class AliVCluster;
class AliESDtrackCuts;
class AliESDEvent;
class THnSparse;
class TClonesArray;
class TString;
class AliCaloPID;
class AliCalorimeterUtils;
class AliCaloTrackReader;


#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskPi0V2 : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskPi0V2();
    AliAnalysisTaskPi0V2(const char *name);
    virtual ~AliAnalysisTaskPi0V2();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    Double_t		GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax) const;
    Double_t		GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const;
    Bool_t		IsWithinFiducialVolume(Short_t id) const;
    Bool_t		IsGoodCluster(const AliVCluster *c) const;
    Bool_t		IsGoodClusterV1(const AliVCluster *c) const;
    Bool_t		IsGoodPion(const TLorentzVector& p1, const TLorentzVector& p2) const;
    void		FillPion(const TLorentzVector& p1, const TLorentzVector& p2, Double_t EPV0r, Double_t EPV0A, Double_t EPV0C, Double_t EPTPC);
    void		FillCluster(const TLorentzVector& p1, Double_t EPV0r, Double_t EPV0A, Double_t EPV0C, Double_t EPTPC, AliVCluster *c, Double_t *vertex);
    void 		GetMom(TLorentzVector& p, const AliVCluster *c, Double_t *vertex);		
    void		SetEventMethod(Double_t e )	{ fEvtSelect  =e ;}
    void		SetVtxCut(Double_t v )	        { fVtxCut     =v ;}
    void		SetClusNcell(Double_t c )	{ fNcellCut   =c ;}
    void		SetClusE(Double_t e )	        { fECut       =e ;}
    void		SetClusEta(Double_t e )	        { fEtaCut     =e ;}
    void		SetClusM02(Double_t m )	        { fM02Cut     =m ;}
    void		SetDrCut(Double_t m )	        { fDrCut      =m ;}
    void		SetPi0Asy(Double_t a )	        { fPi0AsyCut  =a ;}
    void                SetTracksName(const char *n)    { fTracksName =n ;}
    void                SetTrigClass(const char *n)     { fTrigClass  =n ;} 
    void                SetV1ClusName(const char *n)    { fV1ClusName =n ;} 
    void                SetV2ClusName(const char *n)    { fV2ClusName =n ;} 
    Int_t		ConvertToInternalRunNumber(Int_t n);
    void		FillEPQA();
    void		SetIsV1Clus(Bool_t e)		{ isV1Clus   =e  ;}
    

    
 private:
    TList           		*fOutput;	        //! Output list
    AliESDEvent			*fESD;		        //!ESD object
    AliCaloPID                  *fCaloPID;              //!
    AliCalorimeterUtils         *fCaloUtils;            //!
    AliCaloTrackReader          *fReader;               //!
    TString                     fTracksName;	        // name of track collection
    TString                     fV1ClusName;	        // name of V1 Clus collection
    TString                     fV2ClusName;	        // name of V1 Clus collection
    TString                     fTrigClass;	        // trigger class name for event selection
    TClonesArray                *fTracks;		//! pico tracks specific for Skim ESD
    TClonesArray                *fV1Clus;		//! Cluster Array for V1
    TClonesArray                *fV2Clus;		//! Cluster Array for V2
    Int_t			fRunNumber;		//! Run numbers
    Double_t 			fEvtSelect;	  	// 1 = MB+Semi+Central, 2 = MB+Semi, 3 = MB;
    Double_t			fVtxCut;		// vertex cut
    Double_t			fNcellCut;	        // N cells Cut
    Double_t			fECut;			// Cluster E cut
    Double_t			fEtaCut;		// Cluster Eta Cut
    Double_t			fM02Cut;		// Cluster long axis cut
    Double_t			fDrCut;		// Cluster long axis cut
    Bool_t			fPi0AsyCut;		// pion Asymetry cut 0=off 1=on
    Bool_t			isV1Clus;		// pion Asymetry cut 0=off 1=on
    Double_t			fCentrality;	  	//! Centrality
    Double_t			fEPTPC;			//! Evt plane TPC
    Double_t			fEPTPCreso;		//! resolution of TPC method
    Double_t			fEPV0;	    		//! EP V0
    Double_t			fEPV0A;			//! EP V0A
    Double_t			fEPV0C;	  		//! EP V0C
    Double_t			fEPV0Ar;		//! EP V0A reduced
    Double_t			fEPV0Cr;		//! EP V0C reduced
    Double_t			fEPV0r;	  		//! EP V0 reduced
    Double_t			fEPV0AR4;		//! EP V0A ring4 only
    Double_t			fEPV0AR5;		//! EP V0A ring5 only
    Double_t			fEPV0AR6;		//! EP V0A ring6 only
    Double_t			fEPV0AR7;		//! EP V0A ring7 only
    Double_t			fEPV0CR0;		//! EP V0C ring0 only	
    Double_t			fEPV0CR1;		//! EP V0C ring1 only	
    Double_t			fEPV0CR2;		//! EP V0C ring2 only	
    Double_t			fEPV0CR3;		//! EP V0C ring3 only	

    TH1F			*hEvtCount;		//!
    TH1F			*hAllcentV0;		//!
    TH1F			*hAllcentV0r;		//!
    TH1F			*hAllcentV0A;	  	//!
    TH1F			*hAllcentV0C;	    	//!
    TH1F			*hAllcentTPC;	        //!
  
    TH2F			*h2DcosV0r;		//! QA for cos(Phi) V0r vs Run NUmber
    TH2F			*h2DsinV0r;		//! QA for cos(Phi) V0r vs Run NUmber
    TH2F			*h2DcosV0A;		//!
    TH2F			*h2DsinV0A;		//!
    TH2F			*h2DcosV0C;		//!
    TH2F			*h2DsinV0C;		//!
    TH2F			*h2DcosTPC;		//!
    TH2F			*h2DsinTPC;		//!

    TH2F			*hEPTPC;		//! 2-D histo EPTPC  vs cent
    TH2F			*hresoTPC;		//! 2-D histo TPC resolution vs cent
    TH2F			*hEPV0;			//! 2-D histo EPV0   vs cent
    TH2F			*hEPV0A;		//! 2-D histo EPV0A  vs cent
    TH2F			*hEPV0C;		//! 2-D histo EPV0C  vs cent
    TH2F			*hEPV0Ar;		//! 2-D histo EPV0Ar vs cent
    TH2F			*hEPV0Cr;		//! 2-D histo EPV0Cr vs cent
    TH2F			*hEPV0r;		//! 2-D histo EPV0r  vs cent
    TH2F			*hEPV0AR4;		//! 2-D histo EPV0AR4 vs cent
    TH2F			*hEPV0AR7;		//! 2-D histo EPV0AR7 vs cent
    TH2F			*hEPV0CR0;		//! 2-D histo EPV0AR0 vs cent
    TH2F			*hEPV0CR3;		//! 2-D histo EPV0AR3 vs cent

    TH2F			*hdifV0Ar_V0Cr;		//! 2-D histo diff V0Ar, V0Cr vs cent
    TH2F			*hdifV0A_V0CR0;		//! 2-D histo diff V0A, V0CR0 vs cent
    TH2F			*hdifV0A_V0CR3;		//! 2-D histo diff V0A, V0CR3 vs cent
    TH2F			*hdifV0ACR0_V0CR3;	//! 2-D histo diff V0CR0, V0CR3 vs cent
    TH2F			*hdifV0C_V0AR4;		//! 2-D histo diff V0C, V0AR4 vs cent
    TH2F			*hdifV0C_V0AR7;		//! 2-D histo diff V0C, V0AR7 vs cent
    TH2F			*hdifV0AR4_V0AR7;	//! 2-D histo diff V0AR7, V0AR4 vs cent

    TH2F			*hdifV0A_V0C;		//! 2-D histo diff V0A - V0C
    TH2F			*hdifV0A_TPC;		//! 2-D histo diff V0A - TPC
    TH2F			*hdifV0C_TPC;		//! 2-D histo diff V0C - TPC
    TH2F			*hdifV0C_V0A;		//! 2-D histo diff V0C - V0A

    TH2F			*hM02vsPtA;		//! 2-D histo clus M02 vs Pt before cut
    TH2F			*hM02vsPtB;		//! 2-D histo clus M02 vs Pt after cut
    TH2F			*hClusDxDZA;		//! 2-D histo clus Dx vs Dz before
    TH2F			*hClusDxDZB;		//! 2-D histo clus Dx vs Dz after

    TH3F			*hdifEMC_EPV0;		//! 3-D histo dif phi in EMC with EPV0
    TH3F			*hdifEMC_EPV0A;		//! 3-D histo dif phi in EMC with EPV0A
    TH3F			*hdifEMC_EPV0C;		//! 3-D histo dif phi in EMC with EPV0C

    TH3F			*hdifful_EPV0;		//! 3-D histo dif phi in full with EPV0
    TH3F			*hdifful_EPV0A;		//! 3-D histo dif phi in full with EPV0A
    TH3F			*hdifful_EPV0C;		//! 3-D histo dif phi in full with EPV0C

    TH3F			*hdifout_EPV0;		//! 3-D histo dif phi out EMC with EPV0
    TH3F			*hdifout_EPV0A;		//! 3-D histo dif phi out EMC with EPV0A
    TH3F			*hdifout_EPV0C;		//! 3-D histo dif phi out EMC with EPV0C

    TH3F			*hdifEMC_EPTPC;		//! 3-D histo dif phi in EMC with EPTPC
    TH3F			*hdifful_EPTPC;		//! 3-D histo dif phi in full with EPTPC
    TH3F			*hdifout_EPTPC;		//! 3-D histo dif phi out EMC with EPTPC

    THnSparse		        *fClusterPbV0;
    THnSparse		        *fClusterPbV0A;
    THnSparse		        *fClusterPbV0C;
    THnSparse			 *fClusterPbTPC;

    THnSparse                   *fHEPV0r;	        //! Flow 4-D Histo
    THnSparse                   *fHEPV0A;	        //! Flow 4-D Histo
    THnSparse                   *fHEPV0C;	        //! Flow 4-D Histo
    THnSparse                   *fHEPTPC;	        //! Flow 4-D Histo

    
    
    AliAnalysisTaskPi0V2(const AliAnalysisTaskPi0V2&); // not implemented
    AliAnalysisTaskPi0V2& operator=(const AliAnalysisTaskPi0V2&); // not implemented
    
    ClassDef(AliAnalysisTaskPi0V2, 5); // example of analysis
};

#endif

