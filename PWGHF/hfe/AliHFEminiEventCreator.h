/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
//------------------------------------------------------------------------------------------------
//  This is a derived version of  AliHFEreducedEventCreatorESD.cxx/h by M.Fasel <M.Fasel@gsi.de>
//  The contents are much less than  AliHFEreducedEventCreatorESD class.
//  Contributors: Nirbhay K. Behera, Jiyeon Kwon, Jonghan Park
//  At the moment, keep information only for HFE analysis by IP fit method.
//  Works both for ESD and AOD. Collision system: pp, pA and AA.
//------------------------------------------------------------------------------------------------


#ifndef ALIHFEREDUCEDEVENTCREATOR_H
#define ALIHFEREDUCEDEVENTCREATOR_H

#include "AliAnalysisTaskSE.h"

class TString;
class TTree;
class TH1D;
class AliMCEvent;
class AliVEvent;
class AliVParticle;
class AliAODMCHeader;
class TClonesArray;
class AliAnalysisUtils;
class AliPIDResponse;
class AliHFEcuts;
class AliHFEextraCuts;
class AliHFEpidTPC;
class AliHFEsignalCuts;
class AliHFEV0taginfo;
class AliHFEminiEvent;
class AliHFEminiTrack;

class AliHFEminiEventCreator : public AliAnalysisTaskSE{
  public:
    AliHFEminiEventCreator();
    AliHFEminiEventCreator(const char *name);
    virtual ~AliHFEminiEventCreator();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);

    // Setters for cuts
    void SetIsMCEvent(Bool_t IsMCevent){ fHasMCdata = IsMCevent; }
    void SetChi2TPCCut( Double_t chi2TPC ){ fTPCChiSquare = chi2TPC; }
    void SetMinClusterTPC( Int_t TPCCluster ) { fNclustersTPC = TPCCluster; }
    void SetMinClusterTPCPID( Int_t MinTPCclusterPID ){ fNclustersTPCPID = MinTPCclusterPID; }
    void SetClusterRatioTPC( Double_t RatioTPCNcluster ) { fRatioTPCNcluster = RatioTPCNcluster; }
    void SetClusterITS( Int_t MinNclusterITS ) { fNclustersITS = MinNclusterITS; }
    void SetITSLayerStatus( Bool_t layerStatus ) { fITSLayerStatus = layerStatus; }
    void SetKinematicsCuts( Double_t ptMin, Double_t ptMax, Double_t eta ){ fPtMin = ptMin; fPtMax = ptMax; fEta = eta; }
    void SetVertexZCut( Double_t Vz ) { fVz = Vz; }
    void SetDCAcut( Double_t dcaxy, Double_t dcaz ) { fDCAxy = dcaxy; fDCAz = dcaz; }
    void SetProductionVzToSPDVz( Double_t prodVz ) { fProdVz = prodVz; }
    void SetSPDResolution( Double_t spdResolution ) { fSPDResolution = spdResolution; }
    void SetTPCnSigma( Double_t NsigmaTPCLow, Double_t NsigmaTPCHigh ){
      fNsigmaTPClow = NsigmaTPCLow; fNsigmaTPChigh = NsigmaTPCHigh; }
    void SetTOFnSigma( Double_t nSigmaTOF){ fNsigmaTOF = nSigmaTOF; }
    
    void SetRemoveFirstEventFromChunk() { fRemoveFirstEvent = kTRUE; }
    void SetCollisionSystem( TString CollisionSystem ) { fCollisionSystem = CollisionSystem; }  
    void SetHFECuts(AliHFEcuts * const cuts) {fTrackCuts = cuts;}
    
    AliHFEpidTPC *GetTPCResponse() { return fTPCpid; }
    
    Bool_t IsTOFmismatch(const AliVTrack *const track, const AliPIDResponse *const pid) const;
    
 private:
    AliHFEminiEventCreator(const AliHFEminiEventCreator &);
    AliHFEminiEventCreator &operator=(const AliHFEminiEventCreator &);

    void                 SetUpMCEvent();                             //setter for MC events
    void                 SetUpESDAnalysis();                         //setter for ESD analysis
    void                 SetUpAODAnalysis();                         //Setter for AOD analysis
    void                 ExecAODEvent();                             //AOD analysis
    void                 ExecESDEvent();                             //ESD analysis
    Bool_t               RemovePileupEvents(AliVEvent *event);       //Remove pileup events--
    Bool_t               ProperVertex(AliVEvent *event);             //Vz and SPDVz cut
    void                 GetSPDnVZEROMultiplicity(AliVEvent *event); //SPD and VZERO multiplicity
    void                 GetEventCentrality(AliVEvent *event);       //Get Event centrality
    Int_t                GetMCPrimaryNch();                          //Primary Nch from MC 
    
   
    Int_t                fEventNumber;               // Event Number 
    Int_t                fNclustersTPC;              // Min Number of clusters in TPC
    Int_t                fNclustersTPCPID;           // Min Number of clusters for TPC PID
    Int_t                fNclustersITS;              // Min Number of clusters in ITS
    Double_t             fTPCChiSquare;              //TPC chi2 cut
    Double_t             fRatioTPCNcluster;          // TPC ncluster ratio
    Double_t             fPtMin;                     //lower pT cutoff for track
    Double_t             fPtMax;                     // higher pT cutoff for track
    Double_t             fEta;                       // eta cut for track
    Double_t             fVz;                        //event z-vertex cut
    Double_t             fDCAxy;                     //DCAxy cut 
    Double_t             fDCAz;                      //DCAz cut 
    Double_t             fProdVz;                    //Distance of Vz with respect to SPD Vz
    Double_t             fSPDResolution;             //SPDVz resolution 
    Double_t             fNsigmaTPClow;              //nsigmaTPC lower cut
    Double_t             fNsigmaTPChigh;             //nsigmaTPC uppper cut
    Double_t             fNsigmaTOF;                 //symmetric cut for TOF nsigma
    Bool_t               fITSLayerStatus;            // ITS layer status
    Bool_t               fHasMCdata;                 // true when MC events found
    Bool_t               fAODanalysis;               // true for AOD event analysis
    Bool_t               fRemoveFirstEvent;          // Remove first event from chunk
    Bool_t               fSelectSignalOnly;          // Select signal-only tracks
    TString              fCollisionSystem;           //pp or AA (pPb, PbPb)
   
    TList                *fList;                     //List
    TTree                *fHFEtree;                  // HFE tree 
    TH1D                 *fNevents;
    TH1D                 *fNtracks;
    
    //--------------Common classes-------------------
    AliESDEvent          *fESD;                    // for ESD event
    AliAODEvent          *fAOD;                    // for AOD event
    AliVEvent 	         *fVevent;                 // V event
    AliMCEvent           *fMCEvent;                // for MC EVent
    AliStack             *fStack;                  // for MC Stack
    AliPIDResponse       *fPIDResponse;            // for PID response object
    AliAODMCHeader       *fAODMCHeader;            // for MC info AOD
    TClonesArray         *fAODArrayMCInfo;         // for MC info particle AOD
    AliVEventHandler     *fInputHandler;            // for Event handler
    AliMCEventHandler    *fMCEventHandler;         // for MC event handler---
    AliAnalysisUtils     *fAnalysisUtils;          // To remove the first event, pile up

    //---------------HFE classes------------------------
    AliHFEminiEvent      *fHFEevent;               // hfe event
    AliHFEcuts           *fTrackCuts;              // Track
    AliHFEextraCuts      *fExtraCuts;              // HFE IP info
    AliHFEsignalCuts     *fSignalCuts;             // Signal Cuts
    AliHFEpidTPC         *fTPCpid;                 // TPC PID
    AliHFEV0taginfo      *fV0Tagger;               // Tags v0 tracks per Event
    
    
    ClassDef(AliHFEminiEventCreator, 2)
      };
#endif


