/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */
#ifndef AliMCTrackInfoContainer_H
#define AliMCTrackInfoContainer_H

#include "TClonesArray.h"

#include "AliPIDResponse.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"

class AliMCTrackInfoContainer : public TObject
{
public:
    AliMCTrackInfoContainer(); // default constructor
    AliMCTrackInfoContainer(const AliMCTrackInfoContainer &source); // copy constructor 
    virtual ~AliMCTrackInfoContainer();//destructor.
    
    void Fill(const AliAODTrack* primTrack, const AliPIDResponse* fPIDResponse, Double_t lBestPrimaryVtxPos[3], Double_t lMagField, const AliMCEvent* lMCEvent);
    void Reset();
    
    void Copy(AliMCTrackInfoContainer *source);
    
    Float_t GetLengthInActiveZone(const AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b );
    Float_t GetDCAz(AliAODTrack *lTrack);
    
    Bool_t  CheckITSClusterRequirement(AliESDtrackCuts::ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2);

    Bool_t IsSelected( const AliESDtrackCuts* lPrimTrackCuts );
    //setter 
    void SetPassesTrackCuts2010         (Bool_t value) { fTreePrimVarPassesTrackCuts2010            = value ;}
    void SetPassesTrackCuts2011         (Bool_t value) { fTreePrimVarPassesTrackCuts2011            = value ;}
    void SetPassesTrackCutsTPCRefit     (Bool_t value) { fTreePrimVarPassesTrackCutsTPCRefit        = value ;}
    void SetPassesTrackCuts2011Sys      (Bool_t value) { fTreePrimVarPassesTrackCuts2011Sys         = value ;}
    void SetPassesTrackCutsV0           (Bool_t value) { fTreePrimVarPassesTrackCutsV0              = value ;}
    void SetPassesTrackCutsHybrid_kOff  (Bool_t value) { fTreePrimVarPassesTrackCutsHybrid_kOff     = value ;}
    void SetPassesTrackCutsHybrid_kNone (Bool_t value) { fTreePrimVarPassesTrackCutsHybrid_kNone    = value ;}
    
    //getter
    //-----------BASIC-INFO---------------------------
    Int_t GetCharge                             () const { return  fTreePrimVarCharge    ;}
    Double_t GetRapPion                         () const { return  fTreePrimVarRapPion    ;}
    Double_t GetRapProton                       () const { return  fTreePrimVarRapProton    ;}
    Double_t GetRapKaon                         () const { return  fTreePrimVarRapKaon    ;}
    Double_t GetEta                             () const { return  fTreePrimVarEta    ;}
    Double_t GetPtot                            () const { return  fTreePrimVarPtot    ;}
    Double_t GetPt                              () const { return  fTreePrimVarPt    ;}
    
    Double_t GetDCAxyToPV                       () const { return  fTreePrimVarDCAxyToPV    ;}
    Double_t GetDCAzToPV                        () const { return  fTreePrimVarDCAzToPV    ;}
    
    Int_t GetNbrCrossedRows                     () const { return  fTreePrimVarNbrCrossedRows    ;}
    Double_t GetRatioCrossedRowsOverFindable    () const { return  fTreePrimVarRatioCrossedRowsOverFindable    ;}
    Int_t GetITSNbrClusters                     () const { return  fTreePrimVarITSNbrClusters    ;}
    Int_t GetTPCNbrClusters                     () const { return  fTreePrimVarTPCNbrClusters    ;}
    Double_t GetNbrCrossedRowsOverLength        () const { return  fTreePrimVarNbrCrossedRowsOverLength    ;}
    Double_t GetFractionSharedTPCClusters       () const { return  fTreePrimVarFractionSharedTPCClusters    ;}
    Double_t GetChi2perNDF                      () const { return  fTreePrimVarChi2perNDF; }
    Double_t GetITSChi2PerCluster               () const { return  fTreePrimVarITSChi2PerCluster    ;}
    Double_t GetTPCChi2PerCluster               () const { return  fTreePrimVarTPCChi2PerCluster    ;}
    Double_t GetChi2TPCConstainedVsGlobal       () const { return  fTreePrimVarChi2TPCConstrainedVsGlobal;}
    Double_t GetTrackLength                     () const { return  fTreePrimVarTrackLength    ;}
    //---------------PID-TPC-INFO---------------------
    Float_t GetNSigmaPion                       () const { return  fTreePrimVarNSigmaPion    ;}
    Float_t GetNSigmaKaon                       () const { return  fTreePrimVarNSigmaKaon    ;}
    Float_t GetNSigmaProton                     () const { return  fTreePrimVarNSigmaKaon    ;}
    //---------------PID-TOF-INFO---------------------
    Float_t GetTOFNSigmaPion                    () const { return  fTreePrimVarTOFNSigmaPion    ;}
    Float_t GetTOFNSigmaKaon                    () const { return  fTreePrimVarTOFNSigmaKaon    ;}
    Float_t GetTOFNSigmaProton                  () const { return  fTreePrimVarTOFNSigmaProton    ;}
    //---------------PID-ITS-INFO---------------------
    Float_t GetITSNSigmaPion                    () const { return  fTreePrimVarITSNSigmaPion    ;}
    Float_t GetITSNSigmaKaon                    () const { return  fTreePrimVarITSNSigmaKaon    ;}
    Float_t GetITSNSigmaProton                  () const { return  fTreePrimVarITSNSigmaProton    ;}
    //---Raw TPC dEdx + PIDForTracking information----
    Double_t GetdEdx                            () const { return  fTreePrimVardEdx    ;}
    Double_t GetdEdxN                           () const { return  fTreePrimVardEdxN    ;}
    Double_t GetPIDForTracking                  () const { return  fTreePrimVarPIDForTracking    ;}
    //------------------------------------------------
    ULong64_t GetTrackStatus                    () const { return  fTreePrimVarTrackStatus    ;}
    UInt_t GetFilterMap                    		() const { return  fTreePrimVarFilterMap    ;}
    //------------------------------------------------
    Int_t GetLabel                              () const { return  fTreePrimVarLabel    ;}
    //------------FULL-MOMENTUM-INFO------------------
    Double_t GetPx                              () const { return  fTreePrimVarPx    ;}
    Double_t GetPy                              () const { return  fTreePrimVarPy    ;}
    Double_t GetPz                              () const { return  fTreePrimVarPz    ;}
    //---------------CLUSTER-INFO---------------------
    Bool_t GetITSClusters0                      () const { return  fTreePrimVarITSClusters0    ;}
    Bool_t GetITSClusters1                      () const { return  fTreePrimVarITSClusters1    ;}
    Bool_t GetITSClusters2                      () const { return  fTreePrimVarITSClusters2    ;}
    Bool_t GetITSClusters3                      () const { return  fTreePrimVarITSClusters3    ;}
    Bool_t GetITSClusters4                      () const { return  fTreePrimVarITSClusters4    ;}
    Bool_t GetITSClusters5                      () const { return  fTreePrimVarITSClusters5    ;}
    //------------------------------------------------
    Bool_t GetITSSharedClusters0                () const { return  fTreePrimVarITSSharedClusters0    ;}
    Bool_t GetITSSharedClusters1                () const { return  fTreePrimVarITSSharedClusters1    ;}
    Bool_t GetITSSharedClusters2                () const { return  fTreePrimVarITSSharedClusters2    ;}
    Bool_t GetITSSharedClusters3                () const { return  fTreePrimVarITSSharedClusters3    ;}
    Bool_t GetITSSharedClusters4                () const { return  fTreePrimVarITSSharedClusters4    ;}
    Bool_t GetITSSharedClusters5                () const { return  fTreePrimVarITSSharedClusters5    ;}

    //---------------OOB-PILEUP-INFO---------------------
    Double_t GetTOFExpTDiff                     () const { return  fTreePrimVarTOFExpTDiff    ;}
    Double_t GetTOFSignal                       () const { return  fTreePrimVarTOFSignal    ;}
    Int_t   GetTOFBCid                          () const { return  fTreePrimVarTOFBCid    ;}
    
    Double_t GetTOFLength                       () const { return  fTreePrimVarTOFLength    ;}
    Double_t GetTOFDeltaX                       () const { return  fTreePrimVarTOFDeltaX    ;}
    Double_t GetTOFDeltaZ                       () const { return  fTreePrimVarTOFDeltaZ    ;}
    
    //Kink tagging
    Bool_t GetIsKink                            () const { return  fTreePrimVarIsKink    ;}
    
    //------------TRACK-CUTS-CHECK------------------
    Bool_t GetPassesTrackCuts2010               () const { return  fTreePrimVarPassesTrackCuts2010;}
    Bool_t GetPassesTrackCuts2011               () const { return  fTreePrimVarPassesTrackCuts2011;}
    Bool_t GetPassesTrackCutsTPCRefit           () const { return  fTreePrimVarPassesTrackCutsTPCRefit;}
    Bool_t GetPassesTrackCuts2011Sys            () const { return  fTreePrimVarPassesTrackCuts2011Sys;}
    Bool_t GetPassesTrackCutsV0                 () const { return  fTreePrimVarPassesTrackCutsV0;}
    Bool_t GetPassesTrackCutsHybrid_kOff        () const { return  fTreePrimVarPassesTrackCutsHybrid_kOff;}
    Bool_t GetPassesTrackCutsHybrid_kNone       () const { return  fTreePrimVarPassesTrackCutsHybrid_kNone ;}
    
    Double_t GetPositionX                       () const { return  fTreePrimVarPosX     ;}
    Double_t GetPositionY                       () const { return  fTreePrimVarPosY     ;}
    Double_t GetPositionZ                       () const { return  fTreePrimVarPosZ     ;}
    
    //------------MONTE-CARLO-INFO------------------
    Int_t GetPdgCode                            () const { return  fTreePrimVarPdgCode    ;}
    Double_t GetMassMC                          () const { return  fTreePrimVarMassMC    ;}
    Double_t GetRapMC                           () const { return  fTreePrimVarRapMC    ;}
    Double_t GetEtaMC                           () const { return  fTreePrimVarEtaMC    ;}
    Double_t GetPtotMC                          () const { return  fTreePrimVarPtotMC    ;}
    Double_t GetPtMC                            () const { return  fTreePrimVarPtMC    ;}
    Double_t GetPxMC                            () const { return  fTreePrimVarPxMC    ;}
    Double_t GetPyMC                            () const { return  fTreePrimVarPyMC    ;}
    Double_t GetPzMC                            () const { return  fTreePrimVarPzMC    ;}

    
    Int_t GetMotherPdgCode                      () const { return  fTreePrimVarMotherPdgCode                    ;}
    Int_t GetGrandMotherPdgCode                 () const { return  fTreePrimVarGrandMotherPdgCode               ;}
    Int_t GetMotherLabel                        () const { return  fTreePrimVarMotherLabel                      ;}
    Int_t GetGrandMotherLabel                   () const { return  fTreePrimVarGrandMotherLabel                 ;}

    Bool_t GetIsPhysicalPrimary                 () const { return  fTreePrimVarIsPhysicalPrimary                ;}
    Bool_t GetIsPhysicalPrimaryMother           () const { return  fTreePrimVarIsPhysicalPrimaryMother          ;}
    Bool_t GetIsPhysicalPrimaryGrandMother      () const { return  fTreePrimVarIsPhysicalPrimaryGrandMother     ;}
    
private:
    //===========================================================================================
    //   Variables for Primary tracks Tree
    //===========================================================================================
    //-----------BASIC-INFO---------------------------
    Int_t fTreePrimVarCharge; //
    Double_t fTreePrimVarRapPion; //
    Double_t fTreePrimVarRapProton; //
    Double_t fTreePrimVarRapKaon; //
    Double_t fTreePrimVarEta; //
    Double_t fTreePrimVarPtot; //
    Double_t fTreePrimVarPt; //
    
    Double_t fTreePrimVarDCAxyToPV;//
    Double_t fTreePrimVarDCAzToPV;//
    
    Int_t fTreePrimVarNbrCrossedRows;//
    Double_t fTreePrimVarRatioCrossedRowsOverFindable;//
    Int_t fTreePrimVarITSNbrClusters;//
    Int_t fTreePrimVarTPCNbrClusters;//
    Double_t fTreePrimVarNbrCrossedRowsOverLength;//
    Double_t fTreePrimVarFractionSharedTPCClusters;//
    Double_t fTreePrimVarChi2perNDF;//          Max track chi2
    Double_t fTreePrimVarITSChi2PerCluster;//
    Double_t fTreePrimVarTPCChi2PerCluster;//
    Double_t fTreePrimVarChi2TPCConstrainedVsGlobal;//
    Double_t fTreePrimVarTrackLength;  //  
    //---------------PID-TPC-INFO---------------------
    Float_t fTreePrimVarNSigmaPion;//
    Float_t fTreePrimVarNSigmaKaon;//
    Float_t fTreePrimVarNSigmaProton;//
    //---------------PID-TOF-INFO---------------------
    Float_t fTreePrimVarTOFNSigmaPion;//
    Float_t fTreePrimVarTOFNSigmaKaon;//
    Float_t fTreePrimVarTOFNSigmaProton;//
    //---------------PID-ITS-INFO---------------------
    Float_t fTreePrimVarITSNSigmaPion;//
    Float_t fTreePrimVarITSNSigmaKaon;//
    Float_t fTreePrimVarITSNSigmaProton;//
    //---Raw TPC dEdx + PIDForTracking information----
    Double_t fTreePrimVardEdx;//
    Double_t fTreePrimVardEdxN;//
    Double_t fTreePrimVarPIDForTracking;//
    //------------------------------------------------
    ULong64_t fTreePrimVarTrackStatus;//
    UInt_t fTreePrimVarFilterMap;//
    //------------------------------------------------
    Int_t fTreePrimVarLabel;//
    //------------FULL-MOMENTUM-INFO------------------
    Double_t fTreePrimVarPx; //
    Double_t fTreePrimVarPy; //
    Double_t fTreePrimVarPz; //
    //---------------CLUSTER-INFO---------------------
    Bool_t fTreePrimVarITSClusters0;//
    Bool_t fTreePrimVarITSClusters1;//
    Bool_t fTreePrimVarITSClusters2;//
    Bool_t fTreePrimVarITSClusters3;//
    Bool_t fTreePrimVarITSClusters4;//
    Bool_t fTreePrimVarITSClusters5;//
    //------------------------------------------------
    Bool_t fTreePrimVarITSSharedClusters0;//
    Bool_t fTreePrimVarITSSharedClusters1;//
    Bool_t fTreePrimVarITSSharedClusters2;//
    Bool_t fTreePrimVarITSSharedClusters3;//
    Bool_t fTreePrimVarITSSharedClusters4;//
    Bool_t fTreePrimVarITSSharedClusters5;//

    //---------------OOB-PILEUP-INFO---------------------
    Double_t fTreePrimVarTOFExpTDiff; //
    Double_t fTreePrimVarTOFSignal; //
    Int_t   fTreePrimVarTOFBCid; //
    
    Double_t fTreePrimVarTOFLength;//
    Double_t fTreePrimVarTOFDeltaX;//
    Double_t fTreePrimVarTOFDeltaZ;//
    
    //Kink tagging
    Bool_t fTreePrimVarIsKink;//
    
    //------------TRACK-CUTS-CHECK------------------
    Bool_t fTreePrimVarPassesTrackCuts2010;//
    Bool_t fTreePrimVarPassesTrackCuts2011;//
    Bool_t fTreePrimVarPassesTrackCutsTPCRefit;//
    Bool_t fTreePrimVarPassesTrackCuts2011Sys;//
    Bool_t fTreePrimVarPassesTrackCutsV0;//
    Bool_t fTreePrimVarPassesTrackCutsHybrid_kOff;//
    Bool_t fTreePrimVarPassesTrackCutsHybrid_kNone;//
    
    Double_t fTreePrimVarPosX;//
    Double_t fTreePrimVarPosY;//
    Double_t fTreePrimVarPosZ;//
    
    //------------MONTE-CARLO-INFO------------------
    Int_t fTreePrimVarPdgCode;//
    Double_t fTreePrimVarMassMC;//
    Double_t fTreePrimVarRapMC;//
    Double_t fTreePrimVarEtaMC;//
    Double_t fTreePrimVarPtotMC;//
    Double_t fTreePrimVarPtMC;//
    Double_t fTreePrimVarPxMC;//
    Double_t fTreePrimVarPyMC;//
    Double_t fTreePrimVarPzMC;//
    
    Int_t fTreePrimVarMotherPdgCode;//
    Int_t fTreePrimVarGrandMotherPdgCode;//
    Int_t fTreePrimVarMotherLabel;//
    Int_t fTreePrimVarGrandMotherLabel;//

    Bool_t fTreePrimVarIsPhysicalPrimary;//
    Bool_t fTreePrimVarIsPhysicalPrimaryMother;//
    Bool_t fTreePrimVarIsPhysicalPrimaryGrandMother;//
    
    //------------------------------------------------
    ClassDef(AliMCTrackInfoContainer, 1);
};
#endif
