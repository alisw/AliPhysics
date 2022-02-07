/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */
#ifndef AliV0Container_H
#define AliV0Container_H

#include "TClonesArray.h"

#include "AliAODv0.h"
#include "AliPIDResponse.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliV0Result.h"

class AliV0Container : public TObject
{
public:
    AliV0Container(); // default constructor
    virtual ~AliV0Container();//destructor.
    
    void Fill(AliAODv0* v0, AliPIDResponse* fPIDResponse, Double_t lBestPrimaryVtxPos[3], Double_t lMagField);
    void Reset();
    
    Float_t GetLengthInActiveZone( AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b );
    Float_t GetDCAz(AliAODTrack *lTrack);

    Bool_t IsSelected( const AliV0Result* lV0Result );
    //setter 
    
    
    //getter
    //-----------BASIC-INFO---------------------------
    Bool_t   GetOnFlyStatus                         () const { return  fTreeV0VarOnFlyStatus;}
    Double_t GetMassAsK0s                           () const { return  fTreeV0VarMassAsK0s   ;}
    Double_t GetMassAsLambda                        () const { return  fTreeV0VarMassAsLambda   ;}
    Double_t GetMassAsAntiLambda                    () const { return  fTreeV0VarMassAsAntiLambda   ;}
    Double_t GetRapK0Short                          () const { return  fTreeV0VarRapK0Short   ;}
    Double_t GetRapLambda                           () const { return  fTreeV0VarRapLambda   ;}
    Double_t GetV0Eta                               () const { return  fTreeV0VarV0Eta   ;}
    Double_t GetPosEta                              () const { return  fTreeV0VarPosEta   ;}
    Double_t GetNegEta                              () const { return  fTreeV0VarNegEta   ;}
    Double_t GetPtot                                () const { return  fTreeV0VarPtot   ;}
    Double_t GetPt                                  () const { return  fTreeV0VarPt   ;}
    
    //-----------INFO-FOR-CUTS------------------------
    Double_t GetAlpha                               () const { return  fTreeV0VarAlpha   ;}
    Double_t GetPtArm                               () const { return  fTreeV0VarPtArm   ;}
    Double_t GetDCAV0Dau                            () const { return  fTreeV0VarDCAV0Dau   ;}
    Double_t GetDCAV0ToPV                           () const { return  fTreeV0VarDCAV0ToPV   ;}
    Double_t GetDCAPosToPV                          () const { return  fTreeV0VarDCAPosToPV   ;}
    Double_t GetDCANegToPV                          () const { return  fTreeV0VarDCANegToPV   ;}
    Double_t GetCosPA                               () const { return  fTreeV0VarCosPA   ;}
    Double_t GetRadius                              () const { return  fTreeV0VarRadius   ;}
    
    Int_t   GetLeastNbrCrossedRows                  () const { return  fTreeV0VarLeastNbrCrossedRows   ;}
    Int_t   GetLeastNbrClusters                     () const { return  fTreeV0VarLeastNbrClusters   ;}
    Double_t GetLeastRatioCrossedRowsOverFindable   () const { return  fTreeV0VarLeastRatioCrossedRowsOverFindable   ;}
    Double_t GetMaxChi2PerCluster                   () const { return  fTreeV0VarMaxChi2PerCluster   ;}
    Double_t GetMinTrackLength                      () const { return  fTreeV0VarMinTrackLength   ;}
    //-----------DECAY-LENGTH-INFO--------------------
    Double_t GetDistOverTotMom                      () const { return  fTreeV0VarDistOverTotMom   ;}
    //---------------PID-TPC-INFO---------------------
    Float_t GetPosNSigmaProton                      () const { return  fTreeV0VarPosNSigmaProton   ;}
    Float_t GetPosNSigmaPion                        () const { return  fTreeV0VarPosNSigmaPion   ;}
    Float_t GetNegNSigmaProton                      () const { return  fTreeV0VarNegNSigmaProton   ;}
    Float_t GetNegNSigmaPion                        () const { return  fTreeV0VarNegNSigmaPion   ;}
    //---------------PID-TOF-INFO---------------------
    Float_t GetPosTOFNSigmaProton                   () const { return  fTreeV0VarPosTOFNSigmaProton   ;}
    Float_t GetPosTOFNSigmaPion                     () const { return  fTreeV0VarPosTOFNSigmaPion   ;}
    Float_t GetNegTOFNSigmaProton                   () const { return  fTreeV0VarNegTOFNSigmaProton   ;}
    Float_t GetNegTOFNSigmaPion                     () const { return  fTreeV0VarNegTOFNSigmaPion   ;}
    //---------------PID-ITS-INFO---------------------
    Float_t GetPosITSNSigmaProton                   () const { return  fTreeV0VarPosITSNSigmaProton   ;}
    Float_t GetPosITSNSigmaPion                     () const { return  fTreeV0VarPosITSNSigmaPion   ;}
    Float_t GetNegITSNSigmaProton                   () const { return  fTreeV0VarNegITSNSigmaProton   ;}
    Float_t GetNegITSNSigmaPion                     () const { return  fTreeV0VarNegITSNSigmaPion   ;}
    //---Raw TPC dEdx + PIDForTracking information----
    Double_t GetPosdEdx                             () const { return  fTreeV0VarPosdEdx   ;}
    Double_t GetNegdEdx                             () const { return  fTreeV0VarNegdEdx   ;}
    Double_t GetPosdEdxN                            () const { return  fTreeV0VarPosdEdxN   ;}
    Double_t GetNegdEdxN                            () const { return  fTreeV0VarNegdEdxN   ;}
    Double_t GetPosPIDForTracking                   () const { return  fTreeV0VarPosPIDForTracking   ;} 
    Double_t GetNegPIDForTracking                   () const { return  fTreeV0VarNegPIDForTracking   ;}
    //------------------------------------------------
    Double_t GetChi2V0                              () const {return fTreeV0VarChi2V0    ;}
    //------------------------------------------------
    ULong64_t GetPosTrackStatus                     () const { return  fTreeV0VarPosTrackStatus   ;}
    ULong64_t GetNegTrackStatus                     () const { return  fTreeV0VarNegTrackStatus   ;}
    //------------------------------------------------
    Int_t GetPosLabel                               () const { return  fTreeV0VarPosLabel   ;}
    Int_t GetNegLabel                               () const { return  fTreeV0VarNegLabel   ;}
    //------------------------------------------------
    Double_t GetPosDCAz                             () const { return  fTreeV0VarPosDCAz   ;}
    Double_t GetNegDCAz                             () const { return  fTreeV0VarNegDCAz   ;}
    //------------FULL-MOMENTUM-INFO------------------
    Double_t GetPosPx                               () const { return  fTreeV0VarPosPx   ;}
    Double_t GetPosPy                               () const { return  fTreeV0VarPosPy   ;}
    Double_t GetPosPz                               () const { return  fTreeV0VarPosPz   ;}
    Double_t GetNegPx                               () const { return  fTreeV0VarNegPx   ;}
    Double_t GetNegPy                               () const { return  fTreeV0VarNegPy   ;}
    Double_t GetNegPz                               () const { return  fTreeV0VarNegPz   ;}
    //------------------------------------------------
    Double_t GetV0DecayX                            () const { return  fTreeV0VarV0DecayX   ;}
    Double_t GetV0DecayY                            () const { return  fTreeV0VarV0DecayY   ;}
    Double_t GetV0DecayZ                            () const { return  fTreeV0VarV0DecayZ   ;}
    //---------------CLUSTER-INFO---------------------
    Bool_t GetPosITSClusters0                       () const { return  fTreeV0VarPosITSClusters0   ;}
    Bool_t GetPosITSClusters1                       () const { return  fTreeV0VarPosITSClusters1   ;}
    Bool_t GetPosITSClusters2                       () const { return  fTreeV0VarPosITSClusters2   ;}
    Bool_t GetPosITSClusters3                       () const { return  fTreeV0VarPosITSClusters3   ;}
    Bool_t GetPosITSClusters4                       () const { return  fTreeV0VarPosITSClusters4   ;}
    Bool_t GetPosITSClusters5                       () const { return  fTreeV0VarPosITSClusters5   ;}
    
    Bool_t GetNegITSClusters0                       () const { return  fTreeV0VarNegITSClusters0   ;}
    Bool_t GetNegITSClusters1                       () const { return  fTreeV0VarNegITSClusters1   ;}
    Bool_t GetNegITSClusters2                       () const { return  fTreeV0VarNegITSClusters2   ;}
    Bool_t GetNegITSClusters3                       () const { return  fTreeV0VarNegITSClusters3   ;}
    Bool_t GetNegITSClusters4                       () const { return  fTreeV0VarNegITSClusters4   ;}
    Bool_t GetNegITSClusters5                       () const { return  fTreeV0VarNegITSClusters5   ;}
    //------------------------------------------------
    Bool_t GetPosITSSharedClusters0                 () const { return  fTreeV0VarPosITSSharedClusters0   ;}
    Bool_t GetPosITSSharedClusters1                 () const { return  fTreeV0VarPosITSSharedClusters1   ;}
    Bool_t GetPosITSSharedClusters2                 () const { return  fTreeV0VarPosITSSharedClusters2   ;}
    Bool_t GetPosITSSharedClusters3                 () const { return  fTreeV0VarPosITSSharedClusters3   ;}
    Bool_t GetPosITSSharedClusters4                 () const { return  fTreeV0VarPosITSSharedClusters4   ;}
    Bool_t GetPosITSSharedClusters5                 () const { return  fTreeV0VarPosITSSharedClusters5   ;}
    
    Bool_t GetNegITSSharedClusters0                 () const { return  fTreeV0VarNegITSSharedClusters0   ;}
    Bool_t GetNegITSSharedClusters1                 () const { return  fTreeV0VarNegITSSharedClusters1   ;}
    Bool_t GetNegITSSharedClusters2                 () const { return  fTreeV0VarNegITSSharedClusters2   ;}
    Bool_t GetNegITSSharedClusters3                 () const { return  fTreeV0VarNegITSSharedClusters3   ;}
    Bool_t GetNegITSSharedClusters4                 () const { return  fTreeV0VarNegITSSharedClusters4   ;}
    Bool_t GetNegITSSharedClusters5                 () const { return  fTreeV0VarNegITSSharedClusters5   ;}
    //---------------OOB-PILEUP-INFO---------------------
    Double_t GetNegTOFExpTDiff                      () const { return  fTreeV0VarNegTOFExpTDiff   ;}
    Double_t GetPosTOFExpTDiff                      () const { return  fTreeV0VarPosTOFExpTDiff   ;}
    Double_t GetNegTOFSignal                        () const { return  fTreeV0VarNegTOFSignal   ;}
    Double_t GetPosTOFSignal                        () const { return  fTreeV0VarPosTOFSignal   ;}
    Int_t   GetNegTOFBCid                           () const { return  fTreeV0VarNegTOFBCid   ;}
    Int_t   GetPosTOFBCid                           () const { return  fTreeV0VarPosTOFBCid   ;} 
    
    Double_t GetPosTOFLength                        () const { return  fTreeV0VarPosTOFLength   ;}
    Double_t GetNegTOFLength                        () const { return  fTreeV0VarNegTOFLength   ;}
    
    Double_t GetPosTOFDeltaX                        () const { return  fTreeV0VarPosTOFDeltaX   ;}
    Double_t GetNegTOFDeltaX                        () const { return  fTreeV0VarNegTOFDeltaX   ;}
    
    Double_t GetPosTOFDeltaZ                        () const { return  fTreeV0VarPosTOFDeltaZ   ;}
    Double_t GetNegTOFDeltaZ                        () const { return  fTreeV0VarNegTOFDeltaZ   ;}
    
    //Kink tag
    Bool_t GetPosIsKink                             () const { return  fTreeV0VarPosIsKink   ;}
    Bool_t GetNegIsKink                             () const { return  fTreeV0VarNegIsKink   ;}
    
    //Cowboy/sailor studies
    Bool_t GetIsCowboy                              () const { return  fTreeV0VarIsCowboy   ;} 
    
private:
    //===========================================================================================
    //   Variables for V0 Tree
    //===========================================================================================
    //-----------BASIC-INFO---------------------------
    Bool_t   fTreeV0VarOnFlyStatus;// if kTRUE, then this V0 is recontructed "on fly" during the tracking
    Double_t fTreeV0VarMassAsK0s; //
    Double_t fTreeV0VarMassAsLambda; //
    Double_t fTreeV0VarMassAsAntiLambda; //
    Double_t fTreeV0VarRapK0Short; //
    Double_t fTreeV0VarRapLambda; //
    Double_t fTreeV0VarV0Eta; //
    Double_t fTreeV0VarPosEta; //
    Double_t fTreeV0VarNegEta; //
    Double_t fTreeV0VarPtot; //
    Double_t fTreeV0VarPt; //
    
    //-----------INFO-FOR-CUTS------------------------
    Double_t fTreeV0VarAlpha; //
    Double_t fTreeV0VarPtArm;//
    Double_t fTreeV0VarDCAV0Dau; //
    Double_t fTreeV0VarDCAV0ToPV; //
    Double_t fTreeV0VarDCAPosToPV; //
    Double_t fTreeV0VarDCANegToPV; //
    Double_t fTreeV0VarCosPA; //
    Double_t fTreeV0VarRadius; //
    
    Int_t   fTreeV0VarLeastNbrCrossedRows;//
    Int_t   fTreeV0VarLeastNbrClusters;//
    Double_t fTreeV0VarLeastRatioCrossedRowsOverFindable;//
    Double_t fTreeV0VarMaxChi2PerCluster; //
    Double_t fTreeV0VarMinTrackLength; //
    //-----------DECAY-LENGTH-INFO--------------------
    Double_t fTreeV0VarDistOverTotMom;//
    //---------------PID-TPC-INFO---------------------
    Float_t fTreeV0VarPosNSigmaProton;//
    Float_t fTreeV0VarPosNSigmaPion;//
    Float_t fTreeV0VarNegNSigmaProton;//
    Float_t fTreeV0VarNegNSigmaPion;//
    //---------------PID-TOF-INFO---------------------
    Float_t fTreeV0VarPosTOFNSigmaProton;//
    Float_t fTreeV0VarPosTOFNSigmaPion;//
    Float_t fTreeV0VarNegTOFNSigmaProton;//
    Float_t fTreeV0VarNegTOFNSigmaPion;//
    //---------------PID-ITS-INFO---------------------
    Float_t fTreeV0VarPosITSNSigmaProton;//
    Float_t fTreeV0VarPosITSNSigmaPion;//
    Float_t fTreeV0VarNegITSNSigmaProton;//
    Float_t fTreeV0VarNegITSNSigmaPion;//
    //---Raw TPC dEdx + PIDForTracking information----
    Double_t fTreeV0VarPosdEdx;//
    Double_t fTreeV0VarNegdEdx;//
    Double_t fTreeV0VarPosdEdxN;//
    Double_t fTreeV0VarNegdEdxN;//
    Double_t fTreeV0VarPosPIDForTracking; // uses AliPID EParticleType code (0=electron, 1=muon, 2=pion, etc)
    Double_t fTreeV0VarNegPIDForTracking; //
    //------------------------------------------------
    Double_t fTreeV0VarChi2V0;         //
    //------------------------------------------------
    ULong64_t fTreeV0VarPosTrackStatus; //
    ULong64_t fTreeV0VarNegTrackStatus; //
    //------------------------------------------------
    Int_t fTreeV0VarPosLabel;//
    Int_t fTreeV0VarNegLabel;//
    //------------------------------------------------
    Double_t fTreeV0VarPosDCAz; //
    Double_t fTreeV0VarNegDCAz; //
    //------------FULL-MOMENTUM-INFO------------------
    Double_t fTreeV0VarPosPx; //
    Double_t fTreeV0VarPosPy; //
    Double_t fTreeV0VarPosPz; //
    Double_t fTreeV0VarNegPx; //
    Double_t fTreeV0VarNegPy; //
    Double_t fTreeV0VarNegPz; //
    //------------------------------------------------
    Double_t fTreeV0VarV0DecayX; //
    Double_t fTreeV0VarV0DecayY; //
    Double_t fTreeV0VarV0DecayZ; //
    //---------------CLUSTER-INFO---------------------
    Bool_t fTreeV0VarPosITSClusters0;//
    Bool_t fTreeV0VarPosITSClusters1;//
    Bool_t fTreeV0VarPosITSClusters2;//
    Bool_t fTreeV0VarPosITSClusters3;//
    Bool_t fTreeV0VarPosITSClusters4;//
    Bool_t fTreeV0VarPosITSClusters5;//
    
    Bool_t fTreeV0VarNegITSClusters0;//
    Bool_t fTreeV0VarNegITSClusters1;//
    Bool_t fTreeV0VarNegITSClusters2;//
    Bool_t fTreeV0VarNegITSClusters3;//
    Bool_t fTreeV0VarNegITSClusters4;//
    Bool_t fTreeV0VarNegITSClusters5;//
    //------------------------------------------------
    Bool_t fTreeV0VarPosITSSharedClusters0;//
    Bool_t fTreeV0VarPosITSSharedClusters1;//
    Bool_t fTreeV0VarPosITSSharedClusters2;//
    Bool_t fTreeV0VarPosITSSharedClusters3;//
    Bool_t fTreeV0VarPosITSSharedClusters4;//
    Bool_t fTreeV0VarPosITSSharedClusters5;//
    
    Bool_t fTreeV0VarNegITSSharedClusters0;//
    Bool_t fTreeV0VarNegITSSharedClusters1;//
    Bool_t fTreeV0VarNegITSSharedClusters2;//
    Bool_t fTreeV0VarNegITSSharedClusters3;//
    Bool_t fTreeV0VarNegITSSharedClusters4;//
    Bool_t fTreeV0VarNegITSSharedClusters5;//
    //---------------OOB-PILEUP-INFO---------------------
    Double_t fTreeV0VarNegTOFExpTDiff; //
    Double_t fTreeV0VarPosTOFExpTDiff; //
    Double_t fTreeV0VarNegTOFSignal; //
    Double_t fTreeV0VarPosTOFSignal; //
    Int_t   fTreeV0VarNegTOFBCid; //
    Int_t   fTreeV0VarPosTOFBCid; // 
    
    Double_t fTreeV0VarPosTOFLength;//
    Double_t fTreeV0VarNegTOFLength;//
    
    Double_t fTreeV0VarPosTOFDeltaX;//
    Double_t fTreeV0VarNegTOFDeltaX;//
    
    Double_t fTreeV0VarPosTOFDeltaZ;//
    Double_t fTreeV0VarNegTOFDeltaZ;//
    
    //Kink tag
    Bool_t fTreeV0VarPosIsKink;//
    Bool_t fTreeV0VarNegIsKink;//
    
    //Cowboy/sailor studies
    Bool_t fTreeV0VarIsCowboy; // store if V0 is cowboy-like or sailor-like in XY plane
    
    //------------------------------------------------
    ClassDef(AliV0Container, 1);
};
#endif
