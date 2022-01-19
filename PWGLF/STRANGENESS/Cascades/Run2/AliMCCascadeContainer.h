/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */
#ifndef AliMCCascadeContainer_H
#define AliMCCascadeContainer_H

#include "TClonesArray.h"

#include "AliAODcascade.h"
#include "AliPIDResponse.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliCascadeResult.h"

class AliMCCascadeContainer : public TObject
{
public:
    AliMCCascadeContainer(); // default constructor
    AliMCCascadeContainer(const AliMCCascadeContainer& source);
    virtual ~AliMCCascadeContainer();//destructor.
    
    void Fill(AliAODcascade* cascade, AliPIDResponse* fPIDResponse, Double_t lBestPrimaryVtxPos[3], Double_t lMagField, AliMCEvent* lMCEvent);
    void Reset();
    
    Float_t GetLengthInActiveZone( AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b );
    Float_t GetDCAz(AliAODTrack *lTrack);

    Bool_t IsSelected( const AliCascadeResult* lCascadeResult );
    //setter 
    
    
    //getter
    //-----------BASIC-INFO---------------------------
    Int_t GetCharge                                     () const{ return  fTreeCascVarCharge   ;}
    Double_t GetMassAsXi                                () const{ return  fTreeCascVarMassAsXi   ;}
    Double_t GetMassAsOmega                             () const{ return  fTreeCascVarMassAsOmega   ;}
    Double_t GetPtot                                    () const{ return  fTreeCascVarPtot   ;}
    Double_t GetV0Ptot                                  () const{ return  fTreeCascVarV0Ptot   ;}
    Double_t GetPt                                      () const{ return  fTreeCascVarPt   ;}
    Double_t GetV0Pt                                    () const{ return  fTreeCascVarV0Pt   ;}
    Double_t GetRapXi                                   () const{ return  fTreeCascVarRapXi   ;}
    Double_t GetRapOmega                                () const{ return  fTreeCascVarRapOmega   ;}
    Double_t GetCascEta                                 () const{ return  fTreeCascVarCascEta   ;}
    Double_t GetV0Eta                                   () const{ return  fTreeCascVarV0Eta   ;}
    Double_t GetBachEta                                 () const{ return  fTreeCascVarBachEta   ;}
    Double_t GetPosEta                                  () const{ return  fTreeCascVarPosEta   ;}
    Double_t GetNegEta                                  () const{ return  fTreeCascVarNegEta   ;}
    //-----------INFO-FOR-CUTS------------------------
    Double_t GetAlpha                                   () const{ return  fTreeCascVarAlpha   ;}
    Double_t GetPtArm                                   () const{ return  fTreeCascVarPtArm   ;}
    Double_t GetAlphaV0                                 () const{ return  fTreeCascVarAlphaV0   ;}
    Double_t GetPtArmV0                                 () const{ return  fTreeCascVarPtArmV0   ;}
    Double_t GetDCACascDau                              () const{ return  fTreeCascVarDCACascDau   ;}
    Double_t GetDCABachToPV                             () const{ return  fTreeCascVarDCABachToPV   ;}
    Double_t GetDCAV0Dau                                () const{ return  fTreeCascVarDCAV0Dau   ;}
    Double_t GetDCAV0ToPV                               () const{ return  fTreeCascVarDCAV0ToPV   ;}
    Float_t GetDCAPosToPV                               () const{ return  fTreeCascVarDCAPosToPV   ;}
    Float_t GetDCANegToPV                               () const{ return  fTreeCascVarDCANegToPV   ;}
    Double_t GetCascCosPA                               () const{ return  fTreeCascVarCascCosPA   ;}
    Double_t GetCascCosPASpecial                        () const{ return  fTreeCascVarCascCosPASpecial   ;}
    
    Double_t GetCascRadius                              () const{ return  fTreeCascVarCascRadius   ;}
    Double_t GetV0Mass                                  () const{ return  fTreeCascVarV0Mass   ;}
    Double_t GetV0MassAsLambda                          () const{ return  fTreeCascVarV0MassAsLambda   ;}
    Double_t GetV0MassAsAntiLambda                      () const{ return  fTreeCascVarV0MassAsAntiLambda   ;}
    Double_t GetV0Radius                                () const{ return  fTreeCascVarV0Radius   ;}
    Double_t GetV0CosPA                                 () const{ return  fTreeCascVarV0CosPA   ;}
    Double_t GetWrongCosPA                              () const{ return  fTreeCascVarWrongCosPA   ;}
    Double_t GetDCABachToBaryon                         () const{ return  fTreeCascVarDCABachToBaryon   ;}
    Int_t GetLeastNbrCrossedRows                        () const{ return  fTreeCascVarLeastNbrCrossedRows   ;}
    Int_t GetLeastNbrClusters                           () const{ return  fTreeCascVarLeastNbrClusters   ;}
    Double_t GetLeastRatioCrossedRowsOverFindable       () const{ return  fTreeCascVarLeastRatioCrossedRowsOverFindable   ;}
    Double_t GetNbrCrossedRowsOverLength                () const{ return  fTreeCascVarNbrCrossedRowsOverLength   ;}
    Double_t GetMaxChi2PerCluster                       () const{ return  fTreeCascVarMaxChi2PerCluster   ;}
    Double_t GetMinTrackLength                          () const{ return  fTreeCascVarMinTrackLength   ;}  
    //-----------DECAY-LENGTH-INFO--------------------
    Double_t GetDistOverTotMom                          () const{ return  fTreeCascVarDistOverTotMom   ;}
    Double_t GetV0DistOverTotMom                        () const{ return  fTreeCascVarV0DistOverTotMom   ;}
    //------------------------------------------------
    //---------------PID-TPC-INFO---------------------
    Float_t GetBachNSigmaPion                           () const{ return  fTreeCascVarBachNSigmaPion   ;}
    Float_t GetBachNSigmaKaon                           () const{ return  fTreeCascVarBachNSigmaKaon   ;}
    Float_t GetPosNSigmaProton                          () const{ return  fTreeCascVarPosNSigmaProton   ;}
    Float_t GetPosNSigmaPion                            () const{ return  fTreeCascVarPosNSigmaPion   ;}
    Float_t GetNegNSigmaProton                          () const{ return  fTreeCascVarNegNSigmaProton   ;}
    Float_t GetNegNSigmaPion                            () const{ return  fTreeCascVarNegNSigmaPion   ;}
    //---------------PID-TOF-INFO---------------------
    Float_t GetBachTOFNSigmaPion                        () const{ return  fTreeCascVarBachTOFNSigmaPion   ;}
    Float_t GetBachTOFNSigmaKaon                        () const{ return  fTreeCascVarBachTOFNSigmaKaon   ;}
    Float_t GetPosTOFNSigmaProton                       () const{ return  fTreeCascVarPosTOFNSigmaProton   ;}
    Float_t GetPosTOFNSigmaPion                         () const{ return  fTreeCascVarPosTOFNSigmaPion   ;}
    Float_t GetNegTOFNSigmaProton                       () const{ return  fTreeCascVarNegTOFNSigmaProton   ;}
    Float_t GetNegTOFNSigmaPion                         () const{ return  fTreeCascVarNegTOFNSigmaPion   ;}
    //---------------PID-ITS-INFO---------------------
    Float_t GetBachITSNSigmaPion                        () const{ return  fTreeCascVarBachITSNSigmaPion   ;}
    Float_t GetBachITSNSigmaKaon                        () const{ return  fTreeCascVarBachITSNSigmaKaon   ;}
    Float_t GetPosITSNSigmaProton                       () const{ return  fTreeCascVarPosITSNSigmaProton   ;}
    Float_t GetPosITSNSigmaPion                         () const{ return  fTreeCascVarPosITSNSigmaPion   ;}
    Float_t GetNegITSNSigmaProton                       () const{ return  fTreeCascVarNegITSNSigmaProton   ;}
    Float_t GetNegITSNSigmaPion                         () const{ return  fTreeCascVarNegITSNSigmaPion   ;}
    //---Raw TPC dEdx + PIDForTracking information----
    Double_t GetPosdEdx                                 () const{ return  fTreeCascVarPosdEdx   ;}
    Double_t GetNegdEdx                                 () const{ return  fTreeCascVarNegdEdx   ;}
    Double_t GetBachdEdx                                () const{ return  fTreeCascVarBachdEdx   ;}
    Double_t GetPosdEdxN                                () const{ return  fTreeCascVarPosdEdxN   ;}
    Double_t GetNegdEdxN                                () const{ return  fTreeCascVarNegdEdxN   ;}
    Double_t GetBachdEdxN                               () const{ return  fTreeCascVarBachdEdxN   ;}
    Double_t GetPosPIDForTracking                       () const{ return  fTreeCascVarPosPIDForTracking   ;}
    Double_t GetNegPIDForTracking                       () const{ return  fTreeCascVarNegPIDForTracking   ;}
    Double_t GetBachPIDForTracking                      () const{ return  fTreeCascVarBachPIDForTracking   ;}
    //------------------------------------------------
    Double_t GetChi2Cascade                             () const{ return  fTreeCascVarChi2Cascade   ;}
    Double_t GetChi2CascadePerNDF                       () const{ return  fTreeCascVarChi2CascadePerNDF   ;}
    Double_t GetChi2V0                                  () const{ return  fTreeCascVarChi2V0   ;}
    //------------------------------------------------
    ULong64_t GetBachTrackStatus                        () const{ return  fTreeCascVarBachTrackStatus   ;}
    ULong64_t GetPosTrackStatus                         () const{ return  fTreeCascVarPosTrackStatus   ;}
    ULong64_t GetNegTrackStatus                         () const{ return  fTreeCascVarNegTrackStatus   ;}
    //------------------------------------------------
    Float_t GetBachDCAz                                 () const{ return  fTreeCascVarBachDCAz   ;}
    Float_t GetPosDCAz                                  () const{ return  fTreeCascVarPosDCAz   ;}
    Float_t GetNegDCAz                                  () const{ return  fTreeCascVarNegDCAz   ;}
    //------------FULL-MOMENTUM-INFO------------------
    Double_t GetBachPx                                  () const{ return  fTreeCascVarBachPx   ;}
    Double_t GetBachPy                                  () const{ return  fTreeCascVarBachPy   ;}
    Double_t GetBachPz                                  () const{ return  fTreeCascVarBachPz   ;}
    Double_t GetPosPx                                   () const{ return  fTreeCascVarPosPx   ;}
    Double_t GetPosPy                                   () const{ return  fTreeCascVarPosPy   ;}
    Double_t GetPosPz                                   () const{ return  fTreeCascVarPosPz   ;}
    Double_t GetNegPx                                   () const{ return  fTreeCascVarNegPx   ;}
    Double_t GetNegPy                                   () const{ return  fTreeCascVarNegPy   ;}
    Double_t GetNegPz                                   () const{ return  fTreeCascVarNegPz   ;}
    //------------------------------------------------
    Double_t GetCascadeDecayX                           () const{ return  fTreeCascVarCascadeDecayX   ;}
    Double_t GetCascadeDecayY                           () const{ return  fTreeCascVarCascadeDecayY   ;}
    Double_t GetCascadeDecayZ                           () const{ return  fTreeCascVarCascadeDecayZ   ;}
    Double_t GetV0DecayX                                () const{ return  fTreeCascVarV0DecayX   ;}
    Double_t GetV0DecayY                                () const{ return  fTreeCascVarV0DecayY   ;}
    Double_t GetV0DecayZ                                () const{ return  fTreeCascVarV0DecayZ   ;}
    //------------------------------------------------
    Int_t GetBachLabel                                  () const{ return  fTreeCascVarBachLabel  ;}
    Int_t GetPosLabel                                   () const{ return  fTreeCascVarPosLabel   ;}
    Int_t GetNegLabel                                   () const{ return  fTreeCascVarNegLabel   ;}
    //---------------CLUSTER-INFO---------------------
    Bool_t GetBachITSClusters0                          () const{ return  fTreeCascVarBachITSClusters0   ;}
    Bool_t GetBachITSClusters1                          () const{ return  fTreeCascVarBachITSClusters1   ;}
    Bool_t GetBachITSClusters2                          () const{ return  fTreeCascVarBachITSClusters2   ;}
    Bool_t GetBachITSClusters3                          () const{ return  fTreeCascVarBachITSClusters3   ;}
    Bool_t GetBachITSClusters4                          () const{ return  fTreeCascVarBachITSClusters4   ;}
    Bool_t GetBachITSClusters5                          () const{ return  fTreeCascVarBachITSClusters5   ;}
    
    Bool_t GetPosITSClusters0                           () const{ return  fTreeCascVarPosITSClusters0   ;}
    Bool_t GetPosITSClusters1                           () const{ return  fTreeCascVarPosITSClusters1   ;}
    Bool_t GetPosITSClusters2                           () const{ return  fTreeCascVarPosITSClusters2   ;}
    Bool_t GetPosITSClusters3                           () const{ return  fTreeCascVarPosITSClusters3   ;}
    Bool_t GetPosITSClusters4                           () const{ return  fTreeCascVarPosITSClusters4   ;}
    Bool_t GetPosITSClusters5                           () const{ return  fTreeCascVarPosITSClusters5   ;}
    
    Bool_t GetNegITSClusters0                           () const{ return  fTreeCascVarNegITSClusters0   ;}
    Bool_t GetNegITSClusters1                           () const{ return  fTreeCascVarNegITSClusters1   ;}
    Bool_t GetNegITSClusters2                           () const{ return  fTreeCascVarNegITSClusters2   ;}
    Bool_t GetNegITSClusters3                           () const{ return  fTreeCascVarNegITSClusters3   ;}
    Bool_t GetNegITSClusters4                           () const{ return  fTreeCascVarNegITSClusters4   ;}
    Bool_t GetNegITSClusters5                           () const{ return  fTreeCascVarNegITSClusters5   ;}
    
    //------------------------------------------------
    Bool_t GetPosITSSharedClusters0                     () const{ return  fTreeCascVarPosITSSharedClusters0   ;}
    Bool_t GetPosITSSharedClusters1                     () const{ return  fTreeCascVarPosITSSharedClusters1   ;}
    Bool_t GetPosITSSharedClusters2                     () const{ return  fTreeCascVarPosITSSharedClusters2   ;}
    Bool_t GetPosITSSharedClusters3                     () const{ return  fTreeCascVarPosITSSharedClusters3   ;}
    Bool_t GetPosITSSharedClusters4                     () const{ return  fTreeCascVarPosITSSharedClusters4   ;}
    Bool_t GetPosITSSharedClusters5                     () const{ return  fTreeCascVarPosITSSharedClusters5   ;}
    
    Bool_t GetNegITSSharedClusters0                     () const{ return  fTreeCascVarNegITSSharedClusters0   ;}
    Bool_t GetNegITSSharedClusters1                     () const{ return  fTreeCascVarNegITSSharedClusters1   ;}
    Bool_t GetNegITSSharedClusters2                     () const{ return  fTreeCascVarNegITSSharedClusters2   ;}
    Bool_t GetNegITSSharedClusters3                     () const{ return  fTreeCascVarNegITSSharedClusters3   ;}
    Bool_t GetNegITSSharedClusters4                     () const{ return  fTreeCascVarNegITSSharedClusters4   ;}
    Bool_t GetNegITSSharedClusters5                     () const{ return  fTreeCascVarNegITSSharedClusters5   ;}
    
    Bool_t GetBachITSSharedClusters0                    () const{ return  fTreeCascVarBachITSSharedClusters0   ;}
    Bool_t GetBachITSSharedClusters1                    () const{ return  fTreeCascVarBachITSSharedClusters1   ;}
    Bool_t GetBachITSSharedClusters2                    () const{ return  fTreeCascVarBachITSSharedClusters2   ;}
    Bool_t GetBachITSSharedClusters3                    () const{ return  fTreeCascVarBachITSSharedClusters3   ;}
    Bool_t GetBachITSSharedClusters4                    () const{ return  fTreeCascVarBachITSSharedClusters4   ;}
    Bool_t GetBachITSSharedClusters5                    () const{ return  fTreeCascVarBachITSSharedClusters5   ;}

    //---------------OOB-PILEUP-INFO---------------------
    Double_t GetBachTOFExpTDiff                         () const{ return  fTreeCascVarBachTOFExpTDiff   ;}
    Double_t GetPosTOFExpTDiff                          () const{ return  fTreeCascVarPosTOFExpTDiff   ;}
    Double_t GetNegTOFExpTDiff                          () const{ return  fTreeCascVarNegTOFExpTDiff   ;}
    
    Double_t GetBachTOFSignal                           () const{ return  fTreeCascVarBachTOFSignal   ;}
    Double_t GetPosTOFSignal                            () const{ return  fTreeCascVarPosTOFSignal   ;}
    Double_t GetNegTOFSignal                            () const{ return  fTreeCascVarNegTOFSignal   ;}
    
    Int_t   GetBachTOFBCid                              () const{ return  fTreeCascVarBachTOFBCid   ;}
    Int_t   GetPosTOFBCid                               () const{ return  fTreeCascVarPosTOFBCid   ;}
    Int_t   GetNegTOFBCid                               () const{ return  fTreeCascVarNegTOFBCid   ;}
    
    Double_t GetBachTOFLength                           () const{ return  fTreeCascVarBachTOFLength   ;}
    Double_t GetPosTOFLength                            () const{ return  fTreeCascVarPosTOFLength   ;}
    Double_t GetNegTOFLength                            () const{ return  fTreeCascVarNegTOFLength   ;}
    
    Double_t GetBachTOFDeltaX                           () const{ return  fTreeCascVarBachTOFDeltaX   ;}
    Double_t GetPosTOFDeltaX                            () const{ return  fTreeCascVarPosTOFDeltaX   ;}
    Double_t GetNegTOFDeltaX                            () const{ return  fTreeCascVarNegTOFDeltaX   ;}
    
    Double_t GetBachTOFDeltaZ                           () const{ return  fTreeCascVarBachTOFDeltaZ   ;}
    Double_t GetPosTOFDeltaZ                            () const{ return  fTreeCascVarPosTOFDeltaZ   ;}
    Double_t GetNegTOFDeltaZ                            () const{ return  fTreeCascVarNegTOFDeltaZ   ;}
    
    //Kink tagging
    Bool_t GetBachIsKink                                () const{ return  fTreeCascVarBachIsKink   ;}
    Bool_t GetPosIsKink                                 () const{ return  fTreeCascVarPosIsKink   ;}
    Bool_t GetNegIsKink                                 () const{ return  fTreeCascVarNegIsKink   ;}
    
    //Cowboy/sailor studies
    Bool_t  GetIsCowboy                                 () const{ return  fTreeCascVarIsCowboy   ;} 
    Double_t GetCowboyness                              () const{ return  fTreeCascVarCowboyness   ;} 
    Bool_t  GetIsCascadeCowboy                          () const{ return  fTreeCascVarIsCascadeCowboy   ;} 
    Double_t GetCascadeCowboyness                       () const{ return  fTreeCascVarCascadeCowboyness   ;} 
    
    //------------MONTE-CARLO-INFO------------------
    Double_t GetBachEtaMC                               () const{ return  fTreeCascVarBachEtaMC   ;}
    Double_t GetPosEtaMC                                () const{ return  fTreeCascVarPosEtaMC   ;}
    Double_t GetNegEtaMC                                () const{ return  fTreeCascVarNegEtaMC   ;}
    
    Double_t GetBachPxMC                                () const{ return  fTreeCascVarBachPxMC   ;}
    Double_t GetBachPyMC                                () const{ return  fTreeCascVarBachPyMC   ;}
    Double_t GetBachPzMC                                () const{ return  fTreeCascVarBachPzMC   ;}
    Double_t GetPosPxMC                                 () const{ return  fTreeCascVarPosPxMC   ;}
    Double_t GetPosPyMC                                 () const{ return  fTreeCascVarPosPyMC   ;}
    Double_t GetPosPzMC                                 () const{ return  fTreeCascVarPosPzMC   ;}
    Double_t GetNegPxMC                                 () const{ return  fTreeCascVarNegPxMC   ;}
    Double_t GetNegPyMC                                 () const{ return  fTreeCascVarNegPyMC   ;}
    Double_t GetNegPzMC                                 () const{ return  fTreeCascVarNegPzMC   ;}
    
    Int_t GetBachMotherLabel                            () const{ return  fTreeCascVarBachMotherLabel   ;}
    Int_t GetPosMotherLabel                             () const{ return  fTreeCascVarPosMotherLabel   ;}
    Int_t GetNegMotherLabel                             () const{ return  fTreeCascVarNegMotherLabel   ;}
    
    Int_t GetBachGrandMotherLabel                       () const{ return  fTreeCascVarBachGrandMotherLabel   ;}
    Int_t GetPosGrandMotherLabel                        () const{ return  fTreeCascVarPosGrandMotherLabel   ;}
    Int_t GetNegGrandMotherLabel                        () const{ return  fTreeCascVarNegGrandMotherLabel   ;}
    
    Int_t GetBachPdgCode                                () const{ return  fTreeCascVarBachPdgCode   ;}
    Int_t GetPosPdgCode                                 () const{ return  fTreeCascVarPosPdgCode   ;}
    Int_t GetNegPdgCode                                 () const{ return  fTreeCascVarNegPdgCode   ;}
    
    Int_t GetBachMotherPdgCode                          () const{ return  fTreeCascVarBachMotherPdgCode   ;}
    Int_t GetPosMotherPdgCode                           () const{ return  fTreeCascVarPosMotherPdgCode   ;}
    Int_t GetNegMotherPdgCode                           () const{ return  fTreeCascVarNegMotherPdgCode   ;}
    
    Int_t GetBachGrandMotherPdgCode                     () const{ return  fTreeCascVarBachGrandMotherPdgCode   ;}
    Int_t GetPosGrandMotherPdgCode                      () const{ return  fTreeCascVarPosGrandMotherPdgCode   ;}
    Int_t GetNegGrandMotherPdgCode                      () const{ return  fTreeCascVarNegGrandMotherPdgCode   ;}
    
    Int_t GetPdgCode                                    () const{ return  fTreeCascVarPdgCode   ;}
    Double_t GetRapMC                                   () const{ return  fTreeCascVarRapMC   ;}
    Double_t GetPtotMC                                  () const{ return  fTreeCascVarPtotMC   ;}
    Double_t GetPtMC                                    () const{ return  fTreeCascVarPtMC   ;}
    
    Int_t GetV0PdgCode                                  () const{ return  fTreeCascVarV0PdgCode   ;}
    Double_t GetV0PtotMC                                () const{ return  fTreeCascVarV0PtotMC   ;}
    Double_t GetV0PtMC                                  () const{ return  fTreeCascVarV0PtMC   ;}
    
    Double_t GetCascadeDecayXMC                         () const{ return  fTreeCascVarCascadeDecayXMC   ;}
    Double_t GetCascadeDecayYMC                         () const{ return  fTreeCascVarCascadeDecayYMC   ;}
    Double_t GetCascadeDecayZMC                         () const{ return  fTreeCascVarCascadeDecayZMC   ;}
    
    Double_t GetV0DecayXMC                              () const{ return  fTreeCascVarV0DecayXMC   ;}
    Double_t GetV0DecayYMC                              () const{ return  fTreeCascVarV0DecayYMC   ;}
    Double_t GetV0DecayZMC                              () const{ return  fTreeCascVarV0DecayZMC   ;}

    Int_t GetIsPhysicalPrimary                          () const{ return  fTreeCascVarIsPhysicalPrimary   ;}
    
    Bool_t GetIsPhysicalPrimaryBach                     () const{ return  fTreeCascVarIsPhysicalPrimaryBach   ;}
    Bool_t GetIsPhysicalPrimaryPos                      () const{ return  fTreeCascVarIsPhysicalPrimaryPos   ;}
    Bool_t GetIsPhysicalPrimaryNeg                      () const{ return  fTreeCascVarIsPhysicalPrimaryNeg   ;}
    
    Bool_t GetIsPhysicalPrimaryBachMother               () const{ return  fTreeCascVarIsPhysicalPrimaryBachMother   ;}
    Bool_t GetIsPhysicalPrimaryPosMother                () const{ return  fTreeCascVarIsPhysicalPrimaryPosMother   ;}
    Bool_t GetIsPhysicalPrimaryNegMother                () const{ return  fTreeCascVarIsPhysicalPrimaryNegMother   ;}
    
    Bool_t GetIsPhysicalPrimaryBachGrandMother          () const{ return  fTreeCascVarIsPhysicalPrimaryBachGrandMother   ;}
    Bool_t GetIsPhysicalPrimaryPosGrandMother           () const{ return  fTreeCascVarIsPhysicalPrimaryPosGrandMother   ;}
    Bool_t GetIsPhysicalPrimaryNegGrandMother           () const{ return  fTreeCascVarIsPhysicalPrimaryNegGrandMother   ;}
    
private:
    //-----------BASIC-INFO---------------------------
    Int_t fTreeCascVarCharge;//
    Double_t fTreeCascVarMassAsXi;//
    Double_t fTreeCascVarMassAsOmega;//
    Double_t fTreeCascVarPtot;//
    Double_t fTreeCascVarV0Ptot;//
    Double_t fTreeCascVarPt;//
    Double_t fTreeCascVarV0Pt;//
    Double_t fTreeCascVarRapXi;//
    Double_t fTreeCascVarRapOmega;//
    Double_t fTreeCascVarCascEta;//
    Double_t fTreeCascVarV0Eta;//
    Double_t fTreeCascVarBachEta;//
    Double_t fTreeCascVarPosEta;//
    Double_t fTreeCascVarNegEta;//
    //-----------INFO-FOR-CUTS------------------------
    Bool_t  fTreeCascVarOnFlyStatus;//
    Double_t fTreeCascVarAlpha;//
    Double_t fTreeCascVarPtArm;//
    Double_t fTreeCascVarAlphaV0;//
    Double_t fTreeCascVarPtArmV0;//
    Double_t fTreeCascVarDCACascDau;//
    Float_t fTreeCascVarDCABachToPV;//
    Double_t fTreeCascVarDCAV0Dau;//
    Double_t fTreeCascVarDCAV0ToPV;//
    Float_t fTreeCascVarDCAPosToPV;//
    Float_t fTreeCascVarDCANegToPV;//
    Double_t fTreeCascVarCascCosPA;//
    Double_t fTreeCascVarCascCosPASpecial;//
    
    Double_t fTreeCascVarCascRadius;//
    Double_t fTreeCascVarV0Mass;//
    Double_t fTreeCascVarV0MassAsLambda;//
    Double_t fTreeCascVarV0MassAsAntiLambda;//
    Double_t fTreeCascVarV0Radius;//
    Double_t fTreeCascVarV0CosPA;//
    Double_t fTreeCascVarWrongCosPA;//
    Double_t fTreeCascVarDCABachToBaryon;//
    Int_t fTreeCascVarLeastNbrCrossedRows;//
    Int_t fTreeCascVarLeastNbrClusters;//
    Double_t fTreeCascVarLeastRatioCrossedRowsOverFindable;//
    Double_t fTreeCascVarNbrCrossedRowsOverLength;//
    Double_t fTreeCascVarMaxChi2PerCluster;//
    Double_t fTreeCascVarMinTrackLength;  //  
    //-----------DECAY-LENGTH-INFO--------------------
    Double_t fTreeCascVarDistOverTotMom;//
    Double_t fTreeCascVarV0DistOverTotMom;//
    //------------------------------------------------
    //---------------PID-TPC-INFO---------------------
    Float_t fTreeCascVarBachNSigmaPion;//
    Float_t fTreeCascVarBachNSigmaKaon;//
    Float_t fTreeCascVarPosNSigmaProton;//
    Float_t fTreeCascVarPosNSigmaPion;//
    Float_t fTreeCascVarNegNSigmaProton;//
    Float_t fTreeCascVarNegNSigmaPion;//
    //---------------PID-TOF-INFO---------------------
    Float_t fTreeCascVarBachTOFNSigmaPion;//
    Float_t fTreeCascVarBachTOFNSigmaKaon;//
    Float_t fTreeCascVarPosTOFNSigmaProton;//
    Float_t fTreeCascVarPosTOFNSigmaPion;//
    Float_t fTreeCascVarNegTOFNSigmaProton;//
    Float_t fTreeCascVarNegTOFNSigmaPion;//
    //---------------PID-ITS-INFO---------------------
    Float_t fTreeCascVarBachITSNSigmaPion;//
    Float_t fTreeCascVarBachITSNSigmaKaon;//
    Float_t fTreeCascVarPosITSNSigmaProton;//
    Float_t fTreeCascVarPosITSNSigmaPion;//
    Float_t fTreeCascVarNegITSNSigmaProton;//
    Float_t fTreeCascVarNegITSNSigmaPion;//
    //---Raw TPC dEdx + PIDForTracking information----
    Double_t fTreeCascVarPosdEdx;//
    Double_t fTreeCascVarNegdEdx;//
    Double_t fTreeCascVarBachdEdx;//
    Double_t fTreeCascVarPosdEdxN;//
    Double_t fTreeCascVarNegdEdxN;//
    Double_t fTreeCascVarBachdEdxN;//
    Double_t fTreeCascVarPosPIDForTracking;//
    Double_t fTreeCascVarNegPIDForTracking;//
    Double_t fTreeCascVarBachPIDForTracking;//
    //------------------------------------------------
    Double_t fTreeCascVarChi2Cascade;//
    Double_t fTreeCascVarChi2CascadePerNDF;//
    Double_t fTreeCascVarChi2V0;//
    //------------------------------------------------
    ULong64_t fTreeCascVarBachTrackStatus;//
    ULong64_t fTreeCascVarPosTrackStatus; //
    ULong64_t fTreeCascVarNegTrackStatus; //
    //------------------------------------------------
    Float_t fTreeCascVarBachDCAz; //
    Float_t fTreeCascVarPosDCAz; //
    Float_t fTreeCascVarNegDCAz; //
    //------------FULL-MOMENTUM-INFO------------------
    Double_t fTreeCascVarBachPx; //
    Double_t fTreeCascVarBachPy; //
    Double_t fTreeCascVarBachPz; //
    Double_t fTreeCascVarPosPx; //
    Double_t fTreeCascVarPosPy; //
    Double_t fTreeCascVarPosPz; //
    Double_t fTreeCascVarNegPx; //
    Double_t fTreeCascVarNegPy; //
    Double_t fTreeCascVarNegPz; //
    //------------------------------------------------
    Double_t fTreeCascVarCascadeDecayX; //
    Double_t fTreeCascVarCascadeDecayY; //
    Double_t fTreeCascVarCascadeDecayZ; //
    Double_t fTreeCascVarV0DecayX; //
    Double_t fTreeCascVarV0DecayY; //
    Double_t fTreeCascVarV0DecayZ; //
    //------------------------------------------------
    Int_t fTreeCascVarBachLabel; //
    Int_t fTreeCascVarPosLabel; //
    Int_t fTreeCascVarNegLabel; //
    //---------------CLUSTER-INFO---------------------
    Bool_t fTreeCascVarBachITSClusters0;//
    Bool_t fTreeCascVarBachITSClusters1;//
    Bool_t fTreeCascVarBachITSClusters2;//
    Bool_t fTreeCascVarBachITSClusters3;//
    Bool_t fTreeCascVarBachITSClusters4;//
    Bool_t fTreeCascVarBachITSClusters5;//
    
    Bool_t fTreeCascVarPosITSClusters0;//
    Bool_t fTreeCascVarPosITSClusters1;//
    Bool_t fTreeCascVarPosITSClusters2;//
    Bool_t fTreeCascVarPosITSClusters3;//
    Bool_t fTreeCascVarPosITSClusters4;//
    Bool_t fTreeCascVarPosITSClusters5;//
    
    Bool_t fTreeCascVarNegITSClusters0;//
    Bool_t fTreeCascVarNegITSClusters1;//
    Bool_t fTreeCascVarNegITSClusters2;//
    Bool_t fTreeCascVarNegITSClusters3;//
    Bool_t fTreeCascVarNegITSClusters4;//
    Bool_t fTreeCascVarNegITSClusters5;//
    
    //------------------------------------------------
    Bool_t fTreeCascVarPosITSSharedClusters0;//
    Bool_t fTreeCascVarPosITSSharedClusters1;//
    Bool_t fTreeCascVarPosITSSharedClusters2;//
    Bool_t fTreeCascVarPosITSSharedClusters3;//
    Bool_t fTreeCascVarPosITSSharedClusters4;//
    Bool_t fTreeCascVarPosITSSharedClusters5;//
    
    Bool_t fTreeCascVarNegITSSharedClusters0;//
    Bool_t fTreeCascVarNegITSSharedClusters1;//
    Bool_t fTreeCascVarNegITSSharedClusters2;//
    Bool_t fTreeCascVarNegITSSharedClusters3;//
    Bool_t fTreeCascVarNegITSSharedClusters4;//
    Bool_t fTreeCascVarNegITSSharedClusters5;//
    
    Bool_t fTreeCascVarBachITSSharedClusters0;//
    Bool_t fTreeCascVarBachITSSharedClusters1;//
    Bool_t fTreeCascVarBachITSSharedClusters2;//
    Bool_t fTreeCascVarBachITSSharedClusters3;//
    Bool_t fTreeCascVarBachITSSharedClusters4;//
    Bool_t fTreeCascVarBachITSSharedClusters5;//

    //---------------OOB-PILEUP-INFO---------------------
    Double_t fTreeCascVarBachTOFExpTDiff; //
    Double_t fTreeCascVarPosTOFExpTDiff; //
    Double_t fTreeCascVarNegTOFExpTDiff; //
    
    Double_t fTreeCascVarBachTOFSignal; //
    Double_t fTreeCascVarPosTOFSignal; //
    Double_t fTreeCascVarNegTOFSignal; //
    
    Int_t   fTreeCascVarBachTOFBCid; //
    Int_t   fTreeCascVarPosTOFBCid; //
    Int_t   fTreeCascVarNegTOFBCid; //
    
    Double_t fTreeCascVarBachTOFLength;//
    Double_t fTreeCascVarPosTOFLength;//
    Double_t fTreeCascVarNegTOFLength;//
    
    Double_t fTreeCascVarBachTOFDeltaX;//
    Double_t fTreeCascVarPosTOFDeltaX;//
    Double_t fTreeCascVarNegTOFDeltaX;//
    
    Double_t fTreeCascVarBachTOFDeltaZ;//
    Double_t fTreeCascVarPosTOFDeltaZ;//
    Double_t fTreeCascVarNegTOFDeltaZ;//
    
    //Kink tagging
    Bool_t fTreeCascVarBachIsKink;//
    Bool_t fTreeCascVarPosIsKink;//
    Bool_t fTreeCascVarNegIsKink;//
    
    //Cowboy/sailor studies
    Bool_t  fTreeCascVarIsCowboy;   // store if V0 is cowboy-like or sailor-like in XY plane
    Double_t fTreeCascVarCowboyness; // negative -> cowboy, positive -> sailor
    Bool_t  fTreeCascVarIsCascadeCowboy;   // store if V0 is cowboy-like or sailor-like in XY plane
    Double_t fTreeCascVarCascadeCowboyness; // negative -> cowboy, positive -> sailor
    
    //------------MONTE-CARLO-INFO------------------
    Double_t fTreeCascVarBachEtaMC;//
    Double_t fTreeCascVarPosEtaMC;//
    Double_t fTreeCascVarNegEtaMC;//
    
    Double_t fTreeCascVarBachPxMC; //
    Double_t fTreeCascVarBachPyMC; //
    Double_t fTreeCascVarBachPzMC; //
    Double_t fTreeCascVarPosPxMC; //
    Double_t fTreeCascVarPosPyMC; //
    Double_t fTreeCascVarPosPzMC; //
    Double_t fTreeCascVarNegPxMC; //
    Double_t fTreeCascVarNegPyMC; //
    Double_t fTreeCascVarNegPzMC; //
    
    Int_t fTreeCascVarBachMotherLabel;//
    Int_t fTreeCascVarPosMotherLabel;//
    Int_t fTreeCascVarNegMotherLabel;//
    
    Int_t fTreeCascVarBachGrandMotherLabel;//
    Int_t fTreeCascVarPosGrandMotherLabel;//
    Int_t fTreeCascVarNegGrandMotherLabel;//
    
    Int_t fTreeCascVarBachPdgCode;//
    Int_t fTreeCascVarPosPdgCode;//
    Int_t fTreeCascVarNegPdgCode;//
    
    Int_t fTreeCascVarBachMotherPdgCode;//
    Int_t fTreeCascVarPosMotherPdgCode;//
    Int_t fTreeCascVarNegMotherPdgCode;//
    
    Int_t fTreeCascVarBachGrandMotherPdgCode;//
    Int_t fTreeCascVarPosGrandMotherPdgCode;//
    Int_t fTreeCascVarNegGrandMotherPdgCode;//
    
    Int_t fTreeCascVarPdgCode;//
    Double_t fTreeCascVarRapMC;//
    Double_t fTreeCascVarPtotMC;//
    Double_t fTreeCascVarPtMC;//
    
    Int_t fTreeCascVarV0PdgCode;//
    Double_t fTreeCascVarV0PtotMC;//
    Double_t fTreeCascVarV0PtMC;//
    
    Double_t fTreeCascVarCascadeDecayXMC;//
    Double_t fTreeCascVarCascadeDecayYMC;//
    Double_t fTreeCascVarCascadeDecayZMC;//
    
    Double_t fTreeCascVarV0DecayXMC;//
    Double_t fTreeCascVarV0DecayYMC;//
    Double_t fTreeCascVarV0DecayZMC;//

    Int_t fTreeCascVarIsPhysicalPrimary;//
    
    Bool_t fTreeCascVarIsPhysicalPrimaryBach;//
    Bool_t fTreeCascVarIsPhysicalPrimaryPos;//
    Bool_t fTreeCascVarIsPhysicalPrimaryNeg;//
    
    Bool_t fTreeCascVarIsPhysicalPrimaryBachMother;//
    Bool_t fTreeCascVarIsPhysicalPrimaryPosMother;//
    Bool_t fTreeCascVarIsPhysicalPrimaryNegMother;//
    
    Bool_t fTreeCascVarIsPhysicalPrimaryBachGrandMother;//
    Bool_t fTreeCascVarIsPhysicalPrimaryPosGrandMother;//
    Bool_t fTreeCascVarIsPhysicalPrimaryNegGrandMother;//
    
    //------------------------------------------------
    ClassDef(AliMCCascadeContainer, 1);
};
#endif
