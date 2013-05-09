#ifndef ALIEMCALRECOUTILS_H
#define ALIEMCALRECOUTILS_H

/* $Id: AliEMCALRecoUtils.h 33808 2009-07-15 09:48:08Z gconesab $ */

///////////////////////////////////////////////////////////////////////////////
//
// Class AliEMCALRecoUtils
// Some utilities to recalculate the cluster position or energy linearity
//
//
// Author:  Gustavo Conesa (LPSC- Grenoble) 
//          Track matching part: Rongrong Ma (Yale)
///////////////////////////////////////////////////////////////////////////////

//Root includes
#include <TNamed.h>
#include <TMath.h>
class TObjArray;
class TArrayI;
class TArrayF;
#include <TH2I.h>
class TH2F;
#include <TRandom3.h>

//AliRoot includes
class AliVCluster;
class AliVCaloCells;
class AliVEvent;
#include "AliLog.h"

// EMCAL includes
class AliEMCALGeometry;
class AliEMCALPIDUtils;
class AliESDtrack;
class AliExternalTrackParam;

class AliEMCALRecoUtils : public TNamed {
  
public:
  
  AliEMCALRecoUtils();
  AliEMCALRecoUtils(const AliEMCALRecoUtils&); 
  AliEMCALRecoUtils& operator=(const AliEMCALRecoUtils&); 
  virtual ~AliEMCALRecoUtils() ;  
  
  void     InitParameters();
  
  void     Print(const Option_t*) const;

  //enums
  enum     NonlinearityFunctions{kPi0MC=0,kPi0GammaGamma=1,kPi0GammaConversion=2,kNoCorrection=3,kBeamTest=4,kBeamTestCorrected=5,kPi0MCv2=6,kPi0MCv3=7};
  enum     PositionAlgorithms{kUnchanged=-1,kPosTowerIndex=0, kPosTowerGlobal=1};
  enum     ParticleType{kPhoton=0, kElectron=1,kHadron =2, kUnknown=-1};
  enum     { kNCuts = 12 }; //track matching Marcel
  enum     TrackCutsType{kTPCOnlyCut=0, kGlobalCut=1, kLooseCut=2, kITSStandAlone=3};  //Marcel

  //-----------------------------------------------------
  //Position recalculation
  //-----------------------------------------------------

  void     RecalculateClusterPosition               (const AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu); 
  void     RecalculateClusterPositionFromTowerIndex (const AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu); 
  void     RecalculateClusterPositionFromTowerGlobal(const AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu); 
  
  Float_t  GetCellWeight(const Float_t eCell, const Float_t eCluster) const { if (eCell > 0 && eCluster > 0) return TMath::Max( 0., fW0 + TMath::Log( eCell / eCluster )) ;
                                                                              else                           return 0.                                                    ; }
  
  Float_t  GetDepth(const Float_t eCluster, const Int_t iParticle, const Int_t iSM) const ; 
  
  void     GetMaxEnergyCell(const AliEMCALGeometry *geom, AliVCaloCells* cells, const AliVCluster* clu, 
                            Int_t & absId,  Int_t& iSupMod, Int_t& ieta, Int_t& iphi, Bool_t &shared);
  
  Float_t  GetMisalTransShift(const Int_t i)       const { if(i < 15 ) { return fMisalTransShift[i] ; }
                                                           else        { AliInfo(Form("Index %d larger than 15, do nothing\n",i)) ; 
                                                                         return 0.                  ; } }
  Float_t* GetMisalTransShiftArray()                     { return fMisalTransShift ; }

  void     SetMisalTransShift(const Int_t i, const Float_t shift) {
                                                           if(i < 15 ) { fMisalTransShift[i] = shift ; }
                                                           else        { AliInfo(Form("Index %d larger than 15, do nothing\n",i)) ; } }
  void     SetMisalTransShiftArray(Float_t * misal)               { for(Int_t i = 0; i < 15; i++) fMisalTransShift[i] = misal[i]  ; }

  Float_t  GetMisalRotShift(const Int_t i)         const { if(i < 15 ) { return fMisalRotShift[i]    ; }
                                                           else        { AliInfo(Form("Index %d larger than 15, do nothing\n",i)) ; 
                                                                         return 0.                   ; } }
  
  Float_t* GetMisalRotShiftArray()                       { return fMisalRotShift                     ; }
  
  void     SetMisalRotShift(const Int_t i, const Float_t shift) {
                                                           if(i < 15 ) { fMisalRotShift[i] = shift   ; }
                                                           else        { AliInfo(Form("Index %d larger than 15, do nothing\n",i)) ; } }
  
  void     SetMisalRotShiftArray(Float_t * misal)        { for(Int_t i = 0; i < 15; i++)fMisalRotShift[i] = misal[i] ; }
  
  Int_t    GetParticleType()                       const { return  fParticleType    ; }
  void     SetParticleType(Int_t particle)               { fParticleType = particle ; }
  
  Int_t    GetPositionAlgorithm()                  const { return fPosAlgo          ; }
  void     SetPositionAlgorithm(Int_t alg)               { fPosAlgo = alg           ; }
  
  Float_t  GetW0()                                 const { return fW0               ; }
  void     SetW0(Float_t w0)                             { fW0  = w0                ; }

  //-----------------------------------------------------
  // Non Linearity
  //-----------------------------------------------------

  Float_t  CorrectClusterEnergyLinearity(AliVCluster* clu) ;
  
  Float_t  GetNonLinearityParam(const Int_t i)     const { if(i < 7 ){ return fNonLinearityParams[i] ; }
                                                          else      { AliInfo(Form("Index %d larger than 7, do nothing\n",i)) ; 
                                                                       return 0.                     ; } }
  void     SetNonLinearityParam(const Int_t i, const Float_t param) {
                                                          if(i < 7 ){fNonLinearityParams[i] = param ; }
                                                          else { AliInfo(Form("Index %d larger than 7, do nothing\n",i)) ; } }
  void     InitNonLinearityParam();

  Int_t    GetNonLinearityFunction() const               { return fNonLinearityFunction    ; }
  void     SetNonLinearityFunction(Int_t fun)            { fNonLinearityFunction = fun     ; InitNonLinearityParam() ; }

  void     SetNonLinearityThreshold(Int_t threshold)     { fNonLinearThreshold = threshold ; } //only for Alexie's non linearity correction
  Int_t    GetNonLinearityThreshold()              const { return fNonLinearThreshold      ; }
//  
  //-----------------------------------------------------
  // MC clusters energy smearing
  //-----------------------------------------------------
  
  Float_t  SmearClusterEnergy(const AliVCluster* clu) ;
  void     SwitchOnClusterEnergySmearing()               { fSmearClusterEnergy = kTRUE         ; }
  void     SwitchOffClusterEnergySmearing()              { fSmearClusterEnergy = kFALSE        ; }
  Bool_t   IsClusterEnergySmeared()                const { return fSmearClusterEnergy          ; }   
  void     SetSmearingParameters(Int_t i, Float_t param) { if(i < 3){ fSmearClusterParam[i] = param ; }
                                                           else     { AliInfo(Form("Index %d larger than 2, do nothing\n",i)) ; } }
  //-----------------------------------------------------
  // Recalibration
  //-----------------------------------------------------
  Bool_t   AcceptCalibrateCell(const Int_t absId, const Int_t bc,
                               Float_t & amp, Double_t & time, AliVCaloCells* cells) ; // Energy and Time
  void     RecalibrateCells(AliVCaloCells * cells, Int_t bc) ; // Energy and Time
  void     RecalibrateClusterEnergy(const AliEMCALGeometry* geom, AliVCluster* cluster, AliVCaloCells * cells, const Int_t bc=-1) ; // Energy and time
  void     ResetCellsCalibrated()                        { fCellsRecalibrated = kFALSE; }

  // Energy recalibration
  Bool_t   IsRecalibrationOn()                     const { return fRecalibration ; }
  void     SwitchOffRecalibration()                      { fRecalibration = kFALSE ; }
  void     SwitchOnRecalibration()                       { fRecalibration = kTRUE  ; 
                                                           if(!fEMCALRecalibrationFactors)InitEMCALRecalibrationFactors() ; }
  void     InitEMCALRecalibrationFactors() ;
  TObjArray* GetEMCALRecalibrationFactorsArray()   const { return fEMCALRecalibrationFactors ; }

  TH2F *   GetEMCALChannelRecalibrationFactors(Int_t iSM)     const { return (TH2F*)fEMCALRecalibrationFactors->At(iSM) ; }	
  void     SetEMCALChannelRecalibrationFactors(TObjArray *map)      { fEMCALRecalibrationFactors = map                  ; }
  void     SetEMCALChannelRecalibrationFactors(Int_t iSM , TH2F* h) { fEMCALRecalibrationFactors->AddAt(h,iSM)          ; }
  
  Float_t  GetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow) const { 
    if(fEMCALRecalibrationFactors) 
      return (Float_t) ((TH2F*)fEMCALRecalibrationFactors->At(iSM))->GetBinContent(iCol,iRow); 
    else return 1 ; } 
	
  void     SetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
    if(!fEMCALRecalibrationFactors) InitEMCALRecalibrationFactors() ;
    ((TH2F*)fEMCALRecalibrationFactors->At(iSM))->SetBinContent(iCol,iRow,c) ; }
  
  //Recalibrate channels energy with run dependent corrections
  Bool_t   IsRunDepRecalibrationOn()               const { return fUseRunCorrectionFactors ; }

  void     SwitchOffRunDepCorrection()                   { fUseRunCorrectionFactors = kFALSE ; }
  void     SwitchOnRunDepCorrection()                    { fUseRunCorrectionFactors = kTRUE  ; 
                                                           SwitchOnRecalibration()           ; }      
  // Time Recalibration  
  void     RecalibrateCellTime(const Int_t absId, const Int_t bc, Double_t & time) const;
  
  Bool_t   IsTimeRecalibrationOn()                 const { return fTimeRecalibration   ; }
  void     SwitchOffTimeRecalibration()                  { fTimeRecalibration = kFALSE ; }
  void     SwitchOnTimeRecalibration()                   { fTimeRecalibration = kTRUE  ; 
    if(!fEMCALTimeRecalibrationFactors)InitEMCALTimeRecalibrationFactors() ; }
  void     InitEMCALTimeRecalibrationFactors() ;
  TObjArray* GetEMCALTimeRecalibrationFactorsArray() const { return fEMCALTimeRecalibrationFactors ; }

  Float_t  GetEMCALChannelTimeRecalibrationFactor(const Int_t bc, const Int_t absID) const { 
    if(fEMCALTimeRecalibrationFactors) 
      return (Float_t) ((TH1F*)fEMCALTimeRecalibrationFactors->At(bc))->GetBinContent(absID); 
    else return 0 ; } 
	
  void     SetEMCALChannelTimeRecalibrationFactor(const Int_t bc, const Int_t absID, Double_t c = 0) { 
    if(!fEMCALTimeRecalibrationFactors) InitEMCALTimeRecalibrationFactors() ;
    ((TH1F*)fEMCALTimeRecalibrationFactors->At(bc))->SetBinContent(absID,c) ; }  
  
  TH1F *   GetEMCALChannelTimeRecalibrationFactors(const Int_t bc)const       { return (TH1F*)fEMCALTimeRecalibrationFactors->At(bc) ; }	
  void     SetEMCALChannelTimeRecalibrationFactors(TObjArray *map)            { fEMCALTimeRecalibrationFactors = map                 ; }
  void     SetEMCALChannelTimeRecalibrationFactors(const Int_t bc , TH1F* h)  { fEMCALTimeRecalibrationFactors->AddAt(h,bc)          ; }
  
  //-----------------------------------------------------
  // Modules fiducial region, remove clusters in borders
  //-----------------------------------------------------

  Bool_t   CheckCellFiducialRegion(const AliEMCALGeometry* geom, 
                                   const AliVCluster* cluster, 
                                   AliVCaloCells* cells) ;
  void     SetNumberOfCellsFromEMCALBorder(const Int_t n){ fNCellsFromEMCALBorder = n      ; }
  Int_t    GetNumberOfCellsFromEMCALBorder()      const  { return fNCellsFromEMCALBorder   ; }
    
  void     SwitchOnNoFiducialBorderInEMCALEta0()         { fNoEMCALBorderAtEta0 = kTRUE    ; }
  void     SwitchOffNoFiducialBorderInEMCALEta0()        { fNoEMCALBorderAtEta0 = kFALSE   ; }
  Bool_t   IsEMCALNoBorderAtEta0()                 const { return fNoEMCALBorderAtEta0     ; }
  
  //-----------------------------------------------------
  // Bad channels
  //-----------------------------------------------------

  Bool_t   IsBadChannelsRemovalSwitchedOn()        const { return fRemoveBadChannels       ; }
  void     SwitchOffBadChannelsRemoval()                 { fRemoveBadChannels = kFALSE     ; }
  void     SwitchOnBadChannelsRemoval ()                 { fRemoveBadChannels = kTRUE ; 
                                                           if(!fEMCALBadChannelMap)InitEMCALBadChannelStatusMap() ; }
	
  Bool_t   IsDistanceToBadChannelRecalculated()    const { return fRecalDistToBadChannels   ; }
  void     SwitchOffDistToBadChannelRecalculation()      { fRecalDistToBadChannels = kFALSE ; }
  void     SwitchOnDistToBadChannelRecalculation()       { fRecalDistToBadChannels = kTRUE  ; 
                                                           if(!fEMCALBadChannelMap)InitEMCALBadChannelStatusMap() ; }
  
  TObjArray* GetEMCALBadChannelStatusMapArray()     const { return fEMCALBadChannelMap ; }
  void     InitEMCALBadChannelStatusMap() ;
	
  Int_t    GetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow) const { 
    if(fEMCALBadChannelMap) return (Int_t) ((TH2I*)fEMCALBadChannelMap->At(iSM))->GetBinContent(iCol,iRow); 
    else return 0;}//Channel is ok by default
	
  void     SetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
                                                           if(!fEMCALBadChannelMap)InitEMCALBadChannelStatusMap()               ;
                                                           ((TH2I*)fEMCALBadChannelMap->At(iSM))->SetBinContent(iCol,iRow,c)    ; }
	
  TH2I *   GetEMCALChannelStatusMap(Int_t iSM)     const { return (TH2I*)fEMCALBadChannelMap->At(iSM) ; }
  void     SetEMCALChannelStatusMap(TObjArray *map)      { fEMCALBadChannelMap = map                  ; }
  void     SetEMCALChannelStatusMap(Int_t iSM , TH2I* h) { fEMCALBadChannelMap->AddAt(h,iSM)          ; }

  Bool_t   ClusterContainsBadChannel(const AliEMCALGeometry* geom, const UShort_t* cellList, const Int_t nCells);
 
  //-----------------------------------------------------
  // Recalculate other cluster parameters
  //-----------------------------------------------------

  void     RecalculateClusterDistanceToBadChannel (const AliEMCALGeometry * geom, AliVCaloCells* cells, AliVCluster * cluster);
  void     RecalculateClusterShowerShapeParameters(const AliEMCALGeometry * geom, AliVCaloCells* cells, AliVCluster * cluster);
  void     RecalculateClusterShowerShapeParameters(const AliEMCALGeometry * geom, AliVCaloCells* cells, AliVCluster * cluster,
                                                   Float_t & l0,   Float_t & l1,   
                                                   Float_t & disp, Float_t & dEta, Float_t & dPhi,
                                                   Float_t & sEta, Float_t & sPhi, Float_t & sEtaPhi);

  void     RecalculateClusterPID(AliVCluster * cluster);

  AliEMCALPIDUtils * GetPIDUtils() { return fPIDUtils;}


  //----------------------------------------------------
  // Track matching
  //----------------------------------------------------

  void     FindMatches(AliVEvent *event, TObjArray * clusterArr=0x0, const AliEMCALGeometry *geom=0x0);
  Int_t    FindMatchedClusterInEvent(const AliESDtrack *track, const AliVEvent *event, 
                                     const AliEMCALGeometry *geom, Float_t &dEta, Float_t &dPhi);
  Int_t    FindMatchedClusterInClusterArr(const AliExternalTrackParam *emcalParam, 
                                          AliExternalTrackParam *trkParam, 
                                          const TObjArray * clusterArr, 
                                          Float_t &dEta, Float_t &dPhi);
  
  static Bool_t ExtrapolateTrackToEMCalSurface(AliExternalTrackParam *trkParam, 
                                               const Double_t emcalR, const Double_t mass, const Double_t step, 
                                               Float_t &eta, Float_t &phi);
  static Bool_t ExtrapolateTrackToPosition(AliExternalTrackParam *trkParam, const Float_t *clsPos, 
                                           const Double_t mass, const Double_t step, 
                                           Float_t &tmpEta, Float_t &tmpPhi);
  static Bool_t ExtrapolateTrackToCluster (AliExternalTrackParam *trkParam, const AliVCluster *cluster, 
                                           const Double_t mass, const Double_t step,
                                           Float_t &tmpEta, Float_t &tmpPhi);
  Bool_t        ExtrapolateTrackToCluster (AliExternalTrackParam *trkParam, const AliVCluster *cluster, 
                                           Float_t &tmpEta, Float_t &tmpPhi);

  UInt_t   FindMatchedPosForCluster(const Int_t clsIndex) const;
  UInt_t   FindMatchedPosForTrack  (const Int_t trkIndex) const;
  
  void     GetMatchedResiduals       (const Int_t clsIndex, Float_t &dEta, Float_t &dPhi);
  void     GetMatchedClusterResiduals(const Int_t trkIndex, Float_t &dEta, Float_t &dPhi);
  Int_t    GetMatchedTrackIndex(Int_t clsIndex);
  Int_t    GetMatchedClusterIndex(Int_t trkIndex);
  
  Bool_t   IsClusterMatched(const Int_t clsIndex)         const;
  Bool_t   IsTrackMatched  (const Int_t trkIndex)         const;

  void     SetClusterMatchedToTrack (const AliVEvent *event);
  void     SetTracksMatchedToCluster(const AliVEvent *event);  

  void     SwitchOnCutEtaPhiSum()                     { fCutEtaPhiSum      = kTRUE    ; 
                                                        fCutEtaPhiSeparate = kFALSE   ; }
  void     SwitchOnCutEtaPhiSeparate()                { fCutEtaPhiSeparate = kTRUE    ;
                                                        fCutEtaPhiSum      = kFALSE   ; }

  Float_t  GetCutR()                            const { return fCutR                  ; }
  Float_t  GetCutEta()                          const { return fCutEta                ; }
  Float_t  GetCutPhi()                          const { return fCutPhi                ; }
  Double_t GetClusterWindow()                   const { return fClusterWindow         ; }
  void     SetCutR(Float_t cutR)                      { fCutR   = cutR                ; }
  void     SetCutEta(Float_t cutEta)                  { fCutEta = cutEta              ; }
  void     SetCutPhi(Float_t cutPhi)                  { fCutPhi = cutPhi              ; }
  void     SetClusterWindow(Double_t window)          { fClusterWindow = window       ; }
  void     SetCutZ(Float_t cutZ)                      { printf("Obsolete fucntion of cutZ=%1.1f\n",cutZ) ; } //Obsolete

  Double_t GetMass()                            const { return fMass                  ; }
  Double_t GetStep()                            const { return fStepCluster           ; }
  Double_t GetStepSurface()                     const { return fStepSurface           ; }
  void     SetMass(Double_t mass)                     { fMass = mass                  ; }
  void     SetStep(Double_t step)                     { fStepSurface = step           ; }
  void     SetStepCluster(Double_t step)              { fStepCluster = step           ; }
 
  void     SetITSTrackSA(Bool_t isITS)                { fITSTrackSA = isITS           ; } //Special Handle of AliExternTrackParam    
  
  // Exotic cells / clusters
  
  Bool_t   IsExoticCell(const Int_t absId, AliVCaloCells* cells, const Int_t bc =-1) ;
  void     SwitchOnRejectExoticCell()                 { fRejectExoticCells = kTRUE     ; }
  void     SwitchOffRejectExoticCell()                { fRejectExoticCells = kFALSE    ; } 
  Bool_t   IsRejectExoticCell()                 const { return fRejectExoticCells      ; }
   
  Float_t  GetECross(const Int_t absID, const Double_t tcell,
                     AliVCaloCells* cells, const Int_t bc);
  
  Float_t  GetExoticCellFractionCut()           const { return fExoticCellFraction     ; }
  Float_t  GetExoticCellDiffTimeCut()           const { return fExoticCellDiffTime     ; }
  Float_t  GetExoticCellMinAmplitudeCut()       const { return fExoticCellMinAmplitude ; }
  
  void     SetExoticCellFractionCut(Float_t f)        { fExoticCellFraction     = f    ; }
  void     SetExoticCellDiffTimeCut(Float_t dt)       { fExoticCellDiffTime     = dt   ; }
  void     SetExoticCellMinAmplitudeCut(Float_t ma)   { fExoticCellMinAmplitude = ma   ; }
  
  Bool_t   IsExoticCluster(const AliVCluster *cluster, AliVCaloCells* cells, const Int_t bc=0) ;
  void     SwitchOnRejectExoticCluster()              { fRejectExoticCluster = kTRUE   ;
                                                        fRejectExoticCells   = kTRUE   ; }
  void     SwitchOffRejectExoticCluster()             { fRejectExoticCluster = kFALSE  ; }
  Bool_t   IsRejectExoticCluster()              const { return fRejectExoticCluster    ; }
  
  //Cluster cut
  Bool_t   IsGoodCluster(AliVCluster *cluster, const AliEMCALGeometry *geom, 
                         AliVCaloCells* cells, const Int_t bc =-1);

  //Track Cuts 
  Bool_t   IsAccepted(AliESDtrack *track);
  void     InitTrackCuts();
  void     SetTrackCutsType(Int_t type)              { fTrackCutsType = type           ; 
                                                       InitTrackCuts()                 ; }
  Int_t    GetTrackCutsType() const                  { return fTrackCutsType; }

  void     SetAODTrackFilterMask( UInt_t mask)       {fAODFilterMask    = mask  ; }
  void     SwitchOnAODHybridTracksMatch()            {fAODHybridTracks  = kTRUE ; }
  void     SwitchOffAODHybridTracksMatch()           {fAODHybridTracks  = kFALSE ; }
  
  // track quality cut setters  
  void     SetMinTrackPt(Double_t pt=0)              { fCutMinTrackPt           = pt   ; }
  void     SetMinNClustersTPC(Int_t min=-1)          { fCutMinNClusterTPC       = min  ; }
  void     SetMinNClustersITS(Int_t min=-1)          { fCutMinNClusterITS       = min  ; }
  void     SetMaxChi2PerClusterTPC(Float_t max=1e10) { fCutMaxChi2PerClusterTPC = max  ; }
  void     SetMaxChi2PerClusterITS(Float_t max=1e10) { fCutMaxChi2PerClusterITS = max  ; }
  void     SetRequireTPCRefit(Bool_t b=kFALSE)       { fCutRequireTPCRefit      = b    ; }
  void     SetRequireITSRefit(Bool_t b=kFALSE)       { fCutRequireITSRefit      = b    ; }
  void     SetAcceptKinkDaughters(Bool_t b=kTRUE)    { fCutAcceptKinkDaughters  = b    ; }
  void     SetMaxDCAToVertexXY(Float_t dist=1e10)    { fCutMaxDCAToVertexXY     = dist ; }
  void     SetMaxDCAToVertexZ(Float_t dist=1e10)     { fCutMaxDCAToVertexZ      = dist ; }
  void     SetDCAToVertex2D(Bool_t b=kFALSE)         { fCutDCAToVertex2D        = b    ; }
  void     SetRequireITSStandAlone(Bool_t b=kFALSE)    {fCutRequireITSStandAlone = b;} //Marcel
  void     SetRequireITSPureStandAlone(Bool_t b=kFALSE){fCutRequireITSpureSA     = b;}
  // getters								
  Double_t GetMinTrackPt()                     const { return fCutMinTrackPt           ; }
  Int_t    GetMinNClusterTPC()                 const { return fCutMinNClusterTPC       ; }
  Int_t    GetMinNClustersITS()                const { return fCutMinNClusterITS       ; }
  Float_t  GetMaxChi2PerClusterTPC()           const { return fCutMaxChi2PerClusterTPC ; }
  Float_t  GetMaxChi2PerClusterITS()           const { return fCutMaxChi2PerClusterITS ; }
  Bool_t   GetRequireTPCRefit()                const { return fCutRequireTPCRefit      ; }
  Bool_t   GetRequireITSRefit()                const { return fCutRequireITSRefit      ; }
  Bool_t   GetAcceptKinkDaughters()            const { return fCutAcceptKinkDaughters  ; }
  Float_t  GetMaxDCAToVertexXY()               const { return fCutMaxDCAToVertexXY     ; }
  Float_t  GetMaxDCAToVertexZ()                const { return fCutMaxDCAToVertexZ      ; }
  Bool_t   GetDCAToVertex2D()                  const { return fCutDCAToVertex2D        ; }
  Bool_t   GetRequireITSStandAlone()           const { return fCutRequireITSStandAlone ; } //Marcel	

private:  
  //Position recalculation
  Float_t    fMisalTransShift[15];       // Shift parameters
  Float_t    fMisalRotShift[15];         // Shift parameters
  Int_t      fParticleType;              // Particle type for depth calculation
  Int_t      fPosAlgo;                   // Position recalculation algorithm
  Float_t    fW0;                        // Weight0
    
  // Non linearity
  Int_t      fNonLinearityFunction;      // Non linearity function choice
  Float_t    fNonLinearityParams[7];     // Parameters for the non linearity function
  Int_t	     fNonLinearThreshold;        // Non linearity threshold value for kBeamTesh non linearity function 
  
  // Energy smearing for MC
  Bool_t     fSmearClusterEnergy;        // Smear cluster energy, to be done only for simulated data to match real data
  Float_t    fSmearClusterParam[3];      // Smearing parameters
  TRandom3   fRandom;                    // Random generator
    
  // Energy Recalibration 
  Bool_t     fCellsRecalibrated;         // Internal bool to check if cells (time/energy) where recalibrated and not recalibrate them when recalculating different things
  Bool_t     fRecalibration;             // Switch on or off the recalibration
  TObjArray* fEMCALRecalibrationFactors; // Array of histograms with map of recalibration factors, EMCAL
    
  // Time Recalibration 
  Bool_t     fTimeRecalibration;             // Switch on or off the time recalibration
  TObjArray* fEMCALTimeRecalibrationFactors; // Array of histograms with map of time recalibration factors, EMCAL
  
  // Recalibrate with run dependent corrections, energy
  Bool_t     fUseRunCorrectionFactors;   // Use Run Dependent Correction
    
  // Bad Channels
  Bool_t     fRemoveBadChannels;         // Check the channel status provided and remove clusters with bad channels
  Bool_t     fRecalDistToBadChannels;    // Calculate distance from highest energy tower of cluster to closes bad channel
  TObjArray* fEMCALBadChannelMap;        // Array of histograms with map of bad channels, EMCAL

  // Border cells
  Int_t      fNCellsFromEMCALBorder;     // Number of cells from EMCAL border the cell with maximum amplitude has to be.
  Bool_t     fNoEMCALBorderAtEta0;       // Do fiducial cut in EMCAL region eta = 0?
  
  // Exotic cell / cluster
  Bool_t     fRejectExoticCluster;       // Switch on or off exotic cluster rejection
  Bool_t     fRejectExoticCells;         // Remove exotic cells
  Float_t    fExoticCellFraction;        // Good cell if fraction < 1-ecross/ecell
  Float_t    fExoticCellDiffTime;        // If time of candidate to exotic and close cell is too different (in ns), it must be noisy, set amp to 0
  Float_t    fExoticCellMinAmplitude;    // Check for exotic only if amplitud is larger than this value
  
  // PID
  AliEMCALPIDUtils * fPIDUtils;          // Recalculate PID parameters
    
  //Track matching
  UInt_t     fAODFilterMask;             // Filter mask to select AOD tracks. Refer to $ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C
  Bool_t     fAODHybridTracks;           // Match with hybrid
  
  TArrayI  * fMatchedTrackIndex;         // Array that stores indexes of matched tracks      
  TArrayI  * fMatchedClusterIndex;       // Array that stores indexes of matched clusters
  TArrayF  * fResidualEta;               // Array that stores the residual eta
  TArrayF  * fResidualPhi;               // Array that stores the residual phi
  Bool_t     fCutEtaPhiSum;              // Place cut on sqrt(dEta^2+dPhi^2)
  Bool_t     fCutEtaPhiSeparate;         // Cut on dEta and dPhi separately
  Float_t    fCutR;                      // sqrt(dEta^2+dPhi^2) cut on matching
  Float_t    fCutEta;                    // dEta cut on matching
  Float_t    fCutPhi;                    // dPhi cut on matching
  Double_t   fClusterWindow;             // Select clusters in the window to be matched
  Double_t   fMass;                      // Mass hypothesis of the track
  Double_t   fStepSurface;               // Length of step to extrapolate tracks to EMCal surface
  Double_t   fStepCluster;               // Length of step to extrapolate tracks to clusters
  Bool_t     fITSTrackSA;                // If track matching is to be done with ITS tracks standing alone	
 
  // Track cuts  
  Int_t      fTrackCutsType;             // Esd track cuts type for matching
  Double_t   fCutMinTrackPt;             // Cut on track pT
  Int_t      fCutMinNClusterTPC;         // Min number of tpc clusters
  Int_t      fCutMinNClusterITS;         // Min number of its clusters  
  Float_t    fCutMaxChi2PerClusterTPC;   // Max tpc fit chi2 per tpc cluster
  Float_t    fCutMaxChi2PerClusterITS;   // Max its fit chi2 per its cluster
  Bool_t     fCutRequireTPCRefit;        // Require TPC refit
  Bool_t     fCutRequireITSRefit;        // Require ITS refit
  Bool_t     fCutAcceptKinkDaughters;    // Accepting kink daughters?
  Float_t    fCutMaxDCAToVertexXY;       // Track-to-vertex cut in max absolute distance in xy-plane
  Float_t    fCutMaxDCAToVertexZ;        // Track-to-vertex cut in max absolute distance in z-plane
  Bool_t     fCutDCAToVertex2D;          // If true a 2D DCA cut is made.
  Bool_t     fCutRequireITSStandAlone;   // Require ITSStandAlone
  Bool_t     fCutRequireITSpureSA; 	 // ITS pure standalone tracks
  
  
  ClassDef(AliEMCALRecoUtils, 20)
  
};

#endif // ALIEMCALRECOUTILS_H


