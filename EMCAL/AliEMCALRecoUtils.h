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
#include "TNamed.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TH2F.h"

//AliRoot includes
class AliVCluster;
class AliVCaloCells;
class AliVEvent;
#include "AliLog.h"

// EMCAL includes
class AliEMCALGeometry;
class AliEMCALPIDUtils;
class AliESDtrack;

class AliEMCALRecoUtils : public TNamed {
  
public:
  
  AliEMCALRecoUtils();
  AliEMCALRecoUtils(const AliEMCALRecoUtils&); 
  AliEMCALRecoUtils& operator=(const AliEMCALRecoUtils&); 
  virtual ~AliEMCALRecoUtils() ;  
  void Print(const Option_t*) const;

  //enums
  enum NonlinearityFunctions{kPi0MC=0,kPi0GammaGamma=1,kPi0GammaConversion=2,kNoCorrection=3,kBeamTest=4,kBeamTestCorrected=5};
  enum PositionAlgorithms{kUnchanged=-1,kPosTowerIndex=0, kPosTowerGlobal=1};
  enum ParticleType{kPhoton=0, kElectron=1,kHadron =2, kUnknown=-1};
  enum { kNCuts = 11 }; //track matching

  //-----------------------------------------------------
  //Position recalculation
  //-----------------------------------------------------

  void     RecalculateClusterPosition(AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu); 
  void     RecalculateClusterPositionFromTowerIndex (AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu); 
  void     RecalculateClusterPositionFromTowerGlobal(AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu); 
  
  Float_t  GetCellWeight(const Float_t eCell, const Float_t eCluster) const { return TMath::Max( 0., fW0 + TMath::Log( eCell / eCluster ));}
  
  Float_t  GetDepth(const Float_t eCluster, const Int_t iParticle, const Int_t iSM) const ; 
  
  void     GetMaxEnergyCell(AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu, 
                            Int_t & absId,  Int_t& iSupMod, Int_t& ieta, Int_t& iphi, Bool_t &shared);
  
  Float_t  GetMisalTransShift(const Int_t i) const {
    if(i < 15 ){return fMisalTransShift[i]; }
    else { AliInfo(Form("Index %d larger than 15, do nothing\n",i)); return 0.;}
  }
  Float_t* GetMisalTransShiftArray() {return fMisalTransShift; }

  void     SetMisalTransShift(const Int_t i, const Float_t shift) {
    if(i < 15 ){fMisalTransShift[i] = shift; }
    else { AliInfo(Form("Index %d larger than 15, do nothing\n",i));}
  }
  void     SetMisalTransShiftArray(Float_t * misal) 
  { for(Int_t i = 0; i < 15; i++)fMisalTransShift[i] = misal[i]; }

  Float_t  GetMisalRotShift(const Int_t i) const {
    if(i < 15 ){return fMisalRotShift[i]; }
    else { AliInfo(Form("Index %d larger than 15, do nothing\n",i)); return 0.;}
  }
  Float_t* GetMisalRotShiftArray() {return fMisalRotShift; }
  
  void     SetMisalRotShift(const Int_t i, const Float_t shift) {
    if(i < 15 ){fMisalRotShift[i] = shift; }
    else { AliInfo(Form("Index %d larger than 15, do nothing\n",i));}
  }
  void     SetMisalRotShiftArray(Float_t * misal) 
  { for(Int_t i = 0; i < 15; i++)fMisalRotShift[i] = misal[i]; }
  
  Int_t    GetParticleType()        const  { return  fParticleType    ;}
  void     SetParticleType(Int_t particle) { fParticleType = particle ;}
  
  Int_t    GetPositionAlgorithm()   const  { return fPosAlgo          ;}
  void     SetPositionAlgorithm(Int_t alg) { fPosAlgo = alg           ;}
  
  Float_t  GetW0()                 const   { return fW0               ;}
  void     SetW0(Float_t w0)               { fW0  = w0                ;}

  //-----------------------------------------------------
  //Non Linearity
  //-----------------------------------------------------

  Float_t CorrectClusterEnergyLinearity(AliVCluster* clu);
  
  Float_t  GetNonLinearityParam(const Int_t i) const {
    if(i < 6 ){return fNonLinearityParams[i]; }
    else { AliInfo(Form("Index %d larger than 6, do nothing\n",i)); return 0.;}
  }
  void     SetNonLinearityParam(const Int_t i, const Float_t param) {
    if(i < 6 ){fNonLinearityParams[i] = param; }
    else { AliInfo(Form("Index %d larger than 6, do nothing\n",i));}
  }
  
  Int_t GetNonLinearityFunction() const    { return fNonLinearityFunction ;}
  void  SetNonLinearityFunction(Int_t fun) { fNonLinearityFunction = fun  ;}
    
  //-----------------------------------------------------
  //Recalibration
  //-----------------------------------------------------

  void RecalibrateClusterEnergy(AliEMCALGeometry* geom, AliVCluster* cluster, AliVCaloCells * cells);

  Bool_t IsRecalibrationOn()           const { return fRecalibration ; }
  void   SwitchOnRecalibration()             { fRecalibration = kTRUE ; if(!fEMCALRecalibrationFactors)InitEMCALRecalibrationFactors();}
  void   SwitchOffRecalibration()            { fRecalibration = kFALSE ; }
  void   InitEMCALRecalibrationFactors() ;

  //Recalibrate channels with time dependent corrections
  void SwitchOnTimeDepCorrection()          { fUseTimeCorrectionFactors = kTRUE ; SwitchOnRecalibration();}
  void SwitchOffTimeDepCorrection()         { fUseTimeCorrectionFactors = kFALSE;}
  void SetTimeDependentCorrections(Int_t runnumber);
    
  Float_t GetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow) const { 
    if(fEMCALRecalibrationFactors) return (Float_t) ((TH2F*)fEMCALRecalibrationFactors->At(iSM))->GetBinContent(iCol,iRow); 
    else return 1;}
	
  void SetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
    if(!fEMCALRecalibrationFactors) InitEMCALRecalibrationFactors();
    ((TH2F*)fEMCALRecalibrationFactors->At(iSM))->SetBinContent(iCol,iRow,c);}  
  
  TH2F * GetEMCALChannelRecalibrationFactors(Int_t iSM)   const { return (TH2F*)fEMCALRecalibrationFactors->At(iSM) ;}	
  void SetEMCALChannelRecalibrationFactors(TObjArray *map)      { fEMCALRecalibrationFactors = map                  ;}
  void SetEMCALChannelRecalibrationFactors(Int_t iSM , TH2F* h) { fEMCALRecalibrationFactors->AddAt(h,iSM)          ;}

  //-----------------------------------------------------
  //Modules fiducial region, remove clusters in borders
  //-----------------------------------------------------

  Bool_t CheckCellFiducialRegion(AliEMCALGeometry* geom, AliVCluster* cluster, AliVCaloCells* cells) ;
  void   SetNumberOfCellsFromEMCALBorder(Int_t n) { fNCellsFromEMCALBorder = n    ;}
  Int_t  GetNumberOfCellsFromEMCALBorder() const  { return fNCellsFromEMCALBorder ;}
    
  void   SwitchOnNoFiducialBorderInEMCALEta0()    { fNoEMCALBorderAtEta0 = kTRUE  ;}
  void   SwitchOffNoFiducialBorderInEMCALEta0()   { fNoEMCALBorderAtEta0 = kFALSE ;}
  Bool_t IsEMCALNoBorderAtEta0()                  { return fNoEMCALBorderAtEta0   ;}
  
  //-----------------------------------------------------
  // Bad channels
  //-----------------------------------------------------

  Bool_t IsBadChannelsRemovalSwitchedOn()     const { return fRemoveBadChannels       ;}
  void SwitchOnBadChannelsRemoval ()                { fRemoveBadChannels = kTRUE ; if(!fEMCALBadChannelMap)InitEMCALBadChannelStatusMap();}
  void SwitchOffBadChannelsRemoval()                { fRemoveBadChannels = kFALSE     ;}
	
  Bool_t IsDistanceToBadChannelRecalculated() const { return fRecalDistToBadChannels  ;}
  void SwitchOnDistToBadChannelRecalculation()      { fRecalDistToBadChannels = kTRUE  ; if(!fEMCALBadChannelMap)InitEMCALBadChannelStatusMap();}
  void SwitchOffDistToBadChannelRecalculation()     { fRecalDistToBadChannels = kFALSE ;}
  
  void InitEMCALBadChannelStatusMap() ;
	
  Int_t GetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow) const { 
    if(fEMCALBadChannelMap) return (Int_t) ((TH2I*)fEMCALBadChannelMap->At(iSM))->GetBinContent(iCol,iRow); 
    else return 0;}//Channel is ok by default
	
  void SetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
    if(!fEMCALBadChannelMap)InitEMCALBadChannelStatusMap() ;
    ((TH2I*)fEMCALBadChannelMap->At(iSM))->SetBinContent(iCol,iRow,c);}
	
  TH2I * GetEMCALChannelStatusMap(Int_t iSM) const {return (TH2I*)fEMCALBadChannelMap->At(iSM);}
  void   SetEMCALChannelStatusMap(TObjArray *map)  {fEMCALBadChannelMap = map;}
  void   SetEMCALChannelStatusMap(Int_t iSM , TH2I* h) {fEMCALBadChannelMap->AddAt(h,iSM);}

  Bool_t ClusterContainsBadChannel(AliEMCALGeometry* geom, UShort_t* cellList, Int_t nCells);
 
  //-----------------------------------------------------
  // Recalculate other cluster parameters
  //-----------------------------------------------------

  void RecalculateClusterDistanceToBadChannel(AliEMCALGeometry * geom, AliVCaloCells* cells, AliVCluster * cluster);
  void RecalculateClusterPID(AliVCluster * cluster);

  AliEMCALPIDUtils * GetPIDUtils() { return fPIDUtils;}

  void RecalculateClusterShowerShapeParameters(AliEMCALGeometry * geom, AliVCaloCells* cells, AliVCluster * cluster);

  //----------------------------------------------------
  // Track matching
  //----------------------------------------------------

  void    FindMatches(AliVEvent *event, TObjArray * clusterArr=0x0);
  void    GetMatchedResiduals(Int_t index, Float_t &dR, Float_t &dZ);
  Int_t   GetMatchedTrackIndex(Int_t index);
  Bool_t  IsMatched(Int_t index);
  UInt_t  FindMatchedPos(Int_t index) const;

  Float_t GetCutR()       const { return fCutR ;}
  Float_t GetCutZ()       const { return fCutZ ;}
  void    SetCutR(Float_t cutR) { fCutR=cutR   ;}
  void    SetCutZ(Float_t cutZ) { fCutZ=cutZ   ;}

  //Track Cuts 
  Bool_t  IsAccepted(AliESDtrack *track);
  void    InitTrackCuts();

  // track quality cut setters  
  void    SetMinNClustersTPC(Int_t min=-1)          { fCutMinNClusterTPC       = min  ;}
  void    SetMinNClustersITS(Int_t min=-1)          { fCutMinNClusterITS       = min  ;}
  void    SetMaxChi2PerClusterTPC(Float_t max=1e10) { fCutMaxChi2PerClusterTPC = max  ;}
  void    SetMaxChi2PerClusterITS(Float_t max=1e10) { fCutMaxChi2PerClusterITS = max  ;}
  void    SetRequireTPCRefit(Bool_t b=kFALSE)       { fCutRequireTPCRefit      = b    ;}
  void    SetRequireITSRefit(Bool_t b=kFALSE)       { fCutRequireITSRefit      = b    ;}
  void    SetAcceptKinkDaughters(Bool_t b=kTRUE)    { fCutAcceptKinkDaughters  = b    ;}
  void    SetMaxDCAToVertexXY(Float_t dist=1e10)    { fCutMaxDCAToVertexXY     = dist ;}
  void    SetMaxDCAToVertexZ(Float_t dist=1e10)     { fCutMaxDCAToVertexZ      = dist ;}
  void    SetDCAToVertex2D(Bool_t b=kFALSE)         { fCutDCAToVertex2D        = b    ;}

  // getters
  Int_t   GetMinNClusterTPC()               const   { return fCutMinNClusterTPC       ;}
  Int_t   GetMinNClustersITS()              const   { return fCutMinNClusterITS       ;}
  Float_t GetMaxChi2PerClusterTPC()         const   { return fCutMaxChi2PerClusterTPC ;}
  Float_t GetMaxChi2PerClusterITS()         const   { return fCutMaxChi2PerClusterITS ;}
  Bool_t  GetRequireTPCRefit()              const   { return fCutRequireTPCRefit      ;}
  Bool_t  GetRequireITSRefit()              const   { return fCutRequireITSRefit      ;}
  Bool_t  GetAcceptKinkDaughters()          const   { return fCutAcceptKinkDaughters  ;}
  Float_t GetMaxDCAToVertexXY()             const   { return fCutMaxDCAToVertexXY     ;}
  Float_t GetMaxDCAToVertexZ()              const   { return fCutMaxDCAToVertexZ      ;}
  Bool_t  GetDCAToVertex2D()                const   { return fCutDCAToVertex2D        ;}


private:
  
  //Position recalculation
  Float_t    fMisalTransShift[15];       // Shift parameters
  Float_t    fMisalRotShift[15];         // Shift parameters
  Int_t      fNonLinearityFunction;      // Non linearity function choice
  Float_t    fNonLinearityParams[6];     // Parameters for the non linearity function
  Int_t      fParticleType;              // Particle type for depth calculation
  Int_t      fPosAlgo;                   // Position recalculation algorithm
  Float_t    fW0;                        // Weight0
  
  // Recalibration 
  Bool_t     fRecalibration;             // Switch on or off the recalibration
  TObjArray* fEMCALRecalibrationFactors; // Array of histograms with map of recalibration factors, EMCAL

  // Bad Channels
  Bool_t     fRemoveBadChannels;         // Check the channel status provided and remove clusters with bad channels
  Bool_t     fRecalDistToBadChannels;    // Calculate distance from highest energy tower of cluster to closes bad channel
  TObjArray* fEMCALBadChannelMap;        // Array of histograms with map of bad channels, EMCAL

  // Border cells
  Int_t      fNCellsFromEMCALBorder;     // Number of cells from EMCAL border the cell with maximum amplitude has to be.
  Bool_t     fNoEMCALBorderAtEta0;       // Do fiducial cut in EMCAL region eta = 0?
  
  //Track matching 
  TArrayI  * fMatchedTrackIndex;         // Array that stores indexes of matched tracks      
  TArrayI  * fMatchedClusterIndex;       // Array that stores indexes of matched clusters
  TArrayF  * fResidualZ;                 // Array that stores the residual z
  TArrayF  * fResidualR;                 // Array that stores the residual r
  Float_t    fCutR;                      // dR cut on matching
  Float_t    fCutZ;                      // dZ cut on matching
  
  Int_t      fCutMinNClusterTPC;         // Min number of tpc clusters
  Int_t      fCutMinNClusterITS;         // Min number of its clusters  
  Float_t    fCutMaxChi2PerClusterTPC;   // Max tpc fit chi2 per tpc cluster
  Float_t    fCutMaxChi2PerClusterITS;   // Max its fit chi2 per its cluster
  Bool_t     fCutRequireTPCRefit;        // Require TPC refit
  Bool_t     fCutRequireITSRefit;        // Require ITS refit
  Bool_t     fCutAcceptKinkDaughters;    // Accepting kink daughters?
  Float_t    fCutMaxDCAToVertexXY;       // Track-to-vertex cut in max absolute distance in xy-plane
  Float_t    fCutMaxDCAToVertexZ;        // Track-to-vertex cut in max absolute distance in z-plane
  Bool_t     fCutDCAToVertex2D;          // If true a 2D DCA cut is made. Tracks are accepted if sqrt((DCAXY / fCutMaxDCAToVertexXY)^2 + (DCAZ / fCutMaxDCAToVertexZ)^2) < 1 AND sqrt((DCAXY / fCutMinDCAToVertexXY)^2 + (DCAZ / fCutMinDCAToVertexZ)^2) > 1

  //PID
  AliEMCALPIDUtils * fPIDUtils;          // Recalculate PID parameters
  
  //Time Correction
  Bool_t     fUseTimeCorrectionFactors;  // Use Time Dependent Correction
  Bool_t     fTimeCorrectionFactorsSet;  // Time Correction set at leat once
  
  ClassDef(AliEMCALRecoUtils, 7)
  
};

#endif // ALIEMCALRECOUTILS_H


