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
///////////////////////////////////////////////////////////////////////////////

//Root includes
#include "TNamed.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TH2F.h"

//AliRoot includes
class AliVCluster;
class AliVCaloCells;
#include "AliLog.h"
class AliEMCALGeometry;
class AliEMCALPIDUtils;

class AliEMCALRecoUtils : public TNamed {
  
public:
  
  AliEMCALRecoUtils();
  AliEMCALRecoUtils(const AliEMCALRecoUtils&); 
  AliEMCALRecoUtils& operator=(const AliEMCALRecoUtils&); 
  virtual ~AliEMCALRecoUtils() ;
  
  enum NonlinearityFunctions{kPi0MC=0,kPi0GammaGamma=1,kPi0GammaConversion=2,kNoCorrection=3};
  enum PositionAlgorithms{kUnchanged=-1,kPosTowerIndex=0, kPosTowerGlobal=1};
  enum ParticleType{kPhoton=0, kElectron=1,kHadron =2, kUnknown=-1};
  
  //Position recalculation
  void     RecalculateClusterPosition(AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu); 
  void     RecalculateClusterPositionFromTowerIndex (AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu); 
  void     RecalculateClusterPositionFromTowerGlobal(AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu); 
  
  Float_t  GetCellWeight(const Float_t eCell, const Float_t eCluster) const { return TMath::Max( 0., fW0 + TMath::Log( eCell / eCluster ));}
  
  Float_t  GetDepth(const Float_t eCluster, const Int_t iParticle, const Int_t iSM) const ; 
  
  void     GetMaxEnergyCell(AliEMCALGeometry *geom, AliVCaloCells* cells, AliVCluster* clu, 
                            Int_t & absId,  Int_t& iSupMod, Int_t& ieta, Int_t& iphi);
  
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
  
  Int_t    GetParticleType() const         {return  fParticleType    ;}
  void     SetParticleType(Int_t particle) {fParticleType = particle ;}
  
  Int_t    GetPositionAlgorithm() const    {return fPosAlgo;}
  void     SetPositionAlgorithm(Int_t alg) {fPosAlgo = alg ;}
  
  Float_t  GetW0() const     {return fW0;}
  void     SetW0(Float_t w0) {fW0  = w0 ;}

  //Non Linearity
  
  Float_t CorrectClusterEnergyLinearity(AliVCluster* clu);
  
  Float_t  GetNonLinearityParam(const Int_t i) const {
    if(i < 6 ){return fNonLinearityParams[i]; }
    else { AliInfo(Form("Index %d larger than 6, do nothing\n",i)); return 0.;}
  }
  void     SetNonLinearityParam(const Int_t i, const Float_t param) {
    if(i < 6 ){fNonLinearityParams[i] = param; }
    else { AliInfo(Form("Index %d larger than 6, do nothing\n",i));}
  }
  
  Int_t GetNonLinearityFunction() const    {return fNonLinearityFunction;}
  void  SetNonLinearityFunction(Int_t fun) {fNonLinearityFunction = fun ;}
  
  void Print(const Option_t*) const;
  
  //Recalibration
  void RecalibrateClusterEnergy(AliEMCALGeometry* geom, AliVCluster* cluster, AliVCaloCells * cells);

  Bool_t IsRecalibrationOn()  const { return fRecalibration ; }
  void SwitchOnRecalibration()    {fRecalibration = kTRUE ; InitEMCALRecalibrationFactors();}
  void SwitchOffRecalibration()   {fRecalibration = kFALSE ; }
  
  void InitEMCALRecalibrationFactors() ;
  
  Float_t GetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow) const { 
    if(fEMCALRecalibrationFactors) return (Float_t) ((TH2F*)fEMCALRecalibrationFactors->At(iSM))->GetBinContent(iCol,iRow); 
    else return 1;}
	
  void SetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
    if(!fEMCALRecalibrationFactors) InitEMCALRecalibrationFactors();
    ((TH2F*)fEMCALRecalibrationFactors->At(iSM))->SetBinContent(iCol,iRow,c);}  
  
  TH2F * GetEMCALChannelRecalibrationFactors(Int_t iSM) const {return (TH2F*)fEMCALRecalibrationFactors->At(iSM);}	
  void SetEMCALChannelRecalibrationFactors(TObjArray *map)      {fEMCALRecalibrationFactors = map;}
  void SetEMCALChannelRecalibrationFactors(Int_t iSM , TH2F* h) {fEMCALRecalibrationFactors->AddAt(h,iSM);}

  //Modules fiducial region, remove clusters in borders
  Bool_t CheckCellFiducialRegion(AliEMCALGeometry* geom, AliVCluster* cluster, AliVCaloCells* cells) ;
  void   SetNumberOfCellsFromEMCALBorder(Int_t n) {fNCellsFromEMCALBorder = n; }
  Int_t  GetNumberOfCellsFromEMCALBorder() const  {return fNCellsFromEMCALBorder; }
    
  void   SwitchOnNoFiducialBorderInEMCALEta0()  {fNoEMCALBorderAtEta0 = kTRUE; }
  void   SwitchOffNoFiducialBorderInEMCALEta0() {fNoEMCALBorderAtEta0 = kFALSE; }
  Bool_t IsEMCALNoBorderAtEta0()                {return fNoEMCALBorderAtEta0;}
  
  // Bad channels
  Bool_t IsBadChannelsRemovalSwitchedOn()  const { return fRemoveBadChannels ; }
  void SwitchOnBadChannelsRemoval ()  {fRemoveBadChannels = kTRUE  ; InitEMCALBadChannelStatusMap();}
  void SwitchOffBadChannelsRemoval()  {fRemoveBadChannels = kFALSE ; }
	
  void InitEMCALBadChannelStatusMap() ;
	
  Int_t GetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow) const { 
    if(fEMCALBadChannelMap) return (Int_t) ((TH2I*)fEMCALBadChannelMap->At(iSM))->GetBinContent(iCol,iRow); 
    else return 0;}//Channel is ok by default
	
  void SetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow, Double_t c = 1) { 
    if(!fEMCALBadChannelMap)InitEMCALBadChannelStatusMap() ;
    ((TH2I*)fEMCALBadChannelMap->At(iSM))->SetBinContent(iCol,iRow,c);}
	
  TH2I * GetEMCALChannelStatusMap(Int_t iSM) const {return (TH2I*)fEMCALBadChannelMap->At(iSM);}
  void   SetEMCALChannelStatusMap(TObjArray *map)  {fEMCALBadChannelMap = map;}
	
  Bool_t ClusterContainsBadChannel(AliEMCALGeometry* geom, UShort_t* cellList, Int_t nCells);
 
  //Recalculate other cluster parameters
  void RecalculateClusterPID(AliVCluster * cluster);
  AliEMCALPIDUtils * GetPIDUtils() { return fPIDUtils;}

  void RecalculateClusterShowerShapeParameters(AliEMCALGeometry * geom, AliVCaloCells* cells, AliVCluster * cluster);


private:
  
  Float_t fMisalTransShift[15];   // Shift parameters
  Float_t fMisalRotShift[15];     // Shift parameters
  Int_t   fNonLinearityFunction;  // Non linearity function choice
  Float_t fNonLinearityParams[6]; // Parameters for the non linearity function
  Int_t   fParticleType;          // Particle type for depth calculation
  Int_t   fPosAlgo;               // Position recalculation algorithm
  Float_t fW0;                    // Weight0
  
  Bool_t     fRecalibration;             // Switch on or off the recalibration
  TObjArray* fEMCALRecalibrationFactors; // Array of histograms with map of recalibration factors, EMCAL
  Bool_t     fRemoveBadChannels;         // Check the channel status provided and remove clusters with bad channels
  TObjArray* fEMCALBadChannelMap;        // Array of histograms with map of bad channels, EMCAL
  Int_t      fNCellsFromEMCALBorder;     // Number of cells from EMCAL border the cell with maximum amplitude has to be.
  Bool_t     fNoEMCALBorderAtEta0;       // Do fiducial cut in EMCAL region eta = 0?

  AliEMCALPIDUtils * fPIDUtils;               // Recalculate PID parameters
  
  ClassDef(AliEMCALRecoUtils, 5)
  
};

#endif // ALIEMCALRECOUTILS_H


