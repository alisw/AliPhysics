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
#include "TH2F.h";

//AliRoot includes
class AliVCluster;
class AliVCaloCells;
#include "AliLog.h"
class AliEMCALGeometry;

class AliEMCALRecoUtils : public TNamed {
  
public:
  
  AliEMCALRecoUtils();
  AliEMCALRecoUtils(const AliEMCALRecoUtils&); 
  AliEMCALRecoUtils& operator=(const AliEMCALRecoUtils&); 
  virtual ~AliEMCALRecoUtils() ;
  
  enum NonlinearityFunctions{kPi0MC=0,kPi0GammaGamma=1,kPi0GammaConversion=2,kNoCorrection=3};
  enum PositionAlgorithms{kPosTowerIndex=0, kPosTowerGlobal};
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

private:
  
  Float_t fMisalTransShift[15];   // Shift parameters
  Float_t fMisalRotShift[15];     // Shift parameters
  Int_t   fNonLinearityFunction;  // Non linearity function choice
  Float_t fNonLinearityParams[6]; // Parameters for the non linearity function
  Int_t   fParticleType;          // Particle type for depth calculation
  Int_t   fPosAlgo;               // Position recalculation algorithm
  Float_t fW0;                    // Weight0
  Bool_t       fRecalibration;             //  Switch on or off the recalibration
  TObjArray  * fEMCALRecalibrationFactors; // Array of histograms with map of recalibration factors, EMCAL

  ClassDef(AliEMCALRecoUtils, 3)
  
};

#endif // ALIEMCALRECOUTILS_H


