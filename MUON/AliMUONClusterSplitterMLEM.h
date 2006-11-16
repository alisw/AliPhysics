#ifndef AliMUONClusterSplitterMLEM_H
#define AliMUONClusterSplitterMLEM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup
/// \class AliMUONClusterSplitterMLEM
/// \brief
/// 
//  Author Alexander Zinchenko, JINR Dubna; Laurent Aphecetche, SUBATECH
//

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#include "TMatrixDfwd.h"

class AliMUONCluster;
class TH2;
class TObjArray;
class AliMUONPad;
class AliMUONMathieson;

class AliMUONClusterSplitterMLEM : public TObject
{
public:
  AliMUONClusterSplitterMLEM(Int_t detElemId, TObjArray* pixArray);
  
  virtual ~AliMUONClusterSplitterMLEM();

  void AddBin(TH2 *mlem, 
              Int_t ic, Int_t jc, Int_t mode, 
              Bool_t *used, TObjArray *pix);
  
  void AddCluster(Int_t ic, Int_t nclust, 
                  TMatrixD& aijcluclu, 
                  Bool_t *used, Int_t *clustNumb, Int_t &nCoupled);
  
  TObject* BinToPix(TH2 *mlem, Int_t jc, Int_t ic);
  
  Float_t ChargeIntegration(Double_t x, Double_t y, const AliMUONPad& pad);
  
  void Fcn1(const AliMUONCluster& cluster, 
            Int_t & npar, Double_t * gin, 
            Double_t &f, Double_t *par, Int_t iflag);
  
  Int_t Fit(const AliMUONCluster& cluster,
            Int_t iSimple, Int_t nfit,
            Int_t *clustFit, TObjArray **clusters, 
            Double_t *parOk, TObjArray& clusterList);
    
  void Merge(const AliMUONCluster& cluster,
             Int_t nForFit, Int_t nCoupled, 
             Int_t *clustNumb, Int_t *clustFit, 
             TObjArray **clusters, 
             TMatrixD& aijcluclu, TMatrixD& aijclupad);
    
  Double_t MinGroupCoupl(Int_t nCoupled, Int_t *clustNumb, 
                         TMatrixD& aijcluclu, Int_t *minGroup);
      
  Int_t SelectPad(const AliMUONCluster& cluster,
                  Int_t nCoupled, Int_t nForFit, 
                  Int_t *clustNumb, Int_t *clustFit, 
                  TMatrixD& aijclupad);
  
  void Split(const AliMUONCluster& cluster,
               TH2* mlem,
               Double_t* coef, TObjArray& clusterList);
  
  
  void UpdatePads(const AliMUONCluster& cluster, Int_t nfit, Double_t *par);

private:
  // will not be implemented
  AliMUONClusterSplitterMLEM(const AliMUONClusterSplitterMLEM&);
  // will not be implemented
  AliMUONClusterSplitterMLEM& operator=(const AliMUONClusterSplitterMLEM&);
  
private:
  
    static const Double_t fgkCouplMin; ///< threshold on coupling 

  TObjArray* fPixArray; //! < not owner
  AliMUONMathieson* fMathieson; //!< 
  Int_t fDetElemId; //!< detection element we are working on
  Int_t fNpar; //!< number of fit parameters
  Double_t fQtot; //!< total charge
  Int_t fnCoupled; //!< number of coupled pixels ?
  
  ClassDef(AliMUONClusterSplitterMLEM,1) // 
};

#endif
