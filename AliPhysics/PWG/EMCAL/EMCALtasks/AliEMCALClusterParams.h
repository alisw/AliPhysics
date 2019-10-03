#ifndef ALIEMCALCLUSTERPARAMS_H
#define ALIEMCALCLUSTERPARAMS_H

// $Id$

#include <TObject.h>
class AliESDtrack;
class AliESDCaloCluster;
class AliEMCALGeometry;
class AliESDCaloCells;
class TObjArray;

class AliEMCALClusterParams : public TObject 
{
 public: 
  AliEMCALClusterParams(AliESDtrack       *trackin, 
                        AliESDCaloCluster *clusin, 
                        AliEMCALGeometry  *geometryin, 
                        AliESDCaloCells   *cellsin);
  virtual ~AliEMCALClusterParams() {;}

  Double_t   GetPe()                  const; 
  Int_t      IsElectron()             const;  
  void       LoopThroughCells()       const;
  void       PrintClusterParameters() const;

  //==========Un-log-weighted parameters============================================================
  void       GetCentroid(Double_t &xback, Double_t &yback, Double_t &rback) const;
  Double_t   GetR(Double_t x, Double_t y)                const;                //Uses Centroid
  Double_t   GetRfactor()                                const;                //Uses GetR,Centroid
  Double_t   ElectronFraction(Double_t r, Double_t tce)  const;                //Only used locally for K-factor
  Double_t   GetKfactor()                                const;                //Gets K-factor, uses electronfraction
  Double_t   GetDispersionX()                            const;                //Gets unweighted dispersionX
  Double_t   GetDispersionY()                            const;                //Gets unweighted dispersionY
  Double_t   GetDispersionMax()                          const;                //Gets unweighted dispersionMax
  void       GetEllipseParameters(Double_t &param1, Double_t &param2) const;   //Gets M02 M20
  Double_t   GetDispersion()                             const;                //Unweighted Dispersion

  //==========Log-weighted parameters===============================================================
  void       GetWeightedCentroid(Double_t &xback, Double_t &yback, Double_t &rback) const;
  Double_t   GetWeightedR(Double_t x, Double_t y)                             const;       //Weighted r, uses WeightedCentriod
  Double_t   GetWeightedRfactor()                                             const;       //Uses WeightedR
  Double_t   ElectronfractionWeighted(Double_t r, Double_t tce)               const;       //Used locally for K-factor weighted
  Double_t   GetWeightedKfactor()                                             const;       //Gets Weighted K-factor
  Double_t   GetWeightedDispersionX()                                         const;       //Gets Weighted dispersionX
  Double_t   GetWeightedDispersionY()                                         const;       //Gets Weighted dispersionY
  Double_t   GetWeightedDispersionMax()                                       const;       //Gets Weighted dispersionMax
  void       GetWeightedEllipseParameters(Double_t &param1, Double_t &param2) const;       //Gets weighted M02 M20
  Double_t   GetWeightedDispersion(Double_t &dispersionback)                  const;       //Gets Weighted dispersion

  void     RecalculateClusterShowerShapeParameters(Double_t &m02, Double_t &m20, Double_t &dispersion) const; //from Gustavo's code!

 private:
  AliESDtrack        *fTrack;       //!
  AliESDCaloCluster  *fCluster;     //!
  AliEMCALGeometry   *fGeom;        //!
  AliESDCaloCells    *fCells;       //!
	
  AliEMCALClusterParams(const AliEMCALClusterParams & g);
  AliEMCALClusterParams & operator = (const AliEMCALClusterParams & g);

  ClassDef(AliEMCALClusterParams,0) 
};
#endif //ALIEMCALCLUSTERPARAMS_H
