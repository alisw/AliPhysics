#ifndef ALIEMCALCLUSTERIZERV3_H
#define ALIEMCALCLUSTERIZERV3_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliEMCALClusterizerv1.h"
class AliEMCALRecPoint; 
class AliEMCALDigit;

// Define numbers rows/columns for topological representation of cells
const Int_t kNrows    = (24+1)*(6+4);  // 10x supermodule rows (6 for EMCAL, 4 for DCAL). +1 accounts for topological gap between two supermodules
const Int_t kNcolumns = 48*2+1;        // 2x  supermodule columns + 1 empty space in between for DCAL (not used for EMCAL)

struct cellWithE {
  cellWithE() : energy(0.), row(0), column(0) {}
  cellWithE(Float_t e, Int_t r, Int_t c) : energy(e), row(r), column(c) {}
  // std::sort will require operator< to compile.
  bool operator<( cellWithE const& rhs ) const
     { return energy < rhs.energy; }
  Float_t energy;
  Int_t row;
  Int_t column;
};


//_________________________________________________________________________
/// \class AliEMCALClusterFinder
/// \ingroup EMCALrec
/// \brief Meta class for recursive clusterizer
///
///  Implementation of same algorithm version as in AliEMCALClusterizerv2,
///  but optimized.
///
/// \author Rudiger Haake (Yale)
//_________________________________________________________________________
class AliEMCALClusterFinder: public TObject
{
  public:
    AliEMCALClusterFinder(TObjArray* outputArray, AliEMCALGeometry* geometry, Double_t timeCut, Double_t timeMin, Double_t timeMax, Double_t gradientCut, Bool_t doEnergyGradientCut, Double_t thresholdSeedE, Double_t thresholdCellE);
    ~AliEMCALClusterFinder();

    Int_t               FindClusters(TClonesArray* digits);
    TObjArray*          GetFoundClusters() {return fFoundClusters;}
    
  private:
    AliEMCALRecPoint*   GetClusterFromNeighbours(AliEMCALRecPoint* recPoint, Int_t row, Int_t column);
    void                GetTopologicalRowColumn(AliEMCALDigit* digit, Int_t& row, Int_t& column);

    AliEMCALGeometry*   fEMCALGeometry;
    cellWithE           fSeedList[kNrows*kNcolumns];
    AliEMCALDigit*      fDigitMap[kNrows][kNcolumns];
    Bool_t              fCellMask[kNrows][kNcolumns];

    TObjArray*          fFoundClusters;
    Int_t               fNumFoundClusters;

    Double_t            fTimeCut;
    Double_t            fTimeMin;
    Double_t            fTimeMax;
    Double_t            fGradientCut;
    Bool_t              fDoEnergyGradientCut;
    Double_t            fThresholdSeedEnergy;
    Double_t            fThresholdCellEnergy;

    AliEMCALClusterFinder(const AliEMCALClusterFinder&);            // not implemented
    AliEMCALClusterFinder &operator=(const AliEMCALClusterFinder&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEMCALClusterFinder,1);
  /// \endcond
};


//_________________________________________________________________________
/// \class AliEMCALClusterizerv3
/// \ingroup EMCALrec
/// \brief Clusterize neighbour cells, split if several maxima
///
///  Implementation version 3 of the clusterization algorithm                     
///  (optimized version of v2)
///  Performs clusterization (collects neighbouring active cells). 
///  It does not allow having more than one local maxima. 
///  Results are stored in TreeR.
///
/// \author Rudiger Haake (Yale)
//_________________________________________________________________________
class AliEMCALClusterizerv3 : public AliEMCALClusterizerv1
{
  
public:
  AliEMCALClusterizerv3() ;         
  AliEMCALClusterizerv3(AliEMCALGeometry* geometry);
  AliEMCALClusterizerv3(AliEMCALGeometry* geometry, AliEMCALCalibData* calib, 
                        AliEMCALCalibTime * calibt, AliCaloCalibPedestal *pedestal);

  virtual ~AliEMCALClusterizerv3();

 
  virtual const char *Version() const { return "clu-v3";}  

  void                SetDoEnGradCut(Bool_t b) { fDoEnGradCut = b; }

protected:
  virtual void        MakeClusters();            

  Bool_t                  fDoEnGradCut;             ///<  cut on energy gradient
  AliEMCALClusterFinder*  fClusterFinder;           //!<! Cluster finder

private:

  AliEMCALClusterizerv3              (const AliEMCALClusterizerv3 &); //copy ctor
  AliEMCALClusterizerv3 & operator = (const AliEMCALClusterizerv3 &);

  /// \cond CLASSIMP
  ClassDef(AliEMCALClusterizerv3,1);
  /// \endcond

};

#endif // AliEMCALCLUSTERIZERV3_H
