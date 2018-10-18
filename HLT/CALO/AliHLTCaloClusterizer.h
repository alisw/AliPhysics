//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Rudiger Haake                                         *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTCALOCLUSTERIZER_H
#define ALIHLTCALOCLUSTERIZER_H


/**
 * Class does clusterization for EMCAL. Uses same algorithm
 * as implemented in AliEMCALClusterizerv3
 *
 * @file   AliHLTCaloClusterizer.h
 * @author Rudiger Haake (Yale)
 * @date
 * @brief  Clusterizer for CALO HLT
 */

#include "AliHLTCaloRecPointContainerStruct.h"
#include "AliHLTCaloRecPointDataStruct.h"
#include "AliHLTCaloDigitContainerDataStruct.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "TString.h"
#include "AliHLTCaloConstantsHandler.h"

#include "AliHLTLogging.h"

class TClonesArray;
class TString;
class TH2;
class AliHLTEMCALGeometry;


namespace HLTClusterFinder
{
  // Define numbers rows/columns for topological representation of cells
  const UInt_t kNrows    = (24+1)*(6+4);  // 10x supermodule rows (6 for EMCAL, 4 for DCAL). +1 accounts for topological gap between two supermodules
  const UInt_t kNcolumns = 48*2+1;        // 2x  supermodule columns + 1 empty space in between for DCAL (not used for EMCAL)
}

//_________________________________________________________________________
/// \class AliHLTClusterFinder
/// \brief Meta class for recursive clusterizer
///
///  Implementation of same algorithm version as in AliEMCALClusterizerv3
///  optimized for HLT workflow
///
/// \author Rudiger Haake (Yale)
//_________________________________________________________________________
class AliHLTClusterFinder : public AliHLTLogging
{
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

  public:
    AliHLTClusterFinder(AliHLTCaloRecPointDataStruct** outputArray, Int_t maxNumClusters, AliHLTEMCALGeometry* geometry, Double_t timeCut, Double_t timeMin, Double_t timeMax, Double_t gradientCut, Bool_t doEnergyGradientCut, Double_t thresholdSeedE, Double_t thresholdCellE);
    ~AliHLTClusterFinder();

    Int_t               FindClusters(AliHLTCaloDigitDataStruct** digitArray, Int_t numDigits);
    AliHLTCaloRecPointDataStruct**  GetFoundClusters() {return fFoundClusters;}
    void                AddBadCellsForSM(Int_t sm, TH2* hBadCells);
private:
    AliHLTCaloRecPointDataStruct*  GetClusterFromNeighbours(AliHLTCaloRecPointDataStruct* recPoint, Int_t row, Int_t column);
    void                GetTopologicalRowColumn(AliHLTCaloDigitDataStruct* digit, Int_t& row, Int_t& column);
    void                GetTopologicalRowColumn(UInt_t digitID, Int_t& row, Int_t& column);
    void                CalculateCenterOfGravity(AliHLTCaloRecPointDataStruct* recPoint);

    cellWithE           fSeedList[HLTClusterFinder::kNrows*HLTClusterFinder::kNcolumns];      //!<! seed array
    AliHLTCaloDigitDataStruct*      fDigitMap[HLTClusterFinder::kNrows][HLTClusterFinder::kNcolumns];     //!<! topology arrays
    Bool_t              fCellMask[HLTClusterFinder::kNrows][HLTClusterFinder::kNcolumns];     //!<! topology arrays
    Bool_t              fBadCellMask[HLTClusterFinder::kNrows][HLTClusterFinder::kNcolumns];     //!<! topology arrays for bad cells

    AliHLTCaloDigitDataStruct*      fCurrentClusterDigits[500]; //!<! temporary array of digits belonging to a cluster (is used to on-cluster calculations,e.g. center of gravity)
    AliHLTCaloRecPointDataStruct**   fFoundClusters;      //!<! Pointer to found cluster object array
    Int_t               fNumFoundClusters;                ///<  number of found clusters in FindClusters()
    Int_t               fMaxNumClusters;                  ///<  available buffer

    AliHLTEMCALGeometry* fGeometry;                        ///<  HLT geometry object
    Double_t            fTimeCut;                         ///<  maximum time difference between the digits inside EMC cluster
    Double_t            fTimeMin;                         ///<  minimum time of physical signal in a cell/digit
    Double_t            fTimeMax;                         ///<  maximum time of physical signal in a cell/digit 
    Double_t            fGradientCut;                     ///<  minimum energy difference to distinguish local maxima in a cluster
    Bool_t              fDoEnergyGradientCut;             ///<  cut on energy gradient
    Double_t            fThresholdSeedEnergy;             ///<  minimum energy to seed a EC digit in a cluster
    Double_t            fThresholdCellEnergy;             ///<  minimum energy for a digit to be a member of a cluster

    AliHLTClusterFinder(const AliHLTClusterFinder&);            // not implemented
    AliHLTClusterFinder &operator=(const AliHLTClusterFinder&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliHLTClusterFinder,1);
  /// \endcond
};


/** 
 * @class AliHLTCaloClusterizer
 * Wrapper class that uses AliHLTClusterFinder to form digits to cluster
 * for the EMCAL/DCAL
 *
 * @ingroup alihlt_calo
 */


class AliHLTCaloClusterizer : public AliHLTCaloConstantsHandler, public AliHLTLogging
{
  
public:
  
  /** Constructor */
  AliHLTCaloClusterizer(TString det);    

  /** Destructor */
  virtual ~AliHLTCaloClusterizer();
  
  /** Set geometry object */
  void SetGeometry(AliHLTEMCALGeometry* geometry) { fGeometry = geometry; }

  /** Set array with digits */
  void SetDigitArray(AliHLTCaloDigitDataStruct **digitPointerArr)
  { fDigitsPointerArray = digitPointerArr; } 

  /** Set emc clustering threshold */
  void SetEmcClusteringThreshold(Float_t threshold) { fEmcClusteringThreshold = threshold; }

  /** Set emc min energy threshold */
  void SetEmcMinEnergyThreshold(Float_t threshold) { fEmcMinEnergyThreshold = threshold; }

  /** Set emc time gate */
  void SetEmcTimeGate(Float_t gate) { fEmcTimeGate = gate; }

  /** Set cell time min */
  void SetCellTimeMin(Float_t min) { fCellTimeMin = min; }

  /** Set cell time max */
  void SetCellTimeMax(Float_t max) { fCellTimeMax = max; }

  /** Activate usage of gradient cut*/
  void SetUseGradientCut(Bool_t val) {fUseGradientCut = val; }

  /** Set gradient cut value*/
  void SetGradientCut(Float_t val) {fGradientCut = val; }



  
  /** Starts clusterization of the event */ 
  virtual Int_t ClusterizeEvent(Int_t nDigits);
  

  /**
  * Get pointer to the rec points array
  */
  AliHLTCaloRecPointDataStruct** GetRecPoints() const { return fRecPointArray; }

protected:

  /** Array of pointers to the rec point output */
  AliHLTCaloRecPointDataStruct **fRecPointArray; //COMMENT

  /** Size of the rec point array */
  Int_t fArraySize;
  
  /** Energy threshold for starting a cluster for the calorimeter */
  Float_t fEmcClusteringThreshold;                             //COMMENT

  /** Energy threshold for including a crystal in a cluster */
  Float_t fEmcMinEnergyThreshold;                              //COMMENT

  /** Maximum time difference for inclusion in a rec point */
  Float_t fEmcTimeGate;                                        //COMMENT

  /** Min time for cell */
  Float_t fCellTimeMin;                                        //COMMENT

  /** Max time for cell */
  Float_t fCellTimeMax;                                        //COMMENT

  /** Usage of gradient cut */
  Bool_t  fUseGradientCut;                                     //COMMENT

  /** Threshold for gradient cut */
  Float_t fGradientCut;                                        //COMMENT

  

  /** Array of our digits */
  AliHLTCaloDigitDataStruct **fDigitsPointerArray;             //! transient

    /** Array of our digits */
  AliHLTEMCALGeometry*        fGeometry;                       //! transient

  AliHLTClusterFinder*        fClusterFinder;                  //! transient
  
private:

  /** Default constructor, prohibited */
  AliHLTCaloClusterizer();                          // COMMENT
  
  /** Copy constructor, prohibited */
  AliHLTCaloClusterizer (const AliHLTCaloClusterizer &); //COMMENT
  
  /** Assignment operator, prohibited */
  AliHLTCaloClusterizer & operator = (const AliHLTCaloClusterizer &); //COMMENT

  ClassDef(AliHLTCaloClusterizer, 1);

};

#endif
