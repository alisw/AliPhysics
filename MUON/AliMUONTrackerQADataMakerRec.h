#ifndef ALIMUONTRACKERQADATAMAKERREC_H
#define ALIMUONTRACKERQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id: AliMUONTrackerQADataMakerRec.h 35760 2009-10-21 21:45:42Z ivana $

/// \ingroup rec
/// \class AliMUONTrackerQADataMakerRec
/// \brief MUON Quality assurance data maker
///

// --- AliRoot header files ---
#include "AliMUONVQADataMakerRec.h"
#include "AliMUONRecoParam.h"

class AliMUONDigitMaker;
class AliMUONVClusterStore;
class AliMUONVDigitStore;
class AliMUONVStore;
class AliMUONVTrackerData;
class AliMUONVTrackerDataMaker;
class AliMUONCalibrationData;
class AliMUONQAMappingCheck;

class AliMUONTrackerQADataMakerRec: public AliMUONVQADataMakerRec {

public:
  AliMUONTrackerQADataMakerRec(AliQADataMakerRec* master);         
  virtual ~AliMUONTrackerQADataMakerRec();
  
  AliMUONVTrackerData* GetTrackerData() const;

  virtual void InitDigits(); 
  virtual void InitESDs(); 
  virtual void InitRaws(); 
  virtual void InitRecPoints(); 
  
  void EndOfDetectorCycleRaws(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleRecPoints(Int_t specie, TObjArray** list);
  void EndOfDetectorCycleESDs(Int_t specie, TObjArray** list);

  /// Empty implementation 
  void EndOfDetectorCycleDigits(Int_t, TObjArray**) {}
    
  virtual void MakeDigits(TTree* dig); 
  virtual void MakeESDs(AliESDEvent* esd) ;
  virtual void MakeRaws(AliRawReader* rawReader); 
  virtual void MakeRecPoints(TTree* recpo); 
  
public:

  /// Raw histograms indices
  enum ERaw { 
    kTrackerData              = 3,  ///< Accumulated data
    kTrackerBusPatchOccupancy = 4, ///< Bus patch occupancies
    kTrackerBusPatchNofPads   = 5, ///< Number of pads per bus patch
    kTrackerBusPatchNofManus  = 6, ///< Number of manus per bus patch
    kTrackerBusPatchConfig    = 7 ///< Configuration of the tracker
  };
         
  /// Rec points histograms indices
  enum ERecPoints { 
    kTrackerNumberOfClustersPerChamber    = 100, ///< Tracker: number of clusters per chamber
    kTrackerClusterMultiplicityPerChMean  = 101, ///< cluster size per Ch: mean
    kTrackerClusterMultiplicityPerChSigma = 102, ///< cluster size per Ch: dispersion
    kTrackerClusterChargePerChMean        = 103, ///< cluster charge per Ch: mean
    kTrackerClusterChargePerChSigma       = 104, ///< cluster charge per Ch: dispersion

    kTrackerRecPoints = 105, ///< Tracker : tracker data of clusters (all and mono-cathode ones)

    kTrackerClusterMultiplicityPerChamber = 200, ///< Tracker: cluster multiplicity per chamber
    kTrackerClusterChargePerChamber       = 300, ///< Tracker: cluster charge per chamber
    kTrackerClusterHitMapPerChamber       = 400, ///< Tracker: cluster position distribution per chamber
    
    kTrackerNumberOfClustersPerDE        = 1000, ///< Tracker : number of clusters per DE		
    kTrackerClusterMultiplicityPerDEMean = 1001, ///< cluster size per DE: mean
    kTrackerClusterChargePerDEMean       = 1002, ///< cluster charge per DE: mean
    
    kTrackerClusterMultiplicityPerDE = 3000, ///< Tracker : cluster multiplicity per DE		
    kTrackerClusterChargePerDE       = 5000  ///< Tracker : cluster charge per DE
    
  };
  
  /// ESD histograms indices
  enum EESD { 
    kESDnTracks                 = 0,  ///< number of tracks
    kESDMatchTrig               = 1,  ///< number of tracks matched with trigger
    kESDMomentum                = 2,  ///< P distribution
    kESDPt                      = 3,  ///< Pt distribution
    kESDRapidity                = 4,  ///< rapidity distribution
    kESDChi2                    = 5,  ///< normalized chi2 distribution
    kESDProbChi2                = 6,  ///< distribution of probability of chi2
    
    kESDClusterHitMap           = 7,  ///< cluster position distribution in chamber i
    kESDnClustersPerTrack       = 17, ///< number of clusters per track
    kESDnClustersPerCh          = 18, ///< number of clusters per chamber per track
    kESDnClustersPerDE          = 19, ///< number of clusters per DE per track
    kESDClusterChargeInCh       = 20, ///< cluster charge distribution in chamber i
    kESDClusterChargePerChMean  = 30, ///< cluster charge per Ch: mean
    kESDClusterChargePerChSigma = 31, ///< cluster charge per Ch: dispersion
    kESDClusterChargePerDE      = 32, ///< cluster charge per DE: mean
    kESDClusterSizeInCh         = 33, ///< cluster size distribution in chamber i
    kESDClusterSizePerChMean    = 43, ///< cluster size per Ch: mean
    kESDClusterSizePerChSigma   = 44, ///< cluster size per Ch: dispersion
    kESDClusterSizePerDE        = 45, ///< cluster size per DE: mean
    
    kESDResidualXInCh           = 46, ///< cluster-track residual-X distribution in chamber i
    kESDResidualYInCh           = 56, ///< cluster-track residual-Y distribution in chamber i
    kESDResidualXPerChMean      = 66, ///< cluster-track residual-X per Ch: mean
    kESDResidualYPerChMean      = 67, ///< cluster-track residual-Y per Ch: mean
    kESDResidualXPerChSigma     = 68, ///< cluster-track residual-X per Ch: dispersion
    kESDResidualYPerChSigma     = 69, ///< cluster-track residual-Y per Ch: dispersion
    kESDResidualXPerDEMean      = 70, ///< cluster-track residual-X per DE: mean
    kESDResidualYPerDEMean      = 71, ///< cluster-track residual-Y per DE: mean
    kESDResidualXPerDESigma     = 72, ///< cluster-track residual-X per DE: dispersion
    kESDResidualYPerDESigma     = 73, ///< cluster-track residual-Y per DE: dispersion
    kESDLocalChi2XInCh          = 74, ///< local chi2-X distribution in chamber i
    kESDLocalChi2YInCh          = 84, ///< local chi2-Y distribution in chamber i
    kESDLocalChi2XPerChMean     = 94, ///< local chi2-X per Ch: mean
    kESDLocalChi2YPerChMean     = 95, ///< local chi2-Y per Ch: mean
    kESDLocalChi2XPerDEMean     = 96, ///< local chi2-X per DE: mean
    kESDLocalChi2YPerDEMean     = 97, ///< local chi2-Y per DE: mean
    kESDLocalChi2InCh           = 98, ///< local chi2-X distribution in chamber i
    kESDLocalChi2PerChMean      = 108, ///< local chi2 per Ch: mean
    kESDLocalChi2PerDEMean      = 109, ///< local chi2 per DE: mean
    
    kESDThetaX                  = 110, ///< thetaX distribution
    kESDThetaY                  = 111, ///< thetaY distribution
    
    kESDnTotClustersPerCh       = 1000, ///< total number of associated clusters per chamber
    kESDnTotClustersPerDE       = 1001, ///< total number of associated clusters per DE
    kESDnTotFullClustersPerDE   = 1002, ///< total number of associated clusters containing pad info per DE
    kESDSumClusterChargePerDE   = 1003, ///< sum of cluster charge per DE
    kESDSumClusterSizePerDE     = 1004, ///< sum of cluster size per DE
    kESDSumResidualXPerDE       = 1005, ///< sum of cluster-track residual-X per DE
    kESDSumResidualYPerDE       = 1006, ///< sum of cluster-track residual-Y per DE
    kESDSumResidualX2PerDE      = 1007, ///< sum of cluster-track residual-X**2 per DE
    kESDSumResidualY2PerDE      = 1008, ///< sum of cluster-track residual-Y**2 per DE
    kESDSumLocalChi2XPerDE      = 1009, ///< sum of local chi2-X per DE
    kESDSumLocalChi2YPerDE      = 1010, ///< sum of local chi2-Y per DE
    kESDSumLocalChi2PerDE       = 1011  ///< sum of local chi2 per DE
  };
  
private:
  
  void InsertTrackerData(Int_t specie, TObjArray** list, TObject* object, 
                         Int_t indexNumber, Bool_t replace=kFALSE);

private:
  /// Not implemented
  AliMUONTrackerQADataMakerRec(const AliMUONTrackerQADataMakerRec& rhs);
  /// Not implemented
  AliMUONTrackerQADataMakerRec& operator=(const AliMUONTrackerQADataMakerRec& rhs);
  
  AliMUONVDigitStore*   fDigitStore; //!< pointer to digits store
  AliMUONDigitMaker*    fDigitMaker;  //!< pointer to digit maker
  AliMUONVClusterStore* fClusterStore; //!< pointer to cluster store
	
  AliMUONVTrackerDataMaker* fTrackerDataMaker; //!< tracker data accumulation (Raw)
  
  AliMUONQAMappingCheck* fMappingCheckRecPoints; //!< mapping cross-checker (RecPoints)
  
  AliMUONCalibrationData* fCalibrationData; //!< Used to load Local, Regional and Global masks
  
  ClassDef(AliMUONTrackerQADataMakerRec,1)  // MUON Quality assurance data maker

};
#endif
