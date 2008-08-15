// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCSLICETRACKERCOMPONENT_H
#define ALIHLTTPCSLICETRACKERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCSliceTrackerComponent.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  The TPC conformal mapping tracker component.
*/

#include "AliHLTProcessor.h"

class AliHLTTPCConfMapper;
class AliHLTTPCVertex;
class AliHLTTPCInterMerger;

/**
 * @class AliHLTTPCSliceTrackerComponent
 * The TPC conformal mapping tracker component.
 * 
 * Component ID: \b TPCSliceTracker <br>
 * Library: \b libAliHLTTPC.
 *
 * Mandatory arguments: <br>
 * \li -bfield     <i> magnetic field   </i> <br> 
 *      should allways be used. If there is no magnetic field, use 0.00001
 * 
 * Optional arguments: <br>
 * \li -pp-run <br>  
 *      will give fixed trackparameters for p-p runs.
 * \li -Pb-Pb-run <br>  
 *      will give fixed trackparameters for Pb-Pb runs.
 * \li -multiplicity     <i> multiplicity  </i> <br> 
 *      set the multiplicity of your collitions. Used for Pb-Pb events.
 * \li -disable-merger <br>
 *      turns off the intermerging on sectorlevel. Only use in spacial cases.
 * \li -nonvertextracking <br> 
 *      the tracking will do a reconstruction of clusters without vertex constraint.
 * \li -mainvertextrackingoff <br>
 *      will turn off the tracking with vertex constraint. Use together with -nonvertextracking.
 * \li -etarange   <i> eta range   </i> <br> 
 *      sets the eta range for the tracking.
 * \li -etasegment <i> eta segments   </i> <br> 
 *      let you set the number of eta segments in the trscker.
 * \li -chi2cut <i> chi²   </i> <br> 
 *      set the chi² for the tracker.
 * \li -rowscopetracklet <i> number of padrows   </i> <br> 
 *     tells the tracker how many padrows down it should look for the next cluster when making tracklets.
 * \li -rowscopetrack <i> number of padrows   </i> <br>
 *     tells the tracker how many padrows down it should look for the next cluster when making track.
 * \li -trackletlength <i> number of clusters   </i> <br>
 *     the minimum number of clusters to be on a tracklet.
 * \li -tracklength <i> number of clusters   </i> <br>
 *     the minimum number of clusters to be on a track.
 * \li -clusterZ <i> cutoff in z-direction (cm)   </i> <br>
 *     makes the tracker to not read in Clusters with a higher value then the
 *     one given in both directions
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCSliceTrackerComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTTPCSliceTrackerComponent();
  /** destructor */
  virtual ~AliHLTTPCSliceTrackerComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see @ref AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

protected:

  /**
   * Set Tracker parameters
   * Knowledge on that has been lost, so somebody has to insert more
   * documentation.
   */	
  void SetTrackerParam(Int_t phiSegments=50,Int_t etaSegments=100,
		       Int_t trackletlength=3,Int_t tracklength=5,
		       Int_t rowscopetracklet=2,Int_t rowscopetrack=3,
		       Double_t minPtFit=0,Double_t maxangle=1.31,
		       Double_t goodDist=5,Double_t hitChi2Cut=10,
		       Double_t goodHitChi2=20,Double_t trackChi2Cut=50,
		       Int_t maxdist=50,Double_t maxphi=0.1,Double_t maxeta=0.1);

  /**
   * Set Tracker parameters
   * Knowledge on that has been lost, so somebody has to insert more
   * documentation.
   */
  void SetTrackerParam( Bool_t doPP, Bool_t doPbPb, Int_t multiplicity, 
			Double_t bField, Int_t etasegment, Double_t hitchi2cut, 
			Int_t rowscopetracklet, Int_t rowscopetrack, 
			Int_t trackletlength, Int_t tracklength );

  /**
   * Set default tracker parameters
   * Knowledge on that has been lost, so somebody has to insert more
   * documentation.
   */
  void SetTrackerParam1();

  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 

  /** interface function, see @ref AliHLTComponent for description */
  int DoInit( int argc, const char** argv );
  /** interface function, see @ref AliHLTComponent for description */
  int DoDeinit();
  /** interface function, see @ref AliHLTComponent for description */
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
	       AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size, AliHLTComponentBlockDataList& outputBlocks );

  int Reconfigure(const char* cdbEntry, const char* chainId);
  int ReadPreprocessorValues(const char* modules);
  
  using AliHLTProcessor::DoEvent;
  
private:
  /** copy constructor prohibited */
  AliHLTTPCSliceTrackerComponent(const AliHLTTPCSliceTrackerComponent&);
  /** assignment operator prohibited */
  AliHLTTPCSliceTrackerComponent& operator=(const AliHLTTPCSliceTrackerComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */  
  int Configure(const char* arguments);

  /** instance of the tracker */
  AliHLTTPCConfMapper* fTracker;                                   //! transient
  /** vertex object */
  AliHLTTPCVertex* fVertex;                                        //! transient
  /** eta range */
  Float_t fEta[2];                                                 //  see above
  /** switch for subsequent non-vertex tracking */
  Bool_t fDoNonVertex;                                             //  see above
  /** */
  Bool_t  fDoPP;                                                   //  see above
   /** */
  Bool_t  fDoPbPb;                                                 //  see above
  /** multiplicity estimate */
  Int_t fMultiplicity;                                             //  see above
  /** magnetic field */
  Double_t fBField;                                                //  see above

  Bool_t fnonvertextracking;   // enable NONVERTEX Tracking
  Bool_t fmainvertextracking;  // enable MAINVERTEX Tracking

  /** ParameterVariables for the trackeparameters */


  /** phi_segments:     Devide the space into phi_segments */
  Int_t fPhisegment;              //! transient    
  /** ets_segments:     Devide the space into eta_segments */
  Int_t fEtasegment;              //! transient
  /** trackletlength:   Number of hits a tracklet has to have */
  Int_t fTrackletlength;          //! transient
  /** tracklength:      Number of hits a track has to have */
  Int_t fTracklength;             //! transient
  /** rowscopetracklet: Search range of rows for a tracklet */
  Int_t fRowscopetracklet;        //! transient
  /** rowscopetrack:    Search range of rows for a track */
  Int_t fRowscopetrack;           //! transient
  /** min_pt_fit:       Cut for moment fit, use:SetMaxDca(min_pt_fit) */
  Double_t fMinPtFit;             //! transient
  /** maxangle:         AliHLTTPCTransform::Deg2Rad(10), max angle for the three point look aheand */
  Double_t fMaxangle;             //! transient
  /** goodDist:         Threshold distancs between two hits when building tracklets */
  Double_t fGoodDist;             //! transient
  /** hitChi2Cut:       Max chi2 of added hit to track */
  Double_t fHitChi2Cut;           //! transient
  /** goodHitChi2:      Stop looking for next hit to add if chi2 is less then goodHitChi2 */
  Double_t fGoodHitChi2;          //! transient
  /** trackChi2Cut:     Max chi2 for track after final fit */
  Double_t fTrackChi2Cut;         //! transient
  /** maxdist:          Maximum distance between two clusters when forming segments */
  Int_t fMaxdist;                 //! transient
  /** maxphi:           Max phi difference for neighboring hits */
  Double_t fMaxphi;               //! transient
  /** maxeta:           Max eta difference for neighboring hits */
  Double_t fMaxeta;               //! transient 
 
  /** merger object */
  AliHLTTPCInterMerger *fpInterMerger;                             //! transient

  ClassDef(AliHLTTPCSliceTrackerComponent, 0);

};
#endif
