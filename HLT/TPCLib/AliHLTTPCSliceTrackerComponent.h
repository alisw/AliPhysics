// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCSLICETRACKERCOMPONENT_H
#define ALIHLTTPCSLICETRACKERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCSliceTrackerComponent.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  The TPC conformal mapping tracker component.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

class AliHLTTPCConfMapper;
class AliHLTTPCVertex;
class AliHLTTPCInterMerger;

/**
 * @class AliHLTTPCSliceTrackerComponent
 * The TPC conformal mapping tracker component.
 * 
 * The component has the following component arguments:
 * - disable-merger  disable merging of track segments
 * - pp-run          parameter set for pp run
 * - multiplicity    multiplicity to choose parameter set for
 * - bfield          magnatic field
 * 
 * @ingroup alihlt_tpc
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
		       Int_t maxdist=50,Double_t maxphi=0.1,Double_t maxeta=0.1,
		       bool vertexconstraint=true);

  /**
   * Set Tracker parameters
   * Knowledge on that has been lost, so somebody has to insert more
   * documentation.
   */
  void SetTrackerParam( bool doPP, int multiplicity, double bField );

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
  
private:
  /** copy constructor prohibited */
  AliHLTTPCSliceTrackerComponent(const AliHLTTPCSliceTrackerComponent&);
  /** assignment operator prohibited */
  AliHLTTPCSliceTrackerComponent& operator=(const AliHLTTPCSliceTrackerComponent&);

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
  /** multiplicity estimate */
  Int_t fMultiplicity;                                             //  see above
  /** magnetic field */
  Double_t fBField;                                                //  see above

  Bool_t fnonvertextracking;   // enable NONVERTEX Tracking
  Bool_t fmainvertextracking;  // enable MAINVERTEX Tracking

  /** merger object */
  AliHLTTPCInterMerger *fpInterMerger;                             //! transient

  ClassDef(AliHLTTPCSliceTrackerComponent, 0);

};
#endif
