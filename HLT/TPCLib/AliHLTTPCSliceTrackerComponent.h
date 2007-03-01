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

#include "AliHLTProcessor.h"
#include "AliHLTTPCDefinitions.h"

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
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTTPCSliceTrackerComponent(const AliHLTTPCSliceTrackerComponent&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTTPCSliceTrackerComponent& operator=(const AliHLTTPCSliceTrackerComponent&);
  /** destructor */
  virtual ~AliHLTTPCSliceTrackerComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	AliHLTComponent* Spawn();

protected:
	
	void SetTrackerParam(Int_t phiSegments=50,Int_t etaSegments=100,
			     Int_t trackletlength=3,Int_t tracklength=5,
			     Int_t rowscopetracklet=2,Int_t rowscopetrack=3,
			     Double_t minPtFit=0,Double_t maxangle=1.31,
			     Double_t goodDist=5,Double_t hitChi2Cut=10,
			     Double_t goodHitChi2=20,Double_t trackChi2Cut=50,
			     Int_t maxdist=50,Double_t maxphi=0.1,Double_t maxeta=0.1,
			     bool vertexconstraint=true);
	void SetTrackerParam( bool doPP, int multiplicity, double bField );
	void SetTrackerParam1()
		{
		SetTrackerParam( 10, 20, 5, 10, 2,2,
				 0, 1.31, 5, 100,
				 50, 100, 50, 0.1, 0.1,
				 true );
		}

	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
	
private:

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

// BEGINN ############################################## MODIFIY JMT
  Bool_t fnonvertextracking;   // enable NONVERTEX Tracking
  Bool_t fmainvertextracking;  // enable MAINVERTEX Tracking
// END ################################################# MODIFIY JMT

  /** merger object */
  AliHLTTPCInterMerger *fpInterMerger;                             //! transient

  ClassDef(AliHLTTPCSliceTrackerComponent, 0);

};
#endif
