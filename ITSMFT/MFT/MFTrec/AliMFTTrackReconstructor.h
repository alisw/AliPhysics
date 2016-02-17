#ifndef AliMFTTrackReconstructor_H
#define AliMFTTrackReconstructor_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup MFTrec
/// \class AliMFTTrackReconstructor
/// \brief Class Doing MFT Track reconstruction
///
///
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>, IPN-Lyon
/// \date Octerber 1st, 2015

#include "TObject.h"

class AliMFTTrackParam;
class AliMFTTrack;
class AliMFTCACell;
//=============================================================================================

class AliMFTTrackReconstructor : public TObject {
	
public:
	
	AliMFTTrackReconstructor();
	virtual ~AliMFTTrackReconstructor();
	void EventReconstruct(TClonesArray *fMFTTracks);
	Bool_t AddMCSEffect(AliMFTCACell *currentCell,  AliMFTTrackParam *trackParam);

	
protected:
	
	Bool_t TraceTrack(AliMFTTrack *fMFTTracks );
	Double_t RunKalmanFilter(AliMFTTrackParam &trackParamAtCluster);

	
	/// \cond CLASSIMP
	ClassDef(AliMFTTrackReconstructor,1);
	/// \endcond
};

//=============================================================================================

#endif
