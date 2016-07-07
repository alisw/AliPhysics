/// \class AliTrackContainerV0
/// \brief Select tracks based on specific prescriptions of V0s in jets analysis
///
/// This class derives from AliParticleContainer. 
/// It allows to select tracks based on specific prescriptions of analysis
/// of V0 particles associated with jets.
/// Namely, it remove the V0 daughter tracks from track sample prior
/// to jet finding procedure.
///
/// \author Vojtech Pacik <vojtech.pacik@cern.ch>, NPI Czech Academy of Sciences
/// \date Jun 24, 2016

#ifndef ALITRACKCONTAINERV0_H
#define ALITRACKCONTAINERV0_H

/**************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

class AliAODEvent;

#include "AliAODEvent.h"
#include "AliTrackContainer.h"

class AliTrackContainerV0 : public AliTrackContainer {
	public:
		AliTrackContainerV0();
  	AliTrackContainerV0(const char *name);

  	void   SetFilterDaughterTracks(Bool_t bFilter)       { fFilterDaughterTracks = bFilter; }
  	Bool_t GetFilterDaughterTracks()               const { return fFilterDaughterTracks   ; }

  	// reimplementation of inherited methods
  	virtual void    SetArray(const AliVEvent *event);
  	virtual void    NextEvent();
  	virtual Bool_t	ApplyTrackCuts(const AliVTrack* vp, UInt_t &rejectionReason) const;
		
		void ExtractDaughters(AliAODv0* cand);

	protected:	
		Bool_t	IsV0Daughter(const AliAODTrack* track) const;

		Bool_t              fFilterDaughterTracks    ; ///< if the daughter tracks of V0s  candidates should be filtered out
		const AliAODEvent  *fEvent                   ; ///< pointer to current event (pointer stay the same, but the content is changed event-by-event)
		TObjArray          *fV0s                     ; ///< list of V0 candidates
		TObjArray	          fDaughterList            ; ///< list of V0 daughters

	private:
	/// \cond CLASSIMP
  ClassDef(AliTrackContainerV0, 1);
  /// \endcond
};
#endif
