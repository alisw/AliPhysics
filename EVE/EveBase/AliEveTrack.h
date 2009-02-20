// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTrack_H
#define AliEveTrack_H

#include <TEveTrack.h>

//______________________________________________________________________________
// Short description of AliEveTrack
//

class AliEveTrack : public TEveTrack
{
public:
  AliEveTrack();
  AliEveTrack(TParticle* t, Int_t label, TEveTrackPropagator* rs);
  AliEveTrack(TEveMCTrack*  t, TEveTrackPropagator* rs);
  AliEveTrack(TEveRecTrack* t, TEveTrackPropagator* rs);
  AliEveTrack(const AliEveTrack& t);
  virtual ~AliEveTrack();

  virtual void ImportHits();      // *MENU*
  virtual void ImportClusters();  // *MENU*

  void ImportClustersFromIndex(); // *MENU*
  void ImportKine();              // *MENU*
  void ImportKineWithArgs(Bool_t importMother=kTRUE, Bool_t impDaugters=kTRUE,
			  Bool_t colorPdg    =kTRUE, Bool_t recurse    =kTRUE); // *MENU*
  void PrintKineStack();          // *MENU*

protected:

private:
  AliEveTrack& operator=(const AliEveTrack&); // Not implemented

  ClassDef(AliEveTrack, 0); // Short description.
};

#endif
