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

class AliExternalTrackParam;
class AliESDtrack;
class AliAODTrack;

//______________________________________________________________________________
// Short description of AliEveTrack
//

class AliEveTrack : public TEveTrack
{
public:
    AliEveTrack();
    AliEveTrack(TParticle* t, Int_t label, TEveTrackPropagator* prop=0);
    AliEveTrack(TEveMCTrack*  t, TEveTrackPropagator* prop=0);
    AliEveTrack(TEveRecTrack* t, TEveTrackPropagator* prop=0);
    AliEveTrack(AliESDtrack*  t, TEveTrackPropagator* prop=0);
    AliEveTrack(AliAODTrack*  t, TEveTrackPropagator* prop=0);
    AliEveTrack(const AliEveTrack& t);
    virtual ~AliEveTrack();
    
    void SetStartParams(const AliExternalTrackParam* tp);
    
    void ImportHits();      // *MENU*
    
    void ImportClustersFromLabel(); // *MENU*
    void ImportClustersFromIndex(); // *MENU*
    TEvePointSet* ImportClustersFromIndex(Int_t index);
    
    void ImportKine();              // *MENU*
    void ImportKineWithArgs(Bool_t importMother=kTRUE, Bool_t impDaugters=kTRUE,
                            Bool_t colorPdg    =kTRUE, Bool_t recurse    =kTRUE); // *MENU*
    void PrintKineStack();          // *MENU*
    
    virtual void SecSelected(TEveTrack*);        // *SIGNAL*
    virtual void SecSelectedTrack(AliEveTrack*); // *SIGNAL*
    
    AliESDtrack* GetESDTrack() const;
    AliAODTrack* GetAODTrack() const;
    
protected:
    
private:
    AliEveTrack& operator=(const AliEveTrack&); // Not implemented
    
    ClassDef(AliEveTrack, 0); // Short description.
};

#endif
