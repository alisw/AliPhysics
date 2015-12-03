//
//  AliEveKineTracks.h
//
//  Created by Jeremi Niedziela on 3/12/15
//  Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007
//

#ifndef __AliEveKineTracks__
#define __AliEveKineTracks__

#include <AliEveTrack.h>
#include <AliStack.h>

class AliEveKineTracks
{
public:
    AliEveKineTracks(){}
    ~AliEveKineTracks(){}
    
    TEveTrackList* Draw(Double_t min_pt  = 0,     Double_t min_p   = 0,
                        Bool_t   pdg_col = kTRUE, Bool_t   recurse = kTRUE,
                        Bool_t   use_track_refs = kTRUE);

private:
    TEveElement* KineTrack(Int_t  label,
                           Bool_t import_mother = kTRUE, Bool_t import_daughters = kTRUE,
                           Bool_t pdg_col       = kTRUE, Bool_t recurse          = kTRUE,
                           TEveElement* cont = 0);
    
    void KineDaughters(AliEveTrack* parent,  AliStack* stack,
                                         Double_t     min_pt,  Double_t  min_p,
                                         Bool_t       pdg_col, Bool_t    recurse);
    
    void SetTrackColor(AliEveTrack* t, Bool_t pdg_col);
    Color_t GetPDGColor(Int_t pdg);
    void HideNeutrals(TEveElement* el, Int_t level);
    
    AliEveKineTracks(const AliEveKineTracks&);
    AliEveKineTracks& operator=(const AliEveKineTracks&);
};

#endif
