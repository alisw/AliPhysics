/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

#include "AliJJet.h"

AliJJet::AliJJet(): 
    AliJBaseTrack(), 
    fLeadingTrackId(-1),
    fLeadingTrackPt(-1),
    fLeadingTrackE(-1),
    fNConstituent(0),
    fArea(0),
    fConstituents()
{;}


AliJJet::AliJJet(float px,float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge):
    AliJBaseTrack( px, py, pz, e, id, ptype, charge ),
    fLeadingTrackId(-1),
    fLeadingTrackPt(-1),
    fLeadingTrackE(-1),
    fNConstituent(0),
    fArea(0),
    fConstituents()
{;}

AliJJet::AliJJet(const AliJJet& a):
    AliJBaseTrack( a ),
    fLeadingTrackId( a.fLeadingTrackId),
    fLeadingTrackPt( a.fLeadingTrackPt),
    fLeadingTrackE( a.fLeadingTrackE),
    fNConstituent(a.fNConstituent),
    fArea( a.fArea ),
    fConstituents( a.fConstituents )
{;}

AliJJet::AliJJet(const TLorentzVector & a):
    AliJBaseTrack( a ),
    fLeadingTrackId(-1),
    fLeadingTrackPt(-1),
    fLeadingTrackE(-1),
    fNConstituent(0),
    fArea(0),
    fConstituents()
{;}


AliJJet& AliJJet::operator=(const AliJJet& trk){
  //operator =  
  if(this != &trk){
    AliJBaseTrack::operator=(trk);
    fArea = trk.fArea;
    fConstituents = trk.fConstituents;
    fLeadingTrackId = trk.fLeadingTrackId;
    fLeadingTrackPt = trk.fLeadingTrackPt;
    fLeadingTrackE = trk.fLeadingTrackE;
  }
  return *this;
}

void AliJJet::ReSum(){
    TLorentzVector lv;
    int lid = -1;
    double lpt = 1e-10;
    int lidE = -1;
    double lE = 1e-10;
    fNConstituent = fConstituents.GetEntriesFast();
    for( int i=0;i<GetNConstituents();i++ ){
        AliJBaseTrack * trk = (AliJBaseTrack*) fConstituents[i];
        if( !trk ){ 
            cout<<"DEBUG E1 : No trk in "<<endl;
        continue;
        }
        TLorentzVector * v = (TLorentzVector*)fConstituents[i];
        v->SetE( v->Perp() );
        lv += *v;
        double pt = trk->Pt();
        if( pt > lpt ){lpt = pt;lid = i;}
        double e = trk->E();
        if( e > lE ){lE = e;lidE = i;}
    }
	fLeadingTrackId = lid;
    fLeadingTrackPt = lpt;
    fLeadingTrackE = lE;
    TLorentzVector::operator=(lv);
}

ClassImp(AliJJet)
