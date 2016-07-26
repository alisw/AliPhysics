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
    fE2(0),
    fArea(0),
    fConstituents(),
    fConstituentsRef()
{;}


AliJJet::AliJJet(float px,float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge):
    AliJBaseTrack( px, py, pz, e, id, ptype, charge ),
    fLeadingTrackId(-1),
    fLeadingTrackPt(-1),
    fLeadingTrackE(-1),
    fE2(0),
    fArea(0),
    fConstituents(),
    fConstituentsRef()
{;}

AliJJet::AliJJet(const AliJJet& a):
    AliJBaseTrack( a ),
    fLeadingTrackId( a.fLeadingTrackId),
    fLeadingTrackPt( a.fLeadingTrackPt),
    fLeadingTrackE( a.fLeadingTrackE),
    fE2(a.fArea),
    fArea( a.fArea ),
    fConstituents( a.fConstituents),
    fConstituentsRef( a.fConstituentsRef )
{;}

AliJJet::AliJJet(const TLorentzVector & a):
    AliJBaseTrack( a ),
    fLeadingTrackId(-1),
    fLeadingTrackPt(-1),
    fLeadingTrackE(-1),
    fE2(0),
    fArea(0),
    fConstituents(),
    fConstituentsRef()
{;}


AliJJet& AliJJet::operator=(const AliJJet& trk){
  //operator =  
  if(this != &trk){
    AliJBaseTrack::operator=(trk);
    fArea = trk.fArea;
    fConstituents = trk.fConstituents;
    fConstituentsRef = trk.fConstituentsRef;
    fLeadingTrackId = trk.fLeadingTrackId;
    fLeadingTrackPt = trk.fLeadingTrackPt;
    fLeadingTrackE = trk.fLeadingTrackE;
    fE2 = trk.fE2;
  }
  return *this;
}

void AliJJet::ReSum(){
    TLorentzVector lv;
    int lid = -1;
    double lpt = 1e-10;
    int lidE = -1;
    double lE = 1e-10;
    fE2 = E();
    GenConstituentsFromRef();
    if( GetNConstituents() < 1 ) return;
    fNConstituent = fConstituents.GetEntriesFast();
    double E2 = -1;
    for( int i=0;i<GetNConstituents();i++ ){
        AliJBaseTrack * trk = (AliJBaseTrack*) fConstituents.At(i);
        if( !trk ){ 
            continue;
        }
        TLorentzVector * v = (TLorentzVector*)fConstituents.At(i);
        lv += *v;
        double e2 = TMath::Sqrt(TMath::Power(v->P(),2)+TMath::Power(0.1349766,2)) ;
        E2 += e2;
        double pt = trk->Pt();
        if( pt > lpt ){lpt = pt;lid = i;}
        double e = trk->E();
        if( e > lE ){lE = e;lidE = i;}
    }
    fLeadingTrackPt = lpt;
    fLeadingTrackE = lE;
    if( lv.E() < 1e-4 ) return;
    if( lv.E() > 0 ) TLorentzVector::operator=(lv);
    if( E2 < 0 ) E2 = E();
    else E2+=1;
    fE2=E2;
}

ClassImp(AliJJet)
