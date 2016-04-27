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

#include "AliJTrackCut.h"
#include "AliJTrack.h"
#include <TMath.h>

AliJTrackCut::AliJTrackCut()
{
}

/* FOR ESD
bool AliJTrackCut::IsSelected ( AliJTrack *track, int icut )const{
 // Comment needed
    if( ! track ) return false;
    if( 1 ) {// TODO fProdVersion == 10 )
        switch ( icut ){

            case AliJTrackCut::kJTPCOnly :
                return track->IsFiltered( 7 );

            case AliJTrackCut::kJRaa :
                return track->IsFiltered( 10 ); //TODO

            case AliJTrackCut::kJTPCplus :
                return track->IsFiltered( 12 );

            case AliJTrackCut::kJHybridA :
                return track->IsFiltered( 14 );

            case AliJTrackCut::kJHybridB :
                return  ( !track->IsFiltered( 14 ) && track->IsFiltered(15) );

            case AliJTrackCut::kJHybrid :
                return ( IsSelected(track,kJHybridA) || IsSelected(track, kJHybridB ) );

            case AliJTrackCut::kJTPC2A :
                return track->IsFiltered(11);
            case AliJTrackCut::kJRaa2011 :
                return track->IsFiltered(13);
            case AliJTrackCut::kJH22011 :
                return track->IsFiltered(16);
            case AliJTrackCut::kJH32011 :
                return track->IsFiltered(17);

        }
    }
    return false;
}
  */

//==== FOR AOD
bool AliJTrackCut::IsSelected ( AliJTrack *track, int icut )const{
  // comment needed
    if( ! track ) return false;
    if( 1 ) {// TODO fProdVersion == 10 )
        switch ( icut ){

            case AliJTrackCut::kJTPCOnly :
                return track->IsFiltered( 7 );  // 128

            case AliJTrackCut::kJRaa :
                return track->IsFiltered( 10 ); // 1024

            case AliJTrackCut::kJHybrid :
                return track->IsFiltered(8) || track->IsFiltered(9) ; // 768

            case AliJTrackCut::kJHybridAOD86 :
                return track->IsFiltered(4) || track->IsFiltered(8) ; // 272

            case AliJTrackCut::kJGlobalTightDCA :
                return track->IsFiltered( 5 ); // 32

            case AliJTrackCut::kJGlobalDCA :
                return track->IsFiltered( 4 ); // 16

            case AliJTrackCut::kJGlobalSDD :
                return track->IsFiltered( 5 ) || track->IsFiltered( 6 ); // 96
        }
    }
    return false;
}

bool AliJTrackCut::SetMomentum( AliJTrack *track, int icut ) const {
  // comment needed
    if( ! track ) return false;
    if( 1 ) {// TODO fProdVersion == 10 )
        
        switch ( icut ){
        /*
            case AliJTrackCut::kJTPCOnly : 
                if( track->PtTPC() > 1e-9 ){
                    track->SetUseTPCTrack();
                }
                break;
            case AliJTrackCut::kJTPCplus :
                if( track->PtTPC() > 1e-9 ){
                    track->SetUseTPCTrack();
                }
                break;
            case AliJTrackCut::kJHybrid :
                if( (!IsSelected( track, kJHybridA )) && track->PtGCG() > 1e-9 )
                    track->SetUseGCGTrack();
                break;
        */
            default:
                break;
        }
        
    }
    return true;
}

AliJTrackCut& AliJTrackCut::GetSpecialInstance(){
  // comment needed
    static AliJTrackCut instance;
    return instance;
}
const AliJTrackCut& AliJTrackCut::GetInstance(){
  // comment needed
    return GetSpecialInstance();
}
