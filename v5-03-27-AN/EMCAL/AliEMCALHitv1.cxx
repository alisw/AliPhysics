 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: */

//_________________________________________________________________________
//  Hits class for EMCAL    
//  A hit in EMCAL is the usuall hi with full information about xlocal, ylocal, zlocal 
//*-- Author: Aleksei Pavlinov(WSU); Nov 25, 05

#include <Riostream.h>

// --- AliRoot header files ---
#include "AliEMCALHitv1.h"

ClassImp(AliEMCALHitv1)

//______________________________________________________________________
AliEMCALHitv1::AliEMCALHitv1(){
    // Default ctor
   
    fId      = 0;
    fELOS    = 0.0;
    fTime    = 0.0;
    fPrimary = 0;
    fTrack   = 0;
    fX       = 0.0;
    fY       = 0.0;
    fZ       = 0.0;
    fPx      = 0.0;
    fPy      = 0.0;
    fPz      = 0.0;
    fPe      = 0.0;
    fIparent = 0;
    fIenergy = 0.0;
}
//______________________________________________________________________
AliEMCALHitv1::AliEMCALHitv1(const AliEMCALHitv1 & hit) : AliHit(hit){
    // copy ctor
   
    fId      = hit.fId ; 
    fELOS    = hit.fELOS ;
    fPrimary = hit.fPrimary ; 
    fTrack   = hit.fTrack ; 
    fX       = hit.fX;
    fY       = hit.fY;
    fZ       = hit.fZ;
    fPx      = hit.fPx;
    fPy      = hit.fPy;
    fPz      = hit.fPz;
    fPe      = hit.fPe;
    fIparent = hit.fIparent;
    fIenergy = hit.fIenergy;
    fTime    = hit.fTime  ;
}
//______________________________________________________________________
AliEMCALHitv1::AliEMCALHitv1(Int_t shunt, Int_t primary, Int_t track,Int_t iparent, Float_t ienergy, Int_t id,
			 Float_t *hits,Float_t *p):AliHit(shunt, track){
    //
    // Create an EMCAL  hit object
    //
    fX          = hits[0];
    fY          = hits[1];
    fZ          = hits[2];
    fTime       = hits[3] ;
    fId         = id;
    fELOS       = hits[4];
    fPrimary    = primary;
    fPx         = p[0];
    fPy         = p[1];
    fPz         = p[2];
    fPe         = p[3];
    fIparent    = iparent;
    fIenergy    = ienergy;
}

//______________________________________________________________________
Bool_t AliEMCALHitv1::operator==(AliEMCALHitv1 const &rValue) const{ 
  // no identical hits !!
    Bool_t rv = kFALSE;

    if ( (fId == rValue.GetId()));
    //    if ( (fId == rValue.GetId()) && ( fIparent == rValue.GetIparent()))
    //	rv = kTRUE;

    return rv;
}
//______________________________________________________________________
AliEMCALHitv1 AliEMCALHitv1::operator+(const AliEMCALHitv1 &rValue){
    // Add the energy of the hit

    fELOS += rValue.GetEnergy() ;
 
    if(rValue.GetTime() < fTime)
      fTime = rValue.GetTime() ;
 
    return *this;

}
//______________________________________________________________________
ostream& operator << (ostream& out,AliEMCALHitv1& hit){
    // Print out Id and energy

    out << "AliEMCALHitv1:";
    out << "id=" <<  hit.GetId();
    out << ", Eloss=" <<  hit.GetEnergy();
    out << ", Time=" << hit.GetTime();
    out << "GeV , Track no.=" << hit.GetPrimary();
    out << ", (xyz)=(" << hit.X()<< ","<< hit.Y()<< ","<<hit.Z()<<") cm <- local";
    out << ", fTrack=" << hit.GetTrack();
    out << ", P=(" << hit.GetPx() << "," << hit.GetPy() << "," << hit.GetPz()
                  << "," <<hit.GetPe() << ") GeV"  ;
    out << ", Enterring particle ID" << hit.GetIparent();
    out << ", Enterring particle initial energy = " << hit.GetIenergy() << " GeV" ;
    out << endl;

    return out;
}
