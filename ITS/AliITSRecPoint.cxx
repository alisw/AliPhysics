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

/*
$Log$
Revision 1.3.12.1  2002/10/14 13:14:08  hristov
Updating VirtualMC to v3-09-02

Revision 1.4  2002/09/09 17:28:02  nilsen
Added class iostreamer funcionality and Print and Read functions. cleaned
up files.

*/
////////////////////////////////////////////////
//  Reconstructed point class for set:ITS     //
////////////////////////////////////////////////


#include "AliITSRecPoint.h"
ClassImp(AliITSRecPoint)

AliITSRecPoint::AliITSRecPoint() {
    // default creator
    fTracks[0]=fTracks[1]=fTracks[2]=-3; 
    fX=fZ=fQ=fdEdX=0.;
    fSigmaX2=fSigmaZ2=0.;
}
//----------------------------------------------------------------------
void AliITSRecPoint::Print(ostream *os){
    ////////////////////////////////////////////////////////////////////////
    // Standard output format for this class.
    ////////////////////////////////////////////////////////////////////////
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif
 
    fmt = os->setf(ios::fixed);  // set fixed floating point output
    *os << fTracks[0]<< " " << fTracks[1] << " " << fTracks[2] << " ";
    *os << fX << " " << fZ << " " << fQ << " ";
    fmt = os->setf(ios::scientific); // set scientific for dEdX.
    *os << fdEdX << " ";
    fmt = os->setf(ios::fixed); // every fixed
    *os << fSigmaX2 << " " << fSigmaZ2;
    os->flags(fmt); // reset back to old formating.
    return;
}
//----------------------------------------------------------------------
void AliITSRecPoint::Read(istream *is){
////////////////////////////////////////////////////////////////////////
// Standard input format for this class.
////////////////////////////////////////////////////////////////////////
 

    *is >> fTracks[0] >> fTracks[1] >> fTracks[2] >> fX >> fZ >> fQ;
    *is >> fdEdX >> fSigmaX2 >> fSigmaZ2;
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSRecPoint &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////
 
    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSRecPoint &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////
 
    r.Read(&is);
    return is;
}
//----------------------------------------------------------------------
