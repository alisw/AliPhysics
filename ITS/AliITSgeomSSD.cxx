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
Revision 1.8  2000/10/02 16:32:43  barbera
Forward declaration added

Revision 1.2.4.8  2000/10/02 15:53:49  barbera
Forward declaration added

Revision 1.7  2000/07/10 16:07:18  fca
Release version of ITS code

Revision 1.2.4.2  2000/03/04 23:55:59  nilsen
Fixed up the comments/documentation

Revision 1.2.4.1  2000/01/12 19:03:32  nilsen
This is the version of the files after the merging done in December 1999.
See the ReadMe110100.txt file for details

Revision 1.2  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/
#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <TShape.h>
#include <TBRIK.h>

#include "AliITSgeomSSD.h"

ClassImp(AliITSgeomSSD)
AliITSgeomSSD::AliITSgeomSSD(const Float_t *box,Float_t ap,Float_t an,
			     Int_t np,Float_t *p,Int_t nn,Float_t *n){
////////////////////////////////////////////////////////////////////////
//    Standard Constructor. *box={dx,dy,dz}, ap=anode angle, an=cathode angle,
// nn= number of cathodes+1,*n= array of cathode low edges+highest edge,
// np= number of anodes+1, *p= array of anode low edges+lighest edge.
///////////////////////////////////////////////////////////////////////
    Int_t i;

    fShapeSSD = new TBRIK("ActiveSSD","Active volume of SSD","SSD SI DET",
			  box[0],box[1],box[2]);
    if(fLowEdgeP!=0) delete fLowEdgeP;
    if(fLowEdgeN!=0) delete fLowEdgeN;
    fNp = np;
    fNn = nn;
    fAngleP = ap;
    fAngleN = an;
    fLowEdgeP = new Float_t[fNp];
    fLowEdgeN = new Float_t[fNn];
    for(i=0;i<fNp;i++) fLowEdgeP[i] = p[i];
    for(i=0;i<fNn;i++) fLowEdgeN[i] = n[i];
}
//______________________________________________________________________
AliITSgeomSSD::~AliITSgeomSSD(){
    // Destructor.

    delete fLowEdgeP; fLowEdgeP = 0;
    delete fLowEdgeN; fLowEdgeN = 0;
    delete fShapeSSD; fShapeSSD = 0;
    fNp = 0;
    fNn = 0;
    fAngleP = 0.0;
    fAngleN = 0.0;
}
AliITSgeomSSD::AliITSgeomSSD(const AliITSgeomSSD &source){
////////////////////////////////////////////////////////////////////////
//    copy  constructor
////////////////////////////////////////////////////////////////////////
    Int_t i;

    if(this == &source) return;
    this->fShapeSSD = new TBRIK(*(source.fShapeSSD));
    this->fNp = source.fNp;
    this->fNn = source.fNn;
    delete fLowEdgeP;
    delete fLowEdgeN;
    this->fAngleP = source.fAngleP;
    this->fAngleN = source.fAngleN;
    fLowEdgeP = new Float_t[fNp];
    fLowEdgeN = new Float_t[fNn];
    for(i=0;i<fNp;i++) this->fLowEdgeP[i] = source.fLowEdgeP[i];
    for(i=0;i<fNn;i++) this->fLowEdgeN[i] = source.fLowEdgeN[i];
    return;
}  

AliITSgeomSSD& AliITSgeomSSD::operator=(const AliITSgeomSSD &source) {
////////////////////////////////////////////////////////////////////////
//    assignment operator
////////////////////////////////////////////////////////////////////////
    Int_t i;

    if(this == &source) return *this;
    this->fShapeSSD = new TBRIK(*(source.fShapeSSD));
    this->fNp = source.fNp;
    this->fNn = source.fNn;
    delete fLowEdgeP;
    delete fLowEdgeN;
    this->fAngleP = source.fAngleP;
    this->fAngleN = source.fAngleN;
    fLowEdgeP = new Float_t[fNp];
    fLowEdgeN = new Float_t[fNn];
    for(i=0;i<fNp;i++) this->fLowEdgeP[i] = source.fLowEdgeP[i];
    for(i=0;i<fNn;i++) this->fLowEdgeN[i] = source.fLowEdgeN[i];
    return *this;
}
//______________________________________________________________________
void AliITSgeomSSD::Local2Det(Float_t x,Float_t z,Int_t &a,Int_t &c){
    Float_t d,b;
    Int_t i;

    // project on to bonding edges.
    d = x*TMath::Cos(fAngleP)+z*TMath::Sin(fAngleP);
    b = x*TMath::Cos(fAngleN)+z*TMath::Sin(fAngleN);
    if(d<fLowEdgeP[0]) i=-1;
    else for(i=0;i<fNp;i++){
	if(fLowEdgeP[i]<d) break;
    } // end for i
    a = i;
    if(b<fLowEdgeN[0]) i=-1;
    else for(i=0;i<fNn;i++){
	if(fLowEdgeN[i]<b) break;
    } // end for i
    c = i;
    return;
}
//______________________________________________________________________
void AliITSgeomSSD::Det2Local(Int_t a,Int_t c,Float_t &x,Float_t &z){
//    Float_t d,b;
//    Int_t i;

    return;
}
//______________________________________________________________________
void AliITSgeomSSD::Print(ostream *os){
////////////////////////////////////////////////////////////////////////
// Standard output format for this class.
////////////////////////////////////////////////////////////////////////
    ios::fmtflags fmt;
    Int_t i;

    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "TBRIK" << " ";
    *os << setprecision(16) << GetDx() << " ";
    *os << setprecision(16) << GetDy() << " ";
    *os << setprecision(16) << GetDz() << " ";
    *os << fNp << " " << fNn << " ";
    *os << setprecision(16) << fAngleP << " ";
    *os << setprecision(16) << fAngleN << " ";
    for(i=0;i<fNp;i++) *os << setprecision(16) << fLowEdgeP[i] << " ";
    for(i=0;i<fNn;i++) *os << setprecision(16) << fLowEdgeN[i] << " ";
    *os << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSgeomSSD::Read(istream *is){
////////////////////////////////////////////////////////////////////////
// Standard input format for this class.
////////////////////////////////////////////////////////////////////////
    Float_t dx,dy,dz;
    Int_t i;
    char shp[20];

    *is >> shp;
    *is >> dx >> dy >> dz;
    if(fShapeSSD!=0) delete fShapeSSD;
    fShapeSSD = new TBRIK("ActiveSSD","Active volume of SSD","SSD SI DET",
			    dx,dy,dz);
    *is >> fNp >> fNn;
    *is >> fAngleP >> fAngleN;
    if(fLowEdgeP !=0) delete fLowEdgeP;
    if(fLowEdgeN !=0) delete fLowEdgeN;
    fLowEdgeP = new Float_t[fNp];
    fLowEdgeN = new Float_t[fNn];
    for(i=0;0<fNp;i++) *is >> fLowEdgeP[i];
    for(i=0;i<fNn;i++) *is >> fLowEdgeN[i];
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSgeomSSD &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSSD &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}
//======================================================================
/*
$Log$
*/

//#include "AliITSgeomSSD175.h"

ClassImp(AliITSgeomSSD175)

AliITSgeomSSD175::AliITSgeomSSD175() : AliITSgeomSSD(){
////////////////////////////////////////////////////////////////////////
//    default constructor
////////////////////////////////////////////////////////////////////////
    const Float_t kDxyz[] ={3.6500,0.0150,2.000};//cm. (Geant 3.12 units)
    // Size of sensitive detector area x,y(thickness),z
    const Float_t kangle   = 0.0175; // angle in rad. of anode and cathodes
    const Float_t kpitch   = 0.0095;// cm anode separation.
    const Int_t   kNstrips = 768; // number of anode or cathode strips.
    Float_t *leA,*leC; // array of low edges anode and cathorde.
    Int_t i;

    leA = new Float_t[kNstrips+1];
    leC = new Float_t[kNstrips+1];
    leA[0] = -kDxyz[0];
    leA[1] = -kpitch*(0.5*kNstrips-1);
    leC[0] =  kDxyz[0];
    leC[1] =  kpitch*(0.5*kNstrips-1);
    for(i=1;i<kNstrips;i++){
	leA[i+1] = leA[i] + kpitch;
	leC[i+1] = leC[i] - kpitch;
    } // end for i
    leA[kNstrips] =  kDxyz[0];
    leC[kNstrips] = -kDxyz[0];
//    cout << "AliITSgeomSSD175 default creator called: start" << endl;
    AliITSgeomSSD::AliITSgeomSSD(kDxyz,kangle,-kangle,
				 kNstrips+1,leA,kNstrips+1,leC);
    delete leA;
    delete leC;
//    cout << "AliITSgeomSSD175 default creator called: end" << endl;
}
//________________________________________________________________________
ostream &operator<<(ostream &os,AliITSgeomSSD175 &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSSD175 &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}
//======================================================================
/*
$Log$
*/

//#include "AliITSgeomSSD275and75.h"

ClassImp(AliITSgeomSSD275and75)

AliITSgeomSSD275and75::AliITSgeomSSD275and75() : AliITSgeomSSD(){
////////////////////////////////////////////////////////////////////////
//    default constructor
////////////////////////////////////////////////////////////////////////
    const Float_t kDxyz[] ={3.6500,0.0150,2.000};//cm. (Geant 3.12 units)
    // Size of sensitive detector area x,y(thickness),z
    const Float_t kangleA  = 0.0275; // angle in rad. of anode and cathodes
    const Float_t kangleC  = 0.0075; // angle in rad. of anode and cathodes
    const Float_t kpitch   = 0.0095;// cm anode separation.
    const Int_t   kNstrips = 768; // number of anode or cathode strips.
    Float_t *leA,*leC; // array of low edges anode and cathorde.
    Int_t i;

    leA = new Float_t[kNstrips+1];
    leC = new Float_t[kNstrips+1];
    leA[0] = -kDxyz[0];
    leA[1] = -kpitch*(0.5*kNstrips-1);
    leC[0] =  kDxyz[0];
    leC[1] =  kpitch*(0.5*kNstrips-1);
    for(i=1;i<kNstrips;i++){
	leA[i+1] = leA[i] + kpitch;
	leC[i+1] = leC[i] - kpitch;
    } // end for i
    leA[kNstrips] =  kDxyz[0];
    leC[kNstrips] = -kDxyz[0];
//    cout << "AliITSgeomSSD275and75 default creator called: start" << endl;
    AliITSgeomSSD::AliITSgeomSSD(kDxyz,kangleA,kangleC,
				 kNstrips+1,leA,kNstrips+1,leC);
    delete leA;
    delete leC;
//    cout << "AliITSgeomSSD275and75 default creator called: end" << endl;
}
//________________________________________________________________________
ostream &operator<<(ostream &os,AliITSgeomSSD275and75 &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSSD275and75 &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}
//======================================================================
/*
$Log$
*/
//#include "AliITSgeomSSD75and275.h"

ClassImp(AliITSgeomSSD75and275)

AliITSgeomSSD75and275::AliITSgeomSSD75and275() : AliITSgeomSSD(){
////////////////////////////////////////////////////////////////////////
//    default constructor
////////////////////////////////////////////////////////////////////////
    const Float_t kDxyz[] ={3.6500,0.0150,2.000};//cm. (Geant 3.12 units)
    // Size of sensitive detector area x,y(thickness),z
    const Float_t kangleA  = 0.0075; // angle in rad. of anode and cathodes
    const Float_t kangleC  = 0.0275; // angle in rad. of anode and cathodes
    const Float_t kpitch   = 0.0095;// cm anode separation.
    const Int_t   kNstrips = 768; // number of anode or cathode strips.
    Float_t *leA,*leC; // array of low edges anode and cathorde.
    Int_t i;

    leA = new Float_t[kNstrips+1];
    leC = new Float_t[kNstrips+1];
    leA[0] = -kDxyz[0];
    leA[1] = -kpitch*(0.5*kNstrips-1);
    leC[0] =  kDxyz[0];
    leC[1] =  kpitch*(0.5*kNstrips-1);
    for(i=1;i<kNstrips;i++){
	leA[i+1] = leA[i] + kpitch;
	leC[i+1] = leC[i] - kpitch;
    } // end for i
    leA[kNstrips] =  kDxyz[0];
    leC[kNstrips] = -kDxyz[0];
//    cout << "AliITSgeomSSD275and75 default creator called: start" << endl;
    AliITSgeomSSD::AliITSgeomSSD(kDxyz,kangleA,kangleC,
				 kNstrips+1,leA,kNstrips+1,leC);
    delete leA;
    delete leC;
//    cout << "AliITSgeomSSD275and75 default creator called: end" << endl;
}
//________________________________________________________________________
ostream &operator<<(ostream &os,AliITSgeomSSD75and275 &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSSD75and275 &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}
//======================================================================
