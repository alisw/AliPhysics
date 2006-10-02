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

/* $Id$ */

////////////////////////////////////////////////////////////////////////
// This class is for the Silicon Strip Detector, SSD, specific geometry.
// It is being replaced by AliITSsegmentationSSD class. This file also
// constains classes derived from AliITSgeomSSD which do nothing but
// initilize this one with predefined values.
////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <stdlib.h>
#include <TShape.h>
#include <TBRIK.h>

#include "AliITSgeomSSD.h"

ClassImp(AliITSgeomSSD)


AliITSgeomSSD::AliITSgeomSSD():
TObject(),
fName(),
fTitle(),
fMat(),
fDx(0.0),
fDy(0.0),
fDz(0.0),
fNp(0),
fNn(0),
fLowEdgeP(0),
fLowEdgeN(0),
fAngleP(0.0),
fAngleN(0.0){
// Default constructor
    fNp       = 0;
    fNn       = 0;
    fLowEdgeP = 0;
    fLowEdgeN = 0;
    fAngleP   = 0.0;
    fAngleN   = 0.0;
}
//----------------------------------------------------------------------
AliITSgeomSSD::AliITSgeomSSD(const Float_t *box,Float_t ap,Float_t an,
			     Int_t np,Float_t *p,Int_t nn,Float_t *n):
TObject(),
fName(),
fTitle(),
fMat(),
fDx(0.0),
fDy(0.0),
fDz(0.0),
fNp(0),
fNn(0),
fLowEdgeP(0),
fLowEdgeN(0),
fAngleP(0.0),
fAngleN(0.0){
////////////////////////////////////////////////////////////////////////
//    Standard Constructor. *box={dx,dy,dz}, ap=anode angle, an=cathode angle,
// nn= number of cathodes+1,*n= array of cathode low edges+highest edge,
// np= number of anodes+1, *p= array of anode low edges+lighest edge.
///////////////////////////////////////////////////////////////////////
    fNp       = 0;
    fNn       = 0;
    fLowEdgeP = 0;
    fLowEdgeN = 0;
    fAngleP   = 0.0;
    fAngleN   = 0.0;
    ResetSSD(box,ap,an,np,p,nn,n);
}
//----------------------------------------------------------------------
void AliITSgeomSSD::ResetSSD(const Float_t *box,Float_t ap,Float_t an,
			     Int_t np,Float_t *p,Int_t nn,Float_t *n){
////////////////////////////////////////////////////////////////////////
//    Standard Filler. *box={dx,dy,dz}, ap=anode angle, an=cathode angle,
// nn= number of cathodes+1,*n= array of cathode low edges+highest edge,
// np= number of anodes+1, *p= array of anode low edges+lighest edge.
///////////////////////////////////////////////////////////////////////
    Int_t i;

    fName = "ActiveSSD";
    fTitle = "Active volume of SSD";
    fMat = "SSD Si Det";
    fDx  = box[0];
    fDy  = box[1];
    fDz  = box[2];
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

    if(fLowEdgeP) delete [] fLowEdgeP; fLowEdgeP = 0;
    if(fLowEdgeN) delete [] fLowEdgeN; fLowEdgeN = 0;
    fNp = 0;
    fNn = 0;
    fAngleP = 0.0;
    fAngleN = 0.0;
}
//______________________________________________________________________
AliITSgeomSSD::AliITSgeomSSD(const AliITSgeomSSD &source) : TObject(source),
fName(source.fName),fTitle(source.fTitle),fMat(source.fMat),fDx(source.fDx),fDy(source.fDy),fDz(source.fDz),fNp(source.fNp),fNn(source.fNn),fLowEdgeP(0),fLowEdgeN(0),fAngleP(source.fAngleP),fAngleN(source.fAngleN){
////////////////////////////////////////////////////////////////////////
//    copy  constructor
////////////////////////////////////////////////////////////////////////
    Int_t i;

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

  this->~AliITSgeomSSD();
  new(this) AliITSgeomSSD(source);
  return *this;

}
//______________________________________________________________________
void AliITSgeomSSD::Local2Det(Float_t x,Float_t z,Int_t &a,Int_t &c){
    // Given a GEANT detector local coordinate, cm, this function returns
    // the detector specific P and N side strip numbers.
    // Inputs are:
    // Float_t x   Geant detector local x coordinate in cm
    // Float_t z   Geant detector local z coordinate in cm
    // outputs are:
    // Int_t &a    Detector anode strip number (P side)
    // Int_t &c    Detector cathode strip number (N side)
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
    // use AliITSsegmentationSSD.

    x=a;
    z=c;
    Error("Det2Locat","Use AliITSsegmentationSSD");
    return;
}
//______________________________________________________________________
void AliITSgeomSSD::Print(ostream *os) const {
////////////////////////////////////////////////////////////////////////
// Standard output format for this class.
////////////////////////////////////////////////////////////////////////
    Int_t i;
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC || defined __xlC__
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

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
    Int_t i;
    char shp[20];

    *is >> shp;
    *is >> fDx >> fDy >> fDz;
    *is >> fNp >> fNn;
    *is >> fAngleP >> fAngleN;
    if(fLowEdgeP !=0) delete fLowEdgeP;
    if(fLowEdgeN !=0) delete fLowEdgeN;
    fLowEdgeP = new Float_t[fNp];
    fLowEdgeN = new Float_t[fNn];
    for(i=0;i<fNp;i++) *is >> fLowEdgeP[i];
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
    AliITSgeomSSD::ResetSSD(kDxyz,kangle,-kangle,
				 kNstrips+1,leA,kNstrips+1,leC);
    delete leA;
    delete leC;
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
AliITSgeomSSD& AliITSgeomSSD175::operator=(const AliITSgeomSSD &source) {
////////////////////////////////////////////////////////////////////////
//    assignment operator
////////////////////////////////////////////////////////////////////////
  

    if(this == &source) return *this;
    Error("AliITSgeomSSD175","Not allowed to make a = with "
          "AliITSgeomSSD175 Using default creater instead");

    return *this;
}
//======================================================================

ClassImp(AliITSgeomSSD275and75)

AliITSgeomSSD275and75::AliITSgeomSSD275and75() : AliITSgeomSSD(){
////////////////////////////////////////////////////////////////////////
//    default constructor
////////////////////////////////////////////////////////////////////////
}
//----------------------------------------------------------------------
AliITSgeomSSD275and75::AliITSgeomSSD275and75(Int_t npar,Float_t *par) : 
                                                            AliITSgeomSSD(){
    // Default constructor for AliITSgeomSSD with strip angles of
    // 275 miliradians and 75 miliradians. This constructor initlizes
    // AliITSgeomSSD with the correct values. This is the miror image
    // of the AliITSgeomSSD75and275 class.
    const Float_t kDxyz[] ={3.6500,0.0150,2.000};//cm. (Geant 3.12 units)
    // Size of sensitive detector area x,y(thickness),z
    const Float_t kangleA  = 0.0275; // angle in rad. of anode and cathodes
    const Float_t kangleC  = 0.0075; // angle in rad. of anode and cathodes
    const Float_t kpitch   = 0.0095;// cm anode separation.
    const Int_t   kNstrips = 768; // number of anode or cathode strips.
    Float_t *leA,*leC; // array of low edges anode and cathorde.
    Int_t i;

    if(npar<3){
	Error("AliITSgeomSSD275and75",
	      "npar=%d<3. array par must be [3] or larger.",npar);
	return;
    } // end if
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
    AliITSgeomSSD::ResetSSD(par,kangleA,kangleC,
				 kNstrips+1,leA,kNstrips+1,leC);
    delete [] leA;
    delete [] leC;
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
AliITSgeomSSD& AliITSgeomSSD275and75::operator=(const AliITSgeomSSD &source) {
////////////////////////////////////////////////////////////////////////
//    assignment operator
////////////////////////////////////////////////////////////////////////
  

    if(this == &source) return *this;
    Error("AliITSgeomSSD275and75","Not allowed to make a = with "
          "AliITSgeomSSD275and75 Using default creater instead");

    return *this;
}

//======================================================================

ClassImp(AliITSgeomSSD75and275)

AliITSgeomSSD75and275::AliITSgeomSSD75and275() : AliITSgeomSSD(){
////////////////////////////////////////////////////////////////////////
//    default constructor
////////////////////////////////////////////////////////////////////////
}
AliITSgeomSSD75and275::AliITSgeomSSD75and275(Int_t npar,Float_t *par) : 
                                                            AliITSgeomSSD(){
    // Default constructor for AliITSgeomSSD with strip angles of
    // 75 miliradians and 275 miliradians. This constructor initlizes
    // AliITSgeomSSD with the correct values. This is the miror image
    // of the AliITSgeomSSD275and75 class.
    const Float_t kDxyz[] ={3.6500,0.0150,2.000};//cm. (Geant 3.12 units)
    // Size of sensitive detector area x,y(thickness),z
    const Float_t kangleA  = 0.0075; // angle in rad. of anode and cathodes
    const Float_t kangleC  = 0.0275; // angle in rad. of anode and cathodes
    const Float_t kpitch   = 0.0095;// cm anode separation.
    const Int_t   kNstrips = 768; // number of anode or cathode strips.
    Float_t *leA,*leC; // array of low edges anode and cathorde.
    Int_t i;

    if(npar<3){
	Error("AliITSgeomSSD75and275",
	      "npar=%d<3. array par must be [3] or larger.",npar);
	return;
    } // end if
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
    AliITSgeomSSD::ResetSSD(par,kangleA,kangleC,
				 kNstrips+1,leA,kNstrips+1,leC);
    delete leA;
    delete leC;
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
AliITSgeomSSD& AliITSgeomSSD75and275::operator=(const AliITSgeomSSD &source) {
////////////////////////////////////////////////////////////////////////
//    assignment operator
////////////////////////////////////////////////////////////////////////
  

    if(this == &source) return *this;
    Error("AliITSgeomSSD75and275","Not allowed to make a = with "
          "AliITSgeomSSD75and275 Using default creater instead");

    return *this;
}

//======================================================================
