#ifndef ALIITSGEOMSDD_H
#define ALIITSGEOMSDD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
#include <TBRIK.h>

class TShape;

class AliITSgeomSDD: public TObject {
 public:
    AliITSgeomSDD();
    AliITSgeomSDD(const Float_t *box,Float_t per,Float_t vel,
		  Float_t axL,Float_t axR,
		  Int_t nA0,Float_t *le0,Int_t nA1,Float_t *le1);
    AliITSgeomSDD(AliITSgeomSDD &source);
    AliITSgeomSDD& operator=(AliITSgeomSDD &source);
    virtual ~AliITSgeomSDD();
    void ResetSDD(const Float_t *box,Float_t per,Float_t vel,
		  Float_t axL,Float_t axR,
		  Int_t nA0,Float_t *le0,Int_t nA1,Float_t *le1);
    virtual TShape *GetShape() const {return fShapeSDD;}
    virtual Float_t GetDx() const { // Get TBRIK Dx
	if(fShapeSDD!=0) return fShapeSDD->GetDx();
	else return 0.0;}
    virtual Float_t GetDy() const {// Get TBRIK Dy
	if(fShapeSDD!=0) return fShapeSDD->GetDy();
	else return 0.0;}
    virtual Float_t GetDz() const {// Get TBRIK Dz
	if(fShapeSDD!=0) return fShapeSDD->GetDz();
	else return 0.0;}
    virtual Float_t GetAnodeX(Int_t a,Int_t s) const { // returns X position of anode
	if(s==0) return fAnodeXL; else return fAnodeXR;}
    virtual Float_t GetAnodeZ(Int_t a,Int_t s)const { // returns X position of anode
	if(s==0) return 0.5*(fAnodeLowEdgeL[a]+fAnodeLowEdgeL[a+1]);
	else return 0.5*(fAnodeLowEdgeR[a]+fAnodeLowEdgeR[a+1]);}
    virtual void SetNAnodesL(Int_t s)
	{fNAnodesL = s;} // sets the number of anodes on side 0.
    virtual void SetNAnodesR(Int_t s)
	{fNAnodesR = s;} // sets the anodes spacing for side 1.
    virtual void SetSamplingPeriod(Float_t s)
	{fPeriod = s;} // sets the clock sampling period s.
    virtual void SetDriftVelocity(Float_t s)
	{fDvelocity = s;} // sets the SDD Drift velocity cm/s.
    virtual void SetShape(char *name,char *title,char *mat,
			  Float_t dx,Float_t dy,Float_t dz)
	                 {fShapeSDD = new TBRIK(name,title,mat,dx,dy,dz);}
    virtual void Local2Det(Float_t xl,Float_t zl,Int_t &a,Int_t &t,Int_t &s);
    virtual void Det2Local(Int_t a,Int_t t,Int_t s,Float_t &xl,Float_t &zl);
    virtual void Print(ostream *os) const; // Output streamer to standard out.
    virtual void Read(istream *is);   // Input streamer to standard in.
    // or what other or different information that is needed.

 protected:
    // (L)  -+-> x  (R)
    //       |
    //       V z
    Float_t fPeriod; // ADC sampiling period
    Float_t fDvelocity; // Drift velocity
    Int_t   fNAnodesL;  // number of Anodes on size 0
    Int_t   fNAnodesR;  // number of Anodes on size 1
    Float_t fAnodeXL;   // Anode location in x Left side
    Float_t fAnodeXR;   // Anode location in x Right side
    Float_t *fAnodeLowEdgeL; //[fNAnodesL] Anode spacing left edge
    Float_t *fAnodeLowEdgeR; //[fNAnodesR] Anode spacing right edge
    TBRIK *fShapeSDD;     // shape of sensitive volume

    ClassDef(AliITSgeomSDD,1) // ITS SDD detector geometry class

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSDD &source);
istream &operator>>(istream &os,AliITSgeomSDD &source);
#endif
//======================================================================
#ifndef ALIITSGEOMSDD256_H
#define ALIITSGEOMSDD256_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
*/

//#include "AliITSgeomSDD.h"


class AliITSgeomSDD256 : public AliITSgeomSDD {

 public:
    AliITSgeomSDD256();

    // This clas now has version 0 so that it will not be written to a root
    // file. This is good since there are no longer any data members to this
    // class. It is only designed to make it easer to define this standard
    // SDD detector geometry.
    ClassDef(AliITSgeomSDD256,0) // ITS SDD detector geometry class for 256 anodes per side

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSDD256 &source);
istream &operator>>(istream &os,AliITSgeomSDD256 &source);
#endif
//======================================================================
#ifndef ALIITSGEOMSDD300_H
#define ALIITSGEOMSDD300_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
*/

//#include "AliITSgeomSDD.h"


class AliITSgeomSDD300 : public AliITSgeomSDD {

 public:
    AliITSgeomSDD300();

    // This clas now has version 0 so that it will not be written to a root
    // file. This is good since there are no longer any data members to this
    // class. It is only designed to make it easer to define this standard
    // SDD detector geometry.
    ClassDef(AliITSgeomSDD300,0) // ITS SDD detector geometry class for 300 anodes per side

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSDD300 &source);
istream &operator>>(istream &os,AliITSgeomSDD300 &source);
#endif
