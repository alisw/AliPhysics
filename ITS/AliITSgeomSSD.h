#ifndef ALIITSGEOMSSD_H
#define ALIITSGEOMSSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <TBRIK.h>

class TShape;
class ostream;
class istream;

class AliITSgeomSSD : public TObject {

 public:
    AliITSgeomSSD(); // default constructor
    AliITSgeomSSD(const Float_t *box,Float_t ap,Float_t an,
		  Int_t np,Float_t *p,Int_t nn,Float_t *n); // Constructor
    virtual ~AliITSgeomSSD(); // Destructor
    AliITSgeomSSD(const AliITSgeomSSD &source);// copy constructor
    virtual AliITSgeomSSD& operator=(const AliITSgeomSSD &source); // = opt.
    void ResetSSD(const Float_t *box,Float_t ap,Float_t an,
		  Int_t np,Float_t *p,Int_t nn,Float_t *n); // Filler
    virtual TShape *GetShape() const {return fShapeSSD;}// get shape
    virtual Float_t GetDx() const {if(fShapeSSD!=0) return fShapeSSD->GetDx();
                    else return 0.0;}// get Dx
    virtual Float_t GetDy() const {if(fShapeSSD!=0) return fShapeSSD->GetDy();
                    else return 0.0;}// get Dy
    virtual Float_t GetDz() const {if(fShapeSSD!=0) return fShapeSSD->GetDz();
                    else return 0.0;}// get Dz
    virtual Int_t GetNAnodes() const {return fNp-1;}//the number of Anodes "P"
    virtual Int_t GetNCathodess() const {return fNn-1;}//the number of Cathodes "N"
    virtual Float_t GetAnodePitch(Int_t i=0) const { //anode pitch for anode i
	if(i>=0&&i<fNp) return fLowEdgeP[i+1]-fLowEdgeP[i];else return 0.0;}
    virtual Float_t GetCathodePitch(Int_t i=0) const { // cathode pitch for cathode i
	if(i>0&&i<fNn) return fLowEdgeN[1]-fLowEdgeN[0];else return 0.0;}
    virtual Float_t GetAnodeAngle() const {return fAngleP;}//anode strip angle.
    virtual Float_t GetCathodeAngle() const {return fAngleN;}//cathode strip angle.
    virtual void SetShape(char *name,char *title,char *mat,
                          Float_t dx,Float_t dy,Float_t dz){
	// defines TBRIK with given paramters
        fShapeSSD = new TBRIK(name,title,mat,dx,dy,dz);};
    virtual void SetNAnodes(Int_t n) {// sets the number of Anodes "P" and
	// allocates array of low edges.
	fNp=n+1;delete fLowEdgeP;fLowEdgeP = new Float_t[fNp];}
    virtual void SetNCathotess(Int_t n) {// sets the number of Anodes "N" and 
	// allocates array of low edges.
	fNn=n+1;delete fLowEdgeN;fLowEdgeN =new  Float_t[fNn];}
    virtual void SetAnodeLowEdges(Float_t *p){// sets Anode low edges +1.
	for(Int_t i=0;i<fNp;i++)fLowEdgeP[i]=p[i];}
    virtual void SetCathodeLowEdges(Float_t *p){// sets Cathodes low edges +1.
	for(Int_t i=0;i<fNn;i++)fLowEdgeN[i]=p[i];}
    virtual void SetAnodeAngle(Float_t a){fAngleP=a;} //sets anode angle.
    virtual void SetCathodeAngle(Float_t a){fAngleN=a;}//sets cathode  angle.

    virtual void Local2Det(Float_t x,Float_t z,Int_t &a,Int_t &c);
    virtual void Det2Local(Int_t a,Int_t c,Float_t &x,Float_t &z);

    virtual void Print(ostream *os) const;  // Output streamer to standard out.
    virtual void Read(istream *is);   // Input streamer to standard in.

 protected:
    // -+-> x
    //  |
    //  V
    //  z
    TBRIK *fShapeSSD; // comment
    Int_t   fNp;      // Number of Anode strips.
    Int_t   fNn;      // Number of Cathode strips.
    Float_t *fLowEdgeP;  //[fNp] Anode side strip pitch angle==0.
    Float_t *fLowEdgeN;  //[fNn] Cathode side strip pich angle==0.
    Float_t fAngleP;  // Anode side strip angle (rad).
    Float_t fAngleN;  // Cathode side strip angle (rad).
    // or what other or different information that is needed.
    
    ClassDef(AliITSgeomSSD,1) // ITS SSD detector geometry class

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSSD &source);
istream &operator>>(istream &os,AliITSgeomSSD &source);
#endif
//======================================================================
#ifndef ALIITSGEOMSSD175_H
#define ALIITSGEOMSSD175_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
*/

//#include "AliITSgeomSSD.h"

class ostream;
class istream;

class TShape;

class AliITSgeomSSD175 : public AliITSgeomSSD {

 public:
    AliITSgeomSSD175();

    // This clas now has version 0 so that it will not be written to a root
    // file. This is good since there are no longer any data members to this
    // class. It is only designed to make it easer to define this standard
    // SDD detector geometry.
    ClassDef(AliITSgeomSSD175,0) // ITS SSD detector with stips at +- 0.0175 rad.

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSSD175 &source);
istream &operator>>(istream &os,AliITSgeomSSD175 &source);
#endif
//======================================================================
#ifndef ALIITSGEOMSSD27575_H
#define ALIITSGEOMSSD27575_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
*/

//#include "AliITSgeomSSD.h"

class ostream;
class istream;

class TShape;

class AliITSgeomSSD275and75 : public AliITSgeomSSD {

 public:
    AliITSgeomSSD275and75();

    // This clas now has version 0 so that it will not be written to a root
    // file. This is good since there are no longer any data members to this
    // class. It is only designed to make it easer to define this standard
    // SDD detector geometry.
    ClassDef(AliITSgeomSSD275and75,0) // ITS SSD detector with 0.0275 and 0.0075 rad strip angles.

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSSD275and75 &source);
istream &operator>>(istream &os,AliITSgeomSSD275and75 &source);
#endif
//======================================================================
#ifndef ALIITSGEOMSSD75275_H
#define ALIITSGEOMSSD75275_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
*/

//#include "AliITSgeomSSD.h"

class ostream;
class istream;
class TShape;

class AliITSgeomSSD75and275 : public AliITSgeomSSD {

 public:
    AliITSgeomSSD75and275();

    // This clas now has version 0 so that it will not be written to a root
    // file. This is good since there are no longer any data members to this
    // class. It is only designed to make it easer to define this standard
    // SSD detector geometry.
    ClassDef(AliITSgeomSSD75and275,0) // ITS SSD detector geometry class for 0.0075 and 0.0275 rad angled strips.

};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomSSD75and275 &source);
istream &operator>>(istream &os,AliITSgeomSSD75and275 &source);
#endif
