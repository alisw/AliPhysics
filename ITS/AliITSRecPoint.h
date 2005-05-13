#ifndef ALIITSRECPOINT_H
#define ALIITSRECPOINT_H 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
*/

////////////////////////////////////////////////////
//  Reconstructed space point class for set:ITS   //
////////////////////////////////////////////////////

#include <TObject.h>
#include <Riostream.h>


class AliITSRecPoint : public TObject {
 public:
    AliITSRecPoint();
    virtual ~AliITSRecPoint() {}; // distructor
    Bool_t IsSortable() const {return kTRUE;} // allows for sorting
    Int_t   GetLabel(Int_t i) const {return fTracks[i];} // get track label
    Int_t  *GetTracks(){return fTracks;}// Returns pointer to track array
    Int_t   GetNTracks(){return 3;} // returns track array size
    Float_t GetX() const {return fX;} // gets fX
    Float_t GetZ() const {return fZ;} // gets fZ
    Float_t GetQ() const {return fQ;} // gets fQ
    Float_t GetdEdX() const {return fdEdX;} // gets fdEdX
    Float_t GetSigmaX2() const {return fSigmaX2;} // gets fSigmaX2
    Float_t GetSigmaZ2() const {return fSigmaZ2;} // gets fSigmaZ2
    void SetLabel(Int_t i, Int_t lab){fTracks[i]=lab;} // sets track label
    void SetX(Float_t x){fX=x;} // sets fX
    void SetZ(Float_t z){fZ=z;} // sets fZ
    void SetQ(Float_t q){fQ=q;} // sets fQ
    void SetdEdX(Float_t dedx){fdEdX=dedx;} // sets fdEdX
    void SetSigmaX2(Float_t sx2){fSigmaX2=sx2;} // sets fSigmaX2
    void SetSigmaZ2(Float_t sz2){fSigmaZ2=sz2;} // sets fSigmaZ2
    void  Use() { //if fQ<0 cluster is already associated with a track
	fQ=-fQ;}
    Int_t IsUsed() const {return (fQ<0) ? 1 : 0;} // checks Use condision
    Int_t Compare(const TObject *) const {return 0;} //to be defined
    // Prints out the content of this class in ASCII format.
    void Print(ostream *os); 
    // Reads in the content of this class in the format of Print
    void Read(istream *is);
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual Int_t Read(const char *name) {return TObject::Read(name);}

 public:
    Int_t     fTracks[3]; //labels of overlapped tracks
    Float_t   fX ;        //X of cluster
    Float_t   fZ ;        //Z of cluster
    Float_t   fQ ;        //Q of cluster (in ADC counts)
    Float_t   fdEdX;      //dE/dX inside this cluster
    Float_t   fSigmaX2;   //Sigma X square of cluster
    Float_t   fSigmaZ2;   //Sigma Z square of cluster

    ClassDef(AliITSRecPoint,1)  // AliITSRecPoint class
};
// Input and output function for standard C++ input/output.
ostream& operator<<(ostream &os,AliITSRecPoint &source);
istream& operator>>(istream &is,AliITSRecPoint &source);
#endif
