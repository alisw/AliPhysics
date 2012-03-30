// $Id: AliJBaseEventHeader.h,v 1.5 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliJBaseEventHeader.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.5 $
  \date $Date: 2008/05/08 13:44:45 $
  */
////////////////////////////////////////////////////

#ifndef ALIJBASEEVENTHEADER_H
#define ALIJBASEEVENTHEADER_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif
#include <TNamed.h>

class AliJBaseEventHeader : public TNamed {

 public:
  AliJBaseEventHeader();                                    // default constructor
  AliJBaseEventHeader(int eventid, float cent, float vtxz); // constructor
  AliJBaseEventHeader(const AliJBaseEventHeader& a);                    // copy constructor
  virtual ~AliJBaseEventHeader(){;}                         // destructor

  //getter
  int    GetEventID() const {return fEventID;} 
  float  GetCentrality() const {return fCentrality;}
  float  GetXVertex() const {return fVtxX;}
  float  GetYVertex() const {return fVtxY;}
  float  GetZVertex() const {return fVtxZ;}
  float  GetZVertexErr() const {return fVtxZErr;}

  //setter
  void SetEventID(int evid) {fEventID=evid;}
  void SetCentrality(float  cent) {fCentrality=cent;}
  void SetXVertex(float vt) {fVtxX=vt;}
  void SetYVertex(float vt) {fVtxY=vt;}
  void SetZVertex(float vt) {fVtxZ=vt;}
  void SetZVertexErr(float vt) {fVtxZErr=vt;}
  void SetVertex(float x, float y, float z, float err){ fVtxX=x;fVtxY=y;fVtxZ=z;fVtxZErr=err; }

  AliJBaseEventHeader& operator=(const AliJBaseEventHeader& header);

 private:

  Int_t   fEventID;         //event id
  Double32_t fCentrality;   //centrality
  Double32_t fVtxX;         //vertex X
  Double32_t fVtxY;         //vertex Y
  Double32_t fVtxZ;         //vertex Z
  Double32_t fVtxZErr;      //vertex error

  ClassDef(AliJBaseEventHeader,1)

};

#endif
