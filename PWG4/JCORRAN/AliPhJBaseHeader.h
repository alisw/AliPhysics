#ifndef ALIPHJBASEHEADER_H
#define ALIPHJBASEHEADER_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif


// $Id: AliPhJBaseHeader.h,v 1.5 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliPhJBaseHeader.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.5 $
  \date $Date: 2008/05/08 13:44:45 $
*/
////////////////////////////////////////////////////

class AliPhJBaseHeader : public TObject {
  
public:

  AliPhJBaseHeader();                                    // default constructor
  AliPhJBaseHeader(int eventid, short cent, float vtxz); // constructor
  AliPhJBaseHeader(const AliPhJBaseHeader& a);                    // copy constructor
  virtual ~AliPhJBaseHeader(){;}                         // destructor
  
  //getter
  int    GetEventID() const {return fEventID;} 
  short  GetCentrality() const {return fCentrality;}
  float  GetZVertex() const {return fVtxZ;}
  float  GetZVertexErr() const {return fVtxZErr;}
   
 //setter
  void SetEventID(int evid) {fEventID=evid;}
  void SetCentrality(short  cent) {fCentrality=cent;}
  void SetZVertex(float vt) {fVtxZ=vt;}
  void SetZVertexErr(float vt) {fVtxZErr=vt;}

  AliPhJBaseHeader& operator=(const AliPhJBaseHeader& header);

 private:
  
  int   fEventID; //event id
  short fCentrality; //centrality
  float fVtxZ; //vertex Z
  float fVtxZErr; //vertex error
 

  ClassDef(AliPhJBaseHeader,1)

};

#endif
