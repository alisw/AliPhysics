// $Id: AliPhJPhoton.h,v 1.5 2008/05/08 15:19:52 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliPhJPhoton.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.5 $
  \date $Date: 2008/05/08 15:19:52 $
*/
////////////////////////////////////////////////////

#ifndef ALIPHJPHOTON_H
#define ALIPHJPHOTON_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "TMath.h"
#include "JConst.h"
#include "AliPhJBaseTrack.h"

//class TObject;

class AliPhJPhoton : public AliPhJBaseTrack {

public:

  AliPhJPhoton();	    //default constructor
  AliPhJPhoton(const AliPhJPhoton& a); //copy constructor
  ~AliPhJPhoton(){;}		//destructor    
  
  float  GetChi2() const {return fChi2;}
  float  GetTof() const {return fTof;}                   
  float  GetX() const {return fX;}            
  float  GetY() const {return fY;}          
  float  GetZ() const {return fZ;}
  float  GetProbPhot() const {return fProbPhot;}

  void  SetChi2(float chi2) {fChi2=chi2;}
  void  SetTof(float tof) {fTof=tof;}
  void  SetX(float x) {fX=x;}
  void  SetY(float y) {fY=y;}
  void  SetZ(float z) {fZ=z;}
  void  SetProbPhot(float prob) {fProbPhot=prob;}

  AliPhJPhoton& operator=(const AliPhJPhoton& photon);
  
private:

  float  fChi2;      //chi2             
  float  fTof;       //time of flight               
  float  fX, fY, fZ; // x,y,z coordinates              
  float  fProbPhot;  //probability to be a photon

 ClassDef(AliPhJPhoton,1)
 
};

#endif

