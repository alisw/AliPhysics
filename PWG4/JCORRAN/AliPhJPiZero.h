// $Id: AliPhJPiZero.h,v 1.4 2008/05/08 15:19:52 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliPhJPiZero.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.4 $
  \date $Date: 2008/05/08 15:19:52 $
*/
////////////////////////////////////////////////////

#ifndef ALIPHJPIZERO_H
#define ALIPHJPIZERO_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include  "AliPhJBaseTrack.h"
#include  "AliPhJPhoton.h"

//class AliPhJPhoton;
//class TObject;

class AliPhJPiZero : public AliPhJBaseTrack {

public:
  AliPhJPiZero();//{PhotSum.SetPxPyPzE(0,0,0,0);}   //constructor
  AliPhJPiZero(const AliPhJPiZero& a);//{PhotSum.SetPxPyPzE(0,0,0,0);}   //constructor
  ~AliPhJPiZero(){;}		      //destructor

  bool SetMass(AliPhJPhoton* g1, AliPhJPhoton* g2);

  float GetInvMass()  const {return fPizM;}
  float GetAsymm()    const {return fAsymm;}
  
  int   GetMassBin()  const {return fMassBin;}
  void  SetMassBin(int im)  {fMassBin=im;}
  void  ResetToZero(){fPizM=0; fBasePt=0; fBasePhi=0; fBaseTheta=0; fAsymm=0; fMassBin=-1;}
  
  TVector3 GetP() const {return fPi0P;}
  //TLorentzVector &GetSum() {return *PhotSum;}  //I think should be like that
  TLorentzVector GetSum() const {return fPhotSum;} 

  double operator- (const AliPhJPiZero &pi0);

  AliPhJPiZero& operator= (const AliPhJPiZero& piz);

protected:
  TVector3 fV1, fV2; //vector 1 and 2
  TVector3 fPi0P;     // pi0 momentum
  TLorentzVector fPhotSum; //sum of lorentz vectors of two photons
  TLorentzVector fPhoton1, fPhoton2; //decay photons
  float fAsymm, fPizM;//assymtery and inv. mass
  int fMassBin;    //mass bin
  
  ClassDef(AliPhJPiZero,1)
};

#endif

