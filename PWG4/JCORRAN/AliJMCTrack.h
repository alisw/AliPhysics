// $Id: AliJMCTrack.h,v 1.3 2008/05/08 15:19:52 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliJMCTrack.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.3 $
  \date $Date: 2008/05/08 15:19:52 $
*/
////////////////////////////////////////////////////

#ifndef ALIJMCTRACK_H
#define ALIJMCTRACK_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "AliPhJBaseTrack.h"

//class TObject;

class AliJMCTrack : public AliPhJBaseTrack {

public:

  AliJMCTrack();	    //default constructor
  AliJMCTrack(const AliJMCTrack& a);	//copy constructor

  ~AliJMCTrack(){;}		//destructor
  
  Int_t  GetPdgCode()    const {return fPdgCode;}
  Int_t  GetStatusCode()    const {return fStatus;}
  Int_t  GetFlag()       const {return fFlag;}
  Int_t  GetLabel()      const {return fLabel;}
  Int_t  GetMother  (Int_t i) const {return fMother[i];}
  Int_t  GetDaughter(Int_t i) const {return fDaughter[i];}
  Double_t  GetEta()     const {return fEta;}
  bool  IsPrimary()      const {return fPrimary;}
  bool  IsInPHOS()       const {return fPHOS;}
  bool  IsInEMCAL()      const {return fEMCAL;}
  bool  IsInTPC()        const {return fTPC;}
  Double_t  GetPtHard()     const {return fPtHard;}

  void SetPdgCode(Int_t icode) {fPdgCode=icode;}
  void SetStatusCode(Int_t icode) {fStatus=icode;}
  void SetFlag(Int_t i)        {fFlag=i;}
  void SetLabel(Int_t i)       {fLabel=i;}
  void SetMother  (int i, int code){ fMother[i] = code ; }
  void SetDaughter(int i, int code){ fDaughter[i] = code ; }
  void SetEta(Double_t i)      {fEta=i;}
  void SetProductionVertex(Double_t vx, Double_t vy, Double_t vz)
       {fVx=vx; fVy=vy; fVz=vz;}
  void SetIsPrimary(bool in) {fPrimary=in;}
  void SetIsInPHOS(bool in) {fPHOS=in;}
  void SetIsInEMCAL(bool in) {fEMCAL=in;}
  void SetIsInTPC(bool in) {fTPC=in;}
  void SetPtHard(Double_t i)      {fPtHard=i;}


  AliJMCTrack& operator=(const AliJMCTrack& trk);

	
private:

  Int_t            fPdgCode;              // PDG code of the particle
  Int_t            fStatus;               //status flag
  Int_t            fFlag;                 // Flag for indication of primary etc
  Int_t            fLabel;                // Label of the original MCParticle
  Int_t            fMother[2];               // Index of the mother particles
  Int_t            fDaughter[2];          // Indices of the daughter particles

  Double32_t       fVx;                   // [0.,0.,12] x of production vertex
  Double32_t       fVy;                   // [0.,0.,12] y of production vertex
  Double32_t       fVz;                   // [0.,0.,12] z of production vertex
  Double_t         fEta;                  //pseudorapidity
  bool fPrimary;                          //primary particle flag
  bool fPHOS;                             //hit in PHOS
  bool fEMCAL;                            //hit in EMCAL
  bool fTPC;                              //hit in TPC
  Double_t fPtHard;                       //hard particle


 ClassDef(AliJMCTrack,1)
};

#endif
