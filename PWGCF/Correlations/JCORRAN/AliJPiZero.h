/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

// $Id: AliJPiZero.h,v 1.4 2008/05/08 15:19:52 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliJPiZero.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.4 $
  \date $Date: 2008/05/08 15:19:52 $
 */
////////////////////////////////////////////////////

#ifndef ALIJPIZERO_H
#define ALIJPIZERO_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <TVector3.h>
#include <TLorentzVector.h>

class AliJBaseTrack;
class AliJPhoton;
// #include  "AliJBaseTrack.h"
// #include  "AliJPhoton.h"

//class AliJPhoton;
//class TObject;

class AliJPiZero : public AliJBaseTrack {

  public:
    AliJPiZero();
    virtual ~AliJPiZero(){;}          //destructor
    AliJPiZero(const AliJPiZero& a);//{PhotSum.SetPxPyPzE(0,0,0,0);}   //constructor
    AliJPiZero& operator=(const AliJPiZero& pizero);

    bool SetMass(AliJPhoton* g1, AliJPhoton* g2);

    float GetInvMass()  const {return M();}
    float GetAsymm()    const {return fAsymm;}

    int   GetMassBin()  const {return fMassBin;}
    void  SetMassBin(int im)  {fMassBin=im;}
    void  ResetToZero() {SetPxPyPzE(0, 0, 0, 0);fAsymm=0; fMassBin=-1; fIsGood=kFALSE; fConvPlaneDelta=-999; fPtBin=-1;}
    AliJPhoton *GetPhoton( Int_t i ) const { return fPhoton[i]; }
    void SetGood( Bool_t b ) { fIsGood = b;}
    Bool_t GetGood() const { return fIsGood; }
    void SetConvPlaneDelta( Double_t d ) {fConvPlaneDelta = d; }
    Double_t GetConvPlaneDelta() const { return fConvPlaneDelta; }
    void SetPtBin( Int_t b ) { fPtBin = b; }
    Int_t GetPtBin() const { return fPtBin; }
    
    Double_t Theta() { return fPhoton[0]->Angle( fPhoton[1]->Vect() ); }
    Double_t DPhi() { return TMath::Abs( fPhoton[0]->Phi() - fPhoton[1]->Phi() ); }
    Double_t DZ() { return TMath::Abs( fPhoton[0]->GetZ() - fPhoton[1]->GetZ() ); }

    TVector3 GetP() const {return Vect();}
    TLorentzVector GetSum() const {return *(TLorentzVector*)this;} 

    double DeltaConversion();
        
//     double operator- (const AliJPiZero &pi0);
//     AliJPiZero& operator= (const AliJPiZero& piz);

  protected:
    float fAsymm;  //assymtery and inv. mass
    int fMassBin; //mass bin
    Bool_t fIsGood; //!
    AliJPhoton *fPhoton[2]; //!
    Double_t fConvPlaneDelta; //! plane delta form conversion plane
    Int_t fPtBin; //! pt bin

    ClassDef(AliJPiZero,1)
};

#endif

