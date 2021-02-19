/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

// $Id: AliJBaseTrack.h,v 1.5 2008/05/08 15:19:52 djkim Exp $

///////////////////////////////////////////////////
/*
   \file AliJBaseTrack.h
   \brief
   \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
   \email: djkim@jyu.fi
   \version $Revision: 1.5 $
   \date $Date: 2008/05/08 15:19:52 $
   */
///////////////////////////////////////////////////

#ifndef ALIJBASETRACK_H
#define ALIJBASETRACK_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <iostream>
#include <TLorentzVector.h>
#include <TMath.h>
#include  "AliJConst.h"

using namespace std;

class AliJBaseTrack : public TLorentzVector {
    public:
        enum { kIsIsolated,kPrimary,kNFlag };
        AliJBaseTrack();
        AliJBaseTrack(float px,float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge); // constructor
        AliJBaseTrack(const AliJBaseTrack& a);
        AliJBaseTrack(const TLorentzVector & a);
        virtual ~AliJBaseTrack(){;}    //destructor

        double EtaAbs(){ return TMath::Abs(Eta()); }
        float   GetTwoPiPhi() const {return Phi()>-kJPi/3 ? Phi() : kJTwoPi+Phi();} 
        TLorentzVector GetLorentzVector(){ return TLorentzVector(Px(), Py(), Pz(), E());}

        Int_t         GetID()           const { return fID;}
        Int_t         GetLabel()        const { return fLabel; }
        Short_t       GetParticleType() const { return fParticleType;}
        ULong64_t     GetStatus()       const { return fStatus; }
        Short_t       GetCharge()       const { return fCharge; } 
        UInt_t        GetFlags()        const { return fFlags; }
        Bool_t        GetIsIsolated()   const { return IsTrue(kIsIsolated);}

        Int_t         GetTriggBin()     const { return fTriggID; }
        Int_t         GetAssocBin()     const { return fAssocID; }
        Double_t      GetTrackEff()     const { 
            /*if(fTracEff==-1) {  cout<<"AliJBaseTrack: Uninitilized track eff " <<endl;  exit(-1);
            } else return fTracEff;*/ return fTracEff; }
        Bool_t        IsInTriggerBin()  const { return fTriggID>=0; }
        Bool_t        IsInAssocBin()    const { return fAssocID>=0; }
        Double_t      GetWeight()       const { return fWeight;}             
        Int_t         GetMCIndex()      const { return fMCIndex;}

        void SetID      (const int id){fID=id;}
        void SetLabel   (const Int_t label ){ fLabel=label; }
        void SetParticleType(const Short_t ptype){ fParticleType=ptype; }
        void SetStatus  (const ULong64_t status){ fStatus=status; }
        void SetCharge  (const Char_t charge){ fCharge=charge; }
        void SetFlags   (const UInt_t bits ){ fFlags=bits; }        //MC, is primary flag
        void SetIsIsolated(Bool_t tf){ SetFlag( kIsIsolated, tf); }

        void SetTriggBin(const int id){fTriggID = id;}
        void SetAssocBin(const int id){fAssocID = id;}
        void SetTrackEff(const Double_t inEff){fTracEff = inEff;}

        void SetWeight(Double_t weight) { fWeight = weight;}
        void SetMCIndex(Int_t idx) {      fMCIndex = idx;}

        void SetPrimary(Bool_t b=kTRUE){ SetFlag(kPrimary,b);}
        Bool_t IsPrimary(){return IsTrue(kPrimary);}

        virtual void Print(Option_t *option="") const;

        // Handel BitsData
        Bool_t IsTrue(int i ) const { return TESTBIT(fFlags, i); }
        void SetFlag(int i, Bool_t t){ if(t){SETBIT(fFlags,i);}else{CLRBIT(fFlags, i);}}

        // Operators
        AliJBaseTrack& operator=(const AliJBaseTrack& trk);

    protected:
        Int_t         fID;            // Unique track ID
        Int_t         fLabel;         // Unique track label for MC-Data relation
        Short_t       fParticleType;  // ParticleType 
        Char_t        fCharge;        // track charge for real data
        ULong64_t     fStatus;        // reconstruction status flags or MC status 
        UInt_t        fFlags;         // store series of any boolen value.

        Int_t         fTriggID, fAssocID; //!   //id of trigger and assoc particle 
        Double_t      fTracEff;           //!   //track efficiency
        Int_t         fMCIndex;           //!   //index of corresp. MC track
        Double_t      fWeight;            //!   //particle weight

        ClassDef(AliJBaseTrack,2)
};

#endif

