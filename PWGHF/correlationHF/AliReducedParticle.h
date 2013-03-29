#ifndef AliReducedParticle_H
#define AliReducedParticle_H


/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//             Base class for DStar - Hadron Correlations Analysis
//
//-----------------------------------------------------------------------
//          
//
//						   Author S.Bjelogrlic
//                         Utrecht University 
//                      sandro.bjelogrlic@cern.ch
//
//-----------------------------------------------------------------------


#include "AliLog.h"
#include "AliVParticle.h"

// class to get the reduced hadron candidate
class AliReducedParticle : public AliVParticle
{
public:
    AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Int_t mcLabel, Int_t trackid, Double_t impPar, Bool_t checkSoftPi);
    AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Int_t mcLabel, Int_t trackid, Double_t impPar, Bool_t checkSoftPi, Short_t charge);
    AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Int_t mcLabel, Int_t trackid, Double_t impPar, Bool_t checkSoftPi, Short_t charge,Double_t weight);
    AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Int_t mcLabel);
    AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, int charge, int orginmother); 
    AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Double_t invmass, int ptbin, int orginmother=0);

    ~AliReducedParticle();
    
    // kinematics
	
    virtual Double_t Pt()         const { return fpT; }
    virtual Double_t Phi()        const { return fPhi; }
    virtual Double_t Eta()        const { return fEta; }
    virtual Int_t   GetLabel()    const { return fMcLabel; }
    virtual Int_t  GetID()		const{return fid;}
    virtual Double_t GetImpPar() const{return fImpPar;}
    virtual Bool_t CheckSoftPi() const{return fCheckSoftPi;}
    virtual Double_t GetInvMass() const {return fInvMass;}
    virtual int GetPtBin() const  {return fPtBin;}
    virtual int GetOriginMother() const {return fOriginMother;}
    void SetWeight(Double_t weight){fWeight=weight;}
    Double_t GetWeight(){return fWeight;}
    
	// kinematics
    virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
	
    virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
	
    virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }
	
	
    virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
    

    virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }
    
    virtual Short_t Charge()      const { return fCharge;}
   
    // PID
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }      
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
	
	
	
    
private:
    Double_t fEta;       // eta
    Double_t fPhi;       // phi
    Double_t fpT;        // pT
    Int_t fMcLabel;      //mclabel
    Int_t fid;           // track ID
    Double_t fImpPar;    // impact parameter
    Bool_t fCheckSoftPi; // check if the track is compatible with a softpion from D*
    Short_t fCharge;     // charge of the associated track
    Double_t fInvMass;   // Invariant mass of Dmeson
    Double_t fWeight;   // track weight (e.g. 1/efficiency)
    int fPtBin;          // Ptbin of Dmesons
    int fOriginMother;   //  Holds the origin motherquark (process)

    
    ClassDef(AliReducedParticle, 4); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};
#endif
