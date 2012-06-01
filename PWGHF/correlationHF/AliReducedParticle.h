
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

// class to get the reduced hadron candidate
class AliReducedParticle : public AliVParticle
//class AliReducedParticle : public AliAODTrack
{
public:
    AliReducedParticle(Double_t eta, Double_t phi, Double_t pt, Int_t McLabel, Int_t trackid);

    ~AliReducedParticle();
    
    // kinematics
	
    virtual Double_t Pt()         const { return fpT; }
	virtual Double_t Phi()        const { return fPhi; }
	virtual Double_t Eta()        const { return fEta; }
	virtual Int_t   GetLabel()    const { return fMcLabel; }
	virtual Int_t  GetID()		const{return fid;}
	
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
    
    virtual Short_t Charge()      const { AliFatal("Not implemented"); return 0; }
   
    // PID
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }      
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
	
	
	
    
private:
    Double_t fEta;      // eta
    Double_t fPhi;      // phi
    Double_t fpT;       // pT
	Int_t fMcLabel; //mclabel
	Int_t fid;

    
    ClassDef(AliReducedParticle, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};
#endif
