/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

/* $Id: AliRsnPID.h,v 1.5 2007/02/21 14:33:25 pulvir Exp $ */

//-------------------------------------------------------------------------
//                      Class AliRsnPID
//  Simple collection of reconstructed tracks, selected from an ESD event
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNPID_H
#define ALIRSNPID_H

class AliRsnDaughter;
class AliRsnEvent;

class AliRsnPID : public TObject
{
public:
	
    AliRsnPID();
    AliRsnPID(const AliRsnPID& copy);
    AliRsnPID& operator=(const AliRsnPID& copy);
    virtual ~AliRsnPID() {}
    
    // types enum
    enum EType { 
        kElectron = 0, 
        kMuon, 
        kPion, 
        kKaon, 
        kProton, 
        kUnknown, 
        kSpecies = 5 
    };
    
    // identification method enum
    enum EMethod {
        kNone = 0,
        kRealistic,
        kPerfect
    };
    
    // conversions from PDG code to local type 
    static EType         InternalType(Int_t pdgCode);
    
    // retrieve particle informations from internal type
    static Int_t         PDGCode(EType pid);
    static const char *  ParticleName(EType pid);
    static Double_t      ParticleMass(EType pid);
    
    // identification routines
    Bool_t  Identify(AliRsnDaughter *d);
    Bool_t  Identify(AliRsnEvent *e);
    
    // setters
    void SetMethod(EMethod method) {fMethod = method;}
	void SetPriorProbability(EType type, Double_t p);
	void SetMinProb(Double_t p) {fMinProb = p;}
	void SetMaxPt(Double_t p) {fMaxPt = p;}
	
	// getters
	EMethod  GetMethod() {return fMethod;}
    
private:

    // identification routines
    Bool_t  IdentifyPerfect(AliRsnDaughter *d);
    Bool_t  IdentifyRealistic(AliRsnDaughter *d);
    Bool_t  Unidentify(AliRsnDaughter *d);

    EMethod    fMethod;          // PID method
    
    Double_t   fPrior[kSpecies]; // prior probabilities 
    Double_t   fMaxPt;           // realistic PID is done for Pt < this parameter
    Double_t   fMinProb;         // threshold on acceptable largest probability
    
    static const char* fgkParticleName[kSpecies+1];
    static const Int_t fgkParticlePDG[kSpecies+1];

    ClassDef(AliRsnPID,1);
};

#endif
