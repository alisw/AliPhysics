#ifndef ALIGENREADERECALJETS_H
#define ALIGENREADERECALJETS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenReader.h"


class AliGenReaderEcalJets : public AliGenReader
{
 public:
    AliGenReaderEcalJets();
    
    AliGenReaderEcalJets(const AliGenReaderEcalJets &reader):AliGenReader(reader)
	{reader.Copy(*this);}
    virtual ~AliGenReaderEcalJets(){;}
    // Initialise 
    virtual void Init();
    // Read
    virtual Int_t NextEvent();
    virtual TParticle*  NextParticle();
    AliGenReaderEcalJets & operator=(const AliGenReaderEcalJets & rhs);
 private:
    void Copy(AliGenReaderEcalJets&) const;
 protected:
    Int_t           fNcurrent;      // points to the next event
    Int_t           fNparticle;     // points to the next particle 
    Int_t           fNev;           // event number
    Float_t         fX[2];          // 
    Int_t           fXtyp[2];       // parton type
    Int_t           fNpart;         // number of particles  
    Float_t         fXpt[200];      // pt of particle
    Float_t         fXeta[200];     // eta of particle
    Float_t         fXphi[200];     // phi of particle
    Int_t           fXid[200];      // id of particle
    Int_t           fNjet;          // number of jets 
    Float_t         fJet[10];       // E_t of jet
    Float_t         fJeta[10];      // eta of jet
    Float_t         fJphi[10];      // phi of jet
    Int_t           fNsjet;         // number of clusters
    Float_t         fJset[10];      // E_t of cluster 
    Float_t         fJseta[10];     // eta of cluster
    Float_t         fJsphi[10];     // phi of cluster
    Int_t           fNpjet;         // ?
    Float_t         fJpet[10];      // ?
    Float_t         fJpeta[10];     // ?
    Float_t         fJpphi[10];     // ?

    TTree            *fTreeNtuple;    // pointer to the TTree
    //Declaration of leaves types
    ClassDef(AliGenReaderEcalJets,1) // Read particles from cwn-ntuple
};
#endif






