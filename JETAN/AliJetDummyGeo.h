#ifndef ALIJETDUMMYGEO_H
#define ALIJETDUMMYGEO_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include <TObject.h>

class AliJetDummyGeo : public TObject 
{
 public: 
    AliJetDummyGeo(){;}
    virtual ~AliJetDummyGeo(){;}
    static AliJetDummyGeo* GetInstance() {return new AliJetDummyGeo();}
    static AliJetDummyGeo* GetInstance(char* /*name*/, char* /*title*/)
	{return new AliJetDummyGeo();}
    Int_t GetNCells(){ return 0;}
    Float_t GetArm1EtaMin() {return 0.;}
    Float_t GetArm1EtaMax() {return 0.;}
    Float_t GetArm1PhiMin() {return 0.;}
    Float_t GetArm1PhiMax() {return 0.;}
    void    EtaPhiFromIndex(Int_t /*id*/, Float_t& /*eta*/, Float_t& /*phi*/)
	{;}
    Int_t   TowerIndexFromEtaPhi2(Float_t /*eta*/, Float_t /*phi*/) {return 0;}
    void    GetTransformationForSM(){;}
    Float_t GetSampling() {return 0.;}
    ClassDef(AliJetDummyGeo,1)
};
 
#endif
