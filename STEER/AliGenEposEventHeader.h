#ifndef ALIGENEPOSEVENTHEADER_H_
#define ALIGENEPOSEVENTHEADER_H_
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// *
// * 
// * Header for EPOS generated event.
// *
// *      Author: Piotr Ostrowski
// *


#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
class TGenerator;


class AliGenEposEventHeader : public AliGenEventHeader, public AliCollisionGeometry
{
public:
	AliGenEposEventHeader(const char* name);
	AliGenEposEventHeader();
	virtual ~AliGenEposEventHeader() {}

Float_t GetBimevt() const { return fBimevt; }
Float_t GetPhievt() const { return fPhievt; }
Int_t GetKolevt()   const { return fKolevt; }
Int_t GetKoievt()   const { return fKoievt; }
Float_t GetPmxevt() const { return fPmxevt; }
Float_t GetEgyevt() const { return fEgyevt; }
Int_t GetNpjevt()   const { return fNpjevt; }
Int_t GetNtgevt()   const { return fNtgevt; }
Int_t GetNpnevt()   const { return fNpnevt; }
Int_t GetNppevt()   const { return fNppevt; }
Int_t GetNtnevt()   const { return fNtnevt; }
Int_t GetNtpevt()   const { return fNtpevt; }
Int_t GetJpnevt()   const { return fJpnevt; }
Int_t GetJppevt()   const { return fJppevt; }
Int_t GetJtnevt()   const { return fJtnevt; }
Int_t GetJtpevt()   const { return fJtpevt; }
Float_t GetXbjevt() const { return fXbjevt; }
Float_t GetQsqevt() const { return fQsqevt; }
Int_t GetNglevt()   const { return fNglevt; }
Float_t GetZppevt() const { return fZppevt; }
Float_t GetZptevt() const { return fZptevt; }

void SetBimevt(Float_t value) { fBimevt = value; }
void SetPhievt(Float_t value) { fPhievt = value; }
void SetKolevt(Int_t value)   { fKolevt = value; }
void SetKoievt(Int_t value)   { fKoievt = value; }
void SetPmxevt(Float_t value) { fPmxevt = value; }
void SetEgyevt(Float_t value) { fEgyevt = value; }
void SetNpjevt(Int_t value)   { fNpjevt = value; }
void SetNtgevt(Int_t value)   { fNtgevt = value; }
void SetNpnevt(Int_t value)   { fNpnevt = value; }
void SetNppevt(Int_t value)   { fNppevt = value; }
void SetNtnevt(Int_t value)   { fNtnevt = value; }
void SetNtpevt(Int_t value)   { fNtpevt = value; }
void SetJpnevt(Int_t value)   { fJpnevt = value; }
void SetJppevt(Int_t value)   { fJppevt = value; }
void SetJtnevt(Int_t value)   { fJtnevt = value; }
void SetJtpevt(Int_t value)   { fJtpevt = value; }
void SetXbjevt(Float_t value) { fXbjevt = value; }
void SetQsqevt(Float_t value) { fQsqevt = value; }
void SetNglevt(Int_t value)   { fNglevt = value; }
void SetZppevt(Float_t value) { fZppevt = value; }
void SetZptevt(Float_t value) { fZptevt = value; }


protected:

private:
	Float_t fBimevt; //   bimevt ........ absolute value of impact parameter
	Float_t fPhievt; //   phievt ........ angle of impact parameter
	Int_t fKolevt;   //   kolevt ........ number of collisions
	Int_t fKoievt;   //   koievt ........ number of inelastic collisions
	Float_t fPmxevt; //   pmxevt ........ reference momentum
	Float_t fEgyevt; //   egyevt ........ pp cm energy (hadron) or string energy (lepton)
	Int_t fNpjevt;   //   npjevt ........ number of primary projectile participants
	Int_t fNtgevt;   //   ntgevt ........ number of primary target participants
	Int_t fNpnevt;   //   npnevt ........ number of primary projectile neutron spectators
	Int_t fNppevt;   //   nppevt ........ number of primary projectile proton spectators
	Int_t fNtnevt;   //   ntnevt ........ number of primary target neutron spectators
	Int_t fNtpevt;   //   ntpevt ........ number of primary target proton spectators
	Int_t fJpnevt;   //   jpnevt ........ number of absolute projectile neutron spectators
	Int_t fJppevt;   //   jppevt ........ number of absolute projectile proton spectators
	Int_t fJtnevt;   //   jtnevt ........ number of absolute target neutron spectators
	Int_t fJtpevt;   //   jtpevt ........ number of absolute target proton spectators
	Float_t fXbjevt; //   xbjevt ........ bjorken x for dis
	Float_t fQsqevt; //   qsqevt ........ q**2 for dis
	Int_t fNglevt;   //   nglevt ........ number of collisions acc to  Glauber
	Float_t fZppevt; //   zppevt ........ average Z-parton-proj
	Float_t fZptevt; //   zptevt ........ average Z-parton-targ



	ClassDef(AliGenEposEventHeader,2)
};


#endif /* ALIGENEPOSEVENTHEADER_H_ */
