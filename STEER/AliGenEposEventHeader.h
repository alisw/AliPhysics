/*
 * AliGenEposEventHeader.h
 * 
 * Header for EPOS generated event.
 *
 *      Author: Piotr Ostrowski
 */

#ifndef ALIGENEPOSEVENTHEADER_H_
#define ALIGENEPOSEVENTHEADER_H_

#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
class TGenerator;


class AliGenEposEventHeader : public AliGenEventHeader, public AliCollisionGeometry
{
public:
	AliGenEposEventHeader(const char* name);
	AliGenEposEventHeader();
	virtual ~AliGenEposEventHeader() {}

	Float_t GetBimevt();
	Float_t GetPhievt();
	Int_t GetKolevt();
	Int_t GetKoievt();
	Float_t GetPmxevt();
	Float_t GetEgyevt();
	Int_t GetNpjevt();
	Int_t GetNtgevt();
	Int_t GetNpnevt();
	Int_t GetNppevt();
	Int_t GetNtnevt();
	Int_t GetNtpevt();
	Int_t GetJpnevt();
	Int_t GetJppevt();
	Int_t GetJtnevt();
	Int_t GetJtpevt();
	Float_t GetXbjevt();
	Float_t GetQsqevt();
	Int_t GetNglevt();
	Float_t GetZppevt();
	Float_t GetZptevt();

	void FillInternalFields(TGenerator *epos);

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
