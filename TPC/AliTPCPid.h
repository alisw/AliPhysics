#ifndef ALITPCPID_H
#define ALITPCPID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
//#include <TF1.h>
//#include <TClonesArray.h>
#include <TVector.h>
#include "../TPC/AliTPCtrack.h"
//#include "../ITS/AliITSIOTrack.h"
class AliKalmanTrack;
#include <assert.h>

class TF1;
class TClonesArray;
//___________________________________________________________________________
class  AliTPCPid :
  public TObject {

public:
                AliTPCPid(Int_t ntrs=1000);
		virtual ~AliTPCPid(){}
                AliTPCPid( const AliTPCPid& r);
                AliTPCPid &operator = (const AliTPCPid & param); //assignment
	void	SetEdep(Int_t track,Float_t Edep);
	void	SetPmom(Int_t track,Float_t Pmom);
	void	SetPcod(Int_t track,Int_t Pcod);
	void	PrintPID(Int_t track);
	void	Tab(void);
	void    Reset(void);
	void	SetVec(Int_t track,TVector info) const;
	TVector* GetVec(Int_t track) const;
	Int_t	GetPcode(TClonesArray* ,Float_t) const;
	Int_t	GetPcode(Float_t q,Float_t pm);
	Int_t   GetPcode(AliTPCtrack*track);
        //Int_t   GetPcode(AliKalmanTrack* track);
	void	SetCut(Int_t n, Float_t pm, Float_t pilo, Float_t pihi,
			    Float_t klo, Float_t khi, Float_t plo, 
                            Float_t phi);
	void    SetAProb(Int_t ivar,Int_t icut,Float_t apro){ faprob[ivar][icut]=apro; } 
	Float_t GetAProb(Int_t ivar,Int_t icut) const
                        { return faprob[ivar][icut]; } 
	Float_t GetWpi() const {return fWpi;}
	Float_t GetWk() const {return fWk;}
	Float_t GetWp() const {return fWp;}
	Int_t	GetPid() const {return fPcode;};
	Float_t Qcorr(Float_t xc);
	Int_t	Qcomp(Float_t* qa,Float_t* qb) const {return qa[0]>qb[0]?1:0;}
	Float_t Qtrm(Int_t track) const;
	Float_t Qtrm(Float_t qarr[6],Int_t narr);
	Int_t	Wpik(Int_t nc, Float_t q);
	Int_t	Wpikp(Int_t nc, Float_t q);
	Int_t	Pion(){return /*fWpi=1.,*/fPcode=211;}
	Int_t	Kaon(){return /*fWk=1.,*/fPcode=321;}
	Int_t	Proton(){return /*fWp=1.,*/fPcode=2212;}
private:
	TF1 *fCutKa; // function
	TF1 *fCutPr; // function
	Float_t fCutKaTune,fCutPrTune; // tune cuts
	Float_t fSigmin; // sigma min
        Int_t   fSilent; // flag
	Float_t fcut[13][7],faprob[3][8]; //cuts
	Int_t       fmxtrs; // fmxtrs
	TClonesArray *trs; //pointer
	Float_t fqtot; // tot q
	Float_t fWpi,fWk,fWp; // weights
	Float_t fRpik,fRppi,fRpka,fRp; // ratios
	Int_t 	fPcode; //p-code

  ClassDef(AliTPCPid,1) // Class for TPC PID

};

#endif	




