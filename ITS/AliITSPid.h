#ifndef ALIITSPID_H
#define ALIITSPID_H

#include <TObject.h>
#include "AliITSIOTrack.h"

class TClonesArray;
class TVector;
class AliTPCtrack;
class AliITStrackV2;
class TF1;

//___________________________________________________________________________
class  AliITSPid :
  public TObject {

public:
		AliITSPid(Int_t ntrs=1000);
		virtual ~AliITSPid(){}
	void	SetEdep(Int_t track,Float_t Edep);
	void	SetPmom(Int_t track,Float_t Pmom);
	void	SetPcod(Int_t track,Int_t Pcod);
	void	Print(Int_t track);
	void	Tab(void);
	void    Reset(void);
	void	SetVec(Int_t track,TVector info);
	TVector* GetVec(Int_t track);
	Int_t	GetPcode(TClonesArray* rps,Float_t pm);
	Int_t	GetPcode(Float_t p,Float_t pm);
	Int_t   GetPcode(AliTPCtrack* track);
        Int_t   GetPcode(AliITSIOTrack* track); 
        Int_t   GetPcode(AliITStrackV2* track);
	void	SetCut(Int_t n,Float_t pm,Float_t pilo,Float_t pihi,
		       Float_t klo,Float_t khi,Float_t plo,Float_t phi);
	void    SetAProb(Int_t ivar,Int_t icut,Float_t apro){ fAprob[ivar][icut]=apro; } 
	Float_t GetAProb(Int_t ivar,Int_t icut){ return fAprob[ivar][icut]; } 
	Float_t GetWpi(){return fWpi;}
	Float_t GetWk(){return fWk;}
	Float_t GetWp(){return fWp;}
	Int_t	GetPid(){return fPcode;};
protected:
	Float_t fCut[13][7],fAprob[3][8];
	Int_t       fMxtrs;
	TClonesArray *fTrs;
	Float_t fWpi,fWk,fWp;
	Float_t fRpik,fRppi,fRpka,fRp; 
	Int_t 	fPcode;

	int	Qcomp(Float_t* qa,Float_t* qb){return qa[0]>qb[0]?1:0;}
	Float_t Qtrm(Int_t track);
	Float_t Qtrm(Float_t qarr[6],Int_t narr);
	Float_t fSigmin;
	Int_t	Wpik(Float_t pm,Float_t q);
        Int_t	Wpikp(Float_t pm,Float_t q);
	Int_t	Pion(){return fWpi=1.,fPcode=211;}
	Int_t	Kaon(){return fWk=1.,fPcode=321;}
	Int_t	Proton(){return fWp=1.,fPcode=2212;}
        Int_t   fSilent;
	TF1*    fCutKa;
	TF1*    fCutPr;
	TF1*    fGGpi[6];
	TF1*    fGGka[3];
        TF1*    fGGpr[3];
	TF1*    fggpi;
	TF1*    fggka;
	TF1*    fggpr;

  ClassDef(AliITSPid,2) // Class for ITS PID
};

#endif	




