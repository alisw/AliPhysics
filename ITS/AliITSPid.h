#ifndef ALIITSPID_H
#define ALIITSPID_H
#include <TObject.h>
#include <TClonesArray.h>
#include <TVector.h>
#include "../TPC/AliTPCtrack.h"
#include "AliITSIOTrack.h"
#include "AliITStrackV2.h"
#include <TF1.h>
#include <assert.h>
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
	Int_t	GetPcode(TClonesArray*,Float_t);
	Int_t	GetPcode(Float_t,Float_t);
	Int_t   GetPcode(AliTPCtrack*);
        Int_t   GetPcode(AliITSIOTrack*); 
        Int_t   GetPcode(AliITStrackV2*);
	void	SetCut(Int_t,Float_t,Float_t,Float_t,
			    Float_t,Float_t,Float_t,Float_t);
	void    SetAProb(Int_t ivar,Int_t icut,Float_t apro){ aprob[ivar][icut]=apro; } 
	Float_t GetAProb(Int_t ivar,Int_t icut){ return aprob[ivar][icut]; } 
	Float_t GetWpi(){return fWpi;}
	Float_t GetWk(){return fWk;}
	Float_t GetWp(){return fWp;}
	Int_t	GetPid(){return fPcode;};
protected:
public:
	Float_t cut[13][7],aprob[3][8];
	Int_t       mxtrs;
	TClonesArray *trs;
	Float_t qtot;
	Float_t fWpi,fWk,fWp;
	Float_t fRpik,fRppi,fRpka,fRp; 
	Int_t 	fPcode;
//private:
public:
	int	qcomp(Float_t* qa,Float_t* qb){return qa[0]>qb[0]?1:0;}
	Float_t qtrm(Int_t track);
	Float_t qtrm(Float_t qarr[6],Int_t narr);
	Float_t fSigmin;
	Int_t	wpik(Float_t,Float_t);
        Int_t	wpikp(Float_t,Float_t);
	Int_t	pion(){return fWpi=1.,fPcode=211;}
	Int_t	kaon(){return fWk=1.,fPcode=321;}
	Int_t	proton(){return fWp=1.,fPcode=2212;}
        Int_t   fSilent;
	TF1*    fCutKa;
	TF1*    fCutPr;
	TF1*    fGGpi[6];
	TF1*    fGGka[3];
        TF1*    fGGpr[3];
	TF1*    fggpi;
	TF1*    fggka;
	TF1*    fggpr;
public:
  ClassDef(AliITSPid,1) // Class for ITS PID
};

#endif	




