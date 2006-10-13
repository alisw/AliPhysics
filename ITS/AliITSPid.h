#ifndef ALIITSPID_H
#define ALIITSPID_H
/////////////////////////////////////////////////////////////////
// Class for identification of pions,kaons and protons in ITS  //
// Prior particles population (probabilities) are taken from   //
// Hijing event generator.                                     //
/////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TVectorfwd.h>

class TClonesArray;
class AliITSIOTrack;
class AliKalmanTrack;
class AliITStrackV2;
class TF1;


class  AliITSPid : public TObject {

public:
		AliITSPid(Int_t ntrs=1000);
		virtual ~AliITSPid(){}

	void	SetEdep(Int_t track,Float_t Edep);
	void	SetPmom(Int_t track,Float_t Pmom);
	void	SetPcod(Int_t track,Int_t Pcod);
	void	Print(Int_t track);
        virtual void Print(Option_t *option="") const {TObject::Print(option);}
	void	Tab(void);
	void    Reset(void);
	void	SetVec(Int_t track,const TVector& info) const;
	TVector* GetVec(Int_t track) const;
	Int_t	GetPcode(TClonesArray* rps,Float_t pm);
	Int_t	GetPcode(Float_t p,Float_t pm);
	Int_t   GetPcode(AliKalmanTrack* track);
        Int_t   GetPcode(AliITSIOTrack* track); 
        Int_t   GetPcode(AliITStrackV2* track);
	void	SetCut(Int_t n,Float_t pm,Float_t pilo,Float_t pihi,
		       Float_t klo,Float_t khi,Float_t plo,Float_t phi);
	void    SetAProb(Int_t ivar,Int_t icut,Float_t apro){ fAprob[ivar][icut]=apro; } 
	Float_t GetAProb(Int_t ivar,Int_t icut) const { return fAprob[ivar][icut]; } 
	Float_t GetWpi() const {return fWpi;}
	Float_t GetWk() const {return fWk;}
	Float_t GetWp() const {return fWp;}
	Int_t	GetPid() const {return fPcode;};
protected:
	// copy constructor and assignment operator are protected
	// since they are not allowed
        AliITSPid(const AliITSPid &source); // copy constructor. 
        AliITSPid& operator=(const AliITSPid&  source); // = operator.

	int	Qcomp(Float_t* qa,Float_t* qb) const {return qa[0]>qb[0]?1:0;}
	Float_t Qtrm(Int_t track);
	Float_t Qtrm(Float_t qarr[6],Int_t narr) const;
	Int_t	Wpik(Float_t pm,Float_t q);
        Int_t	Wpikp(Float_t pm,Float_t q);
	Int_t	Pion(){return fWpi=1.,fPcode=211;}
	Int_t	Kaon(){return fWk=1.,fPcode=321;}
	Int_t	Proton(){return fWp=1.,fPcode=2212;}
	//================ Data members ========================
	Float_t fCut[13][7],fAprob[3][8]; //Cuts and prior probs tables
	Int_t       fMxtrs; //Maximum tracks limit
	TClonesArray *fTrs; //Tracks set under investigation
	Float_t fWpi,fWk,fWp; //Probabilities for pions,kaons,protons        
	Float_t fRpik,fRppi,fRpka,fRp; //Signal ratios
	Int_t 	fPcode;  //Particle code
	Float_t fSigmin; // Tuning parameter
        Int_t   fSilent; // Output suppresion flag
	TF1*    fCutKa;  // Pions-kaons cut function
	TF1*    fCutPr;  // Kaons-protons cut function
	TF1*    fGGpi[6];// Pions signal parametrization for Hijing
	TF1*    fGGka[3];// Kaons          --//--
        TF1*    fGGpr[3];// Protons        --//--
	TF1*    fggpi;   // Pions signal for given momentum
	TF1*    fggka;   // Kaons          --//--
	TF1*    fggpr;   // Protons        --//--
  ClassDef(AliITSPid,2) // Class for ITS PID
};

#endif	




