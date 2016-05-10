/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#ifndef ALIJCORRAN_H
#define ALIJCORRAN_H

#include <AliJDataManager.h>
#include <AliJConst.h>
#include <TH1D.h>

// jcorran main class
// used in local and grid execution

class AliJCard;
class AliJHistos;
class AliJCorrelations;
class AliJEventHeader;
class AliJEventPool;
class AliJHistos;
class AliJRunHeader;
class AliJEfficiency;
class AliJTrackCounter;
class AliJAcceptanceCorrection;

class TClonesArray;
class TF1;
class AliJRunTable;
class TRandom3;

class AliJCORRAN : public TObject
{
  public:
    AliJCORRAN();
    AliJCORRAN(Bool_t execLocal);
    virtual ~AliJCORRAN();
    AliJCORRAN(const AliJCORRAN& obj);
    AliJCORRAN& operator=(const AliJCORRAN& obj);
    
    void Initialize() const;
    void Init(){ Initialize(); }
    void UserCreateOutputObjects();
    void UserExec();
    void Terminate();
    
    Int_t GetNumberEvents() const { return fnumberEvents; }
    AliJHistos *GetHistos() { return fhistos; }
    AliJCorrelations *GetCorrelations() { return fcorrelations; }
    AliJEventPool *GetAssocPool() { return fassocPool; }
    Int_t GetEventCounter() const { return fEventCounter; }
    AliJCard *GetCard() { return fcard; }

	void FillXtHistos(TClonesArray *inputList, AliJTrackCounter *lpTrackCounter); // Esko

	void SetCard( AliJCard *c ) { fcard = c; }
	void SetTrigger( char* p ) { fjtrigg = GetParticleType(p); }
	void SetAssoc( char* p ) { fjassoc = GetParticleType(p); }
	void SetInclusiveFile( const char *f ){ fInclusiveFile = f; }
	void SetInputFile( char *f ) { finputFile = f; }
	void SetNumberEvents( Int_t n ) { fnumberEvents = n; }

	void SetTrackList( TClonesArray *a ) { fdmg->SetTrackList( a ); }
	void SetPhotonList( TClonesArray *a ) { fdmg->SetPhotonList( a ); }
	void SetCaloCellList( TClonesArray *a ) { fdmg->SetCaloCellList( a ); }
	void SetMCTrackList( TClonesArray *a ) { fdmg->SetMCTrackList( a ); }
	void SetHeaderList( TClonesArray *a ) { fdmg->SetHeaderList( a ); }
	void SetRunHeader( AliJRunHeader *a ) { frunHeader = a; }
	void SetRunInfoList( TList *a ) { fdmg->SetRunInfoList( a ); }
	void SetESDVZERO( TObject *a ) { fdmg->SetESDVZERO( a ); }
 
	//     void SetESDTZERO( AliESDTZERO *a ) { fdmg->SetESDTZERO( a ); }
	//     void SetESDZDC( AliESDZDC *a ) { fdmg->SetESDZDC( a ); }


	double DeltaPhi(double phi1, double phi2);
	particleType  GetParticleType(char *inchar);
	void ScaleNotEquidistantHisto(TH1D *hid, const double sc);

	private:

	Bool_t fExecLocal; // exec mode
	Bool_t fFirstEvent; //!

	particleType fjtrigg; // assoc
	particleType fjassoc; // trigger

	AliJCard *fcard; // card
	char * finputFile; //!
	TString fInclusiveFile; // File for inclusive distributions

	Int_t fevt; // event number
	AliJHistos *fhistos; //!
	AliJCorrelations *fcorrelations; //!
  AliJAcceptanceCorrection *fAcceptanceCorrection; //! Class for acceptance correction
	AliJEventPool *fphotonPool; //!
	AliJEventPool *fassocPool; //! 
	TClonesArray *fphotonList; //! 
	TClonesArray *fchargedHadronList; //!
	TClonesArray *fpizeroList; //!
	TClonesArray *ftriggList; //! 
	TClonesArray *fassocList; //!
	TClonesArray *fpairList; //!
	TClonesArray *fpairCounterList; //!
	TClonesArray *finputList; //!

	AliJDataManager* fdmg; //!
	AliJEventHeader *feventHeader; //!
	AliJRunHeader* frunHeader; //!

	int fnumberEvents; //!
	int fieout; //!
	int fEventCounter; //!
	double fcent; //!
	Int_t fncBin; //!
	Int_t fnPttBin; //!
	bool fbTriggCorrel; //!
	bool fbLPCorrel; //!
	bool fbLPpairCorrel; //!
	double fTrackEtaRange; //!
	double flowerPtAssocBoarder; //!

	TF1 *fCentMultLow;  //!
	TF1 *fCentMultHigh; //!



	Int_t fEventBC; //!
	Float_t fSQRTS; // todo

	AliJEfficiency *fEfficiency; // comment me
	AliJRunTable *fRunTable; // comment needed
	float fIsolationR; // comment1
	int fHadronSelectionCut; /// comment2

	ClassDef(AliJCORRAN, 1); // EMCAL for jcorran

};

#endif
























