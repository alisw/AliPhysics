/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#ifndef ALIJIAAANALYSIS_H
#define ALIJIAAANALYSIS_H

#include "../AliJDataManager.h"
#include "../AliJConst.h"
#include <TH1D.h>

// iaaAnalysis main class
// used in local and grid execution

class AliJCard;
class AliJIaaHistos;
class AliJIaaCorrelations;
class AliJEventHeader;
class AliJEventPool;
class AliJRunHeader;
class AliJEfficiency;
class AliJTrackCounter;
class AliJAcceptanceCorrection;

class TClonesArray;
class TF1;
class AliJRunTable;
class TRandom3;

class AliJIaaAnalysis : public TObject
{
  public:
    AliJIaaAnalysis();
    AliJIaaAnalysis(Bool_t execLocal);
    virtual ~AliJIaaAnalysis();
    AliJIaaAnalysis(const AliJIaaAnalysis& obj);
    AliJIaaAnalysis& operator=(const AliJIaaAnalysis& obj);
    
    void Initialize() const;
    void Init(){ Initialize(); }
    void UserCreateOutputObjects();
    void UserExec();
    void Terminate();
    
    //Int_t GetNumberEvents() const { return fnumberEvents; }
    AliJIaaHistos *GetHistos() { return fhistos; }
    AliJIaaCorrelations *GetCorrelations() { return fcorrelations; }
    AliJEventPool *GetAssocPool() { return fassocPool; }
    //Int_t GetEventCounter() const { return fEventCounter; }
    AliJCard *GetCard() { return fcard; }

    void SetCard( AliJCard *c ) { fcard = c; }
    void SetTrigger( char* p ) { fjtrigg = GetParticleType(p); }
    void SetAssoc( char* p ) { fjassoc = GetParticleType(p); }
    void SetInclusiveFile( const char *f ){ fInclusiveFile = f; }
    void SetInputFile( char *f ) { finputFile = f; }
    //void SetNumberEvents( Int_t n ) { fnumberEvents = n; }

    void SetTrackList( TClonesArray *a ) { fdmg->SetTrackList( a ); }
    void SetPhotonList( TClonesArray *a ) { fdmg->SetPhotonList( a ); }
    void SetCaloCellList( TClonesArray *a ) { fdmg->SetCaloCellList( a ); }
    void SetMCTrackList( TClonesArray *a ) { fdmg->SetMCTrackList( a ); }
    void SetHeaderList( TClonesArray *a ) { fdmg->SetHeaderList( a ); }
    void SetRunHeader( AliJRunHeader *a ) { frunHeader = a; }
    void SetRunInfoList( TList *a ) { fdmg->SetRunInfoList( a ); }
    //void SetESDVZERO( TObject *a ) { fdmg->SetESDVZERO( a ); }

    double DeltaPhi(double phi1, double phi2);
    particleType  GetParticleType(char *inchar);
    //void ScaleNotEquidistantHisto(TH1D *hid, const double sc);

  private:

    Bool_t fExecLocal; // exec mode
    Bool_t fFirstEvent; //!

    particleType fjtrigg; // assoc
    particleType fjassoc; // trigger

    AliJCard *fcard; // card
    char * finputFile; //!
    TString fInclusiveFile; // File for inclusive distributions

    Int_t fevt; // event number
    AliJIaaHistos *fhistos; //!
    AliJIaaCorrelations *fcorrelations; //!

    AliJAcceptanceCorrection *fAcceptanceCorrection; //! Class for acceptance correction
    AliJEventPool *fassocPool; //!
    TClonesArray *fphotonList; //!
    TClonesArray *fchargedHadronList; //!
    TClonesArray *fpizeroList; //!
    TClonesArray *ftriggList; //!
    TClonesArray *fassocList; //!
    TClonesArray *finputList; //!

    AliJDataManager* fdmg; //!
    AliJEventHeader *feventHeader; //!
    AliJRunHeader* frunHeader; //!

    double fcent; //!
    bool fbTriggCorrel; //!
    bool fbLPCorrel; //!

    double fMinimumPt; //!  Minimum pT value for a particle to be still accepted to analysis

    Int_t fEventBC; //!

    AliJEfficiency *fEfficiency;    // Efficience
    AliJRunTable *fRunTable;        // Run Table
    int fHadronSelectionCut;

    ClassDef(AliJIaaAnalysis, 1);

};

#endif // ALIJIAAANALYSIS_H
























