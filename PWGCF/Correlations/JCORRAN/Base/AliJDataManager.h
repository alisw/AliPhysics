// $Id: AliJDataManager.h,v 1.8 2008/02/12 15:51:27 djkim Exp $
////////////////////////////////////////////////////
/*!
  \file AliJDataManager.h
  \brief
  \author J. Rak, D.J.Kim, B.S Chang (University of Jyvaskyla)
  \email: djkim@cc.jyu.fi
  \version $Revision: 1.8 $
  \date $Date: 2008/02/12 15:51:27 $
  */
////////////////////////////////////////////////////

#ifndef ALIJDATAMANAGER_H
#define ALIJDATAMANAGER_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TChain.h>

#include "AliJConst.h"
#include "AliJEventHeader.h"
#include "AliJRunHeader.h"

class TClonesArray;
class TList;
class AliJHistos;
class AliJCard;
class AliJCorrelations;
class AliJDataManager;
class AliJRunHeader;
class AliJTrackCut;
//class AliESDVZERO;TODO



class AliJDataManager  {

public:
  enum  { kTrackCutDefault=0, kTrackCutJFG=1, kTrackCutHBT=2,  kTrackCutRaa=3, kTrackCutHybrid=4, kTrackCutpp=5 };
  AliJDataManager(AliJCard *inCard, AliJHistos *fhistos, AliJCorrelations *corrin, Bool_t execLocal);

  virtual ~AliJDataManager();                        //destructor
  AliJDataManager();

  void ChainInputStream(const char* infileList);
  Double_t LoadEvent( int ievt );
  virtual void RegisterList(TClonesArray* listToFill, TClonesArray* listFromToFill, int cBin, int zBin, 
      particleType whatToFill);
  // corrType whatCorrType);

  virtual bool IsGoodEvent();
  void SetRunHeader(AliJRunHeader *runHeader){ 
	  std::cout<<"DEBUG ########### "<<runHeader<<std::endl;
	  std::cout<<"DEBUG ########### "<<runHeader->GetRunNumber()<<std::endl;
	  fRunHeader = runHeader;
  }

  // GETTER
  TChain * GetChain(){ return fChain; };
  AliJCard * GetCard(){ return fCard; };
  int GetNEvents(){ return fChain->GetEntries(); } 
  AliJRunHeader *GetRunHeader(){ return fRunHeader; }; 
  AliJEventHeader * GetEventHeader(){ return fEventHeader; };
  TClonesArray  *GetEventHeaderList(){ return fEventHeaderList; };
  TClonesArray  *GetTrackList(){ return fTrackList; };
  TClonesArray  *GetPhotonList(){ return fPhotonList; };
  TClonesArray  *GetCellList(){ return fCellList; };
  TClonesArray  *GetPhotonListRecalib(){ return fPhotonListRecalib; };
  TClonesArray  *GetCellListRecalib(){ return fCellListRecalib; };
  TClonesArray  *GetMCTrackList(){ return fMCTrackList; };
  //AliESDVZERO *GetVZERO(){ return fVZEROData; }
  TObject *GetVZERO(){ return fVZEROData; }

  void SetExecLocal( const Bool_t b ) { fExecLocal = b; }

  void SetTrackList( TClonesArray *a ) { fTrackList = a; }
  void SetPhotonList( TClonesArray *a ) { fPhotonList = a; }
  void SetCaloCellList( TClonesArray *a ) { fCellList = a; }
  void SetMCTrackList( TClonesArray *a ) { fMCTrackList = a; }
  void SetHeaderList( TClonesArray *a ) { fEventHeaderList = a; }
  void SetRunInfoList( TList *a ) { fRunInfoList = a; }
  void SetESDVZERO( TObject *a ) { fVZEROData = a; }
  //     void SetESDTZERO( AliESDTZERO *a ) { AliESDTZERO = a; }
  //     void SetESDZDC( AliESDZDC *a ) { AliESDZDC = a; }

  UInt_t GetFilterMap() const { return fFilterMap; }
  void SetFilterMap( UInt_t fm ){ fFilterMap = fm; }
  void AcceptTrackBit( UInt_t bit ){ SETBIT(fFilterMap, bit); } 
  Bool_t TestTrackBit( int i ){ return TESTBIT( GetFilterMap(), i );}

  Bool_t IsSelectedTrigger(unsigned int triggin ) const { return (triggin & fTriggerMask) >0;}

	protected:
  TChain * fChain; // comment me
  AliJCard * fCard; // comment me
  AliJHistos *fhistos; // comment me
  AliJCorrelations* fcorrelations; // comment me

  // Alice specific
  AliJRunHeader *fRunHeader;  // comment me
  AliJEventHeader * fEventHeader; // comment me
  TClonesArray  *fEventHeaderList; // comment me
  TClonesArray  *fTrackList; // comment me
  TClonesArray  *fPhotonList; // comment me
  TClonesArray  *fCellList; // comment me
  TClonesArray  *fPhotonListRecalib; // comment me
  TClonesArray  *fCellListRecalib; // comment me
  TClonesArray  *fMCTrackList; // comment me
  TObject  *fVZEROData; // comment me
  TList *fRunInfoList; // comment me

  int fhadronSelectionCut; // comment me

  UInt_t fFilterMap; // comment me

  TString fFName; // comment me
  Bool_t fExecLocal; // comment me

  int      fTriggerMask; // Trigger mask;


  AliJTrackCut * fTrackCut; // Trigg selection class

  private:
  AliJDataManager(const AliJDataManager& obj) ;
  AliJDataManager& operator=(const AliJDataManager& obj) {return *this;};

};

#endif






















