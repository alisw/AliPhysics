/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//
//     Data container for relative ITS-TPC alignment analysis
//     see info in the implementation file
//
//     Origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

class TObjArray;
class TGraphErrors;
class TNamed;
class TTree;
class TCollection;
class AliESDEvent;
class TBrowser;
class TList;
#include <AliRelAlignerKalman.h>

class AliRelAlignerKalmanArray:public TNamed
{
public:
  AliRelAlignerKalmanArray();
  AliRelAlignerKalmanArray(Int_t t0, Int_t tend, Int_t slotwidth);
  virtual ~AliRelAlignerKalmanArray();
  AliRelAlignerKalmanArray& operator=(const AliRelAlignerKalmanArray& a );
  AliRelAlignerKalmanArray(const AliRelAlignerKalmanArray& a);
  void SetupArray(Int_t t0, Int_t tend, Int_t slotwidth);
  
  AliRelAlignerKalman* GetAligner(UInt_t timestamp);
  AliRelAlignerKalman* GetAligner(AliESDEvent* event);
  AliRelAlignerKalman* GetAlignerTemplate();
  Long64_t Merge( TCollection* list );
  AliRelAlignerKalman* At( Int_t i ) const;
  AliRelAlignerKalman* operator[](Int_t i) const;
  AliRelAlignerKalman*& operator[](Int_t i);
  Int_t GetEntries() const;
  Int_t GetSize() const {return fSize;}
  AliRelAlignerKalman* Last() const;
  UInt_t GetT0() const {return fT0;}
  UInt_t GetTimebinWidth() const {return fTimebinWidth;}
  Int_t Timebin( UInt_t timestamp ) const;
  virtual void Print(Option_t* option="") const;
  void FillTree( TTree* tree )const ;
  TGraphErrors* MakeGraph(Int_t iparam) const;
  AliRelAlignerKalmanArray* MakeSmoothArray() const;
  void SetOutRejSigmaOnMerge(Double_t s) {fOutRejSigmaOnMerge=s;}
  void SetOutRejSigmaOnSmooth(Double_t s) {fOutRejSigmaOnSmooth=s;}
  void Browse(TBrowser *b);

private:
  void ClearContents();
  void PropagateToTime(AliRelAlignerKalman* al, UInt_t timestamp ) const;

  UInt_t fT0;                            //time of first time slot
  Int_t fTimebinWidth;                   //width of the time bin in seconds
  Int_t fSize;                           //size
  Double_t fOutRejSigmaOnMerge;          //how much outlier rejection on merge
  Double_t fOutRejSigmaOnSmooth;          //how much outlier rejection on Smooth
  AliRelAlignerKalman fAlignerTemplate;  //template
  AliRelAlignerKalman** fPArray;         //[fSize] an array of aligners
  TList* fListOfGraphs;                  //!hold the graphs
  
  ClassDef(AliRelAlignerKalmanArray,5);   //AliRelAlignerKalman class
};


