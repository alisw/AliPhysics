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

#ifndef AliRelAlignerKalmanArray_h
#define AliRelAlignerKalmanArray_h


class TString;
class TCollection;
class AliESDEvent;
class TObjArray;
class AliRelAlignerKalman;
class TNamed;

class AliRelAlignerKalmanArray:public TNamed
{
public:
  AliRelAlignerKalmanArray();
  AliRelAlignerKalmanArray(const char* name);
  virtual ~AliRelAlignerKalmanArray();
  AliRelAlignerKalmanArray& operator=(const AliRelAlignerKalmanArray& a );
  AliRelAlignerKalmanArray(const AliRelAlignerKalmanArray& a);
  
  Long64_t Merge( TCollection* list );
  //Bool_t AddESDEvent( AliESDEvent* event );
  Bool_t AddCosmicEvent( AliESDEvent* event );
  void AddLast( AliRelAlignerKalman* al );
  AliRelAlignerKalman* At( Int_t i ) const;
  AliRelAlignerKalman* Last() const;
  Int_t GetEntries() const {return fArray->GetEntriesFast();}
  AliRelAlignerKalman* operator[](Int_t i) const;
  Bool_t SetTimeMatchingTolerance( const UInt_t m );
  Bool_t SetSaveInterval( const UInt_t s );
  UInt_t GetTimeMatchingTolerance() const {return fTimeMatchingTolerance;}
  UInt_t GetSaveInterval() const {return fSaveInterval;}
  UInt_t TimeBin( UInt_t timebin ) const;
  void SetCurrentTimeBin( UInt_t timestamp );
  UInt_t GetCurrentTimeBin() const {return fCurrentTimeBin;}
  Bool_t IsInCurrentTimeBin( UInt_t timestamp ) const;
  AliRelAlignerKalman* GetAligner() const {return fAligner;}
  TObjArray* SortedMerge ( TObjArray* input ); 
  //void SetResetAllAtNewRun( Bool_t s ) {fResetAllAtNewRun = s;}
  //void SetResetTPCAtNewRun( Bool_t s ) {fResetTPCAtNewRun = s;}

private:
  TObjArray* fArray; //an array of aligners
  UInt_t fSaveInterval; //how often to save (in seconds)
  UInt_t fTimeMatchingTolerance; //tolerance for matching timestamps
  UInt_t fCurrentTimeBin; //current timebin
  AliRelAlignerKalman* fAligner;  //aligner object
  //Bool_t fResetAllAtNewRun;
  //Bool_t fResetTPCAtNewRun;
  
  ClassDef(AliRelAlignerKalmanArray,1)     //AliRelAlignerKalman class
};

#endif

