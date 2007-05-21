#ifndef ALIGENINFO_H
#define ALIGENINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliGenInfo                               //
//   collect together MC info for comparison purposes - effieciency studies and so on//                                                                 //
//   marian.ivanov@cern.ch                                                  //
//////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
//
// Start of implementation of the class AliTPCdigitRow
//
////////////////////////////////////////////////////////////////////////

#include <TParticle.h>
#include "AliTrackReference.h"

class TFile;
class AliRunLoader;
class AliStack;
class AliTPCParam;

const Int_t kgRowBytes = 32;

class AliTPCdigitRow: public TObject {
public:
  AliTPCdigitRow();
  virtual ~AliTPCdigitRow(){;}
  void SetRow(Int_t row);
  Bool_t TestRow(Int_t row) const ;
  AliTPCdigitRow & operator=(const AliTPCdigitRow &digOld);
  Int_t RowsOn(Int_t upto=8*kgRowBytes) const;
  Int_t Last() const;
  Int_t First() const ;
  void Reset();

//private:
  UChar_t fDig[kgRowBytes];
  ClassDef(AliTPCdigitRow,1)  // container for digit pattern
};


////////////////////////////////////////////////////////////////////////
//
// Start of implementation of the class AliMCInfo
//
////////////////////////////////////////////////////////////////////////

class AliMCInfo: public TObject {

public:
  AliMCInfo();
  ~AliMCInfo();   
  AliMCInfo(const AliMCInfo& info);
  void Update();


  AliTrackReference  fTrackRef;      // track reference saved in the output tree
  AliTrackReference  fTrackRefOut;   // decay track reference saved in the output tree
  AliTrackReference  fTRdecay;       // track reference at decay point
  //
  Int_t     fPrimPart;               // index of primary particle in TreeH
  TParticle fParticle;               // generated particle 
  Float_t   fMass;                   // mass of the particle
  Float_t   fCharge;                 //
  Int_t     fLabel;                  // track label
  Int_t     fEventNr;                // event number
  Int_t     fMCtracks;               // indication of how many times the track is retuturned back
  Int_t fPdg;                        //pdg code
  Float_t fDecayCoord[3];            // position of particle decay
  Double_t fVDist[4];                //distance of the particle vertex from primary vertex
  Bool_t fTPCdecay;                  //indicates decay in TPC
  Int_t fRowsWithDigitsInn;          // number of rows with digits in the inner sectors
  Int_t fRowsWithDigits;             // number of rows with digits in the outer sectors
  Int_t fRowsTrackLength;            // last - first row with digit
  Float_t fPrim;                     // theoretical dedx in tpc according particle momenta and mass
  AliTPCdigitRow fTPCRow;                  // information about digits row pattern
  Int_t fNTPCRef;                    // tpc references counter
  Int_t fNITSRef;                    // ITS references counter
  Int_t fNTRDRef;                    // TRD references counter
  Int_t fNTOFRef;                    // TOF references counter
  TClonesArray * fTPCReferences;     //containner with all track references -in the TPC
  TClonesArray * fITSReferences;     //container with ITS references
  TClonesArray * fTRDReferences;     //container with TRD references  
  TClonesArray * fTOFReferences;     //container with TRD references  
  //
  ClassDef(AliMCInfo,1)  // container for 
};



class AliGenV0Info: public TObject {
public:
  AliGenV0Info();       //
  AliMCInfo fMCd;       //info about daughter particle - second particle for V0
  AliMCInfo fMCm;       //info about mother particle   - first particle for V0
  TParticle fMotherP;   //particle info about mother particle
  void Update(Float_t vertex[3]);        // put some derived info to special field 
  Double_t    fMCDist1;    //info about closest distance according closest MC - linear DCA
  Double_t    fMCDist2;    //info about closest distance parabolic DCA
  //
  Double_t     fMCPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t     fMCPd[4];     //exact momentum from MC info
  Double_t     fMCX[3];      //exact position of the vertex
  Double_t     fMCXr[3];     //rec. position according helix
  //
  Double_t     fMCPm[3];    //momentum at the vertex mother
  Double_t     fMCAngle[3]; //three angels
  Double_t     fMCRr;       // rec position of the vertex 
  Double_t     fMCR;        //exact r position of the vertex
  Int_t        fPdg[2];   //pdg code of mother and daugter particles
  Int_t        fLab[2];   //MC label of the partecle  
  //
  Double_t       fInvMass;  //reconstructed invariant mass -
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  //
  ClassDef(AliGenV0Info,1)  // container for  
};



class AliGenKinkInfo: public TObject {
public:
  AliGenKinkInfo();          //default cosntructor
  void    Update();          // put some derived info to special field 
  Float_t GetQt();           //

  AliMCInfo   fMCd;          //info about daughter particle - second particle for V0
  AliMCInfo   fMCm;          //info about mother particle   - first particle for V0
  Double_t    fMCDist1;      //info about closest distance according closest MC - linear DCA
  Double_t    fMCDist2;      //info about closest distance parabolic DCA
  //
  Double_t     fMCPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t     fMCPd[4];     //exact momentum from MC info
  Double_t     fMCX[3];      //exact position of the vertex
  Double_t     fMCXr[3];     //rec. position according helix
  //
  Double_t     fMCPm[3];     //momentum at the vertex mother
  Double_t     fMCAngle[3];  //three angels
  Double_t     fMCRr;        // rec position of the vertex 
  Double_t     fMCR;         //exact r position of the vertex
  Int_t        fPdg[2];      //pdg code of mother and daugter particles
  Int_t        fLab[2];      //MC label of the partecle
  ClassDef(AliGenKinkInfo,1) // container for  
};



////////////////////////////////////////////////////////////////////////
// 
// Start of implementation of the class AliGenInfoMaker
//
////////////////////////////////////////////////////////////////////////

class AliGenInfoMaker {

public:
  AliGenInfoMaker();
  AliGenInfoMaker(const char * fnGalice, const char* fnRes    ="genTracks.root",
		   Int_t nEvents=1, Int_t firstEvent=0);
  virtual ~AliGenInfoMaker();
  Int_t Exec();
  Int_t Exec(Int_t nEvents, Int_t firstEventNr);
  void CreateTreeGenTracks();
  void CloseOutputFile();
  Int_t TreeKLoop();
  Int_t TreeTRLoop();
  Int_t TreeDLoop();
  Int_t BuildKinkInfo();  // build information about MC kinks
  Int_t BuildV0Info();  // build information about MC kinks
  Int_t BuildHitLines();  // build information about MC kinks
  void SetFirstEventNr(Int_t i) {fFirstEventNr = i;}
  void SetNEvents(Int_t i) {fNEvents = i;}
  void SetDebug(Int_t level) {fDebug = level;}
  Int_t SetIO(Int_t eventNr);
  Int_t CloseIOEvent();
  Int_t CloseIO();
  Int_t SetIO();
  Float_t TR2LocalX(AliTrackReference *trackRef,
		    AliTPCParam *paramTPC) const;
  AliMCInfo * GetInfo(UInt_t i) const {return (i<fNParticles)? fGenInfo[i]:0;}
  AliMCInfo * MakeInfo(UInt_t i);

private:
  Int_t  fDebug;                   //! debug flag  
  Int_t  fEventNr;                 //! current event number
  Int_t  fLabel;                   //! track label
  Int_t  fNEvents;                 //! number of events to process
  Int_t  fFirstEventNr;            //! first event to process
  UInt_t fNParticles;              //! number of particles in TreeK
  TTree *fTreeGenTracks;          //! output tree with generated tracks
  TTree *fTreeKinks;             //!  output tree with Kinks
  TTree *fTreeV0;                //!  output tree with V0
  TTree *fTreeHitLines;          //! tree with hit lines
  char   fFnRes[1000];             //! output file name with stored tracks
  TFile *fFileGenTracks;          //! output file with stored fTreeGenTracks
  //
  AliRunLoader * fLoader;         //! pointer to the run loader
  TTree * fTreeD;                 //! current tree with digits
  TTree * fTreeTR;                //! current tree with TR
  AliStack *fStack;               //! current stack
  // 
  AliMCInfo **   fGenInfo;    //! array with pointers to gen info
  Int_t   fNInfos;                  //! number of tracks with infos
  //
  AliTPCParam* fParamTPC;         //! AliTPCParam
  Float_t  fVPrim[3];             //! primary vertex position
                                  // the fVDist[3] contains size of the 3-vector
  // cuts
  //
  Double_t fTPCPtCut; // do not store particles with generated pT less than this
  Double_t fITSPtCut; // do not store particles with generated pT less than this
  Double_t fTRDPtCut; // do not store particles with generated pT less than this
  Double_t fTOFPtCut; // do not store particles with generated pT less than this
 
  ClassDef(AliGenInfoMaker,0)    // class which creates and fills tree with TPCGenTrack objects
};





AliTPCParam * GetTPCParam();
Float_t TPCBetheBloch(Float_t bg);

#endif
