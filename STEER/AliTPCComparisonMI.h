

////////////////////////////////////////////////////////////////////////
//
// Start of implementation of the class digitRow
//
////////////////////////////////////////////////////////////////////////
const Int_t kgRowBytes = 32;

class digitRow: public TObject {

public:
  digitRow();
  virtual ~digitRow(){;}
  void SetRow(Int_t row);
  Bool_t TestRow(Int_t row);
  digitRow & operator=(const digitRow &digOld);
  Int_t RowsOn(Int_t upto=8*kgRowBytes);
  Int_t Last();
  Int_t First();
  void Reset();

//private:
  UChar_t fDig[kgRowBytes];

  ClassDef(digitRow,1)  // container for digit pattern
};
ClassImp(digitRow)


////////////////////////////////////////////////////////////////////////
//
// Start of implementation of the class AliTPCGenInfo
//
////////////////////////////////////////////////////////////////////////

class AliTPCGenInfo: public TObject {

public:
  AliTPCGenInfo();
  ~AliTPCGenInfo();

  AliTrackReference fTrackRef;      // track reference saved in the output tree
  AliTrackReference fTrackRefOut;   // decay track reference saved in the output tree
  TParticle fParticle;           // generated particle 
  Int_t fLabel;                   // track label
  Int_t fEventNr;                 // event number

  Float_t fDecayCoord[3];         // position of particle decay
  Double_t fVDist[4];             //distance of the particle vertex from primary vertex
  
  Int_t fRowsWithDigitsInn;    // number of rows with digits in the inner sectors
  Int_t fRowsWithDigits;       // number of rows with digits in the outer sectors
  Int_t fRowsTrackLength;      // last - first row with digit
  Int_t fDigitsInSeed;         // digits in the default seed rows
  Float_t fPrim;               // theoretical dedx in tpc according particle momenta and mass
  digitRow fRow;               // information about digits row pattern
  TClonesArray * fReferences;  //containner with all track references
  ClassDef(AliTPCGenInfo,1)  // container for 
};
ClassImp(AliTPCGenInfo)



class AliTPCGenV0Info: public TObject {
public:
  AliTPCGenInfo fMCd;      //info about daughter particle
  AliTPCGenInfo fMCm;      //info about mother particle
  void Update();        // put some derived info to special field 
  Double_t    fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t    fDist2;    //info about closest distance parabolic DCA
  //
  Double_t     fPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t     fPd[4];     //exact momentum from MC info
  Double_t     fX[3];      //exact position of the vertex
  Double_t     fXr[3];     //rec. position according helix
  //
  Double_t     fPm[3];    //momentum at the vertex mother
  Double_t     fAngle[3]; //three angels
  Double_t     fRr;       // rec position of the vertex 
  Double_t     fR;        //exact r position of the vertex
  Int_t        fPdg[2];   //pdg code of mother and daugter particles
  Int_t        fLab[2];   //MC label of the partecle
  ClassDef(AliTPCGenV0Info,1)  // container for  
};
ClassImp(AliTPCGenV0Info)



void AliTPCGenV0Info::Update()
{
  fPd[0] = fMCd.fParticle.Px();
  fPd[1] = fMCd.fParticle.Py();
  fPd[2] = fMCd.fParticle.Pz();
  fPd[3] = fMCd.fParticle.P();
  fX[0]  = fMCd.fParticle.Vx();
  fX[1]  = fMCd.fParticle.Vy();
  fX[2]  = fMCd.fParticle.Vz();
  fR     = TMath::Sqrt( fX[0]*fX[0]+
				 fX[1]*fX[1]);
  fPdg[0]    = fMCd.fParticle.GetPdgCode();
  fPdg[1]    = fMCm.fParticle.GetPdgCode();
  //
  fLab[0]    = fMCd.fParticle.GetUniqueID();
  fLab[1]    = fMCm.fParticle.GetUniqueID();

}




     /////////////////////////////////////////////////////////////////////////
class AliTPCRecInfo: public TObject {

public:
  AliTPCRecInfo(){fTP = new TClonesArray("AliTPCTrackPoint2",0);}
  ~AliTPCRecInfo(){if (fTP) {fTP->Delete();delete fTP;}}
  //
  AliTPCtrack fTPCTrack;          // tpc track
  Float_t fTRLocalCoord[3];       //local coordinates of the track ref.
  Int_t   fReconstructed;         //flag if track was reconstructed
  Double_t fRecPhi;         // reconstructed phi angle (0;2*kPI)
  Double_t fLambda;         // reconstructed 
  Double_t fRecPt_1;        // reconstructed 
  Float_t fdEdx;           // reconstructed  dEdx      
  Int_t fFake;             // fake track
  Int_t fMultiple;         // number of reconstructions
  TClonesArray *fTP;        //container with track  points 
  void Reset();
  //
  ClassDef(AliTPCRecInfo,1)  // container for 
};
ClassImp(AliTPCRecInfo)

void AliTPCRecInfo::Reset()
{
  fMultiple =0; 
  fFake     =0;
  fReconstructed=0;
  fRecPhi =0;
  fLambda =0;
}


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


class AliTPCRecV0Info: public TObject {
public:
  AliTPCRecInfo  fT1;      //track1
  AliTPCRecInfo  fT2;      //track2  
  Double_t    fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t    fDist2;    //info about closest distance parabolic DCA
  //
  Double_t     fPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t     fXr[3];     //rec. position according helix
  //
  Double_t     fPm[3];    //momentum at the vertex mother
  Double_t     fAngle[3]; //three angles
  Double_t     fRr;       // rec position of the vertex 
  Int_t        fLab[2];   //MC label of the partecle
  ClassDef(AliTPCRecV0Info,1)  // container for  
};

ClassImp(AliTPCRecV0Info)



  


////////////////////////////////////////////////////////////////////////
// 
// Start of implementation of the class TPCFindGenTracks
//
////////////////////////////////////////////////////////////////////////

class TPCFindGenTracks {

public:
  TPCFindGenTracks();
  TPCFindGenTracks(char* fnHits,
		   char* fnDigits ="tpc.digits.root",
		   char* fnRes    ="genTracks.root",
		   Int_t nEvents=1, Int_t firstEvent=0);
  virtual ~TPCFindGenTracks();
  void Reset();
  Int_t Exec();
  Int_t Exec(Int_t nEvents, Int_t firstEventNr);
  void CreateTreeGenTracks();
  void CloseOutputFile();
  Int_t TreeKLoop();
  Int_t TreeTRLoop();
  Int_t TreeDLoop();
  //void  FillInfo(Int_t iParticle);
  void SetFirstEventNr(Int_t i) {fFirstEventNr = i;}
  void SetNEvents(Int_t i) {fNEvents = i;}
  void SetDebug(Int_t level) {fDebug = level;}
  Int_t SetIO(Int_t eventNr);
  Int_t SetIO();
  Float_t TR2LocalX(AliTrackReference *trackRef,
		    AliTPCParam *paramTPC);

public:
  AliTPCGenInfo*  fMCInfo;           //! information writen per particle
  Int_t fDebug;                   //! debug flag  
  Int_t fEventNr;                 //! current event number
  Int_t fLabel;                   //! track label
  Int_t fNEvents;                 //! number of events to process
  Int_t fFirstEventNr;            //! first event to process
  Int_t fNParticles;              //! number of particles in TreeK
  TTree *fTreeGenTracks;          //! output tree with generated tracks
  char *fFnRes;                   //! output file name with stored tracks
  char *fFnHits;                  //! input file name with hits
  char *fFnDigits;                //! input file name with digits
  TFile *fFileGenTracks;          //! output file with stored fTreeGenTracks
  TFile *fFileHits;               //! input file with hits
  TFile *fFileTreeD;              //! input file with digits
  //
  TTree * fTreeD;                 //! current tree with digits
  TTree * fTreeTR;                //! current tree with TR
  AliStack *fStack;               //! current stack
  //
  digitRow *fContainerDigitRow;   //! big container for partial information
  //
  AliTrackReference *fReferences; //! container with track references
  Int_t *fReferenceIndex0;        //! first index for given track
  Int_t *fReferenceIndex1;        //! last  index for given track
  //
  AliTPCParam* fParamTPC;         //! AliTPCParam
  Double_t fVPrim[3];             //! primary vertex position
                                  // the fVDist[3] contains size of the 3-vector

private:

// some constants for the original non-pareller tracking (by Y.Belikov)
  static const Int_t seedRow11 = 158;  // nRowUp - 1
  static const Int_t seedRow12 = 139;  // nRowUp - 1 - (Int_t) 0.125*nRowUp
  static const Int_t seedRow21 = 149;  // seedRow11 - shift
  static const Int_t seedRow22 = 130;  // seedRow12 - shift
  static const Double_t kRaddeg = 180./kPI;

  static const Int_t fgMaxIndexTR = 50000; // maximum number of tracks with a track ref
  static const Int_t fgMaxTR = 1000000; // maximum number of  track refs

  static const Int_t fgMaxParticles = 2000000; // maximum number of generated particles
  static const Double_t fgPtCut = .1; // do not store particles with generated pT less than this
  static const Float_t fgTrackRefLocalXMax = 82.95;
  static const Float_t fgTrackRefLocalXMaxDelta = 5.;

  ClassDef(TPCFindGenTracks,1)    // class which creates and fills tree with TPCGenTrack objects
};
ClassImp(TPCFindGenTracks)



////////////////////////////////////////////////////////////////////////
// 
// Start of implementation of the class TPCCmpTr
//
////////////////////////////////////////////////////////////////////////

class TPCCmpTr {

public:
  TPCCmpTr();
  TPCCmpTr(char* fnRecTracks,
	   char* fnGenTracks   ="genTracks.root",
	   char* fnCmpRes      ="cmpTracks.root", 
	   char* fnGalice      ="galice.root",
	   Int_t nEvents=1, Int_t firstEvent=0);
  virtual ~TPCCmpTr();
  void Reset();
  Int_t Exec();
  Int_t Exec(Int_t nEvents, Int_t firstEventNr);
  void CreateTreeCmp();
  void CloseOutputFile();
  Bool_t ConnectGenTree();
  Int_t TreeGenLoop(Int_t eventNr);
  Int_t TreeTLoop(Int_t eventNr);
  void SetFirstEventNr(Int_t i) {fFirstEventNr = i;}
  void SetNEvents(Int_t i) {fNEvents = i;}
  void SetDebug(Int_t level) {fDebug = level;}

// tmp method, should go to TrackReferenceTPC
  TVector3 TR2Local(AliTrackReference *trackRef,
		    AliTPCParam *paramTPC);

private:

  Int_t fEventNr;                 //! current event number
  Int_t fNEvents;                 //! number of events to process
  Int_t fFirstEventNr;            //! first event to process
  //
  char *fFnCmp;                   //! output file name with cmp tracks
  TFile *fFileCmp;                //! output file with cmp tracks
  TTree *fTreeCmp;                //! output tree with cmp tracks
  //
  char *fFnGenTracks;             //! input file name with gen tracks
  TFile *fFileGenTracks;
  TTree *fTreeGenTracks;
  //
  char *fFnHits;                  //! input file name with gAlice object (needed for B)
  TFile *fFileHits;               //! input file with gAlice
  //
  char *fFnRecTracks;             //! input file name with tpc rec. tracks
  TFile *fFileRecTracks;          //! input file with reconstructed tracks
  TTree *fTreeRecTracks;          //! tree with reconstructed tracks
  TTree *fTreePoints;             //! tree with reconstructed points
  //
  Int_t *fIndexRecTracks;         //! index of particle label in the TreeT_TPC
  Int_t *fFakeRecTracks;          //! number of fake tracks
  Int_t *fMultiRecTracks;         //! number of multiple reconstructions
  //
  TObjArray * fTracks;            //!container with tracks 
  TObjArray * fTrackPoints;       //! container with track points
  //
  AliTPCParam* fParamTPC;         //! AliTPCParam
  Int_t fNParticles;              //! number of particles in the input tree genTracks
  Int_t fDebug;                   //! debug flag  
  Int_t fNextTreeGenEntryToRead;    //! last entry already read from genTracks tree
  //
  AliTPCGenInfo*  fMCInfo;           //! MC information writen per particle
  AliTPCRecInfo*  fRecInfo;          //! Rec. information writen per particle
  //
  AliTPCtrack *fTPCTrack;         //! pointer to TPC track to connect branch

  ClassDef(TPCCmpTr,1)    // class which creates and fills tree with TPCGenTrack objects
};
ClassImp(TPCCmpTr)

