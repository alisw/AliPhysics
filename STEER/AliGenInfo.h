


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
// Start of implementation of the class AliMCInfo
//
////////////////////////////////////////////////////////////////////////

class AliMCInfo: public TObject {

public:
  AliMCInfo();
  ~AliMCInfo();
  void Update();

  AliTrackReference  fTrackRef;      // track reference saved in the output tree
  AliTrackReference  fTrackRefOut;   // decay track reference saved in the output tree
  AliTrackReference  fTRdecay;       // track reference at decay point
  //
  Int_t fPrimPart;                   // index of primary particle in TreeH
  TParticle fParticle;           // generated particle 
  Float_t   fMass;               // mass of the particle
  Float_t fCharge;               //
  Int_t fLabel;                   // track label
  Int_t fEventNr;                 // event number
  Int_t fMCtracks;                // indication of how many times the track is retuturned back
  Int_t fPdg;                     //pdg code
  Float_t fDecayCoord[3];         // position of particle decay
  Double_t fVDist[4];             //distance of the particle vertex from primary vertex
  Bool_t fTPCdecay;               //indicates decay in TPC
  Int_t fRowsWithDigitsInn;    // number of rows with digits in the inner sectors
  Int_t fRowsWithDigits;       // number of rows with digits in the outer sectors
  Int_t fRowsTrackLength;      // last - first row with digit
  Float_t fPrim;               // theoretical dedx in tpc according particle momenta and mass
  digitRow fTPCRow;               // information about digits row pattern
  Int_t fNTPCRef;                 // tpc references counter
  Int_t fNITSRef;                 // ITS references counter
  Int_t fNTRDRef;                  // TRD references counter
  Int_t fNTOFRef;                  // TOF references counter
  TClonesArray * fTPCReferences;  //containner with all track references -in the TPC
  TClonesArray * fITSReferences;  //container with ITS references
  TClonesArray * fTRDReferences;  //container with TRD references  
  TClonesArray * fTOFReferences;  //container with TRD references  
  //
  ClassDef(AliMCInfo,3)  // container for 
};
ClassImp(AliMCInfo)



class AliGenV0Info: public TObject {
public:
  AliMCInfo fMCd;      //info about daughter particle - second particle for V0
  AliMCInfo fMCm;      //info about mother particle   - first particle for V0
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
ClassImp(AliGenV0Info)




class AliGenKinkInfo: public TObject {
public:
  AliMCInfo fMCd;      //info about daughter particle - second particle for V0
  AliMCInfo fMCm;      //info about mother particle   - first particle for V0
  void Update();        // put some derived info to special field 
  Float_t GetQt();      //
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
  ClassDef(AliGenKinkInfo,1)  // container for  
};
ClassImp(AliGenKinkInfo)




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
  void Reset();
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
  //void  FillInfo(Int_t iParticle);
  void SetFirstEventNr(Int_t i) {fFirstEventNr = i;}
  void SetNEvents(Int_t i) {fNEvents = i;}
  void SetDebug(Int_t level) {fDebug = level;}
  Int_t SetIO(Int_t eventNr);
  Int_t CloseIOEvent();
  Int_t CloseIO();
  Int_t SetIO();
  Float_t TR2LocalX(AliTrackReference *trackRef,
		    AliTPCParam *paramTPC);
  AliMCInfo * GetInfo(UInt_t i){return (i<fNParticles)? fGenInfo[i]:0;}
  AliMCInfo * MakeInfo(UInt_t i);

public:
  Int_t fDebug;                   //! debug flag  
  Int_t fEventNr;                 //! current event number
  Int_t fLabel;                   //! track label
  Int_t fNEvents;                 //! number of events to process
  Int_t fFirstEventNr;            //! first event to process
  UInt_t fNParticles;              //! number of particles in TreeK
  TTree * fTreeGenTracks;          //! output tree with generated tracks
  TTree * fTreeKinks;             //!  output tree with Kinks
  TTree * fTreeV0;                //!  output tree with V0
  TTree * fTreeHitLines;          //! tree with hit lines
  char  fFnRes[1000];             //! output file name with stored tracks
  TFile *fFileGenTracks;          //! output file with stored fTreeGenTracks
  //
  AliRunLoader * fLoader;         //! pointer to the run loader
  TTree * fTreeD;                 //! current tree with digits
  TTree * fTreeTR;                //! current tree with TR
  AliStack *fStack;               //! current stack
  // 
  AliMCInfo **   fGenInfo;    //! array with pointers to gen info
  Int_t fNInfos;                  //! number of tracks with infos
  //
  AliTPCParam* fParamTPC;         //! AliTPCParam
  Float_t fVPrim[3];             //! primary vertex position
                                  // the fVDist[3] contains size of the 3-vector

private:

// some constants for the original non-pareller tracking (by Y.Belikov)
  static const Int_t seedRow11 = 158;  // nRowUp - 1
  static const Int_t seedRow12 = 139;  // nRowUp - 1 - (Int_t) 0.125*nRowUp
  static const Int_t seedRow21 = 149;  // seedRow11 - shift
  static const Int_t seedRow22 = 130;  // seedRow12 - shift
  static const Double_t kRaddeg = 180./3.14159265358979312;
  // 
  static const Double_t fgTPCPtCut = 0.03; // do not store particles with generated pT less than this
  static const Double_t fgITSPtCut = 0.2; // do not store particles with generated pT less than this
  static const Double_t fgTRDPtCut = 0.2; // do not store particles with generated pT less than this
  static const Double_t fgTOFPtCut = 0.15; // do not store particles with generated pT less than this
 
  ClassDef(AliGenInfoMaker,1)    // class which creates and fills tree with TPCGenTrack objects
};
ClassImp(AliGenInfoMaker)



class AliComparisonDraw: public TObject{
public:
  AliComparisonDraw(){fPoints=0; fView=0;}
  void InitView();
  TH1F * DrawXY(const char * chx, const char *chy, const char* selection, 
		const char * quality,Int_t nbins, Float_t minx, Float_t maxx, 
		Float_t miny, Float_t maxy, Int_t nBinsRes=30);
  TH1F * DrawLogXY(const char * chx, const char *chy, const char* selection, 
		   const char * quality, Int_t nbins,Float_t minx, Float_t maxx, 
		   Float_t miny, Float_t maxy, Int_t nBinsRes=30); 
  TH1F * Eff(const char *variable, const char* selection, const char * quality, 
	     Int_t nbins,Float_t min, Float_t max); 
  TH1F * EffLog(const char *variable, const char* selection, const char * quality, 
	     Int_t nbins,Float_t min, Float_t max);
  //
  static void   AliLabelAxes(TH1* histo, const char* xAxisTitle, const char* yAxisTitle);
  static Double_t* CreateLogBins(Int_t nBins, Double_t xMin, Double_t xMax);
  //
  static TH1F*  CreateEffHisto(TH1F* hGen, TH1F* hRec);
  static TH1F*  CreateResHisto(TH2F* hRes2, TH1F **phMean, 
				Bool_t drawBinFits = kTRUE,Bool_t overflowBinFits = kFALSE);
  void   DrawFriend2D(const char * chx, const char *chy, const char* selection, TTree * tfriend);
  void   GetPoints3D(const char * label, const char * chpoints, const char* selection, TTree * tpoints, Int_t color=6, Float_t rmin=4.);
  void   Draw3D(Int_t min=0, Int_t max = 10000);
  void SavePoints(const char* name);
 public: 
  TTree * fTree;
  TH1F  * fRes;  //temporary file
  TH1F  * fMean;  //temporary file
  TView * fView;  //3D view
  TCanvas *fCanvas; //canvas
  TObjArray *fPoints;
  ClassDef(AliComparisonDraw,1)
};
ClassImp(AliComparisonDraw)


class AliPointsMI: public TObject{
 public:
  AliPointsMI();
  AliPointsMI(Int_t n, Float_t *x,Float_t *y, Float_t *z);
  void Reset();
  void Reset(AliDetector * det, Int_t particle);  //load points for given particle
  ~AliPointsMI();
  Int_t   fN;  //number of points;
  Float_t *fX; //[fN] pointer to x
  Float_t *fY; //[fN] pointer to y
  Float_t *fZ; //[fN] pointer to Z  
  Int_t   fCapacity; //!allocated size of the x,y,x
  Int_t   fLabel0; //label
  Int_t   fLabel1; //label
  ClassDef(AliPointsMI,1)
};
ClassImp(AliPointsMI)


AliTPCParam * GetTPCParam();

