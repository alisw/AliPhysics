
class AliESDComparisonDraw: public AliComparisonDraw{
public:
  AliESDComparisonDraw(){fTree = 0;}
  void SetIO(const char *fname = "cmpESDTracks.root");  
  ClassDef(AliESDComparisonDraw,1)
};
ClassImp(AliESDComparisonDraw)


/////////////////////////////////////////////////////////////////////////
class AliESDRecInfo: public TObject {

public:
  AliESDRecInfo(){}
  ~AliESDRecInfo(){}
  //
  void Update(AliMCInfo* info,AliTPCParam * par, Bool_t reconstructed);
  Double_t fTPCinR0[5];   //generated position of the track at inner tpc - radius [3] and fi [4]
  Double_t fTPCinR1[5];   //reconstructed postion of the track           - radius [3] and fi [
  Double_t fTPCinP0[5];   //generated position of the track at inner tpc
  Double_t fTPCinP1[5];   //reconstructed postion of the track
  Double_t fTPCAngle0[2]; // generated angle 
  Double_t fTPCAngle1[2]; //refconstructed angle 
  Double_t fTPCDelta[5];  // deltas
  Double_t fTPCPools[5];  // pools
  Double_t fITSinR0[5];   //generated position of the track at inner tpc
  Double_t fITSinR1[5];   //reconstructed postion of the track
  Double_t fITSinP0[5];   //generated position of the track at inner tpc
  Double_t fITSinP1[5];   //reconstructed postion of the track
  Double_t fITSAngle0[2]; // generated angle 
  Double_t fITSAngle1[2]; //refconstructed angle
  Double_t fITSDelta[5];  // deltas
  Double_t fITSPools[5];  // pools
  AliESDtrack fESDTrack;          // tpc track
  AliITStrackV2 fITStrack;        //its track
  AliTRDtrack fTRDtrack;        //its track
  Float_t fTRLocalCoord[3];       //local coordinates of the track ref.
  Int_t   fReconstructed;         //flag if track was reconstructed
  Int_t fFake;             // fake track
  Int_t fMultiple;         // number of reconstructions
  Bool_t fTPCOn;           // TPC refitted inward
  Bool_t fITSOn;           // ITS refitted inward
  Bool_t fTRDOn;           // ITS refitted inward
  Float_t fDeltaP;          //delta of momenta
  void Reset();
  ClassDef(AliESDRecInfo,1)  // container for 
};
ClassImp(AliESDRecInfo)

void AliESDRecInfo::Reset()
{
  fMultiple =0; 
  fFake     =0;
  fReconstructed=0;
}


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


class AliESDRecV0Info: public TObject {
public:
  void Update(Float_t vertex[3]);
  AliESDRecInfo  fT1;      //track1
  AliESDRecInfo  fT2;      //track2  
  Double_t       fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t       fDist2;    //info about closest distance parabolic DCA
  Double_t       fInvMass;  //reconstructed invariant mass -
  //
  Double_t       fPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t       fXr[3];     //rec. position according helix
  //
  Double_t       fPm[3];    //momentum at the vertex mother
  Double_t       fAngle[3]; //three angles
  Double_t       fRr;       // rec position of the vertex 
  Int_t          fLab[2];   //MC label of the partecle
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  ClassDef(AliESDRecV0Info,1)   // container for  
};

ClassImp(AliESDRecV0Info)



////////////////////////////////////////////////////////////////////////
// 
// Start of implementation of the class ESDCmpTr
//
////////////////////////////////////////////////////////////////////////

class ESDCmpTr {

public:
  ESDCmpTr();
  ESDCmpTr(const char* fnGenTracks,
	   const char* fnCmpRes      ="cmpTracks.root", 
	   const char* fnGalice      ="galice.root", Int_t direction=0,
	   Int_t nEvents=1, Int_t firstEvent=0);
  virtual ~ESDCmpTr();
  void Reset();
  Int_t Exec();
  Int_t Exec(Int_t nEvents, Int_t firstEventNr);
  Int_t SetIO();
  Int_t SetIO(Int_t eventNr );
  void CreateTreeCmp();
  void CloseOutputFile();
  Bool_t ConnectGenTree();
  Int_t TreeGenLoop(Int_t eventNr);
  Int_t TreeTLoop();
  void SetFirstEventNr(Int_t i) {fFirstEventNr = i;}
  void SetNEvents(Int_t i) {fNEvents = i;}
  void SetDebug(Int_t level) {fDebug = level;}

// tmp method, should go to TrackReferenceESD
  static TVector3 TR2Local(AliTrackReference *trackRef,
		    AliTPCParam *paramTPC);

private:

  Int_t fEventNr;                 //! current event number
  Int_t fNEvents;                 //! number of events to process
  Int_t fFirstEventNr;            //! first event to process
  //
  char  fFnCmp[1000];                   //! output file name with cmp tracks
  TFile *fFileCmp;                //! output file with cmp tracks
  TTree *fTreeCmp;                //! output tree with cmp tracks
  //
  char  fFnGenTracks[1000];             //! input file name with gen tracks
  TFile *fFileGenTracks;
  TTree *fTreeGenTracks;
  //
  //
  Int_t  fDirection;
  //
  AliRunLoader * fLoader;         //! pointer to the run loader
  //TTree *fTreeRecTracks;          //! tree with reconstructed tracks
  //
  Int_t *fIndexRecTracks;         //! index of particle label in the TreeT_ESD
  Int_t *fFakeRecTracks;          //! number of fake tracks
  Int_t *fMultiRecTracks;         //! number of multiple reconstructions
  //
  TObjArray * fTracks;            //!container with tracks 
  AliESD *fEvent;                 //!event

  //
  AliTPCParam* fParamTPC;         //! AliTPCParam
  Int_t fNParticles;              //! number of particles in the input tree genTracks
  Int_t fDebug;                   //! debug flag  
  Int_t fNextTreeGenEntryToRead;    //! last entry already read from genTracks tree
  //
  AliMCInfo*  fMCInfo;           //! MC information writen per particle
  AliESDRecInfo*  fRecInfo;          //! Rec. information writen per particle
  //

  ClassDef(ESDCmpTr,1)    // class which creates and fills tree with ESDGenTrack objects
};
ClassImp(ESDCmpTr)



