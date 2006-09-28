
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
  void UpdatePoints(AliESDtrack* track);
  void Update(AliMCInfo* info,AliTPCParam * par, Bool_t reconstructed, AliESD *event);
  void Reset();
  Float_t  fTPCPoints[10]; //start , biggest end points,max density .. density at the last 30 pad-rows
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
  AliITStrackMI fITStrack;        //its track
  AliTRDtrack fTRDtrack;        //its track
  Float_t fBestTOFmatch;        //best matching between times
  Float_t fTRLocalCoord[3];       //local coordinates of the track ref.
  Int_t   fReconstructed;         //flag if track was reconstructed
  Int_t fFake;             // fake track
  Int_t fMultiple;         // number of reconstructions
  Bool_t fTPCOn;           // TPC refitted inward
  Int_t  fStatus[4];        // status -0 not found - 1 -only in - 2 -in-out -3 -in -out-refit
  Bool_t fITSOn;           // ITS refitted inward
  Bool_t fTRDOn;           // ITS refitted inward
  Float_t fDeltaP;          //delta of momenta
  Double_t fSign;           // sign
  Int_t fLabels[2];         // labels
  ClassDef(AliESDRecInfo,2)  // container for 
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
  Double_t       fRs[2];     // minimum radius in rphi intersection
  Double_t       fDistMinR; // distance at minimal radius
  Double_t       fPm[3];    //momentum at the vertex mother
  Double_t       fAngle[3]; //three angles
  Double_t       fRr;       // rec position of the vertex 
  Int_t          fLab[2];   //MC label of the partecle
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  Int_t          fV0Status;       // status of the kink
  AliESDV0       fV0tpc;           // Vo information from reconsturction according TPC
  AliESDV0       fV0its;           // Vo information from reconsturction according ITS
  AliESDV0       fV0rec;           // V0 information form the reconstruction
  Int_t          fMultiple;
  Int_t          fV0Multiple;
  Int_t          fRecStatus;    // status form the reconstuction
  ClassDef(AliESDRecV0Info,2)   // container for  
};

ClassImp(AliESDRecV0Info)


class AliESDRecKinkInfo: public TObject {
public:
  void Update();
  AliESDRecInfo  fT1;      //track1
  AliESDRecInfo  fT2;      //track2  
  AliESDkink     fKink;    //kink
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
  Double_t       fMinR;     // minimum radius in rphi intersection
  Double_t       fDistMinR; // distance at minimal radius
  Int_t          fLab[2];   //MC label of the partecle
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  Int_t          fStatus;       //status -tracks 
  Int_t          fRecStatus;    //kink -status- 0 - not found  1-good -  fake
  Int_t          fMultiple;
  Int_t          fKinkMultiple;
  ClassDef(AliESDRecKinkInfo,1)   // container for  
};

ClassImp(AliESDRecKinkInfo)



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
  Int_t BuildKinkInfo0(Int_t eventNr); // build kink info 0
  Int_t BuildV0Info(Int_t eventNr); // build kink info 0
  void  MakePoints(AliESDtrack * track, AliPointsMI &points);
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
  TTree *fTreeCmpKinks;                //! output tree with cmp Kinks
  TTree *fTreeCmpV0;                //! output tree with cmp V0
  //
  char  fFnGenTracks[1000];             //! input file name with gen tracks
  TFile *fFileGenTracks;
  TTree *fTreeGenTracks;
  TTree *fTreeGenKinks;            // tree with gen kinks
  TTree *fTreeGenV0;            // tree with gen V0
  //
  //
  Int_t  fDirection;
  //
  AliRunLoader * fLoader;         //! pointer to the run loader
  //TTree *fTreeRecTracks;          //! tree with reconstructed tracks
  //
  Short_t *fIndexRecTracks;         //! index of particle label in the TreeT_ESD
  Short_t *fFakeRecTracks;          //! number of fake tracks
  Short_t *fMultiRecTracks;         //! number of multiple reconstructions
  //
  Short_t *fIndexRecKinks;         //! index of particle label in treeesd
  Short_t *fMultiRecKinks;         //! number of multiple reconstructions
  Short_t *fSignedKinks;           //! indicator that kink was not fake
  //
  Short_t *fIndexRecV0;         //! index of particle label in treeesd
  Short_t *fMultiRecV0;         //! number of multiple reconstructions
  Short_t *fSignedV0;                //! indicator that kink was not fake
  //
  TObjArray *fRecArray;           // container with rec infos
  AliESD *fEvent;                 //!event
  //
  AliTPCParam* fParamTPC;         //! AliTPCParam
  Int_t fNParticles;              //! number of particles in the input tree genTracks
  Int_t fDebug;                   //! debug flag  
  Int_t fNextTreeGenEntryToRead;    //! last entry already read from genTracks tree
  Int_t fNextKinkToRead;            //! last entry already read from genKinks tree
  Int_t fNextV0ToRead;            //! last entry already read from genV0 tree
  //
  AliMCInfo*  fMCInfo;           //! MC information writen per particle
  AliGenKinkInfo* fGenKinkInfo;      //! MC information writen per Kink
  AliGenV0Info* fGenV0Info;      //! MC information writen per Kink
  AliESDRecInfo*  fRecInfo;          //! Rec. information writen per particle
  AliESDRecKinkInfo* fRecKinkInfo;    //! reconstructed kink info
  AliESDRecV0Info* fRecV0Info;    //! reconstructed kink info
  //

  ClassDef(ESDCmpTr,1)    // class which creates and fills tree with ESDGenTrack objects
};
ClassImp(ESDCmpTr)



