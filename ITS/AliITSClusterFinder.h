#ifndef ALIITSCLUSTERFINDER_H
#define ALIITSCLUSTERFINDER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id$ */

////////////////////////////////////////////////
//  ITS Cluster Finder Class                  //
//                                            //
//                                            //
////////////////////////////////////////////////

#include <TObject.h>
#include <TClonesArray.h>

class AliITSMap;
class AliITSresponse;
class AliITSsegmentation;
class AliITSRawCluster;
class AliITSdigit;
class AliITSRecPoint;
class AliITSDetTypeRec;
class AliRawReader;

//----------------------------------------------------------------------
class AliITSClusterFinder :public TObject{
  public:
    AliITSClusterFinder(); // Default constructor
    // Standard Constructor
    AliITSClusterFinder(AliITSDetTypeRec* dettyp);
    AliITSClusterFinder(AliITSDetTypeRec* dettyp,TClonesArray *digits);// Standard+ Constructor
    virtual ~AliITSClusterFinder(); // Destructor
    //
    // Do the Reconstruction.
    virtual void FindRawClusters(Int_t mod=0); // Finds cluster of digits.
     //
    // Sets the debug flag for debugging output
    void SetDebug(Int_t level=1){fDebug=level;}
    // Clears the debug flag so no debugging output will be generated
    void SetNoDebug(){fDebug=0;}
    // Returns the debug flag value
    Bool_t GetDebug(Int_t level=1)const {return fDebug>=level;}
    // Digit
    virtual void SetDigits(TClonesArray *itsDigits) {// set digits
        fDigits=itsDigits;fNdigits = fDigits->GetEntriesFast();}
    virtual AliITSdigit* GetDigit(Int_t i){ // Returns ith digit
        return (AliITSdigit*) fDigits->UncheckedAt(i);}
    virtual TClonesArray* Digits(){return fDigits;}// Gets fDigits
    virtual Int_t   NDigits() const {return fNdigits;}// Get Number of Digits
    // clulsters
    // Set fClusters up
    virtual void SetClusters(TClonesArray *itsClusters){// set clusters
        fClusters = itsClusters;fNRawClusters = fClusters->GetEntriesFast();}
    // Get fCluters
    virtual TClonesArray* Clusters(){return fClusters;}
    // Get fCluter
    virtual AliITSRawCluster* Cluster(Int_t i){
        return (AliITSRawCluster*)(fClusters->At(i));}
    // Returns the present number of enteries
    virtual Int_t NClusters()const {return fClusters->GetEntriesFast();}
    // returns fNRawClusters
    virtual Int_t GetNRawClusters() const {return fNRawClusters;}
    // Determins if digit i has a neighbor and if so that neighor index is j.
    virtual void AddCluster(Int_t branch,AliITSRawCluster *c);
    virtual void AddCluster(Int_t branch,AliITSRawCluster *c,
                            AliITSRecPoint &rp);
    virtual void   FillCluster(AliITSRawCluster *,Int_t) {}// fiil cluster
    virtual void   FillCluster(AliITSRawCluster *cluster) {// fill cluster
        FillCluster(cluster,1);}

    virtual void SetModule(Int_t module){fModule = module;}// Set module number
    virtual Int_t GetModule()const {return fModule;}// Returns module number

  void SetEvent(Int_t event) { fEvent=event; }
  virtual void RawdataToClusters(AliRawReader* /*rawReader*/,TClonesArray** /*clusters*/) {
    Warning("RawdataToClusters","Method not implemented in this class ");}

    //
    // RecPoints
    // Given a cluster of digits, creates the nessesary RecPoint. May also
    // do some peak separation.
    virtual void CreateRecPoints(TObjArray *,Int_t){};
    // Others
    virtual void  SetMap(AliITSMap *m) {fMap=m;}// map
    AliITSMap* Map(){return fMap;}// map
    virtual Int_t GetNperMax() const {return fNperMax;}// returns fNperMax
    // returns fDeclusterFlag
    virtual Int_t GetDeclusterFlag()const{return fDeclusterFlag;} 
    // returns fClusterSize
    virtual Int_t GetClusterSize() const {return fClusterSize;} 
    virtual Int_t GetNPeaks() const {return fNPeaks;}// returns fNPeaks
    //
    virtual Bool_t IsNeighbor(TObjArray *digs,Int_t i,Int_t j[]) const;
    virtual void Decluster(AliITSRawCluster *) {}// Decluster
    // Set max. Number of cells per local cluster
    virtual void SetNperMax(Int_t npermax=3) {fNperMax = npermax;}
    //Decluster
    virtual void SetDeclusterFlag(Int_t flag=1){fDeclusterFlag=flag;}
        // Set max. cluster size ; bigger clusters will be rejected
    virtual void SetClusterSize(Int_t clsize=3) {fClusterSize = clsize;}
    virtual void CalibrateCOG() {}// Self Calibration of COG 
    virtual void CorrectCOG(){}// correct COG
    virtual Bool_t Centered(AliITSRawCluster *) const {return kTRUE;}// cluster
    //split by local maxima
    virtual void SplitByLocalMaxima(AliITSRawCluster *){}
    // IO functions
    void Print(ostream *os) const; // Class ascii print function
    void Read(istream *os);  // Class ascii read function
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual Int_t Read(const char *name) {return TObject::Read(name);}

    virtual void SetDetTypeRec(AliITSDetTypeRec* dtr) {fDetTypeRec=dtr;}
    AliITSDetTypeRec* GetDetTypeRec() const {return fDetTypeRec;}

    void InitGeometry(); 
 
  protected:
  class Ali1Dcluster {
  public:
    void SetY(Float_t y) {fY=y;}
    void SetQ(Float_t q) {fQ=q;}
    void SetNd(Int_t n)  {fNd=n;}
    void SetLabels(Int_t *lab) {fLab[0]=lab[0];fLab[1]=lab[1];fLab[2]=lab[2];}
    Float_t GetY() const {return fY;}
    Float_t GetQ() const {return fQ;}
    Int_t GetNd()const {return fNd;}
    Int_t GetLabel(Int_t lab) const { return fLab[lab]; }
  protected:
    Float_t fY; //cluster position
    Float_t fQ; //cluster charge
    Int_t fNd;  //number of digits
    Int_t fLab[3]; //track label
  };
  class AliBin {
  public:
    AliBin():fIndex(0),fMask(0xFFFFFFFE),fQ(0){}
    void SetIndex(UInt_t idx) {fIndex=idx;}
    void SetQ(UShort_t q)  {fQ=q;}
    void SetMask(UInt_t m) {fMask=m;}
    void Reset() {fIndex=0; fMask=0xFFFFFFFE; fQ=0;}

    void Use() {fMask&=0xFFFFFFFE;}
    Bool_t IsNotUsed() const {return (fMask&1);}
    Bool_t IsUsed() const {return !(IsNotUsed());}

    UInt_t   GetIndex() const {return fIndex;}
    UShort_t GetQ()     const {return fQ;}
    UInt_t   GetMask()  const {return fMask;}
  protected:
    UInt_t fIndex; //digit index
    UInt_t fMask; //peak mask
    UShort_t fQ;  //signal
  };
  void MakeCluster(Int_t k,Int_t max,AliBin *bins,UInt_t m,AliITSRecPoint &c);
  static Bool_t IsMaximum(Int_t k, Int_t max, const AliBin *bins);
  static void FindPeaks(Int_t k,Int_t m,AliBin*b,Int_t*idx,UInt_t*msk,Int_t&n);
  static void MarkPeak(Int_t k, Int_t max, AliBin *bins, UInt_t m);
  static void FindCluster(Int_t k,Int_t maxz,AliBin *bins,Int_t &n,Int_t *idx);

  static void CheckLabels2(Int_t lab[10]);
  static void AddLabel(Int_t lab[10], Int_t label);      
   // data members       

   Int_t              fDebug;         //! Debug flag/level
   Int_t              fModule;        //! Module number to be reconstuctted
   TClonesArray       *fDigits;       //! digits 
   Int_t              fNdigits;       //! num of digits 
 
   AliITSDetTypeRec* fDetTypeRec; //ITS object for reconstruction
   TClonesArray       *fClusters;     //! Array of clusters
   Int_t              fNRawClusters;  //! in case we split the cluster
   // and want to keep track of 
   // the cluster which was splitted
    AliITSMap          *fMap;          //! map
    Int_t              fNperMax;       //! NperMax
    Int_t              fDeclusterFlag; //! DeclusterFlag
    Int_t              fClusterSize;   //! ClusterSize
    Int_t              fNPeaks;        //! NPeaks  
    // Data members needed to fill AliCluster objects
    Int_t fNdet[2200];           // detector index  
    Int_t fNlayer[2200];         // detector layer

   Int_t fNModules;             // total number of modules    
   Int_t fEvent;                //event number
   
   AliITSClusterFinder(const AliITSClusterFinder &source); // copy constructor
    // assignment operator
    AliITSClusterFinder& operator=(const AliITSClusterFinder &source);
    

    ClassDef(AliITSClusterFinder,8) //Class for clustering and reconstruction of space points
};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSClusterFinder &source);
istream &operator>>(istream &os,AliITSClusterFinder &source);
#endif
