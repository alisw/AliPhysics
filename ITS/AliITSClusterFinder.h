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
class AliITS;
class AliITSdigit;
class AliITSRecPoint;

//----------------------------------------------------------------------
class AliITSClusterFinder :public TObject{
  public:
    AliITSClusterFinder(); // Default constructor
    // Standard Constructor
    AliITSClusterFinder(AliITSsegmentation *seg, AliITSresponse *resp);
    AliITSClusterFinder(AliITSsegmentation *seg, AliITSresponse *resp,
                        TClonesArray *digits);// Standard+ Constructor
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
    //
    // Setters and Getters
    // segmentation
    virtual void SetSegmentation(AliITSsegmentation *segmentation) {
        fSegmentation=segmentation;}
    //Returns fSegmentation
    virtual AliITSsegmentation* GetSegmentation()const{return fSegmentation;}
    // Digit
    virtual void SetDigits(TClonesArray *itsDigits) {// set digits
        fDigits=itsDigits;fNdigits = fDigits->GetEntriesFast();}
    virtual AliITSdigit* GetDigit(Int_t i){ // Returns ith digit
        return (AliITSdigit*) fDigits->UncheckedAt(i);}
    virtual TClonesArray* Digits(){return fDigits;}// Gets fDigits
    virtual Int_t   NDigits() const {return fNdigits;}// Get Number of Digits
    // Response
    //Return Response
    virtual AliITSresponse* GetResponse()const{return fResponse;}
    virtual void SetResponse(AliITSresponse *response) {// set response
        fResponse=response;}
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
    virtual void FindCluster(Int_t,Int_t,AliITSRawCluster *) {}// find cluster
    //
    virtual void SetModule(Int_t module){fModule = module;}// Set module number
    virtual Int_t GetModule()const{return fModule;}// Returns module number
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
    void Print(ostream *os); // Class ascii print function
    void Read(istream *os);  // Class ascii read function

    // Conversion from RecPoints to V2Clusters
    void RecPoints2Clusters(const TClonesArray *points, Int_t idx, TClonesArray *clusters);

 
  protected:
    // methods 
    static void CheckLabels(Int_t lab[3]);
    void Init();                      
    AliITSClusterFinder(const AliITSClusterFinder &source); // copy constructor
    // assignment operator
    AliITSClusterFinder& operator=(const AliITSClusterFinder &source);
   // data members       
    Int_t              fDebug;         //! Debug flag/level
    Int_t              fModule;        //! Module number to be reconstuctted
    TClonesArray       *fDigits;       //! digits
    Int_t              fNdigits;       //! num of digits

    AliITSresponse     *fResponse;     //! response
    AliITSsegmentation *fSegmentation; //! segmentation
    TClonesArray       *fClusters;     //! Array of clusters
    Int_t              fNRawClusters;  //! in case we split the cluster
                                       // and want to keep track of 
                                       // the cluster which was splitted
    AliITSMap          *fMap;          //! map
    AliITS             *fITS;          //! pointer to the ITS      
    Int_t              fNperMax;       //! NperMax
    Int_t              fDeclusterFlag; //! DeclusterFlag
    Int_t              fClusterSize;   //! ClusterSize
    Int_t              fNPeaks;        //! NPeaks  
    // Data members needed to fill AliCluster objects
    Float_t fYshift[2200];       // y-shifts of detector local coor. systems 
    Float_t fZshift[2200];       // z-shifts of detector local coor. systems 
    Int_t fNdet[2200];           // detector index  
    Int_t fNlayer[2200];         // detector layer

    ClassDef(AliITSClusterFinder,4) //Class for clustering and reconstruction of space points
};
// Input and output functions for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSClusterFinder &source);
istream &operator>>(istream &os,AliITSClusterFinder &source);
#endif
