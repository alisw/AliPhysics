#ifndef ALIITS_H
#define ALIITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//           Manager class for set: ITS                               //
////////////////////////////////////////////////////////////////////////

#include <TObjArray.h> // used in inline function GetModule.
#include <TBranch.h>   // used in inline function SetHitsAddressBranch

#include "AliLoader.h"
#include "AliDetector.h"
#include "AliITSDetType.h"

class TString;
class TTree;
class TFile;

//class AliITSDetType;
class AliITSsimulation;
class AliITSClusterFinder;
class AliITSLoader;
class AliITSsegmentation;
class AliITSresponse;
class AliITShit;
class AliITSgeom;
class AliITSpListItem;
class AliITSdigit;
class AliITSRecPoint;
class AliITSRawCluster;
class AliITSmodule;
class AliVertexer;

const Int_t kNTYPES=3;

class AliITS : public AliDetector {

 public:
    //================= Standard Classes ===============================
    AliITS();  // Default creator.
    AliITS(const char *name, const char *title); // standard Creator
    virtual ~AliITS(); // destructor
    AliITS(const AliITS &source); // copy constructor. Not to be used!
    AliITS& operator=(AliITS &source); // = operator. Not to be used!
    virtual Int_t IsVersion() const {return 1;}
    virtual Int_t DistancetoPrimitive(Int_t px, Int_t py) const;

    //===================== Simulation Geometry ========================
    // get geometry version - detailed (major) or coarse (minor)
    virtual Int_t GetMajorVersion() const {return -1;}
    virtual Int_t GetMinorVersion() const {return -1;}
    virtual void  GetGeometryVersion(Int_t &a,Int_t &b) 
	                   {a = GetMajorVersion();b=GetMinorVersion();return;}
    virtual void  SetEUCLID(Bool_t euclid=kTRUE) {fEuclidOut = euclid;}
    //-------------------- Geometry Transformations --------------------
    // ITS geometry functions
    AliITSgeom   *GetITSgeom() const {return fITSgeom;}
    // Sets ITS geometry ! be very careful using this function.
    void   SetITSgeom(AliITSgeom *geom) {fITSgeom = geom;}
    // return pointer to the array of modules
    TObjArray    *GetModules() const {return fITSmodules;}
    // return pointer to a particular module
    AliITSmodule *GetModule(Int_t index) {return (AliITSmodule *)
					      (fITSmodules->At(index));}

    //================ Necessary general Classes =======================
    virtual void Init();
    virtual AliLoader* MakeLoader(const char* topfoldername);
    virtual void SetDefaults();
    virtual void SetDefaultSimulation();
    virtual void SetDefaultClusterFinders();
    virtual void MakeBranch(Option_t *opt=" ");
    virtual void SetTreeAddress();
    // For a given branch from the treeH sets the TClonesArray address.
    virtual void SetHitsAddressBranch(TBranch *b) {b->SetAddress(&fHits);}
    // Return pointer to DetType #id
    AliITSDetType *DetType(Int_t id);
    //Int_t           NDetTypes() {return fNDetTypes;}
    //---------- Configuration Methods (per detector type) -------------
    // Determines which ITS subdetectors will be processed. Effects
    // digitization, and Reconstruction only.
    void SetDetectors(Option_t *opt="All"){fOpt = opt;}
    // Returns the list of ITS subdetectors that will be processed.
    Option_t* GetDetectors(){return fOpt;}
    // Set response 
    virtual void SetResponseModel(Int_t id, AliITSresponse *response);
    // Set segmentation 
    virtual void SetSegmentationModel(Int_t id, AliITSsegmentation *seg);
    // Set simulation - temporary 
    virtual void SetSimulationModel(Int_t id, AliITSsimulation *sim);
    // Set simulation - temporary 
    virtual AliITSsimulation* GetSimulationModel(Int_t id){
	return ((AliITSDetType*)(fDetTypes->At(id)))->GetSimulationModel();}
    // Set reconstruction 
    virtual void SetReconstructionModel(Int_t id, AliITSClusterFinder *rec);
    // Set class names for digit and rec point 
    virtual void SetClasses(Int_t id, const char *digit, const char *cluster);

    //=================== Hits =========================================
    virtual void StepManager() {} // See Step Manager for specific geometry.
    virtual void AddHit(Int_t track, Int_t *vol, Float_t *hits);
    //------------ sort hits by module for Digitisation ----------------
    virtual void InitModules(Int_t size,Int_t &nmodules);  
    virtual void FillModules(TTree *treeH, Int_t mask = 0);
    virtual void FillModules(Int_t evnt,Int_t bgrev,Int_t nmodules,
			     Option_t *opt, const char *filename);
    virtual void ClearModules();

    //===================== Digitisation ===============================
    void MakeBranchS(const char *file);
    void SetTreeAddressS(TTree *treeS);
    TClonesArray * GetSDigits() { return fSDigits; }
    void MakeBranchInTreeD(TTree *treeD,const char *file=0);
    void MakeBranchD(const char *file){
	MakeBranchInTreeD(GetLoader()->TreeD(),file);}
    void SetTreeAddressD(TTree *treeD);
    void Hits2SDigits(); // Turn hits into SDigits
    void Hits2PreDigits(); // Turn hits into SDigits
    AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
    void SDigits2Digits(){SDigitsToDigits("All");} // Turn SDigits to Digits
    void SDigitsToDigits(Option_t *opt="All"); // Turn SDigits to Digits
    void Hits2Digits(); // Turn hits straight into Digits.
    //------------------ Internal functions ----------------------------
    // Standard Hits To SDigits function
    void HitsToSDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                 Option_t *add, Option_t *det, const char *filename);
    // Standard Hits To SDigits function
    void HitsToPreDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                 Option_t *add, Option_t *det, const char *filename);
    // Standard Hits To Digits function
    void HitsToDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                 Option_t *add, Option_t *det, const char *filename);
    void ResetSDigits();                  // Resets the Summable digits.
    void ResetDigits();                   // depending on how the
    void ResetDigits(Int_t branch);       // tree will be filled only
    void AddSumDigit(AliITSpListItem &sdig);
    void AddRealDigit(Int_t branch, Int_t *digits);
    void AddSimDigit(Int_t branch, AliITSdigit *d);
    void AddSimDigit(Int_t branch,Float_t phys,Int_t* digits,
		     Int_t* tracks,Int_t *hits,Float_t* trkcharges);
    // Return pointers to digits 
    TObjArray    *Dtype() {return fDtype;}
    Int_t        *Ndtype() {return fNdtype;}
    TClonesArray *DigitsAddress(Int_t id)
	{return ((TClonesArray *) (*fDtype)[id]);}
    void SelectVertexer(TString sel=" "){fSelectedVertexer = sel;}

    //===================== Raw Data IO ================================
    // Write digits into raw data format
    virtual void Digits2Raw();

    //==================== Clusterization ==============================
    // create separate tree for clusters - declustering refining
    void MakeTreeC(Option_t *option="C");
    void GetTreeC(Int_t event);
    void AddCluster(Int_t branch, AliITSRawCluster *c);
    void ResetClusters();                 // one of the methods in 
    void ResetClusters(Int_t branch);     // the pair will be kept
    // Return pointers to clusters 
    TObjArray    *Ctype() {return fCtype;}
    Int_t        *Nctype() {return fNctype;}
    TClonesArray *ClustersAddress(Int_t id) 
                   {return ((TClonesArray *) (*fCtype)[id]);}

    //=================== Reconstruction ===============================
    void MakeBranchR(const char *file, Option_t *opt=" ");
    void MakeBranchRF(const char *file){MakeBranchR(file,"Fast");}
    void MakeBranchC();
    void SetTreeAddressR(TTree *treeR);
    void AddRecPoint(const AliITSRecPoint &p);
    void HitsToFastRecPoints(Int_t evNumber,Int_t bgrev,Int_t size,
                 Option_t *add, Option_t *det, const char *filename);
    void Digits2Reco();
    void DigitsToRecPoints(Int_t evNumber,Int_t lastEntry,Option_t *det);
    void ResetRecPoints();
    // Return pointer to rec points 
    TClonesArray  *RecPoints()   {return fRecPoints;}

    virtual void Reconstruct() const;
    virtual AliTracker* CreateTracker() const;
    virtual AliVertexer* CreateVertexer() const;
    virtual void FillESD(AliESD* esd) const;
     
 protected:
    //================== Data Members ==================================
    AliITSgeom   *fITSgeom;    // Pointer to ITS geometry
    Bool_t        fEuclidOut;  // Flag to write geometry in euclid format
    TObjArray    *fITSmodules; //! Pointer to ITS modules
    Option_t     *fOpt;        //! Detector option ="All" unless changed.

    Int_t         fIdN;        // the number of layers
    Int_t        *fIdSens;     //[fIdN] layer identifier
    TString      *fIdName;     //[fIdN] layer identifier

    Int_t         fNDetTypes;  // Number of detector types
    TObjArray    *fDetTypes;   // List of detector types

    TClonesArray  *fSDigits;    //! List of Summable digits.
    Int_t         fNSDigits;   // Number of Summable Digits.

    TObjArray    *fDtype;      //! List of digits
    Int_t        *fNdtype;     //[fNDetTypes] Num. of digits per type of det. 

    TObjArray    *fCtype;      //! List of clusters
    Int_t        *fNctype;     //[fNDetTypes] Num. of clust. per type of det.

    TClonesArray *fRecPoints;  //! List of reconstructed points
    Int_t         fNRecPoints; // Number of rec points
    TString fSelectedVertexer; // Vertexer selected in CreateVertexer

    ClassDef(AliITS,4) // Base class for ITS
};

#endif
