#ifndef ALIITS_H
#define ALIITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//           Manager class for set: ITS                               //
////////////////////////////////////////////////////////////////////////

//#define NEWVERSION

#include <TObjArray.h> // used in inline function GetModule.
#include <TBranch.h>   // used in inline function SetHitsAddressBranch

#include "AliRunLoader.h"
#include "AliDetector.h"
#ifndef NEWVERSION
#include "AliITSDetType.h"
#endif
#ifdef NEWVERSION
#include "AliITSDetTypeSim.h"
#include "AliITSDetTypeRec.h"
#endif

class TString;
class TTree;
class TFile;

class AliITSsimulation;
class AliITSClusterFinder;
class AliITSclusterV2;
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
class AliDigitizer;
class AliRunDigitizer;
class AliRawReader;

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
    virtual Int_t DistancetoPrimitive(Int_t,Int_t) const{return 999;};

    //===================== Simulation Geometry ========================
    // get geometry version - detailed (major) or coarse (minor)
    virtual Int_t GetMajorVersion() const {return -1;}
    virtual Int_t GetMinorVersion() const {return -1;}
    virtual void  GetGeometryVersion(Int_t &a,Int_t &b) 
	                   {a = GetMajorVersion();b=GetMinorVersion();return;}
    virtual void  SetEUCLID(Bool_t euclid=kTRUE) {fEuclidOut = euclid;}
    virtual Bool_t GetEUCLID()const {return fEuclidOut;}
    //-------------------- Geometry Transformations --------------------
#ifndef NEWVERSION
    // ITS geometry functions
    AliITSgeom   *GetITSgeom() const {return fITSgeom;}
    // Sets ITS geometry ! be very careful using this function.
    void   SetITSgeom(AliITSgeom *geom) {fITSgeom = geom;}
    // return pointer to the array of modules
    TObjArray    *GetModules() const {return fITSmodules;}
    // return pointer to a particular module
    AliITSmodule *GetModule(Int_t index) {return (AliITSmodule *)
					      (fITSmodules->At(index));}
#endif
#ifdef NEWVERSION
    // ITS geometry functions
    AliITSgeom* GetITSgeom()const{if(fDetTypeSim!=0)
        return fDetTypeSim->GetITSgeom();else if(fDetTypeRec!=0) 
            return fDetTypeRec->GetITSgeom();else return 0;}
    // ITS geometry functions From Simulation
    AliITSgeom* GetITSgeomSim()const{if(fDetTypeSim!=0)
        return fDetTypeSim->GetITSgeom();else return 0;}
    // ITS geometry functions From Reconstruction
    AliITSgeom* GetITSgeomRec()const{if(fDetTypeRec!=0)
        return fDetTypeRec->GetITSgeom();else return 0;}
    // Sets ITS geometry ! be very careful using this function.
    void   SetITSgeom(AliITSgeom *geom) {if(fDetTypeSim!=0)
        fDetTypeSim->SetITSgeom(geom);if(fDetTypeRec!=0) 
            fDetTypeRec->SetITSgeom(geom);}
    // Sets ITS geometry For Simulation ! be very careful using this function.
    void   SetITSgeomSim(AliITSgeom *geom) {if(fDetTypeSim!=0)
        fDetTypeSim->SetITSgeom(geom);}
    // Sets ITS geometry For Reconstruction! be very careful using this fun.
    void   SetITSgeomRec(AliITSgeom *geom) {if(fDetTypeRec!=0)
        fDetTypeRec->SetITSgeom(geom);}
    // return pointer to the array of modules
    TObjArray    *GetModules() const {return fDetTypeSim->GetModules();}
    // return pointer to a particular module
    AliITSmodule *GetModule(Int_t index){return fDetTypeSim->GetModule(index);}
#endif

    //================ Necessary general Classes =======================
    virtual void Init();
    virtual AliLoader* MakeLoader(const char* topfoldername);
    virtual void SetDefaults();
    virtual void SetDefaultSimulation();
    virtual void SetDefaultClusterFinders();
    virtual void SetDefaultClusterFindersV2();
    virtual void MakeBranch(Option_t *opt=" ");
    virtual void SetTreeAddress();
#ifndef NEWVERSION
    // For a given branch from the treeH sets the TClonesArray address.
    virtual void SetHitsAddressBranch(TBranch *b) {b->SetAddress(&fHits);}
    // Return pointer to DetType #id
    AliITSDetType *DetType(Int_t id){
        return ((AliITSDetType*) fDetTypes->At(id));};
    //Int_t           NDetTypes() {return fNDetTypes;}
#endif
    //---------- Configuration Methods (per detector type) -------------
    // Determines which ITS subdetectors will be processed. Effects
    // digitization, and Reconstruction only.
    void SetDetectors(Option_t *opt="All"){fOpt = opt;}
    // Returns the list of ITS subdetectors that will be processed.
    Option_t* GetDetectors(){return fOpt;}
#ifndef NEWVERSION
    // Set response 
    virtual void SetResponseModel(Int_t id, AliITSresponse *response){
        ((AliITSDetType*) fDetTypes->At(id))->ResponseModel(response);};
    // Set segmentation 
    virtual void SetSegmentationModel(Int_t id, AliITSsegmentation *seg){
        ((AliITSDetType*) fDetTypes->At(id))->SegmentationModel(seg);};
    // Set simulation - temporary 
    virtual void SetSimulationModel(Int_t id, AliITSsimulation *sim){
        ((AliITSDetType*) fDetTypes->At(id))->SimulationModel(sim);};
    // Set simulation - temporary 
    virtual AliITSsimulation* GetSimulationModel(Int_t id){
	return ((AliITSDetType*)(fDetTypes->At(id)))->GetSimulationModel();}
    // Set reconstruction 
    virtual void SetReconstructionModel(Int_t id, AliITSClusterFinder *rec){
        ((AliITSDetType*) fDetTypes->At(id))->ReconstructionModel(rec);};
    // Set class names for digit and rec point 
    virtual void SetClasses(Int_t id, const char *digit, const char *cluster){
        ((AliITSDetType*) fDetTypes->At(id))->ClassNames(digit,cluster);};
#endif
#ifdef NEWVERSION
    // Set response
    virtual void SetResponseModel(Int_t module, AliITSresponse *response){
        fDetTypeSim->SetResponseModel(module,response);};
    // Set segmentation for Simulation
    virtual void SetSegmentationModelSim(Int_t id, AliITSsegmentation *seg){
        fDetTypeSim->SetSegmentationModel(id,seg);};
    // Set segmentation for Reconstruction
    virtual void SetSegmentationModelRec(Int_t id, AliITSsegmentation *seg){
        fDetTypeRec->SetSegmentationModel(id,seg);};
    // Set segmentation 
    virtual void SetSegmentationModel(Int_t id, AliITSsegmentation *seg){
        SetSegmentationModelSim(id,seg);SetSegmentationModelRec(id,seg);};
    // Set simulation 
    virtual void SetSimulationModel(Int_t id, AliITSsimulation *sim){
        fDetTypesSim->SetSimulationModel(sim);};
    // Set simulation 
    virtual AliITSsimulation* GetSimulationModel(Int_t id){
	return fDetTypesSim->GetSimulationModel(id);}
    // Set Calibration
    virtual void SetCalibrationModel(Int_t module, AliITSCalibration *cal){
        fDetTypeRec->SetCalibrationModel(module,cal);};
    // Set reconstruction 
    virtual void SetReconstructionModel(Int_t id, AliITSClusterFinder *rec){
        fDetTypesRec->SetReconstructionModel(id,rec);};
    // Set Class name for Hits
    virtual void SetHitClassName(){
        fDetTypeSim->SetHitClassName(this->GetName());}
    // Set Class name for SDigits
    virtual void SetSDigitClassName(const char *sdigit){
        fDetTypeSim->SetSDigitClassName(sdigit);}
    // Set Class name for Digits for simulation
    virtual void SetDigitClassNameSim(const char *digit){
        fDetTypeSim->SetDigitClassName(digit);}///////// Array of names
    // Set Class name for Digits for Reconstruction
    virtual void SetDigitClassNameRec(const char *digit){
        fDetTypeRec->SetDigitClassName(digit);}///////// Array of names
    virtual void SetClusterClassName(const char *digit){
        fDetTypeRec->SetClusterClassName(digit);}///////// Array of names
    virtual void SetRecPointClassName(const char *digit){
        fDetTypeRec->SetRecPointClassName(digit);}///////// Array of names
    // Set class names for digit and rec point 
    virtual void SetClasses(Int_t id, const char *digit, const char *cluster){
        SetDigitClassNameSim(digit);SetDigitClassNameRec(digit);
        SetRecPointClassName(cluster);};
#endif

    //=================== Hits =========================================
    virtual void StepManager() {} // See Step Manager for specific geometry.
    //------------ sort hits by module for Digitisation ----------------
    virtual void FillModules(Int_t evnt,Int_t bgrev,Int_t nmodules,
			     Option_t *opt, const char *filename);
#ifndef NEWVERSION
    virtual void AddHit(Int_t track, Int_t *vol, Float_t *hits);
    virtual void InitModules(Int_t size,Int_t &nmodules);  
    virtual void FillModules(TTree *treeH, Int_t mask = 0);
    virtual void ClearModules(){if(fITSmodules) fITSmodules->Delete();};
#endif
#ifdef NEWVERSION
    virtual void AddHit(Int_t track, Int_t *vol, Float_t *hits){
        if(fDetTypeSim!=0)fDetTypeSim->AddHit(track,vol,hits);};
    virtual void AddHit(AliITShit &hit){
        if(fDetTypeSim!=0)fDetTypeSim->AddHit(hit);};
    virtual void InitModules(Int_t size,Int_t &nmodules){if(fDetTypeSim!=0)
        fDetTypeSim->InitModules(size,nmodules);};
    virtual void FillModules(TTree *treeH, Int_t mask = 0){if(fDetTypeSim!=0)
        fDetTypeSim->FillModules(treeH,mask);};
    virtual void ClearModules(){if(fDetTypeSim!=0) 
        fDetTypeSim->ClearModules();};
#endif

    //===================== Digitisation ===============================
#ifndef NEWVERSION
    void MakeBranchS(const char *file);
    void SetTreeAddressS(TTree *treeS);
    TClonesArray * GetSDigits() { return fSDigits; }
    void MakeBranchInTreeD(TTree *treeD,const char *file=0);
    void MakeBranchD(const char *file){
	MakeBranchInTreeD(GetLoader()->TreeD(),file);}
    void SetTreeAddressD(TTree *treeD);
#endif
#ifndef NEWVERSION
    void Hits2SDigits(); // Turn hits into SDigits
    void Hits2PreDigits(){ // Turn hits into SDigits
        HitsToPreDigits(fLoader->GetRunLoader()->GetEventNumber(),
                        0,-1," ",fOpt," ");};
    AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
    void SDigits2Digits(){SDigitsToDigits("All");} // Turn SDigits to Digits
    void SDigitsToDigits(Option_t *opt="All"); // Turn SDigits to Digits
    void Hits2Digits(); // Turn hits straight into Digits.
    //------------------ Internal functions ----------------------------
    // Standard Hits To SDigits function
    void HitsToSDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                       Option_t *add, Option_t *det, const char *filename)
        {HitsToPreDigits(evNumber,bgrev,size,add,det,filename);};
    // Standard Hits To SDigits function
    void HitsToPreDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                 Option_t *add, Option_t *det, const char *filename);
    // Standard Hits To Digits function
    void HitsToDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                 Option_t *add, Option_t *det, const char *filename);
#endif
#ifdef NEWVERSION
    // Turn hits into SDigits
    void Hits2SDigits(){if(fDetTypeSim)fDetTypeSim->Hits2SDigits();};
    AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
    // Turn SDigits to Digits
    void SDigits2Digits(){if(fDetTypeSim)fDetTypeSim->SDigits2Digits();} 
    // Turn hits straight into Digits.
    void Hits2Digits(){if(fDetTypeSim)fDetTypeSim->Hits2Digits();};
#endif
    // Resets the Summable digits.
#ifndef NEWVERSION
    void ResetSDigits(){if(fSDigits) fSDigits->Clear();fNSDigits = 0;};
    void ResetDigits();                   // depending on how the
    void ResetDigits(Int_t branch);       // tree will be filled only
    void AddSumDigit(AliITSpListItem &sdig);
    void AddRealDigit(Int_t branch, Int_t *digits);
    void AddSimDigit(Int_t branch, AliITSdigit *d);
    void AddSimDigit(Int_t branch,Float_t phys,Int_t* digits,
		     Int_t* tracks,Int_t *hits,Float_t* trkcharges);
#endif
    // Return pointers to digits 
#ifndef NEWVERSION
    TObjArray    *Dtype() {return fDtype;}
    Int_t        *Ndtype() {return fNdtype;}
    TClonesArray *DigitsAddress(Int_t id)
	{return ((TClonesArray *) (*fDtype)[id]);}
#endif
#ifdef NEWVERSION
    TObjArray* GetDigitsSim(){if(fDetTypeSim!=0) 
        return fDetTypeSim->GetDigits();else return 0;}
    TObjArray* GetDigitsRec(){if(fDetTypeRec!=0) 
        return fDetTypeRec->GetDigits();else return 0;}
    TObjArray    *Dtype() {if(fDetTypeSim!=0) return GetDigitsSim();
    else if(fDetTypeRec!=0) return GetDigitsRec(); else return 0;}
    Int_t* GetNDigitArraySim(){if(fDetTypeSim!=0) 
        fDetTypeSim->GetNDigitArray();}
    Int_t* GetNDigitArrayRec(){if(fDetTypeRec!=0) 
        fDetTypeRec->GetNDigitArray();}
    Int_t        *Ndtype() {if(fDetTypeSim!=0) return GetNDigitArraySim();
        else if(fDetTypeRec!=0) return GetNDigitArrayRec(); else return 0;}
    TClonesArray *DigitsAddressSim(Int_t id){if(fDetTypeSim!=0)
	return fDetTypeSim->GetDigitsAddress(id);else return 0;}
    TClonesArray *DigitsAddressRec(Int_t id){if(fDetTypeRec!=0)
	return fDetTypeRec->GetDigitsAddress(id);else return 0;}
    TClonesArray *DigitsAddress(Int_t id){if(fDetTypeSim!=0)
	return DigitsAddressSim(id);else if(fDetTypeRec!=0) DigitsAddressRec(id);
     else return 0;}
#endif
    void SelectVertexer(TString sel=" "){fSelectedVertexer = sel;}
#ifndef NEWVERSION

    //===================== Raw Data IO ================================
    // Write digits into raw data format
    virtual void Digits2Raw();

    //==================== Clusterization ==============================
    // create separate tree for clusters - declustering refining
    void MakeTreeC(Option_t *option="C");
    void GetTreeC(Int_t event);
    void AddCluster(Int_t branch, AliITSRawCluster *c);
    // one of the methods in
    void ResetClusters(){for(Int_t i=0;i<kNTYPES;i++ ) ResetClusters(i);}; 
    void ResetClusters(Int_t branch);     // the pair will be kept
    void MakeBranchC();
    // Return pointers to clusters 
    TObjArray    *Ctype() {return fCtype;}
    Int_t        *Nctype() {return fNctype;}
    TClonesArray *ClustersAddress(Int_t id) 
                   {return ((TClonesArray *) (*fCtype)[id]);}
#endif
#ifdef NEWVERSION

    //===================== Raw Data IO ================================
    // Write digits into raw data format
    virtual void Digits2Raw(){if(fDetTypeSim)fDetTypeSim->Digits2Raw();};

    //==================== Clusterization ==============================
    // create separate tree for clusters - declustering refining
    void MakeTreeC(Option_t *option="C"){fDetTypeRec->MakeTreeC();};
    void GetTreeC(Int_t event){fDetTypeRec->GetTreeC(event);};
    void AddCluster(Int_t branch, AliITSRawCluster *c){
        fDetTypeRec->AddCluster(branch,c);};
    // one of the methods in
    void ResetClusters(){for(Int_t i=0;i<kNTYPES;i++ ) ResetClusters(i);}; 
    void ResetClusters(Int_t branch){fDetTypeRec->ResetCluster(i);};
    void MakeBranchC(){fDetTypeRec->MakeBranchC();};
    // Return pointers to clusters 
    TObjArray    *Ctype() {if(fDetTypeRec!=0)
        return fDetTypeRec->GetClusterArray(); else return 0;}
    Int_t        *Nctype() {if(fDetTypeRec!=0) 
        return fDetTypeRec->GetNClusters(); else return 0;;}
    TClonesArray *ClustersAddress(Int_t id){if(fDetTypeRec!=0) 
        return fDetTypeRec->GetClusterAddress(id]); else return 0;}
#endif

    //=================== Reconstruction ===============================
#ifndef NEWVERSION
    void MakeBranchR(const char *file, Option_t *opt=" ");
    void SetTreeAddressR(TTree *treeR);
    void AddRecPoint(const AliITSRecPoint &p);
    void ResetRecPoints(){if(fRecPoints) fRecPoints->Clear();fNRecPoints = 0;};
    // Return pointer to rec points 
    TClonesArray  *RecPoints()   {return fRecPoints;}

    void AddClusterV2(const AliITSclusterV2 &cl);
    void ResetClustersV2(){if(fClustersV2) fClustersV2->Clear();fNClustersV2=0;} 
    Int_t GetNClustersV2()const {return fNClustersV2;}
// Return pointer to clustersV2
TClonesArray *ClustersV2() {return fClustersV2;}
#endif
#ifdef NEWVERSION
    void MakeBranchR(const char *file, Option_t *opt=" ");
    void SetTreeAddressR(TTree *treeR){fDetTypeRec->SetTreeAddressR(treeR);};
    void AddRecPoint(const AliITSRecPoint &p){fDetTypeRec->AddRecPoint(p);};
    void ResetRecPoints(){if(fDetTypeRec) fDetTypeRec->ResetRecPoints();};
    // Return pointer to rec points 
    TClonesArray* RecPoints() {if(fDetTypeRec!=0)
        return fDetTypeRec->GetRecPoints();else return 0;}
#endif
    void MakeBranchRF(const char *file){MakeBranchR(file,"Fast");}
    void HitsToFastRecPoints(Int_t evNumber,Int_t bgrev,Int_t size,
                 Option_t *add, Option_t *det, const char *filename);
    void Digits2Reco(){
        DigitsToRecPoints(fLoader->GetRunLoader()->GetEventNumber(),0,fOpt);};
    void DigitsToRecPoints(Int_t evNumber,Int_t lastEntry,Option_t *det);
void DigitsToRecPoints(AliRawReader* rawReader);
 protected:
    //================== Data Members ==================================
#ifdef NEWVERSION
    AliITSDetTypeSim *fDetTypeSim; //
    AliITSDetTypeRec *fDetTypeRec; //
#endif
#ifndef NEWVERSION
    AliITSgeom   *fITSgeom;    // Pointer to ITS geometry
#endif
    Bool_t        fEuclidOut;  // Flag to write geometry in euclid format
#ifndef NEWVERSION
    TObjArray    *fITSmodules; //! Pointer to ITS modules
#endif
    Option_t     *fOpt;        //! Detector option ="All" unless changed.

    Int_t         fIdN;        // the number of layers
    Int_t        *fIdSens;     //[fIdN] layer identifier
    TString      *fIdName;     //[fIdN] layer identifier
#ifndef NEWVERSION
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

    TClonesArray *fClustersV2; //!List of reconstructed clusters v2
    Int_t         fNClustersV2;    //Number of clusters v2

#endif
    TString fSelectedVertexer; // Vertexer selected in CreateVertexer
#ifndef NEWVERSION
    ClassDef(AliITS,5) // Base class for ITS
#endif
#ifdef NEWVERSION
    ClassDef(AliITS,5) // Base class for ITS
#endif
};

#endif
