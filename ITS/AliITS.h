#ifndef ALIITS_H
#define ALIITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//           Manager class for set: ITS                               //
////////////////////////////////////////////////////////////////////////


#include "AliDetector.h"


class AliITSDetType;
class AliITSsimulation;
class AliITSClusterFinder;
class AliITSsegmentation;
class AliITSresponse;
class AliITShit;
class AliITSgeom;
class AliITSdigit;
class AliITSRecPoint;
class AliITSRawCluster;
class AliITSmodule;

class TObjArray;
class TString;
class TTree;



class AliITS : public AliDetector {

  public:
                    AliITS();
                    AliITS(const char *name, const char *title);
    virtual        ~AliITS();
    AliITS(AliITS &source);
    AliITS & operator=(AliITS &source);

    virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits);
    virtual void   AddRealDigit(Int_t branch, Int_t *digits);
    virtual void   AddSimDigit(Int_t branch, AliITSdigit *d);
    virtual void   AddSimDigit(Int_t branch,Float_t phys,Int_t* digits,Int_t* tracks,Int_t *hits,Float_t* trkcharges); 
    //
    virtual void   AddCluster(Int_t branch, AliITSRawCluster *c);
    virtual void   AddRecPoint(const AliITSRecPoint &p);

    virtual void   ResetDigits();                   // depending on how the
    virtual void   ResetDigits(Int_t branch);       // tree will be filled only
    virtual void   ResetClusters();                 // one of the methods in 
    virtual void   ResetClusters(Int_t branch);     // the pair will be kept
    virtual void   ResetRecPoints();

    // get geometry version - detailed (major) or coarse (minor)
    void GetGeometryVersion(Int_t &a,Int_t &b) const 
	           {a = fMajorVersion;b=fMinorVersion;return;}
    virtual Int_t  IsVersion() const {return 1;}
            Int_t  DistancetoPrimitive(Int_t px, Int_t py);
    virtual void   Init();
    virtual void   SetDefaults();
    virtual void   SetDefaultSimulation();
    virtual void   SetDefaultClusterFinders();
    // create separate tree for clusters - declustering refining
    virtual  void  MakeTreeC(Option_t *option="C");
    void           GetTreeC(Int_t event);
    virtual void   MakeBranch(Option_t *opt=" ");
    void           SetTreeAddress();
    virtual void   SetEUCLID(Bool_t euclid=1) {fEuclidOut = euclid;}
    virtual void   StepManager() {}
    // sort hits by module
    virtual void   InitModules(Int_t size,Int_t &nmodules);  
    virtual void   FillModules(Int_t evnt,Int_t bgrev,
                       Int_t nmodules,Option_t *opt,Text_t *filename);
    virtual void   ClearModules();
    // Digitisation
    void HitsToDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                 Option_t *add, Option_t *det, Text_t *filename);
    // Reconstruct hits
    void DigitsToRecPoints(Int_t evNumber,Int_t lastEntry,Option_t *det);
    // Fast simulation of space points from hits
    void HitsToFastRecPoints(Int_t evNumber,Int_t bgrev,Int_t size,
                 Option_t *add, Option_t *det, Text_t *filename);

    // Write digits into raw data format
    virtual void Digits2RawData() {}
    // Decode raw data and store digits
    virtual void RawData2Digits() {}

    // Configuration Methods (per detector type )
    // Set response 
    virtual void   SetResponseModel(Int_t id, AliITSresponse *response);
    // Set segmentation 
    virtual void   SetSegmentationModel(Int_t id, AliITSsegmentation *seg);
    // Set simulation - temporary 
    virtual void   SetSimulationModel(Int_t id, AliITSsimulation *sim);
    // Set reconstruction 
    virtual void   SetReconstructionModel(Int_t id, AliITSClusterFinder *rec);
    // Set class names for digit and rec point 
    virtual void   SetClasses(Int_t id, const char *digit, const char *cluster);


    // Getters
    // ITS geometry functions
    virtual AliITSgeom  *GetITSgeom() const {return fITSgeom;}
    // return pointer to the array of modules
    virtual TObjArray   *GetModules() const {return fITSmodules;}
    // return pointer to a particular module
    AliITSmodule *GetModule(Int_t index) {return (AliITSmodule *)
					      (fITSmodules->At(index));}

    // Return pointers to digits 
    TObjArray            *Dtype() {return fDtype;}
    Int_t                *Ndtype() {return fNdtype;}
    virtual TClonesArray *DigitsAddress(Int_t id) 
                   {return ((TClonesArray *) (*fDtype)[id]);}
    // Return pointers to clusters 
    TObjArray            *Ctype() {return fCtype;}
    Int_t                *Nctype() {return fNctype;}
    virtual TClonesArray *ClustersAddress(Int_t id) 
                   {return ((TClonesArray *) (*fCtype)[id]);}

    // Return pointer to rec points 
    TClonesArray  *RecPoints()   {return fRecPoints;}

    // Return pointer to DetType #id
    AliITSDetType  *DetType(Int_t id);
    //Int_t           NDetTypes() {return fNDetTypes;}

    // Return pointer to the tree of clusters
    TTree          *TreeC() {return fTreeC;}


 protected:

    static const Int_t fgkNTYPES;          // Number of detector types

    AliITSgeom  *fITSgeom;                   // Pointer to ITS geometry
    TObjArray   *fITSmodules;                // Pointer to ITS modules

    Bool_t fEuclidOut;                       // Flag to write out geometry 
                                             // in euclid format
    Int_t  fIdN;                             // the number of layers
    Int_t  *fIdSens;                         //[fIdN] layer identifier
    TString *fIdName;                        //[fIdN] layer identifier
    // Geometry and Stepmanager version numbers used.
    Int_t fMajorVersion,fMinorVersion;       //detailed and coarse(minor) versions
    //
    //

    // Int_t             fNDetTypes;           //Number of detector types
    TObjArray            *fDetTypes;           // List of detector types

    TObjArray            *fDtype;              //[fgkNTYPES]  List of digits
    Int_t                *fNdtype;             //[fgkNTYPES] Number of digits per type of
                                               // detector 
    TObjArray            *fCtype;              //[fgkNTYPES] List of clusters
    Int_t                *fNctype;             //[fgkNTYPES] Number of clusters per type
                                               // of detector

    TClonesArray         *fRecPoints;          // List of reconstructed points
    Int_t                 fNRecPoints;         // Number of rec points
    TTree                *fTreeC;              // Tree for raw clusters


    ClassDef(AliITS,1) // Base class for ITS
};

#endif
