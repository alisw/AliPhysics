#ifndef ALIITSDETTYPEREC_H
#define ALIITSDETTYPEREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
$Id$ 
*/

////////////////////////////////////////////////////////////////////////
// This class contains all of the "external" information needed to do //
// detector specific reconstruction for the ITS.                      //
////////////////////////////////////////////////////////////////////////
#include <TObject.h>
#include <TClonesArray.h>
#include <TBits.h>

class TObjArray;
class TTree;
class TBranch;
class AliITSgeom;
class AliITSsegmentation;
class AliITSCalibration;
class AliITSCalibrationSSD;
class AliITSresponseSDD;
class AliITSClusterFinder;
class AliITSRecPoint;
class AliRawReader;
class AliITSGainSSDv2;
class AliITSBadChannelsSSDv2;
class AliITSDDLModuleMapSDD;
class AliITSNoiseSSDv2;
class AliITSTriggerConditions;
class AliITSFOSignalsSPD;
class AliITSRecPointContainer;

class AliITSDetTypeRec : public TObject {
  public:
    AliITSDetTypeRec(); // Default constructor
 
    virtual ~AliITSDetTypeRec(); // Proper Destructor

    virtual AliITSgeom* GetITSgeom() const { return fITSgeom; }
    virtual void SetITSgeom(AliITSgeom* const geom) { fITSgeom = geom; }
    virtual void SetDefaults();
    virtual void SetDefaultClusterFindersV2(Bool_t rawdata=kFALSE);
    virtual void MakeBranch(TTree *tree,Option_t *opt);
    virtual void SetTreeAddressD(TTree* const treeD);

    virtual void SetSegmentationModel(Int_t dettype, AliITSsegmentation *seg);
    virtual void SetCalibrationModel(Int_t iMod, AliITSCalibration *cal);
    virtual void SetSPDDeadModel(Int_t iMod, AliITSCalibration *cal);
    virtual void SetReconstructionModel(Int_t dettype, AliITSClusterFinder *rec);
    virtual Bool_t GetCalibration();
    virtual AliITSsegmentation* GetSegmentationModel(Int_t dettype) const;
    virtual AliITSCalibration* GetCalibrationModel(Int_t iMod) const;
    virtual AliITSCalibration* GetSPDDeadModel(Int_t iMod) const;
    virtual AliITSTriggerConditions* GetTriggerConditions() const;
    virtual AliITSClusterFinder* GetReconstructionModel(Int_t dettype) const;
    virtual AliITSDDLModuleMapSDD* GetDDLModuleMapSDD() const { return fDDLMapSDD;}
    virtual AliITSresponseSDD* GetResponseSDD() const { return fRespSDD;}
    virtual Float_t GetAverageGainSDD() const {
      if(fAveGainSDD>0.) return fAveGainSDD;
      else return 1.;
    }

    virtual void SetDigitClassName(Int_t i,const Char_t *digit) 
      {fkDigClassName[i]=digit;}
    
    virtual void SetLoadOnlySPDCalib(Bool_t opt=kFALSE)
      {fLoadOnlySPDCalib=opt;}

    const Char_t* GetDigitClassName(Int_t i) const {return fkDigClassName[i];}
    
    TObjArray* GetDigits() const {return fDigits;} 
    TClonesArray *DigitsAddress(Int_t id) const {return ((TClonesArray*)(*fDigits)[id]);}

        AliITSFOSignalsSPD* GetFOSignals() const {return fFOSignals;}
    
    TBranch* MakeBranchInTree(TTree* const tree, const char* name, const char *classname, void* address,Int_t size, Int_t splitlevel);

    virtual void ResetDigits();
    virtual void ResetDigits(Int_t branch);

    void MakeBranchR(TTree *treeR,Option_t *opt=" ");
    void SetTreeAddressR(TTree* const treeR);
    void AddRecPoint(const AliITSRecPoint &p);
    void ResetRecPoints(){if(fRecPoints) fRecPoints->Clear();fNRecPoints = 0;};
    // Return pointer to rec points 
    TClonesArray  *RecPoints() const  {return fRecPoints;}
    void MakeBranchRF(TTree *treeR){MakeBranchR(treeR,"Fast");}
    void DigitsToRecPoints(TTree *treeD,TTree *treeR,Int_t lastEntry,Option_t *det, Int_t optCluFind=0);
    void DigitsToRecPoints(AliRawReader* rawReader,TTree *treeR,Option_t *det="All");
    void DigitsToRecPoints(AliRawReader* rawReader,TClonesArray** clusters,Option_t *opt);

    void   SetFastOrFiredMapOnline(UInt_t eq, UInt_t hs, UInt_t chip);
    void   SetFastOrFiredMap(UInt_t chipKey){fFastOrFiredMap.SetBitNumber(chipKey);} 
    TBits  GetFastOrFiredMap() const {return fFastOrFiredMap;}
    TBits  GetFiredChipMap(TClonesArray **clusters) const; // (using SPD RecPoints)
    TBits  GetFiredChipMap(TTree *treeR) const; // (using SPD RecPoints)
    void   ResetFastOrFiredMap(){fFastOrFiredMap.ResetAllBits();}
    void   RemoveFastOrFiredInActive(); // (using Trigger Conditions)
    void   RemoveFastOrFiredFromDead(TBits firedchipmap); // (using SPD RecPoints)
   
    
  private:
    // private methods
    AliITSDetTypeRec(const AliITSDetTypeRec& rec);
    AliITSDetTypeRec& operator=(const AliITSDetTypeRec &source);
 
    //conversion of the old SSD calibration objects tothe new ones
    void ReadOldSSDNoise(const TObjArray *array, 
			 AliITSNoiseSSDv2 *noiseSSD);
    void ReadOldSSDBadChannels(const TObjArray *array, 
			       AliITSBadChannelsSSDv2 *badChannelsSSD);
    void ReadOldSSDGain(const TObjArray *array, 
			AliITSGainSSDv2 *gainSSD);
    virtual Bool_t GetCalibrationSPD(Bool_t cacheStatus);
    virtual Bool_t GetCalibrationSDD(Bool_t cacheStatus);
    virtual Bool_t GetCalibrationSSD(Bool_t cacheStatus);

    //    virtual void SetLoader(AliITSLoader* loader) {fLoader=loader;}
    static const Int_t fgkNdettypes;          // number of det. types
    static const Int_t fgkDefaultNModulesSPD; // Total numbers of SPD modules by default
    static const Int_t fgkDefaultNModulesSDD; // Total numbers of SDD modules by default
    static const Int_t fgkDefaultNModulesSSD; // Total numbers of SSD modules by default
    Int_t *fNMod;     // numbers of modules from different types

    AliITSgeom   *fITSgeom;       //! ITS geometry

    TObjArray    *fReconstruction;//! [NDet]
    TObjArray    *fSegmentation;  //! [NDet]
    TObjArray    *fCalibration;   //! [NMod]
    AliITSCalibrationSSD* fSSDCalibration;  //! SSD calibration object
    TObjArray    *fSPDDead;       //! [fgkDefaultNModulesSPD]
    AliITSTriggerConditions *fTriggerConditions; //! PIT conditions object
    TObjArray    *fDigits;        //! [NMod][NDigits]
    AliITSFOSignalsSPD *fFOSignals; //! Fast-Or signals (used when reconstructing from digits)
    AliITSDDLModuleMapSDD *fDDLMapSDD; //! mapping DDL/module -> SDD module number
    AliITSresponseSDD *fRespSDD;  //! SDD response parameters 
    Float_t       fAveGainSDD;    //! Average gain of SDD good anodes
    const Char_t*       fkDigClassName[3];     //! String with digit class name.


    TClonesArray *fRecPoints;  //! List of reconstructed points
    Int_t         fNRecPoints; // Number of rec points

    Bool_t fFirstcall;         //! flag
    Bool_t fLoadOnlySPDCalib;  //! flag for loading calibrations only fr SPD

    TBits fFastOrFiredMap;     //! Map of FastOr fired chips (after processing of raw signals)

    ClassDef(AliITSDetTypeRec,19) // ITS Reconstruction structure
};

#endif

