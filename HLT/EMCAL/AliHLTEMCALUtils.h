#ifndef ALIHLTEMCALUTILS_H
#define ALIHLTEMCALUTILS_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTEMCALUtils.h
    @author m. ploskon
    @date   
    @brief  Utils for EMCAL-HLT
*/

#include <TObject.h>

/**
 * @class AliHLTEMCALUtils
 */

class AliEMCALGeometry; 
class AliEMCALClusterizerv1;
class AliEMCALRecParam;    
class AliEMCALRawUtils;
class AliRawReader;
class AliESDEvent;

class TString;
class TFile;
class TGeoManager;
class TTree;
class TClonesArray;
class TBranch;

class AliHLTEMCALUtils : public TObject
{
 public:
  AliHLTEMCALUtils();
  AliHLTEMCALUtils(const AliHLTEMCALUtils & /*t*/);
  AliHLTEMCALUtils & operator = (const AliHLTEMCALUtils & /*t*/);
  virtual ~AliHLTEMCALUtils();

  static void                   SetDebug(Int_t idbg) {fgDebug = idbg;}

  static void                   DeleteStaticMembers();
  static void                   Cleanup();
  
  static void                   InitRecParam();
  static AliEMCALRecParam*      GetRecParam();

  static AliEMCALRawUtils*      GetRawUtils(AliEMCALGeometry *pGeometry = 0);
  static AliEMCALClusterizerv1* GetClusterizer();
  static AliEMCALGeometry*      GetGeometry();

  static void                   DeleteGeoManager();
  static Bool_t                 LoadGeoManagerFromFile(const char *fname);
  static TGeoManager*           GetGeoManager() {return fgGeoManager;}

  static Bool_t                 InitFakeCDB(const char *cdbpath = "", Int_t runid = 0);
  static Bool_t                 ConvertDigits(AliRawReader* rawReader, TTree* digitsTree, Option_t* sOption);
  static Bool_t                 Raw2Clusters(AliRawReader* rawReader, TTree* clustersTree, Option_t* sDigitsOption);
  static Bool_t                 RawBuffer2Clusters(UChar_t *buffer, ULong_t buffSize, 
						   Int_t eqID, 
						   TTree* clustersTree, Option_t* sDigitsOption);
  static TTree*                 RawBuffer2Clusters(UChar_t *buffer, ULong_t buffSize, 
						   Int_t eqID, 
						   Option_t* sDigitsOption);

  static TTree*                 GetDigitsTree();
  static TTree*                 GetClustersTree();

  static UChar_t *              FileToMemory(const char *fname, UChar_t *outpbuffer, ULong_t &buffSize);

  static void                   ResetReconstructionTrees();
  static void                   DeleteReconstructionTrees();

  static Bool_t                 FillESD(TTree* digitsTree, TTree* clustersTree, AliESDEvent* esd);
    
 protected:

 private:
  
  static Int_t                  fgDebug; // debug flag

  static AliEMCALGeometry*      fgGeom;        // EMCAL geometry
  static AliEMCALClusterizerv1* fgClusterizer; // ECMAL clusterizer
  static AliEMCALRecParam*      fgRecParam;    // EMCAL reconstruction parameters
  static AliEMCALRawUtils*      fgRawUtils;    // EMCAL raw utilities 
  
  static TFile*                 fgGeometryFile; //! // Pointer to the geom root file
  static TGeoManager*           fgGeoManager; //! Pointer to geometry manager 

  static TTree*                 fgClustersTreeTmp; //! Clusters tree - depending on the reco scheme 
                                                   //  use ResetReconstructionTrees
  static TTree*                 fgDigitsTreeTmp;   //! Digits tree  - reset on each Raw2Clusters call
  static TClonesArray*          fgDigitsArrTmp;    //! Digits array - reset on each Raw2Clusters call
  static TBranch*               fgDigitsBranchTmp; //! Digits branch - reset on each Raw2Clusters call

  ClassDef(AliHLTEMCALUtils, 0)
};

#endif
