#ifndef AliRICH_h
#define AliRICH_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include <TObjArray.h>
#include <AliDetector.h>
#include "AliRICHConst.h"
#include "AliRICHChamber.h"
static const int kNCH=7;

class TFile;

class AliRICHHit;
class AliRICHSDigit;
class AliRICHRawCluster;
class AliRICHRecHit1D;
class AliRICHRecHit3D;
class AliRICHClusterFinder;
class AliRICHDetect;
class AliRICHChamber;
class AliRICHCerenkov;
class AliSegmentation;
class AliRICHResponse;
class AliRICHGeometry;
class AliRICHMerger;

class AliRICH : public AliDetector 
{
public:
                    AliRICH();                                            //default ctor
                    AliRICH(const char *name, const char *title);         //named ctor
  inline            AliRICH(const AliRICH& RICH)                    {;}   //copy ctor  
          virtual  ~AliRICH();                                            //dtor
          
  inline  AliRICH& operator=(const AliRICH& rhs) { return *this;}
          virtual Int_t  IsVersion() const =0;
          
  virtual void  AddHit(Int_t track, Int_t *vol, Float_t *hits);
          void  AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs);
          void  AddSDigit(Int_t *clhits);
          void  AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits);
          void  AddRawCluster(Int_t id, const AliRICHRawCluster& cluster);
          void  AddRecHit1D(Int_t id, Float_t* rechit, Float_t* photons, Int_t* padsx, Int_t* padsy);
          void  AddRecHit3D(Int_t id, Float_t* rechit, Float_t omega, Float_t theta, Float_t phi);
  virtual void  ResetHits();
  virtual void  ResetDigits();
          void  ResetRawClusters();
          void  ResetRecHits1D();
          void  ResetRecHits3D();
  virtual void  FindClusters(Int_t nev,Int_t lastEntry);
  virtual void  Hits2SDigits();
          Int_t Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss,Int_t id, ResponseType res);
  virtual void  SDigits2Digits();
  virtual void  Digits2Reco();

  virtual void    CreateMaterials(); //GEANT materials definition
          Float_t AbsoCH4(Float_t x);
          Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola);
  virtual void    BuildGeometry();   //TNode ROOT variant for event display
  virtual void    CreateGeometry();  //GEANT volumes tree for simulation  
  virtual void    StepManager()=0;
   
  inline Int_t    DistancetoPrimitive(Int_t px, Int_t py)      {return 9999;}
   
  virtual void   MakeBranch(Option_t *opt=" ");
  virtual void   MakeBranchInTreeD(TTree *treeD, const char *file=0);
  virtual void   SetTreeAddress();
   
  
   
   
  AliRICHSDigit* FirstPad(AliRICHHit *hit, TClonesArray *clusters);
  AliRICHSDigit* NextPad(TClonesArray *clusters);
   

  void     SetGeometryModel(Int_t iChamberN, AliRICHGeometry *pRICHGeo)    {       GetChamber(iChamberN)->SetGeometryModel(pRICHGeo);}
  AliRICHGeometry* GetGeometryModel(Int_t iChamberN=0)                        const{return GetChamber(iChamberN)->GetGeometryModel();}    
  void     SetSegmentationModel(Int_t iChamberN, AliSegmentation *pAliSeg) {       GetChamber(iChamberN)->SetSegmentationModel(pAliSeg);}
  AliSegmentation* GetSegmentationModel(Int_t iChamberN=0)                    const{return GetChamber(iChamberN)->GetSegmentationModel();}
  void     SetResponseModel(Int_t iChamberN, AliRICHResponse *pRICHRes)    {       GetChamber(iChamberN)->SetResponseModel(pRICHRes);}
  AliRICHResponse* GetResponseModel(Int_t iChamberN)                          const{return GetChamber(iChamberN)->GetResponseModel();}
  void     SetReconstructionModel(Int_t iChamberN, AliRICHClusterFinder *pRICHReco){GetChamber(iChamberN)->SetReconstructionModel(pRICHReco);}

  virtual void   SetMerger(AliRICHMerger* thisMerger) {fMerger=thisMerger;}  
  AliRICHChamber& Chamber(Int_t id) {return *((AliRICHChamber *) (*fChambers)[id]);}
  AliRICHChamber* GetChamber(Int_t iChamberN)     const{return (AliRICHChamber*) (*fChambers)[iChamberN];}
  
  inline TObjArray     *Dchambers()                     {return fDchambers;}
  inline TObjArray     *RecHits3D()                const{return fRecHits3D;}
  inline TObjArray     *RecHits1D()                const{return fRecHits1D;}
  inline Int_t         *Ndch()                          {return fNdch;}
  inline Int_t         *Nrechits1D()                    {return fNrechits1D;} 
  inline Int_t         *Nrechits3D()                    {return fNrechits3D;} 
  inline TClonesArray  *SDigits()                  const{return fSDigits;}
  inline TClonesArray  *Cerenkovs()                const{return fCerenkovs;}
  inline TClonesArray  *DigitsAddress(Int_t id)         {return ((TClonesArray *) (*fDchambers)[id]);}
  inline TClonesArray  *RecHitsAddress1D(Int_t id) const{return ((TClonesArray *) (*fRecHits1D)[id]);}
  inline TClonesArray  *RecHitsAddress3D(Int_t id) const{return ((TClonesArray *) (*fRecHits3D)[id]);}
  inline TClonesArray  *RawClustAddress(Int_t id)  const{return ((TClonesArray *) (*fRawClusters)[id]);}    

  void DiagnosticsFE(Int_t evNumber1=0,Int_t evNumber2=0);    // Full events
  void DiagnosticsSE(Int_t diaglevel,Int_t evNumber1=0,Int_t evNumber2=0);    // Single events
 
  virtual void Print(Option_t *option)const; // Prints debug information
    
protected:
  TObjArray            *fChambers;           //! List of RICH chambers
  Int_t                 fNSDigits;           //Current number of sdigits
  Int_t                 fNcerenkovs;         //Current number of cerenkovs
  TClonesArray         *fSDigits;            //! List of sdigits
  TObjArray            *fDchambers;          //! Array of lists of digits
  TClonesArray         *fCerenkovs;          //! List of cerenkovs
  Int_t                 fNdch[kNCH];         //Array of current numbers of digits
  TObjArray            *fRawClusters;        // !List of raw clusters
  TObjArray            *fRecHits1D;          // !List of rec. hits
  TObjArray            *fRecHits3D;          // !List of rec. hits
  Int_t                 fNrawch[kNCH];       //Array of current numbers of raw clusters
  Int_t                 fNrechits1D[kNCH];   //Array of current numbers of rec hits 1D
  Int_t                 fNrechits3D[kNCH];   //Array of current numbers of rec hits 3D 

  Int_t fCkovNumber;                         // Number of Cerenkov photons
  Int_t fCkovQuarz;                          // Cerenkovs crossing quartz
  Int_t fCkovGap;                            // Cerenkovs crossing gap
  Int_t fCkovCsi;                            // Cerenkovs crossing csi
  Int_t fLostRfreo;                          // Cerenkovs reflected in freon
  Int_t fLostRquar;                          // Cerenkovs reflected in quartz
  Int_t fLostAfreo;                          // Cerenkovs absorbed in freon 
  Int_t fLostAquarz;                         // Cerenkovs absorbed in quartz
  Int_t fLostAmeta;                          // Cerenkovs absorbed in methane
  Int_t fLostCsi;                            // Cerenkovs below csi quantum efficiency 
  Int_t fLostWires;                          // Cerenkovs lost in wires
  Int_t fFreonProd;                          // Cerenkovs produced in freon
  Float_t fMipx;                             // x coord. of MIP
  Float_t fMipy;                             // y coord. of MIP
  Int_t fFeedbacks;                          // Number of feedback photons
  Int_t fLostFresnel;                        // Cerenkovs lost by Fresnel reflection

// Background event for event mixing
  Text_t *fFileName;                         //! File with background hits
  AliRICHMerger *fMerger;                    //! pointer to merger
    
  ClassDef(AliRICH,2)                        //Main RICH class 
};//class AliRICH
    
#endif
