#ifndef AliRICH_h
#define AliRICH_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////
//  Manager and hits classes for set:RICH     //
////////////////////////////////////////////////

#include <TObjArray.h>
#include <AliDetector.h>
#include <AliRICHConst.h>
#include <AliRICHChamber.h>
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
   
enum EDebugBits {kDebugStart=BIT(0),kDebugParam=BIT(1),kDebugHit=BIT(2),kDebugDigit=BIT(3),kDebugReco=BIT(4)}; // Debug flags
   
public:
             AliRICH();                                    // default ctor
             AliRICH(const char *name, const char *title); // named ctor
             AliRICH(const AliRICH& RICH);                 // copy ctor  
   virtual  ~AliRICH();                     // dtor
// Pure vituls methods     
   virtual Int_t  IsVersion() const =0;
// The following staff is defined in AliRICHChamber.cxx      
   virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits);
   virtual void   AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs);
   virtual void   AddSDigit(Int_t *clhits);
   virtual void   AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits);
   virtual void   AddRawCluster(Int_t id, const AliRICHRawCluster& cluster);
   virtual void   AddRecHit1D(Int_t id, Float_t* rechit, Float_t* photons, Int_t* padsx, Int_t* padsy);
   virtual void   AddRecHit3D(Int_t id, Float_t* rechit, Float_t omega, Float_t theta, Float_t phi);

   virtual void   BuildGeometry();   // TNode ROOT variant for event display
   virtual void   CreateGeometry();  // GEANT volumes tree for simulation
   virtual void   CreateMaterials(); // GEANT materials definition
   virtual Float_t AbsoCH4(Float_t x);
   virtual Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola);
   virtual void   StepManager();
   Int_t          DistancetoPrimitive(Int_t px, Int_t py);
   virtual void   MakeBranch(Option_t *opt=" ");
   virtual void   MakeBranchInTreeD(TTree *treeD, const char *file=0);
   void           SetTreeAddress();
   virtual void   ResetHits();
   virtual void   ResetDigits();
   virtual void   ResetRawClusters();
   virtual void   ResetRecHits1D();
   virtual void   ResetRecHits3D();
   virtual void   FindClusters(Int_t nev,Int_t lastEntry);
// Converters    
   virtual void   Hits2SDigits();
   virtual Int_t  Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss,Int_t id, ResponseType res);
   virtual void   SDigits2Digits();
   virtual void   SDigits2Digits(Int_t nev, Int_t flag);
   virtual void   Digits2Reco();
// Retrieve pad hits for a given Hit
    virtual        AliRICHSDigit* FirstPad(AliRICHHit *hit, TClonesArray *clusters);
    virtual        AliRICHSDigit* NextPad(TClonesArray *clusters);
// inline methods   
   TClonesArray  *SDigits()   const{return fSDigits;}
   TClonesArray  *Cerenkovs() const{return fCerenkovs;}

// Models for chambers
   virtual void     SetGeometryModel(Int_t iChamberN, AliRICHGeometry *pRICHGeo)    {       GetChamber(iChamberN)->SetGeometryModel(pRICHGeo);}
   AliRICHGeometry* GetGeometryModel(Int_t iChamberN=0)                        const{return GetChamber(iChamberN)->GetGeometryModel();}
    
   virtual void     SetSegmentationModel(Int_t iChamberN, AliSegmentation *pAliSeg) {       GetChamber(iChamberN)->SetSegmentationModel(pAliSeg);}
   AliSegmentation* GetSegmentationModel(Int_t iChamberN=0)                    const{return GetChamber(iChamberN)->GetSegmentationModel();}
	 
   virtual void     SetResponseModel(Int_t iChamberN, AliRICHResponse *pRICHRes)    {       GetChamber(iChamberN)->SetResponseModel(pRICHRes);}
   AliRICHResponse* GetResponseModel(Int_t iChamberN)                          const{return GetChamber(iChamberN)->GetResponseModel();}

   virtual void     SetReconstructionModel(Int_t iChamberN, AliRICHClusterFinder *pRICHReco){GetChamber(iChamberN)->SetReconstructionModel(pRICHReco);}
// Debug staff
   void   SetDebugLevel(Int_t level) {fDebugLevel=level;}
   Int_t  GetDebugLevel()       const{return fDebugLevel;}
   
   void   SetDebugStart()     {fDebugLevel+=kDebugStart;}        // Controls debug message at the entring point of methods
   void ResetDebugStart()     {fDebugLevel-=kDebugStart;}        // Controls debug message at the entring point of methods
   Bool_t  IsDebugStart()const{return fDebugLevel&kDebugStart;}  // Controls debug message at the entring point of methods
   
   void   SetDebugParam()     {fDebugLevel+=kDebugParam;}        // Controls debug printout for the parameters
   void ResetDebugParam()     {fDebugLevel-=kDebugParam;}        // Controls debug printout for the parameters
   Bool_t  IsDebugParam()const{return fDebugLevel&kDebugParam;}  // Controls debug printout for the parameters
   
   void   SetDebugHit()       {fDebugLevel+=kDebugHit;}          // Controls debug printout for hits
   void ResetDebugHit()       {fDebugLevel-=kDebugHit;}          // Controls debug printout for hits
   Bool_t  IsDebugHit()  const{return fDebugLevel&kDebugHit;}    // Controls debug printout for hits
   
   void   SetDebugDigit()     {fDebugLevel+=kDebugDigit;}        // Controls debug printout for digits
   void ResetDebugDigit()     {fDebugLevel-=kDebugDigit;}        // Controls debug printout for digits
   Bool_t  IsDebugDigit()const{return fDebugLevel&kDebugDigit;}  // Controls debug printout for digits
   
   void   SetDebugReco()      {fDebugLevel+=kDebugReco;}         // Controls debug printout for reco
   void ResetDebugReco()      {fDebugLevel-=kDebugReco;}         // Controls debug printout for reco
   Bool_t  IsDebugReco() const{return fDebugLevel&kDebugReco;}   // Controls debug printout for reco
   

    virtual void   SetMerger(AliRICHMerger* thisMerger) {fMerger=thisMerger;}  
// Return reference to Chamber #id
    virtual        AliRICHChamber& Chamber(Int_t id) {return *((AliRICHChamber *) (*fChambers)[id]);}
    AliRICHChamber* GetChamber(Int_t iChamberN)     const{return (AliRICHChamber*) (*fChambers)[iChamberN];}
// Return pointers to digits 
    TObjArray     *Dchambers() {return fDchambers;}
    Int_t         *Ndch() {return fNdch;}
    virtual TClonesArray *DigitsAddress(Int_t id) {return ((TClonesArray *) (*fDchambers)[id]);}
// Return pointers to rec. hits
    TObjArray     *RecHits1D()                      const{return fRecHits1D;}
    Int_t         *Nrechits1D()                          {return fNrechits1D;} // returns array
    Int_t         *Nrechits3D()                          {return fNrechits3D;} // returns array
   virtual TClonesArray *RecHitsAddress1D(Int_t id) const{return ((TClonesArray *) (*fRecHits1D)[id]);}
      TObjArray     *RecHits3D()                    const{return fRecHits3D;}
   virtual TClonesArray *RecHitsAddress3D(Int_t id) const{return ((TClonesArray *) (*fRecHits3D)[id]);}
   virtual TClonesArray *RawClustAddress(Int_t id)  const{return ((TClonesArray *) (*fRawClusters)[id]);}    
// Assignment operator
   AliRICH& operator=(const AliRICH& rhs);

   virtual void DiagnosticsFE(Int_t evNumber1=0,Int_t evNumber2=0);    // Full events
   virtual void DiagnosticsSE(Int_t diaglevel,Int_t evNumber1=0,Int_t evNumber2=0);    // Single events
 
   inline virtual void Print(Option_t *option)const; // Prints debug information
    
protected:
   TObjArray            *fChambers;           // !List of RICH chambers aka modules
   Int_t                 fNSDigits;           // Number of clusters
   Int_t                 fNcerenkovs;         // Number of cerenkovs
   TClonesArray         *fSDigits;            // !List of clusters
   TObjArray            *fDchambers;          // !List of digits
   TClonesArray         *fCerenkovs;          // !List of cerenkovs
   Int_t                 fNdch[kNCH];         // Number of digits
   TObjArray            *fRawClusters;        // !List of raw clusters
   TObjArray            *fRecHits1D;          // !List of rec. hits
   TObjArray            *fRecHits3D;          // !List of rec. hits
   Int_t                 fNrawch[kNCH];       // Number of raw clusters
   Int_t                 fNrechits1D[kNCH];   // Number of rec hits 
   Int_t                 fNrechits3D[kNCH];   // Number of rec hits 
   Int_t                 fDebugLevel;         // Source debugging level

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
    
   ClassDef(AliRICH,2)                        // Main RICH class 
};//class AliRICH
    
inline void AliRICH::Print(Option_t *option)const
{
   TObject::Print(option);
   if(IsDebugParam()){
      GetGeometryModel(0)->Print(option);
      GetSegmentationModel(0)->Print(option);
      GetResponseModel(0)->Print(option);
   }
}//inline void AliRICH::Print(Option_t *option)const

#endif
