#ifndef ALIEVE_TOFSector_H
#define ALIEVE_TOFSector_H

#include <Reve/QuadSet.h>
#include <Reve/RenderElement.h>

#include <Reve/RGBAPalette.h>
#include <Reve/FrameBox.h>

#include <TGeoManager.h>
#include <TClonesArray.h>
#include <TTree.h>

#include <AliTOFGeometry.h>

namespace Alieve {
  
  class TOFSector : public Reve::QuadSet
                   
  {
    TOFSector(const TOFSector&);            // Not implemented
    TOFSector& operator=(const TOFSector&); // Not implemented
 
    //Int_t       fSectorID;
  private:

    void LoadQuads();
    
  protected:
    
    AliTOFGeometry *fTOFgeometry;
    
    TClonesArray   *fTOFarray;
    TTree          *fTOFtree;
    
    Int_t fSector;
    //Int_t fPlate;
    //Int_t fStrip;
    
    Float_t  fDx;
    Float_t  fDy;
    Float_t  fDz;
    ///////////////////////////////

    Bool_t      fAutoTrans; 
    Int_t       fMinTime;     
    Int_t       fMaxTime;
    Short_t     fThreshold;
    Int_t       fMaxVal;
    Int_t       fSectorID; 
    Bool_t      *fPlateFlag;

    Bool_t      fPlateFlag0;
    Bool_t      fPlateFlag1;
    Bool_t      fPlateFlag2;
    Bool_t      fPlateFlag3;
    Bool_t      fPlateFlag4;
    
    Color_t     fFrameColor;
    Bool_t      fRnrFrame;
    
    TGeoManager *fGeoManager;
    
  public: 
    // Bool_t       fAutoTrans;
    
    virtual void InitModule();
    virtual void SetTrans(); 
    TOFSector(const Text_t* n="TOFSector", const Text_t* t=0);
    TOFSector(TGeoManager *localGeoManager, Int_t nSector);
    
    TOFSector(TGeoManager *localGeoManager, Int_t nSector,
	      TClonesArray *tofArray);
    TOFSector(TGeoManager *localGeoManager,
	      Int_t nSector, TTree *tofTree);
    virtual ~TOFSector();
        
    static Bool_t    fgStaticInitDone;
    static void      InitStatics();

    void SetSectorID(Int_t id);
    void SetAutoTrans(Bool_t r){fAutoTrans=r;};
    void SetThreshold(Short_t t);
    void SetMaxVal(Int_t mv);
    Bool_t GetPlate(Int_t nPlate) const {return fPlateFlag[nPlate];};
    Short_t GetThreshold() const {return fThreshold;};
    Int_t GetMaxVal() const {return fMaxVal;};
    Bool_t GetAutoTrans() const {return fAutoTrans;};
    Int_t GetSectorID() const {return fSectorID;};
    virtual void DigitSelected(Int_t idx);
    ///////////////////////////////////////////
   
    void SetPlate(Int_t nPlate, Bool_t r);

    static Reve::FrameBox    *fgTOFsectorFrameBox;
    static Reve::RGBAPalette *fgTOFsectorPalette;

  ClassDef(TOFSector, 1);
  }; 
}
#endif
