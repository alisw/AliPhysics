#ifndef ALIITSV11GEOMETRYSPD_H
#define ALIITSV11GEOMETRYSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */
#include <AliITSv11Geometry.h>
class TGeoVolume;

class AliITSv11GeometrySPD : public AliITSv11Geometry {
  public:
    AliITSv11GeometrySPD();
    AliITSv11GeometrySPD(Int_t debug);
    virtual ~AliITSv11GeometrySPD(){};
    //
    virtual TGeoVolume* CenteralSPD(TGeoVolume *moth);
    //
    // Create figures for the documentation of this class
    virtual void CreateFigure0(const Char_t *filepath="",
                               const Char_t *type="gif");
    virtual void CreateFigureLadder(const Char_t*filepath="",
                                    const Char_t* type="gif");
    //
    virtual const char* GetSensitivevolumeName() const{
        //Get Sensitive volume name
        return fSPDSensitiveVolumeName.Data();}
    virtual void SetSensitivevolumeName(const Char_t 
                *n="ITSSPDDetectorSensitiveVolume"){//Set Sensitive volume name
                fSPDSensitiveVolumeName=n;}
  private:
    void SPDsectorShape(Int_t n,const Double_t *xc,const Double_t *yc,
          const Double_t *r,const Double_t *ths,const Double_t *the,Int_t npr,
                        Int_t &m,Double_t **xp,Double_t **yp);
    const char* CreateSensitivevolumeName(const char* app) const{
        //Get Sensitive volume name
        char *a = new char[fSPDSensitiveVolumeName.Length()+strlen(app)+1]; 
        strcpy(a,fSPDSensitiveVolumeName.Data());return strcat(a,app);}
    //
    TGeoVolume* CarbonFiberSector(TGeoVolume *moth);
    TGeoVolume* CreateHalfStaveGroundFoil();
    TGeoVolume* CreateSPDLadder();

  private:
    // Posible Variables
    TString fSPDSensitiveVolumeName; // ITS SPD sensitive volume name
    Double_t fThickDetector; // Detector thickness
    Double_t fThickChip;     // Read out Chip thickness
    // Stave
    // Grounding Foil
    const Double_t fkGrdFoilThick        ;//!  0.05*fgkmm;
    const Double_t fkGrdFoilWidthA       ;//! 15.95*fgkmm;
    const Double_t fkGrdFoilWidthC       ;//!  4.4*fgkmm;
    const Double_t fkGrdFoilLngA         ;//!139.89*fgkmm;
    const Double_t fkGrdFoilLngB         ;//! 11.55*fgkmm;
    const Double_t fkGrdFoilLngC         ;//! 82.0*fgkmm;
    const Int_t    fkGrdFoilNholesAB     ;//!5;
    const Double_t fkGrdFoilHoleCenterAB ;//!  7.8*fgkmm;
    const Double_t fkGrdFoilHoleLengthAB ;//! 12.0*fgkmm;
    const Double_t fkGrdFoilHoleWidthAB  ;//!  7.5*fgkmm;
    const Double_t fkGrdFoilHoleSpacingAB;//! 14.0*fgkmm;
    const Double_t fkGrdFoilHoleStartA   ;//!  1.36*fgkmm;
    const Double_t fkGrdFoilHoleStartB   ;//! 73.08*fgkmm;
    // Ladder
    const Int_t    fkLadNChips        ;//!5;
    const Double_t fkLadChipWidth     ;//!15950.0*fgkmicron;
    const Double_t fkLadChipHight     ;//!  150.0*fgkmicron;
    const Double_t fkLadChipLength    ;//!13490.0*fgkmicron;
    const Double_t fkLadGlue0Thick    ;//!    0.100*fgkmm;
    const Double_t fkLadBumpBondThick ;//!   30.0*fgkmicron;
    const Double_t fkLadDetectorWidth ;//!13700.0*fgkmicron;
    const Double_t fkLadDetectorThick ;//!  200.0*fgkmicron;
    const Double_t fkLadDetectorLength;//!70710.0*fgkmicron;
    const Double_t fkLadSensDetWidth  ;//! 1280.0*fgkmicron;
    const Double_t fkLadSensDetThick  ;//!  200.0*fgkmicron;
    const Double_t fkLadSensDetLength ;//!69490.0*fgkmicron;
    const Double_t fkLadChipSpacing0  ;//!  610.0*fgkmicron;
    const Double_t fkLadChipSpacing1  ;//! (2.*fkLadChipSpacing0+
    //    ((Double_t)fkLadNChips)*fkLadChipLength)/((Double_t)(fkLadNChips-1));
    //

    ClassDef(AliITSv11GeometrySPD,1) // ITS v11 Centeral SPD geometry
};

#endif
