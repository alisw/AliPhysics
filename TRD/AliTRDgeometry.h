#ifndef ALITRDGEOMETRY_H
#define ALITRDGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD geometry class                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliGeometry.h"

#include "TObjArray.h"

class AliRunLoader;
class TGeoHMatrix;

class AliTRDgeometry : public AliGeometry {

 public:

          enum { kNplan =   6
               , kNcham =   5
               , kNsect =  18
               , kNdet  = 540 
               , kNdets =  30 };

  AliTRDgeometry();
  AliTRDgeometry(const AliTRDgeometry &g);
  virtual ~AliTRDgeometry();
  AliTRDgeometry &operator=(const AliTRDgeometry &g);

  virtual void     Init();
  virtual void     CreateGeometry(Int_t *idtmed);
  virtual Int_t    IsVersion()                                         { return 1;               }
  virtual Bool_t   Impact(const TParticle* ) const                     { return kTRUE;           }
  virtual Bool_t   IsHole(Int_t /*p*/, Int_t /*c*/, Int_t /*s*/) const { return kFALSE;          }

  virtual Bool_t   Rotate(Int_t d, Double_t *pos, Double_t *rot) const;
  virtual Bool_t   RotateBack(Int_t d, Double_t *rot, Double_t *pos) const;

          void     GroupChamber(Int_t iplan, Int_t icham, Int_t *idtmed);
          void     CreateFrame(Int_t *idtmed);
          void     CreateServices(Int_t *idtmed);

          Bool_t   ReadGeoMatrices();  

          void     SetSMstatus(Int_t sm, Char_t status)                { fSMstatus[sm] = status; }

  static  AliTRDgeometry* GetGeometry(AliRunLoader *runLoader = NULL);
  
  static  Int_t    GetDetectorSec(Int_t p, Int_t c);
  static  Int_t    GetDetector(Int_t p, Int_t c, Int_t s);
  virtual Int_t    GetPlane(Int_t d) const;
  virtual Int_t    GetChamber(Int_t d) const;
  virtual Int_t    GetSector(Int_t d) const;

  static Float_t   GetTime0(Int_t p)                                   { return fgkTime0[p];     }

          Char_t   GetSMstatus(Int_t sm) const                         { return fSMstatus[sm];   }
          Float_t  GetChamberWidth(Int_t p) const                      { return fCwidth[p];      }
          Float_t  GetChamberLength(Int_t p, Int_t c) const            { return fClength[p][c];  }

  virtual void     GetGlobal(const AliRecPoint*, TVector3&, TMatrixF& ) const { }; 
  virtual void     GetGlobal(const AliRecPoint*, TVector3& ) const            { };
 
  static  Double_t GetAlpha()                                          { return 2.0 
                                                                           * 3.14159265358979324 
                                                                           / fgkNsect;           } 

  static  Int_t    Nsect()                                             { return fgkNsect;        }
  static  Int_t    Nplan()                                             { return fgkNplan;        }
  static  Int_t    Ncham()                                             { return fgkNcham;        }
  static  Int_t    Ndet()                                              { return fgkNdet;         }

  static  Float_t  Cheight()                                           { return fgkCH;           }
  static  Float_t  Cspace()                                            { return fgkVspace;       }
  static  Float_t  CraHght()                                           { return fgkCraH;         }
  static  Float_t  CdrHght()                                           { return fgkCdrH;         }
  static  Float_t  CamHght()                                           { return fgkCamH;         }
  static  Float_t  CroHght()                                           { return fgkCroH;         }
  static  Float_t  CroWid()                                            { return fgkCroW;         }
  static  Float_t  MyThick()                                           { return fgkMyThick;      }
  static  Float_t  DrThick()                                           { return fgkDrThick;      }
  static  Float_t  AmThick()                                           { return fgkAmThick;      }
  static  Float_t  DrZpos()                                            { return fgkDrZpos;       }
  static  Float_t  RpadW()                                             { return fgkRpadW;        }
  static  Float_t  CpadW()                                             { return fgkCpadW;        }

  static  Float_t  Cwidcha()                                           { return (fgkSwidth2 - fgkSwidth1) 
                                                                                / fgkSheight 
                                                                                * (fgkCH + fgkVspace);      }

  TGeoHMatrix     *GetGeoMatrix(Int_t det)                             { return (TGeoHMatrix *) 
                                                                           fMatrixGeo->At(det);             }
  TGeoHMatrix     *GetMatrix(Int_t det)                                { return (TGeoHMatrix *) 
                                                                           fMatrixArray->At(det);           }
  TGeoHMatrix     *GetCorrectionMatrix(Int_t det)                      { return (TGeoHMatrix *) 
                                                                           fMatrixCorrectionArray->At(det); }

 protected:

  static const Int_t    fgkNsect;                            //  Number of sectors in the full detector (18)
  static const Int_t    fgkNplan;                            //  Number of planes of the TRD (6)
  static const Int_t    fgkNcham;                            //  Number of chambers in z-direction (5)
  static const Int_t    fgkNdet;                             //  Total number of detectors (18 * 6 * 5 = 540)

  static const Float_t  fgkSheight;                          //  Height of the TRD-volume in spaceframe (BTRD)
  static const Float_t  fgkSwidth1;                          //  Lower width of the TRD-volume in spaceframe (BTRD)
  static const Float_t  fgkSwidth2;                          //  Upper width of the TRD-volume in spaceframe (BTRD)
  static const Float_t  fgkSlength;                          //  Length of the TRD-volume in spaceframe (BTRD)

  static const Float_t  fgkSMpltT;                           //  Thickness of the super module side plates

  static const Float_t  fgkCraH;                             //  Height of the radiator part of the chambers
  static const Float_t  fgkCdrH;                             //  Height of the drift region of the chambers
  static const Float_t  fgkCamH;                             //  Height of the amplification region of the chambers
  static const Float_t  fgkCroH;                             //  Height of the readout of the chambers
  static const Float_t  fgkCH;                               //  Total height of the chambers

  static const Float_t  fgkVspace;                           //  Vertical spacing of the chambers
  static const Float_t  fgkHspace;                           //  Horizontal spacing of the chambers
  static const Float_t  fgkVrocsm;                           //  Radial distance of the first ROC to the outer SM plates
  static const Float_t  fgkCalT;                             //  Thickness of the lower aluminum frame
  static const Float_t  fgkCalW;                             //  Width of additional aluminum on lower frame
  static const Float_t  fgkCclsT;                            //  Thickness of the lower Wacosit frame sides
  static const Float_t  fgkCclfT;                            //  Thickness of the lower Wacosit frame front
  static const Float_t  fgkCglT;                             //  Thichness of the glue around the radiator
  static const Float_t  fgkCcuT;                             //  Thickness of the upper Wacosit frame
  static const Float_t  fgkCauT;                             //  Thickness of the aluminum frame of the back panel

  static const Float_t  fgkCroW;                             //  Additional width of the readout chamber frames

  static const Float_t  fgkCpadW;                            //  Difference of outer chamber width and pad plane width
  static const Float_t  fgkRpadW;                            //  Difference of outer chamber width and pad plane width

  static const Float_t  fgkMyThick;                          //  Thickness of the mylar-layer
  static const Float_t  fgkRaThick;                          //  Thickness of the radiator
  static const Float_t  fgkXeThick;                          //  Thickness of the gas volume
  static const Float_t  fgkDrThick;                          //  Thickness of the drift region
  static const Float_t  fgkAmThick;                          //  Thickness of the amplification region
  static const Float_t  fgkWrThick;                          //  Thickness of the wire planes
  static const Float_t  fgkCuThick;                          //  Thickness of the pad plane
  static const Float_t  fgkGlThick;                          //  Thickness of the glue layer
  static const Float_t  fgkSuThick;                          //  Thickness of the NOMEX support structure
  static const Float_t  fgkRpThick;                          //  Thickness of the PCB readout boards
  static const Float_t  fgkRcThick;                          //  Thickness of the PCB copper layers
  static const Float_t  fgkRoThick;                          //  Thickness of all other ROB componentes (caps, etc.)

  static const Float_t  fgkRaZpos;                           //  Position of the radiator
  static const Float_t  fgkDrZpos;                           //  Position of the drift region
  static const Float_t  fgkAmZpos;                           //  Position of the amplification region
  static const Float_t  fgkWrZpos;                           //  Position of the wire planes
  static const Float_t  fgkCuZpos;                           //  Position of the pad plane
  static const Float_t  fgkGlZpos;                           //  Position of the glue layer
  static const Float_t  fgkSuZpos;                           //  Position of the HEXCEL+G10 support structure
  static const Float_t  fgkRpZpos;                           //  Position of the PCB readout boards
  static const Float_t  fgkRcZpos;                           //  Position of the PCB copper layers
  static const Float_t  fgkRoZpos;                           //  Position of all other ROB componentes (caps, etc.)

  Char_t                fSMstatus[kNsect];                   //  Super module status byte

  Float_t               fCwidth[kNplan];                     //  Outer widths of the chambers
  Float_t               fClength[kNplan][kNcham];            //  Outer lengths of the chambers

  Float_t               fRotA11[kNsect];                     //  Matrix elements for the rotation
  Float_t               fRotA12[kNsect];                     //  Matrix elements for the rotation
  Float_t               fRotA21[kNsect];                     //  Matrix elements for the rotation
  Float_t               fRotA22[kNsect];                     //  Matrix elements for the rotation

  Float_t               fRotB11[kNsect];                     //  Matrix elements for the backward rotation
  Float_t               fRotB12[kNsect];                     //  Matrix elements for the backward rotation
  Float_t               fRotB21[kNsect];                     //  Matrix elements for the backward rotation
  Float_t               fRotB22[kNsect];                     //  Matrix elements for the backward rotation

  static const Double_t fgkTime0Base;                        //  Base value for calculation of Time-position of pad 0
  static const Float_t  fgkTime0[kNplan];                    //  Time-position of pad 0
  
  Float_t               fChamberUAorig[3*kNdets][3];         //  Volumes origin in
  Float_t               fChamberUDorig[3*kNdets][3];         //  the chamber
  Float_t               fChamberUForig[3*kNdets][3];         //  [3] = x, y, z
  Float_t               fChamberUUorig[3*kNdets][3];         //

  Float_t               fChamberUAboxd[3*kNdets][3];         //  Volumes box
  Float_t               fChamberUDboxd[3*kNdets][3];         //  dimensions (half)
  Float_t               fChamberUFboxd[3*kNdets][3];         //  [3] = x, y, z
  Float_t               fChamberUUboxd[3*kNdets][3];         // 

  TObjArray *           fMatrixArray;                        //! Transformation Global to Local
  TObjArray *           fMatrixCorrectionArray;              //! Transformation Cluster to  Tracking systerm
  TObjArray *           fMatrixGeo;                          //! Geo matrices

  ClassDef(AliTRDgeometry,11)                                //  TRD geometry class

};

#endif
