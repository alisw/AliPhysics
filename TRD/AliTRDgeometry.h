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

#include "TObjArray.h"

#include "AliGeometry.h"

class TGeoHMatrix;

class AliTRDpadPlane;

class AliTRDgeometry : public AliGeometry {

 public:

  enum { kNlayer  =   6
       , kNstack  =   5
       , kNsector =  18
       , kNdet    = 540 
       , kNdets   =  30 };

  AliTRDgeometry();
  AliTRDgeometry(const AliTRDgeometry &g);
  virtual ~AliTRDgeometry();
  AliTRDgeometry &operator=(const AliTRDgeometry &g);

  virtual void     Init();
  virtual void     CreateGeometry(Int_t *idtmed);
  virtual Int_t    IsVersion()                                            { return 1;               }
  virtual Bool_t   Impact(const TParticle* ) const                        { return kTRUE;           }
  virtual Bool_t   IsHole(Int_t la, Int_t st, Int_t se) const;
  virtual Bool_t   IsOnBoundary(Int_t det, Float_t y, Float_t z, Float_t eps = .5) const;
  virtual Bool_t   RotateBack(Int_t det, Double_t *loc, Double_t *glb) const;

          Bool_t   ChamberInGeometry(Int_t det);

          void     AssembleChamber(Int_t ilayer, Int_t istack);
          void     CreateFrame(Int_t *idtmed);
          void     CreateServices(Int_t *idtmed);

          Bool_t   CreateClusterMatrixArray();  
  TGeoHMatrix     *GetClusterMatrix(Int_t det)                           { return (TGeoHMatrix *) 
                                                                             fClusterMatrixArray->At(det); }

          void     SetSMstatus(Int_t sm, Char_t status)                  { fSMstatus[sm] = status;         }

  static  Int_t    GetDetectorSec(Int_t layer, Int_t stack);
  static  Int_t    GetDetector(Int_t layer, Int_t stack, Int_t sector);
  static  Int_t    GetLayer(Int_t det);
  static  Int_t    GetStack(Int_t det);
          Int_t    GetStack(Double_t z, Int_t layer);
  static  Int_t    GetSector(Int_t det);

          void     CreatePadPlaneArray();
  AliTRDpadPlane  *CreatePadPlane(Int_t layer, Int_t stack);
  AliTRDpadPlane  *GetPadPlane(Int_t layer, Int_t stack);
  AliTRDpadPlane  *GetPadPlane(Int_t det)                                { return GetPadPlane(GetLayer(det)
                                                                                             ,GetStack(det)); }
          Int_t    GetRowMax(Int_t layer, Int_t stack, Int_t /*sector*/);
          Int_t    GetColMax(Int_t layer);
          Double_t GetRow0(Int_t layer, Int_t stack, Int_t /*sector*/);
          Double_t GetCol0(Int_t layer);

  static  Float_t  GetTime0(Int_t layer)                                 { return fgkTime0[layer];        }

  static  Double_t GetXtrdBeg()                                          { return fgkXtrdBeg;             }
  static  Double_t GetXtrdEnd()                                          { return fgkXtrdEnd;             }

          Char_t   GetSMstatus(Int_t sm) const                           { return fSMstatus[sm];          }
          Float_t  GetChamberWidth(Int_t layer) const                    { return fCwidth[layer]      ;   }
          Float_t  GetChamberLength(Int_t layer, Int_t stack) const      { return fClength[layer][stack]; }

  virtual void     GetGlobal(const AliRecPoint*, TVector3&, TMatrixF& ) const { }; 
  virtual void     GetGlobal(const AliRecPoint*, TVector3& ) const            { };

  static  Double_t GetAlpha()                                            { return 2.0 
                                                                             * 3.14159265358979324 
                                                                             / fgkNsector;          } 

  static  Int_t    Nsector()                                             { return fgkNsector;       }
  static  Int_t    Nlayer()                                              { return fgkNlayer;        }
  static  Int_t    Nstack()                                              { return fgkNstack;        }
  static  Int_t    Ndet()                                                { return fgkNdet;          }

  static  Float_t  Cheight()                                             { return fgkCH;            }
  static  Float_t  CheightSV()                                           { return fgkCHsv;          }
  static  Float_t  Cspace()                                              { return fgkVspace;        }
  static  Float_t  CraHght()                                             { return fgkCraH;          }
  static  Float_t  CdrHght()                                             { return fgkCdrH;          }
  static  Float_t  CamHght()                                             { return fgkCamH;          }
  static  Float_t  CroHght()                                             { return fgkCroH;          }
  static  Float_t  CsvHght()                                             { return fgkCsvH;          }
  static  Float_t  CroWid()                                              { return fgkCroW;          }

  static  Float_t  AnodePos()                                            { return fgkAnodePos;      }

  static  Float_t  MyThick()                                             { return fgkRMyThick;      }
  static  Float_t  DrThick()                                             { return fgkDrThick;       }
  static  Float_t  AmThick()                                             { return fgkAmThick;       }
  static  Float_t  DrZpos()                                              { return fgkDrZpos;        }
  static  Float_t  RpadW()                                               { return fgkRpadW;         }
  static  Float_t  CpadW()                                               { return fgkCpadW;         }

  static  Float_t  Cwidcha()                                             { return (fgkSwidth2 - fgkSwidth1) 
                                                                                  / fgkSheight 
                                                                                  * (fgkCH + fgkVspace);      }

  static  Int_t    MCMmax()                                              { return fgkMCMmax;        }
  static  Int_t    MCMrow()                                              { return fgkMCMrow;        }
  static  Int_t    ROBmaxC0()                                            { return fgkROBmaxC0;      }
  static  Int_t    ROBmaxC1()                                            { return fgkROBmaxC1;      }
  static  Int_t    ADCmax()                                              { return fgkADCmax;        }
  static  Int_t    TBmax()                                               { return fgkTBmax;         }            
  static  Int_t    Padmax()                                              { return fgkPadmax;        }
  static  Int_t    Colmax()                                              { return fgkColmax;        }
  static  Int_t    RowmaxC0()                                            { return fgkRowmaxC0;      }
  static  Int_t    RowmaxC1()                                            { return fgkRowmaxC1;      }

 protected:

  static const Int_t    fgkNsector;                          //  Number of sectors in the full detector (18)
  static const Int_t    fgkNlayer;                           //  Number of layers of the TRD (6)
  static const Int_t    fgkNstack;                           //  Number of stacks in z-direction (5)
  static const Int_t    fgkNdet;                             //  Total number of detectors (18 * 6 * 5 = 540)

  static const Float_t  fgkTlength;                          //  Length of the TRD-volume in spaceframe (BTRD)

  static const Float_t  fgkSheight;                          //  Height of the supermodule
  static const Float_t  fgkSwidth1;                          //  Lower width of the supermodule
  static const Float_t  fgkSwidth2;                          //  Upper width of the supermodule
  static const Float_t  fgkSlength;                          //  Length of the supermodule

  static const Float_t  fgkFlength;                          //  Length of the service space in front of a supermodule

  static const Float_t  fgkSMpltT;                           //  Thickness of the super module side plates

  static const Float_t  fgkCraH;                             //  Height of the radiator part of the chambers
  static const Float_t  fgkCdrH;                             //  Height of the drift region of the chambers
  static const Float_t  fgkCamH;                             //  Height of the amplification region of the chambers
  static const Float_t  fgkCroH;                             //  Height of the readout of the chambers
  static const Float_t  fgkCsvH;                             //  Height of the services on top of the chambers
  static const Float_t  fgkCH;                               //  Total height of the chambers (w/o services)
  static const Float_t  fgkCHsv;                             //  Total height of the chambers (with services)

  static const Float_t  fgkAnodePos;                         //  Distance of anode wire plane relative to alignabl volume

  static const Float_t  fgkVspace;                           //  Vertical spacing of the chambers
  static const Float_t  fgkHspace;                           //  Horizontal spacing of the chambers
  static const Float_t  fgkVrocsm;                           //  Radial distance of the first ROC to the outer SM plates

  static const Float_t  fgkCalT;                             //  Thickness of the lower aluminum frame
  static const Float_t  fgkCalW;                             //  Width of additional aluminum ledge on lower frame
  static const Float_t  fgkCalH;                             //  Height of additional aluminum ledge on lower frame
  static const Float_t  fgkCalWmod;                          //  Width of additional aluminum ledge on lower frame
  static const Float_t  fgkCalHmod;                          //  Height of additional aluminum ledge on lower frame
  static const Float_t  fgkCwsW;                             //  Width of additional wacosit ledge on lower frame
  static const Float_t  fgkCwsH;                             //  Height of additional wacosit ledge on lower frame
  static const Float_t  fgkCclsT;                            //  Thickness of the lower Wacosit frame sides
  static const Float_t  fgkCclfT;                            //  Thickness of the lower Wacosit frame front
  static const Float_t  fgkCglT;                             //  Thichness of the glue around the radiator
  static const Float_t  fgkCcuTa;                            //  Thickness of the upper Wacosit frame around amp. region
  static const Float_t  fgkCcuTb;                            //  Thickness of the upper Wacosit frame around amp. region
  static const Float_t  fgkCauT;                             //  Thickness of the aluminum frame of the back panel
  static const Float_t  fgkCroW;                             //  Additional width of the readout chamber frames

  static const Float_t  fgkCpadW;                            //  Difference of outer chamber width and pad plane width
  static const Float_t  fgkRpadW;                            //  Difference of outer chamber width and pad plane width

  static const Float_t  fgkXeThick;                          //  Thickness of the gas volume
  static const Float_t  fgkDrThick;                          //  Thickness of the drift region
  static const Float_t  fgkAmThick;                          //  Thickness of the amplification region
  static const Float_t  fgkWrThick;                          //  Thickness of the wire planes

  static const Float_t  fgkPPdThick;                         //  Thickness of copper of the pad plane
  static const Float_t  fgkPPpThick;                         //  Thickness of PCB board of the pad plane
  static const Float_t  fgkPGlThick;                         //  Thickness of the glue layer
  static const Float_t  fgkPCbThick;                         //  Thickness of the carbon layers
  static const Float_t  fgkPHcThick;                         //  Thickness of the honeycomb support structure
  static const Float_t  fgkPPcThick;                         //  Thickness of the PCB readout boards
  static const Float_t  fgkPRbThick;                         //  Thickness of the PCB copper layers
  static const Float_t  fgkPElThick;                         //  Thickness of all other electronics components (caps, etc.)

  static const Float_t  fgkRFbThick;                         //  Thickness of the fiber layers in the radiator
  static const Float_t  fgkRRhThick;                         //  Thickness of the rohacell layers in the radiator
  static const Float_t  fgkRGlThick;                         //  Thickness of the glue layers in the radiator
  static const Float_t  fgkRCbThick;                         //  Thickness of the carbon layers in the radiator
  static const Float_t  fgkRMyThick;                         //  Thickness of the mylar layers in the radiator

  static const Float_t  fgkDrZpos;                           //  Position of the drift region
  static const Float_t  fgkAmZpos;                           //  Position of the amplification region
  static const Float_t  fgkWrZposA;                          //  Position of the wire planes
  static const Float_t  fgkWrZposB;                          //  Position of the wire planes
  static const Float_t  fgkCalZpos;                          //  Position of the additional aluminum ledges

  static const Int_t    fgkMCMmax;                           //  Maximum number of MCMs per ROB
  static const Int_t    fgkMCMrow;                           //  Maximum number of MCMs per ROB Row
  static const Int_t    fgkROBmaxC0;                         //  Maximum number of ROBs per C0 chamber
  static const Int_t    fgkROBmaxC1;                         //  Maximum number of ROBs per C1 chamber
  static const Int_t    fgkADCmax;                           //  Maximum number of ADC channels per MCM
  static const Int_t    fgkTBmax;                            //  Maximum number of Time bins
  static const Int_t    fgkPadmax;                           //  Maximum number of pads per MCM
  static const Int_t    fgkColmax;                           //  Maximum number of pads per padplane row
  static const Int_t    fgkRowmaxC0;                         //  Maximum number of Rows per C0 chamber
  static const Int_t    fgkRowmaxC1;                         //  Maximum number of Rows per C1 chamber

  Float_t               fCwidth[kNlayer];                    //  Outer widths of the chambers
  Float_t               fClength[kNlayer][kNstack];          //  Outer lengths of the chambers

  Float_t               fRotB11[kNsector];                   //  Matrix elements for the backward rotation
  Float_t               fRotB12[kNsector];                   //  Matrix elements for the backward rotation
  Float_t               fRotB21[kNsector];                   //  Matrix elements for the backward rotation
  Float_t               fRotB22[kNsector];                   //  Matrix elements for the backward rotation

  static const Double_t fgkTime0Base;                        //  Base value for calculation of Time-position of pad 0
  static const Float_t  fgkTime0[kNlayer];                   //  Time-position of pad 0

  static const Double_t fgkXtrdBeg;                          //  X-coordinate in tracking system of begin of TRD mother volume
  static const Double_t fgkXtrdEnd;                          //  X-coordinate in tracking system of end of TRD mother volume

  TObjArray            *fClusterMatrixArray;                 //! Transformation matrices loc. cluster to tracking cs
  TObjArray            *fPadPlaneArray;                      //! Array of pad plane objects

  Char_t                fSMstatus[kNsector];                 //  Super module status byte

  ClassDef(AliTRDgeometry,23)                                //  TRD geometry class

};
#endif
