#ifndef THBTPROCESSOR_H
#define THBTPROCESSOR_H
//_____________________________________________________________________________
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class THBTprocessor                                                       //
//                                                                           //
// Wrapper class to HBT processor fortran code.                              //
// For more information see AliGenHBTprocessor class                         //
// HBT processor is written by Lanny Ray                                     //
//                                                                           //
//                                                                           //
// Comunication is done via COMMON BLOCKS declared in HBTprocCOMMON.h        //
// using cfortran.h routines                                                 //
// User should use class AliGenHBTprocessor and all their interface          //
// see there for more description                                            //
//                                                                           //
// Wrapper class written by                                                  //
// Piotr Krzysztof Skowronski (Piotr.Skowronski@cern.ch)                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
/*******************************************************/
//Author Piotr Krzysztof Skowronski e-mial: Piotr.Skowronski@cern.ch
// C++ Wrapper Class for fortran made HBT Processor written by Lanny Ray
// 
// Comunication is done via COMMON BLOCKS declared in HBTprocCOMMON.h
// using cfortran.h routines
// User should use class AliGenHBTprocessor and all their interface
// see there for more description
//

#include <TGenerator.h>

class THBTprocessor: public TGenerator
  {
    public:
      THBTprocessor(); 
      virtual ~THBTprocessor() {};
      
      virtual void  Initialize() const;
      virtual void  GenerateEvent() const;
      virtual void  PrintEvent() const;
      virtual Int_t ImportParticles(TClonesArray *particles, Option_t *option="");
      virtual TObjArray      *ImportParticles(Option_t *option="");

      //Set/Gets
      //comprehensive description off all these methods 
      //can be found in proper methods of AliGenHBTprocessor
      
      virtual void SetTrackRejectionFactor(Float_t trf = 1.0) const;
      virtual void SetRefControl(Int_t rc =2) const;
      virtual void SetPIDs(Int_t pid1 = 8,Int_t pid2 = 9) const;
      virtual void SetNPIDtypes(Int_t npidt = 2)const;
      virtual void SetDeltap(Float_t deltp = 0.1) const;
      virtual void SetMaxIterations(Int_t maxiter = 50) const;
      virtual void SetDelChi(Float_t dc = 0.1)const;
      virtual void SetIRand(Int_t irnd = 76564) const;
      virtual void SetLambda(Float_t lam = 0.6) const;
      virtual void SetR1d(Float_t r = 7.0) const;
      virtual void SetRSide(Float_t rs = 6.0) const;
      virtual void SetROut(Float_t ro = 7.0) const ;
      virtual void SetRLong(Float_t rl = 4.0) const;
      virtual void SetRPerp(Float_t rp = 6.0) const;
      virtual void SetRParallel(Float_t rprl = 4.0) const;
      virtual void SetR0(Float_t r0 = 4.0) const;
      virtual void SetQ0(Float_t q0 = 9.0) const ;
      virtual void SetSwitch1D(Int_t s1d = 3) const ;
      virtual void SetSwitch3D(Int_t s3d = 0) const ;
      virtual void SetSwitchType(Int_t st = 3) const ;
      virtual void SetSwitchCoherence(Int_t sc = 0) const ;
      virtual void SetSwitchCoulomb(Int_t scol = 2) const ;
      virtual void SetSwitchFermiBose(Int_t sfb = 1) const ;
      
      virtual void SetPtRange(Float_t ptmin = 0.1, Float_t ptmax = 0.98)const ; //Pt in GeV/c
      virtual void SetPxRange(Float_t pxmin = -1.0, Float_t pxmax = 1.0) const;
      virtual void SetPyRange(Float_t pymin = -1.0, Float_t pymax = 1.0)const ;
      virtual void SetPzRange(Float_t pzmin = -3.6, Float_t pzmax = 3.6) const;
      virtual void SetPhiRange(Float_t phimin = 0.0, Float_t phimax = 360.0)const; //Angle in degrees
      virtual void SetEtaRange(Float_t etamin = -1.5, Float_t etamax = 1.5)const; //Pseudorapidity  !!!!!!!!!
      virtual void SetNPtBins(Int_t nptbin = 50)const;
      virtual void SetNPhiBins(Int_t nphibin = 50)const;
      virtual void SetNEtaBins(Int_t netabin = 50)const;
      virtual void SetNPxBins(Int_t npxbin = 20)const;
      virtual void SetNPyBins(Int_t npybin = 20)const;
      virtual void SetNPzBins(Int_t npzbin = 70)const;
      virtual void SetNBins1DFineMesh(Int_t n = 10)const;
      virtual void SetBinSize1DFineMesh(Float_t x=0.01)const;
      virtual void SetNBins1DCoarseMesh(Int_t n =2 )const;
      virtual void SetBinSize1DCoarseMesh(Float_t x=0.05)const;
      virtual void SetNBins3DFineMesh(Int_t n = 8)const;
      virtual void SetBinSize3DFineMesh(Float_t x=0.01)const;
      virtual void SetNBins3DCoarseMesh(Int_t n = 2)const;
      virtual void SetBinSize3DCoarseMesh(Float_t x=0.08)const;
      virtual void SetNBins3DFineProjectMesh(Int_t n =3 )const;
      virtual void SetPrintFull(Int_t flag = 1)const;
      
      void         DumpSettings() const;

    ClassDef(THBTprocessor,0)
  };

#endif

/**************************************************************
COMMENTS APROPOS HBTPROCESSOR FORTRAN CODE:

SUBROUTINE INITIALIZE HAS TO BE OMITTED, IT RESETS (ZEROES) COMMON BLOCKS

SUBROUTINE KINEMATICS HAS TO BE RENAMED EITHER IN MEVSIM OR IN PROCESSOR 
AS THEY HAVE IDENTICAL DECLARATIONS AND SOMTIMES IS USED THE BAD ONE IN CASE BOTH LIBRARIES ARE LOADED

****************************************************************/
