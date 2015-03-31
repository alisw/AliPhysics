/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// Revision of includes 07/05/2004
//
/// \ingroup sim
/// \class AliMUONTriggerGeometryBuilder
/// \brief MUON Trigger stations geometry construction class
///
/// \author Philippe Crochet (LPCCFd)

#ifndef ALI_MUON_TRIGGER_GEOMETRY_BUILDER_H
#define ALI_MUON_TRIGGER_GEOMETRY_BUILDER_H

#include "AliMUONVGeometryBuilder.h"

class AliMUON;

class AliMUONTriggerGeometryBuilder : public AliMUONVGeometryBuilder
{
  public:
    AliMUONTriggerGeometryBuilder(AliMUON* muon);
    AliMUONTriggerGeometryBuilder();
    virtual ~AliMUONTriggerGeometryBuilder();
  
    // methods
    virtual void CreateGeometry();
    virtual void SetVolumes();
    virtual void SetTransformations();
    virtual void SetSensitiveVolumes();
    
    /// Do not apply global transformation (geometry is defined in the new ALICE reference frame)
    virtual bool ApplyGlobalTransformation() { return false; }
    
  protected:  
    
  private:
    /// Not implemented
    AliMUONTriggerGeometryBuilder(const AliMUONTriggerGeometryBuilder& rhs);
    /// Not implemented
    AliMUONTriggerGeometryBuilder& operator = (const AliMUONTriggerGeometryBuilder& rhs);

    // methods
    void BuildChamberPrototype(Int_t icount) const;  
    void BuildRPCSupportsVertical(Int_t& iVolNum, Int_t icount) const;    
    void BuildRPCSupportsHorizontal(Int_t icount) const;    
    void BuildAngularSupportForChambers(Int_t icount) const;  
    void BuildGasPipes(Int_t icount) const;         
    void BuildChamberTypeA(Int_t& iVolNum, Int_t icount);  
    void BuildChamberTypeB(Int_t& iVolNum, Int_t icount);  
    void BuildChamberTypeD(Int_t& iVolNum, Int_t icount);  
    void BuildChamberTypeE(Int_t& iVolNum, Int_t icount);  
    void BuildChamberTypeF(Int_t& iVolNum, Int_t icount);  

    // constants
    
    static const Float_t fgkDXZERO; ///<  vertical gap between right and left chambers (kDXZERO*2=4cm)

    // main distances for chamber definition in first plane/first station
    static const Float_t fgkXMIN; ///< xmin distance in first plane/first station    
    static const Float_t fgkXMED; ///< xmed distance in first plane/first station                                 
    static const Float_t fgkXMAX; ///< xmax distance in first plane/first station  

    // 090704 kXMAX changed from 272 to 255.
    // (see fig.2-4 & 2-5 of Local Trigger Board PRR)
    // segmentation updated accordingly
    static const Float_t fgkYMIN; ///< add                             
    static const Float_t fgkYMAX; ///< add                            

    // inner/outer radius of flange between beam shield. and chambers (1/station)
    //    static const Float_t fgkRMIN[2]={50.,50.};
    //    static const Float_t fgkRMAX[2]={64.,68.};
    // z position of the middle of the gas gap in mother vol 
    static const Float_t fgkZm; ///< inner radius of flange between beam shield. and chambers (1/station)
    static const Float_t fgkZp; ///< outer radius of flange between beam shield. and chambers (1/station)
    
    static const Float_t fgkYVSup[4]; ///< y positions of vertical supports

    static const Float_t fgkSizeVSupExt[3]; ///< ext dimensions of vertical supports 
    static const Float_t fgkSizeVSupInt[3]; ///< int dimensions of vertical supports  

    static const Float_t fgkSizeSupport1V[3]; ///< transverse dimensions of 1V angular supports 
    static const Float_t fgkSizeSupport1H[3]; ///< transverse dimensions of 1H angular supports 
       // z should be 1.4 in the installed set-up 
    static const Float_t fgkSizeSupport2V[3]; ///< transverse dimensions of 2V angular supports  
    static const Float_t fgkSizeSupport2H[3]; ///< transverse dimensions of 2H angular supports  
    static const Float_t fgkSizeSupportXV[3]; ///< transverse dimensions of XV angular supports 
    static const Float_t fgkSizeSupportXH[3]; ///< transverse dimensions of XH angular supports 

    static const Float_t fgkSizeSupportCable[3]; /// transverse dimensions of horizontal cable supports
    static const Float_t fgkSizeGasPipe[3]; ///< dimensions of gas pipes (inner and outer radius)

    static const Float_t fgkOffsetGasPipe;  ///< Position of gas pipe with respect to angular support
    static const Float_t fgkAvoidExtrusion; ///<  Small cut on some volumes to avoid extrusion from SC1x
    
    // 
    TString GetVolumeName(const TString& volume, Int_t icount) const; 
    TString GetVolEnvName(Int_t icount, Int_t ienv) const; 
    TString GetVolAluAngSuppName(
                     const TString& type1234X, 
                     const TString& typeHV,
                     Int_t icount) const;                      
    TString GetVolEnvSuppAngName(
                     const TString& type1234X, 
                     const TString& typeHV, 
                     const TString& typeABDEF,
                     Int_t icount, Int_t ivol) const;                      
    TString GetVolEnvInoxGasPipeName(
                     const TString& type12, 
                     const TString& typeABCDEF,
                     Int_t icount, Int_t ivol) const;                      
                           
    
    // data members 
    AliMUON*  fMUON;   ///< the MUON detector class 
    Int_t*    fIdtmed; //!<! tracking media   
    Int_t     fIdAir;  //!<! medium 1
    Int_t     fIdAlu1; //!<! medium 4
    Int_t     fIdInox; //!<! medium 29 Stainless Steel (18%Cr,9%Ni,Fe)
    Float_t   fYEnvPsave; //!<! add
    Float_t   fYEnvMsave; //!<! add
    Float_t   fDYsave;    //!<! add
    Float_t   fDXsave;    //!<! add
    TGeoRotation fRsupportpipe; ///< pipe support rotation 
        
  ClassDef(AliMUONTriggerGeometryBuilder,2) // MUON Trigger stations geometry construction class
};

#endif //ALI_MUON_TRIGGER_GEOMETRY_BUILDER_H
