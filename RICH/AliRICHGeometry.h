#ifndef AliRICHGeometry_h
#define AliRICHGeometry_h


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <iostream.h>
#include <TObject.h>

class AliRICHGeometry :public TObject  
{
 public:
   inline AliRICHGeometry(); // default ctor
   inline void Print();
          
   void    SetGapThickness(Float_t thickness)          {fGapThickness=thickness;}    // Radiator Thickness 
   Float_t GetGapThickness()                      const{return fGapThickness;}       // Radiator thickness

   void    SetProximityGapThickness(Float_t thickness) {fProximityGapThickness=thickness;}
   Float_t GetProximityGapThickness()             const{return fProximityGapThickness;}
    
   void    SetQuartzLength(Float_t length)             {fQuartzLength=length;}
   Float_t GetQuartzLength()                      const{return fQuartzLength;}
   
   void    SetQuartzWidth(Float_t width)               {fQuartzWidth=width;}
   Float_t GetQuartzWidth()                       const{return fQuartzWidth;}

   void    SetQuartzThickness(Float_t thickness)       {fQuartzThickness=thickness;}
   Float_t GetQuartzThickness()                   const{return fQuartzThickness;}
   
   void    SetOuterFreonLength(Float_t length)         {fOuterFreonLength=length;}
   Float_t GetOuterFreonLength()                  const{return fOuterFreonLength;}
   
   void    SetOuterFreonWidth(Float_t width)           {fOuterFreonWidth=width;}
   Float_t GetOuterFreonWidth()                   const{return fOuterFreonWidth;}
   
   void    SetInnerFreonLength(Float_t length)         {fInnerFreonLength=length;}
   Float_t GetInnerFreonLength()                  const{return fInnerFreonLength;}
   
   void    SetInnerFreonWidth(Float_t width)           {fInnerFreonWidth=width;}
   Float_t GetInnerFreonWidth()                   const{return fInnerFreonWidth;}
   
   void    SetFreonThickness(Float_t thickness)        {fFreonThickness=thickness;}
   Float_t GetFreonThickness()                    const{return fFreonThickness;}
   
   void    SetRadiatorToPads(Float_t distance)         {fRadiatorToPads=distance;}
   Float_t GetRadiatorToPads()                    const{return fRadiatorToPads;}   

   void    SetRotationAngle(Float_t rotAngle)          {fRotationAngle=rotAngle;}
   Float_t GetRotationAngle()                    const {return fRotationAngle;} 
          
private:
   Float_t fGapThickness;            // Gap Thickness
   Float_t fProximityGapThickness;   // Proximity Gap Thickness
   Float_t fQuartzLength;            // Quartz Length
   Float_t fQuartzWidth;             // Quartz Width
   Float_t fQuartzThickness;         // Quartz Thickness
   Float_t fOuterFreonLength;        // Outer Freon Length
   Float_t fOuterFreonWidth;         // Outer Freon Width
   Float_t fInnerFreonLength;        // Inner Freon Length
   Float_t fInnerFreonWidth;         // Inner Freon Width
   Float_t fFreonThickness;          // Freon Thickness
   Float_t fRadiatorToPads;          // Distance from radiator to pads
   Float_t fRotationAngle;           // Azimuthal rotation angle in X-Y plane  
   ClassDef(AliRICHGeometry,2)       // Chamber's main geometry parameters
};

inline AliRICHGeometry::AliRICHGeometry()
      	               :TObject()
{
// Define the default values:
      
   fGapThickness           =     8;   // Gap Thickness
   fProximityGapThickness  =   0.4;   // Proximity Gap Thickness
   fQuartzLength           =   133;   // Quartz Length
   fQuartzWidth            = 127.9;   // Quartz Width
   fQuartzThickness        =   0.5;   // Quartz Thickness
   fOuterFreonLength       =   133;   // Outer Freon Length
   fOuterFreonWidth        =  41.3;   // Outer Freon Width
   fInnerFreonLength       =   133;   // Inner Freon Length
   fInnerFreonWidth        =  41.3;   // Inner Freon Width
   fFreonThickness         =   1.5;   // Freon Thickness
   fRadiatorToPads         =     0;   // Distance from radiator to pads
   fRotationAngle          =   -30;   // Azimuthal rotation angle in X-Y plane  
}//AliRICHGeometry::ctor()

inline void AliRICHGeometry::Print()
{
   TObject::Print();
   cout<<"Radiator Gap thickness:  "<<fGapThickness<<endl;
   cout<<"Proximity Gap thickness: "<<fProximityGapThickness<<endl;
   cout<<"Quartz window length:    "<<fQuartzLength<<endl;
   cout<<"Quartz window width:     "<<fQuartzWidth<<endl;
   cout<<"Quartz window thickness: "<<fQuartzThickness<<endl;
   cout<<"Outer freon length:      "<<fOuterFreonLength<<endl;
   cout<<"Outer freon width:       "<<fOuterFreonWidth<<endl;
   cout<<"Inner freon length:      "<<fInnerFreonLength<<endl;
   cout<<"Inner freon width:       "<<fInnerFreonWidth<<endl;
   cout<<"Freon thickness:         "<<fFreonThickness<<endl;
   cout<<"RadiatorToPads:          "<<fRadiatorToPads<<endl;
   cout<<"Rotation angle:          "<<fRotationAngle<<endl;
}//void AliRICHGeometry::Print()

#endif //AliRICHGeometry_h
