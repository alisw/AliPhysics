#ifndef AliRICHGeometry_h
#define AliRICHGeometry_h


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <Riostream.h>
#include <TObject.h>

class AliRICHGeometry :public TObject  
{
public:
// ctor & dtor staff      
   inline AliRICHGeometry(); // default ctor
// inline methods:          
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
       
   void    SetAlphaAngle(Float_t alphaAngle)           {fAlphaAngle=alphaAngle;} // Angle between modules in YZ plane
   Float_t GetAlphaAngle()                       const {return fAlphaAngle;}     // Angle between modules in YZ plane
   
   void    SetBetaAngle(Float_t betaAngle)             {fBetaAngle=betaAngle;}   // Angle between modules in XY plane
   Float_t GetBetaAngle()                        const {return fBetaAngle;}      // Angle between modules in XY plane
   
   void    SetOffset(Float_t offset)                   {fOffset=offset;}         // Modules offset from IP
   Float_t GetOffset()                           const {return fOffset;}         // Modules offset from IP
   
   inline void Print(Option_t *option)const;
          
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
   Float_t fAlphaAngle;              // Angle between modules in YZ plane
   Float_t fBetaAngle;               // Angle between modules in XY plane
   Float_t fOffset;                  // Modules offset from IP 
   ClassDef(AliRICHGeometry,3)       // Chamber's main geometry parameters
};

inline AliRICHGeometry::AliRICHGeometry()
      	               :TObject()
{
// Define the default values:
      
   SetGapThickness         (8);     // Gap Thickness
   SetProximityGapThickness(0.4);   // Proximity Gap Thickness
   SetQuartzLength         (133);   // Quartz Length
   SetQuartzWidth          (127.9); // Quartz Width
   SetQuartzThickness      (0.5);   // Quartz Thickness
   SetOuterFreonLength     (133);   // Outer Freon Length
   SetOuterFreonWidth      (41.3);  // Outer Freon Width
   SetInnerFreonLength     (133);   // Inner Freon Length
   SetInnerFreonWidth      (41.3);  // Inner Freon Width
   SetFreonThickness       (1.5);   // Freon Thickness
   SetRadiatorToPads       (0);     // Distance from radiator to pads
   SetRotationAngle        (-60);   // Azimuthal rotation angle in X-Y plane
   SetAlphaAngle           (19.5);  // Angle between modules in YZ plane
   SetBetaAngle            (20);    // Angle vetween modules in XY plane     
   SetOffset               (490+1.267); // ???1.267??? Modules offset from IP 
}//AliRICHGeometry::ctor()

inline void AliRICHGeometry::Print(Option_t *option)const
{
   TObject::Print(option);
   cout<<"Radiator Gap thickness:  "<<GetGapThickness()          <<endl;
   cout<<"Proximity Gap thickness: "<<GetProximityGapThickness() <<endl;
   cout<<"Quartz window length:    "<<GetQuartzLength()          <<endl;
   cout<<"Quartz window width:     "<<GetQuartzWidth()           <<endl;
   cout<<"Quartz window thickness: "<<GetQuartzThickness()       <<endl;
   cout<<"Outer freon length:      "<<GetOuterFreonLength()      <<endl;
   cout<<"Outer freon width:       "<<GetOuterFreonWidth()       <<endl;
   cout<<"Inner freon length:      "<<GetInnerFreonLength()      <<endl;
   cout<<"Inner freon width:       "<<GetInnerFreonWidth()       <<endl;
   cout<<"Freon thickness:         "<<GetFreonThickness()        <<endl;
   cout<<"RadiatorToPads:          "<<GetRadiatorToPads()        <<endl;
   cout<<"Rotation angle:          "<<GetRotationAngle()         <<endl;
   cout<<"Alpha Angle:             "<<GetAlphaAngle()            <<endl;
   cout<<"Beta angle:              "<<GetBetaAngle()             <<endl;
   cout<<"Modules offset from IP:  "<<GetOffset()                <<endl;
}//void AliRICHGeometry::Print()

//______________________________________________________
// Definition and manipulation with parameters describing RICH parametrised geometry.

#endif //AliRICHGeometry_h
