#ifndef ALIRICHGEOMETRY_H
#define ALIRICHGEOMETRY_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>


class AliRICHGeometry : //Chamber geometry
public TObject {
 public:
    // Radiator Thickness
    void   SetGapThickness(Float_t thickness) {fGapThickness=thickness;} 
    // Proximity Gap Thickness
    void   SetProximityGapThickness(Float_t thickness) {fProximityGapThickness=thickness;}
    // Quartz Length
    void   SetQuartzLength(Float_t length) {fQuartzLength=length;}
    // Quartz Width
    void   SetQuartzWidth(Float_t width) {fQuartzWidth=width;}
    // Quartz Thickness
    void   SetQuartzThickness(Float_t thickness) {fQuartzThickness=thickness;}
    // Freon Length
    void   SetOuterFreonLength(Float_t length) {fOuterFreonLength=length;}
    // Freon Width
    void   SetOuterFreonWidth(Float_t width) {fOuterFreonWidth=width;}
    // Freon Length
    void   SetInnerFreonLength(Float_t length) {fInnerFreonLength=length;}
    // Freon Width
    void   SetInnerFreonWidth(Float_t width) {fInnerFreonWidth=width;}
    // Freon Thickness
    void   SetFreonThickness(Float_t thickness) {fFreonThickness=thickness;}
    // Freon Thickness
    void   SetRadiatorToPads(Float_t distance) {fRadiatorToPads=distance;}

    // Radiator thickness
    Float_t  GetGapThickness() {return fGapThickness;}
    // Proximity Gap thickness
    Float_t  GetProximityGapThickness() {return fProximityGapThickness;}
    // Quartz Length
    Float_t  GetQuartzLength() {return fQuartzLength;}
    // Quartz Width
    Float_t  GetQuartzWidth() {return fQuartzWidth;}
    // Quartz Thickness
    Float_t  GetQuartzThickness() {return fQuartzThickness;}
    // Freon Length
    Float_t  GetOuterFreonLength() {return fOuterFreonLength;}
    // Freon Width
    Float_t  GetOuterFreonWidth() {return fOuterFreonWidth;}
    // Freon Length
    Float_t  GetInnerFreonLength() {return fInnerFreonLength;}
    // Freon Width
    Float_t  GetInnerFreonWidth() {return fInnerFreonWidth;}
    // Freon Thickness
    Float_t  GetFreonThickness() {return fFreonThickness;}
    // Get distance between radiator and pads
    Float_t  GetRadiatorToPads() {return fRadiatorToPads;}   
    
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

    ClassDef(AliRICHGeometry,1)
};
#endif




