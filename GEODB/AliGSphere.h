#ifndef ALIGSPHERE_H
#define ALIGSPHERE_H

#include "AliGShape.h"

class AliGSphere: public AliGShape {

    private:
        Float_t   fAspectRatio; // Relation between asumth and grid size (by default 1.0)
        Double_t* fCoTab;       // Table of cos(fPhimin) .... cos(Phi)
        Double_t* fCoThetaTab;  // Table of sin(gThemin) .... cos(Theta)
        Int_t     fNdiv;        // number of divisions
        Int_t     fNz;          // number of sections
        Double_t* fSiTab;       // Table of sin(fPhimin) .... sin(Phi)

    protected:
        Float_t faX;      // Coeff along Ox
        Float_t faY;      // Coeff along Oy
        Float_t faZ;      // Coeff along Oz
        Float_t fPhimax;  // maximum phi
        Float_t fPhimin;  // minimum phi
        Float_t fRmax;    // maximum radius
        Float_t fRmin;    // minimum radius
        Float_t fThemax;  // maximum theta
        Float_t fThemin;  // minimum theta

        virtual void  MakeTableOfCoSin();  // Create the table of the fSiTab; fCoTab
        virtual void  PaintGLPoints(Float_t *vertex);

    public:
        AliGSphere();
        AliGSphere(Text_t *name, Text_t *title, Float_t rmin, Float_t rmax, Float_t themin, Float_t themax, Float_t phimin, Float_t phimax);
        AliGSphere(Text_t *name, Text_t *title, Float_t rmax);
        AliGSphere(AliGSphere *sphere);
	virtual ~AliGSphere(); // Destructor

        virtual void  Draw(Option_t *option);
        virtual void  DrawShape(Option_t *option); // *MENU*
        virtual Int_t GetNumberOfDivisions () const {if (fNdiv) return fNdiv; else return 0; /*kDiv;*/}
        virtual void  Paint(Option_t *option);
        virtual void  SetEllipse(Float_t *factors);
        virtual void  SetNumberOfDivisions (Int_t p);
        virtual void  SetPoints(Float_t *buff);
                void  Sizeof3D() const;

   ClassDef(AliGSphere,1) //Simple sphere class
};
#endif
