// $Id$
//
// Class AliMUONSt2GeometryBuilder
// -------------------------------
// MUON Station2 coarse geometry construction class.
// Extracted from AliMUONv1
// by Ivana Hrivnacova, IPN Orsay
// Included in AliRoot 2004/01/23

#include <TVirtualMC.h>
#include <TGeoMatrix.h>

#include "AliMUONSt2GeometryBuilder.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONChamberGeometry.h"

ClassImp(AliMUONSt2GeometryBuilder)

//______________________________________________________________________________
AliMUONSt2GeometryBuilder::AliMUONSt2GeometryBuilder(AliMUON* muon)
 : AliMUONVGeometryBuilder(&muon->Chamber(2), &muon->Chamber(3)),
   fMUON(muon)
{
// Standard constructor

}

//______________________________________________________________________________
AliMUONSt2GeometryBuilder::AliMUONSt2GeometryBuilder()
 : AliMUONVGeometryBuilder(),
   fMUON(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONSt2GeometryBuilder::AliMUONSt2GeometryBuilder(const AliMUONSt2GeometryBuilder& rhs)
  : AliMUONVGeometryBuilder(rhs)
{
  Fatal("Copy constructor", 
        "Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONSt2GeometryBuilder::~AliMUONSt2GeometryBuilder() {
//
}

//______________________________________________________________________________
AliMUONSt2GeometryBuilder& 
AliMUONSt2GeometryBuilder::operator = (const AliMUONSt2GeometryBuilder& rhs) 
{
  // check assignement to self
  if (this == &rhs) return *this;

  Fatal("operator=", 
        "Assignment operator is not implemented.");
    
  return *this;  
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONSt2GeometryBuilder::CreateGeometry() 
{
// From AliMUONv1::CreateGeometry()

//
//********************************************************************
//                            Station 2                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers

     // Get tracking medias Ids     
     Int_t *idtmed = fMUON->GetIdtmed()->GetArray()-1099;
     Int_t idAir= idtmed[1100]; // medium 1
     Int_t idAlu1=idtmed[1103]; // medium 4
     Int_t idAlu2=idtmed[1104]; // medium 5
     Int_t idGas=idtmed[1108];  // medium 9 = Ar-CO2 gas (80%+20%)
     Bool_t frameCrosses=kTRUE;     

     // Rotation matrices in the x-y plane  
     // phi= 0 deg
     Int_t irot1;
     fMUON->AliMatrix(irot1,  90.,   0., 90.,  90., 0., 0.);
     // phi= 90 deg
     Int_t irot2;
     fMUON->AliMatrix(irot2,  90.,  90., 90., 180., 0., 0.);

     AliMUONChamber* iChamber = GetChamber(2);
     AliMUONChamber* iChamber1 = iChamber;
     AliMUONChamber* iChamber2 = GetChamber(3);
     
     // Half of the total thickness of frame crosses (including DAlu)
     // for each chamber in stations 1 and 2:
     // 3% of X0 of composite material,
     // but taken as Aluminium here, with same thickness in number of X0
     Float_t dframez = 3. * 8.9 / 100;
     // DGas and DAlu not changed from standard values
     Double_t zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
             // The same parameters are defined in builder for station 1 
     
     // Mother volume
     // Outer excess and inner recess for mother volume radius
     // with respect to ROuter and RInner
     Float_t dframep=.001; // Value for station 3 should be 6 ...
     // Width (RdPhi) of the frame crosses for stations 1 and 2 (cm)
     // Float_t dframep1=.001;
     Float_t phi=2*TMath::Pi()/12/2;
             // The same parameters are defined in builder for station 1 

     Float_t tpar[3];
     Double_t dstation = (-iChamber2->Z()) - (-iChamber1->Z());
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5;

     gMC->Gsvolu("S03M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("S04M", "TUBE", idAir, tpar, 3);

     // CHANGED
     //gMC->Gspos("S03M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     //gMC->Gspos("S04M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
     GetChamber(2)->GetGeometry()->AddEnvelope("S03M", false);
     GetChamber(3)->GetGeometry()->AddEnvelope("S04M", false);
     
     gMC->Gsbool("S03M", "L3DO");
     gMC->Gsbool("S03M", "L3O1");
     gMC->Gsbool("S03M", "L3O2");
     gMC->Gsbool("S04M", "L3DO");
     gMC->Gsbool("S04M", "L3O1");
     gMC->Gsbool("S04M", "L3O2");

// // Aluminium frames
// // Outer frames
//      pgpar[0] = 360/12/2;
//      pgpar[1] = 360.;
//      pgpar[2] = 12.;
//      pgpar[3] =   2;
//      pgpar[4] = -dframez/2;
//      pgpar[5] = iChamber->ROuter();
//      pgpar[6] = pgpar[5]+dframep;
//      pgpar[7] = +dframez/2;
//      pgpar[8] = pgpar[5];
//      pgpar[9] = pgpar[6];
//      gMC->Gsvolu("S03O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gsvolu("S04O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gspos("S03O",1,"S03M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S03O",2,"S03M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("S04O",1,"S04M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S04O",2,"S04M", 0.,0.,+zfpos,  0,"ONLY");
// //
// // Inner frame
//      tpar[0]= iChamber->RInner()-dframep;
//      tpar[1]= iChamber->RInner();
//      tpar[2]= dframez/2;
//      gMC->Gsvolu("S03I", "TUBE", idAlu1, tpar, 3);
//      gMC->Gsvolu("S04I", "TUBE", idAlu1, tpar, 3);

//      gMC->Gspos("S03I",1,"S03M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S03I",2,"S03M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("S04I",1,"S04M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S04I",2,"S04M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     if (frameCrosses) {
         // outside gas
         // security for inside mother volume
 
         // ADDED !! Repeated     
         Float_t dframep1 = 11.0;
         Float_t dframez = 3. * 8.9 / 100;

         Float_t bpar[3];
	 bpar[0] = (iChamber->ROuter() - iChamber->RInner())
	   * TMath::Cos(TMath::ASin(dframep1 /
				   (iChamber->ROuter() - iChamber->RInner())))
	   / 2.0;
	 bpar[1] = dframep1/2;
	 // total thickness will be (4 * bpar[2]) for each chamber,
	 // which has to be equal to (2 * dframez) - DAlu
	 bpar[2] = (2.0 * dframez - iChamber->DAlu()) / 4.0;
	 gMC->Gsvolu("S03B", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("S04B", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("S03B",1,"S03M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S03B",2,"S03M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S03B",3,"S03M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S03B",4,"S03M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S03B",5,"S03M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S03B",6,"S03M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S03B",7,"S03M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S03B",8,"S03M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    irot2,"ONLY");
	 
	 gMC->Gspos("S04B",1,"S04M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S04B",2,"S04M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S04B",3,"S04M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S04B",4,"S04M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S04B",5,"S04M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S04B",6,"S04M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S04B",7,"S04M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S04B",8,"S04M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    irot2,"ONLY");
     }
//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->RInner();
     tpar[1]= iChamber->ROuter();
     tpar[2] = (iChamber->DGas()+iChamber->DAlu())/2;
     gMC->Gsvolu("S03A", "TUBE", idAlu2, tpar, 3);
     gMC->Gsvolu("S04A", "TUBE", idAlu2, tpar, 3);
     gMC->Gspos("S03A", 1, "S03M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("S04A", 1, "S04M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->DGas();
     tpar[2] = iChamber->DGas()/2;
     gMC->Gsvolu("S03G", "TUBE", idGas, tpar, 3);
     gMC->Gsvolu("S04G", "TUBE", idGas, tpar, 3);
     gMC->Gspos("S03G", 1, "S03A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("S04G", 1, "S04A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas 
     // NONE: chambers are sensitive everywhere
//      if (frameCrosses) {

// 	 dr = (iChamber->ROuter() - iChamber->RInner());
// 	 bpar[0] = TMath::Sqrt(dr*dr-dframep1*dframep1/4)/2;
// 	 bpar[1] = dframep1/2;
// 	 bpar[2] = iChamber->DGas()/2;
// 	 gMC->Gsvolu("S03F", "BOX", idAlu1, bpar, 3);
// 	 gMC->Gsvolu("S04F", "BOX", idAlu1, bpar, 3);
	 
// 	 gMC->Gspos("S03F",1,"S03G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S03F",2,"S03G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S03F",3,"S03G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    irot2,"ONLY");
// 	 gMC->Gspos("S03F",4,"S03G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    irot2,"ONLY");
	 
// 	 gMC->Gspos("S04F",1,"S04G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S04F",2,"S04G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S04F",3,"S04G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    irot2,"ONLY");
// 	 gMC->Gspos("S04F",4,"S04G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    irot2,"ONLY");
//      }
/*
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[2];
     iChamber2 =(AliMUONChamber*) (*fChambers)[3];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = zpos2 - zpos1;
     // DGas and DAlu not changed from standard values
     zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5;

     gMC->Gsvolu("S03M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("S04M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("S03M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("S04M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
     gMC->Gsbool("S03M", "L3DO");
     gMC->Gsbool("S03M", "L3O1");
     gMC->Gsbool("S03M", "L3O2");
     gMC->Gsbool("S04M", "L3DO");
     gMC->Gsbool("S04M", "L3O1");
     gMC->Gsbool("S04M", "L3O2");

// // Aluminium frames
// // Outer frames
//      pgpar[0] = 360/12/2;
//      pgpar[1] = 360.;
//      pgpar[2] = 12.;
//      pgpar[3] =   2;
//      pgpar[4] = -dframez/2;
//      pgpar[5] = iChamber->ROuter();
//      pgpar[6] = pgpar[5]+dframep;
//      pgpar[7] = +dframez/2;
//      pgpar[8] = pgpar[5];
//      pgpar[9] = pgpar[6];
//      gMC->Gsvolu("S03O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gsvolu("S04O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gspos("S03O",1,"S03M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S03O",2,"S03M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("S04O",1,"S04M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S04O",2,"S04M", 0.,0.,+zfpos,  0,"ONLY");
// //
// // Inner frame
//      tpar[0]= iChamber->RInner()-dframep;
//      tpar[1]= iChamber->RInner();
//      tpar[2]= dframez/2;
//      gMC->Gsvolu("S03I", "TUBE", idAlu1, tpar, 3);
//      gMC->Gsvolu("S04I", "TUBE", idAlu1, tpar, 3);

//      gMC->Gspos("S03I",1,"S03M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S03I",2,"S03M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("S04I",1,"S04M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S04I",2,"S04M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     if (frameCrosses) {
         // outside gas
         // security for inside mother volume
	 bpar[0] = (iChamber->ROuter() - iChamber->RInner())
	   * TMath::Cos(TMath::ASin(dframep1 /
				   (iChamber->ROuter() - iChamber->RInner())))
	   / 2.0;
	 bpar[1] = dframep1/2;
	 // total thickness will be (4 * bpar[2]) for each chamber,
	 // which has to be equal to (2 * dframez) - DAlu
	 bpar[2] = (2.0 * dframez - iChamber->DAlu()) / 4.0;
	 gMC->Gsvolu("S03B", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("S04B", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("S03B",1,"S03M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S03B",2,"S03M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S03B",3,"S03M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S03B",4,"S03M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S03B",5,"S03M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S03B",6,"S03M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S03B",7,"S03M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S03B",8,"S03M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    irot2,"ONLY");
	 
	 gMC->Gspos("S04B",1,"S04M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S04B",2,"S04M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S04B",3,"S04M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S04B",4,"S04M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S04B",5,"S04M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S04B",6,"S04M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S04B",7,"S04M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S04B",8,"S04M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    irot2,"ONLY");
     }
//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->RInner();
     tpar[1]= iChamber->ROuter();
     tpar[2] = (iChamber->DGas()+iChamber->DAlu())/2;
     gMC->Gsvolu("S03A", "TUBE", idAlu2, tpar, 3);
     gMC->Gsvolu("S04A", "TUBE", idAlu2, tpar, 3);
     gMC->Gspos("S03A", 1, "S03M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("S04A", 1, "S04M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->DGas();
     tpar[2] = iChamber->DGas()/2;
     gMC->Gsvolu("S03G", "TUBE", idGas, tpar, 3);
     gMC->Gsvolu("S04G", "TUBE", idGas, tpar, 3);
     gMC->Gspos("S03G", 1, "S03A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("S04G", 1, "S04A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas 
     // NONE: chambers are sensitive everywhere
//      if (frameCrosses) {

// 	 dr = (iChamber->ROuter() - iChamber->RInner());
// 	 bpar[0] = TMath::Sqrt(dr*dr-dframep1*dframep1/4)/2;
// 	 bpar[1] = dframep1/2;
// 	 bpar[2] = iChamber->DGas()/2;
// 	 gMC->Gsvolu("S03F", "BOX", idAlu1, bpar, 3);
// 	 gMC->Gsvolu("S04F", "BOX", idAlu1, bpar, 3);
	 
// 	 gMC->Gspos("S03F",1,"S03G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S03F",2,"S03G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S03F",3,"S03G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    irot2,"ONLY");
// 	 gMC->Gspos("S03F",4,"S03G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    irot2,"ONLY");
	 
// 	 gMC->Gspos("S04F",1,"S04G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S04F",2,"S04G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S04F",3,"S04G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    irot2,"ONLY");
// 	 gMC->Gspos("S04F",4,"S04G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    irot2,"ONLY");
//      }
     }
*/
}

//______________________________________________________________________________
void AliMUONSt2GeometryBuilder::SetTransformations()
{
// Defines the transformations for the station2 chambers.
// ---

  AliMUONChamber* iChamber1 = GetChamber(2);
  Double_t zpos1 = - iChamber1->Z(); 
  iChamber1->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos1));

  AliMUONChamber* iChamber2 = GetChamber(3);
  Double_t zpos2 = - iChamber2->Z(); 
  iChamber2->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos2));
}

//______________________________________________________________________________
void AliMUONSt2GeometryBuilder::SetSensitiveVolumes()
{
// Defines the sensitive volumes for station2 chambers.
// ---

  GetChamber(2)->GetGeometry()->SetSensitiveVolume("S03G");
  GetChamber(3)->GetGeometry()->SetSensitiveVolume("S04G");
}
