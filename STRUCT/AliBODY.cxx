///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Alice external volume                                                    //
//  This class contains the description of the Alice external volume         //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliBODYClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:andreas.morsch@cern.ch">Andreas Morsch</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliRun.h"
#include "AliBODY.h"

ClassImp(AliBODY)
 
//_____________________________________________________________________________
AliBODY::AliBODY()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliBODY::AliBODY(const char *name, const char *title)
       : AliModule(name,title)
{
  //
  // Standard constructor of the Alice external volume
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliBODY::CreateGeometry()
{
  //
  // Create the geometry of the Alice external body
  //
  //Begin_Html
  /*
    <img src="picts/AliBODYTree.gif">
  */
  //End_Html
  //
  // If the ZDC is present we have an asymmetric box
  // made by a four sides polygone
  //  
  //Begin_Html
  /*
    <img src="picts/AliBODYLarge.gif">
  */
  //End_Html
  //
  // If the ZDC is not present make just a BOX
  //
  //Begin_Html
  /*
    <img src="picts/AliBODYSmall.gif">
  */
  //End_Html

  Float_t DALIC[10];
  Int_t *idtmed = gAlice->Idtmed();
  AliMC *pMC = AliMC::GetMC();
  //
  if(gAlice->GetModule("ZDC")) {
    //
    // If the ZDC is present we have an asymmetric box
    // made by a four sides polygone
    //
    DALIC[0]=45;
    DALIC[1]=360;
    DALIC[2]=4;
    DALIC[3]=2;
    DALIC[4]=-3000;
    DALIC[5]=0;
    DALIC[6]=2000;
    DALIC[7]=15000;
    DALIC[8]=0;
    DALIC[9]=2000;
    pMC->Gsvolu("ALIC","PGON",idtmed[1],DALIC,10);
  } else {
    //
    // If the ZDC is not present make just a BOX
    //
    DALIC[0]=2000;
    DALIC[1]=2000;
    DALIC[2]=3000;
    pMC->Gsvolu("ALIC","BOX ",idtmed[1],DALIC,3);
  }
}
 
//_____________________________________________________________________________
void AliBODY::CreateMaterials()
{
  Int_t ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  //
  AliMaterial(1,"Vacuum  $",1.e-16,1.e-16,1.e-16,1.e16,1.e16);
  AliMaterial(2,"Air     $",14.61,7.3,0.001205,30420,67500);
  AliMaterial(3,"Be      $", 9.01,4 ,1.848   ,35.30,36.70);
  //
  AliMedium(1,"Vacuum  $",1,0,ISXFLD,SXMGMX,10,1,0.1,0.1,10);
  AliMedium(2,"Air     $",2,0,ISXFLD,SXMGMX,10,-1,-0.1,0.1 ,-10);
  AliMedium(3,"Be pipe $",3,0,ISXFLD,SXMGMX,10,0.1,0.1,0.01,0.01);
}
 
//_____________________________________________________________________________
void AliBODY::DrawModule()
{
  //
  // Draw a view of the Alice outside box
  //
  AliMC* pMC = AliMC::GetMC();
  
  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother visible
  pMC->Gsatt("ALIC","SEEN",1);
  //
  // Set the volumes visible
  //
  pMC->Gdopt("hide","off");
  if(gAlice->GetModule("ZDC")) {
    //
    // ZDC is present
    //
    pMC->DefaultRange();
    pMC->Gdraw("alic", 40, 30, 0, 15, 10, .0014, .0014);
    pMC->Gdhead(1111, "Aice Main body with Zero Degree Calorimeter");
  } else {
    //
    // ZDC is not present
    //
    pMC->Gdraw("alic", 40, 30, 0, 10, 9, .0027, .0027);
    pMC->Gdhead(1111, "Aice Main body");
  }
  pMC->Gdman(18, 4, "MAN");
}
 

