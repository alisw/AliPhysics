/*
**********No misalignment*******
i=0 
matrix global_1 - tr=1  rot=1  refl=0  scl=0
  0.766044    0.000000    0.642788    Tx = 310.215729
  0.642788    0.000000   -0.766044    Ty = -369.700684
  0.000000    1.000000    0.000000    Tz =   0.000000
i=1 
matrix global_1 - tr=1  rot=1  refl=0  scl=0
  0.939693    0.000000    0.342020    Tx = 165.062332
  0.342020    0.000000   -0.939693    Ty = -453.505035
  0.000000    1.000000    0.000000    Tz =   0.000000
i=2 
matrix global_1 - tr=1  rot=1  refl=0  scl=0
  1.000000    0.000000    0.000000    Tx =   0.000000
  0.000000    0.000000   -1.000000    Ty = -482.609985
  0.000000    1.000000    0.000000    Tz =   0.000000
i=3 
matrix global_1 - tr=1  rot=1  refl=0  scl=0
  0.939693    0.000000   -0.342020    Tx = -165.062332
 -0.342020    0.000000   -0.939693    Ty = -453.505035
  0.000000    1.000000    0.000000    Tz =   0.000000

*************
=========Should be========

mod 4 (i=0)
matrix  - tr=1  rot=1  refl=0  scl=0
  0.767337   -0.004312    0.641229    Tx = 313.381746
  0.641189   -0.007942   -0.767342    Ty = -372.024764
  0.008402    0.999959   -0.003330    Tz =  -2.354855

mod 3 (i=1)
matrix  - tr=1  rot=1  refl=0  scl=0
  0.940798    0.002136    0.338961    Tx = 165.078518
  0.338959    0.001227   -0.940800    Ty = -454.148853
 -0.002426    0.999997    0.000431    Tz =  -2.801675

mod 3 (i=2)
matrix  - tr=1  rot=1  refl=0  scl=0
  0.999976   -0.000756   -0.006876    Tx =  -0.083626
 -0.006877   -0.001172   -0.999976    Ty = -482.172358
  0.000748    0.999999   -0.001178    Tz =  -1.906440

mod 3 (i=3)
matrix  - tr=1  rot=1  refl=0  scl=0
  0.934365   -0.002092   -0.356311    Tx = -167.233337
 -0.356313   -0.000801   -0.934366    Ty = -460.901723
  0.001670    0.999997   -0.001494    Tz =  -1.821245
  
------------Correct angles, but zero offsets------
  0.767338   -0.004312    0.641229    Tx = 309.463421
  0.641188   -0.007942   -0.767343    Ty = -370.327157
  0.008401    0.999959   -0.003329    Tz =  -1.606798
i=1 
matrix global_1 - tr=1  rot=1  refl=0  scl=0
  0.940798    0.002136    0.338961    Tx = 163.586079
  0.338959    0.001227   -0.940800    Ty = -454.039582
 -0.002425    0.999997    0.000430    Tz =   0.207687
i=2 
matrix global_1 - tr=1  rot=1  refl=0  scl=0
  0.999976   -0.000756   -0.006876    Tx =  -3.318427
 -0.006877   -0.001172   -0.999976    Ty = -482.598242
  0.000748    0.999999   -0.001177    Tz =  -0.568114
i=3 
matrix global_1 - tr=1  rot=1  refl=0  scl=0
  0.934365   -0.002092   -0.356311    Tx = -171.959374
 -0.356314   -0.000801   -0.934366    Ty = -450.934413
  0.001669    0.999997   -0.001494    Tz =  -0.720939

  
  
  
*/

void MakeFinalAlignment(){
  // Create ideal (no misalignment) object for PHOS

  
  const AliPHOSGeometry *phosGeom = AliPHOSGeometry::GetInstance("Run2", "Run2");
  if (!phosGeom) {
    Error(macroName, "Cannot obtain AliPHOSGeometry singleton.\n");
    return;
  }

  //Activate CDB storage and load geometry from CDB
  //[Part of code, taken from ITS version of MakexxxFullMisalignment
  AliCDBManager * cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://./OCDB");
  cdb->SetRun(0);

  AliPHOSEMCAGeometry * emca = phosGeom->GetEMCAGeometry();
  TClonesArray alobj("AliAlignObjParams", 16);
   
  const Double_t dpsi = 0., dtheta = 0., dphi = 0.;
  Int_t iIndex = 0; //let all modules have index=0 in a layer with no LUT
  const AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer, iIndex);
  Int_t i = 0;

  // Alignment for 5 PHOS modules
 Double_t ideal1[9] = {0.766044,    0.000000,   0.642788,  
		       0.642788,    0.000000,  -0.766044,    
		       0.000000,    1.000000,   0.000000} ;    
 Double_t final1[9] = {0.767337,   -0.004312,   0.641229,   
		       0.641189,   -0.007942,  -0.767342, 
		       0.008402,    0.999959,  -0.003330} ;
  TGeoRotation rot1;
  rot1.SetMatrix(ideal1);
  TGeoRotation inv =rot1.Inverse() ;
  TGeoRotation rotF1;
  rotF1.SetMatrix(final1);
  rotF1.MultiplyBy(&inv) ;
  Double_t dX= 313.381746-309.463421; 
  Double_t dY=-372.024764+370.327157; 
  Double_t dZ=-2.354855 +1.606798 ;

 AliAlignObjParams * mod1 = new(alobj[i++]) AliAlignObjParams("PHOS/Module1",volid, dX, dY, dZ, 0., 0., 0., kTRUE);
  mod1->SetRotation(rotF1);

 
  Double_t ideal2[9] = {0.939693,    0.000000,    0.342020,
		        0.342020,    0.000000,   -0.939693,
		        0.000000,    1.000000,   0.000000} ;    
  Double_t final2[9] = {0.940798,    0.002136,    0.338961,
		        0.338959,    0.001227,   -0.940800, 
		       -0.002426,    0.999997,    0.000431} ;
  TGeoRotation rot2;
  rot2.SetMatrix(ideal2);
  TGeoRotation inv2 =rot2.Inverse() ;
  TGeoRotation rotF2;
  rotF2.SetMatrix(final2);
  rotF2.MultiplyBy(&inv2) ;
  dX= 165.078518-163.586079; 
  dY=-454.148853+454.039582; 
  dZ=-2.801675-0.207687 ;

  AliAlignObjParams * mod2 = new(alobj[i++]) AliAlignObjParams("PHOS/Module2",volid, dX, dY, dZ, 0., 0., 0., kTRUE);
  mod2->SetRotation(rotF2);
  
  Double_t ideal3[9] = {1.000000,    0.000000,    0.000000,
		        0.000000,    0.000000,   -1.000000,
		        0.000000,    1.000000,    0.000000} ;    
  Double_t final3[9] = {0.999976,   -0.000756,   -0.006876,
		       -0.006877,   -0.001172,   -0.999976,
		        0.000748,    0.999999,   -0.001178} ;
  TGeoRotation rot3;
  rot3.SetMatrix(ideal3);
  TGeoRotation inv3 =rot3.Inverse() ;
  TGeoRotation rotF3;
  rotF3.SetMatrix(final3);
  rotF3.MultiplyBy(&inv3) ;
  dX= -0.083626+3.318427; 
  dY=-482.172358+482.598242; 
  dZ=-1.906440+0.568114 ;

  AliAlignObjParams * mod3 = new(alobj[i++]) AliAlignObjParams("PHOS/Module3",volid, dX, dY, dZ, 0., 0., 0., kTRUE);
  mod3->SetRotation(rotF3);
  
 Double_t ideal4[9] = {0.939693,    0.000000,  -0.342020,
		      -0.342020,    0.000000,   -0.939693,
		        0.000000,    1.000000,    0.000000} ;    
  Double_t final4[9] = {0.934365,   -0.002092,   -0.356311,
		       -0.356313,   -0.000801,   -0.934366,
		        0.001670,    0.999997,   -0.001494} ;
  TGeoRotation rot4;
  rot4.SetMatrix(ideal4);
  TGeoRotation inv4 =rot4.Inverse() ;
  TGeoRotation rotF4;
  rotF4.SetMatrix(final4);
  rotF4.MultiplyBy(&inv4) ;
  dX=-167.233337+171.959374; 
  dY=-460.901723+450.934413; 
  dZ=-1.821245+0.720939 ;
  
  AliAlignObjParams * mod4 = new(alobj[i++]) AliAlignObjParams("PHOS/Module4",volid, dX, dY, dZ, 0., 0., 0., kTRUE);
  mod4->SetRotation(rotF4);

  new(alobj[i++]) AliAlignObjParams("PHOS/Module5",
	  volid, 0., 0., 0., 0., 0., 0., kTRUE);

  const Double_t dx = 0., dy = 0., dz = 0. ;
  // Alignment of CPV modules
  new(alobj[i++]) AliAlignObjParams("PHOS/Module1/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module2/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module3/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module4/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module5/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
 
  Double_t displacement=0 ;
  // Alignment for PHOS cradle
  new(alobj[i++]) AliAlignObjParams("PHOS/Cradle0",
	  volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Cradle1",
	  volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);

  // Alignment for cradle wheels
  new(alobj[i++]) AliAlignObjParams("PHOS/Wheel0",
	  volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Wheel1",
	  volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Wheel2",
	  volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Wheel3",
	  volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);

  AliCDBMetaData md;
  md.SetResponsible("Dmitri Peressounko");
  md.SetComment("Ideal alignment objects for PHOS");
  md.SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("PHOS/Align/Data",0,AliCDBRunRange::Infinity());
  cdb->Put(&alobj, id, &md);

  alobj.Delete();
}

 