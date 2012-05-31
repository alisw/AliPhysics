
void MakePipeMisalignment(Int_t inputRun=0,
			  const char* ocdbPath="local://$ALICE_ROOT/OCDB",
			  Bool_t verbose=true, 
			  const char* outFileName=0)
{
  // Creates misalignment of the beam pipe, frame and FMD2&3.
  // 
  // The misalignment is based on the ITS global alignment object with
  // the exception of the beam pipe where only the translations are
  // used ignoring the rotation matrix.
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB");
  cdb->SetRun(inputRun);


  // Get the geometry 
  AliGeomManager::LoadGeometry();

  // Get ITS alignment 
  AliCDBEntry* itsEntry = cdb->Get("ITS/Align/Data");
  if (!itsEntry) { 
    Error("MakePipeMisalignment", "Couldn't get ITS alignment data");
    return;
  } 
  TClonesArray* itsArray = static_cast<TClonesArray*>(itsEntry->GetObject());
  if (!itsArray) { 
    Error("MakePipeMisalignment", "No array in CDB entry");
    return;
  }
  AliAlignObjParams* itsAlign = 
    static_cast<AliAlignObjParams*>(itsArray->At(0));
  if (!itsAlign) { 
    Error("MakePipeMisalignment", "No alignment object in CDB entry");
    return;
  }
  
  if (verbose) {
    Info("MakePipeMisalignment", "ITS alignment:");
    itsAlign->Print();
  }

  Double_t itsTrans[3];
  Double_t itsRot[3];
  itsAlign->GetTranslation(itsTrans);
  itsAlign->GetAngles(itsRot);

  Double_t local[] = { 0, 0, 0 };

  if (!gGeoManager->cd("/ALIC_1/CP_1/Cp1_1/CpSupC_1")) return;
  Double_t pCpSupCA[3];
  gGeoManager->LocalToMaster(local, pCpSupCA);

  if (!gGeoManager->cd("/ALIC_1/CP_1/Cp3_1/CpSupC_3")) return;
  Double_t pCpSupCC[3];
  gGeoManager->LocalToMaster(local, pCpSupCC);

  if (!gGeoManager->cd("/ALIC_1/CP_1/Cp1_1")) return;
  Double_t pCp1[3];
  gGeoManager->LocalToMaster(local, pCp1);

  if (!gGeoManager->cd("/ALIC_1/CP_1/Cp1_1/CP1MO_1/CP1AT_1")) return;
  do { 
    TGeoNode* node = gGeoManager->GetCurrentNode();
    if (!node) break;
    TGeoVolume* volume = node->GetVolume();
    if (!volume) break;
    TGeoBBox* bbox = static_cast<TGeoBBox*>(volume->GetShape());
    if (!bbox) break;
    local[2] = -bbox->GetDZ();
  } while (false);
  if (local[2] == 0) { 
    Error("MakePipeMisalignement", "Failed to find end of CP1 section");
    return;
  }
  Double_t pCP1AT[3];
  gGeoManager->LocalToMaster(local, pCP1AT);


  
  // Now, we need to misalign the C-side collar according to the ITS
  // misalignment. 
  TGeoHMatrix* itsMat = new TGeoHMatrix();
  itsAlign->GetMatrix(*itsMat);
  // itsTrans[2] = 0;
  // itsMat->SetTranslation(itsTrans);
  Double_t tCpSupCA[3];
  itsMat->LocalToMaster(pCpSupCA, tCpSupCA);
  Double_t tCpSupCC[3];
  itsMat->LocalToMaster(pCpSupCC, tCpSupCC);

  Info("MakePipeMisalignment", "\n"
       "\tCollar A-side       (%f,%f,%f)\n" 
       "\tCollar C-side       (%f,%f,%f)\n"
       "\tBeam pipe center    (%f,%f,%f)\n"
       "\tEnd of CP1 A-side   (%f,%f,%f)\n"
       "\tCollar A-side after (%f,%f,%f)\n" 
       "\tCollar C-side after (%f,%f,%f)", 
       pCpSupCA[0], pCpSupCA[1], pCpSupCA[2],
       pCpSupCC[0], pCpSupCC[1], pCpSupCC[2],
       pCp1[0],     pCp1[1],     pCp1[2],
       pCP1AT[0],   pCP1AT[1],   pCP1AT[2],
       tCpSupCA[0], tCpSupCA[1], tCpSupCA[2],
       tCpSupCC[0], tCpSupCC[1], tCpSupCC[2]);
  
  // Distance from A-side end to mis-aligned C-side Collar
  Double_t dX = - tCpSupCC[0] + pCP1AT[0] ;
  Double_t dY = - tCpSupCC[1] + pCP1AT[1] ;
  Double_t dZ = - tCpSupCC[2] + pCP1AT[2] ;
  
  // Rotation angle around X and Y
  Double_t rotX = TMath::ATan2(dY, dZ);
  Double_t rotY = TMath::ATan2(dX, dZ);

  
  Double_t dZc = pCp1[2] - pCP1AT[2];
  Double_t dXc = TMath::Tan(rotY) * dZc;
  Double_t dYc = TMath::Tan(rotX) * dZc;


  Double_t dZt = pCpSupCC[2] - pCP1AT[2] + itsTrans[2];

  Info("MakePipeMisalignment", "\n"       
       "\tDistances (%f,%f,%f)\n" 
       "\trotX=%f rotY=%f\n" 
       "\tCP1     displacement (x,y)=(%f,%f) dZ=%f\n"
       "\tCpSubCA displacement (x,y)=(%f,%f) dZ=%f",
       dX, dY, dZ,
       rotX * TMath::RadToDeg(), rotY * TMath::RadToDeg(),
       dXc, dYc, dZc,
       TMath::Tan(rotY)*dZt, TMath::Tan(rotX)*dZt, dZt);


  TGeoRotation* rot = new TGeoRotation();
  rot->RotateX(rotX*TMath::RadToDeg());
  rot->RotateY(rotY*TMath::RadToDeg());
  TGeoCombiTrans* trans = new TGeoCombiTrans(dXc, dYc, 0, rot);
  trans->Print();

  // Output array 
  TClonesArray* structArray = new TClonesArray("AliAlignObjParams", 20);

  //dummy vol id
  UShort_t dvoluid = 
    AliGeomManager::LayerToVolUID(AliGeomManager::kInvalidLayer,0); 

  //base of symbolic name corresponding to base of path "ALIC_1/B077_1/BSEGMO";
  const char* baseSymName = "FRAME/Sector"; 
  for(Int_t sm=0; sm<18; sm++){
    TString symname =  baseSymName;
    symname         += sm;
    new((*structArray)[sm]) AliAlignObjParams(symname.Data(),dvoluid,
					     0.,0.,0.,0.,0.,0.,kTRUE);
  }
   
  // --- Beam pipe ---------------------------------------------------
  AliAlignObjParams* c1Align = 
    new((*structArray)[18]) AliAlignObjParams("CP1",dvoluid,
					      dXc, dYc, 0, rotX, rotY, 0,kTRUE);
  c1Align->Print();

  // --- C side Collar -----------------------------------------------
  AliAlignObjParams* cSupCAlign = 
    new((*structArray)[19]) AliAlignObjParams("CPSUPC",
					      dvoluid, 
					      tCpSupCC[0]-pCpSupCC[0], 
					      tCpSupCC[1]-pCpSupCC[1], 
					      tCpSupCC[2]-pCpSupCC[2], 
					      0, 0, 0, kTRUE);
  cSupCAlign->Print();
  
  // --- save in CDB storage -----------------------------------------
  Info("MakeForwardMisAlignment",
       "Saving alignment objects in CDB storage %s", ocdbPath);
  AliCDBStorage* storage = cdb->GetStorage(ocdbPath);
  if (!storage){
    Error("MakeForwardMisAlignment","Unable to open storage %s\n",ocdbPath);
    return;
  }


  // --- Make objects and store --------------------------------------
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Grosso Raffaele");
  md->SetComment("Misalignment for FRAME and beam pipe");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("GRP/Align/Data",0,AliCDBRunRange::Infinity());
  storage->Put(structArray,id,md);

  gGeoManager->Export("myGeom.root");
}
