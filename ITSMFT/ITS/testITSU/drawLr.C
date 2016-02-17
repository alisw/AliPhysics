// draw the 1st sensor of every ladder of requested layers (axis directions - thin arrows) 
// and the tracking frame of each sensor (axis directions - thick arrows)

void drawLr(int layMin=1,int layMax=1) 
{

  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  gSystem->Load("libITSUpgradeRec");

  AliGeomManager::LoadGeometry("geometry.root");


  // Apply misaligment ... ;-)
  const char *ocdb="local://$ALICE_ROOT/OCDB";
  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  Int_t run = 1;
  AliCDBManager::Instance()->SetRun(run); 
  AliCDBManager::Instance()->SetSpecificStorage("ITS/Align/Data", Form("local://%s",gSystem->pwd()));
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("ITS/Align/Data");
  TClonesArray *array = (TClonesArray*)entry->GetObject();
  AliGeomManager::ApplyAlignObjsToGeom(*array);
  gGeoManager->LockGeometry();

  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE,kTRUE);
  TObjArray* segmArr = gm->GetSegmentations();  
  //
  int nlr = gm->GetNLayers();
  if (layMin<0) layMin = 0;
  if (layMax<0) layMax = nlr-1;
  else if (layMax>=nlr) layMax = nlr-1;
  //
  TH2* hh = new TH2F("hh","hh",100,-70,70,100,-70,70);
  hh->Draw();
  gStyle->SetOptStat(0);
  //
  double rmax = 0;
  TGeoHMatrix *msens=0, *mttr=0;
  double  loc[3]={0,0,0},gloC[3],glo[3],glo1[3],glo2[3],trk[3];
  //
  Int_t mod=0;

  for (Int_t lay=layMin;lay<=layMax;lay++) {
    //
    AliITSMFTSegmentationPix* segm = (AliITSMFTSegmentationPix*)segmArr->At(gm->GetLayerChipTypeID(lay));
    for (int ild=0;ild<gm->GetNLadders(lay);ild++) {
    
      // Sensor Matrices
      printf("Layer %d Ladder %d\n",lay,ild);
      msens = gm->GetMatrixSens(lay,ild,mod);
      printf("Sensor Matrix: ");
      msens->Print();
      mttr = gm->GetMatrixT2L(lay,ild,mod);
      printf("T2L Matrix: ");
      mttr->Print();
      //
      loc[0]=loc[1]=loc[2] = 0;
      msens->LocalToMaster(loc,gloC);
      mttr->MasterToLocal(loc,trk);  
      printf("SensorCenter in Lab: ");
      for (int i=0;i<3;i++) printf("%+9.4f ",gloC[i]); printf("\n");
      printf("SensorCenter in TRK: ");
      for (int i=0;i<3;i++) printf("%+9.4f ",trk[i]); printf("\n");
      //
      double r = TMath::Sqrt(gloC[0]*gloC[0]+gloC[1]*gloC[1]);
      if (rmax<r) rmax = r;
      //
      loc[0]=-segm->Dx()/2;
      msens->LocalToMaster(loc,glo);
      loc[0]= segm->Dx()/2;;
      msens->LocalToMaster(loc,glo1);
      TArrow* linS = new TArrow(glo[0],glo[1],glo1[0],glo1[1],0.012);  // sensor with pos local X direction 
      linS->SetLineColor((ild+1)%4+1); linS->Draw(); 
      //
      loc[0]=0;
      loc[1]=rmax*0.1;
      msens->LocalToMaster(loc,glo2);
      //
      //lin->Print(); 
      linS = new TArrow(gloC[0],gloC[1],glo2[0],glo2[1],0.012);  // pos local Y axis
      linS->SetLineColor((ild+1)%4+1); linS->Draw(); 
      //
      TMarker* mrk = new TMarker(gloC[0],gloC[1],31);
      mrk->SetMarkerColor((ild+1)%4+1);
      mrk->Draw();
      //
      TLatex *latx = new TLatex( gloC[0]*1.2,gloC[1]*1.2,Form("%d",ild));
      latx->SetTextColor((ild+1)%4+1);
      latx->SetTextSize(0.02);
      latx->Draw();
      //
      // Check Tracking to Local Matrix -------
      //
      // Draw sensors tracking frame
      trk[0]=trk[1]=0;
      mttr->LocalToMaster(trk,loc); 
      msens->LocalToMaster(loc,glo1);
      TMarker* mrk = new TMarker(glo1[0],glo1[1],24);
      mrk->SetMarkerColor((ild+1)%4+1);
      mrk->Draw();
      //
      // normal to sensor plane
      TLine* linN = new TLine(0,0,glo1[0],glo1[1]); 
      linN->SetLineWidth(1);
      linN->SetLineStyle(2);
      linN->SetLineColor((ild+1)%4+1); 
      linN->Draw(); 
      //
      // connect tracking and local frame
      TLine* linNP = new TLine(gloC[0],gloC[1],glo1[0],glo1[1]); 
      linNP->SetLineWidth(1);
      linNP->SetLineStyle(2);
      linNP->SetLineColor((ild+1)%4+1); 
      linNP->Draw(); 
      //
       // direction of X axis of the tracking plane
      trk[0]=rmax*0.1;
      mttr->LocalToMaster(trk,loc); 
      msens->LocalToMaster(loc,glo); 
      TArrow* linX = new TArrow(glo1[0],glo1[1],glo[0],glo[1],0.012);
      linX->SetLineWidth(2);
      linX->SetLineStyle(1);
      linX->SetLineColor((ild+1)%4+1); 
      linX->Draw(); 
      //
      trk[0] = 0;
      // direction of tracking Y axis
      double dx = glo[0]-glo1[0];
      double dy = glo[1]-glo1[1];
      double dst = TMath::Sqrt(dx*dx+dy*dy);
      trk[1]=-dst;
      mttr->LocalToMaster(trk,loc); 
      msens->LocalToMaster(loc,glo);
      trk[1]= dst;
      mttr->LocalToMaster(trk,loc); 
      msens->LocalToMaster(loc,glo1);
      //
      TArrow* linT = new TArrow(glo[0],glo[1],glo1[0],glo1[1],0.012); // normal to sensor plane, pox Y of tracking frame
      linT->SetLineWidth(2);
      linT->SetLineStyle(1);
      linT->SetLineColor((ild+1)%4+1); 
      linT->Draw(); 
      //
      printf("Layer %d Ladder %d\n",lay,ild);
      mttr->Print();
    }
  }
  //
  rmax = 1.3*rmax;
  hh->GetXaxis()->SetRangeUser(-rmax,rmax);
  hh->GetYaxis()->SetRangeUser(-rmax,rmax);
  gPad->Modified();
  gPad->Update();
}
