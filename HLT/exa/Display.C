void Display(char *trackfile=0,Bool_t rpoints = kFALSE,Bool_t tpoints = kTRUE)
{
//  gROOT->Reset();
  Int_t minslice=0;
  Int_t maxslice=35;
  /*
    Display tracks and clusters already written to file.
  */
  char fname[256];
  UInt_t npp[36][5];
  AliL3SpacePointData *ppp[36][5];
  AliL3FileHandler *fff[36][5];
  for(int s=minslice;s<=maxslice;s++){
    for(int p=0;p<5;p++){
      ppp[s][p]=0;
      fff[s][p] = new AliL3FileHandler();
      sprintf(fname,"points_%d_%d.raw",s,p);
      if(!fff[s][p]->SetBinaryInput(fname)){
//        printf("no file %s\n",fname);
        delete fff[s][p];
        continue;
      }
      ppp[s][p] = (AliL3SpacePointData *)fff[s][p]->Allocate();
      fff[s][p]->Binary2Memory(npp[s][p],ppp[s][p]);
      fff[s][p]->CloseBinaryInput();
      printf("%s: %d points\n",fname,npp[s][p]);
    }
  }
  if(trackfile){
    AliL3FileHandler *trackf  = new AliL3FileHandler();
    if(!trackf->SetBinaryInput(trackfile)) return;
    AliL3TrackArray *tracks= new AliL3TrackArray();
    trackf->Binary2TrackArray(tracks);
    trackf->CloseBinaryInput();
    delete trackf;
    printf("%s: %d tracks\n",trackfile,tracks->GetNPresent());
  }

  //display the event in 3D
  TCanvas *c1 = new TCanvas("v1","",700,700);
  TView *v = new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);
  //v->SetRange(0,0,-430,430,560,1710);
  c1->Clear();
  c1->SetFillColor(1);
  c1->SetTheta(90.);
  c1->SetPhi(0.);

  if(rpoints){
    for(Int_t s=0;s<36;s++){
      for(int p=0;p<5;p++){
        AliL3SpacePointData *points = ppp[s][p];
        if(!points) continue;
        int npoints = npp[s][p];
        TPolyMarker3D *pm = new TPolyMarker3D(npoints);
  
        float xyz[3];
        for(int i=0; i<npoints; i++){
          xyz[0] = points[i].fX;
          xyz[1] = points[i].fY;
          xyz[2] = points[i].fZ;
           
          pm->SetPoint(i,xyz[0],xyz[1],xyz[2]); 
          if(xyz[0]==0){
            printf("%d %d %d %d\n",npoints,s,p,i);
          }
        }
        pm->SetMarkerColor(2);
        pm->Draw("");
      }
    }
  }

  if(trackfile)
    {
      int ntracks = tracks->GetNTracks();
      TPolyLine3D *line = new TPolyLine3D[ntracks];
      Float_t xcl[174];
      Float_t ycl[174];
      Float_t zcl[174];
      
      printf("ntracks %d\n",ntracks);
      
      for(int j=0; j<ntracks; j++)
	{
	  AliL3Track *gtrack = tracks->GetCheckedTrack(j); 
          if(!gtrack) continue;        
	  Int_t nHits = gtrack->GetNHits();
	  UInt_t *hitnum=gtrack->GetHitNumbers();
	  //printf("nhits %d\n",nHits);
	  if(nHits < 10) continue;
	  TPolyMarker3D *pm = new TPolyMarker3D(nHits);
	  for(int h=0; h<nHits; h++)
	    {
	      UInt_t id=hitnum[h];
              Int_t slice = (id>>25) & 0x7f;
              Int_t patch = (id>>22) & 0x7;
              UInt_t pos = id&0x3fffff;	      
//              printf("slice: %d patch: %d #%d \n",slice,patch,pos);

              AliL3SpacePointData *points = ppp[slice][patch];

              if(!points) {
                printf("*** No Points at slice: %d patch: %d #%d \n",
                      slice,patch,pos);
                continue;
              }
	      if(pos>=npp[slice][patch]) {printf("Error \n"); continue;}
	      //printf("hit %d %d\n",curclust->GetHitNumber(),id);
	      xcl[h] = points[pos].fX;
	      ycl[h] = points[pos].fY;
	      zcl[h] = points[pos].fZ;
	      pm->SetPoint(h,xcl[h],ycl[h],zcl[h]);
	    }
	  pm->SetMarkerColor(4);
	  if(tpoints) pm->Draw();
	  TPolyLine3D *current_line = &(line[j]);
	  current_line = new TPolyLine3D(nHits,xcl,ycl,zcl,"");
	  
	  current_line->SetLineColor(4);
	  current_line->Draw("same");
	}
    }
  TFile *geofile = new TFile("alice.geom","READ");  
  TGeometry *geom=(TGeometry*)geofile->Get("AliceGeom");
  geom->Draw("same");
  
  c1->cd();
  
    c1->x3d();
  // c1->Modified();
  //c1->Update();
  
  return;
}
