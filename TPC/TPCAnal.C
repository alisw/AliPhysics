/*
#include "alles.h"
extern AliTPCParam *gTPCParam;
extern AliTPCClustersArray * gCalcClusters;
extern AliTPCClustersArray * gDifClusters;
*/



TH1F * GetDX(Float_t x1=-2, Float_t x2=3, Float_t a1=-1, Float_t a2=1)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd = new TH1F("delta x", "dx",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);      
      //            if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      //      if ( (cl->fGenerTrack != cl->fTracks[0])) continue;
      if ( (cl->fAngleX<a2) && (cl->fAngleX>a1) )  hd->Fill(cl->fDx);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}

TH1F * GetDY(Float_t x1=-2, Float_t x2=3, Float_t a1=-1, Float_t a2=1)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd = new TH1F("delta y", "dy",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);    
      //if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;    
      //if ( (cl->fGenerTrack != cl->fTracks[0])) continue;
      //	if (cl->fNy!=2) continue;
      if ((cl->fAngleY<a2) && (cl->fAngleY>a1) )  hd->Fill(cl->fDy);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}



TH1F * GetNX(Float_t x1=0, Float_t x2=7, Float_t a1=-1, Float_t a2=1)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd = new TH1F("N x", "N x",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);      
      if ((cl->fAngleX<a2) && (cl->fAngleX>a1) )  hd->Fill(cl->fNx);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}

TH1F * GetNY(Float_t x1=-0, Float_t x2=7, Float_t a1=-1, Float_t a2=1)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd = new TH1F("N y", "N y",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);      
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      if ( (cl->fAngleY<a2) && (cl->fAngleY>a1) )  hd->Fill(cl->fNy);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}




TH1F * GetAngleX(Float_t x1=-2, Float_t x2=3)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hdx = new TH1F("Angle x", "Angle x",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);      
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hdx->Fill(cl->fAngleX);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hdx;
}

TH1F * GetAngleX2(Float_t x1=-2, Float_t x2=3)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hdx = new TH1F("Angle x2", "Angle x2",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);      
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hdx->Fill(cl->fAngleX*cl->fAngleX);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hdx;
}

TH1F * GetAngleY(Float_t x1=-2, Float_t x2=3)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hdx = new TH1F("Angle y", "Angle y",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);      
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hdx->Fill(cl->fAngleY);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hdx;
}

TH1F * GetAngleY2(Float_t x1=-2, Float_t x2=3)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hdx = new TH1F("Angle y2", "Angle y2",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);      
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hdx->Fill(cl->fAngleY*cl->fAngleY);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hdx;
}


TH1F * GetKmax(Float_t x1=0, Float_t x2=200)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd = new TH1F("Max ADC", "Max ADC",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);      
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fMax);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}


TH1F * GetQ(Float_t x1=0, Float_t x2=200)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd = new TH1F("Integral Q", "Integral Q",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);      
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fQ);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}


TH1F * GetCalcQ(Float_t x1=0, Float_t x2=1)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd  = new TH1F("kCalcQ", "kCalcQ",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
 if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fOrigQ);
    } 
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}


TH1F * GetCalcQSQR(Float_t x1=0, Float_t x2=1)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd  = new TH1F("kCalcQSQRt", "kCalcQSQRt",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
 if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(TMath::Sqrt(cl->fOrigQ));
    } 
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}



TH1F * GetQDivQ(Float_t x1=0, Float_t x2=1)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd  = new TH1F("kQdivCalcQ", "kQDivCalcQ",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
      //      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fQ/cl->fOrigQ);
    } 
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}



TH1F * GetKmaxDivQ(Float_t x1=0, Float_t x2=1)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd  = new TH1F("kMaxdivQ", "kMaxDivQ",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fMax/cl->fQ);
    } 
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}








TH2F * GetDXvRow(Float_t x1=0, Float_t x2=20, Float_t y1=0, Float_t y2=3)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH2F * hd = new TH2F("delta x versur row ", "DXvRow",33,x1,x2,33,y1,y2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    Int_t sector,row;
    gTPCParam->AdjustSectorRow(clrow->GetID(),sector,row);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);      
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(row,cl->fDx,1);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}


TH1F * GetSigmaX(Float_t x1=0, Float_t x2=2)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd = new TH1F("SigmaX", "SigmaX",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(TMath::Sqrt(cl->fSigmaX2));
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}


TH1F * GetSigmaY(Float_t x1=0, Float_t x2=2)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd = new TH1F("SigmaY", "SigmaY",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(TMath::Sqrt(cl->fSigmaY2));
    }
    arr->ClearSegment(clrow->GetID());
  }
  
  return hd;
}

TH1F * GetArea(Float_t x1=0, Float_t x2=40)
{
  AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH1F * hd = new TH1F("Cluster area", "cluster area",33,x1,x2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fArea);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}


TH2F * GetSigmaXvX(Float_t x1=0, Float_t x2=500,Float_t y1=0, Float_t y2=1)
{
   AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH2F * hd = new TH2F("SigmaXvX", "SigmaXvX",33,x1,x2,33,y1,y2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fX,TMath::Sqrt(cl->fSigmaX2),1);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}

TH2F * GetSigmaYvX(Float_t x1=0, Float_t x2=500,Float_t y1=0, Float_t y2=1)
{
   AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH2F * hd = new TH2F("SigmaYvX", "SigmaYvX",33,x1,x2,33,y1,y2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
 if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fX,TMath::Sqrt(cl->fSigmaY2),1);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}

TH2F * GetDeltaXvX(Float_t x1=0, Float_t x2=500,Float_t y1=0, Float_t y2=3)
{
   AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH2F * hd = new TH2F("DeltaXvX", "DeltaXvX",33,x1,x2,33,y1,y2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
 if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fX,cl->fDy,1);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}

TH2F * GetDeltaYvX(Float_t x1=0, Float_t x2=500,Float_t y1=0, Float_t y2=0.5)
{
   AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH2F * hd = new TH2F("DeltaYvX", "DeltaYvX",33,x1,x2,33,y1,y2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
 if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fX,cl->fDy,1);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}

TH2F * GetAreavX(Float_t x1=0, Float_t x2=500,Float_t y1=0, Float_t y2=20)
{
   AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH2F * hd = new TH2F("AreavX", "AreavX",33,x1,x2,33,y1,y2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
 if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fX,cl->fArea,1);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}


TH2F * GetKvX(Float_t x1=0, Float_t x2=500,Float_t y1=0, Float_t y2=20)
{
   AliTPCClustersArray *arr=gDifClusters;
  if (arr==0) return 0;  
  TH2F * hd = new TH2F("KvX", "KvX",33,x1,x2,33,y1,y2);
  Int_t nsegment = (Int_t)arr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    AliTPCClustersRow *clrow= (AliTPCClustersRow*)arr->LoadEntry(segment);
    if (clrow==0) continue;
    Int_t ncluster = clrow->GetArray()->GetEntriesFast();
    for (Int_t i = 0; i<ncluster;i++){
      AliDifCluster *cl =(AliDifCluster*)clrow->GetArray()->UncheckedAt(i);
      if ( (cl->fGenerTrack != cl->fTracks[0]) || (cl->fTracks[1]>0)) continue;
      hd->Fill(cl->fX,cl->fMax,1);
    }
    arr->ClearSegment(clrow->GetID());
  }
  return hd;
}

