Int_t AliITSFindClusters() {

  printf("FindClusters\n");

  TFile *in=TFile::Open("galice.root","UPDATE");
  if (!in->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 2;}

  in->ls();

   if (!(gAlice=(AliRun*)in->Get("gAlice"))) {
     cerr<<"gAlice have not been found on galice.root !\n";
     return 2;
   }


   gAlice->GetEvent(0);

   AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS"); 
   if (!ITS) {
     cerr<<"ITSFindClusters.C : AliITS object not found on file\n";
     return 3;
   }
   Int_t ver = ITS->IsVersion(); 
   cerr<<"ITS version "<<ver<<" has been found !\n";

// Set the models for cluster finding
   AliITSgeom *geom = ITS->GetITSgeom();

   // SPD
   AliITSDetType *iDetType=ITS->DetType(0);
   AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
   TClonesArray *dig0  = ITS->DigitsAddress(0);
   TClonesArray *recp0  = ITS->ClustersAddress(0);
   AliITSClusterFinderSPD *rec0=new AliITSClusterFinderSPD(seg0,dig0,recp0);
   ITS->SetReconstructionModel(0,rec0);
   // test
   printf("SPD dimensions %f %f \n",seg0->Dx(),seg0->Dz());
   printf("SPD npixels %d %d \n",seg0->Npz(),seg0->Npx());


   // SDD
   Float_t baseline = 10.;
   Float_t noise = 1.67;
   Int_t thres = baseline+3.*noise;

   AliITSDetType *iDetType=ITS->DetType(1);
   AliITSsegmentationSDD *seg1=(AliITSsegmentationSDD*)iDetType->GetSegmentationModel();
   if (!seg1) seg1 = new AliITSsegmentationSDD(geom);
   AliITSresponseSDD *res1 = (AliITSresponseSDD*)iDetType->GetResponseModel();
   if (!res1) res1=new AliITSresponseSDD();
	
	
   //res1->SetNoiseParam(noise,baseline);

   res1->SetNoiseParam(noise,baseline);
   Float_t magic = res1->MagicValue();
   Float_t top = res1->MaxAdc();
   thres *= top/magic;

   TClonesArray *dig1  = ITS->DigitsAddress(1);
   TClonesArray *recp1  = ITS->ClustersAddress(1);
   AliITSClusterFinderSDD *rec1=new AliITSClusterFinderSDD(seg1,res1,dig1,recp1);
   rec1->SetMinNCells(6);
   rec1->SetTimeCorr(70.);
   rec1->SetCutAmplitude(thres);
   ITS->SetReconstructionModel(1,rec1);


   // SSD
   AliITSDetType *iDetType=ITS->DetType(2);
   AliITSsegmentationSSD *seg2=(AliITSsegmentationSSD*)iDetType->GetSegmentationModel();
   TClonesArray *dig2  = ITS->DigitsAddress(2);
   TClonesArray *recp2  = ITS->ClustersAddress(2);
   AliITSClusterFinderSSD *rec2=new AliITSClusterFinderSSD(seg2,dig2,recp2);
   ITS->SetReconstructionModel(2,rec2);
   // test
   printf("SSD dimensions %f %f \n",seg2->Dx(),seg2->Dz());
   printf("SSD nstrips %d %d \n",seg2->Npz(),seg2->Npx());




   TStopwatch timer;

   switch (ver) {
   case 5:
      cerr<<"Looking for clusters...\n";
      {
	timer.Start();
	ITS->DigitsToRecPoints(0,1,"All");
      }
      break;
   default:
      cerr<<"Invalid ITS version !\n";
      return 5;
   }

   timer.Stop(); timer.Print();

   delete rec0;
   delete rec1;
   delete rec2;


   delete gAlice; gAlice=0;

   in->Close();

   return 0;
}
