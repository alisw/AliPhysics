void scanTemp(const char *fn = "2015/Run235840_235840_v1_s0.root",
	      const int kNumSens = 160, // active number of sensors
	      const int kNSM = 20) 
{
  TFile *file0 = TFile::Open(fn);
  file0->ls();
  AliCDBEntry->Print();
  AliEMCALSensorTempArray *arr = AliCDBEntry->GetObject();
  //  arr->Print();
  cout << " NumSensors " << arr->NumSensors()
       << " GetFirstIdDCS() " << arr->GetFirstIdDCS()
       << " GetLastIdDCS() " << arr->GetLastIdDCS()
       << endl;
  file0->Close();

  // info for each sensor
  int np[kNumSens] = {0};
  double min[kNumSens] = {0};
  double max[kNumSens] = {0};
  int id = 0;

  for (int isensor=0; isensor<kNumSens; isensor++) {
    AliEMCALSensorTemp *o = arr->GetSensor(isensor);
    if (o) {
      id = o->GetIdDCS();
      o->Print();
      cout << " side " << o->GetSide()
	   << " sector " << o->GetSector()
	   << " num " << o->GetNum()
	   << " id " << id
	   << " startTime " << o->GetStartTime()
	   << " endTime " << o->GetEndTime()
	   << endl;
      if (id >= kNumSens) {
	cout << " exiting - invalid id " << id << endl;
	exit();
      }
      
      AliSplineFit *f = o->GetFit();
      np[id] = 0;
      min[id] = 100;
      max[id] = 0;
      
      if (f) {
	np[id] = f->GetKnots();
	cout << " np " << np[id] << endl;
	Double_t *x = f->GetX();
	Double_t *y0 = f->GetY0();
	Double_t *y1 = f->GetY1();
	for (int i=0; i<np[id]; i++) {
	  cout << " i " << i
	       << " x " << x[i]
	       << " y0 " << y0[i]
	       << " y1 " << y1[i]
	       << endl;
	  
	  if (min[id]>y0[i]) min[id]=y0[i];
	  if (max[id]<y0[i]) max[id]=y0[i];
	}
      }
    }
  }

  for (int iSM=0; iSM<kNSM; iSM++) {
    cout << " iSM " << iSM << endl;
    double aveMin = 0;
    double aveMax = 0;
    int nOK = 0;
    for (int is=0; is<8; is++) {
      id = is + iSM*8;
      printf("id %d np %d min %3.2f max %3.2f diff %3.2f\n",
	     id, np[id], min[id], max[id],
	     max[id] - min[id]);

      if (max[id] > 1) {
	aveMin += min[id];
	aveMax += max[id];
	nOK++;
      }
    }

    if (nOK > 0) {
      aveMin /= nOK;
      aveMax /= nOK;
    }
    printf("iSM %d nOK %d average min %3.2f max %3.2f (max+min)/2 %3.2f\n",
	   iSM, nOK, aveMin, aveMax, (aveMin + aveMax)/2.);
  }

}

