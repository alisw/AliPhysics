void raw2treeGrid_collection()
{

  // reading RAW data from test LCS
  // filling histograms
  // fillinf tree
  //  gROOT->LoadMacro("loadlibs.C");
  //  loadlibs();
  Int_t allData[220][5];
  TGrid::Connect("alien://");
  TTree *fT0OutTree=new TTree("t0tree","None here");
 TAlienCollection *collnum = TAlienCollection::Open("wn.xml");
  Int_t numrun;
  collnum->Reset();
  collnum->Next();
  TString buf_runnum;
  TString buf_path = collnum->GetTURL() ;
  for(int i=0; i<buf_path.Length();i++)	{
    if(buf_path(i,4)=="/raw")	{
      buf_runnum = buf_path(i-6,6);
      numrun = buf_runnum.Atoi();
      break;
    }
  }
  
  TString names[220];
  Int_t chvalue[220], meanchvalue[220];
  AliT0LookUpKey* lookkey= new AliT0LookUpKey();
  AliT0LookUpValue*  lookvalue= new AliT0LookUpValue();
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(numrun);
  AliT0Parameters *fParam = AliT0Parameters::Instance();
  fParam->Init(); 
  TMap *lookup = fParam->GetMapLookup();
  TMapIter *iter = new TMapIter(lookup);
  for( Int_t iline=0; iline<212; iline++)
    {
      lookvalue = ( AliT0LookUpValue*) iter->Next();
      lookkey = (AliT0LookUpKey*) lookup->GetValue((TObject*)lookvalue);      
      if(lookkey){
	Int_t key=lookkey->GetKey();
	names[key]=lookkey->GetChannelName();
	fT0OutTree->Branch(Form("%s",names[key].Data()), &chvalue[key]);
      }
      else
	{printf(" no such value %i \n", iline);}
    } 
  Float_t meanCFD[24], meanQT1[24];
 
  for (int ich=0; ich<24; ich++) {
    meanCFD[ich] = fParam->GetCFD(ich);
    meanQT1[ich] = fParam->GetQT1(ich);
  }
  Float_t meanOrA = fParam->GetMeanOrA();
  Float_t meanOrC = fParam->GetMeanOrC();
  Float_t meanTVDC = fParam->GetMeanVertex();
  //new QTC

  Float_t qt01mean[28] = {18712.5, 18487.5, 18487.5, 18537.5, 
  			  18562.5, 18462.5, 18537.5, 18537.5, 
			  18537.5, 18587.5, 18587.5, 18512.5,
			   18512.5, 18512.5, 18487.5, 18562.5,
			    18537.5, 18512.5, 18537.5, 18537.5,
			     18512.5, 18587.5, 18562.5, 18512.5,
		      18358, 18350, 18374, 18362};
  Float_t qt11mean[28] = {18705, 18495, 18465, 18555, 
  			18555, 18435, 18525, 18525, 
			18525, 18585, 18585, 18495, 
			18495, 18525, 18465, 18555, 
			18525, 18495, 18555, 18495, 
			18495, 18585, 18585, 18495,
		      18358, 18350, 18374, 18362};
  Int_t ind[26];
  for (int iii=0; iii<12; iii++) ind[iii]=25;
  for (int iii=12; iii<24; iii++) ind[iii]=57;

  UInt_t event;
   fT0OutTree->Branch("event", &event);
   ULong64_t triggerMask;
   fT0OutTree->Branch("triggers", &triggerMask);
    TAlienCollection *coll = TAlienCollection::Open("wn.xml");
      coll->Reset();
    
   AliRawReader *reader;
while (coll->Next()) {
    TString	fFileName=coll->GetTURL();
     //READ DATA
 //     TString	fFileName=Form("alien:///alice/data/2015/LHC15i/000%i/raw/15000%i028.%i.root", numrun, numrun, chunk);
      reader = new AliRawReaderRoot(fFileName);
     if(!reader) continue;

     reader = new AliRawReaderRoot(fFileName);     
     if(!reader) continue;
//     reader->LoadEquipmentIdsMap("T0map.txt");
     reader->RequireHeader(kTRUE);
     for (Int_t i0=0; i0<220; i0++) {
       chvalue[i0] = 0;
       for (Int_t j0=0; j0<5; j0++)  allData[i0][j0]=0; 
     }
     
     AliT0RawReader *start = new AliT0RawReader(reader);
     while (reader->NextEvent()) {
       start->Next();
       for (Int_t ii=0; ii<211; ii++) {
	 chvalue[ii] = 0;
	 for (Int_t iHit=0; iHit<5; iHit++) 
	   {
	     allData[ii][iHit]= start->GetData(ii,iHit);
	     //   	if(allData[ii][iHit]>0) cout<<ii<<" "<<allData[ii][iHit]<<endl;
	   }
       } 
       
       const  UInt_t type =reader->GetType();
       if(type != 7) continue;
       triggerMask  = reader->GetClassMask();
       
        for (Int_t iHit=0; iHit<5; iHit++) {
         if( allData[50][iHit]>meanTVDC-800 && allData[50][iHit]<meanTVDC+800) {
	 chvalue[50]=allData[50][iHit];
         break;
         }
       } 
       
        for (Int_t in=0; in<24;  in++)
	 {
	   for (Int_t iHit=0; iHit<5; iHit++)  //old QTC C side
	     {
	       if (allData[2*in+ind[in]+1][iHit] > meanQT1[in]-800 &&  
		   allData[2*in+ind[in]+1][iHit] < meanQT1[in]+800 ) {
		 chvalue[2*in+ind[in]+1] = allData[2*in+ind[in]+1][iHit];
		 break;
	       }
	     }
	   for (Int_t iHit=0; iHit<5; iHit++)  //old QTC A side
	     {
	       if( (allData[2*in+ind[in]][iHit] > chvalue[2*in+ind[in]+1]) &&
		   chvalue[2*in+ind[in]+1]>0)
		 {
		   chvalue[2*in+ind[in]] = allData[2*in+ind[in]][iHit];
		//   printf("index %i  pmt %i QTC old start %i stop %i \n", 
		//   2*in+ind[in], in,
		 //  chvalue[2*in+ind[in]+1], chvalue[2*in+ind[in]]);
		   break;
		 }
	     }
	 }
       for (Int_t in=0; in<12; in++)  
	 {
	   chvalue[in+68+1] = allData[in+68+1][0] ;
	   chvalue[in+12+1] = allData[in+12+1][0] ;
	   for (Int_t iHit=0; iHit<5; iHit++)  //CFD C side
	     {
	       if(allData[in+1][iHit] > meanCFD[in]-800 && 
		  allData[in+1][iHit] < meanCFD[in]+800)
		 {
		   chvalue[in+1] = allData[in+1][iHit] ; 
		   break;
		 }
	     }
	   for (Int_t iHit=0; iHit<5; iHit++)  //CFD A side
	     {
	       if(allData[in+1+56][iHit]>0)
		 if(allData[in+1+56][iHit] > meanCFD[in+12]-800 && 
		    allData[in+1+56][iHit] < meanCFD[in+12]+800)
		   {
		     chvalue[in+1+56] = allData[in+56+1][iHit] ;
		     break;
		   }
	     }
	 }
       // new QTC    
       Int_t pmt;
       for (Int_t ik=0; ik<106; ik+=4) {
	   if (ik<48)          pmt=ik/4;
	   if (ik>47 && ik<52) pmt= 24;   
	   if (ik>51 && ik<56) pmt= 25; 
	   if(ik>55)           pmt=(ik-8)/4;
	   for(Int_t iHt = 0; iHt<5; iHt++) { 
	     if(allData[107+ik+1][iHt] > (qt01mean[pmt]-800) &&
		allData[107+ik+1][iHt] < (qt01mean[pmt]+800) ) {
	       chvalue[107+ik+1] = allData[107+ik+1][iHt];
	    //   printf("start newQTC 00 ik %i iHt %i pmt %i  QT00 %i QT01 %i \n", ik, iHt, pmt, allData[107+ik][iHt],  allData[107+ik+1][iHt]);
	       break;
	     }
	   }
	   for(Int_t iHt = 0; iHt<5; iHt++) { 
	     if(allData[107+ik][iHt]>chvalue[107+ik+1] &&
		chvalue[107+ik+1]>0) {
	       chvalue[107+ik]=allData[107+ik][iHt] ;
	    //   printf("stop newQTC 00 ik %i iHt %i pmt %i  QT00 %i QT01 %i \n", ik, iHt, pmt, allData[107+ik][iHt],  allData[107+ik+1][iHt]);
	       break;
	     }
	   }
	   for(Int_t iHt = 0; iHt<5; iHt++) {
	     if( allData[107+ik+3][iHt] > (qt11mean[pmt]-800) &&
		 allData[107+ik+3][iHt] < (qt11mean[pmt]+800) ) {
	       chvalue[107+ik+3] = allData[107+ik+3][iHt];
	       break;
	     }
	   }
	   for(Int_t iHt = 0; iHt<5; iHt++) {
	     if( allData[107+ik+2][iHt] > chvalue[107+ik+3]&&
		 chvalue[107+ik+3]>0 ) {
	       chvalue[107+ik+2] = allData[107+ik+2][iHt];
	   //    printf(" newQTC 11 ik %i iHt %i pmt %i QT10 %i QT11 %i \n", ik, iHt, pmt, allData[107+ik+2][iHt],  allData[107+ik+3][iHt]);
		 break;
	     }
	   }
	 } //end new QTC
	 // Or
         for(Int_t iHt = 0; iHt<5; iHt++) {
	 if(allData[51][iHt]>meanOrA-800 && allData[51][iHt]<meanOrA+800) {
	   chvalue[51]=allData[51][iHt];
	   break;
	  }
         }
        for(Int_t iHt = 0; iHt<5; iHt++) {
	 if(allData[52][iHt]>meanOrC-800 && allData[52][iHt]<meanOrC+800) {
	   chvalue[52]=allData[52][iHt];
	   break;
	  }
         }
	   	    

       event++;
       if(chvalue[50]>0)      fT0OutTree->Fill(); 
       
     } //event
     start->Delete();
     }
     reader->Delete();
	
   TFile *hist = new TFile("T0RAWtree.root","RECREATE");
   hist->cd();
   fT0OutTree ->Write();  
			   
}
