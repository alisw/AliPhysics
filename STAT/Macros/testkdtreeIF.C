// #include <malloc.h>
// #include "TMatrixD.h"
// #include "TRandom.h"
// #include "TGraph.h"
// #include "TStopwatch.h"
// 


Int_t Mem()
{
  // size in 1000 bytes
  static struct mallinfo memdebug; 
  memdebug = mallinfo();
  return memdebug.uordblks/1000;
}



void testindexes(Int_t nloop, Int_t npoints);
void  testkdtreeIF(Int_t npoints=1000, Int_t nloop=1000, Int_t mode = 1, Int_t bsize=9);
void TestSizeIF(Int_t npoints=1000, Int_t nrows = 150, Int_t nsec=18, Int_t mode = 1, Int_t bsize=9);



void testindexes(Int_t nloop, Int_t npoints){
  //
  // test indexing 
  // 
  TKDTree<Int_t, Float_t> *kdtree = new TKDTree<Int_t, Float_t>(0,0,0,0);
  Int_t row =0;
  Int_t collumn =0; 
  TStopwatch timer;
  timer.Start();
  row = 10;
  for (Int_t iloop=0;iloop<nloop;iloop++)
  for (Int_t index=1;index<npoints;index++){
    //row = TMath::Log2(index);
    //row=0;
   //  if (index< (16<<row)) row=0;
//     for (; index>=(32<<row);row+=5);
//     for (; index>=(2<<row);row++);
//     collumn= index-(1<<row);
      TKDTree<Int_t, Float_t>::GetCoord(index,row,collumn);
    //
    //Int_t index2=kdtree->GetIndex(row,collumn);
    //printf("%d\t%d\t%d\t%d\n",index,index2,row,collumn);
    if (kdtree->GetIndex(row,collumn)!=index || collumn<0) {
      printf("Problem\n");
    }
  }
  timer.Stop();
  timer.Print();
}

void TestBuild(Int_t npoints, Int_t bsize){  
  Float_t *data0 =  new Float_t[npoints*2];
  Float_t *data[2];
  data[0] = &data0[0];
  data[1] = &data0[npoints];
  for (Int_t i=0;i<npoints;i++) {
    data[1][i]= gRandom->Rndm();
    data[0][i]= gRandom->Rndm();
    //data[1][i]= 0;
    //data[0][i]= -i;
  }
  Float_t dataf[2];
  TKDTreeIF *kdtree = new TKDTreeIF(npoints, 2, bsize, data);
  //kdtree->Build();
  for (Int_t i=0; i<npoints;i++){
    Int_t index = kdtree->fIndPoints[i];
    //printf("%d\t%d\t%f\n",i,index,data[0][index]);
  }
  for (Int_t i=0; i<kdtree->fNnodes;i++){
    //printf("%d\t%f\n",i,kdtree->fNodes[i].fValue);
  }
  Float_t sumiter = 0;
  for (Int_t i=0;i<npoints;i++){
    dataf[0] = data[0][i];
    dataf[1] = data[1][i];
    Int_t index=-1;
    Int_t iter =0;
    kdtree->FindPoint(dataf,index, iter);
    if (i!=index){
      printf("%d\t%d\t%f\t%f\n",i,index,dataf[0],data[0][index]);
    }
    sumiter+=iter;
  }
  printf("Mean iter = %f\n",float(sumiter)/float(npoints));
}


void TestSizeIF(Int_t npoints, Int_t nrows, Int_t nsec, Int_t mode, Int_t bsize)
{
  //
  // test size to build kdtree
  //
  Int_t before =Mem();
  for (Int_t isec=0; isec<nsec;isec++)
    for (Int_t irow=0;irow<nrows;irow++){
      testkdtreeIF(npoints,0,mode,bsize);
    }
  Int_t after = Mem();
  printf("Memory usage %d\n",after-before);
}


void  testkdtreeIF(Int_t npoints, Int_t nloop, Int_t mode, Int_t bsize)
{
  //
  // test speed and functionality of kdtree
  //
  Float_t rangey  = 100;
  Float_t rangez  = 100;
  Float_t drangey = 0.1;
  Float_t drangez = 0.1;
//   Float_t rangey  = 20;
//   Float_t rangez  = 250;
//   Float_t drangey = 1;
//   Float_t drangez = 1;

  //
  Float_t *data0 =  new Float_t[npoints*2];
  Float_t *data[2];
  data[0] = &data0[0];
  data[1] = &data0[npoints];
  Int_t i;   
  for (i=0; i<npoints; i++){
    data[0][i]          = gRandom->Uniform(-rangey, rangey);
    data[1][i]          = gRandom->Uniform(-rangez, rangez);
    //     data[i+npoints]  = TMath::Nint(gRandom->Uniform(-rangez, rangez)/10.)*10.;
    //printf("%d %f  %f\n", i, data[i], data[i+npoints]);     
  }
  TStopwatch timerbuild;
  TKDTree<Int_t, Float_t> *kdtree = new TKDTree<Int_t, Float_t>(npoints,2,bsize,data);
   kdtree->Build();
   timerbuild.Stop();
   timerbuild.Print();
  TStopwatch timer;

   Float_t countern=0;
   Float_t counteriter  = 0;
   Float_t counterfound = 0;


   if (mode ==2){ 
     if (nloop) timer.Start();
     Int_t res[npoints];
     Int_t nfound = 0;
     for (Int_t kloop = 0;kloop<nloop;kloop++){
       if (kloop==0){
	 counteriter = 0;
	 counterfound= 0; 
	 countern    = 0;
       }
       for (Int_t i=0;i<npoints;i++){
	 Float_t point[2]={data[0][i],data[1][i]};
	 Float_t delta[2]={drangey,drangez};
	 Int_t iter  =0;
	 nfound =0;
	 Int_t bnode =0;
	 //kdtree->FindBNode(point,delta, bnode);
	 //continue;
	 kdtree->FindInRangeA(point,delta,res,nfound,iter,bnode);
	 if (kloop==0){
	   //Bool_t isOK = kTRUE;
	   Bool_t isOK = kFALSE;
 	   for (Int_t ipoint=0;ipoint<nfound;ipoint++)
 	     if (res[ipoint]==i) isOK =kTRUE;
	   counteriter+=iter;
	   counterfound+=nfound;
	   if (isOK) {
	     countern++;
	   }else{
	     printf("Bug\n");
	   }
	 }
       } 
     }
     if (nloop){
       timer.Stop();
       timer.Print();
     }
   }
   delete [] data0;

   counteriter/=npoints;
   counterfound/=npoints;
   if (nloop) printf("Find nearest point:\t%f\t%f\t%f\n",countern, counteriter, counterfound); 
}
