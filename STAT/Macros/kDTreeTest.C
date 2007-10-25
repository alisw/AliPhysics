#include <malloc.h>
#include "TMatrixD.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "../src/TKDTree.h"


Float_t Mem()
{
  // size in KB
  static struct mallinfo memdebug; 
  memdebug = mallinfo();
  return memdebug.uordblks/1024.;
}


void TestBuild(const Int_t npoints = 1000000, const Int_t bsize = 100);
void TestSpeed(const Int_t npower2 = 20, const Int_t bsize = 10);

void testindexes(Int_t nloop, Int_t npoints);
void  testkdtreeIF(Int_t npoints=1000, Int_t bsize=9, Int_t nloop=1000, Int_t mode = 2);
void TestSizeIF(Int_t npoints=1000, Int_t nrows = 150, Int_t nsec=18, Int_t mode = 1, Int_t bsize=9);

//______________________________________________________________________
void kDTreeTest()
{
	printf("\n\tTesting kDTree memory usage ...\n");
	TestBuild();
	
	printf("\n\tTesting kDTree speed ...\n");
	TestSpeed();
}

//______________________________________________________________________
void TestBuild(const Int_t npoints, const Int_t bsize){  
  Float_t *data0 =  new Float_t[npoints*2];
  Float_t *data[2];
  data[0] = &data0[0];
  data[1] = &data0[npoints];
  for (Int_t i=0;i<npoints;i++) {
    data[1][i]= gRandom->Rndm();
    data[0][i]= gRandom->Rndm();
  }
  Float_t before =Mem();
  TKDTreeIF *kdtree = new TKDTreeIF(npoints, 2, bsize, data);
	Float_t after = Mem();
  printf("Memory usage %f KB\n",after-before);
	delete kdtree;
  Float_t end = Mem();
  printf("Memory leak %f KB\n", end-before);
	return;	
}

//______________________________________________________________________
void TestSpeed(const Int_t npower2, const Int_t bsize)
{
	if(npower2 < 10){
		printf("Please specify a power of 2 greater than 10\n");
		return;
	}
	
	Int_t npoints = Int_t(pow(2., npower2))*bsize;
	Float_t *data0 =  new Float_t[npoints*2];
  Float_t *data[2];
  data[0] = &data0[0];
  data[1] = &data0[npoints];
  for (Int_t i=0;i<npoints;i++) {
    data[1][i]= gRandom->Rndm();
    data[0][i]= gRandom->Rndm();
  }

	TGraph *g = new TGraph(npower2-10);
	g->SetMarkerStyle(7);
	TStopwatch timer;
	Int_t tpoints;
  TKDTreeIF *kdtree = 0x0;
	for(int i=10; i<npower2; i++){
		tpoints = Int_t(pow(2., i))*bsize;
		timer.Start(kTRUE);
		kdtree = new TKDTreeIF(tpoints, 2, bsize, data);
		timer.Stop();
		g->SetPoint(i-10, i, timer.CpuTime());
		printf("npoints [%d] nodes [%d] cpu time %f [s]\n", tpoints, kdtree->GetNNodes(), timer.CpuTime());
		//timer.Print("u");
		delete kdtree;
	}
	g->Draw("apl");
	return;
}

//______________________________________________________________________
void TestSizeIF(Int_t npoints, Int_t nrows, Int_t nsec, Int_t mode, Int_t bsize)
{
  //
  // test size to build kdtree
  //
  Float_t before =Mem();
  for (Int_t isec=0; isec<nsec;isec++)
    for (Int_t irow=0;irow<nrows;irow++){
      testkdtreeIF(npoints,0,mode,bsize);
    }
  Float_t after = Mem();
  printf("Memory usage %f\n",after-before);
}


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


void  testkdtreeIF(Int_t npoints, Int_t bsize, Int_t nloop, Int_t mode)
{
//
// Test speed and functionality of 2D kdtree.
// Input parametrs:
// npoints - number of data points
// bsize   - bucket size
// nloop   - number of loops
// mode    - tasks to be performed by the kdTree
//         - 0  : time building the tree
//

 
  Float_t rangey  = 100;
  Float_t rangez  = 100;
  Float_t drangey = 0.1;
  Float_t drangez = 0.1;

  //
  Float_t *data0 =  new Float_t[npoints*2];
  Float_t *data[2];
  data[0] = &data0[0];
  data[1] = &data0[npoints];
  Int_t i;   
  for (i=0; i<npoints; i++){
    data[0][i]          = gRandom->Uniform(-rangey, rangey);
    data[1][i]          = gRandom->Uniform(-rangez, rangez);
  }
  TStopwatch timer;

	// check time build
	printf("building kdTree ...\n");
	timer.Start(kTRUE);
  TKDTreeIF *kdtree = new TKDTreeIF(npoints, 2, bsize, data);
  timer.Stop();
  timer.Print();
	if(mode == 0) return;
	


	Float_t countern=0;
	Float_t counteriter  = 0;
	Float_t counterfound = 0;


	if (mode ==2){
		if (nloop) timer.Start(kTRUE);
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
