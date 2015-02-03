#ifndef __CINT__
#include "alles.h"
#include "AliRunDigitizer.h"
#include "AliTPCDigitizer.h"
#include "AliH2F.h"
#include "TTree.h"


#endif

/// \file AliTPCTestMerge.C
///
/// test of the tpc merging using AliRunDigitizer and
/// TPC  Hits2Digits, Hits2SDigits and SDigits2Digits macros
/// preparation
/// 0. make 2 directorys - ev1 and ev2
/// 1.make hits, digits,sdigits and sdigits-digits in both directories
///  1.a  galice -b -q grun.C and produce hits
///  1.b. cp galice.root galice.root.hits
///  1.c  run AliTPCHits2Digits.C
///  1.d  cp galice.root galice.root.digits
///  1.e  copy back cp galice.root.hits galice.root
///  1.f  run AliTPCSDigits2Digits.C
///  1.g  cp galice.root  galice.root.sdigits
///  1.h  run AliTPCSDigits2Digit.C
///  1.i  cp galice.root galice.root.dig2
/// 2. cp ev1/galice.root/galice.root.sdigit galice.root
/// 3. load this macro and run testmerge()
/// 4. run test function bellow to compare merged digits with original one
/// 5. to be  noticed output ftom testx() function - should be bigger than
///    noise because in proces of digitisation we use different seed
///    of random numbers for diffusion and gas gain
///    -anyway in place where we have the signal should occur the signal in both casses
///    - only the amplitude should be different - by factor of sqrt
///
/// \author Marian Ivanov

void testmerge()
{
  /// merge two example events
  ///
  /// it merge two events -one from current directory -second from directory ev2

  if(gAlice) delete gAlice;
  AliRunDigitizer * manager = new AliRunDigitizer(2,1);
  manager->SetTreeDTPCBaseName("TreeD_75x40_100x60_150x60_");
  manager->SetInputTreeTPCSBaseName("TreeS_75x40_100x60_150x60_");
  manager->SetInputStream(0,"galice.root");
  manager->SetInputStream(1,"ev2/galice.root.sdigits");
  AliTPCDigitizer *dTPC = new AliTPCDigitizer(manager);
  manager->SetNrOfEventsToWrite(1); 
  TStopwatch timer;
  timer.Start();

  manager->Exec(""); 
  timer.Stop();
  timer.Print();

}

void drawmerged(Int_t sec, Int_t row, Int_t x1=-1, Int_t x2=-1, Int_t y1=-1, Int_t y2=-1)
{
  /// if you think that there is memory leak -
  /// you are tru but othervise graphic doesn't work
  /// sec=0; row =0;

  TFile * f = new TFile("galice.root");
  TFile * f1= new TFile("ev1/galice.root.digits");
  TFile * f2= new TFile("ev2/galice.root.digits");
  TTree * tree = (TTree*)f->Get("TreeD_75x40_100x60_150x60_0");
  TTree * tree1 = (TTree*)f1->Get("TreeD_75x40_100x60_150x60_0");
  TTree * tree2 = (TTree*)f2->Get("TreeD_75x40_100x60_150x60_0");
  //
  AliSimDigits *dig=0;
  AliSimDigits *dig1=0;
  AliSimDigits *dig2=0;
  //
  tree->GetBranch("Segment")->SetAddress(&dig);
  tree1->GetBranch("Segment")->SetAddress(&dig1);
  tree2->GetBranch("Segment")->SetAddress(&dig2);
  AliTPCParam * param =(AliTPCParam*) f->Get("75x40_100x60");
  if(param){
    delete param;
    param=new AliTPCParamSR();
  }
  else param=(AliTPCParam*) f->Get("75x40_100x60_150x60");
  Int_t index = param->GetIndex(sec,row);
  tree->GetEvent(index);
  tree1->GetEvent(index);
  tree2->GetEvent(index);

  TCanvas * c = new TCanvas(" Test merged digits", "test",600,900);
  c->Divide(1,3);
  //
  c->cd(1);
  AliH2F * his = dig->DrawDigits("cont1",x1,x2,y1,y2);
  his->SetTitle("MergedDigits");
  his->SetName("Merged Digits");
  his->GetXaxis()->SetTitle("time");
  his->GetYaxis()->SetTitle("pad");
  his->DrawClone("cont1");
  //
  c->cd(2);
  AliH2F * his1 =dig1->DrawDigits("cont1",x1,x2,y1,y2); 
  his1->SetTitle("background");
  his1->SetName("background");
  his1->GetXaxis()->SetTitle("time");
  his1->GetYaxis()->SetTitle("pad");
  his1->DrawClone("cont1"); 
  //
  c->cd(3);
  AliH2F * his2 =dig2->DrawDigits("cont1",x1,x2,y1,y2); 
  his2->SetTitle("signal");
  his2->SetName("signal");
  his2->GetXaxis()->SetTitle("time");
  his2->GetYaxis()->SetTitle("pad");
  his2->DrawClone("cont1");  
}


void drawd(TFile * f, Int_t amp1, Int_t amp2)
{
  TTree * tree = (TTree*)f->Get("TreeD_75x40_100x60_150x60_0");
  AliSimDigits *dig=0;
  tree->GetBranch("Segment")->SetAddress(&dig); 
  TH1F * his = new TH1F("his","his",amp2-amp1,amp1,amp2);
  for (Int_t i=0;i<60;i++){  
    tree->GetEvent(i);   
    dig->ExpandBuffer();
    Int_t nrows = dig->GetNRows();
    Int_t ncols = dig->GetNCols();
    for (Int_t rows=0;rows<nrows; rows++)
      for (Int_t col=0;col<ncols; col++){    
	Int_t d  = dig->GetDigitFast(rows,col);
	his->Fill(d);
      }
  }
  his->Draw();
}

void test1(){
  /// test of the merged digits
  /// compare merged digits with standard digits

  TFile f("galice.root");
  TFile f1("ev1/galice.root.digits");
  TFile f2("ev2/galice.root.digits");
  TTree * tree = (TTree*)f.Get("TreeD_75x40_100x60_150x60_0");
  TTree * tree1 = (TTree*)f1.Get("TreeD_75x40_100x60_150x60_0");
  TTree * tree2 = (TTree*)f2.Get("TreeD_75x40_100x60_150x60_0");
  //
  AliSimDigits *dig=0;
  AliSimDigits *dig1=0;
  AliSimDigits *dig2=0;
  //
  tree->GetBranch("Segment")->SetAddress(&dig);
  tree1->GetBranch("Segment")->SetAddress(&dig1);
  tree2->GetBranch("Segment")->SetAddress(&dig2);
  //
  for (Int_t i=0;i<6000;i++){
    tree->GetEvent(i);
    tree1->GetEvent(i);
    tree2->GetEvent(i);
    dig->ExpandBuffer();
    dig1->ExpandBuffer();
    dig2->ExpandBuffer();
    //
    Int_t nrows = dig->GetNRows();
    Int_t ncols = dig->GetNCols();
    for (Int_t rows=0;rows<nrows; rows++)
      for (Int_t col=0;col<ncols; col++){
	Int_t d  = dig->GetDigitFast(rows,col);
	Int_t d1 = dig1->GetDigitFast(rows,col);
	Int_t d2 = dig2->GetDigitFast(rows,col);
	
	if (abs(d-(d1+d2))>4)
	  printf("%d\t%d\t%d\t%d\t%d\t%d\n",i,rows,col,d,d1,d2);
      }
  }
}

void test5(){
  /// compare merged digits with digits obtained hits2sdig->sdigtodig

  TFile f("galice.root");
  TFile f1("ev1/galice.root.dig2");
  TFile f2("ev2/galice.root.dig2");
  TTree * tree = (TTree*)f.Get("TreeD_75x40_100x60_150x60_0");
  TTree * tree1 = (TTree*)f1.Get("TreeD_75x40_100x60_150x60_0");
  TTree * tree2 = (TTree*)f2.Get("TreeD_75x40_100x60_150x60_0");

  AliSimDigits *dig=0;
  AliSimDigits *dig1=0;
  AliSimDigits *dig2=0;
  
  tree->GetBranch("Segment")->SetAddress(&dig);
  tree1->GetBranch("Segment")->SetAddress(&dig1);
  tree2->GetBranch("Segment")->SetAddress(&dig2);




  for (Int_t i=0;i<6000;i++){
    tree->GetEvent(i);
    tree1->GetEvent(i);
    tree2->GetEvent(i);
    dig->ExpandBuffer();
    dig1->ExpandBuffer();
    dig2->ExpandBuffer();
    
    Int_t nrows = dig->GetNRows();
    Int_t ncols = dig->GetNCols();

    for (Int_t rows=0;rows<nrows; rows++)
      for (Int_t col=0;col<ncols; col++){
	Int_t d  = dig->GetDigitFast(rows,col);
	Int_t d1 = dig1->GetDigitFast(rows,col);
	Int_t d2 = dig2->GetDigitFast(rows,col);
	
	if (abs(d-d1-d2)>4)
	//if (d2>5)
	  printf("%d\t%d\t%d\t%d\t%d\t%d\n",i,rows,col,d,d1,d2);
      }
  }
}

void test3(){
  /// test of the merged digits

  TFile f("galice.root");
  TFile f1("ev1/galice.root.sdigits");
  TFile f2("ev2/galice.root.sdigits");
  TTree * tree = (TTree*)f.Get("TreeD_75x40_100x60_150x60_0");
  TTree * tree1 = (TTree*)f1.Get("TreeS_75x40_100x60_150x60_0");
  TTree * tree2 = (TTree*)f2.Get("TreeS_75x40_100x60_150x60_0");
  //
  AliSimDigits *dig=0;
  AliSimDigits *dig1=0;
  AliSimDigits *dig2=0;
  //  
  tree->GetBranch("Segment")->SetAddress(&dig);
  tree1->GetBranch("Segment")->SetAddress(&dig1);
  tree2->GetBranch("Segment")->SetAddress(&dig2);
  //
  for (Int_t i=0;i<6000;i++){
    tree->GetEvent(i);
    tree1->GetEvent(i);
    tree2->GetEvent(i);
    if ( (dig1->GetID()!=i) || (dig2->GetID()!=i) ) {
      printf("missed segments\n");
    }
    //
    dig->ExpandBuffer();
    dig1->ExpandBuffer();
    dig2->ExpandBuffer();
    //
    Int_t nrows = dig->GetNRows();
    Int_t ncols = dig->GetNCols();
    //
    for (Int_t rows=0;rows<nrows; rows++)
      for (Int_t col=0;col<ncols; col++){
	Int_t d  = dig->GetDigitFast(rows,col);
	Int_t d1 = dig1->GetDigitFast(rows,col)/16;
	Int_t d2 = dig2->GetDigitFast(rows,col)/16;	
	if (abs(d-d1-d2)>4)
	//if (d2>5)
	  printf("%d\t%d\t%d\t%d\t%d\t%d\n",i,rows,col,d,d1,d2);
      }
  }
}


void TestSDigitsDig2(){
  /// test of the digits produced by the Hits2Digits
  /// and Hits2SDigits - SDigits2Digits chain

  TFile f1("galice.root.digits");
  TFile f2("galice.root.dig2");
  //
  TTree * tree1 = (TTree*)f1.Get("TreeD_75x40_100x60_150x60_0");
  TTree * tree2 = (TTree*)f2.Get("TreeD_75x40_100x60_150x60_0");
  //
  AliSimDigits *dig1=0;
  AliSimDigits *dig2=0;
  //  
  tree1->GetBranch("Segment")->SetAddress(&dig1);
  tree2->GetBranch("Segment")->SetAddress(&dig2);
  //
  for (Int_t i=0;i<6000;i++){
    //tree->GetEvent(i);
    tree1->GetEvent(i);
    tree2->GetEvent(i);
    //dig->ExpandBuffer();
    if ( (dig1->GetID()!=i) || (dig2->GetID()!=i) ) {
      printf("miised semgnets\n");
    }
    //        
    dig1->ExpandBuffer();
    dig2->ExpandBuffer();
    dig1->ExpandTrackBuffer();
    dig2->ExpandTrackBuffer();
    //
    Int_t nrows = dig1->GetNRows();
    Int_t ncols = dig1->GetNCols();
    //
    for (Int_t rows=0;rows<nrows; rows++)
      for (Int_t col=0;col<ncols; col++){
	Int_t d1 = dig1->GetDigitFast(rows,col);
	Int_t d2 = dig2->GetDigitFast(rows,col);
	Int_t t1_1 =dig1->GetTrackIDFast(rows,col,0);
	Int_t t1_2 =dig1->GetTrackIDFast(rows,col,1);
	Int_t t1_3 =dig1->GetTrackIDFast(rows,col,2);
	//
	Int_t t2_1 =dig2->GetTrackIDFast(rows,col,0);
	Int_t t2_2 =dig2->GetTrackIDFast(rows,col,1);
	Int_t t2_3 =dig2->GetTrackIDFast(rows,col,2);
	//
	if ( (d2>0) && (d1>0))
	  if ( ( ( d2>2) || (d1>2)) && (t2_1!=t1_1)) {
	    printf("%d\t%d\t%d\t\t%d\t%d\n",i,rows,col,t1_1,t2_1);
	    printf("\t\t\t\t%d\t%d\n",d1,d2);
	  }
      }
  }
}

void TestSDigitsDig1(){
  /// test of the digits produced by the Hits2Digits
  /// and Hits2SDigits - SDigits2Digits chain

  TFile f1("galice.root.digits");
  TFile f2("galice.root.dig2");
  //
  TTree * tree1 = (TTree*)f1.Get("TreeD_75x40_100x60_150x60_0");
  TTree * tree2 = (TTree*)f2.Get("TreeD_75x40_100x60_150x60_0");
  //
  AliSimDigits *dig1=0;
  AliSimDigits *dig2=0;
  //  
  tree1->GetBranch("Segment")->SetAddress(&dig1);
  tree2->GetBranch("Segment")->SetAddress(&dig2);
  //
  for (Int_t i=0;i<6000;i++){
    //tree->GetEvent(i);
    tree1->GetEvent(i);
    tree2->GetEvent(i);
    //dig->ExpandBuffer();
    if ( (dig1->GetID()!=i) || (dig2->GetID()!=i) ) {
      printf("miised semgnets\n");
    }
    //        
    dig1->ExpandBuffer();
    dig2->ExpandBuffer();
    dig1->ExpandTrackBuffer();
    dig2->ExpandTrackBuffer();
    //
    Int_t nrows = dig1->GetNRows();
    Int_t ncols = dig1->GetNCols();
    //
    for (Int_t rows=0;rows<nrows; rows++)
      for (Int_t col=0;col<ncols; col++){
	//	Int_t d  = dig->GetDigitFast(rows,col);
	Int_t d1 = dig1->GetDigitFast(rows,col);
	Int_t d2 = dig2->GetDigitFast(rows,col);
 
	//	
	if ( ((d2>4) || (d1>4))  && abs(d1-d2)>4)
 	  printf("%d\t%d\t%d\t\t%d\t%d\n",i,rows,col,d1,d2);
      }
  }
}



void test4(){
  /// TPC internal test

  TFile f1("galice.root.sdigits");
  TFile f2("galice.root.digits");
  TTree * tree1 = (TTree*)f1.Get("TreeS_75x40_100x60_150x60_0");
  TTree * tree2 = (TTree*)f2.Get("TreeD_75x40_100x60_150x60_0");
  //
  AliSimDigits *dig1=0;
  AliSimDigits *dig2=0;
  //
  tree1->GetBranch("Segment")->SetAddress(&dig1);
  tree2->GetBranch("Segment")->SetAddress(&dig2);
  //
  for (Int_t i=0;i<6000;i++){
    //tree->GetEvent(i);
    tree1->GetEvent(i);
    tree2->GetEvent(i);
    //dig->ExpandBuffer();
    if ( (dig1->GetID()!=i) || (dig2->GetID()!=i) ) {
      printf("miised semgnets\n");
    }
    //        
    dig1->ExpandBuffer();
    dig2->ExpandBuffer();
    //
    Int_t nrows = dig1->GetNRows();
    Int_t ncols = dig1->GetNCols();
    //
    for (Int_t rows=0;rows<nrows; rows++)
      for (Int_t col=0;col<ncols; col++){
	Int_t d1 = dig1->GetDigitFast(rows,col)/16.;
	Int_t d2 = dig2->GetDigitFast(rows,col);		
	if ((d2>5) &&abs(d1-d2)>2) 	  printf("%d\t%d\t%d\t\t%d\t%d\n",i,rows,col,d1,d2);
      }
  }
}

