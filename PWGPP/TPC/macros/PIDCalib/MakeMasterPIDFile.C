#include "TSpline.h"
#include "TFile.h"
#include "TF1.h"
#include "TObjArray.h"
#include "TKey.h"


 
TSpline3* newReducedSpline(TObject *spo);

void MakeMasterPIDFile()
{
  TObjArray arr;

  // TODO: Add the desired files here as shown for the following examples:
  arr.Add(new TNamed("/hera/alice/bhess/gridOutput/pp/8TeV/LHC12a.pass2/splines_12a.pass2.root","")); 
  arr.Add(new TNamed("/hera/alice/bhess/gridOutput/pp/8TeV/LHC12b.pass2/splines_12b.pass2.root",""));
  
  TIter nextFile(&arr);

  TObject *o=0x0;
  TObjArray pid;
  
  while ( (o=nextFile()) ){
    TFile f(o->GetName());
    if (!f.IsOpen()||f.IsZombie()){
      printf("\nERROR opening file: %s\n",o->GetName());
      continue;
    }
    printf("file: %s (%s)\n",o->GetName(),o->GetTitle());
    TList *l=f.GetListOfKeys();
    TIter nextKey(l);
    TKey *key=0x0;
    while ( (key=(TKey*)nextKey()) ) {
      TString name=key->GetName();
      TString cont=o->GetTitle();
      if (!cont.IsNull()) if(cont!="x"&&!name.Contains(cont)) continue;
      TObject *o2=key->ReadObj();
      
      TSpline3 *newsp=newReducedSpline((TSpline3*)o2);
      if (cont=="x") newsp->SetNameTitle(TString(newsp->GetName()).ReplaceAll("LHC10D","LHC10E"),TString(newsp->GetName()).ReplaceAll("LHC10D","LHC10E"));
      
      TString newSpName(newsp->GetName());
      
      if (newSpName.Contains("MC")) newsp->SetNameTitle(newSpName.ReplaceAll("PASS2","PASS1"),newSpName.ReplaceAll("PASS2","PASS1"));
      
      pid.Add(newsp);
      delete o2;
      printf("%s (%d)\n",newsp->GetName(),newsp->GetNp());
    }


    
  }

  printf("splines done!\n");
  
  pid.SetOwner();
  
  TFile f2("TPCPIDResponse.root","recreate");
  pid.Write("TPCPIDResponse",TObject::kSingleKey);
  f2.Close();
  
  printf("done all...\n");
}


TSpline3* newReducedSpline(TObject *spo){

  TSpline3 *sp=0x0;

  if ( (spo->IsA()==TSpline3::Class()) ) sp=(TSpline3*)spo;
  if ( (spo->IsA()==TGraph::Class()) ) {
    sp=new TSpline3(spo->GetName(),(TGraph*)spo);
    sp->SetNameTitle(TString(spo->GetName()).ReplaceAll("TGRAPH","TSPLINE3"),TString(spo->GetName()).ReplaceAll("TGRAPH","TSPLINE3"));
  }
  
  Int_t np=sp->GetNp();
  if (np<=500) return (TSpline3*)sp->Clone(sp->GetName());

  Int_t nth=(Int_t)((Double_t)np/500.);

  TGraph gr;

  Int_t ipnew=0;
  for (Int_t ip=0;ip<np;++ip){
    Double_t x=0,y=0;
    sp->GetKnot(ip,x,y);
    if (!(ip%nth)) gr.SetPoint(ipnew++,x,y);
  }
  TSpline3 *spNew=new TSpline3("xxx",&gr);
  TString title(sp->GetName());
  title.ToUpper();
  spNew->SetNameTitle(title.Data(),title.Data());

  return spNew;
}

