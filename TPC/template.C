void MakeClusterTree(Int_t n)
{
  TFile * f = new TFile("pokus.root","recreate");
  AliTPCClustersArray arr;
  arr.MakeArray(10000);
  arr.MakeTree();
  for (Int_t i=0;i<n;i++) {
    AliTPCClustersRow * row =new AliTPCClustersRow;
    row->SetIndex(i);
    arr.AddSegment(row);
  }
  for (Int_t i=0;i<n;i++) arr.StoreSegment(n-i);
  arr.GetTree()->Write("pokus1");

 
}


void MakeTree(Int_t n)
{
  AliSegmentArray arr;

  TFile * f= new  TFile("pokus.root","recreate");
  arr.MakeArray(10000);
  arr.MakeTree();
  //  for (Int_t i=0;i<n;i++) arr.AddSegment(new AliSegment(Int_t((gRandom->Rndm())*n)));

  //for (Int_t i=0;i<n;i++) arr.StoreSegment(Int_t((gRandom->Rndm())*n));
  for (Int_t i=0;i<n;i++) arr.AddSegment(new AliSegment(i));
  for (Int_t i=0;i<n;i++) arr.StoreSegment(n-i);
  arr.GetTree()->Write("pokus1");
}

void ConnectTree(Int_t n,AliSegmentArray *a)
{
  AliSegmentArray &arr= *a;
  //  TFile * f = new TFile("pokus.root","update");
  TFile * f = new TFile("pokus.root","update");
  arr.MakeArray(10000);
  arr.ConnectTree("pokus1");
  for (Int_t i=0;i<n;i++) arr.LoadSegment(i);  
  for (Int_t i=0;i<n;i++)
    {
      if (arr[i]==0) continue;
      if (arr[i]->GetID()!=i) cout<<i<<"\n";  
    }
}

