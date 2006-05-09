// $Header$

void reve_load_vsd(const Text_t* vsd = "AliVSD.root")
{
  gReve->GetSelector()->LoadVSD(vsd);
  gReve->GetSelector()->SelectHits();
}
