{
   TFile *file = new TFile("galice.root");
    file->ls();
	delete gAlice;
	gAlice = (AliRun*)(file->Get("gAlice"));

	Int_t npart = gAlice->GetEvent(0);
	assert(npart);
	TObjArray*  parray = gAlice->Particles();
	cout<<" load_particles: N part="<<npart<<endl;

}











