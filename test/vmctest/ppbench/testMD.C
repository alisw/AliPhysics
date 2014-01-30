// $Id$
//
// Macro for checking mother-daughter particles relationship in stack.
// By A. Morsch, 24/10/2013

void testMD(Int_t evNumber1 = 0, Int_t evNumber2 = 0) 
{
// Connect the Root Galice file containing Geometry, Kine and Hits

    AliRunLoader* rl = AliRunLoader::Open("galice.root");
//
    TDatabasePDG*  DataBase = new TDatabasePDG();
    
//
//   Loop over events 
//
    rl->LoadKinematics();
    rl->LoadHeader();    
    for (Int_t nev = evNumber1; nev <= evNumber2; nev++) {
	rl->GetEvent(nev);
	AliStack* stack = rl->Stack();
	Int_t npart = stack->GetNtrack();
	if (nev < evNumber1) continue;
	printf("Event %5d \n", nev);
//
// Particle loop
//       
	
	for (Int_t part = 0; part < npart; part++) {
	  TParticle *particle = stack->Particle(part);
	  Int_t pdg  = particle->GetPdgCode();
	  Int_t child1 = particle->GetFirstDaughter();
	  Int_t child2 = particle->GetLastDaughter();
	  if (child1 == -1 || part < 8) continue;
	  
	  for (Int_t i = child1; i<= child2; i++) 
	    {
	      TParticle *daughter = stack->Particle(i);
	      Int_t mo = daughter->GetFirstMother();
	      Int_t dodg = daughter->GetPdgCode();
	      if (mo != part) 
		printf("Particle with index %5d (pdg %5d) has daughters from index %5d to %5d, however, daughter %5d (pdg %5d) points back to %5d \n", 
		       part, pdg, child1, child2, i, dodg, mo);
	    }
	} // primary loop
    }
}







