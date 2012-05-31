void EvalAbso()
{
// Material numbers
    enum {kC=6, kAl=9, kFe=10, kCu=11, kW=12, kPb=13,
	  kNiCuW=21, kVacuum=16, kAir=15, kConcrete=17,
	  kPolyCH2=18, kSteel=10, kInsulation=14, kPolyCc=20};	   

    AliABSO *pABSO  = (AliABSO*) gAlice->GetModule("ABSO");

    Int_t nl[2];
    nl[0] = pABSO->NumberOfLayers(0);
    nl[1] = pABSO->NumberOfLayers(1);    
    
    Int_t   mat[2][15],  mat[2][15];
    Float_t zmin[2][15], zmin[2][15], zmax[2][15], zmax[2][15];
    Int_t i, j;

    for (j=0; j<   2; j++) {
	for (i=0; i< nl[j]; i++) {
	    mat[j][i]  = pABSO->MaterialOfLayer(j,i)-1599;
	    zmax[j][i] = pABSO->ZPositionOfLayer(j,i);
	    if (i+1 < nl[j]) zmin[j][i+1] = pABSO->ZPositionOfLayer(j,i);
	}
	zmin[j][0]=0.;
    }
    Float_t l = zmax[0][nl[0]-1];


//

//
// 1. Limiting angular resolution in the Branson formalism
//
    Float_t f0,f1,f2;
    Float_t a, z, dens, radl, absl;
    char  name[21];

    for (j=0; j< 2; j++) {
	printf("\n                        A            Z           ZPos         DZ            dens        radL        absL");
	printf("\n==================================================================================");
	f0=f1=f2=0;
	for (i=0; i< nl[j]; i++) {
	    pABSO->AliGetMaterial(mat[j][i], name, a, z, dens, radl, absl);
	    Float_t dz = zmax[j][i]-zmin[j][i];
	    f0 += dz/radl;
	    f1 += dz/(2.*radl)*(2.*zmin[j][i]+dz);
	    f2 += dz/(3.*radl)*(3.*zmin[j][i]*zmin[j][i]+3.*dz*zmin[j][i]+dz*dz);	
	    printf("\n %3d %14s %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f",
		   i+1, name, a, z, (zmax[j][i]+zmin[j][i])/2., dz, dens, radl, absl);
	}
//  
	Float_t deltaThetaB = 0.0136  * TMath::Sqrt(f0-f1*f1/f2);
	Float_t lBranson    = f2/f1;
	printf("\n==================================================================================");
	printf("\n Branson plane:                              %12.3f (cm)", lBranson);
	printf("\n Limiting resolution (Branson)             : %12.3f (mrad)", deltaThetaB*1000.);
    }
//
//  2. Limiting Resolution as a function of thickness of extra Fe
//
    Float_t zmin0[15], zmax0[15];
    Float_t x[100], y[100];

    for (i=0; i< nl[0]; i++) {
	zmin0[i] = zmin[0][i];
	zmax0[i] = zmax[0][i];	
    }
    
    Float_t ds = (zmax[0][3]-zmin[0][2])/100.;
    Int_t j;
    
    for (j=0; j< 100; j++) {
	f0=f1=f2=0;
	Float_t zFe = zmin[0][2]+Float_t(j)*ds;
	zmax0[2]=zFe;
	zmin0[3]=zFe;
	for (i=0; i< nl[0]; i++) {
	    Float_t a, z, dens, radl, absl;
	    pABSO->AliGetMaterial(mat[0][i], name, a, z, dens, radl, absl);
	    Float_t dz = zmax0[i]-zmin0[i];
	    f0 += dz/radl;
	    f1 += dz/(2.*radl)*(2.*zmin0[i]+dz);
	    f2 += dz/(3.*radl)*(3.*zmin0[i]*zmin0[i]+3.*dz*zmin0[i]+dz*dz);	
	}
	Float_t deltaThetaB = 0.0136  * TMath::Sqrt(f0-f1*f1/f2);
	x[j] = zmax[0][3]-zFe;
	y[j] = deltaThetaB*1000.;
    }
    TGraph *sGraph = new TGraph(100, x, y);
    TCanvas *c = new TCanvas("c","  ",400,10,600,700);
    c->Divide(2,2);
    c->cd(1);
    sGraph->Draw("AC");
    c->Update();
//
// 3. Limiting Resolution as a function of density of Carbon
//
    Float_t dRho = (2.5-1.5)/100.;

    for (j=0; j< 100; j++) {
	f0=f1=f2=0;
	Float_t rho = 1.5 + Float_t(j)*dRho;
	for (i=0; i< nl[0]; i++) {
	    Float_t a, z, dens, radl, absl;
	    pABSO->AliGetMaterial(mat[0][i], name, a, z, dens, radl, absl);
	    if (i==1) radl*=1.75/rho;
	    Float_t dz = zmax[0][i]-zmin[0][i];
	    f0 += dz/radl;
	    f1 += dz/(2.*radl)*(2.*zmin[0][i]+dz);
	    f2 += dz/(3.*radl)*(3.*zmin[0][i]*zmin[0][i]+3.*dz*zmin[0][i]+dz*dz);	
	}
	deltaThetaB = 0.0136  * TMath::Sqrt(f0-f1*f1/f2);
	x[j] = rho;
	y[j] = deltaThetaB*1000.;
    }
    TGraph *dGraph = new TGraph(100, x, y);
    c->cd(2);
    dGraph->Draw("AC");
    c->Update();

//
//   4. Energy Loss
//
//
    TCanvas *c1 = new TCanvas("c1","  ",400,10,600,700);
    c1->Divide(2,2);
 
    const Float_t kAvo=0.60221367;

    char* tmp;
    tmp = new char[5];
    strncpy(tmp, "LOSS", 4);
    tmp[4]='\0';

    Int_t ixst;
    Float_t ekin[100], dedx[100], de[100], pcut[100], deRad[100];
    Float_t emin =   2;    
    Float_t emax = 200;
    Float_t eps = (emax-emin)/100.;
//  3 < theta < 9
    for (j=0; j< 100; j++) {
	ekin[j] = emin + Float_t(j)*eps;
	de[j]=0;
	deRad[j]=0;
    }
// all losses    
    for (i=0; i< nl[0]; i++) {
	((TGeant3*) gMC)->Gftmat(pABSO->GetMatId(mat[0][i]), 5, tmp, 100, ekin, dedx, pcut, ixst);
	for (j=0; j< 100; j++) de[j]+=dedx[j]*(zmax[0][i]-zmin[0][i])/1000.;
    }
// radative losses
    for (i=0; i< nl[0]; i++) {
	Float_t a, z, dens, radl, absl;
	pABSO->AliGetMaterial(mat[0][i], name, a, z, dens, radl, absl);
	for (Int_t j=0; j<100; j++) {
	    Float_t nor= kAvo*dens/a*(zmax[0][i]-zmin[0][i]);
	    deRad[j] += (((TGeant3*) gMC)->Gbrelm(z,ekin[j],1.e10))*nor;
	    deRad[j] += (((TGeant3*) gMC)->Gprelm(z,ekin[j],1.e10))*nor;
	}
    }
    
    TH2F *hr = new TH2F("hr1","Several graphs in the same pad",2,0,15,2,0,5);
    hr->SetXTitle("X title");
    hr->SetYTitle("Y title");
    c1->cd(1);
    hr->Draw();

    TGraph *eGraph1  = new TGraph(100, ekin, de);
    TGraph *eGraphR1 = new TGraph(100, ekin, deRad);

    eGraph1 ->Draw("LP");
    eGraphR1->Draw("LP");

    hr = new TH2F("hr2","Several graphs in the same pad",2,15,200,2,0,5);
    c1->cd(2);
    hr->Draw();
    eGraph1 ->Draw("LP");
    eGraphR1->Draw("LP");
    

//  2 < theta < 3
    for (j=0; j< 100; j++) {
	ekin[j] = emin + Float_t(j)*eps;
	de[j]=0;
	deRad[j]=0;
    }
// all losses    
    for (i=0; i< nl[1]; i++) {
	((TGeant3*) gMC)->Gftmat(pABSO->GetMatId(mat[1][i]), 5, tmp, 100, ekin, dedx, pcut, ixst);
	for (j=0; j< 100; j++) de[j]+=dedx[j]*(zmax[1][i]-zmin[1][i])/1000.;
    }
// radative losses
    for (i=0; i< nl[1]; i++) {
	Float_t a, z, dens, radl, absl;
	pABSO->AliGetMaterial(mat[1][i], name, a, z, dens, radl, absl);
	for (Int_t j=0; j<100; j++) {
	    Float_t nor= kAvo*dens/a*(zmax[1][i]-zmin[1][i]);
	    deRad[j] += (((TGeant3*) gMC)->Gbrelm(z,ekin[j],1.e10))*nor;
	    deRad[j] += (((TGeant3*) gMC)->Gprelm(z,ekin[j],1.e10))*nor;
	}
    }

    TH2F *hr2 = new TH2F("hr3","Several graphs in the same pad",2,0,15,2,0,5);
    hr2->SetXTitle("X title");
    hr2->SetYTitle("Y title");
    c1->cd(3);
    hr2->Draw();

    TGraph *eGraph2  = new TGraph(100, ekin, de);
    TGraph *eGraphR2 = new TGraph(100, ekin, deRad);

    eGraph2 ->Draw("LP");
    eGraphR2->Draw("LP");

    hr2 = new TH2F("hr4","Several graphs in the same pad",2,15,200,2,0,5);
    c1->cd(4);
    hr2->Draw();
    eGraph2 ->Draw("LP");
    eGraphR2->Draw("LP");


    c->Update();

}




