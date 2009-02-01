void plotField(Int_t iField = 0)
{
//
//  
//  load necessary libraries
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libminicern");
    gSystem->Load("$(ROOTSYS)/lib/libPhysics");
    gSystem->Load("$(ROOTSYS)/lib/libEG");
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libSTEER");
    
//
//
//  create field map

     AliMagF* field = new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG);

//     field-SetL3ConstField(1);
     
//
//  get parameters
     Float_t xMin, xMax, yMin, yMax, zMin, zMax;
     Float_t dX, dY, dZ;
     xMin = -350.;
     xMax =  350.;
     dX   = (field->FieldMap(0))->DelX();
     yMin = -350.;
     yMax =  350.;
     dY   = (field->FieldMap(0))->DelY();
     zMin =   250.;
     zMax =  1450.;
     dZ   = (field->FieldMap(0))->DelZ();
     Int_t nx = (xMax-xMin)/dX;
     Int_t ny = (yMax-yMin)/dY;
     Int_t nz = (zMax-zMin)/dZ;
     
	 
//
// create histogram
     TH2F* hMap = new TH2F("hMap", "Field Map y-z", 
			   nz, zMin, zMax, ny, yMin, yMax);
     TH2F* hMap1 = new TH2F("hMap1", "Field Map y-z", 
			   nz, zMin, zMax, ny, yMin, yMax);
     TH1F* hZ   = new TH1F("hZ", "Field along Z", 
			   nz, zMin, zMax);

     TH2F* hVec = new TH2F("hVec", "Field Map y-z", 
			   nz, zMin, zMax, ny, yMin, yMax);

     TH2F* hCir = new TH2F("hCir", "Field Map y-z", 
			   nz, zMin, zMax, ny, yMin, yMax);
     Float_t bMax = 0.;
     for (Int_t i = 0; i < nz; i++) {
	 for (Int_t j = 0; j < ny; j++) {
	     Float_t x[3];
	     Float_t b[3];
	     
	     x[2] = zMin + i * dZ;
	     x[1] = yMin + j * dY;
	     x[0] = 0.;
	     field->Field(x, b);
	     Float_t bb = TMath::Sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
	     if (bb > bMax) bMax = bb;
	     hMap->Fill(x[2], x[1], bb);
	 }
     }

     for (Int_t i = 0; i < nz; i++) {
	 for (Int_t j = 0; j < ny; j++) {
	     x[2] = zMin + i * dZ;
	     x[1] = yMin + j * dY;
	     x[0] = 0.;
	     field->Field(x, b);
	     Float_t bb = TMath::Sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
	     Float_t db = (bMax-bb)/bMax*100.;
	     
	     hMap1->Fill(x[2], x[1], db);
	 }
     }
     
     for (Int_t i = 0; i < nz; i++) {
	 x[2] = zMin + i * dZ +dZ/2.;
	 x[1] = 0.;
	 x[0] = 0.;
	 field->Field(x, b);
	 Float_t bb = TMath::Sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
	 hZ->Fill(x[2], bb);
     }

     TCanvas *c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
     hMap->Draw();
     TCanvas *c2 = new TCanvas("c2","Canvas 2",400,10,600,700);
     hZ->Draw();
     TCanvas *c3 = new TCanvas("c3","Canvas 3",400,10,600,700);
     hVec->Draw();
     
     Float_t scale1 = 0.9*TMath::Sqrt(dZ*dZ+dY*dY)/bMax/2.;
     Float_t scale2 = 0.005/bMax;
     TArrow* arrow;
     
     
     for (Int_t i = 0; i < nz; i++) {
	 for (Int_t j = 0; j < nx; j++) {
	     x[2] = zMin + i * dZ + dZ/2.;
	     x[0] = xMin + j * dX + dX/2.;
	     x[1] = 0.;
	     field->Field(x, b);
	     Float_t bb = TMath::Sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
	     b[2] *= scale1;
	     b[0] *= scale1;
	     Float_t width = 0.005+scale2*b[1];
	     arrow = new TArrow(x[2], x[0], x[2]+b[2], x[0]+b[0], width,"|>");
	     arrow->SetFillColor(1);
	     arrow->SetFillStyle(1001);
	     arrow->Draw();
	     c3->Modified();
	     c3->cd();
	 }
     }

     TCanvas *c4 = new TCanvas("c4","Canvas 4",400,10,600,700);
     hCir->Draw();
     for (Int_t i = 0; i < nz; i++) {
	 for (Int_t j = 0; j < ny; j++) {
	     x[2] = zMin + i * dZ + dZ/2.;
	     x[1] = yMin + j * dY + dY/2.;
	     x[0] = 0.;
	     field->Field(x, b);
	     Float_t bb = TMath::Sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
	     TEllipse *ellipse; 
	     ellipse= new TEllipse(x[2], x[1], b[1]*scale1/4., b[1]*scale1/4.,0,360,0);
	     ellipse->Draw();
	     c4->Modified();
	     c4->cd();
	 }
     }
}



