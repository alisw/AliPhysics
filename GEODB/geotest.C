/**********************************************/
/*                                            */
/* FILE: geotest.C                            */
/* PURPOSE: To "experiment" with geometric's  */
/*          databases and ROOT.               */
/* LANGUAGE: C++                              */
/* COMPILER: CC for HP-UX 9.x and 10.         */
/* AUTHOR: Joana && David                     */
/* DATE: December 2, 1998                     */
/* ADDRESS: jesanto@cern.ch, dcollado@cern.ch */
/*                                            */
/**********************************************/

geotest () {
    gSystem->Load("./libGEODB");
    AliGBox* box1 = new AliGBox("box1", "box1", 0.2, 0.2, 0.2);
    AliGMaterial *mat1 = new AliGMaterial("mat1", "mat1", 2.3,4.5,6.7);
    AliGNode *node1 = new AliGNode("Node1", 1, "My first node", box1, mat1 );

    AliGBox* box2 = new AliGBox("box2", "box2", 0.3,0.3,0.3);

    AliGTube* tube = new AliGTube("tube", "tube", 0.150,0.200,0.400);
    AliGMaterial *mat2 = new AliGMaterial("mat2", "mat2", 1.3,7.5,8.7);
    AliGNode *node2 = new AliGNode("Node2", 2, "My second node", /*tube*/ box2, mat2);

    //AliGBox* box3 = new AliGBox("box3", "box3", 3.0,2.0,4.0);
    AliGSphere* sphere = new AliGSphere("sphere", "sphere", 0.25);
    AliGMaterial *mat3 = new AliGMaterial("mat3", "mat3", 0.3,3.5,2.7);
    //AliGNode *node3 = new AliGNode("Node3", 3, "My third node", box3, mat3);
    AliGNode *node3 = new AliGNode("Node3", 3, "My third node", sphere, mat3);

    //AliGBox* box4 = new AliGBox("box4", "box4", 4.0,2.0,6.0);
    AliGCone* cone = new AliGCone("cone", "cone", 0.05,0.07,0.12,0.15,0.1);
    AliGMaterial *mat4 = new AliGMaterial("mat4", "mat4", 0.4,3.2,8.7);
    //AliGNode *node4 = new AliGNode("Node4", 4, "My fourth node", box4, mat4);
    AliGNode *node4 = new AliGNode("Node4", 4, "My fourth node", cone, mat4);

    AliGTransform *tran1 = new AliGTransform("tran1","tran1","tra 1. 0. 0.");
    AliGTransform *tran2 = new AliGTransform("tran2","tran2","tra 0.5 0. 0.");
    AliGTransform *tran3 = new AliGTransform("tran3","tran3","tra 0.25 0. 0.");

    node1->Add(node2,tran1);       // TUBE
    node1->Add(node2,tran2);//
    node1->Add(node3,tran1);//
    node1->Add(node3,tran2);       // SPHERE
    node2->Add(node4,tran1);//
    node2->Add(node4,tran2);       // CONE
    node2->Add(node4,tran3);//

    node1->AddConfig("config1","config1","coarse",981201,991201);
    node2->AddConfig("config2","config2","detail",981201,991201);
    node3->AddConfig("config3","config3","",981201,991201);
    node4->AddConfig("config4","config4","",981201,991201);
    
    AliGeometry* geom1 = new AliGeometry( "geom1", "geom1", node1, 0 ); /* Build geometry in memory from Node1 */
    
    TFile* file = new TFile( "AliGeoDB.root", "RECREATE","Alice Geometry Database");
    
    node1->SaveAll( file ); /* Memory tree structure saved in disk */
    
    
    file->Close();

    TFile* file = new TFile( "AliGeoDB.root", "UPDATE" );

    AliGeometry* geom1 = new AliGeometry( "geom1", "geom1", node1, 0 ); /* Build geometry in memory from Node1 */

    geom1->Write(); /* Save geometry in disk */
    file->Write();
   
    file->Close();  /* Everything saved in disk */
    file = new TFile("AliGeoDB.root","read");
     a = TBrowser() ;
    printf( " The end of geotest.C \n" );

    //if( geom1 ) delete geom1;
    //if( node4 ) delete node4;
    //if( node3 ) delete node3;
    //if( node2 ) delete node2;
    //if( node1 ) delete node1;


 
    //AliGeometry* geom3  = new AliGeometry();
 
    /* Retrieve the memory tree structure from AliGeoDB.root and stores it in memory, being new_tree the top node */
    //AliGNode* tree_root = new AliGNode(geom3->FileToMemTree( "AliGeoDB.root", "geom1" ));

    //TFile* file2 = new TFile( "AliGeoDB2.root", "RECREATE","Alice Geometry Database");
    //tree_root->SaveAll( file2 ); /* Memory tree structure saved in disk (second time) */

    
    //TFile* file = new TFile( "AliGeoDB.root", "READ" );
    
    //AliGeometry* geom2 = new AliGeometry((AliGeometry*) file->Get("geom1"));
    //file->Close();
    //file2->cd();
    //geom2.Write();
    //file2->Write();
    //file2->Close();

    //gROOT->Reset();
    //file3 = new TFile( "AliGeoDB.root", "READ" );
    //file4 = new TFile( "AliGeoDB2.root", "READ" );
    //a = TBrowser();

    //gPad->Update();









    /*    AliGNode *obj;

    TIter next(gPad->GetListOfPrimitives());
    while (obj = (AliGNode*)next()) {
        if(obj.GetShape()->Is3D()) {
            obj.GetShape()->Sizeof3D();
	    printf("\n It is a 3D Object\n");
            printf("Total size of x3d primitives:\n");
            Size3D size = obj.GetShape()->GetSize3D();
            printf("     gSize3D.numPoints=%d\n", size.numPoints);
            //printf("     gSize3D.numSegs  =%d\n", gPad->gSize3D.numSegs);
            //printf("     gSize3D.numPolys =%d\n", gSize3D.numPolys);
        }
        else
	    printf("\n It's not a 3D Object\n");
    }

    gPad->x3d();*/

    /*
    c1 = new TCanvas("c1","Common Shapes.",0,0,1202,846);
    c1->Range(0,0,1,1);
    c1->SetFillColor(32); // Light Green
    c1->SetBorderSize(3);
    c1->SetBorderMode(0); // -1 (down) 0 (no) 1 (up)

    // Create a TitleBar in Main Window
    title = new TPaveLabel(0.3,0.93,0.7,0.98,"View Of Different Shapes");
    title->SetFillColor(34); // Grey
    title->Draw();

    // Text for Pads
    TText *t = new TText();
    t->SetTextFont(32); // Light Green
    t->SetTextColor(1); // Black
    t->SetTextSize(0.03);
    t->DrawText( 0.19, 0.89, "AliGBox" );
    t->DrawText( 0.70, 0.89, "AliGSphere" );
    t->DrawText( 0.19, 0.43, "AliGTube" );
    t->DrawText( 0.70, 0.43, "AliGCone" );

    // Define and Draw Pads of Main Window
    Pad1 = new TPad("Pad1","AliGBox Shape",0.05,0.48,0.45,0.88,28);
    Pad1->SetFillColor(34); // Grey
    Pad2 = new TPad("Pad2","AliGSphere Shape",0.55,0.48,0.95,0.88,28);
    Pad2->SetFillColor(34); // Grey
    Pad3 = new TPad("Pad3","AliGTube Shape",0.05,0.02,0.45,0.42,28);
    Pad3->SetFillColor(34); // Grey
    Pad4 = new TPad("Pad4","AliGCone Shape",0.55,0.02,0.95,0.42,28);
    Pad4->SetFillColor(34); // Grey
    Pad1->Draw();
    Pad2->Draw();
    Pad3->Draw();
    Pad4->Draw();

    //  Define some volumes
    AliGTube* tube = new AliGTube("tube", "tube", 150,200,400);
    AliGCone* cone = new AliGCone("cone", "cone", 50,70,120,100,150);
    AliGSphere* sphere = new AliGSphere("sphere", "sphere", 75);
    AliGThorus* thorus = new AliGThorus("thorus", "thorus", 150,200,400);
    
    //  Set shapes attributes
    box1->SetLineColor(2);  // Red
    sphere->SetLineColor(4); // Blue
    tube->SetLineColor(3); // Bright Green
    cone->SetLineColor(5); // Yellow
    thorus->SetLineColor(2);  // Red

    // Draw inside Pad1
    Pad1->cd(); // Set current pad
    box1->Draw();

    // Draw inside Pad2
    Pad2->cd(); //Set current Pad
    sphere->Draw();

    // Draw inside Pad3
    Pad3->cd(); //Set current Pad
    tube->Draw();

    // Draw inside Pad4
    Pad4->cd(); //Set current Pad
    cone->Draw();

    c1->cd();
    c1->Update();  // Update Canvas
    */
}



