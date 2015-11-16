//#include "AliColladaBuffer.h"
//#include "TGeoSphere.h"
//#include "TGeoBBox.h"
//#include "TGeoTube.h"

void colladaBufferExample()
{
    // prepare collada buffer
    AliColladaBuffer *colladaBuffer = new AliColladaBuffer();
    
    int iter=0;
    
    // transformation of shapes
    double translation[3] = {0.,0.,0.};
    double rotation[3] = {0.,0.,0.};
    double scale[3] = {1.,1.,1.};
    
    colladaBuffer->CreateNode("myShapes"); // create one main node
    colladaBuffer->CreateNode("rounded","myShapes"); // add sub-node for rounded shapes
    
    
    // Create sphere and add to Collada buffer
    TGeoSphere *sphere = new TGeoSphere(0.,10.);
    sphere->SetNumberOfDivisions(100);
    
    colladaBuffer->AddShape("sphere",sphere->MakeBuffer3D(),translation,rotation,scale,"rounded");

    // Create a box and add it to Collada buffer
    TGeoBBox *box = new TGeoBBox(5.,10.,10.);
    
    colladaBuffer->AddShape("box",box->MakeBuffer3D(),translation,rotation,scale,"myShapes");

    // Create a tube and add it to Collada buffer
    TGeoTube *tube = new TGeoTube(2.,2.5,30);
    translation[1] = 30;
    
    colladaBuffer->AddShape("tube",tube->MakeBuffer3D(),translation,rotation,scale,"rounded");
    
    
    // save collada buffer to file
    colladaBuffer->SaveAs("exampleBoxes.dae");
    
    // clean up
    delete colladaBuffer;
    return;
    
}


