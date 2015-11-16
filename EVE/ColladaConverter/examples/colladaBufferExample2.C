//#include "AliColladaBuffer.h"

void colladaBufferExample2()
{
    // prepare collada buffer
    AliColladaBuffer *colladaBuffer = new AliColladaBuffer();
    
    int iter=0;
    
    // transformation of shapes
    double translation[3] = {0.,0.,0.};
    double rotation[3] = {0.,0.,0.};
    double scale[3] = {1.,1.,1.};
    
    // parameters of the lambert material
    TColor *emission = new TColor();
    TColor *ambient = new TColor();
    TColor *diffuse = new TColor();
    
    emission->SetAlpha(1.0);
    ambient->SetRGB(0.0,0.0,0.0);
    ambient->SetAlpha(1.0);
    diffuse->SetAlpha(1.0);
    
    // create some shapes in a loop
    for(int x=-5;x<5;x++)
    {
        for(int y=-5;y<5;y++)
        {
            for(int z=-5;z<5;z++)
            {
                translation[0] = x;
                translation[1] = y;
                translation[2] = z;
                
                emission->SetRGB((x+5)/20.,(y+5)/20.,(z+5)/20.);
                diffuse->SetRGB((x+5)/10.,(y+5)/10.,(z+5)/10.);
                
                // add new lambert material
                colladaBuffer->AddNewLambertMaterial(Form("lambert%d",iter),emission,ambient,diffuse,0.5);
             
                // add shape to the buffer
                colladaBuffer->AddBox(Form("cube%d",iter),0.5,0.5,0.5,translation,rotation,scale,"",Form("lambert%d",iter));
                iter++;
            }
        }
    }
    
    // save collada buffer to file
    colladaBuffer->SaveAs("exampleBoxes2.dae");
    
    // clean up
    delete colladaBuffer;
    return;
    
}


