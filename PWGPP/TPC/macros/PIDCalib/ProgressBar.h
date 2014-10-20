#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

static Bool_t printedSomething = kFALSE;
static Int_t progressCounter = 0;
static const TString slash[4] = { "\\", "-", "/", "|" };
void progressbar(Int_t percent)
{
  // If something else has been printed, do not overwrite!
  if (printedSomething) {
    std::cout << std::endl;
    printedSomething = kFALSE;
  }
  std::cout << "\r"; // carriage return back to beginning of line
  for (Int_t i = 0; i < percent; i++)
    std::cout << ".";
    
  std::cout << " " << slash[progressCounter].Data() << " " << percent << " %" << std::flush; // print the bars and percentage
  progressCounter++; // increment to make the slash appear to rotate
  if(progressCounter == 4)
    progressCounter = 0; // reset slash animation
}


#endif
