#include "AliNanoSkimmingV0s.h"
#include <TObject.h>
#include <AliVEvent.h>

AliNanoSkimmingV0s::AliNanoSkimmingV0s()
    : AliAnalysisCuts{}
{
}

bool AliNanoSkimmingV0s::IsSelected(TObject *obj)
{
    fSelected = false;
    AliVEvent *ev = dynamic_cast<AliVEvent *>(obj);
    if (!ev)
    {
        return fSelected;
    }

    fSelected = ev->GetNumberOfV0s() > 0;

    return fSelected;
}


bool AliNanoSkimmingV0s::IsSelected(TList *)
{
    Fatal("AliNanoSkimmingV0s::IsSelected", "Method not implemented for lists");
    return false;
}
