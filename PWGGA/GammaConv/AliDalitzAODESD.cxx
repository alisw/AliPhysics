#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliDalitzData.h"
#include "AliDalitzAODESD.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliLog.h"
#include "AliESDInputHandler.h"
#include "AliExternalTrackParam.h"
#include <iostream>
using namespace std;



ClassImp( AliDalitzAODESD )
//-----------------------------------------------------------------------------------------------
         
     Double_t AliDalitzAODESD::GetPtG(){
        if (fIsESD==kTRUE) return fESDtrack->GetConstrainedParam()->Pt();
        else return fAODtrack->Pt();
     }
     Double_t AliDalitzAODESD::GetPxG(){
        if (fIsESD==kTRUE) return fESDtrack->GetConstrainedParam()->Px();
        else return fAODtrack->Px();
     }
     Double_t AliDalitzAODESD::GetPyG(){
        if (fIsESD==kTRUE) return fESDtrack->GetConstrainedParam()->Py();
        else return fAODtrack->Py();
     }     
     Double_t AliDalitzAODESD::GetPzG(){
        if (fIsESD==kTRUE) return fESDtrack->GetConstrainedParam()->Pz();
        else return fAODtrack->Pz();
     }
    Double_t AliDalitzAODESD::GetPhiG(){
        if (fIsESD==kTRUE) return fESDtrack->GetConstrainedParam()->Phi();
        else return fAODtrack->Phi();
     }
    const AliExternalTrackParam* AliDalitzAODESD::GetParamG(){
        if (fIsESD==kTRUE) return fESDtrack->GetConstrainedParam();
        else{ //AliExternalTrackParam* aodParam;
            const AliExternalTrackParam* aodParam1=0;
           //AliExternalTrackParam etp;
         //   etp.CopyFromVTrack(fAODtrack);

            std::cout<<"definio aodParam"<<std::endl;
           // std::cout<<fAODtrack->Pt()<<endl;
           // AliExternalTrackParam((AliVTrack*)fAODtrack);
           // aodParam->CopyFromVTrack(fAODtrack);
            std::cout<<"Copio"<<std::endl;
           // const AliExternalTrackParam* aodParam1 =static_cast< const AliExternalTrackParam*>(*etp);
            //const AliExternalTrackParam* aodParam1=etp;
            aodParam1 = fAODtrack->GetOuterParam();
            std::cout<<"Lista para retornar"<<std::endl;
            //return aodParam1;
            return aodParam1;
        }
     }
    Int_t AliDalitzAODESD::GetLabelG(){
        if (fIsESD==kTRUE) return fESDtrack->GetLabel();
        else return fAODtrack->GetLabel();
     }
     Bool_t AliDalitzAODESD::GetPxPyPzG(Double_t* p) const{
        if (fIsESD==kTRUE) return fESDtrack->GetConstrainedPxPyPz(p);//NOTE debe ir GetConstrainedParam
        //else return fAODtrack->GetPxPyPz(p);
        else return fAODtrack->PxPyPzAtDCA(p);
     };
     void AliDalitzAODESD::GetImpactParametersG(Float_t* p,Float_t* cov) const{
        p[0]=b[0];
        p[1]=b[1];
        cov[0]=bCov[0];
        cov[1]=bCov[1];
        cov[2]=bCov[2];
     };
     Double_t AliDalitzAODESD::GetSignG(){
        if (fIsESD==kTRUE) return fESDtrack->GetSign();
        else return fAODtrack->Charge();
     };
     Double_t AliDalitzAODESD::GetEtaG(){
        if (fIsESD==kTRUE) return fESDtrack->Eta();
        else return fAODtrack->Eta();
     };
    void AliDalitzAODESD::ComputeImpactParameter(){
          //Selection of Reconstructed electrons
            Float_t bF[2],bCovF[3];
            //fESDtrack->GetImpactParameters(bF,bCovF);//GetImpactParameters(b,cov)
            if(fIsESD) fESDtrack->GetImpactParameters(bF,bCovF);
            else {
                bF[0]=fAODtrack->DCA();
                bF[1]=fAODtrack->ZAtDCA();
                fAODtrack->GetCovarianceXYZPxPyPz((Double_t*) bCovF);
            }
            b[0]=bF[0];//DCAXY
            b[1]=bF[1];//DCAZ
            bCov[0]=bCovF[0];
            bCov[1]=bCovF[1];
            bCov[2]=bCovF[2];
            if (bCov[0]<=0 || bCov[2]<=0) {
     //       AliDebug(1, "Estimated b resolution lower or equal zero!");
            bCov[0]=0; bCov[2]=0;
            }
    }
    void AliDalitzAODESD::ComputeImpactParameter(const AliVVertex* vx,Double_t bmag){
   AliExternalTrackParam etp; etp.CopyFromVTrack(fAODtrack);
   //AliAODVertex *vtxAOD = aod->GetPrimaryVertex();
   // activate the following two lines if you want to check that the event primary vertex is that reconstructed with tracks
   //  TString title=vtxAOD->GetTitle();
   //  if(!title.Contains("VertexerTracks")) { you could decide what to do in case the primary vertex is not reconstructed with tracks }
   //Double_t b=aod->GetMagneticField();
   if(!(etp.PropagateToDCA(vx,bmag,3.0,b,bCov))){
            b[0]=-999;//DCAXY
            b[1]=-999;//DCAZ
            bCov[0]=-999;
            bCov[1]=-999;
            bCov[2]=-999;
   }
   // --> see above what dz and covdz represen
            if (bCov[0]<=0 || bCov[2]<=0) {
        //    AliDebug(1, "Estimated b resolution lower or equal zero!");
            bCov[0]=0; bCov[2]=0;
            }
    }
    Double_t AliDalitzAODESD::GetTPCNclsFG(){
        if (fIsESD==kTRUE) return fESDtrack->GetTPCNclsF();
        else return fAODtrack->GetTPCNclsF();
    }
     Double_t AliDalitzAODESD::GetNclsG(){
        if (fIsESD==kTRUE) return fESDtrack->GetNcls(1);
        else return fAODtrack->GetTPCNcls();
     }
     Double_t AliDalitzAODESD::GetTPCCrossedRowsG(){
        if (fIsESD==kTRUE) return fESDtrack->GetTPCCrossedRows();
        else return fAODtrack->GetTPCNCrossedRows();   
     }
    Int_t AliDalitzAODESD::GetITSclsG(){
        if (fIsESD==kTRUE) return fESDtrack->GetNcls(0);
        else return fAODtrack->GetITSNcls();   
     }
    Bool_t AliDalitzAODESD::HasPointOnITSLayerG(Int_t i){
        if (fIsESD==kTRUE) return fESDtrack->HasPointOnITSLayer(i);
        else return fAODtrack->HasPointOnITSLayer(i);   
     }
    Double_t AliDalitzAODESD::GetDCAxy(){
        return b[0];
    }
    Double_t AliDalitzAODESD::GetDCAz(){
        return b[1];
    }
    
    
