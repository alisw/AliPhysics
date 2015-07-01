#include "AliZMQManager.h"
#include "AliStorageTypes.h"
#include <iostream>
#include <sstream>
#include <ostream>
#include <AliTrackPointArray.h>
#include <TXMLEngine.h>				

using namespace std;

TXMLEngine* xml = new TXMLEngine;

stringstream& getXml(AliESDEvent *event)
{
	// Create main node of document tree
	XMLNodePointer_t mainnode = xml->NewChild(0, 0, "main");
	
	cout<<"tracks:"<<event->GetNumberOfTracks()<<endl;
	XMLNodePointer_t tracks[event->GetNumberOfTracks()];
	
	for(int i=0;i<event->GetNumberOfTracks();i++)
	{
		AliESDtrack *track = event->GetTrack(i);
		tracks[i] = xml->NewChild(mainnode, 0, Form("track%d",i));
		const AliTrackPointArray *array = track->GetTrackPointArray();
		if(array)
		{
			const float *x = array->GetX();
			const float *y = array->GetY();
			const float *z = array->GetZ();
			int n = array->GetNPoints();

			for(int j=0;j<n;j++)
			{
				cout<<"3"<<endl;
				xml->NewChild(tracks[i], 0,Form("point%d",j),Form("%f\t%f\t%f\n",x[j],y[j],z[j]));
			}	
		}
		else cout<<"no array"<<endl;
	}
   
	stringstream streamXml;
	xml->SavePrimitive(streamXml);
	delete xml;
	return streamXml;
}

int main()
{
	AliZMQManager *manager = AliZMQManager::GetInstance();
	AliESDEvent *event;

	while(1)
	{
        manager->Get(event,EVENTS_SERVER_SUB);
		cout<<"sending xml"<<endl;
		manager->SendAsXml(event,XML_PUB);
		cout<<"xml sent"<<endl;
	}
	return 0;
}


