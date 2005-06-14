////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* Dumps the contents of the TriggerSource 'ts' to the screen.
 */
void DumpTriggerSource(AliMUONHLT::TriggerSource* ts)
{
	cout << "File  : " << ts->FileName() << endl;
	cout << "Folder: " << ts->FolderName() << endl;
	if (ts->GetFirstEvent())
	{
		do
		{
			cout << "================= Event: " << ts->CurrentEvent() << " =================" << endl;
			Int_t blocknum = 0;
			if ( ! ts->GetFirstBlock() )
			{
				cout << "No blocks found." << endl;
				continue;
			};
			do
			{
				cout << "Block: " << blocknum++ << endl;
				if ( ts->GetFirstTrigger() == NULL )
				{
					cout << "\tNo trigger records found." << endl;
					continue;
				};
				do
				{
					AliMUONHLT::TriggerRecord* data = ts->GetTrigger();

					cout << "\tTrigger number = " << data->TriggerNumber()
						<< ", Sign = " << data->ParticleSign()
						<< ", Pt = " << data->Pt()
						<< ", ST1 x = " << data->Station1Point().fX
						<< ", ST1 y = " << data->Station1Point().fY
						<< ", ST2 x = " << data->Station2Point().fX
						<< ", ST2 y = " << data->Station2Point().fY
						<< endl;
				} while (ts->GetNextTrigger());
			} while (ts->GetNextBlock());
		} while (ts->GetNextEvent());
	}
	else
	{
		cout << "No events found." << endl;
	};
};

