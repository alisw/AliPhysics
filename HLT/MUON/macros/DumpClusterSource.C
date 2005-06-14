////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* Dumps the contents of the ClusterSource 'cs' to the screen.
 */
void DumpClusterSource(AliMUONHLT::ClusterSource* cs)
{
	cout << "File  : " << cs->FileName() << endl;
	cout << "Folder: " << cs->FolderName() << endl;
	if (cs->GetFirstEvent())
	{
		do
		{
			cout << "================= Event: " << cs->CurrentEvent() << " =================" << endl;
			if ( ! cs->GetFirstBlock() )
			{
				cout << "No blocks found." << endl;
				continue;
			};
			do
			{
				cout << "Block for chamber: " << cs->Chamber() << endl;
				if ( cs->GetFirstCluster() == NULL )
				{
					cout << "\tNo cluster points found." << endl;
					continue;
				};
				do
				{
					Float_t x, y;
					cs->FetchCluster(x, y);
					cout << "\tx = " << x << ", y = " << y << endl;
				} while (cs->GetNextCluster());
			} while (cs->GetNextBlock());
		} while (cs->GetNextEvent());
	}
	else
	{
		cout << "No events found." << endl;
	};
};

