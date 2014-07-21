#!/bin/bash
rm -fr /opt/reco/run197669/
rm -fr /opt/reco/log/run197669.log
rm -fr ~/storedFiles/
mysql -u storage -pstorage123 -e "drop database storage"
$ALICE_ROOT/STORAGE/setupStorageDatabase.sh
