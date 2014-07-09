#!/bin/bash
HOST="localhost"   # IP of machine on which mysql database is located
PORT="5055"
DATABASE="storage"
USER="storage"
PASS="storage123"
TABLE="events"
STORAGE_PATH="/local/home/edis/storedFiles"
MAX_SIZE="5000000"
MAX_OCCUPATION="80"
REMOVE_PERCENT="60"
EVENTS_IN_FILE="7"
EVENT_SERVER="localhost"   # IP of machine running alieventserver
STORAGE_SERVER="localhost"      # IP of machine running alistorage
mysql -u root -pdaq -e "create database if not exists $DATABASE;"
mysql -u root -pdaq -e "grant ALL PRIVILEGES on $DATABASE.* to '$USER'@'$HOST' identified by '$PASS';"
mysql -u root -pdaq -e "use $DATABASE;"
mysql -u root -pdaq -e "CREATE TABLE IF NOT EXISTS $DATABASE.$TABLE(\
run_number int(6) NOT NULL,\
event_number int(6) NOT NULL,\
system text(7) DEFAULT NULL,\
multiplicity int(5) DEFAULT NULL,\
permanent tinyint(1) DEFAULT NULL,\
file_path text(100) DEFAULT NULL,\
PRIMARY KEY(run_number,event_number));"

echo "-----------------------------"
echo "Databases successfuly created"
echo "-----------------------------"






















