#!/bin/bash
HOST="localhost"   # IP of machine on which mysql database is located
PORT="5055"
DATABASE="storage"
USER="storage"
PASS="storage123"
TABLE="events"
STORAGE_PATH="/Users/Jerus/storedFiles"
MAX_SIZE="30000000"
MAX_OCCUPATION="80"
REMOVE_PERCENT="60"
EVENTS_IN_FILE="5"
EVENT_SERVER="pcald39fix"   # IP of machine running alieventserver
EVENT_SERVER_PORT="5024"
STORAGE_SERVER="localhost"      # IP of machine running alistorage
STORAGE_SERVER_PORT="5066"      # server thread communication port
STORAGE_CLIENT_PORT="5088"      # client thread communication port
XML_SERVER_PORT="5099"          # server of xml files
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