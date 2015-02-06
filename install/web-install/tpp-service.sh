#!/bin/bash
while 1 
do
    cd /var/www/erg.biophys.msu.ru/htdocs/tpp/
    ./cron-starter.py > ./111
    sleep 5
done

