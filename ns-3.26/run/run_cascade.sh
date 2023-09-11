#!/bin/bash

#./waf --run "scratch/cascade --ns3::ConfigStore::Mode=Load --ns3::ConfigStore::Filename=input-defaults.txt"
./waf --command-template="%s --ns3::ConfigStore::Filename=input-defaults.txt --ns3::ConfigStore::Mode=Load --ns3::ConfigStore::FileFormat=RawText" --run scratch/cascade
