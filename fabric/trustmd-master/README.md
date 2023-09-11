# TrustMD
Minifab run instructions

```
minifab up -e 7100 -o org1.example.com -n simple -l java -s couchdb -p '"init","a","200","b", "300"'
sudo chown -R ${USER} vars/chaincode/
cd vars/chaincode && git clone https://github.com/DJAcquila/trustmd && cd ../..
minifab ccup -v 1.0.0 -n trustmd -l java -p '"init","ch0","en0","1","0.2","mh0"'
```
