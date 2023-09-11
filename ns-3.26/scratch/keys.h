#include <vector>
#include <utility>
#include <stack>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdlib.h>

using std::map;
using std::pair;

class Key {
private:
	static bool instanceFlag;
	static Key *keyUtil;
	Key() {this->pool_add = -1;}

	map<int, std::string> ECDSA_publicKeyMap;
	map<int, std::string> ECC_publicKeyMap;
	map<int, std::string> ECDSA_privateKeyMap;
	map<int, std::string> ECC_privateKeyMap;

	map<int, map<int, std::string>> ECDSA_pool;
	map<int, map<int, std::string>> ECC_pool;

public:
	int pool_add;
	static Key* getInstance();
	// Public key
	std::string getECDSA_publicKeyFromMap(int nodeID);
	void putECDSA_publicKeyInMap(int nodeID, std::string key);

	std::string getECC_publicKeyFromMap(int nodeID);
	void putECC_publicKeyInMap(int nodeID, std::string key);

	// Private Key
	std::string getECDSA_privateKeyFromMap(int nodeID);
	void putECDSA_privateKeyInMap(int nodeID, std::string key);

	std::string getECC_privateKeyFromMap(int nodeID);
	void putECC_privateKeyInMap(int nodeID, std::string key);


	//pool
	void put_pool(std::string type, int nodeID, int srcNode, std::string key);
	std::string getKeyFromPool(std::string type, int nodeID, int srcNode);
	void printNodeECDSAPool(int nodeID);

	//erase
	void erase(void);

	~Key() {
		instanceFlag = false;
	}
};

bool Key::instanceFlag = false;
Key* Key::keyUtil = NULL;

Key* Key::getInstance() {
	if (!instanceFlag) {
		keyUtil = new Key();
		instanceFlag = true;
	}
	return keyUtil;
}
void Key::erase(void) {
	std::cout << "Keys - ";
	// for(map<int, std::string>::iterator it = ECDSA_publicKeyMap.begin(); it != ECDSA_publicKeyMap.end(); it++) {
	// 	ECDSA_publicKeyMap.erase(it);
	// }
	// for(map<int, std::string>::iterator it = ECC_publicKeyMap.begin(); it != ECC_publicKeyMap.end(); it++) {
	// 	ECC_publicKeyMap.erase(it);
	// }
	// for(map<int, std::string>::iterator it = ECDSA_privateKeyMap.begin(); it != ECDSA_privateKeyMap.end(); it++) {
	// 	ECDSA_privateKeyMap.erase(it);
	// }
	// for(map<int, std::string>::iterator it = ECC_privateKeyMap.begin(); it != ECC_privateKeyMap.end(); it++) {
	// 	ECC_privateKeyMap.erase(it);
	// }
	for(map<int, map<int, std::string>>::iterator it = ECDSA_pool.begin(); it != ECDSA_pool.end(); it++) {
		ECDSA_pool.erase(it);
	}
	for(map<int, map<int, std::string>>::iterator it = ECC_pool.begin(); it != ECC_pool.end(); it++) {
		ECC_pool.erase(it);
	}
	std::cout << "instance erased" << std::endl;
}
//ECDSA: Private
std::string Key::getECDSA_privateKeyFromMap(int nodeID) {
	map<int, std::string>::iterator it;
	it = ECDSA_privateKeyMap.find(nodeID);
	if (it != ECDSA_privateKeyMap.end()) {
		return it->second;
	} else {
		return std::string("-1");
	}
}

//ECDSA: Private
void Key::putECDSA_privateKeyInMap(int nodeID, std::string key) {
	ECDSA_privateKeyMap.insert(pair<int, std::string> (nodeID, key));
}

//ECDSA: Public
std::string Key::getECDSA_publicKeyFromMap(int nodeID) {
	map<int, std::string>::iterator it;
	it = ECDSA_publicKeyMap.find(nodeID);
	if (it != ECDSA_publicKeyMap.end()) {
		return it->second;
	} else {
		return std::string("-1");
	}
}
//ECDSA: Public
void Key::putECDSA_publicKeyInMap(int nodeID, std::string key) {
	ECDSA_publicKeyMap.insert(pair<int, std::string> (nodeID, key));
}
//ECC: Private
std::string Key::getECC_privateKeyFromMap(int nodeID) {
	map<int, std::string>::iterator it;
	it = ECC_privateKeyMap.find(nodeID);
	if (it != ECC_privateKeyMap.end()) {
		return it->second;
	} else {
		return std::string("-1");
	}
}
//ECC: Private
void Key::putECC_privateKeyInMap(int nodeID, std::string key) {
	ECC_privateKeyMap.insert(pair<int, std::string> (nodeID, key));
}

//ECC: Public
std::string Key::getECC_publicKeyFromMap(int nodeID) {
	map<int, std::string>::iterator it;
	it = ECC_publicKeyMap.find(nodeID);
	if (it != ECC_publicKeyMap.end()) {
		return it->second;
	} else {
		return std::string("-1");
	}
}
/*std::vector<int> Key::getPubKeyNodes() {
	map<int, std::string>::iterator it;
	std::vector<int> v;
	for (it = ECDSA_publicKeyMap.begin(); it != ECDSA_publicKeyMap.end(); it++) {
		v.push_back(it->first);
	}
	return v;
}*/
//ECC: Public
void Key::putECC_publicKeyInMap(int nodeID, std::string key) {
	ECC_publicKeyMap.insert(pair<int, std::string> (nodeID, key));
}

void Key::put_pool(std::string type, int nodeID, int srcNode, std::string key) {
	if (type.compare(std::string("ECDSA")) == 0 || type.compare(std::string("ecdsa")) == 0){
		map<int, map<int, std::string>>::iterator it;
		it = ECDSA_pool.find(nodeID);
		if (it != ECDSA_pool.end()){
			map<int, std::string>::iterator it1;
			it1 = it->second.find(srcNode);

			if (it1 != it->second.end()) {
				it->second[srcNode] = key;
			} else {
				it->second.insert(pair<int, std::string>(srcNode, key));
			}
		} else {
			map<int, std::string> tmp;
			tmp.insert(pair<int, std::string>(srcNode, key));
			ECDSA_pool.insert(pair<int, map<int, std::string>>(nodeID, tmp));
		}
	}
	else if (type.compare(std::string("ECC")) == 0 || type.compare(std::string("ecc")) == 0) {
		map<int, map<int, std::string>>::iterator it;
		it = ECC_pool.find(nodeID);
		if (it != ECC_pool.end()){
			map<int, std::string>::iterator it1;
			it1 = it->second.find(srcNode);

			if (it1 != it->second.end()) {
				it->second[srcNode] = key;
			} else {
				it->second.insert(pair<int, std::string>(srcNode, key));
			}
		} else {
			map<int, std::string> tmp;
			tmp.insert(pair<int, std::string>(srcNode, key));
			ECC_pool.insert(pair<int, map<int, std::string>>(nodeID, tmp));
		}
	}
}

std::string Key::getKeyFromPool(std::string type, int nodeID, int srcNode) {
	if (type.compare(std::string("ECDSA")) == 0 || type.compare(std::string("ecdsa")) == 0){
		map<int, map<int, std::string>>::iterator it;
		it = ECDSA_pool.find(nodeID);
		if (it != ECDSA_pool.end()){
			map<int, std::string>::iterator it1;
			it1 = it->second.find(srcNode);
			if (it1 != it->second.end()) {
				return it1->second;
			} else {
				return std::string("-1");
			}
		} else {
			return std::string("keyHit");
		}
	}
	else if (type.compare(std::string("ECC")) == 0 || type.compare(std::string("ecc")) == 0) {
		map<int, map<int, std::string>>::iterator it;
		it = ECC_pool.find(nodeID);
		if (it != ECC_pool.end()){
			map<int, std::string>::iterator it1;
			it1 = it->second.find(srcNode);
			if (it1 != it->second.end()) {
				return it1->second;
			} else {
				return std::string("-1");
			}
		} else {
			return std::string("keyHit");
		}
	}
	return std::string("-1");
}

void Key::printNodeECDSAPool(int nodeID) {
	map<int, map<int, std::string>>::iterator it;
	it = ECDSA_pool.find(nodeID);
	if (it != ECDSA_pool.end()){
		map<int, std::string>::iterator it1;
		std::cout << "  -> keys - Node ("  << nodeID <<  ") ECDSA pool: ";
		for (it1 = it->second.begin(); it1 != it->second.end(); it1++) {
			std::cout << it1->first << " ";
		}
		std::cout << std::endl;
		
	} else {
		std::cout << "  -> keys - Node ("  << nodeID <<  ") ECDSA pool: ";
		std::cout << "Vazia" << std::endl;
	}

	it = ECC_pool.find(nodeID);
	if (it != ECC_pool.end()){
		map<int, std::string>::iterator it1;
		std::cout << "  -> keys - Node ("  << nodeID <<  ") ECDH pool: ";
		for (it1 = it->second.begin(); it1 != it->second.end(); it1++) {
			std::cout << it1->first << " ";
		}
		std::cout << std::endl;
		
	} else {
		std::cout << "  -> keys - Node ("  << nodeID <<  ") ECDH pool: ";
		std::cout << "Vazia" << std::endl;
	}


}