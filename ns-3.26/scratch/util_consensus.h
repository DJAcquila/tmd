#include <algorithm>
#include <sstream> 
#include <iostream>
#include <set>
#include <map>
#include <vector>


class VideoMetadata {
private:
	int valor_;
public:
	VideoMetadata(int valor): valor_(valor) {}
	VideoMetadata(){}

	void setVal(int valor) {
		valor_ = valor;

	}
	int getValor() {
		return valor_;
	}
};


class Transaction {
private:
	int receiverID_;
	Ipv4Address receiverIp_;
	int senderID_;
	Ipv4Address senderIp_;
	int contID_;
	double time_;
	double timeRecv_;
	VideoMetadata met_;
	
	std::string transactionHash_;

public:
	Transaction (int id, Ipv4Address ip, int sender_id, Ipv4Address sender_ip, int cont, double time, double timeRecv, VideoMetadata met): receiverID_(id), receiverIp_(ip), senderID_(sender_id), senderIp_(sender_ip), contID_(cont), time_(time), timeRecv_(timeRecv) {
		met_ = met;	
		generateTxHash();
	}
	Transaction(){}

	Ipv4Address getIP() {
		return this->receiverIp_;
	}

	Ipv4Address getSenderIP() {
		return this->senderIp_;
	}

	int getSenderID() {
		return this->senderID_;
	}

	int getContID() {
		return this->contID_;
	}
	int getRecvID() {
		return this->receiverID_;
	}
	double getTime() {
		return this->time_;
	}
	double getTimeRecv() {
		return this->timeRecv_;
	}
	VideoMetadata getMetadata() {
		return met_;
	}

	std::string generateTxHash() {
		std::string msg  = std::to_string(getRecvID()) 
							+ std::to_string(getSenderID()) 
							+ std::to_string(getMetadata().getValor())
							+ std::to_string(getContID())
							+ std::to_string(getTime())
							+ std::to_string(getTimeRecv());

		CryptoPP::SHA256 hash;
		byte digest[ CryptoPP::SHA256::DIGESTSIZE ];
     	hash.CalculateDigest( digest, reinterpret_cast<byte*>(&msg[0]), msg.length() );

     	CryptoPP::HexEncoder encoder;
     	std::string output;
     	encoder.Attach( new CryptoPP::StringSink( output ) );
     	encoder.Put( digest, sizeof(digest) );
        encoder.MessageEnd();

		this->transactionHash_ = output;
		return this->transactionHash_;
	}
	
	std::string getTxHash() {
		return transactionHash_;
	}

	bool validateTxHash(std::string hash){
		if(hash.compare(this->transactionHash_) == 0)
			return true;
		else
			return false;
	}

	void printTransacion() {
		std::cout << "(" << getIP() 
		<< ", " << getSenderIP()
		<< ", " << getContID()
		<< ", " << getTime()
		<< ", " << getTimeRecv()
		<< ", " << getRecvID()
		<< ", " << getSenderID()
		<< ", " << getMetadata().getValor()
		<< ") " << std::endl;
	}
};


class UpdateCluster{
private:
	std::vector<int> m_members;
	int m_clusterSize;
	int m_chID;
	Ipv4Address m_chAddress;
	bool m_isCH;
	std::vector<Ipv4Address> m_peersAddresses;
	std::map<Ipv4Address, Ptr<Socket>> m_peersSockets;
	std::map<Ipv4Address, Ptr<Socket>> m_clusterPeersSockets;
	Ptr<Socket> m_clusterSocket;

	bool _flag;
public:

	UpdateCluster(){}

	void set(std::vector<int> members, int clusterSize, int chID, Ipv4Address chAddress, bool isCH, std::vector<Ipv4Address> peersAddresses, 
			 std::map<Ipv4Address, Ptr<Socket>> peersSockets, 
			 std::map<Ipv4Address, Ptr<Socket>> clusterPeersSockets, Ptr<Socket> clusterSocket){

		m_members = members;
		m_clusterSize = clusterSize;
		m_chID = chID;
		m_chAddress = chAddress;
		m_isCH = isCH;
		m_peersAddresses = peersAddresses;
		m_peersSockets = peersSockets;
		m_clusterPeersSockets = clusterPeersSockets;
		m_clusterSocket = clusterSocket;
		_flag = true;
	}
	void get(std::vector<int>* members, int* clusterSize, int* chID, Ipv4Address* chAddress, bool* isCH, std::vector<Ipv4Address>* peersAddresses, 
			 std::map<Ipv4Address, Ptr<Socket>>* peersSockets, 
			 std::map<Ipv4Address, Ptr<Socket>>* clusterPeersSockets, Ptr<Socket>* clusterSocket){
		*members = m_members;
		*clusterSize = m_clusterSize;
		*chID = m_chID;
		*chAddress = m_chAddress;
		*isCH = m_isCH;
		*peersAddresses = m_peersAddresses;
		*peersSockets = m_peersSockets;
		*clusterPeersSockets = m_clusterPeersSockets;
		*clusterSocket = m_clusterSocket;
		_flag = false;
	}

	void setMembers(std::vector<int> members){
		m_members = members;
	}
	void setClusterSize(int clusterSize){
		m_clusterSize = clusterSize;
	}
	void setChID(int chID){
		m_chID = chID;
	}
	void setChAddress(Ipv4Address chAddress){
		m_chAddress = chAddress;
	}
	void setIsCH(bool isCH) {
		m_isCH = isCH;
	}
	void setPeersAdress(std::vector<Ipv4Address> peersAddresses){
		m_peersAddresses = peersAddresses;
	}

	void setPeersSockets(std::map<Ipv4Address, Ptr<Socket>> peersSockets){
		m_peersSockets = peersSockets;
	}
	void setClusterPeersSockets(std::map<Ipv4Address, Ptr<Socket>> clusterPeersSockets){
		m_clusterPeersSockets = clusterPeersSockets;
	}
	void setClusterSocket(Ptr<Socket> clusterSocket){
		m_clusterSocket = clusterSocket;
	}

	std::vector<int> getMembers(){
		return m_members;
	}
	int getClusterSize(){
		return m_clusterSize;
	}
	int getChID(){
		return m_chID;
	}
	Ipv4Address getChAddress(){
		return m_chAddress;
	}

	bool getIsCH(){
		return m_isCH;
	}
	std::vector<Ipv4Address> getPeersAddress(){
		return m_peersAddresses;
	}
	std::map<Ipv4Address, Ptr<Socket>> getPeersSockets(){
		return m_peersSockets;
	}
	std::map<Ipv4Address, Ptr<Socket>> getClusterPeersSockets(){
		return m_clusterPeersSockets;
	}
	Ptr<Socket> getClusterSocket(){
		return m_clusterSocket;
	}

	bool hasToChange(){
		return _flag;
	}

	void changed(){
		_flag = false;
	}

};