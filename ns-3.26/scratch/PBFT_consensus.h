#include "ns3/application.h"
#include "ns3/event-id.h"
#include "ns3/ptr.h"
#include "ns3/traced-callback.h"
#include "ns3/address.h"
#include "ns3/boolean.h"
#include "ns3/address.h"
#include "ns3/address-utils.h"
#include "ns3/log.h"
#include "ns3/inet-socket-address.h"
#include "ns3/inet6-socket-address.h"
#include "ns3/node.h"
#include "ns3/socket.h"
#include "ns3/udp-socket.h"
#include "ns3/simulator.h"
#include "ns3/socket-factory.h"
#include "ns3/packet.h"
#include "ns3/trace-source-accessor.h"
#include "ns3/udp-socket-factory.h"
#include "ns3/tcp-socket-factory.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"
#include "ns3/buffer.h"



#include <algorithm>
#include <sstream> 
#include <iostream>
#include <set>

#include "../rapidjson/document.h"
#include "../rapidjson/writer.h"
#include "../rapidjson/stringbuffer.h"

#include "crypto++/cryptlib.h"
#include "crypto++/sha.h"
#include "crypto++/hex.h"

#include "block.h"



using namespace ns3;

/*int main () {
   string myString    = string("abc");
   string myStringSHA = SHA256(myString);
   std::cout << "SHA(abc): " << myStringSHA << std::endl;
}
*/

BasicConfig* basicConfigPBFTInstance = BasicConfig::getInstance();
uint16_t nNodes = basicConfigPBFTInstance->getNumberOfNodes();

struct comp
{
	template<typename T>
	bool operator()(const T& l, const T& r) const
	{
		if (l.first == r.first)
			return l.second > r.second;

		return l.first < r.first;
	}
};

class TransactionPool {
private:
    //static TransactionPool *TransactionPool;
    std::queue <Transaction> pool;
public:
    TransactionPool(){}

    void addTransacion(int id, Ipv4Address ip,  int sender_id, Ipv4Address sender_ip, int cont, double time, double timeRecv, VideoMetadata met) {
        pool.push(Transaction(id, ip, sender_id, sender_ip, cont, time, timeRecv, met));
    } 
    Transaction getNextTransaction() {
    	return pool.front();
    } 
    void resetPool(){

    	std::queue<Transaction> empty;
    	std::swap(pool, empty);
    } 

    void popTransaction() {
    	if (!pool.empty()){
    		pool.pop();
    	}
    }
    bool empty() {
    	return pool.empty();
    } 
    void printQueue() {
    	auto tmp = pool;
    	int tid = 0;
    	while (!tmp.empty())
		{
			std::cout << "T" << tid << ": ";
			Transaction p = tmp.front();
			p.printTransacion();
			tmp.pop();
			tid++;
		}
    }
};

class PingSet {
private:
	std::set<Ipv4Address> m_pingSet;
public:
	PingSet(){}
	void addNode(Ipv4Address node){
		m_pingSet.insert(node);
	}
	void removeNode(Ipv4Address node) {
		m_pingSet.erase(node);
	}
	bool isEmpty() {
		if (m_pingSet.empty()) {
			return true;
		} else {
			return false;
		}
	}
	std::set<Ipv4Address> get(){
		return m_pingSet;
	}

	int count(){
		return (int) m_pingSet.size();
	}
	void reset(){
		m_pingSet.clear();
	}

};

class PingSession {
private:
	std::vector<PingSet> m_pingSet;
public:
	PingSession(){
		for (int i = 0; i < nNodes; i++){
			m_pingSet.push_back(PingSet());
		}
	}

	bool existInSet(int n, Ipv4Address addr){
		for (auto const &p: this->getSet(n)){
			if (p == addr){
				return true;
			}
		}
		return false;
	}

	void addPingSet(int n, Ipv4Address addr) {
		m_pingSet[n].addNode(addr);
	}

	void resetSession(int n) {
		m_pingSet[n].reset();
	}

	std::set<Ipv4Address> getSet(int n) {
		return m_pingSet[n].get();
	}

	void eraseSet(int n, Ipv4Address addr) {
		if (!m_pingSet[n].isEmpty()) {
			m_pingSet[n].removeNode(addr);
		}
	}

	void reset(Ipv4Address addr) {
		for (int i = 0 ; i < nNodes; i++){
			if (!m_pingSet[i].isEmpty()) {
				
				m_pingSet[i].removeNode(addr);
				
			}
		}
	}

	void reset(int n){
		m_pingSet[n].reset();
	}

	bool isEmpty(int n) {
		return m_pingSet[n].isEmpty();
	}
};

PingSession m_pingSet;
PingSession m_pingValidationSet;

class ConclusionSession {
private:
	std::vector<PingSet> m_set;
public:
	ConclusionSession(){
		for (int i = 0; i < nNodes; i++){
			m_set.push_back(PingSet());
		}
	}

	bool existInSet(int n, Ipv4Address addr){
		for (auto const &p: this->getSet(n)){
			if (p == addr){
				return true;
			}
		}
		return false;
	}

	void add(int n, Ipv4Address addr) {
		m_set[n].addNode(addr);
	}

	void resetSession(int n) {
		m_set[n].reset();
	}

	std::set<Ipv4Address> getSet(int n) {
		return m_set[n].get();
	}

	void eraseSet(int n, Ipv4Address addr) {
		if (!m_set[n].isEmpty()) {
			m_set[n].removeNode(addr);
		}
	}

	void reset(Ipv4Address addr) {
		for (int i = 0 ; i < nNodes; i++){
			if (!m_set[i].isEmpty()) {
				
				m_set[i].removeNode(addr);
				
			}
		}
	}

	bool checkConclusion(int n, int number){
		if(m_set[n].count() >= number){
			return true;
		}
		return false;
	}

	void reset(int n){
		m_set[n].reset();
	}

	bool isEmpty(int n) {
		return m_set[n].isEmpty();
	}
};

ConclusionSession m_conclusion;

class ValidationSession {
private:
	std::string sessionTxHash;
	Transaction transaction;
	std::set<std::pair<int, int>, comp>* votingPool;
	int* nodesDecision;
	std::set<Ipv4Address> votingNodes;
	std::set<Ipv4Address> aux_votingNodes;
	std::map<Ipv4Address, int> nodesDecisions; // para os clusterheads
	int decision;
	bool allDecisions = false;

	int clusterSize_penalty = 0;
public:
	bool m_voteComplete = false;
	ValidationSession(){
		votingPool = new std::set<std::pair<int, int>, comp>[nNodes];
		for (int i = 0; i < nNodes; i++) {
			std::set<std::pair<int, int>, comp> set;
			votingPool[i] = set;
		}

		nodesDecision = new int[nNodes];
		for (int i = 0; i < nNodes; i++) {
			nodesDecision[i] = -1;
		}
		clusterSize_penalty = 0;
	}

	void addTransacion(Transaction t){
		this->transaction = t;
	}
	Transaction getTransaction(){
		return this->transaction;
	}

	void newVotingNode(Ipv4Address addr){
		this->votingNodes.insert(addr);
		this->aux_votingNodes.insert(addr);
	}
	std::set<Ipv4Address> getVotingNodes() {
		return this->votingNodes;
	}
	void printAux() {
		std::cout << "aux_votingNodes: ";
		for (auto const &p: this->aux_votingNodes) {
			std::cout << p << " ";
		}
		std::cout << std::endl;
	}

	void newNodeDecision(Ipv4Address addr, int decision){
		std::map<Ipv4Address, int>::iterator it = this->nodesDecisions.find(addr);
		if (it == this->nodesDecisions.end()){
			this->nodesDecisions.insert(std::pair<Ipv4Address, int>(addr, decision));
			// Já recebeu a decisão do nó
			if (!this->aux_votingNodes.empty())
				this->aux_votingNodes.erase(addr);
		} 
	}
	bool decisionAlreadyMade(Ipv4Address addr){
		std::map<Ipv4Address, int>::iterator it = this->nodesDecisions.find(addr);
		if (it != this->nodesDecisions.end()){
			return true;
		}
		return false;
	}

	bool allDecisionsReceived(){
		if (!allDecisions){
			if(this->aux_votingNodes.empty()){
				return true;
				allDecisions = true;
			}
			return false;

		} else{
			return false;
		}

	}

	void setHash(std::string hash) {
		this->sessionTxHash = hash;
	}
	std::string getSessionHash () {
		return this->sessionTxHash;
	}

	void addVotingPool(int srcNode, int nodeId, int vote) {
		this->votingPool[srcNode].insert(std::pair<int, int>(nodeId, vote));
	}
	std::pair<bool, int> canDecide(int srcNode, int clusterTotal) {
		int invalid = 0;
		int valid = 0;

		for (auto const &p: this->votingPool[srcNode]) 
		{
			if (p.second == 1) {
				valid++;
			} else {
				invalid++;
			}
		}
		if (valid/clusterTotal >= 0.8){
			nodesDecision[srcNode] = 1;
			return std::pair<bool, int> (true, 1);
		}	
		if (invalid/clusterTotal >= 0.8){
			nodesDecision[srcNode] = 0;
			return std::pair<bool, int> (true, 0);
		}
		return std::pair<bool, int> (false, -1);
	}

	bool validateTx(int srcNode) {
		int invalid = 0;
		int valid = 0;
		
		for (auto const &p: this->votingPool[srcNode]) 
		{
			if (p.second == 1){
				valid++;
			} else {
				invalid++;
			}
		}
		
		if (valid > invalid){
			std::cout << "Valid: " << valid << "\nInvalid: " << invalid << std::endl;
			m_voteComplete = true;
			nodesDecision[srcNode] = 1;
			return true;
		}
		std::cout << "Valid: " << valid << "\nInvalid: " << invalid << std::endl;
		m_voteComplete = true;
		nodesDecision[srcNode] = 0;
		return false;
		
	}

	int getNodeDecision(int srcNode){
		return nodesDecision[srcNode];
	}
	void desclassifyNode(){
		clusterSize_penalty++;
	}
	
	bool verifyConclusion(int srcNode, std::vector<Ipv4Address> peersAddresses, Ipv4InterfaceContainer ueWiFiIface, int ch) {
		//bool final = true;
		/*for(std::vector<Ipv4Address>::const_iterator i = peersAddresses.begin(); i != peersAddresses.end(); ++i)
		{
			bool check = false;
			for (auto const &p: this->votingPool[srcNode]) 
			{
				if (ueWiFiIface.GetAddress(p.first) == *i || ueWiFiIface.GetAddress(ch) == *i) {
					check = true;
					break;
				}
			}
			if (!check) {
				//std::cout << "Consenso não concluído" << std::endl;
				return false;
			} 
		}

		validateTx(srcNode);
		return true;*/
		
		
		int valid = 0, invalid = 0;
		for (auto const &p: this->votingPool[srcNode]) 
		{
			if (p.second == 1){
				valid++;
			} else {
				invalid++;
			}
		}

		int size = peersAddresses.size() - clusterSize_penalty;
		
		if (valid >= size){
			std::cout << "Valid: " << valid << "\nInvalid: " << invalid << std::endl;
			nodesDecision[srcNode] = 1;
			m_voteComplete = true;
			return true;
		}
		else if ( invalid >= size ) {
			std::cout << "Valid: " << valid << "\nInvalid: " << invalid << std::endl;
			nodesDecision[srcNode] = 0;
			m_voteComplete = true;
			return true;
		}
		//std::cout << "Consenso concluído" << std::endl;
		return false;
		
	}	

	bool voteCompleted(){
		return m_voteComplete;
	}

	void printVotePool(int srcNode) {
		std::cout << "\t-> VotationPool: ";
		for (auto const &p: this->votingPool[srcNode]) 
		{
			std::cout << "(" << p.first << ", " << p.second << ")";
		}
		std::cout << std::endl;
		
	}
	bool validateTxHash(std::string recv_hash) {
		if (recv_hash == this->sessionTxHash)
			return true;
		else
			return false;
	}

	void resetValidationSession(int srcNode) {
		this->votingPool[srcNode].clear();

		//printVotePool(srcNode);
		//this->sessionTxHash.clear();
		this->votingNodes.clear();
		this->aux_votingNodes.clear();
		this->nodesDecisions.clear();
		this->clusterSize_penalty = 0;
		for (int i = 0; i < nNodes; i++) {
			nodesDecision[i] = -1;
		}
		allDecisions = false;
		m_voteComplete = false;
	}

	void setDecision(int d) {
		decision = d;
	}
	int getMineDecision() {
		return decision;
	}
};

class BlockConsenus {
private:
	Block m_block;
	std::map<std::string, Block> m_blockHashSet;
	std::map<std::string, int> m_hashCategories;

	int m_nCluster;
	int m_nblocksReceived;

public:
	BlockConsenus(int nCluster){
		m_nblocksReceived = 0;
		m_nCluster = nCluster;
	}
	BlockConsenus(){
		m_nblocksReceived = 0;
		m_nCluster = 0;
	}
	void setClusterSize(int size){
		m_nCluster = size;
	}
	bool existInSet(std::string hash){
		std::map<std::string, Block>::iterator it;
		for (it = m_blockHashSet.begin(); it != m_blockHashSet.end(); it++){
			if ((it->first).compare(hash) == 0) {
				return true;
			}
		}
		return false;
	}


	void desclassifyNode(){
		m_nCluster--;
	}

	Block getBLockByHash(std::string hash){
		std::map<std::string, Block>::iterator it;
		for (it = m_blockHashSet.begin(); it != m_blockHashSet.end(); it++){
			if ((it->first).compare(hash) == 0) {
				return it->second;
			}
		}
		Transaction t;
		return Block(t, -1, -1, -1, -1, std::string("-1"), -1, -1, -1, Ipv4Address("0.0.0.0"), -1); 
	}

	void update(std::string hash){
		if (existInSet(hash)){
			m_hashCategories[hash] += 1;
		} else {
			m_hashCategories[hash] = 1;
		}
	}

	void addNewBlock(std::string hash, Block block){
		m_nblocksReceived++;
		update(hash);
		m_blockHashSet[hash] = block;
	}
	void print() {
		std::cout << "Faltam: " << m_nCluster - m_nblocksReceived << " blocos" << std::endl;
	}
	bool verifyConclusion(Block* returned){
		if (m_nblocksReceived >= 0.8*(m_nCluster - 1)){
			std::map<std::string, int>::iterator it;
			int greater = 0;
			std::string greatHash;
			for (it = m_hashCategories.begin(); it != m_hashCategories.end(); it++){
				if (it->second > greater) {
					greater = it->second;
					greatHash = it->first;
				}
			}
			*returned = getBLockByHash(greatHash);
			return true;
		}
		return false;
	}
	//void addBlockReceived
	void reset(){
		m_blockHashSet.clear();

		std::cout << "\t-> m_blockHashSet: ";
		std::map<std::string, Block>::iterator it;
		for (it = m_blockHashSet.begin(); it != m_blockHashSet.end(); it++){
			std::cout << "(" << it->first << ", " << (it->second).GetContentID() << "), ";
		}
		std::cout << std::endl;

		m_hashCategories.clear();

		std::cout << "\t-> m_hashCategories: ";
		std::map<std::string, int>::iterator it1;
		for (it1 = m_hashCategories.begin(); it1 != m_hashCategories.end(); it1++){
			std::cout << "(" << it1->first << ", " << it1->second << "), ";
		}
		std::cout << std::endl;

		m_nblocksReceived = 0;

	}
};
//std::set<Ipv4Address> 					m_pingValidationSet;

class TrustVotes{
private:
	std::set<std::pair<int, int>>* trustVotes_pool;

public:
	TrustVotes(){
		trustVotes_pool = new std::set<std::pair<int, int>>[nNodes];
		for (int i = 0; i < nNodes; i++) {
			std::set<std::pair<int, int>> set;
			trustVotes_pool[i] = set;
		}
	}

	void addVotingPool(int srcNode, int nodeId, int vote) {
		this->trustVotes_pool[srcNode].insert(std::pair<int, int>(nodeId, vote));
	}
	
	int getTrustOpinion(int srcNode, int nodeId){
		for (auto const &p: this->trustVotes_pool[srcNode]) 
		{
			if (p.first == nodeId){
				return p.second;
			}
		}
	}

	std::set<std::pair<int, int>> getVotePool(int srcNode) {
		return this->trustVotes_pool[srcNode];
	}

	void print(int srcNode){
		std::cout << "Pool node " << srcNode << ": ";
		for (auto const &p: this->trustVotes_pool[srcNode]) 
		{
			std::cout << "// (" << p.first << "): " << p.second << " //";
		}
		std::cout << std::endl;
	}

	void reset(int srcNode){
		this->trustVotes_pool[srcNode].clear();
	}
};

TrustVotes 									m_trust_pool;

class ConsensusApp : public Application {
protected:
	Ptr<Socket>								m_socket;
	Ptr<Socket>								m_nodeClusterSocket;
	Address 								m_local;
	Address 								m_localSendCluster;
	bool 									m_consensus_open;
	double 									m_consensusStartTime;
	double 									m_consensusEndTime;
	uint32_t 								m_packetSize;
	//uint32_t 								m_nPackets;
	//DataRate 								m_dataRate;
	double 									m_simTime;
	double 									m_lastTime;
	EventId 								m_nextSendEvent;
	EventId 								m_nextPingEvent;
	EventId 								m_nextValidationPingEvent;
	EventId 								m_nextCheckValidationStartEvent;
	EventId 								m_validateNextEvent;
	EventId 								m_nextSendPingValidationEvent;
	bool 									m_running;
	bool 									m_isCH;
	std::map<Address, std::string>			m_bufferedData; 
	std::queue<double>						m_timer;
	std::queue<double>						m_clusterTimeSchedule;
	std::queue<double>						m_updateTimeSchedule;
	TransactionPool 						m_transactionPool;
	std::queue<double>						m_queueSendTransactions;
	std::list<Ptr<Socket>> 					m_socketList;
	//int 									m_numberOfPeers;
	int 									m_CHNodeId;
	Ipv4Address 							m_chAddress;
	std::vector<Ipv4Address> 				m_peersAddresses;
	std::map<Ipv4Address, Ptr<Socket>>	 	m_peersSockets;
	std::map<Ipv4Address, Ptr<Socket>>	 	m_clusterPeersSockets;
	Ptr<Socket> 							m_clusterSocket;
	nodeStatistics 							*m_nodeStats;
	Ipv4InterfaceContainer 					m_ueWiFiIface;
	NodeContainer 							m_ueNodes;
	uint16_t 								m_SendPort;
	uint16_t 								m_sendoToClusterPort;
	ValidationSession 						m_validationSession;
	int 									m_response;
	int 									m_clusterSize;
	std::vector<int> 						m_members;
	//std::set<Ipv4Address> 					m_pingSet;
	rapidjson::StringBuffer 				m_buffer;
	rapidjson::StringBuffer 				m_bufferValidation;
	Blockchain 								m_blockchain;
	std::set<int> 							m_decisionSet;

	const int 								m_trasactionSizeBytes; // 32 bytes
	const int 								m_blockHeaderSizeBytes; // 68 bytes

	BlockConsenus 							m_blockConsenus;
	bool 									m_blockBeingFormed;

	rapidjson::Document 					m_chMessage;
	bool 									m_alreadyScheduled;

	std::set<int>	  						m_waiting;
	int 									m_ackRecv;
	TracedCallback<Ptr<const Packet>, const Address &> m_rxTrace;
public:
	ConsensusApp() : m_socket (0),
    m_local (),
    m_localSendCluster(),
    m_packetSize (0),
    m_lastTime(0.0),
    m_nextSendEvent (),
    m_running (false),
    m_isCH (false),
    m_trasactionSizeBytes(32),
    m_blockHeaderSizeBytes(68)
     {}

	virtual ~ConsensusApp(){}

	static TypeId GetTypeId (void){
		static TypeId tid = TypeId ("ConsensusApp")
		.SetParent<Application> ()
	    .SetGroupName ("Consense")
	    .AddConstructor<ConsensusApp> ()
	    ;
	    return tid;
	}

	// Para organizar o primeiro envio, devo verificar durante o cache, quem eh o nó requisitante e quem é o sender.
	// Preciso fazer esse app para validar o envio de conteúdos para membros do grupo.
	// Preciso de atualizar m_peersSockets e m_peersAddresses de acordo com a evolução do cluster.
	// Essas informações podem ser recuperadas pegando a instância de cache e cluster do sistem, preciso acrescentar algumas variáveis de controle
	//   para ser possível de detectar todas essas informações
	void Setup (uint32_t packetSize, std::vector<double> scheduleTimer, std::vector<double> clusterTimeSchedule, std::vector<double> updateTimeSchedule, Ipv4InterfaceContainer ueWiFiIface, NodeContainer ueNodes, double simTime) {
		
		m_packetSize = packetSize;
		for (std::vector<double>::iterator it = scheduleTimer.begin(); it != scheduleTimer.end(); it++) {
			m_timer.push(*it);
		}

		//for (std::vector<double>::iterator it = clusterTimeSchedule.begin(); it != clusterTimeSchedule.end(); it++) {
			Simulator::Schedule(Seconds(clusterTimeSchedule[0] + 0.1), &ConsensusApp::setPeersAddress, this);
		//}

		for (std::vector<double>::iterator it = updateTimeSchedule.begin(); it != updateTimeSchedule.end(); it++) {
			//Simulator::Schedule(Seconds(*it), &ConsensusApp::setPeersAddress, this);
		}

		//std::cout << std::endl;
		//m_timer = scheduleTimer;
		//m_tid = TcpSocketFactory::GetTypeId ();
		m_ueNodes = ueNodes;
		m_ueWiFiIface = ueWiFiIface;
		m_sendoToClusterPort = 4321;
		m_SendPort = 1234;
		m_simTime = simTime;
		m_chMessage.SetObject();
		m_alreadyScheduled = false;
		m_blockBeingFormed = false;
		m_ackRecv = 0;

		m_consensusEndTime = 0.0;
		m_consensusStartTime = 0.0;
	}

	void ScheduleNextSend(void) {
		double time;
		
		// Escalonar um evento para time segundos no futuro
		time = m_timer.front() - Simulator::Now().GetSeconds(); //- m_lastTime;
		m_lastTime = m_timer.front();

		m_timer.pop();
		//Time tNext (Seconds (100 * 8 / static_cast<double> (m_dataRate.GetBitRate ())))
		//std::cout << "Next send event - " << m_lastTime << std::endl;
		m_nextSendEvent = Simulator::Schedule(Seconds(time), &ConsensusApp::SendTransactionToClusterHead, this);
	}
	
	void ScheduleNextVerifyValidation(void) {
	
		//std::cout << "ScheduleNextVerifyValidation (" << Simulator::Now ().GetSeconds() << ")" << std::endl;

		if (!m_consensus_open && !m_alreadyScheduled){
			Simulator::Schedule(Seconds(1), &ConsensusApp::ScheduleNextVerifyValidation, this);
			if (m_conclusion.checkConclusion((int) GetNode()->GetId(), m_clusterSize)) {
			
				m_alreadyScheduled = true;
				std::cout << "validateNext (" << Simulator::Now ().GetSeconds() << ")" << std::endl;
				m_validateNextEvent = Simulator::ScheduleNow(&ConsensusApp::validateNext, this);
			}
		}
	}

	int getPeerId(Ipv4Address ip) {
		std::vector<int>::iterator it;

		for (it = m_members.begin(); it != m_members.end(); it++) {

			if (m_ueWiFiIface.GetAddress(*it) == ip) {
				return *it;
			}
		}
		return -1;
	}
	
	void schedulePing(){
		//double cont = Simulator::Now ().GetSeconds();
		
		if (m_consensus_open){
			m_nextPingEvent = Simulator::Schedule(Seconds(1), &ConsensusApp::SendPingValidation, this);
			Simulator::ScheduleNow(&ConsensusApp::SendPingVote, this);
		}
	}
	void schedulePingValidation(){
		//double cont = Simulator::Now ().GetSeconds();
		if (m_consensus_open){
			m_nextValidationPingEvent = Simulator::Schedule(Seconds(1), &ConsensusApp::schedulePingValidation, this);
			m_nextSendPingValidationEvent = Simulator::ScheduleNow(&ConsensusApp::SendPingValidation, this);
		}
		
	}

private:

	
	virtual void StartApplication (void){

		srand(time(NULL) + GetNode()->GetId());
		//std::cout << "ConsensusApp iniciado (node " << (int) GetNode()->GetId() << ")" << std::endl;
		//std::cout << "entrei" << std::endl;
		
		m_local = InetSocketAddress (Ipv4Address::GetAny (), m_SendPort);
		m_consensus_open = false;
		if(!m_socket){
			
			//m_socket = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
			m_socket = Socket::CreateSocket (GetNode (), UdpSocketFactory::GetTypeId ());
			m_socket->Bind(m_local);
			m_socket->Listen();
			m_socket->ShutdownSend();
			if (addressUtils::IsMulticast (m_local)){
				Ptr<UdpSocket> udpSocket = DynamicCast<UdpSocket> (m_socket);
				if (udpSocket) {
					udpSocket->MulticastJoinGroup (0, m_local);
				} else {
					std::cout << "Erro: multicast sem o protocolo udp" << std::endl;
				}
			}
		}
		
		m_socket->SetRecvCallback (MakeCallback (&ConsensusApp::HandleRead, this));
		m_socket->SetAcceptCallback (
		 	MakeNullCallback<bool, Ptr<Socket>, const Address &> (),
		 	MakeCallback (&ConsensusApp::HandleAccept, this));
		m_socket->SetCloseCallbacks (
		 	MakeCallback (&ConsensusApp::HandlePeerClose, this),
		 	MakeCallback (&ConsensusApp::HandlePeerError, this));

		m_localSendCluster = InetSocketAddress (Ipv4Address::GetAny (), m_sendoToClusterPort);
		if(!m_nodeClusterSocket){
			
			//m_socket = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
			m_nodeClusterSocket = Socket::CreateSocket (GetNode (), UdpSocketFactory::GetTypeId ());
			m_nodeClusterSocket->Bind(m_localSendCluster);
			m_nodeClusterSocket->Listen();
			m_nodeClusterSocket->ShutdownSend();
			if (addressUtils::IsMulticast (m_localSendCluster)){
				Ptr<UdpSocket> udpSocket = DynamicCast<UdpSocket> (m_nodeClusterSocket);
				if (udpSocket) {
					udpSocket->MulticastJoinGroup (0, m_localSendCluster);
				} else {
					std::cout << "Erro: multicast sem o protocolo udp" << std::endl;
				}
			}
		}
		
		m_nodeClusterSocket->SetRecvCallback (MakeCallback (&ConsensusApp::HandleRead, this));
		m_nodeClusterSocket->SetAcceptCallback (
		 	MakeNullCallback<bool, Ptr<Socket>, const Address &> (),
		 	MakeCallback (&ConsensusApp::HandleAccept, this));
		m_nodeClusterSocket->SetCloseCallbacks (
		 	MakeCallback (&ConsensusApp::HandlePeerClose, this),
		 	MakeCallback (&ConsensusApp::HandlePeerError, this));
		
		ScheduleNextSend();	
		
		//ScheduleNextVerifyValidation();
	}

	void firstAddressSetUp() {
		Cluster* clusterInstance =  Cluster::getInstance();
			
		m_peersAddresses.clear();
		for (std::map<Ipv4Address, Ptr<Socket>>::iterator it = m_peersSockets.begin(); m_peersSockets.end() != it; it++) {
			
			it->second->Close();
		}
		m_peersSockets.clear();

		std::vector<int>::iterator it;
		int chID = -1;
		std::vector<int> members = clusterInstance->getCluster(&chID, GetNode ()->GetId());
		m_clusterSize = (int) members.size();
		
		m_members = members;

		m_CHNodeId = chID;

		if (chID > -1){
			m_chAddress = m_ueWiFiIface.GetAddress(chID);
			//std::cout << "Node " << GetNode ()->GetId() << " clusterhead " << chID << std::endl;
		}

		if (m_CHNodeId == (int) GetNode ()->GetId()) {
			m_isCH = true;
		}

		// Recuperar membros do cluster e colocar em forma de endereços
		///std::cout << "Membros do cluster node " <<  GetNode ()->GetId() << ": ";
		for (it = members.begin(); it != members.end(); it++) {
			//std::cout << *it << " ";
			m_peersAddresses.push_back(m_ueWiFiIface.GetAddress(*it));
		}
		//std::cout << std::endl;

		//std::cout << "Node " <<  GetNode()->GetId() << " antes de criar os sockets " << std::endl;
		for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
		{
			
			m_peersSockets[*i] = Socket::CreateSocket (GetNode (), UdpSocketFactory::GetTypeId ());
			//m_peersSockets[*i] = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
			m_peersSockets[*i]->Connect (InetSocketAddress (*i, m_SendPort));
			
		}
		
		m_clusterSocket = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
		m_clusterSocket->Connect (InetSocketAddress (m_chAddress, m_sendoToClusterPort));
	

		if(m_isCH){
			for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
			{
				
				m_clusterPeersSockets[*i] = Socket::CreateSocket (GetNode (), UdpSocketFactory::GetTypeId ());
				//m_peersSockets[*i] = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
				m_clusterPeersSockets[*i]->Connect (InetSocketAddress (*i, m_sendoToClusterPort));
				
			}
		}
		clusterInstance->toggleClusterChange(false, (int) GetNode ()->GetId());
	}

	bool inWaiting(int c){
		for (auto const &p : m_waiting){
			if (p == c) return true;
		}
		return false;
	}

	void setPeersAddress() {

		Cluster* clusterInstance =  Cluster::getInstance();
		//std::cout << "ConsensusApp - Setting peers address (" << Simulator::Now ().GetSeconds() << " s) testeClusterChanged " << clusterInstance->clusterChange((int) GetNode ()->GetId()) <<std::endl;

		if (clusterInstance->clusterChange((int) GetNode ()->GetId())){
			/*std::cout << "setPeersAddress - Cluster changing!" << std::endl;
			
			std::cout << "Reseting Blockchain..." << std::endl;
			m_blockchain.ResetChain();
			std::cout << "Blockchain Reseted" << std::endl;

			std::cout << "Reseting transactionPool..." << std::endl;
			m_transactionPool.resetPool();
			std::cout << "Blockchain transactionPool" << std::endl;*/



			
			
			m_peersAddresses.clear();
			for (std::map<Ipv4Address, Ptr<Socket>>::iterator it = m_peersSockets.begin(); m_peersSockets.end() != it; it++) {
				
				it->second->Close();
			}
			m_peersSockets.clear();

			std::vector<int>::iterator it;
			int chID = -1;
			std::vector<int> members = clusterInstance->getCluster(&chID, GetNode ()->GetId());
			m_clusterSize = (int) members.size();
			
			m_members = members;

			m_CHNodeId = chID;

			if (chID > -1){
				m_chAddress = m_ueWiFiIface.GetAddress(chID);
				//std::cout << "Node " << GetNode ()->GetId() << " clusterhead " << chID << std::endl;
			}

			if (m_CHNodeId == (int) GetNode ()->GetId()) {
				m_isCH = true;
			}

			// Recuperar membros do cluster e colocar em forma de endereços
			///std::cout << "Membros do cluster node " <<  GetNode ()->GetId() << ": ";
			for (it = members.begin(); it != members.end(); it++) {
				//std::cout << *it << " ";
				m_peersAddresses.push_back(m_ueWiFiIface.GetAddress(*it));
			}
			//std::cout << std::endl;

			//std::cout << "Node " <<  GetNode()->GetId() << " antes de criar os sockets " << std::endl;
			for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
			{
				
				m_peersSockets[*i] = Socket::CreateSocket (GetNode (), UdpSocketFactory::GetTypeId ());
				//m_peersSockets[*i] = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
				m_peersSockets[*i]->Connect (InetSocketAddress (*i, m_SendPort));
				m_peersSockets[*i]->SetAllowBroadcast (true);
				
			}
			

			m_clusterSocket = Socket::CreateSocket (GetNode (), UdpSocketFactory::GetTypeId ());
			m_clusterSocket->Connect (InetSocketAddress (m_chAddress, m_sendoToClusterPort));

			if(m_isCH){
				std::cout << "CH PEERS ";
				for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
				{
					
					m_clusterPeersSockets[*i] = Socket::CreateSocket (GetNode (), UdpSocketFactory::GetTypeId ());
					//m_peersSockets[*i] = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
					m_clusterPeersSockets[*i]->Connect (InetSocketAddress (*i, m_sendoToClusterPort));
					m_clusterPeersSockets[*i]->SetAllowBroadcast(true);

					std::cout << *i << " ";
				}
				std::cout << std::endl;
			}
			
			clusterInstance->toggleClusterChange(false, (int) GetNode ()->GetId());
		}	
	}

	int validateMetadata(VideoMetadata meta) {
		if (meta.getValor() == 1) {
			return 1; // valid tx
		} else {
			return 0; // invalid tx
		}
	}
	
	void HandleRead(Ptr<Socket> socket) {
		Ptr<Packet> packet;
		Address from;
		//double receivedTime = Simulator::Now ().GetSeconds();
		//std::cout << this << socket << std::endl;

		while ((packet = socket->RecvFrom (from)))
		{
			/*std::cout << "ConsensusApp - Node " << GetNode ()->GetId() << 
			" Received one packet from " << InetSocketAddress::ConvertFrom(from).GetIpv4 () << 
			" (" << receivedTime << " s)!" << 
			std::endl;*/

			if (packet->GetSize () == 0)
			{ 
				std::cout << "Empty packet received" << std::endl;
				break;
			}
			if (InetSocketAddress::IsMatchingType (from)) {
				// std::cout << "Recebeu um pacote" << std::endl;
				std::string delimiter = "#";
				std::string parsedPacket;
				size_t pos = 0;
				char *packetInfo = new char[packet->GetSize () + 1];
				std::ostringstream totalStream;

				packet->CopyData (reinterpret_cast<uint8_t*>(packetInfo), packet->GetSize ());
				packetInfo[packet->GetSize ()] = '\0'; // ensure that it is null terminated to avoid bugs
				/**
				* Add the buffered data to complete the packet
				*/
				//m_bufferedData[from] = m_bufferedData[from] + packetInfo;
				totalStream << m_bufferedData[from] << packetInfo;
				std::string totalReceivedData(totalStream.str());

				//std::cout << "Received data: " << totalReceivedData << std::endl;

				while ((pos = totalReceivedData.find(delimiter)) != std::string::npos)
				{		
					parsedPacket = totalReceivedData.substr(0, pos);

					rapidjson::Document d;

					d.Parse(parsedPacket.c_str());
					if(!d.IsObject())
					{
						//std::cout << "The parsed packet is corrupted\n" << std::cout;
						totalReceivedData.erase(0, pos + delimiter.length());
						continue;
					}
					rapidjson::StringBuffer buffer;
					rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
					d.Accept(writer);

					/*std::cout << "At time "  << Simulator::Now ().GetSeconds ()
					<< "s node " << GetNode ()->GetId () << " received "
					<<  packet->GetSize () << " bytes from "
					<< getPeerId(InetSocketAddress::ConvertFrom(from).GetIpv4()) << " port "
					<< InetSocketAddress::ConvertFrom(from).GetPort()
					<< "\n" << buffer.GetString() << std::endl;*/
					// std::cout << "Passou aqui 1" << std::endl;
					ConsensusTimeStatics::getInstance()->addControlPacket();
					// std::cout  << "Número pacotes: " << ConsensusTimeStatics::getInstance()->getControlPackets() << std::endl;
					switch (d["message"].GetInt())
					{
						case SEND_TRANSACTION:
						{
							//std::cout << "Message Type: 'SEND_TRANSACTION'" << std::endl;
							std::string type = d["type"].GetString();
							int contID = d["content"]["id"].GetInt();
							int sender_id = d["content"]["sender"].GetInt();
							int metaVal = d["metadados"]["val"].GetInt();

							std::cout << "(SEND_TRANSACTION) CH " << GetNode ()->GetId () << " received transaction of content " << contID << " from " 
									  << getPeerId(InetSocketAddress::ConvertFrom(from).GetIpv4()) << std::endl;

							VideoMetadata meta(metaVal);
							std::cout << "Adicionando a pull de transações" << std::endl;
							// Adiciona transação na pool
							m_transactionPool.addTransacion(
								getPeerId(InetSocketAddress::ConvertFrom(from).GetIpv4()),
								InetSocketAddress::ConvertFrom(from).GetIpv4(), 
								sender_id,
								m_ueWiFiIface.GetAddress(sender_id),
								contID,
								Simulator::Now ().GetSeconds(),
								Simulator::Now ().GetSeconds(),
								metaVal
							);

							std::cout << "TX HASH: " << m_transactionPool.getNextTransaction().getTxHash() << std::endl;
							if (!m_consensus_open) {
								validateNext();
								//ScheduleNextVerifyValidation();
							} else {
								std::cout << "Validate another time\n" << std::endl;
							}
							//m_transactionPool.printQueue();
							break;
						}

						case VALIDATE_TRANSACTION:
						{
							m_consensus_open = true;
							std::cout << buffer.GetString() << std::endl;

							std::cout << "\n(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") VALIDATE_TRANSACTION" << std::endl;
							std::string type = d["type"].GetString();
							int receiverID = d["transaction"]["receiver"].GetInt(); // Recupera o receptor do conteúdo relativo a essa transação
							int senderID = d["transaction"]["sender"].GetInt(); // Recupera o sender do conteúdo relativo a essa transação
							double timestamp =  d["transaction"]["timestamp"].GetDouble(); // Recupera o timestamp ao qual o clusterhead recebeu a mensagem
							double timeRecv = d["transaction"]["timeRecv"].GetDouble();
							int contID = d["transaction"]["contentID"].GetInt(); // Recupera o ID do conteúdo relativo a essa transaçã
							int metaVal = d["transaction"]["metadata"]["val"].GetInt(); // Recupera o campo de metadados

							int chvote = d["chvote"].GetInt();
							Ipv4Address addr = m_ueWiFiIface.GetAddress(receiverID);
							Ipv4Address sender_addr = m_ueWiFiIface.GetAddress(senderID);
							
							VideoMetadata meta(metaVal);

							Transaction rcv_trans(receiverID, addr, senderID, sender_addr, contID, timestamp, timeRecv, meta);
							std::string rcv_txHash = rcv_trans.getTxHash();
							std::string pckt_hash = d["transaction"]["hash"].GetString();

							m_validationSession.setHash(pckt_hash);
							m_validationSession.addTransacion(rcv_trans);

							//std::cout << "Received Hash: " << rcv_txHash;
							//std::cout << "\nProduced Hash: " << pckt_hash << std::endl;
							
							if (rcv_trans.validateTxHash(pckt_hash)) {
								//std::cout << "(VALIDATE_TRANSACTION) Same hash value!" << std::endl;

								std::string recv_blockHash = d["blockHash"].GetString();

								Transaction transaction = m_validationSession.getTransaction();
								int blockHeight = m_blockchain.GetTheTop()->GetBlockHeight() + 1;
								int parentBlockContentID = m_blockchain.GetTheTop()->GetContentID();
								std::string parentBlockHash = m_blockchain.GetTheTop() ->GetBlockHash();
								int blockSizeBytes = m_blockHeaderSizeBytes + m_trasactionSizeBytes; // 100 bytes
								double timeCreated = transaction.getTime();
								double timeReceived = transaction.getTimeRecv();
								Ipv4Address recvFrom = m_ueWiFiIface.GetAddress(transaction.getRecvID());

								Block newBlock (transaction, 
												chvote, 
												blockHeight, 
												contID, 
												parentBlockContentID, 
												parentBlockHash,
												blockSizeBytes, 
												timeCreated, 
												timeReceived, 
												recvFrom, 
												transaction.getRecvID(), 
												sender_addr, transaction.getSenderID());

								
								//std::cout << "Correct Block Hash" << std::endl;
								m_validationSession.addVotingPool((int) GetNode ()->GetId(), (int) GetNode ()->GetId(), validateMetadata(meta));
								m_validationSession.addVotingPool((int) GetNode ()->GetId(),getPeerId(m_chAddress), chvote);
							

							} else {
								std::cout << "(VALIDATE_TRANSACTION) Wrong hash value!" << std::endl;
								m_validationSession.desclassifyNode();
								//m_validationSession.addVotingPool((int) GetNode ()->GetId(), (int) GetNode ()->GetId(), 0);
								
							}

							
							/*if (rcv_trans.validateTxHash(pckt_hash)) {
								//std::cout << "(VALIDATE_TRANSACTION) Same hash value!" << std::endl;
								d.RemoveMember("type"); // Trocar o tipo de mensagem
								rapidjson::Value 		value;
								value.SetString("vote");
								d.AddMember("type", value, d.GetAllocator());

								rapidjson::Value 		vote(rapidjson::kObjectType);
								vote = validateMetadata(meta);

								m_validationSession.addVotingPool((int) GetNode ()->GetId(), (int) GetNode ()->GetId(), validateMetadata(meta));
								d.AddMember("vote", vote, d.GetAllocator());

							} else {
								std::cout << "(VALIDATE_TRANSACTION) Wrong hash value!" << std::endl;
								d.RemoveMember("type"); // Trocar o tipo de mensagem

								rapidjson::Value 		value;
								value.SetString("vote");
								d.AddMember("type", value, d.GetAllocator());

								rapidjson::Value 		vote(rapidjson::kObjectType);
								vote = 0;

								m_validationSession.addVotingPool((int) GetNode ()->GetId(), (int) GetNode ()->GetId(), 0);
								d.AddMember("vote", vote, d.GetAllocator());
							}*/
							m_pingValidationSet.eraseSet(m_CHNodeId, m_ueWiFiIface.GetAddress(GetNode ()->GetId ()));

							if (!m_pingValidationSet.isEmpty(m_CHNodeId)) {
								std::cout << "(Node "<< GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<")" 
								<< " m_pingValidationSet ACK from " << getPeerId(InetSocketAddress::ConvertFrom(from).GetIpv4()) << " falta: " << std::endl;
								
								for (auto const &p: m_pingValidationSet.getSet(m_CHNodeId)) 
								{
									std::cout << getPeerId(p) << " ";
								}
								std::cout << std::endl;
								//saveCHping(d);
								//scheduleCheckValidationStart();
							} else {
								SendMessage(VALIDATE_TRANSACTION, ACK_CH, d, m_clusterSocket);

								//schedulePing();
							}
							/*m_pingSet.insert(m_chAddress);
							for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
							{
								if (getPeerId(*i) != (int) GetNode()->GetId ())
									m_pingSet.insert(*i);
							}
							
							
							for (auto const &p: m_pingSet)
							{
								m_peersSockets[p]->Send (reinterpret_cast<const uint8_t*>(m_buffer.GetString()), m_buffer.GetSize(), 0);
								m_peersSockets[p]->Send (delimiter, 1, 0);
								
							}		*/	
							
							//m_transactionPool.printQueue();
							
							break;
						}

						case ACK_CH:
						{
							//Simulator::Cancel(m_nextValidationPingEvent);
							//Simulator::Cancel(m_nextSendPingValidationEvent);
							if (m_consensus_open){
								std::cout << "(Node "<< GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") ACK_CH" << std::endl;
								std::cout << buffer.GetString() << std::endl;
								rapidjson::StringBuffer buff;
									
								rapidjson::Writer<rapidjson::StringBuffer> writer(buff);

								d["message"].SetInt(START_VOTE);
								d.Accept(writer);

								if (m_isCH){
									for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
									{
										if (getPeerId(*i) != (int) GetNode()->GetId()){
											const uint8_t delimiter[] = "#";
											m_clusterPeersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(buff.GetString()), buff.GetSize(), 0);
											m_clusterPeersSockets[*i]->Send (delimiter, 1, 0);
										}
									}
								} else{
									for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
									{
										if (getPeerId(*i) != (int) GetNode()->GetId()){
											const uint8_t delimiter[] = "#";
											m_peersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(buff.GetString()), buff.GetSize(), 0);
											m_peersSockets[*i]->Send (delimiter, 1, 0);
										}
									}
								}
							}
								

							break;
						}

						case START_VOTE:
						{
							if (m_consensus_open){
								std::cout << "(Node "<< GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") START_VOTE" << std::endl;
								std::cout << buffer.GetString() << std::endl;

								Transaction transaction = m_validationSession.getTransaction();
								std::string type = d["type"].GetString();
								int metaVal = transaction.getMetadata().getValor();

								VideoMetadata meta(metaVal);

								d.RemoveMember("type"); // Trocar o tipo de mensagem


								rapidjson::Value 		value;
								value.SetString("vote");
								d.AddMember("type", value, d.GetAllocator());

								int nodeMetadataDecision = validateMetadata(meta);

								rapidjson::Value 		vote(rapidjson::kObjectType);
								vote = validateMetadata(meta);
								d.AddMember("vote", vote, d.GetAllocator());

								Trust* trustInstance = Trust::getInstance(nNodes);
								if (trustInstance->hasPreviousContact(GetNode ()->GetId (), transaction.getSenderID())){
									if (trustInstance->is_trustable( GetNode ()->GetId (), transaction.getSenderID())){
										//voto
										rapidjson::Value 		nodeOpinion(rapidjson::kObjectType);
										nodeOpinion = 1;
										d.AddMember("nodeOpinion", nodeOpinion, d.GetAllocator());
									} else {
										rapidjson::Value 		nodeOpinion(rapidjson::kObjectType);
										nodeOpinion = 0;
										d.AddMember("nodeOpinion", nodeOpinion, d.GetAllocator());
									}
								} else {
									rapidjson::Value 		nodeOpinion(rapidjson::kObjectType);
									nodeOpinion = -1;
									d.AddMember("nodeOpinion", nodeOpinion, d.GetAllocator());
								} 
								// if (trustInstance->is_trustable_in_consensus( GetNode ()->GetId (), transaction.getSenderID())){
								// 	//voto
								// 	rapidjson::Value 		nodeOpinion(rapidjson::kObjectType);
								// 	nodeOpinion = 1;
								// 	d.AddMember("nodeOpinion", nodeOpinion, d.GetAllocator());
								// } else if (trustInstance->even_trustable_consensus( GetNode ()->GetId (), transaction.getSenderID())){
								// 	rapidjson::Value 		nodeOpinion(rapidjson::kObjectType);
								// 	nodeOpinion = -1;
								// 	d.AddMember("nodeOpinion", nodeOpinion, d.GetAllocator());
								// } else {
								// 	rapidjson::Value 		nodeOpinion(rapidjson::kObjectType);
								// 	nodeOpinion = 0;
								// 	d.AddMember("nodeOpinion", nodeOpinion, d.GetAllocator());
								// }

								
								int contID = transaction.getContID();
								int blockHeight = m_blockchain.GetTheTop()->GetBlockHeight() + 1;
								int parentBlockContentID = m_blockchain.GetTheTop()->GetContentID();
								std::string parentBlockHash = m_blockchain.GetTheTop() ->GetBlockHash();
								int blockSizeBytes = m_blockHeaderSizeBytes + m_trasactionSizeBytes; // 100 bytes
								double timeCreated = transaction.getTime();
								double timeReceived = transaction.getTimeRecv();
								Ipv4Address recvFrom = m_ueWiFiIface.GetAddress(transaction.getRecvID());
								Ipv4Address sender_addr = m_ueWiFiIface.GetAddress(transaction.getSenderID());

								Block newBlock (transaction, 
												nodeMetadataDecision, 
												blockHeight, 
												contID, 
												parentBlockContentID, 
												parentBlockHash,
												blockSizeBytes, 
												timeCreated, 
												timeReceived, 
												recvFrom, 
												transaction.getRecvID(),
												sender_addr, transaction.getSenderID());

								d.RemoveMember("blockHash");

								value.SetString(newBlock.GetBlockHash().c_str(), newBlock.GetBlockHash().length(), d.GetAllocator());
								d.AddMember("blockHash", value, d.GetAllocator()); // d["transaction"]["hash"]
								
								/*if (rcv_trans.validateTxHash(pckt_hash)) {
									//std::cout << "(VALIDATE_TRANSACTION) Same hash value!" << std::endl;
									d.RemoveMember("type"); // Trocar o tipo de mensagem
									rapidjson::Value 		value;
									value.SetString("vote");
									d.AddMember("type", value, d.GetAllocator());

									rapidjson::Value 		vote(rapidjson::kObjectType);
									vote = validateMetadata(meta);

									//m_validationSession.addVotingPool((int) GetNode ()->GetId(), (int) GetNode ()->GetId(), validateMetadata(meta));
									d.AddMember("vote", vote, d.GetAllocator());
									m_validationSession.addVotingPool((int) GetNode ()->GetId(), (int) GetNode ()->GetId(), 1);

								} else {
									std::cout << "(VALIDATE_TRANSACTION) Wrong hash value!" << std::endl;
									d.RemoveMember("type"); // Trocar o tipo de mensagem

									rapidjson::Value 		value;
									value.SetString("vote");
									d.AddMember("type", value, d.GetAllocator());

									rapidjson::Value 		vote(rapidjson::kObjectType);
									vote = 0;

									//m_validationSession.addVotingPool((int) GetNode ()->GetId(), (int) GetNode ()->GetId(), 0);
									d.AddMember("vote", vote, d.GetAllocator());
									m_validationSession.addVotingPool((int) GetNode ()->GetId(), (int) GetNode ()->GetId(), 0);
								}*/

								rapidjson::StringBuffer buff;
								rapidjson::Writer<rapidjson::StringBuffer> writer(buff);

								d["message"].SetInt(VOTE);
								d.Accept(writer);

								if (m_isCH){
									for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
									{
										if (getPeerId(*i) != (int) GetNode()->GetId()){
											const uint8_t delimiter[] = "#";
											m_clusterPeersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(buff.GetString()), buff.GetSize(), 0);
											m_clusterPeersSockets[*i]->Send (delimiter, 1, 0);
										}
									}
								} else{
									for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
									{
										if (getPeerId(*i) != (int) GetNode()->GetId()){
											const uint8_t delimiter[] = "#";
											m_peersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(buff.GetString()), buff.GetSize(), 0);
											m_peersSockets[*i]->Send (delimiter, 1, 0);
										}
									}
								}
							}
							break;

						}
						case VOTE:
						{
							if (m_consensus_open){
								std::string type = d["type"].GetString();
								// std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") Message Type: 'VOTE'" << std::endl;
								// std::cout << buffer.GetString() << std::endl;
								//Simulator::Cancel(m_nextCheckValidationStartEvent);

								
								int vote = d["vote"].GetInt();
								int nodeOpinion = d["nodeOpinion"].GetInt();
								// std::cout << nodeOpinion << std::endl;
								int receivedFrom = getPeerId(InetSocketAddress::ConvertFrom(from).GetIpv4());
								m_trust_pool.addVotingPool((int) GetNode ()->GetId (), receivedFrom, nodeOpinion);
								//m_trust_pool.print((int) GetNode ()->GetId ());
								int nodeId = getPeerId(InetSocketAddress::ConvertFrom(from).GetIpv4());

								std::string recv_blockHash = d["blockHash"].GetString();

								Transaction transaction = m_validationSession.getTransaction();
								int contID = transaction.getContID();
								int blockHeight = m_blockchain.GetTheTop()->GetBlockHeight() + 1;
								int parentBlockContentID = m_blockchain.GetTheTop()->GetContentID();
								std::string parentBlockHash = m_blockchain.GetTheTop() ->GetBlockHash();
								int blockSizeBytes = m_blockHeaderSizeBytes + m_trasactionSizeBytes; // 100 bytes
								double timeCreated = transaction.getTime();
								double timeReceived = transaction.getTimeRecv();
								Ipv4Address recvFrom = m_ueWiFiIface.GetAddress(transaction.getRecvID());
								Ipv4Address sender_addr = m_ueWiFiIface.GetAddress(transaction.getSenderID());
								Block recv_block (transaction, 
												vote, 
												blockHeight, 
												contID, 
												parentBlockContentID, 
												parentBlockHash,
												blockSizeBytes, 
												timeCreated, 
												timeReceived, 
												recvFrom, 
												transaction.getRecvID(),
												sender_addr, transaction.getSenderID());
								m_validationSession.addVotingPool((int) GetNode ()->GetId(),nodeId, vote);
								// if (recv_block.GetBlockHash().compare(recv_blockHash) == 0){
								// 	//std::cout << "Correct Block Hash" << std::endl;
								// 	m_validationSession.addVotingPool((int) GetNode ()->GetId(),nodeId, vote);
								// } else {
								// 	std::cout << "Incorrect Block Hash" << std::endl;
								// 	m_validationSession.desclassifyNode();
								// 	//m_validationSession.addVotingPool((int) GetNode ()->GetId(), (int) GetNode ()->GetId(), 0);
								// }
								
								
								//m_validationSession.printVotePool((int) GetNode ()->GetId ());
								//std::cout << std::endl;
	/*
								m_pingSet.eraseSet
								(
									getPeerId(InetSocketAddress::ConvertFrom(from).GetIpv4()),
									m_ueWiFiIface.GetAddress (GetNode ()->GetId ())
								);*/

								//SendMessage(VOTE, ACK_VOTE, d, m_peersSockets[InetSocketAddress::ConvertFrom(from).GetIpv4()]);
								
								if (m_validationSession.verifyConclusion((int) GetNode ()->GetId (), m_peersAddresses, m_ueWiFiIface, m_CHNodeId)) {
									
									std::cout << "(Node "<< GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") END OF VOTATION" << std::endl;
									//m_validationSession.validateTx((int) GetNode ()->GetId ());
									m_blockConsenus.setClusterSize((int) m_clusterSize);

									if(m_validationSession.getNodeDecision((int) GetNode ()->GetId ()) == 1) {
										std::cout << "Consider cont " << contID <<  " VALID!\n";
										/*d.RemoveMember("vote");
										rapidjson::Value 	decision;
										decision = 1;
										m_validationSession.setDecision(1);


										d.AddMember("decision", decision, d.GetAllocator());

										rapidjson::StringBuffer buff;					
										rapidjson::Writer<rapidjson::StringBuffer> writer(buff);
										d["message"].SetInt(GET_NEW_BLOCK);
										d.Accept(writer);*/

										std::cout << "Todas as decisões recebidas. O BLOCO PODE SER FORMADO!" << std::endl;

										Block newBlock (transaction, 
														vote, 
														blockHeight, 
														contID, 
														parentBlockContentID, 
														parentBlockHash,
														blockSizeBytes, 
														timeCreated, 
														timeReceived, 
														recvFrom, 
														transaction.getRecvID(),
														sender_addr, transaction.getSenderID());

										if (m_blockchain.HasBlock(contID)){
											//std::cout << "Bloco já adicionado" << std::endl;
											break;
										}
										std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") blockchain insertion" << std::endl;
									
										m_blockchain.AddBlock(newBlock);
										m_blockchain.Print();

										std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<")\n";
										
										m_validationSession.resetValidationSession((int) GetNode ()->GetId ());
										Content* contentInstance =  Content::getInstance();
									
										BlockchainConnection* blockchainConn =  BlockchainConnection::getInstance();

										int nodeId = (int) GetNode()->GetId();
										if(m_blockchain.GetTheTop()->GetDecision() == 1){
											if (m_blockchain.GetTheTop()->GetReceivedFromId() == nodeId){
												
												blockchainConn->insert(nodeId, contID, 1); // blockchainCOnnector to outside

												// contentInstance->tempCache[nodeId].pop();
												// contentInstance->tempCacheSenders[nodeId].pop();
												contentInstance->addCache(nodeId, m_blockchain.GetTheTop()->GetContentID());
											}
										} else {
											// contentInstance->tempCache[nodeId].pop();
											// contentInstance->tempCacheSenders[nodeId].pop();
											blockchainConn->insert(nodeId, contID, 0);
										} 
										m_waiting.erase(m_blockchain.GetTheTop()->GetContentID());

										std::cout << std::endl;
										if (!m_isCH){
											m_consensus_open = false;
											setPeersAddress();
										}

										SendMessage(VOTE, FINISH, d, m_clusterSocket);

									} else {
										std::cout << "Consider cont " << contID <<  " INVALID!\n";
										/*d.RemoveMember("vote");
										rapidjson::Value 	decision;
										decision = 0;
										m_validationSession.setDecision(0);
										d.AddMember("decision", decision, d.GetAllocator());

										rapidjson::StringBuffer buff;					
										rapidjson::Writer<rapidjson::StringBuffer> writer(buff);
										d["message"].SetInt(GET_NEW_BLOCK);
										d.Accept(writer);*/
										std::cout << "Todas as decisões recebidas. O BLOCO PODE SER FORMADO!" << std::endl;
										Block newBlock (transaction, 
														vote, 
														blockHeight, 
														contID, 
														parentBlockContentID, 
														parentBlockHash,
														blockSizeBytes, 
														timeCreated, 
														timeReceived, 
														recvFrom, 
														transaction.getRecvID(),
														sender_addr, transaction.getSenderID());

										if (m_blockchain.HasBlock(contID)){
											//std::cout << "Bloco já adicionado" << std::endl;
											break;
										}
										std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") blockchain insertion" << std::endl;
									
										m_blockchain.AddBlock(newBlock);
										m_blockchain.Print();

										std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<")\n";
										
										m_validationSession.resetValidationSession((int) GetNode ()->GetId ());
										Content* contentInstance =  Content::getInstance();
									
										BlockchainConnection* blockchainConn =  BlockchainConnection::getInstance();

										int nodeId = (int) GetNode()->GetId();
										if(m_blockchain.GetTheTop()->GetDecision() == 1){
											if (m_blockchain.GetTheTop()->GetReceivedFromId() == nodeId){
												
												blockchainConn->insert(nodeId, contID, 1); // blockchainCOnnector to outside

												// contentInstance->tempCache[nodeId].pop();
												// contentInstance->tempCacheSenders[nodeId].pop();
												contentInstance->addCache(nodeId, m_blockchain.GetTheTop()->GetContentID());
												
											}
										} else {
											// contentInstance->tempCache[nodeId].pop();
											// contentInstance->tempCacheSenders[nodeId].pop();
											blockchainConn->insert(nodeId, contID, 0);
										} 
										m_waiting.erase(m_blockchain.GetTheTop()->GetContentID());
										std::cout << std::endl;
										
										if (!m_isCH){
											m_consensus_open = false;
											setPeersAddress();
										}
										SendMessage(VOTE, FINISH, d, m_clusterSocket);
									}
									

									/*m_pingSet.resetSession((int) GetNode()->GetId ());
									m_pingValidationSet.resetSession((int) GetNode()->GetId ());*/
								}
							}
								
							
							/*	
								Disclaimer						
									Aqui devo enviar a mensagem com as informações da blockchain também
							*/
							
							// Temporariamente termina o periodo de consenso aqui

							break;
						}
						case FINISH:
						{
							m_ackRecv++;
							
							if (m_ackRecv >= m_clusterSize){
								setPeersAddress();
								validateNext();
								m_consensusEndTime = (double) Simulator::Now().GetSeconds();

								std::cout << "Consensus finished in " << m_consensusEndTime - m_consensusStartTime << " s\n";
								auto consensusIsntance = ConsensusTimeStatics::getInstance();
								consensusIsntance->add(m_clusterSize, m_consensusEndTime - m_consensusStartTime);

								m_consensusEndTime = 0.0;
								m_consensusStartTime = 0.0;

								m_consensus_open = false;
								std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") Consensus finished !" << "\n";
								//ScheduleNextVerifyValidation();
								// Atualiza confiança indireta
								Trust* trustInstance = Trust::getInstance(nNodes);
								
								for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
								{
									int srcNode = getPeerId(*i);
									int sender = m_blockchain.GetTheTop()->GetSenderId();
									for (auto const &p: m_trust_pool.getVotePool(srcNode)){
										// std::cout << "First -> " << p.first << "\t Second -> " << p.second << std::endl;
										if (p.second == 1){
											trustInstance->addEvidence(srcNode, sender, p.first, true);
										}
										else if (p.second == 0){
											trustInstance->addEvidence(srcNode, sender, p.first, false);
										}
										else {
											// Não teve contato
										}

									}
									m_trust_pool.reset(srcNode);
									std::cout << "Updated indirect trust node (" << srcNode << ") -> ("<< sender <<") " <<  trustInstance->getIndirectTrust(srcNode, sender) << std::endl;
									std::cout << "Updated trust node (" << srcNode << ") -> ("<< sender <<") " <<  trustInstance->getTrust(srcNode, sender) << std::endl;
								}
							}

							break;
						}
						case GET_NEW_BLOCK:
						{
							/*std::cout << "At time "  << Simulator::Now ().GetSeconds ()
							<< "s clusterHead " << GetNode ()->GetId () << " received "
							<< "GET_NEW_BLOCK from "
							<< getPeerId(InetSocketAddress::ConvertFrom(from).GetIpv4())
							<< " with decision " << d["decision"].GetInt() << std::endl;
							//std::cout << buffer.GetString() << std::endl;*/
							if (m_consensus_open){
								m_validationSession.newNodeDecision(InetSocketAddress::ConvertFrom(from).GetIpv4(), d["decision"].GetInt());
								m_validationSession.printAux();

								if (m_validationSession.allDecisionsReceived()){
									//Simulator::Cancel(m_nextSendPingValidationEvent);
									
									std::cout << "Todas as decisões recebidas. O BLOCO PODE SER FORMADO! \n" << std::endl;
									//m_pingValidationSet.resetSession((int) GetNode()->GetId ());

									Transaction transaction = m_validationSession.getTransaction();
									int blockHeight = m_blockchain.GetTheTop()->GetBlockHeight() + 1;
									int contID = transaction.getContID();
									int parentBlockContentID = m_blockchain.GetTheTop()->GetContentID();
									std::string parentBlockHash = m_blockchain.GetTheTop() ->GetBlockHash();
									int blockSizeBytes = m_blockHeaderSizeBytes + m_trasactionSizeBytes; // 100 bytes
									double timeCreated = Simulator::Now ().GetSeconds ();
									double timeReceived = timeCreated;
									Ipv4Address recvFrom = m_ueWiFiIface.GetAddress(transaction.getRecvID());

									Block newBlock (transaction, 
													d["decision"].GetInt(), 
													blockHeight, 
													contID, 
													parentBlockContentID, 
													parentBlockHash,
													blockSizeBytes, 
													timeCreated, 
													timeReceived, 
													recvFrom, 
													transaction.getRecvID());
									std::cout << "Block Hash: " << newBlock.GetBlockHash() << std::endl;

									rapidjson::Document block;
									block.SetObject();

									rapidjson::Value value;
									rapidjson::Value blockInfo(rapidjson::kObjectType);

									value.SetString("block");
									block.AddMember("type", value, block.GetAllocator());

									value = BLOCK;
									block.AddMember("message", value, block.GetAllocator());

									value = timeCreated;
									blockInfo.AddMember("timeCreated", value, block.GetAllocator()); 
									value  = timeReceived;
									blockInfo.AddMember("timeReceived", value, block.GetAllocator()); 

									value.SetString(newBlock.GetBlockHash().c_str(), newBlock.GetBlockHash().length(), block.GetAllocator());
									blockInfo.AddMember("blockHash", value, block.GetAllocator());
									block.AddMember("block", blockInfo, block.GetAllocator());

									rapidjson::StringBuffer buff;
									rapidjson::Writer<rapidjson::StringBuffer> writer(buff);
									block["message"].SetInt(BLOCK); // Manda uma mensagem BLOCK
									block.Accept(writer);

									if (m_isCH){
										for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
										{
											const uint8_t delimiter[] = "#";
											m_clusterPeersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(buff.GetString()), buff.GetSize(), 0);
											m_clusterPeersSockets[*i]->Send (delimiter, 1, 0);
										}
									} else{
										for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
										{
											const uint8_t delimiter[] = "#";
											m_peersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(buff.GetString()), buff.GetSize(), 0);
											m_peersSockets[*i]->Send (delimiter, 1, 0);
										}
									}

									//m_transactionPool.printQueue();	
								}
							}
							break;
						}
						
						case BLOCK:
						{
							/*std::cout << "At time "  << Simulator::Now ().GetSeconds ()
							<< "s node " << GetNode ()->GetId () << " received "
							<< "BLOCK from clusterhead "
							<< getPeerId(InetSocketAddress::ConvertFrom(from).GetIpv4()) << std::endl;*/
							if (m_consensus_open){
								
								//std::cout << "(Node "<< GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") BLOCK" << std::endl;
								//std::cout << buffer.GetString() << std::endl;
								
								std::string type = d["type"].GetString();
								//std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") Message Type: 'BLOCK'" << std::endl;

								std::string recv_blockHash = d["block"]["blockHash"].GetString(); // Recupera o receptor do conteúdo relativo a essa transação
								double timestamp =  d["block"]["timeCreated"].GetDouble(); // Recupera o timestamp ao qual o clusterhead recebeu a mensagem

								Transaction transaction = m_validationSession.getTransaction();
								int blockHeight = m_blockchain.GetTheTop()->GetBlockHeight() + 1;
								int contID = transaction.getContID();
								int parentBlockContentID = m_blockchain.GetTheTop()->GetContentID();
								int blockSizeBytes = m_blockHeaderSizeBytes + m_trasactionSizeBytes; // 100 bytes
								std::string parentBlockHash = m_blockchain.GetTheTop()->GetBlockHash();
								
								//std::cout << "parentBlockHash " << parentBlockHash << std::endl;

								Ipv4Address recvFrom = m_ueWiFiIface.GetAddress(transaction.getRecvID());

								Block newBlock (transaction, 
												m_validationSession.getNodeDecision((int) GetNode ()->GetId ()), 
												blockHeight, 
												contID, 
												parentBlockContentID, 
												parentBlockHash,
												blockSizeBytes, 
												timestamp, 
												timestamp, 
												recvFrom, 
												transaction.getRecvID());

								//std::cout << "Received hash: " << recv_blockHash << "\nBlock hash: " << newBlock.GetBlockHash() << std::endl;

								/*if (newBlock.VerifyHash(recv_blockHash)){
									std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") formed a good block" << std::endl;
									if (newBlock.GetBlockHeight() > m_blockchain.GetBlockchainHeight()){
										m_blockchain.AddBlock(newBlock);
										m_blockchain.Print();
									}
								}*/

								rapidjson::Value 		value;
								//rapidjson::Value      array(rapidjson::kArrayType);
								d.RemoveMember("block");
								d.RemoveMember("message");

								d["type"] = "block_vote";

								value = BLOCK_VOTE;
								d.AddMember("message", value, d.GetAllocator());

								//transaction
								rapidjson::Value 		transInfo(rapidjson::kObjectType);
								value = getPeerId(transaction.getIP());
								transInfo.AddMember("receiver", value, d.GetAllocator()); // final["transaction"]["receiver"]
								value  = transaction.getTime();
								transInfo.AddMember("timestamp", value, d.GetAllocator()); // final["transaction"]["timestamp"]
								value = transaction.getContID();
								transInfo.AddMember("contentID", value, d.GetAllocator()); // final["transaction"]["contentID"]
								value.SetString(transaction.getTxHash().c_str(), transaction.getTxHash().length(), d.GetAllocator());
								transInfo.AddMember("hash", value, d.GetAllocator()); // final["transaction"]["hash"]

								//metadado
								rapidjson::Value 		metadata(rapidjson::kObjectType);
								value = transaction.getMetadata().getValor();
								metadata.AddMember("val", value, d.GetAllocator()); // final["transaction"]["metadata"]["val"]
								transInfo.AddMember("metadata", metadata, d.GetAllocator());
								
								d.AddMember("transaction", transInfo, d.GetAllocator());
								
								rapidjson::Value blockInfo(rapidjson::kObjectType);

								value = blockHeight;
								blockInfo.AddMember("blockHeight", value, d.GetAllocator()); 
								
								value =	contID;
								blockInfo.AddMember("contID", value, d.GetAllocator()); 

								value =	parentBlockContentID;
								blockInfo.AddMember("parentBlockContentID", value, d.GetAllocator()); 

								value.SetString(parentBlockHash.c_str(), parentBlockHash.length(), d.GetAllocator());
								blockInfo.AddMember("parentBlockHash", value, d.GetAllocator()); 
								

								value = blockSizeBytes;
								blockInfo.AddMember("blockSizeBytes", value, d.GetAllocator()); 

								value = timestamp;
								blockInfo.AddMember("timeCreated", value, d.GetAllocator()); 

								value  = timestamp;
								blockInfo.AddMember("timeReceived", value, d.GetAllocator()); 

								// Sem recvFrom -- devemos capturar ao receber a mensagem
								value = transaction.getRecvID();
								blockInfo.AddMember("receivedFromId", value, d.GetAllocator()); 
								
								value.SetString(newBlock.GetBlockHash().c_str(), newBlock.GetBlockHash().length(), d.GetAllocator());
								blockInfo.AddMember("blockHash", value, d.GetAllocator());

								d.AddMember("block", blockInfo, d.GetAllocator());

								rapidjson::StringBuffer buff;
								rapidjson::Writer<rapidjson::StringBuffer> writer(buff);
								d["message"].SetInt(BLOCK_VOTE);
								d.Accept(writer);
							
								if (m_isCH){
									for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
									{
										const uint8_t delimiter[] = "#";
										m_clusterPeersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(buff.GetString()), buff.GetSize(), 0);
										m_clusterPeersSockets[*i]->Send (delimiter, 1, 0);
									}
								} else{
									for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
									{
										const uint8_t delimiter[] = "#";
										m_peersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(buff.GetString()), buff.GetSize(), 0);
										m_peersSockets[*i]->Send (delimiter, 1, 0);
									}
								}
							}
							
							break;
						}

						case BLOCK_VOTE:
						{
							/*std::cout << "At time "  << Simulator::Now ().GetSeconds ()
							<< "s node " << GetNode ()->GetId () << " received "
							<< "BLOCK_VOTE from "
							<< getPeerId(InetSocketAddress::ConvertFrom(from).GetIpv4()) << std::endl;*/
							
							//std::cout << buffer.GetString() << std::endl;
							
							//std::cout << "(Node "<< GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") BLOCK_VOTE" << std::endl;

							int receiverID = d["transaction"]["receiver"].GetInt(); // Recupera o receptor do conteúdo relativo a essa transação
							int senderID = d["transaction"]["sender"].GetInt(); // Recupera o receptor do conteúdo relativo a essa transação
							//std::cout << "receiverID " << receiverID << std::endl;
							double tx_timestamp =  d["transaction"]["timestamp"].GetDouble(); // Recupera o timestamp ao qual o clusterhead recebeu a mensagem
							//std::cout << "tx_timestamp " << tx_timestamp << std::endl;
							double timeRecv = d["transaction"]["timeRecv"].GetDouble();
							int tx_contID = d["transaction"]["contentID"].GetInt(); // Recupera o ID do conteúdo relativo a essa transação
							//std::cout << "tx_contID " << tx_contID << std::endl;
							int metaVal = d["transaction"]["metadata"]["val"].GetInt(); // Recupera o campo de metadados
							//std::cout << "metaVal " << metaVal << std::endl;
							std::string tx_hash = d["transaction"]["hash"].GetString();
							//std::cout << "tx_hash " << tx_hash << std::endl;
							
							Ipv4Address receiver = m_ueWiFiIface.GetAddress(m_validationSession.getTransaction().getRecvID());
							Ipv4Address recvFrom = m_ueWiFiIface.GetAddress(m_validationSession.getTransaction().getSenderID());
							VideoMetadata meta(metaVal);

							//std::cout << "Passou aqui 1" << std::endl;

							Transaction trans(receiverID, receiver, senderID, recvFrom, tx_contID, tx_timestamp, timeRecv, meta);

							std::string recv_blockHash = d["block"]["blockHash"].GetString(); // Recupera o receptor do conteúdo relativo a essa transação

							double timestamp =  d["block"]["timeCreated"].GetDouble(); // Recupera o timestamp ao qual o clusterhead recebeu a mensagem
							int blockHeight = d["block"]["blockHeight"].GetInt();
							int contID = d["block"]["contID"].GetInt();
							int parentBlockContentID = d["block"]["parentBlockContentID"].GetInt();
							int blockSizeBytes = d["block"]["blockSizeBytes"].GetInt();
							std::string parentBlockHash = d["block"]["parentBlockHash"].GetString();
							

							Block newBlock (trans, 
											m_validationSession.getNodeDecision((int) GetNode ()->GetId ()), 
											blockHeight, 
											contID, 
											parentBlockContentID, 
											parentBlockHash,
											blockSizeBytes, 
											timestamp, 
											timestamp, 
											recvFrom, 
											receiverID);

							m_blockBeingFormed = true;
							// Verifica se o hash fornecido pelo bloco compartilhado é o mesmo do hash compartilhado
							if (m_blockchain.HasBlock(contID)){
								//std::cout << "Bloco já adicionado" << std::endl;
								break;
							}
							if (newBlock.VerifyHash(recv_blockHash)){
								m_blockConsenus.addNewBlock(recv_blockHash, newBlock);
								
								//m_blockConsenus.print();
								
								Block b; 
								if (m_blockConsenus.verifyConclusion(&b)) {
									
									std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") blockchain insertion" << std::endl;
									
									m_blockchain.AddBlock(b);
									m_blockchain.Print();

									std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<")\n";
									
									m_validationSession.resetValidationSession((int) GetNode ()->GetId ());
									m_blockConsenus.reset();
									
									m_pingValidationSet.reset((int) GetNode ()->GetId ());
									m_conclusion.add(m_CHNodeId, InetSocketAddress::ConvertFrom(from).GetIpv4());
									
									
									Content* contentInstance =  Content::getInstance();
									
									BlockchainConnection* blockchainConn =  BlockchainConnection::getInstance();

									int nodeId = (int) GetNode()->GetId();
									if(m_blockchain.GetTheTop()->GetDecision() == 1){
										if (m_blockchain.GetTheTop()->GetReceivedFromId() == nodeId){
											
											blockchainConn->insert(nodeId, contID, 1); // blockchainCOnnector to outside

											contentInstance->tempCache[nodeId].pop();
											contentInstance->tempCacheSenders[nodeId].pop();
											contentInstance->addCache(nodeId, m_blockchain.GetTheTop()->GetContentID());
											m_waiting.erase(m_blockchain.GetTheTop()->GetContentID());
										}
									} else {
										contentInstance->tempCache[nodeId].pop();
										contentInstance->tempCacheSenders[nodeId].pop();
										blockchainConn->insert(nodeId, contID, 0);
									} 

									std::cout << std::endl;
									m_consensus_open = false;
									setPeersAddress();
									if (m_isCH){
										validateNext();
										//ScheduleNextVerifyValidation();
									}
									/*m_pingSet.reset((int) GetNode ()->GetId ());
									
									std::cout << "\t-> m_pingSet: ";
									for (auto const &p: m_pingSet.getSet((int)GetNode ()->GetId ())) 
									{
										std::cout << getPeerId(p) << " ";
									}
									std::cout << std::endl;*/
									

									//validateNext();
									
								} else {
									//m_blockConsenus.print();
								}

							} else {
								m_blockConsenus.desclassifyNode();
							}	
						
							break;
						}
						default:{
							break;
						}
					}
					
					totalReceivedData.erase(0, pos + delimiter.length());
				}
				m_bufferedData[from] = totalReceivedData;
				delete[] packetInfo;
			}
			m_rxTrace (packet, from);
		}
	}

	virtual void StopApplication (void){
		m_running = false;
		for (std::vector<Ipv4Address>::iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) //close the outgoing sockets
		{
			m_peersSockets[*i]->Close ();
		}
		if (m_nextSendEvent.IsRunning ())
		{
			Simulator::Cancel (m_nextSendEvent);
		}

		while(!m_socketList.empty ()) //these are accepted sockets, close them
		{
			Ptr<Socket> acceptedSocket = m_socketList.front ();
			m_socketList.pop_front ();
			acceptedSocket->Close ();
		}
		if (m_socket) 
		{
			m_socket->Close ();
			m_socket->SetRecvCallback (MakeNullCallback<void, Ptr<Socket> > ());
		}
	}

	void SendPingVote(void){

		const uint8_t delimiter[] = "#";
		for (auto const &p: m_pingSet.getSet((int)GetNode ()->GetId ())) 
		{
			//std::cout << "SendPingVote - " << p << std::endl;

			m_peersSockets[p]->Send (reinterpret_cast<const uint8_t*>(m_buffer.GetString()), m_buffer.GetSize(), 0);
			m_peersSockets[p]->Send (delimiter, 1, 0);
			
		}
		/*for (auto const &p: m_pingSet) 
		{
			//std::cout << "SendPingVote - " << p << std::endl;
			m_peersSockets[p]->Send (reinterpret_cast<const uint8_t*>(m_buffer.GetString()), m_buffer.GetSize(), 0);
			m_peersSockets[p]->Send (delimiter, 1, 0);
			
		}*/
	}
	void SendPingValidation(void) {
		const uint8_t delimiter[] = "#";
		
		for (auto const &p: m_pingValidationSet.getSet((int)GetNode ()->GetId ())) 
		{
			//std::cout << "SendPingValidation - " << p << std::endl;
			
			m_peersSockets[p]->Send (reinterpret_cast<const uint8_t*>(m_bufferValidation.GetString()), m_bufferValidation.GetSize(), 0);
			m_peersSockets[p]->Send (delimiter, 1, 0);
			
			
		}
	}
	
	void validateNext(void) {
		m_transactionPool.printQueue();
		//Simulator::Cancel(m_validateNextEvent);

		if (!m_transactionPool.empty() && m_isCH) {
			/*m_validationSession.resetValidationSession((int) GetNode ()->GetId ());
			m_blockConsenus.reset();*/
			//m_pingSet.setHash();
			m_ackRecv = 0;
			m_conclusion.reset((int) GetNode ()->GetId());

			m_alreadyScheduled = false;
			m_consensus_open = true;

			m_consensusStartTime = (double) Simulator::Now().GetSeconds();


			std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId << ") validateNext - Validation process at "<< Simulator::Now().GetSeconds() << "s" << std::endl;
			Transaction t = m_transactionPool.getNextTransaction();// Recupera transação da pool


			Transaction transaction = m_transactionPool.getNextTransaction();
			int blockHeight = m_blockchain.GetTheTop()->GetBlockHeight() + 1;
			int contID = transaction.getContID();
			int parentBlockContentID = m_blockchain.GetTheTop()->GetContentID();
			std::string parentBlockHash = m_blockchain.GetTheTop() ->GetBlockHash();
			int blockSizeBytes = m_blockHeaderSizeBytes + m_trasactionSizeBytes; // 100 bytes
			double timeCreated = transaction.getTime();
			double timeReceived = transaction.getTimeRecv();
			Ipv4Address recvFrom = m_ueWiFiIface.GetAddress(transaction.getSenderID());


			m_validationSession.addVotingPool((int) GetNode ()->GetId(), (int) GetNode ()->GetId(), validateMetadata(t.getMetadata()));
			m_validationSession.setHash(t.getTxHash());
			m_validationSession.addTransacion(t);

			std::cout << "Start validation process of transaction " << t.getContID() << " from peer " << getPeerId(t.getIP()) << std::endl;

			m_transactionPool.popTransaction();

			rapidjson::Document 	d;
			rapidjson::Value 		value(VALIDATE_TRANSACTION);
			//rapidjson::Value      array(rapidjson::kArrayType);
			d.SetObject();
			d.AddMember("message", value, d.GetAllocator());

			value.SetString("tx");
			d.AddMember("type", value, d.GetAllocator());

			//transaction
			rapidjson::Value 		transInfo(rapidjson::kObjectType);
			value = getPeerId(t.getIP());
			transInfo.AddMember("receiver", value, d.GetAllocator()); // d["transaction"]["receiver"]
			value = t.getSenderID();
			transInfo.AddMember("sender", value, d.GetAllocator()); // d["transaction"]["sender"]
			value  = t.getTime();
			transInfo.AddMember("timestamp", value, d.GetAllocator()); // d["transaction"]["timestamp"]
			value  = t.getTimeRecv();
			transInfo.AddMember("timeRecv", value, d.GetAllocator()); // d["transaction"]["timestamp"]
			value = t.getContID();
			transInfo.AddMember("contentID", value, d.GetAllocator()); // d["transaction"]["contentID"]
			value.SetString(t.getTxHash().c_str(), t.getTxHash().length(), d.GetAllocator());
			transInfo.AddMember("hash", value, d.GetAllocator()); // d["transaction"]["hash"]

			
			//metadado
			rapidjson::Value 		metadata(rapidjson::kObjectType);
			value = t.getMetadata().getValor();
			metadata.AddMember("val", value, d.GetAllocator()); // d["transaction"]["metadata"]["val"]
			transInfo.AddMember("metadata", metadata, d.GetAllocator());
			
			d.AddMember("transaction", transInfo, d.GetAllocator());

			//voto
			rapidjson::Value 		val_cont(rapidjson::kObjectType);
			val_cont = validateMetadata(t.getMetadata());
			d.AddMember("chvote", val_cont, d.GetAllocator());


			Block newBlock (transaction, 
							validateMetadata(t.getMetadata()), 
							blockHeight, 
							contID, 
							parentBlockContentID, 
							parentBlockHash,
							blockSizeBytes, 
							timeCreated, 
							timeReceived, 
							recvFrom, 
							transaction.getRecvID());
			

			value.SetString(newBlock.GetBlockHash().c_str(), newBlock.GetBlockHash().length(), d.GetAllocator());
			d.AddMember("blockHash", value, d.GetAllocator()); // d["transaction"]["hash"]

			rapidjson::StringBuffer buffer;
			//m_bufferValidation.Clear();
			rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);

			d["message"].SetInt(VALIDATE_TRANSACTION);
			d.Accept(writer);
			
			for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
			{
				if (getPeerId (*i) != (int) GetNode()->GetId()){
				
					m_pingValidationSet.addPingSet((int) GetNode()->GetId (), *i);
				}
				
				m_validationSession.newVotingNode(*i);
			}

			std::cout << "validateNext - m_pingValidationSet ";
			for (auto const &p: m_pingValidationSet.getSet(m_CHNodeId))
			{
				std::cout << getPeerId (p) << " ";
			}	 
			std::cout << std::endl;
/*
			std::cout << "validateNext - m_pingSet ";
			for (auto const &p: m_pingValidationSet.getSet(m_CHNodeId))
			{
				std::cout << p << " ";
			}
			std::cout << std::endl;	 
*/

		
			if (m_isCH){
				for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
				{
					const uint8_t delimiter[] = "#";
					m_clusterPeersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(buffer.GetString()), buffer.GetSize(), 0);
					m_clusterPeersSockets[*i]->Send (delimiter, 1, 0);
				}
			} else{
				for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
				{
					const uint8_t delimiter[] = "#";
					m_peersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(buffer.GetString()), buffer.GetSize(), 0);
					m_peersSockets[*i]->Send (delimiter, 1, 0);
				}
			}
			
			//schedulePingValidation();
		}
	}

	void SendTransactionToClusterHead (void){
		//setPeersAddress();
		ScheduleNextSend();
		Content* instance = Content::getInstance();
		double sendTime =  Simulator::Now ().GetSeconds();
		int node_id = (int)GetNode()->GetId();
		MaliciousNodes* maliciousInstance = MaliciousNodes::getInstance();
				//InetSocketAddress remote = InetSocketAddress (ueWiFiIface.GetAddress (*vecIt, 0), port);
		if (m_CHNodeId >= 0 && node_id != m_CHNodeId){
			if(!instance->tempCache[node_id].isEmpty() && !instance->tempCacheSenders[node_id].isEmpty()){
				// int decision;
				// if (m_blockchain.HasBlock((int) instance->tempCache[node_id].getNext(), &decision) ){
				// 	std::cout << "(Node " << GetNode ()->GetId () << ", Cluster " << m_CHNodeId <<") Content already validated !" << std::endl;

				// 	if (decision == 1){
				// 		int c = instance->tempCache[node_id].pop();
				// 		instance->addCache((int) GetNode()->GetId(), c);
				// 	} else {
				// 		std::cout << "\t-> Invalid !" << std::endl;
				// 		instance->tempCache[node_id].pop();
				// 		instance->tempCacheSenders[node_id].pop();
				// 	}
				// } else
				if(!inWaiting(instance->tempCache[node_id].getNext())){
					Ptr<Packet> packet = Create<Packet> (m_packetSize);
					std::cout << "Cache Senders: ";
					instance->tempCacheSenders[node_id].print();
					std::cout << "Cache Content: ";
					instance->tempCache[node_id].print();
					// Enviar para o ch	
					std::cout << "\n(" << sendTime << " s) ConsensusApp - Node " << GetNode()->GetId() << " enviando transação do conteúdo " 
							  << instance->tempCache[node_id].getNext() << " (originado de " << instance->tempCacheSenders[node_id].getNext() << ") para o CH " << m_CHNodeId
							  << std::endl;

					rapidjson::Document 	d;
					rapidjson::Value 		value(SEND_TRANSACTION);
					//rapidjson::Value      array(rapidjson::kArrayType);
					d.SetObject();
					d.AddMember("message", value, d.GetAllocator());

					value.SetString("transaction");
					d.AddMember("type", value, d.GetAllocator());

					rapidjson::Value 		contInfo(rapidjson::kObjectType);

					
					value = instance->tempCache[node_id].getNext();
					
					m_waiting.insert((int) instance->tempCache[node_id].pop());

					contInfo.AddMember("id",  value, d.GetAllocator());	// Adiciona a informação acerca do conteúdo

					value = instance->tempCacheSenders[node_id].getNext();
					contInfo.AddMember("sender", value, d.GetAllocator());

					d.AddMember("content", contInfo, d.GetAllocator()); // Adiciona campo no arquivo json
					rapidjson::Value 		metaInfo(rapidjson::kObjectType);
					int sender_id = instance->tempCacheSenders[node_id].pop();
					if (maliciousInstance->isMalicious(sender_id) && !basicConfigPBFTInstance->maliciousSwitched(sendTime)){
						// malicious
						value = 0;
					} else {
						// not malicious
						value = 1;
					}
					
					metaInfo.AddMember("val", value, d.GetAllocator()); // Adiciona a informação extraída dos metadados
					d.AddMember("metadados", metaInfo, d.GetAllocator()); // Adiciona campo no arquivo json

					
					//SendMessage(NO_MESSAGE, SEND_TRANSACTION, d, m_peersSockets[*i]);
					SendMessage(NO_MESSAGE, SEND_TRANSACTION, d, m_clusterSocket);
					std::cout << "Passou aqui 2" << std::endl;
					//m_peersSockets[*i]->Send(packet);
						
					//std::cout << "Packet " << m_packetsSent << " sent at " << Simulator::Now().GetSeconds() << std::endl;
					//m_socket->Send (packet);
					instance->lastContAdded = -1; 
					instance->lastNodeCont = -1;
				} 
			} else {
				// std::cout << "Nothing to validate" << std::endl;
			}
		} else {
			//std::cout << "Node não associado a um cluster\n" << std::endl;
		}
		

	}

	void HandlePeerClose (Ptr<Socket> socket)
	{
	  //NS_LOG_FUNCTION (this << socket);
	}

	void HandlePeerError (Ptr<Socket> socket)
	{
	  //NS_LOG_FUNCTION (this << socket);
	}

	void HandleAccept (Ptr<Socket> s, const Address& from)
	{
	  //NS_LOG_FUNCTION (this << s << from);
	  s->SetRecvCallback (MakeCallback (&ConsensusApp::HandleRead, this));
	  m_socketList.push_back (s);
	}


	void SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, rapidjson::Document &d, Ptr<Socket> outgoingSocket)
	{
		//NS_LOG_FUNCTION (this);

		const uint8_t delimiter[] = "#";

		rapidjson::StringBuffer buffer;
		rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);

		d["message"].SetInt(responseMessage);
		d.Accept(writer);

		/*std::cout << "Node " << GetNode ()->GetId () << " got a "
		<< getMessageName(receivedMessage) << " message"
		<< " and sent a " << getMessageName(responseMessage)
		<< " message: " << buffer.GetString() << std::endl;
*/
		outgoingSocket->Send (reinterpret_cast<const uint8_t*>(buffer.GetString()), buffer.GetSize(), 0);
		outgoingSocket->Send (delimiter, 1, 0);
	}



};

