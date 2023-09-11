#include "ns3/address-utils.h"
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "../rapidjson/document.h"
#include "../rapidjson/writer.h"
#include "../rapidjson/stringbuffer.h"
#include "block.h"
#include <map>
#include <vector>

using namespace ns3;

class NodeApp : public Application
{
public:
	static TypeId GetTypeId (void);

	NodeApp();
	NodeApp(Ipv4InterfaceContainer ueWiFiIface);
	NodeApp(Ipv4InterfaceContainer ueWiFiIface, std::string ecdsa_publicKey, std::string ecdsa_privateKey, std::string ecdh_publicKey, std::string ecdh_privateKey);
	NodeApp(Ipv4InterfaceContainer ueWiFiIface, int idCH, std::vector<int> cluster);
	NodeApp(Ipv4InterfaceContainer ueWiFiIface, int idCH, std::vector<int> cluster, std::string ecdsa_publicKey, std::string ecdsa_privateKey, std::string ecdh_publicKey, std::string ecdh_privateKey);

 	virtual ~NodeApp ();

	/**
	  *  \return pointer to listen sockets
	**/
	Ptr<Socket> GetListeningSocket (void) const;

	/**
	  *  \return a list of the clusterMembers addresses
	**/
	std::vector<Ipv4Address> GetMmebersAddresses (void) const;

	/**
	  *  \return the node's ecdsa public key
	**/
	std::string GetECDSAPublickey (void) const;
	/**
	  *  \return the node's ecdsa private key
	**/
	std::string GetECDSAPrivatekey (void) const;

	/**
	  *  \return the node's ecdsa public key
	**/
	std::string GetECDHPublickey (void) const;

	/**
	  *  \return the node's private key
	**/
	std::string GetECDHPrivatekey (void) const;

	/*****PAREEEEI AQUI******/
	/**
	  *  \return the clusterhead public key
	**/
	std::string GetCHecdsaPublickey (void) const;

	/**
	  *  \return the clusterhead public key
	**/
	std::string GetCHecdhPublickey (void) const;

 	/**
 	  *  \brief update and change the node's keys
 	  *  \param new private ecdsa key and public key to be updated.
 	**/
 	void SetNewEcdsaKeys(std::string newPrivKey, std::string newPublicKey);
 	void SetNewEcdsaKeys();

 	/**
 	  *  \brief update and change the node's keys
 	  *  \param new private ecdh key and public key to be updated.
 	**/
 	void SetNewEcdhKeys(std::string newPrivKey, std::string newPublicKey);
 	void SetNewEcdhKeys();

 	/**
 	  *  \brief Update the CH public key
 	  *  \param new public key of the CH
 	*/
 	void SetCHecdsaPublickey (std::string  newCHKey);
 	void SetCHecdsaPublickey ();

	/**
 	  *  \brief Update the CH public key
 	  *  \param new public key of the CH
 	*/
 	void SetCHecdhPublickey (std::string  newCHKey);
 	void SetCHecdhPublickey ();

	/**
	  *  \brief Set the addresses of clusterMembers
	  *  \param peers the reference of a vector containing the Ipv4 addresses of clustermembers
	**/
	void SetMembersAddresses (const std::vector<Ipv4Address> &peers);
	void SetMembersAddresses ();

	/**
	  * \brief Add new member to the set of members
	  * \param newPeer the new member to add
	**/
	void AddMember (Ipv4Address newPeer);

	void AddMember (int newPeer);

	/**
	  *  \brief Update the node cluster
	  *  \param newPeer the new member to bo added
	**/
	void UpdateCluster (int newCH, std::vector<int> members);

	/**
	  *  \brief Set the protocol type(default: STANDARD_PROTOCOL)
	  *  \param protocolType the type of protocol used for advertising 
	**/
	void SetProtocolType (enum ProtocolType protocolType);

	/**
	  *  \brief Sends a message to a node
	  *  \param receivedMessage the type of the received message
	  *  \param responseMessage the type of the response message
	  *  \param d the rapidjson document containing the info of the outgoing message
	  *  \param outgoingSocket the socket of the peer that will receive the message
	**/
	void SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, rapidjson::Document &d, Ptr<Socket> outgoingSocket);

	
	// void SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, rapidjson::Document &d, Address &outgoingAddress);

	
	void SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, std::string packet, Address &outgoingAddress);

	
	void SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, std::string packet, Ipv4Address &outgoingIpv4Address);

	
	void SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, std::string packet, Ptr<Socket> outgoingSocket);

protected:
	virtual void StartApplication (void); 
	virtual void StopApplication (void); 

	virtual void DoDispose (void);

	/**
	  * \brief Handle a packet received by the application
	  * \param socket the receiving socket
	**/
	void HandleRead (Ptr<Socket> socket);
// 	/**
// 	  * \brief Handle an incoming connection
// 	  * \param socket the incoming connection socket
// 	  * \param from the address the connection is from
// 	**/
 	void HandleAccept (Ptr<Socket> socket, const Address& from);

// 	/**
// 	  * \brief Handle an connection close
// 	  * \param socket the connected socket
// 	**/
 	void HandlePeerClose (Ptr<Socket> socket);

// 	/**
// 	  * \brief Handle an connection error
// 	  * \param socket the connected socket
// 	**/
 	void HandlePeerError (Ptr<Socket> socket);

// 	/**
// 	  * \brief Check signature integrity of a message
// 	  * \param message message to check integrity
// 	  * \param sign signature of the message
// 	**/
// 	bool checkSign (std::string message, std::string sign);

// 	/**
// 	  * \brief Check signature integrity of a message
// 	  * \param message message to check integrity
// 	  * \param sign signature of the message
// 	  * \param sender address of the message sender
// 	**/
// 	bool checkSign (std::string message, std::string sign, Ipv4Address sender);

// 	/**
// 	  * \brief decrypt a encrypted message
// 	  * \param message message to be decrypted
// 	**/
// 	std::string decrypt (std::string message);

// 	/**
// 	  * \brief encrypt a message to be sent
// 	  * \param message message to be encrypted
// 	  * \param PublicKey public key of the receiver
// 	**/
// 	std::string encrypt (std::string message,std::string publicKey);

// 	/**
// 	  * \brief encrypt a message to be sent
// 	  * \param message message to be encrypted
// 	  * \param receiver receiver's ip address
// 	  * \param type of encryption - aes or rsa
// 	**/
// 	std::string encrypt (std::string message, Ipv4Address receiver, int type);

// 	/**
// 	  * \brief encrypt a message to be sent
// 	  * \param message message to be encrypted
// 	  * \param aes_key AES session key
// 	  * \param iv IV for the aes
// 	**/
// 	std::string encrypt (std::string message, byte* aes_key, byte* iv);

private:
	Ptr<Socket> m_socket;									// Node scoket listening
	Address m_local;										// Local address to bind to
	TypeId m_tid;											// Protocol TypeId
	int m_numberOfPeers;									// Number of node's peers
	int m_CHNodeId;											// CHNode id
	
	std::string m_CHecdsaPublicKey;							// CH ecdsa Public Key
	std::string m_CHecdhPublicKey;							// CH ecdh Public Key

	Ipv4Address m_CHNodeAddress;							// CHNode address
	std::string m_ECDSAprivateKey;							// Node ECDSA private key
	std::string m_ECDSApublicKey;							// Node ECDSA public key
	std::string m_ECDSACHPublicKey;							// CH ECDSA node public key

	std::string m_ECDHprivateKey;							// Node ECDH private key
	std::string m_ECDHpublicKey;							// Node ECDH public key
	std::string m_ECDHCHPublicKey;							// CH ECDH node public key

	std::vector<Ipv4Address> m_peersAddresses;				// Vector of Ipv4Address 
	std::map<Ipv4Address, std::string> m_publicEcdsaKeys;	// Map of nodes EcdsaPublicKeys
	std::map<Ipv4Address, std::string> m_publicEcdhKeys;	// Map of nodes EcdhPublicKeys
	std::map<Ipv4Address, Ptr<Socket>> m_peersSockets;		// Map holding peerSockets
	std::map<Address, std::string> m_bufferedData;
	enum ProtocolType m_protocolType;						// Protocol type
	const int m_Port;										// Porta de envio
	nodeStatistics *m_nodeStats;

	Ipv4InterfaceContainer m_ueWiFiIface;
	TracedCallback<Ptr<const Packet>, const Address &> m_rxTrace;

};

TypeId NodeApp::GetTypeId (void) {
	static TypeId tid = TypeId ("NodeApp")
						.SetParent<Application> ()
						.SetGroupName("Applications")
						.AddConstructor<NodeApp> ()
						.AddAttribute("Local",
										"The Address on which to Bind the rx socket.",
										AddressValue (),
										MakeAddressAccessor (&NodeApp::m_local),
										MakeAddressChecker ())
						.AddAttribute("Protocol",
										"The Address on which to Bind the rx socket.",
										TypeIdValue (UdpSocketFactory::GetTypeId ()),
										MakeTypeIdAccessor (&NodeApp::m_tid),
										MakeTypeIdChecker ())
						.AddTraceSource("Rx",
										"A packet has been received",
										MakeTraceSourceAccessor (&NodeApp::m_rxTrace),
										"ns3::Packet::AddressTracedCallback")
	;
	return tid;
}
void NodeApp::DoDispose (void)
{
        //NS_LOG_FUNCTION (this);
        m_socket = 0;

        // chain up
        Application::DoDispose ();
}

NodeApp::NodeApp () : m_Port (4321)
{
    //NS_LOG_FUNCTION (this);
    
    m_socket = 0;
    m_CHNodeId = 0;
    m_CHNodeAddress = Ipv4Address ("::1");
    m_numberOfPeers = m_peersAddresses.size();
    // Recupera as chaves ECDSA e ECDH relacionadas ao nó
    SetNewEcdsaKeys();
    SetNewEcdhKeys();
    /*DEV_ECDSA ecdsa;
    std::string privateKey;
    std::string publicKey;
	ecdsa.GenerateKeys(privateKey, publicKey);

    m_ECDSAprivateKey = privateKey;
    m_ECDSApublicKey = publicKey;

    DEV_ECDH ecdh;
    ecdh.GenerateKeys(1024,privateKey, publicKey);


    m_ECDHprivateKey = privateKey;
    m_ECDHpublicKey = publicKey;*/

    // Recupera as chaves públicas do CH
    SetCHecdsaPublickey();
    SetCHecdhPublickey();

}

NodeApp::NodeApp (Ipv4InterfaceContainer ueWiFiIface) : m_Port (4321)
{
    //NS_LOG_FUNCTION (this);
    m_socket = 0;
    m_nodeStats->ch = false;
    m_CHNodeId = 0;
    m_CHNodeAddress = Ipv4Address ("::1");
    m_numberOfPeers = m_peersAddresses.size();

    DEV_ECDSA ecdsa;
    std::string privateKey;
    std::string publicKey;
	ecdsa.GenerateKeys(privateKey, publicKey);

    m_ECDSAprivateKey = privateKey;
    m_ECDSApublicKey = publicKey;

    DEV_ECDH ecdh;
    ecdh.GenerateKeys(1024,privateKey, publicKey);

    m_ECDHprivateKey = privateKey;
    m_ECDHpublicKey = publicKey;
    m_ueWiFiIface = ueWiFiIface;

}

NodeApp::NodeApp(Ipv4InterfaceContainer ueWiFiIface, 
				std::string ecdsa_publicKey, 
				std::string ecdsa_privateKey,
				std::string ecdh_publicKey, 
				std::string ecdh_privateKey) : m_Port (4321)
{
	m_CHNodeId = 0;	
	m_nodeStats->ch = false;

	m_CHNodeAddress = Ipv4Address ("::1");
	m_numberOfPeers = m_peersAddresses.size();
	m_socket = 0;
	m_ECDSAprivateKey = ecdsa_privateKey;
    m_ECDSApublicKey = ecdsa_publicKey;

    m_ECDHprivateKey = ecdh_privateKey;
    m_ECDHpublicKey = ecdh_publicKey;
    m_ueWiFiIface = ueWiFiIface;

}

NodeApp::NodeApp(Ipv4InterfaceContainer ueWiFiIface, int idCH, std::vector<int> cluster) : m_Port (4321) 
{
	m_CHNodeId = idCH;	
	if (m_CHNodeId == (int) GetNode()->GetId()) {
		    m_nodeStats->ch = true;
	}
	m_CHNodeAddress = ueWiFiIface.GetAddress(idCH);
	m_numberOfPeers = 0;
	m_ueWiFiIface = ueWiFiIface;
	m_socket = 0;
	for (std::vector<int>::iterator cluster_it = cluster.begin(); cluster_it != cluster.end(); cluster_it++) {
		m_peersAddresses.push_back(m_ueWiFiIface.GetAddress(*cluster_it));
	}
	m_numberOfPeers = m_peersAddresses.size();
	DEV_ECDSA ecdsa;
    std::string privateKey;
    std::string publicKey;
	ecdsa.GenerateKeys(privateKey, publicKey);

    m_ECDSAprivateKey = privateKey;
    m_ECDSApublicKey = publicKey;

    DEV_ECDH ecdh;
    ecdh.GenerateKeys(1024,privateKey, publicKey);

    m_ECDHprivateKey = privateKey;
    m_ECDHpublicKey = publicKey;
}	

NodeApp::NodeApp(Ipv4InterfaceContainer ueWiFiIface, 
				int idCH, 
				std::vector<int> cluster, 
				std::string ecdsa_publicKey, 
				std::string ecdsa_privateKey,
				std::string ecdh_publicKey, 
				std::string ecdh_privateKey) :  m_Port (4321) 
{
	m_CHNodeId = idCH;	
	if (m_CHNodeId == (int)GetNode()->GetId()) {
		m_nodeStats->ch = true;
	}
	m_CHNodeAddress = ueWiFiIface.GetAddress(idCH);
	m_ueWiFiIface = ueWiFiIface;
	m_socket = 0;
	m_ECDSAprivateKey = ecdsa_privateKey;
    m_ECDSApublicKey = ecdsa_publicKey;

    m_ECDHprivateKey = ecdh_privateKey;
    m_ECDHpublicKey = ecdh_publicKey;

    for (std::vector<int>::iterator cluster_it = cluster.begin(); cluster_it != cluster.end(); cluster_it++) {
		m_peersAddresses.push_back(m_ueWiFiIface.GetAddress(*cluster_it));
	}
	m_numberOfPeers = m_peersAddresses.size();

}

NodeApp::~NodeApp ()
{
    //NS_LOG_FUNCTION (this);
}

Ptr<Socket> NodeApp::GetListeningSocket (void) const {
	//NS_LOG_FUNCTION (this);
	return m_socket;
}

std::vector<Ipv4Address> NodeApp::GetMmebersAddresses (void) const {
	//NS_LOG_FUNCTION (this);
	return m_peersAddresses;
}

std::string NodeApp::GetECDSAPublickey (void) const {
	//NS_LOG_FUNCTION (this);
	return m_ECDSApublicKey;
}

std::string NodeApp::GetECDSAPrivatekey (void) const {
	//NS_LOG_FUNCTION (this);
	return m_ECDSAprivateKey;
}

std::string NodeApp::GetECDHPublickey (void) const {
	//NS_LOG_FUNCTION (this);
	return m_ECDHpublicKey;
}

std::string NodeApp::GetECDHPrivatekey (void) const {
	//NS_LOG_FUNCTION (this);
	return m_ECDHprivateKey;
}

std::string NodeApp::GetCHecdsaPublickey (void) const {
	return m_CHecdsaPublicKey;
}
std::string NodeApp::GetCHecdhPublickey(void) const {
	return m_CHecdhPublicKey;
}

void NodeApp::SetNewEcdsaKeys(std::string newPrivKey, std::string newPublicKey){
	m_ECDSAprivateKey = newPrivKey;
	m_ECDSApublicKey = newPublicKey;
}

void NodeApp::SetNewEcdsaKeys(){
	Key *keyUtil = Key::getInstance();

	m_ECDSAprivateKey = keyUtil->getECDSA_privateKeyFromMap(GetNode()->GetId());
	m_ECDSApublicKey = keyUtil->getECDSA_publicKeyFromMap(GetNode()->GetId());
}

void NodeApp::SetNewEcdhKeys(std::string newPrivKey, std::string newPublicKey){
	m_ECDHprivateKey = newPrivKey;
	m_ECDHpublicKey = newPublicKey;
}

void NodeApp::SetNewEcdhKeys(){
	Key *keyUtil = Key::getInstance();

	m_ECDHprivateKey = keyUtil->getECC_privateKeyFromMap(GetNode()->GetId());
	m_ECDHpublicKey = keyUtil->getECC_publicKeyFromMap(GetNode()->GetId());
}

void NodeApp::SetCHecdsaPublickey (std::string  newCHKey) {
	m_CHecdsaPublicKey = newCHKey;
}

void NodeApp::SetCHecdsaPublickey () {
	Cluster* clusterUtil =  Cluster::getInstance();
	Key *keyUtil = Key::getInstance();
	int newCH;
	clusterUtil->isMember(&newCH, GetNode()->GetId());
	m_CHecdsaPublicKey = keyUtil->getECDSA_publicKeyFromMap(newCH);
}


void NodeApp::SetCHecdhPublickey (std::string  newCHKey) {
	m_CHecdhPublicKey = newCHKey;
}

void NodeApp::SetCHecdhPublickey () {
	Cluster* clusterUtil =  Cluster::getInstance();
	Key *keyUtil = Key::getInstance();
	int newCH;
	clusterUtil->isMember(&newCH, GetNode()->GetId());
	m_CHecdhPublicKey = keyUtil->getECC_publicKeyFromMap(newCH);
}

void NodeApp::SetMembersAddresses (const std::vector<Ipv4Address> &peers) {
	m_peersAddresses = peers;
	m_numberOfPeers = m_peersAddresses.size();
}

void NodeApp::SetMembersAddresses () {
	Cluster* clusterUtil =  Cluster::getInstance();
	std::vector<int> peersId; 
	clusterUtil->getMembers(m_CHNodeId, &peersId);
	for (std::vector<int>::iterator it = peersId.begin(); it != peersId.end(); it++) {
		m_peersAddresses.push_back(m_ueWiFiIface.GetAddress(*it));
	}
}

void NodeApp::AddMember (Ipv4Address newPeer) {
	m_peersAddresses.push_back(newPeer);
	m_numberOfPeers = m_peersAddresses.size();
}

void NodeApp::AddMember (int newPeer) {
	m_peersAddresses.push_back(m_ueWiFiIface.GetAddress (newPeer));
	m_numberOfPeers = m_peersAddresses.size();
}

void NodeApp::UpdateCluster (int newCH, std::vector<int> members) {
	m_CHNodeId = newCH;
	if (m_nodeStats->nodeId == m_CHNodeId) {
		m_nodeStats->ch = true;
	} else {
		m_nodeStats->nodeId = false;
	}
	m_peersAddresses.clear();
	for (std::vector<int>::iterator cluster_it = members.begin(); cluster_it != members.end(); cluster_it++) {
		m_peersAddresses.push_back(m_ueWiFiIface.GetAddress(*cluster_it));
	}
}

void NodeApp::SetProtocolType (enum ProtocolType protocolType) {
	m_protocolType = protocolType;
}

void NodeApp::StartApplication() {
	srand(time(NULL) + GetNode()->GetId());
	std::cout << "Node " << GetNode()->GetId() << ": Cluster Head node = " << m_CHNodeId << std::endl;
	std::cout << "Node " << GetNode()->GetId() << ": My peers are ";

	for (auto it = m_peersAddresses.begin(); it != m_peersAddresses.end(); it++) {
		std::cout << " " << *it;
	}
	std::cout << std::endl;

	if(!m_socket){
		m_socket = Socket::CreateSocket (GetNode (), m_tid);
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
	m_socket->SetRecvCallback (MakeCallback (&NodeApp::HandleRead, this));
	m_socket->SetAcceptCallback (
	 	MakeNullCallback<bool, Ptr<Socket>, const Address &> (),
	 	MakeCallback (&NodeApp::HandleAccept, this));
	m_socket->SetCloseCallbacks (
	 	MakeCallback (&NodeApp::HandlePeerClose, this),
	 	MakeCallback (&NodeApp::HandlePeerError, this));

	std::cout << "Node" <<  GetNode()->GetId() << " before creating sockets " << std::endl;
	for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
	{
		m_peersSockets[*i] = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
		m_peersSockets[*i]->Connect (InetSocketAddress (*i, m_Port));
	}
	std::cout << "Node" <<  GetNode()->GetId() << " after creating sockets " << std::endl;
	m_nodeStats->nodeId = GetNode ()->GetId();
}

void NodeApp::HandleRead(Ptr<Socket> socket) {
	Ptr<Packet> pckt;
	Address from;


	while ((pckt = socket->RecvFrom (from))) {
		if (pckt->GetSize() == 0) {
			std::cout << "Node " << GetNode ()->GetId () << " received empty pckt: " << std::endl;
			break;
		}
		if (InetSocketAddress::IsMatchingType(from)) {
			std::string delimitador = "#";
			std::string parsedPckt;
			size_t pos = 0;
			char *pcktInfo = new char[pckt->GetSize () + 1];
			
			pckt->CopyData (reinterpret_cast<uint8_t*>(pcktInfo), pckt->GetSize ());
			pcktInfo[pckt->GetSize()] = '\0';

			std::ostringstream totalStream;
			totalStream << m_bufferedData[from] << pcktInfo;
			std::string totalReceived(totalStream.str());

			std::cout << "Node " << GetNode ()->GetId () << " received data: " << totalReceived << std::endl;
			// Iterating through  over all data received
			while ((pos = totalReceived.find(delimitador)) != std::string::npos) {
				parsedPckt = totalReceived.substr(0, pos);

				rapidjson::StringBuffer buffer;
				rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
				rapidjson::Document doc;
				doc.Parse(parsedPckt.c_str());
				if (!doc.IsObject()) {
					std::cout << "Corrupted packet: " << parsedPckt << std::endl;
					totalReceived.erase(0, pos + delimitador.length());
					continue;
				}

				std::cout << "At time "  << Simulator::Now ().GetSeconds ()
                       << "s node " << GetNode ()->GetId () << " received "
                       <<  pckt->GetSize () << " bytes from "
                       << InetSocketAddress::ConvertFrom(from).GetIpv4 ()
                       << " port " << InetSocketAddress::ConvertFrom (from).GetPort ()
                       << " with info = " << buffer.GetString() << std::endl;

                switch (doc["message"].GetInt()) {
                	case RECEIVE_PUBLIC_EDSA_KEY:
                	{
                		std::string msgDelimitador = "/";
                		std::string parsed = doc["key"].GetString();
                		size_t keyPos = parsed.find(msgDelimitador);
                		Ipv4Address nodeAddress = InetSocketAddress::ConvertFrom(from).GetIpv4();

                		std::string publicKeyStr = parsed.substr(0,keyPos).c_str();
                		m_publicEcdsaKeys[nodeAddress] = publicKeyStr;
                	}

                	case RECEIVE_PUBLIC_EDH_KEY:
                	{
                		std::string msgDelimitador = "/";
                		std::string parsed = doc["key"].GetString();
                		size_t keyPos = parsed.find(msgDelimitador);
                		Ipv4Address nodeAddress = InetSocketAddress::ConvertFrom(from).GetIpv4();

                		std::string publicKeyStr = parsed.substr(0,keyPos).c_str();
                		m_publicEcdhKeys[nodeAddress] = publicKeyStr;
                	}
                	default:
                		break;
                }
                totalReceived.erase(0, pos + delimitador.length());
			}
			m_bufferedData[from] = totalReceived;
			delete[] pcktInfo;
		} else if( InetSocketAddress::IsMatchingType(from)){
			std::cout << "At time "  << Simulator::Now ().GetSeconds ()
                       << "s node " << GetNode ()->GetId () << " received "
                       <<  pckt->GetSize () << " bytes from "
                       << InetSocketAddress::ConvertFrom(from).GetIpv4 ()
                       << " port " << InetSocketAddress::ConvertFrom (from).GetPort () <<  std::endl;
		}
	}
}

void NodeApp::StopApplication () {
	for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) {
		m_peersSockets[*i]->Close();
	}

	if (m_socket) {
		m_socket->Close();
		m_socket->SetRecvCallback(MakeNullCallback <void, Ptr<Socket> > ());
	}

	std::cout << "Node " << GetNode()->GetId() << " stopped" << std::endl;
}
void NodeApp::SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, rapidjson::Document &d, Ptr<Socket> outgoingSocket) {
	
	const uint8_t delimitador[] = "#";

	rapidjson::StringBuffer buffer;
	rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);

	d["message"].SetInt(responseMessage);
	d.Accept(writer);

	std::cout << "Node " << GetNode()->GetId() << "sent a " << getMessageName(responseMessage) << " message: " << buffer.GetString() << std::endl;

	outgoingSocket->Send(reinterpret_cast<const uint8_t*> (buffer.GetString()), buffer.GetSize(), 0);
	outgoingSocket->Send (delimitador, 1, 0);
}

void NodeApp::SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, std::string packet, Address &outgoingAddress) {
	const uint8_t delimitador[] = "#";

	rapidjson::StringBuffer buffer;
	rapidjson::Document doc;
	rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);

	doc.Parse(packet.c_str());
	doc["message"].SetInt(responseMessage);
	doc.Accept(writer);

	std::cout << "Node " << GetNode()->GetId() << "sent a " << getMessageName(responseMessage) << " message: " << buffer.GetString() << std::endl;
	Ipv4Address outgoingIpv4Address = InetSocketAddress::ConvertFrom(outgoingAddress).GetIpv4();
	std::map<Ipv4Address, Ptr<Socket>>::iterator it = m_peersSockets.find(outgoingIpv4Address);
	if(it == m_peersSockets.end()) {
		m_peersSockets[outgoingIpv4Address] = Socket::CreateSocket (GetNode(), TcpSocketFactory::GetTypeId());
		m_peersSockets[outgoingIpv4Address]->Connect (InetSocketAddress(outgoingIpv4Address, m_Port));
	}

	Ptr<Socket> outgoingSocket = m_peersSockets[outgoingIpv4Address];
	outgoingSocket->Send(reinterpret_cast<const uint8_t*> (buffer.GetString()), buffer.GetSize(), 0);
	outgoingSocket->Send (delimitador, 1, 0);
}

void NodeApp::SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, std::string packet, Ipv4Address &outgoingIpv4Address) {
	const uint8_t delimitador[] = "#";

	rapidjson::StringBuffer buffer;
	rapidjson::Document doc;
	rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);

	doc.Parse(packet.c_str());
	doc["message"].SetInt(responseMessage);
	doc.Accept(writer);

	std::cout << "Node " << GetNode()->GetId() << "sent a " << getMessageName(responseMessage) << " message: " << buffer.GetString() << std::endl;

	std::map<Ipv4Address, Ptr<Socket>>::iterator it = m_peersSockets.find(outgoingIpv4Address);
	if(it == m_peersSockets.end()) {
		m_peersSockets[outgoingIpv4Address] = Socket::CreateSocket (GetNode(), TcpSocketFactory::GetTypeId());
		m_peersSockets[outgoingIpv4Address]->Connect (InetSocketAddress(outgoingIpv4Address, m_Port));
	}
	
	Ptr<Socket> outgoingSocket = m_peersSockets[outgoingIpv4Address];
	outgoingSocket->Send(reinterpret_cast<const uint8_t*> (buffer.GetString()), buffer.GetSize(), 0);
	outgoingSocket->Send (delimitador, 1, 0);
}

void NodeApp::SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, std::string packet, Ptr<Socket> outgoingSocket) {
	const uint8_t delimitador[] = "#";

	rapidjson::StringBuffer buffer;
	rapidjson::Document doc;
	rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);

	doc.Parse(packet.c_str());
	doc["message"].SetInt(responseMessage);
	doc.Accept(writer);

	std::cout << "Node " << GetNode()->GetId() << "sent a " << getMessageName(responseMessage) << " message: " << buffer.GetString() << std::endl;

	outgoingSocket->Send(reinterpret_cast<const uint8_t*> (buffer.GetString()), buffer.GetSize(), 0);
	outgoingSocket->Send (delimitador, 1, 0);
}

void NodeApp::HandlePeerClose (Ptr<Socket> socket)
{
  //NS_LOG_FUNCTION (this << socket);
}

void NodeApp::HandlePeerError (Ptr<Socket> socket)
{
  //NS_LOG_FUNCTION (this << socket);
}

void NodeApp::HandleAccept (Ptr<Socket> s, const Address& from)
{
  //NS_LOG_FUNCTION (this << s << from);
  s->SetRecvCallback (MakeCallback (&NodeApp::HandleRead, this));
}



