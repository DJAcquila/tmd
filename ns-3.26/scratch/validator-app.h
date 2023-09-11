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

#include <algorithm>

#include "../rapidjson/document.h"
#include "../rapidjson/writer.h"
#include "../rapidjson/stringbuffer.h"

//#include "block.h"



using namespace ns3;

class ValidatorApp : public Application {

	
	//NS_OBJECT_ENSURE_REGISTERED (ValidatorApp);

protected:
	Ptr<Socket>									m_socket;
	Address 									m_local;
	TypeId 										m_tid;
	int 										m_numberOfPeers;
	double										m_meanBlockReceiveTime;
	double										m_previousBlockReceiveTime;
	double										m_meanBlockPropagationTime;
	double										m_meanBlockSize;
	Blockchain 									m_blockchain;
	Time 										m_invTimeoutMinutes;
	bool 										m_isCH;
	double 										m_downloadSpeed;
	double 										m_uploadSpeed;
	bool 										m_blockTorrent;
	uint32_t									m_chunkSize;
	bool 										m_spv;
	
	

	std::vector<Ipv4Address> 					m_peersAddresses;
	std::vector<int> 							m_peersId;
	std::map<Ipv4Address, double> 				m_peersDownloadSpeeds;
	std::map<Ipv4Address, double>				m_peersUploadSpeeds;          
	std::map<Ipv4Address, Ptr<Socket>>			m_peersSockets;
	std::map<std::string, std::vector<Address>>	m_queueInv;
	std::map<std::string, std::vector<Address>>	m_queueChunkPeers;
	std::map<std::string, std::vector<int>>		m_queueChunks;
	std::map<std::string, std::vector<int>>		m_receivedChunks;
	std::map<std::string, EventId>				m_invTimeouts;
	std::map<std::string, EventId>				m_chunkTimeouts;
	std::map<Address, std::string>				m_bufferedData;
	std::map<std::string, Block>				m_receivedNotValidated;
	std::map<std::string, Block>				m_onlyHeadersReceived;
	std::vector<double> 						m_sendBlockTimes;
	std::vector<double>							m_sendCompressedBlockTimes;
	std::vector<double>							m_receiveBlockTimes;
	std::vector<double>							m_receiveCompressedBlockTimes;


	const int 									m_blockchainPort;
	const int 									m_secondsPerMin;
	bool 										m_isValidator;
	const int 									m_countBytes;
	const int 									m_blockchainMessageHeader;
	const int 									m_inventorySizeBytes;
	const int 									m_getHeadersSizeBytes;
	int 										m_headersSizeBytes;
	const int 									m_blockHeadersSizeBytes;
	double 										m_averageTransactionSize;
	int 										m_transactionIndexSize;

	nodeStatistics								*m_nodeStats;
	enum ProtocolType 							m_protocolType;

	TracedCallback<Ptr<const Packet>, const Address &> m_rxTrace;



public:
	static TypeId GetTypeId(void) {
		static TypeId tid = TypeId ("ValidatorApp")
		.SetParent<Application> ()
		.SetGroupName("Applications")
		.AddConstructor<ValidatorApp> ()
		.AddAttribute ("Local",
				"The Address on which to Bind the rx socket.",
				AddressValue (),
				MakeAddressAccessor (&ValidatorApp::m_local),
				MakeAddressChecker ())
		.AddAttribute ("Protocol",
				"The type id of the protocol to use for the rx socket.",
				TypeIdValue (UdpSocketFactory::GetTypeId ()),
				MakeTypeIdAccessor (&ValidatorApp::m_tid),
				MakeTypeIdChecker ())
		.AddAttribute ("BlockTorrent",
				"Enable the BlockTorrent protocol",
				BooleanValue (false),
				MakeBooleanAccessor (&ValidatorApp::m_blockTorrent),
				MakeBooleanChecker ())
		.AddAttribute ("SPV",
				"Enable SPV Mechanism",
				BooleanValue (false),
				MakeBooleanAccessor (&ValidatorApp::m_spv),
				MakeBooleanChecker ())
		.AddAttribute ("InvTimeoutMinutes",
				"The timeout of inv messages in minutes",
				TimeValue (Minutes (20)),
				MakeTimeAccessor (&ValidatorApp::m_invTimeoutMinutes),
				MakeTimeChecker())
		.AddAttribute ("ChunkSize",
				"The fixed size of the block chunk",
				UintegerValue (100000),
				MakeUintegerAccessor (&ValidatorApp::m_chunkSize),
				MakeUintegerChecker<uint32_t> ())
		.AddTraceSource ("Rx",
				 "A packet has been received",
				 MakeTraceSourceAccessor (&ValidatorApp::m_rxTrace),
				 "ns3::Packet::AddressTracedCallback")
		;
		return tid;
	}


	ValidatorApp (void):  m_blockchainPort (8333), m_secondsPerMin(60), m_isValidator (true), m_countBytes (4), m_blockchainMessageHeader (90),
                          m_inventorySizeBytes (36), m_getHeadersSizeBytes (72), m_headersSizeBytes (81), m_blockHeadersSizeBytes (81),
                          m_averageTransactionSize (522.4), m_transactionIndexSize (2)
	{
	  m_socket = 0;
	  m_meanBlockReceiveTime = 0;
	  m_previousBlockReceiveTime = 0;
	  m_meanBlockPropagationTime = 0;
	  m_meanBlockSize = 0;
	  m_numberOfPeers = m_peersAddresses.size();

	}

	~ValidatorApp(void) {
	}

	Ptr<Socket> GetListeningSocket (void) const {
		return m_socket;
	}
	
	std::vector<Ipv4Address> GetPeersAddresses (void) const {
		return m_peersAddresses;
	}

	std::vector<int> GetPeersIds (void) const {
		return m_peersId;
	}

	void SetPeersAddresses (const std::vector<Ipv4Address> &peers) {
	  m_peersAddresses = peers;
	  m_numberOfPeers = m_peersAddresses.size();
	}

	void SetPeersIds (const std::vector<int> &peers) {
	  m_peersId = peers;
	  m_numberOfPeers = m_peersId.size();
	}

	void SetPeersDownloadSpeeds (const std::map<Ipv4Address, double> &peersDownloadSpeeds) {
		m_peersDownloadSpeeds = peersDownloadSpeeds;
	}

	void SetPeersUploadSpeeds (const std::map<Ipv4Address, double> &peersUploadSpeeds) {
		m_peersUploadSpeeds = peersUploadSpeeds;
	}

	void SetNodeInternetSpeeds (const nodeInternetSpeeds &internetSpeeds) {
		m_downloadSpeed = internetSpeeds.downloadSpeed * 1000000 / 8 ;
		m_uploadSpeed = internetSpeeds.uploadSpeed * 1000000 / 8 ;
	}
	
	void SetNodeStats (nodeStatistics *nodeStats){
		m_nodeStats = nodeStats;
	}

	void SetProtocolType (enum ProtocolType protocolType) {
		m_protocolType = protocolType;
	}

	void DoDispose (void) {
		m_socket = 0;
		Application::DoDispose ();
	}

	void HandlePeerClose (Ptr<Socket> socket) {
		//NS_LOG_FUNCTION (this << socket);
	}

	void HandlePeerError (Ptr<Socket> socket) {
		//NS_LOG_FUNCTION (this << socket);
	}

	void HandleAccept (Ptr<Socket> s, const Address& from) {
		//NS_LOG_FUNCTION (this << s << from);
		s->SetRecvCallback (MakeCallback (&ValidatorApp::HandleRead, this));
	}



	void StartApplication() {
		srand(time(NULL) + GetNode()->GetId());

		std::cout << "Node " << GetNode()->GetId() << ": download speed = " << m_downloadSpeed << " B/s" << std::endl;
		std::cout << "Node " << GetNode()->GetId() << ": upload speed = " << m_uploadSpeed << " B/s" << std::endl;
		std::cout << "Node " << GetNode()->GetId() << ": m_numberOfPeers = " << m_numberOfPeers << std::endl;
		std::cout << "Node " << GetNode()->GetId() << ": m_invTimeoutMinutes = " << m_invTimeoutMinutes.GetMinutes() << "mins" << std::endl;
		std::cout << "Node " << GetNode()->GetId() << ": m_protocolType = " << getProtocolType(m_protocolType) << std::endl;
		std::cout << "Node " << GetNode()->GetId() << ": m_blockTorrent = " << m_blockTorrent << std::endl;
		std::cout << "Node " << GetNode()->GetId() << ": m_chunkSize = " << m_chunkSize << " Bytes" << std::endl;

		std::cout << "Node " << GetNode()->GetId() << ": My peers are" << std::endl;


		for (auto it = m_peersAddresses.begin(); it != m_peersAddresses.end(); it++)
			std::cout << "\t" << *it << std::endl;

		if (!m_socket) {

			m_socket = Socket::CreateSocket (GetNode (), m_tid);
			m_socket->Bind (m_local);
			m_socket->Listen ();
			m_socket->ShutdownSend ();
			if (addressUtils::IsMulticast (m_local)) {

				Ptr<UdpSocket> udpSocket = DynamicCast<UdpSocket> (m_socket);
				if (udpSocket){

					// equivalent to setsockopt (MCAST_JOIN_GROUP)
					udpSocket->MulticastJoinGroup (0, m_local);
				}
				else {
					std::cout << "Error: joining multicast on a non-UDP socket" << std::endl;
				}
			}
		}
		m_socket->SetRecvCallback (MakeCallback (&ValidatorApp::HandleRead, this));

		m_socket->SetAcceptCallback (
			MakeNullCallback<bool, Ptr<Socket>, const Address &> (),
			MakeCallback (&ValidatorApp::HandleAccept, this)
		);

		m_socket->SetCloseCallbacks (
			MakeCallback (&ValidatorApp::HandlePeerClose, this),
			MakeCallback (&ValidatorApp::HandlePeerError, this)
		);

		std::cout << "Node" <<  GetNode()->GetId() << " before creating sockets " << std::endl;
		for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
		{
			m_peersSockets[*i] = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
			m_peersSockets[*i]->Connect (InetSocketAddress (*i, m_blockchainPort));
		}
		std::cout << "Node" <<  GetNode()->GetId() << " after creating sockets " << std::endl;
	}


	void StopApplication() {
		for (std::vector<Ipv4Address>::iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) //close the outgoing sockets
		{
			m_peersSockets[*i]->Close ();
		}


		if (m_socket)
		{
			m_socket->Close ();
			m_socket->SetRecvCallback (MakeNullCallback<void, Ptr<Socket> > ());
		}

		std::cout << "\n\nBLOCKCHAIN NODE " << GetNode ()->GetId () << ":" << std::endl;
		std::cout << "Current Top Block is:\n" << *(m_blockchain.GetTheTop()) << std::endl;
		std::cout << "Current Blockchain is:\n" << m_blockchain << std::endl;

		std::cout << "Mean Block Receive Time = " << m_meanBlockReceiveTime << " or "
		<< static_cast<int>(m_meanBlockReceiveTime) / m_secondsPerMin << "min and "
		<< m_meanBlockReceiveTime - static_cast<int>(m_meanBlockReceiveTime) / m_secondsPerMin * m_secondsPerMin << "s" << std::endl;
		std::cout << "Mean Block Propagation Time = " << m_meanBlockPropagationTime << "s" << std::endl;
		std::cout << "Mean Block Size = " << m_meanBlockSize << " Bytes" << std::endl;
		std::cout << "Total Blocks = " << m_blockchain.GetTotalBlocks() << std::endl;
		std::cout << "Stale Blocks = " << m_blockchain.GetNoObsoletBlocks() << " (" << 100. * m_blockchain.GetNoObsoletBlocks() / m_blockchain.GetTotalBlocks() << "%)" << std::endl;
		std::cout << "receivedButNotValidated size = " << m_receivedNotValidated.size() << std::endl;
		std::cout << "m_sendBlockTimes size = " << m_sendBlockTimes.size() << std::endl;
		std::cout << "m_receiveBlockTimes size = " << m_receiveBlockTimes.size() << std::endl;
		std::cout << "longest fork = " << m_blockchain.GetTheLongestFork() << std::endl;
		std::cout << "blocks in forks = " << m_blockchain.GetBlocksInFork() << std::endl;

	}

	void HandleRead(Ptr<Socket> socket) {
		Ptr<Packet> packet;
		Address from;

		//double newBlockReceiveTime = Simulator::Now().GetSeconds();

		while ((packet = socket->RecvFrom (from))) {
			if(packet->GetSize() == 0) {
				std::cout << "Empty packet" << std::endl;
				break;
			}
			if (InetSocketAddress::IsMatchingType(from)) {
				std::string delimiter = "#";
				std::string parsedPacket;
				size_t pos = 0;
				char* packetInfo = new char [packet->GetSize() + 1];
				std::ostringstream totalStream;

				packet->CopyData (reinterpret_cast <uint8_t*> (packetInfo), packet->GetSize());
				packetInfo[packet->GetSize()] = '\0';


				totalStream << m_bufferedData[from] << packetInfo;
				std::string totalReceivedData(totalStream.str());
				std::cout << "Node " << GetNode ()->GetId () << " Total Received Data: " << totalReceivedData << std::endl;

				while ((pos = totalReceivedData.find(delimiter)) != std::string::npos) {
					parsedPacket = totalReceivedData.substr(0, pos);
					std::cout << "Node " << GetNode ()->GetId () << " Parsed Packet: " << parsedPacket << std::endl;
					rapidjson::Document d;
					d.Parse(parsedPacket.c_str());
					// Objeto corrompido
					if(!d.IsObject())
					{
						std::cout << "The parsed packet is corrupted" << std::endl;
						totalReceivedData.erase(0, pos + delimiter.length());
						continue;
					}
					rapidjson::StringBuffer buffer;
					rapidjson::Writer <rapidjson::StringBuffer> writer(buffer);
					d.Accept(writer);

					std::cout 	<< "At time "  << Simulator::Now ().GetSeconds ()
								<< "s node " << GetNode ()->GetId () << " received "
								<<  packet->GetSize () << " bytes from "
								<< InetSocketAddress::ConvertFrom(from).GetIpv4 ()
								<< " port " << InetSocketAddress::ConvertFrom (from).GetPort ()
								<< " with info = " << buffer.GetString() << std::endl;
					switch (d["message"].GetInt()) {
						case INV:
						{							
							std::vector<std::string> requestBlocks;
							std::vector<std::string>::iterator block_it;
							m_nodeStats->invReceivedBytes += m_blockchainMessageHeader + m_countBytes +d["ivn"].Size()*m_inventorySizeBytes;

							for (int j = 0; j < (int) d["inv"].Size(); j++) {
								std::string 	invDelimiter = "/";
								std::string 	parsedInv = d["inv"][j].GetString();
								size_t 			invPos = parsedInv.find(invDelimiter);
								EventId 		timeout;

								int height = atoi(parsedInv.substr(0, invPos).c_str());
								int contentId = atoi(parsedInv.substr(invPos+1, parsedInv.size()).c_str());

								// Verifica se o bloco já foi recebido pelo nó ou se foi recebido mas não validado
								if (m_blockchain.HasBlock(height, contentId) || m_blockchain.IsOrphan(height, contentId) || ReceivedButNotValidated(parsedInv)) {
									std::cout 	<< "INV: Blockchain node " << GetNode ()->GetId ()
												<< " has already received the block with height = "
												<< height << " and contentId = " << contentId << std::endl;
								} else {
									std::cout << "INV: Blockchain node " << GetNode ()->GetId ()
				                              << " does not have the block with height = "
				                              << height << " and contentId = " << contentId << std::endl;

				                    // Verifica se o nó não requisitou o bloco ainda
				                    if (m_invTimeouts.find(parsedInv) == m_invTimeouts.end()) {
										std::cout 	<< "INV: Blockchain node " << GetNode ()->GetId () 
													<< " has not requested the block" << std::endl;
										requestBlocks.push_back(parsedInv);
										// Tenho que trocar o período de configuração INV!!!!!
										timeout = Simulator::Schedule(m_invTimeoutMinutes, &ValidatorApp::InvTimeoutExpired, this, parsedInv);
										m_invTimeouts[parsedInv] = timeout;
									} else {
										std::cout 	<< "INV: Blockchain node " << GetNode ()->GetId () 
													<< " has already requested the block" << std::endl;
									}
									m_queueInv[parsedInv].push_back(from);
								}	
							} // end-for

							if (!requestBlocks.empty()) {
								rapidjson::Value 	value;
								rapidjson::Value 	array(rapidjson::kArrayType);
								d.RemoveMember("inv");

								for (block_it = requestBlocks.begin(); block_it < requestBlocks.end(); block_it++) {
									value.SetString(block_it->c_str(), block_it->size(), d.GetAllocator());
									array.PushBack(value, d.GetAllocator());
								}
								d.AddMember("blocks", array, d.GetAllocator());

								SendMessage(INV, GET_HEADERS, d, from);
								SendMessage(INV, GET_DATA, d, from);
								// PAREI AQUI!!!!!!
							}
						} // end-case
					} //end-switch
				} // end-while
			} // end_if
		} // end-while
	} // end-HandleRead

	bool ReceivedButNotValidated(std::string blockHash) {
		if (m_receivedNotValidated.find(blockHash) != m_receivedNotValidated.end()) {
			return true;
		} else {
			return false;
		}
	}

	void InvTimeoutExpired (std::string blockHash) {
		std::string 	invDelimiter = "/";
		size_t			invPos = blockHash.find(invDelimiter);

		int height = atoi(blockHash.substr(0, invPos).c_str());
		int contentId = atoi(blockHash.substr(invPos+1, blockHash.size()).c_str());

		std::cout << "Node" << GetNode ()->GetId() << ": At time " << Simulator::Now ().GetSeconds()
		<< " the timeout for block " << blockHash << " expired" << std::endl; 

		m_nodeStats->blockTimeouts++;

		m_queueInv[blockHash].erase(m_queueInv[blockHash].begin());
		m_invTimeouts.erase(blockHash);

		if (!m_queueInv[blockHash].empty() && !m_blockchain.HasBlock(height, contentId) && !m_blockchain.IsOrphan(height, contentId) && !ReceivedButNotValidated(blockHash))
		{
			rapidjson::Document 	d;
			EventId               	timeout;
			rapidjson::Value 		value(INV);
			rapidjson::Value  		array(rapidjson::kArrayType);

			d.SetObject();
			d.AddMember("message", value, d.GetAllocator());

			value.SetString("block");
			d.AddMember("type", value, d.GetAllocator());


			value.SetString(blockHash.c_str(), blockHash.size(), d.GetAllocator());
			array.PushBack(value, d.GetAllocator());
			d.AddMember("blocks", array, d.GetAllocator());

			int index = rand() % m_queueInv[blockHash].size();
			Address temp = m_queueInv[blockHash][0];
			m_queueInv[blockHash][0] = m_queueInv[blockHash][index];
			m_queueInv[blockHash][index] = temp;

			SendMessage(INV, GET_HEADERS, d, *(m_queueInv[blockHash].begin()));
			SendMessage(INV, GET_DATA, d, *(m_queueInv[blockHash].begin()));


			timeout = Simulator::Schedule (m_invTimeoutMinutes, &ValidatorApp::InvTimeoutExpired, this, blockHash);
			m_invTimeouts[blockHash] = timeout;

			//d.AddMember("message", );
		} else {
			std::cout << "InvTimeoutExpired - else" << std::endl;
			m_queueInv.erase(blockHash);
		}		
	}

	void SendMessage(enum Messages receivedMessage, enum Messages responseMessage, rapidjson::Document &d, Ptr<Socket> outgoingSocket) {
		const uint8_t delimiter[] = "#";

		rapidjson::StringBuffer buffer;
		rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);

		d["message"].SetInt(responseMessage);

		d.Accept(writer);
		std::cout 	<< "Node " << GetNode ()->GetId () << " got a "
					<< getMessageName(receivedMessage) << " message"
					<< " and sent a " << getMessageName(responseMessage)
					<< " message: " << buffer.GetString() << std::endl;

		outgoingSocket->Send (reinterpret_cast<const uint8_t*>(buffer.GetString()), buffer.GetSize(), 0);
		outgoingSocket->Send (delimiter, 1, 0);

		switch (d["message"].GetInt()) {
			case INV:
			{
				m_nodeStats->invSentBytes += m_blockchainMessageHeader + m_countBytes + d["inv"].Size()*m_inventorySizeBytes;
				break;
			}
			case GET_HEADERS:
			{
				m_nodeStats->getHeadersSentBytes += m_blockchainMessageHeader + m_getHeadersSizeBytes;
				break;
			}
			case HEADERS:
			{
				m_nodeStats->headersSentBytes += m_blockchainMessageHeader + m_countBytes + d["blocks"].Size()*m_headersSizeBytes;
				break;
			}
			case BLOCK:
			{
				for(int k = 0; k < (int) d["blocks"].Size(); k++)
				m_nodeStats->blockSentBytes += d["blocks"][k]["size"].GetInt();
				m_nodeStats->blockSentBytes += m_blockchainMessageHeader;
				break;
			}
			case GET_DATA:
			{
				m_nodeStats->getDataSentBytes += m_blockchainMessageHeader + m_countBytes + d["blocks"].Size()*m_inventorySizeBytes;
				break;
			}
		} // end switch
	}

	void SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, rapidjson::Document &d, Address &outgoingAddress) {
		const uint8_t delimiter[] = "#";

		rapidjson::StringBuffer buffer;
		rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);

		d["message"].SetInt(responseMessage);

		d.Accept(writer);
		std::cout 	<< "Node " << GetNode ()->GetId () << " got a "
					<< getMessageName(receivedMessage) << " message"
					<< " and sent a " << getMessageName(responseMessage)
					<< " message: " << buffer.GetString() << std::endl;

		Ipv4Address outgoingIpv4Address = InetSocketAddress::ConvertFrom(outgoingAddress).GetIpv4 ();
		std::map<Ipv4Address, Ptr<Socket>>::iterator it = m_peersSockets.find(outgoingIpv4Address);

		if (it == m_peersSockets.end()) //Create the socket if it doesn't exist
		{
			m_peersSockets[outgoingIpv4Address] = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
			m_peersSockets[outgoingIpv4Address]->Connect (InetSocketAddress (outgoingIpv4Address, m_blockchainPort));
		}
		m_peersSockets[outgoingIpv4Address]->Send (reinterpret_cast<const uint8_t*>(buffer.GetString()), buffer.GetSize(), 0);
		m_peersSockets[outgoingIpv4Address]->Send (delimiter, 1, 0);

		switch (d["message"].GetInt()) {
			case INV:
			{
				m_nodeStats->invSentBytes += m_blockchainMessageHeader + m_countBytes + d["inv"].Size()*m_inventorySizeBytes;
				break;
			}
			case GET_HEADERS:
			{
				m_nodeStats->getHeadersSentBytes += m_blockchainMessageHeader + m_getHeadersSizeBytes;
				break;
			}
			case HEADERS:
			{
				m_nodeStats->headersSentBytes += m_blockchainMessageHeader + m_countBytes + d["blocks"].Size()*m_headersSizeBytes;
				break;
			}
			case BLOCK:
			{
				for(int k = 0; k < (int) d["blocks"].Size(); k++)
				m_nodeStats->blockSentBytes += d["blocks"][k]["size"].GetInt();
				m_nodeStats->blockSentBytes += m_blockchainMessageHeader;
				break;
			}
			case GET_DATA:
			{
				m_nodeStats->getDataSentBytes += m_blockchainMessageHeader + m_countBytes + d["blocks"].Size()*m_inventorySizeBytes;
				break;
			}
		} // end switch
	}

	void SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, std::string packet, Address &outgoingAddress) {
		const uint8_t delimiter[] = "#";
		rapidjson::Document d;

		rapidjson::StringBuffer buffer;
		rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
		d.Parse(packet.c_str());
		d["message"].SetInt(responseMessage);

		d.Accept(writer);
		std::cout 	<< "Node " << GetNode ()->GetId () << " got a "
					<< getMessageName(receivedMessage) << " message"
					<< " and sent a " << getMessageName(responseMessage)
					<< " message: " << buffer.GetString() << std::endl;

		Ipv4Address outgoingIpv4Address = InetSocketAddress::ConvertFrom(outgoingAddress).GetIpv4 ();
		std::map<Ipv4Address, Ptr<Socket>>::iterator it = m_peersSockets.find(outgoingIpv4Address);

		if (it == m_peersSockets.end()) //Create the socket if it doesn't exist
		{
			m_peersSockets[outgoingIpv4Address] = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
			m_peersSockets[outgoingIpv4Address]->Connect (InetSocketAddress (outgoingIpv4Address, m_blockchainPort));
		}
		m_peersSockets[outgoingIpv4Address]->Send (reinterpret_cast<const uint8_t*>(buffer.GetString()), buffer.GetSize(), 0);
		m_peersSockets[outgoingIpv4Address]->Send (delimiter, 1, 0);

		switch (d["message"].GetInt()) {
			case INV:
			{
				m_nodeStats->invSentBytes += m_blockchainMessageHeader + m_countBytes + d["inv"].Size()*m_inventorySizeBytes;
				break;
			}
			case GET_HEADERS:
			{
				m_nodeStats->getHeadersSentBytes += m_blockchainMessageHeader + m_getHeadersSizeBytes;
				break;
			}
			case HEADERS:
			{
				m_nodeStats->headersSentBytes += m_blockchainMessageHeader + m_countBytes + d["blocks"].Size()*m_headersSizeBytes;
				break;
			}
			case BLOCK:
			{
				for(int k = 0; k < (int) d["blocks"].Size(); k++)
				m_nodeStats->blockSentBytes += d["blocks"][k]["size"].GetInt();
				m_nodeStats->blockSentBytes += m_blockchainMessageHeader;
				break;
			}
			case GET_DATA:
			{
				m_nodeStats->getDataSentBytes += m_blockchainMessageHeader + m_countBytes + d["blocks"].Size()*m_inventorySizeBytes;
				break;
			}
		} // end switch
	}

};