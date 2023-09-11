/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
* Copyright (c) 2015
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License version 2 as
* published by the Free Software Foundation;
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
* Author: Fausto Moraes <fausto@inf.ufg.br>
*/
#include "ns3/lte-helper.h"
#include "ns3/epc-helper.h"
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/internet-module.h"
#include "ns3/mobility-module.h"
#include "ns3/lte-module.h"
#include "ns3/applications-module.h"
#include "ns3/packet-sink.h"
#include "ns3/point-to-point-helper.h"
#include "ns3/config-store.h"
#include "ns3/propagation-loss-model.h"
#include "ns3/wifi-module.h"
#include "ns3/flow-monitor-module.h"
#include "ns3/olsr-helper.h"
#include "ns3/energy-module.h"
#include "src/olsr/model/olsr-routing-protocol.h"
#include "ns3/stats-module.h"
#include "trust.h"
//#include "messages.h"
#include "ns3/v4ping-helper.h"
#include "ns3/csma-module.h"

#include "scratch/graph.h"
#include "scratch/keys.h"
#include "scratch/ecc.h"
#include "scratch/evaluate.h"
// #include "scratch/trust.h"
//#include "scratch/myTag.h"
//#include "scratch/node-app.h"
//#include "scratch/validator-app.h"
#include "scratch/PBFT_consensus.h"

#include <limits.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <bits/stdc++.h>
#include <iomanip>
#include <math.h>
#include <string>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/tokenizer.hpp>
#include <vector>
#include <tuple>
#include <random>
#include <vector>
#include <queue> 

#include "ns3/lte-enb-rrc.h"


#include <map>

//#include "ns3/gtk-config-store.h"
using namespace ns3;
using namespace CryptoPP;

struct greater
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};
struct CompareSecond
{
    bool operator()(const std::pair<int, double>& left, const std::pair<int, double>& right) const
    {
        return left.second < right.second;
    }
};

NS_LOG_COMPONENT_DEFINE("EpcFirstExample");

double
lowerIncompleteGamma(double x, void *params);

double
generalisedPareto(double mu, double sigma, double E, double simTime);

double
Weibull(double lambda, double k);

double
closeness(double k, double teta);

void PrintRoutingTable(Ptr<Node> node);

//RoutingTableToMatrix (NodeContainer ueNodes, int numberOfNodes, Graph *g);
void RoutingTableToMatrix(NodeContainer ueNodes, int numberOfNodes, Graph *g, Contact *contact, ContactDistribution *contDist);
void RoutingTableToMatrixTwoHop(NodeContainer ueNodes, int numberOfNodes, TwoHopGraph *tg, TwoHopContact *tcontact, TwoHopContactDistribution *tcontDist);
void showStatistic(Content *cont);
//showStatistic (Content *cont, ContactDistribution *contDist);

void testMatrix(Graph *g);

void searchContent(Graph *sg, int node, int contentNumber,  Contact *contact);

void Cache(uint16_t nextNode, uint16_t nextContent, Contact *contact, Graph *sg);

//void
//Preload(int n, int zipfv, Content *cont, Graph *g);

//searchContent (NodeContainer ueNodes, Ptr<Node> remoteHost, Ipv4InterfaceContainer ueWiFiIface, Ipv4InterfaceContainer ueLteIface, double simTime, int node, int contentNumber, Content *cont, Graph *g);

void nodeSendContent(NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface, Content *cont, int node, Graph *g);

void internetSendContent(unsigned int node, NodeContainer ueNodes, Ptr<Node> remoteHost, Ipv4InterfaceContainer ueLteIfaces, Content *cont);

void sendKeysAndMetadata(NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface, int nodeRcv, int contentNumber);

void exchangeKeys(int nodeSrc, int nodeRcv);
void clusterKeyShare(Cluster* cluster, NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface);

void addContent(Graph *sg, int node, int contentNumber, Content *cont);

void verifyMetadata(int node, int contentNumber, Content *cont);

/// Trace function for remaining energy at node.
void RemainingEnergy(double oldValue, double remainingEnergy)
{
	if (Simulator::Now().GetSeconds() > 1000.05)
		std::cout << Simulator::Now().GetSeconds()
				  << "s Current remaining energy = " << remainingEnergy << "J" << std::endl;
}

/// Trace function for total energy consumption at node.
void TotalEnergy(double oldValue, double totalEnergy)
{
	if (Simulator::Now().GetSeconds() > 1000.05)
		std::cout << Simulator::Now().GetSeconds()
				  << "s Total energy consumed by radio on node = " << totalEnergy << "J" << std::endl;
}
static void SinkRx (Ptr<const Packet> p, const Address &ad)
{
  std::cout << *p << std::endl;
}
//void SendPublicKey(NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface, Content *cont, int node, Graph *g)
void GenerateECDSAKeys(int idx, Key *keyUtil);
void GenerateKeys(int numberOfNodes);

void insertDegMap(std::map<double, std::vector<int>>* m, double id, int v);
void printDegMap(std::map<double, std::vector<int>> m, double id);


void nodePosition(NodeContainer ueNodes, int numberOfNodes);
std::tuple<double, double> ueNodePosition(NodeContainer ueNodes, int node);
std::tuple<double, double> degMapPos(std::map<double, std::vector<int>> m, NodeContainer ueNodes, double id, int node);

void clusterHeadElection(NodeContainer ueNodes, Contact *contact/*, Cluster* cluster*/, Ipv4InterfaceContainer ueWiFiIface, TwoHopContact *tcontact);
void clusterSetUp(NodeContainer ueNodes/*, Cluster* cluster*/, Ipv4InterfaceContainer ueWiFiIface, Contact* g, TwoHopContact *tcontact);
void clusterSetUpPhase2(NodeContainer ueNodes/*, Cluster* cluster*/, Ipv4InterfaceContainer ueWiFiIface, Contact* g,  TwoHopContact *tcontact);
void clusterHeadManage(NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface/*, Cluster* cluster*/, Contact* contact);
void verifyConnectivity(NodeContainer ueNodes, Contact *contact, Ipv4InterfaceContainer ueWiFiIface, TwoHopContact *tcontact);
void sendCLusterNotification(Ptr<Socket> socket, NodeContainer ueNodes, int src, int dest);

void testeClusterElection(Contact *contact, TwoHopContact *tcontact);
int testeMaisProximo(int node, Contact* contact, TwoHopContact *tcontact);
void UpdateQualityLink(int sourceVertice, int destinationVertice, double snr, double sinr, double per, bool received);
void updateNetStats(FlowMonitorHelper* flowmon, Ptr<FlowMonitor> monitor, Ipv4InterfaceContainer ueLteIface);
void receiveClusterNotification(Ptr<Socket> socket);
void updateStatics();
void ReceivePacket (Ptr<Socket> socket)
{
	while (socket->Recv ())
	{
		NS_LOG_UNCOND ("Received one packet!");
	}
}

static void GenerateTraffic (Ptr<Socket> socket, uint32_t pktSize,
                      uint32_t pktCount, Time pktInterval )
{
	if (pktCount > 0)
	{
		socket->Send (Create<Packet> (pktSize));
		Simulator::Schedule (pktInterval, &GenerateTraffic,
		                    socket, pktSize,pktCount - 1, pktInterval);
	}
	else
	{
		socket->Close ();
	}
}
 
void addNodeMap(Ptr<Node> node, int id);
int getNode (Ptr<Node> node);
uint16_t numberOfNodes = 15;

double probV = 0.32;//0.32
int messageLen = 0;

/*********************************** Nós maliciosos ********************************/
double proportion = 0.30;
/***********************************************************************************/


std::map<Ptr<Node>,int> nodeMap;

NodeContainer ueNodes;
Cluster cluster(numberOfNodes);
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> dis(1, INT_MAX);

void ReceivePacketSink (Ptr<const Packet> packet, const Address &address){
	SnrTag tag;
	//std::cout << "SNR: ";
	if (packet->PeekPacketTag(tag) )
	{
		std::cout << tag.Get() << std::endl;
		NS_LOG_DEBUG("Received Packet with SNR ->" << tag.Get());
	}
	
}

RESULTS resultsInstance = {0,0,0.0,0.0,0.0,0,0,0,0,0,0,0,0,0,0, 0};

void cascadeBlockSec(int argc, char *argv[], std::string arg_traceFile)
{

	Time::SetResolution(Time::NS);
	//LogComponentEnable ("BulkSendApplication", LOG_LEVEL_INFO);
	//LogComponentEnable ("PacketSink", LOG_LEVEL_INFO);
	//LogComponentEnable ("OlsrRoutingProtocol", LOG_LEVEL_INFO);

	// DEFAULT VALUES
	
	// Diminuí a quantidade de conteúdos de 100 para 50
	uint16_t contentPopulation = 100;

	double simTime = 100.1;
	double distance = 40.0;
	uint32_t seed = 5;
	double alpha = 10;
	unsigned int cacheSize = 10;
	std::string traceFile = arg_traceFile;
	std::string socialGraphFile = "scratch/srsn.csv";

	BasicConfig* basicConfigInstance = BasicConfig::getInstance();
	basicConfigInstance->printScenarioInfos();

	//double interPacketInterval = 10;

	// Command line arguments
	CommandLine cmd;
	cmd.AddValue("numberOfNodes", "Number of UE nodes", numberOfNodes);
	cmd.AddValue("contentPopulation", "Number of content available", contentPopulation);
	cmd.AddValue("simTime", "Total duration of the simulation [s])", simTime);
	cmd.AddValue("distance", "Distance between nodes [m]", distance);
	cmd.AddValue("seed", "Random number Distribuiton seed", seed);
	cmd.AddValue("traceFile", "Trace File", traceFile);
	cmd.AddValue("socialGraphFile", "Social Graph", socialGraphFile);
	cmd.AddValue("alpha", "Alpha", alpha);
	cmd.AddValue("cacheSize", "Cache Size", cacheSize);
	Config::SetDefault ("ns3::LteEnbRrc::SrsPeriodicity", UintegerValue (80));

	//cmd.AddValue("interPacketInterval", "Inter packet interval [ms])", interPacketInterval);
	cmd.Parse(argc, argv);

	// Create node containers
	
	NodeContainer enbNodes;
	ueNodes.Create(numberOfNodes);
	enbNodes.Create(1);

	// Create LTE
	Ptr<LteHelper> lteHelper = CreateObject<LteHelper>();
	Ptr<PointToPointEpcHelper> epcHelper = CreateObject<PointToPointEpcHelper>();
	lteHelper->SetEpcHelper(epcHelper);

	ConfigStore inputConfig;
	inputConfig.ConfigureDefaults();

	// Parse again so you can override default values from the command line
	cmd.Parse(argc, argv);

	Ptr<Node> pgw = epcHelper->GetPgwNode();

	/// Internet

	// Create a single RemoteHost
	NodeContainer remoteHostContainer;
	remoteHostContainer.Create(1);
	Ptr<Node> remoteHost = remoteHostContainer.Get(0);
	InternetStackHelper internet;
	internet.Install(remoteHostContainer);

	// Create the Internet
	PointToPointHelper p2ph;
	p2ph.SetDeviceAttribute("DataRate", DataRateValue(DataRate("100Gb/s")));
	p2ph.SetDeviceAttribute("Mtu", UintegerValue(1500));
	p2ph.SetChannelAttribute("Delay", TimeValue(Seconds(0.010)));
	NetDeviceContainer internetDevices = p2ph.Install(pgw, remoteHost);
	Ipv4AddressHelper ipv4h;
	ipv4h.SetBase("1.0.0.0", "255.0.0.0");
	Ipv4InterfaceContainer internetIpIfaces = ipv4h.Assign(internetDevices);
	// interface 0 is localhost, 1 is the p2p device
	Ipv4Address remoteHostAddr = internetIpIfaces.GetAddress(0);

	Ipv4StaticRoutingHelper ipv4RoutingHelper;
	Ptr<Ipv4StaticRouting> remoteHostStaticRouting = ipv4RoutingHelper.GetStaticRouting(remoteHost->GetObject<Ipv4>());
	remoteHostStaticRouting->AddNetworkRouteTo(Ipv4Address("7.0.0.0"), Ipv4Mask("255.0.0.0"), Ipv4Address(remoteHostAddr), 1);

	/// Internet fim

	/// Mobility

	// ENB Mobility
	MobilityHelper enbMobility;
	Ptr<ListPositionAllocator> positionAlloc = CreateObject<ListPositionAllocator>();
	positionAlloc->Add(Vector((500.0), (500.0), 0.0));
	enbMobility.SetPositionAllocator(positionAlloc);
	//enbMobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
	enbMobility.Install(enbNodes);

	Ns2MobilityHelper mobility = Ns2MobilityHelper(traceFile);
	mobility.Install();

	/// Mobility FIM

	// Install LTE Devices to the nodes
	NetDeviceContainer enbLteDevs = lteHelper->InstallEnbDevice(enbNodes);
	NetDeviceContainer ueLteDevs = lteHelper->InstallUeDevice(ueNodes);

	// Create WiFi
	YansWifiChannelHelper wifiChannel = YansWifiChannelHelper::Default();
	//std::string phyMode ("DsssRate11Mbps");
	std::string phyMode("ErpOfdmRate24Mbps");

	// range	wifiChannel.SetPropagationDelay("ns3::ConstantSpeedPropagationDelayModel");
	// range	wifiChannel.AddPropagationLoss("ns3::RangePropagationLossModel", "MaxRange", DoubleValue(70.0));
	Ptr<YansWifiChannel> channel = wifiChannel.Create();
	LogDistancePropagationLossModel loss;
	loss.SetPathLossExponent(2.7);
	loss.SetReference(1, 46.6777);
	channel->SetPropagationLossModel(loss.GetObject<LogDistancePropagationLossModel>());
	YansWifiPhyHelper wifiPhy = YansWifiPhyHelper::Default();
	wifiPhy.SetChannel(channel);
	wifiPhy.Set("EnergyDetectionThreshold", DoubleValue(-96.0)); //Default val is -96dBm
	wifiPhy.Set("TxPowerStart", DoubleValue(40)); //!!!!!!
	wifiPhy.Set("TxPowerEnd", DoubleValue(40)); // !!!!!!
	//wifiPhy.Set ("CcaMode1Threshold", DoubleValue (-96 * 0.125892541));

	/// Wifi:
	WifiHelper wifi = WifiHelper::Default();
	// Set standard to 802.11g
	wifi.SetStandard(WIFI_PHY_STANDARD_80211g);
	// Use constant rates for data and control transmissions
	wifi.SetRemoteStationManager("ns3::ConstantRateWifiManager",
								 "MaxSsrc", UintegerValue(0),
								 "RtsCtsThreshold", UintegerValue(10000),
								 "FragmentationThreshold", StringValue("10000"),
								 "DataMode", StringValue(phyMode),
								 "ControlMode", StringValue(phyMode),
								 "NonUnicastMode", StringValue(phyMode));

	NqosWifiMacHelper wifiMac;
	wifiMac = NqosWifiMacHelper::Default();

	// Instal WiFi devices to the nodes
	NetDeviceContainer ueWiFiDevs = wifi.Install(wifiPhy, wifiMac, ueNodes);

	/// Energy Model
	// energy source
	BasicEnergySourceHelper basicSourceHelper;
	// configure energy source
	basicSourceHelper.Set("BasicEnergySourceInitialEnergyJ", DoubleValue(2000.0));
	// install source
	EnergySourceContainer sources = basicSourceHelper.Install(ueNodes);
	// device energy model
	WifiRadioEnergyModelHelper radioEnergyHelper;
	// configure radio energy model
	radioEnergyHelper.SetTxCurrentModel("ns3::LinearWifiTxCurrentModel");
	//radioEnergyHelper.Set ("TxCurrentA", DoubleValue (0.380));
	//radioEnergyHelper.Set ("RxCurrentA", DoubleValue (0.313));
	// install device model
	DeviceEnergyModelContainer deviceModels = radioEnergyHelper.Install(ueWiFiDevs, sources);

	/// connect trace sources
	// all traces are connected to node 1 (Destination)
	// energy source
	for (int j = 0; j < numberOfNodes; j++)
	{
		Ptr<BasicEnergySource> basicSourcePtr = DynamicCast<BasicEnergySource>(sources.Get(j));
		basicSourcePtr->TraceConnectWithoutContext("RemainingEnergy", MakeCallback(&RemainingEnergy));
		// device energy model
		Ptr<DeviceEnergyModel> basicRadioModelPtr =
			basicSourcePtr->FindDeviceEnergyModels("ns3::WifiRadioEnergyModel").Get(0);
		NS_ASSERT(basicRadioModelPtr != 0);
		basicRadioModelPtr->TraceConnectWithoutContext("TotalEnergyConsumption", MakeCallback(&TotalEnergy));
	}

	/// Create stack protocols:
	OlsrHelper olsr;
	Ipv4ListRoutingHelper list;
	list.Add(ipv4RoutingHelper, 0);
	list.Add(olsr, 10);
	//InternetStackHelper stack;
	internet.SetRoutingHelper(list);
	internet.Install(ueNodes);
	//internet.Install (meshNodes);

	//Ipv4InterfaceContainer ueLteIface;
	//ueLteIface = epcHelper->AssignUeIpv4Address (NetDeviceContainer (ueLteDevs));

	// Set IP address for wifi devs
	Ipv4AddressHelper ipv4;
	NS_LOG_INFO("Assign IP Addresses to the WiFi Devices.");
	ipv4.SetBase("10.1.1.0", "255.255.255.0");
	Ipv4InterfaceContainer ueWiFiIface = ipv4.Assign(ueWiFiDevs);

	// Install the IP stack on the UEs
	//internet.Install (ueNodes);
	Ipv4InterfaceContainer ueLteIface;
	ueLteIface = epcHelper->AssignUeIpv4Address(NetDeviceContainer(ueLteDevs));

	// Assign IP address to UEs
	for (uint32_t u = 0; u < ueNodes.GetN(); ++u)
	{
		Ptr<Node> ueNode = ueNodes.Get(u);
		// Set the default gateway for the UE
		Ptr<Ipv4StaticRouting> ueStaticRouting = ipv4RoutingHelper.GetStaticRouting(ueNode->GetObject<Ipv4>());
		ueStaticRouting->SetDefaultRoute(epcHelper->GetUeDefaultGatewayAddress(), 2);
	}

	// Attach one UE per eNodeB
	for (uint16_t i = 0; i < numberOfNodes; i++)
	{
		lteHelper->Attach(ueLteDevs.Get(i), enbLteDevs.Get(0));
		// side effect: the default EPS bearer will be activated
	}

	for (int i = 0; i < numberOfNodes; ++i)
	{
        addNodeMap(ueNodes.Get(i), i);
    }

	//lteHelper->EnableTraces ();
	// Uncomment to enable PCAP tracing
	p2ph.EnablePcapAll("lena-epc-first");

	// Social Network Graph

	//std::string data("srsn.csv");
	std::ifstream in(socialGraphFile.c_str());
	if (!in.is_open())
		std::cout <<  "Social data file do not exist!" << std::endl;
		//return 1;

	typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;

	std::vector<std::string> vec;
	std::string line;
	Graph socialGraph(numberOfNodes);

	while (std::getline(in, line))
	{

		Tokenizer tok(line);
		vec.assign(tok.begin(), tok.end());

		if (vec.size() < 2)
			continue;

		std::vector<int> vecNum;
		for (int i = 0; i < 2; i++)
		{
			int num = std::atoi(vec.at(i).c_str());
			vecNum.push_back(num);
		}
		socialGraph.addEdge((vecNum.at(0) - 1), (vecNum.at(1) - 1));

		//std::copy(vec.begin(), vec.end(),
		//std::ostream_iterator<std::string>(std::cout, "|"));
		//
		//std::cout << "\n----------------------" << std::endl;
	}
	in.close();
	// Número de CHs

	
	//int numberCH = (int) ceil(0.1*numberOfNodes);
	

	// Adjacency check
	//
	Graph graph(numberOfNodes);
	TwoHopGraph tGraph(numberOfNodes);

	Contact contact(numberOfNodes);
	TwoHopContact tcontact(numberOfNodes);
	ContactDistribution contDist(numberOfNodes);
	TwoHopContactDistribution tcontDist(numberOfNodes);

/*RoutingTableToMatrix(NodeContainer ueNodes, int numberOfNodes, Graph *g, TwoHopGraph *tg,  
	Contact *contact, TwoHopContact *tcontact, ContactDistribution *contDist, TwoHopContactDistribution *tcontDist)*/
	for (int update = 2; update < simTime; update += 2)
	{
		Simulator::Schedule(Seconds(update), &RoutingTableToMatrix, ueNodes, numberOfNodes, &graph, &contact, &contDist);
		Simulator::Schedule(Seconds(update), &RoutingTableToMatrixTwoHop, ueNodes, numberOfNodes, &tGraph, &tcontact, &tcontDist);
		// show adjacency Matrix	Simulator::Schedule (Seconds (2.5), &testMatrix, &socialGraph);
	}

	/*
	* Pre Simulation content selection
	*/
	auto blockchainConn = BlockchainConnection::getInstance(numberOfNodes);
	blockchainConn->startInstance();
	auto cont = Content::getInstance(numberOfNodes, contentPopulation, cacheSize);
	//Content cont(numberOfNodes, contentPopulation, cacheSize);
	View fv(numberOfNodes, contentPopulation);
	View fvCopy(numberOfNodes, contentPopulation);
	double shareRate[numberOfNodes];
	memset(shareRate, 0, numberOfNodes*sizeof(uint16_t));
	double shrT[numberOfNodes];
	memset(shrT, 0, numberOfNodes*sizeof(uint16_t));
	double gp = 0.0;
	double drift = 5.0;
	double sTime = drift;
	double simTimeR = simTime * (alpha / 5); // LIMITE PARA O SORTEIO DO TEMPO DE ACESSO
	srand(seed);
	SeedManager::SetSeed(seed);

	//uint32_t N = 100;
	//int neighContent;

	// Zhang ZIPf
	Ptr<ZipfRandomVariable> x = CreateObject<ZipfRandomVariable>();
	x->SetAttribute("N", IntegerValue(contentPopulation));
	x->SetAttribute("Alpha", DoubleValue(1.9));

	// UNIFORM
	Ptr<UniformRandomVariable> shrt = CreateObject<UniformRandomVariable>();
	shrt->SetAttribute("Min", DoubleValue(0.0));
	shrt->SetAttribute("Max", DoubleValue(1.0));

	gsl_rng *r_bn = gsl_rng_alloc(gsl_rng_rand);
	gsl_rng_set(r_bn, seed);

	// Content Schedule

	double gp_tmp;

	for (uint16_t n = 0; n < numberOfNodes; n++)
	{

		int zeros = cont->zeros();

		shareRate[n] = generalisedPareto(-0.227, 0.305, -0.048, simTimeR);
		shrT[n] = shrt->GetValue();
		for (uint16_t c = 0; c < contentPopulation; c++)
		{

			if (cont->getViews(c) == 0)
			{
				int bin = gsl_ran_binomial(r_bn, (alpha / ((n + 1) * zeros)), 1);
				if (bin == 1)
				{
					gp = generalisedPareto(-2654, 6315, -0.669, simTimeR);
					gp_tmp = gp;
					fv.addView(n, c, gp + drift);
					fvCopy.addView(n, c, gp + drift);

					cont->addViews(c);
					if (shareRate[n] > shrT[n])
					{
						for (uint16_t n_shr = 0; n_shr < numberOfNodes; n_shr++)
						{
							if ((((rand() % 1000) / 1000) < probV) && (socialGraph.isEdge(n, n_shr)))
							{
								gp = generalisedPareto(-2654, 6315, -0.669, simTimeR);
								fv.addView(n_shr, c, gp_tmp + gp);
								fvCopy.addView(n_shr, c, gp_tmp + gp);
								cont->addViews(c);
							}
						}
					}
				}
			}
			else
			{
				int pv = rand() % (n);
				if (pv <= cont->getViews(c))
				{
					gp = generalisedPareto(-2654, 6315, -0.669, simTimeR);
					gp_tmp = gp;
					fv.addView(n, c, gp + drift);
					fvCopy.addView(n, c, gp + drift);
					cont->addViews(c);
					if (shareRate[n] > shrT[n])
					{
						for (uint16_t n_shr = 0; n_shr < numberOfNodes; n_shr++)
						{
							if ((((rand() % 1000) / 1000) < probV) && (socialGraph.isEdge(n, n_shr)))
							{
								gp = generalisedPareto(-2654, 6315, -0.669, simTimeR);
								fv.addView(n_shr, c, gp_tmp + gp);
								fvCopy.addView(n_shr, c, gp_tmp + gp);
								cont->addViews(c);
							}
						}
					}
				}
			}
		}
	}

	//// PRELOAD
	//
	//for (uint16_t plc = 0; plc < floor(contentPopulation*0.5) ; plc++){
	//  int randomN = (rand () % 27) ;
	//  uint32_t v = x->GetInteger ();
	//  int zipfv = (short) v;
	//
	//  Simulator::Schedule (Seconds (4.5), &Preload, randomN, zipfv, &cont, &graph);
	//
	//}

	// Aplication Schedule
	//
	//for (uint32_t i = 0; i < flows.size (); i++) {

	Simulator::Schedule (Seconds (0.0), &GenerateKeys, numberOfNodes);
	std::cout << "Salvando os tempos de scheduling" << std::endl;
	
	std::vector<double> vStime;
	double timeScheduling = drift;
	while (timeScheduling < simTime) {
		if (fvCopy.nextContent > contentPopulation)
			fvCopy.firstView();
		else
			fvCopy.nextView(fvCopy.nextNode, fvCopy.nextContent);

		timeScheduling = fvCopy.getView(fvCopy.nextNode, fvCopy.nextContent);
		std::cout << timeScheduling << " ";
		vStime.push_back(timeScheduling);
	}	
	std::cout << std::endl;
   
	Trust* trustInstance = Trust::getInstance(numberOfNodes);
	trustInstance->setNumberOfNodes((int)numberOfNodes);
	if(basicConfigInstance->is_completeScenario()){
		trustInstance->setIsComplete(true);
	} else {
		trustInstance->setIsComplete(false);
	}

	MaliciousNodes* maliciousInstance = MaliciousNodes::getInstance(proportion, numberOfNodes);
	maliciousInstance->getQtd();
	TrustStatistics* trustStaticsInstance = TrustStatistics::getInstance(numberOfNodes);
	NetStats::getInstance(numberOfNodes);
	trustStaticsInstance->setNodeToAnalyze(maliciousInstance->selectANode());
	//trustValues.printValues();

	while (sTime < simTime)
	{

		if (fv.nextContent > contentPopulation)
			fv.firstView();
		else
			fv.nextView(fv.nextNode, fv.nextContent);

		sTime = fv.getView(fv.nextNode, fv.nextContent);



		// Schedule content search on the neighbors cache
		Simulator::Schedule(Seconds(sTime), &searchContent, &socialGraph, fv.nextNode, fv.nextContent, &contact);
				
		//Send Keys trhough neighbors D2D communication

		// Send content through the BS
		Simulator::Schedule(Seconds(sTime), &internetSendContent, fv.nextNode, ueNodes, remoteHost, ueLteIface, cont);

		// Send content through neighbors D2D communication
		Simulator::Schedule(Seconds(sTime), &nodeSendContent, ueNodes, ueWiFiIface, cont, fv.nextNode, &graph);

		//sTime+= incTime;

		// ProSoCaD CACHE - Share probability
		if (shareRate[fv.nextNode] > shrT[fv.nextNode])
		{
			Simulator::Schedule(Seconds(sTime), &Cache, fv.nextNode, fv.nextContent, &contact, &socialGraph);
		}
		

	}

	for (double stime = 0.0; stime < simTime; stime += 2.0) {	
		Simulator::Schedule(Seconds(stime), &updateStatics);
	}
	
	if(basicConfigInstance->is_completeScenario()){
		std::vector<double> clusterTimeSchedule;
		for (double csTime = 4.0; csTime < simTime; csTime += 200.0) {	

			//Simulator::Schedule(Seconds(csTime), &testeClusterElection, &contact, &tcontact);
			Simulator::Schedule(Seconds(csTime), &clusterHeadElection, ueNodes, &contact/*, &cluster*/, ueWiFiIface, &tcontact);
			clusterTimeSchedule.push_back(csTime);	
		}
		
		std::vector<double> updateTimeSchedule;
		
		
		ApplicationContainer ConsensusContainer;
		
		for (int i = 0; i < numberOfNodes; ++i)
		{
			Ptr<ConsensusApp> consensus = CreateObject<ConsensusApp> ();
			consensus->Setup(10, vStime, clusterTimeSchedule, updateTimeSchedule, ueWiFiIface, ueNodes, simTime);
	        ueNodes.Get(i)->AddApplication(consensus);
	        ConsensusContainer.Add (consensus);
	        
	    }
	    ConsensusContainer.Start(Seconds (0));
	  	ConsensusContainer.Stop(Seconds (simTime));
	}
	// csTime = 4.0;
	// while (csTime < simTime - 2.0) {
		
	// 	double f = (double)rand() / RAND_MAX;
 	//double t = csTime + f * (2) - 1;
    	
	// 	std::cout <<"Teste: " << t << std::endl;
	// 	//gp_tmp = gp;
	// 	Simulator::Schedule(Seconds(t), &clusterSetUp, ueNodes, &cluster, ueWiFiIface, &contact);
		
	// 	csTime++;
	// 	//random_number =  .0001*(rand() % 10000);
	// }
	// lteHelper->EnablePhyTraces ();
 //   	lteHelper->EnableMacTraces ();
   	
	// lteHelper->EnableTraces ();

	Simulator::Schedule(Seconds(simTime), &showStatistic, cont);

	// Config::ConnectWithoutContext("/NodeList/*/DeviceList/*/$ns3::WifiNetDevice/Phy/UpdateLinkQuality", MakeCallback(&UpdateQualityLink));


	// Create route file
	Ptr<OutputStreamWrapper> routingStream = Create<OutputStreamWrapper>("wifi-simple-adhoc-grid.routes", std::ios::out);
	olsr.PrintRoutingTableAllEvery(Seconds(2), routingStream);

	FlowMonitorHelper flowmon;
	Ptr<FlowMonitor> monitor = flowmon.Install(ueNodes);
	flowmon.Install(remoteHostContainer);

	for (double stime = 0.0; stime < simTime; stime += 2.0) {	
		Simulator::Schedule(Seconds(stime), &updateNetStats, &flowmon, monitor, ueLteIface);
	}

	Simulator::Stop(Seconds(simTime));
	Simulator::Run();

	/*GtkConfigStore config;
	config.ConfigureAttributes();*/
	std::string fileNameWithNoExtension = "Flowmon";
	std::string graphicsFileName        = fileNameWithNoExtension + ".png";
	std::string plotFileName            = fileNameWithNoExtension + ".plt";
	std::string plotTitle               = "Flow vs Throughput";
	std::string dataTitle               = "Throughput";

	Gnuplot gnuplot (graphicsFileName);
	gnuplot.SetTitle (plotTitle);
	gnuplot.SetTerminal ("png");
	gnuplot.SetLegend ("Flow", "Throughput");
	Gnuplot2dDataset dataset;
	dataset.SetTitle (dataTitle);
	dataset.SetStyle (Gnuplot2dDataset::LINES_POINTS);
	
	// Print per flow statistics
	Ptr<Node> nodeTemp = remoteHostContainer.Get(0);			  // Get pointer to ith node in container
	Ptr<Ipv4> ipv4Temp = nodeTemp->GetObject<Ipv4>();			  // Get Ipv4 instance of the node
	Ipv4Address remoteIP = ipv4Temp->GetAddress(1, 0).GetLocal(); // Get Ipv4InterfaceAddress of xth interface

	monitor->CheckForLostPackets();
	Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier>(flowmon.GetClassifier());
	FlowMonitor::FlowStatsContainer stats = monitor->GetFlowStats();
	double flowTime;
	double Throughput;
	double lostPackets = 0.0;
	double receivedPackets = 0.0;
	double txPackets = 0.0;
	int count = 0;
	for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin(); i != stats.end(); ++i)
	{
		Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow(i->first);
		if (t.protocol != 17){
			std::cout << "\nFlow " << i->first << " (" << t.sourceAddress << " -> " << t.destinationAddress << ")\n";
			std::cout << "  Tx Packets: " << i->second.txPackets << "\n";
			std::cout << "  Tx Bytes:   " << i->second.txBytes << "\n";
			std::cout << "  TxOffered:  " << i->second.txBytes * 8.0 / (i->second.timeLastRxPacket.GetSeconds() - i->second.timeFirstTxPacket.GetSeconds()) / 1000 / 1000 << " Mbps\n";
			std::cout << "  Rx Packets: " << i->second.rxPackets << "\n";
			std::cout << "  Rx Bytes:   " << i->second.rxBytes << "\n";
			double lostPckt = i->second.txPackets - i->second.rxPackets;
			txPackets += i->second.txPackets;
			lostPackets += lostPckt;
			receivedPackets += i->second.rxPackets;
			if (lostPckt > 0.0)
			{
				std::cout << "  Lost Packets:   " << lostPckt << "\n";
			}

			flowTime = i->second.timeLastRxPacket.GetSeconds() - i->second.timeFirstTxPacket.GetSeconds();

			Throughput = i->second.rxBytes * 8.0 / flowTime / 1024 / 1024;
			std::cout << "  Throughput: " << Throughput << " Mbps\n";
			
			if(i->second.rxPackets > 0) {
				resultsInstance.m_delayAvg += i->second.delaySum.GetMilliSeconds() / i->second.rxPackets;
				if (i->second.rxPackets > 1) {
					resultsInstance.m_jitterAvg += i->second.jitterSum.GetMilliSeconds() / (i->second.rxPackets -1);
				}
				count++;
			}
			resultsInstance.m_throughput += (i->second.rxBytes * 8.0 / (i->second.timeLastRxPacket.GetSeconds() - i->second.timeFirstTxPacket.GetSeconds()) / 1024 / 1024);
			resultsInstance.m_rxPackets +=  i->second.rxPackets;
			resultsInstance.m_lostPackets +=  i->second.lostPackets;

			dataset.Add(i->first, Throughput);
			

			if (remoteIP == t.sourceAddress)
			{
				double Pd = Throughput * 51.97 * flowTime;
				std::cout << " Promotion: "
						  << "1210.7 mW\n";
				std::cout << " Download power: " << Pd << " mW\n";
				std::cout << " Flow duration: " << flowTime << " s\n";
			}
			if (remoteIP == t.destinationAddress)
			{
				double Pu = Throughput * 438.39 * flowTime;
				std::cout << " Upload power: " << Pu << " mW\n";
			}
		}
	}
	double packetLossRate = lostPackets/(receivedPackets + lostPackets);
	auto consensusIsntance = ConsensusTimeStatics::getInstance();
	std::cout << "Packets overheaded: " << consensusIsntance->getControlPackets() << "\n";
	std::cout << "Communication Overhead: " << consensusIsntance->getControlPackets()/txPackets << "\n";
	std::cout << "Total Packets: " << txPackets << "\n";
	std::cout << "Packet Loss Rate: " << packetLossRate << "\n";
	std::cout << "Lost Packets: " << resultsInstance.m_lostPackets << "\n";
	std::cout << "Throughput: " << resultsInstance.m_throughput/(numberOfNodes * 75) / 1024 << "Mbps" << "\n";
	//std::cout << "CONSENSE THROUGHPUT: " << resultsInstance.m_throughput_consense << "Mbps" << "\n";
	std::cout << "Average Jitter: " << resultsInstance.m_jitterAvg / count << "\n";
	std::cout << "Average Delay: " << resultsInstance.m_delayAvg / count << "\n";
	std::cout << "Average Consensus Finish Time: " << consensusIsntance->getAverageTime() << "\n";
	consensusIsntance->getTimePerClusterSize();

	gnuplot.AddDataset (dataset);
	std::ofstream plotFile (plotFileName.c_str());
	gnuplot.GenerateOutput (plotFile);

	plotFile.close ();

	flowmon.SerializeToXmlFile("FlowmonResults.xml", true, true);

	Simulator::Destroy();
	//return 0;
}


int main(int argc, char *argv[])
{
	BasicConfig* basicConfigInstance = BasicConfig::getInstance();
	char opt_scenario;
	std::cout << "Select the test scenario\n";
	std::cout << "(1) Field\n";
	std::cout << "(2) Nodes\n";
	std::cout << "(3) Default\n";
	std::cout << "(4) Only Direct Trust\n";
	std::cout << "(5) No Sec\n";
	std::cout << "(6) (dTrust) On/Off behavior\n";
	std::cout << "(7) (dTrust_iTrust) On/Off behavior\n";
	std::cout << "Your option: ";
	std::cin >> opt_scenario;
	switch(opt_scenario){
		case '1':{
			basicConfigInstance->setCompleteScenario();
			basicConfigInstance->setNumberOfNodes(15);
		    basicConfigInstance->setTrustTreshold(0.5);
		    basicConfigInstance->setMaliciousNodesProportion(proportion);
		    numberOfNodes = basicConfigInstance->getNumberOfNodes();
			char op;
			std::cout << "(a) 100x100\n";
			std::cout << "(b) 200x200\n";
			std::cout << "(c) 300x300\n";
			std::cout << "(d) 400x400\n";
			std::cout << "(e) 500x500\n";
			std::cout << "(f) All\n";
			std::cout << "Your option: ";
			std::cin >> op;
			switch(op){
				case 'a':{
					cascadeBlockSec(argc, argv, "scenario/field/100_100_field.ns_movements");
					break;
				}
				case 'b':{
					cascadeBlockSec(argc, argv, "scenario/field/200_200_field.ns_movements");
					break;
				}
				case 'c':{
					cascadeBlockSec(argc, argv, "scenario/field/300_300_field.ns_movements");
					break;
				}
				case 'd':{
					cascadeBlockSec(argc, argv, "scenario/field/400_400_field.ns_movements");
					break;
				}
				case 'e':{
					cascadeBlockSec(argc, argv, "scenario/field/500_500_field.ns_movements");
					break;
				}
				case 'f':{
					cascadeBlockSec(argc, argv, "scenario/field/100_100_field.ns_movements");
					cascadeBlockSec(argc, argv, "scenario/field/200_200_field.ns_movements");
					cascadeBlockSec(argc, argv, "scenario/field/300_300_field.ns_movements");
					cascadeBlockSec(argc, argv, "scenario/field/400_400_field.ns_movements");
					cascadeBlockSec(argc, argv, "scenario/field/500_500_field.ns_movements");
					break;
				}
				default:{
					cascadeBlockSec(argc, argv, "scenario/field/100_100_field.ns_movements");
					break;
				}
			}
			break;
		}
		case '2':{
			basicConfigInstance->setCompleteScenario();
			char op;
			std::cout << "(a) 15 nodes\n";
			std::cout << "(b) 25 nodes\n";
			std::cout << "(c) 35 nodes\n";
			std::cout << "(d) 45 nodes\n";
			std::cout << "(e) 55 nodes\n";
			std::cout << "(f) All\n";
			std::cout << "Your option: ";
			std::cin >> op;
			switch(op){
				case 'a':{
				    basicConfigInstance->setNumberOfNodes(15);
				    basicConfigInstance->setTrustTreshold(0.5);
				    basicConfigInstance->setMaliciousNodesProportion(proportion);
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/15nodes_nNodes.ns_movements");
					break;
				}
				case 'b':{
					basicConfigInstance->setNumberOfNodes(25);
				    basicConfigInstance->setTrustTreshold(0.5);
				    basicConfigInstance->setMaliciousNodesProportion(proportion);
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'c':{
					basicConfigInstance->setNumberOfNodes(35);
				    basicConfigInstance->setTrustTreshold(0.5);
				    basicConfigInstance->setMaliciousNodesProportion(proportion);
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/35nodes_nNodes.ns_movements");
					break;
				}
				case 'd':{
					basicConfigInstance->setNumberOfNodes(45);
				    basicConfigInstance->setTrustTreshold(0.5);
				    basicConfigInstance->setMaliciousNodesProportion(proportion);
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/45nodes_nNodes.ns_movements");
					break;
				}
				case 'e':{
					basicConfigInstance->setNumberOfNodes(55);
				    basicConfigInstance->setTrustTreshold(0.5);
				    basicConfigInstance->setMaliciousNodesProportion(proportion);
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/55nodes_nNodes.ns_movements");
					break;
				}
				case 'f':{
					basicConfigInstance->setNumberOfNodes(15);
				    basicConfigInstance->setTrustTreshold(0.5);
				    basicConfigInstance->setMaliciousNodesProportion(proportion);
				    numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/15nodes_nNodes.ns_movements");
					basicConfigInstance->setNumberOfNodes(25);
				    basicConfigInstance->setTrustTreshold(0.5);
				    basicConfigInstance->setMaliciousNodesProportion(proportion);
				    numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					basicConfigInstance->setNumberOfNodes(35);
				    basicConfigInstance->setTrustTreshold(0.5);
				    basicConfigInstance->setMaliciousNodesProportion(proportion);
				    numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/35nodes_nNodes.ns_movements");
					basicConfigInstance->setNumberOfNodes(45);
				    basicConfigInstance->setTrustTreshold(0.5);
				    basicConfigInstance->setMaliciousNodesProportion(proportion);
				    numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/45nodes_nNodes.ns_movements");
					basicConfigInstance->setNumberOfNodes(55);
				    basicConfigInstance->setTrustTreshold(0.5);
				    basicConfigInstance->setMaliciousNodesProportion(proportion);
				    numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/55nodes_nNodes.ns_movements");
					break;
				}
				default:{
					std::cout << "Wrong option!\nExiting...\n";		
					break;
				}
			}
			break;
		}
		case '3':{
			basicConfigInstance->setNumberOfNodes(25);
		    basicConfigInstance->setTrustTreshold(0.5);
			basicConfigInstance->setCompleteScenario();
			char opB;
			std::cout << "Variate Behavior: " << std::endl;
			std::cout << "(y/Y) Yes\n";
			std::cout << "(n/N) No\n";
			std::cout << "Your option: ";
			std::cin >> opB;
			switch(opB){
				case 'y':{
					basicConfigInstance->variateMaliciousNodes();
					break;
				}
				case 'Y':{
					basicConfigInstance->variateMaliciousNodes();
					break;
				}
				case 'n':{
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				case 'N':{
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				default:{
					std::cout << "Wrong option!\nExiting...\n";
					exit(0);
				}
				break;
			}
			char op;
			std::cout << "Malicious nodes proportion: " << std::endl;
			std::cout << "(a) 0.1 nodes\n";
			std::cout << "(b) 0.2 nodes\n";
			std::cout << "(c) 0.3 nodes\n";
			std::cout << "(d) 0.4 nodes\n";
			std::cout << "(e) 0.5 nodes\n";
			std::cout << "(f) 0.6 nodes\n";
			std::cout << "(g) 0.7 nodes\n";
			std::cout << "Your option: ";
			std::cin >> op;
			switch(op){
				case 'a':{
					basicConfigInstance->setMaliciousNodesProportion(0.1);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'b':{
					basicConfigInstance->setMaliciousNodesProportion(0.2);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'c':{
					basicConfigInstance->setMaliciousNodesProportion(0.3);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'd':{
					basicConfigInstance->setMaliciousNodesProportion(0.4);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'e':{
					basicConfigInstance->setMaliciousNodesProportion(0.5);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'f':{
					basicConfigInstance->setMaliciousNodesProportion(0.6);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'g':{
					basicConfigInstance->setMaliciousNodesProportion(0.7);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				
				default:{
					std::cout << "Wrong option!\nExiting...\n";		
					break;
				}
			}
			break;
		}
		case '4':{
			basicConfigInstance->setNumberOfNodes(25);
		    basicConfigInstance->setTrustTreshold(0.5);
			basicConfigInstance->setDTrustScenario();
			char opB;
			std::cout << "Variate Behavior: " << std::endl;
			std::cout << "(y/Y) Yes\n";
			std::cout << "(n/N) No\n";
			std::cout << "Your option: ";
			std::cin >> opB;
			switch(opB){
				case 'y':{
					basicConfigInstance->variateMaliciousNodes();
					break;
				}
				case 'Y':{
					basicConfigInstance->variateMaliciousNodes();
					break;
				}
				case 'n':{
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				case 'N':{
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				default:{
					std::cout << "Wrong option!\nExiting...\n";
					exit(0);
				}
				break;
			}
			char op;
			std::cout << "Malicious nodes proportion: " << std::endl;
			std::cout << "(a) 0.1 nodes\n";
			std::cout << "(b) 0.2 nodes\n";
			std::cout << "(c) 0.3 nodes\n";
			std::cout << "(d) 0.4 nodes\n";
			std::cout << "(e) 0.5 nodes\n";
			std::cout << "(f) 0.6 nodes\n";
			std::cout << "(g) 0.7 nodes\n";
			std::cout << "Your option: ";
			std::cin >> op;
			switch(op){
				case 'a':{
					basicConfigInstance->setMaliciousNodesProportion(0.1);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'b':{
					basicConfigInstance->setMaliciousNodesProportion(0.2);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'c':{
					basicConfigInstance->setMaliciousNodesProportion(0.3);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'd':{
					basicConfigInstance->setMaliciousNodesProportion(0.4);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'e':{
					basicConfigInstance->setMaliciousNodesProportion(0.5);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'f':{
					basicConfigInstance->setMaliciousNodesProportion(0.6);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'g':{
					basicConfigInstance->setMaliciousNodesProportion(0.7);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				
				default:{
					std::cout << "Wrong option!\nExiting...\n";		
					break;
				}
			}
			break;
		}
		case '5':{
			basicConfigInstance->setNumberOfNodes(15);
		    basicConfigInstance->setTrustTreshold(0.5);
			basicConfigInstance->setNoTrustScenario();
			char opB;
			std::cout << "Variate Behavior: " << std::endl;
			std::cout << "(y/Y) Yes\n";
			std::cout << "(n/N) No\n";
			std::cout << "Your option: ";
			std::cin >> opB;
			switch(opB){
				case 'y':{
					basicConfigInstance->variateMaliciousNodes();
					break;
				}
				case 'Y':{
					basicConfigInstance->variateMaliciousNodes();
					break;
				}
				case 'n':{
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				case 'N':{
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				default:{
					std::cout << "Wrong option!\nExiting...\n";
					exit(0);
				}
				break;
			}
			char op;
			std::cout << "Malicious nodes proportion: " << std::endl;
			std::cout << "(a) 0.1 nodes\n";
			std::cout << "(b) 0.2 nodes\n";
			std::cout << "(c) 0.3 nodes\n";
			std::cout << "(d) 0.4 nodes\n";
			std::cout << "(e) 0.5 nodes\n";
			std::cout << "(f) 0.6 nodes\n";
			std::cout << "(g) 0.7 nodes\n";
			std::cout << "Your option: ";
			std::cin >> op;
			switch(op){
				case 'a':{
					basicConfigInstance->setMaliciousNodesProportion(0.1);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/15nodes_nNodes.ns_movements");
					break;
				}
				case 'b':{
					basicConfigInstance->setMaliciousNodesProportion(0.2);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/15nodes_nNodes.ns_movements");
					break;
				}
				case 'c':{
					basicConfigInstance->setMaliciousNodesProportion(0.3);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/15nodes_nNodes.ns_movements");
					break;
				}
				case 'd':{
					basicConfigInstance->setMaliciousNodesProportion(0.4);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/15nodes_nNodes.ns_movements");
					break;
				}
				case 'e':{
					basicConfigInstance->setMaliciousNodesProportion(0.5);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/15nodes_nNodes.ns_movements");
					break;
				}
				case 'f':{
					basicConfigInstance->setMaliciousNodesProportion(0.6);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/15nodes_nNodes.ns_movements");
					break;
				}
				case 'g':{
					basicConfigInstance->setMaliciousNodesProportion(0.7);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/15nodes_nNodes.ns_movements");
					break;
				}
				
				default:{
					std::cout << "Wrong option!\nExiting...\n";		
					break;
				}
			}
			break;
		}
		case '6':{
			basicConfigInstance->setNumberOfNodes(25);
		    basicConfigInstance->setTrustTreshold(0.5);
			basicConfigInstance->setDTrustScenario();
			char opB;
			std::cout << "Malicious probability: " << std::endl;
			std::cout << "(1) 90\n";
			std::cout << "(2) 80\n";
			std::cout << "(3) 70\n";
			std::cout << "(4) 60\n";
			std::cout << "(5) 50\n";
			std::cout << "Your option: ";
			std::cin >> opB;
			switch(opB){
				case '1':{
					basicConfigInstance->setOnOffBehavior(0.9);
					break;
				}
				case '2':{
					basicConfigInstance->setOnOffBehavior(0.8);
					break;
				}
				case '3':{
					basicConfigInstance->setOnOffBehavior(0.7);
					break;
				}
				case '4':{
					basicConfigInstance->setOnOffBehavior(0.6);
					break;
				}
				case '5':{
					basicConfigInstance->setOnOffBehavior(0.5);
					break;
				}
				default:{
					std::cout << "Wrong option!\nExiting...\n";
					exit(0);
				}
				break;
			}
			char op;
			std::cout << "Malicious nodes proportion: " << std::endl;
			std::cout << "(a) 0.1 nodes\n";
			std::cout << "(b) 0.2 nodes\n";
			std::cout << "(c) 0.3 nodes\n";
			std::cout << "(d) 0.4 nodes\n";
			std::cout << "(e) 0.5 nodes\n";
			std::cout << "(f) 0.6 nodes\n";
			std::cout << "(g) 0.7 nodes\n";
			std::cout << "Your option: ";
			std::cin >> op;
			switch(op){
				case 'a':{
					basicConfigInstance->setMaliciousNodesProportion(0.1);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'b':{
					basicConfigInstance->setMaliciousNodesProportion(0.2);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'c':{
					basicConfigInstance->setMaliciousNodesProportion(0.3);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'd':{
					basicConfigInstance->setMaliciousNodesProportion(0.4);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'e':{
					basicConfigInstance->setMaliciousNodesProportion(0.5);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'f':{
					basicConfigInstance->setMaliciousNodesProportion(0.6);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'g':{
					basicConfigInstance->setMaliciousNodesProportion(0.7);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				
				default:{
					std::cout << "Wrong option!\nExiting...\n";		
					break;
				}
			}
			break;
		}
		case '7':{
			basicConfigInstance->setNumberOfNodes(25);
		    basicConfigInstance->setTrustTreshold(0.5);
			basicConfigInstance->setCompleteScenario();
			char opB;
			std::cout << "Malicious probability: " << std::endl;
			std::cout << "(1) 90\n";
			std::cout << "(2) 80\n";
			std::cout << "(3) 70\n";
			std::cout << "(4) 60\n";
			std::cout << "(5) 50\n";
			std::cout << "Your option: ";
			std::cin >> opB;
			switch(opB){
				case '1':{
					basicConfigInstance->setOnOffBehavior(0.9);
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				case '2':{
					basicConfigInstance->setOnOffBehavior(0.8);
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				case '3':{
					basicConfigInstance->setOnOffBehavior(0.7);
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				case '4':{
					basicConfigInstance->setOnOffBehavior(0.6);
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				case '5':{
					basicConfigInstance->setOnOffBehavior(0.5);
					basicConfigInstance->notVariateMaliciousNodes();
					break;
				}
				default:{
					std::cout << "Wrong option!\nExiting...\n";
					exit(0);
				}
				break;
			}
			char op;
			std::cout << "Malicious nodes proportion: " << std::endl;
			std::cout << "(a) 0.1 nodes\n";
			std::cout << "(b) 0.2 nodes\n";
			std::cout << "(c) 0.3 nodes\n";
			std::cout << "(d) 0.4 nodes\n";
			std::cout << "(e) 0.5 nodes\n";
			std::cout << "(f) 0.6 nodes\n";
			std::cout << "(g) 0.7 nodes\n";
			std::cout << "Your option: ";
			std::cin >> op;
			switch(op){
				case 'a':{
					basicConfigInstance->setMaliciousNodesProportion(0.1);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'b':{
					basicConfigInstance->setMaliciousNodesProportion(0.2);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'c':{
					basicConfigInstance->setMaliciousNodesProportion(0.3);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'd':{
					basicConfigInstance->setMaliciousNodesProportion(0.4);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'e':{
					basicConfigInstance->setMaliciousNodesProportion(0.5);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'f':{
					basicConfigInstance->setMaliciousNodesProportion(0.6);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				case 'g':{
					basicConfigInstance->setMaliciousNodesProportion(0.7);
					proportion = basicConfigInstance->getMaliciousNodesProportion();
					numberOfNodes = basicConfigInstance->getNumberOfNodes();
					cascadeBlockSec(argc, argv, "scenario/nodes/25nodes_nNodes.ns_movements");
					break;
				}
				
				default:{
					std::cout << "Wrong option!\nExiting...\n";		
					break;
				}
			}
			break;
		}
		default:{
			std::cout << "Wrong option!\nExiting...\n";
			break;
		}
	}
	
	return 0;
}

void updateStatics(){
	double sTime = Simulator::Now().GetSeconds();
	
	std::cout << "Update statics " << sTime << "s" << std::endl;

	// Trust* trustInstance = Trust::getInstance(numberOfNodes);

	TrustStatistics* trustStaticsInstance = TrustStatistics::getInstance(numberOfNodes);

	MaliciousNodes* maliciousInstance = MaliciousNodes::getInstance(proportion, numberOfNodes);

	for (int i = 0; i < numberOfNodes; i++){
		if(maliciousInstance->isMalicious(i)){
			trustStaticsInstance->setAvgTrustOfMalNode(i, (int) sTime);
		} else {
			trustStaticsInstance->setAvgTrustOfNormalNode(i, (int) sTime);
		}
		trustStaticsInstance->setGeneralAvgTrustOfNodes(i, (int) sTime);
	}
	trustStaticsInstance->setNodeToAnalyzeTrust((int) sTime);
	trustStaticsInstance->setBadContentDeliveryRate((int) sTime);
	trustStaticsInstance->print((int) sTime);
}

bool isUEnode(Ipv4Address node, Ipv4InterfaceContainer ueLteIface){
	for (int i = 0; i < 15; i++){
		if(ueLteIface.GetAddress(i) == node) {
			return true;
		}
	}
	return false;
}

void updateNetStats(FlowMonitorHelper* flowmon, Ptr<FlowMonitor> monitor, Ipv4InterfaceContainer ueLteIface) {
	NetStats* netStatsInstance = NetStats::getInstance(numberOfNodes);
	double sTime = Simulator::Now().GetSeconds();

	monitor->CheckForLostPackets();
	Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier>(flowmon->GetClassifier());
	FlowMonitor::FlowStatsContainer stats = monitor->GetFlowStats();
	
	double throughput_flows = 0.0;
	double lostPackets = 0.0;
	double receivedPackets = 0.0;
	double delay = 0.0;
	double jitter = 0.0;
	int flow_cont = 0;
	for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin(); i != stats.end(); ++i)
	{
		Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow(i->first);

        if (t.protocol != 17){
        	lostPackets += i->second.txPackets - i->second.rxPackets;
			receivedPackets += i->second.rxPackets;
			double ftime = i->second.timeLastRxPacket.GetSeconds() - i->second.timeFirstTxPacket.GetSeconds();
			throughput_flows += (i->second.rxBytes * 8.0 / ftime / 1024);
			if(i->second.rxPackets > 0) {
				delay += i->second.delaySum.GetMilliSeconds() / i->second.rxPackets;
				if (i->second.rxPackets > 1) {
					jitter += i->second.jitterSum.GetMilliSeconds() / (i->second.rxPackets -1);
				}
			}
        }
        flow_cont++;
	}

	netStatsInstance->updateVazao((int)sTime, throughput_flows/flow_cont);
	netStatsInstance->updatePlr((int)sTime, lostPackets/(lostPackets+receivedPackets));
	netStatsInstance->updateDelay((int)sTime, delay/flow_cont);
	netStatsInstance->updateJitter((int)sTime, jitter/flow_cont);
	netStatsInstance->print((int) sTime );
}

void 
UpdateQualityLink(int sourceVertice, int destinationVertice, double snr, double sinr, double per, bool received) {
	std::cout << "SNR: " << snr << "\nSINR: " << sinr << "\n";
}


double
lowerIncompleteGamma(double x, void *params)
{
	// (t^(k-1))*(exp(-t))

	// The next line recovers alpha from the passed params pointer
	double alpha = *(double *)params;

	return ((pow(x, (alpha - 1))) * (exp(-x)));
}

double
generalisedPareto(double mu, double sigma, double E, double simTime)
{
	//double gp = 0.0;
	//gsl_rng * r1 = gsl_rng_alloc (gsl_rng_rand);
	//double U = gsl_rng_uniform_pos (r1);

	Ptr<UniformRandomVariable> u = CreateObject<UniformRandomVariable>();
	u->SetAttribute("Min", DoubleValue(0.0));
	u->SetAttribute("Max", DoubleValue(1.0));

	double U = 0.0;
	//  double gp = 0.0;
	double gp = std::numeric_limits<int>::max();

	while (gp > simTime)
	{
		U = u->GetValue();
		gp = fabs(mu + sigma * (pow(U, -E) - 1) / E);
	}
	return gp;
}

double
Weibull(double lambda, double k)
{
	//gsl_rng * r2 = gsl_rng_alloc (gsl_rng_rand);
	//double wb = gsl_ran_weibull(x, k, lambda);

	Ptr<WeibullRandomVariable> w = CreateObject<WeibullRandomVariable>();
	w->SetAttribute("Scale", DoubleValue(lambda));
	w->SetAttribute("Shape", DoubleValue(k));

	double wb = w->GetValue();

	return wb;
}

double
closeness(double k, double teta)
{
	gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(1000);

	double lower_limit = 0;			  /* lower limit a */
	double upper_limit = 50.0 / teta; /* upper limit b */
	double abs_error = 1.0e-8;		  /* to avoid round-off problems */
	double rel_error = 1.0e-8;		  /* the result will usually be much better */
	double result;					  /* the result from the integration */
	double error;					  /* the estimated error from the integration */

	double alpha = k; // parameter in integrand

	gsl_function My_function;
	void *params_ptr = &alpha;

	My_function.function = &lowerIncompleteGamma;
	My_function.params = params_ptr;

	gsl_integration_qags(&My_function, lower_limit, upper_limit,
						 abs_error, rel_error, 1000, work_ptr, &result, &error);

	// cout  std::cout.precision (18);
	// cout  std::cout << "result          = " << result << std::endl;
	// cout  std::cout << "estimated error = " << error << std::endl;

	double gamma = gsl_sf_gamma(alpha);

	//  std::cout << "gamma          = " << gamma << std::endl;

	double wij = 1 - (result / gamma);

	//  std::cout << "Closeness (wij) = " << wij << std::endl;

	gsl_integration_workspace_free(work_ptr);

	return (wij);
}

void PrintRoutingTable(Ptr<Node> node)
{

	Ptr<WifiNetDevice> dev = node->GetObject<WifiNetDevice>();
	Ptr<ns3::Ipv4RoutingProtocol> rp(node->GetObject<Ipv4>()->GetRoutingProtocol());
	Ptr<ns3::olsr::RoutingProtocol> proto = Ipv4RoutingHelper::GetRouting<ns3::olsr::RoutingProtocol>(rp);
	std::vector<ns3::olsr::RoutingTableEntry> entry = proto->GetRoutingTableEntries();

	//Ptr<ns3::olsr::RoutingProtocol> routing = node->GetObject<ns3::olsr::RoutingProtocol>();
	// Print  routing table entries for OLSR routing
	//std::vector<ns3::olsr::RoutingTableEntry> entry = routing->GetRoutingTableEntries();
	//std::cout << "Routing table for device: " << Names::FindName(node) << std::endl;
	std::cout << "Routing table for device: " << node->GetId() << std::endl;
	std::cout << "DestinyAddress\t\tNextAddress\t\tInterface\t\tDistance\n";
	for (std::vector<ns3::olsr::RoutingTableEntry>::const_iterator i = entry.begin(); i != entry.end(); i++)
	{
		std::cout << i->destAddr << "\t\t"
				  << i->nextAddr << "\t\t"
				  << i->interface << "\t\t"
				  << i->distance << std::endl;
	}
}

void GenerateECDSAKeys(int idx, Key *keyUtil)
{
	DEV_ECDSA ecdsa;
	std::string publicKey;
	std::string privateKey;
	ecdsa.GenerateKeys(privateKey, publicKey);

	keyUtil->putECDSA_publicKeyInMap(idx, publicKey);
	keyUtil->putECDSA_privateKeyInMap(idx, privateKey);
}

void GenerateECDHKeys(int idx, Key *keyUtil)
{
	DEV_ECDH ecdh;
	std::string publicKey;
	std::string privateKey;
	ecdh.GenerateKeys(1024,privateKey, publicKey);

	keyUtil->putECC_publicKeyInMap(idx, publicKey);
	keyUtil->putECC_privateKeyInMap(idx, privateKey);
}

void GenerateKeys(int numberOfNodes)
{

	Key *keyUtil = Key::getInstance();

	std::cout << "keys - Pares de chaves ECDH e ECDSA";
	for (int node = 0; node < numberOfNodes; node++)
	{
		//std::string privKey, pubKey;
		GenerateECDSAKeys(node, keyUtil);

		GenerateECDHKeys(node, keyUtil);
		//std::cout << "Node (" << node << "): " << " Gerou par de chaves ECDSA e ECDH " << std::endl;
		//std::cout << "\tTempo: " << timestamp << std::endl;
	}
	std::cout << " gerados" << std::endl;
	/*for (int node = 0; node < numberOfNodes; node++) {
		std::std::vector<int> v = keyUtil->getPubKeyNodes();
	}*/
}
//RoutingTableToMatrix (NodeContainer ueNodes, int numberOfNodes, Graph *g)
void RoutingTableToMatrix(NodeContainer ueNodes, int numberOfNodes, Graph *g, Contact *contact, ContactDistribution *contDist)
{
	//std::cout << "Fazendo copia: ";
	Graph *gCopy = new Graph(*g);
	//TwoHopGraph *tgCopy = new TwoHopGraph(*tg);
	//std::cout << "OK " << std::endl;
	//gCopy.printEdge();
	//
	g->zeros();
	//tg->zeros();

	for (int nodeN = 0; nodeN < numberOfNodes; nodeN++)
	{
		Ptr<Node> node = ueNodes.Get(nodeN);

		for (int nodeC = 0; nodeC < numberOfNodes; nodeC++)
		{
			Ptr<Node> nC = ueNodes.Get(nodeC);

			Ptr<WifiNetDevice> dev = node->GetObject<WifiNetDevice>();
			Ptr<ns3::Ipv4RoutingProtocol> rp(node->GetObject<Ipv4>()->GetRoutingProtocol());
			Ptr<ns3::olsr::RoutingProtocol> proto = Ipv4RoutingHelper::GetRouting<ns3::olsr::RoutingProtocol>(rp);
			std::vector<ns3::olsr::RoutingTableEntry> entry = proto->GetRoutingTableEntries();

			for (std::vector<ns3::olsr::RoutingTableEntry>::const_iterator i = entry.begin(); i != entry.end(); i++)
			{
				int nodeInterface = nC->GetObject<Ipv4L3Protocol>()->GetInterfaceForAddress(i->destAddr);
				//std::cout << nodeInterface << std::endl;
				if (nodeInterface != -1)
				{
					g->addEdge(node->GetId(), nC->GetId());

					//std::cout << i->destAddr << std::endl;
					//g->printEdge();
				}
				/*nodeInterface = nC->GetObject<Ipv4L3Protocol>()->GetInterfaceForAddress(i->nextAddr);
				if (nodeInterface != -1)
				{
					tg->addEdge(node->GetId(), nC->GetId());

					//std::cout << i->destAddr << std::endl;
					//g->printEdge();
				}*/
			}
		}
		
	}

	for (int nodeN = 0; nodeN < numberOfNodes; nodeN++)
	{
		for (int nodeC = 0; nodeC < numberOfNodes; nodeC++)
		{
			if ((g->isEdge(nodeN, nodeC)) && !(gCopy->isEdge(nodeN, nodeC)))
			{
				contact->contactBegin(nodeN, nodeC, Simulator::Now().GetSeconds());
				contDist->insertCont(nodeC);
			}
			else if (!(g->isEdge(nodeN, nodeC)) && (gCopy->isEdge(nodeN, nodeC)))
			{
				contDist->insertDur(Simulator::Now().GetSeconds() - contact->getlastContact(nodeN, nodeC));
				contact->contactEnd(nodeN, nodeC, Simulator::Now().GetSeconds());
			}
			else if ((g->isEdge(nodeN, nodeC)) && (gCopy->isEdge(nodeN, nodeC)))
			{
				contact->contactUpdate(nodeN, nodeC, Simulator::Now().GetSeconds());
			}
		}
	}
	//contDist->printContacts();
	/*for (int nodeN = 0; nodeN < numberOfNodes; nodeN++)
	{
		for (int nodeC = 0; nodeC < numberOfNodes; nodeC++)
		{
			if (nodeC != nodeN)
				contact->printContacts(nodeN, nodeC);
		}
	}
*/
	delete gCopy;
}
void RoutingTableToMatrixTwoHop(NodeContainer ueNodes, int numberOfNodes, TwoHopGraph *tg, TwoHopContact *tcontact, TwoHopContactDistribution *tcontDist){
	//std::cout << "Fazendo copia: ";
	TwoHopGraph *tgCopy = new TwoHopGraph(*tg);
	tg->zeros();

	for (int nodeN = 0; nodeN < numberOfNodes; nodeN++)
	{
		Ptr<Node> node = ueNodes.Get(nodeN);

		for (int nodeC = 0; nodeC < numberOfNodes; nodeC++)
		{
			Ptr<Node> nC = ueNodes.Get(nodeC);

			Ptr<WifiNetDevice> dev = node->GetObject<WifiNetDevice>();
			Ptr<ns3::Ipv4RoutingProtocol> rp(node->GetObject<Ipv4>()->GetRoutingProtocol());
			Ptr<ns3::olsr::RoutingProtocol> proto = Ipv4RoutingHelper::GetRouting<ns3::olsr::RoutingProtocol>(rp);
			std::vector<ns3::olsr::RoutingTableEntry> entry = proto->GetRoutingTableEntries();

			for (std::vector<ns3::olsr::RoutingTableEntry>::const_iterator i = entry.begin(); i != entry.end(); i++)
			{
				int nodeInterface = nC->GetObject<Ipv4L3Protocol>()->GetInterfaceForAddress(i->nextAddr);
				
				if (nodeInterface != -1)
				{
					tg->addEdge(node->GetId(), nC->GetId());

					//std::cout << i->nextAddr << std::endl;
					//g->printEdge();
				}
			}
		}
		
	}

	for (int nodeN = 0; nodeN < numberOfNodes; nodeN++)
	{
		for (int nodeC = 0; nodeC < numberOfNodes; nodeC++)
		{
			if ((tg->isEdge(nodeN, nodeC)) && !(tgCopy->isEdge(nodeN, nodeC)))
			{
				tcontact->contactBegin(nodeN, nodeC, Simulator::Now().GetSeconds());
				tcontDist->insertCont(nodeC);
			}
			else if (!(tg->isEdge(nodeN, nodeC)) && (tgCopy->isEdge(nodeN, nodeC)))
			{
				tcontDist->insertDur(Simulator::Now().GetSeconds() - tcontact->getlastContact(nodeN, nodeC));
				tcontact->contactEnd(nodeN, nodeC, Simulator::Now().GetSeconds());
			}
			else if ((tg->isEdge(nodeN, nodeC)) && (tgCopy->isEdge(nodeN, nodeC)))
			{
				tcontact->contactUpdate(nodeN, nodeC, Simulator::Now().GetSeconds());
			}
		}
	}
	//tcontDist->printContacts();
	
	/*for (int nodeN = 0; nodeN < numberOfNodes; nodeN++)
	{
		for (int nodeC = 0; nodeC < numberOfNodes; nodeC++)
		{
			if (nodeC != nodeN)
				tcontact->printContacts(nodeN, nodeC);
		}
	}*/

	delete tgCopy;
}

void addContent(Graph *sg, int node, int contentNumber, Content *cont)
{
	if (cont->nc != -1)
	{
		std::cout << "InNeighborCache(node: " << cont->nc << " )" << std::endl;
		// closeness				teste.test
		//std::cout << " closeness mean: " << contact->getMean(cont->nc,node) << std::endl;
		//std::cout << " closeness var: " << contact->getVariance(cont->nc,node) << std::endl;
		//std::cout << " closeness n: " << contact->getN(cont->nc,node) << std::endl;
		//
		//double closen = closeness(contact->getMean(cont->nc,node), contact->getVariance(cont->nc,node));
		//std::cout << " closeness: " << closen << std::endl;

		cont->addFreq(cont->nc, contentNumber);
		cont->addContent(node, contentNumber);
		cont->addCache(node, contentNumber);
		cont->cacheHit++;
		// register content obtained from friends
		if (sg->isEdge(node, cont->nc))
		{
			std::cout << "\t->FriendCache";
			cont->friendCache++;
		}
		std::cout << std::endl;
	}
}

void
//searchContent (Graph *sg, int node, int contentNumber, Content *cont, Graph *g)
searchContent(Graph *sg, int node, int contentNumber, Contact *contact)
{
	Content* cont = Content::getInstance();
	Trust* trustInstance = Trust::getInstance(numberOfNodes);
	MaliciousNodes* maliciousInstance = MaliciousNodes::getInstance();
	BlockchainConnection* blockchainConn =  BlockchainConnection::getInstance();
	TrustStatistics* trustStaticsInstance = TrustStatistics::getInstance(numberOfNodes);
	BasicConfig* basicConfigInstance = BasicConfig::getInstance();
	double sTime = Simulator::Now().GetSeconds();
	//std::cout << "Test Content Number " << contentNumber << std::endl;
	cont->nc = -1;
	//Key::getInstance()->pool_add = -1;
	//std::cout << " teste " << cont->nc << std::endl;

	std::vector<int> neighborCont;
	neighborCont.reserve(numberOfNodes);

	std::cout << std::endl
			  << "TypeContent: ";
			  
	// Caso o conteúdo seja antigo (já foi transmtido na rede D2D)
	if (cont->oldContent(contentNumber))
	{
		std::cout << "OldContent " << contentNumber << std::endl;

		// Create a vector of neighbors with the content (neighborCont)

		for (unsigned int j = 0; j < numberOfNodes; j++)
		{
			if ((contact->isEdge(node, j)) && cont->hasCache(j, contentNumber))
			{
				neighborCont.push_back(j);
			}
		}

		std::cout << "CacheStatus: ";
		//Conteúdo procurado está no cache local
		//  --> Apenas recupera o conteúdo presente
		if (cont->hasCache(node, contentNumber))
		{
			cont->lc = true;
			cont->addFreq(node, contentNumber);
			std::cout << "InCacheLocal" << std::endl;
			cont->localCache++;
		}

		//Não há vizinhos com o conteúdo em cache
		else if (neighborCont.empty())
		{
			std::cout << "NotInNeighborsCache" << std::endl;
			cont->cacheMiss++;
		}

		else if (neighborCont.size() == 1)
		{
			cont->nc = neighborCont.at(0); 
			if(basicConfigInstance->is_completeScenario()){
				if (maliciousInstance->isMalicious(cont->nc) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
					std::cout << "On/Off - malicious" << std::endl;
					trustStaticsInstance->newInvalidContent();
				}
				else if (maliciousInstance->isMalicious(cont->nc) && !basicConfigInstance->maliciousSwitched(sTime)){
					trustStaticsInstance->newInvalidContent();
				}
				if (trustInstance->is_trustable(node, cont->nc)){
					// lte
					std::cout << "InNeighborCache(node: " << cont->nc << " )" << std::endl;
					cont->addFreq(cont->nc, contentNumber);
					cont->addContent(node, contentNumber);
					// cont->addCache(node, contentNumber);
					cont->addCache(node, cont->nc, contentNumber, true);
					trustStaticsInstance->newContent();
					cont->cacheHit++;
					if (maliciousInstance->isMalicious(cont->nc) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
						std::cout << "On/Off - malicious" << std::endl;
						std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
						trustInstance->badBehavior(node, cont->nc, 1);
						trustStaticsInstance->newBadContentTransmission();
					}
					else if (maliciousInstance->isMalicious(cont->nc) && !basicConfigInstance->maliciousSwitched(sTime)){
						std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
						trustInstance->badBehavior(node, cont->nc, 1);
						trustStaticsInstance->newBadContentTransmission();
					} else {
						std::cout << "Node (" << cont->nc << ") not malicious - valid content" << std::endl;
						trustInstance->goodBehavior(node, cont->nc, 1);
					}
					std::cout << "Direct Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getDirectTrust(node, cont->nc) << std::endl;
					std::cout << "Indirect Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getIndirectTrust(node, cont->nc) << std::endl;
					std::cout << "Trust value of Node (" << node << ") on Node (" << cont->nc << ") -> " << trustInstance->getTrust(node, cont->nc) << std::endl;
					if (sg->isEdge(node, cont->nc))
					{
						std::cout << "\t->FriendCache";
						cont->friendCache++;
					}
				} else {
					std::cout << "Not adding content " << contentNumber << " on cache of node " << node << std::endl;
					std::cout << "Node (" << node << ") do not trust ("<< cont->nc << ") " << trustInstance->getTrust(node, cont->nc) << std::endl;
				}
			} else if (basicConfigInstance->is_DTrustScenario()) {
				if (maliciousInstance->isMalicious(cont->nc) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
					std::cout << "On/Off - malicious" << std::endl;
					trustStaticsInstance->newInvalidContent();
				}
				else if (maliciousInstance->isMalicious(cont->nc) && !basicConfigInstance->maliciousSwitched(sTime)){
					trustStaticsInstance->newInvalidContent();
				}
				if (trustInstance->is_trustable(node, cont->nc)){
					std::cout << "InNeighborCache(node: " << cont->nc << " )" << std::endl;
					cont->addFreq(cont->nc, contentNumber);
					cont->addContent(node, contentNumber);
					cont->addCache(node, contentNumber);
					// cont->addCache(node, cont->nc, contentNumber, true);
					trustStaticsInstance->newContent();
					cont->cacheHit++;
					if (maliciousInstance->isMalicious(cont->nc) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
						std::cout << "On/Off - malicious" << std::endl;
						std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
						trustInstance->badBehavior(node, cont->nc, 1);
						trustStaticsInstance->newBadContentTransmission();
					}
					else if (maliciousInstance->isMalicious(cont->nc) && !basicConfigInstance->maliciousSwitched(sTime)){
						std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
						trustInstance->badBehavior(node, cont->nc, 1);
						trustStaticsInstance->newBadContentTransmission();
					} else {
						std::cout << "Node (" << cont->nc << ") not malicious - valid content" << std::endl;
						trustInstance->goodBehavior(node, cont->nc, 1);
					}
					std::cout << "Direct Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getDirectTrust(node, cont->nc) << std::endl;
					// std::cout << "Indirect Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getIndirectTrust(node, cont->nc) << std::endl;
					// std::cout << "Trust value of Node (" << node << ") on Node (" << cont->nc << ") -> " << trustInstance->getTrust(node, cont->nc) << std::endl;
					if (sg->isEdge(node, cont->nc))
					{
						std::cout << "\t->FriendCache";
						cont->friendCache++;
					}
				} else {
					std::cout << "Not adding content " << contentNumber << " on cache of node " << node << std::endl;
					std::cout << "Node (" << node << ") do not trust ("<< cont->nc << ") " << trustInstance->getTrust(node, cont->nc) << std::endl;
				}
			} else if (basicConfigInstance->is_NoTrustScenario()) {
				std::cout << "InNeighborCache(node: " << cont->nc << " )" << std::endl;
				cont->addFreq(cont->nc, contentNumber);
				cont->addContent(node, contentNumber);
				cont->addCache(node, contentNumber);
				// cont->addCache(node, cont->nc, contentNumber, true);
				trustStaticsInstance->newContent();
				cont->cacheHit++;
				if (maliciousInstance->isMalicious(cont->nc) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
					std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
					trustInstance->badBehavior(node, cont->nc, 1);
					trustStaticsInstance->newBadContentTransmission();
				}
				else if (maliciousInstance->isMalicious(cont->nc) && !basicConfigInstance->maliciousSwitched(sTime)){
					std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
					//trustInstance->badBehavior(node, cont->nc, 1);
					trustStaticsInstance->newInvalidContent();
					trustStaticsInstance->newBadContentTransmission();
				} else {
					std::cout << "Node (" << cont->nc << ") not malicious - valid content" << std::endl;
					// trustInstance->goodBehavior(node, cont->nc, 1);
				}
				// std::cout << "Direct Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getDirectTrust(node, cont->nc) << std::endl;
				// std::cout << "Indirect Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getIndirectTrust(node, cont->nc) << std::endl;
				// std::cout << "Trust value of Node (" << node << ") on Node (" << cont->nc << ") -> " << trustInstance->getTrust(node, cont->nc) << std::endl;
				if (sg->isEdge(node, cont->nc))
				{
					std::cout << "\t->FriendCache";
					cont->friendCache++;
				}
			}
		}

		//Caso contrário
		//  --> Seleciona o melhor nó para distribuir o conteúdo
		else
		{
			// Select the best node to serve de content
			// random
			cont->nc = neighborCont.at(rand() % neighborCont.size()); // lte
			if (basicConfigInstance->is_completeScenario()){ // cond se verifica confiança
				if (maliciousInstance->isMalicious(cont->nc) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
					std::cout << "On/Off - malicious" << std::endl;
					trustStaticsInstance->newInvalidContent();
				}
				else if (maliciousInstance->isMalicious(cont->nc) && !basicConfigInstance->maliciousSwitched(sTime)){
					trustStaticsInstance->newInvalidContent();
				}
				if (trustInstance->is_trustable(node, cont->nc)){
					
					std::cout << "InNeighbor'sCache(randomNode: " << cont->nc << " )" << std::endl;
					// closeness				teste.test
					//std::cout << " closeness mean: " << contact->getMean(cont->nc,node) << std::endl;
					//std::cout << " closeness var: " << contact->getVariance(cont->nc,node) << std::endl;
					//std::cout << " closeness n: " << contact->getN(cont->nc,node) << std::endl;
					//
					//double closen = closeness(contact->getMean(cont->nc,node), contact->getVariance(cont->nc,node));
					//std::cout << " closeness: " << closen << std::endl;
				
					cont->addFreq(cont->nc, contentNumber);
					cont->addContent(node, contentNumber);
					cont->addCache(node, cont->nc, contentNumber, true);
					trustStaticsInstance->newContent();
					// cont->addCache(node, contentNumber);
					cont->cacheHit++;
					if (maliciousInstance->isMalicious(cont->nc) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
						std::cout << "On/Off - malicious" << std::endl;
						std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
						trustInstance->badBehavior(node, cont->nc, 1);
						trustStaticsInstance->newBadContentTransmission();
					}
					else if (maliciousInstance->isMalicious(cont->nc) && !basicConfigInstance->maliciousSwitched(sTime)){
						std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
						trustInstance->badBehavior(node, cont->nc, 1);
						trustStaticsInstance->newBadContentTransmission();
					} else {
						std::cout << "Node (" << cont->nc << ") not malicious - valid content" << std::endl;
						trustInstance->goodBehavior(node, cont->nc, 1);
					}
					std::cout << "Direct Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getDirectTrust(node, cont->nc) << std::endl;
					std::cout << "Indirect Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getIndirectTrust(node, cont->nc) << std::endl;
					std::cout << "Trust value of Node (" << node << ") on Node (" << cont->nc << ") -> " << trustInstance->getTrust(node, cont->nc) << std::endl;				// register content obtained from friends
					if (sg->isEdge(node, cont->nc))
					{
						std::cout << "\t->FriendCache";
						cont->friendCache++;
					}
				} else {
					std::cout << "Not adding content " << contentNumber << " on cache of node " << node << std::endl;
					std::cout << "Node (" << node << ") do not trust ("<< cont->nc << ") " << trustInstance->getTrust(node, cont->nc) << std::endl;
				}
			} else if (basicConfigInstance->is_DTrustScenario()){
				if (maliciousInstance->isMalicious(cont->nc) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
					std::cout << "On/Off - malicious" << std::endl;
					trustStaticsInstance->newInvalidContent();
				}
				else if (maliciousInstance->isMalicious(cont->nc) && !basicConfigInstance->maliciousSwitched(sTime)){
					trustStaticsInstance->newInvalidContent();
				}
				if (trustInstance->is_trustable(node, cont->nc)){
					// cont->nc = neighborCont.at(rand() % neighborCont.size()); // lte
					std::cout << "InNeighbor'sCache(randomNode: " << cont->nc << " )" << std::endl;
					// closeness				teste.test
					//std::cout << " closeness mean: " << contact->getMean(cont->nc,node) << std::endl;
					//std::cout << " closeness var: " << contact->getVariance(cont->nc,node) << std::endl;
					//std::cout << " closeness n: " << contact->getN(cont->nc,node) << std::endl;
					//
					//double closen = closeness(contact->getMean(cont->nc,node), contact->getVariance(cont->nc,node));
					//std::cout << " closeness: " << closen << std::endl;
				
					cont->addFreq(cont->nc, contentNumber);
					cont->addContent(node, contentNumber);
					// cont->addCache(node, cont->nc, contentNumber, true);
					trustStaticsInstance->newContent();
					cont->addCache(node, contentNumber);
					cont->cacheHit++;
					if (maliciousInstance->isMalicious(cont->nc) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
						std::cout << "On/Off - malicious" << std::endl;
						std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
						trustInstance->badBehavior(node, cont->nc, 1);
						trustStaticsInstance->newBadContentTransmission();
					}
					else if (maliciousInstance->isMalicious(cont->nc) && !basicConfigInstance->maliciousSwitched(sTime) ){
						std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
						trustInstance->badBehavior(node, cont->nc, 1);
						// trustStaticsInstance->newInvalidContent();
						trustStaticsInstance->newBadContentTransmission();
					} else {
						std::cout << "Node (" << cont->nc << ") not malicious - valid content" << std::endl;
						trustInstance->goodBehavior(node, cont->nc, 1);
					}
					std::cout << "Direct Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getDirectTrust(node, cont->nc) << std::endl;
					// std::cout << "Indirect Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getIndirectTrust(node, cont->nc) << std::endl;
					// std::cout << "Trust value of Node (" << node << ") on Node (" << cont->nc << ") -> " << trustInstance->getTrust(node, cont->nc) << std::endl;				// register content obtained from friends
					if (sg->isEdge(node, cont->nc))
					{
						std::cout << "\t->FriendCache";
						cont->friendCache++;
					}
				} else {
					std::cout << "Not adding content " << contentNumber << " on cache of node " << node << std::endl;
					std::cout << "Node (" << node << ") do not trust ("<< cont->nc << ") " << trustInstance->getTrust(node, cont->nc) << std::endl;
				}

			} else if (basicConfigInstance->is_NoTrustScenario()) {
				// cont->nc = neighborCont.at(rand() % neighborCont.size()); // lte
				std::cout << "InNeighbor'sCache(randomNode: " << cont->nc << " )" << std::endl;
				// closeness				teste.test
				//std::cout << " closeness mean: " << contact->getMean(cont->nc,node) << std::endl;
				//std::cout << " closeness var: " << contact->getVariance(cont->nc,node) << std::endl;
				//std::cout << " closeness n: " << contact->getN(cont->nc,node) << std::endl;
				//
				//double closen = closeness(contact->getMean(cont->nc,node), contact->getVariance(cont->nc,node));
				//std::cout << " closeness: " << closen << std::endl;
			
				cont->addFreq(cont->nc, contentNumber);
				cont->addContent(node, contentNumber);
				// cont->addCache(node, cont->nc, contentNumber, true);
				trustStaticsInstance->newContent();
				cont->addCache(node, contentNumber);
				cont->cacheHit++;
				if (maliciousInstance->isMalicious(cont->nc) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
					std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
					trustInstance->badBehavior(node, cont->nc, 1);
					trustStaticsInstance->newBadContentTransmission();
				}
				else if (maliciousInstance->isMalicious(cont->nc) && !basicConfigInstance->maliciousSwitched(sTime)){
					std::cout << "Node (" << cont->nc << ") malicious - invalid content" << std::endl;
					//trustInstance->badBehavior(node, cont->nc, 1);
					trustStaticsInstance->newInvalidContent();
					trustStaticsInstance->newBadContentTransmission();
				} else {
					std::cout << "Node (" << cont->nc << ") not malicious - valid content" << std::endl;
					// trustInstance->goodBehavior(node, cont->nc, 1);
				}
				// std::cout << "Direct Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getDirectTrust(node, cont->nc) << std::endl;
				// std::cout << "Indirect Trust: (" << node << ") -> ("<< cont->nc << ") " << trustInstance->getIndirectTrust(node, cont->nc) << std::endl;
				// std::cout << "Trust value of Node (" << node << ") on Node (" << cont->nc << ") -> " << trustInstance->getTrust(node, cont->nc) << std::endl;				// register content obtained from friends
				if (sg->isEdge(node, cont->nc))
				{
					std::cout << "\t->FriendCache";
					cont->friendCache++;
				}
			}
		}
	}

	else
	{
		if (blockchainConn->checkContent(node, contentNumber)){
			std::cout << "NewContent " << contentNumber << std::endl;
			//std::cout<<"sumline_node node = "<<node<< " ";
			//std::cout<<std::endl;
			cont->addContent(node, contentNumber);
			cont->addCache(node, contentNumber);
			cont->cacheMiss++;
		}
	}
	std::cout << std::endl;
}

void exchangeKeys(int nodeSrc, int nodeRcv)
{
	//std::string comp;
	Key::getInstance()->pool_add = 1;
	std::string srcPubKey;

	//comp = Key::getInstance()->getKeyFromPool("ecdsa", nodeRcv, nodeSrc);
	
	srcPubKey = Key::getInstance()->getECDSA_publicKeyFromMap(nodeSrc);
	//std::cout << "Chave do Nó (" << nodeSrc << ") adicionada na pool do ";
	Key::getInstance()->put_pool("ecdsa", nodeRcv, nodeSrc, srcPubKey);
	//std::cout << "Nó (" << nodeRcv << ")" << std::endl;
	//Key::getInstance()->printNodeECDSAPool(nodeRcv);
	

	//comp = Key::getInstance()->getKeyFromPool("ecc", nodeRcv, nodeSrc);
	
	srcPubKey = Key::getInstance()->getECDSA_publicKeyFromMap(nodeSrc);
	//std::cout << "Chave do Nó (" << nodeSrc << ") adicionada na pool do ";
	Key::getInstance()->put_pool("ecc", nodeRcv, nodeSrc, srcPubKey);
	//std::cout << "Nó (" << nodeRcv << ")" << std::endl;
	//Key::getInstance()->printNodeECDSAPool(nodeRcv);
	

}

void clusterKeyShare(Cluster* cluster, NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface) {
	// Na implementação, as chaves serão compartilhadas como se fosse de forma direta entre os nós do cluster, 
	// mas na realidade devo simular com o NS3 os nós do cluster enviando suas chaves para o CH, o qual é responsável por
	// distribuílas entre os membros do cluster.
	//double d_tNext = (uint32_t) 1040 * 8 / static_cast<double> (DataRate("1Mbps").GetBitRate ());
	double d_tNext = 0.0;
	std::cout << "clusterKeyShare - " << Simulator::Now().GetSeconds() << std::endl;
	uint16_t port = 1024 + rand() % 20000;
	std::map<int, std::vector<int>>::iterator it;
	for (it = cluster->clusterMembers.begin(); it != cluster->clusterMembers.end(); it++){
		std::cout << "ClusterHead " << it->first << " NodeSendKey ";
		std::vector<int>::iterator vecIt;
		d_tNext = (uint32_t) 1040 * 8 / static_cast<double> (DataRate("5Mbps").GetBitRate ());//+ Simulator::Now().GetSeconds();
		//std::cout << "clusterKeyShare - Cluster Head " << it->first << " ";
		for (vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++) {
			std::vector<int>::iterator vecIt2;
			for (vecIt2 = it->second.begin(); vecIt2 != it->second.end(); vecIt2++) {
				if (*vecIt != *vecIt2)
					exchangeKeys(*vecIt, *vecIt2);
			}
			if (it->first != *vecIt)
			{
				// Simula o envio das chaves através do clusterHead
				Time tNext (Seconds (d_tNext));
				TypeId tid = TypeId::LookupByName ("ns3::TcpSocketFactory");
				Ptr<Socket> recvSink = Socket::CreateSocket (ueNodes.Get (*vecIt), tid);
				InetSocketAddress local = InetSocketAddress (Ipv4Address::GetAny (), port);
				recvSink->Bind (local);
				recvSink->SetRecvCallback (MakeCallback (&ReceivePacket));

				Ptr<Socket> source = Socket::CreateSocket (ueNodes.Get (it->first), tid);
				InetSocketAddress remote = InetSocketAddress (ueWiFiIface.GetAddress (*vecIt, 0), port);
				source->Connect (remote);
				Simulator::ScheduleNow(&GenerateTraffic, source, 10, 1, tNext);
			}
			//Simulator::Schedule (tNext, &sendKeysAndMetadata, ueNodes, ueWiFiIface, it->first, *vecIt);	
		}
		std::cout << " ok " << std::endl;
		Key::getInstance()->printNodeECDSAPool(it->first);
		//std::cout << std::endl;
	}
}


void verifyMetadata(int node, int contentNumber, Content *cont)
{
	/*
    Inserir verificação dos metadados aqui
  */
}

void sendKeysAndMetadata(NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface, int nodeSrc, int nodeRcv)
{

	NS_LOG_INFO("time = " << Simulator::Now().GetSeconds());
	//if (Key::getInstance()->pool_add != -1)
	//{
		double sTime = Simulator::Now().GetSeconds();
		uint16_t port = 1024 + rand() % 20000;

		//std::cout << sTime << " ";
		
		//std::string signature;
		//signature = ecdsa.SignMessage(Key::getInstance()->getECDSA_privateKeyFromMap(cont->nc), content->metadata);
		PacketSinkHelper packetSinkHelper("ns3::TcpSocketFactory", InetSocketAddress(Ipv4Address::GetAny(), port));
		ApplicationContainer packetSink = packetSinkHelper.Install(ueNodes.Get(nodeRcv));
		packetSink.Start(Seconds(sTime));
		packetSink.Stop(Seconds(10000));

		BulkSendHelper bulkSendHelper("ns3::TcpSocketFactory", InetSocketAddress(ueWiFiIface.GetAddress(nodeRcv), port));
		bulkSendHelper.SetAttribute ("SendSize", UintegerValue (1040));
		//bulkSendHelper.SetAttribute("MaxBytes", UintegerValue(132));
		ApplicationContainer bulkSend = bulkSendHelper.Install(ueNodes.Get(nodeSrc));
		bulkSend.Start(Seconds(sTime));
		bulkSend.Stop(Seconds(10000));

		//std::cout << "Node (" << nodeRcv << ")" << std::endl;
		//std::cout << "Metadata - Content " << contentNumber << " submitted to Node (" << nodeRcv << ")" << std::endl;
		//Key::getInstance()->printNodeECDSAPool(nodeRcv);
	//}
	//Key::getInstance()->pool_add = -1;
}



void nodeSendContent(NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface, Content *cont, int node, Graph *g)
{

	NS_LOG_INFO("time = " << Simulator::Now().GetSeconds());
	//	unsigned int count = 0;
	//	for (unsigned int i = 0; i < numberOfNodes; i++)
	//	{ 
	//		if (g->isEdge(node,i))
	//		{
	//			count++;
	//		}
	//	}
	// Compartilhamento reativo
	if (cont->nc != -1)
	{

		double sTime = Simulator::Now().GetSeconds();
		uint16_t port = 1024 + rand() % 20000;

		std::cout << "NodeSendContentAtTime = " << sTime << std::endl;
		//std::cout << "OldContenton port (" << port << ")" << std::endl;
		std::cout << "SolicitadoPor (" << node << ") -> ";

		//for (std::vector<int>::const_iterator i=neighbors.begin(); i!=neighbors.end(); ++i)
		//{
		//std::cout << "neighbor " << (*i) << std::endl;
		PacketSinkHelper packetSinkHelper("ns3::TcpSocketFactory", InetSocketAddress(Ipv4Address::GetAny(), port));
		ApplicationContainer packetSink = packetSinkHelper.Install(ueNodes.Get(node));
		packetSink.Start(Seconds(sTime));
		packetSink.Stop(Seconds(10000));
		// Vou adicionar mais 32 bytes para o envio da chave de 256 bits.
		BulkSendHelper bulkSendHelper("ns3::TcpSocketFactory", InetSocketAddress(ueWiFiIface.GetAddress(node), port));
		bulkSendHelper.SetAttribute("SendSize", UintegerValue(2048));
		bulkSendHelper.SetAttribute("MaxBytes", UintegerValue(10240000));
		ApplicationContainer bulkSend = bulkSendHelper.Install(ueNodes.Get(cont->nc));
		bulkSend.Start(Seconds(sTime));
		bulkSend.Stop(Seconds(10000));

		std::cout << "ServidoPor (" << cont->nc << ")\n" << std::endl;
		cont->isExpected(node, cont->nc);
		//cont->cacheHit++;//}
	}

	// Compartilhamento proativo
	if (cont->d2dCache != -1)
	{
		double sTime = Simulator::Now().GetSeconds();
		uint16_t port = 1024 + rand() % 20000;
		std::cout << "d2dCache - NodeSendContentAtTime = " << sTime << std::endl;
		std::cout << "d2dCache - OldContentOnPort (" << port << ")" << std::endl;
		std::cout << "d2dCache - SolicitadoPor (" << node << ")" << std::endl;

		PacketSinkHelper packetSinkHelper("ns3::TcpSocketFactory",
										  InetSocketAddress(Ipv4Address::GetAny(), port));
		ApplicationContainer packetSink = packetSinkHelper.Install(ueNodes.Get(cont->d2dCache));
		packetSink.Start(Seconds(sTime));
		packetSink.Stop(Seconds(10000));
		
		BulkSendHelper bulkSendHelper("ns3::TcpSocketFactory", InetSocketAddress(
																   ueWiFiIface.GetAddress(/*cont->d2dCache*/ node), port));
		bulkSendHelper.SetAttribute("SendSize", UintegerValue(2048));
		bulkSendHelper.SetAttribute("MaxBytes", UintegerValue(10240000));

		ApplicationContainer bulkSend = bulkSendHelper.Install(ueNodes.Get(cont->d2dSend));
		bulkSend.Start(Seconds(sTime));
		bulkSend.Stop(Seconds(10000));
		cont->d2dCache = -1;
	}
}



void internetSendContent(unsigned int node, NodeContainer ueNodes, Ptr<Node> remoteHost, Ipv4InterfaceContainer ueLteIface, Content *cont)
{
	if (cont->lc)
	{
		std::cout << "InternetSendContent - localCache" << std::endl;
		cont->lc = false;
		//cont->localCache++;
	}

	else if (cont->nc == -1)
	{

		double sTime = Simulator::Now().GetSeconds();
		uint16_t port = 1024 + rand() % 20000;

		std::cout << "InternetSendContentAtTime = " << sTime << std::endl;
		//std::cout << "New Content on port (" << port << ")" << std::endl;
		std::cout << "SolicitadoPor (" << node << ")" << std::endl;

		PacketSinkHelper packetSinkHelper("ns3::TcpSocketFactory",
										  InetSocketAddress(Ipv4Address::GetAny(), port));
		ApplicationContainer packetSink = packetSinkHelper.Install(ueNodes.Get(node));
		packetSink.Start(Seconds(sTime));
		packetSink.Stop(Seconds(10000));

		BulkSendHelper bulkSendHelper("ns3::TcpSocketFactory", InetSocketAddress(
																   ueLteIface.GetAddress(node), port));
		//bulkSendHelper.SetAttribute ("SendSize", UintegerValue (2048));
		bulkSendHelper.SetAttribute("MaxBytes", UintegerValue(10240000));
		ApplicationContainer bulkSend = bulkSendHelper.Install(remoteHost);
		bulkSend.Start(Seconds(sTime));
		bulkSend.Stop(Seconds(10000));
		std::cout << "\n";
		//cont->cacheMiss++;
	}

	if (cont->nextCache != -1)
	{

		double sTime = Simulator::Now().GetSeconds();
		uint16_t port = 1024 + rand() % 20000;

		std::cout << "proCache - InternetSendContentAtTime = " << sTime << std::endl;
		//std::cout << "proCache - New Content Cache port (" << port << ")" << std::endl;
		std::cout << "proCache - CacheOn  = " << cont->nextCache << std::endl;

		PacketSinkHelper packetSinkHelper("ns3::TcpSocketFactory",
										  InetSocketAddress(Ipv4Address::GetAny(), port));
		ApplicationContainer packetSink = packetSinkHelper.Install(ueNodes.Get(cont->nextCache));
		packetSink.Start(Seconds(sTime));
		packetSink.Stop(Seconds(10000));

		BulkSendHelper bulkSendHelper("ns3::TcpSocketFactory", InetSocketAddress(ueLteIface.GetAddress(node), port));

		//bulkSendHelper.SetAttribute ("SendSize", UintegerValue (2048));
		bulkSendHelper.SetAttribute("MaxBytes", UintegerValue(10240000));
		ApplicationContainer bulkSend = bulkSendHelper.Install(remoteHost);
		bulkSend.Start(Seconds(sTime));
		bulkSend.Stop(Seconds(10000));
		
		std::cout << "\n";
		//cont->cacheMiss++;
		cont->nextCache = -1;
	}
}

void nodePosition(NodeContainer ueNodes, int numberOfNodes) {
	for (int i = 0; i < numberOfNodes; i++) {
		Ptr<MobilityModel> mobility= ueNodes.Get(i)->GetObject<MobilityModel>();
		Vector pos = mobility->GetPosition ();
		std::cout  << "Node " << i+1 << ": POS: x=" << pos.x << ", y=" << pos.y << std::endl;
	}

}

void nodeVelocity(NodeContainer ueNodes, int numberOfNodes) {
	for (int i = 0; i < numberOfNodes; i++) {
		Ptr<MobilityModel> mobility= ueNodes.Get(i)->GetObject<MobilityModel>();
		Vector vel = mobility->GetVelocity ();
		std::cout  << "Node " << i+1 << ": VEL: x=" << vel.x << ", y=" << vel.y << std::endl;
	}

}

std::tuple<double, double> ueNodePosition(NodeContainer ueNodes, int node) {

	Ptr<MobilityModel> mobility= ueNodes.Get(node)->GetObject<MobilityModel>();
	Vector pos = mobility->GetPosition ();
	return std::make_tuple(pos.x, pos.y);
}

std::tuple<double, double> ueNodeVelocity(NodeContainer ueNodes, int node) {

	Ptr<MobilityModel> mobility= ueNodes.Get(node)->GetObject<MobilityModel>();
	Vector vel = mobility->GetVelocity ();
	return std::make_tuple(vel.x, vel.y);
}


void insertDegMap(std::map<double, std::vector<int>>* m, double id, int v) {
    std::map<double, std::vector<int>>::iterator it;
    it = m->find(id);
    if(it != m->end())
    {
        it->second.push_back(v);

    } else {
        std::vector<int> tmp;
        tmp.push_back(v);
        m->insert(std::pair<double, std::vector<int>>(id, tmp));
    }

}

void printDegMap(std::map<double, std::vector<int>> m, double id) {
	std::map<double, std::vector<int>>::iterator it;
	it = m.find(id);
	if (it != m.end()){
		std::vector<int>::iterator it1;
		std::cout << "    -> Degre ("  << id <<  "): ";
		for (it1 = it->second.begin(); it1 != it->second.end(); it1++) {
			std::cout << *it1 << " ";
		}
		std::cout << std::endl;
	} 
}


std::tuple<double, double> degMapPos(std::map<double, std::vector<int>> m, NodeContainer ueNodes, double id, int node) {
	std::map<double, std::vector<int>>::iterator it;
	it = m.find(id);
	if (it != m.end()){

		std::vector<int>::iterator it1;
		for (it1 = it->second.begin(); it1 != it->second.end(); it1++) {
			if(*it1 == node) {
				return ueNodePosition(ueNodes, node);
			}
		}
	}
	return std::make_tuple(-1, -1);
}

/*static void SinkRx (Ptr<const Packet> p, const Address &ad)
{
  std::cout << *p << std::endl;
}*/

void OnOffFlow(NodeContainer ueNodes, int src, int dest, int port, std::string flowRate, double start, double stop) {
	Ipv4InterfaceContainer host;

	PacketSinkHelper sink ("ns3::TcpSocketFactory", InetSocketAddress(Ipv4Address::GetAny (), port));
	ApplicationContainer apps = sink.Install(ueNodes.Get(dest));
	apps.Start(Seconds(start));
	apps.Stop(Seconds(stop));
	// then, print what the packet sink receives.
  	//Config::ConnectWithoutContext ("/NodeList/i/ApplicationList/0/$ns3::PacketSink/Rx", MakeCallback (&SinkRx));

	Ptr<Ipv4> ipv4node = ueNodes.Get(dest)->GetObject<Ipv4>();
	OnOffHelper onoff("ns3::TcpSocketFactory", InetSocketAddress(ipv4node->GetAddress(1,0).GetLocal(), port));
	auto rate = std::string(flowRate);

	onoff.SetAttribute ("DataRate", StringValue(rate));
	onoff.SetAttribute ("PacketSize",StringValue("20")); 
	onoff.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1.0]"));
	onoff.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=1.0]"));

	apps = onoff.Install(ueNodes.Get(src));
	apps.Start(Seconds(start));
	apps.Stop(Seconds(stop));
	Config::ConnectWithoutContext ("/NodeList/i/ApplicationList/0/$ns3::PacketSink/Rx", MakeCallback (&SinkRx));
	
}	


void addNodeMap(Ptr<Node> node, int id) {
	nodeMap.insert(std::pair<Ptr<Node>,int>(node,id));
}

int getNode (Ptr<Node> node) {
	std::map<Ptr<Node>,int>::iterator nodeIt;
	nodeIt = nodeMap.find(node);
	if(nodeIt != nodeMap.end())
		return nodeIt->second;
	else 
		return -1;	
}

void clusterHeadManage(NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface/*, Cluster* cluster*/, Contact* contact) {
	
	auto cluster = Cluster::getInstance();
	bool change = false;
	//std::cout << "clusterHeadSetUp = " << sTime << "s" << std::endl;
	//std::cout << "New Content on port (" << port << ")" << std::endl;
	std::map<int, std::vector<int>>::iterator it;
	//std::cout << "clusterHeadSetup - ";
	//double d_tNext = 0.0;
	for (it = cluster->clusterMembers.begin(); it != cluster->clusterMembers.end(); it++){
		//d_tNext = (uint32_t) 1040 * 8 / static_cast<double> (DataRate("5Mbps").GetBitRate ());//+ Simulator::Now().GetSeconds();
		for (std::vector<int>::iterator member = it->second.begin(); member != it->second.end(); member++) {
			if (*member != it->first && contact->isEdge(*member, it->first) && cluster->isNodeChangeCH(*member)){

				change = true;

				//Time tNext (Seconds (d_tNext));
				/*TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
				Ptr<Socket> recvSink = Socket::CreateSocket (ueNodes.Get (it->first), tid);
				InetSocketAddress local = InetSocketAddress (Ipv4Address::GetAny (), port);
				recvSink->Bind (local);
				recvSink->SetRecvCallback (MakeCallback (&ReceivePacket));

				Ptr<Socket> source = Socket::CreateSocket (ueNodes.Get (*member), tid);
				InetSocketAddress remote = InetSocketAddress (ueWiFiIface.GetAddress (it->first, 0), port);
				source->Connect (remote);
				Simulator::ScheduleNow(&GenerateTraffic, source, 56, 2, Seconds(1.0));*/
			}
			
		}
		//std::cout << std::endl;
	}

	if (change){
		//std::cout << " changes made ok!" << std::endl;
		//Time tNext (Seconds ((uint32_t) 1040 * 8 / static_cast<double> (DataRate("1Mbps").GetBitRate ())));
		Simulator::ScheduleNow(&clusterKeyShare, cluster, ueNodes, ueWiFiIface); 	
	}
	else
		std::cout << " ok!" << std::endl;
	
}

/*void ReceivePacket (Ptr<Socket> socket)
 {
	while (socket->Recv ())
	{
		//Ptr<Packet> packet = socket->Recv();
		//Ptr<Node> recvnode = socket->GetNode();
		//int recvNodeIndex = getNode(recvnode);
		//std::cout << "Node " << recvNodeIndex << " received one packet!" << std::cout;
		//NS_LOG_UNCOND ("Received one packet!");
	}
 }
static void GenerateTraffic (Ptr<Socket> socket, uint32_t pktSize, uint32_t pktCount, Time pktInterval )
{
	if (pktCount > 0)
	{

   		//Time inter = Seconds(interval);
		socket->Send (Create<Packet> (pktSize));
		DataRate dataRate = DataRate("1Mbps");
		Time tNext (Seconds (pktSize * 8 / static_cast<double> (dataRate.GetBitRate ())));
		Simulator::Schedule (pktInterval, &GenerateTraffic,
		                    socket, pktSize,pktCount - 1, tNext);
	}
	else
	{
		socket->Close ();
	}
}*/

void clusterHeadElection(NodeContainer ueNodes, Contact *contact/*, Cluster* cluster*/, Ipv4InterfaceContainer ueWiFiIface, TwoHopContact *tcontact){
	std::cout << "Cluster Head election" <<std::endl;
	//double maxDeg = 0.0;
	//double sTime = Simulator::Now().GetSeconds();
	std::vector<int> candidates;
	std::vector<double> degList;
	std::map <double, std::vector<int>> degMap;
	//double random_number;
	
	//nodeVelocity(ueNodes, numberOfNodes);
	auto cluster = Cluster::getInstance();

	cluster->toggleClusterChange(true, -1);
	// Recupera os graus de centralidade distintos em contact, e coloca no vector degList
	//std::cout << "Graus de centralidade distintos: ";
    for(int i = 0; i < numberOfNodes; i++) {
    	double dc = contact->degreeCent(i);
    	if (std::find(degList.begin(), degList.end(), dc) == degList.end()){
			degList.push_back(dc);
			//std::cout << dc << " ";
		}
    }
    //Ordena o grau de centralidade por ordem crescente
    std::cout << std::endl;
    sort(degList.begin(), degList.end(), std::greater<double>()); 

    //Constroi um mapa de graus de centralidade, onde temos um conjunto de nós para cada grau.
	std::cout << "Cluster: " << std::endl;
  	for (unsigned int i = 0; i < degList.size(); i++) {
  		double dc = degList[i]; 
		for (int j = 0; j < numberOfNodes; j++) {
			double nodeDc = contact->degreeCent(j);
			if (dc == nodeDc) {
				insertDegMap(&degMap, dc, j);
			}
		}
  	}

  	//int cont = 1;
  	for (unsigned int i = 0; i < degList.size(); i++) {
  		double dc = degList[i]; 
  		printDegMap(degMap, dc);
  	}
  	//int contCH = 0;
  	//std::cout << "Clusters da rodada: " << std::endl;
  	//cluster->printClusters();

  	for (unsigned int i = 0; i < degList.size(); i++) {
  		double dc = degList[i]; 
	  	std::map<double, std::vector<int>>::iterator it;
		it = degMap.find(dc);
		if (it != degMap.end()){	
			std::vector<int>::iterator it1;
			it1 = it->second.begin();
			//std::cout << "Size: " << it->second.size() << std::endl;
			if (it->second.size() == 1) {
				//std::cout << *it1 << " (if) ";
				candidates.push_back(*it1);
				
				if(cluster->setClusterHead(/*contCH,*/ *it1) == CH_CHANGE) {
					cluster->ch_change = true;	
				} 

				//contCH++;
			} else {
				std::map<int, double> disMap;

				std::cout << "treshold: " << cluster->getTreshold() << std::endl;

				std::vector<int>::iterator it2;
				std::vector<int> chCandidates = it->second;

				std::vector<int> indices(chCandidates.size());
				std::iota(indices.begin(), indices.end(), 0);
				std::random_shuffle(indices.begin(), indices.end());

				std::random_device rd;
			    std::mt19937 eng(rd());
			    std::uniform_int_distribution<> distr(2, (int) numberOfNodes/3);

				double numberOfCandidates = distr(eng);

				for (int i = 0; i < numberOfCandidates; i++) {
	 				if(cluster->setClusterHead(chCandidates[indices[i]]) == CH_CHANGE) {
	 					//std::cout << *it2 << " ";
						cluster->ch_change = true;
						cluster->setChNumber();
					} 		
				}
			}
		}
	}
	cluster->rnd++;
	cluster->setTreshold();

	
	clusterSetUp(ueNodes/*, cluster*/, ueWiFiIface, contact, tcontact);	
}
void testeClusterElection(Contact *contact, TwoHopContact *tcontact){
	int closest;
	auto cluster = Cluster::getInstance();
	for (int i = 0; i < numberOfNodes; ++i)
	{	
		closest = testeMaisProximo(i, contact, tcontact);
		if(cluster->setClusterHead(closest) == CH_CHANGE){
			cluster->ch_change = true;
			cluster->setChNumber();
		}
		int id;
		// Se o nó pertencer a algum cluster
		if (cluster->isMember(&id, i)) {
			//std::cout << "Node " << node << " membro do cluster " << id << std::endl;
			//Se o mesmo cluster foi escolhido como cluster
			cluster->removeClusterMember(id, i);
			if (id != closest) {
				
				cluster->NodeChangeCH(i);
			} 	
			//cluster->removeClusterHead(id);
		}
		cluster->addClusterMember(closest, i);
	}
	// Agora remove os clusterHeads dos outros clusters
	for (int j = 0; j < numberOfNodes; j++) {
		if (cluster->isClusterHead(j)) {
			int currentCluster;
			// Verifica se já é membro de algum cluster
			if(cluster->isMember(&currentCluster, j)){
				cluster->removeClusterMember(currentCluster, j);
				if (currentCluster != closest) {
					cluster->addMovingNode(j);
					cluster->NodeChangeCH(j);
				} 
				cluster->addClusterMember(j, j);
			}

		}
	}
	
	cluster->printClusters();
}

int testeMaisProximo(int node, Contact* contact, TwoHopContact *tcontact){
	auto posNode = ueNodePosition(ueNodes, node);

	int closest = -1;
	for (int j = 0; j < numberOfNodes; j++) {
		
		if(node != j){ 
			if (contact->isEdge(node,j) || tcontact->isEdge(node,j)) {
				//std::cout << "NodeCH " << j << ": ";
				if (closest != -1) {
					auto posNew = ueNodePosition(ueNodes, j);
					auto posOld = ueNodePosition(ueNodes, closest);
					
					double xDistance = std::get<0>(posNew) - std::get<0>(posNode);
					double yDistance = std::get<1>(posNew) - std::get<1>(posNode);
					double distance_new = sqrt(pow(xDistance, 2) + pow(yDistance, 2));

					xDistance = std::get<0>(posOld) - std::get<0>(posNode);
					yDistance = std::get<1>(posOld) - std::get<1>(posNode);
					double distance_old = sqrt(pow(xDistance, 2) + pow(yDistance, 2));

					if (distance_old > distance_new) {
						//std::cout << " closest" << std::endl;
						closest = j;
					} 
				}
				else if (closest == -1) {
					closest = j;
				}
			}

		} 
	}
	return closest;
		
}

void maisProximo(int node, Contact* contact, TwoHopContact *tcontact) {
	
	auto posNode = ueNodePosition(ueNodes, node);
	auto cluster = Cluster::getInstance();

	int closest = -1;
	for (int j = 0; j < numberOfNodes; j++) {
		if (cluster->isClusterHead(j)) {
			if(node != j && (contact->isEdge(node,j) || tcontact->isEdge(node,j))){
			//std::cout << "NodeCH " << j << ": ";
				if (closest != -1) {
					auto posNew = ueNodePosition(ueNodes, j);
					auto posOld = ueNodePosition(ueNodes, closest);
					
					double xDistance = std::get<0>(posNew) - std::get<0>(posNode);
					double yDistance = std::get<1>(posNew) - std::get<1>(posNode);
					double distance_new = sqrt(pow(xDistance, 2) + pow(yDistance, 2));

					xDistance = std::get<0>(posOld) - std::get<0>(posNode);
					yDistance = std::get<1>(posOld) - std::get<1>(posNode);
					double distance_old = sqrt(pow(xDistance, 2) + pow(yDistance, 2));

					if (distance_old > distance_new) {
						//std::cout << " closest" << std::endl;
						closest = j;
					} 
				}
				else if (closest == -1) {
					closest = j;
				}

			} else if (node == j) {
				if (closest != -1) {
					auto posNew = ueNodePosition(ueNodes, j);
					auto posOld = ueNodePosition(ueNodes, closest);
					
					double xDistance = std::get<0>(posNew) - std::get<0>(posNode);
					double yDistance = std::get<1>(posNew) - std::get<1>(posNode);
					double distance_new = sqrt(pow(xDistance, 2) + pow(yDistance, 2));

					xDistance = std::get<0>(posOld) - std::get<0>(posNode);
					yDistance = std::get<1>(posOld) - std::get<1>(posNode);
					double distance_old = sqrt(pow(xDistance, 2) + pow(yDistance, 2));

					if (distance_old > distance_new) {
						//std::cout << " closest" << std::endl;
						closest = j;
					} 
				}
				else if (closest == -1) {
					closest = j;
				}
			}
		}
		int id;
		// Se o nó pertencer a algum cluster
		if (cluster->isMember(&id, node)) {
			//std::cout << "Node " << node << " membro do cluster " << id << std::endl;
			//Se o mesmo cluster foi escolhido como cluster
			cluster->removeClusterMember(id, node);
			if (id != closest) {
				
				cluster->NodeChangeCH(node);
			} 	
			//cluster->removeClusterHead(id);
		}
		cluster->addClusterMember(closest, node);
	}
	
}
void verifyConnectivity(NodeContainer ueNodes, Contact *contact, Ipv4InterfaceContainer ueWiFiIface, TwoHopContact *tcontact) {
	double sTime = Simulator::Now().GetSeconds();
	auto cluster = Cluster::getInstance();
	std::cout << "\nverifyConnectivity - " << sTime << "s" << std::endl;
	// Todos os nós irão verificar na tabela de roteamento por um clusterhea a um ou dois hops de distância
	for (int node = 0; node < numberOfNodes; node++) {
		std::map<int, std::vector<int>>::iterator clusIt;
		// Calcula a posição do nó
		auto posNode = ueNodePosition(ueNodes, node);
		int currentCluster;
		// Verifica se já é membro de algum cluster
		if(cluster->isMember(&currentCluster, node)){
			
			int closest = currentCluster;
			
			//std::cout << "\t-> Node " << node << " cluster: " << currentCluster << std::endl;
	        for(clusIt = cluster->clusterMembers.begin(); clusIt != cluster->clusterMembers.end(); clusIt++) {
	        	int clusterID = clusIt->first; 
	        	if (clusterID != currentCluster && (contact->isEdge(node,clusterID) || tcontact->isEdge(node,clusterID))) {

					auto posNew = ueNodePosition(ueNodes, clusterID);
					auto posOld = ueNodePosition(ueNodes, currentCluster);
					
					double xDistance = std::get<0>(posNew) - std::get<0>(posNode);
					double yDistance = std::get<1>(posNew) - std::get<1>(posNode);
					double distance_new = sqrt(pow(xDistance, 2) + pow(yDistance, 2));

					xDistance = std::get<0>(posOld) - std::get<0>(posNode);
					yDistance = std::get<1>(posOld) - std::get<1>(posNode);
					double distance_old = sqrt(pow(xDistance, 2) + pow(yDistance, 2));

					if (distance_old > distance_new) {
						//std::cout << " closest" << std::endl;
						closest = clusterID;
					} 
				}	
			}    
			if (closest != currentCluster) {
		    	std::cout << "Node " << node << " cluster changed!" << std::endl;
		    	cluster->removeClusterMember(currentCluster, node);
				if (currentCluster != closest) {
					cluster->addMovingNode(node);
					cluster->NodeChangeCH(node);
				} 
				cluster->addClusterMember(closest, node);
		    }  
	    }  
	    
	    // Caso não esteja associado a nenhum cluster, deve procurar por um
	    else {
	    	std::cout << "Node is not assocated with a cluster head" << std::endl;
	    	int closest = -1;
	    	for(clusIt = cluster->clusterMembers.begin(); clusIt != cluster->clusterMembers.end(); clusIt++) {
	    		int clusterID = clusIt->first; 

	        	if (clusterID != node && (contact->isEdge(node,clusterID) || tcontact->isEdge(node,clusterID))) {
	        		if (closest != -1){
		        		auto posNew = ueNodePosition(ueNodes, clusterID);
						auto posOld = ueNodePosition(ueNodes, closest);
						
						double xDistance = std::get<0>(posNew) - std::get<0>(posNode);
						double yDistance = std::get<1>(posNew) - std::get<1>(posNode);
						double distance_new = sqrt(pow(xDistance, 2) + pow(yDistance, 2));

						xDistance = std::get<0>(posOld) - std::get<0>(posNode);
						yDistance = std::get<1>(posOld) - std::get<1>(posNode);
						double distance_old = sqrt(pow(xDistance, 2) + pow(yDistance, 2));

						if (distance_old > distance_new) {
							//std::cout << " closest" << std::endl;
							closest = clusterID;
						}
					} else {
						closest = clusterID;
					} 
	        	}
	    	}
	    	if (closest != -1) {
		    	//int id;
		    	std::cout << "Node " << node << " go to cluster "<< closest <<"!" << std::endl;
		   		//cluster->removeClusterMember(&id, node);
				cluster->addMovingNode(node);
				cluster->NodeChangeCH(node);
				
				cluster->addClusterMember(closest, node);
			}
	    }
	    
	    
	}
	cluster->printClusters();
}

void clusterSetUp(NodeContainer ueNodes/*, Cluster* cluster*/, Ipv4InterfaceContainer ueWiFiIface, Contact* contact, TwoHopContact *tcontact){

	// Quando os Cluster heads são eleitos, devem fazer um broadcast para os nós vizinhos com uma mensagem informando que são clusterheads
	// Ao receberem estas mensagens, os outros nós (não CHs) devem selecionar o CH mais adequado para se conectar (o mais próximo).
	//cluster->resetNodeChangeCHmap(false);
	//cluster->debugNodeCHangeMap();
	std::cout << "\nCluster Set Up - Fase 1" <<std::endl;
	Key::getInstance()->erase();
	auto cluster = Cluster::getInstance();
	//std::vector<int> leavingMembers;

	for (int i = 0; i < numberOfNodes; ++i)
	{	
		maisProximo(i, contact, tcontact);
	}

	if(cluster->removeEmptyClusters()){
		cluster->resetReintegrate();
	} else {
		cluster->updateHistory();
	}

	cluster->printClusters();
	
	//cluster->ch_change = false;
	//uint32_t packetSize = 1000; // bytes
   	//uint32_t numPackets = 3;

	//Time tNext (Seconds (pktSize * 8 / static_cast<double> (dataRate.GetBitRate ())));
	if (cluster->reintegrateCh) {
		clusterSetUpPhase2(ueNodes/*, cluster*/, ueWiFiIface, contact, tcontact);
	} else {
		if (cluster->nodeChangeMap_Change())
			clusterHeadManage(ueNodes, ueWiFiIface/*, cluster*/, contact);
	}
}


void clusterSetUpPhase2(NodeContainer ueNodes/*, Cluster* cluster*/, Ipv4InterfaceContainer ueWiFiIface, Contact* contact, TwoHopContact *tcontact){
	
	auto cluster = Cluster::getInstance();

	// Quando os Cluster heads são eleitos, devem fazer um broadcast para os nós vizinhos com uma mensagem informando que são clusterheads
	// Ao receberem estas mensagens, os outros nós (não CHs) devem selecionar o CH mais adequado para se conectar (o mais próximo).
	//double sTime = Simulator::Now().GetSeconds();

	
	std::cout << "Cluster Set Up - Fase 2" <<std::endl;
	//std::vector<int> leavingMembers;
	std::vector<int>::iterator it;
	for (it = cluster->CHtoReintegrate.begin(); it != cluster->CHtoReintegrate.end(); ++it)
	{	
		maisProximo(*it, contact, tcontact);
	}
	cluster->setChNumber();
	cluster->printClusters();
	cluster->updateHistory();
	
	//cluster->debugNodeCHangeMap();
	std::cout << std::endl;
	if (cluster->nodeChangeMap_Change())
		clusterHeadManage(ueNodes, ueWiFiIface/*, cluster*/, contact);

	cluster->resetReintegrate();
}

void Cache(uint16_t nextNode, uint16_t nextContent, Contact *contact, Graph *sg)
{
	Content* cont = Content::getInstance();
	Trust* trustInstance = Trust::getInstance(numberOfNodes);
	MaliciousNodes* maliciousInstance = MaliciousNodes::getInstance();
	TrustStatistics* trustStaticsInstance = TrustStatistics::getInstance(numberOfNodes);
	BasicConfig* basicConfigInstance = BasicConfig::getInstance();
	double sTime = Simulator::Now().GetSeconds();
	// BlockchainConnection* blockchainConn =  BlockchainConnection::getInstance();
	//std::cout << "Test Content Number " << contentNumber << std::endl;
	cont->nextCache = -1;
	Contact candidates(numberOfNodes);
	double influence = 0.0;
	int cache = -1;
	//int icont = 0;
	std::vector<unsigned int> influenceN;
	std::vector<unsigned int> neighborCont;

	//std::cout << "Cache - " << nextNode << " " << std::endl;

	
	// Create Candidates graph
	for (unsigned int i = 0; i < numberOfNodes; i++)
	{
		for (unsigned int j = 0; j < numberOfNodes; j++)
		{
			// Nós candidatos a receberem o conteúdo
			if ((contact->isEdge(i, j)) && sg->isEdge(nextNode, j))
			{
				candidates.addEdge(i, j);
			}
		}
	}
	

	for (unsigned int j = 0; j < numberOfNodes; j++)
	{
		//std::cout << "sumline " <<tg.sumLine(j)<< std::endl;
		// Recupera o maior grau de centralidade do grafo dos candidatos, armazenado na variável influence
		if (candidates.degreeCent(j) > influence)
		{
			//std::cout << "PASSEI AQUI" << std::endl;
			influence = candidates.degreeCent(j);
			std::cout << "\ninfluence " << influence << std::endl;
		}
	}

	for (unsigned int j = 0; j < numberOfNodes; j++)
	{
		//Adiciona todos os nós com o maior grau de centralidade encontrado no grafo dos candidatos
		if (candidates.degreeCent(j) == influence)
		{
			//std::cout << "\ninfluence " << influence << std::endl;
			influenceN.push_back(j);
			//std::cout << "Testing - " << "influenceN " <<j<< std::endl;
		}
	}
	//Caso o vetor com os nós de maior grau de centralidade não esteja vazio e a quantidade de elementos é superior a 1
	//  será selecionado um nó aleatório para receber o conteúdo
	if (!influenceN.empty() && influenceN.size() > 1)
	{
		//random
		cache = influenceN.at(rand() % influenceN.size());
		std::cout << "\ninfluenceN rand " << cache << std::endl;
		//closeness
		// Caso apenas um nó tenha o maior grau de centralidade possível, esse nó será o escolhido para fazer o cache
	}
	else if (influenceN.size() == 1)
	{
		cache = influenceN.at(0);
		std::cout << "influenceN 1 " << cache << std::endl;
	}

	// Add expected Cache
	// Se o nó delimitado para fazer o cache já tenha o conteúdo em cache, esse conteúdo não deve ser enviado

	if (cont->hasCache(cache, nextContent))
	{
		cont->lc = true;
		std::cout << "Cascade - Local Cache" << std::endl;
		cont->nextCache = -1;
	}
	// Caso algum nó tenha sido escolhido para fazer cache do conteúdo,
	else if (cache != -1)
	{

		std::cout << "Cascade - ContNeigh  ";
		unsigned int cn = 0;
		unsigned int ce = 0;

		for (unsigned int j = 0; j < numberOfNodes; j++)
		{
			// Verifica se algum vizinho (de acordo com o grafo d2d offline) tem o conteúdo desejado
			if ((contact->isEdge(cache, j)) && cont->hasCache(j, nextContent))
			{
				neighborCont.push_back(j);
				std::cout << j << " ";
				cn++;
			}
			if ((contact->isEdge(cache, j)) && sg->isEdge(j, nextNode))
			{
				cont->addExpected(cache, j);
				ce++;
			}
		}
		std::cout << std::endl;
		std::cout << "Cascade - ContNeigh Sum  " << cn << std::endl;

		// PROCACHE from LTE
		// Caso não tenha nós vizinhos, deve receber o conteúdo da estação base
		/* "O sistema irá realizar um sorteio de uma variável aleatória uniforme, e
		caso o valor sorteado seja menor que pcv o processo de cache proativo irá prosseguir. Se,
		eventualmente, algum de seus vizinhos já possuir o conteúdo em cache, o conteúdo será
		obtido deste vizinho via comunicação D2D, caso contrário o conteúdo será obtido via BS."*/
		if (neighborCont.empty())
		{
			if (((rand() % 1000) / 1000) < (1 - pow((1 - probV), ce)))
			{
				//if ((rand () % 25) < 4){
				cont->nextCache = cache;
				//std::cout<<"sumline_nextCache = "<<cont->nextCache<< " ";
				//std::cout<<std::endl;

				std::cout << "Cascade - cache proativo - ce " << ce << std::endl;

				// cont->addCache(cache, cont->d2dSend, nextContent, true); //
				cont->addCache(cache,nextContent);
				cont->proCache++;
				 
			}
			else
			{
				std::cout << "Cascade - problte fail" << std::endl;
			}
		}
		// PROCACHE FROM D2D
		//                else if (((rand () % 1000)/1000) < (1-pow(0.84,ce))){
		// Caso tenha algum nó vizinho, deve recuperar o conteúdo desse nó
		else
		{
			// Probabilidade de fazer cache proativo de um conteúdo compartilhado é denotado por pcv = 1-pow((1-probV),ce))
			/* "O sistema irá realizar um sorteio de uma variável aleatória uniforme, e
			caso o valor sorteado seja menor que pcv o processo de cache proativo irá prosseguir. Se,
			eventualmente, algum de seus vizinhos já possuir o conteúdo em cache, o conteúdo será
			obtido deste vizinho via comunicação D2D, caso contrário o conteúdo será obtido via BS."*/
			if (((rand() % 1000) / 1000) < (1 - pow((1 - probV), ce)))
			{
				cont->d2dSend = neighborCont.at(rand() % neighborCont.size()); // lte
				cont->d2dCache = cache;
				if (basicConfigInstance->is_completeScenario()){ // cond se verifica confiança
					if (maliciousInstance->isMalicious(cont->d2dSend) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
						std::cout << "On/Off - malicious" << std::endl;
						trustStaticsInstance->newInvalidContent();
					}
					else if (maliciousInstance->isMalicious(cont->d2dSend) && !basicConfigInstance->maliciousSwitched(sTime)){
						trustStaticsInstance->newInvalidContent();
					}
					if (trustInstance->is_trustable(cache, cont->d2dSend)){
						// Select the best node to serve the content
						std::cout << "Cascade - d2d cache " << std::endl;
						cont->addCache(cache, cont->d2dSend, nextContent, true);
						trustStaticsInstance->newContent();
						//cont->addTempCache(cache, nextContent, true);
						if (maliciousInstance->isMalicious(cont->d2dSend) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
							std::cout << "On/Off - malicious" << std::endl;
							std::cout << "Node (" << cont->d2dSend << ") malicious - invalid content" << std::endl;
							trustInstance->badBehavior(cache, cont->d2dSend, 1);
							trustStaticsInstance->newBadContentTransmission();
						}
						else if (maliciousInstance->isMalicious(cont->d2dSend) && !basicConfigInstance->maliciousSwitched(sTime)){
							std::cout << "Node (" << cont->d2dSend << ") malicious - invalid content" << std::endl;
							trustInstance->badBehavior(cache, cont->d2dSend, 1);
							trustStaticsInstance->newBadContentTransmission();
						} else {
							std::cout << "Node (" << cont->d2dSend << ") not malicious - valid content" << std::endl;
							trustInstance->goodBehavior(cache, cont->d2dSend, 1);
						}
						
						std::cout << "Direct Trust: (" << cache << ") -> ("<< cont->d2dSend << ") " << trustInstance->getDirectTrust(cache, cont->d2dSend) << std::endl;
						std::cout << "Indirect Trust: (" << cache << ") -> ("<< cont->d2dSend << ") " << trustInstance->getIndirectTrust(cache, cont->d2dSend) << std::endl;
						std::cout << "Trust value of Node (" << cache << ") on Node (" << cont->d2dSend << ") -> " << trustInstance->getTrust(cache, cont->d2dSend) << std::endl;
						
						cont->proCached2d++;
						if (neighborCont.size() == 1)
							std::cout << "Cascade - only d2d cache " << std::endl;
					} else {
						std::cout << "Not adding content " << nextContent << " on cache of node " << cache << std::endl;
						std::cout << "Node (" << cache << ") do not trust ("<< cont->d2dSend << ") " << trustInstance->getTrust(cache, cont->d2dSend) << std::endl;
					}
				} else if (basicConfigInstance->is_DTrustScenario()) {
					if (maliciousInstance->isMalicious(cont->d2dSend) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
						std::cout << "On/Off - malicious" << std::endl;
						trustStaticsInstance->newInvalidContent();
					}
					else if (maliciousInstance->isMalicious(cont->d2dSend) && !basicConfigInstance->maliciousSwitched(sTime)){
						trustStaticsInstance->newInvalidContent();
					}
					if (trustInstance->is_trustable(cache, cont->d2dSend)){
						// Select the best node to serve the content
						std::cout << "Cascade - d2d cache " << std::endl;

						cont->addCache(cache,nextContent);
						trustStaticsInstance->newContent();
						//cont->addTempCache(cache, nextContent, true);
						if (maliciousInstance->isMalicious(cont->d2dSend) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
							std::cout << "On/Off - malicious" << std::endl;
							std::cout << "Node (" << cont->d2dSend << ") malicious - invalid content" << std::endl;
							trustInstance->badBehavior(cache, cont->d2dSend, 1);
							trustStaticsInstance->newBadContentTransmission();
						}
						else if (maliciousInstance->isMalicious(cont->d2dSend) && !basicConfigInstance->maliciousSwitched(sTime)){
							std::cout << "Node (" << cont->d2dSend << ") malicious - invalid content" << std::endl;
							trustInstance->badBehavior(cache, cont->d2dSend, 1);
							trustStaticsInstance->newBadContentTransmission();
						} else {
							std::cout << "Node (" << cont->d2dSend << ") not malicious - valid content" << std::endl;
							trustInstance->goodBehavior(cache, cont->d2dSend, 1);
						}
						
						std::cout << "Direct Trust: (" << cache << ") -> ("<< cont->d2dSend << ") " << trustInstance->getDirectTrust(cache, cont->d2dSend) << std::endl;
						// std::cout << "Indirect Trust: (" << cache << ") -> ("<< cont->d2dSend << ") " << trustInstance->getIndirectTrust(cache, cont->d2dSend) << std::endl;
						// std::cout << "Trust value of Node (" << cache << ") on Node (" << cont->d2dSend << ") -> " << trustInstance->getTrust(cache, cont->d2dSend) << std::endl;
						
						cont->proCached2d++;
						if (neighborCont.size() == 1)
							std::cout << "Cascade - only d2d cache " << std::endl;
					} else {
						std::cout << "Not adding content " << nextContent << " on cache of node " << cache << std::endl;
						std::cout << "Node (" << cache << ") do not trust ("<< cont->d2dSend << ") " << trustInstance->getTrust(cache, cont->d2dSend) << std::endl;
					}
				} else if (basicConfigInstance->is_NoTrustScenario()) {
					// cont->d2dCache = cache;
					// Select the best node to serve the content
					std::cout << "Cascade - d2d cache " << std::endl;
					cont->addCache(cache,nextContent);
					trustStaticsInstance->newContent();
					//cont->addTempCache(cache, nextContent, true);
					if (maliciousInstance->isMalicious(cont->d2dSend) && basicConfigInstance->haveOnOff() && basicConfigInstance->maliciousOnOff()){
						std::cout << "Node (" << cont->d2dSend << ") malicious - invalid content" << std::endl;
						trustInstance->badBehavior(cache, cont->d2dSend, 1);
						trustStaticsInstance->newBadContentTransmission();
					}
					else if (maliciousInstance->isMalicious(cont->d2dSend) && !basicConfigInstance->maliciousSwitched(sTime)){
						std::cout << "Node (" << cont->d2dSend << ") malicious - invalid content" << std::endl;
						// trustInstance->badBehavior(cache, cont->d2dSend, 1);
						trustStaticsInstance->newInvalidContent();
						trustStaticsInstance->newBadContentTransmission();
					} else {
						std::cout << "Node (" << cont->d2dSend << ") not malicious - valid content" << std::endl;
						// trustInstance->goodBehavior(cache, cont->d2dSend, 1);
					}
					
					// std::cout << "Direct Trust: (" << cache << ") -> ("<< cont->d2dSend << ") " << trustInstance->getDirectTrust(cache, cont->d2dSend) << std::endl;
					// std::cout << "Indirect Trust: (" << cache << ") -> ("<< cont->d2dSend << ") " << trustInstance->getIndirectTrust(cache, cont->d2dSend) << std::endl;
					// std::cout << "Trust value of Node (" << cache << ") on Node (" << cont->d2dSend << ") -> " << trustInstance->getTrust(cache, cont->d2dSend) << std::endl;
					
					cont->proCached2d++;
					if (neighborCont.size() == 1)
						std::cout << "Cascade - only d2d cache " << std::endl;
				}
			}
		}
	}
	
}

//void Preload(int n, int zipfv, Content *cont, Graph *g){
//
//                int neighContent = 0;
//                for (unsigned int j = 0; j < numberOfNodes; j++)
//                {
//                        if ((g->isEdge(n,j)) && cont->hasContent(j,zipfv))
//                                neighContent += 1;
//                }
//                if (neighContent == 0){
//                        cont->addContent(n,zipfv);
//			std::cout << "Cache content " << zipfv << " in node " << n << std::endl;
//		}
//}

void showStatistic(Content *cont)
//showStatistic (Content *cont, ContactDistribution *contDist)
{
	TrustStatistics* trustStaticsInstance = TrustStatistics::getInstance(numberOfNodes);
	NetStats* netStatsInstance = NetStats::getInstance(numberOfNodes);
	trustStaticsInstance->print();
	netStatsInstance->print();
	std::cout << "Cache Hit: " << cont->cacheHit << std::endl;
	std::cout << "Cache Miss: " << cont->cacheMiss << std::endl;
	std::cout << "Local Cache: " << cont->localCache << std::endl;
	std::cout << "Proactive Cache: " << cont->proCache << std::endl;
	std::cout << "Proactive Cache D2D: " << cont->proCached2d << std::endl;
	std::cout << "Cache Eviction: " << cont->cacheEvict << std::endl;
	std::cout << "Friend Cache : " << cont->friendCache << std::endl;
	std::cout << "ExpetedT Cache: " << cont->expectedT << std::endl;
	std::cout << "ExpetedA Cache: " << cont->expectedA << std::endl;

	std::cout << "CacheHit rate: " << cont->getCHrate() << std::endl;
	std::cout << "offload rate: " << cont->getOFFrate() << std::endl;
	//        std::cout << "CacheHIT rate: " << (cont->cacheHit + cont->localCache)/(cont->cacheHit + cont->localCache + cont->cacheMiss) << std::endl;
	//        std::cout << "Offload rate: " << (cont->cacheHit + cont->localCache - cont->proCache)/(cont->cacheHit + cont->localCache + cont->cacheMiss) << std::endl;

	for (unsigned int j = 0; j < numberOfNodes; j++)
	{
		std::cout << "Views of node " << j << " " << cont->sumLine(j) << std::endl;
	}
	//std::cout << "Total Views: " << cont->printViews << std::endl;
	//contDist->printContacts();
}

//void clusterProcedure () {
//	std::vector<unsigned int> neighborCont;
//}

void testMatrix(Graph *g)
{
	g->printEdge();
}
