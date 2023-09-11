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
#include "scratch/graph.h"

#include <fstream>
#include <iostream>
#include <sstream>
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
//#include "ns3/gtk-config-store.h"



using namespace ns3;

/**
 * Sample simulation script for LTE+EPC. It instantiates several eNodeB,
 * attaches one UE per eNodeB starts a flow for each UE to  and from a remote host.
 * It also  starts yet another flow between each UE pair.
 */

NS_LOG_COMPONENT_DEFINE ("EpcFirstExample");

double 
lowerIncompleteGamma (double x, void *params);

double
generalisedPareto (double mu, double sigma, double E, double simTime);

double
Weibull (double lambda, double k);

double
closeness (double k, double teta);

void 
PrintRoutingTable (Ptr<Node> node); 

void
//RoutingTableToMatrix (NodeContainer ueNodes, int numberOfNodes, Graph *g);
RoutingTableToMatrix (NodeContainer ueNodes, int numberOfNodes, Graph *g, Contact *contact, ContactDistribution *contDist);

void
showStatistic (Content *cont);
//showStatistic (Content *cont, ContactDistribution *contDist);

void
testMatrix (Graph *g);

void
searchContent (Graph *sg, int node, int contentNumber, Content *cont, Contact *contact);

void
Cache(uint16_t nextNode, uint16_t nextContent, Content *cont, Contact *contact, Graph *sg);

//void 
//Preload(int n, int zipfv, Content *cont, Graph *g);

//searchContent (NodeContainer ueNodes, Ptr<Node> remoteHost, Ipv4InterfaceContainer ueWiFiIface, Ipv4InterfaceContainer ueLteIface, double simTime, int node, int contentNumber, Content *cont, Graph *g);

void
nodeSendContent (NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface, Content *cont, int node, Graph *g);

void
internetSendContent (unsigned int node,NodeContainer ueNodes, Ptr<Node> remoteHost, Ipv4InterfaceContainer ueLteIfaces, Content *cont);

/// Trace function for remaining energy at node.
void
RemainingEnergy (double oldValue, double remainingEnergy)
{
  if (Simulator::Now () .GetSeconds () > 1000.05)
  std::cout << Simulator::Now ().GetSeconds ()
            << "s Current remaining energy = " << remainingEnergy << "J" << std::endl;
}

/// Trace function for total energy consumption at node.
void
TotalEnergy (double oldValue, double totalEnergy)
{
  if (Simulator::Now () .GetSeconds () > 1000.05)
  std::cout << Simulator::Now ().GetSeconds ()
            << "s Total energy consumed by radio on node = " << totalEnergy << "J" << std::endl;
}

double probV = 0.32;

int
main (int argc, char *argv[])
{

  Time::SetResolution (Time::NS);
  //LogComponentEnable ("BulkSendApplication", LOG_LEVEL_INFO);
  //LogComponentEnable ("PacketSink", LOG_LEVEL_INFO);
  //LogComponentEnable ("OlsrRoutingProtocol", LOG_LEVEL_INFO);

// DEFAULT VALUES

   uint16_t numberOfNodes = 27;
   uint16_t contentPopulation = 100;
   double simTime = 10.1 ;
   double distance = 47.0;
   uint32_t seed = 1;
   double alpha = 10;
   unsigned int cacheSize = 10;
   std::string traceFile = "movements/100m/rpgm1.ns_movements";
   std::string socialGraphFile = "scratch/srsn.csv";

//   double interPacketInterval = 10;
   
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

//  cmd.AddValue("interPacketInterval", "Inter packet interval [ms])", interPacketInterval);
  cmd.Parse(argc, argv);

// Create node containers
  NodeContainer ueNodes;
  NodeContainer enbNodes;
  ueNodes.Create(numberOfNodes);
  enbNodes.Create(1);

// Create LTE
  Ptr<LteHelper> lteHelper = CreateObject<LteHelper> ();
  Ptr<PointToPointEpcHelper>  epcHelper = CreateObject<PointToPointEpcHelper> ();
  lteHelper->SetEpcHelper (epcHelper);

  ConfigStore inputConfig;
  inputConfig.ConfigureDefaults();

  // parse again so you can override default values from the command line
  cmd.Parse(argc, argv);

  Ptr<Node> pgw = epcHelper->GetPgwNode ();

//   Internet
//
//   // Create a single RemoteHost
  NodeContainer remoteHostContainer;
  remoteHostContainer.Create (1);
  Ptr<Node> remoteHost = remoteHostContainer.Get (0);
  InternetStackHelper internet;
  internet.Install (remoteHostContainer);
//
//  // Create the Internet
  PointToPointHelper p2ph;
  p2ph.SetDeviceAttribute ("DataRate", DataRateValue (DataRate ("100Gb/s")));
  p2ph.SetDeviceAttribute ("Mtu", UintegerValue (1500));
  p2ph.SetChannelAttribute ("Delay", TimeValue (Seconds (0.010)));
  NetDeviceContainer internetDevices = p2ph.Install (pgw, remoteHost);
  Ipv4AddressHelper ipv4h;
  ipv4h.SetBase ("1.0.0.0", "255.0.0.0");
  Ipv4InterfaceContainer internetIpIfaces = ipv4h.Assign (internetDevices);
//  // interface 0 is localhost, 1 is the p2p device
  Ipv4Address remoteHostAddr = internetIpIfaces.GetAddress (0);
//
  Ipv4StaticRoutingHelper ipv4RoutingHelper;
  Ptr<Ipv4StaticRouting> remoteHostStaticRouting = ipv4RoutingHelper.GetStaticRouting (remoteHost->GetObject<Ipv4> ());
  remoteHostStaticRouting->AddNetworkRouteTo (Ipv4Address ("7.0.0.0"), Ipv4Mask ("255.0.0.0"), Ipv4Address (remoteHostAddr), 1);

//
//  Internet fim
 
 // Mobility 
 //   ENB Mobility
 //
    MobilityHelper enbMobility;
    Ptr<ListPositionAllocator> positionAlloc = CreateObject<ListPositionAllocator> ();
    positionAlloc->Add (Vector ((500.0), (500.0), 0.0));
    enbMobility.SetPositionAllocator (positionAlloc);
    //enbMobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
    enbMobility.Install (enbNodes);
 
   Ns2MobilityHelper mobility = Ns2MobilityHelper (traceFile);
   mobility.Install ();

 // Mobility FIM


 // Install LTE Devices to the nodes
   NetDeviceContainer enbLteDevs = lteHelper->InstallEnbDevice (enbNodes);
   NetDeviceContainer ueLteDevs = lteHelper->InstallUeDevice (ueNodes);

 // Create WiFi
        YansWifiChannelHelper wifiChannel = YansWifiChannelHelper::Default ();
	//std::string phyMode ("DsssRate11Mbps");
	std::string phyMode ("ErpOfdmRate24Mbps");
// range	wifiChannel.SetPropagationDelay("ns3::ConstantSpeedPropagationDelayModel");
// range	wifiChannel.AddPropagationLoss("ns3::RangePropagationLossModel", "MaxRange",
// range			DoubleValue(70.0));
	Ptr<YansWifiChannel> channel = wifiChannel.Create ();
	LogDistancePropagationLossModel loss;
	loss.SetPathLossExponent (2.7);
	loss.SetReference (1, 46.6777);
	channel->SetPropagationLossModel (loss.GetObject<LogDistancePropagationLossModel> ());
	YansWifiPhyHelper wifiPhy = YansWifiPhyHelper::Default ();
 	wifiPhy.SetChannel (channel);
	wifiPhy.Set("EnergyDetectionThreshold", DoubleValue(-76.0)); //defulat val is -96dBm
	wifiPhy.Set ("TxPowerStart",DoubleValue (22));
	wifiPhy.Set ("TxPowerEnd", DoubleValue (22));
	//wifiPhy.Set ("CcaMode1Threshold", DoubleValue (-96 * 0.125892541));

	/*
	 * Wifi:
	 */
	WifiHelper wifi = WifiHelper::Default ();
	/*
	 * Set standard to 802.11g
	 */
	wifi.SetStandard (WIFI_PHY_STANDARD_80211g);
	/*
	 * Use constant rates for data and control transmissions
	 */
	wifi.SetRemoteStationManager ("ns3::ConstantRateWifiManager",
			/*
			 * The maximum number of retransmission attempts for an RTS. This
			 * value will not have any effect on some rate control algorithms.
			 * Set with class: ns3::UintegerValue
			 * Underlying type: uint32_t 0:4294967295
			 * Initial value: 7
			 */
			"MaxSsrc", UintegerValue (0),
                         /*
			 * Turn off RTS/CTS for frames below 2049 bytes.
			 * If the size of the data packet + LLC header + MAC header + FCS
			 * trailer is bigger than this value, we use an RTS/CTS handshake
			 * before sending the data, as per IEEE Std. 802.11-2007,
			 * Section 9.2.6. This value will not have any effect on some rate
			 * control algorithms.
			 * Set with class: ns3::UintegerValue
			 * Underlying type: uint32_t 0:4294967295
			 * Initial value: 2346
			 */
			"RtsCtsThreshold", UintegerValue (10000),
			/*
			 * Disable fragmentation for frames below 2049 bytes.
			 */
			"FragmentationThreshold", StringValue ("10000"),
			/*
			 * DataMode: The transmission mode to use for every data packet
			 * transmission.
			 * Set with class: WifiModeValue
			 * Underlying type: WifiMode
			 */
			"DataMode", StringValue (phyMode),
			/*
			 * ControlMode: The transmission mode to use for every control
			 * packet transmission.
			 * Set with class: WifiModeValue
			 * Underlying type: WifiMode
			 */
			"ControlMode", StringValue (phyMode),
			/*
			 * NonUnicastMode: Wifi mode used for non-unicast transmissions.
			 * Set with class: WifiModeValue
			 * Underlying type: WifiMode
			 */
			"NonUnicastMode",StringValue (phyMode));

        NqosWifiMacHelper wifiMac;
	wifiMac = NqosWifiMacHelper::Default ();

  // Instal WiFi devices to the nodes
    NetDeviceContainer ueWiFiDevs = wifi.Install (wifiPhy, wifiMac, ueNodes);

  /** Energy Model **/
  /***************************************************************************/
  /* energy source */
  BasicEnergySourceHelper basicSourceHelper;
  // configure energy source
  basicSourceHelper.Set ("BasicEnergySourceInitialEnergyJ", DoubleValue (2000.0));
  // install source
  EnergySourceContainer sources = basicSourceHelper.Install (ueNodes);
  /* device energy model */
  WifiRadioEnergyModelHelper radioEnergyHelper;
  // configure radio energy model
  radioEnergyHelper.SetTxCurrentModel("ns3::LinearWifiTxCurrentModel");
  //radioEnergyHelper.Set ("TxCurrentA", DoubleValue (0.380));
  //radioEnergyHelper.Set ("RxCurrentA", DoubleValue (0.313));
  // install device model
  DeviceEnergyModelContainer deviceModels = radioEnergyHelper.Install (ueWiFiDevs, sources);
  /***************************************************************************/

  /** connect trace sources **/
  /***************************************************************************/
  // all traces are connected to node 1 (Destination)
  // energy source
  for (int j = 0; j < 27; j++){
	  Ptr<BasicEnergySource> basicSourcePtr = DynamicCast<BasicEnergySource> (sources.Get (j));
	  basicSourcePtr->TraceConnectWithoutContext ("RemainingEnergy", MakeCallback (&RemainingEnergy));
	  // device energy model
	  Ptr<DeviceEnergyModel> basicRadioModelPtr =
	    basicSourcePtr->FindDeviceEnergyModels ("ns3::WifiRadioEnergyModel").Get (0);
	  NS_ASSERT (basicRadioModelPtr != 0);
	  basicRadioModelPtr->TraceConnectWithoutContext ("TotalEnergyConsumption", MakeCallback (&TotalEnergy));
  }

  /*
  * Create stack protocols:
  */
    OlsrHelper olsr;
    Ipv4ListRoutingHelper list;
    list.Add (ipv4RoutingHelper, 0);
    list.Add (olsr, 10);
    //InternetStackHelper stack;
    internet.SetRoutingHelper (list);
    internet.Install (ueNodes);
    //internet.Install (meshNodes);
 
    //Ipv4InterfaceContainer ueLteIface;
    //ueLteIface = epcHelper->AssignUeIpv4Address (NetDeviceContainer (ueLteDevs));
 
  // Set IP address for wifi devs
    Ipv4AddressHelper ipv4;
    NS_LOG_INFO ("Assign IP Addresses to the WiFi Devices.");
    ipv4.SetBase ("10.1.1.0", "255.255.255.0");
    Ipv4InterfaceContainer ueWiFiIface = ipv4.Assign (ueWiFiDevs);


// Install the IP stack on the UEs
//   internet.Install (ueNodes);
    Ipv4InterfaceContainer ueLteIface;
    ueLteIface = epcHelper->AssignUeIpv4Address (NetDeviceContainer (ueLteDevs));
 
 // Assign IP address to UEs
   for (uint32_t u = 0; u < ueNodes.GetN (); ++u)
   {
       Ptr<Node> ueNode = ueNodes.Get (u);
       // Set the default gateway for the UE
       Ptr<Ipv4StaticRouting> ueStaticRouting = ipv4RoutingHelper.GetStaticRouting (ueNode->GetObject<Ipv4> ());
       ueStaticRouting->SetDefaultRoute (epcHelper->GetUeDefaultGatewayAddress (), 2);
   }
 
 // Attach one UE per eNodeB
   for (uint16_t i = 0; i < numberOfNodes; i++)
   {
         lteHelper->Attach (ueLteDevs.Get(i), enbLteDevs.Get(0));
         // side effect: the default EPS bearer will be activated
   }

  lteHelper->EnableTraces ();
  // Uncomment to enable PCAP tracing
  p2ph.EnablePcapAll("lena-epc-first");

// Social Network Graph

    //std::string data("srsn.csv");
    std::ifstream in(socialGraphFile.c_str());
    if (!in.is_open()) return 1;

    typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;
    
    std::vector< std::string > vec;
    std::string line;

    Graph socialGraph (numberOfNodes);

    while (std::getline(in,line))
    {

        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());

        if (vec.size() < 2) continue;

        std::vector<int> vecNum;
        for (int i=0; i < 2; i++)
        {
                int num = std::atoi(vec.at(i).c_str());
                vecNum.push_back(num);
        }
        socialGraph.addEdge((vecNum.at(0)-1), (vecNum.at(1)-1));

//        std::copy(vec.begin(), vec.end(),
//             std::ostream_iterator<std::string>(std::cout, "|"));
//
//        std::cout << "\n----------------------" << std::endl;
    }
    in.close();

// Adjacency check 
// 
   Graph graph (numberOfNodes);
   Contact contact (numberOfNodes);
  ContactDistribution contDist(numberOfNodes);

   for (int update=2; update<simTime; update+=2)
   {
        Simulator::Schedule (Seconds (update), &RoutingTableToMatrix, ueNodes, numberOfNodes, &graph, &contact, &contDist);
// show adjacency Matrix	Simulator::Schedule (Seconds (2.5), &testMatrix, &socialGraph);
   }

/*
 * Pre Simulation content selection
 */

   Content cont (numberOfNodes, contentPopulation, cacheSize);
   View fv (numberOfNodes, contentPopulation);
   double shareRate[27] = {0};
   double shrT[27] = {0};
   double gp = 0.0;
   double drift = 5.0;
   double sTime = drift;
   double simTimeR = simTime*(alpha/5); // LIMITE PARA O SORTEIO DO TEMPO DE ACESSO
   srand (seed);
   SeedManager::SetSeed (seed);
   
   //uint32_t N = 100;
   //int neighContent;

   // Zhang ZIPf
   Ptr<ZipfRandomVariable> x = CreateObject<ZipfRandomVariable> ();
   x->SetAttribute ("N", IntegerValue (contentPopulation));
   x->SetAttribute ("Alpha", DoubleValue (1.9));

   // UNIFORM
   Ptr<UniformRandomVariable> shrt = CreateObject<UniformRandomVariable> ();
   shrt->SetAttribute ("Min", DoubleValue (0.0));
   shrt->SetAttribute ("Max", DoubleValue (1.0));

   gsl_rng * r_bn = gsl_rng_alloc (gsl_rng_rand);
   gsl_rng_set (r_bn, seed);

   // Content Schedule

   double gp_tmp;

   for (uint16_t n = 0; n < numberOfNodes; n++){

           int zeros = cont.zeros();

       	   shareRate[n] = generalisedPareto(-0.227, 0.305, -0.048, simTimeR);
           shrT[n] = shrt->GetValue ();
	   for (uint16_t c = 0; c < contentPopulation; c++){

                 if (cont.getViews(c) == 0){
                        int bin = gsl_ran_binomial(r_bn, (alpha/((n+1)*zeros)), 1);
                        if (bin == 1){
                                gp = generalisedPareto(-2654, 6315, -0.669, simTimeR);
				gp_tmp = gp;
                                fv.addView(n, c, gp+drift);
                                cont.addViews(c);
                                if (shareRate[n] > shrT[n]){
                                        for (uint16_t n_shr = 0; n_shr < numberOfNodes; n_shr++){
                                                if ((((rand () % 1000)/1000) < probV) && (socialGraph.isEdge(n,n_shr))){
					                gp = generalisedPareto(-2654, 6315, -0.669, simTimeR);
							fv.addView(n_shr, c, gp_tmp+gp);
                                                        cont.addViews(c);
                                                }
                                        }
                                }
                        }
                 }
                 else {
                        int pv = rand () % (n);
                        if(pv<=cont.getViews(c)){
                                gp = generalisedPareto(-2654, 6315, -0.669, simTimeR);
                                gp_tmp = gp;
                                fv.addView(n, c, gp+drift);
                                cont.addViews(c);
                                if (shareRate[n] > shrT[n]){
                                        for (uint16_t n_shr = 0; n_shr < numberOfNodes; n_shr++){
                                                if ((((rand () % 1000)/1000) < probV) && (socialGraph.isEdge(n,n_shr))){
                                                        gp = generalisedPareto(-2654, 6315, -0.669, simTimeR);
                                                        fv.addView(n_shr, c, gp_tmp+gp);
                                                        cont.addViews(c);
						}
					}
				}
                        }
                }
           }
    }

//// PRELOAD
//
//   for (uint16_t plc = 0; plc < floor(contentPopulation*0.5) ; plc++){
//        int randomN = (rand () % 27) ;
//        uint32_t v = x->GetInteger ();
//        int zipfv = (short) v;
//
//        Simulator::Schedule (Seconds (4.5), &Preload, randomN, zipfv, &cont, &graph);
//
//    }



// Aplication Schedule 
//
//	for (uint32_t i = 0; i < flows.size (); i++) {
	while ( sTime < simTime ) {

		if (fv.nextContent > contentPopulation)
			fv.firstView();
		else
			fv.nextView(fv.nextNode, fv.nextContent);
                  
		sTime = fv.getView(fv.nextNode, fv.nextContent);

                // Schedule content search on the neighbors cache

		Simulator::Schedule (Seconds (sTime), &searchContent, &socialGraph, fv.nextNode, fv.nextContent, &cont, &contact);

                // Send content through the BS

		Simulator::Schedule (Seconds (sTime), &internetSendContent, fv.nextNode, ueNodes, remoteHost, ueLteIface, &cont);

                // Send content through neighbors D2D communication

		Simulator::Schedule (Seconds (sTime), &nodeSendContent, ueNodes, ueWiFiIface, &cont, fv.nextNode, &graph);

        	//sTime+= incTime;

                // ProSoCaD CACHE - Share probability

                if (shareRate[fv.nextNode] > shrT[fv.nextNode]){
                        Simulator::Schedule (Seconds (sTime), &Cache, fv.nextNode, fv.nextContent, &cont, &contact, &socialGraph);
                }
		
	}


  Simulator::Schedule (Seconds (simTime), &showStatistic, &cont);


// Create route file
  Ptr<OutputStreamWrapper> routingStream = Create<OutputStreamWrapper> ("wifi-simple-adhoc-grid.routes", std::ios::out);
  olsr.PrintRoutingTableAllEvery (Seconds (2), routingStream);

  FlowMonitorHelper flowmon;
  Ptr<FlowMonitor> monitor = flowmon.Install (ueNodes);
  flowmon.Install (remoteHostContainer);

  //lteHelper->EnablePdcpTraces ();

  Simulator::Stop(Seconds(simTime));
  Simulator::Run();

  /*GtkConfigStore config;
  config.ConfigureAttributes();*/


  // Print per flow statistics
  Ptr<Node> nodeTemp = remoteHostContainer.Get (0); // Get pointer to ith node in container
  Ptr<Ipv4> ipv4Temp = nodeTemp->GetObject<Ipv4> (); // Get Ipv4 instance of the node
  Ipv4Address remoteIP = ipv4Temp->GetAddress (1, 0).GetLocal (); // Get Ipv4InterfaceAddress of xth interface

  monitor->CheckForLostPackets ();
  Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier> (flowmon.GetClassifier ());
  FlowMonitor::FlowStatsContainer stats = monitor->GetFlowStats ();
  double flowTime;
  double Throughput;
  for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin (); i != stats.end (); ++i)
    {
          Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow (i->first);
          std::cout << "Flow " << i->first  << " (" << t.sourceAddress << " -> " << t.destinationAddress << ")\n";
          std::cout << "  Tx Packets: " << i->second.txPackets << "\n";
          std::cout << "  Tx Bytes:   " << i->second.txBytes << "\n";
          std::cout << "  TxOffered:  " << i->second.txBytes * 8.0 / (i->second.timeLastRxPacket.GetSeconds () - i->second.timeFirstTxPacket.GetSeconds ()) / 1000 / 1000  << " Mbps\n";
          std::cout << "  Rx Packets: " << i->second.rxPackets << "\n";
          std::cout << "  Rx Bytes:   " << i->second.rxBytes << "\n";
	  flowTime = i->second.timeLastRxPacket.GetSeconds () - i->second.timeFirstTxPacket.GetSeconds (); 
	  Throughput = i->second.rxBytes * 8.0 / flowTime / 1000 / 1000;
          std::cout << "  Throughput: " << Throughput  << " Mbps\n";
          if (remoteIP == t.sourceAddress){
		double Pd = Throughput * 51.97 * flowTime;
		std::cout << " Promotion: " << "1210.7 mW\n";
		std::cout << " Download power: " << Pd << " mW\n"; 
		std::cout << " Flow duration: " << flowTime << " s\n";
	  }
          if (remoteIP == t.destinationAddress){
		double Pu = Throughput * 438.39 * flowTime;
                std::cout << " Upload power: " << Pu << " mW\n";
	  }
    }   

  //flowmon.SerializeToXmlFile("test.flowmon", true, true);

  Simulator::Destroy();
  return 0;

}

double
lowerIncompleteGamma (double x, void *params)
{
  // (t^(k-1))*(exp(-t))

  // The next line recovers alpha from the passed params pointer
  double alpha = *(double *) params;

  return ((pow(x,(alpha -1))) * (exp(-x)) );

}

double
generalisedPareto (double mu, double sigma, double E, double simTime)
{
  //double gp = 0.0;
  //gsl_rng * r1 = gsl_rng_alloc (gsl_rng_rand);
  //double U = gsl_rng_uniform_pos (r1);

  Ptr<UniformRandomVariable> u = CreateObject<UniformRandomVariable> ();
  u->SetAttribute ("Min", DoubleValue (0.0));
  u->SetAttribute ("Max", DoubleValue (1.0));

  double U = 0.0;
//  double gp = 0.0;
  double gp = std::numeric_limits<int>::max();

  while (gp > simTime){
  U =  u->GetValue ();
  gp = fabs(mu + sigma*(pow(U,-E) -1)/E);
  }
  return gp;
}

double
Weibull (double lambda, double k)
{
  //gsl_rng * r2 = gsl_rng_alloc (gsl_rng_rand);
  //double wb = gsl_ran_weibull(x, k, lambda);
 
  Ptr<WeibullRandomVariable> w = CreateObject<WeibullRandomVariable> ();
  w->SetAttribute ("Scale", DoubleValue (lambda));
  w->SetAttribute ("Shape", DoubleValue (k));
  
  double wb = w->GetValue ();


  return wb;
}

double
closeness (double k, double teta)
{
 gsl_integration_workspace *work_ptr
    = gsl_integration_workspace_alloc (1000);

  double lower_limit = 0;       /* lower limit a */
  double upper_limit = 50.0/teta;/* upper limit b */
  double abs_error = 1.0e-8;    /* to avoid round-off problems */
  double rel_error = 1.0e-8;    /* the result will usually be much better */
  double result;                /* the result from the integration */
  double error;                 /* the estimated error from the integration */

  double alpha = k;           // parameter in integrand

  gsl_function My_function;
  void *params_ptr = &alpha;

  My_function.function = &lowerIncompleteGamma;
  My_function.params = params_ptr;

  gsl_integration_qags (&My_function, lower_limit, upper_limit,
                        abs_error, rel_error, 1000, work_ptr, &result, &error);

// cout  std::cout.precision (18);
// cout  std::cout << "result          = " << result << std::endl;
// cout  std::cout << "estimated error = " << error << std::endl;

  double gamma = gsl_sf_gamma (alpha);

//  std::cout << "gamma          = " << gamma << std::endl;

  double wij = 1 - (result/gamma);

//  std::cout << "Closeness (wij) = " << wij << std::endl;

  gsl_integration_workspace_free (work_ptr);

  return (wij);
}

void
PrintRoutingTable (Ptr<Node> node)
{
    Ptr<WifiNetDevice> dev = node->GetObject<WifiNetDevice>();
    Ptr<ns3::Ipv4RoutingProtocol> rp(node-> GetObject <Ipv4> () -> GetRoutingProtocol());
    Ptr<ns3::olsr::RoutingProtocol> proto = Ipv4RoutingHelper::GetRouting <ns3::olsr::RoutingProtocol> (rp);
    std::vector<ns3::olsr::RoutingTableEntry> entry =  proto->GetRoutingTableEntries();


//        Ptr<ns3::olsr::RoutingProtocol> routing = node->GetObject<ns3::olsr::RoutingProtocol>();                
// Print  routing table entries for OLSR routing
  //      std::vector<ns3::olsr::RoutingTableEntry> entry = routing->GetRoutingTableEntries();
//        std::cout << "Routing table for device: " << Names::FindName(node) << std::endl;
        std::cout << "Routing table for device: " << node->GetId() << std::endl;
        std::cout << "DestinyAddress\t\tNextAddress\t\tInterface\t\tDistance\n";
        for (std::vector<ns3::olsr::RoutingTableEntry>::const_iterator i=entry.begin(); i!=entry.end(); i++)
        {
                std::cout << i->destAddr << "\t\t"
                                << i->nextAddr << "\t\t"
                                  << i->interface << "\t\t"
                                   << i->distance << std::endl;
        }
}

void  
//RoutingTableToMatrix (NodeContainer ueNodes, int numberOfNodes, Graph *g)
RoutingTableToMatrix (NodeContainer ueNodes, int numberOfNodes, Graph *g, Contact *contact, ContactDistribution *contDist)
{

//    std::cout << "Fazendo copia: ";
    Graph *gCopy = new Graph(*g); 
//    std::cout << "OK " << std::endl;
//    gCopy.printEdge();
//
    g->zeros();

     for (int nodeN=0; nodeN < numberOfNodes; nodeN++)
 //         std::cout << nodeNumber << std::endl;
       {
       Ptr<Node> node = ueNodes.Get(nodeN);
       for (int nodeC=0; nodeC < numberOfNodes; nodeC++)
       {
	  Ptr<Node> nC = ueNodes.Get(nodeC);

          Ptr<WifiNetDevice> dev = node->GetObject<WifiNetDevice>();
          Ptr<ns3::Ipv4RoutingProtocol> rp(node-> GetObject <Ipv4> () -> GetRoutingProtocol());
          Ptr<ns3::olsr::RoutingProtocol> proto = Ipv4RoutingHelper::GetRouting <ns3::olsr::RoutingProtocol> (rp);
          std::vector<ns3::olsr::RoutingTableEntry> entry =  proto->GetRoutingTableEntries();
      
              for (std::vector<ns3::olsr::RoutingTableEntry>::const_iterator i=entry.begin(); i!=entry.end(); i++)
              {
      		 int nodeInterface = nC->GetObject<Ipv4L3Protocol> () -> GetInterfaceForAddress(i->destAddr);
//		 std::cout << nodeInterface << std::endl;
		 if (nodeInterface != -1)
		 { 
			g->addEdge(node->GetId(), nC->GetId());
			
//			std::cout << i->destAddr << std::endl;
//                g->printEdge();

		 }
	      }	
  	}
     }

     for ( int nodeN=0; nodeN < numberOfNodes; nodeN++) {
       for (int nodeC=0; nodeC < numberOfNodes; nodeC++) {
	 if ((g->isEdge(nodeN,nodeC)) && !(gCopy->isEdge(nodeN,nodeC))){
	     contact->contactBegin(nodeN,nodeC,Simulator::Now().GetSeconds());
	     contDist->insertCont(nodeC);
         }
	 else if (!(g->isEdge(nodeN,nodeC)) && (gCopy->isEdge(nodeN,nodeC))){
             contDist->insertDur(Simulator::Now().GetSeconds() - contact->getlastContact(nodeN,nodeC));
             contact->contactEnd(nodeN,nodeC,Simulator::Now().GetSeconds());
         }
         else if ((g->isEdge(nodeN,nodeC)) && (gCopy->isEdge(nodeN,nodeC))){
	     contact->contactUpdate(nodeN,nodeC,Simulator::Now().GetSeconds());	
         }
       }
     }

    delete gCopy;
}


void
//searchContent (Graph *sg, int node, int contentNumber, Content *cont, Graph *g)
searchContent (Graph *sg, int node, int contentNumber, Content *cont, Contact *contact)
{
//                std::cout << "Test Content Number " << contentNumber << std::endl;
		cont->nc = -1;
		//std::cout << " teste " << cont->nc << std::endl;

                std::vector<int> neighborCont;
                neighborCont.reserve(27);

                if (cont->oldContent(contentNumber))
                {	
                        std::cout << "Old Content " << contentNumber << std::endl;

                        // Create a vector of neighbors with the content                        

		        for (unsigned int j = 0; j < 27; j++)
		        {
		                if ((contact->isEdge(node,j)) && cont->hasCache(j,contentNumber)) { 
		                        neighborCont.push_back(j);
		                }
				
		        }

			if (cont->hasCache(node, contentNumber)){
				cont->lc = true;
                                cont->addFreq(node, contentNumber);
				std::cout << "Local Cache" << std::endl;
				cont->localCache++;
			}
			else if (neighborCont.empty())
			{
				std::cout << "Not in neighbors cache" << std::endl;
				cont->cacheMiss++;
			}
			else if (neighborCont.size() == 1)
                        {
                                cont->nc = neighborCont.at(0);   			  // lte
                                std::cout << " Only Neighbor: " << cont->nc << std::endl;           
                                cont->addFreq(cont->nc, contentNumber);
	                        cont->addContent(node,contentNumber);
                                cont->addCache(node,contentNumber);
				cont->cacheHit++;
				if (sg->isEdge(node, cont->nc)){
					std::cout << "Friend Cache " << cont->nc << std::endl;
					cont->friendCache++;
				}
                        }
                        else
			{
                                // Select the best node to serve de content
// random
				cont->nc = neighborCont.at(rand() % neighborCont.size()); // lte
				std::cout << " random node: " << cont->nc << std::endl;
// closeness				teste.test
//                                std::cout << " closeness mean: " << contact->getMean(cont->nc,node) << std::endl;
//                                std::cout << " closeness var: " << contact->getVariance(cont->nc,node) << std::endl;
//                                std::cout << " closeness n: " << contact->getN(cont->nc,node) << std::endl;
//	
//				double closen = closeness(contact->getMean(cont->nc,node), contact->getVariance(cont->nc,node));
//                                std::cout << " closeness: " << closen << std::endl;

                                cont->addFreq(cont->nc, contentNumber);
                                cont->addContent(node,contentNumber);
                                cont->addCache(node,contentNumber);
				cont->cacheHit++;
				// register content obtained from friends
                                if (sg->isEdge(node, cont->nc)){
                                        std::cout << "Friend Cache " << cont->nc << std::endl;
                                        cont->friendCache++;
				}
			}
                }
                else
                {
                        std::cout << "New Content " << contentNumber << std::endl;
			//std::cout<<"sumline_node node = "<<node<< " ";
		         //std::cout<<std::endl;

       			cont->addContent(node,contentNumber);
                        cont->addCache(node,contentNumber);
			cont->cacheMiss++;
                }
}

void
nodeSendContent (NodeContainer ueNodes, Ipv4InterfaceContainer ueWiFiIface, Content *cont, int node, Graph *g)
{

	NS_LOG_INFO("time = " << Simulator::Now().GetSeconds());
//	unsigned int count = 0;
//	for (unsigned int i = 0; i < 27; i++)
//	{
//		if (g->isEdge(node,i))
//		{
//			count++;
//		}
//	}

	if (cont->nc != -1)
	{
		
	        double sTime =  Simulator::Now().GetSeconds();
        	uint16_t port =  1024 + rand () % 20000;

	        std::cout << "app start time = " << sTime << std::endl;
	        std::cout << "Old Content on port " << port << std::endl;
        	std::cout << "nó solicitando = " << node << std::endl;

		//for (std::vector<int>::const_iterator i=neighbors.begin(); i!=neighbors.end(); ++i)
		//{
			//std::cout << "neighbor " << (*i) << std::endl;
                PacketSinkHelper packetSinkHelper ("ns3::TcpSocketFactory",
                                InetSocketAddress (Ipv4Address::GetAny (), port));
                ApplicationContainer packetSink = packetSinkHelper.Install (ueNodes.Get (node));
                packetSink.Start (Seconds (sTime));
                packetSink.Stop (Seconds (10000));
//
                BulkSendHelper bulkSendHelper ("ns3::TcpSocketFactory", InetSocketAddress (
                                ueWiFiIface.GetAddress (node), port));
                //bulkSendHelper.SetAttribute ("SendSize", UintegerValue (2048));
                bulkSendHelper.SetAttribute ("MaxBytes", UintegerValue (10240000));
                ApplicationContainer bulkSend = bulkSendHelper.Install (ueNodes.Get (cont->nc));
                bulkSend.Start (Seconds (sTime));
                bulkSend.Stop (Seconds (10000));

		std::cout << "nó servindo = " << cont->nc << std::endl;
		cont->isExpected(node, cont->nc);
		//cont->cacheHit++;//}
	}

	if (cont->d2dCache != -1){
        	double sTime =  Simulator::Now().GetSeconds();
	        uint16_t port =  1024 + rand () % 20000;

                std::cout << "d2dCache - app start time = " << sTime << std::endl;
                std::cout << "d2dCache - Old Content on port " << port << std::endl;
                std::cout << "d2dCache - nó solicitando = " << node << std::endl;

		PacketSinkHelper packetSinkHelper ("ns3::TcpSocketFactory",
                InetSocketAddress (Ipv4Address::GetAny (), port));
                ApplicationContainer packetSink = packetSinkHelper.Install (ueNodes.Get (cont->d2dCache));
                packetSink.Start (Seconds (sTime));
                packetSink.Stop (Seconds (10000));
//
                BulkSendHelper bulkSendHelper ("ns3::TcpSocketFactory", InetSocketAddress (
                                ueWiFiIface.GetAddress (node), port));
                //bulkSendHelper.SetAttribute ("SendSize", UintegerValue (2048));
                bulkSendHelper.SetAttribute ("MaxBytes", UintegerValue (10240000));
                ApplicationContainer bulkSend = bulkSendHelper.Install (ueNodes.Get (cont->d2dSend));
                bulkSend.Start (Seconds (sTime));
                bulkSend.Stop (Seconds (10000));

		cont->d2dCache = -1;
	}
}

void
internetSendContent (unsigned int node, NodeContainer ueNodes, Ptr<Node> remoteHost, Ipv4InterfaceContainer  ueLteIface, Content *cont)
{
	if (cont->lc){
		std::cout << "ISC - local cache" << std::endl;
		cont->lc = false;
		//cont->localCache++;
	}

        else if (cont->nc == -1)
        {

        double sTime =  Simulator::Now().GetSeconds();
        uint16_t port = 1024 + rand () % 20000;

        std::cout << "app start time = " << sTime << std::endl;
	std::cout << "New Content on port " << port << std::endl;
        std::cout << "nó solicitando = " << node << std::endl;


       PacketSinkHelper packetSinkHelper ("ns3::TcpSocketFactory",
                             InetSocketAddress (Ipv4Address::GetAny (), port));
       ApplicationContainer packetSink = packetSinkHelper.Install (ueNodes.Get (node));
       packetSink.Start (Seconds (sTime));
       packetSink.Stop (Seconds (10000));

       BulkSendHelper bulkSendHelper ("ns3::TcpSocketFactory", InetSocketAddress (
                         ueLteIface.GetAddress (node), port));
         //bulkSendHelper.SetAttribute ("SendSize", UintegerValue (2048));
       bulkSendHelper.SetAttribute ("MaxBytes", UintegerValue (10240000));
       ApplicationContainer bulkSend = bulkSendHelper.Install (remoteHost);
       bulkSend.Start (Seconds (sTime));
       bulkSend.Stop (Seconds (10000));

	//cont->cacheMiss++;

	}

        if (cont->nextCache != -1)
        {

        double sTime =  Simulator::Now().GetSeconds();
        uint16_t port = 1024 + rand () % 20000;

        std::cout << "proCache - app start time = " << sTime << std::endl;
        std::cout << "proCache - New Content Cache port " << port << std::endl;
        std::cout << "proCache - Cache on  = " << cont->nextCache  << std::endl;


       PacketSinkHelper packetSinkHelper ("ns3::TcpSocketFactory",
                             InetSocketAddress (Ipv4Address::GetAny (), port));
       ApplicationContainer packetSink = packetSinkHelper.Install (ueNodes.Get (cont->nextCache));
       packetSink.Start (Seconds (sTime));
       packetSink.Stop (Seconds (10000));

       BulkSendHelper bulkSendHelper ("ns3::TcpSocketFactory", InetSocketAddress (
                         ueLteIface.GetAddress (node), port));
         //bulkSendHelper.SetAttribute ("SendSize", UintegerValue (2048));
       bulkSendHelper.SetAttribute ("MaxBytes", UintegerValue (10240000));
       ApplicationContainer bulkSend = bulkSendHelper.Install (remoteHost);
       bulkSend.Start (Seconds (sTime));
       bulkSend.Stop (Seconds (10000));

        //cont->cacheMiss++;
	cont->nextCache = -1;

        }

}

void
Cache (uint16_t nextNode, uint16_t nextContent, Content *cont, Contact *contact, Graph *sg)
{
//                std::cout << "Test Content Number " << contentNumber << std::endl;
	cont->nextCache = -1;
	Contact candidates (27);
	double influence = 0.0;
        int cache = -1;
	std::vector<unsigned int> influenceN;
	std::vector<unsigned int> neighborCont;

//        std::cout << "Cache - " << nextNode << " " << std::endl;
	
	// Create Candidates graph
	
        for (unsigned int i = 0; i < 27; i++)
        {
		for (unsigned int j = 0; j < 27; j++)
		{
		        if ((contact->isEdge(i,j)) && sg->isEdge(nextNode,j)) { 
		
		                candidates.addEdge(i,j);
		        }
		
		}
	}
	
	for (unsigned int j = 0; j < 27; j++)
        {
//                std::cout << "sumline " <<tg.sumLine(j)<< std::endl;
                if (candidates.degreeCent(j) > influence) {
		      influence = candidates.degreeCent(j);
                      std::cout << "influence " <<influence<< std::endl;
                }

        }
	
        for (unsigned int j = 0; j < 27; j++)
        {
                if (candidates.degreeCent(j) == influence) {
                	influenceN.push_back(j);
//  	                std::cout << "influenceN " <<j<< std::endl;
		}
        }
	
	if (!influenceN.empty() && influenceN.size() > 1){
//random
		cache = influenceN.at(rand() % influenceN.size());
                std::cout << "influenceN rand " <<cache<< std::endl;
//closeness
	}
	else if (influenceN.size() == 1){
		cache = influenceN.at(0);
                std::cout << "influenceN 1 " <<cache<< std::endl;
	}

// Add expected Cache
	

	if (cont->hasCache(cache, nextContent)){
                cont->lc = true;
                std::cout << "Cascade - Local Cache" << std::endl;
                cont->nextCache = -1;

	}
	else if (cache != -1){

                std::cout << "Cascade - ContNeigh  ";
                unsigned int cn = 0;
                unsigned int ce = 0;

                for (unsigned int j = 0; j < 27; j++)
                {
                        if ((contact->isEdge(cache,j)) && cont->hasCache(j,nextContent)) { 
  
                                neighborCont.push_back(j);
				std::cout << j << " ";
				cn++;
                        }

                        if ((contact->isEdge(cache,j)) && sg->isEdge(j,nextNode)) {
				cont->addExpected(cache, j);
				ce++;
			}
                }
                std::cout <<std::endl;
                std::cout << "Cascade - ContNeigh Sum  " << cn << std::endl;

// PROCACHE from LTE
                if (neighborCont.empty() )
                {
			if (((rand () % 1000)/1000) < (1-pow((1-probV),ce))){
//			if ((rand () % 25) < 4){
                		cont->nextCache = cache;
                		//std::cout<<"sumline_nextCache = "<<cont->nextCache<< " ";
                		//std::cout<<std::endl;

                		std::cout << "Cascade - cache proativo - ce "<<ce<<std::endl;

                		cont->addCache(cache,nextContent); //
                		cont->proCache++;       
			}                
			else std::cout << "Cascade - problte fail"<<std::endl; 
                }
// PROCACHE FROM D2D
//                else if (((rand () % 1000)/1000) < (1-pow(0.84,ce))){
	        else {
                	if (((rand () % 1000)/1000) < (1-pow((1-probV),ce))){
 		        	cont->d2dCache = cache;
 		        	// Select the best node to serve the content
 		        	cont->d2dSend = neighborCont.at(rand() % neighborCont.size()); // lte
 		        	std::cout << "Cascade - d2d cache " << std::endl;
 		        	cont->addCache(cache,nextContent);
 		        	cont->proCached2d++;
				if (neighborCont.size() == 1)
				        std::cout << "Cascade - only d2d cache " << std::endl;

 			}
		}
	}
}

//void Preload(int n, int zipfv, Content *cont, Graph *g){
//
//                int neighContent = 0;
//                for (unsigned int j = 0; j < 27; j++)
//                {
//                        if ((g->isEdge(n,j)) && cont->hasContent(j,zipfv))
//                                neighContent += 1;
//                }
//                if (neighContent == 0){
//                        cont->addContent(n,zipfv);
//			std::cout << "Cache content " << zipfv << " in node " << n << std::endl;
//		}
//}


void
showStatistic (Content *cont)
//showStatistic (Content *cont, ContactDistribution *contDist)
{
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

        for (unsigned int j = 0; j < 27; j++){
	        std::cout << "Views of node " << j << " " <<  cont->sumLine(j) << std::endl;
	}
        //std::cout << "Total Views: " << cont->printViews << std::endl;
        //contDist->printContacts();
}

void 
testMatrix (Graph *g)
{
    g->printEdge();
}	
