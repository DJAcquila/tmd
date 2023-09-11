#include <vector>
#include <utility>
#include <fstream>
#include <stack>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include "messages.h"
#include <tuple>
#include <random>
#include <map>
class BasicConfig {
private:
	int numberOfNodes;
	double trustTreshold;
	double maliciousNodesProportion;
	bool bTrust;
	bool dTrust;
	bool noTrust;
	bool onOffBehavior;
	double malProbability;
	int mb;
	static bool instanceFlag;
	static BasicConfig *configUtil;
public:
	static BasicConfig* getInstance();
	static void resetInstance();
	
	BasicConfig(){
		trustTreshold = 0.5;
		bTrust = false;
		dTrust = false;
		noTrust = false;
		onOffBehavior = false;
		mb = 1;
	}

	void variateMaliciousNodes(){
		mb = 0;
	}
	void notVariateMaliciousNodes(){
		mb = 1;
	}

	double getOnOffProbability(){
		return malProbability;
	}

	void setOnOffBehavior(double probability){
		onOffBehavior = true;
		malProbability = probability;
	}

	bool haveOnOff(){
		return onOffBehavior;
	}

	void setNumberOfNodes(int num){
		numberOfNodes = num;
	}

	int getNumberOfNodes(){
		return numberOfNodes;
	}

	void setTrustTreshold(double num){
		trustTreshold = num;
	}

	double getTrustTreshold(){
		return trustTreshold;
	}

	void setMaliciousNodesProportion(double num){
		maliciousNodesProportion = num;
	}

	double getMaliciousNodesProportion(){
		return maliciousNodesProportion;
	}

	void setCompleteScenario(){
		bTrust = true;
	}

	void setDTrustScenario(){
		dTrust = true;
	}

	void setNoTrustScenario(){
		noTrust = true;
	}

	bool is_completeScenario(){
		return bTrust;
	}

	bool is_DTrustScenario(){
		return dTrust;
	}

	bool is_NoTrustScenario(){
		return noTrust;
	}

	bool hasMaliciousVariance() {
		if (mb == 1) {
			return false;
		}
		return true;
	}
	bool maliciousSwitched(double time){
		if (hasMaliciousVariance()){
			if (time >= 0.0 && time < 40){
				return true;
			} 
		} 
		return false;
	}

	bool maliciousOnOff(){
		if (haveOnOff()){
			double valor = ((double) rand() / (RAND_MAX));
			if (valor < getOnOffProbability()){
				std::cout << "(On/Off) Valor sorteado (malicioso): " << valor << std::endl;
				return true;
			}
			std::cout << "(On/Off) Valor sorteado (NÃO malicioso): " << valor << std::endl;
		}
		return false;
	}

	void printScenarioInfos(){
		if (bTrust){
			std::cout << "=============== Complete Scenario ===============" << std::endl;
			std::cout << "\tNumber of nodes: " << numberOfNodes << std::endl;
			std::cout << "\tTrust Treshold: " << trustTreshold << std::endl;
			std::cout << "\tMalicious Nodes Proportion: " << maliciousNodesProportion << std::endl;
			if (hasMaliciousVariance()){
				std::cout << "\tMalicious Nodes Variate Behavior: yes"  << std::endl;
			} else {
				std::cout << "\tMalicious Nodes Variate Behavior: no"  << std::endl;
			}
			if (haveOnOff()){
				std::cout << "\tMalicious Nodes On/Off Behavior: yes"  << std::endl;
			} else {
				std::cout << "\tMalicious Nodes On/Off Behavior: no"  << std::endl;
			}
			std::cout << "=================================================" << std::endl;

		} else if (dTrust) {
			std::cout << "=============== Direct Trust Scenario ===============" << std::endl;
			std::cout << "\tNumber of nodes: " << numberOfNodes << std::endl;
			std::cout << "\tTrust Treshold: " << trustTreshold << std::endl;
			std::cout << "\tMalicious Nodes Proportion: " << maliciousNodesProportion << std::endl;
			if (hasMaliciousVariance()){
				std::cout << "\tMalicious Nodes Variate Behavior: yes"  << std::endl;
			} else {
				std::cout << "\tMalicious Nodes Variate Behavior: no"  << std::endl;
			}
			if (haveOnOff()){
				std::cout << "\tMalicious Nodes On/Off Behavior: yes"  << std::endl;
			} else {
				std::cout << "\tMalicious Nodes On/Off Behavior: no"  << std::endl;
			}
			std::cout << "=====================================================" << std::endl;
		} else if (noTrust) {
			std::cout << "=============== No Trust Scenario ===============" << std::endl;
			std::cout << "\tNumber of nodes: " << numberOfNodes << std::endl;
			std::cout << "\tTrust Treshold: " << trustTreshold << std::endl;
			std::cout << "\tMalicious Nodes Proportion: " << maliciousNodesProportion << std::endl;
			if (hasMaliciousVariance()){
				std::cout << "\tMalicious Nodes Variate Behavior: yes"  << std::endl;
			} else {
				std::cout << "\tMalicious Nodes Variate Behavior: no"  << std::endl;
			}
			if (haveOnOff()){
				std::cout << "\tMalicious Nodes On/Off Behavior: yes"  << std::endl;
			} else {
				std::cout << "\tMalicious Nodes On/Off Behavior: no"  << std::endl;
			}
			std::cout << "====================================================" << std::endl;
		} else {
			std::cout << "Error - scenario is not selected" << std::endl;
			exit(0);
		}
	}
};

bool BasicConfig::instanceFlag = false;
BasicConfig* BasicConfig::configUtil = NULL;

BasicConfig* BasicConfig::getInstance(){
	if (!instanceFlag) {
		configUtil = new BasicConfig();
		instanceFlag = true;
	}
	return configUtil;
}

void BasicConfig::resetInstance(){
	instanceFlag = false;
}

struct AvgDataStats{
	double timestamp;
	double avg;
};


class NetStats {
private:
	int m_nNodes;
	
	static bool instanceFlag;
	static NetStats *netStatsUtil;
	
	AvgDataStats* vazao_avg; //timestamp -> avg
	AvgDataStats* plr_avg; //timestamp -> avg
	AvgDataStats* delay_avg; //timestamp -> avg
	AvgDataStats* jitter_avg; //timestamp -> avg
	double vazao_total = 0.0;
	double plr_total = 0.0;

	int max_t_stamp;


public:
	static NetStats* getInstance(int nNodes);
	static NetStats* getInstance();
	NetStats(int nNodes){
		m_nNodes = nNodes;
		max_t_stamp = 110;
		
		this->vazao_avg = new AvgDataStats[max_t_stamp];
		for (int j = 0; j < max_t_stamp; j+=2) {
			init(j);
		}
		
		
		this->plr_avg = new AvgDataStats[max_t_stamp];
		for (int j = 0; j < max_t_stamp; j+=2) {
			init(j);
		}

		this->delay_avg = new AvgDataStats[max_t_stamp];
		for (int j = 0; j < max_t_stamp; j+=2) {
			init(j);
		}
		
		
		this->jitter_avg = new AvgDataStats[max_t_stamp];
		for (int j = 0; j < max_t_stamp; j+=2) {
			init(j);
		}
		
	}

	void init(int timestamp){
		this->vazao_avg[timestamp].timestamp = (double) timestamp;
		this->vazao_avg[timestamp].avg = 0.0;
	}

	void updateVazao(int timestamp, double avg){
		this->vazao_avg[timestamp].timestamp = (double) timestamp;
		this->vazao_avg[timestamp].avg = avg;
	}

	void updatePlr(int timestamp, double avg){
		this->plr_avg[timestamp].timestamp = (double) timestamp;
		this->plr_avg[timestamp].avg = avg;
	}

	void updateDelay(int timestamp, double avg){
		this->delay_avg[timestamp].timestamp = (double) timestamp;
		this->delay_avg[timestamp].avg = avg;
	}

	void updateJitter(int timestamp, double avg){
		this->jitter_avg[timestamp].timestamp = (double) timestamp;
		this->jitter_avg[timestamp].avg = avg;
	}

	double getVazao(int timestamp){
		return this->vazao_avg[timestamp].avg;
	}

	double getPlr(int timestamp){
		return this->plr_avg[timestamp].avg;
	}


	double getDelay(int timestamp){
		return this->delay_avg[timestamp].avg;
	}

	double getJitter(int timestamp){
		return this->jitter_avg[timestamp].avg;
	}

	void print(){
		
		for (int j = 0; j <= max_t_stamp; j+=2) {
			std::cout << "#--#Timestamp: (" << j << " s)" << std::endl;
			std::cout << "\t-> Average throughput: " << getVazao(j) << std::endl;
			std::cout << "\t-> Average plr: " << getPlr(j) << std::endl;
			std::cout << "\t-> Average delay: " << getDelay(j) << std::endl;
			std::cout << "\t-> Average jitter: " << getJitter(j) << std::endl;
		}
	}

	void print(int timestamp){
		
		std::cout << "#--#Timestamp: (" << timestamp << " s)" << std::endl;
		std::cout << "\t-> Average throughput: " << getVazao(timestamp) << " kbps" << std::endl;
		std::cout << "\t-> Average plr: " << getPlr(timestamp) << std::endl;
		std::cout << "\t-> Average delay: " << getDelay(timestamp) << std::endl;
		std::cout << "\t-> Average jitter: " << getJitter(timestamp) << std::endl;
	
	}
};

bool NetStats::instanceFlag = false;
NetStats* NetStats::netStatsUtil = NULL;

NetStats* NetStats::getInstance(int nNodes){
	if (!instanceFlag) {
		netStatsUtil = new NetStats(nNodes);
		instanceFlag = true;
	}
	return netStatsUtil;
}

NetStats* NetStats::getInstance(){
	if (instanceFlag) 
		return netStatsUtil;
	exit (EXIT_FAILURE);
}


//typedef struct duration{
//  int contCount;
//  double value;
//  struct duration *next;
//} duration;

struct Contacts{
	bool edge;
	int n;
	double mean;
	double variance;
	double M2;
	double lastContact;
};

//struct lteEnergy{
//  double upE;
//  double downE;
//  double activeTime; 
//};
 
class Graph {
private:
	
	bool** adjacencyMatrix;
	int vertexCount;
	int edgeCount;
public:
	Graph(int vertexCount) {
		this->vertexCount = vertexCount;
		this->edgeCount = 0;
		adjacencyMatrix = new bool*[vertexCount];
		for (int i = 0; i < vertexCount; i++) {
			adjacencyMatrix[i] = new bool[vertexCount];
			for (int j = 0; j < vertexCount; j++)
				adjacencyMatrix[i][j] = false;
		}
	}

	Graph(const Graph &m) {
		this->vertexCount = m.vertexCount;
		adjacencyMatrix = new bool*[vertexCount];
		for(int i = 0; i < vertexCount; i++) {
			adjacencyMatrix[i] = new bool[vertexCount];
			for(int j = 0; j < vertexCount; j++)
				adjacencyMatrix[i][j] =  m.adjacencyMatrix[i][j];
		}
	}

	std::vector<int> dominatingSet(){
		std::vector< std::vector<int> > g; 
		int numEdge = 0;
		g.resize(vertexCount);
		for (int i = 0; i < vertexCount; i++) {
			for (int j = 0; j < vertexCount; j++) {
				if (isEdge(i, j)){
					numEdge++;
					g[i].push_back(j);
				}
			}
		}

		bool visit[vertexCount + 1];

		memset(visit, 0, sizeof(visit));

		std::vector<int> set;
		for (int i = 0; i < vertexCount; ++i)
		{
			if(!visit[i]) {
				set.push_back(i);
				visit[i] = true;
				for(int j = 0; j < (int)g[i].size(); j++) {
					if(!visit[g[i][j]]) {
						visit[g[i][j]] = true;
						break;
					}
				}
			}
		}
		return set;
	}

	void addEdge(int i, int j) {
		if (i >= 0 && i < vertexCount && j > 0 && j < vertexCount) {
			adjacencyMatrix[i][j] = true;
			adjacencyMatrix[j][i] = true;
			this->edgeCount += 1;
		}
	}

	void removeEdge(int i, int j) {
		if (i >= 0 && i < vertexCount && j > 0 && j < vertexCount) {
			adjacencyMatrix[i][j] = false;
			adjacencyMatrix[j][i] = false;
			this->edgeCount -= 1;
		}
	}

	bool isEdge(int i, int j) {
		if (i >= 0 && i < vertexCount && j > 0 && j < vertexCount)
			return adjacencyMatrix[i][j];
		else
			return false;
	}

	void zeros () {
		for(int x=0;x<vertexCount;x++)  
		{
			for(int y=0;y<vertexCount;y++)  
			{

			// set current element of the array as false
				adjacencyMatrix[x][y] = false;
			}
		}
	}

	void printEdge () {
		for(int x=0;x<vertexCount;x++)  
		{
			for(int y=0;y<vertexCount;y++)  
			{
				std::cout<<adjacencyMatrix[x][y] << " "; 
			}
			std::cout<<std::endl;  
		}

	}

	unsigned int sumLine (unsigned int x) {
		unsigned int sum = 0;
		for(int y=0;y<vertexCount;y++)
		{
			if (adjacencyMatrix[x][y])
			{
				sum = sum + 1;
			}
		}
		return sum;    
	}

	~Graph() {
		for (int i = 0; i < vertexCount; i++)
			delete[] adjacencyMatrix[i];
		delete[] adjacencyMatrix;
	}
};

class TwoHopGraph {
private:
	
	bool** twoHopAdjacencyMatrix;
	int vertexCount;
	int edgeCount;
public:
	TwoHopGraph(int vertexCount) {
		this->vertexCount = vertexCount;
		this->edgeCount = 0;
		twoHopAdjacencyMatrix = new bool*[vertexCount];
		for (int i = 0; i < vertexCount; i++) {
			twoHopAdjacencyMatrix[i] = new bool[vertexCount];
			for (int j = 0; j < vertexCount; j++)
				twoHopAdjacencyMatrix[i][j] = false;
		}
	}

	TwoHopGraph(const TwoHopGraph &m) {
		this->vertexCount = m.vertexCount;
		twoHopAdjacencyMatrix = new bool*[vertexCount];
		for(int i = 0; i < vertexCount; i++) {
			twoHopAdjacencyMatrix[i] = new bool[vertexCount];
			for(int j = 0; j < vertexCount; j++)
				twoHopAdjacencyMatrix[i][j] =  m.twoHopAdjacencyMatrix[i][j];
		}
	}


	void addEdge(int i, int j) {
		if (i >= 0 && i < vertexCount && j > 0 && j < vertexCount) {
			twoHopAdjacencyMatrix[i][j] = true;
			twoHopAdjacencyMatrix[j][i] = true;
			this->edgeCount += 1;
		}
	}

	void removeEdge(int i, int j) {
		if (i >= 0 && i < vertexCount && j > 0 && j < vertexCount) {
			twoHopAdjacencyMatrix[i][j] = false;
			twoHopAdjacencyMatrix[j][i] = false;
			this->edgeCount -= 1;
		}
	}

	bool isEdge(int i, int j) {
		if (i >= 0 && i < vertexCount && j > 0 && j < vertexCount)
			return twoHopAdjacencyMatrix[i][j];
		else
			return false;
	}

	void zeros () {
		for(int x=0;x<vertexCount;x++)  
		{
			for(int y=0;y<vertexCount;y++)  
			{

			// set current element of the array as false
				twoHopAdjacencyMatrix[x][y] = false;
			}
		}
	}

	void printEdge () {
		for(int x=0;x<vertexCount;x++)  
		{
			for(int y=0;y<vertexCount;y++)  
			{
				std::cout<<twoHopAdjacencyMatrix[x][y] << " "; 
			}
			std::cout<<std::endl;  
		}

	}

	unsigned int sumLine (unsigned int x) {
		unsigned int sum = 0;
		for(int y=0;y<vertexCount;y++)
		{
			if (twoHopAdjacencyMatrix[x][y])
			{
				sum = sum + 1;
			}
		}
		return sum;    
	}

	~TwoHopGraph() {
		for (int i = 0; i < vertexCount; i++)
			delete[] twoHopAdjacencyMatrix[i];
		delete[] twoHopAdjacencyMatrix;
	}
};

class TwoHopContactDistribution {

private:
	std::vector<double> duration;
	unsigned int* count;
	unsigned int nodes;
public:
	TwoHopContactDistribution( unsigned int nodes){
		duration.reserve(nodes);
		this->nodes = nodes;
		count = new unsigned int[nodes];
		for (unsigned int i = 0; i < nodes; i++) {
			count[i] = 0;
		}
	}

	TwoHopContactDistribution(const TwoHopContactDistribution &m){
		//this->duration = m.duration;
		this->nodes = m.nodes;
		count = new unsigned int[nodes];
		for (unsigned int i = 0; i < nodes; i++) {
			count[i] = m.count[i];
		}
	}

	void insertDur(double val){
		duration.push_back(val);
	}

	void insertCont(unsigned int n){
		count[n]++;
	}

	void printContacts(){
		std::cout << "Contats Duration" << std::endl;
		for (std::vector<double>::const_iterator i = duration.begin(); i != duration.end(); ++i)
			std::cout << *i << std::endl;

		std::cout << "Contats Count" << std::endl;
		for (unsigned int n = 0; n < nodes; n++) {
			std::cout << count[n] << std::endl;
		}
	}
};

class ContactDistribution {

private:
	std::vector<double> duration;
	unsigned int* count;
	unsigned int nodes;
public:
	ContactDistribution( unsigned int nodes){
		duration.reserve(nodes);
		this->nodes = nodes;
		count = new unsigned int[nodes];
		for (unsigned int i = 0; i < nodes; i++) {
			count[i] = 0;
		}
	}

	ContactDistribution(const ContactDistribution &m){
		//this->duration = m.duration;
		this->nodes = m.nodes;
		count = new unsigned int[nodes];
		for (unsigned int i = 0; i < nodes; i++) {
			count[i] = m.count[i];
		}
	}

	void insertDur(double val){
		duration.push_back(val);
	}

	void insertCont(unsigned int n){
		count[n]++;
	}

	void printContacts(){
		std::cout << "Contats Duration" << std::endl;
		for (std::vector<double>::const_iterator i = duration.begin(); i != duration.end(); ++i)
			std::cout << *i << std::endl;

		std::cout << "Contats Count" << std::endl;
		for (unsigned int n = 0; n < nodes; n++) {
			std::cout << count[n] << std::endl;
		}
	}
};

class Contact {

private:
	Contacts** adjacencyMatrix;
	int vertexCount;
public:
	Contact(int vertexCount) {
		this->vertexCount = vertexCount;
		adjacencyMatrix = new Contacts*[vertexCount];
		for (int i = 0; i < vertexCount; i++) {
			adjacencyMatrix[i] = new Contacts[vertexCount];
			for (int j = 0; j < vertexCount; j++){
				adjacencyMatrix[i][j].edge = false;
				adjacencyMatrix[i][j].n = 0;
				adjacencyMatrix[i][j].mean = 0.0;
				adjacencyMatrix[i][j].variance = 0.0;
				adjacencyMatrix[i][j].M2 = 0.0;
				adjacencyMatrix[i][j].lastContact = 0.0; 
			}
		}
	}

//  Contact(const Contact &m) {
//      this->vertexCount = m.vertexCount;
//      adjacencyMatrix = new Contacts*[vertexCount];
//      for(int i = 0; i < vertexCount; i++) {
//          adjacencyMatrix[i] = new Contacts[vertexCount];
//          for(int j = 0; j < vertexCount; j++){
//              adjacencyMatrix[i][j].n =  m.adjacencyMatrix[i][j].n;
//              adjacencyMatrix[i][j].mean =  m.adjacencyMatrix[i][j].mean;
//              adjacencyMatrix[i][j].M2 =  m.adjacencyMatrix[i][j].M2;
//              adjacencyMatrix[i][j].lastContact =  m.adjacencyMatrix[i][j].lastContact;
//          }
//      }
//  }

	void addEdge(int i, int j){
		adjacencyMatrix[i][j].edge = true;
	}
	

	void contactBegin(int i, int j, double time){
		adjacencyMatrix[i][j].edge = true; 
		adjacencyMatrix[i][j].n += 1;
		adjacencyMatrix[i][j].lastContact = time;
		//
		adjacencyMatrix[j][i].edge = true;
		adjacencyMatrix[j][i].n += 1;
		adjacencyMatrix[j][i].lastContact = time;
	}


	void contactEnd(int i, int j, double time){
		adjacencyMatrix[i][j].edge = false;
		double x = (time - adjacencyMatrix[i][j].lastContact);
		double delta = x - adjacencyMatrix[i][j].mean;
		adjacencyMatrix[i][j].mean += delta/adjacencyMatrix[i][j].n;
		adjacencyMatrix[i][j].M2 += delta*(x - adjacencyMatrix[i][j].mean);
		if (adjacencyMatrix[i][j].n < 2) {
			adjacencyMatrix[i][j].variance = 0.1;
		}
		else {
			adjacencyMatrix[i][j].variance = adjacencyMatrix[i][j].M2/(adjacencyMatrix[i][j].n - 1);
		}
		//
		adjacencyMatrix[j][i].edge = false;
		adjacencyMatrix[j][i].mean = adjacencyMatrix[i][j].mean;
		adjacencyMatrix[j][i].M2 = adjacencyMatrix[i][j].M2;
		adjacencyMatrix[j][i].variance = adjacencyMatrix[i][j].variance;
	}

	void contactUpdate(int i, int j, double time){
		double x = (time - adjacencyMatrix[i][j].lastContact);
		double delta = x - adjacencyMatrix[i][j].mean;
		adjacencyMatrix[i][j].mean += delta/adjacencyMatrix[i][j].n;
		//adjacencyMatrix[i][j].M2 += delta*(x - adjacencyMatrix[i][j].mean);
		//if (adjacencyMatrix[i][j].n < 2) {
		//  adjacencyMatrix[i][j].variance = 0;
		//}
		//else {
		//  adjacencyMatrix[i][j].variance = adjacencyMatrix[i][j].M2/(adjacencyMatrix[i][j].n - 1);
		//}
		//
		adjacencyMatrix[j][i].mean = adjacencyMatrix[i][j].mean;
		adjacencyMatrix[j][i].M2 = adjacencyMatrix[i][j].M2;
		//adjacencyMatrix[j][i].variance = adjacencyMatrix[i][j].variance;

	}
	double degreeCent (int i) {
		double sum = 0;
		for(int j=0;j<vertexCount;j++)
		{
			if (adjacencyMatrix[i][j].edge)
			{
				sum = sum + adjacencyMatrix[i][j].mean;
			}
		}
		return sum;
	}

	bool isEdge (int i, int j){
		return adjacencyMatrix[i][j].edge;
	}

	std::vector<int> dominatingSet(){
		std::vector< std::vector<int> > g; 
		int numEdge = 0;
		g.resize(vertexCount);
		for (int i = 0; i < vertexCount; i++) {
			for (int j = 0; j < vertexCount; j++) {
				if (isEdge(i, j)){
					numEdge++;
					g[i].push_back(j);
				}
			}
		}

		bool visit[vertexCount + 1];

		memset(visit, 0, sizeof(visit));

		std::vector<int> set;
		for (int i = 0; i < vertexCount; ++i)
		{
			if(!visit[i]) {
				set.push_back(i);
				visit[i] = true;
				for(int j = 0; j < (int)g[i].size(); j++) {
					if(!visit[g[i][j]]) {
						visit[g[i][j]] = true;
						break;
					}
				}
			}
		}
		return set;
	}
	double getlastContact (int i, int j){
		return adjacencyMatrix[i][j].lastContact;
	}

	double getN (int i, int j){
		return adjacencyMatrix[i][j].n;
	}

	double getMean (int i, int j){
		//double mean = adjacencyMatrix[i][j].mean;
		//double x = (time - adjacencyMatrix[i][j].lastContact);
		//double delta = x - mean; 
		//mean += delta/adjacencyMatrix[i][j].n;
		return adjacencyMatrix[i][j].mean; 
	}

	double getVariance (int i, int j){
		return adjacencyMatrix[i][j].variance;
	}

	void printContacts(int i, int j){
		std::cout << "\nContact - Contact Nodes ("<< i << ", " << j <<")" << std::endl;
		std::cout << "edge = " << adjacencyMatrix[i][j].edge << '\n';
		std::cout << "n = " << adjacencyMatrix[i][j].n << '\n';
		std::cout << "mean = " << adjacencyMatrix[i][j].mean << '\n';
		std::cout << "M2 = " << adjacencyMatrix[i][j].M2 << '\n';
		std::cout << "lastContact = " << adjacencyMatrix[i][j].lastContact << '\n';
	}

};

class TwoHopContact {
private:
	Contacts** twoHopAdjacencyMatrix;
	int vertexCount;
public:
	TwoHopContact(int vertexCount) {
		this->vertexCount = vertexCount;
		twoHopAdjacencyMatrix = new Contacts*[vertexCount];
		for (int i = 0; i < vertexCount; i++) {
			twoHopAdjacencyMatrix[i] = new Contacts[vertexCount];
			for (int j = 0; j < vertexCount; j++){
				twoHopAdjacencyMatrix[i][j].edge = false;
				twoHopAdjacencyMatrix[i][j].n = 0;
				twoHopAdjacencyMatrix[i][j].mean = 0.0;
				twoHopAdjacencyMatrix[i][j].variance = 0.0;
				twoHopAdjacencyMatrix[i][j].M2 = 0.0;
				twoHopAdjacencyMatrix[i][j].lastContact = 0.0; 
			}
		}
	}

//  Contact(const Contact &m) {
//      this->vertexCount = m.vertexCount;
//      twoHopAdjacencyMatrix = new Contacts*[vertexCount];
//      for(int i = 0; i < vertexCount; i++) {
//          twoHopAdjacencyMatrix[i] = new Contacts[vertexCount];
//          for(int j = 0; j < vertexCount; j++){
//              twoHopAdjacencyMatrix[i][j].n =  m.twoHopAdjacencyMatrix[i][j].n;
//              twoHopAdjacencyMatrix[i][j].mean =  m.twoHopAdjacencyMatrix[i][j].mean;
//              twoHopAdjacencyMatrix[i][j].M2 =  m.twoHopAdjacencyMatrix[i][j].M2;
//              twoHopAdjacencyMatrix[i][j].lastContact =  m.twoHopAdjacencyMatrix[i][j].lastContact;
//          }
//      }
//  }

	void addEdge(int i, int j){
		twoHopAdjacencyMatrix[i][j].edge = true;
	}
	

	void contactBegin(int i, int j, double time){
		twoHopAdjacencyMatrix[i][j].edge = true; 
		twoHopAdjacencyMatrix[i][j].n += 1;
		twoHopAdjacencyMatrix[i][j].lastContact = time;
		//
		twoHopAdjacencyMatrix[j][i].edge = true;
		twoHopAdjacencyMatrix[j][i].n += 1;
		twoHopAdjacencyMatrix[j][i].lastContact = time;
	}


	void contactEnd(int i, int j, double time){
		twoHopAdjacencyMatrix[i][j].edge = false;
		double x = (time - twoHopAdjacencyMatrix[i][j].lastContact);
		double delta = x - twoHopAdjacencyMatrix[i][j].mean;
		twoHopAdjacencyMatrix[i][j].mean += delta/twoHopAdjacencyMatrix[i][j].n;
		twoHopAdjacencyMatrix[i][j].M2 += delta*(x - twoHopAdjacencyMatrix[i][j].mean);
		if (twoHopAdjacencyMatrix[i][j].n < 2) {
			twoHopAdjacencyMatrix[i][j].variance = 0.1;
		}
		else {
			twoHopAdjacencyMatrix[i][j].variance = twoHopAdjacencyMatrix[i][j].M2/(twoHopAdjacencyMatrix[i][j].n - 1);
		}
		//
		twoHopAdjacencyMatrix[j][i].edge = false;
		twoHopAdjacencyMatrix[j][i].mean = twoHopAdjacencyMatrix[i][j].mean;
		twoHopAdjacencyMatrix[j][i].M2 = twoHopAdjacencyMatrix[i][j].M2;
		twoHopAdjacencyMatrix[j][i].variance = twoHopAdjacencyMatrix[i][j].variance;
	}

	void contactUpdate(int i, int j, double time){
		double x = (time - twoHopAdjacencyMatrix[i][j].lastContact);
		double delta = x - twoHopAdjacencyMatrix[i][j].mean;
		twoHopAdjacencyMatrix[i][j].mean += delta/twoHopAdjacencyMatrix[i][j].n;
		//twoHopAdjacencyMatrix[i][j].M2 += delta*(x - twoHopAdjacencyMatrix[i][j].mean);
		//if (twoHopAdjacencyMatrix[i][j].n < 2) {
		//  twoHopAdjacencyMatrix[i][j].variance = 0;
		//}
		//else {
		//  twoHopAdjacencyMatrix[i][j].variance = twoHopAdjacencyMatrix[i][j].M2/(twoHopAdjacencyMatrix[i][j].n - 1);
		//}
		//
		twoHopAdjacencyMatrix[j][i].mean = twoHopAdjacencyMatrix[i][j].mean;
		twoHopAdjacencyMatrix[j][i].M2 = twoHopAdjacencyMatrix[i][j].M2;
		//twoHopAdjacencyMatrix[j][i].variance = twoHopAdjacencyMatrix[i][j].variance;

	}
	double degreeCent (int i) {
		double sum = 0;
		for(int j=0;j<vertexCount;j++)
		{
			if (twoHopAdjacencyMatrix[i][j].edge)
			{
				sum = sum + twoHopAdjacencyMatrix[i][j].mean;
			}
		}
		return sum;
	}

	bool isEdge (int i, int j){
		return twoHopAdjacencyMatrix[i][j].edge;
	}

	std::vector<int> dominatingSet(){
		std::vector< std::vector<int> > g; 
		int numEdge = 0;
		g.resize(vertexCount);
		for (int i = 0; i < vertexCount; i++) {
			for (int j = 0; j < vertexCount; j++) {
				if (isEdge(i, j)){
					numEdge++;
					g[i].push_back(j);
				}
			}
		}

		bool visit[vertexCount + 1];

		memset(visit, 0, sizeof(visit));

		std::vector<int> set;
		for (int i = 0; i < vertexCount; ++i)
		{
			if(!visit[i]) {
				set.push_back(i);
				visit[i] = true;
				for(int j = 0; j < (int)g[i].size(); j++) {
					if(!visit[g[i][j]]) {
						visit[g[i][j]] = true;
						break;
					}
				}
			}
		}
		return set;
	}
	double getlastContact (int i, int j){
		return twoHopAdjacencyMatrix[i][j].lastContact;
	}

	double getN (int i, int j){
		return twoHopAdjacencyMatrix[i][j].n;
	}

	double getMean (int i, int j){
		//double mean = twoHopAdjacencyMatrix[i][j].mean;
		//double x = (time - twoHopAdjacencyMatrix[i][j].lastContact);
		//double delta = x - mean; 
		//mean += delta/twoHopAdjacencyMatrix[i][j].n;
		return twoHopAdjacencyMatrix[i][j].mean; 
	}

	double getVariance (int i, int j){
		return twoHopAdjacencyMatrix[i][j].variance;
	}

	void printContacts(int i, int j){
		std::cout << "\nTwoHopContact - Contact Nodes ("<< i << ", " << j <<")" << std::endl;
		std::cout << "edge = " << twoHopAdjacencyMatrix[i][j].edge << '\n';
		std::cout << "n = " << twoHopAdjacencyMatrix[i][j].n << '\n';
		std::cout << "mean = " << twoHopAdjacencyMatrix[i][j].mean << '\n';
		std::cout << "M2 = " << twoHopAdjacencyMatrix[i][j].M2 << '\n';
		std::cout << "lastContact = " << twoHopAdjacencyMatrix[i][j].lastContact << '\n';
	}

};
// --- EDITING ----//
class Cluster {
private:
   
	int numClusters;
	int nNodes;
	int* clusterSize;
	double clusterpercent;
	double treshold;
	static bool instanceFlag;
	static Cluster *clusterUtil;
	int chNumber;
	std::vector<bool> testeClusterChanged;
	
	//std::map<int, Graph*> clustersGraph;

public:
	bool ch_change;
	bool reintegrateCh;
	std::map<int, std::vector<int>> clusterMembers; // Para cada clusterID, temos um conjunto de membros
	std::map<int, std::vector<int>> chPerRound;
	std::vector<int> clusterHeads; // Para cada clusterID, temos um CH
	std::vector<int> CHtoReintegrate;
	std::map <int, bool> nodeChangeCH;
	std::vector<int> movingNodes;
	int rnd;
	int cha;

	
	
	Cluster(int nNodes) {
		this->nNodes = nNodes;
		this->ch_change = false;
		this->reintegrateCh = false;
		this->clusterpercent = 0.1;
		this->cha = 0;
		for (int i = 0; i < nNodes; i++) {
			this->testeClusterChanged.push_back(false);
		}
		rnd = 0;

		this->treshold = this->clusterpercent/(1-this->clusterpercent*(rnd % 10));

		for (int i = 0; i < nNodes; i++) {
			nodeChangeCH.insert(std::pair<int, bool>(i, false));
		}
	}
	
	//Definido fora das chaves
	static Cluster* getInstance();
	static void resetInstance();
	void resetReintegrate() {
		if (reintegrateCh){
			reintegrateCh = false;
		} else {
			reintegrateCh = true;
		}
	}
	bool clusterChange(int node){
		if (node != -1)
			return testeClusterChanged[node];
		return false;
	}

	void toggleClusterChange(bool all, int node){

		if (all){
			for (int i = 0; i < nNodes; i++){
				if(testeClusterChanged[i]){
					testeClusterChanged[i] = false;
				} else {
					testeClusterChanged[i] = true;
				}
			}
		} else {
			if (node != -1){
				if(testeClusterChanged[node]){
					testeClusterChanged[node] = false;
				} else {
					testeClusterChanged[node] = true;
				}
			}
		}
	}

	void setChNumber(void) {
		std::map<int, std::vector<int>>::iterator it;
		int cont = 0;
		for(it = clusterMembers.begin(); it != clusterMembers.end(); it++) {
			cont++;
		}
		this->chNumber = cont;
	}
	
	int getChNumber(void) {
		return this->chNumber;
	}

	bool nodeChangeMap_Change() {
		for (std::map <int, bool>::iterator nodeChangeCHIt = nodeChangeCH.begin(); nodeChangeCHIt != nodeChangeCH.end(); nodeChangeCHIt++) {
			if (nodeChangeCHIt->second == true) {
				return true;
			}
		}
		return false;
	}
	void debugNodeCHangeMap() {
		for (std::map <int, bool>::iterator nodeChangeCHIt = nodeChangeCH.begin(); nodeChangeCHIt != nodeChangeCH.end(); nodeChangeCHIt++) {
			std::cout << "Node: " << nodeChangeCHIt->first << " status: " << nodeChangeCHIt->second << std::endl;
		}
	}
	void resetNodeChangeCHmap(bool val) {
		for (std::map <int, bool>::iterator nodeChangeCHIt = nodeChangeCH.begin(); nodeChangeCHIt != nodeChangeCH.end(); nodeChangeCHIt++) {
			nodeChangeCHIt->second = val;
		}
	}
	void NodeChangeCH(int node) {
		std::map<int,bool >::iterator nodeChangeCHIt = nodeChangeCH.find(node);
		if (nodeChangeCHIt != nodeChangeCH.end()) {
			nodeChangeCHIt->second = true;
		}
	}
	void NodeNotChangeCH(int node) {
		std::map<int,bool >::iterator nodeChangeCHIt = nodeChangeCH.find(node);
		if (nodeChangeCHIt != nodeChangeCH.end()) {
			nodeChangeCHIt->second = false;
		}
	}
	bool isNodeChangeCH(int node) {
		std::map<int,bool >::iterator nodeChangeCHIt = nodeChangeCH.find(node);
		if (nodeChangeCHIt != nodeChangeCH.end()) {
			if(nodeChangeCHIt->second) {
				return true;
			}
		}
		return false;
	}

	void addMovingNode(int node){
		this->movingNodes.push_back(node);
	}

	void resetMovingNodes() {
		this->movingNodes.clear();
	}
	bool inMovingNodes(int node){
		for (std::vector<int>::iterator it = this->movingNodes.begin(); it != this->movingNodes.end(); it++){
			if(node == *it) 
				return true;
		}
		return false;
	}

	void setTreshold() {
		this->treshold = this->clusterpercent/(1-this->clusterpercent*(rnd % 10));
	}

	double getTreshold() {
		return this->treshold;
	}

	void updateParams() {
		this->rnd++;
		//std::cout << "Round = " << this->rnd << std::endl;
		setTreshold();
	}

	void addClusterMember(int id, int nodeNum){
		std::map<int, std::vector<int> >::iterator it;
		it = clusterMembers.find(id);
		if(it != clusterMembers.end()) {
			//std:: cout << "Node " << nodeNum << " adicionado ao cluster " << id << std::endl;
			it->second.push_back(nodeNum); 

		} 
	}
	// Verifica se node faz parte de algum cluster e retorna como parâmetro o cluster ao qual ele faz parte
	bool isMember(int *id, int node) {
		std::map<int, std::vector<int>>::iterator it;
		for (it = clusterMembers.begin(); it != clusterMembers.end(); it++) {
			for (unsigned int i = 0; i < it->second.size(); i++) {
				if(it->second.at(i) == node) {
					*id = it->first;
					return true;
				}
			}
		}
		return false;
	}
	// Verifica se node faz parte de algum cluster
	bool isMember(int id,  int node) {
		std::map<int, std::vector<int>>::iterator it;
		it = clusterMembers.find(id);
		if(it != clusterMembers.end()) {
			for (unsigned int i = 0; i < it->second.size(); i++) {
				if (it->second.at(i) == node) {
					return true;
				}
			}
		} 
		return false;
	}

	std::vector<int> getCluster(int* id, int node){
		std::vector<int> members;
		std::map<int, std::vector<int>>::iterator it;

		for (it = clusterMembers.begin(); it != clusterMembers.end(); it++) {
			for (unsigned int i = 0; i < it->second.size(); i++) {
				if(it->second.at(i) == node) {
					// retorna o CH por referência
					*id = it->first;
					getMembers(it->first, &members);
					return members;
				}
			}
		}
		return members;
	}

	void getMembers(int id, std::vector<int> *v) {
		std::map<int, std::vector<int>>::iterator it;
		it = clusterMembers.find(id);
		if(it != clusterMembers.end()) {
			*v = it->second;
		}
		//std::cout << "Cluster ID não identificado" << std::endl;
		
	}

	bool removeClusterMember(int id, int nodeNum) {
		std::map<int, std::vector<int>>::iterator it;
		it = clusterMembers.find(id);
		if(it != clusterMembers.end()) {
			std::vector<int>::iterator vecIt;      
			for (vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++){
				if (*vecIt == nodeNum){
					it->second.erase(vecIt);
					//std::cout << "cluster - (removeClusterMember) Node " << *vecIt << " retirado do cluster " << id << " " << std::endl;
					return true;
					//std::cout << "node "<< nodeNum << " erase " << std::endl; 
				}
			}
			return false;
		}
		std::cout << "removeClusterMember - cluster " << id << " inexistente" << std::endl;
		return false;
	}

	void printClusters(){
		std::map<int, std::vector<int>>::iterator it;
		for(it = clusterMembers.begin(); it != clusterMembers.end(); it++) {
			std::vector<int>::iterator vecIt;  
			std::cout << "Cluster " << it->first << ": ";
			for (vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++){
				std::cout << *vecIt << " ";
			}
			std::cout << std::endl;
		}
	}

	int setClusterHead(int nodeNum) {
		std::map<int, std::vector<int> >::iterator it;
		it = clusterMembers.find(nodeNum);
		
		std::vector<int> tmp;
		int id;
		if (isMember(&id, nodeNum))
			removeClusterMember(id, nodeNum);
		tmp.push_back(nodeNum);
		std:: cout << "Node " << nodeNum << " adicionado ao cluster " << nodeNum << std::endl;
		clusterMembers.insert(std::pair<int, std::vector<int>>(nodeNum, tmp));
		return CH_CHANGE;
		
	}

	bool removeEmptyClusters() {
		std::map<int, std::vector<int>>::iterator it;
		bool hasEmpty = false;
		
		for(it = clusterMembers.begin(); it != clusterMembers.end(); it++) {

			//std::cout << "Cluster " << it->first << " tamanho " << it->second.size() << std::endl;
			if(emptyCluster(it->first, &CHtoReintegrate)) {
				hasEmpty = true;
				
				// if (std::find(CHtoReintegrate.begin(), CHtoReintegrate.end(), it->first) == CHtoReintegrate.end()) {
				//     CHtoReintegrate.push_back(it->first);
				// }
				//std::cout << "Cluster " << it->first << " vazio" << std::endl;
				clusterMembers.erase(it->first);
			}
			//std::cout << " não vazio " ;
		}
		return hasEmpty;
	}

	bool emptyCluster(int id, std::vector<int> *CHtoReintegrate) {
		std::map<int, std::vector<int>>::iterator it;
		it = clusterMembers.find(id);
		if(it != clusterMembers.end()) {
			if (it->second.size() == 2){
				std::vector<int>::iterator vecIt;
				for(vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++){
					if (std::find(CHtoReintegrate->begin(), CHtoReintegrate->end(), *vecIt) == CHtoReintegrate->end()) {
						CHtoReintegrate->push_back(*vecIt);
						clusterMembers.erase(it->first);
					}
				}
				return true;
			}

			if (it->second.size() == 1)
			{
				if(it->second.at(0) == id){
					CHtoReintegrate->push_back(it->second.at(0));
					clusterMembers.erase(it->first);
					return  true;
				}
				else
					return false;
			}

			if(it->second.size() == 0){
				clusterMembers.erase(it->first);
				return true;
			}
			else
				return false;
		}
		//std::cout << "emptyCluster - cluster " << id << " inexistente" << std::endl;   
		return false;
	}

	bool isClusterHead(int id) {
		std::map<int, std::vector<int>>::iterator it;
		it = clusterMembers.find(id);
		if(it != clusterMembers.end()) {
			return true;
		}
		return false;
	}

	void removeClusterHead(int nodeNum) {
		if (!isClusterHead(nodeNum)){
			std::vector<int>::iterator vecIt;
			for(vecIt = clusterHeads.begin(); vecIt != clusterHeads.end(); vecIt++)
			{
				if (*vecIt == nodeNum){
					clusterHeads.erase(vecIt);
					std::cout << "node "<< nodeNum << " erase " << std::endl; 
				}
			}
		}
	}

	void updateHistory() {
		std::map<int, std::vector<int>>::iterator it;
		it = chPerRound.find(this->rnd);
		if(it == chPerRound.end()) {
			std::vector<int> vec;
			std::map<int, std::vector<int>>::iterator itcluster;
			for(itcluster = clusterMembers.begin(); itcluster != clusterMembers.end(); itcluster++){
				vec.push_back(itcluster->first);
			}
			chPerRound.insert(std::pair<int, std::vector<int>>(this->rnd, vec));
		}

		/*for(it = chPerRound.begin(); it != chPerRound.end(); it++) {
			std::vector<int>::iterator vecIt;  
			std::cout << "Round " << it->first << ": ";
			for (vecIt = it->second.begin(); vecIt != it->second.end(); vecIt++){
				std::cout << *vecIt << " ";
			}
			std::cout << std::endl;
		}*/
	}
	bool getLastRoundHistory(int node) {
		std::map<int, std::vector<int>>::iterator it;
		it = chPerRound.find(this->rnd);
		if(it != chPerRound.end()) {
			if (std::find(it->second.begin(), it->second.end(), node) != it->second.end()) {
				//std::cout << "CH " << node << " encontrado na rodada passada" << std::endl;
				return true;
			}
			return false;
		}
		return false;
	}    
};

bool Cluster::instanceFlag = false;
Cluster* Cluster::clusterUtil = NULL;

Cluster* Cluster::getInstance(){
	if (!instanceFlag) {
		BasicConfig* basicConfigInstance = BasicConfig::getInstance();
		uint16_t number_nodes = basicConfigInstance->getNumberOfNodes();

		clusterUtil = new Cluster(number_nodes);
		instanceFlag = true;
	}
	return clusterUtil;
}

void Cluster::resetInstance(){
	instanceFlag = false;
}

class HandleQueue {
private:
	std::queue<int> m_queue;
	std::set<int> aux_set;
	int cont;
public:
	HandleQueue(){cont = 0;}

	void add (int c){
		m_queue.push(c);
		aux_set.insert(c);
		cont ++;
	}
	int getNext(){
		if (!isEmpty())
			return m_queue.front();
		else
			return -1;
	}

	int pop() {
		int r = m_queue.front();
		aux_set.erase(r);
		m_queue.pop();
		cont--;
		return r;
	}

	std::set<int> getQueueElements(){
		return aux_set;
	}

	bool isEmpty(){
		return (cont == 0);
	}

	int count(){
		return (int) cont;
	}
	
	bool hasMember(int c){
		for (auto const &p : aux_set){
			if (p == c) return true;
		}
		return false;
	}

	void print(){
		for (auto const &p : aux_set){
			std::cout << p << " ";
		}
		std::cout << std::endl;
	}
};

class Content {
private:
	bool** contentMatrix;
	int height;
	int width;
	unsigned int cacheSize;
	int** frequency;
	int* views;
	
	std::vector<int>* cache;
	
	std::vector<int>* expectedCache;


	static bool instanceFlag;
	static Content *contUtil;
public:
	std::vector<HandleQueue> tempCache; // Quem enviou, conteúdo
	std::vector<HandleQueue> tempCacheSenders;
	int nc ;
	std::string metadata; // Temporario, enquanto não extraio os metadados
	bool verifyMetadata;
	bool random;
	int nextCache;
	int cacheHit;
	int cacheMiss;
	int localCache;
	int cacheEvict;
	int d2dSend;
	int proCache;
	int proCached2d;
	int d2dCache;
	int friendCache;

	bool lc;
	unsigned int expectedT;
	unsigned int expectedA;

	int lastContAdded;
	int lastNodeCont;

	static Content* getInstance(int height, int width, unsigned int cacheSize);
	static void  resetInstance();
	static Content* getInstance();

	Content(int height, int width, unsigned int cacheSize) {
		this->height = height;
		this->width = width;
		this->cacheSize = cacheSize;
		this->metadata = std::string("metadata");
		this->nc = -1;
		this->verifyMetadata = false;
		this->random = false;
		this->nextCache = -1;
		this->d2dCache = -1;
		this->d2dSend = -1;
		this->cacheHit = 0;
		this->cacheMiss = 0;
		this->localCache = 0;
		this->cacheEvict = 0;
		this->proCache = 0;
		this->proCached2d = 0;
		this->friendCache = 0;
		this->expectedT = 0;
		this->expectedA = 0;
		this->lastContAdded = -1;
		this->lastNodeCont = -1;
		this->lc = false;

		cache = new std::vector<int>[height];
		for (int i = 0; i < height; i++) {
			cache[i].reserve(cacheSize);
		}

		for (int i = 0; i < height; i++){
			tempCache.push_back(HandleQueue());
		}

		for (int i = 0; i < height; i++){
			tempCacheSenders.push_back(HandleQueue());
		}
		
		views = new int[width];
		for (int i = 0; i < width; i++){
			views[i] = 0;
		}
		frequency = new int*[height];
		for (int i = 0; i < height; i++) {
			frequency[i] = new int[width];
		}
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
				frequency[i][j] = 0;
		}

		contentMatrix = new bool*[height];
		for (int i = 0; i < height; i++) {
			contentMatrix[i] = new bool[width];
		}
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
				contentMatrix[i][j] = false;
		}

		expectedCache = new std::vector<int>[height];
		for (int i = 0; i < height; i++) {
			expectedCache[i].reserve(height);
		}
	}

	void addContent(int n, int c) {
		if (n >= 0 && n < height && c > 0 && c < width) {
			contentMatrix[n][c] = true;
		}
	}

	void addCache(int n, int c) {
		if (n >= 0 && n < height && c > 0 && c < width) {
			//if (cache[n].empty()){
			//  cache[n].push_back(c);
			//} else
			if (std::find(cache[n].begin(), cache[n].end(), c) != cache[n].end()){
				std::cout << "cont - Content " << c << " already cached on node " << n << std::endl;
			}
			else if (cache[n].size() < cacheSize){
				std::cout << "cont - Added content " << c << " to node " << n << std::endl;
				cache[n].push_back(c);

				//tempCache[n].push_back(c);
				std::cout << "cont - Cache Size of node "<<n<<": " << cache[n].size() << std::endl;
			}
			else if (cache[n].size() >= cacheSize){
				//std::cout << "cont - CacheSize  " << cache[n].size() << " maxSize " << cache[n].max_size() << std::endl;
				int min_index = 0;

				for(int i = 1; i < width; i++){
					if (frequency[n][i] < frequency[n][min_index] && 
						(std::find(cache[n].begin(), cache[n].end(), i)!= cache[n].end())){
						min_index = i; 
					}
				}
				std::cout << "cont - Replaced content " << cache[n].at(min_index) << " for "<< c << " on node "<< n << std::endl;
				//std::vector<int>::iterator newEnd = std::remove(cache[n].begin(), cache[n].end(), min_index);
				//cache[n].erase(newEnd, cache[n].end());
				cache[n].erase(cache[n].begin()+min_index);
				frequency[n][min_index] = 0;
				cache[n].push_back(c);
				//tempCache[n].push_back(c);      
				cacheEvict++; 
			}
			cacheHit++;

			for (std::vector<int>::const_iterator i = cache[n].begin(); i != cache[n].end(); ++i)
				std::cout << *i << ' ';
			std::cout << std::endl;  
		}
	}

	void addCache(int n, int sender, int c, bool d2d) {
		if (n >= 0 && n < height && c > 0 && c < width) {
			//if (cache[n].empty()){
			//  cache[n].push_back(c);
			//} else
			if (std::find(cache[n].begin(), cache[n].end(), c) != cache[n].end()){
				std::cout << "cont - Content " << c << " already cached on node " << n << std::endl;
			}
			else if (tempCache[n].hasMember(c)) {
				std::cout << "cont - Content " << c << " already on temp cache node " << n << std::endl;
			}
			else if (cache[n].size() < cacheSize){
				std::cout << "cont - Added content "<< c << " to temporary queue of node " << n << std::endl;
				std::cout << "Cont sended by " << sender << std::endl;

				tempCache[n].add(c);
				tempCacheSenders[n].add(sender);
				if (d2d) {
					lastContAdded = c;
					lastNodeCont = n;
					std::cout << "cont - Added content " << c << " from d2d connection " << std::endl;
				}
				//tempCache[n].push_back(c);
				std::cout << "cont - Temp Cache Size of node "<< n <<": " << tempCache[n].count() << std::endl;
			}
			std::cout << "cont - Temp Cache: ";
			tempCache[n].print();
			tempCacheSenders[n].print();
			std::cout << std::endl;
		}
	}

	void removeContent(int i, int j) {
		if (i >= 0 && i < height && j > 0 && j < width) {
			contentMatrix[i][j] = false;
		}
	}

   /* void removeTempCache(int n, int c) {
		if (n >= 0 && n < height && c > 0 && c < width) {
			std::cout << "cont - Temp Cache: ";
			tempCache[n].print();
			
			
			if (tempCache[n].hasMember(c)){
				tempCache[n].remove(c);
				std::cout << "cont - Erased tempCache of node " << n << std::endl;
			}
			std::cout << "Content not in temp cache" << std::endl;
		}
	}*/

	bool hasCache(int n, int c) {
		if (n >= 0 && n < height && c > 0 && c < width){
			if (std::find(cache[n].begin(), cache[n].end(), c) != cache[n].end()){
				return true;
			}
		}
		return false;
	}

	bool oldContent(int j) {
		//if (j >=0 && j < width){
		for (int i = 0; i < height; i++){
			if (contentMatrix[i][j])
				return true;
		}   
		return false;
	}

	unsigned int sumLine (unsigned int x) {
		unsigned int sum = 0;
		//std::cout<<"sumline_metodo X = "<<x<< " "; 
		//std::cout<<std::endl;
		for(int y=0;y<width;y++)
		{
			if (contentMatrix[x][y])
			{
				//std::cout<<"sumline_laco X = "<<x<< " ";
				//std::cout<<std::endl;
				sum = sum + 1;
			}
		}
		return sum;
	}

	int zeros () {
		int sum = 0;
		//std::cout<<"sumline_metodo X = "<<x<< " ";
		//std::cout<<std::endl;
		for(int y=0;y<width;y++)
		{
			if (views[y]==0)
				sum = sum + 1;
		}
		return sum;
	}

	void addViews (int x){
		views[x] += 1;
	}

	int getViews (unsigned int x) {
		return views[x];
	}


	void contentReplace (int x, int v) {
		int min = std::numeric_limits<int>::max();;
		int minPos = -1;
		//int max = 0;
		for(int y=0;y<width;y++)
		{
		//if (contentMatrix[x][y])
			if (frequency[x][y] < min || frequency[x][y] > -1 )
			{
				min = frequency[x][y];
				minPos = y;
			}
		}
		if (minPos != v || minPos > -1){
		//frequency[x][v]=0;
			contentMatrix[x][minPos] = false;
			contentMatrix[x][v] = true;
			//this->removeContent(x,minPos);
			//this->addContent(x,v);
			std::cout << "Replace " << minPos << " for " << v << " on node " << x <<std::endl;
		}
	}

	void addFreq (int n, int c) {
		frequency[n][c] += 1;
		std::cout << "Content " << c << " freq. " << frequency[n][c] << " on node " << n <<std::endl;
	}

	void addExpected (int c, int e) {
		expectedCache[c].push_back(e);
		expectedT++;    
	} 

	void isExpected (int c, int e) {
		if (std::find(expectedCache[c].begin(), expectedCache[c].end(), e) != expectedCache[c].end()){
			expectedA++;
		}
	}

	void printViews () {
		for(int y=0;y<width;y++) 
		// loop for the three elements on the line
		{
			std::cout << views[y] << " "; 
			// display the current element out of the array
		}
		std::cout<<std::endl; 
		// when the inner loop is done, go to a new line
	}

	//void printEdge () {
	//  for(int x=0;x<vertexCount;x++) 
	//  //loop 3 times for three lines
	//  {
	//      for(int y=0;y<vertexCount;y++) 
	//      //loop for the three elements on the line
	//      {
	//          std::cout<<adjacencyMatrix[x][y]; 
	//          //display the current element out of the array
	//      }
	//      std::cout<<std::endl; 
	//      //when the inner loop is done, go to a new line
	//  }
	//}

	double getCHrate () {
		double chRate = (cacheHit + localCache)/(cacheHit + localCache + cacheMiss);
		return chRate;
	}

	double getOFFrate () {
		double offloadRate = (cacheHit + localCache)/(cacheHit + localCache + cacheMiss + proCache);
		return offloadRate;
	}

	~Content() {
		for (int i = 0; i < height; i++){
			delete[] contentMatrix[i];
			delete[] frequency[i];        }
		delete[] contentMatrix;
		//for (int i = 0; i < height; i++)
		//  delete[] frequency[i];
		delete[] frequency;
		delete[] views;
	}
};

bool Content::instanceFlag = false;
Content* Content::contUtil = NULL;

Content* Content::getInstance(int height, int width, unsigned int cacheSize){
	if (!instanceFlag) {
		contUtil = new Content(height, width, cacheSize);
		instanceFlag = true;
	}
	return contUtil;
}

void Content::resetInstance(){
	instanceFlag = true;
}

Content* Content::getInstance(){
	if (!instanceFlag) {
		std::cout << "To call Content::getInstance for the first time /params 'int height, int width, unsigned int cacheSize' needed " << std::endl;
	}
	return contUtil;
}


class View {
private:
	double** viewMatrix;
	int height;
	int width;
public:
	int nextContent;
	int nextNode;

	View(int height, int width) {
		this->height = height;
		this->width = width;
		this->nextContent = width+1;
		this->nextNode = width+1;
		viewMatrix = new double*[height];
		for (int i = 0; i < height; i++) {
			viewMatrix[i] = new double[width];
		}     
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
				viewMatrix[i][j] = INFINITY;
		}           
	} 

	void addView(int i,  int j, double v) { 
		if (i >= 0 && i < height && j > 0 && j < width) {
			viewMatrix[i][j] = v;
		}
	}

	void updateView(int i,  int j, double v) {
		if (i >= 0 && i < height && j > 0 && j < width && v < viewMatrix[i][j]) {
			viewMatrix[i][j] = v;
		}
	}     

	void removeView(int i, int j) {
		if (i >= 0 && i < height && j > 0 && j < width) {
			viewMatrix[i][j] = false;
		}
	}

	double getView(int i, int j) {
		//if (i >= 0 && i < height && j > 0 && j < width) 
		return viewMatrix[i][j];
	}

	void nextView(int i, int j) {
		double lastView = viewMatrix[i][j];
		double nextView = INFINITY;
		for(int i=0;i<height;i++) {
			for(int j=0;j<width;j++) {
				if(viewMatrix[i][j] > lastView && viewMatrix[i][j] < nextView ){
					nextView=viewMatrix[i][j];
					//std::cout << this->nextNode << " " << this->nextContent << std::endl;
					this->nextNode = i;
					this->nextContent = j;
					//std::cout << this->nextNode << " " << this->nextContent << std::endl;
				}
			}
		}
	}

	double getNext(int i, int j) {
		double lastView = viewMatrix[i][j];
		double nextView = INFINITY;
		for(int i=0;i<height;i++) {
			for(int j=0;j<width;j++) {
				if(viewMatrix[i][j] > lastView && viewMatrix[i][j] < nextView ){
					nextView=viewMatrix[i][j];
				}
			}
		}
		return nextView;
	}

	void firstView() {
		double nextView = INFINITY;
		for(int i=0;i<height;i++) {
			//std::cout << i << std::endl;
			for(int j=0;j<width;j++) {
				//std::cout << j << std::endl;
				if(viewMatrix[i][j] < nextView ){
					nextView=viewMatrix[i][j];
					this->nextNode = i;
					this->nextContent = j;
					//std::cout << "first " << this->nextNode << " " << this->nextContent << std::endl;
				}
			}             
		}     
	}

	void printEdge () {
		//std::cout << height << " " << width <<std::endl;
		for(int x=0;x<height;x++) 
		// loop 3 times for three lines
		{
			for(int y=0;y<width;y++) 
			// loop for the three elements on the line
			{
				std::cout<<viewMatrix[x][y] << " "; 
				// display the current element out of the array
			}
			std::cout<<std::endl; 
			// when the inner loop is done, go to a new line
		}
		//return 0;
	}
	
	~View() {
		for (int i = 0; i < height; i++)
			delete[] viewMatrix[i];
		delete[] viewMatrix;
	}
};


//class Energy {
//
//private:
//  nodeEnergy* energyArray;
//  int vertexCount;
//  double duration;
//public:
//  Contact(int vertexCount, double duration) {
//      this->vertexCount = vertexCount;
//      this->duration = duration;
//      energyArray = new nodeEnergy*[vertexCount];
//      for (int i = 0; i < vertexCount; i++) {
//          energyArray[i].upE = 0.0;
//          energyArray[i].downE = 0.0;
//          energyArray[i].transitionE = 0.0;
//          energyArray[i].activeTime = 0.0;
//      }
//  }
//
//  void updateE (int n, double upEn, double downEn, double transition, double flowDuration){
//      energyArray[n].upE += upEn;
//      energyArray[n].downE += downEn;
//      energyArray[n].transitionE += transition;
//      energyArray[n].activeTime += flowDuration;
//  }
//
//  void printTotalEnergy (){
//      for(int i=0;i<vertexCount;i++){
//          
//      } 
//  }
//};


class BlockchainConnection {
private:
	
	std::map<int, std::set<std::pair<int, int>>> m_blockchain;

	static bool instanceFlag;
	static BlockchainConnection *bUtil;
	int start;
public:

	static BlockchainConnection* getInstance(int nNodes);
	static BlockchainConnection* getInstance();
	static void resetInstance();

	 BlockchainConnection(int nNodes){
		for(int i = 0; i < nNodes; i++){
			std::set<std::pair<int, int>> set;
			m_blockchain[i] = set;
		}
		start = 0;
	}
	void startInstance(){
		start = 1;
	}

	void insert(int node, int id, int decision){
	
		m_blockchain[node].insert(std::pair<int, int>(id, decision));
	}

	bool searchContent(int node, int id){
		for (auto &p: m_blockchain[node]) {
			if (p.first == id){
				return true;
			}
			
		}
		return false;
	}

	bool searchContent(int node, int id, int* decision){
		for (auto &p: m_blockchain[node]) {
			if (p.first == id){
				*decision = p.second; 
				return true;
			}
		}
		
		return false;
	}
	int getDecision(int node, int id){
		int decision;
		if (!searchContent(node, id, &decision)){
			return -1;
		}
		return decision;
	}

	bool checkContent(int node, int id){
		int decision = getDecision(node, id);
		if (decision == 1){
			std::cout << "(Node " << node <<  ") Content " << id << " is valid" << std::endl;
			return true;
		} else if (decision == 0){
			std::cout << "(Node " << node <<  ") Content " << id << " is invalid" << std::endl;
			return false;
		} else {
			//std::cout << "(Node " << node <<  ") Content " << id << " not in chain" << std::endl;
			return true;
		}
	}

	void reset(int node){
		m_blockchain[node].clear();
	}

};

bool BlockchainConnection::instanceFlag = false;
BlockchainConnection* BlockchainConnection::bUtil = NULL;

BlockchainConnection* BlockchainConnection::getInstance(int nNodes){
	if (!instanceFlag) {
		bUtil = new BlockchainConnection(nNodes);
		instanceFlag = true;
	}
	return bUtil;
}

void BlockchainConnection::resetInstance(){
	
	instanceFlag = false;
}

BlockchainConnection* BlockchainConnection::getInstance(){
	if (!instanceFlag) {
		std::cout << "To call BlockchainConnection::getInstance for the first time /params 'int height, int width, unsigned int cacheSize' needed " << std::endl;
	}
	return bUtil;
}


class ConsensusTimeStatics{
private:
	std::vector<std::pair<int, double>> t_vec;
	double totMed;
	std::map<int, double> t_medSizeTime;
	std::map<int, std::set<double>> t_sizeTime;
	static bool instanceFlag;
	static ConsensusTimeStatics *consensusTimeUtil;
	double m_controlPackets;
public:
   
	static ConsensusTimeStatics* getInstance();
	static void resetInstance();

	ConsensusTimeStatics(){
		this->totMed = 0.0;
		m_controlPackets = 0.0;
	}

	void add(int c_size, double time){
		std::pair<int, double> pair;
		pair.first = c_size;
		pair.second = time;
		t_vec.push_back(pair);
	}

	void addControlPacket(){
		m_controlPackets = m_controlPackets + 1.0;
	}

	double getControlPackets(){
		return m_controlPackets;
	}

	void calculate(){
		int tv_size = (int) t_vec.size();
		double v_tot = 0.0;
		for(int i = 0; i < tv_size; i++){
			v_tot += t_vec[i].second;
			t_sizeTime[t_vec[i].first].insert(t_vec[i].second);
		}

		this->totMed = v_tot/tv_size;
	}

	double getAverageTime(){
		if(this->totMed == 0.0){
			calculate();
		}
		//std::cout << "Consensus average conclusion time: " << this->totMed << std::endl;
		return this->totMed;
	}

	void getTimePerClusterSize(){
		std::cout << "Time per Cluster Size (clusterSize, time)\n";
		for (std::map<int, std::set<double>>::iterator it = t_sizeTime.begin(); it != t_sizeTime.end(); it++){
			int tot = 0.0;
			for (auto const &p: it->second){
				tot += p;
			}

			t_medSizeTime[it->first] = (tot / (it->second).size());

			std::cout << "\t-> (" << it->first << ", " << t_medSizeTime[it->first] << ")\n";
		}
		std::cout << std::endl;
	}
};

bool ConsensusTimeStatics::instanceFlag = false;
ConsensusTimeStatics* ConsensusTimeStatics::consensusTimeUtil = NULL;

ConsensusTimeStatics* ConsensusTimeStatics::getInstance(){
	if (!instanceFlag) {
		consensusTimeUtil = new ConsensusTimeStatics();
		instanceFlag = true;
	}
	return consensusTimeUtil;
}

void ConsensusTimeStatics::resetInstance(){
	instanceFlag = false;
}
