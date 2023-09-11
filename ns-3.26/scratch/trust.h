#include <vector>
#include <utility>
#include <fstream>
#include <stack>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <random>
#include <map>

#include "scratch/dempstershafer.h"
#include "scratch/bayesian.h"
// Classe de distribuição normal (gausiana) para escolha de nós em geral
class Generator {
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution;
	int min;
	int max;
public:
	Generator(int min, int max):
		distribution(min, max), min(min), max(max)
	{}

	int operator ()() {
		while (true) {
			int number = this->distribution(generator);
			if (number >= this->min && number <= this->max)
				return number;
		}
	}
};

class MaliciousNodes
{
private:
	double m_proportion;
	std::set<int> m_nodes;
	int m_qtd;
	int m_nNodes;
	static bool instanceFlag;
	static MaliciousNodes *maliciousUtil;


	void setMaliciousNodes(){
		std::cout << "Nós maliciosos: " << std::endl;
		Generator g(0, m_nNodes);
		for (int i = 0; i < m_nNodes * m_proportion; i++){
			int no = g();
			m_nodes.insert(no);
			std::cout << no << " ";
		}

		m_qtd = m_nodes.size();

		std::cout << std::endl;
	}

public:
	MaliciousNodes(double proportion, int nNodes) 
		: m_proportion(proportion), m_nNodes(nNodes){
		this->setMaliciousNodes();
	}	
	static MaliciousNodes* getInstance();
	static MaliciousNodes* getInstance(double proportion, int nNodes);
	int getQtd(){
		return m_qtd;
	}
	std::set<int> get(){
		return m_nodes;
	}
	
	bool isMalicious(int in_node){
		for (auto &node: m_nodes)
			if (node == in_node) return true;
		return false;
	}
	int selectANode(){
		for (auto &node: m_nodes){
			return node;	
		}
		return -1;		
	}
};

bool MaliciousNodes::instanceFlag = false;
MaliciousNodes* MaliciousNodes::maliciousUtil = NULL;

MaliciousNodes* MaliciousNodes::getInstance(double proportion, int numberOfNodes){
	if (!instanceFlag) {
		maliciousUtil = new MaliciousNodes(proportion, numberOfNodes);
		instanceFlag = true;
	}
	return maliciousUtil;
}

MaliciousNodes* MaliciousNodes::getInstance(){
	return maliciousUtil;
}


class TrustHistory{
private:
	std::map<int, std::vector<double>> nodesTrust;
	int m_nNodes;
public:
	TrustHistory(){
	}
	void setNumNodes(int nNodes){
		m_nNodes = nNodes;
	}

	TrustHistory(int nNodes): m_nNodes(nNodes){
		for (int i = 0; i < m_nNodes; i++) {
			std::vector<double> new_set;	
			nodesTrust[i] = new_set;
		}
	}

	void update(int node, double t_value){
		nodesTrust[node].push_back(t_value);
	}

	std::vector<double> get(int node){
		std::map<int, std::vector<double>>::iterator it_trusts = nodesTrust.find(node);
		if (it_trusts != nodesTrust.end()) {
			return nodesTrust[node];
		}
	}
	void print(int node){
		std::map<int, std::vector<double>>::iterator it_trusts = nodesTrust.find(node);
		if (it_trusts != nodesTrust.end()) {
			std::cout << "Histórico de confiança do nó " << node << ": ";
			for (auto const &p: nodesTrust[node]){
				std::cout << p << " ";
			}
			std::cout << std::endl;
		}
	}
};

class NodeTrust{
private:
	BayesianInference* directTrust_values;
	DempsterShafer* indirectTurst_values;
	TrustHistory dHistory;
	TrustHistory iHistory;
	TrustHistory history;

	int m_node;
	int m_nNodes;
	//double* directTrust_values; // Valores de confiança

	double* trust_values;	
	double omega;
	
	void init(){
		dHistory.setNumNodes(m_nNodes);
		iHistory.setNumNodes(m_nNodes);
		history.setNumNodes(m_nNodes);
		
		directTrust_values = new BayesianInference[m_nNodes];
		indirectTurst_values = new DempsterShafer[m_nNodes];
		trust_values = new double[m_nNodes];
		omega = 0.4;
		for (int i = 0; i < m_nNodes; i++){
			trust_values[i] = 0.5;
			dHistory.update(i, directTrust_values[i].getTrust());
			iHistory.update(i, indirectTurst_values[i].getTrust());
			history.update(i, trust_values[i]);
		}
	}
public:
	NodeTrust(int node_id, int nNodes) : m_nNodes(nNodes){
		m_node = node_id;
		init();
	}
	double get(int node){
		if (node < m_nNodes)
			return trust_values[node];
		return -1;
	}

	double getDirect(int node){
		if (node < m_nNodes)
			return directTrust_values[node].getTrust();
		return -1;
	}

	double getIndirect(int node){
		if (node < m_nNodes)
			return indirectTurst_values[node].getTrust();
		return -1;
	}
	double addGoodBehavior(int node, double val){
		directTrust_values[node].addGoodBehavior(val);
		dHistory.update(node, directTrust_values[node].getTrust());
		updateFinalTrust(node);
		return this->get(node);
	}
	double addBadBehavior(int node, double val){
		directTrust_values[node].addBadBehavior(val);
		dHistory.update(node, directTrust_values[node].getTrust());
		updateFinalTrust(node);
		return this->get(node);
	}

	double addEvidence(int node, int opinion_node, bool trustable){
		if (trustable){
			indirectTurst_values[node].addTrustEvidence(trust_values[opinion_node]);
		} else {
			indirectTurst_values[node].addUntrustEvidence(trust_values[opinion_node]);
		}
		iHistory.update(node, indirectTurst_values[opinion_node].getTrust());
		updateFinalTrust(node);
		return this->get(node);
	}

	void print(int node){
		if (node < m_nNodes){
			std::cout << "Node " << m_node + 1 << ": ";
			for (int i = 0; i < m_nNodes; i++){
				std::cout << "(" << i + 1 << ") " << get(i) << ", ";
			}
			std::cout << std::endl;
		}
	}
	double updateFinalTrust(int node){
		trust_values[node] = omega*directTrust_values[node].getTrust() + (1-omega)*indirectTurst_values[node].getTrust();
		return trust_values[node];
	}

	std::vector<double> getHistory(int node){
		return  history.get(node);
	}
};

class Trust{
private:
	std::vector<NodeTrust*> vectorOfTrust;
	int** contact;
	int m_nNodes;
	bool m_is_complete;
	static bool instanceFlag;
	static Trust *trustUtil;
	double limiar;
public:
	Trust(int nNodes): m_nNodes(nNodes){
		limiar = 0.5;
		for (int i = 0; i < m_nNodes; i++){
			vectorOfTrust.push_back(new NodeTrust(i, m_nNodes));
		}
		contact = new int*[m_nNodes];

		for (int i = 0; i < m_nNodes; i++){
			contact[i] =  new int[m_nNodes];
			for (int j = 0; j < m_nNodes; j++) {
				contact[i][j] = 0;
			}
		}
	}
	Trust(){}
	static Trust* getInstance(int numberNOfodes);
	static Trust* getInstance();
	void setNumberOfNodes(int nNodes){
		m_nNodes = nNodes;
		for (int i = 0; i < m_nNodes; i++){
			vectorOfTrust.push_back(new NodeTrust(i, m_nNodes));
		}
	}

	void setIsComplete(bool complete){
		m_is_complete = complete;
	}
	
	double goodBehavior(int node, int get_node, double val){
		contact[node][get_node] = 1;
		return vectorOfTrust[node]->addGoodBehavior(get_node, val);
	}
	double badBehavior(int node, int get_node, double val){
		contact[node][get_node] = 1;
		return vectorOfTrust[node]->addBadBehavior(get_node, val);
	}

	double addEvidence(int node, int get_node, int opinion_node, bool trustable){
		contact[node][get_node] = 1;
		return vectorOfTrust[node]->addEvidence(get_node, opinion_node, trustable);
	}

	double getTrust(int node, int get_node){
		if (m_is_complete) {
			return vectorOfTrust[node]->get(get_node);
		} 
		return vectorOfTrust[node]->getDirect(get_node);
		
	}

	bool hasPreviousContact(int node, int get_node){
		return (contact[node][get_node] == 1);
	}

	bool is_trustable(int node, int get_node){
		if (getTrust(node, get_node) >= 0.5){
			return true;
		}
		return false;
	}

	bool is_trustable_in_consensus(int node, int get_node){
		if (getTrust(node, get_node) > 0.5){
			return true;
		}
		return false;
	}

	bool even_trustable_consensus(int node, int get_node){
		if (getTrust(node, get_node) == 0.5){
			return true;
		}
		return false;
	}

	double avgTrust(int node){
		double sum = 0.0;
		for (int i = 0; i < m_nNodes; i++){
			sum += vectorOfTrust[i]->get(node);
		}
		return sum/(m_nNodes);
	}

	double avgContactTrust(int node){
		double sum = 0.0;
		int cont = 0;
		for (int i = 0; i < m_nNodes; i++){
			if(contact[i][node] == 1){
				sum += vectorOfTrust[i]->get(node);
				cont++;
			}
		}
		if (cont > 0)
			return sum/(cont);

		return (double) 0.5;
	}

	double getDirectTrust(int node, int get_node){
		return vectorOfTrust[node]->getDirect(get_node);
	}
	double getIndirectTrust(int node, int get_node){
		return vectorOfTrust[node]->getIndirect(get_node);
	}

	std::vector<double> getHistory(int node, int get_node){
		return vectorOfTrust[node]->getHistory(get_node);
	}
	void printValues(){
		for (int i = 0; i < m_nNodes; i++){
			std::cout << "Node " << i+1 << ": ";
			for (int j = 0; j < m_nNodes; j++)
				std::cout << "(" << j + 1 << ") " << getTrust(i, j) << ", ";

			std::cout << std::endl;
		}
	}
};

bool Trust::instanceFlag = false;
Trust* Trust::trustUtil = NULL;

Trust* Trust::getInstance(int numberOfNodes){
	if (!instanceFlag) {
		trustUtil = new Trust(numberOfNodes);
		instanceFlag = true;
	}
	return trustUtil;
}

Trust* Trust::getInstance(){
	return trustUtil;
}

struct AverageTrust{
	double timestamp;
	double avg_trust;
	double avg_contact_trust;
};

struct BadContent{
	double timestamp;
	int bad_content;
	int total_bad_content;
	double badContentDeliveryRate;
	double badContentDeliverySuccessRate;
};

class TrustStatistics {
private:
	int m_nNodes;
	
	static bool instanceFlag;
	static TrustStatistics *trustStaticsUtil;
	
	AverageTrust** malicious_avg_trust; //timestamp -> trust
	AverageTrust** normal_avg_trust; //timestamp -> trust
	AverageTrust** general_avg_trust;
	AverageTrust*  nodeToAnalyze_avg_Trust;

	int max_t_stamp;
	
	BadContent* maliciousContent;
	int num_cont;
	
	int bad_content;
	int total_bad_content;
	int aux_bad_content;
	int aux_total_bad_content;

	int nodeToAnalyze;
	

public:
	static TrustStatistics* getInstance(int nNodes);
	static TrustStatistics* getInstance();
	TrustStatistics(int nNodes){
		m_nNodes = nNodes;
		max_t_stamp = 110;
		num_cont = 0;
		bad_content = 0;
		total_bad_content = 0;
		aux_bad_content = 0;
		aux_total_bad_content = 0;

		this->malicious_avg_trust = new AverageTrust*[nNodes];
		for (int i = 0; i < nNodes; i++){
			this->malicious_avg_trust[i] = new AverageTrust[max_t_stamp];
			for (int j = 0; j < max_t_stamp; j+=2) {
				setAvgTrustOfMalNode(i, j); // Avg trust of malicioues nodes
			}
		}

		this->normal_avg_trust = new AverageTrust*[nNodes];
		for (int i = 0; i < nNodes; i++){
			this->normal_avg_trust[i] = new AverageTrust[max_t_stamp];
			for (int j = 0; j < max_t_stamp; j+=2) {
				setAvgTrustOfNormalNode(i, j); // Avg trust of normal nodes
			}
		}

		this->general_avg_trust = new AverageTrust*[nNodes];
		for (int i = 0; i < nNodes; i++){
			this->general_avg_trust[i] = new AverageTrust[max_t_stamp];
			for (int j = 0; j < max_t_stamp; j+=2) {
				setGeneralAvgTrustOfNodes(i, j); // Avg trust of normal nodes
			}
		}
		
		this->maliciousContent = new BadContent[max_t_stamp];
		for (int j = 0; j < max_t_stamp; j+=2) {
			maliciousContent[j].timestamp = j;
			maliciousContent[j].bad_content = 0;
			maliciousContent[j].total_bad_content = 0;
			maliciousContent[j].badContentDeliveryRate = ((double)aux_bad_content)/((double)num_cont);
			maliciousContent[j].badContentDeliverySuccessRate = ((double)aux_bad_content)/((double)aux_total_bad_content);
		}
		
	}

	void setAvgTrustOfMalNode(int node, int timestamp){
		Trust* trustInstance = Trust::getInstance(m_nNodes);
		this->malicious_avg_trust[node][timestamp].timestamp = (double) timestamp;
		this->malicious_avg_trust[node][timestamp].avg_trust = trustInstance->avgTrust(node);
		this->malicious_avg_trust[node][timestamp].avg_contact_trust = trustInstance->avgContactTrust(node);
	}
	void setAvgTrustOfNormalNode(int node, int timestamp){
		Trust* trustInstance = Trust::getInstance(m_nNodes);
		this->normal_avg_trust[node][timestamp].timestamp = (double) timestamp;
		this->normal_avg_trust[node][timestamp].avg_trust = trustInstance->avgTrust(node);
		this->normal_avg_trust[node][timestamp].avg_contact_trust = trustInstance->avgContactTrust(node);
	}
	void setGeneralAvgTrustOfNodes(int node, int timestamp){
		Trust* trustInstance = Trust::getInstance(m_nNodes);
		this->general_avg_trust[node][timestamp].timestamp = (double) timestamp;
		this->general_avg_trust[node][timestamp].avg_trust = trustInstance->avgTrust(node);
		this->general_avg_trust[node][timestamp].avg_contact_trust = trustInstance->avgContactTrust(node);
	}


	void setNodeToAnalyzeTrust(int timestamp){
		Trust* trustInstance = Trust::getInstance(m_nNodes);
		this->nodeToAnalyze_avg_Trust[timestamp].timestamp = (double) timestamp;
		this->nodeToAnalyze_avg_Trust[timestamp].avg_trust = trustInstance->avgTrust(nodeToAnalyze);
		this->nodeToAnalyze_avg_Trust[timestamp].avg_contact_trust = trustInstance->avgContactTrust(nodeToAnalyze);
	}

	void setNodeToAnalyze(int node){
		nodeToAnalyze = node;
		this->nodeToAnalyze_avg_Trust = new AverageTrust[max_t_stamp];
		for (int j = 0; j < max_t_stamp; j+=2) {
			setNodeToAnalyzeTrust(j); // Avg trust of normal nodes
		}
	}

	void print(){
		
		for (int j = 0; j <= max_t_stamp; j+=2) {
			std::cout << "#--#Timestamp: (" << j << " s)" << std::endl;
			std::cout << "\t-> Average malicious nodes trust: " << getAvgMaliciousNodesTrust(j) << std::endl;
			std::cout << "\t-> Average normal nodes trust: " << getAvgNormalNodesTrust(j) << std::endl;
			std::cout << "\t-> Average General nodes trust: " << getAvgGeneralNodesTrust(j) << std::endl;
			std::cout << "\t-> Average contacted malicious nodes trust: " << getContactAvgMaliciousNodesTrust(j) << std::endl;
			std::cout << "\t-> Average contacted normal nodes trust: " << getContactAvgNormalNodesTrust(j) << std::endl;
			std::cout << "\t-> Average contacted General nodes trust: " << getContactAvgGeneralNodesTrust(j) << std::endl;
			std::cout << "\t-> Bad content successfully delivered: " << getBadContentNumber(j) << std::endl;
			std::cout << "\t-> Bad content successfully delivered rate: " << getBadContentSuccessRate(j) << std::endl;
			std::cout << "\t-> Bad content totality: " << getTotalBadContentNumber(j) << std::endl;
			std::cout << "\t-> Bad content delivery rate: " << getBadContentDeliveryRate(j) << std::endl;
			std::cout << "\t-> Bad content delivery rate: " << getBadContentDeliveryRate(j) << std::endl;
			std::cout << "\t-> Average malicious node " << nodeToAnalyze << " trust: " << getAvgNodeToAnalyzeTrust(j) << std::endl;
		}
	}

	void print(int timestamp){
		std::cout << "#--#Timestamp: (" << timestamp << " s)" << std::endl;
		std::cout << "\t-> Average malicious nodes trust: " << getAvgMaliciousNodesTrust(timestamp) << std::endl;
		std::cout << "\t-> Average normal nodes trust: " << getAvgNormalNodesTrust(timestamp) << std::endl;
		std::cout << "\t-> Average General nodes trust: " << getAvgGeneralNodesTrust(timestamp) << std::endl;
		std::cout << "\t-> Average contacted malicious nodes trust: " << getContactAvgMaliciousNodesTrust(timestamp) << std::endl;
		std::cout << "\t-> Average contacted normal nodes trust: " << getContactAvgNormalNodesTrust(timestamp) << std::endl;
		std::cout << "\t-> Average contacted General nodes trust: " << getContactAvgGeneralNodesTrust(timestamp) << std::endl;
		std::cout << "\t-> Bad content successfully delivered: " << getBadContentNumber(timestamp) << std::endl;
		std::cout << "\t-> Bad content successfully delivered rate: " << getBadContentSuccessRate(timestamp) << std::endl;
		std::cout << "\t-> Bad content totality: " << getTotalBadContentNumber(timestamp) << std::endl;
		std::cout << "\t-> Bad content delivery rate: " << getBadContentDeliveryRate(timestamp) << std::endl;
		std::cout << "\t-> Average malicious node " << nodeToAnalyze << " trust: " << getAvgNodeToAnalyzeTrust(timestamp) << std::endl;
	}

	int getNodeToAnalyze(){
		return nodeToAnalyze;
	}

	double getAvgNodeToAnalyzeTrust(int timestamp){
		return this->nodeToAnalyze_avg_Trust[timestamp].avg_contact_trust;
	}

	double getAvgNormalNodesTrust(int timestamp){
		double sum = 0.0;
		for (int node = 0; node < m_nNodes; node++){
			sum += this->normal_avg_trust[node][timestamp].avg_trust;
		}
		return sum/m_nNodes;
	}
	
	double getAvgMaliciousNodesTrust(int timestamp){
		double sum = 0.0;
		for (int node = 0; node < m_nNodes; node++){
			sum += this->malicious_avg_trust[node][timestamp].avg_trust;
		}
		return sum/m_nNodes;
	}

	double getAvgGeneralNodesTrust(int timestamp){
		double sum = 0.0;
		for (int node = 0; node < m_nNodes; node++){
			sum += this->general_avg_trust[node][timestamp].avg_trust;
		}
		return sum/m_nNodes;
	}

	double getContactAvgNormalNodesTrust(int timestamp){
		double sum = 0.0;
		for (int node = 0; node < m_nNodes; node++){
			sum += this->normal_avg_trust[node][timestamp].avg_contact_trust;
		}
		return sum/m_nNodes;
	}
	
	double getContactAvgMaliciousNodesTrust(int timestamp){
		double sum = 0.0;
		for (int node = 0; node < m_nNodes; node++){
			sum += this->malicious_avg_trust[node][timestamp].avg_contact_trust;
		}
		return sum/m_nNodes;
	}

	double getContactAvgGeneralNodesTrust(int timestamp){
		double sum = 0.0;
		for (int node = 0; node < m_nNodes; node++){
			sum += this->general_avg_trust[node][timestamp].avg_contact_trust;
		}
		return sum/m_nNodes;
	}

	void newBadContentTransmission(){
		this->bad_content++;	
		this->aux_bad_content++;	
	}
	void newInvalidContent(){
		this->total_bad_content++;
		this->aux_total_bad_content++;
	}

	void resetCounters(){
		this->aux_bad_content = 0;
		this->aux_total_bad_content = 0;
	}
	void setBadContentDeliveryRate(int timestamp){
		this->maliciousContent[timestamp].timestamp = (double) timestamp;
		this->maliciousContent[timestamp].bad_content = aux_bad_content;
		this->maliciousContent[timestamp].total_bad_content = aux_total_bad_content;
		this->maliciousContent[timestamp].badContentDeliveryRate = (double) aux_bad_content/num_cont;
		this->maliciousContent[timestamp].badContentDeliverySuccessRate = (double) aux_bad_content/aux_total_bad_content;
	}
	int getBadContentNumber(int timestamp){
		return this->maliciousContent[timestamp].bad_content;
	}
	int getTotalBadContentNumber(int timestamp){
		return this->maliciousContent[timestamp].total_bad_content;
	}
	double getBadContentDeliveryRate(int timestamp){
		return this->maliciousContent[timestamp].badContentDeliveryRate;
	}
	double getBadContentSuccessRate(int timestamp){
		return this->maliciousContent[timestamp].badContentDeliverySuccessRate;
	}
	void newContent(){
		num_cont++;
	}

};

bool TrustStatistics::instanceFlag = false;
TrustStatistics* TrustStatistics::trustStaticsUtil = NULL;

TrustStatistics* TrustStatistics::getInstance(int nNodes){
	if (!instanceFlag) {
		trustStaticsUtil = new TrustStatistics(nNodes);
		instanceFlag = true;
	}
	return trustStaticsUtil;
}

TrustStatistics* TrustStatistics::getInstance(){
	if (instanceFlag) 
		return trustStaticsUtil;
	exit (EXIT_FAILURE);
}
