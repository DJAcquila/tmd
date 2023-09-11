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

class BayesianInference{
private:
	double m_alpha;
	double m_beta;
	double m_finalTrust;
	double m_punishFactor;

	void calculateTrust(){
		m_finalTrust = m_alpha / (m_alpha + m_punishFactor * m_beta);
	}

public:
	BayesianInference(double punish) : m_alpha(1.0), m_beta(1.0), m_finalTrust(0.5), m_punishFactor(punish) {}
	BayesianInference() : m_alpha(1.0), m_beta(1.0), m_finalTrust(0.5), m_punishFactor(1.0) {}

	void addGoodBehavior(double val){
		m_alpha = m_alpha + val;
		calculateTrust();
	}
	void addBadBehavior(double val){
		m_beta = m_beta + val;
		calculateTrust(); 
	}

	void changePunishFactor(double new_punish){
		m_punishFactor = new_punish;
		calculateTrust();
	}

	double getTrust(){
		return m_finalTrust;
	}
};

