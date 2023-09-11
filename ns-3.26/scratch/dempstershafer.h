#include <set>
#include <list>
#include <bitset>
#include <stdarg.h>
#include <iostream>
#include <cstddef>
#include <stdarg.h>
#include <iostream>
#include <cstddef>
#include <iostream>

#define MAX_HYP 30

using namespace std;

class DSUniverse;
class Evidence;
class DempsterShafer;

string hypothesis_to_string_function(void* element);
string trust("trust");
string untrust("untrust");

typedef struct {
  bitset <MAX_HYP> items;
  double mass;
} FocalSet;


class DSUniverse{
private:
  void* hyps[MAX_HYP]; // Conjunto de hipóteses
  int last_hyp_numb; // Índice da última hipótese adicionada

public:

  DSUniverse();
  void add_hyp(set<void*> &hyps);
  void add_hyp(void* hyp, ...);
  Evidence add_evidence();
  bitset<MAX_HYP> instantiate_bitset(set<void*>& hyp_members);
  bitset<MAX_HYP> instantiate_bitset(void* member, ...);
  friend class Evidence;
};

//universe->last_hyp_numb
class Evidence{
private:
  DSUniverse *universe;
  Evidence(DSUniverse *universe);
  void add_focal_element(FocalSet set);
  list<FocalSet> focal_sets;
public:
  void add_focal_element (double mass, bitset<MAX_HYP> &items);
  void add_focal_element (double mass, set<void*> &items);
  void add_focal_element (double mass, void* item, ...);
  void add_omega_set();
  Evidence operator&(Evidence &otherEvidence);
  double conflict(Evidence &other);
  double belief(bitset<MAX_HYP> &items);
  double belief(set<void*> &items);
  double belief(void *item, ...);
  double plausability(bitset<MAX_HYP> &items);
  double plausability(set<void*> &items);
  double plausability(void *item, ...);
  void* best_hyp();
  friend class DSUniverse;
};



DSUniverse::DSUniverse():
  last_hyp_numb(0){}

Evidence DSUniverse::add_evidence() {
  return Evidence(this);
}

void DSUniverse::add_hyp(set<void*> &hyps) {
  // Percorrer os elementos do conjunto de hipóteses até o último elemento.
  // A hipótese será adicionada no final do conjunto
  for (set<void*>::iterator i = hyps.begin(); i != hyps.end(); i++) {
    if(last_hyp_numb >= MAX_HYP - 1) {
      throw 1;
    }
    this->hyps[last_hyp_numb++] = *i;
  }
}

void DSUniverse::add_hyp(void* hyp, ...) {
  va_list arg_list;
  void *present;
  set<void*> hyps;

  va_start(arg_list, hyp);
  for (present = hyp; present != NULL; present = va_arg(arg_list, void*)) {
    hyps.insert(present);
  }
  va_end(arg_list);

  add_hyp(hyps);
}

bitset<MAX_HYP> DSUniverse::instantiate_bitset(set<void*>& hyp_members) {
  bitset<MAX_HYP> newHyp;
  for (set<void*>::iterator i = hyp_members.begin(); i != hyp_members.end(); i++) {
    for (int j = 0; j < last_hyp_numb; j++) {
      if(hyps[j] == *i) {
        newHyp.set(j, true);
      }
    }
  }
  return newHyp;
}

bitset<MAX_HYP> DSUniverse::instantiate_bitset(void* member, ...) {
  va_list arg_list;
  void *present;
  set<void*> hyp_members;

  va_start(arg_list, member);
  for (present = member; present != NULL; present = va_arg(arg_list, void*)) {
    hyp_members.insert(present);
  }
  va_end(arg_list);

  return instantiate_bitset(hyp_members);
}

Evidence::Evidence(DSUniverse *universe):
  universe(universe){}

void  Evidence::add_focal_element(FocalSet set){
  focal_sets.push_back(set);
}

void Evidence::add_focal_element (double mass, bitset<MAX_HYP> &items) {
  FocalSet set;
  set.mass = mass;
  set.items = items;

  focal_sets.push_back(set);
}

void Evidence::add_focal_element (double mass, set<void*> &items) {
  FocalSet set;
  set.mass = mass;
  set.items = universe->instantiate_bitset(items);

  focal_sets.push_back(set);
}

void Evidence::add_focal_element (double mass, void* item, ...) {
  set<void*> items;
  va_list arg_list;
  void* present;

  va_start(arg_list, item);
  for (present = item; present != NULL; present = va_arg(arg_list, void*)) {
    items.insert(present);
  }
  va_end(arg_list);

  add_focal_element(mass, items);
}

void Evidence::add_omega_set() {
  double total_mass = 0.0;
  for (list<FocalSet>::iterator i = focal_sets.begin(); i != focal_sets.end(); i++) {
    total_mass = total_mass + i->mass;
  }
  if (total_mass < 1.0){
    FocalSet omegaSet;
    omegaSet.mass = 1.0 - total_mass; // Massa restante para o omegaSet
    omegaSet.items.set(); // Passa os bits relativos aos items no bitmap para 1
    focal_sets.push_back(omegaSet);
  }
}

// Sobrepõe o operador & para junção de evidências
Evidence Evidence::operator&(Evidence &otherEvidence) {
  Evidence combinedEvidence(universe);
  double k = 0.0; // Variável de conflito
  list<FocalSet> buffer;

  for (list<FocalSet>::iterator i = focal_sets.begin(); i != focal_sets.end(); i++) {
    for (list<FocalSet>::iterator j = otherEvidence.focal_sets.begin(); j != otherEvidence.focal_sets.end(); j++) {
      FocalSet combine;
      // Combina as massas de acordo com a teoria de Dempster-Shafer
      combine.mass = i->mass * j->mass;
      // Combina os bits no mapa de bits
      combine.items = i->items & j->items;

      /* * Se os bits do mapa de bits se combinarem adicionaremos essa nova evidência
            combinada ao mapa de bits.
         * Caso contrário será adicionado o conflito (k) relativo a combinação das massas
      */

      if (combine.items.any()) {
        buffer.push_back(combine);
      } else {
        k += combine.mass; // Soma o valor de conflito
      }
    }
  }
  for (list<FocalSet>::iterator i = buffer.begin(); i != buffer.end(); i++) {
    i->mass = i->mass * 1.0/(1.0 - k);

    combinedEvidence.add_focal_element(*i);
  }

  return combinedEvidence;
}

double Evidence::conflict(Evidence &otherEvidence) {
  double conflict = 0.0;
  for (list<FocalSet>::iterator i = focal_sets.begin(); i != focal_sets.end(); i++) {
    for (list<FocalSet>::iterator j = otherEvidence.focal_sets.begin(); j != otherEvidence.focal_sets.end(); j++) {
      bitset <MAX_HYP> combined_items = (i->items & j->items); // Retornará um bitset resultante
      if (combined_items.none()) {
        conflict = conflict + (i->mass * j->mass);
      }
    }
  }
  return conflict;
}

double Evidence::belief(bitset<MAX_HYP> &items) {
  double belief = 0.0;
  for (list<FocalSet>::iterator i = focal_sets.begin(); i != focal_sets.end(); i++) {
    bitset <MAX_HYP> subSet = (items & i->items); // Operação binária dos itens passados por parâmetro
    // Caso a quantidade de bits conincidentes seja a mesma em relação aos bits do focal_sets
    // Resumidamente procuramos um sub conjunto dos itens do focal_sets
    if (subSet.count() == i->items.count()) {
      belief += i->mass;
    }
  }
  return belief;
}

double Evidence::belief(set<void*> &items) {

  bitset<MAX_HYP> itemBitset = universe->instantiate_bitset(items);
  return belief(itemBitset);
}

double Evidence::belief(void *item, ...) {
  void *present;
  va_list arg_list;
  set<void*> items;

  va_start(arg_list, item);
  for (present = item; present != NULL; present = va_arg(arg_list, void*)) {
    items.insert(present);
  }
  va_end(arg_list);

  return belief(items);
}

double Evidence::plausability(bitset<MAX_HYP> &items){
    //Verifica os conjuntos focais com intersecção aos itens passados por parâmetro
  double plausability = 0.0;
  for (list<FocalSet>::iterator i = focal_sets.begin(); i != focal_sets.end(); i++) {
    if ((items & i->items).any()){
      plausability = plausability + i->mass;
    }
  }
  return plausability;
}

double Evidence::plausability(set<void*> &items){
  bitset<MAX_HYP> itemBitset = universe->instantiate_bitset(items);

  return plausability(itemBitset);
}

double Evidence::plausability(void *item, ...){
  void *present;
  va_list arg_list;
  set<void*> items;

  va_start(arg_list, item);
  for (present = item; present != NULL; present = va_arg(arg_list, void*)) {
    items.insert(present);
  }
  va_end(arg_list);

  return plausability(items);
}

void* Evidence::best_hyp() {
  double greatBelief = 0.0;
  double greatPlausability = 0.0;
  int bestHyp = 0;
  // Percorre todas as hipóteses a procura dos elementos focais correspondentes
  for (int hyp_it = 0; hyp_it < universe->last_hyp_numb; hyp_it++){
    double belief = 0.0;
    double plausability = 0.0;
    for (list<FocalSet>::iterator i = focal_sets.begin(); i != focal_sets.end(); i++){
      if(i->items.test(hyp_it)) {
        plausability = plausability + i->mass;
        if (i->items.count() == 1)
          belief = belief + i->mass;
      }
    }
    if(((belief == greatBelief) && (plausability > greatPlausability)) || (belief > greatBelief)) {
          greatBelief = belief;
          greatPlausability = plausability;
          bestHyp = hyp_it;
    }
  }
  return universe->hyps[bestHyp];
}


class DempsterShafer{
private:
  // Hipóteses
  DSUniverse universe;
  bitset<MAX_HYP> trustProve;
  bitset<MAX_HYP> untrustProve;
  bitset<MAX_HYP> evenProve;
  std::vector<Evidence> combined_evidences;
  double trust_value;
  bool first_evidence;

public:
  
  DempsterShafer(){
    
    universe.add_hyp(&trust, &untrust, NULL);
    trustProve = universe.instantiate_bitset(&trust, NULL);
    untrustProve = universe.instantiate_bitset(&untrust, NULL);
    evenProve = universe.instantiate_bitset(&trust, &untrust, NULL);
    trust_value = 0.5;
    first_evidence = true;
  }
  
  

  // void addEvidences(std::vector<double> trustValues, std::vector<bool> opinions){
  //   std::vector<Evidence> evidences;
  //   for (int i = 0; i < opinions.size(); i++){
  //     evidences.push_back(universe.add_evidence());
  //   }
  //   for (int i = 0; i < opinions.size(); i++){
  //     if (opinions[i]){
  //       evidences[i].add_focal_element(trustValues[i], trustProve); 
  //       evidences[i].add_focal_element(1.0-trustValues[i], evenProve); 
  //       evidences[i].add_omega_set();
  //     } else {
  //       evidences[i].add_focal_element(trustValues[i], untrustProve); 
  //       evidences[i].add_focal_element(1.0-trustValues[i], evenProve); 
  //       evidences[i].add_omega_set();
  //     }
  //   }
  //   Evidence comb = evidences[0];
  //   for (int i = 1; i < opinions.size(); i++){
  //     comb = comb & evidences[i];
  //   } 
  //   combined_evidences = comb;
  //   std ::cout << combined_evidences.belief(&trust, NULL) << std::endl;
  // }

  void addTrustEvidence(double trustValue){
    // std::cout << "Good Behavior evidence:" << std::endl;
    Evidence new_evidence = universe.add_evidence();
    new_evidence.add_focal_element(trustValue, trustProve); // Trustable prove
    new_evidence.add_focal_element(1.0-trustValue, evenProve); // Universal prove
    new_evidence.add_omega_set();
    // std::cout << "New evidence (trust): " <<  new_evidence.belief(&trust, NULL) << std::endl;
    // std::cout << "New evidence (untrust): " <<  new_evidence.belief(&untrust, NULL) << std::endl;
    combined_evidences.push_back(new_evidence);
    
    updateTrust();
  }

  void addUntrustEvidence(double trustValue){
    // std::cout << "Bad Behavior evidence:" << std::endl;
    Evidence new_evidence = universe.add_evidence();
    new_evidence.add_focal_element(trustValue, untrustProve); // Untrustable prove
    new_evidence.add_focal_element(1.0-trustValue, evenProve); // Universal prove
    new_evidence.add_omega_set();
    // std::cout << "New evidence (trust): " <<  new_evidence.belief(&trust, NULL) << std::endl;
    // std::cout << "New evidence (untrust): " <<  new_evidence.belief(&untrust, NULL) << std::endl;
    combined_evidences.push_back(new_evidence);
  
    updateTrust();
  }

  void updateTrust(){
    std::cout << "Evidencies size: " << combined_evidences.size() << std::endl;
    Evidence comb = combined_evidences[0];
    for (int i = 1; i < (int) combined_evidences.size(); i++){
      comb = comb & combined_evidences[i];
    }
    trust_value = comb.belief(&trust, NULL);
  }

  double getTrust(){
    return trust_value;
  }
};