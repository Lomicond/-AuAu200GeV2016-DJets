#ifndef ANGULARITY_H
#define ANGULARITY_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>

// Třída Angularity, která dědí z fastjet::FunctionOfPseudoJet<double>
// a implementuje generalized angularity s exponenty a (pro delta_R) a kappa (pro pT)
class Angularity : public fastjet::FunctionOfPseudoJet<double> {
public:
  // Konstruktor s parametry:
  //   a     - exponent pro vzdálenost (delta_R)
  //   kappa - exponent pro vážení hybností (pT)
  Angularity(double kappa = 1.0, double a = 1.0, double Radius = 0.4) :  _kappa(kappa), _a(a), _Radius(Radius) {}

  // Metoda, která spočítá generalized angularity pro daný jet
  virtual double result(const fastjet::PseudoJet& jet) const {
  
  double f_result = -999;
  
    double numerator = 0.0;
    double pt_total = jet.pt();
    std::vector<fastjet::PseudoJet> constituents = jet.constituents();
    
    for (size_t i = 0; i < constituents.size(); ++i) {
      double pt_i = constituents[i].pt();
      
      // Výpočet rozdílů v pseudorapiditě a phi (eta-phi prostor)
      double deta = jet.eta() - constituents[i].eta();
      double dphi = jet.phi() - constituents[i].phi();
      
      // Normalizace úhlu phi do intervalu (-pi, pi)
      if (dphi > M_PI)  dphi -= 2 * M_PI;
      if (dphi < -M_PI) dphi += 2 * M_PI;
      
      double dr = std::sqrt(deta * deta + dphi * dphi);
      
      // Přičtení váženého příspěvku každého constituentu
      numerator += std::pow(1.*pt_i/pt_total, _kappa) * std::pow(1.0*dr/_Radius, _a);
    }
    
    
    f_result = numerator;
    
    //cout << numerator << "tohle je to spocitany: " << endl;
    return f_result;
    
  }

  // Popis funkce pro ladění či výpis informací
  virtual std::string description() const {
    std::ostringstream oss;
    oss << "Generalized Angularity with parameters kappa = " << _kappa << " and alpha = " << _a << " and jet R = " << _Radius;
    return oss.str();
  }
  
private:
  double _kappa; // exponent pro vážení hybností (pT)
  double _a;     // exponent pro vzdálenost delta_R
  double _Radius; // jet radius
};

#endif // ANGULARITY_H

