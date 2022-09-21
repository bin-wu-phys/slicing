#ifndef analysis_qT_h
#define analysis_qT_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"
#include "TH1.h"

namespace MA5
{
class qT : public AnalyzerBase
{
  INIT_ANALYSIS(qT,"qT")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

  double DeltaR(const MCParticleFormat*, const MCParticleFormat*);//distance of two particles in (phi, eta)
  double DeltaPhi(const MCParticleFormat*, const MCParticleFormat*);//distance of two particles in (phi, eta)


  double PTVecSum(const MCParticleFormat*, const MCParticleFormat*);//the modulus of the total pT of two particles

  double qTInJet(double, double, double);//In jet qT

  bool selectQ(double pJ, double etaJ);//test selected kinetics.

  void setParas(double, double, double);//set parameters _R, _etaMax, _pTJMin

  void setInitAll();//choose all the channels
  void setInitPartons(int, int);//set the PDF ids of the initial partons: _iA, _iB

  void printSummary(const SampleFormat& summary);//Print out the information in Sample 

  void setLumi(double);//set luminosity in the unit of 1/pb

  void cmpWithNorm(const SampleFormat& summary);//compare with the results of the normal mode.

  void dsdqT(const SampleFormat& summary);//output dsigma/dq_T
  void sigma(const SampleFormat& summary);//output \int d q_T dsigma/dq_T
  void sigma(const SampleFormat& summary, TH1F* hist, const char* fname, const char* xlabel, const char* ylabel);//output \int d q_T dsigma/dq_T

  
 private:
  //Luminosity
  double _L;
  
  //parameters
  double _R; //jet radius
  double _etaMax;//maximum pseudorapidity
  double _pTJMin;//lower jet pT cut

  //swithces
  bool _initPartonsQ;//pick a pair of initial partons?
  int _iA, _iB;//PDG id for the two initial partons
  
  TH1F *_histDphi;//Delta phi histogram
  TH1F *_histqT;//qT histogram
  TH1F *_histqTSJA;//qT histogram in SJA

  void SortpT(const MCParticleFormat**, int *idx);

  bool InJetQ(const MCParticleFormat*, const MCParticleFormat*);//test whether j2 and j3 in the same jet
};
}

#endif
