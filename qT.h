#ifndef analysis_qT_h
#define analysis_qT_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"
#include "TH1.h"

enum { _g =21, _d = 1, _db = -1,  _u =2, _ub = -2};


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
  double rapidity(const MCParticleFormat*, const MCParticleFormat*);//the modulus of the total pT of two particles

  double qTInJet(double, double, double);//In jet qT

  bool selectetaQ(double etaJ);//test selected kinetics.
  bool selectpTQ(double pJ);//test selected kinetics.

  void setParas(double, double, double);//set parameters _R, _etaMax, _pTJMin

  void setInitAll();//choose all the channels
  void setInitPartons(int, int);//set the PDF ids of the initial partons: _iA, _iB

  void printSummary(const SampleFormat& summary);//Print out the information in Sample 

  void setLumi(double);//set luminosity in the unit of 1/pb

  void cmpWithNorm(const SampleFormat& summary);//compare with the results of the normal mode.

  void dsdqT(const SampleFormat& summary);//output dsigma/dq_T
  void dsdphi(const SampleFormat& summary);//output dsigma/ddphi
  void sigma(const SampleFormat& summary, TH1F* hist, const char* fname, const char* xlabel, const char* ylabel);//output \int d q_T dsigma/dq_T

  bool chanelud();//select a chanel
  bool chanelgg2ggg();//select a chanel
  bool chanelall(){return true;};//no constrain
  bool setCut(double *pJ, double *etaJ, unsigned int nJ);

  void sigtot(const SampleFormat& summary);
 private:
  double _aS, _scale;
  unsigned int _pdfIDA, _pdfIDB;

  //make selections
  int _initIDs[2];//make sure there are only two initial-state particles.
  int _finalIDs[3];

  //Luminosity
  double _L;

  //count total # of events: just want to compare with that in the summary
  unsigned _numEvents, _numSelected;
  
  //parameters
  double _R; //jet radius
  double _etaMax;//maximum pseudorapidity
  double _pTJMin;//lower jet pT cut

  //swithces
  bool _initPartonsQ;//pick a pair of initial partons?
  int _iA, _iB;//PDG id for the two initial partons
  
  TH1F *_histDphi;//Delta phi histogram
  TH1F *_histdphi;//delta phi histogram
  TH1F *_histldphi;//log10(dphi) histogram
  TH1F *_histDphiS;//Delta phi histogram with SJA
  TH1F *_histdphiS;//delta phi histogram with SJA
  TH1F *_histldphiS;//log10(dphi) histogram with SJA

  TH1F *_histqT;//qT histogram using WTA
  TH1F *_histqTSJA;//qT histogram

  void SortpT(const MCParticleFormat**, int *idx);

  bool InJetQ(const MCParticleFormat*, const MCParticleFormat*);//test whether j2 and j3 in the same jet

  //total cross section
  double _sigWTA, _sigSJA;
};
}

#endif
