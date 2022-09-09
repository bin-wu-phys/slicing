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

 private:
  TH1F *_histDphi;//Delta phi histogram
  TH1F *_histqT;//qT histogram

  void SortpT(const MCParticleFormat**);
  void SortpT(const MCParticleFormat**, int *idx);
};
}

#endif