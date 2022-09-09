#include "SampleAnalyzer/User/Analyzer/qT.h"
#include "TCanvas.h"
#include "TMath.h"
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool qT::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;
  // Initializing PhysicsService for MC
  PHYSICS->mcConfig().Reset();

  //set parameters
  _R = 0.4;
  
  // Initializing histograms
  _histDphi = new TH1F("#Delta#phi", "#Delta#phi_{12}", 50, 0.0, TMath::Pi());
  _histDphi->GetXaxis()->SetTitle("#Delta#phi");
  _histDphi->GetYaxis()->SetTitle("Events  ( L_{int} = 10 fb^{-1} )");
  
  _histqT = new TH1F("q_{T}", "q_{T}", 50, 0.0, 50);
  _histqT->GetXaxis()->SetTitle("q_{T}[GeV]");
  _histqT->GetYaxis()->SetTitle("Events  ( L_{int} = 10 fb^{-1} )");

  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void qT::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  // Normalization of the histogram: L = 10 fb-1
  double nrm = summary.mc()->xsection() * 10000. /
    static_cast<float>(summary.nevents());
  _histDphi->Scale(nrm); _histqT->Scale(nrm);
  
  //Output
  TCanvas* c1 = new TCanvas("c1","#Delta#phi", 500, 700);
  //_histDphi->SetFillColor(kRed);
  c1->SetLeftMargin(0.14);
  _histDphi->Draw("HIST");
  c1->SaveAs("Dphi.pdf");

  TCanvas* c2 = new TCanvas("c2","q_{T}", 500, 700);
  c2->SetLeftMargin(0.14);
  //_histDphi->SetFillColor(kRed);
  _histqT->Draw("HIST");
  c2->SaveAs("qT.pdf");  
  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool qT::Execute(SampleFormat& sample, const EventFormat& event)
{
  const MCParticleFormat *Js[3]; int iJ = 0;
  for (MAuint32 i=0;i<event.mc()->particles().size();i++){
    const MCParticleFormat* part = &event.mc()->particles()[i];

    if (PHYSICS->Id->IsFinalState(*part)){
      Js[iJ] = part;
      iJ++;
    }
  }
  if(iJ!=3){
    cout << "There are " << iJ << " final-state particles?" << endl;
  }else{
    //Sort the pT
    int idx[3];
    SortpT(Js, idx);

    double pJ, etaJ, qT;
    if(InJetQ(Js[idx[1]], Js[idx[2]])){
      pJ = Js[idx[1]]->pt() + Js[idx[2]]->pt();
      etaJ = Js[idx[1]]->eta(); qT = qTInJet(Js[idx[0]]->pt(), Js[idx[1]]->pt(), Js[idx[2]]->pt());
    }else{
      pJ = Js[idx[0]]->pt(); etaJ = Js[idx[0]]->eta(); qT = PTVecSum(Js[idx[0]], Js[idx[1]]);
    }
    if(selectQ(pJ, etaJ)){
      _histDphi->Fill(DeltaPhi(Js[idx[0]], Js[idx[1]]));
      _histqT->Fill(qT);
    }
  }
  
  return true;
}

void qT::SortpT(const MCParticleFormat** Js, int *idx){
  int iMin=0, iMax=0;
  double pMin = Js[0]->pt(), pMax = Js[0]->pt();
  for(int i=1;i<3;i++){
    if(Js[i]->pt() < pMin){
      iMin = i; pMin = Js[i]->pt();
    }else if(Js[i]->pt() > pMax){
      iMax = i; pMax = Js[i]->pt();
    }
  }

  if(iMax!=iMin){
    idx[0] = iMax; idx[2] = iMin;
    for(int i=0;i<3;i++){
      if((i!=iMax)&&(i!=iMin)) idx[1] = i;
    }
  }else{
    for(int i=0;i<3;i++)
      idx[i] = i;
  }
}

double qT::DeltaR(const MCParticleFormat* p1, const MCParticleFormat* p2){
  double deta = p1->eta()-p2->eta(), dphi = p1->phi()-p2->phi();
  return sqrt(deta*deta + dphi*dphi);
}

double qT::PTVecSum(const MCParticleFormat* p1, const MCParticleFormat* p2){
  double dpx = p1->px()+p2->px(), dpy = p1->py()+p2->py();
  return sqrt(dpx*dpx + dpy*dpy);
}

double qT::qTInJet(double p1, double p2, double p3){
  return sqrt((p1+p2+p3)*(p2+p3-p1)*p3/p2);
}

double qT::DeltaPhi(const MCParticleFormat* p1, const MCParticleFormat* p2){
  /*Because phi() in [-pi, pi], dphi in [-2pi, 2pi].
    We further constrain it to [0, pi], corresponding to DPHI_0_PI.
  */
  double dphi = fabs(p1->phi() - p2->phi());
  if(dphi>TMath::Pi()){
    dphi = 2.0*TMath::Pi() - dphi;
  }
  return dphi;
}

bool qT::InJetQ(const MCParticleFormat* p2, const MCParticleFormat* p3){
  bool Q = false;
  if(DeltaR(p2, p3) < _R) Q = true;
  return Q;
}

bool qT::selectQ(double pJ, double etaJ){
  bool Q = true;
  return Q;
}
