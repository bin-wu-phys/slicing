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
  
  _histqT = new TH1F("q_{T}", "#Delta#phi_{23}", 50, 0.0, TMath::Pi());
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
  TH1F* S1_DPHI_0_PI_0 = new TH1F("S1_DPHI_0_PI_0","S1_DPHI_0_PI_0",50,0.0,3.141593);
  // Content
  S1_DPHI_0_PI_0->SetBinContent(0,0.0); // underflow
  S1_DPHI_0_PI_0->SetBinContent(1,0.0);
  S1_DPHI_0_PI_0->SetBinContent(2,0.0);
  S1_DPHI_0_PI_0->SetBinContent(3,0.0);
  S1_DPHI_0_PI_0->SetBinContent(4,0.0);
  S1_DPHI_0_PI_0->SetBinContent(5,0.0);
  S1_DPHI_0_PI_0->SetBinContent(6,0.0);
  S1_DPHI_0_PI_0->SetBinContent(7,0.0);
  S1_DPHI_0_PI_0->SetBinContent(8,0.0);
  S1_DPHI_0_PI_0->SetBinContent(9,0.0);
  S1_DPHI_0_PI_0->SetBinContent(10,0.0);
  S1_DPHI_0_PI_0->SetBinContent(11,0.0);
  S1_DPHI_0_PI_0->SetBinContent(12,0.0);
  S1_DPHI_0_PI_0->SetBinContent(13,0.0);
  S1_DPHI_0_PI_0->SetBinContent(14,0.0);
  S1_DPHI_0_PI_0->SetBinContent(15,0.0);
  S1_DPHI_0_PI_0->SetBinContent(16,0.0);
  S1_DPHI_0_PI_0->SetBinContent(17,0.0);
  S1_DPHI_0_PI_0->SetBinContent(18,0.0);
  S1_DPHI_0_PI_0->SetBinContent(19,0.0);
  S1_DPHI_0_PI_0->SetBinContent(20,0.0);
  S1_DPHI_0_PI_0->SetBinContent(21,0.0);
  S1_DPHI_0_PI_0->SetBinContent(22,0.0);
  S1_DPHI_0_PI_0->SetBinContent(23,0.0);
  S1_DPHI_0_PI_0->SetBinContent(24,0.0);
  S1_DPHI_0_PI_0->SetBinContent(25,0.0);
  S1_DPHI_0_PI_0->SetBinContent(26,0.0);
  S1_DPHI_0_PI_0->SetBinContent(27,0.0);
  S1_DPHI_0_PI_0->SetBinContent(28,0.0);
  S1_DPHI_0_PI_0->SetBinContent(29,0.0);
  S1_DPHI_0_PI_0->SetBinContent(30,0.0);
  S1_DPHI_0_PI_0->SetBinContent(31,0.0);
  S1_DPHI_0_PI_0->SetBinContent(32,0.0);
  S1_DPHI_0_PI_0->SetBinContent(33,0.0);
  S1_DPHI_0_PI_0->SetBinContent(34,2672747008.399999);
  S1_DPHI_0_PI_0->SetBinContent(35,13268280041.699995);
  S1_DPHI_0_PI_0->SetBinContent(36,21954710069.0);
  S1_DPHI_0_PI_0->SetBinContent(37,27491110086.399982);
  S1_DPHI_0_PI_0->SetBinContent(38,33313880104.699978);
  S1_DPHI_0_PI_0->SetBinContent(39,35127530110.39998);
  S1_DPHI_0_PI_0->SetBinContent(40,33886620106.500008);
  S1_DPHI_0_PI_0->SetBinContent(41,39232110123.29999);
  S1_DPHI_0_PI_0->SetBinContent(42,33027520103.8);
  S1_DPHI_0_PI_0->SetBinContent(43,35986630113.09999);
  S1_DPHI_0_PI_0->SetBinContent(44,34077530107.100006);
  S1_DPHI_0_PI_0->SetBinContent(45,30832050096.9);
  S1_DPHI_0_PI_0->SetBinContent(46,29782040093.599995);
  S1_DPHI_0_PI_0->SetBinContent(47,30259320095.100006);
  S1_DPHI_0_PI_0->SetBinContent(48,33504790105.299976);
  S1_DPHI_0_PI_0->SetBinContent(49,50495830158.69999);
  S1_DPHI_0_PI_0->SetBinContent(50,469639901476.00006);
  S1_DPHI_0_PI_0->SetBinContent(51,0.0); // overflow
  S1_DPHI_0_PI_0->SetEntries(10000);

  c1->SetLeftMargin(0.14);
  _histDphi->SetLineColor(kRed);
  _histDphi->Draw("HIST");
  S1_DPHI_0_PI_0->SetLineColor(kBlue);
  S1_DPHI_0_PI_0->SetLineStyle(2);
  S1_DPHI_0_PI_0->Draw("SAME");
  c1->SaveAs("Dphi12.pdf");

  TCanvas* c2 = new TCanvas("c2","q_{T}", 500, 700);
    TH1F* S2_DPHI_0_PI_0 = new TH1F("S2_DPHI_0_PI_0","S2_DPHI_0_PI_0",50,0.0,3.141593);
  // Content
  S2_DPHI_0_PI_0->SetBinContent(0,0.0); // underflow
  S2_DPHI_0_PI_0->SetBinContent(1,420384925411.0776);
  S2_DPHI_0_PI_0->SetBinContent(2,36654822215.680145);
  S2_DPHI_0_PI_0->SetBinContent(3,23195631402.110195);
  S2_DPHI_0_PI_0->SetBinContent(4,21859251321.329803);
  S2_DPHI_0_PI_0->SetBinContent(5,17181951038.600254);
  S2_DPHI_0_PI_0->SetBinContent(6,14891020900.12002);
  S2_DPHI_0_PI_0->SetBinContent(7,13363740807.800266);
  S2_DPHI_0_PI_0->SetBinContent(8,11931910721.250195);
  S2_DPHI_0_PI_0->SetBinContent(9,11263720680.859999);
  S2_DPHI_0_PI_0->SetBinContent(10,10881900657.780062);
  S2_DPHI_0_PI_0->SetBinContent(11,13841010836.649887);
  S2_DPHI_0_PI_0->SetBinContent(12,10690990646.240091);
  S2_DPHI_0_PI_0->SetBinContent(13,11454630692.39997);
  S2_DPHI_0_PI_0->SetBinContent(14,11550090698.170256);
  S2_DPHI_0_PI_0->SetBinContent(15,9450070571.22999);
  S2_DPHI_0_PI_0->SetBinContent(16,11168260675.089712);
  S2_DPHI_0_PI_0->SetBinContent(17,13268280802.02998);
  S2_DPHI_0_PI_0->SetBinContent(18,13363740807.800266);
  S2_DPHI_0_PI_0->SetBinContent(19,11263720680.859999);
  S2_DPHI_0_PI_0->SetBinContent(20,13936470842.420174);
  S2_DPHI_0_PI_0->SetBinContent(21,15081930911.659988);
  S2_DPHI_0_PI_0->SetBinContent(22,17181951038.600254);
  S2_DPHI_0_PI_0->SetBinContent(23,17372861050.140224);
  S2_DPHI_0_PI_0->SetBinContent(24,17563771061.680195);
  S2_DPHI_0_PI_0->SetBinContent(25,18995601148.230267);
  S2_DPHI_0_PI_0->SetBinContent(26,18231951102.069782);
  S2_DPHI_0_PI_0->SetBinContent(27,24913821505.969917);
  S2_DPHI_0_PI_0->SetBinContent(28,24245641465.580326);
  S2_DPHI_0_PI_0->SetBinContent(29,23672901430.959816);
  S2_DPHI_0_PI_0->SetBinContent(30,26822931621.370213);
  S2_DPHI_0_PI_0->SetBinContent(31,22527441361.719997);
  S2_DPHI_0_PI_0->SetBinContent(32,16131940975.130121);
  S2_DPHI_0_PI_0->SetBinContent(33,9354615565.460005);
  S2_DPHI_0_PI_0->SetBinContent(34,859097351.9300007);
  S2_DPHI_0_PI_0->SetBinContent(35,0.0);
  S2_DPHI_0_PI_0->SetBinContent(36,0.0);
  S2_DPHI_0_PI_0->SetBinContent(37,0.0);
  S2_DPHI_0_PI_0->SetBinContent(38,0.0);
  S2_DPHI_0_PI_0->SetBinContent(39,0.0);
  S2_DPHI_0_PI_0->SetBinContent(40,0.0);
  S2_DPHI_0_PI_0->SetBinContent(41,0.0);
  S2_DPHI_0_PI_0->SetBinContent(42,0.0);
  S2_DPHI_0_PI_0->SetBinContent(43,0.0);
  S2_DPHI_0_PI_0->SetBinContent(44,0.0);
  S2_DPHI_0_PI_0->SetBinContent(45,0.0);
  S2_DPHI_0_PI_0->SetBinContent(46,0.0);
  S2_DPHI_0_PI_0->SetBinContent(47,0.0);
  S2_DPHI_0_PI_0->SetBinContent(48,0.0);
  S2_DPHI_0_PI_0->SetBinContent(49,0.0);
  S2_DPHI_0_PI_0->SetBinContent(50,0.0);
  S2_DPHI_0_PI_0->SetBinContent(51,0.0); // overflow
  S2_DPHI_0_PI_0->SetEntries(10000);

  c2->SetLeftMargin(0.14);
  _histqT->SetLineColor(kRed);
  _histqT->Draw("HIST");
  S2_DPHI_0_PI_0->SetLineColor(kBlue);
  S2_DPHI_0_PI_0->SetLineStyle(2);
  S2_DPHI_0_PI_0->Draw("SAME");
  c2->SaveAs("Dphi23.pdf");  
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
    _histDphi->Fill(DeltaPhi(Js[idx[0]], Js[idx[1]]));
    _histqT->Fill(DeltaPhi(Js[idx[1]], Js[idx[2]]));
    /*
    if(InJetQ(Js[idx[1]], Js[idx[1]])){
    }else{
      _histDphi->Fill(Js[idx[0]]->pt());
      _histqT->Fill(Js[idx[2]]->pt());
    }
    */
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

