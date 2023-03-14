#include "SampleAnalyzer/User/Analyzer/qT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "ResNormal.h"
#include <string>
#include <sstream>

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
  setParas(0.4, 2.0, 50.0);
  //setInitAll();//include all the initial partons
  setInitPartons(2, -2);//pick u and ubar

  //initialize Luminisity
  setLumi(1e4);

  //reset # of events
  _numEvents = 0; _numSelected = 0;
  
  // Initializing histograms
  _histDphi = new TH1F("#Delta#phi", "#Delta#phi_{12}", 50, 2.10, 3.14);//TMath::Pi());//change to sma
  _histdphi = new TH1F("#delta#phi", "#delta#phi_{12}", 50, 1e-4, 1e-2);//TMath::Pi());//change to small delta phi
  _histldphi = new TH1F("#delta#phi", "#delta#phi_{12}", 41, -4.05, 0.05);//TMath::Pi());//change to small delta phi
  //_histdphi = new TH1F("#delta#phi", "#delta#phi_{12}", 100, 2.141593, 3.131593);//TMath::Pi());//change to small delta phi
  _histdphi->SetStats(kFALSE);
  _histdphi->GetXaxis()->SetTitle("#delta#phi");

  _histDphiS = new TH1F("#Delta#phi in SJA", "#Delta#phi_{12}", 50, 2.10, 3.14);//TMath::Pi());//change to sma
  _histdphiS = new TH1F("#delta#phi in SJA", "#delta#phi_{12}", 50, 1e-4, 1e-2);//TMath::Pi());//change to small delta phi
  _histldphiS = new TH1F("#delta#phi in SJA", "#delta#phi_{12}", 41, -4.05, 0.05);//TMath::Pi());//change to small delta phi
  //_histdphi = new TH1F("#delta#phi", "#delta#phi_{12}", 100, 2.141593, 3.131593);//TMath::Pi());//change to small delta phi

  _histqT = new TH1F("q_{T}", "q_{T} in WTA", 50, 5.0, 100);
  _histqT->SetStats(kFALSE);
  _histqT->GetXaxis()->SetTitle("q_{T}[GeV]");

  _histqTSJA = new TH1F("q_{T}", "q_{T} in SJA", 50, 5.0, 100);
  _histqTSJA->SetStats(kFALSE);
  _histqTSJA->GetXaxis()->SetTitle("q_{T}[GeV]");

  _aS = 0.0; _scale = 0.0; _pdfIDA = 0; _pdfIDB = 0;

  //total cross section
  _sigWTA = 0.0; _sigSJA = 0.0;
  
  cout << "END   Initialization" << endl;
  return true;
}

void qT::setLumi(double l){
  _L = l;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void qT::printSummary(const SampleFormat& summary){
  cout << "Summary:" << endl;
  cout << "Total # of events: " << summary.nevents() << endl;
}

void qT::cmpWithNorm(const SampleFormat& summary){
  //scale the number corresponding to _L
  double nrm = summary.mc()->xsection() * _L /static_cast<float>(summary.nevents());

  TGraph g;

  TH1F* S1_DPHI_0_PI_0 = new TH1F("S1_DPHI_0_PI_0","S1_DPHI_0_PI_0", 50, 2.10, 3.14);
  TH1F* S1_DPHI_0_PI_0_S = new TH1F("S1_DPHI_0_PI_0 in SJA","S1_DPHI_0_PI_0", 50, 2.10, 3.14);
  DphiNormal(S1_DPHI_0_PI_0); DphiSJANormal(S1_DPHI_0_PI_0_S);

  cout << "Dphi using WTA:" << endl;
  for(int i=1; i<=_histDphi->GetNbinsX(); i++){//constrained to [1, GetNbinsX()]
    double Dp = _histDphi->GetBinCenter(i);
    double Dn = nrm*_histDphi->GetBinContent(i)-S1_DPHI_0_PI_0->GetBinContent(i);//S1_DPHI_0_PI_0->GetBinContent(_histdphi->GetNbinsX()+1-i);
    if( nrm*_histDphi->GetBinContent(i)!=0.0)
      g.AddPoint(Dp, Dn/(nrm*_histDphi->GetBinContent(i)));
    cout << Dp << " " << S1_DPHI_0_PI_0->GetBinCenter(i);
    cout << " " << nrm*_histDphi->GetBinContent(i) << " " <<  S1_DPHI_0_PI_0->GetBinContent(i) << " " <<  Dn;
    if(_histDphi->GetBinContent(i)!=0.0)
      cout << " " << Dn/(nrm*_histDphi->GetBinContent(i));
    cout << endl;
  }

  cout << "Dphi using SJA:" << endl;
  for(int i=1; i<=_histDphiS->GetNbinsX(); i++){//constrained to [1, GetNbinsX()]
    double Dp = _histDphiS->GetBinCenter(i);
    double Dn = nrm*_histDphiS->GetBinContent(i)-S1_DPHI_0_PI_0_S->GetBinContent(i);//S1_DPHI_0_PI_0->GetBinContent(_histdphi->GetNbinsX()+1-i);
    if( nrm*_histDphiS->GetBinContent(i)!=0.0)
      g.AddPoint(Dp, Dn/(nrm*_histDphiS->GetBinContent(i)));
    cout << Dp << " " << S1_DPHI_0_PI_0_S->GetBinCenter(i);
    cout << " " << nrm*_histDphiS->GetBinContent(i) << " " <<  S1_DPHI_0_PI_0_S->GetBinContent(i) << " " <<  Dn;
    if(_histDphiS->GetBinContent(i)!=0.0)
      cout << " " << Dn/(nrm*_histDphiS->GetBinContent(i));
    cout << endl;
  }

  g.GetYaxis()->SetTitle("Events  ( L_{int} = 10 fb^{-1} )");
  g.GetXaxis()->SetTitle("#Delta#phi");
  
  //Output
  TCanvas* c1 = new TCanvas("c1","#Delta#phi", 500, 700);
  
  c1->SetLeftMargin(0.14);
  g.Draw("AP");
  c1->SaveAs("Dphi.pdf");

  cout << endl;
  
}

void qT::dsdqT(const SampleFormat& summary){
  double nrm = summary.mc()->xsection()/(static_cast<float>(summary.nevents())*_histqT->GetBinWidth(1));
  TH1F hqT = *_histqT;
  TH1F hqTS = *_histqTSJA;
  TH1F *hqTSJA =  new TH1F("S1_PT_0","S1_PT_0",50,5.0,100.0);
  qTSJA(hqTSJA);//read from standard mode results
  
  hqT.GetYaxis()->SetTitle("#frac{d#sigma}{dq_{T}}[pb#bulletGeV^{-1}]");
  hqT.SetLineColor(kRed);
  hqT.Scale(nrm); hqTS.Scale(nrm);

  cout << endl;
  cout << _histqT << " vs " << &hqT << endl;
  
  TCanvas* c = new TCanvas("c","q_{T}", 500, 700);
  c->SetLeftMargin(0.14);
  //_histDphi->SetFillColor(kRed);
  hqT.Draw("HIST");
  hqTSJA->Draw("SAME");
  c->SetGrid();
  c->SetLogy(1);
  c->SaveAs("qT.pdf");

  for(int i=1; i<=hqT.GetNbinsX(); i++){//constrained to [1, GetNbinsX()]
    cout << hqT.GetBinCenter(i) << " " << hqT.GetBinContent(i) << " " << 1e-4*hqTSJA->GetBinContent(i)/hqTSJA->GetBinWidth(i) << " " << hqTS.GetBinContent(i) << endl;
  }

  cout << endl;
}

void qT::dsdphi(const SampleFormat& summary){
  double nrm = summary.mc()->xsection()/(static_cast<float>(summary.nevents())*_histdphi->GetBinWidth(1));
  TH1F hqT = *_histdphi;
  hqT.GetYaxis()->SetTitle("#frac{d#sigma}{d#delta#phi}[pb]");
  hqT.SetLineColor(kRed);
  hqT.Scale(nrm);
  for(int i=1; i<=hqT.GetNbinsX(); i++){//constrained to [1, GetNbinsX()]
    cout << hqT.GetBinCenter(i) << " " << hqT.GetBinContent(i) << endl;
  }  
  cout << endl;
  cout << "#{log(dphi), dsig/ddphi}" << endl;
  cout << "# of events: " << summary.nevents() << endl;
  for(int i=1; i<=_histldphi->GetNbinsX(); i++){//constrained to [1, GetNbinsX()]
    nrm = summary.mc()->xsection()/(static_cast<float>(summary.nevents())*(pow(10.0, _histldphi->GetBinLowEdge(i)+_histldphi->GetBinWidth(i)) - pow(10.0, _histldphi->GetBinLowEdge(i)) ));
    cout << _histldphi->GetBinCenter(i) << " " << nrm*_histldphi->GetBinContent(i)  << " " << nrm*_histldphiS->GetBinContent(i) << endl;
  }

  
  TCanvas* c = new TCanvas("c","dphi_{T}", 500, 700);
  c->SetLeftMargin(0.14);
  //_histDphi->SetFillColor(kRed);
  hqT.Draw("HIST");
  c->SetGrid();
  c->SetLogy(1);
  c->SaveAs("dsigmaodphi.pdf");
}

void qT::sigma(const SampleFormat& summary, TH1F* hist, const char* fname, const char* xlabel, const char* ylabel){
  double nrm = summary.mc()->xsection()/(static_cast<float>(summary.nevents()));
  cout << "Total # of events: " << summary.nevents() << " vs _numEvents: " << _numEvents << endl;
  
  //check the total number
  int numBin = hist->GetNbinsX()+1;
  TGraph g;

  double sig = nrm*hist->GetBinContent(numBin);
  int chk_num = hist->GetBinContent(numBin);
  for(int i=hist->GetNbinsX(); i>0; i--){//constrained to [1, GetNbinsX()]
    sig += nrm*hist->GetBinContent(i); chk_num += hist->GetBinContent(i);
    g.AddPoint(hist->GetBinLowEdge(i), sig);
    cout << hist->GetBinLowEdge(i) << " " << sig << endl;
  }
  cout << "Recounted # = " << chk_num + hist->GetBinContent(0) << " vs " << _numSelected << endl;
  TCanvas* c = new TCanvas("","#sigma", 500, 700);
  c->SetLeftMargin(0.14);
  //_histDphi->SetFillColor(kRed);
  g.Draw();
  g.GetYaxis()->SetTitle(ylabel);
  g.GetXaxis()->SetTitle(xlabel);
  //c->SetLogy(1);
  c->SaveAs(fname);
}


void qT::sigtot(const SampleFormat& summary){
  double nrm = summary.mc()->xsection()/(static_cast<float>(summary.nevents()));
  _sigWTA *= nrm; _sigSJA *= nrm;
  cout << "/nTotal cross section WTA = " << _sigWTA << " pb, SJA = " << _sigSJA << " pb" << ", WTA - SJA = " << _sigWTA - _sigSJA << " pb" << endl;
}

void qT::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;

  //compared with the normal mode
  cmpWithNorm(summary);

  
  //output dsigma/dq_T
  dsdqT(summary);
  dsdphi(summary);

  //output sigma(q_T)
  //cout << "#{log10(dphi), sigma[pb]}" << endl;
  //sigma(summary, _histldphi, "sigmadphi.pdf", "log(#delta#phi)", "#sigma(#delta#phi)[pb]");

  //cout << "#{dphi, sigma[pb]}" << endl;
  //sigma(summary, _histdphi, "sigmadphi.pdf", "#delta#phi", "#sigma(#delta#phi)[pb]");

  //cout << "#{qT, sigma[pb]}" << endl;
  //sigma(summary, _histqT, "sigmaqT.pdf", "q_{T}[GeV]", "#sigma(q_{T})[pb]");
  

  cout << "PDFs: " << _pdfIDA << ", " << _pdfIDB << endl;
  cout << "alpha_s = " << _aS << ", factorization scale = " << _scale << " GeV" << endl;
  //printSummary(summary);

  sigtot(summary);
  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool qT::Execute(SampleFormat& sample, const EventFormat& event)
{//for one event
  _numEvents++;
  if(_aS!=event.mc()->alphaQCD()) _aS = event.mc()->alphaQCD();
  if(_scale!=event.mc()->scale()) _scale = event.mc()->scale();
  if(_pdfIDA!=sample.mc()->beamPDFID().first) _pdfIDA=sample.mc()->beamPDFID().first;
  if(_pdfIDB!=sample.mc()->beamPDFID().second) _pdfIDB=sample.mc()->beamPDFID().second;
  
  //cout << sample.mc()->beamPDFID().second << endl;
  const MCParticleFormat *Js[3]; int iJ = 0;
  int idxI = 0;
  for (MAuint32 i=0;i<event.mc()->particles().size();i++){
    const MCParticleFormat* part = &event.mc()->particles()[i];
    if(part->pdgid()==9) cout << "Warning: pdg id = 9 for gluons!" << endl;
    if (PHYSICS->Id->IsInitialState(*part)){
      _initIDs[idxI] = part->pdgid(); idxI++;
    }else if (PHYSICS->Id->IsFinalState(*part)){
      Js[iJ] = part; _finalIDs[iJ] = part->pdgid(); 
      iJ++;
    }
  }

  if(chanelgg2ggg()){
  if(idxI!=2) cout << "There are " << idxI << " initial-state particles?" << endl;
  if(iJ!=3){
    cout << "There are " << iJ << " final-state particles?" << endl;
  }else{
    //Sort the pT
    int idx[3];
    SortpT(Js, idx);

    double pJ[3], etaJ[3], pJSJA[3], etaJSJA[3], qT, qTSJA=0.0;
    unsigned int nJ;
    bool injet = InJetQ(Js[idx[1]], Js[idx[2]]);
    if(injet){
      nJ = 2;
      //WTA
      pJ[0] = Js[idx[1]]->pt() + Js[idx[2]]->pt(); etaJ[0] = Js[idx[1]]->eta();
      pJ[1] = Js[idx[0]]->pt(); etaJ[1] = Js[idx[0]]->eta();
      qT = qTInJet(Js[idx[0]]->pt(), Js[idx[1]]->pt(), Js[idx[2]]->pt());
      
      //SJA
      pJSJA[0] = Js[idx[0]]->pt(); etaJSJA[0] = Js[idx[0]]->eta();
      pJSJA[1] =  pJSJA[0]; etaJSJA[1] = rapidity(Js[idx[1]], Js[idx[2]]);      
      qTSJA=0.0;
      
    }else{
      nJ = 3;

      for(unsigned int i=0; i<nJ; i++){
	pJ[i] = Js[idx[i]]->pt(); etaJ[i] = Js[idx[i]]->eta();
	pJSJA[i] = pJ[i]; etaJSJA[i] = etaJ[i];
      }
      /*
      pJ = Js[idx[0]]->pt(); etaJ = Js[idx[0]]->eta(); 
      pJ2 = Js[idx[1]]->pt(); etaJ2 = Js[idx[1]]->eta();
      */
      
      qTSJA = Js[idx[2]]->pt(); qT = PTVecSum(Js[idx[0]], Js[idx[1]]);
    }
    
    if(setCut(pJ, etaJ, nJ)){
      //if(selectpTQ(pJ2)&&selectetaQ(etaJ)&&selectetaQ(etaJ2)){
      double Dphi = DeltaPhi(Js[idx[0]], Js[idx[1]]);
      _histldphi->Fill(log10(TMath::Pi()-Dphi));
      _histdphi->Fill(TMath::Pi()-Dphi);
      _histDphi->Fill(Dphi);
      _histqT->Fill(qT);
      /*
      //SJA: only if j2, j3 are not in the same jet, qT>0
      if(!injet){
	_histldphiS->Fill(log10(TMath::Pi()-Dphi));
	_histdphiS->Fill(TMath::Pi()-Dphi);
	_histDphiS->Fill(Dphi);
      }
      _histqTSJA->Fill(qTSJA);
      */
      _sigWTA += 1.0;
    }
   
    if(setCut(pJSJA, etaJSJA, nJ)){
      //if(selectpTQ(pJ2)&&selectetaQ(etaJ)&&selectetaQ(etaJ2)){
      double Dphi = DeltaPhi(Js[idx[0]], Js[idx[1]]);
      /*
      _histldphi->Fill(log10(TMath::Pi()-Dphi));
      _histdphi->Fill(TMath::Pi()-Dphi);
      _histDphi->Fill(Dphi);
      _histqT->Fill(qT);
      */
      //SJA: only if j2, j3 are not in the same jet, qT>0
      if(!injet){
	_histldphiS->Fill(log10(TMath::Pi()-Dphi));
	_histdphiS->Fill(TMath::Pi()-Dphi);
	_histDphiS->Fill(Dphi);
      }
      _histqTSJA->Fill(qTSJA);

      _sigSJA += 1.0;
    }
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

double qT::rapidity(const MCParticleFormat* p1, const MCParticleFormat* p2){
  double px = p1->px()+p2->px(), py = p1->py()+p2->py(), pz = p1->pz()+p2->pz();
  double p = sqrt(px*px + py*py + pz*pz);
  return 0.5*log((p+pz)/(p-pz));
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

bool qT::selectpTQ(double pJ){
  bool Q = false;
  
  if(pJ > _pTJMin) Q = true;

  return Q;
}

bool qT::selectetaQ(double etaJ){
  bool Q = false;
  
  if((etaJ > -_etaMax)&&(etaJ < _etaMax)) Q = true;

  return Q;
}

void qT::setParas(double R, double etaMax, double pTJMin){
  _R = R; _etaMax = etaMax; _pTJMin = pTJMin;
}

void qT::setInitAll(){
  _initPartonsQ = false;
}

void qT::setInitPartons(int iA, int iB){
  _iA = iA; _iB = iB;
  _initPartonsQ = true;
}

bool qT::chanelud(){
  return ((_initIDs[0]==_u&&_initIDs[1]==_d)||(_initIDs[1]==_u&&_initIDs[0]==_d));
}

bool qT::chanelgg2ggg(){
  return (_initIDs[0]==_g&&_initIDs[1]==_g&&_finalIDs[0]==_g&&_finalIDs[1]==_g&&_finalIDs[2]==_g);
}

bool qT::setCut(double *pJ, double *etaJ, unsigned int nJ){
  bool cut = true;
  for(unsigned int i=0; i<nJ; i++){
    cut = cut&&selectpTQ(pJ[i])&&selectetaQ(etaJ[i]);
  }
  
  return cut;
}
