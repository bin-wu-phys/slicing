#include "SampleAnalyzer/User/Analyzer/qT.h"
#include "TCanvas.h"
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


  // Initializing histograms
  _histDphi = new TH1F("_histDphi", "#Delta#phi", 50, 20.0, 200.0);
  _histDphi->GetXaxis()->SetTitle("#Delta#phi");
  
  _histqT = new TH1F("_qT", "q_{T}", 50, 20.0, 200.0);
  _histqT->GetXaxis()->SetTitle("q_{T}[GeV]");

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

    // Creating a new TH1F
  TH1F* S1_PT_0 = new TH1F("S1_PT_0","S1_PT_0",50,20.0,200.0);
  // Content
  S1_PT_0->SetBinContent(0,0.0); // underflow
  S1_PT_0->SetBinContent(1,5059129226.522019);
  S1_PT_0->SetBinContent(2,18900140846.251972);
  S1_PT_0->SetBinContent(3,28445671273.65217);
  S1_PT_0->SetBinContent(4,41809401872.01191);
  S1_PT_0->SetBinContent(5,40377571807.901855);
  S1_PT_0->SetBinContent(6,49445822213.931885);
  S1_PT_0->SetBinContent(7,83618803744.02382);
  S1_PT_0->SetBinContent(8,97841644380.85013);
  S1_PT_0->SetBinContent(9,91159774081.67003);
  S1_PT_0->SetBinContent(10,79227863547.41989);
  S1_PT_0->SetBinContent(11,61377732748.18202);
  S1_PT_0->SetBinContent(12,56318602521.65996);
  S1_PT_0->SetBinContent(13,49159462201.110146);
  S1_PT_0->SetBinContent(14,37323011671.134224);
  S1_PT_0->SetBinContent(15,29495671320.665817);
  S1_PT_0->SetBinContent(16,26822931200.99414);
  S1_PT_0->SetBinContent(17,23863811068.499825);
  S1_PT_0->SetBinContent(18,19568330876.170116);
  S1_PT_0->SetBinContent(19,16513760739.402033);
  S1_PT_0->SetBinContent(20,14891020666.744005);
  S1_PT_0->SetBinContent(21,13650100611.181929);
  S1_PT_0->SetBinContent(22,8400063376.112022);
  S1_PT_0->SetBinContent(23,7731876346.194013);
  S1_PT_0->SetBinContent(24,7063689316.276003);
  S1_PT_0->SetBinContent(25,4963673222.247987);
  S1_PT_0->SetBinContent(26,4772763213.700009);
  S1_PT_0->SetBinContent(27,5154584230.796008);
  S1_PT_0->SetBinContent(28,4772763213.700009);
  S1_PT_0->SetBinContent(29,3245479145.3160133);
  S1_PT_0->SetBinContent(30,3245479145.3160133);
  S1_PT_0->SetBinContent(31,2768202123.9459815);
  S1_PT_0->SetBinContent(32,1145463051.2879968);
  S1_PT_0->SetBinContent(33,2195471098.3020053);
  S1_PT_0->SetBinContent(34,2100016094.0280166);
  S1_PT_0->SetBinContent(35,2195471098.3020053);
  S1_PT_0->SetBinContent(36,859097338.466);
  S1_PT_0->SetBinContent(37,1240918055.5619855);
  S1_PT_0->SetBinContent(38,763642034.1919979);
  S1_PT_0->SetBinContent(39,954552642.7400019);
  S1_PT_0->SetBinContent(40,1145463051.2879968);
  S1_PT_0->SetBinContent(41,572731525.6439984);
  S1_PT_0->SetBinContent(42,190910508.54799947);
  S1_PT_0->SetBinContent(43,668186829.9180005);
  S1_PT_0->SetBinContent(44,954552642.7400019);
  S1_PT_0->SetBinContent(45,668186829.9180005);
  S1_PT_0->SetBinContent(46,477276321.37000096);
  S1_PT_0->SetBinContent(47,381821017.09599894);
  S1_PT_0->SetBinContent(48,668186829.9180005);
  S1_PT_0->SetBinContent(49,286365812.82200146);
  S1_PT_0->SetBinContent(50,95455264.27400018);
  S1_PT_0->SetBinContent(51,0.0); // overflow
  S1_PT_0->SetEntries(10000);
  
  //Output
  TCanvas* c1 = new TCanvas("c1","#Delta#phi", 500, 700);
  _histDphi->Draw("HIST");
  _histDphi->SetFillColor(kRed);
  S1_PT_0->Draw("SAME");
  S1_PT_0->SetFillColorAlpha(kBlue, 0.5);;
  c1->SaveAs("Dphi.pdf");

  TCanvas* c2 = new TCanvas("c2","q_{T}", 500, 700);
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
    _histDphi->Fill(Js[idx[0]]->pt());
    _histqT->Fill(Js[idx[2]]->pt());
  }
  
  return true;
}

void qT::SortpT(const MCParticleFormat** Js){
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
    const MCParticleFormat *j1 = Js[iMax], *j2, *j3 = Js[iMin];
    for(int i=0;i<3;i++){
      if((i!=iMax)&&(i!=iMin)) j2 = Js[i];
    }
    Js[0] = j1; Js[1] = j2; Js[2] = j3;
  }
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
      if((i!=iMax)&&(i!=iMin)) idx[2] = i;
    }
  }else{
    for(int i=0;i<3;i++)
      idx[i] = i;
  }
}
