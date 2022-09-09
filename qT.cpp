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
  _histDphi = new TH1F("#Delta#phi", "#Delta#phi", 50, 0.0, 0.5);
  _histDphi->GetXaxis()->SetTitle("#Delta#phi");
  _histDphi->GetYaxis()->SetTitle("Events  ( L_{int} = 10 fb^{-1} )");
  
  _histqT = new TH1F("q_{T}", "q_{T}", 50, 0.0, 200.0);
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

    TH1F* S1_DELTAR_0 = new TH1F("S1_DELTAR_0","S1_DELTAR_0",50,0.0,0.5);
  // Content
  S1_DELTAR_0->SetBinContent(0,0.0); // underflow
  S1_DELTAR_0->SetBinContent(1,285697597605.6);
  S1_DELTAR_0->SetBinContent(2,34459349711.2);
  S1_DELTAR_0->SetBinContent(3,20713789826.400013);
  S1_DELTAR_0->SetBinContent(4,16036479865.600033);
  S1_DELTAR_0->SetBinContent(5,11931909899.99998);
  S1_DELTAR_0->SetBinContent(6,8590972928.000004);
  S1_DELTAR_0->SetBinContent(7,6300046947.200002);
  S1_DELTAR_0->SetBinContent(8,6681867944.000003);
  S1_DELTAR_0->SetBinContent(9,5631859952.800003);
  S1_DELTAR_0->SetBinContent(10,5250038956.000003);
  S1_DELTAR_0->SetBinContent(11,3722754968.8000016);
  S1_DELTAR_0->SetBinContent(12,4677307960.799998);
  S1_DELTAR_0->SetBinContent(13,4581851961.600004);
  S1_DELTAR_0->SetBinContent(14,3627299969.5999994);
  S1_DELTAR_0->SetBinContent(15,4390941963.2);
  S1_DELTAR_0->SetBinContent(16,2672746977.600003);
  S1_DELTAR_0->SetBinContent(17,3913664967.200006);
  S1_DELTAR_0->SetBinContent(18,2672746977.600003);
  S1_DELTAR_0->SetBinContent(19,2481836979.1999984);
  S1_DELTAR_0->SetBinContent(20,2100015982.399998);
  S1_DELTAR_0->SetBinContent(21,3245478972.799999);
  S1_DELTAR_0->SetBinContent(22,3054567974.400003);
  S1_DELTAR_0->SetBinContent(23,1909104984.000002);
  S1_DELTAR_0->SetBinContent(24,3054567974.400003);
  S1_DELTAR_0->SetBinContent(25,2195470981.6);
  S1_DELTAR_0->SetBinContent(26,1813649984.7999997);
  S1_DELTAR_0->SetBinContent(27,2004559983.200004);
  S1_DELTAR_0->SetBinContent(28,1431828987.9999993);
  S1_DELTAR_0->SetBinContent(29,2100015982.399998);
  S1_DELTAR_0->SetBinContent(30,1145462990.400001);
  S1_DELTAR_0->SetBinContent(31,2100015982.399998);
  S1_DELTAR_0->SetBinContent(32,1622738986.4000037);
  S1_DELTAR_0->SetBinContent(33,1145462990.400001);
  S1_DELTAR_0->SetBinContent(34,2004559983.200004);
  S1_DELTAR_0->SetBinContent(35,1527283987.2000015);
  S1_DELTAR_0->SetBinContent(36,1622738986.4000037);
  S1_DELTAR_0->SetBinContent(37,1145462990.400001);
  S1_DELTAR_0->SetBinContent(38,477276296.00000006);
  S1_DELTAR_0->SetBinContent(39,2195470981.6);
  S1_DELTAR_0->SetBinContent(40,1240917989.6000032);
  S1_DELTAR_0->SetBinContent(41,1336373988.799997);
  S1_DELTAR_0->SetBinContent(42,1050007991.199999);
  S1_DELTAR_0->SetBinContent(43,1622738986.4000037);
  S1_DELTAR_0->SetBinContent(44,954552592.0000001);
  S1_DELTAR_0->SetBinContent(45,2004559983.200004);
  S1_DELTAR_0->SetBinContent(46,286365797.59999985);
  S1_DELTAR_0->SetBinContent(47,1145462990.400001);
  S1_DELTAR_0->SetBinContent(48,859097292.8000004);
  S1_DELTAR_0->SetBinContent(49,1431828987.9999993);
  S1_DELTAR_0->SetBinContent(50,1718194985.5999975);
  S1_DELTAR_0->SetBinContent(51,468971696069.6); // overflow
  S1_DELTAR_0->SetEntries(10000);

  c1->SetLeftMargin(0.14);
  _histDphi->SetLineColor(kBlue);
  _histDphi->Draw("HIST");
  S1_DELTAR_0->SetLineStyle(2);
  S1_DELTAR_0->SetLineColor(kRed);
  S1_DELTAR_0->Draw("SAME");
  c1->SetLogy(1);
  c1->SaveAs("DeltaR23.pdf");

  TCanvas* c2 = new TCanvas("c2","q_{T}", 500, 700);
  TH1F* S2_PT_0 = new TH1F("S2_PT_0","S2_PT_0",50,0.0,200.0);
  // Content
  S2_PT_0->SetBinContent(0,0.0); // underflow
  S2_PT_0->SetBinContent(1,0.0);
  S2_PT_0->SetBinContent(2,0.0);
  S2_PT_0->SetBinContent(3,0.0);
  S2_PT_0->SetBinContent(4,0.0);
  S2_PT_0->SetBinContent(5,0.0);
  S2_PT_0->SetBinContent(6,6300047261.624003);
  S2_PT_0->SetBinContent(7,23959270994.964027);
  S2_PT_0->SetBinContent(8,34554801434.96789);
  S2_PT_0->SetBinContent(9,46295801922.54003);
  S2_PT_0->SetBinContent(10,46009431910.64786);
  S2_PT_0->SetBinContent(11,79896053317.86803);
  S2_PT_0->SetBinContent(12,109009904526.8879);
  S2_PT_0->SetBinContent(13,101182604201.84117);
  S2_PT_0->SetBinContent(14,88677933682.55588);
  S2_PT_0->SetBinContent(15,68345962838.223854);
  S2_PT_0->SetBinContent(16,61664102560.74419);
  S2_PT_0->SetBinContent(17,48873092029.56795);
  S2_PT_0->SetBinContent(18,38182101585.5999);
  S2_PT_0->SetBinContent(19,30736591276.407894);
  S2_PT_0->SetBinContent(20,27872931157.4878);
  S2_PT_0->SetBinContent(21,24150181002.892006);
  S2_PT_0->SetBinContent(22,18900140784.87197);
  S2_PT_0->SetBinContent(23,16609210689.73581);
  S2_PT_0->SetBinContent(24,14604650606.491825);
  S2_PT_0->SetBinContent(25,8877339368.652006);
  S2_PT_0->SetBinContent(26,8590973356.759996);
  S2_PT_0->SetBinContent(27,7254599301.26398);
  S2_PT_0->SetBinContent(28,5250039218.019996);
  S2_PT_0->SetBinContent(29,5727315237.839984);
  S2_PT_0->SetBinContent(30,5154584214.056006);
  S2_PT_0->SetBinContent(31,4104576170.451999);
  S2_PT_0->SetBinContent(32,3436389142.7039905);
  S2_PT_0->SetBinContent(33,2768202114.955982);
  S2_PT_0->SetBinContent(34,1718195071.352016);
  S2_PT_0->SetBinContent(35,2386381099.099983);
  S2_PT_0->SetBinContent(36,2386381099.099983);
  S2_PT_0->SetBinContent(37,1718195071.352016);
  S2_PT_0->SetBinContent(38,1050008043.6040075);
  S2_PT_0->SetBinContent(39,1145463047.567997);
  S2_PT_0->SetBinContent(40,1050008043.6040075);
  S2_PT_0->SetBinContent(41,1145463047.567997);
  S2_PT_0->SetBinContent(42,572731523.7839985);
  S2_PT_0->SetBinContent(43,190910507.92799947);
  S2_PT_0->SetBinContent(44,954552639.6400015);
  S2_PT_0->SetBinContent(45,859097335.6759998);
  S2_PT_0->SetBinContent(46,859097335.6759998);
  S2_PT_0->SetBinContent(47,477276319.82000077);
  S2_PT_0->SetBinContent(48,668186827.7480003);
  S2_PT_0->SetBinContent(49,286365811.8920013);
  S2_PT_0->SetBinContent(50,95455263.96400015);
  S2_PT_0->SetBinContent(51,0.0); // overflow
  S2_PT_0->SetEntries(10000);
  
  c2->SetLeftMargin(0.14);
  _histqT->SetLineColor(kBlue);
  _histqT->Draw("HIST");
  S2_PT_0->SetLineStyle(2);
  S2_PT_0->SetLineColor(kRed);
  S2_PT_0->Draw("SAME");
  c2->SaveAs("pT23.pdf");  
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
    _histDphi->Fill(DeltaR(Js[idx[1]], Js[idx[2]]));
    _histqT->Fill(PTVecSum(Js[idx[1]], Js[idx[2]]));
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
