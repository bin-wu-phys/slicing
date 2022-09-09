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
  
  _histp2 = new TH1F("p2", "p2", 50, 20.0, 200.0);

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
  _histDphi->Scale(nrm); _histp2->Scale(nrm); _histqT->Scale(nrm);

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

  TH1F* S2_PT_0 = new TH1F("S2_PT_0","S2_PT_0",50,20.0,200.0);
  // Content
  S2_PT_0->SetBinContent(0,0.0); // underflow
  S2_PT_0->SetBinContent(1,95455262416.00006);
  S2_PT_0->SetBinContent(2,157692103991.23236);
  S2_PT_0->SetBinContent(3,144232903650.57614);
  S2_PT_0->SetBinContent(4,114737202904.03151);
  S2_PT_0->SetBinContent(5,89918852275.87193);
  S2_PT_0->SetBinContent(6,69300521754.01608);
  S2_PT_0->SetBinContent(7,51259471297.39192);
  S2_PT_0->SetBinContent(8,42954871087.2001);
  S2_PT_0->SetBinContent(9,35413900896.33599);
  S2_PT_0->SetBinContent(10,27395660693.39203);
  S2_PT_0->SetBinContent(11,22909260579.839954);
  S2_PT_0->SetBinContent(12,19091050483.199963);
  S2_PT_0->SetBinContent(13,14795560374.479876);
  S2_PT_0->SetBinContent(14,11359180287.50411);
  S2_PT_0->SetBinContent(15,9927347251.264006);
  S2_PT_0->SetBinContent(16,6490957164.287988);
  S2_PT_0->SetBinContent(17,5536405140.128002);
  S2_PT_0->SetBinContent(18,5154584130.464003);
  S2_PT_0->SetBinContent(19,4009121101.472005);
  S2_PT_0->SetBinContent(20,4486397113.551997);
  S2_PT_0->SetBinContent(21,4295487108.720011);
  S2_PT_0->SetBinContent(22,1622739041.0719905);
  S2_PT_0->SetBinContent(23,1431829036.2400036);
  S2_PT_0->SetBinContent(24,2290926057.9839954);
  S2_PT_0->SetBinContent(25,1813650045.9040027);
  S2_PT_0->SetBinContent(26,1050008026.5760043);
  S2_PT_0->SetBinContent(27,1622739041.0719905);
  S2_PT_0->SetBinContent(28,859097321.7439996);
  S2_PT_0->SetBinContent(29,1240918031.4079912);
  S2_PT_0->SetBinContent(30,1145463028.9919977);
  S2_PT_0->SetBinContent(31,954552624.1600007);
  S2_PT_0->SetBinContent(32,763642019.3279985);
  S2_PT_0->SetBinContent(33,572731514.4959989);
  S2_PT_0->SetBinContent(34,477276312.08000034);
  S2_PT_0->SetBinContent(35,286365807.2480007);
  S2_PT_0->SetBinContent(36,381821009.66399926);
  S2_PT_0->SetBinContent(37,381821009.66399926);
  S2_PT_0->SetBinContent(38,381821009.66399926);
  S2_PT_0->SetBinContent(39,95455262.41600007);
  S2_PT_0->SetBinContent(40,95455262.41600007);
  S2_PT_0->SetBinContent(41,477276312.08000034);
  S2_PT_0->SetBinContent(42,0.0);
  S2_PT_0->SetBinContent(43,95455262.41600007);
  S2_PT_0->SetBinContent(44,95455262.41600007);
  S2_PT_0->SetBinContent(45,0.0);
  S2_PT_0->SetBinContent(46,0.0);
  S2_PT_0->SetBinContent(47,0.0);
  S2_PT_0->SetBinContent(48,0.0);
  S2_PT_0->SetBinContent(49,0.0);
  S2_PT_0->SetBinContent(50,0.0);
  S2_PT_0->SetBinContent(51,0.0); // overflow
  S2_PT_0->SetEntries(10000);


    // Creating a new TH1F
  TH1F* S3_PT_0 = new TH1F("S3_PT_0","S3_PT_0",50,20.0,200.0);
  // Content
  S3_PT_0->SetBinContent(0,0.0); // underflow
  S3_PT_0->SetBinContent(1,462862544746.5709);
  S3_PT_0->SetBinContent(2,219356221205.94727);
  S3_PT_0->SetBinContent(3,109773510612.19629);
  S3_PT_0->SetBinContent(4,58609535665.992584);
  S3_PT_0->SetBinContent(5,38468473718.88438);
  S3_PT_0->SetBinContent(6,22813812205.49249);
  S3_PT_0->SetBinContent(7,14031921356.51582);
  S3_PT_0->SetBinContent(8,9163705885.88809);
  S3_PT_0->SetBinContent(9,6968234673.644067);
  S3_PT_0->SetBinContent(10,3818210369.119997);
  S3_PT_0->SetBinContent(11,2577292249.156022);
  S3_PT_0->SetBinContent(12,1909105184.5599985);
  S3_PT_0->SetBinContent(13,1336374129.1920474);
  S3_PT_0->SetBinContent(14,954552692.2800089);
  S3_PT_0->SetBinContent(15,477276346.14000446);
  S3_PT_0->SetBinContent(16,381821036.9119997);
  S3_PT_0->SetBinContent(17,477276346.14000446);
  S3_PT_0->SetBinContent(18,190910518.45599985);
  S3_PT_0->SetBinContent(19,190910518.45599985);
  S3_PT_0->SetBinContent(20,95455269.2280009);
  S3_PT_0->SetBinContent(21,0.0);
  S3_PT_0->SetBinContent(22,0.0);
  S3_PT_0->SetBinContent(23,95455269.2280009);
  S3_PT_0->SetBinContent(24,0.0);
  S3_PT_0->SetBinContent(25,0.0);
  S3_PT_0->SetBinContent(26,0.0);
  S3_PT_0->SetBinContent(27,0.0);
  S3_PT_0->SetBinContent(28,0.0);
  S3_PT_0->SetBinContent(29,0.0);
  S3_PT_0->SetBinContent(30,0.0);
  S3_PT_0->SetBinContent(31,0.0);
  S3_PT_0->SetBinContent(32,0.0);
  S3_PT_0->SetBinContent(33,0.0);
  S3_PT_0->SetBinContent(34,0.0);
  S3_PT_0->SetBinContent(35,0.0);
  S3_PT_0->SetBinContent(36,0.0);
  S3_PT_0->SetBinContent(37,0.0);
  S3_PT_0->SetBinContent(38,0.0);
  S3_PT_0->SetBinContent(39,0.0);
  S3_PT_0->SetBinContent(40,0.0);
  S3_PT_0->SetBinContent(41,0.0);
  S3_PT_0->SetBinContent(42,0.0);
  S3_PT_0->SetBinContent(43,0.0);
  S3_PT_0->SetBinContent(44,0.0);
  S3_PT_0->SetBinContent(45,0.0);
  S3_PT_0->SetBinContent(46,0.0);
  S3_PT_0->SetBinContent(47,0.0);
  S3_PT_0->SetBinContent(48,0.0);
  S3_PT_0->SetBinContent(49,0.0);
  S3_PT_0->SetBinContent(50,0.0);
  S3_PT_0->SetBinContent(51,0.0); // overflow
  S3_PT_0->SetEntries(10000);

  
  //Output
  TCanvas* c1 = new TCanvas("c1","#Delta#phi", 500, 700);
  _histDphi->SetFillColor(kRed);
  S1_PT_0->SetFillColorAlpha(kBlue, 0.5);;

  _histDphi->Draw("HIST");
  S1_PT_0->Draw("SAME");
  c1->SaveAs("p1.pdf");

  TCanvas* c3 = new TCanvas("c3","p2", 500, 700);
  _histp2->SetFillColor(kRed);
  S2_PT_0->SetFillColorAlpha(kBlue, 0.5);;

  _histp2->Draw("HIST");
  S2_PT_0->Draw("SAME");
  c3->SaveAs("p2.pdf");  

  TCanvas* c2 = new TCanvas("c2","q_{T}", 500, 700);
  _histqT->SetFillColor(kRed);
  S3_PT_0->SetFillColorAlpha(kBlue, 0.5);;

  _histqT->Draw("HIST");
  S3_PT_0->Draw("SAME");
  c2->SaveAs("p3.pdf");  
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
    _histp2->Fill(Js[idx[1]]->pt());
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
      if((i!=iMax)&&(i!=iMin)) idx[1] = i;
    }
  }else{
    for(int i=0;i<3;i++)
      idx[i] = i;
  }
}
