//R__LOAD_LIBRARY(libMcIniData.so)
R__LOAD_LIBRARY(libUniGenFormat.so)

using namespace ROOT;
using namespace ROOT::RDF;

// Kinematic cuts for the multiplicity
const double ETA_CUT = 0.5;      // |eta| < ETA_CUT
const double PT_MIN = 0.2;       // pT > PT_MIN, GeV/c
const bool CHARGED_ONLY = true;  // count only charged particles

// Histogram binning
const int N_B_BINS = 200;
const double B_MIN = 0.;
const double B_MAX = 20.;        // fm
const int N_MULT_BINS = 500;
const double MULT_MIN = 0.;
const double MULT_MAX = 500.;

TChain* makeChain(string& filename, const char* treename) {
  cout << "Adding files to chain:" << endl;
  TChain *chain = new TChain(treename);
  if (filename.rfind(".root") < filename.size())
    chain->Add(filename.data());
  else {
    TFileCollection fc("fc", "", filename.c_str());
    chain->AddFileInfoList((TCollection*)fc.GetList());
  }
  chain->ls();
  return chain;
}

int getCharge(int pdg)
{
  if (abs(pdg) > 1000000000) // nuclei: 10LZZZAAAI
    return pdg/10000%1000;
  auto particle = TDatabasePDG::Instance()->GetParticle(pdg);
  if (!particle)
    return 0;
  return (int) (particle->Charge()/3.); // TParticlePDG charge is in units of |e|/3
}

int getRefMult(const UEvent &event)
{
  int mult = 0;
  for (int i=0; i<event.GetNpa(); i++) {
    auto particle = event.GetParticle(i);
    if (CHARGED_ONLY && getCharge(particle->GetPdg()) == 0) continue;
    double px = particle->Px();
    double py = particle->Py();
    double pz = particle->Pz();
    double pt = sqrt(px*px + py*py);
    if (pt <= PT_MIN) continue;
    double eta = asinh(pz/pt);
    if (abs(eta) >= ETA_CUT) continue;
    mult++;
  }
  return mult;
}

void run_refmult_mcini(string fileIn="", string fileOut="")
{
  TStopwatch timer;
  timer.Start();

  std::cout << "Multiplicity cuts: |eta| < " << ETA_CUT << ", pT > " << PT_MIN << " GeV/c"
            << (CHARGED_ONLY ? ", charged only" : "") << std::endl;

  TChain *chain = makeChain(fileIn, "events");
  ROOT::RDataFrame d(*chain);

  auto dd = d
    .Define("b", [](const UEvent &event){ return event.GetB(); }, {"event"})
    .Define("refMult", getRefMult, {"event"})
  ;

  auto h2_refMult_b = dd.Histo2D({"h2_refMult_b",
                                  Form("Multiplicity (|#eta|<%.1f, p_{T}>%.1f GeV/c) vs b;b, fm;N_{ch}", ETA_CUT, PT_MIN),
                                  N_B_BINS, B_MIN, B_MAX, N_MULT_BINS, MULT_MIN, MULT_MAX},
                                 "b", "refMult");

  TFile fOut(fileOut.c_str(), "recreate");
  fOut.cd();
  h2_refMult_b->Write();
  fOut.Close();

  timer.Stop();
  timer.Print();
}
