#include "FunctionsQa.C"

void runQaMpd_step1(string fileIn="", string fileOut="", std::string cm_energy="2.5", std::string str_nucleus_mass="209")
{
  TStopwatch timer;
  timer.Start();

  const double sNN = std::stod( cm_energy ); // in GeV
  const double M = 0.938; // in GeV/c^2
  const double T = sNN*sNN/(2.*M) - 2.*M;
  const double E = T + M;
  const double P = sqrt( E*E - M*M );
  const double Y_BEAM = 0.25 * log( (E + P) / (E - P) );
  const double nucleus_mass = std::stod(str_nucleus_mass);
  const double NUCLEUS_RADIUS = 1.25 * pow( nucleus_mass, 1.0 / 3.0 );

  std::cout << "sqrtSnn = " << sNN << " GeV; T = " << T << "A GeV; Y_BEAM = " << Y_BEAM << std::endl;
  std::cout << "A = " << nucleus_mass << "; R = " << NUCLEUS_RADIUS << std::endl;

  ROOT::RDataFrame d("t", fileIn.c_str());
  TFile fOut(fileOut.c_str(),"recreate");

  vector <RResultPtr<::TH1D >> hists;
  vector <RResultPtr<::TH2D >> hists2d;
  vector <RResultPtr<::TH3D >> hists3d;
  vector <RResultPtr<::THnD >> histsnd;
  vector <RResultPtr<::TProfile >> profs;
  vector <RResultPtr<::TProfile2D >> profs2d;

  auto dd = d
    .Filter("mcB<20.")
    .Define("simB", "mcB")
    .Define("refMult", getRefMultFxt, {"recoGlobalMom", "recoGlobalNhits"})
    .Define("recDca", getDcaMag, {"recoGlobalDca"})
    .Define("isTOF", "recoGlobalTofFlag")
    .Define("isGoodTrack", isGoodTrack, {"recoGlobalMom", "recoGlobalNhits", "recDca", "isTOF"})
    .Define("isPrimary", isRecPrim, {"recDca"})
    .Define("recTpcDedx", "recoGlobalTpcDedx")
    .Define("recTofMass2", "recoGlobalTofMass2")
    .Define("isProton", "recoGlobalSimPdg==2212")
    .Define("isPionP", "recoGlobalSimPdg==211")
    .Define("isPionM", "recoGlobalSimPdg==-211")
    .Define("isPions", "abs(recoGlobalSimPdg)==211")
    .Define("isKaonP", "recoGlobalSimPdg==321")
    .Define("isKaonM", "recoGlobalSimPdg==-321")
    .Define("isKaons", "abs(recoGlobalSimPdg)==321")
    .Define("recPt", getPt, {"recoGlobalMom"})
    .Define("recP", getP, {"recoGlobalMom"})
    .Define("recRigidity", getRigidity, {"recoGlobalMom", "recoGlobalCharge"})
    .Define("recRigidityReverse", getRigidityReverse, {"recoGlobalMom", "recoGlobalCharge"})
    .Define("recEta", getEta, {"recoGlobalMom"})
    .Define("recPhi", getPhi, {"recoGlobalMom"})
    .Define("recY", [Y_BEAM]( const RVec<int> vec_pdg, vector<fourVector> vec_momentum ){
      RVec<float> vec_y;
      vec_y.reserve( vec_pdg.size() );
      for( int i=0; i<vec_pdg.size(); ++i ){
        auto pdg = vec_pdg.at(i);
        if( pdg == 0 ){
          vec_y.push_back( -999. );
          continue;
        }
        auto pz = vec_momentum.at(i).Pz();
        auto p = vec_momentum.at(i).P();
        TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdg);
        if (!particle) {
          vec_y.push_back( -999. );
          continue;
        }
        auto m = particle->Mass();
        auto E = sqrt( p*p + m*m );
        auto y = 0.5 * log( (E+pz)/(E-pz) ) - Y_BEAM;
        vec_y.push_back(y);
      }
      return vec_y;
    },{"recoGlobalSimPdg", "recoGlobalMom"})
    .Define("isSimPrimary", "simMotherId==-1")
    .Define("isSimProton", "simPdg==2212")
    .Define("isSimPionP", "simPdg==211")
    .Define("isSimPionM", "simPdg==-211")
    .Define("isSimPions", "abs(simPdg)==211")
    .Define("isSimKaonP", "simPdg==321")
    .Define("isSimKaonM", "simPdg==-321")
    .Define("isSimKaons", "abs(simPdg)==321")
    .Define("simPt", getPt, {"simMom"})
    .Define("simEta", getEta, {"simMom"})
    .Define("simPhi", getPhi, {"simMom"})
    .Define("simY", [Y_BEAM]( const RVec<int> vec_pdg, vector<fourVector> vec_momentum ){
      RVec<float> vec_y;
      vec_y.reserve( vec_pdg.size() );
      for( int i=0; i<vec_pdg.size(); ++i ){
        auto pdg = vec_pdg.at(i);
        if( pdg == 0 ){
          vec_y.push_back( -999. );
          continue;
        }
        auto pz = vec_momentum.at(i).Pz();
        auto p = vec_momentum.at(i).P();
        TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdg);
        if (!particle) {
          vec_y.push_back( -999. );
          continue;
        }
        auto m = particle->Mass();
        auto E = sqrt( p*p + m*m );
        auto y = 0.5 * log( (E+pz)/(E-pz) ) - Y_BEAM;
        vec_y.push_back(y);
      }
      return vec_y;
    },{"simPdg", "simMom"})
  ;

  dd.Foreach([](ULong64_t evtId){if (evtId % 100 == 0) cout << "\r" << evtId;}, {"rdfentry_"});

  // Make lists of histograms for QA
  hists.push_back(dd.Histo1D({"h1_recVtx_X","Reconstructed vertex X;x (cm)",500,-1,1}, "recoPrimVtxX")); 
  hists.push_back(dd.Histo1D({"h1_recVtx_Y","Reconstructed vertex Y;y (cm)",500,-1,1}, "recoPrimVtxY")); 
  hists.push_back(dd.Histo1D({"h1_recVtx_Z","Reconstructed vertex Z;z (cm)",500,-1,1}, "recoPrimVtxZ"));
  hists.push_back(dd.Histo1D({"h1_refMult", "Reconstructed N_{ch};N_{ch}",1000,0.,1000.}, "refMult"));
  hists.push_back(dd.Histo1D({"h1_simB", "Simulated b;b, fm",200,0.,20.}, "simB"));
  hists2d.push_back(dd.Histo2D({"h2_recVtx_XY","Reconstructed vertex XY;x (cm);y (cm)",500,-1,1,500,-1,1}, "recoPrimVtxX", "recoPrimVtxY"));
  hists2d.push_back(dd.Histo2D({"h2_simBrefMult", "b vs N_{ch};N_{ch};b, fm",1000,0.,1000.,200,0.,20.}, "refMult", "simB"));
  hists2d.push_back(dd.Histo2D({"h2_recDedxPq", "Reconstructed dEdx vs P/q;p/q, GeV/c;dEdx, a.u.", 1000, -5., 5., 500, 0., 5.e3}, "recRigidity", "recTpcDedx", "isGoodTrack"));
  hists2d.push_back(dd.Histo2D({"h2_recDedxPqReverse", "Reconstructed dEdx vs -P/q;-p/q, GeV/c;dEdx, a.u.", 1000, -5., 5., 500, 0., 5.e3}, "recRigidityReverse", "recTpcDedx", "isGoodTrack"));
  hists2d.push_back(dd.Histo2D({"h2_recDedxPt", "Reconstructed dEdx vs P_{T};p_{T}, GeV/c;dEdx, a.u.", 1000, -5., 5., 500, 0., 5.e3}, "recPt", "recTpcDedx", "isGoodTrack"));
  hists2d.push_back(dd.Histo2D({"h2_recM2Pq", "Reconstructed m^{2}_{TOF} vs P/q;p/q, GeV/c;m^{2}, (GeV/c^{2})^{2}", 1000, -5., 5., 220, -0.2, 2.}, "recRigidity", "recTofMass2", "isGoodTrack"));
  hists2d.push_back(dd.Histo2D({"h2_recM2Pt", "Reconstructed m^{2}_{TOF} vs p_{T};p_{T}, GeV/c;m^{2}, (GeV/c^{2})^{2}", 500, 0., 5., 220, -0.2, 2.}, "recPt", "recTofMass2", "isGoodTrack"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt", "Reconstructed #eta vs p_{T}", 1000, -5., 5., 500, 0., 5.}, "recEta", "recPt"));
  hists3d.push_back(dd.Histo3D({"h3_recDedxM2Pq", "Reconstructed dEdx vs m^{2}_{TOF} vs P/q;m^{2}, (GeV/c^{2})^{2};dEdx, a.u.;p/q, GeV/c", 220, -0.2, 2., 500, 0., 5.e3, 500, -5., 5.}, "recTofMass2", "recTpcDedx", "recRigidity", "isGoodTrack"));
  hists3d.push_back(dd.Histo3D({"h3_recDedxM2PqReverse", "Reconstructed dEdx vs m^{2}_{TOF} vs -P/q;m^{2}, (GeV/c^{2})^{2};dEdx, a.u.;-p/q, GeV/c", 220, -0.2, 2., 500, 0., 5.e3, 500, -5., 5.}, "recTofMass2", "recTpcDedx", "recRigidityReverse", "isGoodTrack"));
  // Write QA histograms to the output file
  fOut.cd();
  for (auto& hist:hists)
    hist->Write();
  for (auto& hist:hists2d)
    hist->Write();
  for (auto& hist:hists3d)
    hist->Write();
  for (auto& hist:profs)
    hist->Write();
  for (auto& hist:profs2d)
    hist->Write();
  fOut.Close();

  std::cout << std::endl;
  timer.Stop();
  timer.Print();
}
