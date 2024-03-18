#include "FunctionsQa.C"
void runQaMpd_step2(string fileIn="", string fileOut="", std::string cm_energy="2.5", std::string str_nucleus_mass="209", std::string fileStep1="")
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

  auto fiStep1 = new TFile(fileStep1.c_str(), "read");
  auto f1_dedx_pip = (TF1*)fiStep1->Get("f1_dedx_pip");
  auto f1_dedx_p = (TF1*)fiStep1->Get("f1_dedx_proton");
  auto f1_dedx_sigm_pip = (TF1*)fiStep1->Get("f1_dedx_sigm_pip");
  auto f1_dedx_sigm_p = (TF1*)fiStep1->Get("f1_dedx_sigm_p");
  auto f1_m2_sigm_pip = (TF1*)fiStep1->Get("f1_m2_sigm_pip");
  auto f1_m2_sigm_p = (TF1*)fiStep1->Get("f1_m2_sigm_p");

  if (!fiStep1){
    std::cerr << "ERROR: Cannot open step1 QA file: " << fileStep1.c_str() << std::endl;
    return;
  }

  if (!f1_dedx_pip){
    std::cerr << "ERROR: Cannot find f1_dedx_pip in the step1 QA file!" << std::endl;
    return;
  }
  if (!f1_dedx_p){
    std::cerr << "ERROR: Cannot find f1_dedx_proton in the step1 QA file!" << std::endl;
    return;
  }
  if (!f1_dedx_sigm_pip){
    std::cerr << "ERROR: Cannot find f1_dedx_sigm_pip in the step1 QA file!" << std::endl;
    return;
  }
  if (!f1_dedx_sigm_p){
    std::cerr << "ERROR: Cannot find f1_dedx_sigm_p in the step1 QA file!" << std::endl;
    return;
  }
  if (!f1_m2_sigm_pip){
    std::cerr << "ERROR: Cannot find f1_m2_sigm_pip in the step1 QA file!" << std::endl;
    return;
  }
  if (!f1_m2_sigm_p){
    std::cerr << "ERROR: Cannot find f1_m2_sigm_p in the step1 QA file!" << std::endl;
    return;
  }

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
    .Define("isGoodTrackProton", isGoodTrack4Protons, {"recoGlobalMom", "recoGlobalNhits", "recDca"})
    .Define("isGoodTrackPionP", isGoodTrack4Pions, {"recoGlobalMom", "recoGlobalNhits", "recDca"})
    .Define("isGoodTrackPionM", isGoodTrack4Pions, {"recoGlobalMom", "recoGlobalNhits", "recDca"})
    .Define("isGoodTrackPions", isGoodTrack4Pions, {"recoGlobalMom", "recoGlobalNhits", "recDca"})
    .Define("isGoodPosEta", isGoodPosEta, {"recoGlobalMom", "isGoodTrack"})
    .Define("isGoodNegEta", isGoodNegEta, {"recoGlobalMom", "isGoodTrack"})
    //.Define("isPrimary", "recoGlobalSimMotherId==-1")
    .Define("isPrimary", isRecPrim, {"recDca"})
    .Define("recTpcDedx", "recoGlobalTpcDedx")
    .Define("recTofMass2", "recoGlobalTofMass2")
    //.Define("isProton", "recoGlobalSimPdg==2212")
    //.Define("isPionP", "recoGlobalSimPdg==211")
    //.Define("isPionM", "recoGlobalSimPdg==-211")
    //.Define("isPions", "abs(recoGlobalSimPdg)==211")
    .Define("isKaonP", "recoGlobalSimPdg==321")
    .Define("isKaonM", "recoGlobalSimPdg==-321")
    .Define("isKaons", "abs(recoGlobalSimPdg)==321")
    .Define("recPt", getPt, {"recoGlobalMom"})
    .Define("recP", getP, {"recoGlobalMom"})
    .Define("recRigidity", getRigidity, {"recoGlobalMom", "recoGlobalCharge"})
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
    .Define("nsigDedxProton", [f1_dedx_p, f1_dedx_sigm_p]( RVec<float> vec_Pq, RVec<float> vec_dedx ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_dedx.size() );
      for( int i=0; i<vec_dedx.size(); ++i){
        auto pq = vec_Pq.at(i);
        auto dedx = vec_dedx.at(i);
        if (pq<0.3)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (dedx - f1_dedx_p->Eval(pq))/(f1_dedx_p->Eval(pq)*f1_dedx_sigm_p->Eval(pq)) );
      }
      return vec_nsig;
    },{"recRigidity", "recTpcDedx"})
    .Define("nsigM2Proton", [f1_m2_sigm_p]( RVec<float> vec_Pq, RVec<float> vec_m2, RVec<bool> vec_tof ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_m2.size() );
      for( int i=0; i<vec_m2.size(); ++i){
        auto tof = vec_tof.at(i);
        auto pq = vec_Pq.at(i);
        auto m2 = vec_m2.at(i);
        auto m2_p = pow(TDatabasePDG::Instance()->GetParticle(2212)->Mass(), 2);
        if (pq<0.3 || !tof)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (m2 - m2_p)/f1_m2_sigm_p->Eval(pq) );
      }
      return vec_nsig;
    },{"recRigidity", "recTofMass2", "isTOF"})
    .Define("nsigDedxPionP", [f1_dedx_pip, f1_dedx_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_dedx ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_dedx.size() );
      for( int i=0; i<vec_dedx.size(); ++i){
        auto pq = vec_Pq.at(i);
        auto dedx = vec_dedx.at(i);
        if (abs(pq)<0.1 || pq<0.)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (dedx - f1_dedx_pip->Eval(pq))/(f1_dedx_pip->Eval(pq)*f1_dedx_sigm_pip->Eval(pq)) );
      }
      return vec_nsig;
    },{"recRigidity", "recTpcDedx"})
    .Define("nsigM2PionP", [f1_m2_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_m2, RVec<bool> vec_tof ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_m2.size() );
      for( int i=0; i<vec_m2.size(); ++i){
        auto tof = vec_tof.at(i);
        auto pq = vec_Pq.at(i);
        auto m2 = vec_m2.at(i);
        auto m2_p = pow(TDatabasePDG::Instance()->GetParticle(211)->Mass(), 2);
        if (abs(pq)<0.1 || pq<0. || !tof)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (m2 - m2_p)/f1_m2_sigm_pip->Eval(pq) );
      }
      return vec_nsig;
    },{"recRigidity", "recTofMass2", "isTOF"})
    .Define("nsigDedxPionM", [f1_dedx_pip, f1_dedx_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_dedx ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_dedx.size() );
      for( int i=0; i<vec_dedx.size(); ++i){
        auto pq = (float)-1.*vec_Pq.at(i);
        auto dedx = vec_dedx.at(i);
        if (abs(pq)<0.1 || pq<0.)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (dedx - f1_dedx_pip->Eval(pq))/(f1_dedx_pip->Eval(pq)*f1_dedx_sigm_pip->Eval(pq)) );
      }
      return vec_nsig;
    },{"recRigidity", "recTpcDedx"})
    .Define("nsigM2PionM", [f1_m2_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_m2, RVec<bool> vec_tof ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_m2.size() );
      for( int i=0; i<vec_m2.size(); ++i){
        auto tof = vec_tof.at(i);
        auto pq = (float)-1.*vec_Pq.at(i);
        auto m2 = vec_m2.at(i);
        auto m2_p = pow(TDatabasePDG::Instance()->GetParticle(211)->Mass(), 2);
        if (abs(pq)<0.1 || pq<0. || !tof)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (m2 - m2_p)/f1_m2_sigm_pip->Eval(pq) );
      }
      return vec_nsig;
    },{"recRigidity", "recTofMass2", "isTOF"})
    .Define("nsigDedxPions", [f1_dedx_pip, f1_dedx_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_dedx ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_dedx.size() );
      for( int i=0; i<vec_dedx.size(); ++i){
        auto pq = abs(vec_Pq.at(i));
        auto dedx = vec_dedx.at(i);
        if (pq<0.1)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (dedx - f1_dedx_pip->Eval(pq))/(f1_dedx_pip->Eval(pq)*f1_dedx_sigm_pip->Eval(pq)) );
      }
      return vec_nsig;
    },{"recRigidity", "recTpcDedx"})
    .Define("nsigM2Pions", [f1_m2_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_m2, RVec<bool> vec_tof ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_m2.size() );
      for( int i=0; i<vec_m2.size(); ++i){
        auto tof = vec_tof.at(i);
        auto pq = abs(vec_Pq.at(i));
        auto m2 = vec_m2.at(i);
        auto m2_p = pow(TDatabasePDG::Instance()->GetParticle(211)->Mass(), 2);
        if (pq<0.1 || !tof)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (m2 - m2_p)/f1_m2_sigm_pip->Eval(pq) );
      }
      return vec_nsig;
    },{"recRigidity", "recTofMass2", "isTOF"})
    .Define("isProton", []( RVec<float> vec_nsigDedx_p, RVec<float> vec_nsigM2_p, RVec<float> vec_nsigDedx_pi, RVec<float> vec_nsigM2_pi ){
      RVec<int> vec_pid;
      vec_pid.reserve( vec_nsigDedx_p.size() );
      for( int i=0; i<vec_nsigDedx_p.size(); ++i){
        auto dedx_p = vec_nsigDedx_p.at(i);
        auto m2_p = vec_nsigM2_p.at(i);
        auto dedx_pi = vec_nsigDedx_pi.at(i);
        auto m2_pi = vec_nsigM2_pi.at(i);
        if ( sqrt(pow(dedx_p,2)+pow(m2_p,2)) < 2. &&
             sqrt(pow(dedx_pi,2)+pow(m2_pi,2)) > 3.)
          vec_pid.push_back(1);
        else
          vec_pid.push_back(0);
      }
      return vec_pid;
    },{"nsigDedxProton", "nsigM2Proton", "nsigDedxPionP", "nsigM2PionP"})
    .Define("isPionP", []( RVec<float> vec_nsigDedx_p, RVec<float> vec_nsigM2_p, RVec<float> vec_nsigDedx_pi, RVec<float> vec_nsigM2_pi ){
      RVec<int> vec_pid;
      vec_pid.reserve( vec_nsigDedx_p.size() );
      for( int i=0; i<vec_nsigDedx_p.size(); ++i){
        auto dedx_p = vec_nsigDedx_p.at(i);
        auto m2_p = vec_nsigM2_p.at(i);
        auto dedx_pi = vec_nsigDedx_pi.at(i);
        auto m2_pi = vec_nsigM2_pi.at(i);
        if ( sqrt(pow(dedx_p,2)+pow(m2_p,2)) > 3. &&
             sqrt(pow(dedx_pi,2)+pow(m2_pi,2)) < 2.)
          vec_pid.push_back(1);
        else
          vec_pid.push_back(0);
      }
      return vec_pid;
    },{"nsigDedxProton", "nsigM2Proton", "nsigDedxPionP", "nsigM2PionP"})
    .Define("isPionM", []( RVec<float> vec_nsigDedx, RVec<float> vec_nsigM2 ){
      RVec<int> vec_pid;
      vec_pid.reserve( vec_nsigDedx.size() );
      for( int i=0; i<vec_nsigDedx.size(); ++i){
        auto dedx = vec_nsigDedx.at(i);
        auto m2 = vec_nsigM2.at(i);
        if ( sqrt(pow(dedx,2)+pow(m2,2)) < 2. )
          vec_pid.push_back(1);
        else
          vec_pid.push_back(0);
      }
      return vec_pid;
    },{"nsigDedxPionM", "nsigM2PionM"})
    .Define("isPions", []( RVec<float> vec_nsigDedx_p, RVec<float> vec_nsigM2_p, RVec<float> vec_nsigDedx_pi, RVec<float> vec_nsigM2_pi ){
      RVec<int> vec_pid;
      vec_pid.reserve( vec_nsigDedx_p.size() );
      for( int i=0; i<vec_nsigDedx_p.size(); ++i){
        auto dedx_p = vec_nsigDedx_p.at(i);
        auto m2_p = vec_nsigM2_p.at(i);
        auto dedx_pi = vec_nsigDedx_pi.at(i);
        auto m2_pi = vec_nsigM2_pi.at(i);
        if ( sqrt(pow(dedx_p,2)+pow(m2_p,2)) > 3. &&
             sqrt(pow(dedx_pi,2)+pow(m2_pi,2)) < 2.)
          vec_pid.push_back(1);
        else
          vec_pid.push_back(0);
      }
      return vec_pid;
    },{"nsigDedxProton", "nsigM2Proton", "nsigDedxPions", "nsigM2Pions"})
    .Define("recPtGoodProton", trGoodPid, {"recPt", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recPtGoodPionP", trGoodPid, {"recPt", "isGoodTrackPionP", "isPrimary", "isPionP"})
    .Define("recPtGoodPionM", trGoodPid, {"recPt", "isGoodTrackPionM", "isPrimary", "isPionM"})
    .Define("recPtGoodPions", trGoodPid, {"recPt", "isGoodTrackPions", "isPrimary", "isPions"})
    .Define("recPtGoodKaonP", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isKaonP"})
    .Define("recPtGoodKaonM", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isKaonM"})
    .Define("recPtGoodKaons", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("recYGoodProton", trGoodPid, {"recY", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recYGoodPionP", trGoodPid, {"recY", "isGoodTrackPionP", "isPrimary", "isPionP"})
    .Define("recYGoodPionM", trGoodPid, {"recY", "isGoodTrackPionM", "isPrimary", "isPionM"})
    .Define("recYGoodPions", trGoodPid, {"recY", "isGoodTrackPions", "isPrimary", "isPions"})
    .Define("recYGoodKaonP", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isKaonP"})
    .Define("recYGoodKaonM", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isKaonM"})
    .Define("recYGoodKaons", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("recEtaGoodProton", trGoodPid,{"recEta", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recEtaGoodPionP", trGoodPid, {"recEta", "isGoodTrackPionP", "isPrimary", "isPionP"})
    .Define("recEtaGoodPionM", trGoodPid, {"recEta", "isGoodTrackPionM", "isPrimary", "isPionM"})
    .Define("recEtaGoodPions", trGoodPid, {"recEta", "isGoodTrackPions", "isPrimary", "isPions"})
    .Define("recEtaGoodKaonP", trGoodPid, {"recEta", "isGoodTrack", "isPrimary", "isKaonP"})
    .Define("recEtaGoodKaonM", trGoodPid, {"recEta", "isGoodTrack", "isPrimary", "isKaonM"})
    .Define("recEtaGoodKaons", trGoodPid, {"recEta", "isGoodTrack", "isPrimary", "isKaons"})
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
    .Define("simPtGoodProton", trPrimPid, {"simPt", "isSimPrimary", "isSimProton"})
    .Define("simPtGoodPionP",  trPrimPid, {"simPt", "isSimPrimary", "isSimPionP"})
    .Define("simPtGoodPionM",  trPrimPid, {"simPt", "isSimPrimary", "isSimPionM"})
    .Define("simPtGoodPions",  trPrimPid, {"simPt", "isSimPrimary", "isSimPions"})
    .Define("simPtGoodKaonP",  trPrimPid, {"simPt", "isSimPrimary", "isSimKaonP"})
    .Define("simPtGoodKaonM",  trPrimPid, {"simPt", "isSimPrimary", "isSimKaonM"})
    .Define("simPtGoodKaons",  trPrimPid, {"simPt", "isSimPrimary", "isSimKaons"})
    .Define("simYGoodProton",  trPrimPid, {"simY",  "isSimPrimary", "isSimProton"})
    .Define("simYGoodPionP",   trPrimPid, {"simY",  "isSimPrimary", "isSimPionP"})
    .Define("simYGoodPionM",   trPrimPid, {"simY",  "isSimPrimary", "isSimPionM"})
    .Define("simYGoodPions",   trPrimPid, {"simY",  "isSimPrimary", "isSimPions"})
    .Define("simYGoodKaonP",   trPrimPid, {"simY",  "isSimPrimary", "isSimKaonP"})
    .Define("simYGoodKaonM",   trPrimPid, {"simY",  "isSimPrimary", "isSimKaonM"})
    .Define("simYGoodKaons",   trPrimPid, {"simY",  "isSimPrimary", "isSimKaons"})
    .Define("simEtaGoodProton", trPrimPid, {"simEta", "isSimPrimary", "isSimProton"})
    .Define("simEtaGoodPionP",  trPrimPid, {"simEta", "isSimPrimary", "isSimPionP"})
    .Define("simEtaGoodPionM",  trPrimPid, {"simEta", "isSimPrimary", "isSimPionM"})
    .Define("simEtaGoodPions",  trPrimPid, {"simEta", "isSimPrimary", "isSimPions"})
    .Define("simEtaGoodKaonP",  trPrimPid, {"simEta", "isSimPrimary", "isSimKaonP"})
    .Define("simEtaGoodKaonM",  trPrimPid, {"simEta", "isSimPrimary", "isSimKaonM"})
    .Define("simEtaGoodKaons",  trPrimPid, {"simEta", "isSimPrimary", "isSimKaons"})
    .Define("simPtGoodFHCalProton",  trPrimFHCalPid, {"simPt", "isSimPrimary", "isSimProton", "simHasHitFHCal"})
    .Define("simPtGoodFHCalPionP",  trPrimFHCalPid, {"simPt", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simPtGoodFHCalPionM",  trPrimFHCalPid, {"simPt", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simPtGoodFHCalKaonP",  trPrimFHCalPid, {"simPt", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("simPtGoodFHCalKaonM",  trPrimFHCalPid, {"simPt", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("simYGoodFHCalProton",  trPrimFHCalPid, {"simY", "isSimPrimary", "isSimProton", "simHasHitFHCal"})
    .Define("simYGoodFHCalPionP",  trPrimFHCalPid, {"simY", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simYGoodFHCalPionM",  trPrimFHCalPid, {"simY", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simYGoodFHCalKaonP",  trPrimFHCalPid, {"simY", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("simYGoodFHCalKaonM",  trPrimFHCalPid, {"simY", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("simEtaGoodFHCalProton",  trPrimFHCalPid, {"simEta", "isSimPrimary", "isSimProton", "simHasHitFHCal"})
    .Define("simEtaGoodFHCalPionP",  trPrimFHCalPid, {"simEta", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simEtaGoodFHCalPionM",  trPrimFHCalPid, {"simEta", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simEtaGoodFHCalKaonP",  trPrimFHCalPid, {"simEta", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("simEtaGoodFHCalKaonM",  trPrimFHCalPid, {"simEta", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("recDPtGoodProton", getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recDPtGoodPionP",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrackPionP",  "isPrimary", "isPionP"}) 
    .Define("recDPtGoodPionM",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrackPionM",  "isPrimary", "isPionM"}) 
    .Define("recDPtGoodPions",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrackPions",  "isPrimary", "isPions"}) 
    .Define("recDPtGoodKaonP",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonP"}) 
    .Define("recDPtGoodKaonM",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonM"}) 
    .Define("recDPtGoodKaons",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("recDYGoodProton", getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recDYGoodPionP",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrackPionP",  "isPrimary", "isPionP"}) 
    .Define("recDYGoodPionM",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrackPionM",  "isPrimary", "isPionM"}) 
    .Define("recDYGoodPions",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrackPions",  "isPrimary", "isPions"}) 
    .Define("recDYGoodKaonP",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonP"}) 
    .Define("recDYGoodKaonM",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonM"}) 
    .Define("recDYGoodKaons",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("recDEtaGoodProton", getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recDEtaGoodPionP",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrackPionP",  "isPrimary", "isPionP"}) 
    .Define("recDEtaGoodPionM",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrackPionM",  "isPrimary", "isPionM"}) 
    .Define("recDEtaGoodPions",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrackPions",  "isPrimary", "isPions"}) 
    .Define("recDEtaGoodKaonP",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonP"}) 
    .Define("recDEtaGoodKaonM",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonM"}) 
    .Define("recDEtaGoodKaons",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("recDPhiGoodProton", getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recDPhiGoodPionP",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrackPionP",  "isPrimary", "isPionP"}) 
    .Define("recDPhiGoodPionM",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrackPionM",  "isPrimary", "isPionM"}) 
    .Define("recDPhiGoodPions",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrackPions",  "isPrimary", "isPions"}) 
    .Define("recDPhiGoodKaonP",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonP"}) 
    .Define("recDPhiGoodKaonM",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonM"}) 
    .Define("recDPhiGoodKaons",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaons"})
  ;

  dd.Foreach([](ULong64_t evtId){if (evtId % 100 == 0) cout << "\r" << evtId;}, {"rdfentry_"});

  // Make lists of histograms for QA
  hists.push_back(dd.Histo1D({"h1_recVtx_X","Reconstructed vertex X;x (cm)",500,-1,1}, "recoPrimVtxX")); 
  hists.push_back(dd.Histo1D({"h1_recVtx_Y","Reconstructed vertex Y;y (cm)",500,-1,1}, "recoPrimVtxY")); 
  hists.push_back(dd.Histo1D({"h1_recVtx_Z","Reconstructed vertex Z;z (cm)",500,-1,1}, "recoPrimVtxZ"));
  hists.push_back(dd.Histo1D({"h1_refMult", "Reconstructed N_{ch};N_{ch}",1000,0.,1000.}, "refMult"));
  hists.push_back(dd.Histo1D({"h1_simB", "Simulated b;b, fm",200,0.,20.}, "simB"));
  profs.push_back(dd.Profile1D({"p1_DPt_recNhits_proton", "Pt-resolution for reconstructed protons vs. Nhits;N_{hits};#deltap_{T}", 50, 0., 50., -998., 1000.}, "recoGlobalNhits", "recDPtGoodProton"));
  profs.push_back(dd.Profile1D({"p1_DPt_recDCA_proton", "Pt-resolution for reconstructed protons vs. DCA;DCA (cm);#deltap_{T}", 50, 0., 5., -998., 1000.}, "recDca", "recDPtGoodProton"));
  profs.push_back(dd.Profile1D({"p1_DPt_recChi2_proton", "Pt-resolution for reconstructed protons vs. Chi2;#chi_{2}/ndf;#deltap_{T}", 5000, 0., 5000., -998., 1000.}, "recoGlobalChi2", "recDPtGoodProton"));
  profs.push_back(dd.Profile1D({"p1_DPt_recNhits_pionP", "Pt-resolution for reconstructed pions (#pi^{+}) vs. Nhits;N_{hits};#deltap_{T}", 50, 0., 50., -998., 1000.}, "recoGlobalNhits", "recDPtGoodPionP"));
  profs.push_back(dd.Profile1D({"p1_DPt_recDCA_pionP", "Pt-resolution for reconstructed pions (#pi^{+}) vs. DCA;DCA (cm);#deltap_{T}", 50, 0., 5., -998., 1000.}, "recDca", "recDPtGoodPionP"));
  profs.push_back(dd.Profile1D({"p1_DPt_recChi2_pionP", "Pt-resolution for reconstructed pions (#pi^{+}) vs. Chi2;#chi_{2}/ndf;#deltap_{T}", 5000, 0., 5000., -998., 1000.}, "recoGlobalChi2", "recDPtGoodPionP"));
  profs.push_back(dd.Profile1D({"p1_DPt_recNhits_pionM", "Pt-resolution for reconstructed pions (#pi^{-}) vs. Nhits;N_{hits};#deltap_{T}", 50, 0., 50., -998., 1000.}, "recoGlobalNhits", "recDPtGoodPionM"));
  profs.push_back(dd.Profile1D({"p1_DPt_recDCA_pionM", "Pt-resolution for reconstructed pions (#pi^{-}) vs. DCA;DCA (cm);#deltap_{T}", 50, 0., 5., -998., 1000.}, "recDca", "recDPtGoodPionM"));
  profs.push_back(dd.Profile1D({"p1_DPt_recChi2_pionM", "Pt-resolution for reconstructed pions (#pi^{-}) vs. Chi2;#chi_{2}/ndf;#deltap_{T}", 5000, 0., 5000., -998., 1000.}, "recoGlobalChi2", "recDPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_recVtx_XY","Reconstructed vertex XY;x (cm);y (cm)",500,-1,1,500,-1,1}, "recoPrimVtxX", "recoPrimVtxY"));
  hists2d.push_back(dd.Histo2D({"h2_simBrefMult", "b vs N_{ch};N_{ch};b, fm",1000,0.,1000.,200,0.,20.}, "refMult", "simB"));
  hists2d.push_back(dd.Histo2D({"h2_recDedxPq", "Reconstructed dEdx vs P/q;p/q, GeV/c;dEdx, a.u.", 1000, -5., 5., 500, 0., 5.e3}, "recRigidity", "recTpcDedx", "isGoodTrack"));
  hists2d.push_back(dd.Histo2D({"h2_recDedxPqProton", "Reconstructed dEdx vs P/q for protons;p/q, GeV/c;dEdx, a.u.", 1000, -5., 5., 500, 0., 5.e3}, "recRigidity", "recTpcDedx", "isProton"));
  hists2d.push_back(dd.Histo2D({"h2_recDedxPqPionP", "Reconstructed dEdx vs P/q for #pi^{+};p/q, GeV/c;dEdx, a.u.", 1000, -5., 5., 500, 0., 5.e3}, "recRigidity", "recTpcDedx", "isPionP"));
  hists2d.push_back(dd.Histo2D({"h2_recDedxPqPionM", "Reconstructed dEdx vs P/q for #pi^{-};p/q, GeV/c;dEdx, a.u.", 1000, -5., 5., 500, 0., 5.e3}, "recRigidity", "recTpcDedx", "isPionM"));
  hists2d.push_back(dd.Histo2D({"h2_recDedxPqPions", "Reconstructed dEdx vs P/q for #pi^{#pm};p/q, GeV/c;dEdx, a.u.", 1000, -5., 5., 500, 0., 5.e3}, "recRigidity", "recTpcDedx", "isPions"));
  hists2d.push_back(dd.Histo2D({"h2_recM2Pq", "Reconstructed m^{2}_{TOF} vs P/q;p/q, GeV/c;m^{2}, (GeV/c^{2})^{2}", 1000, -5., 5., 220, -0.2, 2.}, "recRigidity", "recTofMass2", "isGoodTrack"));
  hists2d.push_back(dd.Histo2D({"h2_recM2PqProton", "Reconstructed m^{2}_{TOF} vs P/q for protons;p/q, GeV/c;m^{2}, (GeV/c^{2})^{2}", 1000, -5., 5., 220, -0.2, 2.}, "recRigidity", "recTofMass2", "isProton"));
  hists2d.push_back(dd.Histo2D({"h2_recM2PqPionP", "Reconstructed m^{2}_{TOF} vs P/q for #pi^{+};p/q, GeV/c;m^{2}, (GeV/c^{2})^{2}", 1000, -5., 5., 220, -0.2, 2.}, "recRigidity", "recTofMass2", "isPionP"));
  hists2d.push_back(dd.Histo2D({"h2_recM2PqPionM", "Reconstructed m^{2}_{TOF} vs P/q for #pi^{-};p/q, GeV/c;m^{2}, (GeV/c^{2})^{2}", 1000, -5., 5., 220, -0.2, 2.}, "recRigidity", "recTofMass2", "isPionM"));
  hists2d.push_back(dd.Histo2D({"h2_recM2PqPions", "Reconstructed m^{2}_{TOF} vs P/q for #pi^{#pm};p/q, GeV/c;m^{2}, (GeV/c^{2})^{2}", 1000, -5., 5., 220, -0.2, 2.}, "recRigidity", "recTofMass2", "isPions"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt", "Reconstructed #eta vs p_{T}", 1000, -5., 5., 500, 0., 5.}, "recEta", "recPt"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_inFHCal_proton", "Simulated protons Ycm-pT in FHCal;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodFHCalProton", "simPtGoodFHCalProton"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_inFHCal_pionP", "Simulated pions (#pi^{+}) Ycm-pT in FHCal;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodFHCalPionP", "simPtGoodFHCalPionP"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_inFHCal_pionM", "Simulated pions (#pi^{+}) Ycm-pT in FHCal;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodFHCalPionM", "simPtGoodFHCalPionM"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_inFHCal_kaonP", "Simulated Kaons (#pi^{+}) Ycm-pT in FHCal;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodFHCalKaonP", "simPtGoodFHCalKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_inFHCal_kaonM", "Simulated Kaons (#pi^{+}) Ycm-pT in FHCal;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodFHCalKaonM", "simPtGoodFHCalKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_inFHCal_proton", "Simulated protons #eta-pT in FHCal;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodFHCalProton", "simPtGoodFHCalProton"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_inFHCal_pionP", "Simulated pions (#pi^{+}) #eta-pT in FHCal;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodFHCalPionP", "simPtGoodFHCalPionP"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_inFHCal_pionM", "Simulated pions (#pi^{+}) #eta-pT in FHCal;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodFHCalPionM", "simPtGoodFHCalPionM"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_inFHCal_kaonP", "Simulated Kaons (#pi^{+}) #eta-pT in FHCal;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodFHCalKaonP", "simPtGoodFHCalKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_inFHCal_kaonM", "Simulated Kaons (#pi^{+}) #eta-pT in FHCal;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodFHCalKaonM", "simPtGoodFHCalKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_proton", "Reconstructed protons Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_pionP", "Reconstructed pions (#pi^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_pionM", "Reconstructed pions (#pi^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_pions", "Reconstructed pions (#pi^{#pm}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPions", "recPtGoodPions"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_kaonP", "Reconstructed kaons (K^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_kaonM", "Reconstructed kaons (K^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_kaons", "Reconstructed kaons (K^{#pm}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaons", "recPtGoodKaons"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_proton", "Reconstructed protons #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodProton", "recPtGoodProton"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_pionP", "Reconstructed pions (#pi^{+}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodPionP", "recPtGoodPionP"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_pionM", "Reconstructed pions (#pi^{-}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodPionM", "recPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_pions", "Reconstructed pions (#pi^{#pm}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodPions", "recPtGoodPions"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_kaonP", "Reconstructed kaons (K^{+}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodKaonP", "recPtGoodKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_kaonM", "Reconstructed kaons (K^{-}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodKaonM", "recPtGoodKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_kaons", "Reconstructed kaons (K^{#pm}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodKaons", "recPtGoodKaons"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_proton", "Simulated protons Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodProton", "simPtGoodProton"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_pionP", "Simulated pions (#pi^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodPionP", "simPtGoodPionP"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_pionM", "Simulated pions (#pi^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodPionM", "simPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_pions", "Simulated pions (#pi^{#pm}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodPions", "simPtGoodPions"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_kaonP", "Simulated kaons (K^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodKaonP", "simPtGoodKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_kaonM", "Simulated kaons (K^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodKaonM", "simPtGoodKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_kaons", "Simulated kaons (K^{#pm}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodKaons", "simPtGoodKaons"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_proton", "Simulated protons #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodProton", "simPtGoodProton"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_pionP", "Simulated pions (#pi^{+}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodPionP", "simPtGoodPionP"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_pionM", "Simulated pions (#pi^{-}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodPionM", "simPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_pions", "Simulated pions (#pi^{#pm}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodPions", "simPtGoodPions"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_kaonP", "Simulated kaons (K^{+}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodKaonP", "simPtGoodKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_kaonM", "Simulated kaons (K^{-}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodKaonM", "simPtGoodKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_kaons", "Simulated kaons (K^{#pm}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodKaons", "simPtGoodKaons"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_proton", "Pt-resolution for reconstructed protons in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton", "recDPtGoodProton"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_pionP", "Pt-resolution for reconstructed pions (#pi^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP", "recDPtGoodPionP"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_pionM", "Pt-resolution for reconstructed pions (#pi^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM", "recDPtGoodPionM"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_pions", "Pt-resolution for reconstructed pions (#pi^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPions", "recPtGoodPions", "recDPtGoodPions"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_kaonP", "Pt-resolution for reconstructed kaons (K^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP", "recDPtGoodKaonP"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_kaonM", "Pt-resolution for reconstructed kaons (K^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM", "recDPtGoodKaonM"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_kaons", "Pt-resolution for reconstructed kaons (K^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaons", "recPtGoodKaons", "recDPtGoodKaons"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_proton", "Phi-resolution for reconstructed protons in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton", "recDPhiGoodProton"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_pionP", "Phi-resolution for reconstructed pions (#pi^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP", "recDPhiGoodPionP"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_pionM", "Phi-resolution for reconstructed pions (#pi^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM", "recDPhiGoodPionM"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_pions", "Phi-resolution for reconstructed pions (#pi^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPions", "recPtGoodPions", "recDPhiGoodPions"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_kaonP", "Phi-resolution for reconstructed kaons (K^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP", "recDPhiGoodKaonP"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_kaonM", "Phi-resolution for reconstructed kaons (K^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM", "recDPhiGoodKaonM"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_kaons", "Phi-resolution for reconstructed kaons (K^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaons", "recPtGoodKaons", "recDPhiGoodKaons"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_proton", "Ycm-resolution for reconstructed protons in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton", "recDYGoodProton"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_pionP", "Ycm-resolution for reconstructed pions (#pi^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP", "recDYGoodPionP"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_pionM", "Ycm-resolution for reconstructed pions (#pi^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM", "recDYGoodPionM"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_pions", "Ycm-resolution for reconstructed pions (#pi^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPions", "recPtGoodPions", "recDYGoodPions"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_kaonP", "Ycm-resolution for reconstructed kaons (K^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP", "recDYGoodKaonP"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_kaonM", "Ycm-resolution for reconstructed kaons (K^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM", "recDYGoodKaonM"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_kaons", "Ycm-resolution for reconstructed kaons (K^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaons", "recPtGoodKaons", "recDYGoodKaons"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_proton", "Eta-resolution for reconstructed protons in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton", "recDEtaGoodProton"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_pionP", "Eta-resolution for reconstructed pions (#pi^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP", "recDEtaGoodPionP"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_pionM", "Eta-resolution for reconstructed pions (#pi^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM", "recDEtaGoodPionM"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_pions", "Eta-resolution for reconstructed pions (#pi^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPions", "recPtGoodPions", "recDEtaGoodPions"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_kaonP", "Eta-resolution for reconstructed kaons (K^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP", "recDEtaGoodKaonP"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_kaonM", "Eta-resolution for reconstructed kaons (K^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM", "recDEtaGoodKaonM"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_kaons", "Eta-resolution for reconstructed kaons (K^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaons", "recPtGoodKaons", "recDEtaGoodKaons"));
  //hists3d.push_back(dd.Histo3D({"h3_recDedxM2Pq", "Reconstructed dEdx vs m^{2}_{TOF} vs P/q;m^{2}, (GeV/c^{2})^{2};dEdx, a.u.;p/q, GeV/c", 220, -0.2, 2., 500, 0., 5.e3, 500, -5., 5.}, "recTofMass2", "recTpcDedx", "recRigidity", "isGoodTrack"));
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
