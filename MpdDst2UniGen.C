#include <Rtypes.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TStopwatch.h>

#include "FairMCEventHeader.h"
#include "MpdMCTrack.h"

R__ADD_INCLUDE_PATH($VMCWORKDIR)

void MpdDst2UniGen(std::string inFileName, std::string outFileName, bool doLab2CmsBoost=false, float snn=0.)
{
  if (inFileName == "") {
    std::cerr << "Input file path is empty! Abort..." << std::endl;
    return;
  }
  if (outFileName == "") {
    std::cerr << "Output file path is empty! Abort..." << std::endl;
    return;
  }

  TStopwatch timer;
  timer.Start();

  float ekin, ec, pc, betaCM, gammaCM;
  float mp = 0.938272f;
  if (doLab2CmsBoost){
    ec = 0.5 * snn;
    pc = sqrt(ec*ec - mp*mp);
    ekin = (snn + 2.*mp)*(snn - 2.*mp)/(2.*mp);
    betaCM = pc/ec;
    gammaCM = ec/mp;
    std::cout << "Lab -> CMS boost is ON. sqrt(sNN) = " << snn << " GeV, E_kin = " << ekin << "A GeV." << std::endl;
    std::cout << "\tbeta_CM = " << betaCM << ", gamma_CM = " << gammaCM << "." << std::endl;
  }

  TChain *dstTree = new TChain("mpdsim");
  dstTree->Add(inFileName.c_str());
  
  // Activate branches
  MpdEvent *event{nullptr};
  dstTree->SetBranchAddress("MPDEvent.", &event);
  FairMCEventHeader *MCHeader{nullptr};
  dstTree->SetBranchAddress("MCEventHeader.", &MCHeader);
  TClonesArray *fMCTracks{nullptr};
  dstTree->SetBranchAddress("MCTrack", &fMCTracks);

  TFile *fo = new TFile(outFileName.c_str(), "recreate");
  URun uheader;
  UEvent *uevent = new UEvent();
  TTree *tree = new TTree("events", "UniGen model file");
  tree->Branch("event", "UEvent", uevent);
  
  long events = dstTree->GetEntries();
  std::cout << "Number of events in the DST file = " << events << std::endl;

  for (long i=0; i<events; i++){
    dstTree->GetEntry(i);

    if (i%100 == 0) std::cout << "\tevent [" << i << "/" << events << "]" << std::endl;

    uevent->Clear();
    float nes = 0.f, nevstep = 0.f, ntstep = 0.f;
    uevent->SetParameters(MCHeader->GetEventID(), MCHeader->GetB(), MCHeader->GetRotZ(), nes, nevstep, ntstep);

    int mctracks = fMCTracks->GetEntries();
    for (int itr=0; itr<mctracks; itr++){
      int tr_counter = 0;
      MpdMCTrack *mctrack = (MpdMCTrack*)fMCTracks->UncheckedAt(itr);
      if (!mctrack) continue;
      if (mctrack->GetMotherId() != -1) continue;

      int status = 0, parent = 0, parentDecay = 0, mate = 0, decay = 0, child[2] = {0, 0};
      double weight = 1.;
      float px, py, pz, en;
      TLorentzVector position(mctrack->GetStartX(), mctrack->GetStartY(), mctrack->GetStartZ(), mctrack->GetStartT());
      TLorentzVector momentum;
      if (!doLab2CmsBoost){
        px = mctrack->GetPx();
        py = mctrack->GetPy();
        pz = mctrack->GetPz();
        en = mctrack->GetEnergy();
      } else {
        px = mctrack->GetPx();
        py = mctrack->GetPy();
        pz = mctrack->GetPz();
        en = mctrack->GetEnergy();
        pz = gammaCM * (pz - betaCM * en);
        en = sqrt(pow(mctrack->GetMass(),2) + px*px + py*py + pz*pz);
      }
      momentum.SetPxPyPzE(px, py, pz, en); 

      UParticle uparticle(tr_counter, mctrack->GetPdgCode(), status, parent, parentDecay, mate, decay, child, momentum, position, weight);
      
      uevent->AddParticle(uparticle);
      tr_counter++;
    }
    tree->Fill();
  }
  uheader.Write();
  tree->Write();

  timer.Stop();
  timer.Print();
}
