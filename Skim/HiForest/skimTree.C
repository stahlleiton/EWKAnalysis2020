#include "skimTree.h"


void skimTree(const std::string& type="muon")
{
  const std::string inDir = "";//"/eos/cms/store/group/phys_heavyions/anstahll/";
  const std::vector<std::string> inFiles = {inDir+"HiForestAOD_1192.root"};


  // prepare multi-threading
  const int nCores = inFiles.size();
  ROOT::EnableImplicitMT();
  ROOT::TProcessExecutor mpe(nCores);

  TH1::AddDirectory(kFALSE);
  auto fillTree = [=](int idx)
  {
    auto c_start = std::clock();
    const auto& inFile = inFiles[idx];
    const auto& pTag = type;
    
    // define variables
    double lepEta, lepIso;
    TLorentzVector met, metNoHF, metChg, metChgPtMin4, lepP4T;
    int lepId, lepCharge, evtID, cent;
    // define tree
    TTree* t = new TTree("AnaTree", "AnaTree");
    auto& tree = *t;
    tree.Branch("lepP4T",&lepP4T);
    tree.Branch("lepEta",&lepEta,"lepEta/D");
    tree.Branch("lepCharge",&lepCharge,"lepCharge/I");
    tree.Branch("lepIso",&lepIso,"lepIso/D");
    tree.Branch("met",&met);
    tree.Branch("metNoHF",&metNoHF);
    tree.Branch("metChg",&metChg);
    tree.Branch("metChgPtMin4",&metChgPtMin4);
    tree.Branch("evtID",&evtID,"evtID/i");
    tree.Branch("cent",&cent,"cent/i");

    // get offline
    RecoReader recoInfo(inFile);

    // get online information
    TriggerReader triggerInfo(inFile);

    // add trigger path to trigger reader
    const std::string path = (pTag=="muon" ? "HLT_HIL3Mu12_v1" : "HLT_HIEle20Gsf_v1");
    triggerInfo.addTrigger(path);

    // initialize offline information
    recoInfo.initBranches(pTag);

    // fill tree
    const auto nEntries = recoInfo.getEntries();
    for (auto iEntry : ROOT::TSeqUL(nEntries)) {
      if ((iEntry%100000)==0) { std::cout << "[INFO] Core " << idx << ":  Processing event " << iEntry << " / " << nEntries << std::endl; }
      recoInfo.setEntry(iEntry, false, true);
      triggerInfo.setEntry(iEntry, false, true);
      // check that event pass event selection
      if (!recoInfo.passEventSelection()) continue;
      // initialize variables
      lepPt=0; lepEta=0; lepPhi=0; lepIso=0; met=TVector2(); metNoHF=TVector2(); metChg=TVector2(); metChgPtMin4=TVector2();
      lepId=0; lepCharge=0; evtID=0; cent=0;
      // extract candidate
      auto candidate = recoInfo.getCandidate(pTag);
      if (candidate.p4().Pt()<=25.) continue;
      // check trigger decision
      if (!triggerInfo.passTrigger(path)) continue;
      // check if trigger matched
      const auto isMatched = (pTag=="muon" ? triggerInfo.isTriggerMatched(candidate.p4(), path) : true);
      if (!isMatched) continue;
      // check if QCD event
      if (candidate.iso() >= (pTag=="muon" ? -0.060000 : 0.080000)) { candidate.setEvtID(2); }
      // fill tree
      lepP4T.SetPtEtaPhiM(candidate.p4().Pt(), 0.0, candidate.p4().Phi(), candidate.p4().M());
      lepEta = candidate.p4().Eta();
      lepCharge = candidate.charge();
      lepIso = candidate.iso();
      met = candidate.met().at("Full");
      metNoHF = candidate.met().at("NoHF");
      metChg = candidate.met().at("Chg");
      metChgPtMin4 = candidate.met().at("ChgPtMin4");
      evtID = candidate.evtID();
      cent = candidate.cent();
      tree.Fill();
    }
    return t;
  };

  const auto& trees = mpe.Map(fillTree, ROOT::TSeqI(nCores));
  TFile ft(Form("output_%s.root", type.c_str()),"RECREATE");
  TList list;
  for (const auto& tree : trees) { list.Add(tree); }
  TTree *newtree = TTree::MergeTrees(&list);
  newtree->SetName("AnaTree"); 
  newtree->Write("AnaTree", 2);
  ft.Close();
}
