#ifndef skimTree_h
#define skimTree_h

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TEfficiency.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <ROOT/TProcessExecutor.hxx>
#include <ROOT/TSeq.hxx>

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <dirent.h>
#include <memory>


class TreeReaderBase
{
 public:
  TreeReaderBase(){};

  TreeReaderBase(const std::string& path)
  {
    setFile(path);
  };

  ~TreeReaderBase()
  {
    if (file_ && file_->IsOpen()) file_->Close();
  };

  // getters
  Long64_t getEntries() const 
  {
    if (reader_.empty()) { return 0; }
    return reader_.begin()->second->GetEntries(false);
  };

  // setters
  void setFile(const std::string& path)
  {
    file_.reset(TFile::Open(path.c_str(), "READ"));
    if (!file_->IsOpen() || file_->IsZombie()) {
      throw std::runtime_error("[ERROR] File "+path+" failed to open!");
    }
  };

  void setEntry(const Long64_t& i, const bool& checkValue, const bool& onDemand)
  {
    if (i<0) { throw std::runtime_error(Form("[ERROR] Invalid index: %lld", i)); }
    onDemand_ = onDemand;
    index_ = i;
    checkValue_ = checkValue;
    if (!onDemand_) {
      for (const auto& r : reader_) {
        loadEntry(r.first);
      }
    }
  };

 protected:
  void checkValue(ROOT::Internal::TTreeReaderValueBase* value) const
  {
    if (value==NULL) {
      throw std::runtime_error("[ERROR] Value pointer is null");
    }
    else if (value->GetSetupStatus() < 0) {
      throw std::runtime_error(Form("[ERROR] Status %d when setting up reader for %s", value->GetSetupStatus(), value->GetBranchName()));
    }
  };

  virtual void checkValues(const std::string& tree) const
  {
    std::cout << tree << std::endl;
  };

  void setTreeReader(const std::string& name, const std::string& dir)
  {
    if (reader_.find(name)!=reader_.end() || !file_) return;
    reader_[name].reset(new TTreeReader(dir.c_str(), file_.get()));
    if (!reader_.at(name)->GetTree()) {
      throw std::runtime_error("[ERROR] Failed to open "+dir+" !");
    }
    currentEntry_[name] = -1;
  };

  void loadTreeEntry(TTreeReader* r) const
  {
    if (!r) return;
    const auto status = r->SetEntry(index_);
    if (status!=TTreeReader::kEntryValid) {
      std::string msg = "";
      if      (status==TTreeReader::kEntryNotLoaded) { msg = "no entry has been loaded yet"; }
      else if (status==TTreeReader::kEntryNoTree) { msg = "the tree does not exist"; }
      else if (status==TTreeReader::kEntryNotFound) { msg = "the tree entry number does not exist"; }
      else if (status==TTreeReader::kEntryChainSetupError) { msg = "problem in accessing a chain element, e.g. file without the tree"; }
      else if (status==TTreeReader::kEntryChainFileError) { msg = "problem in opening a chain's file"; }
      else if (status==TTreeReader::kEntryDictionaryError) { msg = "problem reading dictionary info from tree"; }
      //else if (status==TTreeReader::kEntryLast) { msg = "last entry was reached"; }
      throw std::runtime_error("[ERROR] Invalid entry: "+msg);
    }
  };

  void loadEntry(const std::string& tree) const
  {
    if (index_!=currentEntry_.at(tree)) {
      loadTreeEntry(reader_.at(tree).get());
      if (checkValue_) { checkValues(tree); }
      const_cast<std::map<std::string, Long64_t >* >(&currentEntry_)->at(tree) = index_;
    }
  };

  std::unique_ptr<TFile> file_;
  std::map<std::string, std::unique_ptr<TTreeReader> > reader_;
  bool checkValue_, onDemand_;
  Long64_t index_;
  std::map<std::string, Long64_t > currentEntry_;
};


class TriggerReader : public TreeReaderBase
{
 public:
  TriggerReader(const std::string& path)
  {
    setFile(path);
    setTreeReader("HltTree", "hltanalysis/HltTree");
    evtI_["runN"].reset(new TTreeReaderValue<int>(*reader_.at("HltTree"), "Run"));
    evtUL_["eventN"].reset(new TTreeReaderValue<ULong64_t>(*reader_.at("HltTree"), "Event"));
  };

  // getters
  std::pair<Long64_t, Long64_t> getEventNumber() const
  {
    if (onDemand_) { loadEntry("HltTree"); }
    if (evtUL_.find("eventN")==evtUL_.end()) { return {0, 0}; }
    return {*evtI_.at("runN")->Get(), *evtUL_.at("eventN")->Get()};
  };
    
  bool passTrigger(const std::string& path) const
  {
    // set entry on demand for HltTree
    if (onDemand_) { loadEntry("HltTree"); }
    return *evtI_.at(path)->Get();
  };
    
  bool isTriggerMatched(const TLorentzVector& p4, const std::string& path) const
  {
    // set entry on demand for TriggerObject
    if (onDemand_) { loadEntry(path);      }
    // define delta R threshold
    const auto dR = (path.rfind("HLT_HIL3",0)==0 ? 0.1 : 0.3);
    const auto dPt = (path.rfind("HLT_HIL1",0)==0 ? 2.0 : 10.0);
    // check trigger objects
    bool isMatch = false;
    for (size_t i=0; i<obj_.at(path).at("eta")->Get()->size(); i++) {
      // compare object momenta
      TLorentzVector trigP4; trigP4.SetPtEtaPhiM(obj_.at(path).at("pt")->Get()->at(i), obj_.at(path).at("eta")->Get()->at(i), obj_.at(path).at("phi")->Get()->at(i), p4.M());
      isMatch = (path.rfind("HLT_HIL1",0)==0 ? std::abs(trigP4.Eta()-p4.Eta()) < 0.2 : trigP4.DeltaR(p4) < dR) && std::abs(trigP4.Pt()-p4.Pt())/p4.Pt() < dPt;
      if (isMatch) break;
    }
    return isMatch;
  };

  // setters
  void addTrigger(const std::string& name)
  {
    if (name.rfind("HLT_",0)!=0) {
      throw std::runtime_error("[ERROR] Invalid trigger name "+name+" !");
    }
    if (evtI_.find(name)!=evtI_.end()) return;
    const auto nm = name.substr(0,name.rfind("_v")+2);
    setTreeReader(name, "hltobject/"+nm);
    evtI_[name].reset(new TTreeReaderValue<int>(*reader_.at("HltTree"), name.c_str()));
    for (const auto& var : {"pt", "eta", "phi"}) {
      obj_[name][var].reset(new TTreeReaderValue<std::vector<double> >(*reader_.at(name), var));
    }
  };

  using TreeReaderBase::setEntry;

  bool setEntry(const std::pair<Long64_t, Long64_t>& evtN, const bool& checkValue=true, const bool& onDemand=true)
  {
    const auto index = reader_.at("HltTree")->GetTree()->GetEntryNumberWithIndex(evtN.first, evtN.second);
    if (index<0) { return false; }
    setEntry(index, checkValue, onDemand);
    return true;
  };

 private:
  void checkValues(const std::string& tree) const
  {
    if (tree=="HltTree") {
      for (const auto& r : evtI_ ) { checkValue(r.second.get()); }
      for (const auto& r : evtUL_) { checkValue(r.second.get()); }
    }
    else if (obj_.find(tree)!=obj_.end()) {
      for (const auto& o : obj_.at(tree)) { checkValue(o.second.get()); }
    }
  };

  std::map<std::string, std::unique_ptr<TTreeReaderValue<int> > > evtI_;
  std::map<std::string, std::unique_ptr<TTreeReaderValue<ULong64_t> > > evtUL_;
  std::map<std::string, std::map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<double> > > > > obj_;
};


class RecoReader : public TreeReaderBase
{
 public:
  RecoReader(const std::string& path)
  {
    setFile(path);
    setTreeReader("ggHiNtuplizerGED", "ggHiNtuplizerGED/EventTree");
    setTreeReader("skimanalysis", "skimanalysis/HltTree");
    setTreeReader("hiEvtAnalyzer", "hiEvtAnalyzer/HiTree");
    setTreeReader("pfcandAnalyzer", "pfcandAnalyzer/pfTree");
    setTreeReader("hiFJRhoAnalyzerFinerBins", "hiFJRhoAnalyzerFinerBins/t");
    initEventInfo();
    initPF();
    initRho();
  };

  typedef std::map<std::string, TLorentzVector> metMap;
  
  class Candidate
  {
   public:
    Candidate() {};
    ~Candidate() {};

    typedef TLorentzVector vector;

    // getters
    vector     p4 () const { return p4_;     };
    int    charge () const { return charge_; };
    double    iso () const { return iso_;    };
    int      cent () const { return cent_;   };
    int     evtID () const { return evtId_;  };
    metMap    met () const { return metM_;   };

    // setters
    void setP4    (const vector&    p4) { p4_ = p4;       };
    void setCharge(const int&      chg) { charge_ = chg;  };
    void setIso   (const double&   iso) { iso_ = iso;     };
    void setCent  (const int&      cnt) { cent_ = cnt;    };
    void setEvtID (const int&    evtId) { evtId_ = evtId; };
    void setMET   (const metMap&  metM) { metM_ = metM;   };
    
   private:
    vector p4_;
    int    charge_;
    double iso_;
    int    cent_;
    int    evtId_;
    metMap metM_;
  };
  struct orderByPt { inline bool operator() (const Candidate& i, const Candidate& j) { return (i.p4().Pt() > j.p4().Pt()); } };
  typedef std::vector<Candidate> CandidateCollection;

   // getters
  Candidate getCandidate(const std::string& type) const
  {
    if (onDemand_) { loadEntry("ggHiNtuplizerGED"); }
    if (getSize(type)==0) return Candidate();
    // loop over leptons
    CandidateCollection candidates;
    for (size_t i=0; i<getSize(type); i++) {
      if (passParticleCut(type, i)) {
	candidates.push_back(Candidate());
	const auto& lepChg = (type=="muon" ? objI_.at("muCharge")->Get()->at(i) : objI_.at("eleCharge")->Get()->at(i));
	candidates.back().setCharge(lepChg);
	candidates.back().setP4(getP4(type, i));
	candidates.back().setIso(getIso(type, i));
      }
    }
    if (candidates.empty()) return Candidate();
    std::sort(candidates.begin(), candidates.end(), orderByPt());
    auto& candidate = *candidates.begin();
    if (candidate.p4().Pt()<25.) return Candidate();
    // extract MET
    const auto metM = getMET();
    // extract event info
    const auto cent = getCentrality();
    int evtID = 0;
    if (evtID==0 && candidates.size()>1 && std::next(candidates.begin(), 1)->p4().Pt()>15) { evtID = 1; } // Drell-Yan
    // add information
    candidate.setMET(metM);
    candidate.setCent(cent);
    candidate.setEvtID(evtID);
    // return candidate
    return candidate;
  };

  bool passEventSelection() const
  {
    if (onDemand_) { loadEntry("skimanalysis"); }
    return *skimI_.at("evtSel")->Get();
  };

  int getCentrality() const
  {
    if (onDemand_) { loadEntry("hiEvtAnalyzer"); }
    return *evtI_.at("hiBin")->Get();
  };

  // setters
  void initBranches(const std::string& type)
  {
    if      (type=="electron") { initElectron(); }
    else if (type=="muon"    ) { initMuon();     }
  };

 private:
  void checkValues(const std::string& tree) const
  {
    if (tree=="ggHiNtuplizerGED") {
      for (const auto& o : objI_) { checkValue(o.second.get()); }
      for (const auto& o : objF_) { checkValue(o.second.get()); }
    }
    else if (tree=="hiEvtAnalyzer") {
      for (const auto& o : evtI_ ) { checkValue(o.second.get()); }
      for (const auto& o : evtUI_) { checkValue(o.second.get()); }
      for (const auto& o : evtUL_) { checkValue(o.second.get()); }
    }
    else if (tree=="skimanalysis") {
      for (const auto& o : skimI_ ) { checkValue(o.second.get()); }
    }
    else if (tree=="pfcandAnalyzer") {
      for (const auto& o : objI_) { checkValue(o.second.get()); }
      for (const auto& o : objF_) { checkValue(o.second.get()); }
    }
    else if (tree=="hiFJRhoAnalyzerFinerBins") {
      for (const auto& o : objD_) { checkValue(o.second.get()); }
    } 
  };

  size_t getSize(const std::string& type) const
  {
    if      (type=="electron") { return (objF_.find("elePt")!=objF_.end() ? objF_.at("elePt")->Get()->size() : 0); }
    else if (type=="muon"    ) { return (objF_.find("muPt" )!=objF_.end() ? objF_.at("muPt" )->Get()->size() : 0); }
    return 0;
  };

  TLorentzVector getP4(const std::string& type, const size_t& i) const
  {
    TLorentzVector p4;
    if      (getSize(type)==0) { return p4; }
    else if (type=="electron") { p4.SetPtEtaPhiM(objF_.at("elePt")->Get()->at(i), objF_.at("eleSCEta")->Get()->at(i), objF_.at("eleSCPhi")->Get()->at(i), 0.000511);   }
    else if (type=="muon"    ) { p4.SetPtEtaPhiM(objF_.at("muPt" )->Get()->at(i), objF_.at("muEta"   )->Get()->at(i), objF_.at("muPhi"   )->Get()->at(i), 0.10565837); }
    return p4;
  };

  double getUE(const double& rho, const double& rho0, const double& a, const double& b, const double& c, const double& d, const double& e) const
  {
    const auto logRho0 = std::log(rho0);
    const double f = (a*rho0*rho0 + b*rho0 + c) - (d*logRho0*logRho0 + e*logRho0);
    const auto logRho = std::log(rho);
    return (rho <= rho0 ? a*rho*rho + b*rho + c : d*logRho*logRho + e*logRho + f);
  };

  double getUE(const double& rho, const std::string& type, const double& coneSize=0.2) const
  {
    if (type=="electron" && coneSize==0.3) {
      return getUE(rho, 119.996446, 0.000941, 0.235722, 2.665041, 54.666635, -466.732790);
    }
    else if (type=="electron" && coneSize==0.2) {
      return getUE(rho, 94.772523, 0.000452, 0.101356, 2.222826, 11.076515, -87.827324);
    }
    else if (type=="muon" && coneSize==0.3) {
      return getUE(rho, 108.767334, 0.002178, 0.054611, 7.039144, 0.000000, 68.794928);
    }
    else if (type=="muon" && coneSize==0.2) {
      return getUE(rho, 120.000000, 0.001890, -0.016016, 5.744579, 0.000000, 39.198389);
    }
    throw std::logic_error(Form("Wrong inputs for UE: %s, %g", type.c_str(), coneSize));
    return 0;
  };
  
  double getIso(const std::string& type, const size_t& i, const double& coneSize=0.2, const double& vetoArea=0.015, const double& pTMin=0.5) const
  { 
    // set entry on demand
    if (onDemand_) { loadEntry("pfcandAnalyzer"); }
    if (onDemand_) { loadEntry("hiFJRhoAnalyzerFinerBins"); }
    // get lepton kinematics
    const auto lepP4 = getP4(type, i);
    // loop over PF particles
    double sumPt = 0.0;
    for (size_t iPF=0; iPF<objF_.at("pfPt")->Get()->size(); iPF++) {
      const auto& pfPt  = objF_.at("pfPt" )->Get()->at(iPF);
      if (pfPt < pTMin) continue;
      const auto& pfEta = objF_.at("pfEta")->Get()->at(iPF);
      const auto& pfPhi = objF_.at("pfPhi")->Get()->at(iPF);
      TLorentzVector pfP4;
      pfP4.SetPtEtaPhiM(pfPt, pfEta, pfPhi, 0);
      double deltaR = pfP4.DeltaR(lepP4);
      if ( (deltaR < coneSize) && (deltaR >= vetoArea) ) { sumPt += pfPt; }
    }
    // compute underlying event correction (rho-based)
    double rho = -1.;
    for (size_t iRho=0; iRho<objD_.at("etaMin")->Get()->size(); iRho++) {
      const auto& rhoVal = objD_.at("rho")->Get()->at(iRho);
      const auto& rhoEtaMin = objD_.at("etaMin")->Get()->at(iRho);
      const auto& rhoEtaMax = objD_.at("etaMax")->Get()->at(iRho);
      if (lepP4.Eta()>=rhoEtaMin && lepP4.Eta()<rhoEtaMax) { rho = rhoVal; break; }
    }
    if (rho<0) { throw std::logic_error(Form("Rho not found for lepton eta: %g", lepP4.Eta())); }
    const auto ue = getUE(rho, type, coneSize);
    // return corrected relative lepton isolation
    return ((sumPt - ue)/lepP4.Pt());
  };
  
  metMap getMET() const
  {
    // set entry on demand
    if (onDemand_) { loadEntry("pfcandAnalyzer"); }
    // get MET
    metMap metM;
    const std::vector<std::string> catList({"Full", "Chg", "NoHF", "ChgPtMin4"});
    for (const auto& cat : catList) { metM[cat] = TLorentzVector(0., 0., 0., 0.); }
    for (size_t iPF=0; iPF<objF_.at("pfPt")->Get()->size(); iPF++) {
      const auto& pfId  = objI_.at("pfId")->Get()->at(iPF);
      const auto& pfPt  = objF_.at("pfPt" )->Get()->at(iPF);
      const auto& pfEta = objF_.at("pfEta")->Get()->at(iPF);
      const auto& pfPhi = objF_.at("pfPhi")->Get()->at(iPF);
      const auto& trkNHit = objF_.at("trkNHit")->Get()->at(iPF);
      TLorentzVector vec; vec.SetPtEtaPhiM(pfPt, 0.0, pfPhi, 0.0);
      for (const auto& cat : catList) {
	auto& met = metM.at(cat);
	if (cat=="Full") { met -= vec; }
	else if (cat=="Chg" && trkNHit>0) { met -= vec; }
	else if (cat=="ChgPtMin4" && trkNHit>0 && pfPt>4.) { met -= vec; }
	else if (cat=="NoHF" && pfId<6 && std::abs(pfEta)<3.0) { met -= vec; }
      }
    }
    return metM;
  };

  bool passParticleCut(const std::string& type, const size_t& i) const
  {
    if      (getSize(type)==0) { return false; }
    else if (type=="electron") { return passElectronCut(i); }
    else if (type=="muon"    ) { return passMuonCut(i);     }
    return false;
  };

  bool passElectronCut(const size_t& i) const
  {
    // set entry on demand
    if (onDemand_) { loadEntry("hiEvtAnalyzer"); }
    // check kinematics
    if (objF_.at("elePt")->Get()->at(i) <= 20.) return false;
    if ((std::abs(objF_.at("eleSCEta")->Get()->at(i)) >= 1.4 && std::abs(objF_.at("eleSCEta")->Get()->at(i)) <= 1.6) ||
        std::abs(objF_.at("eleSCEta")->Get()->at(i)) >= 2.1) return false;
    // use Loose ID working points for PbPb 2018
    float cutSigmaIEtaIEta, cutdEtaSeedAtVtx, cutdPhiAtVtx, cutEoverPInv, cutHoverEBc;
    if (*evtI_.at("hiBin")->Get() <= 60) {
      if (std::abs(objF_.at("eleSCEta")->Get()->at(i)) < 1.479) {
        cutSigmaIEtaIEta = 0.013451;
        cutdEtaSeedAtVtx = 0.003814;
        cutdPhiAtVtx = 0.037586;
        cutEoverPInv = 0.017664;
        cutHoverEBc = 0.161613;
      }
      else {
        cutSigmaIEtaIEta = 0.046571;
        cutdEtaSeedAtVtx = 0.006293;
        cutdPhiAtVtx = 0.118592;
        cutEoverPInv = 0.020135;
        cutHoverEBc = 0.131705;
      }
    }
    else {
      if (std::abs(objF_.at("eleSCEta")->Get()->at(i)) < 1.479) {
        cutSigmaIEtaIEta = 0.010867;
        cutdEtaSeedAtVtx = 0.003284;
        cutdPhiAtVtx = 0.020979;
        cutEoverPInv = 0.077633;
        cutHoverEBc = 0.126826;
      }
      else {
        cutSigmaIEtaIEta = 0.033923;
        cutdEtaSeedAtVtx = 0.006708;
        cutdPhiAtVtx = 0.083766;
        cutEoverPInv = 0.019279;
        cutHoverEBc = 0.097703;
      }
    }
    const bool passID = objF_.at("eleSigmaIEtaIEta_2012")->Get()->at(i) < cutSigmaIEtaIEta && std::abs(objF_.at("eledEtaSeedAtVtx")->Get()->at(i)) < cutdEtaSeedAtVtx &&
                        std::abs(objF_.at("eledPhiAtVtx")->Get()->at(i)) < cutdPhiAtVtx && objF_.at("eleEoverPInv")->Get()->at(i) < cutEoverPInv &&
                        objF_.at("eleHoverEBc")->Get()->at(i) < cutHoverEBc && std::abs(objF_.at("eleIP3D")->Get()->at(i)) < 0.03 && objI_.at("eleMissHits")->Get()->at(i) <= 1;
    return passID; 
  };

  bool passMuonCut(const size_t& i) const
  {
    // check kinematics
    const auto& eta = objF_.at("muEta")->Get()->at(i);
    const auto& pt = objF_.at("muPt")->Get()->at(i);
    const bool inAcc = (std::abs(eta)<2.4 && pt>=15.0);
    if (!inAcc) return false;
    // use Tight ID working points for PbPb 2018
    const bool passID = objI_.at("muIsGlobal")->Get()->at(i) && objI_.at("muIsPF")->Get()->at(i) && objI_.at("muStations")->Get()->at(i) > 1 &&
                        objF_.at("muChi2NDF")->Get()->at(i) < 10. && objI_.at("muMuonHits")->Get()->at(i) > 0 &&
                        objI_.at("muTrkLayers")->Get()->at(i) > 5 && objI_.at("muPixelHits")->Get()->at(i) > 0 &&
                        std::abs(objF_.at("muD0")->Get()->at(i)) < 0.2 && std::abs(objF_.at("muDz")->Get()->at(i)) < 0.5;
    return passID;
  };

  void initEventInfo()
  {
    evtUI_["runN"].reset(new TTreeReaderValue<UInt_t>(*reader_.at("hiEvtAnalyzer"), "run"));
    evtUL_["eventN"].reset(new TTreeReaderValue<ULong64_t>(*reader_.at("hiEvtAnalyzer"), "evt"));
    evtI_["hiBin"].reset(new TTreeReaderValue<int>(*reader_.at("hiEvtAnalyzer"), "hiBin"));
    skimI_["evtSel"].reset(new TTreeReaderValue<int>(*reader_.at("skimanalysis"), "collisionEventSelectionAODv2"));
  };

  void initElectron()
  {
    for (const auto& var : {"eleCharge", "eleMissHits"}) {
      objI_[var].reset(new TTreeReaderValue<std::vector<int> >(*reader_.at("ggHiNtuplizerGED"), var));
    }
    for (const auto& var : {"elePt", "eleSCEta", "eleSCPhi", "eleSigmaIEtaIEta_2012", "eledEtaSeedAtVtx", "eledPhiAtVtx", "eleEoverPInv", "eleHoverEBc", "eleIP3D"}) {
      objF_[var].reset(new TTreeReaderValue<std::vector<float> >(*reader_.at("ggHiNtuplizerGED"), var));
    }
  };

  void initMuon()
  {
    for (const auto& var : {"muCharge", "muIsGlobal", "muIsPF", "muIDTight", "muStations", "muTrkLayers", "muPixelHits", "muMuonHits"}) {
      objI_[var].reset(new TTreeReaderValue<std::vector<int> >(*reader_.at("ggHiNtuplizerGED"), var));
    };
    for (const auto& var : {"muPt", "muEta", "muPhi", "muD0", "muDz", "muChi2NDF"}) {
      objF_[var].reset(new TTreeReaderValue<std::vector<float> >(*reader_.at("ggHiNtuplizerGED"), var));
    }
  };

  void initPF()
  {
    for (const auto& var : {"pfId"}) {
      objI_[var].reset(new TTreeReaderValue<std::vector<int> >(*reader_.at("pfcandAnalyzer"), var));
    };
    for (const auto& var : {"pfPt", "pfEta", "pfPhi", "trkNHit"}) {
      objF_[var].reset(new TTreeReaderValue<std::vector<float> >(*reader_.at("pfcandAnalyzer"), var));
    }
  };

  void initRho()
  {
    for (const auto& var : {"etaMin", "etaMax", "rho"}) {
      objD_[var].reset(new TTreeReaderValue<std::vector<double> >(*reader_.at("hiFJRhoAnalyzerFinerBins"), var));
    }
  };
  
  std::map<std::string, std::unique_ptr<TTreeReaderValue<UInt_t> > > evtUI_;
  std::map<std::string, std::unique_ptr<TTreeReaderValue<ULong64_t> > > evtUL_;
  std::map<std::string, std::unique_ptr<TTreeReaderValue<int> > > evtI_, skimI_;
  std::map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<int> > > > objI_;
  std::map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<float> > > > objF_;
  std::map<std::string, std::unique_ptr<TTreeReaderValue<std::vector<double> > > > objD_;
};


bool existDir(const std::string& dir)
{
  bool exist = false;
  const auto& dirp = gSystem->OpenDirectory(dir.c_str());
  if (dirp) { gSystem->FreeDirectory(dirp); exist = true; }
  return exist;
};


void makeDir(const std::string& dir)
{
  if (existDir(dir)==false){
    std::cout << "[INFO] DataSet directory: " << dir << " doesn't exist, will create it!" << std::endl;  
    gSystem->mkdir(dir.c_str(), kTRUE);
  }
};


#endif // ifndef skimTree_h
