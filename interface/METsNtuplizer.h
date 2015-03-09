#ifndef METsNtuplizer_H
#define METsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class TFormula;

class METsNtuplizer : public CandidateNtuplizer {

public:
   METsNtuplizer( std::vector<edm::InputTag> labels, std::vector<std::string> jecLabels, std::vector<std::string> corrFormulas, NtupleBranches* nBranches );
   ~METsNtuplizer( void );
   
   void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
   
   void   addTypeICorr      ( edm::Event const & event );
   void   addSysShiftCorr   ( edm::Event const & event );   
   double getJEC            ( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
   double getJECOffset      ( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
   void   initJetCorrFactors( void );
   
private:      
   edm::InputTag METsRawLabel_           ;
   edm::InputTag METsLabel_              ; 
   edm::InputTag rhoLabel_               ;
   edm::InputTag verticesLabel_          ;
   edm::InputTag verticesForMEtCorrLabel_;
   edm::InputTag jetsLabel_              ;
   edm::InputTag muonsLabel_             ;
   
   std::vector<std::string> jetCorrLabel_   ;
   std::vector<std::string> offsetCorrLabel_;
   std::vector<std::string> corrFormulas_   ;
  
   boost::shared_ptr<FactorizedJetCorrector>   jecAK5_         ;
   boost::shared_ptr<FactorizedJetCorrector>   jecOffset_      ;
  
   edm::Handle< std::vector<pat::MET> >      METsRaw_           ;
   edm::Handle< std::vector<pat::MET> >      METs_              ;
   edm::Handle< double >                     rho_               ;
   edm::Handle< std::vector<reco::Vertex> >  vertices_          ; 
   edm::Handle< std::vector<reco::Vertex> >  verticesForMEtCorr_;
   edm::Handle< std::vector<pat::Jet> >      jets_              ;
   edm::Handle< std::vector<pat::Muon> >     muons_             ;
   
   std::map<std::string,double> SysShiftCorrMap_;
   std::map<std::string,double> TypeICorrMap_   ;
   
   TFormula* corrPx_;
   TFormula* corrPy_;
  
};

#endif // METsNtuplizer_H
