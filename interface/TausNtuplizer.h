#ifndef TausNtuplizer_H
#define TausNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class TausNtuplizer : public CandidateNtuplizer {

public:
   TausNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches );
   ~TausNtuplizer( void );
   
   void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::InputTag tausLabel_   ;
   edm::InputTag eleTausLabel_;
   edm::InputTag muTausLabel_ ;
   edm::InputTag rhoLabel_    ;
   edm::InputTag verticesLabel_    ;

   edm::Handle< std::vector<pat::Tau> > taus_   ;
   edm::Handle< std::vector<pat::Tau> > eleTaus_;
   edm::Handle< std::vector<pat::Tau> > muTaus_ ;
   edm::Handle< double >                rho_     ;
   edm::Handle< std::vector<reco::Vertex> > vertices_;    
};

#endif // TausNtuplizer_H
