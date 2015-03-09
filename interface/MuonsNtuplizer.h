#ifndef MuonsNtuplizer_H
#define MuonsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class MuonsNtuplizer : public CandidateNtuplizer {

public:
   MuonsNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches );
   ~MuonsNtuplizer( void );

   void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );

private:
   edm::InputTag muonsLabel_   ;
   edm::InputTag verticesLabel_;
   edm::InputTag rhoLabel_     ;
   edm::InputTag muonsNotBoostedLabel_   ; 
      
   edm::Handle< std::vector<pat::Muon> >    muons_   ;
   edm::Handle< std::vector<reco::Vertex> > vertices_;
   edm::Handle< double >                    rho_     ;
   edm::Handle< std::vector<pat::Muon> >    muonsNotBoosted_   ;

};

#endif // MuonsNtuplizer_H
