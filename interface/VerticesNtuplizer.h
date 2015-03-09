#ifndef VerticesNtuplizer_H
#define VerticesNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class VerticesNtuplizer : public CandidateNtuplizer {

public:
  VerticesNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches );
  ~VerticesNtuplizer( void );
  
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::InputTag verticesLabel_   ;
   
   edm::Handle< std::vector<reco::Vertex> >  vertices_;
      
};

#endif // VerticesNtuplizer_H
