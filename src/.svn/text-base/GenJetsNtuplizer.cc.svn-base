#include "../interface/GenJetsNtuplizer.h"

//===================================================================================================================        
GenJetsNtuplizer::GenJetsNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches ) 
   : CandidateNtuplizer( nBranches )
   , genJetsAK5Label_( labels[0] )
   , genJetsAK5nonuLabel_( labels[1] )
   , genJetsCA8Label_( labels[2] )
   , genJetsCA8nonuLabel_( labels[3] )
{

}   

//===================================================================================================================        
GenJetsNtuplizer::~GenJetsNtuplizer( void )
{

}

//===================================================================================================================        
void GenJetsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

    event.getByLabel(genJetsAK5Label_    , genJetsAK5_    );	
    event.getByLabel(genJetsAK5nonuLabel_, genJetsAK5nonu_);
    event.getByLabel(genJetsCA8Label_    , genJetsCA8_    );
    event.getByLabel(genJetsCA8nonuLabel_, genJetsCA8nonu_);
    
   /*here we want to save the gen jets info*/
   nBranches_->ngenjetsAK5 = 0;
   for( unsigned j=0; j<genJetsAK5_->size(); ++j ){
      nBranches_->ngenjetsAK5++;
      nBranches_->genJetAK5_pt  .push_back((*genJetsAK5_)[j].pt()    );  
      nBranches_->genJetAK5_eta .push_back((*genJetsAK5_)[j].eta()   ); 
      nBranches_->genJetAK5_mass.push_back((*genJetsAK5_)[j].mass()  );
      nBranches_->genJetAK5_phi .push_back((*genJetsAK5_)[j].phi()   ); 
      nBranches_->genJetAK5_e   .push_back((*genJetsAK5_)[j].energy());   
   }
   
   for( unsigned j=0; j<genJetsAK5nonu_->size(); ++j ){
      nBranches_->genJetAK5nonu_pt  .push_back((*genJetsAK5nonu_)[j].pt()    );  
      nBranches_->genJetAK5nonu_eta .push_back((*genJetsAK5nonu_)[j].eta()   ); 
      nBranches_->genJetAK5nonu_mass.push_back((*genJetsAK5nonu_)[j].mass()  );
      nBranches_->genJetAK5nonu_phi .push_back((*genJetsAK5nonu_)[j].phi()   ); 
      nBranches_->genJetAK5nonu_e   .push_back((*genJetsAK5nonu_)[j].energy());  
   }     

   nBranches_->ngenjetsCA8 = 0;
   for( unsigned j=0; j<genJetsCA8_->size(); ++j ){
      nBranches_->ngenjetsCA8++;
      nBranches_->genJetCA8_pt  .push_back((*genJetsCA8_)[j].pt()    );  
      nBranches_->genJetCA8_eta .push_back((*genJetsCA8_)[j].eta()   ); 
      nBranches_->genJetCA8_mass.push_back((*genJetsCA8_)[j].mass()  );
      nBranches_->genJetCA8_phi .push_back((*genJetsCA8_)[j].phi()   ); 
      nBranches_->genJetCA8_e   .push_back((*genJetsCA8_)[j].energy());   
   }
   
   for( unsigned j=0; j<genJetsCA8nonu_->size(); ++j ){
      nBranches_->genJetCA8nonu_pt  .push_back((*genJetsCA8nonu_)[j].pt()    );  
      nBranches_->genJetCA8nonu_eta .push_back((*genJetsCA8nonu_)[j].eta()   ); 
      nBranches_->genJetCA8nonu_mass.push_back((*genJetsCA8nonu_)[j].mass()  );
      nBranches_->genJetCA8nonu_phi .push_back((*genJetsCA8nonu_)[j].phi()   ); 
      nBranches_->genJetCA8nonu_e   .push_back((*genJetsCA8nonu_)[j].energy());  
   }
              
}

