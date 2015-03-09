#include "../interface/GenParticlesNtuplizer.h"
 
//===================================================================================================================        
GenParticlesNtuplizer::GenParticlesNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches ) 
   : CandidateNtuplizer( nBranches )
   , genParticlesLabel_( labels[0] )
{

}

//===================================================================================================================        
GenParticlesNtuplizer::~GenParticlesNtuplizer( void )
{
}

//===================================================================================================================        
  void GenParticlesNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  
    event.getByLabel(genParticlesLabel_ , genParticles_); 
    
   /* here we want to save  gen particles info*/
    std::vector<int> vDau ;
    std::vector<int> vMoth;
    int nMoth = 0;
    int nDau  = 0;  
    for( unsigned p=0; p<genParticles_->size(); ++p ){
      //if( (*genParticles_)[p].status() != 3 ) continue;
      //std::cout << (*genParticles_)[p].pdgId() << std::endl;
      vDau.clear(); vMoth.clear();
      nDau = 0; nMoth = 0;
      //genParticle_p4.push_back(TLorentzVector((*genParticles_)[p].px(),(*genParticles_)[p].py(),(*genParticles_)[p].pz(),(*genParticles_)[p].energy()));
      nBranches_->genParticle_pt    .push_back((*genParticles_)[p].pt()     );
      nBranches_->genParticle_px    .push_back((*genParticles_)[p].px()     );
      nBranches_->genParticle_py    .push_back((*genParticles_)[p].py()     );
      nBranches_->genParticle_pz    .push_back((*genParticles_)[p].pz()     );
      nBranches_->genParticle_eta   .push_back((*genParticles_)[p].eta()    );
      nBranches_->genParticle_mass  .push_back((*genParticles_)[p].mass()   );
      nBranches_->genParticle_phi   .push_back((*genParticles_)[p].phi()    );
      nBranches_->genParticle_e     .push_back((*genParticles_)[p].energy() );
      nBranches_->genParticle_status.push_back((*genParticles_)[p].status() );
      nBranches_->genParticle_pdgId .push_back((*genParticles_)[p].pdgId()  );

      nBranches_->genParticle_vx    .push_back((*genParticles_)[p].vertex().X()     );
      nBranches_->genParticle_vy    .push_back((*genParticles_)[p].vertex().Y()      );
      nBranches_->genParticle_vz    .push_back((*genParticles_)[p].vertex().Z()     );
      

      for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
	//if((*genParticles_)[p].daughter(d)->status() != 3 ) continue;
        vDau.push_back( (*genParticles_)[p].daughter(d)->pdgId() );
	nDau++;
      }
      for( unsigned int m=0; m<(*genParticles_)[p].numberOfMothers(); ++m ){
	//	if( (*genParticles_)[p].mother(m)->status() != 3 ) continue;
        vMoth.push_back( (*genParticles_)[p].mother(m)->pdgId() );
	//	if ((*genParticles_)[p].mother(m)->pdgId()==15 )std::cout <<"tau mother"<< std::endl;
	nMoth++;
      }
      nBranches_->genParticle_nDau  .push_back( nDau  );
      nBranches_->genParticle_nMoth .push_back( nMoth );      
      nBranches_->genParticle_mother.push_back( vMoth );
      nBranches_->genParticle_dau   .push_back( vDau  );      
    }

    //////Comemnted out for Dibosons : 14/10/14
   //  edm::Handle<LHEEventProduct> product; 
//     event.getByLabel(std::string("source"), product);   
//     lhef::HEPEUP hepeup_ = product->hepeup();
//     std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;

//     bool lCheck=false;
//     bool lbarCheck=false;
//     bool vlCheck=false;
//     bool vlbarCheck=false;
//     nBranches_->lheHT = 0.;
//     nBranches_->lheNj = 0;
//     TLorentzVector l,lbar,vl,vlbar,V_tlv;
    
//     for(unsigned int i=0; i<pup_.size(); ++i){
//        int id=hepeup_.IDUP[i]; //pdgId
//        int status = hepeup_.ISTUP[i];
//        int idabs=TMath::Abs(id);
//        if( status == 1 && ( ( idabs == 21 ) || (idabs > 0 && idabs < 7) ) ){ // gluons and quarks
//           nBranches_->lheHT += TMath::Sqrt( TMath::Power(hepeup_.PUP[i][0],2) + TMath::Power(hepeup_.PUP[i][1],2) ); // first entry is px, second py
//           nBranches_->lheNj++;
//        }	

//        if(id==11){ l.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lCheck=true;}
//        if(id==-11){ lbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lbarCheck=true;}
//        if(id==12){ vl.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlCheck=true;}
//        if(id==-12){ vlbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlbarCheck=true;}

//        if(id==13){ l.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lCheck=true;}
//        if(id==-13){ lbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lbarCheck=true;}
//        if(id==14){ vl.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlCheck=true;}
//        if(id==-14){ vlbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlbarCheck=true;}

//        if(id==15){ l.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lCheck=true;}
//        if(id==-15){ lbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lbarCheck=true;}
//        if(id==16){ vl.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlCheck=true;}
//        if(id==-16){ vlbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlbarCheck=true;}

//     }
    
//     if( lCheck && lbarCheck ) V_tlv = l + lbar; // ZtoLL
//     if( vlCheck && vlbarCheck ) V_tlv = vl + vlbar; // ZtoNuNu
//     if( lCheck && vlbarCheck ) V_tlv = l + vlbar; // WToLNu
//     if( lbarCheck && vlCheck ) V_tlv = lbar + vl; // WToLNu
//     nBranches_->lheV_pt = V_tlv.Pt();
   
}

