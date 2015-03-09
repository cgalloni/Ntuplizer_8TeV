#include "../interface/JetsNtuplizer.h"

//===================================================================================================================        
JetsNtuplizer::JetsNtuplizer( std::vector<edm::InputTag> labels, std::vector<std::string> jecCA8Labels, std::vector<std::string> jecAK5Labels, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , jetsAK5Label_( labels[0] )
   , jetsCA8Label_( labels[1] )
   , jetsCA8prunedLabel_( labels[2] )
   , verticesLabel_( labels[3] )
   , rhoLabel_( labels[4] )
   , flavLabel_( labels[5] )
{
   jecCA8PayloadNames_ = jecCA8Labels;
   jecCA8PayloadNames_.pop_back();
   jecCA8UncName_ = jecCA8Labels.back();
  
   jecAK5PayloadNames_ = jecAK5Labels;
   jecAK5PayloadNames_.pop_back();
   jecAK5UncName_ = jecAK5Labels.back();
  
   initJetCorrFactors();


}

//===================================================================================================================
JetsNtuplizer::~JetsNtuplizer( void )
{
}

//===================================================================================================================
void JetsNtuplizer::initJetCorrFactors( void ){

  std::vector<JetCorrectorParameters> vPar;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecCA8PayloadNames_.begin(), payloadEnd = jecCA8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  
  // Make the FactorizedJetCorrector and Uncertainty
  jecCA8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  jecCA8Unc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecCA8UncName_) );

  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK5PayloadNames_.begin(), payloadEnd = jecAK5PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  
  // Make the FactorizedJetCorrector and Uncertainty
  jecAK5_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  jecAK5Unc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecAK5UncName_) );
        
}

//===================================================================================================================
void JetsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  // std::cout <<"runOnMC "<<runOnMC <<std::endl;
  event.getByLabel(jetsAK5Label_      , jetsAK5_      );
  event.getByLabel(jetsCA8Label_      , jetsCA8_      );
  event.getByLabel(jetsCA8prunedLabel_, jetsCA8pruned_);
  event.getByLabel(rhoLabel_          , rho_          ); 
  event.getByLabel(verticesLabel_     , vertices_     );  
  if ( flavLabel_.label() != "" )
    event.getByLabel(flavLabel_         , theTagByValue );         
 /*here we want to save the jets info*/
   
  nBranches_->njetsAK5 = 0;
  for( unsigned j=0; j<jetsAK5_->size(); ++j ){
  
    nBranches_->njetsAK5++;

    reco::Candidate::LorentzVector uncorrJet = (*jetsAK5_)[j].correctedP4(0);

    jecAK5_->setJetEta( uncorrJet.eta()          );
    jecAK5_->setJetPt ( uncorrJet.pt()           );
    jecAK5_->setJetE  ( uncorrJet.energy()       );
    jecAK5_->setJetA  ( (*jetsAK5_)[j].jetArea() );
    jecAK5_->setRho   ( *(rho_.product())        );
    jecAK5_->setNPV   ( vertices_->size()        );
    double corr = jecAK5_->getCorrection();

    
    jecAK5Unc_->setJetEta( uncorrJet.eta() );
    jecAK5Unc_->setJetPt( corr * uncorrJet.pt() );
    double corrUp = corr * (1 + fabs(jecAK5Unc_->getUncertainty(1)));
    jecAK5Unc_->setJetEta( uncorrJet.eta() );
    jecAK5Unc_->setJetPt( corr * uncorrJet.pt() );
    double corrDown = corr * ( 1 - fabs(jecAK5Unc_->getUncertainty(-1)) );
                
    nBranches_->jetAK5_pt     	    .push_back(corr*uncorrJet.pt());
    nBranches_->jetAK5_eta    	    .push_back((*jetsAK5_)[j].eta());
    nBranches_->jetAK5_mass   	    .push_back((*jetsAK5_)[j].mass());
    nBranches_->jetAK5_phi    	    .push_back((*jetsAK5_)[j].phi());
    nBranches_->jetAK5_e      	    .push_back(corr*uncorrJet.energy());
    nBranches_->jetAK5_jec    	    .push_back(corr);
    nBranches_->jetAK5_jecUp        .push_back(corrUp);
    nBranches_->jetAK5_jecDown      .push_back(corrDown);    
    nBranches_->jetAK5_chf    	    .push_back((*jetsAK5_)[j].chargedHadronEnergyFraction());
    nBranches_->jetAK5_nhf    	    .push_back((*jetsAK5_)[j].neutralHadronEnergyFraction());
    nBranches_->jetAK5_phf    	    .push_back((*jetsAK5_)[j].photonEnergyFraction());
    nBranches_->jetAK5_elf    	    .push_back((*jetsAK5_)[j].electronEnergyFraction());
    nBranches_->jetAK5_muf    	    .push_back((*jetsAK5_)[j].muonEnergyFraction());
    nBranches_->jetAK5_cef    	    .push_back((*jetsAK5_)[j].chargedEmEnergyFraction());
    nBranches_->jetAK5_nef    	    .push_back((*jetsAK5_)[j].neutralEmEnergyFraction());
    nBranches_->jetAK5_chm    	    .push_back((*jetsAK5_)[j].chargedHadronMultiplicity());
    nBranches_->jetAK5_nhm    	    .push_back((*jetsAK5_)[j].neutralHadronMultiplicity());
    nBranches_->jetAK5_npr    	    .push_back((*jetsAK5_)[j].chargedHadronMultiplicity()+(*jetsAK5_)[j].neutralMultiplicity());
    nBranches_->jetAK5_cm     	    .push_back((*jetsAK5_)[j].chargedMultiplicity());
    nBranches_->jetAK5_nm     	    .push_back((*jetsAK5_)[j].neutralMultiplicity());
    nBranches_->jetAK5_nconstituents.push_back((*jetsAK5_)[j].nConstituents());
    nBranches_->jetAK5_flavour	    .push_back((*jetsAK5_)[j].partonFlavour());
    nBranches_->jetAK5_charge 	    .push_back((*jetsAK5_)[j].charge());
    nBranches_->jetAK5_ssv    	    .push_back((*jetsAK5_)[j].bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
    nBranches_->jetAK5_csv    	    .push_back((*jetsAK5_)[j].bDiscriminator("combinedSecondaryVertexBJetTags"));
    nBranches_->jetAK5_tchp   	    .push_back((*jetsAK5_)[j].bDiscriminator("trackCountingHighPurBJetTags"));
    nBranches_->jetAK5_tche   	    .push_back((*jetsAK5_)[j].bDiscriminator("trackCountingHighEffBJetTags"));
    nBranches_->jetAK5_jp     	    .push_back((*jetsAK5_)[j].bDiscriminator("jetProbabilityBJetTags"));
    nBranches_->jetAK5_jbp    	    .push_back((*jetsAK5_)[j].bDiscriminator("jetBProbabilityBJetTags"));
    nBranches_->jetAK5_nSVs   	    .push_back(-99);
    

    //std::cout <<"AK5:  corr " <<  corr <<" corr up " <<   corrUp  << " corr * uncorrJet.pt() " << corr * uncorrJet.pt() << " uncorrJet.pt() "<<  uncorrJet.pt() <<" PAT tuple pt "<< (*jetsAK5_)[j].pt() <<std::endl;


    //const reco::SecondaryVertexTagInfo * svTagInfos = (*jetsAK5_)[j].tagInfoSecondaryVertex("secondaryVertex");
    //nBranches_->jetAK5_nSVs.push_back(svTagInfos->nVertices());
             	       
    //    ===============================
    //        Cesar ID Loose
    //    ===============================
    bool isCJetSelection=false;
    if((fabs((*jetsAK5_)[j].eta())<2.4) &&
       //   ((*jetsAK5_)[j].pt()>30.0) &&
       ( ((*jetsAK5_)[j].chargedMultiplicity() + (*jetsAK5_)[j].neutralMultiplicity())>1  ) && 
       ( ((*jetsAK5_)[j].photonEnergy() / ((*jetsAK5_)[j].jecFactor(0) * (*jetsAK5_)[j].energy() )) <0.99 ) &&
       ( ((*jetsAK5_)[j].neutralHadronEnergyFraction() + (*jetsAK5_)[j].HFHadronEnergyFraction())<0.99  )   &&    
       ( (*jetsAK5_)[j].chargedEmEnergyFraction()<0.99 ) &&
       ( ((*jetsAK5_)[j].muonEnergy() / ( (*jetsAK5_)[j].jecFactor(0) * (*jetsAK5_)[j].energy() ) )<0.99 ) &&
       ( (*jetsAK5_)[j].chargedHadronEnergyFraction()>0 ) )  isCJetSelection=true;
    nBranches_->jetAK5_isCJetSelection .push_back( isCJetSelection);
    
    const reco::GenParticle * genP = (*jetsAK5_)[j].genParton();
    if( genP ) nBranches_->jetAK5_flavour.push_back(abs(genP->pdgId()));
    else nBranches_->jetAK5_flavour.push_back(abs((*jetsAK5_)[j].partonFlavour()));
  }  

  
  nBranches_->njetsCA8 = 0;

  for( unsigned j=0; j<jetsCA8_->size(); ++j ){

    nBranches_->njetsCA8++;	       

    reco::Candidate::LorentzVector uncorrJet = (*jetsCA8_)[j].correctedP4(0);

    jecCA8_->setJetEta( uncorrJet.eta()          );
    jecCA8_->setJetPt ( uncorrJet.pt()           );
    jecCA8_->setJetE  ( uncorrJet.energy()       );
    jecCA8_->setJetA  ( (*jetsCA8_)[j].jetArea() );
    jecCA8_->setRho   ( *(rho_.product())        );
    jecCA8_->setNPV   ( vertices_->size()        );
    double corr = jecCA8_->getCorrection();

    jecCA8Unc_->setJetEta( uncorrJet.eta() );
    jecCA8Unc_->setJetPt( corr * uncorrJet.pt() );
    double corrUp = corr * (1 + fabs(jecCA8Unc_->getUncertainty(1)));
    //  std::cout <<"CA*8:  corr " <<  corr <<" corr up " <<   corrUp  << " corr * uncorrJet.pt() " << corr * uncorrJet.pt() << " uncorrJet.pt() "<<  uncorrJet.pt() <<" PAT tuple pt "<< (*jetsCA8_)[j].pt() <<std::endl;
    jecCA8Unc_->setJetEta( uncorrJet.eta() );
    jecCA8Unc_->setJetPt( corr * uncorrJet.pt() );
    double corrDown = corr * ( 1 - fabs(jecCA8Unc_->getUncertainty(-1)) );
   //  std::cout <<"  corr " <<  corr <<" corr down " <<  corrDown << std::endl;       

    
    nBranches_->jetCA8_pt     	    .push_back(corr*uncorrJet.pt());
    nBranches_->jetCA8_eta    	    .push_back((*jetsCA8_)[j].eta());
    nBranches_->jetCA8_mass   	    .push_back((*jetsCA8_)[j].mass());
    nBranches_->jetCA8_phi    	    .push_back((*jetsCA8_)[j].phi());
    nBranches_->jetCA8_e      	    .push_back(corr*uncorrJet.energy());
    nBranches_->jetCA8_jec    	    .push_back(corr);
    nBranches_->jetCA8_jecUp        .push_back(corrUp);
    nBranches_->jetCA8_jecDown      .push_back(corrDown);
    nBranches_->jetCA8_chf    	    .push_back((*jetsCA8_)[j].chargedHadronEnergyFraction());
    nBranches_->jetCA8_nhf    	    .push_back((*jetsCA8_)[j].neutralHadronEnergyFraction());
    nBranches_->jetCA8_phf    	    .push_back((*jetsCA8_)[j].photonEnergyFraction());
    nBranches_->jetCA8_elf    	    .push_back((*jetsCA8_)[j].electronEnergyFraction());
    nBranches_->jetCA8_muf    	    .push_back((*jetsCA8_)[j].muonEnergyFraction());
    nBranches_->jetCA8_cef    	    .push_back((*jetsCA8_)[j].chargedEmEnergyFraction());
    nBranches_->jetCA8_nef    	    .push_back((*jetsCA8_)[j].neutralEmEnergyFraction());
    nBranches_->jetCA8_chm    	    .push_back((*jetsCA8_)[j].chargedHadronMultiplicity());
    nBranches_->jetCA8_nhm    	    .push_back((*jetsCA8_)[j].neutralHadronMultiplicity());
    nBranches_->jetCA8_npr    	    .push_back((*jetsCA8_)[j].chargedHadronMultiplicity()+(*jetsCA8_)[j].neutralMultiplicity());
    nBranches_->jetCA8_cm     	    .push_back((*jetsCA8_)[j].chargedMultiplicity());
    nBranches_->jetCA8_nm     	    .push_back((*jetsCA8_)[j].neutralMultiplicity());
    nBranches_->jetCA8_nconstituents.push_back((*jetsCA8_)[j].nConstituents());
    nBranches_->jetCA8_flavour	    .push_back((*jetsCA8_)[j].partonFlavour()); 				       
    nBranches_->jetCA8_charge 	    .push_back((*jetsCA8_)[j].charge());					       
    nBranches_->jetCA8_ssv    	    .push_back((*jetsCA8_)[j].bDiscriminator("simpleSecondaryVertexHighPurBJetTags")); 
    nBranches_->jetCA8_csv    	    .push_back((*jetsCA8_)[j].bDiscriminator("combinedSecondaryVertexBJetTags"));      
    nBranches_->jetCA8_tchp   	    .push_back((*jetsCA8_)[j].bDiscriminator("trackCountingHighPurBJetTags"));	       
    nBranches_->jetCA8_tche   	    .push_back((*jetsCA8_)[j].bDiscriminator("trackCountingHighEffBJetTags"));	       
    nBranches_->jetCA8_jp     	    .push_back((*jetsCA8_)[j].bDiscriminator("jetProbabilityBJetTags"));	       
    nBranches_->jetCA8_jbp    	    .push_back((*jetsCA8_)[j].bDiscriminator("jetBProbabilityBJetTags"));
    
    	       
    //    ===============================
    //        Cesar ID Loose
    //    ===============================

    float chf = (*jetsCA8_)[j].chargedHadronEnergyFraction();
    float cemf= (*jetsCA8_)[j].chargedEmEnergyFraction(); 
    float nhf = (*jetsCA8_)[j].neutralHadronEnergyFraction() + (*jetsCA8_)[j].HFHadronEnergyFraction();
    float phf = (*jetsCA8_)[j].photonEnergyFraction();
    float elf = (*jetsCA8_)[j].electronEnergy()/((*jetsCA8_)[j].jecFactor(0) * (*jetsCA8_)[j].energy());
    float muf = (*jetsCA8_)[j].muonEnergyFraction();

    int chm    = (*jetsCA8_)[j].chargedHadronMultiplicity();
    int npr    = (*jetsCA8_)[j].chargedMultiplicity() + (*jetsCA8_)[j].neutralMultiplicity();
    float eta  = fabs((*jetsCA8_)[j].eta());
    float pt   = (*jetsCA8_)[j].pt();
    
    bool idL   = ( (*jetsCA8_)[j].chargedMultiplicity()>0 && npr>1 && phf<0.99 && nhf<0.99 && cemf<0.99 && muf < 0.99 && chf>0); 
    nBranches_->jetCA8_isIdL .push_back(idL);

    bool isCJetSelection=(
			  (!((eta>1.0 && eta<1.5) && ((*jetsCA8_)[j].neutralMultiplicity()==0))) &&
			  (!( (eta>1.0 && eta<1.5) && ((*jetsCA8_)[j].chargedMultiplicity()/(*jetsCA8_)[j].neutralMultiplicity() > 2)))
			  
			  && idL && eta <2.4//  && pt>30.0
			  );
    
    
    nBranches_->jetCA8_isCJetSelection .push_back( isCJetSelection);
    
    //    ===============================
    //       End: Cesar ID Loose CA8
    //    ===============================


    const reco::GenParticle * genP = (*jetsCA8_)[j].genParton();
    if( genP ) nBranches_->jetCA8_flavour.push_back(abs(genP->pdgId()));
    else nBranches_->jetCA8_flavour.push_back(0);

    
    const reco::SecondaryVertexTagInfo * svTagInfos = (*jetsCA8_)[j].tagInfoSecondaryVertex("secondaryVertex");
    nBranches_->jetCA8_nSVs   .push_back(svTagInfos->nVertices());  
    
    nBranches_->jetCA8_tau1.push_back((*jetsCA8_)[j].userFloat("tau1"));	       
    nBranches_->jetCA8_tau2.push_back((*jetsCA8_)[j].userFloat("tau2"));
    nBranches_->jetCA8_tau3.push_back((*jetsCA8_)[j].userFloat("tau3"));	   


    bool isAJetSelection=false;
    float massZ=-9999; float ptZ=-999; bool foundJet=false; float tau21Z=-9999;
   
    float dRmin = 9999.; float mass = 0.;




    for( unsigned int l=0; l<jetsCA8pruned_->size(); ++l ){  
      reco::Candidate::LorentzVector uncorrJet = (*jetsCA8pruned_)[l].correctedP4(0);
      
      jecCA8_->setJetEta( uncorrJet.eta()          );
      jecCA8_->setJetPt ( uncorrJet.pt()           );
      jecCA8_->setJetE  ( uncorrJet.energy()       );
      jecCA8_->setJetA  ( (*jetsCA8pruned_)[l].jetArea() );
      jecCA8_->setRho   ( *(rho_.product())        );
      jecCA8_->setNPV   ( vertices_->size()        );
      double corr = jecCA8_->getCorrection();


      float dRtmp = ROOT::Math::VectorUtil::DeltaR((*jetsCA8_)[j].p4(),uncorrJet);
      if(dRtmp<dRmin && dRtmp<0.8){//matching failed if greater than jet radius
	dRmin=dRtmp;
	mass= corr * uncorrJet.mass();

      }
  
   }
    
    nBranches_->jetCA8_prunedMass.push_back(mass);
    
    if(((*jetsCA8_)[j].muonEnergyFraction()<0.99)  &&
       ((*jetsCA8_)[j].photonEnergyFraction()<0.99)  &&
       ((*jetsCA8_)[j].chargedEmEnergyFraction()<0.99)  &&
       ((*jetsCA8_)[j].neutralHadronEnergyFraction()<0.99)  &&
       ((*jetsCA8_)[j].chargedHadronEnergyFraction()>0.00)  &&
       //  ((*jetsCA8_)[j].pt()>400)  &&
       (fabs((*jetsCA8_)[j].eta())<2.4)  // &&
//        //  ((mass>70 && mass<110))  &&
//        ((*jetsCA8_)[j].userFloat("tau2")/(*jetsCA8_)[j].userFloat("tau1")<=0.75)
       ){
      isAJetSelection=  true;}
    
    nBranches_->jetCA8_isAJetSelection .push_back( isAJetSelection);
    
  }
  
  

  nBranches_->njetsCA8pruned = 0;
  std::vector<float> vSubjetpt     ;
  std::vector<float> vSubjeteta    ;
  std::vector<float> vSubjetmass   ;
  std::vector<float> vSubjetphi    ;
  std::vector<float> vSubjete	   ;
  std::vector<int  > vSubjetcharge ;
  std::vector<int  > vSubjetflavour;
  std::vector<float> vSubjetssv    ;
  std::vector<float> vSubjetcsv    ;
  std::vector<float> vSubjettchp   ;
  std::vector<float> vSubjettche   ; 
  std::vector<float> vSubjetjp     ;
  std::vector<float> vSubjetjbp    ;		
    
  for( unsigned int j=0; j<jetsCA8pruned_->size(); ++j ){

    reco::Candidate::LorentzVector uncorrJet = (*jetsCA8pruned_)[j].correctedP4(0);

    jecCA8_->setJetEta( uncorrJet.eta()          );
    jecCA8_->setJetPt ( uncorrJet.pt()           );
    jecCA8_->setJetE  ( uncorrJet.energy()       );
    jecCA8_->setJetA  ( (*jetsCA8pruned_)[j].jetArea() );
    jecCA8_->setRho   ( *(rho_.product())        );
    jecCA8_->setNPV   ( vertices_->size()        );
    double corr = jecCA8_->getCorrection();
    // std::cout << " AMONG the pruned: uncorrJet.pt() "  <<uncorrJet.pt()  << " corr CA8pruned  " <<corr<< " corr*uncorrJet.mass()" << corr*uncorrJet.mass()<< " Pat tuple: mass"  << (*jetsCA8pruned_)[j].mass()<<std::endl;
    
    jecCA8Unc_->setJetEta( uncorrJet.eta() );
    jecCA8Unc_->setJetPt( corr * uncorrJet.pt() );
    double corrUp = corr * (1 + fabs(jecCA8Unc_->getUncertainty(1)));
    jecCA8Unc_->setJetEta( uncorrJet.eta() );
    jecCA8Unc_->setJetPt( corr * uncorrJet.pt() );
    double corrDown = corr * ( 1 - fabs(jecCA8Unc_->getUncertainty(-1)) );
            
    nBranches_->njetsCA8pruned++;
    vSubjetpt	  .clear();
    vSubjeteta    .clear();    
    vSubjetmass   .clear();   
    vSubjetphi    .clear();    
    vSubjete	  .clear();
    vSubjetcharge .clear();
    vSubjetflavour.clear();
    vSubjetssv    .clear();    
    vSubjetcsv    .clear();    
    vSubjettchp   .clear();   
    vSubjettche   .clear();
    vSubjetjp	  .clear();
    vSubjetjbp    .clear();		 
    nBranches_->jetCA8pruned_pt	    .push_back(corr*uncorrJet.pt());
    nBranches_->jetCA8pruned_eta    .push_back((*jetsCA8pruned_)[j].eta());
    nBranches_->jetCA8pruned_mass   .push_back(corr*uncorrJet.mass());
    nBranches_->jetCA8pruned_phi    .push_back((*jetsCA8pruned_)[j].phi());
    nBranches_->jetCA8pruned_e	    .push_back(corr*uncorrJet.energy());
    nBranches_->jetCA8pruned_jec    .push_back(corr);
    nBranches_->jetCA8pruned_jecUp  .push_back(corrUp);
    nBranches_->jetCA8pruned_jecDown.push_back(corrDown);
    nBranches_->jetCA8pruned_flavour.push_back((*jetsCA8pruned_)[j].partonFlavour()					 );
    nBranches_->jetCA8pruned_charge .push_back((*jetsCA8pruned_)[j].charge()						 );
    nBranches_->jetCA8pruned_ssv    .push_back((*jetsCA8pruned_)[j].bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
    nBranches_->jetCA8pruned_csv    .push_back((*jetsCA8pruned_)[j].bDiscriminator("combinedSecondaryVertexBJetTags"	));
    nBranches_->jetCA8pruned_tchp   .push_back((*jetsCA8pruned_)[j].bDiscriminator("trackCountingHighPurBJetTags"	));
    nBranches_->jetCA8pruned_tche   .push_back((*jetsCA8pruned_)[j].bDiscriminator("trackCountingHighEffBJetTags"	));
    nBranches_->jetCA8pruned_jp	    .push_back((*jetsCA8pruned_)[j].bDiscriminator("jetProbabilityBJetTags"		));
    nBranches_->jetCA8pruned_jbp    .push_back((*jetsCA8pruned_)[j].bDiscriminator("jetBProbabilityBJetTags"		));	 
    nBranches_->nsubjets	    .push_back((*jetsCA8pruned_)[j].numberOfDaughters()  				 );
	   
 
 
    const reco::SecondaryVertexTagInfo * svTagInfos = (*jetsCA8pruned_)[j].tagInfoSecondaryVertex("secondaryVertex");
    int nsubjets = 0;
    nBranches_->jetCA8pruned_nSVs.push_back(svTagInfos->nVertices()); 

    for( unsigned int sj=0; sj<(*jetsCA8pruned_)[j].numberOfDaughters(); ++sj ){

      
      const pat::Jet* subjet = dynamic_cast<const pat::Jet*>((*jetsCA8pruned_)[j].daughter(sj));

      ///Flavor
      if( subjet->pt() < 0.01 ) continue;
      nsubjets++;
      int flavor = 0;
      
      if( flavLabel_.label() != "" ){
      	 TLorentzVector dau; dau.SetPtEtaPhiE(subjet->pt(),subjet->eta(),subjet->phi(),subjet->energy());
      	 float dRmin2 = 9999.; 
      	 for ( reco::JetFlavourMatchingCollection::const_iterator t  = theTagByValue->begin();
      	 						    t != theTagByValue->end();
      	 						    t ++ ) {
      	   edm::RefToBase<reco::Jet> aJet  = (*t).first;
      	   const reco::JetFlavour aFlav = (*t).second;
	   
	   if( aJet.get()->pt() < 0.01 ) continue;

      	   TLorentzVector testJ; testJ.SetPtEtaPhiE(aJet.get()->pt(),aJet.get()->eta(),aJet.get()->phi(),aJet.get()->energy());
	   float dRtmp2 = testJ.DeltaR(dau);
	  
	   if( dRtmp2 < dRmin2 && dRtmp2 < 0.1){
	     dRmin2 = dRtmp2;
	    
	     flavor = aFlav.getFlavour();
	   } 
      	  
      	 }
      }


      vSubjetflavour.push_back(flavor                                                        );
      vSubjetpt     .push_back(subjet->pt()						     ); 
      vSubjeteta    .push_back(subjet->eta()						     ); 
      vSubjetmass   .push_back(subjet->mass()						     ); 
      vSubjetphi    .push_back(subjet->phi()						     ); 
      vSubjete      .push_back(subjet->energy() 					     );
      vSubjetcharge .push_back(subjet->charge() 					     ); 
      vSubjetssv    .push_back(subjet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags")); 
      vSubjetcsv    .push_back(subjet->bDiscriminator("combinedSecondaryVertexBJetTags"     )); 
      vSubjettchp   .push_back(subjet->bDiscriminator("trackCountingHighPurBJetTags"	    )); 
      vSubjettche   .push_back(subjet->bDiscriminator("trackCountingHighEffBJetTags"	    ));
      vSubjetjp     .push_back(subjet->bDiscriminator("jetProbabilityBJetTags"  	    ));
      vSubjetjbp    .push_back(subjet->bDiscriminator("jetBProbabilityBJetTags" 	    ));	      
    }
    nBranches_->subjetCA8pruned_pt     .push_back(vSubjetpt     );  
    nBranches_->subjetCA8pruned_eta    .push_back(vSubjeteta    );  
    nBranches_->subjetCA8pruned_mass   .push_back(vSubjetmass   );  
    nBranches_->subjetCA8pruned_phi    .push_back(vSubjetphi    );  
    nBranches_->subjetCA8pruned_e      .push_back(vSubjete      );
    nBranches_->subjetCA8pruned_charge .push_back(vSubjetcharge ); 
    nBranches_->subjetCA8pruned_flavour.push_back(vSubjetflavour);  
    nBranches_->subjetCA8pruned_ssv    .push_back(vSubjetssv    );  
    nBranches_->subjetCA8pruned_csv    .push_back(vSubjetcsv    );  
    nBranches_->subjetCA8pruned_tchp   .push_back(vSubjettchp   );  
    nBranches_->subjetCA8pruned_tche   .push_back(vSubjettche   );
    nBranches_->subjetCA8pruned_jp     .push_back(vSubjetjp     );
    nBranches_->subjetCA8pruned_jbp    .push_back(vSubjetjbp    );  	      
  }
		
}
