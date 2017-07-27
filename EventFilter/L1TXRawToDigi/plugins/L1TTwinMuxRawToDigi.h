//-------------------------------------------------
//
/**  \class DTTM7FEDReader
 *
 *   L1 DT TwinMux Raw-to-Digi
 *
 *
 *
 *   C. F. Bedoya -- CIEMAT
 *   G. Codispoti -- INFN Bologna
 *   J. Pazzini   -- INFN Padova
 */
//
//--------------------------------------------------
#ifndef L1TXRAWTODIGI_L1TTWINMUXRAWTODIGI_HH
#define L1TXRAWTODIGI_L1TTWINMUXRAWTODIGI_HH

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTTrackContainer.h"

#include <DataFormats/FEDRawData/interface/FEDRawDataCollection.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>

#include "DataFormats/HcalDigi/interface/HOTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDigi/interface/HcalUnpackerReport.h"
#include "EventFilter/L1TXRawToDigi/interface/HOTPUnpacker.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1Trigger/interface/HOTPDigiTwinMux.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/L1Trigger/interface/HOTwinMuxDigiCollection.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"


#include <string>

class L1TTwinMuxRawToDigi : public edm::EDProducer {

public:

  /// Constructor
  L1TTwinMuxRawToDigi( const edm::ParameterSet& pset );

  /// Destructor
  virtual ~L1TTwinMuxRawToDigi();

  /// Produce digis out of raw data
  void produce( edm::Event & e, const edm::EventSetup& c );

  /// HO needful
  std::vector<HOTPDigiTwinMux>* tphoCont; 
  std::auto_ptr<HcalUnpackerReport> report;  

  struct Collections {
    Collections();
    std::vector<HOTPDigiTwinMux>* tphoCont;
  };

  int FindMipFromHcalFeds(int *ieta, int *iphi); 
  bool dodebug =false;

  /// Generate and fill FED raw data for a full event
  bool fillRawData( edm::Event& e,
            L1MuDTChambPhContainer::Phi_Container& phi_data,
            L1MuDTChambThContainer::The_Container& the_data,
		    L1MuDTChambPhContainer::Phi_Container& phi_out_data,
Collections& colls );

  void processFed( int twinmuxfed, int wheel, std::array<short, 12> twinMuxAmcSec,
           edm::Handle<FEDRawDataCollection> data,
           L1MuDTChambPhContainer::Phi_Container& phi_data,
           L1MuDTChambThContainer::The_Container& the_data,
		   L1MuDTChambPhContainer::Phi_Container& phi_out_data, Collections& colls );

private:
  
  bool debug_;
  size_t nfeds_;
  edm::InputTag DTTM7InputTag_;
  std::vector<int> feds_;
  std::vector<int> wheels_;
  std::vector<long long int> amcsecmap_;
  std::vector < std::array<short, 12> > amcsec_;
  
  unsigned char* LineFED_;

  // utilities
  inline void readline( int & lines, long & dataWord )
  { 
    dataWord = *( (long*) LineFED_ );
    LineFED_ += 8;
    ++lines;
  }

  void calcCRC( long word, int & myC );

  edm::InputTag getDTTM7InputTag() { return DTTM7InputTag_; }
  
  edm::EDGetTokenT<FEDRawDataCollection> Raw_token;

  int normBx(int bx_, int bxCnt_);
  int radAngConversion( int radAng_  );
  int benAngConversion( int benAng_  );

  // HO unpacker
  HOTPUnpacker hotpunpacker_;
  std::vector<int> hofedUnpackList_;
  int hofirstFED_;
  int unpackerMode_,expectedOrbitMessageTime_;
  bool silent_, complainEmptyData_;
  edm::FileInPath inputHOLUTs_;
  int mode_ = 0;
  int sourceIdOffset_;
  std::string electronicsMapLabel_;
  bool silent;

  //HO miniElectronicsMap
  int crate = -99, htr=-99, sector=-99;
  int ring=-99, link=-99, indx=-99, eta=-99, phi=99;
  
  struct HOEmap {
  //    HOEmap() {
  int iCrate; //std::vector<int> iCrate;    
  int iHTR; //std::vector<int> iHtr;
  int iSector; //std::vector<int> iSector;
  signed int iWheel;
  int iChan;
  int iEta; //std::vector<int> iEta;
  int iPhi; //std::vector<int> iPhi;
  int iLink;
  int iBitloc;// std::vector<int> iBitloc;
  // }
};
  std::vector<HOEmap> hoemap;
  
  int samples_ = 4;
  int soi_ = 2;
  
  //collections from unpacking HOFEDs
  std::vector<int> hoPhiCol;
  std::vector<int> hoEtaCol;
  std::vector<int> hosoiCol;
  std::vector<int> hosampleCol;
  std::vector<int> hodatabitCol;

};


#endif
