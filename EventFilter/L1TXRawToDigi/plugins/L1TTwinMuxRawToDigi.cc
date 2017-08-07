//-------------------------------------------------
//
//   Class: DTTM7FEDReader
//
//   L1 DT TwinMux Raw-to-Digi
//
//
//
//   Author :
//   C. F. Bedoya  - CIEMAT
//   G. Codispoti -- INFN Bologna
//   J. Pazzini   -- INFN Padova
//   P. Saxena    -- DESY Hamburg

//
//--------------------------------------------------

#include "L1TTwinMuxRawToDigi.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <fstream>
bool dodebug = true;

L1TTwinMuxRawToDigi::L1TTwinMuxRawToDigi(const edm::ParameterSet& pset) :

  debug_( pset.getUntrackedParameter<bool>("debug", false) ), 
  nfeds_(0),
  DTTM7InputTag_( pset.getParameter<edm::InputTag>("DTTM7_FED_Source") ),
  feds_( pset.getUntrackedParameter<std::vector<int> >("feds", std::vector<int>()) ),
  wheels_( pset.getUntrackedParameter<std::vector<int> >("wheels", std::vector<int>())),
  amcsecmap_( pset.getUntrackedParameter<std::vector<long long int> >("amcsecmap", std::vector<long long int>())),
  hotpunpacker_(pset.getUntrackedParameter<int>("HOFirstFED",int(FEDNumbering::MINHCALFEDID)),pset.getParameter<int>("firstSample"),pset.getParameter<int>("lastSample")),
  hofedUnpackList_(pset.getUntrackedParameter<std::vector<int> >("HOFEDs", std::vector<int>())),
  hofirstFED_(pset.getUntrackedParameter<int>("HcalFirstFED", 724)),
  unpackerMode_(pset.getUntrackedParameter<int>("UnpackerMode",0)),
  expectedOrbitMessageTime_(pset.getUntrackedParameter<int>("ExpectedOrbitMessageTime",-1)),
  silent_(pset.getUntrackedParameter<bool>("silent",false)),
  complainEmptyData_(pset.getUntrackedParameter<bool>("ComplainEmptyData",true)),
  inputHOLUTs_(pset.getParameter<edm::FileInPath>("inputHOLUTs")),
  electronicsMapLabel_(pset.getParameter<std::string>("ElectronicsMap"))

{
    
  produces<L1MuDTChambPhContainer>("PhIn").setBranchAlias("PhIn");
  produces<L1MuDTChambThContainer>("ThIn").setBranchAlias("ThIn");
  produces<L1MuDTChambPhContainer>("PhOut").setBranchAlias("PhOut");
  produces<HOTwinMuxDigiCollection>();

  Raw_token = consumes<FEDRawDataCollection> (DTTM7InputTag_);
 
  nfeds_ = feds_.size();
  
  if ( nfeds_ != wheels_.size() )
    throw cms::Exception("TwinMux_unpacker") << "Configuration file error. Size of \'wheels\' and \'feds\' differs.\n";  

  if ( amcsecmap_.size() != wheels_.size() )
    throw cms::Exception("TwinMux_unpacker") << "Configuration file error. Size of \'wheels\' and \'amcsecmap\' differs.\n";  
    
  for (size_t wh_i = 0; wh_i < amcsecmap_.size(); ++wh_i){
    std::array<short, 12> whmap;      
    for (size_t amc_i = 1; amc_i < 13; ++amc_i ){
      short shift = (12-amc_i)*4;
      whmap[amc_i-1] = ( amcsecmap_[wh_i] >> shift ) & 0xF;
    }
    amcsec_.push_back(whmap);
  }
   
  // HO miniElectronicMap                                                                                                                                                  
  edm::LogInfo("TwinMux_unpacker") << "Using ASCII LUTs" << inputHOLUTs_.fullPath() << " for HcalTPGCoderULUT initialization";
  const char* filename = inputHOLUTs_.fullPath().c_str();
  std::ifstream file(filename, std::ios::in);
  assert(file.is_open());

  if(!file.is_open()) {
    std::cout << "Problem opening the LUT file. Program terminating." << std::endl;
    exit(EXIT_FAILURE);
  }

  unsigned int nCol=2160;//180; //2160                                                                                                                                     
  HOEmap hoemap_;
  int extra;
  std::string ex = "new";
  char temp;
  double te;

  if(dodebug)  std::cout<<"will open the HO LUTs file"<<std::endl;
  for(size_t i=0; i < nCol; i++) {

    //  while(!file.eof()) {                                                                                                                                                
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> extra>> extra >> hoemap_.iEta >> hoemap_.iPhi >> extra >> extra >> ex >> ex >> hoemap_.iSector>> extra >> extra >> extra >> extra >> extra >> hoemap_.iChan >>
      temp >> hoemap_.iCrate >> te >> hoemap_.iHTR >> ex >> extra>> extra>> extra>> extra>>extra>>extra >> extra >> extra >> extra >> extra >> extra>> extra >> extra 
	 >> ex>>hoemap_.iLink  >> hoemap_.iBitloc >> hoemap_.iWheel;
    
    hoemap.push_back(hoemap_);
  }

  if (hofedUnpackList_.empty()) {
  
    for (int i=hofirstFED_; i<=FEDNumbering::MAXHCALFEDID; i++) { // 724-731 are HO FEDs                                                                                  
      hofedUnpackList_.push_back(i);  }
  }

  hotpunpacker_.setExpectedOrbitMessageTime(expectedOrbitMessageTime_);
  hotpunpacker_.setMode(unpackerMode_);

  std::ostringstream ss;
  for (unsigned int i=0; i<hofedUnpackList_.size(); i++) {
    ss << hofedUnpackList_[i] << " L1TTwinMuxRawToDigi will unpack HOFEDs ( " << ss.str() << ")";
  }
}

L1TTwinMuxRawToDigi::~L1TTwinMuxRawToDigi(){}

void L1TTwinMuxRawToDigi::produce(edm::Event& e, 
                             const edm::EventSetup& c) {
  std::cout<<"dodebug "<< dodebug << std::endl;
  std::unique_ptr<L1MuDTChambPhContainer> TM7phi_product(new L1MuDTChambPhContainer);
  std::unique_ptr<L1MuDTChambThContainer> TM7the_product(new L1MuDTChambThContainer);
  std::unique_ptr<L1MuDTChambPhContainer> TM7phi_out_product(new L1MuDTChambPhContainer);

  L1MuDTChambPhContainer::Phi_Container phi_data;
  L1MuDTChambThContainer::The_Container the_data;
  L1MuDTChambPhContainer::Phi_Container phi_out_data;

  /////// HO
  hoEtaCol.clear();
  hoPhiCol.clear();
  hosoiCol.clear();
  hosampleCol.clear();
  hodatabitCol.clear();

  //Step-A get input
  edm::Handle<FEDRawDataCollection> rawraw;
  //  e.getByToken(tok_data_,rawraw);
  e.getByToken(Raw_token,rawraw);

  //get the mapping
  edm::ESHandle<HcalDbService> pSetup; 
  c.get<HcalDbRecord>().get( pSetup );
  edm::ESHandle<HcalElectronicsMap> item;
  c.get<HcalElectronicsMapRcd>().get(electronicsMapLabel_, item);
  const HcalElectronicsMap* readoutMap = item.product();
  //  const HcalElectronicsMap*  readoutMap = pSetup->getHcalMapping();
  
  HOTPUnpacker::HOCol col;
  col.etaCont=&hoEtaCol;
  col.phiCont=&hoPhiCol;
  col.sampleCont=&hosampleCol;
  col.soiCont=&hosoiCol;
  col.databitsCont=&hodatabitCol;
  auto report = std::make_unique<HcalUnpackerReport>();

  //Step-B unpack all requested FEDs
  for (std::vector<int>::const_iterator i=hofedUnpackList_.begin(); i!=hofedUnpackList_.end(); i++) {
    const FEDRawData& fed = rawraw->FEDData(*i);
    
    if (fed.size()==0) {
      std::cout<<"EmptyData" << "No data for FED " << *i;
      if (complainEmptyData_) {
	//        if (!silent_) edm::LogWarning("EmptyData") << "No data for FED " << *i;
        if (!silent_) std::cout<<"EmptyData" << "No data for FED " << *i;
        report->addError(*i);
      }
    } else if (fed.size()<8*3) {
      std::cout<<"EmptyData" << "Tiny data " << fed.size() << " for FED " << *i;
      //      if (!silent_) edm::LogWarning("EmptyData") << "Tiny data " << fed.size() << " for FED " << *i;
      if (!silent_) std::cout<<"EmptyData" << "Tiny data " << fed.size() << " for FED " << *i;
      report->addError(*i);
    } else {
      try {
	if(dodebug) std::cout<<"will unpack the HOTP for fed:"<<*i<< " & size: "<<fed.size()<<std::endl;
        hotpunpacker_.unpack(fed,*readoutMap,col, *report,silent_);
        report->addUnpacked(*i);
      } catch (cms::Exception& e) {
	std::cout<<"Unpacking error" << e.what();
	//        if (!silent_) edm::LogWarning("Unpacking error") << e.what();
	if (!silent_) std::cout<<"Unpacking error" << e.what();
        report->addError(*i);
      } catch (...) {
	//  if (!silent_) edm::LogWarning("Unpacking exception");
	if (!silent_) std::cout<<"Unpacking exception";
        report->addError(*i);
      }
    }
  }

  std::vector<HOTPDigiTwinMux> hotp;
  Collections colls;
  colls.tphoCont = &hotp;

  if ( !fillRawData(e, phi_data, the_data, phi_out_data, colls) ) return;

  TM7phi_product->setContainer(phi_data);
  TM7the_product->setContainer(the_data);
  TM7phi_out_product->setContainer(phi_out_data);

  e.put(std::move(TM7phi_product), "PhIn");
  e.put(std::move(TM7the_product), "ThIn");
  e.put(std::move(TM7phi_out_product), "PhOut");

  //HO
  auto hotp_product = std::make_unique<HOTwinMuxDigiCollection>();
  hotp_product->swap_contents(hotp);
  hotp_product->sort();
  e.put(std::move(hotp_product));
}


bool L1TTwinMuxRawToDigi::fillRawData( edm::Event& e,
                                  L1MuDTChambPhContainer::Phi_Container& phi_data,
                                  L1MuDTChambThContainer::The_Container& the_data,
				       L1MuDTChambPhContainer::Phi_Container& phi_out_data,
Collections& colls ) {

  edm::Handle<FEDRawDataCollection> data;
  e.getByToken( Raw_token, data );

  for ( size_t w_i = 0; w_i < nfeds_; ++w_i ) {
    processFed( feds_[w_i], wheels_[w_i], amcsec_[w_i], data, phi_data, the_data, phi_out_data, colls );
  }
  
  return true;
}

int L1TTwinMuxRawToDigi::normBx( int bx_, 
                            int bxCnt_ ){
    
    int bxNorm_ = bx_ - bxCnt_;    
    if ( abs( bxNorm_ ) < 3000 ) return bxNorm_; 
    
    if ( bxNorm_ > 0 ) return bxNorm_ - 3564;
    if ( bxNorm_ < 0 ) return bxNorm_ + 3564;
    
    return -99;
    
}

int L1TTwinMuxRawToDigi::radAngConversion( int radAng_  ) {
    
    if (radAng_>2047) 
        return radAng_-4096;

    return radAng_;
    
}

int L1TTwinMuxRawToDigi::benAngConversion( int benAng_  ) {
    
    if (benAng_>511) 
        return benAng_-1024;

    return benAng_;
    
}



int L1TTwinMuxRawToDigi::FindMipFromHcalFeds(int *ieta, int *iphi) {
  for (unsigned i=0; i<hoEtaCol.size(); i++) {

    if(*ieta ==hoEtaCol[i] &&  *iphi == hoPhiCol[i]) {
      if(dodebug)      std::cout<<"Matching done, (eta,phi,soi,mip) : "<<hoEtaCol[i]<<","<< hoPhiCol[i] << ","<<hosoiCol[i]<< ","<<hodatabitCol[i]<< std::endl;
      if(hodatabitCol[i] > 1  ) return 1;
      else return 0;
    }
  }
  return 0;
}

struct HOUnrolledTP { // parts of an HO trigger primitive, unpacked                
  bool valid,checked;
  int ieta, iphi, mip, sector, index, link, bx;
  unsigned int databits;
  signed int wheel;
  HOUnrolledTP() {
    checked=false;
    ieta=-99;
    iphi=-99;
    bx=-99;
    mip =-99;
    valid=false;    
    wheel=-99;
    sector=-99;
    index=-99;
    link=-99;
  }
  //  void setbit(int i) { databits|=(1<<i); }
};



void L1TTwinMuxRawToDigi::processFed( int twinMuxFed, 
                                 int twinMuxWheel,
                                 std::array<short, 12> twinMuxAmcSec,
                                 edm::Handle<FEDRawDataCollection> data,
                                 L1MuDTChambPhContainer::Phi_Container& phiSegments,
                                 L1MuDTChambThContainer::The_Container& theSegments,
				      L1MuDTChambPhContainer::Phi_Container& phioutSegments,
Collections& colls ) {

  /// Container
  std::vector<long> DTTM7WordContainer;

  /// Debug
  std::ofstream logfile;
  if ( debug_ ) {
    std::ostringstream fname;
    fname << "eventDump_" <<  twinMuxFed << ".txt";
    logfile.open( fname.str() );
  }

  /// Header
  FEDRawData TM7data = data->FEDData(twinMuxFed); 
  if ( TM7data.size() == 0 ) return;

  /// Variables
  LineFED_ = TM7data.data();
  int nline  = 0; // counting already include header
  long dataWord = 0;
  int newCRC = 0xFFFF;

  ///--> Header - line 1 [must start with 0x5]
  readline( nline, dataWord );
  calcCRC( dataWord, newCRC );

  int TM7fedId = ( dataWord >> 8 ) & 0xFFF;  // positions 8 -> 19
  /*** NOT UNPACKED  
  int bunchCnt = ( dataWord >> 20 ) & 0xFFF;  // positions 20 -> 31
  int eventCnt = ( dataWord >> 32 ) & 0xFFFFFF;  // positions 32 -> 55
  ***/
  int BOEevTy  = ( dataWord >> 60 ) & 0xF;  // positions 60 -> 63

  int linecounter = 0;
  if ( debug_ ) logfile << '[' << ++linecounter << "]\t"
                        << std::hex << dataWord << std::dec << "\t|\t"
                        << "BOEevTy " << BOEevTy << '\t'
                        << "TM7fedId "  << TM7fedId << '\n';

  if ( (BOEevTy != 0x5) || ( TM7fedId != twinMuxFed ) ) {
            
    edm::LogWarning("TwinMux_unpacker") << "Not a TM7 of FED " 
                                        << twinMuxFed << " header "
                                        << std::hex << dataWord;
    return;
    
  }

  ///--> Header - line 2
  readline( nline, dataWord );
  calcCRC( dataWord, newCRC );

  std::map<int, int> AMCsizes;
  /*** NOT UNPACKED  
  int orbit = ( dataWord >> 4 ) & 0xFFFFFFFF;  // positions 4 -> 35
  ***/
  int nAMC = ( dataWord >> 52 ) & 0xF;  // positions 52 -> 55

  if ( debug_ ) logfile << '[' << ++linecounter << "]\t" << std::hex
                        << dataWord << std::dec << "\t|\t"
                        << "nAMC " << nAMC << '\n';

  ///--> AMC - line 3 to 3+nAMC
  for ( int j = 0; j < nAMC; ++j ) {
  
    readline( nline, dataWord ); 
    calcCRC( dataWord, newCRC );
   
    int AMCno = (dataWord >> 16 ) & 0xF;  // positions 16 -> 19
    /*** NOT UNPACKED  
    int TM7boardID = dataWord & 0xFFFF;  // positions 0 -> 15
    int bulkno = (dataWord >> 20 ) & 0xFF;  // positions 20 -> 27
    ***/
    if ( (AMCno < 1) || (AMCno > 12) ) {
        edm::LogWarning("TwinMux_unpacker") << "AMCnumber " << std::dec << AMCno 
                                            << " out of range (1-12)";
        return;
    }    

    AMCsizes[AMCno] = ( dataWord >> 32 ) & 0xFFFFFF;  // positions 32 -> 55

    if ( debug_ ) logfile << '[' << ++linecounter << "]\t"
                          << std::hex << dataWord
                          << std::dec << "\t|\t"
                          << "AMCsizes[" << AMCno << "] "
                          << AMCsizes[AMCno]
                          << std::dec << '\n';
  }

  ///--> Store payloads
  std::map<int,int>::iterator AMCiterator = AMCsizes.begin();
  std::map<int,int>::iterator AMCitend = AMCsizes.end();  
  for ( ; AMCiterator != AMCitend; ++AMCiterator ) {
      
    for ( int k=0; k<AMCiterator->second; ++k) {
        
       readline( nline, dataWord );
       calcCRC( dataWord, newCRC);
       DTTM7WordContainer.push_back( dataWord );
    }
  }  

  ///--> Trailer - line 1
  readline( nline, dataWord );
  calcCRC( dataWord, newCRC);

  ///--> Trailer - line 2 [must start with 0xA]

  readline( nline, dataWord );
  calcCRC( dataWord & 0xFFFFFFFF0000FFFF, newCRC); /// needed not to put crc in crc calc

  ///--> AMC trailer - line 2
  int chkEOE = (dataWord >> 60 ) & 0xF;  // positions 60 -> 63
  int CRC = ( dataWord >> 16 ) & 0xFFFF; // positions 17 ->32
  int evtLgth = ( dataWord >> 32 ) & 0xFFFFFF; // positions 33 ->56

  if ( chkEOE != 0xA ) {
      edm::LogWarning("TwinMux_unpacker") << "AMC block closing line " << std::hex << dataWord 
                                          << std::dec << " does not start with 0xA";
      return;
  }    

  if ( debug_ ) logfile << "\tevtLgth " << std::hex
                        << evtLgth << "\tCRC " << CRC << std::dec << '\n';

  if ( nline != evtLgth ) {
    edm::LogWarning("TwinMux_unpacker") << "Number of words read " << std::dec << nline 
                                        << " and event length " << std::dec << evtLgth 
                                        << " differ ";
    return;
  }

  if ( newCRC != CRC ) {
    edm::LogWarning("TwinMux_unpacker") << "Calculated CRC " << std::hex << newCRC 
                                        << " differs from CRC in trailer " << std::hex << CRC;
    return;
  }

  // --> Analyze event 
  std::vector<long>::iterator DTTM7iterator = DTTM7WordContainer.begin();
  std::vector<long>::iterator DTTM7itend = DTTM7WordContainer.end();

  int lcounter = 0;
  for ( ; DTTM7iterator != DTTM7itend; ++DTTM7iterator ) {

    dataWord  = (*DTTM7iterator);
    int dataLenght = (dataWord & 0xFFFFF);         // positions 0 -> 19
    int bxCounter  = (dataWord >> 20 ) & 0xFFF;    // positions 20 -> 31
    int event      = (dataWord >> 32 ) & 0xFFFFFF; // positions 32 -> 55
    int AMC_ID     = (dataWord >> 56 ) & 0xF;      // positions 56 -> 59
    int control    = (dataWord >> 60 ) & 0xF;      // positions 59 -> 63 
    int wheel      = twinMuxWheel;
    
    if( ( AMC_ID < 1 ) or ( AMC_ID > 12 ) ) {
      edm::LogWarning("TwinMux_unpacker") << "%%%%%% AMC_ID OUT OF RANGE \n"
                                          << " TM7fedId "     << TM7fedId
                                          << " AMC_ID "       << AMC_ID;
      break;
    }
    
    int sector     = twinMuxAmcSec[AMC_ID-1];
    
    if( ( sector < 1 ) or ( sector > 12 ) ) {
      if( sector != 15 ) edm::LogWarning("TwinMux_unpacker") << "%%%%%% VALID AMC_ID POINTS TO SECTOR OUT OF RANGE \n"
                                                             << " TM7fedId "     << TM7fedId
                                                             << " AMC_ID "       << AMC_ID
                                                             << " wheel "        << wheel
                                                             << " sector "       << sector;         
      break;
    }

    if ( debug_ ) logfile << '[' << ++lcounter << "]\t"
                          << std::hex << dataWord << std::dec << "\t|\t"
                          << "AMC_ID "     << AMC_ID << '\t' 
                          << "control "    << control   << '\t' 
                          << "event "      << event  << '\t' 
                          << "bxCounter "  << bxCounter  << '\t' 
                          << "dataLenght " << dataLenght << '\n';

    ++DTTM7iterator; // User word empty  /// ==>> increment 2
    if( DTTM7iterator == DTTM7itend ) {
      edm::LogInfo("TwinMux_unpacker") << "TRAILING WORD AS A PAYLOAD END in FED " 
                                       << std::hex << TM7fedId 
                                       << std::hex << dataWord 
                                       << std::dec<< " [it pos " 
                                       << int(DTTM7iterator - DTTM7itend)  << " ]";
      break;
    }

    dataWord = (*DTTM7iterator);
    int boardID   = (dataWord & 0xFFFF); // positions  0 -> 15
    int orbit     = (dataWord >> 16 ) & 0xFFFF; // positions 15 -> 32
    
    if ( DTTM7iterator == DTTM7itend ) {
      edm::LogWarning("TwinMux_unpacker") << "%%%%%% AMC_ID " << AMC_ID
                                          << " control "      << control
                                          << " event "        << event
                                          << " bxCounter "    << bxCounter
                                          << " size "         << dataLenght
                                          << " orbit "        << orbit
                                          << " board "        << boardID
                                          << " AMCsizes "     << AMCsizes[AMC_ID]
                                          << " it pos "       << int(DTTM7iterator - DTTM7itend);
      break;
    }

    if (debug_ ) logfile << '[' << ++lcounter << "]\t" 
                         << std::hex << dataWord 
                         << std::dec << "\t|\t" 
                         << " orbit " << orbit
                         << " board " << boardID << '\n';

    int AMCsize = AMCsizes[AMC_ID] - 1; /// do not consider the trailer
    int bxID =  99;
    int bc0  = -99;
    int bxNr = -99;
    
    /// 2 words already read, last removed because trailer with CRC
    for ( int tm7eventsize = 2; tm7eventsize < AMCsize; ++tm7eventsize ) {
  
      ++DTTM7iterator; /// ==>> increment 3   
      if ( DTTM7iterator == DTTM7itend ) {
          
        edm::LogWarning("TwinMux_unpacker") << "UNEXPECTED END OF PAYLOAD INSIDE CHAMBER DESCRIPTION"
                                          << " [it pos " << int(DTTM7iterator - DTTM7itend)  << " ]" ;
        break;
        
      }

      long dataWordSub = (*DTTM7iterator);
      int selector = ( dataWordSub >> 60 ) & 0xF; // positions 60 -> 63

      if ( selector == 0x4 ) { //TSC word

        bxID = ( dataWordSub >> 48 ) & 0xFFF; // positions 48 -> 60
        bc0  = ( dataWordSub >> 22 ) & 0x1; // positions 22 -> 23
        bxNr = normBx(bxID, bxCounter); /// bx normalized to the bxcounter
		
        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t TSC WORD\t"
                              << "bxID " << bxID << '\t'
                              << "bc0  " << bc0  << '\n';

      }//TSC WORD
 
      else if ( selector == 0x1 ) { //MB1/2 word
 
        int mb2_phi =    ( dataWordSub & 0xFFF);        // positions  0 -> 11
        int mb2_phib =   ( dataWordSub >> 12 ) & 0x3FF; // positions 12 -> 21
        int mb2_qual =   ( dataWordSub >> 22 ) & 0x7;   // positions 22 -> 24
        int mb2_ts2tag = ( dataWordSub >> 28 ) & 0x1;   // positions 28
        /*** NOT UNPACKED  
        int mb2_parity = ( dataWordSub >> 29) & 0x1;    // positions 29
        ***/

        int mb1_phi  =   ( dataWordSub >> 30 ) & 0xFFF; // positions 30 -> 41
        int mb1_phib =   ( dataWordSub >> 42 ) & 0x3FF; // positions 42 -> 51
        int mb1_qual =   ( dataWordSub >> 52 ) & 0x7;   // positions 52 -> 54
        int mb1_ts2tag = ( dataWordSub >> 58 ) & 0x1;   // positions 58
        /*** NOT UNPACKED  
        int mb1_parity = ( dataWordSub >> 59 ) &0x1;    // positions 59
        ***/
        
        int mb1_phi_conv  = radAngConversion(mb1_phi);
        int mb1_phib_conv = benAngConversion(mb1_phib);

        int mb2_phi_conv  = radAngConversion(mb2_phi);
        int mb2_phib_conv = benAngConversion(mb2_phib);

        phiSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                1, mb1_phi_conv, mb1_phib_conv, 
                                                mb1_qual, mb1_ts2tag, bxCounter ) );
        phiSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                2, mb2_phi_conv, mb2_phib_conv, 
                                                mb2_qual, mb2_ts2tag, bxCounter ) );

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t"<< std::hex 
                              << dataWordSub   << std::dec      << "\t|\t"
                              << "mb1_ts2tag " << mb1_ts2tag    << '\t'
                              << "mb1_qual "   << mb1_qual      << '\t'
                              << "mb1_phib "   << mb1_phib_conv << '\t'
                              << "mb1_phi "    << mb1_phi_conv  << '\t'
                              << "mb2_ts2tag " << mb2_ts2tag    << '\t'
                              << "mb2_qual "   << mb2_qual      << '\t'
                              << "mb2_phib "   << mb2_phib_conv << '\t'
                              << "mb2_phi "    << mb2_phi_conv  << '\n';
      }//MB1/2 word
 
      else if ( selector == 0x2 ) {

        int mb4_phi =    ( dataWordSub & 0xFFF);        // positions  0 -> 11
        int mb4_phib =   ( dataWordSub >> 12 ) & 0x3FF; // positions 12 -> 21
        int mb4_qual =   ( dataWordSub >> 22 ) & 0x7;   // positions 22 -> 24
        int mb4_ts2tag = ( dataWordSub >> 28 ) & 0x1;   // positions 28
        /*** NOT UNPACKED  
        int mb4_parity = ( dataWordSub >> 29) & 0x1;    // positions 29
        ***/

        int mb3_phi  =   ( dataWordSub >> 30 ) & 0xFFF; // positions 30 -> 41
        int mb3_phib =   ( dataWordSub >> 42 ) & 0x3FF; // positions 42 -> 51
        int mb3_qual =   ( dataWordSub >> 52 ) & 0x7;   // positions 52 -> 54
        int mb3_ts2tag = ( dataWordSub >> 58 ) & 0x1;   // positions 58
        /*** NOT UNPACKED  
        int mb3_parity = ( dataWordSub >> 59 ) &0x1;    // positions 59
        ***/

        int mb3_phi_conv  = radAngConversion(mb3_phi);
        int mb3_phib_conv = benAngConversion(mb3_phib);

        int mb4_phi_conv  = radAngConversion(mb4_phi);
        int mb4_phib_conv = benAngConversion(mb4_phib);

        phiSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                  3, mb3_phi_conv, mb3_phib_conv, 
                                                  mb3_qual, mb3_ts2tag, bxCounter) );
        phiSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                  4, mb4_phi_conv, mb4_phib_conv, 
                                                  mb4_qual, mb4_ts2tag, bxCounter) );

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t"<< std::hex
                              << dataWordSub   << std::dec      << "\t|\t"
                              << "mb3_ts2tag " << mb3_ts2tag    << '\t'
                              << "mb3_qual "   << mb3_qual      << '\t'
                              << "mb3_phib "   << mb3_phib_conv << '\t'
                              << "mb3_phi "    << mb3_phi_conv  << '\t'
                              << "mb4_ts2tag " << mb4_ts2tag    << '\t'
                              << "mb4_qual "   << mb4_qual      << '\t'
                              << "mb4_phib "   << mb4_phib_conv << '\t'
                              << "mb4_phi "    << mb4_phi_conv  << '\n';

      }//MB3/4 word
  
      else if ( selector == 0x3 ) { //etha word
       
        int posBTI[7], qualBTI[7];
    
        int mb3_eta    = ( dataWordSub & 0xFF );        // positions  0 -> 7
        int mb2_eta    = ( dataWordSub >> 16 ) & 0xFF;  // positions 16 -> 23
        int mb1_eta    = ( dataWordSub >> 40 ) & 0xFF;  // positions 40 -> 47

        int mb3_eta_HQ = ( dataWordSub >> 8  ) & 0xFF;  // positions  8 -> 15
        int mb2_eta_HQ = ( dataWordSub >> 32 ) & 0xFF;  // positions 32 -> 39
        int mb1_eta_HQ = ( dataWordSub >> 48 ) & 0xFF;  // positions 48 -> 55

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex
                              << dataWordSub << std::dec << "\t|\t"
                              << "mb1_eta " << mb1_eta  << '\t'
                              << "mb2_eta " << mb2_eta  << '\t'
                              << "mb3_eta " << mb3_eta  << '\t'
                              << "mb1_eta_HQ " <<  mb1_eta_HQ   << '\t'
                              << "mb2_eta_HQ " <<  mb2_eta_HQ   << '\t'
                              << "mb3_eta_HQ " <<  mb3_eta_HQ   << '\n';
        
        //MB1
        posBTI[0] = (mb1_eta & 0x01);
        posBTI[1] = ((mb1_eta & 0x02)>>1);
        posBTI[2] = ((mb1_eta & 0x04)>>2);
        posBTI[3] = ((mb1_eta & 0x08)>>3);
        posBTI[4] = ((mb1_eta & 0x10)>>4);
        posBTI[5] = ((mb1_eta & 0x20)>>5);
        posBTI[6] = (((mb1_eta & 0x40)>>6) || ((mb1_eta & 0x80)>>7));

        qualBTI[0] = (mb1_eta_HQ & 0x01);
        qualBTI[1] = ((mb1_eta_HQ & 0x02)>>1);
        qualBTI[2] = ((mb1_eta_HQ & 0x04)>>2);
        qualBTI[3] = ((mb1_eta_HQ & 0x08)>>3);
        qualBTI[4] = ((mb1_eta_HQ & 0x10)>>4);
        qualBTI[5] = ((mb1_eta_HQ & 0x20)>>5);
        qualBTI[6] = (((mb1_eta_HQ & 0x40)>>6) || ((mb1_eta_HQ & 0x80)>>7));

        theSegments.push_back( L1MuDTChambThDigi( bxNr, wheel, sector-1, 1, posBTI, qualBTI) );
        
        //MB2
        posBTI[0] = (mb2_eta & 0x01);
        posBTI[1] = ((mb2_eta & 0x02)>>1);
        posBTI[2] = ((mb2_eta & 0x04)>>2);
        posBTI[3] = ((mb2_eta & 0x08)>>3);
        posBTI[4] = ((mb2_eta & 0x10)>>4);
        posBTI[5] = ((mb2_eta & 0x20)>>5);
        posBTI[6] = (((mb2_eta & 0x40)>>6) || ((mb2_eta & 0x80)>>7));

        qualBTI[0] = (mb2_eta_HQ & 0x01);
        qualBTI[1] = ((mb2_eta_HQ & 0x02)>>1);
        qualBTI[2] = ((mb2_eta_HQ & 0x04)>>2);
        qualBTI[3] = ((mb2_eta_HQ & 0x08)>>3);
        qualBTI[4] = ((mb2_eta_HQ & 0x10)>>4);
        qualBTI[5] = ((mb2_eta_HQ & 0x20)>>5);
        qualBTI[6] = (((mb2_eta_HQ & 0x40)>>6) || ((mb2_eta_HQ & 0x80)>>7));

        theSegments.push_back( L1MuDTChambThDigi( bxNr, wheel, sector-1, 2, posBTI, qualBTI) );
        
        //MB3
        posBTI[0] = (mb3_eta & 0x01);
        posBTI[1] = ((mb3_eta & 0x02)>>1);
        posBTI[2] = ((mb3_eta & 0x04)>>2);
        posBTI[3] = ((mb3_eta & 0x08)>>3);
        posBTI[4] = ((mb3_eta & 0x10)>>4);
        posBTI[5] = ((mb3_eta & 0x20)>>5);
        posBTI[6] = (((mb3_eta & 0x40)>>6) || ((mb3_eta & 0x80)>>7));

        qualBTI[0] = (mb3_eta_HQ & 0x01);
        qualBTI[1] = ((mb3_eta_HQ & 0x02)>>1);
        qualBTI[2] = ((mb3_eta_HQ & 0x04)>>2);
        qualBTI[3] = ((mb3_eta_HQ & 0x08)>>3);
        qualBTI[4] = ((mb3_eta_HQ & 0x10)>>4);
        qualBTI[5] = ((mb3_eta_HQ & 0x20)>>5);
        qualBTI[6] = (((mb3_eta_HQ & 0x40)>>6) || ((mb3_eta_HQ & 0x80)>>7));

        theSegments.push_back( L1MuDTChambThDigi( bxNr, wheel, sector-1, 3, posBTI, qualBTI) );
        
      }//etha word
     
      else if ( selector == 0xB ) { //MB1/2 output word
 
        int mb2_phi =    ( dataWordSub & 0xFFF);        // positions  0 -> 11
        int mb2_phib =   ( dataWordSub >> 12 ) & 0x3FF; // positions 12 -> 21
        int mb2_qual =   ( dataWordSub >> 22 ) & 0x7;   // positions 22 -> 24
        int mb2_q3 =     ( dataWordSub >> 25 ) & 0x1;   // positions 25
        int mb2_q4 =     ( dataWordSub >> 26 ) & 0x1;   // positions 26
        int mb2_ts2tag = ( dataWordSub >> 28 ) & 0x1;   // positions 28
        /*** NOT UNPACKED  
        int mb2_parity = ( dataWordSub >> 29) & 0x1;    // positions 29
        ***/

        int mb1_phi  =   ( dataWordSub >> 30 ) & 0xFFF; // positions 30 -> 41
        int mb1_phib =   ( dataWordSub >> 42 ) & 0x3FF; // positions 42 -> 51
        int mb1_qual =   ( dataWordSub >> 52 ) & 0x7;   // positions 52 -> 54
        int mb1_q3 =     ( dataWordSub >> 55 ) & 0x1;   // positions 55
        int mb1_q4 =     ( dataWordSub >> 56 ) & 0x1;   // positions 56
        int mb1_ts2tag = ( dataWordSub >> 58 ) & 0x1;   // positions 58
        /*** NOT UNPACKED  
        int mb1_parity = ( dataWordSub >> 59 ) &0x1;    // positions 59
        ***/

        int mb1_phi_conv  = radAngConversion(mb1_phi);
        int mb1_phib_conv = benAngConversion(mb1_phib);

        int mb2_phi_conv  = radAngConversion(mb2_phi);
        int mb2_phib_conv = benAngConversion(mb2_phib);

        phioutSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                1, mb1_phi_conv, mb1_phib_conv, 
                                                mb1_qual, mb1_ts2tag, bxCounter, mb1_q3 + 2*mb1_q4 ) );
        phioutSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                2, mb2_phi_conv, mb2_phib_conv, 
                                                mb2_qual, mb2_ts2tag, bxCounter, mb2_q3 + 2*mb2_q4 ) );

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t"<< std::hex 
                              << dataWordSub   << std::dec      << "\t|\t"
                              << "mb1_ts2tag_out " << mb1_ts2tag    << '\t'
                              << "mb1_qual_out "   << mb1_qual      << '\t'
                              << "mb1_q3_out "     << mb1_q3        << '\t'
                              << "mb1_q4_out "     << mb1_q4        << '\t'
                              << "mb1_phib_out "   << mb1_phib_conv << '\t'
                              << "mb1_phi_out "    << mb1_phi_conv  << '\t'
                              << "mb2_ts2tag_out " << mb2_ts2tag    << '\t'
                              << "mb2_qual_out "   << mb2_qual      << '\t'
                              << "mb2_q3_out "     << mb2_q3        << '\t'
                              << "mb2_q4_out "     << mb2_q4        << '\t'
                              << "mb2_phib_out "   << mb2_phib_conv << '\t'
                              << "mb2_phi_out "    << mb2_phi_conv  << '\n';

      }//MB1/2 output word
 
 
      else if ( selector == 0xC ) { //MB3/4 output word
 
        int mb4_phi =    ( dataWordSub & 0xFFF);        // positions  0 -> 11
        int mb4_phib =   ( dataWordSub >> 12 ) & 0x3FF; // positions 12 -> 21
        int mb4_qual =   ( dataWordSub >> 22 ) & 0x7;   // positions 22 -> 24
        int mb4_q3 =     ( dataWordSub >> 25 ) & 0x1;   // positions 25
        int mb4_q4 =     ( dataWordSub >> 26 ) & 0x1;   // positions 26
        int mb4_ts2tag = ( dataWordSub >> 28 ) & 0x1;   // positions 28
        /*** NOT UNPACKED  
        int mb4_parity = ( dataWordSub >> 29) & 0x1;    // positions 29
        ***/

        int mb3_phi  =   ( dataWordSub >> 30 ) & 0xFFF; // positions 30 -> 41
        int mb3_phib =   ( dataWordSub >> 42 ) & 0x3FF; // positions 42 -> 51
        int mb3_qual =   ( dataWordSub >> 52 ) & 0x7;   // positions 52 -> 54
        int mb3_q3 =     ( dataWordSub >> 55 ) & 0x1;   // positions 55
        int mb3_q4 =     ( dataWordSub >> 56 ) & 0x1;   // positions 56
        int mb3_ts2tag = ( dataWordSub >> 58 ) & 0x1;   // positions 58
        /*** NOT UNPACKED  
        int mb3_parity = ( dataWordSub >> 59 ) &0x1;    // positions 59
        ***/

        int mb3_phi_conv  = radAngConversion(mb3_phi);
        int mb3_phib_conv = benAngConversion(mb3_phib);

        int mb4_phi_conv  = radAngConversion(mb4_phi);
        int mb4_phib_conv = benAngConversion(mb4_phib);

        phioutSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                3, mb3_phi_conv, mb3_phib_conv, 
                                                mb3_qual, mb3_ts2tag, bxCounter, mb3_q3 + 2*mb3_q4 ) );
        phioutSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                4, mb4_phi_conv, mb4_phib_conv, 
                                                mb4_qual, mb4_ts2tag, bxCounter, mb4_q3 + 2*mb4_q4 ) );

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t"<< std::hex 
                              << dataWordSub   << std::dec      << "\t|\t"
                              << "mb3_ts2tag_out " << mb3_ts2tag    << '\t'
                              << "mb3_qual_out "   << mb3_qual      << '\t'
                              << "mb3_q3_out "     << mb3_q3        << '\t'
                              << "mb3_q4_out "     << mb3_q4        << '\t'
                              << "mb3_phib_out "   << mb3_phib_conv << '\t'
                              << "mb3_phi_out "    << mb3_phi_conv  << '\t'
                              << "mb4_ts2tag_out " << mb4_ts2tag    << '\t'
                              << "mb4_qual_out "   << mb4_qual      << '\t'
                              << "mb4_q3_out "     << mb4_q3        << '\t'
                              << "mb4_q4_out "     << mb4_q4        << '\t'
                              << "mb4_phib_out "   << mb4_phib_conv << '\t'
                              << "mb4_phi_out "    << mb4_phi_conv  << '\n';

      }//MB3/4 output word
 
      else if ( selector == 0xD ) { //etha output word
         
        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t ETHA OUTPUT WORD\n";
    
      }//etha output word
       
      else if ( selector == 0x9 || selector == 0xE ) { //RPC word
          
        edm::LogInfo("TwinMux_unpacker") << "RPC WORD [" << std::dec << tm7eventsize << "] : "
                                         << std::hex << dataWordSub << std::dec
                                         << " it pos " << int(DTTM7iterator - DTTM7itend);
  
        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t RPC WORD\n";          

      }//RPC word
 
      else if ( selector == 0x6 ) { //HO word
          
        edm::LogInfo("TwinMux_unpacker") << "HO WORD [" << std::dec << tm7eventsize << "] : "
                                         << std::hex << dataWordSub << std::dec
                                         << " it pos " << int(DTTM7iterator - DTTM7itend);
  
        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t HO WORD\n";          

	int GoodDataFlag_L6  = ( dataWordSub >> 28 ) & 0x1;  
	int GoodDataFlag_L7  = ( dataWordSub >> 57 ) & 0x1;  
	if(!GoodDataFlag_L6 && !GoodDataFlag_L7 ) continue;
	
	int bx_evt=-99;
	if( tm7eventsize == 75) bx_evt =-1;
	if( tm7eventsize == 76) bx_evt =0;
	if( tm7eventsize == 77) bx_evt = 1;

	if(dodebug) 	std::cout <<"wheel, sector : "<< wheel <<", "<< sector<<", "<< bx_evt<<std::endl;
	
	long int mipl6 = -999;
	long int mipl7 = -999;
	
	std::vector<long int > mipl6_col;
	std::vector<long int > mipl7_col;
	
	int ncol_l6= -99, ncol_l7=-99;
	
	// numbers here refers to the bit-locatioin in Twinmux_payload
	if(wheel ==2 && wheel >0 ) {ncol_l6 = 14, ncol_l7 = 43;}  //// link1:0-14, link2:29-43 -> wheel:2/-2
	if(wheel ==1 && wheel >0 ) {ncol_l6 = 17, ncol_l7 = 46;}  //// link1:0-17, link2:29-46 -> wheel:1/-1   
	if(wheel ==0)              {ncol_l6 = 23, ncol_l7 = 52;}  //// link1:0-23, link2:29-52 -> wheel:0      
	if(wheel ==-1 ) {ncol_l6 = 17, ncol_l7 = 46;} 
	if(wheel ==-2 ) {ncol_l6 = 14, ncol_l7 = 43;}
	
	
	//Link-II
	for(int i =0; i <= ncol_l6; i++) { 
	  mipl6 =   ( dataWordSub >> i ) & 0x1 ;    
	  mipl6_col.push_back(mipl6);
	}
	
	//Link-I
	//	std::cout<<"l2"<<std::endl;
	for(int i =29; i <= ncol_l7; i++) { // bit for link-2 is always starting with bit:29
	  mipl7 =   ( dataWordSub >> i ) & 0x1 ;    
	  mipl7_col.push_back(mipl7);
	}
	
	
	if(dodebug) {	  std::cout<<"i:"; for (int i : mipl6_col) {std::cout <<i;} std::cout<<""<<std::endl; for (int i : mipl7_col) std::cout <<i; std::cout<<std::endl;
	  for (int i : mipl6_col) { if(i >0) std::cout<<"YES"<<std::endl;} for (int i : mipl7_col)  { if(i >0) std::cout<<"YES"<<std::endl;} 
	}
	
	HOUnrolledTP unrolled[2160]; //nChan
	int mipFromHcal = -99;
	
	std::vector<HOEmap>::iterator it;
	int count = 0;
	
	//loop over ieta
	if(dodebug)  std::cout<<"In FindMip(), HCAL FED TP size: "<< hoEtaCol.size()<<std::endl;
	
	if(dodebug) std::cout<<"HCAL HO TPs"<<std::endl;
	for (unsigned i=0; i<hoEtaCol.size(); i++) {
	  if(dodebug && hodatabitCol[i] >0)    std::cout<<"("<<hoEtaCol[i]<<","<< hoPhiCol[i]<<"), "<<hosoiCol[i]<< ","<<hodatabitCol[i]<< std::endl;
	}
	
	  for( it = hoemap.begin( ); it != hoemap.end( ); ++it ) {
	    if((*it).iWheel == wheel && (*it).iSector == sector) { 
	      unrolled[count].checked  = true;
	      unrolled[count].ieta     = (*it).iEta;
	      unrolled[count].iphi     = (*it).iPhi;
	      unrolled[count].bx       = bx_evt;
	      if((*it).iLink ==6 && GoodDataFlag_L6)  {unrolled[count].mip = mipl6_col[(*it).iBitloc];} 
	      else if((*it).iLink ==7 && GoodDataFlag_L7) {unrolled[count].mip = mipl7_col[(*it).iBitloc -29];} 
	      else unrolled[count].mip = -99;
	      unrolled[count].wheel    = (*it).iWheel;
	      unrolled[count].sector   = (*it).iSector;
	      unrolled[count].index    = (*it).iBitloc;
	      unrolled[count].link     = (*it).iLink;
	      mipFromHcal =  FindMipFromHcalFeds(&(*it).iEta, &(*it).iPhi);
	      if((*it).iLink ==6) {
		unrolled[count].valid    = ((mipFromHcal ==  mipl6_col[(*it).iBitloc]) ? true : false); 
	      } 
	      else if((*it).iLink ==7) {
		unrolled[count].valid    = ((mipFromHcal ==  mipl7_col[(*it).iBitloc-29]) ? true : false); 
	      }
	      if(dodebug) std::cout <<(*it).iWheel <<"\t"<<(*it).iSector<<"\t"<< (*it).iEta<<"\t"<<(*it).iPhi<<"\t"<<(*it).iLink<<"\t"<<(*it).iBitloc<<"\t"<<unrolled[count].mip<<"\t"<<mipFromHcal<<"\t"<<unrolled[count].valid<<std::endl;	
	      
	      count++;
	    }
	    
	  }
	  if(dodebug)	std::cout<<"count is :"<< count<< std::endl;
	  
	  
	  //filling the classes
	  for (int i=0; i< count; i++) {
	    if (unrolled[i].mip >0) {
	      
	      colls.tphoCont->push_back(HOTPDigiTwinMux(unrolled[i].ieta,
							unrolled[i].iphi,
							unrolled[i].bx,
							unrolled[i].mip,
							unrolled[i].valid,
							unrolled[i].wheel,
						      unrolled[i].sector,
							unrolled[i].index,
							unrolled[i].link));
	      if(dodebug) std::cout<<"******** filled::(eta, phi, mip, valid):"<<unrolled[i].ieta<<"\t"<<unrolled[i].iphi<<"\t"<<unrolled[i].mip<<"\t"<<unrolled[i].valid<<std::endl;
	      if(dodebug) { if(unrolled[i].valid ==0) std::cout<<"Error!"<< std::endl;}
	    }
	  }
	  

	  
      }//HO word
 
      else if ( selector == 0xF ) { //ERROR word

        edm::LogInfo("TwinMux_unpacker") << "ERROR WORD [" << std::dec << tm7eventsize << "] : "
                                         << std::hex << dataWordSub << std::dec
                                         << " it pos " << int(DTTM7iterator - DTTM7itend);
  
        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t ERROR WORD\n";
      }//ERROR word

      else { //unkown word

        edm::LogInfo("TwinMux_unpacker") << "UNKNOWN WORD received " << std::hex << dataWordSub 
                                           << " in FED " << std::hex << TM7fedId;

   	    if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t UNKNOWN WORD\n";
      }
  
      if( DTTM7iterator == DTTM7itend ) break;
      
    } //end of loop over AMCsize


    /// Trailer of payload with CRC
    ++DTTM7iterator;

    if( DTTM7iterator == DTTM7itend ) break;

  } // end for-loop container content

  return;
}



void L1TTwinMuxRawToDigi::calcCRC( long word, int & myC ) {

  int myCRC[16], D[64], C[16];

  for ( int i = 0; i < 64; ++i ) { D[i]    = (word >> i) & 0x1; }
  for ( int i = 0; i < 16; ++i ) { C[i]    = (myC>>i)  & 0x1; }

  myCRC[0] = ( D[63] + D[62] + D[61] + D[60] + D[55] + D[54] +
               D[53] + D[52] + D[51] + D[50] + D[49] + D[48] +
               D[47] + D[46] + D[45] + D[43] + D[41] + D[40] +
               D[39] + D[38] + D[37] + D[36] + D[35] + D[34] +
               D[33] + D[32] + D[31] + D[30] + D[27] + D[26] +
               D[25] + D[24] + D[23] + D[22] + D[21] + D[20] +
               D[19] + D[18] + D[17] + D[16] + D[15] + D[13] +
               D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  +
               D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
               D[0]  + C[0]  + C[1]  + C[2]  + C[3]  + C[4]  +
               C[5]  + C[6]  + C[7]  + C[12] + C[13] + C[14] +
               C[15] )%2;

  myCRC[1] = ( D[63] + D[62] + D[61] + D[56] + D[55] + D[54] +
	           D[53] + D[52] + D[51] + D[50] + D[49] + D[48] +
	           D[47] + D[46] + D[44] + D[42] + D[41] + D[40] +
	           D[39] + D[38] + D[37] + D[36] + D[35] + D[34] +
    	       D[33] + D[32] + D[31] + D[28] + D[27] + D[26] +
    	       D[25] + D[24] + D[23] + D[22] + D[21] + D[20] +
    	       D[19] + D[18] + D[17] + D[16] + D[14] + D[13] +
    	       D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  +
	           D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
	           C[0]  + C[1]  + C[2]  + C[3]  + C[4]  + C[5]  +
	           C[6]  + C[7]  + C[8]  + C[13] + C[14] + C[15] )%2;

  myCRC[2] = ( D[61] + D[60] + D[57] + D[56] + D[46] + D[42] +
	           D[31] + D[30] + D[29] + D[28] + D[16] + D[14] +
	           D[1]  + D[0]  + C[8]  + C[9]  + C[12] + C[13] )%2;

  myCRC[3] = ( D[62] + D[61] + D[58] + D[57] + D[47] + D[43] +
	           D[32] + D[31] + D[30] + D[29] + D[17] + D[15] +
	           D[2]  + D[1]  + C[9]  + C[10] + C[13] + C[14] )%2;

  myCRC[4] = ( D[63] + D[62] + D[59] + D[58] + D[48] + D[44] +
    	       D[33] + D[32] + D[31] + D[30] + D[18] + D[16] + 
	           D[3]  + D[2]  + C[0]  + C[10] + C[11] + C[14] +
	           C[15] )%2;

  myCRC[5] = ( D[63] + D[60] + D[59] + D[49] + D[45] + D[34] +
	           D[33] + D[32] + D[31] + D[19] + D[17] + D[4]  +
    	       D[3]  + C[1]  + C[11] + C[12] + C[15] )%2;

  myCRC[6] = ( D[61] + D[60] + D[50] + D[46] + D[35] + D[34] +
	           D[33] + D[32] + D[20] + D[18] + D[5]  + D[4]  +
	           C[2]  + C[12] + C[13] )%2;

  myCRC[7] = ( D[62] + D[61] + D[51] + D[47] + D[36] + D[35] +
    	       D[34] + D[33] + D[21] + D[19] + D[6]  + D[5]  +
	           C[3]  + C[13] + C[14] )%2;

  myCRC[8] = ( D[63] + D[62] + D[52] + D[48] + D[37] + D[36] +
	           D[35] + D[34] + D[22] + D[20] + D[7]  + D[6]  +
    	       C[0]  + C[4]  + C[14] + C[15] )%2;

  myCRC[9] = ( D[63] + D[53] + D[49] + D[38] + D[37] + D[36] +
	           D[35] + D[23] + D[21] + D[8]  + D[7]  + C[1]  +
	           C[5]  + C[15] )%2;

  myCRC[10] = ( D[54] + D[50] + D[39] + D[38] + D[37] + D[36] + 
       		    D[24] + D[22] + D[9]  + D[8]  + C[2]  + C[6] )%2;

  myCRC[11] = ( D[55] + D[51] + D[40] + D[39] + D[38] + D[37] +
		        D[25] + D[23] + D[10] + D[9]  + C[3]  + C[7] )%2;

  myCRC[12] = ( D[56] + D[52] + D[41] + D[40] + D[39] + D[38] +
        		D[26] + D[24] + D[11] + D[10] + C[4]  + C[8] )%2;

  myCRC[13] = ( D[57] + D[53] + D[42] + D[41] + D[40] + D[39] +
		        D[27] + D[25] + D[12] + D[11] + C[5]  + C[9] )%2;

  myCRC[14] = ( D[58] + D[54] + D[43] + D[42] + D[41] + D[40] +
        		D[28] + D[26] + D[13] + D[12] + C[6]  + C[10] )%2;

  myCRC[15] = ( D[63] + D[62] + D[61] + D[60] + D[59] + D[54] +
		        D[53] + D[52] + D[51] + D[50] + D[49] + D[48] + 
	        	D[47] + D[46] + D[45] + D[44] + D[42] + D[40] +
        		D[39] + D[38] + D[37] + D[36] + D[35] + D[34] + 
		        D[33] + D[32] + D[31] + D[30] + D[29] + D[26] +
        		D[25] + D[24] + D[23] + D[22] + D[21] + D[20] + 
        		D[19] + D[18] + D[17] + D[16] + D[15] + D[14] +
        		D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  + 
        		D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
        		D[0]  + C[0]  + C[1]  + C[2]  + C[3]  + C[4]  + 
	        	C[5]  + C[6]  + C[11] + C[12] + C[13] + C[14] +
	        	C[15] )%2;

  int tempC = 0x0;  
  for ( int i = 0; i < 16 ; ++i) { tempC = tempC + ( myCRC[i] << i ); }
  myC = tempC;
  return;
}


L1TTwinMuxRawToDigi::Collections::Collections() {
  tphoCont=0;
}


//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TTwinMuxRawToDigi);
