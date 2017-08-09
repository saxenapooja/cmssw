#include "SimCalorimetry/HcalTrigPrimProducers/src/HOTrigPrimDigiProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/HcalDigi/interface/HOTriggerPrimitiveDigi.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/HcalTPGCompressor.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "CondFormats/HcalObjects/interface/HcalLutMetadata.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include <algorithm>

HOTrigPrimDigiProducer::HOTrigPrimDigiProducer(const edm::ParameterSet& ps)
  : 
  theAlgo_(ps.getParameter<std::vector<double> >("weights"),
	   ps.getParameter<int>("latency"),
	   ps.getParameter<int>("numberOfSamples"),
	   ps.getParameter<int>("numberOfPresamples"),
	   ps.getParameter<unsigned int>("muonBitThreshold_high"),
	   ps.getParameter<unsigned int>("muonBitThreshold_low")
	   ),
  inputLabel_(ps.getParameter<std::vector<edm::InputTag> >("inputLabel")),
  inputTagFEDRaw_(ps.getParameter<edm::InputTag> ("InputTagFEDRaw")),
  runFrontEndFormatError_(ps.getParameter<bool>("FrontEndFormatError"))
{

  // register for data access
  if (runFrontEndFormatError_) {
    tok_raw_ = consumes<FEDRawDataCollection>(inputTagFEDRaw_);
  }
  tok_ho_ = consumes<HODigiCollection>(inputLabel_[0]);
  
  
  //typedef edm::SortedCollection<HOTriggerPrimitiveDigi> HOTrigPrimDigiCollection;
  produces<HOTrigPrimDigiCollection>();
  theAlgo_.setPeakFinderAlgorithm(ps.getParameter<int>("PeakFinderAlgorithm"));
}


void HOTrigPrimDigiProducer::produce(edm::Event& iEvent, const edm::EventSetup& eventSetup) {
  //  int e = iEvent.id().event();

  //  if(e == 294807358 || e == 290411280){

 // Step A: get the conditions, for the decoding
  edm::ESHandle<HcalTPGCoder> inputCoder;
  eventSetup.get<HcalTPGRecord>().get(inputCoder);
  
  edm::ESHandle<CaloTPGTranscoder> outTranscoder;
  eventSetup.get<CaloTPGRecord>().get(outTranscoder);

  // edm::ESHandle<HcalLutMetadata> lutMetadata;
  // eventSetup.get<HcalLutMetadataRcd>().get(lutMetadata);
  // float rctlsb = lutMetadata->getRctLsb();

  edm::ESHandle<HcalTrigTowerGeometry> pG;
  eventSetup.get<CaloGeometryRecord>().get(pG);
  
  // Step B: Create empty output
  std::unique_ptr<HOTrigPrimDigiCollection> result(new HOTrigPrimDigiCollection());
 
  edm::Handle<HODigiCollection>   hoDigis;
  iEvent.getByToken(tok_ho_,hoDigis);

  // protect here against missing input collections
  // there is no protection in HcalTriggerPrimitiveAlgo
  
  if (!hoDigis.isValid()) {
    edm::LogInfo("HOTrigPrimDigiProducer")
      << "\nWarning: HODigiCollection with input tag "
      << inputLabel_[2]
      << "\nrequested in configuration, but not found in the event."
      << "\nQuit returning empty product." << std::endl;
    
    // put empty HcalTrigPrimDigiCollection in the event
    iEvent.put(std::move(result));
  
    return;
  }
  
  
  theAlgo_.run(inputCoder.product(), 
	       *hoDigis,
	       *result,
	       &(*pG));
  
  iEvent.put(std::move(result));
  //}
}

