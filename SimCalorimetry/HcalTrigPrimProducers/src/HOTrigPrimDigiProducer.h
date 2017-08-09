#ifndef HOTrigPrimDigiProducer_h
#define HOTrigPrimDigiProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimCalorimetry/HcalTrigPrimAlgos/interface/HcalTriggerPrimitiveAlgo.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HODataFrame.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include <vector>

class HOTrigPrimDigiProducer : public edm::EDProducer {
public:

  explicit HOTrigPrimDigiProducer(const edm::ParameterSet& ps);
  virtual ~HOTrigPrimDigiProducer() {}

  /**Produces the EDM products,*/
  virtual void produce(edm::Event& e, const edm::EventSetup& c);

private:

  HcalTriggerPrimitiveAlgo theAlgo_;

  /// input tags for HCAL digis
  std::vector<edm::InputTag> inputLabel_;
  // this seems a strange way of doing things
  edm::EDGetTokenT<HODigiCollection> tok_ho_;

  /// input tag for FEDRawDataCollection
  edm::InputTag inputTagFEDRaw_;
  edm::EDGetTokenT<FEDRawDataCollection> tok_raw_;

  bool runFrontEndFormatError_;
};

#endif

