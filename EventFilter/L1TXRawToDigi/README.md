L1TXRawToDigi contains CMSSW EDProducers for L1T modules not using the 
common unpacker framework in L1TRawToDigi.

L1TCaloLayer1RawToDigi unpacks Layer-1 calorimeter trigger data obtained from its three FEDs.

The producer makes ECAL and HCAL TPG objects.

These should be identical to the ones read out from ECAL TCCs
and HCAL uHTRs.  However, any link issues could result in 
discrepancies.  



L1TTriwMuxRawToDigi unpacks the TwinMux in the Muon system.


*** For running the HO TPs code setup:
Use config /test/TMReader_Test_collision.py
It has been tested with SingleMuon RAW dataset of 2017B
The required local emap is in /data directory