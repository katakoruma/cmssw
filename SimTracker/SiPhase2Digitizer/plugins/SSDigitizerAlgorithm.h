#ifndef _SimTracker_SiPhase2Digitizer_SSDigitizerAlgorithm_h
#define _SimTracker_SiPhase2Digitizer_SSDigitizerAlgorithm_h

#include <unordered_map>
#include <vector>

#include "CondFormats/DataRecord/interface/SiPhase2OuterTrackerLorentzAngleRcd.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "SimTracker/SiPhase2Digitizer/plugins/Phase2TrackerDigitizerAlgorithm.h"
#include "CondFormats/SiStripObjects/interface/SiStripBadStrip.h"

class SSDigitizerAlgorithm : public Phase2TrackerDigitizerAlgorithm {
public:
  SSDigitizerAlgorithm(const edm::ParameterSet& conf, edm::ConsumesCollector iC);
  ~SSDigitizerAlgorithm() override;

  // initialization that cannot be done in the constructor
  void init(const edm::EventSetup& es) override;
  void initializeEvent(CLHEP::HepRandomEngine& eng) override;

  bool select_hit(const PSimHit& hit, double tCorr, double& sigScale) const override;
  bool isAboveThreshold(const digitizerUtility::SimHitInfo* hitInfo, float charge, float thr) const override;

protected:
  // With common-mode noise the strip noise is split into a part per hybrid, per CBC and per strip
  void add_noise(const Phase2TrackerGeomDetUnit* ph2det) override;
  void add_noisy_channels(const Phase2TrackerGeomDetUnit* ph2det, float channelThreshold) override;

private:
  enum { SquareWindow, SampledMode, LatchedMode, SampledOrLachedMode, HIPFindingMode };
  double cbc3PulsePolarExpansion(double x) const;
  double signalShape(double x) const;
  double getSignalScale(double xval) const;
  void storeSignalShape();
  bool select_hit_sampledMode(const PSimHit& hit, double tCorr, double& sigScale) const;
  bool select_hit_latchedMode(const PSimHit& hit, double tCorr, double& sigScale) const;
  void module_killing_DB(const Phase2TrackerGeomDetUnit* pixdet) override;

  int hitDetectionMode_;
  std::vector<double> pulseShapeVec_;
  std::vector<double> pulseShapeParameters_;
  float deadTime_;
  edm::ESGetToken<SiPhase2OuterTrackerLorentzAngle, SiPhase2OuterTrackerLorentzAngleSimRcd> siPhase2OTLorentzAngleToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  static constexpr float bx_time{25};
  static constexpr size_t interpolationPoints{1000};
  static constexpr int interpolationStep{10};
  edm::ESGetToken<SiStripBadStrip, SiPhase2OuterTrackerBadStripRcd> badChannelToken_;
  const SiStripBadStrip* badChannelPayload_;

  // Common mode in electrons of every CBC of the module (index column * nCBC + CBC), drawn once per
  // module and event, so that both sensors of a 2S module see the same values
  const std::vector<float>& commonMode(const Phase2TrackerGeomDetUnit* ph2det);
  int cbcIndex(int row, int column, int nrows) const { return column * (nrows / stripsPerCBC) + row / stripsPerCBC; }

  // A CBC reads 127 strips of each sensor; a hybrid (one per strip column) holds 8 CBCs
  static constexpr int stripsPerCBC{127};
  const double commonModeCBCFraction_;        // sigma of the CM per CBC / total strip noise
  const double commonModeHybridFraction_;     // sigma of the CM per hybrid / total strip noise
  const double commonModeHybridCorrelation_;  // correlation of the CM of the two hybrids of a module
  const bool addCommonMode_;
  const double independentNoiseFraction_;  // sigma of the remaining noise per strip / total strip noise
  std::unordered_map<uint32_t, std::vector<float>> commonModeCache_;  // per module (stack DetId), cleared per event
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const TrackerTopology* tTopo_{nullptr};
};
#endif
