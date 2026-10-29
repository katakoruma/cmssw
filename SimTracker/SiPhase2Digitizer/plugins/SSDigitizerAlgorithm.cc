#include <algorithm>
#include <array>
#include <iostream>
#include <cmath>
#include <numeric>

#include "SimTracker/SiPhase2Digitizer/plugins/SSDigitizerAlgorithm.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Geometry
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include "CondFormats/DataRecord/interface/SiPhase2OuterTrackerCondDataRecords.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

using namespace edm;

namespace {
  double nFactorial(int n);
  double aScalingConstant(int N, int i);
}  // namespace
namespace {
  double nFactorial(int n) { return std::tgamma(n + 1); }
  double aScalingConstant(int N, int i) {
    return std::pow(-1, (double)i) * nFactorial(N) * nFactorial(N + 2) /
           (nFactorial(N - i) * nFactorial(N + 2 - i) * nFactorial(i));
  }
  // Standard normal deviate above s (as GaussianTailNoiseGenerator: rejection for s < 1, else
  // Marsaglia's tail method, Knuth v2, 3rd ed, pp 123-128)
  double gaussianTail(double s, CLHEP::HepRandomEngine* engine) {
    if (s < 1.) {
      double x;
      do {
        x = CLHEP::RandGaussQ::shoot(engine, 0., 1.);
      } while (x < s);
      return x;
    }
    double u, v, x;
    do {
      u = CLHEP::RandFlat::shoot(engine);
      do {
        v = CLHEP::RandFlat::shoot(engine);
      } while (v == 0.);
      x = std::sqrt(s * s - 2 * std::log(v));
    } while (x * u > s);
    return x;
  }
  // Binomial(n, p) deviate by inversion of the cumulative distribution, exact also for very small p
  // (CLHEP::RandBinomial returns 0 for p <= 1e-5); on average 1 + n p steps
  int binomial(int n, double p, CLHEP::HepRandomEngine* engine) {
    if (p <= 0.)
      return 0;
    if (p >= 1.)
      return n;
    const double u = CLHEP::RandFlat::shoot(engine);
    const double ratio = p / (1. - p);
    double pk = std::exp(n * std::log1p(-p));  // P(0)
    double cdf = pk;
    int k = 0;
    while (u > cdf && k < n) {
      pk *= ratio * (n - k) / (k + 1);
      cdf += pk;
      ++k;
    }
    return k;
  }
}  // namespace

void SSDigitizerAlgorithm::init(const edm::EventSetup& es) {
  if (use_LorentzAngle_DB_)  // Get Lorentz angle from DB record
    siPhase2OTLorentzAngle_ = &es.getData(siPhase2OTLorentzAngleToken_);

  if (use_deadmodule_DB_)  // Get Bad Channel (SiStripBadStrip) from DB
    badChannelPayload_ = &es.getData(badChannelToken_);

  geom_ = &es.getData(geomToken_);
  tTopo_ = &es.getData(topoToken_);
}

void SSDigitizerAlgorithm::initializeEvent(CLHEP::HepRandomEngine& eng) {
  Phase2TrackerDigitizerAlgorithm::initializeEvent(eng);
  commonModeCache_.clear();
}

SSDigitizerAlgorithm::SSDigitizerAlgorithm(const edm::ParameterSet& conf, edm::ConsumesCollector iC)
    : Phase2TrackerDigitizerAlgorithm(conf.getParameter<ParameterSet>("AlgorithmCommon"),
                                      conf.getParameter<ParameterSet>("SSDigitizerAlgorithm"),
                                      iC),
      hitDetectionMode_(conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<int>("HitDetectionMode")),
      pulseShapeParameters_(conf.getParameter<ParameterSet>("SSDigitizerAlgorithm")
                                .getParameter<std::vector<double> >("PulseShapeParameters")),
      deadTime_(conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<double>("CBCDeadTime")),
      geomToken_(iC.esConsumes()),
      commonModeCBCFraction_(
          conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<double>("CommonModeCBCFraction")),
      commonModeHybridFraction_(
          conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<double>("CommonModeHybridFraction")),
      commonModeHybridCorrelation_(
          conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<double>("CommonModeHybridCorrelation")),
      addCommonMode_(commonModeCBCFraction_ > 0. || commonModeHybridFraction_ > 0.),
      independentNoiseFraction_(std::sqrt(std::max(0.,
                                                   1. - commonModeCBCFraction_ * commonModeCBCFraction_ -
                                                       commonModeHybridFraction_ * commonModeHybridFraction_))),
      topoToken_(iC.esConsumes()) {
  if (commonModeCBCFraction_ < 0. || commonModeHybridFraction_ < 0. ||
      commonModeCBCFraction_ * commonModeCBCFraction_ + commonModeHybridFraction_ * commonModeHybridFraction_ >= 1. ||
      std::abs(commonModeHybridCorrelation_) > 1.)
    throw cms::Exception("Configuration")
        << "SSDigitizerAlgorithm: need CommonModeCBCFraction, CommonModeHybridFraction >= 0 with a sum of squares < 1 "
           "and |CommonModeHybridCorrelation| <= 1, got "
        << commonModeCBCFraction_ << ", " << commonModeHybridFraction_ << ", " << commonModeHybridCorrelation_;
  if (addCommonMode_)
    edm::LogInfo("SSDigitizerAlgorithm") << "2S common-mode noise: sigma per hybrid " << commonModeHybridFraction_
                                         << ", per CBC " << commonModeCBCFraction_ << ", per strip "
                                         << independentNoiseFraction_ << " of the strip noise; hybrid correlation "
                                         << commonModeHybridCorrelation_;
  if (use_LorentzAngle_DB_)
    siPhase2OTLorentzAngleToken_ = iC.esConsumes();

  if (use_deadmodule_DB_) {
    std::string badChannelLabel_ = conf.getParameter<ParameterSet>("SSDigitizerAlgorithm")
                                       .getUntrackedParameter<std::string>("BadChannelLabel", "");
    badChannelToken_ = iC.esConsumes(edm::ESInputTag{"", badChannelLabel_});
  }

  pixelFlag_ = false;
  LogDebug("SSDigitizerAlgorithm ") << "SSDigitizerAlgorithm constructed "
                                    << "Configuration parameters:"
                                    << "Threshold/Gain = "
                                    << "threshold in electron Endcap = " << theThresholdInE_Endcap_
                                    << "threshold in electron Barrel = " << theThresholdInE_Barrel_ << " "
                                    << theElectronPerADC_ << " " << theAdcFullScale_ << " The delta cut-off is set to "
                                    << tMax_ << " channel inefficiency flag" << addChannelInefficiency_;
  storeSignalShape();
}
SSDigitizerAlgorithm::~SSDigitizerAlgorithm() { LogDebug("SSDigitizerAlgorithm") << "SSDigitizerAlgorithm deleted"; }
//
// -- Select the Hit for Digitization
//
bool SSDigitizerAlgorithm::select_hit(const PSimHit& hit, double tCorr, double& sigScale) const {
  bool result = false;
  if (hitDetectionMode_ == SSDigitizerAlgorithm::SampledMode)
    result = select_hit_sampledMode(hit, tCorr, sigScale);
  else if (hitDetectionMode_ == SSDigitizerAlgorithm::LatchedMode)
    result = select_hit_latchedMode(hit, tCorr, sigScale);
  else {
    double toa = hit.tof() - tCorr;
    result = (toa > theTofLowerCut_ && toa < theTofUpperCut_);
  }
  return result;
}
//
// -- Select Hits in Sampled Mode
//
bool SSDigitizerAlgorithm::select_hit_sampledMode(const PSimHit& hit, double tCorr, double& sigScale) const {
  double toa = hit.tof() - tCorr;
  double sampling_time = bx_time;

  DetId det_id = DetId(hit.detUnitId());
  float theThresholdInE =
      (det_id.subdetId() == StripSubdetector::TOB) ? theThresholdInE_Barrel_ : theThresholdInE_Endcap_;

  sigScale = getSignalScale(sampling_time - toa);
  return (sigScale * hit.energyLoss() / GeVperElectron_ > theThresholdInE);
}
//
// -- Select Hits in Hit Detection Mode
//
bool SSDigitizerAlgorithm::select_hit_latchedMode(const PSimHit& hit, double tCorr, double& sigScale) const {
  float toa = hit.tof() - tCorr;
  toa -= hit.eventId().bunchCrossing() * bx_time;

  float sampling_time = (-1) * (hit.eventId().bunchCrossing() + 1) * bx_time;

  DetId det_id = DetId(hit.detUnitId());
  float theThresholdInE =
      (det_id.subdetId() == StripSubdetector::TOB) ? theThresholdInE_Barrel_ : theThresholdInE_Endcap_;

  bool lastPulse = true;
  bool aboveThr = false;
  for (float i = deadTime_; i <= bx_time; i++) {
    sigScale = getSignalScale(sampling_time - toa + i);

    aboveThr = (sigScale * hit.energyLoss() / GeVperElectron_ > theThresholdInE);
    if (!lastPulse && aboveThr)
      return true;

    lastPulse = aboveThr;
  }
  return false;
}
double SSDigitizerAlgorithm::cbc3PulsePolarExpansion(double x) const {
  constexpr size_t max_par = 6;
  if (pulseShapeParameters_.size() < max_par)
    return -1;
  double xOffset = pulseShapeParameters_[0];
  double tau = pulseShapeParameters_[1];
  double r = pulseShapeParameters_[2];
  double theta = pulseShapeParameters_[3];
  int nTerms = static_cast<int>(pulseShapeParameters_[4]);

  double fN = 0;
  double xx = x - xOffset;
  if (xx < 0)
    return 0;

  for (int i = 0; i < nTerms; i++) {
    double angularTerm = 0;
    double temporalTerm = 0;
    double rTerm = std::pow(r, i) / (std::pow(tau, 2. * i) * ::nFactorial(i + 2));
    for (int j = 0; j <= i; j++) {
      angularTerm += std::pow(std::cos(theta), (double)(i - j)) * std::pow(std::sin(theta), (double)j);
      temporalTerm += ::aScalingConstant(i, j) * std::pow(xx, (double)(i - j)) * std::pow(tau, (double)j);
    }
    double fi = rTerm * angularTerm * temporalTerm;

    fN += fi;
  }
  return fN;
}
double SSDigitizerAlgorithm::signalShape(double x) const {
  double xOffset = pulseShapeParameters_[0];
  double tau = pulseShapeParameters_[1];
  double maxCharge = pulseShapeParameters_[5];

  double xx = x - xOffset;
  return maxCharge * (std::exp(-xx / tau) * std::pow(xx / tau, 2.) * cbc3PulsePolarExpansion(x));
}
void SSDigitizerAlgorithm::storeSignalShape() {
  for (size_t i = 0; i < interpolationPoints; i++) {
    float val = i / interpolationStep;

    pulseShapeVec_.emplace_back(signalShape(val));
  }
}
double SSDigitizerAlgorithm::getSignalScale(double xval) const {
  double res = 0.0;
  int len = pulseShapeVec_.size();

  if (xval < 0.0 || xval * interpolationStep >= len)
    return res;

  unsigned int lower = std::floor(xval) * interpolationStep;
  unsigned int upper = std::ceil(xval) * interpolationStep;
  for (size_t i = lower + 1; i < upper * interpolationStep; i++) {
    float val = i * 0.1;
    if (val > xval) {
      res = pulseShapeVec_[i - 1];
      break;
    }
  }
  return res;
}
//
// -- Compare Signal with Threshold
//
bool SSDigitizerAlgorithm::isAboveThreshold(const digitizerUtility::SimHitInfo* const hisInfo,
                                            float charge,
                                            float thr) const {
  return (charge >= thr);
}
//
// -- Read Bad Channels from the Condidion DB and kill channels/module accordingly
//
void SSDigitizerAlgorithm::module_killing_DB(const Phase2TrackerGeomDetUnit* ph2det) {
  auto detId = ph2det->geographicalId().rawId();
  auto& theSignal = _signal[detId];  // Caller ensures detId exists

  SiStripBadStrip::Range range = badChannelPayload_->getRange(detId);
  for (std::vector<unsigned int>::const_iterator badChannel = range.first; badChannel != range.second; ++badChannel) {
    const auto& firstStrip = badChannelPayload_->decodePhase2(*badChannel).firstStrip;
    const auto& channelRange = badChannelPayload_->decodePhase2(*badChannel).range;

    for (int index = 0; index < channelRange; index++) {
      for (auto& [channel, sig_data] : theSignal) {
        if (channel == firstStrip + index)
          sig_data.set(0.);
      }
    }
  }
}
//
// -- Common-mode noise
//
// The noise of a strip is split as n = A_h + B_hk + e: A_h is common to the hybrid h (one per strip column),
// B_hk to the CBC k of that hybrid (127 strips of each sensor) and e is independent per strip, with sigmas
// CommonModeHybridFraction, CommonModeCBCFraction and sqrt(1 - both squared) times NoiseInElectrons, so the
// total noise of a strip, and with it the mean occupancy, is unchanged. A_0 and A_1 of a module are
// correlated by CommonModeHybridCorrelation. Both sensors of the module share A and B.
//
const std::vector<float>& SSDigitizerAlgorithm::commonMode(const Phase2TrackerGeomDetUnit* ph2det) {
  auto [it, inserted] = commonModeCache_.try_emplace(tTopo_->stack(ph2det->geographicalId()));
  auto& cm = it->second;
  if (inserted) {
    const Phase2TrackerTopology* topol = &ph2det->specificTopology();
    const int nHybrids = topol->ncolumns();
    const int nCBC = topol->nrows() / stripsPerCBC;
    const double sigmaHybrid = commonModeHybridFraction_ * theNoiseInElectrons_;
    const double sigmaCBC = commonModeCBCFraction_ * theNoiseInElectrons_;
    const double rho = commonModeHybridCorrelation_;
    const double first = CLHEP::RandGaussQ::shoot(rengine_, 0., 1.);
    cm.resize(nHybrids * nCBC);
    for (int hybrid = 0; hybrid < nHybrids; ++hybrid) {
      // correlation rho of every further hybrid with the first one (a 2S module has two)
      const double z =
          hybrid == 0 ? first : rho * first + std::sqrt(1. - rho * rho) * CLHEP::RandGaussQ::shoot(rengine_, 0., 1.);
      for (int cbc = 0; cbc < nCBC; ++cbc)
        cm[hybrid * nCBC + cbc] = sigmaHybrid * z + sigmaCBC * CLHEP::RandGaussQ::shoot(rengine_, 0., 1.);
    }
  }
  return cm;
}
void SSDigitizerAlgorithm::add_noise(const Phase2TrackerGeomDetUnit* ph2det) {
  if (!addCommonMode_) {
    Phase2TrackerDigitizerAlgorithm::add_noise(ph2det);
    return;
  }
  const auto& cm = commonMode(ph2det);
  const int nrows = ph2det->specificTopology().nrows();
  auto& theSignal = _signal[ph2det->geographicalId().rawId()];  // Caller ensures existence of detID
  for (auto& [chan, sig_data] : theSignal) {
    const auto [row, column] = Phase2TrackerDigi::channelToPixel(chan);
    auto noise = independentNoiseFraction_ * gaussDistribution_->fire() + cm[cbcIndex(row, column, nrows)];
    if ((sig_data.ampl() + noise) < 0.)
      sig_data.set(0);
    else
      sig_data += noise;
  }
}
void SSDigitizerAlgorithm::add_noisy_channels(const Phase2TrackerGeomDetUnit* ph2det, float channelThreshold) {
  if (!addCommonMode_ || theNoiseInElectrons_ <= 0.) {
    Phase2TrackerDigitizerAlgorithm::add_noisy_channels(ph2det, channelThreshold);
    return;
  }
  const auto& cm = commonMode(ph2det);
  const Phase2TrackerTopology* topol = &ph2det->specificTopology();
  const int nrows = topol->nrows();
  const int ncolumns = topol->ncolumns();
  auto& theSignal = _signal[ph2det->geographicalId().rawId()];  // Caller ensures existence of detID
  const double sigma = independentNoiseFraction_ * theNoiseInElectrons_;
  const double threshold = channelThreshold * theNoiseInElectrons_;  // in electrons
  std::array<int, stripsPerCBC> strips;
  for (int column = 0; column < ncolumns; ++column) {
    for (int firstStrip = 0; firstStrip < nrows; firstStrip += stripsPerCBC) {
      const int nStrips = std::min(stripsPerCBC, nrows - firstStrip);
      const float cmValue = cm[cbcIndex(firstStrip, column, nrows)];
      // A strip of this CBC fires if its own noise is above the threshold minus the common mode
      const double cut = (threshold - cmValue) / sigma;
      const double probability = 0.5 * std::erfc(cut / std::sqrt(2.));
      const int nNoisy = binomial(nStrips, probability, rengine_);
      // distinct strips by a partial Fisher-Yates shuffle
      std::iota(strips.begin(), strips.begin() + nStrips, firstStrip);
      for (int i = 0; i < nNoisy; ++i) {
        const int j = i + static_cast<int>(CLHEP::RandFlat::shoot(rengine_, nStrips - i));
        std::swap(strips[i], strips[j]);
        const float amplitude = sigma * gaussianTail(cut, rengine_) + cmValue;
        theSignal.try_emplace(Phase2TrackerDigi::pixelToChannel(strips[i], column), amplitude, nullptr, -1.);
      }
    }
  }
}
