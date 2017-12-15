/**
 *  @file   larpandoracontent/LArVertex/TrackDirectionTool.h
 *
 *  @brief  Header file for the candidate vertex creation AlgorithmTool class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_DIRECTION_TOOL_H
#define LAR_TRACK_DIRECTION_TOOL_H 1

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "Pandora/AlgorithmTool.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"

#include <utility>
#include <algorithm>
#include <deque>
#include <unordered_map>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <numeric>

namespace lar_content
{

/**
 *  @brief  TrackDirectionTool::AlgorithmTool class
 */
class TrackDirectionTool : public pandora::AlgorithmTool
{

friend void GetSplitChiSquared(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag);

public:

    /**
     *  @brief  Factory class for instantiating AlgorithmTool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    /**
     *  @brief  Default constructor
     */
    TrackDirectionTool();

    /**
     *  @brief  Default destructor
     */
    ~TrackDirectionTool();

    class HitCharge
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        HitCharge(const pandora::CaloHit* caloHit, float &longitudinalPosition, float &hitWidth, float &hitCharge, float &uncertainty);

        /**
         *  @brief  Constructor
         *
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        HitCharge();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        const pandora::CaloHit* GetCaloHit() const;

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetLongitudinalPosition() const;

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetHitWidth() const;

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetCharge() const;

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetQoverX() const;

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetUncertainty() const;

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        void SetDistanceToNN(float &distance);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetDistanceToNN() const;

        bool                                       m_intails;

    private:
        const pandora::CaloHit*                    m_calohit;
        float                                      m_longitudinalposition;    ///<
        float                                      m_hitwidth;                ///<
        float                                      m_charge;                  ///<
        float                                      m_qoverx;                  ///<
        float                                      m_uncertainty;          ///<
        float                                      m_distancetonearestneighbour;          ///<
    };

    typedef std::vector<HitCharge> HitChargeVector;

    class KinkObject
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        KinkObject(float &xPosition, float &kinkAngle);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetXPosition() const;

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetKinkAngle() const;

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        void SetWLongitudinalPosition(float &rL);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetWLongitudinalPosition() const;

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        void SetNearestHit(HitCharge &hitCharge);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        HitCharge GetNearestHit() const;

    private:
        const float                                      m_xposition;               ///<
        const float                                      m_kinkangle;               ///<
        float                                            m_wrl;
        HitCharge                                        m_hit;
    };

    typedef std::vector<KinkObject> KinkVector;

    class LookupTable
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        LookupTable();

        /**
         *  @brief  Constructor
         *
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        LookupTable(double &initialEnergy, double &binWidth);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        std::map<int, double> GetMap();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        void SetMap(std::map<int, double> &map);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        std::map<double, int> GetReverseMap();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        void SetReverseMap(std::map<double, int> &map);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        double GetInitialEnergy();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        void SetInitialEnergy(double &initialEnergy);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        double GetBinWidth();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        void SetBinWidth(double &binWidth);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        void SetMaxRange(double &maxRange);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        double GetMaxRange();

    private:
        std::map<int, double>                       m_map;
        std::map<double, int>                       m_reversemap;
        double                                      m_binwidth;               ///<
        double                                      m_initialenergy;
        double                                      m_maxrange;
    };

    class FitResult
    {
    public:

        /**
         *  @brief  Constructor
         *
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        FitResult();

        /**
         *  @brief  Constructor
         *
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        FitResult(HitChargeVector &hitChargeVector, int &numberHits, float &meanQoverX, float &forwardsChiSquared, float &backwardsChiSquared);

        /**
         *  @brief  Constructor
         *
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        FitResult(HitChargeVector &hitChargeVector, HitChargeVector &forwardsRecoHits, HitChargeVector &backwardsRecoHits, int &numberHits, float &meanQoverX, float &forwardsChiSquared, float &backwardsChiSquared);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        TrackDirectionTool::HitChargeVector GetHitChargeVector();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        TrackDirectionTool::HitChargeVector GetForwardsRecoHits();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        TrackDirectionTool::HitChargeVector GetBackwardsRecoHits();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetForwardsChiSquared();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetBackwardsChiSquared();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        int GetNHits();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        int GetDirectionEstimate();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetMinChiSquared();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetMinChiSquaredPerHit();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float GetMeanQoverX();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        void SetCoordinates(TwoDSlidingFitResult &slidingFitResult);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        const pandora::CartesianVector GetBeginPoint();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        const pandora::CartesianVector GetEndPoint();

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        float ComputeProbability();

    private:
        HitChargeVector     m_hitchargevector;
        HitChargeVector     m_forwardsrecohits;
        HitChargeVector     m_backwardsrecohits;
        int                 m_nhits;
        float               m_meanqoverx;
        float               m_forwardschisquared;
        float               m_backwardschisquared;
        float               m_probability;
        float               m_begin_x;
        float               m_begin_y;
        float               m_begin_z;
        float               m_end_x;
        float               m_end_y;
        float               m_end_z;
    };

    class JumpObject
    {
        public:

            JumpObject(float &longitudinalPosition, float &jumpValue);

            float GetLongitudinalPosition();

            float GetJumpValue();

        private:
            float       m_longitudinalposition;
            float       m_jumpvalue;

    };

    TrackDirectionTool::FitResult Run(pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pTargetClusterW);

    TrackDirectionTool::FitResult Run(pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject* PFO);

    const pandora::Cluster* GetTargetClusterFromPFO(const pandora::ParticleFlowObject* PFO);

    void WriteLookupTableToTree(LookupTable &lookupTable);

    void ReadLookupTableFromTree(LookupTable &lookupTable);

    void FillHitChargeVector(const pandora::Cluster *const pCluster, HitChargeVector &hitChargeVector);

    void TrackInnerFilter(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector);

    void SetNearestNeighbourValues(HitChargeVector &innerHitChargeVector, int &nNeighboursToConsider);

    void FragmentRemoval(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector);

    void SimpleTrackEndFilter(HitChargeVector &hitChargeVector);

    void TrackEndFilter(HitChargeVector &hitChargeVector);

    void AttemptFragmentRemoval(const HitChargeVector &hitChargeVector, std::vector<JumpObject> &jumpsVector, HitChargeVector &filteredHitChargeVector);

    void FindLargestJumps(const HitChargeVector &hitChargeVector, std::vector<JumpObject> &leftJumps);

    void FindPeakJumps(const HitChargeVector &hitChargeVector, std::vector<JumpObject> &peakJumps);

    void FindTrackEndJumps(const HitChargeVector &hitChargeVector, std::vector<JumpObject> &trackEndJumps);

    void ParticleSplitting(const pandora::Cluster *const pTargetClusterW, HitChargeVector &hitChargeVector, FitResult &fitResult1, FitResult &fitResult2, bool &splitApplied);

    void FindKinkSize(const pandora::Cluster *const pCluster, float &splitPosition, float &kinkSize);

    void CreateCalorimetricSplitHitVector(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions);

    void SplitHitCollectionBySize(const HitChargeVector &hitChargeVector, float &splitPosition, HitChargeVector &smallHitChargeVector, HitChargeVector &largeHitChargeVector);

    void SplitHitCollectionByLeftRight(const HitChargeVector &hitChargeVector, float &splitPosition, HitChargeVector &leftHitChargeVector, HitChargeVector &rightHitChargeVector);

    void GetTrackLength(const HitChargeVector &hitChargeVector, float &trackLength);

    void GetAverageQoverWTrackBody(HitChargeVector &hitChargeVector, float &averageChargeTrackBody);

    void FindKinkSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions);

    void FindPlateauSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions);

    void FindJumpSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions);

    void FindBowlSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions);

    void FitHitChargeVector(const HitChargeVector &hitChargeVector, TrackDirectionTool::FitResult &fitResult, int numberHitsToConsider=1000000);

    void FitHitChargeVector(HitChargeVector &hitChargeVector1, HitChargeVector &hitChargeVector2, TrackDirectionTool::FitResult &fitResult1, TrackDirectionTool::FitResult &fitResult2, int numberHitsToConsider=1000000);

    void SetGlobalMinuitPreliminaries(const HitChargeVector &hitChargeVector);

    void PerformFits(const HitChargeVector &hitChargeVector, HitChargeVector &forwardsFitPoints, HitChargeVector &backwardsFitPoints, int numberHitsToConsider, float &forwardsChiSquared, float &backwardsChiSquared, int &fitStatus1, int &fitStatus2);

    void GetCalorimetricDirection(const pandora::Cluster* pTargetClusterW, FitResult &finalFitResult);

    void AddToSlidingFitCache(const pandora::Cluster *const pCluster);

    const TwoDSlidingFitResult &GetCachedSlidingFitResult(const pandora::Cluster *const pCluster) const;

    void TidyUp();

    void ClearGlobalVariables();

    static bool SortHitChargeVectorByRL(HitCharge &hitCharge1, HitCharge &hitCharge2);
    static bool SortHitChargeVectorByQoverX(HitCharge &hitCharge1, HitCharge &hitCharge2);
    static bool SortByDistanceToNN(HitCharge &hitCharge1, HitCharge &hitCharge2);
    static bool SortJumpVector(JumpObject &jumpObject1, JumpObject &jumpObject2);

    private:

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    //-----------------------------------------------------------------------------------------------

    unsigned int            m_slidingFitWindow;                 ///< The layer window for the sliding linear fits
    TwoDSlidingFitResultMap m_slidingFitResultMap;              ///< The sliding fit result map

    unsigned int            m_minClusterCaloHits;               ///< The min number of hits in base cluster selection method
    float                   m_minClusterLengthSquared;          ///< The min length (squared) in base cluster selection method

    int                     m_targetParticlePDG;
    int                     m_numberTrackEndHits;
    bool                    m_enableFragmentRemoval;
    bool                    m_enableSplitting;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *TrackDirectionTool::Factory::CreateAlgorithmTool() const
{
    return new TrackDirectionTool();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitCharge::HitCharge(const pandora::CaloHit* caloHit, float &longitudinalPosition, float &hitWidth, float &hitCharge, float &uncertainty) :
    m_intails((hitCharge/hitWidth) <= 1.4 ? true : false),
    m_calohit(caloHit),
    m_longitudinalposition(longitudinalPosition),
    m_hitwidth(hitWidth),
    m_charge(hitCharge),
    m_qoverx(hitCharge/hitWidth),
    m_uncertainty(uncertainty)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitCharge::HitCharge() :
    m_intails(false),
    m_calohit(NULL),
    m_longitudinalposition(0.f),
    m_hitwidth(0.f),
    m_charge(0.f),
    m_qoverx(0.f),
    m_uncertainty(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHit* TrackDirectionTool::HitCharge::GetCaloHit() const
{
    return m_calohit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetLongitudinalPosition() const
{
    return m_longitudinalposition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetHitWidth() const
{
    return m_hitwidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetCharge() const
{
    return m_charge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetQoverX() const
{
    return m_qoverx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetUncertainty() const
{
    return m_uncertainty;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetDistanceToNN(float &distance)
{
    m_distancetonearestneighbour = distance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetDistanceToNN() const
{
    return m_distancetonearestneighbour;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::KinkObject::KinkObject(float &xPosition, float &kinkAngle) :
    m_xposition(xPosition),
    m_kinkangle(kinkAngle),
    m_wrl(0.f),
    m_hit(HitCharge())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::KinkObject::GetXPosition() const
{
    return m_xposition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::KinkObject::GetKinkAngle() const
{
    return m_kinkangle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::KinkObject::SetWLongitudinalPosition(float &rL)
{
    m_wrl = rL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::KinkObject::GetWLongitudinalPosition() const
{
    return m_wrl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::KinkObject::SetNearestHit(HitCharge &hitCharge)
{
    m_hit = hitCharge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitCharge TrackDirectionTool::KinkObject::GetNearestHit() const
{
    return m_hit;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::LookupTable::LookupTable()
{
    std::map<int, double> emptyMap;
    std::map<double, int> emptyReverseMap;

    m_map = emptyMap;
    m_reversemap = emptyReverseMap;
    m_binwidth = 0.f;
    m_initialenergy = 0.f,
    m_maxrange = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::LookupTable::LookupTable(double &initialEnergy, double &binWidth)
{
    std::map<int, double> emptyMap;
    std::map<double, int> emptyReverseMap;

    m_map = emptyMap;
    m_reversemap = emptyReverseMap;
    m_binwidth = binWidth;
    m_initialenergy = initialEnergy,
    m_maxrange = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::map<int, double> TrackDirectionTool::LookupTable::GetMap()
{
    return m_map;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetMap(std::map<int, double> &map)
{
    m_map = map;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::map<double, int> TrackDirectionTool::LookupTable::GetReverseMap()
{
    return m_reversemap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetReverseMap(std::map<double, int> &reverseMap)
{
    m_reversemap = reverseMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TrackDirectionTool::LookupTable::GetInitialEnergy()
{
    return m_initialenergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetInitialEnergy(double &initialEnergy)
{
    m_initialenergy = initialEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TrackDirectionTool::LookupTable::GetBinWidth()
{
    return m_binwidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetBinWidth(double &binWidth)
{
    m_binwidth = binWidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetMaxRange(double &maxRange)
{
    m_maxrange = maxRange;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TrackDirectionTool::LookupTable::GetMaxRange()
{
    return m_maxrange;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::FitResult::FitResult()
{
    HitChargeVector emptyVector;
    m_hitchargevector = (emptyVector);
    m_forwardsrecohits = (emptyVector);
    m_backwardsrecohits = (emptyVector);
    m_nhits = 0;
    m_meanqoverx = 0.f;
    m_forwardschisquared = 0.f;
    m_backwardschisquared = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::FitResult::FitResult(HitChargeVector &hitChargeVector, int &numberHits, float &meanQoverX, float &forwardsChiSquared, float &backwardsChiSquared)
{
    HitChargeVector emptyVector;
    m_hitchargevector = hitChargeVector;
    m_forwardsrecohits = emptyVector;
    m_backwardsrecohits = emptyVector;
    m_nhits = numberHits;
    m_meanqoverx = meanQoverX;
    m_forwardschisquared = forwardsChiSquared;
    m_backwardschisquared = backwardsChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::FitResult::FitResult(HitChargeVector &hitChargeVector, HitChargeVector &forwardsRecoHits, HitChargeVector &backwardsRecoHits, int &numberHits, float &meanQoverX, float &forwardsChiSquared, float &backwardsChiSquared) :
    m_hitchargevector(hitChargeVector),
    m_forwardsrecohits(forwardsRecoHits),
    m_backwardsrecohits(backwardsRecoHits),
    m_nhits(numberHits),
    m_meanqoverx(meanQoverX),
    m_forwardschisquared(forwardsChiSquared),
    m_backwardschisquared(backwardsChiSquared)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitChargeVector TrackDirectionTool::FitResult::GetHitChargeVector()
{
    return m_hitchargevector;
}

//------------------------------------------------------------------------------------------------------------------------------------------
inline TrackDirectionTool::HitChargeVector TrackDirectionTool::FitResult::GetForwardsRecoHits()
{
    return m_forwardsrecohits;
}

//------------------------------------------------------------------------------------------------------------------------------------------
inline TrackDirectionTool::HitChargeVector TrackDirectionTool::FitResult::GetBackwardsRecoHits()
{
    return m_backwardsrecohits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::FitResult::GetForwardsChiSquared()
{
    return m_forwardschisquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::FitResult::GetBackwardsChiSquared()
{
    return m_backwardschisquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TrackDirectionTool::FitResult::GetNHits()
{
    return m_nhits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TrackDirectionTool::FitResult::GetDirectionEstimate()
{
    return (m_forwardschisquared <= m_backwardschisquared ? 1 : 0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::FitResult::GetMinChiSquared()
{
    return std::min(m_forwardschisquared, m_backwardschisquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::FitResult::GetMinChiSquaredPerHit()
{
    return (m_nhits != 0 ? std::min(m_forwardschisquared, m_backwardschisquared)/m_nhits : std::numeric_limits<float>::max());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::FitResult::GetMeanQoverX()
{
    return m_meanqoverx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::FitResult::SetCoordinates(TwoDSlidingFitResult &slidingFitResult)
{
    const pandora::CartesianVector minLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
    const pandora::CartesianVector maxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());

    m_begin_x = minLayerPosition.GetX(); 
    m_begin_y = minLayerPosition.GetY(); 
    m_begin_z = minLayerPosition.GetZ(); 

    m_end_x = maxLayerPosition.GetX(); 
    m_end_y = maxLayerPosition.GetY(); 
    m_end_z = maxLayerPosition.GetZ(); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector TrackDirectionTool::FitResult::GetBeginPoint()
{
    const pandora::CartesianVector beginPoint(m_begin_x, m_begin_y, m_begin_z);
    return beginPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector TrackDirectionTool::FitResult::GetEndPoint()
{
    const pandora::CartesianVector endPoint(m_end_x, m_end_y, m_end_z);
    return endPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::FitResult::ComputeProbability()
{
    TFile *f = new TFile("probability.root"); 

    TH1F* forwardsDeltaChiSquared = (TH1F*)f->Get("forwardsDeltaChiSquared"); 
    TH1F* backwardsDeltaChiSquared = (TH1F*)f->Get("backwardsDeltaChiSquared"); 

    float forwardsBinEntry = forwardsDeltaChiSquared->GetBinContent(forwardsDeltaChiSquared->GetBin(m_forwardschisquared/m_nhits));
    float backwardsBinEntry = backwardsDeltaChiSquared->GetBinContent(backwardsDeltaChiSquared->GetBin(m_backwardschisquared/m_nhits));
    float probability(forwardsBinEntry/(forwardsBinEntry + backwardsBinEntry));

    f->Close();

    m_probability = probability;
    return probability; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::JumpObject::JumpObject(float &longitudinalPosition, float&jumpValue) :
    m_longitudinalposition(longitudinalPosition),
    m_jumpvalue(jumpValue)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::JumpObject::GetLongitudinalPosition()
{
    return m_longitudinalposition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::JumpObject::GetJumpValue()
{
    return m_jumpvalue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_TRACK_DIRECTION_TOOL_H
