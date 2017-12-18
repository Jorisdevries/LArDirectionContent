/**
 *  @file   larpandoracontent/LArVertex/TrackDirectionTool.cc
 *
 *  @brief  Implementation of the candidate vertex creation Tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "LArDirection/TrackDirectionTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include <ctime>

using namespace pandora;

//----------------------------------------------------------------------------------------------------------------------------------

//nasty global parameters necessary for TMinuit
lar_content::TrackDirectionTool::HitChargeVector* pMinuitVector_tool = new lar_content::TrackDirectionTool::HitChargeVector;

float globalTotalCharge_tool(0.f);
float globalTrackLength_tool(0.f);
float globalTotalHitWidth_tool(0.f);

double globalBinWidth_tool(0.1);
double globalInitialEnergy_tool(2000.0);

lar_content::TrackDirectionTool::LookupTable globalMuonLookupTable_tool(globalInitialEnergy_tool, globalBinWidth_tool);
std::string treename_tool("lookuptable"), filename_tool("lookuptable.root");

//----------------------------------------------------------------------------------------------------------------------------------

#include "LArDirection/ToolMinuitFunctions.h"

namespace lar_content
{

TrackDirectionTool::TrackDirectionTool() :
    m_slidingFitWindow(20),
    m_minClusterCaloHits(50),
    m_minClusterLengthSquared(30.f * 30.f),
    m_numberTrackEndHits(100000),
    m_enableFragmentRemoval(false),
    m_enableSplitting(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackDirectionTool::~TrackDirectionTool()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), treename_tool.c_str(), filename_tool.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackDirectionTool::FitResult TrackDirectionTool::Run(pandora::Algorithm *const pAlgorithm, const Cluster *const pTargetClusterW)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

       std::cout << "RUNNING TOOL" << std::endl;

    //FillLookupTable_Tool(globalMuonLookupTable_tool, 105.7);
    this->WriteLookupTableToTree(globalMuonLookupTable_tool);
    this->ReadLookupTableFromTree(globalMuonLookupTable_tool);

    try
    {
        FitResult finalFitResult;

        this->AddToSlidingFitCache(pTargetClusterW);
        this->GetCalorimetricDirection(pTargetClusterW, finalFitResult);

        return finalFitResult;
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }

    this->TidyUp();

    //return finalFitResult;
    //return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackDirectionTool::FitResult TrackDirectionTool::Run(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject* PFO)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    //FillLookupTable_Tool(globalMuonLookupTable_tool, 105.7);
    this->WriteLookupTableToTree(globalMuonLookupTable_tool);
    this->ReadLookupTableFromTree(globalMuonLookupTable_tool);

    const Cluster *const pTargetClusterW = this->GetTargetClusterFromPFO(PFO);

    try
    {
        FitResult finalFitResult;

        this->AddToSlidingFitCache(pTargetClusterW);
        this->GetCalorimetricDirection(pTargetClusterW, finalFitResult);

        return finalFitResult;
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }

    this->TidyUp();

    //return finalFitResult;
    //return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster* TrackDirectionTool::GetTargetClusterFromPFO(const ParticleFlowObject* PFO)
{
    ClusterList clusterList(PFO->GetClusterList());
    
    for (const Cluster *const pCluster : clusterList)
    {
        if (LArClusterHelper::GetClusterHitType(pCluster) == TPC_VIEW_W)
            return pCluster;
    }
    
    std::cout << "ERROR: could not find W cluster!" << std::endl;
    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::WriteLookupTableToTree(LookupTable &lookupTable)
{
    std::vector<int> mapVector1, reverseMapVector2;
    std::vector<double> mapVector2, reverseMapVector1;

    for (auto &element : lookupTable.GetMap())
    {
        mapVector1.push_back(element.first);
        mapVector2.push_back(element.second);
    }

    for (auto &element : lookupTable.GetReverseMap())
    {
        reverseMapVector1.push_back(element.first);
        reverseMapVector2.push_back(element.second);
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treename_tool.c_str(), "mapVector1", &mapVector1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treename_tool.c_str(), "mapVector2", &mapVector2));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treename_tool.c_str(), "reverseMapVector1", &reverseMapVector1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treename_tool.c_str(), "reverseMapVector2", &reverseMapVector2));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treename_tool.c_str(), "binWidth", lookupTable.GetBinWidth()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treename_tool.c_str(), "initialEnergy", lookupTable.GetInitialEnergy()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treename_tool.c_str(), "maxRange", lookupTable.GetMaxRange()));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), treename_tool.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::ReadLookupTableFromTree(LookupTable &lookupTable)
{
    TFile *f = TFile::Open(filename_tool.c_str());
    TTree *t1 = (TTree*)f->Get(treename_tool.c_str());

    std::vector<int> mapVector1, reverseMapVector2;
    std::vector<double> mapVector2, reverseMapVector1;

    std::vector<int> *pMapVector1 = 0;
    std::vector<int> *pReverseMapVector2 = 0;
    std::vector<double> *pMapVector2 = 0;
    std::vector<double> *pReverseMapVector1 = 0;

    TBranch *pBranchMapVector1 = 0;
    TBranch *pBranchMapVector2 = 0;
    TBranch *pBranchReverseMapVector1 = 0;
    TBranch *pBranchReverseMapVector2 = 0;

    t1->SetBranchAddress("mapVector1", &pMapVector1, &pBranchMapVector1);
    t1->SetBranchAddress("mapVector2", &pMapVector2, &pBranchMapVector2);
    t1->SetBranchAddress("reverseMapVector1", &pReverseMapVector1, &pBranchReverseMapVector1);
    t1->SetBranchAddress("reverseMapVector2", &pReverseMapVector2, &pBranchReverseMapVector2);

    const auto tEntry = t1->LoadTree(0);
    pBranchMapVector1->GetEntry(tEntry);
    pBranchMapVector2->GetEntry(tEntry);
    pBranchReverseMapVector1->GetEntry(tEntry);
    pBranchReverseMapVector2->GetEntry(tEntry);

    for (int j = 0; j < (int)pMapVector1->size(); ++j)
    {
        mapVector1.push_back(pMapVector1->at(j));
        mapVector2.push_back(pMapVector2->at(j));
        reverseMapVector1.push_back(pReverseMapVector1->at(j));
        reverseMapVector2.push_back(pReverseMapVector2->at(j));
    }

    std::map<int, double> map;
    std::map<double, int> reverseMap;
    double binWidth;
    double initialEnergy;
    double maxRange;

    t1->SetBranchAddress("binWidth", &binWidth);
    t1->SetBranchAddress("initialEnergy", &initialEnergy);
    t1->SetBranchAddress("maxRange", &maxRange);
    t1->GetEntry(0);

    for (int i = 0; i < mapVector1.size(); i++)
        map.insert(std::pair<int,double>(mapVector1.at(i), mapVector2.at(i)));

    for (int i = 0; i < mapVector1.size(); i++)
        reverseMap.insert(std::pair<double, int>(reverseMapVector1.at(i), reverseMapVector2.at(i)));


    lookupTable.SetMap(map);
    lookupTable.SetReverseMap(reverseMap);
    lookupTable.SetBinWidth(binWidth);
    lookupTable.SetInitialEnergy(initialEnergy);
    lookupTable.SetMaxRange(maxRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FillHitChargeVector(const Cluster *const pCluster, HitChargeVector &hitChargeVector)
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    const TwoDSlidingFitResult &slidingFitResult(this->GetCachedSlidingFitResult(pCluster));

    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
    {
        const CaloHit *const pCaloHit(*hitIter);
        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        float hitWidth(pCaloHit->GetCellSize1());

        float caloHitEnergy(pCaloHit->GetInputEnergy());
        caloHitEnergy *= 273.5; //ADC to electron
        caloHitEnergy *= 23.6/1000000; //ionisation energy per electron in MeV
        caloHitEnergy /= 0.62;

        float rL(0.f), rT(0.f);
        slidingFitResult.GetLocalPosition(caloHitPosition, rL, rT);
        if (rL == 0.)
            continue;

        float calibratedUncertainty(std::sqrt((0.00419133 * (caloHitEnergy/hitWidth) * (caloHitEnergy/hitWidth)) + (0.00967141 * (caloHitEnergy/hitWidth)))); //70%
        HitCharge hitCharge(pCaloHit, rL, hitWidth, caloHitEnergy, calibratedUncertainty);
        hitChargeVector.push_back(hitCharge);
    }

    std::sort(hitChargeVector.begin(), hitChargeVector.end(), SortHitChargeVectorByRL);

}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TrackInnerFilter(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector)
{
    //Fill endpoint protected area into filtered vector and put all other hits in a separate vector
    float endpointProtectionRange(0.05);
    filteredHitChargeVector.insert(filteredHitChargeVector.begin(), hitChargeVector.begin(),  hitChargeVector.begin() + endpointProtectionRange * hitChargeVector.size());
    filteredHitChargeVector.insert(filteredHitChargeVector.begin(), hitChargeVector.begin() + (1.0 - endpointProtectionRange) * hitChargeVector.size(), hitChargeVector.end());

    HitChargeVector innerHitChargeVector(hitChargeVector.begin() + endpointProtectionRange * hitChargeVector.size(), hitChargeVector.begin() + (1.0 - endpointProtectionRange) * hitChargeVector.size());

    int nNeighboursToConsider(5);
    this->SetNearestNeighbourValues(innerHitChargeVector, nNeighboursToConsider);

     std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortByDistanceToNN);
     filteredHitChargeVector.insert(filteredHitChargeVector.begin(), innerHitChargeVector.begin(), innerHitChargeVector.begin() + 0.72 * innerHitChargeVector.size()); //lots of testing has been done to optimise percentage
     std::sort(filteredHitChargeVector.begin(), filteredHitChargeVector.end(), SortHitChargeVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetNearestNeighbourValues(HitChargeVector &innerHitChargeVector, int &nNeighboursToConsider)
{
    float trackLength(0.f);
    this->GetTrackLength(innerHitChargeVector, trackLength);

    std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortHitChargeVectorByQoverX);
    float QoverXRange((*innerHitChargeVector.end()).GetQoverX() - (*innerHitChargeVector.begin()).GetQoverX());

    for (HitCharge &hitCharge1 : innerHitChargeVector)
    {
        std::vector<float> distancesToNN;

        for (HitCharge &hitCharge2 : innerHitChargeVector)
        {
            if (&hitCharge1 == &hitCharge2)
                continue;

            float QoverXDistance((trackLength/QoverXRange) * (std::abs(hitCharge1.GetQoverX() - hitCharge2.GetQoverX())));
            float Ldistance(std::abs(hitCharge1.GetLongitudinalPosition() - hitCharge2.GetLongitudinalPosition()));
            float distanceToNN(std::sqrt(QoverXDistance*QoverXDistance + Ldistance*Ldistance));

            distancesToNN.push_back(distanceToNN);
        }

        std::sort(distancesToNN.begin(), distancesToNN.end());
        float nearestNeighboursDistanceSum(std::accumulate(distancesToNN.begin(), distancesToNN.begin() + nNeighboursToConsider, 0.f));
        hitCharge1.SetDistanceToNN(nearestNeighboursDistanceSum);
    }

    std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortHitChargeVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FragmentRemoval(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector)
{
    float trackLength(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    std::vector<JumpObject> jumpsVector;
    this->FindLargestJumps(hitChargeVector, jumpsVector);

    std::vector<JumpObject> peakJumps;
    this->FindPeakJumps(hitChargeVector, jumpsVector);

    this->AttemptFragmentRemoval(hitChargeVector, jumpsVector, filteredHitChargeVector);
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SimpleTrackEndFilter(HitChargeVector &hitChargeVector)
{
    float lowerBound(0.9), upperBound(2.2);

    while ((*(hitChargeVector.begin())).GetQoverX()/(*(std::next(hitChargeVector.begin(), 1))).GetQoverX() <= lowerBound || (*(hitChargeVector.begin())).GetQoverX()/(*(std::next(hitChargeVector.begin(), 1))).GetQoverX() >= upperBound)
        hitChargeVector.erase(hitChargeVector.begin());

    while ((*(std::prev(hitChargeVector.end(), 1))).GetQoverX()/(*(std::prev(hitChargeVector.end(), 2))).GetQoverX() <= lowerBound || (*(std::prev(hitChargeVector.end(), 1))).GetQoverX()/(*(std::prev(hitChargeVector.end(), 2))).GetQoverX() >= upperBound)
        hitChargeVector.pop_back();

    //This piece of logic removes hits that have uncharacteristically high or low Q/w values (in tails of Q/w distribution)
    hitChargeVector.erase(
    std::remove_if(hitChargeVector.begin(), hitChargeVector.end(),
        [](const HitCharge & hitCharge) { return hitCharge.m_intails; }),
    hitChargeVector.end());
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TrackEndFilter(HitChargeVector &hitChargeVector)
{
    float trackLength(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    FitResult beforeFitResult;
    this->FitHitChargeVector(hitChargeVector, beforeFitResult);

    int beforeNumberHits(hitChargeVector.size());
    float bodyQoverW(0.f);
    this->GetAverageQoverWTrackBody(hitChargeVector, bodyQoverW);

    int nHitsToSkip(3);
    float trackEndRange(0.025);

    HitChargeVector filteredHitChargeVector(hitChargeVector);

    for (HitChargeVector::const_iterator iter = filteredHitChargeVector.begin(); iter != filteredHitChargeVector.end(); )
    {
        HitCharge hitCharge(*iter), nextHitCharge(*std::next(iter, 1)), plusNHitCharge(*std::next(iter, nHitsToSkip)), previousHitCharge(*std::prev(iter, 1)), minusNHitCharge(*std::prev(iter, nHitsToSkip));

        if (hitCharge.GetLongitudinalPosition()/trackLength <= trackEndRange || hitCharge.GetLongitudinalPosition()/trackLength >= (1.0 - trackEndRange))
        {
            float nearestRatio(std::max((hitCharge.GetQoverX()/previousHitCharge.GetQoverX()), (hitCharge.GetQoverX()/nextHitCharge.GetQoverX())));
            float plusMinusNRatio(std::max((hitCharge.GetQoverX()/minusNHitCharge.GetQoverX()), (hitCharge.GetQoverX()/plusNHitCharge.GetQoverX())));
            float distanceFromBodyQoverW(std::abs(hitCharge.GetQoverX() - bodyQoverW));

            if (previousHitCharge.GetQoverX() < 0.01)
                nearestRatio = hitCharge.GetQoverX()/nextHitCharge.GetQoverX();
            if (nextHitCharge.GetQoverX() < 0.01)
                nearestRatio = hitCharge.GetQoverX()/previousHitCharge.GetQoverX();
            if (minusNHitCharge.GetQoverX() < 0.01)
                plusMinusNRatio = hitCharge.GetQoverX()/plusNHitCharge.GetQoverX();
            if (plusNHitCharge.GetQoverX() < 0.01)
                plusMinusNRatio = hitCharge.GetQoverX()/minusNHitCharge.GetQoverX();

            if (distanceFromBodyQoverW >= 4.0 || std::abs(1.0 - nearestRatio) >= 0.5 || std::abs(1.0 - plusMinusNRatio) >= 0.5)
                iter = filteredHitChargeVector.erase(iter);
            else
                ++iter;
        }
        else
        {
            ++iter;
        }
    }

    std::sort(filteredHitChargeVector.begin(), filteredHitChargeVector.end(), SortHitChargeVectorByRL);

    FitResult afterFitResult;
    this->FitHitChargeVector(filteredHitChargeVector, afterFitResult);

    float ChiSquaredPerHitChange(beforeFitResult.GetMinChiSquaredPerHit() - afterFitResult.GetMinChiSquaredPerHit());
    float N(beforeNumberHits);
    bool shouldApply(true);

    if (beforeNumberHits < 400 && ChiSquaredPerHitChange < (8.0 - ((N/400) * 7.0)))
        shouldApply = false;
    if (beforeNumberHits >= 400 && ChiSquaredPerHitChange < 1.0)
        shouldApply = false;

    if (shouldApply)
        hitChargeVector = filteredHitChargeVector;
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::AttemptFragmentRemoval(const HitChargeVector &hitChargeVector, std::vector<JumpObject> &jumpsVector, HitChargeVector &filteredHitChargeVector)
{
    FitResult beforeFitResult;
    this->FitHitChargeVector(hitChargeVector, beforeFitResult);

    float beforeRemovalMinchiSquared(beforeFitResult.GetMinChiSquaredPerHit());

    HitChargeVector bestHitChargeVector;
    FitResult bestFitResult(beforeFitResult);

    for (JumpObject &jumpObject : jumpsVector)
    {
        float splitPosition(jumpObject.GetLongitudinalPosition());

        HitChargeVector smallHitCollection, largeHitCollection;
        this->SplitHitCollectionBySize(hitChargeVector, splitPosition, smallHitCollection, largeHitCollection);

        FitResult afterFitResult;
        this->FitHitChargeVector(largeHitCollection, afterFitResult);

        if (afterFitResult.GetMinChiSquaredPerHit() < bestFitResult.GetMinChiSquaredPerHit())
        {
            bestHitChargeVector = largeHitCollection;
            bestFitResult = afterFitResult;
        }
    }

    if (beforeRemovalMinchiSquared - bestFitResult.GetMinChiSquaredPerHit() > 2.0)
        filteredHitChargeVector = bestHitChargeVector;
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindLargestJumps(const HitChargeVector &hitChargeVector, std::vector<JumpObject> &normalJumps)
{
    //HitChargeVector binnedHitChargeVector;
    //BinHitChargeVector_Tool(hitChargeVector, binnedHitChargeVector, 0.5);

    std::vector<HitChargeVector> bothVectors;
    bothVectors.push_back(hitChargeVector);
    //bothVectors.push_back(binnedHitChargeVector);

    for (HitChargeVector &vector : bothVectors)
    {
        int searchRange(0.05 * vector.size());

        for (int jumpRange = 1; jumpRange <= 5; jumpRange++)
        {
            for (int i = 0; i < searchRange; i++)
            {
                float binJump = (std::abs(vector.at(i).GetQoverX() - vector.at(i + jumpRange).GetQoverX()));
                float jumpPosition(vector.at(i + jumpRange).GetLongitudinalPosition());
                JumpObject jumpObject(jumpPosition, binJump);
                normalJumps.push_back(jumpObject);
            }

            for (int j = vector.size() - searchRange; j < vector.size() - jumpRange; j++)
            {
                float binJump = (std::abs(vector.at(j).GetQoverX() - vector.at(j + jumpRange).GetQoverX()));
                float jumpPosition(vector.at(j).GetLongitudinalPosition());
                JumpObject jumpObject(jumpPosition, binJump);

                normalJumps.push_back(jumpObject);
            }
        }
    }

    std::sort(normalJumps.begin(), normalJumps.end(), SortJumpVector);
    if (normalJumps.size() > 3)
        normalJumps.erase(normalJumps.begin() + 3, normalJumps.end());
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindPeakJumps(const HitChargeVector &hitChargeVector, std::vector<JumpObject> &peakJumps)
{
    float jumpPosition(0.f), jumpValue(0.f), currentLargestQoverW(0.f);

    for (const HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetQoverX() > currentLargestQoverW)
        {
            currentLargestQoverW = hitCharge.GetQoverX();
            jumpPosition = hitCharge.GetLongitudinalPosition();
        }
    }

    float jumpPosition1(jumpPosition - 0.5), jumpPosition2(jumpPosition + 0.5);
    JumpObject jumpObject(jumpPosition, jumpValue);
    JumpObject jumpObject1(jumpPosition1, jumpValue);
    JumpObject jumpObject2(jumpPosition2, jumpValue);
    peakJumps.push_back(jumpObject);
    peakJumps.push_back(jumpObject1);
    peakJumps.push_back(jumpObject2);
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindTrackEndJumps(const HitChargeVector &hitChargeVector, std::vector<JumpObject> &trackEndJumps)
{
    float trackLength(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    for (float edge = 0.01; edge <= 0.25; edge += 0.01)
    {
        float jumpPosition1(edge * trackLength), jumpPosition2((1.0 - edge) * trackLength), jumpValue(0.f);
        JumpObject jumpObject1(jumpPosition1, jumpValue);
        JumpObject jumpObject2(jumpPosition2, jumpValue);

        trackEndJumps.push_back(jumpObject1);
        trackEndJumps.push_back(jumpObject2);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::ParticleSplitting(const Cluster* pTargetClusterW, HitChargeVector &hitChargeVector, FitResult &backwardsFitResult, FitResult &forwardsFitResult, bool &splitApplied)
{
    FitResult beforeFitResult;
    this->FitHitChargeVector(hitChargeVector, beforeFitResult);
    FitResult outputBackwardsFitResult(beforeFitResult), outputForwardsFitResult(beforeFitResult);
    backwardsFitResult = beforeFitResult;
    forwardsFitResult = beforeFitResult;

    float afterSplitChiSquared(beforeFitResult.GetMinChiSquaredPerHit());

    std::vector<float> calorimetricSplitPositions;
    this->CreateCalorimetricSplitHitVector(hitChargeVector, calorimetricSplitPositions);

    for (float &splitPosition : calorimetricSplitPositions)
    {
        HitChargeVector backwardsTestHitCollection, forwardsTestHitCollection;
        this->SplitHitCollectionByLeftRight(hitChargeVector, splitPosition, backwardsTestHitCollection, forwardsTestHitCollection);

        FitResult backwardsTestFitResult, forwardsTestFitResult;
        this->FitHitChargeVector(backwardsTestHitCollection, forwardsTestHitCollection, backwardsTestFitResult, forwardsTestFitResult);

        float splitMinChiSquared((backwardsTestFitResult.GetNHits() > 0 ? backwardsTestFitResult.GetBackwardsChiSquared()/backwardsTestFitResult.GetNHits() : 0.f) + (forwardsTestFitResult.GetNHits() > 0 ? forwardsTestFitResult.GetForwardsChiSquared()/forwardsTestFitResult.GetNHits() : 0.f));

        float kinkSize(0.f);
        this->FindKinkSize(pTargetClusterW, splitPosition, kinkSize);

        if (splitMinChiSquared < afterSplitChiSquared)
        {
            afterSplitChiSquared = splitMinChiSquared;
            outputBackwardsFitResult = backwardsTestFitResult;
            outputForwardsFitResult = forwardsTestFitResult;
        }
    }

    if (outputBackwardsFitResult.GetNHits() == 0 || outputForwardsFitResult.GetNHits() == 0)
        return;

    float ChiSquaredPerHitChange(beforeFitResult.GetMinChiSquaredPerHit() - (outputBackwardsFitResult.GetBackwardsChiSquared()/outputBackwardsFitResult.GetNHits() + outputForwardsFitResult.GetForwardsChiSquared()/outputForwardsFitResult.GetNHits()));
    int beforeNumberHits((int)beforeFitResult.GetHitChargeVector().size());
    float N(beforeNumberHits);
    bool shouldApply(true);

    if (beforeNumberHits < 400 && ChiSquaredPerHitChange < (5.0 - ((N/400) * 4.0)))
        shouldApply = false;
    if (beforeNumberHits >= 400 && ChiSquaredPerHitChange < 1.0)
        shouldApply = false;

    if (shouldApply)
    {
        splitApplied = true;
        backwardsFitResult = outputBackwardsFitResult;
        forwardsFitResult = outputForwardsFitResult;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindKinkSize(const Cluster* pCluster, float &splitPosition, float &kinkSize)
{
    try
    {
        const TwoDSlidingFitResult &slidingFitResult(this->GetCachedSlidingFitResult(pCluster));
        const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());
        const int minLayer(layerFitResultMap.begin()->first), maxLayer(layerFitResultMap.rbegin()->first);

        const int nLayersHalfWindow(slidingFitResult.GetLayerFitHalfWindow());
        const int nLayersSpanned(1 + maxLayer - minLayer);

        if (nLayersSpanned <= 2 * nLayersHalfWindow)
            return;

        for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
        {
            const int iLayer(iter->first);

            const float rL(slidingFitResult.GetL(iLayer));
            const float rL1(slidingFitResult.GetL(iLayer - nLayersHalfWindow));
            const float rL2(slidingFitResult.GetL(iLayer + nLayersHalfWindow));

            CartesianVector centralPosition(0.f,0.f,0.f), firstDirection(0.f,0.f,0.f), secondDirection(0.f,0.f,0.f);

            if ((STATUS_CODE_SUCCESS != slidingFitResult.GetGlobalFitPosition(rL, centralPosition)) ||
                (STATUS_CODE_SUCCESS != slidingFitResult.GetGlobalFitDirection(rL1, firstDirection)) ||
                (STATUS_CODE_SUCCESS != slidingFitResult.GetGlobalFitDirection(rL2, secondDirection)))
            {
                continue;
            }

            const float cosTheta(firstDirection.GetDotProduct(secondDirection));
            if (std::abs(splitPosition - rL) <= 3.0 && cosTheta > kinkSize && cosTheta < 1.0)
                kinkSize = (180.0 / 3.1415926535) * std::acos(cosTheta);
        }
    }
    catch (...)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::CreateCalorimetricSplitHitVector(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions)
{
    this->FindKinkSplit(hitChargeVector, splitPositions);
    this->FindJumpSplit(hitChargeVector, splitPositions);
    this->FindPlateauSplit(hitChargeVector, splitPositions);
    this->FindBowlSplit(hitChargeVector, splitPositions);

    sort( splitPositions.begin(), splitPositions.end() );
splitPositions.erase( unique( splitPositions.begin(), splitPositions.end() ), splitPositions.end() );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SplitHitCollectionBySize(const HitChargeVector &hitChargeVector, float &splitPosition, HitChargeVector &smallHitChargeVector, HitChargeVector &largeHitChargeVector)
{
    HitChargeVector leftHitCollection, rightHitCollection;

    for (const  HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() <= splitPosition)
            leftHitCollection.push_back(hitCharge);
        else
            rightHitCollection.push_back(hitCharge);
    }

    if (leftHitCollection.size() <= rightHitCollection.size())
    {
        smallHitChargeVector = leftHitCollection;
        largeHitChargeVector = rightHitCollection;
    }
    else
    {
        smallHitChargeVector = rightHitCollection;
        largeHitChargeVector = leftHitCollection;
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SplitHitCollectionByLeftRight(const HitChargeVector &hitChargeVector, float &splitPosition, HitChargeVector &leftHitChargeVector, HitChargeVector &rightHitChargeVector)
{
    for (const  HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() <= splitPosition)
            leftHitChargeVector.push_back(hitCharge);
        else
            rightHitChargeVector.push_back(hitCharge);
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetTrackLength(const HitChargeVector &hitChargeVector, float &trackLength)
{
    trackLength = 0.f;

    for (const HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() > trackLength)
            trackLength = hitCharge.GetLongitudinalPosition();
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetAverageQoverWTrackBody(HitChargeVector &hitChargeVector, float &averageChargeTrackBody)
{
    //temp vector because I do not want to mess with the sorting of the original vector
    HitChargeVector tempHitChargeVector;

    for (HitCharge &hitCharge : hitChargeVector)
        tempHitChargeVector.push_back(hitCharge);

    std::sort(tempHitChargeVector.begin(), tempHitChargeVector.end(), SortHitChargeVectorByQoverX);

    int nEntries(0);

    for (int q = 0; q < tempHitChargeVector.size() - 1; q++)
    {
        if (q <= 0.1 * tempHitChargeVector.size() || q >= 0.6 * tempHitChargeVector.size())
            continue;

        averageChargeTrackBody += tempHitChargeVector.at(q).GetQoverX();
        nEntries++;
    }

    averageChargeTrackBody /= nEntries;
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindKinkSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions)
{
    HitChargeVector binnedHitChargeVector = hitChargeVector;
    //BinHitChargeVector_Tool(hitChargeVector, binnedHitChargeVector, 1.0);

    std::vector<JumpObject> kinkObjects;

    float minCharge(10000.f), maxCharge(0.f);

    for (HitCharge &hitCharge : binnedHitChargeVector)
    {
        if (hitCharge.GetCharge() < minCharge)
            minCharge = hitCharge.GetCharge();

        if (hitCharge.GetCharge() > maxCharge)
            maxCharge = hitCharge.GetCharge();
    }

    float fullChargeRange(maxCharge - minCharge);

    float chargeHalfWidth(0.1 * fullChargeRange);

    for (HitCharge &bin1 : binnedHitChargeVector)
    {
        HitChargeVector leftHitCollection, rightHitCollection;

        for (HitCharge &vector : binnedHitChargeVector)
        {
            if (vector.GetLongitudinalPosition() <= bin1.GetLongitudinalPosition())
                leftHitCollection.push_back(vector);
            else
                rightHitCollection.push_back(vector);
        }

        if (leftHitCollection.size() == 0 || rightHitCollection.size() == 0)
            continue;

        float bestLeftScore(0.f);

        for (HitCharge &bin2 : leftHitCollection)
        {
            float chargeDifference(bin2.GetCharge() - bin1.GetCharge());
            float positionDifference(bin1.GetLongitudinalPosition() - bin2.GetLongitudinalPosition());
            float slope(chargeDifference/positionDifference);

            int nLeftHits(0);
            for (HitCharge &bin3 : leftHitCollection)
            {
                float lineValue(bin2.GetCharge() - ((bin3.GetLongitudinalPosition() - bin2.GetLongitudinalPosition()) * slope));
                if (std::abs(bin3.GetCharge() - lineValue) < chargeHalfWidth)
                    nLeftHits++;
            }

            float score((float)nLeftHits/leftHitCollection.size());

            if (score > bestLeftScore)
                bestLeftScore = score;
        }


        float bestRightScore(0.f);

        for (HitCharge &bin2 : rightHitCollection)
        {
            float chargeDifference(bin2.GetCharge() - bin1.GetCharge());
            float positionDifference(bin1.GetLongitudinalPosition() - bin2.GetLongitudinalPosition());
            float slope(chargeDifference/positionDifference);

            int nRightHits(0);
            for (HitCharge &bin3 : rightHitCollection)
            {
                float lineValue(bin2.GetCharge() - ((bin3.GetLongitudinalPosition() - bin2.GetLongitudinalPosition()) * slope));
                if (std::abs(bin3.GetCharge() - lineValue) < chargeHalfWidth)
                    nRightHits++;
            }

            float score((float)nRightHits/rightHitCollection.size());

            if (score > bestRightScore)
                bestRightScore = score;
        }

        float kinkPosition(bin1.GetLongitudinalPosition());
        float totalScore(bestLeftScore + bestRightScore);
        JumpObject kinkObject(kinkPosition, totalScore);
        kinkObjects.push_back(kinkObject);
    }

    std::sort(kinkObjects.begin(), kinkObjects.end(), SortJumpVector);

    int cutOff(3), nAdded(0);
    float latestJumpPosition(0.f), latestJumpValue(0.f), range(3.f);

    for (int i = 0; i < kinkObjects.size(); i++)
    {
        if (nAdded >= cutOff)
            break;

        JumpObject kinkObject(kinkObjects.at(i));

        if (kinkObject.GetJumpValue() < 1.5)
            continue;

        if (nAdded == 0)
        {
            splitPositions.push_back(kinkObjects.at(i).GetLongitudinalPosition());
            latestJumpPosition = kinkObjects.at(i).GetLongitudinalPosition();
            latestJumpValue = kinkObjects.at(i).GetJumpValue();
            nAdded++;
        }
        else
        {
            if ((kinkObject.GetLongitudinalPosition() - latestJumpPosition) < range && kinkObject.GetJumpValue() > latestJumpValue)
            {
                splitPositions.pop_back();
                splitPositions.push_back(kinkObjects.at(i).GetLongitudinalPosition());
                latestJumpPosition = kinkObjects.at(i).GetLongitudinalPosition();
                latestJumpValue = kinkObjects.at(i).GetJumpValue();
            }
            else if ((kinkObject.GetLongitudinalPosition() - latestJumpPosition) > range)
            {
                splitPositions.push_back(kinkObjects.at(i).GetLongitudinalPosition());
                latestJumpPosition = kinkObjects.at(i).GetLongitudinalPosition();
                latestJumpValue = kinkObjects.at(i).GetJumpValue();
                nAdded++;
            }
        }
    }
}
//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindPlateauSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions)
{
    float trackLength(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    float averageCharge(0.f);
    this->GetAverageQoverWTrackBody(hitChargeVector, averageCharge);

    std::vector<JumpObject> plateauObjects;


    /////////////CHANGE SETTINGS HERE/////////////

    float positionStepSize(0.1);
    float chargeStep(0.10);
    float trackScanRange(0.05 * trackLength);

    //////////////////////////////////////////////

    for (float currentPosition = 0; currentPosition < trackLength; currentPosition += positionStepSize)
    {
        int totalHitsLeft(0), totalHitsRight(0);

        for (const HitCharge &hitCharge : hitChargeVector)
        {
            if (std::abs(currentPosition - hitCharge.GetLongitudinalPosition()) > trackScanRange)
                continue;

            if (hitCharge.GetLongitudinalPosition() <= currentPosition)
                totalHitsLeft++;
            else
                totalHitsRight++;
        }

        float bestHitFractionLeft(0), bestHitFractionRight(0);
        float bestChargeLeft(0), bestChargeRight(0);

        for (float currentCharge = chargeStep; currentCharge < 10.0; currentCharge += chargeStep)
        {
            float hitCountLeft(0.f), hitCountRight(0.f);

            for (const HitCharge &hitCharge : hitChargeVector)
            {
                if (std::abs(currentPosition - hitCharge.GetLongitudinalPosition()) > trackScanRange)
                    continue;

                if (hitCharge.GetLongitudinalPosition() <= currentPosition && hitCharge.GetCharge() > (currentCharge - chargeStep) && hitCharge.GetCharge() < (currentCharge + chargeStep))
                    hitCountLeft += 1.0;

                if (hitCharge.GetLongitudinalPosition() > currentPosition && hitCharge.GetCharge() > (currentCharge - chargeStep) && hitCharge.GetCharge() < (currentCharge + chargeStep))
                    hitCountRight += 1.0;
            }

            float hitFractionLeft(hitCountLeft/totalHitsLeft), hitFractionRight(hitCountRight/totalHitsRight);

            if (hitFractionLeft > bestHitFractionLeft) // && hitCountLeft >= 5
            {
                bestHitFractionLeft = hitFractionLeft;
                bestChargeLeft = currentCharge;
            }

            if (hitFractionRight > bestHitFractionRight) //&& hitCountRight >= 5
            {
                bestHitFractionRight = hitFractionRight;
                bestChargeRight = currentCharge;
            }
        }

        float chargeRange(std::abs(bestChargeLeft - bestChargeRight));

        float currentScore(bestHitFractionLeft + bestHitFractionRight);
        currentScore *= chargeRange/averageCharge;
        JumpObject plateauObject(currentPosition, currentScore);

        plateauObjects.push_back(plateauObject);
    }

    std::sort(plateauObjects.begin(), plateauObjects.end(), SortJumpVector);
    int cutOff(3), nAdded(0);
    float latestJumpPosition(0.f), latestJumpValue(0.f), range(3.f);

    for (int i = 0; i < plateauObjects.size(); i++)
    {
        if (nAdded >= cutOff)
            break;

        JumpObject plateauObject(plateauObjects.at(i));

        if (plateauObject.GetJumpValue() < 0.1)
            continue;

        if (nAdded == 0)
        {
            splitPositions.push_back(plateauObjects.at(i).GetLongitudinalPosition());
            latestJumpPosition = plateauObjects.at(i).GetLongitudinalPosition();
            latestJumpValue = plateauObjects.at(i).GetJumpValue();
            nAdded++;
        }
        else
        {
            if ((plateauObject.GetLongitudinalPosition() - latestJumpPosition) < range && plateauObject.GetJumpValue() > latestJumpValue)
            {
                splitPositions.pop_back();
                splitPositions.push_back(plateauObjects.at(i).GetLongitudinalPosition());
                latestJumpPosition = plateauObjects.at(i).GetLongitudinalPosition();
                latestJumpValue = plateauObjects.at(i).GetJumpValue();
            }
            else if ((plateauObject.GetLongitudinalPosition() - latestJumpPosition) > range)
            {
                splitPositions.push_back(plateauObjects.at(i).GetLongitudinalPosition());
                latestJumpPosition = plateauObjects.at(i).GetLongitudinalPosition();
                latestJumpValue = plateauObjects.at(i).GetJumpValue();
                nAdded++;
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindJumpSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions)
{
    //HitChargeVector binnedHitChargeVector;
    //BinHitChargeVector_Tool(hitChargeVector, binnedHitChargeVector, 1.0);

    std::vector<JumpObject> normalJumps, binnedJumps;

    std::vector<HitChargeVector> bothVectors;
    bothVectors.push_back(hitChargeVector);
    //bothVectors.push_back(binnedHitChargeVector);

    for (HitChargeVector &vector : bothVectors)
    {
        for (int jumpRange = 1; jumpRange <= 5; jumpRange++)
        {
            for (int i = 0; i <= vector.size() - (jumpRange + 1); i++)
            {
                float combinedUncertainty = vector.at(i).GetUncertainty() + vector.at(i + jumpRange).GetUncertainty();
                float binJump = (std::abs(vector.at(i).GetQoverX() - vector.at(i + jumpRange).GetQoverX()));
                binJump /= combinedUncertainty;
                float jumpPosition(vector.at(i + jumpRange).GetLongitudinalPosition());
                JumpObject jumpObject(jumpPosition, binJump);

                if (vector.size() == hitChargeVector.size())
                    normalJumps.push_back(jumpObject);

                //if (vector.size() == binnedHitChargeVector.size())
                //    binnedJumps.push_back(jumpObject);
            }
        }
    }

    std::sort(normalJumps.begin(), normalJumps.end(), SortJumpVector);
    std::sort(binnedJumps.begin(), binnedJumps.end(), SortJumpVector);

    int cutOff(10), nAdded(0);
    float latestJumpPosition(0.f), latestJumpValue(0.f), range(3.f);

    for (int i = 0; i < normalJumps.size(); i++)
    {
        if (nAdded >= cutOff)
            break;

        JumpObject jumpObject(normalJumps.at(i));

        if (jumpObject.GetJumpValue() < 1.0)
            continue;

        if (nAdded == 0)
        {
            splitPositions.push_back(normalJumps.at(i).GetLongitudinalPosition());
            latestJumpPosition = normalJumps.at(i).GetLongitudinalPosition();
            latestJumpValue = normalJumps.at(i).GetJumpValue();
            nAdded++;
        }
        else
        {
            if ((jumpObject.GetLongitudinalPosition() - latestJumpPosition) < range && jumpObject.GetJumpValue() > latestJumpValue)
            {
                splitPositions.pop_back();
                splitPositions.push_back(normalJumps.at(i).GetLongitudinalPosition());
                latestJumpPosition = normalJumps.at(i).GetLongitudinalPosition();
                latestJumpValue = normalJumps.at(i).GetJumpValue();
            }
            else if ((jumpObject.GetLongitudinalPosition() - latestJumpPosition) > range)
            {
                splitPositions.push_back(normalJumps.at(i).GetLongitudinalPosition());
                latestJumpPosition = normalJumps.at(i).GetLongitudinalPosition();
                latestJumpValue = normalJumps.at(i).GetJumpValue();
                nAdded++;
            }
        }
    }

    nAdded = 0;

    for (int i = 0; i < binnedJumps.size(); i++)
    {
        if (nAdded >= cutOff)
            break;

        JumpObject jumpObject(binnedJumps.at(i));

        if (jumpObject.GetJumpValue() < 1.0)
            continue;

        if (nAdded == 0)
        {
            splitPositions.push_back(binnedJumps.at(i).GetLongitudinalPosition());
            latestJumpPosition = binnedJumps.at(i).GetLongitudinalPosition();
            latestJumpValue = binnedJumps.at(i).GetJumpValue();
            nAdded++;
        }
        else
        {
            if ((jumpObject.GetLongitudinalPosition() - latestJumpPosition) < range && jumpObject.GetJumpValue() > latestJumpValue)
            {
                splitPositions.pop_back();
                splitPositions.push_back(binnedJumps.at(i).GetLongitudinalPosition());
                latestJumpPosition = binnedJumps.at(i).GetLongitudinalPosition();
                latestJumpValue = binnedJumps.at(i).GetJumpValue();
            }
            else if ((jumpObject.GetLongitudinalPosition() - latestJumpPosition) > range)
            {
                splitPositions.push_back(binnedJumps.at(i).GetLongitudinalPosition());
                latestJumpPosition = binnedJumps.at(i).GetLongitudinalPosition();
                latestJumpValue = binnedJumps.at(i).GetJumpValue();
                nAdded++;
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FindBowlSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions)
{
    //HitChargeVector binnedHitChargeVector;
    //BinHitChargeVector_Tool(hitChargeVector, binnedHitChargeVector, 1.0);
    splitPositions.push_back(hitChargeVector.at(hitChargeVector.size()/2).GetLongitudinalPosition());
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FitHitChargeVector(const HitChargeVector &hitChargeVector, TrackDirectionTool::FitResult &fitResult, int numberHitsToConsider)
{
    float particleForwardsChiSquared(0.f), particleBackwardsChiSquared(0.f);
    int numberHits(std::min(2 * numberHitsToConsider, (int)hitChargeVector.size())), particleForwardsFitStatus(-1), particleBackwardsFitStatus(-1);
    HitChargeVector forwardsFitPoints, backwardsFitPoints;
    this->PerformFits(hitChargeVector, forwardsFitPoints, backwardsFitPoints, numberHitsToConsider, particleForwardsChiSquared, particleBackwardsChiSquared, particleForwardsFitStatus, particleBackwardsFitStatus);

    float mean_dEdx(0.f);
    HitChargeVector thisHitChargeVector = hitChargeVector;
    for (HitCharge &hitCharge : thisHitChargeVector)
        mean_dEdx += hitCharge.GetQoverX();
    mean_dEdx /= thisHitChargeVector.size();

    std::sort(thisHitChargeVector.begin(), thisHitChargeVector.end(), SortHitChargeVectorByRL);
    std::sort(forwardsFitPoints.begin(), forwardsFitPoints.end(), SortHitChargeVectorByRL);
    std::sort(backwardsFitPoints.begin(), backwardsFitPoints.end(), SortHitChargeVectorByRL);

    FitResult finalFitResult(thisHitChargeVector, forwardsFitPoints, backwardsFitPoints, numberHits, mean_dEdx, particleForwardsChiSquared, particleBackwardsChiSquared);
    fitResult = finalFitResult;
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FitHitChargeVector(HitChargeVector &hitChargeVector1, HitChargeVector &hitChargeVector2, TrackDirectionTool::FitResult &fitResult1, TrackDirectionTool::FitResult &fitResult2, int numberHitsToConsider)
{
    this->FitHitChargeVector(hitChargeVector1, fitResult1, numberHitsToConsider);
    this->FitHitChargeVector(hitChargeVector2, fitResult2, numberHitsToConsider);
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetGlobalMinuitPreliminaries(const HitChargeVector &hitChargeVector)
{
    this->ClearGlobalVariables();

    for (const HitCharge &hitCharge : hitChargeVector)
        pMinuitVector_tool->push_back(hitCharge);

    for (const HitCharge &hitCharge : *pMinuitVector_tool)
    {
        globalTotalHitWidth_tool += hitCharge.GetHitWidth();
        globalTotalCharge_tool += hitCharge.GetCharge();
    }

    this->GetTrackLength(hitChargeVector, globalTrackLength_tool);
}

//---------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::PerformFits(const HitChargeVector &hitChargeVector, HitChargeVector &forwardsFitPoints, HitChargeVector &backwardsFitPoints, int numberHitsToConsider, float &forwardsChiSquared, float &backwardsChiSquared, int &fitStatus1, int &fitStatus2)
{
    this->SetGlobalMinuitPreliminaries(hitChargeVector);

    //---------------------------------------------------------------------------------------------------
    //Forwards Fit

    double particleMass(105.7), maxScale(globalMuonLookupTable_tool.GetMaxRange()/globalTrackLength_tool);
    LookupTable lookupTable = globalMuonLookupTable_tool;

    const int nParameters = 3;
    const std::string parName[nParameters]   = {"ENDENERGY", "SCALE", "EXTRA"};
    const double vstart[nParameters] = {2.1, 1.0, 1.0};
    const double step[nParameters] = {1.e-1, 1.e-1, 1.e-1};
    const double lowphysbound[nParameters] = {2.0, 0.01, 0.1};
    const double highphysbound[nParameters] = {1.0e3, maxScale, 1.0e1};

    int ierflg(0);

    TMinuit *pMinuit = new TMinuit(nParameters);
    pMinuit->SetPrintLevel(-1);
    pMinuit->SetFCN(GetForwardsChiSquared_Tool);

    for (int j = 0 ; j < nParameters ; ++j)
        pMinuit->mnparm(j, parName[j].c_str(), vstart[j], step[j], lowphysbound[j], highphysbound[j], ierflg);

    double arglist[2];
    arglist[0] = 50000;
    arglist[1] = 2;
    pMinuit->mnexcm("MIGRAD", arglist, 1, fitStatus1);

    double outpar[nParameters], err[nParameters];

    for (int k = 0; k < nParameters; k++)
        pMinuit->GetParameter(k, outpar[k], err[k]);

    delete pMinuit;

    //---------------------------------------------------------------------------------
    //Backwards Fit

    const int nParameters2 = 3;
    const std::string parName2[nParameters2]   = {"ENDENERGY", "SCALE", "EXTRA"};
    const double vstart2[nParameters2] = {2.1, 1.0, 1.0};
    const double step2[nParameters2] = {1.e-1, 1.e-1, 1.e-1};
    const double lowphysbound2[nParameters2] = {2.0, 0.01, 0.1};
    const double highphysbound2[nParameters2] = {1.0e3, maxScale, 1.0e1};

    int ierflg2(0);

    TMinuit *pMinuit2 = new TMinuit(nParameters2);
    pMinuit2->SetPrintLevel(-1);
    pMinuit2->SetFCN(GetBackwardsChiSquared_Tool);

    for (int j = 0 ; j < nParameters2 ; ++j)
        pMinuit2->mnparm(j, parName2[j].c_str(), vstart2[j], step2[j], lowphysbound2[j], highphysbound2[j], ierflg2);

    double arglist2[2];
    arglist2[0] = 50000;
    arglist2[1] = 2;
    pMinuit2->mnexcm("MIGRAD", arglist2, 1, fitStatus2);

    double outpar2[nParameters2], err2[nParameters2];

    for (int k = 0; k < nParameters2; k++)
        pMinuit2->GetParameter(k, outpar2[k], err2[k]);

    delete pMinuit2;

    //--------------------------------------------------------------------------

    double f_Ee(outpar[0]), f_L(outpar[1] * globalTrackLength_tool);
    double f_Le(GetLengthfromEnergy_Tool(lookupTable, f_Ee));
    double f_Ls = f_Le - f_L;

    double f_Es = GetEnergyfromLength_Tool(lookupTable, f_Ls);
    double f_deltaE = f_Es - f_Ee;

    double f_alpha = f_deltaE/globalTotalCharge_tool;
    double f_beta = f_L/globalTotalHitWidth_tool;

    double b_Ee(outpar2[0]), b_L(outpar2[1] * globalTrackLength_tool);
    double b_Le(GetLengthfromEnergy_Tool(lookupTable, b_Ee));
    double b_Ls = b_Le - b_L;

    double b_Es = GetEnergyfromLength_Tool(lookupTable, b_Ls);
    double b_deltaE = b_Es - b_Ee;

    double b_alpha = b_deltaE/globalTotalCharge_tool;
    double b_beta = b_L/globalTotalHitWidth_tool;

    //--------------------------------------------------------------------------

    int nImpureHits(0);

    for (int vectorPosition = 0; vectorPosition < hitChargeVector.size(); vectorPosition++)
    {
        HitCharge hitCharge(hitChargeVector.at(vectorPosition));

        bool isPure(false);

        if (!isPure)
            nImpureHits++;

        double f_L_i = f_Ls + (outpar[1] * hitCharge.GetLongitudinalPosition());
        double f_E_i = GetEnergyfromLength_Tool(lookupTable, f_L_i);
        double f_dEdx_2D = outpar[2] * (f_beta/f_alpha) * BetheBloch_Tool(f_E_i, particleMass);

        double b_L_i = b_Ls + (outpar2[1] * (globalTrackLength_tool - hitCharge.GetLongitudinalPosition()));
        double b_E_i = GetEnergyfromLength_Tool(lookupTable, b_L_i);
        double b_dEdx_2D = outpar2[2] * (b_beta/b_alpha) * BetheBloch_Tool(b_E_i, particleMass);

        double Q_fit_f(f_dEdx_2D * hitCharge.GetHitWidth());
        double Q_fit_b(b_dEdx_2D * hitCharge.GetHitWidth());

        float forwardsDelta(hitCharge.GetQoverX() - f_dEdx_2D), backwardsDelta(hitCharge.GetQoverX() - b_dEdx_2D);

        float f_sigma(std::sqrt((0.00419133 * f_dEdx_2D * f_dEdx_2D) + (0.00967141 * f_dEdx_2D))); //70%
        float b_sigma(std::sqrt((0.00419133 * b_dEdx_2D * b_dEdx_2D) + (0.00967141 * b_dEdx_2D))); //70%

        float lp(hitCharge.GetLongitudinalPosition()), hw(hitCharge.GetHitWidth());
        float f_Q_fit_f(Q_fit_f), f_Q_fit_b(Q_fit_b);
        HitCharge forwardsRecoHitCharge(hitCharge.GetCaloHit(), lp, hw, f_Q_fit_f, f_sigma);
        forwardsFitPoints.push_back(forwardsRecoHitCharge);
        HitCharge backwardsRecoHitCharge(hitCharge.GetCaloHit(), lp, hw, f_Q_fit_b, b_sigma);
        backwardsFitPoints.push_back(backwardsRecoHitCharge);

        float forwardsHitChisquared((forwardsDelta * forwardsDelta)/(f_sigma * f_sigma));
        float backwardsHitChisquared((backwardsDelta * backwardsDelta)/(b_sigma * b_sigma));

        if (!((pMinuitVector_tool->size() >= 2 * numberHitsToConsider) && vectorPosition > numberHitsToConsider && vectorPosition < pMinuitVector_tool->size() - numberHitsToConsider))
        {
            forwardsChiSquared += forwardsHitChisquared;
            backwardsChiSquared += backwardsHitChisquared;
        }
    }

    std::sort(forwardsFitPoints.begin(), forwardsFitPoints.end(), SortHitChargeVectorByRL);
    std::sort(backwardsFitPoints.begin(), backwardsFitPoints.end(), SortHitChargeVectorByRL);

    this->ClearGlobalVariables();
}

//---------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetCalorimetricDirection(const Cluster* pTargetClusterW, FitResult &finalFitResult)
{
    HitChargeVector hitChargeVector;
    this->FillHitChargeVector(pTargetClusterW, hitChargeVector);

    HitChargeVector filteredHitChargeVector;
    this->TrackInnerFilter(hitChargeVector, filteredHitChargeVector);

    this->SimpleTrackEndFilter(filteredHitChargeVector);

    if (m_enableFragmentRemoval)
        this->TrackEndFilter(filteredHitChargeVector);

    std::sort(filteredHitChargeVector.begin(), filteredHitChargeVector.end(), SortHitChargeVectorByRL);

    //Hypothesis one: clean single particle
    FitResult initialFitResult;
    this->FitHitChargeVector(filteredHitChargeVector, initialFitResult);

    bool likelyForwards(initialFitResult.GetDirectionEstimate() == 1 && filteredHitChargeVector.size() >= 400 && initialFitResult.GetForwardsChiSquared()/initialFitResult.GetNHits() <= 1.25);
    bool likelyBackwards(initialFitResult.GetDirectionEstimate() == 0 && filteredHitChargeVector.size() <= 200 && initialFitResult.GetBackwardsChiSquared()/initialFitResult.GetNHits() <= 1.25);

    if (!(likelyForwards || likelyBackwards))
    {
        //Hypothesis two: clean double particle
        FitResult backwardsSplitResult, forwardsSplitResult;
        bool splitApplied(false);

        if (m_enableSplitting)
        {
            this->ParticleSplitting(pTargetClusterW, filteredHitChargeVector, backwardsSplitResult, forwardsSplitResult, splitApplied);
        }

        if (!splitApplied)
        {
            if (m_enableFragmentRemoval)
            {
                //Fragment Removal
                HitChargeVector fragmentlessHitChargeVector;
                this->FragmentRemoval(filteredHitChargeVector, fragmentlessHitChargeVector);

                //Hypothesis three: fragmentless single particle
                FitResult fragmentRemovalFitResult;
                this->FitHitChargeVector(fragmentlessHitChargeVector, fragmentRemovalFitResult);

                bool likelyCorrectFragmentRemoval(initialFitResult.GetDirectionEstimate() != fragmentRemovalFitResult.GetDirectionEstimate() && initialFitResult.GetMinChiSquaredPerHit() - fragmentRemovalFitResult.GetMinChiSquaredPerHit() >= 2.0);

                if (!likelyCorrectFragmentRemoval)
                {
                    std::cout << "Applied Hypothesis #1: After All Tests (Post Desperation Single Particle)" << std::endl;
                    finalFitResult = initialFitResult;
                }
                else
                {
                    std::cout << "Applied Hypothesis #3 (Fragment Removal)" << std::endl;
                    finalFitResult = fragmentRemovalFitResult;
                }
            }
            else
            {
                std::cout << "Applied Hypothesis #1: After All Tests (No Cuts Met: Single Particle)" << std::endl;
                finalFitResult = initialFitResult;
            }
        }
        else
        {
            std::cout << "Applied Hypothesis #2 (First Splitting)" << std::endl;
            finalFitResult = (backwardsSplitResult.GetMeanQoverX() < forwardsSplitResult.GetMeanQoverX() ? backwardsSplitResult: forwardsSplitResult); //TODO
        }
    }
    else
    {
        std::cout << "Applied Hypothesis #1 (Single Clean Particle)" << std::endl;
        finalFitResult = initialFitResult;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::AddToSlidingFitCache(const Cluster *const pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingFitResult &TrackDirectionTool::GetCachedSlidingFitResult(const Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TidyUp()
{
    m_slidingFitResultMap.clear();

    globalTrackLength_tool = (0.f);
    globalTotalHitWidth_tool = (0.f);

    globalBinWidth_tool = (0.1);
    globalInitialEnergy_tool = (2000.0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::ClearGlobalVariables()
{
    pMinuitVector_tool->clear();
    globalTotalCharge_tool = 0.f;
    globalTrackLength_tool = 0.f;
    globalTotalHitWidth_tool = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortHitChargeVectorByRL(HitCharge &hitCharge1, HitCharge &hitCharge2)
{
    return hitCharge1.GetLongitudinalPosition() < hitCharge2.GetLongitudinalPosition();
}

//----------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortHitChargeVectorByQoverX(HitCharge &hitCharge1, HitCharge &hitCharge2)
{
    return hitCharge1.GetQoverX() < hitCharge2.GetQoverX();
}

//----------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortByDistanceToNN(HitCharge &hitCharge1, HitCharge &hitCharge2)
{
    return hitCharge1.GetDistanceToNN() < hitCharge2.GetDistanceToNN();
}

//----------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortJumpVector(JumpObject &jumpObject1, JumpObject &jumpObject2)
{
    return jumpObject1.GetJumpValue() > jumpObject2.GetJumpValue();
}

//----------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackDirectionTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NumberTrackEndHits", m_numberTrackEndHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableFragmentRemoval", m_enableFragmentRemoval));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableSplitting", m_enableSplitting));

    return STATUS_CODE_SUCCESS;
}

//----------------------------------------------------------------------------------------------------------------------------------

} // namespac lar_content