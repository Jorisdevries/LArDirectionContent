/**
 *  @file   ExampleContent/src/ExampleAlgorithms/DirectionAnalysisAlgorithm.cc
 * 
 *  @brief  Implementation of the access lists algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArDirection/DirectionAnalysisAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "LArSpaceChargeHelper.h"

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

std::string clusterTreeName("Cluster"), pfoTreeName("PFO"), hitTreeName("Hit"), vertexTreeName("Vertex"), splittingTreeName("Splitting");

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionAnalysisAlgorithm::DirectionAnalysisAlgorithm():
    m_targetParticlePDG(13),
    m_particleContained(true),
    m_cosmic(true),
    m_data(false),
    m_drawFit(false),
    m_fileIdentifier(0),
    m_eventNumber(-1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionAnalysisAlgorithm::~DirectionAnalysisAlgorithm()
{
    if (m_writeToTree)
    {    
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), clusterTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), pfoTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), hitTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), vertexTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), splittingTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DirectionAnalysisAlgorithm::Run()
{
    m_eventNumber++;

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputHitListName, pCaloHitList));

    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));

    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));

    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    pandora::MCParticleVector mcParticleVector(pMCParticleList->begin(), pMCParticleList->end());
    pandora::ClusterVector clusterVector(pClusterList->begin(), pClusterList->end());
    pandora::PfoVector pfoVector(pPfoList->begin(), pPfoList->end());
    pandora::VertexVector vertexVector(pVertexList->begin(), pVertexList->end());

    //pandora::CaloHitVector hitVector(pCaloHitList->begin(), pCaloHitList->end());
    //for (auto pCaloHit : hitVector)
    //    std::cout << "Hit position: (" << pCaloHit->GetPositionVector().GetX() << "," << pCaloHit->GetPositionVector().GetY() << "," << pCaloHit->GetPositionVector().GetZ() << ")" << std::endl;

    this->WritePfoInformation(pfoVector);
    this->WriteClusterAndHitInformation(clusterVector);

    if (this->CheckEventType(pMCParticleList, pCaloHitList, pfoVector, 1, 0, 1))
        this->WriteVertexInformation(pMCParticleList, pCaloHitList, vertexVector, pfoVector);

    if (this->CheckEventType(pMCParticleList, pCaloHitList, pfoVector, 1, 1, 0)) //0 argument means do not check this count
        this->WriteSplittingInformation(pMCParticleList, pCaloHitList, vertexVector, pfoVector);
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::CheckEventType(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, pandora::PfoVector &pfoVector, int targetNumberMuons, int targetNumberProtons, int targetNumberPfos)
{
    pandora::MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    if (trueNeutrinos.size() != 1)
        return false;

    if (LArMCParticleHelper::GetNuanceCode((*(trueNeutrinos.begin()))) != 1001) 
        return false;

    LArMCParticleHelper::MCContributionMap selectedMCParticles;
    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, selectedMCParticles);

    int nMuons(0), nProtons(0), nOthers(0);

    for (auto mcParticle : selectedMCParticles)
    {
        if (mcParticle.first->GetParticleId() == 13)
            nMuons++;
        else if (mcParticle.first->GetParticleId() == 2212)
            nProtons++;
        else
            nOthers++;
    }

    int nPrimaryPfos(0);

    for (auto pPfo : pfoVector)
    {
        if (!LArPfoHelper::IsFinalState(pPfo))
            continue;

        nPrimaryPfos++;
    }

    if (nMuons != targetNumberMuons || (targetNumberProtons != 0 && nProtons != targetNumberProtons) || (targetNumberPfos != 0 && nPrimaryPfos != targetNumberPfos) || nOthers != 0)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WriteVertexInformation(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, pandora::VertexVector &vertexVector, pandora::PfoVector &pfoVector)
{
    //Reweight vertex scores so that the lowest is 1 and the highest is 5
    this->NormaliseVertexScores(vertexVector);

    float vertexDR(this->GetVertexDR(pMCParticleList, pCaloHitList, vertexVector));
    float minVertexDR(this->GetMinVertexDR(pMCParticleList, pCaloHitList, vertexVector));

    float longestPfoLength(0.f);
    pandora::PfoVector targetPfoVector(pfoVector); 

    for (auto pPfo : pfoVector)
    {   
        if (!LArPfoHelper::IsFinalState(pPfo))
            continue;

        const ClusterList clusterList(pPfo->GetClusterList());

        for (auto pCluster : clusterList)
        {   
            if (LArClusterHelper::GetLength(pCluster) > longestPfoLength)
            {   
                longestPfoLength = LArClusterHelper::GetLength(pCluster);
                targetPfoVector.clear();
                targetPfoVector.push_back(pPfo);
            }   
        }   
    }

    if (targetPfoVector.size() == 0)
        return;

    const ParticleFlowObject *const pTargetPfo(*(targetPfoVector.begin()));
    TrackDirectionTool::DirectionFitObject directionFit = m_pTrackDirectionTool->GetPfoDirection(pTargetPfo);

    pandora::VertexVector vertexVectorReweighted(vertexVector);

    //Linearly reweight vertex score by direction probability
    for (auto pVertex : vertexVectorReweighted)
    {
        float fractionalCloseness(1.0 - ((directionFit.GetBeginpoint() - pVertex->GetPosition()).GetMagnitude())/((directionFit.GetBeginpoint() - directionFit.GetEndpoint()).GetMagnitude()));
        float newScore(pVertex->GetScore() * fractionalCloseness * directionFit.GetProbability());
        pVertex->SetScore(newScore);
    }

    float afterDirectionVertexDR(this->GetVertexDR(pMCParticleList, pCaloHitList, vertexVectorReweighted));

    if (m_drawFit)
    {
        //directionFit.DrawFit();
        //PANDORA_MONITORING_API(Pause(this->GetPandora()));
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "VertexDR", vertexDR));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "MinVertexDR", minVertexDR));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "AfterDirectionVertexDR", afterDirectionVertexDR));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "FileIdentifier", m_fileIdentifier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "EventNumber", m_eventNumber));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), vertexTreeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WriteSplittingInformation(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, pandora::VertexVector &vertexVector, pandora::PfoVector &pfoVector)
{
    float longestClusterLength(0.f);
    pandora::ClusterVector targetClusterVector; 
    pandora::PfoVector targetPfoVector; 

    for (auto pPfo : pfoVector)
    {   
        try
        {
            if (!this->IsGoodPfo(pPfo))
                continue;

            const Cluster *const pCluster = this->GetTargetClusterFromPFO(pPfo); 

            if (LArClusterHelper::GetLength(pCluster) > longestClusterLength)
            {   
                longestClusterLength = LArClusterHelper::GetLength(pCluster);
                targetClusterVector.clear();
                targetClusterVector.push_back(pCluster);
                targetPfoVector.clear();
                targetPfoVector.push_back(pPfo);
            }   
        }
        catch (...)
        {
            continue;
        }
    }

    if (targetClusterVector.size() == 0 || targetPfoVector.size() == 0)
        return;

    const Cluster* const pTargetCluster(*(targetClusterVector.begin()));
    //const ParticleFlowObject *const pTargetPfo(*(targetPfoVector.begin()));

    TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetClusterDirection(pTargetCluster);

    float splitPosition(fitResult.GetSplitObject().GetSplitPosition());

    //const pandora::CaloHit* leftHit(fitResult.GetForwardsFitCharges().back().GetLongitudinalPosition() - fitResult.GetBackwardsFitCharges().front().GetLongitudinalPosition() < fitResult.GetForwardsFitCharges().front().GetLongitudinalPosition() - fitResult.GetBackwardsFitCharges().back().GetLongitudinalPosition() ? fitResult.GetForwardsFitCharges().back().GetCaloHit() : fitResult.GetForwardsFitCharges().front().GetCaloHit());
    //const pandora::CaloHit* rightHit(fitResult.GetForwardsFitCharges().back().GetLongitudinalPosition() - fitResult.GetBackwardsFitCharges().front().GetLongitudinalPosition() < fitResult.GetForwardsFitCharges().front().GetLongitudinalPosition() - fitResult.GetBackwardsFitCharges().back().GetLongitudinalPosition() ? fitResult.GetForwardsFitCharges().front().GetCaloHit() : fitResult.GetForwardsFitCharges().back().GetCaloHit());

    const pandora::CaloHit* leftHit(((fitResult.GetBackwardsFitCharges().back())).GetCaloHit());
    const pandora::CaloHit* rightHit(((fitResult.GetForwardsFitCharges().front())).GetCaloHit());

    const float targetX((leftHit->GetPositionVector().GetX() + rightHit->GetPositionVector().GetX())/2);
    const float lowerZ(leftHit->GetPositionVector().GetZ() < rightHit->GetPositionVector().GetZ() ? leftHit->GetPositionVector().GetZ() : rightHit->GetPositionVector().GetZ());
    const float upperZ(leftHit->GetPositionVector().GetZ() > rightHit->GetPositionVector().GetZ() ? leftHit->GetPositionVector().GetZ() : rightHit->GetPositionVector().GetZ());

    std::cout << "lowerZ: " << lowerZ << std::endl;
    std::cout << "upperZ: " << upperZ << std::endl;
    
    //get closest W calohit
    float closestDistance(1e6);
    const pandora::CaloHit* pTargetCaloHitW(NULL);

    for (auto pCaloHit : *pCaloHitList)
    {
        if (std::abs(pCaloHit->GetPositionVector().GetX() - targetX) < closestDistance && pCaloHit->GetPositionVector().GetZ() > lowerZ && pCaloHit->GetPositionVector().GetZ() < upperZ)
        {
            closestDistance = std::abs(pCaloHit->GetPositionVector().GetX() - targetX); 
            pTargetCaloHitW = pCaloHit; 
        }
    }

    //for (auto pCaloHit : *pCaloHitList)

    /*
    //loop over clusters in pfo to find matching hit in U or V
    ClusterList ClustersV;
    LArPfoHelper::GetClusters(pTargetPfo, TPC_VIEW_V, ClustersV);
    const pandora::CaloHit* pTargetCaloHitV(NULL);
    float closestHitDistance(1e6);

    for (auto pClusterV : ClustersV)
    {
        OrderedCaloHitList orderedCaloHitList(pClusterV->GetOrderedCaloHitList());
        CaloHitList caloHitList;
        orderedCaloHitList.FillCaloHitList(caloHitList);

        for (auto pCaloHitV : caloHitList)
        {
            float deltaX((pCaloHitV->GetPositionVector().GetX() != 0 && pTargetCaloHitW->GetPositionVector().GetX() != 0) ? (std::abs(pCaloHitV->GetPositionVector().GetX() - pTargetCaloHitW->GetPositionVector().GetX())) : 0);
            float deltaY((pCaloHitV->GetPositionVector().GetY() != 0 && pTargetCaloHitW->GetPositionVector().GetY() != 0) ? (std::abs(pCaloHitV->GetPositionVector().GetY() - pTargetCaloHitW->GetPositionVector().GetY())) : 0);
            float deltaZ((pCaloHitV->GetPositionVector().GetZ() != 0 && pTargetCaloHitW->GetPositionVector().GetZ() != 0) ? (std::abs(pCaloHitV->GetPositionVector().GetZ() - pTargetCaloHitW->GetPositionVector().GetZ())) : 0);

            float distance(std::sqrt((deltaX*deltaX) + (deltaY*deltaY) + (deltaZ*deltaZ)));

            if (distance < closestHitDistance)
            {
                closestHitDistance = distance;
                pTargetCaloHitV = pCaloHitV;
            }
        }
    }

    CartesianVector vertexPosition(0.f, 0.f, 0.f);
    float positionChiSquared(0.f);
    LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, pTargetCaloHitV->GetPositionVector(), pTargetCaloHitW->GetPositionVector(), vertexPosition, positionChiSquared);
    */

    //Now find the closest 3D hit
    pandora::CaloHitVector hitVector(pCaloHitList->begin(), pCaloHitList->end());
    const pandora::CaloHit* pTarget3DCaloHit(NULL);
    float closestDistance3D(1e6);

    for (auto pCaloHit : hitVector)
    {
        /*
        if ((pCaloHit->GetPositionVector() - vertexPosition).GetMagnitude() < closestDistance3D)
        {
            closestDistance3D = (pCaloHit->GetPositionVector() - vertexPosition).GetMagnitude();
            pTarget3DCaloHit = pCaloHit;
        }
        */

        if (pCaloHit->GetPositionVector().GetX() == 0 || pCaloHit->GetPositionVector().GetY() == 0 || pCaloHit->GetPositionVector().GetZ() == 0)
            continue;
    
            float deltaX((pCaloHit->GetPositionVector().GetX() != 0 && pTargetCaloHitW->GetPositionVector().GetX() != 0) ? (std::abs(pCaloHit->GetPositionVector().GetX() - pTargetCaloHitW->GetPositionVector().GetX())) : 0);
            float deltaZ((pCaloHit->GetPositionVector().GetZ() != 0 && pTargetCaloHitW->GetPositionVector().GetZ() != 0) ? (std::abs(pCaloHit->GetPositionVector().GetZ() - pTargetCaloHitW->GetPositionVector().GetZ())) : 0);

        float distance(std::sqrt((deltaX*deltaX) + (deltaZ*deltaZ)));
       
        if (distance < closestDistance3D) 
        {
            closestDistance3D = distance; 
            pTarget3DCaloHit = pCaloHit;
        }
    }

    //Test spacecharge helper
    LArSpaceChargeHelper::Configure("SCEoffsets_MicroBooNE_E273.root");
    CartesianVector vertexPosition(pTarget3DCaloHit->GetPositionVector());
    
    std::cout << "Original vertex position: (" << vertexPosition.GetX() << "," << vertexPosition.GetY() << "," << vertexPosition.GetZ() << ")" << std::endl;

    CartesianVector correctedPosition = LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(vertexPosition);

    std::cout << "Spacecharge corrected vertex position: (" << correctedPosition.GetX() << "," << correctedPosition.GetY() << "," << correctedPosition.GetZ() << ")" << std::endl;

    //Create vertex
    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = pTarget3DCaloHit->GetPositionVector();
    parameters.m_vertexLabel = VERTEX_INTERACTION; 
    parameters.m_vertexType = VERTEX_3D; 

    const Vertex *pNewVertex(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));

    float vertexDR(this->GetVertexDR(pMCParticleList, pCaloHitList, vertexVector));

    CartesianVector trueVertexPosition(0.f, 0.f, 0.f);

    for (auto pMCParticle : *pMCParticleList)
    {   
        if (pMCParticle->GetParticleId() == 13 && LArMCParticleHelper::IsPrimary(pMCParticle) && pMCParticle->GetVertex().GetY() != 0)
            trueVertexPosition = pMCParticle->GetVertex();
    }

    /*
    for (auto hitCharge : fitResult.GetBackwardsFitCharges())
    {
        std::cout << "L: " << hitCharge.GetLongitudinalPosition() << std::endl;
        std::cout << "Z: " << hitCharge.GetCaloHit()->GetPositionVector().GetZ() << std::endl;
        std::cout << "------------" << std::endl;
    }

    for (auto hitCharge : fitResult.GetForwardsFitCharges())
    {
        std::cout << "L: " << hitCharge.GetLongitudinalPosition() << std::endl;
        std::cout << "Z: " << hitCharge.GetCaloHit()->GetPositionVector().GetZ() << std::endl;
        std::cout << "------------" << std::endl;
    }
    */

    std::cout << "Split position: " << splitPosition << std::endl;
    std::cout << "True vertex position: (" << trueVertexPosition.GetX() << "," << trueVertexPosition.GetY() << "," << trueVertexPosition.GetZ() << ")" << std::endl;
    std::cout << "New vertex position: (" << pTarget3DCaloHit->GetPositionVector().GetX() << "," << pTarget3DCaloHit->GetPositionVector().GetY() << "," << pTarget3DCaloHit->GetPositionVector().GetZ() << ")" << std::endl;

    //bool extremeChange(std::abs((trueVertexPosition - pNewVertex->GetPosition()).GetMagnitude() - vertexDR) >= 0.9 * LArClusterHelper::GetLength(pTargetCluster));
    float afterSplitVertexDR(fitResult.GetSplitObject().GetSplitApplied() ? (trueVertexPosition - pNewVertex->GetPosition()).GetMagnitude() : vertexDR);

    bool isClusterTwoParticles(false);
    this->IsClusterTwoParticles(pTargetCluster, fitResult.GetForwardsFitCharges(), fitResult.GetBackwardsFitCharges(), isClusterTwoParticles);
    int mcIsClusterTwoParticles(isClusterTwoParticles? 1 : 0);

    std::cout << "Before split vertex DR: " << vertexDR << std::endl;
    std::cout << "After split vertex DR: " << afterSplitVertexDR << std::endl;
    
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), splittingTreeName.c_str(), "MCIsClusterTwoParticles", mcIsClusterTwoParticles));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), splittingTreeName.c_str(), "SplitPosition", splitPosition));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), splittingTreeName.c_str(), "VertexDR", vertexDR));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), splittingTreeName.c_str(), "AfterSplitVertexDR", afterSplitVertexDR));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), splittingTreeName.c_str(), "FileIdentifier", m_fileIdentifier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), splittingTreeName.c_str(), "EventNumber", m_eventNumber));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), splittingTreeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DirectionAnalysisAlgorithm::GetVertexDR(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, pandora::VertexVector &vertexVector)
{
    CartesianVector trueVertexPosition(0.f, 0.f, 0.f);

    for (auto pMCParticle : *pMCParticleList)
    {
        if (pMCParticle->GetParticleId() == 13 && LArMCParticleHelper::IsPrimary(pMCParticle) && pMCParticle->GetVertex().GetY() != 0)
            trueVertexPosition = pMCParticle->GetVertex();
    }

    LArMCParticleHelper::MCContributionMap selectedMCParticles;
    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, selectedMCParticles);

    float maxScore(-1000.f), vertexDR(0.f);
    for (const pandora::Vertex *const pVertex : vertexVector)
    {
        if (pVertex->GetScore() > maxScore)
        {
            maxScore = pVertex->GetScore();
            vertexDR = (pVertex->GetPosition() - trueVertexPosition).GetMagnitude();
        }
    }

    return vertexDR;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DirectionAnalysisAlgorithm::GetMinVertexDR(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, pandora::VertexVector &vertexVector)
{
    CartesianVector trueVertexPosition(0.f, 0.f, 0.f);

    for (auto pMCParticle : *pMCParticleList)
    {
        if (pMCParticle->GetParticleId() == 13 && LArMCParticleHelper::IsPrimary(pMCParticle) && pMCParticle->GetVertex().GetY() != 0)
            trueVertexPosition = pMCParticle->GetVertex();
    }

    LArMCParticleHelper::MCContributionMap selectedMCParticles;
    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, selectedMCParticles);

    float smallestDR(1e6);
    for (const pandora::Vertex *const pVertex : vertexVector)
    {
        if ((pVertex->GetPosition() - trueVertexPosition).GetMagnitude() < smallestDR)
            smallestDR = (pVertex->GetPosition() - trueVertexPosition).GetMagnitude();
    }

    return smallestDR;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::NormaliseVertexScores(pandora::VertexVector &vertexVector)
{
    float lowestScore(1e6), highestScore(0.f);

    for (const pandora::Vertex *const pVertex : vertexVector)
    {
        if (pVertex->GetScore() < lowestScore)
            lowestScore = pVertex->GetScore();
        if (pVertex->GetScore() > highestScore)
            highestScore = pVertex->GetScore();
    }

    float newMinimum(1.0), newMaxiumum(2.0);

    for (const pandora::Vertex *const pVertex : vertexVector)
    {
        float newScore(((newMaxiumum - newMinimum) * (pVertex->GetScore() - lowestScore))/(highestScore - lowestScore) + newMinimum);
        pVertex->SetScore(newScore);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WritePfoInformation(pandora::PfoVector &pfoVector)
{
    for (const pandora::ParticleFlowObject *const pPfo : pfoVector)
    {
        try
        {
            if (!this->IsGoodPfo(pPfo))
                continue;

            const Cluster *const pCluster = this->GetTargetClusterFromPFO(pPfo); 
            TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetPfoDirection(pPfo);
            this->WriteToTree(pCluster, fitResult, pfoTreeName);

            if (m_drawFit)
            {
                std::cout << "Split applied: " << fitResult.GetSplitObject().GetSplitApplied() << std::endl;
                fitResult.DrawFit();
                PANDORA_MONITORING_API(Pause(this->GetPandora()));
            }
        }

        catch (...)
        {
            std::cout << "Skipping..." << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::IsGoodPfo(const pandora::ParticleFlowObject* pPfo)
{
    const MCParticle* const pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));

    if (!LArPfoHelper::IsFinalState(pPfo))
        return false;

    if ((m_particleContained && !m_data && IsParticleContained(pMCParticle)) || (m_particleContained && m_data && !IsRecoParticleContained(pPfo)))
        return false;

    if (pMCParticle->GetParticleId() != m_targetParticlePDG)
        return false;

    if (m_cosmic && !LArMCParticleHelper::IsCosmicRay(pMCParticle))
        return false;

    if (!m_cosmic && LArMCParticleHelper::IsCosmicRay(pMCParticle)) 
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WriteClusterAndHitInformation(pandora::ClusterVector &clusterVector)
{
    for (const pandora::Cluster* const pCluster : clusterVector)
    {
        try
        {
            if (LArClusterHelper::GetClusterHitType(pCluster) != TPC_VIEW_W)
                continue;

            const MCParticle* const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

            if ((m_particleContained && !m_data && IsParticleContained(pMCParticle)) || (m_particleContained && m_data && !IsRecoParticleContained(pCluster)))
                continue;

            if (pMCParticle->GetParticleId() != m_targetParticlePDG)
                continue;

            if (m_cosmic && !LArMCParticleHelper::IsCosmicRay(pMCParticle))
                continue;

            if (!m_cosmic && LArMCParticleHelper::IsCosmicRay(pMCParticle)) 
                continue;

            TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetClusterDirection(pCluster);
            this->WriteToTree(pCluster, fitResult, clusterTreeName);

            for (TrackDirectionTool::HitCharge &hitCharge : fitResult.GetHitChargeVector())
                this->WriteHitToTree(hitCharge, hitTreeName, fitResult.GetMCDirection());

            if (m_drawFit)
            {
                //fitResult.DrawFit();
                //PANDORA_MONITORING_API(Pause(this->GetPandora()));
            }
        }

        catch (...)
        {
            //std::cout << "Skipping..." << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster* DirectionAnalysisAlgorithm::GetTargetClusterFromPFO(const ParticleFlowObject* pPfo)
{
    HitType hitType(TPC_VIEW_W);
    ClusterList clusterListW;
    LArPfoHelper::GetClusters(pPfo, hitType, clusterListW);

    if (clusterListW.size() == 0)
    {    
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }    

    float currentLength(std::numeric_limits<float>::min());
    ClusterList longestClusterList;

    for (const Cluster *const pCluster : clusterListW)
    {    
        if (LArClusterHelper::GetLength(pCluster) > currentLength)
        {    
            currentLength = LArClusterHelper::GetLength(pCluster);
            longestClusterList.clear();
            longestClusterList.push_back(pCluster);
        }    
    }    

    const Cluster *const pCluster(*(longestClusterList.begin()));
    return pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::IsClusterTwoParticles(const Cluster *const pCluster, TrackDirectionTool::HitChargeVector forwardsFitCharges, TrackDirectionTool::HitChargeVector backwardsFitCharges, bool &isTwoParticles)
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    int nContributingParticles(0), contributionThreshold(10);
    int nSecondaryHits(0);
    float chargeContributionThreshold(0.5);
    std::map<int, int> primaryPDGToContributionMap;

    const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

    for (const CaloHit* pCaloHit : caloHitList)
    {    
        MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
        for (MCParticleWeightMap::const_iterator mapIter = mcParticleWeightMap.begin(), mapIterEnd = mcParticleWeightMap.end(); mapIter != mapIterEnd; ++mapIter)
        {    
            int primaryPDG(mapIter->first->GetParticleId());
            if (primaryPDG == 11 && ((pCaloHit->GetPositionVector().GetZ() > pMCParticle->GetVertex().GetZ() && pCaloHit->GetPositionVector().GetZ() < pMCParticle->GetEndpoint().GetZ()) || (pCaloHit->GetPositionVector().GetZ() < pMCParticle->GetVertex().GetZ() && pCaloHit->GetPositionVector().GetZ() > pMCParticle->GetEndpoint().GetZ())))
                primaryPDG = 13;

            float contribution(mapIter->second);

            if (primaryPDGToContributionMap.find(primaryPDG) == primaryPDGToContributionMap.end())
            {    
                if (contribution > chargeContributionThreshold)
                    primaryPDGToContributionMap[primaryPDG] = 1; 
            }    
            else 
            {    
                if (contribution > chargeContributionThreshold)
                    primaryPDGToContributionMap.at(primaryPDG)++;
            }    
        }    
    }    

    for (auto &entry : primaryPDGToContributionMap)
    {    
        if (entry.first != m_targetParticlePDG)
            nSecondaryHits += entry.second;

        if (entry.second >= contributionThreshold)
            nContributingParticles++;
    } 

    float chargeContributionThreshold2(0.25);
    int backwardsNonPrimaryHits(0), backwardsPrimaryHits(0);
    for (TrackDirectionTool::HitCharge &hitCharge : backwardsFitCharges)
    {    
        const pandora::CaloHit* pCaloHit(hitCharge.GetCaloHit());
        MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

        for (MCParticleWeightMap::const_iterator mapIter = mcParticleWeightMap.begin(), mapIterEnd = mcParticleWeightMap.end(); mapIter != mapIterEnd; ++mapIter)
        {    
            int primaryPDG(mapIter->first->GetParticleId());

            float contribution(mapIter->second);
            if (contribution >= chargeContributionThreshold2 && primaryPDG != m_targetParticlePDG)
                backwardsNonPrimaryHits++;
            else if (contribution >= chargeContributionThreshold2 && primaryPDG == m_targetParticlePDG)
                backwardsPrimaryHits++;
        }    
    }    

    int forwardsNonPrimaryHits(0), forwardsPrimaryHits(0);
    for (TrackDirectionTool::HitCharge &hitCharge : forwardsFitCharges)
    {    
        const pandora::CaloHit* pCaloHit(hitCharge.GetCaloHit());
        MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

        for (MCParticleWeightMap::const_iterator mapIter = mcParticleWeightMap.begin(), mapIterEnd = mcParticleWeightMap.end(); mapIter != mapIterEnd; ++mapIter)
        {    
            int primaryPDG(mapIter->first->GetParticleId());

            float contribution(mapIter->second);
            if (contribution >= chargeContributionThreshold2 && primaryPDG != m_targetParticlePDG)
                forwardsNonPrimaryHits++;
            else if (contribution >= chargeContributionThreshold2 && primaryPDG == m_targetParticlePDG)
                forwardsPrimaryHits++;
        }    
    }    

    float forwardsImpurityFraction(((float)forwardsNonPrimaryHits)/forwardsFitCharges.size()), backwardsImpurityFraction(((float)backwardsNonPrimaryHits)/backwardsFitCharges.size());

    if (nSecondaryHits >= 10 || forwardsImpurityFraction > 0.8 || backwardsImpurityFraction > 0.8)
        isTwoParticles= true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::IsParticleContained(const MCParticle* pMCParticle)
{
    const CartesianVector mcVertex(pMCParticle->GetVertex());
    const CartesianVector mcEndpoint(pMCParticle->GetEndpoint());

    if (pMCParticle->GetEnergy() < 0.05)
        return false;

    /*
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    float nHitsInGap(0.f);

    for (const CaloHit* pCaloHit : caloHitList)
    {
        if (LArGeometryHelper::IsInGap(this->GetPandora(), pCaloHit->GetPositionVector(), TPC_VIEW_W, 2.5))
            nHitsInGap += 1.0;
    }

    if (nHitsInGap/caloHitList.size() >= 0.3)
        return false;
    */

    const float eVx(256.35), eVy(233.), eVz(1036.8);
    const float xBorder(10.), yBorder(20.), zBorder(10.);

    if ((mcEndpoint.GetX() < (eVx - xBorder)) && (mcEndpoint.GetX() > xBorder) && (mcEndpoint.GetY() < (eVy / 2. - yBorder)) && (mcEndpoint.GetY() > (-eVy / 2. + yBorder)) && (mcEndpoint.GetZ() < (eVz - zBorder)) && (mcEndpoint.GetZ() > zBorder))
    {
        if (!LArGeometryHelper::IsInGap(this->GetPandora(), mcVertex, TPC_VIEW_W, 0.05*((mcEndpoint - mcVertex).GetMagnitude())) && !LArGeometryHelper::IsInGap(this->GetPandora(), mcEndpoint, TPC_VIEW_W, 0.05*((mcEndpoint - mcVertex).GetMagnitude())))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::IsRecoParticleContained(const ParticleFlowObject* pPfo)
{
    const pandora::Vertex *const pVertex = LArPfoHelper::GetVertex(pPfo);
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    LArTrackStateVector trackStateVector;
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, pVertex, 10, slidingFitPitch, trackStateVector);

    if (trackStateVector.size() < 2)
        return false;

    TrackState firstTrackState(*(trackStateVector.begin())), lastTrackState(trackStateVector.back());
    const pandora::CartesianVector initialPosition(firstTrackState.GetPosition());
    const pandora::CartesianVector endPosition(lastTrackState.GetPosition());

    const float trackLength3D((endPosition - initialPosition).GetMagnitude());

    const pandora::CartesianVector lowYVector(initialPosition.GetY() < endPosition.GetY() ? initialPosition : endPosition);
    const pandora::CartesianVector highYVector(initialPosition.GetY() > endPosition.GetY() ? initialPosition : endPosition);
    const pandora::CartesianVector lowZVector(initialPosition.GetZ() < endPosition.GetZ() ? initialPosition : endPosition);
    const pandora::CartesianVector highZVector(initialPosition.GetZ() > endPosition.GetZ() ? initialPosition : endPosition);

    const float eVx(256.35), eVy(233.), eVz(1036.8);
    const float xBorder(10.), yBorder(20.), zBorder(10.);

    if (!m_cosmic)
    {
        if ((lowZVector.GetX() < (eVx / 2. - xBorder)) && (lowZVector.GetX() > (-eVx / 2. + xBorder)) && (lowZVector.GetY() < (eVy / 2. - yBorder)) && (lowZVector.GetY() > (-eVy / 2. + yBorder)) && (lowYVector.GetZ() < (eVz - zBorder)) && (lowZVector.GetZ() > zBorder) && (highZVector.GetX() < (eVx / 2. - xBorder)) && (highZVector.GetX() > (-eVx / 2. + xBorder)) && (highZVector.GetY() < (eVy / 2. - yBorder)) && (highZVector.GetY() > (-eVy / 2. + yBorder)) && (lowYVector.GetZ() < (eVz - zBorder)) && (highZVector.GetZ() > zBorder))
        {
            if (!LArGeometryHelper::IsInGap3D(this->GetPandora(), lowZVector, TPC_VIEW_W, 0.05*trackLength3D) && !LArGeometryHelper::IsInGap3D(this->GetPandora(), highZVector, TPC_VIEW_W, 0.05*trackLength3D))
                return true;
        }
    }
    else
    { 
        if ((lowYVector.GetX() < (eVx / 2. - xBorder)) && (lowYVector.GetX() > (-eVx / 2. + xBorder)) && (lowYVector.GetY() < (eVy / 2. - yBorder)) && (lowYVector.GetY() > (-eVy / 2. + yBorder)) && (lowYVector.GetZ() < (eVz - zBorder)) && (lowYVector.GetZ() > zBorder))
        {
            if (!LArGeometryHelper::IsInGap3D(this->GetPandora(), lowYVector, TPC_VIEW_W, 0.05*trackLength3D))
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::IsRecoParticleContained(const Cluster* pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFit(pCluster, 20, slidingFitPitch);

    const CartesianVector lowZVector(slidingFit.GetGlobalMinLayerPosition());
    const CartesianVector highZVector(slidingFit.GetGlobalMaxLayerPosition());

    const float eVx(256.35), eVz(1036.8);
    const float xBorder(10.), zBorder(10.);

    if ((lowZVector.GetX() < (eVx / 2. - xBorder)) && (lowZVector.GetX() > (-eVx / 2. + xBorder)) && (lowZVector.GetZ() < (eVz - zBorder)) && (lowZVector.GetZ() > zBorder) && (highZVector.GetX() < (eVx / 2. - xBorder)) && (highZVector.GetX() > (-eVx / 2. + xBorder)) && (lowZVector.GetZ() < (eVz - zBorder)) && (highZVector.GetZ() > zBorder))
    {
        if (!LArGeometryHelper::IsInGap3D(this->GetPandora(), lowZVector, TPC_VIEW_W, 0.05*(highZVector.GetZ() - lowZVector.GetZ())) && !LArGeometryHelper::IsInGap3D(this->GetPandora(), highZVector, TPC_VIEW_W, 0.05*(highZVector.GetZ() - lowZVector.GetZ())))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::IntersectsYFace(const MCParticle* pMCParticle)
{
    CartesianVector mcBeginPoint(pMCParticle->GetVertex()), mcEndpoint(pMCParticle->GetEndpoint());
    CartesianVector lowestPoint(mcBeginPoint.GetY() < mcEndpoint.GetY() ? mcBeginPoint : mcEndpoint), highestPoint(mcBeginPoint.GetY() > mcEndpoint.GetY() ? mcBeginPoint : mcEndpoint);
    float xExtent(highestPoint.GetX() - lowestPoint.GetX()), yExtent(highestPoint.GetY() - lowestPoint.GetY()), zExtent(highestPoint.GetZ() - lowestPoint.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - lowestPoint.GetY());
    CartesianVector yFaceIntersection(lowestPoint.GetX() + xSlope*yDistanceToTravel, 116.5, lowestPoint.GetZ() + zSlope*yDistanceToTravel);

    if (yFaceIntersection.GetX() > 0.0 && yFaceIntersection.GetX() < 128.0 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1037.0)
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::RecoIntersectsYFace(TrackDirectionTool::DirectionFitObject &fitResult)
{
    const pandora::CartesianVector initialPosition(fitResult.GetBeginpoint());
    const pandora::CartesianVector endPosition(fitResult.GetEndpoint());

    const pandora::CartesianVector lowYVector(initialPosition.GetY() < endPosition.GetY() ? initialPosition : endPosition);
    const pandora::CartesianVector highYVector(initialPosition.GetY() > endPosition.GetY() ? initialPosition : endPosition);

    float xExtent(highYVector.GetX() - lowYVector.GetX()), yExtent(highYVector.GetY() - lowYVector.GetY()), zExtent(highYVector.GetZ() - lowYVector.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - lowYVector.GetY());
    CartesianVector yFaceIntersection(lowYVector.GetX() + xSlope*yDistanceToTravel, 116.5, lowYVector.GetZ() + zSlope*yDistanceToTravel);

    if (yFaceIntersection.GetX() > 0.0 && yFaceIntersection.GetX() < 128.0 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1037.0)
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WriteToTree(const Cluster* pCluster, TrackDirectionTool::DirectionFitObject &fitResult, std::string &treeName)
{
    if (!m_writeToTree)
        return;

    CartesianVector xAxis(1.f, 0.f, 0.f), yAxis(0.f, 1.f, 0.f);

    if (!m_data)
    {
        const MCParticle* pMainMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

        CartesianVector mcEndpoint(pMainMCParticle->GetEndpoint());
        CartesianVector mcBeginpoint(pMainMCParticle->GetVertex());
        int neutrinoInduced(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMainMCParticle));

        int mcForwards(mcBeginpoint.GetZ() < mcEndpoint.GetZ() ? 1 : 0);
        int mcDownwards(mcBeginpoint.GetY() > mcEndpoint.GetY() ? 1 : 0);

        CartesianVector mcDirection((mcEndpoint - mcBeginpoint).GetUnitVector());
        float mcPhi(mcDirection.GetOpeningAngle(xAxis)), mcTheta(mcDirection.GetOpeningAngle(yAxis));

        int mcIntersectsYFace(this->IntersectsYFace(MCParticleHelper::GetMainMCParticle(pCluster)) ? 1 : 0);

        bool isClusterTwoParticles(false);
        this->IsClusterTwoParticles(pCluster, fitResult.GetForwardsFitCharges(), fitResult.GetBackwardsFitCharges(), isClusterTwoParticles);
        int mcIsClusterTwoParticles(isClusterTwoParticles? 1 : 0);

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCDownwards", mcDownwards));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCForwards", mcForwards)); 
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCIntersectsYFace", mcIntersectsYFace));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCPhi", mcPhi));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCTheta", mcTheta));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCIsClusterTwoParticles", mcIsClusterTwoParticles));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "NeutrinoInduced", neutrinoInduced));
    }

    CartesianVector recoDirectionVector((fitResult.GetEndpoint() - fitResult.GetBeginpoint()).GetUnitVector());
    float recoPhi(recoDirectionVector.GetOpeningAngle(xAxis)), recoTheta(recoDirectionVector.GetOpeningAngle(yAxis));

    int recoForwards(fitResult.GetEndpoint().GetZ() > fitResult.GetBeginpoint().GetZ() ? 1 : 0); //1 by definition, fill just as sanity check
    int recoDownwards(fitResult.GetEndpoint().GetY() < fitResult.GetBeginpoint().GetY() ? 1 : 0);

    float UpwardsChiSquared(recoDownwards == 0 ? fitResult.GetForwardsChiSquared() : fitResult.GetBackwardsChiSquared()), DownwardsChiSquared(recoDownwards == 1 ? fitResult.GetForwardsChiSquared() : fitResult.GetBackwardsChiSquared());
    float UpwardsChiSquaredPerHit(recoDownwards == 0 ? fitResult.GetForwardsChiSquaredPerHit() : fitResult.GetBackwardsChiSquaredPerHit()), DownwardsChiSquaredPerHit(recoDownwards == 1 ? fitResult.GetForwardsChiSquaredPerHit() : fitResult.GetBackwardsChiSquaredPerHit());
    float UpDownDeltaChiSquaredPerHit(DownwardsChiSquaredPerHit - UpwardsChiSquaredPerHit); 

    float recoLength((fitResult.GetEndpoint() - fitResult.GetBeginpoint()).GetMagnitude());
    int recoIntersectsYFace(this->RecoIntersectsYFace(fitResult) ? 1 : 0);

    //------------------------------------------------------------

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "RecoDownwards", recoDownwards));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "RecoForwards", recoForwards)); 
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "RecoIntersectsYFace", recoIntersectsYFace));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "RecoPhi", recoPhi));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "RecoTheta", recoTheta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "recoLength", recoLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BeginX", fitResult.GetBeginpoint().GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BeginY", fitResult.GetBeginpoint().GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BeginZ", fitResult.GetBeginpoint().GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "EndX", fitResult.GetEndpoint().GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "EndY", fitResult.GetEndpoint().GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "EndZ", fitResult.GetEndpoint().GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "Probability", fitResult.GetProbability()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "Hypothesis", fitResult.GetHypothesis()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MinChiSquaredPerHit", fitResult.GetMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DeltaChiSquaredPerHit", fitResult.GetDeltaChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsChiSquared", fitResult.GetForwardsChiSquared()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsChiSquared", fitResult.GetBackwardsChiSquared()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsChiSquaredPerHit", fitResult.GetForwardsChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsChiSquaredPerHit", fitResult.GetBackwardsChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "UpDownDeltaChiSquaredPerHit", UpDownDeltaChiSquaredPerHit));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "UpwardsChiSquared", UpwardsChiSquared));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DownwardsChiSquared", DownwardsChiSquared));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "UpwardsChiSquaredPerHit", UpwardsChiSquaredPerHit));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DownwardsChiSquaredPerHit", DownwardsChiSquaredPerHit));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "NumberHits", fitResult.GetNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "SplittingBeforeChiSquaredPerHit", fitResult.GetSplitObject().GetBeforeMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "SplittingAfterChiSquaredPerHit", fitResult.GetSplitObject().GetAfterMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "SplittingBeforeNumberHits", fitResult.GetSplitObject().GetBeforeNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "SplittingAfterNumberHits", fitResult.GetSplitObject().GetAfterNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TEFBeforeChiSquaredPerHit", fitResult.GetTEFObject().GetBeforeMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TEFAfterChiSquaredPerHit", fitResult.GetTEFObject().GetAfterMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TEFBeforeNumberHits", fitResult.GetSplitObject().GetBeforeNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TEFAfterNumberHits", fitResult.GetSplitObject().GetAfterNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "FRBeforeChiSquaredPerHit", fitResult.GetFRObject().GetBeforeMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "FRAfterChiSquaredPerHit", fitResult.GetFRObject().GetAfterMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "FRBeforeNumberHits", fitResult.GetSplitObject().GetBeforeNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "FRAfterNumberHits", fitResult.GetSplitObject().GetAfterNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "FileIdentifier", m_fileIdentifier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "EventNumber", m_eventNumber));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), treeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WriteHitToTree(TrackDirectionTool::HitCharge &hitCharge, std::string &treeName, int mcDirection)
{
    if (!m_writeToTree)
        return;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCDirection", mcDirection));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChargeOverWidth", hitCharge.GetChargeOverWidth()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "HitCharge", hitCharge.GetCharge()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "HitWidth", hitCharge.GetHitWidth()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsFitCharge", hitCharge.GetForwardsFitCharge()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsSigma", hitCharge.GetForwardsSigma()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsDelta", hitCharge.GetForwardsDelta()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsChiSquared", hitCharge.GetForwardsChiSquared()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsFitCharge", hitCharge.GetBackwardsFitCharge()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsSigma", hitCharge.GetBackwardsSigma()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsDelta", hitCharge.GetBackwardsDelta()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsChiSquared", hitCharge.GetBackwardsChiSquared()));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), treeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DirectionAnalysisAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputHitListName", m_inputHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListName", m_clusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TargetParticlePDG", m_targetParticlePDG));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ParticleContained", m_particleContained));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Cosmic", m_cosmic));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Data", m_data));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DrawFit", m_drawFit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FileIdentifier", m_fileIdentifier));

    if (m_writeToTree)
    {    
        //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    } 

    AlgorithmTool *pAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackDirection", pAlgorithmTool));

    if (!(this->m_pTrackDirectionTool = dynamic_cast<TrackDirectionTool *>(pAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
