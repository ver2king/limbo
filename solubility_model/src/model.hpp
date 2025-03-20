//| Rewrite solubility model using Eigen
//| 
//| 

#ifndef SOL_MODEL_HPP
#define SOL_MODEL_HPP

#include "../../solubility_model/data/data.hpp"

#pragma once

namespace MODEL
{
    //| List all anemspaces used
    
    using namespace CONST;
    using namespace DATATYPE;
    using namespace UNIT;
    using namespace UTILS;
    using namespace DATA;

    //| Define indexes of phases and component

    DATATYPE::Index const liquidIdx = 1;
    DATATYPE::Index const gasIdx = 0;
    DATATYPE::Index const h2Idx = 0;
    DATATYPE::Index const h2oIdx = 1;

    double const h2Volume = 26.7;
    double const h2oVolume = 18.1;

    //| Define the model params that all need to fit 
    //| the Bayes Optimizer

    struct ModelParams
    {
        double h2ACoeff;
        double h2BCoeff;
        double h2oACoeff;
        double h2oBCoeff;
        
        // Return specific coeffs
        double geth2ACoeff();
        double geth2BCoeff();
        double geth2oACoeff();
        double geth2oBCoeff();

        /// Return A & B 
        Tensor1DFloat64 getMixtureCoefficients(); 
    };

    std::pair< int, int > ConvertToPropertyIndexes(PROPERTY Property);

    double GasPhaseEquilibrium(double T);
    
    double LiquidPhaseEquilibrium(double T);
    
    double GasActivity(double T, double m_s);

    double DetermineRoot(double T, double P,
    Tensor1DFloat64 Roots, Tensor1DFloat64 mixtureCoeffs);
        
    double GasPhaseVolume (double T, double P, Tensor1DFloat64 mixtureCoeffs);
        
    double PhaseFugacity(double T, double P, double V, Tensor1DFloat64 mixtureCoeffs, 
    ModelParams modelParams, Index Idx);
    
    double SolubilityParamA(double K, double Phi, double P, double T);
        
    double SolubilityParamB(double K, double Phi, double P, double T, double GasActivityParam);

    Tensor2DFloat64 ComponentPhaseFractions(double A, double B, double v, double m_s);

    Tensor2DFloat64 SolubilityModel(double T, double P, double m_s, std::string unitT, std::string unitP,
    ModelParams modelParams);

    double ComponentPhaseFractionsWrapper(double T, double P, double m_s, std::string unitT, std::string unitP,
    ModelParams modelParams, PROPERTY Property);

    Tensor2DFloat64 SolubilityModelWrapper(Tensor1DFloat64 & temperatureData, Tensor1DFloat64 & pressureData, 
    ModelParams modelParams, double m_s, std::string temperatureUnit, std::string pressureUnit, 
    PROPERTY Property, PROPERTY_DIMENSION pressureDim, PROPERTY_DIMENSION temperatureDim);

    double SolubilityModelRelativeError(Tensor1DFloat64 & temperatureData, Tensor1DFloat64 & pressureData, 
    ModelParams modelParams, double m_s, std::string temperatureUnit, std::string pressureUnit, 
    Tensor2DFloat64 modelTrueData);
}

#endif