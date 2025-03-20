#include "model.hpp"

namespace MODEL
{   
    double ModelParams::geth2ACoeff() { return h2ACoeff; };
    double ModelParams::geth2BCoeff() { return h2BCoeff; };
    double ModelParams::geth2oACoeff() { return h2oACoeff; };
    double ModelParams::geth2oBCoeff() { return h2oBCoeff; };

    Tensor1DFloat64 ModelParams::getMixtureCoefficients()
    {
        Tensor1DFloat64 Val(2, 0.);
        Val[0] = h2ACoeff;
        Val[1] = h2BCoeff;
        return Val;
    };

    std::pair< int, int > ConvertToPropertyIndexes(PROPERTY Property)
    {
        if ( Property == PROPERTY::H2LIQUID )
        { return std::make_pair(1, 0); }
        else if  ( Property == PROPERTY::H2VAPOR )
        { return std::make_pair(0, 0); }
        else if ( Property == PROPERTY::H2OLIQUID )
        { return std::make_pair(1, 1); }
        else if ( Property == PROPERTY::H2OVAPOR )
        { return std::make_pair(0, 1); }
        else 
        { return std::make_pair(-1, -1); };
    };

    double GasPhaseEquilibrium(double T)
    {   
        double k_gas=2.99475+4.81373*pow(10, -3)*T-5.1773*pow(10, -5)*pow(T, 2)+1.19052*pow(10, -7)*pow(T, 3);
        return pow(10, k_gas);
    };

    double LiquidPhaseEquilibrium(double T)
    {
        double k_liquid=-2.209+3.097*pow(10, -2)*T-1.098*pow(10, -4)*pow(T, 2)+2.048*pow(10, -7)*pow(T, 3);
        return pow(10, k_liquid);
    };
    
    double GasActivity(double T, double m_s)
    {   
        // Since m_s > 0, Heaviaside is 1.
        if ( m_s >= 0. )
        {
            double gas_gamma=2*2.59*0.5+(-0.11+0.0022*T)*m_s;
            return exp(gas_gamma);
        } else
        {
            double gas_gamma=2*2.59*-0.5+(-0.11+0.0022*T)*m_s;
            return exp(gas_gamma);
        }
    };

    double DetermineRoot(double T, double P,
    Tensor1DFloat64 Roots, Tensor1DFloat64 mixtureCoeffs)
    {
        assert( Roots.size() == 3 );
        int numInf = 0;
        for (double r : Roots){ if ( std::isinf(r) == 1 ) { numInf += 1; }; };
            if ( numInf == 2)
            {   
            // Only one root satisfies and it is selected
            for (double r : Roots) 
                { if ( std::isinf(r) == 0 ) 
                    { return r; }; 
                };
            } else
            {
            // Determine gas root and liquid root
            Tensor1DFloat64 positiveRoots;
            for (double r : Roots){ if ( r > 0. ) { positiveRoots.push_back(r); };
            if ( positiveRoots.size() == 1 )
            {
                return positiveRoots[0];
            }; 
                auto minIdx = std::min_element( positiveRoots.begin(), positiveRoots.end() );
                double minRoot = positiveRoots[std::distance( positiveRoots.begin(), minIdx)];
                auto maxIdx = std::max_element( positiveRoots.begin(), positiveRoots.end() );
                double maxRoot = positiveRoots[std::distance( positiveRoots.begin(), maxIdx)];
                // Select gas  root or liquid root depending on w1 & w2 
                double a=mixtureCoeffs[0];
                double b=mixtureCoeffs[1];
                double gasRoot = maxRoot; double liquidRoot = minRoot;
                double w1 = P * ( gasRoot - liquidRoot );
                double w2 = CONST::gasR * T * log( (gasRoot - b) /( liquidRoot - b ));
                w2 += a / ( std::pow(T, 0.5) * b ) * log( ((gasRoot+b)*liquidRoot) / ((liquidRoot+b)*gasRoot) );
            if ( w2 - w1 > 0.) 
            { 
                return gasRoot; 
            } else if ( w2 - w1 < 0.) 
            { 
                return liquidRoot; 
            } else 
            { 
                return gasRoot; 
            };
            };
            };
        return std::numeric_limits<double>::quiet_NaN();
    };

    double GasPhaseVolume (double T, double P, Tensor1DFloat64 mixtureCoeffs)
    {   
        Tensor1DFloat64 cubicCoeffs(4, 0.);
        double mixtureA=mixtureCoeffs[0];
        double mixtureB=mixtureCoeffs[1];
        cubicCoeffs[0] = 1.;
        cubicCoeffs[1] = - CONST::gasR * T / P;
        cubicCoeffs[2] = - ( (CONST::gasR *T*mixtureB) / P - mixtureA/(P*pow(T, 0.5)) + pow(mixtureB, 2) );
        cubicCoeffs[3] = - mixtureA * mixtureB / ( P*pow(T, 0.5)) ;
        Tensor1DFloat64 Roots=CubicEquationSolver(cubicCoeffs);
        double Root = DetermineRoot(T, P, Roots, mixtureCoeffs);
        return Root;
    };
    
    double PhaseFugacity(double T, double P, double V, 
    Tensor1DFloat64 mixtureCoeffs, ModelParams modelParams, Index Idx)
    {
        double mixtureA=mixtureCoeffs[0];
        double mixtureB=mixtureCoeffs[1];
        double componentB = 0.;
        if ( Idx == h2Idx )
        { componentB = modelParams.geth2BCoeff(); } 
        else 
        { componentB = modelParams.geth2oBCoeff(); }
    
        double F1 = log(V / (V-mixtureB)) + componentB / (V-mixtureB);
        double F2 = 0.;

        if ( Idx == h2Idx )
        {
            double tmp = modelParams.geth2ACoeff();
            F2 =-2*tmp/(CONST::gasR * mixtureB*pow(T, 1.5))*log((V+mixtureB)/V);
        } else
        {
            double tmp = modelParams.geth2ACoeff() * modelParams.geth2oACoeff();
            tmp = std::sqrt( tmp ); 
            F2 =-2*tmp/(CONST::gasR * mixtureB*pow(T, 1.5))*log((V+mixtureB)/V);
        };

        double F3=(mixtureA*componentB)/(CONST::gasR*pow(T, 1.5)*pow(mixtureB, 2))*( log((V+componentB)/V) - (mixtureB/(V+mixtureB)) );
        double F4=-log((P*V)/(CONST::gasR*T));
    
        return exp(F1+F2+F3+F4);
    };

    double SolubilityParamA(double K, double Phi, double P, double T) 
    {
        return K/(Phi*P)*exp(((P-CONST::refP) * h2oVolume)/(CONST::gasR * T));
    };

    double SolubilityParamB(double K, double Phi, double P, double T, double GasActivityParam) 
    {
        return (Phi*P)/(55.508*GasActivityParam*K)*exp(-(P-CONST::refP) * h2Volume / (CONST::gasR * T));
    };   

    Tensor2DFloat64 ComponentPhaseFractions(double A, double B, double v, double m_s)
    {   
        /// FIXME: Solid phase fractions not included to check model's error
        double C_s = CONST::saltV * m_s;
        Tensor1DFloat64 gasComponentFractions(2, 0.);
        Tensor1DFloat64 liquidComponentFractions(2, 0.);
        Tensor2DFloat64 Val;

        gasComponentFractions[h2oIdx] = ( 1. - B ) / ( 1./A - B );
        gasComponentFractions[h2Idx] = 1. - gasComponentFractions[h2oIdx];

        liquidComponentFractions[h2Idx] = B * ( 1 - gasComponentFractions[h2oIdx] );
        liquidComponentFractions[h2oIdx] = 1. - liquidComponentFractions[h2Idx];

        Val.push_back(gasComponentFractions); // gasIdx=0
        Val.push_back(liquidComponentFractions); // liquidIdx=1
        return Val;
    };

    Tensor2DFloat64 SolubilityModel(double T, double P, double m_s, std::string unitT, std::string unitP,
    ModelParams modelParams)
    {   
        std::pair<std::string, std::string> tempUnitPair = {unitT, "K"};
        std::pair<std::string, std::string> presUnitPair = {unitT, "Bar"};

        double convP = ConvertPressure(P, presUnitPair);
        double convT = ConvertTemperature(T, tempUnitPair);

        if ( tempUnitPair.first.compare( tempUnitPair.second) == 0 )
        { assert( convT == T ) ; };
        if ( presUnitPair.first.compare( presUnitPair.second) == 0 )
        { assert( convP == P ) ; };      


        Tensor1DFloat64 mixtureCoeffs=modelParams.getMixtureCoefficients();
        double V = GasPhaseVolume(convT, convP, mixtureCoeffs);

        tempUnitPair.first = "K"; tempUnitPair.second = "C";
        double convForKCorrT = ConvertTemperature(convT, tempUnitPair);

        double K_h2=GasPhaseEquilibrium(convForKCorrT); // Require unit temperature C
        double K_h2o=LiquidPhaseEquilibrium(convForKCorrT); // Require unit temperature C
        double Phi_h2=PhaseFugacity(convT, convP, V, mixtureCoeffs, modelParams, gasIdx);
        double Phi_h2o=PhaseFugacity(convT, convP, V, mixtureCoeffs, modelParams, liquidIdx);

        double h2Activity=GasActivity(convT, m_s);
        double A=SolubilityParamA(K_h2o, Phi_h2o, convP, convT);
        double B=SolubilityParamB(K_h2, Phi_h2, convP, convT, h2Activity);

        Tensor2DFloat64 phaseCompFrac=ComponentPhaseFractions(A, B, CONST::saltV, m_s);

        return phaseCompFrac;
    }; 

    double ComponentPhaseFractionsWrapper(double T, double P, double m_s, std::string unitT, std::string unitP,
    ModelParams modelParams, PROPERTY Property)
    {   
        Tensor2DFloat64 phaseCompFrac = SolubilityModel(T, P, m_s, unitT, unitP, modelParams);
        std::pair< int, int > propertyIndexes = ConvertToPropertyIndexes( Property );
        double propertyVal = phaseCompFrac[propertyIndexes.first][propertyIndexes.second];
        return propertyVal;
    };

    Tensor2DFloat64 SolubilityModelWrapper(Tensor1DFloat64 & temperatureData, Tensor1DFloat64 & pressureData, 
    ModelParams modelParams, double m_s, std::string temperatureUnit, std::string pressureUnit, 
    PROPERTY Property, PROPERTY_DIMENSION pressureDim, PROPERTY_DIMENSION temperatureDim)
    {
        int _pressureDim = pressureData.size();
        int _temperatureDim = temperatureData.size();

        Tensor2DFloat64 propertyData(1, Tensor1DFloat64(1, doubleNaN ));

        if ( pressureDim == PROPERTY_DIMENSION::ROW && temperatureDim == PROPERTY_DIMENSION::COL ) 
        { 
            propertyData.resize( _pressureDim, Tensor1DFloat64( _temperatureDim, doubleNaN ) );
        } else if ( pressureDim == PROPERTY_DIMENSION::COL && temperatureDim == PROPERTY_DIMENSION::ROW )
        {
            propertyData.resize( _temperatureDim, Tensor1DFloat64( _pressureDim, doubleNaN ) );
        } else 
        { 
            #warning "Can not resize property data based on given dimension properties."
            exit(1);
        };

        if ( pressureDim == PROPERTY_DIMENSION::ROW && temperatureDim == PROPERTY_DIMENSION::COL ) 
        {   
            for (Index ip = 0; ip < _pressureDim; ++ ip)
            {
            for ( Index it = 0; it < _temperatureDim; ++ it)
            {
            double propertyVal = ComponentPhaseFractionsWrapper(temperatureData[it], pressureData[ip], m_s, 
            temperatureUnit, pressureUnit, modelParams, Property);
            propertyData[ip][it] = propertyVal;
            };
            };
        } else 
        {   
            for (Index it = 0; it < _temperatureDim; ++ it)
            {
            for ( Index ip = 0; ip < _pressureDim; ++ ip)
            {
            double propertyVal = ComponentPhaseFractionsWrapper(temperatureData[it], pressureData[ip], m_s, 
            temperatureUnit, pressureUnit, modelParams, Property);
            propertyData[it][ip] = propertyVal;
            };
            };
        };

        std::cout << propertyData.size() << " , " << propertyData[0].size() << "\n";
        return propertyData;
    };


    double SolubilityModelRelativeError(Tensor1DFloat64 & temperatureData, Tensor1DFloat64 & pressureData, 
    ModelParams modelParams, double m_s, std::string temperatureUnit, std::string pressureUnit, 
    Tensor2DFloat64 modelTrueData)
    {   
        assert(modelTrueData.size() == pressureData.size());
        assert(modelTrueData[0].size() == temperatureData.size());
        //
        int totalValidData = 0;
        int totalOutlierData = 0;
        double totalRelativeError = 0.;
        //
        for (Index ip = 0; ip < pressureData.size(); ++ ip)
        {
        for ( Index it = 0; it < temperatureData.size(); ++ it)
        {
            Tensor2DFloat64 phaseCompFrac = SolubilityModel(temperatureData[it], pressureData[ip], m_s, 
            temperatureUnit, pressureUnit, modelParams);
        
            double yH2O = phaseCompFrac[gasIdx][h2oIdx];
            double yH2 = phaseCompFrac[gasIdx][h2Idx];
            double xH2O = phaseCompFrac[liquidIdx][h2oIdx];
            double xH2 = phaseCompFrac[liquidIdx][h2Idx];

            if ( std::abs( modelTrueData[ip][it] + 1. ) > minDivTol ) 
            {   
                double relErr = std::abs( ( yH2O - modelTrueData[ip][it] ) / modelTrueData[ip][it] );                
                totalRelativeError += relErr;
                totalValidData++;
                // 
                if ( 100 * relErr >= 10.0 ) { 
                    totalOutlierData++; 
                };
            } else { 
                double relErr = -1.; 
            };
        };
        };
        //
        double avgRelativeError = totalRelativeError / totalValidData;
        return avgRelativeError;
    };
}