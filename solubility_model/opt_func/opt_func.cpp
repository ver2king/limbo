#include "opt_func.hpp"

namespace OPT_FUNCTION
{
    ObjectiveFunction::ObjectiveFunction(std::map< PROPERTY, Data > & modelData,
    std::map< PROPERTY, Data > & experimentalData ):
    _modelData( modelData ), _experimentalData ( experimentalData ) 
    { 
        std::vector< PROPERTY > modelProps;
        std::vector< PROPERTY > experimentalProps;
        for ( auto & i : modelData ) { modelProps.push_back( i.first ); };
        for ( auto & j : experimentalData ) { experimentalProps.push_back( j.first ); };
        _modelProps = modelProps;
        _experimentalProps = experimentalProps;
    };

    void ObjectiveFunction::setObjectiveFunctionType( OBJ_FUNC_TYPE & objFuncType )
    { _objFuncType = objFuncType; };

    bool ObjectiveFunction::verifyDimension()
    {       
        if ( _modelProps.size() != _experimentalProps.size() ) 
        { return false; }
        else {
            Tensor2DFloat64 modelPropData = _modelData.at( _modelProps[0] )._propertyData;
            Tensor2DFloat64 expPropData = _experimentalData.at( _experimentalProps[0] )._propertyData;
            if ( modelPropData.size() == expPropData.size() && modelPropData[0].size() == expPropData[0].size() )
            { return true; }
            else { return false; };
        };
    };

    double ObjectiveFunction::computeSingleProperty( PROPERTY & propName )
    {   
        if ( verifyDimension() == false ) 
        { return doubleInf; }
        else if ( FindElementInTensor1D<PROPERTY>( propName, _modelProps ) == false ||
        FindElementInTensor1D<PROPERTY>( propName, _experimentalProps ) == false )
        { return doubleInf;  }
        else 
        {   
            Tensor2DFloat64 modelPropData = _modelData.at( propName )._propertyData;
            Tensor2DFloat64 expPropData = _experimentalData.at( propName )._propertyData;
            //
            int rowDim = modelPropData.size();
            int colDim = modelPropData[0].size();
            //
            int totalValid = 0;
            double totalVal = 0.;
            for ( int i = 0; i < rowDim; ++ i )
            {
            for ( int  j = 0; j < colDim; ++ j )
            {
                double modelVal = modelPropData[i][j];
                double expVal = expPropData[i][j];
                if ( modelVal < 0 || expVal < 0 ) { 
                    totalValid += 0; 
                    totalVal += 0.; 
                } else {
                    totalValid += 1;
                    if ( _objFuncType == OBJ_FUNC_TYPE::MSE )
                    { totalValid += std::pow( ( modelVal - expVal ), 2 ); }
                    else if ( _objFuncType == OBJ_FUNC_TYPE::RMSE )
                    { totalValid += std::abs( modelVal - expVal ); }
                    else if ( _objFuncType == OBJ_FUNC_TYPE::RE )
                    { totalValid += std::abs( ( modelVal - expVal ) / expVal ); }
                    else { totalValid += doubleInf; }
                };
            };
            };
            if ( totalVal == doubleInf )
            { return doubleInf; }
            else if ( totalValid == 0 )
            { return doubleInf; }
            else { return totalVal / totalValid; }; 
        };
    };

    double ObjectiveFunction::computeMultipleProperties( std::vector< PROPERTY > & propNames )
    {
        double objVal = 0;
        for ( int i = 0; i < propNames.size(); ++ i )
        {
            objVal += computeSingleProperty( propNames[i] );
        };
        return objVal;
    }

    ProblemParams::ProblemParams(  std::map< PROBLEM_UNIT, String > _problemUnits,
    Tensor1DFloat64 _temperatureData, Tensor1DFloat64 _pressureData, 
    Tensor1DString _paramNames, std::map< String, Tensor1DFloat64 > _allParamsBounds,
    std::vector< PROPERTY > _propNames )
    {
        problemUnits = _problemUnits;
        temperatureData = _temperatureData;
        pressureData = _pressureData;
        paramNames = _paramNames;
        allParamsBounds = _allParamsBounds;
        propNames = _propNames;
    };

    Tensor1DFloat64 ProblemParams::getTemperatureData()
    { return temperatureData; };

    Tensor1DFloat64 ProblemParams::getPressureData()
    { return pressureData; };

    String ProblemParams::getTemperatureUnit()
    { return problemUnits.at( PROBLEM_UNIT::T ); };

    String ProblemParams::getPressureUnit()
    { return problemUnits.at( PROBLEM_UNIT::P ); };

    Tensor1DFloat64 ProblemParams::getParamBounds(String & paramName )
    { return allParamsBounds.at( paramName ); };

    ModelPrediction::ModelPrediction( ProblemParams & _problemParams ):
    problemParams( _problemParams ) { };
}