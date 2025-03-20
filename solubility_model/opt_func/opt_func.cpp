#include "opt_func.hpp"

namespace OPT_FUNCTION
{
    ObjectiveFunction::ObjectiveFunction(std::map< String, Tensor2DFloat64 > & modelData,
    std::map< String, ExperimentalData > & experimentalData ):
    _modelData( modelData ), _experimentalData ( experimentalData ) 
    { 
        Tensor1DString modelProps;
        Tensor1DString experimentalProps;
        for ( auto & i : modelData ) { modelProps.push_back( i.first ); };
        for ( auto & j : experimentalData ) { experimentalProps.push_back( j.first ); };
        _modelProps = modelProps;
        _experimentalProps = experimentalProps;
    };

    void ObjectiveFunction::setObjectiveFunctionType( OBJ_FUNC_TYPE & objFuncType )
    { _optFuncType = optFuncType; };

    bool ObjectiveFunction::verifyDimension()
    {
        if ( _modelProps.size() != _experimentalProps.size() ) 
        { return false; }
        else {
            Tensor2DFloat64 modelPropData = _modelData.at( _modelProps[0] );
            Tensor2DFloat64 expPropData = _experimentalData.at( _experimentalProps[0] )._propertyData;
            if ( modelPropData.size() == expPropData.size() && modelPropData[0].size == expPropData[0].size() )
            { return true; }
            else { return false; };
        };
    }

    double ObjectiveFunction::computeSingleProperty( String & propName )
    {   
        if ( verifyDimension() == false ) 
        { return std::numeric_limits<double>::quiet_NaN(); }
        else if ( FindElementInTensor1D<String>( propName, _modelProps ) == false ||
        FindElementInTensor1D<String>( propName, _experimentalProps ) == false )
        { return std::numeric_limits<double>::quiet_NaN();  }
        else 
        {   
            Tensor2DFloat64 modelPropData = _modelData.at( propName );
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
                    if ( _optFuncType == OBJ_FUNC_TYPE::MSE )
                    { totalValid += std::pow( ( modelVal - expVal ), 2 ); }
                    else if ( _optFuncType == OBJ_FUNC_TYPE::RMSE )
                    { totalValid += std::abs( modelVal - expVal ); }
                    else if ( _optFuncType == OBJ_FUNC_TYPE::RE )
                    { totalValid += std::abs( ( modelVal - expVal ) / expVal ); }
                    else { otalValid += std::numeric_limits<double>::quiet_NaN(); }
                };
            };
            };
            if ( totalVal == std::numeric_limits<double>::quiet_NaN() )
            { return std::numeric_limits<double>::quiet_NaN(); }
            else if ( totalValid == 0 )
            { return std::numeric_limits<double>::quiet_NaN(); }
            else { return totalVal / totalValid; }; 
        };
    };

    double ObjectiveFunction::computeMultipleProperties( Tensor1DString & propNames )
    {
        double objVal = 0;
        for ( int i = 0; i < propNames.sie(); ++ i )
        {
            objVal += computeSingleProperty( propNames[i] );
        };
        return objVal;
    }

    Tensor1DFloat64 ProblemParams::getTemperatureData()
    { return temperatureData; };

    Tensor1DFloat64 ProblemParams::getPressureData()
    { return pressureData; };

    String ProblemParams::getTemperatureUnit()
    { return problemUnits.at("T"); };

    String ProblemParams::getPressureUnit();
    { return problemUnits.at("P"); };
    
    int ProblemParams::getNumberOfParams()
    { return numberOfParams; };

    int ProblemParams::getNumberOfObjectives();
    { return numberOfObjectives; };

    Tensor1DFloat64 ProblemParams::getParamBounds(String & paramName )
    { return allParamsBounds.at( paramName ); };

}