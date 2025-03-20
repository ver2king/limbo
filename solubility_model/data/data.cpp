#include "data.hpp"

namespace DATA
{
    Data::Data(Tensor1DFloat64 & pressureData, Tensor1DFloat64 & temperatureData, 
    Tensor2DFloat64 & propertyData): _pressureData(pressureData), _temperatureData(temperatureData),
    _propertyData(propertyData) { };
    
    void Data::setPressureDimension( PROPERTY_DIMENSION pressureDim )
    { _pressureDim = pressureDim; };

    void Data::setTemperatureDimension( PROPERTY_DIMENSION temperatureDim )
    { _temperatureDim = temperatureDim; };

    bool Data::verifyDimension()
    { 
        int rowDim = _propertyData.size();
        int colDim = _propertyData[0].size();
        //
        int pressureDim = _pressureData.size();
        int temperatureDim = _temperatureData.size();
        bool rowDimVerify = false;
        bool colDimVerify = false;
        //
        if ( _pressureDim == PROPERTY_DIMENSION::ROW ) { 
            rowDimVerify = ( pressureDim == rowDim );
        } else { 
            colDimVerify = (pressureDim == colDim );
        };
        if ( _temperatureDim == PROPERTY_DIMENSION::ROW ) { 
            colDimVerify = ( temperatureDim == rowDim );
        } else {
            colDimVerify = ( temperatureDim == colDim );
        };
        //
        if ( rowDimVerify == true && colDimVerify == true ) { 
            return true;
        } else { 
            return false;
        };
    };

    Tensor1DFloat64 Data::getPressureData() { return _pressureData; };
    Tensor1DFloat64 Data::getTemperatureData() { return _temperatureData; };

    double Data::getPropertyData( double T, double P ) 
    {   
        //
        double propertyVal = std::numeric_limits<double>::quiet_NaN();
        auto temperatureIt = std::find( _temperatureData.begin(), _temperatureData.end(), T );
        auto pressureIt = std::find( _pressureData.begin(), _pressureData.end(), P );
        if (temperatureIt != _temperatureData.end() && pressureIt != _pressureData.end() ) {
            int temperatureIdx = std::distance( _temperatureData.begin(), temperatureIt );
            int pressureIdx = std::distance( _pressureData.begin(), pressureIt );
            if ( _pressureDim == PROPERTY_DIMENSION::ROW && _temperatureDim == PROPERTY_DIMENSION::COL ) 
            { propertyVal = _propertyData[pressureIdx][temperatureIdx]; } 
            else if ( _pressureDim == PROPERTY_DIMENSION::COL && _temperatureDim == PROPERTY_DIMENSION::ROW ) 
            { propertyVal = _propertyData[temperatureIdx][pressureIdx]; }
            else 
            { propertyVal = std::numeric_limits<double>::quiet_NaN(); };
        } else { 
            propertyVal = std::numeric_limits<double>::quiet_NaN();
        };
        return propertyVal;
    };
}