//| Rewrite solubility model utils using Eigen
//| 
//| 

#ifndef SOL_MODEL_UTILS_HPP
#define SOL_MODEL_UTILS_HPP

#include <cstdlib>

#include <vector>
#include <map>
#include <cmath>
#include <assert.h>
#include <list>

#include <algorithm>
#include <iterator>
#include <utility>
#include <memory>
#include <filesystem>

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include <Eigen/Core>
#include <Eigen/Dense>

#pragma once

namespace CONST
{   
    extern double const minDivTol; 
    extern double const gasR;
    extern double const refP;
    extern double const saltV;

    extern double const floatInf;
    extern double const doubleInf;
    extern double const floatNaN;
    extern double const doubleNaN;
}

namespace DATATYPE
{
    typedef std::vector<int> Tensor1DInt;
    typedef std::vector<float> Tensor1DFloat32;
    typedef std::vector<double> Tensor1DFloat64;
    typedef std::vector<std::vector<int>> Tensor2DInt;
    typedef std::vector<std::vector<float>> Tensor2DFloat32;
    typedef std::vector<std::vector<double>> Tensor2DFloat64;
    typedef std::vector<std::string> Tensor1DString;
    typedef std::vector<std::vector<std::string>> Tensor2DString;

    typedef int Index;
    typedef std::string String; 
}

namespace LOG
{
    using namespace DATATYPE; 

    #define __CHECK_POINT__ std::cerr << "THROW CHECK POINT .\n"

    #define __CHECK_POINT_WITH_MSG__( MSG ) std::cerr << "THROW CHECK POINT WITH MSG: " << MSG << "\n"

    #define __VAR_WITH_EXCEPTION__(VAR, MSG) std::cerr << "THE VARIABLE " << VAR << " HAS EXCEPTION: " << MSG << "\n"

    #define __VAR_WITH_CONDITIONS__( VAR, NAME, MSG, CONDITION1, CONDITION2 ) std::cerr << "THE VARIABLE "  << NAME << ", VALUE: " << VAR << " HAS EXCEPTION: " << MSG <<  " AT (P, T) CONDITIONS: " << CONDITION1 << ", " << CONDITION2 << "\n"

    void PrintTensorString1D(std::vector< String > & Tensor);
    void PrintTensorString2D(std::vector< std::vector< String > > & Tensor);

    template<typename T>
    void PrintTensor1D(std::vector < T > Tensor)
    {
        for (DATATYPE::Index i=0; i<Tensor.size(); ++i)
        {
            printf("%f", Tensor[i]); printf(" , ");
        };
        printf("\n");
    };

    template<typename T>
    void PrintTensor2D(std::vector< std::vector < T > > Tensor)
    {
        for (DATATYPE::Index i=0; i<Tensor.size(); ++i)
        {
            for (DATATYPE::Index j=0; j<Tensor[0].size(); ++j)
            {
                printf("%f", Tensor[i][j]); printf(" , ");
            };
            printf("\n");
        };
        printf("\n");
    };
}

namespace UNIT
{   
    using namespace DATATYPE; 

    int ConvertPressureUnit(std::pair<std::string, std::string> unitPair);
    int ConvertTemperatureUnit(std::pair<std::string, std::string> unitPair);
    double ConvertPressure(double P, std::pair<std::string, std::string> unitPair);
    double ConvertTemperature(double T, std::pair<std::string, std::string> unitPair);
}

namespace UTILS
{   
    using namespace DATATYPE; 

    Tensor2DFloat64 ReadMatrixFile(std::ifstream &f, int &row, int &col);
    Tensor2DFloat64 ReadMatrixFileMultiDelimiters(std::ifstream &f, int &row, int &col);

    template< typename T >
    bool FindElementInTensor1D( T elemVal, std::vector<T> & tensor1D )
    {
        auto elemIt = std::find( tensor1D.begin(), tensor1D.end(), elemVal );
        if (elemIt != tensor1D.end() ) { return true; }
        else { return false; };
    };

    template< typename T>
    double ComputeAverageTensor2D( std::vector< std::vector<T> > & tensor2D )
    {
        double totalVal = 0.;
        double totalValid = 0;
        for (int i = 0; i < tensor2D.size(); ++ i )
        {
        for ( int j = 0; j < tensor2D[i].size(); ++ j)
        {   
            if ( tensor2D[i][j] > 0. ) { 
            totalVal += tensor2D[i][j];
            totalValid ++;
            } else { 
            totalVal += 0.;
            }
        };
        };
        return totalVal / totalValid;
    }
    
    double LinearScaler( double inputVal, double lowerBound, double higherBound );
    Tensor1DFloat64 MakePressureData(double initialPressure, double pressureLimit );
    Tensor1DFloat64 CubicEquationSolver(Tensor1DFloat64 cubicCoeffs);
}

#endif